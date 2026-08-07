## ─────────────────────────────────────────────────────────────────────────────
## helpers_smlm.R  —  SMLM-only helper functions
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   All the computation behind the SMLM (imaging) tab. Unlike the Tracking tab —
##   which asks how molecules MOVE — this tab asks where molecules ARE and how
##   they are ARRANGED: clustering, shape fitting, spatial statistics and
##   two-colour comparison. Nothing here is reactive; every function is a plain
##   input -> output transform testable in the console.
##
## THE DATA HIERARCHY (understand this and the whole file makes sense)
##   localizations          raw points. One BLINK produces many of them.
##      -> tid              trace ID: one continuous blink of one fluorophore
##      -> cluster          "primary" DBSCAN over localizations. Its job is
##                          DE-BLINKING: collapse repeated appearances of the
##                          same molecule into one object.
##      -> cluster2         "median"/secondary clustering, run over the medians
##                          of the level above. Its job is STRUCTURE: group
##                          molecules into the biological object they form.
##   Most functions take a level as input and summarise to the level above, so
##   "which level am I on?" is the first question to ask of any function here.
##
##   WHY THIS MATTERS SCIENTIFICALLY: a fluorophore blinking 10 times looks like
##   10 molecules, and an antibody carries several dyes ~10-40 nm from its
##   target. Feeding RAW localizations into a clustering or spatial statistic
##   therefore measures the LABEL, not the biology. The tiering exists to
##   collapse that artefact before any counting or spatial statistic is applied.
##
## MAP OF THIS FILE (grep the name; line numbers drift, names do not)
##
##   1. AXES & COLOUR
##      get_equal_limits_2d/get_equal_limits ... equal aspect ratio (never let
##                                               plotly distort real geometry)
##      hex_to_rgb, rgb_to_hex, scale_hex_at, colorscale_to_hex, get_palette_hex,
##      palette_colors, plotly_colorscale, map_value_to_scale_hex
##      build_colors_for_pts, color_mapping_for_pts ... per-point colour
##
##   2. TWO-COLOUR / MULTIPLEX
##      assign_two_color_group, auto_detect_dcr_split ... split channels on DCR
##      assign_analysis_group(_from_mplex/_legacy), parse_tc_label, format_tc_id
##      mplex_default_colors
##
##   3. POINT BUILDING & SUMMARISING
##      build_pts_2d, build_pts_3d ..... collapse a level to plottable points
##      median_points_df, medianize_coords, compute_spatial_density
##      remove_unstable_traces, flag_rogue_tids, remove_rogue_traces
##      trace_spread_table(_legacy), sd_* family ... per-trace spread / SD filter
##
##   4. CLUSTERING
##      apply_clustering_pipeline ...... THE main entry point (primary + cluster2)
##      run_cluster2_grouping, .cluster2_core, .prune_sparse_members
##      apply_two_color_clustering_pipeline
##
##   5. NEAREST NEIGHBOUR
##      nn_index_pairs ................. within / between / global channel scope
##      build_nn_edges, pair_passes_scope, compute_batch_nn, compute_cross_nn
##      safe_nn2, safe_k
##
##   6. SHAPE FITTING (the biggest block — see its own banner below)
##      For each shape S in {sphere, ellipsoid, helix, ring, arc}:
##          fit_S_geometry / mle_fit_S ... fit the shape to the points
##          mle_project_S ................ nearest point ON the shape (residuals)
##          mle_sample_S ................. sample the surface (for drawing)
##          make_S_mesh / make_S_tube_mesh  3D geometry for plotly
##          generate_S_from_cluster2 ..... convenience wrapper per cluster
##      calculate_distance_to_fit, add_distance_lines_to_plot
##
##   7. GROUND-TRUTH ESTIMATION (flagged EXPLORATORY in the UI)
##      estimate_gt_from_labeling ...... inverse of the labelling model
##      gt_uncertainty_colors, analyze_cluster2_morphology,
##      mst_contour_length, .grid_thin
##
##   8. SPATIAL STATISTICS
##      compute_cluster_rog ............ radius of gyration per cluster
##      compute_ripley_L, compute_ripley_L12, compute_pcf12,
##      compute_ripley_envelope ........ CSR envelope by randomisation
##      compute_kde_coloc .............. 2D KDE + Pearson co-localization
##
##   9. SURFACES / PLANES
##      fit_poly_surface, create_global_surface, generate_xy_plane,
##      compute_cluster2_plane_distances, estimate_cluster2_curvature,
##      build_curvature_table, add_cluster2_lines_and_angles
##
## UNITS: nanometres throughout. `tim` is seconds in the SMLM tab (the Tracking
## tab converts to ms; this tab does not).
##
## Method provenance and the MINFLUX-appropriateness of each spatial statistic
## are documented in docs/MINFLUX_analysis_methods.md (Part II §6).
## ─────────────────────────────────────────────────────────────────────────────

# ── get_equal_limits_2d() ────────────────────────────────────────────────────
# GOAL: return x/y axis limits that span the SAME distance, so a nanometre on
# screen is the same length on both axes.
# This is not cosmetic: with unequal limits plotly stretches one axis, and a
# round cluster renders as an ellipse — which would corrupt every visual shape
# judgement the user makes. Padding adds a small margin so points do not sit on
# the frame.
get_equal_limits_2d <- function(x, y, padding = 0.05) {
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  
  # Handle degenerate ranges (all points same)
  if (!is.finite(diff(xr)) || diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  if (!is.finite(diff(yr)) || diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
  
  span <- max(diff(xr), diff(yr))
  pad  <- span * padding
  
  x_mid <- mean(xr); y_mid <- mean(yr)
  x_lim <- c(x_mid - span/2 - pad, x_mid + span/2 + pad)
  y_lim <- c(y_mid - span/2 - pad, y_mid + span/2 + pad)
  
  list(x = x_lim, y = y_lim)
}
assign_two_color_group <- function(df, input) {
  if (!"z" %in% names(df)) df$z <- 0
  .dcr_col <- if (isTRUE(input$use_pooled_dcr %||% FALSE) && "dcr_pooled" %in% names(df))
    "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL
  med <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      x_med = stats::median(x, na.rm = TRUE),
      y_med = stats::median(y, na.rm = TRUE),
      z_med = stats::median(z, na.rm = TRUE),
      dcr_med  = if (!is.null(.dcr_col)) stats::median(.data[[.dcr_col]], na.rm = TRUE) else NA_real_,
      fluo_med = if ("fluo" %in% names(df))  stats::median(fluo, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  use_fluo <- isTRUE(input$color_by_fluo %||% FALSE)
  
  if (use_fluo) {
    # FLUO override: expect fluo to be 1/2 (or close); map safely
    med <- med %>%
      dplyr::mutate(
        analysis_group = dplyr::case_when(
          is.finite(fluo_med) & round(fluo_med) == 1L ~ 1L,
          is.finite(fluo_med) & round(fluo_med) == 2L ~ 2L,
          TRUE ~ NA_integer_
        )
      )
  } else {
    # DCR ranges
    g1 <- input$dcr_group1; g2 <- input$dcr_group2
    # fallbacks
    if (is.null(g1) || length(g1)!=2) g1 <- c(0, 0.369)
    if (is.null(g2) || length(g2)!=2) g2 <- c(0.3691, 1)
    # choose vector v
    v <- med$dcr_med
    
    med$analysis_group <- dplyr::case_when(
      is.finite(v) & v >= g1[1] & v <  g1[2] ~ 1L,
      is.finite(v) & v >= g2[1] & v <= g2[2] ~ 2L,
      TRUE ~ NA_integer_
    )
  }
  med
}
# Helper: auto-detect 2-color DCR split via 2-component GMM (mclust)
# Returns a list with $success, $split, $mean1, $mean2, $sd1, $sd2, $n_traces
# On failure returns $success = FALSE and $reason string.
# GOAL: split traces into two colour channels on their DCR (detection channel
# ratio) by fitting a 2-component Gaussian mixture, AND quantify how reliable that
# split is.
#
# WHY THE UNCERTAINTY MATTERS (and why it is not symmetric with a registration
# error): ratiometric MINFLUX has no inter-channel offset — both dyes are excited
# by the same laser through the same beam path, so the two channels share one
# coordinate frame. The residual two-colour error is therefore NOT a shift, it is
# MIS-ASSIGNMENT: the DCR distributions of the two dyes overlap, so a fraction of
# traces lands on the wrong side of the threshold. A mis-assigned molecule sits by
# construction right next to its true neighbours, so crosstalk pushes every
# cross-channel statistic (NN, L12, g12) UP at short r — toward FALSE
# co-localization, exactly where a real interaction would appear. Quantifying it is
# what lets a user judge whether a co-localization result is real.
#
# `col_used` records whether the caller passed raw `dcr` or pooled `dcr_pooled`;
# pooling DCR across the iteration sequence (weighted by photons) usually sharpens
# the two modes considerably, so the split quality is only interpretable alongside it.
auto_detect_dcr_split <- function(dcr_meds, col_used = NA_character_) {
  x <- dcr_meds[is.finite(dcr_meds)]
  if (length(x) < 10)
    return(list(success = FALSE, reason = paste0("Not enough traces (", length(x), "; need \u2265 10)")))

  # Guard: if DCR has no meaningful spread, there is no bimodal distribution to split
  # (e.g. only one detector active → all traces have DCR = 1)
  dcr_range <- diff(range(x))
  if (dcr_range < 0.05)
    return(list(success = FALSE,
                reason = paste0("DCR values have no spread (range = ", round(dcr_range, 4),
                                "). Are both detectors (near + far) active?")))

  fit <- tryCatch(
    mclust::Mclust(x, G = 2, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit))
    return(list(success = FALSE, reason = "GMM fitting failed (data may be unimodal)"))

  # Sort components by mean so group 1 = lower DCR
  ord  <- order(fit$parameters$mean)
  mu1  <- fit$parameters$mean[ord[1]]
  mu2  <- fit$parameters$mean[ord[2]]

  # Guard: if both components landed on nearly the same mean, the distribution is unimodal
  if ((mu2 - mu1) < 0.05)
    return(list(success = FALSE,
                reason = paste0("GMM found two components with nearly identical means (",
                                round(mu1, 3), " and ", round(mu2, 3),
                                "). DCR distribution appears unimodal — is this 2-color data?")))
  pi1  <- fit$parameters$pro[ord[1]]
  pi2  <- fit$parameters$pro[ord[2]]

  # Extract SDs — mclust stores variance differently depending on model type
  vp <- fit$parameters$variance
  if (!is.null(vp$sigmasq)) {
    # "E" (equal) → scalar; "V" (variable) → 2-vector
    sv <- vp$sigmasq
    sd1 <- sqrt(if (length(sv) >= 2) sv[ord[1]] else sv)
    sd2 <- sqrt(if (length(sv) >= 2) sv[ord[2]] else sv)
  } else if (!is.null(vp$sigma)) {
    # 1-D array form: sigma[1,1,k]
    sd1 <- sqrt(vp$sigma[1, 1, ord[1]])
    sd2 <- sqrt(vp$sigma[1, 1, ord[2]])
  } else {
    sd1 <- sd2 <- 0.05  # safe fallback
  }

  # Find the crossing point of the two weighted Gaussians between the two means
  f <- function(v) pi1 * dnorm(v, mu1, sd1) - pi2 * dnorm(v, mu2, sd2)
  split_val <- tryCatch(
    uniroot(f, interval = c(mu1, mu2))$root,
    error = function(e) (mu1 + mu2) / 2
  )

  # ── Mis-assignment quantification at the chosen threshold ──────────────────
  # Component 1 sits below the split, component 2 above, so at threshold s:
  #   traces from comp 1 wrongly called "group 2" = pi1 * P(X1 > s)
  #   traces from comp 2 wrongly called "group 1" = pi2 * P(X2 <= s)
  # err_rate is the total expected wrong-channel fraction; contam1/contam2 are the
  # more useful per-channel numbers — what proportion of everything you CALL
  # group k actually came from the other dye. Those are the values that propagate
  # into a cross-channel co-localization claim.
  p1_up   <- pi1 * stats::pnorm(split_val, mu1, sd1, lower.tail = FALSE)  # comp1 -> grp2
  p2_dn   <- pi2 * stats::pnorm(split_val, mu2, sd2, lower.tail = TRUE)   # comp2 -> grp1
  p1_dn   <- pi1 * stats::pnorm(split_val, mu1, sd1, lower.tail = TRUE)   # comp1 -> grp1 (correct)
  p2_up   <- pi2 * stats::pnorm(split_val, mu2, sd2, lower.tail = FALSE)  # comp2 -> grp2 (correct)
  err_rate <- p1_up + p2_dn
  contam1  <- if ((p1_dn + p2_dn) > 0) p2_dn / (p1_dn + p2_dn) else NA_real_
  contam2  <- if ((p2_up + p1_up) > 0) p1_up / (p2_up + p1_up) else NA_real_

  # Standardised separation of the two modes. Rule of thumb: d' >= 3 is a clean
  # split, ~2 is workable, < 1.5 means the "two colours" heavily overlap and any
  # short-range cross-channel signal is suspect.
  #
  # ⚠ d' MEASURES SEPARATION, NOT WHETHER THE DATA IS TWO-COLOUR AT ALL. Asking for
  # G = 2 always returns two components: a single-dye DCR distribution is skewed, and
  # a 2-component Gaussian mixture fits a skewed unimodal distribution better than one
  # Gaussian. Measured here: on a confirmed SINGLE-colour file this routine still
  # returned a "workable" split (d' = 2.65, 3.6% apparent crosstalk), and a BIC
  # comparison of G=1 vs G=2 does NOT fix it — it also preferred two components on
  # that file (dBIC +16.9) and on a synthetic skewed-unimodal control (+21.1), while
  # missing a genuine two-colour pair at d' = 1 (-9.6). BIC is picking up skewness,
  # not bimodality, so it is deliberately NOT used as a gate here. The only valid
  # basis for splitting is the EXPERIMENTAL DESIGN plus a look at the DCR histogram.
  dprime <- abs(mu2 - mu1) / sqrt((sd1^2 + sd2^2) / 2)

  # Per-trace posterior of the assigned component. Traces whose best posterior is
  # below `post_cut` are genuinely ambiguous — they are the ones that would move if
  # the threshold shifted slightly, and the natural candidates for a strict gate.
  post_cut <- 0.95
  d1 <- pi1 * stats::dnorm(x, mu1, sd1)
  d2 <- pi2 * stats::dnorm(x, mu2, sd2)
  tot <- d1 + d2
  post_max <- ifelse(tot > 0, pmax(d1, d2) / tot, NA_real_)
  ambiguous_frac <- mean(post_max < post_cut, na.rm = TRUE)

  list(
    success  = TRUE,
    split    = round(split_val, 4),
    mean1    = round(mu1, 3),
    mean2    = round(mu2, 3),
    sd1      = round(sd1, 3),
    sd2      = round(sd2, 3),
    pi1      = round(pi1, 3),
    pi2      = round(pi2, 3),
    n_traces = length(x),
    col_used = col_used,
    # Unrounded parameters for arithmetic. The fields above are rounded for display;
    # using them to recompute posteriors made the per-trace edge count disagree with
    # ambiguous_frac (rounding sd to 3 dp is ~1% relative error at DCR sd ≈ 0.05).
    params = list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2,
                  pi1 = pi1, pi2 = pi2, split = split_val),
    # split-quality block
    err_rate       = round(err_rate, 4),
    contam1        = round(contam1, 4),
    contam2        = round(contam2, 4),
    dprime         = round(dprime, 2),
    ambiguous_frac = round(ambiguous_frac, 4),
    post_cut       = post_cut,
    quality = if (!is.finite(dprime)) "unknown"
              else if (dprime >= 3) "clean"
              else if (dprime >= 2) "workable"
              else "poor"
  )
}

# GOAL: per-molecule confidence in its colour-channel assignment.
#
# The GMM split gives one threshold for the whole dataset, but molecules are not
# equally certain: one sitting far out in a mode is safely assigned, one sitting in
# the overlap region could go either way. This turns the fitted mixture into a
# per-trace number so the ambiguous molecules can be identified individually rather
# than only counted in aggregate.
#
#   post_k(v) = pi_k * N(v | mu_k, sd_k) / sum_j pi_j * N(v | mu_j, sd_j)
#
# Returns, per input DCR value:
#   group     1 or 2 — the component with the higher posterior (matches the split
#             threshold by construction, since the threshold IS the crossing point)
#   post      posterior of the assigned component, in [0.5, 1]. 0.5 = the molecule
#             sits exactly on the split and is a coin flip; ~1 = unambiguous.
#   edge      TRUE when post < post_cut, i.e. the molecule is in the overlap band
#             and would switch channel under a small shift of the threshold.
#   log_odds  log(post_assigned / post_other) — the same information on an additive
#             scale, which is easier to threshold consistently across datasets.
#
# `gmm` is the list returned by auto_detect_dcr_split().
dcr_split_confidence <- function(v, gmm, post_cut = NULL) {
  n <- length(v)
  na_out <- list(group = rep(NA_integer_, n), post = rep(NA_real_, n),
                 edge = rep(NA, n), log_odds = rep(NA_real_, n))
  if (is.null(gmm) || !isTRUE(gmm$success)) return(na_out)
  if (is.null(post_cut)) post_cut <- gmm$post_cut %||% 0.95

  # Prefer the unrounded parameters; fall back to the display fields for any older
  # stored result that predates them.
  pr <- gmm$params %||% list(mu1 = gmm$mean1, mu2 = gmm$mean2, sd1 = gmm$sd1,
                             sd2 = gmm$sd2, pi1 = gmm$pi1, pi2 = gmm$pi2,
                             split = gmm$split)
  v  <- as.numeric(v)
  d1 <- pr$pi1 * stats::dnorm(v, pr$mu1, pr$sd1)
  d2 <- pr$pi2 * stats::dnorm(v, pr$mu2, pr$sd2)
  tot <- d1 + d2

  # Underflow guard: far out in the tails both densities can be exactly 0, which
  # would give 0/0. Fall back to the side of the threshold there — the assignment
  # is certain, it is only the arithmetic that degenerates.
  bad <- !is.finite(tot) | tot <= 0
  p1  <- ifelse(bad, ifelse(v <= pr$split, 1, 0), d1 / tot)
  p2  <- 1 - p1

  grp  <- ifelse(p1 >= p2, 1L, 2L)
  post <- pmax(p1, p2)
  lo   <- log(pmax(post, .Machine$double.xmin)) - log(pmax(1 - post, .Machine$double.xmin))
  grp[!is.finite(v)]  <- NA_integer_
  post[!is.finite(v)] <- NA_real_
  lo[!is.finite(v)]   <- NA_real_

  list(group = grp, post = post, edge = post < post_cut, log_odds = lo)
}

# Helper: parse TC label like "A-1" or "B_2" into integer group IDs
parse_tc_label <- function(lbl) {
  if (is.null(lbl) || !nzchar(lbl)) return(NA_integer_)
  lbl <- trimws(lbl)
  
  m <- regexec("^([ABab])\\s*[-_]?\\s*(\\d+)$", lbl)
  g <- regmatches(lbl, m)[[1]]
  if (length(g) != 3) return(NA_integer_)
  
  grp <- toupper(g[2])
  n   <- as.integer(g[3])
  if (!is.finite(n) || n < 1) return(NA_integer_)
  
  if (grp == "A") return(n)
  return(10000L + n)  # B group offset
}
# Helper: assign analysis groups based on DCR or FLUO medians
assign_analysis_group <- function(df, input) {
  .dcr_col <- if (isTRUE(input$use_pooled_dcr %||% FALSE) && "dcr_pooled" %in% names(df))
    "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL
  med_tbl <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      dcr_med  = if (!is.null(.dcr_col)) stats::median(.data[[.dcr_col]], na.rm = TRUE) else NA_real_,
      fluo_med = if ("fluo" %in% names(df)) stats::median(fluo, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  df <- dplyr::left_join(df, med_tbl, by = "tid")
  
  mode <- if (is.null(input$grouping_mode)) "none" else input$grouping_mode
  if (isTRUE(input$color_by_fluo %||% FALSE)) mode <- "fluo"
  if (mode == "fluo") {
    df$analysis_group <- dplyr::case_when(
      is.finite(df$fluo_med) & round(df$fluo_med) == 1 ~ 1L,
      is.finite(df$fluo_med) & round(df$fluo_med) == 2 ~ 2L,
      TRUE ~ NA_integer_
    )
  } else if (mode == "dcr") {
    g1 <- input$dcr_group1; g2 <- input$dcr_group2
    if (!is.null(g1) && !is.null(g2) && length(g1) == 2 && length(g2) == 2 && all(is.finite(c(g1, g2)))) {
      v <- df$dcr_med
      df$analysis_group <- dplyr::case_when(
        v >= g1[1] & v <  g1[2] ~ 1L,
        v >= g2[1] & v <= g2[2] ~ 2L,
        TRUE ~ NA_integer_
      )
    } else {
      df$analysis_group <- NA_integer_
    }
  } else {
    df$analysis_group <- NA_integer_
  }
  df
}
# Helper: safe wrapper around RANN::nn2 to avoid crashes on bad input
build_nn_edges <- function(geo_df, k = 4,
                           scope = c("all","within","between"),
                           dist_range = c(0, Inf)) {
  scope <- match.arg(scope)
  if (!all(c("x","y","z") %in% names(geo_df))) return(data.frame())
  if (nrow(geo_df) < 2) return(data.frame())
  
  coords <- as.matrix(geo_df[, c("x","y","z")])
  tmp    <- safe_nn2(coords, k = k + 1)        # self + k neighbors
  if (is.null(tmp$nn)) return(data.frame())
  
  nn      <- tmp$nn
  usedK   <- max(0, tmp$k_eff - 1)
  if (usedK == 0) return(data.frame())
  
  # collect unique undirected edges (i, j)
  seen <- new.env(parent = emptyenv())
  add_edge <- function(i, j, d) {
    ii <- min(i, j); jj <- max(i, j)
    key <- paste0(ii, "_", jj)
    if (!is.null(seen[[key]])) return(NULL)
    seen[[key]] <- TRUE
    
    g1 <- geo_df$analysis_group[ii] %||% NA_integer_
    g2 <- geo_df$analysis_group[jj] %||% NA_integer_
    if (!pair_passes_scope(g1, g2, scope)) return(NULL)
    if (d < dist_range[1] || d > dist_range[2]) return(NULL)
    
    data.frame(
      i = ii, j = jj,
      x1 = coords[ii,1], y1 = coords[ii,2], z1 = coords[ii,3],
      x2 = coords[jj,1], y2 = coords[jj,2], z2 = coords[jj,3],
      g1 = g1, g2 = g2, dist = d
    )
  }
  
  edges <- list()
  for (i in seq_len(nrow(coords))) {
    for (kk in 1:usedK) {
      j <- nn$nn.idx[i, kk + 1]
      if (is.na(j)) next
      d <- as.numeric(nn$nn.dists[i, kk + 1])
      e <- add_edge(i, j, d)
      if (!is.null(e)) edges[[length(edges) + 1]] <- e
    }
  }
  if (!length(edges)) return(data.frame())
  dplyr::bind_rows(edges)
}
# Helper: check if a pair of groups passes the given scope filter
pair_passes_scope <- function(g1, g2, scope = c("all","within","between")) {
  scope <- match.arg(scope)
  
  # "all" means ignore grouping completely
  if (scope == "all") return(TRUE)
  
  # for within/between, require valid groups
  if (is.na(g1) || is.na(g2)) return(FALSE)
  if (!(g1 %in% c(1L, 2L)) || !(g2 %in% c(1L, 2L))) return(FALSE)
  
  if (scope == "within")  return(g1 == g2)
  # scope == "between"
  g1 != g2
}

# pick the nearest color stop from a Plotly colorscale (no interpolation needed)
hex_to_rgb <- function(hex) {
  hex <- toupper(hex)
  if (!grepl("^#", hex)) hex <- paste0("#", hex)
  grDevices::col2rgb(hex)
}

rgb_to_hex <- function(r, g, b) {
  grDevices::rgb(r, g, b, maxColorValue = 255)
}

scale_hex_at <- function(t01, colorscale) {
  # colorscale: list(list(pos, hex), ...)
  t01 <- max(0, min(1, as.numeric(t01)))
  
  pos  <- vapply(colorscale, function(x) as.numeric(x[[1]]), numeric(1))
  cols <- vapply(colorscale, function(x) as.character(x[[2]]), character(1))
  
  ok <- is.finite(pos) & nzchar(cols)
  pos <- pos[ok]; cols <- cols[ok]
  if (!length(pos)) return("#BDBDBD")
  
  # sort by position
  o <- order(pos)
  pos <- pos[o]; cols <- cols[o]
  
  # clamp ends
  if (t01 <= pos[1]) return(cols[1])
  if (t01 >= pos[length(pos)]) return(cols[length(cols)])
  
  # find interval
  hi <- which(pos >= t01)[1]
  lo <- hi - 1L
  
  p0 <- pos[lo]; p1 <- pos[hi]
  c0 <- cols[lo]; c1 <- cols[hi]
  
  # linear interpolation in RGB
  w <- if (p1 == p0) 0 else (t01 - p0) / (p1 - p0)
  
  rgb0 <- hex_to_rgb(c0)
  rgb1 <- hex_to_rgb(c1)
  
  r <- round((1 - w) * rgb0[1,1] + w * rgb1[1,1])
  g <- round((1 - w) * rgb0[2,1] + w * rgb1[2,1])
  b <- round((1 - w) * rgb0[3,1] + w * rgb1[3,1])
  
  rgb_to_hex(r, g, b)
}
# Map a numeric value v in [cmin,cmax] to a hex color via the given colorscale
map_value_to_scale_hex <- function(v, cmin, cmax, colorscale) {
  if (!is.finite(v) || !is.finite(cmin) || !is.finite(cmax) || cmin == cmax) {
    return(scale_hex_at(0.5, colorscale))
  }
  t01 <- (v - cmin) / (cmax - cmin)
  t01 <- max(0, min(1, t01))
  scale_hex_at(t01, colorscale)
}

# Determine hull color for a given group ID based on the provided cmap and pts table
hull_color_for_gid <- function(gid, id_col, pts, input, cmap) {
  # Fallback
  if (is.null(pts) || !nrow(pts) || is.null(cmap)) return("#bdbdbd")
  
  # Find rows for this group in *pts* (the plotted table)
  sel <- if (is.numeric(pts[[id_col]]) || is.integer(pts[[id_col]])) {
    pts[[id_col]] == suppressWarnings(as.numeric(gid))
  } else {
    as.character(pts[[id_col]]) == as.character(gid)
  }
  sel <- which(sel)
  if (!length(sel)) return("#bdbdbd")
  
  # ----- HEX mode: exact discrete mapping -----
  if (identical(cmap$kind, "hex")) {
    if (!is.null(cmap$hex) && length(cmap$hex) == nrow(pts)) {
      hx <- cmap$hex[sel][1]
      if (!is.null(hx) && nzchar(hx)) return(hx)
    }
    return("#bdbdbd")
  }
  
  # ----- SCALE mode: choose a representative scalar for the group -----
  if (identical(cmap$kind, "scale")) {
    v <- cmap$colorvar
    # In your code, cmap$colorvar is either a numeric vector OR a formula-like mapping.
    # Normalize to a numeric vector:
    if (is.function(v)) {
      vv <- v(pts)
    } else if (is.character(v) && length(v) == 1 && v %in% names(pts)) {
      vv <- pts[[v]]
    } else {
      vv <- v
    }
    
    vv <- suppressWarnings(as.numeric(vv))
    if (length(vv) != nrow(pts)) return("#bdbdbd")
    
    gvals <- vv[sel]
    gvals <- gvals[is.finite(gvals)]
    if (!length(gvals)) return("#bdbdbd")
    
    g_rep <- stats::median(gvals)  # <-- key: representative value for the group
    
    # Convert representative scalar into a hex using the *same* colorscale/cmin/cmax
    return(colorscale_to_hex(
      value = g_rep,
      colorscale = cmap$colorscale,
      cmin = cmap$cmin,
      cmax = cmap$cmax
    ))
  }
  
  "#bdbdbd"
}
# Map a numeric value to a hex color via a Plotly colorscale with interpolation
colorscale_to_hex <- function(value, colorscale, cmin, cmax) {
  if (!is.finite(value) || !is.finite(cmin) || !is.finite(cmax) || cmax <= cmin) return("#bdbdbd")
  
  t <- (value - cmin) / (cmax - cmin)
  t <- max(0, min(1, t))
  
  # colorscale is list of (pos,color) pairs, where pos in [0,1], color like "#RRGGBB"
  # Find segment
  pos <- vapply(colorscale, function(x) as.numeric(x[[1]]), numeric(1))
  col <- vapply(colorscale, function(x) as.character(x[[2]]), character(1))
  
  if (length(pos) < 2) return(col[1])
  
  idx <- findInterval(t, pos, all.inside = TRUE)
  if (idx >= length(pos)) return(col[length(col)])
  
  t0 <- pos[idx]; t1 <- pos[idx + 1]
  c0 <- col[idx]; c1 <- col[idx + 1]
  if (!is.finite(t0) || !is.finite(t1) || t1 <= t0) return(c0)
  
  a <- (t - t0) / (t1 - t0)
  
  # linear interpolate rgb
  rgb0 <- grDevices::col2rgb(c0)
  rgb1 <- grDevices::col2rgb(c1)
  rgbm <- rgb0 + (rgb1 - rgb0) * a
  grDevices::rgb(rgbm[1]/255, rgbm[2]/255, rgbm[3]/255)
}

# Color palette utilities for heatmaps and numeric mapping
NUMERIC_PALETTES <- c(
  "Viridis", "Plasma", "Inferno", "Magma", "Cividis",
  "Turbo",
  "Spectral", "RdBu", "RdYlBu", "RdYlGn",
  "YlOrRd", "Greys"
)
# Return a vector of n hex colors for the given palette name
get_palette_hex <- function(name, n = 256, reverse = FALSE) {
  name <- if (is.null(name) || !nzchar(name)) "Viridis" else as.character(name)
  
  cols <- switch(
    name,
    # viridisLite family
    "Viridis" = viridisLite::viridis(n),
    "Plasma"  = viridisLite::plasma(n),
    "Inferno" = viridisLite::inferno(n),
    "Magma"   = viridisLite::magma(n),
    "Cividis" = viridisLite::cividis(n),
    
    "Turbo" = tryCatch(
      viridisLite::turbo(n),
      error = function(e1) tryCatch(
        grDevices::hcl.colors(n, "turbo"),
        error = function(e2) viridisLite::viridis(n)
      )
    ),
    
    # Brewer (diverging/sequential)
    "Spectral" = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n),
    "RdBu"     = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(n),
    "RdYlBu"   = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n),
    "RdYlGn"   = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(n),
    "YlOrRd"   = colorRampPalette(RColorBrewer::brewer.pal(9,  "YlOrRd"))(n),
    "Greys"    = colorRampPalette(RColorBrewer::brewer.pal(9,  "Greys"))(n),
    
    # default
    viridisLite::viridis(n)
  )
  
  if (isTRUE(reverse)) cols <- rev(cols)
  cols
}

# --- multiplex helpers ---------------------------------------------------------
# n distinct categorical colors for multiplex channels. Uses the same qualitative
# Brewer palettes the rest of the SMLM tab uses for trace/cluster coloring so the
# look is consistent, recycling if more channels than palette entries.
mplex_default_colors <- function(n) {
  if (!is.finite(n) || n < 1) return(character(0))
  qual <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
  base_colors <- unlist(mapply(RColorBrewer::brewer.pal,
                               qual$maxcolors, rownames(qual)))
  base_colors <- unname(base_colors)
  rep(base_colors, length.out = n)[seq_len(n)]
}

# Nearest-neighbour index pairs from a coordinate matrix. Returns a 3-col matrix
# (i, ni, dist); indices are into the original coords rows.
#   channel = NULL          -> global k-NN (self excluded).
#   channel + mode="within" -> per-channel k-NN (neighbours only from same channel).
#   channel + mode="between"-> cross-channel k-NN: each point's nearest neighbours
#                              are drawn ONLY from OTHER channels (never its own).
# These are true searches (not post-hoc filters of a global search).
nn_index_pairs <- function(coords, k, channel = NULL, mode = "within") {
  n <- nrow(coords)
  if (is.null(n) || n < 2L) return(matrix(numeric(0), ncol = 3))
  kk <- as.integer(k)

  # within one index set (self excluded)
  build <- function(idx) {
    if (length(idx) < 2L) return(NULL)
    sub   <- coords[idx, , drop = FALSE]
    k_eff <- min(kk + 1L, length(idx))
    if (k_eff < 2L) return(NULL)
    nn <- RANN::nn2(sub, k = k_eff)
    out <- vector("list", length(idx))
    for (a in seq_along(idx)) {
      pr <- lapply(2:k_eff, function(b) {
        ni <- nn$nn.idx[a, b]
        if (is.na(ni) || ni < 1L) return(NULL)
        c(idx[a], idx[ni], nn$nn.dists[a, b])
      })
      out[[a]] <- do.call(rbind, pr)
    }
    do.call(rbind, out)
  }

  if (is.null(channel)) {
    m <- build(seq_len(n))
    return(if (is.null(m)) matrix(numeric(0), ncol = 3) else m)
  }

  if (identical(mode, "between")) {
    chs <- sort(unique(channel[!is.na(channel)]))
    res <- list()
    for (chv in chs) {
      qi <- which(!is.na(channel) & channel == chv)   # query points (this channel)
      di <- which(!is.na(channel) & channel != chv)   # data points (other channels)
      if (!length(qi) || !length(di)) next
      k_eff <- min(kk, length(di))
      if (k_eff < 1L) next
      nn <- RANN::nn2(coords[di, , drop = FALSE],
                      query = coords[qi, , drop = FALSE], k = k_eff)
      for (a in seq_along(qi)) {
        for (b in seq_len(k_eff)) {
          ni <- nn$nn.idx[a, b]
          if (is.na(ni) || ni < 1L) next
          res[[length(res) + 1L]] <- c(qi[a], di[ni], nn$nn.dists[a, b])
        }
      }
    }
    m <- do.call(rbind, res)
    return(if (is.null(m)) matrix(numeric(0), ncol = 3) else m)
  }

  # mode == "within"
  m <- do.call(rbind, lapply(split(seq_len(n), channel), build))
  if (is.null(m)) matrix(numeric(0), ncol = 3) else m
}

# Map a chosen pair of multiplex channels onto the 2-group analysis_group column
# (1L / 2L) that every cross-channel analysis already consumes. Rows that belong
# to neither chosen channel are dropped. chA/chB are mplex_label strings.
assign_analysis_group_from_mplex <- function(df, chA, chB) {
  if (!"mplex_label" %in% names(df)) return(df[0, , drop = FALSE])
  df$analysis_group <- NA_integer_
  df$analysis_group[df$mplex_label == chA] <- 1L
  df$analysis_group[df$mplex_label == chB] <- 2L
  df[!is.na(df$analysis_group), , drop = FALSE]
}

# Always returns a color vector length nrow(df_pts)
build_colors_for_pts <- function(df_pts, input, default_palette_col = "color") {
  mode <- (input$color_mode %||% "default")
  map_to_gradient <- function(v, pal_name = NULL, reverse = FALSE) {
    v <- as.numeric(v)
    if (!length(v) || all(!is.finite(v))) return(rep("#bdbdbd", length(v)))
    lo <- suppressWarnings(min(v[is.finite(v)], na.rm = TRUE))
    hi <- suppressWarnings(max(v[is.finite(v)], na.rm = TRUE))
    # avoid flat (lo==hi) -> single color
    if (!is.finite(lo) || !is.finite(hi)) return(rep("#bdbdbd", length(v)))
    if (lo == hi) {
      # add tiny epsilon so a scale can still render different bins if needed
      eps <- max(1e-9, abs(lo) * 1e-6)
      v <- v + seq_len(length(v)) * eps
      lo <- min(v, na.rm = TRUE); hi <- max(v, na.rm = TRUE)
    }
    pal <- get_palette_hex(pal_name %||% input$numeric_palette, n = 256, reverse = isTRUE(reverse))
    idx <- round(1 + 255 * (v - lo) / (hi - lo))
    idx[!is.finite(idx)] <- 128
    pal[pmax(1, pmin(256, idx))]
  }
  # -----------------------------
  
  if (mode == "default") {
    # Respect existing categorical colors (trace/cluster) you already computed
    col <- df_pts[[default_palette_col]]
    if (is.null(col)) col <- rep(input$unicolor_hex %||% "#1f77b4", nrow(df_pts))
    col[!nzchar(col)] <- "#bdbdbd"
    col[is.na(col)]   <- "#bdbdbd"
    return(col)
  }
  
  if (mode == "unicolor") {
    return(rep(input$unicolor_hex %||% "#1f77b4", nrow(df_pts)))
  }

  if (mode == "multiplex") {
    col <- df_pts[["mplex_color"]]
    if (is.null(col)) return(rep("#bdbdbd", nrow(df_pts)))
    col[is.na(col) | !nzchar(col)] <- "#bdbdbd"
    return(col)
  }

  # numeric modes: pick the right vector depending on raw vs median
  if (mode == "z") {
    # raw points: 'z' exists; medians: we put 'z_med' in builders
    base <- if ("z" %in% names(df_pts)) df_pts$z else df_pts$z_med
    return(map_to_gradient(base, pal_name = input$numeric_palette, reverse = input$numeric_palette_reverse))
    
  } else if (mode == "tim") {
    # medians: 'tim_med'; raw: fallback to per-point 'tim'
    base <- if ("tim_med" %in% names(df_pts)) df_pts$tim_med else df_pts$tim
    return(map_to_gradient(base, pal_name = input$numeric_palette, reverse = input$numeric_palette_reverse))
    
  }  else if (mode %in% c("density","locdens")) {
    # precomputed in builders: 'dens' (counts per TID/cluster for medians; per-trace for raw)
    base <- if ("dens" %in% names(df_pts)) df_pts$dens else rep(1, nrow(df_pts))
    # if constant, nudge to avoid flat coloring
    if (length(unique(base[is.finite(base)])) <= 1L) {
      eps <- max(1e-9, abs(as.numeric(base[which(is.finite(base))[1]]) %||% 1) * 1e-6)
      base <- as.numeric(base) + seq_len(length(base)) * eps
    }
    return(map_to_gradient(base, pal_name = input$numeric_palette, reverse = input$numeric_palette_reverse))
  }
  
  # fallback
  rep("#bdbdbd", nrow(df_pts))
}
# --- compute_spatial_density ---
compute_spatial_density <- function(x, y, z = NULL, r = 150, k_max = 1000) {
  # Returns integer vector dens, same length as x
  stopifnot(length(x) == length(y))
  n <- length(x)
  if (n == 0) return(integer(0))
  if (n == 1) return(0L)
  
  pts <- if (is.null(z)) cbind(x, y) else cbind(x, y, z)
  
  # frNN returns neighbor lists within eps=r (includes self)
  nn <- dbscan::frNN(pts, eps = r, sort = FALSE)
  
  # density = neighbors within r EXCLUDING self
  dens <- lengths(nn$id) - 1L
  dens[dens < 0] <- 0L
  as.integer(dens)
}
# --- build_pts_2d ---
build_pts_2d <- function(df, input) {
  stopifnot(all(c("x","y","tid") %in% names(df)))
  
  # --- 1) compute per-tid medians once ---
  tid_meds <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z_med   = if ("z"   %in% names(df)) stats::median(z,   na.rm = TRUE) else NA_real_,
      tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
      color = dplyr::first(color),
      hover_text = dplyr::first(hover_text),
      mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
      mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
      mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(x), is.finite(y))

  # --- 2) spatial density of trace medians (2D) ---
  r_nm <- input$dens_radius_nm %||% 150  # add UI input (see below)
  tid_meds$dens <- compute_spatial_density(tid_meds$x, tid_meds$y, r = r_nm)
  
  if (isTRUE(input$show_medians)) {
    
    # If clustering AND you are showing cluster medians:
    # density should probably be computed on *cluster medians* (not tid medians),
    # because that's what you're plotting. We'll keep that consistent.
    if (isTRUE(input$enable_clustering)) {
      pts <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z_med   = if ("z"   %in% names(df)) stats::median(z,   na.rm = TRUE) else NA_real_,
          tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
          color  = dplyr::first(color),
          hover_text = dplyr::first(hover_text),
          mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
          .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(x), is.finite(y))

      pts$dens <- compute_spatial_density(pts$x, pts$y, r = r_nm)
      pts <- pts %>% dplyr::mutate(group_id = cluster)
      
    } else {
      # Plot the tid medians directly (already has dens)
      pts <- tid_meds %>% dplyr::mutate(group_id = tid)
    }
    
  } else {
    # RAW points: join tid-level spatial density onto every localization
    pts <- df %>%
      dplyr::left_join(tid_meds %>% dplyr::select(tid, dens), by = "tid") %>%
      dplyr::mutate(
        z_med   = if ("z"   %in% names(.)) z   else NA_real_,
        tim_med = if ("tim" %in% names(.)) tim else NA_real_,
        group_id = if (isTRUE(input$enable_clustering) && "cluster" %in% names(.))
          ifelse(is.na(cluster) | cluster == 0, -1L, cluster) else tid
      )
  }
  
  # safety columns
  for (nm in c("x","y","z_med","tim_med","dens","color","hover_text","group_id",
               "mplex_color","mplex_channel","mplex_label")) {
    if (!nm %in% names(pts)) pts[[nm]] <- NA
  }
  pts %>% dplyr::filter(is.finite(x), is.finite(y))
}


# --- build_pts_3d ---
build_pts_3d <- function(df, input) {
  stopifnot(all(c("x","y","z","tid") %in% names(df)))
  
  tid_meds <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z = stats::median(z, na.rm = TRUE),
      tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
      color = dplyr::first(color),
      hover_text = dplyr::first(hover_text),
      mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
      mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
      mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(x), is.finite(y), is.finite(z))

  r_nm <- input$dens_radius_nm %||% 150
  tid_meds$dens <- compute_spatial_density(tid_meds$x, tid_meds$y, tid_meds$z, r = r_nm)
  
  if (isTRUE(input$show_medians)) {
    if (isTRUE(input$enable_clustering)) {
      pts <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
          color  = dplyr::first(color),
          hover_text = dplyr::first(hover_text),
          mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
          .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z))

      pts$dens <- compute_spatial_density(pts$x, pts$y, pts$z, r = r_nm)
      pts <- pts %>% dplyr::mutate(group_id = cluster)
      
    } else {
      pts <- tid_meds %>% dplyr::mutate(group_id = tid)
    }
    
  } else {
    pts <- df %>%
      dplyr::left_join(tid_meds %>% dplyr::select(tid, dens), by = "tid") %>%
      dplyr::mutate(
        tim_med = if ("tim" %in% names(.)) tim else NA_real_,
        group_id = if (isTRUE(input$enable_clustering) && "cluster" %in% names(.))
          ifelse(is.na(cluster) | cluster == 0, -1L, cluster) else tid
      )
  }
  
  for (nm in c("x","y","z","tim_med","dens","color","hover_text","group_id",
               "mplex_color","mplex_channel","mplex_label")) {
    if (!nm %in% names(pts)) pts[[nm]] <- NA
  }
  pts %>% dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
}
# --- color_mapping_for_pts ---
color_mapping_for_pts <- function(pts, input) {
  mode     <- (input$color_mode %||% "default")
  pal_name <- input$numeric_palette %||% "Viridis"
  pal_rev  <- isTRUE(input$numeric_palette_reverse)
  
  # Helper: build a numeric colorscale mapping
  build_num_map <- function(v, title) {
    v <- as.numeric(v)
    # All-NA: still return a colorscale shell so Plotly doesn’t choke
    if (!any(is.finite(v))) {
      return(list(
        showscale  = FALSE,
        kind       = "hex",                # fall back to hex (no scale)
        hex        = rep("#bdbdbd", nrow(pts))
      ))
    }
    v_ok <- v[is.finite(v)]
    cmin <- min(v_ok); cmax <- max(v_ok)
    if (!is.finite(cmin) || !is.finite(cmax) || cmin == cmax) {
      # No spread → solid color (pick mid color from scale to be deterministic)
      mid_cols <- plotly_colorscale(pal_name, reverse = pal_rev)
      mid_hex  <- mid_cols[[ceiling(length(mid_cols)/2)]][[2]]
      return(list(
        showscale  = FALSE,
        kind       = "hex",
        hex        = rep(mid_hex, nrow(pts))
      ))
    }
    list(
      showscale  = TRUE,
      kind       = "scale",
      colorvar   = v,
      cmin       = cmin,
      cmax       = cmax,
      colorscale = plotly_colorscale(pal_name, reverse = pal_rev),
      title      = title
    )
  }
  
  if (mode == "default") {
    col <- pts$color
    if (is.null(col)) col <- rep(input$unicolor_hex %||% "#1f77b4", nrow(pts))
    col[!nzchar(col)] <- "#bdbdbd"; col[is.na(col)] <- "#bdbdbd"
    return(list(showscale = FALSE, kind = "hex", hex = col))
  }
  
  if (mode == "unicolor") {
    return(list(
      showscale = FALSE, kind = "hex",
      hex = rep(input$unicolor_hex %||% "#1f77b4", nrow(pts))
    ))
  }

  if (mode == "multiplex") {
    col <- if ("mplex_color" %in% names(pts)) pts$mplex_color else pts$color
    if (is.null(col)) col <- rep("#bdbdbd", nrow(pts))
    col[is.na(col) | !nzchar(col)] <- "#bdbdbd"
    return(list(showscale = FALSE, kind = "hex", hex = col))
  }

  if (mode == "z") {
    zvec <- if ("z" %in% names(pts)) pts$z else if ("z_med" %in% names(pts)) pts$z_med else NA_real_
    return(build_num_map(zvec, "Depth (z)"))
  }
  if (mode == "tim")     return(build_num_map(if ("tim_med" %in% names(pts)) pts$tim_med else pts$tim, "Time"))
  if (mode %in% c("density","locdens")) {
    return(build_num_map(
      if ("dens" %in% names(pts)) pts$dens else rep(1, nrow(pts)),
      "Local trace-median density"
    ))
  }
  
  # Fallback
  list(showscale = FALSE, kind = "hex", hex = rep("#bdbdbd", nrow(pts)))
}
fit_poly_surface <- function(df, degree = 2, grid_n = 50) {
  stopifnot(all(c("x","y","z") %in% names(df)))
  df <- df %>% dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
  if (nrow(df) < (degree + 1L) * (degree + 2L) / 2L) return(NULL) # not enough points
  
  # build formula like: z ~ poly(x,deg,raw=TRUE) + poly(y,deg,raw=TRUE) + cross terms
  # easiest in base R via manual term expansion:
  terms <- c()
  for (i in 0:degree) {
    for (j in 0:(degree - i)) {
      if (i == 0 && j == 0) next
      ti <- if (i == 0) "1" else paste0("I(x^", i, ")")
      tj <- if (j == 0) "1" else paste0("I(y^", j, ")")
      t  <- if (i == 0) tj else if (j == 0) ti else paste0(ti, ":", tj)
      terms <- c(terms, t)
    }
  }
  fml <- stats::as.formula(paste("z ~", paste(terms, collapse = " + ")))
  
  fit <- stats::lm(fml, data = df)
  
  xr <- range(df$x, na.rm = TRUE)
  yr <- range(df$y, na.rm = TRUE)
  
  gx <- seq(xr[1], xr[2], length.out = grid_n)
  gy <- seq(yr[1], yr[2], length.out = grid_n)
  grid <- expand.grid(x = gx, y = gy)
  grid$z <- as.numeric(stats::predict(fit, newdata = grid))
  
  list(
    fit = fit,
    x = matrix(grid$x, nrow = grid_n, ncol = grid_n, byrow = FALSE),
    y = matrix(grid$y, nrow = grid_n, ncol = grid_n, byrow = FALSE),
    z = matrix(grid$z, nrow = grid_n, ncol = grid_n, byrow = FALSE)
  )
}


# --- sd summary and UI update helpers ---
safe_update_sd_range <- function(session, inputId, value) {
  shinyWidgets::updateNumericRangeInput(session, inputId, value = value)
}
# median summary, built on the table above
compute_sd_medians <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT) {
  tbl <- trace_spread_table(df, min_n = min_n)
  if (!nrow(tbl)) {
    # Return a consistent empty result so the callers don't crash
    return(list(x = NA_real_, y = NA_real_, z = if ("sd_z" %in% names(tbl)) NA_real_ else NULL))
  }
  res <- list(
    x = stats::median(tbl$sd_x, na.rm = TRUE),
    y = stats::median(tbl$sd_y, na.rm = TRUE)
  )
  if ("sd_z" %in% names(tbl)) res$z <- stats::median(tbl$sd_z, na.rm = TRUE)
  res
}
trace_spread_table <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT) {
  stopifnot(is.data.frame(df))
  if (!nrow(df)) {
    return(data.frame(
      tid = character(),
      sd_x = numeric(), sd_y = numeric(), sd_z = numeric(),
      spread_rmsd = numeric(), sd_max = numeric()
    ))
  }

  min_n <- max(2L, as.integer(min_n %||% MIN_LOCS_FOR_SPREAD_DEFAULT))

  has_z <- "z" %in% names(df) && any(is.finite(df$z))
  base  <- df %>% dplyr::group_by(tid)

  if (has_z) {
    out <- base %>% dplyr::summarise(
      n_locs = dplyr::n(),
      mx = stats::median(x, na.rm = TRUE),
      my = stats::median(y, na.rm = TRUE),
      mz = stats::median(z, na.rm = TRUE),
      sd_x = stats::sd(x, na.rm = TRUE),
      sd_y = stats::sd(y, na.rm = TRUE),
      sd_z = stats::sd(z, na.rm = TRUE),
      spread_rmsd = sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  } else {
    out <- base %>% dplyr::summarise(
      n_locs = dplyr::n(),
      mx = stats::median(x, na.rm = TRUE),
      my = stats::median(y, na.rm = TRUE),
      sd_x = stats::sd(x, na.rm = TRUE),
      sd_y = stats::sd(y, na.rm = TRUE),
      sd_z = NA_real_,
      spread_rmsd = sqrt(mean((x - mx)^2 + (y - my)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  }
  out %>%
    dplyr::filter(.data$n_locs >= min_n) %>%
    dplyr::mutate(sd_max = pmax(sd_x, sd_y, sd_z, na.rm = TRUE)) %>%
    dplyr::select(-dplyr::any_of(c("mx", "my", "mz", "n_locs")))
}
# Generate histogram breaks for a numeric vector over given number of bins
sd_breaks <- function(x, bins) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    rng <- c(rng[1] - 0.5, rng[2] + 0.5)
  }
  seq(rng[1], rng[2], length.out = max(2, as.integer(bins) + 1L))
}

# Turn a numeric vector into (centers, counts) on given breaks
sd_counts <- function(x, breaks) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (is.null(breaks)) return(list(centers = numeric(0), counts = numeric(0), width = 1))
  
  bw  <- diff(breaks)
  eps <- if (length(bw)) max(1e-9, 1e-6 * max(bw, na.rm = TRUE)) else 1e-6
  brk <- breaks
  brk[1] <- brk[1] - eps
  brk[length(brk)] <- brk[length(brk)] + eps
  
  if (length(x)) {
    lo <- brk[1]; hi <- brk[length(brk)]
    x  <- pmin(hi, pmax(lo, x))
  }
  
  h <- hist(x, breaks = brk, plot = FALSE, include.lowest = TRUE, right = TRUE)
  
  w <- suppressWarnings(as.numeric(diff(h$breaks)[1]))
  if (!is.finite(w) || w <= 0) w <- 1
  
  list(centers = h$mids, counts = h$counts, width = w)
}
# Generate histogram breaks from 0 to max value in raw_vals over given number of bins
sd_breaks_0max_from_raw <- function(raw_vals, bins) {
  v <- as.numeric(raw_vals)
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  
  mx <- max(v, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) mx <- 1
  seq(0, mx, length.out = max(2, as.integer(bins) + 1L))
}
sd_by_tid <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT) {
  need_cols <- intersect(c("tid","x","y","z"), names(df))
  if (!all(c("tid","x","y") %in% need_cols)) return(data.frame())

  min_n <- max(2L, as.integer(min_n %||% MIN_LOCS_FOR_SPREAD_DEFAULT))

  d <- df[, need_cols, drop = FALSE]
  d$x <- as.numeric(d$x); d$y <- as.numeric(d$y)
  has_z <- "z" %in% names(d)
  if (has_z) d$z <- as.numeric(d$z)

  d %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      n_locs = dplyr::n(),
      sd_x = if (dplyr::n() >= min_n) stats::sd(x, na.rm = TRUE) else NA_real_,
      sd_y = if (dplyr::n() >= min_n) stats::sd(y, na.rm = TRUE) else NA_real_,
      sd_z = if (has_z && dplyr::n() >= min_n) stats::sd(z, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}
# =============================================================================
# cov_ellipse_params() — THE ERROR ELLIPSE OF A 2x2 COVARIANCE MATRIX
# =============================================================================
#
# GOAL
#   Turn a covariance matrix into a drawable, ROTATED ellipse plus the two
#   numbers that describe the cloud's shape independently of the axes.
#
# WHY ROTATION MATTERS (the thing σx / σy alone cannot express)
#   An axis-aligned ellipse has two parameters and cannot represent a TILTED
#   cloud. Elongation along a diagonal inflates σx and σy by the SAME amount, so
#   σx ≈ σy and the ellipse renders as a circle over a visibly stretched cloud —
#   all the elongation lives in cov(x,y), which an axis-aligned drawing discards.
#
#   The same blind spot afflicts max(σx,σy)/min(σx,σy) as an anisotropy metric:
#   a 45° cloud scores ~1.0, i.e. "isotropic". Since drift within a trace is
#   DIRECTIONAL, a tilted cloud is the expected signature of drift — precisely
#   the artifact the precision panel exists to catch.
#
#   The eigen-decomposition is the fix: eigenvectors give the true major/minor
#   directions, sqrt(eigenvalues) the semi-axes. Their ratio is rotation-
#   invariant, so it detects elongation at ANY angle.
#
# ARGUMENTS  S : symmetric 2x2 covariance matrix;  k : how many σ to scale to
# RETURNS
#   list(ra, rb, theta, ratio, rho, angle_deg) — ra/rb are semi-major/minor in
#   the data's units, theta the major-axis angle in radians. NULL if unusable.
# =============================================================================
cov_ellipse_params <- function(S, k = 1) {
  if (!is.matrix(S) || any(!is.finite(S))) return(NULL)
  ev <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev)) return(NULL)
  vals <- pmax(ev$values, 0)                      # clamp tiny negative round-off
  if (!all(is.finite(vals)) || vals[1] <= 0) return(NULL)
  th <- atan2(ev$vectors[2, 1], ev$vectors[1, 1]) # direction of the major axis
  den <- sqrt(max(S[1, 1] * S[2, 2], 1e-12))
  list(
    ra    = k * sqrt(vals[1]),
    rb    = k * sqrt(vals[2]),
    theta = th,
    ratio = sqrt(vals[1] / max(vals[2], 1e-12)),  # rotation-invariant anisotropy
    rho   = S[1, 2] / den,                        # correlation, for reference
    angle_deg = ((th * 180 / pi) + 180) %% 180    # fold to [0,180): an axis has no sign
  )
}

# GOAL: an SVG path tracing a rotated ellipse, for plotly shapes.
# plotly's type = "circle" shape can only ever be AXIS-ALIGNED (it is defined by
# a bounding box), so a rotated ellipse has to be drawn as an explicit path.
ellipse_svg_path <- function(ra, rb, theta, cx = 0, cy = 0, n = 181L) {
  if (!all(is.finite(c(ra, rb, theta, cx, cy))) || ra <= 0 || rb <= 0) return(NULL)
  t  <- seq(0, 2 * pi, length.out = n)
  ct <- cos(theta); st <- sin(theta)
  x  <- cx + ra * cos(t) * ct - rb * sin(t) * st
  y  <- cy + ra * cos(t) * st + rb * sin(t) * ct
  paste0("M ", paste(sprintf("%.4f,%.4f", x, y), collapse = " L "), " Z")
}

# =============================================================================
# density_hdr_2d() — SMOOTH DENSITY + HIGHEST-DENSITY-REGION CONTOURS
# =============================================================================
#
# GOAL
#   Estimate the 2D density of a point cloud on a fixed frame and return both
#   (a) a grid smooth enough to draw as a continuous heatmap, and (b) contour
#   lines at chosen PROBABILITY MASSES rather than at arbitrary density levels.
#
# WHY NOT A BINNED HISTOGRAM-CONTOUR (what the precision cloud used to draw)
#   A 2D histogram estimates each bin independently, so every bin carries full
#   Poisson noise. For the precision cloud (~40k points, 60x60 bins) that is
#   ~9% noise at the peak but ~26% at the 2-sigma radius and ~90% at 3 sigma —
#   i.e. exactly the outer contours the eye reads as "the shape of the cloud"
#   are almost pure shot noise, which is what makes the plot look ragged and
#   under-smoothed. Re-binning cannot fix it: coarser bins blur the core,
#   finer bins make the tails worse. A kernel estimate is the fix, because each
#   point contributes to a neighbourhood and the estimate is continuous.
#
# WHY HDR LEVELS
#   A contour at "density = 0.004" means nothing on its own. A contour holding
#   39.3% / 86.5% of the points means a lot, because those are EXACTLY the
#   masses a Gaussian 1-sigma and 2-sigma ellipse contain in 2D. Drawn together
#   with the fitted ellipses, the empirical contour and the theoretical ring
#   describe the same probability mass, so they coincide when the cloud really
#   is Gaussian and separate visibly where it is not. That comparison is the
#   whole diagnostic, and it is unavailable with count-based levels.
#
# BANDWIDTH
#   Scott's rule for 2D (sigma * n^(-1/6)) with a floor at the grid pitch —
#   below the pitch a kernel cannot smooth anything, it just re-creates the
#   binning artefact. sigma is passed in (use a ROBUST sigma) rather than taken
#   from the data, so a fiducial tail cannot inflate the bandwidth.
#
# ARGUMENTS
#   a, b    : coordinates (same units)
#   lim     : frame half-width; the grid spans [-lim, lim] on both axes
#   sigma   : scale for the bandwidth rule (robust sigma of the cloud)
#   n_grid  : grid points per axis
#   probs   : probability masses to contour
#   max_pts : subsample cap for the kernel sum (cost is O(n * n_grid))
# RETURNS
#   list(x, y, z, levels, paths) where paths is a list of
#   list(prob=, d=, path=) SVG paths ready to use as plotly layout shapes, or
#   NULL if the cloud is unusable. A requested mass that extends beyond the
#   frame is skipped rather than drawn wrong.
# =============================================================================
density_hdr_2d <- function(a, b, lim, sigma, n_grid = 160L,
                           probs = c(0.393, 0.865), max_pts = 20000L) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 20L || !is.finite(lim) || lim <= 0) return(NULL)
  if (length(a) > max_pts) {                       # kernel sum is O(n * n_grid)
    s <- sort(sample.int(length(a), max_pts)); a <- a[s]; b <- b[s]
  }
  pitch <- 2 * lim / (n_grid - 1)                  # grid spacing
  if (!is.finite(sigma) || sigma <= 0) sigma <- stats::mad(sqrt(a^2 + b^2))
  if (!is.finite(sigma) || sigma <= 0) return(NULL)
  # Scott's rule, floored at the grid pitch so the kernel always spans >1 cell.
  bw_sd <- max(sigma * length(a)^(-1 / 6), pitch)
  # MASS::kde2d's `h` is 4x the kernel SD (it uses dnorm(sd = h/4) internally).
  k <- tryCatch(MASS::kde2d(a, b, h = c(4 * bw_sd, 4 * bw_sd), n = n_grid,
                            lims = c(-lim, lim, -lim, lim)),
                error = function(e) NULL)
  if (is.null(k) || !all(is.finite(k$z))) return(NULL)

  # HDR level for mass p: walk the grid cells from densest down, accumulating
  # probability (density x cell area) until p is reached. The density there is
  # the level whose contour encloses exactly that mass.
  cell <- (k$x[2] - k$x[1]) * (k$y[2] - k$y[1])
  srt  <- sort(as.vector(k$z), decreasing = TRUE)
  cm   <- cumsum(srt) * cell                       # mass enclosed at each level
  in_frame <- cm[length(cm)]                       # total mass inside the frame
  paths <- list()
  for (p in probs) {
    # If the frame does not even contain p of the cloud, no honest contour for p
    # exists here — skip it rather than draw the frame edge and imply otherwise.
    if (!is.finite(in_frame) || in_frame < p) next
    idx <- which(cm >= p)[1]
    lev <- srt[idx]
    if (!is.finite(lev) || lev <= 0) next
    cls <- tryCatch(grDevices::contourLines(k$x, k$y, k$z, levels = lev),
                    error = function(e) NULL)
    for (cl in cls) {
      if (length(cl$x) < 3L) next
      paths[[length(paths) + 1L]] <- list(
        prob = p, d = lev,
        path = paste0("M ", paste(sprintf("%.4f,%.4f", cl$x, cl$y), collapse = " L "), " Z"))
    }
  }
  list(x = k$x, y = k$y, z = k$z, bw_sd = bw_sd, mass_in_frame = in_frame, paths = paths)
}

# =============================================================================
# robust_median_cov() — ONE VOTE PER TRACE, BUT STILL TILT-AWARE
# =============================================================================
# GOAL
#   Combine many per-trace covariance matrices into one robust matrix by taking
#   the element-wise median, so a single 2000-localization fiducial cannot
#   dominate the shape the way it dominates a pooled estimate.
#
#   Element-wise medians are not guaranteed to stay positive semi-definite, so
#   the off-diagonal is shrunk if it would exceed what the diagonal supports —
#   without that, eigen() can return a negative eigenvalue and the ellipse
#   collapses or inverts.
#
# ARGUMENTS  vxx, vyy : per-trace variances;  vxy : per-trace covariances
# RETURNS    2x2 matrix, or NULL when there is nothing usable.
# =============================================================================
robust_median_cov <- function(vxx, vyy, vxy) {
  m <- function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) }
  a <- m(vxx); b <- m(vyy); c0 <- m(vxy)
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(NULL)
  if (!is.finite(c0)) c0 <- 0
  lim <- 0.99 * sqrt(a * b)                        # keep |rho| < 1
  if (abs(c0) > lim) c0 <- sign(c0) * lim
  matrix(c(a, c0, c0, b), 2L, 2L)
}

# =============================================================================

# --- sigma label helpers ---
sigma_label_med_mad <- function(v, show_mad = FALSE, digits = 2) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  
  med <- stats::median(v, na.rm = TRUE)
  
  if (!isTRUE(show_mad)) {
    return(sprintf("\u03C3\u0303 %.2f nm", round(med, digits)))
  }
  
  mad <- stats::mad(v, constant = 1.4826, na.rm = TRUE)  # MAD in nm
  sprintf("\u03C3\u0303 %.2f nm<br>MAD %.2f nm", round(med, digits), round(mad, digits))
}
sigma_label_med_pm_sd <- function(v, digits = 2) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)

  med <- stats::median(v, na.rm = TRUE)
  mad <- stats::mad(v, constant = 1.4826, na.rm = TRUE)  # normalized MAD (= MAD / 0.67449), matches pyMINFLUX

  if (!is.finite(med)) return(NULL)
  if (!is.finite(mad)) mad <- NA_real_

  # σ̃ is "sigma with tilde" meaning median-ish notation
  if (is.finite(mad)) {
    sprintf("\u03C3\u0303 %.2f \u00B1 %.2f nm", round(med, digits), round(mad, digits))
  } else {
    sprintf("\u03C3\u0303 %.2f nm", round(med, digits))
  }
}

add_sigma_vline <- function(fig, x, y_top, color, label, dash = NULL) {
  if (!is.finite(x) || !is.finite(y_top) || y_top <= 0) return(fig)
  
  # Add the line
  fig <- fig %>% plotly::add_segments(
    x = x, xend = x,
    y = 0, yend = y_top,
    line = c(list(color = color, width = 2), if (!is.null(dash)) list(dash = dash) else list()),
    inherit = FALSE,
    hoverinfo = "skip",
    showlegend = FALSE
  )
  
  # Add the sigma label slightly above the line
  fig %>% plotly::layout(
    annotations = c(
      fig$x$layout$annotations %||% list(),
      list(
        x = x, y = y_top,
        xref = "x", yref = "y",
        text = label,
        showarrow = FALSE,
        xanchor = "left",
        yanchor = "bottom",
        font = list(color = color, size = 12),
        bgcolor = "rgba(255,255,255,0.0)"
      )
    )
  )
}


# convenience for plotly fixed-width bins
compute_xbins <- function(v, bins) {
  rng <- range(v[is.finite(v)], na.rm = TRUE)
  if (!all(is.finite(rng)) || rng[1] == rng[2]) {
    return(list(start = 0, end = 1, size = 1))  # harmless fallback
  }
  w <- (rng[2] - rng[1]) / max(1, bins)
  list(start = rng[1], end = rng[2], size = w)
}


median_points_df <- function(df, use_cluster1 = FALSE) {
  if (use_cluster1 && "cluster" %in% names(df)) {
    # cluster medians
    med <- df %>%
      dplyr::filter(!is.na(cluster) & cluster != 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
        tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
        dens    = dplyr::n(),                 # “localization density” proxy
        .groups = "drop"
      ) %>%
      dplyr::rename(group_id = cluster)
  } else {
    # per-trace medians
    med <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
        tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
        dens    = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::rename(group_id = tid)
  }
  med
}
# 2D circle generator from a single cluster's points
# - Aggregates to per-tid medians (x,y)
# - Robust center via geometric median (fallback: mean)
# - Radius = median distance to center
# - Returns perimeter samples in $virtual_points (N×2)
generate_circle_from_cluster2 <- function(cluster_data, n_points = 2000) {
  # Expect: columns tid, x, y   (z ignored/not required)
  if (!all(c("tid", "x", "y") %in% names(cluster_data))) return(NULL)
  
  coords <- cluster_data %>%
    dplyr::group_by(.data$tid) %>%
    dplyr::summarise(
      x = stats::median(.data$x, na.rm = TRUE),
      y = stats::median(.data$y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::select(.data$x, .data$y)
  
  coords_matrix <- as.matrix(coords[, c("x","y")])
  # Need at least 3 non-collinear points for a stable circle estimate
  if (nrow(coords_matrix) < 3) return(NULL)
  
  # Robust center: geometric median if available; otherwise mean
  center <- tryCatch({
    if (requireNamespace("Gmedian", quietly = TRUE)) {
      as.numeric(Gmedian::geometric.median(coords_matrix)$median)
    } else {
      colMeans(coords_matrix, na.rm = TRUE)
    }
  }, error = function(e) colMeans(coords_matrix, na.rm = TRUE))
  center <- as.numeric(center)[1:2]
  
  # Median radius (robust to outliers)
  distances <- sqrt(rowSums((coords_matrix - matrix(center, nrow(coords_matrix), 2, byrow = TRUE))^2))
  radius <- stats::median(distances, na.rm = TRUE)
  
  if (!is.finite(radius) || radius < 1e-6) return(NULL)
  
  # Sample the circle perimeter deterministically (even angular spacing)
  t <- seq(0, 2*pi, length.out = n_points)
  x <- center[1] + radius * cos(t)
  y <- center[2] + radius * sin(t)
  virtual_points <- cbind(x, y)
  colnames(virtual_points) <- c("x", "y")
  
  list(
    center = center,            # c(cx, cy)
    radius = radius,            # scalar
    virtual_points = virtual_points  # N×2 perimeter points
  )
}
# 2D convex-hull centroid (polygon centroid of the hull)
# Input: coords_matrix (N x 2) with columns [x, y]
# Output: numeric length-2 c(cx, cy)
# --- per-trace SD table from currently filtered data ---
# sd_by_tid is defined above (canonical helper).

# --- generic plotly histogram with a median line ---
sd_hist_plot <- function(vals, title = "", dark = FALSE, bins = 40, show_median = TRUE) {
  vals <- as.numeric(vals)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(plotly::plot_ly())
  
  fill_col <- if (dark) "rgba(0, 150, 255, 0.65)" else "rgba(0, 90, 200, 0.65)"
  line_col <- if (dark) "white" else "black"
  
  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      x = vals,
      type = "histogram",
      nbinsx = max(1L, as.integer(bins)),
      marker = list(color = fill_col, line = list(color = line_col, width = 1)),
      hovertemplate = "Count: %{y}<br>Bin: %{x}<extra></extra>"
    ) %>%
    plotly::layout(
      title = title,
      xaxis = list(title = "SD"),
      yaxis = list(title = "Count"),
      paper_bgcolor = if (dark) "#1e1e1e" else "white",
      plot_bgcolor  = if (dark) "#1e1e1e" else "white",
      font = list(color = if (dark) "white" else "black")
    )
  
  if (isTRUE(show_median)) {
    med <- stats::median(vals, na.rm = TRUE)
    p <- p %>%
      plotly::layout(
        shapes = list(list(
          type = "line",
          x0 = med, x1 = med, xref = "x",
          y0 = 0,   y1 = 1,   yref = "paper",
          line = list(color = if (dark) "orange" else "red", width = 2, dash = "dash")
        )),
        annotations = list(list(
          x = med, y = 1, xref = "x", yref = "paper",
          yanchor = "bottom",
          text = sprintf("median = %.2f", med),
          showarrow = FALSE,
          font = list(color = if (dark) "orange" else "red", size = 12)
        ))
      )
  }
  p
}

# palette helper (uses RColorBrewer that you already load)
palette_colors <- function(name, n = 256) {
  get_palette_hex(name, n = n, reverse = FALSE)
}
plotly_colorscale <- function(name = "Viridis", reverse = FALSE, n = 11) {
  cols <- get_palette_hex(name, n = n, reverse = reverse)
  pos  <- seq(0, 1, length.out = length(cols))
  Map(function(p, c) list(p, c), pos, cols)
}


# Helper: Remove unstable traces instead of correcting
remove_unstable_traces <- function(df, input) {
  min_n <- max(2L, as.integer(input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT))

  # Per-trace spread is only defined for traces with >= min_n localizations;
  # shorter traces get NA (undefined spread) so their noisy SD never drives the
  # filter — they are kept rather than judged on an unreliable spread.
  spread_df <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      sd_x = if (dplyr::n() >= min_n) stats::sd(x, na.rm = TRUE) else NA_real_,
      sd_y = if (dplyr::n() >= min_n) stats::sd(y, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  # Compute fallback thresholds from traces with a defined spread only.
  fallback <- spread_df %>%
    dplyr::summarise(
      x = stats::median(sd_x, na.rm = TRUE),
      y = stats::median(sd_y, na.rm = TRUE)
    )

  target_x <- ifelse(is.na(input$target_sd_x), fallback$x, input$target_sd_x)
  target_y <- ifelse(is.na(input$target_sd_y), fallback$y, input$target_sd_y)

  # Keep traces within threshold; traces with undefined spread (NA) are kept.
  valid_tids <- spread_df %>%
    dplyr::filter(dplyr::coalesce(sd_x, 0) <= target_x,
                  dplyr::coalesce(sd_y, 0) <= target_y) %>%
    dplyr::pull(tid)

  df %>% dplyr::filter(tid %in% valid_tids)
}

# ---- Batch NN helpers -------------------------------------------------------

# Apply current metadata filter inputs to a raw stored df (no bypass logic —
# the ranges are applied as-is so that the same thresholds work across all
# loaded files in batch mode).
apply_meta_filters_batch <- function(df, input) {
  hard_filter <- function(dat, nm, rng_val) {
    if (!nm %in% names(dat)) return(dat)
    if (is.null(rng_val) || length(rng_val) != 2 || any(!is.finite(as.numeric(rng_val)))) return(dat)
    lo <- rng_val[1]; hi <- rng_val[2]
    dplyr::filter(dat, is.na(.data[[nm]]) | (.data[[nm]] >= lo & .data[[nm]] <= hi))
  }
  df <- df %>%
    hard_filter("efo",   input$efo_range)   %>%
    hard_filter("cfr",   input$cfr_range)   %>%
    hard_filter("dcr",   input$dcr_range)   %>%
    hard_filter("eco",   input$eco_range)   %>%
    hard_filter("tim",   input$tim_range)   %>%
    hard_filter("dwell", input$dwell_range) %>%
    hard_filter("fluo",  input$fluo_range)
  # TID length filter
  tlr <- input$tid_length_range
  if ("tid" %in% names(df) && !is.null(tlr) && length(tlr) == 2 && all(is.finite(tlr))) {
    tl   <- table(df$tid)
    keep <- names(tl)[tl >= tlr[1] & tl <= tlr[2]]
    df   <- dplyr::filter(df, tid %in% keep)
  }
  # Trace spread (σ) filter — mirrors remove_unstable_traces() from Display tab
  if (isTRUE(input$view_SD %||% FALSE) && "tid" %in% names(df)) {
    tx <- input$target_sd_x; ty <- input$target_sd_y
    if (!is.null(tx) && !is.null(ty) && is.finite(tx) && is.finite(ty) && tx > 0 && ty > 0) {
      min_n <- max(2L, as.integer(input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT))
      spread <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          sd_x = if (dplyr::n() >= min_n) stats::sd(x, na.rm = TRUE) else NA_real_,
          sd_y = if (dplyr::n() >= min_n) stats::sd(y, na.rm = TRUE) else NA_real_,
          .groups = "drop")
      # Traces with undefined spread (NA, too short) are kept, not filtered.
      valid <- spread$tid[dplyr::coalesce(spread$sd_x, 0) <= tx &
                          dplyr::coalesce(spread$sd_y, 0) <= ty]
      df <- dplyr::filter(df, tid %in% valid)
    }
  }
  df
}

# Compute k-NN distances between unit centroids.
# unit:         "tid"           — one centroid per TID (global)
#               "pcid"          — one centroid per primary cluster (global)
#               "cluster2_tid"  — TID centroids within each cluster2
#               "cluster2_pcid" — primary cluster centroids within each cluster2
# group_filter: NA (all), 1L, or 2L to restrict to an analysis_group
# Returns a data.frame with column `distance`.
compute_batch_nn <- function(df, unit = "tid", k = 1L, group_filter = NA_integer_,
                             cross_groups = NULL) {
  if (!all(c("x", "y", "tid") %in% names(df))) return(data.frame())
  if (!"z" %in% names(df)) df$z <- 0

  # Between-group (cross-channel) NN: query from cross_groups[1] to cross_groups[2]
  if (!is.null(cross_groups) && length(cross_groups) == 2L &&
      "analysis_group" %in% names(df)) {
    get_pts_grp <- function(sub) {
      if (unit == "pcid" && "cluster" %in% names(sub) &&
          any(!is.na(sub$cluster) & sub$cluster > 0L)) {
        pts <- sub %>%
          dplyr::filter(!is.na(cluster), cluster > 0L) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                           y = stats::median(y, na.rm = TRUE),
                           z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
          dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
      } else {
        pts <- sub %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                           y = stats::median(y, na.rm = TRUE),
                           z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
          dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
      }
      pts
    }
    df_a <- df[!is.na(df$analysis_group) & df$analysis_group == cross_groups[1L], ]
    df_b <- df[!is.na(df$analysis_group) & df$analysis_group == cross_groups[2L], ]
    pts_a <- get_pts_grp(df_a); pts_b <- get_pts_grp(df_b)
    if (nrow(pts_a) < 1L || nrow(pts_b) < 1L) return(data.frame())
    k_eff <- min(k, nrow(pts_b))
    nn <- RANN::nn2(as.matrix(pts_b[, c("x","y","z")]),
                    as.matrix(pts_a[, c("x","y","z")]), k = k_eff)
    d <- as.vector(nn$nn.dists[, k_eff])
    return(data.frame(distance = d[is.finite(d) & d > 1e-9]))
  }

  if (!is.na(group_filter) && "analysis_group" %in% names(df))
    df <- df[!is.na(df$analysis_group) & df$analysis_group == as.integer(group_filter), , drop = FALSE]

  run_nn <- function(coords, k_eff) {
    k_eff <- min(k_eff, nrow(coords) - 1L)
    if (k_eff < 1L) return(numeric(0))
    nn <- tryCatch(RANN::nn2(coords, coords, k = k_eff + 1L), error = function(e) NULL)
    if (is.null(nn)) return(numeric(0))
    d <- as.vector(nn$nn.dists[, seq_len(k_eff) + 1L, drop = FALSE])
    d[is.finite(d) & d > 1e-9]
  }

  if (unit == "tid") {
    pts <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(pts) < 2L) return(data.frame())
    return(data.frame(distance = run_nn(as.matrix(pts[, c("x","y","z")]), k)))

  } else if (unit == "pcid") {
    if (!"cluster" %in% names(df) || all(is.na(df$cluster) | df$cluster == 0L))
      return(data.frame())
    pts <- df %>%
      dplyr::filter(!is.na(cluster), cluster > 0L) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(pts) < 2L) return(data.frame())
    return(data.frame(distance = run_nn(as.matrix(pts[, c("x","y","z")]), k)))

  } else if (unit == "cluster2_tid") {
    if (!"cluster2" %in% names(df) || all(is.na(df$cluster2) | df$cluster2 == 0L))
      return(data.frame())
    tid_pts <- df %>%
      dplyr::filter(!is.na(cluster2), cluster2 > 0L) %>%
      dplyr::group_by(cluster2, tid) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(tid_pts) < 2L) return(data.frame())
    rows <- lapply(split(tid_pts, tid_pts$cluster2), function(g) {
      d <- run_nn(as.matrix(g[, c("x","y","z")]), k)
      if (!length(d)) return(NULL)
      data.frame(distance = d)
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    return(if (is.null(out)) data.frame() else out)

  } else if (unit == "cluster2_pcid") {
    if (!all(c("cluster", "cluster2") %in% names(df))) return(data.frame())
    if (all(is.na(df$cluster2) | df$cluster2 == 0L)) return(data.frame())
    if (all(is.na(df$cluster)  | df$cluster  == 0L)) return(data.frame())
    pcid_pts <- df %>%
      dplyr::filter(!is.na(cluster2), cluster2 > 0L,
                    !is.na(cluster),  cluster  > 0L) %>%
      dplyr::group_by(cluster2, cluster) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(pcid_pts) < 2L) return(data.frame())
    rows <- lapply(split(pcid_pts, pcid_pts$cluster2), function(g) {
      d <- run_nn(as.matrix(g[, c("x","y","z")]), k)
      if (!length(d)) return(NULL)
      data.frame(distance = d)
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    return(if (is.null(out)) data.frame() else out)
  }

  data.frame()
}

# ---- 2D aliases
assign_analysis_group_2d            <- assign_analysis_group
get_equal_limits_2d_fn              <- get_equal_limits_2d


#3D Helpers


get_equal_limits <- function(x, y, z = NULL, padding = 0.05) {
  xr <- range(x, na.rm = TRUE); if (!is.finite(diff(xr)) || diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  yr <- range(y, na.rm = TRUE); if (!is.finite(diff(yr)) || diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
  
  if (is.null(z)) {
    span <- max(diff(xr), diff(yr))
    pad  <- span * padding
    xm <- mean(xr); ym <- mean(yr)
    return(list(
      x = c(xm - span/2 - pad, xm + span/2 + pad),
      y = c(ym - span/2 - pad, ym + span/2 + pad)
    ))
  } else {
    zr <- range(z, na.rm = TRUE); if (!is.finite(diff(zr)) || diff(zr) == 0) zr <- zr + c(-0.5, 0.5)
    span <- max(diff(xr), diff(yr), diff(zr))
    pad  <- span * padding
    xm <- mean(xr); ym <- mean(yr); zm <- mean(zr)
    return(list(
      x = c(xm - span/2 - pad, xm + span/2 + pad),
      y = c(ym - span/2 - pad, ym + span/2 + pad),
      z = c(zm - span/2 - pad, zm + span/2 + pad)
    ))
  }
}

safe_which_max <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  which.max(replace(x, is.na(x), -Inf))
}

assign_analysis_group_legacy <- function(df, input) {
  have_dcr  <- "dcr"  %in% names(df)
  have_fluo <- "fluo" %in% names(df)
  
  med_tbl <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      dcr_med  = if (have_dcr)  stats::median(dcr,  na.rm = TRUE) else NA_real_,
      fluo_med = if (have_fluo) stats::median(fluo, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  df <- dplyr::left_join(df, med_tbl, by = "tid")
  mode <- input$grouping_mode %||% "none"
  
  if (mode == "fluo") {
    if (!have_fluo) return(dplyr::mutate(df, analysis_group = NA_integer_))
    fm <- suppressWarnings(as.integer(round(df$fluo_med)))
    df$analysis_group <- dplyr::case_when(
      is.finite(fm) & fm == 1L ~ 1L,
      is.finite(fm) & fm == 2L ~ 2L,
      TRUE                     ~ NA_integer_
    )
    
  } else if (mode == "dcr") {
    v <- df$dcr_med
    v <- as.numeric(v)
    
    g1 <- input$dcr_group1; g2 <- input$dcr_group2
    ok <- have_dcr && !is.null(g1) && !is.null(g2) &&
      length(g1) == 2 && length(g2) == 2 && all(is.finite(c(g1, g2)))
    
    df$analysis_group <- if (ok) {
      dplyr::case_when(
        is.finite(v) & v >= g1[1] & v <  g1[2] ~ 1L,
        is.finite(v) & v >= g2[1] & v <= g2[2] ~ 2L,
        TRUE                                   ~ NA_integer_
      )
    } else NA_integer_
    
  } else {
    df$analysis_group <- NA_integer_
  }
  df
}

trace_spread_table_legacy <- function(df) {
  df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      mx = stats::median(x, na.rm = TRUE),
      my = stats::median(y, na.rm = TRUE),
      mz = stats::median(z, na.rm = TRUE),
      # RMSD in 3D to each trace median (robust overall spread)
      spread_rmsd = sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2, na.rm = TRUE)),
      # Per-axis SDs (catch “needle” traces along one axis)
      sd_x = stats::sd(x, na.rm = TRUE),
      sd_y = stats::sd(y, na.rm = TRUE),
      sd_z = stats::sd(z, na.rm = TRUE),
      n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sd_max = pmax(sd_x, sd_y, sd_z, na.rm = TRUE)
    )
}

# Decide which tids are rogue
flag_rogue_tids <- function(spread_df, rule = "factor_median", factor = 2.5) {
  factor <- as.numeric(factor %||% 2.5)
  factor <- if (is.finite(factor) && factor > 0) factor else 2.5
  
  if (rule == "axis_sd") {
    med_sdmax <- stats::median(spread_df$sd_max, na.rm = TRUE)
    cutoff    <- factor * med_sdmax
    bad <- spread_df$tid[spread_df$sd_max > cutoff]
    attr(bad, "cutoff_desc") <- sprintf("sd_max > %.2f × median(sd_max) = %.3f", factor, cutoff)
    return(bad)
  }
  
  if (rule == "mad_rule") {
    m   <- stats::median(spread_df$spread_rmsd, na.rm = TRUE)
    mad <- stats::mad(spread_df$spread_rmsd, constant = 1.4826, na.rm = TRUE)
    cutoff <- m + factor * mad
    bad <- spread_df$tid[spread_df$spread_rmsd > cutoff]
    attr(bad, "cutoff_desc") <- sprintf("RMSD > median + %.2f × MAD = %.3f", factor, cutoff)
    return(bad)
  }
  
  # default: factor_median
  med <- stats::median(spread_df$spread_rmsd, na.rm = TRUE)
  cutoff <- factor * med
  bad <- spread_df$tid[spread_df$spread_rmsd > cutoff]
  attr(bad, "cutoff_desc") <- sprintf("RMSD > %.2f × median = %.3f", factor, cutoff)
  bad
}

remove_rogue_traces <- function(df, rule = c("factor_median","mad_rule","axis_sd"), factor = 2.5, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT) {
  rule <- match.arg(rule)
  
  spread <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      sd_x = stats::sd(x, na.rm = TRUE),
      sd_y = stats::sd(y, na.rm = TRUE),
      sd_z = if ("z" %in% names(df)) stats::sd(z, na.rm = TRUE) else NA_real_,
      n    = dplyr::n(),
      .groups = "drop"
    )
  
  keep <- rep(TRUE, nrow(spread))
  note <- ""
  
  if (rule == "factor_median") {
    med_all <- c(stats::median(spread$sd_x, na.rm = TRUE),
                 stats::median(spread$sd_y, na.rm = TRUE),
                 if ("z" %in% names(df)) stats::median(spread$sd_z, na.rm = TRUE) else NA)
    thr_x <- factor * med_all[1]
    thr_y <- factor * med_all[2]
    thr_z <- if (!is.na(med_all[3])) factor * med_all[3] else Inf
    keep <- (spread$sd_x <= thr_x) & (spread$sd_y <= thr_y) & (spread$sd_z <= thr_z | is.na(spread$sd_z))
    note <- sprintf("factor_median (×%.2f)", factor)
    
  } else if (rule == "mad_rule") {
    mx <- stats::median(spread$sd_x, na.rm = TRUE); madx <- stats::mad(spread$sd_x, na.rm = TRUE)
    my <- stats::median(spread$sd_y, na.rm = TRUE); mady <- stats::mad(spread$sd_y, na.rm = TRUE)
    if ("z" %in% names(df)) { mz <- stats::median(spread$sd_z, na.rm = TRUE); madz <- stats::mad(spread$sd_z, na.rm = TRUE) } else { mz <- NA; madz <- NA }
    thr_x <- mx + factor * madx; thr_y <- my + factor * mady; thr_z <- if (!is.na(mz)) mz + factor * madz else Inf
    keep <- (spread$sd_x <= thr_x) & (spread$sd_y <= thr_y) & (spread$sd_z <= thr_z | is.na(spread$sd_z))
    note <- sprintf("mad_rule (median + %.2f × MAD)", factor)
    
  } else if (rule == "axis_sd") {
    medx <- stats::median(spread$sd_x, na.rm = TRUE)
    medy <- stats::median(spread$sd_y, na.rm = TRUE)
    medz <- if ("z" %in% names(df)) stats::median(spread$sd_z, na.rm = TRUE) else NA
    thr_x <- factor * medx; thr_y <- factor * medy; thr_z <- if (!is.na(medz)) factor * medz else Inf
    keep <- (spread$sd_x <= thr_x) & (spread$sd_y <= thr_y) & (spread$sd_z <= thr_z | is.na(spread$sd_z))
    note <- sprintf("axis_sd (any axis > %.2f × median(axis SD))", factor)
  }
  
  keep <- keep & (spread$n >= min_n)
  kept_tids <- spread$tid[keep]
  removed_tids <- spread$tid[!keep]
  
  list(
    df = df %>% dplyr::filter(tid %in% kept_tids),
    removed_tids = removed_tids,
    removed_n = length(removed_tids),
    note = note
  )
}





calculate_distance_to_fit <- function(df           = NULL,
                                      coords_matrix = NULL,
                                      fit_type     = c("sphere", "ellipsoid", "helix", "ring", "arc", "best"),
                                      cluster2_id  = NULL,
                                      use_cluster1 = FALSE,
                                      turns        = 3,
                                      helix_n      = 5000,
                                      arc_span_hint = 0) {

  fit_type <- match.arg(fit_type)
  turns    <- max(1, as.numeric(turns))

  ## ── which models to compute (lazy: skip models not required) ─────────────
  need_s <- fit_type %in% c("sphere",    "best")
  need_e <- fit_type %in% c("ellipsoid", "best")
  need_h <- fit_type %in% c("helix",     "best")
  need_r <- fit_type %in% c("ring",      "best")
  need_a <- fit_type %in% c("arc",       "best")
  .empty_dd <- list(rmsd = NA_real_, meandist = NA_real_, sd = NA_real_,
                    dist = numeric(0), nearest = NULL)

  ## ── resolve coordinate matrix (skip re-medianing if already provided) ----
  if (is.null(coords_matrix)) {
    if (is.null(df)) return(NULL)
    if (!is.null(cluster2_id)) df <- dplyr::filter(df, cluster2 == cluster2_id)
    if (nrow(df) == 0) return(NULL)
    medians <- if (use_cluster1 && "cluster" %in% names(df)) {
      df %>% dplyr::group_by(cluster) %>%
        dplyr::summarise(dplyr::across(c(x, y, z), median), .groups = "drop")
    } else {
      df %>% dplyr::group_by(tid) %>%
        dplyr::summarise(dplyr::across(c(x, y, z), median), .groups = "drop")
    }
    if (nrow(medians) == 0) return(NULL)
    coords_matrix <- as.matrix(medians[, c("x", "y", "z")])
  } else {
    medians <- as.data.frame(coords_matrix)
    colnames(medians) <- c("x", "y", "z")
  }

  med_pts <- coords_matrix

  ## ── inner helper: distances → RMSD + mean distance + SD ------------------
  rmsd_sd <- function(d) {
    d <- as.numeric(d); d <- d[is.finite(d) & d >= 0]
    list(
      rmsd     = if (length(d)) sqrt(mean(d^2)) else NA_real_,
      meandist = if (length(d)) mean(d)          else NA_real_,
      sd       = if (length(d) > 1) stats::sd(d) else NA_real_
    )
  }

  ## ── Sphere ───────────────────────────────────────────────────────────────
  dd_sphere <- if (need_s) {
    sph_geo <- fit_sphere_geometry(med_pts)
    pr_s    <- mle_project_sphere(med_pts, sph_geo$center, sph_geo$radius)
    d_s     <- as.numeric(pr_s$dist); rs_s <- rmsd_sd(d_s)
    list(rmsd=rs_s$rmsd, meandist=rs_s$meandist, sd=rs_s$sd, dist=d_s, nearest=pr_s$nearest)
  } else .empty_dd

  ## ── Ellipsoid ─────────────────────────────────────────────────────────────
  dd_ellipsoid <- if (need_e) {
    ell_geo <- fit_ellipsoid_geometry(med_pts)
    if (!is.null(ell_geo)) {
      pr_e <- mle_project_ellipsoid(med_pts, ell_geo$center, ell_geo$axes, ell_geo$R)
      d_e  <- as.numeric(pr_e$dist); rs_e <- rmsd_sd(d_e)
      list(rmsd=rs_e$rmsd, meandist=rs_e$meandist, sd=rs_e$sd, dist=d_e, nearest=pr_e$nearest)
    } else .empty_dd
  } else .empty_dd

  ## ── Helix ─────────────────────────────────────────────────────────────────
  dd_helix <- if (need_h) {
    hel_geo <- fit_helix_geometry(med_pts, turns = turns, n_points = helix_n)
    t_h <- .project_t_helix(med_pts, hel_geo$proj_center, hel_geo$axis,
                             hel_geo$radius, hel_geo$pitch, hel_geo$phi0,
                             check_neighbors = TRUE)
    H_h <- t(vapply(t_h, .helix_point, numeric(3),
                    c=hel_geo$proj_center, u=hel_geo$axis,
                    r=hel_geo$radius, p=hel_geo$pitch, phi0=hel_geo$phi0))
    d_h <- sqrt(rowSums((med_pts - H_h)^2)); rs_h <- rmsd_sd(d_h)
    list(rmsd=rs_h$rmsd, meandist=rs_h$meandist, sd=rs_h$sd, dist=d_h, nearest=H_h)
  } else .empty_dd

  ## ── Ring ──────────────────────────────────────────────────────────────────
  dd_ring <- if (need_r) {
    ring_geo <- fit_ring_geometry(med_pts)
    pr_r <- mle_project_ring(med_pts, ring_geo$center, ring_geo$normal, ring_geo$radius)
    d_r  <- as.numeric(pr_r$dist); rs_r <- rmsd_sd(d_r)
    list(rmsd=rs_r$rmsd, meandist=rs_r$meandist, sd=rs_r$sd, dist=d_r, nearest=pr_r$nearest)
  } else .empty_dd

  ## ── Arc ───────────────────────────────────────────────────────────────────
  dd_arc <- if (need_a) {
    arc_geo <- fit_arc_geometry(med_pts, arc_span_hint = arc_span_hint)
    pr_a <- mle_project_arc(med_pts, arc_geo$center, arc_geo$normal, arc_geo$radius,
                             arc_geo$t_center, arc_geo$arc_half)
    d_a  <- as.numeric(pr_a$dist); rs_a <- rmsd_sd(d_a)
    list(rmsd=rs_a$rmsd, meandist=rs_a$meandist, sd=rs_a$sd, dist=d_a, nearest=pr_a$nearest)
  } else .empty_dd

  ## ── choose which fit's distances to embed in the return dataframe ---------
  if (fit_type == "best") {
    rmsd_vals <- c(sphere = dd_sphere$rmsd, ellipsoid = dd_ellipsoid$rmsd,
                   helix  = dd_helix$rmsd,  ring = dd_ring$rmsd, arc = dd_arc$rmsd)
    rmsd_vals <- rmsd_vals[is.finite(rmsd_vals)]
    if (!length(rmsd_vals)) return(NULL)
    best_name <- names(which.min(rmsd_vals))
  } else {
    best_name <- fit_type
  }

  pick <- switch(best_name,
                 sphere    = dd_sphere,
                 ellipsoid = dd_ellipsoid,
                 helix     = dd_helix,
                 ring      = dd_ring,
                 arc       = dd_arc)
  if (is.null(pick) || !length(pick$dist)) return(NULL)

  ## ── build return dataframe (same columns as before for line-drawing compat)
  n <- nrow(med_pts)
  medians$distance_to_fit <- pick$dist[seq_len(n)]
  if (!is.null(pick$nearest) && nrow(pick$nearest) == n) {
    medians$nearest_x <- pick$nearest[, 1]
    medians$nearest_y <- pick$nearest[, 2]
    medians$nearest_z <- pick$nearest[, 3]
  } else {
    medians$nearest_x <- NA_real_
    medians$nearest_y <- NA_real_
    medians$nearest_z <- NA_real_
  }

  attr(medians, "rmsd_sphere")        <- dd_sphere$rmsd
  attr(medians, "meandist_sphere")    <- dd_sphere$meandist
  attr(medians, "sd_sphere")          <- dd_sphere$sd
  attr(medians, "rmsd_ellipsoid")     <- dd_ellipsoid$rmsd
  attr(medians, "meandist_ellipsoid") <- dd_ellipsoid$meandist
  attr(medians, "sd_ellipsoid")       <- dd_ellipsoid$sd
  attr(medians, "rmsd_helix")         <- dd_helix$rmsd
  attr(medians, "meandist_helix")     <- dd_helix$meandist
  attr(medians, "sd_helix")           <- dd_helix$sd
  attr(medians, "rmsd_ring")          <- dd_ring$rmsd
  attr(medians, "meandist_ring")      <- dd_ring$meandist
  attr(medians, "sd_ring")            <- dd_ring$sd
  attr(medians, "rmsd_arc")           <- dd_arc$rmsd
  attr(medians, "meandist_arc")       <- dd_arc$meandist
  attr(medians, "sd_arc")             <- dd_arc$sd
  attr(medians, "best_fit")           <- best_name

  return(medians)
}
# ---- NN HELPERS  ----
safe_k <- function(n_data, k) {
  # nn2 requires 1 <= k <= nrow(data)
  max(1, min(k, n_data))
}

safe_nn2 <- function(X, k = 4) {
  X <- as.matrix(X)
  k <- max(1, min(k, nrow(X)))
  out <- tryCatch(RANN::nn2(X, k = k), error = function(e) NULL)
  if (is.null(out)) return(list(nn = NULL, k_eff = 0))
  # Drop all-NA rows (shouldn’t happen, but be safe)
  if (any(!is.finite(out$nn.dists))) {
    bad <- !is.finite(rowSums(out$nn.dists))
    if (any(bad)) {
      out$nn.dists <- out$nn.dists[!bad, , drop = FALSE]
      out$nn.idx   <- out$nn.idx  [!bad, , drop = FALSE]
    }
  }
  list(nn = out, k_eff = ncol(out$nn.idx))
}

# --- color helper for fits ---
fit_color_for_subset <- function(subdf, fallback_hex = "#00BFFF") {
  co <- unique(na.omit(as.character(subdf$cluster2_color)))
  if (length(co) && nzchar(co[1])) co[1] else fallback_hex
}
build_center_args <- function(X, input) {
  c_robust <- robust_center(X)
  switch(input$mle_center_mode,
         "fixed" = list(center_fixed = c_robust),
         "reg"   = list(center_init = c_robust,
                        center_prior = c_robust,
                        lambda = input$mle_center_lambda),
         "free"  = list(center_init = c_robust)
  )
}
# ---------- MLE helpers (orthogonal-distance) ----------
# (A) Sphere
# Robust center helper (fallback to geometric median if hull fails)
# --- Robust center guess (hull centroid -> geometric median -> mean) ---
robust_center <- function(X) {
  tryCatch(
    compute_convex_hull_centroid(X),
    error = function(e1) {
      tryCatch(as.numeric(Gmedian::Gmedian(X)), error = function(e2) colMeans(X))
    }
  )
}
# =============================================================================
# SHAPE FITTING — fit a geometric primitive to a cluster of localizations
# =============================================================================
#
# GOAL OF THIS WHOLE BLOCK
#   Given the points belonging to one cluster, ask "what SHAPE is this?" and fit
#   the best sphere / ellipsoid / helix / ring / arc to it. That turns a cloud of
#   localizations into interpretable numbers: a radius, a pitch, an axis length.
#
# THE FITTING PRINCIPLE — ORTHOGONAL-DISTANCE (total least squares)
#   For a shape S with parameters θ, minimise the sum of squared PERPENDICULAR
#   distances from each point to the SURFACE:
#
#        θ* = argmin_θ  Σ_i  || x_i − P_S(x_i; θ) ||²
#
#   where P_S(x) is the nearest point on the shape. This is the right criterion
#   here because localization error is isotropic in space — there is no
#   privileged "vertical" axis to measure residuals along, as ordinary least
#   squares would assume. It is the same reason these are called MLE fits: for
#   isotropic Gaussian error, minimising orthogonal distance IS the maximum
#   likelihood solution.
#
#   Because P_S has no closed form for the curved shapes, several fits iterate:
#   project -> re-fit parameters -> re-project, until the residual stops falling.
#
# THE NAMING CONVENTION (this is the key to navigating ~30 functions)
#   For each shape S in {sphere, ellipsoid, helix, ring, arc} there is a family:
#
#     fit_S_geometry(X) ...... convenience wrapper: fit and return the geometry
#                              in one call (what most callers actually use)
#     mle_fit_S(X, w) ........ the real optimiser. Returns the fitted parameters
#                              plus the residual/RMSD used for model ranking.
#     mle_project_S(X, ...) .. nearest point ON the shape for each input point.
#                              This is what residuals and the "distance to fit"
#                              lines are computed from.
#     mle_sample_S(...) ...... sample points along/over the fitted shape, used
#                              to DRAW it (not to fit it).
#     make_S_mesh() /
#       make_S_tube_mesh() ... turn the sampled shape into plotly 3D geometry.
#                              "tube" variants sweep a small circle along a
#                              curve, because a 1D curve is invisible in 3D.
#     generate_S_from_cluster2(...)  per-cluster wrapper used by the server.
#
#   So: to change how a shape is FITTED edit mle_fit_S; to change how it LOOKS
#   edit mle_sample_S / make_S_*_mesh. They are deliberately independent.
#
# MODEL SELECTION — which shape wins?
#   RMSD always favours the most flexible shape (more parameters always fit
#   better), so ranking uses AICc/BIC, which penalise parameter count. AICc is
#   preferred for the small point counts typical of one cluster, and is left
#   blank when a cluster has too few points to justify that model (n <= params+2).
#
# HONEST CAVEAT (mirrored by the EXPLORATORY badge in the UI)
#   The shape dictionary is IMPOSED by this tool — a cluster is scored only
#   against the five shapes offered, so "best" means best OF THESE, never proof
#   of biological geometry. With few, blinking-inflated localizations several
#   shapes are often statistically indistinguishable. Treat the ranking as a
#   hypothesis, and check the residuals rather than the label.
# =============================================================================

# Triangle index list for a parametric (nu x nv) surface grid — the connectivity
# plotly's mesh3d needs. Splits every grid quad into two triangles. Shared by all
# the surface meshes below, which differ only in their vertex coordinates.
.make_triangles <- function(nu, nv) {
  # grid indices (0-based for plotly mesh3d i/j/k)
  idx <- matrix(seq_len(nu*nv) - 1L, nrow = nu, ncol = nv, byrow = TRUE)
  i <- c(); j <- c(); k <- c()
  for (u in 1:(nu-1)) {
    for (v in 1:(nv-1)) {
      a <- idx[u    , v    ]
      b <- idx[u + 1, v    ]
      c1<- idx[u    , v + 1]
      d <- idx[u + 1, v + 1]
      # two triangles per quad
      i <- c(i, a, b)
      j <- c(j, b, d)
      k <- c(k, c1, c1)
      i <- c(i, a)
      j <- c(j, d)
      k <- c(k, c1)
    }
  }
  list(i = i, j = j, k = k)
}
# Solid color for mesh3d across plotly versions
mesh_color_args <- function(hex, nfaces) {
  hex <- as.character(hex %||% "#bdbdbd")
  if (!nzchar(hex)) hex <- "#bdbdbd"
  list(
    color      = I(hex),                
    facecolor  = rep(hex, nfaces),      
    showscale  = FALSE,
    intensity  = NULL,                  
    colorscale = NULL,
    cmin       = NULL,
    cmax       = NULL,
    flatshading = TRUE
  )
}
make_sphere_mesh <- function(center, radius, nu = 48, nv = 24) {
  u <- seq(0, 2*pi, length.out = nu)
  v <- seq(0,   pi, length.out = nv)
  UU <- rep(u, each = nv)
  VV <- rep(v, times = nu)
  x <- radius * cos(UU) * sin(VV) + center[1]
  y <- radius * sin(UU) * sin(VV) + center[2]
  z <- radius * cos(VV)           + center[3]
  tris <- .make_triangles(nu, nv)
  list(x = x, y = y, z = z, i = tris$i, j = tris$j, k = tris$k)
}

make_ellipsoid_mesh <- function(center, axes, R, nu = 64, nv = 32) {
  u <- seq(0, 2*pi, length.out = nu)
  v <- seq(0,   pi, length.out = nv)
  UU <- rep(u, each = nv)
  VV <- rep(v, times = nu)
  # local parametric ellipsoid
  xl <- axes[1] * cos(UU) * sin(VV)
  yl <- axes[2] * sin(UU) * sin(VV)
  zl <- axes[3] * cos(VV)
  M  <- cbind(xl, yl, zl) %*% t(R)
  x  <- M[,1] + center[1]
  y  <- M[,2] + center[2]
  z  <- M[,3] + center[3]
  tris <- .make_triangles(nu, nv)
  list(x = x, y = y, z = z, i = tris$i, j = tris$j, k = tris$k)
}
make_alpha_mesh3d <- function(pts_xyz, alpha) {
  sh <- alphashape3d::ashape3d(pts_xyz, alpha = alpha)
  tri <- sh$triang
  if (is.null(tri) || !nrow(tri)) return(NULL)
  
  v <- sh$x  # vertices (n x 3)
  
  # plotly mesh3d uses 0-based indices
  list(
    x = v[, 1], y = v[, 2], z = v[, 3],
    i = tri[, 1] - 1L,
    j = tri[, 2] - 1L,
    k = tri[, 3] - 1L
  )
}
# ========= Sphere (free / fixed / regularized) =========
mle_fit_sphere <- function(X, w = NULL,
                           center_init  = NULL,
                           center_fixed = NULL,
                           center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 4) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)
  
  if (!is.null(center_fixed)) {
    c <- as.numeric(center_fixed)
    L <- sqrt(rowSums((X - matrix(c, n, 3, TRUE))^2))
    r <- sum(w * L) / sum(w)
    return(list(center = c, radius = as.numeric(abs(r))))
  }
  
  c0 <- if (!is.null(center_init)) as.numeric(center_init) else colMeans(X)
  d0 <- sqrt(rowSums((X - matrix(c0, n, 3, TRUE))^2))
  r0 <- median(d0)
  # --- NEW: robust upper bound for radius ---
  r_max <- max(stats::quantile(d0, 0.97, na.rm = TRUE), r0 * 2)
  
  par0 <- c(c0, r0)
  # Bound center to an expanded data bounding box and r to [0, r_max].
  # Without centre bounds, the joint (c, r) objective is non-convex: the
  # optimizer can drift c outside the data cloud and inflate r to compensate
  # while keeping all residuals small — a geometrically meaningless solution.
  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)  # ≥10 nm margin
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand
  obj <- function(par) {
    c <- par[1:3]; r <- abs(par[4])
    L <- sqrt(rowSums((X - matrix(c, n, 3, TRUE))^2))
    # Divide by n so the objective is mean squared residual (nm²), independent of
    # cluster size.  This makes lambda directly interpretable across clusters:
    # lambda ≈ 1 means a 1 nm centre displacement costs as much as 1 nm² of mean residual.
    val <- sum(w * (L - r)^2) / n
    if (!is.null(center_prior) && lambda > 0) val <- val + lambda * sum((c - center_prior)^2)
    val
  }
  fit <- optim(
    par0, obj, method = "L-BFGS-B", control = list(maxit = 500),
    lower = c(c_lo, 0),       # center bounded to expanded data box; r >= 0
    upper = c(c_hi, r_max)    # r <= r_max (hard bound replaces soft penalty)
  )
  list(center = fit$par[1:3], radius = abs(fit$par[4]))
}
mle_project_sphere <- function(X, center, radius) {
  V <- sweep(X, 2, center); L <- sqrt(rowSums(V^2))
  scale <- ifelse(L > 0, radius / L, 0)
  Q <- sweep(V, 1, scale, "*")
  nearest <- sweep(Q, 2, center, "+")
  dist <- abs(L - radius)
  list(nearest = nearest, dist = dist)
}
mle_sample_sphere <- function(center, radius, n = 4000) {
  phi <- runif(n, 0, 2*pi)
  u   <- runif(n, -1, 1); theta <- acos(u)
  x <- radius * sin(theta) * cos(phi)
  y <- radius * sin(theta) * sin(phi)
  z <- radius * cos(theta)
  cbind(x, y, z) + matrix(center, n, 3, TRUE)
}

# (B) Ellipsoid  (rotation via yaw/pitch/roll)
.eul2R <- function(yaw, pitch, roll) {
  Ry <- matrix(c(cos(yaw),-sin(yaw),0, sin(yaw),cos(yaw),0, 0,0,1),3,3,TRUE)
  Rp <- matrix(c(cos(pitch),0,sin(pitch), 0,1,0, -sin(pitch),0,cos(pitch)),3,3,TRUE)
  Rr <- matrix(c(1,0,0, 0,cos(roll),-sin(roll), 0,sin(roll),cos(roll)),3,3,TRUE)
  Ry %*% Rp %*% Rr
}
# ========= Ellipsoid (free / fixed / regularized center) =========
# Objective: sum of squared Euclidean distances from each localization to the
# nearest point on the ellipsoid surface — identical units (nm²) to the sphere
# fitter and to the reported RMSD, making fits directly comparable.
#
# Previous version minimized the dimensionless radial ratio sum((s-1)^2) where
# s = ||Xp/axes||.  That objective has a scale degeneracy: the optimizer could
# translate the centre and simultaneously scale the axes while keeping s ≈ 1
# for all points, producing centres far outside the data cloud.  The orthogonal-
# distance objective eliminates this degeneracy because any displacement of the
# centre away from the data increases the distances unambiguously.
#
# Axis bounds (log-space) prevent runaway growth; centre bounds (expanded data
# bounding box) catch any remaining drift for the free-centre case.
mle_fit_ellipsoid <- function(X, w = NULL,
                              center_init  = NULL,
                              center_fixed = NULL,
                              center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 6) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)

  # Shared initialisation: PCA at a given centre + axis/rotation starting values
  .ell_init <- function(cc) {
    Xc <- sweep(X, 2, cc)
    pc <- prcomp(Xc, center = FALSE); R0 <- pc$rotation
    Q  <- Xc %*% R0
    a0 <- pmax(diff(stats::quantile(Q[,1], c(0.05, 0.95))) / 2, 1e-6)
    b0 <- pmax(diff(stats::quantile(Q[,2], c(0.05, 0.95))) / 2, 1e-6)
    c0 <- pmax(diff(stats::quantile(Q[,3], c(0.05, 0.95))) / 2, 1e-6)
    yaw0   <- atan2(R0[2,1], R0[1,1])
    pitch0 <- asin(max(-1, min(1, -R0[3,1])))  # clamp avoids NaN at ±1
    roll0  <- atan2(R0[3,2], R0[3,3])
    # Log-axis bounds: axes must stay in [1% of smallest initial half-extent, 5× RMS spread]
    d_rms  <- sqrt(mean(rowSums(Xc^2)))
    ax_lo  <- log(pmax(min(a0, b0, c0) * 0.01, 1e-3))
    ax_hi  <- log(max(d_rms * 5, 1))
    list(a0 = a0, b0 = b0, c0 = c0, R0 = R0,
         yaw0 = yaw0, pitch0 = pitch0, roll0 = roll0,
         ax_lo = ax_lo, ax_hi = ax_hi)
  }

  if (!is.null(center_fixed)) {
    cc  <- as.numeric(center_fixed)
    ini <- .ell_init(cc)
    par0 <- c(log(c(ini$a0, ini$b0, ini$c0)), ini$yaw0, ini$pitch0, ini$roll0)
    obj <- function(par) {
      axes_now <- exp(par[1:3])
      R  <- .eul2R(par[4], par[5], par[6])
      pr <- mle_project_ellipsoid(X, cc, axes_now, R)
      sum(w * pr$dist^2) / n   # mean squared distance — scale-independent of n
    }
    fit <- optim(par0, obj, method = "L-BFGS-B",
                 lower = c(rep(ini$ax_lo, 3), -Inf, -Inf, -Inf),
                 upper = c(rep(ini$ax_hi, 3),  Inf,  Inf,  Inf),
                 control = list(maxit = 500))
    p <- fit$par
    return(list(center = cc, axes = exp(p[1:3]), R = .eul2R(p[4], p[5], p[6])))
  }

  # Free centre (optionally regularized toward center_prior)
  cc0 <- if (!is.null(center_init)) as.numeric(center_init) else colMeans(X)
  ini <- .ell_init(cc0)
  par0 <- c(cc0, log(c(ini$a0, ini$b0, ini$c0)), ini$yaw0, ini$pitch0, ini$roll0)

  obj <- function(par) {
    cc       <- par[1:3]
    axes_now <- exp(par[4:6])
    R        <- .eul2R(par[7], par[8], par[9])
    pr  <- mle_project_ellipsoid(X, cc, axes_now, R)
    val <- sum(w * pr$dist^2) / n   # mean squared distance — scale-independent of n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((cc - center_prior)^2)
    val
  }

  # Centre bounded to expanded data bounding box to prevent drift
  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand

  fit <- optim(par0, obj, method = "L-BFGS-B",
               lower = c(c_lo, rep(ini$ax_lo, 3), -Inf, -Inf, -Inf),
               upper = c(c_hi, rep(ini$ax_hi, 3),  Inf,  Inf,  Inf),
               control = list(maxit = 500))
  p <- fit$par
  list(center = p[1:3], axes = exp(p[4:6]), R = .eul2R(p[7], p[8], p[9]))
}
mle_project_ellipsoid <- function(X, center, axes, R) {
  # Exact orthogonal nearest-point via vectorised Newton on the KKT condition.
  # Nearest point p to query q (in PCA frame) satisfies:
  #   p_j = q_j * a_j^2 / (a_j^2 + lambda)
  # where lambda solves  F(lambda) = sum_j (q_j*a_j / (a_j^2+lambda))^2 - 1 = 0
  # F'(lambda) = -2 * sum_j q_j^2 * a_j^2 / (a_j^2+lambda)^3
  n    <- nrow(X)
  Xp   <- (X - matrix(center, n, 3, TRUE)) %*% R   # in PCA frame  (n x 3)
  a2   <- axes^2                                     # [a^2, b^2, c^2]
  Xp2  <- Xp^2                                       # precompute  (n x 3)
  r2   <- rowSums(Xp2)                               # squared distance from center

  min_a2 <- min(a2)
  # Detect interior points: sum_j (q_j/a_j)^2 < 1
  inside_ratio <- rowSums(Xp2 / matrix(a2, n, 3, TRUE))
  lam_outside  <- pmax(max(axes), 1) * sqrt(pmax(r2, 1e-20)) - max(a2)
  # Interior points: lambda root lies in (-min_a2, 0); start at midpoint to avoid clamp divergence
  is_interior <- inside_ratio < 1
  lam <- ifelse(is_interior, -min_a2 * 0.5, lam_outside)
  lam <- pmax(lam, -min_a2 + 1e-6)

  for (iter in 1:25) {
    denom <- outer(lam, a2, "+")                      # n x 3 : a_j^2 + lambda_i
    denom <- pmax(denom, 1e-12)
    A2X2  <- matrix(a2, n, 3, TRUE) * Xp2            # a_j^2 * q_j^2
    fval  <- rowSums(A2X2 / denom^2) - 1             # F(lambda)
    fpval <- -2 * rowSums(A2X2 / denom^3)            # F'(lambda)
    step  <- ifelse(abs(fpval) > 1e-30, fval / fpval, 0)
    lam   <- pmax(lam - step, -min_a2 + 1e-6)        # keep denominators positive
    lam[is_interior] <- pmin(lam[is_interior], -1e-12) # interior root < 0; prevent overshoot above 0
    if (max(abs(step)) < 1e-8) break
  }

  denom_f       <- pmax(outer(lam, a2, "+"), 1e-12)  # n x 3
  nearest_local <- matrix(a2, n, 3, TRUE) * Xp / denom_f  # p_j = a_j^2*q_j/(a_j^2+lambda)
  nearest       <- nearest_local %*% t(R) + matrix(center, n, 3, TRUE)
  dist          <- sqrt(rowSums((X - nearest)^2))

  # Guard: verify the projected point actually lies on the ellipsoid surface.
  # A correct KKT solution satisfies sum_j (p_j/a_j)^2 == 1 exactly.
  # If Newton diverged (lambda stuck near -min_a2), nearest_local is astronomically
  # large → surface_check >> 1 → fall back to radial projection.
  surface_check <- rowSums((nearest_local / matrix(axes, n, 3, TRUE))^2)
  runaway <- !is.finite(surface_check) | abs(surface_check - 1) > 0.1 | !is.finite(dist)

  if (any(runaway)) {
    # Radial fallback: project Xp along its own direction to the ellipsoid surface.
    # Point p = Xp * t lands on surface when sum_j (t*q_j/a_j)^2 = 1
    # → t = 1/sqrt(inside_ratio).  Always finite; always on-surface.
    ir_run  <- pmax(inside_ratio[runaway], 1e-20)
    t_rad   <- 1 / sqrt(ir_run)
    Xp_run  <- Xp[runaway, , drop = FALSE]
    xp_norm <- sqrt(rowSums(Xp_run^2))
    nn_local_rad <- Xp_run * t_rad

    # For points at the center (Xp ≈ 0), radial projection is undefined;
    # place the nearest point on the tip of the shortest semi-axis.
    near_ctr <- xp_norm < 1e-12
    if (any(near_ctr)) {
      j_min <- which.min(axes)
      p_short <- numeric(3); p_short[j_min] <- axes[j_min]
      nn_local_rad[near_ctr, ] <- matrix(p_short, sum(near_ctr), 3, TRUE)
    }

    nearest[runaway, ] <- nn_local_rad %*% t(R) + matrix(center, sum(runaway), 3, TRUE)
    dist[runaway]      <- sqrt(rowSums((X[runaway, , drop = FALSE] -
                                        nearest[runaway, , drop = FALSE])^2))
  }

  # Edge case: point at ellipsoid center -> nearest on shortest axis
  at_ctr <- r2 < 1e-20
  if (any(at_ctr)) {
    min_ax            <- which.min(axes)
    p_deg             <- numeric(3); p_deg[min_ax] <- axes[min_ax]
    nearest[at_ctr, ] <- matrix(p_deg %*% t(R) + center, sum(at_ctr), 3, TRUE)
    dist[at_ctr]      <- axes[min_ax]
  }
  list(nearest = nearest, dist = dist)
}
mle_sample_ellipsoid <- function(center, axes, R, n = 4000) {
  phi <- runif(n, 0, 2*pi)
  u   <- runif(n, -1, 1); theta <- acos(u)
  xl <- axes[1]*sin(theta)*cos(phi)
  yl <- axes[2]*sin(theta)*sin(phi)
  zl <- axes[3]*cos(theta)
  (cbind(xl,yl,zl) %*% t(R)) + matrix(center, n, 3, TRUE)
}

# (C) Helix
.ang_to_u <- function(th, ph) c(sin(th)*cos(ph), sin(th)*sin(ph), cos(th))
.rot_to_u <- function(u) {
  ez <- c(0,0,1); v <- pracma::cross(ez, u); s <- sqrt(sum(v*v)); c <- sum(ez*u)
  if (s < 1e-12) return(diag(3))
  Vx <- matrix(c( 0,-v[3], v[2], v[3], 0,-v[1], -v[2], v[1], 0), 3,3, TRUE)
  diag(3) + Vx + Vx%*%Vx*((1-c)/(s^2))
}
.helix_point <- function(t, c, u, r, p, phi0) {
  R <- .rot_to_u(u)
  H <- cbind(r*cos(t+phi0), r*sin(t+phi0), (p/(2*pi))*t)
  drop(H %*% t(R) + c)
}
.helix_newton <- function(t, Xloc, r, phi0, a, iters) {
  for (k in seq_len(iters)) {
    cx  <- r * cos(t + phi0); sx <- r * sin(t + phi0)
    fx  <- cbind(cx, sx, a * t) - Xloc
    d   <- cbind(-r * sin(t + phi0), r * cos(t + phi0), a)
    t   <- t - rowSums(d * fx) / (rowSums(d * d) + 1e-12)
  }
  t
}
.helix_d2 <- function(t, Xloc, r, phi0, a) {
  H <- cbind(r * cos(t + phi0), r * sin(t + phi0), a * t)
  rowSums((H - Xloc)^2)
}
.project_t_helix <- function(X, c, u, r, p, phi0, t0 = NULL, iters = 15,
                              check_neighbors = FALSE) {
  n <- nrow(X); R <- .rot_to_u(u); Xloc <- (X - matrix(c, n, 3, TRUE)) %*% R
  a <- p/(2*pi)
  if (is.null(t0)) {
    # Angular estimate: which angle on any turn?
    t_xy <- atan2(Xloc[,2], Xloc[,1]) - phi0
    if (abs(a) > 1e-12) {
      # Axial estimate: helix z-coordinate = a*t, so t ≈ z/a.
      # Use this to pick which turn (n_turns) and angular position for fine-tuning.
      # This prevents Newton from converging to a local minimum on the wrong turn.
      t_z     <- Xloc[,3] / a
      n_turns <- round((t_z - t_xy) / (2*pi))
      t0 <- t_xy + n_turns * 2*pi
    } else {
      t0 <- t_xy
    }
  }
  t <- .helix_newton(t0, Xloc, r, phi0, a, iters)

  # Check adjacent turns (±2π) and keep whichever gives shortest distance.
  # Handles cases where the axial/angular initial estimate lands on the wrong turn.
  if (check_neighbors && abs(a) > 1e-12) {
    t_m <- .helix_newton(t - 2*pi, Xloc, r, phi0, a, iters)
    t_p <- .helix_newton(t + 2*pi, Xloc, r, phi0, a, iters)
    d0  <- .helix_d2(t,   Xloc, r, phi0, a)
    dm  <- .helix_d2(t_m, Xloc, r, phi0, a)
    dp  <- .helix_d2(t_p, Xloc, r, phi0, a)
    t   <- ifelse(dm < pmin(d0, dp), t_m,
                  ifelse(dp < pmin(d0, dm), t_p, t))
  }
  t
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Coarse half-sphere orientation grid search.
# Sweeps n_theta × n_phi candidate (theta, phi) pairs covering the half-sphere
# (theta 0→π/2; the normal/axis and its negation define the same geometry).
# obj_fn(theta, phi) must return a scalar objective (lower = better).
# Returns the (theta, phi) pair with the lowest objective, or the PCA fallback
# supplied via th_init/ph_init if no grid candidate beats it.
.grid_search_orient <- function(obj_fn, th_init, ph_init, n_theta = 12L, n_phi = 16L) {
  # Strip any R names so the return vector has clean "theta"/"phi" keys.
  th_init <- unname(as.numeric(th_init))
  ph_init <- unname(as.numeric(ph_init))

  # cos-spaced theta rows → roughly equal solid-angle coverage
  th_vals <- acos(seq(1, 0, length.out = n_theta))   # 0 → π/2
  ph_vals <- seq(0, 2*pi, length.out = n_phi + 1L)[-(n_phi + 1L)]

  pca_val  <- tryCatch(obj_fn(th_init, ph_init), error = function(e) Inf)
  best_val <- pca_val
  best_th  <- th_init
  best_ph  <- ph_init

  for (th in th_vals) {
    for (ph in ph_vals) {
      v <- tryCatch(obj_fn(th, ph), error = function(e) Inf)
      if (is.finite(v) && v < best_val) {
        best_val <- v; best_th <- th; best_ph <- ph
      }
    }
  }
  c(theta = best_th, phi = best_ph)
}

mle_fit_helix <- function(X, w = NULL, turns_hint = 3, max_outer = 20,
                          center_init  = NULL,
                          center_fixed = NULL,
                          center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 6) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)
  
  # initiation
  c0 <- if (!is.null(center_init)) as.numeric(center_init) else colMeans(X)
  Xc <- sweep(X, 2, c0)
  pc <- prcomp(Xc, center = FALSE); u0 <- pc$rotation[,1]
  v1 <- pc$rotation[,2]; v2 <- pracma::cross(u0, v1)
  XY <- Xc %*% cbind(v1, v2, u0)
  r0 <- median(sqrt(XY[,1]^2 + XY[,2]^2))
  h  <- diff(range(XY[,3])); p0 <- if (turns_hint > 0) h/turns_hint else h
  par <- c(c0, th = acos(u0[3]), ph = atan2(u0[2], u0[1]), r = r0, p = p0, phi0 = 0)
  t   <- .project_t_helix(X, par[1:3], .ang_to_u(par[4],par[5]), par[6], par[7], par[8])

  # ── Orientation grid search: finds best axis direction before L-BFGS-B ────
  # Cheap radial proxy: distance from candidate axis should equal r0.
  # Avoids Newton projection; captures the right basin for the optimizer.
  c_gs <- if (!is.null(center_fixed)) as.numeric(center_fixed) else par[1:3]
  gs_h <- .grid_search_orient(
    obj_fn = function(th, ph) {
      ax  <- .ang_to_u(th, ph)
      Xc  <- sweep(X, 2, c_gs)
      d   <- as.numeric(Xc %*% ax)
      Xp  <- Xc - outer(d, ax)
      sum(w * (sqrt(rowSums(Xp^2)) - r0)^2) / n
    },
    th_init = par[4], ph_init = par[5]
  )
  if (gs_h[1] != par[4] || gs_h[2] != par[5]) {
    par[4] <- gs_h[1]; par[5] <- gs_h[2]
    t <- .project_t_helix(X, par[1:3], .ang_to_u(par[4], par[5]),
                           par[6], par[7], par[8])
  }

  for (outer in 1:max_outer) {
    obj <- function(pv) {
      # choose center according to mode
      c <- if (!is.null(center_fixed)) as.numeric(center_fixed) else pv[1:3]
      u <- .ang_to_u(pv[4],pv[5]); r <- abs(pv[6]); p <- pv[7]; phi0 <- pv[8]
      H <- t(vapply(t, .helix_point, numeric(3), c=c, u=u, r=r, p=p, phi0=phi0))
      res2 <- rowSums((X - H)^2)
      val <- sum(w * res2) / n   # mean squared distance — scale-independent of n
      if (is.null(center_fixed) && !is.null(center_prior) && lambda > 0) {
        val <- val + lambda * sum((c - center_prior)^2)
      }
      val
    }
    # Optimization dimension depends on whether c is fixed
    if (!is.null(center_fixed)) {
      par_opt <- par[-(1:3)]
      fit <- optim(par_opt, function(q) obj(c(par[1:3], q)),
                   method = "L-BFGS-B", control = list(maxit = 200))
      par[-(1:3)] <- fit$par
    } else {
      fit <- optim(par, obj, method = "L-BFGS-B", control = list(maxit = 200))
      par <- fit$par
    }
    c_now <- if (is.null(center_fixed)) par[1:3] else as.numeric(center_fixed)
    t_new <- .project_t_helix(X, c_now, .ang_to_u(par[4],par[5]), abs(par[6]), par[7], par[8], t0 = t)
    if (mean(abs(t_new - t)) < 1e-6) { t <- t_new; break } else t <- t_new
  }
  c <- if (is.null(center_fixed)) par[1:3] else as.numeric(center_fixed)
  u <- .ang_to_u(par[4],par[5]); r <- abs(par[6]); p <- par[7]; phi0 <- par[8]
  H <- t(vapply(t, .helix_point, numeric(3), c=c, u=u, r=r, p=p, phi0=phi0))
  list(center = c, axis = u, radius = r, pitch = p, phi0 = phi0, t = t, H = H)
}
mle_sample_helix <- function(fit, n = 3000, pad_turns = 0.25) {
  tmin <- min(fit$t) - 2*pi*pad_turns
  tmax <- max(fit$t) + 2*pi*pad_turns
  tt <- seq(tmin, tmax, length.out = n)
  H  <- t(vapply(tt, .helix_point, numeric(3),
                 c = fit$center, u = fit$axis, r = fit$radius, p = fit$pitch, phi0 = fit$phi0))
  as.data.frame(H)
}
# ── (D) Ring ─────────────────────────────────────────────────────────────────
# A 3D ring (circle) lives in a plane whose unit normal is n.
# Orthogonal distance from any point to the ring = sqrt((r_plane - R)^2 + d_ax^2)
# where r_plane = in-plane distance from the ring axis, d_ax = axial offset.

.ring_basis <- function(normal) {
  # Return [u | v] — two orthonormal vectors spanning the plane ⊥ to normal.
  n   <- as.numeric(normal) / sqrt(sum(normal^2))
  ref <- if (abs(n[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0)
  u   <- pracma::cross(ref, n); u <- u / sqrt(sum(u^2))
  v   <- pracma::cross(n, u)           # already unit (n ⊥ u, both unit)
  cbind(u, v)
}

mle_project_ring <- function(X, center, normal, radius) {
  n   <- as.numeric(normal) / sqrt(sum(normal^2))
  Xc  <- sweep(X, 2, center)
  d_ax <- as.numeric(Xc %*% n)
  Xp   <- Xc - outer(d_ax, n)          # in-plane component
  r_pl <- sqrt(rowSums(Xp^2))           # distance from ring axis in-plane

  # Unit in-plane direction; fall back to ring basis x-axis when on the axis
  basis  <- .ring_basis(n)
  safe   <- r_pl > 1e-12
  dir    <- matrix(basis[, 1], nrow(X), 3, TRUE)
  dir[safe, ] <- Xp[safe, , drop = FALSE] / r_pl[safe]

  nearest <- matrix(center, nrow(X), 3, TRUE) + radius * dir
  dist    <- sqrt((r_pl - radius)^2 + d_ax^2)
  list(nearest = nearest, dist = dist)
}

mle_fit_ring <- function(X, w = NULL,
                         center_init  = NULL,
                         center_fixed = NULL,
                         center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 4) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)

  geo  <- fit_ring_geometry(X)
  c0   <- if (!is.null(center_init)) as.numeric(center_init) else geo$center
  nrm0 <- geo$normal
  th0  <- acos(max(-1, min(1, nrm0[3])))
  ph0  <- atan2(nrm0[2], nrm0[1])
  r0   <- geo$radius
  r_max <- max(sqrt(rowSums(sweep(X, 2, c0)^2))) * 2

  # ── Orientation grid search: best normal direction before L-BFGS-B ────────
  c_gs <- if (!is.null(center_fixed)) as.numeric(center_fixed) else c0
  gs_r <- .grid_search_orient(
    obj_fn  = function(th, ph) sum(w * mle_project_ring(X, c_gs, .ang_to_u(th, ph), r0)$dist^2) / n,
    th_init = th0, ph_init = ph0
  )
  th0 <- gs_r[1]; ph0 <- gs_r[2]

  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand

  if (!is.null(center_fixed)) {
    cc   <- as.numeric(center_fixed)
    par0 <- c(th0, ph0, log(max(r0, 1e-3)))
    obj  <- function(par) {
      pr <- mle_project_ring(X, cc, .ang_to_u(par[1], par[2]), exp(par[3]))
      sum(w * pr$dist^2) / n
    }
    fit <- optim(par0, obj, method = "L-BFGS-B",
                 lower = c(-Inf, -Inf, log(1e-3)),
                 upper = c( Inf,  Inf, log(r_max)),
                 control = list(maxit = 500))
    p <- fit$par
    return(list(center = cc, normal = .ang_to_u(p[1], p[2]), radius = exp(p[3])))
  }

  par0 <- c(c0, th0, ph0, log(max(r0, 1e-3)))
  obj  <- function(par) {
    pr  <- mle_project_ring(X, par[1:3], .ang_to_u(par[4], par[5]), exp(par[6]))
    val <- sum(w * pr$dist^2) / n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((par[1:3] - center_prior)^2)
    val
  }
  fit <- optim(par0, obj, method = "L-BFGS-B",
               lower = c(c_lo, -Inf, -Inf, log(1e-3)),
               upper = c(c_hi,  Inf,  Inf, log(r_max)),
               control = list(maxit = 500))
  p <- fit$par
  list(center = p[1:3], normal = .ang_to_u(p[4], p[5]), radius = exp(p[6]))
}

mle_sample_ring <- function(fit, n = 500) {
  basis <- .ring_basis(fit$normal)
  tt    <- seq(0, 2*pi, length.out = n + 1L)[-(n + 1L)]
  pts   <- fit$radius * (outer(cos(tt), basis[, 1]) + outer(sin(tt), basis[, 2])) +
           matrix(fit$center, n, 3, TRUE)
  colnames(pts) <- c("x", "y", "z")
  as.data.frame(pts)
}

# ── (E) Arc ───────────────────────────────────────────────────────────────────
# An arc is a partial ring: same center/normal/radius, plus t_center (midpoint
# angle in the ring-plane frame) and arc_half (half-span, 0 < arc_half <= π).
# The arc runs from t_center − arc_half to t_center + arc_half.

mle_project_arc <- function(X, center, normal, radius, t_center, arc_half) {
  n     <- as.numeric(normal) / sqrt(sum(normal^2))
  basis <- .ring_basis(n)
  u <- basis[, 1]; vv <- basis[, 2]

  Xc   <- sweep(X, 2, center)
  d_ax <- as.numeric(Xc %*% n)
  Xp   <- Xc - outer(d_ax, n)

  # Angle of each in-plane projection in ring frame
  t_ang <- atan2(as.numeric(Xp %*% vv), as.numeric(Xp %*% u))

  # Wrap relative angle to [−π, π] then clamp to the arc half-span
  dt      <- ((t_ang - t_center + pi) %% (2*pi)) - pi
  dt_cl   <- pmax(-arc_half, pmin(arc_half, dt))
  t_near  <- t_center + dt_cl

  nearest <- radius * (outer(cos(t_near), u) + outer(sin(t_near), vv)) +
             matrix(center, nrow(X), 3, TRUE)
  dist    <- sqrt(rowSums((X - nearest)^2))
  list(nearest = nearest, dist = dist)
}

mle_fit_arc <- function(X, w = NULL,
                        center_init  = NULL,
                        center_fixed = NULL,
                        center_prior = NULL, lambda = 0,
                        arc_span_hint = 0) {
  # The arc_half parameter has a monotone incentive to grow (more arc coverage
  # = lower sum-of-squared-distances).  To prevent collapse to a full ring we
  # use a two-stage strategy:
  #   Stage 1 – fix arc_half at the heuristic estimate; optimise center /
  #              normal / radius / t_center.  This anchors the arc geometry.
  #   Stage 2 – release arc_half but bound it to ±40 % of the heuristic value;
  #              fine-tune all parameters together.
  n <- nrow(X); if (n < 4) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)

  geo  <- fit_arc_geometry(X, arc_span_hint = arc_span_hint)
  c0   <- if (!is.null(center_init)) as.numeric(center_init) else geo$center
  nrm0 <- geo$normal
  th0  <- acos(max(-1, min(1, nrm0[3])))
  ph0  <- atan2(nrm0[2], nrm0[1])
  r0   <- geo$radius
  tc0  <- geo$t_center
  ah0  <- geo$arc_half          # heuristic half-span; used as anchor
  r_max <- max(sqrt(rowSums(sweep(X, 2, c0)^2))) * 2

  # ── Orientation grid search: best normal direction before L-BFGS-B ────────
  # Uses the arc projection with tc0/ah0 fixed — only orientation is swept.
  c_gs <- if (!is.null(center_fixed)) as.numeric(center_fixed) else c0
  gs_a <- .grid_search_orient(
    obj_fn  = function(th, ph)
      sum(w * mle_project_arc(X, c_gs, .ang_to_u(th, ph), r0, tc0, ah0)$dist^2) / n,
    th_init = th0, ph_init = ph0
  )
  th0 <- gs_a[1]; ph0 <- gs_a[2]

  # arc_half bounds: ±40 % of heuristic, clamped to [1°, 179°]
  ah_lo_log <- log(pmax(ah0 * 0.60, 0.02))
  ah_hi_log <- log(pmin(ah0 * 1.40, pi - 1e-6))

  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand

  # Shared objective: fixed arc_half supplied externally
  .obj_fixed_ah <- function(par, cc, ah) {
    pr  <- mle_project_arc(X, cc, .ang_to_u(par[1], par[2]), exp(par[3]), par[4], ah)
    sum(w * pr$dist^2) / n
  }
  # Full objective with arc_half as par[5]
  .obj_full_cc <- function(par, cc) {
    ah  <- pmin(exp(par[5]), pi - 1e-6)
    pr  <- mle_project_arc(X, cc, .ang_to_u(par[1], par[2]), exp(par[3]), par[4], ah)
    sum(w * pr$dist^2) / n
  }

  if (!is.null(center_fixed)) {
    cc <- as.numeric(center_fixed)
    # Stage 1
    p1 <- optim(c(th0, ph0, log(max(r0, 1e-3)), tc0),
                .obj_fixed_ah, cc = cc, ah = ah0, method = "L-BFGS-B",
                lower = c(-Inf, -Inf, log(1e-3), -Inf),
                upper = c( Inf,  Inf, log(r_max),  Inf),
                control = list(maxit = 300))$par
    # Stage 2
    p2 <- optim(c(p1, log(ah0)), .obj_full_cc, cc = cc, method = "L-BFGS-B",
                lower = c(-Inf, -Inf, log(1e-3), -Inf, ah_lo_log),
                upper = c( Inf,  Inf, log(r_max),  Inf, ah_hi_log),
                control = list(maxit = 300))$par
    return(list(center = cc, normal = .ang_to_u(p2[1], p2[2]), radius = exp(p2[3]),
                t_center = p2[4], arc_half = pmin(exp(p2[5]), pi - 1e-6)))
  }

  # Free centre ── Stage 1: fix arc_half
  obj1 <- function(par) {
    pr  <- mle_project_arc(X, par[1:3], .ang_to_u(par[4], par[5]),
                           exp(par[6]), par[7], ah0)
    val <- sum(w * pr$dist^2) / n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((par[1:3] - center_prior)^2)
    val
  }
  f1 <- optim(c(c0, th0, ph0, log(max(r0, 1e-3)), tc0), obj1,
              method = "L-BFGS-B",
              lower = c(c_lo, -Inf, -Inf, log(1e-3), -Inf),
              upper = c(c_hi,  Inf,  Inf, log(r_max),  Inf),
              control = list(maxit = 300))$par

  # Stage 2: release arc_half within ±40 % of heuristic
  obj2 <- function(par) {
    ah  <- pmin(exp(par[8]), pi - 1e-6)
    pr  <- mle_project_arc(X, par[1:3], .ang_to_u(par[4], par[5]),
                           exp(par[6]), par[7], ah)
    val <- sum(w * pr$dist^2) / n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((par[1:3] - center_prior)^2)
    val
  }
  f2 <- optim(c(f1, log(ah0)), obj2, method = "L-BFGS-B",
              lower = c(c_lo, -Inf, -Inf, log(1e-3), -Inf, ah_lo_log),
              upper = c(c_hi,  Inf,  Inf, log(r_max),  Inf, ah_hi_log),
              control = list(maxit = 300))$par

  list(center = f2[1:3], normal = .ang_to_u(f2[4], f2[5]), radius = exp(f2[6]),
       t_center = f2[7], arc_half = pmin(exp(f2[8]), pi - 1e-6))
}

mle_sample_arc <- function(fit, n = 300) {
  basis   <- .ring_basis(fit$normal)
  t_start <- fit$t_center - fit$arc_half
  t_end   <- fit$t_center + fit$arc_half
  tt      <- seq(t_start, t_end, length.out = n)
  pts     <- fit$radius * (outer(cos(tt), basis[, 1]) + outer(sin(tt), basis[, 2])) +
             matrix(fit$center, n, 3, TRUE)
  colnames(pts) <- c("x", "y", "z")
  as.data.frame(pts)
}

## ── Tube mesh helpers for helix, ring, and arc ─────────────────────────────
# Internal: given centerline positions and two cross-section frame vectors
# (frame1, frame2, each n_along × 3), build an n_along × n_tube vertex grid
# and return a mesh3d-compatible list (x, y, z, i, j, k — 0-indexed for plotly).
.tube_mesh_from_frames <- function(positions, frame1, frame2, tube_r,
                                   n_tube = 16, wrap_along = FALSE) {
  n_along <- nrow(positions)
  phi <- seq(0, 2*pi, length.out = n_tube + 1L)[-(n_tube + 1L)]

  vert <- matrix(NA_real_, n_along * n_tube, 3L)
  for (a in seq_len(n_along)) {
    base <- (a - 1L) * n_tube
    for (t in seq_len(n_tube)) {
      vert[base + t, ] <- positions[a, ] +
        tube_r * (cos(phi[t]) * frame1[a, ] + sin(phi[t]) * frame2[a, ])
    }
  }

  vidx <- function(a, t) (a - 1L) * n_tube + (t - 1L)   # 0-based for plotly
  n_a_max <- if (wrap_along) n_along else n_along - 1L
  ii <- integer(n_a_max * n_tube * 2L)
  jj <- integer(n_a_max * n_tube * 2L)
  kk <- integer(n_a_max * n_tube * 2L)
  k_out <- 1L
  for (a in seq_len(n_a_max)) {
    a2 <- if (a == n_along) 1L else a + 1L
    for (t in seq_len(n_tube)) {
      t2 <- if (t == n_tube) 1L else t + 1L
      ii[k_out]     <- vidx(a,  t);  jj[k_out]     <- vidx(a2, t);  kk[k_out]     <- vidx(a,  t2)
      ii[k_out + 1] <- vidx(a2, t);  jj[k_out + 1] <- vidx(a2, t2); kk[k_out + 1] <- vidx(a,  t2)
      k_out <- k_out + 2L
    }
  }

  list(x = vert[, 1], y = vert[, 2], z = vert[, 3], i = ii, j = jj, k = kk)
}

# Torus mesh for a ring fit (closed tube around a 3D circle).
make_ring_tube_mesh <- function(center, normal, radius, tube_r,
                                n_ring = 64, n_tube = 16) {
  basis <- .ring_basis(normal)
  u <- basis[, 1]; v <- basis[, 2]
  n_vec <- as.numeric(normal) / sqrt(sum(normal^2))
  theta  <- seq(0, 2*pi, length.out = n_ring + 1L)[-(n_ring + 1L)]
  cos_t  <- cos(theta); sin_t <- sin(theta)
  pos    <- matrix(center, n_ring, 3L, TRUE) +
            radius * (outer(cos_t, u) + outer(sin_t, v))
  radial <- outer(cos_t, u) + outer(sin_t, v)       # outward from ring axis
  axial  <- matrix(n_vec, n_ring, 3L, TRUE)
  .tube_mesh_from_frames(pos, radial, axial, tube_r, n_tube = n_tube, wrap_along = TRUE)
}

# Open-ended tube for an arc fit.
make_arc_tube_mesh <- function(center, normal, radius, t_center, arc_half, tube_r,
                               n_arc = 64, n_tube = 16) {
  basis <- .ring_basis(normal)
  u <- basis[, 1]; v <- basis[, 2]
  n_vec  <- as.numeric(normal) / sqrt(sum(normal^2))
  theta  <- seq(t_center - arc_half, t_center + arc_half, length.out = n_arc)
  cos_t  <- cos(theta); sin_t <- sin(theta)
  pos    <- matrix(center, n_arc, 3L, TRUE) +
            radius * (outer(cos_t, u) + outer(sin_t, v))
  radial <- outer(cos_t, u) + outer(sin_t, v)
  axial  <- matrix(n_vec, n_arc, 3L, TRUE)
  .tube_mesh_from_frames(pos, radial, axial, tube_r, n_tube = n_tube, wrap_along = FALSE)
}

# Open-ended tube for a helix fit.
# center/axis/radius/pitch/phi0 are the helix parameters as used by .helix_point.
# tmin/tmax are the parameter range to render.
make_helix_tube_mesh <- function(center, axis, radius, pitch, phi0, tmin, tmax,
                                 tube_r, n_helix = 200, n_tube = 12) {
  tt  <- seq(tmin, tmax, length.out = n_helix)
  pos <- t(vapply(tt, .helix_point, numeric(3),
                  c = center, u = axis, r = radius, p = pitch, phi0 = phi0))
  ax  <- axis / sqrt(sum(axis^2))

  frame1 <- matrix(NA_real_, n_helix, 3L)
  frame2 <- matrix(NA_real_, n_helix, 3L)
  for (i in seq_len(n_helix)) {
    # Outward radial: component of (pos - center) perpendicular to axis
    diff_i   <- pos[i, ] - center
    p_axial  <- sum(diff_i * ax) * ax
    p_radial <- diff_i - p_axial
    r_len    <- sqrt(sum(p_radial^2))
    f1 <- if (r_len > 1e-12) p_radial / r_len else {
      ref <- if (abs(ax[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0)
      tmp <- pracma::cross(ax, ref); tmp / sqrt(sum(tmp^2))
    }
    frame1[i, ] <- f1
    frame2[i, ] <- pracma::cross(ax, f1)
  }
  .tube_mesh_from_frames(pos, frame1, frame2, tube_r, n_tube = n_tube, wrap_along = FALSE)
}

# Robust meshing from scattered surface points
add_alpha_or_hull_mesh <- function(fig, pts, color, opacity, name, alpha) {
  pts <- as.matrix(pts)
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]
  pts <- unique(pts)
  if (nrow(pts) < 4) return(fig)
  
  # 1) Try alpha-shape
  ashape <- tryCatch(alphashape3d::ashape3d(pts, alpha = alpha), error = function(e) NULL)
  if (!is.null(ashape) && !is.null(ashape$triang) && nrow(ashape$triang) > 0) {
    tri <- as.matrix(ashape$triang[, 1:3, drop = FALSE]) - 1L
    return(
      fig %>% plotly::add_trace(
        type = "mesh3d",
        x = pts[,1], y = pts[,2], z = pts[,3],
        i = tri[,1], j = tri[,2], k = tri[,3],
        opacity = opacity,
        color = I(color), 
        name = name, showlegend = TRUE
      )
    )
  }
  
  # 2) Fallback: convex hull
  hull <- tryCatch(geometry::convhulln(pts, options = "Qt"), error = function(e) NULL)
  if (!is.null(hull) && nrow(hull) > 0 && ncol(hull) >= 3) {
    tri <- as.matrix(hull[, 1:3, drop = FALSE]) - 1L
    fig <- fig %>% plotly::add_trace(
      type = "mesh3d",
      x = pts[,1], y = pts[,2], z = pts[,3],
      i = tri[,1], j = tri[,2], k = tri[,3],
      opacity = opacity,
      color = I(color), 
      name = paste0(name, " (hull)"),
      showlegend = TRUE
    )
  }
  fig
}
add_distance_lines_to_plot <- function(fig, medians, max_dist = Inf) {
  required_cols <- c("x", "y", "z", "nearest_x", "nearest_y", "nearest_z", "distance_to_fit")
  if (!all(required_cols %in% colnames(medians))) return(fig)

  # Drop rows where either endpoint is missing/non-finite; apply length threshold
  valid <- is.finite(medians$x) & is.finite(medians$y) & is.finite(medians$z) &
           is.finite(medians$nearest_x) & is.finite(medians$nearest_y) & is.finite(medians$nearest_z) &
           is.finite(medians$distance_to_fit)
  if (is.finite(max_dist) && max_dist > 0) valid <- valid & (medians$distance_to_fit <= max_dist)
  medians <- medians[valid, , drop = FALSE]
  if (nrow(medians) == 0) return(fig)

  # Build NA-separated x/y/z for a single line trace (each triplet = one segment)
  xs <- as.vector(rbind(medians$x, medians$nearest_x, NA_real_))
  ys <- as.vector(rbind(medians$y, medians$nearest_y, NA_real_))
  zs <- as.vector(rbind(medians$z, medians$nearest_z, NA_real_))

  fig <- fig %>% add_trace(
    x = xs, y = ys, z = zs,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "blue", width = 2),
    showlegend = FALSE,
    inherit = FALSE
  )

  # Text labels at midpoints — one trace for all labels
  fig <- fig %>% add_trace(
    x = (medians$x + medians$nearest_x) / 2,
    y = (medians$y + medians$nearest_y) / 2,
    z = (medians$z + medians$nearest_z) / 2,
    type = "scatter3d",
    mode = "text",
    text = as.character(round(medians$distance_to_fit, 1)),
    textfont = list(color = "blue", size = 10),
    showlegend = FALSE,
    inherit = FALSE
  )

  return(fig)
}

validate_coords <- function(mat) {
  mat[!duplicated(mat), , drop = FALSE]
}

compute_convex_hull_centroid <- function(coords_matrix) {
  if (nrow(coords_matrix) < 4 || ncol(coords_matrix) != 3) {
    warning("Invalid input to convex hull centroid computation.")
    return(colMeans(coords_matrix))
  }
  
  result <- tryCatch({
    hull <- geometry::convhulln(coords_matrix, options = "FA")
    if (is.null(hull$hull)) stop("Hull computation failed")
    
    triangles <- hull$hull
    pts <- coords_matrix
    
    tri_centers <- matrix(0, nrow = nrow(triangles), ncol = 3)
    tri_areas <- numeric(nrow(triangles))
    
    for (i in seq_len(nrow(triangles))) {
      tri_idx <- triangles[i, ]
      if (length(tri_idx) != 3) next  # skip invalid
      tri_pts <- pts[tri_idx, , drop = FALSE]
      if (nrow(tri_pts) != 3 || ncol(tri_pts) != 3) next  # skip if malformed
      
      v1 <- tri_pts[2, ] - tri_pts[1, ]
      v2 <- tri_pts[3, ] - tri_pts[1, ]
      if (length(v1) != 3 || length(v2) != 3) next  # safety
      
      centroid <- colMeans(tri_pts)
      area <- 0.5 * sqrt(sum(crossprod(pracma::cross(v1, v2))))
      
      tri_centers[i, ] <- centroid
      tri_areas[i] <- area
    }
    
    if (sum(tri_areas) == 0) stop("Degenerate triangles")
    
    center <- colSums(tri_centers * tri_areas) / sum(tri_areas)
    center
  }, error = function(e) {
    warning("Convex hull fallback: ", conditionMessage(e))
    return(colMeans(coords_matrix))  # fallback
  })
  
  return(result)
}
medianize_coords <- function(df, use_cluster1 = FALSE) {
  if (use_cluster1 &&
      "cluster" %in% names(df) &&
      any(is.finite(df$cluster) & df$cluster != 0, na.rm = TRUE)) {
    df %>%
      dplyr::filter(is.finite(cluster) & cluster != 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(x = median(x), y = median(y), z = median(z), .groups = "drop") %>%
      dplyr::select(x, y, z)
  } else {
    df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = median(x), y = median(y), z = median(z), .groups = "drop") %>%
      dplyr::select(x, y, z)
  }
}
# ---- lightweight geometry helpers (parameters only, no random sampling) ----
spread_from_group <- function(x, y, z) {
  if (length(x) <= 1) return(0)
  mx <- stats::median(x, na.rm = TRUE)
  my <- stats::median(y, na.rm = TRUE)
  mz <- stats::median(z, na.rm = TRUE)
  sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2, na.rm = TRUE))
}

fit_sphere_geometry <- function(coords_matrix) {
  center <- tryCatch(compute_convex_hull_centroid(coords_matrix), error = function(e) colMeans(coords_matrix))
  r      <- stats::median(sqrt(rowSums(sweep(coords_matrix, 2, center)^2)))
  list(center = center, radius = r)
}

fit_ellipsoid_geometry <- function(coords_matrix) {
  center <- tryCatch(compute_convex_hull_centroid(coords_matrix), error = function(e) colMeans(coords_matrix))
  Xc     <- sweep(coords_matrix, 2, center)   # center at hull centroid before PCA so axes are aligned to ellipsoid center
  pca    <- prcomp(Xc, center = FALSE, scale. = FALSE)
  rot    <- Xc %*% pca$rotation
  d <- diff(stats::quantile(rot[, 1], c(0.05, 0.95))) / 2
  e <- diff(stats::quantile(rot[, 2], c(0.05, 0.95))) / 2
  f <- diff(stats::quantile(rot[, 3], c(0.05, 0.95))) / 2
  sf <- stats::quantile(sqrt(rowSums(rot^2)), 0.90) / max(d, e, f)
  d <- d * sf; e <- e * sf; f <- f * sf
  if (min(d, e, f) < 1e-6) return(NULL)
  list(center = center, axes = c(d, e, f), R = pca$rotation)
}

fit_helix_geometry <- function(coords_matrix, turns = 3, n_points = 5000) {
  center <- colMeans(coords_matrix)
  Xc     <- sweep(coords_matrix, 2, center)
  pca    <- prcomp(Xc, center = FALSE)
  a1 <- pca$rotation[, 1]; a2 <- pca$rotation[, 2]
  a3 <- pracma::cross(a1, a2)
  ht <- diff(range(Xc %*% a1))
  br <- stats::quantile(abs(Xc %*% a2), 0.9, na.rm = TRUE)
  t  <- seq(0, 2 * pi * max(1, turns), length.out = n_points)
  Hloc <- cbind(br * cos(t), br * sin(t),
                seq(-ht / 2, ht / 2, length.out = n_points))
  H <- sweep(Hloc %*% t(cbind(a2, a3, a1)), 2, center, "+")
  colnames(H) <- c("x", "y", "z")

  # Parameters compatible with .project_t_helix / .helix_point so that
  # Newton projection operates on the exact same 3D curve as H above.
  # .rot_to_u(a1) uses a different x/y-axis convention than cbind(a2,a3,a1),
  # so phi0 aligns the phase; proj_center shifts z-origin from -ht/2 to 0.
  #
  # Derivation: .helix_point(t=0) = proj_center + r*(cos(phi0)*R[,1]+sin(phi0)*R[,2])
  # must equal H[1,] = proj_center + r*a2.
  # Projecting a2 onto the {R[,1], R[,2]} Newton frame gives:
  #   cos(phi0) = a2 . R[,1]  →  sin(phi0) = a2 . R[,2]
  R_nwt       <- .rot_to_u(a1)
  phi0        <- atan2(sum(a2 * R_nwt[, 2]), sum(a2 * R_nwt[, 1]))
  proj_center <- center - a1 * (ht / 2)

  list(points = H, center = center,
       axis = a1, radius = br, pitch = ht / max(1, turns),
       proj_center = proj_center, phi0 = phi0)
}

fit_ring_geometry <- function(coords_matrix) {
  center <- tryCatch(compute_convex_hull_centroid(coords_matrix),
                     error = function(e) colMeans(coords_matrix))
  Xc     <- sweep(coords_matrix, 2, center)
  pca    <- prcomp(Xc, center = FALSE)
  # PC3 = direction of smallest variance = normal to the ring plane
  normal <- pca$rotation[, 3]

  d_ax   <- as.numeric(Xc %*% normal)
  Xp     <- Xc - outer(d_ax, normal)      # in-plane projections
  radius <- stats::median(sqrt(rowSums(Xp^2)))
  list(center = center, normal = normal, radius = radius)
}

fit_arc_geometry <- function(coords_matrix, arc_span_hint = 0) {
  geo   <- fit_ring_geometry(coords_matrix)
  basis <- .ring_basis(geo$normal)
  u <- basis[, 1]; vv <- basis[, 2]

  Xc    <- sweep(coords_matrix, 2, geo$center)
  d_ax  <- as.numeric(Xc %*% geo$normal)
  Xp    <- Xc - outer(d_ax, geo$normal)
  t_ang <- atan2(as.numeric(Xp %*% vv), as.numeric(Xp %*% u))

  # Circular mean → arc midpoint angle
  t_center <- atan2(mean(sin(t_ang)), mean(cos(t_ang)))

  arc_half <- if (arc_span_hint > 1) {
    pmin((arc_span_hint / 2) * pi / 180, pi - 1e-6)
  } else {
    # Auto: 97th-percentile angular deviation from the arc midpoint + 5% buffer
    dt <- ((t_ang - t_center + pi) %% (2*pi)) - pi
    pmin(stats::quantile(abs(dt), 0.97, names = FALSE) * 1.05, pi - 1e-6)
  }
  list(center = geo$center, normal = geo$normal, radius = geo$radius,
       t_center = t_center, arc_half = arc_half)
}

generate_ring_from_cluster2 <- function(cluster_data, n_points = 500,
                                        use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  geo <- fit_ring_geometry(as.matrix(coords))
  pts <- mle_sample_ring(geo, n = n_points)
  list(points = pts, center = geo$center, normal = geo$normal, radius = geo$radius)
}

generate_arc_from_cluster2 <- function(cluster_data, n_points = 300,
                                       arc_span_hint = 0, use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  geo <- fit_arc_geometry(as.matrix(coords), arc_span_hint = arc_span_hint)
  pts <- mle_sample_arc(geo, n = n_points)
  list(points = pts, center = geo$center, normal = geo$normal, radius = geo$radius,
       t_center = geo$t_center, arc_half = geo$arc_half)
}

# =============================================================================
# alpha_surface_3d() — WRAP ONE GROUP OF LOCALIZATIONS IN AN ALPHA SHAPE,
#                      RADIUS CHOSEN FROM THE DATA
# =============================================================================
#
# GOAL
#   Given the localizations of ONE group (a trace, a de-blinked molecule, a
#   Median Cluster), return the closed triangular surface that wraps it, ready
#   to hand straight to plotly's mesh3d. Callers only supply points; there is
#   deliberately no radius argument and no UI control for one.
#
# WHY THERE IS NO ALPHA RADIUS CONTROL ANY MORE (this was a real bug)
#   `alphashape3d::ashape3d()` computes ONE Delaunay tetrahedralization and then,
#   for each requested alpha, adds a COLUMN classifying every triangle. The
#   triangle list in `$triang[, 1:3]` is therefore identical for every alpha —
#   alpha lives only in the classification columns. The old code drew
#   `$triang[, 1:3]` wholesale, i.e. the entire tetrahedralization, so the mesh
#   was byte-for-byte the same at every slider position. That is why moving the
#   "Alpha radius" slider visibly did nothing.
#
#   Column layout of `$triang` (verified against alphashape3d):
#     1-3  vertex indices        4  on.ch (on convex hull)   5  attached
#     6-8  rhoT / muT / MuT      9..  one class per alpha:
#          0 = exterior, 1 = singular, 2 = REGULAR, 3 = interior
#   Only REGULAR triangles form the closed boundary surface. (Check: as alpha
#   grows the regular set converges exactly onto the convex-hull faces, while
#   the singular count explodes — singular triangles are dangling slivers and
#   drawing them is what produced the old opaque blob.)
#
# HOW THE RADIUS IS PICKED
#   A ladder of candidate radii is built from the group's OWN spacing (multiples
#   of its median nearest-neighbour distance) and classified in a SINGLE
#   ashape3d call — extra alphas cost one column each, not another
#   tetrahedralization. We then take the SMALLEST radius that passes two tests:
#     (1) every point sits in ONE connected component (nothing shed), and
#     (2) the regular triangles form a CLOSED surface — each edge used by
#         exactly two of them. Test (2) matters because connectivity alone still
#         admits radii whose boundary has holes, and plotly renders a holed mesh
#         as a torn, see-through shell.
#   That is the tightest honest wrap: smaller and the surface fragments, larger
#   and detail is lost for nothing. It is scale-free, so a 30 nm de-blinked
#   molecule and a 500 nm structure are each wrapped correctly without anyone
#   touching a control. Falls back to the convex hull (largest rung, then
#   geometry::convhulln) if no candidate ever qualifies.
#
# ARGUMENTS
#   pts_mat : n x 3 numeric matrix of coordinates (nm)
#   max_pts : subsample cap. Meshing cost grows steeply with n and a surface is
#             not made more truthful by tens of thousands of vertices.
# RETURNS
#   list(x, y, z, i, j, k, alpha, n_tri, n_pts) with i/j/k already 0-BASED for
#   plotly mesh3d and x/y/z the (possibly subsampled) vertices they index, or
#   NULL when the group cannot be wrapped (<4 usable points, or no triangle
#   survives). Callers MUST use the returned x/y/z, not their own copy of the
#   points, because of the subsampling.
# =============================================================================
alpha_surface_3d <- function(pts_mat, max_pts = 2000L) {
  if (!is.matrix(pts_mat)) pts_mat <- as.matrix(pts_mat)     # accept data.frame
  if (ncol(pts_mat) < 3L) return(NULL)
  pts_mat <- pts_mat[, 1:3, drop = FALSE]
  storage.mode(pts_mat) <- "double"
  pts_mat <- pts_mat[stats::complete.cases(pts_mat), , drop = FALSE]
  if (nrow(pts_mat) < 4L) return(NULL)                        # a tetrahedron is the minimum
  # Deterministic subsample (sorted) so the mesh does not flicker between renders.
  if (nrow(pts_mat) > max_pts)
    pts_mat <- pts_mat[sort(sample.int(nrow(pts_mat), max_pts)), , drop = FALSE]

  # A group that is flat in any axis (2D data promoted to 3D, or a trace whose z
  # never moved) has zero volume and makes the tetrahedralization degenerate.
  # Nudge by sub-picometre amounts: invisible at nm scale, enough for the solver.
  for (k in 1:3) {
    sdk <- stats::sd(pts_mat[, k], na.rm = TRUE)
    if (!is.finite(sdk) || sdk < 1e-9)
      pts_mat[, k] <- pts_mat[, k] + seq_len(nrow(pts_mat)) * 1e-6
  }

  # --- candidate radii, expressed in the group's own units of spacing ---------
  nnd <- tryCatch({
    d <- RANN::nn2(pts_mat, pts_mat, k = 2L)$nn.dists[, 2]    # distance to nearest OTHER point
    d[is.finite(d) & d > 0]
  }, error = function(e) numeric(0))
  # Bounding-box diagonal bounds the useful range: an alpha far above it can only
  # ever give the convex hull, so it is the natural top rung.
  diag_nm <- sqrt(sum((apply(pts_mat, 2, max) - apply(pts_mat, 2, min))^2))
  base <- if (length(nnd)) stats::median(nnd) else diag_nm / 10
  if (!is.finite(base) || base <= 0) base <- diag_nm / 10
  if (!is.finite(base) || base <= 0) return(NULL)             # all points coincide
  alphas <- base * c(1, 1.5, 2, 3, 4, 6, 8, 12, 20, 40)
  alphas <- c(alphas[alphas < diag_nm], 10 * diag_nm)         # last rung = convex hull
  alphas <- sort(unique(alphas[is.finite(alphas) & alphas > 0]))
  if (!length(alphas)) return(NULL)

  # Degenerate inputs (n == 4, all-coincident points) make ashape3d emit internal
  # warnings on its way to returning nothing useful; they are not actionable.
  as3d <- suppressWarnings(tryCatch(alphashape3d::ashape3d(pts_mat, alpha = alphas),
                                    error = function(e) NULL))

  # Is this triangle set a closed surface? Every edge of a closed, watertight
  # mesh belongs to exactly two faces; an edge used once is a hole rim.
  is_closed <- function(tt) {
    if (is.null(tt) || nrow(tt) < 4L) return(FALSE)
    e <- rbind(tt[, c(1, 2)], tt[, c(2, 3)], tt[, c(1, 3)])
    e <- cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2]))    # undirected edge key
    all(tabulate(match(paste(e[, 1], e[, 2]), unique(paste(e[, 1], e[, 2])))) == 2L)
  }

  tt <- NULL; alpha_used <- NA_real_
  if (!is.null(as3d) && !is.null(as3d$triang) && nrow(as3d$triang)) {
    tri  <- as.matrix(as3d$triang)
    # The per-alpha classification columns are the LAST length(alphas) columns, in
    # the order the alphas were passed. Indexing from the end keeps this correct
    # if the package ever adds a descriptive column at the front.
    cls0 <- ncol(tri) - length(alphas)
    # components_ashape3d() labels each point with its connected component
    # (1 = largest, < 0 for points isolated at this radius). One component
    # covering every point == "wraps the whole group in one piece".
    comp <- suppressWarnings(tryCatch(
      alphashape3d::components_ashape3d(as3d, indexAlpha = "all"),
      error = function(e) NULL))
    if (!is.null(comp) && !is.list(comp)) comp <- list(comp)
    if (cls0 >= 3L) {
      for (ii in seq_along(alphas)) {
        cc <- if (!is.null(comp) && length(comp) == length(alphas)) comp[[ii]] else NULL
        # Test (1): nothing shed. Skipped only if components could not be computed.
        if (!is.null(cc) && !(length(cc) && all(is.finite(cc)) && all(cc == 1L))) next
        cand <- tri[tri[, cls0 + ii] == 2L, 1:3, drop = FALSE]  # REGULAR = the boundary
        if (!nrow(cand)) next
        # Test (2): watertight. The last rung is the convex hull, always accepted.
        if (is_closed(cand) || ii == length(alphas)) {
          tt <- cand; alpha_used <- alphas[ii]; break
        }
      }
    }
  }
  # Last resort: the convex hull. Correct (if loose) for any point set, and the
  # only thing that works for the n == 4 / near-coincident cases ashape3d rejects.
  if (is.null(tt) || !nrow(tt)) {
    tt <- tryCatch(as.matrix(geometry::convhulln(pts_mat, options = "Qt")[, 1:3, drop = FALSE]),
                   error = function(e) NULL)
    alpha_used <- Inf
  }
  if (is.null(tt) || !nrow(tt)) return(NULL)

  list(x = pts_mat[, 1], y = pts_mat[, 2], z = pts_mat[, 3],
       i = tt[, 1] - 1L, j = tt[, 2] - 1L, k = tt[, 3] - 1L,  # plotly mesh3d is 0-based
       alpha = alpha_used, n_tri = nrow(tt), n_pts = nrow(pts_mat))
}

# NOTE: `virtual_points` is a dense uniform sampling of the fitted SPHERE, not
# data. Its triangulation is the convex hull by construction, so no alpha shape
# is built here — the caller hulls it. (An `ashape` element used to be returned
# alongside; it was never actually drawn.)
generate_sphere_from_cluster2 <- function(cluster_data, n_points = 4000,
                                          use_cluster1 = FALSE) {

  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  coords_matrix <- as.matrix(coords)

  geo <- fit_sphere_geometry(coords_matrix)

  phi <- runif(n_points, 0, 2*pi)
  theta <- acos(runif(n_points, -1, 1))
  virtual_points <- cbind(
    geo$radius * sin(theta) * cos(phi),
    geo$radius * sin(theta) * sin(phi),
    geo$radius * cos(theta)
  ) + matrix(geo$center, n_points, 3, TRUE)

  list(center = geo$center, radius = geo$radius, virtual_points = virtual_points)
}

# Same as generate_sphere_from_cluster2: a dense sampling of the fitted (convex)
# ELLIPSOID, hulled by the caller. No alpha shape involved.
generate_ellipsoid_from_cluster2 <- function(cluster_data, n_points = 4000,
                                             use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  coords_matrix <- as.matrix(coords)

  geo <- fit_ellipsoid_geometry(coords_matrix)
  if (is.null(geo)) return(NULL)

  phi <- runif(n_points, 0, 2*pi)
  theta <- acos(runif(n_points, -1, 1))
  ellipsoid_points <- cbind(
    geo$axes[1] * sin(theta) * cos(phi),
    geo$axes[2] * sin(theta) * sin(phi),
    geo$axes[3] * cos(theta)
  ) %*% t(geo$R) + matrix(geo$center, n_points, 3, TRUE)

  list(center = geo$center, axes = geo$axes, R = geo$R,
       virtual_points = ellipsoid_points)
}
generate_helix_from_cluster2 <- function(cluster_data, n_points = 10000, turns = 3,
                                         use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 5) return(NULL)
  coords_matrix <- as.matrix(coords)

  geo <- fit_helix_geometry(coords_matrix, turns = turns, n_points = n_points)
  list(points = as.data.frame(geo$points), center = geo$center)
}


add_cluster2_lines_and_angles <- function(fig, df, plane,
                                          use_pca_normal = TRUE,
                                          cluster2_table_store = NULL,
                                          helix_turns = 3,
                                          show_plane_angles = FALSE,
                                          show_plane_distance = FALSE,
                                          show_pca_axes = FALSE,
                                          enable_clustering = FALSE,
                                          dark_mode = FALSE,
                                          bending_radius = FALSE,
                                          show_info_labels = FALSE,
                                          compute_rmsd_metrics = FALSE,
                                          selected_fit_type = "sphere",
                                          selected_fit_method = "heuristic",
                                          fit_weighted = FALSE,
                                          mle_center_mode = "free",
                                          mle_center_lambda = 0.1,
                                          arc_span_hint = 0) {
  
  df <- df %>% dplyr::filter(!is.na(cluster2) & cluster2 != 0)
  if (nrow(df) == 0) return(list(fig = fig, log = "No Median Cluster points found."))
  
  log_msgs <- character(0)
  result_table <- data.frame()
  clusters <- df %>% dplyr::group_by(cluster2) %>% dplyr::group_split()
  
  for (cluster_data in clusters) {

  cluster_id <- unique(cluster_data$cluster2)[1]
  
  fmt <- function(x, digits = 2) {
    ifelse(is.na(x) || !is.finite(x), "NA",
           format(round(x, digits), nsmall = digits, trim = TRUE))
  }

  # --- representative points for THIS cluster2 (used for ALL metrics/fit geometry) ---
  has_valid_clusters <- "cluster" %in% names(cluster_data) &&
    any(!is.na(cluster_data$cluster), na.rm = TRUE)

  if (isTRUE(enable_clustering) && has_valid_clusters) {
    medians_df <- cluster_data %>%
      dplyr::filter(!is.na(cluster) & cluster != 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        x        = stats::median(x, na.rm = TRUE),
        y        = stats::median(y, na.rm = TRUE),
        z        = stats::median(z, na.rm = TRUE),
        spread   = spread_from_group(x, y, z),
        n_locs   = dplyr::n(),
        group_id = dplyr::first(cluster),
        .groups  = "drop"
      )
  } else {
    medians_df <- cluster_data %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        x        = stats::median(x, na.rm = TRUE),
        y        = stats::median(y, na.rm = TRUE),
        z        = stats::median(z, na.rm = TRUE),
        spread   = spread_from_group(x, y, z),
        n_locs   = dplyr::n(),
        group_id = dplyr::first(tid),
        .groups  = "drop"
      )
  }

  medians_df <- medians_df %>% dplyr::filter(is.finite(x), is.finite(y), is.finite(z))

  if (nrow(medians_df) < 4) {
    log_msgs <- c(log_msgs, paste("Median Cluster", cluster_id, "- Not enough representative medians (<4)."))
    next
  }
  need_rmsd <- isTRUE(compute_rmsd_metrics)
  
  # ---- representative point matrix for this cluster2 ----
  coords_matrix <- as.matrix(medians_df[, c("x","y","z"), drop = FALSE])
  
  # Weights for MLE fits: proportional to n_locs / spread² (≈ 1/SEM² of the
  # median position).  Traces with a single localization have spread = 0 by
  # construction (variance is undefined for n = 1); they receive the median
  # spread of multi-localization traces as a conservative fallback so they
  # participate at roughly average weight rather than near-infinite weight.
  # RMSD metrics are always reported unweighted.
  w <- if (isTRUE(fit_weighted) && "spread" %in% names(medians_df)) {
    s   <- medians_df$spread
    n_l <- if ("n_locs" %in% names(medians_df)) medians_df$n_locs else rep(1L, nrow(medians_df))
    s_ref <- stats::median(s[n_l > 1 & s > 0], na.rm = TRUE)
    if (!is.finite(s_ref) || s_ref <= 0) s_ref <- max(s[s > 0], na.rm = TRUE)
    if (!is.finite(s_ref) || s_ref <= 0) s_ref <- 1
    s_eff <- ifelse(n_l > 1 & s > 0, s, s_ref)
    n_l / s_eff^2
  } else NULL
  
  # defaults (fast path)
  rmsd_h_sphere <- meandist_h_sphere <- sd_h_sphere <- NA_real_
  rmsd_h_ellipsoid <- meandist_h_ellipsoid <- sd_h_ellipsoid <- NA_real_
  rmsd_h_helix <- meandist_h_helix <- sd_h_helix <- NA_real_
  best_fit_name <- NA_character_

  rmsd_m_sphere <- meandist_m_sphere <- sd_m_sphere <- NA_real_
  rmsd_m_ellipsoid <- meandist_m_ellipsoid <- sd_m_ellipsoid <- NA_real_
  rmsd_m_helix <- meandist_m_helix <- sd_m_helix <- NA_real_
  rmsd_h_ring  <- meandist_h_ring  <- sd_h_ring  <- NA_real_
  rmsd_h_arc   <- meandist_h_arc   <- sd_h_arc   <- NA_real_
  rmsd_m_ring  <- meandist_m_ring  <- sd_m_ring  <- NA_real_
  rmsd_m_arc   <- meandist_m_arc   <- sd_m_arc   <- NA_real_

  # ── Active model RMSD: always computed (regardless of compute_rmsd_metrics) ─
  active_rmsd <- active_meandist <- active_sd <- NA_real_
  if (!is.null(selected_fit_type) && selected_fit_type != "none") {
    .active_ft <- tryCatch(
      calculate_distance_to_fit(
        coords_matrix = coords_matrix,
        fit_type      = selected_fit_type,
        turns         = helix_turns,
        arc_span_hint = arc_span_hint
      ),
      error = function(e) NULL
    )
    if (!is.null(.active_ft)) {
      active_rmsd     <- attr(.active_ft, paste0("rmsd_",     selected_fit_type))
      active_meandist <- attr(.active_ft, paste0("meandist_", selected_fit_type))
      active_sd       <- attr(.active_ft, paste0("sd_",       selected_fit_type))
    }
    # MLE path: override with the MLE projection for the selected model only
    if (identical(selected_fit_method, "mle")) {
      center_args_act <- build_center_args(
        X = coords_matrix,
        input = list(mle_center_mode = mle_center_mode, mle_center_lambda = mle_center_lambda)
      )
      mle_act <- tryCatch(
        switch(selected_fit_type,
          sphere    = do.call(mle_fit_sphere,
                              c(list(X = coords_matrix, w = w), center_args_act)),
          ellipsoid = do.call(mle_fit_ellipsoid,
                              c(list(X = coords_matrix, w = w), center_args_act)),
          helix     = do.call(mle_fit_helix,
                              c(list(X = coords_matrix, w = w, turns_hint = helix_turns), center_args_act)),
          ring      = do.call(mle_fit_ring,
                              c(list(X = coords_matrix, w = w), center_args_act)),
          arc       = do.call(mle_fit_arc,
                              c(list(X = coords_matrix, w = w, arc_span_hint = arc_span_hint), center_args_act)),
          NULL
        ),
        error = function(e) NULL
      )
      if (!is.null(mle_act)) {
        d_act <- switch(selected_fit_type,
          sphere    = as.numeric(mle_project_sphere(coords_matrix, mle_act$center, mle_act$radius)$dist),
          ellipsoid = as.numeric(mle_project_ellipsoid(coords_matrix, mle_act$center, mle_act$axes, mle_act$R)$dist),
          helix = {
            if (!is.null(mle_act$t)) {
              t_a <- .project_t_helix(coords_matrix, mle_act$center, mle_act$axis,
                                       mle_act$radius, mle_act$pitch, mle_act$phi0,
                                       t0 = mle_act$t, check_neighbors = TRUE)
              H_a <- t(vapply(t_a, .helix_point, numeric(3),
                               c = mle_act$center, u = mle_act$axis,
                               r = mle_act$radius, p = mle_act$pitch, phi0 = mle_act$phi0))
              sqrt(rowSums((coords_matrix - H_a)^2))
            } else numeric(0)
          },
          ring = as.numeric(mle_project_ring(coords_matrix, mle_act$center, mle_act$normal, mle_act$radius)$dist),
          arc  = as.numeric(mle_project_arc(coords_matrix, mle_act$center, mle_act$normal,
                                             mle_act$radius, mle_act$t_center, mle_act$arc_half)$dist),
          numeric(0)
        )
        d_act <- d_act[is.finite(d_act)]
        if (length(d_act)) {
          active_rmsd     <- sqrt(mean(d_act^2))
          active_meandist <- mean(d_act)
          active_sd       <- if (length(d_act) > 1) stats::sd(d_act) else NA_real_
        }
      }
    }
  }

  if (need_rmsd) {
    # ---- Heuristic RMSD/SD (passes pre-computed coords_matrix; no re-medianing) ----
    fit_result <- calculate_distance_to_fit(
      coords_matrix = coords_matrix,
      fit_type      = "best",
      turns         = helix_turns
    )

    if (!is.null(fit_result)) {
      rmsd_h_sphere    <- attr(fit_result, "rmsd_sphere");        sd_h_sphere    <- attr(fit_result, "sd_sphere")
      meandist_h_sphere    <- attr(fit_result, "meandist_sphere")
      rmsd_h_ellipsoid <- attr(fit_result, "rmsd_ellipsoid");    sd_h_ellipsoid <- attr(fit_result, "sd_ellipsoid")
      meandist_h_ellipsoid <- attr(fit_result, "meandist_ellipsoid")
      rmsd_h_helix     <- attr(fit_result, "rmsd_helix");        sd_h_helix     <- attr(fit_result, "sd_helix")
      meandist_h_helix     <- attr(fit_result, "meandist_helix")
      rmsd_h_ring      <- attr(fit_result, "rmsd_ring");         sd_h_ring      <- attr(fit_result, "sd_ring")
      meandist_h_ring      <- attr(fit_result, "meandist_ring")
      rmsd_h_arc       <- attr(fit_result, "rmsd_arc");          sd_h_arc       <- attr(fit_result, "sd_arc")
      meandist_h_arc       <- attr(fit_result, "meandist_arc")
    }
    
    # ---- MLE RMSD/SD (only when MLE is the selected method) ----
    if (identical(selected_fit_method, "mle")) {
      center_args <- build_center_args(
        X = coords_matrix,
        input = list(mle_center_mode = mle_center_mode,
                     mle_center_lambda = mle_center_lambda)
      )

      fS <- tryCatch(do.call(mle_fit_sphere, c(list(X = coords_matrix, w = w), center_args)), error = function(e) NULL)
      if (!is.null(fS)) {
        pr <- mle_project_sphere(coords_matrix, fS$center, fS$radius)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_sphere     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_sphere <- if (length(d)) mean(d)         else NA_real_
        sd_m_sphere       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fE <- tryCatch(do.call(mle_fit_ellipsoid, c(list(X = coords_matrix, w = w), center_args)), error = function(e) NULL)
      if (!is.null(fE)) {
        pr <- mle_project_ellipsoid(coords_matrix, fE$center, fE$axes, fE$R)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_ellipsoid     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_ellipsoid <- if (length(d)) mean(d)         else NA_real_
        sd_m_ellipsoid       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fH <- tryCatch(do.call(mle_fit_helix, c(list(X = coords_matrix, w = w, turns_hint = helix_turns), center_args)), error = function(e) NULL)
      if (!is.null(fH) && !is.null(fH$t)) {
        t_rmsd <- .project_t_helix(coords_matrix, fH$center, fH$axis, fH$radius, fH$pitch, fH$phi0,
                                   t0 = fH$t, check_neighbors = TRUE)
        H_rmsd <- t(vapply(t_rmsd, .helix_point, numeric(3),
                           c = fH$center, u = fH$axis, r = fH$radius, p = fH$pitch, phi0 = fH$phi0))
        d <- sqrt(rowSums((coords_matrix - H_rmsd)^2)); d <- d[is.finite(d)]
        rmsd_m_helix     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_helix <- if (length(d)) mean(d)         else NA_real_
        sd_m_helix       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fR <- tryCatch(do.call(mle_fit_ring, c(list(X = coords_matrix, w = w), center_args)), error = function(e) NULL)
      if (!is.null(fR)) {
        pr <- mle_project_ring(coords_matrix, fR$center, fR$normal, fR$radius)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_ring     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_ring <- if (length(d)) mean(d)         else NA_real_
        sd_m_ring       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fA <- tryCatch(do.call(mle_fit_arc, c(list(X = coords_matrix, w = w, arc_span_hint = arc_span_hint), center_args)), error = function(e) NULL)
      if (!is.null(fA)) {
        pr <- mle_project_arc(coords_matrix, fA$center, fA$normal, fA$radius, fA$t_center, fA$arc_half)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_arc     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_arc <- if (length(d)) mean(d)         else NA_real_
        sd_m_arc       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }
    }
  }
    
  # ---- best fit across all 6 models (heuristic + MLE), by RMSD -----------
  if (need_rmsd) {
    all_rmsd <- c(
      sphere_h    = rmsd_h_sphere,    ellipsoid_h    = rmsd_h_ellipsoid,    helix_h    = rmsd_h_helix,
      ring_h      = rmsd_h_ring,      arc_h          = rmsd_h_arc,
      sphere_mle  = rmsd_m_sphere,    ellipsoid_mle  = rmsd_m_ellipsoid,    helix_mle  = rmsd_m_helix,
      ring_mle    = rmsd_m_ring,      arc_mle        = rmsd_m_arc
    )
    finite_rmsd <- all_rmsd[is.finite(all_rmsd)]
    best_fit_name <- if (length(finite_rmsd)) names(which.min(finite_rmsd)) else NA_character_
  }

  # ---- convex hull (compute ONCE) -> center + area + volume ----
  hull <- NULL
  surface_area  <- NA_real_
  convex_volume <- NA_real_
  
  if (nrow(coords_matrix) >= 4) {
    hull <- tryCatch(
      geometry::convhulln(coords_matrix, options = "FA"),
      error = function(e) {
        log_msgs <<- c(log_msgs, paste("Median Cluster", cluster_id, "- Convex hull error:", e$message))
        NULL
      }
    )
  }
  
  if (!is.null(hull)) {
    center <- compute_convex_hull_centroid(coords_matrix)
    surface_area  <- hull$area
    convex_volume <- hull$vol
  } else {
    center <- colMeans(coords_matrix)
  }
    
    distances <- sqrt(rowSums((coords_matrix - matrix(center, nrow = nrow(coords_matrix), ncol = 3, byrow = TRUE))^2))
    radius <- median(distances)
    curvature <- if (radius > 0) 1 / radius^2 else NA
    sphere_area <- if (!is.na(radius) && radius > 0) 4 * pi * radius^2 else NA
    sphere_volume <- if (!is.na(radius) && radius > 0) (4/3) * pi * radius^3 else NA
    
    
    pca <- prcomp(coords_matrix, scale. = FALSE)
    rotated <- (coords_matrix - matrix(center, nrow = nrow(coords_matrix), ncol = 3, byrow = TRUE)) %*% pca$rotation
    d <- diff(stats::quantile(rotated[, 1], c(0.05, 0.95))) / 2
    e <- diff(stats::quantile(rotated[, 2], c(0.05, 0.95))) / 2
    f <- diff(stats::quantile(rotated[, 3], c(0.05, 0.95))) / 2
    scale_factor <- stats::quantile(sqrt(rowSums(rotated^2)), 0.90) / max(c(d, e, f))
    d <- d * scale_factor; e <- e * scale_factor; f <- f * scale_factor
    
    if (min(d, e, f) < 1e-6) {
      ellipsoid_area <- NA; ellipsoid_volume <- NA
    } else {
      # Knud Thomsen approximation (g = 1.6075); exact only for a sphere.
      # Error < ~1.2% for moderate axis ratios; increases for highly elongated shapes.
      g <- 1.6075
      ellipsoid_area <- 4 * pi * ((d^g * e^g + d^g * f^g + e^g * f^g) / 3)^(1/g)
      ellipsoid_volume <- (4/3) * pi * d * e * f
    }
    
    ellipticity <- round((max(d,e,f) - min(d,e,f))/max(d,e,f), 2)
    axes <- sort(c(d, e, f), decreasing = TRUE); a <- axes[1]; b <- axes[2]; c <- axes[3]
    # Eccentricities of the three principal cross-sections (2-D ellipse formula).
    # ecc1: major (PC1) vs intermediate (PC2) plane.
    # ecc2: major (PC1) vs minor (PC3) plane.
    # ecc3: intermediate (PC2) vs minor (PC3) plane.
    # All three are needed to fully characterise a general 3-D ellipsoid.
    ecc1 <- if (is.finite(d) && d > 1e-6) sqrt(max(0, 1 - (e/d)^2)) else NA_real_
    ecc2 <- if (is.finite(d) && d > 1e-6) sqrt(max(0, 1 - (f/d)^2)) else NA_real_
    ecc3 <- if (is.finite(e) && e > 1e-6) sqrt(max(0, 1 - (f/e)^2)) else NA_real_
    # Bending (perceived) radius — ASSUMPTION: the cluster represents a curved arc
    # where PC1 spans the chord half-length and PC2 is the sagitta (bow height).
    # Formula: R = (chord²/2 + sagitta²) / (2 * sagitta).
    # Valid for ring segments (e.g. S. aureus division ring viewed on edge) but
    # not meaningful for spheres, arbitrary ellipsoids, or helices.
    perceived_radius_1 <- ifelse(e > 0, (d^2 + e^2) / (2 * e), NA)
    perceived_radius_2 <- ifelse(f > 0, (d^2 + f^2) / (2 * f), NA)
    perceived_curv_1 <- 1/(perceived_radius_1^2)
    perceived_curv_2 <- 1/(perceived_radius_2^2)
    perceived_area_1 <- if (!is.na(perceived_radius_1) && perceived_radius_1 > 0) 4 * pi * perceived_radius_1^2 else NA
    prog_area <- surface_area/perceived_area_1
    
    coords <- coords_matrix
    tmp <- safe_nn2(coords, k = 4)
    if (is.null(tmp$nn)) {
      nn_idx  <- matrix(integer(0), nrow = 0, ncol = 0)
      nn_dist <- matrix(numeric(0),  nrow = 0, ncol = 0)
    } else {
      nn_res <- tmp$nn
      usedK  <- max(0, tmp$k_eff - 1)
      if (usedK == 0) {
        nn_idx  <- matrix(integer(0), nrow = nrow(coords_matrix), ncol = 0)
        nn_dist <- matrix(numeric(0),  nrow = nrow(coords_matrix), ncol = 0)
      } else {
        nn_idx  <- nn_res$nn.idx[, 2:(usedK + 1), drop = FALSE]
        nn_dist <- nn_res$nn.dists[, 2:(usedK + 1), drop = FALSE]
      }
    }
    
    nn_df <- data.frame()
    for (i in seq_len(nrow(coords_matrix))) {
      origin_id <- medians_df$group_id[i]
      spread_i <- medians_df$spread[i]
      for (j in seq_len(min(3, ncol(nn_idx)))) {
        nbr_idx <- nn_idx[i, j]
        if (is.na(nbr_idx) || nbr_idx < 1) next
        nbr_id <- medians_df$group_id[nbr_idx]
        spread_j <- medians_df$spread[nbr_idx]
        distance <- nn_dist[i, j]
        spread_cutoff <- max(1.5 * max(spread_i, spread_j, na.rm = TRUE), 6)
        if (is.finite(spread_cutoff) && distance > spread_cutoff) {
          nn_df <- rbind(nn_df, data.frame(
            from = origin_id, to = nbr_id, distance = distance,
            spread_i = spread_i, spread_j = spread_j, spread_cutoff = spread_cutoff
          ))
        }
      }
    }
    nn_stats <- list(
      nn_n = nrow(nn_df),
      nn_mean = if (nrow(nn_df)) mean(nn_df$distance) else NA_real_,
      nn_sd   = if (nrow(nn_df))  sd(nn_df$distance)  else NA_real_,
      nn_min  = if (nrow(nn_df)) min(nn_df$distance) else NA_real_,
      nn_max  = if (nrow(nn_df)) max(nn_df$distance) else NA_real_,
      nn_med  = if (nrow(nn_df)) median(nn_df$distance) else NA_real_
    )
    
    
    symmetry_plane_normal <- pca$rotation[,3]
    symmetry_plane_normal <- symmetry_plane_normal / sqrt(sum(symmetry_plane_normal^2))
    vectors_to_center <- t(t(coords) - center)
    distances_to_plane <- vectors_to_center %*% symmetry_plane_normal
    reflected_coords <- coords - 2 * distances_to_plane %*% t(symmetry_plane_normal)
    nn_reflect <- RANN::nn2(data = coords, query = reflected_coords, k = 1)
    mean_nn_distance <- mean(nn_reflect$nn.dists)
    bounding_radius <- max(distances)
    semi_score <- mean_nn_distance / bounding_radius
    
    if (!is.null(plane) && is.list(plane) && !is.null(plane$normal) && !is.null(plane$center)) {
      n_vec <- plane$normal / sqrt(sum(plane$normal^2))
      p0 <- plane$center
      p <- center
      vector_p0p <- p - p0
      distance_to_plane <- sum(n_vec * vector_p0p)
      foot_point <- p - distance_to_plane * n_vec
    } else {
      distance_to_plane <- NA
      foot_point <- NA
    }
    
    fit_center <- center
    if (!is.null(plane) && use_pca_normal && all(c("normal", "center") %in% names(plane))) {
      base_vector <- plane$normal / sqrt(sum(plane$normal^2))
      orange_vector <- if (distance_to_plane < 0) -base_vector else base_vector
    } else {
      orange_vector <- c(0, 0, 1)
    }
    
    angle_deg <- NA
    end_point_orange <- rep(NA_real_, 3)
    end_point_blue   <- rep(NA_real_, 3)
    if (!is.na(distance_to_plane) && abs(distance_to_plane) <= radius) {
      n_vec <- plane$normal / sqrt(sum(plane$normal^2))
      v_perp <- -distance_to_plane * n_vec
      l_parallel <- sqrt(radius^2 - distance_to_plane^2)
      z_up <- c(0, 0, 1)
      if (abs(sum(n_vec * z_up)) > 0.999) z_up <- c(0, 1, 0)
      in_plane_dir <- pracma::cross(n_vec, z_up)
      in_plane_dir <- in_plane_dir / sqrt(sum(in_plane_dir^2))
      blue_vector <- v_perp + l_parallel * in_plane_dir
      blue_vector <- blue_vector / sqrt(sum(blue_vector^2))
      if (sign(sum(blue_vector * n_vec)) == sign(distance_to_plane)) blue_vector <- -blue_vector
      end_point_orange <- center + orange_vector * radius
      end_point_blue   <- center + radius * blue_vector
      dot_product <- sum(orange_vector * blue_vector)
      angle_rad <- acos(min(max(dot_product, -1), 1))
      angle_deg <- round(angle_rad * 180 / pi, 2)
      label_offset <- 0.1 * radius
    }
    
    if (show_plane_angles) {
      fig <- fig %>% plotly::add_trace(
        x = c(center[1], end_point_orange[1]),
        y = c(center[2], end_point_orange[2]),
        z = c(center[3], end_point_orange[3]),
        type = "scatter3d", mode = "lines",
        line = list(color = "orange", width = 4),
        inherit = FALSE, showlegend = FALSE
      ) %>% plotly::add_trace(
        x = c(center[1], end_point_blue[1]),
        y = c(center[2], end_point_blue[2]),
        z = c(center[3], end_point_blue[3]),
        type = "scatter3d", mode = "lines",
        line = list(color = "blue", width = 4),
        inherit = FALSE, showlegend = FALSE
      ) %>% plotly::add_trace(
        x = c(center[1] + label_offset),
        y = c(center[2] + label_offset),
        z = c(center[3] + label_offset),
        type = "scatter3d", mode = "text",
        text = paste0(angle_deg, "°"),
        textfont = list(color = if (dark_mode) "lightblue" else "darkblue", size = 11),
        inherit = FALSE, showlegend = FALSE
      )
    }
    
    if (isTRUE(bending_radius) && !is.na(perceived_radius_1)) {
      end_point_orange <- center + orange_vector * radius
      bend_dir <- center - end_point_orange
      bend_dir <- bend_dir / sqrt(sum(bend_dir^2))
      end_point_perceived <- end_point_orange + bend_dir * perceived_radius_1
      fig <- fig %>% plotly::add_trace(
        x = c(end_point_orange[1], end_point_perceived[1]),
        y = c(end_point_orange[2], end_point_perceived[2]),
        z = c(end_point_orange[3], end_point_perceived[3]),
        type  = "scatter3d", mode  = "lines",
        line  = list(color = "black", width = 4, dash = "dash"),
        name  = "Perceived radius",
        inherit = FALSE, showlegend = TRUE
      )
    }
    
    if (show_plane_distance && !any(is.na(foot_point))) {
      fig <- fig %>% plotly::add_trace(
        x = c(center[1], foot_point[1]),
        y = c(center[2], foot_point[2]),
        z = c(center[3], foot_point[3]),
        type = "scatter3d", mode = "lines",
        line = list(color = if (dark_mode) "#F6CEFC" else "purple", width = 6, dash = 'dot'),
        inherit = FALSE, showlegend = FALSE
      ) %>% plotly::add_trace(
        x = (center[1] + foot_point[1]) / 2,
        y = (center[2] + foot_point[2]) / 2,
        z = (center[3] + foot_point[3]) / 2,
        type = "scatter3d", mode = "text",
        text = paste0(round(abs(distance_to_plane), 1)),
        textfont = list(color = if (dark_mode) "#F6CEFC" else "purple", size = 10),
        inherit = FALSE, showlegend = FALSE
      )
    }
    
    if (show_pca_axes) {
      axis_colors <- c("red", "green", "blue")
      axis_labels <- c("PC1", "PC2", "PC3")
      for (i in 1:3) {
        axis_vector <- pca$rotation[, i] * c(d, e, f)[i]
        fig <- fig %>% plotly::add_trace(
          x = c(center[1], center[1] + axis_vector[1]),
          y = c(center[2], center[2] + axis_vector[2]),
          z = c(center[3], center[3] + axis_vector[3]),
          type = "scatter3d", mode = "lines+text",
          line = list(color = axis_colors[i], width = 6, dash = "dash"),
          text = axis_labels[i], textposition = "top right",
          textfont = list(color = if (dark_mode) "white" else "black", size = 10),
          showlegend = FALSE, inherit = FALSE
        )
      }
    }
    
    if (isTRUE(show_info_labels)) {
      
      
      hover_txt <- paste0(
        "Median cluster: ", cluster_id,
        "<br>CH Surface Area=", fmt(surface_area, 1)
      )
      
      if (identical(selected_fit_type, "sphere")) {
        hover_txt <- paste0(hover_txt,
                            "<br>Radius=", fmt(radius, 1),
                            "<br>Sphere Surface Area=", fmt(sphere_area, 1))
      } else if (identical(selected_fit_type, "ellipsoid")) {
        hover_txt <- paste0(hover_txt,
                            "<br>Ellipticity=", fmt(ellipticity, 2),
                            "<br>Ecc1(PC1-PC2)=", fmt(ecc1, 2),
                            "<br>Ecc2(PC1-PC3)=", fmt(ecc2, 2),
                            "<br>Ecc3(PC2-PC3)=", fmt(ecc3, 2),
                            "<br>Ellipsoid Surface Area (approx)=", fmt(ellipsoid_area, 1))
      } else if (identical(selected_fit_type, "helix")) {
        hover_txt <- paste0(hover_txt, "<br>Helix turns=", helix_turns)
      }
      
      if (!is.na(active_rmsd)) {
        hover_txt <- paste0(
          hover_txt,
          "<br>RMSD(", selected_fit_method, "/", selected_fit_type, ")=",
          fmt(active_rmsd, 2),
          ", mean=", fmt(active_meandist, 2),
          " ± ", fmt(active_sd, 2)
        )
      }
      
      fig <- fig %>% plotly::add_trace(
        type = "scatter3d",
        mode = "markers",
        x = center[1], y = center[2], z = center[3],
        marker = list(size = 6, opacity = 0.001),  # basically invisible but hoverable
        hovertemplate = paste0(hover_txt, "<extra></extra>"),
        hoverlabel = list(font = list(size = 13)),
        showlegend = FALSE,
        inherit = FALSE,
        name = paste0("Median cluster ", cluster_id)
      )
    }

    result_table <- rbind(result_table, data.frame(
      cluster2 = cluster_id,
      center_x = center[1], center_y = center[2], center_z = center[3],
      radius = radius, curvature = curvature,
      surface_area_convex_hull = surface_area,
      surface_area_sphere = sphere_area,
      surface_area_ellipsoid = ellipsoid_area,
      ellipticity = ellipticity,
      eccentricity_1 = ecc1, eccentricity_2 = ecc2, eccentricity_3 = ecc3,
      volume_convex_hull = convex_volume,
      volume_sphere = sphere_volume,
      volume_ellipsoid = ellipsoid_volume,
      best_fit = best_fit_name,
      rmsd_sphere_h = rmsd_h_sphere, meandist_sphere_h = meandist_h_sphere, sd_sphere_h = sd_h_sphere,
      rmsd_ellipsoid_h = rmsd_h_ellipsoid, meandist_ellipsoid_h = meandist_h_ellipsoid, sd_ellipsoid_h = sd_h_ellipsoid,
      rmsd_helix_h = rmsd_h_helix, meandist_helix_h = meandist_h_helix, sd_helix_h = sd_h_helix,
      rmsd_sphere_mle = rmsd_m_sphere, meandist_sphere_mle = meandist_m_sphere, sd_sphere_mle = sd_m_sphere,
      rmsd_ellipsoid_mle = rmsd_m_ellipsoid, meandist_ellipsoid_mle = meandist_m_ellipsoid, sd_ellipsoid_mle = sd_m_ellipsoid,
      rmsd_helix_mle = rmsd_m_helix, meandist_helix_mle = meandist_m_helix, sd_helix_mle = sd_m_helix,
      rmsd_ring_h = rmsd_h_ring, meandist_ring_h = meandist_h_ring, sd_ring_h = sd_h_ring,
      rmsd_arc_h  = rmsd_h_arc,  meandist_arc_h  = meandist_h_arc,  sd_arc_h  = sd_h_arc,
      rmsd_ring_mle = rmsd_m_ring, meandist_ring_mle = meandist_m_ring, sd_ring_mle = sd_m_ring,
      rmsd_arc_mle  = rmsd_m_arc,  meandist_arc_mle  = meandist_m_arc,  sd_arc_mle  = sd_m_arc,
      pc1_length = 2*d, pc2_length = 2*e, pc3_length = 2*f,
      pc2_pc1 = (2*e)/(2*d), pc3_pc1 = (2*f)/(2*d),
      nn_n = nn_stats$nn_n, nn_mean = nn_stats$nn_mean, nn_sd = nn_stats$nn_sd,
      nn_min = nn_stats$nn_min, nn_max = nn_stats$nn_max, nn_med = nn_stats$nn_med,
      bending_radius_1 = perceived_radius_1, bending_radius_2 = perceived_radius_2,
      perceived_curv_1 = perceived_curv_1, perceived_curv_2 = perceived_curv_2,
      perceived_area_1 = perceived_area_1, prog_area = prog_area,
      distance_to_plane = distance_to_plane, closing_angle_deg = angle_deg
    ))
    
    msg <- paste0(
      "Median cluster ", cluster_id,
      ": CH Area=", fmt(surface_area, 2),
      ", CH Vol=",  fmt(convex_volume, 2)
    )
    
    shape_msg <- switch(
      as.character(selected_fit_type),
      "sphere" = paste0(
        ", Radius=", fmt(radius, 2),
        ", Sphere Area=", fmt(sphere_area, 2),
        ", Sphere Vol=", fmt(sphere_volume, 2),
        ", Curv=", ifelse(is.na(curvature), "NA", signif(curvature, 3))
      ),
      "ellipsoid" = paste0(
        ", Ellipticity=", fmt(ellipticity, 2),
        ", Ecc1(PC1-PC2)=", fmt(ecc1, 2),
        ", Ecc2(PC1-PC3)=", fmt(ecc2, 2),
        ", Ecc3(PC2-PC3)=", fmt(ecc3, 2),
        ", Ellipsoid Area(approx)=", fmt(ellipsoid_area, 2),
        ", Ellipsoid Vol=", fmt(ellipsoid_volume, 2)
      ),
      "helix" = paste0(", Helix turns=", helix_turns),
      ""
    )
    
    msg <- paste0(msg, shape_msg)
    
    if (!is.na(active_rmsd)) {
      msg <- paste0(
        msg,
        ", RMSD(", selected_fit_method, "/", selected_fit_type, ")=",
        fmt(active_rmsd, 2),
        ", mean=", fmt(active_meandist, 2),
        " ± ", fmt(active_sd, 2)
      )
    }

    log_msgs <- c(log_msgs, msg)
  }
  if (!is.null(cluster2_table_store)) cluster2_table_store(result_table)
  list(fig = fig, log = paste(log_msgs, collapse = "\n"), table = result_table)
}


compute_cluster2_plane_distances <- function(df) {
  df_medians <- df %>% group_by(tid) %>% summarize(x = median(x), y = median(y), z = median(z), .groups = "drop")
  pca <- prcomp(df_medians[, c("x", "y", "z")], scale. = FALSE)
  normal <- pca$rotation[, 3]
  center <- colMeans(df_medians[, c("x", "y", "z")])
  
  df_cluster2 <- df %>%
    filter(!is.na(cluster2) & cluster2 != 0) %>%
    group_by(cluster2) %>%
    summarize(
      x = median(x), y = median(y), z = median(z),
      trace_ids = paste(unique(tid), collapse = ","),
      .groups = "drop"
    )
  
  df_cluster2$distance_to_plane <- apply(df_cluster2[, c("x", "y", "z")], 1, function(row) {
    sum(normal * (row - center)) / sqrt(sum(normal^2))
  })
  
  df_cluster2
}
format_tc_id <- function(x) {
  ifelse(is.na(x) | x == 0, NA_character_,
         ifelse(x >= 10001L, paste0("B-", x - 10000L), paste0("A-", x)))
}

generate_xy_plane <- function(df, limits, z_offset = 0) {
  medians <- df %>%
    group_by(tid) %>%
    summarize(x = median(x), y = median(y), z = median(z), .groups = "drop")
  
  pca <- prcomp(medians[, c("x", "y", "z")], scale. = FALSE)
  center <- colMeans(medians[, c("x", "y", "z")])
  center[3] <- center[3] + z_offset  # Apply offset
  
  normal <- pca$rotation[, 3]
  grid_size <- max(diff(limits$x), diff(limits$y)) / 2
  grid_x <- seq(center[1] - grid_size, center[1] + grid_size, length.out = 20)
  grid_y <- seq(center[2] - grid_size, center[2] + grid_size, length.out = 20)
  grid <- expand.grid(x = grid_x, y = grid_y)
  
  a <- normal[1]
  b <- normal[2]
  c <- normal[3]
  d <- sum(normal * center)
  grid$z <- (d - a * grid$x - b * grid$y) / c
  
  list(
    x = matrix(grid$x, nrow = length(grid_y)),
    y = matrix(grid$y, nrow = length(grid_y)),
    z = matrix(grid$z, nrow = length(grid_y)),
    normal = normal,  
    center = center          
  )
}

#Helper: surface plot
create_global_surface <- function(df, grid_n = 120, span = 0.25) {
  grid_x <- seq(min(df$x, na.rm = TRUE), max(df$x, na.rm = TRUE), length.out = grid_n)
  grid_y <- seq(min(df$y, na.rm = TRUE), max(df$y, na.rm = TRUE), length.out = grid_n)
  grid   <- expand.grid(x = grid_x, y = grid_y)
  
  fit    <- loess(z ~ x + y, data = df, span = span, control = loess.control(surface = "direct"))
  z_pred <- predict(fit, newdata = grid)
  
  stopifnot(length(z_pred) == nrow(grid))
  list(
    x = matrix(grid$x, nrow = length(grid_y), byrow = TRUE),
    y = matrix(grid$y, nrow = length(grid_y), byrow = TRUE),
    z = matrix(z_pred,  nrow = length(grid_y), byrow = TRUE)
  )
}

# Helper: Estimate curvature and volume for cluster2 using medians
estimate_cluster2_curvature <- function(df) {
  df %>%
    filter(!is.na(cluster2) & cluster2 != 0) %>%
    group_by(cluster2) %>%
    group_modify(~{
      df_medians <- .x %>%
        group_by(tid) %>%
        summarize(x = median(x), y = median(y), z = median(z), .groups = "drop")
      coords <- as.matrix(df_medians[, c("x", "y", "z")])
      
      if (nrow(coords) >= 5) {
        pca <- prcomp(coords, scale. = FALSE)
        eigen_vals <- pca$sdev^2
        curvature_index <- eigen_vals[3] / sum(eigen_vals)
        curvature_deg <- asin(sqrt(curvature_index)) * (180 / pi)
        
        vol <- tryCatch({
          ashape <- alphashape3d::ashape3d(coords, alpha = 10)
          alphashape3d::volume_ashape3d(ashape)
        }, error = function(e) NA)
        
        # fit_sphere_geometry: hull centroid + median orthogonal radius.
        # Consistent with the main pipeline (same centre definition as
        # the heuristic sphere fit shown in the 3-D viewer).
        sphere_geo    <- fit_sphere_geometry(coords)
        radius        <- sphere_geo$radius
        sphere_curvature <- if (is.finite(radius) && radius > 0) 1 / radius else NA
        sphere_deg    <- if (is.finite(radius)) asin(min(1, sqrt(curvature_index))) * (180 / pi) else NA
        
        data.frame(
          curvature = curvature_index,
          curvature_deg = curvature_deg,
          volume = vol,
          sphere_radius = radius,
          sphere_curvature = sphere_curvature,
          sphere_curvature_deg = sphere_deg
        )
      } else {
        data.frame(
          curvature = NA, curvature_deg = NA, volume = NA,
          sphere_radius = NA, sphere_curvature = NA, sphere_curvature_deg = NA
        )
      }
    }) %>%
    ungroup()
}


# Helper: Create downloadable curvature table
build_curvature_table <- function(df) {
  curvature_df <- estimate_cluster2_curvature(df)
  
  merged <- df %>%
    filter(!is.na(cluster2) & cluster2 != 0) %>%
    group_by(cluster2) %>%
    summarize(
      trace_ids = paste(unique(tid), collapse = ", "),
      clusters = paste(unique(cluster), collapse = ", "),
      n_traces = length(unique(tid)),
      n_points = n(),
      x = median(x), y = median(y), z = median(z),
      .groups = "drop"
    )
  
  full <- merged %>% left_join(curvature_df, by = "cluster2")
  return(full)
}

## ── Secondary (cluster2) grouping engine ──────────────────────────────────────
# Runs the secondary grouping on a matrix of centroid coordinates (trace or
# primary-cluster medians) and returns an integer vector of cluster ids using
# the dbscan convention: 0 = noise/unclustered, 1..k = clusters. This keeps the
# output identical in shape regardless of method, so every downstream consumer
# (cluster2 columns, isosurfaces, fitting pipelines, NN, exports) is unaffected.
#
#   method = "dbscan"  -> density clustering with fixed radius `eps` + `minPts`
#   method = "hdbscan" -> hierarchical density clustering, `minPts` only (no eps)
#   method = "optics"  -> OPTICS reachability + Xi extraction; `minPts` + `xi`,
#                         adaptive density, keeps small structures (no big eps)
#
# Optional local-density prefilter (prefilter = TRUE): before the chosen method
# runs, any centroid with fewer than `prefilter_min` neighbours within radius
# `prefilter_r` is forced to noise. This encodes "structures have higher local
# TID density than scattered background" as an explicit, tunable gate, and works
# in front of any method.
#
# `eps` is used only by dbscan; `xi` only by optics. Works in 2D or 3D depending
# on ncol(pts).
run_cluster2_grouping <- function(pts, method = "dbscan", eps = 200, minPts = 15,
                                  xi = 0.05,
                                  prefilter = FALSE, prefilter_r = NULL,
                                  prefilter_min = NULL,
                                  prune = FALSE, prune_tol = 2) {
  if (is.null(pts)) return(integer(0))
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  n <- nrow(pts)
  if (n < 1L) return(integer(0))

  method <- tolower(as.character(method %||% "dbscan"))

  # ── optional local-density prefilter ──
  keep <- rep(TRUE, n)
  if (isTRUE(prefilter) && !is.null(prefilter_r) && isTRUE(prefilter_r > 0)) {
    pf_min <- max(1L, as.integer(round(prefilter_min %||% 1L)))
    nn <- tryCatch(dbscan::frNN(pts, eps = prefilter_r, sort = FALSE),
                   error = function(e) NULL)
    if (!is.null(nn)) keep <- lengths(nn$id) >= pf_min   # neighbours exclude self
  }

  out <- rep(0L, n)
  if (!any(keep)) return(out)
  sub <- pts[keep, , drop = FALSE]
  lab <- .cluster2_core(sub, method, eps, minPts, xi)

  # ── optional sparse-member pruning (post-clustering) ──
  if (isTRUE(prune)) lab <- .prune_sparse_members(sub, lab, tol = prune_tol)

  out[keep] <- lab
  out
}

# Post-clustering refinement: within each cluster, demote to noise (0) any
# member whose local neighbourhood (mean distance to its k nearest fellow
# members) is much sparser than the cluster's typical member. The reference is
# the per-cluster median of that distance, so it is robust on the small member
# counts typical of MINFLUX structures, and scale-free across clusters of
# different tightness. `tol` is the only knob: a member is dropped when its
# local spacing exceeds `tol` x the cluster median (tol ~ 2 = "twice as loose").
.prune_sparse_members <- function(pts, labels, tol = 2) {
  if (!length(labels)) return(labels)
  tol <- suppressWarnings(as.numeric(tol %||% 2))
  if (!is.finite(tol) || tol <= 1) return(labels)   # tol<=1 would trim core too
  out <- labels
  for (cl in unique(labels[labels > 0])) {
    idx <- which(labels == cl)
    nm  <- length(idx)
    if (nm < 4L) next                                # too few to judge sparsity
    sub <- pts[idx, , drop = FALSE]
    k   <- max(1L, min(3L, nm - 1L))
    kd  <- tryCatch(dbscan::kNN(sub, k = k)$dist, error = function(e) NULL)
    if (is.null(kd)) next
    di  <- rowMeans(kd)                              # mean dist to k nearest siblings
    ref <- stats::median(di)
    if (!is.finite(ref) || ref <= 0) next
    out[idx[di > tol * ref]] <- 0L                   # sparse members -> noise
  }
  out
}

# Single-method core; assumes points already prefiltered. Returns 0=noise/1..k.
.cluster2_core <- function(pts, method, eps, minPts, xi) {
  n <- nrow(pts)
  if (n < 1L) return(integer(0))

  if (identical(method, "hdbscan")) {
    # hdbscan needs minPts >= 2 and at least minPts points to form any cluster.
    mp <- max(2L, as.integer(round(minPts)))
    if (n < mp) return(rep(0L, n))
    res <- tryCatch(dbscan::hdbscan(pts, minPts = mp), error = function(e) NULL)
    if (is.null(res) || is.null(res$cluster)) return(rep(0L, n))
    return(as.integer(res$cluster))
  }

  if (identical(method, "optics")) {
    mp <- max(2L, as.integer(round(minPts)))
    if (n < mp) return(rep(0L, n))
    oo <- tryCatch(dbscan::optics(pts, minPts = mp), error = function(e) NULL)
    if (is.null(oo)) return(rep(0L, n))
    xv <- min(max(as.numeric(xi %||% 0.05), 1e-3), 0.999)
    cx <- tryCatch(dbscan::extractXi(oo, xi = xv), error = function(e) NULL)
    if (is.null(cx) || is.null(cx$cluster)) return(rep(0L, n))
    return(as.integer(cx$cluster))
  }

  # default: classic DBSCAN
  as.integer(dbscan::dbscan(pts, eps = eps, minPts = minPts)$cluster)
}


# =============================================================================
# deblink_axis_weights() — PUT X, Y AND Z ON A COMMON PRECISION SCALE
# =============================================================================
#
# GOAL
#   DBSCAN uses a spherical eps, but MINFLUX trace spread is not isotropic — in
#   this app z is additionally multiplied by the axial calibration (default 0.7,
#   a refractive-index mismatch correction, user-set per dataset). A single
#   radius therefore cannot fit all three axes: set it for xy and traces
#   fragment along the tight axis, set it for z and xy over-merges.
#
# WHY RESCALING IS CORRECT *FOR DE-BLINKING* SPECIFICALLY
#   "Are these two localizations the same fluorophore?" is a statistical
#   question, not a geometric one. If the axial scatter is smaller, the same
#   displacement in z is MORE evidence of two distinct molecules than in xy.
#   Weighting each axis by its own spread is the Mahalanobis form of that test.
#
#   This is emphatically NOT valid for measuring physical size — cluster radii,
#   Rg, nearest-neighbour distances must stay in true nm. The weights are for
#   the merge decision only.
#
# UNITS ARE PRESERVED
#   Weights are normalised by the geometric mean of the per-axis spreads, so the
#   transform is volume-preserving and eps keeps its nm meaning. Perfectly
#   isotropic data yields w = (1,1,1) and behaviour is identical to before,
#   which keeps saved settings comparable.
#
# ARGUMENTS
#   df    : data.frame with tid, x, y [, z]
#   min_n : minimum localizations for a trace to contribute a spread estimate
#   use_z : FALSE forces wz = 1 (2D data carries z = 0 and must not be rescaled)
#
# RETURNS
#   list(w = c(wx, wy, wz), sigma = c(sx, sy, sz), has_z, ok)
#   `w` is ALWAYS length 3 so callers can apply it to an (x,y,z) matrix
#   unconditionally. ok = FALSE means the spreads were unusable and w is
#   c(1,1,1) — the caller then behaves exactly as it did before.
# =============================================================================
deblink_axis_weights <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, use_z = TRUE) {
  flat <- list(w = c(1, 1, 1), sigma = c(NA_real_, NA_real_, NA_real_),
               has_z = FALSE, ok = FALSE)
  if (!is.data.frame(df) || !nrow(df)) return(flat)

  has_z <- isTRUE(use_z) && "z" %in% names(df) &&
    any(is.finite(df$z)) && stats::sd(df$z, na.rm = TRUE) > 1e-6

  sds <- tryCatch(sd_by_tid(df, min_n = min_n), error = function(e) NULL)
  if (!is.data.frame(sds) || !nrow(sds)) return(flat)

  med <- function(v) {
    v <- v[is.finite(v) & v > 0]
    if (length(v) < 3L) NA_real_ else stats::median(v)
  }
  s <- c(med(sds$sd_x), med(sds$sd_y), if (has_z) med(sds$sd_z) else NA_real_)

  # xy must both be measurable; z only when the data is 3D.
  if (!all(is.finite(s[1:2]))) return(flat)
  if (has_z && !is.finite(s[3])) has_z <- FALSE

  use <- if (has_z) 1:3 else 1:2
  gm  <- exp(mean(log(s[use])))                    # geometric mean -> volume-preserving
  if (!is.finite(gm) || gm <= 0) return(flat)

  w <- c(1, 1, 1)
  w[use] <- gm / s[use]
  list(w = w, sigma = s, has_z = has_z, ok = TRUE)
}

# GOAL: apply per-axis de-blink weights to an (x,y,z) coordinate matrix.
# Kept as one function so every primary-DBSCAN call site transforms identically
# — the suggestion and the clustering must share a metric or eps is meaningless.
apply_deblink_weights <- function(pts, w) {
  if (is.null(w) || all(w == 1)) return(pts)
  sweep(pts, 2L, w[seq_len(ncol(pts))], "*")
}

# =============================================================================
# deblink_merge_centroids() — MERGE TRACES ON CENTROID AGREEMENT
# =============================================================================
#
# GOAL
#   Decide "same fluorophore?" by comparing where the traces ARE, not whether
#   their point clouds touch. Each trace collapses to its median position and
#   single-linkage joins centroids closer than `thr`.
#
# WHY THIS BEATS RAW-POINT DBSCAN (measured, not assumed)
#   Raw-point single-linkage merges when ANY localization of one cloud comes
#   within eps of ANY localization of another, so the outlying tails do the
#   chaining and the effective resolution degrades as ~eps + 2*sigma*sqrt(2 ln n)
#   — it gets WORSE as traces get longer. Centroids have no tails.
#
#   On ground-truth simulation (360 molecules, 923 traces, sigma 6.5 nm):
#     raw-point DBSCAN   best at 0.62*sigma -> 35 errors (splits + merges)
#     centroid linkage   best at 1.54*sigma -> 10 errors
#   The reason is visible in the distance distributions: two re-blinks of one
#   fluorophore have centroids a median 3.1 nm apart (95th pct 7.1 nm ~ 1.1
#   sigma), while the nearest DIFFERENT molecule's centroid sits ~54 nm away.
#   Those two populations barely overlap; raw-point linkage discards that.
#
# THRESHOLD SCALE
#   ~1.5x the median per-axis trace spread. Note this is NOT the same number as
#   the raw-point eps (~0.6 sigma) — the two methods measure different distances,
#   which is why they have separate inputs rather than sharing one.
#
# NO NOISE BY CONSTRUCTION
#   minPts = 1, so every centroid is a core point and every trace lands in a
#   cluster. That matches the merge-only rule: de-blinking never removes a
#   molecule (see deblink_assign_singletons).
#
# ARGUMENTS
#   sub_df : rows to cluster (needs tid, x, y, z)
#   thr    : centroid-distance threshold in nm (weighted units if w is given)
#   w      : optional per-axis de-blink weights, same metric as everything else
#
# RETURNS  integer cluster id per ROW of sub_df, all > 0.
# =============================================================================
deblink_merge_centroids <- function(sub_df, thr, w = NULL) {
  n <- nrow(sub_df)
  if (!n) return(integer(0))
  thr <- suppressWarnings(as.numeric(thr))
  if (!length(thr) || !is.finite(thr) || thr <= 0) thr <- 1e-9

  tid <- as.character(sub_df$tid)
  P   <- apply_deblink_weights(as.matrix(sub_df[, c("x", "y", "z"), drop = FALSE]), w)

  g <- split(seq_len(n), tid)
  M <- t(vapply(g, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]),
                numeric(ncol(P))))
  # Single linkage at `thr`: DBSCAN with minPts = 1 is exactly that, and it
  # cannot emit noise.
  cl <- tryCatch(dbscan::dbscan(M, eps = thr, minPts = 1L)$cluster,
                 error = function(e) seq_len(nrow(M)))
  lab <- stats::setNames(as.integer(cl), names(g))
  unname(lab[tid])
}

# =============================================================================
# suggest_centroid_threshold() — THE MERGE DISTANCE FROM THE TRACE SPREAD
# =============================================================================
#
# GOAL
#   Set the centroid merge distance from the per-axis trace spread the data
#   actually shows, and hand back the two distance distributions needed to see
#   whether that choice is sane.
#
# THE SCALE, AND WHERE 1.5 COMES FROM
#   Two re-blinks of ONE fluorophore have centroids separated by roughly the
#   distance the label wanders between blinks. Measured on ground-truth
#   simulation, that separation is median 0.48*sigma and 95th percentile
#   1.09*sigma, while the nearest DIFFERENT molecule's centroid sits ~8*sigma
#   away. A threshold of ~1.5*sigma therefore collects essentially every
#   re-blink while still sitting far below the distinct-molecule population, and
#   the simulated error count bottoms out there (10 errors at k = 1.54, 11 at
#   the k = 2.0 default this function ships with).
#
#   THAT COEFFICIENT ASSUMES THE LABEL DOES NOT MOVE BETWEEN BLINKS, and it is
#   the weakest assumption in this whole panel. Re-running the calibration with
#   an extra between-blink displacement shows how fast the optimum runs away:
#
#     between-blink motion   0 nm   2 nm   4 nm   6 nm
#     optimal k              1.4    1.8    2.3    3.0
#
#   Split-half measurement on real files puts the centroid's own motion at ~4.4
#   nm against a 1.5 nm random floor, so a real dataset is NOT in the k = 1.4
#   column. k is therefore a visible control, defaulting to 2.0 rather than 1.5,
#   and the app cannot measure the between-blink term for you — doing so would
#   require knowing which traces are the same molecule, which is the question
#   being asked. Use the count-vs-threshold curve and the diagnostic to judge it.
#
# ANISOTROPY
#   When per-axis weights are supplied the coordinates are whitened first, so
#   every axis contributes at the same precision and sigma collapses to the
#   geometric mean of the per-axis spreads. That is why the returned `sigma` can
#   differ from a plain RMS of sigma_x/y/z.
#
# ARGUMENTS
#   df    : filtered data.frame (tid, x, y [, z])
#   min_n : minimum localizations for a trace to contribute a spread estimate
#   use_z : include z when the data is genuinely 3D
#   w     : optional per-axis de-blink weights (same metric as the clustering)
#   k     : multiplier on sigma (default 2.0, see above)
#
# RETURNS
#   list(sigma_ax, sigma, suggested, within, between, frac_merged, n_traces,
#        dims) — `within` and `between` are the raw distances the plot draws.
#   NULL when the spread cannot be estimated.
# =============================================================================
suggest_centroid_threshold <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT,
                                       use_z = TRUE, w = NULL, k = 2.0) {
  if (!is.data.frame(df) || !nrow(df) || !all(c("tid","x","y") %in% names(df))) return(NULL)
  has_z <- isTRUE(use_z) && has_axial_data(df)
  cols  <- if (has_z) c("x","y","z") else c("x","y")

  keep <- stats::complete.cases(df[, cols, drop = FALSE])
  d    <- df[keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  P   <- as.matrix(d[, cols, drop = FALSE])
  if (!is.null(w)) P <- sweep(P, 2L, w[seq_len(ncol(P))], "*")
  tid <- as.character(d$tid)

  g <- split(seq_len(nrow(P)), tid)
  g <- g[lengths(g) >= max(2L, as.integer(min_n))]
  if (length(g) < 3L) return(NULL)

  # per-axis spread in the metric that will actually be clustered in
  sd_ax <- vapply(seq_len(ncol(P)), function(k2) {
    v <- vapply(g, function(ix) stats::sd(P[ix, k2]), numeric(1))
    stats::median(v[is.finite(v)])
  }, numeric(1))
  if (!all(is.finite(sd_ax)) || any(sd_ax <= 0)) return(NULL)
  sigma <- sqrt(mean(sd_ax^2))

  # within-trace: how far a localization sits from its own trace centroid
  within <- unlist(lapply(g, function(ix) {
    c0 <- matrixStats_colMedian(P[ix, , drop = FALSE])
    sqrt(rowSums(sweep(P[ix, , drop = FALSE], 2L, c0, "-")^2))
  }), use.names = FALSE)

  # between-trace: nearest-neighbour distance between trace centroids. The near
  # population is re-blinks of one fluorophore, the far population is distinct
  # molecules — the threshold has to land between them.
  M  <- t(vapply(g, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]), numeric(ncol(P))))
  nn <- tryCatch(dbscan::kNN(M, k = 1L)$dist[, 1L], error = function(e) numeric(0))
  nn <- nn[is.finite(nn)]

  sug <- k * sigma
  list(sigma_ax   = sd_ax,
       sigma      = sigma,
       suggested  = round(sug, 1),
       k          = k,
       within     = within,
       between    = nn,
       # what fraction of traces would find a partner at this threshold
       frac_merged = if (length(nn)) mean(nn <= sug) else NA_real_,
       n_traces   = length(g),
       dims       = ncol(P))
}

# =============================================================================
# deblink_reblink_excess() — MEASURE k FROM THE DATA, NOT FROM SIMULATION
# =============================================================================
#
# GOAL
#   Find the smallest threshold at which no re-blinking is left over, using only
#   the data. This replaces guessing k from a simulation-calibrated constant.
#
# THE IDEA
#   If merging is incomplete, the leftovers are unmerged re-blinks of the same
#   fluorophore — and those sit ON TOP of each other, a few nm apart. So after
#   merging, look at the molecule nearest-neighbour distances: an excess of
#   pairs closer than ~2 sigma, far above the local density at ordinary
#   separations, means re-blinks are still being counted twice. Raise k until
#   that excess disappears.
#
#   Measured on a real file (288 traces, sigma 5.81 nm):
#     k = 1.0  15.1% of molecules within 2 sigma   78x the reference density
#     k = 1.5   6.6%                               32x
#     k = 2.0   0.0%                                0x   <- clears here
#   That is a direct, per-dataset answer where the 1.5 default was an assumption.
#
# WHAT IT CAN AND CANNOT SEE
#   It detects UNDER-merging sharply. It cannot detect OVER-merging: once the
#   excess is gone, raising k further only fuses genuine neighbours and the test
#   stays silent. So the value it returns is a LOWER BOUND — the smallest safe k
#   — and the count-vs-threshold curve shows what going higher costs.
#
#   Another route was tried and has no power here: two traces overlapping in
#   TIME cannot be the same fluorophore, which would certify wrong merges. On
#   real data, 0 of 239 trace pairs within 60 nm overlapped in time — the
#   instrument does not track two molecules at once — so that test is dead.
#
# ARGUMENTS
#   df, sigma, w : as for suggest_centroid_threshold(); sigma in the same metric
#   ks           : multipliers to try
#   ref_lo/ref_hi: the shell (nm) used as the "ordinary density" reference
#
# RETURNS  list(table, k_clear) — k_clear is NA when the excess never clears.
# =============================================================================
deblink_reblink_excess <- function(df, sigma, w = NULL, ks = seq(1, 4, by = 0.25),
                                   ref_lo = 30, ref_hi = 60) {
  if (!is.data.frame(df) || !nrow(df) || !is.finite(sigma) || sigma <= 0) return(NULL)
  cols <- intersect(c("x", "y", "z"), names(df))
  if (length(cols) < 2L) return(NULL)
  near <- 2 * sigma
  shell <- (ref_hi^length(cols) - ref_lo^length(cols))

  rows <- lapply(ks, function(k) {
    th <- k * sigma
    cl <- tryCatch(deblink_merge_centroids(df, th, w), error = function(e) NULL)
    if (is.null(cl)) return(NULL)
    M <- do.call(rbind, lapply(split(seq_len(nrow(df)), cl), function(ix)
      matrixStats_colMedian(as.matrix(df[ix, cols, drop = FALSE]))))
    if (nrow(M) < 5L) return(NULL)
    nn <- tryCatch(dbscan::kNN(M, k = 1L)$dist[, 1L], error = function(e) NULL)
    if (is.null(nn)) return(NULL)
    frac <- mean(nn < near)
    ref  <- mean(nn >= ref_lo & nn < ref_hi) / shell * near^length(cols)
    data.frame(k = k, threshold = round(th, 2), n_molecules = nrow(M),
               frac_close = round(frac, 4),
               excess = if (ref > 0) round(frac / ref, 2) else NA_real_)
  })
  tb <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(tb) || !nrow(tb)) return(NULL)

  # First k where the short-range pile-up is indistinguishable from ordinary
  # density (excess <= 1) — i.e. nothing looks like a leftover re-blink.
  ok <- which(is.finite(tb$excess) & tb$excess <= 1)
  list(table = tb, k_clear = if (length(ok)) tb$k[ok[1]] else NA_real_,
       near_nm = round(near, 2))
}

# =============================================================================
# deblink_merge_majority_vote() — LET THE LOCALIZATIONS VOTE
# =============================================================================
#
# GOAL
#   Replace "do the point clouds touch anywhere?" with "do MOST of this trace's
#   localizations think they belong to that other trace?".
#
# THE CRITERION
#   For every localization p in trace X, compare two distances:
#     d_own  |p - centroid of X, computed WITHOUT p|
#     d_alt  |p - centroid of the nearest OTHER trace Y|
#   p votes to merge X into Y when d_alt <= rho * d_own. X and Y merge when at
#   least a fraction `f` of X's localizations vote the same way.
#
# WHY LEAVE-ONE-OUT IS NOT OPTIONAL
#   The centroid of X is fitted TO X's own points, so d_own is systematically
#   too small — every point is unfairly close to its own centre. The effect goes
#   as ~1/n, so it hits the short traces hardest, which are exactly the ones
#   most in need of merging. Excluding p from its own centroid removes the bias
#   exactly, and is cheap for a mean: c^(-p) = (n*c - p) / (n - 1). That is why
#   this uses the MEAN centroid rather than the median used elsewhere.
#
# WHY THIS IS BETTER THAN THE RAW-POINT DBSCAN IT REPLACES
#   Single linkage merged when ANY point of one cloud came within eps of ANY
#   point of another, so one stray localization could fuse two molecules, and
#   the effective resolution degraded as traces got longer. Here a stray point
#   casts one vote and is outvoted.
#
#   It is also SCALE-FREE: no nanometre threshold. A diffuse trace automatically
#   demands proportionally more separation than a tight one, so heterogeneous
#   data (mixed labels, varying photon counts) needs no global sigma.
#
# SYMMETRY AND TRANSITIVITY
#   Votes are asymmetric — a short diffuse trace may vote for a long tight one
#   without the reverse. `mutual = TRUE` requires both directions, which is the
#   conservative reading. Accepted pairs then form a graph and connected
#   components become molecules; chaining is still possible but through a far
#   stricter edge test than single linkage.
#
# ARGUMENTS
#   sub_df : rows to cluster (tid, x, y, z)
#   f      : fraction of a trace's localizations that must vote to merge (0.5)
#   rho    : distance tolerance; 1.0 = the other centroid must be strictly
#            closer, >1 loosens
#   mutual : require both traces to vote for each other
#   w      : optional per-axis de-blink weights
#   k_cand : how many candidate centroids to consider per localization
#
# RETURNS  integer cluster id per ROW of sub_df, all > 0.
# =============================================================================
deblink_merge_majority_vote <- function(sub_df, f = 0.5, rho = 1.0, mutual = TRUE,
                                        w = NULL, k_cand = 4L) {
  n <- nrow(sub_df)
  if (!n) return(integer(0))
  f   <- min(max(as.numeric(f   %||% 0.5), 0.01), 1)
  rho <- max(as.numeric(rho %||% 1), 0.01)

  tid  <- as.character(sub_df$tid)
  P    <- apply_deblink_weights(as.matrix(sub_df[, c("x","y","z"), drop = FALSE]), w)
  g    <- split(seq_len(n), tid)
  tn   <- names(g); nT <- length(g)
  if (nT < 2L) return(rep(1L, n))

  nt  <- lengths(g)
  tix <- match(tid, tn)                                   # trace index per row
  C   <- t(vapply(g, function(ix) colMeans(P[ix, , drop = FALSE]), numeric(ncol(P))))

  # leave-one-out own-centroid distance
  Cl    <- (C[tix, , drop = FALSE] * nt[tix] - P) / pmax(nt[tix] - 1L, 1L)
  d_own <- sqrt(rowSums((P - Cl)^2))

  # nearest OTHER trace centroid: ask for a few candidates and skip self
  kq <- min(as.integer(k_cand) + 1L, nT)
  kn <- tryCatch(dbscan::kNN(C, k = kq, query = P), error = function(e) NULL)
  if (is.null(kn)) return(as.integer(tix))
  best_id <- integer(n); best_d <- rep(Inf, n)
  for (j in seq_len(kq)) {
    cand <- kn$id[, j]; dd <- kn$dist[, j]
    take <- !is.finite(best_d) | (best_d == Inf)
    take <- take & cand != tix
    best_id[take] <- cand[take]; best_d[take] <- dd[take]
  }

  vote <- is.finite(best_d) & best_id > 0L & (best_d <= rho * pmax(d_own, 1e-9))
  if (!any(vote)) return(as.integer(tix))

  # Decide on the TOTAL vote, then wire up every sibling that shares it.
  #
  # Each point names only its single nearest other centroid, so a fluorophore
  # that blinked three times has its votes SPLIT across the two siblings and no
  # individual candidate clears f — the trace stays unmerged even though nearly
  # every point voted to merge with something. So the decision uses the total
  # vote fraction, and edges then go to each candidate holding a real share.
  vt <- data.frame(x = tix[vote], y = best_id[vote])
  tot <- tabulate(vt$x, nbins = nT) / nt                   # share of X's points voting
  keep_x <- which(tot >= f)
  if (!length(keep_x)) return(as.integer(tix))

  key   <- paste(vt$x, vt$y, sep = "_")
  cnt   <- table(key)
  parts <- do.call(rbind, strsplit(names(cnt), "_", fixed = TRUE))
  ed <- data.frame(x = as.integer(parts[, 1]), y = as.integer(parts[, 2]),
                   nv = as.integer(cnt))
  ed <- ed[ed$x %in% keep_x, , drop = FALSE]
  # a candidate counts if it holds a real share of the votes cast, not of all
  # points — otherwise siblings cancel each other out again
  ed$share <- ed$nv / pmax(tabulate(vt$x, nbins = nT)[ed$x], 1L)
  ed <- ed[ed$share >= 0.25, , drop = FALSE]

  if (isTRUE(mutual) && nrow(ed)) {
    fwd <- paste(ed$x, ed$y, sep = "_")
    rev <- paste(ed$y, ed$x, sep = "_")
    ed  <- ed[rev %in% fwd, , drop = FALSE]
  }
  # NOTE on f: the vote fraction for two traces of the SAME molecule is ~0.5 by
  # construction (leave-one-out makes d_own and d_alt identically distributed,
  # since both centroids estimate the same position). It falls to ~0.23 at 10 nm
  # and ~0.07 at 20 nm separation. So f must sit well BELOW 0.5 — thresholding
  # at 0.5 tests against the null value and rejects half of all true merges.

  if (!nrow(ed)) return(as.integer(tix))

  # connected components by union-find
  parent <- seq_len(nT)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  for (r in seq_len(nrow(ed))) {
    ra <- find(ed$x[r]); rb <- find(ed$y[r])
    if (ra != rb) parent[rb] <- ra
  }
  comp <- vapply(seq_len(nT), find, integer(1))
  as.integer(factor(comp))[tix]
}

# =============================================================================
# deblink_minpts() — minPts IS THE TRACE-LENGTH FILTER, NOT A FREE KNOB
# =============================================================================
# GOAL
#   Derive the primary DBSCAN minPts from the Trace ID Length filter's lower
#   bound instead of exposing it as a separate control.
#
# WHY IT IS NOT A FREE PARAMETER
#   minPts set independently of the trace-length filter silently contradicts it:
#   a trace long enough to pass filtering could still be called noise purely
#   because minPts was set higher, which is a length-dependent bias applied
#   outside the filter pipeline and invisible in the filter UI.
#
#   Tying them makes the rule self-consistent: any trace that survives the
#   length filter carries enough localizations to seed a molecule on its own,
#   so the only thing that can still leave it unmerged is its own spread —
#   which is exactly what the spread filter governs.
#
# ARGUMENTS  input : the Shiny input list (reads tid_length_range)
# RETURNS    integer >= 2 (DBSCAN needs at least 2)
# =============================================================================
deblink_minpts <- function(input) {
  tlr <- input$tid_length_range
  lo  <- suppressWarnings(as.numeric(tlr[1]))
  if (!length(lo) || !is.finite(lo)) lo <- 4       # filter UI not rendered yet
  max(2L, as.integer(round(lo)))
}

# =============================================================================
# deblink_assign_singletons() — DE-BLINKING MERGES, IT NEVER REMOVES
# =============================================================================
# GOAL
#   Give every trace DBSCAN labelled noise (cluster 0) its own cluster ID, so
#   cluster 0 no longer exists after de-blinking.
#
# WHY
#   A fluorophore that blinked once, or whose localizations were too scattered
#   to link, is still a molecule. Dropping it here would delete molecules on a
#   spatial criterion applied outside the metadata and spread filters, and would
#   bias every count toward multi-blink species. Downstream code universally
#   filters `cluster != 0`, so leaving traces at 0 would remove them silently —
#   promoting them to singletons is what actually keeps them.
#
#   The cost is explicit and preferable: with a too-small radius, counts inflate
#   toward the trace count rather than quietly shrinking.
#
# ARGUMENTS  df : data.frame with `tid` and an integer `cluster` column
# RETURNS    df with cluster > 0 everywhere
# =============================================================================
deblink_assign_singletons <- function(df) {
  if (!nrow(df) || !"cluster" %in% names(df)) return(df)
  noise <- which(df$cluster == 0L | is.na(df$cluster))
  if (!length(noise)) return(df)

  nxt <- suppressWarnings(max(df$cluster, na.rm = TRUE))
  if (!is.finite(nxt) || nxt < 0) nxt <- 0L

  # one new ID per TRACE (not per point) — a noise trace stays one molecule
  df$cluster[noise] <- as.integer(nxt) +
    as.integer(factor(as.character(df$tid[noise])))
  df
}


# =============================================================================
# deblink_count_curve() — MOLECULE COUNT vs RADIUS, AND ITS KNEE   [EXPLORATORY]
# =============================================================================
#
# GOAL
#   Pick the de-blink radius from the data alone, with no model of precision and
#   no ground truth: sweep eps, run the REAL merge at each value, and plot how
#   many molecules come out.
#
# THE SHAPE OF THE CURVE AND WHY THE KNEE IS THE ANSWER
#   small eps  count is inflated — re-blinks of one fluorophore stay separate
#              and long traces fragment. The count falls STEEPLY as eps grows
#              and genuine re-blinks start to merge.
#   the knee   the point where that steep fall stops: the re-blinks have been
#              collected and nothing much is left to merge correctly.
#   large eps  a slow, steady decline — this is DISTINCT molecules being fused.
#              It never plateaus perfectly, because there is always another
#              neighbour to swallow.
#   So the knee is the largest radius that is still merging mostly re-blinks.
#
# WHY THIS NEEDS NO LABEL / TAG INFORMATION
#   Unlike a prior-based estimate (linker length, dye size), the knee is read off
#   the data's own behaviour. The label geometry is already baked into where the
#   knee falls, so asking the user to declare a tag would add an input that
#   changes nothing.
#
# WHY IT IS MARKED EXPLORATORY
#   The knee is a curvature maximum on a noisy, subsampled curve. It is stable
#   when the two regimes are well separated and unreliable when they are not —
#   which is the same crowding condition that makes de-blinking ambiguous in the
#   first place. Treat it as evidence, not as an answer, and look at the curve.
#
# ARGUMENTS
#   df       : filtered data.frame (x, y, z, tid)
#   minpts   : the minPts the real primary DBSCAN uses
#   eps_seq  : radii to try (nm)
#   w        : optional per-axis de-blink weights, so the curve is measured in
#              the same metric the clustering will use
#   max_pts  : subsample cap — the sweep runs DBSCAN once per eps
#
# RETURNS  list(curve = data.frame(eps, n_molecules), knee, n_locs, n_traces)
# =============================================================================
deblink_count_curve <- function(df, minpts = 4L, eps_seq = seq(1, 40, by = 1),
                                w = NULL, max_pts = 60000L, method = "raw") {
  if (!is.data.frame(df) || !nrow(df) || !all(c("x", "y", "tid") %in% names(df)))
    return(NULL)
  cols <- intersect(c("x", "y", "z"), names(df))
  keep <- stats::complete.cases(df[, cols, drop = FALSE])
  d    <- df[keep, , drop = FALSE]
  if (!"z" %in% names(d)) d$z <- 0
  if (nrow(d) < 10L) return(NULL)
  if (nrow(d) > max_pts) {
    set.seed(1L)                                   # reproducible curve
    d <- d[sort(sample.int(nrow(d), max_pts)), , drop = FALSE]
  }
  P   <- apply_deblink_weights(as.matrix(d[, cols, drop = FALSE]), w)
  tid <- as.character(d$tid)
  n_tr <- length(unique(tid))

  # One eps -> the molecule count the pipeline would actually report: DBSCAN,
  # then TID-plurality, then unmerged traces become their own molecule.
  # Both criteria put every trace in a cluster, so the count is just the number
  # of distinct labels. For "vote" the swept quantity is the vote FRACTION, not
  # a distance — the caller supplies the right sequence.
  Mc <- NULL
  if (identical(method, "centroid")) {
    gg <- split(seq_len(nrow(P)), tid)
    Mc <- t(vapply(gg, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]),
                   numeric(ncol(P))))
  }
  # CENTROID ONLY. A vote sweep was tried and withdrawn: it did not return in a
  # usable time, and an unresponsive diagnostic is worse than an absent one.
  # The vote criterion has no distance to sweep anyway — its parameter is a
  # dimensionless fraction whose useful range is already documented on the input.
  if (is.null(Mc)) return(NULL)
  count_at <- function(e) {
    cl <- tryCatch(dbscan::dbscan(Mc, eps = e, minPts = 1L)$cluster,
                   error = function(err) NULL)
    if (is.null(cl)) NA_integer_ else length(unique(cl))
  }
  cnt <- vapply(eps_seq, count_at, numeric(1))
  ok  <- is.finite(cnt)
  if (sum(ok) < 5L) return(NULL)
  es <- eps_seq[ok]; cn <- cnt[ok]

  # Knee = maximum curvature of log(count) against eps. Logs because the fall is
  # multiplicative; curvature because the knee is where the SLOPE changes, and
  # the slope itself never reaches zero.
  knee <- NA_real_
  if (length(es) >= 5L) {
    lc <- log(pmax(cn, 1))
    d1 <- diff(lc) / diff(es)
    d2 <- diff(d1) / diff(es[-1])
    if (any(is.finite(d2))) knee <- es[which.max(replace(d2, !is.finite(d2), -Inf)) + 1L]
  }
  list(curve = data.frame(eps = es, n_molecules = cn),
       knee = knee, n_locs = nrow(d), n_traces = n_tr)
}

# =============================================================================
# apply_clustering_pipeline() — THE TWO-TIER CLUSTERING ENTRY POINT
# =============================================================================
#
# GOAL
#   Turn raw localizations into the two grouping levels the whole SMLM tab is
#   built on. This mirrors exactly what the interactive plot does
#   (`plottable_data()` in server_smlm.R), so batch tools and the on-screen plot
#   can never disagree — that duplication is the reason this function exists.
#
# THE TWO TIERS AND WHY THERE ARE TWO
#   PRIMARY clustering (`cluster`) — DE-BLINKING.
#     Runs DBSCAN over RAW localizations with a SMALL eps (~10-30 nm, i.e. a few
#     times the localization precision). A single fluorophore blinking on and off
#     produces many localizations within a few nm of each other; this collapses
#     them into one object. Without it, one molecule is counted many times and
#     every downstream count and spatial statistic is inflated.
#
#   SECONDARY / "median" clustering (`cluster2`) — STRUCTURE.
#     Runs over the MEDIAN position of each unit from the level above, with a
#     LARGE eps (~100-300 nm). Its job is biological: group molecules into the
#     structure they form (a filament, a ring, a complex).
#
#   The two eps values are doing completely different jobs — primary eps is set
#   by the MICROSCOPE (precision + label size), secondary eps by the BIOLOGY
#   (how far apart molecules in one structure sit). Tuning them together is the
#   most common source of confusion.
#
# WHY DBSCAN (and its alternatives here)
#   DBSCAN needs no pre-set cluster count, finds arbitrary shapes, and labels
#   sparse points as noise (cluster 0) rather than forcing them into a group —
#   all three matter for SMLM data. It does assume roughly uniform density, so
#   the secondary level also offers:
#     HDBSCAN ... adaptive density; only a minimum cluster size, no eps
#     OPTICS .... adaptive; keeps small structures a single eps would merge
#
# MULTIPLEX RULE (enforced, not optional)
#   When a multiplex dataset is active, clustering runs SEPARATELY PER CHANNEL
#   and the IDs are offset. Molecules imaged in different exchange rounds must
#   never be merged into one cluster — that would invent structures that were
#   never simultaneously present.
#
# DE-BLINKING MERGES, IT NEVER REMOVES
#   The primary tier has no "drop unclustered" behaviour and no free minPts:
#     - minPts comes from the Trace ID Length filter (see deblink_minpts), so it
#       cannot contradict the filter the user actually set.
#     - traces DBSCAN cannot merge become their own single-blink molecule (see
#       deblink_assign_singletons), so `cluster` is > 0 everywhere.
#   Removing molecules is the job of the metadata and spread filters, where it
#   is visible and reversible — not a side effect of a clustering radius.
#
# ARGUMENTS
#   df    : an ALREADY-FILTERED data.frame (x, y, z, tid [, mplex_channel])
#   input : the Shiny input list — read directly so batch tools use exactly the
#           settings the user sees. Relevant: enable_clustering, eps,
#           deblink_anisotropic, tid_length_range (-> minPts), enable_cluster2,
#           eps2/minpts2, cluster2_method, grouping_mode, cluster_group_scope
#
# RETURNS
#   the same df with `cluster` and/or `cluster2` integer columns appended.
#   For `cluster2`, ID 0 still means "noise / unassigned" and consumers must
#   exclude it before counting. For `cluster`, 0 no longer occurs — the existing
#   `cluster != 0` filters downstream are harmless no-ops that now keep every
#   trace which survived filtering.
# =============================================================================
apply_clustering_pipeline <- function(df, input) {
  if (!nrow(df)) return(df)

  df <- assign_analysis_group(df, input)
  if (!nrow(df)) return(df)

  scope <- input$cluster_group_scope %||% "all"
  mode  <- input$grouping_mode       %||% "none"

  # ── Primary clustering ────────────────────────────────────────────────────
  if (isTRUE(input$enable_clustering)) {
    df$cluster <- 0L

    # Per-axis metric, computed ONCE over the whole df so every group/channel
    # subset below is clustered in the same space (a per-subset metric would
    # make cluster IDs incomparable across groups).
    aw1 <- if (isTRUE(input$deblink_anisotropic))
      deblink_axis_weights(df, min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT)$w
      else NULL

    # Two merge criteria, both trace-level. "centroid" compares where the traces
    # ARE; "vote" asks whether most of a trace's localizations think they belong
    # to a neighbouring trace. Raw-point DBSCAN was removed: one stray
    # localization could bridge two molecules. Ground-truth errors at each
    # method's own best threshold: raw-point 35, vote 29, centroid 10.
    # (Centroid scores 11 at the k = 2.0 default and 10 at its optimum 1.54 —
    #  quote the basis, not a bare number.)
    method1 <- input$deblink_method %||% "centroid"
    thr1    <- input$eps_centroid   %||% 10

    run_dbscan1 <- function(sub_df) {
      if (identical(method1, "vote"))
        return(deblink_merge_majority_vote(sub_df,
                 f = input$vote_frac %||% 0.2, rho = input$vote_rho %||% 1.2,
                 mutual = !identical(input$vote_link %||% "mutual", "either"), w = aw1))
      deblink_merge_centroids(sub_df, thr1, aw1)
    }

    mplex_on <- "mplex_channel" %in% names(df) && any(!is.na(df$mplex_channel))
    if (mplex_on) {
      # Multiplex: confine clusters to a single channel.
      next_id <- 1L
      for (chv in sort(unique(df$mplex_channel[!is.na(df$mplex_channel)]))) {
        idx <- which(df$mplex_channel == chv)
        if (!length(idx)) next
        cid <- run_dbscan1(df[idx, ])
        cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
        df$cluster[idx] <- cid
        if (any(cid > 0)) next_id <- max(next_id, max(cid) + 1L)
      }
    } else if (scope == "all" || mode == "none") {
      df$cluster <- run_dbscan1(df)
    } else {
      groups_to_run <- switch(scope, g1 = 1L, g2 = 2L, split = c(1L, 2L), integer(0))
      next_id <- 1L
      for (g in groups_to_run) {
        idx <- which(df$analysis_group == g)
        if (!length(idx)) next
        cid <- run_dbscan1(df[idx, ])
        cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
        df$cluster[idx] <- cid
        if (any(cid > 0)) next_id <- max(next_id, max(cid) + 1L)
      }
    }

    # TID-plurality reassignment
    tid_plurality <- df %>%
      dplyr::filter(cluster != 0L) %>%
      dplyr::group_by(tid, cluster) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(tid) %>%
      dplyr::slice_max(n, n = 1L, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(tid, cluster_new = cluster)
    df <- df %>%
      dplyr::left_join(tid_plurality, by = "tid") %>%
      dplyr::mutate(cluster = dplyr::coalesce(cluster_new, 0L)) %>%
      dplyr::select(-cluster_new)

    # Merge-only: unmerged traces become their own molecule rather than being
    # dropped. See deblink_assign_singletons().
    df <- deblink_assign_singletons(df)
  } else {
    df$cluster <- NA_integer_
  }

  # ── Secondary clustering ──────────────────────────────────────────────────
  if (isTRUE(input$enable_cluster2)) {
    eps2    <- input$eps2    %||% 200
    minpts2 <- input$minpts2 %||% 3
    use_primary <- isTRUE(input$enable_clustering) && "cluster" %in% names(df)

    if (use_primary) {
      medians <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0L) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          id = dplyr::first(cluster),
          mx = stats::median(x, na.rm = TRUE),
          my = stats::median(y, na.rm = TRUE),
          mz = stats::median(z, na.rm = TRUE),
          analysis_group = {
            v <- analysis_group[!is.na(analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups = "drop"
        )
    } else {
      medians <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          id = dplyr::first(tid),
          mx = stats::median(x, na.rm = TRUE),
          my = stats::median(y, na.rm = TRUE),
          mz = stats::median(z, na.rm = TRUE),
          analysis_group = dplyr::first(analysis_group),
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups = "drop"
        )
    }

    medians$cluster2 <- 0L
    run_dbscan2 <- function(idx, next_id) {
      if (!length(idx)) return(list(cid = integer(0), next_id = next_id))
      pts <- as.matrix(medians[idx, c("mx", "my", "mz")])
      cid <- run_cluster2_grouping(
        pts,
        method        = input$cluster2_method %||% "dbscan",
        eps           = eps2,
        minPts        = minpts2,
        xi            = input$cluster2_xi %||% 0.05,
        prefilter     = isTRUE(input$cluster2_prefilter),
        prefilter_r   = input$cluster2_prefilter_r,
        prefilter_min = input$cluster2_prefilter_min,
        prune         = isTRUE(input$cluster2_prune),
        prune_tol     = input$cluster2_prune_tol %||% 2
      )
      if (any(cid > 0)) {
        cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
        next_id <- max(next_id, max(cid) + 1L)
      }
      list(cid = cid, next_id = next_id)
    }

    mplex_on2 <- "mplex_channel" %in% names(medians) && any(!is.na(medians$mplex_channel))
    if (mplex_on2) {
      next_id <- 1L
      for (chv in sort(unique(medians$mplex_channel[!is.na(medians$mplex_channel)]))) {
        idx <- which(medians$mplex_channel == chv)
        res <- run_dbscan2(idx, next_id)
        medians$cluster2[idx] <- res$cid
        next_id <- res$next_id
      }
    } else if (scope == "all" || mode == "none") {
      res <- run_dbscan2(seq_len(nrow(medians)), 1L)
      medians$cluster2 <- res$cid
    } else {
      groups_to_run <- switch(scope, g1 = 1L, g2 = 2L, split = c(1L, 2L), integer(0))
      next_id <- 1L
      for (g in groups_to_run) {
        idx <- which(medians$analysis_group == g)
        res <- run_dbscan2(idx, next_id)
        medians$cluster2[idx] <- res$cid
        next_id <- res$next_id
      }
    }

    join_key <- if (use_primary) c("cluster" = "id") else c("tid" = "id")
    df <- df %>% dplyr::left_join(medians[, c("id", "cluster2")], by = join_key)
  } else {
    df$cluster2 <- NA_integer_
  }

  df
}

## ── Inverse GT-position estimator (labeling-aware) ────────────────────────────
# Reverse-engineers the Simulator's forward labeling model to turn observed
# traces/clusters back into CANDIDATE ground-truth molecule positions, an
# approximate molecule count, and a crowding validity flag.
#
# Scientific model (inverts helpers_shared::minflux_label_presets() +
# sim_server's expand_to_fluorophores / split_for_blinking / gen_localizations):
#
#   one GT molecule
#     --(link_nm)-->    ONE common-mode complex offset, shared by all its dyes
#     --(n_fluor)-->    Poisson(n_fluor) fluorophores, each scattered dye_spread_nm
#     --(blink)-->      each fluorophore -> (1 + Poisson(blink_lambda)) TIDs
#     --(sigma_loc)-->  each TID -> many localizations, precision sigma_loc
#
# So the labeling multiplicity at the RAW-TID level is
#     m_tid = n_fluor * (eff/100) * (1 + blink_lambda)      (TIDs per molecule)
# and the PRIMARY-DBSCAN step (eps ~ dye_spread) collapses re-blinks and the
# dyes of one molecule back into ~one cluster => ~one unit per molecule, turning
# an unstable /m_tid calibration into a /~1 count. Primary clusters are therefore
# the preferred unit when present ("auto").
#
# KEY LIMIT: the link_nm offset is COMMON-MODE within a molecule (all its
# localizations share the one draw), so it never averages out — it is the
# irreducible position-uncertainty floor. Reported as `pos_sd_nm`.
#
# Args:
#   df       clustered SMLM frame (needs x,y[,z],tid; uses `cluster`/`cluster2`)
#   preset   a preset list from minflux_label_presets(), OR its name, OR a list
#            with n_fluor/link_nm/dye_spread_nm/footprint_nm/blink_lambda
#   eff_pct  labeling/detection efficiency (%) — folds into the multiplicity
#   level    "auto" (primary clusters when available else TIDs), "primary", "tid"
#   fuse_k   a unit whose per-axis spread exceeds fuse_k x the single-molecule
#            expectation is flagged as fused (multiple GTs merged)
#   minpts   localizations below which a unit is ignored (matches primary DBSCAN)
#   sigma_loc_nm  per-localization precision; estimated from data when NULL
#   scope_cluster2 restrict to points inside a secondary structure (cluster2>0)
#
# Returns list(points, summary, meta) — see field construction below. Empty
# data.frames (never NULL) when nothing is estimable, so callers stay simple.
estimate_gt_from_labeling <- function(df, preset, eff_pct = 70,
                                      level = c("auto", "primary", "tid"),
                                      fuse_k = 2.5, minpts = 3L,
                                      sigma_loc_nm = NULL,
                                      scope_cluster2 = TRUE,
                                      complex_r = NULL,
                                      displace_inward = TRUE,
                                      displace_nm = NULL,
                                      guided_merge = FALSE,
                                      support_strength = 1) {
  level <- match.arg(level)
  empty_pts <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0),
                          pos_sd_nm = numeric(0), n_loc = integer(0),
                          n_tid = integer(0), n_units = integer(0),
                          fused = logical(0), cluster2 = integer(0),
                          mol = integer(0), level = character(0))
  empty_links <- data.frame(mol = integer(0), ux = numeric(0), uy = numeric(0),
                            uz = numeric(0), gx = numeric(0), gy = numeric(0),
                            gz = numeric(0))
  empty <- list(points = empty_pts, links = empty_links,
                summary = data.frame(), meta = list())
  if (is.null(df) || !nrow(df) || !all(c("x", "y", "tid") %in% names(df)))
    return(empty)

  # ── resolve preset parameters ──
  if (is.character(preset)) preset <- minflux_label_presets()[[preset]]
  gp <- function(k, d) { v <- suppressWarnings(as.numeric(preset[[k]])); if (length(v) && is.finite(v)) v else d }
  n_fluor      <- max(1,  gp("n_fluor", 1))
  link_nm      <- max(0,  gp("link_nm", 5))
  dye_spread   <- max(0,  gp("dye_spread_nm", 0))
  footprint    <- max(0,  gp("footprint_nm", link_nm * 1.5))
  blink_lambda <- max(0,  gp("blink_lambda", 0.3))
  eff          <- min(1, max(0.01, eff_pct / 100))
  preset_note  <- if (!is.null(preset$note)) as.character(preset$note) else ""

  # TIDs expected per molecule, and estimator detection efficiency eta =
  # P(a molecule yields >=1 detected unit). For single-dye labels eta = eff;
  # for multivalent labels the Poisson(n_fluor) count makes >=1 detection near-certain.
  m_tid <- n_fluor * eff * (1 + blink_lambda)
  eta   <- if (n_fluor <= 1) eff else 1 - exp(-n_fluor * eff)
  eta   <- min(1, max(1e-3, eta))

  use_z <- "z" %in% names(df) && any(is.finite(df$z) & df$z != 0)

  # ── estimate per-localization precision from data when not supplied ──
  # Each TID's within-trace coordinate SD ~ sigma_loc; take the median over TIDs.
  if (is.null(sigma_loc_nm) || !is.finite(sigma_loc_nm)) {
    per_tid_sd <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(sx = stats::sd(x, na.rm = TRUE),
                       sy = stats::sd(y, na.rm = TRUE),
                       n  = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n >= 3L, is.finite(sx), is.finite(sy))
    sigma_loc_nm <- if (nrow(per_tid_sd))
      stats::median(sqrt((per_tid_sd$sx^2 + per_tid_sd$sy^2) / 2), na.rm = TRUE) else 2
    if (!is.finite(sigma_loc_nm) || sigma_loc_nm <= 0) sigma_loc_nm <- 2
  }
  # expected per-axis spread of a SINGLE molecule's localizations (dye scatter + noise)
  r_single <- sqrt(dye_spread^2 + sigma_loc_nm^2)

  # ── choose the unit level ──
  has_primary <- "cluster" %in% names(df) && any(!is.na(df$cluster) & df$cluster > 0L)
  unit_level  <- switch(level,
                        auto    = if (has_primary) "primary" else "tid",
                        primary = "primary",
                        tid     = "tid")
  if (unit_level == "primary" && !has_primary) {
    if (level == "primary") return(empty)      # explicitly asked for primary, none exist
    unit_level <- "tid"
  }

  # ── optionally scope to points inside a secondary structure ──
  has_c2 <- "cluster2" %in% names(df) && any(!is.na(df$cluster2) & df$cluster2 > 0L)
  work <- df
  if (isTRUE(scope_cluster2) && has_c2)
    work <- work[!is.na(work$cluster2) & work$cluster2 > 0L, , drop = FALSE]
  if (!nrow(work)) return(empty)
  if (!"cluster2" %in% names(work)) work$cluster2 <- 0L
  work$cluster2[is.na(work$cluster2)] <- 0L
  if (unit_level == "primary")
    work <- work[!is.na(work$cluster) & work$cluster > 0L, , drop = FALSE]
  if (!nrow(work)) return(empty)

  unit_col <- if (unit_level == "primary") "cluster" else "tid"

  # ── build candidate units ──
  if (!"z" %in% names(work)) work$z <- 0    # ensure column exists; zeroed later if 2D
  work$z[!is.finite(work$z)] <- 0
  # NOTE: compute the SDs BEFORE the median columns overwrite x/y/z — dplyr
  # summarise evaluates sequentially, so `sd_x` must reference the raw x, not the
  # just-created scalar median (which would give NA).
  units <- work %>%
    dplyr::group_by(.uid = .data[[unit_col]]) %>%
    dplyr::summarise(
      sd_x     = stats::sd(x, na.rm = TRUE),
      sd_y     = stats::sd(y, na.rm = TRUE),
      sd_z     = if (use_z) stats::sd(z, na.rm = TRUE) else 0,
      n_loc    = dplyr::n(),
      n_tid    = dplyr::n_distinct(tid),
      cluster2 = { v <- cluster2[cluster2 > 0L]; if (length(v)) as.integer(names(which.max(table(v)))) else 0L },
      x        = stats::median(x, na.rm = TRUE),
      y        = stats::median(y, na.rm = TRUE),
      z        = if (use_z) stats::median(z, na.rm = TRUE) else 0,
      .groups  = "drop"
    ) %>%
    dplyr::filter(n_loc >= minpts, is.finite(x), is.finite(y))
  if (!nrow(units)) return(empty)

  units$sd_x[!is.finite(units$sd_x)] <- 0
  units$sd_y[!is.finite(units$sd_y)] <- 0
  units$sd_z[!is.finite(units$sd_z)] <- 0

  # observed per-axis spread of each unit (RMS across axes)
  ndim        <- if (use_z) 3 else 2
  units$unit_spread <- sqrt((units$sd_x^2 + units$sd_y^2 +
                             if (use_z) units$sd_z^2 else 0) / ndim)

  # ── MERGE units into candidate molecules ──────────────────────────────────
  # A multivalent label (n_fluor>1, dye_spread>0) scatters several fluorophores,
  # hence several primary clusters / TIDs, around ONE molecule. Pool unit
  # centroids that fall within `complex_r` (the label-complex radius) so those
  # contributors collapse into a single candidate GT — this is what makes
  # "many traces -> one molecule" happen, and it produces the contribution
  # fan-out lines. Merging is confined within each cluster2 so distinct
  # structures never pool together.
  #
  # AUTO radius is PRESET-SPECIFIC: max(footprint, 2.5 x dye_spread). Both terms
  # come from the labeling preset, so a big polyclonal+secondary complex pools
  # over ~45 nm while a nanobody pools over ~3 nm — i.e. the labeling strategy
  # itself sets how many traces collapse into one molecule, as expected.
  if (is.null(complex_r) || !is.finite(complex_r) || complex_r < 0)
    complex_r <- max(footprint, 2.5 * dye_spread)
  coord_cols <- if (use_z) c("x", "y", "z") else c("x", "y")

  # `m_tid` = expected TIDs per molecule from the labeling strategy; the
  # stoichiometric molecule count of a structure is (its TID count) / m_tid.
  tids_per_mol <- max(1e-6, m_tid)

  units$mol <- NA_integer_
  next_mol  <- 1L
  for (sid in unique(units$cluster2)) {
    ix <- which(units$cluster2 == sid)
    if (length(ix) == 1L || complex_r <= 0) {
      units$mol[ix] <- next_mol + seq_along(ix) - 1L
      next_mol <- next_mol + length(ix); next
    }
    cc <- as.matrix(units[ix, coord_cols, drop = FALSE])
    g  <- tryCatch(dbscan::dbscan(cc, eps = complex_r, minPts = 1L)$cluster,
                   error = function(e) seq_along(ix))     # minPts=1 -> no noise
    g  <- as.integer(g)

    # ── stoichiometry-guided merge (opt-in, MODEL-DRIVEN) ──
    # If the spatial merge left MORE groups than the labeling stoichiometry
    # predicts (N_TID / m_tid), collapse further to that target via hierarchical
    # clustering. Only ever REDUCES the count toward stoichiometry — it never
    # splits (that would fabricate molecules). Confined to real structures.
    if (isTRUE(guided_merge) && sid > 0L) {
      n_grp    <- length(unique(g))
      k_stoich <- max(1L, round(sum(units$n_tid[ix]) / tids_per_mol))
      if (k_stoich < n_grp) {
        hc <- tryCatch(stats::hclust(stats::dist(cc), method = "average"),
                       error = function(e) NULL)
        if (!is.null(hc)) g <- as.integer(stats::cutree(hc, k = k_stoich))
      }
    }

    g <- g + next_mol - 1L
    units$mol[ix] <- g
    next_mol <- max(g) + 1L
  }

  # ── aggregate contributing units into molecules ───────────────────────────
  # Member-centroid SDs are computed BEFORE the weighted-mean position columns
  # (sequential-summarise trap), then the n_loc-weighted molecule position.
  # Order matters (dplyr summarise is sequential): SD/spread and the n_loc-
  # weighted positions must be computed while x/y/z and n_loc still hold the raw
  # per-unit vectors; the summed n_tid/n_loc come last so they don't shadow them.
  mols <- units %>%
    dplyr::group_by(mol) %>%
    dplyr::summarise(
      cluster2 = cluster2[1],
      n_units  = dplyr::n(),
      msd_x    = if (dplyr::n() > 1L) stats::sd(x) else 0,
      msd_y    = if (dplyr::n() > 1L) stats::sd(y) else 0,
      msd_z    = if (use_z && dplyr::n() > 1L) stats::sd(z) else 0,
      inner_sp = stats::median(unit_spread, na.rm = TRUE),
      x        = stats::weighted.mean(x, pmax(1, n_loc)),
      y        = stats::weighted.mean(y, pmax(1, n_loc)),
      z        = if (use_z) stats::weighted.mean(z, pmax(1, n_loc)) else 0,
      n_tid    = sum(n_tid),
      n_loc    = sum(n_loc),
      .groups  = "drop"
    )
  mols$msd_x[!is.finite(mols$msd_x)] <- 0
  mols$msd_y[!is.finite(mols$msd_y)] <- 0
  mols$msd_z[!is.finite(mols$msd_z)] <- 0
  mols$inner_sp[!is.finite(mols$inner_sp)] <- 0

  # spread of the fluorophore cloud that formed this molecule (per-axis RMS)
  member_spread <- sqrt((mols$msd_x^2 + mols$msd_y^2 +
                         if (use_z) mols$msd_z^2 else 0) / ndim)

  # Total spatial extent of a molecule = its members' fan-out (member_spread)
  # combined with the internal spread of those members (inner_sp). This catches
  # BOTH failure modes: several units spread wide (crowding merge) AND a single
  # over-wide primary cluster that already chained many molecules together.
  mol_extent <- sqrt(member_spread^2 + mols$inner_sp^2)

  # ── fusion: a molecule far wider than the expected single-molecule label
  # cloud is really several merged molecules (crowding / under-resolution).
  # Expected cloud = the dye-spread scatter plus localization noise; a tightly
  # labeled single molecule stays well under fuse_k x that.
  expected_cloud <- max(dye_spread, r_single, 1e-6)
  mols$fused <- mol_extent > (fuse_k * expected_cloud)

  # ── displace GT inward toward the local structure (linker correction) ─────
  # In the forward model the antibody binds a surface epitope and the linker
  # throws the fluorophore OUTWARD, so a single molecule's trace-cloud centroid
  # sits ~link_nm OUTSIDE the true structure. Push each (non-fused) molecule back
  # toward its structure by `displace_nm` (auto = link_nm). Geometry-aware: for
  # an elongated structure (a filament) "inward" is PERPENDICULAR toward the axis
  # (the axial component is dropped — there is no inside along the length); for a
  # compact / ring structure it points straight at the centroid. FUSED molecules
  # are skipped: their centroid geometry is already corrupted by merging, so the
  # clean outward-offset assumption doesn't hold. Applied only within a defined
  # structure (cluster2 > 0), and it never overshoots the centre.
  disp <- if (is.null(displace_nm) || !is.finite(displace_nm)) link_nm else max(0, displace_nm)
  meta_disp <- 0
  if (isTRUE(displace_inward) && disp > 0 && nrow(mols) >= 2L) {
    mols$disp_step <- 0                          # how far each GT was pushed (nm)
    mols$dir_conf  <- 0.4                         # confidence in the inward direction
    for (sid in unique(mols$cluster2)) {
      if (sid <= 0L) next                       # need a defined structure to aim at
      ix <- which(mols$cluster2 == sid)
      if (length(ix) < 3L) next
      M <- as.matrix(mols[ix, coord_cols, drop = FALSE])
      w <- pmax(1, mols$n_loc[ix])
      C <- colSums(M * w) / sum(w)              # n_loc-weighted structure centre
      axis1 <- NULL; elongated <- FALSE
      cv <- tryCatch(stats::cov(M), error = function(e) NULL)
      if (!is.null(cv)) {
        ev   <- eigen(cv, symmetric = TRUE)
        lam  <- pmax(ev$values, 0)
        rest <- if (length(lam) > 1) mean(lam[-1]) else 0
        if (lam[1] > 0 && rest > 0 && lam[1] / rest > 3) {
          elongated <- TRUE; axis1 <- ev$vectors[, 1]
        }
      }
      # inward direction is well-defined for an elongated structure (perp to a
      # clear axis) but ambiguous for a compact blob -> lower directional conf.
      mols$dir_conf[ix] <- if (elongated) 0.8 else 0.4
      for (k in ix) {
        if (isTRUE(mols$fused[k])) next          # don't displace ambiguous fused GTs
        p  <- as.numeric(mols[k, coord_cols])
        d  <- C - p
        if (elongated && !is.null(axis1)) d <- d - sum(d * axis1) * axis1  # perp to axis
        nd <- sqrt(sum(d^2))
        if (nd < 1e-6) next
        step <- min(disp, nd)                   # cap: never cross the centre
        newp <- p + step * d / nd
        mols$x[k] <- newp[1]; mols$y[k] <- newp[2]
        if (use_z) mols$z[k] <- newp[3]
        mols$disp_step[k] <- step
      }
    }
    meta_disp <- disp
  }
  if (!"disp_step" %in% names(mols)) { mols$disp_step <- 0; mols$dir_conf <- 1 }

  # ── position uncertainty (this drives the green->red colour) ──────────────
  # Independent contributions added in quadrature (all in nm):
  #  s_link    linker floor — irreducible, common-mode per molecule (dominant)
  #  s_sample  how well the label centroid is pinned = fan-out spread averaged
  #            over the contributing units, plus localization precision / #locs
  #  s_support RELIABILITY penalty for few independent labels: a GT from 1 TID
  #            is far less trustworthy than one from many (~link_nm/sqrt(n_tid))
  #  s_disp    the inward displacement guessed a direction; its error scales with
  #            how far we pushed and how ambiguous the direction was (1-dir_conf)
  #  s_fused   full extent when several molecules were merged into one GT
  s_link    <- link_nm
  s_sample  <- sqrt(member_spread^2 / pmax(1, mols$n_units) +
                    sigma_loc_nm^2   / pmax(1, mols$n_loc))
  s_support <- max(0, support_strength) * link_nm / sqrt(pmax(1, mols$n_tid))
  s_disp    <- mols$disp_step * (1 - mols$dir_conf)
  s_fused   <- ifelse(mols$fused, mol_extent, 0)
  mols$pos_sd_nm <- sqrt(s_link^2 + s_sample^2 + s_support^2 + s_disp^2 + s_fused^2)

  points <- data.frame(
    x = mols$x, y = mols$y, z = if (use_z) mols$z else 0,
    pos_sd_nm = mols$pos_sd_nm, n_loc = mols$n_loc, n_tid = mols$n_tid,
    n_units = mols$n_units, fused = mols$fused,
    cluster2 = mols$cluster2, mol = mols$mol, level = unit_level,
    stringsAsFactors = FALSE
  )

  # ── contribution links: one row per member unit -> its molecule centroid ──
  links <- data.frame(mol = units$mol, ux = units$x, uy = units$y,
                      uz = if (use_z) units$z else 0, stringsAsFactors = FALSE)
  links <- merge(links, mols[, c("mol", "x", "y", "z")], by = "mol", all.x = TRUE)
  names(links)[match(c("x", "y", "z"), names(links))] <- c("gx", "gy", "gz")
  # keep only real fan-outs (molecules with >=2 contributing units)
  multi <- mols$mol[mols$n_units >= 2L]
  links <- links[links$mol %in% multi, , drop = FALSE]

  # ── per-molecule count contribution (fused = several TIDs-per-molecule) ────
  mol_est <- ifelse(mols$fused, pmax(1, round(mols$n_tid / tids_per_mol)), 1)

  # ── per-structure summary (grouped by cluster2; 0 = ungrouped/global) ─────
  # Two INDEPENDENT molecule-count estimators per structure:
  #   spatial      = merged GT groups, efficiency-corrected (/eta)
  #   stoichiometric = total TIDs / m_tid  (already efficiency-aware via eff)
  # When they disagree by >2x the merge radius / labeling assumption is off ->
  # flagged. The reported count keeps the spatial value but the CI spans both.
  struct_ids <- sort(unique(mols$cluster2))
  summ <- lapply(struct_ids, function(sid) {
    idx <- which(mols$cluster2 == sid)
    nu  <- length(idx)
    mean_nn <- NA_real_
    if (nu >= 2L) {
      cc <- as.matrix(mols[idx, coord_cols, drop = FALSE])
      nn <- tryCatch(RANN::nn2(cc, cc, k = 2L)$nn.dists[, 2], error = function(e) NULL)
      if (!is.null(nn)) mean_nn <- stats::median(nn[is.finite(nn)])
    }
    sum_mol    <- sum(mol_est[idx])
    n_spatial  <- sum_mol / eta
    n_tid_tot  <- sum(mols$n_tid[idx])
    n_stoich   <- n_tid_tot / tids_per_mol
    tids_per_gt <- n_tid_tot / max(1L, nu)
    ratio      <- n_spatial / max(1e-6, n_stoich)
    mismatch   <- is.finite(ratio) && (ratio > 2 || ratio < 0.5)
    fused_frac <- mean(mols$fused[idx])
    crowded <- (is.finite(mean_nn) && mean_nn < max(footprint, link_nm + dye_spread)) ||
               fused_frac > 0.3
    cnt_se  <- sqrt(pmax(1, sum_mol)) / eta
    if (crowded) cnt_se <- cnt_se * 1.5
    # CI brackets BOTH estimators (honest: count is only known to ~this range)
    lo <- min(n_spatial - 1.96 * cnt_se, n_stoich); lo <- max(1, lo)
    hi <- max(n_spatial + 1.96 * cnt_se, n_stoich)
    data.frame(
      cluster2      = sid,
      n_units       = nu,                # rendered GT groups (spatial)
      n_fused       = sum(mols$fused[idx]),
      fused_frac    = fused_frac,
      n_gt_est      = n_spatial,
      n_stoich      = n_stoich,
      n_tid_total   = n_tid_tot,
      tids_per_gt   = tids_per_gt,
      stoich_mismatch = mismatch,
      n_gt_lo       = lo,
      n_gt_hi       = hi,
      eta           = eta,
      tids_per_mol  = tids_per_mol,
      mean_nn_nm    = mean_nn,
      crowded       = crowded,
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summ)

  meta <- list(
    level        = unit_level,
    n_fluor      = n_fluor, link_nm = link_nm, dye_spread_nm = dye_spread,
    footprint_nm = footprint, blink_lambda = blink_lambda, eff_pct = eff_pct,
    m_tid        = m_tid, eta = eta, sigma_loc_nm = sigma_loc_nm,
    r_single_nm  = r_single, complex_r_nm = complex_r,
    displace_nm  = meta_disp, guided_merge = isTRUE(guided_merge), use_z = use_z,
    preset_note = preset_note, scoped_to_cluster2 = isTRUE(scope_cluster2) && has_c2
  )
  list(points = points, links = links, summary = summary_df, meta = meta)
}

# Map per-candidate position uncertainty (nm) to a green->amber->red hex ramp.
# lo (best achievable, ~link_nm) maps to green; hi (poor) maps to red.
gt_uncertainty_colors <- function(pos_sd, lo, hi) {
  lo <- suppressWarnings(as.numeric(lo)); hi <- suppressWarnings(as.numeric(hi))
  if (!is.finite(lo)) lo <- 0
  if (!is.finite(hi) || hi <= lo) hi <- lo + 1
  tt   <- pmin(1, pmax(0, (pos_sd - lo) / (hi - lo)))
  tt[!is.finite(tt)] <- 1
  ramp <- grDevices::colorRamp(c("#2ecc71", "#f1c40f", "#e74c3c"))
  m    <- ramp(tt)
  grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
}

## ── Phase B: per-structure morphology / shape classification ──────────────────
# Contour length of a set of points via the DIAMETER of their Euclidean minimum
# spanning tree (longest path in the MST). Robust for CURVED filaments where a
# straight principal-axis extent would underestimate. O(n^2) Prim + tree BFS;
# fine for the tens-to-hundreds of backbone points per structure.
mst_contour_length <- function(P) {
  P <- as.matrix(P); n <- nrow(P)
  if (n < 2L) return(0)
  D <- as.matrix(stats::dist(P))
  inX <- logical(n); inX[1] <- TRUE
  best <- D[1, ]; bestfrom <- rep(1L, n); best[1] <- Inf
  from <- integer(n - 1L); to <- integer(n - 1L); w <- numeric(n - 1L)
  for (it in seq_len(n - 1L)) {
    j <- which.min(ifelse(inX, Inf, best))
    inX[j] <- TRUE
    from[it] <- bestfrom[j]; to[it] <- j; w[it] <- best[j]
    upd <- !inX & (D[j, ] < best)
    best[upd] <- D[j, upd]; bestfrom[upd] <- j
    best[inX] <- Inf
  }
  adj <- vector("list", n)
  for (e in seq_along(from)) {
    adj[[from[e]]] <- c(adj[[from[e]]], e)
    adj[[to[e]]]   <- c(adj[[to[e]]], e)
  }
  far <- function(src) {
    dist <- rep(Inf, n); dist[src] <- 0
    visited <- logical(n); visited[src] <- TRUE; queue <- src
    while (length(queue)) {
      u <- queue[1]; queue <- queue[-1]
      for (e in adj[[u]]) {
        v <- if (from[e] == u) to[e] else from[e]
        if (!visited[v]) { visited[v] <- TRUE; dist[v] <- dist[u] + w[e]; queue <- c(queue, v) }
      }
    }
    dist
  }
  d1 <- far(1L); a <- which.max(d1)
  d2 <- far(a); max(d2[is.finite(d2)])
}

# Collapse points onto a coarse grid (cell = h nm) and return cell centroids.
# Used to thin a backbone before the MST so lateral localization noise doesn't
# inflate the contour length with zig-zag, while curvature is preserved.
.grid_thin <- function(P, h) {
  if (h <= 0 || nrow(P) < 3L) return(P)
  key <- apply(round(P / h), 1L, paste, collapse = "_")
  cen <- lapply(split(seq_len(nrow(P)), key),
                function(ix) colMeans(P[ix, , drop = FALSE]))
  do.call(rbind, cen)
}

# Classify each cluster2's shape and report neutral geometric descriptors.
# EXPLORATORY. Shape from the covariance (gyration-tensor) eigenvalues; contour
# length from the backbone's minimum spanning tree. Not tied to any biology —
# the caller supplies an optional generic subunit spacing if a count is wanted.
#
#   df              clustered SMLM frame (x,y[,z], tid, cluster2)
#   backbone_pts    optional data.frame(x,y[,z],cluster2) of representative points
#                   (e.g. estimated GT candidates) used for the contour length and
#                   inter-point spacing; falls back to per-TID medians when NULL.
#   min_pts         structures with fewer representative points are skipped
#   subunit_spacing_nm  optional generic repeat spacing (nm); when finite and >0
#                   a filament's contour length is divided by it to give an
#                   estimated subunit count (no biological meaning imposed).
#
# Shape from the covariance eigenvalues l1>=l2(>=l3):
#   linearity  = (l1-l2)/l1   planarity = (l2-l3)/l1   sphericity = l3/l1
#   filament (linearity high), ring (hollow, ~constant radius), sheet (planar,
#   3D only), punctum (compact), else diffuse.
# Returns one row per cluster2; empty data.frame when nothing is analysable.
analyze_cluster2_morphology <- function(df, backbone_pts = NULL,
                                        min_pts = 4L, contour_smooth_nm = 10,
                                        subunit_spacing_nm = NA_real_) {
  if (is.null(df) || !nrow(df) || !all(c("x", "y", "tid") %in% names(df)))
    return(data.frame())
  if (!"cluster2" %in% names(df)) return(data.frame())
  use_z <- "z" %in% names(df) && any(is.finite(df$z) & df$z != 0)
  cc    <- if (use_z) c("x", "y", "z") else c("x", "y")

  work <- df[!is.na(df$cluster2) & df$cluster2 > 0L, , drop = FALSE]
  if (!nrow(work)) return(data.frame())
  if (!"z" %in% names(work)) work$z <- 0

  sids <- sort(unique(work$cluster2))
  out <- lapply(sids, function(sid) {
    sub <- work[work$cluster2 == sid, , drop = FALSE]
    # representative points for shape: per-TID medians (robust to loc noise)
    reps <- sub %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = stats::median(x), y = stats::median(y),
                       z = if (use_z) stats::median(z) else 0, .groups = "drop")
    P <- as.matrix(reps[, cc, drop = FALSE])
    if (nrow(P) < min_pts) return(NULL)

    ev  <- eigen(stats::cov(P), symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    l1  <- lam[1]; l2 <- if (length(lam) > 1) lam[2] else 0
    l3  <- if (length(lam) > 2) lam[3] else 0
    if (l1 <= 0) return(NULL)
    linearity  <- (l1 - l2) / l1
    planarity  <- (l2 - l3) / l1
    sphericity <- l3 / l1
    rg <- sqrt(sum(lam))

    # ring test: hollow centre + points at ~constant radius in the best-fit plane
    ctr <- colMeans(P)
    rad <- sqrt(rowSums(sweep(P, 2, ctr)^2))
    r_cv    <- if (mean(rad) > 0) stats::sd(rad) / mean(rad) else Inf
    hollow  <- mean(rad < 0.4 * stats::quantile(rad, 0.95)) < 0.12   # few points near centre
    is_ring <- hollow && r_cv < 0.35 && linearity < 0.5

    # contour + straight extent along principal axis
    axis1  <- ev$vectors[, 1]
    proj   <- as.numeric(P %*% axis1)
    extent <- diff(range(proj))                    # straight-line span
    width  <- 2 * sqrt(max(l2, 0))                 # transverse spread
    bpts <- if (!is.null(backbone_pts) && "cluster2" %in% names(backbone_pts)) {
      bp <- backbone_pts[backbone_pts$cluster2 == sid,
                         intersect(cc, names(backbone_pts)), drop = FALSE]
      if (nrow(bp) >= 2L) as.matrix(bp) else P
    } else P
    # thin lateral noise before measuring contour so it doesn't zig-zag-inflate
    thinned <- .grid_thin(bpts, contour_smooth_nm)
    contour <- mst_contour_length(thinned)
    contour <- max(contour, extent)                # never below the straight span

    # Shape class — dimension-aware. In 2D there is no sheet/sphericity axis, so
    # sphericity (l3/l1) is always 0; classify by linearity + ring only.
    shape <- if (is_ring) "ring"
             else if (linearity > 0.6) "filament"
             else if (use_z && planarity > 0.5 && sphericity < 0.25) "sheet"
             else if (use_z && sphericity > 0.45) "punctum"
             else if (!use_z && linearity < 0.35) "punctum"
             else "diffuse"

    # inter-point spacing among the backbone points (generic "distance between
    # GTs") and, for filaments, mean spacing along the contour + linear density.
    n_bb    <- nrow(bpts)
    nn_med  <- NA_real_
    if (n_bb >= 2L) {
      nn <- tryCatch(RANN::nn2(bpts, bpts, k = 2L)$nn.dists[, 2], error = function(e) NULL)
      if (!is.null(nn)) nn_med <- stats::median(nn[is.finite(nn)])
    }
    spacing_along <- if (shape == "filament" && n_bb >= 2L) contour / (n_bb - 1L) else NA_real_
    lin_density   <- if (shape == "filament" && contour > 0) n_bb / contour * 1000 else NA_real_  # per µm
    n_subunit <- if (shape == "filament" && is.finite(subunit_spacing_nm) &&
                     subunit_spacing_nm > 0) contour / subunit_spacing_nm else NA_real_
    data.frame(
      cluster2      = sid,
      shape         = shape,
      n_rep         = nrow(P),
      n_backbone    = n_bb,
      rg_nm         = rg,
      linearity     = linearity,
      planarity     = planarity,
      sphericity    = sphericity,
      extent_nm     = extent,
      contour_nm    = contour,
      width_nm      = width,
      ring_r_nm     = if (is_ring) stats::median(rad) else NA_real_,
      spacing_nn_nm = nn_med,          # median nearest-neighbour distance between points
      spacing_along_nm = spacing_along, # mean spacing along a filament backbone
      lin_density_um   = lin_density,   # points per µm of contour (filament)
      n_subunit_est    = n_subunit,     # only if a generic subunit spacing was given
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out)) data.frame() else out
}

## ── Radius of Gyration ───────────────────────────────────────────────────────
# Returns a data.frame with one row per cluster: cluster_id, n_traces, rog (nm)
# ── compute_cluster_rog() ─────────────────────────────────────────────────────
# GOAL: radius of gyration  Rg = sqrt(<|p - centroid|^2>)  per cluster — the RMS
# size of a point cloud about its centre.
#
# TWO PARAMETERS, TWO QUESTIONS:
#   cluster_col = what a "cluster" IS (the grouping):
#       "cluster2" -> a median cluster -> Rg measures STRUCTURE SIZE.
#       "cluster"  -> one de-blinked molecule -> Rg measures that molecule's own
#                     blink/label FOOTPRINT (spread of its repeat detections).
#   member_col  = the UNIT reduced to one point before measuring spread:
#       "tid"      -> raw traces (one median point per trace).
#       "cluster"  -> de-blinked molecules (one median point per molecule). Use
#                     this for structure size when de-blinking is ON, so a
#                     multivalent label's several traces count as ONE molecule
#                     rather than inflating Rg with several nearby points.
#
# MINFLUX adaptation: points are MEDIAN member positions, never raw localizations.
# A median of N localizations has SD ~ sigma_loc/sqrt(N), so localization noise is
# already suppressed and Rg reflects real geometry, not measurement scatter. (What
# remains is the label footprint — real between-trace spread — which is exactly the
# thing the "cluster" grouping is meant to report and the "cluster2" grouping is
# meant to look past by de-blinking to one point per molecule.)
#
# SMALL-N BIAS (bias_correct): the centroid is estimated from the SAME points whose
# spread is being measured, so it sits inside the cloud and E[Rg²_obs] = Rg²·(N−1)/N.
# That is a 29% underestimate at N=2, 11% at N=5, 1% at N=50 — i.e. an N-DEPENDENT
# bias, which manufactures a spurious "bigger clusters have bigger Rg" trend whenever
# clusters of different member counts are pooled or plotted against size. Multiplying
# the mean squared distance by N/(N−1) (Bessel, exactly as for a sample variance)
# removes it, giving an unbiased estimate of the underlying structure's Rg. Set
# bias_correct = FALSE for the raw configurational Rg of the observed point set.
#
# Requires >= 2 members per cluster (Rg is undefined for one point). Noise ids (0)
# are dropped for both the grouping and, when members are a cluster level, the
# member. Returns data.frame(cluster_id, n_units, rog).  nm throughout.
compute_cluster_rog <- function(df, cluster_col = "cluster2", member_col = "tid",
                                bias_correct = TRUE) {
  if (!all(c(cluster_col, member_col) %in% names(df))) return(data.frame())
  if (!all(c("x", "y") %in% names(df))) return(data.frame())
  use_z <- "z" %in% names(df)

  # 1) one median point per MEMBER, carrying its cluster id.
  mem <- df %>%
    dplyr::group_by(.data[[member_col]]) %>%
    dplyr::summarise(
      mx = median(x, na.rm = TRUE),
      my = median(y, na.rm = TRUE),
      mz = if (use_z) median(z, na.rm = TRUE) else 0,
      cl = dplyr::first(.data[[cluster_col]]),
      .groups = "drop"
    )
  # Drop noise: a member id of 0 (only meaningful when members are a cluster
  # level), and any point whose cluster id is 0 / NA.
  if (member_col != "tid")
    mem <- mem[is.finite(mem[[member_col]]) & mem[[member_col]] != 0L, , drop = FALSE]
  mem <- mem[is.finite(mem$cl) & mem$cl != 0L, , drop = FALSE]
  if (!nrow(mem)) return(data.frame())

  # 2) per cluster: centroid = mean of member points; Rg = RMS distance to it.
  cents <- mem %>%
    dplyr::group_by(cl) %>%
    dplyr::summarise(cx = mean(mx), cy = mean(my), cz = mean(mz),
                     n_units = dplyr::n(), .groups = "drop")
  # Bessel factor n/(n−1) on the MEAN SQUARED distance removes the N-dependent
  # small-N bias described in the header (no-op as n grows).
  bc <- isTRUE(bias_correct)
  out <- dplyr::left_join(mem, cents, by = "cl") %>%
    dplyr::mutate(r2 = (mx - cx)^2 + (my - cy)^2 + (mz - cz)^2) %>%
    dplyr::group_by(cl, n_units) %>%
    dplyr::summarise(
      rog = sqrt(mean(r2, na.rm = TRUE) *
                   (if (bc && dplyr::first(n_units) > 1L)
                      dplyr::first(n_units) / (dplyr::first(n_units) - 1) else 1)),
      .groups = "drop") %>%
    dplyr::filter(n_units >= 2L, is.finite(rog)) %>%
    dplyr::rename(cluster_id = cl)

  as.data.frame(out)
}

## ── Ripley's L(r) ────────────────────────────────────────────────────────────
# Computes L(r) - r for a set of 2-D or 3-D points.
# Returns a data.frame: r, Kr, Lr, Lr_minus_r
# Uses the naive (uncorrected) Ripley K estimator: K(r) = A * 2*pairs(r) / (n*(n-1)).
# No edge correction is applied; K(r) is underestimated at large r due to border effects.
# Restrict interpretation to r << min(window width, window height) / 2.
compute_ripley_L <- function(x, y, z = NULL, r_seq = NULL, use_3d = FALSE,
                             max_subsample = 2000L, window_ranges = NULL) {
  coords <- cbind(x, y)
  if (use_3d && !is.null(z)) coords <- cbind(x, y, z)

  n <- nrow(coords)
  if (n < 3L) return(data.frame(r = numeric(0), Lr_minus_r = numeric(0)))

  # subsample for large datasets
  if (n > max_subsample) {
    idx <- sample.int(n, max_subsample)
    coords <- coords[idx, , drop = FALSE]
    n <- max_subsample
  }

  # observation window (bounding box). A pre-supplied window_ranges (a 2×d matrix
  # of per-axis min/max, as from apply(coords, 2, range)) lets the observed curve
  # and the Monte-Carlo CSR envelope share one identical window so the edge bias
  # cancels between them.
  rng <- if (!is.null(window_ranges)) window_ranges else apply(coords, 2, range)
  area_vol <- prod(rng[2, ] - rng[1, ])
  if (!is.finite(area_vol) || area_vol <= 0) {
    return(data.frame(r = numeric(0), Lr_minus_r = numeric(0)))
  }

  if (is.null(r_seq)) {
    # r_max is capped by the SHORTEST window axis (edge bias grows once r approaches
    # half the smallest extent). In 3D the z extent is usually far smaller than x/y
    # (thin optical section), so it MUST be included in the cap.
    r_max <- min(diff(rng[, 1]), diff(rng[, 2])) / 4
    if (use_3d && ncol(rng) > 2) r_max <- min(r_max, diff(rng[, 3]) / 4)
    r_seq <- seq(0, r_max, length.out = 60)[-1]
  }

  # pairwise distances (only upper triangle)
  dm <- as.matrix(stats::dist(coords))

  d_vals <- dm[lower.tri(dm)]

  Kr <- vapply(r_seq, function(r) {
    cnt <- sum(d_vals <= r)
    # unbiased K: multiply by 2 for both directions, divide by n^2/area
    (area_vol * 2 * cnt) / (n * (n - 1))
  }, numeric(1L))

  Lr <- if (use_3d && !is.null(z)) {
    (3 * Kr / (4 * pi))^(1/3)
  } else {
    sqrt(Kr / pi)
  }

  data.frame(r = r_seq, Kr = Kr, Lr = Lr, Lr_minus_r = Lr - r_seq)
}

## ── Two-color full pipeline (group assignment + tc1/tc2 clustering) ─────────
## Mirrors plottable_data()'s two-color block so batch analysis can use it
## without requiring the user to have pressed Plot first.
## Returns df with analysis_group (1/2), cluster_tc1, cluster_tc2, and
## cluster/cluster2 aliases for downstream compatibility.
apply_two_color_clustering_pipeline <- function(df, input) {
  if (!nrow(df)) return(df)
  if (!"z" %in% names(df)) df$z <- 0

  # ── 1. Assign two-color groups (DCR or fluo) ──

  # Keep pre-existing valid analysis_group (e.g. from file or prior processing)
  has_preassigned <- "analysis_group" %in% names(df) &&
    any(!is.na(df$analysis_group) & df$analysis_group %in% c(1L, 2L))

  if (!has_preassigned) {
    .dcr_col <- if (isTRUE(input$use_pooled_dcr %||% FALSE) && "dcr_pooled" %in% names(df))
      "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL

    med <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        dcr_med  = if (!is.null(.dcr_col)) stats::median(.data[[.dcr_col]], na.rm = TRUE) else NA_real_,
        fluo_med = if ("fluo" %in% names(df)) stats::median(fluo, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )

    use_fluo <- isTRUE(input$color_by_fluo %||% FALSE) && "fluo" %in% names(df)
    if (use_fluo) {
      med$analysis_group <- dplyr::case_when(
        is.finite(med$fluo_med) & round(med$fluo_med) == 1L ~ 1L,
        is.finite(med$fluo_med) & round(med$fluo_med) == 2L ~ 2L,
        TRUE ~ NA_integer_
      )
    } else if (!is.null(.dcr_col)) {
      # Per-dataset GMM split; fall back to manual thresholds if GMM fails
      gmm_res <- tryCatch(auto_detect_dcr_split(med$dcr_med),
                          error = function(e) list(success = FALSE))
      if (isTRUE(gmm_res$success)) {
        sp <- gmm_res$split
        g1 <- c(0, sp); g2 <- c(sp + 1e-4, 1)
      } else {
        g1 <- input$dcr_group1 %||% c(NA_real_, NA_real_)
        g2 <- input$dcr_group2 %||% c(NA_real_, NA_real_)
      }
      if (length(g1) == 2 && length(g2) == 2 && all(is.finite(c(g1, g2)))) {
        med$analysis_group <- dplyr::case_when(
          is.finite(med$dcr_med) & med$dcr_med >= g1[1] & med$dcr_med <  g1[2] ~ 1L,
          is.finite(med$dcr_med) & med$dcr_med >= g2[1] & med$dcr_med <= g2[2] ~ 2L,
          TRUE ~ NA_integer_
        )
      } else {
        med$analysis_group <- NA_integer_
      }
    } else {
      med$analysis_group <- NA_integer_
    }

    df <- df %>%
      dplyr::select(-dplyr::any_of(c("analysis_group", "dcr_med", "fluo_med"))) %>%
      dplyr::left_join(med %>% dplyr::select(tid, analysis_group, dcr_med, fluo_med), by = "tid")
  }

  df$cluster_tc1 <- NA_integer_
  df$cluster_tc2 <- NA_integer_
  use_3d <- isTRUE(input$mode3d == "3D")

  # ── 2. Primary two-color clustering → cluster_tc1 ──
  if (isTRUE(input$enable_two_color_cluster1 %||% FALSE)) {
    # Centroid de-blinking per channel - mirrors the single-colour path exactly.
    # See deblink_merge_centroids(); raw-point DBSCAN was retired to
    # deadcode_smlm.R.
    df$cluster_tc1 <- 0L
    aw_tc <- if (isTRUE(input$deblink_anisotropic))
      deblink_axis_weights(df, min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT)$w
      else NULL
    thr_tc <- suppressWarnings(as.numeric(input$eps_tc1))
    if (!length(thr_tc) || !is.finite(thr_tc) || thr_tc <= 0) thr_tc <- 12
    meth_tc <- input$deblink_method_tc1 %||% "centroid"
    for (g in c(1L, 2L)) {
      idx <- which(df$analysis_group == g & is.finite(df$x) & is.finite(df$y) & is.finite(df$z))
      if (!length(idx)) next
      cid <- if (identical(meth_tc, "vote"))
        deblink_merge_majority_vote(df[idx, , drop = FALSE],
          f = input$vote_frac_tc1 %||% 0.2, rho = input$vote_rho_tc1 %||% 1.2, w = aw_tc)
        else deblink_merge_centroids(df[idx, , drop = FALSE], thr_tc, aw_tc)
      if (g == 2L) cid <- cid + 10000L
      df$cluster_tc1[idx] <- cid
    }
  }

  # ── 3. Secondary two-color clustering → cluster_tc2 ──
  if (isTRUE(input$enable_two_color_cluster2 %||% FALSE)) {
    tc1_on <- isTRUE(input$enable_two_color_cluster1 %||% FALSE) &&
              any(!is.na(df$cluster_tc1) & df$cluster_tc1 != 0L)
    if (tc1_on) {
      med2 <- df %>%
        dplyr::filter(!is.na(cluster_tc1), cluster_tc1 != 0L) %>%
        dplyr::group_by(cluster_tc1) %>%
        dplyr::summarise(
          mx = stats::median(x, na.rm = TRUE),
          my = stats::median(y, na.rm = TRUE),
          mz = stats::median(z, na.rm = TRUE),
          analysis_group = {
            v <- analysis_group[!is.na(analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
          .groups = "drop"
        ) %>%
        dplyr::filter(analysis_group %in% c(1L,2L), is.finite(mx), is.finite(my), is.finite(mz))
      join_key <- "cluster_tc1"
    } else {
      med2 <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          mx = stats::median(x, na.rm = TRUE),
          my = stats::median(y, na.rm = TRUE),
          mz = stats::median(z, na.rm = TRUE),
          analysis_group = dplyr::first(analysis_group),
          .groups = "drop"
        ) %>%
        dplyr::filter(analysis_group %in% c(1L,2L), is.finite(mx), is.finite(my), is.finite(mz))
      join_key <- "tid"
    }
    med2$cluster_tc2 <- 0L
    for (g in c(1L, 2L)) {
      idx <- which(med2$analysis_group == g)
      if (length(idx) < (input$minpts_tc2 %||% 5L)) next
      pts <- if (use_3d) as.matrix(med2[idx, c("mx","my","mz")]) else as.matrix(med2[idx, c("mx","my")])
      cid <- run_cluster2_grouping(
        pts,
        method        = input$cluster2_method_tc2 %||% "dbscan",
        eps           = input$eps_tc2 %||% 100,
        minPts        = input$minpts_tc2 %||% 5L,
        xi            = input$cluster2_xi_tc2 %||% 0.05,
        prefilter     = isTRUE(input$cluster2_prefilter_tc2),
        prefilter_r   = input$cluster2_prefilter_r_tc2,
        prefilter_min = input$cluster2_prefilter_min_tc2,
        prune         = isTRUE(input$cluster2_prune_tc2),
        prune_tol     = input$cluster2_prune_tol_tc2 %||% 2
      )
      if (g == 2L) cid[cid > 0] <- cid[cid > 0] + 10000L
      med2$cluster_tc2[idx] <- cid
    }
    df <- df %>%
      dplyr::select(-dplyr::any_of("cluster_tc2")) %>%
      dplyr::left_join(med2 %>% dplyr::select(dplyr::all_of(join_key), cluster_tc2), by = join_key)
  }

  # ── 4. Alias to cluster/cluster2 for downstream compatibility ──
  df$cluster  <- df$cluster_tc1
  df$cluster2 <- df$cluster_tc2
  df
}

## ── Two-color cross-channel NN ──────────────────────────────────────────────
compute_cross_nn <- function(df, unit = "tid", k = 1L) {
  if (!"analysis_group" %in% names(df)) return(data.frame())
  df1 <- dplyr::filter(df, analysis_group == 1L)
  df2 <- dplyr::filter(df, analysis_group == 2L)
  if (!nrow(df1) || !nrow(df2)) return(data.frame())

  get_pts <- function(d, u) {
    if (u == "tid") {
      d %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(x = median(x, na.rm = TRUE),
                         y = median(y, na.rm = TRUE),
                         z = median(z, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    } else {
      cl_col <- if ("cluster" %in% names(d) && any(!is.na(d$cluster) & d$cluster > 0L)) "cluster" else NULL
      if (is.null(cl_col)) return(data.frame(x = numeric(), y = numeric(), z = numeric()))
      d %>%
        dplyr::filter(!is.na(.data[[cl_col]]) & .data[[cl_col]] > 0L) %>%
        dplyr::group_by(.data[[cl_col]]) %>%
        dplyr::summarise(x = median(x, na.rm = TRUE),
                         y = median(y, na.rm = TRUE),
                         z = median(z, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    }
  }

  pts1 <- get_pts(df1, unit)
  pts2 <- get_pts(df2, unit)
  if (!nrow(pts1) || !nrow(pts2)) return(data.frame())

  k1 <- min(k, nrow(pts2))
  k2 <- min(k, nrow(pts1))

  # Use 2-D (x, y) distances so that the CSR null model (which is 2-D) is consistent.
  run_nn <- function(query, ref, kk) {
    nn <- RANN::nn2(as.matrix(ref[, c("x","y")]),
                    as.matrix(query[, c("x","y")]),
                    k = kk)
    as.vector(nn$nn.dists[, kk])
  }

  bbox_area <- diff(range(c(pts1$x, pts2$x), na.rm = TRUE)) *
               diff(range(c(pts1$y, pts2$y), na.rm = TRUE))
  rbind(
    data.frame(direction = "Group 1 → Group 2", distance = run_nn(pts1, pts2, k1),
               n_ref = nrow(pts2), area = bbox_area),
    data.frame(direction = "Group 2 → Group 1", distance = run_nn(pts2, pts1, k2),
               n_ref = nrow(pts1), area = bbox_area)
  )
}

## ── Bivariate Ripley's L12(r) ───────────────────────────────────────────────
compute_ripley_L12 <- function(x1, y1, z1 = NULL, x2, y2, z2 = NULL,
                               r_seq = NULL, use_3d = FALSE,
                               max_subsample = 2000L, window_ranges = NULL) {
  if (use_3d && (!is.null(z1)) && (!is.null(z2))) {
    pts1 <- cbind(x1, y1, z1)
    pts2 <- cbind(x2, y2, z2)
  } else {
    pts1 <- cbind(x1, y1)
    pts2 <- cbind(x2, y2)
    use_3d <- FALSE
  }

  n1 <- nrow(pts1); n2 <- nrow(pts2)
  if (n1 < 2L || n2 < 2L) return(data.frame())

  if (n1 > max_subsample) pts1 <- pts1[sample(n1, max_subsample), , drop = FALSE]
  if (n2 > max_subsample) pts2 <- pts2[sample(n2, max_subsample), , drop = FALSE]
  n1 <- nrow(pts1); n2 <- nrow(pts2)

  all_pts <- rbind(pts1, pts2)
  # per-axis extents; a pre-supplied window_ranges (2×d min/max matrix) makes the
  # observed curve and its CSR envelope share one identical window.
  rng <- if (!is.null(window_ranges)) (window_ranges[2, ] - window_ranges[1, ])
         else apply(all_pts, 2, function(v) diff(range(v, na.rm = TRUE)))
  if (use_3d) {
    area_vol <- prod(rng)
    if (!is.finite(area_vol) || area_vol == 0) return(data.frame())
  } else {
    area_vol <- prod(rng[1:2])
    if (!is.finite(area_vol) || area_vol == 0) return(data.frame())
  }

  dmat <- as.matrix(dist(rbind(pts1, pts2)))[seq_len(n1), (n1 + 1):(n1 + n2)]

  if (is.null(r_seq) || length(r_seq) == 0) {
    # Cap by the SHORTEST window axis, not the longest: edge bias is set by the
    # dimension that runs out first (in 3D that is almost always z). Matches
    # compute_ripley_L's cap.
    r_max <- min(if (use_3d) rng else rng[1:2]) / 4
    r_seq <- seq(0, r_max, length.out = 60L)[-1]  # drop r=0 (trivially L12(0)-0=0)
  }

  Kr <- vapply(r_seq, function(r) {
    cnt <- sum(dmat <= r)
    area_vol * cnt / (n1 * n2)
  }, numeric(1L))

  Lr <- if (use_3d) (3 * Kr / (4 * pi))^(1/3) else sqrt(Kr / pi)
  data.frame(r = r_seq, Kr = Kr, Lr = Lr, Lr_minus_r = Lr - r_seq)
}

## ── Pair Cross-Correlation Function g₁₂(r) ──────────────────────────────────
## g12(r) = (area × N_pairs_in_annulus) / (n1 × n2 × annulus_area)
## g12 = 1 under CSR; > 1 = attraction; < 1 = repulsion.
compute_pcf12 <- function(x1, y1, z1 = NULL, x2, y2, z2 = NULL,
                           r_seq = NULL, use_3d = FALSE, max_subsample = 2000L,
                           window_ranges = NULL) {
  if (use_3d && !is.null(z1) && !is.null(z2)) {
    pts1 <- cbind(x1, y1, z1); pts2 <- cbind(x2, y2, z2)
  } else {
    pts1 <- cbind(x1, y1);     pts2 <- cbind(x2, y2); use_3d <- FALSE
  }
  n1 <- nrow(pts1); n2 <- nrow(pts2)
  if (n1 < 2L || n2 < 2L) return(data.frame())

  if (n1 > max_subsample) pts1 <- pts1[sample(n1, max_subsample), , drop = FALSE]
  if (n2 > max_subsample) pts2 <- pts2[sample(n2, max_subsample), , drop = FALSE]
  n1 <- nrow(pts1); n2 <- nrow(pts2)

  all_pts <- rbind(pts1, pts2)
  # per-axis extents; a pre-supplied window_ranges (2×d min/max matrix) makes the
  # observed curve and its CSR envelope share one identical window.
  rng     <- if (!is.null(window_ranges)) (window_ranges[2, ] - window_ranges[1, ])
             else apply(all_pts, 2, function(v) diff(range(v, na.rm = TRUE)))
  area_vol <- if (use_3d) prod(rng) else prod(rng[1:2])
  if (!is.finite(area_vol) || area_vol == 0) return(data.frame())

  dmat <- as.matrix(dist(rbind(pts1, pts2)))[seq_len(n1), (n1 + 1):(n1 + n2)]

  if (is.null(r_seq) || !length(r_seq)) {
    # Cap by the SHORTEST window axis (see compute_ripley_L / compute_ripley_L12).
    r_max <- min(if (use_3d) rng else rng[1:2]) / 4
    r_seq <- seq(0, r_max, length.out = 60L)
  }
  dr <- if (length(r_seq) > 1L) mean(diff(r_seq)) else r_seq[1L]

  g12 <- vapply(r_seq, function(r) {
    lo <- max(0, r - dr / 2); hi <- r + dr / 2
    cnt <- sum(dmat > lo & dmat <= hi)
    shell <- if (use_3d) {
      (4 / 3) * pi * (hi^3 - lo^3)
    } else {
      pi * (hi^2 - lo^2)
    }
    if (shell == 0 || !is.finite(shell)) return(NA_real_)
    area_vol * cnt / (n1 * n2 * shell)
  }, numeric(1L))

  data.frame(r = r_seq, g12 = g12)
}

## ── Monte-Carlo CSR envelope for L(r) / L₁₂(r) / g₁₂(r) ──────────────────────
## Generates `nsim` complete-spatial-randomness (CSR) patterns of `n` points
## (and `n2` for the cross statistics) uniformly inside the SAME window as the
## observed pattern (`window_ranges`, a 2×d matrix of per-axis min/max), and runs
## the identical estimator on each. Returns pointwise lo/hi quantiles (default
## 2.5% / 97.5% → a 95% band; a two-sided α≈0.05 test needs nsim ≥ 39) and the
## mean null curve, on the supplied `r_seq`.
##
## Because the sims share the observed window, the (uncorrected) edge bias that
## drags L(r)−r negative at large r is present in the envelope too and cancels
## when the observed curve is read against the band — this envelope is the
## interpretive tool, superseding an explicit edge correction.
compute_ripley_envelope <- function(n, window_ranges, r_seq, use_3d = FALSE,
                                     nsim = 99L, estimator = c("L", "L12", "pcf"),
                                     n2 = NULL, probs = c(0.025, 0.975)) {
  estimator <- match.arg(estimator)
  if (is.null(r_seq) || !length(r_seq)) return(NULL)
  if (is.null(window_ranges) || any(!is.finite(window_ranges))) return(NULL)
  nsim <- max(1L, as.integer(nsim))
  n    <- as.integer(n)
  if (n < 3L) return(NULL)
  if (estimator != "L") {
    if (is.null(n2)) return(NULL)
    n2 <- as.integer(n2)
    if (n2 < 2L) return(NULL)
  }

  d    <- ncol(window_ranges)
  lo_v <- window_ranges[1, ]; hi_v <- window_ranges[2, ]
  is3d <- isTRUE(use_3d) && d >= 3L

  gen <- function(m) {
    list(x = stats::runif(m, lo_v[1], hi_v[1]),
         y = stats::runif(m, lo_v[2], hi_v[2]),
         z = if (is3d) stats::runif(m, lo_v[3], hi_v[3]) else NULL)
  }

  nr  <- length(r_seq)
  mat <- matrix(NA_real_, nrow = nr, ncol = nsim)
  cap <- if (estimator == "L") n else max(n, n2)

  for (s in seq_len(nsim)) {
    p1 <- gen(n)
    col <- if (estimator == "L") {
      cv <- compute_ripley_L(p1$x, p1$y, p1$z, r_seq = r_seq, use_3d = is3d,
                             max_subsample = cap, window_ranges = window_ranges)
      cv$Lr_minus_r
    } else if (estimator == "L12") {
      p2 <- gen(n2)
      cv <- compute_ripley_L12(p1$x, p1$y, p1$z, p2$x, p2$y, p2$z,
                               r_seq = r_seq, use_3d = is3d,
                               max_subsample = cap, window_ranges = window_ranges)
      cv$Lr_minus_r
    } else {
      p2 <- gen(n2)
      cv <- compute_pcf12(p1$x, p1$y, p1$z, p2$x, p2$y, p2$z,
                          r_seq = r_seq, use_3d = is3d,
                          max_subsample = cap, window_ranges = window_ranges)
      cv$g12
    }
    if (length(col) == nr) {
      mat[, s] <- col
    } else if (length(col) > 0L) {
      k <- min(length(col), nr)
      mat[seq_len(k), s] <- col[seq_len(k)]
    }
  }

  lo <- apply(mat, 1, stats::quantile, probs = probs[1], na.rm = TRUE)
  hi <- apply(mat, 1, stats::quantile, probs = probs[2], na.rm = TRUE)
  data.frame(r = r_seq, lo = as.numeric(lo), hi = as.numeric(hi),
             null_mean = rowMeans(mat, na.rm = TRUE))
}

## ── 2D KDE density co-localization map ──────────────────────────────────────
## Returns KDE grids for Ch1 and Ch2 on a shared spatial grid,
## plus Pearson r of the two density surfaces.
## Optional honesty add-ons:
##   bw_sweep    – vector of bandwidths; returns $bw_sweep = data.frame(bw,
##                 pearson_r) so the bandwidth-dependence of Pearson r is visible.
##   costes_nsim – >0 builds a randomization null for Pearson r this many times and
##                 returns a one-sided p-value (fraction of null r ≥ observed) plus
##                 the null mean/SD.
##   null_mode   – WHICH null. Measured false-positive rates are in the note below;
##                 read it before quoting any p-value from this function.
##     "shift" (default) – TOROIDAL RANDOM SHIFT of channel 2 (Lotwick & Silverman
##                 1982). Channel 2 is translated by a random offset with wrap-around,
##                 so its OWN clustering, density and point count are preserved exactly
##                 and only its spatial REGISTRATION with channel 1 is destroyed. This
##                 is the null that isolates cross-channel dependence, and it is the
##                 same random-shift logic already used for the bivariate statistics.
##     "csr"     – channel 2 replaced by uniform random points, i.e. a test against
##                 complete spatial randomness rather than against independent
##                 placement. It is NOT the Costes et al. 2004 test (theirs preserves
##                 the randomized channel's autocorrelation by block randomization).
##
## ⚠ CALIBRATION — measured, not assumed (2026-07-30). Simulated 30× per scenario,
## independent placement (truth = NOT co-localized), so a calibrated test should
## reject at p<0.05 about 5% of the time:
##   • Channels clustered but spread over the WHOLE window:
##         shift 3% · csr 3%   → both fine.
##   • Channels sharing ONE CELL that fills a small part of the window (the realistic
##     SMLM case), placed independently within it:
##         shift 53% · csr 93% → BOTH anti-conservative; csr severely so.
## So "csr" really does mostly report "both channels are non-uniform", but the
## random-shift default is only a PARTIAL fix: wrapping Ch2 moves it outside the cell,
## so the null still contrasts "Ch2 inside the cell" with "Ch2 partly outside it". The
## root cause is that BOTH nulls randomize over the bounding box rather than over the
## actual support (the cell mask). Until randomization is mask-aware, treat a low
## p-value here as WEAK evidence whenever the signal does not fill the field, and get
## significance from g12(r) against its CSR envelope instead. Power is not the problem:
## on a genuinely co-localized positive control the shift null fired 10/10, median
## p = 0.02.
compute_kde_coloc <- function(x1, y1, x2, y2, bandwidth = NULL, n_grid = 100L,
                              bw_sweep = NULL, costes_nsim = 0L,
                              null_mode = c("shift", "csr")) {
  null_mode <- match.arg(null_mode)
  all_x <- c(x1, x2); all_y <- c(y1, y2)
  xlim  <- range(all_x, na.rm = TRUE); ylim <- range(all_y, na.rm = TRUE)
  if (diff(xlim) == 0 || diff(ylim) == 0) return(NULL)

  if (is.null(bandwidth)) {
    n   <- length(all_x)
    bw  <- 0.9 * min(stats::sd(all_x, na.rm = TRUE),
                     stats::IQR(all_x, na.rm = TRUE) / 1.34,
                     na.rm = TRUE) * n^(-0.2)
    bandwidth <- if (is.finite(bw) && bw > 0) bw else diff(xlim) / 20
  }

  k1 <- MASS::kde2d(x1, y1, h = bandwidth, n = n_grid, lims = c(xlim, ylim))
  k2 <- MASS::kde2d(x2, y2, h = bandwidth, n = n_grid, lims = c(xlim, ylim))

  pearson_r <- stats::cor(as.vector(k1$z), as.vector(k2$z))

  out <- list(x = k1$x, y = k1$y, z1 = k1$z, z2 = k2$z,
              pearson_r = pearson_r, bandwidth = bandwidth)

  # ── Bandwidth sweep: r(bw) ─────────────────────────────────────────────────
  if (!is.null(bw_sweep) && length(bw_sweep)) {
    bw_sweep <- sort(unique(bw_sweep[is.finite(bw_sweep) & bw_sweep > 0]))
    if (length(bw_sweep)) {
      sweep_r <- vapply(bw_sweep, function(h) {
        ka <- MASS::kde2d(x1, y1, h = h, n = n_grid, lims = c(xlim, ylim))
        kb <- MASS::kde2d(x2, y2, h = h, n = n_grid, lims = c(xlim, ylim))
        stats::cor(as.vector(ka$z), as.vector(kb$z))
      }, numeric(1L))
      out$bw_sweep <- data.frame(bw = bw_sweep, pearson_r = sweep_r)
    }
  }

  # ── Spatial randomization test (see null_mode in the header) ──────────────
  costes_nsim <- as.integer(costes_nsim %||% 0L)
  if (isTRUE(costes_nsim > 0L)) {
    n2 <- length(x2)
    wx <- diff(xlim); wy <- diff(ylim)
    null_r <- vapply(seq_len(costes_nsim), function(i) {
      if (identical(null_mode, "shift")) {
        # Toroidal shift: wrap channel 2 inside the SAME window. Preserves its own
        # clustering and n exactly; only the cross-channel registration is broken.
        sx <- stats::runif(1, 0, wx); sy <- stats::runif(1, 0, wy)
        rx <- xlim[1] + ((x2 - xlim[1] + sx) %% wx)
        ry <- ylim[1] + ((y2 - ylim[1] + sy) %% wy)
      } else {
        rx <- stats::runif(n2, xlim[1], xlim[2])
        ry <- stats::runif(n2, ylim[1], ylim[2])
      }
      kb <- MASS::kde2d(rx, ry, h = bandwidth, n = n_grid, lims = c(xlim, ylim))
      stats::cor(as.vector(k1$z), as.vector(kb$z))
    }, numeric(1L))
    # one-sided p: how often randomized co-localization matches/exceeds observed
    out$costes_p         <- (sum(null_r >= pearson_r, na.rm = TRUE) + 1) /
                            (sum(is.finite(null_r)) + 1)
    out$costes_null_mean <- mean(null_r, na.rm = TRUE)
    out$costes_null_sd   <- stats::sd(null_r, na.rm = TRUE)
    out$costes_nsim      <- costes_nsim
    out$null_mode        <- null_mode
    out$null_label       <- if (identical(null_mode, "shift"))
      "random toroidal shift of Ch2 (structure-preserving; Lotwick & Silverman 1982)"
    else "uniform CSR resampling of Ch2 (tests against complete spatial randomness)"
  }

  out
}

## ── Grid-based density correlation (multi-scale) ─────────────────────────────
## For each cell size in `cell_sizes`, bins Ch1 and Ch2 molecules into a 2D
## grid and returns Pearson r between per-cell counts — a multi-scale
## co-localization measure.

make_helpers <- function(mode) {
  if (identical(mode, "3D")) {
    list(
      assign_group         = function(df, input) assign_analysis_group(df, input),
      equal_limits         = function(df) { lim <- get_equal_limits(df$x, df$y, df$z); list(x=lim$x, y=lim$y, z=lim$z) },
      distance_to_fit      = function(...) calculate_distance_to_fit(...),         # 3D
      add_distance_lines   = function(fig, med) add_distance_lines_to_plot(fig, med) # 3D
    )
  } else {
    list(
      assign_group         = function(df, input) assign_analysis_group_2d(df, input),
      equal_limits         = function(df) get_equal_limits_2d_fn(df$x, df$y),
      distance_to_fit      = function(...) calculate_distance_to_fit(...),
      add_distance_lines   = function(fig, med) add_distance_lines_to_plot(fig, med) # 3D
    )
  }
}
