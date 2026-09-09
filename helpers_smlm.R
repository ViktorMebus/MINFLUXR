## ─────────────────────────────────────────────────────────────────────────────
## helpers_smlm.R  -  SMLM-only helper functions
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   All the computation behind the SMLM (imaging) tab. Unlike the Tracking tab -
##   which asks how molecules MOVE - this tab asks where molecules ARE and how
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
##      resolve_channels ............... THE one place a channel is defined
##                                       (none / 2-colour split / multiplex)
##      channels_in_scope .............. which channels a run covers
##      deblink_threshold_for_channel .. per-channel merge distance (k x sigma_ch)
##      apply_clustering_pipeline ...... THE main entry point (primary + cluster2),
##                                       loops over channels; never merges across
##      run_cluster2_grouping, .cluster2_core, .prune_sparse_members
##      (membership refinements: core-only / core-distance prune / min members)
##
##   5. NEAREST NEIGHBOUR
##      nn_index_pairs ................. within / between / global channel scope
##      build_nn_edges, pair_passes_scope, compute_batch_nn, compute_cross_nn
##      safe_nn2, safe_k
##
##   6. SHAPE FITTING (the biggest block - see its own banner below)
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
##      compute_ripley_L, compute_ripley_L12, compute_pcf12, compute_pcf,
##      compute_ripley_envelope ........ CSR envelope by randomisation
##      compute_kde_coloc .............. 2D KDE + Pearson co-localization
##
##   8b. LABELLING-DENSITY NORMALISATION (opt-in; see the block banner)
##      unit_centroids ................. one point per trace / de-blinked molecule
##      estimate_support ............... edge- and dark-region-aware support area
##      nn_local_intensity ............. support-free intensity from k-th NN
##      compute_unit_density ........... the per-file density row
##      thin_units_to_p ................ p-thinning at the UNIT level
##
##   8c. LOCAL DENSITY & HOTSPOTS (see the block banner)
##      compute_local_enhancement ...... local density / this cell's mean density
##      compute_dispersion_index ....... variance-to-mean of quadrat counts
##      calibrate_dbscan_eps ........... eps a CSR pattern in this mask would not cluster
##      compute_foci_table ............. one row per detected focus (n, Rg, aniso)
##      hotspot_summary_row ............ the per-file row the table and plots read
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
# round cluster renders as an ellipse - which would corrupt every visual shape
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
  # Thin wrapper. The DCR/fluo split itself lives in resolve_channels(), which
  # is the only place a channel is defined; this used to be a second copy of it
  # (and server_smlm had a third). Callers here want a per-TRACE map, so that is
  # what it returns - one row per tid with analysis_group and the medians.
  if (!is.data.frame(df) || !nrow(df)) return(df[0, , drop = FALSE])
  if (!"z" %in% names(df)) df$z <- 0
  d <- resolve_channels(df, utils::modifyList(as.list(input), list(enable_two_color = TRUE)))
  if (!nrow(d)) return(d[0, c("tid"), drop = FALSE])
  d %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      x_med = stats::median(x, na.rm = TRUE),
      y_med = stats::median(y, na.rm = TRUE),
      z_med = stats::median(z, na.rm = TRUE),
      dcr_med  = if ("dcr_med"  %in% names(d)) dplyr::first(dcr_med)  else NA_real_,
      fluo_med = if ("fluo_med" %in% names(d)) dplyr::first(fluo_med) else NA_real_,
      analysis_group = dplyr::first(analysis_group),
      .groups = "drop")
}
# Helper: auto-detect 2-color DCR split via 2-component GMM (mclust)
# Returns a list with $success, $split, $mean1, $mean2, $sd1, $sd2, $n_traces
# On failure returns $success = FALSE and $reason string.
# GOAL: split traces into two colour channels on their DCR (detection channel
# ratio) by fitting a 2-component Gaussian mixture, AND quantify how reliable that
# split is.
#
# WHY THE UNCERTAINTY MATTERS (and why it is not symmetric with a registration
# error): ratiometric MINFLUX has no inter-channel offset - both dyes are excited
# by the same laser through the same beam path, so the two channels share one
# coordinate frame. The residual two-colour error is therefore NOT a shift, it is
# MIS-ASSIGNMENT: the DCR distributions of the two dyes overlap, so a fraction of
# traces lands on the wrong side of the threshold. A mis-assigned molecule sits by
# construction right next to its true neighbours, so crosstalk pushes every
# cross-channel statistic (NN, L12, g12) UP at short r - toward FALSE
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
                                "). DCR distribution appears unimodal - is this 2-color data?")))
  pi1  <- fit$parameters$pro[ord[1]]
  pi2  <- fit$parameters$pro[ord[2]]

  # Extract SDs - mclust stores variance differently depending on model type
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
  # more useful per-channel numbers - what proportion of everything you CALL
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
  # comparison of G=1 vs G=2 does NOT fix it - it also preferred two components on
  # that file (dBIC +16.9) and on a synthetic skewed-unimodal control (+21.1), while
  # missing a genuine two-colour pair at d' = 1 (-9.6). BIC is picking up skewness,
  # not bimodality, so it is deliberately NOT used as a gate here. The only valid
  # basis for splitting is the EXPERIMENTAL DESIGN plus a look at the DCR histogram.
  dprime <- abs(mu2 - mu1) / sqrt((sd1^2 + sd2^2) / 2)

  # Per-trace posterior of the assigned component. Traces whose best posterior is
  # below `post_cut` are genuinely ambiguous - they are the ones that would move if
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
#   group     1 or 2 - the component with the higher posterior (matches the split
#             threshold by construction, since the threshold IS the crossing point)
#   post      posterior of the assigned component, in [0.5, 1]. 0.5 = the molecule
#             sits exactly on the split and is a coin flip; ~1 = unambiguous.
#   edge      TRUE when post < post_cut, i.e. the molecule is in the overlap band
#             and would switch channel under a small shift of the threshold.
#   log_odds  log(post_assigned / post_other) - the same information on an additive
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
  # would give 0/0. Fall back to the side of the threshold there - the assignment
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

# =============================================================================
# resolve_channels() - THE ONE PLACE A CHANNEL IS DEFINED
# =============================================================================
#
# GOAL
#   Turn whatever grouping the user has set up - nothing, a 2-colour DCR/fluo
#   split, or a multiplex build - into ONE integer `channel` column, so that
#   everything downstream has exactly one notion of "these points must not be
#   grouped together".
#
# THE INVARIANT THIS COLUMN EXISTS TO ENFORCE
#   Grouping never crosses a channel. De-blinking, trace-ID plurality
#   reassignment and median clustering all run strictly inside one channel.
#   DISTANCE may cross a channel, but only through an explicitly pairwise tool
#   (cross-NN, L12, g12, co-localization) - never through a merge.
#
#   "Never merge two colours" and "never merge two multiplex rounds" are the
#   same statement about the same column. Before this function existed they were
#   two separate implementations, and the 2-colour one had silently dropped the
#   multiplex branch - so a dataset with both active de-blinked ACROSS multiplex
#   channels inside each colour.
#
# THE THREE SOURCES
#   none ....... 1 channel. Everything behaves as it did before channels existed.
#   2-colour ... 2 channels, assigned PER TRACE from median DCR (or `fluo` 1/2).
#                A trace is one channel or the other, never split.
#   multiplex .. 2-7 channels, assigned PER LOCALIZATION (file overlay, or TIM
#                window). Carried on `mplex_channel` by the multiplex builder.
#
#   In practice a sample is one of the three. If a DCR split and multiplex are
#   BOTH active the channel becomes their product (colour x round), which is the
#   only reading under which both rules hold at once. That case is defensive, not
#   expected.
#
# WHAT IT WRITES
#   channel ........ integer 1..N, never NA (unassignable rows go to channel 1
#                    when there is nothing to split on; see below)
#   channel_label .. display name ("A"/"B", or the multiplex label)
#   channel_color .. hex, reused by the plot AND by the per-channel histograms
#   analysis_group . 1L/2L/NA - the PAIR VIEW, kept because every cross-channel
#                    statistic (L12, g12, cross-NN) is pairwise by mathematics
#                    and already consumes this column. With N > 2 it is filled
#                    only for the selected pair.
#   two_color_group  alias of `channel` while the source is the 2-colour split,
#                    kept so existing rendering code is untouched.
#
#   attr(df, "channel_table") - one row per channel (channel, label, color, n).
#   attr(df, "channel_source") - "none" | "two_color" | "mplex" | "mplex_x_color"
#
# UNASSIGNED TRACES ARE DROPPED, NOT POOLED
#   A trace whose median DCR falls in neither group range is NA, and NA cannot be
#   given a channel without guessing. Those rows are removed. Keeping them would
#   put them in some channel, and any choice of channel is a merge across a
#   boundary the user asked us to respect.
# =============================================================================
resolve_channels <- function(df, input) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  if (!"z" %in% names(df)) df$z <- 0

  two_color_on <- isTRUE(input$enable_two_color %||% FALSE)
  mplex_on     <- "mplex_channel" %in% names(df) && any(!is.na(df$mplex_channel))

  # ── Colour split (per trace) ────────────────────────────────────────────
  # Reuses assign_analysis_group()'s logic but forced into DCR/fluo mode: when
  # the 2-colour switch is on, the split IS the grouping, so `grouping_mode` has
  # nothing left to decide.
  col_grp <- NULL
  if (two_color_on) {
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
      med$col_grp <- dplyr::case_when(
        is.finite(med$fluo_med) & round(med$fluo_med) == 1L ~ 1L,
        is.finite(med$fluo_med) & round(med$fluo_med) == 2L ~ 2L,
        TRUE ~ NA_integer_
      )
    } else if (!is.null(.dcr_col)) {
      # Per-dataset GMM split, falling back to the manual ranges if it fails.
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
        med$col_grp <- dplyr::case_when(
          is.finite(med$dcr_med) & med$dcr_med >= g1[1] & med$dcr_med <  g1[2] ~ 1L,
          is.finite(med$dcr_med) & med$dcr_med >= g2[1] & med$dcr_med <= g2[2] ~ 2L,
          TRUE ~ NA_integer_
        )
      } else {
        med$col_grp <- NA_integer_
      }
    } else {
      med$col_grp <- NA_integer_
    }

    df <- df %>%
      dplyr::select(-dplyr::any_of(c("dcr_med", "fluo_med", "col_grp"))) %>%
      dplyr::left_join(med, by = "tid")
    col_grp <- df$col_grp
    df$col_grp <- NULL

    # Traces the split could not place are dropped rather than pooled.
    keep <- !is.na(col_grp)
    if (!all(keep)) {
      df      <- df[keep, , drop = FALSE]
      col_grp <- col_grp[keep]
    }
    if (!nrow(df)) return(df)
  }

  # COLOUR IS NOT SET HERE. This function runs inside the cached
  # filter -> thin -> cluster chain, and the colour pickers are not part of what
  # makes a prepared frame: changing a channel's colour must not invalidate the
  # clusters and force every file to be re-clustered. Reading them here also made
  # the cache key dishonest - it did not include them, so a recoloured channel
  # kept serving the old hex. Colour is attached at render time instead, by
  # attach_channel_colors(). See .prep_key_selftest() in server_smlm.R.
  #
  # `mplex_color` is different and stays: it is a DATA column written by the
  # multiplex builder, not a live input.

  # ── Build the channel column ────────────────────────────────────────────
  if (mplex_on && two_color_on) {
    # Defensive product case. Both rules hold simultaneously only here.
    mch <- df$mplex_channel
    mch[is.na(mch)] <- 0L
    key <- interaction(mch, col_grp, drop = TRUE, lex.order = TRUE)
    df$channel <- as.integer(key)
    lab_m <- if ("mplex_label" %in% names(df)) df$mplex_label else as.character(mch)
    df$channel_label <- paste0(lab_m, "/", c("A", "B")[col_grp])
    df$channel_color <- NA_character_
    src <- "mplex_x_color"

  } else if (mplex_on) {
    mch <- df$mplex_channel
    mch[is.na(mch)] <- 0L
    lev <- sort(unique(mch))
    df$channel <- match(mch, lev)
    df$channel_label <- if ("mplex_label" %in% names(df)) as.character(df$mplex_label)
                        else paste0("Ch", df$channel)
    df$channel_color <- if ("mplex_color" %in% names(df)) as.character(df$mplex_color)
                        else mplex_default_colors(length(lev))[df$channel]
    src <- "mplex"

  } else if (two_color_on) {
    df$channel       <- col_grp
    df$channel_label <- c("A", "B")[col_grp]
    df$channel_color <- NA_character_
    src <- "two_color"

  } else {
    df$channel       <- 1L
    df$channel_label <- "All"
    df$channel_color <- NA_character_
    src <- "none"
  }

  # ── Pair view for the cross-channel statistics ──────────────────────────
  # With exactly two channels the pair is unambiguous. With more, the caller
  # picks one (assign_analysis_group_from_mplex), and until it does the pairwise
  # tools have nothing to run on - which is correct, not a failure.
  # A caller that has ALREADY chosen the pair wins - assign_analysis_group_from_mplex()
  # drops every row outside the chosen pair and labels the survivors 1/2, and which
  # of the two is "1" is the user's choice of A and B. Re-deriving it from channel
  # order here would silently swap the channels on any pair whose labels do not sort
  # the same way, which flips the sign of every asymmetric cross-channel statistic.
  preassigned <- "analysis_group" %in% names(df) &&
    all(c(1L, 2L) %in% df$analysis_group[!is.na(df$analysis_group)])

  n_ch <- length(unique(df$channel))
  if (identical(src, "two_color")) {
    df$analysis_group  <- col_grp
    df$two_color_group <- col_grp
  } else if (preassigned) {
    df$two_color_group <- NA_integer_
  } else if (n_ch == 2L) {
    df$analysis_group  <- as.integer(factor(df$channel))
    df$two_color_group <- NA_integer_
  } else {
    if (!"analysis_group" %in% names(df)) df$analysis_group <- NA_integer_
    df$two_color_group <- NA_integer_
  }

  ct <- df %>%
    dplyr::group_by(channel) %>%
    dplyr::summarise(label = dplyr::first(channel_label),
                     color = dplyr::first(channel_color),
                     n     = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(channel)

  attr(df, "channel_table")  <- as.data.frame(ct)
  attr(df, "channel_source") <- src
  df
}

# =============================================================================
# attach_channel_colors() - COLOUR, APPLIED AT RENDER TIME
# =============================================================================
#
# GOAL
#   Fill `channel_color` (and refresh the `channel_table` attribute) from the
#   live colour pickers, for a frame that already has `channel`.
#
# WHY IT IS SEPARATE FROM resolve_channels()
#   resolve_channels() runs inside prepared_df()'s cached
#   filter -> thin -> cluster chain. Anything it reads from `input` becomes part
#   of what defines a prepared frame, and must therefore be in the cache key
#   (.PREP_KEY_INPUTS) or the cache serves stale data.
#
#   Colour should be neither. Putting the pickers in the key would re-cluster
#   every loaded file whenever someone nudges a colour - expensive, and for a
#   purely cosmetic change. Leaving them out of the key while still reading them
#   was the other failure: a recoloured channel kept its old hex until the cache
#   happened to turn over. The way out of that trade-off is not to read them in
#   the cached chain at all.
#
#   So: `channel` and `channel_label` are DATA and are cached; `channel_color` is
#   PRESENTATION and is applied here, cheaply, every render.
#
# ARGUMENTS
#   df    : a frame carrying `channel` (from resolve_channels)
#   input : Shiny input list. Reads color_group1 / color_group2 only.
#
# RETURNS
#   df with channel_color filled and attr "channel_table" refreshed. Rows that
#   already carry a colour (multiplex, whose builder writes `mplex_color` into
#   the data) keep it.
# =============================================================================
attach_channel_colors <- function(df, input) {
  if (!is.data.frame(df) || !nrow(df) || !"channel" %in% names(df)) return(df)
  if (!"channel_color" %in% names(df)) df$channel_color <- NA_character_

  # Multiplex carries its colour as a data column - respect it.
  if ("mplex_color" %in% names(df)) {
    keep <- is.na(df$channel_color) & !is.na(df$mplex_color)
    df$channel_color[keep] <- as.character(df$mplex_color[keep])
  }

  need <- is.na(df$channel_color)
  if (any(need)) {
    ids <- sort(unique(df$channel[!is.na(df$channel)]))
    pal <- if (length(ids) <= 2L)
      c(input$color_group1 %||% "#1f77b4", input$color_group2 %||% "#d62728")
      else mplex_default_colors(length(ids))
    idx <- match(df$channel, ids)
    df$channel_color[need] <- pal[((idx[need] - 1L) %% length(pal)) + 1L]
  }

  ct <- attr(df, "channel_table")
  if (is.data.frame(ct) && nrow(ct)) {
    m <- df$channel_color[match(ct$channel, df$channel)]
    ct$color <- m
    attr(df, "channel_table") <- ct
  }
  df
}

# Channels a clustering run should cover, honouring `cluster_group_scope`.
# There is deliberately NO pooled option: pooling across channels is the thing
# the channel column exists to prevent. "within" runs each channel separately;
# "only:<k>" restricts to one channel for inspection.
channels_in_scope <- function(df, input) {
  if (!"channel" %in% names(df)) return(list(all = TRUE, ids = NA_integer_))
  ids <- sort(unique(df$channel[!is.na(df$channel)]))
  scope <- input$cluster_group_scope %||% "within"
  if (grepl("^only:", scope)) {
    k <- suppressWarnings(as.integer(sub("^only:", "", scope)))
    if (is.finite(k) && k %in% ids) ids <- k
  }
  ids
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
# =============================================================================
# compute_sd_medians() - THE REPORTED SPREAD SUMMARY
# =============================================================================
#
# GOAL
#   Collapse the per-trace spreads into the numbers the UI reports: a centre, a
#   robust width, and the interquartile range, per axis.
#
# REPORT vs DECIDE - WHY robust DEFAULTS TO FALSE HERE
#   The per-trace sigma fed in is the plain n-1 SD (see sd_by_tid), because that
#   is the one number with an outside reference. Robustness enters ACROSS traces
#   (median and MAD), which is where pyMINFLUX puts it too, and where it costs
#   no comparability. The corrected-MAD branch stays reachable for the internal
#   callers that need a breakdown-resistant per-trace sigma.
#
# RETURNS
#   list(x, y[, z]) medians as before, plus mad_x/mad_y[/mad_z] and
#   q1_x/q3_x/... quartiles. NA where an axis cannot be measured.
# =============================================================================
compute_sd_medians <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, robust = FALSE) {
  tbl <- trace_spread_table(df, min_n = min_n, robust = robust)
  has_z <- "sd_z" %in% names(tbl) && any(is.finite(tbl$sd_z))
  if (!nrow(tbl)) {
    # Consistent empty result so the callers don't crash
    return(list(x = NA_real_, y = NA_real_, z = if ("sd_z" %in% names(tbl)) NA_real_ else NULL,
                mad_x = NA_real_, mad_y = NA_real_, mad_z = NA_real_,
                q1_x = NA_real_, q3_x = NA_real_, q1_y = NA_real_, q3_y = NA_real_,
                q1_z = NA_real_, q3_z = NA_real_, n_traces = 0L))
  }

  # One axis -> centre, robust width, quartiles. Every statistic drops NA first
  # so a trace that is too short on one axis does not blank the whole summary.
  ax_stats <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(list(med = NA_real_, mad = NA_real_, q1 = NA_real_, q3 = NA_real_))
    q <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    list(med = stats::median(v),
         # constant = 1.4826 matches scipy's scale = 0.67449 used by pyMINFLUX
         mad = stats::mad(v, constant = 1.4826, na.rm = TRUE),
         q1  = q[1], q3 = q[2])
  }

  sx <- ax_stats(tbl$sd_x); sy <- ax_stats(tbl$sd_y)
  sz <- if (has_z) ax_stats(tbl$sd_z) else NULL

  res <- list(
    x = sx$med, y = sy$med,
    mad_x = sx$mad, mad_y = sy$mad,
    q1_x = sx$q1, q3_x = sx$q3,
    q1_y = sy$q1, q3_y = sy$q3,
    n_traces = nrow(tbl)
  )
  if (!is.null(sz)) {
    res$z <- sz$med; res$mad_z <- sz$mad; res$q1_z <- sz$q1; res$q3_z <- sz$q3
  }
  res
}
trace_spread_table <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, robust = FALSE) {
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

  if (isTRUE(robust)) {
    # Robust branch: per-axis spread comes from sd_by_tid (corrected MAD), and
    # the overall spread is rebuilt from those axes rather than measured as an
    # RMSD. RMSD squares the residuals, so one stray localization would dominate
    # spread_rmsd even though every sd_ column had shrugged it off - the table
    # would then contradict itself. For Gaussian residuals sqrt(sum of per-axis
    # sigma^2) has the same expectation as the RMSD it replaces (sqrt(d)*sigma),
    # so the column keeps its meaning and its scale.
    out <- sd_by_tid(df, min_n = min_n, robust = TRUE)
    if (!nrow(out)) return(data.frame(
      tid = character(), sd_x = numeric(), sd_y = numeric(), sd_z = numeric(),
      spread_rmsd = numeric(), sd_max = numeric()))
    if (!has_z) out$sd_z <- NA_real_
    sq <- out$sd_x^2 + out$sd_y^2 + if (has_z) out$sd_z^2 else 0
    out$spread_rmsd <- sqrt(sq)
    return(out %>%
      dplyr::filter(.data$n_locs >= min_n) %>%
      dplyr::mutate(sd_max = pmax(sd_x, sd_y, sd_z, na.rm = TRUE)) %>%
      dplyr::select(-dplyr::any_of("n_locs")))
  }

  base <- df %>% dplyr::group_by(tid)

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
# =============================================================================
# grp_median() / grp_mad() - PER-GROUP MEDIAN AND MAD IN ONE RADIX SORT
# =============================================================================
#
# GOAL
#   Per-trace median and MAD fast enough to sit in a reactive that refires on
#   every filter change.
#
# WHY NOT stats::mad() UNDER group_by()
#   Measured on 600k localizations / 60k traces, three axes:
#     dplyr + stats::sd    1.5 s
#     dplyr + stats::mad   8.6 s     <- 5.9x, unusable in a reactive
#     this                 0.7 s     <- 0.5x, i.e. FASTER than the sd it replaces
#   stats::mad is a generic that re-sorts for each of its two medians and pays
#   full argument dispatch per group. Sorting ONCE by (group, value) makes every
#   group's values both contiguous and ordered, so the median is an index
#   lookup, not a search. Same trick as the cleaning-pipeline radix fix.
#
#   Verified identical to stats::mad on 60k traces (max abs difference 0).
#
# ARGUMENTS
#   g  : integer group code, 1..ng (NOT character - that is the slow path)
#   v  : values, NAs allowed (dropped per group, like na.rm = TRUE)
#   ng : number of groups
#
# RETURNS  numeric(ng), NA for groups with no non-NA value.
# =============================================================================
grp_median <- function(g, v, ng) {
  ok <- !is.na(v)
  if (!all(ok)) { g <- g[ok]; v <- v[ok] }
  o  <- order(g, v, method = "radix")
  gs <- g[o]; vs <- v[o]

  cnt   <- tabulate(gs, nbins = ng)
  end   <- cumsum(cnt)
  start <- end - cnt + 1L

  res <- rep(NA_real_, ng)
  hit <- cnt > 0L
  # even n -> mean of the two central values; odd n -> lo and hi coincide
  lo  <- start[hit] + (cnt[hit] - 1L) %/% 2L
  hi  <- start[hit] + cnt[hit] %/% 2L
  res[hit] <- (vs[lo] + vs[hi]) / 2
  res
}

grp_mad <- function(g, v, ng) {
  m <- grp_median(g, v, ng)
  1.4826 * grp_median(g, abs(v - m[g]), ng)
}

# =============================================================================
# mad_finite_sample_factor() - MAKE MAD AN UNBIASED sigma AT SHORT TRACE LENGTHS
# =============================================================================
#
# GOAL
#   Remove MAD's small-sample downward bias so a spread threshold in nm means
#   the same thing for a 6-localization trace as for a 40-localization one.
#
# WHY THIS IS NOT OPTIONAL HERE
#   The 1.4826 constant makes MAD consistent for sigma only ASYMPTOTICALLY. At
#   the trace lengths this app actually sees, E[MAD]/sigma is:
#     n = 4  0.735      n = 8  0.887      n = 14 0.940
#     n = 6  0.841      n = 12 0.931      n = 20 0.959
#   Uncorrected, the same nm threshold is therefore far more permissive for
#   SHORT traces - a length-dependent bias applied invisibly inside the filter,
#   which is exactly the failure mode deblink_minpts() documents for minPts.
#
# THE FACTORS
#   1 / E[MAD] under a Gaussian, simulated at 2e5 replicates per n. These
#   reproduce the published Croux & Rousseeuw (1992) correction table to within
#   0.5% (their n = 8: 1.129, here 1.1269), which is the standard reference for
#   this correction. Above n = 20 the usual n / (n - 0.8) form takes over; it is
#   continuous with the table at the join (1.0417 vs 1.0427).
#
# ARGUMENTS  n : integer vector of per-trace localization counts
# RETURNS    numeric multiplier, 1 where n is not usable
# =============================================================================
.MAD_SMALL_N_FACTOR <- c(
  1.1948, 1.4862, 1.3596, 1.2181, 1.1889, 1.1373, 1.1269, 1.0992, 1.0956,
  1.0798, 1.0769, 1.0662, 1.0633, 1.0558, 1.0552, 1.0486, 1.0482, 1.0426,
  1.0427)                                            # n = 2 .. 20

mad_finite_sample_factor <- function(n) {
  n   <- as.integer(n)
  out <- rep(1, length(n))
  small <- !is.na(n) & n >= 2L & n <= 20L
  out[small] <- .MAD_SMALL_N_FACTOR[n[small] - 1L]
  big <- !is.na(n) & n > 20L
  out[big] <- n[big] / (n[big] - 0.8)
  out
}

# =============================================================================
# sd_by_tid() - PER-TRACE, PER-AXIS SPREAD
# =============================================================================
#
# GOAL
#   One spread number per trace per axis, for the Stage-2 spread filter, the
#   spread histograms and the de-blink axis weights.
#
# WHY THE DEFAULT IS THE PLAIN SD
#   This is the number the app REPORTS, and it is the one quantity here that has
#   an outside reference: pyMINFLUX computes exactly this per trace ID
#   (processor/_processor.py, calculate_statistics_on):
#       sx = df_grouped["x"].std()   # pandas, ddof = 1
#   Keeping it identical is what makes MINFLUXR's precision numbers comparable
#   with pyMINFLUX's. pyMINFLUX has NO spread-based trace filter and NO
#   per-localization outlier removal, so nothing downstream of this number is
#   constrained by that comparison - only the number itself.
#
# WHY sd's ZERO BREAKDOWN POINT IS NOT A PROBLEM HERE ANY MORE
#   One stray localization sets an sd. That used to delete whole traces over two
#   bad points. The fix is trace_point_outliers(), which removes the stray
#   localizations BEFORE this runs, rather than swapping the estimator: the
#   contamination is dealt with where it lives instead of being averaged over.
#   With trimming on there is one spread statistic in the whole app.
#
# ARGUMENTS
#   df     : localization table (tid, x, y, optionally z)
#   min_n  : below this, spread is left NA (undefined) rather than guessed
#   robust : FALSE = pyMINFLUX-matching stats::sd (default and user-facing).
#            TRUE  = finite-sample-corrected MAD. INTERNAL machinery only -
#            trace_point_outliers() uses it to set the trimming scale, where a
#            breakdown-resistant sigma is required and comparability is
#            irrelevant. Not exposed as a user choice: two competing readings of
#            "trace spread" in the UI is exactly the confusion to avoid.
#
# RETURNS  data.frame(tid, n_locs, sd_x, sd_y, sd_z), tid in its input type.
# =============================================================================
sd_by_tid <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, robust = FALSE) {
  need_cols <- intersect(c("tid","x","y","z"), names(df))
  if (!all(c("tid","x","y") %in% need_cols)) return(data.frame())

  min_n <- max(2L, as.integer(min_n %||% MIN_LOCS_FOR_SPREAD_DEFAULT))
  has_z <- "z" %in% need_cols

  if (!isTRUE(robust)) {
    d <- df[, need_cols, drop = FALSE]
    d$x <- as.numeric(d$x); d$y <- as.numeric(d$y)
    if (has_z) d$z <- as.numeric(d$z)
    return(d %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        n_locs = dplyr::n(),
        sd_x = if (dplyr::n() >= min_n) stats::sd(x, na.rm = TRUE) else NA_real_,
        sd_y = if (dplyr::n() >= min_n) stats::sd(y, na.rm = TRUE) else NA_real_,
        sd_z = if (has_z && dplyr::n() >= min_n) stats::sd(z, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ))
  }

  # Integer trace codes: grp_median sorts on them, and a character sort here is
  # the ~25x penalty the cleaning pipeline already had to remove once.
  f  <- factor(as.character(df$tid))
  gi <- as.integer(f)
  ng <- nlevels(f)
  if (!ng) return(data.frame())

  n_locs <- tabulate(gi, nbins = ng)
  # The correction is per trace, because it depends on that trace's own n.
  fac <- mad_finite_sample_factor(n_locs)

  ax <- function(nm) {
    if (!nm %in% need_cols) return(rep(NA_real_, ng))
    grp_mad(gi, as.numeric(df[[nm]]), ng) * fac
  }

  # tid must come back in the type it went in as. factor() stringifies it, and
  # returning character where the classic branch returns integer breaks every
  # downstream join on tid ("can't join <integer> with <character>"). The
  # level -> value mapping survives the conversion because as.integer(f) indexes
  # into levels(f), so the pairing is positional, not alphabetical.
  tid_out <- levels(f)
  proto   <- df$tid
  if (is.factor(proto))       tid_out <- factor(tid_out, levels = levels(proto))
  else if (is.integer(proto)) tid_out <- as.integer(tid_out)
  else if (is.numeric(proto)) tid_out <- as.numeric(tid_out)

  out <- data.frame(tid = tid_out, n_locs = n_locs,
                    sd_x = ax("x"), sd_y = ax("y"), sd_z = ax("z"),
                    stringsAsFactors = FALSE)

  # Short traces get NA, not a guess - an unreliable spread must not drive the
  # filter. Matches the non-robust branch above.
  short <- out$n_locs < min_n
  out$sd_x[short] <- NA_real_
  out$sd_y[short] <- NA_real_
  out$sd_z[short] <- NA_real_
  out
}
# =============================================================================
# trace_point_outliers() - FLAG SINGLE LOCALIZATIONS THAT DO NOT FIT THEIR TRACE
# =============================================================================
#
# GOAL
#   Find the individual localizations that sit far from the rest of their OWN
#   trace, so a 12-point trace is not deleted over 2 stray points. Flags rows;
#   it never deletes and never touches whole traces.
#
# THE QUESTION IT ASKS, AND WHY THAT CHOICE MATTERS
#   "Is this point unlike the REST OF THIS TRACE" - not "unlike the dataset".
#   A uniformly broad trace (drift, two molecules, a fiducial) then loses
#   NOTHING here, because all of its points agree with each other; deciding its
#   fate stays with the spread filter and the rogue rules. Only a trace that is
#   tight apart from a few strays loses those strays. Measured on real data
#   (collaborator-cond02-B2, k = 4): a dataset-wide sigma flags 7.6% of all
#   localizations, this rule flags 3.4%, and the difference is almost entirely
#   points inside legitimately broad traces.
#
# THE SCALE: PER-TRACE SIGMA SHRUNK TOWARD THE DATASET
#   A trace's own robust sigma is the right scale but is noisy at n ~ 10: used
#   raw it flags 2.4% of points in CLEAN simulated data, because a trace that
#   happens to draw a small MAD then rejects its own good points. The dataset
#   sigma is stable but wrong for broad traces. So they are blended,
#       sigma_eff = w * sigma_trace + (1 - w) * sigma_pooled,   w = n / (n + 6)
#   which is the usual shrinkage form: short traces borrow the dataset's
#   stability, long traces trust themselves. Measured on clean simulation this
#   costs 0.88% false flags (vs 2.43% per-trace, 0.22% pooled) while still
#   catching 100% of injected contaminants at k = 4.
#
#   sigma is robust (corrected MAD) at BOTH levels on purpose. An sd-based scale
#   would be inflated by the very outliers it is meant to find, so the test
#   would hide them - the estimator is internal machinery, not the reported
#   number (see sd_by_tid).
#
# THE DISTANCE
#   Per-axis whitened, r^2 = sum over axes of ((v - median) / sigma_eff)^2, so
#   the anisotropy of MINFLUX precision (z is not x) does not leak into the
#   threshold. Same metric as apply_deblink_weights.
#
# WHY THE GAUSSIAN k DOES NOT PREDICT THE REJECTION RATE
#   Under a Gaussian null, k = 4 rejects essentially nothing. On real data it
#   flags ~3%. The residuals are genuinely heavier-tailed than Gaussian, so the
#   only honest rejection estimate is the one MEASURED on the loaded data -
#   which is what trim_summary() reports rather than a theoretical number.
#
# PROTECTING SHORT TRACES
#   Below min_trim_n a trace is never trimmed: at n = 5 a single point is 20% of
#   the evidence and the trace median that defines "far" is itself barely
#   determined. Trimming also never takes a trace below keep_n retained points.
#
# ARGUMENTS
#   df         : localization table (tid, x, y, optionally z)
#   k          : reject beyond this many sigma (user-set, default 4)
#   min_trim_n : traces shorter than this are left untouched
#   keep_n     : never trim a trace below this many retained localizations
#   w0         : shrinkage half-weight; n = w0 splits trace and pooled evenly
#
# RETURNS  logical vector, one per ROW of df. TRUE = outlier.
# =============================================================================
trace_point_outliers <- function(df, k = 4, min_trim_n = 6L,
                                 keep_n = MIN_LOCS_FOR_SPREAD_DEFAULT, w0 = 6) {
  n <- nrow(df)
  if (!is.data.frame(df) || !n || !all(c("tid", "x", "y") %in% names(df)))
    return(rep(FALSE, max(0L, n)))
  k <- suppressWarnings(as.numeric(k))
  if (!length(k) || !is.finite(k) || k <= 0) return(rep(FALSE, n))

  axes <- intersect(c("x", "y", "z"), names(df))
  f  <- factor(as.character(df$tid))
  gi <- as.integer(f); ng <- nlevels(f)
  n_tr <- tabulate(gi, nbins = ng)

  # Robust scale at both levels: per trace, and pooled over traces long enough
  # for their own estimate to mean anything.
  sg  <- sd_by_tid(df, min_n = 2L, robust = TRUE)
  sg  <- sg[match(levels(f), as.character(sg$tid)), , drop = FALSE]
  fac <- n_tr / (n_tr + w0)

  r2 <- numeric(n)
  for (ax in axes) {
    v   <- as.numeric(df[[ax]])
    col <- paste0("sd_", ax)
    per <- if (col %in% names(sg)) sg[[col]] else rep(NA_real_, ng)
    pooled <- stats::median(per[n_tr >= max(4L, min_trim_n)], na.rm = TRUE)
    if (!is.finite(pooled) || pooled <= 0) pooled <- stats::median(per, na.rm = TRUE)
    if (!is.finite(pooled) || pooled <= 0) next          # axis carries no spread
    use <- ifelse(is.finite(per) & per > 0, per, pooled)
    sig <- fac * use + (1 - fac) * pooled
    sig[!is.finite(sig) | sig <= 0] <- pooled
    ctr <- grp_median(gi, v, ng)
    r2  <- r2 + ((v - ctr[gi]) / sig[gi])^2
  }
  flag <- sqrt(r2) > k
  flag[is.na(flag)] <- FALSE

  # Short traces are exempt entirely.
  flag[n_tr[gi] < as.integer(min_trim_n)] <- FALSE

  # Never trim a trace below keep_n retained points. A trace that would breach
  # that is not "a good trace with a few strays" - it is a trace this rule
  # cannot describe, so it is left whole for the trace-level filters to judge.
  n_flag <- tabulate(gi[flag], nbins = ng)
  breach <- which(n_tr - n_flag < as.integer(keep_n))
  if (length(breach)) flag[gi %in% breach] <- FALSE

  flag
}

# =============================================================================
# trim_summary() - WHAT THE TRIM ACTUALLY DID, FOR THE INFO BOX
# =============================================================================
# GOAL
#   Report the trim as a measured outcome on THIS dataset, not as a theoretical
#   rejection rate - see trace_point_outliers() on why the Gaussian number is
#   wrong here by a factor of ten.
#
#   The per-trace trimmed FRACTION is the composite-trace diagnostic: a few
#   percent means stray localizations, a large fraction means the trace was
#   never one molecule (simulated, a second molecule contributing 1/12 of the
#   points trims 9%, 2/12 trims 18%, 4/12 trims 36%, 6/12 trims 86%). Reported
#   as a count of traces past a fraction, not applied as a rule, because the
#   fraction is n-dependent and the threshold would be arbitrary.
#
# ARGUMENTS  df : the table the flags were computed on;  flag : that logical
# RETURNS    list of scalars for renderUI
# =============================================================================
trim_summary <- function(df, flag) {
  n <- length(flag)
  if (!n) return(list(n_pts = 0L, n_flag = 0L, pct = 0, n_tr_touched = 0L,
                      pct_tr = 0, n_tr_heavy = 0L, n_tr = 0L))
  tid  <- as.character(df$tid)
  per  <- tapply(flag, tid, mean)
  list(
    n_pts        = n,
    n_flag       = sum(flag),
    pct          = 100 * mean(flag),
    n_tr         = length(per),
    n_tr_touched = sum(per > 0),
    pct_tr       = 100 * mean(per > 0),
    n_tr_heavy   = sum(per >= 0.30)
  )
}

# =============================================================================
# cov_ellipse_params() - THE ERROR ELLIPSE OF A 2x2 COVARIANCE MATRIX
# =============================================================================
#
# GOAL
#   Turn a covariance matrix into a drawable, ROTATED ellipse plus the two
#   numbers that describe the cloud's shape independently of the axes.
#
# WHY ROTATION MATTERS (the thing σx / σy alone cannot express)
#   An axis-aligned ellipse has two parameters and cannot represent a TILTED
#   cloud. Elongation along a diagonal inflates σx and σy by the SAME amount, so
#   σx ≈ σy and the ellipse renders as a circle over a visibly stretched cloud -
#   all the elongation lives in cov(x,y), which an axis-aligned drawing discards.
#
#   The same blind spot afflicts max(σx,σy)/min(σx,σy) as an anisotropy metric:
#   a 45° cloud scores ~1.0, i.e. "isotropic". Since drift within a trace is
#   DIRECTIONAL, a tilted cloud is the expected signature of drift - precisely
#   the artifact the precision panel exists to catch.
#
#   The eigen-decomposition is the fix: eigenvectors give the true major/minor
#   directions, sqrt(eigenvalues) the semi-axes. Their ratio is rotation-
#   invariant, so it detects elongation at ANY angle.
#
# ARGUMENTS  S : symmetric 2x2 covariance matrix;  k : how many σ to scale to
# RETURNS
#   list(ra, rb, theta, ratio, rho, angle_deg) - ra/rb are semi-major/minor in
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

# =============================================================================
# geometric_median() - A ROBUST CENTRE THAT DOES NOT DEPEND ON THE AXES
# =============================================================================
#
# GOAL
#   The point minimising the sum of DISTANCES to a set of points, found by
#   Weiszfeld iteration. Used to summarise a cloud of displacement vectors where
#   the answer is a direction, so the estimator must not care how x and y happen
#   to be oriented.
#
# WHY NOT A PER-AXIS MEDIAN
#   Taking median(x) and median(y) separately is NOT rotation-equivariant:
#   rotate the frame, take the medians, rotate back, and you get a different
#   vector. Measured on displacement clouds like the drift plot's, that
#   disagreement reaches 6.8 nm - which for a plot whose entire content is a
#   DIRECTION means the arrow would depend on the microscope's axis convention.
#   The mean and the geometric median are both exactly equivariant (0.000 nm).
#
# WHY NOT THE MEAN
#   The mean has no breakdown point. On these clouds the contamination is real:
#   a "molecule" can hold two genuinely different fluorophores, which puts a
#   displacement of tens of nm into an otherwise ~10 nm cloud. Measured RMSE
#   against a known 10 nm shift, jitter sd 12 nm, n = 20 per bin:
#
#     contamination        mean    per-axis median   geometric median
#     none                 4.55        5.69              4.95
#     10% wrong molecule   5.70        6.01              5.54
#     25% wrong molecule   7.59        6.81              6.07
#     10% gross outliers  15.65        6.08              5.60
#
#   The geometric median is best or tied everywhere contamination exists and
#   costs ~9% efficiency when there is none, which is the right trade for a
#   panel whose whole purpose is judging a real sample.
#
# ARGUMENTS
#   M    : matrix, one row per point
#   iter : maximum Weiszfeld iterations
#   tol  : stop when the estimate moves less than this
#
# RETURNS  numeric vector of length ncol(M); the mean for degenerate input.
# =============================================================================
geometric_median <- function(M, iter = 200L, tol = 1e-9) {
  if (!is.matrix(M)) M <- as.matrix(M)
  M <- M[stats::complete.cases(M), , drop = FALSE]
  if (!nrow(M)) return(rep(NA_real_, ncol(M)))
  if (nrow(M) <= 2L) return(colMeans(M))
  y <- colMeans(M)
  for (k in seq_len(iter)) {
    d <- sqrt(rowSums(sweep(M, 2L, y)^2))
    # A point sitting exactly on the current estimate would divide by zero; the
    # floor keeps the step finite and costs nothing at these scales.
    d[d < 1e-12] <- 1e-12
    y2 <- colSums(M / d) / sum(1 / d)
    if (!all(is.finite(y2))) return(colMeans(M))
    if (sqrt(sum((y2 - y)^2)) < tol) return(y2)
    y <- y2
  }
  y
}

# GOAL: the outline of a rotated ellipse as plain x/y vectors.
# Drawn as a TRACE rather than a layout shape: a shape cannot be switched off
# from the legend, and these rings are exactly what a reader wants to toggle.
ellipse_xy <- function(ra, rb, theta, cx = 0, cy = 0, n = 181L) {
  if (!all(is.finite(c(ra, rb, theta, cx, cy))) || ra <= 0 || rb <= 0) return(NULL)
  t  <- seq(0, 2 * pi, length.out = n)
  ct <- cos(theta); st <- sin(theta)
  list(x = cx + ra * cos(t) * ct - rb * sin(t) * st,
       y = cy + ra * cos(t) * st + rb * sin(t) * ct)
}


# =============================================================================
# density_hdr_2d() - SMOOTH DENSITY + HIGHEST-DENSITY-REGION CONTOURS
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
#   ~9% noise at the peak but ~26% at the 2-sigma radius and ~90% at 3 sigma -
#   i.e. exactly the outer contours the eye reads as "the shape of the cloud"
#   are almost pure shot noise, which is what makes the plot look ragged and
#   under-smoothed. Re-binning cannot fix it: coarser bins blur the core,
#   finer bins make the tails worse. A kernel estimate is the fix, because each
#   point contributes to a neighbourhood and the estimate is continuous.
#
# WHY HDR LEVELS
#   A contour at "density = 0.004" means nothing on its own. A contour holding
#   half the points means something a reader can act on. Drawn together with a
#   fitted ellipse enclosing the SAME mass, the empirical contour and the
#   theoretical ring coincide when the cloud really is Gaussian and separate
#   visibly where it is not. That comparison is the whole diagnostic, and it is
#   unavailable with count-based levels.
#
#   The masses are round numbers on purpose. Any pair works, so they are chosen
#   to be readable (50%, 90%) and the sigma multiple is derived from them - for
#   a 2D Gaussian the k-sigma ellipse holds 1 - exp(-k^2/2), so 50% sits at
#   1.177 sigma and 90% at 2.146 sigma. Choosing the masses the other way round
#   (0.393, 0.865, i.e. exactly 1 and 2 sigma) makes the model ring easier to
#   label and the plot harder to read; the reader is the one to optimise for.
#
# BANDWIDTH
#   Scott's rule for 2D (sigma * n^(-1/6)) with a floor at the grid pitch -
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
#   list(x, y, z, bw_sd, mass_in_frame, paths) where paths is a list of
#   list(prob=, d=, x=, y=) contour outlines, or NULL if the cloud is unusable.
#   A requested mass that extends beyond the frame is skipped rather than drawn
#   wrong. Draw the outlines with plotly::add_paths(), NOT add_lines() -
#   add_lines() sorts by x and turns a closed contour into a filled zigzag.
# =============================================================================
density_hdr_2d <- function(a, b, lim, sigma, n_grid = 160L,
                           probs = c(0.50, 0.90), max_pts = 20000L) {
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
    # exists here - skip it rather than draw the frame edge and imply otherwise.
    if (!is.finite(in_frame) || in_frame < p) next
    idx <- which(cm >= p)[1]
    lev <- srt[idx]
    if (!is.finite(lev) || lev <= 0) next
    cls <- tryCatch(grDevices::contourLines(k$x, k$y, k$z, levels = lev),
                    error = function(e) NULL)
    for (cl in cls) {
      if (length(cl$x) < 3L) next
      # Plain coordinates: the caller draws these as a TRACE (so the legend can
      # toggle them), which needs x/y, not an SVG path string.
      paths[[length(paths) + 1L]] <- list(prob = p, d = lev, x = cl$x, y = cl$y)
    }
  }
  list(x = k$x, y = k$y, z = k$z, bw_sd = bw_sd, mass_in_frame = in_frame, paths = paths)
}

# =============================================================================
# robust_median_cov() - ONE VOTE PER TRACE, BUT STILL TILT-AWARE
# =============================================================================
# GOAL
#   Combine many per-trace covariance matrices into one robust matrix by taking
#   the element-wise median, so a single 2000-localization fiducial cannot
#   dominate the shape the way it dominates a pooled estimate.
#
#   Element-wise medians are not guaranteed to stay positive semi-definite, so
#   the off-diagonal is shrunk if it would exceed what the diagonal supports -
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

# =============================================================================
# add_sigma_iqr_lines() - THE QUARTILE PAIR THAT FLANKS THE MEDIAN LINE
# =============================================================================
#
# GOAL
#   Draw thin lines at the 25th and 75th percentile of a spread distribution,
#   flanking the median line that add_sigma_vline() already draws.
#
# WHY
#   Matches pyMINFLUX's add_median_line(..., add_iqr = True), which puts dashed
#   quartile lines either side of its median marker on the sx/sy/sz histograms.
#   It also earns its place independently: median +/- nMAD is symmetric, while
#   the distribution of per-trace spread is right-skewed (a spread cannot go
#   below zero but has no ceiling). The quartiles show that asymmetry, which the
#   +/- label cannot.
#
# ARGUMENTS
#   fig    : plotly figure to add to
#   v      : the per-trace spread values the median line was drawn from
#   y_top  : same panel ceiling used for the median line
#   color  : match the median line this pair belongs to
#
# RETURNS  the figure, unchanged when the quartiles are not finite.
# =============================================================================
add_sigma_iqr_lines <- function(fig, v, y_top, color) {
  v <- as.numeric(v); v <- v[is.finite(v)]
  if (length(v) < 4L || !is.finite(y_top) || y_top <= 0) return(fig)
  q <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(q)) || q[1] == q[2]) return(fig)

  for (xq in q) {
    fig <- fig %>% plotly::add_segments(
      x = xq, xend = xq,
      y = 0, yend = y_top,
      line = list(color = color, width = 1, dash = "dot"),
      inherit = FALSE, hoverinfo = "skip", showlegend = FALSE
    )
  }
  fig
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

# --- median marker lines for any plotly histogram -------------------------
# GOAL: put the summary statistic ON the distribution it summarises, so a shifted
# or bimodal histogram can be read against a number instead of by eye.
#
# Draws one dashed vertical line per series, spanning the full panel height
# (yref = "paper", so it does not need the histogram's y scale) with the median
# value as a label at the top. Lines are drawn as layout SHAPES rather than
# traces, so they add nothing to the legend and cannot be toggled off with a
# series - which also means they do NOT follow legend clicks; hide the panel's
# median lines with the panel's own checkbox instead.
#
# Input : fig     - a plotly object
#         meds    - named numeric vector, name = series label, value = median
#         colors  - named character vector of line colors (same names as `meds`)
#         digits  - significant digits in the label
#         unit    - appended to the label ("nm" etc.)
# Output: the plotly object with shapes + annotations appended
add_median_vlines <- function(fig, meds, colors = NULL, digits = 3, unit = "") {
  meds <- meds[is.finite(meds)]
  if (!length(meds)) return(fig)
  nms  <- names(meds) %||% rep("", length(meds))
  # Stagger the labels so two nearby medians do not print on top of each other.
  ytop <- 1 - 0.07 * (seq_along(meds) - 1L)

  shapes <- lapply(seq_along(meds), function(i) {
    col <- if (!is.null(colors) && nms[i] %in% names(colors)) unname(colors[nms[i]]) else "#333333"
    list(type = "line", xref = "x", yref = "paper",
         x0 = meds[[i]], x1 = meds[[i]], y0 = 0, y1 = 1,
         line = list(color = col, width = 2, dash = "dash"))
  })
  anns <- lapply(seq_along(meds), function(i) {
    col <- if (!is.null(colors) && nms[i] %in% names(colors)) unname(colors[nms[i]]) else "#333333"
    list(x = meds[[i]], y = ytop[i], xref = "x", yref = "paper",
         yanchor = "bottom", xanchor = "left",
         text = sprintf("med %s%s", signif(meds[[i]], digits),
                        if (nzchar(unit)) paste0(" ", unit) else ""),
         showarrow = FALSE, font = list(color = col, size = 11))
  })

  fig %>% plotly::layout(
    shapes      = c(fig$x$layout$shapes %||% list(), shapes),
    annotations = c(fig$x$layout$annotations %||% list(), anns))
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
  # shorter traces get NA (undefined spread) so their noisy estimate never
  # drives the filter - they are kept rather than judged on an unreliable
  # spread. Same SD as the Stage-2 filter, on the same trimmed point set.
  if (isTRUE(input$trim_outliers)) {
    fl <- tryCatch(trace_point_outliers(df, k = input$trim_k %||% 4, keep_n = min_n),
                   error = function(e) NULL)
    if (!is.null(fl) && length(fl) == nrow(df)) df <- df[!fl, , drop = FALSE]
  }
  spread_df <- sd_by_tid(df, min_n = min_n)

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

# Apply current metadata filter inputs to a raw stored df (no bypass logic -
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
  # Outlier localizations - Stage 1b, and it MUST come before the length filter
  # below, exactly as in the interactive pipeline: a trimmed point is not a
  # localization of that trace, so it must not count toward its length either.
  if (isTRUE(input$trim_outliers) && "tid" %in% names(df) && nrow(df)) {
    fl <- tryCatch(trace_point_outliers(
            df, k = input$trim_k %||% 4,
            keep_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT),
          error = function(e) NULL)
    if (!is.null(fl) && length(fl) == nrow(df)) df <- df[!fl, , drop = FALSE]
  }
  # TID length filter
  tlr <- input$tid_length_range
  if ("tid" %in% names(df) && !is.null(tlr) && length(tlr) == 2 && all(is.finite(tlr))) {
    tl   <- table(df$tid)
    keep <- names(tl)[tl >= tlr[1] & tl <= tlr[2]]
    df   <- dplyr::filter(df, tid %in% keep)
  }
  # Trace spread (σ) filter - mirrors remove_unstable_traces() from Display tab
  if (isTRUE(input$view_SD %||% FALSE) && "tid" %in% names(df)) {
    tx <- input$target_sd_x; ty <- input$target_sd_y
    if (!is.null(tx) && !is.null(ty) && is.finite(tx) && is.finite(ty) && tx > 0 && ty > 0) {
      min_n <- max(2L, as.integer(input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT))
      # Same estimator as the interactive Stage-2 filter, on the same trimmed
      # point set, or batch and interactive report different trace counts from
      # one setting. Trimming already happened above.
      spread <- sd_by_tid(df, min_n = min_n)
      # Traces with undefined spread (NA, too short) are kept, not filtered.
      valid <- spread$tid[dplyr::coalesce(spread$sd_x, 0) <= tx &
                          dplyr::coalesce(spread$sd_y, 0) <= ty]
      df <- dplyr::filter(df, tid %in% valid)
    }
  }
  df
}

# Compute k-NN distances between unit centroids.
# unit:         "tid"           - one centroid per TID (global)
#               "pcid"          - one centroid per de-blinked molecule (global)
#               "cluster2_tid"  - TID centroids within each cluster2
#               "cluster2_pcid" - de-blinked molecule centroids within each cluster2
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

remove_rogue_traces <- function(df, rule = c("factor_median","mad_rule","axis_sd"),
                                factor = 2.5, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT,
                                robust = TRUE) {
  rule <- match.arg(rule)

  # The per-trace spread MUST use the same estimator as the Stage-2 filter.
  # These rules were robust ACROSS traces (median/MAD of the per-trace values)
  # but the per-trace value itself was a plain sd, so a trace with 2 stray
  # localizations that the spread filter deliberately kept would be deleted here
  # instead - the two stages would silently disagree about what a rogue is.
  spread <- sd_by_tid(df, min_n = 2L, robust = robust)
  if (!nrow(spread)) return(list(df = df, removed_tids = character(0),
                                 removed_n = 0L, note = "no traces"))
  names(spread)[names(spread) == "n_locs"] <- "n"

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
# --- Robust center guess (geometric median -> mean) ---
# Starting centre for the MLE fits, via build_center_args(). Was the convex-hull
# centroid, which is set by the extreme points (breakdown 0); geometric_median()
# has breakdown 0.5, so a loosely-attached member no longer drags the start.
# The MLE optimises the centre from here, so this changes where the search
# begins, not what it converges to - except where the old start was pulled far
# enough out to find a different local optimum.
robust_center <- function(X) {
  tryCatch(geometric_median(X), error = function(e) colMeans(X))
}
# =============================================================================
# SHAPE FITTING - fit a geometric primitive to a cluster of localizations
# =============================================================================
#
# GOAL OF THIS WHOLE BLOCK
#   Given the points belonging to one cluster, ask "what SHAPE is this?" and fit
#   the best sphere / ellipsoid / helix / ring / arc to it. That turns a cloud of
#   localizations into interpretable numbers: a radius, a pitch, an axis length.
#
# THE FITTING PRINCIPLE - ORTHOGONAL-DISTANCE (total least squares)
#   For a shape S with parameters θ, minimise the sum of squared PERPENDICULAR
#   distances from each point to the SURFACE:
#
#        θ* = argmin_θ  Σ_i  || x_i − P_S(x_i; θ) ||²
#
#   where P_S(x) is the nearest point on the shape. This is the right criterion
#   here because localization error is isotropic in space - there is no
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
# MODEL SELECTION - which shape wins?
#   RMSD always favours the most flexible shape (more parameters always fit
#   better), so ranking uses AICc/BIC, which penalise parameter count. AICc is
#   preferred for the small point counts typical of one cluster, and is left
#   blank when a cluster has too few points to justify that model (n <= params+2).
#
# HONEST CAVEAT (mirrored by the EXPLORATORY badge in the UI)
#   The shape dictionary is IMPOSED by this tool - a cluster is scored only
#   against the five shapes offered, so "best" means best OF THESE, never proof
#   of biological geometry. With few, blinking-inflated localizations several
#   shapes are often statistically indistinguishable. Treat the ranking as a
#   hypothesis, and check the residuals rather than the label.
# =============================================================================

# Triangle index list for a parametric (nu x nv) surface grid - the connectivity
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
  # while keeping all residuals small - a geometrically meaningless solution.
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
# nearest point on the ellipsoid surface - identical units (nm²) to the sphere
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
      sum(w * pr$dist^2) / n   # mean squared distance - scale-independent of n
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
    val <- sum(w * pr$dist^2) / n   # mean squared distance - scale-independent of n
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
      val <- sum(w * res2) / n   # mean squared distance - scale-independent of n
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
  # Return [u | v] - two orthonormal vectors spanning the plane ⊥ to normal.
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
  # Uses the arc projection with tc0/ah0 fixed - only orientation is swept.
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
# and return a mesh3d-compatible list (x, y, z, i, j, k - 0-indexed for plotly).
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

  # Text labels at midpoints - one trace for all labels
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

# =============================================================================
# geometric_median() - ROBUST MULTIVARIATE CENTRE (WEISZFELD)
# =============================================================================
#
# GOAL
#   Return the point m minimising the SUM OF EUCLIDEAN DISTANCES to every row of
#   `pts` (the L1 / spatial median). This is the centre used by the heuristic
#   shape fits and as the MLE starting centre.
#
# WHY NOT THE CONVEX-HULL CENTROID (what these fits used before)
#   compute_convex_hull_centroid() is an area-weighted centroid of the hull
#   SURFACE, so it is defined entirely by the extreme points: one stray member
#   drags it, i.e. its breakdown point is 0. The geometric median has breakdown
#   0.5 - half the members must be corrupted before it can be moved without
#   bound - which is what we want on MINFLUX clusters that may still carry a
#   loosely-attached trace or two.
#
# WHY NOT A COORDINATE-WISE median(x), median(y), median(z)
#   That is not rotation-equivariant: rotate the structure and the centre moves
#   relative to it. Since these centres are handed to prcomp() to define ring
#   normals and ellipsoid axes, the centre must rotate WITH the data. The
#   geometric median does; the per-axis median does not.
#
# WHAT THIS IS NOT
#   It is a centre OF THE POINTS, not a centre of curvature. When a structure is
#   only partially labelled (a cap of a vesicle, an arc of a ring) no
#   centre-of-the-points statistic recovers the sphere/ring centre - the
#   geometric median then sits inside the observed patch. Use the MLE pathway
#   (mle_fit_*, which optimises the centre itself) when the centre is the
#   measurement; this is a robust starting point and a heuristic answer.
#
# ALGORITHM
#   Weiszfeld (1937) fixed-point iteration
#       m <- SUM_i (p_i / d_i) / SUM_i (1 / d_i),    d_i = ||p_i - m||
#   started from the coordinate-wise median (already robust, so it lands close
#   and converges in a handful of steps). The classic degeneracy is d_i = 0 when
#   an iterate comes to rest exactly on a data point; d_i is floored at `eps_d`
#   so the weight stays finite.
#
#   Deterministic and exact for the small member counts here (tens to hundreds),
#   which is why this is used rather than Gmedian::Gmedian - that is a stochastic
#   averaged-gradient estimator built for large n.
geometric_median <- function(pts, max_iter = 200L, tol = 1e-7, eps_d = 1e-9) {
  if (is.null(pts)) return(NULL)
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]   # a partial row has no position
  n <- nrow(pts)
  if (n == 0L) return(rep(NA_real_, ncol(pts)))
  if (n <= 2L) return(colMeans(pts))       # 1 point is its own median; for 2 the
                                           # whole segment is optimal - take the midpoint
  m <- apply(pts, 2, stats::median)        # robust starting point
  for (i in seq_len(max_iter)) {
    d <- sqrt(rowSums(sweep(pts, 2, m)^2))
    d[d < eps_d] <- eps_d                  # iterate sitting on a data point
    w <- 1 / d
    m_new <- colSums(pts * w) / sum(w)
    if (max(abs(m_new - m)) < tol) { m <- m_new; break }
    m <- m_new
  }
  unname(m)
}

fit_sphere_geometry <- function(coords_matrix) {
  center <- tryCatch(geometric_median(coords_matrix), error = function(e) colMeans(coords_matrix))
  r      <- stats::median(sqrt(rowSums(sweep(coords_matrix, 2, center)^2)))
  list(center = center, radius = r)
}

fit_ellipsoid_geometry <- function(coords_matrix) {
  center <- tryCatch(geometric_median(coords_matrix), error = function(e) colMeans(coords_matrix))
  Xc     <- sweep(coords_matrix, 2, center)   # center before PCA so axes are aligned to the ellipsoid center
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
  center <- tryCatch(geometric_median(coords_matrix), error = function(e) colMeans(coords_matrix))
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
  center <- tryCatch(geometric_median(coords_matrix),
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
# alpha_surface_3d() - WRAP ONE GROUP OF LOCALIZATIONS IN AN ALPHA SHAPE,
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
#   triangle list in `$triang[, 1:3]` is therefore identical for every alpha -
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
#   the singular count explodes - singular triangles are dangling slivers and
#   drawing them is what produced the old opaque blob.)
#
# HOW THE RADIUS IS PICKED
#   A ladder of candidate radii is built from the group's OWN spacing (multiples
#   of its median nearest-neighbour distance) and classified in a SINGLE
#   ashape3d call - extra alphas cost one column each, not another
#   tetrahedralization. We then take the SMALLEST radius that passes two tests:
#     (1) every point sits in ONE connected component (nothing shed), and
#     (2) the regular triangles form a CLOSED surface - each edge used by
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
# is built here - the caller hulls it. (An `ashape` element used to be returned
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
    # ONE curvature convention app-wide: MEAN curvature H = 1/R (nm^-1), the
    # quantity membrane-bending work is expressed in. This used to be 1/R^2 here
    # (Gaussian curvature K for a sphere) while estimate_cluster2_curvature()
    # reported 1/R, so the same word meant two different things in one codebase.
    # For a fitted sphere both are just reparametrisations of `radius` and carry
    # no information beyond it.
    mean_curvature <- if (radius > 0) 1 / radius else NA
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
    #
    # ecc3 == 0 IS A RESULT, NOT A FAILURE: it is zero exactly when the two
    # SMALLER axes are equal, i.e. the structure is rotationally symmetric about
    # its long axis (rod / cigar / spherocylinder). Verified against built
    # shapes: prolate 400:150:150 -> ecc1 = ecc2 = 0.930, ecc3 = 0.057;
    # triaxial 400:250:100 -> ecc3 = 0.915 against a truth of 0.917.
    #
    # CAVEAT (see docs I-A.6b): PCA orders axes by SAMPLE variance, so the
    # largest sampled direction always becomes PC1 - the ordering manufactures
    # apparent elongation. On a TRUE SPHERE the median reported ecc2 is 0.81 at
    # 10 members, 0.53 at 40, 0.34 at 200, 0.18 at 2000. These numbers are only
    # evidence of elongation when they clearly exceed that floor for the
    # cluster's own member count.
    ecc1 <- if (is.finite(d) && d > 1e-6) sqrt(max(0, 1 - (e/d)^2)) else NA_real_
    ecc2 <- if (is.finite(d) && d > 1e-6) sqrt(max(0, 1 - (f/d)^2)) else NA_real_
    ecc3 <- if (is.finite(e) && e > 1e-6) sqrt(max(0, 1 - (f/e)^2)) else NA_real_
    # Bending (perceived) radius - ASSUMPTION: the cluster represents a curved arc
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
      radius = radius, mean_curvature_inv_nm = mean_curvature,
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
        ", Mean curv=", ifelse(is.na(mean_curvature), "NA", paste0(signif(mean_curvature, 3), " nm^-1"))
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

# =============================================================================
# cluster2_member_col() - WHICH UNIT IS A "MEMBER" OF A MEDIAN CLUSTER
# =============================================================================
#
# GOAL
#   One rule, shared by every per-structure statistic: a member is a DE-BLINKED
#   MOLECULE when de-blinking has run, otherwise a TRACE ID.
#
# WHY IT MATTERS
#   A trace ID is one blinking event. One fluorophore commonly emits several, so
#   counting trace IDs counts blinks, not molecules, and a structure's measured
#   composition then moves with photophysics rather than with biology. This
#   mirrors the server-side smlm_nn_unit() used by the Crowding and
#   NN-vs-size panels, so all per-structure numbers agree on the unit.
cluster2_member_col <- function(df, member_col = NULL) {
  if (!is.null(member_col)) return(member_col)
  if ("cluster" %in% names(df) && any(!is.na(df$cluster) & df$cluster > 0L)) "cluster" else "tid"
}

# =============================================================================
# smlm_loc_precision_nm() - PER-LOCALIZATION PRECISION OF THE CURRENT DATASET
# =============================================================================
#
# GOAL
#   Give the isosurface a PHYSICAL length scale to hang its geometry on, so the
#   drawn size of a structure stops depending on how well it happened to be
#   labelled.
#
# HOW
#   Each trace is one fluorophore sitting still, so the scatter of its own
#   localizations around their centre IS the localization precision. Take the
#   per-trace SD on each axis and report the median trace, matching the
#   convention already used by estimate_gt_positions() (sigma_loc_nm) and by the
#   de-blinking radius.
#
#   The MEDIAN over traces, not the mean: a trace that caught two molecules in
#   sequence has a spread of tens of nm and would drag a mean far off the real
#   single-molecule precision.
#
# WHY NOT THE REPORTED CRLB
#   Not every export carries one, and where it does it is a lower bound on a
#   single localization rather than a measured scatter. The measured spread
#   already contains whatever residual drift and stage vibration survived, which
#   is what actually blurs the rendered structure.
#
# THE AXES ARE NOT EQUAL, AND NOT IN THE DIRECTION CAMERA SMLM TRAINS YOU TO
# EXPECT. MINFLUX z is typically TIGHTER than x/y, the opposite of astigmatism
# or biplane imaging. Measured on the collaborator cond02 files: sx 6.84 / sy 6.81 /
# sz 4.11 nm on B2-010 and 7.07 / 7.16 / 4.17 nm on B3-008, a z:lateral ratio of
# 0.60 and 0.59. Anything that wraps a structure isotropically therefore
# over-wraps it axially by about two thirds, so `per_axis = TRUE` returns the
# three numbers separately and the isosurface uses them.
#
# RETURNS
#   `per_axis = FALSE`: one number in nm (the RMS over the axes in use).
#   `per_axis = TRUE` : c(x, y, z) in nm; z is the lateral mean on 2-D data.
#   `fallback` (recycled) when there are too few multi-localization traces.
smlm_loc_precision_nm <- function(df, min_n = 3L, fallback = 2, per_axis = FALSE) {
  bail <- function() if (isTRUE(per_axis)) rep(fallback, 3L) else fallback
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(bail())
  if (!all(c("tid", "x", "y") %in% names(df))) return(bail())

  # z only counts when the dataset really is 3-D; a 2-D file carries a constant
  # (or all-zero) z whose SD is 0 and would halve the estimate.
  use_z <- "z" %in% names(df) && any(is.finite(df$z) & df$z != 0)

  per_tid <- tryCatch(
    df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        sx = stats::sd(x, na.rm = TRUE),
        sy = stats::sd(y, na.rm = TRUE),
        sz = if (use_z) stats::sd(z, na.rm = TRUE) else NA_real_,
        n  = dplyr::n(), .groups = "drop"
      ) %>%
      dplyr::filter(n >= min_n, is.finite(sx), is.finite(sy)),
    error = function(e) NULL
  )
  if (is.null(per_tid) || !nrow(per_tid)) return(bail())

  # median trace on each axis independently, so a real anisotropy survives
  ax <- c(
    x = stats::median(per_tid$sx, na.rm = TRUE),
    y = stats::median(per_tid$sy, na.rm = TRUE),
    z = if (use_z) stats::median(per_tid$sz, na.rm = TRUE) else NA_real_
  )
  if (!is.finite(ax[["z"]]) || ax[["z"]] <= 0) ax[["z"]] <- mean(ax[c("x", "y")])
  ax[!is.finite(ax) | ax <= 0] <- fallback

  if (isTRUE(per_axis)) return(ax)

  out <- sqrt(mean(ax^2))                       # RMS over the axes
  if (!is.finite(out) || out <= 0) fallback else out
}

# Bandwidth of the isosurface KDE, as a multiple of the member gap. See the
# `sigma` block inside cluster_density_field() for how this was chosen: it is a
# compromise between a surface that fragments (too low) and one that has smoothed
# the structure away (too high), and the mesh smoothing downstream means less
# bandwidth is needed here than when the voxel grid was rendered directly.
SIGMA_GAP_FACTOR <- 0.85

# =============================================================================
# cluster_density_field() - GAUSSIAN KDE OVER A VOXEL GRID + A DENSITY LEVEL
#                           THAT ENCLOSES A STATED FRACTION OF THE LABELLING
# =============================================================================
#
# GOAL
#   Give the 3-D viewer a smooth envelope to wrap a Median Cluster in, built
#   from a field that actually responds to how dense the labelling is.
#
# WHAT IT REPLACED, AND WHY THAT MATTERED
#   The isosurface used to score each voxel by its distance to the SINGLE
#   NEAREST member: val = exp(-d^2 / 2 sigma^2). That is a distance transform,
#   not a density - measured, a probe 60 nm from a lone isolated centroid and a
#   probe 60 nm from a centroid with 30 neighbours within 200 nm both scored
#   0.9231, identical. Visually that means every straggler grows its own
#   detached bubble of exactly the same size as the real structure's rim, which
#   is both misleading and ugly. Summing a kernel over ALL members instead makes
#   the field proportional to local labelling density, so the body of the
#   structure dominates and isolated members fade out.
#
# THE LEVEL IS CALIBRATED TO A SURFACE OFFSET IN NANOMETRES
#   `offset` = "how far outside the member centroids should the skin sit". The
#   level is read off the field at that distance: voxels whose nearest-member
#   distance falls in a band around `offset` vote, and their MEDIAN value is the
#   contour. The surface then follows the offset while density still modulates
#   it - sparse patches pull in slightly, dense ones push out slightly, which is
#   what makes it look like a skin rather than a union of balls.
#
# THE OFFSET IS PHYSICAL, THE BANDWIDTH IS COSMETIC
#   `offset` (AUTO = 3 x localization precision, floored at half the member gap)
#   sets HOW TIGHT the skin is, in nanometres that mean something on the
#   microscope. `sigma` (AUTO = SIGMA_GAP_FACTOR x gap) sets HOW SMOOTH it is.
#   Both are applied PER AXIS, in units of that axis's own precision. They are
#   independent by construction: the level is re-measured at `offset` using the
#   current sigma, so smoothing more does not inflate the structure.
#
#   This is the property that makes the surface comparable between clusters. The
#   earlier rule keyed both dials to the member spacing, so a poorly labelled
#   copy of a structure drew larger than a well labelled one - on an R = 150 nm
#   shell the realised radius ran from 304 nm at 30 members down to 192 nm at
#   400, all of it labelling artifact. Anchored to precision the same four cases
#   land at 166 / 163 / 173 / 169 nm.
#
#   `sampling_limited` in the returned list says which regime you are in: TRUE
#   means the gap floor set the offset, so the surface is as tight as this
#   labelling density permits and the remaining bulk is missing molecules rather
#   than real structure.
#
#   WHY NOT AN ENCLOSED-DENSITY FRACTION (what this replaced). A Gaussian
#   spreads its mass over ~3 sigma in every direction, so "enclose 50% of the
#   labelling density" is already a thick band: measured on an R=150 shell with
#   40 members, frac = 0.50 - the lowest the slider allowed - put the skin
#   56.4 nm out (max 128.5), and frac = 0.85 put it 108.7 nm out. A mass
#   fraction cannot express "hug the points"; a distance can.
#
#   WHY NOT A FRACTION OF THE PEAK either: it collapses too fast. On the same
#   shell, peak fraction 0.4 -> 681 voxels, 0.6 -> 122, 0.8 -> 16, i.e. the
#   surface shatters into fragments before it gets tight.
#
# THE DIALS
#   `offset` sets the wrap and `sigma` the smoothing, both defaulting as above;
#   `pad` adds to or subtracts from the offset and is the coarser tighter/looser
#   control. Grid extent stays derived, so no knob can silently fight another.
#
# BOX
#   `grid_res` is a COST BUDGET of grid_res^3 voxels, not a per-axis count. The
#   voxels are cubic and the box is sized to members + offset + 2 voxels, which
#   is just enough for marching cubes to close the surface without clipping it.
#   Spending the budget per-axis instead gave anisotropic voxels on elongated
#   structures (28 nm along a rod, 10 nm across it, so it broke up lengthwise)
#   and the old offset + 2*sigma margin left 32-70% of the box empty. The box is
#   NOT a user control: when padding was a separate box setting it both overrode
#   a request for a tight box (forced to 3*sigma) and leaked into the contour.
#
# NOT A CALIBRATED DENSITY ESTIMATE
#   The kernel is unnormalised and sampled on a coarse voxel grid, so `value` is
#   not a density in localizations/nm^3. It is a visualization envelope; use the
#   Ripley/PCF panels for density statistics.
cluster_density_field <- function(pts, grid_res = 22L, pad = 0, sigma = NULL,
                                  offset = NULL, precision = NULL,
                                  close_surface = TRUE, sigma_scale = 1) {
  if (is.null(pts)) return(NULL)
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  if (ncol(pts) < 3L) return(NULL)
  pts <- pts[stats::complete.cases(pts), 1:3, drop = FALSE]
  storage.mode(pts) <- "double"
  n <- nrow(pts)
  if (n < 4L) return(NULL)

  # ── work in PRECISION UNITS, not nanometres ──
  # `precision` may be one number or c(x, y, z). Dividing each axis by its own
  # precision makes the space isotropic in the only currency that matters here:
  # how well a molecule's position is known along that axis. Everything below is
  # then a plain isotropic problem, and the mesh is scaled back to nm on the way
  # out. The skin that comes back is an ellipsoidal shell in real space, tighter
  # on whichever axis is measured better.
  #
  # For MINFLUX that axis is usually z (measured z:lateral 0.60 on collaborator cond02),
  # so an isotropic wrap inflates the structure axially by about two thirds. Note
  # the direction: this is the opposite of astigmatic camera SMLM, and assuming
  # the camera convention here would make the error worse rather than better.
  prec_ax <- suppressWarnings(as.numeric(precision %||% NA_real_))
  if (length(prec_ax) == 1L) prec_ax <- rep(prec_ax, 3L)
  if (length(prec_ax) != 3L) prec_ax <- rep(prec_ax[1], 3L)
  prec_ax[!is.finite(prec_ax) | prec_ax <= 0] <- NA_real_
  if (all(is.na(prec_ax))) prec_ax <- rep(2, 3L)
  prec_ax[is.na(prec_ax)] <- mean(prec_ax, na.rm = TRUE)

  # nm <-> precision-unit conversion for the dials the user types in nm. The
  # LATERAL precision is the reference: it is the axis people quote, and it keeps
  # a manual "30 nm wrap" meaning 30 nm laterally.
  prec_ref <- mean(prec_ax[1:2])
  pts_nm   <- pts
  pts      <- sweep(pts, 2, prec_ax, "/")

  # ── two length scales, and they measure different things ──
  # `spacing`  = the TYPICAL member separation (median NN). Reported for
  #              context; it is deliberately NOT what the geometry keys off.
  # `gap`      = the LARGE member separations (90th percentile NN). A closed
  #              surface has to bridge the big gaps, not the typical one, so
  #              this is the scale that decides whether the skin closes.
  #
  # WHY THE DISTINCTION IS NOT PEDANTRY. For a quasi-regular point set the two
  # coincide, but real labelling is Poisson-like, and there the median NN sits
  # far below the gaps that actually need bridging: measured on random points on
  # a 600 x 120 nm rod, median NN was 19.8 nm while the 90th percentile was
  # 39.7 nm, a factor of 2.0. That factor is exactly what the old
  # `offset = 2 x median spacing` rule was silently paying for. Naming the gap
  # directly lets the offset shrink to something physical without the surface
  # shattering.
  nnd <- tryCatch({
    d <- RANN::nn2(pts, pts, k = 2L)$nn.dists[, 2]
    d[is.finite(d) & d > 0]
  }, error = function(e) numeric(0))
  spacing <- if (length(nnd)) stats::median(nnd) else NA_real_
  gap     <- if (length(nnd)) as.numeric(stats::quantile(nnd, 0.90, na.rm = TRUE)) else NA_real_
  if (!is.finite(spacing) || spacing <= 0 || !is.finite(gap) || gap <= 0) {
    diag_nm <- sqrt(sum((apply(pts, 2, max) - apply(pts, 2, min))^2))
    spacing <- if (is.finite(diag_nm) && diag_nm > 0) diag_nm / 10 else 20
    gap     <- spacing
  }

  # In precision units one unit IS one localization sigma, on every axis.
  precision <- 1

  # ── surface offset: how far outside the members the skin sits ──
  # AUTO = 3 x localization precision, i.e. a 3-sigma skin around where each
  # molecule actually is, floored at half the gap.
  #
  # WHY 3 x PRECISION AND NOT 2 x SPACING (what this replaced). The old rule read
  # its scale off the labelling, so the SAME structure drew bigger when it was
  # sampled worse. Measured on an R = 150 nm shell, the realised outer radius was
  # 304 nm at 30 members (+103%), 263 nm at 60 (+76%), 216 nm at 150 (+44%) and
  # 192 nm at 400 (+28%): the drawn size tracked labelling density, and a sparse
  # shell even rendered hollow where a dense one filled in solid. Anchoring to
  # precision instead put the same four cases at 166 / 163 / 173 / 169 nm, which
  # is the truth (150) plus the shell's own thickness plus the offset.
  #
  # WHY THE GAP FLOOR IS NOT OPTIONAL. A continuous surface passing outside a set
  # of points cannot come nearer than half their separation to all of them, so an
  # offset below gap/2 cannot close whatever the bandwidth is. The floor states
  # that limit instead of letting the surface shatter into one bead per molecule.
  # On well-labelled data the precision term wins and the skin is genuinely tight;
  # on sparse data the floor wins and the offset reports the real sampling limit.
  offset <- suppressWarnings(as.numeric(offset %||% NA_real_)) / prec_ref  # nm in
  if (!is.finite(offset) || offset <= 0) offset <- max(3 * precision, 0.5 * gap)

  # `pad` adds to (or, negative, subtracts from) that - the coarser/tighter dial.
  pad <- suppressWarnings(as.numeric(pad %||% 0)) / prec_ref                # nm in
  if (!is.finite(pad)) pad <- 0
  offset <- max(offset + pad, 0.5 * precision)     # floor: never collapse to nothing

  # ── bandwidth: the SMOOTHNESS dial, largely independent of the offset ──
  # AUTO = SIGMA_GAP_FACTOR x gap. Because the contour level is re-measured at
  # `offset` with the current sigma, raising sigma does not inflate the surface
  # the way it did when the two were tied together. Measured on collaborator cond02, the
  # MEDIAN vertex-to-molecule distance moved 25.6 -> 29.0 nm against a 25.3 nm
  # target across a 5x sigma sweep.
  #
  # It is not a perfect invariance and should not be sold as one: the same sweep
  # took the MAXIMUM distance from 39.8 to 65.7 nm, because a wider kernel
  # absorbs small protrusions into the body rather than following them. Sigma
  # holds the typical wrap and trades away fine relief.
  #
  # Measured at fixed offset, sweeping sigma over 0.75 to 2.0 x gap on shells, a
  # solid ball and rods: components fell to 1 by ~1.25x, shells stayed HOLLOW up
  # to 1.5x and collapsed into a solid ball at 2.0x. The factor now sits BELOW
  # that closure point on purpose, because two things downstream did not exist
  # when it was measured: the level search closes the surface without spending
  # bandwidth on it (see `close_surface`), and the mesh is smoothed as geometry
  # rather than as density. Carrying enough bandwidth to close the surface on its
  # own is what made real structures look oversmoothed - it blurs the measurement
  # to fix what is really a topology problem. Smoothing the mesh instead leaves
  # the density where it was measured.
  #
  # It was previously tied to the offset (max(0.5 * offset, spacing)), so
  # tightening the wrap also sharpened the field and re-fragmented the surface -
  # the two dials fought each other.
  # `sigma_scale` tunes the AUTOMATIC bandwidth only; an explicit `sigma` in nm
  # is a full override and is left alone. Keeping them separate means the
  # relative dial stays meaningful no matter what the data does to the auto value.
  sigma_scale <- suppressWarnings(as.numeric(sigma_scale %||% 1))
  if (!is.finite(sigma_scale) || sigma_scale <= 0) sigma_scale <- 1
  sigma <- suppressWarnings(as.numeric(sigma %||% NA_real_)) / prec_ref      # nm in
  if (!is.finite(sigma) || sigma <= 0)
    sigma <- max(SIGMA_GAP_FACTOR * sigma_scale * gap, precision)

  # ── box: cubic voxels, spent where the surface actually is ──
  # `grid_res` is a COST budget (grid_res^3 voxels), not a per-axis count. Two
  # things were wrong with using it per-axis:
  #
  #   1. margin = offset + 2 * sigma padded the box far past anything the surface
  #      could reach - the skin stops at ~offset from the outermost member - so
  #      32-70% of the voxels held empty space and coarsened the rest. Measured
  #      box occupancy on 300 nm structures: 30% at 30 members, 68% at 400.
  #   2. an equal count per axis gives ANISOTROPIC voxels on anything elongated.
  #      A 600 x 120 nm rod got 28.2 nm voxels along x and 10 nm across it, so it
  #      broke up lengthwise while being oversampled across.
  #
  # Cubic voxels under a fixed budget fix both at identical cost: the same rod
  # now gets 11.7 nm voxels on a 59x20x20 grid, and a 300 nm shell goes from 25.6
  # to 18.3 nm. Two voxels of headroom past the offset keeps marching cubes from
  # clipping the surface against the wall of the box.
  grid_res <- max(8L, as.integer(round(grid_res)))
  budget   <- as.numeric(grid_res)^3

  rng <- apply(pts, 2, range)
  ext <- pmax(rng[2, ] - rng[1, ], 1e-9)           # guard a flat/degenerate axis

  # margin depends on the voxel size and vice versa; one fixed-point pass from
  # margin = offset is enough, the correction is a few percent
  vox_for <- function(budget) {                    # solve vox against its own margin
    v <- (prod(ext + 2 * offset) / budget)^(1 / 3)
    for (.i in 1:2) v <- (prod(ext + 2 * (offset + 2 * v)) / budget)^(1 / 3)
    v
  }
  vox <- vox_for(budget)

  # ── resolution floor: the voxel must be able to RESOLVE the skin ──
  # Marching cubes needs the enclosed shell to be a couple of voxels thick or it
  # drops pieces of it, and the shell is only ~`offset` thick. Now that the
  # offset is precision-anchored it is often 10-13 nm, well under the ~19 nm
  # voxel that grid_res = 22 buys - measured, that alone broke a closed shell
  # into 3 fragments purely as a sampling artifact, with nothing wrong with the
  # field. Refine until the skin is ~1.5 voxels thick.
  #
  # This raises the user's grid_res, so it is capped: VOXEL_CAP bounds the redraw
  # at roughly two seconds for the largest clusters. Beyond that the surface is
  # under-resolved and `resolved` comes back FALSE rather than the app quietly
  # spending ten seconds on one cluster.
  VOXEL_CAP <- 2e5
  vox_max   <- offset / 1.5
  resolved  <- TRUE
  if (vox > vox_max) {
    need <- prod(ext + 2 * (offset + 2 * vox_max)) / (vox_max^3)
    if (need <= VOXEL_CAP) {
      vox <- vox_max
    } else {
      vox      <- vox_for(VOXEL_CAP)
      resolved <- vox <= vox_max
    }
  }
  margin <- offset + 2 * vox

  n_ax <- pmax(6L, as.integer(ceiling((ext + 2 * margin) / vox)) + 1L)
  xs <- seq(rng[1, 1] - margin, rng[2, 1] + margin, length.out = n_ax[1])
  ys <- seq(rng[1, 2] - margin, rng[2, 2] + margin, length.out = n_ax[2])
  zs <- seq(rng[1, 3] - margin, rng[2, 3] + margin, length.out = n_ax[3])
  grid <- as.matrix(expand.grid(x = xs, y = ys, z = zs))

  # ── truncated KDE: a Gaussian is numerically dead past 3 sigma, so only the
  #    nearest members can contribute.
  #
  # k IS MEASURED, NOT ASSUMED. It used to be min(n, 512), which sets the cost of
  # every grid point by the SIZE OF THE CLUSTER rather than by how many members
  # can actually reach a probe. The nn2 result is a (grid points x k) matrix, so
  # at k = 512 a 280k-voxel grid allocates 145M distances (~1.2 GB) to multiply
  # almost all of them by a weight of zero: measured 22.8 s for one 1500-member
  # cluster where the real neighbour count within 3 sigma was a few dozen.
  #
  # Instead, count how many members genuinely fall within 3 sigma of a member and
  # size k from that, with headroom. The probe subset is a deterministic stride,
  # never sample() - drawing at random here would advance the GLOBAL RNG stream
  # and silently perturb every downstream random draw (the permutation nulls, the
  # alpha-surface subsample), and a redraw could jitter the surface.
  k_probe <- as.integer(round(seq(1, n, length.out = min(n, 256L))))
  k_probe <- unique(k_probe)
  k_need <- tryCatch({
    kn <- RANN::nn2(pts, pts[k_probe, , drop = FALSE], k = min(n, 512L))
    max(rowSums(kn$nn.dists <= 3 * sigma))
  }, error = function(e) min(n, 512L))
  # 1.5x headroom covers members denser than anywhere in the probe subset; the
  # floor keeps small clusters exact. Truncation past this is not an
  # approximation being traded away - those weights are already numerically 0.
  kmax <- min(n, max(8L, as.integer(ceiling(1.5 * k_need)) + 4L))
  nn <- tryCatch(RANN::nn2(pts, grid, k = kmax), error = function(e) NULL)
  if (is.null(nn)) return(NULL)
  w <- exp(-(nn$nn.dists^2) / (2 * sigma * sigma))
  w[nn$nn.dists > 3 * sigma] <- 0
  val <- as.numeric(rowSums(w))
  if (!any(is.finite(val)) || max(val, na.rm = TRUE) <= 0) return(NULL)

  # ── level calibrated to the offset, measured OFF-GRID ──
  # Probe points are scattered on spheres of radius `offset` around randomly
  # chosen members; the median field value over the probes that really do sit
  # ~offset from the point set is the contour.
  #
  # This is deliberately independent of the voxel grid. Reading the level from
  # grid voxels instead lets the BOX size leak into the surface: one far-flung
  # member stretches the bounding box, the voxels coarsen, the "within one
  # voxel of the offset" band widens to tens of nm, and the median is taken
  # over voxels far from any member - which dropped the level low enough to
  # re-draw the isolated straggler that the density field exists to suppress.
  # Directions come from a Fibonacci sphere and members are cycled, so the probe
  # set is deterministic (a redraw cannot jitter the surface) and evenly spread.
  # Drawing them at random would also have meant calling set.seed() here, which
  # would reset the GLOBAL RNG stream and silently perturb every downstream
  # random draw - the permutation nulls, the alpha-surface subsample.
  n_probe <- 1024L
  kk    <- seq_len(n_probe) - 0.5
  phi   <- acos(1 - 2 * kk / n_probe)
  theta <- pi * (1 + sqrt(5)) * kk
  dir   <- cbind(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi))
  mi    <- rep_len(seq_len(n), n_probe)
  probes <- pts[mi, , drop = FALSE] + dir * offset

  pn <- tryCatch(RANN::nn2(pts, probes, k = kmax), error = function(e) NULL)
  if (is.null(pn)) return(NULL)
  # keep probes that are genuinely at the offset - the rest landed back inside
  # the structure, where the field says nothing about where the skin belongs
  keep_p <- pn$nn.dists[, 1] >= 0.8 * offset
  if (!any(keep_p)) keep_p <- rep(TRUE, nrow(probes))
  pw <- exp(-(pn$nn.dists^2) / (2 * sigma * sigma))
  pw[pn$nn.dists > 3 * sigma] <- 0
  level <- stats::median(rowSums(pw)[keep_p])

  # plotly draws nothing if the level sits at or above the peak
  vmax <- max(val, na.rm = TRUE)
  if (!is.finite(level) || level <= 0) level <- vmax * 1e-3
  if (level >= vmax) level <- vmax * 0.999

  # ── solve the level down until the surface is ONE body ──
  # The offset fixes where the skin belongs; this asks the far cheaper question
  # of whether that level actually yields a connected surface, and if not, how
  # far it has to give. Lowering the level grows the surface, which is the exact
  # inflation the offset rewrite removed, so the search is BOUNDED and the cost
  # is reported rather than absorbed: `closed` says whether it succeeded and
  # `level_relaxed` how many times the level had to step down.
  #
  # The field is computed once; only the threshold moves. Connectivity is tested
  # on the voxel mask by label propagation, which is a few vectorised array
  # shifts over a grid this size.
  lvl_components <- function(thr) {
    mask <- array(val >= thr, n_ax)
    if (!any(mask)) return(.Machine$integer.max)
    lab <- array(0L, n_ax); idx <- which(mask); lab[idx] <- idx
    repeat {
      prev <- lab
      for (ax in 1:3) for (sh in c(-1L, 1L)) {
        sl <- switch(ax,
          lab[if (sh == 1L) c(1L, seq_len(n_ax[1] - 1L)) else c(2:n_ax[1], n_ax[1]), , , drop = FALSE],
          lab[, if (sh == 1L) c(1L, seq_len(n_ax[2] - 1L)) else c(2:n_ax[2], n_ax[2]), , drop = FALSE],
          lab[, , if (sh == 1L) c(1L, seq_len(n_ax[3] - 1L)) else c(2:n_ax[3], n_ax[3]), drop = FALSE])
        hit <- mask & sl > 0L & (lab == 0L | sl < lab)
        lab[hit] <- sl[hit]
      }
      lab[!mask] <- 0L
      if (identical(lab, prev)) break
    }
    length(unique(lab[mask]))
  }

  # NA, not TRUE: when the search is off nothing has checked connectivity, and
  # reporting an untested surface as closed is worse than reporting nothing
  closed <- NA; level_relaxed <- 0L
  if (isTRUE(close_surface)) {
    ncomp <- lvl_components(level)
    # 0.85 per step, at most 6 steps: the level can fall to ~0.38 of the
    # offset-calibrated value and no further. Past that the surface has stopped
    # describing the offset at all and a fragmented honest answer beats a
    # connected invented one.
    while (ncomp > 1L && level_relaxed < 6L) {
      level <- level * 0.85
      level_relaxed <- level_relaxed + 1L
      ncomp <- lvl_components(level)
    }
    closed <- ncomp == 1L
  }

  # ── back to nanometres ──
  # Everything above ran in precision units. `grid` is scaled per axis (so the
  # caller gets real coordinates); the scalar diagnostics are reported against
  # the lateral precision, which is the axis they are quoted in.
  grid_nm <- sweep(grid, 2, prec_ax, "*")

  list(grid = grid_nm, value = val, level = level,
       sigma = sigma * prec_ref, offset = offset * prec_ref,
       spacing = spacing * prec_ref, gap = gap * prec_ref,
       precision = prec_ax, precision_ref = prec_ref,
       margin = margin * prec_ref, voxel = vox * prec_ref,
       # axis scale factors, so a caller can rebuild the anisotropic geometry
       axis_scale = prec_ax, dims = n_ax, n_members = n,
       # grid/level in the isotropic working space, for the mesh builder
       grid_units = grid, breaks_units = list(x = xs, y = ys, z = zs),
       # FALSE when the voxel cap stopped the grid short of resolving the skin;
       # the surface is then under-sampled and may show gaps that are not real
       resolved = resolved,
       # did the level search reach a single connected body, and at what cost
       closed = closed, level_relaxed = level_relaxed,
       # TRUE when the sampling, not the microscope, set the wrap distance -
       # i.e. the surface is as tight as this labelling density allows, no tighter
       sampling_limited = (0.5 * gap) > (3 * precision))
}


# =============================================================================
# density_field_mesh() - MARCHING CUBES OVER A DENSITY FIELD, WELDED + SMOOTHED
# =============================================================================
#
# GOAL
#   Turn the voxel field from cluster_density_field() into a genuine triangle
#   mesh that renders as a smooth skin, instead of handing the raw voxels to
#   plotly's `isosurface` trace.
#
# WHY NOT plotly's OWN isosurface TRACE (what this replaced)
#   plotly runs marching cubes per cell and shades the result FLAT, so every
#   voxel face is visible as its own facet. That is the blocky look: at the
#   default grid the voxels were tens of nm across, so the structure was drawn
#   with tens-of-nm flat plates regardless of how good the underlying data was.
#   It also had to be driven with a hairline isomin/isomax band around the level,
#   which renders a thin SLAB with an inner and an outer wall - under a
#   translucent opacity those two walls shade separately and read as bulk.
#
# THE THREE STEPS THAT MAKE IT LOOK LIKE A SURFACE
#   1. MARCHING CUBES (misc3d::computeContour3d) on the field. This is the same
#      algorithm plotly uses; the point is what happens next.
#   2. WELD. computeContour3d emits independent triangles, so a shared corner
#      appears once per triangle touching it - measured 83% duplicate vertex rows
#      on a test sphere. plotly can only average normals across triangles that
#      literally share a vertex INDEX, so without welding `flatshading = FALSE`
#      has nothing to work with and the mesh shades flat no matter what. Welding
#      is what buys smooth shading; it also cuts the payload about six-fold.
#   3. TAUBIN SMOOTH. Marching cubes leaves a staircase at the voxel scale.
#      Taubin's lambda/mu pass (a positive smoothing step followed by a slightly
#      larger negative one) removes it WITHOUT the shrinkage that plain Laplacian
#      smoothing causes - a Laplacian-only mesh creeps inward every iteration and
#      would quietly undo the offset calibration the field just established.
#
# WHY SMOOTH THE MESH RATHER THAN THE DENSITY
#   Both remove the staircase, but they are not the same operation. Raising the
#   KDE bandwidth blurs the MEASUREMENT: it merges neighbouring molecules and
#   pulls real detail out of the structure. Smoothing the mesh only regularises
#   the ISOSURFACE EXTRACTION, which is a rendering artifact to begin with. That
#   split is why SIGMA_GAP_FACTOR could be lowered from 1.25 to 0.85 here.
#
# NOT A VOLUME ESTIMATOR
#   The mesh is for looking at. Use cluster_alpha_volume() for enclosed volume;
#   the smoothing step alone moves this surface by a few percent.
#
# RETURNS
#   list(x, y, z, i, j, k) with 0-based i/j/k ready for plotly::add_trace(type =
#   "mesh3d"), plus n_vertices / n_triangles, or NULL when the level yields no
#   surface at all.
density_field_mesh <- function(fld, smooth_iter = 8L, lambda = 0.5, mu = -0.53,
                               weld_tol = 1e-6) {
  if (is.null(fld) || is.null(fld$value) || is.null(fld$dims)) return(NULL)
  if (!requireNamespace("misc3d", quietly = TRUE)) return(NULL)

  n_ax <- fld$dims
  vol  <- array(fld$value, n_ax)
  br   <- fld$breaks_units          # isotropic working space (precision units)
  if (is.null(br)) return(NULL)

  # ── 1. marching cubes ──
  tri <- tryCatch(
    misc3d::computeContour3d(vol, level = fld$level, x = br$x, y = br$y, z = br$z),
    error = function(e) NULL
  )
  if (is.null(tri) || !is.matrix(tri) || nrow(tri) < 3L) return(NULL)

  # ── 2. weld ──
  # Quantise to a tolerance relative to the voxel, so vertices that marching
  # cubes computed independently but that ARE the same point collapse together.
  # match() on the packed key is the vectorised way to do this; paste() of three
  # rounded numbers is stable because they came from the same interpolation.
  vox_u <- if (length(br$x) > 1L) abs(br$x[2] - br$x[1]) else 1
  q     <- round(tri / (vox_u * weld_tol))
  key   <- paste(q[, 1], q[, 2], q[, 3], sep = ",")
  first <- !duplicated(key)
  verts <- tri[first, , drop = FALSE]
  idx   <- match(key, key[first])            # every row -> its welded vertex

  faces <- matrix(idx, ncol = 3L, byrow = TRUE)
  # marching cubes can emit slivers where two corners welded together; a triangle
  # with a repeated vertex has no area and no normal, and would poison the
  # averaged normal of every vertex it touches
  keep_f <- faces[, 1] != faces[, 2] & faces[, 2] != faces[, 3] & faces[, 1] != faces[, 3]
  faces  <- faces[keep_f, , drop = FALSE]
  if (!nrow(faces)) return(NULL)

  # ── 3. Taubin smoothing ──
  if (smooth_iter > 0L && nrow(verts) > 3L) {
    # Undirected edge list, each edge once per direction so the neighbour average
    # is symmetric. Built once and reused for every iteration.
    e_from <- c(faces[, 1], faces[, 2], faces[, 3], faces[, 2], faces[, 3], faces[, 1])
    e_to   <- c(faces[, 2], faces[, 3], faces[, 1], faces[, 1], faces[, 2], faces[, 3])
    ord    <- !duplicated(paste(e_from, e_to, sep = ","))
    e_from <- e_from[ord]; e_to <- e_to[ord]
    deg    <- tabulate(e_from, nbins = nrow(verts))
    deg[deg == 0L] <- 1L                     # isolated vertex: leave it alone

    # Neighbour-averaging operator as a sparse matrix, built once. Each
    # iteration is then one sparse matvec instead of three tapply() calls over
    # tens of thousands of edges, which is what made this the slow step.
    A <- Matrix::sparseMatrix(i = e_from, j = e_to, x = 1 / deg[e_from],
                              dims = c(nrow(verts), nrow(verts)))
    for (it in seq_len(smooth_iter)) {
      verts <- verts + lambda * (as.matrix(A %*% verts) - verts)   # shrink
      verts <- verts + mu     * (as.matrix(A %*% verts) - verts)   # push back out
    }
  }

  # ── back to nanometres: undo the per-axis precision scaling ──
  sc <- fld$axis_scale %||% c(1, 1, 1)
  verts <- sweep(verts, 2, sc, "*")

  list(x = verts[, 1], y = verts[, 2], z = verts[, 3],
       i = faces[, 1] - 1L, j = faces[, 2] - 1L, k = faces[, 3] - 1L,   # plotly is 0-based
       n_vertices = nrow(verts), n_triangles = nrow(faces))
}


# GOAL: one hex colour for a set of members, or NA when they carry none.
# Used by the isosurface so a surface is coloured by the SAME channel_color the
# localizations inside it use. Returns NA (not a guess) when the members span
# more than one channel colour, because a two-channel structure has no single
# honest colour and the caller should fall back to its own palette.
iso_channel_color <- function(members) {
  if (!is.data.frame(members) || !nrow(members)) return(NA_character_)
  if (!"channel_color" %in% names(members)) return(NA_character_)
  cc <- unique(members$channel_color[!is.na(members$channel_color)])
  if (length(cc) == 1L) as.character(cc) else NA_character_
}
# =============================================================================
# add_cluster_isosurface() - ONE ISOSURFACE TRACE, USED BY EVERY VIEWING MODE
# =============================================================================
#
# GOAL
#   Wrap one structure in one smooth skin and add it to a plotly 3-D figure.
#   Normal mode, the 2-colour split and multiplex channels all call this; they
#   differ only in which points they hand over and which colour they ask for.
#
# WHY IT EXISTS
#   The cluster2 and tc2 render blocks were near-duplicates that had already
#   drifted apart: only one of them had been switched to de-blinked molecule
#   centroids, and only one honoured channel colours. A wrapping rule that lives
#   in two places is a wrapping rule that will disagree with itself.
#
# COLOUR IS THE CALLER'S DECISION
#   This takes a single hex colour and does not know what it means. The caller
#   picks it from the cluster palette, the channel table or a fixed user colour,
#   which is what keeps "colour by channel" from having to be re-implemented per
#   viewing mode.
#
# RENDERING NOTES
#   The trace is a `mesh3d`, not plotly's `isosurface`. flatshading = FALSE plus
#   the welded vertices from density_field_mesh() is what produces smooth shading;
#   either alone gives flat facets. Lighting is deliberately low-specular: a
#   translucent skin with a hot highlight reads as a plastic bead rather than as
#   a density envelope.
#
# RETURNS
#   The figure, with one trace added. On any failure (too few members, no
#   surface at that level) the figure comes back UNCHANGED rather than erroring,
#   because one bad cluster must not blank the whole plot.
add_cluster_isosurface <- function(fig, pts, color = "#00CED1", opacity = 0.25,
                                   name = NULL, showlegend = FALSE,
                                   grid_res = 22L, pad = 0, sigma = NULL,
                                   offset = NULL, precision = NULL,
                                   smooth_iter = 8L, close_surface = TRUE,
                                   sigma_scale = 1, diag_env = NULL, cache = NULL) {
  # `diag_env` collects one row per surface so the UI can report what the
  # automatic geometry actually chose. An environment rather than a return value
  # because the caller wants the figure back in a pipe; nothing here reads it.
  cache_row <- NULL
  note <- function(...) {
    row <- data.frame(..., stringsAsFactors = FALSE)
    cache_row <<- row
    if (is.null(diag_env)) return(invisible(NULL))
    diag_env$rows <- c(diag_env$rows, list(row))
  }
  if (is.null(pts)) return(fig)
  pts <- as.matrix(pts)
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]
  if (nrow(pts) < 4L || ncol(pts) < 3L) {
    note(name = name %||% NA_character_, n_members = nrow(pts), wrap = NA_real_,
         smoothing = NA_real_, voxel = NA_real_, sampling_limited = NA,
         closed = NA, resolved = NA, drawn = FALSE)
    return(fig)
  }

  # ── geometry cache ──
  # Colour and opacity do not enter the key, because they do not change the mesh.
  # A Shiny render re-runs on ANY of its reactive dependencies, so without this a
  # change of opacity or colour mode rebuilt every field and every mesh from
  # scratch. The caller clears the cache when the committed geometry changes.
  ck <- NULL
  if (!is.null(cache)) {
    ck <- digest::digest(list(pts, grid_res, pad, sigma, offset, precision,
                              smooth_iter, close_surface, sigma_scale))
    hit <- cache[[ck]]
    if (!is.null(hit)) {
      if (!is.null(diag_env) && !is.null(hit$diag))
        diag_env$rows <- c(diag_env$rows, list(hit$diag))
      if (is.null(hit$msh)) return(fig)
      return(iso_mesh_trace(fig, hit$msh, color, opacity, name, showlegend))
    }
  }

  fld <- tryCatch(
    cluster_density_field(pts, grid_res = grid_res, pad = pad, sigma = sigma,
                          offset = offset, precision = precision,
                          close_surface = close_surface, sigma_scale = sigma_scale),
    error = function(e) NULL
  )
  if (is.null(fld)) {
    note(name = name %||% NA_character_, n_members = nrow(pts), wrap = NA_real_,
         smoothing = NA_real_, voxel = NA_real_, sampling_limited = NA,
         closed = NA, resolved = NA, drawn = FALSE)
    return(fig)
  }

  msh <- tryCatch(density_field_mesh(fld, smooth_iter = smooth_iter),
                  error = function(e) NULL)
  note(name = name %||% NA_character_, n_members = nrow(pts),
       wrap = fld$offset, smoothing = fld$sigma, voxel = fld$voxel,
       sampling_limited = isTRUE(fld$sampling_limited),
       closed = fld$closed, resolved = isTRUE(fld$resolved),
       drawn = !is.null(msh))
  if (!is.null(ck)) cache[[ck]] <- list(msh = msh, diag = cache_row)
  if (is.null(msh)) return(fig)

  iso_mesh_trace(fig, msh, color, opacity, name, showlegend)
}

# GOAL: add one finished mesh to a figure. Split out so a cache hit and a fresh
# build draw through exactly the same code.
#
# COLOUR MUST NOT BE PASSED AS `color`. plotly's R wrapper intercepts a `color`
# argument and routes it through its own scale-mapping machinery, which has no
# output slot for mesh3d - so the attribute is silently DROPPED and every surface
# renders in plotly's default. Verified against plotly_build(): both
# `color = "#FF0000"` and `color = I("#FF0000")` vanish from the built trace,
# while intensity/colorscale, facecolor and vertexcolor all survive.
#
# A flat intensity against a two-stop colorscale is the lightest of the three:
# one number per vertex rather than a hex STRING per vertex (vertexcolor) or per
# triangle (facecolor), which matters when a view draws tens of surfaces.
# cmin/cmax are pinned because an all-equal intensity would otherwise autoscale
# to a degenerate cmin == cmax range.
iso_mesh_trace <- function(fig, msh, color, opacity, name, showlegend) {
  nv <- length(msh$x)
  fig %>% plotly::add_trace(
    type = "mesh3d",
    x = msh$x, y = msh$y, z = msh$z,
    i = msh$i, j = msh$j, k = msh$k,
    intensity = rep(0, nv),
    colorscale = list(list(0, color), list(1, color)),
    cmin = 0, cmax = 1, showscale = FALSE,
    opacity = opacity,
    flatshading = FALSE,
    lighting = list(ambient = 0.55, diffuse = 0.8, specular = 0.12,
                    roughness = 0.9, fresnel = 0.1),
    lightposition = list(x = 1e4, y = 1e4, z = 1e4),
    name = name %||% "Isosurface",
    showlegend = isTRUE(showlegend),
    hoverinfo = "skip",
    inherit = FALSE
  )
}

# =============================================================================
# iso_status_ui() - WHAT THE AUTOMATIC ISOSURFACE GEOMETRY ACTUALLY CHOSE
# =============================================================================
#
# GOAL
#   Report, under the isosurface controls, the numbers the surface was built
#   with and anything that stopped it being drawn properly.
#
# WHY THIS IS NOT OPTIONAL POLISH
#   The wrap distance, the smoothing and the grid are all derived from the data,
#   so "Tighten / loosen" adjusts a value the user cannot otherwise see. Without
#   this readout, -20 nm is a blind edit. It also separates the two ways a
#   surface can look wrong that would otherwise be indistinguishable on screen:
#   a structure that is genuinely fragmented, and one that is merely
#   under-resolved by the voxel grid.
#
# THE TWO REGIMES
#   The wrap is 3 x localization precision unless that falls below half the
#   molecule spacing, in which case spacing wins. Which one bound it is the
#   single most useful thing to know about a surface: a spacing-bound surface is
#   as tight as this labelling density permits and its remaining bulk is missing
#   molecules, not structure.
#
# RETURNS
#   A shiny tag, or NULL when nothing has been drawn yet.
iso_status_ui <- function(rows) {
  if (is.null(rows) || !length(rows)) return(NULL)
  d <- if (is.data.frame(rows)) rows else do.call(rbind, rows)
  if (is.null(d) || !nrow(d)) return(NULL)

  ok  <- d[isTRUE_vec(d$drawn), , drop = FALSE]
  n_f <- sum(!isTRUE_vec(d$drawn))
  if (!nrow(ok)) {
    return(shiny::div(class = "iso-status",
      shiny::tags$span(style = "color:#b45309;",
        sprintf("No surface drawn (%d cluster%s too small or too sparse).",
                n_f, if (n_f == 1L) "" else "s"))))
  }

  rng <- function(v, digits = 0) {
    v <- v[is.finite(v)]
    if (!length(v)) return("n/a")
    if (diff(range(v)) < 0.5) sprintf("%.*f nm", digits, v[1])
    else sprintf("%.*f to %.*f nm", digits, min(v), digits, max(v))
  }

  n_samp <- sum(isTRUE_vec(ok$sampling_limited))
  n_prec <- nrow(ok) - n_samp
  regime <- if (n_samp == 0L) "Wrap set by localization precision."
    else if (n_prec == 0L)   "Wrap set by molecule spacing."
    else sprintf("Wrap set by localization precision on %d, by molecule spacing on %d.",
                 n_prec, n_samp)

  warn <- character(0)
  n_open <- sum(!is.na(ok$closed) & !ok$closed)
  if (n_open) warn <- c(warn, sprintf(
    "%d surface%s did not close into one body. The structure is sampled too sparsely to wrap this tightly; loosen the wrap or raise smoothing.",
    n_open, if (n_open == 1L) "" else "s"))
  n_unres <- sum(!isTRUE_vec(ok$resolved))
  if (n_unres) warn <- c(warn, sprintf(
    "%d surface%s hit the voxel cap and %s under-resolved. Gaps in %s may not be real.",
    n_unres, if (n_unres == 1L) "" else "s",
    if (n_unres == 1L) "is" else "are", if (n_unres == 1L) "it" else "them"))
  if (n_f) warn <- c(warn, sprintf(
    "%d cluster%s had too few molecules to wrap.", n_f, if (n_f == 1L) "" else "s"))

  shiny::div(
    class = "iso-status",
    style = paste("margin:6px 0 10px 0; padding:8px 10px; border-radius:4px;",
                  "background:rgba(127,127,127,0.10); font-size:12px; line-height:1.45;"),
    shiny::tags$div(shiny::tags$b(sprintf("%d surface%s drawn.", nrow(ok),
                                          if (nrow(ok) == 1L) "" else "s")),
                    sprintf(" Wrap %s, smoothing %s, voxel %s.",
                            rng(ok$wrap), rng(ok$smoothing), rng(ok$voxel, 1))),
    shiny::tags$div(style = "opacity:0.8;", regime),
    if (length(warn)) shiny::tags$div(
      style = "margin-top:5px; color:#b45309;",
      lapply(warn, function(w) shiny::tags$div(w))
    )
  )
}

# isTRUE() is not vectorised; these columns can carry NA (a diagnostic that was
# never tested), and NA must read as FALSE rather than propagate.
isTRUE_vec <- function(v) !is.na(v) & v
# =============================================================================
# cluster_alpha_volume() - ENCLOSED VOLUME OF ONE STRUCTURE, ALPHA CHOSEN FROM
#                          THAT STRUCTURE'S OWN SPACING
# =============================================================================
#
# GOAL
#   Return list(volume, alpha): the alpha-shape volume enclosed by a structure's
#   member centroids, plus the alpha radius actually used, so the number is
#   reportable rather than dependent on a hidden constant.
#
# WHY THE ALPHA CANNOT BE A CONSTANT
#   Alpha is a probe radius: gaps wider than it are left open, so a value below
#   the members' own spacing returns an empty shape and a volume of ZERO. The
#   previous fixed `alpha = 10` did exactly that - on a 150 nm shell whose
#   members sat 18.5 nm apart it returned 0, as it did for every closed
#   structure tested; that shell needs alpha ~ 200 nm to close. A 30 nm molecule
#   and a 500 nm vesicle cannot share one probe radius, which is the same
#   argument alpha_surface_3d() makes for the rendered surfaces.
#
# HOW ALPHA IS CHOSEN
#   Walk a ladder of radii built from the group's own median nearest-neighbour
#   spacing and take the SMALLEST that both (1) wraps every member in ONE
#   connected component and (2) yields a WATERTIGHT boundary - every edge shared
#   by exactly two faces. Same two tests alpha_surface_3d() uses for rendering.
#
#   Connectivity alone is not enough for SOLID structures. A filled cloud is
#   connected at a very small alpha, but the boundary is then pitted with
#   concavities that eat the volume: measured on a uniformly filled R=150 ball,
#   connectivity-only recovered 21.5% of the true volume (14.3% at n=1500 - it
#   gets WORSE with more points, since denser sampling connects at a still
#   smaller alpha), while adding the watertight test gave 75.1% and 83.8%.
#   Hollow shells are unaffected (96.0% either way), because they only connect
#   once alpha spans the gaps, which is also when they close.
#
# WHAT THE NUMBER MEANS
#   For a hollow structure this is the volume ENCLOSED BY the labels once alpha
#   bridges the gaps between them - the volume of the compartment they bound,
#   NOT the volume of labelled material. It also grows with alpha, so `alpha`
#   is returned alongside and both belong in any report.
cluster_alpha_volume <- function(pts) {
  none <- list(volume = NA_real_, alpha = NA_real_)
  if (is.null(pts)) return(none)
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  if (ncol(pts) < 3L) return(none)
  pts <- pts[stats::complete.cases(pts), 1:3, drop = FALSE]
  storage.mode(pts) <- "double"
  if (nrow(pts) < 4L) return(none)                 # a tetrahedron is the minimum

  # A group flat in any axis makes the tetrahedralization degenerate; nudge by
  # sub-picometre amounts (invisible at nm scale). Same guard as alpha_surface_3d().
  for (k in 1:3) {
    sdk <- stats::sd(pts[, k], na.rm = TRUE)
    if (!is.finite(sdk) || sdk < 1e-9) pts[, k] <- pts[, k] + seq_len(nrow(pts)) * 1e-6
  }

  nnd <- tryCatch({
    d <- RANN::nn2(pts, pts, k = 2L)$nn.dists[, 2]  # distance to nearest OTHER member
    d[is.finite(d) & d > 0]
  }, error = function(e) numeric(0))
  diag_nm <- sqrt(sum((apply(pts, 2, max) - apply(pts, 2, min))^2))
  base <- if (length(nnd)) stats::median(nnd) else diag_nm / 10
  if (!is.finite(base) || base <= 0) base <- diag_nm / 10
  if (!is.finite(base) || base <= 0) return(none)   # all members coincide

  alphas <- c(base * c(1, 1.5, 2, 3, 4, 6, 8, 12, 20, 40), diag_nm)
  alphas <- sort(unique(alphas[is.finite(alphas) & alphas > 0]))
  if (!length(alphas)) return(none)

  as3d <- suppressWarnings(tryCatch(alphashape3d::ashape3d(pts, alpha = alphas),
                                    error = function(e) NULL))
  if (is.null(as3d)) return(none)

  comp <- suppressWarnings(tryCatch(
    alphashape3d::components_ashape3d(as3d, indexAlpha = "all"),
    error = function(e) NULL))
  if (!is.null(comp) && !is.list(comp)) comp <- list(comp)

  # Watertight test: every edge of a closed mesh belongs to exactly two faces;
  # an edge used once is a hole rim. Same criterion as alpha_surface_3d().
  is_closed <- function(tt) {
    if (is.null(tt) || nrow(tt) < 4L) return(FALSE)
    ed  <- rbind(tt[, c(1, 2)], tt[, c(2, 3)], tt[, c(1, 3)])
    ed  <- cbind(pmin(ed[, 1], ed[, 2]), pmax(ed[, 1], ed[, 2]))   # undirected key
    key <- paste(ed[, 1], ed[, 2])
    all(tabulate(match(key, unique(key))) == 2L)
  }
  tri  <- if (!is.null(as3d$triang) && nrow(as3d$triang)) as.matrix(as3d$triang) else NULL
  # per-alpha classification columns are the LAST length(alphas) columns, in the
  # order the alphas were passed
  cls0 <- if (!is.null(tri)) ncol(tri) - length(alphas) else 0L

  for (ii in seq_along(alphas)) {
    cc <- if (!is.null(comp) && length(comp) == length(alphas)) comp[[ii]] else NULL
    # (1) one component covering every member == nothing shed at this radius
    if (!is.null(cc) && !(length(cc) && all(is.finite(cc)) && all(cc == 1L))) next
    if (ii < length(alphas)) {
      # (2) watertight boundary
      if (!is.null(tri) && cls0 >= 3L) {
        cand <- tri[tri[, cls0 + ii] == 2L, 1:3, drop = FALSE]      # REGULAR = boundary
        if (!nrow(cand) || !is_closed(cand)) next
      }
      # (3) the shape must actually CONTAIN the members. Tests (1) and (2) alone
      # admit a degenerate sliver: components_ashape3d() reports connectivity of
      # POINTS, and is_closed() accepts any closed sub-mesh, so a single
      # tetrahedron spanning 4 of the members passes both. Observed on a real
      # 31-member shell - alpha 100.4 gave a 4-triangle shape and a volume of
      # 3,155 nm^3 where the correct rung (alpha 150.5, 58 triangles) gave
      # 10,629,935. Requiring containment is the direct test and works for
      # hollow and solid structures alike, where a "fraction of points on the
      # boundary" rule would not.
      inside <- suppressWarnings(tryCatch(
        alphashape3d::inashape3d(as3d, indexAlpha = ii, points = pts),
        error = function(e) NULL))
      if (!is.null(inside) && !(length(inside) == nrow(pts) && all(inside))) next
    }
    v <- suppressWarnings(tryCatch(
      alphashape3d::volume_ashape3d(as3d, indexAlpha = ii),
      error = function(e) NA_real_))
    if (length(v) && is.finite(v[1]) && v[1] > 0)
      return(list(volume = as.numeric(v[1]), alpha = alphas[ii]))
  }
  none
}

# =============================================================================
# estimate_cluster2_curvature() - PER-STRUCTURE SIZE, CURVATURE AND VOLUME
# =============================================================================
#
# GOAL
#   One row per Median Cluster: how many members, how big (fitted sphere
#   radius), how curved (mean curvature) and what volume its members enclose.
#
# UNIT
#   Members are de-blinked molecules when de-blinking has run, else trace IDs -
#   see cluster2_member_col(). Each member is reduced to its MEDIAN position
#   first, so localization precision enters as sigma/sqrt(n) rather than as
#   the raw point cloud, and a long trace does not outvote a short one.
#
# CURVATURE
#   MEAN curvature H = 1/R in nm^-1, one convention app-wide (see
#   add_cluster2_lines_and_angles). For a fitted sphere this is a
#   reparametrisation of `sphere_radius` and carries no extra information; it is
#   reported because bending is conventionally quoted this way.
#
# WHAT WAS REMOVED AND WHY
#   The old `curvature` / `curvature_deg` columns were lambda3/sum(lambda), the
#   variance fraction on the smallest principal axis. That is a PLANARITY index
#   (the Westin anisotropy family), not curvature: measured, R = 50, 150 and
#   500 nm shells all scored ~0.30, so it could not distinguish a small vesicle
#   from a large one, and a flat disc (0.0002) scored the same as a rod (0.0010),
#   so it could not distinguish a plane from a line. `asin(sqrt(.))` in degrees
#   had no geometric meaning, and `sphere_curvature_deg` was the identical
#   formula under a name implying it came from the sphere fit.
estimate_cluster2_curvature <- function(df, member_col = NULL) {
  mcol <- cluster2_member_col(df, member_col)
  unit_lbl <- if (identical(mcol, "cluster")) "de-blinked molecule" else "trace ID"

  base <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0)
  if (identical(mcol, "cluster"))
    base <- base %>% dplyr::filter(!is.na(cluster) & cluster > 0L)

  base %>%
    dplyr::group_by(cluster2) %>%
    dplyr::group_modify(~{
      med <- .x %>%
        dplyr::group_by(.data[[mcol]]) %>%
        dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                         y = stats::median(y, na.rm = TRUE),
                         z = stats::median(z, na.rm = TRUE), .groups = "drop")
      coords <- as.matrix(med[, c("x", "y", "z")])
      coords <- coords[stats::complete.cases(coords), , drop = FALSE]
      n_mem  <- nrow(coords)

      if (n_mem >= 5) {
        radius <- fit_sphere_geometry(coords)$radius     # geometric-median centre
        av     <- cluster_alpha_volume(coords)
        data.frame(
          member_unit   = unit_lbl,
          n_members     = n_mem,
          sphere_radius = radius,
          mean_curvature_inv_nm = if (is.finite(radius) && radius > 0) 1 / radius else NA_real_,
          volume        = av$volume,
          volume_alpha  = av$alpha
        )
      } else {
        data.frame(
          member_unit   = unit_lbl,
          n_members     = n_mem,
          sphere_radius = NA_real_,
          mean_curvature_inv_nm = NA_real_,
          volume        = NA_real_,
          volume_alpha  = NA_real_
        )
      }
    }) %>%
    dplyr::ungroup()
}


# =============================================================================
# build_curvature_table() - THE EXPORTED PER-STRUCTURE TABLE
# =============================================================================
#
# GOAL
#   One row per Median Cluster, joining membership/identity columns to the
#   geometry from estimate_cluster2_curvature().
#
# COUNTS ARE REPORTED SEPARATELY AND LABELLED
#   n_traces      = blinking events (trace IDs) - photophysics-dependent
#   n_molecules   = de-blinked molecules, NA when de-blinking has not run
#   n_members     = the unit the geometry was actually measured on (see
#                   `member_unit`), so the table never implies the geometry used
#                   a unit it did not.
#   n_points      = raw localizations
build_curvature_table <- function(df) {
  curvature_df <- estimate_cluster2_curvature(df)
  has_mol <- "cluster" %in% names(df) && any(!is.na(df$cluster) & df$cluster > 0L)

  merged <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
    dplyr::group_by(cluster2) %>%
    dplyr::summarize(
      trace_ids   = paste(unique(tid), collapse = ", "),
      molecules   = if (has_mol) paste(unique(cluster[!is.na(cluster) & cluster > 0L]), collapse = ", ") else NA_character_,
      n_traces    = dplyr::n_distinct(tid),
      n_molecules = if (has_mol) dplyr::n_distinct(cluster[!is.na(cluster) & cluster > 0L]) else NA_integer_,
      n_points    = dplyr::n(),
      # Structure position on the same footing as the geometry: the geometric
      # median of the MEMBER centroids, not a coordinate-wise median of every
      # raw localization (which weights long traces and is not rotation-equivariant).
      .groups = "drop"
    )

  mcol <- cluster2_member_col(df)
  centres <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
    { if (identical(mcol, "cluster")) dplyr::filter(., !is.na(cluster) & cluster > 0L) else . } %>%
    dplyr::group_by(cluster2, .data[[mcol]]) %>%
    dplyr::summarise(mx = stats::median(x, na.rm = TRUE),
                     my = stats::median(y, na.rm = TRUE),
                     mz = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(cluster2) %>%
    dplyr::group_modify(~{
      g <- geometric_median(as.matrix(.x[, c("mx", "my", "mz")]))
      data.frame(x = g[1], y = g[2], z = g[3])
    }) %>%
    dplyr::ungroup()

  out <- merged %>%
    dplyr::left_join(centres, by = "cluster2") %>%
    dplyr::left_join(curvature_df, by = "cluster2")

  # ── SMALL-N RELIABILITY FLAG ────────────────────────────────────────────────
  # Every shape metric here is estimated from the cluster's own members, and all of
  # them are biased at low member counts - in DIFFERENT directions, which is why a
  # single "n" column is not enough for a reader to judge them. Measured against
  # known ground truth (see docs and the code notes on estimate_cluster2_curvature):
  #
  #   eccentricity - PCA orders axes by SAMPLE variance, so the ordering itself
  #     manufactures elongation. On a TRUE SPHERE the median ecc2 is 0.81 at n = 10,
  #     0.53 at n = 40, 0.34 at n = 200 and 0.18 at n = 2000. Never read an
  #     eccentricity without comparing it to that floor for the member count.
  #   alpha-shape volume - flat facets under-fill a curved surface: 29.9% of truth
  #     at n = 10, 54.6% at n = 20, 67.6% at n = 31, 86.2% at n = 80, 96.0% at
  #     n = 300. Below ~20 members it is not usable.
  #   sphere_radius - the centre is fitted to the same points, so the radius is
  #     biased low: 147.4 nm at n = 40 against a truth of 150, 149.9 at n = 1000.
  #
  # The flag is advisory and deliberately does not blank anything: which metric the
  # reader cares about decides what counts as "enough", and silently dropping values
  # would hide the size dependence rather than expose it.
  nm_col <- if ("n_members" %in% names(out)) out$n_members else
            if ("n_molecules" %in% names(out) && any(!is.na(out$n_molecules))) out$n_molecules else
            out$n_traces
  out$n_members_used <- nm_col
  out$smallN_flag <- dplyr::case_when(
    is.na(nm_col)  ~ NA_character_,
    nm_col <  20L  ~ "severe (volume <55% of truth; eccentricity floor >0.5)",
    nm_col <  40L  ~ "strong (volume ~68-86% of truth; eccentricity floor ~0.5)",
    nm_col < 200L  ~ "moderate (eccentricity floor ~0.3-0.5)",
    TRUE           ~ "ok"
  )
  # Expected eccentricity of a TRUE SPHERE at this member count, log-interpolated
  # from the measured floor - the number to compare an observed eccentricity against.
  kn <- c(10, 40, 200, 2000); kv <- c(0.81, 0.53, 0.34, 0.18)
  out$ecc_sphere_floor <- ifelse(is.na(nm_col), NA_real_,
    stats::approx(log(kn), kv, xout = log(pmax(nm_col, 2)), rule = 2)$y)
  out
}

## ── Secondary (cluster2) grouping engine ──────────────────────────────────────
# Runs the secondary grouping on a matrix of centroid coordinates (trace or
# de-blinked molecule medians) and returns an integer vector of cluster ids using
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
# Three membership refinements run AFTER the method, in this order:
#
#   core_only    (DBSCAN only) keep only points that are themselves core, i.e.
#                have >= minPts neighbours within eps. Classic DBSCAN also
#                absorbs BORDER points - anything within eps of a core point,
#                however isolated it is itself - which is what attaches lone
#                background centroids to the rim of a real structure. Dropping
#                them is the DBSCAN* cluster definition of Campello et al.
#                (2013), the same definition HDBSCAN uses, so HDBSCAN/OPTICS
#                need no equivalent and ignore this flag.
#
#   prune        graded version of the same idea - see .prune_sparse_members().
#
#   min_members  clusters left with fewer than this many members are dropped
#                whole. Applied LAST so it counts what actually survived, and
#                separate from minPts so the density scale and the "how big is
#                a real structure" decision are not the same knob.
#
# `eps` is used only by dbscan; `xi` only by optics. Works in 2D or 3D depending
# on ncol(pts).
run_cluster2_grouping <- function(pts, method = "dbscan", eps = 200, minPts = 15,
                                  xi = 0.05,
                                  prefilter = FALSE, prefilter_r = NULL,
                                  prefilter_min = NULL,
                                  prune = FALSE, prune_tol = 2,
                                  core_only = FALSE, min_members = 0) {
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

  # ── core-only membership (DBSCAN*) ──
  # Only DBSCAN has border points to drop: HDBSCAN/OPTICS are already defined on
  # core points, and neither has a fixed eps to test against.
  if (isTRUE(core_only) && identical(method, "dbscan")) {
    mp <- max(2L, as.integer(round(minPts)))
    nn <- tryCatch(dbscan::frNN(sub, eps = eps, sort = FALSE), error = function(e) NULL)
    if (!is.null(nn)) {
      is_core <- (lengths(nn$id) + 1L) >= mp        # +1: frNN excludes self, minPts includes it
      lab[!is_core] <- 0L
    }
  }

  # ── optional sparse-member pruning (post-clustering) ──
  if (isTRUE(prune)) lab <- .prune_sparse_members(sub, lab, tol = prune_tol, minPts = minPts)

  # ── minimum members per cluster, on what survived ──
  # as.numeric() first: a restored settings file can deliver this as a string,
  # and round() errors on character rather than coercing.
  mm <- suppressWarnings(as.numeric(min_members %||% 0))
  mm <- if (is.finite(mm)) as.integer(round(mm)) else 0L
  if (is.finite(mm) && mm > 1L) {
    tb <- table(lab[lab > 0])
    small <- as.integer(names(tb)[tb < mm])
    if (length(small)) lab[lab %in% small] <- 0L
  }

  # Renumber survivors to a contiguous 1..k so the documented output contract
  # holds after members and whole clusters have been removed.
  pos <- lab > 0
  if (any(pos)) lab[pos] <- as.integer(factor(lab[pos]))

  out[keep] <- lab
  out
}

# Post-clustering refinement: within each cluster, demote to noise (0) any
# member sitting in a much sparser spot than the members around it.
#
# STATISTIC: each member's CORE DISTANCE - the distance to its (minPts-1)'th
# nearest fellow member, i.e. the radius it needs to reach minPts neighbours.
# Small = deep inside a density maximum, large = on the outskirts. This judges
# density at the SAME scale the clustering used (minPts), and it is a density
# statistic, not a Euclidean distance to a centre or to a mode - so it behaves
# on rings, arcs, rods and helices, where "far from the middle" is a property of
# legitimate members rather than of stragglers.
#
# BASELINE: the median core distance of the member's OWN neighbourhood, not of
# the whole cluster. A single per-cluster median is set by the densest region,
# so on any structure with an internal density gradient it condemns the sparser
# half wholesale even where those members are interior. A local baseline moves
# with the gradient and only flags members that are loose RELATIVE TO THEIR
# SURROUNDINGS.
#
# `tol` is the only knob and the cutoff stays relative (no absolute nm): a
# member is dropped when its core distance exceeds `tol` x its neighbourhood's
# median (tol ~ 2 = "twice as loose as the members around me").
.prune_sparse_members <- function(pts, labels, tol = 2, minPts = 15) {
  if (!length(labels)) return(labels)
  tol <- suppressWarnings(as.numeric(tol %||% 2))
  if (!is.finite(tol) || tol <= 1) return(labels)   # tol<=1 would trim the core too
  mp  <- max(2L, as.integer(round(minPts %||% 15)))
  out <- labels
  for (cl in unique(labels[labels > 0])) {
    idx <- which(labels == cl)
    nm  <- length(idx)
    if (nm < 4L) next                                # no baseline from <4 members
    sub <- pts[idx, , drop = FALSE]

    k  <- max(1L, min(mp - 1L, nm - 1L))             # core distance at the clustering scale
    kd <- tryCatch(dbscan::kNN(sub, k = k)$dist, error = function(e) NULL)
    if (is.null(kd)) next
    cd <- kd[, k]

    m   <- max(4L, min(mp, nm - 1L))                 # baseline neighbourhood
    nbr <- tryCatch(dbscan::kNN(sub, k = m)$id, error = function(e) NULL)
    if (is.null(nbr)) next
    ref <- apply(nbr, 1, function(j) stats::median(cd[j]))   # neighbours', not own

    drop <- is.finite(ref) & ref > 0 & cd > tol * ref
    out[idx[drop]] <- 0L
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
# deblink_axis_weights() - PUT X, Y AND Z ON A COMMON PRECISION SCALE
# =============================================================================
#
# GOAL
#   DBSCAN uses a spherical eps, but MINFLUX trace spread is not isotropic - in
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
#   This is emphatically NOT valid for measuring physical size - cluster radii,
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
#   c(1,1,1) - the caller then behaves exactly as it did before.
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
# Kept as one function so every de-blink call site transforms identically
# - the suggestion and the clustering must share a metric or eps is meaningless.
apply_deblink_weights <- function(pts, w) {
  if (is.null(w)) return(pts)
  # A weight vector shorter than the coordinate matrix used to pad with NA
  # (w[seq_len(3)] on a length-2 w), which poisons every centroid downstream and
  # silently switched de-blinking off. Anything missing or non-finite means
  # "leave that axis alone", which is the isotropic default.
  ww <- rep(1, ncol(pts))
  n  <- min(length(w), ncol(pts))
  if (n) ww[seq_len(n)] <- as.numeric(w)[seq_len(n)]
  ww[!is.finite(ww)] <- 1
  if (all(ww == 1)) return(pts)
  sweep(pts, 2L, ww, "*")
}

# GOAL: the coordinate columns actually present, in canonical order.
# The de-blink family must agree on this: suggest_centroid_threshold() and
# deblink_count_curve() already tolerated 2D input while the two merge functions
# hard-coded c("x","y","z") and errored with "undefined columns selected".
deblink_coord_cols <- function(df) intersect(c("x", "y", "z"), names(df))

# GOAL: per-column median of a coordinate matrix -> one numeric vector.
# Base R only (no matrixStats dependency); the matrices here are 2-3 columns and
# a handful of rows per trace, so a column loop costs nothing. Kept as its own
# function because every de-blink call site must collapse a trace to a centroid
# the same way - median, NA-tolerant, always length ncol even for empty input.
matrixStats_colMedian <- function(m) {
  if (!is.matrix(m)) m <- as.matrix(m)
  storage.mode(m) <- "double"
  if (!nrow(m)) return(rep(NA_real_, ncol(m)))
  vapply(seq_len(ncol(m)), function(j) stats::median(m[, j], na.rm = TRUE),
         numeric(1))
}

# =============================================================================
# link_complete_threshold() - COMPLETE LINKAGE AT A DISTANCE, AT SCALE
# =============================================================================
#
# GOAL
#   Group points so that EVERY pair inside a group is within `thr`, not merely
#   connected through a chain of near neighbours.
#
# WHY NOT JUST hclust()
#   `hclust(dist(M))` needs the full n^2 distance matrix - at 60k traces that is
#   ~1.8e9 pairs and will not run. The way round it is exact, not an
#   approximation: two points in different connected components of the
#   threshold graph are more than `thr` apart, so complete linkage can never
#   merge them. Cutting the global dendrogram at `thr` therefore gives exactly
#   the same partition as cutting each component's dendrogram separately. The
#   components are tiny at de-blinking thresholds, so the n^2 cost is paid only
#   within each of them.
#
# THE GUARD
#   A pathological threshold (or genuinely dense data) can still produce one
#   huge component. Above `max_comp` members the component falls back to its
#   single-linkage answer rather than attempting an allocation that would take
#   the session down - the result is then locally chain-prone, which is the old
#   behaviour, and never an error.
#
# ARGUMENTS
#   M        : numeric matrix, one row per object
#   thr      : distance threshold, same units as M
#   max_comp : component size above which single linkage is used instead
#
# RETURNS  integer group id per ROW of M, all > 0.
# =============================================================================
link_complete_threshold <- function(M, thr, max_comp = 2000L) {
  n <- nrow(M)
  if (!n) return(integer(0))
  if (n == 1L) return(1L)

  # Connected components of the threshold graph, by union-find over frNN.
  fr <- dbscan::frNN(M, eps = thr)
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  for (i in seq_len(n)) {
    nb <- fr$id[[i]]
    if (!length(nb)) next
    for (j in nb) { ra <- find(i); rb <- find(j); if (ra != rb) parent[rb] <- ra }
  }
  comp <- vapply(seq_len(n), find, integer(1))

  out <- integer(n); nxt <- 0L
  for (ix in split(seq_len(n), comp)) {
    if (length(ix) == 1L) { nxt <- nxt + 1L; out[ix] <- nxt; next }
    if (length(ix) > max_comp) {          # too big to resolve - keep the chain
      nxt <- nxt + 1L; out[ix] <- nxt; next
    }
    h <- stats::cutree(stats::hclust(stats::dist(M[ix, , drop = FALSE]),
                                     method = "complete"), h = thr)
    out[ix] <- nxt + as.integer(h); nxt <- nxt + max(h)
  }
  out
}

# =============================================================================
# deblink_merge_centroids() - MERGE TRACES ON CENTROID AGREEMENT
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
#   - it gets WORSE as traces get longer. Centroids have no tails.
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
#   the raw-point eps (~0.6 sigma) - the two methods measure different distances,
#   which is why they have separate inputs rather than sharing one.
#
# NO NOISE BY CONSTRUCTION
#   minPts = 1, so every centroid is a core point and every trace lands in a
#   cluster. That matches the merge-only rule: de-blinking never removes a
#   molecule (see deblink_assign_singletons).
#
# COMPLETE LINKAGE IS THE DEFAULT (measured on the collaborator file, 2026-08-09)
#   Single linkage merges through chains, so a group's diameter is not bounded
#   by `thr`. That is visible in the real data without any ground truth - on
#   250218 collaborator-cond02-B2 at k = 2.5 (thr = 15.2 nm) single linkage produces a
#   group spanning 27.5 nm, 1.8x the threshold and 4.5 sigma across; at k = 5 it
#   spans 73.2 nm. Complete linkage caps the diameter at `thr` by construction
#   (15.1 nm and 30.1 nm for the same two settings).
#
#   Semi-synthetic ground truth built ON that file's own point cloud (real
#   per-trace residuals, real molecule positions, 50% of molecules given a
#   second blink), best F1 over k:
#
#     injected motion   1 nm   2 nm   4 nm   6 nm   9 nm
#     single            0.954  0.944  0.914  0.872  0.786
#     complete          0.955  0.953  0.946  0.920  0.868
#
#   On a fixed, fiducial-corrected sample the motion term is small and the two
#   are nearly tied - the gain at 1-2 nm is within noise. What complete linkage
#   buys at EVERY motion level is tolerance to an over-large threshold, which is
#   the error a user makes when correcting a threshold by hand, which is to
#   overshoot. At 2 nm motion, F1 across k = 3.0 to
#   6.0 falls 0.928 -> 0.715 under single linkage and only 0.952 -> 0.843 under
#   complete. It never scored WORSE at any motion level tested.
#
# ARGUMENTS
#   sub_df  : rows to cluster (needs tid, x, y, z)
#   thr     : centroid-distance threshold in nm (weighted units if w is given)
#   w       : optional per-axis de-blink weights, same metric as everything else
#   linkage : "complete" - every member within thr of every other. The only
#             behaviour the app offers; there is no UI control for this.
#             "single" is the historical chain-through-neighbours rule, kept
#             REACHABLE FROM CODE ONLY so analyses run before 2026-08-09 can be
#             reproduced exactly. It lost to complete linkage at every motion
#             level tested (table above) - do not re-expose it.
#
# RETURNS  integer cluster id per ROW of sub_df, all > 0.
# =============================================================================
deblink_merge_centroids <- function(sub_df, thr, w = NULL,
                                    linkage = c("complete", "single")) {
  linkage <- match.arg(linkage)
  n <- nrow(sub_df)
  if (!n) return(integer(0))
  thr <- suppressWarnings(as.numeric(thr))
  if (!length(thr) || !is.finite(thr) || thr <= 0) thr <- 1e-9

  cols <- deblink_coord_cols(sub_df)
  if (length(cols) < 2L) return(rep(1L, n))
  tid <- as.character(sub_df$tid)
  P   <- apply_deblink_weights(as.matrix(sub_df[, cols, drop = FALSE]), w)

  g <- split(seq_len(n), tid)
  M <- t(vapply(g, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]),
                numeric(ncol(P))))

  # A trace whose centroid is non-finite (e.g. its z is entirely NA) makes
  # dbscan throw. That used to hit the blanket tryCatch below and return
  # "every trace is its own molecule" for the WHOLE dataset - de-blinking
  # silently disabled, with an inflated molecule count and no warning. Instead,
  # cluster the traces that do have a usable centroid and give each unusable one
  # its own id, which is the merge-only rule applied to a trace we cannot place.
  # An axis that is unmeasured across the WHOLE dataset (a z column that is all
  # NA rather than the synthesised 0) is dropped rather than allowed to
  # invalidate every centroid - de-blinking then proceeds in the axes that exist.
  keep_ax <- apply(is.finite(M), 2L, any)
  if (any(keep_ax) && !all(keep_ax)) M <- M[, keep_ax, drop = FALSE]

  ok <- is.finite(rowSums(M))
  cl <- integer(nrow(M))
  if (any(ok)) {
    cl[ok] <- if (identical(linkage, "complete")) {
      tryCatch(link_complete_threshold(M[ok, , drop = FALSE], thr),
               error = function(e) tryCatch(
                 dbscan::dbscan(M[ok, , drop = FALSE], eps = thr, minPts = 1L)$cluster,
                 error = function(e2) seq_len(sum(ok))))
    } else {
      # Single linkage at `thr`: DBSCAN with minPts = 1 is exactly that, and it
      # cannot emit noise.
      tryCatch(dbscan::dbscan(M[ok, , drop = FALSE], eps = thr, minPts = 1L)$cluster,
               error = function(e) seq_len(sum(ok)))
    }
  }
  if (any(!ok)) cl[!ok] <- max(0L, max(cl[ok], 0L)) + seq_len(sum(!ok))

  lab <- stats::setNames(as.integer(cl), names(g))
  unname(lab[tid])
}

# =============================================================================
# suggest_centroid_threshold() - THE MERGE DISTANCE FROM THE TRACE SPREAD
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
#   THE DEFAULT IS k = 2.8, MEASURED AGAINST GROUND TRUTH
#   Ground truth was built on a real file's own point cloud: real per-trace
#   residuals, real molecule positions, half the molecules given a second blink.
#   Scoring which trace PAIRS were correctly grouped, over a range of injected
#   between-blink displacement:
#
#     k        1 nm     2 nm     4 nm      mean
#     2.00    0.928    0.916    0.864     0.903   <- the old default
#     2.355   0.947    0.944    0.897     0.929   <- sigma -> FWHM, still tight
#     2.80    0.954    0.952    0.931     0.946
#     3.00    0.954    0.953    0.942     0.950
#
#   The optimum sits at 2.7-3.1 and the curve is flat across it, so 2.8 is taken
#   as the default rather than the exact peak. Complete linkage bounds a
#   molecule's width by the threshold, so overshooting is far cheaper than
#   undershooting - erring high is the safe direction.
#
#   FWHM (2.355 sigma) was tested as a principled-looking alternative and beats
#   the old 2.0, but is consistently too tight. It also does not describe the
#   right quantity: 2.355 is the width of the single-localization density
#   profile, whereas merging asks whether two CENTROIDS agree.
#
#   Normalising by each trace's own sigma, or by the centroid standard error
#   sigma/sqrt(n), was tested too. Per-trace sigma shrunk toward the global
#   median wins by ~0.004 F1 - real but negligible against the ~0.05 that
#   choosing k buys. sigma/sqrt(n) is markedly WORSE, because it tightens the
#   threshold for long traces while the physical displacement between blinks
#   does not shrink with n.
#
#   The old default assumed the label moves appreciably between blinks. Measured
#   on fiducial-locked SMLM (35 and 105 min runs) there is no coherent sample
#   drift and no dependence of pair separation on the time gap between blinks,
#   so the point-cloud spread - not motion - sets the scale. See
#   sample_drift_estimate(), which flags the rare file where that fails.
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
#        dims) - `within` and `between` are the raw distances the plot draws.
#   NULL when the spread cannot be estimated.
# =============================================================================
suggest_centroid_threshold <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT,
                                       use_z = TRUE, w = NULL, k = 2.8) {
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
  # molecules - the threshold has to land between them.
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
# deblink_threshold_for_channel() - ONE MERGE DISTANCE PER CHANNEL
# =============================================================================
#
# GOAL
#   Return a named list, channel id -> centroid merge distance in nm, so that
#   de-blinking uses a threshold matched to the localization precision of the
#   dye that produced those traces.
#
# WHY NOT ONE SHARED DISTANCE
#   The merge distance is k x the measured trace spread, and trace spread IS the
#   localization precision. Two fluorophores rarely share a photon yield, so they
#   rarely share a precision. A single pooled sigma sits between the two: it
#   over-merges the broader channel and under-merges the tighter one, and the
#   error is invisible because both channels report the same number.
#
#   `k` stays SHARED. It is a modelling choice about how many spreads apart two
#   blinks can be and still be one molecule - a statement about blinking, not
#   about the dye. Only sigma is per channel.
#
# THE FALLBACK IS EXPLICIT, NEVER SILENT
#   A channel with too few long traces to measure a spread (`suggest_centroid_
#   threshold` returns NULL) falls back to the shared `eps_centroid`. The
#   returned attribute `measured` records which channels were measured and which
#   fell back, so the UI can say so rather than showing a number that looks
#   measured but is not.
#
# ARGUMENTS
#   df    : data.frame carrying `channel` (from resolve_channels)
#   input : Shiny input list. Reads eps_centroid, centroid_k, sd_min_locs,
#           deblink_anisotropic, deblink_shared_spread.
#
# RETURNS
#   named list (names = as.character(channel)) of distances in nm, with
#   attr "measured" (named logical) and attr "sigma" (named numeric).
# =============================================================================
deblink_threshold_for_channel <- function(df, input) {
  shared <- suppressWarnings(as.numeric(input$eps_centroid %||% 10))
  if (!length(shared) || !is.finite(shared) || shared <= 0) shared <- 10

  ids <- if ("channel" %in% names(df)) sort(unique(df$channel[!is.na(df$channel)])) else 1L
  out <- stats::setNames(as.list(rep(shared, length(ids))), as.character(ids))
  meas <- stats::setNames(rep(FALSE, length(ids)), as.character(ids))
  sig  <- stats::setNames(rep(NA_real_, length(ids)), as.character(ids))

  # One channel, or the user has asked for a single combined spread: nothing to
  # split, so the shared value stands.
  if (length(ids) < 2L || isTRUE(input$deblink_shared_spread %||% FALSE))
    return(structure(out, measured = meas, sigma = sig))

  kk <- suppressWarnings(as.numeric(input$centroid_k))
  if (!length(kk) || !is.finite(kk) || kk <= 0) kk <- 2.8
  min_n <- input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT

  for (ch in ids) {
    sub <- df[which(df$channel == ch), , drop = FALSE]
    w <- if (isTRUE(input$deblink_anisotropic))
      tryCatch(deblink_axis_weights(sub, min_n = min_n)$w, error = function(e) NULL)
      else NULL
    r <- tryCatch(suggest_centroid_threshold(sub, min_n = min_n, use_z = TRUE,
                                            w = w, k = kk),
                  error = function(e) NULL)
    if (!is.null(r) && is.finite(r$sigma) && r$sigma > 0) {
      out[[as.character(ch)]]  <- kk * r$sigma
      meas[[as.character(ch)]] <- TRUE
      sig[[as.character(ch)]]  <- r$sigma
    }
  }
  structure(out, measured = meas, sigma = sig)
}

# =============================================================================
# sample_drift_estimate() - HOW FAR DID THE SAMPLE MOVE DURING THE RUN?
# =============================================================================
#
# GOAL
#   Put a number, with an interval, on how much the sample moved across the
#   acquisition - so the de-blink radius can be widened on the rare file that
#   drifted, and left on the trace spread everywhere else.
#
# WHY NOT WITHIN-TRACE SLOPES (tried, and it cannot work)
#   The obvious test - regress each axis on time inside a trace and look for a
#   common direction - has no power here. A trace lasts ~0.15 s, so drift of
#   even 0.3 nm/s moves the sample 0.05 nm while it is being recorded, against a
#   localization sigma of ~7 nm. Measured on synthetic data with drift injected
#   from 0 to 0.3 nm/s, the coherence statistic never left its null value.
#   Within-trace slopes see fast jitter only, never slow drift, and scaling such
#   a slope up to the length of the run produced absurd numbers (31 um).
#
# THE LEVER ARM IS THE GAP BETWEEN BLINKS
#   Take pairs of traces whose centroids are close enough to be re-blinks of one
#   molecule. If the sample drifted at velocity v, the later blink of a pair sits
#   at displacement v * dt from the earlier one. Regressing the SIGNED
#   displacement on the SIGNED time gap, through the origin, estimates v with a
#   lever arm of the whole acquisition instead of 0.15 s.
#
#   Pairs that are actually two DIFFERENT molecules do not bias this. Their
#   displacement is isotropic with zero mean, so they widen the interval but
#   cannot tilt the slope. That is why the radius can be generous.
#
# THE ONE REAL BIAS
#   Selecting pairs by separation truncates exactly the large displacements fast
#   drift produces, pulling the estimate toward zero. The radius is therefore set
#   well above the merge threshold, and drift large enough to push same-molecule
#   pairs outside it is UNDERSTATED rather than missed - on a badly drifting file
#   this is a lower bound.
#
# ARGUMENTS
#   df     : localization table (tid, tim, x, y[, z])
#   min_n  : traces shorter than this are ignored
#   radius : pair-selection radius in nm; defaults to 10 x the trace spread
#   conf   : confidence level for the interval on the drift
#
# RETURNS
#   list(rate, rate_lo, rate_hi, drift_run, drift_run_lo, drift_run_hi, dt_typ,
#   drift_typ, jitter, frac_directed, direction, n_units, span, n_pairs, radius,
#   dims, significant) - rates in nm/s, drift_run in nm across the acquisition,
#   drift_typ in nm across a typical gap between blinks, jitter in nm - or NULL
#   when it cannot be measured.
#
#   DRIFT vs JITTER: drift is the displacement that lines up with time and
#   accumulates; jitter is the scatter left once it is removed, which does not.
#   frac_directed says which of the two dominates (0 = pure jitter, 1 = a
#   coherent march). `group` supplies molecule ids (de-blinked IDs) so that a
#   fluorophore seen at several times gives a trajectory instead of one step.
# =============================================================================
sample_drift_estimate <- function(df, min_n = 5L, radius = NULL, conf = 0.95,
                                  group = NULL) {
  if (!is.data.frame(df) || !nrow(df) || !all(c("tid", "tim", "x", "y") %in% names(df)))
    return(NULL)
  cols <- intersect(c("x", "y", "z"), names(df))
  keep <- stats::complete.cases(df[, c("tim", cols), drop = FALSE])
  d <- df[keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  # A synthesised constant z would add a zero-variance axis to every distance.
  cols <- cols[vapply(cols, function(a) stats::sd(d[[a]], na.rm = TRUE) > 0, logical(1))]
  if (length(cols) < 2L) return(NULL)

  gid_row <- if (!is.null(group)) as.character(group)[keep] else NULL

  ord <- split(seq_len(nrow(d)), as.character(d$tid))
  ord <- ord[lengths(ord) >= max(3L, as.integer(min_n))]
  if (length(ord) < 20L) return(NULL)

  M  <- t(vapply(ord, function(ix) matrixStats_colMedian(as.matrix(d[ix, cols, drop = FALSE])),
                 numeric(length(cols))))
  tm <- vapply(ord, function(ix) stats::median(d$tim[ix], na.rm = TRUE), numeric(1))
  # One molecule id per TRACE. Traces never straddle molecules, so the first row
  # of a trace names its molecule.
  tg <- if (!is.null(gid_row)) vapply(ord, function(ix) gid_row[ix[1]], character(1)) else NULL
  ok <- is.finite(rowSums(M)) & is.finite(tm) & (if (is.null(tg)) TRUE else !is.na(tg))
  M  <- M[ok, , drop = FALSE]; tm <- tm[ok]
  if (!is.null(tg)) tg <- tg[ok]
  if (nrow(M) < 20L) return(NULL)

  if (is.null(radius) || !is.finite(radius) || radius <= 0) {
    sg <- stats::median(vapply(ord[ok], function(ix)
      sqrt(mean(apply(as.matrix(d[ix, cols, drop = FALSE]), 2L, stats::var), na.rm = TRUE)),
      numeric(1)), na.rm = TRUE)
    radius <- if (is.finite(sg) && sg > 0) 10 * sg else 60
  }

  # Pair traces that belong to the same molecule. Two routes:
  #   GROUPED - an explicit molecule id (the de-blinked ID). Highest purity, and
  #             a fluorophore seen at three separate times contributes three
  #             pairs, so one long-lived molecule gives a whole mini-trajectory
  #             rather than a single displacement.
  #   SPATIAL - no grouping supplied: pair traces sitting within `radius`.
  # The grouped route is TRUNCATED by whatever threshold produced the grouping:
  # a molecule whose blinks drifted further than that radius was never merged,
  # so it cannot appear here. That is exactly why this is a QC readout and must
  # not feed back into the merge threshold - it would be measuring its own input.
  n_units <- NA_integer_
  if (!is.null(tg)) {
    gs <- split(seq_along(tg), tg)
    gs <- gs[lengths(gs) >= 2L]
    if (!length(gs)) return(NULL)
    cb <- lapply(gs, function(v) utils::combn(v, 2L))
    ii <- unlist(lapply(cb, function(m) m[1, ]), use.names = FALSE)
    jj <- unlist(lapply(cb, function(m) m[2, ]), use.names = FALSE)
    n_units <- length(gs)
  } else {
    nb <- tryCatch(dbscan::frNN(M, eps = radius), error = function(e) NULL)
    if (is.null(nb)) return(NULL)
    ii <- rep(seq_along(nb$id), lengths(nb$id))
    jj <- unlist(nb$id, use.names = FALSE)
    sel <- jj > ii                                 # keep each pair once
    ii <- ii[sel]; jj <- jj[sel]
  }
  if (length(ii) < 20L) return(NULL)

  # Orient every pair the same way in time, so the displacement reads "later
  # minus earlier" and a common drift adds instead of cancelling.
  late  <- ifelse(tm[jj] >= tm[ii], jj, ii)
  early <- ifelse(tm[jj] >= tm[ii], ii, jj)
  dt <- tm[late] - tm[early]
  dP <- M[late, , drop = FALSE] - M[early, , drop = FALSE]
  use <- is.finite(dt) & dt > 0 & is.finite(rowSums(dP))
  dt <- dt[use]; dP <- dP[use, , drop = FALSE]
  if (length(dt) < 20L || sum(dt^2) <= 0) return(NULL)

  # Regression through the origin, per axis: v = sum(dt * dP) / sum(dt^2).
  sdt2  <- sum(dt^2)
  v     <- as.numeric(crossprod(dt, dP)) / sdt2
  resid <- dP - outer(dt, v)
  dfree <- max(length(dt) - 1L, 1L)
  se    <- sqrt(colSums(resid^2) / dfree / sdt2)
  tcrit <- stats::qt(1 - (1 - conf) / 2, dfree)

  span <- diff(range(d$tim, na.rm = TRUE))
  rate <- sqrt(sum(v^2))
  # The displacement that matters is NOT the drift across the whole run - it is
  # the drift across the gap between two blinks of one molecule, which is what
  # the merge radius has to cover. The median gap among candidate pairs is the
  # available proxy for it, and it is typically a fraction of the run.
  dt_typ <- stats::median(dt)
  # Interval on the SPEED: propagate the per-axis errors along the fitted
  # direction, which is the quantity the caller acts on.
  se_r <- if (rate > 0) sqrt(sum((v / rate)^2 * se^2)) else sqrt(sum(se^2))
  lo <- max(0, rate - tcrit * se_r); hi <- rate + tcrit * se_r

  # DRIFT vs JITTER. Drift is the part of the displacement that lines up with
  # time - one direction, growing with the gap. Jitter is what is left once that
  # is removed: how far two blinks of one molecule land apart for reasons that do
  # NOT accumulate. Both are reported because they mean different things about a
  # sample, and because jitter is the term that actually sets how wide a de-blink
  # radius has to be.
  #   jitter here still contains the localization error of the two centroids, so
  #   it is an upper bound on any real positional wobble.
  jitter <- sqrt(mean(rowSums(resid^2)))

  # Displacements FROM EACH MOLECULE'S FIRST BLINK, for the drift plot: every
  # molecule starts at the origin and each later blink is one point. This is a
  # subset of the pairs used for the fit (the fit also uses second-to-third
  # steps, which carry information but have no common origin to plot against),
  # so the plot shows the same phenomenon on a readable frame rather than the
  # exact rows the regression saw.
  if (!is.null(tg)) {
    firsts <- vapply(split(seq_along(tg), tg), function(v) v[which.min(tm[v])], integer(1))
    is_first <- logical(length(tm)); is_first[firsts] <- TRUE
    ssel <- is_first[early[use]]
  } else {
    ssel <- rep(TRUE, length(dt))
  }
  steps <- if (any(ssel)) {
    st <- as.data.frame(dP[ssel, , drop = FALSE])
    names(st) <- paste0("d", cols)
    st$dt <- dt[ssel]
    # Which molecule each step belongs to, so the caller can colour by it and
    # join one molecule's blinks into a trajectory. NA when there is no grouping
    # (spatial pairing has no molecule to name).
    st$mol <- if (!is.null(tg)) tg[early[use]][ssel] else NA_character_
    # WHEN in the acquisition this blink happened, as distinct from dt (how long
    # AFTER its molecule's first blink it happened). Colouring by this shows
    # whether displacements point one way early in the run and another way
    # later - a reversal or a step that a gap-based colour cannot reveal.
    st$t_abs <- tm[late[use]][ssel]
    st[order(st$mol, st$dt), , drop = FALSE]
  } else NULL
  # Directional purity: how much of the total displacement the fitted drift
  # explains. Near 0 = scattered (jitter only), near 1 = a coherent march.
  tot <- sum(dP^2)
  frac_dir <- if (tot > 0) max(0, 1 - sum(resid^2) / tot) else NA_real_

  list(rate = rate, rate_lo = lo, rate_hi = hi,
       drift_run = rate * span, drift_run_lo = lo * span, drift_run_hi = hi * span,
       dt_typ = dt_typ, drift_typ = rate * dt_typ,
       jitter = jitter, frac_directed = frac_dir, n_units = n_units,
       direction = if (rate > 0) v / rate else rep(NA_real_, length(cols)),
       velocity = stats::setNames(v, cols), axes = cols, steps = steps,
       span = span, n_pairs = length(dt), radius = radius, dims = length(cols),
       significant = is.finite(lo) && lo > 0)
}


# =============================================================================
# deblink_merge_majority_vote() - LET THE LOCALIZATIONS VOTE
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
#   too small - every point is unfairly close to its own centre. The effect goes
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
#   Votes are asymmetric - a short diffuse trace may vote for a long tight one
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
deblink_merge_majority_vote <- function(sub_df, f = 0.3, rho = 3.0, mutual = TRUE,
                                        w = NULL, k_cand = 4L) {
  n <- nrow(sub_df)
  if (!n) return(integer(0))
  f   <- min(max(as.numeric(f   %||% 0.3), 0.01), 1)
  rho <- max(as.numeric(rho %||% 3.0), 0.01)

  cols <- deblink_coord_cols(sub_df)
  if (length(cols) < 2L) return(rep(1L, n))
  tid  <- as.character(sub_df$tid)
  P    <- apply_deblink_weights(as.matrix(sub_df[, cols, drop = FALSE]), w)
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
  # individual candidate clears f - the trace stays unmerged even though nearly
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
  # points - otherwise siblings cancel each other out again
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
  # and ~0.07 at 20 nm separation. So f must sit well BELOW 0.5 - thresholding
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
# deblink_minpts() - minPts IS THE TRACE-LENGTH FILTER, NOT A FREE KNOB
# =============================================================================
# GOAL
#   Derive this diagnostic's minPts from the Trace ID Length filter's lower
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
#   so the only thing that can still leave it unmerged is its own spread -
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
# deblink_assign_singletons() - DE-BLINKING MERGES, IT NEVER REMOVES
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
#   filters `cluster != 0`, so leaving traces at 0 would remove them silently -
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

  # One new ID per TRACE (not per point) - a noise trace stays one molecule.
  # Keyed on (tid, channel) for the same reason the plurality step is: under a
  # per-localization channel split, one trace can sit in two channels, and a
  # single shared ID would be a molecule spanning both.
  key <- if ("channel" %in% names(df))
    paste0(as.character(df$tid[noise]), "
", as.character(df$channel[noise]))
    else as.character(df$tid[noise])
  df$cluster[noise] <- as.integer(nxt) + as.integer(factor(key))
  df
}


# =============================================================================
# deblink_count_curve() - MOLECULE COUNT vs RADIUS, AND ITS KNEE   [EXPLORATORY]
# =============================================================================
#
# GOAL
#   Pick the de-blink radius from the data alone, with no model of precision and
#   no ground truth: sweep eps, run the REAL merge at each value, and plot how
#   many molecules come out.
#
# THE SHAPE OF THE CURVE AND WHY THE KNEE IS THE ANSWER
#   small eps  count is inflated - re-blinks of one fluorophore stay separate
#              and long traces fragment. The count falls STEEPLY as eps grows
#              and genuine re-blinks start to merge.
#   the knee   the point where that steep fall stops: the re-blinks have been
#              collected and nothing much is left to merge correctly.
#   large eps  a slow, steady decline - this is DISTINCT molecules being fused.
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
#   when the two regimes are well separated and unreliable when they are not -
#   which is the same crowding condition that makes de-blinking ambiguous in the
#   first place. Treat it as evidence, not as an answer, and look at the curve.
#
# ARGUMENTS
#   df       : filtered data.frame (x, y, z, tid)
#   minpts   : accepted for signature compatibility and currently UNUSED. The
#              sweep counts connected components at each radius (minPts = 1), so
#              there is no density criterion to match. De-blinking itself has no
#              minPts either: it merges on a distance threshold.
#   eps_seq  : radii to try (nm)
#   w        : optional per-axis de-blink weights, so the curve is measured in
#              the same metric the clustering will use
#   max_pts  : subsample cap - the sweep runs DBSCAN once per eps
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
  # a distance - the caller supplies the right sequence.
  Mc <- NULL
  if (identical(method, "centroid")) {
    gg <- split(seq_len(nrow(P)), tid)
    Mc <- t(vapply(gg, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]),
                   numeric(ncol(P))))
  }
  # CENTROID ONLY. A vote sweep was tried and withdrawn: it did not return in a
  # usable time, and an unresponsive diagnostic is worse than an absent one.
  # The vote criterion has no distance to sweep anyway - its parameter is a
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
# apply_clustering_pipeline() - THE TWO-TIER CLUSTERING ENTRY POINT
# =============================================================================
#
# GOAL
#   Turn raw localizations into the two grouping levels the whole SMLM tab is
#   built on. This mirrors exactly what the interactive plot does
#   (`plottable_data()` in server_smlm.R), so batch tools and the on-screen plot
#   can never disagree - that duplication is the reason this function exists.
#
# THE TWO TIERS AND WHY THERE ARE TWO
#   DE-BLINKING (`cluster`).
#     Merges TRACE CENTROIDS - one point per trace, not the raw localizations -
#     when they sit within a threshold of a few times the MEASURED TRACE SPREAD
#     (~10-30 nm). A single fluorophore blinking on and off produces several
#     traces in the same place; this collapses them into one object. Without it,
#     one molecule is counted many times and every downstream count and spatial
#     statistic is inflated.
#
#     ⚠ THIS IS NOT DBSCAN, and it has not been for some time. The two criteria
#     are complete-linkage on the centroids (`deblink_merge_centroids`, the
#     default) and a majority vote of each trace's localizations
#     (`deblink_merge_majority_vote`). Neither has a minPts and neither has a
#     density criterion, so neither can emit noise: every trace comes out with
#     an id and `cluster` has no 0. The single-linkage option calls
#     `dbscan::dbscan(minPts = 1)`, but that is exactly single-linkage connected
#     components at the threshold, used as a components engine rather than as
#     density clustering. Older comments called this stage "primary DBSCAN";
#     they were wrong and have been corrected. The input VALUES "primary" and
#     "pcid" survive in saved settings and still mean the de-blinked level.
#
#   MEDIAN CLUSTERING (`cluster2`) - STRUCTURE.
#     Runs over the MEDIAN position of each unit from the level above, with a
#     LARGE eps (~100-300 nm). Its job is biological: group molecules into the
#     structure they form (a filament, a ring, a complex). THIS is where DBSCAN
#     actually runs, and it is the only density-clustering stage in the pipeline.
#
#   The two distances are doing completely different jobs - the de-blink
#   threshold is set by the MICROSCOPE (precision + label size), cluster2's eps
#   by the BIOLOGY (how far apart molecules in one structure sit). Tuning them
#   together is the most common source of confusion.
#
# WHY DBSCAN FOR THE MEDIAN LEVEL (and its alternatives here)
#   DBSCAN needs no pre-set cluster count, finds arbitrary shapes, and labels
#   sparse points as noise (cluster 0) rather than forcing them into a group -
#   all three matter for SMLM data. It does assume roughly uniform density, so
#   the median level also offers:
#     HDBSCAN ... adaptive density; only a minimum cluster size, no eps
#     OPTICS .... adaptive; keeps small structures a single eps would merge
#
# THE CHANNEL RULE (enforced, not optional) - BOTH TIERS
#   Both tiers run SEPARATELY PER CHANNEL, with IDs offset so they can never
#   collide. Channels come from resolve_channels(): a 2-colour DCR/fluo split, a
#   multiplex build, or a single channel when neither is active.
#
#   Molecules of different colours must never be de-blinked into one molecule,
#   and molecules imaged in different exchange rounds must never be merged into
#   one cluster - the second would invent a structure that was never
#   simultaneously present. These are one rule, so there is one loop.
#
#   There is no pooled option. Pooling across channels is the thing the channel
#   column exists to prevent; `cluster_group_scope` chooses between "within"
#   (every channel, separately) and "only:<k>" (one channel, for inspection).
#
#   PARAMETERS SPLIT WHERE THE PHYSICS SPLITS. The primary merge distance is per
#   channel, because it is k x the trace spread and spread is a property of the
#   dye (see deblink_threshold_for_channel). eps2 is shared, because it is set by
#   the biology - how far apart molecules in one structure sit - which no choice
#   of fluorophore changes.
#
# DE-BLINKING MERGES, IT NEVER REMOVES
#   The primary tier has no "drop unclustered" behaviour and no free minPts:
#     - minPts comes from the Trace ID Length filter (see deblink_minpts), so it
#       cannot contradict the filter the user actually set.
#     - traces DBSCAN cannot merge become their own single-blink molecule (see
#       deblink_assign_singletons), so `cluster` is > 0 everywhere.
#   Removing molecules is the job of the metadata and spread filters, where it
#   is visible and reversible - not a side effect of a clustering radius.
#
# ARGUMENTS
#   df    : an ALREADY-FILTERED data.frame (x, y, z, tid [, mplex_channel])
#   input : the Shiny input list - read directly so batch tools use exactly the
#           settings the user sees. Relevant: enable_clustering, eps_centroid,
#           centroid_k, deblink_method, deblink_anisotropic, deblink_shared_spread,
#           tid_length_range (-> minPts), enable_cluster2, eps2/minpts2,
#           cluster2_method, cluster_group_scope, enable_two_color (-> channels).
#
# RETURNS
#   the same df with `channel`/`channel_label`/`channel_color` (from
#   resolve_channels) and `cluster` and/or `cluster2` integer columns appended.
#   For `cluster2`, ID 0 still means "noise / unassigned" and consumers must
#   exclude it before counting. For `cluster`, 0 no longer occurs - the existing
#   `cluster != 0` filters downstream are harmless no-ops that now keep every
#   trace which survived filtering.
# =============================================================================
apply_clustering_pipeline <- function(df, input) {
  if (!nrow(df)) return(df)

  # One channel definition for the whole run. Everything below loops over it.
  df <- resolve_channels(df, input)
  if (!nrow(df)) return(df)

  chans <- channels_in_scope(df, input)

  # Run `fn` once per channel, offsetting the returned IDs so they can never
  # collide across channels. Running-max, not a fixed constant: the 2-colour
  # path used +10000L, which breaks silently at 3 channels or 10k clusters.
  per_channel <- function(n_rows, subset_idx, fn) {
    out <- integer(n_rows)
    next_id <- 1L
    for (ch in chans) {
      idx <- subset_idx(ch)
      if (!length(idx)) next
      cid <- fn(idx)
      if (!length(cid)) next
      cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
      out[idx] <- cid
      if (any(cid > 0)) next_id <- max(next_id, max(cid, na.rm = TRUE) + 1L)
    }
    out
  }

  # ── De-blinking ─────────────────────────────────────────────────────
  if (isTRUE(input$enable_clustering)) {
    df$cluster <- 0L

    # Per-axis metric, computed ONCE over the whole df so every channel subset
    # below is clustered in the same space (a per-subset metric would make
    # cluster IDs incomparable across channels).
    aw1 <- if (isTRUE(input$deblink_anisotropic))
      deblink_axis_weights(df, min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT)$w
      else NULL

    # Two merge criteria, both trace-level. "centroid" compares where the traces
    # ARE; "vote" asks whether most of a trace's localizations think they belong
    # to a neighbouring trace. Raw-point DBSCAN was removed: one stray
    # localization could bridge two molecules. Ground-truth errors at each
    # method's own best threshold: raw-point 35, vote 29, centroid 10.
    # (Centroid scores 11 at the k = 2.0 default and 10 at its optimum 1.54 -
    #  quote the basis, not a bare number.)
    method1 <- input$deblink_method %||% "centroid"

    # THE MERGE DISTANCE IS PER CHANNEL (see deblink_threshold_for_channel).
    # It is dictated by the trace spread, and trace spread is a property of the
    # dye: two fluorophores with different photon yields have different
    # localization precision, so one shared distance under-serves the tighter
    # channel. The shared value remains available as an explicit fallback.
    thr_by_ch <- deblink_threshold_for_channel(df, input)

    run_primary <- function(idx) {
      sub <- df[idx, , drop = FALSE]
      if (identical(method1, "vote"))
        return(deblink_merge_majority_vote(sub,
                 f = input$vote_frac %||% 0.2, rho = input$vote_rho %||% 1.2,
                 mutual = !identical(input$vote_link %||% "mutual", "either"), w = aw1))
      ch  <- sub$channel[1]
      thr <- thr_by_ch[[as.character(ch)]] %||% (input$eps_centroid %||% 10)
      deblink_merge_centroids(sub, thr, aw1)
    }

    df$cluster <- per_channel(nrow(df),
                              function(ch) which(df$channel == ch),
                              run_primary)

    # TID-plurality reassignment, keyed on (tid, CHANNEL) rather than tid alone.
    #
    # THE KEY MATTERS, and getting it wrong silently undoes the channel loop.
    # A 2-colour split is assigned per TRACE, so there (tid, channel) and tid are
    # the same key and this changes nothing. A MULTIPLEX split is assigned per
    # LOCALIZATION, so a trace whose blinks straddle a TIM breakpoint belongs to
    # two channels at once. Keyed on tid alone, plurality then pulled all of that
    # trace into whichever channel held the majority - re-merging across the
    # boundary the loop above had just respected. Measured on a 3-window TIM split
    # of a real file: 2 de-blinked molecules spanned two channels, after the
    # per-channel clustering had produced none.
    #
    # Splitting the trace instead is the honest reading: the localizations really
    # were recorded in different exchange rounds, so they are not evidence of one
    # simultaneously-present molecule. It does mean a straddling trace is counted
    # twice - see the boundary-splitting note in docs/channel_unification_plan.md.
    tid_plurality <- df %>%
      dplyr::filter(cluster != 0L) %>%
      dplyr::group_by(tid, channel, cluster) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(tid, channel) %>%
      dplyr::slice_max(n, n = 1L, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(tid, channel, cluster_new = cluster)
    df <- df %>%
      dplyr::left_join(tid_plurality, by = c("tid", "channel")) %>%
      dplyr::mutate(cluster = dplyr::coalesce(cluster_new, 0L)) %>%
      dplyr::select(-cluster_new)

    # Merge-only: unmerged traces become their own molecule rather than being
    # dropped. See deblink_assign_singletons().
    df <- deblink_assign_singletons(df)
  } else {
    df$cluster <- NA_integer_
  }

  # ── Secondary clustering (structure) ──────────────────────────────────────
  if (isTRUE(input$enable_cluster2)) {
    eps2    <- input$eps2    %||% 200
    minpts2 <- input$minpts2 %||% 3
    use_primary <- isTRUE(input$enable_clustering) && "cluster" %in% names(df)

    # eps2 is SHARED across channels: it is set by the biology (how far apart
    # molecules in one structure sit), not by the dye. The clustering still runs
    # strictly per channel, so a Median Cluster can never hold molecules from
    # two channels.
    if (use_primary) {
      medians <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0L) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          id = dplyr::first(cluster),
          mx = stats::median(x, na.rm = TRUE),
          my = stats::median(y, na.rm = TRUE),
          mz = stats::median(z, na.rm = TRUE),
          channel = dplyr::first(channel),
          analysis_group = {
            v <- analysis_group[!is.na(analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
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
          channel = dplyr::first(channel),
          analysis_group = dplyr::first(analysis_group),
          .groups = "drop"
        )
    }
    medians <- medians[is.finite(medians$mx) & is.finite(medians$my) &
                       is.finite(medians$mz), , drop = FALSE]

    run_secondary <- function(idx) {
      if (length(idx) < minpts2) return(integer(length(idx)))
      # DIMENSIONALITY IS A PROPERTY OF THE DATA, NOT OF THE VIEW. Median
      # clustering always measures distances in 3D. Following the `mode3d`
      # DISPLAY toggle here meant flipping the view to 2D silently re-clustered
      # the data. A placeholder z (all zeros) contributes nothing to a Euclidean
      # distance, so always-3D is safe for 2D data; on real z it is simply
      # correct. Measured on real MINFLUX files (z sd 112-138 nm), 3D
      # nearest-neighbour distances are 1.46-1.89x the 2D projection, so the
      # choice moved eps2 by nearly a factor of two.
      pts <- as.matrix(medians[idx, c("mx", "my", "mz")])
      run_cluster2_grouping(
        pts,
        method        = input$cluster2_method %||% "dbscan",
        eps           = eps2,
        minPts        = minpts2,
        xi            = input$cluster2_xi %||% 0.05,
        prefilter     = isTRUE(input$cluster2_prefilter),
        prefilter_r   = input$cluster2_prefilter_r,
        prefilter_min = input$cluster2_prefilter_min,
        prune         = isTRUE(input$cluster2_prune),
        prune_tol     = input$cluster2_prune_tol %||% 2,
        core_only     = isTRUE(input$cluster2_core_only),
        min_members   = input$cluster2_min_members %||% 0
      )
    }

    medians$cluster2 <- per_channel(nrow(medians),
                                    function(ch) which(medians$channel == ch),
                                    run_secondary)

    join_key <- if (use_primary) c("cluster" = "id") else c("tid" = "id")
    df <- df %>%
      dplyr::select(-dplyr::any_of("cluster2")) %>%
      dplyr::left_join(medians[, c("id", "cluster2")], by = join_key)
  } else {
    df$cluster2 <- NA_integer_
  }

  df
}

## ── Empirical label-multiplicity calibration (counting standard) ──────────────
# Measures the observed labelling multiplicity from a CALIBRATION STANDARD of
# known stoichiometry, so the count estimator can divide by a measured value
# instead of the preset guess m_tid = n_fluor * eff * (1 + blink_lambda).
#
# The standard is a sparse sample of isolated structures each holding a KNOWN
# number of true molecules (`k_per_structure`; 1 for single, well-separated
# molecules such as sparse DNA docking strands). It MUST be imaged with the same
# label and acquisition as the sample being counted, and clustered with the same
# settings, so each secondary structure (cluster2 > 0) is one standard object.
#
# Observed TIDs per molecule = median over structures of (TIDs in the structure)
# / k. That IS the m_tid the count divides by, measured rather than assumed. The
# per-structure spread (MAD) reports how reproducible the label is.
#
# Returns a list; ok = FALSE (never an error) when nothing is measurable.
measure_label_multiplicity <- function(df, k_per_structure = 1L, minpts = 3L) {
  out <- list(ok = FALSE, n_struct = 0L, k = as.integer(k_per_structure),
              m_tid = NA_real_, m_tid_mad = NA_real_, m_tid_iqr = c(NA_real_, NA_real_),
              m_prim = NA_real_, m_loc = NA_real_, per_struct = data.frame())
  if (is.null(df) || !nrow(df) || !"tid" %in% names(df) || !"cluster2" %in% names(df))
    return(out)
  k <- max(1L, as.integer(k_per_structure))
  has_prim <- "cluster" %in% names(df)
  d <- df[!is.na(df$cluster2) & df$cluster2 > 0L, , drop = FALSE]
  if (!nrow(d)) return(out)
  per <- d %>%
    dplyr::group_by(cluster2) %>%
    dplyr::summarise(
      n_loc  = dplyr::n(),
      n_tid  = dplyr::n_distinct(tid),
      n_prim = if (has_prim) dplyr::n_distinct(cluster[!is.na(cluster) & cluster > 0L]) else NA_integer_,
      .groups = "drop") %>%
    dplyr::filter(n_loc >= minpts)
  if (!nrow(per)) return(out)
  out$ok         <- TRUE
  out$n_struct   <- nrow(per)
  out$k          <- k
  out$m_tid      <- stats::median(per$n_tid, na.rm = TRUE) / k
  out$m_tid_mad  <- stats::mad(per$n_tid / k, na.rm = TRUE)
  out$m_tid_iqr  <- as.numeric(stats::quantile(per$n_tid / k, c(0.25, 0.75), na.rm = TRUE))
  out$m_prim     <- if (has_prim) stats::median(per$n_prim, na.rm = TRUE) / k else NA_real_
  out$m_loc      <- stats::median(per$n_loc, na.rm = TRUE) / k
  out$per_struct <- as.data.frame(per)
  out
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
# and the DE-BLINKING step (merge threshold ~ dye_spread) collapses re-blinks and
# the dyes of one molecule back into ~one object => ~one unit per molecule,
# turning an unstable /m_tid calibration into a /~1 count. De-blinked molecules
# are therefore the preferred unit when present ("auto").
#
# KEY LIMIT: the link_nm offset is COMMON-MODE within a molecule (all its
# localizations share the one draw), so it never averages out - it is the
# irreducible position-uncertainty floor. Reported as `pos_sd_nm`.
#
# Args:
#   df       clustered SMLM frame (needs x,y[,z],tid; uses `cluster`/`cluster2`)
#   preset   a preset list from minflux_label_presets(), OR its name, OR a list
#            with n_fluor/link_nm/dye_spread_nm/footprint_nm/blink_lambda
#   eff_pct  labeling/detection efficiency (%) - folds into the multiplicity
#   level    "auto" (de-blinked molecules when available else TIDs), "primary", "tid"
#            ("primary" is the stored VALUE for the de-blinked level; kept for
#            settings compatibility, it does not mean a primary DBSCAN)
#   fuse_k   a unit whose per-axis spread exceeds fuse_k x the single-molecule
#            expectation is flagged as fused (multiple GTs merged)
#   minpts   localizations below which a unit is ignored
#   sigma_loc_nm  per-localization precision; estimated from data when NULL
#   scope_cluster2 restrict to points inside a secondary structure (cluster2>0)
#
# Returns list(points, summary, meta) - see field construction below. Empty
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
                                      support_strength = 1,
                                      m_tid_override = NULL) {
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
  # NOTE: compute the SDs BEFORE the median columns overwrite x/y/z - dplyr
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
  # hence several de-blinked molecules / TIDs, around ONE molecule. Pool unit
  # centroids that fall within `complex_r` (the label-complex radius) so those
  # contributors collapse into a single candidate GT - this is what makes
  # "many traces -> one molecule" happen, and it produces the contribution
  # fan-out lines. Merging is confined within each cluster2 so distinct
  # structures never pool together.
  #
  # AUTO radius is PRESET-SPECIFIC: max(footprint, 2.5 x dye_spread). Both terms
  # come from the labeling preset, so a big polyclonal+secondary complex pools
  # over ~45 nm while a nanobody pools over ~3 nm - i.e. the labeling strategy
  # itself sets how many traces collapse into one molecule, as expected.
  if (is.null(complex_r) || !is.finite(complex_r) || complex_r < 0)
    complex_r <- max(footprint, 2.5 * dye_spread)
  coord_cols <- if (use_z) c("x", "y", "z") else c("x", "y")

  # `m_tid` = expected TIDs per molecule from the labeling strategy; the
  # stoichiometric molecule count of a structure is (its TID count) / m_tid.
  # A measured multiplicity from a counting standard (m_tid_override) replaces the
  # preset guess when supplied, so the count divides by a calibrated value.
  m_tid_measured <- !is.null(m_tid_override) && is.finite(m_tid_override) && m_tid_override > 0
  tids_per_mol <- if (m_tid_measured) as.numeric(m_tid_override) else max(1e-6, m_tid)

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
    # clustering. Only ever REDUCES the count toward stoichiometry - it never
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
  # over-wide de-blinked molecule that already chained many traces together.
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
  # (the axial component is dropped - there is no inside along the length); for a
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
  #  s_link    linker floor - irreducible, common-mode per molecule (dominant)
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
    tids_per_mol = tids_per_mol,
    m_tid_source = if (m_tid_measured) "measured" else "preset",
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
# length from the backbone's minimum spanning tree. Not tied to any biology -
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

    # Shape class - dimension-aware. In 2D there is no sheet/sphericity axis, so
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
# GOAL: radius of gyration  Rg = sqrt(<|p - centroid|^2>)  per cluster - the RMS
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
# remains is the label footprint - real between-trace spread - which is exactly the
# thing the "cluster" grouping is meant to report and the "cluster2" grouping is
# meant to look past by de-blinking to one point per molecule.)
#
# SMALL-N BIAS (bias_correct): the centroid is estimated from the SAME points whose
# spread is being measured, so it sits inside the cloud and E[Rg²_obs] = Rg²·(N−1)/N.
# That is a 29% underestimate at N=2, 11% at N=5, 1% at N=50 - i.e. an N-DEPENDENT
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
## GOAL: a REPRODUCIBLE subsample for the O(n^2) spatial statistics.
##
## These estimators cap the point count (pairwise distances are quadratic), and the
## cap used a bare `sample.int`, so pressing Run twice on the same data returned
## DIFFERENT curves whenever a file exceeded the cap. That is not acceptable for a
## figure, and it made cached and uncached results compare unequal.
##
## A seeded draw rather than a sorted head: taking the first k rows would be biased
## whenever the table is written in acquisition or spatial order, which MINFLUX
## exports are. The global RNG stream is saved and restored so seeding here cannot
## perturb the CSR envelopes, permutation nulls or thinning draws.
.stable_subsample <- function(n, cap, seed = 20260813L) {
  if (!is.finite(n) || n <= cap) return(seq_len(n))
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  bk  <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", bk, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  sort(sample.int(n, cap))
}

## `support`: optional support object (estimate_support(..., return_mask = TRUE)).
## When given, its measure replaces the bounding-box area in lambda_hat = n / A. See
## the SUPPORT-AWARE OBSERVATION WINDOW banner for why the box is wrong.
compute_ripley_L <- function(x, y, z = NULL, r_seq = NULL, use_3d = FALSE,
                             max_subsample = 2000L, window_ranges = NULL,
                             support = NULL, edge_correct = TRUE) {
  coords <- cbind(x, y)
  if (use_3d && !is.null(z)) coords <- cbind(x, y, z)

  n <- nrow(coords)
  if (n < 3L) return(data.frame(r = numeric(0), Lr_minus_r = numeric(0)))

  # subsample for large datasets
  idx <- .stable_subsample(n, max_subsample)
  if (length(idx) < n) { coords <- coords[idx, , drop = FALSE]; n <- nrow(coords) }

  # observation window (bounding box). A pre-supplied window_ranges (a 2×d matrix
  # of per-axis min/max, as from apply(coords, 2, range)) lets the observed curve
  # and the Monte-Carlo CSR envelope share one identical window so the edge bias
  # cancels between them.
  rng <- if (!is.null(window_ranges)) window_ranges else apply(coords, 2, range)
  # Support measure when supplied, else the bounding box (the historical behaviour).
  area_vol <- if (!is.null(support)) support_area(support, rng) else prod(rng[2, ] - rng[1, ])
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

  # Translation edge correction: each PAIR is weighted by 1/gammabar at its own
  # separation, since K(r) = integral over |h| <= r of g(h) dh accumulates pairs of
  # every separation below r. gammabar = area_vol for all pairs recovers the
  # uncorrected estimator exactly.
  # SINGLE-PASS accumulation. The old form swept every pair distance once PER RADIUS
  # (60 radii x ~2e6 pairs), which dominated the cost; binning each distance once and
  # cumulating is one pass. `left.open` makes bin k hold the pairs whose distance is
  # the smallest r_seq value that is >= d, so cumsum(k) counts exactly d <= r_seq[k].
  nb    <- length(r_seq)
  cnt_b <- tabulate(findInterval(d_vals, r_seq, left.open = TRUE) + 1L, nbins = nb)

  # Edge weight applied PER BIN rather than per pair. gamma-bar is already the
  # direction-averaged (isotropised) set covariance, so resolving it finer than the
  # r bins buys nothing - and interpolating it at ~2e6 individual pair distances per
  # simulation was the remaining bottleneck inside the envelope loop.
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d && !is.null(z)) else NULL

  Kr <- if (is.null(gb)) {
    # unbiased K: multiply by 2 for both directions, divide by n^2/area
    (area_vol * 2 * cumsum(cnt_b)) / (n * (n - 1))
  } else {
    w <- ifelse(is.finite(gb) & gb > 0, area_vol / gb, 1)  # 1 = no area lost
    (area_vol * 2 * cumsum(cnt_b * w)) / (n * (n - 1))
  }

  Lr <- if (use_3d && !is.null(z)) {
    (3 * Kr / (4 * pi))^(1/3)
  } else {
    sqrt(Kr / pi)
  }

  data.frame(r = r_seq, Kr = Kr, Lr = Lr, Lr_minus_r = Lr - r_seq)
}

## Two-colour clustering used to live here as apply_two_color_clustering_pipeline(),
## a near line-by-line copy of apply_clustering_pipeline() with the group loop
## hardcoded to c(1L, 2L). It has been folded into that function: channels come
## from resolve_channels() and the loop covers 1..N, so 2-colour and multiplex
## are the same code path and neither can merge across a channel boundary.


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
      # "cluster2" collapses each Median Cluster to ONE point, so the distance is
      # structure-to-structure ("is my group-1 cluster near a group-2 cluster")
      # rather than molecule-to-molecule. Both channels keep their own clusters -
      # nothing here mixes the two groups.
      want   <- if (u == "cluster2") "cluster2" else "cluster"
      cl_col <- if (want %in% names(d) && any(!is.na(d[[want]]) & d[[want]] > 0L)) want else NULL
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
                               max_subsample = 2000L, window_ranges = NULL,
                               support = NULL, edge_correct = TRUE) {
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

  pts1 <- pts1[.stable_subsample(n1, max_subsample), , drop = FALSE]
  pts2 <- pts2[.stable_subsample(n2, max_subsample, seed = 20260814L), , drop = FALSE]
  n1 <- nrow(pts1); n2 <- nrow(pts2)

  all_pts <- rbind(pts1, pts2)
  # per-axis extents; a pre-supplied window_ranges (2×d min/max matrix) makes the
  # observed curve and its CSR envelope share one identical window.
  rng <- if (!is.null(window_ranges)) (window_ranges[2, ] - window_ranges[1, ])
         else apply(all_pts, 2, function(v) diff(range(v, na.rm = TRUE)))
  if (!is.null(support)) {
    area_vol <- support_area(support)
    if (!is.finite(area_vol) || area_vol <= 0) return(data.frame())
  } else if (use_3d) {
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

  # Translation edge correction, identical in form to the univariate case: divide by
  # the window's isotropised set covariance instead of its area. Applied PER SHELL
  # (gamma-bar is already direction-averaged, so finer resolution buys nothing).
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d) else NULL
  nb <- length(r_seq)
  cnt_b <- tabulate(findInterval(as.vector(dmat), r_seq, left.open = TRUE) + 1L,
                    nbins = nb)
  Kr <- if (is.null(gb)) {
    (area_vol * cumsum(cnt_b)) / (n1 * n2)
  } else {
    w <- ifelse(is.finite(gb) & gb > 0, area_vol / gb, 1)
    (area_vol * cumsum(cnt_b * w)) / (n1 * n2)
  }

  Lr <- if (use_3d) (3 * Kr / (4 * pi))^(1/3) else sqrt(Kr / pi)
  data.frame(r = r_seq, Kr = Kr, Lr = Lr, Lr_minus_r = Lr - r_seq)
}

## ── Pair Cross-Correlation Function g₁₂(r) ──────────────────────────────────
## g12(r) = (area × N_pairs_in_annulus) / (n1 × n2 × annulus_area)
## g12 = 1 under CSR; > 1 = attraction; < 1 = repulsion.
compute_pcf12 <- function(x1, y1, z1 = NULL, x2, y2, z2 = NULL,
                           r_seq = NULL, use_3d = FALSE, max_subsample = 2000L,
                           window_ranges = NULL, support = NULL,
                           edge_correct = TRUE) {
  if (use_3d && !is.null(z1) && !is.null(z2)) {
    pts1 <- cbind(x1, y1, z1); pts2 <- cbind(x2, y2, z2)
  } else {
    pts1 <- cbind(x1, y1);     pts2 <- cbind(x2, y2); use_3d <- FALSE
  }
  n1 <- nrow(pts1); n2 <- nrow(pts2)
  if (n1 < 2L || n2 < 2L) return(data.frame())

  pts1 <- pts1[.stable_subsample(n1, max_subsample), , drop = FALSE]
  pts2 <- pts2[.stable_subsample(n2, max_subsample, seed = 20260814L), , drop = FALSE]
  n1 <- nrow(pts1); n2 <- nrow(pts2)

  all_pts <- rbind(pts1, pts2)
  # per-axis extents; a pre-supplied window_ranges (2×d min/max matrix) makes the
  # observed curve and its CSR envelope share one identical window.
  rng     <- if (!is.null(window_ranges)) (window_ranges[2, ] - window_ranges[1, ])
             else apply(all_pts, 2, function(v) diff(range(v, na.rm = TRUE)))
  area_vol <- if (!is.null(support)) support_area(support)
              else if (use_3d) prod(rng) else prod(rng[1:2])
  if (!is.finite(area_vol) || area_vol <= 0) return(data.frame())

  dmat <- as.matrix(dist(rbind(pts1, pts2)))[seq_len(n1), (n1 + 1):(n1 + n2)]

  if (is.null(r_seq) || !length(r_seq)) {
    # Cap by the SHORTEST window axis (see compute_ripley_L / compute_ripley_L12).
    r_max <- min(if (use_3d) rng else rng[1:2]) / 4
    r_seq <- seq(0, r_max, length.out = 60L)
  }
  dr <- if (length(r_seq) > 1L) mean(diff(r_seq)) else r_seq[1L]

  # See compute_pcf: same translation correction, same single-pass shell binning.
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d) else NULL
  nb    <- length(r_seq)
  edges <- c(r_seq - dr / 2, r_seq[nb] + dr / 2)
  cnt_b <- tabulate(findInterval(as.vector(dmat), edges, left.open = TRUE), nbins = nb)

  g12 <- vapply(seq_along(r_seq), function(k) {
    r <- r_seq[k]
    lo <- max(0, r - dr / 2); hi <- r + dr / 2
    shell <- if (use_3d) (4 / 3) * pi * (hi^3 - lo^3) else pi * (hi^2 - lo^2)
    if (!is.finite(shell) || shell <= 0) return(NA_real_)
    denom <- if (is.null(gb)) area_vol else gb[k]
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    area_vol^2 * cnt_b[k] / (n1 * n2 * denom * shell)
  }, numeric(1L))

  data.frame(r = r_seq, g12 = g12)
}

## ── Univariate Pair-Correlation Function g(r) ────────────────────────────────
## GOAL: the one-channel counterpart of compute_pcf12 - the scale-LOCAL companion
## to Ripley's L(r). Counts ORDERED pairs of the SAME point set falling in the
## annulus/shell (r-dr/2, r+dr/2], normalised by the shell volume and by the
## intensity squared:
##
##     g(r) = A * N_pairs(annulus) / (n(n-1) * shell)      [pairs counted twice]
##
## g = 1 under CSR, > 1 clustering at scale r, < 1 exclusion at scale r.
##
## WHY THIS AND NOT ONLY L(r): K/L are CUMULATIVE, so a short-r label or blinking
## artefact is carried into every larger r and a real cluster shows up as a broad
## rising ramp whose peak position is not the cluster size. g(r) is the derivative
## form, so an aggregate of radius R gives a PEAK near R and a short-r artefact
## stays confined to the first bins. This is the PC-PALM statistic (Sengupta et
## al. 2011; Veatch et al. 2012).
##
## DENSITY NORMALISATION: g is divided by the intensity SQUARED, so - like K/L -
## it is invariant to how densely the sample is labelled, which is what makes it
## the endpoint of choice when comparing cells of differing expression level.
## That invariance is exact only for INDEPENDENT (Poisson) thinning; see the
## LABELLING-DENSITY block below.
##
## Self-pairs (i == j) are excluded via the n(n-1) denominator and by never
## counting the zero diagonal (dr/2 > 0, and the diagonal is dropped explicitly).
## No analytic edge correction - use compute_ripley_envelope(estimator = "g"),
## which carries the identical bias in the same window.
# =============================================================================
# fit_correlation_length() - ONE NUMBER FOR "HOW BIG ARE THE CLUSTERS"
# =============================================================================
#
# GOAL
#   Reduce a g(r) curve to a single comparable number, the correlation length xi,
#   by fitting the excess pair correlation
#
#       g(r) - 1  =  A * exp(-r / xi)
#
#   so that two conditions can be compared with a test on one value per cell
#   instead of by eye across two curve shapes.
#
# WHY AN EXPONENTIAL, AND WHAT xi MEANS
#   g(r) - 1 is the excess probability of finding a second molecule at
#   separation r. For clusters of a characteristic size it decays roughly
#   exponentially, and xi is the distance over which that excess falls by 1/e.
#   Read it as the scale over which molecules "know about" each other - close to
#   a cluster radius for compact clusters, longer for diffuse ones. A is the
#   amplitude at r -> 0, which is NOT a clean measure of clustering strength (see
#   the floor below), so it is reported but should not be the headline.
#
# HOW IT IS FITTED, AND WHY NOT THE OBVIOUS WAY
#   The obvious fit is log-linear: log(g - 1) = log(A) - r/xi is a straight line,
#   so ordinary least squares solves it and cannot fail to converge. It is also
#   BIASED, badly, and the bias is not obvious. Points with g <= 1 have no
#   logarithm and must be dropped - but in a noisy tail those are exactly the
#   downward noise excursions, so dropping them keeps only the upward ones, lifts
#   the tail, flattens the slope and OVERESTIMATES xi. Measured on synthetic
#   curves with a known answer (20 repeats, noise sd 0.05 on g):
#
#       true xi     40      60     120     250
#       log-linear 121      94     115     248      <- unusable below ~100 nm
#       nls         36      56     116     244
#
#   So the fit is nonlinear least squares on g itself, g = 1 + A*exp(-r/xi),
#   which uses every point including those below 1 and needs no logarithm. The
#   log-linear solution is still computed - it seeds nls, and it is the fallback
#   if nls fails to converge, with `method` recording which one produced the
#   number.
#
# THE FLOOR (r_min) IS THE IMPORTANT ARGUMENT
#   Short-r g(r) in SMLM is dominated by LABELLING, not biology: residual blinks
#   of one fluorophore and the several dyes of a multivalent label all pile up
#   within a few tens of nm. Fitting from r = 0 therefore measures the label,
#   not the structure. Start the fit above that floor - a sensible default is a
#   few times the measured trace spread, or the r where a non-clustering control
#   returns to 1.
#
#   Note the floor cuts both ways: if xi turns out to be SMALLER than r_min, the
#   fit began past the decay it was meant to measure and the estimate is not
#   trustworthy. That case is flagged rather than silently returned.
#
# ARGUMENTS
#   r, g   : the curve, as returned by compute_pcf()
#   r_min  : lower bound of the fit (nm). Points below it are ignored.
#   r_max  : upper bound (nm); NULL = use all of the curve above r_min.
#
# RETURNS
#   one-row data.frame: xi, xi_se, amplitude, r2, method, below_floor, n_used,
#   n_avail, r_min, r_max. xi is NA when fewer than 4 usable points remain or the
#   excess does not decay (no decay = no correlation length; reporting one there
#   would be fitting noise).
# =============================================================================
fit_correlation_length <- function(r, g, r_min = 0, r_max = NULL) {
  na_row <- data.frame(xi = NA_real_, xi_se = NA_real_, amplitude = NA_real_,
                       r2 = NA_real_, method = NA_character_, below_floor = FALSE,
                       n_used = 0L, n_avail = 0L,
                       r_min = r_min, r_max = r_max %||% NA_real_,
                       stringsAsFactors = FALSE)
  if (!length(r) || length(r) != length(g)) return(na_row)
  keep <- is.finite(r) & is.finite(g) & r >= r_min
  if (!is.null(r_max) && is.finite(r_max)) keep <- keep & r <= r_max
  na_row$n_avail <- sum(keep)
  if (na_row$n_avail < 4L) return(na_row)

  rr <- r[keep]; gg <- g[keep]
  na_row$n_used <- length(rr)
  na_row$r_max  <- max(rr, na.rm = TRUE)

  # Seed from the log-linear solution, on the excess points only.
  pos <- gg > 1
  if (sum(pos) < 4L) return(na_row)
  lf <- tryCatch(stats::lm(log(gg[pos] - 1) ~ rr[pos]), error = function(e) NULL)
  if (is.null(lf)) return(na_row)
  co <- stats::coef(lf); slope <- unname(co[2])
  if (!is.finite(slope) || slope >= 0) return(na_row)   # no decay to measure
  seed <- list(A = unname(exp(co[1])), xi = -1 / slope)

  fit <- tryCatch(stats::nls(gg ~ 1 + A * exp(-rr / xi), start = seed,
                             control = stats::nls.control(warnOnly = TRUE, maxiter = 200)),
                  error = function(e) NULL)

  if (!is.null(fit) && all(is.finite(stats::coef(fit))) &&
      is.finite(stats::coef(fit)[["xi"]]) && stats::coef(fit)[["xi"]] > 0) {
    cf <- stats::coef(fit)
    se <- tryCatch(summary(fit)$coefficients["xi", "Std. Error"], error = function(e) NA_real_)
    xi <- unname(cf[["xi"]]); amp <- unname(cf[["A"]])
    # Pseudo-R2 against the mean of g: nls has no r.squared, and 1 - RSS/TSS is
    # the comparable quantity.
    rss <- sum(stats::residuals(fit)^2); tss <- sum((gg - mean(gg))^2)
    r2  <- if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_
    meth <- "nls"
  } else {
    # Fallback: the seed itself. Recorded as such - it is the biased estimator,
    # and a reader should know which one produced the number.
    xi  <- seed$xi; amp <- seed$A
    se  <- tryCatch(summary(lf)$coefficients[2, 2] / slope^2, error = function(e) NA_real_)
    r2  <- tryCatch(summary(lf)$r.squared, error = function(e) NA_real_)
    meth <- "log-linear (nls did not converge)"
  }

  data.frame(xi = xi, xi_se = se, amplitude = amp, r2 = r2, method = meth,
             below_floor = is.finite(xi) && xi < r_min,
             n_used = length(rr), n_avail = na_row$n_avail,
             r_min = r_min, r_max = max(rr, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

compute_pcf <- function(x, y, z = NULL, r_seq = NULL, use_3d = FALSE,
                        max_subsample = 2000L, window_ranges = NULL,
                        support = NULL, edge_correct = TRUE) {
  pts <- if (use_3d && !is.null(z)) cbind(x, y, z) else { use_3d <- FALSE; cbind(x, y) }
  n <- nrow(pts)
  if (is.null(n) || n < 3L) return(data.frame())

  idx <- .stable_subsample(n, max_subsample)
  if (length(idx) < n) { pts <- pts[idx, , drop = FALSE]; n <- nrow(pts) }

  # Per-axis extents. A pre-supplied window_ranges (2×d min/max matrix) lets the
  # observed curve and its CSR envelope share one identical window.
  rng      <- if (!is.null(window_ranges)) window_ranges else apply(pts, 2, range)
  ext      <- rng[2, ] - rng[1, ]
  area_vol <- if (!is.null(support)) support_area(support, rng)
              else if (use_3d) prod(ext[1:3]) else prod(ext[1:2])
  if (!is.finite(area_vol) || area_vol <= 0) return(data.frame())

  if (is.null(r_seq) || !length(r_seq)) {
    # Cap by the SHORTEST window axis (matches compute_ripley_L / _L12 / _pcf12):
    # edge bias is set by the dimension that runs out first, in 3D almost always z.
    r_max <- min(if (use_3d) ext[1:3] else ext[1:2]) / 4
    r_seq <- seq(0, r_max, length.out = 60L)
  }
  dr <- if (length(r_seq) > 1L) mean(diff(r_seq)) else r_seq[1L]

  # Upper-triangle distances only: each unordered pair once, so the annulus count
  # is doubled below to get ordered pairs (matching the n(n-1) denominator).
  dm     <- as.matrix(stats::dist(pts))
  d_vals <- dm[lower.tri(dm)]

  # Translation edge correction from the support's set covariance (see
  # set_covariance_iso). Substituting gammabar = area_vol recovers the uncorrected
  # estimator, so `gb = NULL` is exactly the old behaviour.
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d) else NULL

  # SINGLE-PASS shell counts. The shells (r - dr/2, r + dr/2] tile because r_seq is
  # equally spaced, so one findInterval over the shared edges replaces one full sweep
  # of the pair distances per radius.
  nb    <- length(r_seq)
  edges <- c(r_seq - dr / 2, r_seq[nb] + dr / 2)
  cnt_b <- tabulate(findInterval(d_vals, edges, left.open = TRUE), nbins = nb)

  gr <- vapply(seq_along(r_seq), function(k) {
    r <- r_seq[k]
    lo <- max(0, r - dr / 2); hi <- r + dr / 2
    shell <- if (use_3d) (4 / 3) * pi * (hi^3 - lo^3) else pi * (hi^2 - lo^2)
    if (!is.finite(shell) || shell <= 0) return(NA_real_)
    denom <- if (is.null(gb)) area_vol else gb[k]
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    area_vol^2 * 2 * cnt_b[k] / (n * (n - 1) * denom * shell)
  }, numeric(1L))

  out <- data.frame(r = r_seq, g = gr)
  attr(out, "edge_corrected") <- !is.null(gb)
  out
}

## ═══════════════════════════════════════════════════════════════════════════════
## 8b. LABELLING-DENSITY NORMALISATION  (opt-in)
## ═══════════════════════════════════════════════════════════════════════════════
## WHY THIS BLOCK EXISTS
## Transient transfection gives each cell an undefined plasmid load, so labelling
## density varies cell to cell and run to run even within one condition. Almost
## every clustering read-out is density-dependent in the direction that MIMICS
## aggregation:
##   * k-NN distance scales as rho^(-1/d) - under pure CSR in 2D the mean 1-NN
##     distance is 0.5 * rho^(-1/2), so a high expresser "looks clustered".
##   * DBSCAN with fixed (eps, minPts) IS a density threshold: rho_crit =
##     minPts / (pi eps^2) in 2D. Raise rho and more clusters cross it, and they
##     chain into larger ones.
## So two transfections cannot be compared on those statistics as they stand.
##
## THE FIX AND WHY IT IS SOUND
## Independent ("Poisson") p-thinning - keep each unit with probability p,
## independently - sends the intensity lambda -> p*lambda and the second-order
## product density rho2 -> p^2 * rho2. The pair-correlation function
## g = rho2 / lambda^2 is therefore EXACTLY invariant, as is Ripley's K. Real
## clustering survives thinning; a pure density artefact does not. Thinning every
## file down to one common density therefore makes even the fixed-threshold
## DBSCAN read-outs comparable, without erasing the signal being measured.
##
## WHAT THE THINNING TARGET MUST BE - THE ONE REAL TRAP
## The target must be a density that is NOT itself inflated by clustering:
##   * GLOBAL density = n_units / support  is rearrangement-invariant. Clustering
##     moves points around inside the support but does not change n or A, so this
##     is the CORRECT thinning target.
##   * LOCAL NN-based intensity (nn_local_intensity) is CONFOUNDED with
##     clustering: aggregates shorten NN distances, inflating lambda_hat, so a
##     genuinely aggregating cell would be thinned HARDER and part of the signal
##     removed. It is offered as a DIAGNOSTIC and as a support-free cross-check,
##     and is flagged in the UI - it is not the recommended target.
## The high quantiles of the NN distribution are dominated by the sparse
## background rather than the aggregates, which is why nn_local_intensity takes a
## quantile range rather than the median (see `probs`).
## ─────────────────────────────────────────────────────────────────────────────

## One point per de-blinked unit - the point pattern every density statistic and
## every spatial statistic in this app should be measured on. A trace ID is ONE
## BLINK, so trace-ID density confounds expression level with photophysics and
## laser power and is NOT comparable between runs; prefer unit = "cluster"
## (de-blinked molecules) whenever de-blinking is enabled.
##   unit: "tid"     - one median position per trace ID
##         "cluster" - one median position per de-blinked molecule (noise dropped)
## Returns a data.frame(x, y, z) of finite centroids, or NULL if unavailable.
unit_centroids <- function(df, unit = "tid", use_3d = FALSE) {
  if (is.null(df) || !nrow(df) || !all(c("x", "y") %in% names(df))) return(NULL)
  if (!"z" %in% names(df)) df$z <- 0
  if (unit == "cluster") {
    if (!"cluster" %in% names(df) || !any(!is.na(df$cluster) & df$cluster > 0L)) return(NULL)
    pts <- df %>%
      dplyr::filter(!is.na(cluster), cluster > 0L) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop")
  } else {
    if (!"tid" %in% names(df)) return(NULL)
    pts <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop")
  }
  pts <- pts[is.finite(pts$x) & is.finite(pts$y) & (!use_3d | is.finite(pts$z)), , drop = FALSE]
  if (!nrow(pts)) return(NULL)
  as.data.frame(pts[, c("x", "y", "z"), drop = FALSE])
}

## GOAL: estimate the AREA (2D) or VOLUME (3D) actually occupied by the sample,
## so that n_units / support is a labelling density rather than an artefact of how
## much empty frame was scanned.
##
## WHY NOT THE BOUNDING BOX: MINFLUX scan fields routinely include the edge of a
## cell and interior dark regions (organelle exclusions, unlabelled compartments).
## The bounding box counts all of that as "sample", so a cell imaged half-in-frame
## reads at half its true density purely from framing. That is exactly the
## confound this whole block exists to remove.
##
## modes:
##   "grid"  (DEFAULT, edge- and dark-region-aware) - occupancy grid. Bin the unit
##           centroids into cells of side `h_nm`, then MORPHOLOGICALLY CLOSE the
##           occupied set by `close_nm` (dilate then erode) and count. The closing
##           is the knob that separates the two things that both look like holes:
##             - sparse SAMPLING inside the cell -> should be filled in (it is
##               sample), which closing does;
##             - a genuine DARK REGION or the outside of the cell -> should stay
##               excluded, which closing leaves alone as long as close_nm is
##               smaller than the dark region.
##           So close_nm should be a few times the typical unit spacing and well
##           below the size of the dark features you want excluded.
##   "hull"  - convex hull (geometry::convhulln). Follows the outer edge but
##           bridges concavities and dark regions, so it sits between bbox and grid.
##   "bbox"  - axis-aligned bounding box. The implicit choice everywhere else in
##           the app; kept for comparison and as the fallback.
##
## Returns list(support, mode, h, close_nm, n_occupied, bbox, frac_of_bbox, d).
## `support` is in nm^2 (2D) or nm^3 (3D). NULL only if there is nothing to measure.
##
## HONEST LIMITS: the grid estimate is h-dependent (it can only ever be a multiple
## of h^d) and, like any union-of-cells support, is biased upward at the boundary
## by roughly perimeter x h/2. Neither matters much for the intended use, which is
## a RATIO of densities between files measured with the SAME h and close_nm - the
## boundary bias largely cancels between comparable fields. It does matter if you
## quote an absolute molecules/um^2, so both the value and the settings are
## reported, and the fraction of the bounding box is shown so an implausible
## support is visible at a glance.
## `return_mask = TRUE` additionally returns the occupancy array itself (`mask`)
## plus `origin`/`pad`/`dims`/`win`, so the support can be DRAWN over the data -
## which is the only way to see whether the gap-fill radius did something sane.
##
## DARK REGIONS AND THE GAP FILL - the relationship worth knowing, because it means
## no extra control is needed. A morphological closing of radius r can only fill an
## interior hole narrower than about 2r, so the gap fill IS the dark-region control:
##     a dark region wider than ~2 * close_nm stays excluded.
## MEASURED on a 5 um field with a 1500 nm dark region and sparse sampling: at the
## auto plateau (200 nm) the dark region is correctly excluded and the support lands
## within -0.1% of truth. A separate "keep holes bigger than X" knob was built and
## REMOVED - re-excluding interior holes of the raw occupancy made the estimate
## WORSE (-9.0% at 400 nm, -2.0% at 800 nm vs -0.1% without it), because at the raw
## scale sparse-sampling voids are face-connected to genuine dark regions and get
## re-excluded with them. If dark regions are being swallowed, lower close_nm - the
## support map shows immediately whether they are.
estimate_support <- function(pts, use_3d = FALSE, mode = c("grid", "hull", "bbox", "kde"),
                             h_nm = 100, close_nm = 200, return_mask = FALSE,
                             bw_nm = 150, frac = 0.99) {
  mode <- match.arg(mode)
  if (is.null(pts) || !nrow(pts)) return(NULL)
  # KDE / highest-density region: a separate estimator, always 2D (the field is a
  # plan-view density; a 3D HDR would need a volumetric kernel and is not offered).
  if (mode == "kde") {
    k <- kde_hdr_support(pts, res_nm = h_nm, bw_nm = bw_nm, frac = frac,
                         return_mask = return_mask)
    if (is.null(k)) return(NULL)
    rng2 <- apply(as.matrix(pts[, c("x", "y"), drop = FALSE]), 2, range)
    bb2  <- prod(pmax(rng2[2, ] - rng2[1, ], h_nm))
    k$bbox <- bb2
    k$frac_of_bbox <- k$support / bb2
    k$d <- 2L
    # Normalise the returned shape to match the grid/hull/bbox branches, so every
    # caller (density row, map renderer, status readout) can read the same fields.
    k$h        <- k$res
    k$close_nm <- NA_real_
    return(k)
  }
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < 1L) return(NULL)
  d <- ncol(m)

  rng  <- apply(m, 2, range)
  ext  <- rng[2, ] - rng[1, ]
  # A degenerate axis (all points equal) would zero the product; floor it at h so
  # a genuinely flat slab still reports a finite support instead of 0.
  ext_f <- pmax(ext, h_nm)
  bbox  <- prod(ext_f)

  out <- list(mode = mode, h = h_nm, close_nm = close_nm, d = d,
              bbox = bbox, n_occupied = NA_integer_)

  if (mode == "bbox" || nrow(m) < (d + 2L)) {
    out$support <- bbox
    if (mode != "bbox") out$mode <- "bbox (too few units for the chosen mode)"
    out$frac_of_bbox <- 1
    return(out)
  }

  if (mode == "hull") {
    hv <- tryCatch({
      h <- geometry::convhulln(m, options = "FA")
      if (d == 2L) h$vol else h$vol   # convhulln$vol is area in 2D, volume in 3D
    }, error = function(e) NA_real_)
    out$support      <- if (is.finite(hv) && hv > 0) hv else bbox
    out$frac_of_bbox <- out$support / bbox
    return(out)
  }

  # ── occupancy grid + morphological closing ────────────────────────────────
  h <- max(1e-6, as.numeric(h_nm))
  # Integer voxel index per point, 1-based, per axis.
  idx <- sweep(m, 2, rng[1, ], "-") / h
  idx <- floor(idx) + 1L
  dims <- as.integer(pmax(1L, apply(idx, 2, max)))
  # Pad by the closing radius in voxels so a dilation near the border has room.
  rad  <- max(0L, as.integer(floor(close_nm / h)))
  pad  <- rad + 1L
  dimp <- dims + 2L * pad
  idxp <- matrix(as.integer(sweep(idx, 2, rep(pad, d), "+")), ncol = d)

  occ <- array(FALSE, dim = dimp)
  occ[idxp] <- TRUE                      # matrix-index assignment: one row = one cell
  n_raw <- sum(occ)

  if (rad > 0L) {
    # Separable Chebyshev (box) dilation then erosion = morphological closing with
    # a cubic structuring element. Separable => O(rad) passes per axis, not O(rad^d).
    box_max <- function(a, axis, r) {
      res <- a
      for (s in seq_len(r)) {
        res <- res | .array_shift(res, axis, 1L) | .array_shift(res, axis, -1L)
      }
      res
    }
    dil <- occ
    for (ax in seq_len(d)) dil <- box_max(dil, ax, rad)
    # Erosion of X = complement of the dilation of the complement.
    ero <- !dil
    for (ax in seq_len(d)) ero <- box_max(ero, ax, rad)
    closed <- !ero
  } else {
    closed <- occ
  }

  n_occ <- sum(closed)
  out$n_occupied   <- as.integer(n_occ)
  out$n_raw_cells  <- as.integer(n_raw)
  out$support      <- n_occ * h^d
  out$frac_of_bbox <- out$support / bbox
  if (isTRUE(return_mask)) {
    # The mask is PADDED by `pad` cells on every side so the dilation had room, so
    # padded cell index i on an axis starts at world coordinate
    # rng[1, axis] + (i - 1 - pad) * h. `origin` is that world coordinate of the
    # first padded cell, which is what a renderer needs to place the grid.
    out$mask   <- closed
    out$origin <- rng[1, ] - pad * h
    out$pad    <- pad
    out$dims   <- dimp
    out$win    <- rng
  }
  out
}

# Shift a logical array by one cell along `axis` (vacated cells become FALSE).
# Helper for the separable morphological closing in estimate_support; internal.
# Uses `TRUE` as a whole-axis index so the same code works for 2D and 3D.
.array_shift <- function(a, axis, by) {
  dm <- dim(a)
  n  <- dm[axis]
  if (by == 0L) return(a)
  out <- array(FALSE, dim = dm)
  if (n < 2L) return(out)
  src <- if (by > 0L) seq_len(n - 1L) else 2:n
  dst <- if (by > 0L) 2:n             else seq_len(n - 1L)
  i_s <- rep(list(TRUE), length(dm)); i_s[[axis]] <- src
  i_d <- rep(list(TRUE), length(dm)); i_d[[axis]] <- dst
  slice <- do.call("[", c(list(a), i_s, list(drop = FALSE)))
  do.call("[<-", c(list(out), i_d, list(value = slice)))
}

## GOAL: support from a SMOOTHED DENSITY FIELD instead of a binary occupancy grid -
## a kernel density estimate of the localizations, thresholded at the level that
## encloses a given FRACTION of the units (a highest-density region, Hyndman 1996).
##
## WHY THIS FORM AND NOT AN ABSOLUTE DENSITY LEVEL. Thresholding at a fixed
## density level is the obvious version and it is CIRCULAR: the support becomes a
## function of the density it is the denominator for, so a sparsely labelled cell
## falls below the level over most of its area and its measured density is dragged
## toward the threshold. Measured on two scenes with identical geometry and a known
## 4.0x labelling difference, the recovered ratio was **0.77 and 0.55** - not merely
## compressed but INVERTED - and the support of the sparse scene collapsed to zero
## by n = 1000. The HDR form asks instead "where do 99% of the units live", which is
## invariant to how many there are: same scenes, recovered ratio **3.90 and 3.90**
## (occupancy grid: 3.85 / 3.78).
##
## Measured support error vs known truth at n = 3000, bandwidth 150 nm:
##   plain square   grid -0%, KDE-HDR -0.6%
##   L-shape        grid -0%, KDE-HDR +2.3%
##   square + dark  grid -0%, KDE-HDR -0.6%
## Both degrade at very low n (n = 500: grid -15%, KDE-HDR -18%). KDE-HDR holds its
## support better under thinning (to 25% of units: -3.7% vs the grid's -6.2%; to 10%:
## -18.3% vs -26.3%), which matters because thinning is a core feature here.
##
## WHY IT IS THE DEFAULT - measured on CLUSTERED data (60 structures in a 5 um cell,
## 4.0x labelling difference), recovered density ratio:
##   KDE-HDR  bw 150 -> 4.03,  250 -> 3.95,  400 -> 4.03,  600 -> 4.02
##   grid     fill 200 -> 3.20, 400 -> 3.53, 800 -> 3.83, 1200 -> 3.83
## The ratio is nearly INDEPENDENT of bandwidth, and right, even where the support is
## only 30% of the cell - because the level is calibrated per file to enclose the same
## FRACTION of units, so the two supports agree with each other (7.46e6 vs 7.51e6 nm2)
## whether or not either equals the cell. The grid's threshold is an absolute count
## (">=1 unit per cell"), so a sparsely labelled file gets a systematically smaller
## mask (8.05e6 vs 6.45e6, 25% apart) and the ratio is dragged down.
## Consequence for use: bandwidth sets the ABSOLUTE area (over-smoothing spills past
## the cell - 107% of it at bw 600, 123% at 1400) but barely affects the RATIO, which
## is what a between-condition comparison depends on.
##
## `frac` is the fraction of UNITS enclosed; `bw_nm` MUST be held common across the
## samples being compared (it is the smoothing, and a per-file bandwidth would make
## each file's gaps close differently). Dark regions are excluded whenever they are
## wider than roughly 2 * 3 * bw_nm, since no kernel reaches the middle.
##
## Scale-invariance note: the level is a quantile of the SAME field it thresholds, so
## multiplying the field (e.g. more localizations per unit) scales the level equally
## and the mask does not move. That is what makes it thinning-invariant.
##
## The field is built by BINNING then convolving with a separable Gaussian, so cost
## is set by the grid, not by the number of localizations.
kde_hdr_support <- function(pts, res_nm = 50, bw_nm = 150, frac = 0.99,
                           return_mask = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  if (!all(c("x", "y") %in% names(pts))) return(NULL)
  p <- pts[stats::complete.cases(pts[, c("x", "y")]), , drop = FALSE]
  n <- nrow(p)
  if (n < 3L) return(NULL)
  res <- max(1e-6, as.numeric(res_nm)); bw <- max(res / 2, as.numeric(bw_nm))
  pad <- 3 * bw

  xr <- range(p$x); yr <- range(p$y)
  gx <- seq(xr[1] - pad, xr[2] + pad, by = res)
  gy <- seq(yr[1] - pad, yr[2] + pad, by = res)
  nx <- length(gx); ny <- length(gy)
  if (nx < 3L || ny < 3L) return(NULL)

  ix <- pmin(pmax(findInterval(p$x, gx), 1L), nx)
  iy <- pmin(pmax(findInterval(p$y, gy), 1L), ny)
  cnt <- matrix(0, nx, ny)
  tb  <- table(ix, iy)
  cnt[cbind(as.integer(rownames(tb)[row(tb)]), as.integer(colnames(tb)[col(tb)]))] <-
    as.vector(tb)

  # Separable Gaussian smoothing via banded kernel matrices (exact, and independent
  # of n). Kernels are normalised so the field is units per nm^2.
  kern <- function(m) {
    d <- outer(seq_len(m), seq_len(m), function(a, b) (a - b) * res)
    exp(-0.5 * (d / bw)^2)
  }
  z <- kern(nx) %*% cnt %*% kern(ny)
  z <- z / (2 * pi * bw^2)

  # Level = the (1 - frac) quantile of the field AT the unit positions.
  vals <- z[cbind(ix, iy)]
  lvl  <- as.numeric(stats::quantile(vals, 1 - frac, na.rm = TRUE))
  mask <- z >= lvl
  out <- list(support = sum(mask) * res^2, level = lvl, bw = bw, res = res,
              frac = frac, mode = "kde")
  if (isTRUE(return_mask)) {
    out$mask   <- mask
    out$origin <- c(gx[1], gy[1])
    out$win    <- rbind(c(xr[1], yr[1]), c(xr[2], yr[2]))
    out$h      <- res
  }
  out
}

## ═══════════════════════════════════════════════════════════════════════════════
## SUPPORT-AWARE OBSERVATION WINDOW for the spatial statistics
## ═══════════════════════════════════════════════════════════════════════════════
## WHY: Ripley's K/L, the pair correlation g and the co-localization null all need
## an observation window, and all of them used the axis-aligned BOUNDING BOX. For a
## cell that does not fill its box - a concave shape, an interior dark region, a
## scan field larger than the cell - the box overstates the area, so
##     lambda_hat = n / A  is UNDERSTATED,
## and since K = A * (pair fraction), an understated lambda inflates K at every r.
## The pattern then looks clustered because of the framing. The same error makes a
## bounding-box CSR null anti-conservative: the 2026-07-30 audit measured 53%
## (random-shift) and 93% (CSR) false-positive rates for co-localization when both
## channels shared one cell inside a larger window, against a nominal 5%, and
## identified mask-aware randomization as the outstanding fix. This is it.
##
## A `support` here is the list returned by estimate_support(..., return_mask = TRUE):
## $mask (logical array), $origin (world coords of the first cell's low corner),
## $h (cell side) and $support (the measure). `support_area()` reads the measure and
## `sample_in_support()` draws uniformly on it, so the observed statistic and its
## envelope share one geometry.
##
## MEASURED (scratchpad/test_support_window.R). A CSR pattern generated INSIDE an
## L-shaped cell whose support is 63% of its bounding box - truth is "not clustered":
##   Ripley  mean L(r)-r    bbox window +68.9 nm  ->  support window -13.5 nm
##   g(r)    mean g, r>150  bbox window 1.422     ->  support window 0.871   (want 1.0)
## Co-localization false-positive rate, two channels placed INDEPENDENTLY in one cell
## inside a larger window (truth = not co-localized, nominal 5%):
##   random-shift null   bbox 17%  ->  support 0%
##   CSR null            bbox 37%  ->  support 0%
## with power retained (10/10 genuinely co-localized pairs detected, median p 0.025).
##
## ⚠ RESIDUAL, do not oversell: the correction is only as good as the support
## estimate. Here the estimated support came out 6% larger than the true L-shape, and
## that alone leaves the observed curve outside a narrow envelope most of the time
## even though the bias in L-r is down ~80%. And 0/30 rejections cannot distinguish
## "calibrated" from "conservative" - the honest statement is that the support-aware
## null is CONSERVATIVE rather than anti-conservative, which is the safe direction.
support_area <- function(support, fallback_ranges = NULL) {
  if (!is.null(support) && is.finite(support$support %||% NA_real_) && support$support > 0)
    return(as.numeric(support$support))
  if (!is.null(fallback_ranges)) return(prod(fallback_ranges[2, ] - fallback_ranges[1, ]))
  NA_real_
}

## Uniform sample on the support: pick occupied cells with replacement, then a
## uniform position inside each. Exact uniform-on-support up to the cell size, and
## it works for any shape including one with holes - which a rejection sampler on
## the bounding box would also do but far more slowly for a thin or sparse support.
## Returns a list(x, y[, z]) of length n, or NULL if the support has no cells.
sample_in_support <- function(support, n, use_3d = FALSE) {
  if (is.null(support) || is.null(support$mask)) return(NULL)
  m <- support$mask
  d <- length(dim(m))
  occ <- which(m)
  if (!length(occ)) return(NULL)
  h  <- support$h
  o  <- support$origin
  pick <- occ[sample.int(length(occ), n, replace = TRUE)]
  ai <- arrayInd(pick, dim(m))
  # Cell (i,..) spans [origin + (i-1)h, origin + i*h] on each axis.
  jitter_axis <- function(k) o[k] + (ai[, k] - 1L) * h + stats::runif(n, 0, h)
  out <- list(x = jitter_axis(1L), y = jitter_axis(2L))
  out$z <- if (use_3d && d >= 3L) jitter_axis(3L) else NULL
  out
}

## Is a set of points inside the support? Used to restrict a randomized pattern to
## the cell (the co-localization shift null) rather than to the box.
in_support <- function(support, x, y, z = NULL) {
  if (is.null(support) || is.null(support$mask)) return(rep(TRUE, length(x)))
  m <- support$mask; dm <- dim(m); d <- length(dm); h <- support$h; o <- support$origin
  idx <- function(v, k) pmin(pmax(floor((v - o[k]) / h) + 1L, 1L), dm[k])
  ix <- idx(x, 1L); iy <- idx(y, 2L)
  if (d >= 3L && !is.null(z)) m[cbind(ix, iy, idx(z, 3L))] else {
    if (d >= 3L) {
      # 3D mask but 2D query: a column counts if any slice in it is occupied.
      p2 <- apply(m, c(1L, 2L), any); p2[cbind(ix, iy)]
    } else m[cbind(ix, iy)]
  }
}

## GOAL: compute the window's set covariance ONCE and cache it on the support object.
##
## WHY THIS EXISTS: gamma-bar depends only on the WINDOW, never on the point pattern,
## but the estimators are called once per Monte-Carlo simulation. Recomputing an FFT
## plus a million-shift radial binning inside a 99-simulation envelope loop made
## Ripley's L appear to hang - measured 0.78 s per simulation on a 210x210 mask, so
## ~77 s for one file, with the progress bar only advancing per FILE and therefore
## sitting at 0 the whole time. Attach once, reuse for the observed curve and every
## simulation. `gcache_r` is the grid gamma-bar was evaluated on; the estimators
## interpolate off it rather than rebuilding.
attach_gamma_cache <- function(support, r_max, use_3d = FALSE, n_grid = 256L) {
  if (is.null(support) || is.null(support$mask)) return(support)
  if (!is.null(support$gcache) && isTRUE(support$gcache_max >= r_max)) return(support)
  rg <- seq(0, max(r_max, support$h), length.out = as.integer(n_grid))
  gb <- tryCatch(set_covariance_iso(
    if (use_3d) support$mask else mask_project_2d(support$mask), support$h, rg),
    error = function(e) NULL)
  if (is.null(gb)) return(support)
  support$gcache     <- gb
  support$gcache_r   <- rg
  support$gcache_max <- max(rg)
  support
}

# gamma-bar at arbitrary distances, from the cache when present, else computed.
.gamma_at <- function(support, d, use_3d = FALSE) {
  if (is.null(support) || is.null(support$mask) || !length(d)) return(NULL)
  if (!is.null(support$gcache))
    return(stats::approx(support$gcache_r, support$gcache, xout = d, rule = 2)$y)
  rg <- seq(0, max(d, na.rm = TRUE), length.out = 256L)
  gb <- set_covariance_iso(
    if (use_3d) support$mask else mask_project_2d(support$mask), support$h, rg)
  if (is.null(gb)) return(NULL)
  stats::approx(rg, gb, xout = d, rule = 2)$y
}

## GOAL: a rectangular "support" over the bounding box, so the edge correction can be
## applied WITHOUT any support estimate.
##
## This matters because the edge correction and the support window fix two different
## things and only one of them is usually needed. For a cell that FILLS its scan
## field the bounding box already IS the correct window - measured on clustered data
## in a 5 um cell, bbox gave peak L-r +196 against the true window's +202 - so all
## that is missing there is the boundary correction. Requiring a support estimate to
## get it would have made the safe, always-valid fix conditional on the risky one.
bbox_support <- function(pts, use_3d = FALSE, h_nm = 100) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < 3L) return(NULL)
  rng <- apply(m, 2, range)
  ext <- pmax(rng[2, ] - rng[1, ], h_nm)
  dims <- pmax(2L, as.integer(ceiling(ext / h_nm)))
  mask <- array(TRUE, dim = dims)
  list(mask = mask, origin = rng[1, ], h = h_nm,
       support = prod(dims) * h_nm^length(dims),
       win = rng, mode = "bbox", d = length(dims))
}

## GOAL: how many face-connected pieces the window mask has - the diagnostic that
## separates "a cell" from "a pile of structures".
##
## WHY IT IS NEEDED. A density-derived support answers "where is the sample dense",
## which on clustered data is the CLUSTERS, not the cell. Using that as the
## observation window conditions on the very thing a clustering statistic tests -
## inside the clusters the points look locally uniform, so the signal collapses.
## Measured on clustered data in a 5 um cell: the KDE/HDR support came to 29% of the
## cell and crushed peak L(r)-r from +202 to +49 while dragging the tail to -34,
## which is the "straight line into the negative" symptom. A footprint estimator
## (occupancy grid at its closing plateau) gave 79% and peak +163.
## A genuine cell footprint is ONE piece; a structures-hugging mask is many, and that
## is checkable without knowing the true cell.
## Returns list(n, largest_frac) - the piece count AND the share of the window that
## lives in its biggest piece.
##
## ⚠ USE `largest_frac`, NOT `n`. A raw count is the wrong criterion and rejected
## perfectly good real data: measured on four real MINFLUX acquisitions the footprint
## came out in 2, 3, 1 and 1 pieces, because sparse labelling at a cell edge leaves a
## few detached slivers even though 90%+ of the window is one cell. A structures-
## hugging mask is different in kind - many pieces of COMPARABLE size, so no single
## piece dominates. `largest_frac` separates those two cases; a count cannot.
mask_components <- function(mask, h) {
  if (is.null(mask) || !any(mask)) return(list(n = 0L, largest_frac = 0))
  idx <- which(mask)
  if (length(idx) == 1L) return(list(n = 1L, largest_frac = 1))
  sub <- arrayInd(idx, dim(mask))
  cl <- tryCatch(dbscan::dbscan(sub * h, eps = 1.05 * h, minPts = 1L)$cluster,
                 error = function(e) NULL)
  if (is.null(cl)) return(list(n = NA_integer_, largest_frac = NA_real_))
  tb <- table(cl[cl > 0])
  if (!length(tb)) return(list(n = 0L, largest_frac = 0))
  list(n = length(tb), largest_frac = max(tb) / sum(tb))
}

## GOAL: the isotropised SET COVARIANCE of the support, which is what turns Ripley's
## K and the pair correlation g into unbiased estimators in a window of ANY shape.
##
## WHY THIS IS THE REAL FIX, AND WHY THE AREA WAS NOT. Both estimators divide by an
## area, so a wrong area is a MULTIPLICATIVE scale error. But even with the exactly
## correct area they are still biased, because a point near the boundary has fewer
## neighbours available at separation r than one in the middle - nothing about the
## area fixes that. Measured on CSR (truth g = 1) inside an L-shaped cell:
##     bounding-box area, no correction   1.398
##     support area, no correction        0.857
##     EXACT true area, no correction     0.827   <- the area was never the problem
##     exact area + this correction       1.003
##     support area + this correction     1.030
## and on a plain square with the exact area, uncorrected g decays 0.949 -> 0.783
## across r (mean 0.882), corrected 0.984. The bounding box only ever looked closer
## to 1 because its area error and the edge deficit had opposite signs and partly
## cancelled; enabling the support removed one and exposed the other, which is why
## the curves started looking wrong.
##
## THE MATHS. For a window W, gamma(h) = |W ∩ W_h| is the set covariance. The number
## of ordered pairs with separation in a shell S has expectation
##     lambda^2 * integral_S g(h) gamma(h) dh  ~=  lambda^2 g(r) gammabar(r) |S|,
## so the estimator must divide by gammabar(r), not by |W|. Substituting gamma = |W|
## recovers the uncorrected form, i.e. the uncorrected estimator is the special case
## that assumes NO area is lost to the boundary. This is the translation (Ohser)
## correction; see Illian et al. 2008 §4.3.3, Baddeley et al. 2015 §7.4.
##
## COMPUTATION. gamma is the autocorrelation of the mask indicator, so one FFT gives
## it for every shift at once. The array is ZERO-PADDED to twice its size per axis -
## without that the FFT autocorrelation is circular and wraps opposite edges
## together. gamma(0) must come back equal to the support measure; that is asserted
## by the callers' tests.
##
## RADIAL BINNING: shifts live on a lattice of spacing h, so a shell thinner than h
## can contain NO lattice point. Binning at the requested r-spacing therefore leaves
## most shells empty (this bit me in the prototype: NAs at scattered radii, which
## silently dropped those radii from the curve). gammabar is therefore accumulated in
## bins of width h and then INTERPOLATED onto the requested radii.
set_covariance_iso <- function(mask, h, r_out) {
  if (is.null(mask) || !any(mask) || !length(r_out)) return(NULL)
  dm <- dim(mask); nd <- length(dm)
  if (nd < 2L || nd > 3L) return(NULL)
  pad <- vapply(dm, function(k) 2^ceiling(log2(2 * k)), numeric(1))
  # Guard against an absurd FFT (a very fine 3D grid); fall back to no correction.
  if (prod(pad) > 3e7) return(NULL)
  m <- array(0, dim = pad)
  m[slice.index(m, 1) <= dm[1]] <- 0   # (initialised zero; explicit for clarity)
  idx <- lapply(dm, seq_len)
  m <- do.call("[<-", c(list(m), idx, list(value = as.numeric(mask))))

  F  <- stats::fft(m)
  Cc <- Re(stats::fft(F * Conj(F), inverse = TRUE)) / prod(pad)   # overlap in CELLS

  # Shift lattice per axis: index 1 is shift 0, upper half is negative shifts.
  shifts <- lapply(pad, function(k) c(0:(k / 2 - 1), -(k / 2):(-1)) * h)
  d2 <- if (nd == 2L) outer(shifts[[1]]^2, shifts[[2]]^2, "+") else
    outer(outer(shifts[[1]]^2, shifts[[2]]^2, "+"), shifts[[3]]^2, "+")
  dist <- sqrt(d2)

  # Bin at the lattice spacing, then interpolate - see RADIAL BINNING above.
  rmax <- max(r_out) + h
  brk  <- seq(0, rmax + h, by = h)
  nb   <- length(brk) - 1L
  bin  <- findInterval(as.vector(dist), brk, rightmost.closed = TRUE)
  keep <- bin >= 1L & bin <= nb
  cv   <- as.vector(Cc)[keep]; bin <- bin[keep]
  if (!length(bin)) return(NULL)
  # rowsum is C-level; tapply over ~1e6 lattice shifts was a measurable cost.
  rs_  <- rowsum(cv, bin, reorder = TRUE)
  cnt  <- tabulate(bin, nbins = nb)
  ok   <- as.integer(rownames(rs_))
  sums <- as.numeric(rs_[, 1]) / cnt[ok]
  ctr  <- (brk[-1] + brk[-length(brk)]) / 2
  if (length(ok) < 2L) return(NULL)
  g_at <- stats::approx(ctr[ok], as.numeric(sums), xout = r_out, rule = 2)$y
  # (binning above uses rowsum, which is C-level; tapply over ~1e6 lattice shifts
  #  was a measurable share of the cost when this ran inside an envelope loop.)
  cell <- h^nd
  gb   <- g_at * cell                                   # cells -> nm^2 / nm^3
  # gamma is non-increasing in |h| for a convex window and bounded by |W| always;
  # clamp tiny numerical excursions so the correction can never shrink the estimate.
  gb <- pmin(pmax(gb, .Machine$double.eps), sum(mask) * cell)
  gb
}

## GOAL: turn an occupancy mask into drawable OUTLINE segments, so the support can
## be shown as a boundary over the localizations rather than as a number.
##
## An edge is on the boundary when an occupied cell's neighbour is not occupied (or
## is off the grid). Emitting exactly those edges outlines the support including any
## interior holes, which is what makes a swallowed dark region visible at a glance.
## Segments are returned as one x/y pair with NA separators, so the whole outline is
## a SINGLE plotly trace however many disjoint pieces it has.
##
## `mask` is 2D (project a 3D mask first - see mask_project_2d). `origin` is the
## world coordinate of the first cell's low corner, `h` the cell side.
mask_outline_segments <- function(mask, origin, h) {
  if (is.null(mask) || !any(mask)) return(list(x = numeric(0), y = numeric(0)))
  d <- dim(mask)
  occ <- which(mask, arr.ind = TRUE)          # columns: row = x index, col = y index
  ix <- occ[, 1L]; iy <- occ[, 2L]
  # Cell (i,j) spans [x0 + (i-1)h, x0 + i*h] x [y0 + (j-1)h, y0 + j*h].
  x_lo <- origin[1] + (ix - 1L) * h; x_hi <- x_lo + h
  y_lo <- origin[2] + (iy - 1L) * h; y_hi <- y_lo + h

  # Is the neighbour occupied? Off-grid counts as not occupied.
  nb <- function(dx, dy) {
    jx <- ix + dx; jy <- iy + dy
    ok <- jx >= 1L & jx <= d[1] & jy >= 1L & jy <= d[2]
    out <- rep(FALSE, length(ix))
    if (any(ok)) out[ok] <- mask[cbind(jx[ok], jy[ok])]
    out
  }
  xs <- ys <- numeric(0)
  push <- function(x1, y1, x2, y2, sel) {
    if (!any(sel)) return(invisible(NULL))
    xs <<- c(xs, as.numeric(rbind(x1[sel], x2[sel], NA)))
    ys <<- c(ys, as.numeric(rbind(y1[sel], y2[sel], NA)))
  }
  push(x_lo, y_lo, x_hi, y_lo, !nb(0L, -1L))   # bottom edge
  push(x_lo, y_hi, x_hi, y_hi, !nb(0L,  1L))   # top edge
  push(x_lo, y_lo, x_lo, y_hi, !nb(-1L, 0L))   # left edge
  push(x_hi, y_lo, x_hi, y_hi, !nb( 1L, 0L))   # right edge
  list(x = xs, y = ys)
}

# Collapse a 3D occupancy mask onto XY: a column counts as sample if ANY z slice in
# it is occupied. The 2D map is a projection, so a 3D support looks larger in plan
# view than it is - stated in the UI rather than silently.
mask_project_2d <- function(mask) {
  if (is.null(mask)) return(NULL)
  if (length(dim(mask)) < 3L) return(mask)
  apply(mask, c(1L, 2L), any)
}

## GOAL: make the gap-fill (closing) radius CHOOSABLE instead of guessed, by
## sweeping it and reporting the resulting support.
##
## The curve has a characteristic shape and the shape is the whole point:
##   * at close_nm = 0 the support is just the occupied cells, so it tracks the
##     STRUCTURES rather than the cell - strongly n-dependent and not a density
##     denominator you want;
##   * it then RISES steeply as sampling gaps inside the sample get filled;
##   * it PLATEAUS once the sample region is solid - this plateau is the footprint,
##     and its left edge is the radius to use;
##   * it rises again only if the radius grows big enough to bridge genuine dark
##     regions or reach across the frame, which is what you must NOT do.
## Pick the radius at the start of the plateau. A support that is still climbing
## steeply at your chosen radius means the estimate is not yet a footprint; one
## that has climbed to near the bounding box means real dark regions are being
## swallowed.
##
## Returns data.frame(close_nm, support, frac_of_bbox, rel_step) where rel_step is
## the fractional increase from the previous radius - the plateau is where
## rel_step falls near 0. `plateau_close_nm` reads the plateau start off this.
##
## MEASURED, do not re-derive: on two synthetic cells with identical cluster
## geometry and a 4.0x labelling-density difference, a 200 nm radius (the old
## default) gave supports of 30% and 24% of frame and a density ratio of 3.13,
## while the plateau at 600-800 nm gave 71% and 70% - agreeing to 1 percentage
## point - and a ratio of 3.94. Choosing this radius badly is the single largest
## error source in the whole density normalisation, which is why it is swept.
support_closing_sweep <- function(pts, use_3d = FALSE, h_nm = 100,
                                  close_seq = NULL) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  if (is.null(close_seq)) close_seq <- h_nm * c(0, 1, 2, 3, 4, 6, 8, 12, 16, 24)
  close_seq <- sort(unique(pmax(0, as.numeric(close_seq))))
  out <- lapply(close_seq, function(cn) {
    s <- tryCatch(estimate_support(pts, use_3d = use_3d, mode = "grid",
                                   h_nm = h_nm, close_nm = cn),
                  error = function(e) NULL)
    if (is.null(s)) return(NULL)
    data.frame(close_nm = cn, support = s$support, frac_of_bbox = s$frac_of_bbox)
  })
  res <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(res) || !nrow(res)) return(NULL)
  res$rel_step <- c(NA_real_, diff(res$support) / utils::head(res$support, -1))
  res
}

## GOAL: read the plateau start off a closing sweep - the smallest radius at which
## the support has stopped growing appreciably, i.e. sampling gaps are filled but
## dark regions are not yet bridged.
##
## `tol` is the fractional support increase below which the curve counts as flat.
## Guards: never returns a radius whose support already fills more than
## `max_frac` of the bounding box (that is the "swallowed the dark regions" end of
## the curve), and falls back to the largest swept radius below that cap.
plateau_close_nm <- function(sweep, tol = 0.05, max_frac = 0.97) {
  if (is.null(sweep) || !nrow(sweep)) return(NA_real_)
  s <- sweep[order(sweep$close_nm), , drop = FALSE]
  ok <- is.finite(s$frac_of_bbox) & s$frac_of_bbox <= max_frac
  if (!any(ok)) return(NA_real_)
  s <- s[ok, , drop = FALSE]
  hit <- which(is.finite(s$rel_step) & s$rel_step < tol)
  if (length(hit)) return(s$close_nm[hit[1]])
  s$close_nm[nrow(s)]
}

## GOAL: a SUPPORT-FREE intensity estimate, from the k-th nearest-neighbour
## distances of the units themselves.
##
## For a Poisson process of intensity lambda in R^d, the k-th NN distance R_k
## satisfies lambda * omega_d * R_k^d ~ Gamma(k, 1), where omega_d is the volume
## of the unit ball (omega_2 = pi, omega_3 = 4pi/3). Hence
##     E[R_k] = Gamma(k + 1/d) / (Gamma(k) * (lambda * omega_d)^(1/d))
## which inverts to the moment estimator
##     lambda_hat = ( Gamma(k + 1/d) / (Gamma(k) * mean(R_k)) )^d / omega_d.
##
## WHY IT IS USEFUL: it needs no window and no support estimate, and because it
## only measures where points actually ARE, it is automatically blind to the frame
## edge, to the cell boundary and to interior dark regions. That makes it the
## natural cross-check on estimate_support: if the grid support is sane, the two
## density numbers agree to within a factor of order one.
##
## ⚠ WHY IT IS NOT THE THINNING TARGET: it is a LOCAL intensity, so clustering
## inflates it. An aggregating cell has short NN distances and reads as denser
## than it is, and thinning to match it would remove part of the aggregation
## signal. `probs` therefore defaults to the UPPER-MIDDLE of the NN distribution
## (0.5-0.9), which is dominated by the sparse background rather than by the
## aggregates - a partial, not a complete, mitigation. Use the global density for
## thinning; use this to diagnose.
##
## CALIBRATION OF THE QUANTILE BAND (do not remove - it is worth 30%):
## the Gamma relation above gives E[R_k] over ALL points. Averaging R_k over a
## quantile band instead raises the mean, so feeding a band mean into the
## full-distribution formula under-estimates lambda by ~30% at probs = (0.5, 0.9).
## Measured before the fix: ratio to truth 0.64 / 0.70 / 0.74 across three decades
## of lambda in 2D and 0.63 in 3D - a systematic bias, not noise. The band is
## therefore calibrated against its OWN conditional expectation,
##     c(k, d, a, b) = E[ G^(1/d) | G in (qgamma(a, k), qgamma(b, k)) ],  G ~ Gamma(k, 1)
## evaluated numerically, so lambda_hat = (c / mean(R_k | band))^d / omega_d is
## unbiased on CSR for any band. Verified back to ratio ~1.0 (see test 3).
##
## Returns list(lambda, mean_r, median_r, k, d, n) with lambda in units per nm^d.
nn_local_intensity <- function(pts, k = 3L, use_3d = FALSE, probs = c(0.5, 0.9)) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m); d <- ncol(m)
  if (n < 3L) return(NULL)
  k_eff <- min(max(1L, as.integer(k)), n - 1L)
  nn <- tryCatch(RANN::nn2(m, m, k = k_eff + 1L), error = function(e) NULL)
  if (is.null(nn)) return(NULL)
  rk <- nn$nn.dists[, k_eff + 1L]        # column 1 is the point itself (d = 0)
  rk <- rk[is.finite(rk) & rk > 0]
  if (length(rk) < 3L) return(NULL)

  # Restrict to the quantile band that the sparse background dominates.
  pr  <- sort(probs)
  qq  <- stats::quantile(rk, probs = pr, na.rm = TRUE)
  sel <- rk >= qq[1] & rk <= qq[2]
  used_band <- sum(sel) >= 3L
  rb  <- if (used_band) rk[sel] else rk

  omega <- if (d == 3L) 4 * pi / 3 else pi
  mr    <- mean(rb, na.rm = TRUE)
  if (!is.finite(mr) || mr <= 0) return(NULL)

  # Calibration constant for the band actually used (see the note above): the
  # conditional mean of G^(1/d) over the matching Gamma quantile band. Falls back
  # to the closed-form full-distribution constant if the integration fails.
  cst <- if (!used_band || (pr[1] <= 0 && pr[2] >= 1)) {
    gamma(k_eff + 1 / d) / gamma(k_eff)
  } else {
    ga <- stats::qgamma(pr[1], shape = k_eff)
    gb <- stats::qgamma(pr[2], shape = k_eff)
    val <- tryCatch(
      stats::integrate(function(g) g^(1 / d) * stats::dgamma(g, shape = k_eff),
                       lower = ga, upper = gb, rel.tol = 1e-8)$value / (pr[2] - pr[1]),
      error = function(e) NA_real_)
    if (is.finite(val) && val > 0) val else gamma(k_eff + 1 / d) / gamma(k_eff)
  }

  lam <- (cst / mr)^d / omega
  list(lambda = lam, mean_r = mr, median_r = stats::median(rk, na.rm = TRUE),
       k = k_eff, d = d, n = n, probs = pr, band_const = cst, banded = used_band)
}

## GOAL: the PER-UNIT local density, one value per unit, for the distribution plot.
##
## From the same Gamma law as nn_local_intensity but read per point rather than
## averaged: lambda_i = k / (omega_d * R_k,i^d) is the unbiased per-point estimator
## of the local intensity (E[R_k^d] = k / (lambda * omega_d)).
##
## WHY A DISTRIBUTION AND NOT JUST THE MEAN: the single global number n/support says
## how much label there is; the SHAPE says how it is arranged. A cell with genuine
## aggregates is heterogeneous - a sparse background mode plus a high-density tail -
## whereas a uniformly-labelled cell of the same mean density is unimodal. So the
## spread of this distribution is itself clustering information, which is exactly
## why its MEAN must not be used as the thinning target (see nn_local_intensity).
## Returned in units per um^d to match compute_unit_density's `density`.
unit_local_density <- function(pts, k = 3L, use_3d = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(numeric(0))
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m); d <- ncol(m)
  if (n < 3L) return(numeric(0))
  k_eff <- min(max(1L, as.integer(k)), n - 1L)
  nn <- tryCatch(RANN::nn2(m, m, k = k_eff + 1L), error = function(e) NULL)
  if (is.null(nn)) return(numeric(0))
  rk <- nn$nn.dists[, k_eff + 1L]          # column 1 is the point itself
  ok <- is.finite(rk) & rk > 0
  if (!any(ok)) return(numeric(0))
  omega <- if (d == 3L) 4 * pi / 3 else pi
  (k_eff / (omega * rk[ok]^d)) * 1000^d
}

## ═══════════════════════════════════════════════════════════════════════════════
## 8c. LOCAL DENSITY & HOTSPOTS
## ═══════════════════════════════════════════════════════════════════════════════
## WHY THIS BLOCK EXISTS
## Block 8b removes the labelling-density confound by THINNING every file down to a
## common density. That is the right move when you want to compare a curve. It is
## the wrong move when the thing you are looking for is a small number of very
## dense spots, because thinning throws away exactly the points that make them
## dense.
##
## This block takes the other route: leave the data alone and build statistics that
## are DIMENSIONLESS BY CONSTRUCTION, each unit measured against its own cell's
## mean density rather than against an absolute scale. A cell that expressed twice
## as much protein has twice the density everywhere, so a ratio to the cell mean is
## unchanged.
##
## THE FAILURE THIS BLOCK IS ANSWERING
## Mean and median nearest-neighbour distance are CENTRAL statistics. When a
## minority of molecules sit in foci and the rest are diffuse, the diffuse majority
## sets the median and the foci move it hardly at all. Measured on a real
## membrane-protein pair (a self-associating variant against a non-associating one,
## 10 acquisitions): median NN correlated with density at Spearman rho = -1.00
## within one condition and did not separate the conditions at all, while every
## statistic in this block separated them with non-overlapping ranges. NN was not
## noisy - it was measuring density.
##
## SO: read the TAIL of the local-density distribution, not its centre.
##
## WHAT IS HERE
##   compute_local_enhancement .. per-unit local density / that cell's mean density
##   compute_dispersion_index ... variance-to-mean of quadrat counts (1 under CSR)
##   calibrate_dbscan_eps ....... eps that a CSR pattern in THIS cell's own mask
##                                would barely cluster - makes eps comparable
##                                across cells of different density
##   compute_foci_table ......... one row per detected focus (n, Rg, centroid)
##   hotspot_summary_row ........ one row per file, everything above collapsed
##
## EVERY FUNCTION HERE TAKES A `support` (block 8b's footprint) AND USES IT. The
## cell mean density is n / support area, so with the bounding box in place of the
## footprint a cell imaged half in frame reads at half its true density and every
## enhancement doubles. The footprint is not optional for these statistics.

# =============================================================================
# compute_local_enhancement() - LOCAL DENSITY, IN UNITS OF THIS CELL'S MEAN
# =============================================================================
#
# GOAL
#   For every unit, estimate the density of units in its immediate neighbourhood
#   and express it as a MULTIPLE of the density the same cell would have if its
#   units were spread evenly over the footprint. 1.0 means "as dense as this cell
#   is on average"; 20 means "twenty times denser than this cell's own average".
#
# HOW THE LOCAL DENSITY IS ESTIMATED
#   The k-th nearest-neighbour distance r_k. In d dimensions, k units inside a
#   ball of radius r_k gives lambda_i = k / (omega_d * r_k^d), omega_2 = pi,
#   omega_3 = 4pi/3. This is the standard k-NN density estimator.
#
#   WHY k-NN AND NOT A FIXED-RADIUS COUNT: a fixed radius has to be chosen, and the
#   right choice depends on the density you are trying to measure - which is the
#   unknown. k-NN adapts its own scale, so one k works across a 7x density range.
#   WHY NOT A KDE: a KDE needs a bandwidth, same problem, and it smooths away the
#   isolated units that carry the "sparser background" half of the signal.
#
# CHOOSING k
#   k is a smoothing parameter, not a cluster size. Small k is noisy (k = 1 has a
#   ~52% coefficient of variation on a random pattern); large k blurs a focus into
#   its surroundings and cannot resolve a focus holding fewer than about 2k units.
#   k = 5 is the default: stable, and still resolves the ~17-unit foci that the
#   membrane-protein data set showed. Raise it for large diffuse domains, lower it
#   only if you are chasing pairs.
#
# THE BIAS THAT IS LEFT, AND WHICH WAY IT POINTS
#   The footprint estimator shrinks on a sparse cell, which raises lambda_cell and
#   therefore LOWERS enhancement. That acts AGAINST finding hotspots in sparse
#   cells, so a positive result on a sparse cell is conservative. A negative one on
#   a sparse cell is not evidence of absence.
#
# ARGUMENTS
#   pts     : data.frame with x, y (and z when use_3d) in nm, one row per unit
#   k       : neighbour order for the density estimate
#   use_3d  : measure in 3D (spheres) instead of 2D (discs)
#   support : an estimate_support() object; its $support field is the area/volume.
#             If NULL, falls back to the bounding box and flags it in used_bbox.
#
# RETURNS
#   list(enh, lambda_i, lambda_cell, r_k, k, d, n, area, used_bbox, keep)
#   `enh` follows the order of the input rows that had finite coordinates, which
#   `keep` identifies.
# =============================================================================
compute_local_enhancement <- function(pts, k = 5L, use_3d = FALSE, support = NULL) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  keep <- stats::complete.cases(m)
  m <- m[keep, , drop = FALSE]
  n <- nrow(m); d <- ncol(m)
  if (n < 3L) return(NULL)

  k_eff <- min(max(1L, as.integer(k)), n - 1L)
  nn <- tryCatch(RANN::nn2(m, m, k = k_eff + 1L), error = function(e) NULL)
  if (is.null(nn)) return(NULL)
  r_k <- nn$nn.dists[, k_eff + 1L]        # column 1 is the unit itself (d = 0)

  omega <- if (d == 3L) 4 * pi / 3 else pi
  # Guard the degenerate r_k = 0 (two units at identical coordinates, which
  # de-blinking can leave behind): it would divide by zero and produce Inf, which
  # then poisons every quantile. Floor it at half the smallest positive radius.
  pos <- r_k[is.finite(r_k) & r_k > 0]
  floor_r <- if (length(pos)) max(min(pos) / 2, 1e-6) else 1e-6
  r_k <- pmax(r_k, floor_r)
  lambda_i <- k_eff / (omega * r_k^d)

  # The cell's own mean density: units per unit area/volume of the FOOTPRINT.
  used_bbox <- FALSE
  area <- support$support %||% NA_real_
  if (!is.finite(area) || area <= 0) {
    used_bbox <- TRUE
    rng <- apply(m, 2, range)
    area <- prod(pmax(rng[2, ] - rng[1, ], 1))
  }
  lambda_cell <- n / area

  list(enh = lambda_i / lambda_cell, lambda_i = lambda_i, lambda_cell = lambda_cell,
       r_k = r_k, k = k_eff, d = d, n = n, area = area, used_bbox = used_bbox,
       keep = keep)
}

# =============================================================================
# compute_dispersion_index() - VARIANCE-TO-MEAN OF QUADRAT COUNTS
# =============================================================================
#
# GOAL
#   One dimensionless number per cell that is 1.0 under complete spatial
#   randomness BY CONSTRUCTION, so it needs no simulated null to be read.
#
# THE STATISTIC
#   Lay a grid of square (2D) or cubic (3D) quadrats of side `quadrat_nm` over the
#   window, count units in each, and take var(count)/mean(count). A Poisson process
#   has var = mean, so VMR = 1. Clustering pushes it above 1; a regular or inhibited
#   pattern below. This is the index of dispersion (Fisher, Thornton & Mackenzie
#   1922), and the reason to keep it alongside the k-NN work is precisely that it is
#   a DIFFERENT KIND of estimator - counts in fixed cells rather than distances to
#   neighbours - so agreement between the two is real evidence rather than one
#   statistic restated.
#
# THE PART THAT MUST BE RIGHT: THE EMPTY QUADRATS
#   Only quadrats inside the window may be counted, and the empty ones among them
#   MUST be included. Counting only quadrats that contain a unit throws away every
#   zero and drives VMR down towards 0; counting quadrats outside the cell invents
#   zeros and drives it up without limit. Both errors are silent. Hence the window
#   area comes from `support` when there is one, and `used_bbox` records when it
#   did not.
#
# CHOOSING THE QUADRAT SIZE
#   VMR is scale-dependent and that is a feature, not a nuisance - it is asking
#   "is there structure AT THIS SCALE". Set it near the size of the feature you are
#   looking for. On the membrane-protein data set 250 nm quadrats separated the
#   conditions with non-overlapping ranges (0.89-2.07 against 2.46-5.17) while
#   500 nm quadrats overlapped, because the foci are ~200-400 nm across and a
#   500 nm quadrat averages a focus together with its surrounding empty membrane.
#   ALWAYS report the quadrat size with the number.
#
# ARGUMENTS
#   pts        : data.frame with x, y (and z when use_3d), nm
#   quadrat_nm : quadrat side in nm
#   use_3d     : cubic quadrats instead of square
#   support    : estimate_support() object; only its area is used here
#
# RETURNS
#   list(vmr, mean_count, var_count, n_quadrats, quadrat_nm, used_bbox)
# =============================================================================
compute_dispersion_index <- function(pts, quadrat_nm = 250, use_3d = FALSE,
                                     support = NULL) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < 3L) return(NULL)
  q <- max(1, as.numeric(quadrat_nm))
  d <- ncol(m)

  # Index every unit into a quadrat on a lattice anchored at the data minimum.
  o   <- apply(m, 2, min)
  idx <- floor(sweep(m, 2, o, "-") / q)
  key <- apply(idx, 1, paste, collapse = ",")
  occupied <- table(key)

  # How many quadrats does the window hold? Empty ones must be counted, so this is
  # an AREA divided by a quadrat area, not a count of occupied cells.
  used_bbox <- FALSE
  area <- support$support %||% NA_real_
  if (!is.finite(area) || area <= 0) {
    used_bbox <- TRUE
    rng  <- apply(m, 2, range)
    area <- prod(pmax(rng[2, ] - rng[1, ], q))
  }
  n_q <- max(length(occupied), round(area / q^d))

  cnt <- c(as.integer(occupied), rep(0L, max(0L, n_q - length(occupied))))
  mu  <- mean(cnt)
  if (!is.finite(mu) || mu <= 0) return(NULL)
  list(vmr = stats::var(cnt) / mu, mean_count = mu, var_count = stats::var(cnt),
       n_quadrats = length(cnt), quadrat_nm = q, used_bbox = used_bbox)
}

# =============================================================================
# calibrate_dbscan_eps() - ONE eps THAT MEANS THE SAME THING IN EVERY CELL
# =============================================================================
#
# GOAL
#   Choose DBSCAN's eps per cell so that "clustered" means the same thing in a
#   sparse cell and a dense one.
#
# WHY A FIXED eps CANNOT WORK ACROSS CELLS
#   DBSCAN with (eps, minPts) is a DENSITY THRESHOLD: a unit is core when at least
#   minPts units lie within eps, i.e. when the local density exceeds
#   minPts / (omega_d * eps^d). That threshold is absolute, but the cells are not -
#   across a 7x range of labelling density the same eps is crossed trivially in the
#   dense cells and never in the sparse ones. The number of "clusters" then tracks
#   transfection load, which is the confound the whole Analysis tab exists to avoid.
#
# THE CALIBRATION
#   For a candidate eps, scatter n units UNIFORMLY AT RANDOM INSIDE THIS CELL'S OWN
#   FOOTPRINT, run the same DBSCAN, and record what fraction of that random pattern
#   is assigned to a cluster. Repeat `nsim` times and average. Take the LARGEST eps
#   whose random-clustering fraction is still at or below `target_frac`. That eps is
#   the most permissive setting at which chance alone would not produce clusters, so
#   a cluster found at it is a real density excess for THIS cell and the resulting
#   focus statistics are comparable between cells.
#
#   Randomising inside the FOOTPRINT rather than the bounding box is what makes this
#   correct for a cell that does not fill the frame; scattering into empty field
#   would spread the random units too thin and hand back an eps far too permissive.
#
# WHAT IT DOES NOT FIX
#   Cluster COUNT is still not comparable and neither is the fraction of units in a
#   cluster - on real data that fraction reached 31% in a control cell, because at
#   low density the calibrated eps grows large enough to sweep up diffuse
#   background. The statistic that survived was the SIZE OF THE LARGEST focus.
#
# ARGUMENTS
#   pts         : data.frame with x, y (and z when use_3d)
#   support     : estimate_support() object WITH $mask - the whole point is to
#                 randomise inside the cell. Without a mask it falls back to the
#                 bounding box, which weakens the calibration.
#   minPts      : DBSCAN minPts, held fixed across cells
#   eps_grid    : candidate eps values in nm
#   target_frac : the false-cluster rate to calibrate to (0.05 = 5%)
#   nsim        : CSR patterns per candidate. 15 is enough to place the crossing
#                 between two grid rungs; more only sharpens a number the grid
#                 spacing has already quantised.
#
# RETURNS
#   list(eps, curve = data.frame(eps, frac), target_frac, minPts, hit_floor)
#   `hit_floor` is TRUE when even the smallest candidate clustered too much, i.e.
#   the grid did not reach low enough - the returned eps is then the grid minimum
#   and the calibration did NOT converge.
# =============================================================================
calibrate_dbscan_eps <- function(pts, support = NULL, minPts = 5L,
                                 eps_grid = c(40, 60, 80, 100, 125, 150, 175,
                                              200, 250, 300, 400, 500),
                                 target_frac = 0.05, nsim = 15L, use_3d = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m)
  if (n < max(3L, minPts)) return(NULL)
  eps_grid <- sort(unique(eps_grid[is.finite(eps_grid) & eps_grid > 0]))
  if (!length(eps_grid)) return(NULL)
  nsim <- max(3L, as.integer(nsim))

  # Randomiser: inside the footprint when there is a mask, else inside the bbox.
  rng <- apply(m, 2, range)
  draw <- function() {
    p <- if (!is.null(support) && !is.null(support$mask))
      tryCatch(sample_in_support(support, n, use_3d = use_3d), error = function(e) NULL)
    else NULL
    if (is.null(p)) {
      out <- vapply(seq_len(ncol(m)),
                    function(j) stats::runif(n, rng[1, j], rng[2, j]), numeric(n))
      return(matrix(out, nrow = n))
    }
    if (use_3d) cbind(p$x, p$y, p$z %||% rep(0, n)) else cbind(p$x, p$y)
  }

  frac <- vapply(eps_grid, function(e) {
    mean(vapply(seq_len(nsim), function(i) {
      cl <- tryCatch(dbscan::dbscan(draw(), eps = e, minPts = minPts)$cluster,
                     error = function(err) NULL)
      if (is.null(cl)) NA_real_ else mean(cl > 0L)
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))

  ok <- is.finite(frac) & frac <= target_frac
  hit_floor <- !any(ok)
  eps <- if (hit_floor) eps_grid[1] else max(eps_grid[ok])
  list(eps = eps, curve = data.frame(eps = eps_grid, frac = frac),
       target_frac = target_frac, minPts = minPts, hit_floor = hit_floor)
}

# =============================================================================
# compute_foci_table() - ONE ROW PER DETECTED HOTSPOT
# =============================================================================
#
# GOAL
#   Run DBSCAN at a given eps and describe each cluster it finds, so the headline
#   number can be a count of UNITS IN THE BIGGEST FOCUS rather than a curve.
#
# WHY Rg IS BIAS-CORRECTED
#   The plain radius of gyration is biased low at small n by the factor (n-1)/n,
#   because the centroid is estimated from the same points. Uncorrected, that
#   manufactures a "bigger clusters are bigger" trend out of nothing (29% low at
#   n = 2, 11% at n = 5, 1% at n = 50). The Bessel correction is applied here, the
#   same as compute_cluster_rog(). A residual Jensen bias survives the square root
#   at very small n and has no closed form.
#
# WHAT Rg DOES NOT TELL YOU
#   It is a single scalar for a shape that need not be round. Real foci in the
#   membrane-protein data came out as two lobes and as elongated chains, and Rg
#   reports both as "about 150 nm". `aniso` (the ratio of the largest to smallest
#   principal standard deviation) is returned alongside for exactly that reason:
#   read it before treating a focus as a disc.
#
# ARGUMENTS
#   pts    : data.frame with x, y (and z when use_3d)
#   eps    : DBSCAN eps in nm, normally from calibrate_dbscan_eps()
#   minPts : DBSCAN minPts
#
# RETURNS
#   list(table, cluster, keep, eps). `table` is one row per focus (focus, n, rg,
#   aniso, x, y, z) ordered largest first; `cluster` is the per-unit assignment
#   (0 = not in any focus) for the rows `keep` marks.
# =============================================================================
compute_foci_table <- function(pts, eps, minPts = 5L, use_3d = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  keep <- stats::complete.cases(m)
  m <- m[keep, , drop = FALSE]
  if (nrow(m) < max(3L, minPts) || !is.finite(eps) || eps <= 0) return(NULL)

  cl <- tryCatch(dbscan::dbscan(m, eps = eps, minPts = minPts)$cluster,
                 error = function(e) NULL)
  if (is.null(cl)) return(NULL)

  empty <- data.frame(focus = integer(0), n = integer(0), rg = numeric(0),
                      aniso = numeric(0), x = numeric(0), y = numeric(0),
                      z = numeric(0))
  ids <- sort(unique(cl[cl > 0L]))
  if (!length(ids))
    return(list(table = empty, cluster = cl, keep = keep, eps = eps))

  rows <- lapply(ids, function(id) {
    P   <- m[cl == id, , drop = FALSE]
    ctr <- colMeans(P)
    Q   <- sweep(P, 2, ctr, "-")
    nn  <- nrow(Q)
    # Bessel-corrected mean squared distance from the centroid.
    rg <- sqrt(sum(Q^2) / nn * (if (nn > 1L) nn / (nn - 1L) else 1))
    # Principal spread ratio: above ~1.6 the focus is not a disc, so look at it.
    sdv <- tryCatch(stats::prcomp(Q, center = FALSE)$sdev, error = function(e) NA_real_)
    an  <- if (all(is.finite(sdv)) && min(sdv) > 0) max(sdv) / min(sdv) else NA_real_
    data.frame(focus = id, n = nn, rg = rg, aniso = an,
               x = ctr[1], y = ctr[2],
               z = if (length(ctr) >= 3L) ctr[3] else NA_real_)
  })
  tb <- do.call(rbind, rows)
  tb <- tb[order(-tb$n), , drop = FALSE]
  list(table = tb, cluster = cl, keep = keep, eps = eps)
}

# =============================================================================
# csr_cv_local_density() - WHAT THE CV READS WHEN NOTHING IS CLUSTERED
# =============================================================================
#
# GOAL
#   Give the CV of local density a zero point, so it can be read as a multiple of
#   random rather than as a bare number.
#
# THE DERIVATION
#   Under a Poisson process the k-th neighbour distance satisfies
#   lambda*omega_d*r_k^d ~ Gamma(k, 1), and the estimator is
#   lambda_i = k / (omega_d r_k^d), so lambda_i is proportional to 1/Gamma(k, 1).
#   For X ~ Gamma(k, 1): E[1/X] = 1/(k-1) and Var(1/X) = 1/((k-1)^2 (k-2)), so
#
#       CV = sd/mean = 1 / sqrt(k - 2)
#
#   It depends ONLY on k. Not on the density, not on the area, not on the number
#   of units, and not on the dimension. That is what makes the raw CV comparable
#   between cells in the first place.
#
#   k = 3 -> 1.000, k = 5 -> 0.577, k = 10 -> 0.354, k = 20 -> 0.236.
#
# UNDEFINED BELOW k = 3
#   The variance of 1/Gamma(k,1) exists only for k > 2, and the mean only for
#   k > 1. At k = 1 or 2 the sampling distribution has no finite spread and the
#   ratio is meaningless, so this returns NA and the ratio column goes NA with it.
#
# VALIDATED AGAINST SIMULATION, NOT ASSUMED
#   The formula ignores two things the real estimator has: neighbourhoods overlap,
#   so the per-unit estimates are spatially correlated, and units near the edge of
#   the footprint lose neighbours. Both could shift the baseline. Measured on the
#   10 membrane-protein footprints at k = 5, 60 randomisations each, drawing into each file's
#   own footprint at its own n: the simulated median CV came out 0.53-0.59 against
#   the predicted 0.577, and it did NOT track the unit count (Spearman -0.07 over
#   an 86-546 range). So the closed form is usable as-is and no simulation is
#   needed at run time.
# =============================================================================
csr_cv_local_density <- function(k) {
  k <- suppressWarnings(as.numeric(k))
  if (!length(k) || !is.finite(k) || k <= 2) return(NA_real_)
  1 / sqrt(k - 2)
}

# =============================================================================
# hotspot_summary_row() - EVERYTHING ABOVE, COLLAPSED TO ONE ROW PER FILE
# =============================================================================
#
# GOAL
#   The row that the comparison table and the metric-versus-density plot both read.
#
# WHICH COLUMNS TO TRUST, MEASURED
#   On the membrane-protein pair (4 control against 6 self-associating cells, whose
#   densities were NOT significantly different, rank-sum p = 0.92) the columns with
#   completely non-overlapping ranges, ranked by how significant condition stayed
#   after regressing out log density, were:
#     max_focus_n (p 0.0005) > cv_local (0.0007) > enh_q90 (0.0025) >
#     vmr (0.0030) > frac_above (0.0057).
#   `cv_ratio` is `cv_local` divided by a constant that depends only on k, so it
#   ranks identically; it exists to be readable, not to add information. On that
#   data cond01 ran 1.6-2.4x random and cond02 3.6-6.8x.
#   `nn_med` and `clark_evans` are carried in the same row deliberately, so the
#   table shows the two that FAILED next to the ones that worked rather than
#   quietly omitting them.
#
# ARGUMENTS
#   pts        : one row per unit (from unit_centroids())
#   support    : estimate_support() object with $mask
#   k          : k for the local density estimate
#   quadrat_nm : quadrat side for the dispersion index
#   enh_thresh : the enhancement multiple that `frac_above` counts
#   eps_cal    : a calibrate_dbscan_eps() result, or NULL to skip the focus columns
#   minPts     : DBSCAN minPts
#
# RETURNS
#   a one-row data.frame, or NULL when there are too few units.
# =============================================================================
hotspot_summary_row <- function(pts, support = NULL, k = 5L, quadrat_nm = 250,
                                enh_thresh = 5, eps_cal = NULL, minPts = 5L,
                                use_3d = FALSE, dataset = NA_character_) {
  le <- compute_local_enhancement(pts, k = k, use_3d = use_3d, support = support)
  if (is.null(le)) return(NULL)
  di <- compute_dispersion_index(pts, quadrat_nm = quadrat_nm, use_3d = use_3d,
                                 support = support)

  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m); d <- ncol(m)

  # 1-NN distance and the Clark-Evans ratio, carried as the negative controls.
  nn1 <- tryCatch(RANN::nn2(m, m, k = 2L)$nn.dists[, 2], error = function(e) NA_real_)
  lam <- le$lambda_cell
  # Expected mean 1-NN under CSR: 0.5 * lambda^(-1/2) in 2D,
  # gamma(4/3) * (3 / (4 pi lambda))^(1/3) in 3D.
  exp_nn <- if (d == 3L) gamma(4/3) * (3 / (4 * pi * lam))^(1/3) else 0.5 / sqrt(lam)
  ce <- if (is.finite(exp_nn) && exp_nn > 0) mean(nn1, na.rm = TRUE) / exp_nn else NA_real_

  ft <- if (!is.null(eps_cal))
    compute_foci_table(pts, eps = eps_cal$eps, minPts = minPts, use_3d = use_3d) else NULL
  tb <- ft$table
  # Volume in um^3 for 3D, area in um^2 for 2D - the scale factor differs.
  scale_unit <- if (d == 3L) 1e9 else 1e6

  data.frame(
    dataset      = as.character(dataset),
    n_units      = n,
    area_um2     = le$area / scale_unit,
    density      = lam * scale_unit,
    # --- the readouts that worked ---
    cv_local     = stats::sd(le$lambda_i) / mean(le$lambda_i),
    # The same number as a multiple of what a random pattern gives at this k, so
    # 1.0 reads as "random" and 3.0 as "three times as uneven as random".
    cv_ratio     = (stats::sd(le$lambda_i) / mean(le$lambda_i)) /
                     csr_cv_local_density(le$k),
    enh_q90      = unname(stats::quantile(le$enh, 0.90, na.rm = TRUE)),
    enh_q99      = unname(stats::quantile(le$enh, 0.99, na.rm = TRUE)),
    frac_above   = mean(le$enh > enh_thresh, na.rm = TRUE),
    vmr          = di$vmr %||% NA_real_,
    eps_nm       = eps_cal$eps %||% NA_real_,
    n_foci       = if (is.null(tb)) NA_integer_ else nrow(tb),
    max_focus_n  = if (is.null(tb) || !nrow(tb)) 0L else as.integer(tb$n[1]),
    max_focus_rg = if (is.null(tb) || !nrow(tb)) NA_real_ else tb$rg[1],
    max_focus_aniso = if (is.null(tb) || !nrow(tb)) NA_real_ else tb$aniso[1],
    frac_in_foci = if (is.null(ft)) NA_real_ else mean(ft$cluster > 0L),
    # --- the two that FAILED, kept visible on purpose ---
    nn_med       = stats::median(nn1, na.rm = TRUE),
    clark_evans  = ce,
    used_bbox    = isTRUE(le$used_bbox),
    stringsAsFactors = FALSE)
}

## ── Series colours for the batch/spatial panels ──────────────────────────────
## GOAL: two palettes with different jobs.
##
## `qual_palette(n)` - CATEGORICAL colours for one-series-per-file plots. Files are
## nominal, not ordered, so a sequential ramp is the wrong encoding: viridis on 10
## files reads as a blue-green gradient in which neighbouring files are nearly
## indistinguishable. This is a hand-ordered set of maximally separated hues
## (Dark2 then Set1 then darker Paired entries), recycled with a darkening pass if
## more series than colours. DETERMINISTIC, not random - a random draw would
## re-colour every series on every redraw and can collide by chance.
##
## `shade_colors(base, n)` - ORDERED shades of one base colour, for individual
## files inside a group. Here the encoding is right: the group identity is the hue
## and the file is a shade of it, so a red group and a blue group stay readable at a
## glance while individual files remain separable.
QUAL_COLORS <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
  "#a6761d", "#666666", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
  "#ff7f00", "#a65628", "#f781bf", "#1f78b4", "#33a02c", "#6a3d9a",
  "#b15928", "#b2182b"
)

qual_palette <- function(n) {
  n <- max(1L, as.integer(n))
  base <- QUAL_COLORS
  if (n <= length(base)) return(base[seq_len(n)])
  # More series than distinct hues: recycle, darkening each successive lap so a
  # repeat is still visually distinct from its first appearance.
  reps <- ceiling(n / length(base))
  out  <- unlist(lapply(seq_len(reps) - 1L, function(k)
    if (k == 0L) base else .adjust_lum(base, -0.18 * k)))
  out[seq_len(n)]
}

# Shift luminance of hex colours: amt > 0 lightens toward white, < 0 darkens.
.adjust_lum <- function(hex, amt) {
  m <- grDevices::col2rgb(hex) / 255
  m <- if (amt >= 0) m + (1 - m) * amt else m * (1 + amt)
  m <- pmin(pmax(m, 0), 1)
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

shade_colors <- function(base, n) {
  n <- max(1L, as.integer(n))
  if (n == 1L) return(base)
  # Light tint -> base -> dark shade. Kept away from the extremes so no file ends
  # up near-white (invisible) or near-black (indistinguishable from another group).
  grDevices::colorRampPalette(c(.adjust_lum(base, 0.55), base,
                                .adjust_lum(base, -0.30)))(n)
}

## GOAL: the x-axis range for the density-distribution plot, and the clipping span
## that goes with it.
##
## WHY NOT AUTORANGE: a violin's kernel tails run well past the data, and on a log
## axis plotly then pads out to whole decade boundaries - which leaves the panel
## long and mostly empty. So the span is taken from a robust quantile of the local
## densities (`probs`), widened to include the per-file global densities so those
## markers can never fall off frame, and set explicitly.
##
## Values outside `lim` are meant to be CLIPPED to the edge by the caller rather
## than dropped, so each violin ends at the frame instead of trailing off it.
##
## Returns list(lim, range, clip) where `range` is what plotly's xaxis wants
## (log10 of the padded limits when use_log, the padded limits otherwise) and
## `clip` is a function mapping values into `lim`.
density_axis_range <- function(local_vals, global_vals, use_log = TRUE,
                               probs = c(0.005, 0.995)) {
  lv <- local_vals[is.finite(local_vals) & local_vals > 0]
  gv <- global_vals[is.finite(global_vals)]
  # Both empty is reachable (a panel rendering before any measurement); range() on
  # nothing warns and returns c(Inf, -Inf), so short-circuit to a harmless frame.
  if (!length(lv) && !length(gv))
    return(list(lim = c(1, 10), range = if (use_log) c(0, 1) else c(0, 10),
                clip = function(v) v))
  lim <- if (length(lv) >= 10L)
    as.numeric(stats::quantile(lv, probs, na.rm = TRUE)) else range(c(lv, gv))
  lim <- range(c(lim, gv))
  if (!all(is.finite(lim)) || lim[2] <= lim[1])
    lim <- c(max(1e-9, min(lim, na.rm = TRUE)), max(1, max(lim, na.rm = TRUE)) * 10)
  if (use_log) {
    pos <- c(lv, gv[gv > 0])
    if (length(pos)) lim[1] <- max(lim[1], min(pos) / 1.6)
    lim[1] <- max(lim[1], 1e-9)
    rng <- log10(pmax(c(lim[1] / 1.25, lim[2] * 1.25), 1e-12))
  } else {
    sp  <- diff(lim)
    rng <- c(max(0, lim[1] - 0.04 * sp), lim[2] + 0.04 * sp)
  }
  list(lim = lim, range = rng,
       clip = function(v) pmin(pmax(v, lim[1]), lim[2]))
}

## GOAL: the per-file row of the density panel - everything needed to decide
## whether two acquisitions are comparable, and to thin them if they are not.
##
## Returns a one-row data.frame:
##   n_units      units in the file (trace IDs or de-blinked molecules)
##   support_nm   estimated support (nm^2 or nm^3), by `support_mode`
##   density      units per um^2 (2D) or per um^3 (3D)   <- the THINNING TARGET
##   density_nn   the same quantity from nn_local_intensity (support-free, but
##                clustering-confounded - a cross-check, not a target)
##   nn_median    median k-th NN distance (nm), the raw observable
##   frac_of_bbox support / bounding box - a framing/edge diagnostic. Well below
##                1 means a lot of empty frame or dark region was excluded, which
##                is the case where the bounding-box density would have been wrong.
##
## ⚠ `support_override`: THE SUPPORT ESTIMATE IS ITSELF DENSITY-DEPENDENT at low n.
## An occupancy grid built from fewer points finds fewer occupied cells, so
## re-estimating the support after thinning shrinks A and inflates n/A - which
## would make the thinning look like it had failed to equalise anything (measured:
## a 4x labelling difference read as 3.4x before thinning and still 1.16x after,
## purely from this). A cell's footprint does not change when you thin it, so pass
## the PRE-THINNING support here whenever the pattern has been thinned. The app's
## thinning gate does exactly that.
##
## ⚠ `support_from`: WHICH POINTS DEFINE THE FOOTPRINT. The support is a property
## of WHERE SIGNAL WAS DETECTED, so it is best estimated from EVERY LOCALIZATION
## ("locs", the default) rather than from the unit centroids. Using centroids makes
## the support depend on the unit count, which is the very thing being normalised:
## measured on two synthetic cells with identical geometry and a 4x labelling
## difference, centroid-based supports were 5.01e6 vs 3.82e6 nm^2 for the same
## region, so the density ratio read 3.05 instead of 4.0.
compute_unit_density <- function(df, unit = "tid", use_3d = FALSE,
                                support_mode = "grid", h_nm = 100, close_nm = 200,
                                k = 3L, support_override = NULL,
                                support_from = c("locs", "units"),
                                bw_nm = 150, frac = 0.99) {
  support_from <- match.arg(support_from)
  pts <- unit_centroids(df, unit = unit, use_3d = use_3d)
  if (is.null(pts) || nrow(pts) < 3L) return(NULL)
  sup_pts <- if (support_from == "locs") {
    keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
    if (!all(keep %in% names(df))) pts else {
      d2 <- df[, keep, drop = FALSE]
      if (!use_3d) d2$z <- 0
      d2[stats::complete.cases(d2), , drop = FALSE]
    }
  } else pts
  sup <- if (!is.null(support_override) && is.finite(support_override) && support_override > 0) {
    list(support = as.numeric(support_override), mode = "inherited (pre-thinning)",
         h = h_nm, close_nm = close_nm, frac_of_bbox = NA_real_)
  } else {
    estimate_support(sup_pts, use_3d = use_3d, mode = support_mode,
                     h_nm = h_nm, close_nm = close_nm,
                     bw_nm = bw_nm, frac = frac)
  }
  if (is.null(sup) || !is.finite(sup$support) || sup$support <= 0) return(NULL)
  nni <- nn_local_intensity(pts, k = k, use_3d = use_3d)

  # Dimensionality follows the ESTIMATOR, not the 3-D checkbox: the KDE/HDR support
  # is a plan-view area even on 3-D data, so treating it as a volume would divide by
  # nm^3 and report a density three orders of magnitude wrong.
  d      <- if (!is.null(sup$d)) as.integer(sup$d) else if (use_3d) 3L else 2L
  per_um <- 1000^d                       # nm^d -> um^d
  data.frame(
    n_units      = nrow(pts),
    support_nm   = sup$support,
    support_mode = sup$mode,
    h_nm         = sup$h %||% NA_real_,
    close_nm     = sup$close_nm %||% NA_real_,
    frac_of_bbox = sup$frac_of_bbox,
    density      = nrow(pts) / sup$support * per_um,
    density_nn   = if (!is.null(nni)) nni$lambda * per_um else NA_real_,
    nn_median    = if (!is.null(nni)) nni$median_r else NA_real_,
    nn_k         = if (!is.null(nni)) nni$k else NA_integer_,
    dims         = d,
    stringsAsFactors = FALSE
  )
}

## GOAL: independent p-thinning of a localization table, AT THE UNIT LEVEL.
##
## Each unit (trace ID, or de-blinked molecule) is kept or dropped as a WHOLE,
## with all of its localizations. Thinning individual localizations instead would
## corrupt the thing being normalised: it would shorten traces, change their
## spread and precision, and break de-blinking - and it would not thin the point
## pattern at all, since a partially-sampled trace still contributes one centroid.
##
## `p` >= 1 is a no-op (never up-samples - you cannot invent molecules).
## The draw is deterministic given `seed`, and the GLOBAL RNG stream is SAVED AND
## RESTORED around it, so enabling thinning cannot perturb any other randomised
## result in the app (the permutation nulls, the CSR envelopes, the alpha-shape
## subsampling all draw from the same stream).
thin_units_to_p <- function(df, unit = "tid", p = 1, seed = 1L) {
  if (is.null(df) || !nrow(df)) return(df)
  if (!is.finite(p) || p >= 1) return(df)
  if (p <= 0) return(df[0, , drop = FALSE])
  key <- if (unit == "cluster" && "cluster" %in% names(df)) "cluster" else "tid"
  if (!key %in% names(df)) return(df)

  ids <- unique(df[[key]])
  # De-blinking noise (id 0 / NA) is not a unit; it is dropped from the thinned
  # pattern rather than kept wholesale, which would leave a density floor that no
  # thinning could reach.
  if (key == "cluster") ids <- ids[!is.na(ids) & ids > 0L]
  if (!length(ids)) return(df[0, , drop = FALSE])

  had  <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_backup <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", seed_backup, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  set.seed(seed)
  # Fixed-count draw rather than n independent coin flips: same expectation, no
  # binomial noise in the achieved density, which matters because the whole point
  # is to land every file on ONE common density.
  n_keep <- max(1L, as.integer(round(length(ids) * p)))
  keep   <- sample(ids, n_keep)
  df[df[[key]] %in% keep, , drop = FALSE]
}

## ── Monte-Carlo CSR envelope for L(r) / L₁₂(r) / g₁₂(r) / g(r) ───────────────
## Generates `nsim` complete-spatial-randomness (CSR) patterns of `n` points
## (and `n2` for the cross statistics) uniformly inside the SAME window as the
## observed pattern (`window_ranges`, a 2×d matrix of per-axis min/max), and runs
## the identical estimator on each. Returns pointwise lo/hi quantiles (default
## 2.5% / 97.5% → a 95% band; a two-sided α≈0.05 test needs nsim ≥ 39) and the
## mean null curve, on the supplied `r_seq`.
##
## Because the sims share the observed window, the (uncorrected) edge bias that
## drags L(r)−r negative at large r is present in the envelope too and cancels
## when the observed curve is read against the band - this envelope is the
## interpretive tool, superseding an explicit edge correction.
## `support`: when supplied, the CSR patterns are drawn uniformly ON THE SUPPORT
## instead of in the bounding box, and the estimators use the support measure. This
## is the mask-aware randomization the 2026-07-30 audit flagged as outstanding: a
## box-CSR null puts simulated points where the sample never was, so it under-states
## how clustered a random pattern INSIDE the cell already looks, and the observed
## curve leaves the band too easily.
## `progress`: optional function(frac) called after each simulation. The envelope is
## the slow part by a wide margin, so a bar that only advances per FILE looks frozen -
## pass a callback and it advances per simulation instead.
compute_ripley_envelope <- function(n, window_ranges, r_seq, use_3d = FALSE,
                                     nsim = 99L, estimator = c("L", "L12", "pcf", "g"),
                                     n2 = NULL, probs = c(0.025, 0.975),
                                     support = NULL, edge_correct = TRUE,
                                     progress = NULL) {
  estimator <- match.arg(estimator)
  if (is.null(r_seq) || !length(r_seq)) return(NULL)
  if (is.null(window_ranges) || any(!is.finite(window_ranges))) return(NULL)
  nsim <- max(1L, as.integer(nsim))
  n    <- as.integer(n)
  if (n < 3L) return(NULL)
  # The one-sample estimators ("L", "g") need no second point count; the cross
  # statistics do.
  if (!estimator %in% c("L", "g")) {
    if (is.null(n2)) return(NULL)
    n2 <- as.integer(n2)
    if (n2 < 2L) return(NULL)
  }

  d    <- ncol(window_ranges)
  lo_v <- window_ranges[1, ]; hi_v <- window_ranges[2, ]
  is3d <- isTRUE(use_3d) && d >= 3L

  # Draw on the support when one is given, else uniformly in the box.
  use_sup <- !is.null(support) && !is.null(support$mask)
  gen <- function(m) {
    if (use_sup) {
      s <- sample_in_support(support, m, use_3d = is3d)
      if (!is.null(s)) return(s)
    }
    list(x = stats::runif(m, lo_v[1], hi_v[1]),
         y = stats::runif(m, lo_v[2], hi_v[2]),
         z = if (is3d) stats::runif(m, lo_v[3], hi_v[3]) else NULL)
  }

  nr  <- length(r_seq)
  mat <- matrix(NA_real_, nrow = nr, ncol = nsim)
  cap <- if (estimator %in% c("L", "g")) n else max(n, n2)

  # Set covariance ONCE for the whole loop - it is a property of the window, and
  # rebuilding it per simulation is what made this look like a hang.
  if (!is.null(support) && isTRUE(edge_correct)) {
    d_max <- sqrt(sum((hi_v - lo_v)^2))
    support <- attach_gamma_cache(support, r_max = max(max(r_seq), d_max),
                                 use_3d = is3d)
  }

  for (s in seq_len(nsim)) {
    p1 <- gen(n)
    col <- if (estimator == "L") {
      cv <- compute_ripley_L(p1$x, p1$y, p1$z, r_seq = r_seq, use_3d = is3d,
                             max_subsample = cap, window_ranges = window_ranges,
                             support = support, edge_correct = edge_correct)
      cv$Lr_minus_r
    } else if (estimator == "g") {
      cv <- compute_pcf(p1$x, p1$y, p1$z, r_seq = r_seq, use_3d = is3d,
                        max_subsample = cap, window_ranges = window_ranges,
                        support = support, edge_correct = edge_correct)
      cv$g
    } else if (estimator == "L12") {
      p2 <- gen(n2)
      cv <- compute_ripley_L12(p1$x, p1$y, p1$z, p2$x, p2$y, p2$z,
                               r_seq = r_seq, use_3d = is3d,
                               max_subsample = cap, window_ranges = window_ranges,
                               support = support)
      cv$Lr_minus_r
    } else {
      p2 <- gen(n2)
      cv <- compute_pcf12(p1$x, p1$y, p1$z, p2$x, p2$y, p2$z,
                          r_seq = r_seq, use_3d = is3d,
                          max_subsample = cap, window_ranges = window_ranges,
                          support = support)
      cv$g12
    }
    if (length(col) == nr) {
      mat[, s] <- col
    } else if (length(col) > 0L) {
      k <- min(length(col), nr)
      mat[seq_len(k), s] <- col[seq_len(k)]
    }
    if (is.function(progress)) try(progress(1 / nsim), silent = TRUE)
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
## ⚠ CALIBRATION - measured, not assumed (2026-07-30). Simulated 30× per scenario,
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
## `support`: when supplied, BOTH nulls randomize over the sample support instead of
## the bounding box. This is the fix the 2026-07-30 audit named as outstanding - it
## measured 53% (shift) and 93% (CSR) false-positive rates when both channels shared
## one cell inside a larger window, because a box-randomized Ch2 gets scattered into
## regions where no sample exists, decorrelating it far more than a real
## rearrangement inside the cell would. See the measured post-fix rates in
## scratchpad/test_support_window.R and the README.
compute_kde_coloc <- function(x1, y1, x2, y2, bandwidth = NULL, n_grid = 100L,
                              bw_sweep = NULL, costes_nsim = 0L,
                              null_mode = c("shift", "csr"), support = NULL) {
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
    use_sup <- !is.null(support) && !is.null(support$mask)
    null_r <- vapply(seq_len(costes_nsim), function(i) {
      if (identical(null_mode, "shift")) {
        # Toroidal shift: wrap channel 2 inside the SAME window. Preserves its own
        # clustering and n exactly; only the cross-channel registration is broken.
        # With a support, shifts are REJECTION-SAMPLED so the shifted channel still
        # lands inside the cell: a shift that throws most of Ch2 into empty space is
        # not a plausible rearrangement of the sample, and accepting it is exactly
        # what made this null anti-conservative.
        accept <- function(rx, ry) mean(in_support(support, rx, ry)) >= 0.9
        rx <- ry <- NULL
        for (tries in seq_len(if (use_sup) 60L else 1L)) {
          sx <- stats::runif(1, 0, wx); sy <- stats::runif(1, 0, wy)
          cx <- xlim[1] + ((x2 - xlim[1] + sx) %% wx)
          cy <- ylim[1] + ((y2 - ylim[1] + sy) %% wy)
          if (!use_sup || accept(cx, cy)) { rx <- cx; ry <- cy; break }
        }
        if (is.null(rx)) {
          # No acceptable shift found (a cell filling most of its box leaves little
          # room). Fall back to a support-restricted CSR for this replicate rather
          # than silently using an out-of-cell shift.
          s <- sample_in_support(support, n2)
          rx <- s$x; ry <- s$y
        }
      } else if (use_sup) {
        s <- sample_in_support(support, n2)
        rx <- s$x; ry <- s$y
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
    out$null_support     <- use_sup
    out$null_label       <- paste0(
      if (identical(null_mode, "shift"))
        "random toroidal shift of Ch2 (structure-preserving; Lotwick & Silverman 1982)"
      else "uniform CSR resampling of Ch2 (tests against complete spatial randomness)",
      if (use_sup) ", randomized over the sample support" else ", randomized over the bounding box")
  }

  out
}

## ── Grid-based density correlation (multi-scale) ─────────────────────────────
## For each cell size in `cell_sizes`, bins Ch1 and Ch2 molecules into a 2D
## grid and returns Pearson r between per-cell counts - a multi-scale
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
