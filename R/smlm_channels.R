## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_channels.R  -  Channels: two-colour split and multiplex
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Where a channel is defined (resolve_channels), the ratiometric DCR split
##   and its per-trace confidence, the multiplex channel pair, and channel
##   colours. Grouping never crosses a channel; distances may, pairwise only.
##
## CONTENTS (in file order)
##   assign_two_color_group, auto_detect_dcr_split, dcr_split_confidence,
##   assign_analysis_group, resolve_channels, attach_channel_colors,
##   channels_in_scope, pair_passes_scope, mplex_default_colors,
##   assign_analysis_group_from_mplex, format_tc_id, wrap_groups,
##   wrap_group_rows, .mplex_key, .mplex_adopt_csv, channel_parts,
##   split_df_by_channel, .facet_name, get_tc_pts
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

assign_two_color_group <- function(df, settings) {
  # Thin wrapper. The DCR/fluo split itself lives in resolve_channels(), which
  # is the only place a channel is defined; this used to be a second copy of it
  # (and server_smlm had a third). Callers here want a per-TRACE map, so that is
  # what it returns - one row per tid with analysis_group and the medians.
  if (!is.data.frame(df) || !nrow(df)) return(df[0, , drop = FALSE])
  if (!"z" %in% names(df)) df$z <- 0
  d <- resolve_channels(df, smlm_set(settings, enable_two_color = TRUE))
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

# GOAL: tag each row with analysis_group 1/2 from its trace's median DCR or fluo.
#   settings      : smlm_settings or Shiny input (reads use_pooled_dcr,
#                   color_by_fluo, dcr_group1, dcr_group2)
#   grouping_mode : "none", "dcr" or "fluo". The app has no control for it, so
#                   it is "none" there; color_by_fluo switches to "fluo".
assign_analysis_group <- function(df, settings, grouping_mode = "none") {
  s <- as_smlm_settings(settings)
  .dcr_col <- if (isTRUE(s$use_pooled_dcr) && "dcr_pooled" %in% names(df))
    "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL
  med_tbl <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      dcr_med  = if (!is.null(.dcr_col)) stats::median(.data[[.dcr_col]], na.rm = TRUE) else NA_real_,
      fluo_med = if ("fluo" %in% names(df)) stats::median(fluo, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  df <- dplyr::left_join(df, med_tbl, by = "tid")
  
  mode <- grouping_mode %||% "none"
  if (isTRUE(s$color_by_fluo)) mode <- "fluo"
  if (mode == "fluo") {
    df$analysis_group <- dplyr::case_when(
      is.finite(df$fluo_med) & round(df$fluo_med) == 1 ~ 1L,
      is.finite(df$fluo_med) & round(df$fluo_med) == 2 ~ 2L,
      TRUE ~ NA_integer_
    )
  } else if (mode == "dcr") {
    g1 <- s$dcr_group1; g2 <- s$dcr_group2
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
resolve_channels <- function(df, settings) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  if (!"z" %in% names(df)) df$z <- 0
  s <- as_smlm_settings(settings)

  two_color_on <- isTRUE(s$enable_two_color)
  mplex_on     <- "mplex_channel" %in% names(df) && any(!is.na(df$mplex_channel))

  # ── Colour split (per trace) ────────────────────────────────────────────
  # Reuses assign_analysis_group()'s logic but forced into DCR/fluo mode: when
  # the 2-colour switch is on, the split IS the grouping, so `grouping_mode` has
  # nothing left to decide.
  col_grp <- NULL
  if (two_color_on) {
    .dcr_col <- if (isTRUE(s$use_pooled_dcr) && "dcr_pooled" %in% names(df))
      "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL

    med <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        dcr_med  = if (!is.null(.dcr_col)) stats::median(.data[[.dcr_col]], na.rm = TRUE) else NA_real_,
        fluo_med = if ("fluo" %in% names(df)) stats::median(fluo, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )

    use_fluo <- isTRUE(s$color_by_fluo) && "fluo" %in% names(df)
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
        g1 <- s$dcr_group1
        g2 <- s$dcr_group2
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
  # attach_channel_colors(). See the cache key in prepared_df(), server_smlm.R.
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
#   df     : a frame carrying `channel` (from resolve_channels)
#   colors : channel colours, A then B (smlm_channel_colors(); the colour pickers)
#
# RETURNS
#   df with channel_color filled and attr "channel_table" refreshed. Rows that
#   already carry a colour (multiplex, whose builder writes `mplex_color` into
#   the data) keep it.
# =============================================================================
attach_channel_colors <- function(df, colors = smlm_channel_colors()) {
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
      colors
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
#   scope : "within" or "only:<channel id>" (the cluster_group_scope setting)
channels_in_scope <- function(df, scope = SMLM_DEFAULTS$cluster_group_scope) {
  if (!"channel" %in% names(df)) return(list(all = TRUE, ids = NA_integer_))
  ids <- sort(unique(df$channel[!is.na(df$channel)]))
  scope <- scope %||% SMLM_DEFAULTS$cluster_group_scope
  if (grepl("^only:", scope)) {
    k <- suppressWarnings(as.integer(sub("^only:", "", scope)))
    if (is.finite(k) && k %in% ids) ids <- k
  }
  ids
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

# GOAL: label a de-blinked or Median Cluster id with the channel it belongs to,
# e.g. "B-8", for the hover text.
#
# THE LETTER COMES FROM THE CHANNEL COLUMN, NEVER FROM THE ID. This used to
# decode ids >= 10001 as channel B, from when the two-colour path offset channel
# B by a fixed 10000. apply_clustering_pipeline() now offsets by a running
# maximum (channel B starts at max(A) + 1), so every channel-B cluster was being
# labelled "A-". The id is global and unique across channels, so it is shown
# as-is; only the channel prefix is added.
#
# x  : integer cluster ids (0 / NA = none)
# ch : per-row channel label (resolve_channels()' `channel_label`), or NULL
format_tc_id <- function(x, ch = NULL) {
  none <- is.na(x) | x == 0
  if (is.null(ch) || length(ch) != length(x))
    return(ifelse(none, NA_character_, as.character(x)))
  ifelse(none, NA_character_,
         ifelse(is.na(ch) | !nzchar(as.character(ch)), as.character(x),
                paste0(ch, "-", x)))
}

# =============================================================================
# wrap_groups() - THE UNITS AN ALPHA SHAPE OR ISOSURFACE MAY WRAP
# =============================================================================
#
# GOAL
#   Return one row per group a 3-D wrap is allowed to enclose: a Median Cluster
#   WITHIN ONE CHANNEL. Both the alpha-shape wrap and the isosurface iterate over
#   this, so the rule lives in one place.
#
# WHY THE CHANNEL IS PART OF THE KEY
#   apply_clustering_pipeline() already clusters each channel separately and
#   offsets the ids, so today a `cluster2` id never spans two channels and this
#   key adds no rows. It is keyed on (id, channel) anyway so the guarantee holds
#   at the point of drawing, not only as a consequence of how ids happen to be
#   numbered upstream: a wrap that joined two colours would draw one structure
#   out of molecules that were never shown to belong together. With a single
#   channel the result is identical to grouping on the id alone.
#
# ARGUMENTS
#   df     : points carrying `id_col`, optionally `channel`, `channel_label`,
#            `channel_color`
#   id_col : the grouping column, normally "cluster2"
#
# RETURNS
#   data.frame(id, channel, label, color). `channel` is NA and `label` "" when
#   the data has no channel column.
wrap_groups <- function(df, id_col = "cluster2") {
  empty <- data.frame(id = integer(0), channel = integer(0),
                      label = character(0), color = character(0))
  if (!is.data.frame(df) || !nrow(df) || !id_col %in% names(df)) return(empty)
  ids <- df[[id_col]]
  keep <- !is.na(ids) & is.finite(ids) & ids != 0
  if (!any(keep)) return(empty)
  ch  <- if ("channel" %in% names(df)) df$channel[keep] else rep(NA_integer_, sum(keep))
  lab <- if ("channel_label" %in% names(df)) as.character(df$channel_label[keep])
         else rep("", sum(keep))
  col <- if ("channel_color" %in% names(df)) as.character(df$channel_color[keep])
         else rep(NA_character_, sum(keep))
  g <- data.frame(id = ids[keep], channel = ch, label = lab, color = col,
                  stringsAsFactors = FALSE)
  g <- g[!duplicated(g[, c("id", "channel")]), , drop = FALSE]
  g[order(g$id, g$channel), , drop = FALSE]
}

# GOAL: rows of `df` belonging to one wrap_groups() row. Matches on the channel
# as well as the id, which is what makes the (id, channel) key binding.
wrap_group_rows <- function(df, id_col, g) {
  sel <- !is.na(df[[id_col]]) & df[[id_col]] == g$id
  if (!is.na(g$channel) && "channel" %in% names(df))
    sel <- sel & !is.na(df$channel) & df$channel == g$channel
  df[sel, , drop = FALSE]
}


# =============================================================================
# PANEL HELPERS (moved from server_smlm.R; plain functions of their arguments)
# =============================================================================

# Source key of one spec entry (older specs without $src fall back to the channel).
.mplex_key <- function(c) c$src %||% paste0("ch:", c$channel)

# GOAL: recognise a CSV exported from a multiplex and describe it as one.
#   The export writes mplex_channel / mplex_label / mplex_color per row and the
#   total alignment shift already applied (mplex_dx_nm / _dy_nm / _dz_nm). The
#   shift is recorded in the spec as $prior (display and re-export only; the
#   coordinates already carry it) and the columns are dropped so they cannot
#   go stale. New nudges on the loaded file start from zero.
# RETURNS df, with the multiplex attributes when it is a multiplex CSV.
.mplex_adopt_csv <- function(df, fname) {
  if (!is.data.frame(df) || !all(c("mplex_channel", "mplex_label") %in% names(df))) return(df)
  df$mplex_channel <- suppressWarnings(as.integer(df$mplex_channel))
  df <- df[!is.na(df$mplex_channel), , drop = FALSE]
  if (!nrow(df)) return(df)
  if (!"mplex_color" %in% names(df)) df$mplex_color <- NA_character_
  chans <- sort(unique(df$mplex_channel))
  defs  <- mplex_default_colors(length(chans))
  pcol  <- c("mplex_dx_nm", "mplex_dy_nm", "mplex_dz_nm")
  spec <- lapply(seq_along(chans), function(i) {
    rows <- df$mplex_channel == chans[i]
    col  <- df$mplex_color[rows][1]
    if (is.na(col) || !nzchar(col)) { col <- defs[i]; df$mplex_color[rows] <<- col }
    lab  <- as.character(df$mplex_label[rows][1])
    pr   <- vapply(pcol, function(p)
      if (p %in% names(df)) suppressWarnings(as.numeric(df[[p]][rows][1])) else 0, numeric(1))
    pr[!is.finite(pr)] <- 0
    list(channel = chans[i], label = lab, color = col,
         n_loc = sum(rows), n_trace = dplyr::n_distinct(df$tid[rows]),
         src = paste0("csv:", fname, ":", lab),
         prior = stats::setNames(pr, c("x", "y", "z")))
  })
  df <- df[, setdiff(names(df), pcol), drop = FALSE]
  attr(df, "mplex_spec")       <- spec
  attr(df, "mplex_store_name") <- fname
  attr(df, "mplex_off")        <- list()
  df
}

#
# Covers BOTH channel sources - a 2-colour DCR split and a multiplex build -
# because they are one concept now. The older `mplex_split_active()` above is
# multiplex-only and splits the RAW upload; this one splits a frame that has
# already been through resolve_channels(), so the labels and colours match
# everything else on screen.
#
# THE POOLED PARENT COMES BACK WITH EACH CHANNEL, and that is the point.
# Ripley's L and g(r) divide by the observed area. If each channel derived its
# own window, the sparser channel would get a smaller bounding box and
# footprint, hence a higher apparent density, hence a DEFLATED curve - a
# difference between the channels manufactured out of their difference in
# count. The two dyes imaged the SAME field, so they must be measured in the
# same window; then a real difference in curve is a difference in how the two
# targets are arranged in that field, which is the biology being asked about.
# The plain split: one entry per channel, tagged "file · channel". Used by
# every analysis that does NOT divide by an area. Returns a single untagged
# entry when there is nothing to split, so callers need no special case.
channel_parts <- function(df, tag) {
  none <- list(list(df = df, tag = tag))
  if (!is.data.frame(df) || !nrow(df) || !"channel" %in% names(df)) return(none)
  ids <- sort(unique(df$channel[!is.na(df$channel)]))
  if (length(ids) < 2L) return(none)
  lapply(ids, function(ch) {
    sub <- df[which(df$channel == ch), , drop = FALSE]
    list(df = sub, tag = .facet_name(tag, as.character(sub$channel_label[1])))
  })
}

# The same split, plus the pooled parent, for the analyses that DO divide by an
# area (see the note above).
split_df_by_channel <- function(df, tag) {
  lapply(channel_parts(df, tag), function(pt)
    c(pt, list(win_df  = df,    # pooled file - the shared observation window
               win_tag = tag))) # shared support cache key, so it is built once
}

# Append the channel label to a dataset name, so one series per file becomes one
# series per file-and-channel and every existing grouping/colour path downstream
# keeps working unchanged.
.facet_name <- function(fname, label)
  if (nzchar(label %||% "")) paste0(fname, " · ", label) else fname

get_tc_pts <- function(df, grp, unit) {
  sub <- dplyr::filter(df, analysis_group == grp)
  if (!nrow(sub)) return(data.frame(x = numeric(), y = numeric(), z = numeric()))
  if (unit == "tid") {
    sub %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = median(x, na.rm = TRUE), y = median(y, na.rm = TRUE),
                       z = if ("z" %in% names(sub)) median(z, na.rm = TRUE) else 0,
                       .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y))
  } else {
    if (!"cluster" %in% names(sub) || !any(!is.na(sub$cluster) & sub$cluster > 0L))
      return(data.frame(x = numeric(), y = numeric(), z = numeric()))
    sub %>%
      dplyr::filter(!is.na(cluster), cluster > 0L) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(x = median(x, na.rm = TRUE), y = median(y, na.rm = TRUE),
                       z = if ("z" %in% names(sub)) median(z, na.rm = TRUE) else 0,
                       .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y))
  }
}

