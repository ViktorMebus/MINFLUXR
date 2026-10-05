## ============================================================================
## est_calib.R - instrument calibration from an immobile sample
## ----------------------------------------------------------------------------
## The photon-budget model sigma = L / (a sqrt(N)) has one free constant, a. The
## ideal-estimator value follows from the pattern geometry (TCP_PATTERN_GEOM), but
## the positions an instrument exports carry its own estimator, background and
## tracking loop, so a has to be measured for each acquisition sequence on that
## instrument. On an immobile sample all apparent motion is instrumental:
##
##     MSD(k steps) = 4 sigma^2 + 4 D_w * k * step      (lateral, short lags)
##
## The intercept gives sigma, the slope the wander floor D_w, and
##     a = L * sqrt(<1/N>) / sigma
## with <1/N> the mean inverse photon count, since sigma^2 = L^2 <1/N> / a^2.
##
##   calib_cell() ....... one immobile file -> sigma, a, D_w (one row) + its MSD curve
##   calib_summary() .... median and quartiles across cells -> the profile values
##   calib_a_at() ....... a recomputed for another blur coefficient (sensitivity)
##
## Two wander floors are reported. D_w is the slope of the immobile MSD over the
## MSD fit window, the floor that matches the headline (upper-regime) D. D_w_short
## is the slope over the first few steps, where the wander rises fastest.
## ============================================================================

# ── calib_cell() ──────────────────────────────────────────────────────────────
# GOAL: measure the localization error, prefactor a and wander floor of ONE
# immobile-sample file, after the same filter pipeline the live data gets.
# df      : localizations (tid, x, y, tim in ms, eco)
# file    : label for the row
# L       : TCP diameter of the sequence (nm)
# filters : filter settings as filter_pipeline() takes them, without t_c / sw
# lags    : step lags fitted; short lags keep the fit inside the linear regime
# RETURNS one-row data.frame(file, n_locs, n_tracks, t_c_ms, eco_med, msd1, sigma_nm,
#   a, D_w, note). sigma and a are NA when the fitted intercept is not positive.
calib_cell <- function(df, file, L, filters = list(), lags = 1:5,
                       msd = list(window = c(5L, 50L), track_weight = TRUE, min_pts = 50L),
                       curve_max = 200L, R = NULL, n_vert = 6L) {
  row <- data.frame(file = file, n_locs = nrow(df), n_tracks = NA_integer_, t_c_ms = NA_real_,
                    eco_med = NA_real_, msd1 = NA_real_, sigma_nm = NA_real_, a = NA_real_,
                    D_w = NA_real_, D_w_short = NA_real_, R = NA_real_, msd0 = NA_real_,
                    slope1 = NA_real_, inv_N = NA_real_, note = "", stringsAsFactors = FALSE)
  if (!nrow(df) || !all(c("tid", "x", "y", "tim", "eco") %in% names(df))) {
    row$note <- "needs tid, x, y, tim and eco"; return(row)
  }
  cc <- cycle_cluster(df$tid, df$tim)
  if (is.null(cc) || !is.finite(cc$cycle)) { row$note <- "too few intervals for a cycle time"; return(row) }
  row$t_c_ms <- cc$cycle
  fp <- filters; fp$t_c <- cc$cycle; fp$sw <- cc$sw
  df <- filter_pipeline(df, fp)
  if (!nrow(df)) { row$note <- "nothing left after the filters"; return(row) }
  df <- df[order(df$tid, df$tim), ]
  row$n_tracks <- length(unique(df$tid))
  eco <- df$eco[is.finite(df$eco) & df$eco > 0]
  row$eco_med <- stats::median(eco)

  # MSD over exactly k consecutive single-cycle steps: each spans t_c + SW
  step <- cc$cycle + cc$sw
  n <- nrow(df)
  msd_s <- vapply(lags, function(k) {
    m <- n - k; if (m < 2L) return(NA_real_)
    same <- df$tid[(k + 1):n] == df$tid[1:m]
    gap  <- df$tim[(k + 1):n] - df$tim[1:m]
    ok   <- same & abs(gap - k * step) < 0.25 * step
    if (sum(ok) < 50L) return(NA_real_)
    mean((df$x[(k + 1):n][ok] - df$x[1:m][ok])^2 + (df$y[(k + 1):n][ok] - df$y[1:m][ok])^2)
  }, numeric(1))
  row$msd1 <- msd_s[1]

  # The immobile ensemble MSD over many cycle lags (lateral), as the MSD panel
  # computes it; its slope over the fit window is the floor under the headline D.
  win <- as.integer(msd$window %||% c(5L, 50L))
  sums <- msd_track_sums(df, cc$cycle, cc$sw, max(as.integer(curve_max), win[2]),
                         as.integer(msd$min_pts %||% 50L), has_z = FALSE)
  ens <- msd_ensemble(sums, cc$cycle, cc$sw, isTRUE(msd$track_weight %||% TRUE))
  curve <- data.frame(lag = ens$lag, tau = ens$tau_lin, msd = ens$msd, keep = ens$keep_lag)
  attr(row, "curve") <- curve
  okw <- curve$keep & is.finite(curve$msd) & curve$lag >= win[1] & curve$lag <= win[2]
  if (sum(okw) >= 3L)
    row$D_w <- unname(stats::coef(stats::lm(curve$msd[okw] ~ curve$tau[okw]))[2]) / 4   # nm^2/ms

  okl <- is.finite(msd_s)
  if (sum(okl) < 3L) { row$note <- "too few single-cycle runs for the short-lag fit"; return(row) }
  cf <- unname(stats::coef(stats::lm(msd_s[okl] ~ lags[okl])))
  row$D_w_short <- cf[2] / (4 * step)                    # nm^2/ms

  # Motion blur. A position is averaged over its photon collection, so for any
  # residual motion the short-lag MSD is 4 sigma^2 + 4 D (k - 2R) step (Berglund 2010):
  # the intercept sits 2R x slope BELOW 4 sigma^2. Immobile samples still move at
  # this scale (the step correlation comes out positive, the blur signature), so
  # leaving it out reads sigma low and a high. R per step from this file's photon
  # collection window eco/efo, as the main tab does without a sequence file.
  if (is.null(R)) {
    t_coll <- if ("efo" %in% names(df)) stats::median(df$eco / df$efo, na.rm = TRUE) * 1000 else NA_real_
    br <- if (is.finite(t_coll) && t_coll > 0) blur_R_step(step, t_coll, n_vert = n_vert) else NULL
    R  <- if (!is.null(br) && is.finite(br$R)) br$R else BLUR_R_STEP_FALLBACK
  }
  row$R <- R; row$msd0 <- cf[1]; row$slope1 <- cf[2]; row$inv_N <- mean(1 / eco)
  s2 <- (cf[1] + 2 * R * cf[2]) / 4                      # per-axis sigma^2, lateral
  if (!is.finite(s2) || s2 <= 0) { row$note <- "fitted intercept not positive"; return(row) }
  row$sigma_nm <- sqrt(s2)
  row$a <- L * sqrt(row$inv_N) / row$sigma_nm
  row
}

# ── calib_a_at() ──────────────────────────────────────────────────────────────
# GOAL: recompute each cell's a for another blur coefficient (R x f), from the
# stored short-lag fit, to show how much the result leans on R.
calib_a_at <- function(tab, f) {
  s2 <- (tab$msd0 + 2 * tab$R * f * tab$slope1) / 4
  ifelse(is.finite(s2) & s2 > 0, tab$a * tab$sigma_nm / sqrt(s2), NA_real_)
}

# ── calib_summary() ───────────────────────────────────────────────────────────
# GOAL: the profile values from the per-cell rows: median and quartiles of a and
# D_w across cells. The median is the profile value; the quartiles show how well
# the cells agree.
calib_summary <- function(tab) {
  q <- function(v) { v <- v[is.finite(v)]
    if (!length(v)) c(NA_real_, NA_real_, NA_real_)
    else unname(stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)) }
  qa <- q(tab$a); qd <- q(tab$D_w); qs <- q(tab$D_w_short)
  med <- function(v) { v <- v[is.finite(v)]; if (length(v)) stats::median(v) else NA_real_ }
  has_fit <- all(c("msd0", "slope1", "R") %in% names(tab))
  list(n_cells = sum(is.finite(tab$a)), a = qa[2], a_q25 = qa[1], a_q75 = qa[3],
       D_w = qd[2], Dw_q25 = qd[1], Dw_q75 = qd[3], D_w_short = qs[2],
       R = if (has_fit) med(tab$R) else NA_real_,
       a_R_lo = if (has_fit) med(calib_a_at(tab, 0.7)) else NA_real_,   # R 30% lower
       a_R_hi = if (has_fit) med(calib_a_at(tab, 1.3)) else NA_real_)   # R 30% higher
}
