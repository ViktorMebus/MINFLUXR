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
##   calib_cell() ....... one immobile file -> sigma, a, D_w (one row)
##   calib_summary() .... median and quartiles across cells -> the profile values
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
calib_cell <- function(df, file, L, filters = list(), lags = 1:5) {
  row <- data.frame(file = file, n_locs = nrow(df), n_tracks = NA_integer_, t_c_ms = NA_real_,
                    eco_med = NA_real_, msd1 = NA_real_, sigma_nm = NA_real_, a = NA_real_,
                    D_w = NA_real_, note = "", stringsAsFactors = FALSE)
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
  msd <- vapply(lags, function(k) {
    m <- n - k; if (m < 2L) return(NA_real_)
    same <- df$tid[(k + 1):n] == df$tid[1:m]
    gap  <- df$tim[(k + 1):n] - df$tim[1:m]
    ok   <- same & abs(gap - k * step) < 0.25 * step
    if (sum(ok) < 50L) return(NA_real_)
    mean((df$x[(k + 1):n][ok] - df$x[1:m][ok])^2 + (df$y[(k + 1):n][ok] - df$y[1:m][ok])^2)
  }, numeric(1))
  row$msd1 <- msd[1]
  okl <- is.finite(msd)
  if (sum(okl) < 3L) { row$note <- "too few single-cycle runs for the short-lag fit"; return(row) }
  cf <- unname(stats::coef(stats::lm(msd[okl] ~ lags[okl])))
  row$D_w <- cf[2] / (4 * step)                          # nm^2/ms
  if (!is.finite(cf[1]) || cf[1] <= 0) { row$note <- "fitted intercept not positive"; return(row) }
  sigma <- sqrt(cf[1] / 4)                               # per-axis, lateral
  row$sigma_nm <- sigma
  row$a <- L * sqrt(mean(1 / eco)) / sigma
  row
}

# ── calib_summary() ───────────────────────────────────────────────────────────
# GOAL: the profile values from the per-cell rows: median and quartiles of a and
# D_w across cells. The median is the profile value; the quartiles show how well
# the cells agree.
calib_summary <- function(tab) {
  q <- function(v) { v <- v[is.finite(v)]
    if (!length(v)) c(NA_real_, NA_real_, NA_real_)
    else unname(stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)) }
  qa <- q(tab$a); qd <- q(tab$D_w)
  list(n_cells = sum(is.finite(tab$a)), a = qa[2], a_q25 = qa[1], a_q75 = qa[3],
       D_w = qd[2], Dw_q25 = qd[1], Dw_q75 = qd[3])
}
