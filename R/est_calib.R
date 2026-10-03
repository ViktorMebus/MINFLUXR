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
##   calib_cutoff() ..... where the pooled immobile MSD has flattened (suggested
##                        lower cutoff for the MSD fit window)
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
                       curve_max = 200L) {
  row <- data.frame(file = file, n_locs = nrow(df), n_tracks = NA_integer_, t_c_ms = NA_real_,
                    eco_med = NA_real_, msd1 = NA_real_, sigma_nm = NA_real_, a = NA_real_,
                    D_w = NA_real_, D_w_short = NA_real_, note = "", stringsAsFactors = FALSE)
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
  qa <- q(tab$a); qd <- q(tab$D_w); qs <- q(tab$D_w_short)
  list(n_cells = sum(is.finite(tab$a)), a = qa[2], a_q25 = qa[1], a_q75 = qa[3],
       D_w = qd[2], Dw_q25 = qd[1], Dw_q75 = qd[3], D_w_short = qs[2])
}

# ── calib_cutoff() ────────────────────────────────────────────────────────────
# GOAL: suggest where the MSD fit window should start for this sequence: the lag
# after which the immobile MSD has largely stopped rising, so the wander adds
# little slope to a fit that starts there. The cells' curves are averaged (each
# cell counts once); the rise per cycle is smoothed over 5 lags and compared with
# its rate over the first three lags. A heuristic, shown on the plot to be judged.
# curves : list of calib_cell() curves; frac : the fraction of the initial rate
# RETURNS list(lag = suggested cutoff, or NA when the curve does not flatten in its
#   first half, pooled = data.frame(lag, msd))
calib_cutoff <- function(curves, frac = 0.2) {
  curves <- Filter(function(c) !is.null(c) && nrow(c) > 0, curves)
  if (!length(curves)) return(list(lag = NA_integer_, pooled = NULL))
  L <- max(vapply(curves, nrow, integer(1)))
  m <- vapply(curves, function(c) { v <- rep(NA_real_, L); ok <- c$keep & is.finite(c$msd)
                                    v[c$lag[ok]] <- c$msd[ok]; v }, numeric(L))
  m <- matrix(m, nrow = L)
  pooled <- rowMeans(m, na.rm = TRUE); pooled[!is.finite(pooled)] <- NA_real_
  rise <- diff(c(NA_real_, pooled))                      # rise per cycle at each lag
  r0   <- mean(rise[2:4], na.rm = TRUE)                  # initial rate (lags 1 -> 4)
  sm   <- stats::filter(rise, rep(1 / 5, 5), sides = 2)  # 5-lag running mean
  cut  <- if (is.finite(r0) && r0 > 0) which(is.finite(sm) & seq_len(L) >= 3L & sm < frac * r0)[1] else NA_integer_
  # A "flattening" found only in the second half of the curve is a noise dip in a
  # curve that keeps rising; report none rather than a misleading cutoff.
  if (length(cut) && is.finite(cut) && cut > L / 2) cut <- NA_integer_
  list(lag = if (length(cut) && is.finite(cut)) as.integer(cut) else NA_integer_,
       pooled = data.frame(lag = seq_len(L), msd = pooled))
}
