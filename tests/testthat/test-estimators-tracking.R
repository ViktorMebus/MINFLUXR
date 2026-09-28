# Recovery tests for the tracking estimators extracted from the server
# (est_cve.R, est_msd.R, est_jdd.R, est_dmap.R, est_cpd.R), against ground truth
# from sim_ct_state_tracks() (est_hmm.R) at tracking-like acquisition settings:
# t_c 0.192 ms, SW 0.011 ms, dark-cycle mix, per-localization photon σ, blur R from the
# timing model of the standard 2D sequence (MSD_R_FALLBACK, est_blur.R). These tests check the
# algebra at a known R; test-blur.R checks that R against a simulated MINFLUX estimator.
TCP_CRLB_PREFACTOR <- if (exists("TCP_CRLB_PREFACTOR")) TCP_CRLB_PREFACTOR else 4
for (f in c("R/est_blur.R", "R/est_hmm.R", "R/est_cve.R", "R/est_msd.R", "R/est_jdd.R", "R/est_dmap.R", "R/est_cpd.R"))
  source(file.path(.repo_root, f))

.t_c <- 0.192; .sw <- 0.011; .R <- MSD_R_FALLBACK
.free_tracks <- function(D, n_tracks = 200, n_locs = 300, seed = 1, gap_p = c(0.72, 0.2, 0.05, 0.03)) {
  set.seed(seed)
  eco <- pmax(20, rlnorm(1e5, log(36), 0.6))
  sig <- 100.8 / (4 * sqrt(eco))
  Q0  <- matrix(c(-1e-9, 1e-9, 1e-9, -1e-9), 2)          # one effective state
  df  <- sim_ct_state_tracks(n_tracks, n_locs, D = c(D, D * 1.000001), Q = Q0, t_c = .t_c, sw = .sw,
                             gap_p = gap_p, sigma = sig, R_blur = .R, seed = seed)
  attr(df, "sig_rms") <- sqrt(mean(sig^2))
  df
}

# CLAIM (methods doc §I-B.1): the CVE D is blur-independent and unbiased.
test_that("CVE recovers D on blurred, noisy, irregularly sampled tracks", {
  skip_if_not_slow()
  for (D in c(50, 215, 1000)) {
    df <- .free_tracks(D, seed = 21)
    ce <- cve_estimate(df, .t_c)
    expect_lt(abs(ce$D_cve / D - 1), 0.05)
  }
})

# CLAIM (§I-B MSD): OLSF with the blur term recovers σ; R = 0 under-reads it.
test_that("OLSF blur term recovers sigma from the MSD intercept", {
  set.seed(1)
  D <- 215; sig <- 5; dt <- 0.192; nsub <- 24; N <- 1e5
  # Each localization averages the path over a window of 6*R*dt at the start of the
  # cycle, which is exactly blur R (Berglund 2010), so the fit is tested at a known R.
  n_win <- round(6 * .R * nsub)
  sim1 <- function() { h <- dt / nsub
    path <- cumsum(rnorm(N * nsub, 0, sqrt(2 * D * h)))
    rowMeans(matrix(path, ncol = nsub, byrow = TRUE)[, seq_len(n_win), drop = FALSE]) + rnorm(N, 0, sig) }
  x <- sim1(); y <- sim1(); lags <- 1:50
  msd <- sapply(lags, function(k) mean((x[-(1:k)] - x[1:(N - k)])^2 + (y[-(1:k)] - y[1:(N - k)])^2))
  tau <- lags * dt
  f0 <- msd_fit_olsf(tau, msd, rep(1L, 50), k = 50, auto_pmin = TRUE, R = 0, Npts = 1e4, dt_cycle = dt)
  fb <- msd_fit_olsf(tau, msd, rep(1L, 50), k = 50, auto_pmin = TRUE, R = .R, Npts = 1e4, dt_cycle = dt)
  fi <- msd_fit_olsf_iterative(tau, msd, R = .R, dt_cycle = dt, N = 50)
  expect_lt(abs(fb$D / D - 1), 0.03)
  expect_lt(abs(fb$sigma_loc / sig - 1), 0.05)
  expect_lt(abs(fi$sigma_loc / sig - 1), 0.05)
  expect_lt(f0$sigma_loc / sig, 0.9)                        # R = 0 reads ~4.0 nm at the physical R
})

# CLAIM (§I-B.4, fixed 2026-09-17): JDD D is noise- and blur-corrected in both modes.
test_that("JDD recovers D with blur, single-cycle and pooled", {
  skip_if_not_slow()
  for (D in c(50, 1000)) {
    df <- .free_tracks(D, n_tracks = 300, seed = 3)
    for (io in c(TRUE, FALSE)) {
      cj <- jdd_collect_jumps(df, min_pts = 10, integ_only = io, t_c = .t_c, norm_dt = TRUE)
      jf <- jdd_fit(cj$jumps, cj$dt_all, n_comp = 1, sigma_loc = attr(df, "sig_rms"), R_abs = .R * .t_c)
      expect_lt(abs(jf$D_vals / D - 1), 0.05)
    }
  }
})

# CLAIM (§I-B.9, fixed 2026-09-17): the D-map is unbiased on a uniform blurred field.
test_that("D-map recovers a uniform D with blur", {
  skip_if_not_slow()
  df <- .free_tracks(215, n_tracks = 300, seed = 6)
  m  <- dmap_compute(df, 5L, 50L, .t_c, TRUE, attr(df, "sig_rms"), R_abs = .R * .t_c)
  w  <- m$N_mat > 0
  expect_lt(abs(weighted.mean(m$D_mat[w], m$N_mat[w]) / 215 - 1), 0.05)
  m0 <- dmap_compute(df, 5L, 50L, .t_c, TRUE, attr(df, "sig_rms"), R_abs = 0)
  expect_lt(weighted.mean(m0$D_mat[w], m0$N_mat[w]) / 215, 0.88)  # unblurred form ~ −20% at the physical R
})

# CLAIM (§I-B.3): scaled PELT cuts few single-state tracks, and segment D is unbiased.
test_that("change-point: few spurious segments on one free population, D unbiased", {
  skip_if_not_slow()
  df   <- .free_tracks(215, n_tracks = 150, n_locs = 300, seed = 7, gap_p = 1)
  segs <- cpd_segment_tracks(df, min_seg = 8L, noise_sig2 = attr(df, "sig_rms")^2,
                             R_abs = .R * .t_c)
  frac_cut <- mean(tapply(segs$n_states, segs$tid, max) > 1)
  expect_lt(frac_cut, 0.15)
  expect_lt(abs(stats::weighted.mean(segs$D_nm2_ms, segs$n_steps) / 215 - 1), 0.05)
})

test_that("photon CRLB reduces to L/(a sqrt N) at nu = 0 and x = 0", {
  expect_equal(.crlb_doughnut(100.8, 0, c(36, 100), 0, 4), 100.8 / (4 * sqrt(c(36, 100))))
  expect_equal(.crlb_doughnut(100.8, 0, 36, 10, 4), 100.8 / (4 * 6) * (1 + 2 * (10 / 100.8)^2))
})
