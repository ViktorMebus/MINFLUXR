# Tests for est_hmm.R (continuous-time diffusive-state HMM).
source(file.path(.repo_root, "R/est_hmm.R"))

test_that(".expm_small matches the closed form for a two-state generator", {
  a <- 0.3; b <- 0.1; t <- 2.5
  Q <- matrix(c(-a, a, b, -b), 2, byrow = TRUE)
  e <- exp(-(a + b) * t); p1 <- b / (a + b); p2 <- a / (a + b)
  P_exact <- matrix(c(p1 + p2 * e, p2 * (1 - e), p1 * (1 - e), p2 + p1 * e), 2, byrow = TRUE)
  expect_equal(.expm_small(Q * t), P_exact, tolerance = 1e-10)
  expect_equal(.expm_small(Q * 500), matrix(c(p1, p2), 2, 2, byrow = TRUE), tolerance = 1e-8)
  expect_equal(.ctmc_stationary(Q), c(p1, p2), tolerance = 1e-12)
})

test_that("parameter packing round-trips and keeps D ordered", {
  D <- c(20, 150, 900); Q <- matrix(c(0, .1, .2, .3, 0, .4, .5, .6, 0), 3, byrow = TRUE)
  diag(Q) <- -rowSums(Q)
  u <- .hmm_unpack(.hmm_pack(D, Q), 3)
  expect_equal(u$D, D, tolerance = 1e-10); expect_equal(u$Q, Q, tolerance = 1e-10)
  expect_true(all(diff(.hmm_unpack(rnorm(9), 3)$D) > 0))
})

test_that("one-state limit: the emission likelihood is maximised near the true D", {
  set.seed(4)
  Q <- matrix(c(-1e-6, 1e-6, 1e-6, -1e-6), 2)
  df <- sim_ct_state_tracks(60, 200, D = c(200, 200.0001), Q = Q, sigma = 4, seed = 4, R_blur = 0.106)
  pr <- hmm_ct_prepare(df, sigma = 0)
  R_abs <- 0.192 * 0.106                          # the R the tracks were simulated with
  nll <- function(lD) -sum(.hmm_emissions(pr, exp(lD), R_abs, log = TRUE)[!is.na(pr$R2)])
  D_hat <- exp(optimize(nll, c(log(10), log(5000)))$minimum)
  expect_lt(abs(D_hat / 200 - 1), 0.05)          # blur + noise terms are right
})

# CLAIM (methods doc §I-B.5): at 30-cycle dwells and 10-fold D contrast the MLE
# recovers D and dwell times within 15%, where decoded-path counting reads dwells
# more than 40% long. Full map: tests/validation/hmm_recovery_sweep.R.
test_that("two-state recovery: MLE unbiased, decoded-path counting biased long", {
  skip_if_not_slow()
  t_c <- 0.192; dw <- 30 * t_c
  D <- c(30, 300); Q <- matrix(c(0, 1 / dw, 1 / dw, 0), 2); diag(Q) <- -rowSums(Q)
  set.seed(1); eco <- pmax(20, rlnorm(1e5, log(36), 0.6))
  df <- sim_ct_state_tracks(200, 300, D, Q, sigma = 100.8 / (4 * sqrt(eco)), seed = 11)
  pr <- hmm_ct_prepare(df, sigma = 0)
  f  <- hmm_ct_fit(pr, K = 2, R_abs = t_c * 0.106, t_c = t_c)   # the simulator default R
  expect_true(all(abs(f$D / D - 1) < 0.15))
  expect_true(all(abs(f$dwell_ms / dw - 1) < 0.15))
  vs  <- hmm_ct_viterbi(pr, f)
  expect_gt(mean(vs == pr$df$state[pr$steps$row]), 0.8)
  old <- hmm_path_rates(pr, vs, 2, t_c)
  expect_true(all(old$dwell_ms / dw - 1 > 0.4))
})
