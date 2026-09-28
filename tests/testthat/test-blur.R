# Motion-blur coefficient R from the sequence timing (R/est_blur.R).
source(file.path(.repo_root, "R/est_blur.R"))

# GP67H tracking sequence: hexagon, 100 us dwell, CFR centre 16 us, step 202.9 us.
.gp_R <- function(...) blur_R_step(0.2029, 0.100, pattern = "hexagon", cfr = TRUE, ctr_factor = 0.16, ...)

# CLAIM (est_blur.R): the timeline model reproduces the per-axis values of the
# MINFLUX projection model (x 0.103, y 0.097) and the equal-weight staircase (0.106).
test_that("blur_R_step reproduces the per-axis values for the GP67H sequence", {
  b <- .gp_R()
  expect_equal(unname(b$R_axis["x"]), 0.103, tolerance = 0.005)
  expect_equal(unname(b$R_axis["y"]), 0.097, tolerance = 0.005)
  expect_equal(b$R_uniform, 0.106, tolerance = 0.005)
  expect_false(b$overrun)
})

# CLAIM: R follows the sequence timing - shorter dwell and fewer vertices lower it,
# pattern repeats and extra passes spread the collection and raise it - so no single
# constant can stand for all sequences.
test_that("R depends on dwell, vertices, repeats and passes", {
  base <- .gp_R()$R
  expect_lt(blur_R_step(0.12, 0.050, pattern = "triangle")$R, base)          # fast triangle sequence
  expect_gt(.gp_R(pat_repeat = 2)$R, base)
  expect_gt(blur_R_step(0.40, 0.100, pattern = "hexagon", eta = 2, cfr = TRUE, ctr_factor = 0.16)$R, base)
  expect_lt(abs(base - 0.10), 0.01)                                          # the documented fallback
  expect_equal(BLUR_R_STEP_FALLBACK, 0.10)
})

test_that("blur_offsets converts R_step to the estimator forms", {
  o <- blur_offsets(0.10, dt_ms = 0.2029, t_c = 0.1914)
  expect_equal(o$R_abs, 0.02029)
  expect_equal(o$R_msd, 0.02029 / 0.1914)
  expect_equal(blur_offsets(NA, 0.2, 0.19)$R_abs, 0)
})

# CLAIM (non-circular): on tracks produced by a MINFLUX estimator - photon counts at
# the six vertices of a parabolic doughnut while the molecule diffuses, position from
# the linearised MINFLUX estimator - the single-step D corrected with the timing-model
# R recovers the true D, while the constant 1/6.2 (0.152 per step) reads it high.
test_that("timing-model R recovers D from simulated MINFLUX localizations; 1/6.2 does not", {
  skip_if_not_slow()
  set.seed(5)
  dt <- 202.9; K <- 6; dwell <- 100 / K; eod <- 5; Lh <- 50.4; D <- 0.215      # us, nm^2/us
  ang <- (0:(K - 1)) * pi / 3; vx <- Lh * cos(ang); vy <- Lh * sin(ang)
  starts <- cumsum(c(eod, rep(dwell + eod, K - 1)))
  h <- 0.25; per <- round(dt / h); tt <- (seq_len(per) - 1) * h
  idx <- lapply(seq_len(K), function(i) which(tt >= starts[i] & tt < starts[i] + dwell))
  n <- 12000; est <- matrix(NA_real_, n, 2); a <- c(0, 0)
  for (k in seq_len(n)) {
    inc <- cbind(cumsum(rnorm(per, 0, sqrt(2 * D * h))), cumsum(rnorm(per, 0, sqrt(2 * D * h))))
    lam <- vapply(seq_len(K), function(i) mean((inc[idx[[i]], 1] - vx[i])^2 + (inc[idx[[i]], 2] - vy[i])^2), 0)
    p <- lam / sum(lam)
    est[k, ] <- a - c(sum((K * p - 1) * vx), sum((K * p - 1) * vy)) / K        # linearised estimator
    a <- a + inc[per, ]
  }
  d2 <- rowSums(diff(est)^2)                                                   # 2D squared steps, one dt
  single_step_D <- function(R_step) mean(d2) / (4 * dt) / (1 - 2 * R_step)     # no static noise here
  R_model <- blur_R_step(dt / 1000, 0.100, pattern = "hexagon", cfr = TRUE, ctr_factor = 0.16)$R
  expect_lt(abs(single_step_D(R_model) / D - 1), 0.05)
  expect_gt(single_step_D((1 / 6.2) * 191.4 / dt) / D - 1, 0.10)
})

# CLAIM: sequences that cannot describe the lateral blur of the data are recognisable,
# so the app can fall back to the data instead of reporting a wrong R: an axial-only
# final pattern gives no lateral R, and a collection timeline longer than the measured
# step (e.g. an imaging sequence loaded with tracking data) is flagged as an overrun.
test_that("axial patterns give no lateral R and mismatched sequences are flagged", {
  z <- blur_R_step(0.2029, 1.0, pattern = "zline2", pat_repeat = 5)
  expect_false(is.finite(z$R))
  im <- blur_R_step(0.2029, 1.0, pattern = "hexagon", pat_repeat = 5, cfr = TRUE, ctr_factor = 1)
  expect_true(im$overrun)
  expect_false(.gp_R()$overrun)
  expect_null(blur_R_step(0.2029, NA))                  # no dwell: no model, caller uses the data
  expect_true(is.finite(blur_R_step(0.2029, 0.1, pattern = "pentagram", n_vert = NA)$R))  # unknown name
})
