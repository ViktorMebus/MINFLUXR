# Univariate Ripley L(r) and pair correlation g(r). The CSR envelopes, global
# p-value and footprint debiasing are covered in test-spatial-envelope.R.
.W   <- 4000
.win <- rbind(c(0, 0), c(.W, .W))

# CLAIM: K(r) = |A| * 2 * P(r) / (n(n-1)), with P(r) the number of point pairs
# closer than r. Checked by hand on four points in a 10 x 10 window.
test_that("Ripley K equals the hand count on a small pattern", {
  x <- c(0, 1, 0, 10); y <- c(0, 0, 1, 10)
  K <- compute_ripley_L(x, y, r_seq = c(1.2, 1.5, 20), window_ranges = rbind(c(0, 0), c(10, 10)))
  expect_equal(K$Kr, 100 * 2 * c(2, 3, 6) / (4 * 3))    # 2, 3, then all 6 pairs
  expect_equal(K$Lr, sqrt(K$Kr / pi))
  expect_equal(K$Lr_minus_r, K$Lr - K$r)
})

# CLAIM: without a footprint there is no edge correction, so L(r) - r sinks below 0
# as r grows; with a footprint the translation correction removes that bias.
test_that("edge correction brings L(r) - r to 0 on CSR", {
  set.seed(1)
  p   <- data.frame(x = stats::runif(1500, 0, .W), y = stats::runif(1500, 0, .W))
  r   <- seq(25, 500, by = 25)
  raw <- compute_ripley_L(p$x, p$y, r_seq = r, window_ranges = .win, max_subsample = 5000)
  sup <- estimate_support(cbind(p, z = 0), mode = "grid", h_nm = 100, close_nm = 200,
                          return_mask = TRUE)
  sup <- attach_gamma_cache(sup, r_max = 600)
  cor <- compute_ripley_L(p$x, p$y, r_seq = r, window_ranges = .win, support = sup,
                          max_subsample = 5000)
  expect_lt(utils::tail(raw$Lr_minus_r, 1), -15)        # biased low at r = 500
  expect_lt(max(abs(cor$Lr_minus_r)), 8)                # corrected: within a few nm
})

test_that("3D L(r) - r is near 0 on CSR at short range", {
  set.seed(4)
  p3 <- matrix(stats::runif(3 * 1500, 0, 2000), ncol = 3)
  L3 <- compute_ripley_L(p3[, 1], p3[, 2], p3[, 3], r_seq = seq(20, 100, by = 20),
                         use_3d = TRUE, window_ranges = rbind(c(0, 0, 0), c(2000, 2000, 2000)),
                         max_subsample = 5000)
  expect_equal(L3$Lr, (3 * L3$Kr / (4 * pi))^(1 / 3))
  expect_lt(max(abs(L3$Lr_minus_r)), 5)
})

test_that("g(r) is 1 on CSR with the footprint correction", {
  set.seed(2)
  p   <- data.frame(x = stats::runif(1500, 0, .W), y = stats::runif(1500, 0, .W))
  sup <- attach_gamma_cache(estimate_support(cbind(p, z = 0), mode = "grid", h_nm = 100,
                                             close_nm = 200, return_mask = TRUE), r_max = 600)
  g <- compute_pcf(p$x, p$y, r_seq = seq(20, 400, by = 20), window_ranges = .win,
                   support = sup, max_subsample = 5000)
  expect_true(attr(g, "edge_corrected"))
  expect_equal(mean(g$g), 1, tolerance = 0.05)
})

# CLAIM: g(r) is the local form of K, so clusters of size ~s give a peak at short r
# that returns to 1 beyond the cluster scale.
test_that("g(r) peaks inside clusters and returns to 1 outside them", {
  set.seed(3)
  P <- sim_thomas(60, 25, 40, .W)
  g <- compute_pcf(P[, 1], P[, 2], r_seq = seq(10, 600, by = 10), window_ranges = .win,
                   max_subsample = 5000)
  expect_gt(g$g[1], 5)                                  # r = 10 nm
  expect_equal(mean(g$g[g$r >= 300]), 1, tolerance = 0.2)
})

# CLAIM: the subsample above the point cap is seeded, so two runs agree and the
# global RNG stream is left as it was.
test_that("subsampled statistics are reproducible and RNG-neutral", {
  set.seed(5)
  p <- data.frame(x = stats::runif(3000, 0, .W), y = stats::runif(3000, 0, .W))
  set.seed(9); a <- stats::runif(1)
  set.seed(9)
  L1 <- compute_ripley_L(p$x, p$y, r_seq = seq(50, 300, by = 50), max_subsample = 1000)
  b  <- stats::runif(1)
  L2 <- compute_ripley_L(p$x, p$y, r_seq = seq(50, 300, by = 50), max_subsample = 1000)
  expect_identical(L1, L2)
  expect_identical(a, b)
})

# CLAIM: the correlation length is fitted by nls on g itself, which recovers xi on a
# noisy exponential where the log-linear fit is biased high.
test_that("fit_correlation_length recovers xi on a noisy exponential", {
  set.seed(6)
  r <- seq(5, 600, by = 5)
  for (xi in c(40, 120, 250)) {
    est <- replicate(10, {
      g <- 1 + 3 * exp(-r / xi) + stats::rnorm(length(r), 0, 0.05)
      fit_correlation_length(r, g, r_min = 10)$xi
    })
    expect_equal(stats::median(est), xi, tolerance = 0.12)
  }
  flat <- fit_correlation_length(r, 1 + stats::rnorm(length(r), 0, 0.05))
  expect_true(is.na(flat$xi) || flat$xi_se > flat$xi)   # no decay, no trustworthy xi
})

# CLAIM: the panels' per-file computation (spatial_curve_for_points) is exactly the
# estimator on the given window, plus an envelope and global p when asked.
test_that("spatial_curve_for_points gives the estimator's curve and the envelope", {
  set.seed(7)
  p   <- data.frame(x = stats::runif(400, 0, 4000), y = stats::runif(400, 0, 4000), z = 0)
  win <- rbind(c(0, 0), c(4000, 4000)); r <- seq(50, 400, by = 50)
  L <- spatial_curve_for_points(p, win, "L", r_seq = r, tag = "a")
  expect_equal(L$curve$Lr_minus_r,
               compute_ripley_L(p$x, p$y, r_seq = r, window_ranges = win)$Lr_minus_r)
  expect_true(all(L$curve$dataset == "a")); expect_null(L$env)
  g <- spatial_curve_for_points(p, win, "g", r_seq = r, tag = "a")
  expect_equal(g$curve$g, compute_pcf(p$x, p$y, r_seq = r, window_ranges = win)$g)
  expect_true(all(g$curve$n_units == 400)); expect_false("n_units" %in% names(L$curve))
  set.seed(8)
  e <- spatial_curve_for_points(p, win, "L", r_seq = r, envelope = TRUE, nsim = 19, tag = "a")
  expect_true(all(c("lo", "hi", "global_p", "dataset") %in% names(e$env)))
  expect_true(e$env$global_p[1] >= 0 && e$env$global_p[1] <= 1)
  expect_null(spatial_curve_for_points(p[1:2, ], win, "L", r_seq = r))
})
