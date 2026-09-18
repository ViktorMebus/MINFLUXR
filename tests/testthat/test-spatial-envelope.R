# CSR envelopes for the spatial-statistics panels. Full coverage map:
# tests/validation/spatial_envelope_coverage.R (60 reps x 4 panels x 99 sims).
suppressMessages(suppressWarnings(source(file.path(.repo_root, "packages.R"))))

.csr_L <- function(n, seed = NULL) {                    # 8 x 6 um cell, one quadrant empty
  if (!is.null(seed)) set.seed(seed)
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    x <- runif(3 * n, 0, 8000); y <- runif(3 * n, 0, 6000)
    k <- !(x > 4000 & y > 3000)
    out <- rbind(out, data.frame(x = x[k], y = y[k]))
  }
  o <- out[seq_len(n), ]; o$z <- 0; o
}

test_that("the global band and p-value are self-consistent", {
  set.seed(2)
  p <- .csr_L(400, seed = 2)
  sup <- estimate_support(p, mode = "grid", h_nm = 100, close_nm = 400, return_mask = TRUE)
  win <- apply(as.matrix(p[, c("x", "y")]), 2, range)
  r_seq <- seq(50, 800, by = 50)
  env <- compute_ripley_envelope(n = 400, window_ranges = win, r_seq = r_seq, nsim = 39,
                                 estimator = "L", support = sup)
  expect_true(all(c("g_lo", "g_hi", "null_sd") %in% names(env)))
  expect_true(all(env$g_lo <= env$lo + 1e-9))          # global band is wider than pointwise
  expect_true(all(env$g_hi >= env$hi - 1e-9))
  # a curve ON the null mean cannot be significant; a wildly displaced one must be
  expect_equal(envelope_global_p(env, env$null_mean)$p, 1, tolerance = 1e-9)
  expect_lt(envelope_global_p(env, env$null_mean + 50 * env$null_sd)$p, 0.03)
})

# CLAIM (methods doc §I-A.3): the footprint is grown back to the region the points
# fill, because the occupancy grid stops inside the true cell edge at sparse density.
test_that("debias_support_window recovers the true area of a sparse cell", {
  skip_if_not_slow()
  true_area <- 8000 * 6000 - 4000 * 3000
  for (n in c(200, 500)) {
    r <- replicate(4, {
      p <- .csr_L(n)
      sup <- estimate_support(p, mode = "grid", h_nm = 100, close_nm = 400, return_mask = TRUE)
      db  <- debias_support_window(sup, nrow(p))
      c(sup$support / true_area, db$support / true_area)
    })
    expect_lt(mean(r[1, ]), 0.96)                       # raw estimate is biased low
    expect_lt(abs(mean(r[2, ]) - 1), 0.08)              # debiased is close to the truth
  }
})

# CLAIM: with independent channels in one cell, the random-labelling null keeps the
# global test near its nominal level, where the CSR-in-window null does not
# (measured 1.7% vs 21.7% for L12 over 60 replicates).
test_that("random-labelling null is calibrated on independent channels", {
  skip_if_not_slow()
  set.seed(4)
  sig <- replicate(20, {
    a <- .csr_L(400); b <- .csr_L(300)
    pooled <- rbind(a, b)
    sup <- estimate_support(pooled, mode = "grid", h_nm = 100, close_nm = 400, return_mask = TRUE)
    sup <- debias_support_window(sup, nrow(pooled))
    win <- apply(as.matrix(pooled[, c("x", "y")]), 2, range)
    r_seq <- seq(50, 800, by = 50)
    obs <- compute_ripley_L12(a$x, a$y, NULL, b$x, b$y, NULL, r_seq = r_seq,
                              window_ranges = win, support = sup)$Lr_minus_r
    env <- compute_ripley_envelope(n = 400, n2 = 300, window_ranges = win, r_seq = r_seq,
                                   nsim = 99, estimator = "L12", support = sup,
                                   null_mode = "labels", pts1 = a, pts2 = b)
    isTRUE(envelope_global_p(env, obs)$p < 0.05)
  })
  expect_lt(mean(sig), 0.20)                            # 20 reps: nominal 5%, allow noise
})
