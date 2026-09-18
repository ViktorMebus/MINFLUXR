# CLAIM: the confinement-index null reproduces the free-track CI distribution when
# real tracks are motion-blurred and carry a spread of photon-model sigma.
# Measured 2026-09-17 at GP67H values: blurred null 5.5% confined / 4.6% directed on
# true free tracks; the old static, unblurred null gave 4.4% / 7.3%.
test_that("blurred CI null keeps false-positive rates near 5% on blurred free tracks", {
  skip_if_not_slow()
  set.seed(3)
  t_c <- 0.192; D <- 215; N <- 200; R <- 1 / 6.2
  gap <- c(0.72, 0.2, 0.05, 0.03)
  eco <- pmax(20, rlnorm(1e5, log(36), 0.6)); sig <- 100.8 / (4 * sqrt(eco))
  # raw CI samples: same function with the quantile step replaced
  src <- paste(deparse(sim_free_ci_null), collapse = "\n")
  src <- sub("q <- stats::quantile(v, c(0.05, 0.5, 0.95), names = FALSE)", "return(v)", src, fixed = TRUE)
  src <- sub("do.call(rbind, out)", "out", src, fixed = TRUE)
  ci_samples <- eval(parse(text = src))
  truth <- ci_samples(N, D_grid = D, n_sim = 2000, t_c = t_c, sigma = sig, gap_p = gap, R_blur = R)[[1]]
  nul   <- sim_free_ci_null(N, D_grid = D, n_sim = 1000, t_c = t_c, sigma = sig, gap_p = gap, R_blur = R)
  expect_lt(abs(mean(truth < nul$lo05) - 0.05), 0.025)
  expect_lt(abs(mean(truth > nul$hi95) - 0.05), 0.025)
})

test_that("sim_free_ci_null accepts a scalar or a vector sigma and R_blur = 0", {
  set.seed(1)
  a <- sim_free_ci_null(60, D_grid = 200, n_sim = 30, t_c = 0.2, sigma = 5)
  b <- sim_free_ci_null(60, D_grid = 200, n_sim = 30, t_c = 0.2, sigma = runif(100, 3, 6), R_blur = 1 / 6.2)
  expect_equal(nrow(a), 1L); expect_equal(nrow(b), 1L)
  expect_true(all(is.finite(c(a$med, b$med))))
})
