# Precision cloud core (precision_cloud, R/smlm_filters.R).

test_that("precision cloud recovers the simulated per-axis scatter", {
  set.seed(8)
  d <- sim_smlm_traces(grid_molecules(10, 70), blink_mean = 1, n_loc = 20, sigma = c(6, 4))
  d$tim <- sort(stats::runif(nrow(d), 0, 600))
  res <- precision_cloud(list(f = d), min_loc = 5L)
  expect_false(is.null(res))
  expect_equal(res$sigma_x, 6, tolerance = 0.08)
  expect_equal(res$sigma_y, 4, tolerance = 0.08)
  expect_gt(res$anisotropy, 1.3)
  expect_false(res$has_z)
  expect_equal(length(res$dx), res$n_locs)            # under 40k: no subsample
  expect_equal(abs(mean(res$dx)), 0, tolerance = 1e-8)  # recentred at the trace means
})

test_that("precision cloud returns NULL when no trace is long enough", {
  set.seed(8)
  d <- sim_smlm_traces(grid_molecules(3, 70), blink_mean = 1, n_loc = 3, sigma = 5)
  d$tim <- seq_len(nrow(d))
  expect_null(precision_cloud(list(f = d), min_loc = 50L))
})
