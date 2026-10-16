# Fiducial / contaminant scan (fiducial_scan_file and helpers, R/smlm_filters.R).

test_that("robust log z and recurrence", {
  expect_identical(robust_z_log10(rep(5, 4)), rep(0, 4))        # collapsed MAD flags nothing
  z <- robust_z_log10(c(10, 11, 9, 10, 12, 1000))
  expect_gt(z[6], 10)
  expect_identical(trace_recurrence(c(0, 50, 500), c(0, 0, 0), 100), c(1L, 1L, 0L))
  expect_identical(trace_recurrence(0, 0, 100), 0L)
})

test_that("a long-lived bright bead is flagged, ordinary traces are not", {
  set.seed(9)
  d <- sim_smlm_traces(grid_molecules(8, 150), blink_mean = 2, n_loc = 6, sigma = 5)
  d$tim <- d$tid * 2 + stats::runif(nrow(d), 0, 1)
  d$efo <- stats::rlnorm(nrow(d), log(5e4), 0.2)
  bead <- data.frame(tid = max(d$tid) + 1L, x = 2000 + stats::rnorm(60), y = 2000 + stats::rnorm(60),
                     tim = seq(0, 600, length.out = 60), efo = 5e6)
  res <- fiducial_scan_file(dplyr::bind_rows(d, bead), "f")
  expect_identical(res$tid[res$flagged], bead$tid[1])
  expect_identical(res$rule[res$flagged], "long+bright")
  expect_true(all(res$dataset == "f"))
  expect_null(fiducial_scan_file(d[, c("tid", "x")], "f"))
})
