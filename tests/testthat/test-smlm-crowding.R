# Crowding vs detection-time shift (crowding_for_file, crowding_null_test; R/smlm_nn.R).

.crowd_df <- function() {
  set.seed(6)
  # 6 structures of 4..12 traces at 50 nm spread, plus 60 background traces
  n <- c(4, 6, 8, 10, 12, 3)
  ctr <- cbind(seq(0, 2500, 500), 0)
  cl <- rep(seq_along(n), n)
  x <- c(ctr[cl, 1] + stats::rnorm(sum(n), 0, 50), stats::runif(60, -500, 3000))
  y <- c(stats::rnorm(sum(n), 0, 50), stats::runif(60, 500, 2000))
  tid <- seq_along(x)
  data.frame(tid = tid, x = x, y = y, tim = c(stats::rnorm(sum(n), 100, 3), stats::rnorm(60, 40, 3)),
             cluster2 = c(cl, rep(0L, 60)))
}

test_that("crowding rows: one per qualifying structure, shift against the background", {
  out <- crowding_for_file(.crowd_df(), "f", id_col = "tid", min_n = 4L)
  expect_equal(nrow(out$rows), 5L)                    # the 3-trace structure is too small
  expect_equal(out$diag[["too_small"]], 1L)
  expect_equal(out$rows$time_shift, rep(60, 5), tolerance = 0.1)
  expect_identical(out$rows$n_tids, c(4L, 6L, 8L, 10L, 12L))
  expect_true(all(out$rows$unit == "tid"))
  expect_null(out$null)
})

test_that("crowding names why a file gives nothing", {
  d <- .crowd_df()
  expect_equal(crowding_for_file(d[, names(d) != "cluster2"], "f")$diag[["no_c2"]], 1L)
  expect_equal(crowding_for_file(d, "f", use_3d = TRUE)$diag[["no_z"]], 1L)
  expect_equal(crowding_for_file(d[d$cluster2 > 0, ], "f")$diag[["no_background"]], 1L)
})

test_that("crowding null: permutation matrix, band and trend test", {
  set.seed(2)
  out <- crowding_for_file(.crowd_df(), "f", id_col = "tid", min_n = 4L, B = 99L)
  expect_identical(dim(out$null), c(5L, 99L))
  tst <- crowding_null_test(out$rows, out$null)
  expect_true(all(tst$res$null_lo <= tst$res$null_hi))
  expect_true(all(tst$res$time_shift > tst$res$null_hi))   # all structures detected late
  expect_equal(tst$stats$outside, 1)
  expect_true(tst$stats$p >= 1 / 100 && tst$stats$p <= 1)
  bare <- crowding_null_test(out$rows, NULL)
  expect_null(bare$stats)
  expect_true(all(is.na(bare$res$null_med)))
})

test_that("size vs NN: one row per structure with hull size", {
  out <- size_nn_for_file(.crowd_df(), "f", id_col = "tid", min_n = 4L)
  expect_identical(out$rows$size, c(4L, 6L, 8L, 10L, 12L))
  expect_equal(out$diag[["too_small"]], 1L)
  expect_true(all(out$rows$hull_size > 0 & out$rows$nn_median > 0))
  expect_true(all(out$rows$dim == 2L))
  expect_equal(size_nn_for_file(.crowd_df()[, 1:4], "f")$diag[["no_c2"]], 1L)
})
