# Batch NN per-file core (batch_nn_for_file, R/smlm_nn.R) and channel facets.

.bnn_df <- function() {
  set.seed(4)
  sim_smlm_traces(grid_molecules(6, 70), blink_mean = 1, sigma = 5, dcr = rep(c(0.3, 0.62), 18))
}

test_that("channel facets split only when asked and when there are channels", {
  df <- data.frame(channel = c(1L, 1L, 2L), channel_label = c("A", "A", "B"), x = 1:3)
  expect_length(channel_facets(df), 1L)
  f <- channel_facets(df, by_channel = TRUE)
  expect_identical(vapply(f, `[[`, "", "label"), c("A", "B"))
  expect_equal(nrow(f[[2]]$df), 1L)
  expect_length(channel_facets(df[df$channel == 1L, ], by_channel = TRUE), 1L)
})

test_that("batch NN per file: units, groups and facets", {
  df <- .bnn_df()
  s  <- smlm_settings(enable_clustering = TRUE, enable_two_color = TRUE)
  tid <- batch_nn_for_file(df, "f", s, unit = "tid", k = 1)
  mol <- batch_nn_for_file(df, "f", s, unit = "pcid", k = 1)
  expect_named(tid, "f")
  expect_lt(nrow(mol$f), nrow(tid$f))                  # molecules are fewer than traces
  expect_true(all(tid$f$distance > 0))
  btw <- batch_nn_for_file(df, "f", s, unit = "pcid", k = 1, cross_groups = c(1L, 2L))
  expect_gt(nrow(btw$f), 0)
  fac <- batch_nn_for_file(df, "f", s, unit = "pcid", k = 1, by_channel = TRUE)
  expect_setequal(names(fac), c("f · A", "f · B"))
  expect_equal(sum(vapply(fac, nrow, 1L)), nrow(mol$f))  # facets partition the molecules
  # a cross-group run is never faceted
  expect_named(batch_nn_for_file(df, "f", s, unit = "pcid", cross_groups = c(1L, 2L),
                                 by_channel = TRUE), "f")
  expect_length(batch_nn_for_file(df[0, ], "f", s, unit = "tid"), 0L)
})
