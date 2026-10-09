# NN-line overlay (build_nn_segments, .nn_pairs_scoped in R/smlm_nn.R): options come
# from smlm_nn_options(), never from Shiny's input directly.

.nn_df <- function() {
  set.seed(5)
  df <- sim_smlm_traces(grid_molecules(6, 70), blink_mean = 1, sigma = 5,
                        dcr = rep(c(0.3, 0.62), 18))
  df <- apply_clustering_pipeline(df, smlm_settings(enable_clustering = TRUE, enable_two_color = TRUE))
  df <- attach_channel_colors(df)
  df$color <- df$channel_color
  df
}

test_that("NN options are read once with defaults", {
  o <- smlm_nn_options(list())
  expect_equal(o$nn_k, 4); expect_identical(o$nn_mplex_scope, "within")
  expect_false(o$deblinked); expect_equal(o$distance_range, c(0, Inf))
  expect_false(o$within_cluster2); expect_equal(o$label_offset, 5)
  o2 <- smlm_nn_options(list(nn_within_cluster2 = TRUE))   # needs median clustering on too
  expect_false(o2$within_cluster2)
  rv <- shiny::reactiveValues(nn_k = 2, enable_clustering = TRUE)
  expect_equal(shiny::isolate(smlm_nn_options(rv))[c("nn_k", "deblinked")], list(nn_k = 2, deblinked = TRUE))
})

test_that("NN lines honour the channel scope and the distance range", {
  df <- .nn_df()
  base <- list(nn_k = 2, deblinked = TRUE, distance_range = c(0, Inf), label_offset = 5)
  within <- build_nn_segments(df, modifyList(base, list(nn_mplex_scope = "within")))$edges
  across <- build_nn_segments(df, modifyList(base, list(nn_mplex_scope = "across")))$edges
  expect_true(all(within$ch1 == within$ch2))
  expect_true(all(across$ch1 != across$ch2))
  # molecules sit 70 nm apart (99 nm on the diagonal, reached with 6 neighbours):
  # an 85 nm limit keeps the first and drops the second
  full  <- build_nn_segments(df, modifyList(base, list(nn_k = 6, nn_mplex_scope = "all")))$edges
  short <- build_nn_segments(df, modifyList(base, list(nn_k = 6, nn_mplex_scope = "all", distance_range = c(0, 85))))$edges
  expect_true(any(full$dist > 85))
  expect_gt(nrow(short), 0)
  expect_true(all(short$dist <= 85))
  # molecules (de-blinked) give fewer endpoints than traces
  per_mol <- build_nn_segments(df, modifyList(base, list(nn_mplex_scope = "all")))
  per_tid <- build_nn_segments(df, modifyList(base, list(deblinked = FALSE, nn_mplex_scope = "all")))
  expect_lt(nrow(per_mol$edges), nrow(per_tid$edges))
})
