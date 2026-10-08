# The data-chain helpers take `settings` (R/smlm_settings.R). During the migration
# the server still passes Shiny's `input`; as_smlm_settings() adapts it.

test_that("as_smlm_settings passes settings through and wraps everything else", {
  s <- smlm_settings(eps2 = 150)
  expect_identical(as_smlm_settings(s), s)
  expect_identical(as_smlm_settings(NULL)$eps2, 200)
  v <- as_smlm_settings(list(eps2 = "150", minpts2 = NULL, enable_cluster2 = 1))
  expect_s3_class(v, "smlm_settings_view")
  expect_identical(v$eps2, 150)                         # coerced like smlm_settings_from_input
  expect_identical(v[["minpts2"]], 15L)                 # missing -> default
  expect_true(v$enable_cluster2)
  expect_error(v$eps, "not an SMLM setting")            # unknown or partial name
  expect_identical(as_smlm_settings(v), v)
  # values from the app are used as entered; only code-built settings clamp
  expect_identical(as_smlm_settings(list(trim_k = 1))$trim_k, 1)
})

test_that("smlm_set overrides one setting and keeps the rest", {
  v  <- as_smlm_settings(list(eps2 = 120))
  v2 <- smlm_set(v, enable_two_color = TRUE)
  expect_true(v2$enable_two_color)
  expect_false(v$enable_two_color)                      # the original is untouched
  expect_identical(v2$eps2, 120)
  s2 <- smlm_set(smlm_settings(), minpts2 = 9)
  expect_s3_class(s2, "smlm_settings")
  expect_identical(s2$minpts2, 9L)
  expect_error(smlm_set(v, eps_2 = 1), "unknown setting")
})

# CLAIM: the view reads a field only when the helper asks for it, so a reactive
# that calls the pipeline depends on the same inputs as when the helper read
# input$... directly. Settings an inactive branch does not use must not trigger a
# re-run.
test_that("a reactive calling the pipeline depends only on the settings it uses", {
  set.seed(1)
  df <- sim_smlm_traces(grid_molecules(5, 80), blink_mean = 1, sigma = 5)
  rv <- shiny::reactiveValues(enable_clustering = TRUE, enable_cluster2 = FALSE,
                              eps2 = 200, vote_frac = 0.3, deblink_method = "centroid")
  n  <- 0
  r  <- shiny::reactive({ n <<- n + 1; apply_clustering_pipeline(df, rv) })
  step <- function(expr) { force(expr); shiny:::flushReact(); invisible(shiny::isolate(r())) }
  step(NULL);                          expect_equal(n, 1)
  step(rv$eps2 <- 150);                expect_equal(n, 1)   # median clustering is off
  step(rv$vote_frac <- 0.5);           expect_equal(n, 1)   # centroid method selected
  step(rv$enable_cluster2 <- TRUE);    expect_equal(n, 2)
  step(rv$eps2 <- 120);                expect_equal(n, 3)   # now it is used
})

test_that("small helpers take plain arguments", {
  df <- data.frame(channel = c(1L, 2L, 2L, 3L))
  expect_equal(channels_in_scope(df), 1:3)
  expect_equal(channels_in_scope(df, "within"), 1:3)
  expect_equal(channels_in_scope(df, "only:2"), 2L)
  expect_equal(channels_in_scope(df, "only:9"), 1:3)    # unknown channel: all of them
  expect_identical(deblink_minpts(c(7, Inf)), 7L)
  expect_identical(deblink_minpts(NULL), 4L)            # filter not rendered yet
  expect_identical(deblink_minpts(c(1, 10)), 2L)        # DBSCAN needs >= 2
  X <- matrix(stats::rnorm(30), ncol = 3)
  expect_named(build_center_args(X, "fixed", 0.1), "center_fixed")
  expect_equal(build_center_args(X, "reg", 0.3)$lambda, 0.3)
  expect_named(build_center_args(X, "free", NA), "center_init")
})

test_that("channel helpers agree between settings and a plain input list", {
  set.seed(3)
  df <- sim_smlm_traces(grid_molecules(6, 80), blink_mean = 1, sigma = 5,
                        dcr = rep(c(0.3, 0.62), 18))
  inp <- list(enable_two_color = TRUE, use_pooled_dcr = FALSE)
  s   <- smlm_settings(enable_two_color = TRUE)
  expect_identical(resolve_channels(df, inp), resolve_channels(df, s))
  expect_identical(assign_two_color_group(df, list()), assign_two_color_group(df, smlm_settings()))
  expect_identical(assign_analysis_group(df, inp), assign_analysis_group(df, s))
  # grouping_mode has no control in the app; "dcr" applies the DCR group ranges
  g <- assign_analysis_group(df, s, grouping_mode = "dcr")
  expect_setequal(unique(g$analysis_group), c(1L, 2L))
})
