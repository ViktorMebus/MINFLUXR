# The settings-based cache key. The server caches the filter -> thin -> cluster
# chain per file (prepared_df in server_smlm.R) under a key that hashes the same
# settings object the chain runs with (smlm_settings_key), plus the derived state
# the settings do not cover: the data version, the density-thinning
# keep-probabilities and the fiducial exclusion list. These tests check that the
# settings cover every field the chain reads, that the chain's output is a
# function of (data, settings), and that the key is stable across harmless type
# differences and changes when any setting changes. The server side is exercised
# end to end in test-smlm-server-smoke.R.

.server_src <- function() paste(readLines(file.path(.repo_root, "server_smlm.R"), warn = FALSE,
                                          encoding = "UTF-8"), collapse = "\n")

# A source that records every field the chain reads through the adapter.
.reads <- new.env()
registerS3method("[[", "smlm_recorder", function(x, i, ...) {
  assign(i, TRUE, envir = .reads); unclass(x)[[i]]
})
.recorder <- function(...) structure(list(...), class = "smlm_recorder")

.chain_df <- function(seed = 1) {
  set.seed(seed)
  df <- sim_smlm_traces(grid_molecules(8, 60), blink_mean = 1, sigma = 5,
                        dcr = rep(c(0.3, 0.62), 32))
  n <- nrow(df)
  df$efo <- stats::runif(n, 1e4, 1e5); df$cfr <- stats::runif(n); df$eco <- stats::rpois(n, 200)
  df$tim <- seq_len(n) / 100; df$dwell <- stats::runif(n, 1, 3); df$fluo <- 1L
  df
}

# Runs that between them take every branch: filters with trimming and the spread
# filter, two colours, centroid and vote de-blinking, every median-clustering
# option, and a single-channel scope.
.branch_runs <- list(
  list(efo_range = c(0, 2e5), cfr_range = c(0, 1), dcr_range = c(0, 1), eco_range = c(0, 1e3),
       tim_range = c(0, 1e4), dwell_range = c(0, 10), fluo_range = c(0, 5),
       tid_length_range = c(5, Inf), trim_outliers = TRUE, trim_k = 4, sd_min_locs = 5,
       view_SD = TRUE, target_sd_x = 20, target_sd_y = 20,
       enable_two_color = TRUE, use_pooled_dcr = FALSE, color_by_fluo = FALSE,
       enable_clustering = TRUE, deblink_method = "centroid", eps_centroid = 10,
       centroid_k = 2.8, deblink_anisotropic = TRUE, deblink_shared_spread = FALSE,
       enable_cluster2 = TRUE, cluster2_method = "dbscan", eps2 = 120, minpts2 = 4,
       cluster2_prefilter = TRUE, cluster2_prefilter_r = 150, cluster2_prefilter_min = 2,
       cluster2_prune = TRUE, cluster2_prune_tol = 2, cluster2_core_only = TRUE,
       cluster2_min_members = 2, cluster_group_scope = "within"),
  list(enable_two_color = TRUE, enable_clustering = TRUE, deblink_method = "vote",
       vote_frac = 0.3, vote_rho = 3, vote_link = "either", enable_cluster2 = TRUE,
       cluster2_method = "optics", cluster2_xi = 0.05, minpts2 = 4,
       cluster_group_scope = "only:1"),
  list(enable_two_color = TRUE, color_by_fluo = TRUE, enable_clustering = TRUE,
       deblink_shared_spread = TRUE, enable_cluster2 = TRUE, cluster2_method = "hdbscan",
       minpts2 = 4),
  list(enable_two_color = TRUE, dcr_group1 = c(0, 0.45), dcr_group2 = c(0.45, 1))
)

# CLAIM: every field the cached chain reads is a defined setting, so a key built
# from the settings covers everything the chain depends on through its settings.
test_that("the chain reads only defined settings", {
  rm(list = ls(.reads), envir = .reads)
  df <- .chain_df()
  for (run in .branch_runs) {
    s <- do.call(.recorder, run)
    # OPTICS warns when Xi finds no cluster in one channel; that is a result, not a fault
    out <- suppressWarnings(apply_clustering_pipeline(apply_meta_filters_batch(df, s), s))
    expect_gt(nrow(out), 0)
  }
  read <- ls(.reads)
  expect_gt(length(read), 35)                            # the branches really ran
  expect_equal(setdiff(read, names(SMLM_SETTINGS_SPEC)), character(0))
})

# GUARD: the server builds the key from the settings object and runs the chain
# with that same object; the hand-kept input list and its runtime scan are gone.
test_that("prepared_df keys and runs on one settings object", {
  src  <- .server_src()
  body <- regmatches(src, regexpr("prepared_df <- function\\((.|\\n)*?\\n  \\}\\n", src, perl = TRUE))
  expect_match(body, "s <- smlm_settings_from_input(input)", fixed = TRUE)
  expect_match(body, ".prep_key(fname, needs_cluster, two_color, s)", fixed = TRUE)
  expect_match(body, "apply_meta_filters_batch(df, s)", fixed = TRUE)
  expect_match(body, "apply_density_thinning(df, fname, s)", fixed = TRUE)
  expect_match(body, "apply_clustering_pipeline(df, s)", fixed = TRUE)
  expect_false(grepl(".PREP_KEY_INPUTS", src, fixed = TRUE))
  expect_false(grepl(".prep_key_selftest", src, fixed = TRUE))
})

# CLAIM: the chain's output is a function of the data and the settings alone,
# whichever form the settings arrive in.
test_that("output is identical from input, list, view and settings object", {
  df  <- .chain_df(2)
  run <- .branch_runs[[1]]
  rv  <- do.call(shiny::reactiveValues, run)
  a <- shiny::isolate(apply_clustering_pipeline(apply_meta_filters_batch(df, rv), rv))
  b <- apply_clustering_pipeline(apply_meta_filters_batch(df, run), run)
  s <- smlm_settings_from_input(run)
  c <- apply_clustering_pipeline(apply_meta_filters_batch(df, s), s)
  expect_identical(a, b)
  expect_identical(a, c)
})

# CLAIM: a key built on the settings object is stable: values that mean the same
# setting give the same key, so the cache does not miss on a widget returning 15
# where a settings file holds 15L or "15", or on a field the session never set.
test_that("a settings key ignores harmless type differences", {
  key <- function(x) smlm_settings_key(x)
  base <- key(list())
  expect_identical(key(list(minpts2 = 15)), base)
  expect_identical(key(list(minpts2 = 15L)), base)
  expect_identical(key(list(minpts2 = "15")), base)
  expect_identical(key(list(eps2 = 200L, enable_clustering = "FALSE")), base)
  expect_identical(key(list(dcr_group1 = c(0L, 0.369))), base)
  expect_identical(key(list(not_a_setting = 1, color_group1 = "#000000")), base)  # display only
})

# CLAIM: and sensitive: changing any one setting changes the key. Fields an
# inactive branch does not read change it too; that costs a cache miss, never a
# stale result.
test_that("changing any setting changes a settings key", {
  key  <- function(x) smlm_settings_key(x)
  base <- key(list())
  alt  <- list(flag = function(d) !d, num = function(d) if (is.na(d)) 1 else d + 1,
               int = function(d) d + 1L, range = function(d) if (is.null(d)) c(0, 1) else d + 1,
               choice = function(d) "")
  for (id in names(SMLM_SETTINGS_SPEC)) {
    sp <- SMLM_SETTINGS_SPEC[[id]]
    v  <- if (sp$type == "choice") {
      if (id == "cluster_group_scope") "only:1" else setdiff(sp$choices, sp$default)[1]
    } else alt[[sp$type]](sp$default)
    expect_false(identical(key(stats::setNames(list(v), id)), base), info = id)
  }
})
