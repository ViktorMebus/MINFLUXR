# SMLM settings (R/smlm_settings.R): one set of defaults for the UI widgets, the
# helper functions' fallbacks and the reset buttons.

.read_src <- function(f) paste(readLines(file.path(.repo_root, f), warn = FALSE, encoding = "UTF-8"),
                               collapse = "\n")
# Settings fields a function reads (`s$<id>` or `input$<id>`), from its deparsed body
.fields_read <- function(fn) {
  txt <- paste(deparse(get(fn)), collapse = "\n")
  hits <- regmatches(txt, gregexpr("(input|\\bs)[$][A-Za-z0-9_.]+", txt, perl = TRUE))[[1]]
  unique(sub("^(input|s)[$]", "", hits))
}

test_that("defaults are the values a fresh session shows", {
  s <- smlm_settings()
  expect_s3_class(s, "smlm_settings")
  expect_setequal(names(s), names(SMLM_SETTINGS_SPEC))
  # the four that used to disagree between the helpers and the UI
  expect_equal(s$minpts2, 15L)
  expect_equal(s$vote_frac, 0.3)
  expect_equal(s$vote_rho, 3)
  expect_equal(s$eps2, 200)
  expect_null(s$efo_range)                              # no filter until data sets one
  expect_equal(s$sd_min_locs, MIN_LOCS_FOR_SPREAD_DEFAULT)
})

test_that("overrides are coerced, checked and clamped", {
  s <- smlm_settings(eps2 = "150", enable_clustering = 1, cluster_group_scope = "only:2")
  expect_identical(s$eps2, 150)
  expect_true(s$enable_clustering)
  expect_identical(s$cluster_group_scope, "only:2")
  expect_identical(smlm_settings(minpts2 = 7.4)$minpts2, 7L)
  expect_warning(m <- smlm_settings(minpts2 = 2)$minpts2, "clamped")
  expect_identical(m, 4L)                               # the app's floor
  expect_warning(k <- smlm_settings(cluster2_method = "kmeans")$cluster2_method, "not one of")
  expect_identical(k, "dbscan")
  expect_warning(smlm_settings(cluster_group_scope = "all"), "within")
  expect_error(smlm_settings(eps_2 = 150), "unknown setting")   # a typo is an error
  expect_error(smlm_settings(150), "named")
  expect_true(is.na(smlm_settings(target_sd_x = NA)$target_sd_x))
})

# CLAIM: `$` matches names exactly and refuses names that are not settings, so a
# typo or a prefix is an error rather than a longer field or a silent NULL.
test_that("$ on settings matches exactly and refuses unknown names", {
  s <- smlm_settings()
  expect_error(s$cluster2_pre, "not an SMLM setting")
  expect_error(s$eps, "not an SMLM setting")
  expect_identical(s$cluster2_prefilter, FALSE)
})

test_that("settings read from the app fall back to the defaults silently", {
  inp <- list(eps2 = 120, enable_cluster2 = TRUE, minpts2 = NULL, vote_frac = "x")
  expect_silent(s <- smlm_settings_from_input(inp))
  expect_equal(s$eps2, 120)
  expect_true(s$enable_cluster2)
  expect_equal(s$minpts2, 15L)
  expect_equal(s$vote_frac, 0.3)
  rv <- shiny::reactiveValues(eps_centroid = 12, deblink_method = "vote")
  s2 <- shiny::isolate(smlm_settings_from_input(rv))
  expect_equal(s2$eps_centroid, 12)
  expect_identical(s2$deblink_method, "vote")
})

# CLAIM: the settings object drops into any helper that takes `input`, with the
# same result as the plain list it replaces.
test_that("helpers give identical results from a settings object and a plain list", {
  set.seed(1)
  mol <- grid_molecules(8, 60)
  df  <- sim_smlm_traces(mol, blink_mean = 1, sigma = 5, dcr = rep(c(0.3, 0.62), 32))
  s   <- smlm_settings(enable_clustering = TRUE, enable_two_color = TRUE,
                       enable_cluster2 = TRUE, eps2 = 120, minpts2 = 4,
                       tid_length_range = c(5, Inf), trim_outliers = TRUE)
  a <- apply_clustering_pipeline(apply_meta_filters_batch(df, s), s)
  b <- apply_clustering_pipeline(apply_meta_filters_batch(df, unclass(s)), unclass(s))
  expect_identical(a, b)
})

# GUARD: every setting the cached data chain reads is defined here. A field read
# but missing from the spec would fall back to a hidden default and stay out of
# any settings-based cache key.
test_that("the data chain reads only settings that are defined", {
  chain <- c("apply_meta_filters_batch", "resolve_channels",
             "deblink_threshold_for_channel", "channels_in_scope", "deblink_minpts",
             "apply_clustering_pipeline")
  used <- unique(unlist(lapply(chain, .fields_read)))
  expect_equal(setdiff(used, names(SMLM_SETTINGS_SPEC)), character(0))
})

# GUARD: no literal fallback for a defined setting in the helpers or the server;
# they must read SMLM_DEFAULTS. Two exceptions are deliberate: the length-filter
# "no filter" range, and the minPts floor observer, whose 4 is the floor itself.
test_that("no literal defaults for defined settings in helpers or server", {
  ids <- paste(names(SMLM_SETTINGS_SPEC), collapse = "|")
  pat <- sprintf("input\\$(%s)\\s*%%\\|\\|%%\\s*(?!SMLM_DEFAULTS)[^,)\\s]+", ids)
  helper_parts <- setdiff(file.path("R", list.files(file.path(.repo_root, "R"), "^smlm_.*\\.R$")),
                          "R/smlm_settings.R")
  expect_gt(length(helper_parts), 10)
  hits <- unlist(lapply(c(helper_parts, "server_smlm.R"), function(f)
    regmatches(.read_src(f), gregexpr(pat, .read_src(f), perl = TRUE))[[1]]))
  allowed <- c("input$tid_length_range %||% c(0", "input$minpts2 %||% 4L")
  expect_equal(setdiff(hits, allowed), character(0))
})

# GUARD: the UI's static widgets take their initial value from SMLM_DEFAULTS.
test_that("UI widgets take their defaults from SMLM_DEFAULTS", {
  ui <- .read_src("ui_smlm.R")
  for (id in names(SMLM_SETTINGS_SPEC)) {
    m <- regexpr(sprintf('\\(\\s*"%s"[^\\n]*?(\\n[^\\n]*?){0,3}?\\b(value|selected)\\s*=\\s*[^,)]+', id),
                 ui, perl = TRUE)
    if (m < 0) next                                     # widget built in the server
    expect_match(regmatches(ui, m), "SMLM_DEFAULTS\\$", info = id)
  }
})
