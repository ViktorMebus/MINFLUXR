# Settings files (Download / Load settings in the SMLM tab): the smlm_settings block,
# the older layout, and the restore through the widgets.

# Non-default values across every stage, including the awkward types: an NA spread
# target, a NULL (unset) range, a length-2 range, an integer, a single-channel scope.
.custom_settings <- function() smlm_settings(
  efo_range = c(1e4, 2e5), tid_length_range = NULL, trim_outliers = TRUE, trim_k = 5,
  sd_min_locs = 7, view_SD = TRUE, target_sd_x = NA, target_sd_y = 12,
  enable_two_color = TRUE, color_by_fluo = TRUE, dcr_group1 = c(0, 0.4), dcr_group2 = c(0.41, 1),
  enable_clustering = TRUE, deblink_method = "vote", eps_centroid = 12, vote_frac = 0.4,
  enable_cluster2 = TRUE, cluster2_method = "hdbscan", eps2 = 150, minpts2 = 9,
  cluster2_prune = TRUE, cluster2_prune_tol = 2.5, cluster2_min_members = 3,
  cluster_group_scope = "only:1", fid_exclude = TRUE, dens_norm_enable = TRUE,
  dens_norm_unit = "tid")

# The file exactly as the app writes and reads it (same jsonlite options).
.json_round_trip <- function(x) {
  f <- tempfile(fileext = ".json"); on.exit(unlink(f))
  writeLines(jsonlite::toJSON(x, pretty = TRUE, auto_unbox = TRUE, na = "null"), f, useBytes = TRUE)
  jsonlite::fromJSON(f)
}

# CLAIM: a settings file written with the smlm_settings block gives back exactly the
# settings that were saved, every one of them.
test_that("the smlm_settings block survives a JSON round trip", {
  s   <- .custom_settings()
  raw <- .json_round_trip(list(settings_version = SMLM_SETTINGS_FILE_VERSION,
                               smlm_settings = unclass(s), z_scale = 0.8))
  sf  <- smlm_settings_from_file(raw)
  expect_identical(sf$layout, "block")
  expect_setequal(sf$present, names(SMLM_SETTINGS_SPEC))
  expect_identical(smlm_settings_key(sf$settings), smlm_settings_key(s))
  expect_true(is.na(sf$settings$target_sd_x))
  expect_null(sf$settings$tid_length_range)
  expect_identical(sf$settings$minpts2, 9L)
})

# CLAIM: an older file (no block) is read from where it kept each setting, and the
# settings it does not hold are reported as absent, so a loader leaves them alone.
test_that("an older settings file is read from its own layout", {
  raw <- .json_round_trip(list(
    efo_range = c(1e4, 2e5), trim_outliers = TRUE, enable_clustering = TRUE,
    eps_centroid = 12, eps2 = 150, minpts2 = 9,
    two_color = list(enable = TRUE, use_fluo = FALSE, dcr_group1 = c(0, 0.4),
                     dcr_group2 = c(0.41, 1), color_group1 = "#1f77b4"),
    analysis = list(use_pooled_dcr = TRUE, dens_norm_enable = FALSE),
    fiducial = list(fid_exclude = TRUE, fid_efo_z = 3)))
  sf <- smlm_settings_from_file(raw)
  expect_identical(sf$layout, "legacy")
  expect_setequal(sf$present, c("efo_range", "trim_outliers", "enable_clustering",
                                "eps_centroid", "eps2", "minpts2", "enable_two_color",
                                "color_by_fluo", "dcr_group1", "dcr_group2",
                                "use_pooled_dcr", "dens_norm_enable", "fid_exclude"))
  expect_true(sf$settings$enable_two_color)
  expect_true(sf$settings$use_pooled_dcr)
  expect_true(sf$settings$fid_exclude)
  expect_equal(sf$settings$dcr_group2, c(0.41, 1))
  expect_false("cluster2_method" %in% sf$present)       # not in the file: not touched
})

test_that("settings written into the older layout read back the same", {
  s   <- .custom_settings()
  raw <- smlm_settings_to_legacy(list(two_color = list(color_group1 = "#000000")), s,
                                 names(SMLM_SETTINGS_SPEC))
  expect_identical(raw$two_color$color_group1, "#000000")   # neighbours are kept
  expect_true(raw$two_color$enable)
  expect_true(raw$fiducial$fid_exclude)
  expect_true(raw$analysis$dens_norm_enable)
  expect_true("tid_length_range" %in% names(raw))           # NULL kept as a field
  sf <- smlm_settings_from_file(.json_round_trip(raw))
  expect_identical(smlm_settings_key(sf$settings), smlm_settings_key(s))
})

# CLAIM: every setting gets a widget update, and each switch comes after the
# parameters it puts into effect.
test_that("the restore plan covers every setting with switches last", {
  plan <- smlm_settings_restore_plan(.custom_settings(), names(SMLM_SETTINGS_SPEC))
  expect_setequal(plan$id, names(SMLM_SETTINGS_SPEC))
  expect_identical(utils::tail(plan$id, 3),
                   c("enable_two_color", "enable_clustering", "enable_cluster2"))
  w <- stats::setNames(plan$widget, plan$id)
  expect_identical(unname(w[c("trim_outliers", "eps2", "minpts2", "dcr_group1",
                              "deblink_method", "cluster_group_scope", "dens_norm_unit")]),
                   c("checkbox", "numeric", "numeric", "range", "radio", "select", "select"))
  v <- stats::setNames(attr(plan, "values"), plan$id)
  expect_identical(v$cluster2_method, "hdbscan")
  expect_equal(v$dcr_group1, c(0, 0.4))
  expect_identical(smlm_settings_restore_plan(smlm_settings(), c("eps2", "nonsense"))$id, "eps2")
})

# End to end on the real server: session 1 saves with non-default settings, a fresh
# session 2 loads the file, and every setting reaches its widget. Widget updates
# are recorded by stand-ins for shiny's update functions (testServer does not move
# widgets); the main plot trigger is stubbed because the plot path fails under
# testServer for an unrelated reason (missing browser-provided inputs).
test_that("a saved settings file restores every setting in a new session", {
  skip_if_not_slow()
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  suppressMessages(source("global.R"))
  set.seed(1)
  df <- sim_smlm_traces(grid_molecules(6, 60), blink_mean = 1, sigma = 5,
                        dcr = rep(c(0.3, 0.62), 18))
  df$mol <- NULL; df$.dataset_tag <- "t1"
  n <- nrow(df)                                          # metadata the filter widgets read
  df$efo <- stats::runif(n, 1e4, 1e5); df$cfr <- stats::runif(n); df$eco <- stats::rpois(n, 200)
  df$tim <- seq_len(n) / 100; df$dwell <- 1; df$fluo <- 1L
  inp <- list(enable_clustering = TRUE, deblink_method = "vote", eps_centroid = 12,
              vote_frac = 0.4, enable_cluster2 = TRUE, cluster2_method = "hdbscan",
              cluster2_prune = TRUE, cluster2_prune_tol = 2.5, cluster2_min_members = 3,
              minpts2 = 9, eps2 = 150, dens_norm_unit = "tid", cluster_group_scope = "only:1",
              trim_outliers = TRUE, trim_k = 5, z_scale = 0.8, efo_range = c(0, 1e9))
  path <- tempfile(fileext = ".json"); on.exit(unlink(path), add = TRUE)

  shiny::testServer(smlm_server, {
    smlm_file_store(list(f1 = df)); uploaded_df_rv(df)
    do.call(session$setInputs, inp)
    writeLines(jsonlite::toJSON(get_ui_settings(), pretty = TRUE, auto_unbox = TRUE,
                                na = "null"), path, useBytes = TRUE)
  })
  saved <- jsonlite::fromJSON(path)
  expect_setequal(names(saved$smlm_settings), names(SMLM_SETTINGS_SPEC))

  rec <- new.env()
  spies <- c("updateCheckboxInput", "updateNumericInput", "updateRadioButtons", "updateSelectInput")
  for (f in spies) assign(f, local({ fn <- f; function(session, inputId, ...) {
    a <- list(...); assign(inputId, a$value %||% a$selected, envir = rec); invisible(NULL) } }),
    envir = globalenv())
  on.exit(rm(list = spies, envir = globalenv()), add = TRUE)

  # The display panels (histogram sliders, overlay) warn when restored in a session
  # that never rendered them; only the settings are under test here.
  suppressWarnings(shiny::testServer(smlm_server, {
    assign("plot_trigger", function(...) invisible(NULL), envir = session$env)
    smlm_file_store(list(f1 = df)); uploaded_df_rv(df)
    session$setInputs(load_settings = data.frame(name = "s.json", size = file.size(path),
                                                 type = "application/json", datapath = path))
    session$setInputs(.flush = 1)                        # run the onFlushed restore
  }))
  got <- as.list(rec)
  for (id in setdiff(names(inp), "efo_range"))
    expect_equal(got[[id]], inp[[id]], info = id)
  # every setting the plan restores through a checkbox/numeric/radio/select arrived
  plan_ids <- setdiff(names(SMLM_SETTINGS_SPEC), .SETTINGS_FILTER_IDS <- c(
    "efo_range", "cfr_range", "dcr_range", "eco_range", "tim_range", "dwell_range",
    "fluo_range", "tid_length_range", "sd_min_locs", "trim_outliers", "trim_k",
    "view_SD", "target_sd_x", "target_sd_y", "fid_exclude", "use_pooled_dcr",
    "dens_norm_enable", "dens_thin_enable"))
  plan_ids <- setdiff(plan_ids, c("dcr_group1", "dcr_group2"))   # range widgets (shinyWidgets)
  expect_equal(setdiff(plan_ids, names(got)), character(0))
})
