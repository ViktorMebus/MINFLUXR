# End-to-end: the SMLM server's own cached chain (prepared_df) and its
# settings-based cache key, driven through Shiny's `input` as the app does. Loads
# the whole app, so it runs with the slow tests.

.load_app <- function() {
  # global.R sources the app's files into the global environment, as the app does
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  suppressMessages(source("global.R"))
}
.smoke_df <- function() {
  set.seed(1)
  df <- sim_smlm_traces(grid_molecules(8, 60), blink_mean = 1, sigma = 5,
                        dcr = rep(c(0.3, 0.62), 32))
  df$mol <- NULL
  df
}
.smoke_inp <- list(enable_clustering = TRUE, enable_two_color = TRUE, enable_cluster2 = TRUE,
                   eps2 = 120, minpts2 = 4, tid_length_range = c(5, Inf))

test_that("the server's cached chain matches the settings-driven chain", {
  skip_if_not_slow()
  .load_app()
  df <- .smoke_df()
  shiny::testServer(smlm_server, {
    smlm_file_store(list(f1 = df))
    do.call(session$setInputs, .smoke_inp)
    a <- prepared_df("f1")
    expect_equal(length(unique(a$cluster)), 64L)
    s <- smlm_settings_from_input(.smoke_inp)
    expect_identical(a, apply_clustering_pipeline(apply_meta_filters_batch(df, s), s))
    c2 <- cluster2_defaults()                            # Reset = a fresh session
    expect_equal(c(c2$eps2, c2$minpts2), c(SMLM_DEFAULTS$eps2, SMLM_DEFAULTS$minpts2))
  })
})

# CLAIM: the key is stable where it should be (a repeat call is a cache hit; a
# display-only input or the same value in another type leaves it alone) and moves
# where it must (any setting, and each kind of derived state).
test_that("the cache key hits, misses and invalidates as it should", {
  skip_if_not_slow()
  .load_app()
  df <- .smoke_df()
  shiny::testServer(smlm_server, {
    smlm_file_store(list(f1 = df))
    do.call(session$setInputs, .smoke_inp)
    key <- function() .prep_key("f1", TRUE, TRUE)
    a  <- prepared_df("f1")
    k0 <- key()
    expect_length(prep_cache_rv(), 1L)
    expect_identical(prepared_df("f1"), a)               # served from the cache
    expect_length(prep_cache_rv(), 1L)

    session$setInputs(color_group1 = "#000000", color_mode = "z", numeric_palette = "Viridis")
    expect_identical(key(), k0)                          # display only
    session$setInputs(eps2 = 120L, minpts2 = 4L)
    expect_identical(key(), k0)                          # same values, other type

    changes <- list(eps2 = 150, eps_centroid = 12, enable_cluster2 = FALSE,
                    trim_outliers = TRUE, efo_range = c(0, 1), deblink_method = "vote",
                    cluster_group_scope = "only:1", dens_norm_unit = "tid")
    for (id in names(changes)) {
      old <- .smoke_inp[[id]] %||% SMLM_DEFAULTS[[id]]
      do.call(session$setInputs, changes[id])
      expect_false(identical(key(), k0), info = id)
      do.call(session$setInputs, stats::setNames(list(old), id))
      expect_identical(key(), k0, info = paste(id, "restored"))
    }

    # A changed setting re-runs the chain. Within one channel the molecules sit
    # ~85 nm apart, so eps2 = 70 leaves no Median Cluster where 120 found two.
    session$setInputs(eps2 = 70)
    b <- prepared_df("f1")
    expect_length(prep_cache_rv(), 2L)
    expect_equal(sum(b$cluster2 > 0, na.rm = TRUE), 0)
    expect_gt(sum(a$cluster2 > 0, na.rm = TRUE), 0)
    session$setInputs(eps2 = 120)
    expect_identical(prepared_df("f1"), a)              # back to the first entry

    # Derived state 1: the data version (axial calibration rescales z in place)
    z0 <- z_scale_applied(); z_scale_applied(0.5)
    expect_false(identical(key(), k0)); z_scale_applied(z0)
    expect_identical(key(), k0)

    # Derived state 2: fiducial exclusion drops the flagged traces
    fiducial_results(data.frame(dataset = "f1", tid = c(1, 2, 3), flagged = TRUE))
    session$setInputs(fid_exclude = TRUE)
    expect_false(identical(key(), k0))
    f <- prepared_df("f1")
    expect_false(any(f$tid %in% 1:3))
    expect_true(any(a$tid %in% 1:3))
    session$setInputs(fid_exclude = FALSE); fiducial_results(NULL)
    expect_identical(key(), k0)

    # Derived state 3: density thinning to half the density
    # (The density panel's readouts warn on this minimal profile, which lacks the
    # columns their plots draw; only the keep-probability matters here.)
    density_profile_rv(data.frame(dataset = "f1", density = 10, unit = "cluster"))
    suppressWarnings(session$setInputs(dens_norm_enable = TRUE, dens_thin_enable = TRUE,
                                       dens_target_mode = "manual", dens_target_manual = 5))
    expect_false(identical(key(), k0))
    th <- prepared_df("f1")
    expect_equal(length(unique(th$cluster)), 32L)       # half of the 64 molecules
  })
})
