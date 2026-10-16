# Local density & hotspots under complete spatial randomness (CSR) in a footprint.
# Full map: tests/validation/hotspot_csr_null.R (3 shapes x 3 densities x 20 reps).
suppressMessages(suppressWarnings(source(file.path(.repo_root, "packages.R"))))

.rod_csr <- function(n, seed) {
  set.seed(seed)
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    x <- runif(4 * n, 0, 10000); y <- runif(4 * n, 0, 6000)
    k <- ((x - 5000) / 4000)^2 + ((y - 3000) / 1000)^2 <= 1
    out <- rbind(out, data.frame(x = x[k], y = y[k]))
  }
  out <- out[seq_len(n), ]; out$z <- 0; out
}

# CLAIM (fixed 2026-09-17): whole-quadrat VMR reads ~1 on CSR in a narrow footprint
# at any density; the edge-quadrat version read 2.56 at 5000 units in a 2 x 8 µm rod.
test_that("dispersion index is ~1 on CSR in a narrow rod, sparse and dense", {
  skip_if_not_slow()
  for (n in c(200L, 5000L)) {
    pts <- .rod_csr(n, 10 + n)
    sup <- estimate_support(pts, mode = "grid", h_nm = 100, close_nm = 200, return_mask = TRUE)
    di  <- compute_dispersion_index(pts, quadrat_nm = 250, support = sup)
    expect_false(di$edge_quadrats)
    expect_lt(abs(di$vmr - 1), 0.25)
  }
})

# CLAIM (fixed 2026-09-17): eps calibration converges in dense cells (the grid is
# extended below 40 nm) and keeps random units out of foci; the largest focus
# relative to the cell's own random 95th percentile stays near or below 1.
test_that("eps calibration converges on a dense CSR rod and the focus ratio reads chance", {
  skip_if_not_slow()
  pts <- .rod_csr(5000L, 7)
  sup <- estimate_support(pts, mode = "grid", h_nm = 100, close_nm = 200, return_mask = TRUE)
  ec  <- calibrate_dbscan_eps(pts, support = sup, minPts = 5L, target_frac = 0.05, nsim = 15L)
  expect_false(ec$hit_floor)
  expect_true(ec$grid_extended)
  row <- hotspot_summary_row(pts, support = sup, eps_cal = ec)
  expect_lt(row$frac_in_foci, 0.10)
  expect_lt(row$max_focus_ratio, 1.5)
  expect_lt(abs(row$cv_ratio - 1), 0.25)
})

# The per-file core the Run handler calls: one bundle per file, window mode recorded.
test_that("hotspot_for_points bundles enhancement, eps, foci and the summary row", {
  pts <- .rod_csr(400L, 77)
  sup <- estimate_support(pts, mode = "grid", h_nm = 100, close_nm = 200, return_mask = TRUE)
  set.seed(3)
  o <- hotspot_for_points(pts, support = sup, k = 5L, nsim = 5L, tag = "f")
  expect_named(o, c("tag", "pts", "support", "enh", "le", "eps_cal", "foci", "row"))
  expect_identical(o$row$dataset, "f")
  expect_identical(o$row$window_mode, sup$mode %||% NA_character_)
  expect_equal(length(o$enh), sum(o$le$keep))
  expect_null(hotspot_for_points(pts[1:2, ], support = sup))
})
