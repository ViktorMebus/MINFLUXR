# Labelling-density normalisation (section 8b): footprint, per-file density and
# unit-level p-thinning. Cells are an L-shape (4 x 3 um with one quadrant empty).
.l_cell <- function(n) {
  x <- stats::runif(3 * n, 0, 4000); y <- stats::runif(3 * n, 0, 3000)
  k <- !(x > 2000 & y > 1500)
  cbind(x[k], y[k])[seq_len(n), ]
}
.l_area <- 4000 * 3000 - 2000 * 1500

test_that("the grid footprint of a densely sampled rectangle is its area", {
  set.seed(9)
  pts <- data.frame(x = stats::runif(5000, 0, 4000), y = stats::runif(5000, 0, 3000), z = 0)
  sup <- estimate_support(pts, mode = "grid", h_nm = 100, close_nm = 200)
  expect_equal(sup$support / 12e6, 1, tolerance = 0.03)
})

# CLAIM: two cells of the same shape with a 4x labelling difference read a 4x
# density ratio when both footprints are measured on one common gap-fill radius
# (the app's auto radius: the largest plateau start across files).
# Measured at 600 vs 2400 units: ratio 3.73-3.80 (the sparser footprint still runs
# ~2% small). Below ~300 units per cell the footprint shrinks further and the ratio
# compresses toward 1; see the note in the density panel.
test_that("density ratio tracks a 4x labelling difference", {
  set.seed(1)
  lo <- sim_smlm_traces(.l_cell(600),  blink_mean = 0, n_loc = 6, sigma = 5)
  hi <- sim_smlm_traces(.l_cell(2400), blink_mean = 0, n_loc = 6, sigma = 5)
  cl <- max(vapply(list(lo, hi), function(d)
    plateau_close_nm(support_closing_sweep(d[, c("x", "y", "z")], h_nm = 100)), numeric(1)))
  dl <- compute_unit_density(lo, close_nm = cl)
  dh <- compute_unit_density(hi, close_nm = cl)
  expect_equal(dl$support_nm / .l_area, 1, tolerance = 0.05)
  expect_equal(dh$density / dl$density, 4, tolerance = 0.08)
  expect_equal(dh$dims, 2L)
  # density is reported per um^2
  expect_equal(dh$density, dh$n_units / dh$support_nm * 1e6)
})

# CLAIM: thinning the denser cell to p = 0.25 and keeping its PRE-thinning footprint
# lands it on the sparser cell's density.
test_that("p-thinning with the inherited footprint equalises density", {
  set.seed(2)
  lo <- sim_smlm_traces(.l_cell(600),  blink_mean = 0, n_loc = 6, sigma = 5)
  hi <- sim_smlm_traces(.l_cell(2400), blink_mean = 0, n_loc = 6, sigma = 5)
  dl <- compute_unit_density(lo, close_nm = 300)
  dh <- compute_unit_density(hi, close_nm = 300)
  dt <- compute_unit_density(thin_units_to_p(hi, p = 0.25), close_nm = 300,
                             support_override = dh$support_nm)
  expect_equal(dt$support_mode, "inherited (pre-thinning)")
  expect_equal(dt$density / dl$density, 1, tolerance = 0.08)
})

# CLAIM: thinning keeps or drops whole units, hits round(p x n) exactly, is
# reproducible from its seed, and never touches the global RNG stream.
test_that("thin_units_to_p keeps whole units and leaves the RNG alone", {
  set.seed(3)
  df <- sim_smlm_traces(grid_molecules(10, 80), blink_mean = 1, sigma = 5)
  n  <- length(unique(df$tid))
  set.seed(42); before <- stats::runif(1)
  set.seed(42); t1 <- thin_units_to_p(df, p = 0.3, seed = 7); after <- stats::runif(1)
  expect_identical(before, after)
  expect_equal(length(unique(t1$tid)), round(0.3 * n))
  kept <- unique(t1$tid)
  expect_equal(nrow(t1), sum(df$tid %in% kept))      # every kept trace is complete
  expect_identical(t1, thin_units_to_p(df, p = 0.3, seed = 7))
  expect_identical(thin_units_to_p(df, p = 1), df)   # never up-samples
  expect_equal(nrow(thin_units_to_p(df, p = 0)), 0L)
})

test_that("unit = cluster thins de-blinked molecules, not traces", {
  set.seed(4)
  df  <- apply_clustering_pipeline(sim_smlm_traces(grid_molecules(8, 80), blink_mean = 2,
                                                   sigma = 5), smlm_input())
  t1  <- thin_units_to_p(df, unit = "cluster", p = 0.5, seed = 1)
  expect_equal(length(unique(t1$cluster)), round(0.5 * length(unique(df$cluster))))
  expect_equal(nrow(t1), sum(df$cluster %in% unique(t1$cluster)))
  u <- unit_centroids(df, unit = "cluster")
  expect_equal(nrow(u), length(unique(df$cluster)))
})

# CLAIM (methods doc 8b): independent p-thinning leaves g(r) unchanged, so real
# clustering survives normalisation.
test_that("g(r) of a clustered pattern is invariant under p-thinning", {
  skip_if_not_slow()
  W <- 4000; win <- rbind(c(0, 0), c(W, W)); r <- seq(20, 300, by = 20)
  ratio <- replicate(8, {
    P  <- sim_thomas(80, 30, 40, W)
    g0 <- compute_pcf(P[, 1], P[, 2], r_seq = r, window_ranges = win, max_subsample = 1e5)
    k  <- sample(nrow(P), round(0.3 * nrow(P)))
    g1 <- compute_pcf(P[k, 1], P[k, 2], r_seq = r, window_ranges = win, max_subsample = 1e5)
    mean(g1$g[r <= 100]) / mean(g0$g[r <= 100])
  })
  expect_equal(mean(ratio), 1, tolerance = 0.08)
})
