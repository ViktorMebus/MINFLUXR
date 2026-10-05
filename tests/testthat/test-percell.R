# Per-cell comparison (R/est_percell.R) and the timing / cleaning functions it shares
# with the main Tracking pipeline (R/est_clean.R).
for (f in c("R/est_clean.R", "R/est_blur.R", "R/est_cve.R", "R/est_msd.R", "R/est_percell.R")) source(file.path(.repo_root, f))

# Free 2D Brownian tracks sampled the MINFLUX way: each interval spans k cycles,
# dt = k * t_c + sw exactly, with k drawn from k_p. D in nm^2/ms, sigma in nm.
.cell <- function(t_c, sw = 0.011, D = 200, sigma = 5, n_tracks = 60, n_locs = 200,
                  k_p = c(0.7, 0.2, 0.1), seed = 1, prefix = "") {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_tracks), function(j) {
    k  <- sample(seq_along(k_p), n_locs - 1L, replace = TRUE, prob = k_p)
    dt <- k * t_c + sw
    x  <- cumsum(c(0, rnorm(n_locs - 1L, 0, sqrt(2 * D * dt))))
    y  <- cumsum(c(0, rnorm(n_locs - 1L, 0, sqrt(2 * D * dt))))
    data.frame(tid = paste0(prefix, j), tim = cumsum(c(0, dt)),
               x = x + rnorm(n_locs, 0, sigma), y = y + rnorm(n_locs, 0, sigma),
               eco = 30, stringsAsFactors = FALSE)
  }))
}

test_that("cycle_cluster recovers the cycle time and the switching overhead", {
  d  <- .cell(t_c = 0.166, sw = 0.011)
  cc <- cycle_cluster(d$tid, d$tim)
  expect_equal(cc$cycle, 0.166, tolerance = 1e-3)
  expect_equal(cc$sw, 0.011, tolerance = 0.05)
  expect_equal(cc$frac_one_cycle, 0.7, tolerance = 0.05)
})

test_that("split_gaps splits at the threshold, drops short pieces, and is idempotent", {
  d <- data.frame(tid = "a", tim = c(0, 1, 2, 10, 11, 12, 13), x = 0, y = 0)
  s <- split_gaps(d, gap_thresh_ms = 5, min_len = 1L)
  expect_equal(unique(s$tid), c("a_s1", "a_s2"))
  expect_identical(split_gaps(s, 5, 1L), s)                       # splitting twice == once
  expect_equal(unique(split_gaps(d, 5, min_len = 4L)$tid), "a_s2")  # 3-loc piece dropped
})

test_that("adaptive_gap_cycles is median + SD + 1, capped by maxOffTime", {
  tot <- c(rep(1L, 80), rep(2L, 15), rep(3L, 5))
  expect_equal(adaptive_gap_cycles(tot, cyc = 0.2, maxoff_ms = 3),
               as.integer(floor(1 + stats::sd(tot) + 1)))
  expect_equal(adaptive_gap_cycles(c(rep(1L, 10), rep(40L, 10)), cyc = 0.2, maxoff_ms = 1), 5L)
})

test_that("each cell gets its own cycle time, D and sequence-specific sigma", {
  filt <- list(split = list(enable = TRUE, adaptive = FALSE, gap_cycles = 7, min_len = 20L, maxoff_ms = 3))
  slow <- percell_summary(.cell(t_c = 0.166, seed = 2), "slow", L = 100.8, a = 3.6, filters = filt)
  fast <- percell_summary(.cell(t_c = 0.128, seed = 3), "fast", L = 151.2, a = 3.6, filters = filt)
  expect_equal(slow$t_c_ms, 0.166, tolerance = 1e-3)
  expect_equal(fast$t_c_ms, 0.128, tolerance = 1e-3)
  expect_equal(slow$D_um2_s, 0.2, tolerance = 0.1)                # true D 200 nm^2/ms
  expect_equal(fast$D_um2_s, 0.2, tolerance = 0.1)
  # same photons: sigma scales with L (up to the small motion term)
  expect_equal(fast$sigma_nm / slow$sigma_nm, 151.2 / 100.8, tolerance = 0.03)
  expect_equal(slow$note, "")
  # D_w is reported beside D, never folded into it
  w <- percell_summary(.cell(t_c = 0.166, seed = 2), "w", L = 100.8, a = 3.6, D_w = 50, filters = filt)
  expect_equal(w$D_um2_s, slow$D_um2_s)
  expect_equal(w$D_corr_um2_s, slow$D_um2_s - 0.05)
})

test_that("an empty or too-short file returns a row with a note", {
  r <- percell_summary(data.frame(tid = character(0), x = numeric(0), y = numeric(0), tim = numeric(0)),
                       "empty", L = 100.8)
  expect_equal(nrow(r), 1L); expect_equal(r$note, "no localizations")
})

test_that("condition stats use the cell as the unit and skip unassigned cells", {
  tab <- data.frame(condition = c("A", "A", "A", "B", "B", ""), D = c(1, 2, 3, 10, 12, 99))
  st  <- percell_condition_stats(tab, "D")
  expect_equal(st$condition, c("A", "B"))
  expect_equal(st$n_cells, c(3L, 2L))
  expect_equal(st$mean, c(2, 11))
  expect_equal(st$ci_hi[1] - st$mean[1], stats::qt(0.975, 2) * 1 / sqrt(3))
})
