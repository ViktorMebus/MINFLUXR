# Instrument calibration from an immobile sample (R/est_calib.R).
for (f in c("R/est_clean.R", "R/est_calib.R")) source(file.path(.repo_root, f))

# Immobile emitters: the reported position is a slow random walk (the wander,
# D_w in nm^2/ms) plus white noise with sigma_i = L / (a sqrt(N_i)). Intervals are
# k cycles plus the switching overhead, dt = k * t_c + sw.
.immobile <- function(a = 3.9, L = 100.8, D_w = 0, t_c = 0.166, sw = 0.011, n_tracks = 80,
                      n_locs = 300, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_tracks), function(j) {
    k   <- sample(1:3, n_locs - 1L, replace = TRUE, prob = c(0.7, 0.2, 0.1))
    dt  <- k * t_c + sw
    eco <- stats::rpois(n_locs, 30) + 5
    sig <- L / (a * sqrt(eco))
    w   <- cumsum(c(0, rnorm(n_locs - 1L, 0, sqrt(2 * D_w * dt))))
    v   <- cumsum(c(0, rnorm(n_locs - 1L, 0, sqrt(2 * D_w * dt))))
    data.frame(tid = as.character(j), tim = cumsum(c(0, dt)), eco = eco,
               x = w + rnorm(n_locs, 0, sig), y = v + rnorm(n_locs, 0, sig))
  }))
}

test_that("calib_cell recovers a from white localization noise", {
  r <- calib_cell(.immobile(a = 3.9), "sim", L = 100.8)
  expect_equal(r$note, "")
  expect_equal(r$a, 3.9, tolerance = 0.05)
  expect_lt(abs(r$D_w), 5)                                        # no wander simulated
})

test_that("calib_cell separates a from a diffusive wander floor", {
  r <- calib_cell(.immobile(a = 3.9, D_w = 150, seed = 2), "sim", L = 100.8)
  expect_equal(r$a, 3.9, tolerance = 0.08)
  expect_equal(r$D_w, 150, tolerance = 0.15)
})

test_that("calib_cell scales with L, and reports why it cannot measure", {
  r <- calib_cell(.immobile(a = 3.0, L = 151.2, seed = 3), "sim", L = 151.2)
  expect_equal(r$a, 3.0, tolerance = 0.05)
  bad <- calib_cell(data.frame(tid = "1", tim = 1:3, x = 0, y = 0), "x", L = 100)
  expect_match(bad$note, "eco")
})

test_that("calib_summary takes the median and quartiles across cells", {
  sm <- calib_summary(data.frame(a = c(3, 4, 5, NA), D_w = c(100, 200, 300, 400)))
  expect_equal(sm$n_cells, 3L)
  expect_equal(c(sm$a_q25, sm$a, sm$a_q75), c(3.5, 4, 4.5))
  expect_equal(sm$D_w, 250)
})
