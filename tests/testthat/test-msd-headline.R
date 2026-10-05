# The cycle-binned MSD (msd_track_sums / msd_ensemble), the headline upper-regime D
# (msd_headline).
TCP_CRLB_PREFACTOR <- if (exists("TCP_CRLB_PREFACTOR")) TCP_CRLB_PREFACTOR else sqrt(8)
for (f in c("R/est_clean.R", "R/est_blur.R", "R/est_cve.R", "R/est_msd.R", "R/est_calib.R"))
  source(file.path(.repo_root, f))

# Free 2D Brownian tracks, dt = k * t_c + sw with dark cycles; D in nm^2/ms.
.free <- function(D = 200, t_c = 0.2, sw = 0.01, n_tracks = 40, n_locs = 120, sigma = 4, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_tracks), function(j) {
    k  <- sample(1:3, n_locs - 1L, replace = TRUE, prob = c(0.7, 0.2, 0.1))
    dt <- k * t_c + sw
    data.frame(tid = sprintf("t%02d", j), tim = cumsum(c(0, dt)),
               x = cumsum(c(0, rnorm(n_locs - 1L, 0, sqrt(2 * D * dt)))) + rnorm(n_locs, 0, sigma),
               y = cumsum(c(0, rnorm(n_locs - 1L, 0, sqrt(2 * D * dt)))) + rnorm(n_locs, 0, sigma))
  }))
}

test_that("vectorised MSD sums equal a pair-by-pair count", {
  d <- .free(n_tracks = 6, n_locs = 40); t_c <- 0.2; sw <- 0.01; L <- 25L
  s <- msd_track_sums(d, t_c, sw, max_lag = L, min_pts = 5L)
  # naive: every pair within a track, lag = difference of cumulative cycle index
  ref_dr2 <- matrix(0, L, length(s$tk_names)); ref_np <- matrix(0L, L, length(s$tk_names))
  for (k in seq_along(s$tk_names)) {
    tb  <- d[d$tid == s$tk_names[k], ]; tb <- tb[order(tb$tim), ]
    cyc <- c(0L, cumsum(pmax(1L, as.integer(round((diff(tb$tim) - sw) / t_c)))))
    for (i in seq_len(nrow(tb) - 1L)) for (j in (i + 1L):nrow(tb)) {
      g <- cyc[j] - cyc[i]
      if (g >= 1L && g <= L) {
        ref_dr2[g, k] <- ref_dr2[g, k] + (tb$x[j] - tb$x[i])^2 + (tb$y[j] - tb$y[i])^2
        ref_np[g, k]  <- ref_np[g, k] + 1L
      }
    }
  }
  expect_equal(s$dr2, ref_dr2)
  expect_identical(s$np, ref_np)
})

test_that("tracks shorter than min_pts contribute nothing", {
  d <- rbind(.free(n_tracks = 2, n_locs = 30), transform(.free(n_tracks = 1, n_locs = 4, seed = 2), tid = "short"))
  s <- msd_track_sums(d, 0.2, 0.01, max_lag = 10L, min_pts = 10L)
  expect_equal(sum(s$np[, s$tk_names == "short"]), 0L)
})

test_that("the headline D recovers D from the upper regime", {
  d  <- .free(D = 200, n_tracks = 80, n_locs = 200, seed = 3)
  hd <- msd_headline(d, t_c = 0.2, sw = 0.01, window = c(5L, 50L), R = 0, coef = 4, min_pts = 50L)
  expect_equal(hd$D, 200, tolerance = 0.1)
  expect_equal(hd$window, c(5L, 50L))
})
