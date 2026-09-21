# continuity_mobility_table() (helpers_tracking.R): D by continuity quartile.
# It must read ~1 when continuity and mobility are independent, and pick up the
# gradient when fast molecules are the ones collecting dark cycles.
TCP_CRLB_PREFACTOR <- if (exists("TCP_CRLB_PREFACTOR")) TCP_CRLB_PREFACTOR else 4
for (f in c("R/est_hmm.R", "R/est_cve.R")) source(file.path(.repo_root, f))

.t_c <- 0.192
Q0 <- matrix(c(-1e-9, 1e-9, 1e-9, -1e-9), 2)             # one effective state
.pop <- function(D, n, gap_p, seed)
  sim_ct_state_tracks(n, 200, D = c(D, D * 1.000001), Q = Q0, t_c = .t_c, sw = 0.011,
                      gap_p = gap_p, sigma = 5, R_blur = 1 / 6.2, seed = seed)
# Same continuity score as the app's gate: mean over steps of 2^(-(k-1)/H).
.scores <- function(df, H = 3) {
  s <- tapply(df$tim, df$tid, function(v) { k <- pmax(1, round(diff(sort(v)) / .t_c))
                                            mean(2^(-(k - 1) / H)) })
  stats::setNames(as.numeric(s), names(s))
}

test_that("continuity independent of mobility reads a ratio near 1", {
  a <- .pop(215, 60, c(1), 1)                               # every step integrated
  b <- .pop(215, 60, c(0.4, 0.3, 0.2, 0.1), 2)              # gappy, SAME D
  b$tid <- paste0("b", b$tid)
  df <- rbind(a, b)
  tab <- continuity_mobility_table(df, .scores(df), .t_c)
  expect_equal(nrow(tab), 4L)
  expect_lt(abs(log(attr(tab, "ratio_best_rest"))), log(1.2))
})

test_that("fast molecules collecting dark cycles show up as a gradient", {
  slow <- .pop(100, 60, c(1), 3)                            # continuous and slow
  fast <- .pop(400, 60, c(0.4, 0.3, 0.2, 0.1), 4)           # gappy and fast
  fast$tid <- paste0("f", fast$tid)
  df <- rbind(slow, fast)
  tab <- continuity_mobility_table(df, .scores(df), .t_c)
  expect_lt(attr(tab, "ratio_best_rest"), 0.6)              # best quarter = the slow ones
  expect_equal(tab$D_um2_s[4], 0.1, tolerance = 0.1)        # ~100 nm²/ms = 0.1 µm²/s
})

test_that("too few tracks returns NULL instead of thin quarters", {
  df <- .pop(215, 20, c(1), 5)
  expect_null(continuity_mobility_table(df, .scores(df), .t_c))
})
