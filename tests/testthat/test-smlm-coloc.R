# Two-channel co-localization: bivariate L12 / g12, cross-channel nearest
# neighbour, and the KDE Pearson test. The L12/g12 nulls (shift, labels, CSR) are
# covered in test-spatial-envelope.R.
.W   <- 4000
.win <- rbind(c(0, 0), c(.W, .W))
.r   <- seq(20, 400, by = 20)
.thomas_on <- function(par, n, s = 40) {
  k <- sample(nrow(par), n, TRUE)
  par[k, , drop = FALSE] + matrix(stats::rnorm(2 * n, 0, s), ncol = 2)
}

test_that("L12 is ~0 and g12 ~1 for independent CSR channels", {
  set.seed(1)
  a <- cbind(stats::runif(800, 0, .W), stats::runif(800, 0, .W))
  b <- cbind(stats::runif(600, 0, .W), stats::runif(600, 0, .W))
  L <- compute_ripley_L12(a[, 1], a[, 2], NULL, b[, 1], b[, 2], NULL, r_seq = .r,
                          window_ranges = .win, max_subsample = 5000)
  g <- compute_pcf12(a[, 1], a[, 2], NULL, b[, 1], b[, 2], NULL, r_seq = .r,
                     window_ranges = .win, max_subsample = 5000)
  expect_lt(max(abs(L$Lr_minus_r[.r <= 200])), 10)
  expect_equal(mean(g$g12), 1, tolerance = 0.1)        # uncorrected edges read slightly low
})

test_that("channels sharing clusters give g12 >> 1 at short range", {
  set.seed(2)
  par <- cbind(stats::runif(40, 0, .W), stats::runif(40, 0, .W))
  A <- .thomas_on(par, 500); B <- .thomas_on(par, 500)
  g <- compute_pcf12(A[, 1], A[, 2], NULL, B[, 1], B[, 2], NULL, r_seq = .r,
                     window_ranges = .win, max_subsample = 5000)
  L <- compute_ripley_L12(A[, 1], A[, 2], NULL, B[, 1], B[, 2], NULL, r_seq = .r,
                          window_ranges = .win, max_subsample = 5000)
  expect_gt(g$g12[2], 5)                                # r = 40 nm
  expect_gt(L$Lr_minus_r[5], 100)                       # r = 100 nm
})

# CLAIM: the cross K function counts each A-B pair once, so it is symmetric in the
# two channels.
test_that("K12 is symmetric in the two channels", {
  set.seed(3)
  a <- cbind(stats::runif(300, 0, .W), stats::runif(300, 0, .W))
  b <- cbind(stats::runif(200, 0, .W), stats::runif(200, 0, .W))
  K1 <- compute_ripley_L12(a[, 1], a[, 2], NULL, b[, 1], b[, 2], NULL, r_seq = .r, window_ranges = .win)
  K2 <- compute_ripley_L12(b[, 1], b[, 2], NULL, a[, 1], a[, 2], NULL, r_seq = .r, window_ranges = .win)
  expect_equal(K1$Kr, K2$Kr)
})

test_that("cross-channel NN distances match a hand calculation", {
  d <- data.frame(tid = 1:4, x = c(0, 100, 30, 100), y = c(0, 0, 40, 300), z = 0,
                  analysis_group = c(1L, 1L, 2L, 2L))
  nn <- compute_cross_nn(d)
  expect_equal(nn$distance[nn$direction == "Group 1 → Group 2"], c(50, sqrt(70^2 + 40^2)))
  expect_equal(nn$distance[nn$direction == "Group 2 → Group 1"], c(50, 300))
  expect_equal(nrow(compute_cross_nn(d[d$analysis_group == 1L, ])), 0L)
})

# CLAIM: the KDE Pearson r is high for channels sharing structure and ~0 for
# independent structure, and the shift-null p-value separates the two cases.
test_that("KDE co-localization separates shared from independent structure", {
  set.seed(2)
  par  <- cbind(stats::runif(30, 500, 3500), stats::runif(30, 500, 3500))
  par2 <- cbind(stats::runif(30, 500, 3500), stats::runif(30, 500, 3500))
  A <- .thomas_on(par, 400, 50); B <- .thomas_on(par, 400, 50); C <- .thomas_on(par2, 400, 50)
  same <- compute_kde_coloc(A[, 1], A[, 2], A[, 1], A[, 2], bandwidth = 80)
  ab   <- compute_kde_coloc(A[, 1], A[, 2], B[, 1], B[, 2], bandwidth = 80, costes_nsim = 49)
  ac   <- compute_kde_coloc(A[, 1], A[, 2], C[, 1], C[, 2], bandwidth = 80, costes_nsim = 49)
  expect_equal(same$pearson_r, 1)
  expect_gt(ab$pearson_r, 0.4)
  expect_lt(abs(ac$pearson_r), 0.2)
  expect_lt(ab$costes_p, 0.05)
  expect_gt(ac$costes_p, 0.2)
  expect_identical(ab$null_mode, "shift")
  expect_false(ab$null_support)
})

test_that("the KDE null randomises inside a footprint when one is given", {
  set.seed(4)
  par <- cbind(stats::runif(20, 500, 3500), stats::runif(20, 500, 3500))
  A <- .thomas_on(par, 300, 50); B <- .thomas_on(par, 300, 50)
  sup <- estimate_support(data.frame(x = c(A[, 1], B[, 1]), y = c(A[, 2], B[, 2]), z = 0),
                          mode = "grid", h_nm = 100, close_nm = 400, return_mask = TRUE)
  k <- compute_kde_coloc(A[, 1], A[, 2], B[, 1], B[, 2], bandwidth = 80, costes_nsim = 9,
                         support = sup)
  expect_true(k$null_support)
  expect_match(k$null_label, "sample support")
  sw <- compute_kde_coloc(A[, 1], A[, 2], B[, 1], B[, 2], bandwidth = 80, bw_sweep = c(40, 80, 160))
  expect_equal(sw$bw_sweep$bw, c(40, 80, 160))
  expect_equal(sw$bw_sweep$pearson_r[2], sw$pearson_r)  # the sweep reuses the same estimator
})

# CLAIM (methods doc, KDE caveat): with the footprint, the shift null is conservative
# for independent clustered channels sharing one cell inside a larger window.
test_that("KDE shift null with a footprint does not over-reject", {
  skip_if_not_slow()
  set.seed(10)
  cell <- function(n) {                                 # a 2 x 2 um cell in a 4 x 4 um field
    x <- stats::runif(n, 1000, 3000); y <- stats::runif(n, 1000, 3000); cbind(x, y)
  }
  p <- replicate(20, {
    A <- .thomas_on(cell(15), 250, 50); B <- .thomas_on(cell(15), 250, 50)
    sup <- estimate_support(data.frame(x = c(A[, 1], B[, 1]), y = c(A[, 2], B[, 2]), z = 0),
                            mode = "grid", h_nm = 100, close_nm = 400, return_mask = TRUE)
    compute_kde_coloc(A[, 1], A[, 2], B[, 1], B[, 2], bandwidth = 80, costes_nsim = 39,
                      support = sup)$costes_p
  })
  expect_lte(mean(p < 0.05), 0.15)
})

# CLAIM: the L12 / g12 panels' per-file computation (cross_curve_for_points) is the
# estimator on the two channels' shared bounding box, plus the envelope under the
# chosen null when asked.
test_that("cross_curve_for_points matches the estimators and carries the null", {
  set.seed(9)
  a <- data.frame(x = stats::runif(300, 0, 4000), y = stats::runif(300, 0, 4000), z = 0)
  b <- data.frame(x = stats::runif(200, 0, 4000), y = stats::runif(200, 0, 4000), z = 0)
  win <- apply(as.matrix(rbind(a[, 1:2], b[, 1:2])), 2, range); r <- seq(50, 400, by = 50)
  L <- cross_curve_for_points(a, b, "L12", r_seq = r, tag = "f")
  expect_identical(names(L$curve)[1], "dataset")
  expect_equal(L$curve$Lr_minus_r,
               compute_ripley_L12(a$x, a$y, NULL, b$x, b$y, NULL, r_seq = r, window_ranges = win)$Lr_minus_r)
  g <- cross_curve_for_points(a, b, "g12", r_seq = r, tag = "f")
  expect_equal(g$curve$g12,
               compute_pcf12(a$x, a$y, NULL, b$x, b$y, NULL, r_seq = r, window_ranges = win)$g12)
  set.seed(10)
  e <- cross_curve_for_points(a, b, "L12", r_seq = r, envelope = TRUE, nsim = 19,
                              null_mode = "labels", tag = "f")
  expect_true(all(e$env$null_mode == "labels"))
  expect_true(all(c("global_p", "dataset") %in% names(e$env)))
  expect_null(cross_curve_for_points(a[0, ], b, "L12", r_seq = r))
})
