# Local enhancement: k-NN density of each unit divided by the cell's mean density
# (compute_local_enhancement). Dispersion index and eps calibration are covered in
# test-hotspots.R.

# CLAIM: on CSR, pi * lambda * r_k^2 ~ Gamma(k, 1) in any dimension, so the median
# enhancement is k / qgamma(0.5, k) (1.07 at k = 5), not exactly 1.
test_that("median enhancement on CSR matches the k-NN expectation", {
  k <- 5L; expct <- k / stats::qgamma(0.5, k)
  set.seed(1)
  p <- data.frame(x = stats::runif(2000, 0, 5000), y = stats::runif(2000, 0, 5000))
  e <- compute_local_enhancement(p, k = k, support = list(support = 25e6))
  expect_false(e$used_bbox)
  expect_equal(e$lambda_cell, 2000 / 25e6)
  expect_equal(stats::median(e$enh), expct, tolerance = 0.08)
  # 3D: edge effects pull the median down slightly in a cube
  p3 <- data.frame(x = stats::runif(3000, 0, 2000), y = stats::runif(3000, 0, 2000),
                   z = stats::runif(3000, 0, 2000))
  e3 <- compute_local_enhancement(p3, k = k, use_3d = TRUE, support = list(support = 8e9))
  expect_equal(e3$d, 3L)
  expect_equal(stats::median(e3$enh), expct, tolerance = 0.12)
})

test_that("a dense focus stands out against the background", {
  set.seed(7)
  bg <- data.frame(x = stats::runif(1000, 0, 5000), y = stats::runif(1000, 0, 5000))
  fo <- data.frame(x = stats::rnorm(40, 2500, 30), y = stats::rnorm(40, 2500, 30))
  e  <- compute_local_enhancement(rbind(bg, fo), k = 5, support = list(support = 25e6))
  expect_gt(stats::median(e$enh[1001:1040]), 20)
  expect_equal(stats::median(e$enh[1:1000]), 1, tolerance = 0.15)
})

test_that("duplicate coordinates stay finite and a missing footprint is flagged", {
  set.seed(8)
  p  <- data.frame(x = stats::runif(50, 0, 1000), y = stats::runif(50, 0, 1000))
  e  <- compute_local_enhancement(rbind(p, p[1:5, ]), k = 3)
  expect_true(all(is.finite(e$enh)))
  expect_true(e$used_bbox)
  expect_null(compute_local_enhancement(p[1:2, ]))
  # rows with missing coordinates are dropped and reported through `keep`
  p$x[3] <- NA
  e2 <- compute_local_enhancement(p, k = 3)
  expect_equal(sum(e2$keep), 49L)
  expect_length(e2$enh, 49L)
})
