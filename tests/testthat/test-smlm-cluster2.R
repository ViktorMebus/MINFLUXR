# Median clustering (cluster2): run_cluster2_grouping() and its use per channel in
# apply_clustering_pipeline(). Four Gaussian blobs on a sparse uniform background.
.c2_scene <- function(seed = 6) {
  set.seed(seed)
  cent  <- cbind(c(500, 1500, 2500, 1000), c(500, 1500, 500, 2500))
  blobs <- do.call(rbind, lapply(1:4, function(i)
    cbind(stats::rnorm(40, cent[i, 1], 60), stats::rnorm(40, cent[i, 2], 60), 0)))
  bg <- cbind(stats::runif(60, 0, 3000), stats::runif(60, 0, 3000), 0)
  list(P = rbind(blobs, bg), truth = c(rep(1:4, each = 40), rep(0L, 60)))
}

# CLAIM (output contract): 0 = noise, survivors numbered 1..k contiguously, whatever
# the method and refinements.
test_that("labels are 0 for noise and contiguous 1..k otherwise", {
  s <- .c2_scene()
  for (args in list(list(method = "dbscan"), list(method = "hdbscan"),
                    list(method = "optics"), list(method = "dbscan", min_members = 30),
                    list(method = "dbscan", prune = TRUE, core_only = TRUE))) {
    lab <- do.call(run_cluster2_grouping, c(list(s$P, eps = 120, minPts = 8), args))
    expect_length(lab, nrow(s$P))
    k <- max(lab)
    if (k > 0) expect_setequal(unique(lab[lab > 0]), seq_len(k))
  }
})

test_that("DBSCAN finds the four blobs and calls most background noise", {
  s   <- .c2_scene()
  lab <- run_cluster2_grouping(s$P, method = "dbscan", eps = 120, minPts = 8)
  expect_equal(max(lab), 4L)
  expect_true(all(lab[s$truth > 0] > 0))
  expect_gte(mean(lab[s$truth == 0] == 0), 0.9)
  # one label per blob, no label shared between blobs
  tab <- table(s$truth[s$truth > 0], lab[s$truth > 0])
  expect_true(all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1))
})

test_that("HDBSCAN finds the four blobs without an eps", {
  s   <- .c2_scene()
  lab <- run_cluster2_grouping(s$P, method = "hdbscan", minPts = 8)
  expect_equal(max(lab), 4L)
  expect_true(all(lab[s$truth > 0] > 0))
})

# CLAIM: OPTICS adapts to density and may split a blob into sub-structures, but it
# never merges two separate blobs into one cluster.
test_that("OPTICS never merges two separate blobs", {
  s   <- .c2_scene()
  lab <- run_cluster2_grouping(s$P, method = "optics", minPts = 8, xi = 0.05)
  inb <- s$truth > 0 & lab > 0
  tab <- table(lab[inb], s$truth[inb])
  expect_true(all(rowSums(tab > 0) == 1))
})

# CLAIM: the local-density prefilter forces isolated background centroids to noise
# before clustering, without losing the structures.
test_that("the prefilter removes background and keeps the blobs", {
  s   <- .c2_scene()
  lab <- run_cluster2_grouping(s$P, eps = 120, minPts = 8, prefilter = TRUE,
                               prefilter_r = 120, prefilter_min = 5)
  expect_equal(max(lab), 4L)
  expect_equal(mean(lab[s$truth == 0] == 0), 1)
})

# CLAIM: core-only membership (DBSCAN*) drops border points and nothing else, so it
# can only shrink a cluster.
test_that("core_only drops border members only", {
  s  <- .c2_scene()
  l0 <- run_cluster2_grouping(s$P, eps = 120, minPts = 8)
  l1 <- run_cluster2_grouping(s$P, eps = 120, minPts = 8, core_only = TRUE)
  expect_lte(sum(l1 > 0), sum(l0 > 0))
  expect_true(all(l0[l1 > 0] > 0))                   # a core member was a member before
  nn <- dbscan::frNN(s$P, eps = 120)
  expect_true(all(lengths(nn$id)[l1 > 0] + 1L >= 8L))
})

test_that("min_members drops small clusters and renumbers", {
  s    <- .c2_scene()
  keep <- c(which(s$truth != 4), which(s$truth == 4)[1:12])   # blob 4 cut to 12 points
  lab  <- run_cluster2_grouping(s$P[keep, ], eps = 120, minPts = 8, min_members = 20)
  expect_equal(max(lab), 3L)
  expect_true(all(lab[s$truth[keep] == 4] == 0))
})

# CLAIM: pruning judges each member against its OWN neighbourhood, so a uniform ring
# is kept whole while a straggler hanging off it is dropped.
test_that("pruning drops a loose straggler and keeps a uniform ring", {
  set.seed(8)
  th   <- stats::runif(80, 0, 2 * pi)
  ring <- cbind(300 * cos(th), 300 * sin(th)) + matrix(stats::rnorm(160, 0, 8), ncol = 2)
  P    <- rbind(ring, c(300 + 110, 0))               # one point loosely attached
  lab0 <- run_cluster2_grouping(P, eps = 120, minPts = 4)
  expect_gt(lab0[81], 0)                             # plain DBSCAN absorbs it
  lab  <- run_cluster2_grouping(P, eps = 120, minPts = 4, prune = TRUE, prune_tol = 2)
  expect_equal(lab[81], 0L)
  expect_gte(mean(lab[1:80] > 0), 0.95)
})

# CLAIM: grouping never crosses a channel; a Median Cluster holds one channel only,
# even where two channels' molecules interleave in space.
test_that("cluster2 runs strictly per channel", {
  set.seed(9)
  mol  <- grid_molecules(12, 50)
  dcr  <- ifelse(seq_len(nrow(mol)) %% 2 == 0, 0.30, 0.62)
  df   <- sim_smlm_traces(mol, blink_mean = 1, sigma = 5, dcr = dcr)
  out  <- apply_clustering_pipeline(df, smlm_input(enable_two_color = TRUE,
                                                   enable_cluster2 = TRUE, eps2 = 120, minpts2 = 4))
  c2   <- out[!is.na(out$cluster2) & out$cluster2 > 0, ]
  expect_gt(length(unique(c2$cluster2)), 1L)
  expect_true(all(tapply(c2$channel, c2$cluster2, function(v) length(unique(v))) == 1L))
})
