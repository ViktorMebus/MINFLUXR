# De-blinking: complete linkage on trace centroids (helpers_smlm.R, section
# "deblink_*" and apply_clustering_pipeline). Ground truth from sim_smlm_traces().

# CLAIM: link_complete_threshold() gives exactly the partition of a full
# hclust(complete) cut at the threshold, while only clustering inside connected
# components of the threshold graph.
test_that("link_complete_threshold matches hclust complete linkage exactly", {
  set.seed(4)
  M   <- matrix(stats::runif(600, 0, 200), ncol = 2)
  got <- link_complete_threshold(M, 20)
  ref <- stats::cutree(stats::hclust(stats::dist(M), "complete"), h = 20)
  tab <- table(got, ref)
  expect_true(all(rowSums(tab > 0) == 1))            # every group maps to one reference group
  expect_true(all(colSums(tab > 0) == 1))            # and back
  expect_true(all(got > 0))
})

# CLAIM: complete linkage bounds a molecule's diameter by the threshold, so an
# over-large threshold cannot chain neighbours into one long object.
test_that("no de-blinked group is wider than the threshold", {
  set.seed(11)
  df  <- sim_smlm_traces(grid_molecules(8, 25, jitter = 5), blink_mean = 2, sigma = 5)
  thr <- 20
  cl  <- deblink_merge_centroids(df, thr)
  cen <- stats::aggregate(cbind(x, y) ~ tid, data = df, FUN = stats::median)
  cen$cl <- cl[match(cen$tid, df$tid)]
  dia <- tapply(seq_len(nrow(cen)), cen$cl, function(ix)
    if (length(ix) > 1L) max(stats::dist(cen[ix, c("x", "y")])) else 0)
  expect_lte(max(dia), thr + 1e-9)
})

# CLAIM: on molecules spaced well above the trace spread, centroid de-blinking at the
# suggested threshold (k = 2.8 x sigma) merges every re-blink and fuses no molecules.
test_that("de-blinking recovers the molecule count on separable data", {
  set.seed(1)
  mol <- grid_molecules(10, 60, jitter = 10)
  df  <- sim_smlm_traces(mol, blink_mean = 1, sigma = 5, wobble = 1)
  sg  <- suggest_centroid_threshold(df)
  expect_equal(sg$sigma, 5, tolerance = 0.1)         # median per-axis trace spread
  expect_equal(sg$suggested, round(2.8 * sg$sigma, 1))
  cl <- deblink_merge_centroids(df, sg$suggested)
  tr <- unique(data.frame(mol = df$mol, cl = cl))
  s  <- score_grouping(tr$mol, tr$cl)
  expect_equal(unname(s["not_split"]), 1)
  expect_equal(unname(s["not_merged"]), 1)
  expect_equal(unname(s["n_label"]), nrow(mol))
})

# CLAIM: de-blinking merges and never removes. Every trace that enters leaves with a
# cluster id > 0; a trace nothing merges with becomes its own molecule.
test_that("the pipeline keeps every trace and leaves no cluster 0", {
  set.seed(2)
  df  <- sim_smlm_traces(grid_molecules(6, 80), blink_mean = 1, sigma = 5)
  out <- apply_clustering_pipeline(df, smlm_input())
  expect_equal(nrow(out), nrow(df))
  expect_true(all(out$cluster > 0))
  expect_setequal(unique(out$tid), unique(df$tid))
  # one trace, one molecule: the plurality step must not split a trace
  expect_true(all(tapply(out$cluster, out$tid, function(v) length(unique(v))) == 1L))
})

test_that("deblink_assign_singletons gives one new id per noise trace", {
  df <- data.frame(tid = c(1, 1, 2, 2, 3, 3), cluster = c(5L, 5L, 0L, 0L, NA, NA))
  out <- deblink_assign_singletons(df)
  expect_true(all(out$cluster > 0))
  expect_equal(length(unique(out$cluster)), 3L)
  expect_equal(out$cluster[1], 5L)                   # existing ids untouched
  expect_equal(out$cluster[3], out$cluster[4])       # one id per trace, not per point
})

# CLAIM: the per-axis weights put x, y and z on one precision scale and are
# volume-preserving (geometric mean 1), so eps keeps its nm meaning.
test_that("deblink_axis_weights whiten an anisotropic spread and preserve volume", {
  set.seed(5)
  mol <- cbind(grid_molecules(6, 80), stats::runif(36, 0, 300))
  w3  <- deblink_axis_weights(sim_smlm_traces(mol, sigma = c(6, 6, 3)))
  expect_true(w3$ok && w3$has_z)
  expect_equal(prod(w3$w)^(1 / 3), 1, tolerance = 1e-9)
  expect_equal(unname(w3$w[3] / w3$w[1]), 2, tolerance = 0.15)   # z twice as tight
  w2 <- deblink_axis_weights(sim_smlm_traces(grid_molecules(6, 80), sigma = 5))
  expect_false(w2$has_z)
  expect_equal(w2$w[3], 1)                           # 2D data: z left alone
  expect_equal(unname(w2$w[1:2]), c(1, 1), tolerance = 0.1)
})

# CLAIM: the merge distance is k x the trace spread PER CHANNEL, because spread is a
# property of the dye; k is shared.
test_that("each channel gets its own merge distance from its own spread", {
  set.seed(6)
  a <- sim_smlm_traces(grid_molecules(6, 80), sigma = 4);  a$channel <- 1L
  b <- sim_smlm_traces(grid_molecules(6, 80), sigma = 8);  b$channel <- 2L
  b$tid <- b$tid + max(a$tid)
  thr <- deblink_threshold_for_channel(rbind(a, b), smlm_input())
  expect_true(all(attr(thr, "measured")))
  expect_equal(unname(attr(thr, "sigma")), c(4, 8), tolerance = 0.12)
  expect_equal(thr[["2"]] / thr[["1"]], 2, tolerance = 0.15)
  shared <- deblink_threshold_for_channel(rbind(a, b), smlm_input(deblink_shared_spread = TRUE))
  expect_equal(unlist(shared), c(`1` = 10, `2` = 10))
  expect_false(any(attr(shared, "measured")))        # fallback is reported, not hidden
})

# CLAIM: recovery degrades gracefully as molecules move closer, and the error is
# merging (counts fall), never fragmentation of a molecule's own re-blinks.
test_that("re-blinks are never split across a spacing sweep", {
  skip_if_not_slow()
  for (sp in c(60, 40, 25)) {
    r <- replicate(5, {
      df <- sim_smlm_traces(grid_molecules(10, sp, jitter = 3), blink_mean = 1, sigma = 5)
      out <- apply_clustering_pipeline(df, smlm_input())
      tr <- unique(out[, c("mol", "cluster")])
      score_grouping(tr$mol, tr$cluster)[["not_split"]]
    })
    expect_gte(mean(r), 0.98)
  }
})
