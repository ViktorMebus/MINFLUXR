# Ratiometric two-colour separation: GMM split on per-trace median DCR
# (auto_detect_dcr_split), per-trace confidence (dcr_split_confidence) and channel
# assignment (resolve_channels).

# CLAIM: the split is the crossing point of the two weighted Gaussians, and the
# reported mis-assignment rate is the mixture's own tail mass beyond it.
test_that("the GMM split and error rate match the true mixture", {
  set.seed(3)
  v <- c(stats::rnorm(600, 0.30, 0.05), stats::rnorm(400, 0.60, 0.05))
  g <- auto_detect_dcr_split(v)
  f <- function(t) 0.6 * stats::dnorm(t, 0.3, 0.05) - 0.4 * stats::dnorm(t, 0.6, 0.05)
  sp  <- stats::uniroot(f, c(0.3, 0.6))$root
  err <- 0.6 * stats::pnorm(sp, 0.3, 0.05, lower.tail = FALSE) + 0.4 * stats::pnorm(sp, 0.6, 0.05)
  expect_true(g$success)
  expect_equal(g$split, sp, tolerance = 0.02)
  expect_equal(g$err_rate, err, tolerance = 0.002)      # absolute: both are ~0.1 %
  expect_equal(unname(g$dprime), 6, tolerance = 0.1)
  expect_identical(g$quality, "clean")
  expect_equal(g$n_traces, 1000L)
})

test_that("overlapping dyes are reported as a poor split with more crosstalk", {
  set.seed(4)
  v <- c(stats::rnorm(600, 0.45, 0.08), stats::rnorm(400, 0.55, 0.08))
  g <- auto_detect_dcr_split(v)
  skip_if_not(isTRUE(g$success), "GMM merged the two modes")
  expect_identical(g$quality, "poor")
  expect_gt(g$err_rate, 0.1)
  expect_gt(g$ambiguous_frac, 0.3)
})

test_that("the split refuses too few traces and single-detector data", {
  expect_false(auto_detect_dcr_split(stats::runif(5))$success)
  set.seed(5)
  flat <- auto_detect_dcr_split(0.5 + stats::rnorm(50, 0, 0.001))
  expect_false(flat$success)
  expect_match(flat$reason, "no spread")
})

# CLAIM: the per-trace posterior agrees with the threshold side, sits in [0.5, 1],
# flags the overlap band as `edge`, and stays defined far out in the tails.
test_that("per-trace confidence is consistent with the split", {
  set.seed(6)
  g  <- auto_detect_dcr_split(c(stats::rnorm(600, 0.30, 0.05), stats::rnorm(400, 0.60, 0.05)))
  s  <- g$params$split
  v  <- c(-5, 0.2, s - 0.001, s + 0.001, 0.7, 5, NA)
  cf <- dcr_split_confidence(v, g)
  expect_equal(cf$group, c(1L, 1L, 1L, 2L, 2L, 2L, NA))
  expect_true(all(cf$post[1:6] >= 0.5 & cf$post[1:6] <= 1))
  expect_equal(cf$edge[1:6], cf$post[1:6] < 0.95)
  expect_true(cf$edge[3] && cf$edge[4])                 # right at the split: a coin flip
  expect_true(all(is.finite(cf$log_odds[1:6])))         # no 0/0 in the far tails
  expect_true(all(is.na(dcr_split_confidence(v, list(success = FALSE))$group)))
})

# CLAIM: a trace is one channel or the other, never split; channel 1 is the lower
# DCR mode; on a clean split every trace lands in its dye's channel.
test_that("resolve_channels assigns whole traces to the right channel", {
  set.seed(3)
  mol <- grid_molecules(12, 80)
  dm  <- ifelse(seq_len(nrow(mol)) %% 2 == 0, 0.30, 0.62)
  df  <- sim_smlm_traces(mol, blink_mean = 1, sigma = 5, dcr = dm)
  rc  <- resolve_channels(df, smlm_input(enable_two_color = TRUE))
  tr  <- unique(rc[, c("tid", "mol", "channel")])
  expect_false(any(duplicated(tr$tid)))
  expect_equal(nrow(tr), length(unique(df$tid)))
  expect_equal(mean(tr$channel == ifelse(dm[tr$mol] < 0.45, 1L, 2L)), 1)
  expect_identical(attr(rc, "channel_source"), "two_color")
  expect_equal(attr(rc, "channel_table")$label, c("A", "B"))
  expect_setequal(unique(rc$analysis_group), c(1L, 2L))
})

# CLAIM: when the GMM cannot split, the manual DCR ranges apply, and traces outside
# both ranges are dropped rather than pooled into a channel.
test_that("manual ranges are the fallback and unassigned traces are dropped", {
  df <- data.frame(tid = rep(1:4, each = 3), x = 0, y = 0, z = 0,
                   dcr = rep(c(0.20, 0.25, 0.50, 0.90), each = 3))
  inp <- smlm_input(enable_two_color = TRUE, dcr_group1 = c(0, 0.3), dcr_group2 = c(0.8, 1))
  rc  <- resolve_channels(df, inp)
  expect_setequal(unique(rc$tid), c(1L, 2L, 4L))        # trace 3 (0.50) fits neither range
  expect_equal(unique(rc$channel[rc$tid %in% 1:2]), 1L)
  expect_equal(unique(rc$channel[rc$tid == 4]), 2L)
})

# CLAIM: typed ranges replace a GMM that would have succeeded, on whichever DCR
# column is selected; left off, the GMM decides.
test_that("typed ranges override a successful GMM, on raw and pooled DCR", {
  set.seed(9)
  n <- 40L
  meds <- c(stats::rnorm(n, 0.30, 0.02), stats::rnorm(n, 0.62, 0.02))
  df <- data.frame(tid = rep(seq_len(2L * n), each = 3L), x = 0, y = 0, z = 0,
                   dcr = rep(meds, each = 3L), dcr_pooled = rep(meds + 0.2, each = 3L))
  expect_true(auto_detect_dcr_split(meds)$success)

  chan <- function(...) {
    rc <- resolve_channels(df, smlm_input(enable_two_color = TRUE, ...))
    tapply(rc$channel, rc$tid, `[`, 1L)
  }
  # GMM: the two modes separate cleanly, 40 traces each.
  expect_equal(as.numeric(table(chan())), c(n, n))
  # Typed on raw DCR: this boundary puts everything but the very top in channel 1.
  typed <- list(dcr_split_manual = TRUE, dcr_group1 = c(0, 0.9), dcr_group2 = c(0.9001, 1))
  expect_true(all(do.call(chan, typed) == 1L))
  # The same switch applies to the pooled column. The modes sit at 0.50 and 0.82
  # there, so a boundary at 0.7 separates them; on raw DCR (0.30, 0.62) it would not.
  typed_p <- list(dcr_split_manual = TRUE, use_pooled_dcr = TRUE,
                  dcr_group1 = c(0, 0.7), dcr_group2 = c(0.7001, 1.2))
  got <- do.call(chan, typed_p)
  expect_equal(as.numeric(table(got)), c(n, n))
  expect_equal(as.vector(got[seq_len(n)]), rep(1L, n))
  # Ranges left in the boxes do nothing while the GMM is in charge.
  expect_equal(as.numeric(table(chan(dcr_group1 = c(0, 0.9), dcr_group2 = c(0.9001, 1)))), c(n, n))
  # A half-cleared range is not usable, so the GMM decides.
  expect_equal(as.numeric(table(chan(dcr_split_manual = TRUE, dcr_group1 = c(0, NA),
                                     dcr_group2 = c(0.5, 1)))), c(n, n))
})

test_that("dcr_split_column picks the column the split reads", {
  both <- data.frame(dcr = 1, dcr_pooled = 1)
  expect_identical(dcr_split_column(both, TRUE), "dcr_pooled")
  expect_identical(dcr_split_column(both, FALSE), "dcr")
  expect_identical(dcr_split_column(data.frame(dcr = 1), TRUE), "dcr")
  expect_null(dcr_split_column(data.frame(x = 1), FALSE))
})

test_that("fluo-based channels follow the fluo column", {
  df <- data.frame(tid = rep(1:3, each = 2), x = 0, y = 0, z = 0, fluo = rep(c(1, 2, 2), each = 2))
  rc <- resolve_channels(df, smlm_input(enable_two_color = TRUE, color_by_fluo = TRUE))
  expect_equal(rc$channel, c(1L, 1L, 2L, 2L, 2L, 2L))
})

test_that("with no split everything is channel 1", {
  df <- data.frame(tid = 1:3, x = 1:3, y = 1:3)
  rc <- resolve_channels(df, smlm_input())
  expect_true(all(rc$channel == 1L))
  expect_identical(attr(rc, "channel_source"), "none")
})
