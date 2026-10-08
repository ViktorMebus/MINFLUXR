# Inverse ground-truth estimator (estimate_gt_from_labeling), flagged EXPLORATORY
# in the app, and the trace-count model under it (gt_tid_model). Single-dye label
# at 100 % efficiency unless stated, so one resolved molecule = one de-blinked unit.
.gt_preset <- list(n_fluor = 1, link_nm = 5, dye_spread_nm = 0, footprint_nm = 4,
                   blink_lambda = 1)
.gt_run <- function(mol, seed, ...) {
  set.seed(seed)
  df  <- apply_clustering_pipeline(sim_smlm_traces(mol, blink_mean = 1, sigma = 5), smlm_input())
  estimate_gt_from_labeling(df, .gt_preset, eff_pct = 100, scope_cluster2 = FALSE,
                            displace_inward = FALSE, ...)
}

# CLAIM: gt_tid_model is the exact distribution of the Simulator's labelling model
# (max(1, Poisson) dyes, binomial detection, 1 + Poisson blinks per detected dye).
test_that("the trace-count model matches a Monte Carlo of the forward model", {
  set.seed(1)
  n <- 2e5
  for (a in list(c(1, 0.7, 0.3), c(3, 0.7, 0.5), c(9, 0.5, 0.5))) {
    K  <- if (a[1] <= 1) rep(1L, n) else pmax(1L, stats::rpois(n, a[1]))
    Kd <- stats::rbinom(n, K, a[2]); Tm <- Kd + stats::rpois(n, Kd * a[3])
    m  <- gt_tid_model(a[1], a[2], a[3])
    expect_equal(m$eta, mean(Kd > 0), tolerance = 0.01)
    expect_equal(m$mean_u, mean(Tm), tolerance = 0.01)
    expect_equal(m$var_u, stats::var(Tm), tolerance = 0.03)
    expect_equal(m$mean_d, mean(Tm[Kd > 0]), tolerance = 0.01)
    expect_equal(m$var_d, stats::var(Tm[Kd > 0]), tolerance = 0.03)
    expect_equal(m$p_ge(4), mean(Tm[Kd > 0] >= 4), tolerance = 0.03)
  }
  # single dye, no blinking: exactly one trace per detected molecule
  m0 <- gt_tid_model(1, 0.8, 0)
  expect_equal(c(m0$eta, m0$mean_u, m0$mean_d, m0$var_d), c(0.8, 0.8, 1, 0))
  expect_equal(m0$p_ge(c(1, 2)), c(1, 0))
})

test_that("resolved molecules are counted, placed and pass the consistency test", {
  set.seed(1)
  mol <- grid_molecules(8, 100, jitter = 15)
  gt  <- .gt_run(mol, 1)
  s   <- gt$summary
  expect_identical(gt$meta$level, "primary")            # de-blinked units are preferred
  expect_equal(s$n_gt_est, nrow(mol), tolerance = 0.05)
  expect_equal(s$n_stoich, nrow(mol), tolerance = 0.15)
  expect_identical(s$count_verdict, "consistent")
  expect_false(s$crowded || s$stoich_mismatch)
  expect_false(any(gt$points$fused))
  expect_lte(s$n_gt_lo, nrow(mol)); expect_gte(s$n_gt_hi, nrow(mol))
  d <- RANN::nn2(mol, as.matrix(gt$points[, c("x", "y")]), k = 1)$nn.dists[, 1]
  expect_lt(stats::median(d), 5)                        # candidates sit on molecules
})

# CLAIM: the linker offset is common-mode within a molecule, so no candidate is
# reported more precise than link_nm.
test_that("position uncertainty never falls below the linker floor", {
  gt <- .gt_run(grid_molecules(8, 100, jitter = 15), 2)
  expect_true(all(gt$points$pos_sd_nm >= .gt_preset$link_nm))
})

# CLAIM: when molecules sit closer than the precision, de-blinking merges them into
# units no wider than its threshold, so the spread-based fusion test cannot see it.
# The trace total can: the structure is called under-resolved, `crowded` fires, and
# the range (both estimators with their errors) still contains the truth.
test_that("unresolved molecules are caught by the trace-count test", {
  mol <- grid_molecules(8, 6)
  for (seed in 1:3) {
    gt <- .gt_run(mol, seed)
    s  <- gt$summary
    expect_lt(s$n_gt_est, 0.6 * nrow(mol))              # spatial count falls short
    expect_identical(s$count_verdict, "under-resolved")
    expect_gt(s$count_z, 2.58)
    expect_true(s$crowded && s$stoich_mismatch)
    expect_lte(s$n_gt_lo, nrow(mol)); expect_gte(s$n_gt_hi, nrow(mol))
    expect_gt(sum(gt$points$count_fused), 0)            # and the heavy candidates are marked
  }
})

test_that("re-blinks left as separate candidates are called over-split", {
  set.seed(4)
  mol <- grid_molecules(8, 100, jitter = 15)
  df  <- sim_smlm_traces(mol, blink_mean = 1, sigma = 5)
  gt  <- estimate_gt_from_labeling(df, .gt_preset, eff_pct = 100, level = "tid",
                                   complex_r = 0, scope_cluster2 = FALSE,
                                   displace_inward = FALSE)
  s <- gt$summary
  expect_gt(s$n_gt_est, 1.5 * nrow(mol))                # every re-blink counted
  expect_identical(s$count_verdict, "over-split")
  expect_false(s$crowded)
  expect_lte(s$n_gt_lo, nrow(mol))
})

# CLAIM: a measured multiplicity replaces the preset mean, and the blink rate is
# re-solved so the uncertainty model keeps the measured mean.
test_that("a measured multiplicity drives the count", {
  gt <- .gt_run(grid_molecules(8, 100, jitter = 15), 1, m_tid_override = 4)
  expect_identical(gt$meta$m_tid_source, "measured")
  expect_equal(gt$meta$tids_per_mol, 4)
  expect_equal(gt$summary$n_stoich, gt$summary$n_tid_total / 4)
  expect_identical(gt$summary$count_verdict, "over-split") # data carry ~2 per molecule
})

test_that("empty or unusable input returns empty frames, not NULL", {
  e <- estimate_gt_from_labeling(data.frame(), .gt_preset)
  expect_equal(nrow(e$points), 0L)
  expect_true(is.data.frame(e$summary))
  df <- data.frame(tid = 1, x = 0, y = 0)
  expect_equal(nrow(estimate_gt_from_labeling(df, .gt_preset, level = "primary")$points), 0L)
})
