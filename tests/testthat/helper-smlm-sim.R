# Ground-truth simulators shared by the SMLM test files (testthat sources every
# helper-*.R before the tests). Units are nm throughout, as in the SMLM tab.
.smlm_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages(suppressWarnings(source(file.path(.smlm_root, "packages.R"))))

# ── sim_smlm_traces() ────────────────────────────────────────────────────────
# GOAL: a localization table with KNOWN molecules, so de-blinking and counting can
# be scored against the truth.
#   each molecule  -> 1 + Poisson(blink_mean) traces (re-blinks of one dye)
#   each trace     -> its own small offset from the molecule (between-blink wobble)
#   each trace     -> n_loc localizations scattered with per-axis sd sigma
# ARGS  mol    : matrix, one row per molecule, 2 or 3 columns (x, y[, z])
#       sigma  : per-axis localization sd; a length-3 vector gives anisotropic data
#       dcr    : optional per-molecule DCR mean; each trace's dcr is drawn around it
# RETURNS data.frame(tid, x, y, z, mol[, dcr]) - tid unique per trace
sim_smlm_traces <- function(mol, blink_mean = 1, n_loc = 15L, sigma = 5,
                            wobble = 1, dcr = NULL, dcr_sd = 0.03) {
  mol <- as.matrix(mol)
  if (ncol(mol) == 2L) mol <- cbind(mol, 0)
  sig <- rep_len(sigma, 3)
  if (all(mol[, 3] == 0)) sig[3] <- 0               # 2D data carries z = 0
  n_tr <- 1L + stats::rpois(nrow(mol), blink_mean)
  m_of <- rep(seq_len(nrow(mol)), n_tr)             # molecule of each trace
  rows <- lapply(seq_along(m_of), function(t) {
    m   <- m_of[t]
    ctr <- mol[m, ] + c(stats::rnorm(2, 0, wobble), if (sig[3] > 0) stats::rnorm(1, 0, wobble) else 0)
    out <- data.frame(tid = t,
                      x = stats::rnorm(n_loc, ctr[1], sig[1]),
                      y = stats::rnorm(n_loc, ctr[2], sig[2]),
                      z = if (sig[3] > 0) stats::rnorm(n_loc, ctr[3], sig[3]) else 0,
                      mol = m)
    if (!is.null(dcr)) out$dcr <- stats::rnorm(n_loc, stats::rnorm(1, dcr[m], dcr_sd), 0.01)
    out
  })
  do.call(rbind, rows)
}

# GOAL: molecule positions on a jittered grid, so the nearest-neighbour spacing is
# known and well above the localization noise.
grid_molecules <- function(n_side, spacing, jitter = 0) {
  g <- expand.grid(x = seq_len(n_side) * spacing, y = seq_len(n_side) * spacing)
  as.matrix(g) + stats::runif(2 * nrow(g), -jitter, jitter)
}

# GOAL: a Thomas cluster process in a square window - n_par parents, each with
# Poisson(mu) offspring scattered with sd `s`. Returns matrix(x, y).
sim_thomas <- function(n_par, mu, s, W) {
  par <- cbind(stats::runif(n_par, 0, W), stats::runif(n_par, 0, W))
  k   <- stats::rpois(n_par, mu)
  P   <- par[rep(seq_len(n_par), k), , drop = FALSE]
  P   <- P + matrix(stats::rnorm(2 * nrow(P), 0, s), ncol = 2)
  P[P[, 1] >= 0 & P[, 1] <= W & P[, 2] >= 0 & P[, 2] <= W, , drop = FALSE]
}

# GOAL: score a grouping against the truth. Returns the fraction of truth groups
# that came out as exactly one label (no split) and the fraction of labels that hold
# exactly one truth group (no merge), plus the label count.
score_grouping <- function(truth, label) {
  tab <- table(truth, label)
  c(not_split  = mean(rowSums(tab > 0) == 1),
    not_merged = mean(colSums(tab > 0) == 1),
    n_label    = length(unique(label)))
}

# Settings for apply_clustering_pipeline() and friends: the app's defaults
# (smlm_settings, R/smlm_settings.R) with de-blinking switched on, since almost
# every test needs it. Any field can be overridden: smlm_input(enable_cluster2 = TRUE).
smlm_input <- function(...) {
  args <- list(...)
  if (is.null(args$enable_clustering)) args$enable_clustering <- TRUE
  do.call(smlm_settings, args)
}
