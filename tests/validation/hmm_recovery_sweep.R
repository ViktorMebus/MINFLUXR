## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/hmm_recovery_sweep.R  -  where the continuous-time HMM can be trusted
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Recovery sweep for est_hmm.R against ground truth from sim_ct_state_tracks(),
##   which keeps what the model simplifies (switches inside a step, motion blur,
##   dark cycles, per-localization σ). Maps bias and CI coverage of D and dwell
##   times over (dwell time / cycle) × (D contrast), and reports the old
##   decoded-path rate estimate next to the MLE.
##
## RUN   Rscript tests/validation/hmm_recovery_sweep.R [n_reps]
##       (~15-20 min for 3 reps). Writes results/hmm_recovery_sweep.csv and a summary.
## ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args)) as.integer(args[1]) else 3L
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
suppressMessages({ source(file.path(root, "helpers_shared.R")); source(file.path(root, "R/est_hmm.R")) })

# Acquisition matched to the tracking test files (docs/simulation_module.md §9)
t_c <- 0.192; sw <- 0.011; R <- 1 / 6.2; gap_p <- c(0.72, 0.20, 0.05, 0.03)
set.seed(99); eco <- pmax(20, rlnorm(1e5, log(36), 0.6)); sig_pool <- 100.8 / (4 * sqrt(eco))
D_fast <- 300                                              # nm²/ms = 0.3 µm²/s
n_tracks <- 200L; n_locs <- 300L                           # 60k localizations per dataset

grid <- expand.grid(contrast = c(0.03, 0.1, 0.3), dwell_cyc = c(1, 3, 10, 30, 100, 300), rep = seq_len(n_reps))
rows <- vector("list", nrow(grid))
for (g in seq_len(nrow(grid))) {
  ct <- grid$contrast[g]; dw <- grid$dwell_cyc[g] * t_c
  D  <- c(ct * D_fast, D_fast)
  Q  <- matrix(c(0, 1 / dw, 1 / dw, 0), 2); diag(Q) <- -rowSums(Q)
  df <- sim_ct_state_tracks(n_tracks, n_locs, D, Q, t_c = t_c, sw = sw, gap_p = gap_p,
                            sigma = sig_pool, R_blur = R, seed = 1000 * g)
  pr <- hmm_ct_prepare(df, sigma = 0)
  t0 <- Sys.time()
  f  <- tryCatch(hmm_ct_fit(pr, K = 2, R_abs = R * t_c, t_c = t_c), error = function(e) NULL)
  el <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(f)) next
  vs <- hmm_ct_viterbi(pr, f)
  acc <- mean(vs == pr$df$state[pr$steps$row])
  old <- hmm_path_rates(pr, vs, 2, t_c)
  cov_dw <- if (!is.null(f$ci)) as.numeric(f$ci$dwell_ms[1, ] <= dw & f$ci$dwell_ms[2, ] >= dw) else c(NA, NA)
  cov_D  <- if (!is.null(f$ci)) as.numeric(f$ci$D[1, ] <= D & f$ci$D[2, ] >= D) else c(NA, NA)
  rows[[g]] <- data.frame(
    contrast = ct, dwell_cyc = grid$dwell_cyc[g], dwell_ms = dw, rep = grid$rep[g],
    D1_rel = f$D[1] / D[1] - 1, D2_rel = f$D[2] / D[2] - 1,
    dwell1_rel = f$dwell_ms[1] / dw - 1, dwell2_rel = f$dwell_ms[2] / dw - 1,
    covD1 = cov_D[1], covD2 = cov_D[2], cov_dwell1 = cov_dw[1], cov_dwell2 = cov_dw[2],
    ci_dwell1_width_rel = if (!is.null(f$ci)) diff(f$ci$dwell_ms[, 1]) / dw else NA,
    path_dwell1_rel = old$dwell_ms[1] / dw - 1, path_dwell2_rel = old$dwell_ms[2] / dw - 1,
    state_acc = acc, exp_switches = sum(f$expected_switches, na.rm = TRUE),
    converged = f$converged, fit_s = el)
  cat(sprintf("[%d/%d] contrast %.2f dwell %4.0f cyc rep %d: D %+.2f/%+.2f dwell %+.2f/%+.2f (path %+.2f/%+.2f) acc %.2f  %.0fs\n",
              g, nrow(grid), ct, grid$dwell_cyc[g], grid$rep[g],
              rows[[g]]$D1_rel, rows[[g]]$D2_rel, rows[[g]]$dwell1_rel, rows[[g]]$dwell2_rel,
              rows[[g]]$path_dwell1_rel, rows[[g]]$path_dwell2_rel, acc, el))
}
res <- do.call(rbind, rows)
out <- file.path(root, "tests", "validation", "results", "hmm_recovery_sweep.csv")
utils::write.csv(res, out, row.names = FALSE)
agg <- stats::aggregate(cbind(D1_rel, D2_rel, dwell1_rel, dwell2_rel, path_dwell1_rel, path_dwell2_rel,
                              cov_dwell1, cov_dwell2, covD1, covD2, ci_dwell1_width_rel, state_acc, exp_switches) ~ contrast + dwell_cyc,
                        data = res, FUN = mean)
print(format(agg, digits = 2), row.names = FALSE)
utils::write.csv(agg, sub(".csv", "_summary.csv", out, fixed = TRUE), row.names = FALSE)
