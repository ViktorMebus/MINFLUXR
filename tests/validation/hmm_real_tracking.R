## Real-data check of est_hmm.R on the tracking test file (K = 1-3, bootstrap over tracks for K = 2).
## RUN  Rscript tests/validation/hmm_real_tracking.R > tests/validation/results/hmm_real_tracking.log
# The repository root, two folders above this script, so it runs from any machine.
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
source("helpers_shared.R"); source("R/est_hmm.R")
# The test data is local and not part of the repository. Point MINFLUXR_TRACK_PATTERN at a
# filename pattern under "Test data/" to run this on a different tracking acquisition.
f <- list.files("Test data", Sys.getenv("MINFLUXR_TRACK_PATTERN", "track.*_filtered.*csv$"), full.names = TRUE)[1]
stopifnot(!is.na(f))
cat("file:", basename(f), "
")
d <- utils::read.csv(f); d$tim <- d$tim * 1000; d$tid <- as.character(d$tid)
t_c <- 0.1914; R_abs <- (1/6.2) * t_c
d$sig <- 100.8 / (4 * sqrt(pmax(d$eco, 1))) * 1.016
cat(nrow(d), "locs,", length(unique(d$tid)), "tracks; dt quantiles:", round(quantile(diff(d$tim)[d$tid[-1]==d$tid[-nrow(d)]], c(.5,.9,.99,1)),3), "\n")
pr <- hmm_ct_prepare(d, sigma = 4.2, min_pts = 10)
cat("steps", pr$n_steps, "unique dt", length(pr$dt_u), "\n")
fits <- list()
for (K in 1:3) {
  if (K == 1) { # single-state likelihood for BIC reference: D MLE by 1-D optimisation
    nll <- function(lD) { E <- .hmm_emissions(pr, exp(lD), R_abs, log = TRUE); -sum(E[!is.na(pr$R2)]) }
    o <- optimize(nll, c(log(1), log(5000))); cat(sprintf("K=1: D %.3g um2/s  BIC %.0f\n", exp(o$minimum)/1e3, 2*o$objective + log(pr$n_steps)))
    next }
  t0 <- Sys.time(); fk <- hmm_ct_fit(pr, K = K, R_abs = R_abs, t_c = t_c); el <- format(Sys.time() - t0)
  fits[[K]] <- fk
  cat(sprintf("K=%d (%s): D %s um2/s | occupancy %s | dwell %s ms | BIC %.0f | conv %s\n", K, el,
      paste(signif(fk$D/1e3, 3), collapse=" / "), paste(round(100*fk$pi), collapse="/"),
      paste(signif(fk$dwell_ms, 3), collapse=" / "), fk$BIC, fk$converged))
  if (!is.null(fk$ci)) { cat("  D CI:", paste(signif(fk$ci$D/1e3,3), collapse=" "), "\n  dwell CI:", paste(signif(fk$ci$dwell_ms,3), collapse=" "), "\n") }
  vs <- hmm_ct_viterbi(pr, fk); old <- hmm_path_rates(pr, vs, K, t_c)
  cat("  path-count dwell ms:", paste(signif(old$dwell_ms, 3), collapse=" / "), " expected switches:", round(sum(fk$expected_switches, na.rm=TRUE)), "\n")
}
t0 <- Sys.time()
bt <- hmm_ct_bootstrap(pr$df, fits[[2]], B = 20, sigma = 4.2, t_c = t_c)
cat("bootstrap (B=20,", format(Sys.time()-t0), "): D", paste(signif(bt$D/1e3,3), collapse=" "), "| dwell", paste(signif(bt$dwell_ms,3), collapse=" "), "\n")
