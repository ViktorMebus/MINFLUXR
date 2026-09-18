## Diagnostic: is the fast switching on GP67H-014 produced by dark-cycle steps?
## Fits K = 2 using single-cycle steps only (tracks split at every dark cycle), after
## tabulating mean r2 by cycles per step (the compression that would fake a slow state).
## RUN  Rscript tests/validation/hmm_real_gp67h_integrated.R > tests/validation/results/hmm_real_gp67h_integrated.log
setwd(normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "..")))
source("helpers_shared.R"); source("est_hmm.R")
f <- list.files("Test data", "GP67H-track-014_filtered.*csv$", full.names = TRUE)[1]
d <- utils::read.csv(f); d$tim <- d$tim * 1000; d$tid <- as.character(d$tid)
t_c <- 0.1914; sw <- 0.0107; R_abs <- (1/6.2) * t_c
d$sig <- 100.8 / (4 * sqrt(pmax(d$eco, 1))) * 1.016
# 1) mean r2 per cycle count (the known compression), from all steps
pr_all <- hmm_ct_prepare(d, sigma = 4.2, min_pts = 10)
st <- pr_all$steps; st$dt <- pr_all$DT[cbind(st$chunk, st$col)]; st$r2 <- pr_all$R2[cbind(st$chunk, st$col)]
st$k <- pmax(1, round((st$dt - sw) / t_c))
tab <- aggregate(r2 ~ pmin(k, 5), data = st, FUN = mean); names(tab) <- c("k", "mean_r2")
tab$n <- as.vector(table(pmin(st$k, 5))); tab$free_pred <- tab$mean_r2[1] * (tab$k * t_c + sw) / (t_c + sw)
cat("mean r2 (nm2) by cycles per step:\n"); print(tab, row.names = FALSE)
# 2) single-cycle steps only
max_dt <- 1.5 * t_c + sw
pr1 <- hmm_ct_prepare(d, sigma = 4.2, min_pts = 10, max_dt = max_dt)
cat(sprintf("single-cycle steps: %d of %d, chunks %d\n", pr1$n_steps, pr_all$n_steps, nrow(pr1$R2)))
nll <- function(lD) -sum(.hmm_emissions(pr1, exp(lD), R_abs, log = TRUE)[!is.na(pr1$R2)])
o <- optimize(nll, c(log(1), log(5000)))
cat(sprintf("K=1: D %.3g um2/s  BIC %.0f\n", exp(o$minimum)/1e3, 2 * o$objective + log(pr1$n_steps)))
t0 <- Sys.time(); f2 <- hmm_ct_fit(pr1, K = 2, R_abs = R_abs, t_c = t_c)
cat(sprintf("K=2 single-cycle (%s): D %s um2/s | occupancy %s | dwell %s ms (= %s cycles) | BIC %.0f\n",
    format(Sys.time() - t0), paste(signif(f2$D/1e3, 3), collapse = " / "), paste(round(100 * f2$pi), collapse = "/"),
    paste(signif(f2$dwell_ms, 3), collapse = " / "), paste(signif(f2$dwell_ms / t_c, 3), collapse = " / "), f2$BIC))
if (!is.null(f2$ci)) cat("  dwell CI:", paste(signif(f2$ci$dwell_ms, 3), collapse = " "), "\n")
