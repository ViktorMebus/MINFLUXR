## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/hotspot_csr_null.R  -  Local density & hotspots under CSR
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The hotspot panel reports statistics rather than a yes/no call, so its
##   "false-positive" behaviour is the baseline each statistic takes on a pattern
##   with NO clustering: complete spatial randomness inside a cell footprint.
##   For a between-cell comparison to be valid that baseline must not move with
##   labelling density or footprint shape. Measured here:
##     cv_ratio      (CV of local density / CSR expectation; should be ≈ 1)
##     vmr           (quadrat variance-to-mean; ≈ 1 under CSR)
##     frac_in_foci  (should be ≈ the calibration target, 5%)
##     max_focus_n   (largest focus under CSR: the null for the headline metric)
##     frac_above    (units with enhancement > 5)
##     hit_floor     (eps calibration failed to converge)
##   plus a positive control (20% of units in 6 Gaussian foci, σ 40 nm).
##
## Pipeline as in the app's hotspot panel: footprint = estimate_support(grid,
## h 100 nm) at the bridging-jump gap-fill radius of support_closing_sweep();
## eps = calibrate_dbscan_eps(target 5%, minPts 5, nsim 15); hotspot_summary_row(k 5,
## quadrat 250 nm, enh_thresh 5).
##
## RUN   Rscript tests/validation/hotspot_csr_null.R [n_reps]   (~5-15 min)
## ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args)) as.integer(args[1]) else 20L
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings({ source("packages.R"); source("helpers_shared.R"); source("helpers_smlm.R") }))

# Footprints (nm): all inside a 10 x 6 µm frame
in_shape <- function(shape, x, y) switch(shape,
  frame   = rep(TRUE, length(x)),                                         # fills its box
  rod     = ((x - 5000) / 4000)^2 + ((y - 3000) / 1000)^2 <= 1,           # E. coli-like
  lshape  = !(x > 5000 & y > 3000))                                       # concave
draw_csr <- function(shape, n) {
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    x <- stats::runif(4 * n, 0, 10000); y <- stats::runif(4 * n, 0, 6000)
    k <- in_shape(shape, x, y); out <- rbind(out, data.frame(x = x[k], y = y[k]))
  }
  out[seq_len(n), ]
}
draw_clustered <- function(shape, n, frac = 0.2, n_foci = 6, sd_nm = 40) {
  bg <- draw_csr(shape, round(n * (1 - frac)))
  cen <- draw_csr(shape, n_foci)
  m <- round(n * frac); id <- sample.int(n_foci, m, TRUE)
  fo <- data.frame(x = cen$x[id] + stats::rnorm(m, 0, sd_nm), y = cen$y[id] + stats::rnorm(m, 0, sd_nm))
  rbind(bg, fo[in_shape(shape, fo$x, fo$y), ])
}
# Gap-fill radius by the bridging jump (same rule as spatial_resolve_close in server_smlm.R)
bridge_close <- function(sw) {
  if (is.null(sw) || nrow(sw) < 3L) return(NA_real_)
  s  <- sw[order(sw$close_nm), , drop = FALSE]
  s  <- s[is.finite(s$frac_of_bbox) & s$frac_of_bbox <= 0.97, , drop = FALSE]
  if (nrow(s) < 3L) return(NA_real_)
  j <- which.max(s$rel_step[-1]) + 1L
  k <- which(seq_len(nrow(s)) > j & is.finite(s$rel_step) & s$rel_step < 0.10)
  if (length(k)) s$close_nm[k[1]] else s$close_nm[nrow(s)]
}
run_one <- function(pts) {
  sw  <- support_closing_sweep(pts, use_3d = FALSE, h_nm = 100)
  cn  <- bridge_close(sw); if (!is.finite(cn)) cn <- 200
  sup <- estimate_support(pts, use_3d = FALSE, mode = "grid", h_nm = 100, close_nm = cn, return_mask = TRUE)
  ec  <- calibrate_dbscan_eps(pts, support = sup, minPts = 5L, target_frac = 0.05, nsim = 15L)
  r   <- hotspot_summary_row(pts, support = sup, k = 5L, quadrat_nm = 250, enh_thresh = 5,
                             eps_cal = ec, minPts = 5L)
  r$hit_floor <- isTRUE(ec$hit_floor); r$close_nm <- cn; r$frac_of_bbox <- sup$frac_of_bbox
  r
}

set.seed(2026)
grid <- expand.grid(shape = c("frame", "rod", "lshape"), n = c(200L, 1000L, 5000L),
                    pattern = c("csr", "clustered"), rep = seq_len(n_reps), stringsAsFactors = FALSE)
grid <- grid[!(grid$pattern == "clustered" & grid$rep > max(3L, n_reps %/% 4L)), ]   # fewer positive controls
rows <- vector("list", nrow(grid))
t0 <- Sys.time()
for (g in seq_len(nrow(grid))) {
  pts <- if (grid$pattern[g] == "csr") draw_csr(grid$shape[g], grid$n[g]) else draw_clustered(grid$shape[g], grid$n[g])
  pts$z <- 0
  r <- tryCatch(run_one(pts), error = function(e) { message("fail: ", conditionMessage(e)); NULL })
  if (!is.null(r)) rows[[g]] <- cbind(grid[g, ], r[, c("n_units", "cv_ratio", "vmr", "frac_in_foci", "max_focus_n",
                                                      "frac_above", "enh_q99", "eps_nm", "hit_floor", "frac_of_bbox",
                                                      "max_focus_ratio")])
  if (g %% 20 == 0) cat(sprintf("[%d/%d] %.1f min\n", g, nrow(grid), as.numeric(Sys.time() - t0, units = "mins")))
}
res <- do.call(rbind, rows)
dir.create(file.path(root, "tests", "validation", "results"), showWarnings = FALSE)
utils::write.csv(res, file.path(root, "tests", "validation", "results", "hotspot_csr_null.csv"), row.names = FALSE)
agg <- stats::aggregate(cbind(cv_ratio, vmr, frac_in_foci, max_focus_n, max_focus_ratio, frac_above, eps_nm, hit_floor, frac_of_bbox) ~
                          pattern + shape + n, data = res, FUN = mean)
agg <- agg[order(agg$pattern, agg$shape, agg$n), ]
print(format(agg, digits = 3), row.names = FALSE)
# 95th percentile of the CSR null for the headline metrics, per cell
q95 <- stats::aggregate(cbind(cv_ratio, vmr, max_focus_n, max_focus_ratio, frac_in_foci) ~ shape + n, data = res[res$pattern == "csr", ],
                        FUN = function(v) stats::quantile(v, 0.95))
cat("\nCSR 95th percentiles:\n"); print(format(q95, digits = 3), row.names = FALSE)
