## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/hotspot_realdata_recheck.R  -  hotspot metrics on the two-condition set
## ─────────────────────────────────────────────────────────────────────────────
## Re-runs the Local Density & Hotspots statistics on the 10 two-condition acquisitions
## (4 cond01 without IDR, 6 cond02 with IDR) after the 2026-09-17 fixes (whole-
## quadrat VMR, extended eps grid, max_focus_ratio), to check the published
## separation still holds. A headless approximation of the app pipeline: trace
## length >= 5, per-axis spread < 10 nm, trace-ID centroids, 2D, footprint at
## 25 nm pitch with the bridging-jump gap fill shared by the whole set, eps
## calibrated per cell (5%, minPts 5), k 5, quadrat 250 nm.
## RUN  Rscript tests/validation/hotspot_realdata_recheck.R   (needs Test data/)
## ─────────────────────────────────────────────────────────────────────────────
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings({ source("packages.R"); source("helpers_shared.R"); source("helpers_smlm.R") }))
# The acquisitions are local and not part of the repository: point MINFLUXR_SMLM_TESTDIR
# at the folder holding the two-condition .npy files.
files <- list.files(Sys.getenv("MINFLUXR_SMLM_TESTDIR", "Test data"), "[.]npy$",
                    full.names = TRUE, recursive = TRUE)
cond  <- ifelse(grepl("cond01", files), "cond01", "cond02")

units <- lapply(files, function(f) {
  df <- read_minflux_npy_as_df(f)
  df$tid <- as.character(df$tid)
  len <- table(df$tid); df <- df[df$tid %in% names(len)[len >= 5], ]
  sds <- sd_by_tid(df, min_n = 5)
  ok  <- sds$tid[is.finite(sds$sd_x) & is.finite(sds$sd_y) & sds$sd_x < 10 & sds$sd_y < 10]
  df  <- df[df$tid %in% ok, ]
  pts <- unit_centroids(df, unit = "tid", use_3d = FALSE)
  pts$z <- 0
  pts
})
bridge_close <- function(sw) {
  if (is.null(sw) || nrow(sw) < 3L) return(NA_real_)
  s <- sw[order(sw$close_nm), , drop = FALSE]
  s <- s[is.finite(s$frac_of_bbox) & s$frac_of_bbox <= 0.97, , drop = FALSE]
  if (nrow(s) < 3L) return(NA_real_)
  j <- which.max(s$rel_step[-1]) + 1L
  k <- which(seq_len(nrow(s)) > j & is.finite(s$rel_step) & s$rel_step < 0.10)
  if (length(k)) s$close_nm[k[1]] else s$close_nm[nrow(s)]
}
h <- 25
per <- vapply(units, function(p) bridge_close(support_closing_sweep(p, use_3d = FALSE, h_nm = h)), numeric(1))
close_nm <- max(per, na.rm = TRUE)                              # shared radius, as the app does
cat("shared gap fill:", close_nm, "nm\n")
set.seed(1)
rows <- do.call(rbind, lapply(seq_along(units), function(i) {
  p   <- units[[i]]
  sup <- estimate_support(p, use_3d = FALSE, mode = "grid", h_nm = h, close_nm = close_nm, return_mask = TRUE)
  ec  <- calibrate_dbscan_eps(p, support = sup, minPts = 5L, target_frac = 0.05, nsim = 15L)
  r   <- hotspot_summary_row(p, support = sup, k = 5L, quadrat_nm = 250, enh_thresh = 5, eps_cal = ec,
                             minPts = 5L, dataset = basename(files[i]))
  di  <- compute_dispersion_index(p, quadrat_nm = 250, support = sup)
  cbind(cond = cond[i], r, edge_quadrats = di$edge_quadrats, n_quadrats = di$n_quadrats,
        grid_extended = isTRUE(ec$grid_extended))
}))
show <- c("cond", "n_units", "density", "cv_ratio", "enh_q90", "frac_above", "vmr", "n_quadrats", "edge_quadrats",
          "eps_nm", "max_focus_n", "max_focus_null95", "max_focus_ratio", "nn_med")
print(format(rows[order(rows$cond), show], digits = 3), row.names = FALSE)
cat("\nRanges by condition (DISJOINT = no overlap):\n")
for (m in c("density", "cv_ratio", "enh_q90", "frac_above", "vmr", "max_focus_n", "max_focus_ratio", "nn_med")) {
  a <- range(rows[[m]][rows$cond == "cond01"]); b <- range(rows[[m]][rows$cond == "cond02"])
  disj <- a[2] < b[1] || b[2] < a[1]
  p <- suppressWarnings(stats::wilcox.test(rows[[m]][rows$cond == "cond01"], rows[[m]][rows$cond == "cond02"])$p.value)
  cat(sprintf("  %-16s cond01 %8.3g - %-8.3g cond02 %8.3g - %-8.3g %s  rank-sum p %.3f\n", m, a[1], a[2], b[1], b[2],
              if (disj) "DISJOINT" else "overlap ", p))
}
utils::write.csv(rows, file.path(root, "tests", "validation", "results", "hotspot_realdata_recheck.csv"), row.names = FALSE)
