## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/spatial_envelope_coverage.R  -  CSR envelopes: how often random data "leaves"
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The Ripley L, g(r), L12 and g12 panels shade a 2.5-97.5% Monte-Carlo envelope
##   at every r and are read as "significant where the observed curve leaves the
##   band". For data that ARE random inside the cell (one channel: CSR; two
##   channels: independent CSR patterns), this measures
##     pointwise  - fraction of replicates outside the POINTWISE band at one fixed r
##                  (should be ≈ 5% if the envelope is calibrated)
##     anywhere   - fraction of replicates outside the pointwise band at ANY r
##                  (the false-positive rate of the reading the UI invites)
##     global     - fraction outside the GLOBAL band anywhere, and the rate of
##                  global p < 0.05 (both should be ≈ 5%)
##   with the footprint window and edge correction ON, in a concave (L-shaped)
##   cell that does not fill its frame, i.e. the case the window exists for.
##
## RUN   Rscript tests/validation/spatial_envelope_coverage.R [n_reps] [nsim]
##       defaults 100 reps, 99 simulations per envelope (~20-40 min)
## ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100L
nsim   <- if (length(args) >= 2) as.integer(args[2]) else 99L
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings({ source("packages.R"); source("helpers_shared.R"); source("helpers_smlm.R") }))

in_L <- function(x, y) !(x > 4000 & y > 3000)                # 8 x 6 µm frame, one quadrant empty
csr_L <- function(n) {
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    x <- runif(3 * n, 0, 8000); y <- runif(3 * n, 0, 6000)
    k <- in_L(x, y); out <- rbind(out, data.frame(x = x[k], y = y[k]))
  }
  out[seq_len(n), ]
}
bridge_close <- function(sw) {
  if (is.null(sw) || nrow(sw) < 3L) return(NA_real_)
  s <- sw[order(sw$close_nm), , drop = FALSE]
  s <- s[is.finite(s$frac_of_bbox) & s$frac_of_bbox <= 0.97, , drop = FALSE]
  if (nrow(s) < 3L) return(NA_real_)
  j <- which.max(s$rel_step[-1]) + 1L
  k <- which(seq_len(nrow(s)) > j & is.finite(s$rel_step) & s$rel_step < 0.10)
  if (length(k)) s$close_nm[k[1]] else s$close_nm[nrow(s)]
}
window_for <- function(pts) {
  p  <- pts; p$z <- 0
  cn <- bridge_close(support_closing_sweep(p, use_3d = FALSE, h_nm = 100)); if (!is.finite(cn)) cn <- 400
  sup <- estimate_support(p, use_3d = FALSE, mode = "grid", h_nm = 100, close_nm = cn, return_mask = TRUE)
  sup <- debias_support_window(sup, nrow(p))            # as spatial_window_of() does
  win <- apply(as.matrix(pts[, c("x", "y")]), 2, range)
  sup <- attach_gamma_cache(sup, r_max = sqrt(sum((win[2, ] - win[1, ])^2)), use_3d = FALSE)
  list(sup = sup, win = win)
}

r_seq <- seq(25, 1000, by = 25)
r_fix <- 200                                                   # pointwise check at this r
i_fix <- which(r_seq == r_fix)
set.seed(42)
out <- list()
t0 <- Sys.time()
for (est in c("L", "g", "L12", "g12", "L12lab", "g12lab")) {
  rows <- vector("list", n_reps)
  for (rep in seq_len(n_reps)) {
    if (est %in% c("L", "g")) {
      a <- csr_L(500); w <- window_for(a)
      obs <- if (est == "L") compute_ripley_L(a$x, a$y, r_seq = r_seq, window_ranges = w$win, support = w$sup)$Lr_minus_r
             else compute_pcf(a$x, a$y, r_seq = r_seq, window_ranges = w$win, support = w$sup)$g
      env <- compute_ripley_envelope(n = 500, window_ranges = w$win, r_seq = r_seq, nsim = nsim,
                                     estimator = est, support = w$sup)
    } else {
      a <- csr_L(400); b <- csr_L(300); w <- window_for(rbind(a, b))
      is_L <- est %in% c("L12", "L12lab")
      obs <- if (is_L) compute_ripley_L12(a$x, a$y, NULL, b$x, b$y, NULL, r_seq = r_seq,
                                          window_ranges = w$win, support = w$sup)$Lr_minus_r
             else compute_pcf12(a$x, a$y, NULL, b$x, b$y, NULL, r_seq = r_seq,
                                window_ranges = w$win, support = w$sup)$g12
      lab <- grepl("lab$", est)
      env <- compute_ripley_envelope(n = 400, n2 = 300, window_ranges = w$win, r_seq = r_seq, nsim = nsim,
                                     estimator = if (is_L) "L12" else "pcf", support = w$sup,
                                     null_mode = if (lab) "labels" else "csr",
                                     pts1 = if (lab) a else NULL, pts2 = if (lab) b else NULL)
    }
    k <- is.finite(obs) & is.finite(env$lo) & is.finite(env$hi)
    outside <- (obs < env$lo | obs > env$hi) & k
    gk <- is.finite(env$g_lo) & is.finite(env$g_hi) & is.finite(obs)
    g_out <- (obs < env$g_lo | obs > env$g_hi) & gk
    gp <- envelope_global_p(env, obs)
    rows[[rep]] <- data.frame(estimator = est, rep = rep,
                              pointwise = outside[i_fix], anywhere = any(outside),
                              above_anywhere = any(obs > env$hi & k), n_r = sum(k),
                              global_out = any(g_out), global_p = gp$p %||% NA_real_,
                              global_sig = isTRUE((gp$p %||% 1) < 0.05),
                              global_above = any(obs > env$g_hi & gk))
  }
  out[[est]] <- do.call(rbind, rows)
  cat(sprintf("%-4s done (%.1f min)\n", est, as.numeric(Sys.time() - t0, units = "mins")))
}
res <- do.call(rbind, out)
dir.create(file.path(root, "tests", "validation", "results"), showWarnings = FALSE)
utils::write.csv(res, file.path(root, "tests", "validation", "results", "spatial_envelope_coverage.csv"), row.names = FALSE)
cat(sprintf("\n%d reps per estimator, %d sims per envelope, %d radii (%g-%g nm), pointwise at r = %g nm\n",
            n_reps, nsim, length(r_seq), min(r_seq), max(r_seq), r_fix))
agg <- aggregate(cbind(pointwise, anywhere, above_anywhere, global_out, global_sig, global_above) ~ estimator,
                 data = res, FUN = mean)
print(format(agg, digits = 3), row.names = FALSE)
cat("\nglobal p quartiles (uniform under the null):\n")
print(aggregate(global_p ~ estimator, data = res, FUN = function(v) round(stats::quantile(v, c(.25, .5, .75)), 3)))
