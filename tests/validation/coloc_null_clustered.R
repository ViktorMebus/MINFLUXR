## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/coloc_null_clustered.R  -  which two-channel null answers the co-localization question
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The L12 and g12 panels offer three nulls. spatial_envelope_coverage.R checked
##   them on independent CSR channels, the one case where all of them are exact. Real
##   channels are clustered, and there the nulls ask different questions. This runs
##   the global test (studentized max deviation, 5%) on clustered channels in a
##   concave cell that does not fill its frame, footprint window and edge correction
##   on, for three truths:
##     independent - A and B cluster at unrelated places     (co-localized: NO)
##     shared      - A and B drawn from the same clusters    (co-localized: YES)
##     partners    - every B sits within ~8 nm of an A       (co-localized: YES)
##   and reports the rejection rate and its direction per null.
##
## RUN   Rscript tests/validation/coloc_null_clustered.R [n_reps] [nsim]
##       defaults 30 reps, 99 simulations per envelope (~30-60 min)
## OUTPUT tests/validation/results/coloc_null_clustered.csv
## ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 30L
nsim   <- if (length(args) >= 2) as.integer(args[2]) else 99L
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings({ source("packages.R"); source("helpers_shared.R"); source("helpers_smlm.R") }))

in_L <- function(x, y) x >= 0 & x <= 8000 & y >= 0 & y <= 6000 & !(x > 4000 & y > 3000)  # one quadrant empty
# Parent centres uniform inside the cell; offspring Gaussian around them, kept inside.
parents_L <- function(m) {
  out <- NULL
  while (is.null(out) || nrow(out) < m) {
    x <- runif(3 * m, 0, 8000); y <- runif(3 * m, 0, 6000)
    k <- in_L(x, y); out <- rbind(out, cbind(x[k], y[k]))
  }
  out[seq_len(m), , drop = FALSE]
}
thomas_L <- function(P, n, s = 40) {
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    k <- sample(nrow(P), 2 * n, TRUE)
    x <- P[k, 1] + rnorm(2 * n, 0, s); y <- P[k, 2] + rnorm(2 * n, 0, s)
    keep <- in_L(x, y); out <- rbind(out, data.frame(x = x[keep], y = y[keep]))
  }
  out[seq_len(n), ]
}
scenarios <- list(
  independent = function() list(thomas_L(parents_L(40), 400), thomas_L(parents_L(40), 300)),
  shared      = function() { P <- parents_L(40); list(thomas_L(P, 400), thomas_L(P, 300)) },
  partners    = function() {
    a <- thomas_L(parents_L(40), 400); b <- a[sample(nrow(a), 300), ]
    b$x <- b$x + rnorm(300, 0, 8); b$y <- b$y + rnorm(300, 0, 8)
    b <- b[in_L(b$x, b$y), ]; list(a, b)
  })

# Same window construction as spatial_window_of() in the app (bridging-jump closing
# radius, grid footprint, debiased area, cached set covariance).
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
  sup <- debias_support_window(sup, nrow(p))
  win <- apply(as.matrix(pts[, c("x", "y")]), 2, range)
  sup <- attach_gamma_cache(sup, r_max = sqrt(sum((win[2, ] - win[1, ])^2)), use_3d = FALSE)
  list(sup = sup, win = win)
}

r_seq <- seq(25, 600, by = 25)
set.seed(7)
rows <- list(); t0 <- Sys.time()
for (sc in names(scenarios)) for (rep in seq_len(n_reps)) {
  ab <- scenarios[[sc]](); a <- ab[[1]]; b <- ab[[2]]; w <- window_for(rbind(a, b))
  for (est in c("L12", "pcf")) {
    obs <- if (est == "L12")
      compute_ripley_L12(a$x, a$y, NULL, b$x, b$y, NULL, r_seq = r_seq, window_ranges = w$win, support = w$sup)$Lr_minus_r
    else compute_pcf12(a$x, a$y, NULL, b$x, b$y, NULL, r_seq = r_seq, window_ranges = w$win, support = w$sup)$g12
    for (nl in c("shift", "labels", "csr")) {
      env <- compute_ripley_envelope(n = nrow(a), n2 = nrow(b), window_ranges = w$win, r_seq = r_seq,
                                     nsim = nsim, estimator = est, support = w$sup, null_mode = nl,
                                     pts1 = a, pts2 = b)
      gp <- envelope_global_p(env, obs)
      k  <- which.max(abs(obs - env$null_mean) / env$null_sd)
      rows[[length(rows) + 1L]] <- data.frame(scenario = sc, estimator = est, null = nl, rep = rep,
        p = gp$p, above = obs[k] > env$null_mean[k], r_at_max = env$r[k],
        fallbacks = attr(env, "shift_fallbacks"))
    }
  }
  if (rep %% 5 == 0) cat(sprintf("%s rep %d  (%.1f min)\n", sc, rep, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, rows)
dir.create("tests/validation/results", showWarnings = FALSE)
utils::write.csv(res, "tests/validation/results/coloc_null_clustered.csv", row.names = FALSE)
agg <- aggregate(cbind(reject = p < 0.05, called_colocalized = p < 0.05 & above,
                       called_exclusion = p < 0.05 & !above) ~ scenario + estimator + null, res, mean)
agg$fallback_rate <- aggregate(fallbacks ~ scenario + estimator + null, res, function(v) mean(v, na.rm = TRUE) / nsim,
                               na.action = na.pass)$fallbacks
print(agg[order(agg$estimator, agg$scenario, agg$null), ], row.names = FALSE, digits = 3)
