setwd("C:/Users/viktorhm/repos/MINFLUXR")
suppressMessages(suppressWarnings({ source("packages.R"); source("helpers_shared.R"); source("helpers_smlm.R") }))
src <- readLines("tests/validation/spatial_envelope_coverage.R", encoding = "UTF-8")
a <- grep("^in_L <- function", src); b <- grep("^r_seq <- seq", src) - 1
eval(parse(text = paste(src[a:b], collapse = "\n")))
# true L-shape mask on the same 100 nm grid, as an estimate_support-like object
true_support <- function(win) {
  h <- 100; o <- c(0, 0); nx <- 80; ny <- 60
  xc <- (seq_len(nx) - 0.5) * h; yc <- (seq_len(ny) - 0.5) * h
  mk <- outer(xc, yc, function(x, y) in_L(x, y))
  s <- list(mask = mk, origin = o, h = h, support = sum(mk) * h^2, d = 2L, mode = "grid",
            frac_of_bbox = mean(mk), win = win)
  attach_gamma_cache(s, r_max = sqrt(sum((win[2, ] - win[1, ])^2)), use_3d = FALSE)
}
set.seed(7); r_seq <- seq(25, 1000, by = 25); i200 <- which(r_seq == 200)
for (n in c(500)) {
  res <- list()
  for (rep in 1:15) {
    a <- csr_L(n); w <- window_for(a)
    ts <- true_support(w$win)
    cn <- bridge_close(support_closing_sweep(transform(a, z = 0), use_3d = FALSE, h_nm = 100))
    for (kind in c("estimated", "true")) {
      sup <- if (kind == "estimated") w$sup else ts
      obs <- compute_ripley_L(a$x, a$y, r_seq = r_seq, window_ranges = w$win, support = sup)$Lr_minus_r
      env <- compute_ripley_envelope(n = n, window_ranges = w$win, r_seq = r_seq, nsim = 39, estimator = "L", support = sup)
      out <- obs < env$lo | obs > env$hi
      res[[length(res) + 1]] <- data.frame(kind = kind, rep = rep, close = cn, frac_bbox = sup$frac_of_bbox,
                                           area_ratio = sup$support / ts$support,
                                           below200 = obs[i200] < env$lo[i200], above200 = obs[i200] > env$hi[i200],
                                           anywhere = any(out, na.rm = TRUE), obs200 = obs[i200], null200 = env$null_mean[i200])
    }
  }
  r <- do.call(rbind, res)
  cat(sprintf("\n=== n = %d (true area %.1f um2) ===\n", n, true_support(w$win)$support / 1e6))
  print(aggregate(cbind(close, frac_bbox, area_ratio, below200, above200, anywhere, obs200, null200) ~ kind, data = r, FUN = mean), digits = 3)
}
