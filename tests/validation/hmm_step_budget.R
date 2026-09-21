## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/hmm_step_budget.R  -  what the HMM step budget costs
## ─────────────────────────────────────────────────────────────────────────────
## The HMM panel fits at most `Step budget for the fit` steps, drawing whole
## tracks at random when the data hold more. This measures what that costs, by
## fitting simulated tracks with KNOWN D and dwell times at several budgets.
## Fit time scales with the step count; the question is whether the estimate does.
## RUN  Rscript tests/validation/hmm_step_budget.R   (no Test data/ needed)
## OUT  tests/validation/results/hmm_step_budget.csv
## ─────────────────────────────────────────────────────────────────────────────
root <- normalizePath(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
suppressMessages(source(file.path(root, "R/est_hmm.R")))

t_c <- 0.192; sw <- 0.011; R_blur <- 1 / 6.2
D_true <- c(50, 500)                      # nm²/ms
dwell_cycles <- c(20, 10)                 # resolvable: the panel's range starts near 10
dw <- dwell_cycles * t_c
Q <- matrix(c(-1 / dw[1], 1 / dw[2], 1 / dw[1], -1 / dw[2]), 2, byrow = TRUE)

df <- sim_ct_state_tracks(400, 500, D = D_true, Q = Q, t_c = t_c, sw = sw,
                          sigma = 5, R_blur = R_blur, seed = 11)

# The server's rule: whole tracks, drawn at random, until the budget is met.
subsample <- function(d, cap) {
  if (cap <= 0) return(d)
  per <- table(d$tid)
  if (sum(per - 1L) <= cap) return(d)
  set.seed(1L)
  sh <- sample(names(per))
  d[d$tid %in% sh[cumsum(per[sh] - 1L) <= cap], , drop = FALSE]
}

rows <- do.call(rbind, lapply(c(0, 150000, 50000, 20000), function(cap) {
  pr <- hmm_ct_prepare(subsample(df, cap), ndim = 2L, sigma = 5, min_pts = 10)
  el <- system.time(f <- hmm_ct_fit(pr, K = 2L, R_abs = R_blur * t_c, t_c = t_c))[["elapsed"]]
  data.frame(budget = cap, steps = pr$n_steps,
             D1 = f$D[1], D2 = f$D[2], D1_true = D_true[1], D2_true = D_true[2],
             dwell1_cyc = f$dwell_ms[1] / t_c, dwell2_cyc = f$dwell_ms[2] / t_c,
             dwell1_true = dwell_cycles[1], dwell2_true = dwell_cycles[2],
             elapsed_s = el)
}))
rows$D1_err <- rows$D1 / rows$D1_true - 1
rows$D2_err <- rows$D2 / rows$D2_true - 1
print(rows, row.names = FALSE, digits = 3)
out <- file.path(root, "tests/validation/results")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(rows, file.path(out, "hmm_step_budget.csv"), row.names = FALSE)
cat("\nwrote tests/validation/results/hmm_step_budget.csv\n")
