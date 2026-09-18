## ─────────────────────────────────────────────────────────────────────────────
## est_jdd.R  -  jump-distance distribution (JDD) mixture fit
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The estimator behind the Tracking tab's JDD panel, as pure functions. Moved
##   out of tracking_server() on 2026-09-17 (publication audit C1); jdd_bundle()
##   gathers inputs and calls jdd_collect_jumps() then jdd_fit().
##
## METHOD
##   Jump distances |Δr| (optionally √Δt-normalised, optionally only single-cycle
##   steps) are fitted by maximum likelihood with an N-component scaled-chi
##   mixture, the law of |Δr| for an isotropic Gaussian displacement in d dims:
##       p(r) = Σ w_k · C_d · r^(d−1)/s_k^d · exp(−r²/(2 s_k²)),  C_d = 1/(2^(d/2−1)Γ(d/2))
##   d = 2 Rayleigh (Schütz 1997), d = 3 Maxwell. s_k is the per-axis SD, so
##   D_k = (s_k² − noise)/2 in both 2D and 3D.
##
## UNITS: positions nm, tim ms; with √Δt normalisation s is nm/√ms and D nm²/ms.
## ─────────────────────────────────────────────────────────────────────────────


# ── jdd_collect_jumps() ───────────────────────────────────────────────────────
# GOAL: pool jump distances (and their Δt) over all tracks.
#
# df         : data.frame(tid, x, y, tim[, z]), finite coordinates, tim in ms
# min_pts    : minimum localizations per track
# step_mode  : "consec" (every step), "binmed" (median position per time bin of
#              `bin_width` ms), "subsamp" (every `subsamp_n`-th localization)
# min_dt, max_dt : keep steps with min_dt < Δt ≤ max_dt
# integ_only : keep only single-cycle steps, round(Δt/t_c) == 1
# t_c        : cycle time (ms), used by integ_only
# norm_dt    : divide each jump by √Δt
# use_z      : include z (only when z carries real axial data)
# progress   : optional function(fraction)
#
# RETURNS list(jumps, dt_all)
jdd_collect_jumps <- function(df, min_pts = 2L, step_mode = "consec", bin_width = 1,
                              subsamp_n = 4L, min_dt = 0, max_dt = Inf,
                              integ_only = FALSE, t_c = NA_real_, norm_dt = TRUE,
                              use_z = FALSE, progress = NULL) {
  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
  df <- df[order(df$tid, df$tim), ]
  sp <- split(df, df$tid)
  jumps  <- numeric(0)
  dt_all <- numeric(0)
  n_sp <- length(sp)
  for (ti in seq_along(sp)) {
    if (is.function(progress) && ti %% 50L == 0L) progress(50 / n_sp)
    tb <- sp[[ti]]
    if (nrow(tb) < min_pts) next

    # ---- Coarsen the track according to the chosen step mode ----------------
    if (identical(step_mode, "binmed")) {
      # Assign each localization to a time bin and take the median position
      bin_idx     <- floor(tb$tim / bin_width)
      unique_bins <- sort(unique(bin_idx))
      if (length(unique_bins) < 2) next
      bin_rows <- vector("list", length(unique_bins))
      for (bi_i in seq_along(unique_bins)) {
        idx <- which(bin_idx == unique_bins[bi_i])
        row <- data.frame(x = median(tb$x[idx]), y = median(tb$y[idx]),
                          tim = median(tb$tim[idx]))
        if (use_z && "z" %in% names(tb)) row$z <- median(tb$z[idx])
        bin_rows[[bi_i]] <- row
      }
      tb_use <- do.call(rbind, bin_rows)
    } else if (identical(step_mode, "subsamp")) {
      # Take every N-th localization
      idx_keep <- seq(1L, nrow(tb), by = subsamp_n)
      if (length(idx_keep) < 2) next
      tb_use <- tb[idx_keep, , drop = FALSE]
    } else {
      tb_use <- tb   # consecutive steps
    }
    if (nrow(tb_use) < 2) next

    dx <- diff(tb_use$x); dy <- diff(tb_use$y)
    r  <- if (use_z && "z" %in% names(tb_use))
             sqrt(dx^2 + dy^2 + diff(tb_use$z)^2)
          else sqrt(dx^2 + dy^2)
    dt <- diff(tb_use$tim)
    keep <- is.finite(dt) & dt > min_dt & dt <= max_dt & is.finite(r)
    if (isTRUE(integ_only)) keep <- keep & (round(dt / t_c) == 1L)
    if (!any(keep)) next
    if (norm_dt) {
      r    <- r / sqrt(dt)
      keep <- keep & is.finite(r)
    }
    if (!any(keep)) next
    dt_all <- c(dt_all, dt[keep])
    jumps  <- c(jumps,  r[keep])
  }
  list(jumps = jumps, dt_all = dt_all)
}


# ── jdd_fit() ─────────────────────────────────────────────────────────────────
# GOAL: histogram + N-component scaled-chi mixture MLE + D per component.
#
# jumps, dt_all : from jdd_collect_jumps()
# n_comp        : components (1-4)
# use_z         : TRUE = 3D Maxwell law, FALSE = 2D Rayleigh
# norm_dt       : whether the jumps were √Δt-normalised (sets the noise term)
# sigma_loc     : σ_loc (nm) for the noise correction; 0 = none
# R_abs         : motion-blur offset in ms (MSD_BLUR_R · t_c); 0 = no blur correction
#
# NOISE AND BLUR CORRECTION (fixed 2026-09-17). With √Δt normalisation the fitted
# per-axis variance of a component is
#     s² = 2D·(1 − 2·R_abs·⟨1/Δt⟩) + 2σ²·⟨1/Δt⟩
# so D = (s² − 2σ²⟨1/Δt⟩) / (2·(1 − 2·R_abs·⟨1/Δt⟩)). The previous version used
# D = (s² − 2σ²/median Δt)/2: no blur term and the median in place of ⟨1/Δt⟩.
# Simulated (GP67H timings, photon σ, blur R = 1/6.2, D 0.05 / 0.215 / 1 µm²/s):
# old −30% single-cycle and −26% to −54% pooled; corrected +1% in both modes.
#
# RETURNS list(mids, dens, sigs, ws, D_vals, dt0, inv_dt, r_noise, d_fit)
jdd_fit <- function(jumps, dt_all, n_comp = 2L, use_z = FALSE, norm_dt = TRUE,
                    sigma_loc = 0, R_abs = 0, progress = NULL) {
  # Histogram (truncate upper tail at 99.5th pct to avoid outlier stretching)
  q_hi <- quantile(jumps, 0.995, na.rm = TRUE)
  br   <- c(pretty(c(0, q_hi), n = 60), Inf)
  h    <- hist(jumps, breaks = br, plot = FALSE)
  mids <- h$mids[is.finite(h$mids)]
  dens <- h$density[seq_along(mids)]

  # Mixture NLL. Parameters: log(sigma_k), and for k>1: log(w_k / w_K) (softmax).
  d_jdd <- if (isTRUE(use_z)) 3L else 2L
  C_jdd <- 1 / (2^(d_jdd / 2 - 1) * gamma(d_jdd / 2))
  mix_nll <- function(par) {
    sigs <- exp(par[seq_len(n_comp)])
    if (n_comp > 1) {
      raw <- exp(par[n_comp + seq_len(n_comp - 1)])
      ws  <- c(raw, 1) / (1 + sum(raw))
    } else {
      ws <- 1
    }
    if (any(!is.finite(sigs)) || any(!is.finite(ws))) return(1e15)
    lik <- Reduce(`+`, lapply(seq_len(n_comp), function(k)
      ws[k] * C_jdd * jumps^(d_jdd - 1) / sigs[k]^d_jdd *
        exp(-jumps^2 / (2 * sigs[k]^2))))
    lik[lik <= 0] <- .Machine$double.eps
    -sum(log(lik))
  }

  sig0 <- seq(quantile(jumps, 0.05), quantile(jumps, 0.75), length.out = n_comp)
  sig0 <- pmax(sig0, 1e-3)
  par0 <- c(log(sig0), rep(0, max(0, n_comp - 1)))

  if (is.function(progress)) progress(0.4)
  # One component has a closed-form MLE (s² = mean(r²)/d); Nelder-Mead in one
  # dimension is unreliable, so it is only used for mixtures.
  fit <- if (n_comp == 1L) list(par = log(sqrt(mean(jumps^2) / d_jdd)), convergence = 0L)
         else tryCatch(
           optim(par0, mix_nll, method = "Nelder-Mead",
                 control = list(maxit = 8000, reltol = 1e-9)),
           error = function(e) NULL)

  if (!is.null(fit) && fit$convergence == 0) {
    sigs <- exp(fit$par[seq_len(n_comp)])
    if (n_comp > 1) {
      raw <- exp(fit$par[n_comp + seq_len(n_comp - 1)])
      ws  <- c(raw, 1) / (1 + sum(raw))
    } else {
      ws <- 1
    }
  } else {
    sigs <- sig0
    ws   <- rep(1 / n_comp, n_comp)
  }

  # Sort slow -> fast
  ord  <- order(sigs); sigs <- sigs[ord]; ws <- ws[ord]

  # Median Δt across all valid steps used in this fit (ms), reported for display
  dt0 <- median(dt_all, na.rm = TRUE)
  if (!is.finite(dt0) || dt0 <= 0) dt0 <- NA_real_
  # ⟨1/Δt⟩ is what the noise and blur terms average to after √Δt normalisation
  inv_dt <- mean(1 / dt_all[is.finite(dt_all) & dt_all > 0])
  if (!is.finite(inv_dt)) inv_dt <- NA_real_

  # D from the per-axis spread: E[r²] = d·s² = 2·d·D·Δt, so D = s²/2 in 2D and 3D,
  # after removing noise and blur (see the header of this function). Only meaningful
  # with √Δt normalisation; without it s is in nm per step and D keeps that unit.
  if (norm_dt && is.finite(inv_dt)) {
    loc_bias  <- if (sigma_loc > 0) 2 * sigma_loc^2 * inv_dt else 0
    blur_fac  <- max(0.05, 1 - 2 * max(0, R_abs) * inv_dt)   # guard: never divide by ~0
    D_vals <- pmax(0, sigs^2 - loc_bias) / (2 * blur_fac)
  } else {
    D_vals <- sigs^2 / 2
  }

  # Noise floor as a JUMP DISTANCE in the plot's units: RMS √(2d)·σ_loc, ×√⟨1/Δt⟩ if normalised.
  r_noise <- if (sigma_loc > 0) {
    base <- sqrt(2 * d_jdd) * sigma_loc
    if (norm_dt && is.finite(inv_dt)) base * sqrt(inv_dt) else base
  } else NA_real_

  list(mids = mids, dens = dens, sigs = sigs, ws = ws, D_vals = D_vals,
       dt0 = dt0, inv_dt = inv_dt, r_noise = r_noise, d_fit = d_jdd)
}
