## ─────────────────────────────────────────────────────────────────────────────
## est_msd.R  -  MSD model fitting (OLSF linear fit and the non-linear models)
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The fitting core behind the Tracking tab's MSD panels, as pure functions so
##   tests/ can check them against simulated ground truth. Moved verbatim out of
##   tracking_server() on 2026-09-17 (publication audit C1); the server still
##   builds the MSD curves (msd_unified_bundle) and calls these to fit them.
##
## WHAT LIVES HERE
##   MSD_BLUR_R ............... motion-blur coefficient shared by the MSD fits, the
##                              confinement null and the HMM
##   .olsf_sig2() ............. intercept -> σ² with the blur offset added back
##   .michalet_p() ............ Michalet & Berglund optimal point count
##   msd_fit_olsf() ........... anchored OLSF (app default)
##   msd_fit_olsf_iterative() . MinSpt's exact iterative OLSF
##   .msd_model_fit() ......... linear / confined / directed / anomalous dispatcher
##
## UNITS: τ in ms, MSD in nm², D in nm²/ms (÷1000 → µm²/s), σ in nm.
##
## ATTRIBUTION: the OLSF workflow is an R port of MinSpt (Vogler et al. 2025), which is
## GPL-3.0; that is why MINFLUXR as a whole is GPL-3.0. See THIRD-PARTY.md section 1.
## ─────────────────────────────────────────────────────────────────────────────

# --- Shared MSD fit helpers --------------------------------------------------
# ndim = 2 (2D): MSD = 4*D*tau^alpha + 4*sigma^2
# ndim = 3 (3D): MSD = 6*D*tau^alpha + 6*sigma^2
# Units throughout: tau in ms, MSD in nm², D in nm²/ms.
# Display conversion: D [µm²/s] = D [nm²/ms] / 1e3   (1 nm²/ms = 0.001 µm²/s)
#                     MSD [µm²] = MSD [nm²]   / 1e6

# OLSF - faithful port of Vogler's OptimizedLeastSquaresFitting (Michalet &
# Berglund). UNWEIGHTED linear fit MSD = a + b·τ with SEPARATE optimal point
# counts for the slope (→ D) and the intercept (→ σ²), iterated to a fixed
# point.  Edge cases and formulas match his olsf_base_math / olsf_fit_utilities:
#   x   = σ²/(D·dt) − 2R                 (reduced localization error)
#   a   = 2d·σ² − 4d·R·D·dt   (a<0 ⟺ x<0 → p = 2)
#   b   = 2d·D·dt            (b<0 → p = L)
#   f_b = 2 + 1.35·x^0.6 ,  L_b = 0.8 + 0.564·N            (slope)
#   f_a = 2 + 1.6·x^0.51 ,  L_a = 3 + (4.5·N^0.4 − 8.5)^1.2 (intercept)
#   p   = f·L / (f³ + L³)^(1/3)
.michalet_p <- function(D, sig2, dt, N, R, d, which) {
  if (!is.finite(D) || D <= 0 || !is.finite(dt) || dt <= 0) return(NA_real_)
  a <- 2*d*sig2 - 4*d*R*D*dt
  b <- 2*d*D*dt
  if (which == "slope") {
    if (a < 0) return(2); if (b < 0) return(0.8 + 0.564*N)
    x <- sig2/(D*dt) - 2*R; if (x < 0) x <- 0
    fb <- 2 + 1.35*x^0.6; Lb <- 0.8 + 0.564*N
    fb*Lb/((fb^3 + Lb^3)^(1/3))
  } else {
    if (a < 0) return(2); if (b < 0) return(3 + (4.5*N^0.4 - 8.5)^1.2)
    x <- sig2/(D*dt) - 2*R; if (x < 0) x <- 0
    fa <- 2 + 1.6*x^0.51; La <- 3 + (4.5*N^0.4 - 8.5)^1.2
    fa*La/((fa^3 + La^3)^(1/3))
  }
}

# Motion-blur coefficient R for the MSD fits. A MINFLUX localization averages the
# molecule's path over the cycle, which LOWERS every lag's MSD by the same offset
# (Michalet & Berglund 2012):
#     MSD(τ) = coef·D·τ + coef·σ² − 2·coef·R·D·dT        (dT = cycle time)
# so the intercept is coef·(σ² − 2·R·D·dT), not coef·σ². Reading σ² straight off
# the intercept (R = 0) under-reports σ_loc: at D = 0.215 µm²/s, t_c = 0.192 ms and
# σ = 5 nm it returns 3.3 nm (simulated, full-cycle window). D is unaffected.
# 1/6.2 is the value the ported toolbox uses (Vogler et al., citing Schmidt et al.
# 2020); the CVE's instrument-modelled bound [R_lo, R_hi] brackets it, and the
# difference moves σ² by < 2 nm² at these timings. A constant keeps the MSD panel
# independent of whether the CVE has been run.
MSD_BLUR_R <- 1 / 6.2

# Convert a linear fit's intercept to σ² by adding back the blur offset.
# GOAL: one place for the intercept → σ² relation so every OLSF variant agrees.
# intercept_over_coef = cf[1]/coef; D = slope/coef (nm²/ms); dT = cycle time (ms).
.olsf_sig2 <- function(intercept_over_coef, D, R, dT) {
  blur <- if (is.finite(R) && R > 0 && is.finite(D) && D > 0 && is.finite(dT) && dT > 0)
            2 * R * D * dT else 0          # blur offset is only defined for a positive slope
  intercept_over_coef + blur
}

# auto_pmin = FALSE: fit the first k points (caller-chosen window).
# auto_pmin = TRUE : two-stage iterative Michalet optimal point selection.
# `dt_cycle` is the cycle time (for x and the blur offset); `Npts` the trajectory
# length for L(N). `R` is the motion-blur coefficient (MSD_BLUR_R at every call site).
# Returns sigma_loc blur-corrected, and `intercept` (nm², the fitted offset) which is
# what the drawn line uses, since the data sit on the blurred line.
# coef_override: the MSD-vs-D coefficient, when it is NOT simply 2*ndim. A molecule
# confined to a SURFACE moves in 2 degrees of freedom however many coordinates are
# recorded, and a 2D projection of a curved surface is foreshortened by a
# geometry-dependent factor, so the divisor is not the recorded dimensionality (see
# .msd_d_eff). `ndim` still carries the GEOMETRIC dimension, because Michalet's
# optimal point count is derived for real spatial dimensions.
msd_fit_olsf <- function(tau, msd, npairs, k, ndim = 2L,
                         auto_pmin = FALSE, R = 0, Npts = NA_real_,
                         dt_cycle = NA_real_, coef_override = NULL) {
  coef_val <- if (is.null(coef_override)) 2L * ndim else coef_override
  ok_all <- is.finite(tau) & is.finite(msd) & tau > 0
  tau <- tau[ok_all]; msd <- msd[ok_all]
  n_avail <- length(tau)
  if (n_avail < 3L) return(NULL)
  # Cycle time for the Michalet x and the blur offset (fallback: the lag spacing).
  dt <- if (is.finite(dt_cycle) && dt_cycle > 0) dt_cycle
        else if (n_avail >= 2) stats::median(diff(tau)) else tau[1]

  # Unweighted linear fit on the first p points → D, σ² (his curve_fit).
  # sig2 is the blur-corrected σ²; icpt keeps the raw fitted offset for drawing.
  fit_p <- function(p) {
    p  <- max(2L, min(as.integer(p), n_avail))
    s <- tau[seq_len(p)]; y <- msd[seq_len(p)]
    ok <- is.finite(s) & is.finite(y)
    if (sum(ok) < 2L) return(NULL)
    cf <- tryCatch(coef(lm(y[ok] ~ s[ok])), error = function(e) NULL)
    if (is.null(cf) || any(!is.finite(cf))) return(NULL)
    D_p <- cf[2] / coef_val
    list(D = D_p, sig2 = .olsf_sig2(cf[1] / coef_val, D_p, R, dt),
         icpt = unname(cf[1]), p = p)
  }

  if (!auto_pmin) {
    r <- fit_p(max(3L, min(as.integer(k), n_avail)))
    if (is.null(r)) return(NULL)
    D_est <- max(0, r$D); sig2_est <- max(0, r$sig2); p_slope <- r$p; p_int <- r$p
    icpt_est <- r$icpt
  } else {
    N  <- if (is.finite(Npts) && Npts > 0) Npts else n_avail
    d  <- ndim
    # Anchor the reduced localization error x on a STABLE fit over the whole
    # [min_lag, fit_max] window.  Michalet's iteration is otherwise init-
    # dependent: a few-point fit gives a too-small σ² (hence too-small x) and
    # collapses p to the starting guess.  The full-window σ² is the stable
    # baseline that makes the optimal point count well-defined and reproducible.
    r_anc   <- fit_p(n_avail)
    if (is.null(r_anc)) return(NULL)
    D_anc   <- max(1e-9, r_anc$D)
    sig2_anc<- max(0, r_anc$sig2)
    # Separate optimal point counts for the slope (→ D) and intercept (→ σ²),
    # per Michalet & Berglund (Vogler's p_b / p_a), with his +1 shift.
    p_b <- .michalet_p(D_anc, sig2_anc, dt, N, R, d, "slope")
    p_a <- .michalet_p(D_anc, sig2_anc, dt, N, R, d, "intercept")
    p_slope <- if (is.finite(p_b)) max(3L, min(n_avail, as.integer(round(p_b)) + 1L)) else n_avail
    p_int   <- if (is.finite(p_a)) max(3L, min(n_avail, as.integer(round(p_a)) + 1L)) else n_avail
    r_s <- fit_p(p_slope); r_i <- fit_p(p_int)
    D_est    <- if (!is.null(r_s)) max(0, r_s$D)    else max(0, D_anc)
    sig2_est <- if (!is.null(r_i)) max(0, r_i$sig2) else sig2_anc
    icpt_est <- if (!is.null(r_i)) r_i$icpt else r_anc$icpt
  }
  if (!is.finite(D_est)) return(NULL)
  if (!is.finite(sig2_est)) sig2_est <- 0
  if (!is.finite(icpt_est)) icpt_est <- coef_val * sig2_est
  sig_est <- sqrt(sig2_est)
  tau_seq <- seq(min(tau), max(tau), length.out = 200)
  # Draw the line the data actually follow: slope from D, offset = fitted intercept
  # (blurred), so the overlay does not float above the points by the blur term.
  yfit    <- pmax(0, coef_val * D_est * tau_seq + icpt_est)

  list(D = D_est, alpha = 1.0, sigma_loc = sig_est, intercept = icpt_est,
       tau_fit = tau_seq, msd_fit = yfit, n_points = p_slope,
       n_points_intercept = p_int, R_blur = R,
       forced_linear = TRUE, ndim = ndim, method = "OLSF")
}

# Bela's EXACT iterative OLSF (optimized_least_squares_fitting.py). Unlike the
# anchored variant above, it does NOT pre-fit the whole window - it starts at a
# few points and iterates the optimal point count until it repeats, separately
# for slope (→ D) and intercept (→ σ²). This deliberately converges on the
# LOWEST lags (least error-prone), reproducing his numbers. `tau` is the lag
# axis (use the linearized lag×cycle axis for a faithful match); D = slope/(2d).
# coef_override as in msd_fit_olsf.
msd_fit_olsf_iterative <- function(tau, msd, ndim = 2L, R = 0, coef_override = NULL,
                                   dt_cycle = NA_real_, N = NA_real_,
                                   max_attempts = 1000L) {
  coef_val <- if (is.null(coef_override)) 2L * ndim else coef_override
  ok <- is.finite(tau) & is.finite(msd) & tau > 0
  tau <- tau[ok]; msd <- msd[ok]
  n_avail <- length(tau)
  if (n_avail < 3L) return(NULL)
  d  <- ndim
  if (!is.finite(N) || N <= 0) N <- n_avail
  if (!is.finite(dt_cycle) || dt_cycle <= 0)
    dt_cycle <- if (n_avail >= 2) stats::median(diff(tau)) else tau[1]

  # sig2 is blur-corrected (his curve_fit fits σ² inside a = 2dσ² − 4dRDdT);
  # icpt is the raw fitted offset, used to draw the line through the data.
  fit_p <- function(p) {
    p <- max(2L, min(as.integer(p), n_avail))
    s <- tau[seq_len(p)]; y <- msd[seq_len(p)]
    cf <- tryCatch(stats::coef(stats::lm(y ~ s)), error = function(e) NULL)
    if (is.null(cf) || any(!is.finite(cf))) return(NULL)
    D_p <- cf[2] / coef_val
    list(D = D_p, sig2 = .olsf_sig2(cf[1] / coef_val, D_p, R, dt_cycle),
         icpt = unname(cf[1]), p = p)
  }

  p0 <- if (n_avail %/% 10L > 3L) n_avail %/% 10L else 3L
  nps_s <- p0; nps_i <- p0
  hist_s <- p0; hist_i <- p0
  s_done <- FALSE; i_done <- FALSE
  D_opt <- NA_real_; sig2_opt <- NA_real_
  p_s_fin <- p0; p_i_fin <- p0

  for (it in seq_len(max_attempts)) {
    rs <- if (!s_done) fit_p(nps_s) else NULL
    ri <- if (!i_done) fit_p(nps_i) else NULL
    if ((!s_done && is.null(rs)) || (!i_done && is.null(ri))) break
    # Cross-coupled (his exact wiring): σ² from the INTERCEPT fit, D from SLOPE.
    D_cur    <- if (!s_done) rs$D    else D_opt
    sig2_cur <- if (!i_done) ri$sig2 else sig2_opt
    if (!s_done) {
      pb <- .michalet_p(max(1e-9, D_cur), max(0, sig2_cur), dt_cycle, N, R, d, "slope")
      new_s <- if (is.finite(pb)) as.integer(round(pb)) + 1L else n_avail
      if (new_s %in% hist_s) { D_opt <- rs$D; p_s_fin <- nps_s; s_done <- TRUE }
      else { hist_s <- c(hist_s, new_s); nps_s <- max(2L, min(n_avail, new_s)) }
    }
    if (!i_done) {
      pa <- .michalet_p(max(1e-9, D_cur), max(0, sig2_cur), dt_cycle, N, R, d, "intercept")
      new_i <- if (is.finite(pa)) as.integer(round(pa)) + 1L else n_avail
      if (new_i %in% hist_i) { sig2_opt <- ri$sig2; p_i_fin <- nps_i; i_done <- TRUE }
      else { hist_i <- c(hist_i, new_i); nps_i <- max(2L, min(n_avail, new_i)) }
    }
    if (s_done && i_done) break
  }
  if (!is.finite(D_opt))    { r <- fit_p(nps_s); if (!is.null(r)) D_opt <- r$D }
  if (!is.finite(sig2_opt)) { r <- fit_p(nps_i); if (!is.null(r)) sig2_opt <- r$sig2 }
  if (!is.finite(D_opt)) return(NULL)
  D_est <- max(0, D_opt); sig2_est <- max(0, if (is.finite(sig2_opt)) sig2_opt else 0)
  sig_est <- sqrt(sig2_est)
  # Intercept of the converged intercept fit (blurred offset) for the drawn line.
  icpt_est <- { r <- fit_p(p_i_fin); if (!is.null(r)) r$icpt else coef_val * sig2_est }
  tau_seq <- seq(min(tau), max(tau), length.out = 200)
  yfit <- pmax(0, coef_val * D_est * tau_seq + icpt_est)
  list(D = D_est, alpha = 1.0, sigma_loc = sig_est, intercept = icpt_est,
       tau_fit = tau_seq, msd_fit = yfit, n_points = p_s_fin,
       n_points_intercept = p_i_fin, R_blur = R,
       forced_linear = TRUE, ndim = ndim, method = "OLSF-iterative")
}

# Ensemble-MSD model fit with a selectable model. "linear" = the OLSF above
# (Brownian, D from the slope). "confined" = a plateauing corral model (→ local
# D + confinement size L). "directed" = Brownian + drift (→ D + velocity v).
# "anomalous" = power law (→ D + exponent α). The non-linear models fit the WHOLE
# [min_lag, max_lag] window (they need the curvature/plateau, not just low lags).
# Returns a list with $model, $D (nm²/ms), $sigma_loc (nm), $tau_fit, $msd_fit
# (nm²) and model-specific extras (L_conf, tau_c_ms, plateau_um2, v_um_s, alpha).
# Pure model-fit core, reused by the unified plot AND the per-segment analysis.
# s = τ (ms, already windowed), y = MSD (nm²). model ∈ linear/confined/directed/
# anomalous. Returns $model, $D (nm²/ms), $sigma_loc (nm), $tau_fit, $msd_fit (nm²)
# + extras (L_conf, tau_c_ms, plateau_um2, degenerate / v_um_s / alpha).
.msd_model_fit <- function(s, y, npairs, ndim, model, use_iter, t_cycle,
                           coef_override = NULL) {
  coef_v <- if (is.null(coef_override)) 2L * ndim else coef_override
  if (identical(model, "linear")) {
    if (length(s) < 3L) return(NULL)
    return(if (isTRUE(use_iter))
             msd_fit_olsf_iterative(s, y, ndim = ndim, R = MSD_BLUR_R, dt_cycle = t_cycle,
                                    N = length(s), coef_override = coef_override)
           else
             msd_fit_olsf(s, y, npairs, k = length(s), ndim = ndim,
                          auto_pmin = TRUE, R = MSD_BLUR_R, Npts = 1e4, dt_cycle = t_cycle,
                          coef_override = coef_override))
  }
  if (length(s) < 4L) return(NULL)
  # `warnOnly = TRUE` makes nls RETURN a fit instead of erroring when it cannot
  # converge - but it still prints "Convergence failure: singular convergence (7)"
  # to the console, once per attempt, which is noise the user cannot act on (each
  # cluster × each of the three non-linear models is an attempt). The failure is not
  # ignored: it is recorded on the fit and folded into the `degenerate` flag below,
  # which is what the UI already reports. Singular convergence means the parameters
  # are not identifiable over this window - usually no plateau has been reached -
  # so the coefficients that come back must not be presented as a real fit.
  nls_fit <- function(formula, start, lower, upper) {
    f <- tryCatch(suppressWarnings(
      stats::nls(formula, data = data.frame(s = s, y = y), start = start,
                 lower = lower, upper = upper, algorithm = "port",
                 control = stats::nls.control(maxiter = 200, warnOnly = TRUE))),
      error = function(e) NULL)
    if (is.null(f)) return(NULL)
    attr(f, "converged") <- isTRUE(tryCatch(f$convInfo$isConv, error = function(e) NA))
    f
  }
  D0   <- max(1e-6, (max(y) - min(y)) / (coef_v * max(s)))
  C0   <- max(0, min(y, na.rm = TRUE))
  tseq <- seq(min(s), max(s), length.out = 200)
  if (identical(model, "confined")) {
    fit <- nls_fit(y ~ A * (1 - exp(-s / tc)) + C,
                   start = list(A = max(y) * 0.9, tc = max(s) / 3, C = C0),
                   lower = list(A = 0, tc = 1e-4, C = 0),
                   upper = list(A = Inf, tc = Inf, C = Inf))
    if (is.null(fit)) return(NULL)
    cf <- coef(fit); A <- unname(cf["A"]); tc <- unname(cf["tc"]); C <- unname(cf["C"])
    # A non-converged fit is degenerate too: singular convergence means the plateau
    # was never resolved, which is exactly what this flag already communicates.
    degenerate <- (tc > 2 * max(s)) || (A > 3 * max(y, na.rm = TRUE)) ||
                  !isTRUE(attr(fit, "converged"))
    # Confinement size L_c from the fitted plateau A. For a particle uniformly
    # exploring a d-dimensional box of side L_c, each axis contributes
    # Var(Δx) = L_c²/6, so the steady-state MSD is A = d·L_c²/6 and
    #     L_c = sqrt(6·A/d)
    # which gives the familiar sqrt(3A) in 2D and sqrt(2A) in 3D. Verified against
    # reflecting-boundary simulations (d=2 and d=3, L_c = 140/200/300 nm, recovered
    # to <1%). Using the 2D form in 3D would inflate L_c by sqrt(3/2) = 1.22×.
    list(model = "confined", D = A / (coef_v * tc), sigma_loc = sqrt(max(0, C / coef_v)),
         L_conf = sqrt(6 * A / ndim), tau_c_ms = tc, plateau_um2 = A / 1e6, degenerate = degenerate,
         tau_fit = tseq, msd_fit = A * (1 - exp(-tseq / tc)) + C, ndim = ndim, method = "confined")
  } else if (identical(model, "directed")) {
    fit <- nls_fit(y ~ coef_v * D * s + (v * s)^2 + C,
                   start = list(D = D0, v = sqrt(D0) / max(s), C = C0),
                   lower = list(D = 0, v = 0, C = 0),
                   upper = list(D = Inf, v = Inf, C = Inf))
    if (is.null(fit)) return(NULL)
    cf <- coef(fit); D <- unname(cf["D"]); v <- unname(cf["v"]); C <- unname(cf["C"])
    list(model = "directed", D = D, sigma_loc = sqrt(max(0, C / coef_v)), v_um_s = v,
         degenerate = !isTRUE(attr(fit, "converged")),
         tau_fit = tseq, msd_fit = coef_v * D * tseq + (v * tseq)^2 + C, ndim = ndim, method = "directed")
  } else if (identical(model, "anomalous")) {
    fit <- nls_fit(y ~ coef_v * D * s^a + C,
                   start = list(D = D0, a = 1, C = C0),
                   lower = list(D = 0, a = 0.05, C = 0),
                   upper = list(D = Inf, a = 2.5, C = Inf))
    if (is.null(fit)) return(NULL)
    cf <- coef(fit); D <- unname(cf["D"]); a <- unname(cf["a"]); C <- unname(cf["C"])
    list(model = "anomalous", D = D, alpha = a, sigma_loc = sqrt(max(0, C / coef_v)),
         degenerate = !isTRUE(attr(fit, "converged")),
         tau_fit = tseq, msd_fit = coef_v * D * tseq^a + C, ndim = ndim, method = "anomalous")
  } else NULL
}
