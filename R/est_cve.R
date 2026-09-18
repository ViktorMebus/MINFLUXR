## ─────────────────────────────────────────────────────────────────────────────
## est_cve.R  -  covariance-based estimator (CVE) of D and σ_loc, photon-model σ
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The core of the Tracking tab's "Data Quality & Calibration" panel as pure
##   functions. Moved out of trk_cve_bundle() on 2026-09-17 (publication audit
##   C1). The server still builds the diagnostics around it (continuity, per-track
##   σ, linking QC, photon/mobility coupling) from the step arrays it computes.
##
## WHAT LIVES HERE
##   cve_steps() ........ step arrays (dx, dy, dz, dt, single-cycle mask) per row
##   cve_core() ......... lag-1 covariance, R bounds (staircase photon-window
##                        model .blur_R_duty() nested inside), blur-aware σ, D_cve, r1
##   cve_estimate() ..... cve_steps() + cve_core() from a data.frame (tests, scripts)
##   .crlb_doughnut() ... photon-budget σ per localization (Vogler et al. 2025 Eq. 1)
##
## METHOD (Vestergaard, Blainey & Flyvbjerg 2014, with motion blur)
##   Over consecutive single-cycle step pairs, per axis:
##       cov1 = <dx_i dx_(i+1)>,   σ² = R·<dx²> + (2R − 1)·cov1,
##       D = <dx²>/(2Δt) + cov1/Δt  (blur-independent)
##   R is bounded between a staircase model of the photon window (R_lo, used)
##   and a full-cycle window t_c/(6Δt) capped at 1/6 (R_hi).
##
## UNITS: nm, ms, D nm²/ms.
## ─────────────────────────────────────────────────────────────────────────────


# ── cve_steps() ───────────────────────────────────────────────────────────────
# GOAL: per-row step arrays for a data.frame sorted by (tid, tim). Row i holds
# the step arriving at localization i (NA for a track's first localization).
# RETURNS list(df (sorted), n, same, dt, dx, dy, dz, n_cyc, integ_step, n_steps, frac_integ)
cve_steps <- function(df, t_c, has_z = FALSE) {
  df <- df[order(df$tid, df$tim), ]
  n  <- nrow(df)
  same <- c(FALSE, df$tid[-1] == df$tid[-n])                 # step i valid (2nd+ loc)
  dt   <- c(NA_real_, diff(df$tim)); dt[!same] <- NA_real_
  dx   <- c(NA_real_, diff(df$x));   dx[!same] <- NA_real_
  dy   <- c(NA_real_, diff(df$y));   dy[!same] <- NA_real_
  dz   <- if (has_z) { v <- c(NA_real_, diff(df$z)); v[!same] <- NA_real_; v } else rep(NA_real_, n)
  n_cyc <- if (is.finite(t_c) && t_c > 0) round(dt / t_c) else rep(1L, n)
  integ_step <- same & is.finite(n_cyc) & (n_cyc == 1L)
  integ_step[is.na(integ_step)] <- FALSE
  n_steps    <- sum(same, na.rm = TRUE)
  frac_integ <- if (n_steps > 0) sum(integ_step) / n_steps else NA_real_
  list(df = df, n = n, same = same, dt = dt, dx = dx, dy = dy, dz = dz, n_cyc = n_cyc,
       integ_step = integ_step, n_steps = n_steps, frac_integ = frac_integ)
}

# ── cve_core() ────────────────────────────────────────────────────────────────
# GOAL: the CVE numbers from the step arrays.
# dx, dy, dz, dt, integ_step : from cve_steps()
# has_z  : include z (only when z carries real axial data)
# t_c    : cycle time (ms)
# t_coll : photon collection window per localization (ms), median ECO/EFO; NA = unknown
# n_vert : pattern vertices (6 = hexagon)
# RETURNS list(pair, prod_x, prod_y, prod_z, n_pairs, cov1, dt_int, R_lo, R_hi, R_blur,
#   msq, sigma2, sigma, sigma_xy, sigma_z, sig_aniso, D_cve, r1_norm, noise_frac),
#   or NULL when there are fewer than 20 consecutive single-cycle step products.
cve_core <- function(dx, dy, dz, dt, integ_step, has_z = FALSE, t_c = NA_real_,
                     t_coll = NA_real_, n_vert = 6L) {
  n <- length(dx)
  # CVE lag-1 autocovariance over CONSECUTIVE integrated step pairs.
  lead_integ <- c(integ_step[-1], FALSE)
  pair   <- integ_step & lead_integ
  prod_x <- dx * c(dx[-1], NA_real_)
  prod_y <- dy * c(dy[-1], NA_real_)
  vals   <- c(prod_x[pair], prod_y[pair])
  prod_z <- if (has_z) dz * c(dz[-1], NA_real_) else rep(NA_real_, n)
  if (has_z) vals <- c(vals, prod_z[pair])
  vals <- vals[is.finite(vals)]
  if (length(vals) < 20) return(NULL)             # caller reports "too few steps"
  cov1   <- mean(vals)
  # ── Motion blur (Vestergaard's R) ──────────────────────────────────────────
  # A MINFLUX localization is not a snapshot: photons are collected across the
  # whole TCP cycle, so the reported position is a time-average of the molecule's
  # path over that window. Vestergaard's CVE then reads
  #     sigma^2 = R*<dx^2> + (2R - 1)*cov1        (D is blur-INDEPENDENT)
  # with R = t_collect/(6*dt). R = 0 is the instantaneous-exposure special case,
  # where it collapses to sigma^2 = -cov1 - which is what this used to compute,
  # and it under-reports sigma by ~3x on real tracking data (1.5 vs 4.1 nm on the
  # GP67H files). Collection spans the cycle t_c out of a step dt = t_c + SW, so
  # R = t_c/(6*dt); a full-cycle window at SW = 0 gives the familiar 1/6.
  dt_int  <- mean(dt[integ_step], na.rm = TRUE)
  # R is not a single knowable number: it depends on how the cycle's photons are
  # distributed in time, and only part of that is visible in the data. So bound it.
  #
  #   R_hi - collection spans the whole cycle: R = t_c/(6*dt), capped at 1/6.
  #   R_lo - collection is only the measured ECO/EFO window, delivered as a
  #          staircase over the pattern's vertices, with the beam moves between
  #          them and the terminal localization step contributing nothing.
  #          Closed form for K equal dwells of a, moves of m, terminal r:
  #            R = [ sum_i a*(s0(1-s0) + d(1-2*s0)/2 - d^2/3) + sum_i m*(i/K)(1-i/K) ] / dt
  #          with d = 1/K and s0 = (i-1)/K.
  # Vertex count is a weak lever (triangle 0.131 vs hexagon 0.137 at these
  # timings); the terminal step is the strong one (0.165 without it, 0.137 with).
  .blur_R_duty <- function(dt_ms, t_coll_ms, K, t_eod_ms = 0.0053, t_rep_ms = 0.020) {
    if (!is.finite(dt_ms) || dt_ms <= 0 || !is.finite(t_coll_ms) || t_coll_ms <= 0 ||
        !is.finite(K) || K < 1) return(NA_real_)
    t_coll_ms <- min(t_coll_ms, dt_ms)
    a  <- t_coll_ms / K
    m  <- max(0, min(t_eod_ms, (dt_ms - t_coll_ms) / K))
    d  <- 1 / K; i <- seq_len(K); s0 <- (i - 1) / K
    dwell_part <- sum(a * (s0 * (1 - s0) + d * (1 - 2 * s0) / 2 - d^2 / 3))
    flat_part  <- sum(m * (i / K) * (1 - i / K))
    (dwell_part + flat_part) / dt_ms
  }
  R_hi <- if (is.finite(t_c) && t_c > 0 && is.finite(dt_int) && dt_int > 0)
            min(1/6, t_c / (6 * dt_int)) else 1/6
  R_lo <- .blur_R_duty(dt_int, t_coll, n_vert)
  if (!is.finite(R_lo) || R_lo <= 0 || R_lo > R_hi) R_lo <- R_hi
  R_blur <- R_lo   # the physically modelled value; R_hi bounds the other side
  sq_all  <- c(dx[integ_step]^2, dy[integ_step]^2)
  if (has_z) sq_all <- c(sq_all, dz[integ_step]^2)
  sq_all  <- sq_all[is.finite(sq_all)]
  msq     <- if (length(sq_all)) mean(sq_all) else NA_real_
  .sig_from <- function(msq_v, cov_v)
    if (!is.finite(msq_v) || !is.finite(cov_v)) NA_real_
    else sqrt(max(0, R_blur * msq_v + (2 * R_blur - 1) * cov_v))
  sigma2 <- R_blur * msq + (2 * R_blur - 1) * cov1
  sigma  <- sqrt(max(0, sigma2))
  # ── Per-axis-group σ: LATERAL (x,y) vs AXIAL (z) ────────────────────────────
  # `sigma` above pools all axes into ONE number. In 3D MINFLUX the axial precision
  # is typically WORSE than the lateral, so the pooled value over-states σ_xy and
  # under-states σ_z. Several consumers are purely lateral (the D-map is an x–y map;
  # JDD/MSD in 2D mode), and subtracting a z-inflated noise floor from them
  # over-corrects. Report both so the anisotropy is visible and quantified.
  .sig_of <- function(v, sqv) {
    v <- v[is.finite(v)]; sqv <- sqv[is.finite(sqv)]
    if (length(v) < 10L || !length(sqv)) NA_real_ else .sig_from(mean(sqv), mean(v))
  }
  sigma_xy <- .sig_of(c(prod_x[pair], prod_y[pair]),
                      c(dx[integ_step]^2, dy[integ_step]^2))
  sigma_z  <- if (has_z) .sig_of(prod_z[pair], dz[integ_step]^2) else NA_real_
  # Axial/lateral ratio - > ~1.3 means a single pooled σ_loc is a poor description.
  sig_aniso <- if (is.finite(sigma_xy) && sigma_xy > 0 && is.finite(sigma_z))
                 sigma_z / sigma_xy else NA_real_
  # D is blur-independent: D = ⟨Δ²⟩/(2Δt) + cov1/Δt holds for any R.
  sq     <- sq_all
  D_cve  <- if (length(sq) && is.finite(dt_int) && dt_int > 0) mean(sq)/(2*dt_int) + cov1/dt_int else NA_real_

  # ── How well is σ actually resolved? (drives the photon-weighting gate) ─────
  # Normalized lag-1 correlation r1 = cov1/⟨Δ²⟩. For Brownian + noise,
  # r1 = −σ²/(2·D·Δt + 2σ²), so the NOISE FRACTION of the step variance is −2·r1.
  # r1 → −0.5 means noise dominates (σ cleanly resolved); r1 → 0 means σ is a marginal
  # residual sitting on top of motion, and any PER-TRACK or PER-BIN σ then mostly
  # tracks the step variance (mobility), not photons. r1 > 0 means correlated/directed
  # motion is dominating and σ is unresolvable there. The pooled σ above is still the
  # best estimate; this only says how much to trust σ HETEROGENEITY between tracks.
  r1_norm    <- if (length(sq) && mean(sq) > 0) cov1 / mean(sq) else NA_real_
  noise_frac <- -2 * r1_norm
  list(pair = pair, prod_x = prod_x, prod_y = prod_y, prod_z = prod_z, n_pairs = length(vals),
       cov1 = cov1, dt_int = dt_int, R_lo = R_lo, R_hi = R_hi, R_blur = R_blur, msq = msq,
       sigma2 = sigma2, sigma = sigma, sigma_xy = sigma_xy, sigma_z = sigma_z,
       sig_aniso = sig_aniso, D_cve = D_cve, r1_norm = r1_norm, noise_frac = noise_frac)
}

# ── cve_estimate() ────────────────────────────────────────────────────────────
# GOAL: CVE from a data.frame in one call (tests and scripts). Returns cve_core()'s
# list plus n_steps and frac_integ, or NULL.
cve_estimate <- function(df, t_c, has_z = FALSE, t_coll = NA_real_, n_vert = 6L) {
  st <- cve_steps(df, t_c, has_z)
  cc <- cve_core(st$dx, st$dy, st$dz, st$dt, st$integ_step, has_z, t_c, t_coll, n_vert)
  if (is.null(cc)) return(NULL)
  c(cc, list(n_steps = st$n_steps, frac_integ = st$frac_integ))
}

# ── Doughnut CRLB, Vogler et al. 2025 Eq. 1 ─────────────────────────────────
# sigma(x) = L/(4 sqrt(N)) * [1 + 2(x/L)^2] * |cosh(a) - (x/L) sinh(a)|
#            / |1 + ln2 (L/nu)^2 ((2x/L)^2 - 1)| ,   a = 4 ln2 x L / nu^2
#
# L  TCP diameter of the reiterated iteration
# nu FWHM of the excitation doughnut (instrument constant)
# N  detected photons for this localization (eco)
# x  emitter offset from the TCP centre. Zero for an immobile emitter; while
#    tracking, Vogler model it as the 1D diffusive spread over one localization,
#    x = sqrt(2 D t_loc), which is why SPT precision is worse than imaging.
#
# Vectorised over N. Returns L/(4 sqrt(N)) when nu is not usable, i.e. the
# nu -> infinity limit, so a missing constant degrades to the old behaviour.
.crlb_doughnut <- function(L, nu, N, x = 0, pref = TCP_CRLB_PREFACTOR) {
  if (!length(N)) return(numeric(0))
  if (!is.finite(pref) || pref <= 0) pref <- TCP_CRLB_PREFACTOR
  base <- L / (pref * sqrt(N))
  if (!is.finite(L) || L <= 0) return(base)
  if (!is.finite(x) || x < 0) x <- 0
  # nu off (0) = the large-beam limit. The denominator goes to 1 and cosh/sinh to
  # 1/0, so Eq. 1 collapses to base * [1 + 2(x/L)^2] - the moving-emitter term
  # SURVIVES. That term is what makes the estimate mode-aware: x = sqrt(2 D t_loc)
  # is smaller under confinement and larger under directed or fast motion.
  if (!is.finite(nu) || nu <= 0) return(base * (1 + 2 * (x / L)^2))
  a   <- 4 * log(2) * x * L / nu^2
  num <- (1 + 2 * (x / L)^2) * abs(cosh(a) - (x / L) * sinh(a))
  den <- abs(1 + log(2) * (L / nu)^2 * ((2 * x / L)^2 - 1))
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(base)
  base * num / den
}
