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
##                        R from blur_R_step() in est_blur.R), blur-aware σ, D_cve, r1
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
                     t_coll = NA_real_, n_vert = 6L, R_step = NA_real_) {
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
  # tracking test files). Collection spans the cycle t_c out of a step dt = t_c + SW, so
  # R = t_c/(6*dt); a full-cycle window at SW = 0 gives the familiar 1/6.
  dt_int  <- mean(dt[integ_step], na.rm = TRUE)
  # R depends on WHEN the photons that make up a localization are collected inside
  # the step (vertices, dwell, moves, CFR centre, re-centring), not on L; the model
  # and its validation are in est_blur.R. R_hi (collection spanning the whole cycle,
  # t_c/(6*dt), capped at 1/6) is kept as the upper bound.
  R_hi <- if (is.finite(t_c) && t_c > 0 && is.finite(dt_int) && dt_int > 0)
            min(1/6, t_c / (6 * dt_int)) else 1/6
  # R from the shared timing model (est_blur.R) so the CVE, MSD, JDD, D-map,
  # change-point and HMM all use one blur. `R_step` from the caller (built from the
  # sequence file when there is one) wins; otherwise the timeline is built from the
  # measured collection window t_coll = eco/efo and the vertex count.
  R_lo <- if (is.finite(R_step) && R_step > 0) R_step else
            (blur_R_step(dt_int, t_coll, n_vert = n_vert)$R %||% NA_real_)
  if (!is.finite(R_lo) || R_lo <= 0 || R_lo > R_hi) R_lo <- min(R_hi, BLUR_R_STEP_FALLBACK)
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

# ── Emitter offset from the TCP centre, per pattern geometry ───────────────────
# During tracking the emitter moves away from the TCP centre while photons are
# collected, and precision falls off with that offset. Per-axis variance, averaged
# over axes and over the direction of the offset, follows
#     sigma^2(r) = sigma0^2 * (1 + c * <r^2> / L^2)
# with r the radial offset. c comes from the Fisher matrix of the multinomial photon
# split for an ideal doughnut zero, fitted for r/L <= 0.15 (tests/testthat/
# test-estimators-tracking.R recomputes it). The triangle loses precision fastest
# because it samples the fewest directions. `dim` is the number of axes the
# pattern spans, which sets <r^2> = dim * D * t_loc (the offset averaged over the
# exposure, since photons arrive throughout it). The octahedron values assume an
# isotropic 3D zero, which real axial beams are not, so treat them as a guide.
TCP_PATTERN_GEOM <- list(
  hexagon    = list(dim = 2L, a_ideal = sqrt(8),    c_off = 8.6),
  square     = list(dim = 2L, a_ideal = sqrt(8),    c_off = 9.0),
  triangle   = list(dim = 2L, a_ideal = sqrt(8),    c_off = 13.7),
  octahedron = list(dim = 3L, a_ideal = sqrt(16/3), c_off = 9.0)
)
TCP_PATTERN_DEFAULT <- "hexagon"

# ── .tcp_pattern_name() ───────────────────────────────────────────────────────
# GOAL: one canonical pattern name from what a user or a file may write: the name
# itself, a short form, or the 2D vertex count (6 hexagon, 4 square, 3 triangle).
# 6 means the hexagon, since the octahedron also has 6 vertices; 3D needs its name.
# RETURNS a name in TCP_PATTERN_GEOM, "" for blank, or NA when not recognised.
.tcp_pattern_name <- function(x) {
  v <- tolower(trimws(as.character(x)))[1]
  if (is.na(v) || !nzchar(v)) return("")
  aliases <- c("6" = "hexagon", hex = "hexagon", hexagonal = "hexagon",
               "4" = "square", sq = "square",
               "3" = "triangle", tri = "triangle", triangular = "triangle",
               oct = "octahedron", octahedral = "octahedron")
  if (v %in% names(aliases)) v <- unname(aliases[v])
  if (v %in% names(TCP_PATTERN_GEOM)) v else NA_character_
}

# ── .tcp_geom() ───────────────────────────────────────────────────────────────
# GOAL: look up the geometry constants of a TCP pattern (axes spanned, ideal
# prefactor, off-centre coefficient) from TCP_PATTERN_GEOM above.
# Unknown or missing names fall back to the default (hexagon), so a sequence with
# an unusual final pattern still gets a sensible 2D answer.
.tcp_geom <- function(pattern = TCP_PATTERN_DEFAULT) {
  p <- .tcp_pattern_name(pattern %||% TCP_PATTERN_DEFAULT)
  g <- if (is.na(p)) NULL else TCP_PATTERN_GEOM[[p]]
  if (is.null(g)) g <- TCP_PATTERN_GEOM[[TCP_PATTERN_DEFAULT]]
  g
}

# ── .crlb_doughnut() - photon-budget sigma per localization ─────────────────
# GOAL: per-axis localization precision of each localization from its photon
# count, for the pattern the sequence actually ran.
#
#   sigma = L/(a sqrt(N)) * B * sqrt(1 + c <r^2>/L^2)
#
# L        TCP diameter of the reiterated iteration (nm)
# nu       FWHM of the excitation doughnut (nm); 0 switches the beam term off
# N        detected photons for this localization (eco), vectorised
# D, t_loc diffusion coefficient (nm^2/ms) and time to localization (ms). They
#          set the offset of the emitter from the TCP centre during the exposure,
#          <r^2> = dim * D * t_loc (the end-point spread 2*dim*D*t_loc, averaged
#          over the exposure because photons arrive throughout it). Zero for an
#          immobile emitter.
# pref     the effective prefactor a (measured, see packages.R)
# pattern  TCP geometry: hexagon, square, triangle or octahedron. Sets dim and c.
#
# B is the beam-size factor of Vogler et al. 2025 Eq. 1 at the TCP centre,
# 1/|1 - ln2 (L/nu)^2|, which makes sigma grow faster than L as L approaches nu.
# The offset is handled by the pattern's own coefficient rather than by Eq. 1's
# x-dependence, which is written for a single axis.
.crlb_doughnut <- function(L, nu, N, D = 0, t_loc = 0, pref = TCP_CRLB_PREFACTOR,
                           pattern = TCP_PATTERN_DEFAULT) {
  if (!length(N)) return(numeric(0))
  if (!is.finite(pref) || pref <= 0) pref <- TCP_CRLB_PREFACTOR
  base <- L / (pref * sqrt(N))                     # centred emitter, no beam term
  if (!is.finite(L) || L <= 0) return(base)
  g <- .tcp_geom(pattern)
  # mean squared offset over the exposure; motion can only add to it
  r2 <- if (is.finite(D) && D > 0 && is.finite(t_loc) && t_loc > 0) g$dim * D * t_loc else 0
  off <- sqrt(1 + g$c_off * r2 / L^2)
  # beam term, only where the parabolic-zero picture still holds (den well above 0)
  B <- 1
  if (is.finite(nu) && nu > 0) {
    den <- abs(1 - log(2) * (L / nu)^2)
    if (is.finite(den) && den > 0.1) B <- 1 / den
  }
  base * B * off
}
