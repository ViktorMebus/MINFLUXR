## ─────────────────────────────────────────────────────────────────────────────
## est_blur.R  -  motion-blur coefficient R of one MINFLUX tracking localization
## ─────────────────────────────────────────────────────────────────────────────
##
## WHAT THIS FILE HOLDS
##   blur_pattern_vectors() .. unit vectors of the TCP vertices, in visiting order
##   blur_R_step() ........... R per step (and per axis) from the sequence timing
##   blur_offsets() .......... R converted to the forms the estimators consume
##
## WHY R IS NEEDED
##   A localization is not a snapshot: its photons are collected while the molecule
##   moves. Every estimator that reads a step's squared displacement then sees
##       E[dx^2] = 2 D (dt - 2 R dt) + 2 sigma^2           (per axis)
##   (Berglund 2010; Vestergaard, Blainey & Flyvbjerg 2014), so D read without R is
##   too low and D read with too large an R is too high.
##
## WHAT SETS R (and what does not)
##   R = (1/dt) * integral over the step of S(t) (1 - S(t)) dt, where S(t) is the
##   cumulative weight the position estimate gives to the molecule's path up to t.
##   So R is set by WHEN the photons that make up the estimate are collected inside
##   the step: the dwell per vertex, the number and order of vertices, pattern
##   repeats, the passes needed to reach the photon limit (eta), the beam moves
##   between vertices (t_EOD ~ 5 us), the CFR centre dwell and the re-centring
##   overhead (t_REP ~ 20 us) (Vogler et al. 2025, Eq. 4-7). The TCP diameter L
##   does not enter: simulated with a MINFLUX estimator at L = 50, 101 and 302 nm
##   and D = 0.2 and 2 um^2/s, R stayed at 0.084-0.095 (private audit, 2026-09-28).
##
## THE MINFLUX-SPECIFIC PART
##   A vertex only senses the displacement along its own direction, so the linear
##   MINFLUX estimator weights vertex i by u_ia^2 on axis a (hexagon, x: 1, 1/4, 1/4,
##   1, 1/4, 1/4; y: 0, 3/4, 3/4, 0, 3/4, 3/4). R therefore differs per axis and
##   depends on the order the vertices are visited, which is taken as angular order.
##   For the GP67H sequence (hexagon, 100 us dwell, CFR on, 203 us step) this gives
##   x 0.103 / y 0.097 per step, against 0.096 / 0.097 from a full simulation of
##   photon counts and the estimator.
##
## WHY NOT THE CONSTANT 1/6.2
##   MinSpt (Vogler et al.) passes R = 1/6.2 of the cycle, cited to Schmidt et al.
##   2020 with a question mark; that paper gives no blur coefficient, and no MINFLUX
##   sequence field sets one. On the GP67H timing it is 0.152 per step, ~50% above
##   the value the timing gives, and made single-step D read 15-17% high on
##   simulated undamped tracks. The timing model reads within 2-3%.
##
## WHAT IT CANNOT DO
##   Damping (sequence `damping` > 0) smooths the exported positions ACROSS steps.
##   That is not a within-step blur and no R describes it; the caller flags it.
## ─────────────────────────────────────────────────────────────────────────────

# Fallback R per step, used only when neither a sequence file nor eco/efo gives the
# timing: the value blur_R_step() returns for the standard 2D tracking sequence
# (hexagon, 100 us dwell, one pass, ~203 us step).
BLUR_R_STEP_FALLBACK <- 0.10

# =============================================================================
# blur_pattern_vectors() - TCP VERTEX DIRECTIONS, IN VISITING ORDER
# =============================================================================
# GOAL: the unit vector of each vertex, so each can be weighted by the axis it
#   senses. Patterns follow Imspector's names; the order is angular.
# ARGS  pattern : "hexagon", "square", "triangle", "octahedron", or NA
#       n_vert  : vertex count used when the pattern is unknown (regular polygon)
#   "zline*" patterns are axial (vertices on the optical axis).
# RETURNS a K x 3 matrix (x, y, z) of unit vectors.
blur_pattern_vectors <- function(pattern = NA_character_, n_vert = 6L) {
  polygon <- function(k) { a <- (seq_len(k) - 1) * 2 * pi / k; cbind(cos(a), sin(a), 0) }
  p <- tolower(as.character(pattern %||% NA_character_))
  if (identical(p, "hexagon"))  return(polygon(6L))
  if (identical(p, "square"))   return(polygon(4L))
  if (identical(p, "triangle")) return(polygon(3L))
  if (identical(p, "octahedron"))
    return(rbind(c(1, 0, 0), c(0, 1, 0), c(-1, 0, 0), c(0, -1, 0), c(0, 0, 1), c(0, 0, -1)))
  # Axial-only patterns: the vertices lie on the optical axis, so they carry no
  # lateral information and give no lateral R (callers fall back to the data).
  if (!is.na(p) && grepl("^zline", p)) return(rbind(c(0, 0, 1), c(0, 0, -1)))
  k <- suppressWarnings(as.integer(n_vert)); if (!is.finite(k) || k < 2L) k <- 6L
  polygon(k)
}

# =============================================================================
# blur_R_step() - BLUR COEFFICIENT PER STEP FROM THE SEQUENCE TIMING
# =============================================================================
# GOAL: R for one localization, per axis, from how its photons are collected.
#
# THE TIMELINE (one step of length dt_ms, Vogler et al. 2025 Eq. 4)
#   repeated eta times:  for each pattern repeat, for each vertex: a beam move of
#                        t_eod_ms, then a dwell of t_dwell_ms / (K * pat_repeat);
#                        then, if CFR is on, a move and a centre dwell of
#                        ctr_factor * t_dwell_ms (not used for the position)
#   then t_rep_ms for estimation and re-centring, then idle until dt_ms.
#   Where the idle time sits does not change R (S is flat there).
#
# ARGS
#   dt_ms       : mean interval between consecutive single-cycle localizations
#   t_dwell_ms  : pattern dwell time per pass (sequence patDwellTime); when the
#                 sequence is not known, the photon collection time eco/efo
#   n_vert, pattern, pat_repeat, eta : pattern geometry and passes per localization
#   cfr, ctr_factor : CFR centre check on/off and its dwell factor (ctrDwellFactor)
#   t_eod_ms, t_rep_ms : hardware times (Vogler et al. 2025: ~5 us and ~20 us)
#
# RETURNS list(R = mean over x and y, R_axis = c(x, y, z), R_uniform (all vertices
#   weighted equally, the older timing model), t_loc_ms (modelled time to
#   localization), overrun = TRUE when the modelled timeline is longer than dt_ms
#   (then compressed to fit), K), or NULL when the inputs are unusable.
blur_R_step <- function(dt_ms, t_dwell_ms, n_vert = 6L, pattern = NA_character_,
                        pat_repeat = 1L, eta = 1L, cfr = FALSE, ctr_factor = 0,
                        t_eod_ms = 0.005, t_rep_ms = 0.020) {
  if (!is.finite(dt_ms) || dt_ms <= 0 || !is.finite(t_dwell_ms) || t_dwell_ms <= 0) return(NULL)
  U  <- blur_pattern_vectors(pattern, n_vert); K <- nrow(U)
  nr <- max(1L, as.integer(pat_repeat %||% 1L)); if (!is.finite(nr)) nr <- 1L
  ne <- max(1L, as.integer(round(eta %||% 1L))); if (!is.finite(ne)) ne <- 1L
  cf <- if (isTRUE(cfr) && is.finite(ctr_factor) && ctr_factor > 0) ctr_factor else 0
  dw <- t_dwell_ms / (K * nr)                       # dwell per vertex visit

  # Build the collection windows [start, end] and the vertex index of each.
  st <- numeric(0); en <- numeric(0); vid <- integer(0); t <- 0
  for (p in seq_len(ne)) {
    for (r in seq_len(nr)) for (i in seq_len(K)) {
      t <- t + t_eod_ms; st <- c(st, t); en <- c(en, t + dw); vid <- c(vid, i); t <- t + dw
    }
    if (cf > 0) t <- t + t_eod_ms + cf * t_dwell_ms   # centre dwell: time passes, no weight
  }
  t_loc <- t + t_rep_ms
  overrun <- t_loc > dt_ms * 1.02
  if (overrun) { s <- dt_ms / t_loc; st <- st * s; en <- en * s }   # keep the shape, fit the step

  # Per-axis weight of each visit: the share of that axis' sensitivity it carries.
  R_for <- function(w_vert) {
    w <- w_vert[vid]; if (sum(w) <= 0) return(NA_real_)
    w <- w / sum(w)
    # S rises linearly across each window by its weight; exact integral of S(1-S)
    # piecewise: flat parts contribute S0(1-S0)*gap, ramps contribute the integral
    # of (S0 + w u)(1 - S0 - w u) du over the window.
    o <- order(st); st2 <- st[o]; en2 <- en[o]; w2 <- w[o]
    S <- 0; acc <- 0; tprev <- 0
    for (j in seq_along(st2)) {
      acc <- acc + S * (1 - S) * (st2[j] - tprev)                       # flat before the window
      len <- en2[j] - st2[j]; a <- w2[j]
      acc <- acc + len * (S * (1 - S) + a * (1 - 2 * S) / 2 - a^2 / 3)  # ramp
      S <- S + a; tprev <- en2[j]
    }
    acc / dt_ms                                     # S = 1 after the last window: no more area
  }
  sens <- U^2
  R_axis <- c(x = R_for(sens[, 1]), y = R_for(sens[, 2]),
              z = if (any(sens[, 3] > 0)) R_for(sens[, 3]) else NA_real_)
  list(R = mean(R_axis[c("x", "y")], na.rm = TRUE), R_axis = R_axis,
       R_uniform = R_for(rep(1, K)), t_loc_ms = t_loc, overrun = overrun, K = K)
}

# =============================================================================
# blur_offsets() - R IN THE FORMS THE ESTIMATORS USE
# =============================================================================
# GOAL: one conversion, so every estimator gets the same blur.
#   R_step : per step, relative to the single-cycle interval dt_ms
#   R_abs  : blur offset in ms, R_step * dt_ms; each step's effective diffusion
#            time is dt - 2 R_abs (JDD, D-map, change-point, HMM)
#   R_msd  : relative to the cycle time, R_abs / t_c (MSD / OLSF, confinement null)
# ARGS  R_step, dt_ms (single-cycle interval), t_c (cycle time), both ms
blur_offsets <- function(R_step, dt_ms, t_c) {
  ok <- is.finite(R_step) && R_step >= 0 && is.finite(dt_ms) && dt_ms > 0
  R_abs <- if (ok) R_step * dt_ms else 0
  list(R_step = if (ok) R_step else 0, R_abs = R_abs,
       R_msd = if (ok && is.finite(t_c) && t_c > 0) R_abs / t_c else 0)
}
