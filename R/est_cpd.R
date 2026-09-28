## ─────────────────────────────────────────────────────────────────────────────
## est_cpd.R  -  change-point segmentation of tracks by local diffusivity
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The estimator behind the Tracking tab's "Track Segmentation (Change-Point +
##   GMM)" panel, as pure functions. Moved out of tracking_server() on 2026-09-17
##   (publication audit C1); cpd_segments_rv() / cpd_gmm_rv() now only gather
##   inputs and call these.
##
## METHOD
##   Per track, the signal log(r²/Δt) (local log-diffusivity, already Δt-aware) is
##   divided by its fixed SD √trigamma(d/2) and segmented by PELT
##   (Killick, Fearnhead & Eckley 2012; changepoint::cpt.mean, BIC penalty,
##   minseglen). Per segment D = mean(r²/Δt)/(2d) (mean of ratios: no Jensen bias
##   when r² and Δt co-vary), minus the noise floor σ²·mean(1/Δt), clamped at 0,
##   divided by the blur factor 1 − 2R·t_c·mean(1/Δt).
##   The pooled per-segment log10 D is then fitted with a BIC-selected Gaussian
##   mixture (mclust, model "V").
##
## KNOWN LIMITS (measured 2026-09-17, methods doc §I-B.3)
##   - Scaled PELT still cuts 3-9% of single-state 2D tracks.
##   - The mixture treats every segment as equally precise, so segments of very
##     different lengths from ONE population can return two components. The
##     component count is exploratory; the HMM (est_hmm.R) is the tool for
##     "how many states, how fast do they switch".
##
## UNITS: positions nm, tim ms, D nm²/ms (÷1000 → µm²/s).
## ─────────────────────────────────────────────────────────────────────────────


# ── cpd_segment_tracks() ──────────────────────────────────────────────────────
# GOAL: one row per detected segment, for every track long enough to hold two
# segments of `min_seg` steps.
#
# df        : data.frame(tid, x, y, tim[, z]); tim in ms
# min_seg   : minimum steps per segment (PELT minseglen); tracks shorter than
#             2 × min_seg steps are skipped
# use_z     : include z in r² (3D data)
# noise_sig2: σ_loc² (nm²) to subtract as σ²·mean(1/Δt); 0 = no subtraction
# R_abs     : blur offset in ms (blur_offsets()$R_abs, est_blur.R); the segment D is divided by
#             1 − 2·R_abs·mean(1/Δt). 0 = no blur correction. Without it D read
#             about 30% low on blurred simulated tracks (same bias as JDD/D-map).
# progress  : optional function(fraction) called once per track
#
# RETURNS data.frame(tid, seg_idx, n_states, n_steps, D_nm2_ms, D_um2_s), or NULL
# when no track qualifies. A track with any non-positive or non-finite Δt is
# skipped whole (unchanged from the in-server version).
cpd_segment_tracks <- function(df, min_seg = 8L, use_z = FALSE, noise_sig2 = 0,
                               R_abs = 0, progress = NULL) {
  min_seg <- max(3L, as.integer(min_seg))
  ndim    <- if (isTRUE(use_z)) 3L else 2L
  coef_D  <- 2L * ndim                                    # MSD = coef_D * D * dt
  tids    <- split(df, df$tid)
  seg_rows <- lapply(seq_along(tids), function(i) {
    if (is.function(progress)) progress(1 / length(tids))
    td <- tids[[i]][order(tids[[i]]$tim), ]
    n  <- nrow(td)
    if (n < 2L * min_seg) return(NULL)

    dt <- diff(td$tim)                                    # ms
    if (any(!is.finite(dt) | dt <= 0)) return(NULL)
    r2 <- diff(td$x)^2 + diff(td$y)^2                     # nm²
    if (isTRUE(use_z) && "z" %in% names(td) && any(is.finite(td$z)))
      r2 <- r2 + diff(td$z)^2                             # include z for 3D data

    # Step size normalised by time lag -> proportional to local D; log makes the
    # skewed chi-square roughly symmetric and turns multiplicative changes additive.
    r2n     <- pmax(r2 / dt, .Machine$double.eps)
    log_r2n <- log(r2n)

    # PELT on the mean of log(r²/Δt). cpt.mean's Normal cost assumes UNIT variance,
    # so its BIC penalty is only calibrated on a unit-variance signal. log of a
    # scaled chi-square with d degrees of freedom has a fixed variance whatever D
    # and Δt are: trigamma(d/2) = π²/6 ≈ 1.64 in 2D, 0.93 in 3D. Unscaled, the 2D
    # penalty is 1.64× too weak. Simulated single free population (MINFLUX tracking timings,
    # σ 5 nm): unscaled 28% / 35% / 51% of tracks cut at N = 100 / 300 / 1000
    # steps; scaled 8% / 9% / 3%.
    log_r2n_std <- log_r2n / sqrt(trigamma(ndim / 2))
    cpt_obj <- tryCatch(
      changepoint::cpt.mean(log_r2n_std, method = "PELT", penalty = "BIC",
                            minseglen = min_seg),
      error = function(e) NULL)
    bkpts <- if (!is.null(cpt_obj)) changepoint::cpts(cpt_obj) else integer(0)
    seg_s <- c(1L, if (length(bkpts)) bkpts + 1L else integer(0))
    seg_e <- c(bkpts, length(log_r2n))

    tid_val <- as.character(td$tid[1])
    do.call(rbind, Filter(Negate(is.null), lapply(seq_along(seg_s), function(s) {
      idx <- seg_s[s]:seg_e[s]
      if (length(idx) < min_seg) return(NULL)
      # mean(r²/Δt)/coef_D, not mean(r²)/(coef_D·mean(Δt)): the latter is biased
      # when r² and Δt covary (Jensen's inequality).
      D_nm2_ms <- mean(r2[idx] / dt[idx]) / coef_D
      # Noise floor σ²·⟨1/Δt⟩, clamped ≥ 0, so immobile segments read ~0.
      if (noise_sig2 > 0)
        D_nm2_ms <- max(0, D_nm2_ms - noise_sig2 * mean(1 / dt[idx]))
      # Blur shortens the effective diffusion time of each step to Δt − 2·R_abs.
      if (is.finite(R_abs) && R_abs > 0)
        D_nm2_ms <- D_nm2_ms / max(0.05, 1 - 2 * R_abs * mean(1 / dt[idx]))
      data.frame(tid = tid_val, seg_idx = s, n_states = length(seg_s),
                 n_steps = length(idx), D_nm2_ms = D_nm2_ms, D_um2_s = D_nm2_ms / 1e3,
                 stringsAsFactors = FALSE)
    })))
  })
  do.call(rbind, Filter(Negate(is.null), seg_rows))
}

# ── cpd_gmm() ─────────────────────────────────────────────────────────────────
# GOAL: BIC-selected Gaussian mixture over per-segment log10 D (µm²/s).
# Log space because D spans orders of magnitude and is right-skewed. See the
# file header for why the component count is exploratory.
# segs  : output of cpd_segment_tracks()
# max_G : largest component count tried (1-6)
# RETURNS an mclust object, or NULL.
cpd_gmm <- function(segs, max_G = 4L) {
  if (is.null(segs) || nrow(segs) < 3L) return(NULL)
  log10_D <- log10(segs$D_um2_s)
  log10_D <- log10_D[is.finite(log10_D)]
  if (length(log10_D) < 3L) return(NULL)
  max_G <- max(1L, min(6L, as.integer(max_G)))
  tryCatch(mclust::Mclust(log10_D, G = seq_len(max_G), modelNames = "V", verbose = FALSE),
           error = function(e) NULL)
}
