## ─────────────────────────────────────────────────────────────────────────────
## sim_server.R  -  Simulation tab server logic
## Synthetic MINFLUX localizations (shapes OR image-guided)
## ─────────────────────────────────────────────────────────────────────────────

sim_server <- function(input, output, session) {

  # ── Zoom-aware scale bar (mirrors the SMLM tab) ──────────────────────────────
  # The bar lives in paper coordinates with width = length_nm / visible_x_span,
  # so it must be re-pushed whenever the user zooms (the visible span changes).
  sim_axis_x_span <- reactiveVal(NULL)

  # TRUE once each plot has actually rendered WITH event_register(), so the
  # relayout observers never call event_data() before the source is registered.
  # That call would otherwise print a benign "source not registered" warning,
  # deferred via session$onFlushed() so suppressWarnings() can't catch it. They
  # stay FALSE until the user clicks Generate (the render short-circuits on a
  # validate() until then, so the source is never registered before that).
  sim_plot_ready <- reactiveVal(FALSE)
  sim_diff_ready <- reactiveVal(FALSE)

  # Single source of truth: is the FOV / ROI window active?  TRUE in Acquisition
  # mode, or in Normal mode when the user pins an ROI ("Pin axes to a FOV").
  # Used everywhere the FOV matters - axis pinning, background spread, and the
  # confocal window - so they can never disagree.
  sim_use_fov <- function() {
    identical(input$sim_mode %||% "normal", "acquisition") ||
      isTRUE(input$sim_view_fov)
  }

  # Full-view x-span: the FOV X when the FOV is active, else the data extent.
  sim_full_x_span <- function() {
    if (isolate(sim_use_fov()))
      return(max(1, as.numeric(isolate(input$sim_fov_x) %||% 4000)))
    res <- tryCatch(isolate(sim_result()), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    ref  <- if (nrow(res$locs) > 0) res$locs else res$gt
    span <- suppressWarnings(diff(range(ref$x, na.rm = TRUE)))
    if (length(span) == 1 && is.finite(span) && span > 0) span else NULL
  }

  sim_push_scale_bar <- function() {
    proxy <- plotlyProxy("sim_plot", session)
    if (!isTRUE(input$sim_show_scale_bar) ||
        !identical(input$sim_view_mode, "2d")) {
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = list(), annotations = list()))
      return()
    }
    span <- sim_axis_x_span()
    if (length(span) != 1 || !is.finite(span) || span <= 0) {
      span <- sim_full_x_span()
      if (length(span) == 1 && is.finite(span) && span > 0) sim_axis_x_span(span)
    }
    if (length(span) != 1 || !is.finite(span) || span <= 0) return()
    sb_font <- suppressWarnings(as.numeric(input$sim_scale_font_size))
    sb <- make_scale_bar(
      suppressWarnings(as.numeric(input$sim_scale_bar_length)),
      span,
      suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
      dark = isTRUE(input$sim_dark_mode),
      scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
    if (!is.null(sb))
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = sb$shapes, annotations = sb$annotations))
  }

  # Same zoom-aware scale bar for the Confocal view (always white on black).
  sim_diff_x_span <- reactiveVal(NULL)   # zoomed visible x-span
  sim_diff_full   <- reactiveVal(NULL)   # full rendered-image x-span

  sim_push_diff_scale_bar <- function() {
    proxy <- plotlyProxy("sim_diff_plot", session)
    if (!isTRUE(input$sim_show_scale_bar)) {
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = list(), annotations = list()))
      return()
    }
    span <- sim_diff_x_span()
    if (length(span) != 1 || !is.finite(span) || span <= 0) span <- sim_diff_full()
    if (length(span) != 1 || !is.finite(span) || span <= 0) return()
    sb_font <- suppressWarnings(as.numeric(input$sim_scale_font_size))
    sb <- make_scale_bar(
      suppressWarnings(as.numeric(input$sim_scale_bar_length)),
      span,
      suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
      dark = TRUE,
      scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
    if (!is.null(sb))
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = sb$shapes, annotations = sb$annotations))
  }

  # ── Detection-rate SHAPE models ─────────────────────────────────────────────
  # The instantaneous structure detection rate is  struct_rate · g(t),  where
  # g(t) ∈ [0, 1] is a normalized SHAPE and struct_rate is the CHARACTERISTIC
  # (reference) rate - the peak for a decaying run, the plateau for a rising or
  # constant one. Three shapes cover the profiles seen in real MINFLUX runs:
  #
  #   "decay"    g(t) = level + (1-level)·e^{-t/τ}     starts at 1, → level
  #              Blinking dyes: bright early, bleach toward a floor. `level` =
  #              late/peak ratio (0 = decays to zero, 1 = no decay).
  #   "constant" g(t) = 1                              flat (variation is Poisson)
  #   "rise"     g(t) = 1 − (1-level)·e^{-t/τ}         starts at `level`, → 1
  #              Activation-limited: almost no detection early, then rises to a
  #              cap. `level` = initial/plateau ratio (0 = starts from zero).
  #
  # `t` and `tau` must share units (both seconds, or both minutes). level clamped
  # to [0,1]. The decay case is identical to the old hard-coded formula, so
  # existing calibrations are unchanged.
  acq_shape <- function(t, model, tau, level) {
    tau <- max(1e-6, tau); lv <- max(0, min(1, level))
    e   <- exp(-t / tau)
    switch(model %||% "decay",
      constant = rep(1, length(t)),
      rise     = 1 - (1 - lv) * e,
      lv + (1 - lv) * e)                 # decay (default)
  }
  # Closed-form ∫₀ᵀ g(t) dt  (T, tau in the SAME units). Used for expected counts.
  acq_shape_integral <- function(T, model, tau, level) {
    tau <- max(1e-6, tau); lv <- max(0, min(1, level))
    ee  <- tau * (1 - exp(-T / tau))     # ∫₀ᵀ e^{-t/τ} dt
    switch(model %||% "decay",
      constant = T,
      rise     = T - (1 - lv) * ee,
      lv * T + (1 - lv) * ee)            # decay (default)
  }

  # Sample n detection times on [0, T] from the (front- or back-loaded) rate
  # shape, by inverse-CDF sampling of the density g(t). model/level select the
  # profile (see acq_shape); tau_sec is the rate-change time constant.
  sim_sample_detection_times <- function(n, T_sec, tau_sec, floor_frac,
                                         model = "decay") {
    if (n <= 0) return(numeric(0))
    if (!is.finite(T_sec) || T_sec <= 0) return(rep(0, n))
    ng   <- 2000L
    tg   <- seq(0, T_sec, length.out = ng)
    dens <- pmax(1e-9, acq_shape(tg, model, tau_sec, floor_frac))
    cdf  <- cumsum(dens); cdf <- cdf / cdf[ng]
    approx(cdf, tg, runif(n), rule = 2, ties = "ordered")$y
  }

  # Expected number of structure detections over [0, T] for an ABSOLUTE rate
  # r(t) = struct_rate · g(t)  [traces/min]:  N_T = struct_rate · ∫₀ᵀ g dt.
  # (T, τ in min.) For decay this is the old floor·T + (1-floor)·τ·(1-e^{-T/τ}).
  acq_expected_count <- function(struct_rate, T_sec, decay_min, floor_frac,
                                 model = "decay") {
    T_min <- T_sec / 60
    struct_rate * acq_shape_integral(T_min, model, max(1e-6, decay_min), floor_frac)
  }

  # Acquisition timeline from a fixed fluorophore pool + live acquisition
  # parameters.  Operates only on molecule IDs (no positions/localizations), so
  # it is cheap enough to recompute live as rate / duration / background change
  # without re-running the whole structure.  pool_mol = mol_id of every unique
  # fluorophore (length = pool size); n_mol_total = molecules placed (for
  # display).  Returns the same shape of list the acquisition plot/summary use.
  compute_acq_timeline <- function(pool_mol, n_mol_total, struct_rate, T_sec,
                                   blinks, decay_min, floor_frac, bg_rate,
                                   resolve_q, crowd_m = NULL, model = "decay") {
    n_pool <- length(pool_mol)
    if (n_pool == 0 || !is.finite(T_sec) || T_sec <= 0 || struct_rate <= 0)
      return(NULL)
    n_mol_resolvable <- length(unique(pool_mol))

    # Deterministic EXPECTED timeline (no per-realization jitter), so the curves
    # and the time-to-resolve are stable as the user drags rate / duration.
    # struct_rate / bg_rate are TRACES PER HOUR (easier to match real settings);
    # converted to per-minute internally.  Detections follow an ABSOLUTE rate
    # r(t) = X·(floor + (1-floor)·e^{-t/τ}), so how fast the structure resolves
    # depends on the rate and structure - NOT on the chosen duration.
    T_min     <- T_sec / 60
    tau_min   <- max(1e-6, decay_min)
    rate_min  <- struct_rate / 60                    # per hour → per minute
    bg_min    <- bg_rate / 60

    # Photon budget: each fluorophore yields ~blinks detections before bleaching.
    # This caps the TOTAL TIDs (over-counting) only - it does NOT gate resolution
    # (a molecule is resolved by its FIRST detection; the extra blinks just pile
    # more TIDs on top).
    struct_cap <- n_pool * max(1, blinks)
    f_m        <- as.numeric(tapply(rep(1, n_pool), pool_mol, sum))  # dyes / molecule

    ng     <- 240L
    tg_min <- seq(0, T_min, length.out = ng)
    tg_sec <- tg_min * 60
    # Cumulative fluorophore-detection attempts from the absolute rate:
    # ∫₀ᵗ g(s) ds for the chosen shape (decay / constant / rise).
    cum_R <- rate_min * vapply(tg_min,
               function(tt) acq_shape_integral(tt, model, tau_min, floor_frac),
               numeric(1))
    # Photobleaching consumes the photon budget GRADUALLY (fluorophores bleach at
    # different times), so the total saturates smoothly toward struct_cap and the
    # structure rate tapers off - no hard cliff when the budget is reached.
    bleach <- exp(-cum_R / struct_cap)          # fraction of budget still emitting
    N_t    <- struct_cap * (1 - bleach)         # total recorded TIDs (over-counting)
    n_s    <- N_t[ng]

    # Expected distinct molecules resolved by time t: a molecule is found once
    # ANY of its f_m fluorophores is first detected.  Driven by cum_R (rate ×
    # time × pool) - INDEPENDENT of the per-fluorophore blink budget, so
    # "detections before bleaching" never shifts the road to X% resolved.
    # crowd_m (≤ 1, per molecule, aligned to the sorted unique mol_ids) slows
    # discovery where emitters crowd the donut: a crowded molecule's effective
    # detection rate is scaled down, so dense structures resolve slower.
    cw <- if (is.null(crowd_m) || length(crowd_m) != length(f_m))
            rep(1, length(f_m)) else crowd_m
    cum_unique <- vapply(cum_R,
                         function(c) sum(1 - exp(-f_m * cw * c / n_pool)),
                         numeric(1))
    n_mol_det  <- cum_unique[ng]
    denom      <- max(1, n_mol_resolvable)

    # Instantaneous detection rate (traces/HOUR): structure (the chosen shape g(t)
    # × the photobleaching taper) + constant background floor.
    rate_struct <- struct_rate * acq_shape(tg_min, model, tau_min, floor_frac) * bleach
    rate_v <- rate_struct + bg_rate
    n_bg   <- max(0L, min(200000L, as.integer(round(bg_min * T_min))))

    q        <- max(0.5, min(1, resolve_q))
    res_time <- if (n_mol_det >= q * denom)
                  approx(cum_unique, tg_sec, q * denom, ties = "ordered")$y
                else NA_real_

    list(
      T_sec = T_sec, struct_rate = struct_rate, bg_rate = bg_rate,
      n_pool = n_pool, n_mol_total = n_mol_total,
      n_mol_resolvable = n_mol_resolvable, struct_cap = struct_cap, n_bg = n_bg,
      t = tg_sec, cum_unique = cum_unique, rate_t = tg_sec, rate_v = rate_v,
      rate_struct = rate_struct,          # structure-only component of rate_v
      rate_peak  = max(rate_v),
      rate_floor = stats::median(tail(rate_v, max(1L, length(rate_v) %/% 4L))),
      rate_mean  = n_s / max(1e-9, T_min) * 60,      # traces/hour
      n_mol_det = round(n_mol_det), resolved_frac = n_mol_det / denom,
      n_real_tid = round(n_s), det_per_mol = n_s / max(1, n_mol_det),
      resolve_q = q, res_time = res_time,
      overcount_frac = 1 - n_mol_det / max(1, n_s))
  }

  # Calibrate the acquisition model from a REAL MINFLUX run: measure the
  # detection-rate profile (peak / mean / floor / decay τ), the run duration and
  # the FOV extent from the data.  Works on the first detection time of each
  # unique trace (tid).  df = canonical MINFLUX frame (tid, tim [s], x, y[, z] nm).
  analyze_real_run <- function(df) {
    if (!all(c("tid", "tim", "x", "y") %in% names(df))) return(NULL)
    df <- df[is.finite(df$tim) & is.finite(df$x) & is.finite(df$y), , drop = FALSE]
    if (nrow(df) < 10) return(NULL)

    # First detection time + median position of each trace (same tid order).
    t0_tid <- as.numeric(tapply(df$tim, df$tid, min))
    tx     <- as.numeric(tapply(df$x,   df$tid, median))
    ty     <- as.numeric(tapply(df$y,   df$tid, median))
    t0_tid <- t0_tid - min(t0_tid)                    # seconds from start
    n_tr   <- length(t0_tid)
    T_sec  <- max(t0_tid)
    if (!is.finite(T_sec) || T_sec <= 0) return(NULL)

    # Binned detection rate (traces / hour).
    nb     <- max(8L, min(60L, n_tr %/% 10L))
    br     <- seq(0, T_sec, length.out = nb + 1L)
    bin_hr <- max(1e-9, (br[2] - br[1]) / 3600)
    h      <- hist(t0_tid, breaks = br, plot = FALSE)
    rate_t <- h$mids                                  # seconds
    rate_v <- h$counts / bin_hr

    R0     <- mean(head(rate_v, max(1L, nb %/% 10L)))           # early rate
    Rf     <- stats::median(tail(rate_v, max(1L, nb %/% 4L)))   # late rate

    R_mean <- n_tr / (T_sec / 3600)

    # τ from a log-linear fit of |rate − late_rate| toward 0 - works whether the
    # rate falls (dye bleaching) OR rises (e.g. buffer activating over the run).
    tau_min <- 15
    exc     <- abs(rate_v - Rf)
    ok      <- exc > 0.05 * max(exc, 1e-9) & rate_t < T_sec * 0.8
    if (sum(ok) >= 3) {
      fit <- tryCatch(stats::lm(log(exc[ok]) ~ rate_t[ok]), error = function(e) NULL)
      if (!is.null(fit)) {
        slope <- stats::coef(fit)[[2]]
        if (is.finite(slope) && slope < 0) tau_min <- max(0.5, (-1 / slope) / 60)
      }
    }
    # If the rate visibly changed over the run, τ can't exceed it - a noisy
    # (especially rising) profile can fit an absurd τ ≫ run, which would make the
    # simulated change look flat.  Cap at the run length.
    tau_min <- min(tau_min, T_sec / 60)

    # ── Choose the detection-rate MODEL that best fits the profile ─────────────
    # All three are the same family  A + B·e^{-t/τ}:  decay B>0 (high→low),
    # rise B<0 (low→high), constant B≈0. Compare the exponential fit's residual
    # to a flat line; if it barely helps, or the endpoints hardly differ, call it
    # constant. Otherwise the SIGN of (early − late) picks decay vs rise, and the
    # level/reference are set so struct_rate is the peak (decay) or plateau (rise).
    okv       <- is.finite(rate_v)
    rt        <- rate_t[okv]; rv <- rate_v[okv]
    rss_const <- sum((rv - mean(rv))^2)
    pred_exp  <- Rf + (R0 - Rf) * exp(-rt / (tau_min * 60))
    rss_exp   <- sum((rv - pred_exp)^2)
    rel_chg   <- abs(R0 - Rf) / max(1e-9, R_mean)
    # Only fall back to "constant" when the exponential barely helps (< 5% lower
    # residual than a flat line) AND the endpoints differ by < 10% of the mean.
    # (Loosened from 0.9 / 0.15 so a shallow but real decay is not called flat.)
    if (rss_exp >= 0.95 * rss_const || rel_chg < 0.10) {
      det_model <- "constant"; level_frac <- 1;                       ref_rate <- R_mean
    } else if (R0 >= Rf) {
      det_model <- "decay";    level_frac <- max(0, min(1, Rf / max(1e-9, R0))); ref_rate <- R0
    } else {
      det_model <- "rise";     level_frac <- max(0, min(1, R0 / max(1e-9, Rf))); ref_rate <- Rf
    }

    fov_x <- diff(range(df$x, na.rm = TRUE))
    fov_y <- diff(range(df$y, na.rm = TRUE))
    fov_z <- if ("z" %in% names(df)) diff(range(df$z, na.rm = TRUE)) else 0

    # Spatial split: a trace sitting in a locally dense neighbourhood is ON the
    # structure; an isolated trace is background.  Uses the DBSCAN core-point
    # criterion (≥ minPts neighbours within eps) - absolute density, so it works
    # even when the structure itself is densely labeled.  Approximate; the user
    # can override the resulting fraction.  Also bins the two classes separately
    # so the measured structure / background rate curves can be overlaid.
    struct_frac   <- 0.5
    rate_v_struct <- NULL
    rate_v_bg     <- NULL
    if (n_tr >= 20) {
      eps    <- 60          # nm neighbourhood radius
      minPts <- 3L          # core-point threshold
      K      <- min(n_tr, 64L)
      nn     <- RANN::nn2(cbind(tx, ty), cbind(tx, ty), k = K,
                          searchtype = "radius", radius = eps)
      cnt    <- rowSums(nn$nn.idx > 0) - 1L          # neighbours within eps
      is_struct   <- cnt >= minPts
      struct_frac <- max(0, min(1, mean(is_struct)))
      rate_v_struct <- hist(t0_tid[is_struct],  breaks = br, plot = FALSE)$counts / bin_hr
      rate_v_bg     <- hist(t0_tid[!is_struct], breaks = br, plot = FALSE)$counts / bin_hr
    }

    # ── Point-cloud profile ───────────────────────────────────────────────────
    # Localization precision = SD of the repeated final localizations WITHIN a
    # trace (the molecule re-localized while ON), per axis, median over traces
    # with >= 3 localizations.  Photons/loc from the effective count (eco).
    # Effective TCP diameter reproducing this precision: L_eff = 2·σ·√N.
    #
    # CAVEAT (stated in the UI too): a within-trace SD is not pure localization
    # precision — it also contains any real molecular motion, stage drift and
    # sample wobble during the trace.  On a static target that is small; on
    # anything mobile this is an UPPER BOUND on precision, not an estimate of it.
    nl  <- as.numeric(tapply(df$x, df$tid, length))
    sx  <- as.numeric(tapply(df$x, df$tid, stats::sd))
    sy  <- as.numeric(tapply(df$y, df$tid, stats::sd))
    sz  <- if ("z" %in% names(df))
             as.numeric(tapply(df$z, df$tid, stats::sd))
           else rep(NA_real_, length(nl))
    okp <- nl >= 3 & is.finite(sx) & is.finite(sy)

    # Small-n bias correction.  The sample SD under-estimates sigma at small n,
    # which matters now that the user controls the trace-length gate and can let
    # short traces dominate.  Each trace's s is divided by its OWN factor.
    #
    # NOT c4.  The textbook constant c4(n) = E[s]/sigma un-biases the MEAN of the
    # sampling distribution, but this aggregates with a MEDIAN, and median(s)
    # sits BELOW E[s] at small n.  Dividing by c4 under a median therefore leaves
    # about 2.5 % of the bias in place (measured, short-trace-heavy population).
    # The matching constant for a median aggregator is the median of the same
    # chi distribution:  m(n) = sqrt(qchisq(0.5, n-1) / (n-1)).
    #   n = 3:  c4 = 0.886  but  m = 0.833   (raw median is 7.6 % low)
    # Using m recovers sigma to within ~0.2 % of truth.
    # NOTE: either constant CANCELS exactly in the axis ratios below (both axes
    # share the same n), so the ratios do not depend on this choice.
    sd_median_debias <- function(n) {
      n <- pmax(2, n)
      sqrt(stats::qchisq(0.5, n - 1) / (n - 1))
    }
    cc   <- sd_median_debias(nl)
    sx_c <- sx / cc; sy_c <- sy / cc; sz_c <- sz / cc

    # Median and spread (nMAD) of the per-trace sigma.  The SPREAD is the point:
    # real MINFLUX sigma is heteroscedastic (it tracks photon count), and a
    # single median throws exactly that away.
    med_mad <- function(v) {
      v <- v[is.finite(v)]
      if (!length(v)) return(c(NA_real_, NA_real_))
      c(stats::median(v), 1.4826 * stats::mad(v, constant = 1))
    }
    mx <- med_mad(sx_c[okp]); my <- med_mad(sy_c[okp])
    mz <- med_mad(sz_c[okp & is.finite(sz_c)])
    sigma_x <- mx[1]; sigma_y <- my[1]; sigma_z <- mz[1]

    # Axis ratios, as the MEDIAN OF PER-TRACE RATIOS rather than the ratio of
    # medians: within one trace both axes share the same photon count, so the
    # per-trace ratio already divides the photon dependence out and measures
    # cloud SHAPE alone.  This is what an anisotropic simulated cloud needs.
    okr     <- okp & is.finite(sx) & sx > 0
    ratio_yx <- if (any(okr)) stats::median((sy / sx)[okr]) else NA_real_
    okrz     <- okr & is.finite(sz)
    ratio_zx <- if (any(okrz)) stats::median((sz / sx)[okrz]) else NA_real_

    # Distributions kept whole (not just their medians) — 8.1c resamples these
    # to give every simulated localization its own photon count and every
    # simulated trace its own length.  Capped so a huge run cannot bloat the
    # stored profile; resampling from 200 k is indistinguishable from 5 M.
    CAP  <- 200000L
    subs <- function(v) {
      v <- v[is.finite(v)]
      if (length(v) > CAP) v[sample.int(length(v), CAP)] else v
    }
    eco_vec <- if ("eco" %in% names(df)) subs(as.numeric(df$eco)) else numeric(0)
    len_vec <- as.integer(nl[is.finite(nl)])       # ALL gated traces, not just okp

    photons  <- if ("eco" %in% names(df)) stats::median(df$eco, na.rm = TRUE) else NA_real_
    sigma_xy <- mean(c(sigma_x, sigma_y), na.rm = TRUE)
    L_eff    <- if (is.finite(sigma_xy) && is.finite(photons) && photons > 0)
                  2 * sigma_xy * sqrt(photons) else NA_real_

    list(n_traces = n_tr, duration_hr = T_sec / 3600,
         rate_peak = R0, rate_mean = R_mean, rate_floor = Rf,
         det_model = det_model, level_frac = level_frac, ref_rate = ref_rate,
         tau_min = tau_min, struct_frac = struct_frac,
         fov_x = fov_x, fov_y = fov_y, fov_z = fov_z,
         rate_t = rate_t, rate_v = rate_v,
         rate_v_struct = rate_v_struct, rate_v_bg = rate_v_bg,
         sigma_x = sigma_x, sigma_y = sigma_y, sigma_z = sigma_z,
         sigma_x_mad = mx[2], sigma_y_mad = my[2], sigma_z_mad = mz[2],
         ratio_yx = ratio_yx, ratio_zx = ratio_zx,
         n_sigma_traces = sum(okp),
         debias_med = if (any(okp)) stats::median(cc[okp]) else NA_real_,
         eco_vec = eco_vec, eco_n = if ("eco" %in% names(df)) sum(is.finite(df$eco)) else 0L,
         len_vec = len_vec,
         len_med = if (length(len_vec)) stats::median(len_vec) else NA_real_,
         len_q = if (length(len_vec))
                   stats::quantile(len_vec, c(0.25, 0.75), names = FALSE)
                 else c(NA_real_, NA_real_),
         photons = photons, L_eff = L_eff)
  }

  # ╔══════════════════════════════════════════════════════════════════════════╗
  # ║  PHASE 1 - Synthetic localizations                                      ║
  # ╚══════════════════════════════════════════════════════════════════════════╝

  # ── Shape generators ────────────────────────────────────────────────────────

  gen_cluster <- function(n, radius, distribution = "uniform", use_3d = FALSE) {
    sigma <- radius / 2
    if (use_3d) {
      if (distribution == "gaussian") {
        df <- data.frame(x = rnorm(n, 0, sigma),
                         y = rnorm(n, 0, sigma),
                         z = rnorm(n, 0, sigma))
      } else {
        dirs <- matrix(rnorm(n * 3), ncol = 3)
        lens <- sqrt(rowSums(dirs^2)); lens[lens == 0] <- 1
        dirs <- dirs / lens
        r3   <- radius * runif(n)^(1/3)
        df <- data.frame(x = r3 * dirs[,1], y = r3 * dirs[,2], z = r3 * dirs[,3])
      }
    } else {
      r     <- if (distribution == "gaussian") abs(rnorm(n, 0, sigma))
               else radius * sqrt(runif(n))
      theta <- runif(n, 0, 2 * pi)
      df <- data.frame(x = r * cos(theta), y = r * sin(theta))
    }
    df$ldir_type <- NA_character_
    df
  }

  gen_ring <- function(n, diameter, ring_width, z_thickness = 0,
                        use_3d = FALSE) {
    radius <- diameter / 2
    theta  <- runif(n, 0, 2 * pi)
    r      <- rnorm(n, mean = radius, sd = max(0.1, ring_width / 4))
    df     <- data.frame(x = r * cos(theta), y = r * sin(theta))
    if (use_3d) df$z <- rnorm(n, 0, max(0.1, z_thickness / 4))
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  gen_ellipse <- function(n, semi_a, semi_b, semi_c = 1,
                           filled = TRUE, use_3d = FALSE) {
    if (use_3d) {
      phi   <- acos(2 * runif(n) - 1)
      theta <- runif(n, 0, 2 * pi)
      r     <- if (filled) runif(n)^(1/3) else rep(1, n)
      df <- data.frame(x = semi_a * r * sin(phi) * cos(theta),
                       y = semi_b * r * sin(phi) * sin(theta),
                       z = semi_c * r * cos(phi))
    } else {
      theta <- runif(n, 0, 2 * pi)
      if (filled) {
        r  <- sqrt(runif(n))
        df <- data.frame(x = semi_a * r * cos(theta), y = semi_b * r * sin(theta))
      } else {
        df <- data.frame(x = semi_a * cos(theta), y = semi_b * sin(theta))
      }
    }
    df$ldir_type <- NA_character_
    df
  }

  # n molecules spaced at spacing_nm along a straight filament.
  # jitter_lat_nm: noise perpendicular to the filament axis (isotropic in perp plane).
  # jitter_ax_nm:  noise along the filament axis (arc-length perturbation).
  # helix_r / helix_pitch: optional helical arrangement (0 = disabled);
  #   default pitch 36 nm mirrors B-DNA major-groove repeat.
  gen_filament <- function(n, spacing_nm, angle_deg = 0,
                            jitter_lat_nm = 0, jitter_ax_nm = 0,
                            helix_r = 0, helix_pitch = 36,
                            use_3d = FALSE, z_angle_deg = 0) {
    n      <- max(2L, n)
    s      <- (seq_len(n) - (n + 1) / 2) * spacing_nm   # centred arc-lengths
    if (jitter_ax_nm > 0) s <- s + rnorm(n, 0, jitter_ax_nm)

    xy_rad <- angle_deg  * pi / 180
    z_rad  <- z_angle_deg * pi / 180

    # Filament axis unit vector
    ax <- cos(xy_rad) * cos(z_rad)
    ay <- sin(xy_rad) * cos(z_rad)
    az <- sin(z_rad)

    # Two perpendicular unit vectors (analytical, both already unit length):
    # p1 = lateral perpendicular in XY plane; p2 = cross(axis, p1) → toward Z
    p1x <- -sin(xy_rad); p1y <- cos(xy_rad); p1z <- 0
    p2x <- -sin(z_rad) * cos(xy_rad)
    p2y <- -sin(z_rad) * sin(xy_rad)
    p2z <-  cos(z_rad)

    # Lateral jitter as a random walk along the filament: each molecule's
    # perpendicular offset = previous offset + small increment.  This makes
    # the filament curve and meander rather than scatter independently.
    # jitter_lat_nm is the per-step SD (nm); larger values = more curling.
    if (jitter_lat_nm > 0) {
      jit1 <- cumsum(rnorm(n, 0, jitter_lat_nm))
      jit2 <- if (use_3d) cumsum(rnorm(n, 0, jitter_lat_nm)) else rep(0, n)
    } else {
      jit1 <- rep(0, n)
      jit2 <- rep(0, n)
    }

    # Helical displacement: cosine in p1, sine in p2 → right-handed helix
    if (helix_r > 0 && helix_pitch > 0) {
      phi  <- 2 * pi * s / helix_pitch
      hel1 <- helix_r * cos(phi)
      hel2 <- helix_r * sin(phi)
    } else {
      hel1 <- rep(0, n)
      hel2 <- rep(0, n)
    }

    d1 <- jit1 + hel1   # total displacement in p1 direction
    d2 <- jit2 + hel2   # total displacement in p2 direction (0 in 2D XY)

    df <- data.frame(
      x = s * ax + d1 * p1x + d2 * p2x,
      y = s * ay + d1 * p1y + d2 * p2y
    )
    if (use_3d) df$z <- s * az + d1 * p1z + d2 * p2z
    # Store filament axis tangent for direction-aware linker displacement.
    # All molecules share the same global axis direction (straight filament).
    df$ldir_x    <- ax
    df$ldir_y    <- ay
    df$ldir_z    <- if (use_3d) az else 0
    df$ldir_type <- "perp"
    df
  }

  # Hollow sphere - uniform surface sampling with adjustable shell thickness.
  # In 2D collapses to a ring (circle cross-section).
  gen_hollow_sphere <- function(n, radius, shell_thick, use_3d = FALSE) {
    if (!use_3d) {
      r     <- rnorm(n, mean = radius, sd = max(0.1, shell_thick / 4))
      theta <- runif(n, 0, 2 * pi)
      df <- data.frame(x = r * cos(theta), y = r * sin(theta))
    } else {
      phi   <- acos(2 * runif(n) - 1)
      theta <- runif(n, 0, 2 * pi)
      r     <- radius + runif(n, -shell_thick / 2, shell_thick / 2)
      df <- data.frame(x = r * sin(phi) * cos(theta),
                       y = r * sin(phi) * sin(theta),
                       z = r * cos(phi))
    }
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  # Double ring - two parallel rings separated along Z (or offset along Y in 2D).
  # Useful for nuclear pore complex (~35 nm separation), centrioles, etc.
  gen_double_ring <- function(n, diameter, ring_width, separation_nm,
                               z_thickness = 0, use_3d = FALSE) {
    radius <- diameter / 2
    n1 <- round(n / 2); n2 <- n - n1
    make_one <- function(n_pts, z_pos) {
      theta <- runif(n_pts, 0, 2 * pi)
      r     <- rnorm(n_pts, mean = radius, sd = max(0.1, ring_width / 4))
      df    <- data.frame(x = r * cos(theta), y = r * sin(theta))
      if (use_3d)
        df$z <- rnorm(n_pts, z_pos, max(0.1, z_thickness / 4))
      df
    }
    df <- rbind(make_one(n1, -separation_nm / 2),
                make_one(n2,  separation_nm / 2))
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  # Arc - a partial ring spanning arc_extent_deg degrees.
  # arc_start_deg sets the angular position of the start of the arc.
  gen_arc <- function(n, diameter, ring_width, arc_start_deg, arc_extent_deg,
                       z_thickness = 0, use_3d = FALSE) {
    radius     <- diameter / 2
    start_rad  <- arc_start_deg  * pi / 180
    extent_rad <- max(1, arc_extent_deg) * pi / 180
    theta <- runif(n, start_rad, start_rad + extent_rad)
    r     <- rnorm(n, mean = radius, sd = max(0.1, ring_width / 4))
    df    <- data.frame(x = r * cos(theta), y = r * sin(theta))
    if (use_3d)
      df$z <- rnorm(n, 0, max(0.1, z_thickness / 4))
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  # Oligomeric ring stack - fixed, evenly-spaced subunit positions (no randomness).
  # Models protein complexes: hexamers, heptamers,
  # nuclear pore subcomplexes, etc.
  # GT positions are rigidly determined: exactly n_subunits × n_stacks molecules.
  # angle_offset_deg: each successive ring is rotated by this amount (e.g., for
  # symmetry mismatch between stacked assemblies).
  # In 3D, rings are separated along Z by separation_nm.
  # In 2D, all rings project onto XY (stacks overlap); use angle_offset_deg to
  # make the rotational offset visible.
  gen_oligoring_stack <- function(n_subunits, diameter, n_stacks,
                                   separation_nm, angle_offset_deg = 0,
                                   base_angle_deg = 0,
                                   use_3d = FALSE) {
    n_subunits <- max(2L, as.integer(n_subunits))
    n_stacks   <- max(1L, as.integer(n_stacks))
    radius     <- diameter / 2
    off_rad    <- angle_offset_deg * pi / 180
    base_rad   <- base_angle_deg   * pi / 180

    layers <- lapply(seq_len(n_stacks), function(k) {
      z_pos  <- (k - (n_stacks + 1) / 2) * separation_nm
      rot    <- base_rad + (k - 1L) * off_rad
      angles <- rot + 2 * pi * (seq_len(n_subunits) - 1L) / n_subunits
      df     <- data.frame(x = radius * cos(angles),
                           y = radius * sin(angles))
      if (use_3d) df$z <- z_pos
      df
    })
    df <- do.call(rbind, layers)
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    rownames(df) <- NULL
    df
  }

  # ── Image-guided helpers ─────────────────────────────────────────────────────

  read_image_gray <- function(path) {
    img   <- magick::image_read(path)
    img_g <- magick::image_convert(img, colorspace = "Gray")
    arr   <- magick::image_data(img_g, "gray")
    # dim(arr) = [1, width, height]  →  [H x W] matrix
    W <- dim(arr)[2]; H <- dim(arr)[3]
    mat <- t(matrix(as.integer(arr[1,,]), nrow = W, ncol = H)) / 255
    rng <- range(mat)
    if (diff(rng) > .Machine$double.eps)
      mat <- (mat - rng[1]) / diff(rng)
    list(mat = mat, img = img_g)
  }

  img_to_b64 <- function(img_magick, max_px = 800L) {
    info <- magick::image_info(img_magick)
    if (max(info$width, info$height) > max_px)
      img_magick <- magick::image_scale(
        img_magick, paste0(max_px, "x", max_px))
    raw_bytes <- magick::image_write(img_magick, format = "png")
    paste0("data:image/png;base64,", base64enc::base64encode(raw_bytes))
  }

  # Sample proteins from confocal image intensity.
  # XY is always drawn from pixel intensity (proportional to brightness).
  # In 3D mode, Z spread scales with local pixel intensity: bright regions
  # receive more axial depth (z_spread_nm), dim regions stay near Z = 0.
  gen_from_image <- function(img_gray, n_proteins, px_size_nm,
                              threshold, min_trace, max_trace,
                              use_3d = FALSE,
                              z_spread_nm = 200, z_bg_nm = 20) {
    m     <- img_gray; m[m < threshold] <- 0
    probs <- as.vector(m)
    if (sum(probs) == 0)
      stop("All pixels are below the threshold - lower the threshold value.")
    probs <- probs / sum(probs)
    H <- nrow(m); W <- ncol(m)

    idx     <- sample(H * W, n_proteins, replace = TRUE, prob = probs)
    row_idx <- ((idx - 1L) %% H) + 1L
    col_idx <- ((idx - 1L) %/% H) + 1L
    x_nm    <- (col_idx - 0.5 + runif(n_proteins, -0.5, 0.5)) * px_size_nm -
               (W / 2) * px_size_nm
    y_nm    <- ((H - row_idx + 0.5) + runif(n_proteins, -0.5, 0.5)) * px_size_nm -
               (H / 2) * px_size_nm
    local_int  <- probs[idx] / max(probs, na.rm = TRUE)   # 0–1
    n_locs_exp <- min_trace + (max_trace - min_trace) * local_int

    if (!use_3d) {
      df <- data.frame(x = x_nm, y = y_nm, n_locs_expected = n_locs_exp)
    } else {
      # 3D: Z spread = pmax(z_bg_nm, z_spread_nm × local_intensity)
      z_sd  <- pmax(max(0.5, z_bg_nm), z_spread_nm * local_int)
      z_nm  <- rnorm(n_proteins, 0, z_sd)
      z_nm  <- pmax(-1000, pmin(1000, z_nm))
      df <- data.frame(x = x_nm, y = y_nm, z = z_nm, n_locs_expected = n_locs_exp)
    }
    df$ldir_type <- NA_character_
    df
  }

  # ── Labeling presets ─────────────────────────────────────────────────────────
  # n_fluor: mean fluorophores per labeled molecule (Poisson-drawn, min 1).
  # link_nm: mean distance fluorophore→true molecule centre (nm).
  # In MINFLUX each blinking fluorophore = one trace ID.
  # n_fluor is the *theoretical* count; actual detected TIDs are ~60–80 % of
  # this due to incomplete blinking/detection efficiency.
  #
  # Key assumptions:
  #   IgG arm (Fab tip to Fc base) ≈ 11 nm; secondary adds another ≈ 11 nm.
  #   Nanobody body ≈ 2.5 nm → effective linker ~3 nm.
  #   NHS-ester DOL for commercial IgG kits: 2–4 (mean ~3).
  #   Monoclonal: 1 antibody per epitope (bivalent, but typically 1 per target).
  #   Polyclonal: 3–5 antibodies per antigen (multiple epitopes accessible).
  #   Secondary binding stoichiometry: 2–4 secondaries per primary IgG.
  # dye_spread_nm: Gaussian sigma of individual dye positions within the label
  # complex (around the shared complex-centre displaced by link_nm).
  # n_fluor > 1 on the SAME complex → dyes cluster within dye_spread_nm;
  # each molecule's complex is displaced independently by link_nm.
  # Single source of truth lives in helpers_shared.R so the SMLM tab's inverse
  # GT-estimator inverts exactly the numbers this forward model generates with.
  LABEL_PRESETS <- minflux_label_presets()

  # Expand labeled molecules to individual fluorophore positions.
  # eff_pct:       labeling efficiency (% of molecules that get labeled)
  # n_fluor:       mean fluorophores per labeled molecule (Poisson-drawn, min 1)
  # link_nm:       mean distance from label-complex centre to molecule centre (nm)
  # dye_spread_nm: Gaussian SD of individual dye positions within the complex
  #
  # Two-level displacement model:
  #   Level 1 (per molecule): one complex-centre offset of ~link_nm, shared by
  #     all fluorophores on that molecule.
  #   Level 2 (per fluorophore): each dye scatters by dye_spread_nm around the
  #     complex centre. For single-dye labels dye_spread_nm = 0.
  # Each returned row = one fluorophore = one future trace ID.
  expand_to_fluorophores <- function(gt_df, eff_pct, n_fluor,
                                      link_nm, dye_spread_nm = 0,
                                      ldir_aware = TRUE,
                                      protein_r_nm = 0) {
    use_3d <- "z" %in% names(gt_df)

    # Step 1: all GT molecules can be labeled; orig_rows tracks indices into
    # the full gt_df so mol_id can look up res$gt in the plotting code.
    orig_rows <- seq_len(nrow(gt_df))
    n_mol     <- nrow(gt_df)
    if (n_mol == 0) {
      base <- gt_df[, intersect(c("x","y","z","n_locs_expected"), names(gt_df)),
                    drop = FALSE]
      base$mol_id <- integer(0)
      return(base)
    }

    # Step 2: draw fluorophore count per molecule (Poisson, min 1)
    if (n_fluor <= 1L) {
      n_per <- rep(1L, n_mol)
    } else {
      n_per <- pmax(1L, rpois(n_mol, n_fluor))
    }

    # Step 3: expand molecule rows to fluorophore rows.
    # mol_id stores the original GT row index for use in the plotting code.
    mol_idx <- rep(seq_len(n_mol), n_per)
    keep    <- intersect(c("x","y","z","n_locs_expected"), names(gt_df))
    f_df    <- gt_df[mol_idx, keep, drop = FALSE]
    f_df$mol_id <- orig_rows[mol_idx]
    n_f     <- nrow(f_df)

    # Step 3b: per-fluorophore detection filter.
    # Labeling efficiency applies to each individual fluorophore independently -
    # this is the biologically relevant quantity for MINFLUX: what fraction of
    # the theoretical TIDs are actually detected?  Molecules where every
    # fluorophore is missed simply produce no TIDs (equivalent to being
    # unlabeled), which emerges naturally without a separate molecule-level step.
    if (eff_pct < 100 && n_f > 0) {
      detected <- runif(n_f) < (eff_pct / 100)
      f_df  <- f_df[detected, , drop = FALSE]
      n_f   <- nrow(f_df)
      if (n_f == 0) return(f_df)
      mol_idx <- mol_idx[detected]
    }

    # Local helper: draw n_mol structured displacement vectors of given magnitude.
    # Respects ldir_type: "perp" (filament) → lateral only; "radial" (ring/sphere)
    # → outward from structure centre; NA → isotropic.
    # Returns a list with elements dx, dy, and (when use_3d) dz.
    struct_displace <- function(mag_vec) {
      ldt <- if (ldir_aware && "ldir_type" %in% names(gt_df))
               gt_df$ldir_type[1L] else NA_character_

      if (!is.na(ldt) && ldt == "perp" &&
          all(c("ldir_x","ldir_y") %in% names(gt_df))) {
        if (use_3d) {
          tang <- cbind(gt_df$ldir_x, gt_df$ldir_y, gt_df$ldir_z)
          u    <- matrix(rnorm(n_mol * 3L), ncol = 3L)
          u    <- u - rowSums(u * tang) * tang
          u_n  <- sqrt(rowSums(u^2)); u_n[u_n == 0] <- 1
          u    <- u / u_n
          list(dx = mag_vec * u[,1L], dy = mag_vec * u[,2L], dz = mag_vec * u[,3L])
        } else {
          sgn <- sample(c(-1L, 1L), n_mol, replace = TRUE)
          list(dx = mag_vec * sgn * (-gt_df$ldir_y),
               dy = mag_vec * sgn * ( gt_df$ldir_x))
        }
      } else if (!is.na(ldt) && ldt == "radial") {
        cx <- if ("struct_cx" %in% names(gt_df)) gt_df$struct_cx else rep(0, n_mol)
        cy <- if ("struct_cy" %in% names(gt_df)) gt_df$struct_cy else rep(0, n_mol)
        rx <- gt_df$x - cx; ry <- gt_df$y - cy
        if (use_3d) {
          cz <- if ("struct_cz" %in% names(gt_df)) gt_df$struct_cz else rep(0, n_mol)
          rz <- gt_df$z - cz
          rr <- sqrt(rx^2 + ry^2 + rz^2); rr[rr == 0] <- 1
          list(dx = mag_vec * rx/rr, dy = mag_vec * ry/rr, dz = mag_vec * rz/rr)
        } else {
          rr <- sqrt(rx^2 + ry^2); rr[rr == 0] <- 1
          list(dx = mag_vec * rx/rr, dy = mag_vec * ry/rr)
        }
      } else {
        if (use_3d) {
          dirs <- matrix(rnorm(n_mol * 3L), ncol = 3L)
          lens <- sqrt(rowSums(dirs^2)); lens[lens == 0] <- 1
          dirs <- dirs / lens
          list(dx = mag_vec * dirs[,1L], dy = mag_vec * dirs[,2L],
               dz = mag_vec * dirs[,3L])
        } else {
          theta <- runif(n_mol, 0, 2 * pi)
          list(dx = mag_vec * cos(theta), dy = mag_vec * sin(theta))
        }
      }
    }

    # Step 3c (Level 0): epitope offset from protein centroid.
    # The antibody binds a surface epitope, not the geometric centre of the protein.
    # Displacement magnitude is exactly the protein radius (r = 0.66 × MW^(1/3) nm);
    # direction follows the same structure-aware rules as the linker arm.
    if (protein_r_nm > 0 && n_mol > 0) {
      d0 <- struct_displace(rep(protein_r_nm, n_mol))
      f_df$x <- f_df$x + d0$dx[mol_idx]
      f_df$y <- f_df$y + d0$dy[mol_idx]
      if (use_3d) f_df$z <- f_df$z + d0$dz[mol_idx]
    }

    # Step 4 (Level 1): one complex-centre displacement per MOLECULE.
    # All fluorophores on the same molecule share this offset.
    # Direction is constrained by the local structure geometry when ldir_aware=TRUE.
    if (link_nm > 0 && n_mol > 0) {
      sigma_link <- link_nm * sqrt(pi / 2)   # half-normal → mean ≈ link_nm
      d1 <- struct_displace(abs(rnorm(n_mol, 0, sigma_link)))
      f_df$x <- f_df$x + d1$dx[mol_idx]
      f_df$y <- f_df$y + d1$dy[mol_idx]
      if (use_3d) f_df$z <- f_df$z + d1$dz[mol_idx]
    }

    # Step 5 (Level 2): per-fluorophore scatter within the complex body.
    # Skipped when dye_spread_nm = 0 (single-dye labels, custom dye_spread = 0).
    if (dye_spread_nm > 0 && n_f > 0) {
      f_df$x <- f_df$x + rnorm(n_f, 0, dye_spread_nm)
      f_df$y <- f_df$y + rnorm(n_f, 0, dye_spread_nm)
      if (use_3d) f_df$z <- f_df$z + rnorm(n_f, 0, dye_spread_nm)
    }

    rownames(f_df) <- NULL
    f_df
  }

  # ── Blinking / TID splitting ────────────────────────────────────────────────
  # Models dye re-activation: each fluorophore may appear as K additional TIDs,
  # K ~ Poisson(lambda). Extra rows carry the same mol_id and position; each TID
  # generates its own localizations independently via gen_localizations, giving
  # apparent over-counting proportional to (1 + lambda) × n_fluor_per_molecule.
  split_for_blinking <- function(fluoro_df, lambda) {
    if (lambda <= 0 || nrow(fluoro_df) == 0L) return(fluoro_df)
    extra <- rpois(nrow(fluoro_df), lambda)
    if (sum(extra) == 0L) return(fluoro_df)
    rbind(fluoro_df, fluoro_df[rep(seq_len(nrow(fluoro_df)), extra), , drop = FALSE])
  }

  # ── Localization generator ──────────────────────────────────────────────────
  # Uses n_locs_expected column in fluoro_df when present (image-guided mode);
  # otherwise draws from Poisson(trace_length) for all fluorophores.
  # sigma_loc_xy / sigma_loc_z are the per-localization precision in nm.
  #
  # use_geometric:  replace Poisson with Geometric(p = 1/lambda).  Same mean
  #   but heavy right tail: constant per-step bleaching gives many short traces
  #   and a few long ones, matching real MINFLUX trace-length distributions.
  #
  # scale_by_fluor: divide each fluorophore's lambda by the number of
  #   fluorophores on its molecule (requires mol_id column in fluoro_df).
  #   Keeps total locs per molecule ≈ trace_length regardless of labeling
  #   strategy - models a fixed observation budget per molecule.

  gen_localizations <- function(fluoro_df, trace_length,
                                 sigma_loc_xy, sigma_loc_z = sigma_loc_xy,
                                 use_geometric  = FALSE,
                                 scale_by_fluor = FALSE) {
    n      <- nrow(fluoro_df)
    use_3d <- "z" %in% names(fluoro_df)
    if (n == 0) {
      cols  <- if (use_3d) c("x","y","z","tid","tim")
               else        c("x","y","tid","tim")
      empty <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
      names(empty) <- cols
      return(empty)
    }

    lambda <- if ("n_locs_expected" %in% names(fluoro_df))
                fluoro_df$n_locs_expected
              else
                rep(max(1, trace_length), n)

    # Divide each fluorophore's budget by how many siblings it has on its molecule
    if (scale_by_fluor && "mol_id" %in% names(fluoro_df)) {
      mol_counts <- ave(rep(1, n), fluoro_df$mol_id, FUN = sum)
      lambda     <- lambda / pmax(1, mol_counts)
    }
    lambda <- pmax(1, lambda)

    if (use_geometric) {
      # rgeom(n, p) = failures before first success (min 0); +1 → min trace = 1
      # mean = 1/p = lambda, heavy right tail from constant bleach probability
      p     <- pmin(0.999, 1 / lambda)
      n_loc <- pmax(1L, rgeom(n, p) + 1L)
    } else {
      n_loc <- pmax(1L, rpois(n, lambda))
    }

    idx    <- rep(seq_len(n), n_loc)
    total  <- length(idx)
    rep_df <- fluoro_df[idx, intersect(c("x","y","z"), names(fluoro_df)),
                        drop = FALSE]
    rep_df$x <- rep_df$x + rnorm(total, 0, sigma_loc_xy)
    rep_df$y <- rep_df$y + rnorm(total, 0, sigma_loc_xy)
    if (use_3d) {
      rep_df$z <- rep_df$z + rnorm(total, 0, sigma_loc_z)
      rep_df$z <- pmax(-1000, pmin(1000, rep_df$z))
    }
    rep_df$tid <- as.character(idx)
    rep_df$tim <- seq(0, by = 10, length.out = total)
    rownames(rep_df) <- NULL
    rep_df
  }

  # ── Confocal PSF computation from optical parameters ────────────────────────
  # Gaussian approximation of the lateral confocal PSF sigma.
  #
  #   σ_ex   = 0.21 × λ_ex / NA               excitation PSF lateral sigma
  #   σ_em   = 0.21 × λ_em / NA               emission PSF lateral sigma
  #   1 AU   = 1.22 × λ_em / NA  (nm)         one Airy Unit diameter (sample plane)
  #   r_ph   = AU × (1 AU_nm) / 2             pinhole radius in nm
  #   σ_det  = √(σ_em² + (r_ph / 2.355)²)    detection sigma incl. pinhole size
  #   σ_conf = σ_ex × σ_det / √(σ_ex² + σ_det²)   product of two Gaussians
  #
  # For 640 nm / NA 1.4 / 0.83 AU:  σ ≈ 80 nm, FWHM ≈ 188 nm
  # For 560 nm / NA 1.4 / 1.0 AU:   σ ≈ 72 nm, FWHM ≈ 170 nm
  # For 488 nm / NA 1.4 / 1.0 AU:   σ ≈ 63 nm, FWHM ≈ 148 nm
  compute_confocal_sigma <- function(lambda_ex, lambda_em, na, au) {
    sigma_ex  <- 0.21 * lambda_ex / na
    sigma_em  <- 0.21 * lambda_em / na
    au_nm     <- 1.22 * lambda_em / na        # 1 AU diameter in nm (sample plane)
    ph_radius <- au * au_nm / 2               # pinhole radius in nm
    sigma_det <- sqrt(sigma_em^2 + (ph_radius / 2.355)^2)
    max(1, sigma_ex * sigma_det / sqrt(sigma_ex^2 + sigma_det^2))
  }

  # Fraction of a point emitter's light that passes the pinhole - the encircled
  # energy of the Airy pattern within the pinhole radius. v is the optical radius
  # (v = 3.8317 at exactly 1 AU, the first Airy minimum, where ~84 % is enclosed).
  # This drives the real confocal trade-off: closing the pinhole (AU < 1) rejects
  # light - dimmer, noisier - but tightens the PSF and sections better; opening it
  # (AU > 1) barely adds signal yet blurs and lets in out-of-focus haze.
  confocal_collection_efficiency <- function(au) {
    v <- 3.8317 * max(0, au)
    max(0.01, min(1, 1 - besselJ(v, 0)^2 - besselJ(v, 1)^2))
  }


  # ── Diffraction-limited image renderer ──────────────────────────────────────
  # Renders each GT molecule as a 2D Gaussian PSF on a pixel grid, sums them,
  # adds uniform background, then draws Poisson shot noise.
  # Only computes the PSF within a 3-sigma radius per molecule (sparse update)
  # so performance stays acceptable for up to ~5 000 molecules.
  # Grid is capped to 500 × 500 px.
  # px_nm_user: when supplied, that pixel size is used directly (honouring a
  # user-chosen camera/scanner pixel size in nm).  If the structure extent
  # would push the grid over 500 px at that pixel size, px_nm is silently
  # increased just enough to stay within the cap.
  render_diffraction_image <- function(x_nm, y_nm, psf_sigma_nm,
                                        photons_per_mol, bg_photons_per_px,
                                        px_nm_user = NULL,
                                        photons_vec = NULL,
                                        extent = NULL) {
    n <- length(x_nm)
    if (n == 0) return(NULL)
    pad  <- 4 * psf_sigma_nm
    if (!is.null(extent) && length(extent) == 4 && all(is.finite(extent))) {
      # Explicit window (pinned to the FOV) so the confocal view spans the same
      # ROI as the MINFLUX view instead of cropping to the data bounding box.
      x0 <- extent[1]; x1 <- extent[2]; y0 <- extent[3]; y1 <- extent[4]
    } else {
      x0 <- min(x_nm) - pad; x1 <- max(x_nm) + pad
      y0 <- min(y_nm) - pad; y1 <- max(y_nm) + pad
    }
    span <- max(x1 - x0, y1 - y0, 1)
    # Use user pixel size if given; auto-size otherwise.
    # Either way, enforce the 500 px cap.
    px_nm <- if (!is.null(px_nm_user) && px_nm_user > 0)
               max(px_nm_user, span / 500)
             else
               max(psf_sigma_nm / 2, span / 500)
    nx    <- max(10L, min(500L, ceiling((x1 - x0) / px_nm)))
    ny    <- max(10L, min(500L, ceiling((y1 - y0) / px_nm)))
    # Keep pixels square at the requested size (don't rescale to fit)
    py_nm <- px_nm
    xs <- x0 + (seq_len(nx) - 0.5) * px_nm
    ys <- y0 + (seq_len(ny) - 0.5) * py_nm
    img  <- matrix(0, nrow = ny, ncol = nx)
    r_px <- ceiling(3 * psf_sigma_nm / min(px_nm, py_nm))
    for (i in seq_len(n)) {
      xi <- floor((x_nm[i] - x0) / px_nm) + 1L
      yi <- floor((y_nm[i] - y0) / py_nm) + 1L
      ix <- max(1L, xi - r_px):min(nx, xi + r_px)
      iy <- max(1L, yi - r_px):min(ny, yi + r_px)
      if (length(ix) == 0 || length(iy) == 0) next
      gx   <- exp(-0.5 * ((xs[ix] - x_nm[i]) / psf_sigma_nm)^2)
      gy   <- exp(-0.5 * ((ys[iy] - y_nm[i]) / psf_sigma_nm)^2)
      ph_i <- if (!is.null(photons_vec)) photons_vec[i] else photons_per_mol
      img[iy, ix] <- img[iy, ix] + ph_i * outer(gy, gx)
    }
    img <- img + bg_photons_per_px
    img <- matrix(rpois(length(img), as.vector(img)), nrow = ny)
    list(img = img, xs = xs, ys = ys, px_nm = px_nm)
  }

  # ── Linkage presets ──────────────────────────────────────────────────────────
  LINK_PRESETS <- c(
    "Direct tag (HaloTag / SNAP / click-EdU)" = 2,
    "Nanobody"                                 = 4,
    "Primary antibody (IgG)"                  = 10,
    "Primary + secondary antibody"            = 20,
    "Custom"                                   = NA_real_
  )

  # ── Phase 1 simulation (button-gated) ───────────────────────────────────────
  sim_result <- eventReactive(input$sim_run, {
    if (isTRUE(input$sim_fix_seed)) {
      seed <- suppressWarnings(as.integer(input$sim_seed))
      if (!is.finite(seed)) seed <- 42L
    } else {
      seed <- sample.int(.Machine$integer.max, 1L)
    }
    set.seed(seed)

    if (isTRUE(input$sim_fix_label_seed)) {
      label_seed <- suppressWarnings(as.integer(input$sim_label_seed))
      if (!is.finite(label_seed)) label_seed <- 1L
    } else {
      label_seed <- sample.int(.Machine$integer.max, 1L)
    }

    use_3d <- identical(input$sim_view_mode, "3d")
    n_prot <- max(1L, as.integer(input$sim_n_proteins))

    # Acquisition mode adds the FOV, multi-structure placement, timing model and
    # rate-based background; Normal mode is the simple structure + labeling view.
    acq_mode <- identical(input$sim_mode %||% "normal", "acquisition")
    # FOV active → background spreads across the whole FOV (see sim_use_fov()).
    use_fov  <- sim_use_fov()
    fov_x    <- max(1, as.numeric(input$sim_fov_x %||% 4000))
    fov_y    <- max(1, as.numeric(input$sim_fov_y %||% 4000))
    fov_z    <- max(1, as.numeric(input$sim_fov_z %||% 1000))

    img_b64  <- NULL
    img_w_nm <- NULL
    img_h_nm <- NULL

    if (identical(input$sim_shape, "image")) {
      shiny::validate(shiny::need(
        !is.null(input$sim_img_file), "Please upload a confocal image."))
      img_read <- read_image_gray(input$sim_img_file$datapath)
      img_gray <- img_read$mat
      img_b64  <- img_to_b64(img_read$img)
      px_nm    <- max(1, as.numeric(input$sim_img_px_size))
      img_w_nm <- ncol(img_gray) * px_nm
      img_h_nm <- nrow(img_gray) * px_nm
      gt <- gen_from_image(
        img_gray    = img_gray,
        n_proteins  = n_prot,
        px_size_nm  = px_nm,
        threshold   = suppressWarnings(as.numeric(input$sim_img_threshold)),
        min_trace   = max(1L, as.integer(input$sim_img_min_trace)),
        max_trace   = max(1L, as.integer(input$sim_img_max_trace)),
        use_3d      = use_3d,
        z_spread_nm = max(0, as.numeric(input$sim_img_z_spread %||% 200)),
        z_bg_nm     = max(0, as.numeric(input$sim_img_z_bg     %||% 20))
      )
    } else {
      gt <- switch(input$sim_shape,
        cluster = gen_cluster(
          n = n_prot,
          radius       = max(1, as.numeric(input$sim_cluster_radius)),
          distribution = input$sim_cluster_dist,
          use_3d       = use_3d
        ),
        ring = gen_ring(
          n = n_prot,
          diameter    = max(1, as.numeric(input$sim_ring_diam)),
          ring_width  = max(0, as.numeric(input$sim_ring_width)),
          z_thickness = if (use_3d)
                          max(0, as.numeric(input$sim_ring_z_thick))
                        else 0,
          use_3d = use_3d
        ),
        ellipse = gen_ellipse(
          n = n_prot,
          semi_a = max(1, as.numeric(input$sim_ellipse_a)),
          semi_b = max(1, as.numeric(input$sim_ellipse_b)),
          semi_c = if (use_3d) max(1, as.numeric(input$sim_ellipse_c)) else 1,
          filled = identical(input$sim_ellipse_fill, "filled"),
          use_3d = use_3d
        ),
        filament = gen_filament(
          n             = n_prot,
          spacing_nm    = max(0.1, as.numeric(input$sim_fil_spacing)),
          angle_deg     = suppressWarnings(as.numeric(input$sim_fil_angle)),
          jitter_lat_nm = max(0, as.numeric(input$sim_fil_jitter_lat)),
          jitter_ax_nm  = max(0, as.numeric(input$sim_fil_jitter_ax)),
          helix_r       = if (isTRUE(input$sim_fil_helix))
                            max(0, as.numeric(input$sim_fil_helix_r))
                          else 0,
          helix_pitch   = max(0.1, as.numeric(input$sim_fil_helix_pitch)),
          use_3d        = use_3d,
          z_angle_deg   = if (use_3d)
                            suppressWarnings(as.numeric(input$sim_fil_z_angle))
                          else 0
        ),
        hollow_sphere = gen_hollow_sphere(
          n           = n_prot,
          radius      = max(1, as.numeric(input$sim_hsph_radius)),
          shell_thick = max(0, as.numeric(input$sim_hsph_thick)),
          use_3d      = use_3d
        ),
        double_ring = gen_double_ring(
          n              = n_prot,
          diameter       = max(1, as.numeric(input$sim_dring_diam)),
          ring_width     = max(0, as.numeric(input$sim_dring_width)),
          separation_nm  = max(0, as.numeric(input$sim_dring_sep)),
          z_thickness    = if (use_3d)
                             max(0, as.numeric(input$sim_dring_z_thick))
                           else 0,
          use_3d         = use_3d
        ),
        arc = gen_arc(
          n              = n_prot,
          diameter       = max(1, as.numeric(input$sim_arc_diam)),
          ring_width     = max(0, as.numeric(input$sim_arc_width)),
          arc_start_deg  = suppressWarnings(as.numeric(input$sim_arc_start)),
          arc_extent_deg = max(1, as.numeric(input$sim_arc_extent)),
          z_thickness    = if (use_3d)
                             max(0, as.numeric(input$sim_arc_z_thick))
                           else 0,
          use_3d         = use_3d
        ),
        oligoring = gen_oligoring_stack(
          n_subunits       = max(2L, as.integer(input$sim_oligoring_nsub       %||% 6L)),
          diameter         = max(1,  as.numeric(input$sim_oligoring_diam       %||% 12)),
          n_stacks         = max(1L, as.integer(input$sim_oligoring_nstack     %||% 1L)),
          separation_nm    = max(0.5, as.numeric(input$sim_oligoring_sep       %||% 7)),
          angle_offset_deg = suppressWarnings(as.numeric(input$sim_oligoring_angle_off  %||% 0)),
          base_angle_deg   = suppressWarnings(as.numeric(input$sim_oligoring_base_angle %||% 0)),
          use_3d           = use_3d
        )
      )
      # Oligoring GT count is fixed by structure; update n_prot so downstream
      # labeling code (labeling efficiency, steric check) uses the correct N.
      if (identical(input$sim_shape, "oligoring"))
        n_prot <- nrow(gt)
    }

    # Clamp GT Z to ±1000 nm regardless of shape parameters
    if (use_3d && "z" %in% names(gt))
      gt$z <- pmax(-1000, pmin(1000, gt$z))

    # ── Random XY orientation ────────────────────────────────────────────────
    # Rotates the entire GT cloud by a uniformly sampled angle.  Skipped in
    # multi-copy mode (each copy is rotated independently there) and in image
    # mode (the image defines orientation implicitly).
    if (isTRUE(input$sim_rotate) && !identical(input$sim_shape, "image") &&
        !isTRUE(input$sim_multi_copy)) {
      rot_range <- max(0, as.numeric(input$sim_rotate_range %||% 180))
      angle_rad <- runif(1, -rot_range * pi / 180, rot_range * pi / 180)
      cos_a <- cos(angle_rad); sin_a <- sin(angle_rad)
      x_new <- gt$x * cos_a - gt$y * sin_a
      y_new <- gt$x * sin_a + gt$y * cos_a
      gt$x  <- x_new
      gt$y  <- y_new
      # Rotate stored filament tangent with the same angle
      if ("ldir_type" %in% names(gt) && identical(gt$ldir_type[1L], "perp")) {
        lx_new <- gt$ldir_x * cos_a - gt$ldir_y * sin_a
        ly_new <- gt$ldir_x * sin_a + gt$ldir_y * cos_a
        gt$ldir_x <- lx_new; gt$ldir_y <- ly_new
      }
    }

    # ── Multi-structure placement ────────────────────────────────────────────
    # Place N copies of the base structure at random positions in a rectangular
    # FOV (centred at the origin), each independently rotated.  Two drivers:
    #   • Normal mode  → optional "Multi-copy tiling" (square copies FOV).
    #   • Acquisition mode → "Number of structures" placed in the X×Y FOV.
    # Disabled for image mode because the image already fills the full FOV.
    if (acq_mode && !identical(input$sim_shape, "image")) {
      do_tiling <- TRUE
      n_copies  <- max(1L, as.integer(input$sim_acq_n_structures %||% 1L))
      tile_w    <- fov_x; tile_h <- fov_y
      rot_ea    <- TRUE
    } else if (!acq_mode && isTRUE(input$sim_multi_copy) &&
               !identical(input$sim_shape, "image")) {
      do_tiling <- TRUE
      n_copies  <- max(2L, as.integer(input$sim_copies_n   %||% 5L))
      tile_w    <- tile_h <- max(1, as.numeric(input$sim_copies_fov %||% 2000))
      rot_ea    <- isTRUE(input$sim_copies_rotate)
    } else {
      do_tiling <- FALSE
    }

    if (do_tiling && n_copies > 1L) {
      copies <- lapply(seq_len(n_copies), function(i) {
        cp <- gt
        if (rot_ea) {
          a     <- runif(1, 0, 2 * pi)
          ca    <- cos(a); sa <- sin(a)
          x_rot <- cp$x * ca - cp$y * sa
          y_rot <- cp$x * sa + cp$y * ca
          cp$x  <- x_rot; cp$y <- y_rot
          # Rotate filament tangent with each per-copy rotation
          if ("ldir_type" %in% names(cp) && identical(cp$ldir_type[1L], "perp")) {
            lx_r <- cp$ldir_x * ca - cp$ldir_y * sa
            ly_r <- cp$ldir_x * sa + cp$ldir_y * ca
            cp$ldir_x <- lx_r; cp$ldir_y <- ly_r
          }
        }
        # Keep the first structure near the centre (~0-100 nm), like a real
        # MINFLUX ROI is usually centred on a feature; the rest tile randomly.
        if (i == 1L) {
          off_x <- runif(1, -50, 50)
          off_y <- runif(1, -50, 50)
        } else {
          off_x <- runif(1, -tile_w / 2, tile_w / 2)
          off_y <- runif(1, -tile_h / 2, tile_h / 2)
        }
        cp$x <- cp$x + off_x
        cp$y <- cp$y + off_y
        # Track per-copy structure centre for radial displacement
        if ("ldir_type" %in% names(cp) && identical(cp$ldir_type[1L], "radial")) {
          if ("struct_cx" %in% names(cp)) cp$struct_cx <- cp$struct_cx + off_x
          if ("struct_cy" %in% names(cp)) cp$struct_cy <- cp$struct_cy + off_y
        }
        if (use_3d && "z" %in% names(cp))
          cp$z <- pmax(-1000, pmin(1000, cp$z))
        cp
      })
      gt <- do.call(rbind, copies)
      rownames(gt) <- NULL
    }

    # ── Labeling: preset or custom ───────────────────────────────────────────
    set.seed(label_seed)
    label_key <- input$sim_label_type %||% "Custom"
    if (identical(label_key, "Custom")) {
      n_fluor       <- max(1L, as.integer(input$sim_fluor_custom     %||% 1L))
      link_nm       <- max(0,  as.numeric(input$sim_link_custom_nm   %||% 5))
      dye_spread_nm <- max(0,  as.numeric(input$sim_dye_spread_custom %||% 0))
      footprint_nm  <- max(0,  as.numeric(input$sim_footprint_custom  %||% (link_nm * 1.5)))
      blink_lambda  <- max(0,  as.numeric(input$sim_blink_mean        %||% 0.3))
      lbl_note      <- sprintf("%d fluorophore(s)/molecule, %.0f nm linker, %.0f nm dye spread",
                               n_fluor, link_nm, dye_spread_nm)
    } else {
      p             <- LABEL_PRESETS[[label_key]]
      n_fluor       <- p$n_fluor
      link_nm       <- p$link_nm
      dye_spread_nm <- p$dye_spread_nm
      footprint_nm  <- p$footprint_nm
      blink_lambda  <- if (!is.null(p$blink_lambda)) p$blink_lambda else 0.3
      lbl_note      <- p$note
    }
    eff_pct        <- max(1, min(100, as.numeric(input$sim_label_eff)))
    protein_mw_kda <- max(0, as.numeric(input$sim_protein_mw %||% 0))
    protein_r_nm   <- if (protein_mw_kda > 0) 0.66 * protein_mw_kda^(1/3) else 0
    fluoro  <- expand_to_fluorophores(gt, eff_pct, n_fluor, link_nm, dye_spread_nm,
                                       ldir_aware   = isTRUE(input$sim_ldir_aware),
                                       protein_r_nm = protein_r_nm)
    # Unique physical fluorophores (one dye = one position), before any blinking
    # or rate-based over-counting.  The confocal view renders from THIS so a
    # re-detected molecule is not rendered brighter than it physically is.
    fluoro_pool <- fluoro

    # ── Blinking / TID splitting (Normal mode only) ──────────────────────────
    # In Acquisition mode over-counting comes from the detection rate × duration
    # below, so the λ-based blink split is skipped there.
    n_fluoro_base <- nrow(fluoro)   # unique fluorophores (detectable pool)
    if (!acq_mode && isTRUE(input$sim_blink_model)) {
      bl_lambda <- max(0, as.numeric(input$sim_blink_mean %||% blink_lambda))
      fluoro    <- split_for_blinking(fluoro, bl_lambda)
    }

    # ── Localization precision: manual or CRLB-derived ──────────────────────
    # CRLB mode: σ ≈ L / (2√N_eff),  N_eff = N_photons × (1 − bg_fraction).
    # Simplified 2D MINFLUX CRLB (equal pattern arms, ideal PSF, on-axis).
    if (identical(input$sim_prec_mode %||% "manual", "crlb")) {
      crlb_L   <- max(1,    as.numeric(input$sim_crlb_L       %||% 100))
      crlb_ph  <- max(1,    as.numeric(input$sim_crlb_photons %||% 50))
      crlb_bg  <- max(0, min(0.99, as.numeric(input$sim_crlb_bg %||% 0.05)))
      n_eff    <- max(1, crlb_ph * (1 - crlb_bg))
      sigma_loc_xy <- max(0.1, crlb_L / (2 * sqrt(n_eff)))
      sigma_loc_z  <- if (use_3d) {
        z_fac <- max(0.1, as.numeric(input$sim_crlb_z_factor %||% 0.7))
        max(0.1, sigma_loc_xy * z_fac)
      } else sigma_loc_xy
    } else {
      sigma_loc_xy <- max(0.1, as.numeric(input$sim_precision_xy))
      sigma_loc_z  <- if (use_3d)
                        max(0.1, as.numeric(input$sim_precision_z %||% sigma_loc_xy))
                      else sigma_loc_xy
    }
    trace_length <- max(1, as.numeric(input$sim_trace_length))

    # ── Acquisition window + over-counting (Acquisition mode only) ────────────
    # The user sets the duration (hours + minutes) AND the structure detection
    # rate (traces/hour).  Requested structure traces = rate × duration.  Each
    # fluorophore yields a finite number of detections before PERMANENT bleaching
    # (mean = sim_acq_blinks), so the total is capped at that photon budget: a
    # short run under-samples the pool (low % resolved), a longer run resolves
    # more and over-counts, and beyond the budget the structure plateaus (dead)
    # while background keeps accumulating.  Detections are drawn WITHOUT
    # replacement from the per-fluorophore capacity, so molecules with more
    # fluorophores / more blinks are found sooner.
    # (The resolution/saturation curves are computed live in sim_acq_live();
    # here we only need the over-counted detection list for the MINFLUX view.)
    fov_area_nm2 <- fov_x * fov_y
    acq_T_sec    <- 0
    struct_rate  <- max(0, as.numeric(input$sim_struct_rate %||% 20))   # traces/hour
    n_pool       <- nrow(fluoro)        # unique fluorophores (detectable pool)
    if (acq_mode) {
      hrs       <- max(0, as.numeric(input$sim_acq_hours   %||% 0))
      mins      <- max(0, as.numeric(input$sim_acq_minutes %||% 0))
      acq_T_sec <- max(1, (hrs * 60 + mins) * 60)
      blinks    <- max(1, as.numeric(input$sim_acq_blinks %||% 5))   # detections/fluorophore
      dec_min   <- max(0.1, as.numeric(input$sim_acq_decay_min  %||% 15))
      flr_frac  <- max(0, min(1, as.numeric(input$sim_acq_floor_frac %||% 0.25)))
      det_model <- input$sim_acq_det_model %||% "decay"
      n_req     <- as.integer(round(                                # absolute-rate count
                     acq_expected_count(struct_rate / 60, acq_T_sec,  # per hour → per min
                                        dec_min, flr_frac, det_model)))
      if (n_pool > 0 && n_req > 0) {
        k_i        <- 1L + rpois(n_pool, max(0, blinks - 1))         # per-fluorophore budget
        struct_cap <- sum(k_i)
        n_s        <- min(n_req, struct_cap, 200000L)
        det_idx    <- if (struct_cap <= 5e6)
                        sample(rep(seq_len(n_pool), k_i), n_s)        # exact capacity
                      else
                        sample(seq_len(n_pool), n_s, replace = TRUE, prob = k_i)
        fluoro     <- fluoro[det_idx, , drop = FALSE]
        rownames(fluoro) <- NULL
      } else {
        fluoro     <- fluoro[0, , drop = FALSE]
      }
    }

    locs <- gen_localizations(fluoro, trace_length, sigma_loc_xy, sigma_loc_z,
                               use_geometric  = isTRUE(input$sim_tl_geometric),
                               scale_by_fluor = isTRUE(input$sim_tl_scale_fluor))

    # ── Background / noise TIDs ──────────────────────────────────────────────
    # Acquisition mode: count = background detection rate (traces/hour) × duration,
    # mirroring the structure - longer runs accumulate more non-specific TIDs.
    # Normal mode: a plain fixed count.  Background is spread across the FOV when
    # set, else the multi-copy FOV, else the GT bounding box + 10 % padding.
    # Each background detection is a fresh random position. TID numbering
    # continues from the last real TID.
    n_bg   <- 0L
    bg_pos <- NULL
    bg_rate <- max(0, as.numeric(input$sim_bg_rate %||% 5))   # traces/hour
    if (isTRUE(input$sim_bg_noise)) {
      copies_mode <- !acq_mode && isTRUE(input$sim_multi_copy) &&
                     !identical(input$sim_shape, "image")
      if (use_fov) {
        bg_xc <- 0; bg_yc <- 0
        bg_hx <- fov_x / 2; bg_hy <- fov_y / 2
      } else if (copies_mode) {
        fov_half <- max(1, as.numeric(input$sim_copies_fov %||% 2000)) / 2
        bg_xc <- 0; bg_yc <- 0
        bg_hx <- fov_half; bg_hy <- fov_half
      } else {
        pad   <- 0.1
        x_rng <- range(gt$x, na.rm = TRUE); y_rng <- range(gt$y, na.rm = TRUE)
        x_pad <- diff(x_rng) * pad; y_pad <- diff(y_rng) * pad
        bg_xc <- mean(x_rng); bg_yc <- mean(y_rng)
        bg_hx <- diff(x_rng) / 2 + x_pad; bg_hy <- diff(y_rng) / 2 + y_pad
      }

      if (acq_mode) {
        n_bg <- as.integer(round(bg_rate * acq_T_sec / 3600))   # per hour × hours
      } else {
        n_bg <- as.integer(input$sim_bg_n %||% 20L)
      }
      n_bg <- max(0L, min(200000L, n_bg))

      if (n_bg > 0L) {
        bg_df <- data.frame(x = runif(n_bg, bg_xc - bg_hx, bg_xc + bg_hx),
                            y = runif(n_bg, bg_yc - bg_hy, bg_yc + bg_hy))
        if (use_3d && "z" %in% names(gt)) {
          z_rng  <- range(gt$z, na.rm = TRUE)
          # FOV mode → fill the FOV Z slab; otherwise track the GT extent, with
          # a precision-based fallback for degenerate (planar) structures.
          z_half <- if (use_fov)
                      fov_z / 2
                    else if (diff(z_rng) < 1) max(50, sigma_loc_z * 10)
                    else diff(z_rng) / 2 * 1.1
          z_mid  <- if (use_fov) 0 else mean(z_rng)
          bg_df$z <- runif(n_bg, z_mid - z_half, z_mid + z_half)
        }
        bg_df$mol_id <- NA_integer_
        # Save XY(Z) positions so the confocal view can render bg spots
        bg_pos <- bg_df[, intersect(c("x","y","z"), names(bg_df)), drop = FALSE]
        # Each background TID gets its own locs via gen_localizations
        bg_locs <- gen_localizations(bg_df, trace_length, sigma_loc_xy,
                                      sigma_loc_z,
                                      use_geometric  = isTRUE(input$sim_tl_geometric),
                                      scale_by_fluor = FALSE)
        # Offset TID numbers so they don't clash with real TIDs
        bg_locs$tid <- as.character(as.integer(bg_locs$tid) + nrow(fluoro))
        locs <- rbind(locs, bg_locs)
      }
    }

    # ── RMSE: per-TID displacement error (GT position vs trace median) ──────────
    # Only computed for real TIDs (1 … nrow(fluoro)); background TIDs excluded.
    # Pre-split locs by TID for O(1) lookup per trace - essential when multi-
    # copy produces thousands of TIDs.  Capped at 1 000 sampled TIDs so the
    # computation stays fast even with very large multi-copy scenes.
    rmse_vals <- NULL
    if ("mol_id" %in% names(fluoro) && nrow(gt) > 0 && nrow(locs) > 0) {
      n_real    <- nrow(fluoro)
      tid_samp  <- if (n_real > 1000L) sample(n_real, 1000L) else seq_len(n_real)
      real_rows <- which(suppressWarnings(as.integer(locs$tid)) <= n_real &
                         !is.na(suppressWarnings(as.integer(locs$tid))))
      locs_by_tid <- split(real_rows, locs$tid[real_rows])
      errs <- vapply(tid_samp, function(k) {
        rows <- locs_by_tid[[as.character(k)]]
        if (is.null(rows) || length(rows) == 0) return(NA_real_)
        lk <- locs[rows, , drop = FALSE]
        gi <- fluoro$mol_id[k]
        if (is.na(gi) || gi < 1L || gi > nrow(gt)) return(NA_real_)
        dx <- median(lk$x) - gt$x[gi]
        dy <- median(lk$y) - gt$y[gi]
        dz <- if (use_3d && "z" %in% names(lk) && "z" %in% names(gt))
                median(lk$z) - gt$z[gi]
              else 0
        sqrt(dx^2 + dy^2 + dz^2)
      }, numeric(1L))
      rmse_vals <- errs[!is.na(errs)]
    }

    # ── Steric exclusion check ───────────────────────────────────────────────
    # Compute mean nearest-neighbour distance among GT molecule positions and
    # compare it to the label complex footprint.  A footprint larger than the
    # protein spacing means antibodies physically cannot label every molecule.
    steric_warn <- NULL
    gt_nnd      <- NA_real_
    if (nrow(gt) >= 2L) {
      gt_pts <- as.matrix(gt[, intersect(c("x","y","z"), names(gt)), drop = FALSE])
      samp_i <- if (nrow(gt_pts) > 300L) sample(nrow(gt_pts), 300L) else seq_len(nrow(gt_pts))
      d_mat  <- as.matrix(dist(gt_pts[samp_i, , drop = FALSE]))
      diag(d_mat) <- Inf
      gt_nnd <- mean(apply(d_mat, 1, min), na.rm = TRUE)
    }
    # Effective exclusion zone = label footprint + protein diameter.
    # Proteins physically can't overlap, so even before the antibody is considered
    # the centre-to-centre spacing must be ≥ 2 × protein radius.
    eff_exclusion_nm <- max(footprint_nm, 2 * protein_r_nm, na.rm = TRUE)
    if (!is.na(eff_exclusion_nm) && !is.na(gt_nnd)) {
      if (gt_nnd < eff_exclusion_nm * 0.5) {
        steric_warn <- list(
          level = "danger",
          msg   = sprintf(
            paste0("Severe steric exclusion: combined exclusion zone (~%.0f nm: ",
                   "%.0f nm label footprint + %.0f nm protein diameter) exceeds ",
                   "mean protein spacing (~%.0f nm). Full stoichiometric labeling is ",
                   "physically impossible at this density."),
            eff_exclusion_nm, footprint_nm, 2 * protein_r_nm, gt_nnd)
        )
      } else if (gt_nnd < eff_exclusion_nm) {
        steric_warn <- list(
          level = "warning",
          msg   = sprintf(
            paste0("Steric exclusion likely: combined exclusion zone (~%.0f nm) ",
                   "approaches mean protein spacing (~%.0f nm). Some molecules may be ",
                   "inaccessible; effective labeling efficiency could be lower than expected."),
            eff_exclusion_nm, gt_nnd)
        )
      }
    }

    # ── Snapshot detection times (Acquisition mode) ──────────────────────────
    # Assign front-loaded times to THIS realization's TIDs so the MINFLUX Time
    # colour mode and CSV export have realistic timestamps.  The Acquisition-tab
    # saturation/rate curves are computed live in sim_acq_live() (decoupled so
    # they react to rate/duration/background without re-Generating).
    if (acq_mode && nrow(locs) > 0 && acq_T_sec > 0) {
      all_tids   <- unique(locs$tid)
      n_real     <- nrow(fluoro)
      decay_min  <- max(0.1, as.numeric(input$sim_acq_decay_min  %||% 15))
      floor_frac <- max(0, min(1, as.numeric(input$sim_acq_floor_frac %||% 0.25)))
      det_model  <- input$sim_acq_det_model %||% "decay"
      tid_times  <- setNames(rep(NA_real_, length(all_tids)), all_tids)
      real_ids   <- intersect(as.character(seq_len(n_real)), all_tids)
      if (length(real_ids) > 0)
        tid_times[real_ids] <- sim_sample_detection_times(
          length(real_ids), acq_T_sec, decay_min * 60, floor_frac, det_model)
      bg_ids <- setdiff(all_tids, real_ids)
      if (length(bg_ids) > 0)
        tid_times[bg_ids] <- runif(length(bg_ids), 0, acq_T_sec)
      within   <- ave(seq_len(nrow(locs)), locs$tid,
                      FUN = function(i) (seq_along(i) - 1L) * 0.01)
      locs$tim <- unname(tid_times[locs$tid]) + within
    }

    list(gt = gt, fluoro = fluoro, fluoro_pool = fluoro_pool,
         locs = locs,
         link_nm = link_nm, n_fluor = n_fluor, lbl_note = lbl_note,
         dye_spread_nm = dye_spread_nm, blink_lambda = blink_lambda,
         eff_pct = eff_pct,
         protein_r_nm = protein_r_nm,
         sigma_loc_xy = sigma_loc_xy, sigma_loc_z = sigma_loc_z,
         n_prot = n_prot, n_fluoro = nrow(fluoro),
         n_fluoro_base = n_fluoro_base,
         n_locs = nrow(locs), n_bg = n_bg, use_3d = use_3d,
         img_b64 = img_b64,
         img_w_nm = img_w_nm, img_h_nm = img_h_nm,
         rmse_vals = rmse_vals, bg_pos = bg_pos,
         steric_warn = steric_warn,
         footprint_nm = footprint_nm,
         gt_nnd = gt_nnd, acq_mode = acq_mode,
         blink_active = (!acq_mode) && isTRUE(input$sim_blink_model),
         blink_lambda_used = if ((!acq_mode) && isTRUE(input$sim_blink_model))
                               max(0, as.numeric(input$sim_blink_mean %||% blink_lambda))
                             else 0)
  }, ignoreNULL = TRUE)

  # ── Live acquisition timeline ────────────────────────────────────────────────
  # Recomputes the saturation / detection-rate curves from the last-generated
  # fluorophore pool whenever the acquisition parameters change - so the
  # Acquisition tab reacts to rate / duration / background / bleaching live, with
  # no need to re-Generate the structure.  A fixed seed keeps the curve stable
  # for a given parameter set (it only moves when an input actually changes).
  sim_acq_live <- reactive({
    if (!identical(input$sim_mode %||% "normal", "acquisition")) return(NULL)
    res <- sim_result()
    if (is.null(res) || is.null(res$fluoro_pool)) return(NULL)
    pool_mol <- res$fluoro_pool$mol_id
    if (length(pool_mol) == 0) return(NULL)
    hrs   <- max(0, as.numeric(input$sim_acq_hours   %||% 0))
    mins  <- max(0, as.numeric(input$sim_acq_minutes %||% 0))

    # Emitter-crowding weight per molecule: a molecule with many other labeled
    # molecules inside the donut radius is harder to single out, so its effective
    # discovery rate is scaled by 1/(1 + neighbours).  Aligned to the sorted
    # unique mol_ids (= the order compute_acq_timeline builds f_m in).
    crowd_m <- NULL
    if (isTRUE(input$sim_crowding)) {
      u  <- sort(unique(pool_mol))
      gp <- res$gt[u, , drop = FALSE]
      cols <- intersect(c("x", "y", "z"), names(gp))
      if (length(u) >= 2 && length(cols) >= 2) {
        R   <- max(1, as.numeric(input$sim_crowd_radius %||% 50))
        pts <- as.matrix(gp[, cols, drop = FALSE])
        nn  <- RANN::nn2(pts, pts, k = min(nrow(pts), 64L),
                         searchtype = "radius", radius = R)
        neigh   <- rowSums(nn$nn.idx > 0) - 1L
        crowd_m <- 1 / (1 + neigh)
      }
    }

    set.seed(1L)                          # stable curve per parameter set
    compute_acq_timeline(
      pool_mol    = pool_mol,
      n_mol_total = nrow(res$gt),
      struct_rate = max(0, as.numeric(input$sim_struct_rate %||% 20)),
      T_sec       = max(1, (hrs * 60 + mins) * 60),
      blinks      = max(1, as.numeric(input$sim_acq_blinks %||% 5)),
      decay_min   = max(0.1, as.numeric(input$sim_acq_decay_min %||% 15)),
      floor_frac  = max(0, min(1, as.numeric(input$sim_acq_floor_frac %||% 0.25))),
      bg_rate     = max(0, as.numeric(input$sim_bg_rate %||% 5)),
      resolve_q   = as.numeric(input$sim_acq_resolve_q %||% 0.9),
      crowd_m     = crowd_m,
      model       = input$sim_acq_det_model %||% "decay")
  })

  # ── Calibrate from a real MINFLUX run ────────────────────────────────────────
  # Load a .npy / .json / .csv run, measure its detection-rate profile, duration
  # and FOV, pre-fill the Acquisition inputs, and keep the measured rate curve to
  # overlay on the Acquisition plot.
  sim_real_rv  <- reactiveVal(NULL)   # analysis of the GATED run (drives everything)
  sim_real_raw <- reactiveVal(NULL)   # the run exactly as uploaded, never gated

  # ── Calibration gate ────────────────────────────────────────────────────────
  # GOAL: restrict a real run to the localizations / traces the calibration is
  # measured from, and report exactly how much was dropped.
  #
  # Order is deliberate: EFO and CFR gate individual LOCALIZATIONS first, then
  # trace length is recounted on the survivors and the length gate drops whole
  # TRACES.  So "trace length" here means USABLE localizations per trace - gating
  # hard on EFO and then asking for >= 20 locs is a far stricter request than
  # either gate alone, which is why the info box reports both stages.
  #
  # Non-finite values are KEPT (a gate cannot judge a missing number) and counted
  # separately, so leaving a range at its full extent is always a true no-op.
  # A NULL range is a no-op too: the column is absent, or the widget is not built.
  .rng_ok <- function(r)
    length(r) == 2L && all(is.finite(r)) && r[2] >= r[1]

  # A column is usable only if it exists AND carries at least one finite value.
  .col_ok <- function(df, nm)
    nm %in% names(df) &&
      any(is.finite(suppressWarnings(as.numeric(df[[nm]]))))

  gate_real_run <- function(df, efo_rng = NULL, cfr_rng = NULL, len_rng = NULL) {
    if (is.null(df) || !all(c("tid", "tim", "x", "y") %in% names(df))) return(NULL)
    n_loc_in <- nrow(df)
    n_tr_in  <- length(unique(df$tid))

    # ── Stage 1: localization-level gates (EFO, CFR) ──────────────────────────
    keep    <- rep(TRUE, n_loc_in)
    n_na    <- c(efo = 0L, cfr = 0L)   # values a gate could not judge
    narrowed <- character(0)           # gates the user actually tightened
    for (nm in c("efo", "cfr")) {
      rng <- if (identical(nm, "efo")) efo_rng else cfr_rng
      if (!.rng_ok(rng) || !.col_ok(df, nm)) next
      v    <- suppressWarnings(as.numeric(df[[nm]]))
      fin  <- is.finite(v)
      keep <- keep & (!fin | (v >= rng[1] & v <= rng[2]))
      n_na[nm] <- sum(!fin)
      full <- range(v[fin])
      if (rng[1] > full[1] || rng[2] < full[2]) narrowed <- c(narrowed, toupper(nm))
    }
    g <- df[keep, , drop = FALSE]
    n_loc_mid <- nrow(g)
    n_tr_mid  <- length(unique(g$tid))

    # ── Stage 2: trace-level gate on the SURVIVING localization count ─────────
    if (.rng_ok(len_rng) && n_loc_mid > 0) {
      tl   <- table(g$tid)
      ok   <- names(tl)[tl >= len_rng[1] & tl <= len_rng[2]]
      g    <- g[as.character(g$tid) %in% ok, , drop = FALSE]
      full <- range(as.integer(tl))
      if (len_rng[1] > full[1] || len_rng[2] < full[2])
        narrowed <- c(narrowed, "trace length")
    }

    # Traces with >= 3 localizations are the only ones a within-trace SD can be
    # measured from, so they are the sample size that actually sets sigma.
    n_sig <- if (nrow(g) > 0) sum(table(g$tid) >= 3L) else 0L

    list(df = g, narrowed = narrowed, n_na = n_na,
         n_loc_in = n_loc_in, n_loc_mid = n_loc_mid, n_loc_out = nrow(g),
         n_tr_in  = n_tr_in,  n_tr_mid  = n_tr_mid,
         n_tr_out = length(unique(g$tid)), n_tr_sigma = n_sig)
  }

  # The gated run.  Rebuilt whenever the file, a gate range or the axial factor
  # changes; every calibrated number downstream is measured on THIS subset.
  # The stored z is RAW, so the factor is applied here — one place, so sigma_z
  # and fov_z can never disagree with the number shown in the input.
  sim_real_gated <- reactive({
    raw <- sim_real_raw()
    if (is.null(raw)) return(NULL)
    g <- gate_real_run(raw$df,
                       efo_rng = input$sim_real_efo_range,
                       cfr_rng = input$sim_real_cfr_range,
                       len_rng = input$sim_real_len_range)
    if (is.null(g)) return(NULL)
    # Until the widget exists, fall back to the factor detected for this file,
    # so the very first calibration already uses the right one.
    zc <- suppressWarnings(as.numeric(input$sim_real_z_cal %||% raw$cal$factor))
    if (!is.finite(zc) || zc <= 0) zc <- 1
    if ("z" %in% names(g$df)) g$df$z <- g$df$z * zc
    g$z_cal <- zc
    g$z_src <- raw$cal
    g$has_z <- isTRUE(raw$has_z)
    g
  })

  # Axial calibration for the loaded run.  Rebuilt per upload (so a new file can
  # never inherit the previous file's factor) and pre-filled from what the file
  # itself says: .npy/.json are raw and need the correction; a CSV is returned
  # to raw first, then offered the factor its producer had applied.
  # NOTE: this is the refractive-index axial correction, NOT sim_crlb_z_factor
  # (the sigma_z/sigma_xy precision ratio). Both default near 0.7 and they are
  # unrelated quantities.
  output$sim_real_zcal_ui <- renderUI({
    raw <- sim_real_raw()
    if (is.null(raw)) return(NULL)
    cal <- raw$cal
    # A 2D run has no axial structure to calibrate — say nothing about z.
    if (!isTRUE(raw$has_z))
      return(tags$div(style = "font-size:0.78em; color:#999; margin-bottom:6px;",
                      "2D run — no axial calibration needed."))
    tagList(
      numericInput("sim_real_z_cal", "Axial calibration factor (z x)",
                   value = cal$factor, min = 0.05, max = 5, step = 0.05),
      tags$div(style = "font-size:0.78em; color:#666; margin-top:-6px;",
        sprintf("Detected: %s%s. Applied to z before σ_z and the FOV depth are measured.",
                cal$source, if (isTRUE(cal$certain)) "" else " - assumed"))
    )
  })

  # Gate widgets, built from the uploaded run's own ranges so the default value
  # spans the whole data and is therefore a no-op.  Columns the file does not
  # carry (a plain CSV often has no efo / cfr) say so instead of silently
  # passing everything.
  output$sim_real_gate_ui <- renderUI({
    raw <- sim_real_raw()
    if (is.null(raw)) return(NULL)
    df     <- raw$df
    tl_rng <- range(as.integer(table(df$tid)))

    # No `step` here on purpose: numericRangeInput errors on step = NULL, and
    # EFO/CFR have no natural increment - the browser's default suits both.
    num_rng <- function(id, label, nm) {
      if (!.col_ok(df, nm))
        return(tags$div(style = "font-size:0.78em; color:#999; margin-bottom:6px;",
                        sprintf("%s - not in this file, gate unavailable.",
                                toupper(nm))))
      v <- suppressWarnings(as.numeric(df[[nm]]))
      r <- range(v[is.finite(v)])
      shinyWidgets::numericRangeInput(id, label, value = r,
                                      min = r[1], max = r[2])
    }

    tags$details(
      class = "sim-group", open = NA,
      tags$summary(tags$b("Calibration gate")),
      helpText("Restricts which localizations and traces the precision and",
               "trace-length statistics are measured from. EFO and CFR gate",
               "individual localizations; the length gate then drops whole",
               "traces by how many localizations survived."),
      num_rng("sim_real_efo_range", "EFO range (Hz)",  "efo"),
      num_rng("sim_real_cfr_range", "CFR range",       "cfr"),
      shinyWidgets::numericRangeInput(
        "sim_real_len_range", "Localizations per trace",
        value = tl_rng, min = tl_rng[1], max = tl_rng[2], step = 1),
      uiOutput("sim_real_gate_info")
    )
  })

  # Before / after counts for both gate stages.  The gate is only useful if its
  # severity is visible, so this always reports rather than only on trouble.
  output$sim_real_gate_info <- renderUI({
    g <- sim_real_gated()
    if (is.null(g)) return(NULL)
    pct <- function(a, b) if (b > 0) sprintf("%.0f%%", 100 * a / b) else "-"

    too_few <- g$n_loc_out < 10L || g$n_tr_sigma < 3L
    tags$div(
      class = if (too_few) "alert alert-warning" else "well well-sm",
      style = "padding:6px 10px; font-size:0.78em; margin:6px 0;",
      tags$b("Gate"), tags$br(),
      sprintf("Localizations: %d → %d  (%s kept)",
              g$n_loc_in, g$n_loc_out, pct(g$n_loc_out, g$n_loc_in)), tags$br(),
      sprintf("Traces: %d → %d  (%s kept); %d with ≥ 3 locs set σ",
              g$n_tr_in, g$n_tr_out, pct(g$n_tr_out, g$n_tr_in), g$n_tr_sigma),
      if (any(g$n_na > 0)) tagList(tags$br(),
        sprintf("%s value(s) missing and left ungated.",
                paste(sprintf("%s: %d", toupper(names(g$n_na)[g$n_na > 0]),
                              g$n_na[g$n_na > 0]), collapse = ", "))),
      if (too_few) tagList(tags$br(), tags$b(
        "Too little left to calibrate - widen the gate.")),
      # Gating is not neutral: every one of these selections keeps the brighter,
      # better-centred, longer-lived molecules, which are also the most precise.
      if (length(g$narrowed) > 0 && !too_few) tagList(tags$br(), tags$span(
        style = "color:#854d0e;",
        sprintf("%s narrowed - this keeps the brighter / better-centred / longer traces, so the calibrated σ is tighter than the run as a whole.",
                paste(g$narrowed, collapse = ", "))))
    )
  })

  # Split the measured total rate into structure + constant background given the
  # fraction of traces that are on-structure, and push the fitted MODEL plus its
  # struct_rate / bg_rate / level. Background is the constant offset; the
  # structure is the time-varying part on top, whose shape (decay / constant /
  # rise) analyze_real_run already chose. struct_rate is set to the structure's
  # CHARACTERISTIC rate (peak for decay, plateau for rise), and level to the
  # matching [0,1] ratio, so the simulated curve reproduces the measured profile.
  apply_real_rates <- function(a, frac) {
    frac  <- max(0, min(1, frac))
    model <- a$det_model %||% "decay"
    # Background: constant offset, capped below the run minimum so the structure
    # component stays positive at both ends.
    bg <- max(0, min((1 - frac) * a$rate_mean, 0.9 * min(a$rate_peak, a$rate_floor)))
    s_early <- max(0.1, a$rate_peak  - bg)      # structure rate near the start
    s_late  <- max(0,   a$rate_floor - bg)      # structure rate near the end
    if (model == "constant") {
      struct_rate <- max(0.1, a$rate_mean - bg); level <- 1
    } else if (model == "rise") {               # plateau = late, starts low
      struct_rate <- max(0.1, s_late);  level <- max(0, min(1, s_early / struct_rate))
    } else {                                    # decay: peak = early, → floor
      struct_rate <- max(0.1, s_early); level <- max(0, min(1, s_late / struct_rate))
    }
    updateRadioButtons(session, "sim_acq_det_model", selected = model)
    updateNumericInput(session, "sim_struct_rate", value = round(struct_rate))
    updateNumericInput(session, "sim_bg_rate",     value = round(bg))
    updateSliderInput (session, "sim_acq_floor_frac", value = round(level, 2))
  }

  observeEvent(input$sim_real_file, {
    f <- input$sim_real_file
    req(f)
    ext <- get_ext(f$name)
    # Read to RAW, uncorrected z from every source, exactly like the SMLM tab —
    # then ONE axial factor is applied in sim_real_gated(), under user control.
    # Reading with the readers' 0.7 default instead would bake the correction in
    # for .npy/.json while leaving a .csv at whatever its producer applied, so
    # the same run would give a different sigma_z depending on its format.
    got <- tryCatch(
      withProgress(message = "Reading run…", value = 0.2, {
        switch(ext,
          npy  = list(df = read_minflux_npy_as_df(
                             f$datapath, keep_final_only = TRUE,
                             keep_valid_only = TRUE, z_scale = 1.0),
                      cal = list(factor = PYMINFLUX_Z_SCALE_DEFAULT,
                                 source = "raw instrument data (.npy)",
                                 certain = TRUE)),
          json = list(df = read_minflux_json_as_df(
                             f$datapath, keep_final_only = TRUE,
                             keep_valid_only = TRUE, z_scale = 1.0),
                      cal = list(factor = PYMINFLUX_Z_SCALE_DEFAULT,
                                 source = "raw instrument data (.json)",
                                 certain = TRUE)),
          # A CSV arrives with SOME factor already in its z. strip_csv_z_
          # calibration() divides it back out and reports what it found —
          # our own export records the factor explicitly, pyMINFLUX does not.
          csv  = {
            raw <- strip_csv_z_calibration(read.csv(f$datapath, check.names = FALSE))
            list(df = raw, cal = attr(raw, "z_cal"))
          },
          stop("Unsupported file type (use .npy, .json or .csv)."))
      }),
      error = function(e) {
        showNotification(paste("Load failed:", conditionMessage(e)),
                         type = "error", duration = 8); NULL
      })
    if (is.null(got)) return()
    # Same column canonicalization the SMLM tab applies, so alternative headers
    # (trace_id, x_nm, photons, …) are accepted here too.
    df <- canonicalize_minflux_df(got$df)
    if (!all(c("tid", "tim", "x", "y") %in% names(df))) {
      showNotification("Could not read run (need tid, tim, x, y).",
                       type = "error", duration = 6); return()
    }
    # Store the run ungated and at RAW z, together with the detected axial
    # provenance.  Building the widgets from it makes sim_real_gated() fire, and
    # the calibration observer does the rest — so upload, a gate change and an
    # axial-factor change all follow exactly the same path.
    sim_real_raw(list(df = df, cal = got$cal, ext = ext,
                      has_z = has_axial_data(df)))
    showNotification("Run loaded - calibrating from the current gate.",
                     type = "message", duration = 4)
  })

  # Push every measured value into the Acquisition / precision inputs.  Called on
  # upload AND on every gate change, so the inputs always reflect the subset the
  # gate currently selects.
  push_real_calibration <- function(a) {
    # Pre-fill the Acquisition inputs from the measured values.
    hrs  <- floor(a$duration_hr)
    mins <- min(59L, as.integer(round((a$duration_hr - hrs) * 60)))
    updateRadioButtons(session, "sim_mode", selected = "acquisition")
    updateNumericInput(session, "sim_acq_hours",     value = hrs)
    updateNumericInput(session, "sim_acq_minutes",   value = mins)
    updateNumericInput(session, "sim_acq_decay_min", value = round(a$tau_min, 1))
    updateNumericInput(session, "sim_fov_x", value = ceiling(a$fov_x * 1.1 / 100) * 100)
    updateNumericInput(session, "sim_fov_y", value = ceiling(a$fov_y * 1.1 / 100) * 100)
    if (a$fov_z > 1)
      updateNumericInput(session, "sim_fov_z",
                         value = ceiling(a$fov_z * 1.1 / 100) * 100)
    # Structure/background split from the spatial classification (overridable).
    updateSliderInput(session, "sim_real_struct_frac",
                      value = round(a$struct_frac, 2))
    apply_real_rates(a, a$struct_frac)

    # Measured localization precision → Manual σ; photons + effective L → CRLB.
    if (is.finite(a$sigma_x) && is.finite(a$sigma_y))
      updateNumericInput(session, "sim_precision_xy",
                         value = round(mean(c(a$sigma_x, a$sigma_y)), 1))
    if (is.finite(a$sigma_z))
      updateNumericInput(session, "sim_precision_z", value = round(a$sigma_z, 1))
    if (is.finite(a$photons))
      updateNumericInput(session, "sim_crlb_photons", value = round(a$photons))
    if (is.finite(a$L_eff))
      updateNumericInput(session, "sim_crlb_L", value = round(a$L_eff))
  }

  # Re-calibrate whenever the run or the gate changes.  Everything downstream
  # reads sim_real_rv(), so the whole tab follows the gate from this one place.
  # A gate that leaves too little to measure keeps the previous calibration
  # rather than silently resetting the inputs - the gate info box says so.
  observeEvent(sim_real_gated(), {
    g <- sim_real_gated()
    req(g)
    a <- analyze_real_run(g$df)
    if (is.null(a)) return()
    sim_real_rv(a)
    push_real_calibration(a)
  }, ignoreNULL = TRUE)

  # Manual override of the structure/background split (re-derives the rates).
  observeEvent(input$sim_real_struct_frac, {
    a <- sim_real_rv()
    if (!is.null(a)) apply_real_rates(a, input$sim_real_struct_frac)
  }, ignoreInit = TRUE)

  output$sim_real_status <- renderUI({
    a <- sim_real_rv()
    if (is.null(a))
      return(tags$small(style = "color:#888;",
        "Load a .npy / .json MINFLUX run to calibrate the model from real data."))
    # Counts here are those of the GATED subset, not the file - the gate box
    # directly above reports the before/after so the two can never be confused.
    # When the gate is starved, these numbers are the last calibration that
    # succeeded, so they must be labelled as stale rather than as current.
    g     <- sim_real_gated()
    stale <- !is.null(g) && (g$n_loc_out < 10L || g$n_tr_sigma < 3L)
    suffix <- if (stale) " (last valid gate)"
              else if (!is.null(g) && length(g$narrowed) > 0) " (gated)"
              else ""
    tags$div(class = if (stale) "alert alert-warning" else "alert alert-success",
             style = "padding:6px 10px; font-size:0.82em; margin:4px 0;",
      tags$b(sprintf("%d traces over %.1f h%s", a$n_traces, a$duration_hr, suffix)),
      tags$br(),
      sprintf("Total rate: peak %.0f, mean %.0f, floor %.0f /hr  (τ≈%.0f min)",
              a$rate_peak, a$rate_mean, a$rate_floor, a$tau_min),
      tags$br(),
      sprintf("Spatial split: %.0f%% on-structure / %.0f%% background",
              100 * a$struct_frac, 100 * (1 - a$struct_frac)),
      tags$br(),
      sprintf("FOV: %.0f × %.0f%s nm", a$fov_x, a$fov_y,
              if (a$fov_z > 1) sprintf(" × %.0f", a$fov_z) else ""),
      # ── Point-cloud profile ──────────────────────────────────────────────
      # Median ± nMAD per axis: the SPREAD is what makes the cloud MINFLUX-like,
      # so it is shown next to the median rather than collapsed into it.
      if (is.finite(a$sigma_x)) tagList(
        tags$hr(style = "margin:5px 0;"),
        tags$b("Point cloud"),
        sprintf(" (%d traces, ≥3 locs)", a$n_sigma_traces), tags$br(),
        sprintf("σ ± nMAD: x %.1f ± %.1f, y %.1f ± %.1f%s nm",
                a$sigma_x, a$sigma_x_mad, a$sigma_y, a$sigma_y_mad,
                if (is.finite(a$sigma_z))
                  sprintf(", z %.1f ± %.1f", a$sigma_z, a$sigma_z_mad) else ""),
        tags$br(),
        sprintf("Axis ratio: y/x %.2f%s", a$ratio_yx,
                if (is.finite(a$ratio_zx)) sprintf(", z/x %.2f", a$ratio_zx) else ""),
        tags$br(),
        sprintf("Locs per trace: median %.0f (IQR %.0f–%.0f)",
                a$len_med, a$len_q[1], a$len_q[2]),
        if (is.finite(a$photons)) tagList(tags$br(),
          sprintf("%.0f photons/loc (n=%d) → L ≈ %.0f nm", a$photons,
                  a$eco_n, a$L_eff)),
        tags$br(),
        tags$span(style = "color:#666;", sprintf(
          "σ is small-n debiased (×%.3f at the median trace length). Within-trace SD also contains real motion and drift, so it bounds precision rather than measuring it.",
          1 / a$debias_med))))
  })

  # ── Compare-strategies reactive ─────────────────────────────────────────────
  # Applies two labeling strategies to the same GT from sim_result().
  # Geometry is NOT re-run; seeds come from current label_seed + offsets.
  compare_result <- eventReactive(input$sim_compare_run, {
    res <- sim_result()
    shiny::validate(shiny::need(
      !is.null(res),
      "Run the main simulation first to generate a ground truth."))

    gt   <- res$gt
    use_3d <- res$use_3d
    sigma_loc_xy <- res$sigma_loc_xy
    sigma_loc_z  <- res$sigma_loc_z
    trace_length <- max(1, as.numeric(input$sim_trace_length))

    # Base seed: use the fixed label seed if set, else a fresh one per compare run
    base_seed <- if (isTRUE(input$sim_fix_label_seed)) {
      suppressWarnings(as.integer(input$sim_label_seed))
    } else {
      sample.int(.Machine$integer.max, 1L)
    }
    if (!is.finite(base_seed)) base_seed <- 1L

    cmp_protein_mw_kda <- max(0, as.numeric(input$sim_protein_mw %||% 0))
    cmp_protein_r_nm   <- if (cmp_protein_mw_kda > 0) 0.66 * cmp_protein_mw_kda^(1/3) else 0

    run_strategy <- function(label_key, eff_pct, seed_val, suffix) {
      set.seed(seed_val)
      if (identical(label_key, "Custom")) {
        n_fluor       <- max(1L, as.integer(input[[paste0("sim_cmp_", suffix, "_fluor")]] %||% 2L))
        link_nm       <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_link")]]  %||% 10))
        dye_spread_nm <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_spread")]] %||% 0))
        footprint_nm  <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_foot")]]  %||% 8))
        blink_lambda  <- max(0,  as.numeric(input$sim_blink_mean %||% 0.3))
        lbl_note      <- sprintf("Custom: %d fluor / %.0f nm linker",
                                 n_fluor, link_nm)
      } else {
        p             <- LABEL_PRESETS[[label_key]]
        n_fluor       <- p$n_fluor
        link_nm       <- p$link_nm
        dye_spread_nm <- p$dye_spread_nm
        footprint_nm  <- p$footprint_nm
        blink_lambda  <- p$blink_lambda %||% 0.3
        lbl_note      <- p$note
      }

      fluoro <- expand_to_fluorophores(
        gt, eff_pct, n_fluor, link_nm, dye_spread_nm,
        ldir_aware   = isTRUE(input$sim_ldir_aware),
        protein_r_nm = cmp_protein_r_nm)

      if (isTRUE(input$sim_blink_model))
        fluoro <- split_for_blinking(
          fluoro, max(0, as.numeric(input$sim_blink_mean %||% blink_lambda)))

      locs <- gen_localizations(fluoro, trace_length, sigma_loc_xy, sigma_loc_z,
                                 use_geometric  = isTRUE(input$sim_tl_geometric),
                                 scale_by_fluor = isTRUE(input$sim_tl_scale_fluor))

      # RMSE (capped at 500 TIDs for speed)
      rmse_vals <- NULL
      if ("mol_id" %in% names(fluoro) && nrow(gt) > 0 && nrow(locs) > 0) {
        n_samp <- min(nrow(fluoro), 500L)
        real_rows <- which(suppressWarnings(as.integer(locs$tid)) <= nrow(fluoro) &
                           !is.na(suppressWarnings(as.integer(locs$tid))))
        locs_by_tid <- split(real_rows, locs$tid[real_rows])
        errs <- vapply(seq_len(n_samp), function(k) {
          rows <- locs_by_tid[[as.character(k)]]
          if (is.null(rows) || length(rows) == 0L) return(NA_real_)
          lk <- locs[rows, , drop = FALSE]
          gi <- fluoro$mol_id[k]
          if (is.na(gi) || gi < 1L || gi > nrow(gt)) return(NA_real_)
          dx <- median(lk$x) - gt$x[gi]
          dy <- median(lk$y) - gt$y[gi]
          sqrt(dx^2 + dy^2)
        }, numeric(1L))
        rmse_vals <- errs[!is.na(errs)]
      }

      list(fluoro = fluoro, locs = locs, link_nm = link_nm,
           n_fluor = n_fluor, n_fluoro = nrow(fluoro),
           lbl_note = lbl_note, footprint_nm = footprint_nm,
           protein_r_nm = cmp_protein_r_nm,
           rmse_vals = rmse_vals)
    }

    eff_a <- max(1, min(100, as.numeric(input$sim_cmp_eff_a %||% 80)))
    eff_b <- max(1, min(100, as.numeric(input$sim_cmp_eff_b %||% 80)))
    a <- run_strategy(input$sim_cmp_label_a %||%
                        "Direct tag (SNAP / Halo / Click, DOL=1)",
                      eff_a, base_seed, "a")
    b <- run_strategy(input$sim_cmp_label_b %||%
                        "Primary monoclonal IgG + secondary IgG",
                      eff_b, base_seed + 1L, "b")

    list(gt = gt, use_3d = use_3d, a = a, b = b,
         sigma_loc_xy = sigma_loc_xy,
         lbl_a = input$sim_cmp_label_a %||%
                   "Direct tag (SNAP / Halo / Click, DOL=1)",
         lbl_b = input$sim_cmp_label_b %||%
                   "Primary monoclonal IgG + secondary IgG")
  }, ignoreNULL = TRUE)

  # Per-localization local density: number of neighbours within radius_nm
  # (self excluded). Used by the "Local density" colour mode. k is capped so
  # very dense scenes stay fast; counts saturate at k-1 in that regime.
  sim_local_density <- function(locs, radius_nm) {
    cols <- intersect(c("x", "y", "z"), names(locs))
    pts  <- as.matrix(locs[, cols, drop = FALSE])
    n    <- nrow(pts)
    if (n < 2) return(rep(0, n))
    k    <- min(n, 256L)
    nn   <- RANN::nn2(pts, pts, k = k)
    rowSums(nn$nn.dists <= radius_nm) - 1L
  }

  # ── Phase 1 plot ─────────────────────────────────────────────────────────────
  output$sim_plot <- renderPlotly({
    res <- sim_result()
    shiny::validate(
      shiny::need(!is.null(res), "Click 'Generate' to run the simulation."))

    dark    <- isTRUE(input$sim_dark_mode)
    bg_col  <- if (dark) "#111111" else "#FFFFFF"
    fg_col  <- if (dark) "#DDDDDD" else "#444444"
    use_3d  <- isTRUE(res$use_3d)
    tr_type <- if (use_3d) "scatter3d" else "scattergl"
    pt_size <- max(1L, as.integer(input$sim_point_size))
    gc      <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"
    show_grid <- !isTRUE(input$sim_hide_grid)

    # Field of view - display window only, nothing is cropped. Axes are pinned to
    # a real-camera-style FOV; the 3D aspect ratio reflects the true X:Y:Z
    # proportions so Z appears compressed (4000×4000×1000). Always on in
    # Acquisition mode; in Normal mode it's opt-in via "Pin axes to a FOV (ROI)".
    use_fov <- sim_use_fov()
    fov_x   <- max(1, as.numeric(input$sim_fov_x %||% 4000))
    fov_y   <- max(1, as.numeric(input$sim_fov_y %||% 4000))
    fov_z   <- max(1, as.numeric(input$sim_fov_z %||% 1000))

    plt <- plot_ly(source = "sim2d")

    # Helper: add a plain marker layer
    add_layer <- function(p, df, symbol, size, color, name, htxt = NULL) {
      if (nrow(df) == 0) return(p)
      z_val <- if (use_3d && "z" %in% names(df)) df$z else NULL
      args  <- list(
        x = df$x, y = df$y,
        type = tr_type, mode = "markers",
        marker = list(symbol = symbol, size = size, color = color,
                      line = list(color = color, width = 1)),
        name = name,
        hoverinfo = if (is.null(htxt)) "skip" else "text",
        text = htxt
      )
      if (!is.null(z_val)) args$z <- z_val
      do.call(add_trace, c(list(p), args))
    }

    # Build hover text helpers
    hover_xy <- function(df, prefix = "") {
      z_part <- if (use_3d && "z" %in% names(df))
                  paste0("<br>z: ", round(df$z, 1), " nm")
                else rep("", nrow(df))
      paste0(prefix,
             "x: ", round(df$x, 1), " nm",
             "<br>y: ", round(df$y, 1), " nm",
             z_part)
    }

    # GT layer - bright distinct markers so they stand out in both 2D and 3D.
    # In image mode these are the intensity-weighted sampled molecule positions.
    if (isTRUE(input$sim_show_gt) && nrow(res$gt) > 0) {
      gt_clean <- res$gt[, intersect(c("x","y","z"), names(res$gt)), drop = FALSE]
      gt_label <- if (identical(input$sim_shape, "image"))
                    "GT (sampled positions)" else "Ground truth"
      ht <- hover_xy(gt_clean, "GT<br>")
      plt <- add_layer(plt, gt_clean,
                       if (use_3d) "circle-open" else "cross",
                       round(pt_size * 2),
                       if (dark) "#FFEE44" else "#CC2200",
                       gt_label, ht)
    }


    # Localization layer - coloured by TID, Z, time, or local density.
    if (isTRUE(input$sim_show_locs) && nrow(res$locs) > 0) {
      cmode <- input$sim_color_mode %||% "tid"
      if (identical(cmode, "z") && !(use_3d && "z" %in% names(res$locs)))
        cmode <- "tid"                                  # no depth in 2D \u2192 TID
      z_val <- if (use_3d && "z" %in% names(res$locs)) res$locs$z else NULL
      z_ht  <- if (!is.null(z_val))
                  paste0("<br>z: ", round(z_val, 1), " nm")
               else rep("", nrow(res$locs))
      ht <- paste0("TID: ", res$locs$tid,
                   "<br>x: ", round(res$locs$x, 1), " nm",
                   "<br>y: ", round(res$locs$y, 1), " nm",
                   z_ht)

      if (identical(cmode, "tid")) {
        pal    <- tid_color_map(unique(res$locs$tid))
        marker <- list(size = pt_size,
                       color = unname(pal[res$locs$tid]), opacity = 0.85)
      } else {
        cv <- switch(cmode,
          z       = list(v = res$locs$z,   title = "Depth z (nm)"),
          tim     = list(v = res$locs$tim, title = "Time (s)"),
          locdens = list(v = sim_local_density(res$locs,
                              max(1, as.numeric(input$sim_dens_radius_nm %||% 500))),
                         title = "Local density"))
        v   <- as.numeric(cv$v)
        rng <- range(v[is.finite(v)])
        if (!all(is.finite(rng))) rng <- c(0, 1)
        marker <- list(
          size = pt_size, color = v, opacity = 0.85,
          colorscale = plotly_colorscale(
            input$sim_numeric_palette %||% "Spectral",
            reverse = isTRUE(input$sim_palette_reverse)),
          cmin = rng[1], cmax = rng[2], showscale = TRUE,
          colorbar = list(title = cv$title,
                          titlefont = list(color = fg_col),
                          tickfont  = list(color = fg_col)))
      }
      args <- list(
        x = res$locs$x, y = res$locs$y,
        type = tr_type, mode = "markers",
        marker = marker,
        name = sprintf("Localizations (n\u202F=\u202F%d)", nrow(res$locs)),
        hoverinfo = "text", text = ht
      )
      if (!is.null(z_val)) args$z <- z_val
      plt <- do.call(add_trace, c(list(plt), args))
    }

    # Lines: GT molecule position → median of each TID.
    #
    # One line per TID, anchored at the true GT molecule the TID came from.
    # For single-dye labels one short line per molecule shows the small
    # link offset; for polyclonal + secondary many longer lines fan out
    # from the same GT point, reflecting the greater linker distance and
    # the number of fluorophores.  Works in all shape modes including image.
    if (isTRUE(input$sim_show_gt_lines) &&
        nrow(res$gt) > 0 && nrow(res$fluoro) > 0 && nrow(res$locs) > 0 &&
        "mol_id" %in% names(res$fluoro)) {
      mol_ids  <- res$fluoro$mol_id
      has_gz   <- use_3d && "z" %in% names(res$gt)
      has_lz   <- use_3d && "z" %in% names(res$locs)
      seg_col  <- if (dark) "rgba(255,120,80,0.75)" else "rgba(190,50,30,0.75)"
      n_fluoro <- nrow(res$fluoro)

      # One segment per TID: GT[mol_id] → median(locs for that TID)
      segs <- do.call(rbind, lapply(seq_len(n_fluoro), function(k) {
        lk <- res$locs[res$locs$tid == as.character(k), , drop = FALSE]
        if (nrow(lk) == 0) return(NULL)
        gi <- mol_ids[k]
        if (is.na(gi) || gi < 1L || gi > nrow(res$gt)) return(NULL)
        data.frame(
          x0 = res$gt$x[gi],  y0 = res$gt$y[gi],
          x1 = median(lk$x),  y1 = median(lk$y),
          z0 = if (has_gz) res$gt$z[gi]  else NA_real_,
          z1 = if (has_lz) median(lk$z)  else NA_real_
        )
      }))

      if (!is.null(segs) && nrow(segs) > 0) {
        if (use_3d) {
          xs <- c(rbind(segs$x0, segs$x1, NA))
          ys <- c(rbind(segs$y0, segs$y1, NA))
          zs <- c(rbind(segs$z0, segs$z1, NA))
          plt <- do.call(add_trace,
                         list(plt, x = xs, y = ys, z = zs,
                              type = "scatter3d", mode = "lines",
                              line = list(color = seg_col, width = 2),
                              name = "GT \u2192 TID median",
                              hoverinfo = "skip", showlegend = TRUE))
        } else {
          plt <- plt %>%
            add_segments(x = ~x0, xend = ~x1, y = ~y0, yend = ~y1,
                         data = segs, inherit = FALSE,
                         line = list(color = seg_col, width = 1.5),
                         name = "GT \u2192 TID median",
                         hoverinfo = "skip", showlegend = TRUE)
        }
      }
    }

    title_txt <- sprintf("Simulated localizations%s",
                         if (use_3d) " (3D)" else "")

    if (use_3d) {
      ref_pts <- if (nrow(res$locs) > 0) res$locs else res$gt
      if (use_fov) {
        # FOV pins each axis; aspect ratio = true X:Y:Z so Z is compressed.
        xr <- c(-fov_x / 2, fov_x / 2)
        yr <- c(-fov_y / 2, fov_y / 2)
        zr <- c(-fov_z / 2, fov_z / 2)
        aspect <- list(x = 1, y = fov_y / fov_x, z = fov_z / fov_x)
      } else {
        # Equal-scale axes: Z gets the same span as XY so the plot looks cubic.
        # Localizations remain at their actual Z positions near the centre.
        x_rng <- range(ref_pts$x, na.rm = TRUE)
        y_rng <- range(ref_pts$y, na.rm = TRUE)
        half  <- max(diff(x_rng), diff(y_rng)) / 2 * 1.05   # 5% padding
        if (half < 1) half <- 1
        x_ctr <- mean(x_rng); y_ctr <- mean(y_rng)
        z_ctr <- if ("z" %in% names(ref_pts))
                    median(ref_pts$z, na.rm = TRUE) else 0
        xr <- c(x_ctr - half, x_ctr + half)
        yr <- c(y_ctr - half, y_ctr + half)
        zr <- c(z_ctr - half, z_ctr + half)
        aspect <- list(x = 1, y = 1, z = 1)
      }

      plt <- plt %>% layout(
        paper_bgcolor = bg_col,
        scene = list(
          bgcolor = bg_col,
          xaxis = list(title = "X (nm)", color = fg_col, range = xr,
                       showgrid = show_grid, zeroline = show_grid),
          yaxis = list(title = "Y (nm)", color = fg_col, range = yr,
                       showgrid = show_grid, zeroline = show_grid),
          zaxis = list(title = "Z (nm)", color = fg_col, range = zr,
                       showgrid = show_grid, zeroline = show_grid),
          aspectmode  = "manual",
          aspectratio = aspect
        ),
        font = list(color = fg_col), showlegend = TRUE,
        legend = list(x = 0.01, y = 0.99, xanchor = "left", yanchor = "top",
                      bgcolor = if (dark) "rgba(17,17,17,0.6)"
                                else "rgba(255,255,255,0.6)",
                      font = list(size = 11)),
        title = list(text = title_txt, font = list(color = fg_col))
      )
    } else {
      img_layout <- if (!is.null(res$img_b64)) {
        list(list(source  = res$img_b64,
                  xref    = "x", yref = "y",
                  x       = -res$img_w_nm / 2,
                  y       = res$img_h_nm / 2,
                  sizex   = res$img_w_nm, sizey = res$img_h_nm,
                  sizing  = "stretch", opacity = 0.25, layer = "below"))
      } else NULL

      xaxis <- list(title = "X (nm)", color = fg_col,
                    scaleanchor = "y", scaleratio = 1, gridcolor = gc,
                    showgrid = show_grid, zeroline = show_grid)
      yaxis <- list(title = "Y (nm)", color = fg_col, gridcolor = gc,
                    showgrid = show_grid, zeroline = show_grid)
      if (use_fov) {
        xaxis$range <- c(-fov_x / 2, fov_x / 2)
        yaxis$range <- c(-fov_y / 2, fov_y / 2)
      }

      # Scale bar (2D only) - reuses the shared SMLM make_scale_bar helper.
      sb <- NULL
      if (isTRUE(input$sim_show_scale_bar)) {
        x_span <- if (use_fov) fov_x else {
          ref <- if (nrow(res$locs) > 0) res$locs else res$gt
          diff(range(ref$x, na.rm = TRUE))
        }
        sb <- make_scale_bar(
          suppressWarnings(as.numeric(input$sim_scale_bar_length)),
          x_span,
          suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
          dark = dark,
          scale_font_size = suppressWarnings(
            as.numeric(input$sim_scale_font_size %||% 18)))
      }

      lay_args <- list(plt,
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font   = list(color = fg_col),
        xaxis  = xaxis, yaxis = yaxis,
        showlegend = TRUE,
        legend = list(x = 0.01, y = 0.99, xanchor = "left", yanchor = "top",
                      bgcolor = if (dark) "rgba(17,17,17,0.6)"
                                else "rgba(255,255,255,0.6)",
                      font = list(size = 11)),
        title  = list(text = title_txt, font = list(color = fg_col)),
        images = img_layout,
        uirevision = "sim"
      )
      if (!is.null(sb)) {
        lay_args$shapes      <- sb$shapes
        lay_args$annotations <- sb$annotations
      }
      plt <- do.call(layout, lay_args)
    }
    plt <- plt %>% config(displaylogo = FALSE,
                   toImageButtonOptions = list(format = "svg",
                                               filename = "sim_localisations"))
    # Register the relayout event so the "sim2d" source is recognised in both
    # 2D and 3D (avoids the "source not registered" console warning). Flag it as
    # rendered so the relayout observer only reads the event after registration.
    sim_plot_ready(TRUE)
    plt %>% plotly::event_register("plotly_relayout")
  })

  # Capture the visible x-span on zoom/pan and re-push the scale bar so it
  # tracks the current view (matches the SMLM tab's behaviour).
  observe({
    # Only read the relayout event when the plot has rendered + registered AND is
    # on screen (sim tab + MINFLUX subtab + 2D); otherwise plotly warns "source
    # not registered". sim_plot_ready stays FALSE until the first Generate render.
    if (!isTRUE(sim_plot_ready()) ||
        !identical(input$main_navbar, "tab_sim") ||
        !identical(input$sim_main_tabs, "MINFLUX") ||
        !identical(input$sim_view_mode, "2d")) return()
    ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "sim2d"))
    if (!is.list(ed)) return()
    x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
    x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
    if (length(x0) == 1 && length(x1) == 1 &&
        is.finite(x0) && is.finite(x1) && x1 > x0) {
      sim_axis_x_span(x1 - x0)
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      span <- sim_full_x_span()
      if (length(span) == 1 && is.finite(span) && span > 0) sim_axis_x_span(span)
    } else return()
    if (isTRUE(input$sim_show_scale_bar)) sim_push_scale_bar()
  })

  # Confocal view: capture its own zoom and re-push its (white) scale bar.
  observe({
    if (!isTRUE(sim_diff_ready()) ||
        !identical(input$main_navbar, "tab_sim") ||
        !identical(input$sim_main_tabs, "Confocal view")) return()
    ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "sim_diff"))
    if (!is.list(ed)) return()
    x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
    x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
    if (length(x0) == 1 && length(x1) == 1 &&
        is.finite(x0) && is.finite(x1) && x1 > x0) {
      sim_diff_x_span(x1 - x0)
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      sim_diff_x_span(NULL)
    } else return()
    if (isTRUE(input$sim_show_scale_bar)) sim_push_diff_scale_bar()
  })

  observeEvent(list(input$sim_show_scale_bar, input$sim_scale_bar_length,
                    input$sim_scale_bar_thickness, input$sim_scale_font_size,
                    input$sim_dark_mode, input$sim_view_mode,
                    input$sim_mode, input$sim_view_fov, input$sim_fov_x), {
    sim_push_scale_bar()
    sim_push_diff_scale_bar()
  }, ignoreInit = TRUE)

  # ── Phase 1 summary ──────────────────────────────────────────────────────────
  output$sim_summary <- renderUI({
    res <- sim_result()
    if (is.null(res)) return(NULL)
    nn_dist <- if (nrow(res$locs) >= 2) {
      idx  <- if (nrow(res$locs) > 500) sample(nrow(res$locs), 500)
              else seq_len(nrow(res$locs))
      cols <- intersect(c("x","y","z"), names(res$locs))
      pts  <- as.matrix(res$locs[idx, cols])
      d    <- as.matrix(dist(pts)); diag(d) <- Inf
      mean(apply(d, 1, min), na.rm = TRUE)
    } else NA_real_

    row <- function(lbl, val)
      tags$tr(tags$td(style = "padding:2px 8px;", lbl),
              tags$td(style = "padding:2px 8px;", tags$b(val)))

    # Steric warning block
    steric_ui <- if (!is.null(res$steric_warn)) {
      alert_cls <- if (identical(res$steric_warn$level, "danger"))
                     "alert alert-danger" else "alert alert-warning"
      tags$div(class = alert_cls,
               style = "padding:6px 10px; font-size:0.84em; margin-top:6px;",
               tags$b(if (identical(res$steric_warn$level, "danger"))
                        "\u26A0 Severe steric exclusion" else "\u26A0 Steric exclusion warning"),
               tags$br(), res$steric_warn$msg)
    } else NULL

    # Blinking TID inflation row (shown only when blinking is active)
    blink_row <- if (isTRUE(res$blink_active) && !is.null(res$n_fluoro_base)) {
      row("Blinking TID inflation",
          sprintf("%d base \u2192 %d apparent  (\u00D7%.1f, \u03BB\u202F=\u202F%.2f)",
                  res$n_fluoro_base, res$n_fluoro,
                  res$n_fluoro / max(1, res$n_fluoro_base),
                  res$blink_lambda_used))
    } else NULL

    tagList(
      steric_ui,
      tags$div(
        class = "well well-sm", style = "font-size:0.88em; margin-top:10px;",
        tags$b("Simulation summary"), tags$br(),
        tags$table(style = "width:100%; border-collapse:collapse;",
          row("Molecules placed",  res$n_prot),
          row("Labeling strategy", res$lbl_note),
          row("Fluorophores (trace IDs)",
              sprintf("%d  (~%.1f/molecule)", res$n_fluoro,
                      res$n_fluoro / max(1, res$n_prot))),
          blink_row,
          row("Linker distance",   sprintf("~%.0f nm", res$link_nm)),
          if (!is.na(res$footprint_nm))
            row("Label footprint (steric check)", sprintf("~%.0f nm", res$footprint_nm)),
          row("Total localizations", res$n_locs),
          row("Avg locs / trace ID",
              sprintf("%.1f", res$n_locs / max(1, res$n_fluoro))),
          row("Precision XY \u03C3",    sprintf("%.1f nm", res$sigma_loc_xy)),
          if (isTRUE(res$use_3d))
            row("Precision Z \u03C3",   sprintf("%.1f nm", res$sigma_loc_z)),
          if (is.finite(nn_dist))
            row("Mean nearest-neighbour", sprintf("%.1f nm", nn_dist)),
          if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
            row("Mean TID displacement error",
                sprintf("%.1f nm (SD %.1f, n\u202F=\u202F%d)",
                        mean(res$rmse_vals), sd(res$rmse_vals),
                        length(res$rmse_vals))),
          if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
            row("Theoretical combined \u03C3  (\u221A(link\u00B2 + \u03C3\u00B2_xy))",
                sprintf("%.1f nm", sqrt(res$link_nm^2 + res$sigma_loc_xy^2)))
        )
      )
    )
  })

  # ── Labeling strategy info card (live, no Generate needed) ─────────────────
  output$sim_label_info <- renderUI({
    key <- input$sim_label_type %||% "Custom"
    if (identical(key, "Custom")) return(NULL)
    p <- LABEL_PRESETS[[key]]
    if (is.null(p) || is.na(p$n_fluor)) return(NULL)
    tags$div(
      class = "alert alert-info",
      style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
      tags$b(sprintf("%d fluorophore(s) per molecule", p$n_fluor)),
      sprintf(" - linker ~%d nm", round(p$link_nm)),
      tags$br(), p$note
    )
  })

  # ── CRLB precision preview (live - no Generate needed) ──────────────────────
  output$sim_crlb_preview <- renderUI({
    req(identical(input$sim_prec_mode %||% "manual", "crlb"))
    L_nm    <- max(1,    as.numeric(input$sim_crlb_L       %||% 100))
    n_ph    <- max(1,    as.numeric(input$sim_crlb_photons %||% 50))
    bg_frac <- max(0, min(0.99, as.numeric(input$sim_crlb_bg %||% 0.05)))
    n_eff   <- max(1, n_ph * (1 - bg_frac))
    sig_xy  <- max(0.1, L_nm / (2 * sqrt(n_eff)))
    use_3d  <- identical(input$sim_view_mode, "3d")
    z_fac   <- max(0.1, as.numeric(input$sim_crlb_z_factor %||% 0.7))
    sig_z   <- sig_xy * z_fac
    tags$div(
      class = "alert alert-success",
      style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
      tags$b(sprintf("\u03C3_xy \u2248 %.2f nm", sig_xy)),
      if (use_3d) sprintf("   |   \u03C3_z \u2248 %.2f nm", sig_z),
      tags$br(),
      tags$small(sprintf(
        "L\u202F=\u202F%.0f nm,  N_eff\u202F=\u202F%.0f ph",
        L_nm, n_eff))
    )
  })

  # ── Phase 1 per-molecule TID count summary ───────────────────────────────────
  output$sim_tid_summary <- renderUI({
    res <- sim_result()
    if (is.null(res) || !("mol_id" %in% names(res$fluoro)) ||
        nrow(res$fluoro) == 0) return(NULL)

    # Count detected TIDs per molecule (background TIDs have mol_id = NA)
    real_tids  <- res$fluoro[!is.na(res$fluoro$mol_id), , drop = FALSE]
    if (nrow(real_tids) == 0) return(NULL)
    tids_per_mol <- as.integer(table(real_tids$mol_id))
    # Molecules that lost all fluorophores to detection efficiency = 0 TIDs
    n_zero  <- res$n_prot - length(tids_per_mol)
    counts  <- c(tids_per_mol, rep(0L, max(0L, n_zero)))

    row <- function(lbl, val)
      tags$tr(tags$td(style = "padding:2px 8px;", lbl),
              tags$td(style = "padding:2px 8px;", tags$b(val)))

    tagList(tags$div(
      class = "well well-sm",
      style = "font-size:0.88em; margin-top:6px;",
      tags$b("TIDs per molecule"),
      tags$table(style = "width:100%; border-collapse:collapse;",
        row("Mean",    sprintf("%.2f", mean(counts))),
        row("SD",      sprintf("%.2f", sd(counts))),
        row("Median",  sprintf("%d",   as.integer(median(counts)))),
        row("Range",   sprintf("%d \u2013 %d", min(counts), max(counts))),
        row("Molecules with 0 TIDs",
            sprintf("%d (%.0f%%)", n_zero,
                    100 * n_zero / max(1, res$n_prot))),
        if (res$n_bg > 0)
          row("Background TIDs added", res$n_bg)
      )
    ))
  })

  # ── Phase 1 precision budget ─────────────────────────────────────────────────
  # Bar chart breaking down: linker error | localization precision | combined
  # (+ observed RMSE when GT data are available).
  output$sim_precision_budget <- renderUI({
    res <- sim_result()
    if (is.null(res)) return(NULL)

    # ── Error budget, all terms as RADIAL RMS displacement (nm) ───────────────
    # Every component must be expressed in the SAME quantity before it can go in
    # quadrature.  The observed value this is compared against is a radial
    # distance |trace median - GT|, so radial RMS is the common currency.
    #
    #   protein   epitope offset: a FIXED magnitude r        -> r^2
    #   linker    |N(0, s)| with s = link*sqrt(pi/2), so the MEAN magnitude is
    #             link_nm but E[|d|^2] = s^2 = (pi/2)*link^2.  Using link^2 here
    #             (as this panel used to) understates the term by 20 %.
    #   dye       per-AXIS sigma -> radial variance is D * dye_spread^2
    #   loc       per-AXIS sigma, and the TID position is the MEDIAN of n_loc
    #             localizations, so the localization term is DIVIDED DOWN by the
    #             trace length.  Omitting that made "observed vs theoretical"
    #             compare two different quantities.
    #             The textbook Var(median) = (pi/2)*sigma^2/n is ASYMPTOTIC and
    #             runs +17 % at n=3, +9 % at n=5 - exactly the short traces
    #             MINFLUX produces.  (pi/2)/(n + pi/2 - 1) holds to ~2 % from
    #             n=3 up.  It smooths over the odd/even-n zigzag in the exact
    #             value, which is right here because Poisson trace lengths mix
    #             both parities anyway.
    D       <- if (isTRUE(res$use_3d)) 3L else 2L
    link    <- res$link_nm
    dspread <- res$dye_spread_nm %||% 0
    prot_r  <- res$protein_r_nm  %||% 0
    sig_xy  <- res$sigma_loc_xy
    sig_z   <- res$sigma_loc_z %||% sig_xy

    # Localizations per REAL trace (background TIDs excluded).  The WHOLE vector
    # is kept, not just its median: trace lengths are Poisson/geometric and the
    # median-variance factor is convex, so evaluating it at the median length
    # under-states the ensemble mean by 5-7 % (Jensen).  Averaging the factor
    # over the realized lengths removes that exactly.
    n_vec     <- numeric(0)
    n_loc_med <- 1
    if (!is.null(res$locs) && nrow(res$locs) > 0 && !is.null(res$fluoro)) {
      tid_i <- suppressWarnings(as.integer(res$locs$tid))
      real  <- !is.na(tid_i) & tid_i >= 1L & tid_i <= nrow(res$fluoro)
      if (any(real)) {
        n_vec     <- as.numeric(table(tid_i[real]))
        n_loc_med <- max(1, stats::median(n_vec))
      }
    }

    # Variance of the median of n Gaussian draws, per unit sigma^2.
    med_var_f <- function(n) (pi / 2) / (pmax(1, n) + pi / 2 - 1)
    var_loc1  <- if (D == 3L) 2 * sig_xy^2 + sig_z^2 else 2 * sig_xy^2
    v_prot    <- prot_r^2
    v_link    <- (pi / 2) * link^2
    v_dye     <- D * dspread^2
    v_loc     <- var_loc1 * if (length(n_vec)) mean(med_var_f(n_vec))
                            else med_var_f(n_loc_med)

    theory <- sqrt(v_prot + v_link + v_dye + v_loc)
    # What NO amount of acquisition time removes: the per-molecule label offset.
    # dye_spread is excluded because it averages once a molecule's TIDs are
    # merged; protein + linker are drawn once per molecule and never do.
    lab_floor <- sqrt(v_prot + v_link)

    # Observed as an RMS too (it was a mean of magnitudes, which is a different
    # statistic and always smaller), so the two are directly comparable.
    obs <- if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
             sqrt(mean(res$rmse_vals^2)) else NA_real_

    # Per-MOLECULE error after merging that molecule's TIDs — the quantity that
    # decides whether two molecules can be told apart.  dye_spread averages over
    # the dyes detected, sigma_loc over all their localizations.
    n_dye_det <- max(1, res$n_fluor * (res$eff_pct %||% 100) / 100)
    n_tid_mol <- n_dye_det * (1 + (res$blink_lambda %||% 0))
    n_loc_mean <- if (length(n_vec)) mean(n_vec) else n_loc_med
    eps_mol   <- sqrt(v_prot + v_link + D * dspread^2 / n_dye_det +
                      var_loc1 * med_var_f(n_tid_mol * n_loc_mean))

    # Two molecules are separable at ~k x their position error. k = 2 is kept
    # for continuity with previous numbers; strictly, two INDEPENDENT estimates
    # each of error eps separate at 2*sqrt(2)*eps.
    resolve_limit <- 2 * eps_mol
    max_val       <- max(theory, obs, 1, na.rm = TRUE) * 1.05

    make_bar <- function(label, value, color) {
      if (is.na(value) || value < 0) return(NULL)
      w_pct <- sprintf("%.1f%%", 100 * min(value, max_val) / max_val)
      tags$div(style = "margin-bottom:7px;",
        tags$div(
          style = paste0("font-size:0.80em; margin-bottom:2px;",
                         " display:flex; justify-content:space-between;"),
          tags$span(label),
          tags$b(sprintf("%.1f nm", value))
        ),
        tags$div(
          style = "background:#e8e8e8; border-radius:3px; height:11px;",
          tags$div(style = sprintf(
            "background:%s; width:%s; height:11px; border-radius:3px;",
            color, w_pct))
        )
      )
    }

    tags$div(
      class = "well well-sm",
      style = "font-size:0.88em; margin-top:8px;",
      tags$b("Precision Budget"),
      tags$br(),
      tags$small(sprintf(
        "Radial RMS displacement of one TID from its molecule (%dD, %.0f locs/trace):",
        D, n_loc_med)),
      tags$br(), tags$br(),
      if (prot_r > 0) make_bar("Protein size (epitope offset)", prot_r, "#be185d"),
      make_bar("Linker / label error  (\u221A(\u03C0/2)\u00B7link)",
               sqrt(v_link), "#e07b39"),
      if (dspread > 0)
        make_bar("Dye spread in complex", sqrt(v_dye), "#0891b2"),
      make_bar("Localization  (\u03c3 averaged over the trace)",
               sqrt(v_loc), "#3b82c4"),
      make_bar("Combined - theoretical", theory, "#7c3aed"),
      if (!is.na(obs))
        make_bar(sprintf("Observed TID displacement (RMS, n\u202F=\u202F%d)",
                         length(res$rmse_vals)),
                 obs, "#16a34a"),
      tags$small(style = "color:#888;",
        "Observed and theoretical should agree; both are radial RMS."),
      # The label floor is the headline: it is what the acquisition CANNOT fix.
      tags$div(
        style = paste0("margin-top:8px; padding:5px 8px; border-radius:4px;",
                       " font-size:0.82em; background:#f1f5f9; color:#334155;"),
        tags$b(sprintf("Label floor: %.1f nm", lab_floor)),
        tags$br(),
        "Protein offset + linker are drawn once per molecule, so no extra ",
        "localizations or acquisition time reduce them."),
      if (!is.na(res$gt_nnd)) {
        resolvable <- res$gt_nnd >= resolve_limit
        never      <- res$gt_nnd < 2 * lab_floor
        tags$div(
          style = sprintf(
            paste0("margin-top:8px; padding:5px 8px; border-radius:4px;",
                   " font-size:0.82em; background:%s; color:%s;"),
            if (resolvable) "#dcfce7" else "#fef9c3",
            if (resolvable) "#166534" else "#854d0e"),
          tags$b(if (resolvable) "\u2714 Structure resolvable"
                 else if (never) "\u2716 Not resolvable with this label"
                 else "\u26A0 Below resolution limit"),
          tags$br(),
          sprintf(paste0("GT spacing \u2248 %.0f nm vs 2\u00D7 per-molecule error",
                         " (%.1f nm) = %.0f nm.  %s"),
            res$gt_nnd, eps_mol, resolve_limit,
            if (resolvable)
              "Individual molecules can be resolved."
            else if (never)
              paste0("The spacing is below 2\u00D7 the label floor, so this ",
                     "cannot be fixed by imaging longer - a shorter linker ",
                     "is required.")
            else
              paste0("Currently unresolved, but the label floor still allows ",
                     "it - more localizations per molecule would close the gap."))
        )
      }
    )
  })

  # ── Compare strategies: helper to build a stripped-down scatter plot ─────────
  make_cmp_plot <- function(gt, locs, use_3d, title,
                             loc_color = "rgba(220,38,38,0.65)",
                             dark = FALSE) {
    bg_col <- if (dark) "#111111" else "#FFFFFF"
    fg_col <- if (dark) "#DDDDDD" else "#444444"
    tr_type <- if (use_3d) "scatter3d" else "scattergl"
    show_grid <- !isTRUE(input$sim_hide_grid)

    # Compute isotropic axis ranges from GT so both panels share identical scales
    # and 1 nm maps to the same pixel distance on every axis.
    has_z  <- use_3d && "z" %in% names(gt)
    cx     <- mean(range(gt$x, na.rm = TRUE))
    cy     <- mean(range(gt$y, na.rm = TRUE))
    cz     <- if (has_z) mean(range(gt$z, na.rm = TRUE)) else 0
    spans  <- c(diff(range(gt$x, na.rm = TRUE)),
                diff(range(gt$y, na.rm = TRUE)),
                if (has_z) diff(range(gt$z, na.rm = TRUE)) else 0)
    half   <- max(spans, 1) / 2 * 1.15   # 15 % padding, same on all axes
    xrng   <- cx + c(-half, half)
    yrng   <- cy + c(-half, half)
    zrng   <- cz + c(-half, half)

    gt_col <- if (dark) "#FFFFFF" else "#000000"

    # Draw localizations first so GT markers render on top of them.
    plt <- plot_ly() %>%
      add_trace(type = tr_type,
                x = locs$x, y = locs$y,
                z = if (use_3d && "z" %in% names(locs)) locs$z else NULL,
                mode = "markers",
                marker = list(size = if (use_3d) 4 else 7,
                              color = loc_color,
                              opacity = 0.75),
                name = "Localizations",
                hoverinfo = "skip")

    # GT layer drawn on top so it is never hidden by localizations.
    # All shapes use open circle markers (no connecting lines).
    plt <- add_trace(plt,
              type = tr_type,
              x = gt$x, y = gt$y,
              z = if (has_z) gt$z else NULL,
              mode = "markers",
              marker = list(size   = if (use_3d) 10 else 14,
                            color  = gt_col,
                            symbol = "circle-open"),
              name = "GT molecules",
              hoverinfo = "skip")

    plt <- plt %>%
      layout(
        title = list(text = title, font = list(size = 12, color = fg_col)),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        xaxis = list(color = fg_col, title = "x (nm)",
                     scaleanchor = "y", scaleratio = 1,
                     range = xrng, constrain = "domain", showgrid = show_grid, zeroline = show_grid),
        yaxis = list(color = fg_col, title = "y (nm)",
                     range = yrng, constrain = "domain", showgrid = show_grid, zeroline = show_grid),
        showlegend = TRUE,
        legend = list(font = list(size = 10, color = fg_col))
      )
    if (use_3d) plt <- plt %>%
      layout(scene = list(
        aspectmode  = "manual",
        aspectratio = list(x = 1, y = 1, z = 1),
        xaxis = list(title = "x (nm)", color = fg_col, range = xrng,
                     showgrid = show_grid, zeroline = show_grid),
        yaxis = list(title = "y (nm)", color = fg_col, range = yrng,
                     showgrid = show_grid, zeroline = show_grid),
        zaxis = list(title = "z (nm)", color = fg_col, range = zrng,
                     showgrid = show_grid, zeroline = show_grid),
        bgcolor = bg_col))
    plt
  }

  output$sim_cmp_plot_a <- renderPlotly({
    cmp <- compare_result()
    shiny::validate(shiny::need(!is.null(cmp), "Click Compare."))
    make_cmp_plot(cmp$gt, cmp$a$locs, cmp$use_3d,
                  "Strategy A", "rgba(220,38,38,0.75)",
                  dark = isTRUE(input$sim_dark_mode))
  })

  output$sim_cmp_plot_b <- renderPlotly({
    cmp <- compare_result()
    shiny::validate(shiny::need(!is.null(cmp), "Click Compare."))
    make_cmp_plot(cmp$gt, cmp$b$locs, cmp$use_3d,
                  "Strategy B", "rgba(37,99,235,0.75)",
                  dark = isTRUE(input$sim_dark_mode))
  })

  output$sim_cmp_table <- renderUI({
    cmp <- compare_result()
    if (is.null(cmp)) return(NULL)

    sig <- cmp$sigma_loc_xy
    a   <- cmp$a;  b <- cmp$b

    pr_a     <- a$protein_r_nm %||% 0
    pr_b     <- b$protein_r_nm %||% 0
    theory_a <- sqrt(pr_a^2 + a$link_nm^2 + sig^2)
    theory_b <- sqrt(pr_b^2 + b$link_nm^2 + sig^2)
    obs_a    <- if (!is.null(a$rmse_vals) && length(a$rmse_vals) >= 1)
                  sprintf("%.1f nm (n=%d)", mean(a$rmse_vals), length(a$rmse_vals))
                else "-"
    obs_b    <- if (!is.null(b$rmse_vals) && length(b$rmse_vals) >= 1)
                  sprintf("%.1f nm (n=%d)", mean(b$rmse_vals), length(b$rmse_vals))
                else "-"

    cell_style <- "padding:4px 8px; white-space:normal; word-break:break-word;"
    hdr  <- function(...) tags$th(style = paste0(cell_style, " text-align:left;"), ...)
    cel  <- function(...) tags$td(style = cell_style, ...)
    bcel <- function(...) tags$td(style = cell_style, tags$b(...))
    row2 <- function(lbl, va, vb)
      tags$tr(cel(lbl), bcel(va), bcel(vb))

    tags$div(
      class = "well well-sm",
      style = "font-size:0.86em; margin-top:10px; overflow-x:auto;",
      tags$b("Strategy comparison"),
      tags$table(
        style = "width:100%; border-collapse:collapse; margin-top:6px; table-layout:fixed;",
        tags$colgroup(
          tags$col(style = "width:36%;"),
          tags$col(style = "width:32%;"),
          tags$col(style = "width:32%;")
        ),
        tags$thead(tags$tr(
          hdr("Metric"),
          hdr(tags$span(style = "color:#dc2626;", "Strategy A")),
          hdr(tags$span(style = "color:#2563eb;", "Strategy B"))
        )),
        tags$tbody(
          row2("Label",
               a$lbl_note,
               b$lbl_note),
          row2("Fluorophores / molecule (mean)",
               sprintf("%.1f", a$n_fluoro / max(1, nrow(cmp$gt))),
               sprintf("%.1f", b$n_fluoro / max(1, nrow(cmp$gt)))),
          row2("Total TIDs",
               a$n_fluoro, b$n_fluoro),
          row2("Linker distance",
               sprintf("%.0f nm", a$link_nm),
               sprintf("%.0f nm", b$link_nm)),
          row2("Label footprint",
               if (!is.na(a$footprint_nm)) sprintf("~%.0f nm", a$footprint_nm) else "-",
               if (!is.na(b$footprint_nm)) sprintf("~%.0f nm", b$footprint_nm) else "-"),
          if (pr_a > 0 || pr_b > 0)
            row2("Protein radius",
                 sprintf("%.1f nm", pr_a),
                 sprintf("%.1f nm", pr_b)),
          row2("Combined σ",
               sprintf("%.1f nm", theory_a),
               sprintf("%.1f nm", theory_b)),
          row2("Observed TID displacement",
               obs_a, obs_b)
        )
      )
    )
  })

  # ── Confocal view: live PSF preview  ────────────────────
  output$sim_diff_psf_preview <- renderUI({
    na      <- max(0.5, min(1.5, as.numeric(input$sim_diff_NA %||% 1.4)))
    au      <- max(0.1, as.numeric(input$sim_diff_AU %||% 0.83))
    sig     <- compute_confocal_sigma(640, 670, na, au)
    fwhm    <- round(sig * 2.355)
    au_nm   <- round(1.22 * 670 / na)
    ee      <- confocal_collection_efficiency(au)
    tags$div(
      class = "alert alert-info",
      style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
      tags$b(sprintf("\u03C3 \u2248 %.1f nm  (FWHM \u2248 %d nm)  \u00B7  signal \u2248 %.0f%%",
                     sig, fwhm, 100 * ee)),
      tags$br(),
      tags$small(sprintf(
        "1 AU = %d nm  \u00B7  pinhole radius = %.0f nm  \u00B7  closing sharpens + sections but dims",
        au_nm, au * au_nm / 2))
    )
  })

  # ── Confocal view plot ───────────────────────────────────────────────────────
  output$sim_diff_plot <- renderPlotly({
    res <- sim_result()
    shiny::validate(
      shiny::need(!is.null(res) && nrow(res$gt) > 0,
                  "Click 'Generate' first to build the GT positions."))

    dark   <- isTRUE(input$sim_dark_mode)
    bg_col <- if (dark) "#111111" else "#FFFFFF"
    fg_col <- if (dark) "#DDDDDD" else "#444444"
    gc     <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"

    na      <- max(0.5, min(1.5, as.numeric(input$sim_diff_NA %||% 1.4)))
    au      <- max(0.1, as.numeric(input$sim_diff_AU %||% 0.83))
    psf_sig <- compute_confocal_sigma(640, 670, na, au)
    n_ph    <- max(1, as.numeric(input$sim_diff_photons %||% 200))
    bg_ph   <- max(0, as.numeric(input$sim_diff_bg      %||% 5))

    # Pinhole physics: the collected signal is scaled by the encircled energy
    # (closing the pinhole → dimmer, noisier via Poisson), while the out-of-focus
    # haze grows with the pinhole AREA (~AU², opening → washed-out background).
    ee       <- confocal_collection_efficiency(au)
    n_ph_eff <- n_ph * ee                 # signal photons actually collected
    bg_eff   <- bg_ph * au^2              # out-of-focus background admitted

    # Render from the UNIQUE fluorophore pool (one physical dye = one PSF), not
    # the over-counted MINFLUX detection list - re-detecting a molecule over the
    # acquisition does not make its confocal spot brighter.
    #   • Each fluorophore = one PSF at its displaced position, n_ph photons.
    #     Multi-dye labels (e.g. polyclonal + secondary) still stack into a
    #     brighter spot because they are distinct dyes.
    #   • Background TIDs (non-specific binding) render as single n_ph spots.
    fp     <- if (!is.null(res$fluoro_pool)) res$fluoro_pool else res$fluoro
    src_x  <- if (nrow(fp) > 0) fp$x else numeric(0)
    src_y  <- if (nrow(fp) > 0) fp$y else numeric(0)
    src_ph <- rep(n_ph_eff, length(src_x))

    if (!is.null(res$bg_pos) && nrow(res$bg_pos) > 0) {
      src_x  <- c(src_x,  res$bg_pos$x)
      src_y  <- c(src_y,  res$bg_pos$y)
      src_ph <- c(src_ph, rep(n_ph_eff, nrow(res$bg_pos)))
    }

    # Cap at 5 000 sources total to keep rendering fast (sample uniformly)
    if (length(src_x) > 5000L) {
      idx    <- sample(length(src_x), 5000L)
      src_x  <- src_x[idx]; src_y <- src_y[idx]; src_ph <- src_ph[idx]
    }

    shiny::validate(shiny::need(length(src_x) > 0,
      "No fluorophores to render. Adjust labeling efficiency or generate first."))

    px_user  <- max(1, as.numeric(input$sim_diff_px %||% 50))
    # Pin the confocal window to the FOV when it's active so the confocal spans
    # the same ROI as the MINFLUX view rather than cropping to the sources.
    use_fov  <- sim_use_fov()
    fov_ext  <- if (use_fov) {
                  fx <- max(1, as.numeric(input$sim_fov_x %||% 4000))
                  fy <- max(1, as.numeric(input$sim_fov_y %||% 4000))
                  c(-fx / 2, fx / 2, -fy / 2, fy / 2)
                } else NULL
    img_data <- render_diffraction_image(src_x, src_y, psf_sig, n_ph_eff, bg_eff,
                                          px_nm_user = px_user,
                                          photons_vec = src_ph,
                                          extent = fov_ext)
    shiny::validate(shiny::need(!is.null(img_data), "Rendering failed."))

    # FWHM = 2.355 × σ: optical resolution - apparent width of a point source
    # and the minimum resolvable separation between two identical emitters.
    fwhm_nm   <- round(psf_sig * 2.355)
    title_txt <- sprintf(
      "Confocal view  (\u03C3 \u2248 %.0f nm, FWHM \u2248 %d nm, signal %.0f%% @ %.2f AU)",
      psf_sig, fwhm_nm, 100 * ee, au)

    # Hot Red LUT (ImageJ "Red Hot") - fixed regardless of dark mode; the
    # confocal image is always rendered on black. Low → high:
    # black → red → orange → yellow → white.
    cscale <- list(
      list(0.00, "rgb(0,0,0)"),
      list(0.25, "rgb(128,0,0)"),
      list(0.45, "rgb(255,0,0)"),
      list(0.65, "rgb(255,140,0)"),
      list(0.85, "rgb(255,255,0)"),
      list(1.00, "rgb(255,255,255)"))

    # The confocal image is inherently a black-background picture, so render it
    # on black with no gridlines and white axis text regardless of dark mode.
    cf_bg <- "#000000"; cf_fg <- "#DDDDDD"

    # Scale bar - same controls as the MINFLUX view, forced white so it stands
    # out on the black image, measured against the rendered image width.  The
    # full image span is stored so the zoom observer can re-scale the bar.
    x_span <- diff(range(img_data$xs, na.rm = TRUE))
    sim_diff_full(x_span)
    sb <- NULL
    if (isTRUE(input$sim_show_scale_bar)) {
      sb <- make_scale_bar(
        suppressWarnings(as.numeric(input$sim_scale_bar_length)),
        x_span,
        suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
        dark = TRUE,
        scale_font_size = suppressWarnings(
          as.numeric(input$sim_scale_font_size %||% 18)))
    }

    lay_args <- list(
      paper_bgcolor = cf_bg, plot_bgcolor = cf_bg,
      font  = list(color = cf_fg),
      xaxis = list(title = "X (nm)", color = cf_fg,
                   scaleanchor = "y", scaleratio = 1,
                   showgrid = FALSE, zeroline = FALSE),
      yaxis = list(title = "Y (nm)", color = cf_fg,
                   showgrid = FALSE, zeroline = FALSE),
      title = list(text = title_txt, font = list(color = cf_fg)))
    if (!is.null(sb)) {
      lay_args$shapes      <- sb$shapes
      lay_args$annotations <- sb$annotations
    }

    fig <- plot_ly(source = "sim_diff") %>%
      add_trace(
        type       = "heatmap",
        x          = img_data$xs, y = img_data$ys, z = img_data$img,
        colorscale = cscale,
        # Auto-contrast the display (like ImageJ auto B&C) so it looks like a real
        # confocal image. The pinhole still shows through the DATA - closing it
        # collects fewer photons (grainier via Poisson) and tightens the PSF;
        # opening it adds out-of-focus background - not through raw brightness.
        showscale  = FALSE, hoverinfo = "skip",
        name       = "Diffraction-limited"
      )
    sim_diff_ready(TRUE)
    do.call(layout, c(list(fig), lay_args)) %>%
      plotly::event_register("plotly_relayout") %>%
      config(displaylogo = FALSE,
             toImageButtonOptions = list(format = "svg",
                                         filename = "diffraction_limited"))
  })

  # ── Phase 1 CSV export ───────────────────────────────────────────────────────
  # Each row is one localisation.  When GT data exist, mol_id and the GT
  # coordinates (gt_x, gt_y, gt_z) are appended so every loc can be traced
  # back to its source molecule.  Background TIDs get NA in these columns.
  output$sim_dl_csv <- downloadHandler(
    filename = function()
      paste0("sim_localisations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- sim_result()
      shiny::validate(shiny::need(!is.null(res) && nrow(res$locs) > 0,
                                  "No localizations to export."))
      out <- res$locs[, intersect(c("x", "y", "z", "tid", "tim"),
                                  names(res$locs)), drop = FALSE]

      if ("mol_id" %in% names(res$fluoro) && nrow(res$gt) > 0) {
        n_real  <- nrow(res$fluoro)
        tid_idx <- suppressWarnings(as.integer(out$tid))
        valid   <- !is.na(tid_idx) & tid_idx >= 1L & tid_idx <= n_real
        mol_ids <- rep(NA_integer_, nrow(out))
        mol_ids[valid] <- res$fluoro$mol_id[tid_idx[valid]]
        out$mol_id <- mol_ids
        out$gt_x   <- ifelse(is.na(mol_ids), NA_real_, res$gt$x[mol_ids])
        out$gt_y   <- ifelse(is.na(mol_ids), NA_real_, res$gt$y[mol_ids])
        if ("z" %in% names(res$gt))
          out$gt_z <- ifelse(is.na(mol_ids), NA_real_, res$gt$z[mol_ids])
      }

      write.csv(out, file, row.names = FALSE)
    }
  )

  # ── Acquisition timeline: resolution-over-time plot ──────────────────────────
  output$sim_acq_plot <- renderPlotly({
    shiny::validate(shiny::need(
      identical(input$sim_mode %||% "normal", "acquisition"),
      "Switch to Acquisition mode (top of the sidebar) to model an acquisition timeline."))
    a <- sim_acq_live()
    shiny::validate(shiny::need(
      !is.null(a),
      "Click 'Generate' to build the acquisition timeline."))
    dark   <- isTRUE(input$sim_dark_mode)
    bg_col <- if (dark) "#111111" else "#FFFFFF"
    fg_col <- if (dark) "#DDDDDD" else "#444444"
    gc     <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"

    # Work in minutes; lead with the detection rate (traces/hour) since that is
    # the curve directly comparable to a real run (early peak → constant floor).
    rt_min  <- a$rate_t / 60
    t_min   <- a$t / 60
    # % resolved is relative to the LABELED molecules, so the curve plateaus at
    # the fraction actually reachable in this run (reaches 100 % when all
    # labeled molecules are found).
    cum_pct <- 100 * a$cum_unique / max(1, a$n_mol_resolvable)
    rmax    <- max(a$rate_v, 1, na.rm = TRUE)

    p <- plot_ly()
    if (!is.null(a$rate_t))
      p <- p %>% add_trace(x = rt_min, y = a$rate_v, type = "scatter",
                           mode = "lines", fill = "tozeroy",
                           fillcolor = "rgba(224,123,57,0.22)",
                           line = list(color = "#e07b39", width = 1.5),
                           name = "Detection rate", hoverinfo = "x+y")
    # Optional decomposition of the PREDICTED rate into its structure (decaying +
    # bleaching) and constant background components - makes explicit that the
    # orange total = structure + background, and the floor is background.
    # Convention: solid = model, dotted = measured (real run, below).
    if (isTRUE(input$sim_real_show_split) && !is.null(a$rate_struct)) {
      p <- p %>%
        add_trace(x = rt_min, y = a$rate_struct, type = "scatter",
                  mode = "lines",
                  line = list(color = "#2a9d8f", width = 1.5),
                  name = "Structure (model)", hoverinfo = "x+y") %>%
        add_segments(x = 0, xend = max(rt_min),
                     y = a$bg_rate, yend = a$bg_rate,
                     line = list(color = "#b07aa1", width = 1.5),
                     name = "Background (model)", inherit = FALSE)
    }
    # Late-phase floor line (median of the last quarter of the run).
    if (!is.null(a$rate_floor) && is.finite(a$rate_floor))
      p <- p %>% add_segments(x = 0, xend = max(rt_min),
                              y = a$rate_floor, yend = a$rate_floor,
                              line = list(color = "#e07b39", width = 1,
                                          dash = "dot"),
                              name = sprintf("floor ~%.0f/hr", a$rate_floor),
                              inherit = FALSE)
    # Cumulative molecules resolved (%) on the secondary axis.
    p <- p %>% add_trace(x = t_min, y = cum_pct, type = "scatter",
                         mode = "lines", yaxis = "y2",
                         line = list(color = "#3b82c4", width = 2),
                         name = "Molecules resolved", hoverinfo = "x+y")
    # Measured detection rate from a loaded real run (calibration overlay).
    # Optionally decomposed into the spatially-split structure / background.
    real <- sim_real_rv()
    if (!is.null(real)) {
      if (isTRUE(input$sim_real_show_split) && !is.null(real$rate_v_struct)) {
        p <- p %>%
          add_trace(x = real$rate_t / 60, y = real$rate_v_struct,
                    type = "scatter", mode = "lines",
                    line = list(color = "#2a9d8f", width = 1.5, dash = "dot"),
                    name = "Measured structure", hoverinfo = "x+y") %>%
          add_trace(x = real$rate_t / 60, y = real$rate_v_bg,
                    type = "scatter", mode = "lines",
                    line = list(color = "#b07aa1", width = 1.5, dash = "dot"),
                    name = "Measured background", hoverinfo = "x+y")
      } else {
        p <- p %>% add_trace(x = real$rate_t / 60, y = real$rate_v,
                             type = "scatter", mode = "lines",
                             line = list(color = "#777777", width = 1.5,
                                         dash = "dot"),
                             name = "Measured (real run)", hoverinfo = "x+y")
      }
    }
    # Resolution marker only when the threshold is actually reached.
    if (!is.null(a$res_time) && is.finite(a$res_time))
      p <- p %>% add_segments(x = a$res_time / 60, xend = a$res_time / 60,
                              y = 0, yend = rmax,
                              line = list(color = "#16a34a", width = 1.5,
                                          dash = "dash"),
                              name = sprintf("%.0f%% resolved", 100 * a$resolve_q),
                              inherit = FALSE)
    p %>%
      layout(
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font   = list(color = fg_col),
        margin = list(l = 60, r = 95, t = 50, b = 110),
        xaxis  = list(title = "Acquisition time (min)", color = fg_col,
                      gridcolor = gc),
        yaxis  = list(title = "Detection rate (traces/hour)", color = "#e07b39",
                      gridcolor = gc, rangemode = "tozero", automargin = TRUE),
        yaxis2 = list(title = list(text = "Molecules resolved (%)",
                                   standoff = 22),
                      color = "#3b82c4", overlaying = "y", side = "right",
                      range = c(0, 105), showgrid = FALSE, automargin = TRUE),
        title  = list(text = "Detection rate & resolution over time",
                      font = list(color = fg_col), x = 0.02, xanchor = "left",
                      y = 0.97, yanchor = "top"),
        showlegend = TRUE,
        legend = list(orientation = "h", x = 1, xanchor = "right",
                      y = -0.28, yanchor = "top", font = list(size = 10))
      ) %>%
      config(displaylogo = FALSE)
  })

  # ── Acquisition timeline summary ─────────────────────────────────────────────
  output$sim_acq_summary <- renderUI({
    a <- sim_acq_live()
    if (is.null(a)) return(NULL)

    fmt_t <- function(s) if (s >= 60) sprintf("%.1f min", s / 60)
                         else          sprintf("%.1f s", s)
    reached <- !is.null(a$res_time) && is.finite(a$res_time)
    capped  <- !is.null(a$struct_cap) && a$struct_cap > 0 &&
               a$n_real_tid >= a$struct_cap

    row <- function(lbl, val)
      tags$tr(tags$td(style = "padding:2px 8px;", lbl),
              tags$td(style = "padding:2px 8px;", tags$b(val)))

    tags$div(
      class = "well well-sm", style = "font-size:0.88em; margin-top:10px;",
      tags$b("Acquisition timeline"), tags$br(),
      tags$table(style = "width:100%; border-collapse:collapse;",
        row("Total acquisition time", fmt_t(a$T_sec)),
        row("Structure detection rate",
            sprintf("%.0f traces/hour  (mean %.0f)", a$struct_rate, a$rate_mean)),
        row("Structure trace IDs",
            if (capped) sprintf("%d  (photon budget reached)", a$n_real_tid)
            else        sprintf("%d", a$n_real_tid)),
        row("Molecules resolved",
            sprintf("%d / %d labeled  (%.0f%%; %d of %d total)",
                    a$n_mol_det, a$n_mol_resolvable, 100 * a$resolved_frac,
                    a$n_mol_det, a$n_mol_total)),
        row(sprintf("Time to %.0f%% resolved", 100 * a$resolve_q),
            if (reached) fmt_t(a$res_time)
            else sprintf("not reached (%.0f%% max)", 100 * a$resolved_frac)),
        if (reached)
          row("Time after resolved",
              sprintf("%s  (%.0f%% of run, mostly re-detecting)",
                      fmt_t(max(0, a$T_sec - a$res_time)),
                      100 * max(0, a$T_sec - a$res_time) / max(1e-9, a$T_sec))),
        row("Detections per molecule",
            sprintf("%.1f× on average  (%.0f%% of traces are repeats)",
                    a$det_per_mol, 100 * a$overcount_frac)),
        if (!is.null(a$n_bg) && a$n_bg > 0)
          row("Background TIDs",
              sprintf("%d  (@ %.0f/hr)", a$n_bg, a$bg_rate))
      ),
      tags$small(style = "color:#888;",
        if (reached)
          paste("Resolved early, then later traces re-detect the same",
                "molecules - the over-counting you group as repeat blinks.",
                "Increase the rate or duration to resolve more; decrease them",
                "(or shorten the run) to see partial resolution.")
        else
          paste(sprintf("Only %.0f%% of the structure resolved.",
                        100 * a$resolved_frac),
                "Increase the detection rate or duration (or improve labeling)",
                "to reach the resolve threshold."))
    )
  })

  # ── Phase 1 GT positions export ──────────────────────────────────────────────
  output$sim_dl_gt_csv <- downloadHandler(
    filename = function()
      paste0("sim_gt_positions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- sim_result()
      shiny::validate(shiny::need(!is.null(res) && nrow(res$gt) > 0,
                                  "No ground-truth positions to export."))
      cols <- intersect(c("x","y","z"), names(res$gt))
      gt_out <- res$gt[, cols, drop = FALSE]
      gt_out$mol_id <- seq_len(nrow(gt_out))
      write.csv(gt_out, file, row.names = FALSE)
    }
  )

  # ── .sim_fgn() - fractional Gaussian noise (Davies-Harte) ───────────────────
  # GOAL: increments for fractional Brownian motion, the standard model of
  # ANOMALOUS diffusion. Successive increments are correlated (anti-correlated for
  # H < 0.5), which makes MSD grow as t^alpha with alpha = 2H instead of t^1.
  # A cell-wall protein in a crowded peptidoglycan mesh is subdiffusive: the GP67H
  # files measure alpha ~ 0.78.
  # n: number of increments. H: Hurst exponent (alpha/2), 0 < H < 1.
  # returns: numeric(n), unit variance per increment, correlated.
  # Exact via circulant embedding; falls back to white noise at H = 0.5 (Brownian).
  .sim_fgn <- function(n, H) {
    if (n < 2L) return(stats::rnorm(n))
    if (abs(H - 0.5) < 1e-6) return(stats::rnorm(n))
    k <- 0:(n - 1)
    g <- 0.5 * (abs(k + 1)^(2*H) - 2 * abs(k)^(2*H) + abs(k - 1)^(2*H))  # fGn autocovariance
    cc  <- c(g, g[n:2])                       # circulant first row, length m = 2n-2
    m   <- length(cc)
    lam <- Re(stats::fft(cc))
    if (any(lam < 0)) lam <- pmax(lam, 0)     # embedding not non-negative: clip (rare)
    half <- m / 2
    w <- complex(length.out = m)
    w[1]        <- sqrt(lam[1] / m) * stats::rnorm(1)
    w[half + 1] <- sqrt(lam[half + 1] / m) * stats::rnorm(1)
    j  <- 2:half
    zr <- stats::rnorm(length(j)); zi <- stats::rnorm(length(j))
    w[j]         <- sqrt(lam[j] / (2 * m)) * complex(real = zr, imaginary = zi)
    w[m - j + 2] <- Conj(w[j])
    Re(stats::fft(w))[1:n]
  }

  # ══════════════════════════════════════════════════════════════════════════
  # Diffusion / Tracking ground-truth simulator (Phase 1)
  # See docs/simulation_module.md. Free surface Brownian on a sphere or rod,
  # sampled the MINFLUX way (cycle clock + dark cycles + heteroscedastic σ_loc +
  # photobleaching) → MINFLUX-format tracks + a ground-truth D. The CVE round-trip
  # in the summary demonstrates the validation loop in-place.
  # ══════════════════════════════════════════════════════════════════════════
  # --- Calibrate the sampling model from a real MINFLUX file -----------------
  # Reads a real .csv/.npy/.json and auto-fills cycle time, photons, dark-cycle
  # probability and mean track length so the simulator matches a real experiment.
  sim_trk_calib_rv <- reactiveVal(NULL)
  observeEvent(input$sim_trk_calib, {
    fi <- input$sim_trk_calib_file
    shiny::validate(shiny::need(!is.null(fi), "Choose a MINFLUX file first."))
    p <- fi$datapath; nm <- tolower(fi$name)
    df <- withProgress(message = "Reading MINFLUX file…", value = 0, {
      tryCatch({
        if (grepl("\\.npy$",  nm)) read_minflux_npy_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
        else if (grepl("\\.json$", nm)) read_minflux_json_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
        else if (grepl("\\.csv$",  nm)) {
          x <- as.data.frame(data.table::fread(p, showProgress = FALSE))
          if (all(c("tid","tim") %in% names(x))) x else canonicalize_minflux_df(x)
        } else NULL
      }, error = function(e) NULL)
    })
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0 && all(c("tid","tim") %in% names(df)),
      "Could not read tid / tim from that file."))
    df <- df[is.finite(df$tim), ]
    # Calibrate ONLY on analysable traces (≥ MIN_LEN localizations, no upper bound), so
    # the short-junk traces don't drag the statistics down. The representative length is
    # the MEDIAN of the kept traces.
    MIN_LEN  <- max(2L, as.integer(input$sim_trk_calib_minlen %||% 100L))
    len_all  <- table(df$tid)
    keep_tid <- names(len_all)[len_all >= MIN_LEN]
    shiny::validate(shiny::need(length(keep_tid) >= 5, sprintf("Fewer than 5 traces have ≥%d localizations.", MIN_LEN)))
    n_total  <- length(len_all)
    df <- df[as.character(df$tid) %in% keep_tid, ]        # keep only the long traces
    len     <- as.numeric(table(df$tid)); med_len <- stats::median(len)
    photons <- if ("eco" %in% names(df)) stats::median(df$eco[is.finite(df$eco) & df$eco > 0]) else NA_real_
    o    <- order(df$tid, df$tim); df <- df[o, ]
    tid <- df$tid; tim <- df$tim; nrw <- length(tid)
    same <- c(FALSE, tid[-1] == tid[-nrw]); dt <- c(NA_real_, diff(tim)) * 1000; dt[!same] <- NA_real_  # ms
    dtv  <- dt[is.finite(dt) & dt > 0]
    shiny::validate(shiny::need(length(dtv) > 20, "Too few within-track steps to calibrate."))
    # ── Cycle time and SW from the Δt quantisation ────────────────────────────
    # Real Δt obey Δt = k·t_c + SW exactly, so the MODAL Δt is t_c + SW, not t_c.
    # Take the mode as a first guess, assign each Δt its cycle count k, then read
    # t_c off the SPACING between the levels and SW off the intercept. Taking the
    # mode for the cycle time inflates it by SW (~5% on the GP67H files).
    hb    <- hist(dtv[dtv < stats::quantile(dtv, 0.98)], breaks = 400, plot = FALSE)
    modal <- hb$mids[which.max(hb$counts)]                         # modal Δt = t_c + SW
    k0    <- round(dtv / modal)
    ok    <- k0 >= 1 & k0 <= 10 & abs(dtv - k0 * modal) < 0.3 * modal   # on-level Δt only
    cyc <- modal; sw <- 0
    if (sum(ok) > 20) {
      cf <- stats::coef(stats::lm(dtv[ok] ~ k0[ok]))
      if (is.finite(cf[2]) && cf[2] > 0 && cf[1] >= 0 && cf[1] < cf[2]) {
        cyc <- unname(cf[2]); sw <- unname(cf[1])
      }
    }
    ncyc <- pmax(1L, as.integer(round((dtv - sw) / cyc)))
    p_dark <- sum(ncyc - 1L) / sum(ncyc); integ <- mean(ncyc == 1L)
    # ── D and σ_loc from the single-cycle steps (CVE, blur-aware) ─────────────
    # σ² = R⟨Δx²⟩ + (2R−1)·cov₁ per axis with R = 1/6 (full-cycle integration);
    # D itself is blur-independent. Lateral axes only: real 2D tracking files
    # carry a placeholder z.
    dstep <- same & is.finite(dt) & abs(round((dt - sw) / cyc) - 1) < 1e-6
    dstep[is.na(dstep)] <- FALSE
    dx <- c(NA_real_, diff(df$x)); dy <- c(NA_real_, diff(df$y)); dx[!same] <- NA_real_; dy[!same] <- NA_real_
    pair <- dstep & c(dstep[-1], FALSE)
    cov1 <- mean(c(dx[pair] * c(dx[-1], NA_real_)[pair],
                   dy[pair] * c(dy[-1], NA_real_)[pair]), na.rm = TRUE)
    msq  <- mean(c(dx[dstep]^2, dy[dstep]^2), na.rm = TRUE)
    D_meas <- if (is.finite(msq) && is.finite(cov1)) (msq / (2 * (cyc + sw)) + cov1 / (cyc + sw)) / 1000 else NA_real_
    sig_meas <- if (is.finite(msq) && is.finite(cov1)) sqrt(max(0, msq / 6 - (2/3) * cov1)) else NA_real_
    # -- Photon budget: collection window, photon limit, repeats ---------------
    # tau = eco/efo is the collection window, written in the file per localization and
    # quantised in units of one pattern repeat. It gives the blur coefficient and the
    # repeat structure without any assumption.
    pht_lim <- NA_real_; expf <- NA_real_; nrep <- NA_real_
    if (all(c("eco","efo") %in% names(df))) {
      good <- is.finite(df$eco) & is.finite(df$efo) & df$efo > 0 & df$eco > 0
      if (sum(good) > 100) {
        tau_v   <- df$eco[good] / df$efo[good] * 1000                  # ms
        tau0    <- unname(stats::quantile(tau_v, 0.05))                # one pattern repeat
        pht_lim <- unname(stats::quantile(df$eco[good], 0.001))        # hard floor = phtLimit
        expf    <- tau0 / cyc
        nrep    <- unname(stats::quantile(tau_v / tau0, 0.999))
      }
    }
    # -- Brightness: between molecules, within a track, correlation, fade ------
    b_btw <- NA_real_; b_wth <- NA_real_; b_phi <- NA_real_; b_dec <- NA_real_
    if ("efo" %in% names(df)) {
      dd  <- df[is.finite(df$efo) & df$efo > 0, ]
      spl <- split(log(dd$efo), dd$tid)
      spl <- spl[vapply(spl, length, 1L) >= 20L]
      if (length(spl) >= 5L) {
        b_btw <- stats::sd(vapply(spl, stats::median, 0))              # molecule to molecule
        b_wth <- stats::median(vapply(spl, stats::sd, 0))              # inside one track
        b_phi <- stats::median(vapply(spl, function(v)
                   suppressWarnings(stats::cor(v[-length(v)], v[-1])), 0), na.rm = TRUE)
        # Fade: log drop from the first tenth of a track to the last, POOLED over
        # localizations rather than averaged over tracks. The per-track median measures
        # something else - short tracks mostly end by being lost rather than by
        # bleaching, so they barely fade - and calibrating on it under-states the fade
        # by a factor of about seven on this data.
        dec <- unlist(lapply(spl, function(v) {
                 f <- (seq_along(v) - 1) / max(1, length(v) - 1); pmin(9L, floor(f * 10)) }))
        vv  <- unlist(spl, use.names = FALSE)
        b_dec <- stats::median(vv[dec == 0L]) - stats::median(vv[dec == 9L])
      }
    }
    # -- Track-length heterogeneity --------------------------------------------
    # The generator compounds a geometric lifetime with a lognormal mean; mean/median of
    # the resulting lengths is a monotone function of that spread, inverted here.
    ratio  <- mean(len) / stats::median(len)
    lifecv <- stats::approx(c(1.45, 1.68, 2.14, 2.82, 3.60, 5.01), c(0, 0.5, 1, 1.5, 2, 3),
                            xout = ratio, rule = 2)$y
    # -- D-lifetime coupling ----------------------------------------------------
    # Per-track D against track length. Both are noisy, so the correlation the simulator
    # reproduces is weaker than the one measured; the value filled in is de-attenuated.
    rho_meas <- NA_real_; rho_set <- NA_real_
    if (all(c("x","y") %in% names(df))) {
      pt <- data.frame(tid = df$tid, dx = dx, dy = dy,
                       st = dstep, pr = dstep & c(dstep[-1], FALSE))
      sp <- split(pt, pt$tid)
      dd <- vapply(sp, function(z) {
              ii <- which(z$st); jj <- which(z$pr)
              if (length(ii) < 30 || length(jj) < 10) return(NA_real_)
              (mean(c(z$dx[ii]^2, z$dy[ii]^2), na.rm = TRUE) / (2*(cyc+sw)) +
               mean(c(z$dx[jj]*z$dx[jj+1], z$dy[jj]*z$dy[jj+1]), na.rm = TRUE) / (cyc+sw)) / 1000
            }, 0)
      nn <- vapply(sp, nrow, 1L)
      gd <- is.finite(dd) & dd > 0
      if (sum(gd) >= 10) {
        rho_meas <- suppressWarnings(stats::cor(log(nn[gd]), log(dd[gd]), method = "spearman"))
        rho_set  <- max(-0.95, min(0.95, 2 * rho_meas))
      }
    }
    setn <- function(id, v, dig = 2) if (is.finite(v)) updateNumericInput(session, id, value = round(v, dig))
    if (is.finite(cyc))     updateNumericInput(session, "sim_trk_cycle",   value = signif(cyc, 3))
    if (is.finite(sw))      updateNumericInput(session, "sim_trk_sw",      value = signif(sw, 2))
    if (is.finite(p_dark))  updateNumericInput(session, "sim_trk_pdark",   value = round(p_dark, 2))
    if (is.finite(photons)) updateNumericInput(session, "sim_trk_photons", value = round(photons))
    if (is.finite(med_len)) updateNumericInput(session, "sim_trk_meanlen", value = round(med_len))
    if (is.finite(D_meas) && D_meas > 0)
      updateNumericInput(session, "sim_trk_D", value = signif(D_meas, 2))
    if (length(keep_tid) >= 1) updateNumericInput(session, "sim_trk_ntracks", value = length(keep_tid))
    setn("sim_trk_phtlimit", pht_lim, 0); setn("sim_trk_expfrac", expf)
    setn("sim_trk_taurep",   nrep, 0)
    setn("sim_trk_phtbtw",   b_btw); setn("sim_trk_phtcv", b_wth)
    setn("sim_trk_phtphi",   max(0, b_phi)); setn("sim_trk_phtdecay", max(0, b_dec))
    setn("sim_trk_lifecv",   lifecv); setn("sim_trk_Dlife", rho_set)
    sim_trk_calib_rv(list(file = fi$name, n_tracks = length(keep_tid), n_total = n_total,
      min_len = MIN_LEN, n_loc = nrow(df), med_len = med_len, mean_len = mean(len),
      photons = photons, cyc = cyc, sw = sw, modal = modal, p_dark = p_dark, integ = integ,
      D = D_meas, sigma = sig_meas, med_dt = stats::median(dtv),
      pht_lim = pht_lim, expf = expf, nrep = nrep, b_btw = b_btw, b_wth = b_wth,
      b_phi = b_phi, b_dec = b_dec, lifecv = lifecv, rho_meas = rho_meas, rho_set = rho_set))
    showNotification(sprintf("Calibrated on %d traces (≥%d locs) of %d total.",
                             length(keep_tid), MIN_LEN, n_total), type = "message", duration = 6)
  })

  output$sim_trk_calib_report <- renderUI({
    cb <- sim_trk_calib_rv(); if (is.null(cb)) return(NULL)
    rowf <- function(k, v) htmltools::tags$tr(
      htmltools::tags$td(style = "padding:1px 10px 1px 0; color:#666;", k),
      htmltools::tags$td(htmltools::tags$b(v)))
    htmltools::div(style = "font-size:0.8em; color:#333; margin-top:4px;",
      htmltools::tags$div(style = "color:#228833;", sprintf("✓ Calibrated from %s", cb$file)),
      htmltools::tags$table(style = "border-collapse:collapse;",
        rowf(sprintf("Traces used (≥%d locs)", cb$min_len %||% 50L),
             sprintf("%s of %s  ·  %s locs", format(cb$n_tracks, big.mark = ","),
                     format(cb$n_total %||% cb$n_tracks, big.mark = ","), format(cb$n_loc, big.mark = ","))),
        rowf("Cycle time", sprintf("%.4g ms  (SW %.3g ms · modal Δt %.4g)", cb$cyc, cb$sw %||% 0, cb$modal %||% NA)),
        rowf("Photons/loc (median eco)", if (is.finite(cb$photons)) sprintf("%.0f", cb$photons) else "-"),
        rowf("Dark-cycle prob.", sprintf("%.2f  (%.0f%% integrated)", cb$p_dark, 100 * cb$integ)),
        rowf("Track length (median)", sprintf("%.0f · mean %.0f", cb$med_len, cb$mean_len)),
        rowf("D (CVE, lateral)", if (is.finite(cb$D %||% NA)) sprintf("%.3g µm²/s", cb$D) else "-"),
        rowf("σ_loc (CVE, measured R)", if (is.finite(cb$sigma %||% NA)) sprintf("%.2f nm", cb$sigma) else "-"),
        rowf("Collection window", if (is.finite(cb$expf %||% NA))
               sprintf("%.0f%% of the cycle · up to %.0f repeats · photon limit %.0f",
                       100*cb$expf, cb$nrep, cb$pht_lim) else "-"),
        rowf("Brightness", if (is.finite(cb$b_btw %||% NA))
               sprintf("spread %.2f between / %.2f within · correlation %.2f · fades to %.0f%%",
                       cb$b_btw, cb$b_wth, cb$b_phi, 100*exp(-cb$b_dec)) else "-"),
        rowf("Length spread", sprintf("CV %.1f  (mean/median = %.2f)", cb$lifecv, cb$mean_len/cb$med_len)),
        rowf("D vs track length", if (is.finite(cb$rho_meas %||% NA))
               sprintf("%+.2f measured → %+.2f set (noise attenuates it)", cb$rho_meas, cb$rho_set) else "-")),
      htmltools::tags$p(style = "color:#888; margin-top:3px;",
        "Every acquisition field is filled from this file. Set the sample geometry, the motion model ",
        "and the TCP diameter, then Simulate."))
  })

  sim_trk_result <- eventReactive(input$sim_trk_run, {
    set.seed(as.integer(input$sim_trk_seed %||% 1L))
    geom     <- input$sim_trk_geom %||% "rod"
    mode     <- input$sim_trk_mode %||% "surface"      # surface | volume (intracellular)
    mthick   <- max(0, as.numeric(input$sim_trk_mthick %||% 6))              # membrane thickness (nm)
    # Membrane tilt: the sheet is primarily in XY but inclined by mtilt°, so the plane
    # is spanned by e1,e2 (normal n̂ near ẑ). The path is generated flat (local u,v,w)
    # then rotated by R_mem = [e1 e2 n̂] → z varies across the plane. 0° = flat.
    mtilt <- max(0, min(89, as.numeric(input$sim_trk_mtilt %||% 15))) * pi / 180
    maz   <- runif(1, 0, 2 * pi)
    n_hat <- c(sin(mtilt) * cos(maz), sin(mtilt) * sin(maz), cos(mtilt))
    e1    <- { cx <- c(-n_hat[2], n_hat[1], 0); m <- sqrt(sum(cx^2)); if (m < 1e-9) c(1,0,0) else cx / m }  # ẑ×n̂ (in XY)
    e2    <- c(n_hat[2]*e1[3] - n_hat[3]*e1[2], n_hat[3]*e1[1] - n_hat[1]*e1[3], n_hat[1]*e1[2] - n_hat[2]*e1[1])  # n̂×e1
    R_mem <- cbind(e1, e2, n_hat)                                            # local (u,v,w) → world
    # motion DOF (used by the CVE check): free=3, membrane=2 (lateral), else surface=2/volume=3
    d_mot    <- if (identical(geom, "free")) 3L else if (identical(geom, "membrane")) 2L
                else if (identical(mode, "volume")) 3L else 2L
    D_um2_s  <- max(0, as.numeric(input$sim_trk_D %||% 0.05))
    D_nm2ms  <- D_um2_s * 1000                        # µm²/s → nm²/ms
    n_trk    <- max(1L, as.integer(input$sim_trk_ntracks %||% 200L))
    mean_len <- max(5L, as.integer(input$sim_trk_meanlen %||% 80L))
    cyc      <- max(0.01, as.numeric(input$sim_trk_cycle %||% 0.192))   # ms, spacing between cycles
    # Fixed per-localization overhead (settling / window). Real files obey
    # Delta-t = k*cycle + SW exactly, so the modal Delta-t is cycle + SW, not cycle.
    sw       <- max(0, as.numeric(input$sim_trk_sw %||% 0))              # ms
    p_dark   <- min(0.95, max(0, as.numeric(input$sim_trk_pdark %||% 0.28)))
    photons  <- max(1, as.numeric(input$sim_trk_photons %||% 38))
    L_tcp    <- max(1, as.numeric(input$sim_trk_Ltcp %||% 100.8))
    blur     <- isTRUE(input$sim_trk_blur)                               # integrate over the cycle
    # Photon collection window as a fraction of the cycle. tau together with the step
    # duration sets the blur coefficient R = tau/(6*dt): a full-cycle window gives 1/6,
    # a short burst tends to 0.
    expfrac  <- min(1, max(0.05, as.numeric(input$sim_trk_expfrac %||% 0.55)))
    # A localization followed by a dark gap integrates ACROSS that gap (repeated pattern
    # passes) instead of skipping it, up to this many repeats. Measured on real files as
    # tau = eco/efo: 1 repeat before a single-cycle step, 2 before a two-cycle step,
    # saturating near 3. 1 disables it.
    tau_rep  <- max(1L, as.integer(input$sim_trk_taurep %||% 1L))
    # Lognormal spread of D between molecules (0 = every track identical). The input D
    # stays the MEDIAN of the population.
    D_cv     <- max(0, as.numeric(input$sim_trk_Dcv %||% 0))
    D_sdlog  <- if (D_cv > 0) sqrt(log(1 + D_cv^2)) else 0
    # Re-acquisition shrinkage: the first localization after a dark gap is pulled this
    # fraction of the way back toward the previous reported position.
    shrink   <- min(1, max(0, as.numeric(input$sim_trk_shrink %||% 0)))
    # Photon model: the cycle ends once phtLimit photons are in, so eco has a hard
    # floor there; brightness differences between molecules give the long right tail.
    # Rate is lognormal (median = the photons input) and the count Poisson around it,
    # so sdlog = 0 reduces to plain Poisson.
    pht_min  <- max(1L, as.integer(input$sim_trk_phtlimit %||% 1L))
    pht_sd   <- max(0, as.numeric(input$sim_trk_phtcv %||% 0))       # within-track sdlog
    pht_btw  <- max(0, as.numeric(input$sim_trk_phtbtw %||% 0))      # between-molecule sdlog
    pht_phi  <- min(0.95, max(0, as.numeric(input$sim_trk_phtphi %||% 0)))  # AR(1) on log rate
    # Emission fades over a molecule's lifetime: real tracks end at ~44% of their
    # starting rate, and their gap fraction climbs from 18% to 43% as they fade.
    pht_dec  <- max(0, as.numeric(input$sim_trk_phtdecay %||% 0))
    # Photon budget: the tracker repeats the pattern until phtLimit photons are in, and
    # each extra repeat costs a cycle. Dark cycles then FOLLOW from brightness instead of
    # being drawn from p_dark, and eco gets its hard floor at the limit for free.
    budget   <- isTRUE(input$sim_trk_budget)
    # Spread of the bleaching lifetime between molecules. 0 = one geometric lifetime
    # for the whole population; above 0 the per-molecule mean lifetime is lognormal,
    # which makes the recorded-length distribution heavy-tailed like a real file.
    life_cv  <- max(0, as.numeric(input$sim_trk_lifecv %||% 0))
    life_sd  <- if (life_cv > 0) sqrt(log(1 + life_cv^2)) else 0
    # Correlation between a molecule's D and its lifetime. Real files anticorrelate at
    # about −0.4: the slow molecules are the ones tracked longest, which is exactly why
    # a step-weighted pooled D lands below the per-track median.
    D_life   <- min(0.95, max(-0.95, as.numeric(input$sim_trk_Dlife %||% 0)))
    proj2d   <- isTRUE(input$sim_trk_project2d)
    sd_step  <- sqrt(2 * D_nm2ms * cyc)               # per-axis step SD per cycle (nm)
    R    <- max(1, as.numeric(input$sim_trk_R %||% 500))
    rr   <- max(1, as.numeric(input$sim_trk_r %||% 400))
    Lcyl <- max(1, as.numeric(input$sim_trk_Lcyl %||% 2000))
    circ <- 2 * pi * rr
    ea   <- max(1, as.numeric(input$sim_trk_ea %||% 500))
    eb   <- max(1, as.numeric(input$sim_trk_eb %||% 500))
    ec   <- max(1, as.numeric(input$sim_trk_ec %||% 900))
    abc  <- if (identical(geom, "sphere")) c(R, R, R) else c(ea, eb, ec)   # ellipsoid semi-axes
    # Long-axis orientation for the rod: E. coli on a coverslip lie flat → axis in the
    # image plane. Generated in canonical coords (axis ∥ X) then permuted to the choice.
    rod_axis <- input$sim_trk_rod_axis %||% "x"
    rod_perm <- switch(rod_axis, x = c(1L,2L,3L), y = c(2L,1L,3L), z = c(2L,3L,1L))
    # ── FIELD OF VIEW ──────────────────────────────────────────────────────────
    # More than one cell: scatter them over the field at random positions and random
    # in-plane orientations, and give each track a cell. Placement is rejection-based
    # so cells do not sit on top of each other. A single cell stays at the origin, so
    # the estimator-validation view is unchanged.
    n_cells  <- max(1L, as.integer(input$sim_trk_ncells %||% 1L))
    fov_nm   <- max(1, as.numeric(input$sim_trk_fov %||% 10)) * 1000
    cell_len <- if (identical(geom, "rod")) Lcyl + 2 * rr else 2 * max(abc)
    cell_pos <- matrix(0, n_cells, 2); cell_ang <- numeric(n_cells)
    if (n_cells > 1L) {
      for (ci in seq_len(n_cells)) {
        for (try in 1:200) {
          cand <- c(runif(1, -fov_nm/2, fov_nm/2), runif(1, -fov_nm/2, fov_nm/2))
          if (ci == 1L) break
          dmin <- min(sqrt(rowSums((cell_pos[seq_len(ci-1), , drop = FALSE] -
                                    rep(cand, each = ci-1))^2)))
          if (dmin > cell_len) break
        }
        cell_pos[ci, ] <- cand; cell_ang[ci] <- runif(1, 0, 2*pi)
      }
    }

    # ── unified geometry operators (rod = spherocylinder body + hemispherical caps;
    #    else ellipsoid, sphere = a=b=c). project = nearest surface point; reflect =
    #    back inside the volume; normal = unit outward normal (for the tangent step). ─
    is_rod <- identical(geom, "rod")
    geo_normal <- function(p) {
      if (is_rod) { a <- p[1]
        n <- if (a >= 0 && a <= Lcyl) c(0, p[2], p[3]) else { cx <- if (a < 0) 0 else Lcyl; p - c(cx,0,0) }
        m <- sqrt(sum(n^2)); if (m == 0) c(0,1,0) else n / m
      } else { n <- p / abc^2; n / sqrt(sum(n^2)) }
    }
    geo_project <- function(p) {
      if (is_rod) { a <- p[1]
        if (a >= 0 && a <= Lcyl) { rad <- c(p[2],p[3]); rho <- sqrt(sum(rad^2)); rad <- if (rho == 0) c(rr,0) else rad*rr/rho; c(a, rad) }
        else { cx <- if (a < 0) 0 else Lcyl; d <- p - c(cx,0,0); dd <- sqrt(sum(d^2)); d <- if (dd == 0) c(if (a<0) -rr else rr,0,0) else d*rr/dd; c(cx,0,0)+d }
      } else { u <- p/abc; abc * u / sqrt(sum(u^2)) }
    }
    geo_reflect <- function(p) {
      if (is_rod) { a <- p[1]
        if (a >= 0 && a <= Lcyl) { rad <- c(p[2],p[3]); rho <- sqrt(sum(rad^2)); if (rho > rr) rad <- rad*(2*rr - rho)/rho; c(a, rad) }
        else { cx <- if (a < 0) 0 else Lcyl; d <- p - c(cx,0,0); dd <- sqrt(sum(d^2)); if (dd > rr) d <- d*(2*rr - dd)/dd; c(cx,0,0)+d }
      } else { u <- p/abc; nu <- sqrt(sum(u^2)); if (nu > 1) u <- u*(2-nu)/nu; abc*u }
    }
    geo_rand_surf <- function() if (is_rod) { s <- runif(1,0,Lcyl); phi <- runif(1,0,2*pi); c(s, rr*cos(phi), rr*sin(phi)) } else geo_project(rnorm(3) * abc)
    geo_rand_in   <- function() if (is_rod) {
        repeat { p <- c(runif(1,-rr,Lcyl+rr), runif(1,-rr,rr), runif(1,-rr,rr)); a <- p[1]
          if ((if (a >= 0 && a <= Lcyl) (p[2]^2+p[3]^2) <= rr^2 else if (a < 0) sum(p^2) <= rr^2 else sum((p-c(Lcyl,0,0))^2) <= rr^2)) return(p) }
      } else { d <- rnorm(3); d <- d/sqrt(sum(d^2)); abc * d * runif(1)^(1/3) }

    # ── motion modifier (Phase 3): free | confined (OU tether) | directed (drift) |
    #    twostate (mobile↔immobile CTMC). Applied to the proposed 3D point, then the
    #    geometry constraint re-projects/reflects it. ──────────────────────────────
    motion  <- input$sim_trk_motion  %||% "free"
    # Anomalous diffusion: MSD grows as t^alpha via correlated fBm increments. D still
    # sets the single-cycle MSD, so only the growth with lag changes.
    alpha   <- min(1.95, max(0.05, as.numeric(input$sim_trk_alpha %||% 1)))
    is_anom <- identical(motion, "anomalous")
    # Hop diffusion: free inside a corral, reflecting off its barrier, crossing into the
    # neighbour with probability p_perm per encounter. Free at short lags, slower at long
    # ones, so alpha < 1 in between without any long-range memory.
    L_hop   <- max(1, as.numeric(input$sim_trk_Lhop %||% 60))
    p_perm  <- min(1, max(0, as.numeric(input$sim_trk_pperm %||% 0.05)))
    L_c     <- max(1, as.numeric(input$sim_trk_Lc %||% 150))         # confinement radius (nm)
    v_nmms  <- max(0, as.numeric(input$sim_trk_vel %||% 0))          # directed speed (µm/s = nm/ms)
    kswitch <- min(1, max(0, as.numeric(input$sim_trk_kswitch %||% 0.03)))  # per-cycle state-switch prob

    drift0  <- if (is_rod) c(1,0,0)                                          # rod drift ∥ long axis
               else if (identical(geom, "membrane")) { a <- runif(1, 0, 2*pi); c(cos(a), sin(a), 0) }  # in-plane
               else { v <- rnorm(3); v/sqrt(sum(v^2)) }

    is_free <- identical(geom, "free"); is_mem <- identical(geom, "membrane")
    # dur[i] = real elapsed time of the step leaving cycle i (ms). The molecule keeps
    # diffusing through the SW overhead, so SW must enter the step duration and not
    # only the clock; otherwise D is under-recovered by cyc/(cyc+SW) (~5%).
    gen_path <- function(n_cyc, D_trk) {
      surface <- identical(mode, "surface")
      tau_c   <- if (D_trk > 0) L_c^2 / D_trk else Inf                # OU relaxation time (ms)
      p <- if (is_free) c(0,0,0)
           else if (is_mem) c(0, 0, if (mthick > 0) runif(1, -mthick/2, mthick/2) else 0)
           else if (surface) geo_rand_surf() else geo_rand_in()
      home <- p; st <- 1L
      pos <- matrix(NA_real_, n_cyc, 3); pos[1,] <- p; states <- integer(n_cyc); states[1] <- 1L
      # Anomalous motion draws its increments up-front: fBm increments are CORRELATED,
      # so they cannot be generated one step at a time like Brownian ones.
      inc <- if (is_anom)
               cbind(.sim_fgn(n_cyc, alpha/2), .sim_fgn(n_cyc, alpha/2), .sim_fgn(n_cyc, alpha/2))
             else NULL
      # Recorded/dark and the clock are decided inside the loop: dur[i] (the elapsed
      # time of the step leaving cycle i) needs to know whether cycle i was recorded,
      # and a displacement-dependent dark probability needs the step that arrived at
      # cycle i. Deciding them in the loop keeps both exact.
      # ── BRIGHTNESS (photons/ms) over the molecule's life ────────────────────
      # log-rate = median + between-molecule offset − fade over the lifetime + an
      # AR(1) wobble, so brightness is correlated from one localization to the next.
      tau0  <- expfrac * cyc                                  # one pattern repeat (ms)
      lrate <- log(max(1e-6, photons / tau0))
      if (pht_btw > 0) lrate <- lrate + stats::rnorm(1, 0, pht_btw)
      frac  <- if (n_cyc > 1L) (seq_len(n_cyc) - 1L) / (n_cyc - 1L) else 0
      lr    <- lrate - pht_dec * frac
      if (pht_sd > 0) {
        e <- numeric(n_cyc); z <- stats::rnorm(n_cyc)
        e[1] <- pht_sd * z[1]
        if (n_cyc > 1L) for (i in 2:n_cyc)
          e[i] <- pht_phi * e[i-1] + sqrt(1 - pht_phi^2) * pht_sd * z[i]
        lr <- lr + e
      }
      rate <- exp(lr)                                          # photons per ms

      integ <- logical(n_cyc); integ[1] <- TRUE
      dur   <- rep(cyc, n_cyc); dur[1] <- cyc + sw
      eco_i <- rep(NA_real_, n_cyc); rep_i <- rep(NA_integer_, n_cyc)
      if (budget) {
        # Repeat the pattern until the photon limit is reached, up to tau_rep repeats.
        # A window that never reaches the limit is abandoned: those cycles are dark and
        # the gap simply grows. This is what makes eco floor at the limit exactly.
        integ[] <- FALSE
        i <- 1L
        while (i <= n_cyc) {
          pend <- 0; j <- 0L
          while (j < tau_rep && i + j <= n_cyc) {
            j <- j + 1L
            pend <- pend + stats::rpois(1, rate[i + j - 1L] * tau0)
            if (pend >= pht_min) break
          }
          if (pend >= pht_min) { integ[i] <- TRUE; eco_i[i] <- pend; rep_i[i] <- j }
          i <- i + j
        }
        dur <- cyc + sw * integ
      }
      for (i in 2:n_cyc) {
        Dstep <- if (identical(motion,"twostate") && st == 0L) D_trk * 0.01 else D_trk
        dti   <- dur[i - 1L]                                                # this step's duration (ms)
        sdi   <- sqrt(2 * Dstep * dti)
        g0 <- if (is_anom) inc[i, ] * sdi else rnorm(3, 0, sdi)             # correlated when anomalous
        if (is_free)      { p2 <- p + g0 }                                  # unbounded 3D
        else if (is_mem)  { p2 <- p + c(g0[1:2], 0) }                        # 2D lateral (z handled below)
        else if (surface) { nn <- geo_normal(p); p2 <- p + (g0 - sum(g0*nn)*nn) }  # tangent step
        else              { p2 <- p + g0 }                                  # 3D volume
        if (identical(motion,"directed")) p2 <- p2 + drift0 * v_nmms * dti
        if (identical(motion,"confined")) p2 <- home + (p2 - home) * exp(-dti / tau_c)  # OU pull to home
        if (identical(motion,"hop")) {                       # corral with a permeable barrier
          dv <- p2 - home; nd <- sqrt(sum(dv^2))
          if (nd > L_hop) {
            if (runif(1) < p_perm) home <- home + (dv / nd) * (2 * L_hop)   # cross into the neighbour
            else                   p2   <- home + (dv / nd) * (2 * L_hop - nd)  # reflect off the barrier
          }
        }
        p  <- if (is_free) p2
              else if (is_mem) c(p2[1], p2[2],                              # keep z within the membrane thickness
                                 if (mthick > 0) max(-mthick/2, min(mthick/2, p[3] + rnorm(1, 0, sdi * 0.2))) else 0)
              else if (surface) geo_project(p2) else geo_reflect(p2)
        if (identical(motion,"twostate") && runif(1) < kswitch) st <- 1L - st
        pos[i,] <- p; states[i] <- st
        # Was this cycle recorded? Under the photon budget that was already decided from
        # brightness above; otherwise fall back to a flat per-cycle probability.
        if (!budget) {
          integ[i] <- runif(1) > p_dark
          dur[i]   <- cyc + sw * integ[i]
        }
      }
      # ── MOTION BLUR ──────────────────────────────────────────────────────────
      # A MINFLUX localization is the time-average of the path over its cycle, not
      # the position at the start of it. For Brownian motion the average over one
      # cycle, conditioned on both end points, is the midpoint plus bridge noise of
      # variance D·cyc/6 per axis. That is the R = 1/6 blur factor the noise-corrected
      # CVE / MSD estimators are written for. Computed here, in the canonical frame,
      # so the bridge noise can be kept in the surface tangent plane; the last cycle
      # has no successor and stays unblurred.
      bp <- NULL
      if (blur) {
        t_cyc <- c(0, cumsum(dur)[-n_cyc])
        # Each recorded cycle is averaged over its own window, ending at the NEXT
        # recorded cycle. A gap makes both the window (up to tau_rep repeats) and the
        # step longer, which is the structure real files show through tau = eco/efo.
        rec  <- which(integ)
        nxt  <- seq_len(n_cyc); Tstep <- rep(cyc, n_cyc); taui <- rep(expfrac * cyc, n_cyc)
        if (length(rec) > 1L) {
          nxt[rec[-length(rec)]]   <- rec[-1]
          Tstep[rec[-length(rec)]] <- diff(t_cyc[rec])
          reps_used <- if (budget) rep_i[rec[-length(rec)]] else pmin(diff(rec), tau_rep)
          taui[rec[-length(rec)]]  <- reps_used * expfrac * cyc
        }
        mid <- pos + (taui / (2 * Tstep)) * (pos[nxt, , drop = FALSE] - pos)
        Dst <- if (identical(motion, "twostate")) ifelse(states == 0L, D_trk * 0.01, D_trk)
               else rep(D_trk, n_cyc)
        # Time-average over a window tau <= T, conditioned on both end points:
        #   mean = x_i + (tau/(2T))*(x_j - x_i),  var = (2*D*tau/3)*(1 - 3*tau/(4T))
        # which reduces to the midpoint plus D*T/6 when tau = T (R = 1/6).
        br  <- matrix(rnorm(n_cyc * 3), ncol = 3) *
               sqrt(pmax(0, (2 * Dst * taui / 3) * (1 - 3 * taui / (4 * Tstep))))
        if (is_mem) br[, 3] <- 0                                    # in-plane only
        else if (surface && !is_free) {                             # tangent plane only
          nrm <- t(vapply(seq_len(n_cyc), function(i) geo_normal(pos[i, ]), numeric(3)))
          br  <- br - (rowSums(br * nrm)) * nrm
        }
        bp <- mid + br
        bp[n_cyc, ] <- pos[n_cyc, ]                  # no successor: unblurred
      }
      if (is_rod) { pos <- pos[, rod_perm, drop = FALSE]
                    if (!is.null(bp)) bp <- bp[, rod_perm, drop = FALSE] }
      if (is_mem) { pos <- pos %*% t(R_mem)            # tilt the flat local sheet into 3D
                    if (!is.null(bp)) bp <- bp %*% t(R_mem) }
      attr(pos, "states")  <- states
      attr(pos, "blurpos") <- bp
      attr(pos, "integ")   <- integ
      attr(pos, "eco")     <- eco_i
      attr(pos, "reps")    <- rep_i
      # Start time of each cycle: cycles last `cyc`, plus SW once per RECORDED cycle, so
      # two recorded cycles k apart are exactly k·cyc + SW apart.
      attr(pos, "t_cyc")   <- c(0, cumsum(dur)[-n_cyc])
      pos
    }
    # Survivorship: a molecule that moves > R_reacq during a dark gap is lost - the
    # track ends there. Reproduces the dark-cycle bias (long jumps culled → apparent
    # confinement). Off unless enabled.
    surv    <- isTRUE(input$sim_trk_surv)
    R_reacq <- max(1, as.numeric(input$sim_trk_rreacq %||% 300))

    rows <- vector("list", n_trk); gt <- vector("list", n_trk); k <- 0L
    withProgress(message = "Simulating tracks…", value = 0, {
    for (ti in seq_len(n_trk)) {
      if (ti %% 20L == 0L) incProgress(20 / n_trk)
      # mean_len = target MEDIAN recorded localizations/track. Recorded count is ≈
      # geometric, whose median ≈ ln(2)·mean, and only (1−p_dark) of cycles are recorded,
      # so the bleach lifetime in CYCLES has mean = median / (ln2 · (1−p_dark)) to land
      # the recorded MEDIAN on mean_len despite dark-cycling.
      # D and lifetime are drawn as a correlated pair, so a sample where the slow
      # molecules are the ones tracked longest can be reproduced.
      z1 <- stats::rnorm(1)
      z2 <- D_life * z1 + sqrt(1 - D_life^2) * stats::rnorm(1)
      life  <- mean_len / (log(2) * max(0.05, 1 - p_dark))
      if (life_sd > 0) life <- life * exp(life_sd * z2)
      n_cyc <- max(3L, stats::rgeom(1, 1 / max(2, life)) + 1L)
      # Each molecule may carry its own D (lognormal, median = the D input).
      D_trk <- if (D_sdlog > 0) D_nm2ms * exp(D_sdlog * z1) else D_nm2ms
      pos   <- gen_path(n_cyc, D_trk)                   # ground-truth path at cycle resolution
      integ <- attr(pos, "integ"); t_cyc <- attr(pos, "t_cyc")
      idx   <- which(integ)
      if (surv && length(idx) > 1L) {                   # survivorship: cut at first over-range dark-gap jump
        jmp  <- sqrt(rowSums((pos[idx[-1], , drop = FALSE] - pos[idx[-length(idx)], , drop = FALSE])^2))
        lost <- which(jmp > R_reacq)
        if (length(lost)) idx <- idx[seq_len(lost[1])]  # keep up to (and incl.) the pre-jump loc
      }
      if (length(idx) < 5L) next
      # Photons per localization: the budget already counted them (floored at the limit
      # by construction); otherwise draw them from the brightness model directly.
      N <- if (budget) attr(pos, "eco")[idx] else {
             r <- if (pht_sd > 0) photons * exp(stats::rnorm(length(idx), 0, pht_sd)) else photons
             pmax(pht_min, stats::rpois(length(idx), r))
           }
      tau_i <- if (budget) attr(pos, "reps")[idx] * expfrac * cyc else rep(expfrac * cyc, length(idx))
      # MINFLUX CRLB for a hexagonal TCP of diameter L with N detected photons.
      # Cross-checked on real tracking data: L = 100.8 nm (patGeoFactor 0.28) and the
      # median eco of ~38 photons give 4.1 nm, exactly the blur-aware CVE sigma of the
      # GP67H tracking files.
      sig   <- L_tcp / (TCP_CRLB_PREFACTOR * sqrt(N))  # per-localization σ_loc (nm)
      bpos  <- attr(pos, "blurpos")                    # blurred path, or NULL when off
      true_pos <- if (!is.null(bpos)) bpos[idx, , drop = FALSE] else pos[idx, , drop = FALSE]
      obs   <- true_pos + matrix(rnorm(length(idx) * 3), ncol = 3) * sig
      # ── RE-ACQUISITION SHRINKAGE ───────────────────────────────────────────────
      # After a dark gap the molecule sits off the pattern centre, where the estimator
      # is biased back toward it, so the first localization after a gap is reported
      # closer to the previous one than it really is. Sequential: the centre is the
      # PREVIOUS REPORTED position, which is all the instrument has.
      if (shrink > 0 && length(idx) > 1L) {
        gap <- c(FALSE, diff(idx) > 1L)                # first loc after a dark gap
        for (j in which(gap)) obs[j, ] <- obs[j - 1L, ] + (1 - shrink) * (obs[j, ] - obs[j - 1L, ])
      }
      # Place the track on its cell: rotate in-plane, then offset to the cell centre.
      if (n_cells > 1L) {
        ci <- 1L + (ti - 1L) %% n_cells
        ca <- cos(cell_ang[ci]); sa <- sin(cell_ang[ci])
        ctr <- if (identical(geom, "rod")) c(Lcyl/2, 0, 0)[rod_perm] else c(0,0,0)
        ox <- obs[, 1] - ctr[1]; oy <- obs[, 2] - ctr[2]
        obs[, 1] <- ox * ca - oy * sa + cell_pos[ci, 1]
        obs[, 2] <- ox * sa + oy * ca + cell_pos[ci, 2]
        obs[, 3] <- obs[, 3] - ctr[3]
      }
      k <- k + 1L
      rows[[k]] <- data.frame(
        tid = k, tim = t_cyc[idx] / 1000,              # SECONDS (importer ×1000 → ms)
        x = obs[, 1], y = obs[, 2], z = obs[, 3],
        efo = N / (tau_i / 1000), eco = N, stringsAsFactors = FALSE)
      gt[[k]] <- data.frame(tid = k, D_true_um2_s = D_trk / 1000, geometry = geom, mode = mode,
                            motion = motion, L_c_nm = if (identical(motion,"confined")) L_c else NA_real_,
                            v_um_s = if (identical(motion,"directed")) v_nmms else NA_real_,
                            n_locs = length(idx), sigma_loc_med = stats::median(sig),
                            stringsAsFactors = FALSE)
    }
    })
    if (k == 0L) return(NULL)
    df <- dplyr::bind_rows(rows[seq_len(k)])
    if (proj2d) df$z <- NULL
    list(df = df, gt = dplyr::bind_rows(gt[seq_len(k)]),
         D_um2_s = D_um2_s, cyc = cyc, sw = sw, blur = blur, p_dark = p_dark, proj2d = proj2d,
         expfrac = expfrac, D_cv = D_cv, shrink = shrink,
         life_cv = life_cv, pht_min = pht_min, pht_sd = pht_sd, n_cells = n_cells, fov_nm = fov_nm,
         tau_rep = tau_rep, budget = budget, alpha = if (is_anom) alpha else NA_real_,
         pht_dec = pht_dec, L_hop = L_hop, p_perm = p_perm,
         D_life = D_life,
         mode = mode, d_mot = d_mot, geom = geom,
         motion = motion, L_c = L_c, v_nmms = v_nmms, surv = surv, R_reacq = R_reacq)
  })

  # CVE round-trip: recover D from consecutive single-cycle (Δt = cyc) steps, using the
  # TOTAL displacement with the MOTION dimensionality. D = (⟨|Δr|²⟩ + 2⟨Δr_i·Δr_{i+1}⟩)
  # / (2·d·Δt). SURFACE diffusion has d=2 DOF but spreads over 3 Cartesian axes (a
  # per-axis CVE would under-estimate D by ⅓); VOLUME diffusion has d=3. When projected
  # to 2D only n_cart=2 axes are seen, so d = min(motion-DOF, n_cart). n_cart Cartesian
  # axes carry the noise (cov = −n_cart·σ²).
  sim_trk_cve <- reactive({
    res <- sim_trk_result(); if (is.null(res)) return(NULL)
    df <- res$df; cyc <- res$cyc
    dt_step <- cyc + (res$sw %||% 0)                    # true single-cycle Δt = cycle + SW
    axes  <- c("x", "y", if (!is.null(df$z)) "z")
    ncart <- length(axes)
    d_mot <- min(res$d_mot %||% 2L, ncart)              # 2 surface / 3 volume, capped by observed axes
    r2_sum <- 0; r2_n <- 0L; cov_sum <- 0; cov_n <- 0L
    for (tb in split(df, df$tid)) {
      tb  <- tb[order(tb$tim), ]
      dt  <- c(NA_real_, diff(tb$tim) * 1000)           # ms
      one <- is.finite(dt) & abs(dt - dt_step) < cyc * 0.25  # single-cycle step (ending at this row)
      # per-step displacement vector (rows = steps ending at row i)
      dmat <- sapply(axes, function(ax) c(NA_real_, diff(tb[[ax]])))
      if (is.null(dim(dmat))) dmat <- matrix(dmat, ncol = ncart)
      r2  <- rowSums(dmat^2)                             # |Δr|² total
      s1  <- one & is.finite(r2)
      r2_sum <- r2_sum + sum(r2[s1]); r2_n <- r2_n + sum(s1)
      nxt  <- c(s1[-1], FALSE); pair <- s1 & nxt         # consecutive single-cycle pair
      dot  <- rowSums(dmat * rbind(dmat[-1, , drop = FALSE], NA_real_))  # Δr_i · Δr_{i+1}
      cov_sum <- cov_sum + sum(dot[pair], na.rm = TRUE); cov_n <- cov_n + sum(pair)
    }
    if (r2_n < 1L) return(NULL)
    r2m  <- r2_sum / r2_n                                # ⟨|Δr|²⟩ (nm²)
    cov1 <- if (cov_n > 0L) cov_sum / cov_n else 0       # ⟨Δr_i·Δr_{i+1}⟩ = −n_cart·σ²
    # Vestergaard's CVE: D is blur-independent, σ is not. With motion blur on, the
    # localization is a time-average over the collection window τ, so R = τ/(6·Δt) and
    # σ² = R⟨Δx²⟩ + (2R−1)⟨Δx_i Δx_{i+1}⟩ per axis; σ² = −cov is the R = 0 (snapshot)
    # special case. Reading it with the wrong R is what makes a real 4 nm σ_loc look
    # like 1.4 nm.
    Rblur <- if (isTRUE(res$blur)) (res$expfrac %||% 1) * cyc / (6 * dt_step) else 0
    # (single-cycle steps only, so the window is one repeat regardless of tau_rep)
    sig2  <- Rblur * (r2m / ncart) + (2 * Rblur - 1) * (cov1 / ncart)
    list(D_um2_s  = ((r2m + 2 * cov1) / (2 * d_mot * dt_step)) / 1000,   # CVE D (µm²/s)
         sigma_nm = sqrt(max(0, sig2)), R_blur = Rblur, n_pairs = cov_n)
  })

  # ── sim_trk_kcurve() - displacement vs dark-gap length ──────────────────────
  # GOAL: the diagnostic the new sampling knobs exist to drive. Mean per-axis ⟨Δx²⟩
  # against the number of cycles k a step spans, next to the free-diffusion
  # expectation. Free Brownian motion grows linearly in k; the real GP67H files
  # flatten out after k = 2 (see docs/tracking_sigma_and_cve_plan.md §2.3), and
  # re-acquisition shrinkage / speed-dependent dark cycles are the candidate causes
  # this table is meant to discriminate.
  # returns: data.frame(k, n, msq, expect, pct) or NULL.
  sim_trk_kcurve <- reactive({
    res <- sim_trk_result(); if (is.null(res)) return(NULL)
    df  <- res$df; cyc <- res$cyc; sw <- res$sw %||% 0
    axes  <- c("x", "y", if (!is.null(df$z)) "z")
    ncart <- length(axes)
    d_mot <- min(res$d_mot %||% 2L, ncart)
    kk <- list(); dd <- list()
    for (tb in split(df, df$tid)) {
      tb <- tb[order(tb$tim), ]; if (nrow(tb) < 2L) next
      dt <- diff(tb$tim) * 1000
      kk[[length(kk) + 1L]] <- rep(pmax(1L, round((dt - sw) / cyc)), ncart)
      dd[[length(dd) + 1L]] <- unlist(lapply(axes, function(a) diff(tb[[a]])))
    }
    if (!length(kk)) return(NULL)
    k <- unlist(kk); dx <- unlist(dd)
    ok <- is.finite(k) & is.finite(dx) & k >= 1L & k <= 6L
    if (!any(ok)) return(NULL)
    agg <- stats::aggregate(list(msq = dx[ok]^2), by = list(k = k[ok]),
                            FUN = function(v) c(m = mean(v), n = length(v)))
    out <- data.frame(k = agg$k, msq = agg$msq[, "m"], n = agg$msq[, "n"])
    # Free-diffusion expectation per axis. Surface motion has d_mot DOF spread over
    # ncart Cartesian axes, so each axis carries d_mot/ncart of the motion. The blur
    # term −(2/3)·D·τ is subtracted because the simulated positions are time-averages.
    w_gt  <- pmax(0, res$gt$n_locs - 1)
    D_ref <- (if (sum(w_gt) > 0) stats::weighted.mean(res$gt$D_true_um2_s, w_gt) else res$D_um2_s) * 1000
    sig2  <- stats::median(res$gt$sigma_loc_med, na.rm = TRUE)^2
    tau   <- if (isTRUE(res$blur)) (res$expfrac %||% 1) * cyc else 0
    fr    <- d_mot / ncart
    out$expect <- 2 * D_ref * fr * (out$k * cyc + sw) + 2 * sig2 - (2/3) * D_ref * fr * tau
    out$pct    <- 100 * out$msq / out$expect
    out[order(out$k), ]
  })

  output$sim_trk_plot <- renderPlotly({
    res <- sim_trk_result()
    shiny::validate(shiny::need(!is.null(res), "Set parameters and click 'Simulate tracks'."))
    df <- res$df; is3d <- !is.null(df$z)
    show <- head(unique(df$tid), 60L); d <- df[df$tid %in% show, ]
    if (is3d)
      plotly::plot_ly(d, x = ~x, y = ~y, z = ~z, split = ~factor(tid),
                      type = "scatter3d", mode = "lines",
                      line = list(width = 2), showlegend = FALSE) %>%
        plotly::layout(title = "Simulated MINFLUX tracks (first 60)",
          scene = list(xaxis = list(title = "x (nm)"), yaxis = list(title = "y (nm)"),
                       zaxis = list(title = "z (nm)"),
                       aspectmode = "data"))   # equal nm-per-unit on X/Y/Z (no axis stretching)
    else
      plotly::plot_ly(d, x = ~x, y = ~y, split = ~factor(tid),
                      type = "scatter", mode = "lines",
                      line = list(width = 1), showlegend = FALSE) %>%
        plotly::layout(title = "Simulated MINFLUX tracks - 2D projection (first 60)",
          xaxis = list(title = "x (nm)"),
          yaxis = list(title = "y (nm)", scaleanchor = "x"))
  })

  output$sim_trk_summary <- renderUI({
    res <- sim_trk_result(); if (is.null(res)) return(NULL)
    df <- res$df; cve <- sim_trk_cve()
    sig_true <- stats::median(res$gt$sigma_loc_med, na.rm = TRUE)   # ground-truth σ_loc (nm)
    dt_all <- unlist(lapply(split(df$tim, df$tid), function(t) diff(sort(t)) * 1000), use.names = FALSE)
    frac_int <- if (length(dt_all))
                  mean(abs(dt_all - (res$cyc + (res$sw %||% 0))) < res$cyc * 0.25) else NA_real_
    # With a heterogeneous population the pooled CVE cannot recover the median D: it
    # recovers the STEP-WEIGHTED MEAN of the drawn per-track D's, so that is what the
    # round-trip is scored against (the two coincide when the spread is 0).
    w_gt  <- pmax(0, res$gt$n_locs - 1)
    D_ref <- if (sum(w_gt) > 0) stats::weighted.mean(res$gt$D_true_um2_s, w_gt) else res$D_um2_s
    recov <- if (!is.null(cve)) {
      rr  <- cve$D_um2_s / D_ref
      col <- if (is.finite(rr) && abs(rr - 1) < 0.25) "#228833" else "#BB7700"
      htmltools::tags$span(style = sprintf("color:%s; font-weight:bold;", col),
        sprintf("%.4g µm²/s  (%.0f%% of truth)", cve$D_um2_s, 100 * rr))
    } else "-"
    motion_lbl <- switch(res$motion %||% "free",
      anomalous = sprintf("anomalous fBm (α = %.2f)", res$alpha),
      hop = sprintf("hop diffusion (corral %.0f nm, permeability %.2f)", res$L_hop, res$p_perm),
      confined = sprintf("confined (OU tether, L_c = %.0f nm)", res$L_c),
      directed = sprintf("directed (v = %.3g µm/s)", res$v_nmms),
      twostate = "two-state (mobile↔immobile)",
      "free")
    mode_lbl <- if ((res$geom %||% "") %in% c("free","membrane")) ""
                else if (identical(res$mode, "volume")) " · volume (intracellular)" else " · surface"
    geom_lbl <- sprintf("%s%s · %s%s", res$geom %||% "?", mode_lbl,
                        motion_lbl, if (isTRUE(res$surv)) sprintf(" · survivorship R=%.0fnm", res$R_reacq) else "")
    htmltools::div(style = "font-size:0.9em; color:#333; margin-top:8px;",
      htmltools::tags$table(style = "border-collapse:collapse;",
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Model"),
          htmltools::tags$td(htmltools::tags$b(geom_lbl))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Tracks / localizations"),
          htmltools::tags$td(sprintf("%d / %s", length(unique(df$tid)), format(nrow(df), big.mark=",")))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Ground-truth D"),
          htmltools::tags$td(htmltools::tags$b(
            if ((res$D_cv %||% 0) > 0)
              sprintf("%.4g µm²/s median · CV %.2f · step-weighted mean %.4g · per-track 10–90%% %.3g–%.3g",
                      res$D_um2_s, res$D_cv, D_ref,
                      stats::quantile(res$gt$D_true_um2_s, 0.1), stats::quantile(res$gt$D_true_um2_s, 0.9))
            else sprintf("%.4g µm²/s", res$D_um2_s)))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Recovered D (CVE round-trip)"),
          htmltools::tags$td(recov)),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Geometric bias factor"),
          htmltools::tags$td(if (!is.null(cve) && is.finite(cve$D_um2_s) && res$D_um2_s > 0) {
            bf <- cve$D_um2_s / D_ref
            htmltools::tags$span(style = if (abs(bf - 1) > 0.05) "color:#b26a00;font-weight:bold;" else "",
              sprintf("%.2f×  →  divide apparent D by %.2f to correct", bf, bf))
          } else "-")),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "σ_loc  truth / recovered"),
          htmltools::tags$td(sprintf("%.2f nm / %s", sig_true,
            if (!is.null(cve)) sprintf("%.2f nm  (CVE, R = %.2g)", cve$sigma_nm, cve$R_blur %||% 0) else "-"))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Sampling clock"),
          htmltools::tags$td(sprintf("Δt = k × %.4g ms + %.3g ms (SW) · blur %s%s%s",
            res$cyc, res$sw %||% 0,
            if (isTRUE(res$blur))
              sprintf("on (τ = %.0f%% of cycle, R = %.3g)", 100 * (res$expfrac %||% 1),
                      (res$expfrac %||% 1) * res$cyc / (6 * (res$cyc + (res$sw %||% 0))))
            else "off (snapshot)",
            if ((res$shrink %||% 0) > 0) sprintf(" · gap shrinkage %.0f%%", 100 * res$shrink) else "", ""))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Photons / dark cycles"),
          htmltools::tags$td(sprintf("eco median %.0f (limit %d) · %s",
            stats::median(df$eco, na.rm = TRUE), res$pht_min %||% 1L,
            if (isTRUE(res$budget))
              sprintf("dark %.0f%%, emergent from the photon budget", 100 * (1 - frac_int))
            else sprintf("dark %.0f%% as set", 100 * res$p_dark)))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Integrated (single-cycle) fraction"),
          htmltools::tags$td(if (is.finite(frac_int)) sprintf("%.0f%%", 100 * frac_int) else "-")),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Median Δt / median track length"),
          htmltools::tags$td(sprintf("%.3g ms / %.0f locs",
            if (length(dt_all)) stats::median(dt_all) else NA_real_, stats::median(as.numeric(table(df$tid))))))),
      htmltools::tags$p(style = "color:#888; font-size:0.9em; margin-top:6px;",
        if (identical(res$motion, "confined")) {
          tau_c   <- if (res$D_um2_s > 0) (res$L_c^2 / (res$D_um2_s * 1000)) else Inf   # ms
          med_dur <- stats::median(unlist(lapply(split(df$tim, df$tid), function(t) (max(t) - min(t)) * 1000)))
          verdict <- if (is.finite(tau_c) && tau_c < 0.5 * med_dur) "✓ RESOLVABLE" else "⚠ NOT resolvable - looks free"
          sprintf("Confined (OU, L_c=%.0f nm): CVE gives the LOCAL D at single-cycle lags (≈ truth); confinement is an MSD plateau ≈ 2·d·L_c² in the Tracking tab. Resolvability: τ_c = L_c²/D ≈ %.3g ms vs median track duration ≈ %.3g ms → %s. (If τ_c ≳ track duration, free diffusion over a short track mimics the same extent - the key artifact to test.)",
                  res$L_c, tau_c, med_dur, verdict)
        }
        else if (identical(res$motion, "directed"))
          "Directed: the CVE local D reflects the diffusive part; the drift shows as an upward-curving MSD / a sustained drift-velocity in the Tracking tab."
        else if (identical(res$motion, "twostate"))
          "Two-state: mixed mobile/immobile steps - the CVE gives a population-average D; load in the Tracking tab and fit the HMM to recover the two states."
        else
          "Free diffusion: the CVE round-trip should land near 100 % of truth - that is the validation loop.",
        " Download the tracks (CSV) and load them in the Tracking tab to test the full MSD / confinement pipeline against this known ground truth."),
      # Displacement-vs-gap-length table: the k-curve diagnostic of the plan doc.
      local({
        kc <- tryCatch(sim_trk_kcurve(), error = function(e) NULL)
        if (is.null(kc) || !nrow(kc)) return(NULL)
        cell <- function(x, b = FALSE) htmltools::tags$td(
          style = "padding:1px 10px 1px 0; text-align:right;", if (b) htmltools::tags$b(x) else x)
        htmltools::div(style = "margin-top:8px;",
          htmltools::tags$div(style = "color:#666; font-size:0.86em;",
            "Per-axis ⟨Δx²⟩ (nm²) by dark-gap length k"),
          htmltools::tags$table(style = "border-collapse:collapse; font-size:0.84em;",
            htmltools::tags$tr(cell("k"), lapply(kc$k, function(v) cell(v, TRUE))),
            htmltools::tags$tr(cell("steps"), lapply(kc$n, function(v) cell(format(v, big.mark = ",")))),
            htmltools::tags$tr(cell("⟨Δx²⟩"), lapply(kc$msq, function(v) cell(sprintf("%.0f", v)))),
            htmltools::tags$tr(cell("free expectation"), lapply(kc$expect, function(v) cell(sprintf("%.0f", v)))),
            htmltools::tags$tr(cell("% of free"), lapply(kc$pct, function(v)
              htmltools::tags$td(style = sprintf("padding:1px 10px 1px 0; text-align:right; %s",
                if (is.finite(v) && abs(v - 100) > 15) "color:#b26a00;font-weight:bold;" else ""),
                sprintf("%.0f%%", v))))),
          htmltools::tags$p(style = "color:#888; font-size:0.84em; margin-top:3px;",
            "Free diffusion holds this row near 100% at every k. The GP67H tracking files fall to ",
            "~75% at k=2 and ~30% by k=6, which is what the gap-shrinkage and speed-dependent-dark ",
            "controls are for."))
      }),
      # Projection artifact: 2D imaging of a CURVED-surface protein under-measures the
      # azimuthal motion (its out-of-plane component is dropped), so the recovered D is a
      # fixed fraction of truth (~75% for a side-on cylinder) - NOT a bug, and the same
      # bias applies to real 2D data. The factor is itself a useful geometric calibration.
      if ((res$geom %||% "") %in% c("sphere","rod","ellipsoid") && identical(res$mode, "surface") && isTRUE(res$proj2d))
        htmltools::tags$p(style = "color:#b26a00; font-size:0.9em; margin-top:4px;",
          htmltools::HTML("&#9888; "),
          htmltools::tags$b("2D projection of a curved surface:"),
          " recovered D is expected BELOW truth (≈75% for a side-on cylinder) because the ",
          "out-of-plane part of the azimuthal motion is lost - this is a real geometric bias that ",
          "also affects your real 2D data. Use this ratio as the correction factor, or acquire in 3D. ",
          "It applies to free AND confined alike, so it is not the confinement signal."))
  })

  output$sim_trk_dl_csv <- downloadHandler(
    filename = function() paste0("sim_minflux_tracks_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      res <- sim_trk_result()
      shiny::validate(shiny::need(!is.null(res), "Simulate tracks first."))
      write.csv(res$df, file, row.names = FALSE)
    })
  output$sim_trk_dl_gt <- downloadHandler(
    filename = function() paste0("sim_ground_truth_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      res <- sim_trk_result()
      shiny::validate(shiny::need(!is.null(res), "Simulate tracks first."))
      write.csv(res$gt, file, row.names = FALSE)
    })

}
