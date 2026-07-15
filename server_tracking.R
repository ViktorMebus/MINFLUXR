## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
## server_tracking.R   -   Tracking analysis server logic
## All input/output IDs prefixed with trk_ to avoid collisions with SMLM IDs.
## Called from server.R inside the combined server function.

tracking_server <- function(input, output, session) {
  
  
  # --- Reactives / state ---
  rv <- reactiveValues(
    prev_tid_idx          = NULL,
    play_gen              = 0L,   # global playback generation (simultaneous mode)
    manifest_dir          = NULL, # current working chunks (may be filtered)
    original_manifest_dir = NULL, # frozen at upload — Filter & Rechunk always reads from here
    bigfile_notice        = FALSE,
    chunking              = FALSE,
    trk_diag              = NULL  # list(data_clusters, img_blobs) in data coords
  )
  is_playing         <- reactiveVal(FALSE)
  tick               <- reactiveVal(0L)
  slider_guard       <- reactiveVal(0L)
  trace_map          <- reactiveVal(NULL)
  current_tid_idx    <- reactiveVal(1L)
  # Frozen at upload time; never overwritten by chunk-switching or Filter & Rechunk.
  # Provides stable grey-bar reference for the metadata histogram.
  trk_raw_stats_rv      <- reactiveVal(NULL)
  # Updated after each Filter & Rechunk; holds stats across ALL remaining TIDs
  # (not just the current chunk).  NULL means no rechunk has happened yet.
  # Filter UIs use this as their value when non-NULL; Reset Ranges clears it.
  trk_rechunk_stats_rv  <- reactiveVal(NULL)
  # Full spatial extent of all remaining TIDs after rechunk.
  # list(x = c(xmin, xmax), y = c(ymin, ymax)) — used to lock plot axes.
  trk_rechunk_bbox      <- reactiveVal(NULL)
  current_pt_idx     <- reactiveVal(1L)
  # Current plot axis ranges captured from plotly_relayout — used by "Crop to view".
  # list(x = c(x0, x1), y = c(y0, y1)) or NULL before first zoom interaction.
  trk_zoom_ranges    <- reactiveVal(NULL)
  # Active spatial box filter — NULL means no spatial filter is applied.
  # list(x = c(x0, x1), y = c(y0, y1), pct_thr = numeric[0-1])
  trk_spatial_filter <- reactiveVal(NULL)

  # "TSA has been run at least once" gate. The Integrated/Dark bundle is
  # eventReactive-gated on trk_tsa_run, but the sibling display outputs (Δt
  # histogram, stats, quality verdict, per-track progression) read dt_vectors()
  # directly — this flag keeps them dormant until Run is pressed. Cleared when TSA
  # is disabled; stays set across filter changes (re-press Run to refresh, same as
  # the other run-gated analyses).
  tsa_ran <- reactiveVal(FALSE)
  # Run TSA (and Estimate σ_loc below) also APPLY the current cleaning settings, so
  # the raw→analyzed diagnostics reflect split/trim/continuity without a separate
  # Plot Data press. bump_clean() is a no-op when cleaning is unchanged.
  observeEvent(input$trk_tsa_run, { bump_clean(); tsa_ran(TRUE) })
  observeEvent(input$trk_cve_run, { bump_clean() }, ignoreInit = TRUE)
  observeEvent(input$trk_show_tsa, {
    if (!isTRUE(input$trk_show_tsa)) tsa_ran(FALSE)
  }, ignoreInit = TRUE)

  # Same run-gate for the Instantaneous Speed section: both the TID selector
  # (which ranks all tracks by export strength — heavy) and the plot stay dormant
  # until Run is pressed, then remain interactive on TID selection.
  speed_ran <- reactiveVal(FALSE)
  observeEvent(input$trk_speed_run, { speed_ran(TRUE) })
  observeEvent(input$trk_show_speed, {
    if (!isTRUE(input$trk_show_speed)) speed_ran(FALSE)
  }, ignoreInit = TRUE)

  # ── App-wide localization precision σ_loc (the short-lag noise floor) ────────
  # Two independent sources feed the app-wide value, in priority order:
  #   1. manual override (input$trk_sigma_loc_nm > 0)
  #   2. the Calibration panel's CVE estimate (MSD-free; trk_sigma_loc_cve)
  #   3. the unified MSD fit's intercept (trk_sigma_loc_rv)
  # Every short-lag tool (D-map, JDD, confinement, HMM, directionality, speed)
  # reads trk_sigma_loc_eff() so the noise floor is consistent everywhere.
  trk_sigma_loc_rv  <- reactiveVal(NA_real_)   # nm, from the last unified MSD fit
  trk_sigma_loc_cve <- reactiveVal(NA_real_)   # nm, from the CVE calibration
  trk_sigma_loc_eff <- reactive({
    manual <- suppressWarnings(as.numeric(input$trk_sigma_loc_nm %||% 0))
    if (is.finite(manual) && manual > 0) return(manual)
    cve <- trk_sigma_loc_cve()
    if (is.finite(cve) && cve > 0) return(cve)
    s <- trk_sigma_loc_rv()
    if (is.finite(s) && s > 0) s else NA_real_
  })
  # Which source is currently active (for status readouts).
  trk_sigma_loc_src <- reactive({
    manual <- suppressWarnings(as.numeric(input$trk_sigma_loc_nm %||% 0))
    if (is.finite(manual) && manual > 0) "manual override"
    else if (is.finite(trk_sigma_loc_cve())) "calibration (CVE)"
    else if (is.finite(trk_sigma_loc_rv()))  "unified MSD fit"
    else "not set"
  })

  # Mean + SD + t-based 95% CI treating each element (a TRACK) as the independent
  # unit — the honest CI for track-level statistics, since localizations within a
  # track are correlated (the track, not the localization, is the sampling unit).
  # Also returns median + quartiles for a robust summary.
  .mean_ci <- function(x, conf = 0.95) {
    x <- x[is.finite(x)]; n <- length(x)
    if (n < 1L) return(list(mean = NA_real_, sd = NA_real_, se = NA_real_,
                            lo = NA_real_, hi = NA_real_, n = 0L,
                            median = NA_real_, q1 = NA_real_, q3 = NA_real_))
    m  <- mean(x); md <- stats::median(x)
    qs <- stats::quantile(x, c(0.25, 0.75), names = FALSE, na.rm = TRUE)
    if (n < 2L) return(list(mean = m, sd = NA_real_, se = NA_real_, lo = NA_real_,
                            hi = NA_real_, n = n, median = md, q1 = qs[1], q3 = qs[2]))
    s  <- stats::sd(x); se <- s / sqrt(n)
    tc <- stats::qt(1 - (1 - conf) / 2, df = n - 1)
    list(mean = m, sd = s, se = se, lo = m - tc * se, hi = m + tc * se,
         n = n, median = md, q1 = qs[1], q3 = qs[2])
  }

  # Percentile bootstrap CI of a statistic over the sampling units (TRACKS): x is a
  # per-track vector, resampled with replacement B times. This is the honest way to
  # put a CI on an ensemble number, because localizations within a track are
  # correlated — the track is the independent unit, not the localization.
  .boot_ci <- function(x, statfn = stats::median, B = 2000L, conf = 0.95) {
    x <- x[is.finite(x)]; n <- length(x)
    if (n < 3L) return(c(lo = NA_real_, hi = NA_real_))
    bs <- vapply(seq_len(B), function(i) statfn(x[sample.int(n, n, replace = TRUE)]), numeric(1))
    q  <- stats::quantile(bs, c((1 - conf)/2, 1 - (1 - conf)/2), names = FALSE, na.rm = TRUE)
    c(lo = q[1], hi = q[2])
  }

  # Per-track cycle-CONTINUITY score. A step spanning k acquisition cycles (k=1 =
  # integrated, k≥2 = k−1 dark cycles) is weighted w(k) = 2^(−(k−1)/H), H = the
  # dark-cycle "half-life" (tolerance): H small = harsh (only integrated counts),
  # H large = lenient (a few dark cycles barely penalised, so 1–2-cycle steps are
  # forgiven while 10–20-cycle steps are heavily down-weighted). The per-track score
  # is the mean weight over its steps, in [0,1] (1 = every step integrated). Returns
  # a named vector (tid → score). Half-life 0/NA → binary integrated fraction.
  .continuity_scores <- function(df, t_c, half_life = 3) {
    if (is.null(df) || !nrow(df) || !all(c("tid","tim") %in% names(df)) ||
        !is.finite(t_c) || t_c <= 0) return(NULL)
    utid <- unique(df$tid); tf <- match(df$tid, utid)     # integer track codes
    ord  <- order(tf, df$tim, method = "radix")           # radix sort (fast on integers)
    tf   <- tf[ord]; tim <- df$tim[ord]; n <- length(tf)
    step <- c(FALSE, tf[-1] == tf[-n])
    dt   <- c(NA_real_, diff(tim)); dt[!step] <- NA_real_
    k    <- pmax(1L, as.integer(round(dt / t_c)))
    w    <- if (is.finite(half_life) && half_life > 0) 2^(-(k - 1) / half_life)
            else as.numeric(k == 1L)     # half-life 0 → pure integrated fraction
    sc   <- tapply(w[step], tf[step], mean)               # by integer code
    stats::setNames(as.numeric(sc), utid[as.integer(names(sc))])   # map codes → tids
  }

  # ── Data-quality calibration: σ_loc (CVE, no MSD) + integrated/dark score ────
  # Covariance-based estimator (Vestergaard, Blainey & Flyvbjerg 2014): consecutive
  # per-axis steps share one localization's error and are anti-correlated by −σ².
  # With MINFLUX's near-instantaneous dwell (motion-blur R≈0):
  #     σ_loc² = −⟨Δx_i · Δx_{i+1}⟩   (per axis; averaged over x, y[, z]).
  # Uses only INTEGRATED (single-cycle) consecutive step pairs so dark-cycle Δt
  # doesn't bias it. Independent of D — no diffusion fit required. Also returns the
  # dataset integrated-step fraction (the same score the track-quality gate uses).
  trk_cve_bundle <- eventReactive(input$trk_cve_run, {
    df <- tryCatch(trk_gap_split(), error = function(e) NULL)
    req(!is.null(df), nrow(df) > 0, all(c("x","y","tid","tim") %in% names(df)))
    t_c   <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    has_z <- "z" %in% names(df) && sum(is.finite(df$z)) > 0
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

    # Same integrated fraction, but on the FULLY-FILTERED data (post quality-gate /
    # end-trim / length / duration / spatial) — the data the analyses actually use.
    # Lets the panel show raw → after-filters and whether the filters helped.
    frac_integ_of <- function(d) {
      if (is.null(d) || !nrow(d) || !all(c("tid","tim") %in% names(d))) return(NA_real_)
      d <- d[order(d$tid, d$tim), ]; m <- nrow(d)
      sm  <- c(FALSE, d$tid[-1] == d$tid[-m])
      dtt <- c(NA_real_, diff(d$tim)); dtt[!sm] <- NA_real_
      nc  <- if (is.finite(t_c) && t_c > 0) round(dtt / t_c) else rep(1L, m)
      isi <- sm & is.finite(nc) & (nc == 1L); isi[is.na(isi)] <- FALSE
      ns  <- sum(sm, na.rm = TRUE)
      if (ns > 0) sum(isi) / ns else NA_real_
    }
    frac_integ_filt <- frac_integ_of(tryCatch(trk_spatial_filtered(), error = function(e) NULL))
    # RAW reference = TRUE frozen upload data (never filtered / rechunked), NOT df
    # (= trk_gap_split, which reads the current — possibly filtered — manifest). This
    # keeps "raw" at the uploaded value (e.g. 52%) so the raw→analyzed jump is real.
    frac_integ_raw  <- frac_integ_of(tryCatch(trk_raw_original(),      error = function(e) NULL))

    # Dataset-level cycle-CONTINUITY (tolerance-weighted; same weighting as the
    # continuity gate) on TRUE RAW (frozen upload) vs ANALYZED (filtered). This reflects
    # the dark-cycle tolerance, so a dataset with many gentle 1–2-cycle gaps scores
    # higher than its bare integrated fraction.
    cont_of <- function(d, H) {
      if (is.null(d) || !nrow(d) || !all(c("tid","tim") %in% names(d))) return(NA_real_)
      d <- d[order(d$tid, d$tim), ]; m <- nrow(d)
      sm  <- c(FALSE, d$tid[-1] == d$tid[-m])
      dtt <- c(NA_real_, diff(d$tim)); dtt[!sm] <- NA_real_
      kk  <- pmax(1L, as.integer(round(dtt / t_c)))
      ww  <- if (is.finite(H) && H > 0) 2^(-(kk - 1) / H) else as.numeric(kk == 1L)
      if (sum(sm, na.rm = TRUE) > 0) mean(ww[sm], na.rm = TRUE) else NA_real_
    }
    H_disp    <- suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3))
    cont_raw  <- cont_of(tryCatch(trk_raw_original(), error = function(e) NULL), H_disp)  # TRUE raw (frozen at upload)
    cont_filt <- cont_of(tryCatch(trk_spatial_filtered(), error = function(e) NULL), H_disp)  # analyzed (live filtered)

    # CVE lag-1 autocovariance over CONSECUTIVE integrated step pairs.
    lead_integ <- c(integ_step[-1], FALSE)
    pair   <- integ_step & lead_integ
    prod_x <- dx * c(dx[-1], NA_real_)
    prod_y <- dy * c(dy[-1], NA_real_)
    vals   <- c(prod_x[pair], prod_y[pair])
    if (has_z) { prod_z <- dz * c(dz[-1], NA_real_); vals <- c(vals, prod_z[pair]) }
    vals <- vals[is.finite(vals)]
    shiny::validate(shiny::need(length(vals) >= 20,
      "Too few consecutive integrated steps for a CVE estimate — relax filters or check the cycle time."))
    cov1   <- mean(vals)
    sigma2 <- -cov1
    sigma  <- sqrt(max(0, sigma2))
    # Free bonus: bias-corrected D (per axis) = ⟨Δ²⟩/(2Δt) + cov1/Δt.
    sq     <- c(dx[integ_step]^2, dy[integ_step]^2); if (has_z) sq <- c(sq, dz[integ_step]^2)
    sq     <- sq[is.finite(sq)]
    dt_int <- mean(dt[integ_step], na.rm = TRUE)
    D_cve  <- if (length(sq) && is.finite(dt_int) && dt_int > 0) mean(sq)/(2*dt_int) + cov1/dt_int else NA_real_

    # ── Per-track σ_loc heterogeneity ──────────────────────────────────────────
    # One CVE σ per track (tracks with ≥ n_min integrated steps), skipping the
    # first `skip` localizations — MINFLUX iterates onto the molecule, so the
    # opening localizations carry lock-on transient and read worse. The spread of
    # this distribution says whether a single global σ_loc is a good approximation.
    n_min_pt <- max(10L, as.integer(input$trk_sig_min_integ %||% 50L))
    skip_pt  <- max(0L,  as.integer(input$trk_sig_skip %||% 5L))
    pos  <- sequence(rle(df$tid)$lengths)                     # within-track index (df sorted)
    post <- pos > skip_pt                                     # localization kept (past lock-on)
    step_post  <- same & post & c(FALSE, post[-n])            # step between two kept locs
    integ_post <- step_post & is.finite(n_cyc) & (n_cyc == 1L); integ_post[is.na(integ_post)] <- FALSE
    pair_pt <- integ_post & c(integ_post[-1], FALSE)          # consecutive integrated kept pair
    tid_pt  <- df$tid[pair_pt]
    pv <- c(prod_x[pair_pt], prod_y[pair_pt]); tv <- c(tid_pt, tid_pt)
    if (has_z) { pv <- c(pv, prod_z[pair_pt]); tv <- c(tv, tid_pt) }
    okpv <- is.finite(pv); pv <- pv[okpv]; tv <- tv[okpv]
    pertrack_sigma <- if (length(pv)) {
      s2_by  <- tapply(pv, tv, function(v) -mean(v))
      n_by   <- tapply(integ_post, df$tid, sum)[names(s2_by)]
      efo_by <- if ("efo" %in% names(df)) {
                  kl <- post & is.finite(df$efo)
                  tapply(df$efo[kl], df$tid[kl], stats::median)[names(s2_by)]
                } else stats::setNames(rep(NA_real_, length(s2_by)), names(s2_by))
      data.frame(tid = names(s2_by), sigma = sqrt(pmax(0, as.numeric(s2_by))),
                 n_integ = as.integer(n_by), efo_med = as.numeric(efo_by),
                 degenerate = as.numeric(s2_by) <= 0, stringsAsFactors = FALSE)
    } else NULL
    if (!is.null(pertrack_sigma))
      pertrack_sigma <- pertrack_sigma[is.finite(pertrack_sigma$n_integ) &
                                       pertrack_sigma$n_integ >= n_min_pt, , drop = FALSE]
    sig_pt_ok  <- if (!is.null(pertrack_sigma) && nrow(pertrack_sigma))
                    pertrack_sigma$sigma[!pertrack_sigma$degenerate] else numeric(0)
    sig_pt_med <- if (length(sig_pt_ok) >= 5L) stats::median(sig_pt_ok) else NA_real_
    sig_pt_iqr <- if (length(sig_pt_ok) >= 5L)
                    as.numeric(diff(stats::quantile(sig_pt_ok, c(.25, .75)))) else NA_real_
    sig_pt_cv  <- if (is.finite(sig_pt_med) && sig_pt_med > 0 && length(sig_pt_ok) >= 5L)
                    stats::sd(sig_pt_ok) / mean(sig_pt_ok) else NA_real_
    # Precision should improve with photons: negative σ–EFO correlation validates it.
    sig_efo_rho <- if (!is.null(pertrack_sigma) &&
                       sum(is.finite(pertrack_sigma$sigma) & is.finite(pertrack_sigma$efo_med)) >= 8L)
      suppressWarnings(stats::cor(pertrack_sigma$sigma, pertrack_sigma$efo_med,
                                  method = "spearman", use = "complete.obs")) else NA_real_
    # 95% CI on the per-track-median σ_loc (bootstrap over tracks).
    sig_ci <- if (length(sig_pt_ok) >= 5L) .boot_ci(sig_pt_ok, stats::median, B = 1000L) else c(lo = NA_real_, hi = NA_real_)

    # ── Linking / sanity QC ─────────────────────────────────────────────────────
    # (1) Implausible jumps: a step whose |Δr| far exceeds what free diffusion +
    #     localization noise allow for its OWN gap (E[|Δr|²] = 2d·D·Δt + 2d·σ²).
    #     A big jump after one cycle is suspicious; after a long dark gap it isn't
    #     (the Δt term grows), so this naturally accounts for dark cycles. Flags
    #     mis-links / re-acquisition errors. (2) Suspicious tracks: per-track σ_loc
    #     that is degenerate (σ²<0) or a far outlier — candidate bad/mis-linked tracks.
    sig_use <- if (is.finite(sig_pt_med)) sig_pt_med else sigma
    d_dim   <- if (has_z) 3 else 2
    n_implaus <- NA_integer_; frac_implaus <- NA_real_
    if (is.finite(D_cve) && D_cve > 0 && is.finite(sig_use)) {
      dr2  <- dx^2 + dy^2 + if (has_z) dz^2 else 0            # NA where !same
      exp2 <- 2 * d_dim * D_cve * dt + 2 * d_dim * sig_use^2  # expected E[|Δr|²] per step
      ratio2 <- dr2 / exp2
      implaus <- same & is.finite(ratio2) & (ratio2 > 25)     # |Δr| > 5× expected RMS
      n_implaus    <- sum(implaus, na.rm = TRUE)
      frac_implaus <- if (n_steps > 0) n_implaus / n_steps else NA_real_
    }
    n_suspect <- if (!is.null(pertrack_sigma) && nrow(pertrack_sigma) >= 5L) {
      sg  <- pertrack_sigma$sigma
      ok  <- !pertrack_sigma$degenerate & is.finite(sg)
      med0 <- stats::median(sg[ok]); iqr0 <- stats::IQR(sg[ok])
      sum(pertrack_sigma$degenerate | (is.finite(sg) & sg > med0 + 3 * iqr0), na.rm = TRUE)
    } else NA_integer_
    n_tracks_pt <- if (!is.null(pertrack_sigma)) nrow(pertrack_sigma) else NA_integer_

    list(sigma = sigma, sigma2 = sigma2, cov1 = cov1, n_pairs = length(vals),
         D_cve = D_cve, dt_int = dt_int, t_cycle = t_c, has_z = has_z,
         frac_integ = frac_integ_raw, frac_integ_filt = frac_integ_filt,
         cont_raw = cont_raw, cont_filt = cont_filt, cont_H = H_disp,
         n_steps = n_steps, degenerate = sigma2 <= 0,
         pertrack_sigma = pertrack_sigma, sig_pt_med = sig_pt_med,
         sig_pt_iqr = sig_pt_iqr, sig_pt_cv = sig_pt_cv, sig_pt_lo = sig_ci[["lo"]],
         sig_pt_hi = sig_ci[["hi"]], n_pt = length(sig_pt_ok), sig_efo_rho = sig_efo_rho,
         n_implaus = n_implaus, frac_implaus = frac_implaus,
         n_suspect = n_suspect, n_tracks_pt = n_tracks_pt)
  })

  # Feed the CVE σ_loc into the app-wide value: prefer the per-track median (robust
  # to a few bad tracks) when available, else the pooled estimate (skip if
  # degenerate = drift-dominated).
  observeEvent(trk_cve_bundle(), {
    b <- trk_cve_bundle(); if (is.null(b)) return()
    val <- if (is.finite(b$sig_pt_med) && b$sig_pt_med > 0) b$sig_pt_med
           else if (is.finite(b$sigma) && b$sigma > 0 && !isTRUE(b$degenerate)) b$sigma
           else NA_real_
    if (is.finite(val)) trk_sigma_loc_cve(val)
  }, ignoreInit = TRUE)

  output$trk_cve_status <- renderUI({
    b <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    if (is.null(b))
      return(htmltools::div(style = "color:#888; font-size:0.9em;",
        "Press ", htmltools::tags$b("Estimate σ_loc & quality"), " to calibrate from the data (no MSD needed)."))
    # Verdict on the ANALYZED integrated fraction (SAME metric as TSA / cycle-time
    # status, so all quality verdicts agree). The tolerance-weighted continuity is
    # reported alongside (it drives the continuity GATE). "Raw" = the TRUE frozen
    # upload data (never filtered); "Analyzed" = the data the analyses use.
    fi_raw  <- b$frac_integ %||% NA_real_
    fi_filt <- b$frac_integ_filt %||% NA_real_
    c_raw   <- b$cont_raw  %||% NA_real_
    c_filt  <- b$cont_filt %||% NA_real_
    q_use   <- if (is.finite(fi_filt)) fi_filt else fi_raw
    q_verdict <- if (!is.finite(q_use)) c("UNKNOWN", "#888")
                 else if (q_use >= 0.8) c("GOOD", "#1a8a1a")
                 else if (q_use >= 0.5) c("MODERATE", "#b26a00")
                 else c("POOR", "#c03030")
    improved <- is.finite(c_raw) && is.finite(c_filt) && (c_filt - c_raw) > 0.01
    frac_line <- {
      integ_txt <- if (is.finite(fi_raw) && is.finite(fi_filt))
                     sprintf("%.0f%%→%.0f%%", 100 * fi_raw, 100 * fi_filt)
                   else if (is.finite(fi_raw)) sprintf("%.0f%%", 100 * fi_raw) else "—"
      cont_txt  <- if (is.finite(c_raw) && is.finite(c_filt))
                     sprintf("%.2f→%.2f", c_raw, c_filt)
                   else if (is.finite(c_raw)) sprintf("%.2f", c_raw) else "—"
      sprintf("  ·  continuity (raw→analyzed) %s%s  ·  integrated %s",
              cont_txt, if (improved) " ▲" else "", integ_txt)
    }
    sigma_line <- if (isTRUE(b$degenerate) && !is.finite(b$sig_pt_med))
      htmltools::span(style = "color:#c03030;",
        "σ_loc: not resolvable — consecutive steps are positively correlated (drift/directed motion dominates). ",
        "Enter σ_loc manually or use a less directed subset.")
    else if (is.finite(b$sig_pt_med))
      htmltools::span(
        htmltools::tags$b(sprintf("σ_loc = %.2f nm", b$sig_pt_med)),
        if (is.finite(b$sig_pt_lo)) sprintf(" (95%% CI %.2f–%.2f)", b$sig_pt_lo, b$sig_pt_hi) else "",
        sprintf(" · per-track median, n=%d; IQR %.2f nm; CV %.0f%%; pooled %.2f nm",
                b$n_pt, b$sig_pt_iqr %||% NA_real_, 100 * (b$sig_pt_cv %||% NA_real_), b$sigma),
        if (is.finite(b$sig_efo_rho))
          htmltools::span(style = "color:#888;",
            sprintf("  ·  σ–EFO ρ = %.2f%s", b$sig_efo_rho,
                    if (b$sig_efo_rho < -0.2) " ✓ precision improves with photons" else "")) else NULL)
    else
      htmltools::span(htmltools::tags$b(sprintf("σ_loc = %.2f nm", b$sigma)),
        sprintf("  (CVE pooled, %s consecutive integrated pairs)", format(b$n_pairs, big.mark = ",")))
    htmltools::div(style = "font-size:0.92em; margin:4px 0;",
      sigma_line, htmltools::br(),
      htmltools::span(htmltools::tags$b("Data quality (analyzed): "),
        htmltools::span(style = sprintf("color:%s; font-weight:bold;", q_verdict[2]), q_verdict[1]),
        frac_line),
      htmltools::br(),
      htmltools::span(style = "color:#888; font-size:0.94em;",
        sprintf("cycle time = %s  ·  free D (CVE) = %s µm²/s  ·  %s steps",
                if (is.finite(b$t_cycle)) sprintf("%.4g ms", b$t_cycle) else "unknown",
                if (is.finite(b$D_cve))   sprintf("%.3g", b$D_cve / 1e3) else "—",
                format(b$n_steps, big.mark = ","))),
      if (is.finite(b$frac_implaus %||% NA_real_) || is.finite(b$n_suspect %||% NA_real_)) {
        bad <- (is.finite(b$frac_implaus) && b$frac_implaus > 0.01) ||
               (is.finite(b$n_suspect) && is.finite(b$n_tracks_pt) && b$n_tracks_pt > 0 &&
                b$n_suspect / b$n_tracks_pt > 0.05)
        htmltools::div(style = sprintf("font-size:0.9em; margin-top:3px; color:%s;", if (bad) "#b26a00" else "#888"),
          htmltools::tags$b("Linking QC: "),
          if (is.finite(b$n_implaus))
            sprintf("%s implausible jumps (%.2f%% of steps, |Δr|>5× expected for their gap)",
                    format(b$n_implaus, big.mark = ","), 100 * b$frac_implaus) else "—",
          if (is.finite(b$n_suspect))
            sprintf("  ·  %d/%d suspicious tracks (degenerate/outlier σ_loc)", b$n_suspect, b$n_tracks_pt) else "",
          if (bad) "  — inspect at high density / long dark gaps (possible mis-links)." else " — within expectation.")
      } else NULL,
      if (isTRUE(b$degenerate)) NULL else htmltools::div(style = "color:#2a6; font-size:0.9em; margin-top:3px;",
        "→ σ_loc is now the app-wide noise floor for the D-map, JDD, confinement, HMM, directionality and speed. ",
        "Cycle time drives gap-splitting, the quality gate and the D-map.")
    )
  })

  # Per-track σ_loc heterogeneity histogram (σ_loc is NOT a single instrument
  # constant — this shows the spread; a tight peak justifies one global value).
  output$trk_sigma_pertrack_plot <- renderPlotly({
    b  <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    ps <- if (!is.null(b)) b$pertrack_sigma else NULL
    shiny::validate(shiny::need(!is.null(ps) && nrow(ps) >= 5,
      "Run the calibration on data with enough long, integrated tracks to see the per-track σ_loc spread."))
    s  <- ps$sigma[is.finite(ps$sigma) & !ps$degenerate]
    shiny::validate(shiny::need(length(s) >= 5, "Too few resolvable per-track σ_loc values."))
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    med <- stats::median(s)
    plot_ly() %>%
      add_histogram(x = s, nbinsx = min(40L, max(10L, round(length(s)/4))),
                    marker = list(color = "rgba(120,160,220,0.7)",
                                  line = list(color = "rgba(80,120,180,0.9)", width = 0.5)),
                    name = "per-track σ_loc") %>%
      layout(
        xaxis = list(title = "per-track σ_loc (nm)", color = fg, zeroline = FALSE),
        yaxis = list(title = "tracks", color = fg),
        shapes = list(list(type = "line", x0 = med, x1 = med, y0 = 0, y1 = 1, yref = "paper",
                           line = list(color = "#c03030", width = 2, dash = "dash"))),
        annotations = list(list(x = med, y = 1, yref = "paper", yanchor = "bottom",
                                text = sprintf("median %.2f nm", med), showarrow = FALSE,
                                font = list(color = "#c03030", size = 11))),
        paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
        showlegend = FALSE, bargap = 0.04, uirevision = "sigma_pt")
  })

  # Transparency funnel: how many localizations / tracks survive each pipeline
  # stage. Analyzing well-tracked molecules is a legitimate, DISCLOSED choice —
  # this box is that disclosure. Reads debounced stages so it doesn't thrash.
  output$trk_quality_funnel <- renderUI({
    cnt <- function(d) if (is.null(d) || !nrow(d) || !"tid" %in% names(d)) c(0, 0)
                       else c(nrow(d), length(unique(d$tid)))
    rs  <- trk_raw_stats_rv()
    raw <- if (!is.null(rs) && length(rs$tid_length)) c(sum(rs$tid_length), length(rs$tid_length)) else c(NA, NA)
    stg <- list(
      list("Gap-split",         cnt(tryCatch(trk_gap_split_d(),          error = function(e) NULL)), isTRUE(input$trk_split_gaps_enable)),
      list("End-trim",          cnt(tryCatch(trk_end_trimmed_d(),        error = function(e) NULL)), isTRUE(input$trk_endtrim_enable)),
      list("Quality gate",      cnt(tryCatch(trk_quality_filtered_d(),   error = function(e) NULL)), isTRUE(input$trk_quality_enable)),
      list("Metadata filter",   cnt(tryCatch(trk_meta_filtered_d(),      error = function(e) NULL)), TRUE),
      list("Length filter",     cnt(tryCatch(trk_len_filtered_d(),       error = function(e) NULL)), TRUE),
      list("Duration filter",   cnt(tryCatch(trk_dur_filtered_d(),       error = function(e) NULL)), TRUE),
      list("Analyzed (spatial)",cnt(tryCatch(trk_spatial_filtered(),     error = function(e) NULL)), TRUE)
    )
    base <- if (is.finite(raw[1]) && raw[1] > 0) raw[1] else stg[[1]][[2]][1]
    fmt  <- function(v) if (is.finite(v)) format(round(v), big.mark = ",") else "—"
    row  <- function(label, loc, tid, active, pct_base) {
      col <- if (isTRUE(active)) "#333" else "#aaa"
      htmltools::tags$tr(
        htmltools::tags$td(label, style = sprintf("padding:1px 10px 1px 0; color:%s;", col)),
        htmltools::tags$td(fmt(loc), style = sprintf("padding:1px 10px; text-align:right; color:%s;", col)),
        htmltools::tags$td(if (is.finite(pct_base) && pct_base > 0 && is.finite(loc))
                             sprintf("%.0f%%", 100 * loc / pct_base) else "",
                           style = sprintf("padding:1px 10px; text-align:right; color:%s;", col)),
        htmltools::tags$td(fmt(tid), style = sprintf("padding:1px 10px; text-align:right; color:%s;", col)),
        htmltools::tags$td(if (!isTRUE(active)) "(off)" else "",
                           style = "padding:1px 6px; color:#aaa; font-size:0.9em;"))
    }
    rows <- list(if (is.finite(raw[1])) row("Uploaded (raw)", raw[1], raw[2], TRUE, base) else NULL)
    for (s in stg) rows <- c(rows, list(row(s[[1]], s[[2]][1], s[[2]][2], s[[3]], base)))
    htmltools::div(
      htmltools::tags$table(style = "font-size:0.83em; border-collapse:collapse;",
        htmltools::tags$thead(htmltools::tags$tr(lapply(
          c("Stage", "localizations", "% kept", "tracks", ""),
          function(h) htmltools::tags$th(h, style = "padding:1px 10px; text-align:right; border-bottom:1px solid #ccc; color:#666;")))),
        htmltools::tags$tbody(rows)),
      htmltools::tags$p(style = "font-size:0.8em; color:#888; margin:5px 0 0 0;",
        "Analyses run on the last row. Filtering for well-tracked molecules is a disclosed choice — ",
        "report this funnel in your methods. Track count can rise at Gap-split (segments) then fall."))
  })


  # Non-reactive loop state for sequential mode (safe for later())
  is_seq   <- reactive(identical(input$trk_play_mode, "seq"))
  is_simul <- reactive(identical(input$trk_play_mode, "simul"))

  # Bin width (ms) is the master time-resolution control. The frame count is
  # derived from it (window / bin width) rather than entered separately — the
  # old "Frames per playthrough" input was always back-computed into bin_dt.
  # Used by sequential playback (frames per TID) and its step slider; the
  # simultaneous path builds its own bin edges directly from trk_bin_dt.
  n_frames_eff <- reactive({
    tot <- suppressWarnings(as.numeric(input$trk_play_total))
    bw  <- suppressWarnings(as.numeric(input$trk_bin_dt))
    if (!is.finite(tot) || tot <= 0 || !is.finite(bw) || bw <= 0) return(120L)
    max(2L, as.integer(round(tot / bw)))
  })
  loop <- new.env(parent = emptyenv())
  loop$playing     <- FALSE
  loop$sp          <- NULL   # list of data.frames per tid (x, y, t)
  loop$tids        <- NULL   # character vector of tids in play order
  loop$cur_tid     <- 1L
  loop$cur_frame     <- 1L
  loop$ggen      <- 0L      # generation for global (simultaneous) loop
  loop$gplaying  <- FALSE   # playing flag for global loop
  loop$trail_len   <- 0.5
  loop$trail_nmax  <- 50L
  loop$pt_size     <- 6
  loop$col_map     <- NULL
  loop$markers_idx <- 3L
  loop$trail_start <- 2L
  loop$ms_per_unit <- 0
  loop$speed       <- 1    # "Playback speed (×)" multiplier, captured at play start
  loop$gen         <- 0L   # sequential playback generation
  loop$auto_tuning      <- FALSE  # TRUE while apply_auto_tune() is batching input updates
  loop$last_fire_layout <- TRUE   # TRUE = last trigger was fire_layout (full rebuild)

  # Base animation cadence (~16 fps at 1x speed). Replaces the old "Tick interval"
  # input: it drives the simultaneous (prebinned) loop and floors the sequential
  # loop. The "Playback speed (×)" slider scales the per-frame delay around this.
  TICK_INTERVAL_SEC <- 0.06
  # Hard lower bound on per-frame delay (~50 fps). Stops the loop from outpacing
  # the websocket (which is what made playback "refuse to stop"), no matter how
  # high the speed multiplier or how small "ms per time unit" gets.
  MIN_INTERVAL_SEC  <- 0.02

  # Effective per-frame delay (seconds) for the current speed multiplier.
  play_interval_for <- function(speed) {
    s <- suppressWarnings(as.numeric(speed))
    if (!is.finite(s) || s <= 0) s <- 1
    max(MIN_INTERVAL_SEC, TICK_INTERVAL_SEC / s)
  }
  
  # ---- Optional chunk picker (manifest.csv in same folder as upload) ----
  output$trk_chunk_ui <- renderUI({
    req(input$trk_file)
    
    # If we're still preparing chunks, show a placeholder message
    if (isTRUE(rv$chunking)) {
      return(tags$div(
        "Preparing chunks... This can take a bit for very large files.",
        style = "color:#555;"
      ))
    }
    
    # Otherwise, if a manifest exists, show the chunk selector
    dirp <- if (!is.null(rv$manifest_dir)) rv$manifest_dir else dirname(input$trk_file$datapath[1])
    manf <- file.path(dirp, "manifest.csv")
    if (file.exists(manf)) {
      man <- readr::read_csv(manf, show_col_types = FALSE)
      chunk_choices <- setNames(man$file, paste0("Chunk ", man$chunk, " (", man$n_tids, " tids)"))
      all_choices <- c("All TIDs" = "__all__", chunk_choices)
      return(tagList(
        selectInput("trk_chunk_sel", "Playback chunk", choices = all_choices, selected = "__all__"),
        helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
          "Filtering, gap-splitting, analysis and the full-track overlay always use ALL ",
          "TIDs. This selector only limits which TIDs are ANIMATED during playback, to ",
          "keep the animation smooth — \"All TIDs\" animates everything.")
      ))
    }

    NULL
  })
  output$trk_msd_tid_ui <- renderUI({
    req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_show_single))
    df <- trk(); req(nrow(df) > 0)
    min_pts <- max(5L, as.integer(input$trk_msd_min_pts %||% 20L))
    # Only offer TIDs that have enough points to compute MSD
    counts  <- table(df$tid)
    tids    <- sort(names(counts[counts >= min_pts]))
    req(length(tids) > 0)
    selectizeInput(
      "trk_msd_tid", "Choose TID for overlay:",
      choices = tids, selected = tids[1],
      options = list(maxOptions = 5000)
    )
  })
  # View mode toggle: 2D / 3D (only rendered when data has a z column)
  output$trk_view_mode_ui <- renderUI({
    req(isTRUE(has_z()))
    tagList(
      radioButtons("trk_view_mode", NULL,
                   choices  = c("2D" = "2d", "3D" = "3d"),
                   selected = isolate(input$trk_view_mode %||% "2d"),
                   inline   = TRUE),
      numericInput("trk_z_stretch", "Z display scale",
                   value = isolate(input$trk_z_stretch %||% 1.0),
                   min = 0.2, max = 10, step = 0.1)
    )
  })

  # Switching view mode requires a full plot rebuild
  observeEvent(input$trk_view_mode, {
    # Fire only on a GENUINE 2D/3D change. On upload the input is (re)created with the
    # same value the plot was just rendered with (loop$rendered_view), so its creation
    # fire is skipped — no redundant rebuild, no fragile time window.
    if (identical(input$trk_view_mode %||% "2d", loop$rendered_view %||% "2d")) return()
    fire_layout("view_mode_2d3d")
  }, ignoreInit = TRUE)
  # Z stretch: debounce so rapid typing / arrow-key changes don't hammer fire_layout
  trk_z_stretch_d <- debounce(reactive({ input$trk_z_stretch }), 500)
  observeEvent(trk_z_stretch_d(), {
    val <- trk_z_stretch_d()
    req(is.numeric(val), is.finite(val), val >= 0.1)
    if (isTRUE(all.equal(val, loop$rendered_zs %||% 1.0))) return()  # unchanged since last render → skip creation fire
    fire_layout("z_stretch")
  }, ignoreInit = TRUE)

  # ---- Big-file detection & chunking/compression on upload ---------------------
  observeEvent(input$trk_file, {
    req(input$trk_file)
    rv$bigfile_notice        <- FALSE
    rv$manifest_dir          <- NULL
    rv$original_manifest_dir <- NULL

    # input$trk_file may hold several rows (multiple files selected for merging)
    n_files   <- nrow(input$trk_file)
    paths     <- input$trk_file$datapath
    names_v   <- tolower(input$trk_file$name)
    size_bytes <- sum(suppressWarnings(as.numeric(input$trk_file$size)), na.rm = TRUE)
    if (!is.finite(size_bytes) || size_bytes <= 0) {
      info <- tryCatch(file.info(paths), error = function(e) NULL)
      if (!is.null(info)) size_bytes <- sum(info$size, na.rm = TRUE)
    }
    size_mb <- round(as.numeric(size_bytes) / (1024^2), 1)

    # If a single file is uploaded next to an existing manifest, adopt it and bail.
    # (Skipped for multi-file merges — those always re-chunk the combined data.)
    if (n_files == 1L) {
      dirp <- dirname(paths[1])
      if (file.exists(file.path(dirp, "manifest.csv"))) {
        rv$manifest_dir          <- dirp
        rv$original_manifest_dir <- dirp
        trk_data_key(as.character(as.numeric(Sys.time())))
        rv$chunking <- FALSE
        if (!rv$bigfile_notice) {
          showNotification("Found existing chunk manifest; using chunked loading.", type = "message", duration = 6)
          rv$bigfile_notice <- TRUE
        }
        return(invisible(NULL))
      }
    }

    # Always chunk/compress every upload for efficient chunked loading
    rv$chunking <- TRUE
    on.exit({ rv$chunking <- FALSE }, add = TRUE)
    safe_base <- gsub("[^A-Za-z0-9_-]+", "_", tools::file_path_sans_ext(basename(names_v[1])))
    if (n_files > 1L) safe_base <- paste0(safe_base, "_merged", n_files)
    out_dir <- file.path(tempdir(), paste0("chunks_", safe_base, "_", as.integer(Sys.time())))

    # Loader for one file → data.frame (or NULL).  Returns attr "json_itr".
    load_one_trk <- function(p, nm) {
      d <- NULL; jitr <- NA_integer_
      if (grepl("\\.csv(\\.gz)?$", nm)) {
        d <- data.table::fread(p, showProgress = FALSE, nThread = max(1, parallel::detectCores()-1)) |>
             as.data.frame()
        d <- canonicalize_minflux_df(d)
      } else if (grepl("\\.qs$", nm)) {
        d <- qs::qread(p)
      } else if (grepl("\\.fst$", nm)) {
        d <- fst::read_fst(p) |> as.data.frame()
      } else if (grepl("\\.rds$", nm)) {
        d <- readRDS(p)
        if (!is.data.frame(d)) return(NULL)
        d <- canonicalize_minflux_df(d)
      } else if (grepl("\\.json$", nm)) {
        d <- read_minflux_json_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
        jitr <- if ("itr" %in% names(d) && nrow(d) > 0) max(d$itr, na.rm = TRUE) else NA_integer_
      } else if (grepl("\\.npy$", nm)) {
        d <- read_minflux_npy_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
        jitr <- if ("itr" %in% names(d) && nrow(d) > 0) max(d$itr, na.rm = TRUE) else NA_integer_
      } else {
        return(NULL)
      }
      attr(d, "json_itr") <- jitr
      d
    }

    # Load all files; prefix TIDs per file so IDs never collide on merge.
    df       <- NULL
    json_itr <- NA_integer_
    withProgress(message = sprintf("Loading %d file(s) (%.1f MB)...", n_files, size_mb), value = 0, {
      dfs <- vector("list", n_files)
      for (fi in seq_len(n_files)) {
        incProgress(1 / n_files, detail = basename(names_v[fi]))
        di <- tryCatch(load_one_trk(paths[fi], names_v[fi]), error = function(e) NULL)
        if (is.null(di) || !is.data.frame(di) || !nrow(di)) next
        ji <- attr(di, "json_itr"); if (!is.na(ji)) json_itr <- ji
        if (n_files > 1L && "tid" %in% names(di))
          di$tid <- paste0("F", fi, "_", as.character(di$tid))
        dfs[[fi]] <- di
      }
      dfs <- Filter(Negate(is.null), dfs)
      df  <- if (length(dfs)) dplyr::bind_rows(dfs) else NULL
    })
    if (n_files > 1L && !is.null(df))
      showNotification(sprintf("Merged %d files into one dataset (%s localizations).",
                               n_files, format(nrow(df), big.mark = ",")),
                       type = "message", duration = 7)

    if (!is.null(df) && is.data.frame(df)) {
      if (!all(c("tid","x","y","tim") %in% names(df))) {
        showNotification("File is missing required columns (tid, x, y, tim).", type = "error", duration = 8)
        return(invisible(NULL))
      }
      df$tid <- as.character(df$tid)

      # Keep only the columns tracking needs  -  reduces chunk file size significantly
      keep_cols <- intersect(c("tid","x","y","z","tim","itr","efo","cfr","dcr","eco","dwell","fluo"), names(df))
      df <- df[, keep_cols, drop = FALSE]

      # NOTE: no automatic CFR pre-filter on upload. CFR is not a strong artefact
      # indicator for tracking, so all localizations are kept. The user can still
      # filter CFR manually via the metadata filter (CFR range) if desired.

      # Show JSON iteration notification (mirrors SMLM tab behaviour)
      if (!is.na(json_itr)) {
        showNotification(
          sprintf(".json or .npy import: using final localization (itr = %d) for x/y.", as.integer(json_itr)),
          type = "message", duration = 8
        )
      }

      # Compute per-TID summary (duration in ms, n_pts) and sort longest-first
      # tim * 1000 matches the ms conversion trk_raw_df() applies after chunk loading
      tid_sum <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          duration_ms = (max(tim, na.rm = TRUE) - min(tim, na.rm = TRUE)) * 1000,
          n_pts        = dplyr::n(),
          .groups      = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(duration_ms))

      # Reorder df rows so chunks end up sorted by duration (longest first)
      df <- df[order(match(df$tid, tid_sum$tid)), , drop = FALSE]

      # Freeze raw stats for histogram grey bars and filter reset.
      # Done once per upload — never overwritten by chunk switching or rechunking.
      # tim stored in ms (raw file seconds; trk_raw_df multiplies by 1000).
      loc_vec <- function(col) if (col %in% names(df)) df[[col]][is.finite(df[[col]])] else numeric(0)
      trk_raw_stats_rv(list(
        tid_length   = as.numeric(table(df$tid)),
        tid_duration = tid_sum$duration_ms,
        efo          = loc_vec("efo"),
        cfr          = loc_vec("cfr"),
        dcr          = loc_vec("dcr"),
        tim          = loc_vec("tim") * 1000,
        eco          = loc_vec("eco"),
        dwell        = loc_vec("dwell"),
        fluo         = loc_vec("fluo")
      ))
      trk_rechunk_stats_rv(NULL)  # clear any leftover rechunk stats from prior dataset
      trk_rechunk_bbox(NULL)      # clear locked spatial extent from prior dataset

      # use user-entered chunk value if provided, else the default
      user_tpc <- suppressWarnings(as.integer(isolate(input$trk_tids_per_chunk)))
      if (!is.finite(user_tpc) || user_tpc < 50) user_tpc <- TIDS_PER_CHUNK_DEFAULT

      manf <- split_by_tid_to_dir(
        df, out_dir,
        tids_per_chunk = user_tpc,
        fmt = CHUNK_FORMAT
      )
      if (file.exists(manf)) {
        rv$manifest_dir          <- dirname(manf)
        rv$original_manifest_dir <- rv$manifest_dir
        trk_data_key(as.character(as.numeric(Sys.time())))
        # Save per-TID summary alongside the manifest
        readr::write_csv(tid_sum, file.path(rv$manifest_dir, "tid_summary.csv"))
        man <- readr::read_csv(manf, show_col_types = FALSE)
        rv$chunking <- FALSE
        avg_dur <- round(mean(tid_sum$duration_ms, na.rm = TRUE))
        showNotification(
          paste0("Split into ", nrow(man), " chunks @ ~", user_tpc,
                 " TIDs/chunk. Avg trace lifetime: ", avg_dur, " ms."),
          type = "message", duration = 8
        )
      } else {
        # Fallback: compress to a single .qs chunk
        out_qs <- file.path(out_dir, paste0(safe_base, ".qs"))
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        qs::qsave(df, out_qs, preset = "high")
        man <- data.frame(chunk = 1L, n_tids = length(unique(df$tid)), file = out_qs)
        readr::write_csv(man, file.path(out_dir, "manifest.csv"))
        readr::write_csv(tid_sum, file.path(out_dir, "tid_summary.csv"))
        rv$manifest_dir          <- out_dir
        rv$original_manifest_dir <- out_dir
        trk_data_key(as.character(as.numeric(Sys.time())))
        rv$chunking <- FALSE
        showNotification("File compressed to a single .qs chunk.", type = "message", duration = 8)
      }
    } else {
      showNotification("File type not supported for chunking; loading as-is.", type = "warning", duration = 6)
    }
  }, ignoreInit = TRUE)
  # Auto-plot fallback: fires only if chunking failed (unsupported file type).
  # For all supported types, the chunk_sel observer handles the first plot after
  # chunk_ui renders.
  observeEvent(input$trk_file, {
    if (is.null(rv$manifest_dir)) {
      bump_clean()
      fire_layout("file_nomanifest")
    }
  }, ignoreInit = TRUE)
  
  # ---- Load CSV/fast formats OR a pre-split chunk ----
  trk_raw_df <- reactive({
    req(input$trk_file)
    required <- c("tid","x","y","tim")
    
    # Source for the WHOLE filter pipeline. The filters and gap-splitting must
    # operate on ALL TIDs, never a single chunk, so we always load every chunk of
    # the CURRENT manifest (rv$manifest_dir) — the chunk selector does not restrict
    # the data. After Filter & Rechunk the current manifest is the rechunked set;
    # use Reset Ranges (which restores rv$original_manifest_dir) to return to the
    # full original, un-split data before re-splitting.
    df <- NULL
    if (!is.null(rv$manifest_dir)) {
      manf <- file.path(rv$manifest_dir, "manifest.csv")
      if (file.exists(manf)) {
        man <- tryCatch(readr::read_csv(manf, show_col_types = FALSE), error = function(e) NULL)
        if (!is.null(man) && nrow(man) > 0 && "file" %in% names(man)) {
          read_chunk <- function(f) tryCatch(
            if (grepl("\\.qs$",  f, TRUE)) qs::qread(f)
            else if (grepl("\\.fst$", f, TRUE)) fst::read_fst(f) |> as.data.frame()
            else NULL,
            error = function(e) NULL
          )
          dfs <- lapply(man$file, read_chunk)
          df  <- dplyr::bind_rows(Filter(Negate(is.null), dfs))
        }
      }
    }
    if (!is.null(df)) {
      # df already loaded (all chunks) — skip to post-processing below
    } else {
      path <- input$trk_file$datapath[1]
      name <- tolower(input$trk_file$name[1])
      df <- if (grepl("\\.qs$", name)) {
        qs::qread(path)
      } else if (grepl("\\.fst$", name)) {
        fst::read_fst(path) |> as.data.frame()
      } else if (grepl("\\.rds$", name)) {
        obj <- readRDS(path); shiny::validate(shiny::need(is.data.frame(obj), "RDS must contain a data.frame")); obj
      } else if (grepl("\\.csv(\\.gz)?$", name)) {
        data.table::fread(path, showProgress = FALSE, nThread = max(1, parallel::detectCores()-1)) |> as.data.frame()
      } else if (grepl("\\.json$", name)) {
        withProgress(message = "Reading MINFLUX JSON...", value = 0, {
          out <- read_minflux_json_as_df(path, keep_final_only = TRUE, keep_valid_only = TRUE)
          if (is.data.frame(out) && nrow(out) > 0) {
            itr_used <- if ("itr" %in% names(out)) max(out$itr, na.rm = TRUE) else NA_integer_
            if (!is.na(itr_used)) {
              showNotification(
                sprintf("JSON import: using final localization (itr = %d) for x/y.", as.integer(itr_used)),
                type = "message", duration = 8
              )
            }
          }
          out
        })
      } else if (grepl("\\.npy$", name)) {
        withProgress(message = "Reading MINFLUX .npy...", value = 0, {
          read_minflux_npy_as_df(path, keep_final_only = TRUE, keep_valid_only = TRUE)
        })
      } else {
        shiny::validate(shiny::need(FALSE, "Unsupported file type"))
      }
    }
    
    # --- common post-processing (applies to both chunk & single-file) ---
    miss <- setdiff(required, names(df))
    shiny::validate(shiny::need(length(miss) == 0, paste("Missing columns:", paste(miss, collapse = ", "))))
    
    df$tid <- as.character(df$tid)
    df$x   <- as.numeric(df$x)
    df$y   <- as.numeric(df$y)
    df$tim <- as.numeric(df$tim)
    
    # tim is always in seconds  -  convert to milliseconds
    df$tim <- df$tim * 1000

    df
  })

  # TRUE, frozen raw — the ORIGINAL uploaded data, NEVER touched by the filters or by
  # Filter & Rechunk. trk_raw_df() reads the CURRENT manifest (which becomes the
  # FILTERED set after a rechunk), so it can't be used as the "raw" reference. This
  # reactive always reads rv$original_manifest_dir (frozen at upload, restored by
  # Reset Ranges), so "raw" stays the uploaded data (e.g. 52% integrated) regardless
  # of any filtering/rechunk. Used for the RAW side of the σ_loc estimator and the
  # TSA before/after. Only recomputes on a genuinely new upload.
  trk_raw_original <- reactive({
    dir <- rv$original_manifest_dir
    # Fallback: unchunked / manifest-less loads → best-effort current raw
    if (is.null(dir) || !file.exists(file.path(dir, "manifest.csv")))
      return(tryCatch(trk_raw_df(), error = function(e) NULL))
    man <- tryCatch(readr::read_csv(file.path(dir, "manifest.csv"), show_col_types = FALSE),
                    error = function(e) NULL)
    if (is.null(man) || !nrow(man) || !"file" %in% names(man)) return(NULL)
    read_chunk <- function(f) tryCatch(
      if (grepl("\\.qs$",  f, TRUE)) qs::qread(f)
      else if (grepl("\\.fst$", f, TRUE)) fst::read_fst(f) |> as.data.frame()
      else NULL, error = function(e) NULL)
    df <- dplyr::bind_rows(Filter(Negate(is.null), lapply(man$file, read_chunk)))
    if (is.null(df) || !nrow(df) || !all(c("tid","tim") %in% names(df))) return(NULL)
    df$tid <- as.character(df$tid)
    for (cc in intersect(c("x","y","tim"), names(df))) df[[cc]] <- as.numeric(df[[cc]])
    df$tim <- df$tim * 1000                  # s → ms (match trk_raw_df post-processing)
    df
  })


  # ---- Dynamic range UIs ----
  # Bounds come from trk_raw_stats_rv() (frozen at upload, whole dataset) so
  # Sticky bar label: filename + current chunk indicator
  output$trk_bar_label <- renderUI({
    req(input$trk_file)
    nm <- if (nrow(input$trk_file) > 1L)
            sprintf("%d files merged (%s, …)", nrow(input$trk_file), input$trk_file$name[1])
          else input$trk_file$name[1]
    chunk_info <- {
      sel <- input$trk_chunk_sel
      if (identical(sel, "__all__"))        " — All chunks"
      else if (isTruthy(sel))              paste0(" — ", basename(sel))
      else                                 ""
    }
    tags$span(style = "font-size:13px; color:#555;", nm, chunk_info)
  })

  # sliders reflect the full data extent regardless of which chunk is selected.
  # Re-renders on: upload (trk_raw_stats_rv() changes) OR Reset Ranges press
  # (reset_ranges_trigger increments).  Never re-renders on chunk switches, so
  # user-set filter values survive chunk navigation.
  make_range_ui <- function(col, label) {
    renderUI({
      reset_ranges_trigger()        # re-render on Reset Ranges press
      rs_orig    <- trk_raw_stats_rv()
      rs_rechunk <- trk_rechunk_stats_rv()
      # Bounds from original upload so the slider never permanently shrinks
      orig_v <- if (!is.null(rs_orig) && col %in% names(rs_orig)) rs_orig[[col]] else {
        df <- trk_raw_df(); if (!col %in% names(df)) return(NULL); df[[col]]
      }
      orig_v <- orig_v[is.finite(orig_v)]; if (!length(orig_v)) return(NULL)
      bounds <- range(orig_v, na.rm = TRUE)
      # Value: rechunk stats when available (all remaining TIDs), else full original range
      val <- if (!is.null(rs_rechunk) && col %in% names(rs_rechunk)) {
        curr_v <- rs_rechunk[[col]][is.finite(rs_rechunk[[col]])]
        if (length(curr_v)) range(curr_v, na.rm = TRUE) else bounds
      } else bounds
      numericRangeInput(paste0(col, "_rng"), label, value = val, min = bounds[1], max = bounds[2])
    })
  }
  output$trk_efo_ui   <- make_range_ui("efo",   "EFO range")
  output$trk_cfr_ui   <- make_range_ui("cfr",   "CFR range")
  output$trk_dcr_ui   <- make_range_ui("dcr",   "DCR range")
  output$trk_tim_ui   <- make_range_ui("tim",   "TIM range (ms)")
  output$trk_eco_ui   <- make_range_ui("eco",   "ECO range")
  output$trk_dwell_ui <- make_range_ui("dwell", "Dwell range")
  output$trk_fluo_ui  <- make_range_ui("fluo",  "FLUO range")

  # Stage 2 (length) auto-detection. The slider BOUNDS reflect the metadata-filtered
  # (post-split) population, but update at PLOT cadence only — keyed on data_ver()
  # (bumps on upload / Plot Data / Reset Ranges), NOT live while other sliders move.
  # So the range always describes the currently PLOTTED data and stops jumping as you
  # drag. The VALUE is PRESERVED (clamped into the new bounds) so a manual minimum
  # (e.g. 50) is not wiped; it snaps to full only on new upload / Reset Ranges /
  # post-rechunk, signalled by trk_len_reset_flag.
  output$trk_tid_len_ui <- renderUI({
    reset_ranges_trigger()
    data_ver()                                                           # update at PLOT cadence (upload / Plot Data), not live
    seg_df <- tryCatch(isolate(trk_meta_filtered()), error = function(e) NULL)  # isolate: no live re-fire while dragging
    lv <- if (!is.null(seg_df) && nrow(seg_df) > 0 && "tid" %in% names(seg_df))
            as.numeric(table(seg_df$tid)) else numeric(0)
    if (!length(lv)) {                          # pipeline not ready → rechunk stats, else upload
      rs_re <- trk_rechunk_stats_rv(); rs_or <- trk_raw_stats_rv()
      lv <- if (!is.null(rs_re) && length(rs_re$tid_length)) rs_re$tid_length
            else if (!is.null(rs_or) && length(rs_or$tid_length)) rs_or$tid_length else numeric(0)
    }
    if (!length(lv)) return(NULL)
    bounds <- range(lv, na.rm = TRUE)
    curr <- isolate(input$trk_tid_len)
    snap <- isTRUE(isolate(trk_len_reset_flag())) || is.null(curr) ||
            length(curr) != 2 || any(!is.finite(as.numeric(curr)))
    val  <- if (snap) bounds
            else c(max(bounds[1], min(as.numeric(curr))), min(bounds[2], max(as.numeric(curr))))
    isolate(trk_len_reset_flag(FALSE))
    numericRangeInput("trk_tid_len", "Trace ID length (localizations)",
                      value = val, min = bounds[1], max = bounds[2])
  })

  # Histogram zoom range — always anchored to the ORIGINAL upload stats so the
  # range never shrinks after Filter & Rechunk.  If it used tid_summary_rv() or
  # trk_raw_df() (both reflect the current, possibly filtered, manifest), the zoom
  # window would clip the grey bars to the filtered range, making them look
  # identical to the blue bars.
  output$trk_hist_range_ui <- renderUI({
    req(input$trk_hist_param)
    param <- input$trk_hist_param
    rs <- trk_raw_stats_rv()
    if (param == "tid_length") {
      vals <- if (!is.null(rs) && length(rs$tid_length)) rs$tid_length else {
        df <- tryCatch(trk_raw_df(), error = function(e) NULL)
        req(!is.null(df)); as.numeric(table(df$tid))
      }
    } else if (param == "tid_duration") {
      vals <- if (!is.null(rs) && length(rs$tid_duration)) rs$tid_duration else {
        ts <- tid_summary_rv()
        if (!is.null(ts) && nrow(ts) > 0) ts$duration_ms else {
          df <- tryCatch(trk_raw_df(), error = function(e) NULL)
          req(!is.null(df))
          df %>% dplyr::group_by(tid) %>%
            dplyr::summarise(d = max(tim) - min(tim), .groups = "drop") %>% dplyr::pull(d)
        }
      }
    } else {
      vals <- if (!is.null(rs) && param %in% names(rs)) rs[[param]] else {
        df <- tryCatch(trk_raw_df(), error = function(e) NULL)
        req(!is.null(df), param %in% names(df)); df[[param]]
      }
      vals <- vals[is.finite(vals)]
    }
    req(length(vals) > 0)
    rng <- range(vals, na.rm = TRUE)
    numericRangeInput("trk_hist_range", "Histogram zoom range:",
                      value = rng, min = rng[1], max = rng[2])
  })

  # Reset all filter inputs to the original upload's full range, and restore
  # rv$manifest_dir to the original upload so trk_raw_df sees all TIDs again.
  # All filter UIs (metadata, tid_len, tid_dur) are renderUI-wrapped — incrementing
  # reset_ranges_trigger forces them all to re-render with value = full rng from
  # trk_raw_stats_rv().  updateNumericRangeInput does not reliably update
  # renderUI-wrapped widgets and is intentionally not used here.
  reset_tracking_filters <- function() {
    if (!is.null(rv$original_manifest_dir))
      rv$manifest_dir <- rv$original_manifest_dir
    trk_rechunk_stats_rv(NULL)  # clear rechunk stats so UIs revert to original upload ranges
    trk_rechunk_bbox(NULL)      # clear locked extent so axes revert to per-TID bbox
    trk_spatial_filter(NULL)    # clear any active spatial crop
    trk_len_reset_flag(TRUE)    # snap length + duration sliders to full range
    reset_dur_filter_flag(TRUE)
    reset_ranges_trigger(isolate(reset_ranges_trigger()) + 1L)
  }

  observeEvent(input$trk_reset_filters, {
    reset_tracking_filters()
  }, ignoreInit = TRUE)


  # When TRUE, Stage 3 (trk_dur_filtered) bypasses input$trk_dur_filter entirely
  # so the stale value from a previous dataset never touches new data. Cleared by
  # trk_dur_filter_ui once it re-renders with the correct new range.
  # Must only be set on NEW UPLOADS (trk_raw_stats_rv changes), NOT on rechunks.
  # trk_data_key fires on both upload and rechunk; trk_raw_stats_rv only on upload.
  reset_dur_filter_flag <- reactiveVal(FALSE)
  # Sibling flag for the length slider: when TRUE the slider snaps to the full
  # detected range; otherwise it PRESERVES the user's current value (clamped to the
  # new bounds). Set on new upload / Reset Ranges / after Filter & Rechunk.
  trk_len_reset_flag    <- reactiveVal(FALSE)
  # Incrementing this forces trk_dur_filter_ui to re-render with value = full rng.
  # Used by reset_tracking_filters() as a reliable alternative to updateNumericRangeInput,
  # which does not reliably update renderUI-wrapped widgets.
  reset_ranges_trigger  <- reactiveVal(0L)

  # On new data upload, snap the TID-level filters to the new full range.
  # trk_raw_stats_rv() is frozen at upload and never updated by rechunking, so
  # this observer fires only on genuine new uploads.
  observeEvent(trk_raw_stats_rv(), {
    rs <- trk_raw_stats_rv()
    if (is.null(rs)) return()
    if (length(rs$tid_length))    trk_len_reset_flag(TRUE)
    if (length(rs$tid_duration))  reset_dur_filter_flag(TRUE)
  }, ignoreInit = TRUE)

  # ---- Three-stage filtering pipeline ----------------------------------------
  # Pipeline: Stage 1 (metadata) ÃÂ¢Ã¢ÂÂ Ã¢ÂÂ Stage 2 (trace ID length) ÃÂ¢Ã¢ÂÂ Ã¢ÂÂ Stage 3 (duration).
  # Each stage uses metrics computed from the PREVIOUS stage's output so that
  # upstream filters propagate correctly (e.g. metadata filtering shrinks trace
  # localization counts, which then drives the length filter).
  # All filter inputs use NULL pass-through (input not yet rendered ÃÂ¢Ã¢ÂÂ Ã¢ÂÂ all data
  # passes) and a stale-range guard (filter range outside data range ÃÂ¢Ã¢ÂÂ Ã¢ÂÂ all data
  # passes).  TID values are coerced to character to guard against type mismatches
  # when readr parses numeric TIDs from tid_summary.csv.

  # Stage 1: per-localization metadata filters (efo, cfr, dcr, tim, eco, dwell, fluo)
  # Stage 4 (metadata filter): now runs AFTER track cleaning (split → end-trim →
  # quality gate), so efo/cfr/dcr/… filters apply to the cleaned tracks. Track
  # cleaning is done on the RAW timing (dark-cycle detection is not distorted by
  # metadata-removed localizations), which is the correct order.
  trk_meta_filtered <- reactive({
    df <- trk_quality_filtered()
    shiny::validate(shiny::need(all(c("x","y","tid","tim") %in% names(df)),
                                "Data must contain x, y, tid, tim."))
    rng_filter <- function(dat, nm) {
      rng <- input[[paste0(nm, "_rng")]]
      if (is.null(rng) || !nm %in% names(dat)) return(dat)
      v <- dat[[nm]]; v <- v[is.finite(v)]
      # Stale guard: if the stored range has no overlap with the data, pass through.
      if (length(v) && (rng[1] > max(v) || rng[2] < min(v))) return(dat)
      dplyr::filter(dat, is.na(.data[[nm]]) | dplyr::between(.data[[nm]], rng[1], rng[2]))
    }
    df %>%
      rng_filter("efo")   %>% rng_filter("cfr")   %>% rng_filter("dcr") %>%
      rng_filter("tim")   %>% rng_filter("eco")   %>% rng_filter("dwell") %>%
      rng_filter("fluo")
  })

  # ===========================================================================
  #  TRACKING FILTER PIPELINE  (read top-to-bottom; each stage feeds the next)
  #    raw → split → end-trim → continuity gate → metadata → length → duration
  #  The cycle-time + TSA machinery lives further down, "CYCLE-TIME & TSA SUPPORT".
  # ===========================================================================

  # Cleaning-apply trigger. Split / end-trim / continuity-gate are HEAVY (per-track
  # passes over ALL localizations), so they must not recompute on every slider
  # touch. They read their parameters via isolate() and recompute only when this
  # bumps — on Plot Data, a new upload, or Filter & Rechunk. (trk_raw_df is still a
  # live dependency, so a new manifest also refreshes them.) Metadata / length /
  # duration filters stay LIVE (they read the cached cleaning output).
  clean_trigger <- reactiveVal(0L)
  bump_clean    <- function() clean_trigger(isolate(clean_trigger()) + 1L)

  # Stage 1 (gap splitting) — the FIRST cleaning stage, on the RAW data, so dark-
  # cycle detection sees the true acquisition timing (metadata filtering, which
  # can remove localizations and widen gaps, now happens LATER). Splits tracks
  # wherever an interval is longer than the threshold (in acquisition cycles). New
  # TIDs "<orig>_s1", "<orig>_s2", … In ADAPTIVE (Vogler/Bela) mode the threshold
  # is data-driven (median+SD+1 of the total-cycle counts, capped by maxOffTime)
  # and segments shorter than the min-length are dropped.
  trk_gap_split <- reactive({
    clean_trigger()                        # apply only on Plot Data / upload / rechunk
    df <- trk_raw_df()
    if (!isTRUE(isolate(input$trk_split_gaps_enable)) || nrow(df) == 0) return(df)
    adaptive   <- isTRUE(isolate(input$trk_split_adaptive))
    gap_cycles <- if (adaptive) isolate(trk_adaptive_gap_cycles())
                  else suppressWarnings(as.numeric(isolate(input$trk_split_gap_cycles)))
    modal_dt   <- isolate(trk_cycle_time())
    if (!is.finite(gap_cycles) || gap_cycles <= 0 ||
        !is.finite(modal_dt) || modal_dt <= 0) return(df)
    gap_thresh <- gap_cycles * modal_dt   # ms
    # Min-length applies in BOTH modes now (was adaptive-only): drop short fragments
    # so a fixed 7-cycle split doesn't manufacture 2–3-loc stubs and flood the plot.
    min_len <- max(1L, as.integer(isolate(input$trk_split_min_len %||% 50L)))
    # Idempotency: strip our own "_s<N>" suffix (only when present) before splitting
    # so it always operates on the ORIGINAL track. After a Filter & Rechunk bakes
    # "_s" segments into the saved TIDs, trk_raw_df reloads them and (split still on)
    # they'd otherwise be split a SECOND time (251_s1 → 251_s1_s1 …). Stripping
    # recombines siblings in time order → splitting twice == splitting once.
    if (any(grepl("_s", df$tid, fixed = TRUE)))
      df$tid <- sub("_s\\d+$", "", as.character(df$tid), perl = TRUE)
    # Fully VECTORISED split (no per-track split()/loop) with an INTEGER-coded RADIX
    # sort. Sorting the character TIDs directly (base order() on character) was the
    # real bottleneck on full raw data — ~25× slower than radix-sorting integer codes.
    tf  <- match(df$tid, unique(df$tid))              # integer code per track
    ord <- order(tf, df$tim, method = "radix")
    df  <- df[ord, , drop = FALSE]; tf <- tf[ord]
    n   <- nrow(df); tid <- as.character(df$tid)
    newtrk <- c(TRUE, tf[-1] != tf[-n])               # first row of each track
    dt     <- c(NA_real_, diff(df$tim)); dt[newtrk] <- NA_real_
    is_gap <- is.finite(dt) & (dt > gap_thresh)       # split point at this row
    grp    <- cumsum(newtrk)                          # track id 1..K
    cg     <- cumsum(is_gap)
    seg    <- cg - cg[newtrk][grp] + 1L               # 1-based segment # within track
    nseg   <- (as.integer(rowsum(as.integer(is_gap), grp)) + 1L)[grp]  # #segments in track
    new_tid <- tid; spl <- nseg > 1L
    if (any(spl)) new_tid[spl] <- paste0(tid[spl], "_s", seg[spl])
    if (min_len > 1L) {                               # drop segments shorter than min_len
      u    <- unique(new_tid); mi <- match(new_tid, u)
      keep <- tabulate(mi, length(u))[mi] >= min_len
      df   <- df[keep, , drop = FALSE]; new_tid <- new_tid[keep]
    }
    if (!nrow(df)) return(df[0, , drop = FALSE])
    df$tid <- new_tid
    df
  })

  # Stage 2.45: end-trim (optional; pass-through when disabled).
  # Photobleaching degrades the ENDS of a track: fewer photons (lower EFO), worse
  # precision, and blinky dark-gap-separated localizations. These inflate Rmax /
  # MSD tails and produce spurious "exports at track ends". Peel leading & trailing
  # localizations that are reached across a dark gap (≥ the dark-cycle threshold) or
  # are dim (EFO < frac × the track's median), from each end inward, until a solid,
  # well-connected localization is reached. Idempotent — a second pass stops at the
  # same boundary. Runs per segment (after gap-split), so it also cleans mid-track
  # segment ends, not just the bleaching tail.
  trk_end_trimmed <- reactive({
    df <- trk_gap_split()                  # inherits the clean_trigger recompute
    if (!isTRUE(isolate(input$trk_endtrim_enable)) || nrow(df) == 0 || !"tim" %in% names(df))
      return(df)
    use_dark <- isTRUE(isolate(input$trk_endtrim_dark))
    use_efo  <- isTRUE(isolate(input$trk_endtrim_efo)) && "efo" %in% names(df)
    if (!use_dark && !use_efo) return(df)
    t_c      <- tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_)
    dark_cyc <- max(2L, as.integer(isolate(input$trk_endtrim_dark_cycles %||% 3L)))
    efo_frac <- suppressWarnings(as.numeric(isolate(input$trk_endtrim_efo_frac %||% 0.5)))
    dark_gap <- if (use_dark && is.finite(t_c) && t_c > 0) dark_cyc * t_c else Inf
    tf  <- match(df$tid, unique(df$tid))                  # integer track codes (fast sort/group)
    ord <- order(tf, df$tim, method = "radix")
    df  <- df[ord, , drop = FALSE]
    sp  <- split(seq_len(nrow(df)), tf[ord])              # group by integer, not character
    keep <- logical(nrow(df))
    for (idx in sp) {
      m <- length(idx)
      if (m < 2L) { keep[idx] <- TRUE; next }
      tim <- df$tim[idx]
      d   <- diff(tim)                                   # d[j] = gap idx[j]→idx[j+1]
      efo_thr <- if (use_efo) efo_frac * stats::median(df$efo[idx], na.rm = TRUE) else -Inf
      efo_v   <- if (use_efo) df$efo[idx] else rep(Inf, m)
      lo <- 1L
      while (lo < m &&
             ((use_efo && is.finite(efo_v[lo]) && efo_v[lo] < efo_thr) ||
              (use_dark && d[lo]      >= dark_gap))) lo <- lo + 1L
      hi <- m
      while (hi > lo &&
             ((use_efo && is.finite(efo_v[hi]) && efo_v[hi] < efo_thr) ||
              (use_dark && d[hi - 1L] >= dark_gap))) hi <- hi - 1L
      if (hi >= lo) keep[idx[lo:hi]] <- TRUE
    }
    df[keep, , drop = FALSE]
  })

  # Stage 2.5: track-CONTINUITY gate (optional; pass-through when disabled).
  # MINFLUX tracks vary in quality: a well-tracked molecule is localized nearly
  # every cycle, whereas a dark-cycle-heavy track has long gaps → larger positional
  # uncertainty, mis-linking risk, a bleaching-degraded tail. Rather than a hard
  # integrated/dark cut, each step is weighted by its dark-cycle count with a
  # tunable tolerance (half-life): a few 1–2-cycle gaps barely hurt, long dark gaps
  # do. Keeps tracks whose mean continuity score ≥ the threshold.
  trk_quality_filtered <- reactive({
    df <- trk_end_trimmed()                # inherits the clean_trigger recompute
    if (!isTRUE(isolate(input$trk_quality_enable)) || nrow(df) == 0 || !"tim" %in% names(df))
      return(df)
    t_c <- tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_)
    if (!is.finite(t_c) || t_c <= 0) return(df)
    thr <- max(0, min(1, suppressWarnings(as.numeric(isolate(input$trk_quality_min_integ) %||% 0.6))))
    H   <- suppressWarnings(as.numeric(isolate(input$trk_quality_darktol) %||% 3))
    sc  <- .continuity_scores(df, t_c, H)
    if (is.null(sc)) return(df)
    keep <- names(sc)[is.finite(sc) & sc >= thr]
    keep <- union(keep, setdiff(unique(df$tid), names(sc)))  # keep unjudged single-loc tracks
    df[df$tid %in% keep, , drop = FALSE]
  })

  # Stage 5: trace ID length — localization count per TID, measured AFTER cleaning
  # (split/trim/quality) and the metadata filter.
  trk_len_filtered <- reactive({
    df  <- trk_meta_filtered()
    rng <- input$trk_tid_len
    if (!is.null(rng) && nrow(df) > 0 && "tid" %in% names(df)) {
      tl <- table(df$tid)
      lv <- as.numeric(tl)
      # Stale guard: pass through if range has no overlap with current lengths.
      if (!(rng[1] > max(lv) || rng[2] < min(lv))) {
        keep <- names(tl)[lv >= rng[1] & lv <= rng[2]]
        df   <- df[df$tid %in% keep, , drop = FALSE]
      }
    }
    df
  })

  # Stage 4: trace duration — max(tim)-min(tim) per TID measured AFTER stages 1,2 and 3.
  trk_dur_filtered <- reactive({
    df  <- trk_len_filtered()
    # On upload the flag is TRUE: bypass duration filtering entirely so the stale
    # input$trk_dur_filter value (from the previous dataset) is never applied to
    # the new data. The flag is cleared by trk_dur_filter_ui after it re-renders
    # with the new full range; from that point normal filtering resumes.
    # Read with isolate() — Stage 3 must NOT take a reactive dep on the flag.
    if (shiny::isolate(reset_dur_filter_flag())) return(df)
    rng <- input$trk_dur_filter
    if (!is.null(rng) && nrow(df) > 0 && "tim" %in% names(df)) {
      lo <- rng[1]; hi <- rng[2]
      dur_tbl <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(dur = max(tim) - min(tim), .groups = "drop")
      dv <- dur_tbl$dur
      # Apply unless the range already covers the full extent (a no-op). An extreme
      # bound (e.g. lower = 200 ms when most tracks are shorter) is HONOURED — it is
      # not treated as "stale". Dataset switches are handled by the reset flag above.
      if (length(dv) > 0 && !(lo <= min(dv) && hi >= max(dv))) {
        keep_tids <- dur_tbl$tid[dv >= lo & dv <= hi]
        df <- df[df$tid %in% keep_tids, , drop = FALSE]
      }
    }
    df
  })

  # Stage 5: spatial crop — keeps TIDs where >= pct_thr fraction of their
  # localizations fall inside the user-defined view box. Passes through unchanged
  # when no spatial filter is active. This is the LAST stage; trk() and the plot
  # read from here.
  trk_spatial_filtered <- reactive({
    df <- trk_dur_filtered()
    sf <- trk_spatial_filter()
    if (is.null(sf) || nrow(df) == 0) return(df)
    in_box <- df$x >= sf$x[1] & df$x <= sf$x[2] &
              df$y >= sf$y[1] & df$y <= sf$y[2]
    tids   <- unique(df$tid)
    keep   <- vapply(tids, function(tid) {
      idx <- df$tid == tid
      mean(in_box[idx], na.rm = TRUE) >= sf$pct_thr
    }, logical(1))
    df[df$tid %in% tids[keep], , drop = FALSE]
  })

  # Debounced views of the filter chain — UIs and the histogram use these so they
  # don't recompute on every slider tick; they settle with fire_frames (400 ms).
  # The length slider auto-detects from the split output (gap_split_d); the
  # duration slider auto-detects from the length-filtered output (len_filtered_d) —
  # so each slider's range reflects the LIVE previous stage (incl. split segments).
  trk_gap_split_d        <- debounce(trk_gap_split,        350)
  trk_end_trimmed_d      <- debounce(trk_end_trimmed,      350)
  trk_quality_filtered_d <- debounce(trk_quality_filtered, 350)
  trk_meta_filtered_d    <- debounce(trk_meta_filtered,    350)
  trk_len_filtered_d     <- debounce(trk_len_filtered,     350)
  trk_dur_filtered_d     <- debounce(trk_dur_filtered,     350)
  # ===========================  END FILTER PIPELINE  =========================

  # ===========================================================================
  #  CYCLE-TIME & TSA SUPPORT  (read by Stage 2 above and the TSA section)
  #  Cycle-time estimate (cluster / hardware / manual), Time-Signature Analysis
  #  readouts, sequence-JSON import, and the adaptive gap-split threshold.
  # ===========================================================================

  # Half-sample mode (Robertson–Cryer / Bickel): a bin-free Δt mode
  # estimator.  Keeps the densest half-sample until ≤3 points remain.
  half_sample_mode <- function(x) {
    x <- sort(x[is.finite(x)]); n <- length(x)
    if (n == 0L) return(NA_real_)
    while (n > 3L) {
      h <- (n + 1L) %/% 2L
      w <- x[h:n] - x[1:(n - h + 1L)]
      j <- which.min(w)
      x <- x[j:(j + h - 1L)]; n <- length(x)
    }
    if (n == 1L) return(x[1])
    if (n == 2L) return(mean(x))
    d1 <- x[2] - x[1]; d2 <- x[3] - x[2]
    if (d1 <= d2) mean(x[1:2]) else mean(x[2:3])
  }

  # Typical Δt (ms) of the RAW data — the empirical acquisition interval, bin-free
  # via the half-sample mode. Computed on trk_raw_df (upstream of ALL cleaning) so
  # gap-split / end-trim / quality can use the cycle time with no circular
  # dependency, and so metadata-removed localizations don't distort the Δt.
  trk_pipeline_modal_dt <- reactive({
    df <- tryCatch(trk_raw_df(), error = function(e) NULL)
    if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df))) return(NA_real_)
    df  <- df[order(df$tid, df$tim), ]
    dts <- unlist(tapply(df$tim, df$tid, function(t) diff(sort(t))), use.names = FALSE)
    dts <- dts[is.finite(dts) & dts > 0]
    if (!length(dts)) return(NA_real_)
    half_sample_mode(dts)
  })

  # ---- Data-driven cycle estimate (a light TSA) ------------------------------
  # Groups consecutive Δt into integer-cycle CLUSTERS and takes the median
  # spacing between adjacent cluster centres as the cycle time.  Because the
  # cluster *spacing* is the cycle time regardless of which cluster is most
  # populated, this recovers the true (atomic) cycle even when the data are poor
  # and the typical interval is itself 2–3 cycles (dense regions / frequent
  # dropouts) — exactly where the half-sample mode would mis-report the cycle.
  # Returns list(cycle, sw, n_clusters, centers, modal, frac_one_cycle).
  # Reads trk_raw_df (upstream of all cleaning → no circular dependency; raw Δt is
  # the true acquisition timing).
  trk_cluster_cycle <- reactive({
    df <- tryCatch(trk_raw_df(), error = function(e) NULL)
    if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df))) return(NULL)
    df  <- df[order(df$tid, df$tim), ]
    dts <- unlist(tapply(df$tim, df$tid, function(t) diff(sort(t))), use.names = FALSE)
    dts <- dts[is.finite(dts) & dts > 0]
    if (length(dts) < 50L) return(NULL)
    dts <- sort(dts)
    modal <- half_sample_mode(dts)
    rough <- as.numeric(quantile(dts, 0.02))          # rough smallest-cycle scale
    if (!is.finite(rough) || rough <= 0) return(NULL)
    thr  <- 0.4 * rough                                # gap threshold for clustering
    grp  <- cumsum(c(0, diff(dts) > thr))
    cl   <- split(dts, grp)
    centers <- vapply(cl, mean, numeric(1))
    counts  <- vapply(cl, length, integer(1))
    keep    <- counts >= max(5L, as.integer(0.005 * length(dts)))
    centers <- sort(centers[keep])
    if (length(centers) < 2L) {
      cyc <- if (length(centers) == 1L) centers[1] else rough
      return(list(cycle = cyc, sw = 0, n_clusters = length(centers),
                  centers = centers, modal = modal, frac_one_cycle = NA_real_))
    }
    cyc <- stats::median(diff(centers))                # spacing = atomic cycle
    sw  <- max(0, centers[1] - cyc)                    # intercept (one-time offset)
    # fraction of intervals that are a single (atomic) cycle
    n_cyc_each <- pmax(1L, as.integer(round((dts - sw) / cyc)))
    list(cycle = cyc, sw = sw, n_clusters = length(centers), centers = centers,
         modal = modal, frac_one_cycle = mean(n_cyc_each == 1L))
  })

  # Regression-based cycle/SW estimate — Bela's EOD/SW fit, reduced to a line fit.
  # The cluster centres lie at ~ k·cycle + SW for integer cycle counts k, so a
  # least-squares fit of centre ~ k gives slope = cycle time and intercept = SW.
  # This is the same quantity Bela obtains by fitting EOD_OVERHEAD / SW_OVERHEAD to
  # the time signal (EOD is just a reparametrisation of the slope through the
  # hardware vertices). Using ALL centres at once is more precise than the median
  # spacing and is robust to an occasional missing cycle. Returns list(cycle, sw,
  # n, r2, eod) — eod filled only when the hardware vertices/dwell are known.
  trk_regression_cycle <- reactive({
    cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    if (is.null(cc) || length(cc$centers) < 2L ||
        !is.finite(cc$cycle) || cc$cycle <= 0) return(NULL)
    centers <- sort(as.numeric(cc$centers))
    idx <- round(centers / cc$cycle)                   # integer cycle index per centre
    keep <- is.finite(idx) & idx >= 1L & !duplicated(idx)
    centers <- centers[keep]; idx <- idx[keep]
    if (length(unique(idx)) < 2L) return(NULL)
    fit <- tryCatch(stats::lm(centers ~ idx), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    cf <- stats::coef(fit)
    slope <- unname(cf[2]); intercept <- unname(cf[1])
    if (!is.finite(slope) || slope <= 0) return(NULL)
    r2 <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
    # Optional EOD per move, if the hardware geometry is supplied: cycle = dwell +
    # EOD × vertices  →  EOD = (cycle − dwell) / vertices  (µs).
    eod <- NA_real_
    if (isTRUE(input$trk_hw_enable)) {
      dwell_us <- suppressWarnings(as.numeric(input$trk_hw_dwell_us))
      pat <- input$trk_hw_pattern %||% "6"
      nv  <- if (identical(pat, "other")) suppressWarnings(as.integer(input$trk_hw_vertices))
             else suppressWarnings(as.integer(pat))
      if (is.finite(dwell_us) && is.finite(nv) && nv >= 1L)
        eod <- (slope * 1000 - dwell_us) / nv
    }
    list(cycle = slope, sw = max(0, intercept), n = length(idx), r2 = r2, eod = eod)
  })

  # Cycle time (ms) from the hardware model — TSA formula derived from Vogler:
  #   cycle = T_DWELL + EOD_OVERHEAD × NUM_VERTICES

  trk_hw_cycle_time <- reactive({
    if (!isTRUE(input$trk_hw_enable)) return(NA_real_)
    dwell_us <- suppressWarnings(as.numeric(input$trk_hw_dwell_us))
    eod_us   <- suppressWarnings(as.numeric(input$trk_hw_eod_us))
    pat      <- input$trk_hw_pattern %||% "6"
    nv <- if (identical(pat, "other"))
            suppressWarnings(as.integer(input$trk_hw_vertices))
          else suppressWarnings(as.integer(pat))
    if (!is.finite(dwell_us) || dwell_us <= 0 || !is.finite(nv) || nv < 1) return(NA_real_)
    if (!is.finite(eod_us) || eod_us < 0) eod_us <- 0
    (dwell_us + eod_us * nv) / 1000   # µs → ms
  })

  # ---- SINGLE SOURCE OF TRUTH for the acquisition cycle time (ms) -------------
  # Priority: (1) hardware model if enabled & valid, (2) manual override in the
  # Advanced MSD panel if set > 0, (3) data-driven CLUSTER cycle (robust atomic
  # cycle), (4) half-sample-mode modal Δt as a last resort.  Everything that
  # needs a cycle time reads this so they can never disagree.
  trk_cycle_time <- reactive({
    hw <- trk_hw_cycle_time()
    if (is.finite(hw) && hw > 0) return(hw)
    man <- suppressWarnings(as.numeric(input$trk_msd_cycle_ms))
    if (is.finite(man) && man > 0) return(man)
    # Optional: Bela-style regression cycle time (toggle), preferred over the
    # cluster median when enabled and valid.
    if (isTRUE(input$trk_cycle_use_regression)) {
      rc <- tryCatch(trk_regression_cycle(), error = function(e) NULL)
      if (!is.null(rc) && is.finite(rc$cycle) && rc$cycle > 0) return(rc$cycle)
    }
    cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    if (!is.null(cc) && is.finite(cc$cycle) && cc$cycle > 0) return(cc$cycle)
    trk_pipeline_modal_dt()
  })

  # Standard-MSD fit window bounds, entered in CYCLES, converted to ms via the
  # cycle time (Vogler's fit_limits, default 5–50).  0 means "off".
  msd_min_lag_ms <- reactive({
    cyc <- suppressWarnings(as.numeric(input$trk_msd_min_lag_cyc))
    tc  <- trk_cycle_time()
    if (length(cyc) == 0 || !is.finite(cyc) || cyc <= 0 || !is.finite(tc) || tc <= 0) return(0)
    cyc * tc
  })
  msd_fit_max_ms <- reactive({
    cyc <- suppressWarnings(as.numeric(input$trk_msd_fit_max_cyc))
    tc  <- trk_cycle_time()
    if (length(cyc) == 0 || !is.finite(cyc) || cyc <= 0 || !is.finite(tc) || tc <= 0) return(0)
    cyc * tc
  })

  # Status / cross-check line shown in the Acquisition Parameters panel
  output$trk_cycle_time_status <- renderUI({
    req(input$trk_file)
    mdt <- trk_pipeline_modal_dt()
    hw  <- trk_hw_cycle_time()
    tc  <- trk_cycle_time()
    rc0 <- tryCatch(trk_regression_cycle(), error = function(e) NULL)
    reg_active <- isTRUE(input$trk_cycle_use_regression) &&
                  !is.null(rc0) && is.finite(rc0$cycle) && rc0$cycle > 0
    src <- if (is.finite(hw) && hw > 0) "hardware model"
           else {
             man <- suppressWarnings(as.numeric(input$trk_msd_cycle_ms))
             if (is.finite(man) && man > 0) "manual override"
             else if (reg_active) "regression fit (EOD/SW)"
             else "cluster cycle (data)"
           }
    if (!is.finite(tc) || tc <= 0)
      return(tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
                    "Cycle time: not yet determined (load data)."))
    rows <- list(
      tags$p(style = "font-size:0.84em; color:#2070a0; margin:2px 0;",
             sprintf("Cycle time in use: %.4g ms  (source: %s)", tc, src))
    )
    if (is.finite(hw) && hw > 0 && is.finite(mdt) && mdt > 0) {
      # The atomic cycle (model) is the MINIMUM interval, so it should sit
      # slightly BELOW the typical/modal Δt (which includes small per-cycle
      # variability).  ratio = model/modal ≈ 0.8–1.0 is the healthy regime.
      ratio <- hw / mdt
      verdict <- if (ratio >= 0.75 && ratio <= 1.05)
                   "consistent ✓ (atomic cycle just below typical Δt, as expected)"
                 else if (ratio > 1.05)
                   "⚠ model exceeds typical Δt - dwell/overhead likely too high"
                 else
                   "⚠ model well below typical Δt - check dwell, or many dark cycles inflate Δt"
      col <- if (ratio >= 0.75 && ratio <= 1.05) "#2a8a3e"
             else if (ratio > 1.05) "#c0392b" else "#b07000"
      msg <- sprintf("Atomic cycle (model) %.4g ms vs typical Δt %.4g ms (ratio %.2f) — %s",
                     hw, mdt, ratio, verdict)
      rows <- c(rows, list(tags$p(style = sprintf("font-size:0.82em; color:%s; margin:2px 0;", col), msg)))
    } else if (is.finite(mdt) && mdt > 0) {
      rows <- c(rows, list(tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
             sprintf("Empirical modal Δt = %.4g ms.", mdt))))
    }
    # Data-driven cluster cycle + data-quality verdict (a light TSA). The verdict
    # uses the SINGLE-CYCLE (integrated) fraction with the SAME thresholds as the
    # TSA / Data Quality panels (≥80% GOOD, ≥50% MODERATE), so all three agree. The
    # modal/cluster ratio is shown alongside only as a cycle-cadence diagnostic.
    # This is the RAW fraction (cluster_cycle reads raw data).
    cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    if (!is.null(cc) && is.finite(cc$cycle) && cc$cycle > 0 && is.finite(cc$modal)) {
      ratio <- cc$modal / cc$cycle                # typical interval / atomic cycle (diagnostic)
      fo    <- cc$frac_one_cycle
      verdict <- if (is.finite(fo)) {
          if (fo >= 0.8)
            sprintf("data quality (raw): GOOD — %.0f%% single-cycle (integrated) steps.", fo * 100)
          else if (fo >= 0.5)
            sprintf("data quality (raw): MODERATE — %.0f%% single-cycle steps; some dark-cycle dropouts.", fo * 100)
          else
            sprintf("data quality (raw): POOR — %.0f%% single-cycle steps; frequent dark cycles.", fo * 100)
        } else if (is.finite(ratio)) {
          if (ratio < 1.4) "data quality (raw): GOOD (typical interval ≈ 1 cycle)."
          else if (ratio < 2.3) sprintf("data quality (raw): MODERATE (typical interval ≈ %.1f cycles).", ratio)
          else sprintf("data quality (raw): POOR (typical interval ≈ %.1f cycles).", ratio)
        } else ""
      good <- if (is.finite(fo)) fo >= 0.8 else is.finite(ratio) && ratio < 1.4
      mod  <- if (is.finite(fo)) fo >= 0.5 else is.finite(ratio) && ratio < 2.3
      qcol <- if (good) "#2a8a3e" else if (mod) "#b07000" else "#c0392b"
      rows <- c(rows, list(
        tags$p(style = "font-size:0.82em; color:#555; margin:2px 0;",
               sprintf("Cluster (atomic) cycle = %.4g ms  |  typical interval (HSM) = %.4g ms  |  %d cycle-clusters detected",
                       cc$cycle, cc$modal, cc$n_clusters)),
        tags$p(style = sprintf("font-size:0.82em; color:%s; margin:2px 0;", qcol), verdict)))
    }
    # Regression (EOD/SW) cycle — Bela-style line fit, shown for side-by-side
    # comparison with the cluster median (and used as the source when toggled on).
    if (!is.null(rc0) && is.finite(rc0$cycle) && rc0$cycle > 0) {
      eod_txt <- if (is.finite(rc0$eod)) sprintf("  |  EOD ≈ %.2f µs/move", rc0$eod) else ""
      r2_txt  <- if (is.finite(rc0$r2)) sprintf("  |  R² = %.4f", rc0$r2) else ""
      rows <- c(rows, list(
        tags$p(style = sprintf("font-size:0.82em; color:%s; margin:2px 0;",
                               if (reg_active) "#2070a0" else "#777"),
               sprintf("%s Regression fit: cycle = %.4g ms  |  SW = %.4g ms%s  (%d cluster centres)%s",
                       if (reg_active) "▶" else "○",
                       rc0$cycle, rc0$sw, eod_txt, rc0$n, r2_txt))))
    }
    do.call(tags$div, rows)
  })

  # ms-equivalent hint shown under the cycle-based gap threshold
  output$trk_split_gap_ms_hint <- renderUI({
    req(isTRUE(input$trk_split_gaps_enable))
    cyc <- suppressWarnings(as.numeric(input$trk_split_gap_cycles))
    mdt <- trk_cycle_time()
    if (!is.finite(cyc) || cyc <= 0 || !is.finite(mdt) || mdt <= 0) return(NULL)
    tags$p(style = "font-size:0.81em; color:#2070a0; margin:2px 0;",
           sprintf("≈ %.4g ms  (at cycle time = %.4g ms)", cyc * mdt, mdt))
  })

  # Adaptive (Bela) threshold readout: the data-driven cycle threshold plus how
  # many tracks/localizations survive the split + min-length filter.
  output$trk_split_adaptive_hint <- renderUI({
    req(isTRUE(input$trk_split_gaps_enable), isTRUE(input$trk_split_adaptive))
    thr <- trk_adaptive_gap_cycles()
    tot <- trk_tot_cycles()
    if (!is.finite(thr) || length(tot) < 2)
      return(tags$p(style = "color:#888; font-size:0.81em;", "Need data to compute the adaptive threshold."))
    med <- as.integer(stats::median(tot)); sdv <- stats::sd(tot)
    cyc <- trk_cycle_time()
    n_split <- sum(tot >= thr)
    df_in    <- trk_raw_df()   # split is the first stage now → its input is raw
    df_split <- trk_gap_split()
    df_final <- tryCatch(trk_spatial_filtered(), error = function(e) df_split)
    n_tid_in    <- if (nrow(df_in)    > 0) length(unique(df_in$tid))    else 0L
    n_tid_split <- if (nrow(df_split) > 0) length(unique(df_split$tid)) else 0L
    n_tid_final <- if (nrow(df_final) > 0) length(unique(df_final$tid)) else 0L
    tags$div(style = "font-size:0.81em; color:#2a6; background:#f0f7f0; border-radius:4px; padding:5px 7px; margin:3px 0;",
      tags$b(sprintf("Auto threshold = %d cycles", thr)),
      sprintf("  (median %d + SD %.2f + 1)", med, sdv), tags$br(),
      sprintf("≈ %.4g ms · %s of %s steps exceed it",
              thr * cyc, format(n_split, big.mark = ","), format(length(tot), big.mark = ",")), tags$br(),
      sprintf("Tracks: %s in → %s after split+min-length → %s displayed",
              format(n_tid_in, big.mark = ","), format(n_tid_split, big.mark = ","),
              format(n_tid_final, big.mark = ",")),
      if (n_tid_final < n_tid_split)
        tags$span(style = "color:#888;", tags$br(),
          "(length / duration / spatial filters remove the rest — their ranges now expand to cover segments)"))
  })

  # End-trim funnel hint (localizations removed from track ends).
  output$trk_endtrim_hint <- renderUI({
    req(isTRUE(input$trk_endtrim_enable))
    din  <- tryCatch(trk_gap_split_d(),   error = function(e) NULL)
    dout <- tryCatch(trk_end_trimmed_d(), error = function(e) NULL)
    if (is.null(din) || is.null(dout) || !nrow(din))
      return(tags$p(style = "color:#888; font-size:0.81em;", "Need data to compute end-trim."))
    n_in <- nrow(din); n_out <- nrow(dout)
    tags$div(style = "font-size:0.81em; color:#2a6; background:#f0f7f0; border-radius:4px; padding:5px 7px; margin:3px 0;",
      tags$b(sprintf("Localizations: %s → %s kept", format(n_in, big.mark = ","), format(n_out, big.mark = ","))),
      sprintf("  (%s trimmed from track ends, %.1f%%)",
              format(n_in - n_out, big.mark = ","), 100 * (n_in - n_out) / max(1L, n_in)))
  })

  # Track-quality gate funnel hint.
  output$trk_quality_hint <- renderUI({
    req(isTRUE(input$trk_quality_enable))
    din  <- tryCatch(trk_gap_split_d(),        error = function(e) NULL)
    dout <- tryCatch(trk_quality_filtered_d(), error = function(e) NULL)
    if (is.null(din) || is.null(dout) || !nrow(din))
      return(tags$p(style = "color:#888; font-size:0.81em;", "Need data to compute track quality."))
    n_in  <- length(unique(din$tid)); n_out <- length(unique(dout$tid))
    t_c   <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    thr   <- max(0, min(1, suppressWarnings(as.numeric(input$trk_quality_min_integ %||% 0.6))))
    H     <- suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3))
    tags$div(style = "font-size:0.81em; color:#2a6; background:#f0f7f0; border-radius:4px; padding:5px 7px; margin:3px 0;",
      tags$b(sprintf("Tracks: %s → %s kept (%.0f%%)",
                     format(n_in, big.mark = ","), format(n_out, big.mark = ","),
                     100 * n_out / max(1L, n_in))),
      sprintf("  at continuity ≥ %.2f (tolerance %.1f)", thr, H),
      if (is.finite(t_c)) sprintf("  ·  cycle = %.4g ms", t_c) else "")
  })

  # ── TSA: MINFLUX sequence (.json) import ────────────────────────────────────
  # Parses an Imspector imaging-sequence JSON. The acquisition cycle time is still
  # recovered from the DATA (the robust cluster method) — the sequence is used to
  # cross-check the scan geometry and, most usefully, to set the gap-split
  # threshold from the hardware's own maxOffTime ("particle lost" timeout).
  .pattern_vertices <- function(p) {
    if (is.null(p) || is.na(p)) return(NA_integer_)
    switch(tolower(p),
           hexagon = 6L, square = 4L, triangle = 3L,
           zline = 2L, zline2 = 2L, line = 2L,
           NA_integer_)
  }

  trk_seq_params <- reactive({
    req(input$trk_seq_json)
    fp <- input$trk_seq_json$datapath
    seq <- tryCatch(jsonlite::fromJSON(fp, simplifyVector = FALSE),
                    error = function(e) NULL)
    if (is.null(seq) || is.null(seq$Itr) || length(seq$Itr) == 0)
      return(list(ok = FALSE, msg = "Not a recognised MINFLUX sequence (no 'Itr' array)."))
    itrs <- seq$Itr
    getp  <- function(it, k) { v <- it[[k]]; if (is.null(v)) NA else v }
    getm  <- function(it, k) { m <- it[["Mode"]]; if (is.null(m)) return(NA); v <- m[[k]]; if (is.null(v)) NA else v }
    pats  <- vapply(itrs, function(it) { v <- getm(it, "pattern"); if (is.na(v)) NA_character_ else as.character(v) }, character(1))
    dims  <- suppressWarnings(vapply(itrs, function(it) as.integer(getm(it, "dim")), integer(1)))
    moff  <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "maxOffTime")))   # seconds
    prep  <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "patRepeat")))
    pht   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "phtLimit")))
    pdw   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "patDwellTime"))) # seconds (timeout)
    final_pat <- tail(pats[!is.na(pats)], 1L)
    if (length(final_pat) == 0) final_pat <- NA_character_
    dim_final <- suppressWarnings(as.integer(tail(dims[!is.na(dims)], 1L)))
    list(
      ok = TRUE,
      n_itr = length(itrs),
      patterns = pats,
      final_pattern = final_pat,
      vertices = .pattern_vertices(final_pat),
      dim = if (length(dim_final) && is.finite(dim_final)) dim_final else NA_integer_,
      maxOffTime_ms = if (is.finite(moff)) moff * 1000 else NA_real_,
      patRepeat = prep,
      phtLimit = pht,
      patDwellTime_ms = if (is.finite(pdw)) pdw * 1000 else NA_real_
    )
  })

  # Apply the parsed sequence: set the gap-split threshold from maxOffTime, enable
  # gap splitting, and (when the lateral pattern is unambiguous) pre-fill the
  # hardware scan pattern. The cycle time stays data-driven.
  observeEvent(input$trk_seq_json, {
    sp <- trk_seq_params()
    if (is.null(sp) || !isTRUE(sp$ok)) return(invisible(NULL))
    # gap-split threshold = round(maxOffTime / cycle time)  [in cycles]
    cyc <- isolate(trk_cycle_time())
    if (is.finite(sp$maxOffTime_ms) && is.finite(cyc) && cyc > 0) {
      n_cyc <- max(2, round(sp$maxOffTime_ms / cyc))
      updateNumericInput(session, "trk_split_gap_cycles", value = n_cyc)
      updateCheckboxInput(session, "trk_split_gaps_enable", value = TRUE)
    }
    # pre-fill the hardware pattern if it maps to a clean lateral polygon
    nv <- sp$vertices
    if (is.finite(nv) && nv %in% c(3L, 4L, 6L)) {
      updateSelectInput(session, "trk_hw_pattern", selected = as.character(nv))
    } else if (is.finite(nv)) {
      updateSelectInput(session, "trk_hw_pattern", selected = "other")
      updateNumericInput(session, "trk_hw_vertices", value = nv)
    }
  }, ignoreInit = TRUE)

  output$trk_seq_json_status <- renderUI({
    if (is.null(input$trk_seq_json)) return(NULL)
    sp <- trk_seq_params()
    if (is.null(sp) || !isTRUE(sp$ok))
      return(tags$p(style = "color:#b00; font-size:0.85em;", sp$msg %||% "Could not parse file."))
    cyc <- trk_cycle_time()
    gap_cyc <- if (is.finite(sp$maxOffTime_ms) && is.finite(cyc) && cyc > 0)
                 max(2, round(sp$maxOffTime_ms / cyc)) else NA
    upat <- unique(sp$patterns[!is.na(sp$patterns)])
    tags$div(style = "font-size:0.84em; background:#eef5ee; border:1px solid #cdd; border-radius:4px; padding:6px 8px; margin-top:-2px;",
      tags$b("Sequence loaded "), sprintf("(%d iterations)", sp$n_itr), tags$br(),
      sprintf("Patterns: %s", paste(upat, collapse = ", ")), tags$br(),
      sprintf("Tracking iteration: %s  ·  dim = %s",
              if (length(sp$final_pattern) && !is.na(sp$final_pattern)) sp$final_pattern else "?",
              if (is.finite(sp$dim)) sprintf("%dD", sp$dim) else "?"), tags$br(),
      if (is.finite(sp$maxOffTime_ms))
        sprintf("maxOffTime = %.3g ms  →  gap-split @ %s cycles",
                sp$maxOffTime_ms, ifelse(is.na(gap_cyc), "?", as.character(gap_cyc)))
      else "maxOffTime: not found",
      if (is.finite(sp$phtLimit)) tags$span(tags$br(), sprintf("Photon limit = %g  ·  patRepeat = %g",
                                                               sp$phtLimit, sp$patRepeat)) )
  })

  # ── TSA: top-line data-quality verdict ──────────────────────────────────────
  # Compares the misleading modal Δt against the robust cluster cycle. A large
  # ratio means dense regions are over-populated with multi-cycle gaps, biasing
  # any naive cycle estimate — a direct read on whether the data are clean.
  output$trk_tsa_quality <- renderUI({
    req(isTRUE(input$trk_show_tsa), tsa_ran())
    cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    if (is.null(cc) || !is.finite(cc$cycle) || cc$cycle <= 0)
      return(tags$p(style = "color:#888;", "Upload tracking data to run the time-signature analysis."))
    modal <- cc$modal
    ratio <- if (is.finite(modal) && cc$cycle > 0) modal / cc$cycle else NA_real_
    frac1 <- cc$frac_one_cycle                       # RAW (metadata-filtered)
    # Verdict on the POST-FILTERED single-cycle fraction — the SAME data the
    # Integrated/Dark charts below use (dt_vectors → trk(), i.e. after the quality
    # gate / end-trim / length / duration / spatial filters). frac1 (raw) is kept as
    # a diagnostic row so the filter improvement is visible here too. Same metric
    # and thresholds as the Calibration panel, so every verdict agrees.
    ccb       <- tryCatch(cycle_class_bundle(), error = function(e) NULL)
    frac_filt <- if (!is.null(ccb) && is.finite(ccb$frac_integ)) ccb$frac_integ else NA_real_
    frac_use  <- if (is.finite(frac_filt)) frac_filt else frac1
    verdict <- if (is.finite(frac_use)) {
                 if (frac_use >= 0.8) c("GOOD", "#1a8a1a")
                 else if (frac_use >= 0.5) c("MODERATE", "#c08000")
                 else c("POOR", "#b00")
               } else if (is.finite(ratio)) {
                 if (ratio < 1.4) c("GOOD", "#1a8a1a")
                 else if (ratio < 3.0) c("MODERATE", "#c08000")
                 else c("POOR", "#b00")
               } else c("UNKNOWN", "#888")
    cyc_in_use <- trk_cycle_time()
    tags$div(style = sprintf("border-left:5px solid %s; padding:6px 10px; background:#fafafa;", verdict[2]),
      tags$div(style = sprintf("font-size:1.05em; font-weight:bold; color:%s;", verdict[2]),
               sprintf("Data quality (%s): %s",
                       if (is.finite(frac_filt)) "filtered" else "raw", verdict[1])),
      tags$table(style = "font-size:0.86em; margin-top:4px;",
        tags$tr(tags$td(style="padding-right:12px;", "Atomic cycle (cluster)"),
                tags$td(tags$b(sprintf("%.4g ms", cc$cycle)))),
        tags$tr(tags$td("Modal Δt (raw)"),
                tags$td(sprintf("%.4g ms", modal))),
        tags$tr(tags$td("Modal / cluster ratio"),
                tags$td(tags$b(sprintf("%.2f×", ratio)))),
        if (is.finite(frac_filt)) tags$tr(tags$td("Single-cycle fraction (filtered)"),
                tags$td(tags$b(sprintf("%.0f%%", 100 * frac_filt)))),
        if (is.finite(frac1)) tags$tr(tags$td("Single-cycle fraction (raw)"),
                tags$td(sprintf("%.0f%%", 100 * frac1))),
        tags$tr(tags$td("Cycle time used downstream"),
                tags$td(tags$b(sprintf("%.4g ms", cyc_in_use))))),
      tags$p(style = "font-size:0.8em; color:#666; margin:5px 0 0 0;",
        if (identical(verdict[1], "GOOD"))
          "Most steps are single-cycle (integrated): gap-free tracking, well-conditioned for diffusion analysis."
        else if (identical(verdict[1], "MODERATE"))
          "A sizeable minority of steps span 2+ cycles (dark cycles). Enable the track-quality gate / end-trim (Display tab) to clean it — the Calibration panel shows the after-filter improvement."
        else if (identical(verdict[1], "POOR"))
          "Most steps span multiple cycles — frequent dropouts. Trust the cluster cycle (not the modal Δt), and filter dark-heavy tracks/ends before trusting D."
        else "Insufficient structure to judge. Modal/cluster ratio is used as a fallback cadence check."))
  })

  # Per-interval TOTAL acquisition-cycle counts across all meta-filtered tracks:
  # TOT_i = round((Δt_i − SW_overhead) / cycle_time).  This is Bela's "TOT" signal
  # (integrated + dark cycles spanned by each consecutive step).
  trk_tot_cycles <- reactive({
    df  <- trk_raw_df()   # split operates on raw → estimate its threshold on raw
    cyc <- trk_cycle_time()
    if (nrow(df) == 0 || !is.finite(cyc) || cyc <= 0) return(integer(0))
    cc  <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    sw  <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0
    df  <- df[order(df$tid, df$tim), ]
    # Bela computes suggest_split_thresholds AFTER filter_by_length(50), so the
    # threshold is estimated only from tracks that meet the minimum length.
    # Mirror that here when adaptive splitting is on, so our med+SD+1 matches his.
    if (isTRUE(input$trk_split_adaptive)) {
      min_len <- max(1L, as.integer(input$trk_split_min_len %||% 50L))
      tl <- table(df$tid); keep_tid <- names(tl)[tl >= min_len]
      df <- df[df$tid %in% keep_tid, , drop = FALSE]
      if (nrow(df) == 0) return(integer(0))
    }
    dts <- unlist(lapply(split(df$tim, df$tid),
                         function(t) if (length(t) > 1L) diff(t) else numeric(0)),
                  use.names = FALSE)
    tot <- round((dts - sw) / cyc)
    tot <- tot[is.finite(tot) & tot >= 1]
    as.integer(tot)
  })

  # Bela's adaptive split threshold (in cycles): med_th = floor(median + std + 1)
  # of the TOT signal, capped by the hardware "particle lost" ceiling
  # high_th = round(maxOffTime / cycle).  maxOffTime comes from the uploaded
  # sequence JSON if present, else Bela's 3 ms default.
  trk_adaptive_gap_cycles <- reactive({
    tot <- trk_tot_cycles()
    if (length(tot) < 2L) return(NA_integer_)
    med <- as.integer(stats::median(tot))
    sdv <- stats::sd(tot); if (!is.finite(sdv)) sdv <- 0
    med_th <- as.integer(floor(med + sdv + 1))
    cyc <- trk_cycle_time()
    moff_ms <- tryCatch({ sp <- trk_seq_params(); if (!is.null(sp) && isTRUE(sp$ok) &&
                            is.finite(sp$maxOffTime_ms)) sp$maxOffTime_ms else 3.0 },
                        error = function(e) 3.0)
    high_th <- if (is.finite(cyc) && cyc > 0) as.integer(round(moff_ms / cyc)) else NA_integer_
    if (is.finite(high_th) && high_th >= 2L) min(med_th, high_th) else med_th
  })

  # (The filter pipeline stages 1.5–4 are defined together near Stage 1, under
  # the "TRACKING FILTER PIPELINE" banner. trk_tot_cycles / trk_adaptive_gap_cycles
  # above are the cycle-time support those stages read.)

  output$trk_chunk_stat <- renderUI({
    req(input$trk_chunk_sel)
    if (input$trk_chunk_sel == "__all__") {
      return(tags$div(style = "color:#555;", "All TIDs selected  -  loading full file."))
    }
    df <- tryCatch({
      if (grepl("\\.qs$", input$trk_chunk_sel, TRUE)) qs::qread(input$trk_chunk_sel) else
        if (grepl("\\.fst$", input$trk_chunk_sel, TRUE)) fst::read_fst(input$trk_chunk_sel) |> as.data.frame() else NULL
    }, error = function(e) NULL)
    if (is.null(df)) return(NULL)
    tids <- length(unique(as.character(df$tid)))
    tags$div(style="color:#555;",
             paste0("Selected chunk: ", format(nrow(df), big.mark=","), " rows * ",
                    tids, " TIDs"))
  })
  
  # ---- Per-TID summary (from tid_summary.csv saved at chunking time) ----------
  tid_summary_rv <- reactive({
    req(!is.null(rv$manifest_dir))
    f <- file.path(rv$manifest_dir, "tid_summary.csv")
    if (!file.exists(f)) return(NULL)
    readr::read_csv(f, show_col_types = FALSE)
  })

  # Lifetime summary stat shown above the metadata filter section
  output$trk_lifetime_stat <- renderUI({
    ts <- tid_summary_rv()
    if (is.null(ts) || !nrow(ts)) return(NULL)
    avg <- mean(ts$duration_ms, na.rm = TRUE)
    med <- stats::median(ts$duration_ms, na.rm = TRUE)
    n   <- nrow(ts)
    tags$div(
      style = "color:#555; margin-bottom:4px; font-size:0.9em;",
      tags$b("Trace lifetimes \u2014"),
      sprintf("%d traces \u00B7 mean %.0f ms \u00B7 median %.0f ms", n, avg, med)
    )
  })

  # Stage 3 (duration) auto-detection. BOUNDS reflect the length-filtered population
  # but update at PLOT cadence only — keyed on data_ver() (upload / Plot Data / Reset
  # Ranges), NOT live. So the range describes the currently PLOTTED data and stops
  # jumping as sliders move. The VALUE is PRESERVED (clamped into the new bounds) so a
  # manual minimum (e.g. 50 ms) sticks; it snaps to full only on new upload / Reset
  # Ranges / post-rechunk, signalled by reset_dur_filter_flag.
  output$trk_dur_filter_ui <- renderUI({
    reset_ranges_trigger()
    data_ver()                                                           # update at PLOT cadence (upload / Plot Data), not live
    df <- tryCatch(isolate(trk_len_filtered()), error = function(e) NULL)
    dv <- if (!is.null(df) && nrow(df) > 0 && all(c("tid","tim") %in% names(df)))
            (df %>% dplyr::group_by(tid) %>%
               dplyr::summarise(dur = max(tim) - min(tim), .groups = "drop"))$dur
          else numeric(0)
    if (!length(dv)) {                          # pipeline not ready → rechunk stats, else upload
      rs_re <- trk_rechunk_stats_rv(); rs_or <- trk_raw_stats_rv()
      dv <- if (!is.null(rs_re) && length(rs_re$tid_duration)) rs_re$tid_duration
            else if (!is.null(rs_or) && length(rs_or$tid_duration)) rs_or$tid_duration else numeric(0)
    }
    if (!length(dv)) return(NULL)
    bounds <- range(dv, na.rm = TRUE)
    curr <- isolate(input$trk_dur_filter)
    snap <- isTRUE(isolate(reset_dur_filter_flag())) || is.null(curr) ||
            length(curr) != 2 || any(!is.finite(as.numeric(curr)))
    val  <- if (snap) bounds
            else c(max(bounds[1], min(as.numeric(curr))), min(bounds[2], max(as.numeric(curr))))
    shiny::isolate(reset_dur_filter_flag(FALSE))
    numericRangeInput("trk_dur_filter", "Trace duration filter (ms)",
                      value = val, min = bounds[1], max = bounds[2])
  })

  # ---- Filter & Rechunk -------------------------------------------------------
  # Always reads from original_manifest_dir (frozen at upload) so repeated
  # rechunks are idempotent — Reset Ranges + Rechunk restores the full dataset.
  # Writes filtered result to a new dir that becomes the new manifest_dir.
  observeEvent(input$trk_rechunk_filtered, {
    req(!is.null(rv$original_manifest_dir))
    manf <- file.path(rv$original_manifest_dir, "manifest.csv")
    req(file.exists(manf))

    withProgress(message = "Filter & Rechunk: loading chunks...", value = 0, {
      man <- readr::read_csv(manf, show_col_types = FALSE)
      n_chunks <- nrow(man)

      # Split parameters computed ONCE so every chunk uses the same threshold and
      # they MATCH the live trk_gap_split (adaptive → median+SD+1 cap maxOffTime,
      # else the manual cycle count; the same cycle time; and min-length drop in
      # adaptive mode). Each TID lives wholly in one chunk, so per-chunk splitting
      # with this shared threshold reproduces the live result exactly.
      g_split_on  <- isolate(isTRUE(input$trk_split_gaps_enable))
      g_adaptive  <- isolate(isTRUE(input$trk_split_adaptive))
      g_cyc       <- if (g_adaptive) isolate(trk_adaptive_gap_cycles())
                     else isolate(suppressWarnings(as.numeric(input$trk_split_gap_cycles)))
      g_cycletime <- isolate(trk_cycle_time())
      g_minlen    <- max(1L, isolate(as.integer(input$trk_split_min_len %||% 50L)))  # both modes now
      g_thr       <- if (g_split_on && is.finite(g_cyc) && g_cyc > 0 &&
                         is.finite(g_cycletime) && g_cycletime > 0) g_cyc * g_cycletime else NA_real_

      # End-trim + quality-gate params (match live trk_end_trimmed / trk_quality_filtered).
      et_on       <- isolate(isTRUE(input$trk_endtrim_enable))
      et_dark     <- isolate(isTRUE(input$trk_endtrim_dark))
      et_efo      <- isolate(isTRUE(input$trk_endtrim_efo))
      et_dark_cyc <- max(2L, isolate(as.integer(input$trk_endtrim_dark_cycles %||% 3L)))
      et_efo_frac <- isolate(suppressWarnings(as.numeric(input$trk_endtrim_efo_frac %||% 0.5)))
      et_dark_gap <- if (et_dark && is.finite(g_cycletime) && g_cycletime > 0) et_dark_cyc * g_cycletime else Inf
      q_on        <- isolate(isTRUE(input$trk_quality_enable))
      q_thr       <- max(0, min(1, isolate(suppressWarnings(as.numeric(input$trk_quality_min_integ %||% 0.6)))))
      q_H         <- isolate(suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3)))
      trim_ends_one <- function(tb) {   # kept row indices of one (time-sorted) track
        m <- nrow(tb); if (m < 2L) return(seq_len(m))
        d <- diff(tb$tim)
        efo_thr <- if (et_efo && "efo" %in% names(tb)) et_efo_frac * stats::median(tb$efo, na.rm = TRUE) else -Inf
        efo_v   <- if (et_efo && "efo" %in% names(tb)) tb$efo else rep(Inf, m)
        lo <- 1L
        while (lo < m && ((et_efo && is.finite(efo_v[lo]) && efo_v[lo] < efo_thr) ||
                          (et_dark && d[lo]      >= et_dark_gap))) lo <- lo + 1L
        hi <- m
        while (hi > lo && ((et_efo && is.finite(efo_v[hi]) && efo_v[hi] < efo_thr) ||
                           (et_dark && d[hi - 1L] >= et_dark_gap))) hi <- hi - 1L
        if (hi >= lo) lo:hi else integer(0)
      }

      # --- Pass 1: per chunk, in the LIVE pipeline order:
      #     split → end-trim → quality gate → metadata → length ---
      all_dfs <- vector("list", n_chunks)
      for (i in seq_len(n_chunks)) {
        incProgress(i / (n_chunks * 2), detail = paste("Reading chunk", i))
        f <- man$file[i]
        ch <- tryCatch(
          if (grepl("\\.qs$",  f, TRUE)) qs::qread(f)
          else if (grepl("\\.fst$", f, TRUE)) fst::read_fst(f) |> as.data.frame()
          else NULL,
          error = function(e) NULL
        )
        if (is.null(ch) || !nrow(ch)) next
        # Same post-processing as trk_raw_df
        ch$tid <- as.character(ch$tid)
        ch$x   <- as.numeric(ch$x)
        ch$y   <- as.numeric(ch$y)
        ch$tim <- as.numeric(ch$tim) * 1000  # seconds -> ms

        # Stage 1: gap splitting — shared threshold/min-length (matches live).
        if (g_split_on && is.finite(g_thr) && g_thr > 0 && nrow(ch) > 0) {
          ch  <- ch[order(ch$tid, ch$tim), ]
          spc <- split(ch, ch$tid)
          new_segs <- vector("list", length(spc) * 2L)
          ks <- 0L
          for (snm in names(spc)) {
            tb2 <- spc[[snm]]
            if (nrow(tb2) < 2L) {
              if (nrow(tb2) >= g_minlen) { ks <- ks+1L; new_segs[[ks]] <- tb2 }
              next
            }
            gp <- which(diff(tb2$tim) > g_thr)
            if (!length(gp)) {
              if (nrow(tb2) >= g_minlen) { ks <- ks+1L; new_segs[[ks]] <- tb2 }
              next
            }
            st2 <- c(1L, gp+1L); en2 <- c(gp, nrow(tb2))
            for (ii in seq_along(st2)) {
              seg2 <- tb2[st2[ii]:en2[ii], , drop=FALSE]
              if (nrow(seg2) < g_minlen) next
              seg2$tid <- paste0(snm, "_s", ii)
              ks <- ks+1L; new_segs[[ks]] <- seg2
            }
          }
          ch <- if (ks > 0L) dplyr::bind_rows(new_segs[seq_len(ks)]) else ch[0, , drop = FALSE]
        }
        # Stage 2: end-trim (peel degraded track ends).
        if (et_on && (et_dark || et_efo) && nrow(ch) > 0) {
          ch <- ch[order(ch$tid, ch$tim), ]
          spt <- split(seq_len(nrow(ch)), ch$tid)
          keep_et <- logical(nrow(ch))
          for (idx in spt) {
            k <- trim_ends_one(ch[idx, , drop = FALSE])
            if (length(k)) keep_et[idx[k]] <- TRUE
          }
          ch <- ch[keep_et, , drop = FALSE]
        }
        # Stage 3: continuity gate (drop dark-cycle-heavy tracks) — same score as live.
        if (q_on && is.finite(g_cycletime) && g_cycletime > 0 && nrow(ch) > 0) {
          sc <- .continuity_scores(ch, g_cycletime, q_H)
          if (!is.null(sc)) {
            keep_q <- names(sc)[is.finite(sc) & sc >= q_thr]
            keep_q <- union(keep_q, setdiff(unique(ch$tid), names(sc)))
            ch <- ch[ch$tid %in% keep_q, , drop = FALSE]
          }
        }
        # Stage 4: metadata filters (AFTER cleaning, matching the live order).
        for (nm in c("efo","cfr","dcr","tim","eco","dwell","fluo")) {
          rng_m <- isolate(input[[paste0(nm, "_rng")]])
          if (!is.null(rng_m) && nm %in% names(ch) && nrow(ch) > 0) {
            v <- ch[[nm]]; v <- v[is.finite(v)]
            if (length(v) && !(rng_m[1] > max(v) || rng_m[2] < min(v))) {
              ch <- ch[is.na(ch[[nm]]) | (ch[[nm]] >= rng_m[1] & ch[[nm]] <= rng_m[2]), , drop = FALSE]
            }
          }
        }
        # Stage 5: TID length filter.
        rng_len <- isolate(input$trk_tid_len)
        if (!is.null(rng_len) && "tid" %in% names(ch) && nrow(ch) > 0) {
          tl <- table(ch$tid); tl_v <- as.numeric(tl)
          if (!(rng_len[1] <= min(tl_v) && rng_len[2] >= max(tl_v))) {  # apply unless it covers all
            keep <- names(tl)[tl_v >= rng_len[1] & tl_v <= rng_len[2]]
            ch <- ch[ch$tid %in% keep, , drop = FALSE]
          }
        }
        if (nrow(ch) > 0) all_dfs[[i]] <- ch
      }

      df_all <- dplyr::bind_rows(Filter(Negate(is.null), all_dfs))
      if (nrow(df_all) == 0) {
        showNotification("No data survived the filters  -  rechunk aborted.", type = "error", duration = 8)
        return(invisible(NULL))
      }

      # Duration filter on combined data
      rng_dur <- isolate(input$trk_dur_filter)
      if (!is.null(rng_dur)) {
        tdur <- df_all %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(duration_ms = max(tim) - min(tim), .groups = "drop")
        dv <- tdur$duration_ms
        if (length(dv) > 0 && !(rng_dur[1] <= min(dv) && rng_dur[2] >= max(dv))) {  # apply unless covers all
          keep_tids <- tdur$tid[dv >= rng_dur[1] & dv <= rng_dur[2]]
          df_all <- df_all[df_all$tid %in% keep_tids, , drop = FALSE]
        }
      }
      if (nrow(df_all) == 0) {
        showNotification("No TIDs survived the duration filter  -  rechunk aborted.", type = "error", duration = 8)
        return(invisible(NULL))
      }

      # Spatial filter — apply if active, using the same pct_thr threshold
      sf <- isolate(trk_spatial_filter())
      if (!is.null(sf) && "x" %in% names(df_all) && "y" %in% names(df_all)) {
        in_box   <- df_all$x >= sf$x[1] & df_all$x <= sf$x[2] &
                    df_all$y >= sf$y[1] & df_all$y <= sf$y[2]
        all_tids <- unique(df_all$tid)
        keep_sp  <- vapply(all_tids, function(tid) {
          idx <- df_all$tid == tid
          mean(in_box[idx], na.rm = TRUE) >= sf$pct_thr
        }, logical(1))
        df_all <- df_all[df_all$tid %in% all_tids[keep_sp], , drop = FALSE]
        if (nrow(df_all) == 0) {
          showNotification("No TIDs survived the spatial filter  -  rechunk aborted.", type = "error", duration = 8)
          return(invisible(NULL))
        }
      }

      # Convert tim back to seconds (trk_raw_df will multiply by 1000 on load)
      df_all$tim <- df_all$tim / 1000

      # Compute new tid_sum and sort longest-first
      tid_sum_new <- df_all %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          duration_ms = (max(tim) - min(tim)) * 1000,
          n_pts = dplyr::n(), .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(duration_ms))
      df_all <- df_all[order(match(df_all$tid, tid_sum_new$tid)), , drop = FALSE]

      incProgress(0.6, detail = "Writing new chunks...")
      out_dir  <- file.path(tempdir(), paste0("rechunked_", as.integer(Sys.time())))
      user_tpc <- suppressWarnings(as.integer(isolate(input$trk_tids_per_chunk)))
      if (!is.finite(user_tpc) || user_tpc < 50) user_tpc <- TIDS_PER_CHUNK_DEFAULT

      manf_new <- split_by_tid_to_dir(df_all, out_dir, tids_per_chunk = user_tpc, fmt = CHUNK_FORMAT)
      if (!file.exists(manf_new)) {
        showNotification("Rechunking failed (split_by_tid_to_dir returned no manifest).", type = "error", duration = 8)
        return(invisible(NULL))
      }

      readr::write_csv(tid_sum_new, file.path(out_dir, "tid_summary.csv"))
      rv$manifest_dir <- out_dir
      trk_data_key(as.character(as.numeric(Sys.time())))

      # Compute stats across ALL remaining TIDs so filter UIs auto-update to the
      # new data extent (not just the first chunk). tim is currently in seconds;
      # multiply by 1000 to convert to ms format used in the filter UIs and stats displays.
      loc_vec_r <- function(col) if (col %in% names(df_all)) df_all[[col]][is.finite(df_all[[col]])] else numeric(0)
      trk_rechunk_stats_rv(list(
        tid_length   = as.numeric(table(df_all$tid)),
        tid_duration = tid_sum_new$duration_ms,
        efo          = loc_vec_r("efo"),
        cfr          = loc_vec_r("cfr"),
        dcr          = loc_vec_r("dcr"),
        tim          = loc_vec_r("tim") * 1000,
        eco          = loc_vec_r("eco"),
        dwell        = loc_vec_r("dwell"),
        fluo         = loc_vec_r("fluo")
      ))
      # Lock plot axes to the full spatial extent of all remaining TIDs so the
      # view stays stable as individual TIDs are played through chunks.
      rx <- df_all$x[is.finite(df_all$x)]; ry <- df_all$y[is.finite(df_all$y)]
      if (length(rx) && length(ry))
        trk_rechunk_bbox(list(x = range(rx), y = range(ry)))

      # The split (if any) is now BAKED into the saved segment TIDs. Live splitting
      # can stay ON: trk_gap_split strips the "_s<N>" suffix before splitting, so the
      # reloaded segments are recombined and re-split to the SAME segments rather than
      # split a second time. (That double split is what used to shrink the duration
      # upper bound on every rechunk.) Pressing Filter & Rechunk again is idempotent.
      # Flags stay TRUE so the length/duration sliders snap to the new full range
      # on their next (data-driven) re-render — no forced re-render here, which
      # would fire during the async checkbox round-trip while the split is still on.
      trk_len_reset_flag(TRUE)
      reset_dur_filter_flag(TRUE)

      n_tids <- length(unique(df_all$tid))
      n_new  <- nrow(readr::read_csv(manf_new, show_col_types = FALSE))
      showNotification(
        sprintf("Rechunked: %d TIDs in %d chunk(s).", n_tids, n_new),
        type = "message", duration = 8
      )
    })
  }, ignoreInit = TRUE)

  # SMLM-style metadata histogram: raw (grey) + filtered (blue) overlay
  output$trk_meta_hist <- renderPlotly({
    param <- input$trk_hist_param
    req(!is.null(param), param != "none")
    n_bins   <- max(5L, as.integer(input$trk_hist_bins %||% 60L))
    show_raw <- isTRUE(input$trk_hist_show_raw)
    zoom     <- input$trk_hist_range  # may be NULL before UI renders

    # Grey bars: frozen at upload time, unaffected by filtering or rechunking
    raw_stats <- trk_raw_stats_rv()
    raw_vec   <- if (!is.null(raw_stats) && param %in% names(raw_stats))
                   raw_stats[[param]] else numeric(0)

    # Blue bars: fully-filtered data — debounced so the histogram only recomputes
    # after the user stops moving sliders, not on every tick.
    filt_df  <- tryCatch(trk_dur_filtered_d(), error = function(e) NULL)
    filt_vec <- if (param == "tid_length") {
      if (!is.null(filt_df) && nrow(filt_df) > 0 && "tid" %in% names(filt_df))
        as.numeric(table(filt_df$tid))
      else numeric(0)
    } else if (param == "tid_duration") {
      if (!is.null(filt_df) && nrow(filt_df) > 0 && all(c("tid","tim") %in% names(filt_df)))
        filt_df %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(dur = max(tim) - min(tim), .groups = "drop") %>%
          dplyr::pull(dur)
      else numeric(0)
    } else {
      if (!is.null(filt_df) && nrow(filt_df) > 0 && param %in% names(filt_df))
        filt_df[[param]][is.finite(filt_df[[param]])]
      else numeric(0)
    }

    # Subset both vectors to the zoom range (same approach as SMLM tab)
    if (!is.null(zoom) && length(zoom) == 2 && all(is.finite(zoom)) && zoom[1] < zoom[2]) {
      raw_vec  <- raw_vec[raw_vec  >= zoom[1] & raw_vec  <= zoom[2]]
      filt_vec <- filt_vec[filt_vec >= zoom[1] & filt_vec <= zoom[2]]
    }
    shiny::validate(shiny::need(
      length(filt_vec) > 0 || (show_raw && length(raw_vec) > 0),
      "No data to plot for the current filters/zoom range."
    ))

    xlab <- switch(param,
      tid_length   = "Trace ID length (localizations)",
      tid_duration = "Trace ID duration (ms)",
      toupper(param)
    )
    ylab <- if (param %in% c("tid_length", "tid_duration")) "Number of traces" else "Count"

    fig <- plot_ly()

    if (show_raw && length(raw_vec) > 0) {
      fig <- fig %>% add_histogram(
        x        = raw_vec,
        nbinsx   = n_bins,
        bingroup = "x",
        name     = "Raw data",
        opacity  = 0.35,
        marker   = list(color = "lightgray")
      )
    }
    if (length(filt_vec) > 0) {
      fig <- fig %>% add_histogram(
        x        = filt_vec,
        nbinsx   = n_bins,
        bingroup = "x",
        name     = "Filtered data",
        opacity  = 0.85,
        marker   = list(color = "steelblue")
      )
    }
    fig %>% layout(
      barmode = "overlay",
      xaxis   = list(title = xlab),
      yaxis   = list(title = ylab),
      legend  = list(orientation = "h", x = 0, y = 1.1)
    )
  })

  # ---- Background image overlay (tracking tab) --------------------------------
  # All IDs are trk_-prefixed; logic mirrors the SMLM tab exactly.

  trk_confocal_meta <- reactive({
    file <- input$trk_bg_image
    if (is.null(file) || !isTRUE(input$trk_show_bg)) return(NULL)
    img <- tryCatch(magick::image_read(file$datapath), error = function(e) NULL)
    if (is.null(img)) return(NULL)
    rot <- suppressWarnings(as.integer(input$trk_img_rotation %||% 0))
    if (is.finite(rot) && (rot %% 360) != 0) img <- magick::image_rotate(img, rot)
    if (isTRUE(input$trk_img_flip_x)) img <- magick::image_flop(img)
    img     <- apply_background_lut(img, input$trk_img_lut %||% "original",
                                    custom_color = input$trk_img_lut_custom,
                                    brightness   = input$trk_img_brightness %||% 100,
                                    contrast     = input$trk_img_contrast   %||% 0)
    info    <- magick::image_info(img)
    raw_png <- magick::image_write(img, format = "png")
    list(
      uri       = paste0("data:image/png;base64,", base64enc::base64encode(raw_png)),
      width_px  = info$width,
      height_px = info$height
    )
  })

  trk_image_overlay <- reactive({
    meta <- trk_confocal_meta()
    if (is.null(meta)) return(NULL)
    ps <- suppressWarnings(as.numeric(input$trk_px_size))
    if (!is.finite(ps) || ps <= 0) return(NULL)
    x0 <- suppressWarnings(as.numeric(input$trk_img_x0)); if (!is.finite(x0)) x0 <- 0
    y0 <- suppressWarnings(as.numeric(input$trk_img_y0)); if (!is.finite(y0)) y0 <- 0
    sizex <- meta$width_px  * ps
    sizey <- meta$height_px * ps
    if (isTRUE(input$trk_img_invert_y)) {
      y_anchor <- "top";    y_val <- y0 + sizey
    } else {
      y_anchor <- "bottom"; y_val <- y0
    }
    list(
      source  = meta$uri,
      xref = "x", yref = "y",
      x = x0,     y = y_val,
      sizex = sizex, sizey = sizey,
      xanchor = "left", yanchor = y_anchor,
      sizing  = "stretch",
      opacity = if (is.null(input$trk_img_opacity)) 0.6 else input$trk_img_opacity,
      layer   = "below"
    )
  })

  # ---- Filter settings save/load ----
  output$trk_download_filter_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$trk_file$name)) tools::file_path_sans_ext(input$trk_file$name[1]) else "tracking"
      paste0(base, "_filter_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content = function(file) {
      s <- list(
        efo_rng   = input$efo_rng,   cfr_rng  = input$cfr_rng,
        dcr_rng   = input$dcr_rng,   tim_rng  = input$tim_rng,
        eco_rng   = input$eco_rng,   dwell_rng = input$dwell_rng,
        fluo_rng  = input$fluo_rng,
        trk_tid_len    = input$trk_tid_len,
        trk_dur_filter = input$trk_dur_filter
      )
      writeLines(jsonlite::toJSON(s, pretty = TRUE, auto_unbox = TRUE, na = "null"), file)
    }
  )
  observeEvent(input$trk_load_filter_settings, {
    req(input$trk_load_filter_settings$datapath)
    s <- jsonlite::fromJSON(input$trk_load_filter_settings$datapath)
    upd_rng <- function(id) {
      v <- s[[id]]
      if (!is.null(v) && length(v) == 2)
        shinyWidgets::updateNumericRangeInput(session, id, value = v)
    }
    for (id in c("efo_rng","cfr_rng","dcr_rng","tim_rng","eco_rng","dwell_rng","fluo_rng"))
      upd_rng(id)
    if (!is.null(s$trk_tid_len) && length(s$trk_tid_len) == 2)
      shinyWidgets::updateNumericRangeInput(session, "trk_tid_len", value = s$trk_tid_len)
    if (!is.null(s$trk_dur_filter) && length(s$trk_dur_filter) == 2)
      shinyWidgets::updateNumericRangeInput(session, "trk_dur_filter", value = s$trk_dur_filter)
    showNotification("Filter settings loaded.", type = "message", duration = 3)
  })

  # ---- Playback settings save/load ----
  output$trk_download_playback_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$trk_file$name)) tools::file_path_sans_ext(input$trk_file$name[1]) else "tracking"
      paste0(base, "_playback_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content = function(file) {
      s <- list(
        trk_play_mode      = input$trk_play_mode,
        trk_play_total     = input$trk_play_total,
        trk_bin_dt         = input$trk_bin_dt,
        trk_play_speed     = input$trk_play_speed,
        trk_target_markers = input$trk_target_markers,
        trk_trail_len      = input$trk_trail_len,
        trk_ms_per_unit    = input$trk_ms_per_unit,
        trk_per_tid_cap    = input$trk_per_tid_cap,
        trk_trail_nmax     = input$trk_trail_nmax,
        trk_max_pts_tid    = input$trk_max_pts_tid
      )
      writeLines(jsonlite::toJSON(s, pretty = TRUE, auto_unbox = TRUE, na = "null"), file)
    }
  )
  observeEvent(input$trk_load_playback_settings, {
    req(input$trk_load_playback_settings$datapath)
    s <- jsonlite::fromJSON(input$trk_load_playback_settings$datapath)
    upd_num <- function(id) {
      v <- s[[id]]
      if (!is.null(v) && is.finite(as.numeric(v)))
        updateNumericInput(session, id, value = as.numeric(v))
    }
    if (!is.null(s$trk_play_mode))
      updateSelectInput(session, "trk_play_mode", selected = s$trk_play_mode)
    if (!is.null(s$trk_play_speed) && is.finite(as.numeric(s$trk_play_speed)))
      updateSliderInput(session, "trk_play_speed", value = as.numeric(s$trk_play_speed))
    for (id in c("trk_play_total","trk_bin_dt","trk_target_markers",
                 "trk_trail_len","trk_ms_per_unit",
                 "trk_per_tid_cap","trk_trail_nmax","trk_max_pts_tid"))
      upd_num(id)
    showNotification("Playback settings loaded.", type = "message", duration = 3)
  })

  output$trk_download_overlay_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$trk_file$name)) tools::file_path_sans_ext(input$trk_file$name[1]) else "tracking"
      paste0(base, "_image_alignment_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content  = function(file) {
      ov <- list(
        show_confocal = input$trk_show_confocal,
        show_bg       = input$trk_show_bg,
        bg_image_name = if (!is.null(input$trk_bg_image)) input$trk_bg_image$name else NULL,
        px_size       = input$trk_px_size,
        img_x0        = input$trk_img_x0,
        img_y0        = input$trk_img_y0,
        img_invert_y  = input$trk_img_invert_y,
        img_flip_x    = input$trk_img_flip_x,
        img_opacity   = input$trk_img_opacity,
        img_rotation  = input$trk_img_rotation,
        img_lut       = input$trk_img_lut,
        img_lut_custom = input$trk_img_lut_custom,
        img_brightness = input$trk_img_brightness,
        img_contrast   = input$trk_img_contrast
      )
      writeLines(jsonlite::toJSON(ov, pretty = TRUE, auto_unbox = TRUE, na = "null"), file)
    }
  )

  observeEvent(input$trk_load_overlay_settings, {
    req(input$trk_load_overlay_settings)
    ov <- tryCatch(jsonlite::fromJSON(input$trk_load_overlay_settings$datapath),
                   error = function(e) NULL)
    if (is.null(ov)) return()
    if (!is.null(ov$show_confocal)) updateCheckboxInput(session, "trk_show_confocal", value = ov$show_confocal)
    if (!is.null(ov$show_bg))       updateCheckboxInput(session, "trk_show_bg",       value = ov$show_bg)
    if (!is.null(ov$px_size))       updateNumericInput(session,  "trk_px_size",        value = ov$px_size)
    if (!is.null(ov$img_x0))        updateNumericInput(session,  "trk_img_x0",         value = ov$img_x0)
    if (!is.null(ov$img_y0))        updateNumericInput(session,  "trk_img_y0",         value = ov$img_y0)
    if (!is.null(ov$img_invert_y))  updateCheckboxInput(session, "trk_img_invert_y",  value = ov$img_invert_y)
    if (!is.null(ov$img_flip_x))    updateCheckboxInput(session, "trk_img_flip_x",    value = ov$img_flip_x)
    if (!is.null(ov$img_opacity))   updateSliderInput(session,   "trk_img_opacity",   value = ov$img_opacity)
    if (!is.null(ov$img_rotation))  updateSelectInput(session,   "trk_img_rotation",  selected = as.character(ov$img_rotation))
    if (!is.null(ov$img_lut))       updateSelectInput(session,   "trk_img_lut",       selected = as.character(ov$img_lut))
    if (!is.null(ov$img_lut_custom)) colourpicker::updateColourInput(session, "trk_img_lut_custom", value = as.character(ov$img_lut_custom))
    if (!is.null(ov$img_brightness)) updateSliderInput(session,   "trk_img_brightness", value = ov$img_brightness)
    if (!is.null(ov$img_contrast))   updateSliderInput(session,   "trk_img_contrast",   value = ov$img_contrast)
    if (!is.null(ov$bg_image_name)) {
      cur <- if (!is.null(input$trk_bg_image)) input$trk_bg_image$name else NULL
      if (is.null(cur) || !identical(cur, ov$bg_image_name))
        showNotification(paste0("Re-upload '", ov$bg_image_name, "' to reproduce alignment."),
                         type = "message", duration = 6)
    }
  })

  output$trk_download_filtered_csv <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$trk_file$name))
        tools::file_path_sans_ext(input$trk_file$name[1]) else "tracking"
      paste0(base, "_filtered_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      df <- trk_spatial_filtered()
      shiny::validate(shiny::need(nrow(df) > 0, "No filtered data to export."))
      drop_cols <- intersect(c("col", "t"), names(df))
      if (length(drop_cols)) df <- df[, !names(df) %in% drop_cols, drop = FALSE]
      # tim is stored internally in ms; convert back to seconds for export
      if ("tim" %in% names(df)) df$tim <- df$tim / 1000
      write.csv(df, file, row.names = FALSE)
    }
  )
  observeEvent(input$trk_download_filtered_csv, {
    showNotification("Exported filtered tracking data.", type = "message", duration = 3)
  })

  observeEvent(input$trk_reset_overlay, {
    updateCheckboxInput(session, "trk_show_confocal", value = FALSE)
    updateCheckboxInput(session, "trk_show_bg",       value = FALSE)
    updateNumericInput(session,  "trk_px_size",        value = 10)
    updateNumericInput(session,  "trk_img_x0",         value = 0)
    updateNumericInput(session,  "trk_img_y0",         value = 0)
    updateCheckboxInput(session, "trk_img_invert_y",  value = TRUE)
    updateCheckboxInput(session, "trk_img_flip_x",    value = FALSE)
    updateSliderInput(session,   "trk_img_opacity",   value = 0.6)
    updateSelectInput(session,   "trk_img_rotation",  selected = "0")
    updateSelectInput(session,   "trk_img_lut",       selected = "original")
    colourpicker::updateColourInput(session, "trk_img_lut_custom", value = "#00FF00")
    updateSliderInput(session,   "trk_img_brightness", value = 100)
    updateSliderInput(session,   "trk_img_contrast",   value = 0)
    shinyjs::reset("trk_overlay_form")
  })

  # Overlay changes: push via proxy relayout (no DOM/frames rebuild needed).
  observeEvent(
    list(input$trk_show_bg, input$trk_bg_image, input$trk_px_size,
         input$trk_img_x0, input$trk_img_y0, input$trk_img_invert_y,
         input$trk_img_flip_x, input$trk_img_opacity, input$trk_img_rotation,
         input$trk_img_lut, input$trk_img_lut_custom,
         input$trk_img_brightness, input$trk_img_contrast),
    {
      req(!is.na(layout_trigger()))   # plot must already exist
      ovl <- tryCatch(isolate(trk_image_overlay()), error = function(e) NULL)
      plotlyProxyInvoke(
        plotlyProxy("trk_track2D", session),
        "relayout",
        list(images = if (!is.null(ovl)) list(ovl) else list())
      )
    },
    ignoreInit = TRUE
  )

  # ---- Auto-align overlay helpers (shared image-loading logic) ----
  .load_trk_img <- function() {
    f <- input$trk_bg_image
    if (is.null(f)) return(NULL)
    img <- tryCatch(magick::image_read(f$datapath), error = function(e) NULL)
    if (is.null(img)) return(NULL)
    rot <- suppressWarnings(as.integer(input$trk_img_rotation %||% 0))
    if (is.finite(rot) && (rot %% 360) != 0) img <- magick::image_rotate(img, rot)
    if (isTRUE(input$trk_img_flip_x)) img <- magick::image_flop(img)
    img
  }

  observeEvent(input$trk_auto_foci, {
    img <- .load_trk_img(); req(!is.null(img))
    df  <- isolate(trk()); req(nrow(df) > 0)
    ps  <- suppressWarnings(as.numeric(input$trk_px_size))
    req(is.finite(ps), ps > 0)
    res <- tryCatch(
      auto_align_overlay(img, df, ps, mode = "foci",
                         invert_y    = isTRUE(input$trk_img_invert_y),
                         n_peaks        = as.integer(input$trk_foci_n_peaks        %||% 5L),
                         bright_pct     = (input$trk_foci_bright_pct     %||% 95)  / 100,
                         border_pct     = (input$trk_foci_border_pct     %||% 10)  / 100,
                         eps_pct        = (input$trk_foci_eps_pct        %||% 3)   / 100,
                         aspect_max     =  input$trk_foci_aspect_max     %||% 4.0,
                         density_thr    = (input$trk_foci_density_thr    %||% 20)  / 100,
                         density_cell   =  input$trk_foci_density_cell   %||% 100,
                         density_smooth =  input$trk_foci_density_smooth %||% 100),
      error = function(e) {
        showNotification(paste("Foci alignment error:", conditionMessage(e)),
                         type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(res)) {
      showNotification("Foci alignment: could not detect matching peaks. Check image and pixel size.",
                       type = "warning", duration = 6)
      return()
    }
    # Store diagnostic peaks (img_peaks are in nm assuming x0=0; add x0/y0 for plot coords)
    n_blobs    <- if (!is.null(res$img_peaks))  nrow(res$img_peaks)  else 0L
    n_clusters <- if (!is.null(res$data_peaks)) nrow(res$data_peaks) else 0L
    if (n_blobs > 0L) {
      rv$trk_diag <- list(
        data_clusters = res$data_peaks,
        img_blobs     = data.frame(x = res$x0 + res$img_peaks$x,
                                   y = res$y0 + res$img_peaks$y)
      )
    }
    updateCheckboxInput(session, "trk_show_bg", value = TRUE)
    updateNumericInput(session, "trk_img_x0", value = res$x0)
    updateNumericInput(session, "trk_img_y0", value = res$y0)
    if (!is.na(layout_trigger())) fire_layout("foci_align")
    showNotification(
      paste0("Foci alignment applied (", n_blobs, " confocal blobs, ",
             n_clusters, " MINFLUX clusters detected). ",
             "Enable 'Show detection diagnostics' to visualise."),
      type = "message", duration = 6)
  }, ignoreInit = TRUE)

  observeEvent(input$trk_show_diag, {
    if (!is.na(layout_trigger())) fire_layout("show_diag")
  }, ignoreInit = TRUE)

  observeEvent(input$trk_auto_boundary, {
    img <- .load_trk_img(); req(!is.null(img))
    df  <- isolate(trk()); req(nrow(df) > 0)
    ps  <- suppressWarnings(as.numeric(input$trk_px_size))
    req(is.finite(ps), ps > 0)
    res <- tryCatch(
      auto_align_overlay(img, df, ps, mode = "boundary",
                         invert_y = isTRUE(input$trk_img_invert_y)),
      error = function(e) {
        showNotification(paste("Boundary alignment error:", conditionMessage(e)),
                         type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(res)) {
      showNotification("Boundary alignment: could not compute centroid. Check image and pixel size.",
                       type = "warning", duration = 6)
      return()
    }
    updateCheckboxInput(session, "trk_show_bg", value = TRUE)
    updateNumericInput(session, "trk_img_x0", value = res$x0)
    updateNumericInput(session, "trk_img_y0", value = res$y0)
    showNotification("Boundary alignment applied — nudge to fine-tune.",
                     type = "message", duration = 4)
  }, ignoreInit = TRUE)

  # ---- Nudge buttons for X and Y origin ----
  local({
    deltas <- c(m1000 = -1000, m100 = -100, m10 = -10,
                p10   =    10, p100 =  100, p1000 = 1000)
    for (sfx in names(deltas)) {
      local({
        s <- sfx          # bind into this env so the reactive expression
        d <- deltas[[s]]  # resolves the correct string after loop ends
        observeEvent(input[[paste0("trk_nx_", s)]], {
          v <- suppressWarnings(as.numeric(input$trk_img_x0))
          updateNumericInput(session, "trk_img_x0",
                             value = (if (is.finite(v)) v else 0) + d)
        }, ignoreInit = TRUE)
        observeEvent(input[[paste0("trk_ny_", s)]], {
          v <- suppressWarnings(as.numeric(input$trk_img_y0))
          updateNumericInput(session, "trk_img_y0",
                             value = (if (is.finite(v)) v else 0) + d)
        }, ignoreInit = TRUE)
      })
    }
  })


  output$t_window_ui <- renderUI({
    req(!is.na(layout_trigger()))
    df <- trk(); req(nrow(df) > 0)
    tcol <- if (isTRUE(is_simul())) "tim0" else "tim"
    v <- df[[tcol]]; v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    rng <- range(v)
    numericRangeInput(
      "t_window",
      paste0("Time window on ", tcol, " [", signif(rng[1], 4), ", ", signif(rng[2], 4), "]"),
      value = rng, min = rng[1], max = rng[2]
    )
  })
  
  
  # ---- Two-level trigger system -----------------------------------------------
  # layout_trigger: fires only when the Plotly CHART STRUCTURE changes (number of
  #   traces, colors, static lines, heatmap layer, play mode).  renderPlotly
  #   depends on this  -  reconstructing the full WebGL chart is expensive and
  #   blanks the display, so it must happen as rarely as possible.
  #
  # plot_trigger: fires whenever FRAME DATA needs rebuilding (bin_dt, n_frames,
  #   trail params, marker quotas, auto-tune).  frames_store/trk/seq_data depend
  #   on this.  renderPlotly does NOT depend on plot_trigger; after a frames
  #   rebuild the new frame 1 is pushed to the chart via draw_frame() proxy.
  layout_trigger    <- reactiveVal(NA)
  plot_trigger      <- reactiveVal(NA)
  # TRUE once the track plot has actually rendered WITH event_register(), so the
  # relayout observer never reads event_data() before the "track" source is
  # registered (that prints a benign, deferred "source not registered" warning
  # suppressWarnings can't catch). Stays FALSE until a plot with data renders —
  # layout_trigger alone isn't enough (it can be set while frames_store is empty).
  trk_plot_ready    <- reactiveVal(FALSE)
  conf_color_map_stored    <- reactiveVal(NULL)  # pre-computed tid->hex when conf coloring is on
  cluster_color_map_stored <- reactiveVal(NULL)  # pre-computed tid->hex when MSD cluster coloring is on

  # TRUE when the loaded data has a usable z column (enables 3D toggle)
  has_z <- reactive({
    df <- trk()
    "z" %in% names(df) && sum(is.finite(df$z)) > 0
  })

  # TRUE when the user has selected 3D mode and data supports it
  use_3d <- reactive({
    isTRUE(has_z()) && identical(input$trk_view_mode %||% "2d", "3d")
  })
  trk_axis_x_span   <- reactiveVal(NULL)
  trk_locked_ranges <- reactiveVal(NULL)   # list(x = c(x0,x1), y = c(y0,y1)) — held fixed during playback
  # Stable key that only changes on new upload/rechunk — used as uirevision so
  # display-only re-renders (show/hide traces, dark mode) do NOT reset user zoom.
  trk_data_key      <- reactiveVal("init")
  trk_default_axis_state <- function() {
    pre <- isolate(frames_store())
    if (!is.list(pre) || !length(pre$frames)) return(NULL)

    f1 <- pre$frames[[1]]
    bx <- f1$markers_x[is.finite(f1$markers_x)]
    by <- f1$markers_y[is.finite(f1$markers_y)]
    if (!length(bx) || !length(by)) {
      sp <- pre$sp
      bx <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
      by <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
      bx <- bx[is.finite(bx)]
      by <- by[is.finite(by)]
    }
    if (!length(bx) || !length(by)) return(NULL)

    xr   <- range(bx)
    yr   <- range(by)
    xpad <- max(1, diff(xr) * 0.04)
    ypad <- max(1, diff(yr) * 0.04)
    list(
      xr = xr,
      yr = yr,
      xpad = xpad,
      ypad = ypad,
      x_span = diff(xr) + 2 * xpad
    )
  }

  # data_ver: bumped whenever the filtered dataset itself changes (filter inputs OR new
  # upload/layout).  trk() is eventReactive on this so prep_tracking + tid_color_map
  # are skipped on rebuild_tick-only (playback param) plot_trigger fires.
  data_ver <- reactiveVal(0L)

  # src is a human-readable label of the trigger (kept for self-documenting call sites).
  fire_layout <- function(src = "?") {
    t <- Sys.time()
    loop$last_fire_layout <- TRUE
    layout_trigger(t)
    plot_trigger(t)
    data_ver(isolate(data_ver()) + 1L)   # layout = new data shape; rerun trk()
  }
  fire_frames <- function(src = "?") {
    loop$last_fire_layout <- FALSE
    plot_trigger(Sys.time())
  }

  # Plot Data: full mode re-renders the static plot (no frames); playback proxy-pushes.
  # Plot Data = a full, reliable rebuild in BOTH modes (re-reads the filtered data
  # and re-renders the chart, like Filter & Rechunk). fire_layout bumps data_ver
  # (reruns trk()) and rebuilds the WebGL chart from scratch, so the plot always
  # reflects the current filters when pressed.
  observeEvent(input$trk_plot_data, {
    bump_clean()      # (re)apply split / end-trim / continuity gate with current settings
    fire_layout("plot_data_btn")
  })
  observeEvent(list(
    input$trk_play_mode,
    input$trk_view_mode_main, input$trk_full_connect, input$trk_full_raw,
    input$trk_show_full_tracks_faint, input$trk_show_moving_trails
    # trk_show_heatmap deliberately NOT here — the dwell-time heatmap is heavy, so it
    # applies on the next Plot Data (heatmap_store is plot_trigger-gated). Toggling the
    # checkbox alone no longer forces a full rebuild that showed nothing until Plot Data.
  ), {
    if (!is.na(layout_trigger()) && !isTRUE(is_playing())) fire_layout("display_opts")
  }, ignoreInit = TRUE)

  # Frame-data-only changes -> frames rebuild, chart stays
  rebuild_tick <- reactiveVal(0L)
  observeEvent(list(
    input$trk_bin_dt,
    input$trk_target_markers, input$trk_per_tid_cap,
    input$trk_trail_len, input$trk_trail_nmax,
    input$trk_max_pts_tid, input$trk_play_total
  ), {
    if (!isTRUE(is_playing()) && !is.na(plot_trigger()) && !isTRUE(loop$auto_tuning)) {
      rebuild_tick(isolate(rebuild_tick()) + 1L)
    }
  }, ignoreInit = TRUE)

  observeEvent(debounce(rebuild_tick(), 200), {
    fire_frames()
  }, ignoreInit = TRUE)

  # Filter changes do NOT live-update the track plot (the old filter_tick auto-replot
  # was removed). The metadata / length / duration filters — like the split / trim /
  # quality-gate cleaning stages — apply to the plotted data ONLY when the user presses
  # Plot Data (fire_layout → bumps data_ver → reruns trk()) or Filter & Rechunk. Live,
  # filter changes still drive the metadata HISTOGRAM (reads the debounced filtered
  # stages directly) and the length/duration range read-outs (re-keyed on data_ver, so
  # they describe the plotted data) — but the WebGL track render stays put until an
  # explicit Plot Data / Rechunk. See observeEvent(input$trk_plot_data) above.


  # ---- Prepare data once per filter-data change ----
  # trk() is keyed on data_ver, NOT plot_trigger.  This means it skips
  # prep_tracking + tid_color_map for rebuild_tick-only (playback param) triggers
  # where the filtered dataframe is identical to the previous run.
  # frames_store() still runs on every plot_trigger; it reads isolate(trk()) and
  # gets the cached result when data_ver did not change.
  trk <- eventReactive(data_ver(), {
    df <- isolate(trk_spatial_filtered()); req(nrow(df) > 0)
    df$tid <- as.character(df$tid)
    df <- prep_tracking(df, normalize_per_tid = isolate(isTRUE(is_simul())))
    cmap <- tid_color_map(df$tid)

    # For split-track segments (TIDs matching "<parent>_s<N>"), derive a darker
    # shade of the first segment's colour so sibling segments share a hue family.
    # _s1 keeps the base colour; _s2 is 80% brightness, _s3 is 65%, etc.
    split_ids <- grep("_s\\d+$", names(cmap), value = TRUE, perl = TRUE)
    if (length(split_ids)) {
      parents <- unique(sub("_s\\d+$", "", split_ids, perl = TRUE))
      for (par in parents) {
        segs <- sort(as.integer(sub(paste0("^", par, "_s"), "", grep(
          paste0("^", par, "_s\\d+$"), names(cmap), value = TRUE, perl = TRUE
        ), perl = TRUE)))
        if (!length(segs)) next
        base_nm  <- paste0(par, "_s", segs[1])
        base_hex <- cmap[[base_nm]]
        for (sn in segs[-1]) {
          nm_sn <- paste0(par, "_s", sn)
          if (!nm_sn %in% names(cmap)) next
          rgb_m <- grDevices::col2rgb(base_hex) / 255
          hsv_m <- grDevices::rgb2hsv(rgb_m[1], rgb_m[2], rgb_m[3])
          v_new <- max(0.55, hsv_m["v", 1] * (1 - 0.09 * (sn - 1L)))
          cmap[[nm_sn]] <- grDevices::hsv(hsv_m["h", 1], hsv_m["s", 1], v_new)
        }
      }
    }

    df$col <- unname(cmap[as.character(df$tid)])
    attr(df, "col_map") <- cmap
    df
  }, ignoreNULL = TRUE)

  # Sequential data (always order by raw 'tim'; decimate for speed)
  # ── Chunk-scoped PLAYBACK ───────────────────────────────────────────────────
  # Chunking exists only to keep the animation smooth: the filter pipeline,
  # analysis, and the full-track static overlay all use ALL TIDs, but the moving
  # markers/trails are limited to the selected chunk's TIDs. Returns the chunk's
  # ORIGINAL TIDs (or NULL = animate everything when "All TIDs" is selected).
  trk_playback_chunk_tids <- reactive({
    sel <- input$trk_chunk_sel
    if (!isTruthy(sel) || identical(sel, "__all__") || !file.exists(sel)) return(NULL)
    ch <- tryCatch(
      if (grepl("\\.qs$",  sel, TRUE)) qs::qread(sel)
      else if (grepl("\\.fst$", sel, TRUE)) fst::read_fst(sel, columns = "tid") |> as.data.frame()
      else NULL,
      error = function(e) NULL
    )
    if (is.null(ch) || !"tid" %in% names(ch) || !nrow(ch)) return(NULL)
    unique(as.character(ch$tid))
  })

  # Subset a display data.frame to the chunk's TIDs. Split segments ("42_s1") are
  # matched by their base TID ("42") so a chunk still animates its split children.
  # Falls back to the full frame if the scope would leave nothing to animate.
  scope_playback <- function(d, chunk_tids) {
    if (is.null(chunk_tids) || !nrow(d) || !"tid" %in% names(d)) return(d)
    base <- sub("_s\\d+$", "", as.character(d$tid), perl = TRUE)
    keep <- base %in% chunk_tids
    if (!any(keep)) return(d)            # chunk fully filtered out — don't blank playback
    d[keep, , drop = FALSE]
  }

  seq_data <- eventReactive(plot_trigger(), {
    df <- isolate(trk()); req(nrow(df) > 0)
    df <- scope_playback(df, isolate(trk_playback_chunk_tids()))   # animate selected chunk only
    # decimate on raw 'tim' to keep timing shape but lighten the load
    df_d <- decimate_per_tid(df, max_pts = isolate(input$trk_max_pts_tid), tcol = "tim")
    df_d <- df_d[order(df_d$tid, df_d$t), ]  # 't' is raw tim here
    sp   <- split(df_d, df_d$tid)
    tids <- names(sp)
    # play TIDs in order of earliest time
    ord  <- order(sapply(sp, function(tb) min(tb$t, na.rm = TRUE)))
    list(sp = sp[tids[ord]], tids = tids[ord])
  }, ignoreNULL = TRUE)
  frames_store <- eventReactive(plot_trigger(), {
    df <- isolate(trk()); req(nrow(df) > 0)

    # Optional faint background of ALL raw (pre-filter) tracks — shown in BOTH full
    # and playback modes so you can see what the filters removed. Derived from the
    # POST-SPLIT raw object (trk_gap_split): when splitting is on the faint tracks
    # are split too (so they line up with the analyzed segments); when it's off,
    # gap_split returns the raw data unchanged. Computed once here (decimated), only
    # when the toggle is on.
    sp_raw <- if (isTRUE(isolate(input$trk_full_raw))) {
      mpb <- suppressWarnings(as.integer(isolate(input$trk_max_pts_tid)))
      if (!is.finite(mpb) || mpb < 50) mpb <- 400L
      rdf <- tryCatch(isolate(trk_gap_split()), error = function(e) NULL)
      if (!is.null(rdf) && nrow(rdf) && all(c("tid","tim","x","y") %in% names(rdf))) {
        rd <- decimate_per_tid(rdf, max_pts = mpb, tcol = "tim")
        split(rd, rd$tid)
      } else NULL
    } else NULL

    # ── FULL-TRACKS MODE: skip ALL playback frame/bin/trail machinery ──────────
    # The static render (view_mode_main == "full") reads only sp_static / col_map /
    # tids — never the per-frame data. Building frames here is what made length,
    # duration and split edits feel slow in the default view, so short-circuit:
    # compute just the full-track overlay + colours and return one empty placeholder
    # frame (keeps downstream req(length(frames) > 0) happy). Frames are built only
    # when the user switches to Playback, which fires fire_layout() → full rebuild.
    if (identical(isolate(input$trk_view_mode_main), "full")) {
      use3d_f <- isolate(isTRUE(use_3d()))
      max_pts <- suppressWarnings(as.integer(isolate(input$trk_max_pts_tid)))
      if (!is.finite(max_pts) || max_pts < 50) max_pts <- 400L

      # Progress feedback: building the per-track overlay (decimate + split into one
      # trace per track) is the visible cost of each Plot Data / upload-settle render.
      full_build <- withProgress(message = "Building track view…", value = 0.2, {
        df_static <- decimate_per_tid(df, max_pts = max_pts, tcol = "tim")
        setProgress(0.5, detail = sprintf("%s tracks", format(length(unique(df_static$tid)), big.mark = ",")))
        sp_static <- split(df_static, df_static$tid)
        if (use3d_f) {
          zs <- suppressWarnings(as.numeric(isolate(input$trk_z_stretch %||% 1.0)))
          if (is.finite(zs) && zs >= 0.1 && zs != 1.0) {
            zscale <- function(spl) lapply(spl, function(tb) { if ("z" %in% names(tb)) tb$z <- tb$z * zs; tb })
            sp_static <- zscale(sp_static); if (!is.null(sp_raw)) sp_raw <- zscale(sp_raw)
          }
        }
        tids_all     <- names(sp_static)
        base_col_map <- if (!is.null(attr(df, "col_map"))) attr(df, "col_map") else tid_color_map(tids_all)
        cluster_cm   <- isolate(cluster_color_map_stored())
        conf_cm      <- isolate(conf_color_map_stored())
        override_cm  <- if (!is.null(cluster_cm)) cluster_cm else conf_cm
        if (!is.null(override_cm))
          for (nm in names(override_cm)) if (nm %in% names(base_col_map)) base_col_map[[nm]] <- override_cm[[nm]]
        setProgress(1)
        list(sp_static = sp_static, tids_all = tids_all,
             base_col_map = base_col_map, sp_raw = sp_raw)
      })
      sp_static    <- full_build$sp_static
      tids_all     <- full_build$tids_all
      base_col_map <- full_build$base_col_map
      sp_raw       <- full_build$sp_raw

      return(list(
        static_lines = NULL,
        frames       = list(f1 = list(
          trails_per_tid = list(), trails_x = numeric(0), trails_y = numeric(0),
          markers_x = numeric(0), markers_y = numeric(0), markers_z = NULL,
          marker_cols = character(0), step_dt = NA_real_, bin_width = NA_real_
        )),
        t_grid     = NA_real_,
        col_map    = base_col_map,
        tids       = tids_all,
        sp         = sp_static,    # axis fallback / faint overlay read this
        sp_static  = sp_static,
        sp_raw     = sp_raw,       # optional faint background of ALL raw tracks
        trail_len  = 0, trail_nmax = 1L, bin_width = NA_real_,
        using_bins = FALSE, step_dt = NA_real_, has_3d = use3d_f
      ))
    }
    # ── PLAYBACK MODE: full windowing + per-frame binning below ────────────────

    # df         = ALL TIDs → used for the static full-track overlay (sp_static)
    # df_play    = selected chunk only → used for the animated markers/trails
    df_play <- scope_playback(df, isolate(trk_playback_chunk_tids()))
    use3d <- isolate(isTRUE(use_3d()))

    # Playback time basis
    tcol_play <- if (isTRUE(isolate(is_simul()))) "tim0" else "tim"

    # Fixed total playback window (in the chosen basis)
    total <- suppressWarnings(as.numeric(isolate(input$trk_play_total)))
    if (!is.finite(total) || total <= 0) total <- 1000
    
    # Determine [tmin, tmax] in the playback basis (over the animated chunk)
    rng_play_all <- range(df_play[[tcol_play]], na.rm = TRUE)
    if (isTRUE(isolate(is_simul()))) {
      tmin <- 0
      tmax <- min(tmin + total, max(rng_play_all, na.rm = TRUE))
    } else {
      tmin <- min(rng_play_all, na.rm = TRUE)
      tmax <- min(tmin + total, max(rng_play_all, na.rm = TRUE))
    }
    if (!is.finite(tmin) || !is.finite(tmax) || tmax <= tmin) {
      tmin <- rng_play_all[1]; tmax <- rng_play_all[2]
    }
    
    # Filter the ANIMATED data (chunk-scoped) to the fixed window
    dwin <- df_play[df_play[[tcol_play]] >= tmin & df_play[[tcol_play]] <= tmax, , drop = FALSE]
    req(nrow(dwin) > 0)
    
    # Light pre-decimation for TRAILS only
    max_pts_tid  <- suppressWarnings(as.integer(isolate(input$trk_max_pts_tid)))
    if (!is.finite(max_pts_tid) || max_pts_tid < 50) max_pts_tid <- 400
    dwin_tr_dec  <- decimate_per_tid(dwin, max_pts = max_pts_tid, tcol = tcol_play)
    dwin_tr_dec  <- dwin_tr_dec[order(dwin_tr_dec$tid, dwin_tr_dec$t), ]
    
    # Split full-res for markers, and decimated for trails
    sp_full <- split(dwin[order(dwin$tid, dwin[[tcol_play]]), ], dwin$tid)  # FULL (markers)
    sp_tr   <- split(dwin_tr_dec, dwin_tr_dec$tid)                          # DECIMATED (trails)

    # Full-range (unwindowed) spatial data for the static full-tracks trace.
    # Using df (not dwin) so tracks outside the playback window still appear.
    df_static    <- decimate_per_tid(df, max_pts = max_pts_tid, tcol = "tim")
    sp_static    <- split(df_static, df_static$tid)

    # Z display scale (3D only) — applied to display data only; analysis reactives
    # (MSD, JDD, confinement) read from trk() which is unaffected.
    if (use3d) {
      z_stretch <- suppressWarnings(as.numeric(isolate(input$trk_z_stretch %||% 1.0)))
      if (!is.finite(z_stretch) || z_stretch < 0.1) z_stretch <- 1.0
      if (z_stretch != 1.0) {
        scale_z <- function(tb) { if ("z" %in% names(tb)) tb$z <- tb$z * z_stretch; tb }
        sp_full   <- lapply(sp_full,   scale_z)
        sp_tr     <- lapply(sp_tr,     scale_z)
        sp_static <- lapply(sp_static, scale_z)
        if (!is.null(sp_raw)) sp_raw <- lapply(sp_raw, scale_z)
      }
    }

    # Colors: start with random per-TID colors, then overlay confinement or cluster map.
    # Priority: cluster > confinement > base random.
    tids         <- names(sp_full)
    base_col_map <- if (!is.null(attr(df, "col_map"))) attr(df, "col_map") else tid_color_map(tids)
    cluster_cm  <- isolate(cluster_color_map_stored())
    conf_cm     <- isolate(conf_color_map_stored())
    override_cm <- if (!is.null(cluster_cm)) cluster_cm else conf_cm
    col_map <- if (!is.null(override_cm)) {
      cm <- base_col_map
      for (nm in names(override_cm)) if (nm %in% names(cm)) cm[[nm]] <- override_cm[[nm]]
      cm
    } else base_col_map
    
    # Frame grid from bin width OR #frames over [tmin, tmax]
    bin_dt <- suppressWarnings(as.numeric(isolate(input$trk_bin_dt)))
    if (!is.finite(bin_dt) || bin_dt <= 0) bin_dt <- (tmax - tmin) / 120  # fallback
    edges <- seq(tmin, tmax, by = bin_dt)
    if (length(edges) < 2) edges <- c(tmin, tmax)
    if (tail(edges, 1) < (tmax - 1e-12)) edges <- c(edges, tmax)
    
    n_frames <- length(edges) - 1L
    bin_width <- if (length(edges) >= 2) (edges[2] - edges[1]) else NA_real_
    req(n_frames >= 1L)
    
    # Quota params
    target_markers <- suppressWarnings(as.integer(isolate(input$trk_target_markers)))
    if (!is.finite(target_markers) || target_markers < 0) target_markers <- 0L
    per_tid_cap    <- suppressWarnings(as.integer(isolate(input$trk_per_tid_cap)))
    if (!is.finite(per_tid_cap)    || per_tid_cap    < 0) per_tid_cap    <- 0L
    
    # Trail params
    trail_len  <- suppressWarnings(as.numeric(isolate(input$trk_trail_len)));  if (!is.finite(trail_len)  || trail_len  < 0) trail_len  <- 0
    trail_nmax <- suppressWarnings(as.integer(isolate(input$trk_trail_nmax))); if (!is.finite(trail_nmax) || trail_nmax < 1) trail_nmax <- 1L
    
    # ----- Precompute frame membership for markers (O(N_total) once) -----
    # For each TID, split its row-indices by bin index so the per-frame loop
    # does a cheap list-lookup instead of scanning all points every frame.
    # findInterval returns 1..n_frames for points inside [tmin,tmax]; 0 or
    # n_frames+1 for out-of-window points (excluded via the 'valid' mask).
    frame_rows_per_tid <- lapply(tids, function(tid) {
      tb  <- sp_full[[tid]]
      bi  <- findInterval(tb[[tcol_play]], edges)
      valid <- bi >= 1L & bi <= n_frames
      if (!any(valid)) return(list())
      split(which(valid), bi[valid])   # named list: "k" -> integer row-indices
    })

    # ----- Precompute sorted trail-time vectors for O(log n) window lookups -----
    sp_tr_t <- lapply(tids, function(tid) {
      tb <- sp_tr[[tid]]
      if (is.null(tb) || !nrow(tb)) numeric(0) else tb$t
    })

    # Build frames (combined trail + markers)
    frames <- vector("list", n_frames)
    names(frames) <- paste0("f", seq_len(n_frames))

    for (k in seq_len(n_frames)) {
      e0    <- edges[k]
      k_str <- as.character(k)

      # --- Markers (FULL), with quotas per bin ---
      # Fast O(1) lookup from precomputed membership instead of O(n_pts) which()
      per_tid_idx <- lapply(frame_rows_per_tid, function(fr) fr[[k_str]] %||% integer(0))
      bin_counts  <- lengths(per_tid_idx)

      total_in_bin <- sum(bin_counts)
      if (target_markers > 0 && total_in_bin > target_markers) {
        # proportional allocation + largest remainder
        prop   <- bin_counts / total_in_bin
        raw_q  <- prop * target_markers
        q      <- pmin(bin_counts, floor(raw_q))
        if (per_tid_cap > 0) q <- pmin(q, per_tid_cap)
        left   <- max(0L, as.integer(round(target_markers - sum(q))))
        if (left > 0L) {
          gain <- raw_q - q
          add_order <- order(gain, decreasing = TRUE)
          j <- 1L
          while (left > 0L && length(add_order) > 0L) {
            i <- add_order[j]
            if (q[i] < bin_counts[i] && (per_tid_cap == 0L || q[i] < per_tid_cap)) {
              q[i] <- q[i] + 1L; left <- left - 1L
            }
            j <- if (j < length(add_order)) j + 1L else 1L
          }
        }
        quotas <- q
      } else {
        quotas <- bin_counts
        if (per_tid_cap > 0) quotas <- pmin(quotas, per_tid_cap)
      }

      mx <- my <- numeric(0); mz <- numeric(0); mcols <- character(0)

      for (i in seq_along(tids)) {
        w <- per_tid_idx[[i]]
        if (!length(w)) next
        kq <- quotas[i]
        if (kq >= length(w)) {
          keep <- w
        } else if (kq > 0) {
          keep <- w[unique(round(seq(1, length(w), length.out = kq)))]
        } else {
          keep <- integer(0)
        }
        if (length(keep)) {
          tb <- sp_full[[tids[i]]]
          mx    <- c(mx, tb$x[keep])
          my    <- c(my, tb$y[keep])
          if (use3d) mz <- c(mz, tb$z[keep])
          mcols <- c(mcols, rep(unname(as.character(col_map[[tids[i]]])), length(keep)))
        }
      }

      # --- Per-TID trail data: binary search via findInterval (O(log n) per TID) ---
      trails_per_tid <- vector("list", length(tids))
      tx <- ty <- numeric(0)  # combined, used only for dedup check
      if (trail_len > 0) {
        for (i in seq_along(tids)) {
          tt <- sp_tr_t[[i]]
          if (!length(tt)) {
            trails_per_tid[[i]] <- list(x = numeric(0), y = numeric(0))
            next
          }
          # Binary search: points with t in (e0 - trail_len, e0]
          lo_idx <- findInterval(e0 - trail_len, tt) + 1L
          hi_idx <- findInterval(e0, tt, rightmost.closed = TRUE)
          if (hi_idx >= lo_idx) {
            wtr <- seq.int(lo_idx, hi_idx)
            if (length(wtr) > trail_nmax) wtr <- tail(wtr, trail_nmax)
            tb_tr <- sp_tr[[tids[i]]]
            tr_entry <- list(x = tb_tr$x[wtr], y = tb_tr$y[wtr])
            if (use3d && "z" %in% names(tb_tr)) tr_entry$z <- tb_tr$z[wtr]
            trails_per_tid[[i]] <- tr_entry
            tx <- c(tx, tb_tr$x[wtr], NA_real_)
            ty <- c(ty, tb_tr$y[wtr], NA_real_)
          } else {
            trails_per_tid[[i]] <- list(x = numeric(0), y = numeric(0))
          }
        }
      } else {
        for (i in seq_along(tids))
          trails_per_tid[[i]] <- list(x = numeric(0), y = numeric(0))
      }

      frames[[k]] <- list(
        trails_per_tid = trails_per_tid,
        trails_x       = tx,   # combined, used only for dedup
        trails_y       = ty,
        markers_x      = mx,
        markers_y      = my,
        markers_z      = if (use3d) mz else NULL,
        marker_cols    = mcols,
        step_dt        = if (length(edges) >= 2) (edges[2] - edges[1]) else NA_real_,
        bin_width      = bin_width
      )
    }
    
    # --- Collapse identical consecutive frames (same markers + same trails) ---
    same_vec <- function(a, b, tol = 1e-12) {
      length(a) == length(b) && (length(a) == 0 || max(abs(a - b)) <= tol)
    }
    keep <- logical(length(frames)); prev <- NULL
    for (k in seq_along(frames)) {
      fr <- frames[[k]]
      if (is.null(prev)) {
        keep[k] <- TRUE
      } else {
        dup <- same_vec(fr$markers_x, prev$markers_x) &&
          same_vec(fr$markers_y, prev$markers_y) &&
          same_vec(fr$trails_x,  prev$trails_x)  &&
          same_vec(fr$trails_y,  prev$trails_y)
        keep[k] <- !dup
      }
      prev <- fr
    }
    t_grid_full <- edges[-length(edges)]  # start-of-bin time for each frame
    frames  <- frames[keep]
    t_grid_full <- t_grid_full[keep]     # match filtered frames
    
    list(
      static_lines = NULL,
      frames      = frames,
      t_grid      = t_grid_full,
      col_map     = col_map,
      tids        = tids,
      sp          = sp_tr,      # trail animation uses windowed decimated data
      sp_static   = sp_static,  # full-tracks overlay uses unwindowed data
      sp_raw      = sp_raw,     # optional faint background of ALL raw tracks
      trail_len   = trail_len,
      trail_nmax  = trail_nmax,
      bin_width   = bin_width,
      using_bins  = TRUE,
      step_dt     = if (length(edges) >= 2) (edges[2] - edges[1]) else NA_real_,
      has_3d      = use3d      # TRUE when z column is present and 3D mode is active
    )
  }, ignoreNULL = TRUE)
  
  fmt_num <- function(x) ifelse(is.finite(x), format(signif(x, 4), trim = TRUE), "NA")
  hex_to_rgba <- function(hex, alpha = 0.5) {
    rgb_vals <- grDevices::col2rgb(hex)
    sprintf("rgba(%d,%d,%d,%.2f)", rgb_vals[1L], rgb_vals[2L], rgb_vals[3L], alpha)
  }
  # --- show frame dt and time basis ----------------------------------------
  output$trk_frame_meta <- renderUI({
    req(!is.na(plot_trigger()))  # plot_trigger: updates after auto-tune frame rebuilds too
    pre <- frames_store(); req(is.list(pre), length(pre$frames) > 0)
    tcol <- if (isTRUE(is_simul())) "tim0" else "tim"
    dt_show <- pre$step_dt
    if (!is.finite(dt_show) || dt_show <= 0) dt_show <- pre$bin_width
    htmltools::div(
      style = "color:#555;margin-top:6px;margin-bottom:6px;",
      htmltools::tags$b("Time basis: "), tcol, htmltools::HTML("&nbsp;&middot;&nbsp;"),
      htmltools::tags$b("Frames: "), length(pre$frames), htmltools::HTML("&nbsp;&middot;&nbsp;"),
      htmltools::tags$b("Frame dt ~ "), fmt_num(dt_show), " ms"
    )
  })
  
  outputOptions(output, "trk_frame_meta", suspendWhenHidden = FALSE)

  # Wrap plotlyOutput in a renderUI so that every fire_layout() destroys and
  # recreates the DOM element. This forces Plotly to call newPlot() (not react()),
  # which properly initialises the WebGL context on chunk switches and structural
  # changes. Without this, Plotly's in-place react() path blanks scattergl charts.
  output$trk_track2D_wrapper <- renderUI({
    req(!is.na(layout_trigger()))
    plotlyOutput("trk_track2D", height = "1000px")
  })

  # ---- Initial render ----
  output$trk_track2D <- renderPlotly({
    req(!is.na(layout_trigger()))
    # Record the EFFECTIVE view-mode / z-stretch / chunk this render used (with the
    # same defaults the render applies). The view_mode / z_stretch / chunk_sel
    # observers compare against these and skip firing when unchanged — so their
    # INPUT CREATION after an upload (value == the just-rendered default) never
    # triggers a redundant rebuild, without any fragile time window.
    loop$rendered_view  <- isolate(input$trk_view_mode  %||% "2d")
    loop$rendered_zs    <- isolate(input$trk_z_stretch  %||% 1.0)
    loop$rendered_chunk <- isolate(input$trk_chunk_sel  %||% "__all__")
    pre  <- isolate(frames_store()); req(is.list(pre), length(pre$frames) > 0)
    f1    <- pre$frames[[1]]
    tids  <- pre$tids
    cmap  <- pre$col_map
    sp    <- pre$sp
    use3d <- isTRUE(pre$has_3d)

    fig <- plot_ly(source = "track")

    # Faint background of ALL raw tracks (pre$sp_raw, when the Display toggle is on).
    # Drawn first so it sits behind everything; consolidated into one trace.
    draw_raw_faint <- function(fig) {
      spr <- pre$sp_raw
      if (is.null(spr) || !length(spr)) return(fig)
      dk  <- isTRUE(isolate(input$trk_dark_mode))
      col <- if (dk) "rgba(255,255,255,0.10)" else "rgba(120,120,120,0.12)"
      nms <- names(spr)
      rx  <- unlist(lapply(nms, function(nm) c(spr[[nm]]$x, NA_real_)), use.names = FALSE)
      ry  <- unlist(lapply(nms, function(nm) c(spr[[nm]]$y, NA_real_)), use.names = FALSE)
      tr  <- list(x = rx, y = ry, type = if (use3d) "scatter3d" else "scattergl", mode = "lines",
                  line = list(width = 1, color = col), hoverinfo = "skip",
                  showlegend = TRUE, name = "All raw tracks")
      if (use3d) tr$z <- unlist(lapply(nms, function(nm) c(spr[[nm]]$z, NA_real_)), use.names = FALSE)
      do.call(add_trace, c(list(fig), tr))
    }

    # Dwell-time heatmap (2D only) — shared by FULL-tracks and PLAYBACK modes so the
    # toggle works in both. heatmap_store() is computed from trk() on plot_trigger
    # (not playback-specific), so it has data in either mode. Added as a background
    # layer (before the tracks) so paths draw on top.
    hm_data <- if (!use3d && isTRUE(isolate(input$trk_show_heatmap)))
                 tryCatch(isolate(heatmap_store()), error = function(e) NULL) else NULL
    draw_heatmap <- function(fig) {
      if (is.null(hm_data)) return(fig)
      z_proc <- apply_heatmap_contrast(hm_data$z,
                  clip = isolate(input$trk_hm_clip_pct), gamma = isolate(input$trk_hm_gamma))
      fig %>% add_heatmap(x = hm_data$x, y = hm_data$y, z = z_proc, zmin = 0, zmax = 1,
                          showscale = TRUE, opacity = isolate(input$trk_hm_opacity),
                          colorscale = hm_colorscale(isolate(input$trk_hm_colors)),
                          hoverinfo = "skip")
    }

    # Ã¢ÂÂÃ¢ÂÂ Full-tracks static mode Ã¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂ
    if (identical(isolate(input$trk_view_mode_main), "full")) {
      connect  <- isTRUE(isolate(input$trk_full_connect))
      sp_all   <- if (!is.null(pre$sp_static)) pre$sp_static else sp
      all_tids <- names(sp_all)
      dark     <- isTRUE(isolate(input$trk_dark_mode))
      bg_col   <- if (dark) "#000000" else "#FFFFFF"
      fg_col   <- if (dark) "#FFFFFF" else "#63666A"

      fig <- draw_heatmap(fig)     # dwell-time heatmap background (2D only; no-op otherwise)
      fig <- draw_raw_faint(fig)   # all-raw-tracks background (behind the filtered tracks)

      if (connect) {
        # One line+marker trace per TID (allows per-TID color on lines + hover)
        tr_type <- if (use3d) "scatter3d" else "scattergl"
        for (i in seq_along(all_tids)) {
          nm  <- all_tids[i]
          pts <- sp_all[[nm]]
          if (!nrow(pts)) next
          pts <- pts[order(pts$t), ]
          col_i    <- as.character(cmap[[nm]])
          tim_col  <- if ("tim" %in% names(pts)) pts$tim else pts$t
          hover_z  <- if (use3d && "z" %in% names(pts))
                        sprintf("<br>z: %.1f nm", pts$z) else
                        rep("", nrow(pts))
          hover_txt <- sprintf("TID: %s<br>t: %.1f ms<br>x: %.1f nm<br>y: %.1f nm%s",
                               nm, tim_col, pts$x, pts$y, hover_z)
          tr_args <- list(
            x = pts$x, y = pts$y,
            type = tr_type, mode = "lines+markers",
            line   = list(width = 1, color = col_i),
            marker = list(size = isolate(input$trk_pt_size), color = col_i),
            text = hover_txt, hovertemplate = "%{text}<extra></extra>",
            showlegend = (i == 1L), legendgroup = "Localizations",
            name = if (i == 1L) "Localizations" else paste0("tid_", nm)
          )
          if (use3d) tr_args$z <- pts$z
          fig <- do.call(add_trace, c(list(fig), tr_args))
        }
      } else {
        # Single fast scattergl trace with per-point colors + hover text
        tr_type <- if (use3d) "scatter3d" else "scattergl"
        x_v <- c(); y_v <- c(); z_v <- c(); col_v <- c(); txt_v <- c()
        for (nm in all_tids) {
          pts   <- sp_all[[nm]]
          if (!nrow(pts)) next
          col_i   <- as.character(cmap[[nm]])
          tim_col <- if ("tim" %in% names(pts)) pts$tim else pts$t
          hover_z <- if (use3d && "z" %in% names(pts))
                       sprintf("<br>z: %.1f nm", pts$z) else
                       rep("", nrow(pts))
          txt_v <- c(txt_v, sprintf("TID: %s<br>t: %.1f ms<br>x: %.1f nm<br>y: %.1f nm%s",
                                    nm, tim_col, pts$x, pts$y, hover_z))
          x_v   <- c(x_v,   pts$x)
          y_v   <- c(y_v,   pts$y)
          if (use3d) z_v <- c(z_v, pts$z)
          col_v <- c(col_v, rep(col_i, nrow(pts)))
        }
        tr_args <- list(
          x = x_v, y = y_v,
          type = tr_type, mode = "markers",
          marker = list(size = isolate(input$trk_pt_size), color = col_v),
          text = txt_v, hovertemplate = "%{text}<extra></extra>",
          showlegend = TRUE, name = "Localizations", legendgroup = "Localizations"
        )
        if (use3d) tr_args$z <- z_v
        fig <- do.call(add_trace, c(list(fig), tr_args))
      }

      # Axis layout (reuse shared logic below)
      all_x <- unlist(lapply(sp_all, `[[`, "x"), use.names = FALSE)
      all_y <- unlist(lapply(sp_all, `[[`, "y"), use.names = FALSE)
      all_x <- all_x[is.finite(all_x)]; all_y <- all_y[is.finite(all_y)]
      hide_grid <- isTRUE(isolate(input$trk_hide_grid))
      grid_col  <- if (dark) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
      axis_base <- list(showgrid = !hide_grid, zeroline = !hide_grid,
                        showticklabels = !hide_grid,
                        gridcolor = grid_col, color = fg_col)
      ovl <- tryCatch(isolate(trk_image_overlay()), error = function(e) NULL)

      if (use3d) {
        all_z <- unlist(lapply(sp_all, `[[`, "z"), use.names = FALSE)
        lim3  <- get_equal_limits(all_x, all_y, all_z)
        ax3   <- c(list(backgroundcolor = bg_col, autorange = FALSE), axis_base)
        fig   <- fig %>% layout(
          paper_bgcolor = bg_col,
          showlegend = TRUE,                      # always show the "Localizations" toggle
          legend = list(font = list(color = fg_col), groupclick = "togglegroup"),
          scene = list(
            bgcolor = bg_col,
            xaxis = c(list(title = if (hide_grid) "" else "X (nm)", range = lim3$x), ax3),
            yaxis = c(list(title = if (hide_grid) "" else "Y (nm)", range = lim3$y), ax3),
            zaxis = c(list(title = if (hide_grid) "" else "Z (nm)", range = lim3$z), ax3)
          )
        )
      } else {
        xspan <- if (length(all_x) > 1) diff(range(all_x)) else 1000
        yspan <- if (length(all_y) > 1) diff(range(all_y)) else 1000
        half  <- max(xspan, yspan) / 2 * 1.04
        xctr  <- if (length(all_x) > 1) mean(range(all_x)) else 0
        yctr  <- if (length(all_y) > 1) mean(range(all_y)) else 0
        xr    <- c(xctr - half, xctr + half)
        yr    <- c(yctr - half, yctr + half)
        trk_locked_ranges(list(x = xr, y = yr))
        fig <- fig %>% layout(
          paper_bgcolor = bg_col, plot_bgcolor = bg_col,
          xaxis = c(list(title = if (hide_grid) "" else "X (nm)", range = xr), axis_base),
          yaxis = c(list(title = if (hide_grid) "" else "Y (nm)", range = yr,
                         scaleanchor = "x", scaleratio = 1), axis_base),
          images = if (!is.null(ovl)) list(ovl) else list(),
          # showlegend=TRUE forces the legend to appear even with a single track trace
          # (plotly hides it otherwise, and the heatmap adds only a colorbar) so the
          # "Localizations" entry is always clickable to reveal the dwell-time heatmap.
          # groupclick=togglegroup so one click hides ALL track traces (connect mode
          # draws one trace per TID).
          showlegend = TRUE,
          legend = list(font = list(color = fg_col), groupclick = "togglegroup"),
          uirevision = isolate(trk_data_key())
        )
        # Diagnostic markers (same inline logic as the main render path)
        diag <- isolate(rv$trk_diag)
        if (!is.null(diag) && isTRUE(isolate(input$trk_show_diag))) {
          if (!is.null(diag$data_clusters) && nrow(diag$data_clusters) > 0)
            fig <- fig %>% plotly::add_trace(
              x = diag$data_clusters$x, y = diag$data_clusters$y,
              type = "scatter", mode = "markers",
              marker = list(symbol = "cross", size = 16, color = "red",
                            line = list(color = "red", width = 3)),
              name = "MINFLUX clusters", showlegend = TRUE, hoverinfo = "skip")
          if (!is.null(diag$img_blobs) && nrow(diag$img_blobs) > 0)
            fig <- fig %>% plotly::add_trace(
              x = diag$img_blobs$x, y = diag$img_blobs$y,
              type = "scatter", mode = "markers",
              marker = list(symbol = "circle-open", size = 20, color = "#00FF80",
                            line = list(color = "#00FF80", width = 3)),
              name = "Confocal blobs", showlegend = TRUE, hoverinfo = "skip")
        }
        # BAKE the scale bar into the figure (paper-coord shapes), exactly like the
        # SMLM tab — reading inputs ISOLATED so toggling the checkbox does NOT invalidate
        # this render (previously trk_push_scale_bar()'s non-isolated input reads made the
        # plot depend on the checkbox → a rebuild that wiped the proxy-added bar). The
        # checkbox observer still proxy-pushes it live; a genuine rebuild re-bakes it here
        # so it never vanishes. Width uses the rendered x-span; zoom updates via the
        # plotly_relayout observer. trk_axis_x_span is seeded so the live proxy matches.
        trk_axis_x_span(diff(xr))
        if (isTRUE(isolate(input$trk_show_scale_bar))) {
          sb_len   <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_length)))
          sb_thick <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_thickness)))
          sb_font  <- suppressWarnings(as.numeric(isolate(input$trk_scale_font_size)))
          sb <- make_scale_bar(sb_len, diff(xr), sb_thick, dark = dark,
                               scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
          if (!is.null(sb))
            fig <- fig %>% plotly::layout(shapes = sb$shapes, annotations = sb$annotations)
        }
      }
      trk_plot_ready(TRUE)
      return(fig %>% plotly::config(scrollZoom = TRUE, displayModeBar = TRUE) %>%
               plotly::event_register('plotly_relayout'))
    }
    # Ã¢ÂÂÃ¢ÂÂ End full-tracks static mode Ã¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂÃ¢ÂÂ

    # ---- Heatmap (optional, 2D only) — shared helper (see draw_heatmap above) ----
    fig    <- draw_heatmap(fig)
    n_heat <- if (!is.null(hm_data)) 1L else 0L

    # ---- Faint static traces (optional) ----
    # Consolidated into ONE trace with NA separators instead of one per TID.
    # This drops per-TID colouring (opaciOkaty 0.15 grey is visually equivalent),
    # but reduces WebGL trace count from O(N_TIDs) to O(1), which dramatically
    # shortens the initial render time and eliminates the proxy-during-init blank.
    show_faint <- isTRUE(isolate(input$trk_show_full_tracks_faint))
    if (show_faint) {
      # Faint background = the full (post-filtered, post-split) paths of exactly the
      # TIDs being ANIMATED (chunk-scoped), so it never shows tracks that aren't
      # playing. With "All TIDs" selected this is every filtered track.
      #
      # pre$tids is windowed to the playback window [tmin, tmax], but sequential
      # mode plays EVERY chunk TID (seq_data$tids), not just the windowed subset.
      # Use the full played set in sequential mode so every coloured trail that
      # runs has a faint full-track underneath it.
      faint_tids <- pre$tids
      if (isolate(is_seq())) {
        sd_faint <- isolate(seq_data())
        if (is.list(sd_faint) && length(sd_faint$tids)) faint_tids <- sd_faint$tids
      }
      sp_st    <- if (!is.null(pre$sp_static))
                    pre$sp_static[intersect(faint_tids, names(pre$sp_static))] else sp
      all_tids <- names(sp_st)
      static_x <- unlist(lapply(all_tids, function(nm) c(sp_st[[nm]]$x, NA_real_)), use.names = FALSE)
      static_y <- unlist(lapply(all_tids, function(nm) c(sp_st[[nm]]$y, NA_real_)), use.names = FALSE)
      static_z <- if (use3d) unlist(lapply(all_tids, function(nm) c(sp_st[[nm]]$z, NA_real_)), use.names = FALSE) else NULL
      tr_args <- list(
        x = static_x, y = static_y,
        type = if (use3d) "scatter3d" else "scattergl", mode = "lines",
        line = list(width = 1, color = if (isTRUE(isolate(input$trk_dark_mode))) {
          if (use3d) "rgba(255,255,255,0.55)" else "rgba(255,255,255,0.25)"
        } else {
          if (use3d) "rgba(80,80,80,0.55)" else "rgba(80,80,80,0.15)"
        }),
        hoverinfo = "skip", showlegend = TRUE, name = "Full tracks"
      )
      if (use3d) tr_args$z <- static_z
      fig <- do.call(add_trace, c(list(fig), tr_args))
    }
    n_static <- if (show_faint) 1L else 0L

    # ---- Faint ALL-RAW-tracks background (optional; one consolidated trace) ----
    has_raw_bg <- !is.null(pre$sp_raw) && length(pre$sp_raw) > 0
    if (has_raw_bg) fig <- draw_raw_faint(fig)
    n_raw <- if (has_raw_bg) 1L else 0L

    # IMPORTANT: total background traces BEFORE dynamic ones
    n_bg <- n_heat + n_static + n_raw

    # ---- Dynamic setup ----
    show_trails <- isTRUE(isolate(input$trk_show_moving_trails))

    if (isolate(is_seq())) {
      # Sequential playback = the moving trail + a leading marker sweeping through
      # each TID in chronological order. The trail IS the playback here (it's how
      # you watch a TID light up), so it is always present regardless of the
      # simultaneous-mode "Animate moving trails" toggle. The full-track context
      # is drawn faintly underneath (see the faint-tracks block above).
      first_tid <- tids[1]; col1 <- as.character(cmap[[first_tid]])
      tr_type   <- if (use3d) "scatter3d" else "scattergl"
      trail_args  <- list(x = numeric(0), y = numeric(0), type = tr_type, mode = "lines",
                          line = list(width = 3, color = col1), hoverinfo = "skip",
                          showlegend = TRUE, legendgroup = "Localizations", name = "Localizations")
      marker_args <- list(x = numeric(0), y = numeric(0), type = tr_type, mode = "markers",
                          marker = list(size = isolate(input$trk_pt_size), color = col1),
                          hoverinfo = "skip", showlegend = FALSE,
                          legendgroup = "Localizations", name = "marker_seq")
      if (use3d) { trail_args$z <- numeric(0); marker_args$z <- numeric(0) }
      fig <- do.call(add_trace, c(list(fig), trail_args))
      fig <- do.call(add_trace, c(list(fig), marker_args))

      trace_map(list(
        trails_mode = "sequential",
        trail_start = n_bg,
        n_tid       = 0L,
        markers_idx = n_bg + 1L,
        static_idx  = if (show_faint) n_heat else NA_integer_
      ))

    } else {
      if (show_trails) {
        # One line trace per TID, allows per-TID color without per-point tricks.
        # Markers are hidden while trails are active.
        tr_type <- if (use3d) "scatter3d" else "scattergl"
        for (i in seq_along(tids)) {
          col_i   <- as.character(cmap[[tids[i]]])
          tr_args <- list(x = numeric(0), y = numeric(0), type = tr_type, mode = "lines",
                          line = list(width = 3, color = col_i), hoverinfo = "skip",
                          showlegend = (i == 1L), legendgroup = "Localizations",
                          name = if (i == 1L) "Localizations" else paste0("trail_", i))
          if (use3d) tr_args$z <- numeric(0)
          fig <- do.call(add_trace, c(list(fig), tr_args))
        }
        trail_mode  <- "per_tid"
        trail_start <- n_bg
        markers_idx <- n_bg + length(tids)
      } else {
        trail_mode  <- "none"
        trail_start <- NA_integer_
        markers_idx <- n_bg
      }

      # markers ÃÂ¯Â¿Â½ empty when trails are shown (lines replace dots)
      mk_type <- if (use3d) "scatter3d" else "scattergl"
      mk_args <- list(
        x = if (show_trails) numeric(0) else f1$markers_x,
        y = if (show_trails) numeric(0) else f1$markers_y,
        type = mk_type, mode = "markers",
        marker = list(size = isolate(input$trk_pt_size),
                      color = if (show_trails) "rgba(0,0,0,0)" else as.character(f1$marker_cols)),
        hoverinfo = "skip",
        showlegend = !show_trails,
        legendgroup = "Localizations",
        name = "Localizations"
      )
      if (use3d) mk_args$z <- if (show_trails) numeric(0) else (f1$markers_z %||% numeric(0))
      fig <- do.call(add_trace, c(list(fig), mk_args))

      trace_map(list(
        trails_mode = trail_mode,
        trail_start = trail_start,
        n_tid       = if (identical(trail_mode, "per_tid")) length(tids) else 0L,
        markers_idx = markers_idx,
        static_idx  = if (show_faint) n_heat else NA_integer_
      ))
    }
    
    # Bounding box for axis locking.  After Filter & Rechunk, use the stored
    # full spatial extent of all remaining TIDs so the view stays stable while
    # playing through individual TIDs across chunks.  Without a rechunk bbox,
    # fall back to frame 1's markers (or the whole current chunk if frame 1 is empty).
    rechunk_bb <- isolate(trk_rechunk_bbox())
    if (!is.null(rechunk_bb)) {
      bbox_x <- rechunk_bb$x
      bbox_y <- rechunk_bb$y
    } else if (isolate(is_seq())) {
      # Sequential plays one TID at a time, so frame 1 is just the earliest TID's
      # first bin — using it would zoom the view into a tiny region and push the
      # faint full-track context off-screen. Use the full spatial extent of every
      # PLAYED track (seq_data$tids, matching the faint overlay) so every track
      # stays visible while stepping through TIDs.
      seq_tids <- pre$tids
      sd_bb <- isolate(seq_data())
      if (is.list(sd_bb) && length(sd_bb$tids)) seq_tids <- sd_bb$tids
      src <- if (!is.null(pre$sp_static))
               pre$sp_static[intersect(seq_tids, names(pre$sp_static))] else sp
      bbox_x <- unlist(lapply(src, `[[`, "x"), use.names = FALSE)
      bbox_y <- unlist(lapply(src, `[[`, "y"), use.names = FALSE)
      bbox_x <- bbox_x[is.finite(bbox_x)]; bbox_y <- bbox_y[is.finite(bbox_y)]
      if (!length(bbox_x) || !length(bbox_y)) {
        bbox_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
        bbox_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
        bbox_x <- bbox_x[is.finite(bbox_x)]; bbox_y <- bbox_y[is.finite(bbox_y)]
      }
    } else {
      bbox_x <- f1$markers_x[is.finite(f1$markers_x)]
      bbox_y <- f1$markers_y[is.finite(f1$markers_y)]
      if (!length(bbox_x) || !length(bbox_y)) {
        bbox_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
        bbox_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
        bbox_x <- bbox_x[is.finite(bbox_x)]; bbox_y <- bbox_y[is.finite(bbox_y)]
      }
    }
    have_bbox <- length(bbox_x) > 0 && length(bbox_y) > 0
    xr <- yr <- NULL
    if (have_bbox) {
      xr   <- range(bbox_x); yr <- range(bbox_y)
      xpad <- max(1, diff(xr) * 0.04); ypad <- max(1, diff(yr) * 0.04)
      if (!use3d) {
        # 2D only — locked ranges and scale bar depend on these
        trk_axis_x_span(diff(xr) + 2 * xpad)
        trk_locked_ranges(list(x = c(xr[1]-xpad, xr[2]+xpad),
                               y = c(yr[1]-ypad, yr[2]+ypad)))
      }
    }

    ovl       <- tryCatch(isolate(trk_image_overlay()), error = function(e) NULL)
    hide_grid <- isTRUE(isolate(input$trk_hide_grid))
    dark      <- isTRUE(isolate(input$trk_dark_mode))
    bg_col    <- if (dark) "#000000" else "#FFFFFF"
    fg_col    <- if (dark) "#FFFFFF" else "#63666A"
    grid_col  <- if (dark) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
    axis_grid <- list(
      showgrid       = !hide_grid, zeroline   = !hide_grid,
      showticklabels = !hide_grid,
      gridcolor  = grid_col, color = fg_col
    )
    if (use3d) {
      # 3D scene layout — uses scene instead of xaxis/yaxis
      ax3 <- list(color = fg_col, gridcolor = grid_col,
                  showgrid = !hide_grid, zeroline = !hide_grid,
                  showticklabels = !hide_grid, backgroundcolor = bg_col)

      # Equal axis ranges: use the same helper as the SMLM 3D plot.
      # get_equal_limits() finds the largest span across X, Y, Z and applies it
      # to all three axes (centred on each axis's own midpoint) — identical nm/tick.
      all_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
      all_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
      all_z <- unlist(lapply(sp, `[[`, "z"), use.names = FALSE)
      lim3  <- get_equal_limits(all_x, all_y, all_z)

      # autorange = FALSE is required — Plotly 3D ignores range unless autorange is
      # explicitly disabled.
      ax3_x <- c(list(title = if (hide_grid) "" else "X (nm)", range = lim3$x, autorange = FALSE), ax3)
      ax3_y <- c(list(title = if (hide_grid) "" else "Y (nm)", range = lim3$y, autorange = FALSE), ax3)
      ax3_z <- c(list(title = if (hide_grid) "" else "Z (nm)", range = lim3$z, autorange = FALSE), ax3)

      lay <- list(
        title        = list(text = "t = frame 1", font = list(color = fg_col)),
        paper_bgcolor = bg_col,
        font          = list(color = fg_col),
        scene = list(
          bgcolor     = bg_col,
          xaxis       = ax3_x,
          yaxis       = ax3_y,
          zaxis       = ax3_z,
          aspectmode  = "manual",
          aspectratio = list(x = 1, y = 1, z = 1)
        ),
        uirevision  = isolate(trk_data_key()),
        transition  = list(duration = 0)
      )
    } else {
      lay <- list(
        title = list(text = "t = frame 1", font = list(color = fg_col)),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font = list(color = fg_col),
        xaxis = c(
          if (have_bbox) list(title = if (hide_grid) "" else "X (nm)",
                              range = c(xr[1]-xpad, xr[2]+xpad), autorange = FALSE)
          else           list(title = if (hide_grid) "" else "X (nm)"),
          axis_grid
        ),
        yaxis = c(
          if (have_bbox) list(title = if (hide_grid) "" else "Y (nm)",
                              range = c(yr[1]-ypad, yr[2]+ypad), autorange = FALSE,
                              scaleanchor = "x", scaleratio = 1)
          else           list(title = if (hide_grid) "" else "Y (nm)",
                              scaleanchor = "x", scaleratio = 1),
          axis_grid
        ),
        uirevision = isolate(trk_data_key()),
        transition = list(duration = 0)
      )
      if (!is.null(ovl)) lay$images <- list(ovl)
      # ---- Scale bar (2D only) ----
      if (have_bbox && isTRUE(isolate(input$trk_show_scale_bar))) {
        sb_len   <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_length)))
        sb_thick <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_thickness)))
        sb_font  <- suppressWarnings(as.numeric(isolate(input$trk_scale_font_size)))
        sb <- make_scale_bar(sb_len, isolate(trk_axis_x_span()), sb_thick, dark = dark,
                             scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
        if (!is.null(sb)) {
          lay$shapes      <- sb$shapes
          lay$annotations <- sb$annotations
        }
      }
    }
    # ---- Diagnostic: detected cluster centres (red +) and confocal blobs (green ÃÂ¢—Ã¢ÂÂ¹) ----
    # Use type="scatter" (SVG), not scattergl — scattergl lacks circle-open / cross symbols.
    diag <- isolate(rv$trk_diag)
    if (!use3d && isTRUE(isolate(input$trk_show_diag)) && !is.null(diag)) {
      if (!is.null(diag$data_clusters) && nrow(diag$data_clusters) > 0) {
        fig <- fig %>% plotly::add_trace(
          x = diag$data_clusters$x, y = diag$data_clusters$y,
          type = "scatter", mode = "markers",
          marker = list(symbol = "cross", size = 16, color = "red",
                        line = list(color = "red", width = 3)),
          name = "MINFLUX clusters", showlegend = TRUE, hoverinfo = "skip"
        )
      }
      if (!is.null(diag$img_blobs) && nrow(diag$img_blobs) > 0) {
        fig <- fig %>% plotly::add_trace(
          x = diag$img_blobs$x, y = diag$img_blobs$y,
          type = "scatter", mode = "markers",
          marker = list(symbol = "circle-open", size = 20, color = "#00FF80",
                        line = list(color = "#00FF80", width = 3)),
          name = "Confocal blobs", showlegend = TRUE, hoverinfo = "skip"
        )
      }
    }

    trk_plot_ready(TRUE)
    do.call(plotly::layout, c(list(fig), lay)) %>% plotly::event_register('plotly_relayout')
  })

  # After every layout rebuild, restore the locked axis ranges via proxy (2D only).
  # lr is read INSIDE the delay so it captures the value written by renderPlotly
  # (which runs in the same flush cycle as layout_trigger but finishes before the
  # 120ms callback fires).  Reading it outside would capture the stale pre-render value.
  observeEvent(layout_trigger(), {
    if (isTRUE(isolate(use_3d()))) return()   # 3D camera manages its own view
    shinyjs::delay(120, {
      lr <- isolate(trk_locked_ranges())
      if (!is.null(lr)) {
        plotlyProxyInvoke(
          plotlyProxy("trk_track2D", session), "relayout",
          list(`xaxis.range[0]` = lr$x[1], `xaxis.range[1]` = lr$x[2],
               `yaxis.range[0]` = lr$y[1], `yaxis.range[1]` = lr$y[2],
               `xaxis.autorange` = FALSE,   `yaxis.autorange` = FALSE)
        )
      }
    })
  }, ignoreInit = TRUE)

  cap_for <- function(n) {
    if (n <= 200) 0L else ceiling(sqrt(n) * 6)
  }
  
  
  # ---- Draw helpers ----
  draw_frame <- function(k, force = FALSE) {
    if (identical(isolate(input$trk_view_mode_main), "full")) {
      return(invisible(NULL))
    }
    # Pull current prebuilt frames and mapping safely
    # Use if/return guards (not req())  -  req() throws inside later() callbacks
    pre <- isolate(frames_store())
    if (!is.list(pre) || !length(pre$frames)) { return(invisible(NULL)) }
    tm  <- isolate(trace_map())
    if (is.null(tm)) { return(invisible(NULL)) }
    pt_size <- suppressWarnings(isolate(input$trk_pt_size))

    # Respect play/pause unless explicitly forced
    if (!force && !(isTRUE(loop$gplaying) || isTRUE(loop$playing))) { return(invisible(NULL)) }
    
    # Clamp k and extract frame
    k  <- max(1L, min(as.integer(k), length(pre$frames)))
    fr <- pre$frames[[k]]
    
    proxy <- plotlyProxy("trk_track2D", session)

    is3d <- isTRUE(pre$has_3d)

    # ----- Trails (per-TID line traces, batched into ONE proxy call) -----
    if (identical(tm$trails_mode, "per_tid")) {
      tpt <- if (!is.null(fr$trails_per_tid)) fr$trails_per_tid else list()
      n   <- tm$n_tid
      xs  <- vector("list", n); ys <- vector("list", n); zs <- vector("list", n)
      for (i in seq_len(n)) {
        td      <- if (i <= length(tpt)) tpt[[i]] else NULL
        xs[[i]] <- if (!is.null(td) && length(td$x)) td$x else NA_real_
        ys[[i]] <- if (!is.null(td) && length(td$y)) td$y else NA_real_
        if (is3d) zs[[i]] <- if (!is.null(td) && length(td$z)) td$z else NA_real_
      }
      restyle_data <- list(x = xs, y = ys)
      if (is3d) restyle_data$z <- zs
      plotlyProxyInvoke(proxy, "restyle", restyle_data,
        seq(tm$trail_start, tm$trail_start + n - 1L))
      # Markers stay empty while trails are shown
      empty_mk <- list(x = list(NA_real_), y = list(NA_real_))
      if (is3d) empty_mk$z <- list(NA_real_)
      plotlyProxyInvoke(proxy, "restyle", empty_mk, tm$markers_idx)
    } else {
      # ----- Markers (no trails active) -----
      if (length(fr$markers_x)) {
        mk_data <- list(
          x = list(fr$markers_x),
          y = list(fr$markers_y),
          `marker.size`  = if (is.finite(pt_size)) pt_size else 6,
          `marker.color` = list(as.character(fr$marker_cols))
        )
        if (is3d) mk_data$z <- list(fr$markers_z %||% numeric(0))
        plotlyProxyInvoke(proxy, "restyle", mk_data, tm$markers_idx)
      } else {
        empty_data <- list(x = list(NA_real_), y = list(NA_real_),
                           `marker.size` = if (is.finite(pt_size)) pt_size else 6)
        if (is3d) empty_data$z <- list(NA_real_)
        plotlyProxyInvoke(proxy, "restyle", empty_data, tm$markers_idx)
      }
    }

    # ----- Title (+ range re-pin while trails animate) -----
    # A title relayout re-runs plotly's layout solver, which \u2014 with scaleanchor on \u2014
    # re-fits the axes and jitters the FOV. When trails are animating we therefore fold
    # the current view's ranges into the SAME (atomic) relayout, every frame, so the view
    # can't drift. We pin to the USER'S CURRENT VIEW (trk_zoom_ranges, updated live by the
    # plotly_relayout observer when they zoom/pan) and fall back to the full extent only
    # when they haven't zoomed \u2014 so zooming in DURING playback sticks instead of snapping
    # back. Without trails the markers are one stable-extent trace, so the plain title
    # relayout (throttled) is safe. (3D manages its own camera via aspectmode.)
    playing_now <- isTRUE(loop$gplaying) || isTRUE(loop$playing)
    trails_on   <- identical(tm$trails_mode, "per_tid")
    if (trails_on || !playing_now || (k %% 3L) == 1L) {
      n_pts <- length(fr$markers_x)
      rl <- list(title = list(text = sprintf("frame %d \u00B7 markers: %d", k, n_pts)))
      if (trails_on && !is3d) {
        lr <- isolate(trk_zoom_ranges())              # live user view (zoom/pan)
        if (is.null(lr)) lr <- isolate(trk_locked_ranges())   # else full extent
        if (!is.null(lr)) rl <- c(rl, list(
          `xaxis.range[0]` = lr$x[1], `xaxis.range[1]` = lr$x[2],
          `yaxis.range[0]` = lr$y[1], `yaxis.range[1]` = lr$y[2],
          `xaxis.autorange` = FALSE,  `yaxis.autorange` = FALSE))
      }
      plotlyProxyInvoke(proxy, "relayout", rl)
    }

    # ----- Slider sync (throttled to every 5 frames during playback) -----
    if (!is.null(isolate(input$trk_frame_idx))) {
      if (!isTRUE(loop$gplaying) && !isTRUE(loop$playing) || (k %% 5L) == 1L) {
        slider_guard(isolate(slider_guard()) + 1L)
        updateSliderInput(session, "trk_frame_idx", value = k, max = length(pre$frames))
      }
    }
    
    invisible(k)
  }
  
  smart_defaults <- function(df, mode = c("simul","seq"), target_frames = 200, target_fps = 15) {
    mode <- match.arg(mode)
    tcol <- if (mode == "simul") "tim0" else "tim"
    
    # robust time span
    rng <- range(df[[tcol]], na.rm = TRUE)
    total_units <- if (mode == "simul") max(df$tim0, na.rm = TRUE) else diff(rng)
    if (!is.finite(total_units) || total_units <= 0) total_units <- 1
    
    # round up to a nice number for the playback window
    nice_total <- max(pretty(c(0, total_units), n = 5))
    
    # frame count / width
    n_frames <- target_frames
    bin_dt   <- nice_total / n_frames
    
    # trail ~5% of window, at least one frame wide
    trail_len <- max(nice_total * 0.05, bin_dt)
    
    # marker quota ~ avg points per frame (*1.5 safety), clipped
    n_pts <- nrow(df)
    target_markers <- max(200, min(8000, round((n_pts / n_frames) * 1.5)))
    
    # sequential timing: choose ms_per_unit to get ~target_fps
    # use median per-TID duration for stability
    durs <- df |>
      dplyr::group_by(tid) |>
      dplyr::summarise(dur = max(.data[[tcol]]) - min(.data[[tcol]]), .groups = "drop") |>
      dplyr::pull(dur)
    med_dur <- suppressWarnings(stats::median(durs[is.finite(durs) & durs > 0]))
    if (!is.finite(med_dur) || med_dur <= 0) med_dur <- total_units
    step_units   <- med_dur / n_frames           # units per frame (seq)
    frame_ms     <- 1000 / target_fps            # ~15 fps
    ms_per_unit  <- max(1, round(frame_ms / max(step_units, .Machine$double.eps)))
    
    list(
      nice_total      = nice_total,
      n_frames        = n_frames,
      bin_dt          = bin_dt,
      trail_len       = trail_len,
      target_markers  = target_markers,
      ms_per_unit     = ms_per_unit
    )
  }  
  
  # Draw one step of the CURRENT TID by point index (1..nrow(tb))
  draw_tid_point <- function(tid, i) {
    # Read reactives safely (no dependencies)
    pre <- isolate(frames_store()); if (is.null(pre)) return(invisible(NULL))
    tm  <- isolate(trace_map());    if (is.null(tm))  return(invisible(NULL))
    sd  <- isolate(seq_data());     if (is.null(sd) || !(tid %in% sd$tids)) return(invisible(NULL))
    
    tb <- sd$sp[[tid]]; if (is.null(tb) || !nrow(tb)) return(invisible(NULL))
    n  <- nrow(tb)
    i  <- max(1L, min(as.integer(i), n))
    
    tnow  <- tb$t[i]
    is3d  <- isTRUE(pre$has_3d)
    has_z_col <- is3d && "z" %in% names(tb)
    w <- which(tb$t >= (tnow - pre$trail_len) & tb$t <= tnow)
    if (length(w) > pre$trail_nmax) w <- tail(w, pre$trail_nmax)
    xx <- if (length(w)) c(tb$x[w], NA_real_) else numeric(0)
    yy <- if (length(w)) c(tb$y[w], NA_real_) else numeric(0)
    zz <- if (has_z_col && length(w)) c(tb$z[w], NA_real_) else numeric(0)

    col <- as.character(pre$col_map[[tid]])
    proxy <- plotlyProxy("trk_track2D", session)

    # Trail (single trace in sequential mode) — only when the trail trace exists
    # (i.e. "Animate moving trails" is on). Otherwise it's a plain moving point.
    if (!is.na(tm$trail_start)) {
      trail_data <- list(x = list(xx), y = list(yy), `line.color` = list(col))
      if (is3d) trail_data$z <- list(zz)
      plotlyProxyInvoke(proxy, "restyle", trail_data, tm$trail_start)
    }

    # Marker
    mk_data <- list(x = list(tb$x[i]), y = list(tb$y[i]),
                    `marker.size`  = isolate(input$trk_pt_size),
                    `marker.color` = list(col))
    if (has_z_col) mk_data$z <- list(tb$z[i])
    plotlyProxyInvoke(proxy, "restyle", mk_data, tm$markers_idx)
    
    plotlyProxyInvoke(proxy, "relayout",
                      list(title = paste0("tid = ", tid,
                                          " (", which(pre$tids == tid), "/", length(pre$tids), ")",
                                          ", idx = ", i,
                                          ", t = ", signif(tnow, 6))))
    
    # Keep the slider in sync (ok to use isolate here)
    if (!is.null(isolate(input$trk_frame_idx))) {
      slider_guard(isolate(slider_guard()) + 1L)
      updateSliderInput(session, "trk_frame_idx", value = i, max = n)
    }
  }
  
  
  #--------HEATMAP-------
  hm_colorscale <- function(name = "fire") {
    switch(name,
      magenta = list(
        c(0.00, "#000000"), c(0.40, "#3b0f3f"),
        c(0.70, "#a100a1"), c(1.00, "#ffe6ff")
      ),
      fire = list(
        c(0.00, "#000000"), c(0.25, "#4d0000"),
        c(0.50, "#b30000"), c(0.75, "#ff8c00"),
        c(1.00, "#ffff00")
      ),
      # all NUMERIC_PALETTES entries fall through to plotly_colorscale
      plotly_colorscale(name)
    )
  }
  apply_heatmap_contrast <- function(zmat, clip = c(1, 99), gamma = 1) {
    zvec <- as.numeric(zmat)
    zvec <- zvec[is.finite(zvec)]
    if (!length(zvec)) return(zmat)
    
    lo <- stats::quantile(zvec, clip[1] / 100, names = FALSE, type = 7)
    hi <- stats::quantile(zvec, clip[2] / 100, names = FALSE, type = 7)
    if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
      lo <- min(zvec, na.rm = TRUE)
      hi <- max(zvec, na.rm = TRUE) + 1e-9
    }
    
    zclip <- pmin(pmax(zmat, lo), hi)
    z01   <- (zclip - lo) / (hi - lo)
    
    if (is.finite(gamma) && gamma > 0 && abs(gamma - 1) > 1e-9) {
      # gamma < 1 brightens faint values; gamma > 1 compresses lows
      z01 <- z01 ^ gamma
    }
    z01
  }
  
  build_dwell_heatmap <- function(df, nbins = 80) {
    # df needs columns x, y, tid, tim
    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
    if (!nrow(df)) return(NULL)
    
    df <- df[order(df$tid, df$tim), ]
    # dwell weight = time until next sample within same tid
    df <- df %>%
      group_by(tid) %>%
      mutate(dt = lead(tim) - tim) %>%
      ungroup()
    df$dt[!is.finite(df$dt) | df$dt < 0] <- 0
    
    xr <- range(df$x, na.rm = TRUE); if (xr[1] == xr[2]) xr <- xr + c(-0.5, 0.5)
    yr <- range(df$y, na.rm = TRUE); if (yr[1] == yr[2]) yr <- yr + c(-0.5, 0.5)
    
    xb <- seq(xr[1], xr[2], length.out = nbins + 1)
    yb <- seq(yr[1], yr[2], length.out = nbins + 1)
    ix <- findInterval(df$x, xb, all.inside = TRUE)
    iy <- findInterval(df$y, yb, all.inside = TRUE)
    
    # full grid with zeros
    df$ix <- factor(pmax(1, pmin(nbins, ix)), levels = 1:nbins)
    df$iy <- factor(pmax(1, pmin(nbins, iy)), levels = 1:nbins)
    
    z <- xtabs(df$dt ~ df$iy + df$ix)     # rows=y, cols=x
    z <- as.matrix(z)
    
    xcent <- 0.5 * (xb[-1] + xb[-length(xb)])
    ycent <- 0.5 * (yb[-1] + yb[-length(yb)])
    list(x = xcent, y = ycent, z = z)
  }
  auto_hm_bins <- function(npts) {
    if (npts < 2e5) 300 else if (npts < 1e6) 200 else 120
  }
  
  heatmap_store <- eventReactive(plot_trigger(), {
    df <- isolate(trk()); req(nrow(df) > 0)
    bins <- if (isTruthy(isolate(input$trk_hm_bins))) isolate(input$trk_hm_bins) else auto_hm_bins(nrow(df))
    build_dwell_heatmap(df[, c("x","y","tid","tim")], nbins = bins)
  }, ignoreInit = TRUE)
  
  # --- dt extraction per TID (based on raw 'tim') --------------------------
  dt_vectors <- reactive({
    df <- trk(); req(nrow(df) > 0)
    df <- df[is.finite(df$tim), ]
    df <- df[order(df$tid, df$tim), ]
    sp  <- split(df$tim, df$tid)
    lapply(sp, function(tt) {
      if (length(tt) < 2) return(numeric(0))
      dt <- diff(tt)
      dt[is.finite(dt) & dt >= 0]
    })
  })
  # --- dt histogram ---------------------------------------------------------

  # --- ÃÂt histogram in analysis tab (richer version with stats) ---------------
  output$trk_dt_hist_jdd <- renderPlotly({
    req(isTRUE(input$trk_show_tsa), tsa_ran())
    dts  <- dt_vectors()
    vec  <- unlist(dts, use.names = FALSE)
    vec  <- vec[is.finite(vec) & vec > 0]
    shiny::validate(shiny::need(length(vec) > 0,
      "No Δt values — need at least 2 time points per TID."))
    n_bins  <- max(10L, as.integer(input$trk_dt_hist_bins %||% 80L))
    dark    <- isTRUE(isolate(input$trk_dark_mode))
    bg_col  <- if (dark) "#000000" else "#FFFFFF"
    fg_col  <- if (dark) "#FFFFFF" else "#333333"

    dt_med  <- median(vec, na.rm = TRUE)
    # Bin-FREE typical Δt (half-sample mode) — independent of the bin slider.
    dt_hsm  <- half_sample_mode(vec)
    h_tmp   <- hist(vec, breaks = n_bins, plot = FALSE)   # for display + y-scaling only

    fig <- plot_ly(x = vec, type = "histogram", nbinsx = n_bins,
                   marker = list(color = "rgba(100,150,220,0.7)",
                                 line  = list(color = "rgba(70,110,180,0.9)",
                                              width = 0.5))) %>%
      add_trace(x = c(dt_med, dt_med), y = c(0, max(h_tmp$counts) * 1.1),
                type = "scatter", mode = "lines",
                line   = list(color = "#E03020", width = 2, dash = "dash"),
                marker = list(size = 0),
                name = sprintf("Median %.3g ms", dt_med)) %>%
      add_trace(x = c(dt_hsm, dt_hsm), y = c(0, max(h_tmp$counts) * 1.1),
                type = "scatter", mode = "lines",
                line   = list(color = "#20A040", width = 2, dash = "dashdot"),
                marker = list(size = 0),
                name = sprintf("Typical Δt (bin-free) %.3g ms", dt_hsm)) %>%
      layout(
        xaxis  = list(title = "Δt between consecutive localizations (ms)",
                      color = fg_col),
        yaxis  = list(title = "Count", color = fg_col),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font   = list(color = fg_col),
        legend = list(orientation = "h", x = 0, y = 1.12,
                      font = list(color = fg_col)),
        bargap = 0.05,
        uirevision = "dt_jdd"
      )
    fig
  })

  output$trk_dt_stats_ui <- renderUI({
    req(isTRUE(input$trk_show_tsa), tsa_ran())
    dts <- dt_vectors()
    vec <- unlist(dts, use.names = FALSE)
    vec <- vec[is.finite(vec) & vec > 0]
    req(length(vec) > 0)
    dt_hsm  <- half_sample_mode(vec)      # bin-free typical Δt
    dt_med  <- median(vec)
    dt_p10  <- quantile(vec, 0.10)
    dt_p25  <- quantile(vec, 0.25)
    dt_p75  <- quantile(vec, 0.75)
    hw      <- trk_hw_cycle_time()
    tags$div(style = "font-size:0.84em; color:#555; margin:4px 0;",
      tags$div(
        tags$b("Typical Δt (bin-free, half-sample mode): "),
        tags$span(style = "color:#20803a;", sprintf("%.4g ms", dt_hsm)),
        sprintf("   ·   median %.4g   ·   p10 %.4g   ·   p25 %.4g   ·   p75 %.4g   ·   n = %s",
                dt_med, dt_p10, dt_p25, dt_p75, format(length(vec), big.mark = ","))
      ),
      if (is.finite(hw) && hw > 0)
        tags$div(style = "color:#888; margin-top:2px;",
          sprintf("Hardware atomic cycle = %.4g ms. If the typical Δt is ≈2× this, dark cycles dominate (the molecule is often briefly unseen between localizations).", hw)),
      tags$div(style = "color:#888; margin-top:2px; font-size:0.95em;",
        "This estimate is independent of the histogram bin slider (which is display-only).")
    )
  })

  #-----MSD time binned
  # Returns per-bin mean square displacement, pair count, AND the sum of the
  # actual Δt of all pairs in each bin.  The Δt sum lets the caller use the true
  # mean lag time per bin (Σdt / n_pairs) as the τ value, instead of the geometric
  # bin centre — important for MINFLUX where pairs within a wide bin are NOT
  # centred on the midpoint (far more short-lag pairs than long-lag pairs).
  msd_time_binned_one <- function(x, y, t, tau_edges, z = NULL) {
    n     <- length(t); nb <- length(tau_edges) - 1L
    use_z <- !is.null(z) && length(z) == n
    if (n < 2 || nb < 1)
      return(list(msd = rep(NA_real_, nb), npairs = integer(nb), dt_sum = double(nb)))
    msd_acc <- double(nb); np_acc <- integer(nb); dt_acc <- double(nb)
    tau_lo <- tau_edges[1L]; tau_hi <- tau_edges[nb + 1L]

    for (i in seq_len(n - 1L)) {
      # t is sorted; binary-search for last j where dt <= tau_hi
      dt_all <- t[seq.int(i + 1L, n)] - t[i]
      j_hi   <- findInterval(tau_hi, dt_all, rightmost.closed = TRUE)
      if (j_hi < 1L) next
      dt <- dt_all[seq_len(j_hi)]
      ok <- dt >= tau_lo          # dt <= tau_hi already ensured by j_hi cutoff
      if (!any(ok)) next
      dt_ok <- dt[ok]
      b     <- findInterval(dt_ok, tau_edges, rightmost.closed = TRUE)
      in_b  <- b >= 1L & b <= nb
      if (!any(in_b)) next
      js2    <- seq.int(i + 1L, i + j_hi)[ok][in_b]
      b2     <- b[in_b]
      dt_in  <- dt_ok[in_b]
      dx  <- x[js2] - x[i]; dy <- y[js2] - y[i]
      dr2 <- dx * dx + dy * dy
      if (use_z) { dz <- z[js2] - z[i]; dr2 <- dr2 + dz * dz }
      # rowsum() is a fast C-level grouped sum
      rs  <- rowsum(dr2, b2, reorder = FALSE)
      idx <- as.integer(rownames(rs))
      msd_acc[idx] <- msd_acc[idx] + rs[, 1L]
      rt  <- rowsum(dt_in, b2, reorder = FALSE)
      idt <- as.integer(rownames(rt))
      dt_acc[idt] <- dt_acc[idt] + rt[, 1L]
      np_acc <- np_acc + tabulate(b2, nbins = nb)
    }
    list(msd = ifelse(np_acc > 0L, msd_acc / np_acc, NA_real_),
         npairs = np_acc, dt_sum = dt_acc)
  }

  # Per-track diffusion-sign filter (Bela's TAMSD D>0 rule). Given a per-track MSD
  # matrix (rows = lags/bins, cols = tracks) and the lag-time vector, fit each
  # track's MSD ~ a + b·tau by OLS over its finite points and return a logical
  # mask of columns whose slope b > 0 (positive diffusion). Tracks with < 2 finite
  # points are kept (cannot judge). Used to drop non-diffusing / artefact tracks
  # before the ensemble in ALL three MSD pipelines.
  .track_positive_D <- function(MSD_mat, tau) {
    if (is.null(dim(MSD_mat))) MSD_mat <- matrix(MSD_mat, ncol = 1L)
    apply(MSD_mat, 2L, function(m) {
      ok <- is.finite(m) & is.finite(tau)
      if (sum(ok) < 2L) return(TRUE)
      x <- tau[ok]; y <- m[ok]
      vx <- sum((x - mean(x))^2)
      if (vx <= 0) return(TRUE)
      b <- sum((x - mean(x)) * (y - mean(y))) / vx
      is.finite(b) && b > 0
    })
  }

  # ---  build per-TID MSD and aggregate -------------------------------------
  msd_bundle <- reactive({
    # Only recompute when the "Plot MSD" button is clicked.
    # All settings are read inside isolate() so changing them alone does not
    # trigger a recompute — the user must click the button again.
    req(isTRUE(input$trk_show_msd))
    req(input$trk_msd_compute > 0)

    df <- isolate(trk()); req(nrow(df) > 0)
    
    # raw chronological order for diffusion
    # Detect genuine 3D data: z column must exist with meaningful variation (SD > 1 nm).
    # A z column of zeros or NAs (2D experiment) must NOT trigger 3D MSD.
    has_real_z <- "z" %in% names(df) && {
      zv <- df$z[is.finite(df$z)]
      length(zv) > 10 && sd(zv) > 1
    }
    # For 3D, also require finite z so NA z never poisons pair accumulation.
    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim) &
               if (has_real_z) is.finite(df$z) else TRUE, ]
    df <- df[order(df$tid, df$tim), ]
    sp  <- split(df, df$tid)

    # Gap splitting is now handled in the main filter pipeline (trk_gap_split,
    # Stage 1.5).  trk() already contains split segments when the filter is on,
    # so no additional splitting is needed here.
    min_pts <- isolate(max(5L, as.integer(input$trk_msd_min_pts)))

    # Build per-track MSDs (time-binned only ÃÂ¯Â¿Â½ correct for irregular MINFLUX dt)
    per <- vector("list", length(sp))
    names(per) <- names(sp)

    tau_max <- isolate(suppressWarnings(as.numeric(input$trk_msd_tau_max)))
    if (!is.finite(tau_max) || tau_max <= 0) tau_max <- 2000

    # τ-bin width: auto-detect from modal Δt (cycle-aware) when enabled, so the
    # default on MINFLUX data is cycle resolution rather than a coarse 10 ms bin
    # that pools many lags into the first point and biases the intercept.
    auto_bin <- isolate(isTRUE(input$trk_msd_tau_bin_auto))
    if (auto_bin) {
      vmod <- unlist(isolate(dt_vectors()), use.names = FALSE)
      vmod <- vmod[is.finite(vmod) & vmod > 0]
      d_tau <- if (length(vmod)) half_sample_mode(vmod) else NA_real_
      if (!is.finite(d_tau) || d_tau <= 0) d_tau <- 10
    } else {
      d_tau <- isolate(suppressWarnings(as.numeric(input$trk_msd_tau_bin)))
      if (!is.finite(d_tau) || d_tau <= 0) d_tau <- 10
    }
    tau_edges   <- seq(0, tau_max, by = d_tau)
    if (length(tau_edges) < 2) tau_edges <- c(0, tau_max)
    tau_centers <- 0.5 * (tau_edges[-1] + tau_edges[-length(tau_edges)])

    withProgress(message = "Computing MSD...", value = 0, {
      for (k in seq_along(sp)) {
        incProgress(1/length(sp), detail = names(sp)[k])
        tb <- sp[[k]]
        if (nrow(tb) < min_pts) next
        out <- msd_time_binned_one(
          tb$x, tb$y, tb$tim, tau_edges,
          z = if (has_real_z) tb$z else NULL
        )
        per[[k]] <- list(msd = out$msd, npairs = out$npairs, dt_sum = out$dt_sum)
      }
    })

    keep <- !vapply(per, is.null, logical(1))
    if (!any(keep)) return(list(tau = numeric(0), msd = numeric(0), dt0 = NA_real_))
    per <- per[keep]

    nb  <- length(tau_centers)
    MSD <- do.call(cbind, lapply(per, function(z) {
      v <- z$msd; if (length(v) < nb) c(v, rep(NA_real_, nb - length(v))) else v
    }))
    NP  <- do.call(cbind, lapply(per, function(z) {
      v <- z$npairs; if (length(v) < nb) c(v, rep(0L, nb - length(v))) else v
    }))
    DTS <- do.call(cbind, lapply(per, function(z) {
      v <- z$dt_sum; if (length(v) < nb) c(v, rep(0, nb - length(v))) else v
    }))

    # True mean Δt per bin (pair-weighted across all tracks). Use this as the τ
    # value instead of the geometric bin centre — removes the systematic τ bias
    # that occurs when a bin is wide relative to the cycle time. Empty bins fall
    # back to the geometric centre.
    dt_sum_bin <- rowSums(DTS, na.rm = TRUE)
    np_sum_bin <- rowSums(NP,  na.rm = TRUE)
    tau_actual <- ifelse(np_sum_bin > 0, dt_sum_bin / np_sum_bin, tau_centers)

    # Drop per-track TAMSDs with non-positive diffusion (negative slope) — Bela's
    # D>0 filter. Removes non-diffusing / artefact tracks before the ensemble.
    # Guarded so it never blanks the plot (only fires when some tracks are +D and
    # some are −D).
    posD <- .track_positive_D(MSD, tau_actual)
    n_dropped_negD <- sum(!posD)
    if (any(posD) && !all(posD)) {
      MSD <- MSD[, posD, drop = FALSE]; NP <- NP[, posD, drop = FALSE]; DTS <- DTS[, posD, drop = FALSE]
      dt_sum_bin <- rowSums(DTS, na.rm = TRUE)
      np_sum_bin <- rowSums(NP,  na.rm = TRUE)
      tau_actual <- ifelse(np_sum_bin > 0, dt_sum_bin / np_sum_bin, tau_centers)
    } else {
      n_dropped_negD <- 0L
    }

    # Always pair-weighted (smoother, robust): long tracks contribute more at
    # longer lags.  Track-weighted is computed only where needed for ergodicity.
    {
      W <- NP; W[W == 0] <- NA
      w   <- rowSums(W, na.rm = TRUE)
      num <- rowSums(MSD * W, na.rm = TRUE)
      agg <- ifelse(w > 0, num / w, NA_real_)
    }

    pair_sum         <- rowSums(NP, na.rm = TRUE)
    n_tracks_per_bin <- rowSums(!is.na(MSD))

    # Mask ensemble mean where too few tracks contribute.
    # Threshold: user-set fraction of the peak track count, with a hard floor of 2.
    min_trk_pct  <- isolate(suppressWarnings(as.numeric(input$trk_msd_min_trk_pct)))
    if (!is.finite(min_trk_pct) || min_trk_pct <= 0) min_trk_pct <- 10
    min_trk_abs  <- max(2L, floor((min_trk_pct / 100) * max(n_tracks_per_bin, na.rm = TRUE)))
    agg[n_tracks_per_bin < min_trk_abs] <- NA_real_

    dts <- isolate(dt_vectors())
    dt0 <- median(unlist(dts, use.names = FALSE), na.rm = TRUE)
    if (!is.finite(dt0) || dt0 <= 0) dt0 <- 1

    # Per-track IQR band — same threshold as the ensemble mean for consistency
    msd_q25 <- apply(MSD, 1L, quantile, probs = 0.25, na.rm = TRUE)
    msd_q75 <- apply(MSD, 1L, quantile, probs = 0.75, na.rm = TRUE)
    msd_q25[n_tracks_per_bin < min_trk_abs] <- NA_real_
    msd_q75[n_tracks_per_bin < min_trk_abs] <- NA_real_

    ndim_sel <- isolate(input$trk_msd_ndim %||% "auto")
    ndim_msd <- if (identical(ndim_sel, "3")) 3L else if (identical(ndim_sel, "2")) 2L else if (has_real_z) 3L else 2L

    list(tau = tau_actual, msd = agg, dt0 = dt0, npairs = pair_sum,
         tau_centers = tau_centers,     # geometric bin centres (kept for reference)
         d_tau = d_tau, tau_bin_auto = auto_bin,
         msd_q25 = msd_q25, msd_q75 = msd_q75,
         n_tracks = ncol(MSD),          # tracks kept AFTER the D>0 filter
         n_dropped_negD = n_dropped_negD,
         MSD_matrix = MSD,              # nrow=nb, ncol=n_valid_tracks (for spaghetti + histograms)
         n_tracks_per_bin = n_tracks_per_bin,
         ndim = ndim_msd)               # 2 or 3 — controls fit coefficient (4D vs 6D)
  })
  #---- MSD export ----
  msd_export <- reactive({
    ub <- msd_unified_bundle(); req(ub)
    linear_tau <- isTRUE(input$trk_uni_linear_tau)
    tau_all <- if (linear_tau) ub$tau_lin else ub$tau
    ok  <- is.finite(tau_all) & is.finite(ub$msd) & ub$keep_lag
    fr  <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
    td  <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)

    data <- data.frame(
      cycle_lag = ub$lag[ok],
      tau_ms    = tau_all[ok],
      MSD_um2   = ub$msd[ok] / 1e6,
      n_pairs   = if (!is.null(ub$npairs)) ub$npairs[ok] else NA_integer_,
      n_tracks  = ub$n_trk_lag[ok],
      stringsAsFactors = FALSE
    )

    params <- list(
      mode               = "unified MSD (cumulative cycle-index binning)",
      time_unit          = "ms",
      space_unit         = "µm",
      tau_axis           = if (linear_tau) "linearized lag × cycle" else "real mean Δt per bin",
      ensemble_weighting = if (isTRUE(ub$track_weighted)) "track-weighted (eMSD)" else "pair-weighted",
      fit_method         = if (isTRUE(input$trk_uni_iterative)) "iterative OLSF" else "anchored OLSF + Michalet",
      fit_window_cycles  = sprintf("%s-%s", input$trk_uni_min_lag %||% 5, input$trk_uni_max_lag %||% 50),
      cycle_time_ms      = ub$t_cycle,
      msd_ndim           = ub$ndim,
      eDIFF_um2_s        = if (!is.null(fr)) fr$D / 1e3 else NA_real_,
      tDIFF_um2_s        = if (!is.null(td)) td$tDIFF else NA_real_,
      n_tracks           = ub$n_tracks,
      msd_formula        = sprintf("MSD = %dD*tau + %d*sigma^2", 2L * ub$ndim, 2L * ub$ndim)
    )

    list(data = data, params = params)
  })
  msd_single_for_tb <- function(tb, tau_edges, min_pts = 5L) {
    if (nrow(tb) < min_pts) return(NULL)
    out <- msd_time_binned_one(
      tb$x, tb$y, tb$tim,
      tau_edges = tau_edges,
      z = if ("z" %in% names(tb)) tb$z else NULL
    )
    tau <- 0.5 * (tau_edges[-1] + tau_edges[-length(tau_edges)])
    data.frame(tau = tau, msd = out$msd, npairs = out$npairs)
  }
  msd_single <- reactive({
    req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_show_single), nzchar(input$trk_msd_tid))
    df <- trk(); req(nrow(df) > 0)
    has_z_s <- "z" %in% names(df) && {
      zv <- df$z[is.finite(df$z)]; length(zv) > 10 && sd(zv) > 1
    }
    cols_s <- if (has_z_s) c("x","y","z","tim") else c("x","y","tim")
    tb <- df[df$tid == input$trk_msd_tid, cols_s]
    if (has_z_s) tb <- tb[is.finite(tb$z), ]
    tau_max <- suppressWarnings(as.numeric(input$trk_msd_tau_max))
    if (!is.finite(tau_max) || tau_max <= 0) tau_max <- 2000
    # Match the ensemble bin: auto = modal Δt, else the manual value
    if (isTRUE(input$trk_msd_tau_bin_auto)) {
      vmod <- unlist(dt_vectors(), use.names = FALSE)
      vmod <- vmod[is.finite(vmod) & vmod > 0]
      d_tau <- if (length(vmod)) half_sample_mode(vmod) else NA_real_
      if (!is.finite(d_tau) || d_tau <= 0) d_tau <- 10
    } else {
      d_tau <- suppressWarnings(as.numeric(input$trk_msd_tau_bin))
      if (!is.finite(d_tau) || d_tau <= 0) d_tau <- 10
    }
    tau_edges <- seq(0, tau_max, by = d_tau)
    if (length(tau_edges) < 2) tau_edges <- c(0, tau_max)
    msd_single_for_tb(tb, tau_edges, min_pts = max(5L, input$trk_msd_min_pts))
  })
  # Per-track D and alpha fits — computed when histogram OR clustering is open
  # Per-track fits feeding the D histogram, the D/CI clustering and the per-track
  # CSV export \u2014 all derived from the SAME Unified MSD calculation (per-track OLSF
  # D + confinement index), via the cached trk_unified_tdiff reactive.
  msd_pertk_fits <- reactive({
    req(isTRUE(input$trk_show_msd))
    req(isTRUE(input$trk_msd_show_hist) || isTRUE(input$trk_msd_cluster))
    td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(td) && !is.null(td$pertrack) && nrow(td$pertrack) > 0,
      "Enable the Unified MSD and click 'Plot MSD' to compute per-track diffusion."))
    td$pertrack
  })

  # --- MSD cluster analysis ---------------------------------------------------
  # k-means on (log10 D, log10 CI); clusters sorted ascending by mean log10(D) so that
  # cluster 1 = most confined / slowest, cluster k = most mobile / fastest.
  MSD_CLUSTER_PAL <- c("#3060C0","#C03030","#2CA02C","#FF7F0E","#9467BD","#8C564B")

  msd_cluster_rv <- reactive({
    req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
    fits <- msd_pertk_fits()
    req(!is.null(fits) && nrow(fits) >= 3L)

    k <- max(2L, min(6L, as.integer(input$trk_msd_n_clusters %||% 3L)))
    req(nrow(fits) >= k)

    # Feature matrix: alpha (primary) + log10(CI) (secondary).
    # After z-scoring, alpha is doubled so it is the dominant axis.
    # With k=2 the split is primarily along alpha (confined vs free/directed);
    # CI adds refinement for k >= 3.
    req(!is.null(fits$conf_index))
    ci_safe   <- pmax(fits$conf_index, 1e-4)
    logD_safe <- log10(pmax(fits$D_um2_s, 1e-9))
    feat <- scale(cbind(logD_safe, log10(ci_safe)))
    feat[!is.finite(feat)] <- 0
    feat[, 1] <- feat[, 1] * 2  # log10(D) is the primary separator

    set.seed(42L)
    km <- tryCatch(
      kmeans(feat, centers = k, nstart = 25L, iter.max = 200L),
      error = function(e) NULL
    )
    req(!is.null(km))

    # Sort clusters by mean log10(D) ascending: cluster 1 = most confined, cluster k = most mobile
    mean_logD <- tapply(logD_safe, km$cluster, mean, na.rm = TRUE)
    rank_ord  <- order(mean_logD)
    new_lbl  <- integer(nrow(fits))
    for (r in seq_len(k)) new_lbl[km$cluster == rank_ord[r]] <- r

    pal     <- MSD_CLUSTER_PAL[seq_len(k)]
    col_map <- setNames(pal[new_lbl], as.character(fits$tid))

    fits$cluster    <- new_lbl
    lbl_names <- c("confined","intermediate","mobile","fast","faster","fastest")
    fits$cluster_lbl <- paste0("Cluster ", new_lbl, " (", lbl_names[new_lbl], ")")

    list(col_map = col_map, palette = pal, n_clusters = k,
         fits = fits, cluster_ids = setNames(new_lbl, as.character(fits$tid)))
  })

  # Observer: keep cluster_color_map_stored in sync
  observe({
    if (!isTRUE(input$trk_msd_cluster_color_tracks)) {
      cluster_color_map_stored(NULL)
      return()
    }
    cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
    cluster_color_map_stored(if (!is.null(cl)) cl$col_map else NULL)
  })

  # Rebuild frames whenever stored cluster map changes
  observeEvent(cluster_color_map_stored(), {
    fire_frames()
  }, ignoreInit = TRUE, ignoreNULL = FALSE)

  # Also rebuild when checkbox is toggled or k changes
  observeEvent(input$trk_msd_cluster_color_tracks, {
    fire_frames()
  }, ignoreInit = TRUE)
  observeEvent(input$trk_msd_n_clusters, {
    if (isTRUE(input$trk_msd_cluster_color_tracks)) fire_frames()
  }, ignoreInit = TRUE)

  # ---------------------------------------------------------------------------
  # ---- Change-point segmentation + GMM clustering ---------------------------
  # Splits each track into diffusion-state segments via PELT on log(r²/ÃÂt),
  # then fits a GMM (mclust, BIC-selected G) to the per-segment D distribution.
  # ---------------------------------------------------------------------------
  CPD_PAL <- c("#3060C0","#C03030","#2CA02C","#FF7F0E","#9467BD","#8C564B")

  cpd_segments_rv <- eventReactive(input$trk_cpd_compute, {
    req(isTRUE(input$trk_show_cpd))
    df <- tryCatch(trk_spatial_filtered(), error = function(e) NULL)
    req(!is.null(df), nrow(df) > 0,
        all(c("x","y","tid","tim") %in% names(df)))

    min_seg    <- max(3L, as.integer(input$trk_cpd_min_seg %||% 8L))
    use_z_cpd  <- isTRUE(has_z())
    ndim       <- if (use_z_cpd) 3L else 2L
    coef_D     <- 2L * ndim   # 4 (2D) or 6 (3D): MSD = coef_D * D * dt

    tids <- split(df, df$tid)

    seg_rows <- withProgress(message = "Segmenting tracks…", value = 0, {
      lapply(seq_along(tids), function(i) {
        incProgress(1 / length(tids))
        td <- tids[[i]][order(tids[[i]]$tim), ]
        n  <- nrow(td)
        if (n < 2L * min_seg) return(NULL)

        dt <- diff(td$tim)                       # ms (already converted)
        if (any(!is.finite(dt) | dt <= 0)) return(NULL)
        r2 <- diff(td$x)^2 + diff(td$y)^2       # nm²
        if (use_z_cpd && "z" %in% names(td) && any(is.finite(td$z)))
          r2 <- r2 + diff(td$z)^2               # include z for 3D data

        # Normalise step size by time lag Ã¢ÂÂ value Ã¢ÂÂ local D
        # Log-transform Ã¢ÂÂ approximately normal; mean shifts between D-states
        r2n <- pmax(r2 / dt, .Machine$double.eps)
        log_r2n <- log(r2n)

        # PELT change-point on the mean of log(r²/ÃÂt)
        cpt_obj <- tryCatch(
          changepoint::cpt.mean(log_r2n, method = "PELT", penalty = "BIC",
                                minseglen = min_seg),
          error = function(e) NULL
        )
        bkpts <- if (!is.null(cpt_obj)) changepoint::cpts(cpt_obj) else integer(0)
        seg_s <- c(1L, if (length(bkpts)) bkpts + 1L else integer(0))
        seg_e <- c(bkpts, length(log_r2n))

        tid_val <- as.character(td$tid[1])
        do.call(rbind, Filter(Negate(is.null), lapply(seq_along(seg_s), function(s) {
          idx <- seg_s[s]:seg_e[s]
          if (length(idx) < min_seg) return(NULL)
          # Use mean(r²/ÃÂt) / coef_D, not mean(r²) / (coef_D * mean(ÃÂt)):
          # the latter is biased when r² and ÃÂt covary (Jensen's inequality).
          D_nm2_ms <- mean(r2[idx] / dt[idx]) / coef_D
          data.frame(
            tid      = tid_val,
            seg_idx  = s,
            n_states = length(seg_s),     # total segments in this track
            n_steps  = length(idx),
            D_nm2_ms = D_nm2_ms,
            D_um2_s  = D_nm2_ms / 1e3,
            stringsAsFactors = FALSE
          )
        })))
      })
    })

    result <- do.call(rbind, Filter(Negate(is.null), seg_rows))
    shiny::validate(shiny::need(
      !is.null(result) && nrow(result) >= 3L,
      "No segments found — try lowering Min segment length."
    ))
    result
  }, ignoreNULL = TRUE)

  cpd_gmm_rv <- reactive({
    req(isTRUE(input$trk_show_cpd))
    segs <- tryCatch(cpd_segments_rv(), error = function(e) NULL)
    req(!is.null(segs), nrow(segs) >= 3L)
    log10_D <- log10(segs$D_um2_s)
    log10_D <- log10_D[is.finite(log10_D)]
    req(length(log10_D) >= 3L)
    max_G <- max(1L, min(6L, as.integer(input$trk_cpd_max_g %||% 4L)))
    tryCatch(
      mclust::Mclust(log10_D, G = seq_len(max_G), modelNames = "V", verbose = FALSE),
      error = function(e) NULL
    )
  })

  output$trk_cpd_hist <- renderPlotly({
    req(isTRUE(input$trk_show_cpd))
    segs <- tryCatch(cpd_segments_rv(), error = function(e) {
      shiny::validate(shiny::need(FALSE, conditionMessage(e))); NULL
    })
    req(!is.null(segs))
    gmm <- tryCatch(cpd_gmm_rv(), error = function(e) NULL)

    dark   <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"

    log10_D_fin <- log10(segs$D_um2_s)
    log10_D_fin <- log10_D_fin[is.finite(log10_D_fin)]
    n_bins <- 40L

    fig <- plot_ly(
      x    = log10_D_fin,
      type = "histogram",
      nbinsx = n_bins,
      marker = list(color = "rgba(100,120,210,0.55)",
                    line  = list(color = "rgba(80,100,190,0.8)", width = 0.5)),
      name = "Segments"
    )

    if (!is.null(gmm) && gmm$G >= 1L) {
      x_seq  <- seq(min(log10_D_fin) - 0.5, max(log10_D_fin) + 0.5, length.out = 400)
      props  <- gmm$parameters$pro
      means  <- gmm$parameters$mean
      sds    <- sqrt(gmm$parameters$variance$sigmasq)
      bw     <- diff(range(log10_D_fin)) / n_bins
      scale  <- length(log10_D_fin) * bw     # histogram Ã¢ÂÂ density scale factor

      total_y <- numeric(length(x_seq))
      for (g in seq_len(gmm$G)) {
        y_g      <- props[g] * dnorm(x_seq, means[g], sds[g]) * scale
        total_y  <- total_y + y_g
        lbl      <- sprintf("Pop. %d  D =%.3g µm²/s  (%.1f%%)",
                            g, 10^means[g], props[g] * 100)
        fig <- fig %>% add_trace(
          x = x_seq, y = y_g, type = "scatter", mode = "lines",
          line = list(color = CPD_PAL[g], width = 2.5),
          name = lbl
        )
      }
      fig <- fig %>% add_trace(
        x = x_seq, y = total_y, type = "scatter", mode = "lines",
        line = list(color = "rgba(200,30,30,0.75)", width = 1.5, dash = "dash"),
        name = "GMM total"
      )
    }

    fig %>% layout(
      barmode       = "overlay",
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font          = list(color = fg_col),
      xaxis = list(title = "log₁₀(D)  (µm²/s)", color = fg_col,
                   gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
      yaxis = list(title = "Segments", color = fg_col,
                   gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
      legend = list(orientation = "h", x = 0, y = 1.12, font = list(size = 11))
    )
  })

  output$trk_cpd_summary <- renderUI({
    req(isTRUE(input$trk_show_cpd))
    segs <- tryCatch(cpd_segments_rv(), error = function(e) NULL)
    req(!is.null(segs))
    gmm  <- tryCatch(cpd_gmm_rv(),      error = function(e) NULL)

    n_tracks   <- length(unique(segs$tid))
    n_segs     <- nrow(segs)
    n_multi    <- sum(tapply(segs$n_states, segs$tid, max) > 1L)

    gmm_line <- if (!is.null(gmm)) {
      G <- gmm$G
      paste0("GMM selected ", G, " component",
             if (G > 1L) "s" else "", " by BIC:  ",
             paste0(sapply(seq_len(G), function(g)
               sprintf("Pop.%d D =%.3g µm²/s (%.1f%%)",
                       g, 10^gmm$parameters$mean[g],
                       gmm$parameters$pro[g] * 100)),
               collapse = "  |  "))
    } else "GMM not yet computed."

    tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
      tags$b("Segmentation:"), " ",
      n_tracks, " tracks  •  ",
      n_segs, " segments  •  ",
      n_multi, " tracks with ≥2 states detected", tags$br(),
      tags$b("GMM: "), gmm_line
    )
  })

  output$trk_dl_cpd_csv <- downloadHandler(
    filename = function() paste0("cpd_segments_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      segs <- tryCatch(cpd_segments_rv(), error = function(e) NULL)
      gmm  <- tryCatch(cpd_gmm_rv(),      error = function(e) NULL)
      if (!is.null(segs) && !is.null(gmm)) {
        log10_D   <- log10(segs$D_um2_s)
        segs$gmm_component <- ifelse(
          is.finite(log10_D),
          mclust::predict.Mclust(gmm, log10_D)$classification,
          NA_integer_
        )
      }
      write.csv(segs, file, row.names = FALSE)
    }
  )

  # --- Shared MSD fit helpers --------------------------------------------------
  # ndim = 2 (2D): MSD = 4*D*tau^alpha + 4*sigma^2
  # ndim = 3 (3D): MSD = 6*D*tau^alpha + 6*sigma^2
  # Units throughout: tau in ms, MSD in nm², D in nm²/ms.
  # Display conversion: D [µm²/s] = D [nm²/ms] / 1e3   (1 nm²/ms = 0.001 µm²/s)
  #                     MSD [ÃÂµm²] = MSD [nm²]   / 1e6

  # OLSF — faithful port of Vogler's OptimizedLeastSquaresFitting (Michalet &
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

  # auto_pmin = FALSE: fit the first k points (caller-chosen window).
  # auto_pmin = TRUE : two-stage iterative Michalet optimal point selection.
  # `dt_cycle` is the cycle time (for x); `Npts` the trajectory length for L(N).
  msd_fit_olsf <- function(tau, msd, npairs, k, ndim = 2L,
                           auto_pmin = FALSE, R = 0, Npts = NA_real_,
                           dt_cycle = NA_real_) {
    coef_val <- 2L * ndim  # 4 (2D) or 6 (3D)
    ok_all <- is.finite(tau) & is.finite(msd) & tau > 0
    tau <- tau[ok_all]; msd <- msd[ok_all]
    n_avail <- length(tau)
    if (n_avail < 3L) return(NULL)

    # Unweighted linear fit on the first p points → D, σ² (his curve_fit).
    fit_p <- function(p) {
      p  <- max(2L, min(as.integer(p), n_avail))
      s <- tau[seq_len(p)]; y <- msd[seq_len(p)]
      ok <- is.finite(s) & is.finite(y)
      if (sum(ok) < 2L) return(NULL)
      cf <- tryCatch(coef(lm(y[ok] ~ s[ok])), error = function(e) NULL)
      if (is.null(cf) || any(!is.finite(cf))) return(NULL)
      list(D = cf[2] / coef_val, sig2 = cf[1] / coef_val, p = p)
    }

    if (!auto_pmin) {
      r <- fit_p(max(3L, min(as.integer(k), n_avail)))
      if (is.null(r)) return(NULL)
      D_est <- max(0, r$D); sig2_est <- max(0, r$sig2); p_slope <- r$p; p_int <- r$p
    } else {
      dt <- if (is.finite(dt_cycle) && dt_cycle > 0) dt_cycle
            else if (n_avail >= 2) stats::median(diff(tau)) else tau[1]
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
    }
    if (!is.finite(D_est)) return(NULL)
    if (!is.finite(sig2_est)) sig2_est <- 0
    sig_est <- sqrt(sig2_est)
    tau_seq <- seq(min(tau), max(tau), length.out = 200)
    yfit    <- pmax(0, coef_val * D_est * tau_seq + coef_val * sig_est^2)

    list(D = D_est, alpha = 1.0, sigma_loc = sig_est,
         tau_fit = tau_seq, msd_fit = yfit, n_points = p_slope,
         n_points_intercept = p_int,
         forced_linear = TRUE, ndim = ndim, method = "OLSF")
  }

  # Bela's EXACT iterative OLSF (optimized_least_squares_fitting.py). Unlike the
  # anchored variant above, it does NOT pre-fit the whole window — it starts at a
  # few points and iterates the optimal point count until it repeats, separately
  # for slope (→ D) and intercept (→ σ²). This deliberately converges on the
  # LOWEST lags (least error-prone), reproducing his numbers. `tau` is the lag
  # axis (use the linearized lag×cycle axis for a faithful match); D = slope/(2d).
  msd_fit_olsf_iterative <- function(tau, msd, ndim = 2L, R = 0,
                                     dt_cycle = NA_real_, N = NA_real_,
                                     max_attempts = 1000L) {
    coef_val <- 2L * ndim
    ok <- is.finite(tau) & is.finite(msd) & tau > 0
    tau <- tau[ok]; msd <- msd[ok]
    n_avail <- length(tau)
    if (n_avail < 3L) return(NULL)
    d  <- ndim
    if (!is.finite(N) || N <= 0) N <- n_avail
    if (!is.finite(dt_cycle) || dt_cycle <= 0)
      dt_cycle <- if (n_avail >= 2) stats::median(diff(tau)) else tau[1]

    fit_p <- function(p) {
      p <- max(2L, min(as.integer(p), n_avail))
      s <- tau[seq_len(p)]; y <- msd[seq_len(p)]
      cf <- tryCatch(stats::coef(stats::lm(y ~ s)), error = function(e) NULL)
      if (is.null(cf) || any(!is.finite(cf))) return(NULL)
      list(D = cf[2] / coef_val, sig2 = cf[1] / coef_val, p = p)
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
    tau_seq <- seq(min(tau), max(tau), length.out = 200)
    yfit <- pmax(0, coef_val * D_est * tau_seq + coef_val * sig_est^2)
    list(D = D_est, alpha = 1.0, sigma_loc = sig_est,
         tau_fit = tau_seq, msd_fit = yfit, n_points = p_s_fin,
         n_points_intercept = p_i_fin,
         forced_linear = TRUE, ndim = ndim, method = "OLSF-iterative")
  }

  msd_fit_result <- function(tau, msd, k, anomalous = TRUE, ndim = 2L) {
    df <- data.frame(tau = tau[seq_len(k)], msd = msd[seq_len(k)])
    df <- df[is.finite(df$tau) & is.finite(df$msd) & df$tau > 0, ]
    if (nrow(df) < 3L) return(NULL)
    coef_val <- 2L * ndim   # 4 for 2D, 6 for 3D

    D0     <- max(1e-6, (max(df$msd) - min(df$msd)) / (coef_val * max(df$tau)))
    sig0   <- max(0, sqrt(abs(df$msd[1]) / coef_val))
    alpha0 <- 1.0
    alpha_se <- NA_real_   # parametric SE on alpha from the ensemble fit (set below)

    fm_anom <- as.formula(sprintf("msd ~ %d*D*tau^alpha + %d*sig2", coef_val, coef_val))
    fm_lin  <- as.formula(sprintf("msd ~ %d*D*tau + %d*sig2",       coef_val, coef_val))

    if (anomalous) {
      fit <- tryCatch(
        nls(fm_anom, data = df,
            start   = list(D = D0, alpha = alpha0, sig2 = sig0^2),
            lower   = list(D = 0,  alpha = 0.05,   sig2 = 0),
            upper   = list(D = Inf, alpha = 2.5,   sig2 = Inf),
            algorithm = "port",
            control = nls.control(maxiter = 200, warnOnly = TRUE)),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        fit <- tryCatch(
          nls(fm_lin, data = df,
              start = list(D = D0, sig2 = sig0^2),
              lower = list(D = 0, sig2 = 0),
              algorithm = "port",
              control = nls.control(maxiter = 200, warnOnly = TRUE)),
          error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        cf <- coef(fit)
        D_est <- unname(cf["D"]); alpha_est <- 1.0
        sig_est <- sqrt(max(0, unname(cf["sig2"]))); forced_linear <- TRUE
      } else {
        cf <- coef(fit)
        D_est <- unname(cf["D"]); alpha_est <- unname(cf["alpha"])
        sig_est <- sqrt(max(0, unname(cf["sig2"]))); forced_linear <- FALSE
        # Parametric SE on alpha from the fit covariance. NB this treats the averaged
        # MSD points as independent, so it is a LOWER bound on the true uncertainty
        # (MSD lags are correlated) — surfaced in the UI with that caveat.
        alpha_se <- tryCatch(unname(sqrt(diag(vcov(fit)))["alpha"]), error = function(e) NA_real_)
      }
    } else {
      fit <- tryCatch(
        nls(fm_lin, data = df,
            start = list(D = D0, sig2 = sig0^2),
            lower = list(D = 0, sig2 = 0),
            algorithm = "port",
            control = nls.control(maxiter = 200, warnOnly = TRUE)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      cf <- coef(fit)
      D_est <- unname(cf["D"]); alpha_est <- 1.0
      sig_est <- sqrt(max(0, unname(cf["sig2"]))); forced_linear <- FALSE
    }

    tau_seq <- seq(min(df$tau), max(df$tau), length.out = 200)
    yfit    <- pmax(0, coef_val * D_est * tau_seq^alpha_est + coef_val * sig_est^2)

    list(D = D_est, alpha = alpha_est, alpha_se = alpha_se, sigma_loc = sig_est,
         tau_fit = tau_seq, msd_fit = yfit,
         forced_linear = if (anomalous) forced_linear else FALSE,
         ndim = ndim)
  }

  # Shared fit-window selector — applied consistently by the ensemble estimates,
  # the plotted fit line, and the per-track fits so they always agree.
  # Lower bound: min_lag_ms.  Upper bound: fit_max_ms (absolute) when > 0, else
  # fit_pct % of the remaining points.  Returns a logical mask over `tau`.
  msd_fit_mask <- function(tau, min_lag_ms = 0, fit_max_ms = 0, fit_pct = 30) {
    keep <- is.finite(tau)
    if (min_lag_ms > 0) keep <- keep & tau >= min_lag_ms
    idx  <- which(keep)
    if (!length(idx)) return(keep)
    if (fit_max_ms > 0) {
      keep <- keep & tau <= fit_max_ms
    } else {
      k <- max(3L, floor((fit_pct / 100) * length(idx)))
      cutoff <- sort(tau[idx])[min(k, length(idx))]
      keep <- keep & tau <= cutoff
    }
    keep
  }

  # --- MSD plot ---------------------------------------------------------------
  output$trk_msd_plot <- renderPlotly({
    mb <- msd_bundle()
    shiny::validate(shiny::need(length(mb$tau) > 0 && length(mb$msd) > 0,
      "No MSD to display yet \u2014 try lowering 'Min points per TID' or widening tau bins."))

    ok  <- is.finite(mb$tau) & is.finite(mb$msd)
    tau <- mb$tau[ok]; msd <- mb$msd[ok]
    shiny::validate(shiny::need(length(tau) >= 2, "Not enough finite MSD points to plot."))
    tau_max_shown <- max(tau)   # spaghetti clipped to this — beyond it the mean is masked

    fit_pct   <- 30
    min_lag_ms_f <- msd_min_lag_ms()
    fit_max_ms_f <- msd_fit_max_ms()
    t_cyc_f      <- trk_cycle_time()
    anomalous <- isTRUE(input$trk_msd_anomalous)
    use_olsf   <- !anomalous
    auto_pmin_f <- use_olsf
    R_blur_f    <- 0
    # Fit window [min_lag, fit_max] (Vogler fit_limits); OLSF picks points within.
    fmask_p <- rep(TRUE, length(tau))
    if (min_lag_ms_f > 0) fmask_p <- fmask_p & tau >= min_lag_ms_f
    if (fit_max_ms_f > 0) fmask_p <- fmask_p & tau <= fit_max_ms_f
    if (!use_olsf) fmask_p <- fmask_p & msd_fit_mask(tau, min_lag_ms_f, fit_max_ms_f, fit_pct)
    tau_fit_v <- tau[fmask_p]; msd_fit_v <- msd[fmask_p]
    npr_fit_v <- (mb$npairs %||% rep(1L, length(mb$msd)))[ok][fmask_p]
    k         <- max(3L, length(tau_fit_v))
    show_fit  <- isTRUE(input$trk_msd_show_fit)
    show_band <- isTRUE(input$trk_msd_show_band)
    show_spag  <- isTRUE(input$trk_msd_show_spaghetti)
    show_ntau  <- isTRUE(input$trk_msd_show_n_tau)
    norm_tau   <- FALSE

    # Scale MSD from nm\u00B2 to \u00B5m\u00B2 for display (1 nm\u00B2 = 1e-6 \u00B5m\u00B2)
    msd_scale <- 1e-6

    # Helper: divide y by x, element-wise, skipping zeros/NAs
    norm_y <- function(y_vec, x_vec) ifelse(is.finite(x_vec) & x_vec > 0, y_vec / x_vec, NA_real_)

    y_label <- if (norm_tau) "MSD(\u03C4)/\u03C4  [\u00B5m\u00B2/ms]" else "MSD(\u03C4) [\u00B5m\u00B2]"

    plt <- plot_ly()

    # 1. Spaghetti — individual per-track TAMSDs
    #    When MSD cluster coloring is active: one colored trace per cluster.
    #    Otherwise: single grey concatenated trace (NA breaks between tracks).
    if (show_spag && !is.null(mb$MSD_matrix) && ncol(mb$MSD_matrix) > 0) {
      MSD_mat   <- mb$MSD_matrix
      max_lines <- 200L
      ncols     <- ncol(MSD_mat)

      # Try to get cluster assignments (tid -> cluster label)
      cl_rv <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
      use_cluster_spag <- isTRUE(input$trk_msd_cluster) && !is.null(cl_rv)

      set.seed(42L)
      idx <- if (ncols > max_lines) sort(sample(ncols, max_lines)) else seq_len(ncols)

      if (use_cluster_spag) {
        # One trace per cluster, colored; tracks not in fits shown as grey
        k_cls  <- cl_rv$n_clusters
        pal    <- cl_rv$palette
        tid_cl <- cl_rv$cluster_ids  # named integer: tid -> cluster label (1..k)

        # Get column TIDs from MSD_matrix column names (set during msd_bundle)
        col_tids <- colnames(MSD_mat)

        for (cl in seq_len(k_cls)) {
          spag_x <- numeric(0); spag_y <- numeric(0)
          n_cl   <- 0L
          for (i in idx) {
            tid_i <- if (!is.null(col_tids)) col_tids[i] else NA_character_
            if (is.na(tid_i) || !(tid_i %in% names(tid_cl))) next
            if (tid_cl[tid_i] != cl) next
            msd_i <- MSD_mat[, i]
            ok_i  <- is.finite(msd_i) & mb$tau <= tau_max_shown
            if (sum(ok_i) < 2L) next
            tau_i <- mb$tau[ok_i]
            y_i   <- msd_scale * (if (norm_tau) norm_y(msd_i[ok_i], tau_i) else msd_i[ok_i])
            spag_x <- c(spag_x, tau_i,  NA_real_)
            spag_y <- c(spag_y, y_i,    NA_real_)
            n_cl   <- n_cl + 1L
          }
          if (length(spag_x) > 0L && n_cl > 0L) {
            cl_hex  <- pal[cl]
            cl_rgba <- hex_to_rgba(cl_hex, alpha = 0.30)
            plt <- plt %>% add_trace(
              x = spag_x, y = spag_y,
              type = "scattergl", mode = "lines",
              line = list(color = cl_rgba, width = 1),
              name = sprintf("Cluster %d (n=%d)", cl, n_cl),
              hoverinfo = "skip"
            )
          }
        }

        # Grey trace for tracks with no cluster assignment in sampled set
        spag_x_oth <- numeric(0); spag_y_oth <- numeric(0)
        for (i in idx) {
          tid_i <- if (!is.null(col_tids)) col_tids[i] else NA_character_
          if (!is.na(tid_i) && tid_i %in% names(tid_cl)) next  # has cluster, already drawn
          msd_i <- MSD_mat[, i]
          ok_i  <- is.finite(msd_i) & mb$tau <= tau_max_shown
          if (sum(ok_i) < 2L) next
          tau_i <- mb$tau[ok_i]
          y_i   <- msd_scale * (if (norm_tau) norm_y(msd_i[ok_i], tau_i) else msd_i[ok_i])
          spag_x_oth <- c(spag_x_oth, tau_i,  NA_real_)
          spag_y_oth <- c(spag_y_oth, y_i,    NA_real_)
        }
        if (length(spag_x_oth) > 0L) {
          plt <- plt %>% add_trace(
            x = spag_x_oth, y = spag_y_oth,
            type = "scattergl", mode = "lines",
            line = list(color = "rgba(140,140,140,0.18)", width = 1),
            name = "Unassigned TAMSDs", hoverinfo = "skip"
          )
        }

      } else {
        # Plain grey spaghetti — single concatenated trace
        spag_x <- numeric(0); spag_y <- numeric(0)
        for (i in idx) {
          msd_i <- MSD_mat[, i]
          ok_i  <- is.finite(msd_i) & mb$tau <= tau_max_shown
          if (sum(ok_i) < 2L) next
          tau_i <- mb$tau[ok_i]
          y_i   <- msd_scale * (if (norm_tau) norm_y(msd_i[ok_i], tau_i) else msd_i[ok_i])
          spag_x <- c(spag_x, tau_i,  NA_real_)
          spag_y <- c(spag_y, y_i,    NA_real_)
        }
        if (length(spag_x) > 0L) {
          plt <- plt %>% add_trace(
            x = spag_x, y = spag_y,
            type = "scattergl", mode = "lines",
            line = list(color = "rgba(120,120,120,0.18)", width = 1),
            name = sprintf("Individual TAMSDs (n=%d)", length(idx)),
            hoverinfo = "skip"
          )
        }
      }
    }

    # 2. IQR shaded band
    if (show_band && !is.null(mb$msd_q25)) {
      ok_b <- ok & is.finite(mb$msd_q25) & is.finite(mb$msd_q75)
      if (any(ok_b)) {
        tau_b <- mb$tau[ok_b]
        q25_b <- msd_scale * (if (norm_tau) norm_y(mb$msd_q25[ok_b], tau_b) else mb$msd_q25[ok_b])
        q75_b <- msd_scale * (if (norm_tau) norm_y(mb$msd_q75[ok_b], tau_b) else mb$msd_q75[ok_b])
        plt <- plt %>%
          add_trace(
            x = c(tau_b, rev(tau_b)),
            y = c(q75_b, rev(q25_b)),
            type = "scatter", mode = "none", fill = "toself",
            fillcolor = "rgba(100,160,220,0.20)",
            line = list(width = 0),
            name = "IQR (25\u201375%)", showlegend = TRUE, hoverinfo = "skip"
          )
      }
    }

    # 3. Ensemble mean MSD
    msd_plot <- msd_scale * (if (norm_tau) norm_y(msd, tau) else msd)
    plt <- plt %>%
      add_trace(x = tau, y = msd_plot,
                type = "scattergl", mode = "markers",
                marker = list(size = 5, color = "#3060C0"),
                name = if (norm_tau) "Mean MSD(\u03C4)/\u03C4" else "Mean MSD(\u03C4)")

    # 4. Fit curve — uses the same fit window as the estimates panel
    ndim <- mb$ndim %||% 2L
    if (show_fit && length(tau_fit_v) >= 3L) {
      fr <- if (use_olsf)
        msd_fit_olsf(tau_fit_v, msd_fit_v, npr_fit_v, k, ndim = ndim,
                     auto_pmin = auto_pmin_f, R = R_blur_f, Npts = 1e4, dt_cycle = t_cyc_f)
      else
        msd_fit_result(tau_fit_v, msd_fit_v, k, anomalous, ndim = ndim)
      if (!is.null(fr)) {
        fit_y    <- msd_scale * (if (norm_tau) norm_y(fr$msd_fit, fr$tau_fit) else fr$msd_fit)
        coef_lbl <- sprintf("%dD", 2L * ndim)
        lbl <- if (!use_olsf && anomalous && !fr$forced_linear)
          sprintf("Fit: %s\u03C4\u1D45  (\u03B1=%.2f)", coef_lbl, fr$alpha)
        else if (use_olsf)
          sprintf("OLSF fit: %s\u03C4 (Brownian, opt. weighted)", coef_lbl)
        else sprintf("Fit: %s\u03C4 (linear)", coef_lbl)
        plt <- plt %>%
          add_trace(x = fr$tau_fit, y = fit_y,
                    type = "scattergl", mode = "lines",
                    line = list(color = "#E03020", width = 2),
                    name = lbl)
      }
    }

    # 5. Single-TID overlay
    if (isTRUE(input$trk_msd_show_single)) {
      one <- tryCatch(msd_single(), error = function(e) NULL)
      if (!is.null(one)) {
        ok1 <- is.finite(one$tau) & is.finite(one$msd)
        if (any(ok1)) {
          one_y <- msd_scale * (if (norm_tau) norm_y(one$msd[ok1], one$tau[ok1]) else one$msd[ok1])
          plt <- plt %>%
            add_trace(x = one$tau[ok1], y = one_y,
                      type = "scattergl", mode = "lines+markers",
                      name = paste0("TID ", input$trk_msd_tid),
                      line   = list(color = "#2CA02C", width = 2, dash = "dot"),
                      marker = list(color = "#2CA02C", size = 6))
        }
      }
    }

    # 6. N(tau) — number of tracks contributing at each lag (right y-axis)
    has_ntau <- show_ntau && !is.null(mb$n_tracks_per_bin) &&
                any(is.finite(mb$tau) & mb$n_tracks_per_bin > 0L)
    if (has_ntau) {
      ok_n <- is.finite(mb$tau) & mb$n_tracks_per_bin > 0L
      plt <- plt %>% add_trace(
        x = mb$tau[ok_n], y = mb$n_tracks_per_bin[ok_n],
        type = "scatter", mode = "lines",
        line = list(color = "rgba(180,100,40,0.70)", width = 1.5, dash = "dot"),
        name = "N(\u03C4) tracks",
        yaxis = "y2"
      )
    }

    # Lower-regime shading: shade \u03C4 < min_lag_ms to flag the noise-dominated zone
    min_lag_ms_plot <- min_lag_ms_f
    lower_regime_shapes <- if (min_lag_ms_plot > 0 && min_lag_ms_plot < max(tau)) {
      list(list(
        type = "rect", xref = "x", yref = "paper",
        x0 = 0, x1 = min_lag_ms_plot, y0 = 0, y1 = 1,
        fillcolor = "rgba(220,160,80,0.10)",
        line = list(color = "rgba(200,140,60,0.40)", width = 1, dash = "dot")
      ))
    } else list()

    # Layout (yaxis2 added only when N(tau) trace is present)
    lyt <- list(
      xaxis      = list(title = "\u03C4 (ms)"),
      yaxis      = list(title = y_label, rangemode = "tozero"),
      showlegend = TRUE,
      shapes     = lower_regime_shapes
    )
    if (has_ntau) {
      lyt$yaxis2 <- list(title = "N(\u03C4) tracks", overlaying = "y",
                         side = "right", showgrid = FALSE, rangemode = "tozero",
                         tickfont = list(color = "rgba(180,100,40,0.85)"),
                         titlefont = list(color = "rgba(180,100,40,0.85)"))
    }
    # Annotation explaining the lower-regime shading
    if (min_lag_ms_plot > 0) {
      lyt$annotations <- list(list(
        x = min_lag_ms_plot * 0.5, y = 0.97, xref = "x", yref = "paper",
        text = "lower\nregime", showarrow = FALSE,
        font = list(size = 9, color = "rgba(180,120,40,0.70)"),
        xanchor = "center", yanchor = "top"
      ))
    }

    do.call(plotly::layout, c(list(plt), lyt)) %>%
      config(displaylogo = FALSE,
             toImageButtonOptions = list(format = "svg", filename = "msd_plot"))
  })

  output$trk_msd_estimates <- renderUI({
    mb <- msd_bundle()
    shiny::validate(shiny::need(length(mb$tau) > 0 && length(mb$msd) > 0, NULL))

    ok  <- is.finite(mb$tau) & is.finite(mb$msd)
    tau <- mb$tau[ok]; msd <- mb$msd[ok]
    shiny::validate(shiny::need(length(tau) >= 3, NULL))

    fit_pct    <- 30                      # internal default (anomalous path only)
    min_lag_ms <- msd_min_lag_ms()
    fit_max_ms <- msd_fit_max_ms()         # Vogler fit_limits upper bound (lag 50)
    t_cyc      <- trk_cycle_time()
    anomalous  <- isTRUE(input$trk_msd_anomalous)
    use_olsf   <- !anomalous               # OLSF always, unless anomalous α requested
    auto_pmin  <- use_olsf                 # Michalet optimal point count always on
    R_blur     <- 0                        # motion blur off (MINFLUX)
    ndim_est   <- mb$ndim %||% 2L

    # Restrict to the fit window [min_lag, fit_max] (Vogler's fit_limits), then
    # OLSF picks the optimal point count within it.
    keep <- rep(TRUE, length(tau))
    if (min_lag_ms > 0) keep <- keep & tau >= min_lag_ms
    if (fit_max_ms > 0) keep <- keep & tau <= fit_max_ms
    if (!use_olsf) keep <- keep & msd_fit_mask(tau, min_lag_ms, fit_max_ms, fit_pct)
    tau   <- tau[keep]; msd <- msd[keep]
    npr   <- (mb$npairs %||% rep(1L, length(mb$msd)))[ok][keep]
    shiny::validate(shiny::need(length(tau) >= 3,
      "Fewer than 3 \u03c4 points in the fit window \u2014 widen the lag bounds."))
    k <- length(tau)

    fr <- if (use_olsf)
      msd_fit_olsf(tau, msd, npr, k, ndim = ndim_est,
                   auto_pmin = auto_pmin, R = R_blur, Npts = 1e4, dt_cycle = t_cyc)
    else
      msd_fit_result(tau, msd, k, anomalous, ndim = ndim_est)
    if (is.null(fr)) return(htmltools::div("Fit did not converge \u2014 try a wider fit range or more tau bins."))

    # D conversion: 1 nm\u00b2/ms = 0.001 \u00b5m\u00b2/s
    D_um2_s   <- fr$D / 1e3
    # Cycle-equivalent context \u2014 uses the ACTUAL bin width used by the bundle
    t_c_ctx   <- suppressWarnings(as.numeric(isolate(input$trk_msd_cycle_ms)))
    tau_bin_v <- mb$d_tau %||% NA_real_
    cycle_ctx <- if (is.finite(t_c_ctx) && t_c_ctx > 0 && is.finite(tau_bin_v) && tau_bin_v > 0) {
      bins_per_cyc <- round(tau_bin_v / t_c_ctx, 2)
      sprintf("\u03c4 bin = %.4g ms = %.3g\u00d7 cycle%s  |  fit: %.4g\u2013%.4g ms (%.1f\u2013%.1f cycles)",
              tau_bin_v, bins_per_cyc,
              if (isTRUE(mb$tau_bin_auto)) " (auto)" else "",
              tau[1], tau[length(tau)],
              round(tau[1]/t_c_ctx, 1), round(tau[length(tau)]/t_c_ctx, 1))
    } else NULL

    n_trk_alpha <- if (!is.null(mb$MSD_matrix)) ncol(mb$MSD_matrix) else NA_integer_
    alpha_line <- if (!use_olsf && anomalous) {
      label <- if (fr$forced_linear)
        "1.00 (anomalous fit failed \u2014 fell back to linear)"
      else {
        ci_txt <- if (is.finite(fr$alpha_se %||% NA_real_))
          sprintf(" (95%% CI %.2f\u2013%.2f)",
                  fr$alpha - 1.96 * fr$alpha_se, fr$alpha + 1.96 * fr$alpha_se) else ""
        sprintf("%.3f%s%s", fr$alpha, ci_txt,
          if (fr$alpha < 0.9) "  \u2192 subdiffusion / confinement"
          else if (fr$alpha > 1.1) "  \u2192 directed / superdiffusion"
          else "  \u2192 normal Brownian")
      }
      # Honesty caveat: alpha is an ENSEMBLE estimate; single-track / short-track alpha
      # is unreliable (Kepten et al. 2015), and the CI above is a lower bound (MSD lags
      # correlated). Irregular MINFLUX \u0394t further limits usable lag bins.
      caveat <- htmltools::tags$span(style = "color:#888; font-size:0.85em;",
        sprintf("ensemble over %s tracks; single-track \u03B1 is unreliable at short lengths (Kepten 2015). CI is a lower bound.",
                if (is.finite(n_trk_alpha)) format(n_trk_alpha, big.mark = ",") else "?"))
      htmltools::tagList(
        htmltools::tags$b("\u03B1 (anomalous exponent): "), label, htmltools::br(),
        caveat, htmltools::br()
      )
    }

    # D-ceiling check (shared helper with the D_max panel): warn when the fitted D
    # approaches the trackable limit, where MINFLUX under-estimates D (Vogler 2025 /
    # SimuFLUX). Isolated so it doesn't add a reactive dependency.
    dmax_warn <- {
      de2 <- tryCatch(dmax_estimate(isolate(trk()),
                       suppressWarnings(as.numeric(isolate(input$trk_dmax_L))),
                       suppressWarnings(as.numeric(isolate(input$trk_dmax_t_hw)))),
                      error = function(e) NULL)
      if (!is.null(de2) && is.finite(de2$D_max) && de2$D_max > 0 && is.finite(D_um2_s) &&
          (D_um2_s / de2$D_max) > 0.5) {
        rr <- D_um2_s / de2$D_max
        htmltools::div(style = sprintf("color:%s; font-size:0.85em; margin-top:3px;",
                                        if (rr > 0.8) "#CC2200" else "#BB7700"),
          sprintf("\u26A0 D is %.0f%% of the trackable ceiling (D_max \u2248 %.3g \u00B5m\u00B2/s) \u2014 above it MINFLUX under-estimates D (Vogler 2025; SimuFLUX). See the D_max panel.",
                  100 * rr, de2$D_max))
      } else NULL
    }

    htmltools::div(
      style = "color:#444;",
      htmltools::tags$b("D (diffusion coefficient): "),
        sprintf("%.4g \u00B5m\u00B2/s", D_um2_s),
        if (use_olsf) htmltools::tags$span(style="color:#888;font-size:0.9em;", " (OLSF, Brownian)"),
      htmltools::br(),
      dmax_warn,
      if (anomalous) alpha_line,
      htmltools::tags$b("\u03C3_loc (localization error): "),
        sprintf("%.2f nm", fr$sigma_loc),
        htmltools::tags$span(
          style = "color:#888; font-size:0.88em;",
          " \u2014 how wrong each position estimate can be (smaller = better)"
        ),
        htmltools::br(),
        htmltools::tags$span(
          style = "color:#2070a0; font-size:0.88em;",
          sprintf("[\u2192 enter %.2f nm in the JDD \u03C3 correction field  |  MSD y-intercept \u2248 %.1f nm\u00B2]",
                  fr$sigma_loc, (2L * ndim_est) * fr$sigma_loc^2)
        ),
      htmltools::br(),
      htmltools::tags$b("MSD formula: "),
        sprintf("%dD diffusion \u2014 MSD = %d\u00B7D\u00B7\u03C4",
                ndim_est, 2L * ndim_est),
        htmltools::tags$span(style = "color:#888; font-size:0.9em;",
          if (ndim_est == 3L && !identical(isolate(input$trk_msd_ndim %||% "auto"), "3"))
            " (auto-detected from z data)"
          else if (ndim_est == 2L && !identical(isolate(input$trk_msd_ndim %||% "auto"), "2"))
            " (auto-detected: no z variation)"
          else " (user-selected)"
        ),
      htmltools::br(),
      htmltools::tags$b("Fit range: "),
        if (auto_pmin)
          sprintf("Michalet optimal: %d points from %.4g ms%s",
                  fr$n_points %||% k, tau[1],
                  if (min_lag_ms > 0) " (lower regime skipped)" else "")
        else if (fit_max_ms > 0)
          sprintf("%.4g \u2013 %.4g ms absolute (%d points)", tau[1], tau[length(tau)], k)
        else if (min_lag_ms > 0)
          sprintf("%.4g ms \u2013 %d%% of \u03C4 (%d points, lower regime skipped)",
                  min_lag_ms, as.integer(fit_pct), k)
        else
          sprintf("first %d%% of \u03C4 (%d points)", as.integer(fit_pct), k),
      htmltools::br(),
      # Cycle-equivalent context \u2014 explains \u03C3_loc vs the Advanced (cycle) MSD
      if (!is.null(cycle_ctx)) htmltools::tagList(
        htmltools::tags$span(style = "color:#888; font-size:0.85em;", cycle_ctx),
        htmltools::br(),
        if (is.finite(tau_bin_v) && is.finite(t_c_ctx) && tau_bin_v > 2 * t_c_ctx)
          htmltools::tags$span(
            style = "color:#c05020; font-size:0.85em;",
            sprintf("\u26A0 \u03C4 bin (%.4g ms) is much wider than the cycle time (%.4g ms): the first bin pools many noise-dominated short-lag pairs, inflating \u03C3_loc. Click 'Sync \u03C4-bin & fit-lag to cycle time' for a cleaner intercept.",
                    tau_bin_v, t_c_ctx)
          )
        else NULL,
        htmltools::br()
      ),
      # ---- Ergodicity ---------------------------------------------------------
      {
        ergo_row <- NULL
        MSD_mat  <- mb$MSD_matrix
        tau_all  <- mb$tau
        if (!is.null(MSD_mat) && ncol(MSD_mat) > 1 &&
            is.finite(D_um2_s) && D_um2_s > 0) {
          # Same fit window as the ensemble fit, for a consistent comparison
          ok_e    <- msd_fit_mask(tau_all, min_lag_ms, fit_max_ms, fit_pct)
          tau_e   <- tau_all[ok_e]
          MSD_e   <- MSD_mat[ok_e, , drop = FALSE]
          coef_v  <- 2L * ndim_est
          per_D <- apply(MSD_e, 2L, function(mv) {
            ok_c <- is.finite(mv)
            if (sum(ok_c) < 3L) return(NA_real_)
            cf <- tryCatch(coef(lm(mv[ok_c] ~ tau_e[ok_c])),
                           error = function(e) NULL)
            if (is.null(cf) || !is.finite(cf[2]) || cf[2] < 0) return(NA_real_)
            cf[2] / coef_v / 1e3
          })
          mean_ta_D <- mean(per_D, na.rm = TRUE)
          if (is.finite(mean_ta_D) && mean_ta_D > 0) {
            ergo <- mean_ta_D / D_um2_s
            ergo_row <- htmltools::tagList(
              htmltools::tags$b("Ergodicity (\u27E8D\u27E9_TA / D_EA): "),
              sprintf("%.3f", ergo),
              htmltools::tags$span(
                style = "color:#888; font-size:0.88em;",
                if (abs(ergo - 1) < 0.15)
                  "  \u2014 consistent with ergodic diffusion"
                else if (ergo < 0.85)
                  "  \u2014 time-averaged D < ensemble D (heterogeneous / confined sub-population present)"
                else
                  "  \u2014 time-averaged D > ensemble D (a few fast tracks dominate the ensemble)"
              ),
              htmltools::br()
            )
          }
        }
        ergo_row
      },
      htmltools::tags$b("n tracks: "), mb$n_tracks,
      if (!is.null(mb$n_dropped_negD) && mb$n_dropped_negD > 0)
        sprintf("  (\u2212%d non-diffusing, D\u22640)", mb$n_dropped_negD),
      "  |  median \u0394t: ", fmt_num(mb$dt0), " ms",
      htmltools::br(),
      # Track-length percentiles — helps calibrate tau_max and min_pts
      {
        df_pts <- trk()
        if (!is.null(df_pts) && nrow(df_pts) > 0) {
          n_pts <- as.integer(table(df_pts$tid))
          q <- quantile(n_pts, probs = c(0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
          htmltools::tagList(
            htmltools::tags$b("Track length (pts):"),
            sprintf(" median=%d  |  75th pct=%d  |  90th pct=%d",
                    as.integer(q[2]), as.integer(q[3]), as.integer(q[4])),
            htmltools::br(),
            htmltools::tags$span(style = "color:#888; font-size:0.9em;",
              "(Tip: set \u03C4 max \u2264 median track duration to keep N(\u03C4) high. ",
              "See histogram \u2192 'Trace ID length' for full distribution.)")
          )
        }
      }
    )
  })
  # Shared trackable-diffusion ceiling (D_max) estimate — used by the D_max panel AND
  # surfaced as a warning next to the MSD/eDIFF result (single source of truth). Above
  # ~D_max the MINFLUX localization loop can no longer keep the doughnut centred, so D
  # is systematically UNDER-estimated (Vogler et al. 2025; SimuFLUX, Marin & Ries 2025).
  # Pure function of the localization table + TCP diameter L (nm) & hardware dead-time
  # t_hw (ms). Returns D_max in um^2/s (Vogler form preferred, simple form fallback).
  dmax_estimate <- function(df, L = 400, t_hw = 0.05) {
    if (is.null(df) || !nrow(df) || !all(c("tid","tim") %in% names(df))) return(NULL)
    if (!is.finite(L)    || L    <= 0) L    <- 400
    if (!is.finite(t_hw) || t_hw <  0) t_hw <- 0.05
    efo_med <- if ("efo" %in% names(df)) { v <- df$efo[is.finite(df$efo) & df$efo > 0]; if (length(v)) median(v) else NA_real_ } else NA_real_
    eco_med <- if ("eco" %in% names(df)) { v <- df$eco[is.finite(df$eco) & df$eco > 0]; if (length(v)) median(v) else NA_real_ } else NA_real_
    n_pl <- eco_med   # ECO is the photon count per localization (= N_PL directly)
    t_loc_eco <- if (is.finite(eco_med) && is.finite(efo_med) && efo_med > 0)
      (eco_med / efo_med) * 1000 else NA_real_   # ms to collect ECO photons at rate EFO
    # Median within-track dt fallback (fast integer-coded, not slow character tapply).
    tf  <- match(df$tid, unique(df$tid))
    ord <- order(tf, df$tim, method = "radix")
    tf2 <- tf[ord]; tim2 <- df$tim[ord]; n <- length(tf2)
    same <- c(FALSE, tf2[-1] == tf2[-n])
    dtd  <- c(NA_real_, diff(tim2))[same]; dtd <- dtd[is.finite(dtd) & dtd > 0]
    t_loc_dt <- if (length(dtd)) median(dtd) else NA_real_
    t_loc <- if (is.finite(t_loc_eco)) t_loc_eco else t_loc_dt
    D_max_simple <- if (is.finite(t_loc) && t_loc > 0) (L^2 / (16 * t_loc)) / 1e3 else NA_real_
    D_max_vogler <- if (is.finite(efo_med) && is.finite(n_pl) && is.finite(t_hw)) {
      psi <- efo_med; t_hw_s <- t_hw / 1000
      ((L^2 / 16) * psi / (n_pl + psi * t_hw_s)) / 1e6
    } else NA_real_
    list(efo_med = efo_med, n_pl = n_pl, t_loc_eco = t_loc_eco, t_loc_dt = t_loc_dt,
         t_loc = t_loc, D_max_simple = D_max_simple, D_max_vogler = D_max_vogler,
         D_max = if (is.finite(D_max_vogler)) D_max_vogler else D_max_simple)
  }

  # --- Diffusion limit (D_max) panel -----------------------------------------
  output$trk_dmax_panel <- renderUI({
    req(isTRUE(input$trk_show_dmax))
    df  <- trk(); req(nrow(df) > 0)
    ub  <- tryCatch(msd_unified_bundle(), error = function(e) NULL)

    L    <- suppressWarnings(as.numeric(input$trk_dmax_L))
    t_hw <- suppressWarnings(as.numeric(input$trk_dmax_t_hw))
    if (!is.finite(L)    || L    <= 0) L    <- 400
    if (!is.finite(t_hw) || t_hw <  0) t_hw <- 0.05

    de <- dmax_estimate(df, L, t_hw)            # shared helper (single source of truth)
    efo_med      <- de$efo_med;      n_pl         <- de$n_pl
    t_loc_eco    <- de$t_loc_eco;    t_loc_dt     <- de$t_loc_dt
    t_loc        <- de$t_loc
    D_max_simple <- de$D_max_simple; D_max_vogler <- de$D_max_vogler

    # --- Measured D from the Unified MSD (eDIFF, ensemble fit) -----------------
    D_meas <- if (!is.null(ub)) {
      fr <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
      if (!is.null(fr) && is.finite(fr$D)) fr$D / 1e3 else NA_real_  # nm2/ms -> um2/s
    } else NA_real_

    # --- Build output ---------------------------------------------------------
    fmt <- function(x, digits = 3)
      if (is.finite(x)) sprintf("%.*g", digits, x) else "—"

    ratio_simple <- if (is.finite(D_max_simple) && D_max_simple > 0 && is.finite(D_meas))
      D_meas / D_max_simple else NA_real_
    ratio_vogler <- if (is.finite(D_max_vogler) && D_max_vogler > 0 && is.finite(D_meas))
      D_meas / D_max_vogler else NA_real_

    flag_col <- function(ratio) {
      if (!is.finite(ratio)) return("#444")
      if (ratio > 0.8)  return("#CC2200")   # red: likely affected
      if (ratio > 0.5)  return("#BB7700")   # amber: caution
      return("#228833")                      # green: safe
    }
    flag_msg <- function(ratio) {
      if (!is.finite(ratio)) return("")
      if (ratio > 0.8)  return(" ⚠ D near/above D_max — the tracking loop loses lock, so D is likely UNDER-estimated (Vogler 2025; SimuFLUX)")
      if (ratio > 0.5)  return(" ⚠ D approaching D_max — D may be under-estimated; interpret with caution")
      return(" ✓ D is well below D_max")
    }

    use_ratio <- if (is.finite(ratio_vogler)) ratio_vogler else ratio_simple
    main_dmax <- if (is.finite(D_max_vogler)) D_max_vogler else D_max_simple

    htmltools::div(
      style = "font-size:0.88em; color:#333; padding:6px 0;",

      htmltools::tags$table(
        style = "width:100%; border-collapse:collapse;",
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "TCP diameter L"),
          htmltools::tags$td(style="padding:3px 6px; font-weight:bold;", sprintf("%g nm", L))
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "t_loc from ECO/EFO"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(t_loc_eco)) sprintf("%.3f ms", t_loc_eco) else "— (needs ECO & EFO)")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median Δt within tracks"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(t_loc_dt)) sprintf("%.3f ms", t_loc_dt) else "—")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median EFO (ψ_MFX)"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(efo_med)) sprintf("%.0f Hz", efo_med) else "— (no EFO column)")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median ECO (N_PL, photons/loc)"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(n_pl)) sprintf("%.1f", n_pl) else "— (no ECO column)")
        ),
        htmltools::tags$tr(style="border-top:1px solid #ddd;",
          htmltools::tags$td(style="padding:4px 6px; color:#666;", "D_max simple (L² / 16·t_loc)"),
          htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_max_simple)) sprintf("%.3g µm²/s", D_max_simple) else "—")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:4px 6px; color:#666;", "D_max Vogler (EFO/N_PL)"),
          htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_max_vogler)) sprintf("%.3g µm²/s", D_max_vogler) else "— (needs ECO & EFO)")
        ),
        htmltools::tags$tr(style="border-top:1px solid #ddd;",
          htmltools::tags$td(style="padding:4px 6px; color:#666;", "Measured D (MSD fit)"),
          htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_meas)) sprintf("%.3g µm²/s", D_meas) else "— (run MSD first)")
        ),
        if (is.finite(use_ratio)) htmltools::tags$tr(
          htmltools::tags$td(style="padding:4px 6px; color:#666;", "D / D_max"),
          htmltools::tags$td(
            style = sprintf("padding:4px 6px; font-weight:bold; color:%s;", flag_col(use_ratio)),
            sprintf("%.4f%s", use_ratio, flag_msg(use_ratio))
          )
        )
      ),

      htmltools::tags$div(
        style = "margin-top:8px; font-size:0.82em; color:#666;",
        "D_max simple = L² / (16 · t_loc). ",
        "D_max Vogler = (L²/16) · ψ_MFX / (N_PL + ψ_MFX · t_hw), from Vogler et al. (2025, Commun. Biol.). ",
        "t_loc = ECO / EFO × 1000 (ms) — the time to collect ECO photons at rate EFO. ",
        "ECO is the photon count per localization (= N_PL directly). ",
        "Median Δt is shown as a cross-check. ",
        htmltools::tags$b("Note:"), " this is the THEORETICAL ceiling (Vogler 2025 report ≈10–11 µm²/s), ",
        "but tracking has only been DEMONSTRATED to ≈2.5 µm²/s (lipids) — so treat D above ~2–3 µm²/s ",
        "as suspect even when the flag is still green."
      )
    )
  })

  # Per-track D histogram (log10 or linear scale)
  output$trk_msd_hist_plot <- renderPlotly({
    req(isTRUE(input$trk_msd_show_hist))
    fits <- msd_pertk_fits()
    shiny::validate(shiny::need(!is.null(fits) && nrow(fits) >= 3L,
      "Need at least 3 successfully fitted tracks. Try lowering Min points per TID."))

    dark    <- isTRUE(isolate(input$trk_dark_mode))
    bg_col  <- if (dark) "#000000" else "#FFFFFF"
    txt_col <- if (dark) "#CCCCCC" else "#333333"
    use_log <- isTRUE(input$trk_msd_hist_log %||% TRUE)

    # Exclude D ≤ 0 (noise-clamped, unphysical); keep only positive fits
    D_vals <- fits$D_um2_s[is.finite(fits$D_um2_s) & fits$D_um2_s > 0]
    shiny::validate(shiny::need(length(D_vals) >= 3L, "No tracks with positive D after fitting."))

    if (use_log) {
      x_data  <- log10(D_vals)
      x_title <- "log10(D) [µm²/s]"
      ann_txt <- sprintf("n = %d tracks  |  median log10(D) = %.2f  |  median D = %.3g µm²/s",
                         length(D_vals), median(x_data), median(D_vals))
    } else {
      x_data  <- D_vals
      x_title <- "D [µm²/s]"
      ann_txt <- sprintf("n = %d tracks  |  median D = %.3g µm²/s",
                         length(D_vals), median(D_vals))
    }

    plot_ly() %>%
      add_trace(x = x_data, type = "histogram", nbinsx = 35,
                marker = list(color = "rgba(48,96,192,0.75)",
                              line = list(color = "rgba(30,60,140,0.9)", width = 0.5))) %>%
      layout(
        showlegend = FALSE,
        plot_bgcolor  = bg_col, paper_bgcolor = bg_col,
        font = list(color = txt_col),
        xaxis = list(title = x_title, color = txt_col,
                     gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
        yaxis = list(title = "Count", color = txt_col,
                     gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
        bargap = 0.05,
        annotations = list(list(
          text = ann_txt,
          xref = "paper", yref = "paper", x = 0.5, y = 1.06,
          xanchor = "center", showarrow = FALSE,
          font = list(size = 11, color = txt_col)
        ))
      ) %>%
      config(displaylogo = FALSE,
             toImageButtonOptions = list(format = "svg", filename = "msd_pertk_hist"))
  })
  output$trk_dl_msd_pertk_csv <- downloadHandler(
    filename = function() paste0("msd_pertk_fits_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) write.csv(msd_pertk_fits(), file, row.names = FALSE)
  )

  # MSD cluster scatter: D (log10 or linear) vs CI, colored by cluster
  output$trk_msd_cluster_plot <- renderPlotly({
    req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
    cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(cl),
      "Clustering not ready. Need at least 3 fitted tracks and k \u2264 n tracks."))
    fits    <- cl$fits
    pal     <- cl$palette
    use_log <- isTRUE(input$trk_msd_hist_log %||% TRUE)

    # Exclude D \u2264 0 as in the histogram
    ok_d <- is.finite(fits$D_um2_s) & fits$D_um2_s > 0
    fits <- fits[ok_d, , drop = FALSE]
    shiny::validate(shiny::need(nrow(fits) >= 3L, "No tracks with positive D."))

    x_data  <- if (use_log) log10(fits$D_um2_s) else fits$D_um2_s
    x_title <- if (use_log) "log10(D) [\u00B5m\u00B2/s]" else "D [\u00B5m\u00B2/s]"
    x_min   <- min(x_data, na.rm = TRUE); x_max <- max(x_data, na.rm = TRUE)

    plt <- plot_ly()
    for (c_lbl in seq_len(cl$n_clusters)) {
      idx <- fits$cluster == c_lbl
      if (!any(idx)) next
      hover_x <- if (use_log)
        paste0("log10(D): ", round(x_data[idx], 3))
      else
        paste0("D: ", round(fits$D_um2_s[idx], 6), " \u00B5m\u00B2/s")
      plt <- plt %>% add_trace(
        x    = x_data[idx],
        y    = fits$conf_index[idx],
        type = "scattergl", mode = "markers",
        marker = list(color = pal[c_lbl], size = 6, opacity = 0.75,
                      line = list(color = "#ffffff", width = 0.4)),
        name = sprintf("Cluster %d (n=%d)", c_lbl, sum(idx)),
        text = paste0("TID: ", fits$tid[idx], "<br>", hover_x,
                      "<br>D: ", round(fits$D_um2_s[idx], 6), " \u00B5m\u00B2/s",
                      "<br>CI: ", round(fits$conf_index[idx], 3)),
        hoverinfo = "text"
      )
    }
    plt %>%
      layout(
        xaxis = list(title = x_title),
        yaxis = list(title = "Confinement index (CI)", type = "log"),
        showlegend = TRUE,
        shapes = list(
          list(type = "line", x0 = x_min, x1 = x_max,
               y0 = 1, y1 = 1, xref = "x", yref = "y",
               line = list(color = "rgba(80,80,80,0.35)", dash = "dash", width = 1))
        ),
        annotations = list(
          list(x = 0.5, y = 1.02, xref = "paper", yref = "paper", showarrow = FALSE,
               text = "CI = 1: free diffusion boundary",
               font = list(size = 10, color = "rgba(100,100,100,0.6)"),
               xanchor = "center", yanchor = "bottom")
        )
      ) %>%
      config(displaylogo = FALSE,
             toImageButtonOptions = list(format = "svg", filename = "msd_cluster_scatter"))
  })

  # ---- ms hint for the standard MSD cycle-based min-lag input -----------------
  # Input is in cycles; show the ms equivalent at the current cycle time.
  output$trk_msd_min_lag_cycle_hint <- renderUI({
    t_c <- trk_cycle_time()
    v   <- suppressWarnings(as.numeric(input$trk_msd_min_lag_cyc))
    if (!is.finite(t_c) || t_c <= 0 || !is.finite(v) || v <= 0) return(NULL)
    tags$p(style = "font-size:0.81em; color:#2070a0; margin:2px 0;",
           sprintf("\u2248 %.4g ms  (at %.4g ms/cycle)", v * t_c, t_c))
  })

  # Sync button: make the standard MSD reproduce the Advanced (cycle) MSD as
  # closely as possible.  Auto-detects the cycle time from the modal \u0394t (same
  # logic as the Advanced "Auto" button), then aligns every parameter that
  # affects the result: \u03c4-bin, fit-from-lag, \u03c4-max, weighting mode, and OLSF.
  observeEvent(input$trk_msd_sync_cycle, {
    # 1. Cycle time from the single source of truth (hardware \u2192 manual \u2192 modal \u0394t)
    t_c <- trk_cycle_time()
    if (!is.finite(t_c) || t_c <= 0) {
      showNotification("Could not determine cycle time (no \u0394t data). Set it manually or via the hardware model.",
                       type = "warning", duration = 6)
      return()
    }

    max_lag  <- suppressWarnings(as.integer(input$trk_msd_cycle_max_lag))
    if (!is.finite(max_lag) || max_lag < 5L) max_lag <- 50L
    min_lag  <- suppressWarnings(as.integer(input$trk_msd_cycle_min_lag))
    if (!is.finite(min_lag) || min_lag < 1L) min_lag <- 4L

    # 2. Push the same cycle time to the Advanced panel so both use one value.
    updateNumericInput(session, "trk_msd_cycle_ms",  value = round(t_c, 4))
    # 3. Align standard-MSD binning & fit range (fit window in cycles).
    updateNumericInput(session, "trk_msd_tau_bin",     value = round(t_c, 4))
    updateNumericInput(session, "trk_msd_tau_max",     value = round(max_lag * t_c, 4))
    updateNumericInput(session, "trk_msd_min_lag_cyc", value = min_lag)

    showNotification(
      sprintf(paste0("Synced to cycle = %.4g ms:  \u03c4-bin = 1 cycle  |  \u03c4-max = %d cycles  |  ",
                     "fit from lag %d (%.4g ms).  Both use OLSF + pair-weighting. ",
                     "Click 'Plot MSD' to recompute."),
              t_c, max_lag, min_lag, min_lag * t_c),
      type = "message", duration = 9
    )
  }, ignoreInit = TRUE)

  # ---- Advanced MINFLUX MSD: cycle-based (Vogler et al. 2025) -----------------
  # Resolved cycle-time status shown inside the Advanced MSD panel
  output$trk_msd_cycle_status <- renderUI({
    req(isTRUE(input$trk_msd_advanced))
    tc  <- trk_cycle_time()
    if (!is.finite(tc) || tc <= 0)
      return(tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
                    "Cycle time: not yet determined."))
    hw  <- trk_hw_cycle_time()
    man <- suppressWarnings(as.numeric(input$trk_msd_cycle_ms))
    src <- if (is.finite(hw) && hw > 0) "hardware model"
           else if (is.finite(man) && man > 0) "manual override"
           else "modal \u0394t (data)"
    tags$p(style = "font-size:0.83em; color:#2070a0; margin:2px 0;",
           sprintf("Cycle time in use: %.4g ms (source: %s)", tc, src))
  })

  # Cycle-based MSD (Vogler et al. 2025): round each lag to the nearest integer
  # multiple of t_cycle, then compute per-track TAMSDs and the track-unweighted
  # ensemble mean (each track contributes equally, matching Vogler's approach).
  msd_cycle_bundle <- eventReactive(input$trk_msd_compute, {
    req(isTRUE(input$trk_show_msd), isolate(isTRUE(input$trk_msd_advanced)))

    t_c     <- isolate(trk_cycle_time())   # single source of truth
    cc_sw   <- isolate({ cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
                         if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0 })
    max_lag <- isolate(max(5L, as.integer(input$trk_msd_cycle_max_lag %||% 50L)))
    min_pts <- isolate(max(5L, as.integer(input$trk_msd_min_pts)))
    shiny::validate(shiny::need(is.finite(t_c) && t_c > 0, "Could not determine cycle time."))

    df <- isolate(trk()); req(nrow(df) > 0)
    has_z <- "z" %in% names(df) && {
      zv <- df$z[is.finite(df$z)]; length(zv) > 10 && sd(zv) > 1
    }
    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim) &
               if (has_z) is.finite(df$z) else TRUE, ]
    df <- df[order(df$tid, df$tim), ]
    sp <- split(df, df$tid)

    # Per-track TAMSD matrix: rows=lags, cols=tracks
    # Each cell = time-averaged r\u00b2 for that track at that lag
    n_tracks  <- length(sp)
    tk_acc    <- matrix(0,  max_lag, n_tracks)
    tk_np     <- matrix(0L, max_lag, n_tracks)
    tk_names  <- names(sp)

    withProgress(message = "Computing cycle MSD (per-track)\u2026", value = 0, {
      for (k in seq_along(sp)) {
        incProgress(1 / n_tracks)
        tb <- sp[[k]]
        n  <- nrow(tb)
        if (n < min_pts) next
        # Cumulative cycle index (Bela's TOT): round EACH interval to an integer
        # number of cycles, then cumsum. The pair lag = cycidx[j] - cycidx[i].
        # Per-interval rounding is robust to a small cycle-time error, so the
        # integer lag never drifts the way round(total/t_c) does (which skips and
        # doubles bins → the sawtooth).
        dcyc   <- pmax(1L, as.integer(round((diff(tb$tim) - cc_sw) / t_c)))
        cycidx <- c(0L, cumsum(dcyc))
        for (i in seq_len(n - 1L)) {
          lags  <- cycidx[seq.int(i + 1L, n)] - cycidx[i]
          ok    <- lags >= 1L & lags <= max_lag
          if (!any(ok)) next
          js2   <- seq.int(i + 1L, n)[ok]
          b2    <- lags[ok]
          dr2   <- (tb$x[js2] - tb$x[i])^2 + (tb$y[js2] - tb$y[i])^2
          if (has_z) dr2 <- dr2 + (tb$z[js2] - tb$z[i])^2
          rs    <- rowsum(dr2, b2, reorder = FALSE)
          ix    <- as.integer(rownames(rs))
          tk_acc[ix, k] <- tk_acc[ix, k] + rs[, 1L]
          tk_np[, k]    <- tk_np[, k] + tabulate(b2, nbins = max_lag)
        }
      }
    })

    # Per-track TAMSD: NA when no pairs at that lag for that track
    track_msd <- ifelse(tk_np > 0L, tk_acc / tk_np, NA_real_)

    # Drop non-positive-D tracks (Bela's TAMSD D>0 filter) before the ensemble.
    posD <- .track_positive_D(track_msd, seq_len(max_lag) * t_c)
    n_dropped_negD <- if (any(posD) && !all(posD)) sum(!posD) else 0L
    if (any(posD) && !all(posD)) {
      track_msd <- track_msd[, posD, drop = FALSE]
      tk_acc    <- tk_acc[, posD, drop = FALSE]
      tk_np     <- tk_np[, posD, drop = FALSE]
      tk_names  <- tk_names[posD]
    }

    n_trk_lag <- rowSums(!is.na(track_msd))   # N(\u03c4): tracks per lag
    msd_q25   <- apply(track_msd, 1L, quantile, probs = 0.25, na.rm = TRUE)
    msd_q75   <- apply(track_msd, 1L, quantile, probs = 0.75, na.rm = TRUE)

    # Pair count
    np_acc <- rowSums(tk_np)
    tk_acc_total <- rowSums(tk_acc)

    # Always pair-weighted (smoother, robust)
    weight_eq <- FALSE
    msd_ens <- if (weight_eq)
      rowMeans(track_msd, na.rm = TRUE)           # track-unweighted (Vogler)
    else
      ifelse(np_acc > 0L, tk_acc_total / np_acc, NA_real_)  # pair-weighted (smoother)

    # Mask lags where too few tracks contribute (floor: 10% of max or 2)
    min_trk <- max(2L, floor(0.10 * max(n_trk_lag)))
    msd_ens[n_trk_lag  < min_trk] <- NA_real_
    msd_q25[n_trk_lag  < min_trk] <- NA_real_
    msd_q75[n_trk_lag  < min_trk] <- NA_real_

    list(lag         = seq_len(max_lag),
         tau_ms      = seq_len(max_lag) * t_c,
         msd         = msd_ens,
         msd_q25     = msd_q25,
         msd_q75     = msd_q75,
         n_trk_lag   = n_trk_lag,
         track_msd   = track_msd,
         tk_names    = tk_names,
         npairs      = np_acc,
         weight_equal = weight_eq,
         n_dropped_negD = n_dropped_negD,
         t_cycle     = t_c,
         has_z       = has_z,
         ndim        = if (has_z) 3L else 2L)
  })

  output$trk_msd_cycle_plot <- renderPlotly({
    req(isTRUE(input$trk_msd_advanced))
    cb <- msd_cycle_bundle()
    shiny::validate(shiny::need(!is.null(cb), "Click 'Plot MSD' to compute."))

    dark        <- isTRUE(isolate(input$trk_dark_mode))
    bg_col      <- if (dark) "#000000" else "#FFFFFF"
    fg_col      <- if (dark) "#FFFFFF" else "#333333"
    show_iqr    <- isTRUE(input$trk_cycle_show_band)
    show_spag   <- isTRUE(input$trk_cycle_show_spag)
    show_ntau   <- isTRUE(input$trk_cycle_show_ntau)

    min_lag <- max(1L, as.integer(input$trk_msd_cycle_min_lag %||% 4L))
    fit_pct <- 30

    tau_ms  <- cb$tau_ms
    msd_um2 <- cb$msd / 1e6
    ok_data <- is.finite(msd_um2)

    # OLSF on upper regime — auto-p_min picks the optimal point count from min_lag
    auto_pmin_c <- TRUE
    R_blur_c    <- 0
    ok_fit  <- ok_data & cb$lag >= min_lag
    tau_fit <- tau_ms[ok_fit]; msd_fit_nm2 <- cb$msd[ok_fit]
    k_fit   <- max(3L, floor((fit_pct / 100) * length(tau_fit)))
    fr <- if (length(tau_fit) >= 3L)
      tryCatch(msd_fit_olsf(tau_fit, msd_fit_nm2,
                            cb$npairs[ok_fit], k_fit, ndim = cb$ndim,
                            auto_pmin = auto_pmin_c, R = R_blur_c, Npts = 1e4,
                            dt_cycle = cb$t_cycle),
               error = function(e) NULL)
    else NULL

    fig <- plot_ly()

    # Individual track TAMSDs (spaghetti) \u2014 show up to 80 tracks
    if (show_spag && !is.null(cb$track_msd)) {
      n_show  <- min(80L, ncol(cb$track_msd))
      idx_show <- sample.int(ncol(cb$track_msd), n_show)
      for (j in idx_show) {
        y_j <- cb$track_msd[, j] / 1e6
        ok_j <- is.finite(y_j)
        if (sum(ok_j) < 2L) next
        fig <- fig %>% add_trace(
          x = tau_ms[ok_j], y = y_j[ok_j],
          type = "scatter", mode = "lines",
          line = list(color = "rgba(100,160,220,0.18)", width = 0.8),
          showlegend = FALSE, hoverinfo = "skip"
        )
      }
    }

    # IQR band
    if (show_iqr) {
      q25 <- cb$msd_q25 / 1e6; q75 <- cb$msd_q75 / 1e6
      ok_q <- is.finite(q25) & is.finite(q75)
      if (any(ok_q)) {
        fig <- fig %>%
          add_trace(x = c(tau_ms[ok_q], rev(tau_ms[ok_q])),
                    y = c(q75[ok_q], rev(q25[ok_q])),
                    type = "scatter", mode = "none",
                    fill = "toself",
                    fillcolor = "rgba(70,130,200,0.18)",
                    line = list(color = "transparent"),
                    name = "IQR (track spread)", showlegend = TRUE)
      }
    }

    # Ensemble mean
    fig <- fig %>% add_trace(
      x = tau_ms[ok_data], y = msd_um2[ok_data],
      type = "scatter", mode = "lines+markers",
      marker = list(size = 5, color = "steelblue", opacity = 0.9),
      line   = list(color = "steelblue", width = 2),
      name   = "Ensemble mean"
    )

    # OLSF fit line
    if (!is.null(fr)) {
      coef_v <- 2L * cb$ndim
      t_seq  <- seq(tau_fit[1], tau_fit[min(k_fit, length(tau_fit))], length.out = 200)
      fit_y  <- pmax(0, coef_v * fr$D * t_seq + coef_v * fr$sigma_loc^2) / 1e6
      fig <- fig %>% add_trace(
        x = t_seq, y = fit_y,
        type = "scatter", mode = "lines",
        line = list(color = "#E03020", width = 2.5, dash = "dash"),
        name = sprintf("OLSF  D=%.3g \u00B5m\u00B2/s", fr$D / 1e3)
      )
    }

    # Regime boundary
    regime_ms <- min_lag * cb$t_cycle
    y_top     <- max(msd_um2[ok_data], na.rm = TRUE)
    fig <- fig %>% add_trace(
      x = c(regime_ms, regime_ms), y = c(0, y_top * 1.18),
      type = "scatter", mode = "lines",
      marker = list(size = 0),
      line = list(color = "rgba(180,100,20,0.6)", width = 1.5, dash = "dot"),
      name = sprintf("Regime boundary (lag %d)", min_lag)
    )

    # N(\u03C4): show as a line on secondary axis only when it actually varies.
    # When constant (all tracks contribute at all lags \u2014 the ideal case),
    # show a text annotation instead of a flat, uninformative line.
    ntau_varies  <- show_ntau && diff(range(cb$n_trk_lag, na.rm = TRUE)) > 0
    ntau_const   <- show_ntau && !ntau_varies
    ntau_val     <- if (ntau_const) cb$n_trk_lag[1L] else NA_integer_

    if (ntau_varies) {
      fig <- fig %>% add_trace(
        x = tau_ms, y = cb$n_trk_lag,
        type = "scatter", mode = "lines",
        line = list(color = "rgba(160,160,160,0.5)", width = 1, dash = "dot"),
        yaxis = "y2", name = "N(\u03C4) tracks", showlegend = TRUE
      )
    }

    y2_layout <- if (ntau_varies)
      list(overlaying = "y", side = "right", title = "N(\u03C4) tracks",
           color = "rgba(150,150,150,0.7)", showgrid = FALSE)
    else NULL

    lay <- list(
      xaxis = list(title = "Time lag \u03C4 (ms)", color = fg_col, zeroline = FALSE),
      yaxis = list(title = "MSD (\u00B5m\u00B2)",  color = fg_col, zeroline = FALSE),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col),
      legend = list(x = 0.02, y = 0.98, xanchor = "left", yanchor = "top",
                    font = list(color = fg_col, size = 11)),
      annotations = {
        anns <- list(list(
          x = 0.98, y = 0.02, xref = "paper", yref = "paper", showarrow = FALSE,
          xanchor = "right", yanchor = "bottom",
          text = sprintf("t_cycle = %.4g ms  |  lower regime: lags 1\u2013%d",
                         cb$t_cycle, min_lag - 1L),
          font = list(size = 10, color = "rgba(130,130,130,0.85)")
        ))
        if (ntau_const)
          anns <- c(anns, list(list(
            x = 0.02, y = 0.02, xref = "paper", yref = "paper", showarrow = FALSE,
            xanchor = "left", yanchor = "bottom",
            text = sprintf("N(\u03c4) = %d tracks at all lags (all tracks contribute)", ntau_val),
            font = list(size = 10, color = "rgba(130,130,130,0.75)")
          )))
        anns
      },
      uirevision = "msd_cycle"
    )
    if (!is.null(y2_layout)) lay$yaxis2 <- y2_layout

    do.call(plotly::layout, c(list(fig), lay))
  })

  output$trk_msd_cycle_estimates <- renderUI({
    req(isTRUE(input$trk_msd_advanced))
    cb <- tryCatch(msd_cycle_bundle(), error = function(e) NULL)
    if (is.null(cb)) return(NULL)

    min_lag <- max(1L, as.integer(input$trk_msd_cycle_min_lag %||% 4L))
    fit_pct <- 30
    auto_pmin_c <- TRUE
    R_blur_c    <- 0
    ok_fit <- is.finite(cb$msd) & cb$lag >= min_lag
    tau_f  <- cb$tau_ms[ok_fit]; msd_f <- cb$msd[ok_fit]
    k_fit  <- max(3L, floor((fit_pct / 100) * length(tau_f)))
    if (length(tau_f) < 3L)
      return(htmltools::div(style="color:#888;font-size:0.86em;",
                            "Not enough upper-regime lag points to fit."))

    fr <- tryCatch(msd_fit_olsf(tau_f, msd_f, cb$npairs[ok_fit], k_fit, ndim = cb$ndim,
                                auto_pmin = auto_pmin_c, R = R_blur_c, Npts = 1e4,
                                dt_cycle = cb$t_cycle),
                   error = function(e) NULL)
    if (is.null(fr))
      return(htmltools::div(style="color:#888;font-size:0.86em;", "OLSF fit did not converge."))

    weight_lbl <- if (isTRUE(cb$weight_equal)) "track-weighted" else "pair-weighted"
    npts_lbl   <- if (auto_pmin_c) sprintf("Michalet optimal %d pts", fr$n_points %||% k_fit)
                  else sprintf("%d pts", fr$n_points %||% k_fit)
    htmltools::div(
      style = "font-size:0.87em; color:#444; margin:6px 0;",
      htmltools::tags$b("Cycle-based OLSF (upper regime only):"),
      htmltools::br(),
      sprintf("D = %.4g \u00B5m\u00B2/s  |  \u03C3_loc = %.2f nm  |  cycle time = %.4g ms  |  %s  |  %s",
              fr$D / 1e3, fr$sigma_loc, cb$t_cycle, weight_lbl, npts_lbl),
      if (!is.null(cb$n_dropped_negD) && cb$n_dropped_negD > 0)
        sprintf("  |  \u2212%d non-diffusing (D\u22640) tracks dropped", cb$n_dropped_negD),
      htmltools::br(),
      htmltools::tags$span(
        style = "color:#888; font-size:0.88em;",
        paste0("\u03C3_loc here is the pure-diffusion intercept \u2014 extrapolated from the upper regime only. ",
               "It excludes linker-dynamics and is the best value to use in the JDD correction. ",
               "The standard MSD \u03C3_loc (fitting from lag 1) is higher because it includes ",
               "early-lag linker wobble.")
      ),
      htmltools::br(),
      htmltools::tags$span(
        style = "color:#2070a0; font-size:0.9em;",
        sprintf("[\u2192 enter %.2f nm in the JDD \u03C3 correction field]", fr$sigma_loc)
      )
    )
  })

  # ============================================================================
  # UNIFIED MSD (proposed consolidation) — one correct-by-default estimator.
  #   • Cycle-resolution binning: lag = round(Δt / cycle_time)
  #   • Each bin placed at its TRUE mean Δt → D is independent of the cycle-time
  #     estimate (combines the Advanced cycle binning with the Standard real-time
  #     τ axis).
  #   • Pair-weighted ensemble for D; track-weighted computed only for ergodicity.
  #   • OLSF + Michalet optimal n_points, always on. R = 0 (MINFLUX).
  #   • Only user choices: lower-regime cutoff (cycles) and 2D/3D (auto).
  # ============================================================================
  msd_unified_bundle <- eventReactive(input$trk_msd_compute, {
    req(isTRUE(input$trk_show_msd), isolate(isTRUE(input$trk_msd_unified)))
    t_c <- isolate(trk_cycle_time())
    cc_sw <- isolate({ cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
                       if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0 })
    shiny::validate(shiny::need(is.finite(t_c) && t_c > 0, "Could not determine cycle time."))
    # Number of cycle-lags to compute = the requested upper fit limit (so the
    # confined-model plateau is reachable), floored at 60 and capped at 5000 to
    # avoid runaway. Re-click Plot MSD after changing the upper limit to recompute.
    max_lag <- max(60L, min(5000L, as.integer(isolate(input$trk_uni_max_lag) %||% 50L)))
    min_pts <- isolate(max(5L, as.integer(input$trk_msd_min_pts)))

    df <- isolate(trk()); req(nrow(df) > 0)
    # Dimensionality: honour the user's 2D/3D selector (auto = 3D only if z varies).
    ndim_sel   <- isolate(input$trk_msd_ndim %||% "auto")
    has_z_auto <- "z" %in% names(df) && { zv <- df$z[is.finite(df$z)]; length(zv) > 10 && sd(zv) > 1 }
    has_z <- if (identical(ndim_sel, "2")) FALSE
             else if (identical(ndim_sel, "3")) ("z" %in% names(df))
             else has_z_auto
    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim) &
               if (has_z) is.finite(df$z) else TRUE, ]
    df <- df[order(df$tid, df$tim), ]
    sp <- split(df, df$tid)
    n_tracks <- length(sp)
    tk_names <- names(sp)                      # per-track TIDs (column order)

    tk_dr2 <- matrix(0,  max_lag, n_tracks)   # per-track Σ dr²
    tk_np  <- matrix(0L, max_lag, n_tracks)   # per-track pair counts
    tk_dt  <- matrix(0,  max_lag, n_tracks)   # per-track Σ Δt (for mean-Δt-per-bin)

    withProgress(message = "Computing Unified MSD…", value = 0, {
      for (k in seq_along(sp)) {
        incProgress(1 / n_tracks)
        tb <- sp[[k]]; n <- nrow(tb)
        if (n < min_pts) next
        # Cumulative cycle index (Bela's TOT) for the integer lag — robust to
        # cycle-time error, removing the round(total/t_c) aliasing / sawtooth.
        # The real Δt (dt_v) is still kept for the mean-Δt-per-bin x-axis.
        dcyc   <- pmax(1L, as.integer(round((diff(tb$tim) - cc_sw) / t_c)))
        cycidx <- c(0L, cumsum(dcyc))
        for (i in seq_len(n - 1L)) {
          dt_v <- tb$tim[seq.int(i + 1L, n)] - tb$tim[i]
          lag  <- cycidx[seq.int(i + 1L, n)] - cycidx[i]
          ok   <- lag >= 1L & lag <= max_lag
          if (!any(ok)) next
          js   <- seq.int(i + 1L, n)[ok]; b <- lag[ok]; dtk <- dt_v[ok]
          dr2  <- (tb$x[js] - tb$x[i])^2 + (tb$y[js] - tb$y[i])^2
          if (has_z) dr2 <- dr2 + (tb$z[js] - tb$z[i])^2
          rs <- rowsum(dr2, b, reorder = FALSE); ix <- as.integer(rownames(rs))
          tk_dr2[ix, k] <- tk_dr2[ix, k] + rs[, 1L]
          rt <- rowsum(dtk, b, reorder = FALSE); it <- as.integer(rownames(rt))
          tk_dt[it, k]  <- tk_dt[it, k] + rt[, 1L]
          tk_np[, k]    <- tk_np[, k] + tabulate(b, nbins = max_lag)
        }
      }
    })

    track_msd <- ifelse(tk_np > 0L, tk_dr2 / tk_np, NA_real_)   # per-track TAMSD (nm²)

    # Drop non-positive-D tracks before the ensemble. OFF by default to match
    # Bela's eDIFF, which keeps EVERY track (his D≤0 filter applies only to the
    # per-track tDIFF, which we still enforce there). Toggle on for the previous
    # QC behaviour.
    n_dropped_negD <- 0L
    if (isTRUE(input$trk_uni_drop_negd)) {
      posD <- .track_positive_D(track_msd, seq_len(max_lag) * t_c)
      n_dropped_negD <- if (any(posD) && !all(posD)) sum(!posD) else 0L
      if (any(posD) && !all(posD)) {
        track_msd <- track_msd[, posD, drop = FALSE]
        tk_dr2    <- tk_dr2[, posD, drop = FALSE]
        tk_np     <- tk_np[, posD, drop = FALSE]
        tk_dt     <- tk_dt[, posD, drop = FALSE]
        tk_names  <- tk_names[posD]
      }
    }
    colnames(track_msd) <- tk_names

    np_lag    <- rowSums(tk_np)
    dr2_lag   <- rowSums(tk_dr2)
    dt_lag    <- rowSums(tk_dt)
    # Ensemble weighting: pair-weighted (default, smoother) or track-weighted
    # (Vogler/Bela eMSD = mean of per-track TAMSDs).
    track_weight <- isolate(isTRUE(input$trk_uni_track_weight))
    msd_pair  <- if (track_weight) rowMeans(track_msd, na.rm = TRUE)
                 else ifelse(np_lag > 0, dr2_lag / np_lag, NA_real_)
    tau_real  <- ifelse(np_lag > 0, dt_lag  / np_lag, seq_len(max_lag) * t_c)  # true mean Δt
    n_trk_lag <- rowSums(!is.na(track_msd))
    msd_q25   <- apply(track_msd, 1L, quantile, probs = 0.25, na.rm = TRUE)
    msd_q75   <- apply(track_msd, 1L, quantile, probs = 0.75, na.rm = TRUE)

    min_trk <- max(2L, floor(0.10 * max(n_trk_lag)))
    keep_lag <- n_trk_lag >= min_trk
    list(lag = seq_len(max_lag), tau = tau_real, tau_lin = seq_len(max_lag) * t_c,
         msd = msd_pair,
         npairs = np_lag, n_trk_lag = n_trk_lag, keep_lag = keep_lag,
         track_msd = track_msd, tk_names = tk_names,
         msd_q25 = msd_q25, msd_q75 = msd_q75,
         t_cycle = t_c, has_z = has_z, ndim = if (has_z) 3L else 2L,
         track_weighted = track_weight, n_dropped_negD = n_dropped_negD,
         n_tracks = sum(colSums(tk_np) > 0))
  })

  # Fit helper shared by the unified plot + estimates (OLSF + Michalet p_min)
  msd_unified_fit <- function(ub) {
    min_lag <- max(1L, as.integer(isolate(input$trk_uni_min_lag) %||% 4L))
    # Upper fit limit (Vogler fit_limits, default 50). Read WITHOUT isolate so the
    # plot/estimates that call this re-render when the user changes it.
    max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    tau_use <- if (isTRUE(isolate(input$trk_uni_linear_tau))) ub$tau_lin else ub$tau
    ok  <- is.finite(tau_use) & is.finite(ub$msd) & ub$keep_lag &
           ub$lag >= min_lag & ub$lag <= max_lag
    if (sum(ok) < 3L) return(NULL)
    if (isTRUE(input$trk_uni_iterative)) {
      # Bela's exact iterative loop — converges on the lowest, least-noisy lags.
      msd_fit_olsf_iterative(tau_use[ok], ub$msd[ok], ndim = ub$ndim, R = 0,
                             dt_cycle = ub$t_cycle, N = sum(ok))
    } else {
      # Robust full-window-anchored Michalet (default).
      msd_fit_olsf(tau_use[ok], ub$msd[ok], ub$npairs[ok],
                   k = sum(ok), ndim = ub$ndim,
                   auto_pmin = TRUE, R = 0, Npts = 1e4, dt_cycle = ub$t_cycle)
    }
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
  .msd_model_fit <- function(s, y, npairs, ndim, model, use_iter, t_cycle) {
    coef_v <- 2L * ndim
    if (identical(model, "linear")) {
      if (length(s) < 3L) return(NULL)
      return(if (isTRUE(use_iter))
               msd_fit_olsf_iterative(s, y, ndim = ndim, R = 0, dt_cycle = t_cycle, N = length(s))
             else
               msd_fit_olsf(s, y, npairs, k = length(s), ndim = ndim,
                            auto_pmin = TRUE, R = 0, Npts = 1e4, dt_cycle = t_cycle))
    }
    if (length(s) < 4L) return(NULL)
    nls_fit <- function(formula, start, lower, upper) tryCatch(
      stats::nls(formula, data = data.frame(s = s, y = y), start = start,
                 lower = lower, upper = upper, algorithm = "port",
                 control = stats::nls.control(maxiter = 200, warnOnly = TRUE)),
      error = function(e) NULL)
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
      degenerate <- (tc > 2 * max(s)) || (A > 3 * max(y, na.rm = TRUE))
      list(model = "confined", D = A / (coef_v * tc), sigma_loc = sqrt(max(0, C / coef_v)),
           L_conf = sqrt(3 * A), tau_c_ms = tc, plateau_um2 = A / 1e6, degenerate = degenerate,
           tau_fit = tseq, msd_fit = A * (1 - exp(-tseq / tc)) + C, ndim = ndim, method = "confined")
    } else if (identical(model, "directed")) {
      fit <- nls_fit(y ~ coef_v * D * s + (v * s)^2 + C,
                     start = list(D = D0, v = sqrt(D0) / max(s), C = C0),
                     lower = list(D = 0, v = 0, C = 0),
                     upper = list(D = Inf, v = Inf, C = Inf))
      if (is.null(fit)) return(NULL)
      cf <- coef(fit); D <- unname(cf["D"]); v <- unname(cf["v"]); C <- unname(cf["C"])
      list(model = "directed", D = D, sigma_loc = sqrt(max(0, C / coef_v)), v_um_s = v,
           tau_fit = tseq, msd_fit = coef_v * D * tseq + (v * tseq)^2 + C, ndim = ndim, method = "directed")
    } else if (identical(model, "anomalous")) {
      fit <- nls_fit(y ~ coef_v * D * s^a + C,
                     start = list(D = D0, a = 1, C = C0),
                     lower = list(D = 0, a = 0.05, C = 0),
                     upper = list(D = Inf, a = 2.5, C = Inf))
      if (is.null(fit)) return(NULL)
      cf <- coef(fit); D <- unname(cf["D"]); a <- unname(cf["a"]); C <- unname(cf["C"])
      list(model = "anomalous", D = D, alpha = a, sigma_loc = sqrt(max(0, C / coef_v)),
           tau_fit = tseq, msd_fit = coef_v * D * tseq^a + C, ndim = ndim, method = "anomalous")
    } else NULL
  }

  # Ensemble-MSD model fit for the unified plot (selectable model). Non-linear
  # models fit the whole [min_lag, max_lag] window (they need curvature/plateau).
  msd_unified_modelfit <- function(ub) {
    model <- input$trk_uni_model %||% "linear"
    if (identical(model, "linear")) return(msd_unified_fit(ub))
    min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
    max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    tau_use <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
    ok <- is.finite(tau_use) & is.finite(ub$msd) & ub$keep_lag &
          ub$lag >= min_lag & ub$lag <= max_lag
    if (sum(ok) < 4L) return(NULL)
    .msd_model_fit(tau_use[ok], ub$msd[ok], ub$npairs[ok], ub$ndim, model,
                   isTRUE(input$trk_uni_iterative), ub$t_cycle)
  }

  # Cached per-track analysis from the UNIFIED MSD — the single source for the
  # time-averaged D (tDIFF), the per-track D histogram, the D/CI clustering and the
  # per-track CSV export. Runs the same per-track OLSF as the ensemble (iterative
  # when toggled, else anchored) on every track's TAMSD over the fit window, and
  # computes a confinement index (CI) from the same windowed curve. Cached so the
  # heavy loop only re-runs when the bundle or the window/method changes.
  # Returns list(tDIFF µm²/s, n_tD, per_D, pertrack = data.frame).
  trk_unified_tdiff <- reactive({
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
    if (is.null(ub) || is.null(ub$track_msd) || ncol(ub$track_msd) == 0L) return(NULL)
    min_lag   <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
    max_lag_e <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    use_iter  <- isTRUE(input$trk_uni_iterative)
    tau_sel   <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
    win   <- ub$lag >= min_lag & ub$lag <= max_lag_e & is.finite(tau_sel)
    tau_w <- tau_sel[win]; t_c_e <- ub$t_cycle
    tids  <- if (!is.null(ub$tk_names)) ub$tk_names else as.character(seq_len(ncol(ub$track_msd)))
    rows <- lapply(seq_len(ncol(ub$track_msd)), function(j) {
      col <- ub$track_msd[, j]
      mv  <- col[win]
      okc <- is.finite(mv); if (sum(okc) < 3L) return(NULL)
      tt <- tau_w[okc]; yy <- mv[okc]; N_tr <- sum(is.finite(col))
      fr_tr <- if (use_iter)
        msd_fit_olsf_iterative(tt, yy, ndim = ub$ndim, R = 0, dt_cycle = t_c_e, N = N_tr)
      else
        msd_fit_olsf(tt, yy, rep(1L, length(tt)), k = length(tt), ndim = ub$ndim,
                     auto_pmin = TRUE, R = 0, Npts = N_tr, dt_cycle = t_c_e)
      if (is.null(fr_tr) || !is.finite(fr_tr$D) || fr_tr$D <= 0) return(NULL)
      # Confinement index: actual MSD growth vs the free-diffusion expectation,
      # using the median of the first/last 15% of the windowed bins (robust to a
      # noisy endpoint). CI < 1 confined, ~1 Brownian, > 1 directed.
      n_ci <- length(tt); n_edge <- max(2L, floor(0.15 * n_ci))
      tf <- median(tt[seq_len(n_edge)]); tl <- median(tt[seq(n_ci - n_edge + 1L, n_ci)])
      mf <- median(yy[seq_len(n_edge)]); ml <- median(yy[seq(n_ci - n_edge + 1L, n_ci)])
      ci <- if (is.finite(tf) && tf > 0 && is.finite(mf) && mf > 0 &&
                is.finite(tl) && tl > tf && is.finite(ml)) (ml / mf) / (tl / tf) else NA_real_
      data.frame(tid = tids[j], D_um2_s = fr_tr$D / 1e3, alpha = 1,
                 sigma_loc_nm = fr_tr$sigma_loc, conf_index = ci,
                 n_bins = sum(okc), stringsAsFactors = FALSE)
    })
    pertrack <- do.call(rbind, Filter(Negate(is.null), rows))
    per_D <- if (!is.null(pertrack)) pertrack$D_um2_s else numeric(0)
    ci <- .mean_ci(per_D)   # track-level SD + 95% CI on the per-track-mean D
    list(tDIFF = ci$mean, n_tD = ci$n, sd_tD = ci$sd, ci_lo = ci$lo, ci_hi = ci$hi,
         med_tD = ci$median, per_D = per_D, pertrack = pertrack)
  })

  # Per-cluster diffusion summary (cached): eDIFF from the cluster eMSD fit (same
  # window + method as the overall ensemble) and tDIFF from the cluster's per-track
  # D's. Used by BOTH the unified plot legend and the estimates breakdown so the
  # cluster eMSD is fit only once. Returns a data.frame(cluster, eDIFF, tDIFF, n,
  # color) or NULL.
  trk_unified_cluster_diff <- reactive({
    if (!isTRUE(input$trk_msd_cluster)) return(NULL)
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
    cl <- tryCatch(msd_cluster_rv(),    error = function(e) NULL)
    td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
    if (is.null(ub) || is.null(cl) || is.null(ub$track_msd)) return(NULL)
    tids_u    <- if (!is.null(ub$tk_names)) ub$tk_names else as.character(seq_len(ncol(ub$track_msd)))
    cids      <- cl$cluster_ids[tids_u]
    min_lag   <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
    max_lag_c <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    tau_sel   <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
    use_iter  <- isTRUE(input$trk_uni_iterative)
    pt <- if (!is.null(td)) td$pertrack else NULL
    do.call(rbind, lapply(seq_len(cl$n_clusters), function(clk) {
      cols_clk <- which(!is.na(cids) & cids == clk)
      if (!length(cols_clk)) return(NULL)
      emsd <- rowMeans(ub$track_msd[, cols_clk, drop = FALSE], na.rm = TRUE)  # nm²
      winc <- ub$lag >= min_lag & ub$lag <= max_lag_c & is.finite(tau_sel) & is.finite(emsd)
      eD <- if (sum(winc) >= 3L) {
        frc <- tryCatch(
          if (use_iter) msd_fit_olsf_iterative(tau_sel[winc], emsd[winc], ndim = ub$ndim,
                                               R = 0, dt_cycle = ub$t_cycle, N = sum(winc))
          else msd_fit_olsf(tau_sel[winc], emsd[winc], rep(1L, sum(winc)), k = sum(winc),
                            ndim = ub$ndim, auto_pmin = TRUE, R = 0, Npts = 1e4, dt_cycle = ub$t_cycle),
          error = function(e) NULL)
        if (!is.null(frc) && is.finite(frc$D)) frc$D / 1e3 else NA_real_
      } else NA_real_
      tD <- if (!is.null(pt)) {
        dv <- pt$D_um2_s[!is.na(cl$cluster_ids[pt$tid]) & cl$cluster_ids[pt$tid] == clk]
        if (length(dv)) mean(dv, na.rm = TRUE) else NA_real_
      } else NA_real_
      data.frame(cluster = clk, eDIFF = eD, tDIFF = tD, n = length(cols_clk),
                 color = cl$palette[clk], stringsAsFactors = FALSE)
    }))
  })

  output$trk_msd_unified_plot <- renderPlotly({
    req(isTRUE(input$trk_msd_unified))
    ub <- msd_unified_bundle()
    shiny::validate(shiny::need(!is.null(ub), "Click 'Plot MSD' to compute."))
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    show_band <- isTRUE(input$trk_uni_show_band); show_ntau <- isTRUE(input$trk_uni_show_ntau)
    min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
    linear_tau <- isTRUE(input$trk_uni_linear_tau)
    tau_full <- if (linear_tau) ub$tau_lin else ub$tau

    ok <- is.finite(tau_full) & is.finite(ub$msd) & ub$keep_lag
    tau <- tau_full[ok]; msd <- ub$msd[ok] / 1e6   # µm²
    fr  <- msd_unified_modelfit(ub)
    # Regime-cutoff x-position: exact on the linearized axis; the real τ at the
    # cutoff lag on the mean-Δt axis.
    regime_x <- if (linear_tau) min_lag * ub$t_cycle else {
      idx <- which(ub$lag == min_lag)
      if (length(idx) && is.finite(ub$tau[idx])) ub$tau[idx] else min_lag * ub$t_cycle
    }

    fig <- plot_ly()

    # Cluster map (tid -> cluster / colour) when clustering is enabled.
    cl <- if (isTRUE(input$trk_msd_cluster)) tryCatch(msd_cluster_rv(), error = function(e) NULL) else NULL
    tids_u <- if (!is.null(ub$tk_names)) ub$tk_names else as.character(seq_len(ncol(ub$track_msd)))

    # Per-track TAMSD spaghetti (optional, drawn first = behind). Coloured by
    # cluster when clustering is on, else faint grey. Sampled for performance.
    if (isTRUE(input$trk_uni_show_spag) && !is.null(ub$track_msd) && ncol(ub$track_msd) > 0) {
      tm <- ub$track_msd
      ncols <- ncol(tm); max_lines <- 200L
      set.seed(42L)
      idx_s <- if (ncols > max_lines) sort(sample(ncols, max_lines)) else seq_len(ncols)
      if (!is.null(cl)) {
        cids <- cl$cluster_ids[tids_u]      # cluster per column (NA if not fit)
        for (clk in seq_len(cl$n_clusters)) {
          sx <- numeric(0); sy <- numeric(0); ncl <- 0L
          for (j in idx_s) {
            if (is.na(cids[j]) || cids[j] != clk) next
            yv <- tm[ok, j] / 1e6; okj <- is.finite(yv)
            if (sum(okj) < 2L) next
            sx <- c(sx, tau[okj], NA_real_); sy <- c(sy, yv[okj], NA_real_); ncl <- ncl + 1L
          }
          if (length(sx))
            fig <- fig %>% add_trace(x = sx, y = sy, type = "scattergl", mode = "lines",
              line = list(color = hex_to_rgba(cl$palette[clk], 0.22), width = 1),
              name = sprintf("Cluster %d tracks (n=%d)", clk, ncl),
              hoverinfo = "skip", showlegend = TRUE)
        }
      } else {
        sx <- numeric(0); sy <- numeric(0)
        for (j in idx_s) {
          yv <- tm[ok, j] / 1e6; okj <- is.finite(yv)
          if (sum(okj) < 2L) next
          sx <- c(sx, tau[okj], NA_real_); sy <- c(sy, yv[okj], NA_real_)
        }
        if (length(sx))
          fig <- fig %>% add_trace(x = sx, y = sy, type = "scattergl", mode = "lines",
            line = list(color = "rgba(140,140,140,0.22)", width = 1),
            name = "Per-track TAMSD", hoverinfo = "skip", showlegend = TRUE)
      }
    }

    if (show_band) {
      q25 <- ub$msd_q25[ok] / 1e6; q75 <- ub$msd_q75[ok] / 1e6
      okq <- is.finite(q25) & is.finite(q75)
      if (any(okq)) fig <- fig %>% add_trace(
        x = c(tau[okq], rev(tau[okq])), y = c(q75[okq], rev(q25[okq])),
        type = "scatter", mode = "none", fill = "toself",
        fillcolor = "rgba(70,130,200,0.18)", line = list(color = "transparent"),
        name = "IQR", showlegend = TRUE)
    }
    fig <- fig %>% add_trace(x = tau, y = msd, type = "scatter", mode = "lines+markers",
        marker = list(size = 5, color = "steelblue"), line = list(color = "steelblue", width = 2),
        name = "Unified MSD")
    if (!is.null(fr)) {
      # Model fit curve (from fr$tau_fit / fr$msd_fit; same for every model).
      fx <- if (!is.null(fr$tau_fit)) fr$tau_fit else seq(regime_x, max(tau), length.out = 100)
      fy <- if (!is.null(fr$msd_fit)) pmax(0, fr$msd_fit) / 1e6 else {
        coef_v <- 2L * ub$ndim
        pmax(0, coef_v * fr$D * fx + coef_v * fr$sigma_loc^2) / 1e6
      }
      fit_name <- switch(fr$model %||% "linear",
        confined  = if (isTRUE(fr$degenerate))
                      sprintf("Confined: no plateau in window (≈linear, D=%.3g)", fr$D / 1e3)
                    else sprintf("Confined: D=%.3g µm²/s, L=%.0f nm", fr$D / 1e3, fr$L_conf),
        directed  = sprintf("Directed: D=%.3g µm²/s, v=%.3g µm/s", fr$D / 1e3, fr$v_um_s),
        anomalous = sprintf("Anomalous: D=%.3g, α=%.2f", fr$D / 1e3, fr$alpha),
        sprintf("OLSF D=%.3g µm²/s (%d pts)", fr$D / 1e3, fr$n_points %||% NA))
      fig <- fig %>% add_trace(x = fx, y = fy, type = "scatter", mode = "lines",
        line = list(color = "#E03020", width = 2.5, dash = "dash"), name = fit_name)
    }

    # Per-cluster eMSD: track-weighted mean TAMSD of each cluster's tracks, labelled
    # with the cluster eDIFF + tDIFF (from the cached trk_unified_cluster_diff).
    if (!is.null(cl) && !is.null(ub$track_msd)) {
      cdiff <- tryCatch(trk_unified_cluster_diff(), error = function(e) NULL)
      cids  <- cl$cluster_ids[tids_u]
      for (clk in seq_len(cl$n_clusters)) {
        cols_clk <- which(!is.na(cids) & cids == clk)
        if (length(cols_clk) < 1L) next
        emsd_clk <- rowMeans(ub$track_msd[, cols_clk, drop = FALSE], na.rm = TRUE)[ok] / 1e6
        okc <- is.finite(emsd_clk); if (sum(okc) < 2L) next
        eD <- tD <- NA_real_
        if (!is.null(cdiff)) { r <- cdiff[cdiff$cluster == clk, ]; if (nrow(r)) { eD <- r$eDIFF[1]; tD <- r$tDIFF[1] } }
        nm <- if (is.finite(eD) || is.finite(tD))
                sprintf("Cluster %d eMSD (eDIFF %s, tDIFF %s)", clk,
                        if (is.finite(eD)) sprintf("%.3g", eD) else "—",
                        if (is.finite(tD)) sprintf("%.3g", tD) else "—")
              else sprintf("Cluster %d eMSD", clk)
        fig <- fig %>% add_trace(x = tau[okc], y = emsd_clk[okc], type = "scatter", mode = "lines",
          line = list(color = cl$palette[clk], width = 2.8), name = nm)
      }
    }
    regime <- regime_x
    # N(τ) on the right axis — same warm style as the standard plot.
    show_n2 <- show_ntau && diff(range(ub$n_trk_lag[ok])) > 0
    if (show_n2)
      fig <- fig %>% add_trace(x = tau, y = ub$n_trk_lag[ok], type = "scatter", mode = "lines",
        line = list(color = "rgba(180,100,40,0.70)", width = 1.5, dash = "dot"),
        yaxis = "y2", name = "N(τ) tracks")
    # Lower-regime shading (standard-plot style): shade τ < cutoff as the
    # noise / linker-dominated zone excluded from the fit.
    regime_shapes <- if (is.finite(regime) && regime > 0 && regime < max(tau))
      list(list(type = "rect", xref = "x", yref = "paper",
                x0 = 0, x1 = regime, y0 = 0, y1 = 1,
                fillcolor = "rgba(220,160,80,0.10)",
                line = list(color = "rgba(200,140,60,0.40)", width = 1, dash = "dot")))
      else list()
    lay <- list(xaxis = list(title = if (linear_tau)
                               "τ (ms, linearized = lag × cycle)" else
                               "τ (ms, true mean Δt per bin)", color = fg),
                yaxis = list(title = "MSD (µm²)", color = fg, rangemode = "tozero"),
                paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
                # Legend beside the plot (right), not overlaying it.
                legend = list(font = list(color = fg, size = 10)),
                margin = list(l = 55, r = 10, t = 10, b = 45),
                shapes = regime_shapes, uirevision = "msd_unified")
    if (is.finite(regime) && regime > 0)
      lay$annotations <- list(list(
        x = regime * 0.5, y = 0.97, xref = "x", yref = "paper",
        text = "lower\nregime", showarrow = FALSE,
        font = list(size = 9, color = "rgba(180,120,40,0.70)"),
        xanchor = "center", yanchor = "top"))
    if (show_n2)
      lay$yaxis2 <- list(overlaying = "y", side = "right", title = "N(τ) tracks",
                         showgrid = FALSE, rangemode = "tozero",
                         tickfont = list(color = "rgba(180,100,40,0.85)"),
                         titlefont = list(color = "rgba(180,100,40,0.85)"))
    do.call(plotly::layout, c(list(fig), lay))
  })

  # Cache the unified MSD's σ_loc (from the OLSF intercept) whenever it recomputes,
  # so the D-map and other short-lag tools can use it as the noise floor.
  observeEvent(msd_unified_bundle(), {
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); if (is.null(ub)) return()
    fr <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
    if (!is.null(fr) && is.finite(fr$sigma_loc) && fr$sigma_loc >= 0)
      trk_sigma_loc_rv(fr$sigma_loc)
  }, ignoreInit = TRUE)

  output$trk_sigma_loc_status <- renderUI({
    eff <- trk_sigma_loc_eff()
    txt <- if (is.finite(eff))
             sprintf("Effective σ_loc = %.2f nm (%s). Used as the noise floor everywhere.",
                     eff, trk_sigma_loc_src())
           else paste0("σ_loc not set — estimate it in Analysis → Data Quality & Calibration ",
                       "(no MSD needed), or enter it manually. Noise-floor corrections stay ",
                       "inactive until then.")
    htmltools::div(style = "font-size:0.8em; color:#666; margin-top:-4px;", txt)
  })

  # eDIFF 95% CI by CLUSTER-BOOTSTRAP over tracks: resample track TAMSD columns with
  # replacement, rebuild the (track-weighted) ensemble MSD, refit — the honest CI
  # given within-track correlation. Gated behind a checkbox (B refits are slower);
  # recomputes only when Plot MSD is pressed with the box ticked.
  trk_unified_boot <- eventReactive(input$trk_msd_compute, {
    req(isTRUE(input$trk_uni_bootstrap))
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); req(!is.null(ub))
    tm <- ub$track_msd; req(!is.null(tm), ncol(tm) >= 5L)
    B  <- 300L; nc <- ncol(tm)
    boot_D <- withProgress(message = "Bootstrapping eDIFF 95% CI…", value = 0, {
      vapply(seq_len(B), function(i) {
        if (i %% 30L == 0L) incProgress(30 / B)
        ub2 <- ub
        ub2$msd <- rowMeans(tm[, sample.int(nc, nc, replace = TRUE), drop = FALSE], na.rm = TRUE)
        fr <- tryCatch(msd_unified_fit(ub2), error = function(e) NULL)
        if (is.null(fr) || !is.finite(fr$D)) NA_real_ else fr$D / 1e3   # µm²/s
      }, numeric(1))
    })
    q <- stats::quantile(boot_D, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    list(lo = q[1], hi = q[2], B = sum(is.finite(boot_D)))
  })

  output$trk_msd_unified_estimates <- renderUI({
    req(isTRUE(input$trk_msd_unified))
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); if (is.null(ub)) return(NULL)
    fr <- msd_unified_modelfit(ub)
    if (is.null(fr)) return(htmltools::div(style="color:#888;font-size:0.86em;", "Fit not available."))
    # eDIFF = fr$D is the ensemble-fit D. tDIFF = the time-averaged D (mean of
    # per-track OLSF D, D>0) — computed in the cached trk_unified_tdiff reactive so
    # the heavy per-track fit doesn't rerun on every estimates refresh.
    model <- fr$model %||% "linear"
    linear_tau <- isTRUE(input$trk_uni_linear_tau)
    td    <- trk_unified_tdiff()
    tDIFF <- if (!is.null(td)) td$tDIFF else NA_real_
    n_tD  <- if (!is.null(td)) td$n_tD  else 0L
    tD_ci <- if (!is.null(td) && is.finite(td$ci_lo))
               sprintf(" [95%% CI %.3g–%.3g, SD %.3g]", td$ci_lo, td$ci_hi, td$sd_tD) else ""
    eb    <- tryCatch(trk_unified_boot(), error = function(e) NULL)
    eD_ci <- if (!is.null(eb) && is.finite(eb$lo)) sprintf(" (95%% CI %.3g–%.3g)", eb$lo, eb$hi) else ""
    sig_eff <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)   # noise floor for sanity flags
    ergo  <- if (is.finite(tDIFF) && is.finite(fr$D) && fr$D > 0) tDIFF / (fr$D / 1e3) else NA_real_
    extra <- switch(model,
      confined  = if (isTRUE(fr$degenerate))
                    "  |  ⚠ no plateau reached in window — L unconstrained (MSD still rising; lengthen tracks / raise the upper limit, or this population isn't confined on this timescale)"
                  else sprintf("  |  confinement L = %.0f nm%s  ·  plateau = %.3g µm²  ·  τ_c = %.3g ms",
                               fr$L_conf,
                               if (is.finite(sig_eff) && is.finite(fr$L_conf) && fr$L_conf < 2 * sig_eff)
                                 sprintf(" ⚠ below precision (2·σ_loc=%.0f nm — unresolved)", 2 * sig_eff) else "",
                               fr$plateau_um2, fr$tau_c_ms),
      directed  = {
        vnoise <- if (is.finite(sig_eff) && is.finite(ub$t_cycle) && ub$t_cycle > 0)
                    sqrt(2 * ub$ndim) * sig_eff / ub$t_cycle else NA_real_   # nm/ms = µm/s
        sprintf("  |  velocity v = %.3g µm/s%s", fr$v_um_s,
                if (is.finite(vnoise) && is.finite(fr$v_um_s) && fr$v_um_s < vnoise)
                  sprintf(" ⚠ at/below noise floor (%.3g µm/s — consistent with noise)", vnoise) else "")
      },
      anomalous = sprintf("  |  α = %.2f (%s)", fr$alpha,
                          if (fr$alpha < 0.9) "sub-diffusive / confined"
                          else if (fr$alpha > 1.1) "super-diffusive / directed" else "≈ Brownian"),
      "")
    method_lbl <- switch(model,
      linear = if (identical(fr$method, "OLSF-iterative")) "iterative OLSF" else "anchored Michalet",
      confined = "confined model", directed = "directed model", anomalous = "anomalous model", "model")
    D_lbl <- if (identical(model, "linear")) "eDIFF (ensemble)" else "D (local, ensemble)"
    htmltools::div(style = "font-size:0.87em; color:#444; margin:6px 0;",
      htmltools::tags$b("MSD: "),
      sprintf("%s = %.4g%s  |  tDIFF (per-track mean, n=%d) = %.4g µm²/s%s  |  σ_loc = %.2f nm  |  %s  |  ergodicity = %.2f%s",
              D_lbl, fr$D / 1e3, eD_ci, n_tD, tDIFF, tD_ci, fr$sigma_loc, method_lbl, ergo, extra),
      htmltools::br(),
      htmltools::tags$span(style = "color:#888; font-size:0.9em;",
        sprintf("n = %d tracks%s  |  %s  |  cycle = %.4g ms  |  τ axis = %s",
                ub$n_tracks,
                if (!is.null(ub$n_dropped_negD) && ub$n_dropped_negD > 0)
                  sprintf(" (−%d D≤0)", ub$n_dropped_negD) else "",
                if (isTRUE(ub$track_weighted)) "track-weighted (eMSD)" else "pair-weighted",
                ub$t_cycle,
                if (linear_tau) "linearized lag×cycle" else "real mean Δt (cycle-time-independent)")),
      # Confinement honesty. Bela Vogler's MINFLUX MSD workflow (and Vogler 2025) is built for
      # FREE Brownian diffusion; anomalous/confined MSD is explicitly "under development / implement
      # on your own" in his toolbox. So when fitting Linear, remind the user that in a confined
      # space (e.g. E. coli) the reported D is only an apparent/local D, and point to the CI +
      # Confined/Anomalous models — flagged as beyond-standard-for-MINFLUX.
      if (identical(model, "linear"))
        htmltools::div(style = "color:#8a6d00; font-size:0.82em; margin-top:5px; border-top:1px dashed #e0d0a0; padding-top:4px;",
          htmltools::HTML("&#9432; "),
          htmltools::tags$b("Assumes free Brownian diffusion."),
          " In a confined space (e.g. inside E. coli) the MSD plateaus and this Linear D is only an ",
          "apparent/local D. Check the confinement index (CI, in the Cluster panel) or switch ",
          htmltools::tags$b("Fit model"), " to ", htmltools::tags$b("Confined"),
          " (plateau → domain size L) or ", htmltools::tags$b("Anomalous"),
          " (α<1 = sub-diffusive). Confined/anomalous MSD is beyond the standard MINFLUX workflow ",
          "(Vogler 2025 covers free diffusion) — treat those fits as exploratory."),
      # Per-cluster eDIFF / tDIFF breakdown (when clustering is enabled).
      local({
        cdiff <- tryCatch(trk_unified_cluster_diff(), error = function(e) NULL)
        if (is.null(cdiff) || !nrow(cdiff)) return(NULL)
        rows <- lapply(seq_len(nrow(cdiff)), function(i) {
          r <- cdiff[i, ]
          htmltools::tags$div(style = sprintf("color:%s; font-size:0.85em; margin-left:8px;", r$color),
            sprintf("● Cluster %d: eDIFF = %s, tDIFF = %s µm²/s  (n = %d)", r$cluster,
                    if (is.finite(r$eDIFF)) sprintf("%.3g", r$eDIFF) else "—",
                    if (is.finite(r$tDIFF)) sprintf("%.3g", r$tDIFF) else "—", r$n))
        })
        htmltools::tagList(htmltools::br(),
          htmltools::tags$b(style = "font-size:0.85em;", "Per-cluster diffusion:"), rows)
      }))
  })

  # Head-to-head comparison table: Standard vs Advanced vs Unified
  output$trk_msd_compare_table <- renderUI({
    req(isTRUE(input$trk_msd_unified))
    getD <- function(expr) tryCatch(expr, error = function(e) NA_real_)
    # Unified
    ub <- getD(msd_unified_bundle()); fr_u <- if (!is.null(ub)) msd_unified_fit(ub) else NULL
    D_u <- if (!is.null(fr_u)) fr_u$D / 1e3 else NA_real_
    # Advanced (cycle)
    cb <- getD(msd_cycle_bundle())
    D_a <- if (!is.null(cb)) {
      ml <- max(1L, as.integer(isolate(input$trk_msd_cycle_min_lag) %||% 4L))
      okf <- is.finite(cb$msd) & cb$lag >= ml
      fr_a <- if (sum(okf) >= 3L) tryCatch(msd_fit_olsf(cb$tau_ms[okf], cb$msd[okf],
                cb$npairs[okf], sum(okf), ndim = cb$ndim, auto_pmin = TRUE, R = 0,
                Npts = 1e4, dt_cycle = cb$t_cycle),
                error = function(e) NULL) else NULL
      if (!is.null(fr_a)) fr_a$D / 1e3 else NA_real_
    } else NA_real_
    row <- function(lbl, d) htmltools::tags$tr(
      htmltools::tags$td(lbl, style = "padding:2px 10px;"),
      htmltools::tags$td(if (is.finite(d)) sprintf("%.4g µm²/s", d) else "—",
                         style = "padding:2px 10px; font-weight:bold;"))
    htmltools::div(style = "margin-top:8px;",
      htmltools::tags$table(style = "font-size:0.86em; border-collapse:collapse;",
        htmltools::tags$thead(htmltools::tags$tr(
          htmltools::tags$th("Estimator", style="padding:2px 10px;border-bottom:1px solid #ccc;"),
          htmltools::tags$th("D", style="padding:2px 10px;border-bottom:1px solid #ccc;"))),
        htmltools::tags$tbody(
          row("Unified (cycle-bin + mean-Δt + OLSF/p_min)", D_u),
          row("Advanced (cycle MSD)", D_a))),
      htmltools::tags$p(style = "font-size:0.8em; color:#888; margin-top:4px;",
        "Standard MSD D is shown in its own panel above. Compare all three; ",
        "the Unified should track the Advanced while being cycle-time-independent."))
  })

  # ---- Integrated / Dark cycle classification --------------------------------
  # Approximates Vogler's TSA cycle decomposition: each consecutive interval is
  # mapped to an integer number of acquisition cycles, n = round(Δt / t_cycle).
  # n = 1 → Integrated (localized every cycle).  n ≥ 2 → (n−1) Dark cycles passed
  # before the next localization.  Uses only the Δt distribution — no hardware
  # simulation — so the cycle time is the modal Δt rather than a fitted overhead.
  cycle_class_bundle <- eventReactive(input$trk_tsa_run, {
    req(isTRUE(input$trk_show_tsa))
    withProgress(message = "Time Signature Analysis…", value = 0.5, {
    t_c <- trk_cycle_time()   # single source of truth
    dts <- dt_vectors()
    vec <- unlist(dts, use.names = FALSE)
    vec <- vec[is.finite(vec) & vec > 0]
    shiny::validate(shiny::need(length(vec) > 0, "No Δt data for cycle classification."))
    shiny::validate(shiny::need(is.finite(t_c) && t_c > 0, "Could not determine cycle time."))
    hw_used <- { hw <- trk_hw_cycle_time(); is.finite(hw) && hw > 0 }

    n_cyc      <- pmax(1L, as.integer(round(vec / t_c)))
    n_total    <- length(n_cyc)
    n_integ    <- sum(n_cyc == 1L)
    n_dark_int <- sum(n_cyc >= 2L)           # intervals containing ≥1 dark cycle
    dark_cycles_total <- sum(n_cyc - 1L)     # total dark cycles across all intervals
    frac_integ <- n_integ / n_total

    max_disp <- max(5L, as.integer(input$trk_cycle_dark_max %||% 20L))
    tab      <- tabulate(pmin(n_cyc, max_disp), nbins = max_disp)
    list(t_cycle = t_c, n_cyc = n_cyc, counts = tab, max_disp = max_disp,
         n_total = n_total, n_integ = n_integ, n_dark_int = n_dark_int,
         dark_cycles_total = dark_cycles_total, frac_integ = frac_integ,
         hw_used = hw_used)
    })  # withProgress
  })

  output$trk_cycle_class_plot <- renderPlotly({
    cc <- cycle_class_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"

    x_lab  <- seq_len(cc$max_disp)
    frac   <- cc$counts / cc$n_total
    # Integrated (n=1) green, Dark (n>=2) amber→red gradient
    cols <- ifelse(x_lab == 1L, "#2CA02C",
              grDevices::colorRampPalette(c("#E0A030", "#C03030"))(cc$max_disp)[x_lab])
    x_txt <- as.character(x_lab); x_txt[cc$max_disp] <- paste0("≥", cc$max_disp)

    plot_ly(x = x_lab, y = frac, type = "bar",
            marker = list(color = cols),
            text = sprintf("%d cycle%s: %.1f%%", x_lab, ifelse(x_lab==1,"","s"), frac*100),
            hovertemplate = "%{text}<extra></extra>") %>%
      layout(
        xaxis = list(title = "Cycles per interval (1 = Integrated, ≥2 = Dark)",
                     tickmode = "array", tickvals = x_lab, ticktext = x_txt,
                     color = fg_col),
        yaxis = list(title = "Fraction of intervals", color = fg_col),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font = list(color = fg_col), showlegend = FALSE,
        bargap = 0.1, uirevision = "cycle_class"
      )
  })

  output$trk_cycle_class_summary <- renderUI({
    cc <- cycle_class_bundle()
    htmltools::div(
      style = "font-size:0.87em; color:#444; margin:6px 0;",
      htmltools::tags$b("Cycle classification "),
      sprintf("(cycle time = %.4g ms, %s):", cc$t_cycle,
              if (isTRUE(cc$hw_used)) "snapped to hardware model" else "from modal Δt"),
      htmltools::br(),
      sprintf("Integrated (1 cycle): %.1f%% of %s intervals  |  with Dark cycles (≥2): %.1f%%",
              cc$frac_integ * 100, format(cc$n_total, big.mark = ","),
              (1 - cc$frac_integ) * 100),
      htmltools::br(),
      sprintf("Total Dark cycles: %s  |  mean cycles per interval: %.2f",
              format(cc$dark_cycles_total, big.mark = ","),
              mean(cc$n_cyc)),
      htmltools::br(),
      htmltools::tags$span(
        style = "color:#888; font-size:0.9em;",
        if (cc$frac_integ > 0.9)
          "Excellent tracking quality — almost every cycle yielded a localization."
        else if (cc$frac_integ > 0.7)
          "Good tracking — most cycles localized, modest dark fraction."
        else
          "Frequent dropouts — consider that Dark cycles add temporal gaps; the gap-split filter (Display tab) may help downstream analysis."
      )
    )
  })

  # Cycle stats on an arbitrary df (Integrated/Dark) at a given cycle time.
  compute_cycle_stats <- function(df, t_c) {
    if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df)) ||
        !is.finite(t_c) || t_c <= 0) return(NULL)
    tf  <- match(df$tid, unique(df$tid))           # integer codes → radix sort (fast)
    ord <- order(tf, df$tim, method = "radix")
    tf  <- tf[ord]; tim <- df$tim[ord]; n <- length(tf)
    same <- c(FALSE, tf[-1] == tf[-n])             # consecutive pair within same track
    dts  <- c(NA_real_, diff(tim))[same]           # within-track gaps only
    dts <- dts[is.finite(dts) & dts > 0]
    if (!length(dts)) return(NULL)
    n_cyc <- pmax(1L, as.integer(round(dts / t_c)))
    list(n_tracks   = length(unique(tf)),
         n_total    = length(n_cyc),
         frac_integ = mean(n_cyc == 1L),
         dark_total = sum(n_cyc - 1L),
         mean_cyc   = mean(n_cyc),
         max_cyc    = max(n_cyc))
  }

  # Raw vs analyzed cycle-stats: the NET effect of all track cleaning (gap-split +
  # end-trim + quality gate + metadata + length/duration) on the Integrated/Dark
  # balance. "Raw" = the uploaded data; "Analyzed" = what the analyses actually use.
  # They differ whenever any cleaning is active; identical means nothing is cleaning.
  output$trk_cycle_class_beforeafter <- renderUI({
    req(isTRUE(input$trk_show_tsa), tsa_ran())
    t_c    <- trk_cycle_time()
    before <- compute_cycle_stats(tryCatch(trk_raw_original(),     error = function(e) NULL), t_c)
    after  <- compute_cycle_stats(tryCatch(trk_spatial_filtered(), error = function(e) NULL), t_c)
    if (is.null(before)) return(NULL)
    any_clean <- isTRUE(input$trk_split_gaps_enable) || isTRUE(input$trk_endtrim_enable) ||
                 isTRUE(input$trk_quality_enable)

    cell <- function(x) htmltools::tags$td(x, style = "padding:2px 10px; text-align:right;")
    row  <- function(lbl, b, a) htmltools::tags$tr(
      htmltools::tags$td(lbl, style = "padding:2px 10px;"), cell(b), cell(a))
    fmt  <- function(v, pct = FALSE) if (is.null(v)) "—" else
      if (pct) sprintf("%.1f%%", v * 100) else format(round(v, 2), big.mark = ",")

    htmltools::div(style = "font-size:0.86em; color:#444; margin:8px 0;",
      htmltools::tags$b("Cycle stats — raw vs analyzed (after cleaning)"),
      if (!any_clean) htmltools::tags$span(style = "color:#888;",
        " (no cleaning active — enable Split / Trim / Quality gate in the Display tab to see the effect)"),
      htmltools::tags$table(style = "border-collapse:collapse; margin-top:4px;",
        htmltools::tags$thead(htmltools::tags$tr(
          htmltools::tags$th("", style = "padding:2px 10px; border-bottom:1px solid #ccc;"),
          htmltools::tags$th("Raw", style = "padding:2px 10px; border-bottom:1px solid #ccc;"),
          htmltools::tags$th("Analyzed", style = "padding:2px 10px; border-bottom:1px solid #ccc;"))),
        htmltools::tags$tbody(
          row("Tracks / segments", fmt(before$n_tracks), fmt(after$n_tracks)),
          row("Integrated (1 cycle)", fmt(before$frac_integ, TRUE), fmt(after$frac_integ, TRUE)),
          row("With Dark cycles (≥2)", fmt(1 - before$frac_integ, TRUE), fmt(1 - after$frac_integ, TRUE)),
          row("Total Dark cycles", fmt(before$dark_total), fmt(after$dark_total)),
          row("Mean cycles / interval", fmt(before$mean_cyc), fmt(after$mean_cyc)),
          row("Max cycles in one gap", fmt(before$max_cyc), fmt(after$max_cyc)))),
      htmltools::tags$p(style = "color:#888; font-size:0.95em; margin-top:4px;",
        "Cleaning (split at dark gaps, trim degraded ends, drop dark-heavy tracks) should raise ",
        "the Integrated fraction, cut total Dark cycles, and shrink the max gap — the analyzed ",
        "intervals are clean, near-consecutive localizations."))
  })

  # Per-track cycle progression (plot_cycle_progression_for_ID rewrite to R)
  output$trk_cycle_prog_tid_ui <- renderUI({
    req(isTRUE(input$trk_show_tsa))
    df <- tryCatch(trk(), error = function(e) NULL)
    req(!is.null(df), nrow(df) > 0)
    counts <- table(df$tid)
    tids   <- names(counts[counts >= 5L])
    req(length(tids) > 0)
    # Offer longest tracks first
    tids <- tids[order(as.numeric(counts[tids]), decreasing = TRUE)]
    selectizeInput("trk_cycle_prog_tid", "Track ID:",
                   choices = tids, selected = tids[1],
                   options = list(maxOptions = 5000))
  })

  output$trk_cycle_prog_plot <- renderPlotly({
    req(isTRUE(input$trk_show_tsa), tsa_ran())
    req(nzchar(input$trk_cycle_prog_tid))
    df  <- trk(); req(nrow(df) > 0)
    t_c <- trk_cycle_time()
    if (!is.finite(t_c) || t_c <= 0) t_c <- 1
    tb <- df[df$tid == input$trk_cycle_prog_tid, , drop = FALSE]
    tb <- tb[order(tb$tim), ]
    shiny::validate(shiny::need(nrow(tb) >= 2, "Track too short."))

    # Cluster (atomic) cycle + its SW offset for the cycle-count mapping
    cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    sw <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0

    dt_v  <- diff(tb$tim)
    n_cyc <- pmax(1L, as.integer(round((dt_v - sw) / t_c)))  # total cycles to next loc
    dark  <- n_cyc - 1L                                       # dark cycles (0 = Integrated)
    step  <- seq_along(n_cyc)
    is_int <- dark == 0L

    dm <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dm) "#000000" else "#FFFFFF"
    fg_col <- if (dm) "#FFFFFF" else "#333333"
    hov <- function(i) sprintf("step %d: %d cycle%s = %s (%.3g ms)",
              step[i], n_cyc[i], ifelse(n_cyc[i]==1,"","s"),
              ifelse(is_int[i], "Integrated", sprintf("%d Dark", dark[i])), dt_v[i])

    plot_ly() %>%
      add_trace(x = step, y = dark, type = "scatter", mode = "lines",
                line = list(color = "rgba(120,120,120,0.4)", width = 1),
                showlegend = FALSE, hoverinfo = "skip") %>%
      add_trace(x = step[is_int], y = dark[is_int], type = "scatter", mode = "markers",
                marker = list(color = "#2070D0", size = 6),
                name = "Integrated (0 dark)",
                text = vapply(which(is_int), hov, character(1)),
                hovertemplate = "%{text}<extra></extra>") %>%
      add_trace(x = step[!is_int], y = dark[!is_int], type = "scatter", mode = "markers",
                marker = list(color = "#C03030", size = 6),
                name = "Dark cycles present",
                text = vapply(which(!is_int), hov, character(1)),
                hovertemplate = "%{text}<extra></extra>") %>%
      layout(
        xaxis = list(title = "Localization index along track", color = fg_col),
        yaxis = list(title = "Dark cycles before this localization", color = fg_col,
                     rangemode = "tozero"),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font = list(color = fg_col),
        legend = list(orientation = "h", x = 0, y = 1.12, font = list(color = fg_col)),
        uirevision = "cycle_prog"
      )
  })

  # Plain CSV: tidy columns for Excel / Prism
  output$trk_dl_msd_csv <- downloadHandler(
    filename = function() {
      ex <- msd_export(); mode <- ex$params$mode
      paste0("msd_", mode %||% "time", ".csv")
    },
    content = function(file) {
      ex <- msd_export()
      utils::write.csv(ex$data, file, row.names = FALSE)
    }
  )
  
  # ZIP: includes msd_data.csv + params.txt (easy to read in Excel/Notepad/Prism)
  output$trk_dl_msd_zip <- downloadHandler(
    filename = function() {
      ex <- msd_export(); mode <- ex$params$mode
      paste0("msd_", mode %||% "time", "_export.zip")
    },
    content = function(file) {
      ex <- msd_export()
      td <- tempdir()
      f_csv <- file.path(td, "msd_data.csv")
      f_txt <- file.path(td, "params.txt")
      
      # data
      utils::write.csv(ex$data, f_csv, row.names = FALSE)
      
      # human-readable params
      lines <- c(
        "MSD export parameters:",
        sprintf("mode: %s", ex$params$mode),
        sprintf("time_unit: %s", ex$params$time_unit),
        sprintf("space_unit: %s", ex$params$space_unit),
        sprintf("tau_value: %s", ex$params$tau_value),
        sprintf("msd_min_pts: %s", ex$params$msd_min_pts),
        sprintf("msd_weighting: %s", ex$params$msd_weighting),
        sprintf("msd_fit_method: %s", ex$params$msd_fit_method),
        sprintf("msd_tau_max: %s", ex$params$msd_tau_max),
        sprintf("msd_tau_bin_used_ms: %s", ex$params$msd_tau_bin_used),
        sprintf("msd_tau_bin_auto: %s", ex$params$msd_tau_bin_auto),
        sprintf("msd_fit_min_lag_cycles: %s", ex$params$msd_fit_min_lag_cycles),
        sprintf("msd_fit_min_lag_ms: %s", ex$params$msd_fit_min_lag_ms),
        sprintf("msd_ndim: %s", ex$params$msd_ndim),
        sprintf("msd_formula: %s", ex$params$msd_formula),
        sprintf("dt0_ms (median dt): %s", ex$params$dt0_ms)
      )
      writeLines(lines, f_txt)
      
      # zip it
      zip::zipr(zipfile = file, files = c(f_csv, f_txt), root = td)
    }
  )
  
  # ---- Jump Distance Distribution ------------------------------------------
  jdd_bundle <- eventReactive(input$trk_jdd_run, {
    req(isTRUE(input$trk_show_jdd))
    df <- trk(); req(nrow(df) > 0)
    min_pts   <- max(2L, as.integer(input$trk_jdd_min_pts))
    norm_dt   <- isTRUE(input$trk_jdd_norm_dt)
    n_comp    <- max(1L, min(4L, as.integer(input$trk_jdd_n_comp)))
    min_dt    <- suppressWarnings(as.numeric(input$trk_jdd_min_dt))
    max_dt    <- suppressWarnings(as.numeric(input$trk_jdd_max_dt))
    sigma_loc <- suppressWarnings(as.numeric(input$trk_jdd_sigma_loc))
    step_mode <- input$trk_jdd_step_mode %||% "consec"
    bin_width <- suppressWarnings(as.numeric(input$trk_jdd_bin_width))
    subsamp_n <- suppressWarnings(as.integer(input$trk_jdd_subsamp_n))
    if (!is.finite(min_dt)    || min_dt    <= 0) min_dt    <- 0
    if (!is.finite(max_dt)    || max_dt    <= 0) max_dt    <- Inf
    if (!is.finite(sigma_loc) || sigma_loc <= 0) {   # 0 = auto from app-wide σ_loc
      s_app <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
      sigma_loc <- if (is.finite(s_app) && s_app > 0) s_app else 0
    }
    if (!is.finite(bin_width) || bin_width <= 0) bin_width <- 1
    if (!is.finite(subsamp_n) || subsamp_n <  2) subsamp_n <- 4L

    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
    df <- df[order(df$tid, df$tim), ]
    sp <- split(df, df$tid)

    use_z_jdd <- "z" %in% names(df) && sum(is.finite(df$z)) > 0
    jumps  <- numeric(0)
    dt_all <- numeric(0)
    withProgress(message = "Computing JDD…", value = 0, {
    n_sp_jdd <- length(sp)
    for (ti_jdd in seq_along(sp)) {
      if (ti_jdd %% 50L == 0L) incProgress(50 / n_sp_jdd, detail = "collecting steps")
      tb <- sp[[ti_jdd]]
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
          if (use_z_jdd && "z" %in% names(tb)) row$z <- median(tb$z[idx])
          bin_rows[[bi_i]] <- row
        }
        tb_use <- do.call(rbind, bin_rows)
      } else if (identical(step_mode, "subsamp")) {
        # Take every N-th localization
        idx_keep <- seq(1L, nrow(tb), by = subsamp_n)
        if (length(idx_keep) < 2) next
        tb_use <- tb[idx_keep, , drop = FALSE]
      } else {
        tb_use <- tb   # consecutive steps — original behaviour
      }
      if (nrow(tb_use) < 2) next

      dx <- diff(tb_use$x); dy <- diff(tb_use$y)
      r  <- if (use_z_jdd && "z" %in% names(tb_use))
               sqrt(dx^2 + dy^2 + diff(tb_use$z)^2)
            else sqrt(dx^2 + dy^2)
      dt <- diff(tb_use$tim)
      keep <- is.finite(dt) & dt > min_dt & dt <= max_dt & is.finite(r)
      if (!any(keep)) next
      if (norm_dt) {
        r    <- r / sqrt(dt)
        keep <- keep & is.finite(r)
      }
      if (!any(keep)) next
      dt_all <- c(dt_all, dt[keep])
      jumps  <- c(jumps,  r[keep])
    }
    shiny::validate(shiny::need(length(jumps) >= 10,
      "Too few jumps — lower 'Min points per TID' or relax filters."))

    # Histogram (truncate upper tail at 99.5th pct to avoid outlier stretching)
    q_hi <- quantile(jumps, 0.995, na.rm = TRUE)
    br   <- c(pretty(c(0, q_hi), n = 60), Inf)
    h    <- hist(jumps, breaks = br, plot = FALSE)
    mids <- h$mids[is.finite(h$mids)]
    dens <- h$density[seq_along(mids)]

    # --- Rayleigh mixture MLE ------------------------------------------------
    # p(r) = sum_k w_k * (r/s_k^2) * exp(-r^2 / (2*s_k^2))
    # Parameters: log(sigma_k), and for k>1: log(w_k / w_K) (softmax)
    rayleigh_mix_nll <- function(par) {
      sigs <- exp(par[seq_len(n_comp)])
      if (n_comp > 1) {
        raw  <- exp(par[n_comp + seq_len(n_comp - 1)])
        ws   <- c(raw, 1) / (1 + sum(raw))
      } else {
        ws <- 1
      }
      if (any(!is.finite(sigs)) || any(!is.finite(ws))) return(1e15)
      lik <- Reduce(`+`, lapply(seq_len(n_comp), function(k)
        ws[k] * (jumps / sigs[k]^2) * exp(-jumps^2 / (2 * sigs[k]^2))
      ))
      lik[lik <= 0] <- .Machine$double.eps
      -sum(log(lik))
    }

    sig0 <- seq(quantile(jumps, 0.05), quantile(jumps, 0.75), length.out = n_comp)
    sig0 <- pmax(sig0, 1e-3)
    par0 <- c(log(sig0), rep(0, max(0, n_comp - 1)))

    incProgress(0.4, detail = "fitting Rayleigh mixture")
    fit <- tryCatch(
      optim(par0, rayleigh_mix_nll, method = "Nelder-Mead",
            control = list(maxit = 8000, reltol = 1e-9)),
      error = function(e) NULL
    )

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

    # Median ÃÂt across all valid steps used in this fit (ms)
    dt0 <- median(dt_all, na.rm = TRUE)
    if (!is.finite(dt0) || dt0 <= 0) dt0 <- NA_real_

    # D from sigma (when norm_dt = TRUE, ÃÂ is in nm/Ã¢ÂÂms so ÃÂ²/2 is nm²/ms):
    #   2D (Rayleigh): E[r²] = 2ÃÂ²; 2D MSD = 4D·ÃÂt  Ã¢ÂÂ  D = ÃÂ²/2
    #   3D (Maxwell):  E[r²] = 3ÃÂ²; 3D MSD = 6D·ÃÂt  Ã¢ÂÂ  D = ÃÂ²/3
    # When z is included the fit is still Rayleigh-shaped but ÃÂ reflects 3D spread.
    divisor <- if (use_z_jdd) 3 else 2

    # ÃÂ_loc correction (only meaningful when norm_dt = TRUE):
    # The measured ÃÂ_fit² = 2D + 2ÃÂ_loc²/ÃÂt_typical (per-coordinate loc noise adds
    # 2ÃÂ_loc² to the raw step variance; after Ã¢ÂÂÃÂt normalisation this becomes 2ÃÂ_loc²/ÃÂt).
    # We subtract the bias term before computing D, clamping to ≥ 0.
    loc_bias <- if (norm_dt && sigma_loc > 0 && is.finite(dt0))
                  2 * sigma_loc^2 / dt0
                else 0
    D_vals <- pmax(0, sigs^2 - loc_bias) / divisor

    # Localization-noise floor as a JUMP DISTANCE in the plot's units. A pure-noise
    # jump has RMS √(2·d)·σ_loc (d = 2/3); after √Δt normalisation, divide by √dt0.
    # Peaks/steps left of it are noise-dominated, not real motion.
    d_dim_jdd <- if (use_z_jdd) 3 else 2
    r_noise <- if (sigma_loc > 0) {
      base <- sqrt(2 * d_dim_jdd) * sigma_loc
      if (norm_dt && is.finite(dt0) && dt0 > 0) base / sqrt(dt0) else base
    } else NA_real_

    list(mids = mids, dens = dens, sigs = sigs, ws = ws, D_vals = D_vals,
         n_comp = n_comp, norm_dt = norm_dt, n_jumps = length(jumps),
         use_z = use_z_jdd, dt0 = dt0, sigma_loc = sigma_loc, r_noise = r_noise,
         step_mode = step_mode,
         min_dt_used = if (min_dt > 0) min_dt else NA_real_,
         max_dt_used = if (is.finite(max_dt)) max_dt else NA_real_)
    })  # withProgress
  })

  output$trk_jdd_plot <- renderPlotly({
    bd <- jdd_bundle()
    shiny::validate(shiny::need(!is.null(bd), "No jump data."))
    dark   <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"

    r_seq <- seq(0, max(bd$mids) * 1.15, length.out = 500)

    rayleigh_comp <- function(r, sig, w)
      w * (r / sig^2) * exp(-r^2 / (2 * sig^2))

    total_dens <- Reduce(`+`, lapply(seq_len(bd$n_comp), function(k)
      rayleigh_comp(r_seq, bd$sigs[k], bd$ws[k])
    ))

    comp_cols <- c("#E07040", "#4080D0", "#50B060", "#B040B0")

    fig <- plot_ly() %>%
      add_bars(
        x = bd$mids, y = bd$dens,
        name = "Jumps",
        marker = list(color = "rgba(120,160,220,0.6)",
                      line  = list(color = "rgba(80,120,180,0.8)", width = 0.5))
      ) %>%
      add_lines(
        x = r_seq, y = total_dens,
        name = "Fit (total)",
        line = list(color = "#E03020", width = 2.5)
      )

    if (bd$n_comp > 1) {
      for (k in seq_len(bd$n_comp)) {
        fig <- fig %>% add_lines(
          x = r_seq,
          y = rayleigh_comp(r_seq, bd$sigs[k], bd$ws[k]),
          name = sprintf("Comp %d  D=%.2f", k, bd$D_vals[k]),
          line = list(color = comp_cols[k], width = 1.5, dash = "dash")
        )
      }
    }

    if (is.finite(bd$r_noise %||% NA_real_) && bd$r_noise > 0 && bd$r_noise < max(r_seq)) {
      ymax <- max(c(bd$dens, total_dens), na.rm = TRUE)
      fig <- fig %>% add_lines(
        x = c(bd$r_noise, bd$r_noise), y = c(0, ymax),
        name = "loc-noise floor",
        line = list(color = "#888888", width = 1.5, dash = "dot"),
        hovertemplate = sprintf("loc-noise floor \u2248 %.1f (left = noise-dominated)<extra></extra>", bd$r_noise)
      )
    }

    x_lab <- if (bd$norm_dt) "Step size / \u221a\u0394t  (nm / \u221ams)"
             else            "Step size  (nm)"

    fig %>% layout(
      xaxis = list(title = x_lab,   color = fg_col, zeroline = FALSE),
      yaxis = list(title = "Density", color = fg_col, zeroline = FALSE),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font    = list(color = fg_col),
      legend  = list(font = list(color = fg_col)),
      bargap  = 0.05,
      uirevision = "jdd"
    )
  })

  output$trk_jdd_estimates <- renderUI({
    bd <- jdd_bundle()
    if (is.null(bd)) return(NULL)

    show_D    <- isTRUE(bd$norm_dt)
    d_formula <- if (isTRUE(bd$use_z)) "D = (\u03C3\u00B2 \u2212 2\u03C3\u1D38\u1D52\u1D9C\u00B2/\u27E8\u0394t\u27E9)/3" else
                                        "D = (\u03C3\u00B2 \u2212 2\u03C3\u1D38\u1D52\u1D9C\u00B2/\u27E8\u0394t\u27E9)/2"
    dt0_str   <- if (is.finite(bd$dt0)) sprintf("%.2g ms", bd$dt0) else "n/a"
    max_dt_str <- if (is.finite(bd$max_dt_used)) sprintf("%.4g ms", bd$max_dt_used) else "none"

    sigma_note <- if (show_D && bd$sigma_loc > 0 && is.finite(bd$dt0))
      tags$p(sprintf("\u03C3\u1D38\u1D52\u1D9C = %.1f nm  |  \u27E8\u0394t\u27E9 = %s  |  bias subtracted: %.2g nm\u00B2/ms",
                     bd$sigma_loc, dt0_str, 2 * bd$sigma_loc^2 / bd$dt0),
             style = "font-size:0.82em; color:#888; margin:2px 0;")
    else NULL

    col_headers <- if (show_D) c("", "Weight", "Width (\u03C3)", "Diffusivity (D, nm\u00B2/ms)")
                   else        c("", "Weight", "Width (\u03C3)")

    rows <- lapply(seq_len(bd$n_comp), function(k) {
      cells <- list(
        tags$td(sprintf("Component %d", k),          style = "padding:2px 8px;"),
        tags$td(sprintf("%.1f%%", bd$ws[k] * 100),   style = "padding:2px 8px;"),
        tags$td(sprintf("\u03C3 = %.1f nm/\u221Ams", bd$sigs[k]), style = "padding:2px 8px;")
      )
      if (show_D) {
        d_val <- bd$D_vals[k]
        d_str <- if (d_val <= 0 || d_val < 1e-6)
          tags$span("D ≈ 0", tags$span(style="color:#c05020; font-size:0.88em;",
            " (within loc. noise — step sizes explained by σₗₒₓ alone)"))
        else
          sprintf("D = %.3g", d_val)
        cells <- c(cells, list(tags$td(d_str, style = "padding:2px 8px;")))
      }
      do.call(tags$tr, cells)
    })

    no_d_warn <- if (!show_D)
      tags$p("\u26A0 D values are not shown: without \u221A\u0394t normalisation \u03C3 has units nm, not nm/\u221Ams, so \u03C3\u00B2/2 is not a diffusion coefficient.",
             style = "font-size:0.82em; color:#c05020; margin:4px 0;")
    else NULL

    min_dt_str <- if (is.finite(bd$min_dt_used) && bd$min_dt_used > 0)
                     sprintf("%.4g ms", bd$min_dt_used) else "none"
    mode_str <- switch(bd$step_mode %||% "consec",
      binmed  = "bin-then-median",
      subsamp = "subsampled",
      "consecutive steps"
    )
    tags$div(
      tags$p(sprintf("n = %d jumps  |  %s  |  %s  |  \u0394t window: [%s, %s]",
                     bd$n_jumps, d_formula, mode_str, min_dt_str, max_dt_str),
             style = "font-size:0.85em; color:#888; margin:4px 0;"),
      sigma_note,
      no_d_warn,
      tags$table(
        style = "font-size:0.88em; border-collapse:collapse;",
        tags$thead(tags$tr(
          lapply(col_headers, function(h)
            tags$th(h, style = "padding:2px 8px; border-bottom:1px solid #ccc;"))
        )),
        tags$tbody(rows)
      )
    )
  })

  output$trk_dl_jdd_csv <- downloadHandler(
    filename = function() paste0("jdd_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      bd <- jdd_bundle()
      shiny::validate(shiny::need(!is.null(bd), "No JDD computed."))
      rayleigh_comp <- function(r, sig, w)
        w * (r / sig^2) * exp(-r^2 / (2 * sig^2))
      out <- data.frame(mid_nm = bd$mids, density = bd$dens)
      out$fit_total <- Reduce(`+`, lapply(seq_len(bd$n_comp), function(k)
        rayleigh_comp(bd$mids, bd$sigs[k], bd$ws[k])
      ))
      for (k in seq_len(bd$n_comp))
        out[[paste0("comp_", k)]] <- rayleigh_comp(bd$mids, bd$sigs[k], bd$ws[k])
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ---- Confinement analysis ------------------------------------------------
  conf_bundle <- eventReactive(input$trk_conf_run, {
    req(isTRUE(input$trk_show_conf))
    df <- trk(); req(nrow(df) > 0)
    min_pts <- max(4L, as.integer(input$trk_conf_min_pts))

    use_z_conf <- "z" %in% names(df) && sum(is.finite(df$z)) > 0
    df <- df[is.finite(df$x) & is.finite(df$y) & (if (use_z_conf) is.finite(df$z) else TRUE), ]
    sp <- split(df, df$tid)

    # Optional σ_loc correction: localization noise inflates Rg. Measured per-axis
    # variance = true variance + σ_loc², so Rg²_meas = Rg²_true + d·σ_loc²·(1−1/n)
    # (d = 2 or 3; the (1−1/n) accounts for estimating the centroid from the same
    # noisy points). Subtract it, clamped ≥ 0. Rmax (an extreme value) is not
    # noise-corrected and left as-is. No-op when σ_loc is unset.
    d_dim     <- if (use_z_conf) 3 else 2
    sig_loc   <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    sub_noise <- isTRUE(input$trk_conf_sub_noise) && is.finite(sig_loc) && sig_loc > 0

    rows <- withProgress(message = "Computing confinement…", value = 0,
      lapply(seq_along(sp), function(i) {
      if (i %% 50L == 0L) incProgress(50 / length(sp))
      nm <- names(sp)[i]
      tb <- sp[[nm]]
      if (nrow(tb) < min_pts) return(NULL)
      cx <- mean(tb$x); cy <- mean(tb$y)
      cz <- if (use_z_conf) mean(tb$z) else 0
      d2 <- (tb$x - cx)^2 + (tb$y - cy)^2 + if (use_z_conf) (tb$z - cz)^2 else 0
      rg2 <- mean(d2)
      if (sub_noise) rg2 <- rg2 - d_dim * sig_loc^2 * (1 - 1 / nrow(tb))
      Rg   <- sqrt(max(0, rg2))
      Rmax <- max(sqrt(d2))
      if (!is.finite(Rmax) || Rmax == 0) return(NULL)
      ratio <- Rg / Rmax
      data.frame(tid = nm, Rg_nm = Rg, Rmax_nm = Rmax,
                 ratio = ratio, n_pts = nrow(tb),
                 stringsAsFactors = FALSE)
    }))
    rows <- rows[!sapply(rows, is.null)]
    shiny::validate(shiny::need(length(rows) > 0,
      "Too few tracks — lower 'Min points per TID'."))
    res <- do.call(rbind, rows)
    attr(res, "sub_noise") <- sub_noise
    attr(res, "sigma_loc") <- if (sub_noise) sig_loc else NA_real_
    res
  })

  output$trk_conf_plot <- renderPlotly({
    res  <- conf_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"

    # colour by ratio (confined = red, free = blue)
    plot_ly(res,
      x = ~Rmax_nm, y = ~Rg_nm,
      color = ~ratio, colors = c("#C03030","#808080","#3060C0"),
      text  = ~sprintf("TID: %s<br>R\u2098: %.1f nm<br>R\u2090: %.1f nm<br>ratio: %.3f<br>n: %d",
                       tid, Rg_nm, Rmax_nm, ratio, n_pts),
      hovertemplate = "%{text}<extra></extra>",
      type = "scatter", mode = "markers",
      marker = list(size = 6, opacity = 0.75)
    ) %>% layout(
      xaxis  = list(title = "R\u2098\u2090\u02E3 (nm)", color = fg_col),
      yaxis  = list(title = "R\u2092 (nm)",              color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col),
      annotations = list(list(
        text = "Confined \u2192 R\u2092/R\u2098\u2090\u02E3 \u2248 0.25  |  Free \u2192 \u2248 0.5",
        xref = "paper", yref = "paper", x = 0.5, y = 1.05,
        showarrow = FALSE, font = list(size = 11, color = fg_col)
      )),
      uirevision = "conf"
    )
  })

  output$trk_conf_summary <- renderUI({
    res <- conf_bundle()
    rc <- .mean_ci(res$ratio)                       # track-level median/IQR + mean 95% CI
    med_ratio <- rc$median
    pct_conf  <- mean(res$ratio < 0.35, na.rm = TRUE) * 100
    pct_free  <- mean(res$ratio > 0.45, na.rm = TRUE) * 100
    ci_note <- if (is.finite(rc$q1))
      tags$span(style = "color:#555;",
        sprintf("  |  IQR %.3f–%.3f  |  mean %.3f (95%% CI %.3f–%.3f)",
                rc$q1, rc$q3, rc$mean, rc$lo, rc$hi))
    else NULL
    corr_note <- if (isTRUE(attr(res, "sub_noise")))
      tags$span(style = "color:#888;",
        sprintf("  |  Rg σ_loc-corrected (σ_loc = %.2f nm)", attr(res, "sigma_loc")))
    else
      tags$span(style = "color:#888;", "  |  Rg not noise-corrected (set σ_loc to enable)")
    tags$div(style = "font-size:0.88em; margin:4px 0;",
      tags$b(sprintf("n = %d tracks", nrow(res))),
      tags$span(sprintf(
        "  |  median R\u2092/R\u2098\u2090\u02E3 = %.3f  |  confined (<0.35): %.0f%%  |  free (>0.45): %.0f%%",
        med_ratio, pct_conf, pct_free
      )),
      ci_note,
      corr_note
    )
  })

  output$trk_dl_conf_csv <- downloadHandler(
    filename = function() paste0("confinement_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) write.csv(conf_bundle(), file, row.names = FALSE)
  )

  # Observer: keep conf_color_map_stored in sync with confinement results.
  # Runs BEFORE the fire_frames observer so frames_store sees the fresh map.
  observe({
    if (!isTRUE(input$trk_conf_color_tracks)) {
      conf_color_map_stored(NULL)
      return()
    }
    res <- conf_bundle()
    if (is.null(res) || nrow(res) == 0) { conf_color_map_stored(NULL); return() }
    ramp    <- colorRamp(c("#C03030", "#808080", "#3060C0"))
    r       <- pmax(0, pmin(1, res$ratio))
    rgb_mat <- ramp(r)
    hex_cols <- rgb(rgb_mat[, 1], rgb_mat[, 2], rgb_mat[, 3], maxColorValue = 255)
    conf_color_map_stored(setNames(hex_cols, as.character(res$tid)))
  })

  # Rebuild frames (baking new colors in) when the confinement color toggle changes.
  observeEvent(input$trk_conf_color_tracks, {
    fire_frames()
  }, ignoreInit = TRUE)

  # ---- Directionality / turn-angle autocorrelation -------------------------
  dir_bundle <- eventReactive(input$trk_dir_run, {
    req(isTRUE(input$trk_show_dir))
    df <- trk(); req(nrow(df) > 0)
    min_pts <- max(4L, as.integer(input$trk_dir_min_pts))

    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
    df <- df[order(df$tid, df$tim), ]
    sp <- split(df, df$tid)

    # Pre-allocate angle storage as a list then unlist once at the end
    angle_list <- vector("list", length(sp))
    max_lag    <- 20L
    ac_sums    <- rep(0,  max_lag)
    ac_n       <- rep(0L, max_lag)
    # Localization-noise floor on a 2D step magnitude: RMS pure-noise step = 2·σ_loc.
    # Steps below it have noise-dominated directions (inflate the 180° turn-angle
    # peak and the lag-1 anti-correlation). Count them to quantify how much of the
    # directionality signal is noise. Uses the app-wide σ_loc.
    sig_dir <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    r_noise <- if (is.finite(sig_dir) && sig_dir > 0) 2 * sig_dir else NA_real_
    n_steps_tot <- 0L; n_below <- 0L

    withProgress(message = "Computing directionality…", value = 0, {
    for (ii in seq_along(sp)) {
      if (ii %% 50L == 0L) incProgress(50 / length(sp))
      tb <- sp[[ii]]
      if (nrow(tb) < min_pts) next
      dx <- diff(tb$x); dy <- diff(tb$y)
      if (length(dx) < 2) next

      # --- Turn angles: fully vectorised ---
      # mag of each displacement vector
      mag <- sqrt(dx^2 + dy^2)
      ok  <- mag > 0
      n_steps_tot <- n_steps_tot + length(mag)
      if (is.finite(r_noise)) n_below <- n_below + sum(mag < r_noise)
      # consecutive pairs where both vectors are non-zero
      ok_pair <- ok[-length(ok)] & ok[-1]
      if (any(ok_pair)) {
        dot   <- dx[-length(dx)] * dx[-1] + dy[-length(dy)] * dy[-1]
        denom <- mag[-length(mag)] * mag[-1]
        cos_a <- dot[ok_pair] / denom[ok_pair]
        cos_a <- pmax(-1, pmin(1, cos_a))
        angle_list[[ii]] <- acos(cos_a) * 180 / pi
      }

      # --- Velocity autocorrelation ---
      # unit vectors (only for non-zero steps)
      vx <- ifelse(ok, dx / mag, NA_real_)
      vy <- ifelse(ok, dy / mag, NA_real_)
      n_v <- length(vx)
      for (lag in seq_len(min(max_lag, n_v - 1L))) {
        i1   <- seq_len(n_v - lag)
        i2   <- i1 + lag
        keep <- !is.na(vx[i1]) & !is.na(vx[i2])
        if (!any(keep)) next
        dots <- vx[i1[keep]] * vx[i2[keep]] + vy[i1[keep]] * vy[i2[keep]]
        ac_sums[lag] <- ac_sums[lag] + sum(dots)
        ac_n[lag]    <- ac_n[lag]    + length(dots)
      }
    }

    all_angles <- unlist(angle_list, use.names = FALSE)

    shiny::validate(shiny::need(length(all_angles) >= 10,
      "Too few steps — lower 'Min points per TID' or relax filters."))

    ac_vals <- ifelse(ac_n > 0, ac_sums / ac_n, NA_real_)
    list(angles = all_angles, ac = ac_vals, lags = seq_len(max_lag),
         r_noise = r_noise, sigma_loc = if (is.finite(r_noise)) sig_dir else NA_real_,
         frac_below_noise = if (n_steps_tot > 0) n_below / n_steps_tot else NA_real_)
    })  # withProgress
  })

  output$trk_dir_angle_plot <- renderPlotly({
    bd   <- dir_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"

    h <- hist(bd$angles, breaks = seq(0, 180, by = 10), plot = FALSE)
    plot_ly() %>%
      add_bars(
        x = h$mids, y = h$density,
        marker = list(color = "rgba(80,140,210,0.65)",
                      line  = list(color = "rgba(50,100,170,0.8)", width = 0.5))
      ) %>%
      layout(
        xaxis = list(title = "Turn angle (\u00B0)", range = c(0, 180),
                     tickvals = seq(0, 180, 30), color = fg_col),
        yaxis = list(title = "Density", color = fg_col),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font = list(color = fg_col), showlegend = FALSE,
        uirevision = "dir_angle"
      )
  })

  output$trk_dir_autocorr_plot <- renderPlotly({
    bd   <- dir_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"

    plot_ly() %>%
      add_lines(x = bd$lags, y = bd$ac,
                line = list(color = "#D06030", width = 2)) %>%
      add_markers(x = bd$lags, y = bd$ac,
                  marker = list(color = "#D06030", size = 5)) %>%
      add_lines(x = range(bd$lags), y = c(0, 0),
                line = list(color = fg_col, width = 1, dash = "dot"),
                showlegend = FALSE) %>%
      layout(
        xaxis  = list(title = "Lag (steps)", color = fg_col),
        yaxis  = list(title = "Velocity autocorrelation", range = c(-1.05, 1.05),
                      color = fg_col),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font = list(color = fg_col), showlegend = FALSE,
        uirevision = "dir_ac"
      )
  })

  output$trk_dir_summary <- renderUI({
    bd <- dir_bundle()
    med_ang <- median(bd$angles, na.rm = TRUE)
    pct_fwd  <- mean(bd$angles <  60, na.rm = TRUE) * 100
    pct_back <- mean(bd$angles > 120, na.rm = TRUE) * 100
    ac1 <- if (length(bd$ac) >= 1 && is.finite(bd$ac[1])) round(bd$ac[1], 3) else NA
    noise_note <- if (is.finite(bd$frac_below_noise %||% NA_real_))
      tags$span(style = "color:#b26a00;", tags$br(),
        sprintf("⚠ %.0f%% of steps are below the localization-noise floor (2·σ_loc = %.1f nm) — their directions are noise-dominated and inflate the backward (180°) peak / lag-1 anti-correlation.",
                100 * bd$frac_below_noise, bd$r_noise))
    else
      tags$span(style = "color:#888;", tags$br(),
        "Set σ_loc (run the unified MSD, or enter it in Acquisition Parameters) to see how much of the turn-angle signal is localization noise.")
    tags$div(style = "font-size:0.88em; margin:4px 0;",
      tags$b(sprintf("n = %d turn angles", length(bd$angles))),
      tags$span(sprintf(
        "  |  median angle = %.1f\u00B0  |  forward (<60\u00B0): %.0f%%  |  backward (>120\u00B0): %.0f%%  |  AC(lag 1) = %s",
        med_ang, pct_fwd, pct_back, if (is.na(ac1)) "—" else ac1
      )),
      noise_note
    )
  })

  output$trk_dl_dir_csv <- downloadHandler(
    filename = function() paste0("directionality_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      bd <- dir_bundle()
      # Two sheets packed into one CSV with a blank-row separator
      ac_df  <- data.frame(lag_steps = bd$lags, velocity_autocorr = bd$ac)
      ang_df <- data.frame(turn_angle_deg = bd$angles)
      out <- rbind(
        cbind(sheet = "autocorrelation", ac_df,
              data.frame(turn_angle_deg = rep(NA, nrow(ac_df)))),
        data.frame(sheet = NA, lag_steps = NA, velocity_autocorr = NA,
                   turn_angle_deg = NA),
        cbind(sheet = "angles", lag_steps = NA, velocity_autocorr = NA,
              turn_angle_deg = ang_df$turn_angle_deg)
      )
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ---- Instantaneous speed (per-step v = |Δr|/Δt) -------------------------
  # Per-track speed time-series + rolling median. nm/ms = µm/s directly.
  trk_speed_one <- function(tb, smooth = 5L) {
    tb <- tb[order(tb$tim), ]
    t  <- as.numeric(tb$tim)
    has_z <- "z" %in% names(tb) && any(is.finite(tb$z))
    dt <- diff(t)
    dr <- sqrt(diff(tb$x)^2 + diff(tb$y)^2 + (if (has_z) diff(tb$z)^2 else 0))
    ok <- is.finite(dt) & dt > 0 & is.finite(dr)
    if (!any(ok)) return(NULL)
    spd   <- dr[ok] / dt[ok]            # nm/ms = µm/s
    t_rel <- t[-1][ok] - t[1]           # ms relative to track start
    k <- max(1L, as.integer(smooth)); if (k %% 2L == 0L) k <- k + 1L
    med <- if (length(spd) >= k && k > 1L) as.numeric(stats::runmed(spd, k)) else spd
    data.frame(t_ms = t_rel, speed_um_s = spd, speed_med = med)
  }

  # ---- Rolling drift (transport) velocity --------------------------------------
  # v_drift(t) = |net displacement over a time window W| / W.  Unlike per-step
  # path speed, this measures COHERENT transport: diffusion's net displacement
  # grows only as sqrt(t), so its apparent velocity ~ sqrt(4D/W) shrinks with W,
  # while directed transport gives a large, sustained v_drift. Windows are defined
  # by TIME (ms), robust to MINFLUX's irregular Δt, with a min-localizations guard.
  # nm/ms = µm/s directly. Returns midpoint-time, raw v_drift, and a light median.
  # min_span_frac guards the track ends: near the end the forward window can no
  # longer extend, so its DURATION collapses toward 0 while displacement noise
  # stays — inflating v = dr/span into a spurious spike (which otherwise reads as a
  # false "export at the very end"). Requiring the window to span >= frac·W drops
  # those unreliable end windows. Cost: the last <frac·W of each track is not
  # measured — you cannot estimate a velocity over near-zero time.
  trk_drift_velocity_one <- function(tb, win_ms = 200, min_win_pts = 5L, smooth = 5L,
                                     min_span_frac = 0.5) {
    tb <- tb[order(tb$tim), ]
    has_z <- "z" %in% names(tb) && any(is.finite(tb$z))
    # Drop non-finite rows: trk() may carry NA localizations, and the net-displacement
    # endpoints are used directly, so a single NA endpoint would poison v with NA
    # (and crash the downstream if(has_export)). Matches the rest of the tracking code.
    keep <- is.finite(tb$tim) & is.finite(tb$x) & is.finite(tb$y) &
            (if (has_z) is.finite(tb$z) else TRUE)
    tb <- tb[keep, , drop = FALSE]
    t  <- as.numeric(tb$tim)
    n  <- length(t)
    min_win_pts <- max(2L, as.integer(min_win_pts))
    if (n < min_win_pts || !is.finite(win_ms) || win_ms <= 0) return(NULL)
    min_span <- max(0, min_span_frac) * win_ms
    x <- tb$x; y <- tb$y; z <- if (has_z) tb$z else rep(0, n)
    t0 <- t[1]
    out_t <- numeric(n); out_v <- numeric(n); m <- 0L
    j <- 1L
    for (i in seq_len(n)) {
      if (j < i) j <- i
      while (j < n && t[j + 1L] <= t[i] + win_ms) j <- j + 1L   # last pt within [t_i, t_i+W]
      span <- t[j] - t[i]; npts <- j - i + 1L
      if (span >= min_span && span > 0 && npts >= min_win_pts) {
        dr <- sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2 + (z[j] - z[i])^2)
        m <- m + 1L
        out_t[m] <- (t[i] + t[j]) / 2 - t0   # centre of the window
        out_v[m] <- dr / span                # nm/ms = µm/s
      }
    }
    if (m < 1L) return(NULL)
    out_t <- out_t[seq_len(m)]; out_v <- out_v[seq_len(m)]
    k <- max(1L, as.integer(smooth)); if (k %% 2L == 0L) k <- k + 1L
    vmed <- if (length(out_v) >= k && k > 1L) as.numeric(stats::runmed(out_v, k)) else out_v
    data.frame(t_ms = out_t, v_drift = out_v, v_drift_med = vmed)
  }

  # ---- Local anomalous exponent α(t) -------------------------------------------
  # Rolling log-log slope of a short-window MSD: α≈1 Brownian, α→2 ballistic/
  # directed. This is the per-track analog of the ensemble TAMSD "shooting upward".
  # Windowed by STEP COUNT (needs contiguous lag structure), but the log-log fit uses
  # the REAL mean elapsed time per step-lag (τ = mean Δt over that lag), NOT the step
  # index — so MINFLUX's irregular Δt / dark cycles don't bias α (a window straddling a
  # dark gap would otherwise stretch a "1-step" lag over many cycles and flatten the
  # slope). Assigned at the window-centre time. Returns NULL if a track is too short.
  # sigma_nm > 0 subtracts the static localization-noise offset (2·d·σ²) from each
  # MSD lag before the log-log fit. Without it, noise flattens the short-lag MSD and
  # biases α downward (directed motion can read as α≈1 instead of ≈2).
  trk_local_alpha_one <- function(tb, win_steps = 21L, max_lag = 6L, sigma_nm = 0) {
    tb <- tb[order(tb$tim), ]
    has_z <- "z" %in% names(tb) && any(is.finite(tb$z))
    keep <- is.finite(tb$tim) & is.finite(tb$x) & is.finite(tb$y) &
            (if (has_z) is.finite(tb$z) else TRUE)
    tb <- tb[keep, , drop = FALSE]
    n  <- nrow(tb)
    win_steps <- max(7L, as.integer(win_steps)); if (win_steps %% 2L == 0L) win_steps <- win_steps + 1L
    if (n < win_steps) return(NULL)
    x <- tb$x; y <- tb$y; z <- if (has_z) tb$z else rep(0, n)
    t <- as.numeric(tb$tim); t0 <- t[1]
    d_dim     <- if (has_z) 3 else 2
    noise_off <- if (is.finite(sigma_nm) && sigma_nm > 0) 2 * d_dim * sigma_nm^2 else 0
    half <- win_steps %/% 2L
    L    <- max(2L, min(as.integer(max_lag), win_steps %/% 3L))
    centres <- (half + 1L):(n - half)
    out_t <- numeric(length(centres)); out_a <- rep(NA_real_, length(centres)); li <- 0L
    lags  <- seq_len(L)
    for (cc in centres) {
      idx <- (cc - half):(cc + half)
      xi <- x[idx]; yi <- y[idx]; zi <- z[idx]; ti <- t[idx]; mlen <- length(idx)
      msd <- rep(NA_real_, L); tau <- rep(NA_real_, L)   # MSD and REAL mean time-lag per step-lag
      for (lg in lags) {
        if (lg >= mlen) next
        dx  <- xi[(lg + 1L):mlen] - xi[1:(mlen - lg)]
        dy  <- yi[(lg + 1L):mlen] - yi[1:(mlen - lg)]
        dz  <- zi[(lg + 1L):mlen] - zi[1:(mlen - lg)]
        dtt <- ti[(lg + 1L):mlen] - ti[1:(mlen - lg)]    # actual elapsed time of each lag-lg pair
        v   <- mean(dx^2 + dy^2 + dz^2) - noise_off
        tt  <- mean(dtt[is.finite(dtt) & dtt > 0])       # mean real time spanned at this step-lag
        if (is.finite(v) && v > 0 && is.finite(tt) && tt > 0) { msd[lg] <- v; tau[lg] <- tt }
      }
      okl <- is.finite(msd) & is.finite(tau)
      li <- li + 1L
      out_t[li] <- t[cc] - t0
      if (sum(okl) >= 2L) {   # slope of log(MSD) vs log(REAL τ) = time-based α
        fit <- suppressWarnings(stats::lm.fit(cbind(1, log(tau[okl])), log(msd[okl])))
        out_a[li] <- unname(fit$coefficients[2])
      }
    }
    ok <- is.finite(out_a)
    if (!any(ok)) return(NULL)
    data.frame(t_ms = out_t[ok], alpha = out_a[ok])
  }

  # Effective threshold for calling a drift-velocity window "directed":
  #   absolute onset_thr if the user set one, else fold × the track's own diffusive
  #   baseline (median drift velocity ~ sqrt(4D/W)), always floored at the
  #   localization-noise floor √2σ/W. Fold-over-baseline (not median+k·MAD of the
  #   signal) is what keeps pure-diffusion tracks from being flagged: validated on
  #   synthetic data, fold≈3 gives ~5% false positives with full sensitivity.
  trk_export_thr <- function(v, onset_thr, fold, vfloor) {
    base <- stats::median(v, na.rm = TRUE)
    thr  <- if (is.finite(onset_thr) && onset_thr > 0) onset_thr else fold * base
    max(thr, vfloor, na.rm = TRUE)
  }

  # Effective σ_loc for the speed/transport noise floor: the user's field if set
  # (> 0), else the app-wide value (unified MSD / manual), else 0. Keeps the drift
  # noise floor √2σ/W and the α de-bias consistent with the rest of the app.
  trk_speed_sigma_eff <- reactive({
    v <- suppressWarnings(as.numeric(input$trk_speed_sigma_nm %||% 0))
    if (is.finite(v) && v > 0) return(v)
    s <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    if (is.finite(s) && s > 0) s else 0
  })

  # ---- Per-track export quantification -----------------------------------------
  # One summary row per qualifying track. An "export" = a run of >= min_run
  # consecutive drift-velocity windows above the effective threshold (max of the
  # user threshold / auto median+3·MAD, and the localization-noise floor √2σ/W).
  # Requiring a sustained run — not a single peak — avoids flagging noise spikes.
  # onset_frac = position of the export onset along the track (0 start, 1 end),
  # which is what reveals the "exports cluster at track ends" pattern.
  trk_export_summary <- eventReactive(input$trk_speed_run, {
    req(isTRUE(input$trk_show_speed))
    df <- trk(); req(nrow(df) > 0)
    min_pts   <- max(5L, as.integer(input$trk_speed_min_pts %||% 20L))
    smooth    <- as.integer(input$trk_speed_smooth %||% 5L)
    drift_win <- as.numeric(input$trk_speed_drift_win %||% 200)
    drift_mp  <- as.integer(input$trk_speed_drift_minpts %||% 5L)
    sigma_nm  <- trk_speed_sigma_eff()
    onset_thr <- as.numeric(input$trk_speed_onset_thr %||% 0)
    fold      <- as.numeric(input$trk_speed_fold %||% 3)
    min_run   <- max(1L, as.integer(input$trk_speed_min_run %||% 3L))
    guard_win <- max(0, as.numeric(input$trk_speed_edge_guard %||% 1))
    a_win     <- as.integer(input$trk_speed_alpha_win %||% 21L)
    a_maxlag  <- as.integer(input$trk_speed_alpha_maxlag %||% 6L)
    vfloor    <- if (is.finite(sigma_nm) && sigma_nm > 0 && drift_win > 0) sqrt(2) * sigma_nm / drift_win else 0

    counts <- table(df$tid); tids <- names(counts[counts >= min_pts])
    shiny::validate(shiny::need(length(tids) > 0, "No tracks with enough points — lower 'Min points per TID'."))

    rows <- withProgress(message = "Speed & transport…", value = 0,
      lapply(seq_along(tids), function(ti_sp) {
      if (ti_sp %% 25L == 0L) incProgress(25 / length(tids), detail = "per-track drift/α")
      t_i <- tids[ti_sp]
      d  <- df[df$tid == t_i, , drop = FALSE]
      dv <- trk_drift_velocity_one(d, drift_win, drift_mp, smooth)
      if (is.null(dv) || nrow(dv) < max(2L, min_run)) return(NULL)
      al <- trk_local_alpha_one(d, a_win, a_maxlag, sigma_nm)
      v  <- dv$v_drift_med
      thr_eff <- trk_export_thr(v, onset_thr, fold, vfloor)
      # Edge guard: the first/last window of a track is unreliable (initial
      # localization jitter at track start, collapsing span at the end), which
      # otherwise piles false "export at frac 0/1" events at the edges. Restrict
      # detection to interior windows >= guard_win·W from each end.
      tmin <- min(dv$t_ms); tmax <- max(dv$t_ms); t_span <- tmax - tmin
      guard    <- guard_win * drift_win
      interior <- dv$t_ms >= (tmin + guard) & dv$t_ms <= (tmax - guard)
      n_int    <- sum(interior)
      above <- (v >= thr_eff) & interior
      r <- rle(above); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
      seg_ok <- r$values & r$lengths >= min_run
      has_exp <- isTRUE(any(seg_ok))          # NA-proof
      n_events <- sum(seg_ok, na.rm = TRUE)   # number of directed episodes (hops/exports)
      # "Sits still" fraction: interior windows whose net motion is at/below the
      # localization-noise floor — i.e. confined/immobile at this timescale. High
      # immobile_frac + several n_events = the hop-diffusion (dwell → hop → dwell) look.
      immobile_frac <- if (n_int > 0) sum(interior & v <= max(vfloor, 0), na.rm = TRUE) / n_int else NA_real_
      onset_ms <- NA_real_; onset_frac <- NA_real_; peak_frac <- NA_real_
      export_v <- NA_real_; export_dur <- NA_real_
      if (has_exp) {
        first    <- which(seg_ok)[1]
        onset_ms <- dv$t_ms[starts[first]]
        onset_frac <- if (is.finite(t_span) && t_span > 0) (onset_ms - tmin) / t_span else NA_real_
        keep <- unlist(mapply(function(s, e) s:e, starts[seg_ok], ends[seg_ok], SIMPLIFY = FALSE))
        # Export POSITION = the strongest directed window (peak drift velocity),
        # a far less start-biased locator than the first threshold crossing.
        peak_idx  <- keep[which.max(v[keep])]
        peak_frac <- if (is.finite(t_span) && t_span > 0) (dv$t_ms[peak_idx] - tmin) / t_span else NA_real_
        export_v   <- stats::median(v[keep])
        export_dur <- if (n_int > 0) sum(r$lengths[seg_ok]) / n_int * t_span else NA_real_
      }
      data.frame(
        tid = t_i, n_loc = nrow(d), track_dur_ms = t_span,
        v_drift_peak = max(v), v_drift_median = stats::median(v),
        directed_frac = if (n_int > 0) sum(above, na.rm = TRUE) / n_int else 0,
        immobile_frac = immobile_frac,
        alpha_max = if (!is.null(al)) max(al$alpha) else NA_real_,
        has_export = has_exp, n_events = n_events,
        export_onset_ms = onset_ms, export_onset_frac = onset_frac,
        export_pos_frac = peak_frac,
        export_v = export_v, export_dur_ms = export_dur,
        stringsAsFactors = FALSE
      )
    }))
    out <- do.call(rbind, rows)
    shiny::validate(shiny::need(!is.null(out) && nrow(out) > 0, "No tracks produced a drift-velocity curve."))
    out
  })

  output$trk_speed_tid_ui <- renderUI({
    req(isTRUE(input$trk_show_speed), speed_ran())
    s <- tryCatch(trk_export_summary(), error = function(e) NULL)
    req(!is.null(s), nrow(s) > 0)
    # Rank tracks by export strength so the picker doubles as a "top exporters"
    # follow-up list: exporters first, then by peak drift velocity (descending).
    s <- s[order(-as.integer(s$has_export %in% TRUE), -s$v_drift_peak), , drop = FALSE]
    labs <- sprintf("%s%s  (peak %.2g µm/s)",
                    ifelse(s$has_export %in% TRUE, "★ ", "  "),
                    s$tid, s$v_drift_peak)
    choices <- stats::setNames(s$tid, labs)
    sel <- if (!is.null(input$trk_speed_tid) && input$trk_speed_tid %in% s$tid)
             input$trk_speed_tid else s$tid[1]
    selectInput("trk_speed_tid", "Track (TID) — ★ = flagged exporter, ranked by peak drift",
                choices = choices, selected = sel)
  })

  output$trk_speed_plot <- renderPlotly({
    req(isTRUE(input$trk_show_speed), speed_ran())
    df <- trk(); req(nrow(df) > 0)
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    min_pts <- max(5L, as.integer(input$trk_speed_min_pts %||% 20L))
    smooth  <- as.integer(input$trk_speed_smooth %||% 5L)
    counts  <- table(df$tid); tids <- names(counts[counts >= min_pts])
    shiny::validate(shiny::need(length(tids) > 0, "No tracks with enough points — lower 'Min points per TID'."))
    tid_sel <- input$trk_speed_tid
    if (is.null(tid_sel) || !(tid_sel %in% tids)) tid_sel <- tids[1]

    # New detectors (default on): rolling drift velocity + local exponent α(t).
    show_path  <- isTRUE(input$trk_speed_show_path  %||% TRUE)
    show_drift <- isTRUE(input$trk_speed_show_drift %||% TRUE)
    show_alpha <- isTRUE(input$trk_speed_show_alpha %||% TRUE)
    drift_win  <- as.numeric(input$trk_speed_drift_win %||% 200)
    drift_mp   <- as.integer(input$trk_speed_drift_minpts %||% 5L)
    sigma_nm   <- trk_speed_sigma_eff()
    a_win      <- as.integer(input$trk_speed_alpha_win %||% 21L)
    a_maxlag   <- as.integer(input$trk_speed_alpha_maxlag %||% 6L)
    onset_thr  <- as.numeric(input$trk_speed_onset_thr %||% 0)
    fold       <- as.numeric(input$trk_speed_fold %||% 3)

    fig <- plot_ly()

    dsel <- df[df$tid == tid_sel, , drop = FALSE]
    spd  <- trk_speed_one(dsel, smooth)
    dv   <- trk_drift_velocity_one(dsel, drift_win, drift_mp, smooth)
    al   <- trk_local_alpha_one(dsel, a_win, a_maxlag, sigma_nm)
    shiny::validate(shiny::need(!is.null(spd) && nrow(spd) >= 2, "Selected track has too few valid steps."))

    # Path speed (per-step) — kept for reference; it is dominated by diffusion+noise
    if (show_path) {
      fig <- fig %>%
        add_trace(x = spd$t_ms, y = spd$speed_um_s, type = "scattergl", mode = "markers",
                  marker = list(size = 4, color = "rgba(140,140,140,0.30)"),
                  name = "path speed (raw)", hoverinfo = "x+y") %>%
        add_trace(x = spd$t_ms, y = spd$speed_med, type = "scattergl", mode = "lines",
                  line = list(color = "rgba(120,120,120,0.85)", width = 1.5, dash = "dot"),
                  name = "path speed (median)")
    }

    # Rolling drift (transport) velocity — the export-detecting curve
    if (show_drift && !is.null(dv) && nrow(dv) >= 2) {
      fig <- fig %>%
        add_trace(x = dv$t_ms, y = dv$v_drift_med, type = "scattergl", mode = "lines",
                  line = list(color = "#1B9E5A", width = 3), name = "drift velocity")
      vfloor <- if (is.finite(sigma_nm) && sigma_nm > 0 && is.finite(drift_win) && drift_win > 0)
                  sqrt(2) * sigma_nm / drift_win else 0
      # Localization-noise floor √2·σ/W (µm/s): drift bumps below this are noise.
      if (vfloor > 0) {
        xr <- range(dv$t_ms)
        fig <- fig %>% add_trace(x = xr, y = c(vfloor, vfloor), type = "scatter", mode = "lines",
          line = list(color = "rgba(200,120,0,0.9)", width = 1.2, dash = "dash"),
          name = "noise floor √2σ/W", hoverinfo = "name")
      }
      # Onset overlay: same threshold as the population classifier (fold × baseline,
      # or absolute if set), highlight windows above it.
      thr <- trk_export_thr(dv$v_drift_med, onset_thr, fold, vfloor)
      if (is.finite(thr) && thr > 0) {
        xr <- range(dv$t_ms)
        above <- dv$v_drift_med >= thr
        fig <- fig %>% add_trace(x = xr, y = c(thr, thr), type = "scatter", mode = "lines",
          line = list(color = "rgba(224,48,32,0.85)", width = 1.2, dash = "dash"),
          name = "onset threshold", hoverinfo = "name")
        if (any(above))
          fig <- fig %>% add_trace(x = dv$t_ms[above], y = dv$v_drift_med[above],
            type = "scattergl", mode = "markers",
            marker = list(size = 6, color = "#E03020"),
            name = "above threshold", hoverinfo = "x+y")
      }
    }

    # Local exponent α(t) on a secondary axis (0..2): 1 = Brownian, 2 = directed
    if (show_alpha && !is.null(al) && nrow(al) >= 2) {
      fig <- fig %>% add_trace(x = al$t_ms, y = al$alpha, type = "scattergl", mode = "lines",
        line = list(color = "#7B4FBE", width = 2), name = "local α (right axis)",
        yaxis = "y2", hoverinfo = "x+y")
    }

    # α reference lines (1 = Brownian, 2 = ballistic) on the secondary axis
    a_shapes <- list()
    if (show_alpha && !is.null(al) && nrow(al) >= 2) {
      for (av in c(1, 2))
        a_shapes[[length(a_shapes) + 1L]] <- list(type = "line", xref = "paper",
          x0 = 0, x1 = 1, yref = "y2", y0 = av, y1 = av,
          line = list(color = "rgba(123,79,190,0.35)", width = 1, dash = "dot"))
    }

    fig %>% layout(
      xaxis = list(title = "time within track (ms)", color = fg),
      yaxis = list(title = "velocity (µm/s)", color = fg, rangemode = "tozero"),
      yaxis2 = list(title = "local α", overlaying = "y", side = "right",
                    range = c(0, 2.1), color = "#7B4FBE", showgrid = FALSE),
      shapes = a_shapes,
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(font = list(color = fg, size = 10), orientation = "h", y = -0.2),
      title = list(text = sprintf("TID %s", tid_sel), font = list(color = fg, size = 12)),
      uirevision = "speed")
  })

  # ---- Population view: quantify export across all tracks -----------------------
  output$trk_export_summary_txt <- renderUI({
    req(isTRUE(input$trk_show_speed))
    s <- trk_export_summary()
    n  <- nrow(s); ne <- sum(s$has_export, na.rm = TRUE)
    med_peak <- stats::median(s$v_drift_peak, na.rm = TRUE)
    exp_rows <- s[s$has_export %in% TRUE, , drop = FALSE]
    med_onset <- if (nrow(exp_rows)) stats::median(exp_rows$export_pos_frac, na.rm = TRUE) else NA_real_
    med_expv  <- if (nrow(exp_rows)) stats::median(exp_rows$export_v, na.rm = TRUE) else NA_real_
    n_multi   <- if (nrow(exp_rows)) sum(exp_rows$n_events >= 2, na.rm = TRUE) else 0
    htmltools::div(style = "font-size:0.9em; line-height:1.5;",
      htmltools::HTML(sprintf(
        "<b>%d</b> tracks analysed &middot; <b style='color:#1B9E5A;'>%d (%.0f%%)</b> flagged as exporters%s<br>%s%smedian peak drift = <b>%.3g µm/s</b>",
        n, ne, if (n > 0) 100 * ne / n else 0,
        if (n_multi > 0) sprintf(" &middot; <b>%d</b> with multiple episodes (hop/multi-step)", n_multi) else "",
        if (is.finite(med_expv))  sprintf("median export velocity = <b>%.3g µm/s</b> &middot; ", med_expv) else "",
        if (is.finite(med_onset)) sprintf("median export at <b>%.0f%%</b> of track length &middot; ", 100 * med_onset) else "",
        med_peak)))
  })

  # Distribution of a per-track metric (the "histogram of speeds" across tracks)
  output$trk_export_hist <- renderPlotly({
    req(isTRUE(input$trk_show_speed))
    s <- trk_export_summary()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    metric <- input$trk_export_hist_metric %||% "v_drift_peak"
    lab <- switch(metric,
      v_drift_peak  = "peak drift velocity (µm/s)",
      export_v      = "export-segment velocity (µm/s)",
      directed_frac = "fraction of track time directed",
      immobile_frac = "fraction of track time immobile/confined",
      n_events      = "number of transport episodes per track",
      alpha_max     = "max local α",
      "value")
    # export_v exists only for exporters; other metrics span all tracks
    vals_all <- s[[metric]]
    exp_mask <- s$has_export %in% TRUE
    fig <- plot_ly()
    if (metric %in% c("v_drift_peak", "directed_frac", "immobile_frac", "n_events", "alpha_max")) {
      # non-exporters vs exporters stacked so the split is visible
      ve <- vals_all[exp_mask];  vn <- vals_all[!exp_mask]
      ve <- ve[is.finite(ve)];   vn <- vn[is.finite(vn)]
      if (length(vn)) fig <- fig %>% add_histogram(x = vn, name = "non-exporter",
        marker = list(color = "rgba(150,150,150,0.55)"), nbinsx = 30)
      if (length(ve)) fig <- fig %>% add_histogram(x = ve, name = "exporter",
        marker = list(color = "rgba(27,158,90,0.75)"), nbinsx = 30)
      fig <- fig %>% layout(barmode = "stack")
    } else {
      v <- vals_all[exp_mask]; v <- v[is.finite(v)]
      shiny::validate(shiny::need(length(v) > 0, "No exporter tracks yet — loosen the threshold / min-run."))
      fig <- fig %>% add_histogram(x = v, name = "exporters",
        marker = list(color = "rgba(27,158,90,0.75)"), nbinsx = 30)
    }
    fig %>% layout(
      xaxis = list(title = lab, color = fg),
      yaxis = list(title = "tracks", color = fg),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(font = list(color = fg, size = 10), orientation = "h", y = -0.25),
      title = list(text = "Per-track distribution", font = list(color = fg, size = 12)),
      uirevision = "exporthist")
  })

  # Where along the track is transport strongest? (0 = start, 1 = end)
  # Uses the peak-drift-velocity position (not first crossing, which is start-biased).
  output$trk_export_onset_hist <- renderPlotly({
    req(isTRUE(input$trk_show_speed))
    s <- trk_export_summary()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    v <- s$export_pos_frac[s$has_export %in% TRUE]; v <- v[is.finite(v)]
    shiny::validate(shiny::need(length(v) > 0, "No exporter tracks yet — loosen the threshold / min-run."))
    plot_ly(x = v, type = "histogram", nbinsx = 20,
            marker = list(color = "rgba(123,79,190,0.75)")) %>%
      layout(
        xaxis = list(title = "export position (0 = start … 1 = end of track)",
                     range = c(0, 1), color = fg),
        yaxis = list(title = "exporter tracks", color = fg),
        paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
        title = list(text = sprintf("Export location (n=%d) — pile-up near 1 = exports at track ends", length(v)),
                     font = list(color = fg, size = 12)),
        uirevision = "onsethist")
  })

  # Per-track summary table download (the quantification deliverable)
  output$trk_dl_export_csv <- downloadHandler(
    filename = function() paste0("export_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) write.csv(trk_export_summary(), file, row.names = FALSE)
  )

  output$trk_dl_speed_csv <- downloadHandler(
    filename = function() paste0("speed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      df <- tryCatch(trk(), error = function(e) NULL); req(!is.null(df))
      min_pts   <- max(5L, as.integer(input$trk_speed_min_pts %||% 20L))
      smooth    <- as.integer(input$trk_speed_smooth %||% 5L)
      drift_win <- as.numeric(input$trk_speed_drift_win %||% 200)
      drift_mp  <- as.integer(input$trk_speed_drift_minpts %||% 5L)
      a_win     <- as.integer(input$trk_speed_alpha_win %||% 21L)
      a_maxlag  <- as.integer(input$trk_speed_alpha_maxlag %||% 6L)
      sigma_nm  <- trk_speed_sigma_eff()
      counts  <- table(df$tid); tids <- names(counts[counts >= min_pts])
      # Long/tidy format: all three metrics stacked (they have different time bases).
      out <- do.call(rbind, lapply(tids, function(t_i) {
        d <- df[df$tid == t_i, , drop = FALSE]
        sp <- trk_speed_one(d, smooth)
        dv <- trk_drift_velocity_one(d, drift_win, drift_mp, smooth)
        al <- trk_local_alpha_one(d, a_win, a_maxlag, sigma_nm)
        rows <- list()
        if (!is.null(sp)) rows[[length(rows)+1L]] <- data.frame(tid = t_i, metric = "path_speed_um_s",
          t_ms = sp$t_ms, value = sp$speed_um_s, value_smoothed = sp$speed_med)
        if (!is.null(dv)) rows[[length(rows)+1L]] <- data.frame(tid = t_i, metric = "drift_velocity_um_s",
          t_ms = dv$t_ms, value = dv$v_drift, value_smoothed = dv$v_drift_med)
        if (!is.null(al)) rows[[length(rows)+1L]] <- data.frame(tid = t_i, metric = "local_alpha",
          t_ms = al$t_ms, value = al$alpha, value_smoothed = NA_real_)
        if (!length(rows)) return(NULL)
        do.call(rbind, rows)
      }))
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ---- Diffusion State HMM (Baum-Welch EM + Viterbi) ----------------------
  hmm_bundle <- eventReactive(input$trk_hmm_run, {
    req(isTRUE(input$trk_show_hmm))
    df <- trk(); req(nrow(df) > 0)
    K        <- max(2L, min(4L, as.integer(input$trk_hmm_n_states)))
    min_pts  <- max(4L, as.integer(input$trk_hmm_min_pts))
    n_ex     <- max(1L, as.integer(input$trk_hmm_n_ex))
    norm_dt  <- isTRUE(input$trk_hmm_norm_dt)

    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
    df <- df[order(df$tid, df$tim), ]
    sp <- split(df, df$tid)
    sp <- sp[sapply(sp, nrow) >= min_pts]
    shiny::validate(shiny::need(length(sp) > 0,
      "No tracks meet the minimum length — lower 'Min localizations per TID'."))

    # Build per-track observation sequences.
    # When norm_dt = TRUE each step r is divided by Ã¢ÂÂÃÂt so that states correspond
    # to diffusion coefficients rather than raw step sizes (important for irregular
    # MINFLUX sampling where variable ÃÂt within a track would otherwise create
    # spurious apparent state switches).
    obs_list <- Filter(function(r) length(r) >= 2L,
      lapply(sp, function(tb) {
        r  <- sqrt(diff(tb$x)^2 + diff(tb$y)^2)
        if (norm_dt) {
          dt <- diff(tb$tim)
          ok <- is.finite(dt) & dt > 0 & is.finite(r) & r > 0
          r  <- ifelse(ok, r / sqrt(dt), NA_real_)
          r  <- r[!is.na(r)]
        } else {
          r <- r[is.finite(r) & r > 0]
        }
        r
      })
    )
    all_r <- unlist(obs_list, use.names = FALSE)
    shiny::validate(shiny::need(length(all_r) >= 20, "Too few valid steps."))

    # ---- HMM helpers --------------------------------------------------------
    lse <- function(x) {
      x <- x[is.finite(x)]; if (!length(x)) return(-Inf)
      m <- max(x); m + log(sum(exp(x - m)))
    }
    # Fast logsumexp over columns of a KÃÃxK matrix; returns -Inf (not NaN) for all-(-Inf) columns
    # Column-wise max of a matrix without apply() (K small -> 1-3 pmax calls).
    col_max <- function(mat) {
      nr <- nrow(mat); m <- mat[1L, ]
      if (nr > 1L) for (r in 2L:nr) m <- pmax(m, mat[r, ])
      m
    }
    # Column-wise logsumexp; returns -Inf (not NaN) for all-(-Inf) columns.
    # Vectorised replacement for the old apply()/sweep() version (identical result).
    lse_cols <- function(mat) {
      nr  <- nrow(mat)
      m   <- col_max(mat)
      fin <- is.finite(m)
      out <- rep(-Inf, length(m))
      if (any(fin)) {
        mm <- m; mm[!fin] <- 0
        cs <- colSums(exp(mat - rep(mm, each = nr)))
        out[fin] <- m[fin] + log(cs[fin])
      }
      out
    }
    # Row-wise logsumexp for a T x K matrix (replaces apply(mat, 1L, lse)).
    lse_rows <- function(mat) {
      nc  <- ncol(mat); m <- mat[, 1L]
      if (nc > 1L) for (cc in 2L:nc) m <- pmax(m, mat[, cc])
      fin <- is.finite(m)
      out <- rep(-Inf, nrow(mat))
      if (any(fin)) {
        mm <- m; mm[!fin] <- 0
        rs <- rowSums(exp(mat - mm))       # mm (length T) recycles down columns
        out[fin] <- m[fin] + log(rs[fin])
      }
      out
    }
    # Elementwise log(exp(a) + exp(b)) for equally-shaped matrices/vectors.
    log_add <- function(a, b) {
      m   <- pmax(a, b)
      res <- m + log1p(exp(-abs(a - b)))
      res[!is.finite(m)] <- -Inf          # both -Inf -> -Inf (not NaN)
      res
    }
    # Log Rayleigh emission: TÃÃxK matrix
    lem <- function(r, sigma)
      outer(r, sigma, function(ri, si) log(ri) - 2*log(si) - ri^2/(2*si^2))

    fwd <- function(le, lA, lpi) {
      T <- nrow(le); la <- matrix(-Inf, T, K)
      la[1L, ] <- lpi + le[1L, ]
      for (t in seq(2L, T))
        la[t, ] <- lse_cols(lA + la[t-1L, ]) + le[t, ]
      la
    }
    bwd <- function(le, lA) {
      T <- nrow(le); lb <- matrix(0, T, K)
      tlA <- t(lA)                          # hoist transpose out of the loop
      for (t in seq(T-1L, 1L)) {
        v <- le[t+1L, ] + lb[t+1L, ]
        lb[t, ] <- lse_cols(tlA + v)   # tlA[j,k]=lA[k,j]; +v[j] via recycle -> col k = lb[t,k]
      }
      lb
    }
    safe_which_max <- function(x) {
      # Replace NaN/NA with -Inf so which.max never returns integer(0)
      x[!is.finite(x)] <- -Inf
      w <- which.max(x)
      if (!length(w)) 1L else w
    }
    viterbi_path <- function(obs, lA, lpi, sigma) {
      le <- lem(obs, sigma)
      # Replace any -Inf rows with uniform log-prob so decoding never collapses
      bad_rows <- !is.finite(rowSums(le))
      if (any(bad_rows)) le[bad_rows, ] <- log(1 / K)
      T <- nrow(le)
      delta <- matrix(-Inf, T, K); psi <- matrix(1L, T, K)
      delta[1L, ] <- lpi + le[1L, ]
      for (t in seq(2L, T)) {
        mat        <- lA + delta[t-1L, ]
        mat[!is.finite(mat)] <- -Inf
        # column argmax via max.col(t(mat)); ties.method="first" matches which.max
        psi[t, ]   <- max.col(t(mat), ties.method = "first")
        delta[t, ] <- col_max(mat) + le[t, ]
      }
      path <- integer(T); path[T] <- safe_which_max(delta[T, ])
      for (t in seq(T-1L, 1L)) path[t] <- psi[t+1L, path[t+1L]]
      path
    }

    # ---- Initialisation -----------------------------------------------------
    # Sigma for a Rayleigh is close to the median / sqrt(ln 4) ÃÂ¢Ã¢ÂÂ°ÃÂ median / 1.18.
    # Spread K initial sigmas evenly across 10thÃÂ¢Ã¢ÂÂ¬Ã¢ÂÂ90th percentile of step sizes.
    # Floor at 1% of the overall median (at least 0.5 nm) to prevent collapse.
    r_floor <- max(0.5, median(all_r) * 0.01)
    sigma_init_pts <- quantile(all_r, seq(0.10, 0.90, length.out = K))
    sigma <- pmax(sort(as.numeric(sigma_init_pts) / sqrt(log(4))), r_floor)
    pi_v <- rep(1/K, K)
    A    <- matrix(0.1 / max(K-1L, 1L), K, K); diag(A) <- 0.9

    # Fit on at most max_fit_trk longest tracks to keep runtime bounded
    max_fit_trk <- max(10L, as.integer(isolate(input$trk_hmm_max_fit_trk) %||% 100L))
    fit_obs <- obs_list[order(sapply(obs_list, length), decreasing = TRUE)[seq_len(min(max_fit_trk, length(obs_list)))]]

    # Optional per-track step cap for FITTING only (0 = use every step). The EM
    # cost is linear in the total number of steps, and the longest tracks are the
    # slowest; a contiguous prefix of a long track already contains ample emission
    # + transition statistics, so capping barely moves the fitted σ/A but can cut
    # runtime several-fold. A CONTIGUOUS prefix is used (not every-Nth) because the
    # observation is a per-step displacement — subsampling localizations would
    # inflate the step sizes and bias the diffusion estimate. The Viterbi decode
    # below always runs on the FULL tracks, so per-localization state resolution
    # is unaffected.
    fit_step_cap <- as.integer(isolate(input$trk_hmm_max_fit_steps) %||% 0L)
    if (is.finite(fit_step_cap) && fit_step_cap >= 10L)
      fit_obs <- lapply(fit_obs, function(o)
        if (length(o) > fit_step_cap) o[seq_len(fit_step_cap)] else o)

    # ---- Baum-Welch EM -------------------------------------------------------
    # Speed: (1) ≤30 iterations with an early-stop convergence threshold;
    # (2) all per-step log-space math is vectorised (col_max/lse_cols/lse_rows/
    #     log_add) — no apply()/sweep() and no K×K lse() loop in the hot path;
    # (3) fitting is bounded by max_fit_trk (tracks) and the optional step cap
    #     above. The remaining cost is the inherently sequential forward/backward
    #     recursion (one iteration per localization); shrink tracks/steps to speed
    #     it up further, or move the recursion to Rcpp for a step change.
    prev_ll <- -Inf
    withProgress(message = "Fitting HMM (Baum-Welch)…", value = 0, {
      for (iter in seq_len(30L)) {
        setProgress(iter / 30)
        lA <- log(A); lpi <- log(pi_v)
        sg1  <- rep(0, K); sg <- rep(0, K); sgr2 <- rep(0, K); total_ll <- 0
        log_sx <- matrix(-Inf, K, K)

        for (obs in fit_obs) {
          le  <- lem(obs, sigma)
          la  <- fwd(le, lA, lpi)
          lb  <- bwd(le, lA)
          ll  <- lse(la[nrow(la), ])
          if (!is.finite(ll)) next
          total_ll <- total_ll + ll

          lg  <- la + lb; lg <- lg - lse_rows(lg)
          gam <- exp(lg)
          bad <- !is.finite(rowSums(gam))
          if (any(bad)) gam[bad, ] <- 1 / K
          sg1  <- sg1  + gam[1L, ]
          sg   <- sg   + colSums(gam)
          sgr2 <- sgr2 + colSums(gam * obs^2)

          # Vectorised log-xi: broadcast v2_t[,j] over all K columns of la_t,
          # then colwise logsumexp — replaces the inner K*K lse() loop.
          T_obs <- length(obs)
          la_t  <- la[seq_len(T_obs - 1L), , drop = FALSE]
          v2_t  <- le[seq(2L, T_obs), , drop = FALSE] +
                   lb[seq(2L, T_obs), , drop = FALSE]
          log_inner <- matrix(-Inf, K, K)
          nrt <- nrow(la_t)
          for (j in seq_len(K)) {
            m     <- la_t + v2_t[, j]          # (T-1) x K broadcast
            m_max <- col_max(m)
            ok    <- is.finite(m_max)
            if (any(ok)) {
              cs <- colSums(exp(m[, ok, drop = FALSE] - rep(m_max[ok], each = nrt)))
              log_inner[ok, j] <- m_max[ok] + log(pmax(cs, .Machine$double.eps))
            }
          }
          log_xi_new <- lA + log_inner - ll
          # Accumulate in log-space, vectorised over the whole K x K matrix.
          log_sx <- log_add(log_sx, log_xi_new)
        }

        if (any(!is.finite(c(sg, sgr2)))) break
        pi_v <- sg1 / sum(sg1)
        for (i in seq_len(K)) {
          row_lse <- lse(log_sx[i, ])
          A[i, ]  <- if (is.finite(row_lse)) exp(log_sx[i, ] - row_lse) else rep(1/K, K)
        }
        sigma <- sqrt(pmax(sgr2 / (2 * pmax(sg, .Machine$double.eps)), r_floor^2))

        if (is.finite(total_ll) && is.finite(prev_ll) &&
            abs(total_ll - prev_ll) < 1e-3 * abs(prev_ll) + 1e-5) break
        if (is.finite(total_ll)) prev_ll <- total_ll
      }
    })

    # Sort states slow -> fast (by sigma)
    ord   <- order(sigma)
    sigma <- sigma[ord]; pi_v <- pi_v[ord]; A <- A[ord, ][, ord]

    # ---- Viterbi decode all tracks ------------------------------------------
    # Observations must use the same normalisation as training.
    # viterbi_path patches degenerate emission rows with uniform log-prob,
    # so NA/Inf entries (from dt Ã¢ÂÂ¤ 0 when norm_dt=TRUE) decode as the most
    # likely state given context rather than crashing.
    lA <- log(A); lpi <- log(pi_v)
    seg_list <- lapply(names(obs_list), function(nm) {
      tb  <- sp[[nm]]
      if (norm_dt) {
        r_raw <- sqrt(diff(tb$x)^2 + diff(tb$y)^2)
        dt_v  <- diff(tb$tim)
        ok    <- is.finite(dt_v) & dt_v > 0 & is.finite(r_raw)
        obs   <- ifelse(ok, r_raw / sqrt(dt_v), NA_real_)
      } else {
        obs <- sqrt(diff(tb$x)^2 + diff(tb$y)^2)
      }
      # Replace NA/non-positive with a small floor so viterbi_path can patch them
      obs[!is.finite(obs) | obs <= 0] <- NA_real_
      path <- viterbi_path(obs, lA, lpi, sigma)
      data.frame(tid = nm, x = tb$x, y = tb$y, tim = tb$tim,
                 state = c(path, tail(path, 1L)), stringsAsFactors = FALSE)
    })

    all_states <- unlist(lapply(seg_list, `[[`, "state"), use.names = FALSE)
    fracs      <- tabulate(all_states, nbins = K) / length(all_states)

    # Median ÃÂt across all tracks (ms)
    all_dt <- unlist(lapply(sp, function(tb) diff(tb$tim)), use.names = FALSE)
    dt_med <- median(all_dt[all_dt > 0], na.rm = TRUE)

    # D computation:
    # norm_dt=TRUE:  ÃÂ is in nm/Ã¢ÂÂms  Ã¢ÂÂ  D = ÃÂ²/2  in nm²/ms (proper diffusivity)
    # norm_dt=FALSE: ÃÂ is in nm      Ã¢ÂÂ  D = ÃÂ²/2  in nm²/step (not nm²/ms)
    D_vals <- sigma^2 / 2

    # Flag states that are indistinguishable from the localization-noise floor.
    # A pure-noise 2D step has per-coordinate variance 2σ_loc², i.e. an apparent
    # D of σ_loc²/Δt (normalised) or σ_loc² per step (un-normalised). A fitted
    # state at or below this is a noise/immobile-limited state, not a resolved
    # diffusive one. σ_loc is the app-wide value (unified MSD fit or manual).
    sig_loc <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    D_noise <- if (is.finite(sig_loc) && sig_loc > 0) {
      if (norm_dt && is.finite(dt_med) && dt_med > 0) sig_loc^2 / dt_med  # nm²/ms
      else sig_loc^2                                                       # nm²/step
    } else NA_real_
    noise_state <- if (is.finite(D_noise)) D_vals <= D_noise else rep(FALSE, K)

    list(K = K, sigma = sigma, A = A, pi = pi_v,
         D_vals = D_vals, fracs = fracs,
         norm_dt = norm_dt,
         seg_list = seg_list, n_tracks = length(obs_list),
         dt_med_ms = dt_med,
         sigma_loc_nm = sig_loc, D_noise = D_noise, noise_state = noise_state)
  })

  output$trk_hmm_params_plot <- renderPlotly({
    bd   <- hmm_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"
    state_cols <- c("#C03030","#D08020","#2080C0","#7030A0")[seq_len(bd$K)]

    d_unit  <- if (isTRUE(bd$norm_dt)) "nm\u00B2/ms" else "nm\u00B2/step"
    sig_unit <- if (isTRUE(bd$norm_dt)) "nm/\u221Ams" else "nm"
    is_noise <- if (!is.null(bd$noise_state)) bd$noise_state else rep(FALSE, bd$K)
    labels  <- sprintf("State %d\n%.0f%%", seq_len(bd$K), bd$fracs * 100)
    labels  <- ifelse(is_noise, paste0(labels, "\n\u26A0 noise-limited"), labels)
    bar_cols <- ifelse(is_noise, "#9e9e9e", state_cols)   # grey out noise/immobile-limited states
    plot_ly() %>%
      add_bars(
        x = bd$D_vals, y = labels,
        orientation = "h",
        marker = list(color = bar_cols),
        text  = sprintf("\u03C3=%.1f %s  D=%.2g %s%s", bd$sigma, sig_unit, bd$D_vals, d_unit,
                        ifelse(is_noise,
                               sprintf("  (\u2264 noise floor %.2g %s \u2014 immobile/noise-limited)", bd$D_noise, d_unit),
                               "")),
        hovertemplate = "%{text}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(title = paste("D (", d_unit, ")"), color = fg_col),
        yaxis = list(title = "", autorange = "reversed", color = fg_col),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font = list(color = fg_col), showlegend = FALSE,
        uirevision = "hmm_params"
      )
  })

  output$trk_hmm_trans_plot <- renderPlotly({
    bd   <- hmm_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"
    K      <- bd$K
    labels <- paste0("S", seq_len(K))
    plot_ly(
      x = labels, y = rev(labels),
      z = bd$A[rev(seq_len(K)), ],
      type = "heatmap",
      colorscale = list(list(0, "#08306b"), list(0.5, "#4393c3"), list(1, "#ffffff")),
      reversescale = FALSE, zmin = 0, zmax = 1,
      text  = matrix(sprintf("%.3f", bd$A[rev(seq_len(K)), ]), K, K),
      texttemplate = "%{text}",
      hovertemplate = "From %{y} \u2192 To %{x}: %{text}<extra></extra>"
    ) %>% layout(
      xaxis = list(title = "To state",   color = fg_col, side = "top"),
      yaxis = list(title = "From state", color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col), uirevision = "hmm_trans"
    )
  })

  output$trk_hmm_example_plot <- renderPlotly({
    bd     <- hmm_bundle()
    n_ex   <- max(1L, as.integer(isolate(input$trk_hmm_n_ex)))
    dark   <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"
    state_cols <- c("#C03030","#D08020","#2080C0","#7030A0")[seq_len(bd$K)]

    ex <- head(bd$seg_list, n_ex)
    # Lay out tracks in a grid — compute a consistent offset spacing from
    # the 90th percentile of Rmax across example tracks so the grid never
    # overlaps but also isn't wastefully spread out.
    spans <- sapply(ex, function(seg) {
      cx <- mean(seg$x); cy <- mean(seg$y)
      r  <- sqrt((seg$x - cx)^2 + (seg$y - cy)^2)
      max(r, na.rm = TRUE)
    })
    spacing <- max(50, quantile(spans, 0.90, na.rm = TRUE) * 2.5)
    n_cols  <- max(1L, ceiling(sqrt(length(ex))))

    fig <- plot_ly()
    is_first <- TRUE
    for (i in seq_along(ex)) {
      seg   <- ex[[i]]
      col_i <- ((i - 1L) %% n_cols)
      row_i <- ((i - 1L) %/% n_cols)
      cx <- mean(seg$x); cy <- mean(seg$y)
      xc <- seg$x - cx + col_i * spacing
      yc <- seg$y - cy - row_i * spacing   # negative so row 0 is at top
      # Grey track line
      fig <- fig %>% add_lines(
        x = xc, y = yc,
        line = list(color = "rgba(150,150,150,0.3)", width = 1),
        hoverinfo = "skip", showlegend = FALSE
      )
      # Points coloured by state
      for (k in seq_len(bd$K)) {
        idx <- seg$state == k
        if (!any(idx)) next
        fig <- fig %>% add_markers(
          x = xc[idx], y = yc[idx],
          marker = list(size = 4, color = state_cols[k], opacity = 0.8),
          name = paste0("State ", k),
          legendgroup = paste0("state", k),
          showlegend = is_first,
          hovertemplate = paste0("TID: ", seg$tid[1L], "  State ", k, "<extra></extra>")
        )
      }
      is_first <- FALSE
    }
    fig %>% layout(
      xaxis = list(title = "position (nm, centred + offset)", color = fg_col,
                   zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", color = fg_col, zeroline = FALSE,
                   showticklabels = FALSE, scaleanchor = "x", scaleratio = 1),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col),
      uirevision = "hmm_ex"
    )
  })

  output$trk_hmm_summary <- renderUI({
    bd       <- hmm_bundle()
    K        <- bd$K
    norm_dt  <- isTRUE(bd$norm_dt)
    dt       <- bd$dt_med_ms
    d_unit   <- if (norm_dt) "nm\u00B2/ms"
                else if (is.finite(dt) && dt > 0)
                  sprintf("nm\u00B2/step  (\u00F7%.2g ms \u2192 nm\u00B2/ms)", dt)
                else "nm\u00B2/step"
    sig_unit <- if (norm_dt) "nm/\u221Ams" else "nm"

    is_noise <- if (!is.null(bd$noise_state)) bd$noise_state else rep(FALSE, K)
    rows <- lapply(seq_len(K), function(k)
      tags$tr(
        tags$td(sprintf("State %d", k),                                     style = "padding:2px 8px;"),
        tags$td(sprintf("%.1f%%",  bd$fracs[k] * 100),                     style = "padding:2px 8px;"),
        tags$td(sprintf("\u03C3 = %.1f %s", bd$sigma[k], sig_unit),        style = "padding:2px 8px;"),
        tags$td(sprintf("D = %.3g %s", bd$D_vals[k], d_unit),              style = "padding:2px 8px;"),
        tags$td(sprintf("self-trans = %.3f", bd$A[k,k]),                   style = "padding:2px 8px;"),
        tags$td(if (isTRUE(is_noise[k]))
                  tags$span("\u26A0 noise/immobile-limited", style = "color:#b26a00;")
                else "",                                                    style = "padding:2px 8px;")
      )
    )
    tags$div(style = "font-size:0.88em; margin:4px 0;",
      tags$p(sprintf("n = %d tracks  |  median \u0394t = %.2g ms  |  \u221A\u0394t normalisation: %s",
                     bd$n_tracks, bd$dt_med_ms, if (norm_dt) "on" else "off"),
             style = "color:#888; margin:2px 0;"),
      tags$table(
        style = "border-collapse:collapse;",
        tags$thead(tags$tr(lapply(
          c("", "Fraction", "Width (\u03C3)", "Diffusivity (D)", "Self-transition", ""),
          function(h) tags$th(h, style = "padding:2px 8px; border-bottom:1px solid #ccc;")
        ))),
        tags$tbody(rows)
      ),
      if (is.finite(bd$D_noise %||% NA_real_))
        tags$p(style = "color:#888; font-size:0.92em; margin-top:4px;",
          sprintf("Noise floor D = %.2g %s (\u03C3_loc = %.2f nm). States at or below it are marked \u26A0 \u2014 they are indistinguishable from localization noise, not resolved diffusive states.",
                  bd$D_noise, d_unit, bd$sigma_loc_nm))
      else
        tags$p(style = "color:#888; font-size:0.92em; margin-top:4px;",
          "Set \u03C3_loc (run the unified MSD, or enter it in Acquisition Parameters) to flag noise/immobile-limited states.")
      ,
      # MINFLUX caveat: emissions are \u221A\u0394t-normalised (state \u03C3/D handle variable \u0394t),
      # but the transition matrix A is DISCRETE-TIME (per-step), so self-transition /
      # switching probabilities assume a roughly uniform \u0394t. Under heavy dark-gapping a
      # step spanning several skipped cycles is still counted as one transition
      # opportunity, so transition RATES are approximate. Read state \u03C3/D as robust,
      # transition probabilities as indicative. (See docs/MINFLUX_analysis_methods.md.)
      tags$p(style = "color:#888; font-size:0.85em; margin-top:6px; border-top:1px dashed #ddd; padding-top:4px;",
        "\u2139 Transition probabilities are per-step (discrete-time) and assume ~uniform \u0394t; ",
        "treat them as indicative under heavy dark-gapping. State \u03C3/D use \u221A\u0394t-normalised steps and are robust.")
    )
  })

  output$trk_dl_hmm_csv <- downloadHandler(
    filename = function() paste0("hmm_states_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      bd  <- hmm_bundle()
      out <- do.call(rbind, bd$seg_list)
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ---- Per-segment MSD by diffusion state (uses HMM Viterbi segments) -------
  # Splits each track into contiguous same-state runs, pools each state's segment
  # displacement-pairs into a per-state eMSD (same cumulative cycle-index binning
  # as the Unified MSD), and fits it with the Unified model + window settings.
  # This is the payoff for mixed-phase data: the confined model on the slow/wall
  # state gives a real confinement size L; the directed model on the fast/exported
  # state gives a velocity v — instead of one washed-out whole-track MSD.
  SEG_MSD_PAL <- c("#3060C0","#C03030","#2CA02C","#FF7F0E","#9467BD","#8C564B")

  seg_msd_rv <- eventReactive(input$trk_seg_msd_run, {
    req(isTRUE(input$trk_show_seg_msd))
    bd <- tryCatch(hmm_bundle(), error = function(e) NULL)
    shiny::validate(shiny::need(
      !is.null(bd) && !is.null(bd$seg_list) && length(bd$seg_list) > 0,
      "Enable the Diffusion State HMM (above) first — its Viterbi states define the segments."))
    K   <- bd$K
    t_c <- trk_cycle_time(); req(is.finite(t_c) && t_c > 0)
    cc  <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    sw  <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0
    min_seg <- max(5L, as.integer(input$trk_seg_msd_min_pts %||% 15L))
    has_z <- any(vapply(bd$seg_list, function(d) "z" %in% names(d) && any(is.finite(d$z)), logical(1)))
    ndim_sel <- input$trk_msd_ndim %||% "auto"
    ndim <- if (identical(ndim_sel, "3")) 3L else if (identical(ndim_sel, "2")) 2L else if (has_z) 3L else 2L
    max_lag <- max(60L, min(5000L, as.integer(input$trk_uni_max_lag %||% 50L)))

    dr2 <- matrix(0,  max_lag, K); np <- matrix(0L, max_lag, K); dts <- matrix(0, max_lag, K)
    n_seg <- integer(K); n_pts <- integer(K)

    withProgress(message = "Per-segment MSD…", value = 0, {
      for (di in seq_along(bd$seg_list)) {
        incProgress(1 / length(bd$seg_list))
        d <- bd$seg_list[[di]]; d <- d[order(d$tim), ]
        st <- d$state; n <- nrow(d); if (n < 2L) next
        runs <- rle(st); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
        for (r in seq_along(runs$lengths)) {
          sstate <- runs$values[r]
          if (is.na(sstate) || sstate < 1L || sstate > K) next
          idx <- starts[r]:ends[r]; if (length(idx) < min_seg) next
          tb <- d[idx, , drop = FALSE]
          n_seg[sstate] <- n_seg[sstate] + 1L; n_pts[sstate] <- n_pts[sstate] + length(idx)
          dcyc <- pmax(1L, as.integer(round((diff(tb$tim) - sw) / t_c)))
          cyc  <- c(0L, cumsum(dcyc)); nn <- nrow(tb)
          for (i in seq_len(nn - 1L)) {
            dt_v <- tb$tim[seq.int(i + 1L, nn)] - tb$tim[i]
            lag  <- cyc[seq.int(i + 1L, nn)] - cyc[i]
            okk  <- lag >= 1L & lag <= max_lag
            if (!any(okk)) next
            js <- seq.int(i + 1L, nn)[okk]; b <- lag[okk]; dtk <- dt_v[okk]
            r2 <- (tb$x[js] - tb$x[i])^2 + (tb$y[js] - tb$y[i])^2
            if (has_z && "z" %in% names(tb)) r2 <- r2 + (tb$z[js] - tb$z[i])^2
            rs <- rowsum(r2, b, reorder = FALSE); ix <- as.integer(rownames(rs))
            dr2[ix, sstate] <- dr2[ix, sstate] + rs[, 1L]
            rt <- rowsum(dtk, b, reorder = FALSE); it <- as.integer(rownames(rt))
            dts[it, sstate] <- dts[it, sstate] + rt[, 1L]
            np[, sstate] <- np[, sstate] + tabulate(b, nbins = max_lag)
          }
        }
      }
    })

    model    <- input$trk_uni_model %||% "linear"
    min_lag  <- max(1L, as.integer(input$trk_uni_min_lag %||% 5L))
    max_fit  <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    use_iter <- isTRUE(input$trk_uni_iterative)
    linear_tau <- isTRUE(input$trk_uni_linear_tau)
    lag_v <- seq_len(max_lag)
    states <- Filter(Negate(is.null), lapply(seq_len(K), function(k) {
      npk <- np[, k]; if (sum(npk) == 0L) return(NULL)
      msd <- ifelse(npk > 0L, dr2[, k] / npk, NA_real_)
      tau_use <- if (linear_tau) lag_v * t_c else ifelse(npk > 0L, dts[, k] / npk, lag_v * t_c)
      okf <- is.finite(tau_use) & is.finite(msd) & lag_v >= min_lag & lag_v <= max_fit & npk > 0L
      fr  <- if (sum(okf) >= 3L)
               tryCatch(.msd_model_fit(tau_use[okf], msd[okf], npk[okf], ndim, model, use_iter, t_c),
                        error = function(e) NULL) else NULL
      list(state = k, n_seg = n_seg[k], n_pts = n_pts[k], D_hmm = bd$D_vals[k] / 1e3,
           lag = lag_v, tau = tau_use, msd = msd, npairs = npk, fr = fr)
    }))
    shiny::validate(shiny::need(length(states) > 0, "No segments long enough — lower 'Min points per segment' or the HMM min points."))
    list(states = states, K = K, ndim = ndim, t_cycle = t_c, model = model, linear_tau = linear_tau)
  })

  output$trk_seg_msd_plot <- renderPlotly({
    req(isTRUE(input$trk_show_seg_msd))
    sm <- seg_msd_rv()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    fig <- plot_ly()
    for (st in sm$states) {
      ok <- is.finite(st$tau) & is.finite(st$msd) & st$npairs > 0L
      if (sum(ok) < 2L) next
      col <- SEG_MSD_PAL[((st$state - 1L) %% length(SEG_MSD_PAL)) + 1L]
      fig <- fig %>% add_trace(x = st$tau[ok], y = st$msd[ok] / 1e6, type = "scatter",
        mode = "lines+markers", marker = list(size = 4, color = col),
        line = list(color = col, width = 2), name = sprintf("State %d eMSD (%d seg)", st$state, st$n_seg))
      fr <- st$fr
      if (!is.null(fr)) {
        lbl <- switch(fr$model %||% "linear",
          confined  = if (isTRUE(fr$degenerate)) sprintf("State %d: ≈linear (no plateau)", st$state)
                      else sprintf("State %d confined: D=%.3g, L=%.0f nm", st$state, fr$D / 1e3, fr$L_conf),
          directed  = sprintf("State %d directed: D=%.3g, v=%.3g µm/s", st$state, fr$D / 1e3, fr$v_um_s),
          anomalous = sprintf("State %d: D=%.3g, α=%.2f", st$state, fr$D / 1e3, fr$alpha),
          sprintf("State %d linear: D=%.3g µm²/s", st$state, fr$D / 1e3))
        fig <- fig %>% add_trace(x = fr$tau_fit, y = pmax(0, fr$msd_fit) / 1e6, type = "scatter",
          mode = "lines", line = list(color = col, width = 2.5, dash = "dash"), name = lbl)
      }
    }
    fig %>% layout(
      xaxis = list(title = if (sm$linear_tau) "τ (ms, lag × cycle)" else "τ (ms, mean Δt)", color = fg),
      yaxis = list(title = "MSD (µm²)", color = fg, rangemode = "tozero"),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(font = list(color = fg, size = 10)),
      margin = list(l = 55, r = 10, t = 10, b = 45), uirevision = "segmsd")
  })

  output$trk_seg_msd_summary <- renderUI({
    req(isTRUE(input$trk_show_seg_msd))
    sm <- tryCatch(seg_msd_rv(), error = function(e) NULL); if (is.null(sm)) return(NULL)
    rows <- lapply(sm$states, function(st) {
      fr <- st$fr
      param <- if (is.null(fr)) "—" else switch(fr$model %||% "linear",
        confined  = if (isTRUE(fr$degenerate)) "no plateau (≈linear)" else sprintf("L = %.0f nm", fr$L_conf),
        directed  = sprintf("v = %.3g µm/s", fr$v_um_s),
        anomalous = sprintf("α = %.2f", fr$alpha), "—")
      D_msd <- if (!is.null(fr) && is.finite(fr$D)) sprintf("%.3g", fr$D / 1e3) else "—"
      col <- SEG_MSD_PAL[((st$state - 1L) %% length(SEG_MSD_PAL)) + 1L]
      htmltools::tags$tr(
        htmltools::tags$td(htmltools::tags$span(style = sprintf("color:%s;font-weight:bold;", col), sprintf("State %d", st$state))),
        htmltools::tags$td(st$n_seg), htmltools::tags$td(format(st$n_pts, big.mark = ",")),
        htmltools::tags$td(D_msd), htmltools::tags$td(param),
        htmltools::tags$td(sprintf("%.3g", st$D_hmm)))
    })
    htmltools::div(style = "font-size:0.85em; margin-top:6px;",
      htmltools::tags$b(sprintf("Per-state MSD (%s model):", sm$model)),
      htmltools::tags$table(class = "table table-condensed", style = "width:auto; margin-top:4px;",
        htmltools::tags$thead(htmltools::tags$tr(lapply(
          c("State", "segments", "points", "D (MSD, µm²/s)", "model param", "D (HMM, µm²/s)"),
          htmltools::tags$th))),
        htmltools::tags$tbody(rows)),
      htmltools::tags$div(style = "color:#888; font-size:0.92em;",
        "States are the HMM diffusion states (sorted slow→fast). Switch the Unified MSD ",
        "'Fit model' to Confined (for a slow/wall state → size L) or Directed (fast/exported → velocity v)."))
  })

  output$trk_dl_seg_msd_csv <- downloadHandler(
    filename = function() paste0("segment_msd_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      sm <- seg_msd_rv()
      out <- do.call(rbind, lapply(sm$states, function(st) {
        ok <- st$npairs > 0L
        data.frame(state = st$state, cycle_lag = st$lag[ok], tau_ms = st$tau[ok],
                   MSD_um2 = st$msd[ok] / 1e6, n_pairs = st$npairs[ok])
      }))
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ---- Spatially resolved D-map --------------------------------------------
  # Uses a single-lag estimator: D_local = mean(r²) / (4 * median_dt)
  # This is far more robust than fitting a full MSD curve per cell, because
  # in small cells (e.g. E. coli) there are rarely enough steps per tau-bin
  # to constrain a multi-parameter fit. The single-lag estimator is the short-lag
  # limit of the covariance-based estimator (Vestergaard, Blainey & Flyvbjerg 2014,
  # Phys. Rev. E 89:022726); for MINFLUX-specific diffusion-parameter validity see
  # Vogler et al. 2025 (Commun. Biol. 8:1573). [Prev. mis-cited to Ostersehlt 2022,
  # which is a DNA-PAINT imaging paper with no diffusion mapping — corrected.]
  dmap_bundle <- eventReactive(input$trk_dmap_run, {
    req(isTRUE(input$trk_show_dmap))
    df <- trk(); req(nrow(df) > 0)

    n_bins  <- max(5L, as.integer(input$trk_dmap_bins))
    min_pts <- max(3L, as.integer(input$trk_dmap_min_pts))

    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
    df <- df[order(df$tid, df$tim), ]
    sp <- split(df, df$tid)

    withProgress(message = "Computing D-map…", value = 0.1, {
    # Build step table: midpoint position, dt, r²
    step_rows <- lapply(sp, function(tb) {
      if (nrow(tb) < 2L) return(NULL)
      n  <- nrow(tb)
      i1 <- seq_len(n - 1L); i2 <- i1 + 1L
      dt <- tb$tim[i2] - tb$tim[i1]
      data.frame(
        mx = 0.5 * (tb$x[i1] + tb$x[i2]),
        my = 0.5 * (tb$y[i1] + tb$y[i2]),
        dt = dt,
        r2 = (tb$x[i2] - tb$x[i1])^2 + (tb$y[i2] - tb$y[i1])^2,
        stringsAsFactors = FALSE
      )
    })
    steps <- do.call(rbind, step_rows[!sapply(step_rows, is.null)])
    steps <- steps[is.finite(steps$dt) & steps$dt > 0 & is.finite(steps$r2), ]
    shiny::validate(shiny::need(nrow(steps) >= min_pts * 2,
      "Too few steps \u2014 relax filters or lower 'Min steps per cell'."))

    # Optional: keep only INTEGRATED (single-cycle) steps. A step spanning dark
    # cycles has a large \u0394t and a spatially meaningless midpoint (the molecule may
    # have moved far while unseen), so mixing them in smears the map. Classified by
    # cycle count with the same cycle time the TSA uses.
    t_c_dmap <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    integrated_only <- isTRUE(input$trk_dmap_integrated_only) &&
                       is.finite(t_c_dmap) && t_c_dmap > 0
    if (integrated_only) {
      keep_step <- pmax(1L, as.integer(round(steps$dt / t_c_dmap))) == 1L
      steps <- steps[keep_step, , drop = FALSE]
      shiny::validate(shiny::need(nrow(steps) >= min_pts * 2,
        "Too few integrated (single-cycle) steps \u2014 untick 'Integrated steps only' or relax filters."))
    }

    # Global median dt used to normalise all cells to the same lag
    dt_med <- median(steps$dt, na.rm = TRUE)
    shiny::validate(shiny::need(is.finite(dt_med) && dt_med > 0,
      "Cannot determine median \u0394t."))

    # Per-step apparent 2D diffusion: D_raw = r² / (4·Δt).
    steps$D_raw <- steps$r2 / (4 * steps$dt)
    # Optional: subtract the static localization-error floor. A single 2D step has
    # E[r²_noise] = 4σ² → an apparent D of σ²/Δt; without this, immobile/confined
    # regions show a spurious diffusion floor. Clamp at 0. σ_loc is the app-wide
    # value (unified MSD fit or manual override); no-op when it is unset.
    sig_loc <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    noise_subtracted <- isTRUE(input$trk_dmap_sub_noise) &&
                        is.finite(sig_loc) && sig_loc > 0
    if (noise_subtracted)
      steps$D_raw <- pmax(0, steps$D_raw - sig_loc^2 / steps$dt)

    x_rng    <- range(steps$mx); y_rng <- range(steps$my)
    x_breaks <- seq(x_rng[1], x_rng[2], length.out = n_bins + 1L)
    y_breaks <- seq(y_rng[1], y_rng[2], length.out = n_bins + 1L)
    cell_w   <- diff(x_rng) / n_bins
    cell_h   <- diff(y_rng) / n_bins

    steps$ix <- pmin(findInterval(steps$mx, x_breaks), n_bins)
    steps$iy <- pmin(findInterval(steps$my, y_breaks), n_bins)

    # Per-cell D: median of per-step D_raw (median is robust to outlier jumps)
    D_mat    <- matrix(NA_real_, n_bins, n_bins)
    N_mat    <- matrix(0L,       n_bins, n_bins)

    cell_idx <- steps$ix + (steps$iy - 1L) * n_bins
    by_cell  <- split(steps$D_raw, cell_idx)

    for (nm in names(by_cell)) {
      vals <- by_cell[[nm]][is.finite(by_cell[[nm]])]
      if (length(vals) < min_pts) next
      idx  <- as.integer(nm)
      ix   <- ((idx - 1L) %% n_bins) + 1L
      iy   <- ((idx - 1L) %/% n_bins) + 1L
      D_mat[ix, iy] <- median(vals)
      N_mat[ix, iy] <- length(vals)
    }

    x_centers <- 0.5 * (x_breaks[-1] + x_breaks[-length(x_breaks)])
    y_centers <- 0.5 * (y_breaks[-1] + y_breaks[-length(y_breaks)])

    list(D_mat = D_mat, N_mat = N_mat,
         x_centers = x_centers, y_centers = y_centers,
         cell_w_nm = cell_w, cell_h_nm = cell_h,
         n_bins = n_bins, n_steps = nrow(steps),
         dt_med_ms = dt_med,
         integrated_only = integrated_only,
         noise_subtracted = noise_subtracted,
         sigma_loc_nm = if (noise_subtracted) sig_loc else NA_real_)
    })  # withProgress
  })

  output$trk_dmap_plot <- renderPlotly({
    bd     <- dmap_bundle()
    dark   <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"

    # D_mat[ix, iy]: ix = x-column, iy = y-row
    # Plotly heatmap z[row, col] => transpose so y is rows, x is cols
    z     <- t(bd$D_mat)
    z_log <- log10(z); z_log[!is.finite(z_log)] <- NA_real_

    finite_vals <- z_log[is.finite(z_log)]
    shiny::validate(shiny::need(length(finite_vals) > 0,
      "No cells met the minimum step threshold \u2014 lower 'Min steps per cell' or reduce grid resolution."))

    # Clip to 2ndÃÂ¢Ã¢ÂÂ¬Ã¢ÂÂ98th percentile so one extreme cell doesn't dominate the scale
    zmin <- floor(quantile(finite_vals, 0.02, na.rm = TRUE))
    zmax <- ceiling(quantile(finite_vals, 0.98, na.rm = TRUE))
    if (zmin >= zmax) { zmin <- zmin - 1; zmax <- zmax + 1 }

    cs <- list(
      list(0,    "#08306b"),
      list(0.25, "#2171b5"),
      list(0.5,  "#f7f7f7"),
      list(0.75, "#f4a582"),
      list(1,    "#b2182b")
    )

    # Build hover text matrix (D in nm²/ms + step count)
    n_row <- nrow(z_log); n_col <- ncol(z_log)
    N_t   <- t(bd$N_mat)
    hover <- matrix("", n_row, n_col)
    for (r in seq_len(n_row))
      for (c in seq_len(n_col)) {
        d <- bd$D_mat[c, r]
        n <- bd$N_mat[c, r]
        hover[r, c] <- if (is.finite(d))
          sprintf("x: %.0f nm<br>y: %.0f nm<br>D: %.3g nm\u00B2/ms<br>n steps: %d",
                  bd$x_centers[c], bd$y_centers[r], d, n)
        else sprintf("x: %.0f nm<br>y: %.0f nm<br>(no data)", bd$x_centers[c], bd$y_centers[r])
      }

    fig <- plot_ly(
      x = bd$x_centers, y = bd$y_centers, z = z_log,
      type = "heatmap",
      colorscale = cs, zmin = zmin, zmax = zmax,
      text = hover, hovertemplate = "%{text}<extra></extra>",
      colorbar = list(
        title = list(text = "log\u2081\u2080 D<br>(nm\u00B2/ms)"),
        tickvals = seq(zmin, zmax, by = 1),
        ticktext = paste0("10^", seq(zmin, zmax, by = 1)),
        len = 0.6
      )
    )

    # Optional: overlay step count as text annotations
    if (isTRUE(isolate(input$trk_dmap_show_n))) {
      ann <- list()
      for (r in seq_len(n_row))
        for (c in seq_len(n_col))
          if (bd$N_mat[c, r] > 0)
            ann <- c(ann, list(list(
              x = bd$x_centers[c], y = bd$y_centers[r],
              text = as.character(bd$N_mat[c, r]),
              showarrow = FALSE, xref = "x", yref = "y",
              font = list(size = 8, color = "rgba(0,0,0,0.6)")
            )))
    } else { ann <- list() }

    fig %>% layout(
      xaxis = list(title = "x (nm)", color = fg_col,
                   scaleanchor = "y", scaleratio = 1),
      yaxis = list(title = "y (nm)", color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col),
      annotations = ann,
      uirevision = "dmap"
    )
  })

  output$trk_dmap_summary <- renderUI({
    bd <- dmap_bundle()
    D_vals <- bd$D_mat[is.finite(bd$D_mat) & bd$D_mat > 0]
    if (!length(D_vals)) return(NULL)
    corr_lbl <- paste0(
      if (isTRUE(bd$integrated_only)) "integrated steps only" else "all steps",
      "  |  ",
      if (isTRUE(bd$noise_subtracted))
        sprintf("noise floor subtracted (\u03C3_loc %.2f nm)", bd$sigma_loc_nm)
      else "no noise-floor subtraction")
    tags$div(style = "font-size:0.88em; margin:4px 0;",
      sprintf(
        "Filled cells: %d / %d  |  D: %.3g \u2013 %.3g nm\u00B2/ms  |  median D: %.3g nm\u00B2/ms  |  cell size: %.0f \u00D7 %.0f nm  |  median \u0394t: %.2g ms  |  total steps: %d",
        length(D_vals), bd$n_bins^2,
        min(D_vals), max(D_vals), median(D_vals),
        bd$cell_w_nm, bd$cell_h_nm,
        bd$dt_med_ms, bd$n_steps
      ),
      tags$br(),
      tags$span(style = "color:#888;", corr_lbl)
    )
  })

  output$trk_dl_dmap_csv <- downloadHandler(
    filename = function() paste0("dmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      bd  <- dmap_bundle()
      rows <- lapply(seq_len(bd$n_bins), function(ix)
        lapply(seq_len(bd$n_bins), function(iy)
          data.frame(x_nm = bd$x_centers[ix], y_nm = bd$y_centers[iy],
                     D_nm2_per_ms = bd$D_mat[ix, iy],
                     n_steps = bd$N_mat[ix, iy],
                     stringsAsFactors = FALSE)
        )
      )
      out <- do.call(rbind, do.call(c, rows))
      out <- out[is.finite(out$D_nm2_per_ms), ]
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ---- Playback engine ----
  # --- Play / Pause ---
  observeEvent(input$trk_btn_play, {
    if (isTRUE(is_playing())) return(invisible(NULL))
    if (is_seq()) {
      # validate data BEFORE mutating loop state  -  req() would leave state inconsistent
      sd  <- isolate(seq_data())
      pre <- isolate(frames_store())
      tm  <- isolate(trace_map())
      if (!is.list(sd) || !length(sd$tids)) return(invisible(NULL))
      if (!is.list(pre)) return(invisible(NULL))
      if (is.null(tm))   return(invisible(NULL))
      # start/restart sequential
      loop$gen       <- loop$gen + 1L
      loop$playing   <- TRUE
      
      loop$sp            <- sd$sp
      loop$tids          <- sd$tids
      loop$trail_len     <- pre$trail_len
      loop$trail_nmax    <- pre$trail_nmax
      loop$pt_size       <- input$trk_pt_size
      loop$col_map       <- pre$col_map
      loop$markers_idx   <- tm$markers_idx
      loop$trail_start   <- tm$trail_start
      loop$frames_per_tid <- max(2L, n_frames_eff())
      loop$speed         <- as.numeric(input$trk_play_speed)
      loop$interval_sec  <- play_interval_for(loop$speed)
      loop$ms_per_unit   <- as.numeric(input$trk_ms_per_unit)
      loop$times <- lapply(loop$tids, function(nm) {
        tb <- loop$sp[[nm]]; r <- range(tb$t, na.rm = TRUE)
        if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) r <- c(0, 1)
        seq(r[1], r[2], length.out = loop$frames_per_tid)
      })
      loop$cur_tid   <- 1L
      loop$cur_frame <- 1L
      is_playing(TRUE)
      
      run_sequential_loop_plain(loop$gen)
      
    } else {
      # start/restart global
      loop$ggen <- loop$ggen + 1L
      loop$gplaying <- TRUE 
      is_playing(TRUE)
      start_global_loop(loop$ggen)
    }
  }, ignoreInit = TRUE)
  start_global_loop <- function(id) {
    pre <- isolate(frames_store())
    if (!is.list(pre) || !length(pre$frames)) return(invisible(NULL))
    # capture maxk once; the per-frame delay is read live from loop$speed so the
    # "Playback speed (×)" slider takes effect mid-playback (see loop_once below).
    maxk     <- length(pre$frames)
    step     <- 1L
    loop$speed    <- as.numeric(isolate(input$trk_play_speed))
    loop$gplaying <- TRUE
    draw_frame(1L, force = TRUE)
    
    k_local <- 0L   # private counter, don't read tick()
    
    loop_once <- function() {
      # Wrap entire body so an unexpected error never silently kills the loop
      tryCatch({
        if (!isTRUE(loop$gplaying) || id != loop$ggen) return(invisible(NULL))
        
        k_local <<- k_local + step
        if (k_local < 1L) k_local <<- 1L
        if (k_local > maxk) {
          # Playback finished naturally  -  reset state inside domain so reactive
          # writes flush correctly from a later() callback context
          loop$gplaying <- FALSE
          loop$ggen <- loop$ggen + 1L
          shiny::withReactiveDomain(session$domain, {
            isolate({ is_playing(FALSE); slider_guard(0L) })
          })
          return(invisible(NULL))
        }

        # run Shiny ops inside domain, and reactive reads inside isolate()
        shiny::withReactiveDomain(session$domain, {
          isolate({
            tick(k_local)
            draw_frame(k_local)
          })
        })
        
        if (isTRUE(loop$gplaying) && id == loop$ggen) {
          later::later(loop_once, play_interval_for(loop$speed))
        }
      }, error = function(e) {
        # On error: log to console and keep the loop alive if still playing
        message("draw_frame error (frame ", k_local, "): ", conditionMessage(e))
        if (isTRUE(loop$gplaying) && id == loop$ggen) {
          later::later(loop_once, play_interval_for(loop$speed))
        }
      })
    }
    
    # kick off immediately
    later::later(loop_once, 0)
  }
  observeEvent(input$trk_chunk_sel, {
    # Skip the creation fire after upload/rechunk (value == the just-rendered chunk);
    # data_key already rebuilt. Only a genuine chunk switch gets through.
    if (identical(input$trk_chunk_sel %||% "__all__", loop$rendered_chunk %||% "__all__")) return()
    if (isTruthy(input$trk_chunk_sel)) fire_layout("chunk_sel")
  }, ignoreInit = TRUE)
  # Fires whenever new data is loaded (upload or rechunk), regardless of whether
  # trk_chunk_sel changed. Ensures the chart always rebuilds with the new dataset.
  # This is the SINGLE authoritative first-plot trigger on load; the chunk_sel /
  # view_mode / z_stretch observers compare against loop$rendered_* (set by the render)
  # so their input-CREATION after load does not fire extra redundant rebuilds.
  observeEvent(trk_data_key(), {
    if (!identical(trk_data_key(), "init")) fire_layout("data_key")
  }, ignoreInit = TRUE)
  observeEvent(input$trk_btn_pause, {
    loop$gplaying <- FALSE   # stop global loop (belt-and-suspenders with ggen check)
    loop$playing  <- FALSE   # stop sequential loop
    loop$ggen <- loop$ggen + 1L
    loop$gen  <- loop$gen  + 1L
    is_playing(FALSE)
    slider_guard(0L)
  }, ignoreInit = TRUE)
  
  observeEvent(input$trk_btn_stop, {
    is_playing(FALSE)
    loop$playing  <- FALSE
    loop$gplaying <- FALSE
    loop$ggen     <- loop$ggen + 1L
    loop$gen      <- loop$gen + 1L
    slider_guard(0L)
    
    loop$cur_tid   <- 1L
    loop$cur_frame <- 1L
    
    if (isTRUE(isolate(is_seq()))) {
      sd <- isolate(seq_data())
      if (!is.null(sd) && length(sd$tids) > 0) draw_tid_point(sd$tids[1], 1L)
    } else {
      pre <- isolate(frames_store())
      if (!is.null(pre) && length(pre$frames) > 0) {
        tick(1L)
        draw_frame(1L, force = TRUE)   # <-- force a redraw while "not playing"
      }
    }
  })
  apply_auto_tune <- function() {
    df <- isolate(trk()); req(nrow(df) > 0)
    mode <- if (isTRUE(isolate(is_simul()))) "simul" else "seq"
    s <- smart_defaults(df, mode = mode, target_frames = 200, target_fps = 15)

    # Set the plain (non-reactive) flag BEFORE queuing updateNumericInput calls.
    # Browser round-trips from those calls arrive ~100-400ms later; rebuild_tick
    # must stay suppressed until AFTER fire_frames() has consumed the new values.
    # fire_frames() is delayed 500 ms so the new input values are in place when
    # frames_store() re-evaluates.  auto_tuning resets 1.5 s after that.
    loop$auto_tuning <- TRUE
    later::later(function() { loop$auto_tuning <- FALSE }, delay = 1.5)

    updateNumericInput(session, "trk_play_total",      value = s$nice_total)
    updateNumericInput(session, "trk_trail_len",       value = signif(s$trail_len, 3))
    updateNumericInput(session, "trk_target_markers",  value = s$target_markers)
    updateNumericInput(session, "trk_ms_per_unit",     value = s$ms_per_unit)
    updateNumericInput(session, "trk_bin_dt",          value = s$nice_total / s$n_frames)

    # Delay fire_frames until after input round-trips arrive (~200-400 ms),
    # so frames_store() reads the newly-tuned bin_dt / n_frames / total values.
    later::later(function() { fire_frames() }, delay = 0.5)
  }
  
  observeEvent(input$trk_auto_tune, {
    apply_auto_tune()
  }, ignoreInit = TRUE)
  # --- Global (prebinned) mode stays on invalidateLater (fast) ---
  run_sequential_loop_plain <- function(id) {
    tryCatch({
      # cancel if paused or superseded
      if (id != loop$gen || !isTRUE(loop$playing)) return(invisible(NULL))
      
      tids <- loop$tids; sp <- loop$sp
      if (is.null(tids) || is.null(sp) || !length(tids)) return(invisible(NULL))
      
      ti <- loop$cur_tid
      fi <- loop$cur_frame
      if (ti < 1L || ti > length(tids)) ti <- 1L
      if (fi < 1L || fi > loop$frames_per_tid) fi <- 1L
      nm <- tids[ti]; tb <- sp[[nm]]
      
      if (is.null(tb) || !nrow(tb)) {
        if (ti >= length(tids)) {
          loop$playing <- FALSE
          loop$gen <- loop$gen + 1L
          shiny::withReactiveDomain(session$domain, {
            isolate({ is_playing(FALSE); slider_guard(0L) })
          })
          return(invisible(NULL))
        }
        loop$cur_tid <- ti + 1L
        loop$cur_frame <- 1L
        later::later(function() { run_sequential_loop_plain(id) }, 0.05)
        return(invisible(NULL))
      }
      
      # current frame time + where we're going next
      tnow    <- loop$times[[ti]][fi]
      next_fi <- if (fi < loop$frames_per_tid) fi + 1L else fi
      next_ti <- ti
      if (fi >= loop$frames_per_tid && ti < length(tids)) {
        next_ti <- ti + 1L
        next_fi <- 1L
      }
      
      # draw this frame
      w <- which(tb$t >= (tnow - loop$trail_len) & tb$t <= tnow)
      if (length(w) > loop$trail_nmax) w <- tail(w, loop$trail_nmax)
      xx <- if (length(w)) c(tb$x[w], NA_real_) else numeric(0)
      yy <- if (length(w)) c(tb$y[w], NA_real_) else numeric(0)
      
      idx <- nearest_idx(tb$t, tnow)
      mx  <- if (is.na(idx)) numeric(0) else tb$x[idx]
      my  <- if (is.na(idx)) numeric(0) else tb$y[idx]
      
      col   <- as.character(loop$col_map[[nm]])
      proxy <- plotlyProxy("trk_track2D", session)
      if (!is.na(loop$trail_start))
        plotlyProxyInvoke(proxy, "restyle", list(x = list(xx), y = list(yy), `line.color` = list(col)), loop$trail_start)
      plotlyProxyInvoke(proxy, "restyle", list(x = list(mx), y = list(my), `marker.size` = loop$pt_size, `marker.color` = list(col)), loop$markers_idx)
      # Skip relayout during playback to avoid scaleanchor jitter
      
      # advance state
      if (fi < loop$frames_per_tid) {
        loop$cur_frame <- fi + 1L
      } else if (ti < length(tids)) {
        loop$cur_tid   <- ti + 1L
        loop$cur_frame <- 1L
      } else {
        loop$playing <- FALSE
        loop$gen <- loop$gen + 1L
        shiny::withReactiveDomain(session$domain, {
          isolate({ is_playing(FALSE); slider_guard(0L) })
        })
        return(invisible(NULL))
      }
      
      # compute delay (ms_per_unit -> variable; else fixed play_interval)
      # The "Playback speed (×)" multiplier scales every branch: higher speed =>
      # shorter delay. Everything is floored at MIN_INTERVAL_SEC so the loop can
      # never outpace the websocket (the old 0.01s firehose is what made playback
      # "refuse to stop" — it starved the Stop click + reactive flush).
      speed     <- if (is.finite(loop$speed) && loop$speed > 0) loop$speed else 1
      floor_sec <- play_interval_for(speed)
      if (is.finite(loop$ms_per_unit) && loop$ms_per_unit > 0) {
        tnext    <- loop$times[[ti]][min(next_fi, loop$frames_per_tid)]
        dt_units <- if (next_ti != ti) 0 else (tnext - tnow)
        delay_sec <- max(floor_sec, abs(dt_units) * (loop$ms_per_unit / 1000) / speed)
        if (next_ti != ti) delay_sec <- max(delay_sec, 0.15 / speed)
      } else {
        delay_sec <- floor_sec
      }
      
      if (id == loop$gen && isTRUE(loop$playing)) {
        later::later(function() { run_sequential_loop_plain(id) }, delay_sec)
      }
    }, error = function(e) {
      message("Sequential loop error (tid=", loop$cur_tid, " frame=", loop$cur_frame, "): ", conditionMessage(e))
      if (id == loop$gen && isTRUE(loop$playing)) {
        later::later(function() { run_sequential_loop_plain(id) }, loop$interval_sec)
      }
    })
  }
  
  # Live speed updates: both loops read loop$speed on each scheduled tick, so
  # dragging the slider during playback changes the rate without restarting.
  observeEvent(input$trk_play_speed, {
    loop$speed <- as.numeric(input$trk_play_speed)
    if (isTRUE(loop$playing)) loop$interval_sec <- play_interval_for(loop$speed)
  }, ignoreInit = TRUE)

  observeEvent(input$trk_play_mode, {
    is_playing(FALSE); loop$playing <- FALSE; loop$gplaying <- FALSE
    loop$ggen <- loop$ggen + 1L
    loop$gen    <- loop$gen + 1L
    slider_guard(0L)
  })
  
  # ---- Slider (seek without auto-pause) ----
  output$trk_frame_slider_ui <- renderUI({
    req(!is.na(layout_trigger()))
    if (is_seq()) {
      sd <- seq_data(); req(is.list(sd), length(sd$tids) > 0)
      ti <- current_tid_idx(); if (ti < 1L || ti > length(sd$tids)) ti <- 1L
      nm <- sd$tids[ti]
      sliderInput(
        "trk_frame_idx", paste0("Step (tid = ", nm, ")"),
        min = 1, max = max(2L, n_frames_eff()), value = 1, step = 1, width = "100%"
      )
    } else {
      pre <- frames_store(); req(length(pre$frames) > 0)
      sliderInput("trk_frame_idx", "Frame",
                  min = 1, max = length(pre$frames), value = 1, step = 1, width = "100%")
    }
  })
  
  
  observeEvent(input$trk_frame_idx, {
    req(!is.na(layout_trigger()))
    if (isTRUE(is_playing())) return(invisible(NULL))
    guard <- isolate(slider_guard())
    if (guard > 0L) {
      slider_guard(guard - 1L)
      return(invisible(NULL))
    }
    
    if (is_seq()) {
      # unchanged: sequential seeks by time within the current TID
      sd <- seq_data(); req(is.list(sd), length(sd$tids) > 0)
      ti <- current_tid_idx(); if (ti < 1L || ti > length(sd$tids)) ti <- 1L
      nm <- sd$tids[ti]; tb <- sd$sp[[nm]]; req(nrow(tb) > 0)
      
      k <- max(1L, min(as.integer(input$trk_frame_idx), max(2L, n_frames_eff())))
      tnow <- if (!is.null(loop$times) && length(loop$times) >= ti) {
        loop$times[[ti]][k]
      } else {
        r <- range(tb$t, na.rm = TRUE)
        seq(r[1], r[2], length.out = max(2L, n_frames_eff()))[k]
      }
      i <- nearest_idx(tb$t, tnow)
      draw_tid_point(nm, i)
      
    } else {
      # simultaneous seeks into the prebuilt frames
      pre <- frames_store(); req(length(pre$frames) > 0)
      k <- max(1L, min(as.integer(input$trk_frame_idx), length(pre$frames)))
      tick(k)                     # optional: keep any tick-based UI in sync
      draw_frame(k, force = TRUE) # <- force draw while paused
    }
  })

  # ---- Reset state + push frame 1 on data rebuild --------------------------------
  # For fire_layout: renderPlotly handles the full rebuild including frame 1 data.
  #   The plot_trigger observer just resets play state and absorbs the slider
  #   round-trip from trk_frame_slider_ui re-rendering.
  # For fire_frames (Plot Data, rebuild_tick, auto-tune): renderPlotly does NOT
  #   re-fire (frames_store is now isolated there). After frames_store re-evaluates,
  #   we push frame 1 data + new axis ranges via proxy so the chart updates in-place
  #   with no WebGL reinit and no blank.
  observeEvent(plot_trigger(), {
    is_playing(FALSE); tick(0L)
    loop$playing  <- FALSE
    loop$gplaying <- FALSE
    loop$ggen     <- loop$ggen + 1L
    loop$gen      <- loop$gen  + 1L
    was_layout <- isTRUE(loop$last_fire_layout)
    slider_guard(if (was_layout) 3L else 0L)

    if (!was_layout) {
      # fire_frames path: push new frame 1 + axes via proxy after frames_store settles.
      # Skip entirely in full-tracks mode — the static plot needs no proxy update.
      if (identical(isolate(input$trk_view_mode_main), "full")) return()
      # Capture generation so a stale callback from a prior session cannot jump the
      # display to frame 1 after a new upload/playback has already begun.
      snap_ggen <- loop$ggen
      snap_gen  <- loop$gen
      later::later(function() {
        shiny::withReactiveDomain(session$domain, {
          tryCatch({
            # Abort if a new upload or play session has superseded this callback.
            if (loop$ggen != snap_ggen || loop$gen != snap_gen) return()
            pre <- isolate(frames_store())
            if (!is.list(pre) || !length(pre$frames)) return()
            tm <- isolate(trace_map())
            if (is.null(tm)) return()

            # push frame 1 marker data (colors already baked into frames)
            draw_frame(1L, force = TRUE)

            # If per-TID trail traces exist, update their line colors too
            # (trails are set at renderPlotly time; conf coloring needs a proxy push)
            if (identical(tm$trails_mode, "per_tid") && tm$n_tid > 0) {
              line_cols <- vector("list", tm$n_tid)
              for (i in seq_len(tm$n_tid)) {
                tid_i <- if (i <= length(pre$tids)) pre$tids[i] else ""
                line_cols[[i]] <- unname(pre$col_map[[as.character(tid_i)]])
              }
              plotlyProxyInvoke(plotlyProxy("trk_track2D", session), "restyle",
                list(`line.color` = line_cols),
                seq(tm$trail_start, tm$trail_start + tm$n_tid - 1L))
            }

            # push updated axis ranges so filtered/new data is fully visible
            f1     <- pre$frames[[1]]
            bx     <- f1$markers_x[is.finite(f1$markers_x)]
            by     <- f1$markers_y[is.finite(f1$markers_y)]
            if (!length(bx) || !length(by)) {
              sp <- pre$sp
              bx <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE); bx <- bx[is.finite(bx)]
              by <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE); by <- by[is.finite(by)]
            }
            if (length(bx) > 0 && length(by) > 0) {
              xr <- range(bx); yr <- range(by)
              xpad <- max(1, diff(xr) * 0.04); ypad <- max(1, diff(yr) * 0.04)
              trk_axis_x_span(diff(xr) + 2 * xpad)
              trk_locked_ranges(list(x = c(xr[1]-xpad, xr[2]+xpad),
                                    y = c(yr[1]-ypad, yr[2]+ypad)))
              plotlyProxyInvoke(plotlyProxy("trk_track2D", session), "relayout", list(
                `xaxis.range[0]` = xr[1]-xpad, `xaxis.range[1]` = xr[2]+xpad,
                `xaxis.autorange` = FALSE,
                `yaxis.range[0]` = yr[1]-ypad, `yaxis.range[1]` = yr[2]+ypad,
                `yaxis.autorange` = FALSE
              ))
            }
          }, error = function(e) NULL)
        })
      }, delay = 0.2)
    } else {
      # fire_layout path: renderPlotly rebuilds the chart fresh. Push a grid/colour
      # relayout shortly after so the axis styling is consistent regardless of what
      # Plotly does internally during the DOM recreate + newPlot sequence.
      later::later(function() {
        shiny::withReactiveDomain(session$domain, {
          tryCatch(trk_axis_relayout(), error = function(e) NULL)
        })
      }, delay = 0.35)
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  observeEvent(input$trk_play_mode, { loop$playing <- FALSE })   # switching modes stops loop
  
  # ---- Live marker size update ----
  observeEvent(input$trk_pt_size, {
    req(!is.na(layout_trigger()))
    tm <- isolate(trace_map()); req(!is.null(tm))
    plotlyProxyInvoke(plotlyProxy("trk_track2D", session), "restyle",
                      list(`marker.size` = input$trk_pt_size), tm$markers_idx)
  })

  trk_axis_relayout <- function() {
    req(!is.na(layout_trigger()))
    hide  <- isTRUE(input$trk_hide_grid)
    dark  <- isTRUE(input$trk_dark_mode)
    is3d  <- isTRUE(use_3d())
    fg    <- if (dark) "#FFFFFF" else "#63666A"
    gc    <- if (dark) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
    bg    <- if (dark) "#000000" else "#FFFFFF"

    if (is3d) {
      # Recompute equal ranges so the relayout doesn't reset autorange.
      pre  <- isolate(frames_store())
      sp   <- if (is.list(pre)) pre$sp else list()
      all_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
      all_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
      all_z <- unlist(lapply(sp, `[[`, "z"), use.names = FALSE)
      lim3  <- get_equal_limits(all_x, all_y, all_z)

      ax3_style <- list(
        showgrid = !hide, zeroline = !hide, showticklabels = !hide,
        gridcolor = gc, color = fg, backgroundcolor = bg
      )
      make_ax3 <- function(title_txt, lim) {
        c(list(title = title_txt, range = lim, autorange = FALSE), ax3_style)
      }
      plotlyProxyInvoke(
        plotlyProxy("trk_track2D", session),
        "relayout",
        list(
          paper_bgcolor = bg,
          font = list(color = fg),
          scene = list(
            bgcolor = bg,
            xaxis = make_ax3(if (hide) "" else "X (nm)", lim3$x),
            yaxis = make_ax3(if (hide) "" else "Y (nm)", lim3$y),
            zaxis = make_ax3(if (hide) "" else "Z (nm)", lim3$z),
            aspectmode  = "manual",
            aspectratio = list(x = 1, y = 1, z = 1)
          )
        )
      )
    } else {
      ovl <- tryCatch(isolate(trk_image_overlay()), error = function(e) NULL)
      # Include locked ranges so the xaxis/yaxis object-form relayout doesn't
      # discard range + autorange that were set by the initial layout or proxy.
      lr         <- isolate(trk_locked_ranges())
      x_range_kv <- if (!is.null(lr)) list(range = lr$x, autorange = FALSE) else list()
      y_range_kv <- if (!is.null(lr)) list(range = lr$y, autorange = FALSE) else list()
      plotlyProxyInvoke(
        plotlyProxy("trk_track2D", session),
        "relayout",
        list(
          paper_bgcolor = bg, plot_bgcolor = bg,
          font = list(color = fg),
          xaxis = c(list(
            title = if (hide) "" else "X (nm)",
            showgrid = !hide, zeroline = !hide, showticklabels = !hide,
            gridcolor = gc, color = fg
          ), x_range_kv),
          yaxis = c(list(
            title = if (hide) "" else "Y (nm)",
            showgrid = !hide, zeroline = !hide, showticklabels = !hide,
            gridcolor = gc, color = fg,
            scaleanchor = "x", scaleratio = 1
          ), y_range_kv),
          images = if (!is.null(ovl)) list(ovl) else list()
        )
      )
    }

    # Restyle faint static trace color when dark mode toggles
    tm <- isolate(trace_map())
    if (!is.null(tm) && !is.null(tm$static_idx) && !is.na(tm$static_idx)) {
      faint_color <- if (dark) {
        if (is3d) "rgba(255,255,255,0.55)" else "rgba(255,255,255,0.25)"
      } else {
        if (is3d) "rgba(80,80,80,0.55)" else "rgba(80,80,80,0.15)"
      }
      plotlyProxyInvoke(
        plotlyProxy("trk_track2D", session),
        "restyle",
        list(`line.color` = faint_color),
        list(tm$static_idx)
      )
    }
  }

  observeEvent(input$trk_hide_grid,  { trk_axis_relayout() }, ignoreInit = TRUE)
  observeEvent(input$trk_dark_mode,  { trk_axis_relayout() }, ignoreInit = TRUE)

  # Scale bar: push via proxy whenever any scale-bar input changes.
  trk_push_scale_bar <- function() {
    req(!is.na(layout_trigger()))
    proxy <- plotlyProxy("trk_track2D", session)
    if (!isTRUE(input$trk_show_scale_bar)) {
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = list(), annotations = list()))
      return()
    }
    sb_len   <- suppressWarnings(as.numeric(input$trk_scale_bar_length))
    sb_thick <- suppressWarnings(as.numeric(input$trk_scale_bar_thickness))
    sb_font  <- suppressWarnings(as.numeric(input$trk_scale_font_size))
    sb <- make_scale_bar(sb_len, trk_axis_x_span() %||% 1000, sb_thick, dark = isTRUE(input$trk_dark_mode),
                         scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
    if (!is.null(sb))
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = sb$shapes, annotations = sb$annotations))
  }

  observe({
    # Only read the relayout event once the track plot has rendered + registered
    # AND the Tracking tab is on screen, so plotly's deferred "source not
    # registered" warning is never scheduled. trk_plot_ready stays FALSE until a
    # plot with data renders (layout_trigger alone can be set with no frames yet).
    if (!isTRUE(trk_plot_ready()) ||
        !identical(input$main_navbar, "tab_tracking")) return()
    req(!is.na(layout_trigger()))
    ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "track"))
    if (!is.list(ed)) return()

    x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
    x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
    if (length(x0) == 1 && length(x1) == 1 && is.finite(x0) && is.finite(x1) && x1 > x0) {
      y0 <- suppressWarnings(as.numeric(ed[["yaxis.range[0]"]]))
      y1 <- suppressWarnings(as.numeric(ed[["yaxis.range[1]"]]))
      if (length(y0) == 1 && length(y1) == 1 && is.finite(y0) && is.finite(y1) && y1 > y0) {
        # Capture only GENUINE user zoom/pan. During trail playback draw_frame re-pins
        # the view every frame (an identical-range relayout) and scaleanchor can nudge
        # it sub-% — both echo back as plotly_relayout events. Capturing those would
        # churn trk_zoom_ranges and let the re-pin race a fresh user zoom (snapping it
        # back). So only update when the range moved by > 1% of its span.
        cur <- isolate(trk_zoom_ranges())
        changed <- is.null(cur) || {
          sx <- max(1e-9, x1 - x0); sy <- max(1e-9, y1 - y0)
          abs(cur$x[1] - x0) > 0.01 * sx || abs(cur$x[2] - x1) > 0.01 * sx ||
          abs(cur$y[1] - y0) > 0.01 * sy || abs(cur$y[2] - y1) > 0.01 * sy
        }
        if (changed) { trk_axis_x_span(x1 - x0); trk_zoom_ranges(list(x = c(x0, x1), y = c(y0, y1))) }
      }
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      axis_state <- trk_default_axis_state()
      if (!is.null(axis_state)) trk_axis_x_span(axis_state$x_span)
      trk_zoom_ranges(NULL)
    } else {
      return()
    }

    if (isTRUE(input$trk_show_scale_bar)) trk_push_scale_bar()
  })

  # ---- Spatial crop: "Crop to view" / "Clear" --------------------------------
  observeEvent(input$trk_crop_to_view, {
    zr  <- trk_zoom_ranges()
    pct <- suppressWarnings(as.numeric(input$trk_crop_pct))
    if (is.null(zr)) {
      showNotification("Zoom into a region first, then press Crop to view.", type = "warning", duration = 5)
      return()
    }
    if (!is.finite(pct) || pct <= 0 || pct > 100) pct <- 50
    trk_spatial_filter(list(x = zr$x, y = zr$y, pct_thr = pct / 100))
    data_ver(isolate(data_ver()) + 1L)
    fire_frames()
    showNotification(
      sprintf("Spatial crop active: x [%.0f, %.0f] nm, y [%.0f, %.0f] nm (≥%.0f%% per TID).",
              zr$x[1], zr$x[2], zr$y[1], zr$y[2], pct),
      type = "message", duration = 5
    )
  }, ignoreInit = TRUE)

  observeEvent(input$trk_clear_spatial, {
    trk_spatial_filter(NULL)
    data_ver(isolate(data_ver()) + 1L)
    fire_frames()
    showNotification("Spatial crop cleared.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  # Info UI: show active box or a hint
  output$trk_spatial_info <- renderUI({
    sf <- trk_spatial_filter()
    if (is.null(sf)) {
      tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
             "No spatial filter active. Zoom into a region and press Crop to view.")
    } else {
      tags$p(style = "font-size:0.82em; color:#2070a0; margin:2px 0;",
             sprintf("Active: x [%.0f, %.0f] nm  |  y [%.0f, %.0f] nm  |  ≥%.0f%% per TID",
                     sf$x[1], sf$x[2], sf$y[1], sf$y[2], sf$pct_thr * 100))
    }
  })

  observeEvent(list(input$trk_show_scale_bar, input$trk_scale_bar_length,
                    input$trk_scale_bar_thickness, input$trk_scale_font_size,
                    input$trk_dark_mode),
    { trk_push_scale_bar() }, ignoreInit = TRUE)
} # end tracking_server



