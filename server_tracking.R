## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
## server_tracking.R   -   Tracking analysis server logic
## All input/output IDs prefixed with trk_ to avoid collisions with SMLM IDs.
## Called from server.R inside the combined server function.
## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##
## GOAL OF THIS FILE
##   Everything behind the "MINFLUXR - Tracking" tab: load a MINFLUX
##   single-particle-tracking dataset, clean it, animate it, and run the
##   diffusion analyses (MSD, CVE, JDD, HMM, confinement, D-map, ...).
##   It is one large Shiny server function, `tracking_server()`.
##
## HOW TO READ A SHINY SERVER FILE
##   Reactives are LAZY: the order they appear in the file is NOT the order they
##   run. What matters is which reactive calls which. So to follow the data,
##   follow the CALLS. The dependency chain is drawn below.
##
## ── 1. THE DATA FLOW (the single most important thing in this file) ──────────
##
##   upload (.json / .csv / .npy)
##        |  split into per-TID chunk files on disk + manifest.csv
##        |  (helpers_tracking.R::split_by_tid_to_dir - keeps RAM bounded)
##        v
##   trk_raw_df ................ loads EVERY chunk of the current manifest
##        v
##   ---- THE FILTER PIPELINE (each stage feeds the next; ~line 1171) ----------
##   trk_gap_split ............. Stage 1  cut tracks at dark-cycle gaps
##   trk_end_trimmed ........... Stage 2  peel dim / gap-reached track ends
##   trk_quality_filtered ...... Stage 3  drop dark-cycle-heavy tracks
##   trk_meta_filtered ......... Stage 4  efo / cfr / dcr / eco / dwell ranges
##   trk_len_filtered .......... Stage 5  localizations per track
##   trk_dur_filtered .......... Stage 6  track duration
##   trk_spatial_filtered ...... Stage 7  spatial crop
##        v
##   trk() ..................... THE analysis-ready dataset. Almost every
##                               analysis below starts by calling trk().
##
##   Cleaning stages are HEAVY, so they read their parameters via isolate() and
##   only recompute when the "apply cleaning" trigger fires 
##
## ── 2. THE OTHER FOUNDATION: TIME ────────────────────────────────────────────
##   MINFLUX has no fixed clock. Almost every estimator here needs the base
##   acquisition CYCLE TIME to convert irregular timestamps into integer cycle
##   counts. That is what makes the analyses comparable.
##     trk_cycle_time .......... the resolved cycle time (hardware > manual >
##                               data-driven cluster > regression)
##     trk_cluster_cycle ....... data-driven estimate (robust to dark cycles)
##     n_cyc = round(dt / t_c) . cycles spanned by a step; 1 = "integrated"
##   A step with n_cyc > 1 skipped (n_cyc - 1) "dark" cycles.
##
## ── 2b. THE MATH IN ONE PLACE (every estimator below rests on these) ─────────
##
##   NOTATION   r = |Δr| step length (nm) | Δt = elapsed time (ms)
##              d = dimensions (2 or 3)   | D = diffusion coefficient (nm²/ms)
##              σ = σ_loc localization precision (nm)
##
##   (1) FREE DIFFUSION.  Squared displacement grows LINEARLY with time:
##
##           <r²> = 2·d·D·Δt
##
##       so the displacement itself grows as the SQUARE ROOT of time:
##           r ~ sqrt(2·d·D·Δt)
##
##   (2) WHY WE DIVIDE BY sqrt(Δt)
##       MINFLUX steps span different Δt, so before pooling them we must remove
##       the Δt dependence. From (1) that dependence is sqrt(Δt):
##
##           r / sqrt(Δt)  ->  scale sqrt(2dD)   CONSTANT in Δt   (what we use)
##           r / Δt        ->  a VELOCITY, scales as sqrt(2dD/Δt)
##                             -> DIVERGES as Δt -> 0             (never use)
##
##       Brownian motion is continuous but nowhere differentiable: it has no
##       instantaneous velocity. A "speed" r/Δt therefore measures your sampling
##       rate as much as the molecule. (This is exactly why the Speed panel's
##       path speed cannot separate diffusion from transport.)
##
##       Equivalent forms - the same operation, two representations:
##           (r / sqrt(Δt))² = r² / Δt  ∝ D
##       So: JDD + HMM normalise LENGTHS by sqrt(Δt) (Rayleigh-distributed),
##           D-map + change-point normalise SQUARED lengths by Δt (local D).
##
##       CAVEAT: this invariance assumes FREE diffusion. Under anomalous motion
##       <r²> ~ Δt^α the correct exponent is Δt^(α/2), so sqrt(Δt) over-corrects
##       for α<1. And localization noise does NOT scale with sqrt(Δt) (it is
##       constant), leaving a Δt-dependent noise term - hence the σ_loc
##       corrections, and why "integrated steps only" (one fixed Δt, no
##       normalisation needed) is the cleanest mode where it is available.
##
##   (3) LOCALIZATION NOISE adds a constant offset:
##
##           <r²>_measured = 2·d·D·Δt  +  2·d·σ²
##
##       -> MSD fit: SLOPE gives D, INTERCEPT gives σ.
##       -> An immobile molecule still shows apparent D = σ²/Δt. That is the
##          "noise floor" the D-map and change-point subtract.
##       -> A pure-noise step has RMS length sqrt(2d)·σ.
##
##   (4) CYCLE COUNTS. MINFLUX has no fixed clock, so time is measured in whole
##       acquisition cycles: n_cyc = round(Δt / t_c). n_cyc = 1 is "integrated"
##       (no gap); n_cyc > 1 skipped (n_cyc − 1) dark cycles.
##
##   UNITS: nm and ms internally. D in nm²/ms; divide by 1e3 for µm²/s.
##
## ── 3. MAP OF THE FILE (grep the name; line numbers drift, names do not) ─────
##   Setup / state ............. reactiveValues + reactiveVal block, top
##   sigma_loc (noise floor) ... trk_sigma_loc_eff, trk_sigma_loc_src
##   Calibration (CVE) ......... trk_cve_bundle          ~154
##   Upload & chunking ......... trk_raw_df, trk_raw_original ~870
##   FILTER PIPELINE ........... banner at                ~1171
##   Cycle time & TSA .......... banner at                ~1382
##   Playback / animation ...... seq_data, frames_store   ~2925
##   MSD (main engine) ......... msd_unified_bundle       ~6202
##     per-track D / ergodicity  trk_unified_tdiff        ~6410
##   Change-point + GMM ........ cpd_segments_rv          ~4616
##   Jump-distance (JDD) ....... jdd_bundle               ~7078
##   Confinement index (MSD) ... conf_ci_bundle           ~8060
##     (the Rg/Rmax shape test was retired 2026-08-03 -> retired_tracking_server.R)
##   Directionality ............ dir_bundle               ~7627
##   Speed & transport ......... trk_export_summary       ~7978
##   Diffusion-state HMM ....... hmm_bundle               ~8298
##   Per-state MSD ............. seg_msd_rv               ~8826
##   Spatial D-map ............. dmap_bundle              ~8982
##
## ── 4. CONVENTIONS & GOTCHAS (things that have bitten us) ────────────────────
##   * `tim` is MILLISECONDS in here (seconds on disk, x1000 right after load).
##     Coordinates are nanometres.
##   * "TID" = trace ID = one tracked molecule. Always character.
##   * `*_bundle` / `*_rv` reactives are eventReactive: they only recompute when
##     their Run button is pressed. That is deliberate - they are expensive.
##   * Chunking is PLAYBACK-only. The analyses always run on all TIDs; the
##     chunk selector just limits what gets animated so playback stays smooth.
##   * jsonlite is not attached - always write jsonlite:: (its `validate` would
##     shadow shiny::validate).
##   * Some sprintf() format strings contain literal \u escapes (e.g. °)
##     rather than the character itself. Match them exactly when editing.
##
##   Method provenance and the MINFLUX-appropriateness of every estimator are
##   documented in docs/MINFLUX_analysis_methods.md.
## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

tracking_server <- function(input, output, session) {
  
  
  # --- Reactives / state ---
  rv <- reactiveValues(
    prev_tid_idx          = NULL,
    play_gen              = 0L,   # global playback generation (simultaneous mode)
    manifest_dir          = NULL, # current working chunks (may be filtered)
    original_manifest_dir = NULL, # frozen at upload - the rechunk always reads from here
    last_rechunk_dir      = NULL, # previous rechunk output, deleted once superseded
    bigfile_notice        = FALSE,
    chunking              = FALSE,
    trk_diag              = NULL  # list(data_clusters, img_blobs) in data coords
  )
  # ── Chunk size actually used ────────────────────────────────────────────────
  # GOAL: turn the "TIDs per chunk" input into the value the chunker gets, for
  # both the upload path and Filter & Rechunk, so the two cannot disagree.
  #
  # A DELIBERATELY SMALL SETTING IS HONOURED. Small chunks are the whole point of
  # the control: playback animates one chunk at a time, so a modest machine that
  # cannot render 100 tracks needs 30. Only a blank or non-numeric entry falls
  # back to the default. The floor is TPC_MIN rather than the default, so a value
  # below it is clamped up by a little instead of being thrown back to 100.
  #
  # WHY THERE IS A FLOOR AT ALL: each chunk becomes one file on disk, so a value
  # near 1 would emit thousands of files for a large dataset and make every
  # manifest read slow.
  TPC_MIN <- 5L
  .tids_per_chunk <- function() {
    v <- suppressWarnings(as.integer(isolate(input$trk_tids_per_chunk)))
    if (!length(v) || !is.finite(v[1]) || v[1] < 1L) return(TIDS_PER_CHUNK_DEFAULT)   # unset / blank / garbage
    max(TPC_MIN, v[1])
  }

  # ── Which chunk the selector lands on when it is rebuilt ────────────────────
  # "__all__" = animate every TID; "first" = land on chunk 1.
  #
  # Rebuilding the selector (upload, Filter & Rechunk, Reset ranges) has to pick a
  # value, and the previous selection cannot simply be carried over: re-chunking at
  # a different size renumbers the chunks, so "chunk 5" may not exist afterwards and
  # would not hold the same tracks if it did. What IS meaningful is whether the user
  # was scoped to a chunk at all, so that is what is preserved - someone who was on a
  # chunk stays on a chunk (the first one, which always exists), and someone on
  # "All TIDs" stays on all. Landing a machine that could not render everything back
  # on "All TIDs" is what this avoids.
  trk_chunk_pref <- reactiveVal("__all__")

  # Plot Data asks for a rechunk through this counter rather than calling the
  # routine directly, so the 260-line filter/rechunk observer stays exactly where
  # it is and keeps its own progress bar and error handling.
  rechunk_request <- reactiveVal(0L)

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
  # list(x = c(xmin, xmax), y = c(ymin, ymax)) - used to lock plot axes.
  trk_rechunk_bbox      <- reactiveVal(NULL)
  current_pt_idx     <- reactiveVal(1L)
  # Current plot axis ranges captured from plotly_relayout - used by "Crop to view".
  # list(x = c(x0, x1), y = c(y0, y1)) or NULL before first zoom interaction.
  trk_zoom_ranges    <- reactiveVal(NULL)
  # Active spatial box filter - NULL means no spatial filter is applied.
  # list(x = c(x0, x1), y = c(y0, y1), pct_thr = numeric[0-1])
  trk_spatial_filter <- reactiveVal(NULL)

  # "TSA has been run at least once" gate. The Integrated/Dark bundle is
  # eventReactive-gated on trk_tsa_run, but the sibling display outputs (Δt
  # histogram, stats, quality verdict, per-track progression) read dt_vectors()
  # directly - this flag keeps them dormant until Run is pressed. Cleared when TSA
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
  # (which ranks all tracks by export strength - heavy) and the plot stay dormant
  # until Run is pressed, then remain interactive on TID selection.
  speed_ran <- reactiveVal(FALSE)
  observeEvent(input$trk_speed_run, { speed_ran(TRUE) })
  observeEvent(input$trk_show_speed, {
    if (!isTRUE(input$trk_show_speed)) speed_ran(FALSE)
  }, ignoreInit = TRUE)

  # ── App-wide localization precision σ_loc (the short-lag noise floor) ────────
  # Two independent sources feed the app-wide value, in priority order:
  #   1. manual override (input$trk_sigma_loc_nm > 0)
  #   2. the photon-budget estimate, L_TCP/(4*sqrt(eco))  (trk_sigma_loc_photon)
  #   3. the blur-aware CVE estimate                      (trk_sigma_loc_cve)
  #   4. the unified MSD fit's intercept                  (trk_sigma_loc_rv)
  # Every short-lag tool (D-map, JDD, confinement, HMM, directionality, speed)
  # reads trk_sigma_loc_eff() so the noise floor is consistent everywhere.
  #
  # WHY THE PHOTON MODEL OUTRANKS THE CVE. The CVE σ is dominated by its
  # R*⟨Δx²⟩ blur term (~90% of σ² at R ≈ 1/6), so it mostly reports sqrt(R) x rms
  # step - sensitive to the assumed R, and inflated by population heterogeneity and
  # drift. L/(4*sqrt(N)) is a per-localization instrument model that assumes nothing
  # about the motion. The CVE is kept as an independent cross-check on L.
  trk_sigma_loc_rv     <- reactiveVal(NA_real_)   # nm, from the last unified MSD fit
  trk_sigma_loc_cve    <- reactiveVal(NA_real_)   # nm, blur-aware CVE (cross-check)
  trk_sigma_loc_photon <- reactiveVal(NA_real_)   # nm, median L/(4*sqrt(eco))
  trk_sigma_loc_eff <- reactive({
    manual <- suppressWarnings(as.numeric(input$trk_sigma_loc_nm %||% 0))
    if (is.finite(manual) && manual > 0) return(manual)
    ph <- trk_sigma_loc_photon()
    if (is.finite(ph) && ph > 0) return(ph)
    cve <- trk_sigma_loc_cve()
    if (is.finite(cve) && cve > 0) return(cve)
    s <- trk_sigma_loc_rv()
    if (is.finite(s) && s > 0) s else NA_real_
  })
  # Which source is currently active (for status readouts).
  trk_sigma_loc_src <- reactive({
    manual <- suppressWarnings(as.numeric(input$trk_sigma_loc_nm %||% 0))
    if (is.finite(manual) && manual > 0) "manual override"
    else if (is.finite(trk_sigma_loc_photon())) "photon budget"
    else if (is.finite(trk_sigma_loc_cve())) "CVE (cross-check)"
    else if (is.finite(trk_sigma_loc_rv()))  "MSD fit"
    else "not set"
  })

  # Mean + SD + t-based 95% CI treating each element (a TRACK) as the independent
  # unit - the honest CI for track-level statistics, since localizations within a
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
  # correlated - the track is the independent unit, not the localization.
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
  .continuity_scores <- continuity_scores                  # R/est_clean.R

  # Does the z column carry REAL axial data, or is it a placeholder?
  # 2D MINFLUX exports routinely ship an all-zero (or constant) z column, and 0 is
  # finite - so the obvious test `any(is.finite(df$z))` returns TRUE for 2D data and
  # every dimension-aware estimator downstream silently switches to 3D. That pools a
  # ZERO-VARIANCE axis into the maths: it diluted the pooled σ_loc by √(2/3) (−18%)
  # and D_cve by 2/3 (−33%), and it would have the JDD fit a 3D Maxwell law to 2D
  # distances. Requiring real spread is the criterion msd_unified_bundle and
  # trk_disp_bundle already used; this makes it the single shared definition so the
  # dimension decision cannot drift between tools again.
  # Tracking's 3D test, kept as a named wrapper because it has many call sites.
  # Delegates to the shared predicate so the axial-calibration panel and the
  # 2D/3D view toggle can never disagree about whether a dataset has z.
  .has_real_z <- function(df, min_sd = 1) {
    has_axial_data(df, min_sd = min_sd, min_n = 11L)
  }

  # ============================================================================
  # trk_cve_bundle() - LOCALIZATION PRECISION (σ_loc) + DATA-QUALITY CALIBRATION
  # ============================================================================
  #
  # GOAL
  #   Measure how precisely the microscope located each molecule (σ_loc, nm) and
  #   score how clean the acquisition timing is. σ_loc becomes the app-wide
  #   "noise floor" that the D-map, JDD, confinement, HMM, directionality and
  #   speed tools all subtract. Run this FIRST - everything else calibrates off it.
  #
  # THE IDEA (CVE - Vestergaard, Blainey & Flyvbjerg 2014)
  #   Localization error appears in TWO consecutive steps with OPPOSITE sign: the
  #   noise at point i is the END of step i-1 and the START of step i. That makes
  #   neighbouring steps anti-correlated by exactly −σ², while real diffusion
  #   increments are independent and contribute nothing to that correlation. So
  #
  #       σ_loc² = −< Δx_i · Δx_{i+1} >          (per axis, averaged over x,y[,z])
  #
  #   That is the R = 0 (instantaneous exposure) case. A MINFLUX localization
  #   integrates photons across the TCP cycle, so the code below uses the blur-aware
  #   form σ² = R·<Δx²> + (2R − 1)·cov₁ (see "Motion blur" in the bundle).
  #   Bonus: D falls out too, noise-corrected, with no MSD fit at all:
  #       D = <Δ²>/(2Δt) + cov₁/Δt
  #
  # THE MINFLUX ADAPTATION
  #   The algebra above assumes a CONSTANT Δt. MINFLUX Δt is irregular. So the
  #   covariance is taken ONLY over pairs where BOTH steps are single-cycle
  #   ("integrated", n_cyc == 1) and share a midpoint inside one track. That
  #   reconstructs constant-Δt sub-runs, at the cost of discarding steps next to
  #   a dark gap.
  #
  # OBJECTS IT CALLS (inputs)
  #   trk_gap_split() ......... data after gap-splitting (cleaning stage 1)
  #   trk_spatial_filtered() .. the fully filtered data, for the "after filters" score
  #   trk_raw_original() ...... the frozen upload, for the "raw" reference score
  #   trk_cycle_time() ........ base cycle time, to classify integrated vs dark
  #   .continuity_scores(), .boot_ci()
  #   inputs: trk_cve_run (Run trigger), trk_sig_min_integ, trk_sig_skip,
  #           trk_quality_darktol
  #
  # STEP BY STEP
  #   1. Build per-step Δx, Δy[, Δz] and Δt within each track (never across tracks).
  #   2. Classify each step: n_cyc = round(Δt / t_c); n_cyc == 1 is "integrated".
  #   3. Data-quality scores: fraction of integrated steps, raw vs after-filters,
  #      plus a dark-tolerance-weighted continuity score.
  #   4. POOLED σ_loc: average Δx_i·Δx_{i+1} over consecutive integrated pairs,
  #      negate, square-root. Also the bias-corrected D_cve.
  #   5. PER-TRACK σ_loc: repeat per track (tracks with enough integrated steps),
  #      skipping the first few localizations - MINFLUX needs a moment to lock on
  #      and those opening points read worse. Gives median, IQR, CV, bootstrap CI,
  #      and a σ-vs-EFO correlation (precision SHOULD improve with photons; a
  #      negative correlation validates the whole estimate).
  #   6. PHOTON WEIGHTING (optional): fit σ ∝ EFO^p across tracks, assign every
  #      localization a precision from its own EFO, and inverse-variance-weight
  #      the pooled σ and D so bright, well-localized steps count for more.
  #   7. QC flags: implausible jumps (|Δr| > 5× expected for that gap - likely
  #      mis-links) and suspicious tracks (degenerate or outlier σ).
  #
  # RETURNS a list. The fields consumers actually use:
  #   sigma ......... pooled σ_loc (nm)
  #   sig_pt_med .... per-track median σ_loc    D_cve ..... noise-corrected D
  #   degenerate .... TRUE if σ² came out <= 0 (drift dominates -> unresolvable)
  #   frac_integ, cont_raw, cont_filt ... data-quality scores
  #   pertrack_sigma .... per-track table (tid, sigma, n_integ, phot_med, var_step)
  #   r1_norm / noise_frac ... how much of the step variance is localization noise
  #   phot_mob_rho ... ρ[photons, step variance]: negative = photons track MOBILITY
  #   n_implaus, n_suspect, ...
  #
  # WHEN IT FAILS: strong directed motion makes consecutive steps POSITIVELY
  # correlated, so σ² goes negative and σ_loc is genuinely unresolvable. That is
  # reported as `degenerate` - enter σ_loc manually in that case.
  # ============================================================================
  trk_cve_bundle <- eventReactive(input$trk_cve_run, {
    # The filtered data of the active dataset, as every other analysis reads it.
    df <- tryCatch(trk_spatial_filtered(), error = function(e) NULL)
    req(!is.null(df), nrow(df) > 0, all(c("x","y","tid","tim") %in% names(df)))
    t_c   <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    # Real axial data only - a placeholder all-zero z column would dilute both σ_loc
    # and D_cve (see .has_real_z).
    has_z <- .has_real_z(df)
    # Step arrays from est_cve.R, shared with cve_core() and the diagnostics below.
    st <- cve_steps(df, t_c, has_z)
    df <- st$df; n <- st$n; same <- st$same; dt <- st$dt
    dx <- st$dx; dy <- st$dy; dz <- st$dz; n_cyc <- st$n_cyc
    integ_step <- st$integ_step; n_steps <- st$n_steps; frac_integ <- st$frac_integ

    # Same integrated fraction, but on the FULLY-FILTERED data (post quality-gate /
    # end-trim / length / duration / spatial) - the data the analyses actually use.
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
    # (= trk_gap_split, which reads the current - possibly filtered - manifest). This
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

    # Collection window from the photon budget: t_coll = ECO / EFO.
    .med_of <- function(cn) if (cn %in% names(df)) {
                              v <- df[[cn]][is.finite(df[[cn]]) & df[[cn]] > 0]
                              if (length(v)) stats::median(v) else NA_real_
                            } else NA_real_
    eco_m <- .med_of("eco"); efo_m <- .med_of("efo")
    t_coll <- if (is.finite(eco_m) && is.finite(efo_m) && efo_m > 0) eco_m / efo_m * 1000 else NA_real_
    n_vert <- {
      sp <- tryCatch(if (is.null(input$trk_seq_json)) NULL else trk_seq_params(), error = function(e) NULL)
      v  <- if (!is.null(sp) && isTRUE(sp$ok)) sp$vertices else NA_integer_
      if (!is.finite(v)) suppressWarnings(as.integer(input$trk_hw_vertices %||% 6L)) else v
    }
    if (!is.finite(n_vert) || n_vert < 1L) n_vert <- 6L
    # CVE core (covariance, blur bounds, σ, D) lives in est_cve.R.
    cc <- cve_core(dx, dy, dz, dt, integ_step, has_z, t_c, t_coll, n_vert,
                   R_step = trk_blur()$R_step)
    shiny::validate(shiny::need(!is.null(cc),
      "Too few consecutive integrated steps for a CVE estimate - relax filters or check the cycle time."))
    pair <- cc$pair; prod_x <- cc$prod_x; prod_y <- cc$prod_y; prod_z <- cc$prod_z
    vals <- rep(NA_real_, cc$n_pairs)                       # only its length is read below
    cov1 <- cc$cov1; dt_int <- cc$dt_int; R_lo <- cc$R_lo; R_hi <- cc$R_hi; R_blur <- cc$R_blur
    msq <- cc$msq; sigma2 <- cc$sigma2; sigma <- cc$sigma
    sigma_xy <- cc$sigma_xy; sigma_z <- cc$sigma_z; sig_aniso <- cc$sig_aniso
    D_cve <- cc$D_cve; r1_norm <- cc$r1_norm; noise_frac <- cc$noise_frac

    # ── Per-track σ_loc heterogeneity ──────────────────────────────────────────
    # One CVE σ per track (tracks with ≥ n_min integrated steps), skipping the
    # first `skip` localizations - MINFLUX iterates onto the molecule, so the
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
    # ECO = photon counts (correct for σ ∝ 1/√N); EFO = emission frequency (fallback).
    phot_col <- if ("eco" %in% names(df) && any(is.finite(df$eco) & df$eco > 0)) "eco"
                else if ("efo" %in% names(df) && any(is.finite(df$efo) & df$efo > 0)) "efo"
                else NA_character_
    # Per-track σ carries the same blur correction. Note what it can and cannot
    # say: at R ≈ 1/6 the R*⟨Δx²⟩ term is ~90% of σ², so per-track differences here
    # track the per-track STEP SIZE (mobility) far more than the photon budget.
    # Kept for the outlier/degenerate QC below, not as a precision measurement -
    # the per-localization photon estimate is what the panel plots.
    sq_x_pt <- dx^2; sq_y_pt <- dy^2; sq_z_pt <- if (has_z) dz^2 else rep(NA_real_, n)
    msq_by_tid <- {
      v <- c(sq_x_pt[integ_post], sq_y_pt[integ_post])
      t2 <- c(df$tid[integ_post], df$tid[integ_post])
      if (has_z) { v <- c(v, sq_z_pt[integ_post]); t2 <- c(t2, df$tid[integ_post]) }
      ok <- is.finite(v); tapply(v[ok], t2[ok], mean)
    }
    pertrack_sigma <- if (length(pv)) {
      cov_by <- tapply(pv, tv, mean)                    # per-track cov1
      s2_by  <- vapply(names(cov_by), function(nm) {
                  mq <- msq_by_tid[[nm]]
                  if (is.null(mq) || !is.finite(mq)) -1
                  else R_blur * mq + (2 * R_blur - 1) * cov_by[[nm]]
                }, numeric(1))
      n_by   <- tapply(integ_post, df$tid, sum)[names(s2_by)]
      # Photon predictor for the precision model. σ_loc ≈ L_TCP/(2√N) depends on the
      # photon COUNT N, so ECO (effective counts) is the physically correct variable.
      # EFO is a FREQUENCY (emission rate); it only tracks N when the dwell time is
      # constant, which MINFLUX tracking violates (dwell varies, especially around dark
      # cycles), so it is a fallback only. Using ECO also gives a testable prediction:
      # the log-log slope should be ≈ -0.5.
      phot_by <- if (!is.na(phot_col)) {
                   kl <- post & is.finite(df[[phot_col]]) & df[[phot_col]] > 0
                   tapply(df[[phot_col]][kl], df$tid[kl], stats::median)[names(s2_by)]
                 } else stats::setNames(rep(NA_real_, length(s2_by)), names(s2_by))
      # Per-track step variance - used to test whether the photon channel is really a
      # MOBILITY proxy (bright = slow), which would make photon weighting bias D.
      v_by <- tapply(dx[integ_post]^2 + dy[integ_post]^2,
                     df$tid[integ_post], mean, na.rm = TRUE)[names(s2_by)]
      data.frame(tid = names(s2_by), sigma = sqrt(pmax(0, as.numeric(s2_by))),
                 n_integ = as.integer(n_by), phot_med = as.numeric(phot_by),
                 var_step = as.numeric(v_by),
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
    # 95% CI on the per-track-median σ_loc (bootstrap over tracks).
    sig_ci <- if (length(sig_pt_ok) >= 5L) .boot_ci(sig_pt_ok, stats::median, B = 1000L) else c(lo = NA_real_, hi = NA_real_)

    # ── Linking / sanity QC ─────────────────────────────────────────────────────
    # (1) Implausible jumps: a step whose |Δr| far exceeds what free diffusion +
    #     localization noise allow for its OWN gap (E[|Δr|²] = 2d·D·Δt + 2d·σ²).
    #     A big jump after one cycle is suspicious; after a long dark gap it isn't
    #     (the Δt term grows), so this naturally accounts for dark cycles. Flags
    #     mis-links / re-acquisition errors. (2) Suspicious tracks: per-track σ_loc
    #     that is degenerate (σ²<0) or a far outlier - candidate bad/mis-linked tracks.
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

    # ── Photon channel: DIAGNOSTIC ONLY (photon weighting has been REMOVED) ─────
    # There used to be an inverse-variance "photon-weighted" CVE here (fit σ ∝ photons^p,
    # then weight steps by the predicted σ). It is gone, because on MINFLUX TRACKING data
    # the photon channel is structurally coupled to MOBILITY, not just brightness: a fast
    # molecule leaves the beam pattern sooner, so it collects FEWER photons per
    # localization. Weighting by photons therefore up-weights slow molecules and biases D
    # downward - the estimator was measuring the thing it was supposed to correct.
    # We keep the two quantities that revealed it, as reported diagnostics:
    #   noise_frac  - how much of the step variance is localization noise (−2·r₁). Small
    #                 ⇒ per-track σ differences are motion, not photons.
    #   phot_mob_rho - Spearman ρ[per-track photons, per-track step variance]. Strongly
    #                 negative ⇒ the photon channel is a mobility proxy.
    # The pooled σ_loc / D_cve above are unaffected: they never used the photon model.
    phot_mob_rho <- if (!is.null(pertrack_sigma) &&
                        sum(is.finite(pertrack_sigma$phot_med) &
                            is.finite(pertrack_sigma$var_step)) >= 8L)
      suppressWarnings(stats::cor(pertrack_sigma$phot_med, pertrack_sigma$var_step,
                                  method = "spearman", use = "complete.obs")) else NA_real_

    # ── σ_loc from the photon budget - the value the app actually uses ─────────
    # σ_i = L_TCP / (a*sqrt(N_i)) per localization, N = eco, a the measured effective
    # prefactor, plus the pattern's off-centre term (.crlb_doughnut, est_cve.R). It is
    # per-localization rather than pooled, and does not inherit the step-size
    # inflation that population heterogeneity and drift push into the CVE's ⟨Δx²⟩
    # term. The CVE σ above is kept as a diagnostic only.
    # ── Are the increments consistent with ANY motion model? ───────────────────
    # A subdiffusive walk already has anti-correlated increments before any noise:
    # for fBm with Hurst H = alpha/2, fGn gives cov1/var = 2^(alpha-1) - 1.
    # Localization noise and blur can only push r1 MORE negative, so
    #     r1_observed <= 2^(alpha-1) - 1
    # must hold. A violation means something is adding POSITIVE increment
    # correlation, i.e. the trajectory was smoothed before it was written - galvo
    # damping (DMP > 0) being the obvious candidate. No estimator can recover
    # sigma from filtered increments, so this is a data-quality flag, not a knob.
    sm <- smoothing_check(df, t_c, dt_int)              # R/est_cve.R
    alpha_short <- sm$alpha_short; r1_fgn <- sm$r1_fgn

    L_tcp <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
    if (!is.finite(L_tcp) || L_tcp <= 0) L_tcp <- DMAX_L_TCP_DEFAULT
    ecov  <- if ("eco" %in% names(df)) df$eco[is.finite(df$eco) & df$eco > 0] else numeric(0)
    # The emitter drifts off the TCP centre while photons are collected; D from the
    # CVE above and the single-cycle interval set how far. The pattern sets how
    # fast precision falls off with that offset (TCP_PATTERN_GEOM).
    nu_beam <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
    if (!is.finite(nu_beam) || nu_beam <= 0) nu_beam <- TCP_BEAM_FWHM_NM
    tcp_pat <- input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT
    D_off   <- if (is.finite(D_cve) && D_cve > 0) D_cve else 0
    t_off   <- if (is.finite(dt_int) && dt_int > 0) dt_int else 0
    r_off   <- sqrt(.tcp_geom(tcp_pat)$dim * D_off * t_off)   # rms offset over the exposure, nm
    a_pref <- .a_in_use(tcp_pat)
    sig_phot_i   <- if (length(ecov)) .crlb_doughnut(L_tcp, nu_beam, ecov, D = D_off, t_loc = t_off,
                                                     pref = a_pref, pattern = tcp_pat) else numeric(0)
    sig_phot_med <- if (length(sig_phot_i)) stats::median(sig_phot_i) else NA_real_
    sig_phot_iqr <- if (length(sig_phot_i) >= 5L)
                      as.numeric(diff(stats::quantile(sig_phot_i, c(.25, .75)))) else NA_real_
    eco_med      <- if (length(ecov)) stats::median(ecov) else NA_real_
    # Invert the blur-aware CVE for L and compare against the sequence's L. The two
    # routes share no assumptions, so agreement is a real check on the prefactor.
    # The prefactor itself is calibrated on an immobile sample; see docs/sigma_loc_model.md.

    list(sigma = sigma, sigma2 = sigma2, cov1 = cov1, n_pairs = length(vals),
         R_blur = R_blur, msq = msq, n_vert = n_vert, t_coll = t_coll,
         sig_phot = sig_phot_i, sig_phot_med = sig_phot_med, sig_phot_iqr = sig_phot_iqr,
         eco_med = eco_med, L_tcp = L_tcp, n_phot = length(ecov),
         nu_beam = nu_beam, r_off = r_off, tcp_pattern = tcp_pat,
         sig_simple = if (length(ecov)) stats::median(L_tcp / (a_pref * sqrt(ecov))) else NA_real_,
         a_pref = a_pref,
         alpha_short = alpha_short, r1_fgn = r1_fgn,
         sigma_xy = sigma_xy, sigma_z = sigma_z, sig_aniso = sig_aniso,
         D_cve = D_cve, dt_int = dt_int, t_cycle = t_c, has_z = has_z,
         frac_integ = frac_integ_raw, frac_integ_filt = frac_integ_filt,
         cont_raw = cont_raw, cont_filt = cont_filt, cont_H = H_disp,
         n_steps = n_steps, degenerate = sigma2 <= 0,
         pertrack_sigma = pertrack_sigma, sig_pt_med = sig_pt_med,
         sig_pt_iqr = sig_pt_iqr, sig_pt_cv = sig_pt_cv, sig_pt_lo = sig_ci[["lo"]],
         sig_pt_hi = sig_ci[["hi"]], n_pt = length(sig_pt_ok),
         n_implaus = n_implaus, frac_implaus = frac_implaus,
         n_suspect = n_suspect, n_tracks_pt = n_tracks_pt,
         phot_col = phot_col,   # photon column used for the diagnostics
         r1_norm = r1_norm, noise_frac = noise_frac, phot_mob_rho = phot_mob_rho)
  })

  # Publish the calibration result as the app-wide σ_loc.
  #
  # The PHOTON estimate is the one used. The CVE σ is not: even with the blur
  # correction applied, ~90% of its σ² comes from the R*⟨Δx²⟩ term, so it is
  # essentially sqrt(R) x rms step size - a function of the assumed R and of the
  # step distribution, which population heterogeneity and drift both inflate. The
  # data's own contribution, cov1, is a few-percent residual. It stays on screen as
  # a cross-check on L, not as a precision measurement.
  observeEvent(trk_cve_bundle(), {
    b <- trk_cve_bundle(); if (is.null(b)) return()
    if (is.finite(b$sig_phot_med) && b$sig_phot_med > 0) trk_sigma_loc_photon(b$sig_phot_med)
    val <- if (is.finite(b$sigma) && b$sigma > 0 && !isTRUE(b$degenerate)) b$sigma else NA_real_
    if (is.finite(val)) trk_sigma_loc_cve(val)
  }, ignoreInit = TRUE)

  # GOAL: report the CVE calibration as numbers plus short verdicts.
  # The reasoning that used to be printed here on every run (why photons are not
  # used as weights, what sigma_loc feeds downstream, why only integrated pairs
  # qualify) now lives in the panel's collapsible Method box and in README
  # "CVE (Covariance-Based Estimator)". On-screen text stays numeric.
  output$trk_cve_status <- renderUI({
    b <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    if (is.null(b))
      return(htmltools::div(style = "color:#888; font-size:0.9em;",
        "Press ", htmltools::tags$b("Estimate ", sigloc(), " & quality"), " to calibrate."))

    AMB <- "#B8860B"; RED <- "#c03030"; DIM <- "#888"
    row  <- function(...) htmltools::div(style = "margin:1px 0;", ...)
    dim_ <- function(...) htmltools::span(style = sprintf("color:%s;", DIM), ...)
    warn <- function(...) htmltools::span(style = sprintf("color:%s;", AMB), ...)
    num  <- function(x, fmt, na = "-") if (is.finite(x)) sprintf(fmt, x) else na

    # ── sigma_loc ──────────────────────────────────────────────────────────────
    # The PHOTON estimate is the value the app uses; the CVE follows as a check.
    sigma_line <- if (is.finite(b$sig_phot_med))
      row(htmltools::tags$b(sigloc(), sprintf(" %.2f nm", b$sig_phot_med)),
          dim_(sprintf("  photon budget  \u00b7  %s  \u00b7  L = %.1f nm  \u00b7  a = %.2f%s  \u00b7  offset %.1f nm  \u00b7  median ECO %.0f  \u00b7  %s locs",
                       b$tcp_pattern %||% TCP_PATTERN_DEFAULT, b$L_tcp, b$a_pref %||% NA_real_,
                       if (is.finite(b$nu_beam %||% NA_real_) && b$nu_beam > 0)
                         sprintf("  \u00b7  \u03bd = %.0f nm", b$nu_beam) else "",
                       b$r_off %||% 0,
                       b$eco_med %||% NA_real_, format(b$n_phot %||% 0L, big.mark = ","))),
          if (is.finite(b$sig_simple %||% NA_real_))
            dim_(sprintf("  \u00b7  %.2f before the motion term", b$sig_simple)) else NULL)
    else
      row(warn(htmltools::HTML("&#9888; "), "no ECO column, so ", sigloc(),
               " cannot come from the photon budget; the CVE value below is used instead."))

    # Consistency of the increments with ANY motion model. Subdiffusion alone gives
    # r\u2081 = 2^(\u03b1\u22121) \u2212 1; noise and blur only push r\u2081 further negative. If the observed
    # r\u2081 sits ABOVE that, the trajectory was smoothed before it was written (galvo
    # damping), and no estimator can recover \u03c3 from it.
    smooth_line <- if (is.finite(b$r1_fgn %||% NA_real_) && is.finite(b$r1_norm %||% NA_real_)) {
      viol <- b$r1_norm > b$r1_fgn + 0.01
      row(htmltools::span(style = sprintf("color:%s;", if (viol) AMB else DIM),
        sprintf("increments: \u03b1 = %.2f \u2192 r\u2081 should be \u2264 %+.3f, observed %+.3f",
                b$alpha_short %||% NA_real_, b$r1_fgn, b$r1_norm),
        if (viol) htmltools::HTML(sprintf(
          "  &#9888; %.1f\u00d7 too weak: the positions carry correlated error (galvo damping, or the tracking loop when motion is small)",
          b$r1_fgn / b$r1_norm)) else "  \u2713 consistent"))
    } else NULL

    # Motion blur used by every estimator (trk_blur(), est_blur.R), and where it came
    # from. Damping in the sequence smooths positions ACROSS steps, which no blur model
    # describes, so it is flagged rather than corrected.
    bl <- tryCatch(trk_blur(), error = function(e) NULL)
    blur_line <- if (!is.null(bl)) row(
      dim_(sprintf("motion blur R = %.3f per step, from %s", bl$R_step, bl$source)),
      if (!is.null(bl$note)) warn(sprintf("  (sequence file not used: %s)", bl$note)) else NULL,
      if (is.finite(bl$damping %||% NA_real_) && bl$damping > 0)
        warn(htmltools::HTML(sprintf(
          "  &#9888; sequence damping = %g: positions are smoothed across steps, so no blur model applies and D from single steps is approximate",
          bl$damping))) else NULL) else NULL
    # evidence that the photon prefactor is right; a large gap means one is wrong.
    # CVE sigma, kept only as a sanity read on the covariance. It is NOT a precision
    # measurement: ~90% of it is the assumed blur term, and on immobile samples the
    # increments come out positively correlated so it fails outright. The prefactor
    # it used to cross-check is now measured directly on a fixed sample, so the old
    # implied-L / R-range readout was retired.
    cve_line <- if (isTRUE(b$degenerate))
      row(dim_("CVE σ: "), warn("unresolvable (consecutive steps are positively correlated)"))
    else if (is.finite(b$sigma))
      row(dim_(sprintf("CVE σ %.2f nm (blur-aware, R = %.3f), %s",
                       b$sigma, b$R_blur %||% NA_real_,
                       if (is.finite(b$noise_frac) && b$noise_frac < 0.15)
                         "the covariance carries little of it" else "not the app-wide value")))
    else NULL

    # ── lateral vs axial (3D only) - the pooled sigma averages x, y and z, so when
    #    the axial precision differs from the lateral one it misstates the floor the
    #    x-y-only tools (D-map, 2D JDD/MSD) subtract. The axial spread is usually the
    #    TIGHTER one in MINFLUX, so the warning is two-sided: an axial-tighter file
    #    (ratio < 1/1.3) under-subtracts in those tools just as an axial-wider one
    #    over-subtracts.
    aniso_line <- if (isTRUE(b$has_z) && is.finite(b$sigma_xy) && is.finite(b$sigma_z))
      row(dim_(sprintf("lateral %.2f  \u00b7  axial %.2f nm", b$sigma_xy, b$sigma_z)),
          if (is.finite(b$sig_aniso) && (b$sig_aniso > 1.3 || b$sig_aniso < 1 / 1.3))
            warn(htmltools::HTML(sprintf(
              "  &#9888; axial is %.1f\u00d7 %s. Enter %.2f nm (lateral) for the D-map and 2D JDD.",
              if (b$sig_aniso > 1) b$sig_aniso else 1 / b$sig_aniso,
              if (b$sig_aniso > 1) "wider" else "tighter",
              b$sigma_xy))) else NULL)
    else NULL

    # ── how well sigma is resolved (r1 -> 0 means the per-track spread is motion) ──
    noise_line <- if (is.finite(b$noise_frac))
      row(dim_(sprintf("noise %.0f%% of step variance (r\u2081 = %.3f)", 100 * max(0, b$noise_frac), b$r1_norm)),
          if (b$noise_frac < 0.15) warn(htmltools::HTML("  &#9888; "), sigloc(),
                                          " spread is dominated by motion, so read it as indicative only") else NULL)
    else NULL

    # ── photon channel = mobility proxy, so an ECO/EFO filter is a mobility filter ──
    eco_line <- if (is.finite(b$phot_mob_rho) && b$phot_mob_rho <= -0.35)
      row(warn(htmltools::HTML(sprintf(
        "&#9888; %s correlates with mobility (&rho; = %.2f), so filtering on %s removes the faster molecules and lowers D",
        toupper(b$phot_col %||% "photon"), b$phot_mob_rho, toupper(b$phot_col %||% "photon")))))
    else NULL

    # ── data quality: same integrated-fraction metric as TSA / cycle-time status ──
    fi_raw  <- b$frac_integ %||% NA_real_;      fi_filt <- b$frac_integ_filt %||% NA_real_
    c_raw   <- b$cont_raw   %||% NA_real_;      c_filt  <- b$cont_filt      %||% NA_real_
    q_use   <- if (is.finite(fi_filt)) fi_filt else fi_raw
    q_verd  <- if (!is.finite(q_use)) c("UNKNOWN", DIM)
               else if (q_use >= 0.7) c("GOOD", "#1a8a1a")
               else if (q_use >= 0.5) c("MODERATE", "#b26a00")
               else c("POOR", RED)
    pair_txt <- function(a, z, fmt) if (is.finite(a) && is.finite(z)) sprintf(paste0(fmt, "\u2192", fmt), a, z)
                                    else if (is.finite(a)) sprintf(fmt, a) else "-"
    # D in the MOTION's degrees of freedom: a per-axis CVE spreads 2-DOF motion over
    # 3 recorded axes, so surface motion tracked in 3D reads 2/3 of the truth.
    d_obs <- if (isTRUE(b$has_z)) 3L else 2L
    d_mot <- { mo <- input$trk_motion %||% "volume"
               if (identical(mo, "filament")) 1L else if (identical(mo, "surface")) min(2L, d_obs) else d_obs }
    d_use <- b$D_cve * d_obs / d_mot
    qual_line <- row(
      htmltools::tags$b("Quality "),
      htmltools::span(style = sprintf("color:%s; font-weight:bold;", q_verd[2]), q_verd[1]),
      dim_(sprintf("  \u00b7  integrated %s  \u00b7  continuity %s%s",
                   pair_txt(100 * fi_raw, 100 * fi_filt, "%.0f%%"),
                   pair_txt(c_raw, c_filt, "%.2f"),
                   if (is.finite(c_raw) && is.finite(c_filt) && (c_filt - c_raw) > 0.01) " \u25b2" else "")))
    # Where a came from: a calibration profile, the ideal-estimator fallback
    # (uncalibrated), or a value typed by hand.
    prof_line <- if (is.finite(b$sig_phot_med %||% NA_real_)) {
      a_used <- b$a_pref %||% NA_real_
      ideal  <- .tcp_geom(b$tcp_pattern %||% TCP_PATTERN_DEFAULT)$a_ideal
      prof   <- trk_cal_profiles()[[input$trk_cal_profile %||% ""]]
      if (!is.null(prof) && isTRUE(abs(a_used - prof$a) < 1e-3))
        row(dim_(sprintf("calibrated: %s  \u00b7  a = %.2f (quartiles %.2f\u2013%.2f over %d cells)",
                         prof$name, prof$a, prof$a_q25 %||% NA_real_, prof$a_q75 %||% NA_real_,
                         as.integer(prof$n_cells %||% NA))))
      else if (isTRUE(abs(a_used - ideal) < 1e-9))
        row(warn(sprintf(paste0("uncalibrated: a = %.2f is the ideal-estimator value for this pattern; ",
                                "measure it for this sequence under Instrument Calibration"), a_used)))
      else
        row(dim_(sprintf("a = %.2f entered by hand  \u00b7  at the ideal-estimator value %.2f, %s would be %.2f nm",
                         a_used, ideal, "\u03c3", b$sig_phot_med * a_used / ideal)))
    } else NULL

    calib_line <- row(dim_(sprintf(
      "cycle %s  \u00b7  D %s \u00b5m\u00b2/s%s  \u00b7  %s steps",
      num(b$t_cycle, "%.4g ms", "unknown"), num(d_use / 1e3, "%.3g"),
      if (d_mot != d_obs) sprintf(" (%d/%d for %s)", d_obs, d_mot, input$trk_motion %||% "volume") else "",
      format(b$n_steps, big.mark = ","))))

    # Tracker wander floor: the apparent diffusion the loop produces on an immobile
    # target. Reported beside D, never subtracted silently - D_w scales with the
    # photon budget, so a single constant would over-correct dim tracks.
    dw_line <- {
      dw <- suppressWarnings(as.numeric(input$trk_dw_floor %||% 0))
      if (!is.finite(dw) || dw <= 0 || !is.finite(d_use)) NULL else {
        frac <- dw / d_use
        row(dim_(sprintf("wander floor %.3g µm²/s  →  D − Dw = %.3g",
                         dw / 1e3, max(0, d_use - dw) / 1e3)),
            if (frac > 0.5) warn(htmltools::HTML(sprintf(
              "  &#9888; floor is %.0f%% of D (not resolved above the tracker)", 100 * frac)))
            else if (frac > 0.2) warn(sprintf("  (%.0f%% of D)", 100 * frac))
            else dim_(sprintf("  (%.0f%% of D)", 100 * frac)))
      }
    }

    # ── linking QC: implausible jumps and per-track sigma outliers ──────────────
    qc_line <- if (is.finite(b$frac_implaus %||% NA_real_) || is.finite(b$n_suspect %||% NA_real_)) {
      bad <- (is.finite(b$frac_implaus) && b$frac_implaus > 0.01) ||
             (is.finite(b$n_suspect) && is.finite(b$n_tracks_pt) && b$n_tracks_pt > 0 &&
              b$n_suspect / b$n_tracks_pt > 0.05)
      txt <- paste(c(
        if (is.finite(b$n_implaus)) sprintf("%s implausible jumps (%.2f%%)", format(b$n_implaus, big.mark = ","), 100 * b$frac_implaus),
        if (is.finite(b$n_suspect)) sprintf("%d/%d suspicious tracks", b$n_suspect, b$n_tracks_pt)), collapse = "  \u00b7  ")
      row(htmltools::span(style = sprintf("color:%s;", if (bad) AMB else DIM),
        htmltools::tags$b("Linking QC "), txt, if (bad) ". Check for mis-links at high density or across long gaps." else ""))
    } else NULL

    htmltools::div(style = "font-size:0.92em; margin:4px 0; line-height:1.45;",
      sigma_line, prof_line, smooth_line, blur_line, cve_line, aniso_line, noise_line, eco_line, qual_line, calib_line,
      dw_line, qc_line)
  })

  # Per-track σ_loc heterogeneity histogram (σ_loc is NOT a single instrument
  # constant - this shows the spread; a tight peak justifies one global value).
  # Distribution of the sigma_loc the app actually uses: one value per localization,
  # sigma_i = L/(a*sqrt(N_i)). Plotting this rather than the per-track CVE sigma keeps
  # the curve and the headline number the same quantity - the dashed median IS the
  # figure reported above, so the two can never disagree.
  output$trk_sigma_pertrack_plot <- renderPlotly({
    b <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    sg <- if (!is.null(b)) b$sig_phot else NULL
    sg <- sg[is.finite(sg) & sg > 0]
    shiny::validate(shiny::need(length(sg) >= 20,
      "Needs an ECO column: sigma_loc is computed per localization from the photon count."))
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    med <- stats::median(sg); q <- stats::quantile(sg, c(.25, .75))
    plot_ly() %>%
      add_histogram(x = sg, nbinsx = 60,
                    marker = list(color = "rgba(120,160,220,0.7)",
                                  line = list(color = "rgba(80,120,180,0.9)", width = 0.4)),
                    hovertemplate = "%{x:.2f} nm<br>%{y} localizations<extra></extra>",
                    name = "sigma") %>%
      layout(
        xaxis = list(title = "σ<sub>loc</sub> per localization (nm)", color = fg, zeroline = FALSE),
        yaxis = list(title = "localizations", color = fg),
        shapes = list(
          list(type = "rect", x0 = q[1], x1 = q[2], y0 = 0, y1 = 1, yref = "paper",
               fillcolor = "rgba(120,160,220,0.13)", line = list(width = 0), layer = "below"),
          list(type = "line", x0 = med, x1 = med, y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#c03030", width = 2, dash = "dash"))),
        annotations = list(list(
          x = med, y = 1, yref = "paper", yanchor = "bottom", showarrow = FALSE,
          text = sprintf("median %.2f nm  ·  IQR %.2f–%.2f", med, q[1], q[2]),
          font = list(color = "#c03030", size = 11))),
        paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
        showlegend = FALSE, bargap = 0.04, uirevision = "sigma_pl")
  })

  # Transparency funnel: how many localizations / tracks survive each pipeline
  # stage. Analyzing well-tracked molecules is a legitimate, DISCLOSED choice -
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
    fmt  <- function(v) if (is.finite(v)) format(round(v), big.mark = ",") else "-"
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
        "Analyses run on the last row. Filtering for well-tracked molecules is a disclosed choice - ",
        "report this funnel in your methods. Track count can rise at Gap-split (segments) then fall."))
  })

  # GOAL: show whether track continuity follows mobility in the analysed data, i.e.
  # whether keeping "well-tracked" molecules means keeping slow ones.
  # A molecule that moves fast relative to the TCP is lost for a cycle more often, so
  # its track collects dark cycles; a continuity threshold can then remove fast
  # molecules under the name of quality. The table ranks the analysed tracks by the
  # same continuity score the gate uses (same half-life) and gives each quarter its
  # own pooled CVE D. It is computed after "Estimate σ_loc & quality", on the data
  # the analyses see, so it describes the population actually being reported.
  # The ±20% band for the verdict is a reading aid, not a significance test.
  output$trk_continuity_mobility <- renderUI({
    b <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    if (is.null(b))
      return(htmltools::tags$p(style = "font-size:0.83em; color:#888;",
        "Run Estimate σ_loc & quality first."))
    df  <- tryCatch(trk(), error = function(e) NULL)
    t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    H   <- suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3))
    sc  <- .continuity_scores(df, t_c, H)
    tab <- continuity_mobility_table(df, sc, t_c, has_z = .has_real_z(df))
    if (is.null(tab))
      return(htmltools::tags$p(style = "font-size:0.83em; color:#888;",
        "Too few tracks to split into quarters (at least 40 needed)."))

    # Two readings, kept apart. r_mol: the gradient across molecules (best quarter vs
    # the other three), which is what a per-track statistic sees. r_pool: what keeping
    # only the best quarter would do to the POOLED D, which weights tracks by their
    # steps; when that quarter holds most localizations, r_mol can be large while the
    # reported number barely moves.
    r_mol  <- attr(tab, "ratio_best_rest")
    r_pool <- attr(tab, "ratio_best_all")
    share4 <- tab$loc_share[4]
    coupled <- is.finite(r_mol) && abs(log(r_mol)) > log(1.2)
    verdict <- if (!is.finite(r_mol)) NULL
      else if (!coupled)
        htmltools::tags$p(style = "font-size:0.85em; color:#1a8a1a; margin:4px 0;",
          sprintf("D is within 20%% across continuity (most continuous quarter %.2f\u00d7 the rest): the continuity gate is not selecting on speed.", r_mol))
      else
        htmltools::tags$div(style = "font-size:0.85em; margin:4px 0;",
          htmltools::tags$p(style = "color:#B8860B; margin:0 0 3px 0;",
            htmltools::HTML(sprintf(
              "&#9888; <b>Across molecules</b>, the most continuous quarter is %.0f%% %s than the rest (D %.3g vs %.3g \u00b5m\u00b2/s): per-track statistics depend on how strictly tracks are gated.",
              100 * abs(r_mol - 1), if (r_mol < 1) "slower" else "faster",
              tab$D_um2_s[4], attr(tab, "D_rest_um2_s")))),
          htmltools::tags$p(style = "color:#555; margin:0;",
            htmltools::HTML(sprintf(
              "<b>On the pooled D</b>, that quarter holds %.0f%% of the localizations, so keeping only it moves D from %.3g to %.3g \u00b5m\u00b2/s (%+.0f%%).",
              100 * share4, attr(tab, "D_all_um2_s"), tab$D_um2_s[4], 100 * (r_pool - 1)))))

    fmt <- function(v, d) ifelse(is.finite(v), formatC(v, format = "f", digits = d), "-")
    rows <- lapply(seq_len(nrow(tab)), function(i) htmltools::tags$tr(
      lapply(c(tab$quartile[i], tab$n_tracks[i],
               sprintf("%s-%s", fmt(tab$score_lo[i], 2), fmt(tab$score_hi[i], 2)),
               paste0(fmt(100 * tab$loc_share[i], 1), "%"),
               paste0(fmt(100 * tab$integ_frac[i], 0), "%"),
               fmt(tab$D_um2_s[i], 3), fmt(tab$sigma_nm[i], 2)),
             function(v) htmltools::tags$td(v, style = "padding:1px 8px; text-align:right;"))))
    htmltools::tagList(
      verdict,
      htmltools::tags$table(style = "font-size:0.83em; border-collapse:collapse;",
        htmltools::tags$thead(htmltools::tags$tr(lapply(
          c("Continuity", "tracks", "score", "of locs", "integrated steps", "D (µm²/s)", "σ (nm)"),
          function(h) htmltools::tags$th(h, style = "padding:1px 8px; text-align:right; border-bottom:1px solid #ccc; color:#666;")))),
        htmltools::tags$tbody(rows)),
      htmltools::tags$p(style = "font-size:0.8em; color:#888; margin:5px 0 0 0;",
        sprintf("Four equal groups by continuity score (half-life %.3g cycles); D is the pooled CVE of each.", H)))
  })


  # Non-reactive loop state for sequential mode (safe for later())
  is_seq   <- reactive(identical(input$trk_play_mode, "seq"))
  is_simul <- reactive(identical(input$trk_play_mode, "simul"))

  # Bin width (ms) is the master time-resolution control. The frame count is
  # derived from it (window / bin width) rather than entered separately - the
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
    req(trk_has_data())
    
    # If we're still preparing chunks, show a placeholder message
    if (isTRUE(rv$chunking)) {
      return(tags$div(
        "Preparing chunks... This can take a bit for very large files.",
        style = "color:#555;"
      ))
    }
    
    # Otherwise, if a manifest exists, show the chunk selector
    req(!is.null(rv$manifest_dir))
    dirp <- rv$manifest_dir
    manf <- file.path(dirp, "manifest.csv")
    if (file.exists(manf)) {
      man <- readr::read_csv(manf, show_col_types = FALSE)
      chunk_choices <- setNames(man$file, paste0("Chunk ", man$chunk, " (", man$n_tids, " tids)"))
      all_choices <- c("All TIDs" = "__all__", chunk_choices)
      # Chunk 1 if the user was scoped to a chunk before this rebuild, else all TIDs.
      sel_now <- if (identical(trk_chunk_pref(), "first") && nrow(man) >= 1L)
                   man$file[1] else "__all__"
      return(tagList(
        selectInput("trk_chunk_sel", "Playback chunk", choices = all_choices, selected = sel_now),
        helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
          "Filtering, gap-splitting, analysis and the full-track overlay always use ALL ",
          "TIDs. This selector only limits which TIDs are ANIMATED during playback, to ",
          "keep the animation smooth - \"All TIDs\" animates everything.")
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
  # 2D/3D toggle and the XY/XZ/YZ plane selector are static inputs (ui_tracking.R).
  # They used to be renderUIs that VANISHED on planar data. Two reasons that is
  # worse than greying them:
  #   - a control that disappears reads as a bug, not as "not applicable";
  #   - the plane selector keyed off trk(), which only exists after Plot Data, so
  #     it was absent on upload and appeared later.
  # The observer below disables the z-dependent choices instead, using the SAME
  # planar test as the axial calibration panel.
  observe({
    planar <- identical(trk_z_state(), "planar")

    shinyjs::toggleState(selector = "input[name='trk_view_mode'][value='3d']",  condition = !planar)
    shinyjs::toggleState(selector = "input[name='trk_view_plane'][value='XZ']", condition = !planar)
    shinyjs::toggleState(selector = "input[name='trk_view_plane'][value='YZ']", condition = !planar)

    # Disabling a radio does not deselect it - force the view back, or a planar
    # upload while in 3D/XZ leaves the plot in a mode whose control is greyed out.
    if (planar) {
      if (!identical(isolate(input$trk_view_mode) %||% "2d", "2d"))
        updateRadioButtons(session, "trk_view_mode", selected = "2d")
      if (!identical(isolate(input$trk_view_plane) %||% "XY", "XY"))
        updateRadioButtons(session, "trk_view_plane", selected = "XY")
    }
  })

  output$trk_z_axes_note <- renderUI({
    if (!identical(trk_z_state(), "planar")) return(NULL)
    helpText(tags$b("2D dataset: "), "3D mode and the XZ / YZ side views are unavailable.")
  })

  # Switching view mode requires a full plot rebuild
  observeEvent(input$trk_view_mode, {
    # Compare against what was last RENDERED, so a programmatic reset to "2d"
    # (planar guard above) or a value that already matches the plot does not
    # trigger a redundant rebuild.
    if (identical(input$trk_view_mode %||% "2d", loop$rendered_view %||% "2d")) return()
    fire_layout("view_mode_2d3d")
  }, ignoreInit = TRUE)
  # A plane change is a coordinate remap baked in frames_store, so it needs a rebuild -
  # same creation-fire guard as the 2D/3D toggle (compare to the plane last rendered).
  observeEvent(input$trk_view_plane, {
    if (identical(input$trk_view_plane %||% "XY", loop$rendered_plane %||% "XY")) return()
    fire_layout("view_plane")
  }, ignoreInit = TRUE)
  # ── Axial calibration: pending vs applied ────────────────────────────────
  # input$trk_z_scale is only a PENDING value; trk_z_applied() is what the data
  # is actually scaled by. Committing on Plot Data matches how every other
  # tracking setting behaves (see the note above observeEvent(input$trk_plot_data))
  # and avoids re-running the whole pipeline on every keystroke - the reason the
  # old control needed a debounce.
  trk_z_applied <- reactiveVal(TRK_Z_SCALE_DEFAULT)

  observeEvent(input$trk_plot_data, {
    val <- suppressWarnings(as.numeric(input$trk_z_scale))
    if (!length(val) || !is.finite(val) || val <= 0) return()
    if (isTRUE(all.equal(val, trk_z_applied()))) return()
    trk_z_applied(val)
    # Record it either way, but only announce it when there is axial data for it
    # to act on - a 2D dataset should never surface a z message.
    if (identical(trk_z_state(), "3d"))
      showNotification(sprintf("Axial calibration applied: %.3f.", val),
                       type = "message", duration = 5)
    # No fire_layout here on purpose: the generic Plot Data handler below already
    # calls it, and priority = 10 guarantees this runs first - so that single
    # rebuild picks up the new factor. Firing here too would rebuild twice.
  }, ignoreInit = TRUE, priority = 10)

  # Three states, because "nothing uploaded yet" must not look like "this data is
  # 2D". Reads the RAW loaded data, not trk(), so the panel settles as soon as a
  # file is uploaded rather than waiting for the first Plot Data.
  trk_z_state <- reactive({
    df <- tryCatch(trk_raw_loaded(), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return("nodata")
    if (isTRUE(.has_real_z(df))) "3d" else "planar"
  })

  # Dim (not hide) for planar data - see the SMLM equivalent for the reasoning.
  observe({
    ok <- !identical(trk_z_state(), "planar")
    shinyjs::toggleState("trk_z_scale", condition = ok)
    shinyjs::toggleClass(id = "trk_z_scale_group", class = "dimmed", condition = !ok)
  })

  output$trk_z_scale_status <- renderUI({
    if (identical(trk_z_state(), "planar"))
      return(helpText(tags$b("Not applicable: "),
                      "the uploaded dataset has no axial (z) data."))

    pend <- suppressWarnings(as.numeric(input$trk_z_scale))
    appl <- trk_z_applied()
    if (!length(pend) || !is.finite(pend) || pend <= 0)
      return(helpText(sprintf("Applied: %.3f. Enter a positive factor.", appl)))
    if (isTRUE(all.equal(pend, appl)))
      helpText(sprintf("Applied: %.3f.", appl))
    else
      helpText(tags$b(sprintf("Pending: %.3f", pend)),
               sprintf(" (applied: %.3f) - press Plot Data to apply.", appl))
  })

  # ============================================================================
  # DATASETS - uploads add files to a store; one VIEW of it is active
  # ============================================================================
  # Every uploaded file is read once and kept on disk (trk_store). The active
  # view is either one file or all loaded files merged (track ids then prefixed
  # F1_, F2_, ... so they cannot collide). Activating a view runs the ingest
  # below - chunking, frozen raw statistics, the Filter & Rechunk reset point -
  # exactly as an upload used to, so everything downstream sees a fresh dataset.
  # Removing a file drops it from the store and activates another one.
  trk_store     <- reactiveVal(list())          # name -> list(path, n, json_itr)
  trk_view      <- reactiveVal(character(0))    # names of the files in the active view
  trk_store_dir <- file.path(tempdir(), paste0("trk_store_", session$token))
  .trk_last_active <- ""                        # last ingested view (guards the selector echo)
  .trk_ingest_n    <- 0L                        # makes every ingest's chunk dir unique

  # TRUE once there is a dataset to show (any view, or an adopted chunk manifest).
  trk_has_data <- reactive(length(trk_view()) > 0L || !is.null(rv$manifest_dir))

  # Name of the active view for labels and download file names.
  trk_view_label <- reactive({
    v <- trk_view()
    if (!length(v)) "tracking" else if (length(v) == 1L) v
    else sprintf("%d files merged (%s, …)", length(v), v[1])
  })
  trk_view_base <- reactive({
    v <- trk_view()
    if (!length(v)) return("tracking")
    b <- tools::file_path_sans_ext(v[1])
    if (length(v) > 1L) paste0(b, "_merged", length(v)) else b
  })

  # Loader for one file → data.frame (or NULL).  Returns attr "json_itr".
  load_one_trk <- function(p, nm) {
    d <- NULL; jitr <- NA_integer_
    if (grepl("\\.csv(\\.gz)?$", nm)) {
      d <- data.table::fread(p, showProgress = FALSE, nThread = max(1, parallel::detectCores()-1)) |>
           as.data.frame()
      # A CSV arrives with some axial factor already applied; return z to raw so
      # trk_raw_df() applies the user's calibration once (it was applied twice).
      d <- strip_csv_z_calibration(d)
      d <- canonicalize_minflux_df(d)
    } else if (grepl("\\.(qs|fst)$", nm)) {
      # read_chunk_file reports a missing 'qs' rather than throwing, so a machine
      # without it degrades to "this file was skipped" instead of a stack trace.
      d <- read_chunk_file(p)
      if (is.null(d)) return(NULL)
    } else if (grepl("\\.rds$", nm)) {
      d <- readRDS(p)
      if (!is.data.frame(d)) return(NULL)
      d <- canonicalize_minflux_df(d)
    } else if (grepl("\\.json$", nm)) {
      # z_scale = 1: chunks store RAW z. The axial calibration is applied once,
      # in trk_raw_df(), so every format lands on the same footing.
      d <- read_minflux_json_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE,
                                   z_scale = 1.0)
      jitr <- if ("itr" %in% names(d) && nrow(d) > 0) max(d$itr, na.rm = TRUE) else NA_integer_
    } else if (grepl("\\.npy$", nm)) {
      d <- read_minflux_npy_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE,
                                  z_scale = 1.0)
      jitr <- if ("itr" %in% names(d) && nrow(d) > 0) max(d$itr, na.rm = TRUE) else NA_integer_
    } else {
      return(NULL)
    }
    attr(d, "json_itr") <- jitr
    d
  }

  # Choices for the Active dataset selector: every file, plus the merged view.
  .trk_view_choices <- function(store) {
    nm <- names(store)
    ch <- stats::setNames(nm, nm)
    if (length(nm) >= 2L)
      ch <- c(ch, stats::setNames("__merged__", sprintf("All %d loaded files, merged", length(nm))))
    ch
  }

  # ── .trk_ingest() ──────────────────────────────────────────────────────────
  # GOAL: make one data.frame the active dataset: per-track summary, frozen raw
  # statistics (histogram backgrounds, filter reset), chunks on disk, manifest.
  # This is the body the upload observer used to run on its merged data.
  .trk_ingest <- function(df, label) {
    rv$bigfile_notice <- FALSE
    rv$chunking <- TRUE
    on.exit({ rv$chunking <- FALSE }, add = TRUE)
    .trk_ingest_n <<- .trk_ingest_n + 1L
    safe_base <- gsub("[^A-Za-z0-9_-]+", "_", label)
    out_dir <- file.path(tempdir(), sprintf("chunks_%s_%d_%d", safe_base, as.integer(Sys.time()), .trk_ingest_n))

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
    # Done once per dataset - never overwritten by chunk switching or rechunking.
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
    trk_spatial_filter(NULL)    # a crop drawn on another dataset does not apply here
    trk_chunk_pref("__all__")   # a new dataset does not inherit a chunk scoping

    user_tpc <- .tids_per_chunk()

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
      rv$chunking <- FALSE
    } else {
      # Fallback: compress to a single chunk, in whatever format the app writes.
      # This used to call qs::qsave() unconditionally, which crashed the upload on
      # any machine without qs even though CHUNK_FORMAT had already fallen back.
      out_one <- file.path(out_dir, paste0(safe_base, ".", CHUNK_FORMAT))
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      if (identical(CHUNK_FORMAT, "qs") && requireNamespace("qs", quietly = TRUE))
        qs::qsave(df, out_one, preset = "high")
      else
        fst::write_fst(as.data.frame(df), out_one, compress = 50)
      man <- data.frame(chunk = 1L, n_tids = length(unique(df$tid)), file = out_one)
      readr::write_csv(man, file.path(out_dir, "manifest.csv"))
      readr::write_csv(tid_sum, file.path(out_dir, "tid_summary.csv"))
      rv$manifest_dir          <- out_dir
      rv$original_manifest_dir <- out_dir
      trk_data_key(as.character(as.numeric(Sys.time())))
      rv$chunking <- FALSE
    }
    invisible(NULL)
  }

  # ── .trk_activate() ────────────────────────────────────────────────────────
  # GOAL: build the data of one view from the store and ingest it.
  # sel : a file name, or "__merged__" for all loaded files
  .trk_activate <- function(sel) {
    store <- isolate(trk_store())
    files <- if (identical(sel, "__merged__")) names(store) else intersect(sel, names(store))
    if (!length(files)) return(invisible(NULL))
    dfs <- lapply(seq_along(files), function(i) {
      d <- read_chunk_file(store[[files[i]]]$path)
      if (is.null(d)) return(NULL)
      if (length(files) > 1L) d$tid <- paste0("F", i, "_", as.character(d$tid))
      d
    })
    df <- dplyr::bind_rows(Filter(Negate(is.null), dfs))
    if (!nrow(df)) {
      showNotification("The selected dataset could not be read.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    .trk_last_active <<- sel
    trk_view(files)
    withProgress(message = sprintf("Preparing %s...", if (length(files) > 1L)
                                     sprintf("%d merged files", length(files)) else files),
                 value = 0.5, .trk_ingest(df, isolate(trk_view_base())))
    if (length(files) > 1L)
      showNotification(sprintf("Merged %d files into one dataset (%s localizations).",
                               length(files), format(nrow(df), big.mark = ",")),
                       type = "message", duration = 6)
    invisible(NULL)
  }

  # Point the selector at a view and ingest it directly. Setting the selector to
  # the value it already holds does not fire its observer (e.g. re-uploading a
  # file of the same name), so the ingest cannot be left to the observer.
  .trk_set_active <- function(sel) {
    updateSelectInput(session, "trk_active", choices = .trk_view_choices(isolate(trk_store())),
                      selected = sel)
    .trk_activate(sel)
  }

  # ---- Upload: read every file into the store -------------------------------
  observeEvent(input$trk_file, {
    req(input$trk_file)
    n_files   <- nrow(input$trk_file)
    paths     <- input$trk_file$datapath
    names_o   <- as.character(input$trk_file$name)          # as uploaded, for display
    names_v   <- tolower(names_o)                           # for extension checks
    size_bytes <- sum(suppressWarnings(as.numeric(input$trk_file$size)), na.rm = TRUE)
    if (!is.finite(size_bytes) || size_bytes <= 0) {
      info <- tryCatch(file.info(paths), error = function(e) NULL)
      if (!is.null(info)) size_bytes <- sum(info$size, na.rm = TRUE)
    }
    size_mb <- round(as.numeric(size_bytes) / (1024^2), 1)

    # If a single file is uploaded next to an existing manifest, adopt it and bail.
    # (Skipped for multi-file uploads.) The adopted chunks bypass the store.
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

    dir.create(trk_store_dir, showWarnings = FALSE, recursive = TRUE)
    store <- trk_store(); added <- character(0); json_itr <- NA_integer_
    keep_all <- c("tid","x","y","z","tim","itr","efo","cfr","dcr","eco","dwell","fluo")
    withProgress(message = sprintf("Loading %d file(s) (%.1f MB)...", n_files, size_mb), value = 0, {
      for (fi in seq_len(n_files)) {
        incProgress(1 / n_files, detail = names_o[fi])
        di <- tryCatch(load_one_trk(paths[fi], names_v[fi]), error = function(e) NULL)
        if (is.null(di) || !is.data.frame(di) || !nrow(di)) {
          showNotification(sprintf("%s could not be read and was skipped.", names_o[fi]),
                           type = "warning", duration = 8)
          next
        }
        if (!all(c("tid","x","y","tim") %in% names(di))) {
          showNotification(sprintf("%s is missing required columns (tid, x, y, tim) and was skipped.",
                                   names_o[fi]), type = "error", duration = 8)
          next
        }
        ji <- attr(di, "json_itr"); if (!is.na(ji)) json_itr <- ji
        di$tid <- as.character(di$tid)
        # Keep only the columns tracking needs  -  reduces stored file size significantly.
        # No automatic CFR pre-filter: CFR is not a strong artefact indicator for
        # tracking, and the metadata filter offers a CFR range when it is wanted.
        di <- di[, intersect(keep_all, names(di)), drop = FALSE]
        old <- store[[names_o[fi]]]
        if (!is.null(old)) unlink(old$path)                    # same name: replaced
        f <- file.path(trk_store_dir, sprintf("f%03d_%s.fst", length(store) + length(added) + 1L,
                                              gsub("[^A-Za-z0-9_-]+", "_", tools::file_path_sans_ext(names_o[fi]))))
        fst::write_fst(di, f, compress = 50)
        store[[names_o[fi]]] <- list(path = f, n = nrow(di), json_itr = ji)
        added <- c(added, names_o[fi])
      }
    })
    if (!length(added)) return(invisible(NULL))
    trk_store(store)

    # Show JSON iteration notification (mirrors SMLM tab behaviour)
    if (!is.na(json_itr))
      showNotification(sprintf(".json or .npy import: using final localization (itr = %d) for x/y.",
                               as.integer(json_itr)), type = "message", duration = 8)
    if (length(store) > 1L)
      showNotification(sprintf("%d datasets loaded. Showing %s; switch or merge them under Active dataset.",
                               length(store), added[1]), type = "message", duration = 7)
    .trk_set_active(added[1])                                 # one file at a time by default
  }, ignoreInit = TRUE)

  # User picks another view.
  observeEvent(input$trk_active, {
    sel <- input$trk_active
    req(nzchar(sel %||% ""))
    if (identical(sel, .trk_last_active)) return(invisible(NULL))   # already active
    .trk_activate(sel)
  }, ignoreInit = TRUE)

  # Remove the selected file from the store and show another one.
  observeEvent(input$trk_remove, {
    sel <- input$trk_active %||% ""
    store <- trk_store()
    if (identical(sel, "__merged__")) {
      showNotification("Select a single dataset to remove it.", type = "warning", duration = 5)
      return(invisible(NULL))
    }
    if (!nzchar(sel) || is.null(store[[sel]])) return(invisible(NULL))
    unlink(store[[sel]]$path)
    store[[sel]] <- NULL
    trk_store(store)
    if (length(store)) {
      .trk_set_active(names(store)[1])
    } else {
      .trk_last_active <<- ""
      trk_view(character(0))
      rv$manifest_dir <- NULL; rv$original_manifest_dir <- NULL
      trk_raw_stats_rv(NULL)
      updateSelectInput(session, "trk_active", choices = character(0))
      trk_data_key(as.character(as.numeric(Sys.time())))
    }
    showNotification(sprintf("%s removed.", sel), type = "message", duration = 4)
  })

  # ── Adopt an axial factor recorded in an uploaded CSV ─────────────────────
  # Reads a small head of each file rather than the whole thing, so this stays
  # cheap even for a multi-GB export. A sample rather than a single row because
  # z VARIATION is what decides whether the calibration applies at all, and one
  # row cannot show variation - a planar file must be passed over in silence.
  #
  # Done here rather than inside trk_raw_loaded() because that is a reactive:
  # writing to an input from inside one is a side effect in the wrong place.
  #
  # Only `certain` factors are adopted - i.e. a real z_scale_applied column. The
  # 0.7 assumed for pyMINFLUX and legacy exports is a guess, and a guess must not
  # silently become the session's calibration.
  observeEvent(input$trk_file, {
    req(input$trk_file)
    adopt <- NULL
    for (i in seq_len(nrow(input$trk_file))) {
      if (!grepl("\\.csv(\\.gz)?$", tolower(input$trk_file$name[i]))) next
      head_df <- tryCatch(
        as.data.frame(data.table::fread(input$trk_file$datapath[i], nrows = 2000L,
                                        showProgress = FALSE)),
        error = function(e) NULL)
      if (is.null(head_df) || !nrow(head_df)) next
      if (!has_axial_data(head_df, min_sd = 1)) next   # planar: nothing to adopt
      cal <- detect_csv_z_calibration(head_df)
      if (isTRUE(cal$certain) && is.finite(cal$factor) && cal$factor > 0) adopt <- cal$factor
    }
    if (is.null(adopt) || isTRUE(all.equal(adopt, trk_z_applied()))) return()

    updateNumericInput(session, "trk_z_scale", value = adopt)
    # Committed straight away: an upload is itself an explicit action, and the
    # pipeline is being rebuilt anyway, so there is nothing to defer to Plot Data.
    trk_z_applied(adopt)
    showNotification(sprintf("Axial calibration %.3f adopted from the uploaded file.", adopt),
                     type = "message", duration = 6)
  }, ignoreInit = TRUE, priority = 20)   # before the loaders read the data
  
  # ---- Load CSV/fast formats OR a pre-split chunk ----
  trk_raw_loaded <- reactive({
    req(trk_has_data())
    required <- c("tid","x","y","tim")
    
    # Source for the WHOLE filter pipeline. The filters and gap-splitting must
    # operate on ALL TIDs, never a single chunk, so we always load every chunk of
    # the CURRENT manifest (rv$manifest_dir) - the chunk selector does not restrict
    # the data. After Filter & Rechunk the current manifest is the rechunked set;
    # use Reset Ranges (which restores rv$original_manifest_dir) to return to the
    # full original, un-split data before re-splitting.
    df <- NULL
    if (!is.null(rv$manifest_dir)) {
      manf <- file.path(rv$manifest_dir, "manifest.csv")
      if (file.exists(manf)) {
        man <- tryCatch(readr::read_csv(manf, show_col_types = FALSE), error = function(e) NULL)
        if (!is.null(man) && nrow(man) > 0 && "file" %in% names(man)) {
          dfs <- lapply(man$file, read_chunk_file)
          df  <- dplyr::bind_rows(Filter(Negate(is.null), dfs))
        }
      }
    }
    # Every dataset is ingested into chunks (.trk_ingest), so there is always a
    # manifest once anything is loaded.
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No dataset loaded."))
    
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

  # ── Physical Z correction factor (applied at the pipeline ROOT) ──────────────
  # A refractive-index / mounting-media mismatch scales the apparent axial (z)
  # coordinate by a constant. Applying it in a thin layer directly over the chunk
  # loader means EVERY downstream reactive - track cleaning, gap-splitting, the CVE
  # σ_loc estimator, trk(), and the 2D/3D/side-view display - reads corrected z from
  # this single point, so the correction reaches all z-dependent parameters at once.
  #
  # Kept SEPARATE from trk_raw_loaded() so changing the factor re-runs only this
  # cheap rescale, never the expensive chunk reload. The z dependency is intentionally
  # NOT isolated: it must invalidate the whole filter pipeline so the cached stages
  # (which trk() reads via isolate()) recompute with the new z rather than serving a
  # stale cached value. z is not a filter/histogram variable, so rescaling it here has
  # no effect on the metadata panels or the raw-stats ranges.
  trk_raw_df <- reactive({
    df <- trk_raw_loaded()
    # trk_raw_loaded() is RAW, uncorrected z for every format (readers use
    # z_scale = 1; CSVs are stripped on load), so this is the ONE place the
    # axial calibration enters the tracking pipeline.
    zf <- trk_z_applied()
    if (is.finite(zf) && zf > 0 && zf != 1.0 && "z" %in% names(df) && any(is.finite(df$z)))
      df$z <- df$z * zf
    df
  })

  # TRUE, frozen raw - the ORIGINAL uploaded data, NEVER touched by the filters or by
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
    df <- dplyr::bind_rows(Filter(Negate(is.null), lapply(man$file, read_chunk_file)))
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
    req(trk_has_data())
    nm <- trk_view_label()
    # Chunks scope the playback animation only, so the label mentions them only in
    # Playback mode - in Full tracks the whole dataset is always shown.
    chunk_info <- if (!identical(input$trk_view_mode_main %||% "full", "playback")) "" else {
      sel <- input$trk_chunk_sel
      if (identical(sel, "__all__"))        " - All chunks"
      else if (isTruthy(sel))              paste0(" - ", basename(sel))
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
      # Value: the stored filter, clamped to this dataset's range (full range if none).
      val <- .filt_value(col, bounds)
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
  # (post-split) population, but update at PLOT cadence only - keyed on data_ver()
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
    val <- .filt_value("tid_len", bounds)      # stored filter, clamped (full range if none)
    isolate(trk_len_reset_flag(FALSE))
    numericRangeInput("trk_tid_len", "Trace ID length (localizations)",
                      value = val, min = bounds[1], max = bounds[2])
  })

  # Histogram zoom range - always anchored to the ORIGINAL upload stats so the
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
  # All filter UIs (metadata, tid_len, tid_dur) are renderUI-wrapped - incrementing
  # reset_ranges_trigger forces them all to re-render with value = full rng from
  # trk_raw_stats_rv().  updateNumericRangeInput does not reliably update
  # renderUI-wrapped widgets and is intentionally not used here.
  reset_tracking_filters <- function() {
    if (!is.null(rv$original_manifest_dir))
      rv$manifest_dir <- rv$original_manifest_dir
    trk_rechunk_stats_rv(NULL)  # clear rechunk stats so UIs revert to original upload ranges
    trk_rechunk_bbox(NULL)      # clear locked extent so axes revert to per-TID bbox
    trk_spatial_filter(NULL)    # clear any active spatial crop
    for (k in names(.FILT_IDS)) trk_filt[[k]] <- NULL   # no range filters
    trk_len_reset_flag(TRUE)    # snap length + duration sliders to full range
    reset_dur_filter_flag(TRUE)
    reset_ranges_trigger(isolate(reset_ranges_trigger()) + 1L)
  }

  observeEvent(input$trk_reset_filters, {
    reset_tracking_filters()
  }, ignoreInit = TRUE)


  # Legacy re-render flags for the length / duration widgets. The filters themselves
  # live in trk_filt (FILTER SETTINGS above), so a new dataset can no longer receive
  # a stale widget value, and these no longer change what is filtered.
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
  # Pipeline: Stage 1 (metadata) → Stage 2 (trace ID length) → Stage 3 (duration).
  # Each stage uses metrics computed from the PREVIOUS stage's output so that
  # upstream filters propagate correctly (e.g. metadata filtering shrinks trace
  # localization counts, which then drives the length filter).
  # All filter inputs use NULL pass-through (input not yet rendered → all data
  # passes) and a stale-range guard (filter range outside data range → all data
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
    filter_ranges(df, trk_filt_ranges())                    # R/est_clean.R
  })

  # ===========================================================================
  #  FILTER SETTINGS - stored as the user set them, independent of any dataset
  # ===========================================================================
  # A range filter ("EFO between 50k and 200k") applies to every loaded dataset
  # alike, so it is kept here as the user entered it. The pipeline reads these
  # values; the widgets only DISPLAY them, clamped to the active dataset's range.
  # Storing the clamped widget values instead would silently change the filter
  # every time a dataset with a narrower range was shown. Reset Ranges clears them.
  .FILT_IDS <- c(efo = "efo_rng", cfr = "cfr_rng", dcr = "dcr_rng", tim = "tim_rng",
                 eco = "eco_rng", dwell = "dwell_rng", fluo = "fluo_rng",
                 tid_len = "trk_tid_len", dur = "trk_dur_filter")
  trk_filt    <- reactiveValues()                    # key -> c(lo, hi); NULL = no filter
  .filt_shown <- new.env(parent = emptyenv())        # key -> list(val, bounds) last displayed

  # The value a widget shows: the stored filter clamped into the dataset's bounds,
  # or the full bounds when there is no filter. Remembered, so the echo of our own
  # render is not mistaken for a user edit.
  .filt_value <- function(key, bounds) {
    sp  <- isolate(trk_filt[[key]])
    val <- if (is.null(sp)) bounds
           else c(min(max(sp[1], bounds[1]), bounds[2]), max(min(sp[2], bounds[2]), bounds[1]))
    assign(key, list(val = val, bounds = bounds), envir = .filt_shown)
    val
  }

  # A user edit becomes the stored filter; covering the full range clears it.
  for (.k in names(.FILT_IDS)) local({
    key <- .k; id <- .FILT_IDS[[.k]]
    observeEvent(input[[id]], {
      v <- suppressWarnings(as.numeric(input[[id]]))
      if (length(v) != 2L || any(!is.finite(v))) return(invisible(NULL))
      sh <- if (exists(key, envir = .filt_shown, inherits = FALSE)) get(key, envir = .filt_shown) else NULL
      if (!is.null(sh) && isTRUE(all.equal(v, sh$val))) return(invisible(NULL))   # our own render
      b <- sh$bounds
      trk_filt[[key]] <- if (!is.null(b) && v[1] <= b[1] && v[2] >= b[2]) NULL else sort(v)
      if (!is.null(sh)) assign(key, list(val = v, bounds = b), envir = .filt_shown)
    }, ignoreInit = TRUE)
  })
  rm(.k)

  # The prefactor a in use: the value entered (or set by a calibration profile),
  # else the ideal-estimator value of the TCP pattern, which is uncalibrated.
  .a_in_use <- function(pattern = input$trk_tcp_pattern) {
    v <- suppressWarnings(as.numeric(input$trk_crlb_pref))
    if (length(v) && is.finite(v[1]) && v[1] > 0) v[1]
    else .tcp_geom(pattern %||% TCP_PATTERN_DEFAULT)$a_ideal
  }

  # The per-localization range filters, as filter_ranges() takes them.
  trk_filt_ranges <- reactive(stats::setNames(
    lapply(c("efo", "cfr", "dcr", "tim", "eco", "dwell", "fluo"), function(k) trk_filt[[k]]),
    c("efo", "cfr", "dcr", "tim", "eco", "dwell", "fluo")))

  # Every filter setting, for running the whole pipeline on any dataset
  # (filter_pipeline() in R/est_clean.R): Plot Data's rechunk and Compare Cells.
  # t_c / sw belong to the dataset and are filled in by the caller.
  trk_filter_params <- function() isolate(list(
    split   = list(enable = isTRUE(input$trk_split_gaps_enable), adaptive = isTRUE(input$trk_split_adaptive),
                   gap_cycles = input$trk_split_gap_cycles, min_len = input$trk_split_min_len %||% 50L,
                   maxoff_ms = tryCatch({ sp <- trk_seq_params()
                                          if (!is.null(sp) && isTRUE(sp$ok) && is.finite(sp$maxOffTime_ms))
                                            sp$maxOffTime_ms else 3.0 }, error = function(e) 3.0)),
    trim    = list(enable = isTRUE(input$trk_endtrim_enable), dark = isTRUE(input$trk_endtrim_dark),
                   dark_cycles = input$trk_endtrim_dark_cycles %||% 3L, efo = isTRUE(input$trk_endtrim_efo),
                   efo_frac = suppressWarnings(as.numeric(input$trk_endtrim_efo_frac %||% 0.5))),
    gate    = list(enable = isTRUE(input$trk_quality_enable),
                   min_score = suppressWarnings(as.numeric(input$trk_quality_min_integ %||% 0.6)),
                   half_life = suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3))),
    ranges  = trk_filt_ranges(), tid_len = trk_filt$tid_len, dur = trk_filt$dur))

  # ===========================================================================
  #  TRACKING FILTER PIPELINE  (read top-to-bottom; each stage feeds the next)
  #    raw → split → end-trim → continuity gate → metadata → length → duration
  #  The cycle-time + TSA machinery lives further down, "CYCLE-TIME & TSA SUPPORT".
  # ===========================================================================

  # Cleaning-apply trigger. Split / end-trim / continuity-gate are HEAVY (per-track
  # passes over ALL localizations), so they must not recompute on every slider
  # touch. They read their parameters via isolate() and recompute only when this
  # bumps - on Plot Data, a new upload, or Filter & Rechunk. (trk_raw_df is still a
  # live dependency, so a new manifest also refreshes them.) Metadata / length /
  # duration filters stay LIVE (they read the cached cleaning output).
  clean_trigger <- reactiveVal(0L)
  bump_clean    <- function() clean_trigger(isolate(clean_trigger()) + 1L)

  # Stage 1 (gap splitting) - the FIRST cleaning stage, on the RAW data, so dark-
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
    # The split itself (idempotent on an existing "_s<N>" suffix, radix-sorted,
    # short fragments dropped) is split_gaps() in R/est_clean.R, shared with the
    # per-cell comparison.
    split_gaps(df, gap_thresh, min_len)
  })

  # Stage 2.45: end-trim (optional; pass-through when disabled).
  # Photobleaching degrades the ENDS of a track: fewer photons (lower EFO), worse
  # precision, and blinky dark-gap-separated localizations. These inflate Rmax /
  # MSD tails and produce spurious "exports at track ends". Peel leading & trailing
  # localizations that are reached across a dark gap (≥ the dark-cycle threshold) or
  # are dim (EFO < frac × the track's median), from each end inward, until a solid,
  # well-connected localization is reached. Idempotent - a second pass stops at the
  # same boundary. Runs per segment (after gap-split), so it also cleans mid-track
  # segment ends, not just the bleaching tail.
  trk_end_trimmed <- reactive({
    df <- trk_gap_split()                  # inherits the clean_trigger recompute
    if (!isTRUE(isolate(input$trk_endtrim_enable)) || nrow(df) == 0 || !"tim" %in% names(df))
      return(df)
    # Length floor: gap splitting promises segments of at least `min_len`
    # localizations, so trimming stops there rather than breaking that promise.
    floor_len <- if (isTRUE(isolate(input$trk_split_gaps_enable)))
                   max(1L, as.integer(isolate(input$trk_split_min_len) %||% 1L)) else 1L
    trim_track_ends(df, tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_),
                    use_dark    = isTRUE(isolate(input$trk_endtrim_dark)),
                    dark_cycles = isolate(input$trk_endtrim_dark_cycles %||% 3L),
                    use_efo     = isTRUE(isolate(input$trk_endtrim_efo)),
                    efo_frac    = suppressWarnings(as.numeric(isolate(input$trk_endtrim_efo_frac %||% 0.5))),
                    floor_len   = floor_len)                 # R/est_clean.R
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
    gate_continuity(df, tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_),
                    min_score = suppressWarnings(as.numeric(isolate(input$trk_quality_min_integ) %||% 0.6)),
                    half_life = suppressWarnings(as.numeric(isolate(input$trk_quality_darktol) %||% 3)))
  })

  # Stage 5: trace ID length - localization count per TID, measured AFTER cleaning
  # (split/trim/quality) and the metadata filter.
  trk_len_filtered <- reactive({
    filter_tid_length(trk_meta_filtered(), trk_filt$tid_len)   # R/est_clean.R
  })

  # Stage 4: trace duration - max(tim)-min(tim) per TID measured AFTER stages 1,2 and 3.
  trk_dur_filtered <- reactive({
    filter_duration(trk_len_filtered(), trk_filt$dur)          # R/est_clean.R
  })

  # Stage 5: spatial crop - keeps TIDs where >= pct_thr fraction of their
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

  # Debounced views of the filter chain - UIs and the histogram use these so they
  # don't recompute on every slider tick; they settle with fire_frames (400 ms).
  # The length slider auto-detects from the split output (gap_split_d); the
  # duration slider auto-detects from the length-filtered output (len_filtered_d) -
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
  # half_sample_mode() lives in R/est_clean.R.

  # Typical Δt (ms) of the RAW data - the empirical acquisition interval, bin-free
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
  # dropouts) - exactly where the half-sample mode would mis-report the cycle.
  # Returns list(cycle, sw, n_clusters, centers, modal, frac_one_cycle).
  # Reads trk_raw_df (upstream of all cleaning → no circular dependency; raw Δt is
  # the true acquisition timing).
  trk_cluster_cycle <- reactive({
    df <- tryCatch(trk_raw_df(), error = function(e) NULL)
    if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df))) return(NULL)
    cycle_cluster(df$tid, df$tim)                     # R/est_clean.R
  })

  # Regression-based cycle/SW estimate - Bela's EOD/SW fit, reduced to a line fit.
  # The cluster centres lie at ~ k·cycle + SW for integer cycle counts k, so a
  # least-squares fit of centre ~ k gives slope = cycle time and intercept = SW.
  # This is the same quantity Bela obtains by fitting EOD_OVERHEAD / SW_OVERHEAD to
  # the time signal (EOD is just a reparametrisation of the slope through the
  # hardware vertices). Using ALL centres at once is more precise than the median
  # spacing and is robust to an occasional missing cycle. Returns list(cycle, sw,
  # n, r2, eod) - eod filled only when the hardware vertices/dwell are known.
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

  # Cycle time (ms) from the hardware model - TSA formula derived from Vogler:
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
  # Upper end of the same fit window, also entered in CYCLES and converted to ms.
  # Same "0 means off" convention as msd_min_lag_ms above.
  msd_fit_max_ms <- reactive({
    cyc <- suppressWarnings(as.numeric(input$trk_msd_fit_max_cyc))   # user value, in cycles
    tc  <- trk_cycle_time()                                          # ms per cycle
    # Any missing/nonsensical input (or unknown cycle time) disables the bound.
    if (length(cyc) == 0 || !is.finite(cyc) || cyc <= 0 || !is.finite(tc) || tc <= 0) return(0)
    cyc * tc                                                         # cycles -> milliseconds
  })

  # Status / cross-check line shown in the Acquisition Parameters panel
  output$trk_cycle_time_status <- renderUI({
    req(trk_has_data())
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
      msg <- sprintf("Atomic cycle (model) %.4g ms vs typical Δt %.4g ms (ratio %.2f) - %s",
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
          if (fo >= 0.7)
            sprintf("data quality (raw): GOOD, %.0f%% single-cycle (integrated) steps.", fo * 100)
          else if (fo >= 0.5)
            sprintf("data quality (raw): MODERATE, %.0f%% single-cycle steps; some dark-cycle dropouts.", fo * 100)
          else
            sprintf("data quality (raw): POOR, %.0f%% single-cycle steps; frequent dark cycles.", fo * 100)
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
    # Regression (EOD/SW) cycle - Bela-style line fit, shown for side-by-side
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
          "(length / duration / spatial filters remove the rest - their ranges now expand to cover segments)"))
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
  # recovered from the DATA (the robust cluster method) - the sequence is used to
  # cross-check the scan geometry and, most usefully, to set the gap-split
  # threshold from the hardware's own maxOffTime ("particle lost" timeout).
  .pattern_vertices <- function(p) {
    if (is.null(p) || is.na(p)) return(NA_integer_)
    switch(tolower(p),
           hexagon = 6L, square = 4L, triangle = 3L, octahedron = 6L,
           zline = 2L, zline2 = 2L, line = 2L,
           NA_integer_)
  }

  # .crlb_doughnut() (photon-budget σ, pattern-aware) and .tcp_geom() live in est_cve.R.

  # TCP diameter of one iteration, in nm. Sequence files store the pattern size as
  # the dimensionless `patGeoFactor`; L = patGeoFactor x TCP_L_PER_GEOFACTOR_NM.
  # A 3D iteration may carry a [lateral, axial] pair (e.g. zline gives [0.8, 4.0]);
  # the LATERAL element is the one that governs escape from the pattern, and it is
  # also the L the D_max form is written for, so the first element is used.
  .itr_L_nm <- function(it) {
    pg <- it[["patGeoFactor"]]
    if (is.null(pg)) return(NA_real_)
    pg <- suppressWarnings(as.numeric(unlist(pg)))
    pg <- pg[is.finite(pg) & pg > 0]
    if (!length(pg)) return(NA_real_)
    base <- suppressWarnings(as.numeric(input$trk_geofactor_base %||% TCP_L_PER_GEOFACTOR_NM))
    if (!is.finite(base) || base <= 0) base <- TCP_L_PER_GEOFACTOR_NM
    pg[1] * base
  }

  # ── trk_seq_params() - read the Imspector imaging sequence JSON ─────────────
  # GOAL: pull the ACQUISITION SETTINGS straight from the microscope's sequence
  # file, so the app does not have to infer everything from the data alone.
  # Used as a cross-check on the data-driven cycle time, and to pre-fill the
  # gap-split threshold (from maxOffTime) and the hardware pattern.
  #
  # The file's `Itr` array is one entry per MINFLUX iteration; the last one is
  # the tracking iteration we care about.
  # returns: list(ok = FALSE, msg = ...) on any problem, else the parsed fields.
  #          Callers must always check $ok first.
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
    # maxOffTime may sit on the final iteration or at the top level, and is often
    # the string "unspecified" on stock tracking sequences. Prefer the iteration
    # value, fall back to the top-level one, else leave NA for the 3 ms default.
    moff  <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "maxOffTime")))   # seconds
    if (!is.finite(moff)) moff <- suppressWarnings(as.numeric(getp(seq, "maxOffTime")))
    prep  <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "patRepeat")))
    pht   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "phtLimit")))
    pdw   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "patDwellTime"))) # seconds (timeout)
    # Timing inputs for the blur model (est_blur.R): the CFR check adds a centre dwell
    # of ctrDwellFactor x dwell when ccrLimit > 0, and `damping` > 0 smooths the
    # exported positions across steps, which no blur model describes.
    ccr   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "ccrLimit")))
    cdf   <- suppressWarnings(as.numeric(getp(seq, "ctrDwellFactor")))
    dmp   <- suppressWarnings(as.numeric(getp(seq, "damping")))
    final_pat <- tail(pats[!is.na(pats)], 1L)
    if (length(final_pat) == 0) final_pat <- NA_character_
    dim_final <- suppressWarnings(as.integer(tail(dims[!is.na(dims)], 1L)))
    # TCP diameter per iteration, and the final one - the iteration a tracking
    # sequence reiterates once locked on, which is the L the D_max ceiling needs.
    # The whole ladder is returned too, so a sequence that loops an earlier
    # iteration can be read off and entered by hand.
    L_by_itr <- vapply(itrs, .itr_L_nm, numeric(1))
    L_final  <- if (length(L_by_itr)) L_by_itr[length(L_by_itr)] else NA_real_
    # `headstart` is the loop-back offset: -1 (every tracking sequence seen) means
    # the last iteration is the one repeated. Reported, not acted on.
    hs <- suppressWarnings(as.integer(getp(seq, "headstart")))
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
      patDwellTime_ms = if (is.finite(pdw)) pdw * 1000 else NA_real_,
      cfr_on = is.finite(ccr) && ccr > 0,
      ctrDwellFactor = if (is.finite(cdf)) cdf else NA_real_,
      damping = if (is.finite(dmp)) dmp else NA_real_,
      L_by_itr = L_by_itr,
      L_tcp_nm = L_final,
      headstart = if (is.finite(hs)) hs else NA_integer_
    )
  })

  # ── trk_blur() - ONE MOTION-BLUR COEFFICIENT FOR EVERY TRACKING ESTIMATOR ───
  # GOAL: the blur R of a localization, from the photon-collection timeline of the
  # acquisition (est_blur.R), converted once into the forms the estimators use:
  # R_abs (ms, JDD / D-map / change-point / HMM), R_msd (relative to the cycle, MSD
  # fits and the confinement null) and R_step (the CVE).
  # Inputs, best first: the sequence JSON (pattern, dwell, pattern repeats, CFR
  # centre dwell); otherwise the data (collection window eco/efo, the hardware
  # pattern setting). The step dt is the mean single-cycle interval of the data.
  # eta (passes per localization) = eco / (efo * dwell) when the dwell is known.
  trk_blur <- reactive({
    df  <- trk(); t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    out <- list(R_step = BLUR_R_STEP_FALLBACK, R_abs = NA_real_, R_msd = MSD_R_FALLBACK,
                source = "fallback (standard 2D sequence)", damping = NA_real_, detail = NULL)
    if (is.null(df) || !nrow(df) || !is.finite(t_c) || t_c <= 0) return(out)
    d  <- df[order(df$tid, df$tim), c("tid", "tim", intersect(c("eco", "efo"), names(df)))]
    m  <- nrow(d); same <- c(FALSE, d$tid[-1] == d$tid[-m])
    dt <- c(NA_real_, diff(d$tim)); dt[!same] <- NA_real_
    one <- is.finite(dt) & round(dt / t_c) == 1
    dt_int <- if (any(one)) mean(dt[one]) else t_c
    med <- function(v) { v <- v[is.finite(v) & v > 0]; if (length(v)) stats::median(v) else NA_real_ }
    t_coll <- if (all(c("eco", "efo") %in% names(d))) med(d$eco / d$efo) * 1000 else NA_real_   # ms
    sp <- tryCatch(if (is.null(input$trk_seq_json)) NULL else trk_seq_params(), error = function(e) NULL)
    has_seq <- !is.null(sp) && isTRUE(sp$ok) && is.finite(sp$patDwellTime_ms)
    b <- NULL; note <- NULL
    if (has_seq) {
      eta <- if (is.finite(t_coll)) max(1, round(t_coll / sp$patDwellTime_ms)) else 1
      b <- blur_R_step(dt_int, sp$patDwellTime_ms, pattern = sp$final_pattern,
                       n_vert = sp$vertices, pat_repeat = sp$patRepeat %||% 1, eta = eta,
                       cfr = isTRUE(sp$cfr_on), ctr_factor = sp$ctrDwellFactor %||% 0)
      src <- sprintf("sequence file (%s, %.0f µs dwell%s%s)", sp$final_pattern %||% "?",
                     1000 * sp$patDwellTime_ms, if (isTRUE(sp$cfr_on)) ", CFR" else "",
                     if (eta > 1) sprintf(", %d passes", eta) else "")
      # A sequence whose collection alone outlasts the measured step does not describe
      # this data (e.g. an imaging sequence loaded with tracking data); one whose final
      # pattern is axial has no lateral R. Both fall back to the data, and say so.
      if (is.null(b) || !is.finite(b$R)) {
        note <- sprintf("the sequence's final pattern (%s) gives no lateral blur", sp$final_pattern %||% "?")
        b <- NULL
      } else if (isTRUE(b$overrun)) {
        note <- sprintf("the sequence needs %.0f µs per localization but the data step is %.0f µs, so it does not match this data",
                        1000 * b$t_loc_ms, 1000 * dt_int)
        b <- NULL
      }
    }
    if (is.null(b) && is.finite(t_coll)) {
      # Lateral polygon only: an axial pattern loaded from a sequence can set the
      # vertex control to 2, which says nothing about the lateral collection.
      nv <- suppressWarnings(as.integer(input$trk_hw_vertices %||% 6L))
      if (!is.finite(nv) || nv < 3L) nv <- 6L
      b <- blur_R_step(dt_int, t_coll, n_vert = nv)
      src <- sprintf("data (collection window eco/efo %.0f µs, %d vertices)", 1000 * t_coll, nv)
    }
    if (!is.null(b) && is.finite(b$R)) {
      o <- blur_offsets(b$R, dt_int, t_c)
      out <- c(o, list(source = src, detail = b))
    } else {
      out$R_abs <- blur_offsets(out$R_step, dt_int, t_c)$R_abs
    }
    out$note    <- note
    out$damping <- if (!is.null(sp) && isTRUE(sp$ok)) sp$damping else NA_real_
    out$dt_int  <- dt_int
    out
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
    # D_max ceiling: take L from the sequence rather than a default. The ceiling
    # scales as L^2, so this is the input that most needs to come from the file.
    if (is.finite(sp$L_tcp_nm) && sp$L_tcp_nm > 0)
      updateNumericInput(session, "trk_L_tcp", value = round(sp$L_tcp_nm, 1))
    # TCP geometry for the photon-budget sigma (off-centre coefficient, axes
    # spanned). Only the patterns with a modelled geometry are pre-selected; any
    # other final pattern leaves the current choice alone.
    fp <- tolower(sp$final_pattern %||% "")
    if (fp %in% names(TCP_PATTERN_GEOM))
      updateSelectInput(session, "trk_tcp_pattern", selected = fp)
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
      sprintf("Tracking iteration: %s  ·  dim = %s  ·  L_TCP = %s",
              if (length(sp$final_pattern) && !is.na(sp$final_pattern)) sp$final_pattern else "?",
              if (is.finite(sp$dim)) sprintf("%dD", sp$dim) else "?",
              if (is.finite(sp$L_tcp_nm)) sprintf("%.1f nm", sp$L_tcp_nm) else "?"), tags$br(),
      if (any(is.finite(sp$L_by_itr)))
        tags$span(sprintf("Iteration L: %s nm",
                          paste(ifelse(is.finite(sp$L_by_itr), sprintf("%.0f", sp$L_by_itr), "?"),
                                collapse = " · ")), tags$br()) else NULL,
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
  # any naive cycle estimate - a direct read on whether the data are clean.
  output$trk_tsa_quality <- renderUI({
    req(isTRUE(input$trk_show_tsa), tsa_ran())
    cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    if (is.null(cc) || !is.finite(cc$cycle) || cc$cycle <= 0)
      return(tags$p(style = "color:#888;", "Upload tracking data to run the time-signature analysis."))
    modal <- cc$modal
    ratio <- if (is.finite(modal) && cc$cycle > 0) modal / cc$cycle else NA_real_
    frac1 <- cc$frac_one_cycle                       # RAW (metadata-filtered)
    # Verdict on the POST-FILTERED single-cycle fraction - the SAME data the
    # Integrated/Dark charts below use (dt_vectors → trk(), i.e. after the quality
    # gate / end-trim / length / duration / spatial filters). frac1 (raw) is kept as
    # a diagnostic row so the filter improvement is visible here too. Same metric
    # and thresholds as the Calibration panel, so every verdict agrees.
    ccb       <- tryCatch(cycle_class_bundle(), error = function(e) NULL)
    frac_filt <- if (!is.null(ccb) && is.finite(ccb$frac_integ)) ccb$frac_integ else NA_real_
    frac_use  <- if (is.finite(frac_filt)) frac_filt else frac1
    verdict <- if (is.finite(frac_use)) {
                 if (frac_use >= 0.7) c("GOOD", "#1a8a1a")
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
        if (identical(verdict[1], "GOOD")) "Mostly single-cycle steps."
        else if (identical(verdict[1], "MODERATE"))
          "Many steps span 2+ cycles - try the continuity gate / end-trim (Display tab)."
        else if (identical(verdict[1], "POOR"))
          "Frequent dropouts - filter dark-heavy tracks and ends before trusting D."
        else "Too little structure to judge; falling back on the modal/cluster ratio."))
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
    # Bela computes suggest_split_thresholds AFTER filter_by_length(50), so in
    # adaptive mode only tracks meeting the minimum length count.
    min_len <- if (isTRUE(input$trk_split_adaptive))
                 max(1L, as.integer(input$trk_split_min_len %||% 50L)) else NULL
    total_cycles(df$tid, df$tim, cyc, sw, min_len)    # R/est_clean.R
  })

  # Bela's adaptive split threshold (in cycles): med_th = floor(median + std + 1)
  # of the TOT signal, capped by the hardware "particle lost" ceiling
  # high_th = round(maxOffTime / cycle).  maxOffTime comes from the uploaded
  # sequence JSON if present, else Bela's 3 ms default.
  trk_adaptive_gap_cycles <- reactive({
    moff_ms <- tryCatch({ sp <- trk_seq_params(); if (!is.null(sp) && isTRUE(sp$ok) &&
                            is.finite(sp$maxOffTime_ms)) sp$maxOffTime_ms else 3.0 },
                        error = function(e) 3.0)
    adaptive_gap_cycles(trk_tot_cycles(), trk_cycle_time(), moff_ms)   # R/est_clean.R
  })

  # (The filter pipeline stages 1.5–4 are defined together near Stage 1, under
  # the "TRACKING FILTER PIPELINE" banner. trk_tot_cycles / trk_adaptive_gap_cycles
  # above are the cycle-time support those stages read.)

  output$trk_chunk_stat <- renderUI({
    req(input$trk_chunk_sel)
    if (input$trk_chunk_sel == "__all__") {
      return(tags$div(style = "color:#555;", "All TIDs selected  -  loading full file."))
    }
    df <- read_chunk_file(input$trk_chunk_sel)
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
      tags$b("Trace lifetimes -"),
      sprintf("%d traces \u00B7 mean %.0f ms \u00B7 median %.0f ms", n, avg, med)
    )
  })

  # Stage 3 (duration) auto-detection. BOUNDS reflect the length-filtered population
  # but update at PLOT cadence only - keyed on data_ver() (upload / Plot Data / Reset
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
    val <- .filt_value("dur", bounds)          # stored filter, clamped (full range if none)
    shiny::isolate(reset_dur_filter_flag(FALSE))
    numericRangeInput("trk_dur_filter", "Trace duration filter (ms)",
                      value = val, min = bounds[1], max = bounds[2])
  })

  # ---- Filter + rechunk (driven by Plot Data) ---------------------------------
  # Always reads from original_manifest_dir (frozen at upload) so repeated
  # rechunks are idempotent - Reset Ranges + Rechunk restores the full dataset.
  # Writes filtered result to a new dir that becomes the new manifest_dir.
  observeEvent(rechunk_request(), {
    req(!is.null(rv$original_manifest_dir))
    manf <- file.path(rv$original_manifest_dir, "manifest.csv")
    req(file.exists(manf))
    # Read BEFORE anything is rebuilt: was the user scoped to a chunk, or on all TIDs?
    sel_before <- isolate(input$trk_chunk_sel)

    withProgress(message = "Applying filters: loading chunks...", value = 0, {
      man <- readr::read_csv(manf, show_col_types = FALSE)
      n_chunks <- nrow(man)

      # The SAME pipeline the live view runs (filter_pipeline, R/est_clean.R), on all
      # of the original data at once and with the stored filters, so what Plot Data
      # bakes into the new chunks cannot drift from what the analyses read. (This
      # used to be a second, hand-written copy of every stage, which had drifted.)
      incProgress(0.1, detail = "Reading chunks")
      df_all <- dplyr::bind_rows(Filter(Negate(is.null), lapply(man$file, read_chunk_file)))
      if (!nrow(df_all)) {
        showNotification("No data to filter  -  rechunk aborted.", type = "error", duration = 8)
        return(invisible(NULL))
      }
      # Same post-processing as trk_raw_df
      df_all$tid <- as.character(df_all$tid)
      df_all$x   <- as.numeric(df_all$x)
      df_all$y   <- as.numeric(df_all$y)
      df_all$tim <- as.numeric(df_all$tim) * 1000           # seconds -> ms
      fp <- trk_filter_params()
      fp$t_c <- isolate(trk_cycle_time())
      cc_r   <- tryCatch(isolate(trk_cluster_cycle()), error = function(e) NULL)
      fp$sw  <- if (!is.null(cc_r) && is.finite(cc_r$sw)) cc_r$sw else 0
      incProgress(0.3, detail = "Applying filters")
      df_all <- filter_pipeline(df_all, fp)
      if (nrow(df_all) == 0) {
        showNotification("No data survived the filters  -  rechunk aborted.", type = "error", duration = 8)
        return(invisible(NULL))
      }

      # Spatial filter - apply if active, using the same pct_thr threshold
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
      user_tpc <- .tids_per_chunk()

      manf_new <- split_by_tid_to_dir(df_all, out_dir, tids_per_chunk = user_tpc, fmt = CHUNK_FORMAT)
      if (!file.exists(manf_new)) {
        showNotification("Rechunking failed (split_by_tid_to_dir returned no manifest).", type = "error", duration = 8)
        return(invisible(NULL))
      }

      readr::write_csv(tid_sum_new, file.path(out_dir, "tid_summary.csv"))
      # Carry the scoping, not the index: on a chunk before -> chunk 1 of the new
      # chunking; on "All TIDs" before -> still all. Set before the manifest swap so
      # the selector rebuild below already sees it.
      trk_chunk_pref(if (isTruthy(sel_before) && !identical(sel_before, "__all__"))
                       "first" else "__all__")
      old_rechunk         <- rv$last_rechunk_dir   # superseded by out_dir; removed below
      rv$last_rechunk_dir <- out_dir
      rv$manifest_dir     <- out_dir
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

      # Auto-fit the total playback window to the longest remaining trace so the whole
      # longest track plays without being cut off. tid_sum_new is sorted longest-first,
      # so [1] is the max duration (ms). The user can still lower it afterwards.
      max_dur <- suppressWarnings(tid_sum_new$duration_ms[1])
      if (length(max_dur) == 1 && is.finite(max_dur) && max_dur > 0)
        updateNumericInput(session, "trk_play_total", value = signif(max_dur, 4))

      # The split (if any) is now BAKED into the saved segment TIDs. Live splitting
      # can stay ON: trk_gap_split strips the "_s<N>" suffix before splitting, so the
      # reloaded segments are recombined and re-split to the SAME segments rather than
      # split a second time. (That double split is what used to shrink the duration
      # upper bound on every rechunk.) Pressing Filter & Rechunk again is idempotent.
      # Flags stay TRUE so the length/duration sliders snap to the new full range
      # on their next (data-driven) re-render - no forced re-render here, which
      # would fire during the async checkbox round-trip while the split is still on.
      trk_len_reset_flag(TRUE)
      reset_dur_filter_flag(TRUE)

      n_tids <- length(unique(df_all$tid))
      n_new  <- nrow(readr::read_csv(manf_new, show_col_types = FALSE))
      showNotification(
        sprintf("%s TIDs kept, in %d chunk(s) @ %d TIDs/chunk. Showing %s.",
                format(n_tids, big.mark = ","), n_new, user_tpc,
                if (identical(isolate(trk_chunk_pref()), "first") && n_new >= 1)
                  "chunk 1" else "all TIDs"),
        type = "message", duration = 8
      )

      # Each press writes a fresh chunk set to tempdir(). Now that this runs on
      # every Plot Data, the superseded one is removed rather than left to
      # accumulate copies of the dataset. Done last, so nothing is deleted unless
      # the new set was written and swapped in successfully. The frozen upload set
      # is never a candidate.
      if (!is.null(old_rechunk) && !identical(old_rechunk, out_dir) &&
          !identical(old_rechunk, rv$original_manifest_dir) && dir.exists(old_rechunk))
        unlink(old_rechunk, recursive = TRUE)
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

    # Blue bars: fully-filtered data - debounced so the histogram only recomputes
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

  # ── trk_image_overlay() - place a confocal image behind the tracks ──────────
  # GOAL: turn an uploaded confocal/widefield image into a plotly background
  # layer positioned in DATA COORDINATES (nm), so tracks overlay the cell image.
  #
  # The image is in pixels; the tracks are in nm. Converting between them needs
  # the pixel size and an origin offset, both user-supplied:
  #   nm = pixel_index * trk_px_size + (trk_img_x0 / trk_img_y0)
  # returns NULL (draw nothing) if there is no image or the pixel size is unset -
  # every caller treats NULL as "no overlay".
  trk_image_overlay <- reactive({
    meta <- trk_confocal_meta()                                  # image + its dimensions
    if (is.null(meta)) return(NULL)                              # no image uploaded
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
      base <- trk_view_base()
      paste0(base, "_filter_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content = function(file) {
      sf <- trk_spatial_filter()
      s <- list(
        # ── Metadata (per-localization) ranges ──────────────────────────────
        # The stored filters (null = no filter), not the widgets' clamped display.
        efo_rng   = trk_filt$efo %||% NA,   cfr_rng   = trk_filt$cfr %||% NA,
        dcr_rng   = trk_filt$dcr %||% NA,   tim_rng   = trk_filt$tim %||% NA,
        eco_rng   = trk_filt$eco %||% NA,   dwell_rng = trk_filt$dwell %||% NA,
        fluo_rng  = trk_filt$fluo %||% NA,
        trk_tid_len    = trk_filt$tid_len %||% NA,
        trk_dur_filter = trk_filt$dur %||% NA,

        # ── Track CLEANING - these change the tracks themselves ─────────────
        # Cleaning runs before every metadata filter (split → trim → gate), so
        # omitting it meant a reloaded settings file reproduced the ranges but
        # silently analysed uncleaned tracks: different track IDs, lengths and
        # durations, hence different MSD/D. These matter more than the ranges.
        trk_split_gaps_enable   = isTRUE(input$trk_split_gaps_enable),
        trk_split_adaptive      = isTRUE(input$trk_split_adaptive),
        trk_split_gap_cycles    = input$trk_split_gap_cycles,
        trk_split_min_len       = input$trk_split_min_len,

        trk_endtrim_enable      = isTRUE(input$trk_endtrim_enable),
        trk_endtrim_dark        = isTRUE(input$trk_endtrim_dark),
        trk_endtrim_dark_cycles = input$trk_endtrim_dark_cycles,
        trk_endtrim_efo         = isTRUE(input$trk_endtrim_efo),
        trk_endtrim_efo_frac    = input$trk_endtrim_efo_frac,

        trk_quality_enable      = isTRUE(input$trk_quality_enable),
        trk_quality_darktol     = input$trk_quality_darktol,
        trk_quality_min_integ   = input$trk_quality_min_integ,

        # ── Spatial crop (also a data filter) ───────────────────────────────
        trk_crop_pct   = input$trk_crop_pct,
        trk_spatial_x  = if (!is.null(sf)) unname(sf$x) else NULL,
        trk_spatial_y  = if (!is.null(sf)) unname(sf$y) else NULL,
        trk_spatial_pct_thr = if (!is.null(sf)) sf$pct_thr else NULL,

        # ── Histogram panel display state ───────────────────────────────────
        trk_show_metadata_hist = isTRUE(input$trk_show_metadata_hist),
        trk_hist_param         = input$trk_hist_param,
        trk_hist_bins          = input$trk_hist_bins,
        trk_hist_show_raw      = isTRUE(input$trk_hist_show_raw)
      )
      writeLines(jsonlite::toJSON(s, pretty = TRUE, auto_unbox = TRUE, na = "null"), file)
    }
  )
  observeEvent(input$trk_load_filter_settings, {
    req(input$trk_load_filter_settings$datapath)
    s <- jsonlite::fromJSON(input$trk_load_filter_settings$datapath)

    # Every restore is guarded on the key being PRESENT, so an older settings file
    # (saved before cleaning was persisted) leaves the current cleaning setup alone
    # instead of silently switching split/trim/gate off.
    num <- function(v) suppressWarnings(as.numeric(unlist(v)))
    upd_rng <- function(id) {
      v <- s[[id]]
      if (!is.null(v) && length(unlist(v)) == 2)
        shinyWidgets::updateNumericRangeInput(session, id, value = num(v))
    }
    upd_num <- function(id) {
      v <- s[[id]]
      if (!is.null(v)) { n <- num(v)[1]; if (is.finite(n)) updateNumericInput(session, id, value = n) }
    }
    upd_sld <- function(id) {
      v <- s[[id]]
      if (!is.null(v)) { n <- num(v)[1]; if (is.finite(n)) updateSliderInput(session, id, value = n) }
    }
    upd_chk <- function(id) if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))

    # Range filters go straight into the store (a missing key leaves a filter as it
    # is, null clears it); the widgets then re-render from it.
    for (k in names(.FILT_IDS)) {
      id <- .FILT_IDS[[k]]
      if (!id %in% names(s)) next
      v <- num(s[[id]])
      trk_filt[[k]] <- if (length(v) == 2L && all(is.finite(v))) sort(v) else NULL
    }
    reset_ranges_trigger(isolate(reset_ranges_trigger()) + 1L)

    # cleaning: enables first, then their parameters
    for (id in c("trk_split_gaps_enable","trk_split_adaptive",
                 "trk_endtrim_enable","trk_endtrim_dark","trk_endtrim_efo",
                 "trk_quality_enable",
                 "trk_show_metadata_hist","trk_hist_show_raw"))
      upd_chk(id)
    for (id in c("trk_split_gap_cycles","trk_split_min_len",
                 "trk_endtrim_dark_cycles","trk_endtrim_efo_frac","trk_crop_pct"))
      upd_num(id)
    for (id in c("trk_quality_darktol","trk_quality_min_integ","trk_hist_bins"))
      upd_sld(id)
    if (!is.null(s$trk_hist_param))
      updateSelectInput(session, "trk_hist_param", selected = as.character(s$trk_hist_param))

    # Spatial crop: restore the stored box directly (it is a reactiveVal, not an
    # input) and bump the data version so the filter is actually re-applied.
    if (!is.null(s$trk_spatial_x) && !is.null(s$trk_spatial_y) &&
        length(unlist(s$trk_spatial_x)) == 2 && length(unlist(s$trk_spatial_y)) == 2) {
      pct <- num(s$trk_spatial_pct_thr)[1]
      if (!is.finite(pct) || pct <= 0 || pct > 1) pct <- 0.5
      trk_spatial_filter(list(x = num(s$trk_spatial_x), y = num(s$trk_spatial_y), pct_thr = pct))
      data_ver(isolate(data_ver()) + 1L)
      fire_frames()
    }

    cleaning_on <- any(vapply(c("trk_split_gaps_enable","trk_endtrim_enable","trk_quality_enable"),
                              function(k) isTRUE(s[[k]]), logical(1)))
    showNotification(
      paste0("Filter settings loaded.",
             if (cleaning_on) " Cleaning settings restored - press Plot Data to re-apply them." else "",
             if (!is.null(s$trk_spatial_x)) " Spatial crop restored." else ""),
      type = "message", duration = 6)
  })

  # ---- Playback settings save/load ----
  output$trk_download_playback_settings <- downloadHandler(
    filename = function() {
      base <- trk_view_base()
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
        trk_per_tid_cap    = input$trk_per_tid_cap,
        trk_trail_nmax     = input$trk_trail_nmax,
        trk_max_pts_tid    = input$trk_max_pts_tid,
        # Panel-visibility and trail toggles. The advanced numerics above were
        # already saved but their gating checkboxes were not, so a reloaded file
        # restored values the user could not see and trails stayed off.
        trk_show_playback_params = isTRUE(input$trk_show_playback_params),
        trk_playback_ctrl        = isTRUE(input$trk_playback_ctrl),
        trk_show_moving_trails   = isTRUE(input$trk_show_moving_trails),
        trk_show_advanced        = isTRUE(input$trk_show_advanced),
        # View mode lives just above the playback panel and decides whether any of
        # the settings below are even used (playback vs static full tracks), so it
        # belongs with them rather than in the filter file.
        trk_view_mode_main       = input$trk_view_mode_main,
        trk_full_connect         = isTRUE(input$trk_full_connect)
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
                 "trk_trail_len",
                 "trk_per_tid_cap","trk_trail_nmax","trk_max_pts_tid"))
      upd_num(id)
    # Gating checkboxes - restore only when present, so an older playback file
    # does not collapse the panels the user currently has open.
    for (id in c("trk_show_playback_params","trk_playback_ctrl",
                 "trk_show_moving_trails","trk_show_advanced","trk_full_connect"))
      if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))
    if (!is.null(s$trk_view_mode_main))
      updateRadioButtons(session, "trk_view_mode_main", selected = as.character(s$trk_view_mode_main))
    showNotification("Playback settings loaded.", type = "message", duration = 3)
  })

  # Live guidance for the drift-speed window W. W must sit well BELOW a typical track
  # duration: once W >= a track's duration every point in it shares the same window
  # (lo=1, hi=n), so the whole track renders as ONE flat colour and the within-track
  # gradient the mode exists for disappears. Suggest ~1/10 of the median duration, which
  # leaves ~10 independent window positions along a typical track.
  output$trk_speed_win_note <- renderUI({
    if (!identical(input$trk_color_mode, "speed")) return(NULL)
    rs <- trk_rechunk_stats_rv(); ro <- trk_raw_stats_rv()
    dur <- if (!is.null(rs) && length(rs$tid_duration)) rs$tid_duration
           else if (!is.null(ro) && length(ro$tid_duration)) ro$tid_duration else NULL
    dur <- dur[is.finite(dur) & dur > 0]
    if (!length(dur)) return(NULL)
    med  <- stats::median(dur)
    W    <- suppressWarnings(as.numeric(input$trk_color_speed_win))
    sugg <- signif(med / 10, 2)
    flat <- if (is.finite(W)) mean(dur <= W) else 0
    warn <- is.finite(W) && W > med / 3
    div(style = sprintf("color:%s; font-size:0.8em; line-height:1.35; margin-top:-6px;",
                        if (warn) "#b8860b" else "#666"),
        if (warn) "⚠ " else "",
        sprintf("Median track duration %s ms. ", format(signif(med, 3), trim = TRUE)),
        if (flat > 0) sprintf(
          "At W = %s ms, %.0f%% of tracks are shorter than the window - those render as one flat colour. ",
          format(W, trim = TRUE), 100 * flat) else "",
        sprintf("Suggested W ≈ %s ms (~1/10 of a typical track).", format(sugg, trim = TRUE))
    )
  })

  # Warn when the trail is shorter than one bin: each frame advances by one bin width,
  # so a trail shorter than a bin can't bridge consecutive frames - it skips ahead and
  # looks jumpy. Smooth, continuous trails need trail length >= bin width.
  output$trk_trail_warn <- renderUI({
    if (!isTRUE(input$trk_show_moving_trails)) return(NULL)
    tl <- suppressWarnings(as.numeric(input$trk_trail_len))
    bw <- suppressWarnings(as.numeric(input$trk_bin_dt))
    if (!is.finite(tl) || !is.finite(bw) || bw <= 0 || tl >= bw) return(NULL)
    div(style = "color:#b8860b; font-size:12px; line-height:1.3; margin:-4px 0 6px 0;",
        sprintf(paste0("⚠ Trail length (%g ms) is below the bin width (%g ms). The trail ",
                       "advances one bin per frame, so a trail shorter than a bin skips ahead ",
                       "each frame and looks jumpy. Set trail length ≥ bin width for smooth trails."),
                tl, bw))
  })

  output$trk_download_overlay_settings <- downloadHandler(
    filename = function() {
      base <- trk_view_base()
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
      base <- trk_view_base()
      paste0(base, "_filtered_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      df <- trk_spatial_filtered()
      shiny::validate(shiny::need(nrow(df) > 0, "No filtered data to export."))
      drop_cols <- intersect(c("col", "t"), names(df))
      if (length(drop_cols)) df <- df[, !names(df) %in% drop_cols, drop = FALSE]
      # tim is stored internally in ms; convert back to seconds for export
      if ("tim" %in% names(df)) df$tim <- df$tim / 1000
      # Axial provenance: z below is at this factor. Re-importing (here or in the
      # SMLM tab) divides it out instead of falling back to the 0.7 guess in
      # detect_csv_z_calibration(), so a round-trip is exact.
      df$z_scale_applied <- trk_z_applied()
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
    showNotification("Boundary alignment applied - nudge to fine-tune.",
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
  # suppressWarnings can't catch). Stays FALSE until a plot with data renders -
  # layout_trigger alone isn't enough (it can be set while frames_store is empty).
  trk_plot_ready    <- reactiveVal(FALSE)
  # Always NULL since the Rg/Rmax panel was retired (2026-08-03) - the observer that
  # filled it moved to retired_tracking_server.R. Kept so the playback colour override
  # below needs no change, and so reinstating the retired panel just works.
  conf_color_map_stored    <- reactiveVal(NULL)  # pre-computed tid->hex when conf coloring is on
  cluster_color_map_stored <- reactiveVal(NULL)  # pre-computed tid->hex when MSD cluster coloring is on

  # TRUE when the loaded data has a usable z column (enables 3D toggle). Requires real
  # axial spread, so a placeholder all-zero z column does not offer a 3D mode that
  # would silently pool a zero-variance axis (see .has_real_z).
  has_z <- reactive({
    .has_real_z(trk())
  })

  # TRUE when the user has selected 3D mode and data supports it
  use_3d <- reactive({
    isTRUE(has_z()) && identical(input$trk_view_mode %||% "2d", "3d")
  })
  trk_axis_x_span   <- reactiveVal(NULL)
  trk_locked_ranges <- reactiveVal(NULL)   # list(x = c(x0,x1), y = c(y0,y1)) - held fixed during playback
  # Stable key that only changes on new upload/rechunk - used as uirevision so
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

  # ── Plot Data - the single apply action ─────────────────────────────────────
  # Applies everything the user has changed: cleaning (split / end-trim / continuity
  # gate), the metadata / length / duration / spatial filters, re-chunks the
  # survivors at the current TIDs-per-chunk setting, and re-renders.
  #
  # There is no separate rechunk button. The rechunk path already ended in
  # trk_data_key(), which fires the layout, so the two actions were doing the same
  # work with the rechunk merely persisting the result - splitting them into two
  # presses only raised the question of which order to use them in.
  #
  # NO fire_layout() on the rechunk branch: trk_data_key() fires it, and calling
  # both would render the chart twice. The fallback branch covers a dataset with no
  # manifest (chunking is skipped for some file types), where there is nothing to
  # rechunk and a plain rebuild is the whole job.
  observeEvent(input$trk_plot_data, {
    bump_clean()      # (re)apply split / end-trim / continuity gate with current settings
    manf <- if (!is.null(rv$original_manifest_dir))
              file.path(rv$original_manifest_dir, "manifest.csv") else NA_character_
    if (!is.na(manf) && file.exists(manf)) rechunk_request(isolate(rechunk_request()) + 1L)
    else fire_layout("plot_data_btn")
  })
  observeEvent(list(
    input$trk_play_mode,
    input$trk_view_mode_main, input$trk_full_connect, input$trk_full_raw,
    input$trk_show_full_tracks_faint, input$trk_show_moving_trails, input$trk_color_mode,
    input$trk_dark_show, input$trk_color_speed_win, input$trk_chunk_scope_full
    # trk_show_heatmap deliberately NOT here - the dwell-time heatmap is heavy, so it
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
  # was removed). The metadata / length / duration filters - like the split / trim /
  # quality-gate cleaning stages - apply to the plotted data ONLY when the user presses
  # Plot Data (fire_layout → bumps data_ver → reruns trk()) or Filter & Rechunk. Live,
  # filter changes still drive the metadata HISTOGRAM (reads the debounced filtered
  # stages directly) and the length/duration range read-outs (re-keyed on data_ver, so
  # they describe the plotted data) - but the WebGL track render stays put until an
  # explicit Plot Data / Rechunk. See observeEvent(input$trk_plot_data) above.


  # ---- Prepare data once per filter-data change ----
  # trk() is keyed on data_ver, NOT plot_trigger.  This means it skips
  # prep_tracking + tid_color_map for rebuild_tick-only (playback param) triggers
  # where the filtered dataframe is identical to the previous run.
  # frames_store() still runs on every plot_trigger; it reads isolate(trk()) and
  # gets the cached result when data_ver did not change.
  # ── Dark-cycle colouring ─────────────────────────────────────────────────────
  # A step spanning k = round(dt / cycle_time) acquisition cycles has k-1 DARK cycles
  # (k=1 = integrated, i.e. a localization every cycle). Each localization inherits the
  # dark-cycle count of the step LEADING INTO it. Colour scheme: integrated = blue,
  # 1+ dark cycles = light→deep red gradient (capped at DARK_NMAX; higher = deeper red).
  DARK_NMAX <- 6L
  .dark_hex <- function(nmax = DARK_NMAX) {
    reds <- grDevices::colorRampPalette(c("#fca9a4", "#e34a33", "#7f000d"))(nmax)
    c("#2c7fb8", reds)                       # length nmax+1; index = n_dark + 1
  }
  dark_cycle_cols <- function(n_dark, nmax = DARK_NMAX) {
    pal <- .dark_hex(nmax)
    i <- pmin(pmax(as.integer(n_dark), 0L), nmax) + 1L
    i[!is.finite(i)] <- 1L
    pal[i]
  }
  dark_cycle_colorscale <- function(nmax = DARK_NMAX) {
    pal <- .dark_hex(nmax)                   # colours for values 0..nmax
    lapply(0:nmax, function(k) list(k / nmax, pal[k + 1]))
  }

  # ── Speed colouring (drift velocity) ─────────────────────────────────────────
  # Per-step speed |Δr|/Δt is NOT usable here: Brownian motion has no instantaneous
  # velocity, so r/Δt diverges as Δt→0 and mostly measures the sampling rate (see the
  # header note, and why the Speed panel's path speed can't separate diffusion from
  # transport). We therefore colour by the same COHERENT-transport measure the Speed &
  # Transport panel uses: v_drift = |net displacement| / span over a CENTRED window W.
  # Diffusion's net displacement grows only as sqrt(t) so its apparent v ~ sqrt(4D/W)
  # stays low, while sustained directed motion reads high - which is what makes a
  # confined stretch read blue and a directed excursion read red.
  # Vectorised per track via findInterval on the window edges → O(N) overall.
  trk_point_drift_speed <- function(df, win_ms = 200, min_pts = 3L) {
    n_all <- nrow(df)
    v <- rep(NA_real_, n_all)
    if (!is.finite(win_ms) || win_ms <= 0 || n_all < 2L) return(v)
    has_z <- "z" %in% names(df) && any(is.finite(df$z))
    hw <- win_ms / 2
    for (ii in split(seq_len(n_all), df$tid)) {
      n <- length(ii); if (n < 2L) next
      t <- df$tim[ii]; xx <- df$x[ii]; yy <- df$y[ii]
      zz <- if (has_z) df$z[ii] else NULL
      lo <- pmin(pmax(findInterval(t - hw, t) + 1L, 1L), n)   # first pt >= t-W/2
      hi <- pmin(pmax(findInterval(t + hw, t),      1L), n)   # last  pt <= t+W/2
      span <- t[hi] - t[lo]
      dr <- sqrt((xx[hi] - xx[lo])^2 + (yy[hi] - yy[lo])^2 +
                 (if (has_z) (zz[hi] - zz[lo])^2 else 0))
      okw <- is.finite(span) & span > 0 & is.finite(dr) & (hi - lo + 1L) >= min_pts
      vv <- rep(NA_real_, n); vv[okw] <- dr[okw] / span[okw]   # nm/ms = µm/s
      v[ii] <- vv
    }
    v
  }
  # Blue (slow / confined) → red (fast / directed). NA (window can't be formed, e.g.
  # at track ends) stays GREY rather than being painted a misleading "slow" blue.
  .speed_hex <- function(n = 64L)
    grDevices::colorRampPalette(c("#2c7fb8", "#41b6c4", "#fecc5c", "#fd8d3c", "#e31a1c"))(n)
  speed_cols <- function(v, vmax) {
    pal <- .speed_hex()
    out <- rep("#9e9e9e", length(v))
    ok  <- is.finite(v) & isTRUE(is.finite(vmax) && vmax > 0)
    if (any(ok)) {
      i <- 1L + as.integer(round(pmin(pmax(v[ok] / vmax, 0), 1) * (length(pal) - 1L)))
      out[ok] <- pal[i]
    }
    out
  }
  speed_colorscale <- function() {
    pal <- .speed_hex(9L)
    lapply(seq_along(pal), function(k) list((k - 1) / (length(pal) - 1), pal[k]))
  }

  # Colorbar-only trace. Per-point colours are baked as hex strings (needed by the
  # restyle-driven playback path), and hex arrays carry no colourbar - so we append one
  # invisible 2-point trace whose sole job is to render the scale. Added LAST so it can
  # never shift trail_start / markers_idx in trace_map.
  # ax/ay/az anchor the trace on REAL coordinates: plotly (scattergl especially) drops a
  # trace whose points are all NA, and a dropped trace draws no colourbar. The markers
  # are made invisible with opacity 0 + size 0 instead.
  add_colorbar_trace <- function(fig, mode, fg_col, vmax = NA_real_, use3d = FALSE,
                                 ax = 0, ay = 0, az = 0) {
    if (!mode %in% c("dark", "speed")) return(fig)
    if (!is.finite(ax) || !is.finite(ay)) { ax <- 0; ay <- 0 }
    if (!is.finite(az)) az <- 0
    if (identical(mode, "dark")) {
      cs <- dark_cycle_colorscale(); cmin <- 0; cmax <- DARK_NMAX
      cb <- list(title = list(text = "dark<br>cycles", font = list(color = fg_col)),
                 tickmode = "array", tickvals = 0:DARK_NMAX,
                 ticktext = c(as.character(0:(DARK_NMAX - 1L)), paste0(DARK_NMAX, "+")),
                 tickfont = list(color = fg_col), len = 0.5, thickness = 14)
    } else {
      if (!is.finite(vmax) || vmax <= 0) return(fig)
      cs <- speed_colorscale(); cmin <- 0; cmax <- vmax
      cb <- list(title = list(text = "drift speed<br>(µm/s)", font = list(color = fg_col)),
                 tickfont = list(color = fg_col), len = 0.5, thickness = 14)
    }
    args <- list(x = c(ax, ax), y = c(ay, ay),
                 type = if (use3d) "scatter3d" else "scattergl", mode = "markers",
                 marker = list(size = 0, opacity = 0, color = c(cmin, cmax),
                               colorscale = cs, cmin = cmin, cmax = cmax,
                               showscale = TRUE, colorbar = cb),
                 hoverinfo = "skip", showlegend = FALSE, name = "scale")
    if (use3d) args$z <- c(az, az)
    do.call(add_trace, c(list(fig), args))
  }

  # ── trk() - THE analysis-ready dataset ──────────────────────────────────────
  # GOAL: the single hand-off point between the filter pipeline and everything
  # that analyses or draws the data. Almost every bundle below starts with
  # `df <- trk()`. If a filter change is not showing up in an analysis, this is
  # the first place to look.
  #
  # It does three things:
  #   1. takes the fully-filtered data (last pipeline stage, trk_spatial_filtered)
  #   2. canonicalises it via prep_tracking() - sorted by (tid, time), non-finite
  #      rows dropped, per-track `tim0` clock added
  #   3. assigns each track a stable colour, with split segments of the same
  #      parent track sharing a hue family
  #
  # WHY eventReactive(data_ver()): recomputing on every upstream twitch would be
  # far too expensive. `data_ver()` is a version counter bumped only when the
  # data genuinely changes, so trk() recomputes once per real change. Everything
  # inside is read with isolate() so only data_ver() triggers it.
  trk <- eventReactive(data_ver(), {
    df <- isolate(trk_spatial_filtered()); req(nrow(df) > 0)   # last pipeline stage
    df$tid <- as.character(df$tid)                             # TIDs are character everywhere
    # normalize_per_tid: in simultaneous playback every track must start at t=0
    # so they animate together; in sequential mode real timestamps are kept.
    df <- prep_tracking(df, normalize_per_tid = isolate(isTRUE(is_simul())))
    cmap <- tid_color_map(df$tid)                              # stable colour per TID

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

    # Per-localization dark-cycle count + colour (df is sorted by tid,tim here).
    # dt uses RAW tim (ms); normalization only offsets per-TID so diffs are real.
    t_c <- tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_)
    n2   <- nrow(df)
    same <- if (n2 > 1) c(FALSE, df$tid[-1] == df$tid[-n2]) else FALSE
    dt   <- if (n2 > 1) c(NA_real_, diff(df$tim)) else NA_real_
    dt[!same] <- NA_real_
    if (is.finite(t_c) && t_c > 0) {
      ndk <- pmax(0L, as.integer(round(dt / t_c)) - 1L)   # k-1 dark cycles; k=round(dt/t_c)
    } else ndk <- rep(0L, n2)
    ndk[is.na(ndk)] <- 0L                                 # first loc of each TID → integrated
    df$n_dark   <- ndk
    df$col_dark <- dark_cycle_cols(ndk)
    # (Z is already corrected at the pipeline root in trk_raw_df() by the physical
    # Z factor, so it flows in here through the filter chain - not re-applied.)
    df
  }, ignoreNULL = TRUE)

  # Sequential data (always order by raw 'tim'; decimate for speed)
  # ── Chunk-scoped PLAYBACK ───────────────────────────────────────────────────
  # Chunking exists only to keep the animation smooth: the filter pipeline,
  # analysis, and the full-track static overlay all use ALL TIDs, but the moving
  # markers/trails are limited to the selected chunk's TIDs. Returns the TIDs exactly
  # as the chunk file stores them - original ("42") for upload chunks, segment
  # ("42_s1") once a rechunk has baked the split in; scope_playback() below accepts
  # either. NULL = animate everything ("All TIDs" selected).
  trk_playback_chunk_tids <- reactive({
    sel <- input$trk_chunk_sel
    if (!isTruthy(sel) || identical(sel, "__all__") || !file.exists(sel)) return(NULL)
    ch <- read_chunk_file(sel, columns = "tid")
    if (is.null(ch) || !"tid" %in% names(ch) || !nrow(ch)) return(NULL)
    unique(as.character(ch$tid))
  })

  # Subset a display data.frame to the chunk's TIDs.
  #
  # A chunk file holds one of TWO TID conventions, and the match has to accept both:
  #   * UPLOAD chunks hold the original TIDs ("42"). Gap-splitting happens live
  #     afterwards, so the display frame carries segments ("42_s1") that must still
  #     be matched to their parent chunk - hence the base-TID comparison.
  #   * RE-CHUNKED sets have the split BAKED IN (the rechunk writes "42_s1" as the
  #     TID), so the display frame and the chunk file agree exactly and the base
  #     comparison alone never matches - "42" is not in {"42_s1", ...}. Every row
  #     then failed, the guard below returned the whole frame, and selecting a
  #     chunk appeared to do nothing.
  # Testing both forms covers each case without over-matching: a chunk holding
  # "42_s1" still does not pull in "42_s2", which may belong to a different chunk.
  #
  # Falls back to the full frame if the scope would leave nothing to animate.
  scope_playback <- function(d, chunk_tids) {
    if (is.null(chunk_tids) || !nrow(d) || !"tid" %in% names(d)) return(d)
    tid_d <- as.character(d$tid)
    base  <- sub("_s\\d+$", "", tid_d, perl = TRUE)
    keep  <- tid_d %in% chunk_tids | base %in% chunk_tids
    if (!any(keep)) return(d)            # chunk fully filtered out - don't blank playback
    d[keep, , drop = FALSE]
  }

  # ── seq_data() - prepare tracks for PLAYBACK ────────────────────────────────
  # GOAL: turn the analysis dataset into a per-track structure the animation can
  # step through cheaply. Playback must stay smooth, so this does the expensive
  # work once (scope to chunk, thin, split by track) instead of per frame.
  #
  # Two deliberate cost controls:
  #   * only the SELECTED playback chunk is animated (analyses still use all TIDs)
  #   * each track is decimated to at most `trk_max_pts_tid` points
  # Both are DISPLAY-only; no estimator ever sees this thinned data.
  #
  # Recomputes only when plot_trigger() fires, not on every input change.
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

    # Drift-speed colouring is computed HERE, not in trk(): trk() is keyed on data_ver
    # only, so it would not refresh when the colour mode or window changes. frames_store
    # re-runs on plot_trigger, which fire_layout("display_opts") bumps. Only computed
    # when the speed mode is active, since it is an O(N) pass.
    speed_max <- NA_real_
    if (identical(isolate(input$trk_color_mode %||% "tid"), "speed")) {
      w_ms <- suppressWarnings(as.numeric(isolate(input$trk_color_speed_win %||% 200)))
      if (!is.finite(w_ms) || w_ms <= 0) w_ms <- 200
      df$v_drift <- trk_point_drift_speed(df, win_ms = w_ms)
      # Robust upper limit so a few outlier windows don't flatten the whole scale.
      speed_max <- suppressWarnings(as.numeric(stats::quantile(df$v_drift, 0.95, na.rm = TRUE)))
      if (!is.finite(speed_max) || speed_max <= 0)
        speed_max <- suppressWarnings(max(df$v_drift, na.rm = TRUE))
      if (!is.finite(speed_max) || speed_max <= 0) speed_max <- 1
      df$col_speed <- speed_cols(df$v_drift, speed_max)
    }

    # ── 2D viewing plane (XY / XZ / YZ) ────────────────────────────────────────
    # Bake the plane choice in as a COORDINATE REMAP, so the whole 2D pipeline
    # (full-tracks + playback frames + axes) plots the chosen plane without any
    # per-branch z handling - the same trick the SMLM Display tab uses. Only in 2D
    # mode and only when z is real; 3D mode always keeps native XY(Z). Placed AFTER
    # the drift-speed pass so v_drift stays a true in-plane physical speed.
    plane <- if (!isTRUE(isolate(use_3d()))) (isolate(input$trk_view_plane) %||% "XY") else "XY"
    # Fall back to XY unless z is present AND actually varies (defensive: mirrors the
    # plane-selector gate, so a stale XZ/YZ selection on planar data can't flat-line it).
    z_fin <- if ("z" %in% names(df)) df$z[is.finite(df$z)] else numeric(0)
    if (plane %in% c("XZ", "YZ") && !(length(z_fin) > 0 && diff(range(z_fin)) > 1e-6)) plane <- "XY"
    remap_plane <- function(d) {
      if (is.null(d) || !("z" %in% names(d))) return(d)
      if (plane == "XZ")      d$y <- d$z              # X horizontal, Z vertical
      else if (plane == "YZ") { d$x <- d$y; d$y <- d$z }  # Y horizontal, Z vertical
      d
    }
    df <- remap_plane(df)

    # Optional faint background of ALL raw (pre-filter) tracks - shown in BOTH full
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
        rd <- remap_plane(rd)   # keep the faint raw tracks in the same plane as the data
        split(rd, rd$tid)
      } else NULL
    } else NULL

    # ── FULL-TRACKS MODE: skip ALL playback frame/bin/trail machinery ──────────
    # The static render (view_mode_main == "full") reads only sp_static / col_map /
    # tids - never the per-frame data. Building frames here is what made length,
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
      # Full tracks normally draws every filtered TID - the chunk is a playback
      # cost control. When the user asks for it, scope this view to the selected
      # chunk as well, using the same TIDs the animation would show. Analysis and
      # filtering are untouched either way: they always read all TIDs.
      df_full <- if (isTRUE(isolate(input$trk_chunk_scope_full)))
                   scope_playback(df, isolate(trk_playback_chunk_tids())) else df

      full_build <- withProgress(message = "Building track view…", value = 0.2, {
        df_static <- decimate_per_tid(df_full, max_pts = max_pts, tcol = "tim",
                                      extra_cols = c("col_dark", "n_dark",
                                                     "col_speed", "v_drift"))
        setProgress(0.5, detail = sprintf("%s tracks", format(length(unique(df_static$tid)), big.mark = ",")))
        sp_static <- split(df_static, df_static$tid)
        # (Z is already corrected in trk() by the physical Z factor - no display-only
        # stretch here, otherwise it would be applied twice.)
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
        using_bins = FALSE, step_dt = NA_real_, has_3d = use3d_f,
        plane      = plane,        # 2D viewing plane the coordinates were remapped to
        speed_max  = speed_max     # drift-speed colour-scale limit (NA unless speed mode)
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
    # Carry the per-point colour columns so trails can be coloured by the active mode.
    dwin_tr_dec  <- decimate_per_tid(dwin, max_pts = max_pts_tid, tcol = tcol_play,
                                     extra_cols = c("col_dark", "col_speed"))
    dwin_tr_dec  <- dwin_tr_dec[order(dwin_tr_dec$tid, dwin_tr_dec$t), ]
    
    # Split full-res for markers, and decimated for trails
    sp_full <- split(dwin[order(dwin$tid, dwin[[tcol_play]]), ], dwin$tid)  # FULL (markers)
    sp_tr   <- split(dwin_tr_dec, dwin_tr_dec$tid)                          # DECIMATED (trails)

    # Full-range (unwindowed) spatial data for the static full-tracks trace.
    # Using df (not dwin) so tracks outside the playback window still appear.
    df_static    <- decimate_per_tid(df, max_pts = max_pts_tid, tcol = "tim")
    sp_static    <- split(df_static, df_static$tid)

    # (Z is already corrected in trk() by the physical Z factor, so the per-TID
    # spatial lists inherit it - no separate display stretch is applied here.)

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

    # Per-point colour column for the active mode (NULL = one colour per TID).
    pcol <- switch(isolate(input$trk_color_mode %||% "tid"),
                   dark = "col_dark", speed = "col_speed", NULL)
    color_by_dark <- !is.null(pcol)   # any per-point mode: markers carry the colour
    # Display filter by integration status: "all" | "integrated" (n_dark==0) | "dark" (n_dark>=1)
    dark_show <- isolate(input$trk_dark_show %||% "all")

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
      # Integration-status display filter (markers only; trails keep the full path).
      if (dark_show != "all" && "n_dark" %in% names(tb)) {
        valid <- valid & (if (dark_show == "integrated") tb$n_dark == 0L else tb$n_dark >= 1L)
      }
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
          if (!is.null(pcol) && pcol %in% names(tb)) {
            mcols <- c(mcols, as.character(tb[[pcol]][keep]))
          } else {
            mcols <- c(mcols, rep(unname(as.character(col_map[[tids[i]]])), length(keep)))
          }
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
            # Per-point colour mode: colour the whole trail by its LEADING (most recent)
            # point. plotly 2D lines take a single colour per trace - per-vertex line
            # colour is not supported - but a trail only spans trail_len, so its leading
            # value is "the value right now" and the trail recolours as playback advances.
            if (!is.null(pcol) && pcol %in% names(tb_tr))
              tr_entry$col <- as.character(tb_tr[[pcol]][wtr[length(wtr)]])
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
      has_3d      = use3d,     # TRUE when z column is present and 3D mode is active
      plane       = plane,     # 2D viewing plane the coordinates were remapped to
      color_by_dark = color_by_dark,  # per-point colour mode: markers carry it (overlaid even in trails mode)
      speed_max   = speed_max         # robust upper limit for the drift-speed colour scale
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
    # Live frame + marker count (moved off the plot title: updating a title needs a
    # relayout, which re-solves scaleanchor and rescales the FOV every frame).
    k_now <- suppressWarnings(as.integer(tick()))
    cur   <- if (is.finite(k_now) && k_now >= 1L && k_now <= length(pre$frames))
               pre$frames[[k_now]] else NULL
    htmltools::div(
      style = "color:#555;margin-top:6px;margin-bottom:6px;",
      htmltools::tags$b("Time basis: "), tcol, htmltools::HTML("&nbsp;&middot;&nbsp;"),
      htmltools::tags$b("Frames: "), length(pre$frames), htmltools::HTML("&nbsp;&middot;&nbsp;"),
      htmltools::tags$b("Frame dt ~ "), fmt_num(dt_show), " ms",
      # Own line + nowrap so changing digit counts can't reflow the sidebar mid-playback.
      if (!is.null(cur)) htmltools::div(
        style = "white-space:nowrap;",
        htmltools::tags$b("Frame: "), k_now, " / ", length(pre$frames),
        htmltools::HTML("&nbsp;&middot;&nbsp;"),
        htmltools::tags$b("Markers: "), length(cur$markers_x)
      )
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
    # same defaults the render applies). The view_mode / z_scale / chunk_sel
    # observers compare against these and skip firing when unchanged - so their
    # INPUT CREATION after an upload (value == the just-rendered default) never
    # triggers a redundant rebuild, without any fragile time window.
    loop$rendered_view  <- isolate(input$trk_view_mode  %||% "2d")
    loop$rendered_zs    <- isolate(trk_z_applied())
    loop$rendered_chunk <- isolate(input$trk_chunk_sel  %||% "__all__")
    pre  <- isolate(frames_store()); req(is.list(pre), length(pre$frames) > 0)
    f1    <- pre$frames[[1]]
    tids  <- pre$tids
    cmap  <- pre$col_map
    sp    <- pre$sp
    use3d <- isTRUE(pre$has_3d)

    # 2D viewing plane the frames were built in (XY / XZ / YZ). The coordinates are
    # already remapped in frames_store, so here we only derive axis labels and decide
    # whether the XY-physical overlays (dwell heatmap, background image) still apply.
    plane     <- pre$plane %||% "XY"
    loop$rendered_plane <- plane
    plane_lab <- switch(plane,
      XZ = list(x = "X (nm)", y = "Z (nm)"),
      YZ = list(x = "Y (nm)", y = "Z (nm)"),
      list(x = "X (nm)", y = "Y (nm)"))
    plane_xy  <- identical(plane, "XY")

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

    # Dwell-time heatmap (2D only) - shared by FULL-tracks and PLAYBACK modes so the
    # toggle works in both. heatmap_store() is computed from trk() on plot_trigger
    # (not playback-specific), so it has data in either mode. Added as a background
    # layer (before the tracks) so paths draw on top.
    # The dwell-time heatmap is computed in real XY space, so it is only meaningful in
    # the XY plane - suppress it (like the background image) in the XZ/YZ side views.
    hm_data <- if (!use3d && plane_xy && isTRUE(isolate(input$trk_show_heatmap)))
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

    # ── Full-tracks static mode ───────────────────────────────────────────────
    if (identical(isolate(input$trk_view_mode_main), "full")) {
      connect  <- isTRUE(isolate(input$trk_full_connect))
      sp_all   <- if (!is.null(pre$sp_static)) pre$sp_static else sp
      all_tids <- names(sp_all)
      dark     <- isTRUE(isolate(input$trk_dark_mode))
      bg_col   <- if (dark) "#000000" else "#FFFFFF"
      fg_col   <- if (dark) "#FFFFFF" else "#63666A"
      cmode_st <- isolate(input$trk_color_mode %||% "tid")
      pcol_st  <- switch(cmode_st, dark = "col_dark", speed = "col_speed", NULL)
      color_by_dark <- !is.null(pcol_st) &&
                       any(vapply(sp_all, function(p) pcol_st %in% names(p), logical(1)))
      dark_line_col <- if (dark) "rgba(200,200,200,0.30)" else "rgba(120,120,120,0.30)"
      # Integration-status display filter: keep all / integrated (n_dark==0) / dark (>=1).
      dark_show <- isolate(input$trk_dark_show %||% "all")
      dark_keep <- function(pts) {
        if (dark_show == "all" || !("n_dark" %in% names(pts))) return(pts)
        m <- if (dark_show == "integrated") pts$n_dark == 0L else pts$n_dark >= 1L
        pts[m, , drop = FALSE]
      }

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
          pts <- dark_keep(pts)
          if (!nrow(pts)) next
          col_i    <- as.character(cmap[[nm]])
          tim_col  <- if ("tim" %in% names(pts)) pts$tim else pts$t
          hover_z  <- if (use3d && "z" %in% names(pts))
                        sprintf("<br>z: %.1f nm", pts$z) else
                        rep("", nrow(pts))
          hover_d  <- if ("n_dark" %in% names(pts))
                        sprintf("<br>dark cycles: %d", pts$n_dark) else rep("", nrow(pts))
          hover_txt <- sprintf("TID: %s<br>t: %.1f ms<br>x: %.1f nm<br>y: %.1f nm%s%s",
                               nm, tim_col, pts$x, pts$y, hover_z, hover_d)
          mk_col <- if (color_by_dark && pcol_st %in% names(pts)) as.character(pts[[pcol_st]]) else col_i
          ln_col <- if (color_by_dark) dark_line_col else col_i
          tr_args <- list(
            x = pts$x, y = pts$y,
            type = tr_type, mode = "lines+markers",
            line   = list(width = 1, color = ln_col),
            marker = list(size = isolate(input$trk_pt_size), color = mk_col),
            text = hover_txt, hovertemplate = "%{text}<extra></extra>",
            showlegend = (i == 1L && !color_by_dark), legendgroup = "Localizations",
            name = if (i == 1L) "Localizations" else paste0("tid_", nm)
          )
          if (use3d) tr_args$z <- pts$z
          fig <- do.call(add_trace, c(list(fig), tr_args))
        }
        # Connect mode draws one trace per TID with hex colours, so it gets no automatic
        # colourbar - append the scale explicitly.
        if (color_by_dark) {
          a1 <- sp_all[[all_tids[1]]]
          fig <- add_colorbar_trace(fig, cmode_st, fg_col,
                                    vmax = suppressWarnings(as.numeric(pre$speed_max)),
                                    use3d = use3d,
                                    ax = a1$x[1], ay = a1$y[1],
                                    az = if (use3d && "z" %in% names(a1)) a1$z[1] else 0)
        }
      } else {
        # Single fast scattergl trace with per-point colors + hover text
        tr_type <- if (use3d) "scatter3d" else "scattergl"
        x_v <- c(); y_v <- c(); z_v <- c(); col_v <- c(); nd_v <- c(); vd_v <- c(); txt_v <- c()
        for (nm in all_tids) {
          pts   <- sp_all[[nm]]
          if (!nrow(pts)) next
          pts <- dark_keep(pts)
          if (!nrow(pts)) next
          col_i   <- as.character(cmap[[nm]])
          tim_col <- if ("tim" %in% names(pts)) pts$tim else pts$t
          hover_z <- if (use3d && "z" %in% names(pts))
                       sprintf("<br>z: %.1f nm", pts$z) else
                       rep("", nrow(pts))
          hover_d <- if ("n_dark" %in% names(pts))
                       sprintf("<br>dark cycles: %d", pts$n_dark) else rep("", nrow(pts))
          hover_v <- if ("v_drift" %in% names(pts))
                       ifelse(is.finite(pts$v_drift),
                              sprintf("<br>drift speed: %.2f µm/s", pts$v_drift),
                              "<br>drift speed: n/a") else rep("", nrow(pts))
          txt_v <- c(txt_v, sprintf("TID: %s<br>t: %.1f ms<br>x: %.1f nm<br>y: %.1f nm%s%s%s",
                                    nm, tim_col, pts$x, pts$y, hover_z, hover_d, hover_v))
          x_v   <- c(x_v,   pts$x)
          y_v   <- c(y_v,   pts$y)
          if (use3d) z_v <- c(z_v, pts$z)
          col_v <- c(col_v, rep(col_i, nrow(pts)))
          if ("n_dark" %in% names(pts))  nd_v <- c(nd_v, pts$n_dark)
          if ("v_drift" %in% names(pts)) vd_v <- c(vd_v, pts$v_drift)
        }
        # Dark-cycle colouring uses a numeric colour array + colorscale so plotly draws
        # a self-documenting colorbar (0 = integrated/blue → deeper red = more dark cycles).
        if (identical(cmode_st, "dark") && color_by_dark && length(nd_v) == length(x_v)) {
          mk <- list(size = isolate(input$trk_pt_size), color = nd_v,
                     colorscale = dark_cycle_colorscale(), cmin = 0, cmax = DARK_NMAX,
                     showscale = TRUE,
                     colorbar = list(title = list(text = "dark<br>cycles", font = list(color = fg_col)),
                                     tickmode = "array", tickvals = 0:DARK_NMAX,
                                     ticktext = c(as.character(0:(DARK_NMAX - 1L)), paste0(DARK_NMAX, "+")),
                                     tickfont = list(color = fg_col), len = 0.5, thickness = 14))
          show_lg <- FALSE
        } else if (identical(cmode_st, "speed") && color_by_dark && length(vd_v) == length(x_v)) {
          vmax_st <- suppressWarnings(as.numeric(pre$speed_max))
          if (!is.finite(vmax_st) || vmax_st <= 0)
            vmax_st <- suppressWarnings(stats::quantile(vd_v, 0.95, na.rm = TRUE))
          if (!is.finite(vmax_st) || vmax_st <= 0) vmax_st <- 1
          mk <- list(size = isolate(input$trk_pt_size), color = vd_v,
                     colorscale = speed_colorscale(), cmin = 0, cmax = vmax_st,
                     showscale = TRUE,
                     colorbar = list(title = list(text = "drift speed<br>(µm/s)",
                                                 font = list(color = fg_col)),
                                     tickfont = list(color = fg_col), len = 0.5, thickness = 14))
          show_lg <- FALSE
        } else {
          mk <- list(size = isolate(input$trk_pt_size), color = col_v)
          show_lg <- TRUE
        }
        tr_args <- list(
          x = x_v, y = y_v,
          type = tr_type, mode = "markers",
          marker = mk,
          text = txt_v, hovertemplate = "%{text}<extra></extra>",
          showlegend = show_lg, name = "Localizations", legendgroup = "Localizations"
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
      ovl <- if (plane_xy) tryCatch(isolate(trk_image_overlay()), error = function(e) NULL) else NULL

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
          xaxis = c(list(title = if (hide_grid) "" else plane_lab$x, range = xr), axis_base),
          yaxis = c(list(title = if (hide_grid) "" else plane_lab$y, range = yr,
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
        # SMLM tab - reading inputs ISOLATED so toggling the checkbox does NOT invalidate
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
    # ── End full-tracks static mode ───────────────────────────────────────────

    # ---- Heatmap (optional, 2D only) - shared helper (see draw_heatmap above) ----
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

      # markers - empty when trails are shown (lines replace dots)
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

      # Colour scale for playback. Added AFTER trace_map is set so the extra trace can
      # never disturb trail_start / markers_idx (draw_frame restyles by index).
      # fg_col is local to the static branch (which returns above), so derive it here.
      cmode_pb <- isolate(input$trk_color_mode %||% "tid")
      if (cmode_pb %in% c("dark", "speed")) {
        fg_pb <- if (isTRUE(isolate(input$trk_dark_mode))) "#FFFFFF" else "#63666A"
        # Anchor on a real localization (frame 1, else the first track) - an all-NA
        # trace is dropped by plotly and would draw no colourbar.
        a1x <- f1$markers_x[1]; a1y <- f1$markers_y[1]
        a1z <- if (use3d) (f1$markers_z %||% numeric(0))[1] else 0
        if (!length(a1x) || !is.finite(a1x %||% NA_real_)) {
          fb  <- if (length(tids)) sp[[tids[1]]] else NULL
          a1x <- if (!is.null(fb)) fb$x[1] else 0
          a1y <- if (!is.null(fb)) fb$y[1] else 0
          a1z <- if (!is.null(fb) && use3d && "z" %in% names(fb)) fb$z[1] else 0
        }
        fig <- add_colorbar_trace(fig, cmode_pb, fg_pb,
                                  vmax = suppressWarnings(as.numeric(pre$speed_max)),
                                  use3d = use3d, ax = a1x, ay = a1y, az = a1z)
      }
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
      # first bin - using it would zoom the view into a tiny region and push the
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
        # 2D only - locked ranges and scale bar depend on these
        trk_axis_x_span(diff(xr) + 2 * xpad)
        trk_locked_ranges(list(x = c(xr[1]-xpad, xr[2]+xpad),
                               y = c(yr[1]-ypad, yr[2]+ypad)))
      }
    }

    ovl       <- if (plane_xy) tryCatch(isolate(trk_image_overlay()), error = function(e) NULL) else NULL
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
      # 3D scene layout - uses scene instead of xaxis/yaxis
      ax3 <- list(color = fg_col, gridcolor = grid_col,
                  showgrid = !hide_grid, zeroline = !hide_grid,
                  showticklabels = !hide_grid, backgroundcolor = bg_col)

      # Equal axis ranges: use the same helper as the SMLM 3D plot.
      # get_equal_limits() finds the largest span across X, Y, Z and applies it
      # to all three axes (centred on each axis's own midpoint) - identical nm/tick.
      all_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
      all_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
      all_z <- unlist(lapply(sp, `[[`, "z"), use.names = FALSE)
      lim3  <- get_equal_limits(all_x, all_y, all_z)

      # autorange = FALSE is required - Plotly 3D ignores range unless autorange is
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
        # No frame text in the title: updating it would need a relayout, which re-solves
        # scaleanchor and jitters the FOV. The live frame readout lives in trk_frame_meta.
        title = list(text = "", font = list(color = fg_col)),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font = list(color = fg_col),
        xaxis = c(
          if (have_bbox) list(title = if (hide_grid) "" else plane_lab$x,
                              range = c(xr[1]-xpad, xr[2]+xpad), autorange = FALSE)
          else           list(title = if (hide_grid) "" else plane_lab$x),
          axis_grid
        ),
        yaxis = c(
          if (have_bbox) list(title = if (hide_grid) "" else plane_lab$y,
                              range = c(yr[1]-ypad, yr[2]+ypad), autorange = FALSE,
                              scaleanchor = "x", scaleratio = 1)
          else           list(title = if (hide_grid) "" else plane_lab$y,
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
    # ---- Diagnostic: detected cluster centres (red +) and confocal blobs (green o) ----
    # Use type="scatter" (SVG), not scattergl - scattergl lacks circle-open / cross symbols.
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
      # Per-point colour mode: recolour each trail LINE from its leading point, so the
      # trails carry the same colour meaning the markers would (blue→red speed, or the
      # dark-cycle gradient) and no dots are needed on top.
      per_pt <- isTRUE(pre$color_by_dark)
      lc <- if (per_pt) vector("list", n) else NULL
      for (i in seq_len(n)) {
        td      <- if (i <= length(tpt)) tpt[[i]] else NULL
        xs[[i]] <- if (!is.null(td) && length(td$x)) td$x else NA_real_
        ys[[i]] <- if (!is.null(td) && length(td$y)) td$y else NA_real_
        if (is3d) zs[[i]] <- if (!is.null(td) && length(td$z)) td$z else NA_real_
        if (per_pt) {
          # Empty trail -> transparent (its x/y are NA so nothing is drawn anyway).
          lc[[i]] <- if (!is.null(td) && !is.null(td$col) && nzchar(td$col)) td$col else "rgba(0,0,0,0)"
        }
      }
      restyle_data <- list(x = xs, y = ys)
      if (is3d) restyle_data$z <- zs
      if (per_pt) restyle_data$`line.color` <- lc
      plotlyProxyInvoke(proxy, "restyle", restyle_data,
        seq(tm$trail_start, tm$trail_start + n - 1L))
      # Markers stay empty while trails are shown - the trails themselves now carry the
      # per-point colour, so overlaying dots would just double-draw the same information.
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

    # ----- NO RELAYOUT HERE, EVER -----
    # draw_frame is RESTYLE-ONLY in every path. Any relayout - even a title-only one -
    # makes plotly re-solve the scaleanchor aspect constraint and visibly rescales the
    # FOV (the sequential loop skips relayout for the same reason; see
    # run_sequential_loop_plain). Gating it on "not playing" was not enough: draw_frame
    # is also driven by the slider observer and by stale slider echoes that arrive after
    # playback ends (slider_guard is reset to 0 there), so the flags can read "paused"
    # while frames are still advancing - which produced a per-frame FOV jitter after
    # stop-then-replay. The frame/marker readout now lives in trk_frame_meta (driven by
    # tick()), which costs the plot nothing. uirevision keeps the view - and any user
    # zoom - stable across the restyles.

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

    # Trail (single trace in sequential mode) - only when the trail trace exists
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
  
  # ── heatmap_store() - dwell-time heatmap ────────────────────────────────────
  # GOAL: a 2D map of how much TIME molecules spent in each region of the field,
  # rather than how many localizations landed there. Time-weighting matters
  # because a slow/trapped molecule produces many localizations in one spot; a
  # raw count map would exaggerate it purely as a sampling artefact.
  # Bin count is user-set, or auto-scaled from the dataset size.
  heatmap_store <- eventReactive(plot_trigger(), {
    df <- isolate(trk()); req(nrow(df) > 0)
    # Auto bin count keeps the map readable across very different dataset sizes.
    bins <- if (isTruthy(isolate(input$trk_hm_bins))) isolate(input$trk_hm_bins) else auto_hm_bins(nrow(df))
    build_dwell_heatmap(df[, c("x","y","tid","tim")], nbins = bins)   # needs tim for dwell weighting
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

  # --- Δt histogram in analysis tab (richer version with stats) ---------------
  output$trk_dt_hist_jdd <- renderPlotly({
    req(isTRUE(input$trk_show_tsa), tsa_ran())
    dts  <- dt_vectors()
    vec  <- unlist(dts, use.names = FALSE)
    vec  <- vec[is.finite(vec) & vec > 0]
    shiny::validate(shiny::need(length(vec) > 0,
      "No Δt values - need at least 2 time points per TID."))
    n_bins  <- 80L    # only used by the no-cycle-time fallback binning below
    dark    <- isTRUE(isolate(input$trk_dark_mode))
    bg_col  <- if (dark) "#000000" else "#FFFFFF"
    fg_col  <- if (dark) "#FFFFFF" else "#333333"

    dt_med  <- median(vec, na.rm = TRUE)
    # Bin-FREE typical Δt (half-sample mode) - independent of the bin slider.
    dt_hsm  <- half_sample_mode(vec)
    rng     <- range(vec)                                 # x-extent, for the cycle markers

    # The app's best cycle-time estimate anchors the bins and the cycle markers. Read
    # reactively so both follow a cycle change.
    t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    # SW offset: peaks sit at SW + k·t_c, so the bin grid has to be shifted by it.
    # Prefer the fitted value; fall back to (typical Δt − t_c), which is the same
    # quantity read straight off the 1-cycle peak.
    cc_sw <- {
      cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
      s  <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else NA_real_
      if (!is.finite(s) && is.finite(dt_hsm) && is.finite(t_c)) max(0, dt_hsm - t_c) else s
    }
    if (!is.finite(cc_sw)) cc_sw <- 0

    # ── Binning: one bar per cycle count ──────────────────────────────────────
    # Δt is not continuous: it takes values SW + k·t_c for integer k. Fixed-width bins
    # alias against that comb - at 80 bins the grid runs ~1.63 bins per cycle, so each
    # peak lands at a different offset inside its bin and the bars drift visibly away
    # from the cycle markers, with empty bins interleaved. Instead the bin EDGES sit on
    # the classification boundaries SW + (k ± ½)·t_c, so every bar is exactly one cycle
    # count centred on its peak, and the picture is identical to the
    # round((Δt − SW)/t_c) the analysis performs. Falls back to fixed-width bins only
    # when no cycle time is available.
    use_cycle_bins <- is.finite(t_c) && t_c > 0
    if (use_cycle_bins) {
      kv    <- round((vec - cc_sw) / t_c)
      kv    <- kv[is.finite(kv) & kv >= 1]
      kmax  <- max(1L, min(as.integer(max(kv)), 400L))
      kcnt  <- tabulate(pmin(kv, kmax), nbins = kmax)
      bar_x <- cc_sw + seq_len(kmax) * t_c      # bar centre = the peak itself
      bar_y <- kcnt
      bar_w <- t_c * 0.9
      keepb <- bar_y > 0
      bar_x <- bar_x[keepb]; bar_y <- bar_y[keepb]
    } else {
      h_tmp <- hist(vec, breaks = n_bins, plot = FALSE)
      bar_x <- h_tmp$mids; bar_y <- h_tmp$counts
      bar_w <- diff(h_tmp$breaks)[1]
    }

    # Counts span ~100× between the 1-cycle peak and the dark tail, so the tail is
    # invisible on a linear axis; log y makes the whole decay readable.
    log_y <- isTRUE(input$trk_dt_log_y %||% FALSE)
    ymax  <- if (log_y) 10^(log10(max(bar_y, na.rm = TRUE)) * 1.06)
             else max(bar_y, na.rm = TRUE) * 1.08
    ybase <- if (log_y) 0.7 else 0        # log axes cannot start at 0

    fig <- plot_ly() %>%
      add_bars(x = bar_x, y = bar_y, width = bar_w,
               marker = list(color = "rgba(100,150,220,0.7)",
                             line  = list(color = "rgba(70,110,180,0.9)", width = 0.5)),
               hovertemplate = if (use_cycle_bins)
                 "Δt %{x:.4g} ms<br>%{y:,} steps<extra></extra>"
               else "Δt %{x:.4g} ms<br>%{y:,}<extra></extra>",
               showlegend = FALSE)
    fig <- fig %>%
      add_trace(x = c(dt_med, dt_med), y = c(ybase, ymax),
                type = "scatter", mode = "lines",
                line   = list(color = "#E03020", width = 2, dash = "dash"),
                name = sprintf("Median %.3g ms", dt_med)) %>%
      add_trace(x = c(dt_hsm, dt_hsm), y = c(ybase, ymax),
                type = "scatter", mode = "lines",
                line   = list(color = "#20A040", width = 2, dash = "dashdot"),
                name = sprintf("Typical Δt (bin-free) %.3g ms", dt_hsm))

    # ── Cycle-multiple markers ────────────────────────────────────────────────
    # A k-cycle interval sits at Δt = SW + k·t_c, i.e. modal + (k−1)·t_c (the peaks
    # are NOT at k·t_c - that omits the once-per-step SW overhead and lands every
    # marker short of its bell). So anchor marker 1 on the observed modal Δt (the
    # tallest = 1-cycle peak, = the green Typical-Δt line) and step by the SELECTED
    # cycle time t_c. Consequence: if the selected cycle time matches the data, all
    # markers land on the bells; if a manual/hardware override is wrong, the higher
    # markers visibly drift off - a built-in check that the chosen cycle time is right.
    # Layout shapes (not traces) so they stay out of the legend.
    shp <- list(); ann <- list()
    if (isTRUE(input$trk_dt_mark_cycles %||% TRUE) &&
        is.finite(t_c) && t_c > 0 && is.finite(dt_hsm)) {
      n_cap <- suppressWarnings(as.integer(input$trk_cycle_dark_max %||% 10L))
      if (!is.finite(n_cap) || n_cap < 1L) n_cap <- 10L
      Nmax  <- max(1L, min(n_cap, 12L, as.integer(floor((rng[2] - dt_hsm) / t_c)) + 1L))
      grid_c <- if (dark) "rgba(200,200,200,0.35)" else "rgba(90,90,90,0.35)"
      for (nn in seq_len(Nmax)) {
        xc <- dt_hsm + (nn - 1L) * t_c        # SW + nn·t_c, anchored on the modal peak
        shp[[length(shp) + 1L]] <- list(type = "line", x0 = xc, x1 = xc,
          y0 = ybase, y1 = ymax,
          line = list(color = grid_c, width = 1, dash = "dot"), layer = "below")
        ann[[length(ann) + 1L]] <- list(x = xc, y = ymax, text = paste0(nn, "×"),
          showarrow = FALSE, yanchor = "bottom",
          font = list(size = 9, color = fg_col))
      }
    }

    fig %>% layout(
      xaxis  = list(title = "Δt between consecutive localizations (ms)",
                    color = fg_col),
      yaxis  = if (log_y)
                 list(title = "Count (log)", color = fg_col, type = "log",
                      range = c(log10(ybase), log10(ymax)))
               else list(title = "Count", color = fg_col, range = c(0, ymax)),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font   = list(color = fg_col),
      # Vertical legend to the RIGHT of the plot (x > 1 = outside the axes) so it no
      # longer overlaps the bars / cycle labels along the top.
      legend = list(orientation = "v", x = 1.02, xanchor = "left",
                    y = 1, yanchor = "top", font = list(color = fg_col)),
      margin = list(r = 10),
      bargap = 0.05,
      shapes = shp,
      annotations = ann,
      uirevision = "dt_jdd"
    )
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
    # Peak-based cycle time: the Δt peaks sit at k·cycle + SW, so a least-squares fit
    # of the peak centres against the integer cycle index k gives slope = cycle time,
    # intercept = SW. This is the "fit the 1×,2×,3×… peaks and read the cycle time off
    # the spacing" idea, and is the same quantity Bela/Vogler recover from the EOD/SW
    # time-signature fit. Computed here regardless of the acquisition-tab toggle so it
    # is always visible next to the density peaks it comes from.
    rc <- tryCatch(trk_regression_cycle(), error = function(e) NULL)
    cc <- tryCatch(trk_cluster_cycle(),    error = function(e) NULL)
    tags$div(style = "font-size:0.84em; color:#555; margin:4px 0;",
      tags$div(
        tags$b("Typical Δt (bin-free, half-sample mode): "),
        tags$span(style = "color:#20803a;", sprintf("%.4g ms", dt_hsm)),
        sprintf("   ·   median %.4g   ·   p10 %.4g   ·   p25 %.4g   ·   p75 %.4g   ·   n = %s",
                dt_med, dt_p10, dt_p25, dt_p75, format(length(vec), big.mark = ","))
      ),
      if (!is.null(rc) && is.finite(rc$cycle) && rc$cycle > 0)
        tags$div(style = "color:#3a5a8a; margin-top:2px;",
          tags$b("Peak-fit cycle time: "),
          tags$span(style = "font-weight:bold;", sprintf("%.4g ms", rc$cycle)),
          sprintf("   ·   SW offset %.3g ms   ·   R² = %.4f   ·   fit on %d cycle peaks",
                  rc$sw, rc$r2 %||% NA_real_, rc$n),
          if (!is.null(cc) && is.finite(cc$frac_one_cycle))
            sprintf("   ·   %.0f%% are single (integrated) cycles", 100 * cc$frac_one_cycle) else "",
          # The cycle time is smaller than every gap you can actually observe, which
          # reads as a discrepancy unless the decomposition is spelled out.
          if (is.finite(rc$cycle) && is.finite(rc$sw))
            tags$div(style = "color:#666; font-weight:normal; margin-top:1px;",
              sprintf("Cycle time %.4g + SW %.3g = %.4g ms, the typical \u0394t between consecutive localizations. The cycle time itself is the repeat unit, not an observable gap.",
                      rc$cycle, rc$sw, rc$cycle + rc$sw))
          else NULL)
      else if (!is.null(cc) && is.finite(cc$cycle) && cc$cycle > 0)
        tags$div(style = "color:#3a5a8a; margin-top:2px;",
          tags$b("Peak-spacing cycle time: "),
          tags$span(style = "font-weight:bold;", sprintf("%.4g ms", cc$cycle)),
          sprintf("   ·   from %d Δt peaks (need ≥2 for the regression fit)", length(cc$centers %||% numeric(0)))),
      if (is.finite(hw) && hw > 0)
        tags$div(style = "color:#888; margin-top:2px;",
          sprintf("Hardware cycle = %.4g ms. A typical \u0394t near twice this means dark cycles dominate.", hw)),
      tags$div(style = "color:#888; margin-top:2px; font-size:0.95em;",
        "This estimate is independent of the histogram bin slider (which is display-only).",
        " The peak-fit cycle time regresses the Δt density peaks (1×, 2×, 3×… the cycle) against",
        " their cycle index - enable ", tags$em("Fit cycle time from data (EOD/SW regression)"),
        " in Acquisition Parameters to use it everywhere.")
    )
  })


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
      mode               = "MSD (cumulative cycle-index binning)",
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
      msd_formula        = sprintf("MSD = %.4gD*tau + %.4g*sigma^2 (d_eff = %.4g)",
                                   ub$coef_eff %||% (2L * ub$ndim),
                                   ub$coef_eff %||% (2L * ub$ndim),
                                   ub$d_eff %||% ub$ndim)
    )

    list(data = data, params = params)
  })
  # Per-track D and alpha fits - computed when histogram OR clustering is open
  # Per-track fits feeding the D histogram, the D/CI clustering and the per-track
  # CSV export - all derived from the SAME Unified MSD calculation (per-track OLSF
  # D + confinement index), via the cached trk_unified_tdiff reactive.
  msd_pertk_fits <- reactive({
    req(isTRUE(input$trk_show_msd))
    req(isTRUE(input$trk_msd_show_hist) || isTRUE(input$trk_msd_cluster))
    td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(td) && !is.null(td$pertrack) && nrow(td$pertrack) > 0,
      "Enable the MSD and click 'Plot MSD' to compute per-track diffusion."))
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

  # ── Per-cluster motion model ────────────────────────────────────────────────
  # The k-means split above groups tracks by D and confinement index, but it does not
  # say WHAT KIND of motion each group shows. Here each cluster's own ensemble MSD
  # (mean of its members' TAMSDs) is fitted with all four models and the best one is
  # selected by AICc, so a cluster can be reported as free / confined / directed /
  # anomalous rather than just "slow" or "fast".
  #
  # AICc (not raw residuals) because the models differ in parameter count: linear has
  # 2, the others 3. n here is the number of MSD lag points in the fit window, which is
  # small, so the small-sample correction matters. Lags are correlated, so this is a
  # ranking heuristic, not a likelihood test - reported as such.
  msd_cluster_models <- reactive({
    cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
    if (is.null(cl) || is.null(ub) || is.null(ub$track_msd)) return(NULL)
    min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
    max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    tau_all <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
    win <- ub$lag >= min_lag & ub$lag <= max_lag & ub$keep_lag & is.finite(tau_all)
    if (sum(win) < 5L) return(NULL)

    ids <- cl$cluster_ids
    rows <- lapply(seq_len(cl$n_clusters), function(cc) {
      tids_c <- names(ids)[ids == cc]
      cols   <- which(colnames(ub$track_msd) %in% tids_c)
      if (!length(cols)) return(NULL)
      # Cluster ensemble MSD = mean of member TAMSDs (same convention as the eMSD).
      y <- rowMeans(ub$track_msd[, cols, drop = FALSE], na.rm = TRUE)[win]
      s <- tau_all[win]
      okp <- is.finite(s) & is.finite(y) & y > 0
      if (sum(okp) < 5L) return(NULL)
      s <- s[okp]; y <- y[okp]; np <- ub$npairs[win][okp]

      aicc_of <- function(fr, k_par) {
        if (is.null(fr) || is.null(fr$tau_fit)) return(NA_real_)
        yhat <- stats::approx(fr$tau_fit, fr$msd_fit, xout = s, rule = 2)$y
        rss  <- sum((y - yhat)^2); n <- length(y)
        if (!is.finite(rss) || rss <= 0 || n <= k_par + 2L) return(NA_real_)
        n * log(rss / n) + 2 * k_par + (2 * k_par * (k_par + 1)) / (n - k_par - 1)
      }
      # The free model is fitted here by plain OLS over the SAME window as the other
      # three, NOT by the OLSF used for the headline D. OLSF deliberately fits only the
      # first p (Michalet-optimal) lags, so its residuals over the full window are not
      # comparable with full-window non-linear fits - scoring it that way made every
      # free cluster look like a 3-parameter model (free won 0/10 on simulated free
      # curves; 9/10 once fitted on the same window). The headline D is unaffected;
      # this OLS is used only to rank the models against each other.
      cv <- 2L * ub$ndim
      free_fr <- tryCatch({
        lf <- stats::lm(y ~ s)
        list(D = unname(stats::coef(lf)[2]) / cv,
             sigma_loc = sqrt(max(0, unname(stats::coef(lf)[1])) / cv),
             tau_fit = s, msd_fit = unname(stats::fitted(lf)))
      }, error = function(e) NULL)
      cand <- list(
        free      = list(fr = free_fr, k = 2L),
        confined  = list(fr = tryCatch(.msd_model_fit(s, y, np, ub$ndim, "confined",
                                        FALSE, ub$t_cycle, coef_override = ub$coef_eff, R = trk_blur()$R_msd),
                                       error = function(e) NULL), k = 3L),
        directed  = list(fr = tryCatch(.msd_model_fit(s, y, np, ub$ndim, "directed",
                                        FALSE, ub$t_cycle, coef_override = ub$coef_eff, R = trk_blur()$R_msd),
                                       error = function(e) NULL), k = 3L),
        anomalous = list(fr = tryCatch(.msd_model_fit(s, y, np, ub$ndim, "anomalous",
                                        FALSE, ub$t_cycle, coef_override = ub$coef_eff, R = trk_blur()$R_msd),
                                       error = function(e) NULL), k = 3L))
      aic <- vapply(names(cand), function(nm) aicc_of(cand[[nm]]$fr, cand[[nm]]$k), numeric(1))
      if (!any(is.finite(aic))) return(NULL)
      auto_best <- names(cand)[which.min(replace(aic, !is.finite(aic), Inf))]
      # Per-cluster override: the whole point is to be able to say "Cluster 1 is free,
      # Cluster 2 is confined" and see both fits, so an explicit choice wins over AICc.
      # "auto" (default) keeps the AICc pick.
      pick <- input[[paste0("trk_cmodel_", cc)]] %||% "auto"
      if (!identical(pick, "auto") && !is.null(cand[[pick]]$fr)) best <- pick else best <- auto_best
      bfr  <- cand[[best]]$fr
      forced <- !identical(pick, "auto") && !identical(best, auto_best)
      # A confined fit that never reached its plateau is really just a line.
      if (identical(best, "confined") && isTRUE(bfr$degenerate)) best <- "free (no plateau)"
      detail <- switch(sub(" .*", "", best),
        confined  = sprintf("L_c = %.0f nm", bfr$L_conf),
        directed  = sprintf("v = %.3g µm/s", bfr$v_um_s),
        anomalous = sprintf("α = %.2f", bfr$alpha),
        "")
      d_aic <- sort(aic[is.finite(aic)])
      list(
        row = data.frame(
          cluster = cc, n_tracks = length(tids_c),
          model = best, auto_model = auto_best, forced = forced,
          D_um2_s = if (is.finite(bfr$D)) bfr$D / 1e3 else NA_real_,
          detail = detail,
          # Gap to the runner-up: small = the models are not really distinguishable here.
          dAICc = if (length(d_aic) >= 2L) unname(d_aic[2] - d_aic[1]) else NA_real_,
          stringsAsFactors = FALSE),
        # Observed cluster eMSD + the fitted curve, for the overlay plot (µm², ms).
        curve = list(cluster = cc, tau = s, msd = y,
                     fit_tau = bfr$tau_fit, fit_msd = bfr$msd_fit))
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(NULL)
    list(table  = do.call(rbind, lapply(rows, `[[`, "row")),
         curves = lapply(rows, `[[`, "curve"))
  })

  # One model selector per cluster. Auto = the AICc pick; anything else forces that
  # model on that cluster, so different populations can be modelled differently.
  output$trk_msd_cluster_model_ui <- renderUI({
    req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
    cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
    if (is.null(cl)) return(NULL)
    ch <- c("Auto (best AICc)" = "auto", "Free / Brownian" = "free",
            "Confined" = "confined", "Directed" = "directed", "Anomalous" = "anomalous")
    tags$div(style = "margin-top:4px;",
      tags$div(style = "font-size:0.85em; color:#555; margin-bottom:2px;",
               "Model per cluster"),
      lapply(seq_len(cl$n_clusters), function(cc) {
        div(style = "display:inline-block; width:49%; padding-right:4px; vertical-align:top;",
            selectInput(paste0("trk_cmodel_", cc),
                        label = HTML(sprintf(
                          "<span style='color:%s;font-weight:bold;'>Cluster %d</span>",
                          cl$palette[cc], cc)),
                        choices = ch, selected = "auto", width = "100%"))
      }))
  })

  output$trk_msd_cluster_models <- renderUI({
    req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
    res <- tryCatch(msd_cluster_models(), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    tb  <- res$table
    pal <- tryCatch(msd_cluster_rv()$palette, error = function(e) NULL)
    tags$div(style = "font-size:0.85em; margin-top:6px;",
      tags$b("Motion model per cluster"),
      tags$table(style = "width:100%; border-collapse:collapse; margin-top:3px;",
        tags$tr(style = "color:#666; text-align:left;",
          tags$th("Cluster"), tags$th("Tracks"), tags$th("Model"),
          tags$th("D (µm²/s)"), tags$th(""), tags$th("ΔAICc")),
        lapply(seq_len(nrow(tb)), function(i) {
          r <- tb[i, ]
          tags$tr(style = "border-top:1px solid #eee;",
            tags$td(tags$span(style = sprintf("color:%s; font-weight:bold;",
                                              if (!is.null(pal)) pal[r$cluster] else "#333"),
                              sprintf("%d", r$cluster))),
            tags$td(r$n_tracks),
            tags$td(tags$b(r$model),
                    # Make it obvious when a forced choice differs from the AICc pick.
                    if (isTRUE(r$forced))
                      tags$span(style = "color:#b26a00;",
                                sprintf(" (set; AICc: %s)", r$auto_model)) else NULL),
            tags$td(if (is.finite(r$D_um2_s)) sprintf("%.3g", r$D_um2_s) else "-"),
            tags$td(r$detail),
            tags$td(if (is.finite(r$dAICc)) sprintf("%.1f", r$dAICc) else "-"))
        })),
      tags$div(style = "color:#888; margin-top:3px;",
        "Each cluster's ensemble MSD is fitted with every model; Auto picks the lowest AICc. ",
        tags$b("\u0394AICc"), " below about 2: the two best models are not distinguishable. It ranks models; it is not a test.")
    )
  })

  # Cluster eMSD curves with each cluster's fitted model overlaid.
  output$trk_msd_cluster_curves <- renderPlotly({
    req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
    res <- tryCatch(msd_cluster_models(), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(res), "Click 'Plot MSD' to compute cluster models."))
    cl   <- msd_cluster_rv()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    tb <- res$table
    fig <- plot_ly()
    for (cv in res$curves) {
      col <- cl$palette[cv$cluster]
      lbl <- tb$model[match(cv$cluster, tb$cluster)]
      fig <- fig %>%
        add_trace(x = cv$tau, y = cv$msd / 1e6, type = "scatter", mode = "markers",
                  marker = list(size = 5, color = col),
                  name = sprintf("Cluster %d eMSD", cv$cluster))
      if (!is.null(cv$fit_tau))
        fig <- fig %>%
          add_trace(x = cv$fit_tau, y = pmax(0, cv$fit_msd) / 1e6,
                    type = "scatter", mode = "lines",
                    line = list(color = col, width = 2, dash = "dash"),
                    name = sprintf("Cluster %d: %s", cv$cluster, lbl))
    }
    fig %>% layout(
      xaxis = list(title = "τ (ms)", color = fg),
      yaxis = list(title = "MSD (µm²)", color = fg, rangemode = "tozero"),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(font = list(color = fg, size = 10), orientation = "h",
                    x = 0, xanchor = "left", y = -0.34, yanchor = "top"),
      margin = list(l = 55, r = 20, t = 10, b = 110),
      uirevision = "msd_cluster_curves")
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
  # Splits each track into diffusion-state segments via PELT on log(r²/Δt),
  # then fits a GMM (mclust, BIC-selected G) to the per-segment D distribution.
  # ---------------------------------------------------------------------------
  # ==========================================================================
  # cpd_segments_rv() + cpd_gmm_rv() - CHANGE-POINT SEGMENTATION + STATE COUNT
  # ==========================================================================
  #
  # GOAL
  #   Find the MOMENTS a molecule changes its mobility, without assuming in
  #   advance how many states exist. This is the model-free counterpart to the
  #   HMM: the HMM fits K states globally, whereas this first detects abrupt
  #   changes within each track and only then asks how many groups the
  #   resulting segments fall into.
  #
  # THE MATH
  #   (a) THE SIGNAL. From <r²> = 2dDΔt (header §2b), the per-step quantity
  #                        r² / Δt   ∝   local D
  #       is an (extremely noisy) instantaneous diffusivity. Dividing by Δt is
  #       what makes irregular MINFLUX sampling comparable - this is the
  #       SQUARED-displacement form of the sqrt(Δt) rule.
  #
  #       We segment log(r²/Δt), not r²/Δt, because:
  #         - single-step D is roughly exponential and strongly right-skewed;
  #           the log makes it ~Gaussian, which a mean-shift test assumes;
  #         - a change in D is MULTIPLICATIVE (e.g. "10x slower"), and in log
  #           space that becomes an ADDITIVE mean shift - what PELT detects.
  #
  #   (b) PELT (Killick, Fearnhead & Eckley 2012) finds the change-point set
  #       minimising  Σ_segments cost(segment) + β·(n change-points),
  #       exactly and in O(n). The BIC penalty β stops it inventing a change for
  #       every noisy step; `minseglen` sets the shortest allowed segment.
  #
  #   (c) PER-SEGMENT D uses the mean of ratios, NOT the ratio of means:
  #            D = mean(r²/Δt) / (2d)          (correct)
  #        not D = mean(r²) / (2d·mean(Δt))    (biased when r² and Δt co-vary -
  #                                             Jensen's inequality)
  #       Then subtract the noise floor at the D level:
  #            D_true = max(0, D − σ²·mean(1/Δt))
  #       so an immobile segment reads ~0 rather than a small positive D.
  #
  #   (d) HOW MANY STATES: pool the per-segment D values and fit a Gaussian
  #       mixture in log10(D), component count chosen by BIC (`cpd_gmm_rv`).
  #       Exploratory: segments of different lengths from one population can
  #       return two components (est_cpd.R header).
  #
  # OBJECTS IT CALLS (inputs)
  #   trk_spatial_filtered() (last pipeline stage), has_z(), trk_sigma_loc_eff()
  #   inputs: trk_cpd_compute (Run trigger), trk_cpd_min_seg, trk_cpd_max_g,
  #           trk_cpd_sub_noise
  #
  # STEP BY STEP
  #   1. Per track build r², Δt, and log(r²/Δt).
  #   2. Run PELT -> change-point indices -> segment boundaries.
  #   3. Per segment: D via mean-of-ratios, minus the noise floor.
  #   4. Pool all segments; fit the BIC-selected mixture in log10 D.
  #
  # RETURNS
  #   cpd_segments_rv : data.frame(tid, seg_idx, n_states, n_steps,
  #                                D_nm2_ms, D_um2_s) - one row per segment
  #   cpd_gmm_rv      : fitted mclust object; its components ARE the inferred
  #                     diffusive states
  #
  # CAVEAT: single-step D is very noisy, so PELT can over-segment (`min_seg` is
  # the guard). Cross-check the state count against the HMM - same question,
  # different route, they should broadly agree.
  # ==========================================================================
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

    # Optional σ_loc noise-floor subtraction on the per-segment D. A pure-noise step
    # has E[|Δr|²] = coef_D·σ² (each of the two localizations contributes σ² per axis),
    # so the apparent D of an immobile segment is σ²·⟨1/Δt⟩. Subtract it at the D level
    # (not per step - that would rectify and bias upward), clamped ≥ 0. Same floor the
    # D-map uses, so immobile segments no longer read as a small non-zero D. No-op when
    # σ_loc is unset or the toggle is off.
    sig_cpd    <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    sub_noise_cpd <- isTRUE(input$trk_cpd_sub_noise) && is.finite(sig_cpd) && sig_cpd > 0
    noise_D_off   <- if (sub_noise_cpd) sig_cpd^2 else 0   # × ⟨1/Δt⟩ per segment below

    t_c_cpd <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    # Motion blur needs the cycle time. Without it the fit would silently run with no
    # blur correction, which reads D about a third low at MINFLUX cycle times, so the
    # panel asks for the calibration instead of guessing.
    shiny::validate(shiny::need(is.finite(t_c_cpd) && t_c_cpd > 0,
      "Cycle time unknown - run Data Quality & Calibration first (the blur correction needs it)."))
    # Segmentation core lives in est_cpd.R (pure, tested).
    result <- withProgress(message = "Segmenting tracks…", value = 0,
      cpd_segment_tracks(df, min_seg = min_seg, use_z = use_z_cpd,
                         noise_sig2 = noise_D_off,
                         R_abs = trk_blur()$R_abs,
                         progress = function(f) incProgress(f)))
    shiny::validate(shiny::need(
      !is.null(result) && nrow(result) >= 3L,
      "No segments found - try lowering Min segment length."
    ))
    result
  }, ignoreNULL = TRUE)

  # ── cpd_gmm_rv() - how many diffusive states are there? ─────────────────────
  # GOAL: take the per-segment diffusion coefficients produced by the
  # change-point step and ask whether they fall into distinct populations.
  # A Gaussian mixture is fitted and the number of components chosen by BIC, so
  # the state count is selected by the data rather than asserted by the user.
  #
  # Fitted in log10(D) space: diffusion coefficients span orders of magnitude and
  # are right-skewed, so log space is where the populations look Gaussian.
  cpd_gmm_rv <- reactive({
    req(isTRUE(input$trk_show_cpd))
    segs <- tryCatch(cpd_segments_rv(), error = function(e) NULL)   # per-segment D values
    req(!is.null(segs), nrow(segs) >= 3L)
    cpd_gmm(segs, max_G = input$trk_cpd_max_g %||% 4L)              # est_cpd.R
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

    # add_histogram, not inline plot_ly(type="histogram", nbinsx=): the GMM overlay
    # scatter lines added below would otherwise inherit nbinsx from the plot_ly level
    # and warn "'scatter' objects don't have these attributes: 'nbinsx'".
    fig <- plot_ly() %>% add_histogram(
      x    = log10_D_fin,
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
      scale  <- length(log10_D_fin) * bw     # histogram → density scale factor

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

  # MSD fit core (blur R from trk_blur(), .olsf_sig2, .michalet_p, msd_fit_olsf,
  # msd_fit_olsf_iterative) lives in est_msd.R since 2026-09-17.

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
        # (MSD lags are correlated) - surfaced in the UI with that caveat.
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

  # Shared fit-window selector - applied consistently by the ensemble estimates,
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

  # Shared trackable-diffusion ceiling (D_max) estimate, backing the D_max panel.
  # Above ~D_max the MINFLUX localization loop can no longer keep the doughnut centred,
  # so D is systematically UNDER-estimated (Vogler et al. 2025; SimuFLUX, Marin & Ries 2025).
  #
  # WHICH L THIS IS. L_TCP is the diameter of the FINAL iteration - the one that is
  # reiterated continuously while tracking - because that is the pattern the particle
  # has to stay inside between consecutive fixes. The coarse search iterations that
  # precede lock-on are irrelevant here, and so is the much tighter final iteration an
  # IMAGING sequence would use: SPT deliberately stops at a larger L, since precision
  # improves only linearly as L shrinks (sigma_CRLB ~ L/(4*sqrt(N))) while the escape
  # risk grows quadratically (D_max ~ L^2). 150 nm is the tracking-iteration diameter
  # in Vogler et al. 2025 and matches the simulator's L default (sim_server.R).
  #
  # Pure function of the localization table + L_TCP (nm) & hardware dead-time t_hw (ms).
  # Returns D_max in um^2/s (Vogler form preferred, simple form fallback).
  dmax_estimate <- function(df, L = DMAX_L_TCP_DEFAULT, t_hw = 0.05) {
    if (is.null(df) || !nrow(df) || !all(c("tid","tim") %in% names(df))) return(NULL)
    if (!is.finite(L)    || L    <= 0) L    <- DMAX_L_TCP_DEFAULT
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

    L    <- suppressWarnings(as.numeric(input$trk_L_tcp))
    t_hw <- suppressWarnings(as.numeric(input$trk_dmax_t_hw))
    if (!is.finite(L)    || L    <= 0) L    <- DMAX_L_TCP_DEFAULT
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
      if (is.finite(x)) sprintf("%.*g", digits, x) else "-"

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
    # Short verdicts - the mechanism (the loop loses lock above D_max, so D reads low)
    # is in README "Fast-Diffusion Ceiling" and the panel footnote below.
    flag_msg <- function(ratio) {
      if (!is.finite(ratio)) return("")
      if (ratio > 0.8)  return(" ⚠ at the ceiling - D is likely under-estimated")
      if (ratio > 0.5)  return(" ⚠ approaching the ceiling")
      return(" ✓ safe")
    }

    use_ratio <- if (is.finite(ratio_vogler)) ratio_vogler else ratio_simple
    main_dmax <- if (is.finite(D_max_vogler)) D_max_vogler else D_max_simple

    htmltools::div(
      style = "font-size:0.88em; color:#333; padding:6px 0;",

      htmltools::tags$table(
        style = "width:100%; border-collapse:collapse;",
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "L_TCP (final iteration)"),
          htmltools::tags$td(style="padding:3px 6px; font-weight:bold;", sprintf("%g nm", L),
            # D_max goes as L^2, so the wrong iteration moves the ceiling by an order
            # of magnitude silently. Only flag clearly imaging-range diameters: fine
            # patterns are legitimate for slow targets, so the threshold sits well
            # below the 100-300 nm a tracking iteration normally uses.
            if (is.finite(L) && L < 60) htmltools::span(style = "font-weight:normal; color:#B8860B;",
              htmltools::HTML(" &#9888; tight for tracking. Confirm this is the reiterated iteration."))
            else NULL)
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "t_loc from ECO/EFO"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(t_loc_eco)) sprintf("%.3f ms", t_loc_eco) else "- (needs ECO & EFO)")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median Δt within tracks"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(t_loc_dt)) sprintf("%.3f ms", t_loc_dt) else "-")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median EFO (ψ_MFX)"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(efo_med)) sprintf("%.0f Hz", efo_med) else "- (no EFO column)")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median ECO (N_PL, photons/loc)"),
          htmltools::tags$td(style="padding:3px 6px;",
            if(is.finite(n_pl)) sprintf("%.1f", n_pl) else "- (no ECO column)")
        ),
        htmltools::tags$tr(style="border-top:1px solid #ddd;",
          htmltools::tags$td(style="padding:4px 6px; color:#666;", htmltools::HTML("D<sub>max</sub> simple (L_TCP² / 16·t_loc)")),
          htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_max_simple)) sprintf("%.3g µm²/s", D_max_simple) else "-")
        ),
        htmltools::tags$tr(
          htmltools::tags$td(style="padding:4px 6px; color:#666;", htmltools::HTML("D<sub>max</sub> Vogler (EFO/N_PL)")),
          htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_max_vogler)) sprintf("%.3g µm²/s", D_max_vogler) else "- (needs ECO & EFO)")
        ),
        htmltools::tags$tr(style="border-top:1px solid #ddd;",
          htmltools::tags$td(style="padding:4px 6px; color:#666;", "Measured D (MSD fit)"),
          htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_meas)) sprintf("%.3g µm²/s", D_meas) else "- (run MSD first)")
        ),
        if (is.finite(use_ratio)) htmltools::tags$tr(
          htmltools::tags$td(style="padding:4px 6px; color:#666;", htmltools::HTML("D / D<sub>max</sub>")),
          htmltools::tags$td(
            style = sprintf("padding:4px 6px; font-weight:bold; color:%s;", flag_col(use_ratio)),
            htmltools::HTML(sprintf("%.4f%s", use_ratio, flag_msg(use_ratio)))
          )
        )
      ),

      htmltools::tags$div(
        style = "margin-top:8px; font-size:0.82em; color:#666;",
        htmltools::HTML("L<sub>TCP</sub>²/(16·t_loc) &nbsp;·&nbsp; Vogler: (L<sub>TCP</sub>²/16)·ψ/(N_PL + ψ·t_hw), "),
        "with L the final iteration reiterated during tracking. The ceiling scales as L². ",
        htmltools::tags$b("Theoretical ceiling"), " (Vogler 2025 report ≈10–11 µm²/s); tracking is only ",
        "demonstrated to ≈2.5 µm²/s, so treat D above ~2–3 µm²/s as suspect even on a green flag."
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


  # ============================================================================
  # msd_unified_bundle() - THE MSD ENGINE
  # ============================================================================
  #
  # GOAL
  #   Build the mean-squared-displacement curve, MSD(τ), from which the
  #   diffusion coefficient D is fitted. This is the single MSD estimator in the
  #   app; every per-track D, the D histogram, the clustering and the CSV export
  #   all derive from this one bundle.
  #
  #   MSD(τ) = < |r(t+τ) − r(t)|² >  averaged over all pairs separated by τ.
  #   For free Brownian motion MSD = 2·d·D·τ + 2·d·σ_loc², so the SLOPE gives D
  #   and the INTERCEPT gives the localization error.
  #
  # THE MINFLUX PROBLEM THIS SOLVES
  #   MSD needs pairs grouped by a common time lag τ. MINFLUX has no fixed clock:
  #   consecutive localizations are separated by integer multiples of the base
  #   cycle, with "dark" cycles skipped at random. Binning by STEP INDEX would
  #   put physically different Δt into the same bin and corrupt the curve.
  #   Fix: bin by CUMULATIVE CYCLE INDEX (Vogler et al. 2025) - see step 3.
  #
  # OBJECTS IT CALLS (inputs)
  #   trk() ................. the analysis-ready dataset (see its header)
  #   trk_cycle_time() ...... base cycle time t_c, ms
  #   trk_cluster_cycle()$sw  software/settling overhead subtracted before
  #                           converting an interval into a cycle count
  #   .track_positive_D() ... optional per-track D>0 quality filter
  #   inputs: trk_msd_compute (the Run trigger), trk_uni_max_lag,
  #           trk_msd_min_pts, trk_msd_ndim, trk_uni_drop_negd,
  #           trk_uni_track_weight
  #
  # STEP BY STEP
  #   1. Resolve the cycle time t_c and the max cycle-lag to compute.
  #   2. Take trk(), pick 2D/3D, drop non-finite rows, split into one
  #      data.frame per track (`sp`).
  #   3. For each track, convert its irregular timestamps into a CUMULATIVE
  #      CYCLE INDEX:
  #         dcyc   = round((Δt − sw) / t_c)   cycles spanned by each step (>= 1)
  #         cycidx = cumsum(dcyc)             position on a linearized timeline
  #      The lag between two localizations is then the INTEGER difference of
  #      their cycidx - exact regardless of how many cycles were skipped.
  #      (Per-interval rounding then cumsum is deliberate: it is robust to a
  #      small error in t_c, whereas round(total_Δt / t_c) accumulates that
  #      error and produces a sawtooth in the curve.)
  #   4. For every pair (i, j) within a track, accumulate into three
  #      lag × track matrices:
  #         tk_dr2  Σ squared displacement    tk_np  pair count
  #         tk_dt   Σ real elapsed Δt  (kept so the x-axis can be TRUE mean Δt,
  #                                     not an assumed multiple of t_c)
  #   5. Per-track TAMSD = tk_dr2 / tk_np. Optionally drop tracks whose TAMSD
  #      slope is <= 0 (non-diffusing / artefact).
  #   6. Collapse to the ensemble curve, two weightings available:
  #         pair-weighted  Σdr² / Σn   (default; smoother, long tracks dominate)
  #         track-weighted mean of per-track TAMSDs (each track counts once -
  #                        this is Vogler's eMSD convention)
  #   7. Compute the IQR band across tracks and drop lags supported by fewer
  #      than 10% of the tracks (`keep_lag`) - the long-lag tail is built from
  #      very few, highly correlated pairs and is not trustworthy.
  #
  # RETURNS a list (this is the contract every consumer relies on):
  #   lag ........... integer cycle lags 1..max_lag
  #   tau ........... TRUE mean elapsed Δt per lag bin (ms)  <- honest x-axis
  #   tau_lin ....... lag × t_c + SW (linearized x-axis = Vogler's exact T_LIN)
  #   msd ........... ensemble MSD per lag (nm²)
  #   npairs ........ pairs contributing to each lag
  #   n_trk_lag ..... tracks contributing to each lag
  #   keep_lag ...... logical, TRUE where enough tracks support the lag
  #   track_msd ..... lag × track matrix of per-track TAMSDs (nm²)
  #   tk_names ...... TIDs, in the column order of track_msd
  #   msd_q25/q75 ... inter-quartile band across tracks
  #   t_cycle, has_z, ndim, track_weighted, n_dropped_negD, n_tracks
  #
  # NOTE: units are nm² and ms throughout. Consumers convert to µm²/s at the
  # point of display (divide nm²/ms by 1e3).
  # ============================================================================
  # ── Effective MSD divisor ───────────────────────────────────────────────────
  # MSD = 2 · d_eff · D · τ. d_eff is NOT simply the number of recorded coordinates:
  # it is set by how many directions the molecule can actually move in, and by how
  # that motion is foreshortened when projected onto the image plane.
  #
  #   motion   tracked        d_eff   why
  #   volume   3D              3      3 DOF on 3 axes - the plain case
  #   volume   2D              2      dropping z leaves 2 independent axes
  #   surface  3D              2      2 DOF however many coordinates are recorded
  #   surface  2D              1 + ⟨n_z²⟩, geometry-dependent (below)
  #
  # The last row is the subtle one: a molecule on a curved surface steps in its local
  # tangent plane, and the fraction of that step that survives projection depends on
  # how the surface is tilted. Averaging over the shape gives ⟨n_z²⟩ = 1/3 for a
  # sphere, 1/2 for a cylinder whose axis lies in the image plane, and 1 for a flat
  # membrane parallel to it - so d_eff = 4/3, 3/2 and 2. Verified by simulation to
  # three decimals. A membrane tilted by θ gives 1 + cos²θ, between the rod and flat
  # values; it is not offered separately.
  #
  # Getting this wrong is a first-order error on D, not a refinement: using the
  # recorded dimensionality for a surface-bound molecule under-reports D by 1/3 when
  # tracked in 3D, and by 1/3 (sphere) or 1/4 (rod) when tracked in 2D.
  # Extended for 1-DOF motion (a motor on a filament). Recorded in 3D the answer is
  # just the motion's own dimensionality; recorded in 2D it is reduced by whatever
  # fraction of a step points out of the image plane, which is ⟨n_z²⟩ (surface) or
  # ⟨u_z²⟩ (filament) for the shape in use:
  #   surface  on a sphere    ⟨n_z²⟩ = 1/3 → d_eff = 4/3
  #            on a rod       ⟨n_z²⟩ = 1/2 → d_eff = 3/2   (axis in the image plane)
  #            on a membrane  ⟨n_z²⟩ = 1   → d_eff = 2     (tilt θ gives 1 + cos²θ)
  #   filament in-plane                    → d_eff = 1
  # An unspecified shape ("none") is taken to lie in the image plane, the common
  # flat-mount case and the least surprising default - it applies no reduction.
  # The simulated free-motion nulls understand closed/flat surfaces only. "None" is
  # genuinely unbounded, and a filament has no surface to diffuse on, so both fall back
  # to the unbounded case - the null then reduces to the infinite-plane one, which is
  # the right answer rather than a fudge.
  .geom_for_null <- function(g) if (g %in% c("sphere", "rod", "ellipsoid", "membrane")) g else "free"

  # Human-readable name of the shared geometry, for panels that only echo it.
  .geom_label <- function(g) switch(g %||% "none",
    none = "None - unbounded", sphere = "Sphere", rod = "Rod / spherocylinder",
    ellipsoid = "Ellipsoid", membrane = "Flat membrane", filament = "Filament / linear track", g)

  .msd_d_eff <- function(d_obs) {
    motion <- input$trk_motion %||% "volume"
    if (identical(motion, "volume")) return(as.numeric(d_obs))
    d_mot <- if (identical(motion, "filament")) 1 else 2
    if (d_obs >= 3L) return(d_mot)               # all coordinates recorded: no projection loss
    if (identical(motion, "filament")) return(1) # a track lying in the image plane
    switch(input$trk_geom_type %||% "none",
           sphere = 4 / 3, ellipsoid = 4 / 3, rod = 3 / 2, membrane = 2, 2)
  }

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
    has_z_auto <- .has_real_z(df)
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

    withProgress(message = "Computing MSD…", value = 0, {
      for (k in seq_along(sp)) {
        incProgress(1 / n_tracks)
        tb <- sp[[k]]; n <- nrow(tb)
        if (n < min_pts) next
        # Cumulative cycle index (Bela's TOT) for the integer lag - robust to
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
    # Median of the same per-track TAMSDs, so the plotted band has a central line
    # that is commensurable with it (the eMSD is a mean of a skewed distribution).
    msd_med   <- apply(track_msd, 1L, stats::median, na.rm = TRUE)

    min_trk <- max(2L, floor(0.10 * max(n_trk_lag)))
    keep_lag <- n_trk_lag >= min_trk
    # Linearized (Vogler/MinSpt) time axis = k·t_c + SW, matching his T_LIN exactly
    # (simulate_cycles(k) = k·t_c + SW_OVERHEAD; get_cycle_time subtracts the SW).
    # The +SW is a CONSTANT offset, so it does NOT change the fitted slope / D - it
    # only shifts the curve right by SW and aligns the OLSF intercept (σ_loc) with his.
    list(lag = seq_len(max_lag), tau = tau_real, tau_lin = seq_len(max_lag) * t_c + cc_sw,
         msd = msd_pair,
         npairs = np_lag, n_trk_lag = n_trk_lag, keep_lag = keep_lag,
         track_msd = track_msd, tk_names = tk_names,
         msd_q25 = msd_q25, msd_q75 = msd_q75, msd_med = msd_med,
         t_cycle = t_c, has_z = has_z, ndim = if (has_z) 3L else 2L,
         # MSD = 2*d_eff*D*tau. d_eff follows the MOTION and its projection, not the
         # recorded coordinate count -- see .msd_d_eff.
         d_eff = .msd_d_eff(if (has_z) 3L else 2L),
         coef_eff = 2 * .msd_d_eff(if (has_z) 3L else 2L),
         track_weighted = track_weight, n_dropped_negD = n_dropped_negD,
         n_tracks = sum(colSums(tk_np) > 0))
  })

  # ── Free-diffusion-on-a-geometry MSD reference ──────────────────────────────
  # On a closed cell the MSD of a FREE molecule bends and eventually saturates, so a
  # bend in the measured curve is not by itself evidence of trapping. This simulates
  # free diffusion on the chosen shape at the fitted D and returns the curve it would
  # produce, on the same lag grid and the same observed coordinates as the data.
  #
  # The geometry is taken from the Displacement vs Diffusion panel rather than
  # duplicated here - one cell definition, used by both tests.
  #
  # D fed in is the FITTED D from this panel, which is already corrected for motion
  # degrees of freedom and projection (see .msd_d_eff), i.e. a true surface D. That is
  # why no projection re-calibration is applied, unlike the net-displacement null which
  # receives a raw per-axis D.
  msd_geom_null <- reactive({
    req(isTRUE(input$trk_msd_geom))
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); req(!is.null(ub))
    fr <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
    req(!is.null(fr), is.finite(fr$D), fr$D > 0)
    geom <- .geom_for_null(input$trk_geom_type %||% "none")
    params <- list(
      R    = max(1, as.numeric(input$trk_disp_R    %||% 500)),
      r    = max(1, as.numeric(input$trk_disp_r    %||% 400)),
      Lcyl = max(1, as.numeric(input$trk_disp_L    %||% 2000)),
      a    = max(1, as.numeric(input$trk_disp_a    %||% 500)),
      b    = max(1, as.numeric(input$trk_disp_b    %||% 500)),
      c    = max(1, as.numeric(input$trk_disp_c    %||% 900)),
      mtilt = max(0, as.numeric(input$trk_disp_mtilt %||% 0)))
    # Anything that is not volume-filling sits on the boundary (same mapping the
    # displacement null uses), so filament and surface both simulate on the surface.
    mode <- if (identical(input$trk_motion %||% "volume", "volume")) "volume" else "surface"
    keep <- ub$keep_lag & is.finite(ub$msd)
    lags <- ub$lag[keep]
    if (!length(lags)) return(NULL)
    # Lag grid: log-spaced across the whole range so the overlay reads well on any
    # scale, PLUS an even spread inside the fit window. The fit-window points are what
    # the D-shortfall estimate below is regressed on, and a purely log grid leaves that
    # regression resting on a handful of long lags, which made it noisy.
    lo_w <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
    hi_w <- max(lo_w + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    hi_w <- min(hi_w, max(lags))
    grid <- unique(c(
      round(exp(seq(log(1), log(max(lags)), length.out = 40))),
      round(seq(lo_w, hi_w, length.out = min(30L, max(3L, hi_w - lo_w + 1L))))))
    grid <- sort(grid[grid >= 1L])
    # Cost scales with the LARGEST lag, not the number of points, so keep the ensemble
    # large - the long-lag points are the noisy ones and they carry the signal.
    M <- 4000L
    res <- withProgress(message = "Simulating free diffusion on the geometry…", value = 0.4,
      sim_free_geom_msd(grid, D_surf = fr$D, geom = geom, params = params,
                        cyc_ms = ub$t_cycle, proj2d = (ub$ndim == 2L),
                        mode = mode, M = M))
    if (is.null(res) || !nrow(res)) return(NULL)
    res$tau_lin <- res$lag * ub$t_cycle
    attr(res, "geom") <- geom; attr(res, "mode") <- mode
    attr(res, "dims") <- .disp_geom_dims(geom, params)
    attr(res, "D") <- fr$D
    attr(res, "coef") <- ub$coef_eff %||% (2 * ub$ndim)
    res
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
      # Bela's exact iterative loop - converges on the lowest, least-noisy lags.
      msd_fit_olsf_iterative(tau_use[ok], ub$msd[ok], ndim = ub$ndim, R = trk_blur()$R_msd,
                             dt_cycle = ub$t_cycle, N = sum(ok),
                             coef_override = ub$coef_eff)
    } else {
      # Robust full-window-anchored Michalet (default).
      msd_fit_olsf(tau_use[ok], ub$msd[ok], ub$npairs[ok],
                   k = sum(ok), ndim = ub$ndim,
                   auto_pmin = TRUE, R = trk_blur()$R_msd, Npts = 1e4, dt_cycle = ub$t_cycle,
                   coef_override = ub$coef_eff)
    }
  }

  # .msd_model_fit() (linear / confined / directed / anomalous) lives in est_msd.R.

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
                   isTRUE(input$trk_uni_iterative), ub$t_cycle,
                   coef_override = ub$coef_eff, R = trk_blur()$R_msd)
  }

  # Cached per-track analysis from the UNIFIED MSD - the single source for the
  # time-averaged D (tDIFF), the per-track D histogram, the D/CI clustering and the
  # per-track CSV export. Runs the same per-track OLSF as the ensemble (iterative
  # when toggled, else anchored) on every track's TAMSD over the fit window, and
  # computes a confinement index (CI) from the same windowed curve. Cached so the
  # heavy loop only re-runs when the bundle or the window/method changes.
  # Returns list(tDIFF µm²/s, n_tD, per_D, pertrack = data.frame).
  #
  # WHY TWO DIFFUSION NUMBERS - eDIFF vs tDIFF (they answer different questions)
  #   eDIFF ("ensemble")   fit ONE line to the pooled ensemble MSD curve.
  #                        -> the population's typical D. Computed by
  #                           msd_unified_fit() on msd_unified_bundle()$msd.
  #   tDIFF ("time-average") fit EACH track separately, then average the D values.
  #                        -> the average individual molecule. Computed here.
  #
  #   For a homogeneous, ergodic sample the two agree. When they DISAGREE that is
  #   itself the result - it means the population is heterogeneous or
  #   non-ergodic (a few fast molecules pull the ensemble fit, while the typical
  #   molecule is slower, or vice versa). The ratio is reported as
  #
  #        ergodicity = <D>_time-average / D_ensemble      (1.0 = ergodic)
  #
  #   Only tracks with D > 0 enter tDIFF (a fitted negative slope is noise, not a
  #   molecule diffusing backwards).
  #
  #   CONFINEMENT INDEX (CI): the ratio of the MSD curve's late-window rise to
  #   its early-window rise. Free diffusion rises linearly so CI ~ 1; a confined
  #   molecule's MSD flattens towards a plateau, so CI < 1.
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
        msd_fit_olsf_iterative(tt, yy, ndim = ub$ndim, R = trk_blur()$R_msd, dt_cycle = t_c_e, N = N_tr,
                               coef_override = ub$coef_eff)
      else
        msd_fit_olsf(tt, yy, rep(1L, length(tt)), k = length(tt), ndim = ub$ndim,
                     auto_pmin = TRUE, R = trk_blur()$R_msd, Npts = N_tr, dt_cycle = t_c_e,
                     coef_override = ub$coef_eff)
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
                                               R = trk_blur()$R_msd, dt_cycle = ub$t_cycle, N = sum(winc),
                                               coef_override = ub$coef_eff)
          else msd_fit_olsf(tau_sel[winc], emsd[winc], rep(1L, sum(winc)), k = sum(winc),
                            ndim = ub$ndim, auto_pmin = TRUE, R = trk_blur()$R_msd, Npts = 1e4,
                            dt_cycle = ub$t_cycle, coef_override = ub$coef_eff),
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
        # mode = "lines" with a transparent line, NOT mode = "none": passing a line
        # object with mode="none" makes plotly print "A line object has been specified,
        # but lines is not in the mode / Adding lines to the mode..." to the console.
        type = "scatter", mode = "lines", fill = "toself",
        fillcolor = "rgba(70,130,200,0.18)",
        line = list(color = "rgba(0,0,0,0)", width = 0),
        # Named for what it is: the spread of per-track TAMSDs, NOT a confidence
        # band on the plotted ensemble value. The line is a MEAN of that same
        # distribution, and per-track TAMSD is strongly right-skewed, so the mean
        # sits high in the band - near Q3 - with Q1 far below it.
        hoverinfo = "skip", name = "IQR across tracks", showlegend = TRUE)
    }
    fig <- fig %>% add_trace(x = tau, y = msd, type = "scatter", mode = "lines+markers",
        marker = list(size = 5, color = "steelblue"), line = list(color = "steelblue", width = 2),
        # Name the curve for what it actually is: eMSD = mean of the per-track TAMSDs
        # (track-weighted), otherwise the pair-weighted ensemble average.
        name = if (isTRUE(ub$track_weighted)) "eMSD (track-weighted)"
               else "Ensemble MSD (pair-weighted)")
    # Median of the same per-track TAMSDs. The eMSD above is a MEAN of a strongly
    # right-skewed distribution, so it sits high inside the IQR band - near Q3, with
    # Q1 far below. Without the median next to it the line looks like it is sitting
    # outside its own band. Legend-toggleable, like every other trace.
    if (show_band && !is.null(ub$msd_med)) {
      mm <- ub$msd_med[ok] / 1e6
      okm <- is.finite(mm)
      if (any(okm)) fig <- fig %>% add_trace(
        x = tau[okm], y = mm[okm], type = "scatter", mode = "lines",
        line = list(color = "rgba(70,130,200,0.95)", width = 1.6, dash = "dot"),
        name = "median track")
    }
    # Free-diffusion-on-the-cell reference. Drawn in the same green as the geometry
    # band in Displacement vs Diffusion, since it is the same idea applied to the MSD.
    gm <- tryCatch(msd_geom_null(), error = function(e) NULL)
    if (!is.null(gm) && nrow(gm)) {
      gx <- if (linear_tau) gm$tau_lin else gm$tau
      fig <- fig %>% add_trace(
        x = gx, y = gm$msd / 1e6, type = "scatter", mode = "lines",
        line = list(color = "#1B9E5A", width = 2, dash = "dash"),
        name = sprintf("free on %s (%s)", attr(gm, "geom"), attr(gm, "mode")),
        hovertemplate = "free-on-geometry<br>τ %{x:.4g} ms<br>MSD %{y:.4g} µm²<extra></extra>")
    }
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
                    else sprintf("Confined: D=%.3g µm²/s, L_c=%.0f nm", fr$D / 1e3, fr$L_conf),
        directed  = sprintf("Directed: D=%.3g µm²/s, v=%.3g µm/s", fr$D / 1e3, fr$v_um_s),
        anomalous = sprintf("Anomalous: D=%.3g, α=%.2f", fr$D / 1e3, fr$alpha),
        sprintf("OLSF D=%.3g µm²/s (%d pts)", fr$D / 1e3, fr$n_points %||% NA))
      # Draw the fit SOLID only across the lags it was fitted to, and faint beyond.
      # OLSF picks an optimal (usually small) number of low lags, but the fitted
      # curve was being drawn across the whole axis - on a bending MSD that puts the
      # line far above the data at long lags and reads as a failed fit rather than an
      # extrapolation. Splitting the trace makes the fitted range explicit.
      # OLSF fits the first n_points of the WINDOWED curve, i.e. lags starting at the
      # lower cutoff (regime_x in tau), so the fitted range runs from there. The
      # non-linear models use the whole window, so they fall back to the upper limit.
      fit_hi <- {
        v_lag  <- suppressWarnings(as.numeric(input$trk_uni_max_lag))
        tau_hi <- if (is.finite(v_lag) && v_lag > 0) v_lag * ub$t_cycle else max(fx, na.rm = TRUE)
        n_used <- suppressWarnings(as.integer(fr$n_points %||% NA))
        # The fitted range ends at the n_used-th point OF THE WINDOWED CURVE, so read
        # it off the tau grid itself. Adding n_used * t_cycle to the cutoff is only
        # right on the linearized axis; on the mean-Delta-t axis the bins are spaced
        # wider than one cycle (they carry dark cycles), so that arithmetic lands
        # short and the solid segment renders far shorter than the range it stands for.
        tw <- sort(tau[is.finite(tau)])
        if (is.finite(n_used) && n_used > 0 && length(tw))
          min(tau_hi, tw[min(n_used, length(tw))])
        else tau_hi
      }
      in_fit <- is.finite(fx) & fx <= fit_hi
      # One curve, two segments: solid-dashed over the lags actually fitted, lighter
      # dotted beyond. Extrapolating an OLSF fit past its optimal point count is what
      # the method warns against, so the split stays - but both halves share a legend
      # group, so they read (and toggle) as a single trace.
      if (any(in_fit)) fig <- fig %>% add_trace(
        x = fx[in_fit], y = fy[in_fit], type = "scatter", mode = "lines",
        legendgroup = "msdfit",
        line = list(color = "#E03020", width = 2.5, dash = "dash"), name = fit_name)
      if (any(!in_fit)) fig <- fig %>% add_trace(
        x = fx[!in_fit], y = fy[!in_fit], type = "scatter", mode = "lines",
        legendgroup = "msdfit", showlegend = FALSE,
        line = list(color = "rgba(224,48,32,0.55)", width = 2, dash = "dot"),
        name = "fit, extrapolated",
        hovertemplate = "extrapolated beyond fitted lags<br>tau %{x:.4g} ms<br>MSD %{y:.4g} um2<extra></extra>")
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
                        if (is.finite(eD)) sprintf("%.3g", eD) else "-",
                        if (is.finite(tD)) sprintf("%.3g", tD) else "-")
              else sprintf("Cluster %d eMSD", clk)
        fig <- fig %>% add_trace(x = tau[okc], y = emsd_clk[okc], type = "scatter", mode = "lines",
          line = list(color = cl$palette[clk], width = 2.8), name = nm)
      }
    }
    regime <- regime_x
    # N(τ) on the right axis - same warm style as the standard plot.
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
                # Legend BELOW the plot: the right-hand side is occupied by the N(τ)
                # axis, and plotly's default right-hand legend sits on top of it.
                # y is in paper coordinates (a fraction of the plotting area), so it
                # is derived from the element height rather than hard-coded: the gap
                # below the axes then stays put when TRK_MSD_PLOT_PX changes.
                legend = list(font = list(color = fg, size = 10),
                              orientation = "h", x = 0, xanchor = "left",
                              y = -TRK_MSD_LEGEND_GAP /
                                   max(60, TRK_MSD_PLOT_PX - TRK_MSD_MARGIN_T - TRK_MSD_MARGIN_B),
                              yanchor = "top"),
                margin = list(l = 55, r = 60, t = TRK_MSD_MARGIN_T, b = TRK_MSD_MARGIN_B),
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

  # ── Cycle-binning diagnostic ────────────────────────────────────────────────
  # Re-inspects whether displacements are being binned correctly by cycle. Reads
  # the existing bundle (no recompute). Two signals per cycle-lag:
  #
  #   • BLUE, pairs per lag (left, log) - how much data backs each lag. Jagged at low
  #     lag because dark gaps make some cycle-lags reachable more often than others;
  #     then falls away, because a track of N localizations offers only ~(N−k) pairs at
  #     lag k and short tracks stop contributing at all. Where the count collapses the
  #     MSD goes noisy/saw-toothed - hence the "counts drop past N×" marker.
  #
  #   • RED, real elapsed Δt / linearized axis (right) - a rise to a PLATEAU is the
  #     EXPECTED, benign shape, not a fault. Per step, real elapsed = dcyc·t_c + SW,
  #     so a pair spanning k cycles in m steps takes k·t_c + m·SW, while the linearized
  #     axis (Vogler's T_LIN) charges SW just once: k·t_c + SW. Hence
  #
  #         ratio(k) = (k·t_c + m·SW) / (k·t_c + SW)
  #
  #     which is exactly 1 at lag 1 (m = 1) and climbs to a ceiling
  #
  #         plateau = 1 + (m/k)·SW/t_c ,   m/k → fraction of cycles that start a step
  #
  #     i.e. ≤ 1 + SW/t_c, hit only if every step is integrated. Verified in silico
  #     (t_c = 0.5, SW = 0.08, 35% dark: predicted 1.1041, observed 1.1033). So the
  #     PLATEAU HEIGHT reads out the overhead, and its excess over 1 is the absolute-D
  #     difference between the two axis conventions (it cancels in A-vs-B comparisons).
  #     What would actually be wrong: a ratio that keeps climbing without flattening,
  #     or one that saw-tooths - that is a mis-estimated t_c, not overhead.
  trk_msd_bindiag <- reactive({
    ub <- msd_unified_bundle(); req(!is.null(ub))
    lag   <- ub$lag; tc <- ub$t_cycle
    npair <- ub$npairs; ntrk <- ub$n_trk_lag
    keep  <- ub$keep_lag & is.finite(npair) & npair > 0
    # Compare real elapsed against the ACTUAL linearized axis (ub$tau_lin = k·t_c + SW,
    # Vogler's T_LIN), not the bare k·t_c - that's the axis the D is fit on. The two
    # coincide at lag 1 (both carry one SW) and diverge as real time accumulates an
    # extra SW per step; that divergence is the real-vs-linearized D difference.
    ratio <- ifelse(npair > 0 & is.finite(ub$tau) & is.finite(ub$tau_lin) & ub$tau_lin > 0,
                    ub$tau / ub$tau_lin, NA_real_)
    # Trust boundary: last kept lag whose pair count is still ≥ 20% of the peak.
    pk <- if (any(keep)) max(npair[keep], na.rm = TRUE) else NA_real_
    trust_lag <- if (is.finite(pk)) {
      good <- which(keep & npair >= 0.20 * pk)
      if (length(good)) max(lag[good]) else NA_integer_
    } else NA_integer_
    # Observed ceiling of the red curve: median ratio over the top decile of kept lags.
    # Its excess over 1 is (steps per cycle)·SW/t_c - the overhead read-out.
    kl  <- lag[keep & is.finite(ratio)]
    plateau <- if (length(kl) >= 10L)
      stats::median(ratio[keep & is.finite(ratio) &
                          lag >= stats::quantile(kl, 0.90, names = FALSE)], na.rm = TRUE)
    else NA_real_
   
    # UI says so. Use the TSA / cycle-time panel to actually determine t_c.
    sw_eff  <- if (length(ub$tau_lin) && is.finite(ub$tau_lin[1]) && is.finite(tc))
                 ub$tau_lin[1] - tc else NA_real_
    ceiling_ratio <- if (is.finite(sw_eff) && is.finite(tc) && tc > 0)
                       1 + sw_eff / tc else NA_real_
    tc_flag <- if (!is.finite(plateau) || !is.finite(ceiling_ratio)) NA_character_
               else if (plateau < 0.995) "too_large"
               else if (plateau > ceiling_ratio * 1.005) "too_small"
               else "ok"
    list(lag = lag, npair = npair, ntrk = ntrk, ratio = ratio, keep = keep,
         tc = tc, trust_lag = trust_lag, plateau = plateau,
         sw_eff = sw_eff, ceiling_ratio = ceiling_ratio, tc_flag = tc_flag,
         min_lag = max(1L, as.integer(input$trk_uni_min_lag %||% 4L)),
         max_lag = max(3L, as.integer(input$trk_uni_max_lag %||% 50L)))
  })

  output$trk_msd_bindiag_plot <- renderPlotly({
    req(isTRUE(input$trk_msd_unified), isTRUE(input$trk_msd_show_bindiag))
    d <- tryCatch(trk_msd_bindiag(), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(d), "Click 'Plot MSD' first."))
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    k  <- d$keep
    fig <- plot_ly() %>%
      add_trace(x = d$lag[k], y = d$npair[k], type = "scatter", mode = "lines+markers",
                line = list(color = "rgba(70,110,180,0.9)", width = 2),
                marker = list(size = 3, color = "rgba(70,110,180,0.9)"),
                name = "pairs / lag", yaxis = "y") %>%
      add_trace(x = d$lag[k], y = d$ratio[k], type = "scatter", mode = "lines",
                line = list(color = "#C0392B", width = 2),
                name = "real Δt / linearized axis", yaxis = "y2")
    # Reference: real Δt == the linearized axis (ratio = 1) on the right axis.
    shp <- list(list(type = "line", x0 = min(d$lag[k]), x1 = max(d$lag[k]),
                     y0 = 1, y1 = 1, yref = "y2",
                     line = list(color = if (dark) "rgba(200,200,200,0.4)" else "rgba(90,90,90,0.4)",
                                 width = 1, dash = "dot")))
    # Shade the fit window, and mark the trust boundary.
    shp[[length(shp)+1L]] <- list(type = "rect", x0 = d$min_lag, x1 = d$max_lag,
      y0 = 0, y1 = 1, yref = "paper", fillcolor = "rgba(120,160,120,0.10)",
      line = list(width = 0), layer = "below")
    ann <- list()
    if (is.finite(d$trust_lag)) {
      shp[[length(shp)+1L]] <- list(type = "line", x0 = d$trust_lag, x1 = d$trust_lag,
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "#E67E22", width = 1.5, dash = "dash"))
      ann[[1]] <- list(x = d$trust_lag, y = 1, yref = "paper", yanchor = "bottom",
        text = sprintf("counts drop past %d×", d$trust_lag), showarrow = FALSE,
        font = list(size = 9, color = "#E67E22"))
    }
    fig %>% layout(
      xaxis = list(title = "cycle lag", color = fg),
      yaxis = list(title = "pairs / lag", type = "log", color = fg,
                   gridcolor = if (dark) "rgba(255,255,255,0.12)" else "rgba(0,0,0,0.10)"),
      yaxis2 = list(title = "real Δt / linearized axis", overlaying = "y", side = "right",
                    color = "#C0392B", zeroline = FALSE),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(orientation = "v", x = 1.08, xanchor = "left", y = 1,
                    font = list(color = fg)),
      margin = list(r = 10),
      shapes = shp, annotations = ann,
      uirevision = "msd_bindiag"
    )
  })

  output$trk_msd_bindiag_note <- renderUI({
    req(isTRUE(input$trk_msd_unified), isTRUE(input$trk_msd_show_bindiag))
    d <- tryCatch(trk_msd_bindiag(), error = function(e) NULL); req(!is.null(d))
    k <- d$keep & d$lag >= d$min_lag & d$lag <= d$max_lag
    mean_ratio <- if (any(k)) mean(d$ratio[k], na.rm = TRUE) else NA_real_
    bias_pct   <- if (is.finite(mean_ratio)) 100 * (mean_ratio - 1) else NA_real_
    linear_tau <- isTRUE(input$trk_uni_linear_tau)
    tags$div(style = "font-size:0.82em; color:#555; margin:4px 0;",
      if (is.finite(mean_ratio))
        tags$div(
          tags$b("Linearized vs real-Δt axis: "),
          sprintf("real Δt averages %+.1f%% vs the linearized axis (lag × cycle + SW) over the fit window. ", bias_pct),
          if (is.finite(bias_pct) && abs(bias_pct) > 1)
            tags$span(
              "This shifts absolute D by about that much and cancels in comparisons between datasets. ",
              if (linear_tau)
                "Linearized axis in use, the published convention (Vogler et al. 2025)."
              else
                "Real-\u0394t axis in use, as a cross-check; the linearized axis is the published convention.")
          else tags$span(style = "color:#3a7d3a;", "The two axes agree here.")),
      if (is.finite(d$trust_lag))
        tags$div(style = "margin-top:2px;",
          tags$b("Pair-count trust boundary: "),
          sprintf("counts stay \u2265 20%% of peak up to lag %d (\u2248 %.3g ms); keep the upper fit limit at or below it.",
                  d$trust_lag, d$trust_lag * d$tc)),
      tags$div(style = "color:#888; margin-top:2px;",
        tags$b("Red curve: "),
        "real time over lag \u00d7 cycle. It should rise from 1 to a plateau between 1 and 1 + SW/t_c. ",
        if (is.finite(d$plateau))
          tags$span("Plateau here: ", tags$b(sprintf("%.3f", d$plateau)),
                    sprintf(" (+%.1f%%)", 100 * (d$plateau - 1)),
                    if (is.finite(d$ceiling_ratio))
                      sprintf(", band 1.000\u2013%.3f (SW \u2248 %.3g ms). ",
                              d$ceiling_ratio, d$sw_eff) else ". ") else NULL,
        tags$b("Cycle-time check: "),
        if (identical(d$tc_flag, "too_large"))
          tags$span(style = "color:#c03030;",
            sprintf("\u26a0 %.3f is below 1, so the cycle time is too large. Re-check it in TSA.",
                    d$plateau))
        else if (identical(d$tc_flag, "too_small"))
          tags$span(style = "color:#c03030;",
            sprintf("\u26a0 %.3f is above the %.3f ceiling, so the cycle time is too small. Re-check it in TSA.",
                    d$plateau, d$ceiling_ratio))
        else if (identical(d$tc_flag, "ok"))
          tags$span(style = "color:#3a7d3a;", "Within the band, consistent with the cycle time in use.")
        else NULL,
        " A mildly over-estimated cycle time can still land inside the band.",
        tags$br(),
        tags$b("Blue curve: "),
        "pairs per lag (log axis). Keep the upper fit limit where it is well populated.")
    )
  })

  # Cache the unified MSD's σ_loc (the LINEAR OLSF intercept - never a non-linear
  # model's offset, which is not separable) as the FALLBACK noise floor, used only when
  # neither a manual value nor the CVE estimate is available.
  #
  # This intercept is an extrapolation of the fit window back to τ = 0, so it absorbs
  # anything that is not linear-in-τ at short lags, not just localization noise. On
  # confined or otherwise bending data it therefore reads HIGH - simulated confined
  # tracks with a true σ_loc of 2.4 nm returned 3.7 nm from the linear intercept. The
  # CVE estimate is model-free (it reads the anti-correlation of consecutive steps,
  # which is a pure noise signature) and is preferred wherever available; that ordering
  # is enforced in trk_sigma_loc_eff.
  #
  # Require strictly > 0: a fit whose intercept collapses to exactly zero has not
  # measured anything, and caching 0 would silently disable every noise correction.
  observeEvent(msd_unified_bundle(), {
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); if (is.null(ub)) return()
    fr <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
    if (!is.null(fr) && is.finite(fr$sigma_loc) && fr$sigma_loc > 0)
      trk_sigma_loc_rv(fr$sigma_loc)
  }, ignoreInit = TRUE)

  output$trk_sigma_loc_status <- renderUI({
    eff <- trk_sigma_loc_eff()
    txt <- if (is.finite(eff))
             htmltools::tagList("Effective ", sigloc(),
               sprintf(" = %.2f nm (%s). Used as the noise floor everywhere.",
                       eff, trk_sigma_loc_src()))
           else htmltools::tagList(sigloc(),
                       " not set - estimate it in Analysis → Data Quality & Calibration ",
                       "(no MSD needed), or enter it manually. Noise-floor corrections stay ",
                       "inactive until then.")
    htmltools::div(style = "font-size:0.8em; color:#666; margin-top:-4px;", txt)
  })

  # eDIFF 95% CI by CLUSTER-BOOTSTRAP over tracks: resample track TAMSD columns with
  # replacement, rebuild the (track-weighted) ensemble MSD, refit - the honest CI
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
    # per-track OLSF D, D>0) - computed in the cached trk_unified_tdiff reactive so
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
    # Say which MSD-to-D relation is in force whenever it is not the plain one, since
    # it scales D directly.
    dof_note <- {
      de <- ub$d_eff %||% ub$ndim
      if (abs(de - ub$ndim) > 1e-9)
        htmltools::div(style = "color:#555; font-size:0.86em; margin-top:2px;",
          htmltools::HTML(sprintf(
            "Restricted motion: MSD = %.4g·D·τ (d<sub>eff</sub> = %.4g) on %dD coordinates%s. Reading it as %dD motion would report D a factor %.3g too low.",
            ub$coef_eff, de, ub$ndim,
            if (ub$ndim == 2L) sprintf(", projected from a %s", .geom_label(input$trk_geom_type)) else "",
            ub$ndim, de / ub$ndim)))
      else NULL
    }
    extra <- switch(model,
      confined  = if (isTRUE(fr$degenerate))
                    "  |  \u26a0 no plateau in the window, so L_c is unconstrained. Raise the upper limit, or the population is not confined on this timescale."
                  else sprintf("  |  confinement size L_c = %.0f nm%s  ·  plateau = %.3g µm²  ·  τ_c = %.3g ms",
                               fr$L_conf,
                               if (is.finite(sig_eff) && is.finite(fr$L_conf) && fr$L_conf < 2 * sig_eff)
                                 # extra[] flows through sprintf into an escaped string, so no
                                 # HTML here - use plain 2σ notation instead of an underscore.
                                 sprintf(" ⚠ below precision (2σ = %.0f nm - unresolved)", 2 * sig_eff) else "",
                               fr$plateau_um2, fr$tau_c_ms),
      directed  = {
        vnoise <- if (is.finite(sig_eff) && is.finite(ub$t_cycle) && ub$t_cycle > 0)
                    sqrt(2 * ub$ndim) * sig_eff / ub$t_cycle else NA_real_   # nm/ms = µm/s
        sprintf("  |  velocity v = %.3g µm/s%s%s", fr$v_um_s,
                if (is.finite(vnoise) && is.finite(fr$v_um_s) && fr$v_um_s < vnoise)
                  sprintf(" ⚠ at/below noise floor (%.3g µm/s - consistent with noise)", vnoise) else "",
                if (isTRUE(fr$degenerate)) " ⚠ fit did not converge - parameters not identifiable over this window" else "")
      },
      anomalous = sprintf("  |  α = %.2f (%s)%s", fr$alpha,
                          if (fr$alpha < 0.9) "sub-diffusive / confined"
                          else if (fr$alpha > 1.1) "super-diffusive / directed" else "≈ Brownian",
                          if (isTRUE(fr$degenerate)) " ⚠ fit did not converge - parameters not identifiable over this window" else ""),
      "")
    method_lbl <- switch(model,
      linear = if (identical(fr$method, "OLSF-iterative")) "iterative OLSF" else "anchored Michalet",
      confined = "confined model", directed = "directed model", anomalous = "anomalous model", "model")
    D_lbl <- if (identical(model, "linear")) "eDIFF (ensemble)" else "D (local, ensemble)"
    # How much is the CELL bending the curve inside the fit window? Compares the
    # simulated free-on-geometry MSD to a straight line at the same D over the same
    # lags: a large shortfall means the fitted D is being pulled down by the geometry,
    # not by biology, and the window should be shortened (or the bend interpreted as
    # the cell).
    geom_note <- {
      gm <- tryCatch(msd_geom_null(), error = function(e) NULL)
      if (is.null(gm) || !nrow(gm)) NULL else {
        lo <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
        hi <- max(lo + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
        w  <- gm$lag >= lo & gm$lag <= hi
        # Fit the SELECTED model to the simulated free-on-geometry curve, over the same
        # window, and report what it returns. Nothing is confined, drifting or anomalous
        # in that simulation, so whatever comes back is the shape alone - the baseline
        # the measured parameter has to beat. This generalises the D-shortfall check to
        # every model: on a closed shape the confined model reports a spurious corral
        # and the anomalous model a sub-diffusive exponent, purely from the ceiling.
        fg <- if (sum(w) < 4L) NULL else
          tryCatch(.msd_model_fit(gm$tau[w], gm$msd[w], rep(1000L, sum(w)), ub$ndim,
                                  model, FALSE, ub$t_cycle, coef_override = attr(gm, "coef"), R = trk_blur()$R_msd),
                   error = function(e) NULL)
        if (is.null(fg)) NULL else {
          Dtrue <- attr(gm, "D")
          rat   <- if (is.finite(fg$D) && Dtrue > 0) fg$D / Dtrue else NA_real_
          # Model-specific line: what does pure geometry hand back?
          say <- switch(model,
            confined = if (is.finite(fg$L_conf))
                sprintf("would itself fit as a %s of L_c &asymp; %.0f nm%s",
                        if (isTRUE(fg$degenerate)) "non-converged corral" else "corral",
                        fg$L_conf,
                        if (is.finite(fr$L_conf))
                          sprintf(", against your fitted %.0f nm", fr$L_conf) else "")
              else NULL,
            directed = if (is.finite(fg$v_um_s))
                sprintf("would itself fit as v &asymp; %.3g &micro;m/s%s", fg$v_um_s,
                        if (is.finite(fr$v_um_s)) sprintf(", against your fitted %.3g", fr$v_um_s) else "")
              else NULL,
            anomalous = if (is.finite(fg$alpha))
                sprintf("would itself fit as &alpha; &asymp; %.2f%s", fg$alpha,
                        if (is.finite(fr$alpha)) sprintf(", against your fitted %.2f", fr$alpha) else "")
              else NULL,
            if (is.finite(rat)) sprintf("fitted the same way returns about %.0f%% of its true D",
                                        5 * round(100 * rat / 5)) else NULL)
          if (is.null(say)) NULL else {
            # Warn when the shape alone already explains most of what was measured.
            bad <- switch(model,
              confined  = is.finite(fg$L_conf) && is.finite(fr$L_conf) && fg$L_conf < 2 * fr$L_conf,
              anomalous = is.finite(fg$alpha) && fg$alpha < 0.95,
              directed  = FALSE,
              is.finite(rat) && rat < 0.95)
            col <- if (isTRUE(bad)) "#b26a00" else "#3a7d3a"
            htmltools::div(style = sprintf("color:%s; font-size:0.86em; margin-top:2px;", col),
              htmltools::HTML(sprintf(
                "Simulated free diffusion on this %s (%s, %s) %s, %s",
                attr(gm, "geom"), attr(gm, "dims"), attr(gm, "mode"), say,
                if (isTRUE(bad))
                  "geometry, not biology: shorten the upper fit limit."
                else "so the shape does not explain your value over this window.")))
          }
        }
      }
    }
    htmltools::div(style = "font-size:0.87em; color:#444; margin:6px 0;",
      dof_note, geom_note,
      htmltools::tags$b("MSD: "),
      sprintf("%s = %.4g%s  |  tDIFF (per-track mean, n=%d) = %.4g µm²/s%s  |  ",
              D_lbl, fr$D / 1e3, eD_ci, n_tD, tDIFF, tD_ci),
      sigloc(),
      # The offset is only separable in the LINEAR model. In the non-linear models it
      # competes with the plateau / exponent for the same part of the curve and is
      # bounded below at zero, so it routinely collapses to exactly 0 - which is a
      # failure to measure it, not a measurement of zero. Say so rather than print 0.00.
      if (identical(model, "linear") || (is.finite(fr$sigma_loc) && fr$sigma_loc > 0.01))
        sprintf(" = %.2f nm", fr$sigma_loc)
      else " not separable in this model",
      sprintf("  |  %s  |  ergodicity = %.2f%s", method_lbl, ergo, extra),
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
      # Confined/Anomalous models - flagged as beyond-standard-for-MINFLUX.
      # One line, not a lecture: the reasoning is in the panel's Method box.
      if (identical(model, "linear"))
        htmltools::div(style = "color:#8a6d00; font-size:0.82em; margin-top:4px;",
          htmltools::HTML("&#9432; "),
          htmltools::tags$b("Free-Brownian fit"),
          " - under confinement this D is apparent/local. Cross-check the confinement index."),
      # Per-cluster eDIFF / tDIFF breakdown (when clustering is enabled).
      local({
        cdiff <- tryCatch(trk_unified_cluster_diff(), error = function(e) NULL)
        if (is.null(cdiff) || !nrow(cdiff)) return(NULL)
        rows <- lapply(seq_len(nrow(cdiff)), function(i) {
          r <- cdiff[i, ]
          htmltools::tags$div(style = sprintf("color:%s; font-size:0.85em; margin-left:8px;", r$color),
            sprintf("● Cluster %d: eDIFF = %s, tDIFF = %s µm²/s  (n = %d)", r$cluster,
                    if (is.finite(r$eDIFF)) sprintf("%.3g", r$eDIFF) else "-",
                    if (is.finite(r$tDIFF)) sprintf("%.3g", r$tDIFF) else "-", r$n))
        })
        htmltools::tagList(htmltools::br(),
          htmltools::tags$b(style = "font-size:0.85em;", "Per-cluster diffusion:"), rows)
      }))
  })


  # ==========================================================================
  # cycle_class_bundle() - INTEGRATED vs DARK CYCLE CLASSIFICATION (TSA core)
  # ==========================================================================
  #
  # GOAL
  #   Reconstruct WHAT THE MICROSCOPE WAS DOING between consecutive
  #   localizations. MINFLUX runs on a fixed acquisition cycle, but it does not
  #   report a position every cycle - if the molecule was dark (blinking) or too
  #   few photons arrived, that cycle yields nothing and the next localization
  #   appears one or more cycles later. This turns the raw timestamps into that
  #   cycle-level story, which is the foundation every Δt-aware estimator uses.
  #
  #   Adapted from Bela Vogler's MINFLUX Time-Signature Analysis (MFX-TSA),
  #   Vogler et al. 2025 (Commun. Biol.); reference implementation
  #   Toolbox_for_MINFLUX_SPT (cf. its `print_cycle_stats`).
  #
  # THE MATH (deliberately simple, and that is the point)
  #   For each consecutive interval within a track:
  #
  #        n = round(Δt / t_cycle)          cycles spanned by the step
  #
  #        n = 1  -> INTEGRATED: localized on the very next cycle, no gap.
  #        n >= 2 -> (n − 1) DARK cycles passed before the molecule reappeared.
  #
  #   Everything downstream is built on this one integer:
  #     - integrated fraction = the headline data-quality score
  #     - the gap-split threshold is expressed in dark cycles
  #     - the MSD's cumulative cycle index is cumsum(n)
  #     - CVE, JDD and directionality all filter on n == 1
  #
  # DIFFERENCE FROM VOGLER'S IMPLEMENTATION (stated honestly)
  #   His TSA SIMULATES the hardware timing (dwell time, pattern vertices, EOD
  #   and software overheads) and fits those parameters to the data. We infer
  #   t_cycle from the Δt distribution instead - no hardware simulation. That is
  #   why the app's cycle time is a data-driven estimate (cluster cycle, or an
  #   optional regression) rather than a fitted overhead model, and why the
  #   sequence-JSON import exists as a cross-check.
  #
  #   IMPORTANT: the MODAL Δt is a poor estimate of t_cycle when dark cycles are
  #   common, because the mode can land on a 2-cycle gap. The cluster-cycle
  #   estimator (smallest robust Δt cluster) is used instead.
  #
  # OBJECTS IT CALLS (inputs)
  #   dt_vectors() (per-track Δt), trk_cycle_time(), trk_cluster_cycle()
  #   inputs: trk_tsa_run (Run trigger), trk_cycle_dark_max
  #
  # RETURNS the per-interval cycle counts plus the summary statistics the TSA
  # panel displays: integrated fraction, dark-cycle distribution, per-track
  # cycle progression, and the before/after comparison for gap splitting.
  #
  # CAVEAT: descriptive, not inferential. If t_cycle is wrong, EVERY cycle count
  # is wrong - check the data-quality verdict and the sequence JSON before
  # trusting anything downstream.
  # ==========================================================================
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
          "Excellent tracking quality - almost every cycle yielded a localization."
        else if (cc$frac_integ > 0.7)
          "Good tracking - most cycles localized, modest dark fraction."
        else
          "Frequent dropouts. Splitting tracks at gaps (Display tab) keeps long dark gaps out of the analyses."
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
    fmt  <- function(v, pct = FALSE) if (is.null(v)) "-" else
      if (pct) sprintf("%.1f%%", v * 100) else format(round(v, 2), big.mark = ",")

    htmltools::div(style = "font-size:0.86em; color:#444; margin:8px 0;",
      htmltools::tags$b("Cycle stats - raw vs analyzed (after cleaning)"),
      if (!any_clean) htmltools::tags$span(style = "color:#888;",
        " (no cleaning active - enable Split / Trim / Quality gate in the Display tab to see the effect)"),
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
        "the Integrated fraction, cut total Dark cycles, and shrink the max gap - the analyzed ",
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
  
  # ============================================================================
  # jdd_bundle() - JUMP-DISTANCE DISTRIBUTION (how many diffusing populations?)
  # ============================================================================
  #
  # GOAL
  #   The MSD gives ONE average D for the whole dataset. But a sample often holds
  #   several populations at once (e.g. a bound slow fraction + a free fast one).
  #   This asks a different question: what does the DISTRIBUTION of single-step
  #   jump lengths look like, and how many populations does it take to explain it?
  #
  # THE MATH
  #   For 2D free diffusion each coordinate step is Gaussian,
  #        Δx, Δy ~ N(0, s²)  with  s² = 2·D·Δt
  #   so the step LENGTH r = sqrt(Δx² + Δy²) follows a RAYLEIGH distribution:
  #
  #        p(r) = (r / s²) · exp( −r² / (2 s²) ),        s² = 2·D·Δt
  #
  #   Its mean square is <r²> = 2s² = 4DΔt = 2dDΔt (d = 2), consistent with the
  #   MSD relation in the header. Inverting for D:
  #
  #        D = s² / (2·Δt)        and after sqrt(Δt) normalisation:  D = s'² / d
  #
  #   N co-existing populations -> a WEIGHTED MIXTURE of N chi-distributed components,
  #        p(r) = Σ_k w_k · C_d · r^(d−1) / s_k^d · exp( −r² / (2 s_k²) ),   Σ w_k = 1
  #   fitted here by maximum likelihood. Each fitted s_k gives that population's
  #   D_k, and w_k its abundance. d = spatial dimensions: d=2 is the classic Rayleigh
  #   (Schütz 1997), d=3 the Maxwell form used when z is included. s_k is the PER-AXIS
  #   spread either way, so D_k = s_k²/2 regardless of d.
  #
  # THE MINFLUX ADAPTATION (two routes, pick one)
  #   Jumps must be pooled at a COMMON Δt or the mixture is meaningless.
  #   (a) sqrt(Δt) normalisation (default): divide every jump by sqrt of its own
  #       Δt, which removes the Δt dependence for free diffusion (see header §2b).
  #   (b) "Integrated steps only": keep single-cycle steps, so every jump already
  #       shares one true Δt and NO normalisation assumption is needed. Cleaner -
  #       prefer it when the motion may not be free.
  #
  # OBJECTS IT CALLS (inputs)
  #   trk(), trk_cycle_time(), trk_sigma_loc_eff()
  #   inputs: trk_jdd_run (Run trigger), trk_jdd_min_pts, trk_jdd_n_comp,
  #           trk_jdd_norm_dt, trk_jdd_integrated_only, trk_jdd_min_dt/max_dt,
  #           trk_jdd_sigma_loc, trk_jdd_step_mode (+ bin_width / subsamp_n)
  #
  # STEP BY STEP
  #   1. Split into tracks; optionally coarsen each track first
  #      (consecutive steps | bin-then-median | every N-th localization).
  #   2. Per track, form step lengths r and their Δt.
  #   3. Keep steps inside [min_dt, max_dt]; if "integrated only", additionally
  #      require round(Δt/t_c) == 1.
  #   4. If normalising, r <- r / sqrt(Δt).
  #   5. Histogram the pooled jumps (upper tail truncated at the 99.5th pct so a
  #      few outliers do not stretch the axis).
  #   6. Fit the N-component chi-d mixture by MLE (Rayleigh in 2D, Maxwell in 3D).
  #      Parameters are optimised as log(s_k) and softmax weights, so widths stay
  #      positive and weights sum to 1 without constrained optimisation.
  #   7. Convert widths to D, subtracting the localization-noise contribution
  #      (`loc_bias`), and clamp at 0 - a component at the floor is noise, not motion.
  #   8. Compute r_noise = sqrt(2d)·σ_loc, the RMS length of a pure-noise step:
  #      anything left of it on the plot is noise-dominated.
  #
  # RETURNS a list:
  #   mids, dens ...... histogram (x, density) for plotting
  #   sigs[k], ws[k] .. fitted per-axis widths and mixture weights
  #   D_vals[k] ....... per-component D (noise-corrected, clamped >= 0)
  #   n_comp, norm_dt, integrated_only, t_cycle, dt0, use_z, d_fit (2 or 3)
  #   n_jumps, sigma_loc, r_noise, min_dt_used, max_dt_used, step_mode
  #
  # CAVEAT: the mixture assumes every population is FREELY diffusing. Confined or
  # anomalous motion does not follow this law, and extra components may then be
  # fitting curve shape rather than real populations. Keep N small.
  # ============================================================================
  jdd_bundle <- eventReactive(input$trk_jdd_run, {
    req(isTRUE(input$trk_show_jdd))
    df <- trk(); req(nrow(df) > 0)
    min_pts   <- max(2L, as.integer(input$trk_jdd_min_pts))
    norm_dt   <- isTRUE(input$trk_jdd_norm_dt)
    t_c_jdd   <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    # Integrated-only: keep ONLY single-cycle steps (round(Δt/t_c)==1). All jumps are
    # then at the same true Δt, so the fit needs no √Δt-pooling assumption at all - the
    # rigorous MINFLUX-native mode (at the cost of dropping dark-gap-adjacent steps).
    integ_only <- isTRUE(input$trk_jdd_integrated_only) && is.finite(t_c_jdd) && t_c_jdd > 0
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

    # Real axial data only. A placeholder all-zero z column would make this fit the 3D
    # Maxwell law to what are really 2D distances (dz ≡ 0), biasing D by −33%.
    # Motion blur needs the cycle time. Without it the fit would silently run with no
    # blur correction, which reads D about a third low at MINFLUX cycle times, so the
    # panel asks for the calibration instead of guessing.
    shiny::validate(shiny::need(is.finite(t_c_jdd) && t_c_jdd > 0,
      "Cycle time unknown - run Data Quality & Calibration first (the blur correction needs it)."))
    use_z_jdd <- .has_real_z(df)
    # Collection and fit live in est_jdd.R (pure, tested).
    withProgress(message = "Computing JDD…", value = 0, {
    cj <- jdd_collect_jumps(df, min_pts = min_pts, step_mode = step_mode, bin_width = bin_width,
                            subsamp_n = subsamp_n, min_dt = min_dt, max_dt = max_dt,
                            integ_only = integ_only, t_c = t_c_jdd, norm_dt = norm_dt,
                            use_z = use_z_jdd,
                            progress = function(f) incProgress(f, detail = "collecting steps"))
    shiny::validate(shiny::need(length(cj$jumps) >= 10,
      "Too few jumps - lower 'Min points per TID' or relax filters."))
    jf <- jdd_fit(cj$jumps, cj$dt_all, n_comp = n_comp, use_z = use_z_jdd, norm_dt = norm_dt,
                  sigma_loc = sigma_loc,
                  R_abs = trk_blur()$R_abs,
                  progress = function(f) incProgress(f, detail = "fitting jump-distance mixture"))

    list(mids = jf$mids, dens = jf$dens, sigs = jf$sigs, ws = jf$ws, D_vals = jf$D_vals,
         n_comp = n_comp, norm_dt = norm_dt, n_jumps = length(cj$jumps),
         use_z = use_z_jdd, d_fit = jf$d_fit, dt0 = jf$dt0, inv_dt = jf$inv_dt, sigma_loc = sigma_loc,
         r_noise = jf$r_noise, R_abs = trk_blur()$R_abs,
         step_mode = step_mode, integrated_only = integ_only, t_cycle = t_c_jdd,
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

    # Same chi-d PDF the MLE fitted (Rayleigh in 2D, Maxwell in 3D) - the drawn curve
    # must be the fitted curve, or the overlay silently disagrees with the table.
    .d_ov <- bd$d_fit %||% 2L
    .C_ov <- 1 / (2^(.d_ov / 2 - 1) * gamma(.d_ov / 2))
    rayleigh_comp <- function(r, sig, w)
      w * .C_ov * r^(.d_ov - 1) / sig^.d_ov * exp(-r^2 / (2 * sig^2))

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
    # \u03C3 is the per-axis spread in both 2D and 3D (chi-d PDF), so the /2 is dimension-
    # independent; only the fitted PDF shape changes (Rayleigh 2D / Maxwell 3D).
    d_formula <- sprintf("D = (\u03C3\u00B2 \u2212 2\u03C3\u1D38\u1D52\u1D9C\u00B2\u27E81/\u0394t\u27E9) / (2(1 \u2212 2R\u00B7t\u1D9C\u27E81/\u0394t\u27E9))   [%s]",
                         if (isTRUE(bd$use_z)) "3D, Maxwell" else "2D, Rayleigh")
    dt0_str   <- if (is.finite(bd$dt0)) sprintf("%.2g ms", bd$dt0) else "n/a"
    max_dt_str <- if (is.finite(bd$max_dt_used)) sprintf("%.4g ms", bd$max_dt_used) else "none"

    sigma_note <- if (show_D && bd$sigma_loc > 0 && is.finite(bd$dt0))
      tags$p(sprintf("\u03C3\u1D38\u1D52\u1D9C = %.1f nm  |  median \u0394t = %s  |  noise subtracted: %.2g nm\u00B2/ms  |  blur factor %.2f",
                     bd$sigma_loc, dt0_str, 2 * bd$sigma_loc^2 * (bd$inv_dt %||% (1 / bd$dt0)),
                     1 - 2 * (bd$R_abs %||% 0) * (bd$inv_dt %||% (1 / bd$dt0))),
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
            " (within loc. noise - step sizes explained by σₗₒₓ alone)"))
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
      # chi-d PDF as fitted (Rayleigh 2D / Maxwell 3D) - keep the export in step
      # with the plot and the fitted parameters.
      .d_ov <- bd$d_fit %||% 2L
      .C_ov <- 1 / (2^(.d_ov / 2 - 1) * gamma(.d_ov / 2))
      rayleigh_comp <- function(r, sig, w)
        w * .C_ov * r^(.d_ov - 1) / sig^.d_ov * exp(-r^2 / (2 * sig^2))
      out <- data.frame(mid_nm = bd$mids, density = bd$dens)
      out$fit_total <- Reduce(`+`, lapply(seq_len(bd$n_comp), function(k)
        rayleigh_comp(bd$mids, bd$sigs[k], bd$ws[k])
      ))
      for (k in seq_len(bd$n_comp))
        out[[paste0("comp_", k)]] <- rayleigh_comp(bd$mids, bd$sigs[k], bd$ws[k])
      write.csv(out, file, row.names = FALSE)
    }
  )




  # ============================================================================
  # conf_ci_bundle() - PER-TRACK CONFINEMENT FROM THE MSD SHAPE (CI + simulated null)
  # ============================================================================
  # The most sensitive per-track confinement test in the app. Uses the confinement
  # index already computed per track by the MSD:
  #
  #     CI = (MSD_late / MSD_early) / (tau_late / tau_early)
  #
  # WHY THIS ONE. It compares a track's own SHORT lags against its own LONG lags, so
  # it needs no external D and cannot be biased by D varying between molecules - the
  # confound that inflates the Displacement-vs-Diffusion split. It is not circular
  # (short and long lags are different information) and it inherits the cycle-index
  # binning, so irregular Δt is already handled.
  #
  # WHY IT NEEDS A SIMULATED NULL. For ideal free motion CI = 1, but the OBSERVED free
  # value is not: the 2·d·σ_loc² intercept sits on every MSD point and lifts the EARLY
  # window more than the late one, pushing free CI below 1. How far below depends on
  # σ_loc, D, cycle time, fit window and track length, so there is no universal
  # threshold - the null is simulated at the data's own conditions and matched per
  # track on localization count. Simulated free CI medians run ≈0.67–0.77 with the 5th
  # percentile moving 0.28 → 0.63 from N=50 to N=800, which is why matching on N matters.
  #
  # VALIDATION (in silico, MINFLUX-like: t_c 0.191 ms, σ_loc 10 nm, D 0.2 µm²/s, N=200):
  #   free → 7.3% confined / 89.2% free-like / 3.5% directed  (target 5/90/5)
  #   confined L=100 nm → 89.5% detected;  L=250 nm → 25.0% detected
  #   directed v=0.4 µm/s → NOT detected (6.2%) - see the caveat below.
  #
  # LIMITS. Detection is set by how far a free molecule would travel in the same time
  # relative to the confinement size: the resolution floor √(2·d·D·T) is reported so a
  # null result can be read as "no confinement below this scale" rather than "free".
  # CI is weak for DIRECTED motion (a modest drift barely lifts the late/early ratio) -
  # use Directionality or Speed & Transport for transport.
  conf_ci_bundle <- eventReactive(input$trk_conf_ci_run, {
    td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(td) && !is.null(td$pertrack) && nrow(td$pertrack) > 0,
      "Run the MSD first (Plot MSD) - the confinement index comes from the per-track MSD."))
    pt <- td$pertrack
    pt <- pt[is.finite(pt$conf_index), , drop = FALSE]
    shiny::validate(shiny::need(nrow(pt) >= 5L, "Too few tracks with a usable confinement index."))

    df  <- trk(); t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    shiny::validate(shiny::need(is.finite(t_c) && t_c > 0, "Cycle time unknown."))
    ub  <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
    ndim <- if (!is.null(ub)) ub$ndim else 2L
    sig <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    if (!is.finite(sig) || sig < 0) sig <- 0
    b   <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    D_use <- if (!is.null(b) && is.finite(b$D_cve) && b$D_cve > 0) b$D_cve else {
      dm <- stats::median(pt$D_um2_s[is.finite(pt$D_um2_s) & pt$D_um2_s > 0], na.rm = TRUE)
      if (is.finite(dm)) dm * 1e3 else NA_real_
    }
    shiny::validate(shiny::need(is.finite(D_use) && D_use > 0,
      "Need a D - run Data Quality & Calibration, or the MSD."))

    # Per-track localization count and the dataset's real dark-gap distribution, so the
    # null tracks are sampled exactly like the data.
    npts <- table(df$tid)
    pt$n_loc <- as.integer(npts[as.character(pt$tid)])
    pt <- pt[is.finite(pt$n_loc) & pt$n_loc >= 5L, , drop = FALSE]
    shiny::validate(shiny::need(nrow(pt) >= 5L, "Too few tracks after matching localization counts."))
    d2 <- df[order(df$tid, df$tim), ]
    same <- c(FALSE, d2$tid[-1] == d2$tid[-nrow(d2)])
    kk <- pmax(1L, as.integer(round(c(NA, diff(d2$tim))[same] / t_c)))
    kk <- kk[is.finite(kk) & kk >= 1L]
    gap_p <- if (length(kk) > 50L) { tb <- tabulate(pmin(kk, 15L), nbins = 15L); tb / sum(tb) }
             else c(0.55, 0.22, 0.09, 0.05, 0.04, 0.05)

    min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
    max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
    Ns <- pt$n_loc
    grid <- if (length(unique(Ns)) <= 8L) sort(unique(Ns))
            else sort(unique(pmax(5L, round(stats::quantile(Ns, seq(0, 1, length.out = 8L),
                                                            names = FALSE)))))
    nsim <- max(50L, min(500L, as.integer(input$trk_conf_ci_nsim %||% 200L)))

    # ── D grid ────────────────────────────────────────────────────────────────
    # The free CI depends on x = σ_loc²/(D·t_c), so the null is matched on each
    # track's OWN D as well as its length - otherwise the slow half of a
    # heterogeneous population is measured against a null built for the fast half and
    # comes out "confined". The grid collapses to a single point (the cheap path, and
    # identical to matching on N alone) when that dependence cannot bite: either the
    # noise is negligible at the fastest tracks, or the per-track D barely varies.
    Dtr <- pt$D_um2_s * 1e3                       # per-track D, nm²/ms
    ok_D <- is.finite(Dtr) & Dtr > 0
    Dtr[!ok_D] <- D_use
    x_of <- function(D) if (sig > 0) sig^2 / (D * t_c) else 0
    D_rng <- if (any(ok_D)) stats::quantile(Dtr[ok_D], c(0.05, 0.95), names = FALSE)
             else c(D_use, D_use)
    needs_D_grid <- sig > 0 && is.finite(D_rng[1]) && D_rng[1] > 0 &&
                    x_of(D_rng[2]) > 0.3 &&        # noise matters somewhere
                    (D_rng[2] / D_rng[1]) > 2      # and D actually varies
    D_grid <- if (needs_D_grid)
                exp(seq(log(D_rng[1]), log(D_rng[2]), length.out = 4L))
              else D_use
    ncell <- length(grid) * length(D_grid)
    # Noise for the null: when σ_loc comes from the photon model, resample the
    # per-localization values (the MSD offset is set by mean σ², which a long-tailed
    # photon count pushes above median σ²). A manual or CVE/MSD σ stays static.
    sig_null <- if (identical(trk_sigma_loc_src(), "photon budget") && !is.null(b) &&
                    sum(is.finite(b$sig_phot)) >= 50L) b$sig_phot[is.finite(b$sig_phot)]
                else sig
    # The simulation is ~85% of the wait, so it owns 0 -> 0.85 of the bar and the
    # per-track matching that follows owns the rest. sim_free_ci_null() reports a
    # CUMULATIVE fraction, so it is set, not added: incrementing would run the bar
    # to 100% at the end of the simulation and leave the matching looking hung.
    nul <- withProgress(message = "Confinement: simulating the free null…", value = 0,
      sim_free_ci_null(grid, D_grid = D_grid, n_sim = nsim, t_c = t_c, sigma = sig_null,
                       ndim = ndim, gap_p = gap_p, min_lag = min_lag, max_lag = max_lag,
                       progress = function(f) setProgress(0.85 * f,
                         detail = sprintf("%d of %d (length, D) cells", round(f * ncell), ncell)),
                       R_blur = trk_blur()$R_msd))   # same blur coefficient as the MSD fits
    # nearest cell in (log N, log D)
    idx <- withProgress(message = "Confinement: matching tracks to the null…", value = 0.85, {
      ln <- log(nul$n); lD <- log(nul$D)
      vapply(seq_len(nrow(pt)), function(i) {
        which.min((ln - log(Ns[i]))^2 + (lD - log(Dtr[i]))^2)
      }, integer(1))
    })
    pt$null_D  <- nul$D[idx]
    pt$ci_lo05 <- nul$lo05[idx]; pt$ci_med <- nul$med[idx]; pt$ci_hi95 <- nul$hi95[idx]
    pt$ci_class <- ifelse(pt$conf_index < pt$ci_lo05, "confined",
                   ifelse(pt$conf_index > pt$ci_hi95, "directed", "free-like"))
    # Resolution floor: how far a free molecule travels over the median track duration.
    durs <- tapply(df$tim, df$tid, function(v) diff(range(v)))
    T_med <- stats::median(as.numeric(durs), na.rm = TRUE)
    attr(pt, "resolv") <- sqrt(2 * ndim * D_use * T_med)
    attr(pt, "D_use") <- D_use; attr(pt, "sigma") <- sig
    attr(pt, "T_med") <- T_med; attr(pt, "null") <- nul; attr(pt, "nsim") <- nsim
    attr(pt, "D_matched") <- length(D_grid) > 1L
    attr(pt, "x_noise") <- x_of(D_use)
    pt
  })

  output$trk_conf_ci_summary <- renderUI({
    req(isTRUE(input$trk_show_conf))
    pt <- tryCatch(conf_ci_bundle(), error = function(e) NULL)
    if (is.null(pt) || !nrow(pt)) return(NULL)
    f_c <- mean(pt$ci_class == "confined"); f_f <- mean(pt$ci_class == "free-like")
    f_d <- mean(pt$ci_class == "directed")
    rv2 <- attr(pt, "resolv")
    tags$div(style = "font-size:0.9em; margin:6px 0; padding:6px 8px; border-left:3px solid #1a8a1a; background:rgba(26,138,26,0.05);",
      tags$b("MSD confinement index vs a simulated free null"),
      tags$span(style = "color:#666; font-weight:normal;", " - the primary call on confinement"), tags$br(),
      tags$b(style = "color:#1a8a1a; font-size:1.02em;",
        sprintf("Confined %.0f%%  ·  free-like %.0f%%  ·  directed %.0f%%",
                100 * f_c, 100 * f_f, 100 * f_d)),
      tags$div(style = "color:#666; margin-top:2px;",
        sprintf("n = %d tracks · null: σ_loc = %.2f nm, D = %.3g µm²/s, %d free tracks per grid point, matched on localization count%s",
                nrow(pt), attr(pt, "sigma"), attr(pt, "D_use") / 1e3, attr(pt, "nsim"),
                if (isTRUE(attr(pt, "D_matched"))) " and per-track D" else "")),
      tags$div(style = "color:#666; margin-top:2px;",
        htmltools::HTML(sprintf(
          "Resolution floor <b>%.0f nm</b> in the median track (%.1f ms) - confinement larger than this is invisible.",
          rv2, attr(pt, "T_med")))),
      tags$div(style = "color:#888; margin-top:2px;",
        "~5% per tail by construction; weak for directed motion - use Directionality or Speed & Transport."))
  })

  output$trk_dl_conf_ci_csv <- downloadHandler(
    filename = function() paste0("confinement_index_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) write.csv(conf_ci_bundle(), file, row.names = FALSE)
  )

  # CI against localization count, with the simulated free band behind it. Plotting
  # against N makes the N-matching visible: the band narrows as tracks get longer,
  # which is exactly why a single fixed threshold would be wrong.
  output$trk_conf_ci_plot <- renderPlotly({
    req(isTRUE(input$trk_show_conf))
    pt <- tryCatch(conf_ci_bundle(), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(pt) && nrow(pt) > 0,
      "Run the MSD, then click 'Run confinement index'."))
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    nul <- attr(pt, "null")
    # When the null is matched on D too there is a band per D value; show the one
    # nearest the ensemble D so the plot stays a single readable band.
    if (!is.null(nul) && length(unique(nul$D)) > 1L) {
      Dref <- unique(nul$D)[which.min(abs(log(unique(nul$D)) - log(attr(pt, "D_use"))))]
      nul <- nul[nul$D == Dref, , drop = FALSE]
    }
    cols <- c(confined = "#2c7fb8", `free-like` = "#9e9e9e", directed = "#e31a1c")
    fig <- plot_ly()
    if (!is.null(nul) && nrow(nul) > 1) {
      o <- order(nul$n)
      fig <- fig %>%
        add_trace(x = nul$n[o], y = nul$hi95[o], type = "scatter", mode = "lines",
                  line = list(color = "rgba(120,120,120,0)"), showlegend = FALSE,
                  hoverinfo = "skip") %>%
        add_trace(x = nul$n[o], y = nul$lo05[o], type = "scatter", mode = "lines",
                  fill = "tonexty", fillcolor = "rgba(150,150,150,0.18)",
                  line = list(color = "rgba(120,120,120,0)"),
                  name = if (isTRUE(attr(pt, "D_matched")))
                           "free 5–95% (at the ensemble D)" else "free 5–95% (simulated)",
                  hoverinfo = "skip") %>%
        add_trace(x = nul$n[o], y = nul$med[o], type = "scatter", mode = "lines",
                  line = list(color = "#888", dash = "dash", width = 2),
                  name = "free median", hoverinfo = "skip")
    }
    for (cl in names(cols)) {
      idx <- pt$ci_class == cl
      if (!any(idx)) next
      fig <- fig %>% add_trace(
        x = pt$n_loc[idx], y = pt$conf_index[idx], type = "scatter", mode = "markers",
        marker = list(size = 6, color = cols[[cl]], opacity = 0.75,
                      line = list(width = 0.3, color = "rgba(0,0,0,0.3)")),
        name = sprintf("%s (%d)", cl, sum(idx)),
        text = sprintf("TID %s<br>%d locs · CI %.3f<br>%s",
                       pt$tid[idx], pt$n_loc[idx], pt$conf_index[idx], cl),
        hovertemplate = "%{text}<extra></extra>")
    }
    fig %>% layout(
      xaxis = list(title = "localizations in track", color = fg, type = "log"),
      yaxis = list(title = "confinement index (CI)", color = fg),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(font = list(color = fg, size = 10), orientation = "h",
                    x = 0, xanchor = "left", y = -0.34, yanchor = "top"),
      margin = list(l = 55, r = 20, t = 10, b = 110),
      uirevision = "conf_ci")
  })




  # ============================================================================
  # trk_disp_bundle() - NET DISPLACEMENT vs DIFFUSION  (Step 1 of the geometry test)
  # ----------------------------------------------------------------------------
  # Per track: the NET start→end displacement R_net = |r(T) − r(0)|, the max excursion
  # from the start, the (sampling-dependent) path length, and the duration T. Each is
  # compared to the FREE-diffusion expectation at the ensemble D (from the CVE panel):
  # for 2D Brownian motion ⟨R_net²⟩ = 4·D·T, and R_net itself is Rayleigh-distributed
  # (scale √(2DT)), so a single track has a wide spread even when free.
  #
  # WHY 2D + WHY THIS IS ONLY STEP 1: displacement is measured in the 2D image plane
  # (x, y). This is the INFINITE-PLANE null - it does NOT account for a finite/curved
  # surface, which itself caps how far a FREE molecule can get. So "sub-diffusive" here
  # means "confined OR geometrically constrained"; Step 2 (simulate free diffusion ON a
  # chosen geometry at this D) separates the two. Path length is reported but NEVER used
  # to classify - like path speed it scales with sampling rate, not the molecule.
  #
  # NOISE: endpoints carry σ_loc, so ⟨R_net²⟩ gains 4σ² (2 endpoints × 2 axes). That
  # (and the extreme-value analogue for Rmax) is subtracted when σ_loc is set.
  # ============================================================================
  trk_disp_bundle <- eventReactive(input$trk_disp_run, {
    req(isTRUE(input$trk_show_disp))
    df <- trk(); req(nrow(df) > 0)
    min_pts <- max(3L, as.integer(input$trk_disp_min_pts %||% 10L))

    # Dimensionality: net displacement (and the free null) in 2D or 3D. "Auto" uses z when
    # it carries real signal. 2D data → 2D; 3D data → 3D (Rayleigh vs Maxwell null). The
    # geometry null (Step 2) then projects (2D) or stays 3D to match.
    z_has <- .has_real_z(df)
    dmode <- input$trk_disp_dims %||% "auto"
    use3d <- (identical(dmode, "auto") && z_has) || identical(dmode, "3d")
    if (use3d && !z_has) use3d <- FALSE                       # asked for 3D but no z signal
    d_obs <- if (use3d) 3L else 2L
    # ── Degrees of freedom of the MOTION vs of the OBSERVATION ─────────────────
    # d_obs is how many coordinates are recorded; d_mot is how many the molecule can
    # actually move in. They differ for a molecule confined to a SURFACE (a membrane,
    # a cell envelope) but tracked in 3D: 2 DOF observed on 3 axes.
    # Two consequences, both handled here:
    #   1. A per-axis estimator (CVE, or MSD with ndim = 3) spreads 2-DOF motion over
    #      3 axes, so it reports D_meas = D_true · d_mot/d_obs - exactly 2/3 for a
    #      surface in 3D (simulated: 0.675). Undo it: D_mot = D_meas · d_obs/d_mot.
    #   2. R_net² follows χ² with d_mot degrees of freedom, not d_obs.
    # Getting (1) right but not (2) leaves the MEAN correct while the BAND is the wrong
    # shape: on simulated free-on-surface data the old d_obs-only form gave
    # 11.4% sub / 81.4% free / 7.2% super instead of 5/90/5; with d_mot it gives
    # 5.2/90.0/4.8. When d_mot == d_obs (2D data, or true volume diffusion) both terms
    # cancel and this reduces exactly to the previous expression.
    # Degrees of freedom of the motion, from the shared Sample geometry setting.
    d_mot <- {
      mo <- input$trk_motion %||% "volume"
      if (identical(mo, "filament")) 1L
      else if (identical(mo, "surface")) min(2L, d_obs)
      else d_obs
    }

    keep <- is.finite(df$x) & is.finite(df$y) & is.finite(df$tim) &
            (if (use3d) is.finite(df$z) else TRUE)
    df <- df[keep, ]

    # ── Which D sets the free null? ────────────────────────────────────────────
    # The null is √(2·d·D·T), so D directly sets the verdict: a larger D raises the band
    # and flags MORE tracks sub-diffusive. It must be the SHORT-TIME (local) mobility.
    #   CVE      - short-lag, noise-corrected (Vestergaard), model-free. DEFAULT, because
    #              it is measured before any confinement bends the MSD.
    #   MSD fit  - fitted over the lag window. If the population IS confined the MSD has
    #              already bent down inside that window, so this D is DEFLATED and the
    #              test loses power (it under-calls confinement). Offered for comparison.
    #   Manual   - a literature/known value, or for sensitivity testing.
    dsrc  <- input$trk_disp_dsrc %||% "cve"
    b     <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    D_cve <- if (!is.null(b)) b$D_cve else NA_real_
    D_msd <- {
      ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
      fr <- if (!is.null(ub)) tryCatch(msd_unified_modelfit(ub), error = function(e) NULL) else NULL
      if (!is.null(fr) && is.finite(fr$D)) fr$D else NA_real_
    }
    D_man <- suppressWarnings(as.numeric(input$trk_disp_dman %||% NA_real_)) * 1e3  # µm²/s → nm²/ms
    D_use <- switch(dsrc, msd = D_msd, manual = D_man, D_cve)
    d_lab <- switch(dsrc, msd = "MSD fit", manual = "manual", "CVE")
    # In-plane (motion) D. CVE and the MSD fit are per-axis over d_obs recorded axes, so
    # for 2-DOF motion tracked in 3D they read low by d_mot/d_obs and are rescaled here.
    # A MANUAL D is taken to be the true motion D already, so it is passed through.
    D_mot <- if (identical(dsrc, "manual")) D_use else D_use * d_obs / d_mot
    shiny::validate(shiny::need(is.finite(D_use) && D_use > 0, switch(dsrc,
      msd    = "No MSD D yet - press Plot MSD in the MSD panel (or switch the D source to CVE).",
      manual = "Enter a positive D (µm²/s) for the manual D source.",
      "Need an ensemble D first - run Analysis → Data Quality & Calibration (Estimate σloc & quality).")))
    sig    <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    # Endpoint noise: 2 endpoints × d_obs axes × σ².
    noise2 <- if (is.finite(sig) && sig > 0) 2 * d_obs * sig^2 else 0

    sp <- split(df, df$tid)
    rows <- withProgress(message = "Displacement vs diffusion…", value = 0,
      lapply(seq_along(sp), function(i) {
        if (i %% 100L == 0L) incProgress(100 / length(sp))
        tb <- sp[[i]]; if (nrow(tb) < min_pts) return(NULL)
        tb <- tb[order(tb$tim), ]; ne <- nrow(tb)
        T_ms <- tb$tim[ne] - tb$tim[1]
        if (!is.finite(T_ms) || T_ms <= 0) return(NULL)
        dz2 <- if (use3d) (tb$z[ne] - tb$z[1])^2 else 0
        Rnet <- sqrt(max(0, (tb$x[ne] - tb$x[1])^2 + (tb$y[ne] - tb$y[1])^2 + dz2 - noise2))
        ex2  <- (tb$x - tb$x[1])^2 + (tb$y - tb$y[1])^2 + if (use3d) (tb$z - tb$z[1])^2 else 0
        Rmax <- sqrt(max(0, max(ex2) - noise2))
        pl   <- sum(sqrt(diff(tb$x)^2 + diff(tb$y)^2 + if (use3d) diff(tb$z)^2 else 0))
        # ⟨R_net²⟩ = 2·d_mot·D_mot·T. With d_mot = d_obs this is √(4DT) in 2D and
        # √(6DT) in 3D; for a surface tracked in 3D it becomes √(4·D_mot·T).
        free_rms <- sqrt(2 * d_mot * D_mot * T_ms)
        data.frame(tid = names(sp)[i], n_pts = ne, dur_ms = T_ms,
                   Rnet_nm = Rnet, Rmax_nm = Rmax, path_nm = pl,
                   free_rms_nm = free_rms,
                   rho = if (free_rms > 0) Rnet / free_rms else NA_real_,
                   stringsAsFactors = FALSE)
      }))
    rows <- rows[!vapply(rows, is.null, logical(1))]
    shiny::validate(shiny::need(length(rows) > 0, "Too few tracks - lower 'Min points per TID'."))
    res <- do.call(rbind, rows)

    # Per-track free (infinite-plane) 5/95 band and two-tailed classification. R_net² for
    # free diffusion is 2·D_mot·T·χ²(d_mot) - the χ² degrees of freedom follow the MOTION,
    # not the number of recorded coordinates (see the d_mot note above).
    sc2 <- 2 * D_mot * res$dur_ms
    res$free_lo05 <- sqrt(sc2 * stats::qchisq(0.05, d_mot))
    res$free_hi95 <- sqrt(sc2 * stats::qchisq(0.95, d_mot))
    res$class <- ifelse(res$Rnet_nm < res$free_lo05, "sub-diffusive",
                 ifelse(res$Rnet_nm > res$free_hi95, "super-diffusive", "free-consistent"))
    # Self-contained columns so the CSV carries everything Step 2 needs.
    res$D_ens_um2_s  <- D_use / 1e3
    res$D_motion_um2_s <- D_mot / 1e3
    res$sigma_loc_nm <- if (is.finite(sig)) sig else NA_real_
    res$dims         <- d_obs
    res$motion_dof   <- d_mot
    res$D_source     <- d_lab
    attr(res, "D_use") <- D_use; attr(res, "sigma") <- sig; attr(res, "d_obs") <- d_obs
    attr(res, "d_mot") <- d_mot; attr(res, "D_mot") <- D_mot
    attr(res, "d_lab") <- d_lab
    attr(res, "D_cve") <- D_cve; attr(res, "D_msd") <- D_msd
    res
  })

  # Echo the shared geometry inside the displacement panel, so the Step-2 band states
  # which shape it is using without owning a second copy of the control.
  output$trk_disp_geom_echo <- renderUI({
    g <- input$trk_geom_type %||% "none"
    params <- list(R = input$trk_disp_R, r = input$trk_disp_r, Lcyl = input$trk_disp_L,
                   a = input$trk_disp_a, b = input$trk_disp_b, c = input$trk_disp_c,
                   mtilt = input$trk_disp_mtilt)
    dims <- tryCatch(.disp_geom_dims(.geom_for_null(g), params), error = function(e) "")
    htmltools::div(style = "font-size:0.85em; margin:2px 0 4px;",
      htmltools::tags$b("Geometry: "), .geom_label(g),
      if (nzchar(dims)) htmltools::span(style = "color:#666;", sprintf("  (%s)", dims)) else NULL,
      if (identical(.geom_for_null(g), "free"))
        htmltools::div(style = "color:#b26a00;",
          "No bounding shape, so this band matches the unbounded null. Set a shape in Sample geometry to test a closed cell.")
      else NULL)
  })

  # ── Step 2: GEOMETRY-AWARE free null (free diffusion ON the chosen geometry) ───────
  # Replaces the infinite-plane band with a band simulated from free diffusion on the
  # selected shape (sphere / rod / ellipsoid / membrane / free - the Diffusion Simulator's
  # geometries), so "sub-diffusive" can be split into TRULY confined vs merely
  # geometrically constrained. The measured 2D D is a PROJECTED D (< the true surface D on a
  # curved shape), so sim_free_geom_rnet empirically calibrates that per-geometry factor and
  # drives the null at the matching surface D → both bands agree at short T and the
  # closed-surface band diverges only where finite size + curvature bite. Same physics as
  # the Diffusion Simulator.
  # Labels for the current geometry's size parameters (used in the summary).
  .disp_geom_dims <- function(geom, p) switch(geom,
    sphere    = sprintf("R = %.0f nm", p$R),
    rod       = sprintf("r = %.0f nm, L = %.0f nm", p$r, p$Lcyl),
    ellipsoid = sprintf("semi-axes %.0f / %.0f / %.0f nm", p$a, p$b, p$c),
    membrane  = sprintf("tilt = %.0f°", p$mtilt),
    free      = "unbounded", "")
  trk_disp_geom_null <- reactive({
    if (!isTRUE(input$trk_disp_geom)) return(NULL)
    res <- trk_disp_bundle(); req(nrow(res) > 0)
    D_use <- attr(res, "D_use")
    geom  <- .geom_for_null(input$trk_geom_type %||% "none")
    params <- list(
      R    = max(1, as.numeric(input$trk_disp_R    %||% 500)),
      r    = max(1, as.numeric(input$trk_disp_r    %||% 400)),
      Lcyl = max(1, as.numeric(input$trk_disp_L    %||% 2000)),
      a    = max(1, as.numeric(input$trk_disp_a    %||% 500)),
      b    = max(1, as.numeric(input$trk_disp_b    %||% 500)),
      c    = max(1, as.numeric(input$trk_disp_c    %||% 900)),
      mtilt = max(0, as.numeric(input$trk_disp_mtilt %||% 0)))
    cyc   <- tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_)
    if (!is.finite(cyc) || cyc <= 0) cyc <- stats::median(diff(sort(res$dur_ms))) # crude fallback
    if (!is.finite(cyc) || cyc <= 0) return(NULL)
    # Match the real DURATION distribution (R_net depends on T) over its FULL range.
    # This used to be truncated at the 98th percentile to bound the per-cycle loop, which
    # left the longest tracks - the ones carrying the most information about confinement -
    # judged against a flat-extrapolated band. Cost scales with the number of simulation
    # steps (longest duration / step), so instead of dropping those tracks the SIMULATION
    # STEP is coarsened: Brownian statistics are step-size invariant, so a larger step with
    # SD = √(2·D·dt) gives the same R_net distribution provided the step stays small
    # compared with the geometry. Measured on a rod (r = 400 nm, D = 200 nm²/ms): stepping
    # at 1×/5×/10×/25× the cycle gave median R_net 292/305/306/302 nm and p95 753/782/787/786
    # while the run time fell from 2.6 s to 0.09 s.
    # Simulate on a LOG-SPACED GRID spanning the observed durations, with replicates at
    # each point, rather than resampling the empirical durations. Resampling plus
    # equal-count binning put the last bin's MEDIAN far below the longest track, so the
    # band visibly stopped short of the data and the few longest tracks - the most
    # informative ones - were classified against a flat extrapolation. On a grid the
    # band is defined at the maximum duration by construction.
    dv <- res$dur_ms[is.finite(res$dur_ms) & res$dur_ms > 0]
    if (!length(dv)) return(NULL)
    d_rng <- range(dv)
    if (d_rng[1] <= 0 || !is.finite(d_rng[2])) return(NULL)
    n_grid_t <- 24L; n_rep <- 200L
    grid_t <- if (d_rng[2] / d_rng[1] > 1.05)
                exp(seq(log(d_rng[1]), log(d_rng[2]), length.out = n_grid_t))
              else rep(d_rng[2], 1L)
    durs  <- rep(grid_t, each = n_rep)
    T_max <- max(durs, na.rm = TRUE)
    step_budget <- 4000L
    dt_sim <- max(cyc, T_max / step_budget)
    # Keep the step well below the smallest geometry dimension, or the tangent-step-and-
    # reproject operator loses accuracy. Unbounded shapes have no such limit.
    gmin <- switch(geom, sphere = params$R, rod = params$r,
                   ellipsoid = min(params$a, params$b, params$c), Inf)
    if (is.finite(gmin) && gmin > 0 && is.finite(D_use) && D_use > 0) {
      dt_geom <- (0.05 * gmin)^2 / (2 * D_use)     # step SD ≤ 5% of the smallest dimension
      dt_sim  <- min(dt_sim, max(cyc, dt_geom))
    }
    # 2D data → simulate the 3D geometry and PROJECT (the cell is 3D even if tracking is 2D);
    # 3D data → keep all three axes.
    proj2d <- (attr(res, "d_obs") %||% 2L) == 2L
    # Same "Molecule moves" choice that sets the Step-1 degrees of freedom.
    # "surface" and "filament" both mean the molecule is on the boundary; only a
    # volume-filling molecule reflects around inside it.
    gmode  <- if (identical(input$trk_motion %||% "volume", "volume")) "volume" else "surface"
    rnet <- withProgress(message = "Simulating free-on-geometry null…", value = 0.4,
                         sim_free_geom_rnet(durs, D_use, geom, params, dt_sim, proj2d = proj2d,
                                            mode = gmode))
    if (!length(rnet)) return(NULL)
    D_surf <- attr(rnet, "D_surf") %||% D_use
    # One 5/50/95 summary per grid duration - the band therefore spans the full range of
    # observed track durations, ending exactly at the longest track.
    rs  <- as.numeric(rnet)
    gid <- rep(seq_along(grid_t), each = n_rep)[seq_along(rs)]
    bt  <- grid_t[sort(unique(gid))]
    blo <- as.numeric(tapply(rs, gid, stats::quantile, 0.05))
    bmd <- as.numeric(tapply(rs, gid, stats::quantile, 0.50))
    bhi <- as.numeric(tapply(rs, gid, stats::quantile, 0.95))
    # Classify each REAL track against the geometry null at its own duration.
    lo_r <- stats::approx(bt, blo, xout = res$dur_ms, rule = 2)$y
    hi_r <- stats::approx(bt, bhi, xout = res$dur_ms, rule = 2)$y
    gclass <- ifelse(res$Rnet_nm < lo_r, "confined",
              ifelse(res$Rnet_nm > hi_r, "directed", "free-on-geometry"))
    list(bt = bt, blo = blo, bmd = bmd, bhi = bhi, gclass = gclass, D_surf = D_surf,
         geom = geom, dims = .disp_geom_dims(geom, params), cyc = cyc, gmode = gmode)
  })

  # Blue = sub-diffusive (confined/geometric), grey = free-consistent, red = super-
  # diffusive (directed) - matching the Speed colouring (blue slow, red fast).
  .disp_class_col <- function(cls) c(`sub-diffusive` = "#2c7fb8",
                                     `free-consistent` = "#9e9e9e",
                                     `super-diffusive` = "#e31a1c",
                                     `confined` = "#2c7fb8",
                                     `free-on-geometry` = "#9e9e9e",
                                     `directed` = "#e31a1c")[cls]

  output$trk_disp_plot <- renderPlotly({
    req(isTRUE(input$trk_show_disp))
    res <- trk_disp_bundle(); req(nrow(res) > 0)
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    D_use <- attr(res, "D_use"); d_obs <- attr(res, "d_obs") %||% 2L
    # Band uses the MOTION dof / motion D (see trk_disp_bundle); identical to the old
    # expression whenever d_mot == d_obs.
    d_mot <- attr(res, "d_mot") %||% d_obs
    D_mot <- attr(res, "D_mot") %||% D_use
    tg  <- seq(max(1e-6, min(res$dur_ms)), max(res$dur_ms), length.out = 120)
    sc2 <- 2 * D_mot * tg   # R_net² = 2·D_mot·T·χ²(d_mot)
    med <- sqrt(sc2 * stats::qchisq(0.5, d_mot))
    lo  <- sqrt(sc2 * stats::qchisq(0.05, d_mot))
    hi  <- sqrt(sc2 * stats::qchisq(0.95, d_mot))
    gn  <- trk_disp_geom_null()   # NULL unless the Step-2 geometry option is on
    pt_class <- if (!is.null(gn)) gn$gclass else res$class
    fig <- plot_ly() %>%
      add_trace(x = tg, y = hi, type = "scatter", mode = "lines",
                line = list(color = "rgba(120,120,120,0)"), showlegend = FALSE, hoverinfo = "skip") %>%
      add_trace(x = tg, y = lo, type = "scatter", mode = "lines", fill = "tonexty",
                fillcolor = "rgba(150,150,150,0.15)", line = list(color = "rgba(120,120,120,0)"),
                name = "free 5–95% (unbounded)", hoverinfo = "skip") %>%
      add_trace(x = tg, y = med, type = "scatter", mode = "lines",
                line = list(color = "#888", dash = "dash", width = 2),
                name = "free median (unbounded)", hoverinfo = "skip")
    if (!is.null(gn)) {
      # Free-on-E.coli band (green): the geometry-aware null.
      fig <- fig %>%
        add_trace(x = gn$bt, y = gn$bhi, type = "scatter", mode = "lines",
                  line = list(color = "rgba(27,158,90,0)"), showlegend = FALSE, hoverinfo = "skip") %>%
        add_trace(x = gn$bt, y = gn$blo, type = "scatter", mode = "lines", fill = "tonexty",
                  fillcolor = "rgba(27,158,90,0.18)", line = list(color = "rgba(27,158,90,0)"),
                  name = "free 5–95% (on geometry)", hoverinfo = "skip") %>%
        add_trace(x = gn$bt, y = gn$bmd, type = "scatter", mode = "lines",
                  line = list(color = "#1B9E5A", width = 2),
                  name = "free median (on geometry)", hoverinfo = "skip")
    }
    fig <- fig %>%
      add_trace(x = res$dur_ms, y = res$Rnet_nm, type = "scatter", mode = "markers",
                marker = list(size = 6, color = .disp_class_col(pt_class),
                              line = list(width = 0.3, color = "rgba(0,0,0,0.3)")),
                text = sprintf("TID %s<br>dur %.1f ms · %d pts<br>R_net %.0f nm · ρ %.2f<br>%s",
                               res$tid, res$dur_ms, res$n_pts, res$Rnet_nm, res$rho, pt_class),
                hovertemplate = "%{text}<extra></extra>", name = "tracks")
    fig %>%
      layout(
        xaxis = list(title = "track duration T (ms)", color = fg, zeroline = FALSE),
        yaxis = list(title = "net displacement R_net (nm)", color = fg, zeroline = FALSE),
        paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
        legend = list(font = list(color = fg, size = 10), orientation = "h", y = -0.2),
        title = list(text = "Net displacement vs free-diffusion expectation",
                     font = list(color = fg, size = 12)),
        uirevision = "disp")
  })

  output$trk_disp_summary <- renderUI({
    req(isTRUE(input$trk_show_disp))
    res <- trk_disp_bundle(); req(nrow(res) > 0)
    D_use <- attr(res, "D_use"); sig <- attr(res, "sigma")
    n <- nrow(res)
    f_sub <- mean(res$class == "sub-diffusive"); f_free <- mean(res$class == "free-consistent")
    f_sup <- mean(res$class == "super-diffusive")
    T_med <- stats::median(res$dur_ms)
    d_obs <- attr(res, "d_obs") %||% 2L
    resolv <- sqrt(2 * d_obs * D_use * T_med)   # RMS free extent; confinement below it is invisible
    htmltools::div(style = "font-size:0.9em; line-height:1.55;",
      htmltools::HTML(sprintf(
        "<b>%d</b> tracks · <b>%dD</b> net displacement · D = <b>%.3g µm²/s</b> (%s)%s · median duration <b>%.1f ms</b>.",
        n, attr(res, "d_obs") %||% 2L, D_use / 1e3, attr(res, "d_lab") %||% "CVE",
        if (is.finite(sig)) sprintf(" · σ<sub>loc</sub> = %.2f nm", sig) else "", T_med)),
      {
        # Only worth stating when the motion dof differs from the recorded coordinates.
        dm <- attr(res, "d_mot") %||% d_obs; Dm <- attr(res, "D_mot") %||% D_use
        if (dm != d_obs)
          htmltools::div(style = "color:#555; margin-top:3px;",
            htmltools::HTML(sprintf(
              "Surface motion: %d DOF on %d axes &rarr; in-plane D = <b>%.3g µm²/s</b>, band uses χ²(%d).",
              dm, d_obs, Dm / 1e3, dm)))
        else NULL
      },
      {
        # Cross-check: if the two estimators disagree, the verdict is D-sensitive.
        Dc <- attr(res, "D_cve"); Dm <- attr(res, "D_msd")
        if (is.finite(Dc) && is.finite(Dm) && Dc > 0 && Dm > 0 && abs(log(Dm / Dc)) > log(1.3))
          htmltools::div(style = "color:#B8860B; margin-top:3px; font-size:0.93em;",
            htmltools::HTML(sprintf(
              "&#9888; CVE D = %.3g vs MSD D = %.3g \u00b5m\u00b2/s (%.1f\u00d7 apart), so the split below depends on which D is used.",
              Dc / 1e3, Dm / 1e3, max(Dm, Dc) / min(Dm, Dc))))
        else NULL
      },
      htmltools::br(),
      htmltools::HTML(sprintf(
        "Infinite-plane null: <b style='color:#2c7fb8;'>%.0f%% sub-diffusive</b> · <b style='color:#9e9e9e;'>%.0f%% free</b> · <b style='color:#e31a1c;'>%.0f%% super-diffusive</b> (~5%% per tail by construction)",
        100 * f_sub, 100 * f_free, 100 * f_sup)),
      {
        # This panel corroborates the confinement index; it does not compete with it.
        # The two ask different questions of the same track - the CI reads the shape of
        # the MSD curve across lags, this reads how far the track actually travelled -
        # so their percentages are not meant to match. What is informative is how far
        # they agree track by track, which is only shown once both have been run.
        pt <- tryCatch(conf_ci_bundle(), error = function(e) NULL)
        if (is.null(pt) || !nrow(pt) || !all(c("tid", "ci_class") %in% names(pt)))
          htmltools::div(style = "color:#555; margin-top:4px;",
            htmltools::HTML("Run the <b>confinement index</b>, the primary call, to see how far the two agree."))
        else {
          m <- merge(res[, c("tid", "class")], pt[, c("tid", "ci_class")], by = "tid")
          if (!nrow(m)) NULL else {
            ci_conf <- m$ci_class == "confined"
            both    <- sum(ci_conf & m$class == "sub-diffusive")
            htmltools::div(style = "color:#555; margin-top:4px;",
              htmltools::HTML(sprintf(
                "Of the <b>%d</b> tracks the confinement index calls confined, <b>%d</b> (%.0f%%) also fall below this band (%d tracks scored by both).",
                sum(ci_conf), both,
                if (sum(ci_conf)) 100 * both / sum(ci_conf) else 0, nrow(m))))
          }
        }
      },
      htmltools::div(style = "color:#555; margin-top:5px;",
        htmltools::HTML(sprintf(
          "Resolution floor: <b>%.0f nm</b> RMS free extent. Confinement on a larger scale is invisible.", resolv))),
      {
        # ONE ensemble D sets the band for every track, so a population with a wide
        # spread of per-track D reads as sub-diffusive at the slow end even when every
        # molecule is free. Quantified on simulated all-free populations: the
        # sub-diffusive call runs ~5.5% at no spread, ~6.9% at sd(ln D)=0.6 and ~12%
        # at sd(ln D)=1.0. Warn using the per-track D spread when the MSD has been run.
        td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
        Dv <- if (!is.null(td) && !is.null(td$pertrack)) td$pertrack$D_um2_s else NULL
        Dv <- Dv[is.finite(Dv) & Dv > 0]
        sdl <- if (length(Dv) >= 20L) stats::sd(log10(Dv)) else NA_real_
        if (is.finite(sdl) && sdl > 0.35)
          htmltools::div(style = "color:#B8860B; margin-top:5px; font-size:0.95em;",
            htmltools::HTML(sprintf(
              "&#9888; Wide per-track D (sd log<sub>10</sub>D = %.2f, ~%.0f\u00d7 across the 5\u201395%% range): slow free tracks fall below the band, so the sub-diffusive %% is an upper bound.",
              sdl, 10^(2 * 1.645 * sdl))))
        else NULL
      },
      {
        gn <- trk_disp_geom_null()
        if (is.null(gn))
          htmltools::div(style = "color:#B8860B; margin-top:5px; font-size:0.95em;",
            htmltools::HTML("&#9888; Infinite-plane null: free diffusion on a curved cell already reads mildly sub-diffusive. <b>Step 2</b> separates the two. A low directed % does not rule out transport."))
        else {
          gc <- gn$gclass
          fc <- mean(gc == "confined"); fg2 <- mean(gc == "free-on-geometry"); fd <- mean(gc == "directed")
          geom_lab <- c(sphere = "sphere", rod = "spherocylinder (rod)", ellipsoid = "ellipsoid",
                        membrane = "membrane plane", free = "free (unbounded)")[gn$geom]
          unbounded <- gn$geom %in% c("free", "membrane")
          if (!unbounded)
            geom_lab <- paste0(geom_lab, if (identical(gn$gmode, "volume"))
                               " - inside the volume" else " - on the surface")
          htmltools::tagList(
            htmltools::div(style = "margin-top:7px; padding-top:5px; border-top:1px solid #ddd;",
              htmltools::HTML(sprintf(
                "<b>Step 2, free diffusion on %s</b> (%s; surface D = %.3g \u00b5m\u00b2/s): <b style='color:#2c7fb8;'>%.0f%% confined</b> · <b style='color:#1B9E5A;'>%.0f%% free-on-geometry</b> · <b style='color:#e31a1c;'>%.0f%% directed</b>.",
                geom_lab, gn$dims, gn$D_surf / 1e3, 100 * fc, 100 * fg2, 100 * fd))),
            htmltools::div(style = "color:#555; margin-top:4px;",
              htmltools::HTML(
                if (unbounded)
                  sprintf("Unbounded geometry, so this matches the infinite-plane null (%.0f%% confined). Pick a closed surface to test finite size.", 100 * fc)
                else sprintf(
                  "Of the %.0f%% called sub-diffusive on the plane, geometry keeps <b>%.0f%% as truly confined</b>%s",
                  100 * f_sub, 100 * fc,
                  if (fc <= 0.08) ", at the ~5% free baseline: consistent with free surface diffusion."
                  else " - the excess over the ~5% free baseline."))),
            htmltools::div(style = "color:#555; margin-top:4px; font-size:0.94em;",
              htmltools::HTML(
                if ((attr(res, "d_obs") %||% 2L) == 2L)
                  "2D data: the geometry is simulated in 3D and projected to the image plane."
                else "3D data: the geometry null is evaluated in 3D, no projection.")),
            htmltools::div(style = "color:#B8860B; margin-top:4px; font-size:0.94em;",
              htmltools::HTML("&#9888; The null reflects at the surface but does not model dark-cycle survivorship."))
          )
        }
      }
    )
  })

  output$trk_dl_disp_csv <- downloadHandler(
    filename = function() paste0("displacement_vs_D_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) write.csv(trk_disp_bundle(), file, row.names = FALSE)
  )

  # ============================================================================
  # dir_bundle() - DIRECTIONALITY (is the motion persistent, random, or reversing?)
  # ============================================================================
  #
  # GOAL
  #   Free diffusion has no memory of direction: each step is independent, so the
  #   angle between consecutive steps is completely random. Any DEVIATION from
  #   that is informative - persistence points to directed transport, systematic
  #   reversal points to confinement or tethering.
  #
  # THE MATH
  #   (a) TURN ANGLE between consecutive displacement vectors:
  #
  #          cos θ = (Δr_i · Δr_{i+1}) / (|Δr_i| · |Δr_{i+1}|),   θ = acos(...)
  #
  #       Reading the distribution of θ:
  #          flat/uniform  -> free Brownian (no preferred direction)
  #          peak near 0°  -> persistent / directed motion
  #          peak near 180 -> anti-persistent: confinement, tethering, OR noise
  #
  #       Crucially this null is Δt-INDEPENDENT: the angle between two
  #       independent isotropic vectors is uniform no matter how long each step
  #       lasted. So the turn-angle test needs no Δt normalisation at all.
  #
  #   (b) VELOCITY AUTOCORRELATION over unit step vectors v̂ = Δr/|Δr|:
  #
  #          C(τ) = < v̂(t) · v̂(t + τ) >
  #
  #       C > 0 and slowly decaying = persistent motion; C ~ 0 = free diffusion.
  #
  #   (c) THE NOISE FLOOR (why a negative first lag is EXPECTED, not a finding).
  #       A pure-noise step has RMS length sqrt(2d)·σ ( = 2σ in 2D), because each
  #       of the two localizations contributes σ per axis. More importantly the
  #       SHARED localization error at the midpoint enters the two neighbouring
  #       steps with OPPOSITE sign - the very same effect the CVE exploits. That
  #       makes lag-1 anti-correlated even for a perfectly free molecule:
  #
  #          C(first lag) ≈ − σ² / MSD(one cycle)
  #
  #       So a small negative first lag is the signature of measurement noise.
  #       Only a value substantially more negative than the following lags is
  #       evidence of real anti-persistence. The fraction of steps below the
  #       noise floor is reported so you can judge how much is noise.
  #
  # THE MINFLUX ADAPTATION
  #   The autocorrelation lag must be ELAPSED TIME, not step index - on
  #   irregularly sampled data a step-index lag would draw unevenly spaced points
  #   as if evenly spaced and smear the decay. So pairs are binned by the
  #   difference of their CUMULATIVE CYCLE INDEX, and the axis is cycles × t_c.
  #   Steps spanning more than `trk_dir_max_dark` dark cycles are dropped
  #   entirely: their direction is unreliable and a long gap may even bridge two
  #   different molecules.
  #
  # OBJECTS IT CALLS (inputs)
  #   trk(), trk_cycle_time(), trk_sigma_loc_eff()
  #   inputs: trk_dir_run (Run trigger), trk_dir_min_pts, trk_dir_max_dark
  #
  # STEP BY STEP
  #   1. Split into tracks; per track build steps, their Δt and cycle counts.
  #   2. Mark "usable" steps: non-zero length AND <= max_dark dark cycles.
  #   3. Turn angles: only between two CONSECUTIVE usable steps (they share a
  #      localization, so the angle is meaningful).
  #   4. Autocorrelation: for each pair of usable steps in a track, bin the dot
  #      product of their unit vectors by the elapsed cycle lag between them.
  #   5. Count steps below the noise floor to quantify how much of the signal is
  #      localization error rather than motion.
  #
  # RETURNS a list:
  #   angles ............. all turn angles, degrees (pooled across tracks)
  #   ac[lag], lags ...... autocorrelation per integer cycle lag
  #   lag_time, have_cyc . the lag axis in ms (falls back to step index if the
  #                        cycle time is unknown - and says so on the plot)
  #   t_cycle, max_dark, r_noise, sigma_loc
  #   frac_below_noise ... fraction of steps shorter than the noise floor
  #   frac_dropped_dark .. fraction excluded for spanning too many dark cycles
  #
  # CAVEAT: single short tracks carry little power - these are pooled statistics.
  # ============================================================================
  dir_bundle <- eventReactive(input$trk_dir_run, {
    req(isTRUE(input$trk_show_dir))
    df <- trk(); req(nrow(df) > 0)
    min_pts  <- max(4L, as.integer(input$trk_dir_min_pts))
    max_dark <- max(0L, as.integer(input$trk_dir_max_dark %||% 2L))

    df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
    df <- df[order(df$tid, df$tim), ]
    sp <- split(df, df$tid)

    # Cycle time: the autocorrelation lag is measured in ELAPSED TIME, expressed as
    # integer cycle counts (like the MSD), not raw step index. On irregular MINFLUX
    # sampling a "step" is not a fixed Δt, so a step-index lag would misrepresent
    # unevenly-spaced points as evenly spaced. When the cycle time is unknown we fall
    # back to step-index lags and say so.
    t_c       <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    have_cyc  <- is.finite(t_c) && t_c > 0
    max_lag   <- 20L                       # cycles (or steps in fallback)

    # Pre-allocate angle storage as a list then unlist once at the end
    angle_list <- vector("list", length(sp))
    ac_sums    <- rep(0,  max_lag)
    ac_n       <- rep(0L, max_lag)
    # Localization-noise floor on a 2D step magnitude: RMS pure-noise step = 2·σ_loc.
    # Steps below it have noise-dominated directions (inflate the 180° turn-angle
    # peak and the lag-1 anti-correlation). Count them to quantify how much of the
    # directionality signal is noise. Uses the app-wide σ_loc.
    sig_dir <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    r_noise <- if (is.finite(sig_dir) && sig_dir > 0) 2 * sig_dir else NA_real_
    n_steps_tot <- 0L; n_below <- 0L; n_dropped_dark <- 0L

    withProgress(message = "Computing directionality…", value = 0, {
    for (ii in seq_along(sp)) {
      if (ii %% 50L == 0L) incProgress(50 / length(sp))
      tb <- sp[[ii]]
      if (nrow(tb) < min_pts) next
      dx  <- diff(tb$x); dy <- diff(tb$y)
      dtt <- diff(tb$tim)
      if (length(dx) < 2) next

      mag <- sqrt(dx^2 + dy^2)
      # Per-step dark-cycle count. A "usable" step spans ≤ max_dark dark cycles
      # (integrated = 0 dark) and has non-zero length. Steps across long gaps are
      # dropped: the direction over a big gap is unreliable and may bridge a merge
      # with another molecule.
      n_cyc  <- if (have_cyc) pmax(1L, as.integer(round(dtt / t_c))) else rep(1L, length(dtt))
      usable <- mag > 0 & is.finite(n_cyc) & ((n_cyc - 1L) <= max_dark)
      n_steps_tot    <- n_steps_tot + length(mag)
      n_dropped_dark <- n_dropped_dark + sum(mag > 0 & (n_cyc - 1L) > max_dark, na.rm = TRUE)
      if (is.finite(r_noise)) n_below <- n_below + sum(mag < r_noise & usable, na.rm = TRUE)

      # --- Turn angles: only between two CONSECUTIVE usable steps (share a loc) ---
      ok_pair <- usable[-length(usable)] & usable[-1]
      if (any(ok_pair)) {
        dot   <- dx[-length(dx)] * dx[-1] + dy[-length(dy)] * dy[-1]
        denom <- mag[-length(mag)] * mag[-1]
        cos_a <- dot[ok_pair] / denom[ok_pair]
        cos_a <- pmax(-1, pmin(1, cos_a))
        angle_list[[ii]] <- acos(cos_a) * 180 / pi
      }

      # --- Velocity autocorrelation, binned by ELAPSED-TIME (cycle-index) lag ---
      # Each usable step's reference position on the track's cumulative cycle
      # timeline is the cycle index of its start localization. Pairs of usable steps
      # are accumulated into the bin = difference of those cycle indices, so a lag of
      # k means k cycles of real time, however many localizations were skipped.
      cyc_start <- c(0L, cumsum(n_cyc))          # per-loc cumulative cycle index
      cyc_start <- cyc_start[-length(cyc_start)] # start-of-step index (length = #steps)
      ux <- ifelse(usable, dx / mag, NA_real_)
      uy <- ifelse(usable, dy / mag, NA_real_)
      idx_u <- which(usable)
      for (a in idx_u) {
        ca <- cyc_start[a]
        for (b in idx_u[idx_u > a]) {
          lag <- cyc_start[b] - ca
          if (lag < 1L) next
          if (lag > max_lag) break            # cyc_start is increasing → no closer pairs beyond
          ac_sums[lag] <- ac_sums[lag] + ux[a]*ux[b] + uy[a]*uy[b]
          ac_n[lag]    <- ac_n[lag]    + 1L
        }
      }
    }

    all_angles <- unlist(angle_list, use.names = FALSE)

    shiny::validate(shiny::need(length(all_angles) >= 10,
      "Too few usable steps - lower 'Min points per TID', raise 'Max dark cycles', or relax filters."))

    ac_vals  <- ifelse(ac_n > 0, ac_sums / ac_n, NA_real_)
    lag_time <- if (have_cyc) seq_len(max_lag) * t_c else NA_real_   # ms, real elapsed time
    list(angles = all_angles, ac = ac_vals, lags = seq_len(max_lag),
         lag_time = lag_time, have_cyc = have_cyc, t_cycle = t_c, max_dark = max_dark,
         r_noise = r_noise, sigma_loc = if (is.finite(r_noise)) sig_dir else NA_real_,
         frac_below_noise = if (n_steps_tot > 0) n_below / n_steps_tot else NA_real_,
         frac_dropped_dark = if (n_steps_tot > 0) n_dropped_dark / n_steps_tot else NA_real_)
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

    # x-axis is real elapsed time (cycle lag × cycle time) when the cycle time is
    # known; step-index only as a labelled fallback.
    xv    <- if (isTRUE(bd$have_cyc)) bd$lag_time else bd$lags
    xlab  <- if (isTRUE(bd$have_cyc)) "Time lag (ms)" else "Lag (steps - cycle time unknown)"
    plot_ly() %>%
      add_lines(x = xv, y = bd$ac,
                line = list(color = "#D06030", width = 2)) %>%
      add_markers(x = xv, y = bd$ac,
                  marker = list(color = "#D06030", size = 5)) %>%
      add_lines(x = range(xv, na.rm = TRUE), y = c(0, 0),
                line = list(color = fg_col, width = 1, dash = "dot"),
                showlegend = FALSE) %>%
      layout(
        xaxis  = list(title = xlab, color = fg_col),
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
    ac1_lbl <- if (isTRUE(bd$have_cyc)) sprintf("AC(%.0f ms)", bd$t_cycle) else "AC(lag 1 step)"
    dark_note <- if (is.finite(bd$frac_dropped_dark %||% NA_real_) && (bd$frac_dropped_dark %||% 0) > 0)
      tags$span(style = "color:#888;", tags$br(),
        sprintf("%.0f%% of steps span more than %d dark cycle(s) and were left out.",
                100 * bd$frac_dropped_dark, bd$max_dark))
    else NULL
    noise_note <- if (is.finite(bd$frac_below_noise %||% NA_real_))
      tags$span(style = "color:#b26a00;", tags$br(),
        sprintf("⚠ %.0f%% of steps are below the localization-noise floor (2·", 100 * bd$frac_below_noise),
        sigloc(),
        sprintf(" = %.1f nm) - their directions are noise-dominated and inflate the backward (180°) peak / lag-1 anti-correlation.",
                bd$r_noise))
    else
      tags$span(style = "color:#888;", tags$br(),
        "Set ", sigloc(), " (run the MSD, or enter it in Acquisition Parameters) to see how much of the turn-angle signal is localization noise.")
    tags$div(style = "font-size:0.88em; margin:4px 0;",
      tags$b(sprintf("n = %d turn angles", length(bd$angles))),
      tags$span(sprintf(
        "  |  median angle = %.1f\u00B0  |  forward (<60\u00B0): %.0f%%  |  backward (>120\u00B0): %.0f%%  |  %s = %s",
        med_ang, pct_fwd, pct_back, ac1_lbl, if (is.na(ac1)) "-" else ac1
      )),
      dark_note,
      noise_note
    )
  })

  output$trk_dl_dir_csv <- downloadHandler(
    filename = function() paste0("directionality_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      bd <- dir_bundle()
      # Two sheets packed into one CSV with a blank-row separator
      lag_ms <- if (isTRUE(bd$have_cyc)) bd$lag_time else rep(NA_real_, length(bd$lags))
      ac_df  <- data.frame(lag_cycles = bd$lags, lag_ms = lag_ms, velocity_autocorr = bd$ac)
      ang_df <- data.frame(turn_angle_deg = bd$angles)
      out <- rbind(
        cbind(sheet = "autocorrelation", ac_df,
              data.frame(turn_angle_deg = rep(NA, nrow(ac_df)))),
        data.frame(sheet = NA, lag_cycles = NA, lag_ms = NA, velocity_autocorr = NA,
                   turn_angle_deg = NA),
        cbind(sheet = "angles", lag_cycles = NA, lag_ms = NA, velocity_autocorr = NA,
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
  # stays - inflating v = dr/span into a spurious spike (which otherwise reads as a
  # false "export at the very end"). Requiring the window to span >= frac·W drops
  # those unreliable end windows. Cost: the last <frac·W of each track is not
  # measured - you cannot estimate a velocity over near-zero time.
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
  # index - so MINFLUX's irregular Δt / dark cycles don't bias α (a window straddling a
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

  # ==========================================================================
  # trk_export_summary() - DIRECTED TRANSPORT ("export") DETECTION, per track
  # ==========================================================================
  #
  # GOAL
  #   Detect when a molecule stops diffusing and starts being CARRIED - a
  #   sustained, directional excursion (e.g. motor-driven transport, or a
  #   molecule leaving a compartment). Produces one summary row per track: did
  #   it export, how fast, and WHERE along the track it happened.
  #
  # THE CENTRAL DIFFICULTY: DIFFUSION LOOKS FAST TOO
  #   You cannot detect transport from speed. From <r²> = 2dDΔt, the per-step
  #   speed of a purely diffusing molecule is
  #
  #        v_step = r/Δt ~ sqrt(2dD/Δt)    ->  DIVERGES as Δt -> 0
  #
  #   so a freely diffusing molecule shows an arbitrarily high instantaneous
  #   "speed" if you sample it finely enough (header §2b: Brownian motion has no
  #   instantaneous velocity). Any threshold on path speed therefore measures
  #   your sampling rate, not transport.
  #
  #   THE FIX - NET DISPLACEMENT OVER A WINDOW. Diffusion is unbiased: over a
  #   window W its net displacement grows only as sqrt(W), because the steps
  #   cancel. Directed motion is biased: its net displacement grows LINEARLY
  #   with W. So the drift velocity
  #
  #        v_drift(t) = |r(t+W) − r(t)| / W
  #
  #   stays small for diffusion but stays high for transport, and the two
  #   separate more strongly the longer the window. That is what is thresholded.
  #
  # WHAT COUNTS AS AN EXPORT (three guards against false positives)
  #   1. THRESHOLD = max(user or auto threshold, localization-noise floor).
  #      Auto = median + 3·MAD (robust, so a few fast tracks don't set the bar).
  #      Noise floor = sqrt(2)·σ_loc / W - the apparent drift velocity produced
  #      by localization error alone over that window.
  #   2. SUSTAINED RUN of >= min_run consecutive windows above threshold. A
  #      single spike is noise; transport persists.
  #   3. Tracks too short for a full window are skipped rather than guessed at.
  #
  #   onset_frac = where the export starts along the track (0 = start, 1 = end).
  #   This is the diagnostic that reveals the "exports cluster at track ends"
  #   pattern - which is usually an ARTEFACT (photobleaching / loss of lock at
  #   the end of a track), not biology. Check this histogram before believing
  #   an export result.
  #
  # OBJECTS IT CALLS (inputs)
  #   trk(), trk_speed_sigma_eff() (σ_loc for the noise floor),
  #   trk_drift_velocity_one() (the per-track drift-velocity curve)
  #   inputs: trk_speed_drift_win (W), trk_speed_drift_minpts,
  #           trk_speed_min_pts, threshold + min-run settings
  #
  # RETURNS a data.frame, one row per qualifying track:
  #   tid, has_export, v_max / v_export (µm/s), export_pos_frac (onset),
  #   run length, and the threshold actually used.
  #
  # CAVEAT: this reports SUSTAINED DIRECTIONAL BIAS. Confirm against the MSD
  # directed model (which fits a velocity) and the HMM/per-state MSD before
  # calling it active transport - confinement release and flow can look similar.
  # ==========================================================================
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
    shiny::validate(shiny::need(length(tids) > 0, "No tracks with enough points - lower 'Min points per TID'."))

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
      # localization-noise floor - i.e. confined/immobile at this timescale. High
      # immobile_frac + several n_events = the hop-diffusion (dwell → hop → dwell) look.
      immobile_frac <- if (n_int > 0) sum(interior & v <= max(vfloor, 0), na.rm = TRUE) / n_int else NA_real_
      onset_ms <- NA_real_; onset_frac <- NA_real_; peak_frac <- NA_real_
      export_v <- NA_real_; export_dur <- NA_real_
      # Per-event corroboration (artifact vs signal). A GENUINE directed excursion is
      # made of INTEGRATED localizations with normal photon counts - the molecule is
      # present and tracked throughout. A tracker SEARCH after losing the molecule has
      # the same gross signature (fast, outward, track ends) but its localizations are
      # dark-cycle-dominated and photon-starved. These fields come from dark cycles and
      # photons, NOT from drift velocity, so they do not restate what defined the event.
      ev_integ <- NA_real_; ev_eco_ratio <- NA_real_; ev_cfr <- NA_real_
      ev_nloc  <- NA_integer_; ev_at_end <- NA; ev_support <- NA_character_
      t0_abs <- suppressWarnings(min(d$tim, na.rm = TRUE))
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

        # Corroborate the STRONGEST event (the one export_v / export_pos_frac report):
        # find the run containing the peak window, widen by W/2 to cover the
        # localizations the windows were actually built from, then score that span.
        wrun <- which(seg_ok)
        prun <- wrun[which(peak_idx >= starts[wrun] & peak_idx <= ends[wrun])[1]]
        if (!is.na(prun)) {
          ev_lo <- dv$t_ms[starts[prun]] - drift_win / 2
          ev_hi <- dv$t_ms[ends[prun]]   + drift_win / 2
          rel   <- d$tim - t0_abs
          insid <- is.finite(rel) & rel >= ev_lo & rel <= ev_hi
          ev_nloc <- sum(insid)
          if (ev_nloc > 0) {
            if ("n_dark" %in% names(d))
              ev_integ <- mean(d$n_dark[insid] == 0L, na.rm = TRUE)
            if ("eco" %in% names(d)) {
              e_ev <- stats::median(d$eco[insid], na.rm = TRUE)
              e_tr <- stats::median(d$eco,        na.rm = TRUE)   # track's own baseline
              if (is.finite(e_ev) && is.finite(e_tr) && e_tr > 0) ev_eco_ratio <- e_ev / e_tr
            }
            if ("cfr" %in% names(d)) ev_cfr <- stats::median(d$cfr[insid], na.rm = TRUE)
          }
          ev_at_end <- isTRUE((tmax - dv$t_ms[ends[prun]]) <= drift_win)
          ev_support <-
            if (is.finite(ev_integ) && ev_integ < 0.5) "dark-dominated"
            else if (is.finite(ev_eco_ratio) && ev_eco_ratio < 0.6) "low-photon"
            else if (is.finite(ev_integ) && ev_integ >= 0.5) "supported"
            else "unscored"
        }
      }
      # Track-level TAIL quality: integrated fraction over the final W of the track.
      # Directly answers "is the outward tail the molecule, or the tracker hunting?"
      tail_integ <- NA_real_
      if ("n_dark" %in% names(d) && is.finite(t0_abs)) {
        rel <- d$tim - t0_abs
        rmx <- suppressWarnings(max(rel, na.rm = TRUE))
        tl  <- is.finite(rel) & rel >= (rmx - drift_win)
        if (any(tl)) tail_integ <- mean(d$n_dark[tl] == 0L, na.rm = TRUE)
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
        # --- corroboration for the strongest event (independent of drift velocity) ---
        event_support     = ev_support,      # supported | dark-dominated | low-photon
        event_integ_frac  = ev_integ,        # 1 = every loc in the event was integrated
        event_eco_ratio   = ev_eco_ratio,    # event photons / this track's own median
        event_cfr_med     = ev_cfr,
        event_n_loc       = ev_nloc,
        event_at_track_end = ev_at_end,
        tail_integ_frac   = tail_integ,      # integrated fraction in the track's final W
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
    # Rank by segment strength so the picker doubles as a follow-up list:
    # flagged tracks first, then by peak drift velocity (descending).
    s <- s[order(-as.integer(s$has_export %in% TRUE), -s$v_drift_peak), , drop = FALSE]
    labs <- sprintf("%s%s  (peak %.2g µm/s)",
                    ifelse(s$has_export %in% TRUE, "★ ", "  "),
                    s$tid, s$v_drift_peak)
    choices <- stats::setNames(s$tid, labs)
    sel <- if (!is.null(input$trk_speed_tid) && input$trk_speed_tid %in% s$tid)
             input$trk_speed_tid else s$tid[1]
    selectInput("trk_speed_tid", "Track (TID) - ★ = candidate directed segment, ranked by peak drift",
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
    shiny::validate(shiny::need(length(tids) > 0, "No tracks with enough points - lower 'Min points per TID'."))
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

    # Path speed (per-step) - kept for reference; it is dominated by diffusion+noise
    if (show_path) {
      fig <- fig %>%
        add_trace(x = spd$t_ms, y = spd$speed_um_s, type = "scattergl", mode = "markers",
                  marker = list(size = 4, color = "rgba(140,140,140,0.30)"),
                  name = "path speed (raw)", hoverinfo = "x+y") %>%
        add_trace(x = spd$t_ms, y = spd$speed_med, type = "scattergl", mode = "lines",
                  line = list(color = "rgba(120,120,120,0.85)", width = 1.5, dash = "dot"),
                  name = "path speed (median)")
    }

    # Rolling drift (transport) velocity - the segment-detecting curve
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
    # Corroboration breakdown: how many candidate events are actually supported by
    # evidence independent of the drift threshold that found them.
    sup  <- if (nrow(exp_rows)) sum(exp_rows$event_support %in% "supported")      else 0L
    dkd  <- if (nrow(exp_rows)) sum(exp_rows$event_support %in% "dark-dominated") else 0L
    lowp <- if (nrow(exp_rows)) sum(exp_rows$event_support %in% "low-photon")     else 0L
    med_ti <- stats::median(s$tail_integ_frac, na.rm = TRUE)
    htmltools::div(style = "font-size:0.9em; line-height:1.5;",
      htmltools::HTML(sprintf(
        "<b>%d</b> tracks analysed &middot; <b>%d (%.0f%%)</b> with a candidate directed segment%s<br>%s%smedian peak drift = <b>%.3g µm/s</b>",
        n, ne, if (n > 0) 100 * ne / n else 0,
        if (n_multi > 0) sprintf(" &middot; <b>%d</b> with multiple episodes (hop/multi-step)", n_multi) else "",
        if (is.finite(med_expv))  sprintf("median segment velocity = <b>%.3g µm/s</b> &middot; ", med_expv) else "",
        if (is.finite(med_onset)) sprintf("median segment at <b>%.0f%%</b> of track length &middot; ", 100 * med_onset) else "",
        med_peak)),
      # Corroboration = evidence from dark cycles and photons, i.e. NOT from the drift
      # threshold that found the segment. The circularity argument itself is in the
      # panel's Method box, not reprinted on every run.
      if (ne > 0) htmltools::HTML(sprintf(
        "<div style='margin-top:6px;'>Corroboration (dark cycles &amp; photons, independent of drift): <b style='color:#1B9E5A;'>%d supported</b> &middot; <b style='color:#C0392B;'>%d dark-dominated</b> &middot; <b style='color:#B8860B;'>%d low-photon</b></div>",
        sup, dkd, lowp)),
      if (is.finite(med_ti)) htmltools::HTML(sprintf(
        "<div style='margin-top:4px; color:#555;'>Median tail integrated fraction = <b>%.2f</b>%s</div>",
        med_ti, if (med_ti < 0.5) " - track ends are search-dominated" else "")),
      htmltools::HTML(
        "<div style='margin-top:6px; color:#B8860B; font-size:0.95em;'>&#9888; Segment velocity is <b>not independent evidence</b>: the flag that found it was a speed. Judge on corroboration and a matched negative control.</div>")
    )
  })

  # Distribution of a per-track metric (the "histogram of speeds" across tracks)
  output$trk_export_hist <- renderPlotly({
    req(isTRUE(input$trk_show_speed))
    s <- trk_export_summary()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    metric <- input$trk_export_hist_metric %||% "v_drift_peak"
    lab <- switch(metric,
      v_drift_peak     = "peak drift velocity (µm/s)",
      export_v         = "segment velocity (µm/s)",
      directed_frac    = "fraction of track time directed",
      immobile_frac    = "fraction of track time immobile/confined",
      n_events         = "number of transport episodes per track",
      alpha_max        = "max local α",
      event_integ_frac = "integrated fraction within the event (1 = all integrated)",
      event_eco_ratio  = "event photons / track median photons",
      tail_integ_frac  = "integrated fraction in the track's final window",
      "value")
    # Event-scoped metrics exist only for flagged tracks; the rest span all tracks.
    event_only <- c("export_v", "event_integ_frac", "event_eco_ratio")
    vals_all <- s[[metric]]
    exp_mask <- s$has_export %in% TRUE
    fig <- plot_ly()
    if (!(metric %in% event_only)) {
      # unflagged vs flagged stacked so the split is visible
      ve <- vals_all[exp_mask];  vn <- vals_all[!exp_mask]
      ve <- ve[is.finite(ve)];   vn <- vn[is.finite(vn)]
      if (length(vn)) fig <- fig %>% add_histogram(x = vn, name = "unflagged",
        marker = list(color = "rgba(150,150,150,0.55)"), nbinsx = 30)
      if (length(ve)) fig <- fig %>% add_histogram(x = ve, name = "flagged",
        marker = list(color = "rgba(27,158,90,0.75)"), nbinsx = 30)
      fig <- fig %>% layout(barmode = "stack")
    } else {
      v <- vals_all[exp_mask]; v <- v[is.finite(v)]
      shiny::validate(shiny::need(length(v) > 0, "No flagged tracks yet - loosen the threshold / min-run."))
      fig <- fig %>% add_histogram(x = v, name = "flagged",
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
  # ── Is a late-track excursion real, or the tracker searching? ────────────────
  # WHERE an event sits along the track cannot answer that on its own: a molecule
  # leaving and a tracker losing its molecule both produce fast outward motion at the
  # end. What separates them is whether the localizations during the event were still
  # being found on schedule. So plot position (x) against the event's integrated
  # fraction (y): a molecule that is genuinely present keeps being localized and sits
  # high; a tracker searching accumulates dark cycles and sinks, typically at the far
  # right where tracks end. y is independent of the drift threshold that found the
  # event, which is what makes it evidence rather than restatement.
  output$trk_export_corrob <- renderPlotly({
    req(isTRUE(input$trk_show_speed))
    s <- trk_export_summary()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    e <- s[s$has_export %in% TRUE, , drop = FALSE]
    e <- e[is.finite(e$export_onset_frac) & is.finite(e$event_integ_frac), , drop = FALSE]
    shiny::validate(shiny::need(nrow(e) > 0,
      "No candidate segments yet - press 'Run speed & transport', or relax the threshold."))
    pal <- c(supported = "#1B9E5A", `dark-dominated` = "#C0392B", `low-photon` = "#B8860B")
    fig <- plot_ly()
    for (lv in names(pal)) {
      i <- which(e$event_support %in% lv)
      if (!length(i)) next
      fig <- fig %>% add_trace(
        x = e$export_onset_frac[i], y = e$event_integ_frac[i],
        type = "scatter", mode = "markers",
        marker = list(size = 8, color = pal[[lv]], opacity = 0.75,
                      line = list(width = 0.4, color = "rgba(0,0,0,0.35)")),
        name = sprintf("%s (%d)", lv, length(i)),
        # hoverinfo="text" rather than a hovertemplate: no template parsing, so the
        # tooltip cannot fall back to printing the placeholder itself.
        text = sprintf("TID %s<br>starts at %.0f%% of track<br>integrated fraction %.2f<br>photons %.2f x track median<br>%s",
                       as.character(e$tid[i]), 100 * e$export_onset_frac[i],
                       e$event_integ_frac[i], e$event_eco_ratio[i],
                       as.character(e$event_support[i])),
        hoverinfo = "text")
    }
    fig %>% layout(
      xaxis = list(title = "where the segment starts along the track (0 = start, 1 = end)",
                   range = c(-0.02, 1.02), color = fg),
      yaxis = list(title = "integrated fraction during the segment",
                   range = c(-0.02, 1.02), color = fg),
      shapes = list(list(type = "line", x0 = -0.02, x1 = 1.02, y0 = 0.5, y1 = 0.5,
                         line = list(color = "rgba(120,120,120,0.45)", dash = "dot", width = 1))),
      annotations = list(list(x = 0.02, y = 0.52, xref = "x", yref = "y",
        text = "above: molecule still being found · below: dark-cycle dominated",
        showarrow = FALSE, xanchor = "left", yanchor = "bottom",
        font = list(size = 9, color = "rgba(120,120,120,0.9)"))),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(font = list(color = fg, size = 10), orientation = "h",
                    x = 0, xanchor = "left", y = -0.24, yanchor = "top"),
      margin = list(l = 60, r = 20, t = 10, b = 90),
      uirevision = "export_corrob")
  })

  output$trk_export_onset_hist <- renderPlotly({
    req(isTRUE(input$trk_show_speed))
    s <- trk_export_summary()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
    v <- s$export_pos_frac[s$has_export %in% TRUE]; v <- v[is.finite(v)]
    shiny::validate(shiny::need(length(v) > 0, "No flagged tracks yet - loosen the threshold / min-run."))
    plot_ly(x = v, type = "histogram", nbinsx = 20,
            marker = list(color = "rgba(123,79,190,0.75)")) %>%
      layout(
        xaxis = list(title = "segment position (0 = start … 1 = end of track)",
                     range = c(0, 1), color = fg),
        yaxis = list(title = "flagged tracks", color = fg),
        paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
        # Neutral wording: a late pile-up is equally "molecule left" and "tracker lost it".
        title = list(text = sprintf("Segment location (n=%d) - pile-up near 1 is ambiguous: check corroboration", length(v)),
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

  # ============================================================================
  # hmm_bundle() - DIFFUSION-STATE SEGMENTATION (continuous-time HMM, est_hmm.R)
  # ============================================================================
  #
  # GOAL
  #   Answer "does this population switch between diffusive states, and how fast?"
  #   with K states (D_k) and switching RATES estimated by maximum likelihood over
  #   all tracks at once. The estimator itself lives in est_hmm.R (pure functions,
  #   tested against simulated ground truth); this reactive only gathers inputs,
  #   calls it and packages the result for the plots.
  #
  # THE MODEL (full statement in the est_hmm.R header)
  #   Hidden state per step, a continuous-time Markov chain with generator Q, so a
  #   step of any Δt has transition matrix expm(Q·Δt): dark cycles are handled
  #   exactly. Emission: r² over x,y is Gamma(1, 2·v) with the per-axis variance
  #       v = 2·D_s·(Δt − 2·R·t_c) + σ_j² + σ_{j+1}²
  #   i.e. diffusion, motion blur (trk_blur()) and the noise of both localizations.
  #   D is therefore noise- and blur-corrected by construction.
  #
  # WHAT CHANGED AND WHY (2026-09-17)
  #   The old version fitted a per-step matrix by Baum-Welch and counted switches
  #   along the Viterbi path. Path counting absorbs short dwells: on simulated data
  #   with dwells of 10 / 20 ms it read 15.6 / 30.2 ms (+50%), while this MLE read
  #   9.9 / 18.8 ms. The Viterbi path is still computed, but ONLY to colour tracks
  #   and to feed the per-state MSD.
  #
  # OBJECTS IT CALLS (inputs)
  #   trk() ................. analysis-ready data (tid, x, y, tim in ms, eco)
  #   trk_cycle_time() ...... t_c, for the blur offset and starting rates
  #   trk_sigma_loc_eff() / trk_sigma_loc_src() ... the noise model: per-localization
  #                           photon σ when the photon model is active, else one σ
  #   trk_cve_bundle() ...... optional, supplies the emitter offset x for the photon σ
  #   inputs: trk_hmm_run, trk_hmm_n_states, trk_hmm_min_pts, trk_hmm_n_ex,
  #           trk_hmm_boot (bootstrap resamples, 0 = off)
  #
  # RETURNS a list (fields read by the plots, the summary and seg_msd_rv):
  #   K, D_vals[K] (nm²/ms), fracs[K] (stationary occupancy), P_cycle (K×K
  #   transition probabilities over one cycle), Qrate (per ms), dwell_ms[K],
  #   fit (the full est_hmm result incl. CIs), boot (bootstrap CIs or NULL),
  #   seg_list (per-track data.frame(tid, x, y, tim, state)), n_tracks, n_steps,
  #   dt_med_ms, t_cycle_ms, sigma_loc_nm, sigma_mode, D_noise, noise_state[K],
  #   resolv (per-state flags from the validated range), bic_other (BIC at K±1)
  # ============================================================================
  hmm_bundle <- eventReactive(input$trk_hmm_run, {
    req(isTRUE(input$trk_show_hmm))
    df <- trk(); req(nrow(df) > 0)
    K       <- max(2L, min(4L, as.integer(input$trk_hmm_n_states)))
    min_pts <- max(4L, as.integer(input$trk_hmm_min_pts))
    B       <- max(0L, min(200L, as.integer(input$trk_hmm_boot %||% 0L)))
    t_c     <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    shiny::validate(shiny::need(is.finite(t_c) && t_c > 0,
      "Cycle time unknown - run Data Quality & Calibration first."))
    R_abs   <- trk_blur()$R_abs                         # blur offset in ms, shared (est_blur.R)

    # ── Noise model ──────────────────────────────────────────────────────────
    # The emission needs each step's noise variance. With the photon model active,
    # every localization gets its own σ from its photon count (the same model the
    # rest of the tab uses); otherwise one σ for all. Without any σ the slow state
    # would absorb the noise and read too fast, so that case is refused.
    sig_eff <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    src     <- tryCatch(trk_sigma_loc_src(), error = function(e) "not set")
    # If Data Quality & Calibration has not published a σ yet, derive the photon-budget
    # value here from the same model rather than refusing: every input it needs (ECO,
    # L_TCP, the prefactor) is already on this tab.
    if (!is.finite(sig_eff) || sig_eff <= 0) {
      if ("eco" %in% names(df) && any(is.finite(df$eco) & df$eco > 0)) {
        L0 <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
        if (!is.finite(L0) || L0 <= 0) L0 <- DMAX_L_TCP_DEFAULT
        n0 <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
        if (!is.finite(n0) || n0 <= 0) n0 <- TCP_BEAM_FWHM_NM
        a0 <- .a_in_use()
        ok0 <- is.finite(df$eco) & df$eco > 0
        sig_eff <- stats::median(.crlb_doughnut(L0, n0, df$eco[ok0], pref = a0,
                                                pattern = input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT),
                                 na.rm = TRUE)
        src <- "photon budget"
      }
    }
    shiny::validate(shiny::need(is.finite(sig_eff) && sig_eff > 0,
      "Localization precision unknown: the file carries no photon counts, so enter σ_loc on this tab."))
    d_in <- df[, intersect(c("tid", "x", "y", "tim", "eco"), names(df)), drop = FALSE]

    # ── Step budget ──────────────────────────────────────────────────────────
    # Fit time scales with the number of steps, and the optimiser evaluates the
    # likelihood over every step many times. The spread of D across molecules, not
    # the step count, is what limits how precisely D is known, so beyond a budget
    # the extra steps buy time rather than certainty. Whole TRACKS are drawn, never
    # parts of them, so no track is cut in half.
    cap <- suppressWarnings(as.integer(input$trk_hmm_max_steps %||% 0L))
    n_sub_tracks <- NA_integer_
    if (is.finite(cap) && cap > 0) {
      per_tid <- table(d_in$tid)
      if (sum(per_tid - 1L) > cap) {
        set.seed(1L)                                     # same subset on a re-run
        shuf <- sample(names(per_tid))
        keep <- shuf[cumsum(per_tid[shuf] - 1L) <= cap]
        if (length(keep) >= 2L) {
          d_in <- d_in[d_in$tid %in% keep, , drop = FALSE]
          n_sub_tracks <- length(per_tid)
        }
      }
    }

    sigma_mode <- "one value for all localizations"
    if (identical(src, "photon budget") && "eco" %in% names(d_in)) {
      cve <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
      L_tcp  <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
      if (!is.finite(L_tcp) || L_tcp <= 0) L_tcp <- DMAX_L_TCP_DEFAULT
      nu     <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
      if (!is.finite(nu) || nu <= 0) nu <- TCP_BEAM_FWHM_NM
      a_pref <- .a_in_use()
      D_off  <- if (!is.null(cve) && is.finite(cve$D_cve) && cve$D_cve > 0) cve$D_cve else 0
      okc    <- is.finite(d_in$eco) & d_in$eco > 0
      d_in$sig <- NA_real_                               # rows without eco fall back to sig_eff
      if (any(okc)) d_in$sig[okc] <- .crlb_doughnut(L_tcp, nu, d_in$eco[okc], D = D_off, t_loc = t_c,
                                                    pref = a_pref,
                                                    pattern = input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT)
      sigma_mode <- "per localization, from its photon count"
    }

    pr <- hmm_ct_prepare(d_in, ndim = 2L, sigma = sig_eff, min_pts = min_pts)
    shiny::validate(shiny::need(!is.null(pr) && pr$n_steps >= 50L,
      "Too few valid steps - lower 'Min localizations per TID'."))

    fit <- withProgress(message = "Fitting continuous-time HMM…", value = 0.2, {
      f <- hmm_ct_fit(pr, K = K, R_abs = R_abs, t_c = t_c)
      incProgress(0.5)
      f
    })
    shiny::validate(shiny::need(!is.null(fit), "The HMM fit failed - try fewer states."))

    # Neighbouring model sizes, so the summary can show whether the data support K.
    # Two further fits, so it is opt-in: it roughly triples the run time.
    bic_other <- if (!isTRUE(input$trk_hmm_bic)) NULL else
      withProgress(message = "Comparing state counts (BIC)…", value = 0.3, {
      ks <- setdiff(c(K - 1L, K + 1L), 1L); ks <- ks[ks <= 4L]
      stats::setNames(vapply(ks, function(k) {
        f2 <- tryCatch(hmm_ct_fit(pr, K = k, R_abs = R_abs, t_c = t_c, hessian = FALSE,
                                  n_starts = 2L), error = function(e) NULL)
        if (is.null(f2)) NA_real_ else f2$BIC
      }, numeric(1)), ks)
    })

    boot <- if (B >= 10L) withProgress(message = sprintf("Bootstrap over tracks (%d refits)…", B), value = 0,
      hmm_ct_bootstrap(pr$df, fit, B = B, ndim = 2L, sigma = sig_eff, t_c = t_c,
                       min_pts = min_pts, progress = function(p) setProgress(p)))
      else NULL

    # Display-only decode under the fitted model.
    step_states <- hmm_ct_viterbi(pr, fit)
    lab <- hmm_ct_states_by_loc(pr, step_states)
    lab <- lab[is.finite(lab$state), , drop = FALSE]
    seg_list <- lapply(split(lab, lab$tid), function(d)
      data.frame(tid = d$tid, x = d$x, y = d$y, tim = d$tim, state = as.integer(d$state),
                 stringsAsFactors = FALSE))
    # Longest tracks first, so the example panel shows tracks with room to switch.
    seg_list <- seg_list[order(-vapply(seg_list, nrow, integer(1)))]

    dt_med  <- stats::median(pr$DT, na.rm = TRUE)
    D_noise <- sig_eff^2 / dt_med                        # D equivalent of the noise per step
    # Validated range (tests/validation/hmm_recovery_sweep.R, HMM_RESOLV_* in packages.R).
    exp_sw  <- fit$expected_switches
    resolv <- data.frame(
      state        = seq_len(K),
      dwell_cycles = fit$dwell_ms / t_c,
      too_short    = fit$dwell_ms / t_c < HMM_RESOLV_MIN_DWELL_CYC,
      few_switches = rowSums(exp_sw, na.rm = TRUE) < HMM_RESOLV_MIN_SWITCHES,
      stringsAsFactors = FALSE)
    # Two states with the same D differ only in how long they last: the model is
    # describing a non-exponential dwell distribution of ONE diffusivity, not a
    # separate state. Seen on the tracking test file at K = 3 (D 0.026 / 0.026 / 0.49 µm²/s,
    # dwells 17 / 0.26 / 0.17 ms). Flag adjacent states within 5% of each other.
    same_D <- rep(FALSE, K)
    for (k in seq_len(K - 1L)) if (abs(fit$D[k + 1L] / fit$D[k] - 1) < 0.05) same_D[c(k, k + 1L)] <- TRUE
    resolv$same_D <- same_D

    list(K = K, D_vals = fit$D, fracs = fit$pi,
         P_cycle = .expm_small(fit$Q * t_c), Qrate = fit$Q, dwell_ms = fit$dwell_ms,
         fit = fit, boot = boot, bic_other = bic_other, n_tracks_all = n_sub_tracks,
         seg_list = seg_list, n_tracks = fit$n_tracks, n_steps = fit$n_steps,
         dt_med_ms = dt_med, t_cycle_ms = t_c, R_abs = R_abs,
         sigma_loc_nm = sig_eff, sigma_mode = sigma_mode,
         D_noise = D_noise, noise_state = fit$D <= D_noise, resolv = resolv,
         norm_dt = TRUE)
  })

  output$trk_hmm_params_plot <- renderPlotly({
    bd   <- hmm_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"
    state_cols <- c("#C03030","#D08020","#2080C0","#7030A0")[seq_len(bd$K)]
    D_um <- bd$D_vals / 1e3                               # µm²/s
    ci   <- if (!is.null(bd$boot)) bd$boot$D else bd$fit$ci$D   # bootstrap CI preferred
    lo   <- if (!is.null(ci)) ci[1, ] / 1e3 else rep(NA_real_, bd$K)
    hi   <- if (!is.null(ci)) ci[2, ] / 1e3 else rep(NA_real_, bd$K)
    is_noise <- bd$noise_state %||% rep(FALSE, bd$K)
    labels <- sprintf("State %d\n%.0f%% of time", seq_len(bd$K), bd$fracs * 100)
    labels <- ifelse(is_noise, paste0(labels, "\n⚠ at noise level"), labels)
    bar_cols <- ifelse(is_noise, "#9e9e9e", state_cols)
    p <- plot_ly() %>%
      add_bars(
        x = D_um, y = labels, orientation = "h",
        marker = list(color = bar_cols),
        error_x = if (all(is.finite(c(lo, hi))))
                    list(type = "data", symmetric = FALSE, array = hi - D_um,
                         arrayminus = D_um - lo, color = fg_col) else NULL,
        text = sprintf("D = %.3g µm²/s  (95%% CI %.3g to %.3g)", D_um, lo, hi),
        hovertemplate = "%{text}<extra></extra>")
    p %>% layout(
      xaxis = list(title = "D (µm²/s)", color = fg_col),
      yaxis = list(title = "", autorange = "reversed", color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col), showlegend = FALSE, uirevision = "hmm_params")
  })

  output$trk_hmm_trans_plot <- renderPlotly({
    bd   <- hmm_bundle()
    dark <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"
    K      <- bd$K
    labels <- paste0("S", seq_len(K))
    P      <- bd$P_cycle                                  # expm(Q·t_c): one cycle
    plot_ly(
      x = labels, y = rev(labels),
      z = P[rev(seq_len(K)), , drop = FALSE],
      type = "heatmap",
      colorscale = list(list(0, "#08306b"), list(0.5, "#4393c3"), list(1, "#ffffff")),
      zmin = 0, zmax = 1,
      text  = matrix(sprintf("%.3f", P[rev(seq_len(K)), ]), K, K),
      texttemplate = "%{text}",
      hovertemplate = "From %{y} → To %{x} within one cycle: %{text}<extra></extra>"
    ) %>% layout(
      xaxis = list(title = "To state",   color = fg_col, side = "top"),
      yaxis = list(title = "From state", color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col), uirevision = "hmm_trans")
  })

  output$trk_hmm_example_plot <- renderPlotly({
    bd     <- hmm_bundle()
    n_ex   <- max(1L, as.integer(isolate(input$trk_hmm_n_ex)))
    dark   <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"
    state_cols <- c("#C03030","#D08020","#2080C0","#7030A0")[seq_len(bd$K)]

    ex <- head(bd$seg_list, n_ex)
    # Lay out tracks in a grid - compute a consistent offset spacing from
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
    bd  <- hmm_bundle()
    K   <- bd$K; fit <- bd$fit
    fmt_time <- function(ms) if (!is.finite(ms)) "-"
                             else if (ms >= 1000) sprintf("%.2f s", ms / 1000)
                             else if (ms >= 1) sprintf("%.1f ms", ms)
                             else sprintf("%.2f ms", ms)
    ci_src <- if (!is.null(bd$boot)) bd$boot else fit$ci   # bootstrap preferred
    ci_lbl <- if (!is.null(bd$boot)) sprintf("bootstrap over tracks, %d refits", bd$boot$n_ok)
              else "from the likelihood curvature; too narrow when molecules differ"
    rng <- function(m, k, scale = 1, f = function(v) sprintf("%.3g", v))
      if (is.null(m)) "" else sprintf(" (%s to %s)", f(m[1, k] * scale), f(m[2, k] * scale))
    rs  <- bd$resolv
    rows <- lapply(seq_len(K), function(k) {
      flag <- c(if (isTRUE(bd$noise_state[k])) "D at the noise level",
                if (isTRUE(rs$too_short[k])) "dwell too short to resolve reliably",
                if (isTRUE(rs$few_switches[k])) "too few switches to constrain the rate",
                if (isTRUE(rs$same_D[k])) "same D as a neighbouring state (one diffusivity with non-exponential dwells)")
      tags$tr(
        tags$td(sprintf("State %d", k), style = "padding:2px 8px;"),
        tags$td(sprintf("%.0f%%", bd$fracs[k] * 100), style = "padding:2px 8px;"),
        tags$td(sprintf("%.3g", bd$D_vals[k] / 1e3), rng(ci_src$D, k, 1e-3), style = "padding:2px 8px;"),
        tags$td(fmt_time(bd$dwell_ms[k]), rng(ci_src$dwell_ms, k, 1, fmt_time),
                sprintf(" ≈ %.0f cycles", rs$dwell_cycles[k]), style = "padding:2px 8px;"),
        tags$td(if (length(flag)) tags$span(paste("⚠", paste(flag, collapse = "; ")),
                                            style = "color:#b26a00;") else "",
                style = "padding:2px 8px;"))
    })
    rate_txt <- {
      off <- character(0)
      for (a in seq_len(K)) for (b in seq_len(K)) if (a != b)
        off <- c(off, sprintf("S%d→S%d %.3g s⁻¹%s", a, b, bd$Qrate[a, b] * 1000,
                              if (!is.null(ci_src$rates_lo))
                                sprintf(" (%.3g to %.3g)", ci_src$rates_lo[a, b], ci_src$rates_hi[a, b])
                              else ""))
      paste(off, collapse = "  ·  ")
    }
    bic_txt <- {
      bo <- bd$bic_other
      if (length(bo) && any(is.finite(bo))) {
        all_bic <- c(stats::setNames(fit$BIC, K), bo)
        best <- as.integer(names(all_bic)[which.min(all_bic)])
        sprintf("BIC prefers %d states (%s).", best,
                paste(sprintf("K=%s: %+.0f", names(all_bic), all_bic - min(all_bic, na.rm = TRUE)),
                      collapse = ", "))
      } else ""
    }
    # JDD cross-check: at short lags the JDD component weights estimate the same
    # occupancies (per step rather than per unit time), with no segmentation at all.
    jdd_txt <- {
      jb <- tryCatch(jdd_bundle(), error = function(e) NULL)
      if (!is.null(jb) && isTRUE(jb$n_comp == K) && length(jb$ws) == K) {
        o <- order(jb$D_vals)
        sprintf("JDD (%d components): D %s µm²/s, weights %s.", K,
                paste(sprintf("%.3g", jb$D_vals[o] / 1e3), collapse = " / "),
                paste(sprintf("%.0f%%", 100 * jb$ws[o]), collapse = " / "))
      } else NULL
    }
    tags$div(style = "font-size:0.88em; margin:4px 0;",
      tags$p(style = "color:#888; margin:2px 0;",
        sprintf("n = %d tracks%s, %d steps  |  cycle %.3g ms  |  σ_loc %s (median %.2f nm)  |  log-likelihood %.1f",
                bd$n_tracks,
                if (is.finite(bd$n_tracks_all %||% NA)) sprintf(" of %d, drawn at random to meet the step budget",
                                                                bd$n_tracks_all) else "",
                bd$n_steps, bd$t_cycle_ms, bd$sigma_mode, bd$sigma_loc_nm, fit$loglik)),
      tags$table(style = "border-collapse:collapse;",
        tags$thead(tags$tr(lapply(
          c("", "Occupancy", "D (µm²/s, 95% CI)", "Mean dwell (95% CI)", ""),
          function(h) tags$th(h, style = "padding:2px 8px; border-bottom:1px solid #ccc;")))),
        tags$tbody(rows)),
      tags$p(style = "color:#444; margin:6px 0 2px;", tags$b("Switching rates: "), rate_txt),
      tags$p(style = "color:#888; font-size:0.9em; margin:2px 0;",
        sprintf("CIs %s. %s", ci_lbl, bic_txt)),
      if (!is.null(jdd_txt))
        tags$p(style = "color:#888; font-size:0.9em; margin:2px 0;", tags$b("Cross-check. "), jdd_txt),
      if (!isTRUE(fit$converged))
        tags$p(style = "color:#b26a00; margin:2px 0;", "The optimiser did not report convergence; treat the estimates with caution."))
  })

  output$trk_dl_hmm_csv <- downloadHandler(
    filename = function() paste0("hmm_states_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      bd  <- hmm_bundle()
      out <- do.call(rbind, bd$seg_list)
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # seg_msd_rv() - PER-STATE MSD (one MSD curve per HMM diffusive state)
  # ==========================================================================
  #
  # GOAL
  #   Fix the biggest weakness of a whole-track MSD: AVERAGING OVER PHASES.
  #   If a molecule spends part of its track confined at a wall and part of it
  #   being transported, one MSD curve mixes both and describes neither. The
  #   fitted "D" is then an artefact of the mixing ratio, not a property of the
  #   molecule.
  #
  #   This splits every track by the HMM's state labels FIRST, then builds a
  #   separate MSD for each state. Now the physics-specific models finally mean
  #   something:
  #       Confined model on the slow/wall state  -> a real confinement size L
  #       Directed model on the fast state       -> a real velocity v
  #
  # THE MATH
  #   Identical to the main MSD engine - same cumulative cycle-index binning,
  #   same fit window, same model set - but pairs are pooled PER STATE instead of
  #   per dataset. Only pairs where BOTH localizations carry the same state label
  #   contribute, so a curve never mixes phases.
  #
  #   Being an ensemble over many short same-state runs (rather than a few long
  #   tracks) is actually an advantage here: it is the vbSPT-style trick of
  #   pooling statistical power across many short segments.
  #
  # OBJECTS IT CALLS (inputs)
  #   hmm_bundle()$seg_list ... per-localization state labels - THE dependency;
  #                             this analysis cannot run until the HMM has run
  #   trk_cycle_time(), trk_cluster_cycle()$sw
  #   .msd_model_fit() ........ the shared fitter, so the model/window/linear-τ/
  #                             iterative settings are INHERITED from the main
  #                             MSD panel rather than duplicated here
  #   inputs: trk_seg_msd_run (Run trigger), trk_seg_msd_min_pts,
  #           trk_msd_ndim, trk_uni_max_lag (+ the inherited fit settings)
  #
  # STEP BY STEP
  #   1. Take the HMM's Viterbi labels; states are already ordered slow -> fast.
  #   2. Within each track, run-length-split into CONTIGUOUS same-state runs and
  #      discard runs shorter than `min_pts` (too few pairs to be meaningful).
  #   3. For each run, accumulate displacement pairs into that STATE's bins,
  #      using the same cumulative cycle index as the main MSD.
  #   4. Per state: ensemble MSD = Σdr² / Σn per lag.
  #   5. Fit each state's curve with the inherited model and window.
  #   6. Cross-check: the fitted D per state should broadly agree with the HMM's
  #      own D_k = σ_k²/2. Large disagreement means the state assignment is
  #      unreliable - treat both with suspicion.
  #
  # RETURNS a list with, per state: the MSD curve (lag, tau, msd, npairs), the
  # fit result and its model-specific extras (L for confined, v for directed,
  # α for anomalous), segment and point counts.
  #
  # CAVEAT: everything here is downstream of the HMM. A poor state decode makes
  # these curves meaningless no matter how good the fit looks. Confined/anomalous
  # MSD on MINFLUX is also flagged exploratory (see the MSD model note).
  # ==========================================================================
  SEG_MSD_PAL <- c("#3060C0","#C03030","#2CA02C","#FF7F0E","#9467BD","#8C564B")

  seg_msd_rv <- eventReactive(input$trk_seg_msd_run, {
    req(isTRUE(input$trk_show_seg_msd))
    bd <- tryCatch(hmm_bundle(), error = function(e) NULL)
    shiny::validate(shiny::need(
      !is.null(bd) && !is.null(bd$seg_list) && length(bd$seg_list) > 0,
      "Enable the Diffusion State HMM (above) first - its Viterbi states define the segments."))
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
               # Same motion/projection correction as the main MSD panel, so a
               # per-state D is on the same footing as the ensemble D.
               tryCatch(.msd_model_fit(tau_use[okf], msd[okf], npk[okf], ndim, model, use_iter, t_c,
                                       coef_override = 2 * .msd_d_eff(ndim), R = trk_blur()$R_msd),
                        error = function(e) NULL) else NULL
      list(state = k, n_seg = n_seg[k], n_pts = n_pts[k], D_hmm = bd$D_vals[k] / 1e3,
           lag = lag_v, tau = tau_use, msd = msd, npairs = npk, fr = fr)
    }))
    shiny::validate(shiny::need(length(states) > 0, "No segments long enough - lower 'Min points per segment' or the HMM min points."))
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
                      else sprintf("State %d confined: D=%.3g, L_c=%.0f nm", st$state, fr$D / 1e3, fr$L_conf),
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
      param <- if (is.null(fr)) "-" else switch(fr$model %||% "linear",
        confined  = if (isTRUE(fr$degenerate)) "no plateau (≈linear)" else sprintf("L_c = %.0f nm", fr$L_conf),
        directed  = sprintf("v = %.3g µm/s", fr$v_um_s),
        anomalous = sprintf("α = %.2f", fr$alpha), "-")
      D_msd <- if (!is.null(fr) && is.finite(fr$D)) sprintf("%.3g", fr$D / 1e3) else "-"
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
        "States are the HMM diffusion states (sorted slow→fast). Switch the MSD ",
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

  # ==========================================================================
  # dmap_bundle() - SPATIALLY RESOLVED DIFFUSION MAP
  # ==========================================================================
  #
  # GOAL
  #   Every other estimator here collapses the data to ONE number (or a few
  #   states) for the whole field of view. This asks instead: does D vary with
  #   LOCATION? It bins the field into cells and estimates a local D in each,
  #   revealing slow zones (membranes, walls, condensates) vs fast interior.
  #
  # THE MATH - a SINGLE-LAG estimator, not a fitted MSD curve
  #   Rearranging <r²> = 2dDΔt for one step gives the local diffusivity directly:
  #
  #        D_local = r² / (2·d·Δt)        (2d = 4 in 2D, 6 in 3D)
  #
  #   Each step is assigned to the cell containing its MIDPOINT, and the cell's D
  #   is the MEAN of its steps' D, with the noise floor subtracted once per cell:
  #
  #        D_cell = max(0, mean(r²/(4Δt)) − σ²·mean(1/Δt)) / (1 − 2R·t_c·mean(1/Δt))
  #
  #   (the denominator is the motion-blur correction added 2026-09-17; without it
  #   the map read −30% on blurred simulated tracks. Estimator: est_dmap.R)
  #
  #   WHY THE MEAN AND NOT THE MEDIAN. The relation <r²> = 4DΔt + 4σ² holds for
  #   the MEAN only. A single 2D step's r² is exponentially distributed, so its
  #   median sits at ln2 ≈ 0.69 of the mean. The old median form, with the noise
  #   floor then subtracted per step and clamped at 0, measured (σ 5 nm, Δt 0.2 ms)
  #   −50% at D = 0.2 µm²/s and read exactly 0 at D ≤ 0.05. The mean form is
  #   unbiased at every D, whatever the step distribution. The price is
  #   robustness: one mis-linked jump inflates its cell, so the per-cell standard
  #   error is returned alongside D and upstream track cleaning carries the load.
  #   Subtracting per cell (not per step) avoids clamping each step at 0, which
  #   would rectify the noise and bias slow cells upward.
  #
  #   WHY SINGLE-LAG RATHER THAN A PER-CELL MSD FIT: a cell contains few steps.
  #   Fitting a multi-parameter MSD curve per cell would be badly under-
  #   constrained, whereas the single-lag form needs only first moments. It is
  #   the short-lag limit of the covariance-based estimator (Vestergaard,
  #   Blainey & Flyvbjerg 2014, Phys. Rev. E 89:022726); MINFLUX-specific
  #   validity in Vogler et al. 2025 (Commun. Biol. 8:1573).
  #   [Previously mis-cited to Ostersehlt 2022 - a DNA-PAINT imaging paper with
  #   no diffusion mapping. Corrected.]
  #
  #   NOISE FLOOR. An immobile molecule still yields apparent D = σ²/Δt (header
  #   §2b eq. 3), so without correction every cell reads as diffusing. Subtracting
  #   that offset is what lets genuinely immobile regions show as ~0.
  #
  # THE MINFLUX ADAPTATION
  #   "Integrated steps only" drops steps spanning dark cycles. Two reasons: a
  #   long-gap step has a large Δt whose midpoint is not a meaningful location
  #   (the molecule wandered between the two ends), and it would smear the map.
  #
  # OBJECTS IT CALLS (inputs)
  #   trk(), trk_cycle_time(), trk_sigma_loc_eff()
  #   inputs: trk_dmap_run (Run trigger), trk_dmap_bins, trk_dmap_min_pts,
  #           trk_dmap_integrated_only, trk_dmap_sub_noise
  #
  # STEP BY STEP
  #   1. Build per-step r², Δt and midpoint (x,y) within each track.
  #   2. Optionally keep only single-cycle (integrated) steps.
  #   3. Bin midpoints onto a regular grid.
  #   4. Per cell: D = [mean(r²/(4Δt)) − σ²·mean(1/Δt)]₊ / (1 − 2R·t_c·mean(1/Δt)), plus its SE.
  #   5. Blank cells with fewer than `min_pts` steps - too few to trust.
  #
  # CAVEAT: single-lag D is an INSTANTANEOUS estimate; it says nothing about
  # confinement or anomaly. And where sampling is spatially uneven, the map's
  # blank/noisy regions reflect coverage, not biology.
  # ==========================================================================
  dmap_bundle <- eventReactive(input$trk_dmap_run, {
    req(isTRUE(input$trk_show_dmap))
    df <- trk(); req(nrow(df) > 0)
    t_c_dmap <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
    # Motion blur needs the cycle time. Without it the fit would silently run with no
    # blur correction, which reads D about a third low at MINFLUX cycle times, so the
    # panel asks for the calibration instead of guessing.
    shiny::validate(shiny::need(is.finite(t_c_dmap) && t_c_dmap > 0,
      "Cycle time unknown - run Data Quality & Calibration first (the blur correction needs it)."))
    # σ_loc floor is optional (toggle); the app-wide value is used when set.
    sig_loc  <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    sig_use  <- if (isTRUE(input$trk_dmap_sub_noise) && is.finite(sig_loc) && sig_loc > 0) sig_loc else NA_real_
    # Estimator lives in est_dmap.R (pure, tested); blur offset as in the MSD fits.
    withProgress(message = "Computing D-map…", value = 0.1, {
      res <- dmap_compute(df, n_bins = input$trk_dmap_bins, min_pts = input$trk_dmap_min_pts,
                          t_c = t_c_dmap, integrated_only = isTRUE(input$trk_dmap_integrated_only),
                          sigma_loc = sig_use,
                          R_abs = trk_blur()$R_abs)
    })
    shiny::validate(shiny::need(is.null(res$error), res$error %||% ""))
    res
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
      "No cells met the minimum step threshold - lower 'Min steps per cell' or reduce grid resolution."))

    # Clip to 2nd–98th percentile so one extreme cell doesn't dominate the scale
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
        d  <- bd$D_mat[c, r]
        se <- bd$SE_mat[c, r]
        n  <- bd$N_mat[c, r]
        hover[r, c] <- if (is.finite(d))
          sprintf("x: %.0f nm<br>y: %.0f nm<br>D: %.3g \u00B1 %.2g nm\u00B2/ms (SE)<br>n steps: %d",
                  bd$x_centers[c], bd$y_centers[r], d, se, n)
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
                     D_se_nm2_per_ms = bd$SE_mat[ix, iy],
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
      loop$ms_per_unit   <- 0   # sequential mode uses a fixed frame interval scaled by the Playback speed slider
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
  # view_mode / z_scale observers compare against loop$rendered_* (set by the render)
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
      # "refuse to stop" - it starved the Stop click + reactive flush).
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
    # Consume programmatic echoes FIRST, before the is_playing() gate. draw_frame bumps
    # slider_guard for every updateSliderInput it pushes; if the is_playing() gate ran
    # first those credits were never spent during playback, so the counter accumulated
    # and was then wiped by slider_guard(0L) at stop/end - after which late-arriving
    # echoes from the last frames sailed through and re-drove draw_frame.
    guard <- isolate(slider_guard())
    if (guard > 0L) {
      slider_guard(guard - 1L)
      return(invisible(NULL))
    }
    if (isTRUE(is_playing())) return(invisible(NULL))

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
      # Skip entirely in full-tracks mode - the static plot needs no proxy update.
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
            # (trails are set at renderPlotly time; conf coloring needs a proxy push).
            # Skipped in per-point colour modes (dark cycles / speed): there draw_frame
            # above already coloured each trail from its leading point, and pushing the
            # per-TID colours here would clobber it back to one colour per track.
            if (identical(tm$trails_mode, "per_tid") && tm$n_tid > 0 &&
                !isTRUE(pre$color_by_dark)) {
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
    # Keep the 2D axis titles matching the active viewing plane so a dark-mode /
    # grid relayout does not silently revert an XZ/YZ side view to "X/Y".
    plane2 <- if (!is3d) (input$trk_view_plane %||% "XY") else "XY"
    plab   <- switch(plane2,
      XZ = list(x = "X (nm)", y = "Z (nm)"),
      YZ = list(x = "Y (nm)", y = "Z (nm)"),
      list(x = "X (nm)", y = "Y (nm)"))

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
      ovl <- if (plane2 == "XY") tryCatch(isolate(trk_image_overlay()), error = function(e) NULL) else NULL
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
            title = if (hide) "" else plab$x,
            showgrid = !hide, zeroline = !hide, showticklabels = !hide,
            gridcolor = gc, color = fg
          ), x_range_kv),
          yaxis = c(list(
            title = if (hide) "" else plab$y,
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
      trk_axis_x_span(x1 - x0)
      y0 <- suppressWarnings(as.numeric(ed[["yaxis.range[0]"]]))
      y1 <- suppressWarnings(as.numeric(ed[["yaxis.range[1]"]]))
      if (length(y0) == 1 && length(y1) == 1 && is.finite(y0) && is.finite(y1) && y1 > y0)
        trk_zoom_ranges(list(x = c(x0, x1), y = c(y0, y1)))
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
    # The spatial crop filters the real x,y. In an XZ/YZ side view the plotted
    # vertical axis is Z, so the zoom box's "y" is a z-range - cropping on it would
    # silently filter the wrong axis. Require the XY plane for cropping.
    if (!identical(input$trk_view_plane %||% "XY", "XY")) {
      showNotification("Crop to view only works in the XY plane. Switch the viewing plane to XY first.",
                       type = "warning", duration = 6)
      return()
    }
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


  # ============================================================================
  # INSTRUMENT CALIBRATION - a and D_w from an immobile sample, kept as profiles
  # ============================================================================
  # GOAL: measure the photon-budget prefactor a and the wander floor D_w for one
  # acquisition sequence (calib_cell, R/est_calib.R), keep the result as a named
  # profile that records the sequence it belongs to, and let profiles be saved,
  # loaded and applied, globally or per file in Compare Cells. A profile is the
  # unit a lab measures once per sequence and reuses; without one, a falls back to
  # the ideal-estimator value of the pattern and the readout says "uncalibrated".
  trk_cal_profiles <- reactiveVal(list())   # name -> profile (list)
  trk_cal_result   <- reactiveVal(NULL)     # per-cell table of the last measurement

  # Loaded files to choose the immobile ones from.
  observe({
    fs  <- names(trk_store())
    sel <- isolate(input$trk_cal_files)
    updateSelectInput(session, "trk_cal_files", choices = fs, selected = intersect(sel, fs))
  })

  # A loaded sequence file fills in the descriptors the profile records.
  observeEvent(input$trk_seq_json, {
    sp <- tryCatch(trk_seq_params(), error = function(e) NULL)
    if (is.null(sp) || !isTRUE(sp$ok)) return(invisible(NULL))
    if (is.finite(sp$patDwellTime_ms %||% NA)) updateNumericInput(session, "trk_cal_dwell", value = round(sp$patDwellTime_ms * 1000))
    if (is.finite(sp$phtLimit %||% NA))        updateNumericInput(session, "trk_cal_pl", value = sp$phtLimit)
    if (is.finite(sp$damping %||% NA))         updateNumericInput(session, "trk_cal_damping", value = sp$damping)
  }, ignoreInit = TRUE, priority = -1)

  observeEvent(input$trk_cal_run, {
    files <- input$trk_cal_files
    if (!length(files)) {
      showNotification("Choose the immobile-sample files first.", type = "warning", duration = 5)
      return(invisible(NULL))
    }
    L   <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
    if (!is.finite(L) || L <= 0) L <- DMAX_L_TCP_DEFAULT
    pat  <- .tcp_pattern_name(input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT)
    if (is.na(pat) || !nzchar(pat)) pat <- TCP_PATTERN_DEFAULT
    filt <- trk_filter_params()
    rows <- withProgress(message = "Measuring each immobile cell...", value = 0, {
      lapply(files, function(f) {
        incProgress(1 / length(files), detail = f)
        d <- .cmp_read_file(f)
        if (is.null(d)) d <- data.frame(tid = character(0), x = numeric(0), y = numeric(0),
                                        tim = numeric(0), eco = numeric(0))
        tryCatch(calib_cell(d, f, L, filt), error = function(e) {
          r <- calib_cell(d[0, , drop = FALSE], f, L); r$note <- conditionMessage(e); r })
      })
    })
    tab <- do.call(rbind, rows)
    attr(tab, "L") <- L; attr(tab, "pattern") <- pat
    trk_cal_result(tab)
  })

  output$trk_cal_summary <- renderUI({
    tab <- trk_cal_result(); req(!is.null(tab))
    sm  <- calib_summary(tab)
    AMB <- "#B8860B"
    ideal <- .tcp_geom(attr(tab, "pattern"))$a_ideal
    tagList(
      tags$p(style = "margin:4px 0;",
        tags$b(sprintf("a = %.2f", sm$a)),
        sprintf("  (quartiles %.2f–%.2f, %d cell%s)  ·  ", sm$a_q25, sm$a_q75, sm$n_cells,
                if (sm$n_cells == 1L) "" else "s"),
        tags$b(sprintf("D_w = %.0f nm²/ms", sm$D_w)),
        sprintf("  (%.0f–%.0f)", sm$Dw_q25, sm$Dw_q75)),
      helpText(style = "font-size:0.82em; color:#666; margin-top:0;",
        sprintf("L %.1f nm, %s. The ideal-estimator value for this pattern is %.2f.",
                attr(tab, "L"), attr(tab, "pattern"), ideal)),
      if (sm$n_cells < 3L) helpText(style = sprintf("color:%s; font-size:0.82em;", AMB),
        "Fewer than 3 cells: the quartiles say little about how well a is known."),
      if (any(nzchar(tab$note))) helpText(style = sprintf("color:%s; font-size:0.82em;", AMB),
        sprintf("%d file(s) gave no value; see the note column.", sum(nzchar(tab$note)))))
  })

  output$trk_cal_results <- DT::renderDT({
    tab <- trk_cal_result(); req(!is.null(tab))
    show <- tab[, c("file", "n_tracks", "t_c_ms", "eco_med", "msd1", "sigma_nm", "a", "D_w", "note")]
    num  <- vapply(show, is.numeric, logical(1))
    show[num] <- lapply(show[num], function(v) signif(v, 3))
    DT::datatable(show, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE, scrollX = TRUE))
  })

  observeEvent(input$trk_cal_save, {
    tab <- trk_cal_result()
    nm  <- trimws(input$trk_cal_name %||% "")
    if (is.null(tab)) { showNotification("Measure first.", type = "warning", duration = 5); return(invisible(NULL)) }
    if (!nzchar(nm))  { showNotification("Give the profile a name.", type = "warning", duration = 5); return(invisible(NULL)) }
    sm <- calib_summary(tab)
    if (!sm$n_cells) { showNotification("No cell gave a value of a; nothing to save.", type = "error", duration = 6); return(invisible(NULL)) }
    num <- function(x) { v <- suppressWarnings(as.numeric(x)); if (length(v) && is.finite(v[1])) v[1] else NA_real_ }
    profs <- trk_cal_profiles()
    profs[[nm]] <- list(name = nm, a = sm$a, a_q25 = sm$a_q25, a_q75 = sm$a_q75,
                        D_w = sm$D_w, Dw_q25 = sm$Dw_q25, Dw_q75 = sm$Dw_q75, n_cells = sm$n_cells,
                        pattern = attr(tab, "pattern"), L_nm = attr(tab, "L"),
                        dwell_us = num(input$trk_cal_dwell), photon_limit = num(input$trk_cal_pl),
                        damping = num(input$trk_cal_damping), note = trimws(input$trk_cal_note %||% ""),
                        date = format(Sys.Date()), files = tab$file[is.finite(tab$a)])
    trk_cal_profiles(profs)
    showNotification(sprintf("Profile '%s' saved (a = %.2f, D_w = %.0f nm²/ms).", nm, sm$a, sm$D_w),
                     type = "message", duration = 6)
  })

  # Keep the Acquisition Parameters selector in step with the profiles.
  observeEvent(trk_cal_profiles(), {
    nms <- names(trk_cal_profiles())
    sel <- isolate(input$trk_cal_profile %||% "")
    updateSelectInput(session, "trk_cal_profile",
                      choices = c("None (uncalibrated)" = "", stats::setNames(nms, nms)),
                      selected = if (sel %in% nms) sel else "")
  }, ignoreInit = TRUE)

  output$trk_cal_profiles_tab <- DT::renderDT({
    profs <- trk_cal_profiles(); req(length(profs) > 0)
    f <- function(x, d = 2) if (is.null(x) || !is.finite(x)) "" else formatC(x, digits = d, format = "f")
    tab <- do.call(rbind, lapply(profs, function(p) data.frame(
      profile = p$name, a = sprintf("%s (%s–%s)", f(p$a), f(p$a_q25), f(p$a_q75)),
      D_w = f(p$D_w, 0), cells = p$n_cells %||% NA, pattern = p$pattern %||% "", L_nm = f(p$L_nm, 1),
      dwell_us = f(p$dwell_us, 0), photon_limit = f(p$photon_limit, 0), damping = f(p$damping, 0),
      note = p$note %||% "", date = p$date %||% "", stringsAsFactors = FALSE)))
    DT::datatable(tab, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE, scrollX = TRUE))
  })

  output$trk_cal_profiles_dl <- downloadHandler(
    filename = function() paste0("calibration_profiles_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"),
    content  = function(file) writeLines(jsonlite::toJSON(unname(trk_cal_profiles()), pretty = TRUE,
                                                          auto_unbox = TRUE, na = "null", digits = NA), file))

  observeEvent(input$trk_cal_profiles_ul, {
    up <- tryCatch(jsonlite::fromJSON(input$trk_cal_profiles_ul$datapath, simplifyVector = FALSE),
                   error = function(e) NULL)
    if (!is.null(up) && !is.null(up$name)) up <- list(up)          # a single profile
    ok <- Filter(function(p) is.list(p) && is.character(p$name) && nzchar(p$name) &&
                   is.numeric(p$a) && is.finite(p$a) && p$a > 0, up %||% list())
    if (!length(ok)) {
      showNotification("No valid calibration profile in that file.", type = "error", duration = 6)
      return(invisible(NULL))
    }
    num <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x))
    profs <- trk_cal_profiles()
    for (p in ok) {
      for (k in c("a", "a_q25", "a_q75", "D_w", "Dw_q25", "Dw_q75", "L_nm", "dwell_us", "photon_limit", "damping"))
        p[[k]] <- num(p[[k]])
      p$files <- unlist(p$files)
      profs[[p$name]] <- p
    }
    trk_cal_profiles(profs)
    showNotification(sprintf("%d profile(s) loaded.", length(ok)), type = "message", duration = 5)
  })

  # Applying a profile sets a, D_w, the pattern and L; "None" returns to uncalibrated.
  observeEvent(input$trk_cal_profile, {
    nm <- input$trk_cal_profile %||% ""
    if (!nzchar(nm)) {
      updateNumericInput(session, "trk_crlb_pref", value = NA)
      updateNumericInput(session, "trk_dw_floor", value = 0)
      return(invisible(NULL))
    }
    p <- trk_cal_profiles()[[nm]]; if (is.null(p)) return(invisible(NULL))
    updateNumericInput(session, "trk_crlb_pref", value = round(p$a, 3))
    if (is.finite(p$D_w %||% NA)) updateNumericInput(session, "trk_dw_floor", value = round(max(0, p$D_w)))
    pat <- .tcp_pattern_name(p$pattern %||% "")
    if (!is.na(pat) && nzchar(pat)) updateSelectInput(session, "trk_tcp_pattern", selected = pat)
    if (is.finite(p$L_nm %||% NA)) updateNumericInput(session, "trk_L_tcp", value = p$L_nm)
  }, ignoreInit = TRUE)

  # ============================================================================
  # COMPARE CELLS - one value per file (= one cell), conditions compared across cells
  # ============================================================================
  # GOAL: label every loaded file with its condition and the settings of the
  # sequence it was acquired with, analyse each file on its own (percell_summary,
  # R/est_percell.R: own cycle time, own gap split, own L / pattern / a / D_w), and
  # compare conditions with the cell as the unit of replication. Pooling tracks
  # across cells would treat correlated tracks as independent and let one busy
  # cell dominate a condition.

  # Every loaded file, whichever view is active.
  cmp_file_names <- reactive(names(trk_store()))

  # One stored file, prepared the way trk_raw_loaded() / trk_raw_df() prepare the
  # active view: character track ids, times in ms, the axial factor applied.
  .cmp_read_file <- function(name) {
    st <- isolate(trk_store())[[name]]
    d  <- if (is.null(st)) NULL else read_chunk_file(st$path)
    if (is.null(d) || !nrow(d)) return(NULL)
    d$tid <- as.character(d$tid)
    for (cc in intersect(c("x", "y", "tim"), names(d))) d[[cc]] <- as.numeric(d[[cc]])
    d$tim <- d$tim * 1000                                   # s -> ms
    zf <- isolate(trk_z_applied())
    if (is.finite(zf) && zf > 0 && zf != 1 && "z" %in% names(d) && any(is.finite(d$z))) d$z <- d$z * zf
    d
  }

  cmp_tab     <- reactiveVal(NULL)   # per-file settings: file, condition, L_nm, pattern, a, D_w
  cmp_results <- reactiveVal(NULL)   # per-cell numbers from the last run, with the settings used

  # New upload: one row per file. Sequence settings start BLANK, and a blank means
  # "the current Acquisition Parameters value", so a dataset acquired with one
  # sequence follows the sequence file whenever it is loaded, and only settings
  # that differ between files need entering. Adding or removing a file keeps the
  # other rows.
  observeEvent(cmp_file_names(), {
    nm  <- cmp_file_names()
    old <- isolate(cmp_tab())
    if (!length(nm)) { cmp_tab(NULL); cmp_results(NULL); return(invisible(NULL)) }
    tb <- data.frame(file = nm, condition = "", calibration = "", L_nm = NA_real_, pattern = "",
                     a = NA_real_, D_w = NA_real_, stringsAsFactors = FALSE)
    if (!is.null(old)) {                                    # adding/removing keeps other rows
      m <- match(tb$file, old$file); k <- which(!is.na(m))
      if (length(k)) tb[k, -1] <- old[m[k], -1]
    }
    cmp_tab(tb)
    r <- isolate(cmp_results())                             # drop results of removed files
    if (!is.null(r)) cmp_results(if (any(r$file %in% nm)) r[r$file %in% nm, , drop = FALSE] else NULL)
  }, ignoreNULL = FALSE)

  # The table with every blank filled from the current Acquisition Parameters:
  # the values a run actually uses.
  cmp_resolved <- reactive({
    tb <- cmp_tab(); if (is.null(tb)) return(NULL)
    num <- function(x, d) { v <- suppressWarnings(as.numeric(x)); if (length(v) && is.finite(v[1])) v[1] else d }
    pat <- tolower(input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT)
    if (!pat %in% names(TCP_PATTERN_GEOM)) pat <- TCP_PATTERN_DEFAULT
    # A named calibration profile fills its row's blanks first (it belongs to one
    # sequence), then the Acquisition Parameters fill what is still blank.
    profs <- trk_cal_profiles()
    for (i in seq_len(nrow(tb))) {
      pr <- profs[[tb$calibration[i]]]
      if (is.null(pr)) next
      if (!is.finite(tb$L_nm[i]) && is.finite(pr$L_nm %||% NA)) tb$L_nm[i] <- pr$L_nm
      if (!nzchar(tb$pattern[i])) { pp <- .tcp_pattern_name(pr$pattern %||% ""); if (!is.na(pp)) tb$pattern[i] <- pp }
      if (!is.finite(tb$a[i])   && is.finite(pr$a %||% NA))   tb$a[i]   <- pr$a
      if (!is.finite(tb$D_w[i]) && is.finite(pr$D_w %||% NA)) tb$D_w[i] <- max(0, pr$D_w)
    }
    tb$L_nm[!is.finite(tb$L_nm)] <- num(input$trk_L_tcp, DMAX_L_TCP_DEFAULT)
    tb$pattern[!nzchar(tb$pattern)] <- pat
    a_in <- suppressWarnings(as.numeric(input$trk_crlb_pref))
    for (i in which(!is.finite(tb$a)))                       # entered a, else the pattern's ideal
      tb$a[i] <- if (length(a_in) && is.finite(a_in[1]) && a_in[1] > 0) a_in[1]
                 else .tcp_geom(tb$pattern[i])$a_ideal
    tb$D_w[!is.finite(tb$D_w)]   <- num(input$trk_dw_floor, 0)
    tb
  })

  # One value typed for its column, or NULL when it does not fit (the caller then
  # leaves the table unchanged and says so).
  # An empty value clears a setting back to blank (= Acquisition Parameters).
  .cmp_coerce <- function(col, val) {
    if (is.null(val) || (length(val) == 1L && is.na(val))) val <- ""
    if (col %in% c("L_nm", "a", "D_w")) {
      if (!nzchar(trimws(as.character(val)))) return(NA_real_)
      v <- suppressWarnings(as.numeric(val))
      if (!length(v) || !is.finite(v) || v < 0 || (col != "D_w" && v == 0)) return(NULL)
      v
    } else if (col == "pattern") {
      v <- .tcp_pattern_name(val)                      # names, short forms, or 3 / 4 / 6
      if (is.na(v)) return(NULL)
      v
    } else trimws(as.character(val))
  }

  output$trk_cmp_table <- DT::renderDT({
    tb <- cmp_tab(); req(!is.null(tb))
    DT::datatable(tb, rownames = FALSE, selection = "none",
                  editable = list(target = "cell", disable = list(columns = 0)),
                  options = list(dom = "t", paging = FALSE, scrollY = "260px",
                                 scrollCollapse = TRUE, ordering = FALSE))
  })

  # Cell edits from the table. DT reports 0-based columns without row names.
  observeEvent(input$trk_cmp_table_cell_edit, {
    tb <- cmp_tab(); info <- input$trk_cmp_table_cell_edit
    req(!is.null(tb), NROW(info) > 0)
    for (k in seq_len(nrow(info))) {
      col <- names(tb)[info$col[k] + 1L]
      v   <- .cmp_coerce(col, info$value[k])
      if (is.null(v)) {
        showNotification(sprintf("'%s' is not a valid value for %s.", info$value[k], col),
                         type = "warning", duration = 5)
        next
      }
      tb[info$row[k], col] <- v
    }
    cmp_tab(tb)
  })

  # Bulk fill: one column, every file whose name contains the text.
  observeEvent(input$trk_cmp_fill_apply, {
    tb <- cmp_tab(); req(!is.null(tb))
    col <- input$trk_cmp_fill_col %||% "condition"
    v   <- .cmp_coerce(col, input$trk_cmp_fill_value %||% "")
    if (is.null(v)) {
      ok <- if (col == "pattern") paste0(paste(names(TCP_PATTERN_GEOM), collapse = ", "),
                                         ", or 3 / 4 / 6 vertices") else "a positive number"
      showNotification(sprintf("%s needs %s.", col, ok), type = "warning", duration = 5)
      return(invisible(NULL))
    }
    pat <- tolower(trimws(input$trk_cmp_fill_match %||% ""))
    hit <- if (nzchar(pat)) grepl(pat, tolower(tb$file), fixed = TRUE) else rep(TRUE, nrow(tb))
    if (!any(hit)) {
      showNotification("No file name contains that text.", type = "warning", duration = 5)
      return(invisible(NULL))
    }
    tb[hit, col] <- v
    cmp_tab(tb)
    showNotification(sprintf("%s set for %d file(s).", col, sum(hit)), type = "message", duration = 4)
  })

  # Save / load the settings table, so a dataset's conditions are entered once and
  # travel with the analysis.
  output$trk_cmp_table_dl <- downloadHandler(
    filename = function() paste0("cell_conditions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) utils::write.csv(cmp_tab(), file, row.names = FALSE, na = ""))

  observeEvent(input$trk_cmp_table_ul, {
    tb <- cmp_tab(); req(!is.null(tb))
    up <- tryCatch(utils::read.csv(input$trk_cmp_table_ul$datapath, stringsAsFactors = FALSE,
                                   check.names = FALSE, colClasses = "character"),
                   error = function(e) NULL)
    if (is.null(up) || !"file" %in% names(up)) {
      showNotification("The table needs a 'file' column.", type = "error", duration = 6)
      return(invisible(NULL))
    }
    m <- match(tb$file, up$file); n_bad <- 0L
    for (col in intersect(c("condition", "calibration", "L_nm", "pattern", "a", "D_w"), names(up)))
      for (i in which(!is.na(m))) {
        v <- .cmp_coerce(col, up[[col]][m[i]])
        if (is.null(v)) n_bad <- n_bad + 1L else tb[i, col] <- v
      }
    cmp_tab(tb)
    showNotification(sprintf("Settings loaded for %d of %d file(s)%s.", sum(!is.na(m)), nrow(tb),
                             if (n_bad) sprintf("; %d invalid value(s) skipped", n_bad) else ""),
                     type = if (any(!is.na(m))) "message" else "warning", duration = 6)
  })

  # Run: analyse every loaded file on its own, straight from the store, so the
  # comparison does not depend on which view is active.
  observeEvent(input$trk_cmp_run, {
    tb <- cmp_resolved(); req(!is.null(tb), nrow(tb) > 0)
    filt <- trk_filter_params()                             # the Display tab's filters
    nu <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
    if (!is.finite(nu) || nu < 0) nu <- 0
    rows <- withProgress(message = "Analysing each cell...", value = 0, {
      lapply(seq_len(nrow(tb)), function(i) {
        incProgress(1 / nrow(tb), detail = tb$file[i])
        di <- .cmp_read_file(tb$file[i])
        if (is.null(di)) di <- data.frame(tid = character(0), x = numeric(0), y = numeric(0), tim = numeric(0))
        r  <- tryCatch(percell_summary(di, tb$file[i], L = tb$L_nm[i], pattern = tb$pattern[i],
                                       a = tb$a[i], D_w = tb$D_w[i], nu = nu, filters = filt),
                       error = function(e) {
                         x <- percell_summary(di[0, , drop = FALSE], tb$file[i], tb$L_nm[i])
                         x$note <- conditionMessage(e); x })
        if (nzchar(tb$calibration[i]) && is.null(trk_cal_profiles()[[tb$calibration[i]]]))
          r$note <- trimws(paste(r$note, sprintf("calibration '%s' is not loaded", tb$calibration[i])))
        cbind(r, calibration = tb$calibration[i], L_nm = tb$L_nm[i], pattern = tb$pattern[i],
              a = tb$a[i], D_w = tb$D_w[i], stringsAsFactors = FALSE)
      })
    })
    cmp_results(do.call(rbind, rows))
  })

  # Results joined to the CURRENT conditions, so relabelling needs no re-run; flags
  # results whose sequence settings were edited after the run.
  cmp_view <- reactive({
    r <- cmp_results(); tb <- cmp_resolved()
    if (is.null(r) || is.null(tb)) return(NULL)
    m <- match(r$file, tb$file)
    r$condition <- ifelse(is.na(m), "", tb$condition[m])
    stale <- any(!is.na(m) & (abs(r$L_nm - tb$L_nm[m]) > 1e-9 | r$pattern != tb$pattern[m] |
                                abs(r$a - tb$a[m]) > 1e-9 | abs(r$D_w - tb$D_w[m]) > 1e-9))
    attr(r, "stale") <- isTRUE(stale)
    r
  })

  output$trk_cmp_plot <- renderPlotly({
    r <- cmp_view(); req(!is.null(r))
    met <- input$trk_cmp_metric %||% "D_um2_s"
    d   <- r[nzchar(trimws(r$condition)) & is.finite(r[[met]]), , drop = FALSE]
    shiny::validate(shiny::need(nrow(d) > 0,
      "Give the files a condition in the table above to compare them."))
    st <- percell_condition_stats(d, met)
    xi <- match(d$condition, st$condition)
    set.seed(1L); jit <- (stats::runif(nrow(d)) - 0.5) * 0.3     # horizontal spread only
    dark   <- isTRUE(isolate(input$trk_dark_mode))
    bg_col <- if (dark) "#000000" else "#FFFFFF"
    fg_col <- if (dark) "#FFFFFF" else "#333333"
    ci_up  <- ifelse(is.finite(st$ci_hi), st$ci_hi - st$mean, 0)
    ci_dn  <- ifelse(is.finite(st$ci_lo), st$mean - st$ci_lo, 0)
    plot_ly() %>%
      add_trace(x = xi + jit, y = d[[met]], type = "scatter", mode = "markers",
                marker = list(size = 9, color = "rgba(70,130,200,0.75)",
                              line = list(color = "rgba(40,90,160,0.9)", width = 1)),
                text = d$file, name = "one cell",
                hovertemplate = "%{text}<br>%{y:.4g}<extra></extra>") %>%
      add_trace(x = seq_len(nrow(st)), y = st$mean, type = "scatter", mode = "markers",
                marker = list(symbol = "line-ew-open", size = 34, color = fg_col,
                              line = list(width = 3, color = fg_col)),
                error_y = list(type = "data", symmetric = FALSE, array = ci_up, arrayminus = ci_dn,
                               color = fg_col, thickness = 1.5, width = 10),
                name = "mean ± 95% CI across cells",
                hovertemplate = "mean %{y:.4g}<extra></extra>") %>%
      layout(
        xaxis = list(title = "", tickmode = "array", tickvals = seq_len(nrow(st)),
                     ticktext = sprintf("%s<br>(%d cell%s)", st$condition, st$n_cells,
                                        ifelse(st$n_cells == 1L, "", "s")),
                     range = c(0.4, nrow(st) + 0.6), color = fg_col, zeroline = FALSE),
        yaxis = list(title = PERCELL_METRICS[[met]], color = fg_col,
                     rangemode = if (met == "r1_excess") "normal" else "tozero"),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col, font = list(color = fg_col),
        legend = list(orientation = "h", y = -0.22), uirevision = "trk_cmp")
  })

  # Condition table for the chosen metric, plus the flags that change how the
  # numbers should be read.
  output$trk_cmp_stats <- renderUI({
    r <- cmp_view(); req(!is.null(r))
    met <- input$trk_cmp_metric %||% "D_um2_s"
    st  <- percell_condition_stats(r, met)
    fmt <- function(x) ifelse(is.finite(x), formatC(x, digits = 3, format = "g"), "-")
    AMB <- "#B8860B"
    asg <- nzchar(trimws(r$condition))
    n_smooth <- sum(asg & is.finite(r$r1_excess) & r$r1_excess > 0.01)
    n_dw     <- sum(asg & is.finite(r$Dw_frac) & r$Dw_frac > 0.5)
    tagList(
      if (isTRUE(attr(r, "stale"))) helpText(style = sprintf("color:%s; font-size:0.85em;", AMB),
        "Sequence settings were edited after the last run: press Analyse each cell again."),
      if (nrow(st)) tags$table(class = "table table-condensed",
        style = "font-size:0.85em; margin:6px 0; width:auto;",
        tags$thead(tags$tr(lapply(c("condition", "cells", "mean", "SD", "95% CI"), tags$th))),
        tags$tbody(lapply(seq_len(nrow(st)), function(i) tags$tr(
          tags$td(st$condition[i]), tags$td(st$n_cells[i]), tags$td(fmt(st$mean[i])),
          tags$td(fmt(st$sd[i])), tags$td(sprintf("%s to %s", fmt(st$ci_lo[i]), fmt(st$ci_hi[i]))))))),
      if (nrow(st) && any(st$n_cells < 3L)) helpText(style = "font-size:0.82em; color:#666;",
        "A condition needs at least 3 cells for a meaningful CI."),
      if (n_smooth) helpText(style = sprintf("color:%s; font-size:0.82em;", AMB), sprintf(
        "%d cell(s) have increments correlated beyond any motion model (r₁ above the bound). Galvo damping causes this, and so does the tracking loop when motion is small, so compare D and σ only between cells acquired with the same damping setting.",
        n_smooth)),
      if (n_dw) helpText(style = sprintf("color:%s; font-size:0.82em;", AMB), sprintf(
        "In %d cell(s) the wander floor is more than half of the measured D.", n_dw)),
      if (sum(!asg)) helpText(style = "font-size:0.82em; color:#666;", sprintf(
        "%d file(s) without a condition are listed below and left out of the comparison.", sum(!asg))))
  })

  output$trk_cmp_results <- DT::renderDT({
    r <- cmp_view(); req(!is.null(r))
    cols <- c("condition", "file", "calibration", "n_tracks", "len_med", "t_c_ms", "frac_single", "eco_med",
              "sigma_nm", "a", "D_um2_s", "D_corr_um2_s", "Dw_frac", "alpha_short", "r1_excess", "note")
    show <- r[, cols, drop = FALSE]
    num  <- vapply(show, is.numeric, logical(1))
    show[num] <- lapply(show[num], function(v) signif(v, 3))
    DT::datatable(show, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE, scrollX = TRUE,
                                 scrollY = "300px", scrollCollapse = TRUE))
  })

  output$trk_cmp_results_dl <- downloadHandler(
    filename = function() paste0("per_cell_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      r <- cmp_view(); if (is.null(r)) r <- data.frame()
      utils::write.csv(r, file, row.names = FALSE)
    })
} # end tracking_server


