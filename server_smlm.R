## ─────────────────────────────────────────────────────────────────────────────
## server_smlm.R  -  SMLM analysis server logic
## Called from server.R inside the combined server function.
## ─────────────────────────────────────────────────────────────────────────────
##
## GOAL OF THIS FILE
##   Everything behind the SMLM (imaging) tab: load localization data, filter it,
##   group it into molecules and structures, and run the spatial analyses
##   (nearest neighbour, Ripley, clustering, shape fitting, two-colour
##   co-localization, ground-truth estimation). One large Shiny server function.
##
##   The Tracking tab asks how molecules MOVE. This tab asks where they ARE and
##   how they are ARRANGED.
##
## HOW TO READ A SHINY SERVER FILE
##   Reactives are LAZY - file order is NOT execution order. Follow the CALLS.
##   The dependency chain is drawn below.
##
## ── 1. THE DATA FLOW ─────────────────────────────────────────────────────────
##
##   upload (.json / .csv / .npy)
##        |   readers + canonicalisation live in helpers_shared.R
##        v
##   uploaded_df_rv -> uploaded_df() ..... the raw, unfiltered dataset
##        v
##   ---- THE FILTER PIPELINE (three stages, each feeds the next) --------------
##   metadata_filtered_df .... Stage 1  per-LOCALIZATION quality ranges
##                                      (efo / cfr / dcr / eco / tim ...).
##                                      Drops individual points.
##   spread_filtered_df ...... Stage 2  per-TRACE spread (SD) filter.
##                                      Drops whole traces that are too diffuse
##                                      - i.e. a "trace" that is really drift or
##                                      two molecules merged.
##   filtered_df ............. Stage 3  per-TRACE localization count.
##                                      Drops traces that are too short/long.
##        v
##   plottable_data() ........ THE analysis-ready dataset. Adds the trimmer, the
##                             2D crop, and - importantly - runs the CLUSTERING.
##                             Almost every analysis starts by calling this.
##
##   ORDER MATTERS AND IS DELIBERATE: rows are filtered first (Stage 1), then
##   whole traces are judged on what survived (Stages 2-3). Filtering traces
##   before cleaning their points would judge them on data you already rejected.
##
## ── 2. THE GROUPING HIERARCHY (the key concept of this tab) ──────────────────
##   localizations -> tid -> cluster (primary) -> cluster2 (secondary)
##
##     tid ....... one continuous blink of one fluorophore
##     cluster ... DE-BLINKING. Merges repeat appearances of one molecule into
##                 one object, at a threshold of a few times the measured trace
##                 spread (~10-30 nm). Complete-linkage on trace centroids, or
##                 a majority vote of the localizations - NOT DBSCAN. It merges
##                 and never removes, so this column has NO id 0.
##     cluster2 .. MEDIAN CLUSTERING over the medians above, LARGE eps
##                 (~100-300 nm), by DBSCAN / HDBSCAN / OPTICS. Job: STRUCTURE -
##                 group molecules into the biological object they form.
##
##   WHY IT EXISTS: a fluorophore blinks many times, and an antibody carries
##   several dyes 10-40 nm from its target. Statistics run on RAW localizations
##   therefore measure the LABEL, not the biology. In `cluster2`, id 0 means
##   noise/unassigned and must be excluded before counting; `cluster` has no 0
##   because de-blinking assigns every trace (see deblink_assign_singletons).
##   (Implementation: apply_clustering_pipeline() in helpers_smlm.R.)
##
## ── 3. MAP OF THE ANALYSES (grep the name; line numbers drift) ───────────────
##   Upload & settings ......... uploaded_df, get_ui_settings, raw_df, df_full_raw
##   Filter pipeline ........... metadata_filtered_df / spread_filtered_df /
##                               filtered_df   (+ the sd_* trace-spread helpers)
##   Analysis-ready data ....... plottable_data
##   Nearest neighbour ......... nn_table, nn_table_cluster2, secondary_nn_lines,
##                               nn_batch_filtered (batch across files/channels)
##   Two-colour / multiplex .... mplex_split_active, mplex_cross_pair,
##                               cross_dcr_nn_table
##   Ground-truth estimate ..... gt_estimate, gt_morphology   [EXPLORATORY]
##   Cluster viewer / shapes ... viewer_ids, viewer_current, viewer_fit_core,
##                               viewer_fit                    [EXPLORATORY]
##   Dataset groups ............ file_group_map, group_color_map
##   Image overlay ............. confocal_meta, image_overlay
##
##   The heavy spatial statistics themselves (Ripley L/L12, PCF, KDE
##   co-localization, RoG, batch NN) live in helpers_smlm.R - this file wires
##   them to inputs and outputs.
##
## ── 4. CONVENTIONS & GOTCHAS ─────────────────────────────────────────────────
##   * Units: nanometres. `tim` stays in SECONDS here (the Tracking tab converts
##     to ms; this tab does not).
##   * plottable_data is an eventReactive on plot_trigger(): analyses refresh
##     when Plot is pressed, not on every keystroke.
##   * Multiplex: when active, clustering is forced to run PER CHANNEL so
##     molecules from different exchange rounds are never merged.
##   * Retired code is in deadcode_smlm.R (NOT sourced). Several reactives were
##     previously defined TWICE at this scope - the later definition silently
##     shadowed the earlier, so the first never ran. Before adding a reactive,
##     check the name is not already used further down the file.
##
##   Method provenance and MINFLUX-appropriateness of each statistic:
##   docs/MINFLUX_analysis_methods.md (Part II §6).
## ─────────────────────────────────────────────────────────────────────────────

smlm_server <- function(input, output, session) {
  smlm_axis_x_span  <- reactiveVal(NULL)
  smlm_zoom_ranges  <- reactiveVal(NULL)   # list(x=c(x0,x1), y=c(y0,y1)) from 2D plot zoom
  smlm_2d_crop_rv   <- reactiveVal(NULL)   # active 2D crop box: list(x=c(x0,x1), y=c(y0,y1))
  # Committed 3D crop region. The box/sphere inputs are only a PREVIEW; the crop
  # is applied when the user presses "Apply crop", which commits the region here.
  # list(mode="box", x=c(lo,hi), y=..., z=...) or list(mode="sphere", c=c(cx,cy,cz), r=r).
  trim_applied_rv   <- reactiveVal(NULL)

  smlm_default_x_span <- function() {
    df <- isolate(plottable_data())
    if (!is.data.frame(df) || !nrow(df)) return(NULL)

    lim <- get_equal_limits_2d(df$x, df$y)
    ovl <- tryCatch(isolate(image_overlay()), error = function(e) NULL)
    lim <- expand_limits_with_image(lim, ovl, invert_y = isTRUE(isolate(input$img_invert_y)))

    span <- suppressWarnings(as.numeric(diff(lim$x)))
    if (length(span) == 1 && is.finite(span) && span > 0) span else NULL
  }

  smlm_push_scale_bar <- function() {
    proxy <- plotlyProxy("scatter2D", session)
    if (!isTRUE(input$show_scale_bar) || !identical(input$mode3d, "2D")) {
      plotlyProxyInvoke(proxy, "relayout", list(shapes = list(), annotations = list()))
      return()
    }

    span <- smlm_axis_x_span()
    if (length(span) != 1 || !is.finite(span) || span <= 0) {
      span <- smlm_default_x_span()
      if (length(span) == 1 && is.finite(span) && span > 0) smlm_axis_x_span(span)
    }

    sb_len   <- suppressWarnings(as.numeric(input$scale_bar_length))
    sb_thick <- suppressWarnings(as.numeric(input$scale_bar_thickness))
    sb_font  <- suppressWarnings(as.numeric(input$scale_font_size))
    sb <- make_scale_bar(sb_len, span, sb_thick, dark = isTRUE(input$dark_mode),
                         scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
    if (!is.null(sb)) {
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = sb$shapes, annotations = sb$annotations))
    }
  }

  # ── 3D auto-rotate ──────────────────────────────────────────────────────────
  observeEvent(
    list(input$auto_rotate_3d, input$rotate_speed_3d, input$rotate_elev_3d, input$rotate_dist_3d),
    ignoreNULL = FALSE, {
      if (isTRUE(input$auto_rotate_3d)) {
        speed <- input$rotate_speed_3d %||% 2
        elev  <- input$rotate_elev_3d  %||% 25
        dist  <- input$rotate_dist_3d  %||% 2.5
        shinyjs::runjs(sprintf("startRotate3D(%s, %s, %s);", speed, elev, dist))
      } else {
        shinyjs::runjs("stopRotate3D();")
      }
    }
  )

  # ── 3D recording ────────────────────────────────────────────────────────────
  observeEvent(input$record_start_3d, {
    w       <- input$record_width_3d   %||% 800
    h       <- input$record_height_3d  %||% 600
    fps     <- input$record_fps_3d     %||% 20
    bitrate <- input$record_bitrate_3d %||% 10
    shinyjs::runjs(sprintf("startRecord3D(%d, %d, %d, %s);", as.integer(w), as.integer(h), as.integer(fps), bitrate))
  })

  observeEvent(input$record_stop_3d, {
    shinyjs::runjs("stopRecord3D();")
  })

  #1 ) UPLOAD AND PREPROCESS DATA ----------
  uploaded_df_rv  <- reactiveVal(NULL)
  smlm_file_store <- reactiveVal(list())   # fname → data frame

  # ── Axial (z) calibration ────────────────────────────────────────────────
  # MINFLUX z is multiplied by a scaling factor at import to correct the
  # refractive-index mismatch between immersion medium and sample. pyMINFLUX
  # ships 0.7 as a default; the physically correct value depends on the actual
  # index pair and imaging depth, and 0.67 is also in use - so it is exposed
  # rather than hard-coded.
  #
  # It multiplies EVERY axial distance downstream: z trace spread, cluster
  # extent in z, 3D nearest-neighbour distances, Rg. Treat it as a calibration
  # to be measured on a known standard, not a display preference.
  SMLM_Z_SCALE_DEFAULT <- 0.7
  # The factor currently baked into smlm_file_store() - the single source of
  # truth. input$z_scale is only a PENDING value until Plot Data commits it, so
  # nothing may read the input directly to decide what the data is scaled by.
  z_scale_applied <- reactiveVal(SMLM_Z_SCALE_DEFAULT)

  # Helper: load one file from a temp path, return df or NULL
  .load_one_file <- function(path, fname) {
    ext <- get_ext(fname)
    tryCatch({
      withProgress(message = paste0("Loading ", fname, "…"), value = 0, {
        incProgress(0.03, detail = "Preparing")
        loaded <- switch(
          ext,
          "csv" = {
            incProgress(0.20, detail = "Reading CSV...")
            read.csv(path, check.names = FALSE)
          },

          "json" = {
            incProgress(0.10, detail = "Reading JSON...")
            read_minflux_json_as_df(
              path,
              keep_final_only = TRUE,
              keep_valid_only = TRUE,
              final_itr = NA_integer_,
              cfr_strategy = "last_nonzero",
              # Read RAW: .npy/.json come straight off the instrument with no
              # axial correction applied. The app-wide factor is applied once,
              # below, for every source alike.
              z_scale = 1.0,
              dwell_time_ms = 1.0,
              pool_dcr = isTRUE(input$pool_dcr_json)
            )
          },
          "npy" = {
            incProgress(0.10, detail = "Reading .npy (binary)")
            read_minflux_npy_as_df(
              path,
              keep_final_only = TRUE,
              keep_valid_only = TRUE,
              final_itr = NA_integer_,
              cfr_strategy = "last_nonzero",
              # Read RAW: .npy/.json come straight off the instrument with no
              # axial correction applied. The app-wide factor is applied once,
              # below, for every source alike.
              z_scale = 1.0,
              dwell_time_ms = 1.0,
              pool_dcr = isTRUE(input$pool_dcr_json)
            )
          },
          stop("Unsupported file type: ", ext)
        )
        if (!is.data.frame(loaded) || !nrow(loaded)) stop("No rows loaded from file.")
        incProgress(0.15, detail = "Canonicalizing columns")
        loaded <- canonicalize_minflux_df(loaded)
        if (!"z" %in% names(loaded)) loaded$z <- 0

        # ── Normalize the axial calibration across all three sources ──────────
        # Everything reaching this point is reduced to RAW (uncorrected) z, then
        # the single app-wide factor is applied. Without this, .npy/.json would
        # follow the calibration control while CSVs stayed frozen at whatever
        # their producer baked in - and a CSV would get 0.7 applied twice.
        #
        #   .npy / .json  raw already (read above with z_scale = 1)
        #   .csv          divide out the factor its producer applied
        if (identical(ext, "csv")) {
          loaded  <- strip_csv_z_calibration(loaded)     # -> raw
          src_cal <- attr(loaded, "z_cal")
        } else {
          src_cal <- list(factor = 1.0, source = "raw instrument data", certain = TRUE)
        }
        # Applied at the COMMITTED factor, not the pending input value, so the
        # store is never left at a factor the user has not confirmed.
        loaded$z <- loaded$z * z_scale_applied()         # -> app calibration
        attr(loaded, "z_cal") <- src_cal

        # Say nothing about the axial calibration for planar data - a 2D upload
        # should not produce messages about z. has_axial_data() tests VARIATION,
        # not presence: 2D exports still carry a z column, filled with zeros.
        if (identical(ext, "csv") && isTRUE(src_cal$has_z)) {
          showNotification(
            sprintf("%s: axial calibration %.3f detected (%s)%s.",
                    fname, src_cal$factor, src_cal$source,
                    if (isTRUE(src_cal$certain)) "" else " - assumed"),
            type = "message", duration = 8)
        }
        incProgress(0.07, detail = "Finalizing")
        loaded
      })
    }, error = function(e) {
      showNotification(paste0("Could not load ", fname, ": ", conditionMessage(e)),
                       type = "error", duration = 15)
      NULL
    })
  }

  observeEvent(input$file, {
    req(input$file)
    store     <- smlm_file_store()
    new_names <- character(0)

    adopt_factor <- NULL   # a recorded (not guessed) factor from an uploaded file

    for (i in seq_len(nrow(input$file))) {
      fname <- input$file$name[i]
      if (fname %in% names(store)) next   # already loaded - skip duplicate

      df <- .load_one_file(input$file$datapath[i], fname)
      if (is.null(df)) next

      # A file of ours records the factor it was written at. Adopt it so the
      # control reflects the data being loaded and z comes back exactly as
      # exported. Only `certain` factors are adopted - the 0.7 assumed for
      # pyMINFLUX and legacy exports is a guess and must not silently become
      # the session's calibration.
      # ...and only when the file actually has axial data: adopting a factor for
      # a planar dataset would change a number that affects nothing, while
      # announcing it.
      zc <- attr(df, "z_cal")
      if (!is.null(zc) && isTRUE(zc$certain) && isTRUE(zc$has_z) &&
          is.finite(zc$factor) && zc$factor > 0)
        adopt_factor <- zc$factor

      df$.dataset_tag <- paste0(sample(c(0:9, letters[1:6]), 4, replace = TRUE), collapse = "")
      store[[fname]] <- df
      new_names <- c(new_names, fname)

      if (identical(get_ext(fname), "json")) {
        final_used <- if ("itr" %in% names(df)) max(df$itr, na.rm = TRUE) else NA_integer_
        showNotification(
          sprintf("%s - JSON import: using final localization (itr=%s); CFR carried over per TID.",
                  fname, final_used),
          type = "message", duration = 8)
      }
    }

    if (!length(new_names)) return()
    smlm_file_store(store)

    # Directly update uploaded_df_rv here so plotting works even when
    # input$active_file doesn't change (e.g. re-uploading the same filename
    # after removal, where the client-side value may still equal new_names[1]
    # and the observeEvent(input$active_file) observer never fires).
    selected_df <- store[[new_names[1]]]
    if (!is.null(selected_df)) {
      uploaded_df_rv(selected_df)
      hist_base <- c("efo", "tid", "cfr", "tim", "eco", "dcr", "dwell", "fluo")
      hist_choices <- if ("dcr_pooled" %in% names(selected_df))
        c(hist_base, "dcr_pooled") else hist_base
      updateSelectInput(session, "hist_variable", choices = hist_choices,
                        selected = isolate(input$hist_variable) %||% "efo")
    }

    updateSelectInput(session, "active_file",
                      choices  = names(store),
                      selected = new_names[1])

    # Adopt a recorded axial factor, once, at upload. commit_z_scale() rescales
    # everything already in the store to match, so the newly loaded file ends up
    # at exactly the z it was exported with.
    if (!is.null(adopt_factor)) {
      updateNumericInput(session, "z_scale", value = adopt_factor)
      commit_z_scale(adopt_factor, reason = "Axial calibration from file")
    }
  }, ignoreInit = TRUE)

  # Switch active dataset when the file selector changes
  observeEvent(input$active_file, {
    req(input$active_file, nzchar(input$active_file))
    df <- smlm_file_store()[[input$active_file]]
    if (is.null(df)) return()

    # Snapshot sd_filt BEFORE the data change so we can restore it after the
    # reactive graph settles.  This is the authoritative lock mechanism for SD
    # ranges - more reliable than the timing-based sd_skip_arming approach.
    snap_sd <- if (isTRUE(isolate(input$lock_filters)))
      isolate(list(x = sd_filt$x, y = sd_filt$y, z = sd_filt$z))
    else NULL

    uploaded_df_rv(df)

    # Update histogram variable choices - show dcr_pooled only when present
    hist_base <- c("efo", "tid", "cfr", "tim", "eco", "dcr", "dwell", "fluo")
    hist_choices <- if ("dcr_pooled" %in% names(df))
      c(hist_base, "dcr_pooled") else hist_base
    updateSelectInput(session, "hist_variable", choices = hist_choices,
                      selected = isolate(input$hist_variable) %||% "efo")

    if (!is.null(snap_sd)) {
      session$onFlushed(function() {
        sd_skip_arming(TRUE)
        sd_filt$x <- snap_sd$x
        sd_filt$y <- snap_sd$y
        sd_filt$z <- snap_sd$z
        # Update the widgets to display the restored values
        shinyWidgets::updateNumericRangeInput(session, "sd_x_range",
          value = c(max(0, snap_sd$x[1]), if (is.infinite(snap_sd$x[2])) 1e6 else snap_sd$x[2]))
        shinyWidgets::updateNumericRangeInput(session, "sd_y_range",
          value = c(max(0, snap_sd$y[1]), if (is.infinite(snap_sd$y[2])) 1e6 else snap_sd$y[2]))
        shinyWidgets::updateNumericRangeInput(session, "sd_z_range",
          value = c(max(0, snap_sd$z[1]), if (is.infinite(snap_sd$z[2])) 1e6 else snap_sd$z[2]))
        later::later(function() sd_skip_arming(FALSE), delay = 0.5)
      }, once = TRUE)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$flush_file, {
    fname <- isolate(input$active_file)
    if (!nzchar(fname %||% "")) return()
    store <- smlm_file_store()
    store[[fname]] <- NULL
    smlm_file_store(store)
    remaining <- names(store)
    if (length(remaining) == 0) {
      uploaded_df_rv(NULL)
      updateSelectInput(session, "active_file", choices = character(0))
    } else {
      updateSelectInput(session, "active_file",
                        choices  = remaining,
                        selected = remaining[1])
    }
  })

  # ── MULTIPLEXED MINFLUX ─────────────────────────────────────────────────────
  # Build an N-channel combined dataset (overlay loaded files, or split the active
  # file by TIM breakpoints) carrying mplex_channel / mplex_label / mplex_color.
  # Registered as one synthetic entry in smlm_file_store so the plot and every
  # Analysis tool operate on it unchanged. See plan: graceful-marinating-pudding.
  mplex_name   <- "★ Multiplex"
  MPLEX_TID_OFF <- 1e10            # per-channel tid offset (> u4 tid max ~4.3e9)
  mplex_spec_rv <- reactiveVal(NULL)  # list(channel,label,color,n_loc,n_trace)

  MPLEX_MAX_BREAKS <- 6L              # up to 6 splits → 7 time-window channels
  MPLEX_MAX_WIN    <- 7L              # up to 7 explicit windows (= 7 channels)
  .mplex_parse_breaks <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(numeric(0))
    v <- suppressWarnings(as.numeric(strsplit(txt, "[,;[:space:]]+")[[1]]))
    v <- sort(unique(v[is.finite(v)]))
    utils::head(v, MPLEX_MAX_BREAKS)
  }
  # Parse "start-end, start-end, …" into a sorted N×2 matrix of included windows.
  .mplex_parse_windows <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(matrix(numeric(0), ncol = 2))
    parts <- strsplit(txt, "[,;]")[[1]]
    rows <- lapply(parts, function(p) {
      v <- suppressWarnings(as.numeric(strsplit(trimws(p), "[-–:[:space:]]+")[[1]]))
      v <- v[is.finite(v)]
      if (length(v) >= 2) sort(v[1:2]) else NULL
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(matrix(numeric(0), ncol = 2))
    m <- do.call(rbind, rows)
    m <- m[order(m[, 1]), , drop = FALSE]
    utils::head(m, MPLEX_MAX_WIN)
  }
  # Assign a tim vector to channel indices per the active TIM method.
  # Returns list(ch = integer channel per tim (NA = removed/wash), nch, method,
  # breaks/wins). Used by the colors UI, the preview histogram and the builder.
  .mplex_tim_assign <- function(tim) {
    method <- input$mplex_tim_method %||% "breaks"
    if (identical(method, "windows")) {
      wins <- .mplex_parse_windows(input$mplex_tim_windows)
      ch <- rep(NA_integer_, length(tim))
      for (w in seq_len(nrow(wins))) {
        inw <- is.na(ch) & tim >= wins[w, 1] & tim <= wins[w, 2]
        ch[inw] <- w
      }
      list(ch = ch, nch = nrow(wins), method = "windows", wins = wins)
    } else {
      br <- .mplex_parse_breaks(input$mplex_tim_breaks)
      list(ch = findInterval(tim, br) + 1L, nch = length(br) + 1L,
           method = "breaks", breaks = br)
    }
  }

  # Keep the overlay file picker in sync with loaded files (excluding the
  # synthetic multiplex entry itself).
  observeEvent(smlm_file_store(), {
    nm <- setdiff(names(smlm_file_store()), mplex_name)
    updateSelectInput(session, "mplex_files", choices = nm,
                      selected = intersect(isolate(input$mplex_files), nm))
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  output$mplex_tim_info <- renderUI({
    df <- uploaded_df_rv()
    if (is.null(df) || !"tim" %in% names(df) || !nrow(df))
      return(helpText("Load a file with a tim column to split by time."))
    rng <- range(df$tim, na.rm = TRUE)
    helpText(sprintf("Active file tim range: %.2f – %.2f s (span %.2f s).",
                     rng[1], rng[2], diff(rng)))
  })

  output$mplex_files_colors_ui <- renderUI({
    fs <- input$mplex_files
    if (is.null(fs) || !length(fs))
      return(helpText("Pick one or more files above; each becomes a channel."))
    defs <- mplex_default_colors(length(fs))
    lapply(seq_along(fs), function(j) {
      lid <- paste0("mplex_flab_", j); cid <- paste0("mplex_fcol_", j)
      cur_l <- isolate(input[[lid]]); cur_c <- isolate(input[[cid]])
      fluidRow(
        column(7, textInput(lid, NULL,
          value = if (!is.null(cur_l) && nzchar(cur_l)) cur_l else fs[j])),
        column(5, colourpicker::colourInput(cid, NULL,
          value = if (!is.null(cur_c) && nzchar(cur_c)) cur_c else defs[j]))
      )
    })
  })

  output$mplex_tim_colors_ui <- renderUI({
    method <- input$mplex_tim_method %||% "breaks"
    df  <- isolate(uploaded_df_rv())
    has_tim <- !is.null(df) && "tim" %in% names(df) && nrow(df)
    if (identical(method, "windows")) {
      wins <- .mplex_parse_windows(input$mplex_tim_windows)
      nch  <- nrow(wins)
      deflab_fn <- function(j) sprintf("Ch%d (%.0f–%.0f s)", j, wins[j, 1], wins[j, 2])
    } else {
      br  <- .mplex_parse_breaks(input$mplex_tim_breaks)
      nch <- length(br) + 1L
      edges <- c(if (has_tim) min(df$tim, na.rm = TRUE) else 0, br,
                 if (has_tim) max(df$tim, na.rm = TRUE) else NA_real_)
      deflab_fn <- function(j) sprintf("Ch%d (%.0f–%.0f s)", j, edges[j], edges[j + 1])
    }
    if (nch < 1L) return(helpText("Define at least one window / breakpoint above."))
    defs <- mplex_default_colors(nch)
    lapply(seq_len(nch), function(j) {
      lid <- paste0("mplex_tlab_", j); cid <- paste0("mplex_tcol_", j)
      cur_l <- isolate(input[[lid]]); cur_c <- isolate(input[[cid]])
      fluidRow(
        column(7, textInput(lid, NULL,
          value = if (!is.null(cur_l) && nzchar(cur_l)) cur_l else deflab_fn(j))),
        column(5, colourpicker::colourInput(cid, NULL,
          value = if (!is.null(cur_c) && nzchar(cur_c)) cur_c else defs[j]))
      )
    })
  })

  # Traces-over-time preview: one point per trace (median tim) binned over time,
  # bars colored by the time-window channel they fall in, dashed lines at splits.
  # Helps the user place breakpoints. Live on the active file + breakpoints input.
  output$mplex_tim_hist <- plotly::renderPlotly({
    df <- uploaded_df_rv()
    shiny::validate(shiny::need(!is.null(df) && "tim" %in% names(df) && nrow(df),
                                "Load a file with a tim column."))
    tt <- df %>% dplyr::group_by(tid) %>%
      dplyr::summarise(tim = stats::median(tim, na.rm = TRUE), .groups = "drop")
    tt <- tt$tim[is.finite(tt$tim)]
    shiny::validate(shiny::need(length(tt) > 1 && diff(range(tt)) > 0,
                                "Not enough time spread to plot."))
    sp   <- .mplex_tim_assign(tt)
    nch  <- max(1L, sp$nch)
    cols <- mplex_default_colors(nch)
    for (j in seq_len(nch)) {
      v <- input[[paste0("mplex_tcol_", j)]]
      if (!is.null(v) && nzchar(v)) cols[j] <- v
    }
    rng   <- range(tt)
    nbins <- max(2L, as.integer(input$mplex_tim_bins %||% 60L))
    edges <- seq(rng[1], rng[2], length.out = nbins + 1L)
    h     <- graphics::hist(tt, breaks = edges, plot = FALSE)
    binch <- .mplex_tim_assign(h$mids)$ch         # channel per bin (NA = wash)
    barcols <- ifelse(is.na(binch), "#cccccc", cols[pmin(binch, nch)])
    # boundary lines: window edges (windows mode) or breakpoints (breaks mode)
    bounds <- if (identical(sp$method, "windows")) as.numeric(t(sp$wins)) else sp$breaks
    shapes <- lapply(bounds, function(b) list(
      type = "line", x0 = b, x1 = b, y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "#333333", width = 2, dash = "dash")))
    plotly::plot_ly(
      x = h$mids, y = h$counts, type = "bar",
      marker = list(color = barcols),
      hovertemplate = "t = %{x:.0f} s<br>%{y} traces<extra></extra>"
    ) %>% plotly::layout(
      bargap = 0,
      xaxis = list(title = "time (s)"),
      yaxis = list(title = "# trace IDs"),
      showlegend = FALSE,
      margin = list(t = 10, b = 36, l = 48, r = 8),
      shapes = shapes
    )
  })

  observeEvent(input$mplex_build, {
    mode <- input$mplex_mode %||% "files"
    spec <- list()

    if (mode == "files") {
      fs <- input$mplex_files
      if (is.null(fs) || length(fs) < 1) {
        showNotification("Pick at least one file to overlay.", type = "warning"); return()
      }
      store <- smlm_file_store()
      defs  <- mplex_default_colors(length(fs))
      parts <- list()
      for (j in seq_along(fs)) {
        d <- store[[fs[j]]]
        if (is.null(d) || !nrow(d)) next
        lab <- input[[paste0("mplex_flab_", j)]]
        col <- input[[paste0("mplex_fcol_", j)]]
        if (is.null(lab) || !nzchar(lab)) lab <- fs[j]
        if (is.null(col) || !nzchar(col)) col <- defs[j]
        d <- d %>% dplyr::select(-dplyr::any_of(
          c("mplex_channel", "mplex_label", "mplex_color", ".dataset_tag")))
        d$tid <- suppressWarnings(as.numeric(d$tid)) + j * MPLEX_TID_OFF
        d$mplex_channel <- j; d$mplex_label <- lab; d$mplex_color <- col
        parts[[length(parts) + 1L]] <- d
        spec[[length(spec) + 1L]] <- list(channel = j, label = lab, color = col,
          n_loc = nrow(d), n_trace = dplyr::n_distinct(d$tid))
      }
      if (!length(parts)) {
        showNotification("Selected files contained no data.", type = "error"); return()
      }
      combined <- dplyr::bind_rows(parts)

    } else {
      d <- isolate(uploaded_df_rv())
      if (is.null(d) || !nrow(d)) {
        showNotification("No active file to split.", type = "warning"); return()
      }
      if (!"tim" %in% names(d)) {
        showNotification("Active file has no tim column.", type = "error"); return()
      }
      d <- d %>% dplyr::select(-dplyr::any_of(
        c("mplex_channel", "mplex_label", "mplex_color")))
      sp  <- .mplex_tim_assign(d$tim)
      nch <- sp$nch
      if (nch < 1L || (identical(sp$method, "breaks") && nch < 2L)) {
        showNotification(
          if (identical(sp$method, "windows"))
            "Enter at least one time window (start-end)."
          else "Enter at least one TIM breakpoint (seconds).",
          type = "warning"); return()
      }
      ch   <- sp$ch
      keep <- !is.na(ch)            # drop wash localizations (NA = outside windows)
      if (!any(keep)) {
        showNotification("No localizations fall inside the selected windows.",
                         type = "error"); return()
      }
      d  <- d[keep, , drop = FALSE]
      ch <- ch[keep]
      defs <- mplex_default_colors(nch)
      labs <- vapply(seq_len(nch), function(j) {
        v <- input[[paste0("mplex_tlab_", j)]]
        if (is.null(v) || !nzchar(v)) paste0("Ch", j) else v }, character(1))
      cols <- vapply(seq_len(nch), function(j) {
        v <- input[[paste0("mplex_tcol_", j)]]
        if (is.null(v) || !nzchar(v)) defs[j] else v }, character(1))
      d$mplex_channel <- ch
      d$mplex_label   <- labs[ch]
      d$mplex_color   <- cols[ch]
      d$tid <- suppressWarnings(as.numeric(d$tid)) + ch * MPLEX_TID_OFF
      combined <- d
      spec <- lapply(seq_len(nch), function(j) {
        sub <- d[d$mplex_channel == j, , drop = FALSE]
        list(channel = j, label = labs[j], color = cols[j],
             n_loc = nrow(sub), n_trace = dplyr::n_distinct(sub$tid))
      })
    }

    if (!"z" %in% names(combined)) combined$z <- 0
    combined$.dataset_tag <- paste0(
      sample(c(0:9, letters[1:6]), 4, replace = TRUE), collapse = "")

    st <- smlm_file_store(); st[[mplex_name]] <- combined; smlm_file_store(st)
    mplex_spec_rv(spec)
    uploaded_df_rv(combined)
    updateSelectInput(session, "active_file", choices = names(st), selected = mplex_name)
    updateSelectInput(session, "color_mode", selected = "multiplex")
    plot_trigger(Sys.time())
    showNotification(sprintf("Built multiplex: %d channel(s), %d localizations.",
                             length(spec), nrow(combined)),
                     type = "message", duration = 6)
  }, ignoreInit = TRUE)

  observeEvent(input$mplex_clear, {
    st <- smlm_file_store()
    if (!is.null(st[[mplex_name]])) { st[[mplex_name]] <- NULL; smlm_file_store(st) }
    mplex_spec_rv(NULL)
    rem <- names(smlm_file_store())
    if (length(rem)) {
      updateSelectInput(session, "active_file", choices = rem, selected = rem[1])
      uploaded_df_rv(smlm_file_store()[[rem[1]]])
    } else {
      updateSelectInput(session, "active_file", choices = character(0))
      uploaded_df_rv(NULL)
    }
    if (identical(isolate(input$color_mode), "multiplex"))
      updateSelectInput(session, "color_mode", selected = "default")
  })

  output$mplex_summary <- renderText({
    sp <- mplex_spec_rv()
    if (is.null(sp) || !length(sp)) return("No multiplex built yet.")
    paste(vapply(sp, function(c)
      sprintf("Ch%d  %s:  %s loc,  %s traces",
              c$channel, c$label,
              format(c$n_loc, big.mark = ","), format(c$n_trace, big.mark = ",")),
      character(1)), collapse = "\n")
  })

  # Split the active multiplex dataset into a named list (label → channel df) so
  # the per-file Analysis loops produce one result per channel. NULL when the
  # active dataset is not a multiplex. Used by the *_by_channel toggles.
  mplex_split_active <- reactive({
    df <- uploaded_df_rv()
    if (is.null(df) || !"mplex_label" %in% names(df) || !nrow(df)) return(NULL)
    ord  <- order(df$mplex_channel)
    labs <- unique(df$mplex_label[ord])
    stats::setNames(
      lapply(labs, function(l) df[df$mplex_label == l, , drop = FALSE]), labs)
  })

  # ── Split one prepared file into its channels, for the spatial statistics ──
  #
  # Covers BOTH channel sources - a 2-colour DCR split and a multiplex build -
  # because they are one concept now. The older `mplex_split_active()` above is
  # multiplex-only and splits the RAW upload; this one splits a frame that has
  # already been through resolve_channels(), so the labels and colours match
  # everything else on screen.
  #
  # THE POOLED PARENT COMES BACK WITH EACH CHANNEL, and that is the point.
  # Ripley's L and g(r) divide by the observed area. If each channel derived its
  # own window, the sparser channel would get a smaller bounding box and
  # footprint, hence a higher apparent density, hence a DEFLATED curve - a
  # difference between the channels manufactured out of their difference in
  # count. The two dyes imaged the SAME field, so they must be measured in the
  # same window; then a real difference in curve is a difference in how the two
  # targets are arranged in that field, which is the biology being asked about.
  # The plain split: one entry per channel, tagged "file · channel". Used by
  # every analysis that does NOT divide by an area. Returns a single untagged
  # entry when there is nothing to split, so callers need no special case.
  channel_parts <- function(df, tag) {
    none <- list(list(df = df, tag = tag))
    if (!is.data.frame(df) || !nrow(df) || !"channel" %in% names(df)) return(none)
    ids <- sort(unique(df$channel[!is.na(df$channel)]))
    if (length(ids) < 2L) return(none)
    lapply(ids, function(ch) {
      sub <- df[which(df$channel == ch), , drop = FALSE]
      list(df = sub, tag = .facet_name(tag, as.character(sub$channel_label[1])))
    })
  }

  # The same split, plus the pooled parent, for the analyses that DO divide by an
  # area (see the note above).
  split_df_by_channel <- function(df, tag) {
    lapply(channel_parts(df, tag), function(pt)
      c(pt, list(win_df  = df,    # pooled file - the shared observation window
                 win_tag = tag))) # shared support cache key, so it is built once
  }

  # TRUE when the active dataset is a multiplex (drives the *_by_channel UIs).
  output$mplex_active <- reactive({
    df <- uploaded_df_rv()
    !is.null(df) && "mplex_label" %in% names(df)
  })
  outputOptions(output, "mplex_active", suspendWhenHidden = FALSE)

  # How many channels the loaded data currently resolves to. Drives every
  # conditionalPanel that only makes sense once there is more than one - a scope
  # control on single-channel data is a control with one legal answer.
  # Deliberately counts SOURCES rather than calling resolve_channels(): this runs
  # before Plot and must not depend on the clustering chain.
  output$n_channels <- reactive({
    df <- uploaded_df_rv()
    if (is.null(df)) return(1L)
    n_m <- if ("mplex_channel" %in% names(df))
      length(unique(df$mplex_channel[!is.na(df$mplex_channel)])) else 1L
    n_c <- if (isTRUE(input$enable_two_color)) 2L else 1L
    max(1L, n_m) * n_c
  })
  outputOptions(output, "n_channels", suspendWhenHidden = FALSE)

  # Scope picker, built from the channels that actually exist. There is no
  # "pooled" entry: pooling across channels is exactly what the channel column
  # exists to prevent, so offering it would be offering a wrong answer.
  output$channel_scope_ui <- renderUI({
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    ct <- NULL
    if (!is.null(df) && nrow(df))
      ct <- tryCatch(attr(attach_channel_colors(resolve_channels(df, input), input),
                          "channel_table"), error = function(e) NULL)
    if (is.null(ct) || nrow(ct) < 2L) return(NULL)
    ch <- stats::setNames(paste0("only:", ct$channel),
                          paste0("Only ", ct$label))
    selectInput("cluster_group_scope", "Channel scope",
                choices  = c("Within each channel (recommended)" = "within", ch),
                selected = input$cluster_group_scope %||% "within")
  })

  # Channel pair (A, B) for the pairwise cross-channel analyses, or NULL. When
  # set, the cross tools map channel A→group 1 / B→group 2 and reuse the entire
  # existing 2-color machinery (which honours a pre-assigned analysis_group).
  mplex_cross_pair <- reactive({
    df <- uploaded_df_rv()
    if (is.null(df) || !"mplex_label" %in% names(df)) return(NULL)
    a <- input$mplex_cross_a; b <- input$mplex_cross_b
    if (is.null(a) || is.null(b) || !nzchar(a) || !nzchar(b) || identical(a, b))
      return(NULL)
    c(a, b)
  })

  # Keep the channel selectors (cross-channel A/B + Batch-NN X/Y) populated.
  observeEvent(mplex_spec_rv(), {
    sp <- mplex_spec_rv()
    labs <- if (is.null(sp)) character(0) else vapply(sp, function(c) c$label, character(1))
    sel1 <- labs[1] %||% NULL
    sel2 <- if (length(labs) >= 2) labs[2] else NULL
    updateSelectInput(session, "mplex_cross_a", choices = labs, selected = sel1)
    updateSelectInput(session, "mplex_cross_b", choices = labs, selected = sel2)
    updateSelectInput(session, "nn_batch_chx",  choices = labs, selected = sel1)
    updateSelectInput(session, "nn_batch_chy",  choices = labs, selected = sel2)
  }, ignoreNULL = FALSE)


  uploaded_df <- reactive({
    df <- uploaded_df_rv()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data loaded yet."))
    df
  })

  # =========================================================================
  # commit_z_scale() - APPLY THE AXIAL CALIBRATION TO DATA ALREADY IN MEMORY
  # =========================================================================
  # GOAL
  #   Rescale every loaded dataset from the currently-applied axial factor to a
  #   new one. Rescaling by the ratio is exact and avoids re-reading files;
  #   z_scale_applied() tracks what is baked into the store so the ratio is
  #   always right.
  #
  # WHY IT IS NOT WIRED TO input$z_scale DIRECTLY
  #   Applying it replaces the active dataset, which fires the "new dataset"
  #   reset below. Committing only on Plot Data (or on an upload that records
  #   its own factor) keeps the change deliberate and lets the guard below tell
  #   a rescale apart from a genuine dataset change.
  #
  # ARGUMENTS  new : the factor to apply;  reason : shown in the notification
  # RETURNS    TRUE if anything was rescaled
  # =========================================================================
  z_rescale_guard <- reactiveVal(FALSE)

  commit_z_scale <- function(new, reason = "Axial calibration") {
    new <- suppressWarnings(as.numeric(new))
    old <- z_scale_applied()
    if (!length(new) || !is.finite(new) || new <= 0) return(invisible(FALSE))
    if (!is.finite(old) || old <= 0) { z_scale_applied(new); return(invisible(FALSE)) }
    if (isTRUE(all.equal(new, old))) return(invisible(FALSE))

    ratio <- new / old
    store <- smlm_file_store()

    # Nothing axial loaded: record the factor so a later 3D upload and any export
    # label are right, but do no work and stay silent. Rescaling zeros and
    # announcing it is the noise this guard exists to prevent.
    if (!any(vapply(store, has_axial_data, logical(1)))) {
      z_scale_applied(new)
      return(invisible(FALSE))
    }

    if (length(store)) {
      # Tell the reset observer this is a rescale of the SAME data, not a new
      # dataset - otherwise every metadata filter the user set is wiped.
      z_rescale_guard(TRUE)
      store <- lapply(store, function(d) {
        if (is.data.frame(d) && "z" %in% names(d)) d$z <- d$z * ratio
        d
      })
      smlm_file_store(store)
      cur <- isolate(input$active_file)
      if (!is.null(cur) && nzchar(cur) && !is.null(store[[cur]])) uploaded_df_rv(store[[cur]])
      # The z spread filter is expressed in nm, and nm along z just changed
      # meaning - carrying the old bound over would filter on a stale scale.
      sd_filt$z <- c(0, Inf)
    }
    z_scale_applied(new)
    showNotification(sprintf("%s set to %.3f - z distances rescaled by %.4f. Axial spread filter cleared.",
                             reason, new, ratio), type = "message", duration = 6)
    invisible(TRUE)
  }

  # Commit on Plot Data, so typing a value does not rescale mid-keystroke.
  observeEvent(input$plot_data, {
    commit_z_scale(input$z_scale)
  }, ignoreInit = TRUE)

  # Does the ACTIVE dataset actually carry axial data? Three states, because
  # "nothing loaded yet" must not look like "this data is 2D": an empty app
  # showing a greyed-out panel reads as broken rather than as not-applicable.
  smlm_z_state <- reactive({
    df <- uploaded_df_rv()
    if (is.null(df) || !nrow(df)) return("nodata")
    if (!"z" %in% names(df) || !any(is.finite(df$z)) ||
        stats::sd(df$z, na.rm = TRUE) <= 1e-6) return("planar")
    "3d"
  })

  # Dim (not hide) the panel for planar data. Hiding it would raise "where did
  # the setting go?"; dimmed-and-disabled shows it exists and says why it is off.
  observe({
    ok <- !identical(smlm_z_state(), "planar")
    shinyjs::toggleState("z_scale", condition = ok)
    shinyjs::toggleClass(id = "z_scale_group", class = "dimmed", condition = !ok)
  })

  # Every axis-bearing control follows the SAME planar test as the calibration
  # panel, so a dataset can never be simultaneously "2D" for one control and
  # "3D" for another. 3D mode and the XZ/YZ side views on planar data do not
  # error - they render a flat line - which is exactly the kind of quietly wrong
  # plot worth making unreachable.
  observe({
    planar <- identical(smlm_z_state(), "planar")

    shinyjs::toggleState(selector = "input[name='mode3d'][value='3D']",       condition = !planar)
    shinyjs::toggleState(selector = "input[name='view_plane_2d'][value='XZ']", condition = !planar)
    shinyjs::toggleState(selector = "input[name='view_plane_2d'][value='YZ']", condition = !planar)

    # Disabling a radio does not deselect it. Without this, loading planar data
    # while sitting in 3D/XZ would leave the view stuck in a mode whose control
    # is greyed out - unreachable and unfixable from the UI.
    if (planar) {
      if (!identical(isolate(input$mode3d) %||% "2D", "2D"))
        updateRadioButtons(session, "mode3d", selected = "2D")
      if (!identical(isolate(input$view_plane_2d) %||% "XY", "XY"))
        updateRadioButtons(session, "view_plane_2d", selected = "XY")
    }
  })

  output$z_axes_note <- renderUI({
    if (!identical(smlm_z_state(), "planar")) return(NULL)
    helpText(tags$b("2D dataset: "), "3D mode and the XZ / YZ side views are unavailable.")
  })

  # Shows whether the typed value is live yet.
  output$z_scale_status <- renderUI({
    if (identical(smlm_z_state(), "planar"))
      return(helpText(tags$b("Not applicable: "),
                      "the active dataset has no axial (z) data."))

    pend <- suppressWarnings(as.numeric(input$z_scale))
    appl <- z_scale_applied()
    if (!length(pend) || !is.finite(pend) || pend <= 0)
      return(helpText(sprintf("Applied: %.3f. Enter a positive factor.", appl)))
    if (isTRUE(all.equal(pend, appl)))
      helpText(sprintf("Applied: %.3f.", appl))
    else
      helpText(tags$b(sprintf("Pending: %.3f", pend)),
               sprintf(" (applied: %.3f) - press Plot Data to apply.", appl))
  })


  # --------- 2) RESET FILTERS ON NEW DATASET ----------
  observeEvent(uploaded_df(), {
    # An axial rescale replaces the active dataset object but it is the SAME
    # data - resetting every filter here would discard the user's whole setup
    # just because they changed the z factor. commit_z_scale() has already
    # cleared the one filter whose units actually changed (axial spread).
    if (isTRUE(z_rescale_guard())) { z_rescale_guard(FALSE); return() }

    # Respect the lock: when the user has checked "Lock filters", keep the
    # current filter values when switching between loaded files.
    # We still need to arm sd_skip_arming so the SD range observer doesn't fire
    # during the widget rebuild (sd_ranges_ui re-renders because sd_grey_tbl()
    # depends on uploaded_df()).  Without this, the brief NULL window while
    # widgets are being recreated causes the observer to write c(0,Inf) into
    # sd_filt, resetting the locked spread filter.
    if (isTRUE(isolate(input$lock_filters))) {
      sd_skip_arming(TRUE)
      later::later(function() sd_skip_arming(FALSE), delay = 0.8)
      return()
    }

    # Freeze meta_full_rng to current metadata widget values so the all.equal()
    # bypass in metadata_filtered_df passes during the browser round-trip window
    # between upload and the browser processing the new renderUI widgets.
    freeze_if_valid <- function(val) {
      if (!is.null(val) && length(val) == 2 && all(is.finite(as.numeric(val)))) val else NULL
    }
    cur_efo   <- freeze_if_valid(isolate(input$efo_range))
    cur_cfr   <- freeze_if_valid(isolate(input$cfr_range))
    cur_tim   <- freeze_if_valid(isolate(input$tim_range))
    cur_dcr   <- freeze_if_valid(isolate(input$dcr_range))
    cur_eco   <- freeze_if_valid(isolate(input$eco_range))
    cur_dwell <- freeze_if_valid(isolate(input$dwell_range))
    cur_fluo  <- freeze_if_valid(isolate(input$fluo_range))
    if (!is.null(cur_efo))   meta_full_rng$efo   <- cur_efo
    if (!is.null(cur_cfr))   meta_full_rng$cfr   <- cur_cfr
    if (!is.null(cur_tim))   meta_full_rng$tim   <- cur_tim
    if (!is.null(cur_dcr))   meta_full_rng$dcr   <- cur_dcr
    if (!is.null(cur_eco))   meta_full_rng$eco   <- cur_eco
    if (!is.null(cur_dwell)) meta_full_rng$dwell <- cur_dwell
    if (!is.null(cur_fluo))  meta_full_rng$fluo  <- cur_fluo

    # Reset SD filter synchronously - no JS round-trip needed.
    sd_filt$x <- c(0, Inf); sd_filt$y <- c(0, Inf); sd_filt$z <- c(0, Inf)
    # Signal tid_length_range_ui to reset to the new data's full range instead of
    # carrying forward the previous dataset's range.
    tid_reset_flag(TRUE)

    # Whenever a NEW dataset is loaded, reset filters to its full ranges
    shinyjs::reset("filters_form")

    vals <- filters_defaults_from_data(uploaded_df())

    session$onFlushed(function() {
      # Suppress SD range observer fires triggered by apply_filters updating the
      # SD widgets (only fires when the SD panel is already open).
      sd_skip_arming(TRUE)
      apply_filters(vals)
      meta_full_rng$efo   <- vals$efo_range
      meta_full_rng$cfr   <- vals$cfr_range
      meta_full_rng$tim   <- vals$tim_range
      meta_full_rng$dcr   <- vals$dcr_range
      meta_full_rng$eco   <- vals$eco_range
      meta_full_rng$dwell <- vals$dwell_range
      meta_full_rng$fluo  <- vals$fluo_range
      later::later(function() sd_skip_arming(FALSE), delay = 0.8)
    }, once = TRUE)

    # Also clear any loaded filter-settings JSON file input (optional)
    shinyjs::runjs("
    const el = document.getElementById('load_settings');
    if (el) { el.value=''; $(el).trigger('change'); }
  ")

  }, ignoreInit = TRUE)
  
  # --------- 3) COLLECT CURRENT UI SETTINGS ----------
  get_ui_settings <- reactive({
    req(uploaded_df_rv())
    list(
  
      efo_range        = input$efo_range,
      cfr_range        = input$cfr_range,
      tim_range        = input$tim_range,
      dcr_range        = input$dcr_range,
      eco_range        = input$eco_range,
      tid_length_range = input$tid_length_range,
      dwell_range      = input$dwell_range,
      fluo_range       = input$fluo_range,
      # SD filter - save from sd_filt (source of truth)
      # which is unreliable when the SD panel is closed (renderUI removes the widgets).
      sd_x_range       = sd_filt$x,
      sd_y_range       = sd_filt$y,
      sd_z_range       = sd_filt$z,
      sd_min_locs      = input$sd_min_locs,
      # Outlier trimming REMOVES LOCALIZATIONS, so it changes every per-trace
      # statistic downstream. Omitting it would make a reloaded file silently
      # analyse a different point set than the one the ranges were set on.
      trim_outliers    = isTRUE(input$trim_outliers),
      trim_k           = input$trim_k,
      make_positive    = input$make_positive,

      # ── Rogue-trace removal - THIS CHANGES THE DATA ────────────────────────
      # Dropped traces are removed before grouping, clustering, fits and surfaces,
      # so omitting these from the saved file meant reloading a settings file
      # silently reverted rogue removal to OFF and changed what was analysed.
      remove_rogues = isTRUE(input$remove_rogues),
      rogue_rule    = input$rogue_rule,
      rogue_factor  = input$rogue_factor,

      # ── Histogram + spread-panel display state ─────────────────────────────
      # apply_filters() already restores hist_variable/hist_bins/sd_* , but they
      # were never written out, so every load reset them to the data defaults.
      show_metadata_hist      = isTRUE(input$show_metadata_hist),
      hist_variable           = input$hist_variable,
      hist_bins               = input$hist_bins,
      hist_show_raw           = isTRUE(input$hist_show_raw),
      view_SD                 = isTRUE(input$view_SD),
      sd_bins                 = input$sd_bins,
      sd_hist_shared          = isTRUE(input$sd_hist_shared),
      sd_hist_range           = input$sd_hist_range,
      sd_overlay_all          = isTRUE(input$sd_overlay_all),
      sd_show_mean_sd         = isTRUE(input$sd_show_mean_sd),
      sd_show_overall_mean_sd = isTRUE(input$sd_show_overall_mean_sd),
      sd_use_q95              = isTRUE(input$sd_use_q95),
      target_sd_x             = input$target_sd_x,
      target_sd_y             = input$target_sd_y,
      target_sd_z             = input$target_sd_z,

      z_scale             = input$z_scale,

      enable_clustering   = input$enable_clustering,
      deblink_method      = input$deblink_method,
      eps_centroid        = input$eps_centroid,
      vote_frac           = input$vote_frac,
      vote_rho            = input$vote_rho,
      vote_link           = input$vote_link,
      centroid_k          = input$centroid_k,
      # part of the merge metric: eps only means the same thing alongside it
      deblink_anisotropic = input$deblink_anisotropic,
      # Channel handling. cluster_group_scope decides WHICH channels a run
      # covers; deblink_shared_spread decides whether the merge distance is
      # measured per channel or once over all of them. Both change what a
      # de-blinked molecule IS, so neither may be left to a default on reload.
      cluster_group_scope   = input$cluster_group_scope,
      deblink_shared_spread = isTRUE(input$deblink_shared_spread),
      # minPts is derived from tid_length_range and exclude_unclustered no longer
      # exists, so neither is saved - both would be stale on reload.
      
      enable_cluster2  = input$enable_cluster2,
      eps2             = input$eps2,
      minpts2          = input$minpts2,
      
      two_color = list(
        enable        = input$enable_two_color,
        use_fluo      = input$color_by_fluo,
        dcr_group1    = input$dcr_group1,
        dcr_group2    = input$dcr_group2,
        color_group1  = input$color_group1,
        color_group2  = input$color_group2
      ),
      
      # background overlay settings 
      overlay = list(
        show_confocal = input$show_confocal,
        show_bg       = input$show_bg,
        px_size   = input$px_size,
        img_x0    = input$img_x0,
        img_y0    = input$img_y0,
        img_invert_y = input$img_invert_y,
        img_flip_x   = input$img_flip_x,
        img_opacity  = input$img_opacity,
        img_rotation = input$img_rotation,
        img_lut      = input$img_lut,
        img_lut_custom = input$img_lut_custom,
        img_brightness = input$img_brightness,
        img_contrast   = input$img_contrast
      ),
      # heatmap settings
      heatmap = list(
        show_heatmap      = input$show_heatmap,
        heatmap_n         = input$heatmap_n,
        heatmap_opacity   = input$heatmap_opacity,
        heatmap_palette   = input$heatmap_palette,
        heatmap_bw_factor = input$heatmap_bw_factor,
        heatmap_bw_abs    = input$heatmap_bw_abs,
        heatmap_floor     = input$heatmap_floor,
        heatmap_use_raw   = input$heatmap_use_raw,
        heatmap_clip      = input$heatmap_clip,
        heatmap_gamma     = input$heatmap_gamma,
        heatmap_log       = input$heatmap_log,
        heatmap_showlegend= input$heatmap_showlegend
      ),

      # ── Spatial crop / trimmer - CHANGES THE DATA (the point set) ──────────
      # Both the 2D "Crop to view" box and the 3D box/sphere crop
      # remove localizations, so a reload without them analyses a wider region
      # than the ranges were set on. Saved so a round-trip keeps the same crop.
      crop = {
        c2d <- smlm_2d_crop_rv()
        list(
          crop2d_x       = if (!is.null(c2d)) c2d$x else NULL,
          crop2d_y       = if (!is.null(c2d)) c2d$y else NULL,
          enable_trimmer = isTRUE(input$enable_trimmer),
          trim_mode      = input$trim_mode,
          trim_x_range   = input$trim_x_range,
          trim_y_range   = input$trim_y_range,
          trim_z_range   = input$trim_z_range,
          trim_cx = input$trim_cx, trim_cy = input$trim_cy,
          trim_cz = input$trim_cz, trim_r  = input$trim_r,
          # the committed region is what actually crops (inputs above are the editor)
          trim_applied   = trim_applied_rv()
        )
      },

      # ── Fiducial-scan exclusion - CHANGES THE DATA when on ─────────────────
      # The thresholds are saved too so re-running the scan reproduces the flags.
      fiducial = list(
        fid_exclude   = isTRUE(input$fid_exclude),
        fid_efo_z     = input$fid_efo_z,
        fid_dur_x     = input$fid_dur_x,
        fid_min_locs  = input$fid_min_locs,
        fid_len_z     = input$fid_len_z,
        fid_recur_r   = input$fid_recur_r,
        fid_recur_min = input$fid_recur_min,
        fid_all_files = isTRUE(input$fid_all_files)
      ),

      # ── Analysis options that change results ───────────────────────────────
      analysis = list(
        use_pooled_dcr       = isTRUE(input$use_pooled_dcr),
        spatial_edge_correct = isTRUE(input$spatial_edge_correct),
        spatial_use_support  = isTRUE(input$spatial_use_support),
        dens_support_from    = input$dens_support_from,
        dens_norm_enable     = isTRUE(input$dens_norm_enable),
        dens_thin_enable     = isTRUE(input$dens_thin_enable)
      ),

      # ── Local density & hotspots ────────────────────────────────────
      # Every one of these changes a reported number, not just a picture: k and the
      # quadrat size set the scale the density is measured at, the footprint pair
      # sets the area every statistic divides by, and the eps calibration decides
      # what counts as a focus. A reload that silently reverted any of them would
      # reproduce the same plot from different data.
      hotspot = list(
        hot_unit        = input$hot_unit,
        hot_use_3d      = isTRUE(input$hot_use_3d),
        hot_k           = input$hot_k,
        hot_enh_thresh  = input$hot_enh_thresh,
        hot_quadrat     = input$hot_quadrat,
        hot_minpts      = input$hot_minpts,
        hot_target_frac = input$hot_target_frac,
        hot_nsim        = input$hot_nsim,
        hot_palette     = input$hot_palette,
        hot_show_foci   = isTRUE(input$hot_show_foci),
        hot_box_top     = input$hot_box_top,
        hot_metric      = input$hot_metric,
        hot_curve_mode  = input$hot_curve_mode
      )
    )
  })
  # --------- 4) DEFAULTS saved + APPLY functions  ----------
  filters_defaults_from_data <- function(df) {
    rng_of <- function(v) { v <- v[is.finite(v)]; if (!length(v)) c(NA,NA) else range(v, na.rm=TRUE) }
    tid_rng <- {
      if ("tid" %in% names(df)) {
        tl <- as.numeric(table(df$tid)); if (length(tl)) range(tl, na.rm=TRUE) else c(NA,NA)
      } else c(NA,NA)
    }
    # SD defaults (computed on uploaded data)
    # Lower bound is always 0 to match the widget initial value; the bypass check
    # in spread_filtered_df compares input values against meta_full_rng which is
    # also initialized with c(0, max), so they must be consistent.
    max_or1 <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return(1); mx <- max(v, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) 1 else mx }
    sda <- sd_by_tid(df, min_n = input$sd_min_locs)
    rx  <- if (nrow(sda)) c(0, max_or1(sda$sd_x)) else c(NA, NA)
    ry  <- if (nrow(sda)) c(0, max_or1(sda$sd_y)) else c(NA, NA)
    rz  <- if (nrow(sda) && "sd_z" %in% names(sda) && any(is.finite(sda$sd_z))) c(0, max_or1(sda$sd_z)) else c(NA, NA)
    
    list(
      hist_variable="efo", hist_bins=200,
      efo_range = if ("efo"%in%names(df)) rng_of(df$efo) else c(NA,NA),
      cfr_range = if ("cfr"%in%names(df)) rng_of(df$cfr) else c(NA,NA),
      tim_range = if ("tim"%in%names(df)) rng_of(df$tim) else c(NA,NA),
      dcr_range = if ("dcr"%in%names(df)) rng_of(df$dcr) else c(NA,NA),
      eco_range = if ("eco"%in%names(df)) rng_of(df$eco) else c(NA,NA),
      dwell_range = if ("dwell"%in%names(df)) rng_of(df$dwell) else c(NA,NA),
      fluo_range  = if ("fluo"%in%names(df)) rng_of(df$fluo) else c(NA,NA),
      tid_length_range = tid_rng,
      #SD summary controls
      target_sd_x=NA_real_, target_sd_y=NA_real_, target_sd_z=NA_real_,
      sd_bins = 100,
      sd_overlay_all = TRUE,
      sd_show_mean_sd = TRUE,
      sd_show_overall_mean_sd = FALSE,
      
      # --- SD filter ranges ---
      sd_x_range = rx, sd_y_range = ry, sd_z_range = rz
    )
  }
  overlay_defaults <- reactive({
    list(show_confocal=FALSE, show_bg=FALSE, px_size=10, img_x0=0, img_y0=0,
         img_invert_y=TRUE, img_flip_x=FALSE, img_opacity=0.6, img_rotation="0",
         img_lut="original", img_lut_custom="#00FF00", img_brightness=100,
         img_contrast=0, bg_image_name=NULL)
  })
  heatmap_defaults <- reactive({
    list(show_heatmap=FALSE, heatmap_n=200, heatmap_opacity=0.35, heatmap_palette="YlOrRd",
         heatmap_bw_factor=1, heatmap_bw_abs=NA_real_, heatmap_floor=0, heatmap_use_raw=FALSE,
         heatmap_clip=0.99, heatmap_gamma=1, heatmap_log=FALSE, heatmap_showlegend=FALSE)
  })
  cluster2_defaults <- function() {
    list(
      cluster2_method = "dbscan",
      cluster2_xi = 0.05,
      cluster2_prefilter = FALSE,
      cluster2_prefilter_r = 200,
      cluster2_prefilter_min = 3,
      cluster2_prune = FALSE,
      cluster2_prune_tol = 2,
      cluster2_core_only = FALSE,
      cluster2_min_members = 0,
      eps2 = 150,
      minpts2 = 8,
      
      # 3D wrappers
      show_cluster2_shape = FALSE,
      show_cluster2_isosurface = FALSE,
      cluster2_iso_selected_only = TRUE,
      cluster2_iso_id = 1,
      cluster2_iso_color = "#00CED1",
      cluster2_iso_unicolor = FALSE,
      cluster2_iso_by_channel = TRUE,
      cluster2_iso_grid_res = 22,
      cluster2_iso_pad = 0,
      cluster2_iso_manual_offset = FALSE,
      cluster2_iso_offset = 30,
      cluster2_iso_opacity = 0.25,
      cluster2_iso_manual_sigma = FALSE,
      cluster2_iso_smooth_rel = 1,
      cluster2_iso_sigma = 30,
      
      # 2D wrapper
      show_cluster2_shape_2D = FALSE,
      
      # metrics toggles
      show_cluster2_metrics = FALSE,
      
      # 2D metrics
      show_cluster2_pca_axes = FALSE,   # note: you reuse same id in 2D & 3D panels
      show_fit_distances_2d = FALSE,
      fit_color_2d = "#00BFFF",
      fit_line_width_2d = 3,
      
      # 3D metrics
      show_cluster2_labels = FALSE,
      show_cluster2_plane_angle = FALSE,
      show_cluster2_plane_distance = FALSE,
      show_bending_radius = FALSE,
      fit_method = "heuristic",
      fit_weighted = FALSE,
      mle_center_mode = "free",
      mle_center_lambda = 0.1,
      fit_type = "none",
      helix_turns = 3,
      show_fit_tube = FALSE,
      fit_tube_radius = 5,
      show_fit_distances = FALSE,
      fit_opacity = 0.4,
      
      # NN within cluster2
      show_cluster2_nn_lines = FALSE,
      nn_k_cluster2 = 3,
      distance_range_cluster2 = c(0, 1250)
    )
  }
  two_color_defaults <- function() {
    list(
      color_by_fluo = FALSE,
      dcr_group1 = c(0, 0.369),
      dcr_group2 = c(0.3691, 1),
      
      color_group1 = "#1f77a4",
      color_group2 = "#d62728",
      
      
      # CLUSTERING PARAMETERS ARE NOT HERE ANY MORE. There is one clustering
      # pipeline, so there is one set of de-blink / median-clustering controls,
      # and they live in the main settings file. This file now carries only what
      # DEFINES the channels and how they are drawn. Old files that still hold
      # the retired *_tc1 / *_tc2 keys are migrated on load (see below).
      show_tc1_alpha = FALSE,
      tc1_alpha = 120,
      tc1_alpha_opacity = 0.15,
      tc1_alpha_min_pts = 5
      # The Median Cluster wrap and isosurface are no longer duplicated here;
      # they are saved with the Median Clustering controls. Old files carrying
      # show_tc2_alpha / show_tc2_isosurface are migrated on load.
    )
  }
  
  apply_two_color_settings <- function(session, s) {
    getv <- function(name, default = NULL) if (!is.null(s[[name]])) s[[name]] else default
    
    # checkboxes
    cb_ids <- c(
      # master toggle first, so the restored 2-color parameters take effect
      "enable_two_color",
      "color_by_fluo", "use_pooled_dcr",
      "show_tc1_alpha"
    )
    for (id in cb_ids) if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))
    
    # numeric range inputs (shinyWidgets)
    if (!is.null(getv("dcr_group1")) && length(getv("dcr_group1")) == 2)
      shinyWidgets::updateNumericRangeInput(session, "dcr_group1", value = as.numeric(getv("dcr_group1")))
    if (!is.null(getv("dcr_group2")) && length(getv("dcr_group2")) == 2)
      shinyWidgets::updateNumericRangeInput(session, "dcr_group2", value = as.numeric(getv("dcr_group2")))
    
    # colors
    if (!is.null(getv("color_group1"))) colourpicker::updateColourInput(session, "color_group1", value = as.character(getv("color_group1")))
    if (!is.null(getv("color_group2"))) colourpicker::updateColourInput(session, "color_group2", value = as.character(getv("color_group2")))
    
    # NN
    # Retired 2-colour NN controls migrate onto the shared NN panel: there is one
    # NN section now, and enabling channels adapts it rather than adding a second.
    if (!is.null(getv("nn_k_two_color")) && is.null(s[["nn_k"]]))
      updateNumericInput(session, "nn_k", value = as.integer(getv("nn_k_two_color")))
    if (!is.null(getv("tc_label_font_size")) && is.null(s[["nn_label_font_size"]]))
      updateNumericInput(session, "nn_label_font_size", value = as.numeric(getv("tc_label_font_size")))
    if (!is.null(getv("tc_label_offset")) && is.null(s[["nn_label_offset"]]))
      updateNumericInput(session, "nn_label_offset", value = as.numeric(getv("tc_label_offset")))
    rng <- getv("distance_range_two_color")
    if (!is.null(rng) && length(rng) == 2 && is.null(s[["distance_range"]]))
      updateSliderInput(session, "distance_range", value = as.numeric(rng))
    # "between" / "within" map onto the shared three-valued channel scope.
    if (!is.null(getv("two_color_scope")) && is.null(s[["nn_mplex_scope"]]))
      updateRadioButtons(session, "nn_mplex_scope",
        selected = if (identical(as.character(getv("two_color_scope")), "between")) "across" else "within")
    if (!is.null(getv("draw_nn_two_color")) && isTRUE(getv("draw_nn_two_color")))
      updateCheckboxInput(session, "show_nn_lines", value = TRUE)
    
    # ── LEGACY MIGRATION: retired *_tc1 / *_tc2 keys ────────────────────────
    # Settings files written before the pipelines were unified carry a second,
    # 2-colour-only copy of every clustering parameter. There is one pipeline
    # now, so those values are applied to the single set of controls rather than
    # dropped - silently ignoring them would reload a saved analysis with default
    # radii, which is the failure mode this migration exists to prevent.
    #
    # Only fills a control the new file did not already set, so a settings file
    # containing both old and new keys keeps the new ones.
    tc_legacy <- list(
      eps_tc1                    = "eps_centroid",
      deblink_method_tc1         = "deblink_method",
      vote_frac_tc1              = "vote_frac",
      vote_rho_tc1               = "vote_rho",
      eps_tc2                    = "eps2",
      minpts_tc2                 = "minpts2",
      cluster2_method_tc2        = "cluster2_method",
      cluster2_xi_tc2            = "cluster2_xi",
      cluster2_prefilter_tc2     = "cluster2_prefilter",
      cluster2_prefilter_r_tc2   = "cluster2_prefilter_r",
      cluster2_prefilter_min_tc2 = "cluster2_prefilter_min",
      cluster2_prune_tc2         = "cluster2_prune",
      cluster2_prune_tol_tc2     = "cluster2_prune_tol",
      cluster2_core_only_tc2     = "cluster2_core_only",
      cluster2_min_members_tc2   = "cluster2_min_members",
      enable_two_color_cluster1  = "enable_clustering",
      enable_two_color_cluster2  = "enable_cluster2"
    )
    migrated <- character(0)
    for (old_id in names(tc_legacy)) {
      v <- getv(old_id)
      if (is.null(v)) next
      new_id <- tc_legacy[[old_id]]
      if (!is.null(s[[new_id]])) next          # new key wins
      if (is.logical(v) || old_id %in% c("enable_two_color_cluster1", "enable_two_color_cluster2",
                                         "cluster2_prefilter_tc2", "cluster2_prune_tc2",
                                         "cluster2_core_only_tc2")) {
        updateCheckboxInput(session, new_id, value = isTRUE(v))
      } else if (is.character(v)) {
        updateRadioButtons(session, new_id, selected = as.character(v))
      } else {
        updateNumericInput(session, new_id, value = as.numeric(v))
      }
      migrated <- c(migrated, old_id)
    }
    if (length(migrated))
      showNotification(sprintf(
        "Loaded an older 2-colour settings file: %d retired clustering setting(s) were applied to the unified de-blink / Median Clustering controls.",
        length(migrated)), type = "message", duration = 8)

    # 2-colour display settings (these are still 2-colour-specific)
    if (!is.null(getv("tc1_alpha"))) updateNumericInput(session, "tc1_alpha", value = as.numeric(getv("tc1_alpha")))
    if (!is.null(getv("tc1_alpha_opacity"))) updateSliderInput(session, "tc1_alpha_opacity", value = as.numeric(getv("tc1_alpha_opacity")))
    if (!is.null(getv("tc1_alpha_min_pts"))) updateNumericInput(session, "tc1_alpha_min_pts", value = as.integer(getv("tc1_alpha_min_pts")))
    
    # RETIRED: the two-colour Median Cluster wrap and isosurface. They now live
    # under Median Clustering and draw per channel. An older file that had them
    # switched on turns on the unified controls instead, carrying the isosurface
    # geometry across so the surface looks as it did.
    tc2_moved <- character(0)
    if (isTRUE(getv("show_tc2_alpha"))) {
      updateCheckboxInput(session, "show_cluster2_shape", value = TRUE)
      tc2_moved <- c(tc2_moved, "alpha-shape wrap")
    }
    if (isTRUE(getv("show_tc2_isosurface"))) {
      updateCheckboxInput(session, "show_cluster2_isosurface", value = TRUE)
      updateCheckboxInput(session, "cluster2_iso_by_channel", value = TRUE)
      for (pr in list(c("tc2_iso_grid_res",   "cluster2_iso_grid_res"),
                      c("tc2_iso_pad",        "cluster2_iso_pad"),
                      c("tc2_iso_smooth_rel", "cluster2_iso_smooth_rel"),
                      c("tc2_iso_offset",     "cluster2_iso_offset"),
                      c("tc2_iso_sigma",      "cluster2_iso_sigma"),
                      c("tc2_iso_opacity",    "cluster2_iso_opacity"))) {
        v <- suppressWarnings(as.numeric(getv(pr[1])))
        if (length(v) == 1L && is.finite(v)) updateSliderInput(session, pr[2], value = v)
      }
      for (pr in list(c("tc2_iso_manual_offset", "cluster2_iso_manual_offset"),
                      c("tc2_iso_manual_sigma",  "cluster2_iso_manual_sigma")))
        if (!is.null(getv(pr[1])))
          updateCheckboxInput(session, pr[2], value = isTRUE(getv(pr[1])))
      tc2_moved <- c(tc2_moved, "isosurface")
    }
    if (length(tc2_moved))
      showNotification(sprintf(
        "This settings file used the two-colour Median Cluster %s. That now lives under Median Clustering, drawn separately for each channel, and has been switched on there.",
        paste(tc2_moved, collapse = " and ")), type = "message", duration = 10)
  }



  apply_overlay_settings <- function(ov, session, current_bg_name=NULL) {  # (same as your 2D)
    if (!is.null(ov$show_confocal)) updateCheckboxInput(session,"show_confocal", value=ov$show_confocal)
    if (!is.null(ov$show_bg))       updateCheckboxInput(session,"show_bg",       value=ov$show_bg)
    if (!is.null(ov$px_size))      updateNumericInput(session,"px_size", value=ov$px_size)
    if (!is.null(ov$img_x0))       updateNumericInput(session,"img_x0",  value=ov$img_x0)
    if (!is.null(ov$img_y0))       updateNumericInput(session,"img_y0",  value=ov$img_y0)
    if (!is.null(ov$img_invert_y)) updateCheckboxInput(session,"img_invert_y", value=ov$img_invert_y)
    if (!is.null(ov$img_flip_x))   updateCheckboxInput(session,"img_flip_x",   value=ov$img_flip_x)
    if (!is.null(ov$img_opacity))  updateSliderInput(session,"img_opacity",    value=ov$img_opacity)
    if (!is.null(ov$img_rotation)) updateSelectInput(session,"img_rotation",   selected=as.character(ov$img_rotation))
    if (!is.null(ov$img_lut))      updateSelectInput(session,"img_lut",        selected=as.character(ov$img_lut))
    if (!is.null(ov$img_lut_custom)) colourpicker::updateColourInput(session,"img_lut_custom", value=as.character(ov$img_lut_custom))
    if (!is.null(ov$img_brightness)) updateSliderInput(session,"img_brightness", value=ov$img_brightness)
    if (!is.null(ov$img_contrast))   updateSliderInput(session,"img_contrast",   value=ov$img_contrast)
    if (!is.null(ov$bg_image_name)) {
      cur <- if (!is.null(input$bg_image)) input$bg_image$name else NULL
      if (is.null(cur) || !identical(cur, ov$bg_image_name)) {
        showNotification(paste0("Overlay loaded. Re-upload '", ov$bg_image_name, "' to reproduce alignment."),
                         type="message", duration=6)
      }
    }
  }
  apply_filters <- function(vals) {
    # jsonlite::fromJSON may parse JSON arrays as lists rather than numeric vectors;
    # coerce all range values defensively so updateNumericRangeInput never receives a list.
    rng <- function(v) if (is.null(v)) NULL else as.numeric(unlist(v))
    updateSelectInput(session,"hist_variable", selected=vals$hist_variable)
    updateSliderInput(session,"hist_bins", value=vals$hist_bins)
    shinyWidgets::updateNumericRangeInput(session,"efo_range",  value=rng(vals$efo_range))
    shinyWidgets::updateNumericRangeInput(session,"cfr_range",  value=rng(vals$cfr_range))
    shinyWidgets::updateNumericRangeInput(session,"tim_range",  value=rng(vals$tim_range))
    shinyWidgets::updateNumericRangeInput(session,"dcr_range",  value=rng(vals$dcr_range))
    shinyWidgets::updateNumericRangeInput(session,"eco_range",  value=rng(vals$eco_range))
    shinyWidgets::updateNumericRangeInput(session,"dwell_range",value=rng(vals$dwell_range))
    shinyWidgets::updateNumericRangeInput(session,"fluo_range", value=rng(vals$fluo_range))
    shinyWidgets::updateNumericRangeInput(session,"tid_length_range", value=rng(vals$tid_length_range))
    updateNumericInput(session,"target_sd_x", value=vals$target_sd_x)
    updateNumericInput(session,"target_sd_y", value=vals$target_sd_y)
    updateNumericInput(session,"target_sd_z", value=vals$target_sd_z)
    updateSliderInput(session, "sd_bins", value = vals$sd_bins)
    updateCheckboxInput(session, "sd_hist_shared", value = isTRUE(vals$sd_hist_shared))
    if (!is.null(vals$sd_hist_range))
      shinyWidgets::updateNumericRangeInput(session, "sd_hist_range",
                                            value = rng(vals$sd_hist_range))
    updateCheckboxInput(session, "sd_overlay_all", value = vals$sd_overlay_all)
    updateCheckboxInput(session, "sd_show_mean_sd", value = vals$sd_show_mean_sd %||% TRUE)
    updateCheckboxInput(session, "sd_show_overall_mean_sd", value = vals$sd_show_overall_mean_sd %||% FALSE)
    if (!is.null(vals$sd_x_range)) shinyWidgets::updateNumericRangeInput(session, "sd_x_range", value = rng(vals$sd_x_range))
    if (!is.null(vals$sd_y_range)) shinyWidgets::updateNumericRangeInput(session, "sd_y_range", value = rng(vals$sd_y_range))
    if (!is.null(vals$sd_z_range)) shinyWidgets::updateNumericRangeInput(session, "sd_z_range", value = rng(vals$sd_z_range))

    # ── Rogue removal + panel display state ───────────────────────────────────
    # Only restore when the key is actually present, so an OLD settings file (saved
    # before these were persisted) leaves the current UI alone instead of silently
    # switching rogue removal off.
    if (!is.null(vals$remove_rogues)) updateCheckboxInput(session, "remove_rogues", value = isTRUE(vals$remove_rogues))
    if (!is.null(vals$rogue_rule))    updateRadioButtons(session, "rogue_rule", selected = as.character(vals$rogue_rule))
    if (!is.null(vals$rogue_factor))  updateNumericInput(session, "rogue_factor", value = as.numeric(vals$rogue_factor))
    if (!is.null(vals$show_metadata_hist)) updateCheckboxInput(session, "show_metadata_hist", value = isTRUE(vals$show_metadata_hist))
    if (!is.null(vals$hist_show_raw))      updateCheckboxInput(session, "hist_show_raw", value = isTRUE(vals$hist_show_raw))
    if (!is.null(vals$sd_use_q95))         updateCheckboxInput(session, "sd_use_q95", value = isTRUE(vals$sd_use_q95))
  }
  apply_overlay  <- function(vals) { apply_overlay_settings(vals, session) }
  apply_heatmap  <- function(vals) {
    updateCheckboxInput(session,"show_heatmap", value=vals$show_heatmap)
    updateNumericInput(session,"heatmap_n", value=vals$heatmap_n)
    updateSliderInput(session,"heatmap_opacity", value=vals$heatmap_opacity)
    updateSelectInput(session,"heatmap_palette", selected=vals$heatmap_palette)
    updateSliderInput(session,"heatmap_bw_factor", value=vals$heatmap_bw_factor)
    updateNumericInput(session,"heatmap_bw_abs", value=vals$heatmap_bw_abs)
    updateSliderInput(session,"heatmap_floor", value=vals$heatmap_floor)
    updateCheckboxInput(session,"heatmap_use_raw", value=vals$heatmap_use_raw)
    updateSliderInput(session,"heatmap_clip", value=vals$heatmap_clip)
    updateSliderInput(session,"heatmap_gamma", value=vals$heatmap_gamma)
    updateCheckboxInput(session,"heatmap_log", value=vals$heatmap_log)
    updateCheckboxInput(session,"heatmap_showlegend", value=vals$heatmap_showlegend)
  }
  apply_cluster2_settings <- function(session, s) {
    getv <- function(name, default = NULL) {
      if (!is.null(s[[name]])) s[[name]] else default
    }
    
    # method + numeric inputs
    if (!is.null(getv("cluster2_method"))) updateRadioButtons(session, "cluster2_method", selected = as.character(getv("cluster2_method")))
    if (!is.null(getv("cluster2_xi"))) updateNumericInput(session, "cluster2_xi", value = as.numeric(getv("cluster2_xi")))
    if (!is.null(getv("cluster2_prefilter"))) updateCheckboxInput(session, "cluster2_prefilter", value = isTRUE(getv("cluster2_prefilter")))
    if (!is.null(getv("cluster2_prefilter_r"))) updateNumericInput(session, "cluster2_prefilter_r", value = as.numeric(getv("cluster2_prefilter_r")))
    if (!is.null(getv("cluster2_prefilter_min"))) updateNumericInput(session, "cluster2_prefilter_min", value = as.numeric(getv("cluster2_prefilter_min")))
    if (!is.null(getv("cluster2_prune"))) updateCheckboxInput(session, "cluster2_prune", value = isTRUE(getv("cluster2_prune")))
    if (!is.null(getv("cluster2_prune_tol"))) updateNumericInput(session, "cluster2_prune_tol", value = as.numeric(getv("cluster2_prune_tol")))
    if (!is.null(getv("cluster2_core_only"))) updateCheckboxInput(session, "cluster2_core_only", value = isTRUE(getv("cluster2_core_only")))
    if (!is.null(getv("cluster2_min_members"))) updateNumericInput(session, "cluster2_min_members", value = as.numeric(getv("cluster2_min_members")))
    if (!is.null(getv("eps2")))    updateNumericInput(session, "eps2", value = as.numeric(getv("eps2")))
    if (!is.null(getv("minpts2"))) updateNumericInput(session, "minpts2", value = as.numeric(getv("minpts2")))

    # ── PRIMARY (de-blinking) settings ──
    # These were written to the settings file but never restored, so a reloaded
    # analysis silently ran with the default radius. minpts / exclude_unclustered
    # are deliberately absent: minPts now derives from tid_length_range, and
    # exclude_unclustered no longer exists (de-blinking never removes molecules).
    for (nmv in c("vote_frac", "vote_rho"))
      if (!is.null(getv(nmv))) updateNumericInput(session, nmv, value = as.numeric(getv(nmv)))
    if (!is.null(getv("vote_link")))
      updateRadioButtons(session, "vote_link", selected = as.character(getv("vote_link")))
    if (!is.null(getv("eps_centroid")))
      updateNumericInput(session, "eps_centroid", value = as.numeric(getv("eps_centroid")))
    if (!is.null(getv("centroid_k")))
      updateNumericInput(session, "centroid_k", value = as.numeric(getv("centroid_k")))
    # Restored before enable_clustering, so the criterion is in place when it fires.
    if (!is.null(getv("deblink_method")))
      updateRadioButtons(session, "deblink_method", selected = as.character(getv("deblink_method")))
    if (!is.null(getv("deblink_anisotropic")))
      updateCheckboxInput(session, "deblink_anisotropic", value = isTRUE(getv("deblink_anisotropic")))
    if (!is.null(getv("cluster_group_scope")))
      updateSelectInput(session, "cluster_group_scope", selected = as.character(getv("cluster_group_scope")))
    if (!is.null(getv("deblink_shared_spread")))
      updateCheckboxInput(session, "deblink_shared_spread", value = isTRUE(getv("deblink_shared_spread")))
    # Restored last so the parameters above are in place before it takes effect.
    if (!is.null(getv("enable_clustering")))
      updateCheckboxInput(session, "enable_clustering", value = isTRUE(getv("enable_clustering")))

    # Axial calibration: rescales every z distance, so it must survive a reload.
    if (!is.null(getv("z_scale"))) {
      zs <- as.numeric(getv("z_scale"))
      updateNumericInput(session, "z_scale", value = zs)
      # Loading settings is itself an explicit action, so commit rather than
      # leaving the restored calibration pending behind a Plot Data press.
      commit_z_scale(zs, reason = "Axial calibration from settings")
    }


    # checkboxes
    cb_ids <- c(
      # master toggle first, so the restored parameters actually take effect
      "enable_cluster2", "show_cluster2_knn",
      "show_cluster2_shape", "show_cluster2_isosurface", "cluster2_iso_selected_only",
      "cluster2_iso_by_channel",
      "show_cluster2_shape_2D", "show_cluster2_metrics",
      "show_cluster2_pca_axes",
      "show_cluster2_labels", "show_cluster2_plane_angle", "show_cluster2_plane_distance",
      "show_bending_radius", "fit_weighted",
      "show_fit_distances", "show_fit_tube",
      "show_cluster2_nn_lines",
      "compute_rmsd_metrics", "membrane_protein_options", "show_c2_nn_labels"
    )
    for (id in cb_ids) {
      if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))
    }
    
    # selects / radios
    if (!is.null(getv("fit_method")))     updateRadioButtons(session, "fit_method", selected = as.character(getv("fit_method")))
    if (!is.null(getv("mle_center_mode")))updateRadioButtons(session, "mle_center_mode", selected = as.character(getv("mle_center_mode")))
    if (!is.null(getv("fit_type")))       updateSelectInput(session, "fit_type", selected = as.character(getv("fit_type")))
    # sliders / numeric misc
    if (!is.null(getv("mle_center_lambda"))) updateSliderInput(session, "mle_center_lambda", value = as.numeric(getv("mle_center_lambda")))
    if (!is.null(getv("helix_turns")))       updateSliderInput(session, "helix_turns", value = as.numeric(getv("helix_turns")))
    if (!is.null(getv("fit_tube_radius")))   updateNumericInput(session, "fit_tube_radius", value = as.numeric(getv("fit_tube_radius")))
    if (!is.null(getv("fit_opacity")))       updateNumericInput(session, "fit_opacity", value = as.numeric(getv("fit_opacity")))
    
    # isosurface controls
    if (!is.null(getv("cluster2_iso_id")))       updateNumericInput(session, "cluster2_iso_id", value = as.integer(getv("cluster2_iso_id")))
    if (!is.null(getv("cluster2_iso_grid_res"))) updateSliderInput(session, "cluster2_iso_grid_res", value = as.numeric(getv("cluster2_iso_grid_res")))
    if (!is.null(getv("cluster2_iso_pad")))      updateSliderInput(session, "cluster2_iso_pad", value = as.numeric(getv("cluster2_iso_pad")))
    if (!is.null(getv("cluster2_iso_offset")))   updateSliderInput(session, "cluster2_iso_offset", value = as.numeric(getv("cluster2_iso_offset")))
    if (!is.null(getv("cluster2_iso_manual_offset"))) updateCheckboxInput(session, "cluster2_iso_manual_offset", value = isTRUE(getv("cluster2_iso_manual_offset")))
    if (!is.null(getv("cluster2_iso_opacity")))  updateSliderInput(session, "cluster2_iso_opacity", value = as.numeric(getv("cluster2_iso_opacity")))
    if (!is.null(getv("cluster2_iso_sigma")))    updateSliderInput(session, "cluster2_iso_sigma", value = as.numeric(getv("cluster2_iso_sigma")))
    if (!is.null(getv("cluster2_iso_manual_sigma"))) updateCheckboxInput(session, "cluster2_iso_manual_sigma", value = isTRUE(getv("cluster2_iso_manual_sigma")))
    if (!is.null(getv("cluster2_iso_smooth_rel"))) updateSliderInput(session, "cluster2_iso_smooth_rel", value = as.numeric(getv("cluster2_iso_smooth_rel")))
    if (!is.null(getv("cluster2_iso_color")))    colourpicker::updateColourInput(session, "cluster2_iso_color", value = as.character(getv("cluster2_iso_color")))
    if (!is.null(getv("cluster2_iso_unicolor"))) updateCheckboxInput(session, "cluster2_iso_unicolor", value = isTRUE(getv("cluster2_iso_unicolor")))
    if (!is.null(getv("cluster2_iso_by_channel"))) updateCheckboxInput(session, "cluster2_iso_by_channel", value = isTRUE(getv("cluster2_iso_by_channel")))
    
    # color
    if (!is.null(getv("fit_color_2d"))) colourpicker::updateColourInput(session, "fit_color_2d", value = as.character(getv("fit_color_2d")))
    
    # fit/metric options added after this block was last updated
    if (!is.null(getv("arc_span_deg")))      updateSliderInput(session, "arc_span_deg", value = as.numeric(getv("arc_span_deg")))
    if (!is.null(getv("max_dist_lines_nm"))) updateNumericInput(session, "max_dist_lines_nm", value = as.numeric(getv("max_dist_lines_nm")))

    # NN controls
    if (!is.null(getv("nn_k_cluster2"))) updateNumericInput(session, "nn_k_cluster2", value = as.integer(getv("nn_k_cluster2")))
    if (!is.null(getv("c2_label_font_size"))) updateNumericInput(session, "c2_label_font_size", value = as.numeric(getv("c2_label_font_size")))
    if (!is.null(getv("c2_label_offset")))    updateNumericInput(session, "c2_label_offset", value = as.numeric(getv("c2_label_offset")))
    rng <- getv("distance_range_cluster2")
    if (!is.null(rng) && length(rng) == 2) {
      updateSliderInput(session, "distance_range_cluster2", value = as.numeric(rng))
    }
  }
  
  
  # --------- 5) RANGE UI ( ----------
  # Store committed default ranges for metadata + spread filters (used to avoid "filtering on load").
  # Updated in two steps per upload: (1) frozen to current inputs at the top of the upload observer
  # so the all.equal bypass passes during the browser round-trip window, then (2) updated to the
  # new data's ranges inside onFlushed after apply_filters runs.
  # Bypass check for metadata filters: NULL = no bypass; finite c(lo,hi) = committed
  # full-range value that the all.equal check in metadata_filtered_df compares against.
  meta_full_rng <- reactiveValues(
    efo = NULL, cfr = NULL, tim = NULL, dcr = NULL, eco = NULL, dwell = NULL, fluo = NULL
  )
  # Single source of truth for the SD (spread) filter.
  # c(0, Inf) = "no filter" - passes all traces regardless of spread.
  # Set to Inf on upload/reset; set to a finite upper bound only when the user
  # explicitly narrows the slider. Using Inf avoids all floating-point edge cases
  # that arise when comparing widget values against data maxima.
  sd_filt <- reactiveValues(x = c(0, Inf), y = c(0, Inf), z = c(0, Inf))
  # When TRUE, suppress the SD range observer (widget initialisation fires that
  # should not count as user-driven narrowing).
  sd_skip_arming <- reactiveVal(FALSE)
  # When TRUE, tid_length_range_ui resets to [min, max] on next render instead
  # of preserving the current value.  Set on new upload; cleared after one render.
  tid_reset_flag <- reactiveVal(FALSE)
  
  
  # ==========================================================================
  # FILTER PIPELINE - Stage 1 of 3: per-LOCALIZATION quality
  # ==========================================================================
  # GOAL: drop individual localizations that fail the microscope's own quality
  # metrics, before any trace-level judgement is made.
  #
  # These are the per-point columns MINFLUX reports (see the canonical schema in
  # helpers_shared.R):
  #     efo  photon rate      -> low efo = few photons = poor precision
  #     cfr  centre-freq ratio-> molecule poorly centred in the doughnut
  #     eco  photon count     dcr  channel ratio (the 2-colour axis)
  #     tim  acquisition time -> lets you cut the drifty tail of a long run
  #
  # KEY DESIGN CHOICE - `rng_filter` only applies a filter when the user's range
  # actually DIFFERS from the dataset's full range. A slider sitting at its
  # default must be a true no-op: comparing floating-point widget values against
  # data maxima otherwise drops boundary points for no reason.
  #
  # INPUT : uploaded_df()   OUTPUT: same columns, fewer rows
  # ==========================================================================
  metadata_filtered_df <- reactive({
    # Drop Fiducial-scan-flagged traces first (only when the exclusion is on;
    # reversible - see qc_exclude_from). Active file = input$active_file.
    df <- qc_exclude_from(uploaded_df(), input$active_file)

    # Apply range filter ONLY if user range differs from the dataset's full range
    rng_filter <- function(dat, nm, rng_val, full_rng) {
      if (!nm %in% names(dat)) return(dat)
      if (is.null(rng_val) || length(rng_val) != 2 || any(!is.finite(as.numeric(rng_val)))) return(dat)
      
      v <- suppressWarnings(as.numeric(dat[[nm]]))
      if (!length(v)) return(dat)
      
      # Prefer stored full range; fall back to computed full finite range if needed
      if (is.null(full_rng) || length(full_rng) != 2 || any(!is.finite(as.numeric(full_rng)))) {
        vv <- v[is.finite(v)]
        if (!length(vv)) return(dat)
        full_rng <- range(vv, na.rm = TRUE)
      }
      
      # If unchanged from defaults, treat as "no filter"
      if (isTRUE(all.equal(as.numeric(rng_val), as.numeric(full_rng), tolerance = 0))) return(dat)
      
      # When filtering is active: keep NA rows (so they don't vanish unexpectedly)
      dplyr::filter(dat, is.na(.data[[nm]]) | dplyr::between(.data[[nm]], rng_val[1], rng_val[2]))
    }
    
    df %>%
      rng_filter("efo",   input$efo_range,   meta_full_rng$efo) %>%
      rng_filter("cfr",   input$cfr_range,   meta_full_rng$cfr) %>%
      rng_filter("tim",   input$tim_range,   meta_full_rng$tim) %>%
      rng_filter("dcr",   input$dcr_range,   meta_full_rng$dcr) %>%
      rng_filter("eco",   input$eco_range,   meta_full_rng$eco) %>%
      rng_filter("dwell", input$dwell_range, meta_full_rng$dwell) %>%
      rng_filter("fluo",  input$fluo_range,  meta_full_rng$fluo)
  })

  # ===========================================================================
  #  SMLM FILTER PIPELINE  (read top-to-bottom; each stage feeds the next)
  #    Stage 1  metadata  metadata_filtered_df   (defined just above)
  #    Stage 2  spread    spread_filtered_df
  #    Stage 3  length    filtered_df            -> feeds plotting / clustering
  #  The Trace-ID-Length slider UI (tid_length_range_ui) is rendered right after
  #  Stage 3 because it derives its range from the spread-filtered data.
  # ===========================================================================

  # ==========================================================================
  # FILTER PIPELINE - Stage 2 of 3: per-TRACE spread (the SD filter)
  # ==========================================================================
  # GOAL: remove whole traces whose localizations are too SCATTERED to be one
  # stationary molecule.
  #
  # THE REASONING. In SMLM a `tid` should be one fluorophore sitting still while
  # blinking, so its localizations should scatter only by the localization
  # precision - a few nm. A trace with a much larger spread is not a molecule:
  # it is usually stage drift during a long blink, or two nearby molecules whose
  # localizations were merged into one trace ID. Either way its median position
  # is meaningless, and the median is exactly what the clustering and every
  # spatial statistic downstream will use.
  #
  # So the filter is per-axis SD of each trace's localizations, computed AFTER
  # Stage 1 (a trace should be judged on its surviving good points, not on
  # points already rejected as low quality).
  #
  # `sd_filt` holds the committed bounds. An upper bound of Inf means "no
  # filter" - an explicit sentinel, chosen so no floating-point comparison
  # against a data maximum is needed to express "keep everything".
  #
  # INPUT : metadata_filtered_df()   OUTPUT: whole traces removed
  # ==========================================================================
  # ==========================================================================
  # FILTER PIPELINE - Stage 1b of 3: outlier LOCALIZATIONS inside a trace
  # ==========================================================================
  # GOAL: drop the individual localizations that do not belong to their trace,
  # BEFORE any per-trace statistic is computed from them.
  #
  # WHY IT SITS HERE. Every later stage measures a trace: its spread (Stage 2),
  # its length (Stage 3), its median position (clustering). All of those are
  # computed from the surviving points, so a stray localization has to be gone
  # before them or it corrupts each one in turn. Placing it after the metadata
  # filters is deliberate for the same reason those come first - a point already
  # rejected on quality should not get a vote on whether its neighbours are
  # outliers.
  #
  # This is the ONLY place localizations are removed for spatial reasons.
  # Trace-level judgement stays in Stage 2 and the rogue rules; the two are kept
  # apart so that trimming can never quietly rescue a trace that is not one
  # molecule - the trim FRACTION is reported precisely so that case stays visible.
  #
  # OFF by default: opening the spread panel shows the untrimmed data first, and
  # resetting the filters returns to it.
  #
  # INPUT : metadata_filtered_df()   OUTPUT: outlier localizations removed
  # ==========================================================================
  trim_flags <- reactive({
    df <- metadata_filtered_df()
    if (!isTRUE(input$trim_outliers) || !is.data.frame(df) || !nrow(df))
      return(NULL)
    tryCatch(
      trace_point_outliers(df, k = input$trim_k %||% 4,
                           keep_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT),
      error = function(e) NULL)
  })

  trimmed_df <- reactive({
    df <- metadata_filtered_df()
    fl <- trim_flags()
    if (is.null(fl) || !is.data.frame(df) || length(fl) != nrow(df)) return(df)
    df[!fl, , drop = FALSE]
  })

  # Measured, not theoretical - the Gaussian rejection rate for a given k is
  # wrong by about a factor of ten on real MINFLUX residuals.
  output$trim_summary_ui <- renderUI({
    fl <- trim_flags()
    df <- metadata_filtered_df()
    if (is.null(fl) || !is.data.frame(df) || length(fl) != nrow(df))
      return(helpText("No data to trim yet."))
    s <- trim_summary(df, fl)

    # Trimming shortens traces, so a trace can drop below the Trace-ID-Length
    # minimum and be removed by Stage 3 even though it passed before. That is
    # the correct reading - a trace with 4 valid localizations IS a 4-point
    # trace - but it must not happen silently, so it is counted here rather
    # than prevented. Raising the trim guard to the length minimum instead
    # would make WHICH POINTS ARE OUTLIERS depend on how long traces are
    # required to be, and would keep known-bad localizations purely to satisfy
    # a count - corrupting the very median and spread the length filter exists
    # to protect.
    lo <- suppressWarnings(as.numeric((input$tid_length_range %||% c(0, Inf))[1]))
    demoted <- NA_integer_
    if (is.finite(lo) && lo > 0) {
      nb <- table(as.character(df$tid))
      na <- table(as.character(df$tid[!fl]))
      na <- na[match(names(nb), names(na))]; na[is.na(na)] <- 0
      demoted <- sum(nb >= lo & na < lo)
    }

    helpText(HTML(paste0(sprintf(
      "Rejected <b>%s</b> of %s localizations (<b>%.2f%%</b>) in %s of %s traces (%.0f%%).<br/>%s",
      format(s$n_flag, big.mark = ","), format(s$n_pts, big.mark = ","), s$pct,
      format(s$n_tr_touched, big.mark = ","), format(s$n_tr, big.mark = ","), s$pct_tr,
      if (s$n_tr_heavy > 0)
        sprintf(paste0("<b>%s trace(s) lost ≥ 30%% of their localizations</b> - that is ",
                       "not a clean trace with a few strays, and is worth inspecting before ",
                       "trusting its position."), format(s$n_tr_heavy, big.mark = ","))
      else "No trace lost 30% or more of its localizations."),
      if (!is.na(demoted) && demoted > 0)
        sprintf(paste0("<br/><b>%s trace(s) fell below the trace-length minimum of %g</b> ",
                       "because of trimming and are now removed by that filter."),
                format(demoted, big.mark = ","), lo)
      else ""
    )))
  })

  # Stage 2: spread filter - per-trace SD applied after metadata row filters
  spread_filtered_df <- reactive({
    df <- trimmed_df()
    if (!is.data.frame(df) || !nrow(df) || !"tid" %in% names(df)) return(df)

    # sd_filt holds the committed filter bounds set by the user.
    # c(0, Inf) upper bound = "no filter" - all traces pass regardless of spread.
    # This avoids every floating-point comparison edge case that arose from
    # matching widget values against data maxima computed from different subsets.
    rx <- c(sd_filt$x[1], sd_filt$x[2])
    ry <- c(sd_filt$y[1], sd_filt$y[2])
    rz <- c(sd_filt$z[1], sd_filt$z[2])

    # "No filter on this axis" needs BOTH ends open, not just the upper one. An
    # upper bound of Inf with a lower bound still set means "at least this much,
    # no ceiling" - treating Inf alone as no-filter would silently discard the
    # lower bound the user typed.
    open <- function(r) is.infinite(r[2]) && r[1] <= 0
    if (open(rx) && open(ry) && open(rz)) return(df)

    # per-tid SDs computed from CURRENT df (after metadata filters)
    sda <- sd_by_tid(df, min_n = input$sd_min_locs)
    if (!is.data.frame(sda) || !nrow(sda)) return(df[0, , drop = FALSE])

    if (!"sd_z" %in% names(sda)) sda$sd_z <- NA_real_

    # A trace whose spread is undefined on an axis (too few localizations) is
    # never judged on that axis - unchanged behaviour.
    in_rng <- function(v, r) is.na(v) | (v >= r[1] & (is.infinite(r[2]) | v <= r[2]))
    pass_x <- open(rx) | in_rng(sda$sd_x, rx)
    pass_y <- open(ry) | in_rng(sda$sd_y, ry)
    pass_z <- open(rz) | in_rng(sda$sd_z, rz)

    keep_tids <- sda$tid[ pass_x & pass_y & pass_z ]

    df %>% dplyr::filter(tid %in% keep_tids)
  })
  # ==========================================================================
  # FILTER PIPELINE - Stage 3 of 3: per-TRACE localization count
  # ==========================================================================
  # GOAL: keep only traces with a sensible number of localizations. This is the
  # last stage; its output is what the analyses consume.
  #
  # BOTH ENDS OF THE RANGE ARE MEANINGFUL:
  #   too FEW  localizations -> the trace's median position is poorly determined,
  #                             and it may be a spurious detection rather than a
  #                             real fluorophore.
  #   too MANY localizations -> often not a single molecule at all: a persistent
  #                             bright contaminant, or several molecules sharing
  #                             one trace ID. Such a trace would also dominate
  #                             any per-localization statistic by sheer count.
  #
  # Counting happens AFTER Stages 1-2 for the same reason as before: the count
  # that matters is of localizations that actually passed quality control.
  #
  # INPUT : spread_filtered_df()   OUTPUT: the analysis-ready trace set
  # ==========================================================================
  # Stage 3: trace ID length filter - applied after spread filter
  filtered_df <- reactive({
    df <- spread_filtered_df()

    # tid length filter
    if ("tid" %in% names(df) && !is.null(input$tid_length_range)) {
      tl <- table(df$tid)
      keep <- names(tl)[tl >= input$tid_length_range[1] & tl <= input$tid_length_range[2]]
      df <- df %>% dplyr::filter(tid %in% keep)
    }
    # Make positive
    if (isTRUE(input$make_positive)) {
      min_x <- suppressWarnings(min(df$x, na.rm = TRUE))
      min_y <- suppressWarnings(min(df$y, na.rm = TRUE))
      if (is.finite(min_x) && min_x < 0) df$x <- df$x - min_x
      if (is.finite(min_y) && min_y < 0) df$y <- df$y - min_y
      if ("z" %in% names(df)) {
        min_z <- suppressWarnings(min(df$z, na.rm = TRUE))
        if (is.finite(min_z) && min_z < 0) df$z <- df$z - min_z
      }
    }

    df
  })

  # Trace-ID-Length slider UI - range derived from the spread-filtered data
  # (Stage 2 output), so it sits at the end of the pipeline. Stage 3 (filtered_df)
  # reads input$tid_length_range from here.
  output$tid_length_range_ui <- renderUI({
    df <- spread_filtered_df()
    req("tid" %in% names(df))

    tl <- table(df$tid)
    shiny::validate(shiny::need(length(tl) > 0, "No trace lengths available (after metadata + spread filters)."))

    min_len <- min(tl)
    max_len <- max(tl)

    # Clamp the user's current value against the FULL uploaded data's max trace length,
    # not the spread-filtered max.  Without this, applying a spread filter (which removes
    # long traces) would silently lower the upper bound here.  On reset the clamped value
    # would be re-used, permanently excluding any trace whose length exceeded the old
    # spread-filtered maximum - even after the spread filter is cleared.
    df_full     <- uploaded_df()
    tl_full     <- if (!is.null(df_full) && "tid" %in% names(df_full)) table(df_full$tid) else tl
    max_len_full <- if (length(tl_full) > 0) max(tl_full) else max_len

    # On a fresh upload tid_reset_flag is TRUE → always reset to the new data's
    # full range.  On any other re-render (spread filter change, settings load) →
    # preserve the user's current value so they don't lose their TID filter.
    # The flag is read with isolate() so it does NOT create a reactive dependency
    # here; it is cleared immediately after use so the next render (e.g. triggered
    # by a spread filter tweak) falls through to the preserve-curr branch.
    if (isTRUE(isolate(tid_reset_flag()))) {
      isolate(tid_reset_flag(FALSE))
      val <- c(min_len, max_len)
    } else {
      curr <- isolate(input$tid_length_range)
      if (!is.null(curr) && length(curr) == 2 && all(is.finite(as.numeric(curr)))) {
        v1 <- max(min_len, min(as.numeric(curr[1]), max_len_full))
        v2 <- max(min_len, min(as.numeric(curr[2]), max_len_full))
        val <- c(min(v1, v2), max(v1, v2))
      } else {
        val <- c(min_len, max_len)
      }
    }

    shinyWidgets::numericRangeInput(
      "tid_length_range", "Trace ID Length Range",
      value = val, min = min_len, max = max_len
    )
  })
  # ===========================  END SMLM FILTER PIPELINE  =====================

  # --------- 7) SETTINGS LOAD / RESET ----------
  output$download_ui_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
      paste0(base, "_filter_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content  = function(file) {
      settings <- get_ui_settings()
      writeLines(jsonlite::toJSON(settings, pretty=TRUE, auto_unbox=TRUE, na="null"), file, useBytes=TRUE)
    }
  )
  output$download_cluster2_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
      paste0(base, "_cluster2_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content = function(file) {
      s <- list(
        # Master toggle - without it, loading a Median Cluster settings file left
        # clustering switched OFF and none of the restored parameters took effect.
        enable_cluster2 = isTRUE(input$enable_cluster2),
        show_cluster2_knn = isTRUE(input$show_cluster2_knn),

        cluster2_method = input$cluster2_method %||% "dbscan",
        cluster2_xi = input$cluster2_xi,
        cluster2_prefilter = isTRUE(input$cluster2_prefilter),
        cluster2_prefilter_r = input$cluster2_prefilter_r,
        cluster2_prefilter_min = input$cluster2_prefilter_min,
        cluster2_prune = isTRUE(input$cluster2_prune),
        cluster2_prune_tol = input$cluster2_prune_tol,
        cluster2_core_only = isTRUE(input$cluster2_core_only),
        cluster2_min_members = input$cluster2_min_members,
        eps2 = input$eps2,
        minpts2 = input$minpts2,

        show_cluster2_shape = isTRUE(input$show_cluster2_shape),
        show_cluster2_isosurface = isTRUE(input$show_cluster2_isosurface),
        cluster2_iso_selected_only = isTRUE(input$cluster2_iso_selected_only),
        cluster2_iso_id = input$cluster2_iso_id,
        cluster2_iso_color = input$cluster2_iso_color,
        cluster2_iso_unicolor = isTRUE(input$cluster2_iso_unicolor),
        cluster2_iso_by_channel = isTRUE(input$cluster2_iso_by_channel),
        cluster2_iso_grid_res = input$cluster2_iso_grid_res,
        cluster2_iso_pad = input$cluster2_iso_pad,
        cluster2_iso_manual_offset = isTRUE(input$cluster2_iso_manual_offset),
        cluster2_iso_offset = input$cluster2_iso_offset,
        cluster2_iso_opacity = input$cluster2_iso_opacity,
        cluster2_iso_sigma = input$cluster2_iso_sigma,
        cluster2_iso_manual_sigma = isTRUE(input$cluster2_iso_manual_sigma),
        cluster2_iso_smooth_rel = input$cluster2_iso_smooth_rel,
        
        show_cluster2_shape_2D = isTRUE(input$show_cluster2_shape_2D),
        
        show_cluster2_metrics = isTRUE(input$show_cluster2_metrics),
        
        show_cluster2_pca_axes = isTRUE(input$show_cluster2_pca_axes),
        
        show_cluster2_labels = isTRUE(input$show_cluster2_labels),
        show_cluster2_plane_angle = isTRUE(input$show_cluster2_plane_angle),
        show_cluster2_plane_distance = isTRUE(input$show_cluster2_plane_distance),
        show_bending_radius = isTRUE(input$show_bending_radius),
        
        fit_method = input$fit_method,
        fit_weighted = isTRUE(input$fit_weighted),
        mle_center_mode = input$mle_center_mode,
        mle_center_lambda = input$mle_center_lambda,
        
        fit_type = input$fit_type,
        helix_turns = input$helix_turns,
        show_fit_tube = isTRUE(input$show_fit_tube),
        fit_tube_radius = input$fit_tube_radius,
        show_fit_distances = isTRUE(input$show_fit_distances),

        fit_opacity = input$fit_opacity,

        # Fit/metric options that were added after this save list was last updated
        compute_rmsd_metrics     = isTRUE(input$compute_rmsd_metrics),
        membrane_protein_options = isTRUE(input$membrane_protein_options),
        arc_span_deg             = input$arc_span_deg,
        max_dist_lines_nm        = input$max_dist_lines_nm,

        show_cluster2_nn_lines = isTRUE(input$show_cluster2_nn_lines),
        nn_k_cluster2 = input$nn_k_cluster2,
        # NN label styling
        show_c2_nn_labels  = isTRUE(input$show_c2_nn_labels),
        c2_label_font_size = input$c2_label_font_size,
        c2_label_offset    = input$c2_label_offset,
        distance_range_cluster2 = unname(input$distance_range_cluster2) # avoid jsonlite named-vector warning
      )
      
      jsonlite::write_json(s, path = file, auto_unbox = TRUE, pretty = TRUE)
    }
  )

  # ── kNN-distance diagnostic for tuning cluster2 eps ──
  # Plots, for the same centroid set cluster2 operates on (trace medians, or
  # de-blinked molecule medians when de-blinking is on), the sorted distance
  # to the (MinPts-1)-th neighbour. The "elbow" separates dense structures from
  # scattered background = a sensible eps. Red dashed line = current eps2.
  output$cluster2_knn_plot <- renderPlot({
    req(isTRUE(input$enable_cluster2))
    df <- plottable_data()
    req(!is.null(df), nrow(df) > 0)

    use_deblinked <- isTRUE(input$enable_clustering) && "cluster" %in% names(df)
    meds <- if (use_deblinked) {
      df %>%
        dplyr::filter(!is.na(.data$cluster) & .data$cluster != 0) %>%
        dplyr::group_by(.data$cluster) %>%
        dplyr::summarise(mx = stats::median(x, na.rm = TRUE),
                         my = stats::median(y, na.rm = TRUE),
                         mz = stats::median(z, na.rm = TRUE), .groups = "drop")
    } else {
      df %>%
        dplyr::group_by(.data$tid) %>%
        dplyr::summarise(mx = stats::median(x, na.rm = TRUE),
                         my = stats::median(y, na.rm = TRUE),
                         mz = stats::median(z, na.rm = TRUE), .groups = "drop")
    }

    # Must match what median clustering actually measures (always 3D), or the elbow
    # you read eps2 off is on the wrong scale - measured 1.46-1.89x too small on real
    # 3D MINFLUX data when this followed the 2D display toggle.
    cols <- c("mx","my","mz")
    pts <- as.matrix(meds[, cols, drop = FALSE])
    pts <- pts[stats::complete.cases(pts), , drop = FALSE]
    validate(need(nrow(pts) >= 3, "Not enough centroids for a kNN-distance plot."))

    k <- max(1L, as.integer(input$minpts2 %||% 5L) - 1L)
    k <- min(k, nrow(pts) - 1L)
    dbscan::kNNdistplot(pts, k = k)
    graphics::title(main = sprintf("kNN-distance (k = %d, n = %d centroids)", k, nrow(pts)))
    if (isTRUE(input$cluster2_method == "dbscan")) {
      e <- suppressWarnings(as.numeric(input$eps2))
      if (length(e) == 1 && is.finite(e)) graphics::abline(h = e, col = "red", lty = 2, lwd = 2)
    }
  })


  # ── Centroid merge: threshold suggestion + diagnostic ─────────────────────
  # Same shape as the raw-point pair, but measured on CENTROID distances, which
  # is what that criterion actually compares.
  centroid_suggestion <- reactiveVal(NULL)
  # Last measured trace spread, in whatever metric is active. Cached so that
  # editing k can convert to nanometres without re-measuring on every keystroke.
  centroid_sigma_rv <- reactiveVal(NA_real_)

  compute_centroid_suggestion <- function() {
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(NULL)
    w <- if (isTRUE(input$deblink_anisotropic))
      deblink_axis_weights(df, min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT)$w
      else NULL
    r <- suggest_centroid_threshold(
      df,
      min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT,
      use_z = TRUE,          # decided from the data, never from the plot mode
      w     = w,
      k     = { kk <- suppressWarnings(as.numeric(input$centroid_k))
                if (!length(kk) || !is.finite(kk) || kk <= 0) 2.8 else kk })
    if (!is.null(r) && is.finite(r$sigma)) centroid_sigma_rv(r$sigma)
    r
  }

  # ── k drives the distance ────────────────────────────────────────────────
  # The merge itself only ever reads the DISTANCE. k is a convenience for
  # expressing that distance in units of the measured trace spread, so it has to
  # actually move the distance - otherwise editing k changes nothing and the
  # control is decorative. Debounced so typing does not re-measure per keystroke.
  centroid_k_d <- debounce(reactive(input$centroid_k), 600)
  observeEvent(centroid_k_d(), {
    req(isTRUE(input$enable_clustering),
        identical(input$deblink_method %||% "centroid", "centroid"))
    kk <- suppressWarnings(as.numeric(centroid_k_d()))
    if (!length(kk) || !is.finite(kk) || kk <= 0) return()
    sg <- centroid_sigma_rv()
    if (!is.finite(sg)) {                     # never measured yet - measure now
      r <- tryCatch(compute_centroid_suggestion(), error = function(e) NULL)
      if (is.null(r)) return()
      sg <- r$sigma
    }
    newd <- round(kk * sg, 1)
    cur  <- suppressWarnings(as.numeric(input$eps_centroid))
    if (length(cur) == 1 && is.finite(cur) && isTRUE(all.equal(cur, newd))) return()
    updateNumericInput(session, "eps_centroid", value = newd)
  }, ignoreInit = TRUE)

  observeEvent(input$suggest_eps_centroid, {
    r <- compute_centroid_suggestion()
    if (is.null(r)) {
      centroid_suggestion(NULL)
      showNotification("Not enough traces with enough localizations to measure the spread.",
                       type = "warning", duration = 6); return()
    }
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    w  <- if (isTRUE(input$deblink_anisotropic))
      deblink_axis_weights(df, min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT)$w
      else NULL
    # Raw per-axis spread, UNWHITENED. The whitened numbers are identical on
    # every axis by construction, so quoting them three times says nothing.
    r$sigma_ax_raw <- tryCatch(
      suggest_centroid_threshold(df, min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT,
                                 use_z = TRUE, w = NULL, k = r$k)$sigma_ax,
      error = function(e) NULL)
    r$axis_w <- w

    # DELIBERATELY NO DRIFT TERM HERE. Drift is measured, but in the precision
    # QC panel, not in this suggestion - see sample_drift_estimate(). Feeding it
    # back into the merge radius would be circular: the only molecules whose
    # drift can be observed are the ones already merged, i.e. the ones that
    # drifted LESS than the very threshold being set. The threshold is the trace
    # spread times k, and nothing else.
    updateNumericInput(session, "centroid_k", value = r$k)
    centroid_suggestion(r)
    updateNumericInput(session, "eps_centroid", value = r$suggested)
    showNotification(sprintf("Centroid threshold %.1f nm (k = %.2f x σ %.2f nm).",
                             r$suggested, r$k, r$sigma),
                     type = "message", duration = 7)
  })

  output$centroid_suggest_summary <- renderUI({
    r <- centroid_suggestion(); if (is.null(r)) return(NULL)
    axr <- r$sigma_ax_raw %||% r$sigma_ax
    # Derived from the LIVE input, not from the value the suggestor proposed:
    # if the user has since typed a distance by hand, this must say so.
    cur_thr <- suppressWarnings(as.numeric(input$eps_centroid))
    if (!length(cur_thr) || !is.finite(cur_thr)) cur_thr <- r$suggested
    shape <- if (isTRUE(input$deblink_anisotropic) && !is.null(r$axis_w))
      # Whitened per-axis spreads are equal by construction, so the useful
      # number is the real-space merge volume, which is an ellipsoid.
      sprintf("ellipsoid %s nm", paste(sprintf("%.1f", cur_thr / r$axis_w[seq_len(r$dims)]),
                                       collapse = " / "))
      else sprintf("sphere r = %.1f nm", cur_thr)
    tagList(
      helpText(sprintf("Trace spread %s nm  |  threshold %.1f nm = %.2f x spread (%s)  |  %s traces, %.0f%% merge.",
                       paste(sprintf("%.1f", axr), collapse = " / "),
                       cur_thr, if (is.finite(r$sigma) && r$sigma > 0) cur_thr / r$sigma else NA_real_,
                       shape, format(r$n_traces, big.mark = ","),
                       100 * (r$frac_merged %||% NA))),
      helpText("Sample drift is reported in the Localization Precision Cloud QC panel; it is not folded into this threshold.")
    )
  })

  centroid_gap_data <- reactive({
    req(isTRUE(input$view_centroid_gap), isTRUE(input$enable_clustering))
    compute_centroid_suggestion()
  })

  # Drawn by BOTH the single-colour and two-colour panels. Factored so the two
  # can never drift apart visually or in what the lines mean.
  .draw_centroid_gap <- function(r, cur_thr) {
    wi <- r$within[is.finite(r$within) & r$within > 0]
    bw <- r$between[is.finite(r$between) & r$between > 0]
    validate(need(length(wi) >= 2, "Not enough within-trace distances to plot."))

    COL_W <- "#1f78b4"; COL_B <- "#e08214"; COL_S <- "#1a9850"; COL_E <- "#d7191c"
    # LOG x-axis. The re-blink hump sits at a few nm and the distinct-molecule
    # hump at hundreds; on a linear axis clipped to a quantile of the pooled
    # distances the second hump fell off the plot entirely, which defeats the
    # point - the threshold is supposed to be placed in the GAP between them.
    xhi <- max(stats::quantile(bw, 0.99, names = FALSE), stats::quantile(wi, 0.99, names = FALSE),
               3 * r$suggested, na.rm = TRUE)
    xlo <- max(0.5, min(stats::quantile(c(wi, bw), 0.005, names = FALSE), r$sigma / 4))
    scal <- function(v) {
      v <- v[v > 0]
      d <- stats::density(log10(v), n = 512L, from = log10(xlo), to = log10(xhi))
      d$x <- 10^d$x; d$y <- d$y / max(d$y); d
    }
    # Legend goes in the BOTTOM margin. Note xpd is NOT set in par(): doing so
    # let abline() draw its vertical markers out through the margins too, which
    # is what made the lines run past the plot box.
    op <- graphics::par(mar = c(8.2, 4.5, 3.0, 1.2))
    on.exit(graphics::par(op), add = TRUE)
    dw <- scal(wi)
    graphics::plot(dw$x, dw$y, type = "l", col = COL_W, lwd = 2.5, ylim = c(0, 1.05),
                   xlim = c(xlo, xhi), log = "x",
                   xlab = "distance (nm, log scale)", ylab = "density (scaled to peak)",
                   main = "Within-trace scatter vs between-trace centroid spacing")
    if (length(bw) >= 2) { db <- scal(bw); graphics::lines(db$x, db$y, col = COL_B, lwd = 2.5) }
    graphics::abline(v = r$sigma,     col = COL_W, lwd = 2, lty = 3)
    graphics::abline(v = r$suggested, col = COL_S, lwd = 2, lty = 4)
    if (length(cur_thr) == 1 && is.finite(cur_thr))
      graphics::abline(v = cur_thr, col = COL_E, lwd = 2.5, lty = 2)
    # Placed in device coords so the log x-axis cannot distort it; xpd = NA is
    # scoped to this call alone.
    graphics::legend(x = graphics::grconvertX(0.02, "ndc", "user"),
                     y = graphics::grconvertY(0.14, "ndc", "user"),
      bty = "n", cex = 0.85, ncol = 2, xpd = NA,
      legend = c("within trace (loc to its centroid)", "between traces (centroid NN)",
                 "sigma", "suggested", "current threshold"),
      col = c(COL_W, COL_B, COL_W, COL_S, COL_E),
      lty = c(1, 1, 3, 4, 2), lwd = c(2.5, 2.5, 2, 2, 2.5))
    invisible(NULL)
  }

  output$centroid_gap_plot <- renderPlot({
    r <- centroid_gap_data()
    validate(need(!is.null(r), "Not enough traces to measure the spread."))
    validate(need(sum(is.finite(r$within) & r$within > 0) >= 2,
                  "Not enough within-trace distances to plot."))
    .draw_centroid_gap(r, suppressWarnings(as.numeric(input$eps_centroid)))
  })

  # ── PER-CHANNEL TRACE SPREAD ──────────────────────────────────────────────
  # The 2-colour panel used to run this same estimator POOLED over both channels
  # and set a separate eps_tc1 from it. Its own comment admitted the problem:
  # "two dyes can genuinely differ in precision, and if they do this under-serves
  # the tighter channel". They do differ - photon yield differs, so localization
  # precision differs, so trace spread differs - and the merge distance is k x
  # that spread. So it is measured per channel now, and the pooled number is a
  # labelled fallback rather than the default.
  #
  # k stays SHARED. It says how many spreads apart two blinks can be and still be
  # one molecule, which is a statement about blinking, not about the dye.
  channel_spread_table <- reactive({
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- tryCatch(attach_channel_colors(resolve_channels(df, input), input),
                   error = function(e) NULL)
    if (is.null(df) || !nrow(df) || !"channel" %in% names(df)) return(NULL)
    thr  <- tryCatch(deblink_threshold_for_channel(df, input), error = function(e) NULL)
    if (is.null(thr)) return(NULL)
    ct   <- attr(df, "channel_table")
    meas <- attr(thr, "measured"); sig <- attr(thr, "sigma")
    ct$sigma    <- as.numeric(sig[as.character(ct$channel)])
    ct$distance <- as.numeric(unlist(thr[as.character(ct$channel)]))
    ct$measured <- as.logical(meas[as.character(ct$channel)])
    ct
  })

  # ── WHAT THE CHANNEL RULE COSTS ───────────────────────────────────────────
  # "Never merge across a channel" is the right rule, and it is not free. Both
  # costs below bias the de-blinked molecule count DOWNWARD, and neither was
  # visible anywhere before this readout.
  #
  # 1. AMBIGUOUS DCR. A trace is assigned to a channel by its MEDIAN DCR, so two
  #    blinks of one molecule that land on opposite sides of the split can never
  #    be merged - the molecule is counted twice, once per channel. The size of
  #    the effect is already computed as the GMM posterior behind the split
  #    confidence plot; it is just never stated where the merge is configured.
  #
  #    Note this is the SAME root cause as the crosstalk warning on the split
  #    itself, but a different consequence. Crosstalk pushes cross-channel
  #    statistics UP at short r (false co-localization); here it pushes molecule
  #    counts DOWN. One mis-assignment, two directions.
  #
  # 2. STRADDLING TRACES (multiplex only). A multiplex channel is assigned per
  #    LOCALIZATION, so a trace blinking across a TIM breakpoint sits in two
  #    channels and is split into two molecules. Merging it would invent a
  #    molecule that was never simultaneously present, so the split is correct -
  #    but it is a real count, not a rounding error, and it belongs on screen.
  output$channel_integrity_summary <- renderUI({
    items <- list()

    # 1. ambiguous DCR band
    if (isTRUE(input$enable_two_color)) {
      tb <- tryCatch(dcr_conf_tbl(), error = function(e) NULL)
      if (is.data.frame(tb) && nrow(tb)) {
        n_edge <- sum(tb$on_edge, na.rm = TRUE)
        if (n_edge > 0)
          items <- c(items, list(tags$div(
            sprintf("%s of %s traces sit in the ambiguous DCR band (%.1f%%). A molecule whose blinks straddle the split is counted once in each channel, so the molecule count is biased low by about that fraction.",
                    format(n_edge, big.mark = ","), format(nrow(tb), big.mark = ","),
                    100 * n_edge / nrow(tb)))))
      }
    }

    # 2. traces straddling a multiplex boundary
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    if (is.data.frame(df) && nrow(df) && "mplex_channel" %in% names(df) &&
        any(!is.na(df$mplex_channel))) {
      per_tid <- tapply(df$mplex_channel, df$tid, function(v) length(unique(v[!is.na(v)])))
      n_str <- sum(per_tid > 1L, na.rm = TRUE)
      if (n_str > 0)
        items <- c(items, list(tags$div(style = "margin-top:4px;",
          sprintf("%s of %s traces span more than one multiplex channel and are split into one molecule per channel. Merging them would claim a molecule was present in two exchange rounds at once.",
                  format(n_str, big.mark = ","), format(length(per_tid), big.mark = ",")))))
    }

    if (!length(items)) return(NULL)
    tags$div(style = "margin-top:8px; padding:6px 9px; border-left:3px solid #e0c98a; background:#fdfaf2; font-size:0.82em; color:#5a4a20;",
             tags$b("What the channel rule costs here"), tags$br(), items)
  })

  output$channel_spread_summary <- renderUI({
    ct <- tryCatch(channel_spread_table(), error = function(e) NULL)
    if (is.null(ct) || nrow(ct) < 2L) return(NULL)

    rows <- lapply(seq_len(nrow(ct)), function(i) {
      col <- ct$color[i]; if (is.na(col)) col <- "#888888"
      tags$div(style = "display:flex; align-items:center; gap:8px; font-size:0.86em; margin:2px 0;",
        tags$span(style = sprintf("width:10px;height:10px;border-radius:2px;background:%s;display:inline-block;", col)),
        tags$b(as.character(ct$label[i])),
        tags$span(if (isTRUE(ct$measured[i]))
                    sprintf("spread %.1f nm  ->  merge at %.1f nm", ct$sigma[i], ct$distance[i])
                  else sprintf("spread not measurable  ->  shared %.1f nm", ct$distance[i])),
        tags$span(style = "color:#7b8794;", sprintf("(%s locs)", format(ct$n[i], big.mark = ","))))
    })

    # The spread RATIO is the number that decides whether one shared distance is
    # defensible, so state it rather than making the user divide two numbers.
    sg <- ct$sigma[is.finite(ct$sigma)]
    note <- if (length(sg) >= 2) {
      rr <- max(sg) / min(sg)
      if (isTRUE(input$deblink_shared_spread))
        helpText(style = "color:#7a5c00;", sprintf(
          "Using ONE shared distance while the channels differ in spread by %.2fx. The broader channel is under-merged and the tighter one over-merged.", rr))
      else if (rr >= 1.15)
        helpText(sprintf("Channels differ in spread by %.2fx, so they are merged at different distances.", rr))
      else
        helpText(sprintf("Channels agree in spread to within %.2fx - a shared distance would give nearly the same result here.", rr))
    } else NULL

    tagList(tags$div(style = "margin:4px 0 2px; font-weight:600; font-size:0.86em;",
                     "Merge distance per channel"), rows, note)
  })

  # ── EXPLORATORY: molecule count vs radius ─────────────────────────────────
  # Gated on its own checkbox: the sweep runs one DBSCAN per radius, so it must
  # not fire on every filter change.
  deblink_knee_data <- reactive({
    req(isTRUE(input$view_deblink_knee), isTRUE(input$enable_clustering))
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    req(!is.null(df), nrow(df) > 0)
    # Same metric the clustering uses, or the radii on the x-axis mean nothing.
    w <- if (isTRUE(input$deblink_anisotropic))
      deblink_axis_weights(df, min_n = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT)$w
      else NULL
    # Sweep whatever criterion is ACTIVE - a curve over raw-point eps says
    # nothing about a centroid threshold, and vice versa.
    # Sweep range follows the DATA, not a hard-coded 40 nm: on the centroid
    # criterion the interesting structure runs well past that.
    # Centroid criterion only - see deblink_count_curve().
    req(identical(input$deblink_method %||% "centroid", "centroid"))
    cs  <- tryCatch(compute_centroid_suggestion(), error = function(e) NULL)
    top <- if (!is.null(cs))
      max(40, ceiling(stats::quantile(cs$between, 0.75, names = FALSE)), 4 * cs$suggested)
      else 40
    deblink_count_curve(df, minpts = deblink_minpts(input),
                        eps_seq = seq(1, min(top, 400), length.out = 60),
                        w = w, method = "centroid")
  })

  output$deblink_knee_plot <- renderPlot({
    r <- deblink_knee_data()
    validate(need(!is.null(r), "Not enough data to sweep the radius."))
    cv <- r$curve
    graphics::plot(cv$eps, cv$n_molecules, type = "l", lwd = 2.5, col = "#1f78b4",
                   xlab = if (identical(input$deblink_method %||% "centroid", "vote"))
                            "vote fraction f" else "centroid distance threshold (nm)",
                   ylab = "de-blinked molecules",
                   main = "Molecule count vs merge threshold")
    graphics::points(cv$eps, cv$n_molecules, pch = 16, cex = 0.5, col = "#1f78b4")
    graphics::abline(h = r$n_traces, col = "#969696", lty = 3, lwd = 2)   # no merging at all
    if (is.finite(r$knee)) graphics::abline(v = r$knee, col = "#e08214", lwd = 2, lty = 4)
    e_cur <- suppressWarnings(as.numeric(
      if (identical(input$deblink_method %||% "centroid", "vote")) input$vote_frac
      else input$eps_centroid))
    if (length(e_cur) == 1 && is.finite(e_cur))
      graphics::abline(v = e_cur, col = "#d7191c", lwd = 2.5, lty = 2)
    graphics::legend("topright", bty = "n", cex = 0.9,
      legend = c("molecules", "traces (no merging)", "knee", "current threshold"),
      col = c("#1f78b4", "#969696", "#e08214", "#d7191c"),
      lty = c(1, 3, 4, 2), lwd = c(2.5, 2, 2, 2.5))
  })

  output$deblink_knee_summary <- renderUI({
    r <- tryCatch(deblink_knee_data(), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    cv <- r$curve
    at <- function(e) { i <- which.min(abs(cv$eps - e)); cv$n_molecules[i] }
    e_cur <- suppressWarnings(as.numeric(
      if (identical(input$deblink_method %||% "centroid", "vote")) input$vote_frac
      else input$eps_centroid))
    tagList(
      helpText(sprintf("%s traces → %s molecules at the knee (%.0f nm)%s.",
                       format(r$n_traces, big.mark = ","),
                       if (is.finite(r$knee)) format(at(r$knee), big.mark = ",") else "-",
                       r$knee,
                       if (length(e_cur) == 1 && is.finite(e_cur))
                         sprintf("; %s at your current %.0f nm",
                                 format(at(e_cur), big.mark = ","), e_cur) else "")),
      if (r$n_locs < 60000L) NULL else
        helpText("Curve computed on a 60,000-localization subsample."),
      # Measured: the knee detector is calibrated for the raw-point curve, whose
      # steep initial drop comes from traces FRAGMENTING at small eps. The
      # centroid curve has no such regime - it simply starts at "no merging" and
      # declines - so the curvature maximum lands mid-slope and reads far too
      # low (synthetic: knee 5 nm -> 422 molecules, against 360 truth and 360 at
      # the σ-based suggestion). Do not tune the centroid threshold from it.
      if (identical(input$deblink_method %||% "centroid", "centroid"))
        helpText(tags$b("Unreliable for this criterion: "),
                 "the knee is calibrated for the raw-point curve shape and reads too low ",
                 "on centroid distances. Use 'Suggest threshold from trace spread' instead; ",
                 "the curve itself is still worth reading.")
      else NULL
    )
  })


  output$download_two_color_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
      paste0(base, "_channel_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content = function(file) {
      s <- list(
        # Master toggle - as with cluster2, without it the restored 2-color
        # parameters had no effect because 2-color stayed switched off.
        enable_two_color = isTRUE(input$enable_two_color),
        color_by_fluo = isTRUE(input$color_by_fluo),
        use_pooled_dcr = isTRUE(input$use_pooled_dcr),
        dcr_group1 = unname(input$dcr_group1),
        dcr_group2 = unname(input$dcr_group2),
        color_group1 = input$color_group1,
        color_group2 = input$color_group2,

        
        # Clustering parameters are NOT written here. One pipeline means one set
        # of de-blink / Median Clustering controls, and those are saved by the
        # main settings file. Writing a second copy here is what let the two
        # drift apart in the first place.
        show_tc1_alpha = isTRUE(input$show_tc1_alpha),
        tc1_alpha = input$tc1_alpha,
        tc1_alpha_opacity = input$tc1_alpha_opacity,
        tc1_alpha_min_pts = input$tc1_alpha_min_pts
        # Median Cluster wrap / isosurface settings are saved with Median
        # Clustering, not here.
      )
      
      jsonlite::write_json(s, path = file, auto_unbox = TRUE, pretty = TRUE)
    }
  )
  output$download_filtered_csv <- downloadHandler(
    filename = function() {
      paste0(
        tools::file_path_sans_ext(input$active_file %||% "smlm"),
        "_filtered_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".csv"
      )
    },
    
    content = function(file) {
      df <- plottable_data()

      shiny::validate(
        shiny::need(!is.null(df) && nrow(df) > 0, "No data to export.")
      )

      # Core MINFLUX columns - always included when present
      core_cols <- c("x", "y", "z", "tid", "efo", "cfr", "eco", "fbg", "dcr")
      # Optional columns included only when they exist and carry real data
      # dcr_pooled is exported so the pooled value SURVIVES a round-trip: it is
      # derived from the per-iteration dcr1, which a CSV does not carry, so
      # without this it could never be recomputed after re-import.
      optional_cols <- c("fluo", "tim", "itr", "dwell", "dcr_pooled",
                         "pcid", "scid", "n_in_cluster2",
                         ".dataset_tag", "analysis_group")

      keep <- c(
        intersect(core_cols, names(df)),
        intersect(optional_cols, names(df))
      )
      # Drop optional columns that are entirely NA (e.g. clustering not run)
      keep <- keep[sapply(keep, function(col) !all(is.na(df[[col]])))]

      df <- df[, keep, drop = FALSE]

      # Friendly export names
      names(df)[names(df) == ".dataset_tag"]  <- "dataset_tag"
      names(df)[names(df) == "n_in_cluster2"] <- "n_in_scid"

      # Axial provenance: z below is written at this factor. Re-importing divides
      # it out before applying the then-current calibration, so a round-trip is
      # exact instead of relying on the 0.7 assumption in
      # detect_csv_z_calibration(). Written even for 2D data - it costs one
      # column and removes the guess for every future reader.
      # The COMMITTED factor - that is what the exported z is actually at. Using
      # the pending input value would mislabel the file whenever the control has
      # been typed into but not applied.
      df$z_scale_applied <- z_scale_applied()

      write.csv(df, file, row.names = FALSE)
    }
  )
  observeEvent(input$download_filtered_csv, {
    showNotification("Exported filtered dataset.", type = "message", duration = 3)
  })
  observeEvent(input$load_settings, {
    req(input$load_settings$datapath)
    raw <- jsonlite::fromJSON(input$load_settings$datapath)
    
    base_fil <- filters_defaults_from_data(uploaded_df())
    fil <- modifyList(base_fil, raw, keep.null = TRUE)
    ov  <- modifyList(overlay_defaults(), raw$overlay %||% list(), keep.null = TRUE)
    hm  <- modifyList(heatmap_defaults(),  raw$heatmap %||% list(), keep.null = TRUE)
    
    # 1) Push everything EXCEPT tid_length_range first
    fil_no_tid <- fil
    fil_no_tid$tid_length_range <- NULL
    
    session$onFlushed(function() {
      apply_filters(fil_no_tid)
      # Restore view_SD from the saved file (apply_filters no longer touches it).
      if (!is.null(fil$view_SD)) updateCheckboxInput(session, "view_SD", value = isTRUE(fil$view_SD))
      # Restore the min-localizations-for-spread setting.
      if (!is.null(fil$sd_min_locs)) {
        v <- suppressWarnings(as.integer(fil$sd_min_locs))
        if (length(v) == 1 && is.finite(v)) updateNumericInput(session, "sd_min_locs", value = v)
      }

      # Restore trimming BEFORE the ranges below: it changes which localizations
      # exist, so the nm bounds have to land on the same point set they were set
      # on. Files written before this setting existed carry no field and were
      # saved with no trimming.
      updateCheckboxInput(session, "trim_outliers", value = isTRUE(fil$trim_outliers))
      if (!is.null(fil$trim_k)) {
        v <- suppressWarnings(as.numeric(fil$trim_k))
        if (length(v) == 1 && is.finite(v)) updateNumericInput(session, "trim_k", value = v)
      }

      # ── Restore spatial crop / trimmer (changes the point set) ─────────────
      cr <- raw$crop %||% list()
      # 2D crop box: a server-side reactiveVal, so no DOM round-trip needed.
      cx <- suppressWarnings(as.numeric(unlist(cr$crop2d_x)))
      cy <- suppressWarnings(as.numeric(unlist(cr$crop2d_y)))
      if (length(cx) == 2 && length(cy) == 2 && all(is.finite(c(cx, cy))))
        smlm_2d_crop_rv(list(x = cx, y = cy)) else smlm_2d_crop_rv(NULL)
      # 3D box/sphere crop: static inputs restore directly.
      updateCheckboxInput(session, "enable_trimmer", value = isTRUE(cr$enable_trimmer))
      if (!is.null(cr$trim_mode)) updateRadioButtons(session, "trim_mode", selected = cr$trim_mode)
      for (nm in c("trim_cx", "trim_cy", "trim_cz", "trim_r")) {
        v <- suppressWarnings(as.numeric(cr[[nm]]))
        if (length(v) == 1 && is.finite(v)) updateNumericInput(session, nm, value = v)
      }
      # Re-activate the committed crop region so the reloaded file is cropped as saved.
      ta <- cr$trim_applied
      if (!is.null(ta) && !is.null(ta$mode)) {
        if (identical(ta$mode, "sphere")) {
          cc <- suppressWarnings(as.numeric(unlist(ta$c))); rr <- suppressWarnings(as.numeric(ta$r))
          trim_applied_rv(if (length(cc) == 3 && all(is.finite(cc)) && is.finite(rr))
            list(mode = "sphere", c = cc, r = rr) else NULL)
        } else {
          num2 <- function(v) { v <- suppressWarnings(as.numeric(unlist(v))); if (length(v) == 2 && all(is.finite(v))) v else NULL }
          trim_applied_rv(list(mode = "box", x = num2(ta$x), y = num2(ta$y), z = num2(ta$z)))
        }
      } else trim_applied_rv(NULL)
      # Box ranges live in a renderUI (trim_box_ui); set them after it rebuilds.
      if (isTRUE(cr$enable_trimmer) && identical(cr$trim_mode %||% "box", "box")) {
        bx <- suppressWarnings(as.numeric(unlist(cr$trim_x_range)))
        by <- suppressWarnings(as.numeric(unlist(cr$trim_y_range)))
        bz <- suppressWarnings(as.numeric(unlist(cr$trim_z_range)))
        later::later(function() {
          if (length(bx) == 2 && all(is.finite(bx))) shinyWidgets::updateNumericRangeInput(session, "trim_x_range", value = bx)
          if (length(by) == 2 && all(is.finite(by))) shinyWidgets::updateNumericRangeInput(session, "trim_y_range", value = by)
          if (length(bz) == 2 && all(is.finite(bz))) shinyWidgets::updateNumericRangeInput(session, "trim_z_range", value = bz)
          plot_trigger(Sys.time())
        }, delay = 0.4)
      }

      # ── Restore fiducial-scan exclusion + its thresholds ───────────────────
      fd <- raw$fiducial %||% list()
      updateCheckboxInput(session, "fid_exclude",   value = isTRUE(fd$fid_exclude))
      updateCheckboxInput(session, "fid_all_files", value = isTRUE(fd$fid_all_files))
      for (nm in c("fid_efo_z", "fid_dur_x", "fid_min_locs", "fid_len_z", "fid_recur_r", "fid_recur_min")) {
        v <- suppressWarnings(as.numeric(fd[[nm]]))
        if (length(v) == 1 && is.finite(v)) updateNumericInput(session, nm, value = v)
      }

      # ── Restore analysis options that change results ───────────────────────
      an <- raw$analysis %||% list()
      for (nm in c("use_pooled_dcr", "spatial_edge_correct", "spatial_use_support",
                   "dens_norm_enable", "dens_thin_enable"))
        if (!is.null(an[[nm]])) updateCheckboxInput(session, nm, value = isTRUE(an[[nm]]))
      # Which points draw the observation window changes every area on the tab.
      if (!is.null(an$dens_support_from))
        updateRadioButtons(session, "dens_support_from", selected = an$dens_support_from)

      # ── Restore Local density & hotspots ──────────────────────────────
      hs <- raw$hotspot %||% list()
      for (nm in c("hot_k", "hot_enh_thresh", "hot_quadrat", "hot_minpts",
                   "hot_target_frac", "hot_nsim", "hot_box_top")) {
        v <- suppressWarnings(as.numeric(hs[[nm]]))
        if (length(v) == 1 && is.finite(v)) updateNumericInput(session, nm, value = v)
      }
      for (nm in c("hot_unit", "hot_palette", "hot_metric"))
        if (!is.null(hs[[nm]])) updateSelectInput(session, nm, selected = hs[[nm]])
      if (!is.null(hs$hot_curve_mode))
        updateRadioButtons(session, "hot_curve_mode", selected = hs$hot_curve_mode)
      for (nm in c("hot_use_3d", "hot_show_foci"))
        if (!is.null(hs[[nm]])) updateCheckboxInput(session, nm, value = isTRUE(hs[[nm]]))

      # Update sd_filt so the actual spread filter reflects the loaded settings.
      # NA in the upper bound means the filter was at full range when saved
      # (Inf serialises to null/NA in JSON); restore it as Inf.
      load_sd_rng <- function(v) {
        if (is.null(v)) return(NULL)
        r <- as.numeric(unlist(v))
        if (length(r) != 2) return(NULL)
        if (is.na(r[2]) || is.infinite(r[2])) r[2] <- Inf
        r
      }
      loaded_x <- load_sd_rng(fil$sd_x_range)
      loaded_y <- load_sd_rng(fil$sd_y_range)
      loaded_z <- load_sd_rng(fil$sd_z_range)
      if (!is.null(loaded_x)) sd_filt$x <- loaded_x
      if (!is.null(loaded_y)) sd_filt$y <- loaded_y
      if (!is.null(loaded_z)) sd_filt$z <- loaded_z

      # Set tid_length_range after the filter UI is in the DOM; coerce to numeric
      # in case jsonlite parsed the JSON array as a list.
      if (!is.null(fil$tid_length_range) && length(fil$tid_length_range) == 2) {
        df <- isolate(metadata_filtered_df())
        tl <- table(df$tid)
        if (length(tl)) {
          min_len <- min(tl); max_len <- max(tl)
          v <- as.numeric(unlist(fil$tid_length_range))
          v[1] <- max(min_len, min(v[1], max_len))
          v[2] <- max(min_len, min(v[2], max_len))
          shinyWidgets::updateNumericRangeInput(session, "tid_length_range", value = v)
        }
      }

      apply_overlay(ov)
      apply_heatmap(hm)
      # Re-run the plotted/exported pipeline so a restored 2D crop or trimmer
      # takes effect (plottable_data is gated on plot_trigger).
      plot_trigger(Sys.time())
    }, once = TRUE)
  })
  cluster2_file_reset_key <- reactiveVal(0)
  
  output$cluster2_settings_file_ui <- renderUI({
    cluster2_file_reset_key()  
    fileInput("cluster2_settings_file", "Load Median Cluster settings (.json)", accept = c(".json"))
  })
  
  observeEvent(input$cluster2_settings_file, {
    req(input$cluster2_settings_file$datapath)
    
    s <- jsonlite::read_json(
      input$cluster2_settings_file$datapath,
      simplifyVector = TRUE
    )
    
    apply_cluster2_settings(session, s)
    
    showNotification(
      paste0(
        "Median Cluster settings loaded: ",
        input$cluster2_settings_file$name
      ),
      type = "message",
      duration = 3
    )
  })
  two_color_file_reset_key <- reactiveVal(0)
  
  output$two_color_settings_file_ui <- renderUI({
    two_color_file_reset_key()
    fileInput("two_color_settings_file", "Load channel settings (.json)", accept = c(".json"))
  })
  observeEvent(input$two_color_settings_file, {
    req(input$two_color_settings_file$datapath)
    s <- jsonlite::read_json(input$two_color_settings_file$datapath, simplifyVector = TRUE)
    apply_two_color_settings(session, s)
    
    showNotification(
      paste0("2-color settings loaded: ", input$two_color_settings_file$name),
      type = "message", duration = 3
    )
  })
  observeEvent(input$minpts2 < 4L, {
    mp <- as.integer(input$minpts2 %||% 4L)
    if (!is.na(mp) && mp < 4L) {
      updateNumericInput(session, "minpts2", value = 4L)
      showNotification("Median Cluster minPts raised to 4 to keep geometry stable (convex hull / fits need ≥4 points).",
                       type = "warning", duration = 4)
    }
  })
  
  observeEvent(input$reset_filters, {
    req(input$file)
    was_view_SD <- isTRUE(isolate(input$view_SD))
    vals <- filters_defaults_from_data(uploaded_df())
    # Reset SD filter synchronously - Inf upper bound means "pass all".
    sd_filt$x <- c(0, Inf); sd_filt$y <- c(0, Inf); sd_filt$z <- c(0, Inf)
    # Trimming removes localizations, so "reset filters" has to put them back -
    # otherwise the spread the user returns to is not the raw spread.
    updateCheckboxInput(session, "trim_outliers", value = FALSE)
    shinyjs::reset("filters_form")
    session$onFlushed(function() {
      sd_skip_arming(TRUE)
      apply_filters(vals)
      if (was_view_SD) updateCheckboxInput(session, "view_SD", value = TRUE)
      # Update meta_full_rng to the full-range values NOW so the all.equal bypass
      # in metadata_filtered_df passes when the plot re-evaluates after the delay.
      meta_full_rng$efo   <- vals$efo_range
      meta_full_rng$cfr   <- vals$cfr_range
      meta_full_rng$tim   <- vals$tim_range
      meta_full_rng$dcr   <- vals$dcr_range
      meta_full_rng$eco   <- vals$eco_range
      meta_full_rng$dwell <- vals$dwell_range
      meta_full_rng$fluo  <- vals$fluo_range
      # Fire plot_trigger after a short delay to let the browser finish processing
      # all the widget updates sent by apply_filters above.  Without this delay,
      # input$tid_length_range / input$efo_range etc. still hold old filtered
      # values when the plot re-evaluates.
      later::later(function() {
        sd_skip_arming(FALSE)
        if (!is.null(isolate(input$file))) plot_trigger(Sys.time())
      }, delay = 0.5)
    }, once = TRUE)
    shinyjs::runjs("const el = document.getElementById('load_settings'); if (el) { el.value=''; $(el).trigger('change'); }")
    showNotification("Filters reset to data ranges.", type = "message", duration = 3)
  })
  # When the SD panel opens its renderUI creates fresh widgets, which fire their
  # initial values as input changes. Suppress those for 800 ms so they are not
  # mistaken for a user-driven narrowing of the filter.
  observeEvent(input$view_SD, {
    if (isTRUE(input$view_SD)) {
      sd_skip_arming(TRUE)
      later::later(function() sd_skip_arming(FALSE), delay = 0.8)
    }
  }, ignoreInit = TRUE)
  # Commit the user's chosen range into sd_filt - the single source of truth.
  # Uses Inf for any axis whose widget is at the full-range maximum so that
  # spread_filtered_df always passes those traces without a fragile value comparison.
  observeEvent(list(input$sd_x_range, input$sd_y_range, input$sd_z_range), {
    req(input$view_SD)
    if (isolate(sd_skip_arming())) return()
    sdr   <- isolate(sd_grey_tbl())
    get_mx <- function(col) { v <- sdr[[col]]; v <- v[is.finite(v)]; if (!length(v)) Inf else max(v) }
    # True when the widget value is at (or indistinguishably close to) the data maximum.
    # A tiny relative tolerance handles JSON float round-trip precision loss, which can
    # cause rx[2] to come back as mx * (1 - ~1e-15), making the strict >= check fail and
    # setting sd_filt to a finite value that then excludes the highest-spread traces.
    at_max <- function(val, mx) is.infinite(val) || (is.finite(mx) && val >= mx - 1e-9 * max(1, abs(mx)))
    # Clearing a box makes numericRangeInput report NA, not NULL, so `%||%` never
    # sees it and the NA propagates through at_max() into `if`, which stops the
    # session with "missing value where TRUE/FALSE is required". An empty box is
    # read as "no bound on that side": blank lower = 0, blank upper = no upper
    # limit, which is the same thing the widget means when it sits at the data
    # maximum. A reversed pair is swapped rather than refused, so dragging the
    # upper below the lower cannot wedge the panel either.
    san <- function(v) {
      v  <- suppressWarnings(as.numeric(v))
      lo <- if (length(v) >= 1L && is.finite(v[1])) max(0, v[1]) else 0
      hi <- if (length(v) >= 2L && !is.na(v[2]))    v[2]         else Inf
      if (is.finite(hi) && hi < lo) { tmp <- lo; lo <- hi; hi <- tmp }
      c(lo, hi)
    }
    mx_x  <- get_mx("sd_x"); mx_y <- get_mx("sd_y")
    rx <- san(input$sd_x_range)
    ry <- san(input$sd_y_range)
    rz <- san(input$sd_z_range)
    sd_filt$x <- if (rx[1] <= 0 && at_max(rx[2], mx_x)) c(0, Inf) else rx
    sd_filt$y <- if (ry[1] <= 0 && at_max(ry[2], mx_y)) c(0, Inf) else ry
    if ("sd_z" %in% names(sdr)) {
      mx_z <- get_mx("sd_z")
      sd_filt$z <- if (rz[1] <= 0 && at_max(rz[2], mx_z)) c(0, Inf) else rz
    }
  }, ignoreInit = TRUE)
  observeEvent(input$reset_cluster2_settings, {
    apply_cluster2_settings(session, cluster2_defaults())
    
    # Force fileInput to be recreated
    cluster2_file_reset_key(cluster2_file_reset_key() + 1)
    
    showNotification("Median Cluster settings reset (and file cleared).", type = "message", duration = 3)
  })
  
  # Uncheck "within cluster2" when cluster2 DBSCAN is disabled
  observeEvent(input$enable_cluster2, {
    if (!isTRUE(input$enable_cluster2) && isTRUE(input$nn_within_cluster2)) {
      updateCheckboxInput(session, "nn_within_cluster2", value = FALSE)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$reset_two_color_settings, {
    apply_two_color_settings(session, two_color_defaults())
    
    # Force file input to be recreated
    two_color_file_reset_key(two_color_file_reset_key() + 1)
    
    showNotification("2-color settings reset to defaults.", type = "message", duration = 3)
  })
  
  
  observeEvent(input$load_overlay_settings, {
    req(input$load_overlay_settings)
    ov <- jsonlite::fromJSON(input$load_overlay_settings$datapath)
    cur_name <- if (!is.null(input$bg_image)) input$bg_image$name else NULL
    apply_overlay_settings(ov, session, current_bg_name = cur_name)
    showNotification("Overlay settings applied.", type = "message", duration = 4)
  })
  observeEvent(input$reset_overlay, {
    ov <- overlay_defaults()
    apply_overlay_settings(ov, session, current_bg_name = NULL)
    shinyjs::reset("overlay_form")
    showNotification("Overlay reset to defaults.", type = "message", duration = 3)
  })
  # --------- 7) Trace spread inspection  ----------
  sd_grey_tbl <- reactive({
    sd_by_tid(raw_df(), min_n = input$sd_min_locs)
  })

  sd_black_tbl <- reactive({
    sd_by_tid(filtered_df(), min_n = input$sd_min_locs)
  })

  output$sd_summary_ui <- renderUI({
    req(input$view_SD)
    df  <- filtered_df()    
    med <- compute_sd_medians(df, min_n = input$sd_min_locs)
    
    # One axis per line, in the pyMINFLUX reporting form: median +/- MAD with the
    # interquartile range after it. The MAD is the width of the trace-to-trace
    # distribution, NOT an error on the median - the IQR is shown next to it so
    # the asymmetry of that distribution stays visible rather than being hidden
    # behind a symmetric +/-.
    fmt_axis <- function(axis, med_v, mad_v, q1, q3) {
      if (!is.finite(med_v %||% NA_real_)) return(NULL)
      txt <- sprintf("Median trace spread %s: %.2f ± %.2f nm", axis, med_v, mad_v %||% NA_real_)
      if (is.finite(q1 %||% NA_real_) && is.finite(q3 %||% NA_real_))
        txt <- sprintf("%s  (IQR %.2f–%.2f)", txt, q1, q3)
      helpText(txt)
    }

    els <- list(
      fmt_axis("X", med$x, med$mad_x, med$q1_x, med$q3_x),
      fmt_axis("Y", med$y, med$mad_y, med$q1_y, med$q3_y)
    )
    if (!is.null(med$z)) {
      els <- append(els, list(fmt_axis("Z", med$z, med$mad_z, med$q1_z, med$q3_z)))
    }
    els <- Filter(Negate(is.null), els)
    if (length(els) && is.finite(med$n_traces %||% NA_real_))
      els <- append(els, list(helpText(sprintf("Across %d traces with ≥ %d localizations.",
                                               as.integer(med$n_traces),
                                               as.integer(input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT)))))
    do.call(tagList, els)
    })
  output$sd_ranges_ui <- renderUI({
    req(input$view_SD)
    sdr <- sd_grey_tbl()

    get_max <- function(v) {
      v <- v[is.finite(v)]
      if (!length(v)) return(1)
      mx <- max(v, na.rm = TRUE)
      if (!is.finite(mx) || mx <= 0) 1 else mx
    }
    # Initialise widget from sd_filt (the committed filter state) so that a
    # previously loaded or user-set range is shown when the panel is reopened.
    # Falls back to [0, mx] when the filter is at full range (Inf upper bound).
    # When "lock filters" is active, preserve the locked upper bound exactly and
    # extend the widget max to fit it - prevents clamping from silently resetting
    # sd_filt via the sd_x_range observer when the new file has a smaller max SD.
    locked <- isTRUE(input$lock_filters)
    init_val <- function(filt_rng, mx) {
      lo <- max(0, filt_rng[1])
      hi <- if (is.infinite(filt_rng[2])) mx
            else if (locked) filt_rng[2]
            else min(mx, filt_rng[2])
      c(lo, hi)
    }
    eff_max <- function(mx, hi) if (locked && is.finite(hi) && hi > mx) hi else mx

    mx_x <- get_max(sdr$sd_x)
    mx_y <- get_max(sdr$sd_y)
    # Use isolate() so that changes to sd_filt do NOT cause sd_ranges_ui to
    # re-render.  Without isolate, setting sd_filt inside the SD range observer
    # would invalidate sd_ranges_ui → new widgets fire → observer fires again →
    # loop.  With isolate, sd_filt is read once when the panel opens and the
    # output is not tracked as a reactive dependency of sd_filt.
    cur <- isolate(list(x = sd_filt$x, y = sd_filt$y, z = sd_filt$z))
    vx  <- init_val(cur$x, mx_x)
    vy  <- init_val(cur$y, mx_y)

    tagList(
      shinyWidgets::numericRangeInput("sd_x_range", "Allowable trace spread in X",
        value = vx, min = 0, max = eff_max(mx_x, vx[2])),
      shinyWidgets::numericRangeInput("sd_y_range", "Allowable trace spread in Y",
        value = vy, min = 0, max = eff_max(mx_y, vy[2])),
      if ("sd_z" %in% names(sdr) && any(is.finite(sdr$sd_z))) {
        mx_z <- get_max(sdr$sd_z)
        vz   <- init_val(cur$z, mx_z)
        shinyWidgets::numericRangeInput("sd_z_range", "Allowable trace spread in Z",
          value = vz, min = 0, max = eff_max(mx_z, vz[2]))
      },
      helpText("An empty box means no bound on that side: blank lower reads as 0, ",
               "blank upper as no upper limit.")
    )
  })
  observeEvent(input$view_SD, {
    req(input$file)
    df  <- filtered_df()
    med <- compute_sd_medians(df, min_n = input$sd_min_locs)
    isolate({
      updateNumericInput(session, "target_sd_x", value = round(med$x %||% NA_real_, 2))
      updateNumericInput(session, "target_sd_y", value = round(med$y %||% NA_real_, 2))
      if (!is.null(med$z)) {
        updateNumericInput(session, "target_sd_z", value = round(med$z %||% NA_real_, 2))
      }
    })
  })
  # If user toggles q95 checkbox, immediately refresh the preset ranges
  observeEvent(input$sd_use_q95, {
    req(input$view_SD)
    sdr <- sd_grey_tbl()
    req(is.data.frame(sdr), nrow(sdr) > 0)
    
    set_preset <- function(col, inputId) {
      if (!col %in% names(sdr)) return()
      v <- as.numeric(sdr[[col]])
      v <- v[is.finite(v)]
      if (!length(v)) return()
      
      mx  <- max(v, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1
      top <- if (isTRUE(input$sd_use_q95)) stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE) else mx
      if (!is.finite(top) || top <= 0) top <- mx
      
      shinyWidgets::updateNumericRangeInput(session, inputId, value = c(0, round(top, 3)))
    }
    
    set_preset("sd_x", "sd_x_range")
    set_preset("sd_y", "sd_y_range")
    if ("sd_z" %in% names(sdr) && any(is.finite(sdr$sd_z))) set_preset("sd_z", "sd_z_range")
  }, ignoreInit = TRUE)
  has_sd_z <- reactive({
    sda <- sd_all_tbl()
    is.data.frame(sda) && "sd_z" %in% names(sda) && any(is.finite(sda$sd_z))
  })
  sd_all_tbl <- reactive({
    df0 <- uploaded_df()
    shiny::validate(shiny::need(!is.null(df0) && nrow(df0) > 0, "No data"))
    sd_by_tid(df0, min_n = input$sd_min_locs)            # RAW
  })
  
  sd_filtered_tbl <- reactive({
    dfF <- filtered_df()
    shiny::validate(shiny::need(!is.null(dfF) && nrow(dfF) > 0, "No data after current filters"))
    sd_by_tid(dfF, min_n = input$sd_min_locs)
  })

  # Effective spread display ranges: committed filter bounds with Inf replaced
  # by the actual data maximum.  Mirrors exactly what the spread sliders show,
  # even when the "Inspect trace spread" panel has never been opened.
  smlm_sd_display_ranges <- reactive({
    sda    <- tryCatch(sd_all_tbl(), error = function(e) NULL)
    mx_fn  <- function(v) { v <- v[is.finite(v)]; if (length(v)) max(v) else NA_real_ }
    eff    <- function(filt_val, mx) if (is.infinite(filt_val)) mx else filt_val
    mx_x   <- if (is.data.frame(sda)) mx_fn(sda$sd_x) else NA_real_
    mx_y   <- if (is.data.frame(sda)) mx_fn(sda$sd_y) else NA_real_
    mx_z   <- if (is.data.frame(sda) && "sd_z" %in% names(sda)) mx_fn(sda$sd_z) else NA_real_
    list(
      x = c(sd_filt$x[1], eff(sd_filt$x[2], mx_x)),
      y = c(sd_filt$y[1], eff(sd_filt$y[2], mx_y)),
      z = c(sd_filt$z[1], eff(sd_filt$z[2], mx_z))
    )
  })
  observe({
    for (id in c("sd_hist_x", "sd_hist_y", "sd_hist_z")) {
      outputOptions(output, id, suspendWhenHidden = FALSE)
    }
  })
  safe_update_sd_range_local <- function(session, inputId, value) {
    shinyWidgets::updateNumericRangeInput(session, inputId, value = value)
  }
  
  
  # Dynamic container that shows the right set of histograms (X,Y and Z if present)
  output$sd_hist_ui <- renderUI({
    req(input$view_SD)
    # Decide if Z is present from uploaded data
    df_u <- uploaded_df()
    shiny::validate(shiny::need(!is.null(df_u) && nrow(df_u) > 0, "No data uploaded"))
    has_z <- "z" %in% names(df_u) && any(is.finite(df_u$z))
    
    els <- list(
      plotlyOutput("sd_hist_x", height = 220),
      plotlyOutput("sd_hist_y", height = 220)
    )
    if (isTRUE(has_sd_z())) {
      els <- append(els, list(plotlyOutput("sd_hist_z", height = 220)))
    }
    do.call(tagList, els)
  })
  
  # ==========================================================================
  # Shared histogram ranges across sigma X / Y / Z - OPT IN
  # ==========================================================================
  # GOAL: put the three spread histograms on ONE set of bin edges and ONE count
  # axis, so the panels can be read against each other and dropped into a figure
  # as a set.
  #
  # WHY IT IS NOT THE DEFAULT. Per-axis ranges are the better view for reading a
  # single axis: MINFLUX z precision is genuinely different from xy, so a shared
  # range spends most of the x panel's width on empty space. Sharing is a
  # presentation choice, so it is opt-in and changes nothing else.
  #
  # The bin EDGES are shared too, not just the axis limits - histograms on
  # different edges are not comparable however the axes are labelled, because
  # the bar heights answer different questions.
  #
  # RETURNS NULL when off, so each renderer falls back to its own range.
  # ==========================================================================
  # Values outside a MANUAL range must be dropped, not drawn. sd_counts() clamps
  # anything beyond the outermost break into the edge bin, which is right for an
  # auto range (nothing is outside it) but would put every broad trace into one
  # tall bar at the right edge of a user-set range - a spike that is pure
  # artifact and would be read as a real population in a figure.
  sd_hist_clip <- function(v, sh) {
    if (is.null(sh) || !isTRUE(sh$manual)) return(v)
    v <- as.numeric(v)
    v[is.finite(v) & v >= sh$xrange[1] & v <= sh$xrange[2]]
  }

  # Median +/- normalized MAD, ONE LINE PER CHANNEL, each in its channel's colour.
  #
  # The bars were split by channel before this was, which made the panel say two
  # different things at once: coloured distributions with a single pooled median
  # drawn across them. That pooled median is not a property of either channel -
  # it lands between them, weighted by whichever channel happens to have more
  # traces, so on an unbalanced sample (which is normal: one target is usually
  # denser than the other) it sits close to the bigger channel and describes the
  # smaller one badly. It is also the number the spread FILTER is set from, and
  # that filter is shared across channels, so reading it per channel is what
  # tells you whether one threshold cuts them evenly.
  #
  # Falls back to the pooled line when there is only one channel.
  .sd_add_channel_medians <- function(fig, cur, colname, y_top, show_iqr) {
    cm <- tryCatch(sd_channel_map(), error = function(e) NULL)
    if (is.null(cm) || !"tid" %in% names(cur)) return(NULL)
    j <- dplyr::left_join(cur, cm, by = "tid")
    if (!"channel" %in% names(j) || all(is.na(j$channel))) return(NULL)
    chs <- sort(unique(j$channel[!is.na(j$channel)]))
    if (length(chs) < 2L) return(NULL)

    # Stagger the label heights so two nearby medians do not overprint. Channels
    # that differ by only a few percent are exactly the case worth reading.
    n <- length(chs)
    for (i in seq_along(chs)) {
      sub <- j[which(j$channel == chs[i]), , drop = FALSE]
      v   <- sub[[colname]]
      if (!length(v[is.finite(v)])) next
      col <- sub$channel_color[1]; if (is.na(col)) col <- "#666666"
      lab <- sigma_label_med_pm_sd(v)
      if (!is.null(lab)) lab <- paste0(sub$channel_label[1], "  ", lab)
      yl  <- y_top * (1 - 0.10 * (i - 1) / max(1, n - 1))
      fig <- add_sigma_vline(fig, stats::median(v, na.rm = TRUE), yl, col, lab)
      if (isTRUE(show_iqr)) fig <- add_sigma_iqr_lines(fig, v, yl, col)
    }
    fig
  }

  # ── TRACE SPREAD BY CHANNEL ───────────────────────────────────────────────
  # Trace spread is the ONE metadata distribution that is split by channel, and
  # deliberately the only one: EFO, CFR, DCR, ECO and dwell are properties of the
  # acquisition, shared by every channel in it, so splitting them would add noise
  # without adding information. Spread is different - it IS the localization
  # precision, precision follows photon yield, and photon yield is a property of
  # the dye. It is also what sets the de-blink merge distance, so seeing the two
  # distributions apart is what tells you whether one shared distance is honest.
  sd_channel_map <- reactive({
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(NULL)
    d <- tryCatch(attach_channel_colors(resolve_channels(df, input), input),
                  error = function(e) NULL)
    if (is.null(d) || !nrow(d) || !"channel" %in% names(d)) return(NULL)
    if (length(unique(d$channel[!is.na(d$channel)])) < 2L) return(NULL)
    d %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(channel = dplyr::first(channel),
                       channel_label = dplyr::first(channel_label),
                       channel_color = dplyr::first(channel_color),
                       .groups = "drop")
  })

  # One bar series per channel, in that channel's own colour, replacing the single
  # black "Filtered" series. Bars are stacked, so the silhouette still reads as the
  # same total distribution the black bars showed.
  .sd_add_channel_bars <- function(fig, cur, colname, sh, brk, A, in_rng) {
    cm <- tryCatch(sd_channel_map(), error = function(e) NULL)
    if (is.null(cm) || !"tid" %in% names(cur)) return(NULL)
    j <- dplyr::left_join(cur, cm, by = "tid")
    if (!"channel" %in% names(j) || all(is.na(j$channel))) return(NULL)
    chs <- sort(unique(j$channel[!is.na(j$channel)]))
    if (length(chs) < 2L) return(NULL)

    # STACKED, BUT DRAWN IN OVERLAY MODE - and the distinction is the whole fix.
    #
    # plotly's barmode is a LAYOUT property, so it applies to EVERY bar trace on
    # the panel. Asking for "stack" to stack the channels also stacked the grey
    # raw-data overlay, which put the grey in FRONT of the colours. The panel
    # needs both behaviours at once: grey behind, channels on top of each other.
    #
    # So the stack is built by hand. Draw the cumulative total in the last
    # channel's colour, then the running total one channel shorter over it, down
    # to the first channel. Each bar starts at zero and covers the taller one
    # behind it, so the visible bands are exactly the per-channel counts - a real
    # stack - while barmode stays "overlay" and the grey drawn earlier stays put.
    # Bars are opaque: a translucent stack would blend with the grey behind it and
    # neither colour would be the channel's own any more.
    per_ch <- lapply(chs, function(ch) {
      sub <- j[which(j$channel == ch), , drop = FALSE]
      cnt <- sd_counts(sd_hist_clip(sub[[colname]], sh), brk)
      cnt$counts[!in_rng] <- 0
      cc <- sub$channel_color[1]
      list(counts = cnt$counts,
           color  = if (is.na(cc)) "#666666" else cc,
           label  = sub$channel_label[1])
    })

    for (i in rev(seq_along(per_ch))) {
      cum <- Reduce(`+`, lapply(per_ch[seq_len(i)], `[[`, "counts"))
      fig <- fig %>% plotly::add_bars(
        x = A$centers, y = cum, width = A$width,
        marker = list(color = per_ch[[i]]$color, line = list(width = 0)),
        name = paste0(per_ch[[i]]$label, " (filtered)"), hoverinfo = "skip")
    }
    fig
  }

  sd_hist_common <- reactive({
    if (!isTRUE(input$sd_hist_shared)) return(NULL)
    all <- tryCatch(sd_all_tbl(),      error = function(e) NULL)
    cur <- tryCatch(sd_filtered_tbl(), error = function(e) NULL)
    if (!is.data.frame(all) || !nrow(all)) return(NULL)

    cols <- c("sd_x", "sd_y")
    if (isTRUE(has_sd_z()) && "sd_z" %in% names(all)) cols <- c(cols, "sd_z")
    cols <- intersect(cols, names(all))
    if (!length(cols)) return(NULL)

    bins <- as.integer(input$sd_bins %||% 60)

    # A user range wins when it is usable; otherwise fall back to the data's own
    # shared span, so clearing the boxes cannot leave the panels blank.
    ur     <- suppressWarnings(as.numeric(input$sd_hist_range))
    manual <- length(ur) == 2 && all(is.finite(ur)) && ur[2] > ur[1]

    brk <- if (manual) seq(ur[1], ur[2], length.out = max(2L, bins + 1L))
           else        sd_breaks(unlist(all[cols], use.names = FALSE), bins)
    if (is.null(brk) || length(brk) < 2) return(NULL)
    sh0 <- list(manual = manual, xrange = range(brk, na.rm = TRUE))

    # Counts have to be recomputed on the SHARED edges to find the common
    # ceiling - the per-axis maxima are on different edges and would not
    # transfer. Raw counts are included whether or not the grey overlay is on,
    # so toggling the overlay cannot rescale the shared axis under the user.
    ymax <- 0
    n_out <- 0L; n_tot <- 0L
    for (cc in cols) {
      v <- all[[cc]]; v <- v[is.finite(v)]
      n_tot <- n_tot + length(v)
      n_out <- n_out + sum(v < sh0$xrange[1] | v > sh0$xrange[2])
      ymax <- max(ymax, sd_counts(sd_hist_clip(all[[cc]], sh0), brk)$counts, na.rm = TRUE)
      if (is.data.frame(cur) && cc %in% names(cur))
        ymax <- max(ymax, sd_counts(sd_hist_clip(cur[[cc]], sh0), brk)$counts, na.rm = TRUE)
    }
    list(brk    = brk,
         manual = manual,
         xrange = sh0$xrange,
         yrange = c(0, if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1),
         n_out  = n_out, n_tot = n_tot)
  })

  # A figure must not hide data silently - say how much the chosen range leaves out.
  output$sd_hist_range_note <- renderUI({
    sh <- sd_hist_common()
    if (is.null(sh) || !isTRUE(sh$manual)) return(NULL)
    if (sh$n_out <= 0)
      return(helpText("All per-trace spreads fall inside this range."))
    helpText(HTML(sprintf(
      "<b>%s of %s per-trace spread values (%.1f%%) fall outside this range</b> and are not drawn. They are still in the data - this is a display range only, and it does not filter anything.",
      format(sh$n_out, big.mark = ","), format(sh$n_tot, big.mark = ","),
      100 * sh$n_out / max(1L, sh$n_tot))))
  })

  output$sd_hist_x <- renderPlotly({
    req(input$view_SD)
    all <- sd_all_tbl(); cur <- sd_filtered_tbl()
    shiny::validate(shiny::need(nrow(all) > 0, "No SD data"))

    bins <- as.integer(input$sd_bins %||% 60)
    sh   <- sd_hist_common()
    brk  <- if (!is.null(sh)) sh$brk else sd_breaks(all$sd_x, bins)
    A    <- sd_counts(sd_hist_clip(all$sd_x, sh), brk)
    B    <- sd_counts(sd_hist_clip(cur$sd_x, sh), brk)
    
    # active filter range
    r <- input$sd_x_range %||% c(-Inf, Inf)
    in_rng <- A$centers >= r[1] & A$centers <= r[2]
    
    # clamp filtered counts to 0 outside range (so black bars disappear there)
    B$counts[!in_rng] <- 0
    
    fig <- plotly::plot_ly()
    
    if (isTRUE(input$sd_overlay_all)) {
      fig <- fig %>% plotly::add_bars(
        x = A$centers, y = A$counts, width = A$width,
        # Drawn FIRST so it sits BEHIND the filtered bars, and lightly, so what
        # the filter removed reads as data disappearing rather than as a series
        # in its own right. Same idiom as the metadata histogram.
        marker = list(color = "rgba(160,160,160,0.35)"),
        name = "Raw data", hoverinfo = "skip"
      )
    }
    
    # Black foreground: filtered (B already clamped)
    # One series per channel when the data has channels, one black series when it
    # does not. Same totals either way.
    figc <- .sd_add_channel_bars(fig, cur, "sd_x", sh, brk, A, in_rng)
    fig <- if (!is.null(figc)) figc else fig %>% plotly::add_bars(
      x = A$centers, y = B$counts, width = A$width,
      marker = list(color = "rgba(0,0,0,0.85)"),
      name = "Filtered", hoverinfo = "skip"
    )
    ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
    # Under a shared axis the median label has to sit inside the SHARED
    # ceiling, or it is drawn off the top of the panel.
    y_top <- if (!is.null(sh)) sh$yrange[2]
             else if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
    
    # Filtered median ± nMAD (BLUE)
    if (isTRUE(input$sd_show_mean_sd)) {
      # One median per channel when the data has channels; the pooled median
      # otherwise. A pooled median drawn over channel-coloured bars would be a
      # number belonging to neither of them.
      figm <- .sd_add_channel_medians(fig, cur, "sd_x", y_top,
                                      show_iqr = TRUE)
      if (!is.null(figm)) {
        fig <- figm
      } else {
        muF  <- stats::median(cur$sd_x, na.rm = TRUE)
        labF <- sigma_label_med_pm_sd(cur$sd_x)
        fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
      # Quartiles flank the median only for the filtered set - drawing them for
      # the raw overlay too would put four thin lines on one panel.
        fig  <- add_sigma_iqr_lines(fig, cur$sd_x, y_top, "steelblue")
      }
    }

    # Raw median ± nMAD (TEAL, dashed)
    if (isTRUE(input$sd_show_overall_mean_sd)) {
      muA  <- stats::median(all$sd_x, na.rm = TRUE)
      labA <- sigma_label_med_pm_sd(all$sd_x)
      fig  <- add_sigma_vline(fig, muA, y_top, "cadetblue", labA, dash = "dash")
    }

    fig %>% plotly::layout(
      barmode = "overlay",
      title = "Trace spread \u03C3 X (nm)",
      xaxis = c(list(title = "\u03C3 X (nm)"), if (!is.null(sh)) list(range = sh$xrange)),
      yaxis = c(list(title = "Count"),         if (!is.null(sh)) list(range = sh$yrange))
    )
  })
  
  
  output$sd_hist_y <- renderPlotly({
    req(input$view_SD)
    all <- sd_all_tbl(); cur <- sd_filtered_tbl()
    shiny::validate(shiny::need(nrow(all) > 0, "No SD data"))
    
    bins <- as.integer(input$sd_bins %||% 60)
    sh   <- sd_hist_common()
    brk  <- if (!is.null(sh)) sh$brk else sd_breaks(all$sd_y, bins)
    A    <- sd_counts(sd_hist_clip(all$sd_y, sh), brk)
    B    <- sd_counts(sd_hist_clip(cur$sd_y, sh), brk)
    
    r <- input$sd_y_range %||% c(-Inf, Inf)
    in_rng <- A$centers >= r[1] & A$centers <= r[2]
    B$counts[!in_rng] <- 0
    
    fig <- plotly::plot_ly()
    
    if (isTRUE(input$sd_overlay_all)) {
      fig <- fig %>% plotly::add_bars(
        x = A$centers, y = A$counts, width = A$width,
        # Drawn FIRST so it sits BEHIND the filtered bars, and lightly, so what
        # the filter removed reads as data disappearing rather than as a series
        # in its own right. Same idiom as the metadata histogram.
        marker = list(color = "rgba(160,160,160,0.35)"),
        name = "Raw data", hoverinfo = "skip"
      )
    }
    
    # Black foreground: filtered (B already clamped)
    # One series per channel when the data has channels, one black series when it
    # does not. Same totals either way.
    figc <- .sd_add_channel_bars(fig, cur, "sd_y", sh, brk, A, in_rng)
    fig <- if (!is.null(figc)) figc else fig %>% plotly::add_bars(
      x = A$centers, y = B$counts, width = A$width,
      marker = list(color = "rgba(0,0,0,0.85)"),
      name = "Filtered", hoverinfo = "skip"
    )
    
    ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
    # Under a shared axis the median label has to sit inside the SHARED
    # ceiling, or it is drawn off the top of the panel.
    y_top <- if (!is.null(sh)) sh$yrange[2]
             else if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
    
    # Filtered median ± nMAD (BLUE)
    if (isTRUE(input$sd_show_mean_sd)) {
      # One median per channel when the data has channels; the pooled median
      # otherwise. A pooled median drawn over channel-coloured bars would be a
      # number belonging to neither of them.
      figm <- .sd_add_channel_medians(fig, cur, "sd_y", y_top,
                                      show_iqr = TRUE)
      if (!is.null(figm)) {
        fig <- figm
      } else {
        muF  <- stats::median(cur$sd_y, na.rm = TRUE)
        labF <- sigma_label_med_pm_sd(cur$sd_y)
        fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
        fig  <- add_sigma_iqr_lines(fig, cur$sd_y, y_top, "steelblue")
      }
    }

    # Raw median ± nMAD (TEAL, dashed)
    if (isTRUE(input$sd_show_overall_mean_sd)) {
      muA  <- stats::median(all$sd_y, na.rm = TRUE)
      labA <- sigma_label_med_pm_sd(all$sd_y)
      fig  <- add_sigma_vline(fig, muA, y_top, "cadetblue", labA, dash = "dash")
    }

    fig %>% plotly::layout(
      barmode = "overlay",
      title = "Trace spread \u03C3 Y (nm)",
      xaxis = c(list(title = "\u03C3 Y (nm)"), if (!is.null(sh)) list(range = sh$xrange)),
      yaxis = c(list(title = "Count"),         if (!is.null(sh)) list(range = sh$yrange))
    )
  })
  
  
  output$sd_hist_z <- renderPlotly({
    req(input$view_SD, has_sd_z())
    all <- sd_all_tbl(); cur <- sd_filtered_tbl()
    shiny::validate(shiny::need("sd_z" %in% names(all) && any(is.finite(all$sd_z)), "No Z in dataset"))
    
    bins <- as.integer(input$sd_bins %||% 60)
    sh   <- sd_hist_common()
    brk  <- if (!is.null(sh)) sh$brk else sd_breaks(all$sd_z, bins)
    A    <- sd_counts(sd_hist_clip(all$sd_z, sh), brk)
    B    <- sd_counts(sd_hist_clip(cur$sd_z, sh), brk)
    
    r <- input$sd_z_range %||% c(-Inf, Inf)
    in_rng <- A$centers >= r[1] & A$centers <= r[2]
    B$counts[!in_rng] <- 0
    
    fig <- plotly::plot_ly()
    
    if (isTRUE(input$sd_overlay_all)) {
      fig <- fig %>% plotly::add_bars(
        x = A$centers, y = A$counts, width = A$width,
        # Drawn FIRST so it sits BEHIND the filtered bars, and lightly, so what
        # the filter removed reads as data disappearing rather than as a series
        # in its own right. Same idiom as the metadata histogram.
        marker = list(color = "rgba(160,160,160,0.35)"),
        name = "Raw data", hoverinfo = "skip"
      )
    }
    
    # Black foreground: filtered (B already clamped)
    # One series per channel when the data has channels, one black series when it
    # does not. Same totals either way.
    figc <- .sd_add_channel_bars(fig, cur, "sd_z", sh, brk, A, in_rng)
    fig <- if (!is.null(figc)) figc else fig %>% plotly::add_bars(
      x = A$centers, y = B$counts, width = A$width,
      marker = list(color = "rgba(0,0,0,0.85)"),
      name = "Filtered", hoverinfo = "skip"
    )
    
    ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
    # Under a shared axis the median label has to sit inside the SHARED
    # ceiling, or it is drawn off the top of the panel.
    y_top <- if (!is.null(sh)) sh$yrange[2]
             else if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
    
    # Filtered median ± nMAD (BLUE)
    if (isTRUE(input$sd_show_mean_sd)) {
      # One median per channel when the data has channels; the pooled median
      # otherwise. A pooled median drawn over channel-coloured bars would be a
      # number belonging to neither of them.
      figm <- .sd_add_channel_medians(fig, cur, "sd_z", y_top,
                                      show_iqr = TRUE)
      if (!is.null(figm)) {
        fig <- figm
      } else {
        muF  <- stats::median(cur$sd_z, na.rm = TRUE)
        labF <- sigma_label_med_pm_sd(cur$sd_z)
        fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
        fig  <- add_sigma_iqr_lines(fig, cur$sd_z, y_top, "steelblue")
      }
    }

    # Raw median ± nMAD (TEAL, dashed)
    if (isTRUE(input$sd_show_overall_mean_sd)) {
      muA  <- stats::median(all$sd_z, na.rm = TRUE)
      labA <- sigma_label_med_pm_sd(all$sd_z)
      fig  <- add_sigma_vline(fig, muA, y_top, "cadetblue", labA, dash = "dash")
    }

    fig %>% plotly::layout(
      barmode = "overlay",
      title = "Trace spread \u03C3 Z (nm)",
      xaxis = c(list(title = "\u03C3 Z (nm)"), if (!is.null(sh)) list(range = sh$xrange)),
      yaxis = c(list(title = "Count"),         if (!is.null(sh)) list(range = sh$yrange))
    )
  })

  .sanitize_range <- function(r, v) {
    v <- v[is.finite(v)]
    if (is.null(r) || length(r) != 2 || any(!is.finite(r))) {
      if (length(v)) {
        # default to [min, 95th percentile] so you don't nuke the data
        c(min(v), stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE))
      } else {
        c(-Inf, Inf)
      }
    } else r
  }
  
  # --------- 8) TRIGGERS & STORES ----------
  plot_trigger <- reactiveVal(NULL)
  smlm_diag      <- reactiveVal(NULL)   # list(data_clusters, img_blobs) in data coords
  smlm_diag_mats <- reactiveVal(NULL)  # list(raw=mat_c, smooth=mat_cs) for debug plot
  # One row per isosurface drawn in the last 3-D render, so the controls can
  # report the geometry they picked. See iso_status_ui().
  cluster2_iso_diag <- reactiveVal(NULL)

  # ==========================================================================
  # ISOSURFACE: COMMITTED GEOMETRY + MESH CACHE
  # ==========================================================================
  # Wrapping every Median Cluster is the most expensive thing the 3-D view does,
  # and in 2-colour mode it runs over both channels. Two separate problems made
  # that painful, and they need different fixes:
  #
  #   1. DRAGGING A SLIDER re-ran the whole wrap on every intermediate value. The
  #      geometry sliders are therefore COMMITTED, exactly like the clustering
  #      parameters: they are snapshotted when Plot Data is pressed and the view
  #      ignores them until then.
  #   2. ANY re-render rebuilt the meshes, because a renderPlotly body re-runs on
  #      all of its reactive dependencies - so changing point size or colour mode,
  #      which have nothing to do with the surfaces, still paid for them. The
  #      cache keys each mesh on its geometry alone, so those re-renders reuse it.
  #
  # The split follows the control panel: Shape is committed, Appearance is live.
  # Opacity and colour stay instant because they are cache hits, not rebuilds.
  iso_mesh_cache <- new.env(parent = emptyenv())

  iso_geom <- eventReactive(plot_trigger(), {
    # A new commit invalidates every cached mesh. Clearing here rather than on a
    # timer keeps the cache bounded to one settings epoch.
    rm(list = ls(iso_mesh_cache), envir = iso_mesh_cache)
    list(
      c2_grid_res     = input$cluster2_iso_grid_res,
      c2_pad          = input$cluster2_iso_pad,
      c2_smooth_rel   = input$cluster2_iso_smooth_rel,
      c2_manual_off   = isTRUE(input$cluster2_iso_manual_offset),
      c2_offset       = input$cluster2_iso_offset,
      c2_manual_sigma = isTRUE(input$cluster2_iso_manual_sigma),
      c2_sigma        = input$cluster2_iso_sigma
    )
  })
  output$cluster2_iso_status <- renderUI(iso_status_ui(cluster2_iso_diag()))
  cluster2_table_store <- reactiveVal(NULL)
  # Cache of clustered data per file - populated whenever a file is plotted.
  # Enables batch analysis across all previously-plotted files.
  smlm_plotted_store <- reactiveVal(list())
  shinyjs::hide("scatter2D")
  observeEvent(plot_trigger(), {
    shinyjs::show("scatter2D")
  }, ignoreNULL = TRUE, once = TRUE)
  observeEvent(input$plot_data, { plot_trigger(Sys.time()) })

  output$plot_bar_file_label <- renderUI({
    nm <- input$active_file %||% ""
    if (!nzchar(nm)) return(NULL)
    tags$span(style = "font-size: 13px; color: #555;", nm)
  })

  # view_plane_2d is now a static input (ui_smlm.R). The old renderUI dropped
  # XZ/YZ from the choices when `any(is.finite(df$z))` was FALSE - which never
  # fired for a 2D MINFLUX file, because those carry z = 0 and 0 IS finite. So
  # side views were offered on planar data and produced a flat line. The planar
  # observer below disables them on z VARIATION instead.
  observeEvent(input$enable_clustering, {
    req(input$file)
    plot_trigger(Sys.time())
  }, ignoreInit = TRUE)
  fit_log_rv <- reactiveVal("")

  # Keep smlm_plotted_store in sync: every time a file is plotted, cache its
  # clustered data so batch analysis can iterate all plotted files.
  observeEvent(plottable_data(), {
    df  <- tryCatch(plottable_data(), error = function(e) NULL)
    nm  <- isolate(input$active_file) %||% isolate(input$file)[[1]]
    if (is.null(df) || !nrow(df) || is.null(nm)) return()
    store <- smlm_plotted_store()
    store[[nm]] <- df
    smlm_plotted_store(store)
  }, ignoreNULL = TRUE)


  output$fit_log <- renderText({
    fit_log_rv()
  })
  
  # clear the box whenever metrics are disabled
  observeEvent(
    list(input$show_cluster2_metrics, input$enable_cluster2),
    {
      if (!isTRUE(input$enable_cluster2) || !isTRUE(input$show_cluster2_metrics)) {
        fit_log_rv("")
      }
    },
    ignoreInit = TRUE
  )
  observeEvent(
    list(input$remove_rogues, input$view_SD),
    {
      if (!isTRUE(input$remove_rogues) || !isTRUE(input$view_SD)) {
        # only clear if the current message is the rogue-removal one
        cur <- fit_log_rv()
        if (is.character(cur) && length(cur) == 1 && grepl("^Rogue removal:", cur)) {
          fit_log_rv("")
        }
      }
    },
    ignoreInit = TRUE
  )
  observeEvent(input$enable_two_color, {
    if (isTRUE(input$enable_two_color) && !is.null(input$color_mode) && input$color_mode != "default") {
      updateSelectInput(session, "color_mode", selected = "default")
      showNotification(
        "2-color mode uses Default coloring. If you switch to a numeric color mode (z/time/density), that will override the 2-color colors.",
        type = "message", duration = 6
      )
    }
  }, ignoreInit = TRUE)

  # Show "Use Pooled DCR" checkbox only when dcr_pooled is present in the data
  output$use_pooled_dcr_ui <- renderUI({
    df <- uploaded_df_rv()
    if (!is.null(df) && "dcr_pooled" %in% names(df)) {
      checkboxInput("use_pooled_dcr",
                    "Use Pooled DCR for color splitting",
                    value = FALSE)
    }
  })

  # ------ 2-color: auto-detect DCR split via GMM
  dcr_gmm_result <- reactiveVal(NULL)

  observeEvent(input$auto_detect_dcr, {
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    has_dcr <- "dcr" %in% names(df) || "dcr_pooled" %in% names(df)
    if (is.null(df) || !has_dcr || !"tid" %in% names(df)) {
      showNotification("No DCR data available for auto-detection.", type = "warning", duration = 5)
      dcr_gmm_result(list(success = FALSE, reason = "No DCR column in data"))
      return()
    }
    # Use dcr_pooled when requested and available; otherwise per-trace median of raw dcr
    .col <- if (isTRUE(input$use_pooled_dcr) && "dcr_pooled" %in% names(df)) "dcr_pooled" else "dcr"
    dcr_meds <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(dcr_med = stats::median(.data[[.col]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::pull(dcr_med)

    res <- auto_detect_dcr_split(dcr_meds, col_used = .col)
    dcr_gmm_result(res)

    if (!isTRUE(res$success)) {
      showNotification(paste("GMM auto-detect failed:", res$reason), type = "error", duration = 7)
      return()
    }

    split <- res$split
    shinyWidgets::updateNumericRangeInput(session, "dcr_group1", value = c(0, split))
    shinyWidgets::updateNumericRangeInput(session, "dcr_group2", value = c(split + 1e-4, 1))
    showNotification(
      paste0("Split set at DCR = ", split,
             "  |  G1 ", res$mean1, " (SD ", res$sd1, ")",
             "  |  G2 ", res$mean2, " (SD ", res$sd2, ")",
             "  |  d' = ", res$dprime, " (", res$quality, "), ",
             sprintf("%.1f%% expected crosstalk", 100 * res$err_rate),
             if (!identical(.col, "dcr_pooled"))
               "  |  using RAW dcr - try pooled DCR for a cleaner split" else ""),
      type = if (identical(res$quality, "poor")) "warning" else "message", duration = 9
    )
  })

  # Per-trace split confidence: which individual molecules sit on the edge of the
  # DCR split. The GMM summary only gives an aggregate ambiguous %, which cannot
  # tell you WHICH molecules are uncertain - and those are exactly the ones that
  # manufacture false short-range cross-channel pairs.
  dcr_conf_tbl <- reactive({
    gmm <- dcr_gmm_result()
    if (is.null(gmm) || !isTRUE(gmm$success)) return(NULL)
    df <- tryCatch(filtered_df(), error = function(e) NULL)
    if (is.null(df) || !"tid" %in% names(df)) return(NULL)
    .col <- if (identical(gmm$col_used, "dcr_pooled") && "dcr_pooled" %in% names(df))
      "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL
    if (is.null(.col)) return(NULL)

    med <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(dcr_med = stats::median(.data[[.col]], na.rm = TRUE),
                       n_locs  = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(is.finite(dcr_med))
    if (!nrow(med)) return(NULL)

    cf <- dcr_split_confidence(med$dcr_med, gmm)
    med$group      <- cf$group
    med$posterior  <- round(cf$post, 4)
    med$log_odds   <- round(cf$log_odds, 3)
    med$on_edge    <- cf$edge
    med$dcr_column <- .col
    med[order(med$posterior), , drop = FALSE]   # least certain first
  })

  output$dcr_conf_plot <- renderPlotly({
    tb  <- dcr_conf_tbl()
    gmm <- dcr_gmm_result()
    shiny::validate(shiny::need(!is.null(tb) && nrow(tb) > 0,
      "Run 'Auto-detect split (GMM)' to score each trace's assignment confidence."))
    cols <- ifelse(tb$on_edge, "#e67e22", ifelse(tb$group == 1L, "#1f77b4", "#d62728"))
    plotly::plot_ly(
      x = tb$dcr_med, y = tb$posterior, type = "scatter", mode = "markers",
      marker = list(color = cols, size = 6, opacity = 0.75,
                    line = list(color = "white", width = 0.4)),
      text = paste0("tid ", tb$tid, "<br>DCR ", round(tb$dcr_med, 4),
                    "<br>group ", tb$group, "<br>posterior ", round(tb$posterior, 3),
                    ifelse(tb$on_edge, "<br><b>on edge</b>", "")),
      hoverinfo = "text") %>%
      plotly::layout(
        title = list(text = sprintf(
          "Per-trace split confidence - %d of %d traces on the edge (posterior < %.2f), orange",
          sum(tb$on_edge, na.rm = TRUE), nrow(tb), gmm$post_cut), font = list(size = 12)),
        margin = list(t = 46),
        xaxis = list(title = paste0("Trace median ", gmm$col_used)),
        yaxis = list(title = "Posterior of assigned group", range = c(0.45, 1.02)),
        shapes = list(
          list(type = "line", x0 = gmm$split, x1 = gmm$split, y0 = 0.45, y1 = 1.02,
               line = list(color = "#666666", width = 1.5, dash = "dash")),
          list(type = "line", x0 = min(tb$dcr_med), x1 = max(tb$dcr_med),
               y0 = gmm$post_cut, y1 = gmm$post_cut,
               line = list(color = "#e67e22", width = 1, dash = "dot"))),
        showlegend = FALSE)
  })

  output$download_dcr_conf <- downloadHandler(
    filename = function() paste0("dcr_split_confidence_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      tb <- dcr_conf_tbl()
      shiny::validate(shiny::need(!is.null(tb) && nrow(tb) > 0, "Run the GMM auto-detect first."))
      write.csv(tb, file, row.names = FALSE)
    }
  )

  # Split summary + an explicit split-QUALITY verdict. The crosstalk numbers are
  # what a cross-channel co-localization result has to be read against, so they are
  # shown next to the split rather than buried: an apparent co-localization signal
  # at short r is only meaningful if it is large compared with the contamination.
  output$dcr_gmm_info <- renderUI({
    res <- dcr_gmm_result()
    if (is.null(res)) return(NULL)
    if (!isTRUE(res$success))
      return(tags$div(style = "font-size:0.85em; color:#c0392b;",
                      paste("GMM failed:", res$reason)))

    col_lbl  <- if (identical(res$col_used, "dcr_pooled")) "pooled DCR" else "raw DCR"
    pooled   <- identical(res$col_used, "dcr_pooled")
    has_pooled_col <- {
      d <- tryCatch(filtered_df(), error = function(e) NULL)
      !is.null(d) && "dcr_pooled" %in% names(d)
    }
    qcol <- switch(res$quality, clean = "#1a8a1a", workable = "#b26a00", "#c0392b")
    qtxt <- switch(res$quality,
      clean    = "clean split",
      workable = "workable split",
      poor     = "POOR split - channels overlap heavily",
      "separation unknown")

    tagList(
      tags$div(style = "font-size:0.85em; margin-top:4px;",
        tags$b("Split: "), res$split, " on ", col_lbl,
        tags$br(),
        "G1 mean ", res$mean1, " (SD ", res$sd1, ", ", round(100 * res$pi1), "%)",
        " | G2 mean ", res$mean2, " (SD ", res$sd2, ", ", round(100 * res$pi2), "%)",
        " | n = ", res$n_traces, " traces"),
      tags$div(style = paste0("font-size:0.85em; margin-top:3px; color:", qcol, "; font-weight:600;"),
        "d' = ", res$dprime, " - ", qtxt),
      tags$div(style = "font-size:0.82em; color:#555; margin-top:2px;",
        tags$b("Expected crosstalk: "),
        sprintf("%.1f%% of traces mis-assigned overall; group 1 is ~%.1f%% contaminated by dye 2, group 2 ~%.1f%% by dye 1. %.1f%% of traces are ambiguous (max posterior < %.2f).",
                100 * res$err_rate, 100 * res$contam1, 100 * res$contam2,
                100 * res$ambiguous_frac, res$post_cut)),
      tags$div(style = "font-size:0.8em; color:#7a5c00; margin-top:3px;",
        HTML(paste0(
          "Crosstalk pushes cross-channel NN / L₁₂ / g₁₂ <b>up at short r</b>, toward false ",
          "co-localization. Treat a co-localization signal as real only if it is large ",
          "compared with these numbers.",
          # Only suggest pooling where pooling is actually possible: it happens in the
          # .npy/.json readers. A .csv is a pyMINFLUX export of the final iteration and
          # has already been through their pipeline, so there is nothing to pool here.
          if (!pooled && has_pooled_col)
            paste0(" <b>This file has a pooled DCR column you are not using</b> - tick ",
                   "<i>Use pooled DCR</i> and re-upload; pooled DCR usually separates the dyes better.")
          else if (!pooled)
            paste0(" Splitting on raw DCR. If this came from <b>.npy/.json</b>, re-uploading with ",
                   "<i>Pool DCR</i> enabled (photon-weighted across the iteration sequence) usually ",
                   "sharpens the split; a <b>.csv</b> is already a pyMINFLUX final-iteration export, ",
                   "so there is nothing further to pool.")
          else ""))),
      # Keep this proportionate: the GMM is a reasonable estimator and the split is
      # what the data allows. The point is only that split QUALITY is not evidence of
      # two dyes existing - the experiment design is.
      tags$div(style = "margin-top:5px; padding:5px 8px; border-left:3px solid #b26a00; background:rgba(178,106,0,0.06); font-size:0.8em; color:#7f8c8d; line-height:1.45;",
        HTML(paste0(
          "<b>d' measures separation, not whether two dyes exist.</b>",
          " A two-component GMM will always produce 2 groups, even if there is only one dye."))))
    })

  # ------ HISTOGRAM for filters
  output$metadata_histogram <- renderPlotly({
    req(input$show_metadata_hist, input$hist_variable, uploaded_df_rv())
    
    varname <- input$hist_variable
    
    get_vals <- function(df, var) {
      if (var == "tid") {
        tlens <- as.numeric(table(df$tid))
        tlens[is.finite(tlens)]
      } else {
        v <- df[[var]]
        v[is.finite(v)]
      }
    }
    
    df_full <- uploaded_df()
    
    df_keep <- filtered_df()
    
    vals_full <- get_vals(df_full, varname)
    vals_keep <- get_vals(df_keep, varname)
    
    if (!is.null(input$hist_range)) {
      lo <- input$hist_range[1]; hi <- input$hist_range[2]
      vals_full <- vals_full[vals_full >= lo & vals_full <= hi]
      vals_keep <- vals_keep[vals_keep >= lo & vals_keep <= hi]
    }
    
    if (length(vals_full) == 0 && length(vals_keep) == 0) return(NULL)
    
    bg_col   <- if (isTRUE(input$dark_mode)) "#1e1e1e" else "white"
    font_col <- if (isTRUE(input$dark_mode)) "white" else "black"
    
    fig <- plot_ly()
    
    # Optional raw overlay (grey)
    if (isTRUE(input$hist_show_raw)) {
      fig <- fig %>%
        add_histogram(
          x = vals_full, nbinsx = input$hist_bins, name = "Raw data",
          opacity = 0.35, marker = list(color = "lightgray"), bingroup = "x"
        )
    }
    
    # Always show the "filtered" bars (blue)
    fig %>%
      add_histogram(
        x = vals_keep, nbinsx = input$hist_bins, name = "Filtered data",
        opacity = 0.85, marker = list(color = "steelblue"), bingroup = "x"
      ) %>%
      layout(
        barmode = "overlay",
        title = paste("Histogram of", if (varname == "tid") "Trace lengths" else varname),
        xaxis = list(title = if (varname == "tid") "Trace length (# per TID)" else varname),
        yaxis = list(title = "Count"),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col, font = list(color = font_col),
        legend = list(orientation = "h", x = 0, y = 1.1)
      )
  })
  #2-color NN helpers
  # One-row-per-entity table: either per TID or per cluster, with optional 2-color grouping.
  # - If z doesn't exist, we use z=0.
  # - If enable_clustering is TRUE, we use cluster medians (cluster != 0) and propagate group by majority vote.
  build_geo_medians <- function(df, input, use_two_color = FALSE) {
    if (!"z" %in% names(df)) df$z <- 0
    
    # Always guarantee the column exists (even if we can't assign groups)
    if (!"analysis_group" %in% names(df)) df$analysis_group <- NA_integer_
    
    if (isTRUE(use_two_color)) {
      # assign_two_color_group() returns one row per tid with analysis_group
      gmap <- tryCatch(assign_two_color_group(df, input), error = function(e) NULL)
      
      if (is.data.frame(gmap) && nrow(gmap) && "analysis_group" %in% names(gmap)) {
        df <- df %>%
          dplyr::left_join(gmap %>% dplyr::select(tid, analysis_group), by = "tid")
      }
      
      # If join didn’t happen for any reason, still keep the column
      if (!"analysis_group" %in% names(df)) df$analysis_group <- NA_integer_
      
    } else {
      # your general grouping (may set NA_integer_ when disabled)
      df2 <- tryCatch(assign_analysis_group(df, input), error = function(e) df)
      df <- df2
      if (!"analysis_group" %in% names(df)) df$analysis_group <- NA_integer_
    }
    
    if (isTRUE(input$enable_clustering) && "cluster" %in% names(df)) {
      geo <- df %>%
        dplyr::filter(!is.na(.data$cluster) & .data$cluster != 0) %>%
        dplyr::group_by(.data$cluster) %>%
        dplyr::summarise(
          id = dplyr::first(.data$cluster),
          x  = stats::median(.data$x, na.rm = TRUE),
          y  = stats::median(.data$y, na.rm = TRUE),
          z  = stats::median(.data$z, na.rm = TRUE),
          analysis_group = {
            v <- .data$analysis_group[!is.na(.data$analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
          .groups = "drop"
        ) %>%
        dplyr::mutate(id_type = "cluster")
    } else {
      geo <- df %>%
        dplyr::group_by(.data$tid) %>%
        dplyr::summarise(
          id = dplyr::first(.data$tid),
          x  = stats::median(.data$x, na.rm = TRUE),
          y  = stats::median(.data$y, na.rm = TRUE),
          z  = stats::median(.data$z, na.rm = TRUE),
          analysis_group = {
            v <- .data$analysis_group[!is.na(.data$analysis_group)]
            if (length(v)) as.integer(round(stats::median(v))) else NA_integer_
          },
          .groups = "drop"
        ) %>%
        dplyr::mutate(id_type = "tid")
    }
    
    geo
  }
  
  
  safe_nn2_local <- function(coords, k) {
    # Use your existing one if you already have it
    out <- tryCatch({
      if (!is.matrix(coords) || nrow(coords) < 1) return(list(nn = NULL, k_eff = 0))
      k_eff <- max(1, min(k, nrow(coords)))
      nn <- RANN::nn2(coords, k = k_eff)
      list(nn = nn, k_eff = k_eff)
    }, error = function(e) list(nn = NULL, k_eff = 0))
    out
  }
  
  # Shift midpoint of segment (p1→p2) perpendicularly in XY plane by offset_nm.
  # Works for both 2-element (x,y) and 3-element (x,y,z) vectors.
  perp_offset_xy <- function(p1, p2, offset_nm) {
    p1 <- as.numeric(p1); p2 <- as.numeric(p2)
    mid <- (p1 + p2) / 2
    if (!is.finite(offset_nm) || offset_nm <= 0) return(mid)
    dx <- p2[1] - p1[1]; dy <- p2[2] - p1[2]
    len <- sqrt(dx^2 + dy^2)
    if (len < 1e-10) return(mid)
    perp <- c(-dy / len, dx / len)
    if (length(mid) == 3) perp <- c(perp, 0)
    mid + perp * offset_nm
  }

  # NN index pairs honouring the channel scope radio. Reads `channel`, so one
  # control covers a 2-colour split and a multiplex build alike; single-channel
  # data falls back to a global k-NN. "within" = same-channel only, "across" =
  # strictly between different channels (never within), "all" = both.
  # Used by every NN-line / NN-table path.
  .nn_pairs_scoped <- function(df_geo, coords, input) {
    k <- input$nn_k %||% 4
    ch <- if ("channel" %in% names(df_geo)) df_geo$channel
          else if ("mplex_channel" %in% names(df_geo)) df_geo$mplex_channel else NULL
    if (is.null(ch) || length(unique(ch[!is.na(ch)])) < 2L)
      return(nn_index_pairs(coords, k))
    scope <- input$nn_mplex_scope %||% "within"
    if (identical(scope, "all")) return(nn_index_pairs(coords, k))
    nn_index_pairs(coords, k, channel = ch,
                   mode = if (identical(scope, "across")) "between" else "within")
  }

  # core builder used by both 2D/3D
  build_nn_segments <- function(df, input, is3d = TRUE) {
    if (!nrow(df)) return(NULL)
    
    # If analysis_group isn't present, create NA column (scope=="all" will still work)
    if (!"analysis_group" %in% names(df)) df$analysis_group <- NA_integer_
    
    scope <- "all"
    
    # 1) Choose medians: by cluster if enabled, else by tid
    if (isTRUE(input$enable_clustering) && "cluster" %in% names(df)) {
      df_geo <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
          color = dplyr::first(color),
          cluster2 = if ("cluster2" %in% names(df)) dplyr::first(cluster2) else NA_integer_,
          analysis_group = suppressWarnings(dplyr::first(analysis_group)),
          channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
          channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
          channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups = "drop"
        )
    } else {
      df_geo <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
          color = dplyr::first(color),
          cluster2 = if ("cluster2" %in% names(df)) dplyr::first(cluster2) else NA_integer_,
          analysis_group = suppressWarnings(dplyr::first(analysis_group)),
          channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
          channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
          channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups = "drop"
        )
    }

    if (nrow(df_geo) < 2) return(NULL)

    coords <- as.matrix(df_geo[, c("x","y","z")])
    pm <- .nn_pairs_scoped(df_geo, coords, input)
    if (!nrow(pm)) return(NULL)

    seg_x <- c(); seg_y <- c(); seg_z <- c()
    lab_x <- c(); lab_y <- c(); lab_z <- c(); labs <- c()
    e_i <- c(); e_j <- c(); e_d <- c()      # per-edge record, for channel colouring
    seen  <- new.env(parent = emptyenv())

    dmin <- (input$distance_range %||% c(0, Inf))[1]
    dmax <- (input$distance_range %||% c(0, Inf))[2]

    for (r in seq_len(nrow(pm))) {
      i <- pm[r, 1]; ni <- pm[r, 2]
      {
        # dedupe undirected
        key <- paste0(min(i, ni), "_", max(i, ni))
        if (!is.null(seen[[key]])) next
        seen[[key]] <- TRUE

        # within/between/all
        g1 <- df_geo$analysis_group[i]
        g2 <- df_geo$analysis_group[ni]
        if (!pair_passes_scope(g1, g2, scope)) next

        dist <- pm[r, 3]
        if (dist < dmin || dist > dmax) next

        # within-cluster2 filter
        if (isTRUE(input$nn_within_cluster2) && isTRUE(input$enable_cluster2) && "cluster2" %in% names(df_geo)) {
          c2i <- df_geo$cluster2[i]; c2ni <- df_geo$cluster2[ni]
          if (is.na(c2i) || is.na(c2ni) || c2i == 0 || c2ni == 0 || c2i != c2ni) next
        }

        seg_x <- c(seg_x, coords[i,1], coords[ni,1], NA_real_)
        seg_y <- c(seg_y, coords[i,2], coords[ni,2], NA_real_)
        seg_z <- c(seg_z, coords[i,3], coords[ni,3], NA_real_)

        nn_off <- suppressWarnings(as.numeric(input$nn_label_offset %||% 5))
        if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
        mid <- perp_offset_xy(coords[i,], coords[ni,], nn_off)
        lab_x <- c(lab_x, mid[1]); lab_y <- c(lab_y, mid[2]); lab_z <- c(lab_z, mid[3])
        labs  <- c(labs, round(dist, 2))
        e_i <- c(e_i, i); e_j <- c(e_j, ni); e_d <- c(e_d, dist)
      }
    }

    if (!length(seg_x)) return(NULL)

    # Per-edge table, so the plot can colour a line by the channels it joins.
    # This is what the separate 2-colour NN panel used to exist for: a
    # within-channel line drawn in that channel's colour, a cross-channel line in
    # magenta. It is the one thing that panel did which the shared one could not,
    # so it moved here rather than being dropped with it.
    edges <- data.frame(
      x1 = coords[e_i, 1], y1 = coords[e_i, 2], z1 = coords[e_i, 3],
      x2 = coords[e_j, 1], y2 = coords[e_j, 2], z2 = coords[e_j, 3],
      dist = e_d,
      ch1 = df_geo$channel[e_i], ch2 = df_geo$channel[e_j],
      lab1 = df_geo$channel_label[e_i], lab2 = df_geo$channel_label[e_j],
      stringsAsFactors = FALSE
    )
    same <- !is.na(edges$ch1) & !is.na(edges$ch2) & edges$ch1 == edges$ch2
    col_same <- df_geo$channel_color[e_i]
    edges$color <- ifelse(same,
                          ifelse(is.na(col_same), "#666666", col_same),
                          "magenta")

    list(
      seg_x = seg_x, seg_y = seg_y, seg_z = seg_z,
      lab_x = lab_x, lab_y = lab_y, lab_z = lab_z, labs = labs,
      edges = edges
    )
  }
  
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # Always return a correctly-shaped empty edge table
  empty_edges <- function() {
    data.frame(
      x1 = numeric(0), y1 = numeric(0), z1 = numeric(0),
      x2 = numeric(0), y2 = numeric(0), z2 = numeric(0),
      dist = numeric(0),
      g1 = integer(0), g2 = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  # 1) Median coordinates depending on "enable_clustering"
  #    - If clustering = TRUE, use cluster medians (cluster != 0)
  #    - Else, use trace medians by tid
  #    Returns: data.frame(group_id, x,y,z?, optional aux columns you might add upstream)
  nn_medians <- function(df, enable_clustering = FALSE) {
    if (enable_clustering) {
      out <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          analysis_group = suppressWarnings(dplyr::first(analysis_group)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(key_id = paste0("C", cluster))
    } else {
      out <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          analysis_group = suppressWarnings(dplyr::first(analysis_group)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(key_id = paste0("T", tid))
    }
    out
  }
  
  
  # 2) Attach analysis_group (1/2) to the *median table*
  #    Uses your existing assign_analysis_group(df, input) which tags raw rows.
  #    We then propagate group to the medians (by tid or cluster mapping).
  attach_group_to_medians <- function(df_med, df_raw, input, enable_clustering = FALSE) {
    if (!nrow(df_med)) return(dplyr::mutate(df_med, analysis_group = NA_integer_))
    
    tagged <- assign_analysis_group(df_raw, input)
    # Build a map from tid or cluster to analysis_group at the *trace level*.
    # For cluster medians, we derive the cluster->group via the majority vote of member traces.
    if (enable_clustering && "cluster" %in% names(tagged)) {
      # majority group per cluster (using tid-level group assignment)
      tid_groups <- tagged %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(analysis_group = suppressWarnings(as.integer(round(stats::median(analysis_group, na.rm = TRUE)))),
                         cluster = suppressWarnings(as.integer(round(stats::median(cluster, na.rm = TRUE)))),
                         .groups = "drop")
      cl_map <- tid_groups %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          analysis_group = {
            v <- analysis_group[!is.na(analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
          .groups = "drop"
        ) %>%
        dplyr::rename(group_id = cluster)
      
      out <- df_med %>%
        dplyr::left_join(cl_map, by = "group_id")
    } else {
      # simple map: tid -> group (median over rows for that tid)
      tid_map <- tagged %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          analysis_group = suppressWarnings(as.integer(round(stats::median(analysis_group, na.rm = TRUE)))),
          .groups = "drop"
        ) %>%
        dplyr::rename(group_id = tid)
      
      out <- df_med %>%
        dplyr::left_join(tid_map, by = "group_id")
    }
    
    if (!"analysis_group" %in% names(out)) out$analysis_group <- NA_integer_
    out
  }
  
  
  #--- 2D-Only server parts:
  
  # --- 2D helpers --------------------------------------------------------------
  
  # Build a dataURI for the uploaded background image and optionally rotate (0/90/180/270)
  
  
  
  # ---- NN Table for Cluster2 (2D) ----
  nn_table_cluster2 <- reactive({
    df <- plottable_data()
    req(isTRUE(input$enable_cluster2))
    
    df_geo <- df %>%
      dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
      dplyr::group_by(cluster2) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        color = dplyr::first(color),    
        .groups = "drop"
      )
    if (nrow(df_geo) < 2) return(df_geo[0, ])
    
    coords <- as.matrix(df_geo[, c("x","y")])
    k_eff  <- min(input$nn_k_cluster2 + 1, nrow(coords))
    nn     <- RANN::nn2(coords, k = k_eff)
    
    rows <- list()
    for (i in 1:nrow(coords)) {
      for (j in 2:k_eff) {
        ni   <- nn$nn.idx[i, j]
        dist <- as.numeric(nn$nn.dists[i, j])
        if (dist >= input$distance_range_cluster2[1] && dist <= input$distance_range_cluster2[2]) {
          rows[[length(rows) + 1]] <- data.frame(
            source_cluster2 = df_geo$cluster2[i],
            target_cluster2 = df_geo$cluster2[ni],
            distance = round(dist, 2),
            source_x = df_geo$x[i], source_y = df_geo$y[i],
            target_x = df_geo$x[ni], target_y = df_geo$y[ni]
          )
        }
      }
    }
    if (length(rows)) dplyr::bind_rows(rows) else df_geo[0, ]
  })
  
  compute_kde <- function(df, lim, n = 200, bw = NA_real_, bw_factor = 1) {
    if (!all(c("x","y") %in% names(df)) || nrow(df) < 2) return(NULL)
    n <- max(25, min(1000, as.integer(n)))
    if (any(!is.finite(c(lim$x, lim$y)))) return(NULL)
    
    if (is.finite(bw) && bw > 0) {
      # ABSOLUTE bandwidth (same for x and y)
      h <- c(bw, bw)
    } else {
      # Scott rule (isotropic), then scale by bw_factor
      vx  <- stats::var(df$x, na.rm = TRUE)
      vy  <- stats::var(df$y, na.rm = TRUE)
      sdx <- sqrt(max(vx, 0)); sdy <- sqrt(max(vy, 0))
      sx  <- min(sdx, stats::IQR(df$x, na.rm = TRUE) / 1.34)
      sy  <- min(sdy, stats::IQR(df$y, na.rm = TRUE) / 1.34)
      s   <- sqrt((sx^2 + sy^2) / 2)
      npts <- nrow(df)
      hiso <- 1.06 * s * (npts ^ (-1/6))
      if (!is.finite(hiso) || hiso <= 0) {
        rngx <- diff(range(df$x, na.rm = TRUE))
        rngy <- diff(range(df$y, na.rm = TRUE))
        hiso <- max((rngx + rngy)/2 * 0.02, 1e-6)
      }
      hiso <- hiso * if (is.finite(bw_factor) && bw_factor > 0) bw_factor else 1
      h <- c(hiso, hiso)
    }
    
    MASS::kde2d(
      df$x, df$y, n = n,
      lims = c(lim$x[1], lim$x[2], lim$y[1], lim$y[2]),
      h = h
    )
  }
  
  # ---- Reset filters to full ranges ----
  observeEvent(uploaded_df(), {
    # When filters are locked, sd_ranges_ui already initialises the widgets from
    # sd_filt (via isolate(sd_filt) in the renderUI).  Calling updateNumericRangeInput
    # here would override those correct locked values with the new file's max range.
    if (isTRUE(isolate(input$lock_filters))) return()

    df0 <- uploaded_df()
    req(df0)

    sda <- sd_by_tid(df0, min_n = input$sd_min_locs)
    if (!is.data.frame(sda) || !nrow(sda)) return()

    max_or1 <- function(v) {
      v <- v[is.finite(v)]
      if (!length(v)) return(1)
      mx <- max(v, na.rm = TRUE)
      if (!is.finite(mx) || mx <= 0) 1 else mx
    }

    shinyWidgets::updateNumericRangeInput(session, "sd_x_range", value = c(0, max_or1(sda$sd_x)))
    shinyWidgets::updateNumericRangeInput(session, "sd_y_range", value = c(0, max_or1(sda$sd_y)))

    if ("sd_z" %in% names(sda) && any(is.finite(sda$sd_z))) {
      shinyWidgets::updateNumericRangeInput(session, "sd_z_range", value = c(0, max_or1(sda$sd_z)))
    } else {
      # if z doesn't exist, set a wide range so it never filters anything
      shinyWidgets::updateNumericRangeInput(session, "sd_z_range", value = c(0, Inf))
    }
  }, ignoreInit = TRUE)
  
  
  raw_df <- reactive({
    df <- uploaded_df()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data uploaded"))
    df
  })
  
  observeEvent(input$hist_variable, {
    df <- filtered_df()
    if (input$hist_variable == "tid") {
      tlens <- as.numeric(table(df$tid))
      if (!length(tlens)) return()
      rng <- range(tlens, na.rm = TRUE)
      if (all(is.finite(rng))) {
        shinyWidgets::updateNumericRangeInput(session, "hist_range", value = rng)
      }
    } else {
      vals <- df[[input$hist_variable]]
      if (is.null(vals) || !is.numeric(vals) || !length(vals)) return()
      rng <- range(vals, na.rm = TRUE)
      shinyWidgets::updateNumericRangeInput(session, "hist_range", value = rng)
    }
  })
  # Use trace medians (or de-blinked molecule medians when de-blinking is on),
  # unless the user explicitly requests raw points.
  heatmap_input_df <- function(df, use_raw = FALSE, clustering = FALSE) {
    if (isTRUE(use_raw)) {
      return(df %>% dplyr::select(x, y))
    }
    if (isTRUE(clustering) && "cluster" %in% names(df) && any(df$cluster > 0, na.rm = TRUE)) {
      return(
        df %>%
          dplyr::filter(is.finite(cluster) & cluster > 0) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            .groups = "drop"
          )
      )
    }
    if ("tid" %in% names(df)) {
      return(
        df %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            .groups = "drop"
          )
      )
    }
    df %>% dplyr::select(x, y)
  }
  
  
  # Reads, optionally rotates, and prepares a data-URI + pixel size
  confocal_meta <- reactive({
    file <- input$bg_image
    if (is.null(file) || !isTRUE(input$show_bg)) return(NULL)
    
    img <- tryCatch(magick::image_read(file$datapath), error = function(e) NULL)
    if (is.null(img)) return(NULL)
    
    rot <- suppressWarnings(as.integer(if (is.null(input$img_rotation)) 0 else input$img_rotation))
    if (is.finite(rot) && (rot %% 360) != 0) {
      img <- magick::image_rotate(img, rot)
    }
    # 2) Optional: mirror horizontally (left↔right)
    if (isTRUE(input$img_flip_x)) {
      img <- magick::image_flop(img)   # horizontal mirror
    }
    # 3) Optional: recolour through a single-colour LUT (Green, Magenta, ...)
    #    plus brightness / contrast adjustment.
    img <- apply_background_lut(img, input$img_lut %||% "original",
                                custom_color = input$img_lut_custom,
                                brightness   = input$img_brightness %||% 100,
                                contrast     = input$img_contrast   %||% 0)

    info <- magick::image_info(img)
    raw_png <- magick::image_write(img, format = "png")
    list(
      uri = paste0("data:image/png;base64,", base64enc::base64encode(raw_png)),
      width_px  = info$width,
      height_px = info$height
    )
  })
  
  image_overlay <- reactive({
    meta <- confocal_meta()
    if (is.null(meta)) return(NULL)
    
    ps <- suppressWarnings(as.numeric(input$px_size));  if (!is.finite(ps) || ps <= 0) return(NULL)
    x0 <- suppressWarnings(as.numeric(input$img_x0));   if (!is.finite(x0)) x0 <- 0
    y0 <- suppressWarnings(as.numeric(input$img_y0));   if (!is.finite(y0)) y0 <- 0
    
    sizex <- meta$width_px  * ps
    sizey <- meta$height_px * ps
    
    if (isTRUE(input$img_invert_y)) {
      y_anchor <- "top";    y_val <- y0 + sizey
    } else {
      y_anchor <- "bottom"; y_val <- y0
    }
    
    list(
      source  = meta$uri,
      xref    = "x",   yref = "y",
      x       = x0,    y    = y_val,
      sizex   = sizex, sizey = sizey,
      xanchor = "left", yanchor = y_anchor,
      sizing  = "stretch",
      opacity = if (is.null(input$img_opacity)) 0.6 else input$img_opacity,
      layer   = "below"
    )
  })
  
  # helper: expand axes to include the image bbox; safe if overlay is NULL
  expand_limits_with_image <- function(lim, overlay, invert_y = TRUE) {
    if (is.null(overlay)) return(lim)
    x0 <- overlay$x;  w <- overlay$sizex
    y  <- overlay$y;  h <- overlay$sizey
    xr <- c(x0, x0 + w)
    yr <- if (invert_y) c(y - h, y) else c(y, y + h)
    list(x = range(c(lim$x, xr), na.rm = TRUE),
         y = range(c(lim$y, yr), na.rm = TRUE))
  }

  # ---- Auto-align overlay helpers ----
  .load_smlm_img <- function() {
    f <- input$bg_image
    if (is.null(f)) return(NULL)
    img <- tryCatch(magick::image_read(f$datapath), error = function(e) NULL)
    if (is.null(img)) return(NULL)
    rot <- suppressWarnings(as.integer(
      if (is.null(input$img_rotation)) 0 else input$img_rotation))
    if (is.finite(rot) && (rot %% 360) != 0) img <- magick::image_rotate(img, rot)
    if (isTRUE(input$img_flip_x)) img <- magick::image_flop(img)
    img
  }

  observeEvent(input$auto_foci, {
    img <- .load_smlm_img(); req(!is.null(img))
    df  <- isolate(plottable_data()); req(nrow(df) > 0)
    ps  <- suppressWarnings(as.numeric(input$px_size))
    req(is.finite(ps), ps > 0)
    res <- tryCatch(
      auto_align_overlay(img, df, ps, mode = "foci",
                         invert_y    = isTRUE(input$img_invert_y),
                         n_peaks        = as.integer(input$foci_n_peaks        %||% 5L),
                         bright_pct     = (input$foci_bright_pct     %||% 95)  / 100,
                         border_pct     = (input$foci_border_pct     %||% 10)  / 100,
                         eps_pct        = (input$foci_eps_pct        %||% 3)   / 100,
                         aspect_max     =  input$foci_aspect_max     %||% 4.0,
                         density_thr    = (input$foci_density_thr    %||% 20)  / 100,
                         density_cell   =  input$foci_density_cell   %||% 100,
                         density_smooth =  input$foci_density_smooth %||% 100),
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
    if (!is.null(res$img_peaks)) {
      smlm_diag(list(
        data_clusters = res$data_peaks,
        img_blobs     = data.frame(x = res$x0 + res$img_peaks$x,
                                   y = res$y0 + res$img_peaks$y)
      ))
    }
    updateCheckboxInput(session, "show_bg", value = TRUE)
    updateNumericInput(session, "img_x0", value = res$x0)
    updateNumericInput(session, "img_y0", value = res$y0)
    smlm_diag_mats(list(raw = res$img_mat_debug, smooth = res$img_mat_smooth))
    if (isTRUE(input$show_diag)) {
      pv <- res$img_peaks_debug
      mr <- res$mat_range_debug
      dd <- res$img_dim_debug
      dim_str <- if (!is.null(dd))
        paste0(" | src ", dd[1], "×", dd[2], " → proc ", dd[3], "×", dd[4],
               " px, raw bytes: ", dd[5], " (", dd[6], ")")
      else ""
      showNotification(
        paste0("Confocal detection - gray range: ", mr[1], "–", mr[2],
               " | peaks: ", paste(round(as.numeric(pv), 0), collapse = ", "),
               dim_str),
        type = "message", duration = 20)
    }
    showNotification("Foci alignment applied - nudge to fine-tune.",
                     type = "message", duration = 4)
  }, ignoreInit = TRUE)

  observeEvent(input$auto_boundary, {
    img <- .load_smlm_img(); req(!is.null(img))
    df  <- isolate(plottable_data()); req(nrow(df) > 0)
    ps  <- suppressWarnings(as.numeric(input$px_size))
    req(is.finite(ps), ps > 0)
    res <- tryCatch(
      auto_align_overlay(img, df, ps, mode = "boundary",
                         invert_y = isTRUE(input$img_invert_y)),
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
    updateCheckboxInput(session, "show_bg", value = TRUE)
    updateNumericInput(session, "img_x0", value = res$x0)
    updateNumericInput(session, "img_y0", value = res$y0)
    showNotification("Boundary alignment applied - nudge to fine-tune.",
                     type = "message", duration = 4)
  }, ignoreInit = TRUE)

  # ---- Nudge buttons for X and Y origin ----
  local({
    deltas <- c(m1000 = -1000, m100 = -100, m10 = -10,
                p10   =    10, p100 =  100, p1000 = 1000)
    for (sfx in names(deltas)) {
      local({
        s <- sfx
        d <- deltas[[s]]
        observeEvent(input[[paste0("nx_", s)]], {
          v <- suppressWarnings(as.numeric(input$img_x0))
          updateNumericInput(session, "img_x0",
                             value = (if (is.finite(v)) v else 0) + d)
        }, ignoreInit = TRUE)
        observeEvent(input[[paste0("ny_", s)]], {
          v <- suppressWarnings(as.numeric(input$img_y0))
          updateNumericInput(session, "img_y0",
                             value = (if (is.finite(v)) v else 0) + d)
        }, ignoreInit = TRUE)
      })
    }
  })

  output$download_overlay_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
      paste0(base, "_image_alignment_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content  = function(file) {
      ov <- list(
        show_confocal = input$show_confocal,
        show_bg       = input$show_bg,
        bg_image_name = if (!is.null(input$bg_image)) input$bg_image$name else NULL,
        px_size       = input$px_size,
        img_x0        = input$img_x0,
        img_y0        = input$img_y0,
        img_invert_y  = input$img_invert_y,
        img_flip_x    = input$img_flip_x,
        img_opacity   = input$img_opacity,
        img_rotation  = input$img_rotation,
        img_lut       = input$img_lut,
        img_lut_custom = input$img_lut_custom,
        img_brightness = input$img_brightness,
        img_contrast   = input$img_contrast
      )
      writeLines(jsonlite::toJSON(ov, pretty = TRUE, auto_unbox = TRUE, na = "null"),
                 file, useBytes = TRUE)
    }
  )
  # --- Foci detection debug plot (shown when Show diagnostics is checked) ----
  output$foci_debug_plot <- renderPlot({
    req(isTRUE(input$show_diag))
    mats <- smlm_diag_mats(); req(!is.null(mats))
    raw    <- mats$raw;    req(!is.null(raw))
    smooth <- mats$smooth; req(!is.null(smooth))
    par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))
    image(t(raw)[, nrow(raw):1],    col = grDevices::gray.colors(256),
          main = "Raw gray (border masked)", axes = FALSE)
    image(t(smooth)[, nrow(smooth):1], col = grDevices::gray.colors(256),
          main = "After smoothing (detection input)", axes = FALSE)
  }, bg = "#222")

  # --- Main 2D plot -----------------------------------------------------------

  # Shared "exploratory + references" footer for the GT / shape read-outs.
  gt_exploratory_note <- function() {
    div(style = "margin-top:8px;padding:6px 8px;border-left:3px solid #e67e22;background:rgba(230,126,34,0.06);font-size:11.5px;color:#7f8c8d;line-height:1.4;",
        HTML(paste0(
          "<b style='color:#e67e22;'>Exploratory.</b> This is a model-based estimate, not a direct measurement; validate against your own controls. ",
          "Methods: gyration-tensor shape anisotropy (Rudnick &amp; Gaspari, <i>J. Phys. A</i> <b>19</b>, L191, 1986); ",
          "linear/planar/spherical anisotropy coefficients (Westin et&nbsp;al., <i>MICCAI/Med. Image Anal.</i>, 1997/2002); ",
          "minimum-spanning-tree backbone (Prim, 1957). ",
          "The labelling inverse is a heuristic of this app, checked against the simulator rather than a known structure. Cite it as exploratory.")))
  }

  # ── GT counting-standard calibration ───────────────────────────────────────
  # Populate the standard-dataset selector from loaded files.
  observeEvent(smlm_file_store(), {
    fs  <- names(smlm_file_store())
    sel <- isolate(input$gt_calib_file)
    updateSelectInput(session, "gt_calib_file", choices = fs,
      selected = if (!is.null(sel) && sel %in% fs) sel else (isolate(input$active_file) %||% fs[1]))
  }, ignoreNULL = FALSE)

  # Measure observed TIDs-per-molecule from the chosen standard: apply the SAME
  # metadata filters + clustering settings the sample uses, then measure per
  # Median Cluster (each = one standard structure of known stoichiometry).
  gt_calib <- reactiveVal(NULL)
  observeEvent(input$gt_calib_run, {
    store <- smlm_file_store(); fn <- input$gt_calib_file
    if (is.null(fn) || !nzchar(fn %||% "") || !fn %in% names(store)) {
      showNotification("Pick a loaded standard dataset first.", type = "warning"); return()
    }
    k <- max(1L, as.integer(input$gt_calib_k %||% 1L))
    res <- tryCatch({
      df <- apply_clustering_pipeline(apply_meta_filters_batch(store[[fn]], input), input)
      measure_label_multiplicity(df, k_per_structure = k)
    }, error = function(e) { message("[GT-calib] ", conditionMessage(e)); NULL })
    if (is.null(res) || !isTRUE(res$ok)) {
      gt_calib(NULL)
      showNotification("No Median Clusters in the standard. Enable de-blinking and Median Clustering, then re-run.",
                       type = "warning", duration = 6); return()
    }
    res$file <- fn; gt_calib(res)
    showNotification(sprintf("Standard '%s': %.1f TIDs/molecule over %d structure(s).",
                             fn, res$m_tid, res$n_struct), type = "message", duration = 5)
  })

  output$gt_calib_readout <- renderUI({
    r <- gt_calib()
    if (is.null(r) || !isTRUE(r$ok))
      return(helpText(style = "color:#888;", "Measure a standard to calibrate the count."))
    active <- isTRUE(input$gt_use_calib)
    div(style = sprintf("font-size:12px; margin:4px 0; padding:5px 8px; border-left:3px solid %s; background:%s;",
                        if (active) "#1a8a1a" else "#b0b8c0", if (active) "rgba(26,138,26,0.06)" else "#f6f7f8"),
      HTML(sprintf("Measured <b>%.1f</b> TIDs / molecule (MAD %.1f, IQR %.0f&ndash;%.0f) from <b>%d</b> structure%s in <b>%s</b>%s. %s",
        r$m_tid, r$m_tid_mad, r$m_tid_iqr[1], r$m_tid_iqr[2], r$n_struct,
        if (r$n_struct == 1) "" else "s", r$file,
        if (r$k > 1) sprintf(" (÷ %d molecules each)", r$k) else "",
        if (active) "<b style=\"color:#1a8a1a\">In use for the count.</b>" else "Tick the box below to use it.")))
  })

  # ── Inverse GT-position estimate (labeling-aware) ───────────────────────────
  # Shared by the 2D/3D overlays and the readout panel. Returns NULL when the
  # feature is off or there is nothing to estimate.
  gt_estimate <- reactive({
    if (!isTRUE(input$show_gt_estimate)) return(NULL)
    df <- plottable_data()
    if (is.null(df) || !nrow(df)) return(NULL)

    presets     <- minflux_label_presets()
    preset_name <- input$gt_label_preset %||% "Primary monoclonal IgG + secondary IgG"
    preset      <- presets[[preset_name]]
    if (is.null(preset)) return(NULL)

    # Custom: build the preset from the numeric inputs.
    if (identical(preset_name, "Custom")) {
      preset <- list(
        n_fluor       = max(1, suppressWarnings(as.numeric(input$gt_fluor_custom   %||% 1))),
        link_nm       = max(0, suppressWarnings(as.numeric(input$gt_link_custom     %||% 5))),
        dye_spread_nm = max(0, suppressWarnings(as.numeric(input$gt_spread_custom   %||% 0))),
        footprint_nm  = max(0, suppressWarnings(as.numeric(input$gt_footprint_custom %||% 8))),
        blink_lambda  = max(0, suppressWarnings(as.numeric(input$gt_blink_custom    %||% 0.3))),
        note          = "Custom labeling parameters"
      )
    }

    eff <- suppressWarnings(as.numeric(input$gt_label_eff)); if (!is.finite(eff)) eff <- 70
    lvl <- input$gt_model_level %||% "auto"
    # Merge radius: blank / NA / negative -> auto (NULL); 0 -> no merge; >0 -> override.
    cr_in <- suppressWarnings(as.numeric(input$gt_complex_r))
    cr    <- if (length(cr_in) && is.finite(cr_in) && cr_in >= 0) cr_in else NULL
    # Inward displacement magnitude: blank -> auto (link_nm).
    dn_in <- suppressWarnings(as.numeric(input$gt_displace_nm))
    dn    <- if (length(dn_in) && is.finite(dn_in) && dn_in >= 0) dn_in else NULL
    tryCatch(
      estimate_gt_from_labeling(
        df, preset, eff_pct = eff, level = lvl, complex_r = cr,
        scope_cluster2  = isTRUE(input$gt_scope_cluster2 %||% TRUE),
        displace_inward = isTRUE(input$gt_displace_inward %||% TRUE),
        displace_nm     = dn,
        guided_merge    = isTRUE(input$gt_guided_merge),
        support_strength = { v <- suppressWarnings(as.numeric(input$gt_support_strength)); if (is.finite(v)) v else 1 },
        m_tid_override  = if (isTRUE(input$gt_use_calib)) {
          g <- gt_calib(); if (!is.null(g) && isTRUE(g$ok)) g$m_tid else NULL
        } else NULL
      ),
      error = function(e) { message("[GT-est] ", conditionMessage(e)); NULL }
    )
  })

  # Add the hollow-circle GT candidate overlay to a 2D figure (view-plane aware).
  add_gt_overlay_2d <- function(fig, vp) {
    ge <- gt_estimate()
    if (is.null(ge) || !nrow(ge$points)) return(fig)
    pts <- ge$points
    remap_x <- function(x, y, z) switch(vp, XZ = x, YZ = y, x)
    remap_y <- function(x, y, z) switch(vp, XZ = z, YZ = z, y)
    px <- remap_x(pts$x, pts$y, pts$z); py <- remap_y(pts$x, pts$y, pts$z)
    lo  <- ge$meta$link_nm %||% 0
    hi  <- max(ge$meta$footprint_nm %||% 0, 2 * lo, stats::quantile(pts$pos_sd_nm, 0.9, na.rm = TRUE))
    col <- gt_uncertainty_colors(pts$pos_sd_nm, lo, hi)
    sz  <- suppressWarnings(as.numeric(input$gt_circle_size)); if (!is.finite(sz) || sz <= 0) sz <- 15
    struct_lbl <- ifelse(pts$cluster2 > 0, paste0("Structure ", pts$cluster2), "ungrouped")
    txt <- sprintf("<b>%s</b><br>GT candidate · ± %.0f nm%s<br>%d contributing unit(s) · %d TID / %d loc",
                   struct_lbl, pts$pos_sd_nm, ifelse(pts$fused, " (fused: multiple GTs)", ""),
                   pts$n_units, pts$n_tid, pts$n_loc)

    # Contribution fan-out lines: each member unit centroid -> its GT centroid.
    if (isTRUE(input$gt_show_lines %||% TRUE) && !is.null(ge$links) && nrow(ge$links)) {
      lk <- ge$links
      lx <- as.vector(rbind(remap_x(lk$ux, lk$uy, lk$uz),
                            remap_x(lk$gx, lk$gy, lk$gz), NA_real_))
      ly <- as.vector(rbind(remap_y(lk$ux, lk$uy, lk$uz),
                            remap_y(lk$gx, lk$gy, lk$gz), NA_real_))
      line_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.35)" else "rgba(0,0,0,0.35)"
      fig <- fig %>% plotly::add_trace(
        x = lx, y = ly, type = "scatter", mode = "lines",
        line = list(color = line_col, width = 1), name = "GT contributions",
        showlegend = TRUE, inherit = FALSE, hoverinfo = "skip"
      )
    }
    fig %>% plotly::add_trace(
      x = px, y = py, type = "scatter", mode = "markers",
      marker = list(symbol = "circle-open", size = sz,
                    color = col, line = list(color = col, width = 3)),
      name = "Estimated GT", showlegend = TRUE, inherit = FALSE,
      text = txt, hoverinfo = "text"
    )
  }

  # Text/HTML readout of the GT estimate: counts, uncertainty floor, crowding.
  output$gt_estimate_readout <- renderUI({
    if (!isTRUE(input$show_gt_estimate))
      return(helpText("Enable the estimate to see candidate GT positions and counts."))
    ge <- gt_estimate()
    if (is.null(ge) || !nrow(ge$summary))
      return(helpText("No estimable structures (need clustered data; try enabling de-blinking and Median clustering)."))
    s   <- ge$summary
    m   <- ge$meta
    tot_units  <- sum(s$n_units)
    tot_gt     <- sum(s$n_gt_est)
    tot_stoich <- sum(s$n_stoich)
    tot_lo     <- sum(s$n_gt_lo);  tot_hi <- sum(s$n_gt_hi)
    any_crowd  <- any(s$crowded)
    any_mm     <- any(s$stoich_mismatch)
    lvl_lbl    <- if (identical(m$level, "primary")) "de-blinked molecules" else "raw trace IDs"

    rows <- lapply(seq_len(nrow(s)), function(i) {
      r <- s[i, ]
      lab <- if (r$cluster2 == 0) "All (ungrouped)" else paste("Structure", r$cluster2)
      warn <- if (isTRUE(r$crowded))
        span(style = "color:#e67e22;font-weight:600;", " ⚠ crowded / under-resolved") else ""
      mm <- if (isTRUE(r$stoich_mismatch))
        span(style = "color:#c0392b;font-weight:600;",
             sprintf(" ⚠ spatial vs stoichiometry disagree (%.0f obs TIDs/GT vs %.0f expected)",
                     r$tids_per_gt, m$m_tid)) else ""
      tags$li(HTML(sprintf(
        "<b>%s</b>: <b>%.0f</b> GT (spatial) vs <b>%.0f</b> (stoichiometry) · CI %.0f–%.0f · %d TID / %d unit%s%s%s%s",
        lab, r$n_gt_est, r$n_stoich, r$n_gt_lo, r$n_gt_hi, r$n_tid_total, r$n_units,
        ifelse(r$n_units == 1, "", "s"),
        ifelse(r$n_fused > 0, sprintf(", %d fused", r$n_fused), ""),
        as.character(warn), as.character(mm))))
    })

    tagList(
      div(style = "font-size:15px;margin-bottom:6px;",
          HTML(sprintf("Estimated total: <b>%.0f</b> GT molecules (spatial) · <b>%.0f</b> (stoichiometry) · range %.0f–%.0f",
                  tot_gt, tot_stoich, tot_lo, tot_hi))),
      tags$ul(style = "margin:4px 0 8px 0;padding-left:18px;", rows),
      if (any_mm) div(style = "font-size:12px;color:#c0392b;margin-bottom:6px;",
          HTML("The two independent counts disagree in ≥1 structure - the merge radius or the labeling assumption is off. Try the stoichiometry-guided merge, adjust the merge radius, or re-check the preset/efficiency.")),
      div(style = "font-size:12px;color:#7f8c8d;line-height:1.45;",
          HTML(sprintf(
            "Modeled from %s. Units pooled within <b>%.0f nm</b> (label-complex radius)%s into one molecule%s. Position floor <b>±%.0f nm</b> (linker, irreducible) · <b>%.1f TIDs/molecule (%s)</b> · detection η=%.0f%% · σ<sub>loc</sub>≈%.1f nm.<br>Marker colour = position uncertainty, combining the linker floor, localization sampling (σ<sub>loc</sub>/√n<sub>loc</sub>), a <b>support penalty</b> (fewer TIDs → less certain), inward-displacement direction error, and a fusion penalty.<br><i>Counts are conditional on the labeling calibration; positions cannot be resolved below the linker floor. The two counts bracket the honest range.%s</i>",
            lvl_lbl, m$complex_r_nm %||% 0,
            if (isTRUE(m$guided_merge)) " (stoichiometry-guided)" else "",
            if ((m$displace_nm %||% 0) > 0) sprintf(", then displaced <b>%.0f nm</b> inward toward the structure", m$displace_nm) else "",
            m$link_nm, m$tids_per_mol %||% m$m_tid,
            if (identical(m$m_tid_source, "measured")) "measured from a standard" else "from preset",
            100 * m$eta, m$sigma_loc_nm,
            ifelse(any_crowd, " Crowded structures undercount - treat those counts as lower bounds.", ""))),
      gt_exploratory_note())
    )
  })

  # ── Phase B: per-structure shape / polymer morphology ───────────────────────
  gt_morphology <- reactive({
    if (!isTRUE(input$gt_show_shape)) return(NULL)
    df <- plottable_data()
    if (is.null(df) || !nrow(df)) return(NULL)
    sp <- suppressWarnings(as.numeric(input$gt_subunit_spacing))
    if (!is.finite(sp) || sp <= 0) sp <- NA_real_
    # use the estimated GT candidates as a clean backbone when available
    bpts <- tryCatch({ ge <- gt_estimate(); if (!is.null(ge)) ge$points else NULL },
                     error = function(e) NULL)
    tryCatch(
      analyze_cluster2_morphology(df, backbone_pts = bpts, subunit_spacing_nm = sp),
      error = function(e) { message("[GT-morph] ", conditionMessage(e)); NULL }
    )
  })

  output$gt_morphology_readout <- renderUI({
    if (!isTRUE(input$gt_show_shape))
      return(helpText("Enable to classify each cluster's shape (filament / ring / punctum / sheet / diffuse) and report neutral size descriptors."))
    mo <- gt_morphology()
    if (is.null(mo) || !nrow(mo))
      return(helpText("No analysable structures (need Median Clusters with ≥4 traces)."))
    icon_for <- function(s) switch(s, filament = "――", ring = "○",
                                   punctum = "●", sheet = "▧", "∴")
    rows <- lapply(seq_len(nrow(mo)), function(i) {
      r <- mo[i, ]
      spc <- if (is.finite(r$spacing_nn_nm)) sprintf(", spacing %.0f nm", r$spacing_nn_nm) else ""
      detail <- if (r$shape == "filament")
        sprintf("contour <b>%.0f nm</b> (straight %.0f, width %.0f)%s, %.1f points/µm%s",
                r$contour_nm, r$extent_nm, r$width_nm, spc, r$lin_density_um,
                if (is.finite(r$n_subunit_est)) sprintf(" → <b>~%.0f subunits</b>", r$n_subunit_est) else "")
      else if (r$shape == "ring")
        sprintf("radius <b>%.0f nm</b> (circumference %.0f nm), R<sub>g</sub> %.0f%s", r$ring_r_nm, r$contour_nm, r$rg_nm, spc)
      else
        sprintf("extent <b>%.0f nm</b>, width %.0f, R<sub>g</sub> %.0f%s", r$extent_nm, r$width_nm, r$rg_nm, spc)
      tags$li(HTML(sprintf("<b>Structure %d</b> - <b>%s</b> %s: %s",
                           r$cluster2, r$shape, icon_for(r$shape), detail)))
    })
    tagList(
      tags$ul(style = "margin:4px 0 8px 0;padding-left:18px;", rows),
      div(style = "font-size:12px;color:#7f8c8d;line-height:1.45;",
          HTML("Shape is <b>inferred</b> from the gyration-tensor eigenvalues and a hollow-centre test for rings. Contour length follows the minimum spanning tree, so it tracks curves; spacing is the median nearest-neighbour distance."),
          gt_exploratory_note())
    )
  })

  output$scatter2D <- renderPlotly({
    req(input$mode3d == "2D")
    if (is.null(input$file) || is.null(plot_trigger())) {
      return(plotly::plot_ly(source = "smlm2d", type = "scatter", mode = "markers",
                             x = numeric(0), y = numeric(0)) %>%
               plotly::event_register("plotly_relayout"))
    }

    df <- plottable_data()
    shiny::validate(shiny::need(nrow(df) > 0, "No points to display - please check filters or upload a valid file."))

    # Remap x/y columns for XZ or YZ viewing planes
    vp <- input$view_plane_2d %||% "XY"
    axis_labels <- list(x = "X (nm)", y = "Y (nm)")
    if (vp == "XZ" && "z" %in% names(df)) {
      df <- df %>% dplyr::mutate(y = z)
      axis_labels <- list(x = "X (nm)", y = "Z (nm)")
    } else if (vp == "YZ" && "z" %in% names(df)) {
      df <- df %>% dplyr::mutate(x = y, y = z)
      axis_labels <- list(x = "Y (nm)", y = "Z (nm)")
    }

    df_raw <- uploaded_df()
    # Ensure x/y exist on raw
    df_raw <- df_raw %>% dplyr::mutate(x = .data[["x"]], y = .data[["y"]])

    # Base limits from the data
    lim <- get_equal_limits_2d(df$x, df$y)
    
    # Try to build an overlay (may be NULL)
    ovl <- image_overlay()
    
    # Expand limits to include image bbox
    lim <- expand_limits_with_image(lim, ovl, invert_y = isTRUE(input$img_invert_y))
    
    # -------------------------
    # 0) Start empty figure
    # -------------------------
    fig    <- plotly::plot_ly(source = "smlm2d")
    fg_col <- if (isTRUE(input$dark_mode)) "#FFFFFF" else "#63666A"

    # -------------------------
    # 1) HEATMAP (below points)
    # -------------------------
    if (isTRUE(input$show_heatmap)) {
      df_used <- heatmap_input_df(
        df,
        use_raw    = isTRUE(input$heatmap_use_raw),
        clustering = isTRUE(input$enable_clustering)
      )
      lim_hm <- get_equal_limits_2d(df_used$x, df_used$y)
      ovl_for_hm <- image_overlay()
      lim_hm <- expand_limits_with_image(lim_hm, ovl_for_hm, invert_y = isTRUE(input$img_invert_y))
      
      est <- compute_kde(
        df_used,
        lim = lim_hm,
        n   = if (is.null(input$heatmap_n)) 200 else input$heatmap_n,
        bw  = if (is.null(input$heatmap_bw_abs)) NA_real_ else input$heatmap_bw_abs,
        bw_factor = if (is.null(input$heatmap_bw_factor)) 1 else input$heatmap_bw_factor
      )
      
      if (!is.null(est)) {
        dens <- est$z
        floor_p <- if (is.null(input$heatmap_floor)) 0 else input$heatmap_floor
        if (is.finite(floor_p) && floor_p > 0) {
          thr <- stats::quantile(dens, probs = floor_p, na.rm = TRUE)
          dens[dens < thr] <- 0
          if (max(dens, na.rm = TRUE) > 0) dens <- dens / max(dens, na.rm = TRUE)
        }
        dens[!is.finite(dens)] <- 0
        
        clip_p <- if (is.null(input$heatmap_clip)) 0.99 else input$heatmap_clip
        cap    <- stats::quantile(dens, probs = clip_p, na.rm = TRUE)
        dens   <- pmin(dens, cap)
        
        if (max(dens, na.rm = TRUE) > 0) dens <- dens / max(dens, na.rm = TRUE)
        
        if (isTRUE(input$heatmap_log)) {
          dens <- log1p(dens)
          dens <- dens / max(dens, na.rm = TRUE)
        }
        
        gm <- if (is.null(input$heatmap_gamma)) 1 else input$heatmap_gamma
        if (is.finite(gm) && gm > 0) {
          dens <- dens ^ gm
          if (max(dens, na.rm = TRUE) > 0) dens <- dens / max(dens, na.rm = TRUE)
        }
        
        fig <- fig %>% plotly::add_heatmap(
          x = est$x,
          y = est$y,
          z = t(dens),
          zmin = 0, zmax = 1,
          colors = palette_colors(input$heatmap_palette),
          showscale = isTRUE(input$heatmap_showlegend),
          colorbar = list(
            titlefont = list(color = fg_col),
            tickfont  = list(color = fg_col)
          ),
          opacity = if (is.null(input$heatmap_opacity)) 0.35 else input$heatmap_opacity,
          hoverinfo = "skip",
          inherit = FALSE
        )
        
        # lock ranges to heatmap grid
        lim <- lim_hm
      }
    }
    
    # -------------------------
    # 2) BACKGROUND IMAGE (below points)
    # -------------------------
    if (!is.null(ovl)) {
      fig <- fig %>% layout(images = list(ovl))
    }
    
    # -----------------------------------
    # 3) POINTS (build pts + marker once)
    # -----------------------------------
    # --- build the points once (same idea as 3D) ---
    pts <- build_pts_2d(df, input)
    pts$color <- build_colors_for_pts(pts, input, default_palette_col = "color")
    cmap <- color_mapping_for_pts(pts, input)
    
    # bail out early if empty; still keep overlay + axes
    if (!nrow(pts)) {
      if (!is.null(ovl)) fig <- fig %>% layout(images = list(ovl))
      return(
        fig %>% layout(
          title = "2D Scatter Plot",
          xaxis = list(title = axis_labels$x, range = lim$x, scaleanchor = "y", scaleratio = 1),
          yaxis = list(title = axis_labels$y, range = lim$y),
          paper_bgcolor = if (isTRUE(input$dark_mode)) "#1e1e1e" else "white",
          plot_bgcolor  = if (isTRUE(input$dark_mode)) "#1e1e1e" else "white"
        )
      )
    }
    
    # --- draw the points once ---
    mplex_draw <- identical(input$color_mode %||% "default", "multiplex") &&
                  "mplex_label" %in% names(pts) && any(!is.na(pts$mplex_label))
    if (mplex_draw) {
      # one trace per channel so each gets a named, colored legend entry
      ptso  <- pts[order(pts$mplex_channel), , drop = FALSE]
      chans <- unique(ptso$mplex_channel[!is.na(ptso$mplex_channel)])
      for (chv in chans) {
        sub <- ptso[!is.na(ptso$mplex_channel) & ptso$mplex_channel == chv, , drop = FALSE]
        if (!nrow(sub)) next
        fig <- fig %>% plotly::add_trace(
          data = sub, x = ~x, y = ~y,
          type = "scatter", mode = "markers",
          name = sub$mplex_label[1],
          marker = list(size = input$median_point_size %||% 4, color = sub$mplex_color[1]),
          text = ~hover_text,
          hovertemplate = "%{text}<extra></extra>",
          hoverinfo = "text",
          showlegend = TRUE, inherit = FALSE
        )
      }
    } else if (isTRUE(cmap$showscale) && identical(cmap$kind, "scale")) {
      fig <- fig %>% plotly::add_trace(
        data = pts,
        x = ~x, y = ~y,
        type = "scatter", mode = "markers",
        name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
        marker = list(
          size = input$median_point_size %||% 4,
          color = cmap$colorvar,
          colorscale = cmap$colorscale,
          cmin = cmap$cmin, cmax = cmap$cmax,
          showscale = TRUE,
          colorbar = list(
            title     = cmap$title,
            titlefont = list(color = fg_col),
            tickfont  = list(color = fg_col)
          )
        ),
        text = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        hoverinfo = "text",
        showlegend = TRUE, inherit = FALSE
      )
    } else {
      hex <- if (!is.null(cmap$hex)) cmap$hex else rep("#bdbdbd", nrow(pts))
      fig <- fig %>% plotly::add_trace(
        data = pts,
        x = ~x, y = ~y,
        type = "scatter", mode = "markers",
        name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
        marker = list(size = input$median_point_size %||% 4, color = hex),
        text = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        hoverinfo = "text",
        showlegend = TRUE, inherit = FALSE
      )
    }
    
    # ------------------------------------
    # 5) CLUSTER2 NN / FITS (unchanged)
    # ------------------------------------
    if (isTRUE(input$show_cluster2_nn_lines) && isTRUE(input$enable_cluster2) && !is.null(df$cluster2)) {
      lines <- secondary_nn_lines()
      if (length(lines)) {
        seg_x <- c(); seg_y <- c(); lab_x <- c(); lab_y <- c(); labs <- c()
        for (ln in lines) {
          seg_x <- c(seg_x, ln$x, NA_real_)
          seg_y <- c(seg_y, ln$y, NA_real_)
          lab_x <- c(lab_x, ln$midpoint[1]); lab_y <- c(lab_y, ln$midpoint[2]); labs <- c(labs, ln$text)
        }
        fig <- fig %>% plotly::add_trace(
          x = seg_x, y = seg_y, type = "scatter", mode = "lines",
          line = list(color = "green", width = 2),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
        if (length(labs) && isTRUE(input$show_c2_nn_labels)) {
          nn_lfs <- suppressWarnings(as.numeric(input$c2_label_font_size %||% 12))
          if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
          fig <- fig %>% plotly::add_trace(
            x = lab_x, y = lab_y, type = "scatter", mode = "text",
            text = labs, textfont = list(color = "green", size = nn_lfs),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      }
    }


    # ------------------------------------
    # 6) GENERAL NN + 2-color NN (your blocks)
    # ------------------------------------
    if (isTRUE(input$show_nn_lines) && nrow(df) > 1) {
      nn <- build_nn_segments(df, input, is3d = FALSE)
      if (!is.null(nn)) {
        ln_col <- if (isTRUE(input$dark_mode)) "white" else "black"
        # With channels, colour each line by the channels it joins: a
        # within-channel line takes that channel's colour, a cross-channel line
        # is magenta. Without channels it is one plain series, as before.
        ed <- nn$edges
        if (!is.null(ed) && nrow(ed) && length(unique(stats::na.omit(ed$ch1))) +
            length(unique(stats::na.omit(ed$ch2))) > 2L) {
          for (cc in unique(ed$color)) {
            sel <- which(ed$color == cc)
            fig <- fig %>% plotly::add_trace(
              x = as.numeric(t(cbind(ed$x1[sel], ed$x2[sel], NA))),
              y = as.numeric(t(cbind(ed$y1[sel], ed$y2[sel], NA))),
              type = "scatter", mode = "lines",
              line = list(color = cc, width = 2),
              hoverinfo = "none", showlegend = FALSE, inherit = FALSE
            )
          }
        } else {
        fig <- fig %>% plotly::add_trace(
          x = nn$seg_x, y = nn$seg_y,
          type = "scatter", mode = "lines",
          line = list(color = ln_col, width = 1),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
        }
        if (isTRUE(input$show_nn_labels) && length(nn$labs)) {
          lab_col <- if (isTRUE(input$dark_mode)) "orange" else "red"
          nn_lfs <- suppressWarnings(as.numeric(input$nn_label_font_size %||% 12))
          if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
          fig <- fig %>% plotly::add_trace(
            x = nn$lab_x, y = nn$lab_y,
            type = "scatter", mode = "text",
            text = nn$labs, textfont = list(color = lab_col, size = nn_lfs),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      }
    }
    
    # The separate 2-colour NN block that used to live here is gone. NN lines are
    # drawn once, by the shared panel above, which now colours each line by the
    # channels it joins - so enabling channels ADAPTS the Nearest Neighbour
    # controls rather than adding a parallel set of them.

    # ------------------------------------
    # 7) AXES + THEMING (equal aspect)
    # ------------------------------------
    bg_col   <- if (isTRUE(input$dark_mode)) "#000000" else "#FFFFFF"
    fg_col   <- if (isTRUE(input$dark_mode)) "#FFFFFF" else "#63666A"
    show_g   <- !isTRUE(input$hide_grid)
    grid_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
    zero_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.25)" else "rgba(0,0,0,0.25)"
    
    fig <- fig %>%
      plotly::layout(
        title = "2D Scatter Plot",
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        xaxis = list(
          title = axis_labels$x,
          range = lim$x,
          showgrid = show_g, gridcolor = grid_col, gridwidth = 0.4,
          zeroline = show_g, zerolinecolor = zero_col, zerolinewidth = 0.6,
          showticklabels = show_g, tickfont = list(color = fg_col),
          color = fg_col
        ),
        yaxis = list(
          title = axis_labels$y,
          range = lim$y,
          showgrid = show_g, gridcolor = grid_col, gridwidth = 0.4,
          zeroline = show_g, zerolinecolor = zero_col, zerolinewidth = 0.6,
          showticklabels = show_g, tickfont = list(color = fg_col),
          color = fg_col,
          scaleanchor = "x", scaleratio = 1
        ),
        showlegend = TRUE,
        font = list(color = fg_col),
        legend = list(
          font = list(color = fg_col),
          orientation = "h",
          x = 0, y = 0.97,
          xanchor = "left", yanchor = "bottom"
        )
      )

    smlm_axis_x_span(diff(lim$x))

    # --- 2D Cluster2 convex hulls (optional) ---
    if (!is.null(input$show_cluster2_shape_2D) &&
        isTRUE(input$show_cluster2_shape_2D) &&
        isTRUE(input$enable_cluster2) &&
        !is.null(df$cluster2)) {
      
      hull_opacity <- if (!is.null(input$hull_opacity)) input$hull_opacity else 0.2
      
      for (cl2 in unique(df$cluster2[df$cluster2 != 0])) {
        group_df <- df %>% dplyr::filter(cluster2 == cl2)
        if (nrow(group_df) >= 3) {
          hidx <- tryCatch(grDevices::chull(group_df$x, group_df$y), error = function(e) integer(0))
          if (length(hidx) >= 3) {
            hull_xy <- group_df[c(hidx, hidx[1]), c("x","y")]
            base_col <- if ("cluster2_color" %in% names(group_df) && !is.na(group_df$cluster2_color[1]))
              group_df$cluster2_color[1] else "#00BFFF"
            fill_col <- grDevices::adjustcolor(base_col, alpha.f = hull_opacity)
            
            fig <- fig %>%
              plotly::add_trace(
                x = hull_xy$x, y = hull_xy$y,
                type = "scatter", mode = "lines",
                fill = "toself",
                fillcolor = fill_col,
                line = list(color = base_col, width = 2),
                name = paste("Median Cluster", cl2),
                showlegend = TRUE,
                inherit = FALSE
              )
          }
        }
      }
    }

    # ---- Scale bar ----
    if (isTRUE(input$show_scale_bar) && input$mode3d == "2D") {
      sb_len   <- suppressWarnings(as.numeric(input$scale_bar_length))
      sb_thick <- suppressWarnings(as.numeric(input$scale_bar_thickness))
      sb_font  <- suppressWarnings(as.numeric(input$scale_font_size))
      sb <- make_scale_bar(sb_len, diff(lim$x), sb_thick, dark = isTRUE(input$dark_mode),
                           scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
      if (!is.null(sb))
        fig <- fig %>% plotly::layout(shapes = sb$shapes,
                                      annotations = sb$annotations)
    }

    # ---- Diagnostic: detected cluster centres (red +) and confocal blobs (green ○) ----
    # type="scatter" (SVG) - scattergl lacks circle-open / cross symbol support.
    diag <- smlm_diag()   # non-isolated: toggling show_diag or new alignment re-renders
    if (isTRUE(input$show_diag) && !is.null(diag)) {
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

    # ---- Estimated ground-truth candidates (hollow circles, green->red) ----
    fig <- add_gt_overlay_2d(fig, vp)

    fig %>% plotly::event_register('plotly_relayout')
  })

  observe({
    # Only read the relayout event while the SMLM tab is on screen, so plotly's
    # deferred "source not registered" warning is never scheduled.
    if (!identical(input$main_navbar, "tab_smlm")) return()
    req(input$file, plot_trigger(), identical(input$mode3d, "2D"))
    ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "smlm2d"))
    if (!is.list(ed)) return()

    x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
    x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
    if (length(x0) == 1 && length(x1) == 1 && is.finite(x0) && is.finite(x1) && x1 > x0) {
      smlm_axis_x_span(x1 - x0)
      y0 <- suppressWarnings(as.numeric(ed[["yaxis.range[0]"]]))
      y1 <- suppressWarnings(as.numeric(ed[["yaxis.range[1]"]]))
      if (length(y0) == 1 && length(y1) == 1 && is.finite(y0) && is.finite(y1) && y1 > y0)
        smlm_zoom_ranges(list(x = c(x0, x1), y = c(y0, y1)))
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      span <- smlm_default_x_span()
      if (length(span) == 1 && is.finite(span) && span > 0) smlm_axis_x_span(span)
      smlm_zoom_ranges(NULL)
    } else {
      return()
    }

    if (isTRUE(input$show_scale_bar)) smlm_push_scale_bar()
  })

  observeEvent(list(input$show_scale_bar, input$scale_bar_length,
                    input$scale_bar_thickness, input$scale_font_size,
                    input$dark_mode, input$mode3d), {
    smlm_push_scale_bar()
  }, ignoreInit = TRUE)

  # ---- 2D crop-to-view → feeds the existing box trimmer -----------------------
  # When the user zooms in 2D and presses "Crop to view", we enable the box
  # trimmer and push the current axis ranges into trim_x_range / trim_y_range.
  # Because those inputs are created by renderUI (trim_box_ui) which only renders
  # after enable_trimmer=TRUE + trim_mode="box", we use a short later() delay to
  # let the DOM catch up before sending the update.  The crop then persists
  # automatically when the user switches to 3D.
  observeEvent(input$smlm_crop_to_view, {
    zr <- smlm_zoom_ranges()
    if (is.null(zr)) {
      showNotification("Zoom into a region first, then press Crop to view.",
                       type = "warning", duration = 5)
      return()
    }
    smlm_2d_crop_rv(zr)
    plot_trigger(Sys.time())
    showNotification(
      sprintf("Crop active: x [%.0f, %.0f] nm  |  y [%.0f, %.0f] nm (retained in 3D).",
              zr$x[1], zr$x[2], zr$y[1], zr$y[2]),
      type = "message", duration = 6
    )
  }, ignoreInit = TRUE)

  observeEvent(input$smlm_clear_crop, {
    smlm_2d_crop_rv(NULL)
    plot_trigger(Sys.time())
    showNotification("Spatial crop cleared.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  spatial_info_ui <- reactive({
    crop <- smlm_2d_crop_rv()
    if (is.null(crop))
      tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
             "No 2D crop active.")
    else
      tags$p(style = "font-size:0.82em; color:#2070a0; margin:2px 0;",
             sprintf("Active: x [%.0f, %.0f]  |  y [%.0f, %.0f] nm",
                     crop$x[1], crop$x[2], crop$y[1], crop$y[2]))
  })
  output$smlm_spatial_info    <- renderUI({ spatial_info_ui() })
  output$smlm_spatial_info_3d <- renderUI({ spatial_info_ui() })

  observeEvent(input$smlm_clear_crop_3d, {
    smlm_2d_crop_rv(NULL)
    plot_trigger(Sys.time())
    showNotification("Spatial crop cleared.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  #--- 3D-Only server parts:
  
  
  df_trim_ref <- reactive({
    df <- filtered_df()
    df
  })
  
  output$trim_box_ui <- renderUI({
    df <- df_trim_ref(); shiny::validate(shiny::need(nrow(df) > 0, "No data"))
    rx <- range(df$x, na.rm = TRUE); ry <- range(df$y, na.rm = TRUE); rz <- range(df$z, na.rm = TRUE)
    tagList(
      shinyWidgets::numericRangeInput("trim_x_range", "Keep X range", value = rx, min = rx[1], max = rx[2]),
      shinyWidgets::numericRangeInput("trim_y_range", "Keep Y range", value = ry, min = ry[1], max = ry[2]),
      shinyWidgets::numericRangeInput("trim_z_range", "Keep Z range", value = rz, min = rz[1], max = rz[2])
    )
  })
  
  output$trim_sphere_ui <- renderUI({
    df <- df_trim_ref(); shiny::validate(shiny::need(nrow(df) > 0, "No data"))
    c0 <- colMeans(df[,c("x","y","z")], na.rm = TRUE)
    r0 <- stats::quantile(sqrt(rowSums((as.matrix(df[,c("x","y","z")]) - 
                                          matrix(c0, nrow(df), 3, TRUE))^2)), probs = 0.95, na.rm = TRUE)
    tagList(
      numericInput("trim_cx", "Center X", value = round(c0[1], 2)),
      numericInput("trim_cy", "Center Y", value = round(c0[2], 2)),
      numericInput("trim_cz", "Center Z", value = round(c0[3], 2)),
      numericInput("trim_r",  "Radius",   value = round(r0, 2), min = 0)
    )
  })
  
  # When "Lock filters" is on, preserve the current widget value across file switches.
  # The widget is still rebuilt (so bounds reflect the new file's data range), but
  # the displayed value comes from the last-set filter rather than the new data's range.
  .locked_val <- function(cur_input, rng) {
    if (isTRUE(input$lock_filters)) isolate(cur_input) %||% rng else rng
  }
  output$efo_range_ui <- renderUI({
    df  <- uploaded_df()
    rng <- range(df$efo, na.rm = TRUE)
    shinyWidgets::numericRangeInput("efo_range", "EFO Range",
      value = .locked_val(input$efo_range, rng), min = rng[1], max = rng[2])
  })
  output$cfr_range_ui <- renderUI({
    df  <- uploaded_df()
    rng <- range(df$cfr, na.rm = TRUE)
    shinyWidgets::numericRangeInput("cfr_range", "CFR Range",
      value = .locked_val(input$cfr_range, rng), min = rng[1], max = rng[2])
  })
  output$tim_range_ui <- renderUI({
    df  <- uploaded_df()
    rng <- range(df$tim, na.rm = TRUE)
    shinyWidgets::numericRangeInput("tim_range", "TIM Range",
      value = .locked_val(input$tim_range, rng), min = rng[1], max = rng[2])
  })
  output$dcr_range_ui <- renderUI({
    df  <- uploaded_df()
    rng <- range(df$dcr, na.rm = TRUE)
    shinyWidgets::numericRangeInput("dcr_range", "DCR Range",
      value = .locked_val(input$dcr_range, rng), min = rng[1], max = rng[2])
  })
  output$dwell_range_ui <- renderUI({
    df  <- uploaded_df()
    if (!"dwell" %in% names(df)) return(NULL)
    rng <- range(df$dwell, na.rm = TRUE)
    shinyWidgets::numericRangeInput("dwell_range", "DWELL Range",
      value = .locked_val(input$dwell_range, rng), min = rng[1], max = rng[2])
  })
  output$fluo_range_ui <- renderUI({
    df  <- uploaded_df()
    if (!"fluo" %in% names(df)) return(NULL)
    rng <- range(df$fluo, na.rm = TRUE)
    shinyWidgets::numericRangeInput("fluo_range", "FLUO Range",
      value = .locked_val(input$fluo_range, rng), min = rng[1], max = rng[2])
  })
  output$eco_range_ui <- renderUI({
    df  <- uploaded_df()
    rng <- range(df$eco, na.rm = TRUE)
    shinyWidgets::numericRangeInput("eco_range", "ECO Range",
      value = .locked_val(input$eco_range, rng), min = rng[1], max = rng[2])
  })
  
  data <- reactive({
    req(input$file)
    df <- filtered_df()
    df
  })
  cross_dcr_nn_table <- reactive({
    # Only compute when the toggle is on
    if (!isTRUE(input$nn_cross_dcr)) return(data.frame())
    
    df <- plottable_data()
    if (!nrow(df)) return(data.frame())
    
    # We need a DCR column to define groups 1/2
    if (!("dcr" %in% names(df))) return(data.frame())
    
    # DCR ranges must exist and be finite
    g1 <- input$dcr_group1; g2 <- input$dcr_group2
    if (is.null(g1) || is.null(g2) || length(g1) != 2 || length(g2) != 2 ||
        !all(is.finite(c(g1, g2)))) return(data.frame())
    
    # Per-trace medians (geometry + dcr); use dcr_pooled when requested
    .col <- if (isTRUE(input$use_pooled_dcr) && "dcr_pooled" %in% names(df)) "dcr_pooled" else "dcr"
    med <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = stats::median(z, na.rm = TRUE),
        dcr_med = stats::median(.data[[.col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z), is.finite(dcr_med))
    
    if (!nrow(med)) return(data.frame())
    
    # Label groups strictly by DCR (this feature is explicitly “cross DCR”)
    med$grp <- dplyr::case_when(
      med$dcr_med >= g1[1] & med$dcr_med <  g1[2] ~ 1L,
      med$dcr_med >= g2[1] & med$dcr_med <= g2[2] ~ 2L,
      TRUE ~ NA_integer_
    )
    med <- med %>% dplyr::filter(grp %in% c(1L, 2L))
    if (!nrow(med) || length(unique(med$grp)) < 2) return(data.frame())
    
    # Split groups
    A <- med %>% dplyr::filter(grp == 1L)
    B <- med %>% dplyr::filter(grp == 2L)
    if (!nrow(A) || !nrow(B)) return(data.frame())
    
    # For each A, nearest neighbor in B (3D)
    coordsA <- as.matrix(A[, c("x","y","z")])
    coordsB <- as.matrix(B[, c("x","y","z")])
    nn <- tryCatch(RANN::nn2(data = coordsB, query = coordsA, k = 1),
                   error = function(e) NULL)
    if (is.null(nn)) return(data.frame())
    
    rng <- input$distance_range %||% c(0, Inf)
    
    out <- lapply(seq_len(nrow(A)), function(i) {
      j <- nn$nn.idx[i, 1]
      d <- nn$nn.dists[i, 1]
      if (!is.finite(d) || d < rng[1] || d > rng[2]) return(NULL)
      data.frame(
        from_x = A$x[i], from_y = A$y[i], from_z = A$z[i],
        to_x   = B$x[j], to_y   = B$y[j], to_z   = B$z[j],
        distance = d
      )
    })
    
    out <- do.call(rbind, out)
    if (is.null(out)) data.frame() else out
  })
  
    
  # Applies the COMMITTED 3D crop region (from trim_applied_rv, set by "Apply
  # crop"), gated by the master "Crop to region (3D)" checkbox. The live box/
  # sphere inputs drive only the preview, never the data, so tweaking coordinates
  # never crops until the user applies.
  apply_trimmer <- function(df, input) {
    if (!isTRUE(input$enable_trimmer) || nrow(df) == 0) return(df)
    ar <- trim_applied_rv()
    if (is.null(ar)) return(df)
    keep <- rep(TRUE, nrow(df))
    if (identical(ar$mode, "sphere")) {
      cc <- ar$c; r <- ar$r
      if (length(cc) == 3 && all(is.finite(c(cc, r))) && r >= 0)
        keep <- (df$x - cc[1])^2 + (df$y - cc[2])^2 + (df$z - cc[3])^2 <= r^2
    } else {
      if (!is.null(ar$x)) keep <- keep & df$x >= ar$x[1] & df$x <= ar$x[2]
      if (!is.null(ar$y)) keep <- keep & df$y >= ar$y[1] & df$y <= ar$y[2]
      if (!is.null(ar$z)) keep <- keep & df$z >= ar$z[1] & df$z <= ar$z[2]
    }
    df[keep, , drop = FALSE]
  }

  # Commit the current box/sphere inputs as the active crop region.
  trim_region_from_inputs <- function() {
    mode <- input$trim_mode %||% "box"
    num2 <- function(v) { v <- suppressWarnings(as.numeric(unlist(v))); if (length(v) == 2 && all(is.finite(v))) v else NULL }
    if (identical(mode, "sphere")) {
      cc <- suppressWarnings(as.numeric(c(input$trim_cx, input$trim_cy, input$trim_cz)))
      r  <- suppressWarnings(as.numeric(input$trim_r))
      if (length(cc) == 3 && all(is.finite(cc)) && is.finite(r) && r >= 0)
        list(mode = "sphere", c = cc, r = r) else NULL
    } else {
      list(mode = "box", x = num2(input$trim_x_range),
           y = num2(input$trim_y_range), z = num2(input$trim_z_range))
    }
  }
  observeEvent(input$trim_apply, {
    reg <- trim_region_from_inputs()
    if (is.null(reg)) { showNotification("Set the crop coordinates first.", type = "warning"); return() }
    trim_applied_rv(reg); plot_trigger(Sys.time())
    showNotification(if (identical(reg$mode, "sphere"))
      sprintf("3D crop applied: sphere r %.0f nm at (%.0f, %.0f, %.0f).", reg$r, reg$c[1], reg$c[2], reg$c[3])
      else "3D crop applied (box).", type = "message", duration = 4)
  }, ignoreInit = TRUE)
  observeEvent(input$trim_clear, {
    trim_applied_rv(NULL); plot_trigger(Sys.time())
    showNotification("3D crop cleared.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  # Redraw the preview shape a short moment after the coordinates change (only
  # while previewing in 3D), so the low-opacity box/sphere tracks the inputs
  # without re-rendering on every keystroke.
  trim_preview_tick <- shiny::debounce(reactive({
    list(input$trim_mode, input$trim_x_range, input$trim_y_range, input$trim_z_range,
         input$trim_cx, input$trim_cy, input$trim_cz, input$trim_r, input$trim_preview)
  }), 400)
  observeEvent(trim_preview_tick(), {
    if (isTRUE(input$enable_trimmer) && isTRUE(input$trim_preview) && identical(input$mode3d, "3D"))
      plot_trigger(Sys.time())
  }, ignoreInit = TRUE)

  # Add the low-opacity crop-region preview to a 3D figure from the LIVE inputs.
  add_3d_crop_preview <- function(fig) {
    if (!isTRUE(input$enable_trimmer) || !isTRUE(input$trim_preview)) return(fig)
    reg <- trim_region_from_inputs(); if (is.null(reg)) return(fig)
    col <- "#2c7fb8"
    if (identical(reg$mode, "sphere")) {
      cc <- reg$c; r <- reg$r
      nu <- 20L; nv <- 14L
      th <- seq(0, 2 * pi, length.out = nu); ph <- seq(0, pi, length.out = nv)
      g  <- expand.grid(iu = seq_len(nu), iv = seq_len(nv))
      vx <- cc[1] + r * sin(ph[g$iv]) * cos(th[g$iu])
      vy <- cc[2] + r * sin(ph[g$iv]) * sin(th[g$iu])
      vz <- cc[3] + r * cos(ph[g$iv])
      idx <- function(iu, iv) (iv - 1L) * nu + iu
      ii <- jj <- kk <- integer(0)
      for (iv in seq_len(nv - 1L)) for (iu in seq_len(nu - 1L)) {
        a <- idx(iu, iv); b <- idx(iu + 1L, iv); c2 <- idx(iu, iv + 1L); d <- idx(iu + 1L, iv + 1L)
        ii <- c(ii, a - 1L, b - 1L); jj <- c(jj, b - 1L, d - 1L); kk <- c(kk, c2 - 1L, c2 - 1L)
      }
      return(fig %>% plotly::add_trace(type = "mesh3d", x = vx, y = vy, z = vz,
        i = ii, j = jj, k = kk, opacity = 0.15, color = I(col), facecolor = rep(col, length(ii)),
        flatshading = TRUE, hoverinfo = "skip", showlegend = FALSE, showscale = FALSE, name = "crop preview"))
    }
    x <- reg$x; y <- reg$y; z <- reg$z
    if (is.null(x) || is.null(y) || is.null(z)) return(fig)
    vx <- c(x[1], x[2], x[2], x[1], x[1], x[2], x[2], x[1])
    vy <- c(y[1], y[1], y[2], y[2], y[1], y[1], y[2], y[2])
    vz <- c(z[1], z[1], z[1], z[1], z[2], z[2], z[2], z[2])
    ii <- c(0,0,4,4,0,0,3,3,0,0,1,1); jj <- c(1,2,5,6,1,5,2,6,3,7,2,6); kk <- c(2,3,6,7,5,4,6,7,7,4,6,5)
    fig %>% plotly::add_trace(type = "mesh3d", x = vx, y = vy, z = vz,
      i = ii, j = jj, k = kk, opacity = 0.15, color = I(col), facecolor = rep(col, length(ii)),
      flatshading = TRUE, hoverinfo = "skip", showlegend = FALSE, showscale = FALSE, name = "crop preview")
  }
  
  assign_analysis_group_local <- function(df, input) {
    # guard columns
    if (!("tid" %in% names(df))) stop("Missing 'tid' column in data.")
    .dcr_col <- if (isTRUE(input$use_pooled_dcr %||% FALSE) && "dcr_pooled" %in% names(df))
      "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL
    med_tbl <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        dcr_med  = if (!is.null(.dcr_col)) stats::median(.data[[.dcr_col]], na.rm = TRUE) else NA_real_,
        fluo_med = if ("fluo" %in% names(df)) stats::median(fluo, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
    df <- dplyr::left_join(df, med_tbl, by = "tid")
    
    mode <- (input$grouping_mode %||% "none")
    if (isTRUE(input$color_by_fluo)) mode <- "fluo"  # force fluo when override is on
    
    if (mode == "fluo") {
      df$analysis_group <- dplyr::case_when(
        is.finite(df$fluo_med) & round(df$fluo_med) == 1 ~ 1L,
        is.finite(df$fluo_med) & round(df$fluo_med) == 2 ~ 2L,
        TRUE ~ NA_integer_
      )
    } else if (mode == "dcr") {
      g1 <- input$dcr_group1; g2 <- input$dcr_group2
      if (!is.null(g1) && !is.null(g2) && length(g1) == 2 && length(g2) == 2 && all(is.finite(c(g1,g2)))) {
        v <- df$dcr_med
        df$analysis_group <- dplyr::case_when(
          is.finite(v) & v >= g1[1] & v <  g1[2] ~ 1L,
          is.finite(v) & v >= g2[1] & v <= g2[2] ~ 2L,
          TRUE ~ NA_integer_
        )
      } else {
        df$analysis_group <- NA_integer_
      }
    } else {
      df$analysis_group <- NA_integer_
    }
    df
  }
  # Put these with your other helpers ----
  normalize_scope <- function(s) {
    s <- as.character(s %||% "")
    if (s %in% c("within","between","all")) s else "between"
  }
  
  pair_passes_scope_2c <- function(g1, g2, scope = c("within","between")) {
    scope <- match.arg(scope)
    if (is.na(g1) || is.na(g2)) return(FALSE)
    if (!(g1 %in% c(1L,2L)) || !(g2 %in% c(1L,2L))) return(FALSE)
    if (scope == "within") return(g1 == g2)
    return(g1 != g2)  # "between"
  }
  
  pair_passes_scope_all <- function(g1, g2, scope = c("within","between","all")) {
    scope <- match.arg(scope)
    if (is.na(g1) || is.na(g2)) return(FALSE)
    if (!(g1 %in% c(1L,2L)) || !(g2 %in% c(1L,2L))) return(FALSE)
    if (scope == "within")  return(g1 == g2)
    if (scope == "between") return(g1 != g2)
    TRUE  # "all"
  }
  # make_two_color_tid_map() was a third copy of the DCR/fluo trace split, kept
  # alongside assign_two_color_group() and assign_analysis_group(). All three are
  # now resolve_channels() in helpers_smlm.R, which is the only place a channel
  # is defined.
  
  
  # two_color_nn_edges() lived here: a second NN builder that existed only to
  # colour lines by channel pair. build_nn_segments() returns that colouring now,
  # so the shared NN panel does the whole job and this one is gone.
  
  
  
  
  # ==========================================================================
  # plottable_data() - THE ANALYSIS-READY DATASET
  # ==========================================================================
  #
  # GOAL
  #   The single hand-off point between the filter pipeline and everything that
  #   analyses or draws the data. Nearly every analysis below begins with
  #   `df <- plottable_data()`. If a filter change is not reaching an analysis,
  #   this is the first place to look.
  #
  #   It is the SMLM counterpart of trk() in the tracking tab, but it does more:
  #   as well as delivering the filtered data it performs the CLUSTERING, so
  #   every consumer sees the same cluster / cluster2 assignment. That is
  #   deliberate - if each analysis clustered independently they could disagree.
  #
  # WHAT IT ADDS ON TOP OF filtered_df()
  #   1. apply_trimmer() ..... optional outlier/edge trimming
  #   2. 2D crop ............. applies the zoom rectangle SERVER-SIDE, so the
  #                            crop affects the analyses too and behaves
  #                            identically in 2D and 3D (rather than relying on
  #                            a browser round-trip that only changes the view)
  #   3. CLUSTERING .......... de-blinking and/or median clustering
  #                            (structure) via apply_clustering_pipeline().
  #                            Multiplex datasets are clustered PER CHANNEL.
  #   4. colour columns ...... per-point colour for the chosen colour mode
  #
  # WHY eventReactive(plot_trigger())
  #   Clustering a large localization table is expensive. Recomputing on every
  #   keystroke would make the app unusable, so the analyses refresh when Plot is
  #   pressed. Consequence worth knowing: after changing a clustering parameter
  #   you must press Plot before the numbers update.
  #
  # RETURNS the filtered data.frame plus (when enabled):
  #   cluster ..... de-blinked molecule ID (no 0: de-blinking assigns every trace)
  #   cluster2 .... secondary/structure ID (0 = noise)
  #   colour columns, and the multiplex channel/label columns when active
  # ==========================================================================
  plottable_data <- eventReactive(plot_trigger(), {
    df <- filtered_df()
    df <- apply_trimmer(df, input)
    # 2D zoom crop - applied server-side so it works in both 2D and 3D without
    # relying on the browser round-trip from update* calls.
    crop2d <- isolate(smlm_2d_crop_rv())
    if (!is.null(crop2d) && nrow(df) > 0)
      df <- df[df$x >= crop2d$x[1] & df$x <= crop2d$x[2] &
               df$y >= crop2d$y[1] & df$y <= crop2d$y[2], , drop = FALSE]
    if (nrow(df) == 0) return(df)
    # Remove rogue traces
    if (isTRUE(input$remove_rogues)) {
      rr <- remove_rogue_traces(
        df,
        rule   = input$rogue_rule %||% "factor_median",
        factor = input$rogue_factor %||% 2.5,
        min_n  = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT,
      )
      df <- rr$df
      #log text
      fit_log_rv(
        paste0(
          "Rogue removal: removed ", rr$removed_n, " trace(s). Rule: ", rr$note
        )
      )
    }
    
    if (nrow(df) == 0) return(df)
    
    df <- assign_analysis_group(df, input)
    
    # Safety check
    if (nrow(df) == 0) return(df)
    
    # --- Precompute per-trace median DCR (used by unmixer + optional color scale) ---
    if ("dcr" %in% names(df)) {
      tid_dcr <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(tid_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
      df <- df %>% dplyr::left_join(tid_dcr, by = "tid")
    }
    
    # 2) Base color assignment (use a single name everywhere)
    qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
    base_colors   <- unlist(mapply(RColorBrewer::brewer.pal,
                                   qual_col_pals$maxcolors,
                                   rownames(qual_col_pals)))
    
    # ── CLUSTERING: one call, shared with every batch tool ──────────────────
    # This block used to inline both pipelines - the single-colour one and a
    # second copy for 2-colour - which is why the old comments kept insisting
    # the two "must agree". They cannot disagree any more: the plot and the
    # batch tools now run the SAME function, and channels (none / 2-colour split
    # / multiplex) come from resolve_channels() inside it.
    df <- apply_clustering_pipeline(df, input)
    # Colour is presentation, so it is applied here rather than inside the cached
    # chain - see attach_channel_colors(). This is what keeps the colour pickers
    # out of the prepared-data cache key.
    df <- attach_channel_colors(df, input)
    n_ch <- length(unique(df$channel[!is.na(df$channel)]))

    # ── Colours ─────────────────────────────────────────────────────────────
    # Colouring is a display concern, so it stays here rather than moving into
    # the pipeline.
    if (isTRUE(input$enable_clustering)) {
      pal <- rep(base_colors, length.out = max(1, max(df$cluster, na.rm = TRUE) + 1))
      df$color <- ifelse(df$cluster == 0, "#bdbdbd",
                         pal[(df$cluster %% length(pal)) + 1])
    } else {
      extended_palette_tid <- rep(base_colors, length.out = nrow(df))
      df$color <- extended_palette_tid[
        (as.numeric(as.factor(df$tid)) %% length(extended_palette_tid)) + 1]
    }

    if (isTRUE(input$enable_cluster2)) {
      # Members per Median Cluster, counted in whatever unit fed the clustering.
      id_col <- if (isTRUE(input$enable_clustering)) "cluster" else "tid"
      c2_counts <- df %>%
        dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
        dplyr::group_by(cluster2) %>%
        dplyr::summarise(n_in_cluster2 = dplyr::n_distinct(.data[[id_col]]), .groups = "drop")
      df <- df %>%
        dplyr::select(-dplyr::any_of("n_in_cluster2")) %>%
        dplyr::left_join(c2_counts, by = "cluster2")
      df$n_in_cluster2 <- as.integer(df$n_in_cluster2)

      max_c2 <- suppressWarnings(max(df$cluster2, na.rm = TRUE))
      if (is.finite(max_c2) && max_c2 > 0) {
        pal2 <- rep(base_colors, length.out = max_c2 + 1)
        df$cluster2_color <- ifelse(is.na(df$cluster2) | df$cluster2 == 0, NA_character_,
                                    pal2[(df$cluster2 %% length(pal2)) + 1])
      } else {
        df$cluster2_color <- NA_character_
      }
    } else {
      df$cluster2_color <- NA_character_
      df$n_in_cluster2  <- NA_integer_
    }
    df$curvature <- NA_real_

    # Channel colours override the per-cluster palette, because when channels
    # exist the thing the user is looking for is which channel a point is in.
    # Numeric colour modes (z / time / density) still override this downstream.
    if (n_ch > 1L && "channel_color" %in% names(df) && any(!is.na(df$channel_color)))
      df$color <- ifelse(is.na(df$channel_color), df$color, df$channel_color)

    if (!"two_color_group" %in% names(df)) df$two_color_group <- NA_integer_
    # cluster_tc1 / cluster_tc2 are gone: they were aliases of cluster / cluster2
    # kept while the two pipelines were being merged, and every reader now uses
    # the real columns.
    
    
    # If requested, (re)attach TID-median DCR for hover text
    if (isTRUE(input$color_tid_by_median_dcr) && "dcr" %in% names(df)) {
      tid_dcr <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(tid_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
      df <- df %>% dplyr::left_join(tid_dcr, by = "tid")
    }
    # median DCR to the summarized df (trace or cluster) so hover is stable
    if (input$show_medians) {
      if (input$enable_clustering) {
        # cluster-level medians for DCR
        cl_dcr <- df %>%
          dplyr::filter(!is.na(cluster) & cluster != 0) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(cluster_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
        df <- df %>% dplyr::left_join(cl_dcr, by = "cluster")
      } else {
        # trace-level medians for DCR
        tid_dcr <- df %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(tid_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
        df <- df %>% dplyr::left_join(tid_dcr, by = "tid")
      }
    }
    # 7) Hover text 
    show_dcr_hover  <- isTRUE(input$enable_two_color) && ("dcr" %in% names(df))
    show_fluo_hover <- isTRUE(input$color_by_fluo) && ("fluo" %in% names(df))

    # A-#/B-# labels come straight off cluster/cluster2 now - there is one
    # clustering, so the "2-colour cluster" was never a different object.
    ch_lbl  <- if ("channel_label" %in% names(df)) df$channel_label else NULL
    tc1_lbl <- format_tc_id(df$cluster,  ch_lbl)
    tc2_lbl <- format_tc_id(df$cluster2, ch_lbl)
    
    df <- df %>% dplyr::mutate(
      hover_text = paste0(
        if ("tid" %in% names(.)) paste0("Trace ID: ", tid) else "",
        "<br>X: ", round(x, 2),
        "<br>Y: ", round(y, 2),
        if ("z" %in% names(.)) paste0("<br>Z: ", round(z, 2)) else "",
        if (!all(is.na(cluster)))  paste0("<br>De-blinked ID: ", cluster) else "",
        if (!all(is.na(cluster2))) paste0(
          "<br>Median Cluster: ", cluster2,
          ifelse(!is.na(n_in_cluster2), paste0(" (n=", n_in_cluster2, ")"), "")
        ) else "",

        # show group label only in the mode that defines it
        if (show_fluo_hover && "fluo_med" %in% names(.))
          paste0("<br>FLUO group: ", suppressWarnings(as.integer(round(fluo_med))))
        else "",
        
        # DCR hover text only when two-color is enabled
        if (show_dcr_hover && "tid_median_dcr" %in% names(.))
          paste0("<br>TID median DCR: ", round(tid_median_dcr, 2))
        else "",
        
        if (show_dcr_hover && "cluster_median_dcr" %in% names(.))
          paste0("<br>Cluster median DCR: ", round(cluster_median_dcr, 2))
        else "",
        
        if (show_dcr_hover && "dcr" %in% names(.))
          paste0("<br>Raw DCR: ", round(dcr, 2))
        else "",
        ifelse(
          isTRUE(input$enable_two_color) & !is.na(tc1_lbl),
          paste0("<br>Tc-De-blinked ID: ", tc1_lbl),
          ""),
        ifelse(
          isTRUE(input$enable_two_color) & !is.na(tc2_lbl),
          paste0("<br>Tc-Median Cluster: ", tc2_lbl),
          ""
        )
        )
    )

    # Multiplex: append the channel label to hover when present
    if ("mplex_label" %in% names(df))
      df$hover_text <- paste0(df$hover_text, "<br>Channel: ", df$mplex_label)

    # Derive globally-unique pcid / scid using the per-file random tag
    .tag <- if (".dataset_tag" %in% names(df)) df$.dataset_tag[[1L]] else "0000"
    df$pcid <- ifelse(
      !is.na(df$cluster) & df$cluster > 0L,
      paste0(.tag, "-", sprintf("%04d", df$cluster)),
      NA_character_
    )
    df$scid <- ifelse(
      !is.na(df$cluster2) & df$cluster2 > 0L,
      paste0(.tag, "-", sprintf("%04d", df$cluster2)),
      NA_character_
    )

    df
  })
  output$hist_range_ui <- renderUI({
    req(input$file, input$hist_variable)
    
    # Use RAW data to define min/max so the slider doesn't shrink/grow with filters
    df <- uploaded_df()
    varname <- input$hist_variable
    
    # Always reset to the full data range on re-render (new upload or variable change).
    # Preserving a stale curr value across uploads causes the histogram to appear pre-filtered.
    if (varname == "tid") {
      tlens <- as.numeric(table(df$tid))
      if (!length(tlens)) return(NULL)
      min_v <- min(tlens, na.rm = TRUE)
      max_v <- max(tlens, na.rm = TRUE)
      shinyWidgets::numericRangeInput(
        "hist_range", "Histogram Range:",
        value = c(min_v, max_v), min = min_v, max = max_v
      )
    } else {
      vals <- df[[varname]]
      if (is.null(vals) || !is.numeric(vals) || !length(vals)) return(NULL)
      vals <- vals[is.finite(vals)]
      if (!length(vals)) return(NULL)
      min_v <- min(vals, na.rm = TRUE)
      max_v <- max(vals, na.rm = TRUE)
      shinyWidgets::numericRangeInput(
        "hist_range", "Histogram Range:",
        value = c(min_v, max_v), min = min_v, max = max_v
      )
    }
  })
  # Enhanced NN table with all required metadata
  nn_table <- reactive({
    req(plot_trigger())
    df <- plottable_data()
    
    scope <- if (is.null(input$nn_group_scope)) "all" else input$nn_group_scope
    
    if (isTRUE(input$enable_clustering)) {
      # cluster medians
      df_geo <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          id             = dplyr::first(cluster),
          tid            = dplyr::first(tid),
          x              = stats::median(x, na.rm = TRUE),
          y              = stats::median(y, na.rm = TRUE),
          z              = stats::median(z, na.rm = TRUE),
          color = dplyr::first(color),     
          cluster2_color = dplyr::first(cluster2_color),
          cluster2       = dplyr::first(cluster2),
          analysis_group = dplyr::first(analysis_group),
          channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
          channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
          channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups = "drop"
        )
      id_label <- "cluster"
    } else {
      # tid medians
      df_geo <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          id             = dplyr::first(tid),
          x              = stats::median(x, na.rm = TRUE),
          y              = stats::median(y, na.rm = TRUE),
          z              = stats::median(z, na.rm = TRUE),
          color = dplyr::first(color),
          cluster2_color = dplyr::first(cluster2_color),
          cluster        = dplyr::first(cluster),
          cluster2       = dplyr::first(cluster2),
          analysis_group = dplyr::first(analysis_group),
          channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
          channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
          channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups = "drop"
        )
      id_label <- "tid"
    }

    if (nrow(df_geo) < 2) return(df_geo[0, ])

    coords <- as.matrix(df_geo[, c("x","y","z")])
    pm <- .nn_pairs_scoped(df_geo, coords, input)
    if (!nrow(pm)) return(df_geo[0, ])

    rows <- list()
    for (r in seq_len(nrow(pm))) {
      i <- pm[r, 1]; ni <- pm[r, 2]
      {
        # ---- HERE: apply within/between/all filter
        g1 <- df_geo$analysis_group[i]
        g2 <- df_geo$analysis_group[ni]
        if (!pair_passes_scope(g1, g2, scope)) next

        dist <- pm[r, 3]
        if (dist < input$distance_range[1] || dist > input$distance_range[2]) next

        # within-cluster2 filter
        if (isTRUE(input$nn_within_cluster2) && isTRUE(input$enable_cluster2)) {
          c2i <- df_geo$cluster2[i]; c2ni <- df_geo$cluster2[ni]
          if (is.na(c2i) || is.na(c2ni) || c2i == 0 || c2ni == 0 || c2i != c2ni) next
        }

        rows[[length(rows) + 1]] <- data.frame(
          source_id  = df_geo$id[i],
          target_id  = df_geo$id[ni],
          distance   = round(dist, 2),
          source_x   = df_geo$x[i],  source_y = df_geo$y[i],  source_z = df_geo$z[i],
          target_x   = df_geo$x[ni], target_y = df_geo$y[ni], target_z = df_geo$z[ni],
          source_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[i] else NA,
          target_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[ni] else NA,
          source_cluster2 = df_geo$cluster2[i],
          target_cluster2 = df_geo$cluster2[ni],
          stringsAsFactors = FALSE
        )
      }
    }
    out <- if (length(rows)) dplyr::bind_rows(rows) else df_geo[0, ]
    names(out)[names(out) == "source_id"] <- paste0("source_", id_label)
    names(out)[names(out) == "target_id"] <- paste0("target_", id_label)
    out
  })
  
  output$download_nn_table <- downloadHandler(
    filename = function() paste0("nn_distances_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      out <- list(distances = nn_table())
      # Second sheet when channels exist: one row per drawn pair, carrying which
      # channel each end belongs to. This replaces the separate "2-color NN table"
      # export, which was the same thing computed by a second code path.
      ed <- tryCatch(build_nn_segments(plottable_data(), input, is3d = TRUE)$edges,
                     error = function(e) NULL)
      if (is.data.frame(ed) && nrow(ed) &&
          length(unique(stats::na.omit(c(ed$ch1, ed$ch2)))) > 1L)
        out$channel_pairs <- ed[, c("x1","y1","z1","x2","y2","z2",
                                    "dist","ch1","ch2","lab1","lab2")]
      openxlsx::write.xlsx(out, file, row.names = FALSE)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  # Secondary cluster median-to-median NN lines (green)
  secondary_nn_lines <- reactive({
    df <- plottable_data()
    req(isTRUE(input$enable_cluster2))
    
    # majority group per cluster2 (for within/between filtering)
    mode_int <- function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_integer_)
      as.integer(names(which.max(table(x))))
    }
    
    df_geo <- df %>%
      dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
      dplyr::group_by(cluster2) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = stats::median(z, na.rm = TRUE),
        color = dplyr::first(color),    
        analysis_group = mode_int(analysis_group),
        .groups = "drop"
      )
    if (nrow(df_geo) < 2) return(list())
    
    scope  <- if (is.null(input$nn_group_scope)) "all" else input$nn_group_scope
    coords <- as.matrix(df_geo[, c("x","y","z")])
    tmp    <- safe_nn2(coords, k = input$nn_k_cluster2 + 1)
    if (is.null(tmp$nn)) return(list())
    
    nn    <- tmp$nn
    usedK <- max(0, tmp$k_eff - 1)
    
    lines <- list()
    seen  <- new.env(parent = emptyenv())
    for (i in seq_len(nrow(coords))) {
      for (j in 1:usedK) {
        ni <- nn$nn.idx[i, j + 1]
        if (is.na(ni)) next
        
        # de-duplicate undirected pairs
        key <- paste0(min(i, ni), "_", max(i, ni))
        if (!is.null(seen[[key]])) next
        seen[[key]] <- TRUE
        
        # apply within/between/all filter
        g1 <- df_geo$analysis_group[i]
        g2 <- df_geo$analysis_group[ni]
        if (!pair_passes_scope(g1, g2, scope)) next
        
        dist <- as.numeric(nn$nn.dists[i, j + 1])
        if (dist < input$distance_range_cluster2[1] || dist > input$distance_range_cluster2[2]) next

        nn_off <- suppressWarnings(as.numeric(input$c2_label_offset %||% 5))
        if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
        midpoint <- perp_offset_xy(coords[i, ], coords[ni, ], nn_off)
        lines[[length(lines) + 1]] <- list(
          x = c(coords[i, 1], coords[ni, 1]),
          y = c(coords[i, 2], coords[ni, 2]),
          z = c(coords[i, 3], coords[ni, 3]),
          text = round(dist, 2),
          midpoint = midpoint
        )
      }
    }
    lines
  })
  

  
  # 3) Plane distance table ------------ 
  output$download_plane_distance_table <- downloadHandler(
    filename = function() paste0("cluster2_xy_plane_distances_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      df <- plottable_data()
      tbl <- compute_cluster2_plane_distances(df)
      openxlsx::write.xlsx(
        x         = tbl,
        file      = file,
        sheetName = "Plane distances",
        rowNames  = FALSE
      )
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  # 4) Cluster2 metrics ---------------- 
  output$download_cluster2_table <- downloadHandler(
    filename = function() paste0("cluster2_fit_metrics_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      tbl <- cluster2_table_store()
      if (is.null(tbl) || nrow(tbl) == 0) stop("No Median Cluster metrics to download.")
      openxlsx::write.xlsx(
        x         = tbl,
        file      = file,
        sheetName = "Median Cluster metrics",
        rowNames  = FALSE
      )
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  # 5) Cluster2 NN table --------------- 
  output$download_nn_table2 <- downloadHandler(
    filename = function() paste0("nn_distances_cluster2_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      openxlsx::write.xlsx(
        x         = nn_table_cluster2(),
        file      = file,
        sheetName = "NN Median Cluster",
        rowNames  = FALSE
      )
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  # 6) Curvature table -----------------
  output$download_curvature_table <- downloadHandler(
    filename = function() paste0("cluster2_curvature_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      df <- plottable_data()
      curvature_table <- build_curvature_table(df)
      if (nrow(curvature_table) == 0) stop("No curvature data to export.")
      openxlsx::write.xlsx(
        x         = curvature_table,
        file      = file,
        sheetName = "Curvature",
        rowNames  = FALSE
      )
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  df_full_raw <- reactive({
    data()
  })
  
  output$scatter3D <- renderPlotly({
    req(input$file)
    req(plot_trigger(), input$mode3d == "3D")
    
    df <- plottable_data()
    shiny::validate(shiny::need(nrow(df) > 0, "No points to display - please check filters or upload a valid file."))
    
    df_raw <- uploaded_df()
    
    # Match coordinate shift if make_positive was applied
    if (isolate(input$make_positive)) {
      # Use same shift as used in plottable_data()
      min_x <- min(df_raw$x, na.rm = TRUE)
      min_y <- min(df_raw$y, na.rm = TRUE)
      min_z <- min(df_raw$z, na.rm = TRUE)
      
      df_raw <- df_raw %>%
        mutate(
          x = x - min_x,
          y = y - min_y,
          z = z - min_z
        )
    }
    
    # Use consistent axis limits for visual stability
    limits <- get_equal_limits(df$x, df$y, df$z)
    
    
    df$color <- as.character(df$color)
    
    # Axis theming in 3D
    bg_col   <- if (isTRUE(input$dark_mode)) "#000000" else "#FFFFFF"
    fg_col   <- if (isTRUE(input$dark_mode)) "#FFFFFF" else "#63666A"
    show_g   <- !isTRUE(input$hide_grid)
    grid_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
    zero_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.25)" else "rgba(0,0,0,0.25)"
    
    # Build the point table we actually plot
    pts <- build_pts_3d(df, input)
    
    # Compute colors for that table (default uses your existing palette in pts$color)
    pts$color <- build_colors_for_pts(pts, input, default_palette_col = "color")
    
    cmap <- color_mapping_for_pts(pts, input)
    
    # Plot the points ONCE
    fig <- plotly::plot_ly()
    
    id_col <- if (isTRUE(input$enable_clustering)) "cluster" else "tid"
    
    if (identical(id_col, "cluster")) {
      if (!("cluster" %in% names(pts)) || all(is.na(pts$cluster) | pts$cluster == 0)) {
        id_col <- "tid"
      }
    }
    
    # Key vector for points (use character keys ALWAYS)
    gid_key_pts <- as.character(pts[[id_col]])
    marker_hex  <- as.character(cmap$hex)
    # drop cluster 0
    keep <- !is.na(gid_key_pts) & nzchar(gid_key_pts)
    if (identical(id_col, "cluster")) keep <- keep & gid_key_pts != "0"
    
    gid_key_pts <- gid_key_pts[keep]
    marker_hex  <- marker_hex[keep]
    
    gid_hex_map <- tapply(marker_hex, gid_key_pts, function(x) x[1], simplify = TRUE)
    gid_hex_map <- unlist(gid_hex_map, use.names = TRUE)
    
    group_ids <- names(gid_hex_map)
    if (!length(group_ids)) return(fig)
    
    
  
    
    mplex_draw <- identical(input$color_mode %||% "default", "multiplex") &&
                  "mplex_label" %in% names(pts) && any(!is.na(pts$mplex_label))
    if (mplex_draw) {
      # one trace per channel so each gets a named, colored legend entry
      ptso  <- pts[order(pts$mplex_channel), , drop = FALSE]
      chans <- unique(ptso$mplex_channel[!is.na(ptso$mplex_channel)])
      for (chv in chans) {
        sub <- ptso[!is.na(ptso$mplex_channel) & ptso$mplex_channel == chv, , drop = FALSE]
        if (!nrow(sub)) next
        fig <- fig %>% plotly::add_trace(
          data = sub, x = ~x, y = ~y, z = ~z,
          type = "scatter3d", mode = "markers",
          name = sub$mplex_label[1],
          marker = list(size = input$median_point_size %||% 4, color = sub$mplex_color[1]),
          text = ~hover_text,
          hovertemplate = "%{text}<extra></extra>",
          hoverinfo = "text",
          showlegend = TRUE, inherit = FALSE
        )
      }
    } else if (isTRUE(cmap$showscale) && identical(cmap$kind, "scale")) {
      fig <- fig %>% plotly::add_trace(
        data = pts,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
        marker = list(
          size = input$median_point_size %||% 4,
          color = cmap$colorvar,
          colorscale = cmap$colorscale,
          cmin = cmap$cmin, cmax = cmap$cmax,
          showscale = TRUE,
          colorbar = list(
            title     = cmap$title,
            titlefont = list(color = fg_col),
            tickfont  = list(color = fg_col)
          )
        ),
        text = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        hoverinfo = "text",
        showlegend = TRUE, inherit = FALSE
      )
    } else {
      # hex colors (default / unicolor / degenerate numeric)
      hex <- if (!is.null(cmap$hex)) cmap$hex else rep("#bdbdbd", nrow(pts))
      fig <- fig %>% plotly::add_trace(
        data = pts,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        marker = list(size = input$median_point_size %||% 4, color = hex),
        text = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        hoverinfo = "text",
        name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
        showlegend = TRUE, inherit = FALSE
      )
    }
    fig <- fig %>% plotly::layout(
      title = "3D Scatter Plot",
      scene = list(
        xaxis = list(title = "X", range = limits$x),
        yaxis = list(title = "Y", range = limits$y),
        zaxis = list(title = "Z", range = limits$z),
        aspectmode = "manual",
        aspectratio = list(x = 1, y = 1, z = 1)
      ),
      showlegend = TRUE,
      legend = list(
        font = list(color = fg_col),
        orientation = "h",
        x = 0, y = 0.97,
        xanchor = "left", yanchor = "bottom"
      )
    )

    # Low-opacity preview of the 3D crop region (box / sphere) from the live inputs.
    fig <- add_3d_crop_preview(fig)

    cluster2_medians <- NULL
    if (isTRUE(input$enable_cluster2) && "cluster2" %in% names(df)) {
      # Member unit follows de-blinking, matching the two-colour isosurface and
      # every other per-structure statistic (cluster2_member_col / smlm_nn_unit).
      # This used to be hardcoded to `tid`, so with de-blinking on the envelope
      # was built from blink events rather than molecules.
      iso_unit <- cluster2_member_col(df)
      cluster2_medians <- df %>%
        dplyr::filter(!is.na(cluster2), cluster2 != 0) %>%
        { if (identical(iso_unit, "cluster"))
            dplyr::filter(., !is.na(cluster), cluster > 0L) else . } %>%
        # Grouped by (unit, channel), not the unit alone. A de-blinked molecule
        # never spans channels, but a trace ID can in multiplex mode (that split
        # is per localization), and one median over both would sit between them.
        dplyr::group_by(.data[[iso_unit]],
                        .ch = if ("channel" %in% names(.)) channel else NA_integer_) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          cluster2 = dplyr::first(cluster2),
          # carried so the wrap can be keyed and coloured by channel (see
          # wrap_groups()); constant within a group by construction
          channel_label = if ("channel_label" %in% names(df))
                            dplyr::first(channel_label) else NA_character_,
          channel_color = if ("channel_color" %in% names(df))
                            dplyr::first(channel_color) else NA_character_,
          .groups = "drop"
        ) %>%
        dplyr::rename(channel = .ch) %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z), is.finite(cluster2), cluster2 != 0)
    }
    
    
    # PRIMARY ALPHA SHAPE RENDERING
    if (isTRUE(input$show_convex_hull) && nrow(df) > 3) {
      
      id_col <- if (isTRUE(input$enable_clustering)) "cluster" else "tid"
      if (!id_col %in% names(pts)) return(fig)   # since we now hull-wrap pts
      
      # group_ids + gid_hex_map must already exist (built from pts + cmap above)
      if (!exists("group_ids") || !length(group_ids) ||
          !exists("gid_hex_map") || !length(gid_hex_map)) {
        return(fig)
      }
      
      # No alpha radius is asked for: alpha_surface_3d() derives one per group
      # from that group's own point spacing. See its header in helpers_smlm.R for
      # why a single global radius was both meaningless (a 30 nm molecule and a
      # 500 nm structure need different radii) and inert (the old code drew the
      # whole tetrahedralization, which is identical at every alpha).
      for (gid_key in as.character(group_ids)) {

        subset_df <- pts[as.character(pts[[id_col]]) == gid_key, , drop = FALSE]
        if (nrow(subset_df) < 4) next

        srf <- alpha_surface_3d(as.matrix(subset_df[, c("x","y","z"), drop = FALSE]))
        if (is.null(srf)) next

        # name-based lookup (NOW guaranteed)
        hull_col <- unname(gid_hex_map[[as.character(gid_key)]])

        if (is.null(hull_col) || !nzchar(hull_col)) hull_col <- "#bdbdbd"

        fig <- fig %>% plotly::add_trace(
          type = "mesh3d",
          # x/y/z come from srf, NOT from subset_df: alpha_surface_3d may have
          # subsampled, and i/j/k index the points it returned.
          x = srf$x, y = srf$y, z = srf$z,
          i = srf$i, j = srf$j, k = srf$k,
          opacity = input$hull_opacity %||% 0.2,

          # Force solid mesh color (Plotly sometimes ignores `color` for mesh3d)
          color       = I(hull_col),
          facecolor   = rep(hull_col, srf$n_tri),         # per-triangle color
          vertexcolor = rep(hull_col, srf$n_pts),         # per-vertex color

          flatshading = TRUE,
          lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
          name = paste("Hull", gid_key),
          hoverinfo = "skip",
          showlegend = FALSE,
          inherit = FALSE
        )
      }
    }
    
    
    
    
    # --- 4B) Primary NN lines on plot (uses same k/range as the download table) ---
    # General NN
    # --- General NN (trace medians OR cluster medians) ---
    # 3D NN lines. This used to be a second, inline copy of build_nn_segments() -
    # the same geometry, pair selection, dedupe and range filtering written out
    # again - which is how the 2D and 3D views could disagree about which pairs
    # were neighbours. One builder now serves both.
    if (isTRUE(input$show_nn_lines) && nrow(df) > 1) {
      nn <- build_nn_segments(df, input, is3d = TRUE)
      if (!is.null(nn)) {
        ln_col <- if (isTRUE(input$dark_mode)) "white" else "black"
        ed <- nn$edges
        multi_ch <- !is.null(ed) && nrow(ed) &&
          length(unique(stats::na.omit(c(ed$ch1, ed$ch2)))) > 1L
        if (multi_ch) {
          # Line colour says which channels the pair spans: within-channel takes
          # that channel's colour, cross-channel is magenta.
          for (cc in unique(ed$color)) {
            sel <- which(ed$color == cc)
            fig <- fig %>% plotly::add_trace(
              x = as.numeric(t(cbind(ed$x1[sel], ed$x2[sel], NA))),
              y = as.numeric(t(cbind(ed$y1[sel], ed$y2[sel], NA))),
              z = as.numeric(t(cbind(ed$z1[sel], ed$z2[sel], NA))),
              type = "scatter3d", mode = "lines",
              line = list(color = cc, width = 3),
              hoverinfo = "none", showlegend = FALSE, inherit = FALSE
            )
          }
        } else {
          fig <- fig %>% plotly::add_trace(
            x = nn$seg_x, y = nn$seg_y, z = nn$seg_z,
            type = "scatter3d", mode = "lines",
            line = list(color = ln_col, width = 2),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
        if (length(nn$labs) && isTRUE(input$show_nn_labels)) {
          lab_col <- if (isTRUE(input$dark_mode)) "orange" else "red"
          nn_lfs <- suppressWarnings(as.numeric(input$nn_label_font_size %||% 12))
          if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
          fig <- fig %>% plotly::add_trace(
            x = nn$lab_x, y = nn$lab_y, z = nn$lab_z,
            type = "scatter3d", mode = "text",
            text = nn$labs, textfont = list(color = lab_col, size = nn_lfs),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      }
    }
    # ── 2-color alpha shapes (de-blinked molecules) ─────────────────────────
    # Only the per-molecule wrap lives here. The Median Cluster wrap and
    # isosurface are drawn by the Median Clustering controls, per channel (see
    # wrap_groups()); a second copy here was redundant and had drifted.
    if (isTRUE(input$enable_two_color) && isTRUE(input$show_tc1_alpha) &&
        all(c("x","y","z") %in% names(df))) {
      
      add_alpha_for <- function(cluster_col, alpha, opacity, min_pts) {
        ids <- sort(unique(df[[cluster_col]]))
        ids <- ids[is.finite(ids) & ids != 0]
        if (!length(ids)) return(fig)
        
        for (cid in ids) {
          sub <- df[df[[cluster_col]] == cid, , drop = FALSE]
          if (nrow(sub) < min_pts) next
          
          pts <- as.matrix(sub[, c("x","y","z")])
          pts <- pts[stats::complete.cases(pts), , drop = FALSE]
          if (nrow(pts) < min_pts) next
          
          msh <- make_alpha_mesh3d(pts, alpha = alpha)
          if (is.null(msh)) next
          
          fig <- fig %>% plotly::add_trace(
            type = "mesh3d",
            x = msh$x, y = msh$y, z = msh$z,
            i = msh$i, j = msh$j, k = msh$k,
            opacity = opacity,
            hoverinfo = "skip",
            showscale = FALSE,
            inherit = FALSE
          )
        }
        fig
      }
      
      if (isTRUE(input$show_tc1_alpha) && "cluster" %in% names(df)) {
        fig <- add_alpha_for(
          cluster_col = "cluster",
          alpha   = input$tc1_alpha %||% 120,
          opacity = input$tc1_alpha_opacity %||% 0.15,
          min_pts = input$tc1_alpha_min_pts %||% 30
        )
      }
      
    }
    if (!isTRUE(input$show_cluster2_isosurface)) cluster2_iso_diag(NULL)

    # ── cluster2 isosurface (3D) ────────────────────────────────────────────
    if (isTRUE(input$show_cluster2_isosurface) &&
        isTRUE(input$enable_cluster2) &&
        !is.null(cluster2_medians) &&
        nrow(cluster2_medians) > 0 &&
        isTRUE(input$mode3d == "3D")) {
      
      # Geometry from the committed snapshot (Plot Data); opacity stays live.
      ig       <- iso_geom()
      grid_res <- ig$c2_grid_res %||% 22
      pad      <- ig$c2_pad %||% 0
      op       <- input$cluster2_iso_opacity %||% 0.25
      # NULL = auto, i.e. 3 x localization precision floored at half the member
      # gap. Anchoring to the microscope rather than to the labelling is what
      # keeps two clusters comparable; see cluster_density_field().
      iso_off  <- if (!isTRUE(ig$c2_manual_off)) NULL else {
        o <- suppressWarnings(as.numeric(ig$c2_offset %||% 30))
        if (!is.finite(o) || o <= 0) NULL else o
      }
      sigma    <- if (!isTRUE(ig$c2_manual_sigma)) NULL else {
        s <- suppressWarnings(as.numeric(ig$c2_sigma %||% 30))
        if (!is.finite(s) || s <= 0) NULL else s
      }

      # See the tc2 block: one dataset-wide precision, so every Median Cluster in
      # the view is wrapped on the same physical ruler.
      # Per-axis: MINFLUX z is usually TIGHTER than x/y (measured z:lateral 0.60
      # on the a real imaging file files), so an isotropic wrap over-wraps axially.
      iso_precision <- smlm_loc_precision_nm(df, per_axis = TRUE)
      c2_diag <- new.env(parent = emptyenv()); c2_diag$rows <- list()

      # One wrap per Median Cluster PER CHANNEL, so a surface can never enclose
      # molecules of two colours. See wrap_groups().
      groups <- wrap_groups(cluster2_medians, "cluster2")
      # Decided on the WHOLE dataset, before any single cluster is picked out:
      # a selected channel-B surface must still say "(B)" and still be red.
      multi_ch <- length(unique(stats::na.omit(groups$channel))) > 1L

      if (isTRUE(input$cluster2_iso_selected_only)) {
        sel <- .int_input(input$cluster2_iso_id, 1L, 1L)
        groups <- groups[groups$id == sel, , drop = FALSE]
      }

      if (nrow(groups)) {

      # build per-cluster color map using the same palette as cluster2 points
      qual_pals_iso <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
      base_colors_iso <- unlist(mapply(RColorBrewer::brewer.pal, qual_pals_iso$maxcolors, rownames(qual_pals_iso)))
      pal_iso <- rep(base_colors_iso, length.out = max(groups$id) + 1)

      for (gi in seq_len(nrow(groups))) {
        g   <- groups[gi, , drop = FALSE]
        cid <- g$id
        sub <- wrap_group_rows(cluster2_medians, "cluster2", g)
        if (nrow(sub) < 4) next

        pts <- as.matrix(sub[, c("x","y","z")])
        pts <- pts[stats::complete.cases(pts), , drop = FALSE]
        if (nrow(pts) < 4) next

        # Gaussian KDE over the member centroids, contoured at the level that
        # sits `offset` nm outside them. See cluster_density_field() for why the
        # field is a density rather than a nearest-distance transform, and why
        # the offset is anchored to localization precision rather than spacing.
        # Colour precedence: an explicit single colour the user picked, then the
        # channel the members belong to, then the cluster palette. Channel colour
        # is what makes this usable in multiplex mode, where the cluster palette
        # says nothing about which target a structure is.
        # By-channel colouring only means something with more than one channel;
        # with one it would paint every surface the same colour.
        iso_col <- if (isTRUE(input$cluster2_iso_unicolor)) {
          input$cluster2_iso_color %||% "#00CED1"
        } else if (multi_ch && isTRUE(input$cluster2_iso_by_channel)) {
          cc <- iso_channel_color(sub)
          if (is.na(cc)) pal_iso[(cid %% length(pal_iso)) + 1] else cc
        } else {
          pal_iso[(cid %% length(pal_iso)) + 1]
        }

        fig <- add_cluster_isosurface(
          fig, pts, color = iso_col, opacity = op,
          name = if (multi_ch && nzchar(g$label %||% ""))
                   paste0("Median Cluster ", cid, " (", g$label, ")")
                 else paste0("Median Cluster ", cid),
          grid_res = grid_res, pad = pad, sigma = sigma, offset = iso_off,
          precision = iso_precision, diag_env = c2_diag,
          sigma_scale = ig$c2_smooth_rel %||% 1, cache = iso_mesh_cache
        )
      }
      } # end else(!length(ids))
      cluster2_iso_diag(c2_diag$rows)
    }


    if (isTRUE(input$show_xy_plane)) {
      df_raw <- isolate(data())
      df_full <- plottable_data()  # or df_full <- data()
      limits <- get_equal_limits(df_full$x, df_full$y, df_full$z)
      df_for_plane <- if (isTRUE(input$plane_from_unclustered)) {
        df_full %>% filter(is.na(cluster2) | cluster2 == 0)
      } else {
        df_full
      }
      
      plane <- generate_xy_plane(df_for_plane, limits, z_offset = input$plane_offset)
      fig <- fig %>% add_surface(
        x = plane$x,
        y = plane$y,
        z = plane$z,
        opacity = 0.3,
        showscale = FALSE,
        colorscale = list(list(0, "gray"), list(1, "gray")),
        name = "XY Tilted Plane"
      )
    }
    if (isTRUE(input$enable_cluster2) && isTRUE(input$show_cluster2_metrics)) {
      df_full <- plottable_data()
      limits <- get_equal_limits(df_full$x, df_full$y, df_full$z)
      df_for_plane <- if (isTRUE(input$plane_from_unclustered)) {
        df_full %>% filter(is.na(cluster2) | cluster2 == 0)
      } else {
        df_full
      }
      
      plane <- generate_xy_plane(df_for_plane, limits, z_offset = input$plane_offset)
      
      membrane_mode <- isTRUE(input$membrane_protein_options)
      
      result <- add_cluster2_lines_and_angles(
        fig,
        df_full,
        plane = plane,
        use_pca_normal = input$show_xy_plane,
        cluster2_table_store = cluster2_table_store,
        enable_clustering = input$enable_clustering,
        helix_turns = input$helix_turns,
        dark_mode = input$dark_mode,
        bending_radius      = membrane_mode && isTRUE(input$show_bending_radius),
        show_plane_angles   = membrane_mode && isTRUE(input$show_cluster2_plane_angle),
        show_plane_distance = membrane_mode && isTRUE(input$show_cluster2_plane_distance),
        show_pca_axes = input$show_cluster2_pca_axes,
        show_info_labels = input$show_cluster2_labels,
        selected_fit_type   = input$fit_type,
        selected_fit_method = input$fit_method,
        compute_rmsd_metrics = isTRUE(input$compute_rmsd_metrics),
        fit_weighted        = isTRUE(input$fit_weighted),
        mle_center_mode     = input$mle_center_mode,
        mle_center_lambda   = input$mle_center_lambda,
        arc_span_hint       = input$arc_span_deg %||% 0
      )
      
      fig <- result$fig
      fit_log_rv(result$log)
    }
    # ---- GLOBAL SURFACE (from representative points) ----
    if (isTRUE(input$view_surface)) {
      df_surface <- plottable_data()
      
      # 1) optionally use only unclustered (must happen on raw points BEFORE aggregation)
      if (isTRUE(input$plane_from_unclustered)) {
        df_surface <- df_surface %>% dplyr::filter(is.na(cluster2) | cluster2 == 0)
        # If you meant "unclustered in cluster1", use:
        # df_surface <- df_surface %>% dplyr::filter(is.na(cluster) | cluster == 0)
      }
      
      # 2) exclude TIDs (still on raw points, before aggregation)
      ex_txt <- input$surface_exclude_tids %||% ""
      ex_ids <- unique(unlist(strsplit(ex_txt, "[,;\\s]+")))
      ex_ids <- ex_ids[nzchar(ex_ids)]
      
      if (length(ex_ids) && "tid" %in% names(df_surface)) {
        df_surface <- df_surface %>%
          dplyr::mutate(tid_chr = as.character(tid)) %>%
          dplyr::filter(!(.data$tid_chr %in% ex_ids))
      }

      # 2b) filter by DCR/fluo group (group 1 or 2 only)
      surface_group_sel <- input$surface_color_group %||% "all"
      if (surface_group_sel != "all") {
        # If plottable_data's analysis_group is all NA (e.g. grouping_mode="none")
        # but two-color is active, derive group tags via assign_two_color_group.
        has_valid_group <- "analysis_group" %in% names(df_surface) &&
                           any(!is.na(df_surface$analysis_group))
        if (!has_valid_group && isTRUE(input$enable_two_color)) {
          gmap <- tryCatch(assign_two_color_group(df_surface, input), error = function(e) NULL)
          if (is.data.frame(gmap) && "analysis_group" %in% names(gmap)) {
            df_surface <- df_surface %>%
              dplyr::select(-dplyr::any_of("analysis_group")) %>%
              dplyr::left_join(gmap %>% dplyr::select(tid, analysis_group), by = "tid")
          }
        }
        g_int <- suppressWarnings(as.integer(surface_group_sel))
        if ("analysis_group" %in% names(df_surface) && is.finite(g_int)) {
          df_surface <- df_surface %>%
            dplyr::filter(.data$analysis_group == g_int)
        }
      }

      # 3) build representative point set:
      #    - if de-blinking enabled -> de-blinked molecule medians
      #    - else -> tid medians
      if (isTRUE(input$enable_clustering) && "cluster" %in% names(df_surface)) {
        df_surface_pts <- df_surface %>%
          dplyr::filter(!is.na(cluster) & cluster != 0) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            z = stats::median(z, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        df_surface_pts <- df_surface %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            z = stats::median(z, na.rm = TRUE),
            .groups = "drop"
          )
      }
      
      df_surface_pts <- df_surface_pts %>%
        dplyr::mutate(across(c(x,y,z), as.numeric)) %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
      
      method <- input$surface_method %||% "blnk"
      if (method != "blnk" && nrow(df_surface_pts) >= 10) {
        
        if (method == "loess") {
          if (nrow(df_surface_pts) > 100000) df_surface_pts <- dplyr::slice_sample(df_surface_pts, n = 100000)
          
          surf <- tryCatch(
            create_global_surface(
              df_surface_pts,
              grid_n = input$surface_grid_n %||% 120,
              span   = input$surface_span   %||% 0.25
            ),
            error = function(e) NULL
          )
          
          if (!is.null(surf)) {
            fig <- fig %>% plotly::add_surface(
              x = surf$x, y = surf$y, z = surf$z,
              opacity = input$surface_opacity %||% 0.35,
              showscale = TRUE,
              colorscale = plotly_colorscale(
                input$surface_palette %||% "Viridis",
                reverse = isTRUE(input$surface_palette_reverse)
              ),
              colorbar = list(
                title     = "Surface z (nm)",
                titlefont = list(color = fg_col),
                tickfont  = list(color = fg_col)
              ),
              name = "Surface (LOESS)",
              hoverinfo = "skip",
              hovertemplate = ""
            )
          }
          
        } else if (method == "poly") {
          poly <- fit_poly_surface(
            df_surface_pts,
            degree = as.integer(input$poly_degree %||% 2L),
            grid_n  = as.integer(input$surface_grid_n %||% 120L)
          )
          if (!is.null(poly)) {
            fig <- fig %>% plotly::add_surface(
              x = poly$x, y = poly$y, z = poly$z,
              opacity = input$surface_opacity %||% 0.35,
              showscale = TRUE,
              colorscale = plotly_colorscale(
                input$surface_palette %||% "Viridis",
                reverse = isTRUE(input$surface_palette_reverse)
              ),
              colorbar = list(
                title     = "Surface z (nm)",
                titlefont = list(color = fg_col),
                tickfont  = list(color = fg_col)
              ),
              name = "Surface (Polynomial)",
              hoverinfo = "skip",
              hovertemplate = ""
            )
          }
          
        } else if (method == "tps") {
          if (nrow(df_surface_pts) > 50000) df_surface_pts <- dplyr::slice_sample(df_surface_pts, n = 50000)
          
          grid_n <- as.integer(input$surface_grid_n %||% 120L)
          gx <- seq(min(df_surface_pts$x), max(df_surface_pts$x), length.out = grid_n)
          gy <- seq(min(df_surface_pts$y), max(df_surface_pts$y), length.out = grid_n)
          
          gr <- expand.grid(x = gx, y = gy)
          
          fit2s <- fields::Tps(
            x = as.matrix(df_surface_pts[, c("x","y")]),
            Y = df_surface_pts$z,
            lambda = input$tps_lambda %||% 0.01
          )
          
          zhat <- predict(fit2s, as.matrix(gr[, c("x","y")]))
          

          zmat <- matrix(zhat, nrow = length(gx), ncol = length(gy), byrow = FALSE)
          zmat <- t(zmat)  # now rows = y, cols = x
          
          surf <- list(x = gx, y = gy, z = zmat)
          
          if (!is.null(surf)) {
            fig <- fig %>% plotly::add_surface(
              x = surf$x, y = surf$y, z = surf$z,
              opacity = input$surface_opacity %||% 0.35,
              showscale = TRUE,
              colorscale = plotly_colorscale(
                input$surface_palette %||% "Viridis",
                reverse = isTRUE(input$surface_palette_reverse)
              ),
              colorbar = list(
                title     = "Surface z (nm)",
                titlefont = list(color = fg_col),
                tickfont  = list(color = fg_col)
              ),
              name = "Surface (TPS)",
              hoverinfo = "skip",
              hovertemplate = ""
            )
          }
        }
        
      }
    }
    
    
    if (isTRUE(input$show_cluster2_nn_lines) && isTRUE(input$enable_cluster2) && !is.null(df$cluster2)) {
      lines <- secondary_nn_lines()
      if (length(lines)) {
        seg_x <- c(); seg_y <- c(); seg_z <- c()
        lab_x <- c(); lab_y <- c(); lab_z <- c(); labs <- c()
        for (ln in lines) {
          seg_x <- c(seg_x, ln$x, NA_real_)
          seg_y <- c(seg_y, ln$y, NA_real_)
          seg_z <- c(seg_z, ln$z, NA_real_)
          lab_x <- c(lab_x, ln$midpoint[1]); lab_y <- c(lab_y, ln$midpoint[2]); lab_z <- c(lab_z, ln$midpoint[3]); labs <- c(labs, ln$text)
        }
        fig <- fig %>% plotly::add_trace(
          x = seg_x, y = seg_y, z = seg_z,
          type = "scatter3d", mode = "lines",
          line = list(color = "green", width = 2),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
        if (length(labs) && isTRUE(input$show_c2_nn_labels)) {
          nn_lfs <- suppressWarnings(as.numeric(input$c2_label_font_size %||% 12))
          if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
          fig <- fig %>% plotly::add_trace(
            x = lab_x, y = lab_y, z = lab_z,
            type = "scatter3d", mode = "text",
            text = labs, textfont = list(color = "green", size = nn_lfs),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      }
    }
    # ==== Cluster2 fit drawing (Alpha or Smooth) =====================
    if (isTRUE(input$enable_cluster2) && !is.null(input$fit_type) && input$fit_type != "none") {
      df_full <- plottable_data()
      df_cluster2 <- df_full %>% dplyr::filter(!is.na(cluster2) & cluster2 != 0)
      clusters <- df_cluster2 %>% dplyr::group_by(cluster2) %>% dplyr::group_split()
      if (!length(clusters)) {
        fit_log_rv("No Median Cluster structures available for fitting. Check Median Clustering DBSCAN settings (eps2/minpts2).")
      }
      cap_radius_for_cluster <- function(radius, X, cap_mult = 2.5) {
        if (!is.finite(radius) || radius <= 0) return(NA_real_)
        ctr <- colMeans(X, na.rm = TRUE)
        d <- sqrt(rowSums((X - matrix(ctr, nrow(X), 3, byrow = TRUE))^2))
        d <- d[is.finite(d)]
        if (!length(d)) return(radius)
        cap <- as.numeric(stats::quantile(d, 0.95, na.rm = TRUE)) * cap_mult
        if (!is.finite(cap) || cap <= 0) cap <- max(d, na.rm = TRUE)
        min(radius, cap)
      }
      cap_axes_for_cluster <- function(axes, X, cap_mult = 2.5) {
        ax <- as.numeric(axes)
        if (length(ax) != 3) return(rep(NA_real_, 3))
        ctr <- colMeans(X, na.rm = TRUE)
        d <- sqrt(rowSums((X - matrix(ctr, nrow(X), 3, byrow = TRUE))^2))
        d <- d[is.finite(d)]
        if (!length(d)) return(ax)
        cap <- as.numeric(stats::quantile(d, 0.95, na.rm = TRUE)) * cap_mult
        if (!is.finite(cap) || cap <= 0) cap <- max(d, na.rm = TRUE)
        pmax(1e-6, pmin(ax, cap))
      }
      # GOAL: triangulate the surface of a FIT (sphere / ellipsoid), given the
      # dense uniform point sampling of it that the fit functions return.
      # Those surfaces are convex by construction, so the convex hull IS the
      # exact triangulation - no alpha shape is needed or wanted here. (This
      # used to build an ashape3d first and select `triang[, 4] >= 2`; column 4
      # is `on.ch`, a 0/1 flag, so that test never passed and every call fell
      # through to this hull anyway, after paying for a 4000-point
      # tetrahedralization. Data clouds, which are NOT convex, are wrapped by
      # alpha_surface_3d() instead.)
      surface_triangles <- function(pts) {
        tryCatch({
          hh <- geometry::convhulln(pts, options = "Qt")
          as.matrix(hh[, 1:3, drop = FALSE])
        }, error = function(e) NULL)
      }
      build_tid_fit_medians <- function(d) {
        d %>%
          dplyr::group_by(.data$tid) %>%
          dplyr::summarise(
            x = stats::median(.data$x, na.rm = TRUE),
            y = stats::median(.data$y, na.rm = TRUE),
            z = stats::median(.data$z, na.rm = TRUE),
            spread = stats::mad(
              sqrt(
                (.data$x - stats::median(.data$x, na.rm = TRUE))^2 +
                (.data$y - stats::median(.data$y, na.rm = TRUE))^2 +
                (.data$z - stats::median(.data$z, na.rm = TRUE))^2
              ),
              na.rm = TRUE
            ),
            n_locs = dplyr::n(),
            .groups = "drop"
          ) %>%
          dplyr::filter(is.finite(.data$x), is.finite(.data$y), is.finite(.data$z))
      }
      
      # Cache MLE fit objects keyed by cluster_id so the distance-lines loop
      # below can reuse them without refitting the same data a second time.
      mle_fit_cache <- list()

      for (cluster_data in clusters) {
        cluster_id <- unique(cluster_data$cluster2)
        if (is.na(cluster_id) || nrow(cluster_data) < 4) next

        # --- choose color for this cluster’s fits ---
        fit_col <- {
          co <- unique(na.omit(as.character(cluster_data$cluster2_color)))
          if (length(co) && nzchar(co[1])) co[1] else "#00BFFF"   # default fallback
        }
        
        # medians used by MLE path weights, unchanged
        if (isTRUE(input$enable_clustering) && "cluster" %in% names(cluster_data)) {
          med <- cluster_data %>%
            dplyr::filter(!is.na(.data$cluster) & .data$cluster != 0) %>%
            dplyr::group_by(.data$cluster) %>%
            dplyr::summarise(
              mx     = stats::median(.data$x, na.rm = TRUE),
              my     = stats::median(.data$y, na.rm = TRUE),
              mz     = stats::median(.data$z, na.rm = TRUE),
              spread = if (dplyr::n() > 1)
                sqrt(mean((.data$x - mx)^2 + (.data$y - my)^2 + (.data$z - mz)^2, na.rm = TRUE)) else 0,
              n_locs = dplyr::n(),
              .groups = "drop"
            ) %>%
            dplyr::rename(x = mx, y = my, z = mz)
        } else {
          med <- cluster_data %>%
            dplyr::group_by(.data$tid) %>%
            dplyr::summarise(
              mx     = stats::median(.data$x, na.rm = TRUE),
              my     = stats::median(.data$y, na.rm = TRUE),
              mz     = stats::median(.data$z, na.rm = TRUE),
              spread = if (dplyr::n() > 1)
                sqrt(mean((.data$x - mx)^2 + (.data$y - my)^2 + (.data$z - mz)^2, na.rm = TRUE)) else 0,
              n_locs = dplyr::n(),
              .groups = "drop"
            ) %>%
            dplyr::rename(x = mx, y = my, z = mz)
        }

        min_fit_n <- if (identical(input$fit_type, "helix")) 5L else 4L
        if (nrow(med) < min_fit_n) {
          med <- build_tid_fit_medians(cluster_data)
        }
        if (nrow(med) < min_fit_n) next
        X <- as.matrix(med[, c("x","y","z")])
        w <- if (isTRUE(input$fit_weighted)) {
          s   <- med$spread
          n_l <- if ("n_locs" %in% names(med)) med$n_locs else rep(1L, nrow(med))
          s_ref <- stats::median(s[n_l > 1 & s > 0], na.rm = TRUE)
          if (!is.finite(s_ref) || s_ref <= 0) s_ref <- max(s[s > 0], na.rm = TRUE)
          if (!is.finite(s_ref) || s_ref <= 0) s_ref <- 1
          s_eff <- ifelse(n_l > 1 & s > 0, s, s_ref)
          n_l / s_eff^2
        } else NULL
        
        if (input$fit_method == "heuristic") {
          # --- Heuristic path ---
          if (input$fit_type == "sphere") {
            fit <- generate_sphere_from_cluster2(
              cluster_data,
              use_cluster1 = isTRUE(input$enable_clustering)
            )
            if (!is.null(fit)) {
              fit$radius <- cap_radius_for_cluster(fit$radius, X)
              if (!is.finite(fit$radius) || fit$radius <= 0) next
              pts <- fit$virtual_points
              surf_tris <- surface_triangles(pts)
              if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
                fig <- fig %>% add_trace(
                  type = "mesh3d",
                  x = pts[,1], y = pts[,2], z = pts[,3],
                  i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                  opacity = input$fit_opacity,
                  color = I(fit_col),
                  facecolor   = rep(fit_col, nrow(surf_tris)),
                  vertexcolor = rep(fit_col, nrow(pts)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Spherical Fit", cluster_id),
                  hoverinfo = "skip",
                  hovertemplate = "",
                  showlegend = TRUE
                )
              }
            }

          } else if (input$fit_type == "ellipsoid") {
            fit <- generate_ellipsoid_from_cluster2(
              cluster_data,
              use_cluster1 = isTRUE(input$enable_clustering)
            )
            if (!is.null(fit)) {
              fit$axes <- cap_axes_for_cluster(fit$axes, X)
              if (!all(is.finite(fit$axes)) || any(fit$axes <= 0)) next
              pts <- fit$virtual_points
              surf_tris <- surface_triangles(pts)
              if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
                fig <- fig %>% add_trace(
                  type = "mesh3d",
                  x = pts[,1], y = pts[,2], z = pts[,3],
                  i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                  opacity = input$fit_opacity,
                  color = I(fit_col),
                  facecolor   = rep(fit_col, nrow(surf_tris)),
                  vertexcolor = rep(fit_col, nrow(pts)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Ellipsoid Fit", cluster_id),
                  hoverinfo = "skip",
                  hovertemplate = "",
                  showlegend = TRUE
                )
              }
            }
            
          } else if (input$fit_type == "helix") {
            fit <- generate_helix_from_cluster2(
              cluster_data,
              n_points = 3000,
              turns = input$helix_turns,
              use_cluster1 = isTRUE(input$enable_clustering)
            )
            if (!is.null(fit)) {
              if (isTRUE(input$show_fit_tube)) {
                geo_h  <- fit_helix_geometry(X, turns = input$helix_turns)
                tube_r <- max(1, input$fit_tube_radius %||% 5)
                tmax_h <- 2 * pi * max(1, input$helix_turns)
                msh <- tryCatch(make_helix_tube_mesh(
                  center = geo_h$proj_center, axis = geo_h$axis,
                  radius = geo_h$radius, pitch = geo_h$pitch, phi0 = geo_h$phi0,
                  tmin = 0, tmax = tmax_h, tube_r = tube_r
                ), error = function(e) NULL)
                if (!is.null(msh)) {
                  fig <- fig %>% add_trace(
                    type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                    i = msh$i, j = msh$j, k = msh$k,
                    opacity = input$fit_opacity, color = I(fit_col),
                    facecolor = rep(fit_col, length(msh$i)),
                    vertexcolor = rep(fit_col, length(msh$x)),
                    flatshading = TRUE,
                    lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                    name = paste("Helix Fit", cluster_id),
                    showlegend = TRUE, inherit = FALSE,
                    hoverinfo = "skip", hovertemplate = ""
                  )
                }
              } else {
                hp <- fit$points
                fig <- fig %>% add_trace(
                  x = hp$x, y = hp$y, z = hp$z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = fit_col, width = 4),
                  name = paste("Helix Fit", cluster_id),
                  showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            }

          } else if (input$fit_type == "ring") {
            fit <- generate_ring_from_cluster2(
              cluster_data,
              n_points = 500,
              use_cluster1 = isTRUE(input$enable_clustering)
            )
            if (!is.null(fit)) {
              if (isTRUE(input$show_fit_tube)) {
                tube_r <- max(1, input$fit_tube_radius %||% 5)
                msh <- tryCatch(make_ring_tube_mesh(
                  fit$center, fit$normal, fit$radius, tube_r
                ), error = function(e) NULL)
                if (!is.null(msh)) {
                  fig <- fig %>% add_trace(
                    type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                    i = msh$i, j = msh$j, k = msh$k,
                    opacity = input$fit_opacity, color = I(fit_col),
                    facecolor = rep(fit_col, length(msh$i)),
                    vertexcolor = rep(fit_col, length(msh$x)),
                    flatshading = TRUE,
                    lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                    name = paste("Ring Fit", cluster_id),
                    showlegend = TRUE, inherit = FALSE,
                    hoverinfo = "skip", hovertemplate = ""
                  )
                }
              } else {
                rp <- fit$points
                fig <- fig %>% add_trace(
                  x = rp$x, y = rp$y, z = rp$z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = fit_col, width = 4),
                  name = paste("Ring Fit", cluster_id),
                  showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            }

          } else if (input$fit_type == "arc") {
            fit <- generate_arc_from_cluster2(
              cluster_data,
              n_points = 300,
              arc_span_hint = input$arc_span_deg %||% 0,
              use_cluster1 = isTRUE(input$enable_clustering)
            )
            if (!is.null(fit)) {
              if (isTRUE(input$show_fit_tube)) {
                tube_r <- max(1, input$fit_tube_radius %||% 5)
                msh <- tryCatch(make_arc_tube_mesh(
                  fit$center, fit$normal, fit$radius,
                  fit$t_center, fit$arc_half, tube_r
                ), error = function(e) NULL)
                if (!is.null(msh)) {
                  fig <- fig %>% add_trace(
                    type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                    i = msh$i, j = msh$j, k = msh$k,
                    opacity = input$fit_opacity, color = I(fit_col),
                    facecolor = rep(fit_col, length(msh$i)),
                    vertexcolor = rep(fit_col, length(msh$x)),
                    flatshading = TRUE,
                    lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                    name = paste("Arc Fit", cluster_id),
                    showlegend = TRUE, inherit = FALSE,
                    hoverinfo = "skip", hovertemplate = ""
                  )
                }
              } else {
                rp <- fit$points
                fig <- fig %>% add_trace(
                  x = rp$x, y = rp$y, z = rp$z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = fit_col, width = 4),
                  name = paste("Arc Fit", cluster_id),
                  showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            }
          }

        } else {
          # --------- MLE path with center handling ---------
          center_args <- build_center_args(X, input)
          
          if (input$fit_type == "sphere") {
            f <- do.call(mle_fit_sphere, c(list(X = X, w = w), center_args))
            mle_fit_cache[[as.character(cluster_id)]] <- f   # cache for distance-lines loop
            if (!is.null(f)) {
              f$radius <- cap_radius_for_cluster(f$radius, X)
              if (!is.finite(f$radius) || f$radius <= 0) next
              pts <- mle_sample_sphere(f$center, f$radius, n = 4000)
              surf_tris <- surface_triangles(pts)
              if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
                fig <- fig %>% add_trace(
                  type = "mesh3d",
                  x = pts[,1], y = pts[,2], z = pts[,3],
                  i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                  opacity = input$fit_opacity,
                  color   = I(fit_col),
                  facecolor   = rep(fit_col, nrow(surf_tris)),
                  vertexcolor = rep(fit_col, nrow(pts)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Sphere (MLE)", cluster_id),
                  hoverinfo = "skip",
                  hovertemplate = "",
                  showlegend = TRUE
                )
              }
            }

          } else if (input$fit_type == "ellipsoid") {
            f <- do.call(mle_fit_ellipsoid, c(list(X = X, w = w), center_args))
            mle_fit_cache[[as.character(cluster_id)]] <- f   # cache for distance-lines loop
            if (!is.null(f)) {
              f$axes <- cap_axes_for_cluster(f$axes, X)
              if (!all(is.finite(f$axes)) || any(f$axes <= 0)) next
              pts <- mle_sample_ellipsoid(f$center, f$axes, f$R, n = 4000)
              surf_tris <- surface_triangles(pts)
              if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
                fig <- fig %>% add_trace(
                  type = "mesh3d",
                  x = pts[,1], y = pts[,2], z = pts[,3],
                  i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                  opacity = input$fit_opacity,
                  color   = I(fit_col),
                  facecolor   = rep(fit_col, nrow(surf_tris)),
                  vertexcolor = rep(fit_col, nrow(pts)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Ellipsoid (MLE)", cluster_id),
                  showlegend = TRUE,
                  hoverinfo = "skip",
                  hovertemplate = ""
                )
              }
            }
            
          } else if (input$fit_type == "helix") {
            f <- do.call(mle_fit_helix, c(list(X = X, w = w, turns_hint = input$helix_turns), center_args))
            mle_fit_cache[[as.character(cluster_id)]] <- f
            if (!is.null(f)) {
              if (isTRUE(input$show_fit_tube)) {
                tube_r <- max(1, input$fit_tube_radius %||% 5)
                pad    <- 0.25
                tmin_h <- min(f$t) - 2*pi*pad
                tmax_h <- max(f$t) + 2*pi*pad
                msh <- tryCatch(make_helix_tube_mesh(
                  center = f$center, axis = f$axis,
                  radius = f$radius, pitch = f$pitch, phi0 = f$phi0,
                  tmin = tmin_h, tmax = tmax_h, tube_r = tube_r
                ), error = function(e) NULL)
                if (!is.null(msh)) {
                  fig <- fig %>% add_trace(
                    type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                    i = msh$i, j = msh$j, k = msh$k,
                    opacity = input$fit_opacity, color = I(fit_col),
                    facecolor = rep(fit_col, length(msh$i)),
                    vertexcolor = rep(fit_col, length(msh$x)),
                    flatshading = TRUE,
                    lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                    name = paste("Helix (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                    hoverinfo = "skip", hovertemplate = ""
                  )
                }
              } else {
                hp <- mle_sample_helix(f, n = 3000)
                fig <- fig %>% add_trace(
                  x = hp$V1, y = hp$V2, z = hp$V3,
                  type = "scatter3d", mode = "lines",
                  line = list(color = fit_col, width = 4),
                  name = paste("Helix (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            }

          } else if (input$fit_type == "ring") {
            f <- do.call(mle_fit_ring, c(list(X = X, w = w), center_args))
            mle_fit_cache[[as.character(cluster_id)]] <- f
            if (!is.null(f)) {
              if (isTRUE(input$show_fit_tube)) {
                tube_r <- max(1, input$fit_tube_radius %||% 5)
                msh <- tryCatch(make_ring_tube_mesh(
                  f$center, f$normal, f$radius, tube_r
                ), error = function(e) NULL)
                if (!is.null(msh)) {
                  fig <- fig %>% add_trace(
                    type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                    i = msh$i, j = msh$j, k = msh$k,
                    opacity = input$fit_opacity, color = I(fit_col),
                    facecolor = rep(fit_col, length(msh$i)),
                    vertexcolor = rep(fit_col, length(msh$x)),
                    flatshading = TRUE,
                    lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                    name = paste("Ring (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                    hoverinfo = "skip", hovertemplate = ""
                  )
                }
              } else {
                rp <- mle_sample_ring(f, n = 500)
                fig <- fig %>% add_trace(
                  x = rp$x, y = rp$y, z = rp$z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = fit_col, width = 4),
                  name = paste("Ring (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            }

          } else if (input$fit_type == "arc") {
            f <- do.call(mle_fit_arc, c(list(X = X, w = w, arc_span_hint = input$arc_span_deg %||% 0), center_args))
            mle_fit_cache[[as.character(cluster_id)]] <- f
            if (!is.null(f)) {
              if (isTRUE(input$show_fit_tube)) {
                tube_r <- max(1, input$fit_tube_radius %||% 5)
                msh <- tryCatch(make_arc_tube_mesh(
                  f$center, f$normal, f$radius,
                  f$t_center, f$arc_half, tube_r
                ), error = function(e) NULL)
                if (!is.null(msh)) {
                  fig <- fig %>% add_trace(
                    type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                    i = msh$i, j = msh$j, k = msh$k,
                    opacity = input$fit_opacity, color = I(fit_col),
                    facecolor = rep(fit_col, length(msh$i)),
                    vertexcolor = rep(fit_col, length(msh$x)),
                    flatshading = TRUE,
                    lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                    name = paste("Arc (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                    hoverinfo = "skip", hovertemplate = ""
                  )
                }
              } else {
                rp <- mle_sample_arc(f, n = 300)
                fig <- fig %>% add_trace(
                  x = rp$x, y = rp$y, z = rp$z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = fit_col, width = 4),
                  name = paste("Arc (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            }
          }
        }
      }
    }
    
    if (!is.null(input$fit_type) &&
        input$fit_type != "none" &&
        isTRUE(input$show_fit_distances)) {
      
      df_full <- plottable_data()
      df_cluster2 <- df_full %>% dplyr::filter(!is.na(cluster2) & cluster2 != 0)
      clusters <- df_cluster2 %>% dplyr::group_by(cluster2) %>% dplyr::group_split()
      
      for (cluster_data in clusters) {
        # medians for distance computation
        if (isTRUE(input$enable_clustering)) {
          # Use de-blinked molecule medians
          med <- cluster_data %>%
            dplyr::filter(!is.na(cluster) & cluster != 0) %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarise(
              mx     = median(x), my = median(y), mz = median(z),
              spread = if (dplyr::n() > 1)
                sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2)) else 0,
              n_locs = dplyr::n(),
              .groups = "drop"
            ) %>%
            dplyr::rename(x = mx, y = my, z = mz)
        } else {
          # Default: TID medians
          med <- cluster_data %>%
            dplyr::group_by(tid) %>%
            dplyr::summarise(
              mx     = median(x), my = median(y), mz = median(z),
              spread = if (dplyr::n() > 1)
                sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2)) else 0,
              n_locs = dplyr::n(),
              .groups = "drop"
            ) %>%
            dplyr::rename(x = mx, y = my, z = mz)
        }
        if (nrow(med) < 4) next
        X <- as.matrix(med[, c("x","y","z")])
        w <- if (isTRUE(input$fit_weighted)) {
          s   <- med$spread
          n_l <- if ("n_locs" %in% names(med)) med$n_locs else rep(1L, nrow(med))
          s_ref <- stats::median(s[n_l > 1 & s > 0], na.rm = TRUE)
          if (!is.finite(s_ref) || s_ref <= 0) s_ref <- max(s[s > 0], na.rm = TRUE)
          if (!is.finite(s_ref) || s_ref <= 0) s_ref <- 1
          s_eff <- ifelse(n_l > 1 & s > 0, s, s_ref)
          n_l / s_eff^2
        } else NULL
        
        if (input$fit_method == "heuristic") {
          if (input$fit_type == "sphere") {
            g <- fit_sphere_geometry(X)
            r <- cap_radius_for_cluster(g$radius, X)
            if (is.finite(r) && r > 0) {
              pr <- mle_project_sphere(X, g$center, r)
              med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
              med$distance_to_fit <- pr$dist
              fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
            }
          } else if (input$fit_type == "ellipsoid") {
            g <- fit_ellipsoid_geometry(X)
            if (!is.null(g)) {
              g$axes <- cap_axes_for_cluster(g$axes, X)
              if (all(is.finite(g$axes)) && all(g$axes > 0)) {
                pr <- mle_project_ellipsoid(X, g$center, g$axes, g$R)
                med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
                med$distance_to_fit <- pr$dist
                fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
              }
            }
          } else if (input$fit_type == "helix") {
            g <- fit_helix_geometry(X, turns = input$helix_turns, n_points = 3000)
            t_h <- .project_t_helix(X, g$proj_center, g$axis, g$radius, g$pitch, g$phi0,
                                     check_neighbors = TRUE)
            H_h <- t(vapply(t_h, .helix_point, numeric(3),
                            c = g$proj_center, u = g$axis, r = g$radius, p = g$pitch, phi0 = g$phi0))
            med$nearest_x <- H_h[,1]; med$nearest_y <- H_h[,2]; med$nearest_z <- H_h[,3]
            med$distance_to_fit <- sqrt(rowSums((X - H_h)^2))
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)

          } else if (input$fit_type == "ring") {
            g  <- fit_ring_geometry(X)
            pr <- mle_project_ring(X, g$center, g$normal, g$radius)
            med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
            med$distance_to_fit <- pr$dist
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)

          } else if (input$fit_type == "arc") {
            g  <- fit_arc_geometry(X, arc_span_hint = input$arc_span_deg %||% 0)
            pr <- mle_project_arc(X, g$center, g$normal, g$radius, g$t_center, g$arc_half)
            med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
            med$distance_to_fit <- pr$dist
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
          }

        } else {
          # --- MLE: reuse cached fit from the drawing loop above (same data, same settings).
          # Falls back to refitting only if the cache is missing (e.g. show_fit_distances
          # was enabled on a render where the shape loop was skipped).
          cl_key      <- as.character(unique(cluster_data$cluster2))
          center_args <- build_center_args(X, input)

          if (input$fit_type == "sphere") {
            f <- mle_fit_cache[[cl_key]] %||%
                   do.call(mle_fit_sphere, c(list(X = X, w = w), center_args))
            if (!is.null(f)) {
              f$radius <- cap_radius_for_cluster(f$radius, X)
              if (!is.finite(f$radius) || f$radius <= 0) next
              pr <- mle_project_sphere(X, f$center, f$radius)
              med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
              med$distance_to_fit <- pr$dist
              fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
            }

          } else if (input$fit_type == "ellipsoid") {
            f <- mle_fit_cache[[cl_key]] %||%
                   do.call(mle_fit_ellipsoid, c(list(X = X, w = w), center_args))
            if (!is.null(f)) {
              f$axes <- cap_axes_for_cluster(f$axes, X)
              if (!all(is.finite(f$axes)) || any(f$axes <= 0)) next
              pr <- mle_project_ellipsoid(X, f$center, f$axes, f$R)
              med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
              med$distance_to_fit <- pr$dist
              fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
            }

          } else if (input$fit_type == "helix") {
            f <- mle_fit_cache[[cl_key]] %||%
                   do.call(mle_fit_helix, c(list(X = X, w = w, turns_hint = input$helix_turns), center_args))
            if (!is.null(f)) {
              t_mle <- .project_t_helix(X, f$center, f$axis, f$radius, f$pitch, f$phi0,
                                        t0 = f$t, check_neighbors = TRUE)
              H_mle <- t(vapply(t_mle, .helix_point, numeric(3),
                                c = f$center, u = f$axis, r = f$radius, p = f$pitch, phi0 = f$phi0))
              med$nearest_x <- H_mle[,1]; med$nearest_y <- H_mle[,2]; med$nearest_z <- H_mle[,3]
              med$distance_to_fit <- sqrt(rowSums((X - H_mle)^2))
              fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
            }

          } else if (input$fit_type == "ring") {
            f <- mle_fit_cache[[cl_key]] %||%
                   do.call(mle_fit_ring, c(list(X = X, w = w), center_args))
            if (!is.null(f)) {
              pr <- mle_project_ring(X, f$center, f$normal, f$radius)
              med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
              med$distance_to_fit <- pr$dist
              fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
            }

          } else if (input$fit_type == "arc") {
            f <- mle_fit_cache[[cl_key]] %||%
                   do.call(mle_fit_arc, c(list(X = X, w = w, arc_span_hint = input$arc_span_deg %||% 0), center_args))
            if (!is.null(f)) {
              pr <- mle_project_arc(X, f$center, f$normal, f$radius, f$t_center, f$arc_half)
              med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
              med$distance_to_fit <- pr$dist
              fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
            }
          }
        }
      }
    }

    if (input$dark_mode) {
      fig <- fig %>%
        layout(
          paper_bgcolor = 'black',
          plot_bgcolor = 'black',
          scene = list(
            bgcolor = "black",
            xaxis = list(
              color = "white", gridcolor = "white", zerolinecolor = "white",
              titlefont = list(color = "white"), tickfont = list(color = "white")
            ),
            yaxis = list(
              color = "white", gridcolor = "white", zerolinecolor = "white",
              titlefont = list(color = "white"), tickfont = list(color = "white")
            ),
            zaxis = list(
              color = "white", gridcolor = "white", zerolinecolor = "white",
              titlefont = list(color = "white"), tickfont = list(color = "white")
            )
          )
        )
    }
    grid_color <- if (input$hide_grid) "rgba(0,0,0,0)" else if (input$dark_mode) "white" else "black"
    
    fig <- fig %>%
      layout(
        paper_bgcolor = if (input$dark_mode) "black" else "white",
        plot_bgcolor  = if (input$dark_mode) "black" else "white",
        scene = list(
          bgcolor = if (input$dark_mode) "black" else "white",
          xaxis = list(
            showgrid = !input$hide_grid,
            zeroline = !input$hide_grid,
            showticklabels = !input$hide_grid,
            title = if (input$hide_grid) "" else "X",
            color = if (input$dark_mode) "white" else "black"
          ),
          yaxis = list(
            showgrid = !input$hide_grid,
            zeroline = !input$hide_grid,
            showticklabels = !input$hide_grid,
            title = if (input$hide_grid) "" else "Y",
            color = if (input$dark_mode) "white" else "black"
          ),
          zaxis = list(
            showgrid = !input$hide_grid,
            zeroline = !input$hide_grid,
            showticklabels = !input$hide_grid,
            title = if (input$hide_grid) "" else "Z",
            color = if (input$dark_mode) "white" else "black"
          )
        )
      )
    
    if (isTRUE(input$show_cluster2_shape) && isTRUE(input$enable_cluster2) && !is.null(df$cluster2)) {
      # One wrap per Median Cluster PER CHANNEL, the same rule the isosurface
      # follows, so an alpha shape can never join two colours (wrap_groups()).
      a_groups   <- wrap_groups(df, "cluster2")
      a_multi_ch <- length(unique(stats::na.omit(a_groups$channel))) > 1L
      for (gi in seq_len(nrow(a_groups))) {
        g   <- a_groups[gi, , drop = FALSE]
        cl2 <- g$id
        group_df <- wrap_group_rows(df, "cluster2", g)
        if (nrow(group_df) >= 4) {
          # Radius derived per cluster from its own point spacing - see
          # alpha_surface_3d() in helpers_smlm.R.
          srf <- alpha_surface_3d(as.matrix(group_df[, c("x", "y", "z")]))
          if (!is.null(srf)) {

            # pick a single solid color: first non-NA cluster2 color (fallback to light gray)
            col_hex <- group_df$cluster2_color[which(!is.na(group_df$cluster2_color) & nzchar(group_df$cluster2_color))][1]
            if (is.na(col_hex) || !nzchar(col_hex)) col_hex <- "#bdbdbd"
            # With several channels, the shared "Color wraps by channel" option
            # takes the channel's colour instead, as the isosurface does.
            if (a_multi_ch && isTRUE(input$cluster2_iso_by_channel) &&
                !is.na(g$color) && nzchar(g$color)) col_hex <- g$color

            fig <- fig %>% plotly::add_trace(
              type = "mesh3d",
              x = srf$x, y = srf$y, z = srf$z,   # srf may be subsampled; i/j/k index these
              i = srf$i, j = srf$j, k = srf$k,
              opacity = input$hull_opacity,
              # Colour as a flat intensity + one-colour scale, the encoding
              # iso_mesh_trace() uses. `color = I(hex)` on a mesh3d is DROPPED by
              # plotly_build() (with or without I()), so this wrap used to render
              # in plotly's default mesh colour whatever col_hex said.
              intensity = rep(0, length(srf$x)),
              colorscale = list(list(0, col_hex), list(1, col_hex)),
              showscale = FALSE,
              name = if (a_multi_ch && nzchar(g$label %||% ""))
                       paste0("Median Cluster ", cl2, " (", g$label, ")")
                     else paste("Median Cluster", cl2),
              hoverinfo = "skip",
              hovertemplate = "",
              showlegend = FALSE,
              inherit = FALSE
            )
          }
        }
      }
    }

    # ---- Estimated ground-truth candidates (hollow spheres, green->red) ----
    ge_3d <- gt_estimate()
    if (!is.null(ge_3d) && nrow(ge_3d$points)) {
      gpts <- ge_3d$points
      lo3  <- ge_3d$meta$link_nm %||% 0
      hi3  <- max(ge_3d$meta$footprint_nm %||% 0, 2 * lo3,
                  stats::quantile(gpts$pos_sd_nm, 0.9, na.rm = TRUE))
      col3 <- gt_uncertainty_colors(gpts$pos_sd_nm, lo3, hi3)
      sz3  <- suppressWarnings(as.numeric(input$gt_circle_size))
      sz3  <- if (!is.finite(sz3) || sz3 <= 0) 7 else max(3, sz3 * 0.5)
      struct3 <- ifelse(gpts$cluster2 > 0, paste0("Structure ", gpts$cluster2), "ungrouped")
      txt3 <- sprintf("<b>%s</b><br>GT candidate · ± %.0f nm%s<br>%d contributing unit(s) · %d TID / %d loc",
                      struct3, gpts$pos_sd_nm, ifelse(gpts$fused, " (fused: multiple GTs)", ""),
                      gpts$n_units, gpts$n_tid, gpts$n_loc)

      # Contribution fan-out lines in 3D.
      if (isTRUE(input$gt_show_lines %||% TRUE) && !is.null(ge_3d$links) && nrow(ge_3d$links)) {
        lk <- ge_3d$links
        lx <- as.vector(rbind(lk$ux, lk$gx, NA_real_))
        ly <- as.vector(rbind(lk$uy, lk$gy, NA_real_))
        lz <- as.vector(rbind(lk$uz, lk$gz, NA_real_))
        line_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.35)" else "rgba(0,0,0,0.35)"
        fig <- fig %>% plotly::add_trace(
          x = lx, y = ly, z = lz, type = "scatter3d", mode = "lines",
          line = list(color = line_col, width = 2), name = "GT contributions",
          showlegend = TRUE, inherit = FALSE, hoverinfo = "skip"
        )
      }
      fig <- fig %>% plotly::add_trace(
        x = gpts$x, y = gpts$y, z = gpts$z,
        type = "scatter3d", mode = "markers",
        marker = list(symbol = "circle-open", size = sz3,
                      color = col3, line = list(color = col3, width = 2)),
        name = "Estimated GT", showlegend = TRUE, inherit = FALSE,
        text = txt3, hoverinfo = "text"
      )
    }

    fig %>% plotly::event_register('plotly_relayout')
  })

  ## ── Dataset Groups (shared across Analysis-tab analyses) ──────────────────
  # Keep the 6 group file-selectors populated with the loaded file names.
  observeEvent(smlm_file_store(), {
    fs <- names(smlm_file_store())
    for (i in 1:6) {
      sel <- isolate(input[[paste0("grp_files_", i)]])
      updateSelectInput(session, paste0("grp_files_", i),
                        choices = fs, selected = intersect(sel, fs))
    }
  }, ignoreNULL = FALSE)

  ## GOAL: read a numericInput as a whole number inside a range, and survive the
  ## user emptying the box.
  ##
  ## WHY `%||%` IS NOT ENOUGH. An emptied numericInput does not send NULL, it sends
  ## **NA**, and `%||%` only substitutes for NULL. NA then survives as.integer(),
  ## and min()/max() propagate it rather than clamping it, so
  ## `max(1L, min(6L, NA))` is NA and the loop below died on `seq_len(NA)` the
  ## moment the group-count box was cleared. Anything that SIZES something from a
  ## user-typed number needs this, not `%||%`.
  .int_input <- function(x, default, lo = -Inf, hi = Inf) {
    v <- suppressWarnings(as.integer(x))
    if (!length(v) || is.na(v[1])) v <- as.integer(default)
    as.integer(max(lo, min(hi, v[1])))
  }

  # file -> group-name map (only the first n_groups slots, only assigned files)
  file_group_map <- reactive({
    ng <- .int_input(input$n_groups, 6L, 1L, 6L)
    m <- character(0)
    for (i in seq_len(ng)) {
      files <- input[[paste0("grp_files_", i)]]
      if (!length(files)) next
      nm <- input[[paste0("grp_name_", i)]]
      if (is.null(nm) || !nzchar(trimws(nm))) nm <- paste("Group", i)
      m[files] <- trimws(nm)
    }
    m
  })

  # group-name -> colour (from the per-group colour pickers)
  .GROUP_DEFAULT_COLORS <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf")
  group_color_map <- reactive({
    ng <- .int_input(input$n_groups, 6L, 1L, 6L)
    cols <- character(0)
    for (i in seq_len(ng)) {
      nm <- input[[paste0("grp_name_", i)]]
      if (is.null(nm) || !nzchar(trimws(nm))) nm <- paste("Group", i)
      cols[trimws(nm)] <- input[[paste0("grp_color_", i)]] %||% .GROUP_DEFAULT_COLORS[i]
    }
    cols
  })

  # Colour palette for a set of series labels. In group mode, use the group
  # colour pickers (so every analysis is consistent); otherwise the auto palette.
  series_palette <- function(labels) {
    labels <- unique(as.character(labels))
    if ((input$analysis_by %||% "file") == "group") {
      gm  <- group_color_map()
      out <- unname(gm[labels])
      miss <- is.na(out)
      if (any(miss)) {
        fb <- if (sum(miss) <= 8) RColorBrewer::brewer.pal(max(3L, sum(miss)), "Set2")[seq_len(sum(miss))]
              else viridisLite::viridis(sum(miss))
        out[miss] <- fb
      }
      return(stats::setNames(as.character(out), labels))
    }
    # File mode: CATEGORICAL colours. Files are nominal, so a sequential ramp
    # (the old Set2-then-viridis) was the wrong encoding - it made every file a
    # shade of blue-green with neighbours hard to tell apart.
    stats::setNames(qual_palette(length(labels)), labels)
  }

  # Per-FILE colours that respect the grouping: in group mode each file becomes a
  # shade of its group's colour (so the hue reads as the condition and the shade as
  # the replicate); in file mode each file gets its own categorical colour.
  # Returns a named vector indexed by FILE name.
  file_palette <- function(files) {
    files <- unique(as.character(files))
    if (!length(files)) return(character(0))
    if ((input$analysis_by %||% "file") != "group")
      return(stats::setNames(qual_palette(length(files)), files))
    gm <- file_group_map()
    gc <- group_color_map()
    out <- stats::setNames(rep(NA_character_, length(files)), files)
    for (g in unique(unname(gm[intersect(files, names(gm))]))) {
      fs   <- intersect(files, names(gm)[gm == g])
      base <- unname(gc[g])
      if (is.na(base)) base <- QUAL_COLORS[1]
      out[fs] <- shade_colors(base, length(fs))
    }
    miss <- is.na(out)
    if (any(miss)) out[miss] <- .adjust_lum(qual_palette(sum(miss)), 0.25)
    out
  }

  # ── CHANNEL FACET for univariate analyses ─────────────────────────────────
  #
  # Splits an already-prepared per-file frame into one subset per channel, so a
  # distribution can be read per channel instead of pooled. Returns a list of
  # (label, df); with the facet off, or on single-channel data, that list is one
  # entry with an empty label and the caller behaves exactly as before.
  #
  # WHERE THIS IS SAFE, AND WHERE IT IS NOT.
  #   Safe for DISTRIBUTIONS that do not divide by an area: NN distances, radius
  #   of gyration, molecules per Median Cluster. Subsetting changes what is in
  #   the sample but not what the statistic means.
  #
  #   NOT safe, as written, for Ripley's L or g(r). Those divide by the observed
  #   area, and the observation window is currently derived from the points that
  #   are passed in - so a channel with fewer localizations would get a smaller
  #   footprint, a higher apparent density, and a DEFLATED curve. The window has
  #   to be built once from the pooled file and then shared by both channels,
  #   which is a change to spatial_gate_close()'s plumbing rather than a call to
  #   this function. Until that exists the spatial statistics stay pooled: a
  #   pooled curve is honest, a per-channel curve with a per-channel window is
  #   not.
  channel_facets <- function(df) {
    none <- list(list(label = "", df = df))
    if (!isTRUE(input$analysis_by_channel %||% FALSE)) return(none)
    if (!is.data.frame(df) || !nrow(df) || !"channel" %in% names(df)) return(none)
    ids <- sort(unique(df$channel[!is.na(df$channel)]))
    if (length(ids) < 2L) return(none)
    lapply(ids, function(ch) {
      sub <- df[which(df$channel == ch), , drop = FALSE]
      list(label = as.character(sub$channel_label[1]), df = sub)
    })
  }

  # Append the channel label to a dataset name, so one series per file becomes one
  # series per file-and-channel and every existing grouping/colour path downstream
  # keeps working unchanged.
  .facet_name <- function(fname, label)
    if (nzchar(label %||% "")) paste0(fname, " · ", label) else fname

  # Relabel a results data frame's `dataset` column to its group. In group mode,
  # files not assigned to any group are dropped. Analyses call this at read time,
  # so toggling files<->groups updates plots without re-running.
  apply_dataset_groups <- function(res, col = "dataset") {
    if (is.null(res) || !nrow(res) || !col %in% names(res)) return(res)
    if ((input$analysis_by %||% "file") != "group") return(res)
    m <- file_group_map()
    if (!length(m)) return(res[0, , drop = FALSE])
    res <- res[res[[col]] %in% names(m), , drop = FALSE]
    if (nrow(res)) res[[col]] <- unname(m[res[[col]]])
    res
  }

  ## ── Per-file summary export (unit of replication = the file / cell) ────────
  ## GOAL: one CSV with per-file summaries of every distribution analysis that has
  ## been run, tagged with each file's group, so between-group statistics use one
  ## value per cell. Group-mode histograms pool units across files, which is fine
  ## for looking but pseudo-replicates any test run on them; see per_file_summary()
  ## in helpers_smlm.R. Reads the stored results at download time, so it reflects
  ## the last Run of each panel. The result stores are defined further down this
  ## function; they exist by the time the handler can fire.
  output$download_per_file_summary <- downloadHandler(
    filename = function() paste0("per_file_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      rv <- function(f) tryCatch(f(), error = function(e) NULL)   # a store may be unset
      cross <- rv(cross_nn_results)
      if (is.data.frame(cross) && nrow(cross) && "direction" %in% names(cross))
        cross$dataset <- paste(cross$dataset, cross$direction, sep = " · ")  # one row per direction
      results <- list(
        batch_nn         = list(df = rv(nn_batch_results),  metrics = "distance",   unit = "NN distance (nm)"),
        units_per_median_cluster = list(df = rv(tid_count_results), metrics = "n_tids", unit = "Median Cluster"),
        radius_of_gyration = list(df = rv(rog_results),     metrics = c("rog", "n_units"), unit = "cluster"),
        crowding         = list(df = rv(crowd_results),     metrics = c("nn_median", "time_shift", "n_tids"), unit = "Median Cluster"),
        within_cluster_nn_vs_size = list(df = rv(size_nn_results), metrics = c("nn_median", "size"), unit = "Median Cluster"),
        cross_channel_nn = list(df = cross,                  metrics = "distance",   unit = "cross-channel NN distance (nm)"))
      out <- per_file_summary(results, group_map = file_group_map(),
                              files = names(smlm_file_store() %||% list()))
      if (!nrow(out))
        out <- data.frame(note = "No distribution analysis has been run yet. Run a panel, then download again.")
      utils::write.csv(out, file, row.names = FALSE)
    })

  ## ══════════════════════════════════════════════════════════════════════════
  ## LABELLING-DENSITY NORMALISATION  (opt-in, default OFF)
  ## ══════════════════════════════════════════════════════════════════════════
  ## GOAL: make clustering read-outs comparable between acquisitions whose
  ## LABELLING DENSITY differs for reasons that are not biology - transient
  ## transfection giving each cell an undefined plasmid load, different expression
  ## level, run-to-run labelling efficiency.
  ##
  ## Two independent opt-ins, both off by default because most experiments do not
  ## need them (a stable line, or one acquisition, has no such confound):
  ##   1. `density_profile_on`  - just MEASURE and report per-file density.
  ##   2. `density_thin_on`     - additionally THIN every file to one common
  ##                              density before any analysis sees it.
  ##
  ## The maths, the choice of thinning target and the trap in using NN distances
  ## as that target are all documented at the "8b. LABELLING-DENSITY
  ## NORMALISATION" banner in helpers_smlm.R. Read that before changing anything
  ## here. In one line: g(r) and K(r) are exactly invariant under independent
  ## thinning, so thinning removes a density artefact without removing real
  ## clustering - but the TARGET must be the global density n/support, which
  ## clustering cannot inflate, not the NN-based local intensity, which it can
  ## (measured 4.25x on a clustered pattern).
  density_profile_rv <- reactiveVal(NULL)   # one row per file (see compute_unit_density)
  density_sweep_rv   <- reactiveVal(NULL)   # support vs closing radius, per file
  density_local_rv   <- reactiveVal(NULL)   # long df(dataset, local) - per-UNIT density

  # Which unit the density is measured and thinned on. De-blinked molecules
  # whenever de-blinking is on: a trace ID is ONE BLINK, so trace-ID density
  # mixes expression level with photophysics and laser power and is not
  # comparable between runs.
  density_unit <- reactive({
    u <- input$dens_norm_unit %||% "auto"
    if (u != "auto") return(u)
    if (isTRUE(input$enable_clustering)) "cluster" else "tid"
  })

  # Master gates. Thinning implies profiling (it needs the per-file densities).
  density_profile_on <- reactive(isTRUE(input$dens_norm_enable))
  density_thin_on    <- reactive(isTRUE(input$dens_norm_enable) && isTRUE(input$dens_thin_enable))

  ## ── Measure: build the per-file density profile ───────────────────────────
  observeEvent(input$run_density, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }

    unit    <- density_unit()
    use_3d  <- isTRUE(input$dens_norm_3d)
    smode   <- input$dens_support_mode %||% "grid"
    h_nm    <- max(10,  as.numeric(input$dens_grid_h     %||% 100))
    close_n <- max(0,   as.numeric(input$dens_close_nm   %||% 200))
    kk      <- max(1L,  as.integer(input$dens_nn_k       %||% 3L))
    kde_bw  <- max(10,  as.numeric(input$dens_kde_bw     %||% 150))
    kde_fr  <- min(1, max(0.5, as.numeric(input$dens_kde_frac %||% 0.99)))
    needs_cluster <- unit == "cluster"

    auto_close <- isTRUE(input$dens_close_auto) && smode == "grid"
    want_sweep <- auto_close || (isTRUE(input$dens_show_sweep) && smode == "grid")

    rows <- list(); skipped <- character(0); sweeps <- list()
    prepped <- list(); locals <- list()
    withProgress(message = "Measuring labelling density...", value = 0, {
      nf <- length(store)

      ## ── Pass 1: filter (and de-blink) each file, and sweep the closing radius ──
      # The sweep must finish before any density is computed, because the radius it
      # picks has to be THE SAME for every file - different denominators measuring
      # different kinds of region is precisely what makes a density ratio wrong.
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(0.6 / nf, detail = fname)
        df <- qc_exclude_from(store[[fname]], fname)
        if (!is.data.frame(df) || !nrow(df)) { skipped <- c(skipped, fname); next }
        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) { skipped <- c(skipped, fname); next }
        # De-blinking must run before molecule density can be measured. NOT thinned
        # here - this is the pre-thinning reference the whole panel is built on.
        if (needs_cluster) {
          df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) df)
        }
        prepped[[fname]] <- df
        # Sweep on the same points the support is built from (all localizations),
        # so the plateau the user reads is the plateau the density actually used.
        if (want_sweep) {
          keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
          if (all(keep %in% names(df))) {
            sp <- df[, keep, drop = FALSE]
            if (!use_3d) sp$z <- 0
            sw <- tryCatch(support_closing_sweep(sp[stats::complete.cases(sp), , drop = FALSE],
                                                 use_3d = use_3d, h_nm = h_nm),
                           error = function(e) NULL)
            if (!is.null(sw)) { sw$dataset <- fname; sweeps[[fname]] <- sw }
          }
        }
      }

      ## ── Resolve the one common gap-fill radius ──
      # Each file's plateau start is read off its own sweep, then the LARGEST is
      # taken: a radius below a file's plateau leaves that file's support tracking
      # its structures, and it is better to be slightly generous on the files that
      # plateaued early than wrong on the one that did not.
      if (auto_close && length(sweeps)) {
        pls <- vapply(sweeps, function(s) plateau_close_nm(s), numeric(1))
        pls <- pls[is.finite(pls)]
        if (length(pls)) close_n <- max(pls)
      }

      ## ── Pass 2: densities, all on that one radius ──
      for (fname in names(prepped)) {
        incProgress(0.4 / max(1L, length(prepped)), detail = fname)
        df <- prepped[[fname]]
        r <- tryCatch(compute_unit_density(df, unit = unit, use_3d = use_3d,
                                          support_mode = smode, h_nm = h_nm,
                                          close_nm = close_n, k = kk,
                                          support_from = input$dens_support_from %||% "units",
                                          bw_nm = kde_bw, frac = kde_fr),
                      error = function(e) NULL)
        if (is.null(r)) { skipped <- c(skipped, fname); next }
        r$dataset    <- fname
        r$unit       <- unit
        r$close_auto <- auto_close
        # Freeze the settings this row was measured with, so the status line and
        # the density plot describe the run rather than the live sliders.
        r$run_mode   <- smode
        r$run_bw     <- kde_bw
        r$run_frac   <- kde_fr
        r$run_3d     <- use_3d
        rows[[fname]] <- r
        # Per-unit local density for the distribution plot. Subsampled - the plot
        # only needs the shape, and violins on 10^5 points are slow to render.
        lv <- tryCatch(unit_local_density(unit_centroids(df, unit, use_3d),
                                         k = kk, use_3d = use_3d),
                       error = function(e) numeric(0))
        lv <- lv[is.finite(lv) & lv > 0]
        if (length(lv) > 4000L) lv <- lv[seq(1, length(lv), length.out = 4000L)]
        if (length(lv)) locals[[fname]] <- data.frame(dataset = fname, local = lv)
      }
    })
    res <- do.call(rbind, rows)
    if (is.null(res) || !nrow(res)) {
      showNotification("Could not measure density on any file (need >=3 units per file).",
                       type = "warning")
      density_profile_rv(NULL); return()
    }
    density_profile_rv(res)
    density_sweep_rv(if (length(sweeps)) do.call(rbind, sweeps) else NULL)
    density_local_rv(if (length(locals)) do.call(rbind, locals) else NULL)
    rr <- range(res$density, na.rm = TRUE)
    showNotification(sprintf(
      "Density measured on %d file(s): %.1f - %.1f units/%s (%.1fx spread).",
      nrow(res), rr[1], rr[2], if (use_3d) "um3" else "um2",
      if (rr[1] > 0) rr[2] / rr[1] else NA_real_),
      type = "message", duration = 5)
  })

  ## ── Target: the one density every file is thinned to ──────────────────────
  # Returned in the same units as density_profile_rv()$density. The target must be
  # REACHABLE by thinning alone: nothing can be up-sampled, so a file already below
  # the target is left untouched and flagged rather than silently misrepresented.
  density_target <- reactive({
    prof <- density_profile_rv()
    if (is.null(prof) || !nrow(prof)) return(NULL)
    mode <- input$dens_target_mode %||% "min"
    if (mode == "manual") {
      v <- suppressWarnings(as.numeric(input$dens_target_manual %||% NA))
      return(if (is.finite(v) && v > 0) v else NULL)
    }
    if (mode == "control" && (input$analysis_by %||% "file") == "group") {
      gm <- file_group_map()
      ctl <- input$dens_target_group %||% ""
      if (nzchar(ctl) && length(gm)) {
        fs <- names(gm)[gm == ctl]
        v  <- stats::median(prof$density[prof$dataset %in% fs], na.rm = TRUE)
        if (is.finite(v) && v > 0) return(v)
      }
      return(NULL)
    }
    # "min" (default): the lowest density in the comparison set. Requires no
    # up-sampling of any file, so it is the only target that is always reachable.
    scope <- prof
    if ((input$analysis_by %||% "file") == "group") {
      gm <- file_group_map()
      if (length(gm)) scope <- prof[prof$dataset %in% names(gm), , drop = FALSE]
    }
    if (!nrow(scope)) return(NULL)
    v <- min(scope$density, na.rm = TRUE)
    if (is.finite(v) && v > 0) v else NULL
  })

  # Per-file keep-probability p = target / observed, capped at 1.
  density_p_map <- reactive({
    prof <- density_profile_rv(); tgt <- density_target()
    if (is.null(prof) || is.null(tgt)) return(NULL)
    p <- pmin(1, tgt / prof$density)
    stats::setNames(p, prof$dataset)
  })

  ## ── The gate: thin one file's localizations ───────────────────────────────
  # Mirrors qc_exclude_from's contract - call it on a df together with the file it
  # came from, and it returns the df unchanged whenever the feature is off. Must be
  # applied AFTER the metadata filters and BEFORE clustering, because making the
  # fixed-(eps, minPts) DBSCAN comparable is the entire point.
  #
  # The seed is per-file and fixed, so a given file always thins to the same subset
  # (results are reproducible across re-runs), and thin_units_to_p saves/restores
  # the global RNG stream so no other randomised analysis is perturbed.
  apply_density_thinning <- function(df, fname) {
    if (!density_thin_on() || is.null(df) || !nrow(df)) return(df)
    pm <- density_p_map()
    if (is.null(pm) || !fname %in% names(pm)) return(df)
    p <- pm[[fname]]
    if (!is.finite(p) || p >= 1) return(df)
    unit <- density_unit()
    # Molecule-level thinning needs the de-blinked IDs to exist. De-blinking is a
    # deterministic DBSCAN, so deriving it here and again downstream on the
    # surviving points gives the same molecule assignment.
    if (unit == "cluster" && !"cluster" %in% names(df)) {
      df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) df)
    }
    seed <- as.integer(sum(utf8ToInt(fname)) %% 100000L) + 1L
    thin_units_to_p(df, unit = unit, p = p, seed = seed)
  }

  ## ── Readouts ──────────────────────────────────────────────────────────────
  output$dens_norm_active <- reactive(density_thin_on())
  outputOptions(output, "dens_norm_active", suspendWhenHidden = FALSE)

  output$density_status <- renderUI({
    if (!density_profile_on())
      return(helpText("Off - analyses run on the data as acquired."))
    prof <- density_profile_rv()
    if (is.null(prof) || !nrow(prof))
      return(helpText("Press 'Measure density' to build the per-file profile."))
    un  <- if (isTRUE(input$dens_norm_3d)) "units/µm³" else "units/µm²"
    ul  <- if (prof$unit[1] == "cluster") "de-blinked molecules" else "trace IDs"
    rr  <- range(prof$density, na.rm = TRUE)
    tgt <- density_target()
    spread <- if (rr[1] > 0) rr[2] / rr[1] else NA_real_

    # Headline: the spread IS the confound. Everything else is a flag, shown only
    # when it fires, so the panel stays quiet when the measurement is clean.
    lines <- list(tags$div(style = "font-weight:600;", sprintf(
      "%d file(s) · %.1f–%.1f %s · %.1f× spread · unit = %s",
      nrow(prof), rr[1], rr[2], un, spread, ul)))
    # From the PROFILE, so the description matches the measurement.
    smode_now <- (if ("run_mode" %in% names(prof)) prof$run_mode[1] else NULL) %||%
                 (input$dens_support_mode %||% "kde")
    if (smode_now == "grid")
      lines <- c(lines, list(tags$div(style = "color:#777;", sprintf(
        "Support: %g nm grid, gap fill %g nm (%s) - dark regions wider than ~%g nm stay excluded.",
        prof$h_nm[1], prof$close_nm[1],
        if (isTRUE(prof$close_auto[1])) "auto" else "manual", 2 * prof$close_nm[1]))))
    if (smode_now == "kde")
      lines <- c(lines, list(tags$div(style = "color:#777;", sprintf(
        "Support: smoothed density, region enclosing %.0f%% of units, %g nm bandwidth on a %g nm grid - dark regions wider than ~%g nm fall below the level.",
        100 * ((if ("run_frac" %in% names(prof)) prof$run_frac[1] else NULL) %||% 0.99),
        (if ("run_bw" %in% names(prof)) prof$run_bw[1] else NULL) %||% 150,
        prof$h_nm[1],
        6 * ((if ("run_bw" %in% names(prof)) prof$run_bw[1] else NULL) %||% 150)))))

    warn <- function(txt) tags$div(style = "color:#b26a00;", txt)
    if (is.finite(spread) && spread < 1.5)
      lines <- c(lines, list(tags$div(style = "color:#1a8a1a;",
        "Already comparable (< 1.5×) - thinning would only cost points.")))

    # THE failure mode of the support estimate: with too small a gap-fill radius the
    # occupancy grid outlines the STRUCTURES instead of the cell, so the density
    # denominator absorbs the clustering it is meant to be independent of - and by a
    # different amount per file. Measured: 200 nm gave a 3.13 ratio where truth was
    # 4.0; the plateau gave 3.82.
    if (smode_now %in% c("grid", "kde")) {
      fb <- prof$frac_of_bbox
      if (any(is.finite(fb) & fb < 0.5))
        lines <- c(lines, list(warn(sprintf(
          "Support is only %.0f–%.0f%% of the frame. If your sample fills more than that, the smoothing is too tight and the support is outlining your structures, not the cell - check the support map.",
          100 * min(fb, na.rm = TRUE), 100 * max(fb, na.rm = TRUE)))))
      if (any(is.finite(fb) & fb > 0.97))
        lines <- c(lines, list(warn(
          "Support has reached the bounding box - dark regions are being counted as sample. Tighten the smoothing (lower bandwidth, or a lower enclosed fraction).")))
      if (sum(is.finite(fb)) >= 2 && diff(range(fb, na.rm = TRUE)) > 0.25)
        lines <- c(lines, list(warn(
          "Support fraction differs by > 25 points between files, so the denominators are not the same kind of region.")))
    }
    if (density_thin_on()) {
      pm <- density_p_map()
      if (is.null(pm) || is.null(tgt)) {
        lines <- c(lines, list(warn("Thinning on but no target resolved - check the target setting.")))
      } else {
        worst <- if (length(pm)) min(pm, na.rm = TRUE) else NA_real_
        lines <- c(lines, list(tags$div(style = "font-weight:600; color:#c03030;", sprintf(
          "Thinning ON → %.1f %s · %d/%d thinned · densest keeps %.0f%%",
          tgt, un, sum(pm < 1), length(pm), 100 * worst))))
        if (is.finite(worst) && worst < 0.2)
          lines <- c(lines, list(warn(sprintf(
            "Only %.0f%% of the densest file kept. g(r) survives this, but count-based read-outs lose precision and small clusters can fall below minPts.",
            100 * worst))))
      }
    }
    do.call(tags$div, c(list(style = "font-size:0.85em; line-height:1.5;"), lines))
  })

  output$density_table <- DT::renderDataTable({
    prof <- density_profile_rv()
    shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0,
                                "Press 'Measure density'."))
    pm  <- density_p_map()
    d3  <- isTRUE(input$dens_norm_3d)
    out <- data.frame(
      File           = prof$dataset,
      Units          = prof$n_units,
      `Support (um2/3)` = signif(prof$support_nm / 1000^prof$dims, 4),
      `Support/bbox` = round(prof$frac_of_bbox, 3),
      Density        = round(prof$density, 2),
      `Density (NN)` = round(prof$density_nn, 2),
      `Median NN (nm)` = round(prof$nn_median, 1),
      `Keep p`       = if (!is.null(pm)) round(unname(pm[prof$dataset]), 3) else NA_real_,
      check.names    = FALSE, stringsAsFactors = FALSE)
    DT::datatable(out, rownames = FALSE, options = list(dom = "t", paging = FALSE,
                  scrollX = TRUE), class = "compact stripe") %>%
      DT::formatStyle("Keep p", color = DT::styleInterval(c(0.2), c("#b26a00", "#333")))
  })

  ## Per-file density DISTRIBUTION: one horizontal violin per dataset over the
  ## per-unit local density, with that file's global n/support marked on it.
  ##
  ## Why a distribution and not one dot per file: the global number says how much
  ## label there is, the shape says how it is arranged. A cell with real aggregates
  ## is heterogeneous (sparse background mode + high-density tail); a uniformly
  ## labelled cell at the same mean is narrow and unimodal. So the width here is
  ## itself clustering information - and it is also the visual reason the mean of
  ## this distribution must not be the thinning target (see nn_local_intensity).
  ##
  ## Series are identified by the y-axis tick labels, so no legend is needed - which
  ## is what keeps the panel readable with several files.
  ##
  ## AXIS: the x range is set EXPLICITLY from a robust quantile span of the data
  ## (plus the global-density markers), never left to plotly's autorange - a violin's
  ## KDE tails run well past the data and, on a log axis, autorange then pads out to
  ## whole decades, which is what left the panel long and mostly empty. Log is
  ## offered because local density is a ratio-scale quantity that can span decades
  ## inside one cell, and density comparisons are multiplicative (a 4x difference
  ## should look the same anywhere on the axis) - but it is a toggle, since a sample
  ## with a narrow spread reads better linear.
  build_density_plot_fig <- function() {
    prof <- density_profile_rv()
    shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0, "Press 'Measure density'."))
    loc  <- density_local_rv()
    grouped <- (input$analysis_by %||% "file") == "group"
    # dims from the profile: the 2D/3D switch is a compute parameter.
    un   <- if (isTRUE((prof$dims[1] %||% 2L) == 3L)) "units/µm³" else "units/µm²"

    pr <- apply_dataset_groups(prof)
    shiny::validate(shiny::need(nrow(pr) > 0, "No files assigned to groups (group mode)."))
    lc <- if (!is.null(loc)) apply_dataset_groups(loc) else NULL

    # Order rows by the group/file median so the panel reads bottom-to-top.
    ord <- stats::aggregate(density ~ dataset, data = pr, FUN = stats::median)
    ord <- ord$dataset[order(ord$density)]
    pal <- series_palette(ord)

    use_log <- isTRUE(input$dens_plot_log)
    ax      <- density_axis_range(
      if (!is.null(lc)) lc$local else numeric(0), pr$density, use_log = use_log)
    xrange  <- ax$range
    clip    <- ax$clip

    fig <- plotly::plot_ly()
    for (g in ord) {
      col <- unname(pal[g])
      if (!is.null(lc) && any(lc$dataset == g)) {
        fig <- fig %>% plotly::add_trace(
          x = clip(lc$local[lc$dataset == g]), y = rep(g, sum(lc$dataset == g)),
          type = "violin", orientation = "h", name = g, legendgroup = g,
          showlegend = FALSE, side = "positive", points = FALSE,
          line = list(color = col, width = 1),
          fillcolor = sprintf("rgba(%s,0.35)", paste(grDevices::col2rgb(col), collapse = ",")),
          hoverinfo = "x", meanline = list(visible = FALSE), box = list(visible = FALSE))
      }
      sub <- pr[pr$dataset == g, , drop = FALSE]
      fig <- fig %>% plotly::add_trace(
        x = sub$density, y = rep(g, nrow(sub)), type = "scatter", mode = "markers",
        name = g, legendgroup = g, showlegend = FALSE,
        marker = list(color = col, size = 11, symbol = "diamond",
                      line = list(color = "#222", width = 1)),
        text = sprintf("%s<br><b>%.1f %s</b> (global n/support)<br>%d units · support/bbox %.2f",
                       sub$dataset, sub$density, un, sub$n_units, sub$frac_of_bbox),
        hoverinfo = "text")
    }

    tgt <- density_target()
    shp <- list(); ann <- list()
    if (!is.null(tgt) && density_thin_on()) {
      shp <- list(list(type = "line", x0 = tgt, x1 = tgt, y0 = 0, y1 = 1, yref = "paper",
                       line = list(color = "#c03030", dash = "dash", width = 2)))
      # An annotation's x is in AXIS coordinates, which on a log axis means log10 of
      # the value - a shape's x0/x1 are not. Getting this wrong parks the label off
      # frame.
      ann <- list(list(x = if (use_log) log10(max(tgt, 1e-9)) else tgt,
                       y = 1, yref = "paper",
                       text = " target", showarrow = FALSE, xanchor = "left",
                       font = list(color = "#c03030", size = 11)))
    }
    fig %>% plotly::layout(
      shapes = shp, annotations = ann, showlegend = FALSE,
      xaxis = list(title = sprintf("Density (%s%s)", un, if (use_log) ", log" else ""),
                   type = if (use_log) "log" else "linear",
                   range = xrange, zeroline = FALSE),
      yaxis = list(title = "", categoryorder = "array", categoryarray = ord,
                   automargin = TRUE),
      margin = list(l = 10, r = 24, t = 28, b = 44),
      title = list(text = "◆ global density   ▏shaded = spread of per-unit local density",
                   font = list(size = 11, color = "#666"), x = 0),
      paper_bgcolor = "white", plot_bgcolor = "white")
  }
  output$density_plot <- renderPlotly(build_density_plot_fig())

  ## Support vs gap-fill radius - the diagnostic that makes the closing radius a
  ## choice rather than a guess. See support_closing_sweep for how to read it.
  build_density_sweep_plot_fig <- function() {
    sw <- density_sweep_rv()
    shiny::validate(shiny::need(!is.null(sw) && nrow(sw) > 0,
      "Tick 'Show gap-fill sweep' and press 'Measure density'."))
    pal  <- series_palette(unique(sw$dataset))
    # Mark the radius ACTUALLY used. With auto on that is the resolved plateau, not
    # the (hidden) manual input, which would otherwise draw a misleading marker.
    prof <- density_profile_rv()
    cur  <- if (!is.null(prof) && nrow(prof) && is.finite(prof$close_nm[1])) prof$close_nm[1]
            else max(0, as.numeric(input$dens_close_nm %||% 200))
    fig <- plotly::plot_ly()
    for (nm in unique(sw$dataset)) {
      s <- sw[sw$dataset == nm, ]; s <- s[order(s$close_nm), ]
      fig <- fig %>% plotly::add_trace(
        x = s$close_nm, y = s$frac_of_bbox, type = "scatter", mode = "lines+markers",
        name = nm, line = list(color = unname(pal[nm]), width = 2),
        marker = list(color = unname(pal[nm]), size = 6),
        text = sprintf("%s<br>fill %.0f nm<br>support %.3g nm^%d (%.0f%% of bbox)<br>step %+.1f%%",
                       nm, s$close_nm, s$support, if (isTRUE(input$dens_norm_3d)) 3 else 2,
                       100 * s$frac_of_bbox, 100 * s$rel_step),
        hoverinfo = "text")
    }
    fig %>% plotly::layout(
      shapes = list(list(type = "line", x0 = cur, x1 = cur, y0 = 0, y1 = 1, yref = "paper",
                         line = list(color = "#c03030", dash = "dot", width = 2))),
      annotations = list(list(x = cur, y = 1, yref = "paper",
                              text = sprintf(" in use: %g nm", cur),
                              showarrow = FALSE, xanchor = "left",
                              font = list(color = "#c03030", size = 10))),
      xaxis = list(title = "Gap fill radius (nm)"),
      yaxis = list(title = "Support / frame", range = c(0, 1.04)),
      # The curve rises left-to-right, so the top-left of the panel is always empty -
      # putting the legend there instead of below the axis buys back the vertical
      # space a horizontal legend was taking from the plot itself.
      margin = list(l = 58, r = 24, t = 26, b = 46),
      title = list(text = "Pick a radius where the curve has flattened",
                   font = list(size = 11, color = "#666"), x = 0),
      legend = list(x = 0.02, y = 0.98, xanchor = "left", yanchor = "top",
                    font = list(size = 10), bgcolor = "rgba(255,255,255,0.75)",
                    bordercolor = "#ddd", borderwidth = 1),
      paper_bgcolor = "white", plot_bgcolor = "white")
  }
  output$density_sweep_plot <- renderPlotly(build_density_sweep_plot_fig())

  ## ── Support map: what the density denominator actually counted ─────────────
  ## GOAL: draw one dataset in 2D with the estimated support outlined in GREEN and
  ## everything inside the imaged frame that was DISREGARDED shaded in RED.
  ##
  ## This is the panel that makes the support estimate falsifiable. The numbers say
  ## "support is 30% of the frame"; only the picture says whether that 30% is the
  ## cell or is the outline of your own structures - which is the difference between
  ## a density and a disguised clustering measure.
  ##
  ## Recomputed on demand for the SELECTED file only (masks for every file would be
  ## needlessly heavy), using the same filter chain and the same resolved gap-fill
  ## radius the profile used, so what is drawn is what was measured.
  observe({
    fs <- names(smlm_file_store())
    sel <- isolate(input$dens_map_file)
    updateSelectInput(session, "dens_map_file", choices = fs,
                      selected = if (!is.null(sel) && sel %in% fs) sel else
                                 (isolate(input$active_file) %||% fs[1]))
  })

  build_density_map_fig <- function() {
    prof <- density_profile_rv()
    shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0, "Press 'Measure density'."))
    smode <- input$dens_support_mode %||% "grid"
    store <- smlm_file_store()
    fname <- input$dens_map_file %||% (input$active_file %||% names(store)[1])
    shiny::validate(shiny::need(!is.null(fname) && fname %in% names(store), "Pick a dataset."))

    use_3d <- isTRUE(input$dens_norm_3d)
    unit   <- density_unit()
    h_nm   <- max(10, as.numeric(input$dens_grid_h %||% 100))
    # The radius the PROFILE resolved (auto plateau or manual), not the raw input -
    # otherwise the picture could disagree with the measured number.
    close_n <- if ("close_nm" %in% names(prof) && is.finite(prof$close_nm[1]))
      prof$close_nm[1] else max(0, as.numeric(input$dens_close_nm %||% 200))

    df <- qc_exclude_from(store[[fname]], fname)
    df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No localizations after filtering."))
    if (unit == "cluster")
      df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) df)

    keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
    shiny::validate(shiny::need(all(keep %in% names(df)), "Missing coordinates."))
    sp <- df[, keep, drop = FALSE]
    if (!use_3d) sp$z <- 0
    sp <- sp[stats::complete.cases(sp), , drop = FALSE]
    shiny::validate(shiny::need(nrow(sp) >= 3, "Too few localizations."))

    sup <- estimate_support(sp, use_3d = use_3d, mode = smode, h_nm = h_nm,
                            close_nm = close_n, return_mask = TRUE,
                            bw_nm = max(10, as.numeric(input$dens_kde_bw %||% 150)),
                            frac = min(1, max(0.5, as.numeric(input$dens_kde_frac %||% 0.99))))
    shiny::validate(shiny::need(!is.null(sup), "Could not estimate the support."))

    win  <- if (!is.null(sup$win)) sup$win else apply(as.matrix(sp[, c("x","y")]), 2, range)
    npts <- nrow(sp)
    sp_d <- if (npts > 30000L) sp[seq(1, npts, length.out = 30000L), , drop = FALSE] else sp

    # Red = disregarded, but only INSIDE the data bounding box. Two reasons: the grid's
    # padding ring exists solely to give the dilation room, and a localization table
    # carries no record of the real scan field - so the widest defensible "frame" is
    # the extent of the data itself. Consequence, stated in the title: red marks
    # interior dark regions and concave margins, never area outside the data.
    fig <- plotly::plot_ly()
    shp <- list()

    # Grid and KDE both produce a cell mask, so they share one renderer; only hull
    # and bbox need the polygon path. `sup$h` is the cell size the estimator
    # actually used, which is NOT necessarily the h_nm input.
    if (smode %in% c("grid", "kde")) {
      shiny::validate(shiny::need(!is.null(sup$mask), "Could not build the support mask."))
      m2 <- mask_project_2d(sup$mask)
      d  <- dim(m2)
      cell <- sup$h %||% h_nm
      x0 <- sup$origin[1]; y0 <- sup$origin[2]
      # Cell CENTRES for the heatmap (plotly places a heatmap cell on its centre);
      # the outline below uses cell EDGES. Mixing the two shifts the overlay by h/2.
      xc <- x0 + (seq_len(d[1]) - 0.5) * cell
      yc <- y0 + (seq_len(d[2]) - 0.5) * cell
      in_frame <- outer(xc >= win[1, 1] & xc <= win[2, 1],
                        yc >= win[1, 2] & yc <= win[2, 2], "&")
      zred <- ifelse(!m2 & in_frame, 1, NA_real_)
      fig <- fig %>% plotly::add_heatmap(
        x = xc, y = yc, z = t(zred), showscale = FALSE, hoverinfo = "skip",
        colorscale = list(list(0, "rgba(200,48,48,0.22)"), list(1, "rgba(200,48,48,0.22)")))
      outl <- mask_outline_segments(m2, c(x0, y0), cell)
    } else {
      # Hull and bbox have no interior structure, so they are drawn exactly as
      # POLYGONS rather than discretised: a red rectangle over the whole data box,
      # then the support polygon filled opaque white on top of it. Points and the
      # green outline are traces, which plotly draws above all shapes.
      rect_path <- sprintf("M %f,%f L %f,%f L %f,%f L %f,%f Z",
                           win[1,1], win[1,2], win[2,1], win[1,2],
                           win[2,1], win[2,2], win[1,1], win[2,2])
      if (smode == "hull") {
        hidx <- tryCatch(grDevices::chull(sp$x, sp$y), error = function(e) integer(0))
        shiny::validate(shiny::need(length(hidx) >= 3, "Could not build a convex hull."))
        hx <- sp$x[hidx]; hy <- sp$y[hidx]
        hull_path <- paste0("M ", paste(sprintf("%f,%f", hx, hy), collapse = " L "), " Z")
        shp <- list(
          list(type = "path", path = rect_path, fillcolor = "rgba(200,48,48,0.22)",
               line = list(width = 0), layer = "below"),
          list(type = "path", path = hull_path, fillcolor = "#ffffff",
               line = list(width = 0), layer = "below"))
        outl <- list(x = c(hx, hx[1]), y = c(hy, hy[1]))
      } else {
        # bbox: the support IS the frame, so nothing is disregarded.
        shp  <- list()
        outl <- list(x = c(win[1,1], win[2,1], win[2,1], win[1,1], win[1,1]),
                     y = c(win[1,2], win[1,2], win[2,2], win[2,2], win[1,2]))
      }
    }

    fig <- fig %>% plotly::add_trace(
      x = sp_d$x, y = sp_d$y, type = "scattergl", mode = "markers",
      marker = list(size = 2, color = "rgba(40,40,40,0.55)"),
      name = "localizations", hoverinfo = "skip", showlegend = FALSE)
    if (length(outl$x))
      fig <- fig %>% plotly::add_trace(
        x = outl$x, y = outl$y, type = "scatter", mode = "lines",
        line = list(color = "#1a8a1a", width = 1.6), connectgaps = FALSE,
        name = "counted as sample", hoverinfo = "skip", showlegend = FALSE)
    # Frame outline, so "how much of the field was disregarded" is readable.
    fig <- fig %>% plotly::add_trace(
      x = c(win[1,1], win[2,1], win[2,1], win[1,1], win[1,1]),
      y = c(win[1,2], win[1,2], win[2,2], win[2,2], win[1,2]),
      type = "scatter", mode = "lines",
      line = list(color = "#999", width = 1, dash = "dot"),
      name = "frame", hoverinfo = "skip", showlegend = FALSE)

    mode_lbl <- switch(smode,
      kde  = sprintf("smoothed density, %.0f%% of units enclosed, %g nm bandwidth",
                     100 * (sup$frac %||% 0.99), sup$bw %||% NA_real_),
      grid = sprintf("%g nm grid, %g nm gap fill (dark regions wider than ~%g nm stay excluded)",
                     h_nm, close_n, 2 * close_n),
      hull = "convex hull - bridges every concavity and dark region by construction",
      "bounding box - nothing is excluded")
    ttl <- sprintf(
      "%s - <span style='color:#1a8a1a'>green = counted as sample (%.0f%% of the data's bounding box)</span>, <span style='color:#c03030'>red = disregarded</span> · %s%s",
      fname, 100 * sup$frac_of_bbox, mode_lbl,
      if (use_3d && smode == "grid") " · plan view of a 3-D support (a column counts if any z slice does)" else "")
    fig %>% plotly::layout(
      shapes = shp,
      title = list(text = ttl, font = list(size = 11, color = "#555"), x = 0),
      xaxis = list(title = "x (nm)", zeroline = FALSE, constrain = "domain"),
      # Equal aspect: 1 nm in x is 1 nm in y, or the support's shape is a lie.
      yaxis = list(title = "y (nm)", zeroline = FALSE,
                   scaleanchor = "x", scaleratio = 1),
      margin = list(l = 60, r = 20, t = 34, b = 44),
      paper_bgcolor = "white", plot_bgcolor = "white") %>%
      plotly::config(doubleClick = "reset")
  }
  output$density_map <- renderPlotly(build_density_map_fig())

  output$download_density <- downloadHandler(
    filename = function() paste0("labelling_density_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      prof <- density_profile_rv()
      shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0, "No results."))
      pm <- density_p_map()
      if (!is.null(pm)) prof$keep_p <- unname(pm[prof$dataset])
      prof$target_density <- density_target() %||% NA_real_
      write.csv(prof, file, row.names = FALSE)
    }
  )

  # Keep the control-group selector populated from the group names.
  observe({
    gm <- file_group_map()
    gs <- unique(unname(gm))
    updateSelectInput(session, "dens_target_group", choices = gs,
                      selected = isolate(input$dens_target_group) %||% gs[1])
  })

  ## ── Support-aware observation window for the spatial statistics ────────────
  ## One switch, because it is one methodological choice: whether Ripley's L, g(r),
  ## L12, g12 and the co-localization null measure intensity over the SAMPLE SUPPORT
  ## or over the axis-aligned bounding box. See the SUPPORT-AWARE OBSERVATION WINDOW
  ## banner in helpers_smlm.R for the measured effect (a CSR pattern in an L-shaped
  ## cell read as clustered by +68.9 nm of L-r with the box, -13.5 nm with the
  ## support; co-localization false positives 17-37% -> 0%).
  ##
  ## The support settings are the Labelling Density panel's; if that panel has never
  ## been opened, these defaults apply so the switch works on its own.
  ## TWO SEPARATE CONCERNS, deliberately not bundled:
  ##
  ## `spatial_edge_correct` - always safe, and the bigger win. The bounding box is
  ##   a perfectly valid window; it is just missing its boundary correction. For a
  ##   cell that FILLS the scan field the box IS the cell, so this alone is the whole
  ##   fix (measured: clustered data in a 5 um cell, bbox peak L-r +196 vs the true
  ##   window's +202).
  ##
  ## `spatial_use_support` - only for a cell that genuinely does NOT fill the frame
  ##   (concave, or partly out of field). Measured on CSR in an L-shaped cell, the box
  ##   invents clustering (peak L-r +116 where the truth is 0) and the footprint
  ##   window removes it (+0.7).
  ##
  ## ⚠ The window must be a FOOTPRINT, never the density HDR. A density-derived
  ## support is "where the sample is dense", which on clustered data is the CLUSTERS,
  ## and using that as the window conditions on the very thing being tested: measured
  ## 29% of the cell, crushing peak L-r from +202 to +49 and dragging the tail
  ## negative - the "straight line into the negative" symptom. So this deliberately
  ## does NOT follow the Labelling Density panel's estimator; it always uses the
  ## occupancy grid at its closing plateau, which is a footprint by construction.
  ## ── ONE CLOSING RADIUS FOR THE WHOLE COMPARISON SET ────────────────────────
  ## The gap-fill radius used to be chosen per file, from that file's own closing
  ## sweep, by plateau_close_nm() = the FIRST sweep step growing under 5%. Both
  ## halves of that are wrong for a between-file comparison:
  ##
  ##   PER FILE. Two files of the same biology get different windows, so their
  ##   areas mean different things and K (which is proportional to the area) is not
  ##   comparable. Measured on one cell with 26 identical structures, varying ONLY
  ##   the localizations per structure: the auto radius came out 400, 1600 or
  ##   2400 nm and the window 12.6-36.8 um^2, a 2.9x spread on identical geometry.
  ##
  ##   FIRST-DIP. A patchy sample's closing curve is not monotone-then-flat, it is
  ##   bumpy: the area flattens once the STRUCTURES fill in, long before the gaps
  ##   BETWEEN them bridge. Whether a file dips under the tolerance early is luck
  ##   of its density, and the file that does gets a structure-hugging window and
  ##   lands next to the 50% rejection threshold (measured: largest_frac 54%).
  ##
  ## So the radius is now resolved ONCE for the whole set: the smallest swept
  ## radius at which EVERY file resolves to one dominant piece, still capped by
  ## frac_of_bbox so it cannot swallow the frame. Rejection then stops being
  ## something to fall back from - it is impossible by construction - and if no
  ## radius satisfies every file the comparison is REFUSED rather than run with
  ## mixed windows, which is not a comparison at all.
  SPATIAL_MAX_BBOX_FRAC <- 0.97   # same cap plateau_close_nm() uses

  ## Cost note: this runs ONE closing sweep per file - the same work the per-file
  ## path already did - and builds NO masks and runs NO connected-components. An
  ## earlier version searched a radius grid with a mask + component test per file
  ## per radius, which ran before any progress bar and looked like a hang.
  ##
  ## Why max-of-plateaus is the right shared radius: closing is MONOTONE, so a
  ## larger radius can never fragment a window that was already connected. Taking
  ## the largest per-file plateau therefore leaves every file at least as connected
  ## as it would have been on its own, while giving them all one radius.
  ## `force` and `h_nm` exist for the Local Density & Hotspots panel, which needs a
  ## footprint unconditionally (its statistics are all "per unit of THIS cell's own
  ## area", so a bounding box does not merely lose power, it changes the number) and
  ## renders at a finer grid pitch than the curve panels. Everything else calls this
  ## with the defaults and behaves exactly as before.
  spatial_resolve_close <- function(df_list, use_3d = FALSE, force = FALSE,
                                    h_nm = NULL) {
    if (!isTRUE(input$spatial_use_support) && !isTRUE(force))
      return(list(ok = TRUE, close_nm = NULL))
    h    <- max(10, as.numeric(h_nm %||% input$dens_grid_h %||% 100))
    smode <- input$dens_support_mode %||% "grid"
    keep <- if (use_3d) c("x", "y", "z") else c("x", "y")

    prep <- function(df) {
      if (is.null(df) || !nrow(df) || !all(keep %in% names(df))) return(NULL)
      sp <- df[, keep, drop = FALSE]
      if (!use_3d) sp$z <- 0
      sp <- sp[stats::complete.cases(sp), , drop = FALSE]
      if (nrow(sp) < 3L) NULL else sp
    }
    items <- Filter(function(it) !is.null(it$sp),
                    lapply(df_list, function(it)
                      list(sp = prep(it$df), tag = as.character(it$tag %||% ""))))
    if (!length(items)) return(list(ok = TRUE, close_nm = NULL))

    # Radius read off the sweep by the BRIDGING JUMP, not by plateau_close_nm().
    # plateau_close_nm() takes the first step growing under its tolerance, which on
    # a patchy sample is the moment the STRUCTURES fill in - long before the gaps
    # BETWEEN them close. It is a knife-edge: over 20 independent 6-file sets from
    # one generator it picked 400 nm in 8 of them and 2400 nm in the other 12, and
    # every 400 nm set left all six files fragmented. The largest relative jump in
    # the sweep IS the gaps bridging, so take the first settled radius after it -
    # stable at 2400 nm in 20/20, with nothing fragmented.
    bridge_close <- function(sw) {
      if (is.null(sw) || nrow(sw) < 3L) return(NA_real_)
      s  <- sw[order(sw$close_nm), , drop = FALSE]
      ok <- is.finite(s$frac_of_bbox) & s$frac_of_bbox <= SPATIAL_MAX_BBOX_FRAC
      s  <- s[ok, , drop = FALSE]
      if (nrow(s) < 3L) return(NA_real_)
      j <- which.max(s$rel_step[-1]) + 1L               # the bridging jump
      k <- which(seq_len(nrow(s)) > j & is.finite(s$rel_step) & s$rel_step < 0.10)
      if (length(k)) s$close_nm[k[1]] else s$close_nm[nrow(s)]
    }
    per <- vapply(items, function(it) tryCatch({
      p <- bridge_close(support_closing_sweep(it$sp, use_3d = use_3d, h_nm = h))
      if (is.finite(p)) p else NA_real_
    }, error = function(e) NA_real_), numeric(1))

    if (!any(is.finite(per)))
      return(list(ok = TRUE, close_nm = max(0, as.numeric(input$dens_close_nm %||% 200))))

    shared <- max(per, na.rm = TRUE)
    list(ok = TRUE, close_nm = shared,
         per_file = per,
         tags = vapply(items, function(it) it$tag, character(1)),
         # the "your files disagree about their own area scale" signal
         spread = if (sum(is.finite(per)) > 1)
                    max(per, na.rm = TRUE) / max(min(per[is.finite(per)]), 1e-9) else 1)
  }

  ## ONE WINDOW PER FILE, FOR EVERY PANEL THAT DIVIDES BY AN AREA.
  ## The estimator settings all come from "Density & area normalization" at the top
  ## of the tab, so Ripley's L, g(r) and Local Density & Hotspots measure the same
  ## cell the same way. Three panels resolving three different areas for one cell
  ## cannot be reconciled afterwards: K(r) scales with the window area and every
  ## hotspot statistic is per unit of it, so the numbers would silently disagree.
  ##
  ## `force` is for the hotspot panel, which needs a footprint UNCONDITIONALLY.
  ## Elsewhere the footprint is an opt-in that improves an estimate; there it is
  ## the definition of the quantity, since a bounding box makes a half-framed cell
  ## read at half its density and doubles every enhancement.
  ## WHICH POINTS DRAW THE OUTLINE is `dens_support_from`, shared with the density
  ## panel so all four windows match. Neither choice is free of a bias:
  ##
  ##   "units"  (default) - trace / molecule centroids. What the analyses actually
  ##            count, so the outline is the region those units occupy. Its bias is
  ##            that a sparsely labelled cell yields a raggeder, smaller footprint,
  ##            so the area carries some of the unit count (measured on synthetic
  ##            pairs: 5.01e6 vs 3.82e6 nm2 for identical geometry at a 4x labelling
  ##            difference, turning a 4.0x density ratio into 3.05x).
  ##
  ##   "locs" - every localization. A denser sample of the same region, so a
  ##            smoother outline, but trace LENGTH is photophysics (how long one
  ##            fluorophore stayed on), and it enters the closing sweep. Measured on
  ##            the two-condition set: the localization outline inflated the area enough to
  ##            loosen the CSR-calibrated eps, and the largest-focus size went from
  ##            cleanly separating the conditions (6-12 against 19-51 units) to
  ##            overlapping (10-24 against 19-52), by sweeping diffuse background
  ##            into the control's foci.
  ##
  ## Default is "units": hotspots are hotspots OF CENTROIDS, and the failure above
  ## is a false negative on the headline metric, which is worse than the density
  ## ratio compression.
  spatial_window_of <- function(df, use_3d = FALSE, close_nm = NULL, force = FALSE) {
    if (!isTRUE(input$spatial_edge_correct) && !isTRUE(input$spatial_use_support) &&
        !isTRUE(force))
      return(NULL)
    keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
    if (is.null(df) || !nrow(df) || !all(keep %in% names(df))) return(NULL)
    sp <- df[, keep, drop = FALSE]
    if (!use_3d) sp$z <- 0
    sp <- sp[stats::complete.cases(sp), , drop = FALSE]
    if (nrow(sp) < 3L) return(NULL)
    h <- max(10, as.numeric(input$dens_grid_h %||% 100))
    smode <- input$dens_support_mode %||% "grid"

    if (!isTRUE(input$spatial_use_support) && !isTRUE(force))
      return(bbox_support(sp, use_3d = use_3d, h_nm = h))

    # Shared radius when the caller resolved one for the set; the per-file plateau
    # remains only as the single-file / un-resolved fallback.
    cn <- if (!is.null(close_nm) && is.finite(close_nm)) close_nm else tryCatch({
      sw <- support_closing_sweep(sp, use_3d = use_3d, h_nm = h)
      p  <- plateau_close_nm(sw)
      if (is.finite(p)) p else max(0, as.numeric(input$dens_close_nm %||% 200))
    }, error = function(e) max(0, as.numeric(input$dens_close_nm %||% 200)))
    # The estimator mode is the density block's, not a hard-coded "grid". `hull`
    # and `bbox` return no mask, and the randomisation-based callers (the CSR
    # envelopes, the calibrated eps) cannot work without one, so those two fall
    # back to the bounding box here rather than half-working downstream.
    fp <- tryCatch(estimate_support(sp, use_3d = use_3d, mode = smode, h_nm = h,
                                   close_nm = cn, return_mask = TRUE,
                                   bw_nm = max(10, as.numeric(input$dens_kde_bw %||% 150)),
                                   frac  = min(1, max(0.5, as.numeric(input$dens_kde_frac %||% 0.99)))),
                   error = function(e) NULL)
    if (is.null(fp) || is.null(fp$mask)) return(bbox_support(sp, use_3d = use_3d, h_nm = h))

    # GUARD: reject a window that is many pieces of COMPARABLE size (structures),
    # not one that merely has a few detached slivers (a real cell edge). Measured on
    # four real MINFLUX acquisitions the footprint came out in 2/3/1/1 pieces with
    # >=93% of the window in the largest, so a count-based rule rejected good data.
    # The occupancy footprint stops inside the real cell edge at sparse labelling, so
    # a null simulated inside it is denser than the data and random points read as
    # regular (measured: L-r = -6.9 nm at r = 200 on CSR, 33% below the pointwise
    # band). Grow it back to the region the points actually fill. Not applied to the
    # hotspot panel (`force`), whose eps calibration randomises inside this same mask
    # and is therefore already self-consistent.
    if (!isTRUE(force)) {
      fp <- tryCatch(debias_support_window(fp, nrow(sp)), error = function(e) fp)
    }
    cmp <- mask_components(if (use_3d) fp$mask else mask_project_2d(fp$mask), fp$h)
    fp$n_components  <- cmp$n
    fp$largest_frac  <- cmp$largest_frac
    if (is.finite(cmp$largest_frac) && cmp$largest_frac < 0.5) {
      showNotification(sprintf(
        "Footprint rejected: %d pieces, %.0f%% of the area in the largest, so it traces the structures rather than the cell. Using the bounding box instead; check the support map.",
        cmp$n, 100 * cmp$largest_frac), type = "warning", duration = 12)
      bb <- bbox_support(sp, use_3d = use_3d, h_nm = h)
      bb$rejected_footprint <- TRUE
      return(bb)
    }
    fp
  }
  # Back-compatible name used by the analysis panels. Also records what each file's
  # window came out as, so the readout below can show it rather than the user having
  # to infer it from a warning that may or may not have fired.
  spatial_window_log <- reactiveVal(NULL)
  # The radius the current run resolved for the whole set, shown once in the readout.
  spatial_common_close <- reactiveVal(NULL)
  spatial_support_of <- function(df, use_3d = FALSE, tag = NULL, close_nm = NULL,
                                 force = FALSE) {
    w <- spatial_window_of(df, use_3d, close_nm = close_nm, force = force)
    if (!is.null(w) && !is.null(tag)) {
      # Edge correction on a bounding box is only right when the sample fills the box.
      # When the footprint window is off, measure (without using it) how much of the
      # box a footprint would cover, so the readout can prompt for the footprint
      # window. Uses the footprint at the LARGEST swept gap fill (24 grid steps), the
      # most generous outline: clustered data filling its frame reads ~80% there
      # (closing cannot fill indentations at the box edge), so only a box that is
      # mostly empty at the micrometre scale falls below the prompt threshold.
      # Plan view only; ~0.1 s per file.
      fp_frac <- NA_real_
      if (identical(w$mode, "bbox") && !isTRUE(input$spatial_use_support) && !isTRUE(force) &&
          isTRUE(input$spatial_edge_correct) && !is.null(df) && all(c("x", "y") %in% names(df))) {
        fp_frac <- tryCatch({
          sp <- df[stats::complete.cases(df[, c("x", "y")]), c("x", "y"), drop = FALSE]
          sp$z <- 0                                             # plan view
          h  <- max(10, as.numeric(input$dens_grid_h %||% 100))
          sw <- support_closing_sweep(sp, use_3d = FALSE, h_nm = h)
          if (is.null(sw) || !nrow(sw)) NA_real_
          else sw$frac_of_bbox[which.max(sw$close_nm)]          # most generous footprint
        }, error = function(e) NA_real_)
      }
      row <- data.frame(
        dataset      = tag,
        window       = if (identical(w$mode, "bbox"))
                         (if (isTRUE(w$rejected_footprint)) "bounding box (footprint rejected)"
                          else "bounding box") else "footprint",
        area_nm2     = w$support %||% NA_real_,
        frac_of_bbox = w$frac_of_bbox %||% 1,
        gap_fill_nm  = w$close_nm %||% NA_real_,
        pieces       = w$n_components %||% NA_integer_,
        largest_frac = w$largest_frac %||% NA_real_,
        footprint_frac = fp_frac,        # diagnostic only, filled when the box is in use
        stringsAsFactors = FALSE)
      cur <- spatial_window_log()
      cur <- if (is.null(cur)) row else rbind(cur[cur$dataset != tag, , drop = FALSE], row)
      spatial_window_log(cur)
    }
    w
  }

  ## Resolves the shared radius and stores it for the readout. It never blocks and
  ## it does NOT warn about the areas differing between files - different cells
  ## genuinely have different areas, and that is what dividing by the area is for.
  ## The thing that actually breaks a comparison is files being measured by
  ## DIFFERENT METHODS, which is only knowable once the windows are built; see
  ## spatial_report_windows(), called after the run.
  ##
  ## Shown under its own progress message: this is a couple of seconds of sweeping
  ## before any curve exists, and in silence it reads as a hang.
  spatial_gate_close <- function(df_list, use_3d = FALSE, force = FALSE,
                                 h_nm = NULL) {
    cc <- if (isTRUE(input$spatial_use_support) || isTRUE(force))
      withProgress(message = "Measuring sample area...", value = 0.5,
                   spatial_resolve_close(df_list, use_3d, force = force, h_nm = h_nm))
      else spatial_resolve_close(df_list, use_3d, force = force, h_nm = h_nm)
    spatial_common_close(cc)
    cc
  }

  ## Called AFTER a run, with the tags it produced. Warns only when the set was not
  ## measured the same way - i.e. some files got the footprint and others fell back
  ## to the bounding box. That mix is the real problem: the bounding box counts
  ## empty field as sample, K is proportional to the area, so the fallback file's
  ## curve is inflated relative to the rest and the sets are not on one scale.
  ## Deliberately posted after the run so it cannot sit on top of the progress bar.
  spatial_report_windows <- function(tags) {
    if (!isTRUE(input$spatial_use_support)) return(invisible(NULL))
    lg <- spatial_window_log()
    if (is.null(lg) || !nrow(lg)) return(invisible(NULL))
    lg <- lg[lg$dataset %in% tags, , drop = FALSE]
    if (nrow(lg) < 2L) return(invisible(NULL))
    fell_back <- grepl("bounding box", lg$window, fixed = TRUE)
    if (!any(fell_back) || all(fell_back)) return(invisible(NULL))
    showNotification(HTML(sprintf(
      "<b>These files were not measured the same way.</b> %s fell back to the bounding box while the rest used the footprint, so their curves are not on the same scale. Check the support map, or turn the footprint window off for all files.",
      paste(utils::head(lg$dataset[fell_back], 4), collapse = ", "))),
      type = "warning", duration = NULL)
    invisible(NULL)
  }

  output$spatial_window_readout <- renderUI({
    if (!isTRUE(input$spatial_edge_correct) && !isTRUE(input$spatial_use_support))
      return(helpText("Off - spatial statistics use the bounding box with no boundary correction."))
    lg <- spatial_window_log()
    if (is.null(lg) || !nrow(lg))
      return(helpText("Run Ripley's L or g(r) to see the window each file resolved to."))
    rows <- lapply(seq_len(nrow(lg)), function(i) {
      warn <- grepl("rejected", lg$window[i], fixed = TRUE)
      tags$tr(
        tags$td(style = "padding:1px 12px 1px 0;", lg$dataset[i]),
        tags$td(style = paste0("padding:1px 12px 1px 0;",
                               if (warn) "color:#b26a00;font-weight:600;" else ""), lg$window[i]),
        tags$td(style = "padding:1px 12px 1px 0;", sprintf("%.3g nm2", lg$area_nm2[i])),
        tags$td(style = "padding:1px 12px 1px 0;", sprintf("%.0f%% of box", 100 * lg$frac_of_bbox[i])),
        tags$td(style = "padding:1px 12px 1px 0; color:#777;",
                if (is.na(lg$pieces[i])) "-" else
                  sprintf("%d piece(s), %.0f%% in largest", lg$pieces[i], 100 * lg$largest_frac[i])))
    })
    cc <- spatial_common_close()
    small <- which(is.finite(lg$frac_of_bbox) & lg$frac_of_bbox < 0.40 &
                     !grepl("bounding box", lg$window, fixed = TRUE))
    # Edge-corrected bounding box on a sample that leaves much of the box empty: the
    # correction then amplifies the framing error, so prompt for the footprint window.
    # 0.70 is a prompt threshold (a sample filling 70% of its box inflates K by
    # ~1.4x). Synthetic check: CSR filling the frame 100%, clustered filling the frame
    # 79%, L-shaped cell 75%, so it fires on half-framed cells and diagonal rods.
    empty_box <- if ("footprint_frac" %in% names(lg))
                   which(is.finite(lg$footprint_frac) & lg$footprint_frac < 0.70 &
                         grepl("bounding box", lg$window, fixed = TRUE))
                 else integer(0)
    tags$div(style = "font-size:0.85em; margin-top:6px;",
      tags$div(style = "font-weight:600; margin-bottom:3px;", "Window used per file:"),
      tags$table(style = "border-collapse:collapse;", rows),
      if (!is.null(cc) && isTRUE(cc$ok) && !is.null(cc$close_nm) && is.finite(cc$close_nm))
        tags$p(style = "margin:5px 0 0; color:#2c3e50;",
          tags$b(sprintf("Gap fill: %.0f nm, shared by every file in this run.", cc$close_nm)),
          " One radius means every file's area is measured the same way. The areas themselves are",
          " expected to differ - that is what dividing by the area is for."),
      if (length(unique(lg$window)) > 1L)
        tags$p(style = "margin:5px 0 0; color:#b26a00;",
          tags$b("Mixed methods above. "),
          "Some files used the footprint and others the bounding box. The bounding box counts empty",
          " field as sample and K scales with the window area, so those curves are inflated relative",
          " to the rest and the two groups are not on one scale."),
      if (length(empty_box))
        tags$p(style = "margin:5px 0 0; color:#b26a00;",
          tags$b("Turn on the footprint window. "),
          sprintf("The sample covers about %s of the bounding box in %s. ",
                  paste(sprintf("%.0f%%", 100 * lg$footprint_frac[empty_box]), collapse = ", "),
                  paste(lg$dataset[empty_box], collapse = ", ")),
          "The box counts that empty field as sample, which reads as clustering, and edge",
          " correction on such a box makes the error larger."),
      if (length(small))
        tags$p(style = "margin:5px 0 0; color:#b26a00;",
          tags$b("Check the support map. "),
          sprintf("%s came out under 40%% of the bounding box. ",
                  paste(lg$dataset[small], collapse = ", ")),
          "That is expected if the cell genuinely occupies a small part of the field, but a window",
          " that small can instead be tracing the structures - which understates clustering,",
          " because a window drawn round the clusters conditions on the thing being measured."),
      tags$p(style = "margin:5px 0 0; color:#7a5c00;",
        "Each file's window is built from its own localizations, so the area estimation method is",
        " identical across files."))
  })

  ## ── Small-N reliability advisory for the per-structure geometry panels ──────
  ## Every shape metric is estimated from a cluster's own members and is biased at
  ## low member counts, in DIFFERENT directions, so a bare member count does not tell
  ## a reader how much to trust a value. Measured against known ground truth:
  ##   eccentricity  a TRUE SPHERE reads median ecc 0.81 at n = 10, 0.53 at n = 40,
  ##                 0.34 at n = 200, 0.18 at n = 2000 - PCA orders axes by sample
  ##                 variance, so the ordering itself manufactures elongation
  ##   alpha volume  29.9% of truth at n = 10, 54.6% at n = 20, 67.6% at n = 31,
  ##                 86.2% at n = 80, 96.0% at n = 300 - unusable below ~20
  ##   sphere radius 147.4 nm at n = 40 against a truth of 150; 149.9 at n = 1000
  ##   radius of gyration  Bessel-corrected here, but ~12% low at n = 2 from Jensen
  ## Advisory only: it never blanks a value, because which metric matters decides
  ## what counts as "enough" and hiding the numbers would hide the size dependence.
  .small_n_note <- function(nvec, what = "structures") {
    n <- nvec[is.finite(nvec)]
    if (!length(n)) return(NULL)
    br <- c(sum(n < 20), sum(n >= 20 & n < 40), sum(n >= 40 & n < 200), sum(n >= 200))
    pct <- round(100 * br / length(n))
    sev <- br[1] / length(n)
    col <- if (sev > 0.5) "#b26a00" else if (sev > 0.2) "#7a5c00" else "#777"
    tags$div(style = paste0("font-size:0.82em; color:", col,
                            "; margin-top:4px; line-height:1.45;"),
      tags$b(sprintf("Members per %s: ", what)),
      sprintf("%d%% under 20, %d%% 20-39, %d%% 40-199, %d%% 200+ (median %.0f).",
              pct[1], pct[2], pct[3], pct[4], stats::median(n)),
      if (sev > 0.2) tags$span(
        tags$br(),
        "Shape measurements get less reliable the fewer members a structure has, and they do not ",
        "all fail in the same direction. Enclosed volume comes out too small, and is not usable ",
        "below about 20 members. Eccentricity comes out too large: even a perfect sphere reads as ",
        "elongated when it is built from few points, because the axes are ordered by the scatter of ",
        "those points. Compare an eccentricity against the sphere floor for its own member count, ",
        "given in the exported table, rather than against zero.")
      else NULL)
  }

  ## ── Per-file prepared-data cache ────────────────────────────────────────────
  ## GOAL: run the filter + de-blink + median-clustering chain ONCE per file per
  ## settings combination, instead of once per file per PANEL.
  ##
  ## WHY: every batch analysis independently repeated
  ##   qc_exclude_from -> apply_meta_filters_batch -> apply_density_thinning
  ##   -> apply_clustering_pipeline
  ## on every loaded file. Measured on 4 real acquisitions (61,737 localizations),
  ## the panels that cluster took 6.1-11.3 s each while the three that do NOT cluster
  ## (trace detection 1.1 s, fiducial 0.8 s, precision cloud 4.0 s) were fast - so
  ## roughly 5 s of every clustering panel was re-deriving identical clusters.
  ##
  ## KEY DESIGN, and the thing to preserve if you touch this: the cache key is built
  ## from the FULL SET of inputs the chain consumes, enumerated in .PREP_KEY_INPUTS
  ## below. That list was derived by scanning apply_meta_filters_batch,
  ## apply_clustering_pipeline and resolve_channels for `input$` reads - it is not
  ## from memory. Under-keying is the dangerous failure (stale clusters served
  ## silently), so the list includes everything that DEFINES a channel as well as
  ## everything that clusters one, and `.prep_key_selftest()` re-scans those
  ## functions at runtime and warns if any input is consumed but missing.
  .PREP_KEY_INPUTS <- c(
    # apply_meta_filters_batch
    "efo_range", "cfr_range", "dcr_range", "eco_range", "tim_range", "dwell_range",
    "fluo_range", "tid_length_range", "trim_outliers", "trim_k", "sd_min_locs",
    "view_SD", "target_sd_x", "target_sd_y",
    # apply_clustering_pipeline (+ assign_analysis_group)
    "enable_clustering", "deblink_method", "eps_centroid", "deblink_anisotropic",
    "vote_frac", "vote_rho", "vote_link",
    "enable_cluster2", "eps2", "minpts2", "cluster2_method", "cluster2_xi",
    "cluster2_prefilter", "cluster2_prefilter_r", "cluster2_prefilter_min",
    "cluster2_prune", "cluster2_prune_tol", "cluster2_core_only",
    "cluster2_min_members", "cluster_group_scope", "centroid_k",
    "deblink_shared_spread",
    # resolve_channels - what DEFINES a channel. Getting one of these wrong
    # repartitions the data, so they are as important as the eps values.
    "enable_two_color", "color_by_fluo", "dcr_group1", "dcr_group2",
    "use_pooled_dcr",
    # the two opt-in gates that rewrite the data before clustering
    "fid_exclude", "dens_norm_enable", "dens_thin_enable", "dens_norm_unit"
  )

  prep_cache_rv <- reactiveVal(list())   # key -> prepared data.frame

  # Everything the prepared data depends on, collapsed to one comparable object.
  # Includes the thinning keep-probabilities and the fiducial exclusion list, which
  # are derived state rather than raw inputs and would otherwise be missed.
  .prep_key <- function(fname, needs_cluster, two_color) {
    vals <- lapply(.PREP_KEY_INPUTS, function(k) input[[k]])
    names(vals) <- .PREP_KEY_INPUTS
    # density_p_map() is NULL until the density panel has been run; signif() on NULL
    # errors, which silently made every prepared_df() call fail and every panel
    # return no data.
    pm <- tryCatch(density_p_map(), error = function(e) NULL)
    pm_txt <- if (is.numeric(pm) && length(pm))
      paste(names(pm), signif(unname(pm), 8), collapse = ",") else ""

    fr <- tryCatch(fiducial_results(), error = function(e) NULL)
    ex <- if (!is.null(fr) && all(c("flagged", "tid", "dataset") %in% names(fr)))
      sort(fr$tid[fr$flagged %in% TRUE & fr$dataset == fname]) else NULL
    paste(fname, needs_cluster, two_color,
          paste(utils::capture.output(utils::str(vals, max.level = 2)), collapse = "|"),
          pm_txt,
          paste(ex, collapse = ","), sep = "##")
  }

  ## Prepared data for one file: cached filter -> thin -> cluster chain.
  ## `needs_cluster = FALSE` skips clustering entirely (the panels that only need
  ## trace positions), and is cached separately so it cannot be served a clustered
  ## frame or vice versa.
  prepared_df <- function(fname, needs_cluster = TRUE, two_color = NULL) {
    store <- smlm_file_store()
    if (!fname %in% names(store)) return(NULL)
    if (is.null(two_color)) two_color <- isTRUE(input$enable_two_color)
    key <- .prep_key(fname, needs_cluster, two_color)
    cache <- prep_cache_rv()
    if (!is.null(cache[[key]])) return(cache[[key]])

    df <- qc_exclude_from(store[[fname]], fname)
    if (!is.data.frame(df) || !nrow(df)) return(NULL)
    df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- apply_density_thinning(df, fname)
    if (is.null(df) || !nrow(df)) return(NULL)
    if (isTRUE(needs_cluster)) {
      # ONE pipeline. `two_color` survives only as part of the cache key, because
      # enabling the 2-colour split repartitions the data and so must invalidate
      # a cached frame - it no longer selects between two implementations.
      df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) return(NULL)
    }
    # Bounded: one entry per (file x settings). Settings churn would otherwise grow
    # this without limit over a long session.
    if (length(cache) > 60L) cache <- list()
    cache[[key]] <- df
    prep_cache_rv(cache)
    df
  }

  # Warn (once) if the chain reads an input the cache key does not cover - the
  # under-keying failure mode, which would serve stale clusters silently.
  .prep_key_selftest <- function() {
    # attach_channel_colors() is deliberately NOT scanned: it is applied at
    # render time, outside the cached chain, precisely so its inputs (the colour
    # pickers) do not have to be in the key.
    fns <- c("apply_meta_filters_batch", "apply_clustering_pipeline",
             "resolve_channels", "deblink_threshold_for_channel")
    used <- unique(unlist(lapply(fns, function(f) {
      if (!exists(f)) return(character(0))
      txt <- paste(deparse(get(f)), collapse = "
")
      m <- gregexpr("input[$][A-Za-z0-9_.]+", txt)
      sub("^input[$]", "", unlist(regmatches(txt, m)))
    })))
    miss <- setdiff(used, .PREP_KEY_INPUTS)
    if (length(miss))
      warning("prepared_df cache key is missing inputs the clustering chain reads: ",
              paste(miss, collapse = ", "), call. = FALSE)
    invisible(miss)
  }
  observe({ .prep_key_selftest() })

  ## ── "Results are stale" indicator ───────────────────────────────────────────
  ## Compute controls deliberately take effect only on Run, which is what stops an
  ## accidental scroll through r max from recomputing anything. The cost is that a
  ## changed control looks ignored, so each panel exposes a small badge comparing the
  ## live compute settings against the ones its results were computed with.
  ## Display-only controls (loess smoothing, bin counts, log axis, trend lines,
  ## per-file overlays) are NOT listed here - those stay live on purpose.
  .stale_badge <- function(live, used) {
    if (is.null(used)) return(NULL)
    live <- live[!vapply(live, is.null, logical(1))]
    used <- used[!vapply(used, is.null, logical(1))]
    keys <- intersect(names(live), names(used))
    diff <- keys[!vapply(keys, function(k)
      isTRUE(all.equal(live[[k]], used[[k]])), logical(1))]
    if (!length(diff)) return(NULL)
    tags$div(style = paste0("margin:4px 0; padding:4px 9px; border-left:4px solid #e67e22;",
                            "background:#fff8f0; font-size:0.82em; color:#7a5c00;"),
      tags$b("Settings changed since this was computed"), " - press Run to apply: ",
      paste(diff, collapse = ", "), ".")
  }

  output$ripley_stale <- renderUI(.stale_badge(
    list(unit = input$ripley_unit, `r max` = input$ripley_r_max,
         `r steps` = input$ripley_r_steps, `3-D` = isTRUE(input$ripley_use_3d),
         `all files` = isTRUE(input$ripley_all_files),
         envelope = isTRUE(input$ripley_envelope), simulations = input$ripley_nsim,
         `edge correction` = isTRUE(input$spatial_edge_correct),
         `support window` = isTRUE(input$spatial_use_support)),
    ripley_params()))

  output$pcf_stale <- renderUI(.stale_badge(
    list(unit = input$pcf_unit, `r max` = input$pcf_r_max,
         `r steps` = input$pcf_r_steps, `3-D` = isTRUE(input$pcf_use_3d),
         `all files` = isTRUE(input$pcf_all_files),
         envelope = isTRUE(input$pcf_envelope), simulations = input$pcf_nsim,
         `edge correction` = isTRUE(input$spatial_edge_correct),
         `support window` = isTRUE(input$spatial_use_support)),
    pcf_params()))

  output$nn_batch_stale <- renderUI(.stale_badge(
    list(unit = input$nn_batch_unit, k = input$nn_batch_k), nn_batch_params_live()))

  ## ── De-blinking is a GLOBAL setting: resolve every panel's unit against it ──
  ## The per-panel "Unit" selectors used to be independent of `enable_clustering`,
  ## which produced two silent failures:
  ##   1. de-blinking ON but a panel left on trace IDs -> that panel quietly ignored
  ##      de-blinking, so different panels reported on different units;
  ##   2. a molecule unit chosen while de-blinking is OFF -> apply_clustering_pipeline
  ##      never creates a `cluster` column, every file yields no points, and the
  ##      panel came back empty with no explanation. (Verified on two datasets:
  ##      Ripley, g(r) and Batch NN all returned nothing.)
  ## `resolve_unit()` maps the "auto" choice - now the default everywhere - onto the
  ## molecule unit when de-blinking is on and trace IDs otherwise, so switching
  ## de-blinking on once applies it to every analysis and every dataset. An explicit
  ## molecule unit with de-blinking off is downgraded to trace IDs and warned about
  ## rather than silently returning nothing.
  ## ONE SWITCH PAIR. There used to be two: the single-colour pair
  ## (enable_clustering / enable_cluster2) and a 2-colour pair
  ## (enable_two_color_cluster1 / enable_two_color_cluster2), with prepared_df()
  ## running whichever pipeline matched the mode. That meant a panel gating on
  ## `enable_cluster2` refused to run for a 2-colour user who had 2-colour Median
  ## Clustering enabled, and told them to tick a checkbox that would have had no
  ## effect. There is one pipeline now, so there is one pair of switches and
  ## nothing left to resolve between them.
  two_color_on <- reactive(isTRUE(input$enable_two_color))
  deblink_on   <- reactive(isTRUE(input$enable_clustering))
  cluster2_on  <- reactive(isTRUE(input$enable_cluster2))

  deblink_label  <- reactive("de-blinking")
  cluster2_label <- reactive("Median Clustering")

  resolve_unit <- function(sel, mol_value, tid_value = "tid", panel = "This panel") {
    sel <- sel %||% "auto"
    if (identical(sel, "auto")) return(if (deblink_on()) mol_value else tid_value)
    if (identical(sel, mol_value) && !deblink_on()) {
      showNotification(sprintf(
        "%s: de-blinked molecules requested but de-blinking is off - using trace IDs. Enable %s in the Display tab to group blinks into molecules.",
        panel, deblink_label()), type = "warning", duration = 7)
      return(tid_value)
    }
    sel
  }

  # ------ Batch NN Distance Histogram ----------------------------------------
  nn_batch_results <- reactiveVal(NULL)
  nn_batch_params  <- reactiveVal(NULL)   # RESOLVED compute params, frozen at Run
  nn_batch_params_live <- reactiveVal(NULL)  # raw selector values, for the stale badge

  # Two-color group filter - only visible when two-color is enabled
  output$nn_batch_group_ui <- renderUI({
    if (!isTRUE(input$enable_two_color)) return(NULL)
    selectInput("nn_batch_group", "Analysis group",
                choices = c("All"               = "all",
                            "Within Group 1"    = "within_1",
                            "Within Group 2"    = "within_2",
                            "Between Ch1 → Ch2" = "between_12",
                            "Between Ch2 → Ch1" = "between_21"),
                selected = "all")
  })

  observeEvent(input$run_batch_nn, {
    unit  <- resolve_unit(input$nn_batch_unit, "pcid", "tid", "Batch NN")
    k     <- max(1L, as.integer(input$nn_batch_k %||% 1L))
    gf    <- input$nn_batch_group %||% "all"
    group_filter <- switch(gf, "within_1" = 1L, "within_2" = 2L, NA_integer_)
    cross_groups <- switch(gf,
      "between_12" = c(1L, 2L), "between_21" = c(2L, 1L), NULL)

    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }

    # ── Multiplex NN modes ──────────────────────────────────────────────────
    adf      <- uploaded_df_rv()
    mplex_on <- !is.null(adf) && "mplex_label" %in% names(adf)
    mmode    <- if (mplex_on) (input$nn_batch_mplex_mode %||% "split") else "off"
    chx      <- input$nn_batch_chx %||% ""
    chy      <- input$nn_batch_chy %||% ""
    mplex_groups <- FALSE   # analysis_group derived from channel membership (between)

    if (mmode == "split") {
      cs <- mplex_split_active()
      if (length(cs)) store <- cs
    } else if (mmode == "within") {
      if (!nzchar(chx)) { showNotification("Pick Channel X.", type = "warning"); return() }
      sub <- adf[!is.na(adf$mplex_label) & adf$mplex_label == chx, , drop = FALSE]
      if (!nrow(sub)) { showNotification("Channel X has no data.", type = "warning"); return() }
      store <- stats::setNames(list(sub), chx)
      group_filter <- NA_integer_; cross_groups <- NULL
    } else if (mmode == "between") {
      if (!nzchar(chx) || !nzchar(chy) || identical(chx, chy)) {
        showNotification("Pick two distinct channels (X and Y).", type = "warning"); return()
      }
      store <- stats::setNames(list(adf), paste0(chx, " → ", chy))
      group_filter <- NA_integer_; cross_groups <- c(1L, 2L)
      mplex_groups <- TRUE
    }

    needs_clustering <- unit != "tid"
    two_color_mode   <- isTRUE(input$enable_two_color)
    needs_grouping   <- !is.na(group_filter) || !is.null(cross_groups)

    results   <- list()
    n_skipped <- 0L
    withProgress(message = "Running batch NN...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- qc_exclude_from(store[[fname]], fname)
        if (!is.data.frame(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }

        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }
        df <- apply_density_thinning(df, fname)   # opt-in labelling-density match
        if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }

        if (mplex_groups) {
          # Between-channel: cluster (within-channel) if needed, then derive
          # analysis_group from channel membership (X→1, Y→2) for compute_batch_nn.
          if (needs_clustering)
            df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) df)
          if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }
          df$analysis_group <- dplyr::case_when(
            df$mplex_label == chx ~ 1L,
            df$mplex_label == chy ~ 2L,
            TRUE ~ NA_integer_)
        } else if (two_color_mode && needs_grouping) {
          # Same pipeline as every other branch; resolve_channels() inside it
          # derives analysis_group from the 2-colour split.
          df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
          if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }
        } else if (needs_clustering) {
          df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
          if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }
        } else if (needs_grouping) {
          df <- tryCatch(assign_analysis_group(df, input), error = function(e) df)
        }

        # Per-channel facet. A CROSS-channel run (mplex_groups, or a group
        # filter) is never faceted: it is a statistic ABOUT the pair, so
        # splitting it by channel would ask a question that has no answer.
        facets <- if (mplex_groups || needs_grouping) list(list(label = "", df = df))
                  else channel_facets(df)
        got <- FALSE
        for (fc in facets) {
          d <- tryCatch(
            compute_batch_nn(fc$df, unit = unit, k = k, group_filter = group_filter,
                             cross_groups = cross_groups),
            error = function(e) data.frame()
          )
          if (!nrow(d)) next
          nm2 <- .facet_name(fname, fc$label)
          d$dataset <- nm2
          d$tag     <- if (".dataset_tag" %in% names(df)) df$.dataset_tag[[1L]] else fname
          results[[nm2]] <- d
          got <- TRUE
        }
        if (!got) { n_skipped <- n_skipped + 1L; next }
      }
    })

    if (!length(results)) {
      hint <- switch(unit,
        pcid          = "Check that Primary point DBSCAN Clustering is enabled and eps/minPts produce clusters.",
        cluster2_tid  = "Check that Median Clustering is enabled.",
        cluster2_pcid = "Check that both clustering levels are enabled.",
        "Check filters - all rows may have been filtered out."
      )
      showNotification(paste0("No distances computed. ", hint), type = "warning")
      nn_batch_results(NULL)
      return()
    }
    nn_batch_results(do.call(rbind, results))
    # Snapshot the compute parameters so the plot labels the data it is actually
    # showing. Reading them live meant changing the Unit selector relabelled an
    # old plot without recomputing it - and with unit = "auto" the axis read
    # "[auto]" because the label switch has no such case.
    nn_batch_params(list(unit = unit, k = k))
    nn_batch_params_live(list(unit = input$nn_batch_unit, k = input$nn_batch_k))
    msg <- sprintf("NN computed across %d of %d dataset(s).", length(results), length(store))
    if (n_skipped > 0)
      msg <- paste0(msg, sprintf(" %d file(s) yielded no distances (filtered out or no clusters).", n_skipped))
    showNotification(msg, type = "message", duration = 5)
  })

  # Shared reactive: stored results filtered by the distance range cutoff.
  # Both the histogram and the stats table read from this so they stay in sync.
  nn_batch_filtered <- reactive({
    res <- nn_batch_results()
    if (is.null(res) || !nrow(res)) return(NULL)
    dr <- input$nn_batch_dist_range %||% c(0, Inf)
    lo <- if (is.finite(dr[1])) dr[1] else 0
    hi <- if (is.finite(dr[2])) dr[2] else Inf
    res <- res[res$distance >= lo & res$distance <= hi, , drop = FALSE]
    # Keep the source FILE name alongside the (possibly group-relabelled) dataset,
    # so the plot can draw one series per file and colour it by its group.
    if (nrow(res)) res$file <- res$dataset
    apply_dataset_groups(res)   # relabel dataset -> group in group mode
  })

  build_nn_batch_hist_fig <- function() {
    res <- nn_batch_filtered()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Press 'Run Batch NN' to generate histogram (or widen the distance filter)."))

    bins       <- max(5L, as.integer(input$nn_batch_bins %||% 60L))
    density    <- identical(input$nn_batch_yaxis %||% "count", "density")
    error_type <- input$nn_batch_error %||% "iqr"
    aggregate  <- isTRUE(input$nn_batch_aggregate)
    col_g1     <- input$color_group1 %||% "#1f77b4"
    col_g2     <- input$color_group2 %||% "#d62728"
    group_col  <- switch(input$nn_batch_group %||% "all",
      within_1   = col_g1,
      within_2   = col_g2,
      between_12 = col_g1,
      between_21 = col_g2,
      "#2166ac"
    )

    # shared bins across all series so overlaid histograms + indicator lines align
    rng_all <- range(res$distance, na.rm = TRUE)
    if (!all(is.finite(rng_all)) || diff(rng_all) <= 0) rng_all <- c(rng_all[1], rng_all[1] + 1)
    bsize   <- diff(rng_all) / bins
    brk     <- seq(rng_all[1], rng_all[2], length.out = bins + 1L)
    xbins   <- list(start = rng_all[1], end = rng_all[2], size = bsize)
    bin_hts <- function(x) {
      h <- graphics::hist(x, breaks = brk, plot = FALSE, include.lowest = TRUE)
      if (density) h$density else h$counts
    }
    subsets <- if (aggregate) list(res$distance) else split(res$distance, res$dataset)
    ytop <- suppressWarnings(max(vapply(subsets,
              function(x) if (length(x)) max(bin_hts(x)) else 0, numeric(1)), na.rm = TRUE))
    if (!is.finite(ytop) || ytop <= 0) ytop <- 1

    # median + IQR/MAD as vertical LINE TRACES sharing the histogram's legendgroup
    # (so they toggle away with the group when its legend entry is clicked)
    add_spread <- function(fig, d, col, lg) {
      med <- stats::median(d)
      sv  <- switch(error_type,
        iqr = c(stats::quantile(d, 0.25, names = FALSE), stats::quantile(d, 0.75, names = FALSE)),
        mad = { mv <- stats::mad(d, constant = 1); c(max(0, med - mv), med + mv) },
        NULL)
      fig <- fig %>% plotly::add_trace(x = c(med, med), y = c(0, ytop),
        type = "scatter", mode = "lines", legendgroup = lg, showlegend = FALSE,
        hoverinfo = "skip", opacity = 0.9, line = list(color = col, width = 2))
      for (sv_i in sv) fig <- fig %>% plotly::add_trace(x = c(sv_i, sv_i), y = c(0, ytop),
        type = "scatter", mode = "lines", legendgroup = lg, showlegend = FALSE,
        hoverinfo = "skip", opacity = 0.65, line = list(color = col, width = 1.2, dash = "dot"))
      fig
    }

    fig <- plotly::plot_ly()

    if (aggregate) {
      fig <- fig %>% plotly::add_histogram(
        x        = res$distance,
        name     = sprintf("All datasets (N = %d)", nrow(res)),
        legendgroup = "all",
        marker   = list(color = paste0(group_col, "99")),
        xbins    = xbins,
        histnorm = if (density) "probability density" else "",
        opacity  = 0.75
      )
      fig <- add_spread(fig, res$distance, group_col, "all")
    } else {
      grouped_mode <- (input$analysis_by %||% "file") == "group"
      # ONE SERIES PER FILE in both modes. Colour encodes what distinguishes them:
      #   file mode  - a categorical colour per file (the old code shaded ONE blue
      #                by alpha, which made every file a shade of blue);
      #   group mode - the hue is the group and the shade is the file within it, so
      #                the conditions separate at a glance while replicates stay
      #                individually visible.
      files <- unique(res$file %||% res$dataset)
      fpal  <- file_palette(files)
      # Draw grouped by condition so the legend lists a group's files together.
      if (grouped_mode) files <- files[order(match(files, names(file_group_map())))]

      for (fl in files) {
        sub <- res[(res$file %||% res$dataset) == fl, , drop = FALSE]
        if (!nrow(sub)) next
        col_base <- unname(fpal[fl])
        if (is.na(col_base)) col_base <- group_col
        lgroup   <- if (grouped_mode) sub$dataset[[1L]] else fl
        lbl      <- if (grouped_mode) sprintf("%s - %s", sub$dataset[[1L]], fl)
                    else if ("tag" %in% names(sub)) sub$tag[[1L]] else fl
        fig <- fig %>% plotly::add_histogram(
          x        = sub$distance,
          name     = lbl,
          legendgroup = lgroup,
          marker   = list(color = paste0(col_base, "99")),
          xbins    = xbins,
          histnorm = if (density) "probability density" else "",
          opacity  = 0.75
        )
        fig <- add_spread(fig, sub$distance, col_base, lgroup)
      }
    }

    indicator_note <- switch(error_type,
      iqr  = "- median  ··· Q25 / Q75",
      mad  = "- median  ··· ±MAD",
      "- median"
    )
    # From the frozen snapshot (see nn_batch_params), never the live selector.
    npar <- nn_batch_params()
    u_used <- npar$unit %||% "tid"
    k_used <- npar$k    %||% 1L
    unit_lbl <- switch(u_used,
      tid           = "TID centroid",
      pcid          = "De-blinked molecule centroid",
      cluster2_tid  = "TID-within-Median-Cluster",
      cluster2_pcid = "De-blinked molecule centroid within Median Cluster",
      u_used
    )
    fig %>% plotly::layout(
      barmode = "overlay",
      xaxis   = list(title = sprintf("k=%d NN distance (nm)  [%s]",
                                     as.integer(k_used), unit_lbl)),
      yaxis   = list(title = if (density) "Density" else "Count"),
      legend  = list(orientation = "h", y = -0.18, groupclick = "togglegroup"),
      annotations = list(list(
        xref = "paper", yref = "paper", x = 0.01, y = 0.99,
        xanchor = "left", yanchor = "top", showarrow = FALSE,
        text = indicator_note, font = list(size = 11, color = "grey40")
      )),
      paper_bgcolor = "white",
      plot_bgcolor  = "white"
    )
  }
  output$nn_batch_hist <- renderPlotly(build_nn_batch_hist_fig())

  output$nn_batch_stats <- renderTable({
    res        <- nn_batch_filtered()
    error_type <- input$nn_batch_error %||% "iqr"
    aggregate  <- isTRUE(input$nn_batch_aggregate)
    if (is.null(res) || !nrow(res)) return(NULL)

    make_stats_row <- function(d, file_lbl, ds_lbl) {
      med <- stats::median(d)
      spread <- switch(error_type,
        iqr = {
          q25 <- stats::quantile(d, 0.25, names = FALSE)
          q75 <- stats::quantile(d, 0.75, names = FALSE)
          data.frame(Q25 = round(q25, 1), Q75 = round(q75, 1),
                     IQR = round(q75 - q25, 1), check.names = FALSE)
        },
        mad = {
          mv <- stats::mad(d, constant = 1)
          data.frame(`Median - MAD` = round(med - mv, 1),
                     `Median + MAD` = round(med + mv, 1),
                     MAD = round(mv, 1), check.names = FALSE)
        },
        data.frame()
      )
      cbind(data.frame(File = file_lbl, Dataset = ds_lbl, N = length(d),
                       `Median (nm)` = round(med, 1), check.names = FALSE), spread)
    }

    if (aggregate) {
      return(make_stats_row(res$distance, "All", "Aggregated"))
    }

    datasets <- unique(res$dataset)
    rows <- lapply(datasets, function(ds) {
      sub <- res[res$dataset == ds, , drop = FALSE]
      lbl <- if ("tag" %in% names(sub)) sub$tag[[1L]] else ds
      make_stats_row(sub$distance, ds, lbl)
    })
    do.call(rbind, rows)
  }, striped = TRUE, hover = TRUE, spacing = "s", digits = 1)

  output$download_batch_nn <- downloadHandler(
    filename = function()
      paste0("batch_nn_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- isolate(nn_batch_filtered())
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No NN results."))
      keep <- intersect(c("dataset", "tag", "distance"), names(res))
      write.csv(res[, keep, drop = FALSE], file, row.names = FALSE)
    }
  )

  # ------ Traces per Cluster -------------------------------------------------
  tid_count_results <- reactiveVal(NULL)
  tid_count_params  <- reactiveVal(NULL)   # compute params, frozen at Run

  observeEvent(input$run_tid_count, {
    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }

    # The GROUPING is ALWAYS the median cluster (cluster2). The UNIT counted
    # inside each median cluster is either raw trace IDs ("tid") or primary
    # de-blinked molecules ("pcid" - needs de-blinking enabled).
    unit      <- input$tid_count_unit %||% "tid"
    count_col <- if (unit == "pcid") "cluster" else "tid"

    all_counts <- list()
    withProgress(message = "Counting units per median cluster...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
        if (is.null(df) || !nrow(df)) next
        # Need median clusters (cluster2) to group by.
        if (!"cluster2" %in% names(df) || all(is.na(df$cluster2) | df$cluster2 == 0L)) next
        # Counting de-blinked molecules needs the de-blink (primary `cluster`) column present.
        if (count_col == "cluster" &&
            (!"cluster" %in% names(df) || all(is.na(df$cluster) | df$cluster == 0L))) next

        # Split after preparation, so a 2-colour split counts as channels too.
        # Every Median Cluster already belongs to exactly one channel (the
        # clustering ran inside channels), so this partitions the clusters
        # without changing any count.
        parts <- if (isTRUE(input$tid_count_by_channel)) channel_parts(df, fname)
                 else list(list(df = df, tag = fname))

        for (pt in parts) {
          sub <- pt$df %>% dplyr::filter(!is.na(cluster2) & cluster2 > 0L)
          # Legacy guard: de-blinking assigns every trace, so cluster 0 should not occur.
          if (count_col == "cluster") sub <- sub %>% dplyr::filter(!is.na(cluster) & cluster > 0L)
          if (!nrow(sub)) next
          cnt <- sub %>%
            dplyr::group_by(cluster2) %>%
            dplyr::summarise(n_tids = dplyr::n_distinct(.data[[count_col]]), .groups = "drop") %>%
            dplyr::mutate(dataset = pt$tag)
          all_counts[[pt$tag]] <- cnt
        }
      }
    })
    counts <- do.call(rbind, Filter(Negate(is.null), all_counts))

    if (is.null(counts) || !nrow(counts)) {
      hint <- if (count_col == "cluster")
        sprintf("Enable BOTH %s and %s in the Display tab, then press Plot for each file.",
                deblink_label(), cluster2_label())
      else
        sprintf("Enable %s in the Display tab, then press Plot for each file.", cluster2_label())
      showNotification(paste0("No median clusters found. ", hint), type = "warning")
      tid_count_results(NULL)
      return()
    }
    tid_count_results(counts)
    # Compute parameters frozen at Run: reading them live relabelled an existing
    # plot without recomputing it, so the axis could describe a different unit
    # than the data on screen.
    tid_count_params(list(unit = unit))
    unit_lbl <- if (count_col == "cluster") "de-blinked molecules" else "trace IDs"
    showNotification(
      sprintf("Counted %s in %d median clusters across %d file(s) (median = %g).",
              unit_lbl, nrow(counts), length(unique(counts$dataset)),
              stats::median(counts$n_tids)),
      type = "message", duration = 4)
  })

  build_tid_count_hist_fig <- function() {
    counts <- tid_count_results()
    shiny::validate(shiny::need(!is.null(counts) && nrow(counts) > 0,
                                "Press 'Run' to count traces per cluster."))
    counts <- apply_dataset_groups(counts)
    shiny::validate(shiny::need(nrow(counts) > 0, "No files assigned to groups (group mode)."))

    bins      <- max(5L, as.integer(input$tid_count_bins %||% 30L))
    density   <- identical(input$tid_count_yaxis %||% "count", "density")
    unit      <- (tid_count_params()$unit) %||% (input$tid_count_unit %||% "tid")
    # Grouping is always the median cluster; the unit is what we counted inside it.
    unit_lbl  <- if (unit == "pcid") "De-blinked molecules" else "Trace IDs"
    aggregate <- isTRUE(input$tid_count_aggregate)

    med_n   <- stats::median(counts$n_tids)
    grouped <- (input$analysis_by %||% "file") == "group" && !aggregate
    # y-extent for per-group median line-traces (approx histogram top)
    rng_c <- range(counts$n_tids, na.rm = TRUE)
    if (!all(is.finite(rng_c)) || diff(rng_c) <= 0) rng_c <- c(rng_c[1], rng_c[1] + 1)
    brk_c <- seq(rng_c[1], rng_c[2], length.out = bins + 1L)
    ytop_c <- suppressWarnings(max(vapply(
      if (grouped) split(counts$n_tids, counts$dataset) else list(counts$n_tids),
      function(x) { h <- graphics::hist(x, breaks = brk_c, plot = FALSE, include.lowest = TRUE)
                    max(if (density) h$density else h$counts) }, numeric(1)), na.rm = TRUE))
    if (!is.finite(ytop_c) || ytop_c <= 0) ytop_c <- 1
    fig   <- plotly::plot_ly()

    if (aggregate) {
      fig <- fig %>% plotly::add_histogram(
        x        = counts$n_tids,
        name     = sprintf("All datasets (N clusters = %d)", nrow(counts)),
        nbinsx   = bins,
        histnorm = if (density) "probability density" else "",
        marker   = list(color = "#4C8FBF99", line = list(color = "white", width = 0.5)),
        opacity  = 0.75
      )
    } else {
      datasets <- if ("dataset" %in% names(counts)) unique(counts$dataset) else NULL
      n_ds     <- length(datasets)
      pal <- if (n_ds >= 1) series_palette(datasets) else NULL

      if (!is.null(pal)) {
        for (ds in datasets) {
          sub <- counts[counts$dataset == ds, , drop = FALSE]
          fig <- fig %>% plotly::add_histogram(
            x        = sub$n_tids,
            name     = ds,
            legendgroup = ds,
            nbinsx   = bins,
            histnorm = if (density) "probability density" else "",
            marker   = list(color = paste0(pal[ds], "99"),
                            line  = list(color = pal[ds], width = 0.5)),
            opacity  = 0.75
          )
        }
        # per-group median line-traces (toggle away with their group)
        if (grouped) {
          for (ds in datasets) {
            m <- stats::median(counts$n_tids[counts$dataset == ds], na.rm = TRUE)
            fig <- fig %>% plotly::add_trace(x = c(m, m), y = c(0, ytop_c),
              type = "scatter", mode = "lines", legendgroup = ds, showlegend = FALSE,
              hoverinfo = "skip", line = list(color = pal[ds], width = 2, dash = "dash"))
          }
        }
      } else {
        fig <- fig %>% plotly::add_histogram(
          x        = counts$n_tids,
          nbinsx   = bins,
          histnorm = if (density) "probability density" else "",
          marker   = list(color = "#4C8FBF", line = list(color = "white", width = 0.5)),
          name     = "Clusters"
        )
      }
    }

    # overall median line/label only when not grouped (grouped shows per-group lines)
    ov_shapes <- if (!grouped) list(list(
        type = "line", x0 = med_n, x1 = med_n, y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "#E53935", dash = "dash", width = 1.5))) else list()
    ov_annot  <- if (!grouped) list(list(
        x = med_n, y = 1, yref = "paper", xanchor = "left",
        text = paste0(" median = ", round(med_n, 1)), showarrow = FALSE,
        font = list(color = "#E53935", size = 12))) else list()
    fig %>%
      plotly::layout(
        barmode = "overlay",
        shapes = ov_shapes,
        annotations = ov_annot,
        xaxis = list(title = paste0(unit_lbl, " per Median Cluster")),
        yaxis = list(title = if (density) "Density" else "Number of median clusters"),
        legend = list(orientation = "h", y = -0.18, groupclick = "togglegroup"),
        paper_bgcolor = "white",
        plot_bgcolor  = "white"
      )
  }
  output$tid_count_hist <- renderPlotly(build_tid_count_hist_fig())

  output$download_tid_count <- downloadHandler(
    filename = function()
      paste0("tid_per_cluster_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      counts <- tid_count_results()
      shiny::validate(shiny::need(!is.null(counts) && nrow(counts) > 0, "No results."))
      write.csv(counts, file, row.names = FALSE)
    }
  )

  ## ── Cluster2 Detection over Time ──────────────────────────────────────────
  # Each TID is one detection event placed at its median acquisition time
  # (tim_med), tagged in_c2 = whether it belongs to a Cluster2 structure. The
  # plot contrasts structure vs background detection over the acquisition.
  det_time_results <- reactiveVal(NULL)

  observeEvent(input$run_det_time, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    if (!cluster2_on()) {
      showNotification(sprintf("Enable %s in the Display tab first.", cluster2_label()),
                       type = "warning"); return()
    }
    if (!isTRUE(input$det_time_all_files)) {
      af <- input$active_file
      if (!is.null(af) && af %in% names(store)) store <- store[af]
    }

    all_tid <- list()
    withProgress(message = "Computing detection over time...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
        if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) next

        # Per channel when asked. This is the panel the multiplex boundary
        # question is read from: with a TIM-split build, the detection rate per
        # channel is what shows whether a channel is losing traces over its own
        # window rather than over the acquisition as a whole.
        parts <- if (isTRUE(input$det_time_by_channel)) channel_parts(df, fname)
                 else list(list(df = df, tag = fname))

        for (pt in parts) {
          per_tid <- pt$df %>%
            dplyr::group_by(tid) %>%
            dplyr::summarise(
              tim_med = stats::median(tim, na.rm = TRUE),
              in_c2   = any(!is.na(cluster2) & cluster2 > 0L),
              .groups = "drop"
            ) %>%
            dplyr::filter(is.finite(tim_med)) %>%
            dplyr::mutate(dataset = pt$tag)
          if (nrow(per_tid)) all_tid[[pt$tag]] <- per_tid
        }
      }
    })
    res <- do.call(rbind, Filter(Negate(is.null), all_tid))
    if (is.null(res) || !nrow(res)) {
      showNotification(sprintf("No traces with a tim column and Median cluster found. Enable %s and press Plot for each file.", cluster2_label()),
                       type = "warning")
      det_time_results(NULL); return()
    }
    det_time_results(res)
    showNotification(
      sprintf("Detection over time: %d traces across %d file(s); %.0f%% in Median Cluster",
              nrow(res), length(unique(res$dataset)), 100 * mean(res$in_c2)),
      type = "message", duration = 4)
  })

  build_det_time_plot_fig <- function() {
    res <- det_time_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Press 'Run' to compute detection over time."))
    res <- apply_dataset_groups(res)
    shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
    grouped <- (input$analysis_by %||% "file") == "group"

    bins <- max(5L, as.integer(input$det_time_bins %||% 40L))
    mode <- input$det_time_yaxis %||% "count"
    normalized <- mode %in% c("norm", "norm_cumulative")

    rng <- range(res$tim_med, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 0) rng <- c(rng[1], rng[1] + 1)
    breaks  <- seq(rng[1], rng[2], length.out = bins + 1L)
    ctr     <- (breaks[-1] + breaks[-(bins + 1L)]) / 2
    res$bin <- pmin(pmax(findInterval(res$tim_med, breaks, rightmost.closed = TRUE), 1L), bins)

    # per-subset y-values (respecting the selected mode) + KS gap
    yvals <- function(cnt_c2, cnt_bg) {
      cum_c2 <- cumsum(cnt_c2); cum_bg <- cumsum(cnt_bg)
      tot_c2 <- sum(cnt_c2);    tot_bg <- sum(cnt_bg)
      sc2 <- if (tot_c2 > 0) tot_c2 else 1; sbg <- if (tot_bg > 0) tot_bg else 1
      ks  <- if (tot_c2 > 0 && tot_bg > 0) max(abs(cum_c2 / tot_c2 - cum_bg / tot_bg)) else NA_real_
      v <- switch(mode,
        cumulative      = list(c2 = cum_c2,      bg = cum_bg,      yt = "Cumulative traces detected"),
        norm            = list(c2 = cnt_c2 / sc2, bg = cnt_bg / sbg, yt = "Detection fraction per bin"),
        norm_cumulative = list(c2 = cum_c2 / sc2, bg = cum_bg / sbg, yt = "Cumulative detection fraction"),
                          list(c2 = cnt_c2,      bg = cnt_bg,      yt = "Traces detected per bin"))
      c(v, list(ks = ks))
    }

    fig    <- plotly::plot_ly()
    ytitle <- "Traces detected per bin"
    ttl    <- NULL

    if (!grouped) {
      cnt_c2 <- tabulate(res$bin[res$in_c2],  nbins = bins)
      cnt_bg <- tabulate(res$bin[!res$in_c2], nbins = bins)
      yy <- yvals(cnt_c2, cnt_bg); ytitle <- yy$yt
      fig <- fig %>%
        plotly::add_trace(x = ctr, y = yy$c2, type = "scatter", mode = "lines+markers",
                          line = list(color = "#2C7FB8"), marker = list(color = "#2C7FB8"), name = "In Median Cluster") %>%
        plotly::add_trace(x = ctr, y = yy$bg, type = "scatter", mode = "lines+markers",
                          line = list(color = "#9E9E9E"), marker = list(color = "#9E9E9E"), name = "Unclustered")
      if (normalized && is.finite(yy$ks))
        ttl <- list(text = sprintf("Detection-timing discrepancy (KS gap) = %.3f", yy$ks), font = list(size = 13))
    } else {
      gs <- unique(res$dataset); ng <- length(gs)
      pal <- series_palette(gs)
      for (g in gs) {
        sub    <- res[res$dataset == g, , drop = FALSE]
        cnt_c2 <- tabulate(sub$bin[sub$in_c2],  nbins = bins)
        cnt_bg <- tabulate(sub$bin[!sub$in_c2], nbins = bins)
        yy <- yvals(cnt_c2, cnt_bg); ytitle <- yy$yt
        kslab <- if (normalized && is.finite(yy$ks)) sprintf(" (KS %.2f)", yy$ks) else ""
        fig <- fig %>%
          plotly::add_trace(x = ctr, y = yy$c2, type = "scatter", mode = "lines",
                            line = list(color = pal[g], width = 2),
                            name = paste0(g, " - In Median Cluster", kslab)) %>%
          plotly::add_trace(x = ctr, y = yy$bg, type = "scatter", mode = "lines",
                            line = list(color = pal[g], width = 2, dash = "dash"),
                            name = paste0(g, " - Unclustered"))
      }
    }
    fig %>% plotly::layout(
      title  = ttl,
      xaxis  = list(title = "Acquisition time (tim)"),
      yaxis  = list(title = ytitle),
      legend = list(orientation = "h")
    )
  }
  output$det_time_plot <- renderPlotly(build_det_time_plot_fig())

  output$download_det_time <- downloadHandler(
    filename = function()
      paste0("cluster2_detection_over_time_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- det_time_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # Trace-ID detection over time - unique traces vs acquisition time
  # ==========================================================================
  # GOAL: show how many UNIQUE TRACE IDs are detected as a function of
  # acquisition time. A trace ID is one continuous detection of a blinking
  # fluorophore, so this is the detection-event RATE over the run - a direct QC
  # of whether imaging stayed stable (a declining rate flags photobleaching,
  # exhausted blinking, drift, or label/sample depletion).
  #
  # DIFFERS from the Median-Cluster detection plot above: no clustering is
  # required (works on raw traces), and it draws ONE SERIES PER FILE - or per
  # group in group mode - rather than a clustered/unclustered split.
  #
  # Each trace is counted ONCE, at its median tim, so a long trace with many
  # localizations does not dominate the early bins. Per-localization metadata
  # filters apply (via apply_meta_filters_batch); the trace-length filter is
  # deliberately NOT applied, because short blinking events are real detections.
  #
  # Each dataset has its OWN tim timeline, so each series is drawn only across
  # its own acquisition window (see the plot): the line STOPS where that
  # dataset's tim ends, rather than dropping to a misleading zero.
  #
  # RESULT stored: data.frame(tid, tim_med, dataset) - one row per trace.
  tid_time_results <- reactiveVal(NULL)

  observeEvent(input$run_tid_time, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    # Active-file-only unless "all files" is ticked.
    if (!isTRUE(input$tid_time_all_files)) {
      af <- input$active_file
      if (!is.null(af) && af %in% names(store)) store <- store[af]
    }

    # Trace-ID length prefilter (localizations per trace), the SAME range applied
    # to every dataset. len_max = 0 means "no upper limit". Read once here.
    len_min <- max(1L, suppressWarnings(as.integer(input$tid_time_len_min %||% 1L)))
    len_max <- suppressWarnings(as.integer(input$tid_time_len_max %||% 0L))
    if (is.na(len_min)) len_min <- 1L
    if (is.na(len_max)) len_max <- 0L

    all_tid <- list()
    withProgress(message = "Computing trace detection over time...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- qc_exclude_from(store[[fname]], fname)
        if (!is.data.frame(df) || !nrow(df)) next
        # Row-level metadata filters only (efo/cfr/dcr/...); NOT spread/length.
        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df))) next
        df <- apply_density_thinning(df, fname)   # opt-in labelling-density match
        if (is.null(df) || !nrow(df)) next
        # Trace-length prefilter: keep only trace IDs whose localization count is
        # in [len_min, len_max]. Counts are per-file (a tid is within one file).
        tl <- table(df$tid)
        keep_tids <- names(tl)[tl >= len_min & (len_max <= 0L | tl <= len_max)]
        df <- df[df$tid %in% keep_tids, , drop = FALSE]
        if (!nrow(df)) next
        # One row per trace, placed at its median acquisition time.
        per_tid <- df %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(tim_med = stats::median(tim, na.rm = TRUE), .groups = "drop") %>%
          dplyr::filter(is.finite(tim_med)) %>%
          dplyr::mutate(dataset = fname)
        if (nrow(per_tid)) all_tid[[fname]] <- per_tid
      }
    })
    res <- do.call(rbind, Filter(Negate(is.null), all_tid))
    if (is.null(res) || !nrow(res)) {
      showNotification("No traces with a tim column found. Check the metadata filters.",
                       type = "warning")
      tid_time_results(NULL); return()
    }
    tid_time_results(res)
    showNotification(
      sprintf("Trace detection over time: %d traces across %d file(s).",
              nrow(res), length(unique(res$dataset))),
      type = "message", duration = 4)
  })

  build_tid_time_plot_fig <- function() {
    res <- tid_time_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Press 'Run' to compute trace detection over time."))
    # Relabel file -> group when in group mode (drops unassigned files).
    res <- apply_dataset_groups(res)
    shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))

    bins <- max(5L, as.integer(input$tid_time_bins %||% 40L))
    mode <- input$tid_time_yaxis %||% "count"

    # Shared time axis (union of all datasets) so files/groups line up on x.
    rng <- range(res$tim_med, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 0) rng <- c(rng[1], rng[1] + 1)
    breaks  <- seq(rng[1], rng[2], length.out = bins + 1L)
    ctr     <- (breaks[-1] + breaks[-(bins + 1L)]) / 2
    res$bin <- pmin(pmax(findInterval(res$tim_med, breaks, rightmost.closed = TRUE), 1L), bins)

    # Per-series y-values honouring the selected mode. `active` marks the bins
    # that fall inside THIS dataset's own tim window; bins outside it are set to
    # NA so the line STOPS at the dataset's end instead of dropping to zero
    # (zero would wrongly read as "detected nothing" rather than "acquisition
    # over"). A genuine empty bin WITHIN the window stays 0 - a real dip.
    yvals <- function(cnt, active) {
      tot <- sum(cnt); sc <- if (tot > 0) tot else 1
      y <- switch(mode,
        cumulative      = cumsum(cnt),
        norm            = cnt / sc,
        norm_cumulative = cumsum(cnt) / sc,
                          cnt)
      y[!active] <- NA_real_
      yt <- switch(mode,
        cumulative      = "Cumulative unique traces",
        norm            = "Detection fraction per bin",
        norm_cumulative = "Cumulative detection fraction",
                          "Unique traces per bin")
      list(y = y, yt = yt)
    }

    series <- unique(res$dataset)
    pal    <- series_palette(series)
    fig    <- plotly::plot_ly()
    ytitle <- "Unique traces per bin"
    for (s in series) {
      tm  <- res$tim_med[res$dataset == s]
      cnt <- tabulate(res$bin[res$dataset == s], nbins = bins)
      # Bins spanned by this dataset's own [min, max] tim.
      b_lo   <- pmin(pmax(findInterval(min(tm), breaks, rightmost.closed = TRUE), 1L), bins)
      b_hi   <- pmin(pmax(findInterval(max(tm), breaks, rightmost.closed = TRUE), 1L), bins)
      active <- seq_len(bins) >= b_lo & seq_len(bins) <= b_hi
      yy  <- yvals(cnt, active); ytitle <- yy$yt
      fig <- fig %>% plotly::add_trace(
        x = ctr, y = yy$y, type = "scatter", mode = "lines+markers",
        line = list(color = pal[s], width = 2), marker = list(color = pal[s], size = 5),
        connectgaps = FALSE, name = s)
    }
    fig %>% plotly::layout(
      xaxis  = list(title = "Acquisition time (tim)"),
      yaxis  = list(title = ytitle),
      legend = list(orientation = "h")
    )
  }
  output$tid_time_plot <- renderPlotly(build_tid_time_plot_fig())

  output$download_tid_time <- downloadHandler(
    filename = function()
      paste0("trace_detection_over_time_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- tid_time_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # QC 1 - Localization precision cloud (within-trace scatter about the centroid)
  # ==========================================================================
  # GOAL: the direct measure of localization PRECISION. Pool every localization
  # minus its OWN trace centroid, recentre at (0,0), and look at the scatter.
  # Its RMS radius is σ_loc. The SHAPE flags problems: elongation = drift within a
  # trace or direction-dependent estimator precision; heavy tails = mislinks /
  # double events; skew = a systematic bias. Uses raw uploaded data (before
  # metadata filters), so it reports the precision as acquired.
  #
  # NOT a PSF width. MINFLUX fits no PSF - the doughnut is the EXCITATION
  # pattern and position comes from photon-count ratios about its central zero,
  # so this scatter is not a spot width and does not scale like one.
  prec_cloud_results <- reactiveVal(NULL)
  observeEvent(input$run_prec_cloud, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    # Precision is ALWAYS assessed on the single ACTIVE dataset - never pooled across
    # files (an aggregate σ would blur per-acquisition precision). Restrict to it here.
    af <- input$active_file
    if (is.null(af) || !af %in% names(store)) af <- names(store)[[1]]
    store <- store[af]
    min_loc <- max(3L, as.integer(input$prec_cloud_min_loc %||% 5L))
    use_z   <- isTRUE(input$prec_cloud_use_z)
    dxall <- dyall <- dzall <- numeric(0)          # pooled residuals (length-weighted)
    # Per-TRACE VARIANCES (one vote each, robust) and per-TRACE covariances. The
    # cross terms are needed because a cloud tilted off the x/y axes has its
    # elongation entirely there - per-axis spreads alone cannot see it.
    vxt <- vyt <- vzt <- numeric(0)
    cxyt <- cxzt <- cyzt <- numeric(0)
    nlt <- integer(0)                              # localizations per trace
    # Pooled sums of squares + their exact degrees of freedom. Kept separately
    # from the residual vectors because the correct divisor is sum(n_i - 1),
    # NOT the number of localizations: subtracting each trace's own centre
    # spends one degree of freedom PER TRACE. Dividing by N - 1 instead (what
    # sd() on the pooled residuals does) reads ~5% low at 5 localizations per
    # trace and ~1% low at 15 - a real, n-dependent underestimate of sigma.
    ssx_all <- ssy_all <- ssz_all <- 0
    dof_all <- dofz_all <- 0L
    per_file <- list()
    pooled   <- list()                             # per-loc residuals (dx/dy/dz) for the plot
    # Median-of-chi-square correction. s^2 ~ sigma^2 * chisq(n-1)/(n-1), which is
    # right-skewed, so the MEDIAN of the per-trace variances sits BELOW sigma^2 -
    # by 16% at n = 5, i.e. sigma reads 8% low. Dividing by that distribution's
    # own median removes it exactly (same idea as the 1.4826 that makes the MAD
    # consistent for sigma). NOTE the usual c4 factor is the WRONG correction
    # here: c4 corrects E[s], and this estimator takes a median, not a mean.
    med_fac <- function(n) ifelse(n > 1L, stats::qchisq(0.5, n - 1L) / (n - 1L), NA_real_)
    withProgress(message = "Computing precision cloud...", value = 0, {
      nf <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]; incProgress(1 / nf, detail = fname)
        # Assess precision on the SAME localizations the Display tab plots: apply the
        # metadata / trace-length / trace-spread filters (and any fiducial exclusion),
        # so the cloud reflects the post-filter data rather than the raw acquisition.
        df <- tryCatch(apply_meta_filters_batch(qc_exclude_from(store[[fname]], fname), input),
                       error = function(e) store[[fname]])
        if (!is.data.frame(df) || !nrow(df) || !all(c("tid","x","y") %in% names(df))) next
        hasz <- use_z && has_axial_data(df)   # z must vary; a planar z = 0 column is not 3D
        # Residuals are taken about each trace's MEAN, not its median. The mean is
        # the maximum-likelihood centre under the Gaussian error model this panel
        # assumes, and it is the only centre for which the degree-of-freedom
        # corrections below are exact - a median centre has no matching
        # correction. It is also the honest choice: a robust centre partly
        # ABSORBS the mislinked localizations that this panel exists to reveal.
        # Outlier resistance is supplied instead by the per-trace median further
        # down, which gives every trace one vote regardless of length.
        dd <- df %>%
          dplyr::group_by(tid) %>%
          dplyr::mutate(nloc = dplyr::n()) %>%
          dplyr::filter(nloc >= min_loc) %>%
          dplyr::mutate(dx = x - mean(x, na.rm = TRUE),
                        dy = y - mean(y, na.rm = TRUE),
                        dz = if (hasz) z - mean(z, na.rm = TRUE) else NA_real_) %>%
          dplyr::ungroup()
        if (!nrow(dd)) next
        pooled[[fname]] <- data.frame(dx = dd$dx, dy = dd$dy,
                                      dz = if (hasz) dd$dz else NA_real_)
        dxall <- c(dxall, dd$dx); dyall <- c(dyall, dd$dy)
        if (hasz) dzall <- c(dzall, dd$dz[is.finite(dd$dz)])
        # One robust σ per trace (each trace counts once, so a 2000-loc fiducial
        # doesn't dominate the precision estimate the way the pooled cloud does).
        safe_cov <- function(u, v) {
          ok <- is.finite(u) & is.finite(v)
          if (sum(ok) < 3L) return(NA_real_)
          stats::cov(u[ok], v[ok])
        }
        pts <- dd %>% dplyr::group_by(tid) %>%
          dplyr::summarise(nl  = dplyr::n(),
                           ssx = sum(dx^2, na.rm = TRUE),
                           ssy = sum(dy^2, na.rm = TRUE),
                           ssz = if (hasz) sum(dz^2, na.rm = TRUE) else NA_real_,
                           vx  = stats::var(dx, na.rm = TRUE),
                           vy  = stats::var(dy, na.rm = TRUE),
                           vz  = if (hasz) stats::var(dz, na.rm = TRUE) else NA_real_,
                           cxy = safe_cov(dx, dy),
                           cxz = if (hasz) safe_cov(dx, dz) else NA_real_,
                           cyz = if (hasz) safe_cov(dy, dz) else NA_real_,
                           .groups = "drop")
        # ONE factor applied to every second moment of a trace, so the ellipse
        # SHAPE (and hence the anisotropy and tilt) is untouched - only its size
        # is de-biased.
        mf <- med_fac(pts$nl)
        vxt <- c(vxt, pts$vx / mf); vyt <- c(vyt, pts$vy / mf)
        # Keep vzt the SAME LENGTH as vxt/vyt/cxzt by carrying NA through rather
        # than dropping non-finite entries here. robust_median_cov() drops them
        # itself; filtering only this one vector meant var_z and cov_xz were
        # medians over different subsets of traces, so the robust xz/yz matrices
        # mixed two populations.
        if (hasz) vzt <- c(vzt, pts$vz / mf)
        cxyt <- c(cxyt, pts$cxy / mf)
        if (hasz) { cxzt <- c(cxzt, pts$cxz / mf); cyzt <- c(cyzt, pts$cyz / mf) }
        nlt <- c(nlt, pts$nl)
        # Pooled sums of squares with their exact degrees of freedom.
        dof_f  <- sum(pts$nl - 1L)
        dofz_f <- sum((pts$nl - 1L)[is.finite(pts$ssz)])
        ssx_all <- ssx_all + sum(pts$ssx, na.rm = TRUE)
        ssy_all <- ssy_all + sum(pts$ssy, na.rm = TRUE)
        ssz_all <- ssz_all + sum(pts$ssz, na.rm = TRUE)
        dof_all <- dof_all + dof_f; dofz_all <- dofz_all + dofz_f
        per_file[[fname]] <- data.frame(
          dataset  = fname,
          sigma_x  = sqrt(sum(pts$ssx, na.rm = TRUE) / max(dof_f, 1L)),
          sigma_y  = sqrt(sum(pts$ssy, na.rm = TRUE) / max(dof_f, 1L)),
          sigma_z  = if (dofz_f > 0) sqrt(sum(pts$ssz, na.rm = TRUE) / dofz_f) else NA_real_,
          n_traces = nrow(pts), n_locs = nrow(dd),
          stringsAsFactors = FALSE)
      }
    })
    if (!length(dxall)) {
      showNotification("No traces with enough localizations (lower 'Min localizations / trace').",
                       type = "warning"); prec_cloud_results(NULL); return()
    }
    mom <- function(v, p) { v <- v[is.finite(v)]; m <- mean(v); s <- stats::sd(v)
                            if (s <= 0) NA_real_ else mean((v - m)^p) / s^p }
    # Pooled sigma, exact degrees of freedom (see the accumulator comment above).
    sx <- sqrt(ssx_all / max(dof_all, 1L)); sy <- sqrt(ssy_all / max(dof_all, 1L))
    sz <- if (dofz_all > 0) sqrt(ssz_all / dofz_all) else NA_real_
    sx_med <- sqrt(stats::median(vxt, na.rm = TRUE))
    sy_med <- sqrt(stats::median(vyt, na.rm = TRUE))
    sz_med <- if (length(vzt)) sqrt(stats::median(vzt, na.rm = TRUE)) else NA_real_
    # subsample rows for the plot (dx/dy/dz/is_fid kept aligned); σ computed from all
    P <- do.call(rbind, pooled)

    # ── Covariance matrices: the shape of the cloud, tilt included ────────────
    # POOLED describes exactly the population the heatmap draws (every
    # localization, so long traces weigh more). ROBUST gives each trace one vote.
    # Both are kept and both are drawn: when they disagree the gap IS the
    # "a few long/wide traces dominate" diagnostic, made visible.
    pair_cov <- function(u, v) {
      ok <- is.finite(u) & is.finite(v)
      if (sum(ok) < 3L) return(NULL)
      stats::cov(cbind(u[ok], v[ok]))
    }
    # cov() divides by (rows - 1); the honest divisor is sum(n_i - 1). Rescaling
    # by that ratio keeps the drawn pooled ellipse consistent with the pooled
    # sigma reported above, instead of ~5% tighter at short traces.
    dof_scale <- if (dof_all > 0) (nrow(P) - 1) / dof_all else 1
    cov_pooled <- lapply(list(
      xy = pair_cov(P$dx, P$dy),
      xz = if (length(dzall)) pair_cov(P$dx, P$dz) else NULL,
      yz = if (length(dzall)) pair_cov(P$dy, P$dz) else NULL),
      function(S) if (is.null(S)) NULL else S * dof_scale)
    cov_trace <- list(
      xy = robust_median_cov(vxt, vyt, cxyt),
      xz = if (length(vzt)) robust_median_cov(vxt, vzt, cxzt) else NULL,
      yz = if (length(vzt)) robust_median_cov(vyt, vzt, cyzt) else NULL)

    # Headline anisotropy comes from the ROBUST matrix and is the eigenvalue
    # ratio, not max(σx,σy)/min(σx,σy) - the latter reads ~1.0 on a 45°-tilted
    # cloud and would report drifted data as "isotropic".
    ell_r <- cov_ellipse_params(cov_trace$xy)

    # ── How heavy SHOULD the tails be? ────────────────────────────────────────
    # The pooled cloud mixes traces that do NOT share one sigma: MINFLUX
    # precision scales as 1/sqrt(photons) and photon counts vary a lot between
    # traces. A mixture of Gaussians of DIFFERENT widths is leptokurtic no
    # matter how clean the data are - for X|v ~ N(0, v), kurtosis =
    # 3·E[v²]/E[v]². So comparing the observed kurtosis against 3 (which this
    # panel used to do, warning "possible mislinks / double events") flags
    # perfectly good data: at a realistic spread of per-trace sigma the
    # kurtosis of artifact-free simulated data reaches 7-15.
    # What is diagnostic is the EXCESS over what the observed sigma spread
    # already explains. Traces are weighted by n_i because that is how they
    # enter the pooled cloud. The per-trace variances are themselves noisy -
    # v̂ ~ v·χ²(n-1)/(n-1) gives E[v̂²] = v²(1 + 2/(n-1)) - so E[v²] is recovered
    # by the EXACT factor (n-1)/(n+1). Without it, short traces manufacture a
    # heterogeneity that is not there, the prediction comes out too low, and
    # the panel goes back to crying artifact on clean data.
    vv <- (vxt + vyt) / 2                          # per-trace variance, de-biased
    ok_v <- is.finite(vv) & vv > 0 & is.finite(nlt)
    kurt_exp <- NA_real_; sigma_spread <- NA_real_
    if (sum(ok_v) >= 5L) {
      v <- vv[ok_v]; nl <- nlt[ok_v]; w <- nl / sum(nl)
      Ev  <- sum(w * v)
      Ev2 <- sum(w * v^2 * (nl - 1L) / (nl + 1L))
      kurt_exp <- 3 * max(Ev2, Ev^2) / Ev^2
      # Spread of the per-trace precision itself - a real property of the
      # acquisition, reported in its own right rather than mislabelled as tails.
      sigma_spread <- stats::IQR(sqrt(v)) / stats::median(sqrt(v))
    }

    # ── Sample stability: drift and jitter ────────────────────────────────────
    # Belongs here rather than in the de-blink panel: it describes the SAMPLE,
    # and folding it into the merge radius would be circular (the only molecules
    # whose motion is observable are the ones already merged, i.e. those that
    # moved LESS than the radius in question).
    # Measured per file and pooled, because time is only comparable within one
    # acquisition - two files start their clocks independently.
    drift_by_file <- lapply(names(store), function(fname) {
      dfx <- tryCatch(apply_meta_filters_batch(qc_exclude_from(store[[fname]], fname), input),
                      error = function(e) NULL)
      if (!is.data.frame(dfx) || !nrow(dfx) || !"tim" %in% names(dfx)) return(NULL)
      # Group traces into molecules first, so a fluorophore seen at several
      # separate times contributes a trajectory rather than one displacement.
      # A generous radius is used deliberately: pairs that are really two
      # different molecules displace in no particular direction, so they widen
      # the interval without biasing the drift.
      sgq <- tryCatch(suggest_centroid_threshold(dfx, min_n = min_loc)$sigma,
                      error = function(e) NULL)
      gid <- if (!is.null(sgq) && is.finite(sgq))
        tryCatch(deblink_merge_centroids(dfx, 5 * sgq, NULL), error = function(e) NULL) else NULL
      tryCatch(sample_drift_estimate(dfx, min_n = min_loc, group = gid),
               error = function(e) NULL)
    })
    names(drift_by_file) <- names(store)
    drift_by_file <- Filter(Negate(is.null), drift_by_file)

    keep <- if (nrow(P) > 40000L) sort(sample.int(nrow(P), 40000L)) else seq_len(nrow(P))
    Pp <- P[keep, , drop = FALSE]
    prec_cloud_results(list(
      drift = drift_by_file,
      dx = Pp$dx, dy = Pp$dy, dz = Pp$dz,
      sigma_x = sx, sigma_y = sy, sigma_z = sz,
      sigma_x_med = sx_med, sigma_y_med = sy_med, sigma_z_med = sz_med,
      sigma_xy   = sqrt((sx^2 + sy^2) / 2),
      # robust, one-vote-per-trace precision (immune to fiducials / long traces)
      sigma_xy_med = sqrt((sx_med^2 + sy_med^2) / 2),
      cov_pooled = cov_pooled, cov_trace = cov_trace,
      # Rotation-invariant shape of the robust cloud, and WHICH WAY it leans.
      aniso_ev  = if (!is.null(ell_r)) ell_r$ratio     else NA_real_,
      tilt_deg  = if (!is.null(ell_r)) ell_r$angle_deg else NA_real_,
      rho_xy    = if (!is.null(ell_r)) ell_r$rho       else NA_real_,
      # Kept for the per-axis read-out; NOT used for the warning any more, since
      # it cannot see elongation that runs diagonally.
      anisotropy = max(sx, sy) / max(1e-9, min(sx, sy)),
      skew_x = mom(dxall, 3), skew_y = mom(dyall, 3),
      kurt_x = mom(dxall, 4), kurt_y = mom(dyall, 4),
      # What the per-trace sigma spread alone predicts, and that spread itself.
      kurt_exp = kurt_exp, sigma_spread = sigma_spread,
      n_locs = length(dxall), n_traces = length(vxt), has_z = length(dzall) > 0,
      per_file = do.call(rbind, per_file)))
    showNotification(sprintf("Precision cloud: %s localizations; σ_xy = %.2f nm.",
                             format(length(dxall), big.mark = ","), sqrt((sx^2 + sy^2) / 2)),
                     type = "message", duration = 4)
  })

  output$prec_cloud_plot <- renderPlotly({
    r <- prec_cloud_results()
    shiny::validate(shiny::need(!is.null(r), "Press 'Run' to compute the precision cloud."))

    # Each axis pair gets its own SQUARE panel: density heatmap + 1σ/2σ ellipses, a
    # faint zero crosshair through (0,0), no background grid, and a shared ±L nm frame.
    # The panels are built manually on explicit independent axes (x/y, x2/y2, x3/y3)
    # rather than plotly::subplot() - subplot does not reliably re-map per-panel
    # scaleanchor, which cross-linked the zoom across panels. Here each y_i is scale-
    # anchored to its own x_i, so aspect is equal (1 nm x == 1 nm y) AND zooming one
    # panel leaves the others untouched. Ellipse semi-axes are the per-axis robust σ.
    # Colour encodes the PROBABILITY MASS, line style encodes data-vs-model:
    # solid = where that fraction of the localizations actually is, dotted =
    # where it would be if the cloud were Gaussian. The two coincide on clean
    # data and separate where the model fails, which is the entire diagnostic.
    c_50 <- "#c03030"; c_90 <- "#2c7fb8"   # 50% red, 90% blue
    c_rb <- "#1a9850"                      # per-trace robust fit = green
    # TWO KINDS OF RING, EACH LABELLED AS WHAT IT IS
    #   SOLID  = where the localizations actually are: the contours enclosing
    #            50% and 90% of them. Round masses because a reader can act on
    #            them.
    #   DOTTED = the fitted Gaussian at 1 and 2 sigma. Semi-axes come from the
    #            eigen-decomposition of the covariance, so each axis gets its own
    #            sigma and the ring is a rotated ELLIPSE, not a circle - which is
    #            what makes its tilt readable as drift.
    #
    #   The two families are NOT the same probability mass: a 2D Gaussian holds
    #   1 - exp(-k^2/2) inside k sigma, i.e. 39.3% at 1 sigma and 86.5% at 2. So
    #   on a perfectly Gaussian cloud the solid 50% ring sits ~18% outside the
    #   dotted 1 sigma ring, and the solid 90% ~7% outside the dotted 2 sigma.
    #   Judge SHAPE and TILT by these rings, not a millimetre-exact overlap; the
    #   tail comparison is carried numerically by the kurtosis in the stats line.
    P_IN <- c(0.50, 0.90)
    # Legend names, also used as legendgroup keys - keep them short and literal.
    NM_50  <- "50 % of locs";  NM_1S <- "1σ (Gaussian fit)"
    NM_90  <- "90 % of locs";  NM_2S <- "2σ (Gaussian fit)"
    NM_TR  <- "1σ per trace"
    has_z <- isTRUE(r$has_z) && any(is.finite(r$dz))
    # Sp = pooled covariance (the population the heatmap shows)
    # Sr = robust per-trace covariance (the headline σ)
    panels <- if (has_z) list(
      list(a = r$dx, b = r$dy, Sp = r$cov_pooled$xy, Sr = r$cov_trace$xy, la = "Δx (nm)", lb = "Δy (nm)"),
      list(a = r$dx, b = r$dz, Sp = r$cov_pooled$xz, Sr = r$cov_trace$xz, la = "Δx (nm)", lb = "Δz (nm)"),
      list(a = r$dy, b = r$dz, Sp = r$cov_pooled$yz, Sr = r$cov_trace$yz, la = "Δy (nm)", lb = "Δz (nm)")
    ) else list(
      list(a = r$dx, b = r$dy, Sp = r$cov_pooled$xy, Sr = r$cov_trace$xy,
           la = "Δx from trace centroid (nm)", lb = "Δy from trace centroid (nm)"))
    # ── Frame half-width, derived from the data ───────────────────────────────
    # This used to be a hard-coded 30 nm, which only frames the plot correctly
    # when σ happens to be ~10 nm. On any dataset with a wider cloud the heatmap
    # filled the whole panel and the 1σ/2σ rings ran clean off it, leaving stray
    # arcs cutting across the corners with no visible centre - the plot looked
    # broken and over-zoomed because it WAS over-zoomed.
    # The frame follows the CORE - the robust, one-vote-per-trace σ - at a plain
    # 3σ. It must not be driven by the POOLED covariance: σ varies between
    # traces (photon count), so the pooled ellipse is pulled outward by a
    # minority of wide traces and sizing the frame to it leaves the visible
    # cloud small in the middle with only the zero crosshair reaching the edges.
    # A clean Gaussian therefore lands at exactly ±3σ and fills the panel; the
    # 50%/90% contours sit at 1.18σ/2.15σ, comfortably inside.
    # The second term is the safety valve: the radius that actually contains the
    # outer mass, so a genuinely wide or heavy-tailed cloud still gets its outer
    # contour inside the frame instead of having it
    # silently dropped. It is a quantile, not a moment, so outliers cannot run
    # away with it. ONE L is shared by all panels so they stay comparable.
    # Each panel is framed on ITS OWN pair, not on a single shared L. In 3D
    # σ_z is routinely 2-3x σ_xy, and one shared frame is then sized by z -
    # which leaves the Δx-Δy cloud as a small dot in the middle of an empty
    # panel. The axis numbers, which are printed per panel, carry the
    # cross-panel comparison instead.
    # NOTE no rounding to "nice" numbers: plotly chooses tick POSITIONS inside
    # whatever range it is given, so rounding the range buys nothing and only
    # overshoots (36 -> 40 pushes the cloud back toward the middle third).
    ring_reach <- function(S, k) {          # semi-MAJOR bounds the ring on both axes
      e <- cov_ellipse_params(S, k = k)
      if (is.null(e)) NA_real_ else e$ra
    }
    panel_lim <- function(P) {
      need <- c(3.0 * ring_reach(P$Sr, 1),
                1.15 * stats::quantile(sqrt(P$a^2 + P$b^2), P_IN[2],
                                       na.rm = TRUE, names = FALSE))
      need <- need[is.finite(need)]
      core <- max(if (length(need)) signif(max(need), 3) else 30, 5)  # never over-zoom
      # The 2σ model ring may reach past a frame sized to the CORE, and when it
      # does it was being dropped and the panel silently lost a legend entry.
      # Let the frame grow to admit it, but CAP that growth: sizing the frame to
      # the pooled fit without a cap is what once collapsed the cloud to a dot
      # and washed the heatmap out, because a minority of wide traces can
      # inflate the pooled covariance without limit. Beyond the cap the ring is
      # still dropped - a ring that far out cannot be shown without wrecking the
      # view of the data it is meant to describe.
      r2 <- ring_reach(P$Sp, 2)
      if (is.finite(r2) && r2 > core) core <- min(1.05 * r2, 1.4 * core)
      core
    }
    Lv <- vapply(panels, panel_lim, numeric(1))

    np  <- length(panels)
    gap <- if (np > 1) 0.06 else 0          # horizontal gap between panels (paper frac)
    w   <- (1 - gap * (np - 1)) / np        # each panel's width (paper frac)

    fig <- plotly::plot_ly()
    shapes <- list(); axis_args <- list()
    for (i in seq_len(np)) {
      P <- panels[[i]]
      L <- Lv[i]                            # this panel's own half-range
      xid <- if (i == 1) "x" else paste0("x", i)
      yid <- if (i == 1) "y" else paste0("y", i)
      keep <- is.finite(P$a) & is.finite(P$b)
      # ── Density: kernel estimate, not a binned histogram ───────────────────
      # A 2D histogram estimates every bin independently, so each carries its
      # own full Poisson noise: at 60x60 bins and ~40k points that is ~9% at the
      # peak but ~26% at the 2σ radius and ~90% at 3σ. The outer contours - the
      # ones the eye reads as "the shape of the cloud" - were therefore mostly
      # shot noise, which is exactly the ragged, under-smoothed look. Re-binning
      # cannot fix that (coarser blurs the core, finer worsens the tails); a
      # kernel estimate can, because each point informs a neighbourhood.
      # zsmooth = "best" then interpolates the grid so no cell edges show.
      s_rob <- if (!is.null(P$Sr)) sqrt(mean(diag(P$Sr))) else NA_real_
      dh <- density_hdr_2d(P$a[keep], P$b[keep], lim = L, sigma = s_rob, probs = P_IN)
      if (!is.null(dh)) {
        fig <- fig %>% plotly::add_heatmap(
          x = dh$x, y = dh$y, z = t(dh$z),      # kde2d is z[x, y]; plotly wants z[y, x]
          xaxis = xid, yaxis = yid, colorscale = "Hot", reversescale = TRUE,
          showscale = FALSE, zsmooth = "best", hoverinfo = "skip", showlegend = FALSE)
      }
      # Ellipses are drawn from explicit outline coordinates: plotly's "circle"
      # shape is defined by a bounding box and can only ever be axis-aligned,
      # which made a diagonally elongated cloud look circular.
      # Rings are TRACES, not layout shapes, so the legend can switch them off.
      # legendgroup ties the three panels together (one click hides a ring
      # everywhere) and only panel 1 contributes the legend entry.
      # add_PATHS, never add_lines: add_lines() SORTS the points by x before
      # connecting them. These are closed curves, so sorting turns every ring
      # into a zigzag that fills its own interior - the rings render as solid
      # blobs and the dashed ones as scattered stipple. add_paths() keeps the
      # order the points were traced in.
      ring <- function(x, y, colr, dsh, nm) {
        fig <<- plotly::add_paths(fig, x = x, y = y, xaxis = xid, yaxis = yid,
          name = nm, legendgroup = nm, showlegend = (i == 1L),
          hoverinfo = "skip", inherit = FALSE,
          line = c(list(color = colr, width = 2), if (nzchar(dsh)) list(dash = dsh)))
      }
      # A ring whose major axis leaves the frame is DROPPED, not drawn. Drawn, it
      # appears only as arcs cutting the corners with no visible centre, which
      # reads as a broken plot; and the frame must stay sized to the core, so
      # growing it to fit is not the alternative. The stats line under the plot
      # already names this case ("pooled σ >> per-trace σ").
      ell <- function(S, k, colr, dsh, nm) {
        e <- cov_ellipse_params(S, k = k)
        if (is.null(e) || !is.finite(e$ra) || e$ra > L) return(invisible(NULL))
        xy <- ellipse_xy(e$ra, e$rb, e$theta)
        if (is.null(xy)) return(invisible(NULL))
        ring(xy$x, xy$y, colr, dsh, nm)
      }
      axl <- function(x0, x1, y0, y1) list(type = "line", xref = xid, yref = yid,
        x0 = x0, x1 = x1, y0 = y0, y1 = y1, line = list(color = "rgba(0,0,0,0.4)", width = 1))
      # Empirical highest-density contours, drawn SOLID at the same two
      # probability masses the dotted Gaussian rings enclose. Same mass, same
      # colour, different line style: wherever a solid line sits outside its
      # dotted partner, the cloud has more localizations out there than a
      # Gaussian of the fitted width would put there. A count-level contour
      # could not be compared to anything.
      # Match each empirical contour to its mass by VALUE, not by a hard-coded
      # midpoint - with the inner mass at 0.50 a "prob < 0.5" test would colour
      # it as the outer one. A mass can yield several disjoint contours; they
      # share one legendgroup so they toggle together, and only the first
      # carries the legend entry.
      seen <- character(0)
      for (q in dh$paths %||% list()) {
        inner <- q$prob <= P_IN[1] + 1e-9
        nm <- if (inner) NM_50 else NM_90
        show <- (i == 1L) && !(nm %in% seen); seen <- c(seen, nm)
        fig <- plotly::add_paths(fig, x = q$x, y = q$y, xaxis = xid, yaxis = yid,
          name = nm, legendgroup = nm, showlegend = show, hoverinfo = "skip",
          inherit = FALSE,
          line = list(color = if (inner) c_50 else c_90, width = 2))
      }
      ell(P$Sp, 1, c_50, "dot",  NM_1S)   # 1σ ellipse, pooled fit
      ell(P$Sp, 2, c_90, "dot",  NM_2S)   # 2σ ellipse, pooled fit
      # Same 1σ, but fitted with one vote per trace instead of one per
      # localization: where the green sits well inside the dotted red, the
      # pooled width is being carried by a minority of wide or long traces.
      ell(P$Sr, 1, c_rb, "dash", NM_TR)
      shapes <- c(shapes, list(axl(0, 0, -L, L), axl(-L, L, 0, 0)))
      dom <- c((i - 1) * (w + gap), (i - 1) * (w + gap) + w)
      # Anchor the X axis to Y and let X give up DOMAIN (width) to satisfy the
      # 1:1 scale. Anchoring the other way round let plotly satisfy the ratio by
      # widening X's RANGE instead, which produced panels spanning x -20..40
      # against y -20..20: asymmetric, off-centre, and 1.5x wider than tall. An
      # ellipse read against axes like that looks rotated even when it is not.
      # Both axes now show exactly -L..L, so the data window really is square
      # and a tilt on screen is a tilt in the data.
      axis_args[[if (i == 1) "xaxis" else paste0("xaxis", i)]] <- list(
        title = P$la, range = c(-L, L), zeroline = FALSE, showgrid = FALSE,
        domain = dom, anchor = yid, scaleanchor = yid, scaleratio = 1,
        constrain = "domain")
      # Y carries the range and the full height; X is the axis that yields (above).
      # Only ONE of the pair may declare scaleanchor - declaring it on both makes
      # the constraint circular and plotly resolves it by widening a range.
      axis_args[[if (i == 1) "yaxis" else paste0("yaxis", i)]] <- list(
        title = P$lb, range = c(-L, L), zeroline = FALSE, showgrid = FALSE,
        domain = c(0, 1), anchor = xid)
    }
    # The rings are real traces now, so the legend entries are the rings
    # themselves: clicking one hides it across every panel (they share a
    # legendgroup), and double-click isolates it. Solid is where the
    # localizations are, dotted is where a Gaussian of the fitted width would
    # put the same fraction, dashed gives every trace one vote instead of every
    # localization.
    p <- do.call(plotly::layout, c(list(fig, shapes = shapes,
      legend = list(orientation = "h", x = 0.5, xanchor = "center",
                    y = -0.24, yanchor = "top",
                    itemclick = "toggle", itemdoubleclick = "toggleothers"),
      margin = list(t = 10, b = 100)), axis_args))
    # Autoscale-to-data interacts with scaleanchor to blow the ranges up (±100) and
    # desync the panel sizes. Make double-click / reset return to the fixed ±L frame
    # instead, and drop the autoscale button so the only "fit" action is that reset -
    # the panels then stay identical squares. Box-zoom (per panel) still works.
    plotly::config(p, doubleClick = "reset", displaylogo = FALSE,
                   modeBarButtonsToRemove = list("autoScale2d"))
  })

  output$prec_cloud_stats <- renderUI({
    r <- prec_cloud_results(); if (is.null(r)) return(NULL)
    warn <- character(0)
    # Keyed off the EIGENVALUE ratio, not max(σx,σy)/min(σx,σy): the latter scores
    # ~1.0 on a 45°-tilted cloud and would call drifted data isotropic. The angle
    # is reported because drift has a direction and that is the useful clue.
    if (is.finite(r$aniso_ev) && r$aniso_ev > 1.3)
      warn <- c(warn, sprintf("elongated %.2f:1 along %.0f° (drift within traces, or direction-dependent estimator precision)",
                              r$aniso_ev, r$tilt_deg))
    # Tails are judged against what the per-trace σ SPREAD already predicts, not
    # against 3. A pooled cloud of Gaussians with different widths is
    # leptokurtic by construction, and MINFLUX σ varies with photon count, so
    # "kurtosis > 5" fires on clean data (7-15 on artifact-free simulations at a
    # realistic σ spread). Only the excess over the predicted value is evidence.
    # Both an absolute and a relative margin must be cleared. The relative one
    # carries the very-short-trace case: at 5 localizations per trace the
    # per-trace variances have 4 degrees of freedom, so the prediction is itself
    # noisy and a bare "excess > 2" still misfires on clean simulated data.
    kobs <- suppressWarnings(max(r$kurt_x, r$kurt_y, na.rm = TRUE))
    if (is.finite(kobs) && is.finite(r$kurt_exp) &&
        kobs - r$kurt_exp > 2 && kobs > 1.5 * r$kurt_exp)
      warn <- c(warn, sprintf("heavy tails beyond what the σ spread explains (kurtosis %.1f vs %.1f expected - possible mislinks / double events)",
                              kobs, r$kurt_exp))
    # If the length-weighted pooled σ is much larger than the per-trace median, a
    # few long/wide traces (fiducials!) dominate - flag it and point at the scan.
    if (is.finite(r$sigma_xy_med) && r$sigma_xy > 1.5 * r$sigma_xy_med)
      warn <- c(warn, sprintf("pooled σ (%.1f) ≫ per-trace σ (%.1f): the pooled number is carried by a minority of wide or long traces%s",
                              r$sigma_xy, r$sigma_xy_med,
                              # Both a fiducial AND ordinary photon-count spread
                              # produce this gap; naming only fiducials sent the
                              # user hunting for contamination that need not exist.
                              if (is.finite(r$sigma_spread) && r$sigma_spread > 0.5)
                                " - σ also varies widely between traces (photon-count spread), so check that before assuming contamination"
                              else " - run the Fiducial scan"))
    htmltools::div(style = "font-size:0.9em; margin-top:4px;",
      htmltools::tags$b(sprintf("σ_xy = %.2f nm (median per trace, robust)", r$sigma_xy_med)),
      if (r$has_z && is.finite(r$sigma_z_med)) sprintf(" · σ_z %.2f nm", r$sigma_z_med) else "",
      sprintf("  |  elongation %.2f:1 at %.0f° (ρ %.2f)  |  pooled σ_xy %.2f (σ_x %.2f · σ_y %.2f · axis ratio %.2f)  |  σ varies %.0f%% across traces (IQR/median)  |  skew x/y %.2f/%.2f  ·  kurtosis x/y %.1f/%.1f (%.1f expected from that σ spread)  |  %s locs, %s traces",
              r$aniso_ev, r$tilt_deg, r$rho_xy,
              r$sigma_xy, r$sigma_x, r$sigma_y, r$anisotropy,
              100 * (r$sigma_spread %||% NA_real_), r$skew_x, r$skew_y,
              r$kurt_x, r$kurt_y, r$kurt_exp %||% NA_real_,
              format(r$n_locs, big.mark = ","),
              format(r$n_traces, big.mark = ",")),
      if (length(warn)) htmltools::div(style = "color:#b26a00; margin-top:2px;",
        "⚠ ", paste(warn, collapse = "; "), ".") else
        htmltools::div(style = "color:#2a6; margin-top:2px;", "Cloud is compact and near-isotropic."),
      # ── Sample stability, one line per file ─────────────────────────────────
      # Drift accumulates with time and has a direction; jitter does not. They
      # are reported separately because they say different things about a
      # sample, and because jitter - not drift - is what sets how far apart two
      # blinks of one molecule land.
      if (length(r$drift %||% list())) {
        htmltools::tagList(
          htmltools::div(style = "margin-top:6px; font-weight:600;", "Sample stability"),
          lapply(names(r$drift), function(fn) {
            d <- r$drift[[fn]]
            txt <- sprintf("%s: drift %.1f nm over %.0f min (95%% CI %.1f-%.1f), %.1f nm across a typical %.0f s gap between blinks; jitter %.1f nm; %.0f%% of pair displacement is directional. From %s pairs%s.",
                           fn, d$drift_run, d$span / 60, d$drift_run_lo, d$drift_run_hi,
                           d$drift_typ, d$dt_typ, d$jitter, 100 * (d$frac_directed %||% NA_real_),
                           format(d$n_pairs, big.mark = ","),
                           if (is.finite(d$n_units)) sprintf(" in %s multi-blink molecules",
                                                             format(d$n_units, big.mark = ",")) else "")
            # Highlighted only when the drift is big enough to MATTER: a few tens
            # of nm across a whole run is excellent stability, and colouring it
            # as a warning would train the user to ignore the line. The bar is
            # the drift across a typical gap between blinks exceeding the
            # localization spread itself.
            htmltools::div(style = if (isTRUE(d$significant) &&
                                       isTRUE(d$drift_typ > (r$sigma_xy_med %||% Inf)))
                             "color:#b26a00;" else "color:#555;", txt)
          }),
          htmltools::div(style = "color:#777; font-size:0.85em;",
            "Only molecules whose blinks were merged can be seen here, so drift large enough to break a merge is understated. This is a sample readout and is not used to set the de-blink threshold.")
        )
      } else NULL)
  })


  # ==========================================================================
  # SAMPLE DRIFT PLOT - displacement from each molecule's FIRST blink
  # ==========================================================================
  # GOAL: show, rather than assert, whether the sample moved. Every molecule
  # that blinked more than once is placed at the origin at its first blink; each
  # later blink is one point, positioned relative to that origin and coloured by
  # the gap since it.
  #
  # HOW TO READ IT: a cloud centred on the origin is a still sample - the
  # scatter is jitter, which does not accumulate. Points marching off in one
  # direction, further the longer the gap, are drift. The two are hard to tell
  # apart from a single number, which is the whole reason for the plot.
  #
  # The binned means carry the reading. Averaging kills the isotropic part, so
  # under pure jitter every bin sits on the origin whatever its time gap, while
  # under drift the bins step outward in order. That is a far stronger cue than
  # the raw scatter, which drift barely changes when jitter dominates.
  output$prec_drift_plot <- renderPlotly({
    r <- prec_cloud_results()
    validate(need(!is.null(r), "Run the precision cloud first."))
    dl <- r$drift %||% list()
    validate(need(length(dl) > 0, "Sample drift could not be measured (needs traces with times, and molecules that blinked more than once)."))
    st <- do.call(rbind, lapply(names(dl), function(fn) {
      s <- dl[[fn]]$steps; if (is.null(s) || !nrow(s)) return(NULL)
      s$file <- fn
      # Molecule ids are per-file integers, so two files would otherwise share
      # ids and be drawn as one molecule.
      if (!is.null(s$mol)) s$mol <- paste(fn, s$mol, sep = " / ")
      s
    }))
    validate(need(!is.null(st) && nrow(st) >= 5, "Too few multi-blink molecules to plot."))

    ax  <- intersect(c("dx", "dy", "dz"), names(st))
    lab <- c(dx = "Δx", dy = "Δy", dz = "Δz")
    pairs_ax <- if (length(ax) >= 3) list(c("dx","dy"), c("dx","dz"), c("dy","dz")) else list(ax[1:2])

    # One frame for all panels so the eye compares them directly; driven by the
    # data rather than the fit, and by a quantile so a single wild pair cannot
    # zoom the plot out to nothing.
    L <- max(stats::quantile(abs(unlist(st[ax])), 0.99, na.rm = TRUE, names = FALSE) * 1.15, 5)

    # Only draw a fitted drift arrow when there is exactly ONE file: drift has a
    # direction, and two acquisitions do not share one, so a single arrow over
    # pooled files would be meaningless.
    one <- if (length(dl) == 1L) dl[[1]] else NULL
    cmode  <- input$drift_plot_color %||% "time"
    # Colouring by molecule needs a molecule to name. Spatial pairing (no
    # de-blinked IDs available) leaves `mol` NA, so fall back to time.
    by_mol <- identical(cmode, "mol") && "mol" %in% names(st) && any(!is.na(st$mol))
    # Acquisition time needs t_abs; older cached results will not carry it.
    by_abs <- identical(cmode, "abs") && "t_abs" %in% names(st) && any(is.finite(st$t_abs))
    cvals  <- (if (by_abs) st$t_abs else st$dt) / 60
    cpal   <- if (by_abs) "Turbo" else "Viridis"
    ctitle <- if (by_abs) "run time (min)" else "gap (min)"

    # ── Geometric median displacement per time-gap bin ───────────────────────
    # Estimated ONCE in full dimensionality and then projected into each panel,
    # so all three panels show the same estimated point viewed three ways. A
    # per-panel 2D fit would give the XY marker a different underlying point
    # from the XZ marker, which cannot be read as one displacement.
    #
    # GEOMETRIC median, not the mean and not a per-axis median. The mean has no
    # breakdown point and these bins really are contaminated - a merged
    # "molecule" occasionally holds two different fluorophores, dropping a
    # tens-of-nm displacement into an otherwise ~10 nm cloud. A per-axis median
    # is robust but NOT rotation-equivariant, so its marker would move if x and
    # y were relabelled, which is fatal for a plot whose content is a direction.
    # The geometric median is both robust and equivariant.
    binc <- NULL
    nb <- min(5L, max(2L, floor(nrow(st) / 8)))
    br <- unique(stats::quantile(st$dt, seq(0, 1, length.out = nb + 1), na.rm = TRUE))
    if (length(br) > 2) {
      bin  <- cut(st$dt, br, include.lowest = TRUE)
      idx  <- split(seq_len(nrow(st)), bin)
      idx  <- idx[lengths(idx) > 0]
      if (length(idx) >= 2) {
        cen <- t(vapply(idx, function(ix)
          geometric_median(as.matrix(st[ix, ax, drop = FALSE])), numeric(length(ax))))
        colnames(cen) <- ax
        binc <- list(cen = cen,
                     t = vapply(idx, function(ix) stats::median(st$dt[ix], na.rm = TRUE), numeric(1)))
      }
    }

    np  <- length(pairs_ax)
    gap <- if (np > 1) 0.07 else 0
    w   <- (1 - gap * (np - 1)) / np
    fig <- plotly::plot_ly(); shapes <- list(); axis_args <- list()

    for (i in seq_len(np)) {
      a <- pairs_ax[[i]][1]; b <- pairs_ax[[i]][2]
      xid <- if (i == 1) "x" else paste0("x", i)
      yid <- if (i == 1) "y" else paste0("y", i)

      if (by_mol) {
        # One colour per molecule, cycled: with hundreds of molecules a colour
        # cannot IDENTIFY one, but it does group a molecule's blinks so the eye
        # can follow them. Each molecule is joined from the origin in time order,
        # which turns every multi-blink fluorophore into its own small track -
        # under drift those tracks all lean the same way.
        pal <- grDevices::hcl.colors(12, "Dark 3")
        mm  <- split(seq_len(nrow(st)), st$mol)
        # Cap the number drawn. Each molecule is its own plotly trace, times the
        # number of panels, so an unbounded count stalls the browser rather than
        # informing anyone. The molecules with the MOST blinks are kept, since a
        # three- or four-blink track shows a direction and a single step cannot.
        cap <- 200L
        n_mol_all <- length(mm)
        if (n_mol_all > cap) mm <- mm[order(-lengths(mm))][seq_len(cap)]
        for (mi in seq_along(mm)) {
          ix <- mm[[mi]]
          fig <- fig %>% plotly::add_trace(
            x = c(0, st[[a]][ix]), y = c(0, st[[b]][ix]), xaxis = xid, yaxis = yid,
            type = "scatter", mode = "lines+markers",
            line = list(color = pal[(mi - 1L) %% length(pal) + 1L], width = 1),
            marker = list(size = 6, opacity = 0.75,
                          color = pal[(mi - 1L) %% length(pal) + 1L]),
            name = if (n_mol_all > cap)
                     sprintf("molecule track (%s of %s shown)", cap,
                             format(n_mol_all, big.mark = ","))
                   else "molecule track",
            legendgroup = "mol",
            showlegend = (i == 1 && mi == 1), inherit = FALSE,
            text = paste("de-blinked ID", st$mol[ix][1]),
            hovertemplate = paste0("%{text}<br>", lab[[a]], " %{x:.1f} nm<br>",
                                   lab[[b]], " %{y:.1f} nm<extra></extra>"))
        }
      } else {
        # Two continuous colourings, both drawn from the app's own palette set:
        #   time - how long AFTER its molecule's first blink (Viridis)
        #   abs  - WHEN in the acquisition it happened (Turbo, blue -> red)
        # They answer different questions. The gap shows accumulation; the
        # acquisition time shows whether early and late parts of the run pull in
        # different directions, which a gap colouring cannot separate.
        fig <- fig %>% plotly::add_markers(
          x = st[[a]], y = st[[b]], xaxis = xid, yaxis = yid,
          marker = list(size = 6, opacity = 0.55, color = cvals,
                        colorscale = plotly_colorscale(cpal),
                        showscale = (i == 1),
                        colorbar = list(title = ctitle, len = 0.9, thickness = 12)),
          name = "blink displacement", legendgroup = "pts", showlegend = FALSE,
          inherit = FALSE, text = sprintf("%.1f min", cvals),
          hovertemplate = paste0(lab[[a]], " %{x:.1f} nm<br>", lab[[b]],
                                 " %{y:.1f} nm<br>", ctitle,
                                 " %{text}<extra></extra>"))
      }

      # The per-bin centres, projected into this panel (computed once above).
      if (!is.null(binc)) {
        mx <- binc$cen[, a]; my <- binc$cen[, b]; mt <- binc$t
        keep <- is.finite(mx) & is.finite(my)
        # add_trace, not add_paths: add_paths forces mode = "lines" and plotly
        # then warns and re-adds markers on every call. Neither sorts by x, so
        # the point order is preserved either way (add_LINES is the one to
        # avoid - it sorts, and would scramble this trajectory).
        if (sum(keep) >= 2)
          fig <- fig %>% plotly::add_trace(
            x = c(0, mx[keep]), y = c(0, my[keep]), xaxis = xid, yaxis = yid,
            type = "scatter", mode = "lines+markers",
            line = list(color = "#c03030", width = 2),
            marker = list(size = 10, color = "#c03030", symbol = "circle"),
            name = "geometric median per time bin",
            legendgroup = "bin", showlegend = (i == 1), inherit = FALSE,
            text = c("start", sprintf("%.0f min", mt[keep] / 60)),
            hovertemplate = "%{text}<br>%{x:.1f}, %{y:.1f} nm<extra></extra>")
      }

      # Fitted drift, extended over the observed range of gaps.
      if (!is.null(one) && is.finite(one$rate) && one$rate > 0) {
        v <- one$velocity; tmax <- max(st$dt, na.rm = TRUE)
        va <- v[[sub("^d", "", a)]]; vb <- v[[sub("^d", "", b)]]
        if (is.finite(va) && is.finite(vb))
          fig <- fig %>% plotly::add_paths(
            x = c(0, va * tmax), y = c(0, vb * tmax), xaxis = xid, yaxis = yid,
            line = list(color = "#2c7fb8", width = 2, dash = "dash"),
            name = "fitted drift", legendgroup = "fit", showlegend = (i == 1),
            hoverinfo = "skip", inherit = FALSE)
      }

      axl <- function(x0, x1, y0, y1) list(type = "line", xref = xid, yref = yid,
        x0 = x0, x1 = x1, y0 = y0, y1 = y1, line = list(color = "rgba(0,0,0,0.35)", width = 1))
      shapes <- c(shapes, list(axl(0, 0, -L, L), axl(-L, L, 0, 0)))

      dom <- c((i - 1) * (w + gap), (i - 1) * (w + gap) + w)
      # X yields domain to satisfy the 1:1 scale; only one axis of the pair may
      # declare scaleanchor or plotly widens a range to resolve the loop.
      axis_args[[if (i == 1) "xaxis" else paste0("xaxis", i)]] <- list(
        title = paste(lab[[a]], "from first blink (nm)"), range = c(-L, L),
        zeroline = FALSE, showgrid = FALSE, domain = dom, anchor = yid,
        scaleanchor = yid, scaleratio = 1, constrain = "domain")
      axis_args[[if (i == 1) "yaxis" else paste0("yaxis", i)]] <- list(
        title = paste(lab[[b]], "(nm)"), range = c(-L, L),
        zeroline = FALSE, showgrid = FALSE, domain = c(0, 1), anchor = xid)
    }

    do.call(plotly::layout, c(list(fig, shapes = shapes,
      showlegend = TRUE,
      legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.22,
                    yanchor = "top"),
      margin = list(t = 10, b = 90)), axis_args)) %>%
      plotly::config(displaylogo = FALSE)
  })

  output$download_prec_cloud <- downloadHandler(
    filename = function() paste0("precision_cloud_perfile_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      r <- prec_cloud_results(); shiny::validate(shiny::need(!is.null(r), "No results."))
      write.csv(r$per_file, file, row.names = FALSE)
    })

  # ==========================================================================
  # QC 3 - Fiducial / bright-contaminant flagging (report only)
  # ==========================================================================
  # GOAL: flag traces that behave like gold fiducials / bright aggregates rather
  # than single blinking molecules. Learned from a labelled pair of real datasets
  # (one with a fiducial section, one clean), fiducials separate on TWO robust
  # signatures, NOT on spread (spread does not distinguish them):
  #   • BRIGHT   - efo is an extreme upper outlier. A gold bead reflects ~100× a
  #                fluorophore (efo ~1e7 vs ~1e5), so log10(efo) sits many robust
  #                MADs above the bulk. This is the single strongest signal.
  #   • PERSISTENT + RECURRENT - the bead is re-localized far more than a blinking
  #                dye (long trace) AND many such traces pile up at the SAME spot
  #                (fiducials recur in one place; real labels are spread out).
  # Scores are per-file and ROBUST (median + k·MAD on log; MAD tolerates up to ~50%
  # contamination, so a heavy fiducial section does not hide itself). On the clean
  # control file this flags 0 traces; on the fiducial file it flags a coherent,
  # spatially-clustered set. REPORT ONLY unless "Exclude flagged" is ticked.
  fiducial_results <- reactiveVal(NULL)
  observeEvent(input$run_fiducial, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    if (!isTRUE(input$fid_all_files)) {
      af <- input$active_file; if (!is.null(af) && af %in% names(store)) store <- store[af]
    }
    k_efo     <- max(2, input$fid_efo_z    %||% 4)      # brightness outlier (robust z)
    k_dur     <- max(1.5, input$fid_dur_x  %||% 5)      # duration, as a multiple of the file's q99
    min_locs  <- max(1, input$fid_min_locs %||% 5)      # min localizations for the brightness rule
    k_len     <- max(2, input$fid_len_z    %||% 4)      # length outlier (robust z), reported only
    recur_r   <- max(1, input$fid_recur_r  %||% 100)    # recurrence radius (nm), reported only
    recur_min <- max(1, input$fid_recur_min %||% 4)     # co-located traces, reported only
    # Robust z of log10(v): (log10 v − median) / MAD. 0 when MAD collapses (too few
    # / identical values) so a tiny or degenerate file flags nothing.
    rzlog <- function(v) {
      x <- log10(pmax(as.numeric(v), 1e-9))
      m <- stats::median(x, na.rm = TRUE); s <- stats::mad(x, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) rep(0, length(x)) else (x - m) / s
    }
    # For each trace centroid, count OTHER centroids within radius R (spatial
    # recurrence). Grid-bucketed so it stays O(n) on thousands of traces.
    recurrence <- function(cx, cy, R) {
      n <- length(cx); if (n < 2) return(integer(n))
      ox <- floor(cx / R); oy <- floor(cy / R)
      idx <- split(seq_len(n), paste(ox, oy)); cnt <- integer(n)
      for (i in seq_len(n)) {
        nb <- integer(0)
        for (a in -1:1) for (b in -1:1) {
          k <- paste(ox[i] + a, oy[i] + b); j <- idx[[k]]
          if (!is.null(j)) nb <- c(nb, j)
        }
        d <- sqrt((cx[nb] - cx[i])^2 + (cy[nb] - cy[i])^2)
        cnt[i] <- sum(d <= R) - 1L
      }
      cnt
    }
    rows <- list()
    withProgress(message = "Scanning for fiducials / contaminants...", value = 0, {
      nf <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]; incProgress(1 / nf, detail = fname)
        # A single odd file must not abort the whole scan - skip it and continue.
        rows[[fname]] <- tryCatch({
          df <- store[[fname]]
          if (!is.data.frame(df) || !nrow(df) || !all(c("tid","x","y") %in% names(df)))
            stop("no coords")
          # Shared predicate: z counts only if it VARIES (planar exports carry z = 0).
          hasz  <- has_axial_data(df)
          has_e <- "efo" %in% names(df); has_c <- "cfr" %in% names(df); has_t <- "tim" %in% names(df)
          feat <- df %>%
            dplyr::group_by(tid) %>%
            dplyr::summarise(
              n_loc    = dplyr::n(),
              duration = if (has_t) max(tim, na.rm = TRUE) - min(tim, na.rm = TRUE) else NA_real_,
              efo_med  = if (has_e) stats::median(efo, na.rm = TRUE) else NA_real_,
              cfr_med  = if (has_c) stats::median(cfr, na.rm = TRUE) else NA_real_,
              cx       = stats::median(x, na.rm = TRUE),
              cy       = stats::median(y, na.rm = TRUE),
              spread   = sqrt(mean((x - stats::median(x, na.rm = TRUE))^2 +
                                   (y - stats::median(y, na.rm = TRUE))^2 +
                                   (if (hasz) (z - stats::median(z, na.rm = TRUE))^2 else 0), na.rm = TRUE)),
              .groups  = "drop") %>%
            dplyr::mutate(dataset = fname)
          if (!nrow(feat)) stop("no traces")
          # Robust outlier scores + spatial recurrence (all per file).
          feat$z_efo <- if (has_e) rzlog(feat$efo_med) else rep(NA_real_, nrow(feat))
          feat$z_len <- rzlog(feat$n_loc)
          feat$recur <- recurrence(feat$cx, feat$cy, recur_r)
          # ── HOW LONG THE TRACKER STAYED LOCKED ON, RELATIVE TO THIS FILE ──────
          # A fiducial does not bleach, so the sequence keeps re-finding it and its
          # trace outlives every real fluorophore in the same acquisition. That is
          # the primary signature; brightness is secondary, and on some samples the
          # bead is not even the brightest thing in the file.
          # The reference is the file's own 99th percentile of the durations of
          # traces that live at all (duration > 0) - a file can be mostly
          # single-localization traces, which would collapse a median/MAD scale.
          live <- feat$duration[is.finite(feat$duration) & feat$duration > 0]
          dur_q99 <- if (length(live) >= 20L) unname(stats::quantile(live, 0.99)) else NA_real_
          feat$dur_x_q99 <- if (is.finite(dur_q99) && dur_q99 > 0) feat$duration / dur_q99 else NA_real_
          long_lived <- is.finite(feat$dur_x_q99) & feat$dur_x_q99 >= k_dur
          # Brightness outlier, but only for a trace that persisted: a single very
          # bright localization is a blip, not a bead. Measured on the labelled set
          # (Test data/Fiducial test): every brightness-only false positive in the
          # bead-free files was a 1-3 localization trace at 850-5700 kHz.
          bright <- has_e & is.finite(feat$z_efo) & feat$z_efo >= k_efo & feat$n_loc >= min_locs
          feat$flagged <- long_lived | bright
          feat$rule <- ifelse(long_lived & bright, "long+bright",
                       ifelse(long_lived, "long-lived", ifelse(bright, "bright", "")))
          feat
        }, error = function(e) NULL)
      }
    })
    rows <- Filter(Negate(is.null), rows)
    res <- if (length(rows)) do.call(rbind, rows) else NULL
    if (is.null(res) || !nrow(res)) {
      # Not an error - just nothing to look at. Store an empty (0-row) result so
      # the outputs show a calm "nothing found" state rather than a red message.
      fiducial_results(res %||% data.frame())
      showNotification("Scan complete - no traces to assess in the selected files.",
                       type = "message", duration = 5); return()
    }
    fiducial_results(res)
    nfl <- sum(res$flagged)
    showNotification(
      if (nfl == 0)
        sprintf("Scan complete - no fiducial/contaminant candidates found in %d trace(s) across %d file(s).",
                nrow(res), length(unique(res$dataset)))
      else
        sprintf("Flagged %d of %d traces (%.1f%%) as fiducial/contaminant candidates across %d file(s).",
                nfl, nrow(res), 100 * mean(res$flagged), length(unique(res$dataset))),
      type = "message", duration = 5)
  })

  output$fiducial_plot <- renderPlotly({
    res <- fiducial_results()
    # Calm empty state (no red validate) before a run or when nothing was found.
    if (is.null(res) || !nrow(res) || !"flagged" %in% names(res))
      # Declare the trace type (scatter/markers) so plotly does not print a
      # "No trace type specified" console warning for this empty placeholder.
      return(plotly::plot_ly(type = "scatter", mode = "markers") %>% plotly::layout(
        xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
        annotations = list(list(text = "Press 'Run' to scan for fiducials / contaminants.",
                                showarrow = FALSE, font = list(color = "#888", size = 13)))))
    # The two discriminating axes: brightness (efo, log - gold reflects ~100×) and
    # persistence (localizations per trace, log). Flagged traces sit in the bright
    # and/or long extremes. Falls back to spread on the x-axis when efo is absent.
    has_e <- "efo_med" %in% names(res) && any(is.finite(res$efo_med))
    xv <- if (has_e) pmax(res$efo_med, 1) else res$spread
    xl <- if (has_e) "Trace median EFO (log)" else "Trace spread (nm)"
    ok <- !res$flagged
    ht <- function(m) paste0("TID ", res$tid[m], "<br>", res$dataset[m],
      "<br>", res$n_loc[m], " locs · spread ", round(res$spread[m], 1), " nm",
      if (has_e) paste0("<br>efo ", format(round(res$efo_med[m]), big.mark = ","),
                        " (z ", round(res$z_efo[m], 1), ")") else "",
      "<br>len z ", round(res$z_len[m], 1),
      if ("recur" %in% names(res)) paste0(" · ", res$recur[m], " neighbours") else "")
    # Add each layer ONLY when it has points - an add_markers() with empty x/y is
    # what makes plotly warn "No trace type / no positional attributes specified".
    p <- plotly::plot_ly()
    if (any(ok))
      p <- p %>% plotly::add_markers(x = xv[ok], y = res$n_loc[ok], name = "traces",
        marker = list(color = "rgba(120,120,120,0.5)", size = 5,
                      line = list(color = "white", width = 0.3)),
        text = ht(ok), hoverinfo = "text")
    if (any(!ok))
      p <- p %>% plotly::add_markers(x = xv[!ok], y = res$n_loc[!ok], name = "flagged",
        marker = list(color = "#d62728", size = 8, symbol = "x"),
        text = ht(!ok), hoverinfo = "text")
    # If nothing was flagged, say so on the plot (the points are still shown).
    anns <- if (!any(!ok)) list(list(text = "No fiducials / contaminants detected",
      xref = "paper", yref = "paper", x = 0.5, y = 1, showarrow = FALSE,
      font = list(color = "#2a6", size = 13))) else NULL
    p %>% plotly::layout(
      xaxis = list(title = xl, type = if (has_e) "log" else "linear"),
      yaxis = list(title = "Localizations per trace", type = "log"),
      legend = list(orientation = "h"), annotations = anns)
  })

  output$fiducial_table <- renderTable({
    res <- fiducial_results()
    if (is.null(res) || !nrow(res) || !"flagged" %in% names(res)) return(NULL)
    f <- res[isTRUE(res$flagged) | res$flagged %in% TRUE, , drop = FALSE]
    if (!nrow(f)) return(data.frame(Message = "No traces flagged with the current thresholds."))
    f <- f[order(-f$z_efo, -f$n_loc), , drop = FALSE]
    utils::head(data.frame(
      dataset = f$dataset, tid = f$tid, n_loc = f$n_loc,
      efo_med = round(f$efo_med, 0), z_efo = round(f$z_efo, 1),
      z_len = round(f$z_len, 1), neighbours = f$recur,
      spread_nm = round(f$spread, 1), cfr_med = round(f$cfr_med, 3),
      check.names = FALSE), 200L)
  }, striped = TRUE, spacing = "xs", width = "100%")

  output$download_fiducial <- downloadHandler(
    filename = function() paste0("fiducial_scan_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- fiducial_results(); shiny::validate(shiny::need(!is.null(res), "No results."))
      write.csv(res, file, row.names = FALSE)
    })

  # Optional exclusion: when "Exclude flagged traces" is on, drop the traces the
  # Fiducial scan flagged (for that file) from a data frame before it is analysed.
  # NON-DESTRUCTIVE and reversible - the loaded data is untouched; unticking the box
  # (or re-running the scan) restores everything. tid excludes are matched PER FILE
  # (tids are only unique within a file), so `fname` must be the file the df is from.
  # Applied both in the active-file pipeline and in every batch analysis.
  qc_exclude_from <- function(df, fname) {
    if (!isTRUE(input$fid_exclude) || is.null(df) || !nrow(df) || !"tid" %in% names(df))
      return(df)
    fr <- tryCatch(fiducial_results(), error = function(e) NULL)
    if (is.null(fr) || !all(c("flagged","tid","dataset") %in% names(fr))) return(df)
    ex <- fr$tid[fr$flagged %in% TRUE & fr$dataset == fname]
    if (length(ex)) df[!(df$tid %in% ex), , drop = FALSE] else df
  }

  # ── Shared within-cluster NN statistics ─────────────────────────────────────
  # SINGLE SOURCE OF TRUTH for "how close are the members of a median cluster",
  # used by BOTH the Crowding panel and the Within-cluster-NN-vs-size panel so the
  # two report the SAME number. Both panels build their member centroids on the
  # same unit (de-blinked molecule when de-blinking is on, else trace ID; see
  # smlm_nn_unit below) and pass the member coordinates here.
  #
  # For each member, mean distance to its k nearest members (k clamped to n-1;
  # k=1 => nearest-neighbour spacing), then the MEDIAN over members with a robust
  # spread (MAD and the inter-quartile range) as the ± error.
  # pts = numeric matrix of member coordinates (one row per member).
  member_nn_stats <- function(pts, k = 1L) {
    n <- nrow(pts)
    if (is.null(n) || n < 2L) return(NULL)
    k_eff <- min(max(1L, as.integer(k)), n - 1L)
    nn <- tryCatch(dbscan::kNN(pts, k = k_eff)$dist, error = function(e) NULL)
    if (is.null(nn)) return(NULL)
    spacing <- if (k_eff == 1L) as.numeric(nn) else rowMeans(nn, na.rm = TRUE)
    qs <- as.numeric(stats::quantile(spacing, c(0.25, 0.75), na.rm = TRUE))
    list(median = stats::median(spacing, na.rm = TRUE),
         mad    = stats::mad(spacing, na.rm = TRUE),
         q1 = qs[1], q3 = qs[2], k_used = k_eff, n = n)
  }
  # Which unit the within-cluster NN is measured between: de-blinked molecules
  # (primary `cluster`) when de-blinking is enabled and present, else trace IDs.
  smlm_nn_unit <- function(df) {
    if (deblink_on() && "cluster" %in% names(df) &&
        any(!is.na(df$cluster) & df$cluster > 0L)) "cluster" else "tid"
  }

  ## ── Median Cluster (cluster2) Crowding vs Detection-Time Shift ─────────────────────────────
  # One row per Cluster2 structure (built from its member centroids - de-blinked
  # molecules when de-blinking is on, else trace IDs):
  #   nn_median  = median k-NN distance among members (smaller = denser) ± MAD
  #   time_shift = median(member tim_med) - median(unclustered background tim_med),
  #                computed per file (positive = structure detected later)
  # Tests whether denser (more crowded) structures are detected later.
  crowd_results    <- reactiveVal(NULL)
  crowd_null_stats <- reactiveVal(NULL)   # list(rho_obs, p, B, outside) from the permutation null

  observeEvent(input$run_crowd, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    if (!cluster2_on()) {
      showNotification(sprintf("Enable %s in the Display tab first.", cluster2_label()),
                       type = "warning"); return()
    }
    if (!isTRUE(input$crowd_all_files)) {
      af <- input$active_file
      if (!is.null(af) && af %in% names(store)) store <- store[af]
    }

    min_n  <- max(3L, as.integer(input$crowd_min_members %||% 4L))
    kk     <- max(1L, as.integer(input$crowd_k %||% 1L))
    use_3d <- isTRUE(input$crowd_use_3d)
    cols   <- if (use_3d) c("x", "y", "z") else c("x", "y")

    # ── no-crowding null (label permutation) ──────────────────────────────────
    # The panel's y-axis has a trend built in even with NO crowding effect: denser
    # clusters hold more traces, and the MEDIAN of a larger sample drawn from a
    # decaying detection-rate distribution lands later purely by sampling. To make
    # the plot readable we need the "zero crowding" expectation drawn on it.
    #
    # WHY A PERMUTATION AND NOT THE SIMULATOR: shuffling the detection times among
    # the traces of a file keeps the real cluster membership, the real sizes, the
    # real geometry (so every point's x stays exactly where it is) and the real
    # detection-time distribution of that acquisition - it destroys ONLY the
    # association being tested. A simulator null would instead require assuming a
    # bleaching/detection model, which is itself calibrated from data of this kind;
    # that would make the test partly circular. The permutation is exact,
    # assumption-free and cheaper.
    do_null <- isTRUE(input$crowd_null)
    B       <- if (do_null) max(49L, min(2000L, as.integer(input$crowd_nsim %||% 200L))) else 0L

    rows      <- list()
    null_rows <- list()   # one B-length vector of null time_shifts per cluster row
    # Why-it-found-nothing counters. An empty result has several distinct causes and
    # the old generic message named only one of them, so a stale filter or a missing
    # background was indistinguishable from "no structures".
    diag <- c(no_tim = 0L, filtered_out = 0L, no_c2 = 0L, no_z = 0L,
              no_clusters = 0L, no_background = 0L, too_small = 0L)
    withProgress(message = "Computing crowding vs detection-time shift...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- qc_exclude_from(store[[fname]], fname)
        if (!is.data.frame(df) || !nrow(df) || !"tim" %in% names(df)) { diag["no_tim"] <- diag["no_tim"] + 1L; next }
        df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
        # NOTE: the column is `cluster2`. "Median cluster" is only ever the *display*
        # name (hover labels); testing for it here skipped every file unconditionally,
        # so this panel could never return a result.
        if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) { diag["no_c2"] <- diag["no_c2"] + 1L; next }
        if (use_3d && !"z" %in% names(df)) { diag["no_z"] <- diag["no_z"] + 1L; next }

        # Members = de-blinked molecules (when de-blinking is on) or trace IDs -
        # the SAME unit the Within-cluster-NN panel uses, so the two agree. When
        # using molecules, primary noise (cluster 0) is dropped; each molecule's
        # detection time is the median tim of its traces.
        id_col <- smlm_nn_unit(df)
        base   <- if (id_col == "cluster") dplyr::filter(df, !is.na(cluster) & cluster > 0L) else df
        per_tid <- base %>%
          dplyr::group_by(.data[[id_col]]) %>%
          dplyr::summarise(
            tim_med = stats::median(tim, na.rm = TRUE),
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
            cl2 = { v <- cluster2[!is.na(cluster2) & cluster2 > 0L]; if (length(v)) v[1] else 0L },
            .groups = "drop"
          ) %>%
          dplyr::filter(is.finite(tim_med))
        if (!nrow(per_tid)) next

        cids <- unique(per_tid$cl2[per_tid$cl2 > 0L])
        if (!length(cids)) { diag["no_clusters"] <- diag["no_clusters"] + 1L; next }

        bg_med <- stats::median(per_tid$tim_med[per_tid$cl2 == 0L], na.rm = TRUE)
        if (!is.finite(bg_med)) {        # need an unclustered background to reference
          diag["no_background"] <- diag["no_background"] + 1L; next
        }

        keep_idx <- list()   # row indices of each qualifying cluster, for the permutation
        for (cid in cids) {
          which_c <- which(per_tid$cl2 == cid)
          sub <- per_tid[which_c, , drop = FALSE]
          n   <- nrow(sub)
          if (n < min_n) { diag["too_small"] <- diag["too_small"] + 1L; next }
          pts <- as.matrix(sub[, cols, drop = FALSE])
          pts <- pts[stats::complete.cases(pts), , drop = FALSE]
          if (nrow(pts) < 2L) next
          st <- member_nn_stats(pts, k = kk)
          if (is.null(st)) next

          rows[[length(rows) + 1L]] <- data.frame(
            dataset    = fname,
            cluster2   = cid,
            n_tids     = n,
            nn_median  = st$median,
            nn_mad     = st$mad,
            unit       = if (id_col == "cluster") "molecule" else "tid",
            time_shift = stats::median(sub$tim_med) - bg_med,
            stringsAsFactors = FALSE
          )
          keep_idx[[length(keep_idx) + 1L]] <- which_c
        }

        # Permute detection times WITHIN this file; geometry and membership untouched.
        if (do_null && length(keep_idx)) {
          tv     <- per_tid$tim_med
          is_bg  <- per_tid$cl2 == 0L
          nullm  <- matrix(NA_real_, nrow = length(keep_idx), ncol = B)
          for (b in seq_len(B)) {
            tv_b   <- tv[sample.int(length(tv))]          # shuffle times across all traces
            bg_b   <- stats::median(tv_b[is_bg], na.rm = TRUE)
            nullm[, b] <- vapply(keep_idx,
                                 function(ix) stats::median(tv_b[ix], na.rm = TRUE) - bg_b,
                                 numeric(1))
          }
          for (r in seq_len(nrow(nullm))) null_rows[[length(null_rows) + 1L]] <- nullm[r, ]
        }
      }
    })

    res <- if (length(rows)) do.call(rbind, rows) else NULL
    if (is.null(res) || !nrow(res)) {
      why <- c(
        if (diag["no_tim"])        sprintf("%d file(s) have no `tim` column", diag["no_tim"]),
        if (diag["filtered_out"])  sprintf("%d file(s) left empty by the Display-tab metadata filters", diag["filtered_out"]),
        if (diag["no_c2"])         sprintf("%d file(s) produced no cluster2 column (is Median Clustering enabled?)", diag["no_c2"]),
        if (diag["no_z"])          sprintf("%d file(s) have no z (uncheck 'Use 3-D coordinates')", diag["no_z"]),
        if (diag["no_clusters"])   sprintf("%d file(s) yielded no Median Cluster structures at all", diag["no_clusters"]),
        if (diag["no_background"]) sprintf("%d file(s) had no UNCLUSTERED traces to use as the time reference", diag["no_background"]),
        if (diag["too_small"])     sprintf("%d structure(s) had fewer than %d traces", diag["too_small"], min_n))
      showNotification(
        HTML(paste0("<b>No structures to plot.</b><br>",
                    if (length(why)) paste(why, collapse = "<br>") else
                      "No qualifying Median Cluster structures found.")),
        type = "warning", duration = 12)
      crowd_results(NULL); crowd_null_stats(NULL); return()
    }
    # ── attach the null band + the global trend test ──────────────────────────
    crowd_null_stats(NULL)
    if (do_null && length(null_rows) == nrow(res)) {
      NM <- do.call(rbind, null_rows)                     # clusters x B
      qs <- t(apply(NM, 1, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
      res$null_lo  <- qs[, 1]; res$null_med <- qs[, 2]; res$null_hi <- qs[, 3]

      # Global test on the TREND, not the offset: Spearman rho between crowding and
      # time shift, observed vs its permutation distribution. Rank-based, so it is
      # unaffected by a constant clustered-vs-background offset and by the skew of
      # the shift distribution - it asks only "does the shift depend on crowding?".
      rho_of <- function(y) suppressWarnings(
        stats::cor(res$nn_median, y, method = "spearman", use = "complete.obs"))
      rho_obs  <- rho_of(res$time_shift)
      rho_null <- apply(NM, 2, rho_of)
      rho_null <- rho_null[is.finite(rho_null)]
      p_two <- if (is.finite(rho_obs) && length(rho_null))
        (1 + sum(abs(rho_null) >= abs(rho_obs))) / (length(rho_null) + 1) else NA_real_
      outside <- mean(res$time_shift < res$null_lo | res$time_shift > res$null_hi, na.rm = TRUE)
      crowd_null_stats(list(rho_obs = rho_obs, p = p_two, B = length(rho_null),
                            outside = outside))
    } else {
      res$null_lo <- NA_real_; res$null_med <- NA_real_; res$null_hi <- NA_real_
    }

    crowd_results(res)
    ns <- crowd_null_stats()
    showNotification(
      sprintf("Crowding vs detection-time shift: %d structures across %d file(s).%s",
              nrow(res), length(unique(res$dataset)),
              if (!is.null(ns) && is.finite(ns$p))
                sprintf(" Trend vs no-crowding null: rho = %.3f, p = %.3f (%d permutations).",
                        ns$rho_obs, ns$p, ns$B) else ""),
      type = "message", duration = 6)
  })

  build_crowd_plot_fig <- function() {
    res <- crowd_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Press 'Run' to compute crowding vs detection-time shift."))
    n_all <- nrow(res)
    res <- apply_dataset_groups(res)
    shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
    # Group mode can drop files; the stored rho/p were computed over the full run,
    # so only quote them when every structure is still on screen.
    stats_match <- nrow(res) == n_all
    fig <- plotly::plot_ly()
    datasets <- unique(res$dataset)
    n_ds     <- length(datasets)
    pal <- series_palette(datasets)

    # ── no-crowding null band (drawn first, so points sit on top) ─────────────
    # Each cluster has its own null interval (its size is fixed), so the band is
    # per-cluster quantiles ordered along x and lightly smoothed for readability.
    if ("null_lo" %in% names(res) && any(is.finite(res$null_lo))) {
      nb  <- res[is.finite(res$null_lo) & is.finite(res$nn_median), , drop = FALSE]
      nb  <- nb[order(nb$nn_median), , drop = FALSE]
      smooth_band <- function(y) {
        if (nrow(nb) < 8) return(y)
        tr <- tryCatch(stats::loess(y ~ nb$nn_median, span = 0.9, degree = 1), error = function(e) NULL)
        if (is.null(tr)) return(y)
        p <- stats::predict(tr); ifelse(is.finite(p), p, y)
      }
      lo <- smooth_band(nb$null_lo); hi <- smooth_band(nb$null_hi); md <- smooth_band(nb$null_med)
      fig <- fig %>%
        plotly::add_lines(x = nb$nn_median, y = hi, name = "No-crowding null (95%)",
                          line = list(color = "rgba(120,120,120,0)", width = 0),
                          hoverinfo = "skip", showlegend = FALSE) %>%
        plotly::add_lines(x = nb$nn_median, y = lo, name = "No-crowding null (95%)",
                          fill = "tonexty", fillcolor = "rgba(120,120,120,0.18)",
                          line = list(color = "rgba(120,120,120,0)", width = 0),
                          hoverinfo = "skip") %>%
        plotly::add_lines(x = nb$nn_median, y = md, name = "Null median",
                          line = list(color = "#777777", width = 1.5, dash = "dot"),
                          hoverinfo = "skip")
    }

    for (ds in datasets) {
      sub <- res[res$dataset == ds, , drop = FALSE]
      ex  <- if ("nn_mad" %in% names(sub)) list(type = "data", array = sub$nn_mad,
                                                color = pal[ds], thickness = 1, width = 0,
                                                opacity = 0.4) else NULL
      fig <- fig %>% plotly::add_markers(
        x = sub$nn_median, y = sub$time_shift,
        name = if (n_ds > 1) ds else "Structures",
        error_x = ex,   # ± MAD of the within-cluster NN spacing
        marker = list(color = pal[ds], size = 7, opacity = 0.7,
                      line = list(color = "white", width = 0.5)),
        text = paste0("Median cluster ", sub$cluster2, "<br>n = ", sub$n_tids,
                      " ", if ("unit" %in% names(sub))
                             ifelse(sub$unit == "molecule", "de-blinked molecules", "trace IDs") else "",
                      "<br>NN median = ", round(sub$nn_median, 1), " ± ",
                      round(if ("nn_mad" %in% names(sub)) sub$nn_mad else 0, 1), " nm",
                      "<br>time shift = ", round(sub$time_shift, 2)),
        hoverinfo = "text"
      )
    }

    # optional loess trend - per group when grouped, else over the pooled cloud
    add_crowd_loess <- function(fig, d, col, nm) {
      if (!(isTRUE(input$crowd_trend) && nrow(d) >= 8)) return(fig)
      ord <- order(d$nn_median)
      tr  <- tryCatch(stats::loess(time_shift ~ nn_median, data = d,
                       span = 0.9, degree = 1, family = "symmetric"), error = function(e) NULL)
      if (is.null(tr)) return(fig)
      xx <- d$nn_median[ord]; yy <- stats::predict(tr, newdata = data.frame(nn_median = xx)); ok <- is.finite(yy)
      if (sum(ok) >= 2) fig <- fig %>% plotly::add_lines(x = xx[ok], y = yy[ok], name = nm,
                                 line = list(color = col, width = 2))
      fig
    }
    if ((input$analysis_by %||% "file") == "group" && length(datasets) > 1) {
      for (ds in datasets) fig <- add_crowd_loess(fig, res[res$dataset == ds, , drop = FALSE], unname(pal[ds]), paste(ds, "loess"))
    } else {
      fig <- add_crowd_loess(fig, res, "#D7301F", "Trend (loess)")
    }

    # Verdict line: the permutation test on the TREND (rho), plus how many points
    # left their own null interval. Stated explicitly so the slope is never read
    # against zero - the null band, not zero, is the correct reference here.
    ns  <- crowd_null_stats()
    ttl <- if (!is.null(ns) && is.finite(ns$p) && stats_match) {
      list(text = sprintf(
        "Crowding vs detection-time shift  |  trend vs no-crowding null: Spearman rho = %.3f, p = %.3f (%d permutations)  |  %.0f%% of structures outside their null band (5%% expected)",
        ns$rho_obs, ns$p, ns$B, 100 * ns$outside), font = list(size = 12))
    } else if (!is.null(ns) && is.finite(ns$p)) {
      list(text = "Crowding vs detection-time shift  |  null band shown per structure; trend test not shown (group filter changed the set)",
           font = list(size = 12))
    } else NULL

    fig %>% plotly::layout(
      title  = ttl,
      margin = list(l = 70, t = if (is.null(ttl)) 40 else 60),
      xaxis  = list(title = "Median NN distance among members (nm)  -  smaller = more crowded"),
      yaxis  = list(title = "Detection-time shift (tim)",
                    zeroline = TRUE, zerolinecolor = "#888888", zerolinewidth = 1),
      legend = list(orientation = "h")
    )
  }
  output$crowd_plot <- renderPlotly(build_crowd_plot_fig())

  output$download_crowd <- downloadHandler(
    filename = function()
      paste0("cluster2_crowding_vs_detection_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- crowd_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  ## ── Within-Cluster2 NN vs Cluster Size ────────────────────────────────────
  # One row per Cluster2: size = # member units (trace IDs, or de-blinked molecules
  # if de-blinking is on); nn_median = median 1-NN distance among the
  # member centroids. Tests whether spacing is constant with size or "balloons".
  size_nn_results <- reactiveVal(NULL)
  size_nn_params  <- reactiveVal(NULL)   # compute params, frozen at Run

  observeEvent(input$run_size_nn, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    if (!cluster2_on()) {
      showNotification(sprintf("Enable %s in the Display tab first.", cluster2_label()),
                       type = "warning"); return()
    }
    if (!isTRUE(input$size_nn_all_files)) {
      af <- input$active_file
      if (!is.null(af) && af %in% names(store)) store <- store[af]
    }
    min_n       <- max(2L, as.integer(input$size_nn_min_members %||% 3L))
    kk          <- max(1L, as.integer(input$size_nn_k %||% 1L))
    use_3d      <- isTRUE(input$size_nn_use_3d)
    cols        <- if (use_3d) c("x", "y", "z") else c("x", "y")
    use_primary <- deblink_on()

    rows <- list()
    # Same why-it-found-nothing counters as the crowding panel: an empty result has
    # several distinct causes and the old message named only the size threshold.
    diag <- c(empty = 0L, filtered_out = 0L, no_c2 = 0L, no_z = 0L,
              no_clusters = 0L, too_small = 0L, degenerate = 0L)
    withProgress(message = "Computing within-cluster NN vs size...", value = 0, {
      nf <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]; incProgress(1 / nf, detail = fname)
        df <- qc_exclude_from(store[[fname]], fname)
        if (!is.data.frame(df) || !nrow(df)) { diag["empty"] <- diag["empty"] + 1L; next }
        df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
        # NOTE: the column is `cluster2`. "Median Cluster" is only the DISPLAY name;
        # testing for it here skipped every file, so this panel could never return a
        # result. (Same defect as the crowding panel, which used "Median cluster".)
        if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) { diag["no_c2"] <- diag["no_c2"] + 1L; next }
        if (use_3d && !"z" %in% names(df)) { diag["no_z"] <- diag["no_z"] + 1L; next }

        id_col <- smlm_nn_unit(df)   # shared unit rule: molecule if de-blinking on, else tid
        cent <- df %>%
          dplyr::filter(!is.na(cluster2) & cluster2 > 0L) %>%
          { if (id_col == "cluster") dplyr::filter(., !is.na(cluster) & cluster > 0L) else . } %>%
          dplyr::group_by(cluster2, .data[[id_col]]) %>%
          dplyr::summarise(
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
            .groups = "drop")
        if (!nrow(cent)) { diag["no_clusters"] <- diag["no_clusters"] + 1L; next }

        for (cid in unique(cent$cluster2)) {
          sub <- cent[cent$cluster2 == cid, , drop = FALSE]
          n   <- nrow(sub)
          if (n < min_n) { diag["too_small"] <- diag["too_small"] + 1L; next }
          pts <- as.matrix(sub[, cols, drop = FALSE])
          pts <- pts[stats::complete.cases(pts), , drop = FALSE]
          if (nrow(pts) < 2L) { diag["degenerate"] <- diag["degenerate"] + 1L; next }
          # Shared k-NN spacing (median over members ± MAD) - identical definition
          # to the Crowding panel, so the two report the same within-cluster NN.
          st <- member_nn_stats(pts, k = kk)
          if (is.null(st)) next
          # convex-hull size: area (2D) / volume (3D); NA if degenerate
          hull <- tryCatch(geometry::convhulln(pts, options = "FA")$vol, error = function(e) NA_real_)
          rows[[length(rows) + 1L]] <- data.frame(
            dataset   = fname, cluster2 = cid, size = n,
            nn_median = st$median,
            nn_mad    = st$mad,
            hull_size = if (is.finite(hull) && hull > 0) hull else NA_real_,
            dim       = length(cols),
            unit      = if (id_col == "cluster") "molecule" else "tid",
            stringsAsFactors = FALSE)
        }
      }
    })

    res <- if (length(rows)) do.call(rbind, rows) else NULL
    if (is.null(res) || !nrow(res)) {
      why <- c(
        if (diag["empty"])        sprintf("%d file(s) were empty", diag["empty"]),
        if (diag["filtered_out"]) sprintf("%d file(s) left empty by the Display-tab metadata filters", diag["filtered_out"]),
        if (diag["no_c2"])        sprintf("%d file(s) produced no cluster2 column (is Median Clustering enabled?)", diag["no_c2"]),
        if (diag["no_z"])         sprintf("%d file(s) have no z (uncheck 'Use 3-D coordinates')", diag["no_z"]),
        if (diag["no_clusters"])  sprintf("%d file(s) yielded no Median Cluster structures at all", diag["no_clusters"]),
        if (diag["too_small"])    sprintf("%d structure(s) had fewer than %d members", diag["too_small"], min_n),
        if (diag["degenerate"])   sprintf("%d structure(s) had fewer than 2 usable coordinates", diag["degenerate"]))
      showNotification(
        HTML(paste0("<b>No structures to plot.</b><br>",
                    if (length(why)) paste(why, collapse = "<br>") else
                      "No qualifying Median Cluster structures found.")),
        type = "warning", duration = 12)
      size_nn_results(NULL); return()
    }
    size_nn_results(res)
    # Compute parameters frozen at Run: reading them live relabelled an existing
    # plot without recomputing it, so the axis could describe a different unit
    # than the data on screen.
    size_nn_params(list(k = kk, xmode = input$size_nn_xmode %||% "count"))
    rho <- suppressWarnings(stats::cor(res$size, res$nn_median, method = "spearman"))
    showNotification(
      sprintf("Within-cluster NN vs size: %d structures across %d file(s); Spearman rho = %.2f.",
              nrow(res), length(unique(res$dataset)), rho),
      type = "message", duration = 4)
  })

  build_size_nn_plot_fig <- function() {
    res <- size_nn_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Press 'Run' to compute within-cluster NN vs size."))
    res <- apply_dataset_groups(res)
    shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
    # res$unit holds "molecule" / "tid" (written by smlm_nn_unit above). The old
    # test compared against "primary_cluster", a name from when de-blinking was
    # raw-point DBSCAN, so it never matched and the label always read "trace IDs"
    # even when the panel was plotting de-blinked molecules.
    unit_lbl <- if (all(res$unit == "molecule")) "de-blinked molecules" else "trace IDs"
    xmode    <- (size_nn_params()$xmode) %||% (input$size_nn_xmode %||% "count")

    if (identical(xmode, "hull")) {
      res <- res[is.finite(res$hull_size), , drop = FALSE]
      shiny::validate(shiny::need(nrow(res) > 0,
        "No clusters with a valid convex hull (need >= 3 members in 2D / >= 4 in 3D)."))
      res$xv  <- res$hull_size
      is3d    <- all(res$dim == 3)
      xlab    <- if (is3d) "Convex-hull volume (nm^3)" else "Convex-hull area (nm^2)"
      xname   <- if (is3d) "hull vol" else "hull area"
      slope_u <- if (is3d) "nm per nm^3" else "nm per nm^2"
    } else {
      res$xv  <- res$size
      xlab    <- sprintf("Cluster size (# %s)", unit_lbl)
      xname   <- "size"
      slope_u <- "nm/member"
    }

    rho <- suppressWarnings(stats::cor(res$xv, res$nn_median, method = "spearman"))
    sl  <- tryCatch(unname(stats::coef(stats::lm(nn_median ~ xv, data = res))[2]),
                    error = function(e) NA_real_)

    fig <- plotly::plot_ly()
    datasets <- unique(res$dataset); n_ds <- length(datasets)
    pal <- series_palette(datasets)

    for (ds in datasets) {
      sub <- res[res$dataset == ds, , drop = FALSE]
      ey  <- if ("nn_mad" %in% names(sub)) list(type = "data", array = sub$nn_mad,
                                                color = pal[ds], thickness = 1, width = 0,
                                                opacity = 0.4) else NULL
      fig <- fig %>% plotly::add_markers(x = sub$xv, y = sub$nn_median,
        name = if (n_ds > 1) ds else "Structures",
        error_y = ey,   # ± MAD of the within-cluster NN spacing
        marker = list(color = pal[ds], size = 7, opacity = 0.7, line = list(color = "white", width = 0.5)),
        text = paste0("Median cluster ", sub$cluster2, "<br>", xname, " ", signif(sub$xv, 4),
                      "<br>median NN ", round(sub$nn_median, 1), " ± ",
                      round(if ("nn_mad" %in% names(sub)) sub$nn_mad else 0, 1), " nm<br>", sub$dataset),
        hoverinfo = "text")
    }

    grouped <- (input$analysis_by %||% "file") == "group" && length(datasets) > 1

    # dashed linear fit (always) + optional loess, over a data subset
    add_trends <- function(fig, d, col, prefix) {
      if (nrow(d) < 3) return(fig)
      sl_d  <- tryCatch(unname(stats::coef(stats::lm(nn_median ~ xv, data = d))[2]),
                        error = function(e) NA_real_)
      rho_d <- suppressWarnings(stats::cor(d$xv, d$nn_median, method = "spearman"))
      xr    <- range(d$xv, na.rm = TRUE)
      yr    <- tryCatch(stats::predict(stats::lm(nn_median ~ xv, data = d),
                        newdata = data.frame(xv = xr)), error = function(e) NULL)
      if (!is.null(yr)) {
        nm <- if (nzchar(prefix))
                sprintf("%s fit (slope %s, rho %.2f)", prefix,
                        if (is.finite(sl_d)) signif(sl_d, 3) else "NA", rho_d)
              else "Linear fit"
        fig <- fig %>% plotly::add_lines(x = xr, y = yr, name = nm,
                 line = list(color = col, width = 2, dash = "dash"))
      }
      if (isTRUE(input$size_nn_trend) && nrow(d) >= 8) {
        ord <- order(d$xv)
        tr  <- tryCatch(stats::loess(nn_median ~ xv, data = d, span = 0.9, degree = 1, family = "symmetric"),
                        error = function(e) NULL)
        if (!is.null(tr)) {
          xx <- d$xv[ord]; yy <- stats::predict(tr, newdata = data.frame(xv = xx)); ok <- is.finite(yy)
          if (sum(ok) >= 2)
            fig <- fig %>% plotly::add_lines(x = xx[ok], y = yy[ok],
                     name = if (nzchar(prefix)) paste(prefix, "loess") else "Trend (loess)",
                     line = list(color = col, width = 2))
        }
      }
      fig
    }

    if (grouped) {
      for (ds in datasets) fig <- add_trends(fig, res[res$dataset == ds, , drop = FALSE], unname(pal[ds]), ds)
      ttl <- "Within-cluster NN vs size - per-group fits (slope & rho in legend)"
    } else {
      fig <- add_trends(fig, res, "#555555", "")
      ttl <- sprintf("Spearman rho = %.2f   |   linear slope = %s %s",
                     rho, if (is.finite(sl)) signif(sl, 3) else "NA", slope_u)
    }
    kk   <- max(1L, as.integer((size_nn_params()$k) %||% (input$size_nn_k %||% 1L)))
    ylab <- if (kk == 1L) "Within-cluster median NN distance (nm)"
            else sprintf("Within-cluster median spacing (mean of %d nearest, nm)", kk)
    fig %>% plotly::layout(
      title  = list(text = ttl, font = list(size = 13)),
      xaxis  = list(title = xlab,
                    type = if (isTRUE(input$size_nn_logx)) "log" else "linear"),
      yaxis  = list(title = ylab),
      legend = list(orientation = "h"))
  }
  output$size_nn_plot <- renderPlotly(build_size_nn_plot_fig())

  output$size_nn_small_n <- renderUI({
    res <- size_nn_results()
    if (is.null(res) || !nrow(res)) return(NULL)
    # `size` is this panel's member count (see the run_size_nn result frame).
    nv <- if ("size" %in% names(res)) res$size else
          if ("n_tids" %in% names(res)) res$n_tids else NULL
    if (is.null(nv)) return(NULL)
    .small_n_note(nv, "structure")
  })

  output$download_size_nn <- downloadHandler(
    filename = function() paste0("within_cluster2_nn_vs_size_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- size_nn_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  ## ── Cluster Viewer (2D / 3D + shape fit) ──────────────────────────────────
  # Inspect one Median Cluster structure at a time. A Run button snapshots the
  # clusters; Prev/Next and view options then work on that snapshot, so it does
  # not refresh on filter/plot changes until Run is pressed again.
  # MLE best-fit via mle_fit_*; overlay from the fitted parameters.
  #
  # SPANS EVERY LOADED FILE, not just the plotted one. Median Clustering is a
  # global setting that applies to all loaded datasets, so the viewer was the odd
  # one out in reading only plottable_data(): with eight files loaded and the
  # active one a non-clustering control, it reported "no Median Clusters" and
  # refused to open - while seven other files were full of them.
  #
  # Cluster IDs are only unique WITHIN a file, so structures are keyed by
  # file + id (`.vkey`). Without that, cluster 3 of file A and cluster 3 of file B
  # would be the same entry and the viewer would show their points merged.
  viewer_snapshot <- reactiveVal(NULL)   # list(df, ids, keys) captured on Run
  viewer_idx      <- reactiveVal(1L)

  observeEvent(input$run_viewer, {
    all_files <- !isTRUE(input$viewer_active_only)
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning")
                          viewer_snapshot(NULL); return() }

    nms <- names(store)
    if (!all_files) {
      af <- input$active_file
      nms <- if (!is.null(af) && af %in% nms) af else nms[1]
    }

    parts <- list(); skipped <- character(0)
    withProgress(message = "Loading Median Clusters...", value = 0, {
      for (i in seq_along(nms)) {
        nm <- nms[[i]]
        incProgress(1 / length(nms), detail = nm)
        d <- tryCatch(prepared_df(nm, needs_cluster = TRUE), error = function(e) NULL)
        if (is.null(d) || !nrow(d) || !"cluster2" %in% names(d)) { skipped <- c(skipped, nm); next }
        d <- d[!is.na(d$cluster2) & d$cluster2 > 0L, , drop = FALSE]
        if (!nrow(d)) { skipped <- c(skipped, nm); next }
        d$.vfile <- nm
        d$.vkey  <- paste0(nm, "#", d$cluster2)
        parts[[nm]] <- d
      }
    })

    if (!length(parts)) {
      showNotification(sprintf("No Median Clusters in any loaded file. Enable %s and check the filters.",
                               cluster2_label()), type = "warning")
      viewer_snapshot(NULL); return()
    }
    df <- dplyr::bind_rows(parts)
    keys <- df %>%
      dplyr::distinct(.vkey, .vfile, cluster2) %>%
      dplyr::arrange(.vfile, cluster2)

    viewer_snapshot(list(df = df, ids = keys$.vkey, keys = keys))
    viewer_idx(1L)
    msg <- sprintf("Loaded %d Median Cluster structure(s) from %d file(s).",
                   nrow(keys), length(parts))
    if (length(skipped))
      msg <- paste0(msg, sprintf(" %d file(s) had none.", length(skipped)))
    showNotification(msg, type = "message", duration = 4)
  })

  viewer_ids <- reactive({ s <- viewer_snapshot(); if (is.null(s)) integer(0) else s$ids })

  observeEvent(input$viewer_prev, {
    n <- length(viewer_ids()); if (!n) return()
    viewer_idx(((viewer_idx() - 2L) %% n) + 1L)   # wrap-around
  })
  observeEvent(input$viewer_next, {
    n <- length(viewer_ids()); if (!n) return()
    viewer_idx((viewer_idx() %% n) + 1L)
  })

  viewer_current <- reactive({
    s   <- viewer_snapshot()
    shiny::validate(shiny::need(!is.null(s) && length(s$ids) > 0,
      "Press 'Run / Load clusters' to load Median Cluster structures from the loaded files."))
    ids <- s$ids; n <- length(ids)
    idx <- min(max(viewer_idx(), 1L), n)
    key <- ids[idx]
    df  <- s$df
    # Subset on the file+id key, never on cluster2 alone - see the note on Run.
    sub <- df[df$.vkey == key, , drop = FALSE]
    id  <- sub$cluster2[1]
    src <- sub$.vfile[1]
    # THE UNIT FOLLOWS DE-BLINKING, like every other panel.
    #
    # This used to group by `tid` unconditionally, so the viewer showed - and the
    # shape fit consumed - trace-ID medians even with de-blinking switched on.
    # That is the same defect the tc2 isosurface had: one fluorophore reappears as
    # several trace IDs, so a molecule that blinked six times contributed six
    # points to the surface and one that blinked once contributed one. The fitted
    # shape then tracks BLINKING as much as geometry, and the long-lived molecules
    # pull hardest. De-blinking exists precisely to collapse that, and the viewer
    # was ignoring it.
    unit_col <- if (deblink_on() && "cluster" %in% names(sub) &&
                    any(!is.na(sub$cluster) & sub$cluster > 0L)) "cluster" else "tid"
    csub <- if (identical(unit_col, "cluster"))
      sub[!is.na(sub$cluster) & sub$cluster > 0L, , drop = FALSE] else sub

    # Photon budget summed over the UNIT (all localizations of the trace, or of the
    # whole de-blinked molecule). MINFLUX localization variance scales as
    # 1/photons and the unit's centroid averages its localizations, so the
    # centroid's variance goes as 1/sum(eco) - the inverse-variance weight for
    # that centroid IS its total photon count.
    cent <- csub %>% dplyr::group_by(.data[[unit_col]]) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = if ("z" %in% names(csub)) stats::median(z, na.rm = TRUE) else 0,
        w_eco = if ("eco" %in% names(csub)) sum(eco, na.rm = TRUE) else NA_real_,
        .groups = "drop")
    w <- cent$w_eco
    if (!length(w) || all(!is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
      w <- NULL                                   # no usable eco -> unweighted
    } else {
      w[!is.finite(w) | w <= 0] <- stats::median(w[is.finite(w) & w > 0])
    }
    zc    <- if ("z" %in% names(sub)) sub$z else rep(0, nrow(sub))
    cloud <- cbind(x = sub$x, y = sub$y, z = zc)
    cloud <- cloud[stats::complete.cases(cloud), , drop = FALSE]
    list(idx = idx, n = n, id = id, file = src, sub = sub,
         unit = unit_col,
         mat = as.matrix(cent[, c("x", "y", "z")]),
         w_eco = w,
         cloud = cloud)
  })

  # Full MLE fit (optim orthogonal-distance) of all 5 shapes to the trace-median
  # centroids - same pathway the Display tab uses when fit_method = "mle", not
  # the heuristic moment-based geometry. Returns RMSDs, best (min-RMSD) shape and
  # the fitted parameters (used to draw the overlay).
  # Heavy part: the 5 MLE fits + RMSD + AIC/AICc/BIC. Independent of the ranking
  # criterion, so switching AICc<->BIC<->RMSD does NOT re-run the optimisers.
  viewer_fit_core <- reactive({
    if (!isTRUE(input$viewer_enable_fit)) return(NULL)   # inspect-only unless fitting enabled
    cur <- viewer_current()
    X   <- cur$mat
    if (nrow(X) < 4) return(NULL)
    turns <- input$helix_turns %||% 3
    cargs <- tryCatch(build_center_args(
      X, list(mle_center_mode   = input$mle_center_mode   %||% "free",
              mle_center_lambda = input$mle_center_lambda %||% 0.1)),
      error = function(e) list())

    # Inverse-variance (photon) weighting. The default iid-isotropic assumption is
    # wrong in detail for MINFLUX: precision varies with photon count, so a
    # 50-photon trace should not pull the surface as hard as a 5000-photon one.
    # w_i = total eco of trace i (see viewer_current). n_eff is Kish's effective
    # sample size - it exposes the cost of weighting, i.e. whether a couple of very
    # bright traces are now carrying the fit.
    W <- if (isTRUE(input$viewer_weight_eco)) cur$w_eco else NULL
    n_eff <- if (!is.null(W)) sum(W)^2 / sum(W^2) else nrow(X)

    fitS <- tryCatch(do.call(mle_fit_sphere,    c(list(X = X, w = W), cargs)), error = function(e) NULL)
    fitE <- tryCatch(do.call(mle_fit_ellipsoid, c(list(X = X, w = W), cargs)), error = function(e) NULL)
    fitH <- tryCatch(do.call(mle_fit_helix,     c(list(X = X, w = W, turns_hint = turns), cargs)), error = function(e) NULL)
    fitR <- tryCatch(do.call(mle_fit_ring,      c(list(X = X, w = W), cargs)), error = function(e) NULL)
    fitA <- tryCatch(do.call(mle_fit_arc,       c(list(X = X, w = W), cargs)), error = function(e) NULL)

    # Weighted RMSD: sqrt(sum(w*d^2)/sum(w)). With W = NULL this is the plain RMSD.
    rmsd_of <- function(d) {
      ok <- is.finite(d)
      if (!any(ok)) return(NA_real_)
      ww <- if (is.null(W)) rep(1, length(d)) else W
      sqrt(sum(ww[ok] * d[ok]^2) / sum(ww[ok]))
    }
    rS <- if (!is.null(fitS)) rmsd_of(as.numeric(mle_project_sphere(X, fitS$center, fitS$radius)$dist)) else NA_real_
    rE <- if (!is.null(fitE)) rmsd_of(as.numeric(mle_project_ellipsoid(X, fitE$center, fitE$axes, fitE$R)$dist)) else NA_real_
    rR <- if (!is.null(fitR)) rmsd_of(as.numeric(mle_project_ring(X, fitR$center, fitR$normal, fitR$radius)$dist)) else NA_real_
    rA <- if (!is.null(fitA)) rmsd_of(as.numeric(mle_project_arc(X, fitA$center, fitA$normal, fitA$radius, fitA$t_center, fitA$arc_half)$dist)) else NA_real_
    rH <- NA_real_
    if (!is.null(fitH) && !is.null(fitH$t)) {
      t_h <- .project_t_helix(X, fitH$center, fitH$axis, fitH$radius, fitH$pitch, fitH$phi0,
                              t0 = fitH$t, check_neighbors = TRUE)
      H_h <- t(vapply(t_h, .helix_point, numeric(3),
                      c = fitH$center, u = fitH$axis, r = fitH$radius, p = fitH$pitch, phi0 = fitH$phi0))
      rH <- rmsd_of(sqrt(rowSums((X - H_h)^2)))
    }
    rmsd <- c(sphere = rS, ellipsoid = rE, helix = rH, ring = rR, arc = rA)

    # ── model selection: AIC / AICc / BIC under isotropic Gaussian residuals ─────
    # sigma^2_hat = RMSD^2, so -2logL = n*log(2*pi*sigma2) + n. k = geometric DOF + 1
    # (the variance). AICc adds the small-sample correction and -> Inf when a model
    # has too few points to estimate (n <= k+1), which refuses over-parameterised fits.
    # When photon-weighted, point i has variance sigma^2/w_i and -2logL gains a
    # -sum(log w_i) term - identical for every shape (same weights), so it cancels in
    # every AIC/AICc/BIC *difference*. The criteria stay comparable across shapes;
    # only their absolute values shift, so never compare them across weighting modes.
    n   <- nrow(X)
    pk  <- c(sphere = 4, ellipsoid = 9, helix = 8, ring = 6, arc = 8)  # geometric DOF
    shp <- names(rmsd)
    M <- t(vapply(shp, function(s) {
      r <- rmsd[[s]]
      if (!is.finite(r)) return(c(rmsd = r, aic = NA_real_, aicc = NA_real_, bic = NA_real_, params = pk[[s]]))
      sigma2 <- max(r, 1e-6)^2
      k    <- pk[[s]] + 1L                       # + variance parameter
      n2ll <- n * log(2 * pi) + n * log(sigma2) + n
      aic  <- n2ll + 2 * k
      aicc <- if (n - k - 1L > 0) aic + (2 * k * (k + 1)) / (n - k - 1L) else Inf
      bic  <- n2ll + k * log(n)
      c(rmsd = r, aic = aic, aicc = aicc, bic = bic, params = pk[[s]])
    }, numeric(5)))
    rownames(M) <- shp

    list(rmsd = rmsd, M = M, n = n,
         weighted = !is.null(W), n_eff = n_eff,
         fits = list(sphere = fitS, ellipsoid = fitE, helix = fitH, ring = fitR, arc = fitA))
  })

  # Cheap part: pick the best shape for the chosen criterion (just a which.min).
  viewer_fit <- reactive({
    core <- viewer_fit_core()
    if (is.null(core)) return(NULL)
    M <- core$M
    crit_col <- switch(input$viewer_rank %||% "aicc",
                       aicc = "aicc", bic = "bic", rmsd = "rmsd", "aicc")
    pick_best <- function(col) {
      v <- M[, col]; v[!is.finite(v)] <- Inf
      if (all(!is.finite(v))) return(NA_character_)
      rownames(M)[which.min(v)]
    }
    best <- pick_best(crit_col)
    if (is.na(best) && crit_col == "aicc") best <- pick_best("aic")   # tiny clusters
    if (is.na(best)) best <- pick_best("rmsd")
    c(core, list(best = best, criterion = crit_col))
  })

  # Curve overlays (ring/arc/helix): ordered points to draw as a connected line.
  .viewer_curve_points <- function(f, shape) {
    pts <- tryCatch(switch(shape,
      ring = { p <- as.matrix(mle_sample_ring(f, n = 240)); rbind(p, p[1, ]) },  # close loop
      arc  = as.matrix(mle_sample_arc(f, n = 240)),
      helix = {
        if (!is.null(f$t)) {
          tr <- range(f$t, na.rm = TRUE)
          ts <- seq(tr[1], tr[2], length.out = 600L)
          t(vapply(ts, .helix_point, numeric(3),
                   c = f$center, u = f$axis, r = f$radius, p = f$pitch, phi0 = f$phi0))
        } else if (!is.null(f$H)) as.matrix(f$H) else NULL
      },
      NULL), error = function(e) NULL)
    if (is.null(pts) || !nrow(pts)) return(NULL)
    pts <- as.matrix(pts)[, 1:3, drop = FALSE]; colnames(pts) <- c("x", "y", "z"); pts
  }

  # Surface overlays (sphere/ellipsoid): a (theta,phi) vertex grid + triangle
  # faces for a translucent mesh3d in 3D; the projected convex hull is the exact
  # orthographic silhouette outline for these convex bodies in 2D.
  .viewer_surface_mesh <- function(f, shape, nu = 48L, nv = 28L) {
    res <- tryCatch({
      uu <- rep(seq(0, 2*pi, length.out = nu), times = nv)
      vv <- rep(seq(0, pi,  length.out = nv), each  = nu)
      base <- cbind(sin(vv) * cos(uu), sin(vv) * sin(uu), cos(vv))
      verts <- switch(shape,
        sphere    = base * f$radius,
        ellipsoid = sweep(base, 2, f$axes, "*") %*% t(f$R),
        NULL)
      if (is.null(verts)) stop("unsupported")
      verts <- sweep(verts, 2, f$center, "+")
      iu <- rep(seq_len(nu - 1L), times = nv - 1L)
      iv <- rep(seq_len(nv - 1L), each  = nu - 1L)
      a <- (iv - 1L) * nu + (iu - 1L); b <- (iv - 1L) * nu + iu     # 0-based vertex ids
      cc <- iv * nu + iu;              d <- iv * nu + (iu - 1L)
      list(verts = verts, i = c(a, a), j = c(b, cc), k = c(cc, d))
    }, error = function(e) NULL)
    res
  }

  output$viewer_header <- renderText({
    cur  <- viewer_current()
    vf   <- viewer_fit()
    best <- if (!is.null(vf)) vf$best else NA_character_
    rb   <- if (!is.null(vf) && !is.na(best)) as.numeric(vf$rmsd[[best]]) else NA_real_
    crit <- toupper(if (!is.null(vf)) vf$criterion else "aicc")
    wtxt <- if (!is.null(vf) && isTRUE(vf$weighted))
      sprintf("   |   photon-weighted (n_eff %.1f of %d)", vf$n_eff, vf$n) else ""
    unit_lbl <- if (identical(cur$unit, "cluster")) "de-blinked molecules" else "trace IDs"
    sprintf("Cluster %d / %d   |   %s   |   id %d   |   %d %s%s%s",
            cur$idx, cur$n, cur$file %||% "?", cur$id, nrow(cur$mat), unit_lbl,
            if (is.finite(rb)) sprintf("   |   best MLE fit (%s): %s (%sRMSD %.1f nm)", crit, best,
                                       if (!is.null(vf) && isTRUE(vf$weighted)) "weighted " else "", rb) else "",
            wtxt)
  })

  output$viewer_rmsd_table <- renderTable({
    vf <- viewer_fit()
    shiny::validate(shiny::need(!is.null(vf), "Need >= 4 traces for shape fitting."))
    M      <- vf$M
    shapes <- rownames(M)
    fnum   <- function(v) ifelse(is.finite(v), round(v, 1), NA_real_)
    data.frame(
      Shape        = shapes,
      params       = as.integer(M[, "params"]),
      `RMSD (nm)`  = round(M[, "rmsd"], 2),
      AICc         = fnum(M[, "aicc"]),
      BIC          = fnum(M[, "bic"]),
      Best         = ifelse(shapes == vf$best, "<-- best", ""),
      check.names  = FALSE
    )
  }, striped = TRUE, spacing = "xs")

  output$viewer_plot <- renderPlotly({
    cur   <- viewer_current()
    mat   <- cur$mat
    cloud <- cur$cloud
    # Series name follows the unit, so the legend cannot claim "trace medians"
    # while showing de-blinked molecule centroids.
    unit_nm <- if (identical(cur$unit, "cluster")) "De-blinked molecules" else "Trace medians"
    mode  <- input$viewer_mode %||% "3d"

    vf    <- viewer_fit()
    sel   <- input$viewer_fit_shape %||% "best"
    shape <- if (sel == "none") NA_character_
             else if (sel == "best") { if (!is.null(vf)) vf$best else NA_character_ }
             else sel
    f          <- if (!is.na(shape) && !is.null(vf)) vf$fits[[shape]] else NULL
    is_surface <- !is.na(shape) && shape %in% c("sphere", "ellipsoid")
    is_curve   <- !is.na(shape) && shape %in% c("ring", "arc", "helix")
    fit_col    <- "#D7301F"

    if (identical(mode, "pca")) {
      shiny::validate(shiny::need(nrow(mat) >= 3, "Need >= 3 units for a PCA projection."))
      # PCA basis from the UNIT centroids so unit length (localization count)
      # does not weight the axes; the full cloud is then projected onto that basis.
      ctr  <- colMeans(mat)
      pca  <- stats::prcomp(sweep(mat, 2, ctr), center = FALSE, scale. = FALSE)
      proj <- function(P) sweep(P, 2, ctr) %*% pca$rotation
      cl   <- proj(cloud); md <- proj(mat)
      fig  <- plotly::plot_ly() %>%
        plotly::add_markers(x = cl[, 1], y = cl[, 2],
          marker = list(color = "#2C7FB8", size = 3, opacity = 0.30),
          name = "Localizations") %>%
        plotly::add_markers(x = md[, 1], y = md[, 2],
          marker = list(color = "#08306B", size = 7, line = list(color = "white", width = 0.5)),
          name = unit_nm)
      if (!is.null(f) && is_curve) {
        cp <- .viewer_curve_points(f, shape)
        if (!is.null(cp)) { ps <- proj(cp)
          fig <- fig %>% plotly::add_trace(x = ps[, 1], y = ps[, 2],
            type = "scatter", mode = "lines",
            line = list(color = fit_col, width = 3), name = paste("Fit:", shape)) }
      } else if (!is.null(f) && is_surface) {
        mh <- .viewer_surface_mesh(f, shape)
        if (!is.null(mh)) { pv <- proj(mh$verts); h <- grDevices::chull(pv)
          poly <- pv[c(h, h[1]), , drop = FALSE]                # closed silhouette
          fig <- fig %>% plotly::add_trace(x = poly[, 1], y = poly[, 2],
            type = "scatter", mode = "lines", fill = "toself",
            fillcolor = "rgba(215,48,31,0.12)",
            line = list(color = fit_col, width = 3), name = paste("Fit:", shape)) }
      }
      fig %>% plotly::layout(
        xaxis  = list(title = "PC1 (nm)", scaleanchor = "y", scaleratio = 1),
        yaxis  = list(title = "PC2 (nm)"),
        legend = list(orientation = "h"))
    } else {
      fig <- plotly::plot_ly() %>%
        plotly::add_trace(x = cloud[, 1], y = cloud[, 2], z = cloud[, 3],
          type = "scatter3d", mode = "markers",
          marker = list(color = "#2C7FB8", size = 2, opacity = 0.30), name = "Localizations") %>%
        plotly::add_trace(x = mat[, 1], y = mat[, 2], z = mat[, 3],
          type = "scatter3d", mode = "markers",
          marker = list(color = "#08306B", size = 4), name = unit_nm)
      if (!is.null(f) && is_curve) {
        cp <- .viewer_curve_points(f, shape)
        if (!is.null(cp))
          fig <- fig %>% plotly::add_trace(x = cp[, 1], y = cp[, 2], z = cp[, 3],
            type = "scatter3d", mode = "lines",
            line = list(color = fit_col, width = 5), name = paste("Fit:", shape))
      } else if (!is.null(f) && is_surface) {
        mh <- .viewer_surface_mesh(f, shape)
        if (!is.null(mh))
          fig <- fig %>% plotly::add_trace(type = "mesh3d",
            x = mh$verts[, 1], y = mh$verts[, 2], z = mh$verts[, 3],
            i = mh$i, j = mh$j, k = mh$k,
            facecolor = rep(fit_col, length(mh$i)), opacity = 0.20,
            flatshading = FALSE, hoverinfo = "skip",
            name = paste("Fit:", shape), showlegend = TRUE)
      }
      fig %>% plotly::layout(
        scene = list(xaxis = list(title = "x (nm)"),
                     yaxis = list(title = "y (nm)"),
                     zaxis = list(title = "z (nm)"),
                     aspectmode = "data"),
        legend = list(orientation = "h"))
    }
  })

  ## ── Radius of Gyration ────────────────────────────────────────────────────
  rog_results <- reactiveVal(NULL)
  rog_params  <- reactiveVal(NULL)   # compute params, frozen at Run

  observeEvent(input$run_rog, {
    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }
    # mode drives BOTH what we group by and what a "point" is:
    #   "cluster2"  structure size - group by median cluster, one point per
    #               de-blinked MOLECULE when de-blinking is on, else one point per
    #               trace ID. (The auto fallback is the whole point of this mode.)
    #   "cluster"   label footprint - group by de-blinked molecule, one point per trace.
    mode        <- input$rog_cluster_col %||% "cluster2"
    cluster_col <- if (mode == "cluster") "cluster" else "cluster2"
    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }
    all_res <- list()
    withProgress(message = "Computing RoG...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
        if (is.null(df) || !nrow(df) || !cluster_col %in% names(df)) next

        # Split AFTER preparation, so the channels come from resolve_channels()
        # and the split covers a 2-colour DCR split as well as a multiplex build.
        # A cluster belongs to exactly one channel by construction, so this only
        # partitions the results - it never recomputes them differently.
        parts <- if (isTRUE(input$rog_by_channel)) channel_parts(df, fname)
                 else list(list(df = df, tag = fname))

        for (pt in parts) {
          # Only the plain structure-size mode collapses to one point per de-blinked
          # molecule (when de-blinking is on); the other two use raw traces as points.
          member_col <- if (mode == "cluster2" && deblink_on() &&
                            "cluster" %in% names(pt$df) && any(pt$df$cluster > 0, na.rm = TRUE))
                          "cluster" else "tid"
          r <- tryCatch(compute_cluster_rog(pt$df, cluster_col = cluster_col,
                                            member_col = member_col),
                        error = function(e) NULL)
          if (is.null(r) || !nrow(r)) next
          r$dataset <- pt$tag
          all_res[[pt$tag]] <- r
        }
      }
    })
    res <- do.call(rbind, Filter(Negate(is.null), all_res))

    if (is.null(res) || !nrow(res)) {
      hint <- if (cluster_col == "cluster")
        sprintf("Enable %s and press Plot for each file.", deblink_label())
      else
        sprintf("Enable %s and press Plot for each file.", cluster2_label())
      showNotification(paste0("No RoG computed. ", hint), type = "warning")
      rog_results(NULL); return()
    }
    rog_results(res)
    # Compute parameters frozen at Run: reading them live relabelled an existing
    # plot without recomputing it, so the axis could describe a different unit
    # than the data on screen.
    rog_params(list(mode = mode))
    showNotification(
      sprintf("RoG computed for %d cluster(s) across %d file(s).",
              nrow(res), length(unique(res$dataset))),
      type = "message", duration = 3)
  })

  build_rog_hist_fig <- function() {
    res <- rog_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run RoG first."))
    res <- apply_dataset_groups(res)
    shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))

    density_mode <- identical(input$rog_yaxis, "density")
    nbins     <- max(5L, as.integer(input$rog_bins %||% 30L))
    aggregate <- isTRUE(input$rog_aggregate)
    med_v     <- stats::median(res$rog, na.rm = TRUE)

    fig    <- plotly::plot_ly()
    shapes <- list()

    if (aggregate) {
      col <- "rgba(31,119,180,0.65)"
      fig <- fig %>% plotly::add_histogram(
        x        = res$rog,
        name     = sprintf("All datasets (N clusters = %d)", nrow(res)),
        nbinsx   = nbins,
        histnorm = if (density_mode) "probability density" else "",
        marker   = list(color = col, line = list(color = "rgba(31,119,180,1)", width = 1)),
        opacity  = 0.75
      )
      shapes <- list(list(
        type = "line", x0 = med_v, x1 = med_v, y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "firebrick", width = 2, dash = "dash")
      ))
    } else {
      datasets <- if ("dataset" %in% names(res)) unique(res$dataset) else NULL
      n_ds     <- length(datasets)
      pal <- if (n_ds >= 1) series_palette(datasets) else NULL

      if (!is.null(pal)) {
        for (ds in datasets) {
          sub     <- res[res$dataset == ds, , drop = FALSE]
          col_hex <- pal[ds]
          med_ds  <- stats::median(sub$rog, na.rm = TRUE)
          fig <- fig %>% plotly::add_histogram(
            x        = sub$rog,
            name     = ds,
            nbinsx   = nbins,
            histnorm = if (density_mode) "probability density" else "",
            marker   = list(color = paste0(col_hex, "99"),
                            line  = list(color = col_hex, width = 0.5)),
            opacity  = 0.75
          )
          shapes <- c(shapes, list(list(
            type = "line", x0 = med_ds, x1 = med_ds, y0 = 0, y1 = 1, yref = "paper",
            line = list(color = col_hex, width = 2, dash = "dash"), opacity = 0.9
          )))
        }
      } else {
        fig <- fig %>% plotly::add_histogram(
          x        = res$rog,
          nbinsx   = nbins,
          histnorm = if (density_mode) "probability density" else "",
          marker   = list(color = "rgba(31,119,180,0.65)",
                          line  = list(color = "rgba(31,119,180,1)", width = 1)),
          name     = "RoG"
        )
        shapes <- list(list(
          type = "line", x0 = med_v, x1 = med_v, y0 = 0, y1 = 1, yref = "paper",
          line = list(color = "firebrick", width = 2, dash = "dash")
        ))
      }
    }

    fig <- fig %>%
      plotly::layout(
        barmode = "overlay",
        xaxis = list(title = if (identical((rog_params()$mode) %||% (input$rog_cluster_col %||% "cluster2"), "cluster"))
                               "Label footprint Rg (nm)" else "Median-cluster Rg - structure size (nm)"),
        yaxis = list(title = if (density_mode) "Density" else "Count"),
        shapes = shapes,
        annotations = list(
          list(x = med_v, y = 1, yref = "paper",
               text = sprintf(" overall median: %.1f nm", med_v),
               showarrow = FALSE, xanchor = "left",
               font = list(color = "firebrick", size = 11))
        ),
        legend = list(orientation = "h", y = -0.18),
        paper_bgcolor = "white", plot_bgcolor = "white"
      )
    fig
  }
  output$rog_hist <- renderPlotly(build_rog_hist_fig())

  output$rog_small_n <- renderUI({
    res <- rog_results()
    if (is.null(res) || !nrow(res)) return(NULL)
    nv <- if ("n_units" %in% names(res)) res$n_units else
          if ("n_traces" %in% names(res)) res$n_traces else NULL
    if (is.null(nv)) return(NULL)
    .small_n_note(nv, "structure")
  })

  output$download_rog <- downloadHandler(
    filename = function() paste0("rog_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- rog_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  ## ── Ripley's L(r) ─────────────────────────────────────────────────────────
  ## ── Global (whole-curve) CSR test: band + caption ──────────────────────────
  ## The shaded envelope is POINTWISE: it holds at each r on its own, so a curve
  ## that leaves it somewhere is not significant at 5% - on simulated random data
  ## in a cell footprint that happens in 96-100% of runs. `g_lo`/`g_hi` (from
  ## compute_ripley_envelope) are the studentized maximum-deviation band, which the
  ## curve crosses anywhere with probability 5% under the null, and `global_p` is
  ## the matching p-value for the whole r range. Both are drawn/reported per file.
  ## `smooth_fn` is the builder's own smoother (each panel defines its own from its
  ## smoothing controls), so the band gets exactly the same treatment as the curve.
  .add_global_band <- function(fig, e, col, nm, smooth_fn = NULL) {
    if (is.null(e$g_lo) || all(!is.finite(e$g_lo))) return(fig)
    sm <- if (is.function(smooth_fn)) smooth_fn else function(r, y) y
    gl <- sm(e$r, e$g_lo); gh <- sm(e$r, e$g_hi)
    fig %>%
      plotly::add_trace(x = e$r, y = gh, type = "scatter", mode = "lines",
        line = list(color = col, width = 1, dash = "dash"), legendgroup = nm,
        name = paste0(nm, " global 95% band"), showlegend = TRUE, hoverinfo = "skip") %>%
      plotly::add_trace(x = e$r, y = gl, type = "scatter", mode = "lines",
        line = list(color = col, width = 1, dash = "dash"), legendgroup = nm,
        showlegend = FALSE, hoverinfo = "skip")
  }
  ## One line naming each file's global p, for the plot caption. Files beyond the
  ## fourth are summarised so the caption cannot grow without bound.
  .global_p_caption <- function(env_res) {
    if (is.null(env_res) || !"global_p" %in% names(env_res)) return(NULL)
    d <- unique(env_res[, c("dataset", "global_p")])
    d <- d[is.finite(d$global_p), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    fmt <- function(p) if (p <= 1 / 100) sprintf("p \u2264 %.2g", p) else sprintf("p = %.2f", p)
    txt <- paste(sprintf("%s %s", d$dataset, vapply(d$global_p, fmt, character(1))), collapse = "  \u00b7  ")
    if (nrow(d) > 4) txt <- paste0(paste(sprintf("%s %s", d$dataset[1:4], vapply(d$global_p[1:4], fmt, character(1))),
                                         collapse = "  \u00b7  "),
                                   sprintf("  \u00b7  (+%d more, %d with p < 0.05)", nrow(d) - 4,
                                           sum(d$global_p < 0.05, na.rm = TRUE)))
    paste("Global test over the whole r range:", txt)
  }

  ## Caption as a plot annotation (NULL when no envelope was computed).
  .global_p_annot <- function(env_res) {
    cap <- .global_p_caption(env_res)
    if (is.null(cap)) return(NULL)
    list(list(text = cap, xref = "paper", yref = "paper", x = 0, y = 1.06,
              xanchor = "left", showarrow = FALSE, font = list(size = 11, color = "#555")))
  }

  ripley_results     <- reactiveVal(NULL)
  ripley_params      <- reactiveVal(NULL)   # compute settings, frozen at Run
  ripley_env_results <- reactiveVal(NULL)  # Monte-Carlo CSR envelope (per dataset)

  observeEvent(input$run_ripley, {
    unit      <- resolve_unit(input$ripley_unit, "cluster", "tid", "Ripley's L")
    use_3d    <- isTRUE(input$ripley_use_3d)
    r_max_in  <- as.numeric(input$ripley_r_max  %||% 0)
    r_steps   <- max(10L, as.integer(input$ripley_r_steps %||% 60L))
    all_files <- isTRUE(input$ripley_all_files)
    want_env  <- isTRUE(input$ripley_envelope)
    nsim      <- max(19L, min(199L, as.integer(input$ripley_nsim %||% 99L)))
    env_file_cap <- 6L  # cap envelope work in all-files mode (nsim × files is O(n²)-heavy)

    needs_cluster <- unit == "cluster"

    get_pts_flat <- function(df) {
      if (unit == "tid") {
        df %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(x=median(x,na.rm=T), y=median(y,na.rm=T),
                           z=if("z"%in%names(df)) median(z,na.rm=T) else 0, .groups="drop") %>%
          dplyr::filter(is.finite(x), is.finite(y))
      } else {
        if (!"cluster" %in% names(df) || !any(!is.na(df$cluster) & df$cluster > 0L)) return(NULL)
        df %>%
          dplyr::filter(!is.na(cluster), cluster > 0L) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(x=median(x,na.rm=T), y=median(y,na.rm=T),
                           z=if("z"%in%names(df)) median(z,na.rm=T) else 0, .groups="drop") %>%
          dplyr::filter(is.finite(x), is.finite(y))
      }
    }

    compute_ripley_for_df <- function(df, tag, r_seq_use, do_env, close_nm = NULL,
                                      win_df = NULL, win_tag = NULL) {
      pts <- get_pts_flat(df)
      if (is.null(pts) || nrow(pts) < 2L) return(NULL)
      # THE WINDOW COMES FROM win_df, NOT FROM pts. Normally they are the same
      # frame. They differ when this file has been split into channels: there
      # win_df is the pooled file, so both channels are measured in one window
      # (see split_df_by_channel). Deriving it from the channel's own points
      # would shrink the window around the sparser channel and deflate its curve.
      wdf  <- if (is.data.frame(win_df) && nrow(win_df)) win_df else df
      wpts <- if (identical(wdf, df)) pts else get_pts_flat(wdf)
      if (is.null(wpts) || !nrow(wpts)) wpts <- pts
      # Freeze the observation window so the curve and its CSR envelope share one
      # identical window (edge bias then cancels).
      axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
      win  <- apply(as.matrix(wpts[, axes, drop = FALSE]), 2, range)
      # Support window (opt-in): the observed curve and the envelope MUST get the
      # same support, or the band no longer references the same geometry. Keyed on
      # win_tag so channels of one file share a single cached support.
      # Same shared choice of outline points as the hotspot panel, so one cell
      # gets one area whichever panel measured it.
      sup  <- spatial_support_of(spatial_outline_pts(wdf, wpts), use_3d,
                                 tag = win_tag %||% tag, close_nm = close_nm)
      # One set-covariance build per file, reused by the curve AND every simulation.
      if (!is.null(sup) && isTRUE(input$spatial_edge_correct))
        sup <- attach_gamma_cache(sup, r_max = sqrt(sum((win[2, ] - win[1, ])^2)),
                                  use_3d = use_3d)
      curve <- compute_ripley_L(pts$x, pts$y, if (use_3d) pts$z else NULL,
                                r_seq = r_seq_use, use_3d = use_3d,
                                window_ranges = win, support = sup,
                                edge_correct = isTRUE(input$spatial_edge_correct))
      if (!nrow(curve)) return(NULL)
      curve$dataset <- tag
      env <- NULL
      if (do_env) {
        env <- compute_ripley_envelope(n = min(nrow(pts), 2000L),
                                       window_ranges = win, r_seq = curve$r,
                                       use_3d = use_3d, nsim = nsim, estimator = "L",
                                       support = sup,
                                       edge_correct = isTRUE(input$spatial_edge_correct),
                                       progress = env_prog)
        # Global (studentized max-deviation) test, computed here because rbind drops
        # the attributes the envelope carries. One p per file for the whole r range.
        if (!is.null(env)) {
          gp <- envelope_global_p(env, curve$Lr_minus_r)
          env$global_p <- gp$p %||% NA_real_
          env$dataset <- tag
        }
      }
      list(curve = curve, env = env)
    }

    get_df_list <- function() {
      if (all_files) {
        store <- smlm_file_store()
        shiny::validate(shiny::need(length(store) > 0, "No files loaded."))
        # Prepare each FILE once, then split it into channels if asked. Splitting
        # after preparation (rather than handing the pipeline pre-split raw data,
        # as this used to) is what lets both channels keep the pooled file as
        # their shared observation window - and it means the split now covers a
        # 2-colour DCR split as well as a multiplex build.
        items <- lapply(names(store), function(nm) {
          # Cached chain. NOTE this also brings run_ripley in line with every other
          # panel by applying qc_exclude_from, which it previously skipped (a no-op
          # unless 'Exclude flagged traces' is on).
          df <- prepared_df(nm, needs_cluster = needs_cluster)
          if (is.null(df)) df <- store[[nm]][0, , drop = FALSE]
          list(df = df, tag = nm)
        })
        if (isTRUE(input$ripley_by_channel))
          items <- unlist(lapply(items, function(it) split_df_by_channel(it$df, it$tag)),
                          recursive = FALSE)
        else
          items <- lapply(items, function(it) c(it, list(win_df = it$df, win_tag = it$tag)))
        items
      } else {
        df <- plottable_data()
        shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "Press Plot first."))
        nm <- input$active_file %||% "Active file"
        if (isTRUE(input$ripley_by_channel)) split_df_by_channel(df, nm)
        else list(list(df = df, tag = nm, win_df = df, win_tag = nm))
      }
    }

    df_list <- get_df_list()

    # One gap-fill radius for the whole set (see spatial_resolve_close()): a
    # per-file radius makes the windows - and therefore K, which scales with the
    # window area - incomparable between files. Never blocks.
    # The sweep runs on the same points the window will be drawn from.
    cc <- spatial_gate_close(lapply(df_list, function(it) {
      u <- tryCatch(get_pts_flat(it$df), error = function(e) NULL)
      list(df = spatial_outline_pts(it$df, u), tag = it$tag)
    }), use_3d)
    common_close <- cc$close_nm

    if (want_env && length(df_list) > env_file_cap) {
      showNotification(sprintf(
        "CSR envelope computed for the first %d of %d files only (cost is nsim × files).",
        env_file_cap, length(df_list)), type = "message")
    }

    # The envelope dominates the runtime, so the bar advances per SIMULATION.
    n_env_files <- if (want_env) min(length(df_list), env_file_cap) else 0L
    env_prog <- if (want_env)
      function(frac) incProgress(frac / max(1L, n_env_files)) else NULL

    run_all <- function() {
      lapply(seq_along(df_list), function(i) {
        item      <- df_list[[i]]
        r_seq_use <- if (r_max_in > 0) seq(0, r_max_in, length.out = r_steps + 1L)[-1] else NULL
        do_env    <- want_env && i <= env_file_cap
        compute_ripley_for_df(item$df, item$tag, r_seq_use, do_env, common_close,
                              win_df = item$win_df, win_tag = item$win_tag)
      })
    }
    all_out <- if (want_env) {
      withProgress(message = sprintf("Ripley L: %d CSR simulations x %d file(s)...",
                                     nsim, n_env_files), value = 0, run_all())
    } else {
      # Without an envelope there is no per-simulation tick, but building the
      # windows still takes seconds - show a bar rather than appear frozen.
      withProgress(message = sprintf("Ripley L: %d file(s)...", length(df_list)),
                   value = 0.5, run_all())
    }
    spatial_report_windows(vapply(df_list, function(it) as.character(it$tag), character(1)))

    all_out    <- Filter(Negate(is.null), all_out)
    all_curves <- do.call(rbind, lapply(all_out, `[[`, "curve"))
    shiny::validate(shiny::need(!is.null(all_curves) && nrow(all_curves) > 0, "Ripley computation returned no data."))
    ripley_results(all_curves)
    ripley_params(list(unit = input$ripley_unit, `r max` = input$ripley_r_max,
                       `r steps` = input$ripley_r_steps, `3-D` = use_3d,
                       `all files` = all_files, envelope = want_env,
                       simulations = input$ripley_nsim,
                       `edge correction` = isTRUE(input$spatial_edge_correct),
                       `support window` = isTRUE(input$spatial_use_support)))

    envs <- Filter(Negate(is.null), lapply(all_out, `[[`, "env"))
    ripley_env_results(if (length(envs)) do.call(rbind, envs) else NULL)

  })

  build_ripley_plot_fig <- function() {
    res <- ripley_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run Ripley's L first."))

    do_smooth <- isTRUE(input$ripley_smooth)
    span_val  <- as.numeric(input$ripley_smooth_span %||% 0.2)
    smooth_vec <- function(r, y) {
      if (!do_smooth || length(y) < 5L) return(y)
      tryCatch(stats::predict(stats::loess(y ~ r, data = data.frame(r = r, y = y), span = span_val)),
               error = function(e) y)
    }
    hex2rgba <- function(hex, a) {
      cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
    }
    grouped <- (input$analysis_by %||% "file") == "group"
    env_res <- ripley_env_results()
    # Declared up front because the axis code below reads them in both branches.
    rmin_c <- NA_real_; rmax_c <- NA_real_
    show_files <- isTRUE(input$ripley_show_files)
    fig <- plotly::plot_ly()

    if (!grouped) {
      datasets <- unique(res$dataset)
      pal <- series_palette(datasets)
      for (nm in datasets) {
        sub <- res[res$dataset == nm, ]
        # CSR envelope band (drawn first so the observed curve sits on top)
        if (!is.null(env_res) && nm %in% env_res$dataset) {
          # The BAND MUST GET THE SAME SMOOTHING AS THE CURVE. Comparing a smoothed
          # observed curve against a raw envelope compares two different statistics:
          # smoothing damps the curve's excursions but left the band at its full
          # unsmoothed width, so the comparison silently became conservative (and the
          # band looked jagged next to a smooth curve). Smoothing the quantiles is a
          # close approximation to taking quantiles of smoothed simulations - the two
          # differ only in how a monotone smoother interacts with the quantile, which
          # is negligible here and keeps the band consistent with whatever smoothing
          # is displayed (smoothing is a live display control, the envelope is not).
          e <- env_res[env_res$dataset == nm, ]
          e <- e[order(e$r), ]
          e_hi <- smooth_vec(e$r, e$hi); e_lo <- smooth_vec(e$r, e$lo)
          e_nm <- smooth_vec(e$r, e$null_mean)
          fig <- fig %>%
            plotly::add_trace(x = e$r, y = e_hi, type = "scatter", mode = "lines",
              line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
            plotly::add_trace(x = e$r, y = e_lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
              legendgroup = nm, name = paste0(nm, " CSR 95% envelope"),
              showlegend = TRUE, hoverinfo = "skip") %>%
            plotly::add_trace(x = e$r, y = e_nm, type = "scatter", mode = "lines",
              line = list(color = unname(pal[nm]), width = 1, dash = "dot"),
              legendgroup = nm, showlegend = FALSE, hoverinfo = "skip")
          fig <- .add_global_band(fig, e, unname(pal[nm]), nm, smooth_vec)
        }
        fig <- fig %>% plotly::add_trace(
          x = sub$r, y = smooth_vec(sub$r, sub$Lr_minus_r),
          type = "scatter", mode = "lines", name = nm, legendgroup = nm,
          line = list(color = unname(pal[nm]), width = 2), hoverinfo = "skip")
      }
    } else {
      gm <- file_group_map()
      shiny::validate(shiny::need(length(gm) > 0, "No files assigned to groups (group mode)."))
      res <- res[res$dataset %in% names(gm), , drop = FALSE]
      shiny::validate(shiny::need(nrow(res) > 0, "No group files in the Ripley results."))
      res$group <- unname(gm[res$dataset])
      pal  <- series_palette(unique(res$group))
      fpal <- file_palette(unique(res$dataset))
      band <- input$ripley_group_band %||% "iqr"

      # common r-grid up to the smallest per-file r max (avoid extrapolation),
      # then average each group's per-file curves + between-file spread band.
      file_maxr <- tapply(res$r, res$dataset, function(v) max(v, na.rm = TRUE))
      rmax_c <- suppressWarnings(min(file_maxr, na.rm = TRUE))
      # Common grid spans only the range EVERY file actually measured. It used to
      # start at r = 0 while the stored curves start at r_max/steps, so approx(rule=2)
      # flat-extrapolated the first real value back to r = 0 - a false non-zero point
      # at a radius where L(r)-r is 0 by definition, injected into every group curve
      # and its band.
      file_minr <- tapply(res$r, res$dataset, function(v) min(v, na.rm = TRUE))
      rmin_c <- suppressWarnings(max(file_minr, na.rm = TRUE))
      if (!is.finite(rmin_c) || rmin_c < 0) rmin_c <- 0
      # Grid density from the STORED curve, not the live control: r steps is a
      # compute parameter and only takes effect on Run.
      steps  <- max(10L, length(unique(res$r)))
      shiny::validate(shiny::need(is.finite(rmax_c) && rmax_c > 0, "Could not form a common r-grid."))
      rr <- seq(rmin_c, rmax_c, length.out = steps)

      for (g in unique(res$group)) {
        files <- unique(res$dataset[res$group == g])
        mat <- vapply(files, function(f) {
          d <- res[res$dataset == f, ]
          stats::approx(d$r, d$Lr_minus_r, xout = rr, rule = 2)$y
        }, numeric(length(rr)))
        if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
        meanc <- smooth_vec(rr, rowMeans(mat, na.rm = TRUE))
        col   <- unname(pal[g]); nfl <- length(files)

        # ── BAND MUST MATCH THE DRAWN MEAN ──────────────────────────────────
        # Three faults were fixed here together:
        #  1. the mean was smoothed and the band was not, so a smooth line sat in a
        #     jagged band (measured: band edge 16-18x the roughness of the mean).
        #     Same argument as the CSR envelope above - a band and a curve that got
        #     different treatment are not comparable.
        #  2. the SD band was built around the RAW rowMeans while the SMOOTHED mean
        #     was drawn, so the line was not centred in its own band (up to 22 nm
        #     off). It is now centred on the curve actually shown.
        #  3. quartiles need >= 4 values. On 2 files the IQR band is the middle half
        #     of a 2-point range and contains NEITHER file (measured 0% coverage),
        #     which draws a reassuringly tight band over maximal disagreement. Under
        #     4 files the full min-max range is drawn instead, and the legend says so.
        # ── CSR ENVELOPE IN GROUP MODE ──────────────────────────────────────
        # This used to be omitted entirely: envelopes were computed per file but
        # the grouped branch never drew them, so ticking "Show CSR envelope" did
        # nothing as soon as you switched to groups - and the between-file spread
        # band was easy to mistake for it.
        #
        # The group curve is the mean of its files' curves, so its CSR reference
        # is the mean of those files' envelopes, interpolated onto the same grid
        # and smoothed the same way. Drawn as a dashed outline rather than a fill,
        # so it cannot be confused with the solid spread band: they answer
        # different questions - the spread band is "how much do my cells differ",
        # the envelope is "what would randomness look like".
        env_n <- 0L
        if (!is.null(env_res)) {
          ef <- intersect(files, unique(env_res$dataset))
          if (length(ef)) {
            elo <- vapply(ef, function(f) {
              e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
              stats::approx(e$r, e$lo, xout = rr, rule = 2)$y }, numeric(length(rr)))
            ehi <- vapply(ef, function(f) {
              e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
              stats::approx(e$r, e$hi, xout = rr, rule = 2)$y }, numeric(length(rr)))
            if (is.null(dim(elo))) { elo <- matrix(elo, ncol = 1); ehi <- matrix(ehi, ncol = 1) }
            env_n <- length(ef)
            e_lo <- smooth_vec(rr, rowMeans(elo, na.rm = TRUE))
            e_hi <- smooth_vec(rr, rowMeans(ehi, na.rm = TRUE))
            # TWO lines per group - the upper and lower bound of the same band -
            # so N groups draw 2N dotted lines. They are given their own legend
            # entry and their own legendgroup, so it is obvious what they are and
            # they can be switched off independently of the group curve.
            fig <- fig %>% plotly::add_trace(
              x = rr, y = e_hi, type = "scatter", mode = "lines",
              line = list(color = col, width = 1, dash = "dot"),
              name = sprintf("%s – CSR 95%% band", g),
              legendgroup = paste0(g, "__csr"), showlegend = TRUE, hoverinfo = "skip")
            fig <- fig %>% plotly::add_trace(
              x = rr, y = e_lo, type = "scatter", mode = "lines",
              line = list(color = col, width = 1, dash = "dot"),
              legendgroup = paste0(g, "__csr"), showlegend = FALSE, hoverinfo = "skip")
          }
        }

        band_lbl <- ""
        if (band != "none" && nfl >= 2) {
          if (band == "iqr" && nfl >= 4L) {
            lo <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.25, na.rm = TRUE))
            hi <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.75, na.rm = TRUE))
            band_lbl <- ", IQR"
          } else if (band == "iqr") {
            lo <- smooth_vec(rr, apply(mat, 1, min, na.rm = TRUE))
            hi <- smooth_vec(rr, apply(mat, 1, max, na.rm = TRUE))
            band_lbl <- ", range"
          } else {
            sdv <- smooth_vec(rr, apply(mat, 1, stats::sd, na.rm = TRUE))
            lo <- meanc - sdv; hi <- meanc + sdv
            band_lbl <- ", ±1 SD"
          }
          fig <- fig %>%
            plotly::add_trace(x = rr, y = hi, type = "scatter", mode = "lines",
              line = list(width = 0), legendgroup = g, showlegend = FALSE, hoverinfo = "skip") %>%
            plotly::add_trace(x = rr, y = lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = hex2rgba(col, 0.18), line = list(width = 0),
              legendgroup = g, showlegend = FALSE, hoverinfo = "skip")
        }
        # Individual files as shades of the group colour, drawn thin and behind the
        # group mean, so a group difference can be judged against its replicates.
        if (show_files && nfl >= 2) {
          for (f in files) {
            d <- res[res$dataset == f, ]; d <- d[order(d$r), ]
            fc <- unname(fpal[f]); if (is.na(fc)) fc <- col
            fig <- fig %>% plotly::add_trace(
              x = d$r, y = d$Lr_minus_r, type = "scatter", mode = "lines",
              name = f, legendgroup = g, showlegend = FALSE,
              line = list(color = fc, width = 1), opacity = 0.75,
              hoverinfo = "skip")
          }
        }
        # Say how many files the envelope averages. With more files than the
        # envelope cap it covers only some of them, and a reader has to know that
        # before reading the dotted lines as the group's null.
        env_lbl <- if (env_n > 0L && env_n < nfl) sprintf(", CSR n=%d", env_n)
                   else if (env_n > 0L) ", CSR" else ""
        fig <- fig %>% plotly::add_trace(x = rr, y = meanc, type = "scatter", mode = "lines",
          name = sprintf("%s (n=%d%s%s)", g, nfl, band_lbl, env_lbl), legendgroup = g,
          line = list(color = col, width = 2.5), hoverinfo = "skip")
      }
    }

    # ── The x extent that is actually DRAWN ───────────────────────────────
    # In group mode the mean and its band live on a common grid capped at the
    # SMALLEST per-file r max, so that no file is extrapolated. The reference line
    # used to be drawn over range(res$r) instead - the range across ALL files - so
    # one file with a larger r max stretched it far past every curve, plotly
    # autoranged to the line, and the plot ended in a wide empty strip. The
    # reference line now spans exactly what is plotted, and nothing else reaches
    # beyond it unless the individual per-file curves are switched on, which draw
    # at their own full range and legitimately extend the axis.
    x_drawn <- if (grouped) {
      c(rmin_c, if (isTRUE(show_files)) max(res$r, na.rm = TRUE) else rmax_c)
    } else range(res$r, na.rm = TRUE)

    csr_label <- if (!is.null(env_res) && !grouped)
      "theoretical CSR (uncorrected - see envelope)" else "CSR (L−r = 0)"
    fig %>%
      plotly::add_trace(
        x = x_drawn, y = c(0, 0), type = "scatter", mode = "lines",
        name = csr_label,
        line = list(color = "black", width = 1, dash = "dash"), showlegend = TRUE,
        hoverinfo = "skip") %>%
      plotly::layout(
        xaxis  = list(title = .rmax_axis_title(grouped, res, rmax_c, show_files)),
        yaxis  = list(title = "L(r) − r (nm)"),
        annotations = .global_p_annot(env_res),
        # Hover deliberately OFF on this plot: the curve is read as a shape against
        # the envelope, and a unified tooltip listing every series' L−r at each r
        # was noise that also covered the plot.
        hovermode = FALSE,
        legend = list(groupclick = "togglegroup"),
        paper_bgcolor = "white", plot_bgcolor  = "white")
  }
  output$ripley_plot <- renderPlotly(build_ripley_plot_fig())

  output$download_ripley <- downloadHandler(
    filename = function() paste0("ripley_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- ripley_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  ## ── Univariate pair-correlation g(r) ──────────────────────────────────────
  ## GOAL: the scale-LOCAL clustering statistic, and the recommended endpoint for
  ## comparing conditions whose labelling density differs.
  ##
  ## WHY IT IS THE RIGHT ENDPOINT HERE (see also the LABELLING-DENSITY block):
  ##   * g is normalised by the intensity SQUARED, so it is invariant to how
  ##     densely the sample happens to be labelled - verified invariant to within
  ##     1-2% under thinning down to 10% of the units.
  ##   * being a derivative rather than a cumulative statistic, an aggregate of
  ##     radius R shows up as a PEAK near R instead of the broad rising ramp L(r)
  ##     gives, and a short-r blinking/label artefact stays in the first few bins
  ##     instead of contaminating every larger r.
  ## Runs on de-blinked units, same as Ripley's L.
  pcf_results     <- reactiveVal(NULL)
  pcf_params      <- reactiveVal(NULL)   # compute settings, frozen at Run
  pcf_env_results <- reactiveVal(NULL)

  observeEvent(input$run_pcf, {
    unit      <- resolve_unit(input$pcf_unit, "cluster", "tid", "g(r)")
    use_3d    <- isTRUE(input$pcf_use_3d)
    r_max_in  <- as.numeric(input$pcf_r_max %||% 0)
    r_steps   <- max(10L, as.integer(input$pcf_r_steps %||% 50L))
    all_files <- isTRUE(input$pcf_all_files)
    want_env  <- isTRUE(input$pcf_envelope)
    nsim      <- max(19L, min(199L, as.integer(input$pcf_nsim %||% 99L)))
    env_file_cap <- 6L

    needs_cluster <- unit == "cluster"

    compute_for_df <- function(df, tag, r_seq_use, do_env, close_nm = NULL,
                               win_df = NULL, win_tag = NULL) {
      pts <- unit_centroids(df, unit = unit, use_3d = use_3d)
      if (is.null(pts) || nrow(pts) < 3L) return(NULL)
      # Window from win_df, not from pts - see compute_ripley_for_df(). When this
      # file is split into channels, win_df is the pooled file so both channels
      # are measured in one window; g(r) divides by density, so a per-channel
      # window would invent a difference out of the count imbalance alone.
      wdf  <- if (is.data.frame(win_df) && nrow(win_df)) win_df else df
      wpts <- if (identical(wdf, df)) pts else unit_centroids(wdf, unit = unit, use_3d = use_3d)
      if (is.null(wpts) || !nrow(wpts)) wpts <- pts
      # Freeze the window so the curve and its envelope share one window and the
      # (uncorrected) edge bias cancels between them.
      axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
      win  <- apply(as.matrix(wpts[, axes, drop = FALSE]), 2, range)
      # Same shared choice of outline points as the hotspot panel, so one cell
      # gets one area whichever panel measured it.
      sup  <- spatial_support_of(spatial_outline_pts(wdf, wpts), use_3d,
                                 tag = win_tag %||% tag, close_nm = close_nm)
      if (!is.null(sup) && isTRUE(input$spatial_edge_correct))
        sup <- attach_gamma_cache(sup, r_max = sqrt(sum((win[2, ] - win[1, ])^2)),
                                  use_3d = use_3d)
      curve <- compute_pcf(pts$x, pts$y, if (use_3d) pts$z else NULL,
                           r_seq = r_seq_use, use_3d = use_3d, window_ranges = win,
                           support = sup,
                           edge_correct = isTRUE(input$spatial_edge_correct))
      if (is.null(curve) || !nrow(curve)) return(NULL)
      curve$dataset <- tag
      curve$n_units <- nrow(pts)
      env <- NULL
      if (do_env) {
        env <- compute_ripley_envelope(n = min(nrow(pts), 2000L), window_ranges = win,
                                       r_seq = curve$r, use_3d = use_3d,
                                       nsim = nsim, estimator = "g", support = sup,
                                       edge_correct = isTRUE(input$spatial_edge_correct),
                                       progress = env_prog)
        if (!is.null(env)) {
          gp <- envelope_global_p(env, curve$g)
          env$global_p <- gp$p %||% NA_real_
          env$dataset <- tag
        }
      }
      list(curve = curve, env = env)
    }

    get_df_list <- function() {
      if (all_files) {
        store <- smlm_file_store()
        shiny::validate(shiny::need(length(store) > 0, "No files loaded."))
        # Prepare each FILE once, then split into channels - so both channels
        # share the pooled file as their observation window.
        items <- lapply(names(store), function(nm) {
          df <- prepared_df(nm, needs_cluster = needs_cluster)
          if (is.null(df)) df <- store[[nm]][0, , drop = FALSE]
          list(df = df, tag = nm)
        })
        if (isTRUE(input$pcf_by_channel))
          unlist(lapply(items, function(it) split_df_by_channel(it$df, it$tag)),
                 recursive = FALSE)
        else
          lapply(items, function(it) c(it, list(win_df = it$df, win_tag = it$tag)))
      } else {
        df <- plottable_data()
        shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "Press Plot first."))
        nm <- input$active_file %||% "Active file"
        if (isTRUE(input$pcf_by_channel)) split_df_by_channel(df, nm)
        else list(list(df = df, tag = nm, win_df = df, win_tag = nm))
      }
    }

    df_list <- get_df_list()

    # One gap-fill radius for the whole set (see spatial_resolve_close()), swept on
    # the same points the window will be drawn from.
    cc <- spatial_gate_close(lapply(df_list, function(it) {
      u <- tryCatch(unit_centroids(it$df, unit = unit, use_3d = use_3d),
                    error = function(e) NULL)
      list(df = spatial_outline_pts(it$df, u), tag = it$tag)
    }), use_3d)
    common_close <- cc$close_nm

    if (want_env && length(df_list) > env_file_cap)
      showNotification(sprintf(
        "CSR envelope computed for the first %d of %d files only (cost is simulations x files).",
        env_file_cap, length(df_list)), type = "message")

    n_env_files <- if (want_env) min(length(df_list), env_file_cap) else 0L
    env_prog <- if (want_env)
      function(frac) incProgress(frac / max(1L, n_env_files)) else NULL

    run_all <- function() lapply(seq_along(df_list), function(i) {
      r_seq_use <- if (r_max_in > 0) seq(0, r_max_in, length.out = r_steps) else NULL
      compute_for_df(df_list[[i]]$df, df_list[[i]]$tag, r_seq_use,
                     want_env && i <= env_file_cap, common_close,
                     win_df = df_list[[i]]$win_df, win_tag = df_list[[i]]$win_tag)
    })
    all_out <- if (want_env) withProgress(
                 message = sprintf("g(r): %d CSR simulations x %d file(s)...", nsim, n_env_files),
                 value = 0, run_all())
               else withProgress(message = sprintf("g(r): %d file(s)...", length(df_list)),
                                 value = 0.5, run_all())
    spatial_report_windows(vapply(df_list, function(it) as.character(it$tag), character(1)))

    all_out <- Filter(Negate(is.null), all_out)
    curves  <- do.call(rbind, lapply(all_out, `[[`, "curve"))
    shiny::validate(shiny::need(!is.null(curves) && nrow(curves) > 0,
                               "g(r) returned no data (need >=3 units per file)."))
    pcf_results(curves)
    pcf_params(list(unit = input$pcf_unit, `r max` = input$pcf_r_max,
                    `r steps` = input$pcf_r_steps, `3-D` = use_3d,
                    `all files` = all_files, envelope = want_env,
                    simulations = input$pcf_nsim,
                    `edge correction` = isTRUE(input$spatial_edge_correct),
                    `support window` = isTRUE(input$spatial_use_support)))
    envs <- Filter(Negate(is.null), lapply(all_out, `[[`, "env"))
    pcf_env_results(if (length(envs)) do.call(rbind, envs) else NULL)
  })

  # ── Correlation length xi, one row per curve ──────────────────────────────
  # Fits each curve the plot actually draws: per FILE in file mode, and per GROUP
  # in group mode - where the group curve is the mean of its files' curves on a
  # common r-grid, exactly as plotted, so the number matches the line on screen.
  #
  # Per-file rows are kept in group mode too. A single xi per group hides whether
  # the group is tight or whether one cell is carrying it, and that spread is the
  # thing a test between conditions has to work with.
  pcf_xi_results <- reactive({
    req(isTRUE(input$pcf_fit_xi))
    res <- pcf_results()
    if (is.null(res) || !nrow(res)) return(NULL)
    r_min <- suppressWarnings(as.numeric(input$pcf_fit_rmin))
    if (!length(r_min) || !is.finite(r_min) || r_min < 0) r_min <- 0

    per_file <- do.call(rbind, lapply(split(res, res$dataset), function(d) {
      d <- d[order(d$r), ]
      f <- fit_correlation_length(d$r, d$g, r_min = r_min)
      cbind(data.frame(curve = d$dataset[1], level = "file",
                       stringsAsFactors = FALSE), f)
    }))

    out <- per_file
    if ((input$analysis_by %||% "file") == "group") {
      gm <- file_group_map()
      if (length(gm)) {
        res$group <- unname(gm[res$dataset])
        res <- res[!is.na(res$group), , drop = FALSE]
        if (nrow(res)) {
          # Same common-grid averaging the plot uses - never pool coordinates.
          per_grp <- do.call(rbind, lapply(split(res, res$group), function(d) {
            files <- unique(d$dataset)
            rmin_c <- max(tapply(d$r, d$dataset, min), na.rm = TRUE)
            rmax_c <- min(tapply(d$r, d$dataset, max), na.rm = TRUE)
            if (!is.finite(rmin_c) || !is.finite(rmax_c) || rmax_c <= rmin_c) return(NULL)
            rr <- seq(rmin_c, rmax_c, length.out = max(10L, length(unique(d$r))))
            mat <- vapply(files, function(f) {
              dd <- d[d$dataset == f, ]; dd <- dd[order(dd$r), ]
              stats::approx(dd$r, dd$g, xout = rr, rule = 2)$y
            }, numeric(length(rr)))
            if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
            f <- fit_correlation_length(rr, rowMeans(mat, na.rm = TRUE), r_min = r_min)
            cbind(data.frame(curve = d$group[1], level = "group",
                             stringsAsFactors = FALSE), f)
          }))
          if (!is.null(per_grp) && nrow(per_grp))
            out <- rbind(per_grp, per_file[per_file$curve %in% names(gm), , drop = FALSE])
        }
      }
    }
    out
  })

  output$pcf_xi_table <- renderTable({
    x <- pcf_xi_results()
    shiny::validate(shiny::need(!is.null(x) && nrow(x), "Run g(r) first."))
    data.frame(
      Curve   = x$curve,
      Level   = x$level,
      `xi (nm)` = ifelse(is.na(x$xi), "-", sprintf("%.0f", x$xi)),
      `± SE`  = ifelse(is.na(x$xi_se), "-", sprintf("%.0f", x$xi_se)),
      A       = ifelse(is.na(x$amplitude), "-", sprintf("%.2f", x$amplitude)),
      `R²`    = ifelse(is.na(x$r2), "-", sprintf("%.3f", x$r2)),
      `fit range (nm)` = sprintf("%.0f-%.0f", x$r_min, x$r_max),
      n       = x$n_used,
      check.names = FALSE, stringsAsFactors = FALSE)
  }, striped = TRUE, spacing = "xs", width = "100%")

  output$pcf_xi_note <- renderUI({
    x <- pcf_xi_results(); if (is.null(x) || !nrow(x)) return(NULL)
    msgs <- list()
    if (any(is.na(x$xi)))
      msgs <- c(msgs, sprintf("%d curve(s) show no decaying excess over the fit range - no correlation length is reported for them, rather than a number fitted to noise.",
                              sum(is.na(x$xi))))
    if (any(isTRUE(x$below_floor)) || any(x$below_floor %in% TRUE))
      msgs <- c(msgs, "At least one ξ is SMALLER than the fit's starting radius, so the fit began past the decay it was measuring - lower 'Fit from r' for those curves before quoting the value.")
    if (any(grepl("log-linear", x$method %||% "")))
      msgs <- c(msgs, "One or more fits fell back to the log-linear estimator (nls did not converge); that estimator over-estimates ξ on noisy tails.")
    if (!length(msgs)) return(NULL)
    tags$div(style = "font-size:0.82em; color:#7a5c00; margin-top:4px;",
             lapply(msgs, function(m) tags$div(m)))
  })

  build_pcf_plot_fig <- function() {
    res <- pcf_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run g(r) first."))

    do_smooth <- isTRUE(input$pcf_smooth)
    span_val  <- as.numeric(input$pcf_smooth_span %||% 0.25)
    smooth_vec <- function(r, y) {
      ok <- is.finite(y)
      if (!do_smooth || sum(ok) < 5L) return(y)
      tryCatch(stats::predict(stats::loess(y ~ r, data = data.frame(r = r, y = y),
                                           span = span_val, na.action = stats::na.exclude),
                              newdata = data.frame(r = r)),
               error = function(e) y)
    }
    hex2rgba <- function(hex, a) {
      cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
    }
    grouped <- (input$analysis_by %||% "file") == "group"
    env_res <- pcf_env_results()
    # Declared up front because the axis code below reads them in both branches.
    rmin_c <- NA_real_; rmax_c <- NA_real_
    show_files <- isTRUE(input$pcf_show_files)
    fig <- plotly::plot_ly()

    if (!grouped) {
      datasets <- unique(res$dataset)
      pal <- series_palette(datasets)
      for (nm in datasets) {
        sub <- res[res$dataset == nm, ]; sub <- sub[order(sub$r), ]
        if (!is.null(env_res) && nm %in% env_res$dataset) {
          # The BAND MUST GET THE SAME SMOOTHING AS THE CURVE. Comparing a smoothed
          # observed curve against a raw envelope compares two different statistics:
          # smoothing damps the curve's excursions but left the band at its full
          # unsmoothed width, so the comparison silently became conservative (and the
          # band looked jagged next to a smooth curve). Smoothing the quantiles is a
          # close approximation to taking quantiles of smoothed simulations - the two
          # differ only in how a monotone smoother interacts with the quantile, which
          # is negligible here and keeps the band consistent with whatever smoothing
          # is displayed (smoothing is a live display control, the envelope is not).
          e <- env_res[env_res$dataset == nm, ]; e <- e[order(e$r), ]
          e_hi <- smooth_vec(e$r, e$hi); e_lo <- smooth_vec(e$r, e$lo)
          fig <- fig %>%
            plotly::add_trace(x = e$r, y = e_hi, type = "scatter", mode = "lines",
              line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
            plotly::add_trace(x = e$r, y = e_lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
              legendgroup = nm, name = paste0(nm, " CSR 95% envelope"), hoverinfo = "skip")
          fig <- .add_global_band(fig, e, unname(pal[nm]), nm, smooth_vec)
        }
        fig <- fig %>% plotly::add_trace(
          x = sub$r, y = smooth_vec(sub$r, sub$g), type = "scatter", mode = "lines",
          name = sprintf("%s (n=%d)", nm, sub$n_units[1]), legendgroup = nm,
          line = list(color = unname(pal[nm]), width = 2), hoverinfo = "skip")
      }
    } else {
      gm <- file_group_map()
      shiny::validate(shiny::need(length(gm) > 0, "No files assigned to groups (group mode)."))
      res <- res[res$dataset %in% names(gm), , drop = FALSE]
      shiny::validate(shiny::need(nrow(res) > 0, "No group files in the g(r) results."))
      res$group <- unname(gm[res$dataset])
      pal  <- series_palette(unique(res$group))
      fpal <- file_palette(unique(res$dataset))
      band <- input$pcf_group_band %||% "iqr"

      # Same contract as the Ripley group mode: interpolate each FILE's curve onto a
      # common r-grid and average within the group. Coordinates are never pooled -
      # pooling would mix windows and densities and invent structure.
      file_maxr <- tapply(res$r, res$dataset, function(v) max(v, na.rm = TRUE))
      rmax_c <- suppressWarnings(min(file_maxr, na.rm = TRUE))
      # Common grid spans only the range EVERY file actually measured. It used to
      # start at r = 0 while the stored curves start at r_max/steps, so approx(rule=2)
      # flat-extrapolated the first real value back to r = 0 - a false non-zero point
      # at a radius where L(r)-r is 0 by definition, injected into every group curve
      # and its band.
      file_minr <- tapply(res$r, res$dataset, function(v) min(v, na.rm = TRUE))
      rmin_c <- suppressWarnings(max(file_minr, na.rm = TRUE))
      if (!is.finite(rmin_c) || rmin_c < 0) rmin_c <- 0
      shiny::validate(shiny::need(is.finite(rmax_c) && rmax_c > 0, "Could not form a common r-grid."))
      # See the Ripley note: grid density comes from the stored curve.
      rr <- seq(rmin_c, rmax_c, length.out = max(10L, length(unique(res$r))))

      for (g in unique(res$group)) {
        files <- unique(res$dataset[res$group == g])
        mat <- vapply(files, function(f) {
          d <- res[res$dataset == f, ]; d <- d[order(d$r), ]
          stats::approx(d$r, d$g, xout = rr, rule = 2)$y
        }, numeric(length(rr)))
        if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
        meanc <- smooth_vec(rr, rowMeans(mat, na.rm = TRUE))
        col <- unname(pal[g]); nfl <- length(files)
        # ── BAND MUST MATCH THE DRAWN MEAN ──────────────────────────────────
        # Three faults were fixed here together:
        #  1. the mean was smoothed and the band was not, so a smooth line sat in a
        #     jagged band (measured: band edge 16-18x the roughness of the mean).
        #     Same argument as the CSR envelope above - a band and a curve that got
        #     different treatment are not comparable.
        #  2. the SD band was built around the RAW rowMeans while the SMOOTHED mean
        #     was drawn, so the line was not centred in its own band (up to 22 nm
        #     off). It is now centred on the curve actually shown.
        #  3. quartiles need >= 4 values. On 2 files the IQR band is the middle half
        #     of a 2-point range and contains NEITHER file (measured 0% coverage),
        #     which draws a reassuringly tight band over maximal disagreement. Under
        #     4 files the full min-max range is drawn instead, and the legend says so.
        band_lbl <- ""
        if (band != "none" && nfl >= 2) {
          if (band == "iqr" && nfl >= 4L) {
            lo <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.25, na.rm = TRUE))
            hi <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.75, na.rm = TRUE))
            band_lbl <- ", IQR"
          } else if (band == "iqr") {
            lo <- smooth_vec(rr, apply(mat, 1, min, na.rm = TRUE))
            hi <- smooth_vec(rr, apply(mat, 1, max, na.rm = TRUE))
            band_lbl <- ", range"
          } else {
            sdv <- smooth_vec(rr, apply(mat, 1, stats::sd, na.rm = TRUE))
            lo <- meanc - sdv; hi <- meanc + sdv
            band_lbl <- ", ±1 SD"
          }
          fig <- fig %>%
            plotly::add_trace(x = rr, y = hi, type = "scatter", mode = "lines",
              line = list(width = 0), legendgroup = g, showlegend = FALSE, hoverinfo = "skip") %>%
            plotly::add_trace(x = rr, y = lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = hex2rgba(col, 0.18), line = list(width = 0),
              legendgroup = g, showlegend = FALSE, hoverinfo = "skip")
        }
        # Individual files as shades of the group colour, behind the group mean.
        if (show_files && nfl >= 2) {
          for (f in files) {
            d <- res[res$dataset == f, ]; d <- d[order(d$r), ]
            fc <- unname(fpal[f]); if (is.na(fc)) fc <- col
            fig <- fig %>% plotly::add_trace(
              x = d$r, y = d$g, type = "scatter", mode = "lines",
              name = f, legendgroup = g, showlegend = FALSE,
              line = list(color = fc, width = 1), opacity = 0.75, hoverinfo = "skip")
          }
        }
        # CSR envelope for the group - same reasoning as the Ripley panel: the
        # group curve is the mean of its files' curves, so its null is the mean of
        # those files' envelopes. Dotted, so it reads as a reference rather than
        # as the solid between-file spread band.
        env_n <- 0L
        if (!is.null(env_res)) {
          ef <- intersect(files, unique(env_res$dataset))
          if (length(ef)) {
            elo <- vapply(ef, function(f) {
              e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
              stats::approx(e$r, e$lo, xout = rr, rule = 2)$y }, numeric(length(rr)))
            ehi <- vapply(ef, function(f) {
              e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
              stats::approx(e$r, e$hi, xout = rr, rule = 2)$y }, numeric(length(rr)))
            if (is.null(dim(elo))) { elo <- matrix(elo, ncol = 1); ehi <- matrix(ehi, ncol = 1) }
            env_n <- length(ef)
            # Two lines per group: the upper and lower bound of one band. Own
            # legend entry and legendgroup so they are identifiable and can be
            # toggled off without hiding the group curve.
            fig <- fig %>% plotly::add_trace(
              x = rr, y = smooth_vec(rr, rowMeans(ehi, na.rm = TRUE)),
              type = "scatter", mode = "lines",
              line = list(color = col, width = 1, dash = "dot"),
              name = sprintf("%s – CSR 95%% band", g),
              legendgroup = paste0(g, "__csr"), showlegend = TRUE, hoverinfo = "skip")
            fig <- fig %>% plotly::add_trace(
              x = rr, y = smooth_vec(rr, rowMeans(elo, na.rm = TRUE)),
              type = "scatter", mode = "lines",
              line = list(color = col, width = 1, dash = "dot"),
              legendgroup = paste0(g, "__csr"), showlegend = FALSE, hoverinfo = "skip")
          }
        }
        env_lbl <- if (env_n > 0L && env_n < nfl) sprintf(", CSR n=%d", env_n)
                   else if (env_n > 0L) ", CSR" else ""
        fig <- fig %>% plotly::add_trace(x = rr, y = meanc, type = "scatter", mode = "lines",
          name = sprintf("%s (n=%d%s%s)", g, nfl, band_lbl, env_lbl), legendgroup = g,
          line = list(color = col, width = 2.5), hoverinfo = "skip")
      }
    }

    # Same as the Ripley plot: the reference line spans what is actually drawn, so
    # a file with a larger r max cannot stretch the axis past every curve.
    x_drawn <- if (grouped) {
      c(rmin_c, if (isTRUE(show_files)) max(res$r, na.rm = TRUE) else rmax_c)
    } else range(res$r, na.rm = TRUE)

    fig %>%
      plotly::add_trace(x = x_drawn, y = c(1, 1), type = "scatter",
        mode = "lines", name = "CSR (g = 1)",
        line = list(color = "black", width = 1, dash = "dash"), hoverinfo = "skip") %>%
      plotly::layout(
        xaxis = list(title = .rmax_axis_title(grouped, res, rmax_c, show_files)),
        yaxis = list(title = "g(r)"),
        annotations = .global_p_annot(env_res),
        # Hover off, as on the Ripley plot: the curve is read as a shape against the
        # envelope, and a tooltip per series at every r covered the plot.
        hovermode = FALSE,
        legend = list(groupclick = "togglegroup"),
        paper_bgcolor = "white", plot_bgcolor = "white")
  }
  output$pcf_plot <- renderPlotly(build_pcf_plot_fig())

  ## Summary numbers - the actual endpoint you would compare between conditions.
  ## Reports the MEAN g over a user-set window rather than only the peak: a max
  ## over ~50 bins is upward-biased at low unit counts, so peak g is not safely
  ## comparable between files of different density even though g itself is.
  ## Per-FILE scalars first, then the group summary - never scalars from pooled rows.
  ##
  ## WHY THIS ORDER MATTERS. A group's rows are several files' curves concatenated, so
  ## reducing them together makes the CELL stop being the unit of replication: the peak
  ## becomes a max over files, and a denser file contributes more rows and so more
  ## weight. Reducing each file first and then averaging the scalars keeps one vote per
  ## cell, which is the standard unit of replication in SMLM, and it yields a spread
  ## across cells - the only thing that says whether a between-condition difference is
  ## bigger than the between-cell variability.
  ##
  ## Measured on the four real a real imaging file acquisitions (mean g over 30-300 nm, common r grid):
  ## cond01 1.24 / 1.50, cond02 2.31 / 4.38. Between-condition 2.44x against a
  ## within-condition spread of 1.21x and 1.89x, and no overlap - so averaging is
  ## justified there. But at r = 20 nm the two cond02 files read 5.03 and 14.32, so the
  ## MEAN CURVE (9.68) describes neither: read strengths and sizes off the per-file
  ## rows, and use the group row only for the comparison.
  output$pcf_summary <- renderUI({
    res <- pcf_results()
    if (is.null(res) || !nrow(res)) return(NULL)
    r1 <- as.numeric(input$pcf_win_lo %||% 0)
    r2 <- as.numeric(input$pcf_win_hi %||% 0)
    if (!is.finite(r1) || !is.finite(r2) || r2 <= r1) {
      rmx <- max(res$r, na.rm = TRUE); r1 <- 0.05 * rmx; r2 <- 0.4 * rmx
    }
    grouped <- (input$analysis_by %||% "file") == "group"

    # One row per FILE, always computed on the file's own curve.
    per_file <- do.call(rbind, lapply(split(res, res$dataset), function(d) {
      d  <- d[order(d$r), ]
      ok <- is.finite(d$g) & d$r >= r1 & d$r <= r2
      pk <- if (any(is.finite(d$g))) which.max(replace(d$g, !is.finite(d$g), -Inf)) else NA_integer_
      # Size readout: first radius AFTER the peak where g returns to ~1. NA when the
      # curve has not settled inside the r range, which is itself informative -
      # the structure is larger than the window can resolve.
      rz <- NA_real_
      if (!is.na(pk)) {
        ti <- which(seq_along(d$r) > pk & is.finite(d$g) & d$g <= 1.05)
        if (length(ti)) rz <- d$r[ti[1]]
      }
      data.frame(file = d$dataset[1],
                 g_win = if (any(ok)) mean(d$g[ok]) else NA_real_,
                 excess = if (any(ok)) sum((d$g[ok] - 1) * c(diff(d$r[ok]), 0)) else NA_real_,
                 r_settle = rz, rmax = max(d$r, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
    gm <- if (grouped) file_group_map() else character(0)
    per_file$group <- if (grouped) unname(gm[per_file$file]) else per_file$file
    per_file <- per_file[!is.na(per_file$group), , drop = FALSE]
    if (!nrow(per_file)) return(NULL)

    num <- function(v, dg = 2) if (!is.finite(v)) "-" else formatC(v, format = "f", digits = dg)
    hdr <- function(cols) tags$tr(lapply(cols, function(h)
      tags$th(style = "text-align:left; padding:2px 12px 2px 0; border-bottom:1px solid #ddd;", h)))
    row_td <- function(...) tags$tr(lapply(list(...), function(v)
      tags$td(style = "padding:1px 12px 1px 0;", v)))

    file_tbl <- tags$table(style = "border-collapse:collapse;",
      hdr(c("File", if (grouped) "Group" else NULL,
            sprintf("mean g (%.0f-%.0f nm)", r1, r2), "∫(g-1)dr (nm)",
            "r where g→1 (nm)", "r max")),
      lapply(seq_len(nrow(per_file)), function(i) {
        args <- list(per_file$file[i])
        if (grouped) args <- c(args, list(per_file$group[i]))
        args <- c(args, list(num(per_file$g_win[i]),
                             num(per_file$excess[i], 0),
                             if (is.finite(per_file$r_settle[i])) num(per_file$r_settle[i], 0)
                             else sprintf("> %.0f", per_file$rmax[i]),
                             num(per_file$rmax[i], 0)))
        do.call(row_td, args)
      }))

    out <- list(tags$div(style = "font-weight:600; margin:2px 0 3px;",
                         "Per file (one row per cell - this is the unit of replication)"),
                file_tbl)

    if (grouped) {
      ag <- do.call(rbind, lapply(split(per_file, per_file$group), function(d) data.frame(
        group = d$group[1], n = nrow(d),
        mean_g = mean(d$g_win, na.rm = TRUE),
        sd_g   = if (nrow(d) > 1) stats::sd(d$g_win, na.rm = TRUE) else NA_real_,
        lo = min(d$g_win, na.rm = TRUE), hi = max(d$g_win, na.rm = TRUE),
        stringsAsFactors = FALSE)))
      grp_tbl <- tags$table(style = "border-collapse:collapse;",
        hdr(c("Group", "n cells", sprintf("mean g (%.0f-%.0f nm)", r1, r2),
              "SD across cells", "range")),
        lapply(seq_len(nrow(ag)), function(i) row_td(
          ag$group[i], ag$n[i], num(ag$mean_g[i]),
          if (is.finite(ag$sd_g[i])) num(ag$sd_g[i]) else "-",
          sprintf("%s - %s", num(ag$lo[i]), num(ag$hi[i])))))
      # The decisive comparison: is the between-group gap bigger than the within-group
      # spread? With no overlap it is, and that does not need a test to state.
      verdict <- NULL
      if (nrow(ag) == 2 && all(is.finite(c(ag$lo, ag$hi)))) {
        sep <- ag$lo[which.max(ag$mean_g)] > ag$hi[which.min(ag$mean_g)]
        ratio <- max(ag$mean_g) / min(ag$mean_g)
        verdict <- tags$p(style = paste0("margin:5px 0 0;",
                            if (sep) "color:#1a8a1a;" else "color:#b26a00;"),
          sprintf("Group means differ %.2f×; per-cell ranges %s.", ratio,
                  if (sep) "do NOT overlap, so the difference exceeds between-cell variability"
                  else "OVERLAP, so this difference is within between-cell variability"),
          if (min(ag$n) < 3) " With fewer than 3 cells per group this is descriptive only - no test is meaningful." else "")
      }
      out <- c(out, list(tags$div(style = "font-weight:600; margin:10px 0 3px;",
                                  "Per group (mean of the per-cell values)"),
                        grp_tbl, verdict))
    }

    out <- c(out, list(tags$p(style = "margin:8px 0 0; color:#7a5c00;",
      "g = 1 is complete randomness. Read STRENGTH from mean g and SIZE from the radius where g returns to 1. ",
      "A value shown as '> r max' means the structure is larger than this window can resolve - set a fixed r max so files are comparable. ",
      "Do not read a size off a group mean curve: if the cells differ, the average describes none of them.")))
    do.call(tags$div, c(list(style = "font-size:0.85em;"), out))
  })

  output$download_pcf <- downloadHandler(
    filename = function() paste0("pcf_g_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- pcf_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  ## ══════════════════════════════════════════════════════════════════════════
  ## LOCAL DENSITY & HOTSPOTS
  ## ══════════════════════════════════════════════════════════════════════════
  ## Computes, per file: the cell footprint, every unit's local density as a
  ## MULTIPLE OF THAT CELL'S OWN MEAN, a CSR-calibrated DBSCAN eps, and the foci
  ## that eps finds. See the block banner "8c. LOCAL DENSITY & HOTSPOTS" in
  ## helpers_smlm.R for why these statistics and not nearest-neighbour distance.
  ##
  ## THE FOOTPRINT IS NOT OPTIONAL HERE, which is the one place this panel departs
  ## from the rest of the tab. Elsewhere `spatial_use_support` is an opt-in that
  ## improves an estimate; here every statistic is "per unit of this cell's own
  ## area", so a bounding box does not lose power, it silently changes the number:
  ## a cell imaged half in frame reads at half its true density and every
  ## enhancement doubles. The panel therefore resolves its own footprint with
  ## `force = TRUE` and says so in the readout.

  hotspot_rows   <- reactiveVal(NULL)   # one summary row per file
  hotspot_detail <- reactiveVal(NULL)   # per-file pts / support / enh / foci, for the map
  hotspot_params <- reactiveVal(NULL)   # resolved params, frozen at Run
  hotspot_unit   <- reactiveVal(NULL)   # the unit the stored results were built on

  ## GOAL: say on the axis why a group curve stops where it does.
  ##
  ## In group mode the curves are averaged on a grid capped at the SMALLEST
  ## per-file r max, so if one cell was measured to 400 nm and another to 1500 the
  ## group curve ends at 400. That is deliberate - averaging past it would compare
  ## a mean of five files against a mean of one - but with the axis now ending
  ## there too, nothing on screen says the longer files were truncated. This puts
  ## it in the axis title, where it costs no space and cannot be missed.
  .rmax_axis_title <- function(grouped, res, rmax_c, show_files) {
    base <- "r (nm)"
    if (!isTRUE(grouped) || isTRUE(show_files)) return(base)
    rmax_all <- suppressWarnings(max(res$r, na.rm = TRUE))
    if (!is.finite(rmax_c) || !is.finite(rmax_all) || rmax_all <= rmax_c * 1.02)
      return(base)
    sprintf("r (nm) - group curves stop at %.0f nm, the smallest r max in the set (longest file reaches %.0f)",
            rmax_c, rmax_all)
  }

  ## GOAL: pick the point set that draws the observation window, per the shared
  ## `dens_support_from` setting. `locs` is the localization table, `units` the
  ## one-point-per-trace/molecule centroids. See spatial_window_of() for the
  ## measured trade-off; the analyses themselves always run on the units.
  spatial_outline_pts <- function(locs_df, unit_pts) {
    if (identical(input$dens_support_from %||% "units", "locs") &&
        is.data.frame(locs_df) && nrow(locs_df)) return(locs_df)
    if (is.data.frame(unit_pts) && nrow(unit_pts)) unit_pts else locs_df
  }

  ## Look a name up in a named vector WITHOUT the `[[` trap. `v[["absent"]]` raises
  ## "subscript out of bounds" rather than returning NULL, so `%||%` never sees it
  ## and cannot supply the default - which is how an ungrouped file took the whole
  ## enhancement curve down. Single-bracket indexing yields NA for a missing name,
  ## which this turns into the default.
  .name_or <- function(v, key, default) {
    if (is.null(v) || !length(v) || is.null(key) || !length(key)) return(default)
    out <- unname(v[as.character(key)[1]])
    if (length(out) != 1L || is.na(out)) default else out
  }

  ## THE UNIT IS THE WHOLE MEASUREMENT, so it is named in words on every axis and
  ## title rather than left to the control panel.
  ##
  ## RAW LOCALIZATIONS ARE NOT OFFERED, and this is not a convenience choice.
  ## Localizations inside one trace sit a few nm apart, so the k-th neighbour
  ## distance collapses to the localization precision and the local density
  ## measures HOW LONG THAT FLUOROPHORE STAYED ON, not how the protein is
  ## arranged. Measured on the two-condition pair: on localizations the enhancement runs
  ## to 3,000-17,000x the cell mean and NEITHER the CV nor the 90th percentile
  ## separates the conditions (CV 1.46-1.74 against 1.23-2.54, overlapping), while
  ## the same statistics on trace-ID centroids are cleanly disjoint. Trace length
  ## is photophysics; the spatial arrangement of separate traces is the biology.
  hotspot_unit_label <- function(u) switch(u %||% "tid",
    tid     = "trace ID centroids",
    cluster = "de-blinked molecules",
    "units")

  observeEvent(input$run_hotspot, {
    unit    <- resolve_unit(input$hot_unit, "cluster", "tid", "Local density & hotspots")
    ulab    <- hotspot_unit_label(unit)
    use_3d  <- isTRUE(input$hot_use_3d)
    k       <- max(1L, as.integer(input$hot_k %||% 5L))
    quad    <- max(10, as.numeric(input$hot_quadrat %||% 250))
    thresh  <- max(1, as.numeric(input$hot_enh_thresh %||% 5))
    minpts  <- max(2L, as.integer(input$hot_minpts %||% 5L))
    target  <- min(0.5, max(0.001, as.numeric(input$hot_target_frac %||% 0.05)))
    nsim    <- max(3L, as.integer(input$hot_nsim %||% 15L))
    # The window is NOT this panel's to define. Pitch, closing radius and estimator
    # mode all come from "Density & area normalization", so this panel measures the
    # same cell area that Ripley's L and g(r) measure. See spatial_window_of().
    h_nm    <- max(10, as.numeric(input$dens_grid_h %||% 100))
    needs_cluster <- unit == "cluster"

    store <- smlm_file_store()
    shiny::validate(shiny::need(length(store) > 0, "No files loaded."))

    # ---- 1. prepare every file, then take one point per unit -------------------
    items <- lapply(names(store), function(nm) {
      df <- prepared_df(nm, needs_cluster = needs_cluster)
      if (is.null(df) || !nrow(df)) return(NULL)
      pts <- unit_centroids(df, unit = unit, use_3d = use_3d)
      if (is.null(pts) || nrow(pts) < max(3L, minpts)) return(NULL)
      list(tag = nm, df = df, pts = pts)
    })
    items <- Filter(Negate(is.null), items)
    shiny::validate(shiny::need(length(items) > 0,
      "No file had enough units. Check the filters, or lower minPts."))

    # ---- 2. ONE gap-fill radius for the whole set -----------------------------
    # Shared, for the same reason the curve panels share theirs: the footprint sets
    # the area, the area sets the density, and a per-file radius would put each
    # cell on its own density scale - which is precisely the comparison this panel
    # exists to make. Resolved on the UNIT positions (not the raw localizations),
    # because the units are what the density is measured on.
    # THE FOOTPRINT MUST NOT BE BUILT FROM THE UNITS, so the sweep runs on every
    # LOCALIZATION, exactly as the curve panels and the density panel do. The
    # support is a property of where signal was detected; estimating it from unit
    # centroids makes the area depend on the unit count, which is the very thing
    # being normalised out. Measured (see compute_unit_density): two synthetic
    # cells with identical geometry and a 4x labelling difference gave
    # centroid-based supports of 5.01e6 against 3.82e6 nm^2, turning a true 4.0x
    # density ratio into 3.05x.
    cc <- spatial_gate_close(lapply(items, function(it)
            list(df = spatial_outline_pts(it$df, it$pts), tag = it$tag)),
            use_3d, force = TRUE, h_nm = h_nm)
    close_nm <- if (!is.null(cc$close_nm) && is.finite(cc$close_nm)) cc$close_nm else
      max(0, as.numeric(input$dens_close_nm %||% 200))

    # ---- 3. per file: footprint -> enhancement -> calibrated eps -> foci -------
    out <- withProgress(message = sprintf("Local density: %d file(s)...", length(items)),
                        value = 0, {
      lapply(seq_along(items), function(i) {
        it <- items[[i]]
        incProgress(1 / length(items), detail = it$tag)
        # The SAME call the curve panels make, on the SAME point set (all
        # localizations - see the note on the radius above), so the window is
        # logged once and the "these files were not measured the same way" check
        # covers this panel too. The unit centroids are what the STATISTICS run
        # on; they are not what defines the cell.
        sup <- tryCatch(spatial_support_of(spatial_outline_pts(it$df, it$pts),
                                           use_3d, tag = it$tag,
                                           close_nm = close_nm, force = TRUE),
                        error = function(e) NULL)
        le <- compute_local_enhancement(it$pts, k = k, use_3d = use_3d, support = sup)
        if (is.null(le)) return(NULL)
        ec <- tryCatch(calibrate_dbscan_eps(it$pts, support = sup, minPts = minpts,
                                            target_frac = target, nsim = nsim,
                                            use_3d = use_3d),
                       error = function(e) NULL)
        ft <- if (!is.null(ec))
          compute_foci_table(it$pts, eps = ec$eps, minPts = minpts, use_3d = use_3d) else NULL
        row <- hotspot_summary_row(it$pts, support = sup, k = k, quadrat_nm = quad,
                                   enh_thresh = thresh, eps_cal = ec, minPts = minpts,
                                   use_3d = use_3d, dataset = it$tag)
        # bbox_support() returns an all-TRUE mask and a finite area, so a bounding
        # box does not fail here - it silently measures the imaged FRAME instead of
        # the cell. Record the mode so the run can say so.
        row$window_mode <- sup$mode %||% NA_character_
        list(tag = it$tag, pts = it$pts, support = sup, enh = le$enh, le = le,
             eps_cal = ec, foci = ft, row = row)
      })
    })
    out <- Filter(Negate(is.null), out)
    shiny::validate(shiny::need(length(out) > 0, "Local density returned no data."))
    names(out) <- vapply(out, function(o) o$tag, character(1))

    rows <- do.call(rbind, lapply(out, `[[`, "row"))
    rows$close_nm <- close_nm
    hotspot_rows(rows)
    hotspot_unit(unit)
    hotspot_detail(out)
    hotspot_params(list(unit = input$hot_unit, k = input$hot_k,
                        `quadrat (nm)` = input$hot_quadrat,
                        `enhancement threshold` = input$hot_enh_thresh,
                        minPts = input$hot_minpts,
                        `false-cluster target` = input$hot_target_frac,
                        `CSR draws` = input$hot_nsim,
                        `3-D` = use_3d))

    # The map file picker only ever offers files that actually produced a result.
    spatial_report_windows(vapply(items, function(it) as.character(it$tag),
                                  character(1)))

    updateSelectInput(session, "hot_map_file", choices = names(out),
                      selected = if (!is.null(input$hot_map_file) &&
                                     input$hot_map_file %in% names(out))
                        input$hot_map_file else names(out)[1])

    # A bounding-box window makes every number here a per-FRAME density rather than
    # a per-CELL one, which is a different quantity, not a noisier version of the
    # same one. It happens two ways: the estimator mode is set to bbox or hull (no
    # mask), or the footprint guard rejected the window for that file.
    fb <- rows$dataset[grepl("bbox", rows$window_mode %||% "", fixed = TRUE)]
    if (length(fb) == nrow(rows))
      showNotification(HTML(sprintf(
        "<b>Every window is a bounding box</b>, so densities are per imaged frame, not per cell. Set the estimator to <b>grid</b> or <b>kde</b> under Density &amp; area normalization. (%s)",
        rows$window_mode[1])), type = "warning", duration = NULL)
    else if (length(fb))
      showNotification(HTML(sprintf(
        "<b>No footprint for %s.</b> Those files fell back to the bounding box and are not on the same scale as the rest. Raise the gap-fill radius or check the support map.",
        paste(utils::head(fb, 4), collapse = ", "))), type = "warning", duration = NULL)
  })

  output$hotspot_stale <- renderUI(.stale_badge(
    list(unit = input$hot_unit, k = input$hot_k, `quadrat (nm)` = input$hot_quadrat,
         `enhancement threshold` = input$hot_enh_thresh, minPts = input$hot_minpts,
         `false-cluster target` = input$hot_target_frac, `CSR draws` = input$hot_nsim,
         `3-D` = isTRUE(input$hot_use_3d)),
    hotspot_params()))
  # hot_curve_mode and hot_metric are deliberately NOT in the stale badge: they
  # re-draw stored results and never change what was computed.

  # ── Readout: what the run resolved, in one line ────────────────────────────
  output$hotspot_status <- renderUI({
    rows <- hotspot_rows()
    if (is.null(rows) || !nrow(rows))
      return(helpText("Press Run to measure local density and detect hotspots."))
    det <- hotspot_detail()
    floors <- sum(vapply(det, function(o) isTRUE(o$eps_cal$hit_floor), logical(1)))
    tags$div(style = paste0("margin:4px 0 8px; padding:6px 11px; border-left:4px solid #cfe0ee;",
                            "background:#f6f9fc; font-size:0.86em; color:#445;"),
      sprintf("%d file(s), measured on %s. Window: %s, gap fill %.0f nm, pitch %.0f nm (from Density & area normalization). Footprint %.1f-%.1f um2, density %.1f-%.1f per um2.",
              nrow(rows), hotspot_unit_label(hotspot_unit()),
              input$dens_support_mode %||% "grid",
              rows$close_nm[1], max(10, as.numeric(input$dens_grid_h %||% 100)),
              min(rows$area_um2), max(rows$area_um2),
              min(rows$density), max(rows$density)),
      # eps is reported separately because a file whose calibration failed leaves
      # NA here, and min()/max() over an all-NA column returns Inf with a warning.
      if (any(is.finite(rows$eps_nm)))
        sprintf(" Calibrated eps %.0f-%.0f nm.",
                min(rows$eps_nm, na.rm = TRUE), max(rows$eps_nm, na.rm = TRUE))
      else " eps was not calibrated on any file.",
      if (floors > 0) tags$div(style = "color:#a04000; margin-top:4px;",
        sprintf("%d file(s) hit the bottom of the eps grid - a random pattern still clustered at the smallest eps tried, so their focus numbers are not calibrated.", floors)))
  })

  # ── The map: footprint + units coloured by local-density enhancement ───────
  build_hotspot_map_fig <- function() {
    det <- hotspot_detail()
    shiny::validate(shiny::need(!is.null(det) && length(det) > 0, "Run this panel first."))
    fname <- input$hot_map_file %||% names(det)[1]
    shiny::validate(shiny::need(fname %in% names(det), "Pick a dataset."))
    o <- det[[fname]]

    dark   <- isTRUE(input$dark_mode)
    bg     <- if (dark) "#1e1e1e" else "white"
    fgcol  <- if (dark) "#e6e6e6" else "#333a45"
    maskc  <- if (dark) "rgba(255,255,255,0.09)" else "rgba(120,132,148,0.13)"

    pts <- o$pts[o$le$keep, , drop = FALSE]
    enh <- o$enh
    # Log colour scale: enhancement spans two orders of magnitude in a clustered
    # cell, so a linear scale renders every background unit the same colour and
    # only the single densest spot is visible.
    lo  <- log10(max(min(enh, na.rm = TRUE), 0.1))
    hi  <- log10(max(enh, na.rm = TRUE))
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) { lo <- -1; hi <- 1 }

    fig <- plotly::plot_ly()
    map_shapes <- list(); map_annots <- list()
    # Footprint underneath, as a flat mask - this is the "crop to the cell" the
    # density numbers are actually divided by, so it must be visible, not implied.
    if (!is.null(o$support) && !is.null(o$support$mask)) {
      m2   <- mask_project_2d(o$support$mask)
      cell <- o$support$h
      d    <- dim(m2)
      xc <- o$support$origin[1] + (seq_len(d[1]) - 0.5) * cell
      yc <- o$support$origin[2] + (seq_len(d[2]) - 0.5) * cell
      zin <- ifelse(m2, 1, NA_real_)
      fig <- fig %>% plotly::add_heatmap(
        x = xc, y = yc, z = t(zin), showscale = FALSE, hoverinfo = "skip",
        colorscale = list(list(0, maskc), list(1, maskc)))
    }

    fig <- fig %>% plotly::add_markers(
      x = pts$x, y = pts$y, type = "scatter", mode = "markers",
      marker = list(size = 5, color = log10(pmax(enh, 0.1)),
                    colorscale = plotly_colorscale(input$hot_palette %||% "Turbo"),
                    cmin = lo, cmax = hi, line = list(width = 0),
                    colorbar = list(title = list(text = "local density\n/ cell mean"),
                                    tickvals = log10(c(0.1, 0.3, 1, 3, 10, 30, 100)),
                                    ticktext = c("0.1x", "0.3x", "1x", "3x", "10x", "30x", "100x"),
                                    len = 0.8)),
      text = sprintf("%.1fx cell mean", enh), hoverinfo = "text",
      name = hotspot_unit_label(hotspot_unit()))

    # ── Box the detected foci, at their real size ─────────────────────────
    # A marker was the wrong mark for this. `marker$size` is in PIXELS, so a ring
    # sat at a fixed screen size whatever the focus measured and whatever the zoom:
    # on foci spanning 100-490 nm it read as a dot in the middle of the group
    # rather than as its outline. These are rectangle SHAPES in data coordinates,
    # so each box is the size of the focus it encloses and stays correct on zoom.
    #
    # Every focus is boxed, not only the biggest, and each carries its count. That
    # is deliberate: "largest" ranks by member count while the eye reads DENSITY,
    # and the two disagree in half the two-condition files - a big loose focus outranks a
    # small tight one. Boxing all of them with their counts removes the question of
    # which one the panel picked.
    if (isTRUE(input$hot_show_foci) && !is.null(o$foci) && nrow(o$foci$table)) {
      tb  <- o$foci$table
      # Boxing every focus is unreadable on a busy cell - one two-condition file yields 23.
      # The table is ordered largest first, so taking the head keeps the ones worth
      # looking at. 0 means "all", for when the spatial arrangement of the small
      # ones is the point.
      ntop <- as.integer(input$hot_box_top %||% 5L)
      if (is.finite(ntop) && ntop > 0L) tb <- utils::head(tb, ntop)
      fp  <- o$pts[o$foci$keep, , drop = FALSE]
      cl  <- o$foci$cluster
      ring <- if (dark) "#ffd166" else "#c0392b"
      shp <- list(); ann <- list()
      for (ii in seq_len(nrow(tb))) {
        m <- as.matrix(fp[cl == tb$focus[ii], c("x", "y"), drop = FALSE])
        if (!nrow(m)) next
        rx <- range(m[, 1]); ry <- range(m[, 2])
        # A little air, floored so a very tight focus still gets a visible box.
        pad <- max(0.10 * max(diff(rx), diff(ry)), 20)
        # The biggest focus gets the heavier line, so the table's headline number
        # is findable on the map without counting.
        lw <- if (ii == 1L) 2.2 else 1.2
        shp[[length(shp) + 1L]] <- list(
          type = "rect", xref = "x", yref = "y",
          x0 = rx[1] - pad, x1 = rx[2] + pad,
          y0 = ry[1] - pad, y1 = ry[2] + pad,
          line = list(color = ring, width = lw), fillcolor = "rgba(0,0,0,0)",
          layer = "above")
        ann[[length(ann) + 1L]] <- list(
          x = mean(rx), y = ry[2] + pad, xref = "x", yref = "y",
          text = sprintf("%d", tb$n[ii]), showarrow = FALSE,
          yanchor = "bottom", font = list(size = 11, color = ring))
      }
      map_shapes <- shp
      map_annots <- ann
      # An invisible trace over the focus centroids carries the hover detail; the
      # shapes themselves cannot hold a tooltip.
      fig <- fig %>% plotly::add_markers(
        x = tb$x, y = tb$y, type = "scatter", mode = "markers",
        marker = list(size = 12, color = "rgba(0,0,0,0.01)",
                      line = list(width = 0)),
        text = sprintf("focus %d: %d %s, Rg %.0f nm, aniso %.2f",
                       tb$focus, tb$n, hotspot_unit_label(hotspot_unit()),
                       tb$rg, tb$aniso),
        hoverinfo = "text", name = "foci", showlegend = FALSE)
    }

    plotly::layout(fig,
      paper_bgcolor = bg, plot_bgcolor = bg,
      font = list(color = fgcol),
      showlegend = FALSE,
      title = list(text = sprintf("%s - %d %s, %.1f um2, %.1f/um2",
                                  fname, o$row$n_units, hotspot_unit_label(hotspot_unit()),
                                  o$row$area_um2, o$row$density),
                   font = list(size = 13)),
      shapes = map_shapes, annotations = map_annots,
      xaxis = list(title = "x (nm)", zeroline = FALSE),
      # scaleanchor keeps the cell shape honest; a stretched footprint reads as a
      # different cell.
      yaxis = list(title = "y (nm)", zeroline = FALSE,
                   scaleanchor = "x", scaleratio = 1))
  }
  output$hotspot_map <- renderPlotly(build_hotspot_map_fig())

  # ── The enhancement distribution, one curve per file ──────────────────────
  ## Cumulative rather than a histogram: the conditions are indistinguishable over
  ## the bulk of the distribution and separate only in the upper tail, so a
  ## cumulative curve shows where they separate without a bin-width choice.
  ##
  ## TWO THINGS ABOUT THE x AXIS, both learned from the real data.
  ##
  ## 1. THE LOW TAIL IS REAL BUT UNINFORMATIVE, AND IT DOMINATES A LOG AXIS.
  ##    An isolated unit in a sparse corner has a large k-th neighbour distance, so
  ##    its enhancement runs down to ~0.05x. Left alone the axis spans nearly four
  ##    decades and about 40% of the width carries the bottom 1% of the units,
  ##    which is where the plot "wants to go far below anything reasonable". The
  ##    lower limit is therefore the smallest 1st percentile across the drawn
  ##    files, floored at 0.1x. Nothing is dropped from the curves - only the view
  ##    is clipped - and the floor is stated on the axis title.
  ##
  ## 2. THE SEPARATION IS ENTIRELY IN THE UPPER TAIL, which a cumulative curve
  ##    compresses into its last few percent. The "upper tail" mode plots the
  ##    SURVIVAL function 1 - F on log-log instead, where "1 unit in 20 is at 10x
  ##    or denser" is a readable distance between curves rather than two lines
  ##    converging on 1.
  build_hotspot_ecdf_fig <- function() {
    det <- hotspot_detail()
    shiny::validate(shiny::need(!is.null(det) && length(det) > 0, "Run this panel first."))
    dark <- isTRUE(input$dark_mode)
    bg   <- if (dark) "#1e1e1e" else "white"
    fg   <- if (dark) "#e6e6e6" else "#333a45"
    tail_mode <- !identical(input$hot_curve_mode %||% "tail", "cdf")
    ulab <- hotspot_unit_label(hotspot_unit())

    grouped <- (input$analysis_by %||% "file") == "group"
    gm  <- file_group_map()
    # Resolve every file's label FIRST, "Ungrouped" included, and build the palette
    # from that complete set. Building it from the grouped files only left an
    # ungrouped file asking for a colour that had never been created.
    labs <- vapply(names(det), function(nm)
      if (grouped) .name_or(gm, nm, "Ungrouped") else nm, character(1))
    pal <- series_palette(unique(labs))
    thresh <- max(1, as.numeric(input$hot_enh_thresh %||% 5))

    # One shared x window for every curve, so the files are visually comparable.
    p1s <- vapply(det, function(o) {
      e <- o$enh[is.finite(o$enh) & o$enh > 0]
      if (length(e) < 3L) NA_real_ else unname(stats::quantile(e, 0.01))
    }, numeric(1))
    x_lo <- max(0.1, min(p1s, na.rm = TRUE))
    x_hi <- max(vapply(det, function(o) max(o$enh[is.finite(o$enh)], na.rm = TRUE),
                       numeric(1)), na.rm = TRUE)
    if (!is.finite(x_lo) || x_lo <= 0) x_lo <- 0.1
    if (!is.finite(x_hi) || x_hi <= x_lo) x_hi <- x_lo * 100

    fig <- plotly::plot_ly()
    seen <- character(0)
    for (nm in names(det)) {
      e <- sort(det[[nm]]$enh)
      e <- e[is.finite(e) & e > 0]
      if (length(e) < 3L) next
      lab <- unname(labs[nm])
      col <- .name_or(pal, lab, "#888888")
      F_e <- seq_along(e) / length(e)
      if (tail_mode) {
        # Survival, and drop the final point: 1 - F is exactly 0 there, which a log
        # axis cannot draw and which would clip the curve short.
        y <- 1 - F_e
        k <- y > 0
        xx <- e[k]; yy <- y[k]
        htm <- paste0(nm, "<br>%{y:.3f} of ", ulab,
                      " at %{x:.2f}x or denser<extra></extra>")
      } else {
        xx <- e; yy <- F_e
        htm <- paste0(nm, "<br>%{y:.2f} of ", ulab,
                      " below %{x:.2f}x<extra></extra>")
      }
      if (!length(xx)) next
      fig <- fig %>% plotly::add_lines(
        x = xx, y = yy, name = lab,
        legendgroup = lab, showlegend = !(lab %in% seen),
        line = list(color = col, width = 2), hovertemplate = htm)
      seen <- c(seen, lab)
    }

    # 1x (this cell's own average) and the reporting threshold, as references.
    vline <- function(at, col, dash) list(
      type = "line", xref = "x", yref = "paper",
      x0 = log10(at), x1 = log10(at), y0 = 0, y1 = 1,
      line = list(dash = dash, width = 1, color = col))

    # EXPLICIT TICKS. Plotly's own log ticks come out as 5e-02 / 5e-01 / 5e+00,
    # which reads as scientific notation rather than as "half the cell average" /
    # "five times the cell average" - the quantity this axis carries. A 1-2-5
    # ladder labelled with an x is what the number means.
    ladder <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
    tk <- ladder[ladder >= x_lo * 0.999 & ladder <= x_hi * 1.15]
    if (length(tk) < 3L) tk <- ladder
    xax <- list(type = "log", range = c(log10(x_lo), log10(x_hi * 1.15)),
                tickvals = tk,
                ticktext = ifelse(tk < 1, sprintf("%.1fx", tk), sprintf("%gx", tk)),
                title = sprintf("local density / cell mean (view from %.1fx)", x_lo))

    # In tail mode the y axis is also a fraction, so give it the same treatment.
    yfr <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
    yax <- if (tail_mode)
      list(type = "log", tickvals = yfr, ticktext = sprintf("%g", yfr),
           title = sprintf("fraction of %s at least this dense", ulab))
    else
      list(title = sprintf("cumulative fraction of %s", ulab), range = c(0, 1))

    plotly::layout(fig,
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      title = list(text = sprintf("Local-density enhancement (%s)", ulab),
                   font = list(size = 13)),
      legend = list(orientation = "h", y = -0.22),
      xaxis = xax, yaxis = yax,
      shapes = list(vline(1, "#888888", "dot"), vline(thresh, "#c0392b", "dash")))
  }
  output$hotspot_ecdf <- renderPlotly(build_hotspot_ecdf_fig())

  # ── One metric per file, as a dot plot by group ────────────────────────────
  .HOTSPOT_METRICS <- c(
    "Density unevenness (x random)"      = "cv_ratio",
    "CV of local density (raw)"          = "cv_local",
    "Largest focus (x random 95th pct)"  = "max_focus_ratio",
    "Largest focus (n)"                  = "max_focus_n",
    "Enhancement, 90th percentile"       = "enh_q90",
    "Enhancement, 99th percentile"       = "enh_q99",
    "Fraction above the threshold"       = "frac_above",
    "Index of dispersion (quadrats)"     = "vmr",
    "Number of foci"                     = "n_foci",
    "Fraction in foci"                   = "frac_in_foci",
    "Largest focus Rg (nm)"              = "max_focus_rg",
    "Calibrated eps (nm)"                = "eps_nm",
    "Median NN distance (nm)"            = "nn_med",
    "Clark-Evans R"                      = "clark_evans",
    "Density (per um2)"                  = "density")

  build_hotspot_metric_fig <- function() {
    rows <- hotspot_rows()
    shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "Run this panel first."))
    mkey <- input$hot_metric %||% "cv_ratio"
    shiny::validate(shiny::need(mkey %in% names(rows), "Unknown metric."))
    lab  <- names(.HOTSPOT_METRICS)[match(mkey, .HOTSPOT_METRICS)]
    if (!length(lab) || is.na(lab)) lab <- mkey
    dark <- isTRUE(input$dark_mode)
    bg   <- if (dark) "#1e1e1e" else "white"
    fg   <- if (dark) "#e6e6e6" else "#333a45"

    gm  <- file_group_map()
    grp <- unname(gm[rows$dataset]); grp[is.na(grp)] <- "Ungrouped"
    lv  <- unique(grp)
    pal <- series_palette(lv)
    v   <- rows[[mkey]]

    fig <- plotly::plot_ly()
    for (g in lv) {
      s <- grp == g
      # Jitter only on x, so the value stays exactly readable off the y axis.
      fig <- fig %>% plotly::add_markers(
        x = match(g, lv) + stats::runif(sum(s), -0.09, 0.09), y = v[s],
        name = g, marker = list(size = 11, color = .name_or(pal, g, "#888"),
                                line = list(width = 1, color = bg)),
        text = rows$dataset[s],
        hovertemplate = paste0("%{text}<br>", lab, ": %{y:.3g}<extra></extra>"))
      # Median bar: one line, no box - with 4-6 cells a boxplot draws quartiles
      # from too few points to mean anything.
      md <- stats::median(v[s], na.rm = TRUE)
      if (is.finite(md))
        fig <- fig %>% plotly::add_segments(
          x = match(g, lv) - 0.22, xend = match(g, lv) + 0.22, y = md, yend = md,
          line = list(color = fg, width = 3), showlegend = FALSE, hoverinfo = "skip")
    }
    # Three of these metrics have a known value for an unclustered pattern, so
    # draw it. Without the line the reader has no way to tell "clustered" from
    # "large number".
    ref <- switch(mkey, cv_ratio = 1, vmr = 1, clark_evans = 1, max_focus_ratio = 1,
                  cv_local = csr_cv_local_density(input$hot_k %||% 5), NULL)
    if (!is.null(ref) && is.finite(ref))
      fig <- fig %>% plotly::add_segments(
        x = 0.4, xend = length(lv) + 0.6, y = ref, yend = ref,
        line = list(color = "#888888", width = 1, dash = "dash"),
        showlegend = FALSE, hoverinfo = "skip") %>%
        plotly::add_annotations(x = length(lv) + 0.55, y = ref, text = "random",
          showarrow = FALSE, xanchor = "right", yanchor = "bottom",
          font = list(size = 10, color = "#888888"))
    plotly::layout(fig,
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      title = list(text = lab, font = list(size = 13)), showlegend = FALSE,
      xaxis = list(title = "", tickvals = seq_along(lv), ticktext = lv,
                   range = c(0.4, length(lv) + 0.6)),
      yaxis = list(title = lab))
  }
  output$hotspot_metric <- renderPlotly(build_hotspot_metric_fig())

  # ── The metric-versus-density plot ─────────────────────────────────────────
  ## THE POINT OF THIS PLOT: with transient transfection every cell carries a
  ## different plasmid load, so a metric that rises with density is reporting
  ## expression, not organisation. Plot it and look. If the groups separate
  ## VERTICALLY at matched density the effect is real; if both groups lie on one
  ## line the metric is a density readout. Median NN does exactly that, which is
  ## why it is in the metric list - it is the worked example of the failure.
  build_hotspot_density_fig <- function() {
    rows <- hotspot_rows()
    shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "Run this panel first."))
    mkey <- input$hot_metric %||% "cv_ratio"
    # A metric restored from an older settings file may name a column this build no
    # longer has; without the check the plot dies on rows[[mkey]] instead of saying so.
    shiny::validate(shiny::need(mkey %in% names(rows), "Unknown metric."))
    lab  <- names(.HOTSPOT_METRICS)[match(mkey, .HOTSPOT_METRICS)]
    if (!length(lab) || is.na(lab)) lab <- mkey
    dark <- isTRUE(input$dark_mode)
    bg   <- if (dark) "#1e1e1e" else "white"
    fg   <- if (dark) "#e6e6e6" else "#333a45"

    gm  <- file_group_map()
    grp <- unname(gm[rows$dataset]); grp[is.na(grp)] <- "Ungrouped"
    lv  <- unique(grp); pal <- series_palette(lv)
    v   <- rows[[mkey]]

    fig <- plotly::plot_ly()
    for (g in lv) {
      s <- grp == g
      fig <- fig %>% plotly::add_markers(
        x = rows$density[s], y = v[s], name = g,
        marker = list(size = 11, color = .name_or(pal, g, "#888"),
                      line = list(width = 1, color = bg)),
        text = rows$dataset[s],
        hovertemplate = "%{text}<br>%{x:.1f}/um2<br>%{y:.3g}<extra></extra>")
      # Trend line per group, drawn only with 3+ cells - with two points a line is
      # not a trend, it is a line through two points.
      if (sum(s & is.finite(v)) >= 3L) {
        ld <- log10(rows$density[s]); yy <- v[s]
        ok <- is.finite(ld) & is.finite(yy)
        if (sum(ok) >= 3L && stats::sd(ld[ok]) > 0) {
          cf <- stats::coef(stats::lm(yy[ok] ~ ld[ok]))
          xr <- range(rows$density[s][ok])
          xs <- seq(xr[1], xr[2], length.out = 20)
          fig <- fig %>% plotly::add_lines(
            x = xs, y = cf[1] + cf[2] * log10(xs), showlegend = FALSE,
            line = list(color = .name_or(pal, g, "#888"), width = 1, dash = "dash"),
            hoverinfo = "skip")
        }
      }
    }
    ref <- switch(mkey, cv_ratio = 1, vmr = 1, clark_evans = 1, max_focus_ratio = 1,
                  cv_local = csr_cv_local_density(input$hot_k %||% 5), NULL)
    if (!is.null(ref) && is.finite(ref)) {
      xr <- range(rows$density[is.finite(rows$density)])
      fig <- fig %>% plotly::add_lines(
        x = xr, y = c(ref, ref), showlegend = FALSE, hoverinfo = "skip",
        line = list(color = "#888888", width = 1, dash = "dash"))
    }
    plotly::layout(fig,
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      title = list(text = sprintf("%s vs labelling density", lab), font = list(size = 13)),
      xaxis = list(title = sprintf("%s per um2 of footprint",
                                   hotspot_unit_label(hotspot_unit())), type = "log"),
      yaxis = list(title = lab))
  }
  output$hotspot_vs_density <- renderPlotly(build_hotspot_density_fig())

  # ── The per-file table ─────────────────────────────────────────────────────
  output$hotspot_table <- renderTable({
    rows <- hotspot_rows()
    shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "Run this panel first."))
    gm <- file_group_map()
    data.frame(
      Dataset   = rows$dataset,
      Group     = ifelse(is.na(gm[rows$dataset]), "-", gm[rows$dataset]),
      N         = rows$n_units,
      `Area um2` = round(rows$area_um2, 2),
      `Density`  = round(rows$density, 1),
      `Unevenness` = round(rows$cv_ratio, 2),
      `CV local` = round(rows$cv_local, 2),
      `Enh q90`  = round(rows$enh_q90, 1),
      `Frac above` = sprintf("%.1f%%", 100 * rows$frac_above),
      `Dispersion` = round(rows$vmr, 2),
      `eps nm`   = round(rows$eps_nm),
      `Foci`     = rows$n_foci,
      `Largest`  = rows$max_focus_n,
      `Largest / random` = round(rows$max_focus_ratio %||% NA_real_, 2),
      `Rg nm`    = round(rows$max_focus_rg),
      `Aniso`    = round(rows$max_focus_aniso, 2),
      `NN med`   = round(rows$nn_med),
      `Clark-Evans` = round(rows$clark_evans, 2),
      check.names = FALSE, stringsAsFactors = FALSE)
  }, striped = TRUE, hover = TRUE, spacing = "xs", width = "100%")

  output$download_hotspot <- downloadHandler(
    filename = function() paste0("local_density_hotspots_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      rows <- hotspot_rows()
      shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "No results."))
      gm <- file_group_map()
      rows$group <- unname(gm[rows$dataset])
      write.csv(rows, file, row.names = FALSE)
    })

  # Per-focus detail: one row per detected hotspot across every file, which is what
  # you need to plot a focus-size distribution rather than one number per cell.
  output$download_hotspot_foci <- downloadHandler(
    filename = function() paste0("hotspot_foci_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      det <- hotspot_detail()
      shiny::validate(shiny::need(!is.null(det) && length(det) > 0, "No results."))
      gm <- file_group_map()
      out <- do.call(rbind, lapply(names(det), function(nm) {
        tb <- det[[nm]]$foci$table
        if (is.null(tb) || !nrow(tb)) return(NULL)
        tb$dataset <- nm
        tb$group   <- unname(gm[nm]) %||% NA_character_
        tb$eps_nm  <- det[[nm]]$eps_cal$eps %||% NA_real_
        tb
      }))
      shiny::validate(shiny::need(!is.null(out) && nrow(out) > 0, "No foci detected."))
      write.csv(out, file, row.names = FALSE)
    })

  ## ══════════════════════════════════════════════════════════════════════════
  ## MULTI-PANEL VIEW  (3 x 2 slots, one analysis per slot)
  ## ══════════════════════════════════════════════════════════════════════════
  ## GOAL
  ##   Put any six of the Analysis-tab plots side by side, so a set of files can be
  ##   read across several statistics at once instead of by scrolling between
  ##   panels. Each slot has its own picker; the slot draws whatever that picker
  ##   names.
  ##
  ## HOW IT AVOIDS BEING A SECOND IMPLEMENTATION
  ##   It does not recompute anything and it owns no analysis of its own. Every
  ##   entry in the registry below is the SAME builder function the full-size panel
  ##   renders (`build_<name>_fig()`), so a slot shows exactly what its panel shows,
  ##   including its "Run this first" message when the panel has not been run. The
  ##   builders were factored out of the existing renderPlotly() calls for this;
  ##   their bodies were not otherwise touched.
  ##
  ##   CONSEQUENCE, and the thing to remember when using it: a slot is NOT
  ##   independent of its panel. It reads that panel's controls and its stored
  ##   results, so changing a setting in the panel changes the slot, and a panel
  ##   that has not been run shows a placeholder rather than an empty plot. There is
  ##   deliberately no per-slot copy of the settings - two copies of one control is
  ##   how the same statistic ends up displayed twice with different parameters.

  # id -> label shown in the picker, and the builder that draws it. The `section`
  # entry only groups the dropdown, so a long list stays navigable.
  .PANEL_SOURCES <- list(
    none              = list(label = "(empty)",                          section = NULL,             fig = NULL),
    hotspot_map       = list(label = "Hotspot map (footprint + enhancement)", section = "Local density", fig = function() build_hotspot_map_fig()),
    hotspot_ecdf      = list(label = "Local-density enhancement curve",   section = "Local density",  fig = function() build_hotspot_ecdf_fig()),
    hotspot_metric    = list(label = "Hotspot metric per cell",           section = "Local density",  fig = function() build_hotspot_metric_fig()),
    hotspot_density   = list(label = "Hotspot metric vs density",         section = "Local density",  fig = function() build_hotspot_density_fig()),
    ripley_plot       = list(label = "Ripley's L(r) - r",                 section = "Spatial",        fig = function() build_ripley_plot_fig()),
    pcf_plot          = list(label = "Pair correlation g(r)",             section = "Spatial",        fig = function() build_pcf_plot_fig()),
    nn_batch_hist     = list(label = "Batch NN distance histogram",       section = "Spatial",        fig = function() build_nn_batch_hist_fig()),
    crowd_plot        = list(label = "Crowding",                          section = "Spatial",        fig = function() build_crowd_plot_fig()),
    density_plot      = list(label = "Labelling density profile",         section = "Density",        fig = function() build_density_plot_fig()),
    density_map       = list(label = "Support map",                       section = "Density",        fig = function() build_density_map_fig()),
    density_sweep_plot = list(label = "Support closing sweep",            section = "Density",        fig = function() build_density_sweep_plot_fig()),
    rog_hist          = list(label = "Radius of gyration",                section = "Median cluster", fig = function() build_rog_hist_fig()),
    size_nn_plot      = list(label = "Cluster size vs NN",                section = "Median cluster", fig = function() build_size_nn_plot_fig()),
    tid_count_hist    = list(label = "Units per median cluster",          section = "Median cluster", fig = function() build_tid_count_hist_fig()),
    det_time_plot     = list(label = "Detections over time",              section = "Time",           fig = function() build_det_time_plot_fig()),
    tid_time_plot     = list(label = "Trace IDs over time",               section = "Time",           fig = function() build_tid_time_plot_fig()),
    cross_nn_hist     = list(label = "Cross-channel NN",                  section = "Two colour",     fig = function() build_cross_nn_hist_fig()),
    coloc_curve       = list(label = "Co-localization curve",             section = "Two colour",     fig = function() build_coloc_curve_fig()),
    ripley12_plot     = list(label = "Bivariate Ripley's L12",            section = "Two colour",     fig = function() build_ripley12_plot_fig()),
    pcf12_plot        = list(label = "Bivariate g12(r)",                  section = "Two colour",     fig = function() build_pcf12_plot_fig()),
    kde_coloc_plot    = list(label = "KDE co-localization",               section = "Two colour",     fig = function() build_kde_coloc_plot_fig())
  )

  # Choices for the slot pickers, grouped by section. "(empty)" stays ungrouped and
  # first, so a slot can always be cleared.
  panel_source_choices <- function() {
    secs <- unique(Filter(Negate(is.null),
                          lapply(.PANEL_SOURCES, function(s) s$section)))
    out <- list(`(empty)` = "none")
    for (sc in secs) {
      ids <- names(.PANEL_SOURCES)[vapply(.PANEL_SOURCES,
              function(s) identical(s$section, sc), logical(1))]
      grp <- stats::setNames(as.list(ids),
                             vapply(ids, function(i) .PANEL_SOURCES[[i]]$label, character(1)))
      out[[sc]] <- grp
    }
    out
  }

  output$panel_slot_pickers <- renderUI({
    ch <- panel_source_choices()
    # Defaults spell out the workflow this panel was built for: see the pattern,
    # then its distribution, then the two statistics that read its tail.
    defaults <- c("hotspot_map", "hotspot_ecdf", "pcf_plot",
                  "hotspot_metric", "hotspot_density", "ripley_plot")
    tagList(lapply(1:6, function(i)
      column(2, selectInput(paste0("panel_slot_", i), paste("Panel", i),
                            choices = ch, selected = defaults[i], width = "100%"))))
  })

  ## One renderer per slot. Written as an explicit loop with `local()` so each
  ## iteration captures its own `i`; without it all six outputs would close over
  ## the final value of the loop variable and every slot would draw panel 6.
  for (.slot in 1:6) local({
    i  <- .slot
    id <- paste0("panel_slot_", i)
    output[[paste0(id, "_plot")]] <- renderPlotly({
      sel <- input[[id]] %||% "none"
      src <- .PANEL_SOURCES[[sel]]
      shiny::validate(shiny::need(!is.null(src) && !is.null(src$fig),
                                  sprintf("Panel %d: pick an analysis.", i)))
      src$fig()
    })
  })
  rm(.slot)

  ## The 3 x 2 grid itself. Built server-side rather than written out in the UI so
  ## the row height stays ONE control instead of six copies of the same number.
  output$panel_grid <- renderUI({
    hgt <- paste0(max(200, min(900, as.numeric(input$panel_row_height %||% 340))), "px")
    row_of <- function(idx) fluidRow(lapply(idx, function(i)
      column(4, plotlyOutput(paste0("panel_slot_", i, "_plot"), height = hgt))))
    tagList(row_of(1:3), tags$div(style = "height:10px;"), row_of(4:6))
  })

  # ── Two-Color: Cross-Channel NN ──────────────────────────────────────────
  cross_nn_results <- reactiveVal(NULL)

  observeEvent(input$run_cross_nn, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

    unit <- input$cross_nn_unit %||% "tid"
    k    <- max(1L, as.integer(input$cross_nn_k %||% 1L))
    needs_cluster  <- unit %in% c("cluster", "cluster2_cluster")
    needs_cluster2 <- unit %in% c("cluster2_tid", "cluster2_cluster", "cluster2_centroid")
    base_unit <- if (unit == "cluster2_centroid") "cluster2"
                 else if (unit %in% c("tid", "cluster2_tid")) "tid" else "cluster"

    two_color_mode <- isTRUE(input$enable_two_color)
    mplex_pair     <- mplex_cross_pair()

    results <- list()
    withProgress(message = "Cross-channel NN...", value = 0, {
      for (i in seq_along(store)) {
        incProgress(1 / length(store))
        nm <- names(store)[i]
        df <- apply_meta_filters_batch(store[[nm]], input)
        if (!nrow(df)) next

        # A chosen multiplex pair is collapsed to analysis_group 1/2 FIRST;
        # resolve_channels() then preserves that choice rather than re-deriving
        # it from channel order. cluster2 units need the chain run even when the
        # counted unit is the raw trace, so it is not conditional here.
        if (!is.null(mplex_pair))
          df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
        if (nrow(df)) df <- apply_clustering_pipeline(df, input)
        if (!nrow(df)) next
        if (!all(c(1L, 2L) %in% df$analysis_group)) next

        # ── MEDIAN-CLUSTER MODES RESTRICT, THEY DO NOT MIX ──────────────────
        # A Median Cluster belongs to ONE channel by design - the 2-colour
        # pipeline clusters each group separately, and that separation is the
        # point of the 2-colour split. So these modes do NOT pair molecules
        # inside a shared structure (no such structure exists, and building one
        # would merge the two fluorophore populations). They keep every distance
        # strictly between channels and drop the unclustered background from
        # BOTH sides, so the question asked is "how far is a group-1 STRUCTURE
        # from group-2 material" rather than "...from any stray group-2 point".
        if (needs_cluster2) {
          if (!"cluster2" %in% names(df) || !any(!is.na(df$cluster2) & df$cluster2 > 0L)) next
          df <- dplyr::filter(df, !is.na(cluster2) & cluster2 > 0L)
          if (!nrow(df) || !all(c(1L, 2L) %in% df$analysis_group)) next
        }
        d <- compute_cross_nn(df, unit = base_unit, k = k)
        if (!is.null(d) && nrow(d)) results[[nm]] <- cbind(dataset = nm, d)
      }
    })

    if (!length(results)) {
      hint <- if (needs_cluster2)
        sprintf(" One or both channels have no Median Cluster structures - check that %s produces clusters (eps/MinPts).", cluster2_label())
      else if (needs_cluster)
        sprintf(" The molecule units need %s to be enabled.", deblink_label())
      else ""
      showNotification(paste0(
        "No cross-channel data found. Check that both channels are present and analysis_group is set correctly.",
        hint), type = "warning")
      cross_nn_results(NULL)
      return()
    }
    cross_nn_results(do.call(rbind, results))
  })

  build_cross_nn_hist_fig <- function() {
    res <- cross_nn_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run cross-channel NN first."))

    aggregate   <- isTRUE(input$cross_nn_aggregate)
    hist_type   <- input$cross_nn_hist_type %||% "density"
    nbins       <- as.integer(input$cross_nn_bins %||% 40L)
    histnorm    <- if (hist_type == "density") "probability density" else ""
    col_g1      <- input$color_group1 %||% "#1f77b4"
    col_g2      <- input$color_group2 %||% "#d62728"
    grp_colors  <- c("Group 1 → Group 2" = col_g1, "Group 2 → Group 1" = col_g2)

    fig <- plotly::plot_ly()

    if (aggregate) {
      for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
        d <- res$distance[res$direction == dir]
        if (!length(d)) next
        fig <- fig %>% plotly::add_histogram(
          x = d, name = dir, histnorm = histnorm, nbinsx = nbins,
          marker = list(color = grp_colors[dir], opacity = 0.7),
          autobinx = FALSE
        )
      }
    } else {
      datasets <- unique(res$dataset)
      n_ds     <- length(datasets)
      alphas   <- if (n_ds == 1L) "cc"
                  else sprintf("%02x", as.integer(seq(204, 80, length.out = n_ds)))
      names(alphas) <- datasets
      for (nm in datasets) {
        for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
          d <- res$distance[res$dataset == nm & res$direction == dir]
          if (!length(d)) next
          base_col <- grp_colors[dir]
          fig <- fig %>% plotly::add_histogram(
            x = d, name = sprintf("%s | %s", nm, dir),
            histnorm = histnorm, nbinsx = nbins,
            marker = list(color = paste0(base_col, alphas[nm]), opacity = 0.75),
            autobinx = FALSE
          )
        }
      }
    }

    # Median marker per drawn series. The cross-NN distribution is strongly
    # right-skewed, so the median (not the mode of the tallest bin) is the number
    # that summarises it.
    # NOTE when not aggregating: one line per (dataset x direction) would be
    # unreadable on many files, so the line marks the POOLED median of each
    # direction over the files shown - which is not the median of the per-file
    # medians, and is weighted by how many points each file contributed.
    if (isTRUE(input$cross_nn_median_lines %||% TRUE)) {
      meds <- vapply(c("Group 1 → Group 2", "Group 2 → Group 1"), function(dir)
        stats::median(res$distance[res$direction == dir], na.rm = TRUE), numeric(1))
      fig <- add_median_vlines(fig, meds, grp_colors, unit = "nm")
    }

    fig %>% plotly::layout(
      barmode = "overlay",
      xaxis = list(title = "Distance (nm)"),
      yaxis = list(title = if (hist_type == "density") "Density" else "Count"),
      title = "Cross-Channel Nearest-Neighbour Distances",
      paper_bgcolor = "white", plot_bgcolor = "white"
    )
  }
  output$cross_nn_hist <- renderPlotly(build_cross_nn_hist_fig())

  build_coloc_curve_fig <- function() {
    res <- cross_nn_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run cross-channel NN first."))

    aggregate   <- isTRUE(input$cross_nn_aggregate)
    col_ch1     <- input$color_group1 %||% "#1f77b4"
    col_ch2     <- input$color_group2 %||% "#d62728"
    ch_colors   <- c("Group 1 → Group 2" = col_ch1, "Group 2 → Group 1" = col_ch2)

    fig <- plotly::plot_ly()

    add_ecdf_trace <- function(d, nm, col, dash = "solid", show_leg = TRUE) {
      if (!length(d)) return(invisible(NULL))
      ord  <- sort(d)
      frac <- seq_along(ord) / length(ord)
      fig <<- fig %>% plotly::add_trace(
        x = ord, y = frac * 100,
        type = "scatter", mode = "lines", name = nm,
        showlegend = show_leg,
        line = list(color = col, width = 2, dash = dash))
    }

    add_null_trace <- function(sub, col) {
      if (!all(c("n_ref", "area") %in% names(sub))) return(invisible(NULL))
      n_ref <- sub$n_ref[1L]; area <- sub$area[1L]
      if (!is.finite(area) || area <= 0 || n_ref < 1L) return(invisible(NULL))
      r_seq <- seq(0, max(sub$distance, na.rm = TRUE), length.out = 200L)
      null_frac <- (1 - pmax(0, 1 - pi * r_seq^2 / area)^n_ref) * 100
      fig <<- fig %>% plotly::add_trace(
        x = r_seq, y = null_frac,
        type = "scatter", mode = "lines", name = "CSR null",
        showlegend = FALSE,
        line = list(color = col, width = 1, dash = "dash"),
        opacity = 0.45)
    }

    if (aggregate) {
      for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
        sub <- res[res$direction == dir, ]
        add_ecdf_trace(sub$distance, dir, ch_colors[dir])
        add_null_trace(sub, ch_colors[dir])
      }
    } else {
      datasets <- unique(res$dataset)
      pal <- if (length(datasets) <= 8) RColorBrewer::brewer.pal(max(3, length(datasets)), "Set2")[seq_along(datasets)]
             else scales::hue_pal()(length(datasets))
      names(pal) <- datasets
      for (nm in datasets) {
        for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
          sub  <- res[res$dataset == nm & res$direction == dir, ]
          dash <- if (dir == "Group 1 → Group 2") "solid" else "dot"
          add_ecdf_trace(sub$distance, sprintf("%s | %s", nm, dir), pal[nm], dash)
          add_null_trace(sub, pal[nm])
        }
      }
    }

    fig %>% plotly::layout(
      xaxis = list(title = "Distance threshold (nm)"),
      yaxis = list(title = "% molecules with neighbour within d", range = c(0, 100)),
      title = "Co-localization Fraction  (solid = observed, dashed = CSR null)",
      hovermode = "x unified",
      paper_bgcolor = "white", plot_bgcolor = "white"
    )
  }
  output$coloc_curve <- renderPlotly(build_coloc_curve_fig())

  output$download_cross_nn <- downloadHandler(
    filename = function() paste0("cross_nn_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- cross_nn_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ── Two-Color: Bivariate Ripley's L12 ────────────────────────────────────
  ripley12_results     <- reactiveVal(NULL)
  ripley12_env_results <- reactiveVal(NULL)  # Monte-Carlo CSR envelope (per dataset)

  observeEvent(input$run_ripley12, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

    unit      <- input$ripley12_unit %||% "tid"
    use_3d    <- isTRUE(input$ripley12_use_3d)
    r_max_raw <- as.numeric(input$ripley12_r_max %||% 0)
    r_steps   <- as.integer(input$ripley12_r_steps %||% 60L)
    want_env  <- isTRUE(input$ripley12_envelope)
    nsim      <- max(19L, min(199L, as.integer(input$ripley12_nsim %||% 99L)))
    env_file_cap <- 6L
    needs_cluster <- unit == "cluster"

    get_pts_for_grp <- function(d, grp) {
      sub <- dplyr::filter(d, analysis_group == grp)
      if (!nrow(sub)) return(data.frame(x=numeric(), y=numeric(), z=numeric()))
      if (unit == "tid") {
        sub %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(x=median(x,na.rm=T), y=median(y,na.rm=T),
                           z=if("z"%in%names(sub)) median(z,na.rm=T) else 0, .groups="drop") %>%
          dplyr::filter(is.finite(x), is.finite(y))
      } else {
        if (!"cluster" %in% names(sub) || !any(!is.na(sub$cluster) & sub$cluster > 0L))
          return(data.frame(x=numeric(), y=numeric(), z=numeric()))
        sub %>%
          dplyr::filter(!is.na(cluster), cluster > 0L) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(x=median(x,na.rm=T), y=median(y,na.rm=T),
                           z=if("z"%in%names(sub)) median(z,na.rm=T) else 0, .groups="drop") %>%
          dplyr::filter(is.finite(x), is.finite(y))
      }
    }

    two_color_mode <- isTRUE(input$enable_two_color)
    mplex_pair     <- mplex_cross_pair()

    # Prep every file FIRST, so one gap-fill radius can be resolved for the whole
    # set before any curve is computed (see spatial_resolve_close()). This also
    # avoids prepping twice.
    prepped <- Filter(Negate(is.null), lapply(seq_along(store), function(i) {
      nm <- names(store)[i]
      df <- apply_meta_filters_batch(store[[nm]], input)
      if (!nrow(df)) return(NULL)
      if (!is.null(mplex_pair))
        df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
      if (nrow(df)) df <- apply_clustering_pipeline(df, input)
      if (!nrow(df) || !all(c(1L, 2L) %in% df$analysis_group)) return(NULL)
      list(df = df, tag = nm)
    }))
    if (!length(prepped)) {
      showNotification("No file has both channels after filtering.", type = "warning"); return()
    }
    cc <- spatial_gate_close(prepped, use_3d)
    common_close <- cc$close_nm

    results  <- list()
    env_list <- list()
    withProgress(message = "Bivariate Ripley L12...", value = 0, {
      for (i in seq_along(prepped)) {
        incProgress(1 / length(prepped))
        nm <- prepped[[i]]$tag
        df <- prepped[[i]]$df

        r_seq <- if (r_max_raw > 0) seq(0, r_max_raw, length.out = r_steps) else NULL

        pts1 <- get_pts_for_grp(df, 1L)
        pts2 <- get_pts_for_grp(df, 2L)
        if (!nrow(pts1) || !nrow(pts2)) next
        # Shared window for observed curve and CSR envelope (bounding box of both
        # channels' points), so the edge bias cancels in the envelope.
        axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
        win  <- apply(as.matrix(rbind(pts1[, axes, drop = FALSE],
                                      pts2[, axes, drop = FALSE])), 2, range)
        sup  <- spatial_support_of(df, use_3d, tag = nm, close_nm = common_close)
        lr <- compute_ripley_L12(
          x1=pts1$x, y1=pts1$y, z1=if(use_3d) pts1$z else NULL,
          x2=pts2$x, y2=pts2$y, z2=if(use_3d) pts2$z else NULL,
          r_seq=r_seq, use_3d=use_3d, window_ranges = win,
          support = sup, edge_correct = isTRUE(input$spatial_edge_correct)
        )
        if (is.null(lr) || !nrow(lr)) next
        results[[nm]] <- cbind(dataset = nm, lr)

        if (want_env && length(env_list) < env_file_cap) {
          env <- compute_ripley_envelope(
            n = min(nrow(pts1), 2000L), n2 = min(nrow(pts2), 2000L),
            window_ranges = win, r_seq = lr$r, use_3d = use_3d,
            nsim = nsim, estimator = "L12", support = sup,
            edge_correct = isTRUE(input$spatial_edge_correct),
            null_mode = input$coloc_null_mode %||% "labels", pts1 = pts1, pts2 = pts2)
          if (!is.null(env)) {
            gp <- envelope_global_p(env, lr$Lr_minus_r)
            env$global_p <- gp$p %||% NA_real_
            env$dataset <- nm; env_list[[nm]] <- env
          }
        }
      }
    })

    if (!length(results)) {
      showNotification("No results. Ensure both channels are present.", type = "warning")
      ripley12_results(NULL)
      ripley12_env_results(NULL)
      return()
    }
    ripley12_results(do.call(rbind, results))
    ripley12_env_results(if (length(env_list)) do.call(rbind, env_list) else NULL)
    spatial_report_windows(vapply(prepped, function(it) as.character(it$tag), character(1)))
  })

  build_ripley12_plot_fig <- function() {
    res <- ripley12_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run Bivariate Ripley L12 first."))

    aggregate <- isTRUE(input$ripley12_aggregate)
    do_smooth <- isTRUE(input$ripley12_smooth)
    span_val  <- as.numeric(input$ripley12_smooth_span %||% 0.2)
    env_res   <- ripley12_env_results()
    smooth_y  <- function(sub) {
      if (!do_smooth || nrow(sub) < 5L) return(sub$Lr_minus_r)
      tryCatch(predict(loess(Lr_minus_r ~ r, data = sub, span = span_val)),
               error = function(e) sub$Lr_minus_r)
    }
    hex2rgba <- function(hex, a) {
      cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
    }
    fig <- plotly::plot_ly()

    if (aggregate) {
      agg <- res %>%
        dplyr::group_by(r) %>%
        dplyr::summarise(Lr_minus_r = mean(Lr_minus_r, na.rm = TRUE), .groups = "drop")
      fig <- fig %>% plotly::add_trace(
        x = agg$r, y = smooth_y(agg),
        type = "scatter", mode = "lines",
        name = "All datasets (mean)",
        line = list(color = "#7B2D8B", width = 2)
      )
    } else {
      datasets <- unique(res$dataset)
      pal <- if (length(datasets) <= 8) RColorBrewer::brewer.pal(max(3, length(datasets)), "Set2")[seq_along(datasets)]
             else scales::hue_pal()(length(datasets))
      names(pal) <- datasets
      for (nm in datasets) {
        sub <- res[res$dataset == nm, ]
        if (!is.null(env_res) && nm %in% env_res$dataset) {
          e <- env_res[env_res$dataset == nm, ]; e <- e[order(e$r), ]
          fig <- fig %>%
            plotly::add_trace(x = e$r, y = e$hi, type = "scatter", mode = "lines",
              line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
            plotly::add_trace(x = e$r, y = e$lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
              legendgroup = nm, name = paste0(nm, " CSR 95% envelope"),
              showlegend = TRUE, hoverinfo = "skip") %>%
            plotly::add_trace(x = e$r, y = e$null_mean, type = "scatter", mode = "lines",
              line = list(color = unname(pal[nm]), width = 1, dash = "dot"),
              legendgroup = nm, showlegend = FALSE, hoverinfo = "skip")
          fig <- .add_global_band(fig, e, unname(pal[nm]), nm)
        }
        fig <- fig %>% plotly::add_trace(
          x = sub$r, y = smooth_y(sub),
          type = "scatter", mode = "lines", name = nm, legendgroup = nm,
          line = list(color = pal[nm], width = 2)
        )
      }
    }

    csr_label <- if (!is.null(env_res) && !aggregate)
      "theoretical CSR (uncorrected - see envelope)" else "Random (L−r = 0)"
    fig %>%
      plotly::add_trace(
        x = range(res$r), y = c(0, 0),
        type = "scatter", mode = "lines",
        name = csr_label,
        line = list(color = "black", width = 1, dash = "dash"),
        showlegend = TRUE
      ) %>%
      plotly::layout(
        xaxis = list(title = "r (nm)"),
        yaxis = list(title = "L₁₂(r) − r (nm)"),
        annotations = .global_p_annot(env_res),
        title = "Bivariate Ripley’s L₁₂(r)",
        hovermode = "x unified",
        paper_bgcolor = "white", plot_bgcolor = "white"
      )
  }
  output$ripley12_plot <- renderPlotly(build_ripley12_plot_fig())

  output$download_ripley12 <- downloadHandler(
    filename = function() paste0("ripley_L12_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- ripley12_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ── shared helper: build per-group point sets for two-color analyses ─────
  get_tc_pts <- function(df, grp, unit) {
    sub <- dplyr::filter(df, analysis_group == grp)
    if (!nrow(sub)) return(data.frame(x = numeric(), y = numeric(), z = numeric()))
    if (unit == "tid") {
      sub %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(x = median(x, na.rm = TRUE), y = median(y, na.rm = TRUE),
                         z = if ("z" %in% names(sub)) median(z, na.rm = TRUE) else 0,
                         .groups = "drop") %>%
        dplyr::filter(is.finite(x), is.finite(y))
    } else {
      if (!"cluster" %in% names(sub) || !any(!is.na(sub$cluster) & sub$cluster > 0L))
        return(data.frame(x = numeric(), y = numeric(), z = numeric()))
      sub %>%
        dplyr::filter(!is.na(cluster), cluster > 0L) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(x = median(x, na.rm = TRUE), y = median(y, na.rm = TRUE),
                         z = if ("z" %in% names(sub)) median(z, na.rm = TRUE) else 0,
                         .groups = "drop") %>%
        dplyr::filter(is.finite(x), is.finite(y))
    }
  }

  # ── shared helper: prepare data for one file in two-color batch ───────────
  prep_tc_df <- function(raw, input, needs_cluster) {
    df <- apply_meta_filters_batch(raw, input)
    if (!nrow(df)) return(NULL)
    mplex_pair <- mplex_cross_pair()
    if (!is.null(mplex_pair))
      df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
    if (nrow(df)) df <- apply_clustering_pipeline(df, input)
    if (!nrow(df) || !all(c(1L, 2L) %in% df$analysis_group)) return(NULL)
    df
  }

  # ── Pair Cross-Correlation g₁₂(r) ────────────────────────────────────────
  pcf12_results     <- reactiveVal(NULL)
  pcf12_env_results <- reactiveVal(NULL)  # Monte-Carlo CSR envelope (per dataset)

  observeEvent(input$run_pcf12, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

    unit          <- input$pcf12_unit %||% "tid"
    use_3d        <- isTRUE(input$pcf12_use_3d)
    r_max_raw     <- as.numeric(input$pcf12_r_max %||% 0)
    r_steps       <- as.integer(input$pcf12_r_steps %||% 60L)
    want_env      <- isTRUE(input$pcf12_envelope)
    nsim          <- max(19L, min(199L, as.integer(input$pcf12_nsim %||% 99L)))
    env_file_cap  <- 6L
    needs_cluster <- unit == "cluster"

    # Prep first so one gap-fill radius covers the whole set (spatial_resolve_close()).
    prepped <- Filter(Negate(is.null), lapply(seq_along(store), function(i) {
      nm <- names(store)[i]
      df <- prep_tc_df(store[[nm]], input, needs_cluster)
      if (is.null(df)) return(NULL)
      list(df = df, tag = nm)
    }))
    if (!length(prepped)) {
      showNotification("No file usable after filtering.", type = "warning"); return()
    }
    cc <- spatial_gate_close(prepped, use_3d)
    common_close <- cc$close_nm

    results  <- list()
    env_list <- list()
    withProgress(message = "Pair cross-correlation g₁₂...", value = 0, {
      for (i in seq_along(prepped)) {
        incProgress(1 / length(prepped))
        nm <- prepped[[i]]$tag
        df <- prepped[[i]]$df
        pts1 <- get_tc_pts(df, 1L, unit)
        pts2 <- get_tc_pts(df, 2L, unit)
        if (!nrow(pts1) || !nrow(pts2)) next
        r_seq <- if (r_max_raw > 0) seq(0, r_max_raw, length.out = r_steps) else NULL
        axes  <- if (use_3d) c("x", "y", "z") else c("x", "y")
        win   <- apply(as.matrix(rbind(pts1[, axes, drop = FALSE],
                                       pts2[, axes, drop = FALSE])), 2, range)
        sup   <- spatial_support_of(df, use_3d, tag = nm, close_nm = common_close)
        g <- compute_pcf12(pts1$x, pts1$y, if (use_3d) pts1$z else NULL,
                           pts2$x, pts2$y, if (use_3d) pts2$z else NULL,
                           r_seq = r_seq, use_3d = use_3d, window_ranges = win,
                           support = sup,
                           edge_correct = isTRUE(input$spatial_edge_correct))
        if (!nrow(g)) next
        results[[nm]] <- cbind(dataset = nm, g)
        if (want_env && length(env_list) < env_file_cap) {
          env <- compute_ripley_envelope(
            n = min(nrow(pts1), 2000L), n2 = min(nrow(pts2), 2000L),
            window_ranges = win, r_seq = g$r, use_3d = use_3d,
            nsim = nsim, estimator = "pcf", support = sup,
            edge_correct = isTRUE(input$spatial_edge_correct),
            null_mode = input$coloc_null_mode %||% "labels", pts1 = pts1, pts2 = pts2)
          if (!is.null(env)) {
            gp <- envelope_global_p(env, g$g12)
            env$global_p <- gp$p %||% NA_real_
            env$dataset <- nm; env_list[[nm]] <- env
          }
        }
      }
    })

    if (!length(results)) {
      showNotification("No results. Check that both channels are present.", type = "warning")
      pcf12_results(NULL); pcf12_env_results(NULL); return()
    }
    pcf12_results(do.call(rbind, results))
    pcf12_env_results(if (length(env_list)) do.call(rbind, env_list) else NULL)
    spatial_report_windows(vapply(prepped, function(it) as.character(it$tag), character(1)))
  })

  build_pcf12_plot_fig <- function() {
    res <- pcf12_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run g₁₂ first."))

    aggregate <- isTRUE(input$pcf12_aggregate)
    env_res   <- pcf12_env_results()
    hex2rgba  <- function(hex, a) {
      cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
    }
    fig <- plotly::plot_ly()

    if (aggregate) {
      agg <- res %>%
        dplyr::group_by(r) %>%
        dplyr::summarise(g12 = mean(g12, na.rm = TRUE), .groups = "drop")
      fig <- fig %>% plotly::add_trace(x = agg$r, y = agg$g12,
        type = "scatter", mode = "lines", name = "All datasets (mean)",
        line = list(color = "#7B2D8B", width = 2))
    } else {
      datasets <- unique(res$dataset)
      pal <- if (length(datasets) <= 8) RColorBrewer::brewer.pal(max(3, length(datasets)), "Set2")[seq_along(datasets)]
             else scales::hue_pal()(length(datasets))
      names(pal) <- datasets
      for (nm in datasets) {
        sub <- res[res$dataset == nm, ]
        if (!is.null(env_res) && nm %in% env_res$dataset) {
          e <- env_res[env_res$dataset == nm, ]; e <- e[order(e$r), ]
          fig <- fig %>%
            plotly::add_trace(x = e$r, y = e$hi, type = "scatter", mode = "lines",
              line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
            plotly::add_trace(x = e$r, y = e$lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
              legendgroup = nm, name = paste0(nm, " CSR 95% envelope"),
              showlegend = TRUE, hoverinfo = "skip") %>%
            plotly::add_trace(x = e$r, y = e$null_mean, type = "scatter", mode = "lines",
              line = list(color = unname(pal[nm]), width = 1, dash = "dot"),
              legendgroup = nm, showlegend = FALSE, hoverinfo = "skip")
          fig <- .add_global_band(fig, e, unname(pal[nm]), nm)
        }
        fig <- fig %>% plotly::add_trace(x = sub$r, y = sub$g12,
          type = "scatter", mode = "lines", name = nm, legendgroup = nm,
          line = list(color = pal[nm], width = 2))
      }
    }

    fig %>%
      plotly::add_trace(x = range(res$r, na.rm = TRUE), y = c(1, 1),
        type = "scatter", mode = "lines", name = "CSR (g = 1)",
        line = list(color = "black", width = 1, dash = "dash"), showlegend = TRUE) %>%
      plotly::layout(
        xaxis = list(title = "r (nm)"),
        yaxis = list(title = "g₁₂(r)"),
        annotations = .global_p_annot(env_res),
        title = "Pair Cross-Correlation Function g₁₂(r)",
        hovermode = "x unified",
        paper_bgcolor = "white", plot_bgcolor = "white"
      )
  }
  output$pcf12_plot <- renderPlotly(build_pcf12_plot_fig())

  output$download_pcf12 <- downloadHandler(
    filename = function() paste0("pcf12_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      res <- pcf12_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ── 2D KDE Density Co-localization Map (active file only) ───────────────
  kde_coloc_results <- reactiveVal(NULL)

  output$kde_plane_ui <- renderUI({
    nm  <- input$active_file %||% ""
    raw <- if (nzchar(nm)) smlm_file_store()[[nm]] else NULL
    has_z <- !is.null(raw) && "z" %in% names(raw) && any(is.finite(raw$z))
    choices <- c("XY plane" = "xy")
    if (has_z) choices <- c(choices, "XZ plane" = "xz", "YZ plane" = "yz")
    selectInput("kde_plane", "Projection plane", choices = choices, selected = "xy")
  })

  observeEvent(input$run_kde_coloc, {
    nm <- input$active_file %||% ""
    if (!nzchar(nm)) { showNotification("No active file selected.", type = "warning"); return() }
    store <- smlm_file_store()
    raw   <- store[[nm]]
    if (is.null(raw) || !nrow(raw)) { showNotification("Active file has no data.", type = "warning"); return() }

    unit          <- input$kde_unit %||% "tid"
    bw_raw        <- as.numeric(input$kde_bandwidth %||% 0)
    n_grid        <- as.integer(input$kde_n_grid %||% 100L)
    do_sweep      <- isTRUE(input$kde_bw_sweep)
    costes_n      <- max(0L, as.integer(input$kde_costes_nsim %||% 0L))
    needs_cluster <- unit == "cluster"

    df <- prep_tc_df(raw, input, needs_cluster)
    if (is.null(df)) {
      showNotification("Could not assign two-color groups. Check 2-color settings.", type = "warning")
      kde_coloc_results(NULL); return()
    }

    pts1 <- get_tc_pts(df, 1L, unit)
    pts2 <- get_tc_pts(df, 2L, unit)
    if (!nrow(pts1) || !nrow(pts2)) {
      showNotification("One or both channels empty after filtering.", type = "warning")
      kde_coloc_results(NULL); return()
    }

    plane <- input$kde_plane %||% "xy"
    ax <- switch(plane,
      xz = c("x", "z"), yz = c("y", "z"), c("x", "y"))
    has_ax <- function(pts) all(ax %in% names(pts)) &&
      any(is.finite(pts[[ax[1]]])) && any(is.finite(pts[[ax[2]]]))
    if (!has_ax(pts1) || !has_ax(pts2)) {
      showNotification("Selected projection plane requires z coordinates which are absent.", type = "warning")
      kde_coloc_results(NULL); return()
    }

    bw  <- if (bw_raw > 0) bw_raw else NULL

    # Bandwidth sweep grid: log-spaced fractions of the field span (0.4%–20%),
    # so the bandwidth-dependence of Pearson r is visible.
    bw_sweep <- NULL
    if (do_sweep) {
      a1 <- c(pts1[[ax[1]]], pts2[[ax[1]]]); a2 <- c(pts1[[ax[2]]], pts2[[ax[2]]])
      span <- min(diff(range(a1, na.rm = TRUE)), diff(range(a2, na.rm = TRUE)))
      if (is.finite(span) && span > 0)
        bw_sweep <- span * exp(seq(log(0.004), log(0.2), length.out = 12))
    }

    # Single-file panel (active file only), so there is no set to share a gap-fill
    # radius with and the per-file plateau is the right choice here.
    kde_sup <- spatial_support_of(df, FALSE)
    res <- withProgress(message = "KDE co-localization...", value = 0.3,
      compute_kde_coloc(pts1[[ax[1]]], pts1[[ax[2]]],
                        pts2[[ax[1]]], pts2[[ax[2]]], bw, n_grid,
                        bw_sweep = bw_sweep, costes_nsim = costes_n,
                        null_mode = input$kde_null_mode %||% "shift",
                        support = kde_sup))
    if (is.null(res)) {
      showNotification("KDE computation failed.", type = "warning")
      kde_coloc_results(NULL); return()
    }
    res$dataset <- nm
    res$ax      <- ax
    # Whether the randomisation ran inside the sample footprint. With edge correction
    # on and the footprint switch off, the window passed in is still the bounding
    # box, so the switch (and an accepted footprint) is what decides the reading.
    res$footprint_null <- isTRUE(input$spatial_use_support) && !is.null(kde_sup) &&
                          !isTRUE(kde_sup$rejected_footprint)
    kde_coloc_results(res)
  })

  build_kde_coloc_plot_fig <- function() {
    r <- kde_coloc_results()
    shiny::validate(shiny::need(!is.null(r), "Run KDE first."))

    ax     <- r$ax %||% c("x", "y")
    xlabel <- paste0(ax[1], " (nm)")
    ylabel <- paste0(ax[2], " (nm)")

    style   <- input$kde_style    %||% "lines"
    cpair   <- input$kde_colorpair %||% "redblue"
    opacity <- as.numeric(input$kde_opacity %||% 0.75)

    color_pairs <- list(
      redblue      = list(ch1 = list(c(0,"white"), c(1,"#e41a1c")),
                          ch2 = list(c(0,"white"), c(1,"#377eb8"))),
      orangeteal   = list(ch1 = list(c(0,"white"), c(1,"#f46d43")),
                          ch2 = list(c(0,"white"), c(1,"#1a9988"))),
      magentagreen = list(ch1 = list(c(0,"white"), c(1,"#d01c8b")),
                          ch2 = list(c(0,"white"), c(1,"#4dac26"))),
      yellowpurple = list(ch1 = list(c(0,"white"), c(1,"#fc8d59")),
                          ch2 = list(c(0,"white"), c(1,"#542788")))
    )
    pal <- color_pairs[[cpair]] %||% color_pairs$redblue

    op <- if (style == "lines") 1.0 else opacity

    plotly::plot_ly() %>%
      plotly::add_contour(x = r$x, y = r$y, z = t(r$z1),
        name = "Ch1", opacity = op,
        colorscale = pal$ch1,
        showscale = FALSE, contours = list(coloring = style),
        line = list(width = if (style == "lines") 1.5 else 0.5)) %>%
      plotly::add_contour(x = r$x, y = r$y, z = t(r$z2),
        name = "Ch2", opacity = op,
        colorscale = pal$ch2,
        showscale = FALSE, contours = list(coloring = style),
        line = list(width = if (style == "lines") 1.5 else 0.5)) %>%
      plotly::layout(
        xaxis = list(title = xlabel),
        yaxis = list(title = ylabel, scaleanchor = "x", scaleratio = 1),
        title = sprintf("KDE - %s  [%s]  (Ch1 / Ch2)",
                        r$dataset, paste(ax, collapse = "/")),
        paper_bgcolor = "white", plot_bgcolor = "white"
    )
  }
  output$kde_coloc_plot <- renderPlotly(build_kde_coloc_plot_fig())

  output$kde_pearson_display <- renderUI({
    r <- kde_coloc_results()
    if (is.null(r)) return(NULL)
    costes_line <- if (!is.null(r$costes_p)) {
      sig <- if (r$costes_p < 0.05) "significant co-localization" else "not significant"
      tagList(
        tags$p(sprintf("Spatial-randomization p = %.3f  (%s; %d randomizations, null r = %.3f ± %.3f)",
                       r$costes_p, sig, r$costes_nsim, r$costes_null_mean, r$costes_null_sd),
               style = "font-size:13px; margin:2px 0; color:#555;"),
        tags$p(sprintf("Null: %s", r$null_label %||% "-"),
               style = "font-size:12px; margin:0 0 2px; color:#888;"),
        # Measured false-positive rates (simulation, see compute_kde_coloc header):
        # signal spread over the whole window -> shift 3% / csr 3%; signal confined to
        # one cell inside a larger window -> shift 53% / csr 93%, vs the 5% a calibrated
        # test should give over the bounding box. Drawn inside the footprint instead,
        # the test is conservative, so the caveat depends on which window was used.
        if (identical(r$null_mode, "csr")) tags$p(
          HTML("&#9888; The uniform-CSR null tests against complete spatial randomness, so it largely reports &ldquo;both channels are non-uniform&rdquo; and rejects far too readily. Prefer the <b>random-shift</b> null."),
          style = "font-size:12px; margin:0; color:#b26a00;")
        else tags$p(
          HTML(if (isTRUE(r$footprint_null))
                 "Null drawn inside the sample footprint, so the test is conservative: p &gt; 0.05 can be quoted as not significant. Take a positive claim from g&#8321;&#8322;(r)."
               else
                 "&#9888; Null drawn over the bounding box, so for a structure filling part of the view this rejects too readily. Switch the footprint window on, or confirm with g&#8321;&#8322;(r)."),
          style = paste0("font-size:12px; margin:0; color:",
                         if (isTRUE(r$footprint_null)) "#555;" else "#b26a00;")))
    } else NULL
    tags$div(
      style = "margin-top:6px;",
      tags$b("Density-density Pearson r:"),
      tags$p(sprintf("%.3f  (bandwidth = %.1f nm)", r$pearson_r, r$bandwidth),
             style = "font-size:13px; margin:2px 0;"),
      costes_line,
      tags$p(HTML("Pearson r is bandwidth-dependent. For a bandwidth-free, scale-resolved colocalization measure use Pair Cross-Correlation g₁₂(r)."),
             style = "font-size:12px; margin:4px 0 0; color:#888; font-style:italic;")
    )
  })

  output$kde_bw_sweep_plot <- renderPlotly({
    r <- kde_coloc_results()
    shiny::validate(shiny::need(!is.null(r) && !is.null(r$bw_sweep),
                                "Enable 'Bandwidth sweep' and run KDE to see r(bandwidth)."))
    sw <- r$bw_sweep
    plotly::plot_ly() %>%
      plotly::add_trace(x = sw$bw, y = sw$pearson_r, type = "scatter", mode = "lines+markers",
        line = list(color = "#7B2D8B", width = 2), marker = list(size = 5),
        name = "Pearson r") %>%
      plotly::add_trace(x = r$bandwidth, y = r$pearson_r, type = "scatter", mode = "markers",
        marker = list(color = "#e41a1c", size = 10, symbol = "diamond"),
        name = "chosen bandwidth") %>%
      plotly::layout(
        xaxis = list(title = "KDE bandwidth (nm)"),
        yaxis = list(title = "Density-density Pearson r"),
        title = "Bandwidth dependence of Pearson r",
        hovermode = "x unified",
        paper_bgcolor = "white", plot_bgcolor = "white")
  })

  output$download_kde_coloc <- downloadHandler(
    filename = function() paste0("kde_coloc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      r <- kde_coloc_results()
      shiny::validate(shiny::need(!is.null(r), "No results."))
      tbl <- expand.grid(x = r$x, y = r$y) %>%
        dplyr::mutate(dataset     = r$dataset,
                      density_ch1 = as.vector(t(r$z1)),
                      density_ch2 = as.vector(t(r$z2)),
                      pearson_r   = r$pearson_r)
      write.csv(tbl, file, row.names = FALSE)
    }
  )

  list(
    sd_filt           = sd_filt,
    sd_display_ranges = smlm_sd_display_ranges,
    file_store        = smlm_file_store,
    plottable_data    = plottable_data,
    uploaded_df       = uploaded_df,
    meta_full_rng     = meta_full_rng,
    sd_filtered_tbl   = sd_filtered_tbl,
    sd_all_tbl        = sd_all_tbl
  )
} # end smlm_server





