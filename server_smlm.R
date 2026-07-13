## ─────────────────────────────────────────────────────────────────────────────
## server_smlm.R  —  SMLM analysis server logic
## Called from server.R inside the combined server function.
## ─────────────────────────────────────────────────────────────────────────────

smlm_server <- function(input, output, session) {
  smlm_axis_x_span  <- reactiveVal(NULL)
  smlm_zoom_ranges  <- reactiveVal(NULL)   # list(x=c(x0,x1), y=c(y0,y1)) from 2D plot zoom
  smlm_2d_crop_rv   <- reactiveVal(NULL)   # active 2D crop box: list(x=c(x0,x1), y=c(y0,y1))

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

  # Helper: load one file from a temp path, return df or NULL
  .load_one_file <- function(path, fname) {
    ext <- get_ext(fname)
    tryCatch({
      withProgress(message = paste0("Loading ", fname, "…"), value = 0, {
        incProgress(0.05, detail = "Preparing")
        loaded <- switch(
          ext,
          "csv" = {
            incProgress(0.20, detail = "Reading CSV")
            read.csv(path, check.names = FALSE)
          },
          "json" = {
            incProgress(0.10, detail = "Reading JSON (can take a while)")
            read_minflux_json_as_df(
              path,
              keep_final_only = TRUE,
              keep_valid_only = TRUE,
              final_itr = NA_integer_,
              cfr_strategy = "last_nonzero",
              z_scale = 0.7,
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
              z_scale = 0.7,
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
        incProgress(0.05, detail = "Finalizing")
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

    for (i in seq_len(nrow(input$file))) {
      fname <- input$file$name[i]
      if (fname %in% names(store)) next   # already loaded — skip duplicate

      df <- .load_one_file(input$file$datapath[i], fname)
      if (is.null(df)) next

      df$.dataset_tag <- paste0(sample(c(0:9, letters[1:6]), 4, replace = TRUE), collapse = "")
      store[[fname]] <- df
      new_names <- c(new_names, fname)

      if (identical(get_ext(fname), "json")) {
        final_used <- if ("itr" %in% names(df)) max(df$itr, na.rm = TRUE) else NA_integer_
        showNotification(
          sprintf("%s — JSON import: using final localization (itr=%s); CFR carried over per TID.",
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
  }, ignoreInit = TRUE)

  # Switch active dataset when the file selector changes
  observeEvent(input$active_file, {
    req(input$active_file, nzchar(input$active_file))
    df <- smlm_file_store()[[input$active_file]]
    if (is.null(df)) return()

    # Snapshot sd_filt BEFORE the data change so we can restore it after the
    # reactive graph settles.  This is the authoritative lock mechanism for SD
    # ranges — more reliable than the timing-based sd_skip_arming approach.
    snap_sd <- if (isTRUE(isolate(input$lock_filters)))
      isolate(list(x = sd_filt$x, y = sd_filt$y, z = sd_filt$z))
    else NULL

    uploaded_df_rv(df)

    # Update histogram variable choices — show dcr_pooled only when present
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

  # TRUE when the active dataset is a multiplex (drives the *_by_channel UIs).
  output$mplex_active <- reactive({
    df <- uploaded_df_rv()
    !is.null(df) && "mplex_label" %in% names(df)
  })
  outputOptions(output, "mplex_active", suspendWhenHidden = FALSE)

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
  
  # --------- 2) RESET FILTERS ON NEW DATASET ----------
  observeEvent(uploaded_df(), {
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

    # Reset SD filter synchronously — no JS round-trip needed.
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
      # SD filter — save from sd_filt (source of truth), not from input$sd_x_range
      # which is unreliable when the SD panel is closed (renderUI removes the widgets).
      sd_x_range       = sd_filt$x,
      sd_y_range       = sd_filt$y,
      sd_z_range       = sd_filt$z,
      sd_min_locs      = input$sd_min_locs,
      make_positive    = input$make_positive,
      
      enable_clustering   = input$enable_clustering,
      eps                 = input$eps,
      minpts              = input$minpts,
      exclude_unclustered = input$exclude_unclustered,
      
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
      eps2 = 150,
      minpts2 = 8,
      
      # 3D wrappers
      show_cluster2_shape = FALSE,
      show_cluster2_isosurface = FALSE,
      cluster2_iso_selected_only = TRUE,
      cluster2_iso_id = 1,
      cluster2_iso_color = "#00CED1",
      cluster2_iso_unicolor = FALSE,
      cluster2_iso_grid_res = 22,
      cluster2_iso_pad = 100,
      cluster2_iso_level_q = 0.85,
      cluster2_iso_opacity = 0.25,
      cluster2_iso_sigma = 150,
      
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
      
      draw_nn_two_color = FALSE,
      nn_k_two_color = 4,
      show_nn_labels_two_color = FALSE,
      distance_range_two_color = c(0, 500),
      two_color_scope = "between",
      
      enable_two_color_cluster1 = FALSE,
      eps_tc1 = 14,
      minpts_tc1 = 6,
      show_tc1_alpha = FALSE,
      tc1_alpha = 120,
      tc1_alpha_opacity = 0.15,
      tc1_alpha_min_pts = 5,
      exclude_unclustered_tc1 = FALSE,
      
      enable_two_color_cluster2 = FALSE,
      cluster2_method_tc2 = "dbscan",
      cluster2_xi_tc2 = 0.05,
      cluster2_prefilter_tc2 = FALSE,
      cluster2_prefilter_r_tc2 = 150,
      cluster2_prefilter_min_tc2 = 3,
      cluster2_prune_tc2 = FALSE,
      cluster2_prune_tol_tc2 = 2,
      eps_tc2 = 150,
      minpts_tc2 = 8,
      show_tc2_alpha = FALSE,
      tc2_alpha = 150,
      tc2_alpha_opacity = 0.18,
      tc2_alpha_min_pts = 15,
      
      show_tc2_isosurface = FALSE,
      tc2_iso_selected_only = TRUE,
      tc2_iso_cluster_label = "A-1",
      tc2_iso_grid_res = 22,
      tc2_iso_pad = 100,
      tc2_iso_level_q = 0.85,
      tc2_iso_sigma = 150,
      tc2_iso_opacity = 0.25,
      tc2_iso_eps = 1e-4
    )
  }
  
  apply_two_color_settings <- function(session, s) {
    getv <- function(name, default = NULL) if (!is.null(s[[name]])) s[[name]] else default
    
    # checkboxes
    cb_ids <- c(
      "color_by_fluo",
      "draw_nn_two_color", "show_nn_labels_two_color",
      "enable_two_color_cluster1", "show_tc1_alpha", "exclude_unclustered_tc1",
      "enable_two_color_cluster2", "show_tc2_alpha",
      "show_tc2_isosurface", "tc2_iso_selected_only"
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
    if (!is.null(getv("nn_k_two_color"))) updateNumericInput(session, "nn_k_two_color", value = as.integer(getv("nn_k_two_color")))
    rng <- getv("distance_range_two_color")
    if (!is.null(rng) && length(rng) == 2) updateSliderInput(session, "distance_range_two_color", value = as.numeric(rng))
    if (!is.null(getv("two_color_scope"))) updateSelectInput(session, "two_color_scope", selected = as.character(getv("two_color_scope")))
    
    # clustering tc1
    if (!is.null(getv("eps_tc1"))) updateNumericInput(session, "eps_tc1", value = as.numeric(getv("eps_tc1")))
    if (!is.null(getv("minpts_tc1"))) updateNumericInput(session, "minpts_tc1", value = as.integer(getv("minpts_tc1")))
    if (!is.null(getv("tc1_alpha"))) updateNumericInput(session, "tc1_alpha", value = as.numeric(getv("tc1_alpha")))
    if (!is.null(getv("tc1_alpha_opacity"))) updateSliderInput(session, "tc1_alpha_opacity", value = as.numeric(getv("tc1_alpha_opacity")))
    if (!is.null(getv("tc1_alpha_min_pts"))) updateNumericInput(session, "tc1_alpha_min_pts", value = as.integer(getv("tc1_alpha_min_pts")))
    
    # clustering tc2
    if (!is.null(getv("cluster2_method_tc2"))) updateRadioButtons(session, "cluster2_method_tc2", selected = as.character(getv("cluster2_method_tc2")))
    if (!is.null(getv("cluster2_xi_tc2"))) updateNumericInput(session, "cluster2_xi_tc2", value = as.numeric(getv("cluster2_xi_tc2")))
    if (!is.null(getv("cluster2_prefilter_tc2"))) updateCheckboxInput(session, "cluster2_prefilter_tc2", value = isTRUE(getv("cluster2_prefilter_tc2")))
    if (!is.null(getv("cluster2_prefilter_r_tc2"))) updateNumericInput(session, "cluster2_prefilter_r_tc2", value = as.numeric(getv("cluster2_prefilter_r_tc2")))
    if (!is.null(getv("cluster2_prefilter_min_tc2"))) updateNumericInput(session, "cluster2_prefilter_min_tc2", value = as.numeric(getv("cluster2_prefilter_min_tc2")))
    if (!is.null(getv("cluster2_prune_tc2"))) updateCheckboxInput(session, "cluster2_prune_tc2", value = isTRUE(getv("cluster2_prune_tc2")))
    if (!is.null(getv("cluster2_prune_tol_tc2"))) updateNumericInput(session, "cluster2_prune_tol_tc2", value = as.numeric(getv("cluster2_prune_tol_tc2")))
    if (!is.null(getv("eps_tc2"))) updateNumericInput(session, "eps_tc2", value = as.numeric(getv("eps_tc2")))
    if (!is.null(getv("minpts_tc2"))) updateNumericInput(session, "minpts_tc2", value = as.integer(getv("minpts_tc2")))
    if (!is.null(getv("tc2_alpha"))) updateNumericInput(session, "tc2_alpha", value = as.numeric(getv("tc2_alpha")))
    if (!is.null(getv("tc2_alpha_opacity"))) updateSliderInput(session, "tc2_alpha_opacity", value = as.numeric(getv("tc2_alpha_opacity")))
    if (!is.null(getv("tc2_alpha_min_pts"))) updateNumericInput(session, "tc2_alpha_min_pts", value = as.integer(getv("tc2_alpha_min_pts")))
    
    # isosurface
    if (!is.null(getv("tc2_iso_cluster_label"))) updateTextInput(session, "tc2_iso_cluster_label", value = as.character(getv("tc2_iso_cluster_label")))
    if (!is.null(getv("tc2_iso_grid_res"))) updateSliderInput(session, "tc2_iso_grid_res", value = as.numeric(getv("tc2_iso_grid_res")))
    if (!is.null(getv("tc2_iso_pad"))) updateSliderInput(session, "tc2_iso_pad", value = as.numeric(getv("tc2_iso_pad")))
    if (!is.null(getv("tc2_iso_level_q"))) updateSliderInput(session, "tc2_iso_level_q", value = as.numeric(getv("tc2_iso_level_q")))
    if (!is.null(getv("tc2_iso_sigma"))) updateSliderInput(session, "tc2_iso_sigma", value = as.numeric(getv("tc2_iso_sigma")))
    if (!is.null(getv("tc2_iso_opacity"))) updateSliderInput(session, "tc2_iso_opacity", value = as.numeric(getv("tc2_iso_opacity")))
    if (!is.null(getv("tc2_iso_eps"))) updateSliderInput(session, "tc2_iso_eps", value = as.numeric(getv("tc2_iso_eps")))
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
    updateCheckboxInput(session, "sd_overlay_all", value = vals$sd_overlay_all)
    updateCheckboxInput(session, "sd_show_mean_sd", value = vals$sd_show_mean_sd %||% TRUE)
    updateCheckboxInput(session, "sd_show_overall_mean_sd", value = vals$sd_show_overall_mean_sd %||% FALSE)
    if (!is.null(vals$sd_x_range)) shinyWidgets::updateNumericRangeInput(session, "sd_x_range", value = rng(vals$sd_x_range))
    if (!is.null(vals$sd_y_range)) shinyWidgets::updateNumericRangeInput(session, "sd_y_range", value = rng(vals$sd_y_range))
    if (!is.null(vals$sd_z_range)) shinyWidgets::updateNumericRangeInput(session, "sd_z_range", value = rng(vals$sd_z_range))
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
    if (!is.null(getv("eps2")))    updateNumericInput(session, "eps2", value = as.numeric(getv("eps2")))
    if (!is.null(getv("minpts2"))) updateNumericInput(session, "minpts2", value = as.numeric(getv("minpts2")))
    
    # checkboxes
    cb_ids <- c(
      "show_cluster2_shape", "show_cluster2_isosurface", "cluster2_iso_selected_only",
      "show_cluster2_shape_2D", "show_cluster2_metrics",
      "show_cluster2_pca_axes",
      "show_cluster2_labels", "show_cluster2_plane_angle", "show_cluster2_plane_distance",
      "show_bending_radius", "fit_weighted",
      "show_fit_distances", "show_fit_tube",
      "show_cluster2_nn_lines"
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
    if (!is.null(getv("cluster2_iso_level_q")))  updateSliderInput(session, "cluster2_iso_level_q", value = as.numeric(getv("cluster2_iso_level_q")))
    if (!is.null(getv("cluster2_iso_opacity")))  updateSliderInput(session, "cluster2_iso_opacity", value = as.numeric(getv("cluster2_iso_opacity")))
    if (!is.null(getv("cluster2_iso_sigma")))    updateSliderInput(session, "cluster2_iso_sigma", value = as.numeric(getv("cluster2_iso_sigma")))
    if (!is.null(getv("cluster2_iso_color")))    colourpicker::updateColourInput(session, "cluster2_iso_color", value = as.character(getv("cluster2_iso_color")))
    if (!is.null(getv("cluster2_iso_unicolor"))) updateCheckboxInput(session, "cluster2_iso_unicolor", value = isTRUE(getv("cluster2_iso_unicolor")))
    
    # color
    if (!is.null(getv("fit_color_2d"))) colourpicker::updateColourInput(session, "fit_color_2d", value = as.character(getv("fit_color_2d")))
    
    # NN controls
    if (!is.null(getv("nn_k_cluster2"))) updateNumericInput(session, "nn_k_cluster2", value = as.integer(getv("nn_k_cluster2")))
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
  # c(0, Inf) = "no filter" — passes all traces regardless of spread.
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
  
  
  metadata_filtered_df <- reactive({
    df <- uploaded_df()
    
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

  # --------- 6) FILTER CORE (shared by 2D+3D) ----------
  # Stage 2: spread filter — per-trace SD applied after metadata row filters
  spread_filtered_df <- reactive({
    df <- metadata_filtered_df()
    if (!is.data.frame(df) || !nrow(df) || !"tid" %in% names(df)) return(df)

    # sd_filt holds the committed filter bounds set by the user.
    # c(0, Inf) upper bound = "no filter" — all traces pass regardless of spread.
    # This avoids every floating-point comparison edge case that arose from
    # matching widget values against data maxima computed from different subsets.
    rx <- c(sd_filt$x[1], sd_filt$x[2])
    ry <- c(sd_filt$y[1], sd_filt$y[2])
    rz <- c(sd_filt$z[1], sd_filt$z[2])

    if (is.infinite(rx[2]) && is.infinite(ry[2]) && is.infinite(rz[2])) return(df)

    # per-tid SDs computed from CURRENT df (after metadata filters)
    sda <- sd_by_tid(df, min_n = input$sd_min_locs)
    if (!is.data.frame(sda) || !nrow(sda)) return(df[0, , drop = FALSE])

    if (!"sd_z" %in% names(sda)) sda$sd_z <- NA_real_

    pass_x <- is.infinite(rx[2]) | is.na(sda$sd_x) | (sda$sd_x >= rx[1] & sda$sd_x <= rx[2])
    pass_y <- is.infinite(ry[2]) | is.na(sda$sd_y) | (sda$sd_y >= ry[1] & sda$sd_y <= ry[2])
    pass_z <- is.infinite(rz[2]) | is.na(sda$sd_z) | (sda$sd_z >= rz[1] & sda$sd_z <= rz[2])

    keep_tids <- sda$tid[ pass_x & pass_y & pass_z ]

    df %>% dplyr::filter(tid %in% keep_tids)
  })
  # Stage 3: trace ID length filter — applied after spread filter
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

  # Trace-ID-Length slider UI — range derived from the spread-filtered data
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
    # spread-filtered maximum — even after the spread filter is cleared.
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
        cluster2_method = input$cluster2_method %||% "dbscan",
        cluster2_xi = input$cluster2_xi,
        cluster2_prefilter = isTRUE(input$cluster2_prefilter),
        cluster2_prefilter_r = input$cluster2_prefilter_r,
        cluster2_prefilter_min = input$cluster2_prefilter_min,
        cluster2_prune = isTRUE(input$cluster2_prune),
        cluster2_prune_tol = input$cluster2_prune_tol,
        eps2 = input$eps2,
        minpts2 = input$minpts2,

        show_cluster2_shape = isTRUE(input$show_cluster2_shape),
        show_cluster2_isosurface = isTRUE(input$show_cluster2_isosurface),
        cluster2_iso_selected_only = isTRUE(input$cluster2_iso_selected_only),
        cluster2_iso_id = input$cluster2_iso_id,
        cluster2_iso_color = input$cluster2_iso_color,
        cluster2_iso_unicolor = isTRUE(input$cluster2_iso_unicolor),
        cluster2_iso_grid_res = input$cluster2_iso_grid_res,
        cluster2_iso_pad = input$cluster2_iso_pad,
        cluster2_iso_level_q = input$cluster2_iso_level_q,
        cluster2_iso_opacity = input$cluster2_iso_opacity,
        cluster2_iso_sigma = input$cluster2_iso_sigma,
        
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
        
        show_cluster2_nn_lines = isTRUE(input$show_cluster2_nn_lines),
        nn_k_cluster2 = input$nn_k_cluster2,
        distance_range_cluster2 = unname(input$distance_range_cluster2) # avoid jsonlite named-vector warning
      )
      
      jsonlite::write_json(s, path = file, auto_unbox = TRUE, pretty = TRUE)
    }
  )

  # ── kNN-distance diagnostic for tuning cluster2 eps ──
  # Plots, for the same centroid set cluster2 operates on (trace medians, or
  # primary-cluster medians when primary clustering is on), the sorted distance
  # to the (MinPts-1)-th neighbour. The "elbow" separates dense structures from
  # scattered background = a sensible eps. Red dashed line = current eps2.
  output$cluster2_knn_plot <- renderPlot({
    req(isTRUE(input$enable_cluster2))
    df <- plottable_data()
    req(!is.null(df), nrow(df) > 0)

    use_primary <- isTRUE(input$enable_clustering) && "cluster" %in% names(df)
    meds <- if (use_primary) {
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

    cols <- if (isTRUE(input$mode3d == "2D")) c("mx","my") else c("mx","my","mz")
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

  output$download_two_color_settings <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
      paste0(base, "_two_color_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
    },
    content = function(file) {
      s <- list(
        color_by_fluo = isTRUE(input$color_by_fluo),
        dcr_group1 = unname(input$dcr_group1),
        dcr_group2 = unname(input$dcr_group2),
        color_group1 = input$color_group1,
        color_group2 = input$color_group2,
        
        draw_nn_two_color = isTRUE(input$draw_nn_two_color),
        nn_k_two_color = input$nn_k_two_color,
        show_nn_labels_two_color = isTRUE(input$show_nn_labels_two_color),
        distance_range_two_color = unname(input$distance_range_two_color),
        two_color_scope = input$two_color_scope,
        
        enable_two_color_cluster1 = isTRUE(input$enable_two_color_cluster1),
        eps_tc1 = input$eps_tc1,
        minpts_tc1 = input$minpts_tc1,
        show_tc1_alpha = isTRUE(input$show_tc1_alpha),
        tc1_alpha = input$tc1_alpha,
        tc1_alpha_opacity = input$tc1_alpha_opacity,
        tc1_alpha_min_pts = input$tc1_alpha_min_pts,
        exclude_unclustered_tc1 = isTRUE(input$exclude_unclustered_tc1),
        
        enable_two_color_cluster2 = isTRUE(input$enable_two_color_cluster2),
        cluster2_method_tc2 = input$cluster2_method_tc2 %||% "dbscan",
        cluster2_xi_tc2 = input$cluster2_xi_tc2,
        cluster2_prefilter_tc2 = isTRUE(input$cluster2_prefilter_tc2),
        cluster2_prefilter_r_tc2 = input$cluster2_prefilter_r_tc2,
        cluster2_prefilter_min_tc2 = input$cluster2_prefilter_min_tc2,
        cluster2_prune_tc2 = isTRUE(input$cluster2_prune_tc2),
        cluster2_prune_tol_tc2 = input$cluster2_prune_tol_tc2,
        eps_tc2 = input$eps_tc2,
        minpts_tc2 = input$minpts_tc2,
        show_tc2_alpha = isTRUE(input$show_tc2_alpha),
        tc2_alpha = input$tc2_alpha,
        tc2_alpha_opacity = input$tc2_alpha_opacity,
        tc2_alpha_min_pts = input$tc2_alpha_min_pts,
        
        show_tc2_isosurface = isTRUE(input$show_tc2_isosurface),
        tc2_iso_selected_only = isTRUE(input$tc2_iso_selected_only),
        tc2_iso_cluster_label = input$tc2_iso_cluster_label,
        tc2_iso_grid_res = input$tc2_iso_grid_res,
        tc2_iso_pad = input$tc2_iso_pad,
        tc2_iso_level_q = input$tc2_iso_level_q,
        tc2_iso_sigma = input$tc2_iso_sigma,
        tc2_iso_opacity = input$tc2_iso_opacity,
        tc2_iso_eps = input$tc2_iso_eps
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

      # Core MINFLUX columns — always included when present
      core_cols <- c("x", "y", "z", "tid", "efo", "cfr", "eco", "fbg", "dcr")
      # Optional columns included only when they exist and carry real data
      optional_cols <- c("fluo", "tim", "itr", "dwell",
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
    }, once = TRUE)
  })
  cluster2_file_reset_key <- reactiveVal(0)
  
  output$cluster2_settings_file_ui <- renderUI({
    cluster2_file_reset_key()  
    fileInput("cluster2_settings_file", "Load Cluster2 settings (.json)", accept = c(".json"))
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
        "Cluster2 settings loaded: ",
        input$cluster2_settings_file$name
      ),
      type = "message",
      duration = 3
    )
  })
  two_color_file_reset_key <- reactiveVal(0)
  
  output$two_color_settings_file_ui <- renderUI({
    two_color_file_reset_key()
    fileInput("two_color_settings_file", "Load 2-color settings (.json)", accept = c(".json"))
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
      showNotification("Cluster2 minPts raised to 4 to keep geometry stable (convex hull / fits need ≥4 points).",
                       type = "warning", duration = 4)
    }
  })
  
  observeEvent(input$reset_filters, {
    req(input$file)
    was_view_SD <- isTRUE(isolate(input$view_SD))
    vals <- filters_defaults_from_data(uploaded_df())
    # Reset SD filter synchronously — Inf upper bound means "pass all".
    sd_filt$x <- c(0, Inf); sd_filt$y <- c(0, Inf); sd_filt$z <- c(0, Inf)
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
  # Commit the user's chosen range into sd_filt — the single source of truth.
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
    mx_x  <- get_mx("sd_x"); mx_y <- get_mx("sd_y")
    rx <- input$sd_x_range %||% c(0, Inf)
    ry <- input$sd_y_range %||% c(0, Inf)
    rz <- input$sd_z_range %||% c(0, Inf)
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
    
    els <- list(
      helpText(paste("Median Trace spread X:", round(med$x %||% NA_real_, 2))),
      helpText(paste("Median Trace spread Y:", round(med$y %||% NA_real_, 2)))
    )
    if (!is.null(med$z)) {
      els <- append(els, list(helpText(paste("Median Trace spread Z:", round(med$z %||% NA_real_, 2)))))
    }
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
    # extend the widget max to fit it — prevents clamping from silently resetting
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
      }
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
  
  output$sd_hist_x <- renderPlotly({
    req(input$view_SD)
    all <- sd_all_tbl(); cur <- sd_filtered_tbl()
    shiny::validate(shiny::need(nrow(all) > 0, "No SD data"))
    
    bins <- as.integer(input$sd_bins %||% 60)
    brk  <- sd_breaks(all$sd_x, bins)
    A    <- sd_counts(all$sd_x, brk)
    B    <- sd_counts(cur$sd_x, brk)
    
    # active filter range
    r <- input$sd_x_range %||% c(-Inf, Inf)
    in_rng <- A$centers >= r[1] & A$centers <= r[2]
    
    # clamp filtered counts to 0 outside range (so black bars disappear there)
    B$counts[!in_rng] <- 0
    
    fig <- plotly::plot_ly()
    
    if (isTRUE(input$sd_overlay_all)) {
      fig <- fig %>% plotly::add_bars(
        x = A$centers, y = A$counts, width = A$width,
        marker = list(color = "rgba(160,160,160,0.6)"),
        name = "Raw data", hoverinfo = "skip"
      )
    }
    
    # Black foreground: filtered (B already clamped)
    fig <- fig %>% plotly::add_bars(
      x = A$centers, y = B$counts, width = A$width,
      marker = list(color = "rgba(0,0,0,0.85)"),
      name = "Filtered", hoverinfo = "skip"
    )
    ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
    y_top <- if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
    
    # Filtered median ± nMAD (BLUE)
    if (isTRUE(input$sd_show_mean_sd)) {
      muF  <- stats::median(cur$sd_x, na.rm = TRUE)
      labF <- sigma_label_med_pm_sd(cur$sd_x)
      fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
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
      xaxis = list(title = "\u03C3 X (nm)"),
      yaxis = list(title = "Count")
    )
  })
  
  
  output$sd_hist_y <- renderPlotly({
    req(input$view_SD)
    all <- sd_all_tbl(); cur <- sd_filtered_tbl()
    shiny::validate(shiny::need(nrow(all) > 0, "No SD data"))
    
    bins <- as.integer(input$sd_bins %||% 60)
    brk  <- sd_breaks(all$sd_y, bins)
    A    <- sd_counts(all$sd_y, brk)
    B    <- sd_counts(cur$sd_y, brk)
    
    r <- input$sd_y_range %||% c(-Inf, Inf)
    in_rng <- A$centers >= r[1] & A$centers <= r[2]
    B$counts[!in_rng] <- 0
    
    fig <- plotly::plot_ly()
    
    if (isTRUE(input$sd_overlay_all)) {
      fig <- fig %>% plotly::add_bars(
        x = A$centers, y = A$counts, width = A$width,
        marker = list(color = "rgba(160,160,160,0.6)"),
        name = "Raw data", hoverinfo = "skip"
      )
    }
    
    # Black foreground: filtered (B already clamped)
    fig <- fig %>% plotly::add_bars(
      x = A$centers, y = B$counts, width = A$width,
      marker = list(color = "rgba(0,0,0,0.85)"),
      name = "Filtered", hoverinfo = "skip"
    )
    
    ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
    y_top <- if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
    
    # Filtered median ± nMAD (BLUE)
    if (isTRUE(input$sd_show_mean_sd)) {
      muF  <- stats::median(cur$sd_y, na.rm = TRUE)
      labF <- sigma_label_med_pm_sd(cur$sd_y)
      fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
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
      xaxis = list(title = "\u03C3 Y (nm)"),
      yaxis = list(title = "Count")
    )
  })
  
  
  output$sd_hist_z <- renderPlotly({
    req(input$view_SD, has_sd_z())
    all <- sd_all_tbl(); cur <- sd_filtered_tbl()
    shiny::validate(shiny::need("sd_z" %in% names(all) && any(is.finite(all$sd_z)), "No Z in dataset"))
    
    bins <- as.integer(input$sd_bins %||% 60)
    brk  <- sd_breaks(all$sd_z, bins)
    A    <- sd_counts(all$sd_z, brk)
    B    <- sd_counts(cur$sd_z, brk)
    
    r <- input$sd_z_range %||% c(-Inf, Inf)
    in_rng <- A$centers >= r[1] & A$centers <= r[2]
    B$counts[!in_rng] <- 0
    
    fig <- plotly::plot_ly()
    
    if (isTRUE(input$sd_overlay_all)) {
      fig <- fig %>% plotly::add_bars(
        x = A$centers, y = A$counts, width = A$width,
        marker = list(color = "rgba(160,160,160,0.6)"),
        name = "Raw data", hoverinfo = "skip"
      )
    }
    
    # Black foreground: filtered (B already clamped)
    fig <- fig %>% plotly::add_bars(
      x = A$centers, y = B$counts, width = A$width,
      marker = list(color = "rgba(0,0,0,0.85)"),
      name = "Filtered", hoverinfo = "skip"
    )
    
    ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
    y_top <- if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
    
    # Filtered median ± nMAD (BLUE)
    if (isTRUE(input$sd_show_mean_sd)) {
      muF  <- stats::median(cur$sd_z, na.rm = TRUE)
      labF <- sigma_label_med_pm_sd(cur$sd_z)
      fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
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
      xaxis = list(title = "\u03C3 Z (nm)"),
      yaxis = list(title = "Count")
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
  cluster2_table_store <- reactiveVal(NULL)
  # Cache of clustered data per file — populated whenever a file is plotted.
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

  output$view_plane_2d_ui <- renderUI({
    df <- uploaded_df()
    has_z <- !is.null(df) && "z" %in% names(df) && any(is.finite(df$z))
    choices <- if (has_z) c("XY", "XZ", "YZ") else c("XY")
    radioButtons("view_plane_2d", "Viewing plane",
                 choices = choices, inline = TRUE, selected = "XY")
  })
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

    res <- auto_detect_dcr_split(dcr_meds)
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
             "  |  Group 1 mean: ", res$mean1, " (SD ", res$sd1, ")",
             "  |  Group 2 mean: ", res$mean2, " (SD ", res$sd2, ")"),
      type = "message", duration = 8
    )
  })

  output$dcr_gmm_info <- renderText({
    res <- dcr_gmm_result()
    if (is.null(res)) return("")
    if (!isTRUE(res$success)) return(paste("GMM failed:", res$reason))
    paste0(
      "Split: ", res$split,
      " | G1 mean: ", res$mean1, " (SD ", res$sd1, ")",
      " | G2 mean: ", res$mean2, " (SD ", res$sd2, ")",
      " | n = ", res$n_traces, " traces"
    )
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

  # NN index pairs honouring the multiplex scope radio. Non-multiplex datasets
  # fall back to a global k-NN. "within" = same-channel only, "across" = strictly
  # between different channels (never within). Used by all NN-line / NN-table paths.
  .nn_pairs_scoped <- function(df_geo, coords, input) {
    k <- input$nn_k %||% 4
    mplex_on <- "mplex_channel" %in% names(df_geo) && any(!is.na(df_geo$mplex_channel))
    if (!mplex_on) return(nn_index_pairs(coords, k))
    scope <- input$nn_mplex_scope %||% "within"
    nn_index_pairs(coords, k, channel = df_geo$mplex_channel,
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
      }
    }

    if (!length(seg_x)) return(NULL)

    list(
      seg_x = seg_x, seg_y = seg_y, seg_z = seg_z,
      lab_x = lab_x, lab_y = lab_y, lab_z = lab_z, labs = labs
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
  # ---- Enhanced NN table (2D) with metadata ----
  nn_table <- reactive({
    req(plot_trigger())
    df <- plottable_data()
    
    if (isTRUE(input$enable_clustering)) {
      # Use cluster medians
      df_geo <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          id        = dplyr::first(cluster),
          tid       = dplyr::first(tid),
          x         = stats::median(x, na.rm = TRUE),
          y         = stats::median(y, na.rm = TRUE),
          cluster2  = dplyr::first(cluster2),
          color = dplyr::first(color),
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups   = "drop"
        )
      id_label <- "cluster"
    } else {
      # Use tid medians
      df_geo <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          id        = dplyr::first(tid),
          x         = stats::median(x, na.rm = TRUE),
          y         = stats::median(y, na.rm = TRUE),
          cluster   = dplyr::first(cluster),
          cluster2  = dplyr::first(cluster2),
          color = dplyr::first(color),
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups   = "drop"
        )
      id_label <- "tid"
    }

    if (nrow(df_geo) < 2) return(df_geo[0, ])

    coords <- as.matrix(df_geo[, c("x","y")])
    pm <- .nn_pairs_scoped(df_geo, coords, input)

    rows <- list()
    for (r in seq_len(nrow(pm))) {
      i <- pm[r, 1]; ni <- pm[r, 2]; dist <- pm[r, 3]
      if (dist >= input$distance_range[1] && dist <= input$distance_range[2]) {
        rows[[length(rows) + 1]] <- data.frame(
          source_id = df_geo$id[i],
          target_id = df_geo$id[ni],
          distance  = round(dist, 2),
          source_x  = df_geo$x[i],
          source_y  = df_geo$y[i],
          target_x  = df_geo$x[ni],
          target_y  = df_geo$y[ni],
          source_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[i] else NA,
          source_cluster2 = df_geo$cluster2[i],
          target_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[ni] else NA,
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
  
  # ---- Secondary cluster median-to-median NN lines (green, 2D) ----
  secondary_nn_lines <- reactive({
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
    if (nrow(df_geo) < 2) return(list())
    
    coords <- as.matrix(df_geo[, c("x","y")])
    k_eff  <- min(input$nn_k_cluster2 + 1, nrow(coords))
    nn     <- RANN::nn2(coords, k = k_eff)
    
    lines <- list()
    for (i in 1:nrow(coords)) {
      for (j in 2:k_eff) {
        ni   <- nn$nn.idx[i, j]
        dist <- as.numeric(nn$nn.dists[i, j])
        if (dist >= input$distance_range_cluster2[1] && dist <= input$distance_range_cluster2[2]) {
          nn_off <- suppressWarnings(as.numeric(input$c2_label_offset %||% 5))
          if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
          midpoint <- perp_offset_xy(coords[i, ], coords[ni, ], nn_off)
          lines[[length(lines) + 1]] <- list(
            x = c(coords[i, 1], coords[ni, 1]),
            y = c(coords[i, 2], coords[ni, 2]),
            text = round(dist, 2),
            midpoint = midpoint
          )
        }
      }
    }
    lines
  })
  
  # ---- Downloads ----
  output$download_nn_table <- downloadHandler(
    filename = function() paste0("nn_distances_", Sys.Date(), ".xlsx"),
    content  = function(file) openxlsx::write.xlsx(nn_table(), file, row.names = FALSE),
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  
  output$download_nn_table2 <- downloadHandler(
    filename = function() paste0("nn_distances_cluster2_", Sys.Date(), ".xlsx"),
    content  = function(file) openxlsx::write.xlsx(nn_table_cluster2(), file, row.names = FALSE),
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
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
  output$download_cluster2_table <- downloadHandler(
    filename = function() paste0("cluster2_metrics_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      tbl <- cluster2_table_store()
      if (is.null(tbl) || nrow(tbl) == 0) stop("No Cluster2 metrics to download.")
      
      membrane_mode <- isTRUE(input$membrane_protein_options)
      
      if (!membrane_mode) {
        # Drop membrane-specific columns if they exist
        drop_cols <- c(
          "perceived_radius", "bending_radius",
          "distance_to_plane", "angle_to_plane",
          "plane_normal_x", "plane_normal_y", "plane_normal_z"
        )
        tbl <- tbl[, setdiff(names(tbl), drop_cols), drop = FALSE]
      }
      
      openxlsx::write.xlsx(
        x         = tbl,
        file      = file,
        sheetName = "Cluster2 metrics",
        rowNames  = FALSE
      )
    }
  )
  
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
  # Use trace medians (or cluster medians when primary clustering is enabled),
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
        paste0("Confocal detection — gray range: ", mr[1], "–", mr[2],
               " | peaks: ", paste(round(as.numeric(pv), 0), collapse = ", "),
               dim_str),
        type = "message", duration = 20)
    }
    showNotification("Foci alignment applied — nudge to fine-tune.",
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
    showNotification("Boundary alignment applied — nudge to fine-tune.",
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
          "The labeling-displacement inverse (candidate-GT estimation, complex-radius merge, inward linkage displacement, stoichiometric cross-check) is a heuristic developed in this tool and validated in&nbsp;silico against the simulator — treat quantitatively with care and cite as exploratory.")))
  }

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
        guided_merge    = isTRUE(input$gt_guided_merge)
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
    txt <- sprintf("GT candidate<br>± %.0f nm%s<br>%d contributing unit(s) · %d TID / %d loc",
                   pts$pos_sd_nm, ifelse(pts$fused, " (fused: multiple GTs)", ""),
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
      return(helpText("No estimable structures (need clustered data; try enabling primary + secondary clustering)."))
    s   <- ge$summary
    m   <- ge$meta
    tot_units  <- sum(s$n_units)
    tot_gt     <- sum(s$n_gt_est)
    tot_stoich <- sum(s$n_stoich)
    tot_lo     <- sum(s$n_gt_lo);  tot_hi <- sum(s$n_gt_hi)
    any_crowd  <- any(s$crowded)
    any_mm     <- any(s$stoich_mismatch)
    lvl_lbl    <- if (identical(m$level, "primary")) "primary clusters (de-blinked)" else "raw trace IDs"

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
          HTML("The two independent counts disagree in ≥1 structure — the merge radius or the labeling assumption is off. Try the stoichiometry-guided merge, adjust the merge radius, or re-check the preset/efficiency.")),
      div(style = "font-size:12px;color:#7f8c8d;line-height:1.45;",
          HTML(sprintf(
            "Modeled from %s. Units pooled within <b>%.0f nm</b> (label-complex radius)%s into one molecule%s. Position floor <b>±%.0f nm</b> (linker, irreducible) · ~%.1f TIDs/molecule expected · detection η=%.0f%% · σ<sub>loc</sub>≈%.1f nm.<br>Marker colour = position uncertainty, combining the linker floor, localization sampling (σ<sub>loc</sub>/√n<sub>loc</sub>), a <b>support penalty</b> (fewer TIDs → less certain), inward-displacement direction error, and a fusion penalty.<br><i>Counts are conditional on the labeling calibration; positions cannot be resolved below the linker floor. The two counts bracket the honest range.%s</i>",
            lvl_lbl, m$complex_r_nm %||% 0,
            if (isTRUE(m$guided_merge)) " (stoichiometry-guided)" else "",
            if ((m$displace_nm %||% 0) > 0) sprintf(", then displaced <b>%.0f nm</b> inward toward the structure", m$displace_nm) else "",
            m$link_nm, m$m_tid, 100 * m$eta, m$sigma_loc_nm,
            ifelse(any_crowd, " Crowded structures undercount — treat those counts as lower bounds.", ""))),
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
      return(helpText("No analysable structures (need median-clusters / cluster2 with ≥4 traces)."))
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
      tags$li(HTML(sprintf("<b>Structure %d</b> — <b>%s</b> %s: %s",
                           r$cluster2, r$shape, icon_for(r$shape), detail)))
    })
    tagList(
      tags$ul(style = "margin:4px 0 8px 0;padding-left:18px;", rows),
      div(style = "font-size:12px;color:#7f8c8d;line-height:1.45;",
          HTML("Shape is <b>inferred from the point distribution</b> (you don't declare it): from the gyration-tensor eigenvalues (linearity / planarity / sphericity) plus a hollow-centre test for rings. Contour length = diameter of the backbone's minimum spanning tree, so it follows curves where the straight extent underestimates. Spacing = median nearest-neighbour distance between estimated points."),
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
    shiny::validate(shiny::need(nrow(df) > 0, "No points to display — please check filters or upload a valid file."))

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
    if (isTRUE(input$show_nn_lines) && !isTRUE(input$draw_nn_two_color) && nrow(df) > 1) {
      nn <- build_nn_segments(df, input, is3d = FALSE)
      if (!is.null(nn)) {
        ln_col <- if (isTRUE(input$dark_mode)) "white" else "black"
        fig <- fig %>% plotly::add_trace(
          x = nn$seg_x, y = nn$seg_y,
          type = "scatter", mode = "lines",
          line = list(color = ln_col, width = 1),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
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
    
    if (isTRUE(input$draw_nn_two_color) && isTRUE(input$enable_two_color)) {
      edges <- two_color_nn_edges()
      if (is.data.frame(edges) && nrow(edges) > 0) {
        col_g1 <- input$color_group1 %||% "#1f77b4"
        col_g2 <- input$color_group2 %||% "#d62728"
        edge_col <- ifelse(edges$g1 == edges$g2,
                           ifelse(edges$g1 == 1L, col_g1, col_g2),
                           "magenta")
        for (cc in unique(edge_col)) {
          sel <- which(edge_col == cc)
          fig <- fig %>% plotly::add_trace(
            x = as.numeric(t(cbind(edges$x1[sel], edges$x2[sel], NA))),
            y = as.numeric(t(cbind(edges$y1[sel], edges$y2[sel], NA))),
            type = "scatter", mode = "lines",
            line = list(color = cc, width = 2),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
        if (isTRUE(input$show_nn_labels_two_color)) {
          nn_off <- suppressWarnings(as.numeric(input$tc_label_offset %||% 5))
          if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
          nn_lfs <- suppressWarnings(as.numeric(input$tc_label_font_size %||% 12))
          if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
          offsets <- t(mapply(function(x1, y1, x2, y2)
            perp_offset_xy(c(x1, y1), c(x2, y2), nn_off),
            edges$x1, edges$y1, edges$x2, edges$y2))
          mids <- data.frame(
            mx = offsets[, 1],
            my = offsets[, 2],
            lab = sprintf("%.1f", edges$dist)
          )
          fig <- fig %>% plotly::add_trace(
            x = mids$mx, y = mids$my,
            type = "scatter", mode = "text",
            text = mids$lab,
            textfont = list(color = if (isTRUE(input$dark_mode)) "white" else "black", size = nn_lfs),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      }
    }

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
                name = paste("Cluster2", cl2),
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
    # type="scatter" (SVG) — scattergl lacks circle-open / cross symbol support.
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
  observeEvent(input$draw_nn_two_color, {
    df <- plottable_data()
    g  <- assign_two_color_group(df, input)
    e <- two_color_nn_edges()
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
  
    
  apply_trimmer <- function(df, input) {
    if (!isTRUE(input$enable_trimmer) || nrow(df) == 0) return(df)
    
    mode <- if (is.null(input$trim_mode)) "box" else input$trim_mode
    keep <- rep(TRUE, nrow(df))
    
    if (mode == "box") {
      xr <- input$trim_x_range; yr <- input$trim_y_range; zr <- input$trim_z_range
      if (!is.null(xr)) keep <- keep & df$x >= xr[1] & df$x <= xr[2]
      if (!is.null(yr)) keep <- keep & df$y >= yr[1] & df$y <= yr[2]
      if (!is.null(zr)) keep <- keep & df$z >= zr[1] & df$z <= zr[2]
      
    } else if (mode == "sphere") {
      cx <- input$trim_cx; cy <- input$trim_cy; cz <- input$trim_cz; r <- input$trim_r
      if (all(is.finite(c(cx,cy,cz,r))) && r >= 0) {
        d2 <- (df$x - cx)^2 + (df$y - cy)^2 + (df$z - cz)^2
        keep <- keep & (d2 <= r^2)
      }
      
    } else if (mode == "ellip") {
      X <- as.matrix(df[,c("x","y","z")])
      c0 <- colMeans(X, na.rm = TRUE)
      Xc <- sweep(X, 2, c0)
      pc <- tryCatch(prcomp(Xc, center = FALSE), error = function(e) NULL)
      if (!is.null(pc)) {
        Q  <- Xc %*% pc$rotation
        sdv <- apply(Q, 2, stats::sd, na.rm = TRUE)
        sdv[sdv <= 1e-12] <- 1e-12
        d2 <- rowSums((Q / matrix(sdv, nrow(Q), 3, TRUE))^2)
        thr <- stats::qchisq(ifelse(is.null(input$trim_ellip_q), 0.95, input$trim_ellip_q), df = 3)
        keep <- keep & (d2 <= thr)
      }
    }
    df[keep, , drop = FALSE]
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
  make_two_color_tid_map <- function(df, input) {
    if (!"z" %in% names(df)) df$z <- 0

    .dcr_col <- if (isTRUE(input$use_pooled_dcr %||% FALSE) && "dcr_pooled" %in% names(df))
      "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL
    med <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        dcr_med  = if (!is.null(.dcr_col)) stats::median(.data[[.dcr_col]], na.rm = TRUE) else NA_real_,
        fluo_med = if ("fluo" %in% names(df)) stats::median(fluo, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
    
    use_fluo <- isTRUE(input$color_by_fluo) && ("fluo" %in% names(df))
    
    if (use_fluo) {
      # pyMINFLUX-style: fluo is 1/2 per trace; median + round makes it robust
      med$analysis_group <- dplyr::case_when(
        is.finite(med$fluo_med) & round(med$fluo_med) == 1L ~ 1L,
        is.finite(med$fluo_med) & round(med$fluo_med) == 2L ~ 2L,
        TRUE ~ NA_integer_
      )
    } else {
      # user-defined DCR thresholds
      g1 <- input$dcr_group1 %||% c(NA_real_, NA_real_)
      g2 <- input$dcr_group2 %||% c(NA_real_, NA_real_)
      
      if (length(g1) != 2 || length(g2) != 2 || !all(is.finite(c(g1, g2)))) {
        med$analysis_group <- NA_integer_
      } else {
        med$analysis_group <- dplyr::case_when(
          is.finite(med$dcr_med) & med$dcr_med >= g1[1] & med$dcr_med <  g1[2] ~ 1L,
          is.finite(med$dcr_med) & med$dcr_med >= g2[1] & med$dcr_med <= g2[2] ~ 2L,
          TRUE ~ NA_integer_
        )
      }
    }
    
    med %>% dplyr::select(tid, analysis_group, dcr_med, fluo_med)
  }
  
  
  two_color_nn_edges <- reactive({
    req(isTRUE(input$draw_nn_two_color), isTRUE(input$enable_two_color))
    
    df <- plottable_data()
    if (!nrow(df)) return(data.frame())
    
    if (!"z" %in% names(df)) df$z <- 0
    
    # Use the canonical map (DCR by default; FLUO override if enabled + present)
    map <- make_two_color_tid_map(df, input)
    
    df <- df %>%
      dplyr::select(-dplyr::any_of("analysis_group")) %>%
      dplyr::left_join(map %>% dplyr::select(tid, analysis_group), by = "tid") %>%
      dplyr::filter(.data$analysis_group %in% c(1L, 2L))
    
    if (!nrow(df)) return(data.frame())
    
    # geometry medians: tc1 cluster > general cluster > tid
    if (isTRUE(input$enable_two_color_cluster1) && "cluster_tc1" %in% names(df)) {
      geo <- df %>%
        dplyr::filter(!is.na(cluster_tc1) & cluster_tc1 != 0L) %>%
        dplyr::group_by(cluster_tc1) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          analysis_group = {
            v <- analysis_group[!is.na(analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
          .groups = "drop"
        )
    } else if (isTRUE(input$enable_clustering) && "cluster" %in% names(df)) {
      geo <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          analysis_group = {
            v <- analysis_group[!is.na(analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
          .groups = "drop"
        )
    } else {
      geo <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          analysis_group = as.integer(round(stats::median(analysis_group, na.rm = TRUE))),
          .groups = "drop"
        )
    }
    
    geo <- geo %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z), analysis_group %in% c(1L, 2L))
    if (nrow(geo) < 2) return(data.frame())
    
    scope <- match.arg(input$two_color_scope %||% "between", c("within", "between"))
    k <- input$nn_k_two_color %||% 4
    rng   <- input$distance_range_two_color %||% c(0, Inf)
    
    build_nn_edges(geo_df = geo, k = k, scope = scope, dist_range = rng)
  })
  
  
  
  
  plottable_data <- eventReactive(plot_trigger(), {
    df <- filtered_df()
    df <- apply_trimmer(df, input)
    # 2D zoom crop — applied server-side so it works in both 2D and 3D without
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
        min_n  = input$sd_min_locs %||% MIN_LOCS_FOR_SPREAD_DEFAULT
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
    
   # Primary clustering or TID-based coloring
    if (isTRUE(input$enable_clustering)) {
      scope <- if (is.null(input$cluster_group_scope)) "all" else input$cluster_group_scope
      mode  <- if (is.null(input$grouping_mode)) "none" else input$grouping_mode
      
      df$cluster <- 0L

      mplex_on <- "mplex_channel" %in% names(df) && any(!is.na(df$mplex_channel))
      if (mplex_on) {
        # Multiplex: clusters are confined to a single channel (run DBSCAN per
        # channel, offsetting IDs so they never merge across channels).
        next_id <- 1L
        for (chv in sort(unique(df$mplex_channel[!is.na(df$mplex_channel)]))) {
          idx <- which(df$mplex_channel == chv)
          if (!length(idx)) next
          pts <- as.matrix(df[idx, c("x","y","z")])
          cl  <- dbscan::dbscan(pts, eps = input$eps, minPts = input$minpts)
          cid <- cl$cluster
          cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
          df$cluster[idx] <- cid
          if (any(cid > 0)) next_id <- max(cid, na.rm = TRUE) + 1L
        }
      } else if (scope == "all" || mode == "none") {
        pts <- as.matrix(df[, c("x","y","z")])
        cl  <- dbscan::dbscan(pts, eps = input$eps, minPts = input$minpts)
        df$cluster <- cl$cluster
      } else {
        # run DBSCAN per requested group subset
        groups_to_run <- switch(scope,
                                g1 = 1L, g2 = 2L, split = c(1L, 2L),
                                all = integer(0))
        next_id <- 1L
        for (g in groups_to_run) {
          idx <- which(df$analysis_group == g)
          if (length(idx) >= 1) {
            pts <- as.matrix(df[idx, c("x","y","z")])
            cl  <- dbscan::dbscan(pts, eps = input$eps, minPts = input$minpts)
            cid <- cl$cluster
            # offset to keep unique IDs across groups
            cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
            df$cluster[idx] <- cid
            if (any(cid > 0)) next_id <- max(next_id, max(cid, na.rm = TRUE) + 1L)
          }
        }
      }
      
      # Preserve-TID mode: reassign all points of each TID to the cluster that
      # holds the plurality of that TID's non-noise DBSCAN points. TIDs where
      # every point is noise remain noise. This prevents TID splitting and merges
      # TIDs that overlap spatially into the same cluster.
      if (TRUE) {
        tid_plurality <- df %>%
          dplyr::filter(cluster != 0L) %>%
          dplyr::group_by(tid, cluster) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(tid) %>%
          dplyr::slice_max(n, n = 1L, with_ties = FALSE) %>%
          dplyr::ungroup() %>%
          dplyr::select(tid, cluster_new = cluster)
        df <- df %>%
          dplyr::left_join(tid_plurality, by = "tid") %>%
          dplyr::mutate(cluster = dplyr::coalesce(cluster_new, 0L)) %>%
          dplyr::select(-cluster_new)
      }

      pal <- rep(base_colors, length.out = max(1, max(df$cluster, na.rm = TRUE) + 1))
      df$color <- ifelse(df$cluster == 0, "#bdbdbd",
                         pal[(df$cluster %% length(pal)) + 1])

      if (isTRUE(input$exclude_unclustered)) {
        df <- df %>% dplyr::filter(cluster != 0)
      }
    } else {
      # color by TID
      extended_palette_tid <- rep(base_colors, length.out = nrow(df))
      df <- df %>%
        dplyr::mutate(
          color = extended_palette_tid[(as.numeric(as.factor(tid)) %% length(extended_palette_tid)) + 1]
        )
      df$cluster <- NA_integer_
    }
    
    # Secondary clustering on medians  — group-aware
    if (isTRUE(input$enable_cluster2)) {
      scope <- if (is.null(input$cluster_group_scope)) "all" else input$cluster_group_scope
      mode  <- if (is.null(input$grouping_mode)) "none" else input$grouping_mode
      
      use_primary <- isTRUE(input$enable_clustering) && "cluster" %in% names(df)
      
      if (use_primary) {
        # per-primary-cluster medians
        medians <- df %>%
          dplyr::filter(!is.na(.data$cluster) & .data$cluster != 0) %>%
          dplyr::group_by(.data$cluster) %>%
          dplyr::summarise(
            id = dplyr::first(.data$cluster),
            mx = stats::median(.data$x, na.rm = TRUE),
            my = stats::median(.data$y, na.rm = TRUE),
            mz = stats::median(.data$z, na.rm = TRUE),
            analysis_group = {
              v <- .data$analysis_group[!is.na(.data$analysis_group)]
              if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
            },
            mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(.data$mplex_channel) else NA_integer_,
            .groups = "drop"
          ) %>%
          dplyr::mutate(id_type = "cluster")
      } else {
        # per-trace medians
        medians <- df %>%
          dplyr::group_by(.data$tid) %>%
          dplyr::summarise(
            id = dplyr::first(.data$tid),
            mx = stats::median(.data$x, na.rm = TRUE),
            my = stats::median(.data$y, na.rm = TRUE),
            mz = stats::median(.data$z, na.rm = TRUE),
            analysis_group = dplyr::first(.data$analysis_group),
            mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(.data$mplex_channel) else NA_integer_,
            .groups = "drop"
          ) %>%
          dplyr::mutate(id_type = "tid")
      }
      
      medians$cluster2 <- 0L
      
      # helper to run dbscan on a subset (3D) and keep IDs unique across groups
      run_dbscan_subset <- function(idx, next_id) {
        if (length(idx) < 1) return(list(cid = integer(0), next_id = next_id))
        pts <- as.matrix(medians[idx, c("mx","my","mz")])
        cid <- run_cluster2_grouping(
          pts,
          method        = input$cluster2_method %||% "dbscan",
          eps           = input$eps2,
          minPts        = input$minpts2,
          xi            = input$cluster2_xi %||% 0.05,
          prefilter     = isTRUE(input$cluster2_prefilter),
          prefilter_r   = input$cluster2_prefilter_r,
          prefilter_min = input$cluster2_prefilter_min,
          prune         = isTRUE(input$cluster2_prune),
          prune_tol     = input$cluster2_prune_tol %||% 2
        )
        if (any(cid > 0)) {
          cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
          next_id <- max(next_id, max(cid, na.rm = TRUE) + 1L)
        }
        list(cid = cid, next_id = next_id)
      }
      
      mplex_on2 <- "mplex_channel" %in% names(medians) && any(!is.na(medians$mplex_channel))
      if (mplex_on2) {
        # Multiplex: secondary clusters confined to a single channel.
        next_id <- 1L
        for (chv in sort(unique(medians$mplex_channel[!is.na(medians$mplex_channel)]))) {
          idx <- which(medians$mplex_channel == chv)
          res <- run_dbscan_subset(idx, next_id)
          medians$cluster2[idx] <- res$cid
          next_id <- res$next_id
        }
      } else if (scope == "all" || mode == "none") {
        res <- run_dbscan_subset(seq_len(nrow(medians)), next_id = 1L)
        medians$cluster2 <- res$cid
      } else {
        groups_to_run <- switch(scope,
                                g1 = 1L, g2 = 2L, split = c(1L, 2L), all = integer(0)
        )
        next_id <- 1L
        for (g in groups_to_run) {
          idx <- which(medians$analysis_group == g)
          res <- run_dbscan_subset(idx, next_id)
          medians$cluster2[idx] <- res$cid
          next_id <- res$next_id
        }
      }
      
      # Attach cluster2 back to raw points (join key depends on id_type)
      if (use_primary) {
        df <- df %>% dplyr::left_join(medians[, c("id","cluster2")], by = c("cluster" = "id"))
      } else {
        df <- df %>% dplyr::left_join(medians[, c("id","cluster2")], by = c("tid" = "id"))
      }
      
      # count members per cluster2 (number of unique primary units)
      id_col <- if (use_primary) "cluster" else "tid"
      c2_counts <- df %>%
        dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
        dplyr::group_by(cluster2) %>%
        dplyr::summarise(n_in_cluster2 = dplyr::n_distinct(.data[[id_col]]), .groups = "drop")
      df <- df %>% dplyr::left_join(c2_counts, by = "cluster2")
      df$n_in_cluster2 <- as.integer(df$n_in_cluster2)

      # consistent colors for cluster2
      max_c2 <- suppressWarnings(max(df$cluster2, na.rm = TRUE))
      if (is.finite(max_c2) && max_c2 > 0) {
        pal2 <- rep(base_colors, length.out = max_c2 + 1)
        df$cluster2_color <- ifelse(is.na(df$cluster2) | df$cluster2 == 0, NA_character_,
                                    pal2[(df$cluster2 %% length(pal2)) + 1])
      } else {
        df$cluster2_color <- NA_character_
      }

      df$curvature <- NA_real_
    } else {
      df$cluster2 <- NA_integer_
      df$cluster2_color <- NA_character_
      df$n_in_cluster2 <- NA_integer_
      df$curvature <- NA_real_
    }

    
    # 5) Two-color grouping + optional 2-color clustering
    if (isTRUE(input$enable_two_color)) {
      
      # Make sure columns exist
      if (!"z" %in% names(df)) df$z <- 0
      df$cluster_tc1 <- NA_integer_
      df$cluster_tc2 <- NA_integer_
      
      # --- two-color group map (DCR default; FLUO override if enabled+present) ---
      map <- make_two_color_tid_map(df, input)
      
      # Preserve your existing analysis_group from assign_analysis_group()
      # (optional but safer, since you currently drop it in the join)
      df$analysis_group_main <- df$analysis_group
      
      df <- df %>%
        dplyr::select(-dplyr::any_of(c("two_color_group", "analysis_group", "dcr_med", "fluo_med"))) %>%
        dplyr::left_join(map, by = "tid")
      
      df$two_color_group <- df$analysis_group   # this is the two-color group (1/2)
      df$analysis_group  <- df$analysis_group_main
      df$analysis_group_main <- NULL
      
      # --- 2-color primary clustering (raw points within each group, TIDs kept atomic) ---
      if (isTRUE(input$enable_two_color_cluster1)) {

        df$cluster_tc1 <- 0L

        for (g in c(1L, 2L)) {
          idx <- which(df$two_color_group == g & is.finite(df$x) & is.finite(df$y) & is.finite(df$z))
          if (length(idx) < (input$minpts_tc1 %||% 5L)) next

          pts <- if (isTRUE(input$mode3d == "3D")) {
            as.matrix(df[idx, c("x","y","z")])
          } else {
            as.matrix(df[idx, c("x","y")])
          }

          cl  <- dbscan::dbscan(pts, eps = input$eps_tc1 %||% 20, minPts = input$minpts_tc1 %||% 5L)
          cid <- cl$cluster
          if (g == 2L) cid[cid > 0] <- cid[cid > 0] + 10000L
          df$cluster_tc1[idx] <- cid
        }

        # Lock TIDs: every localization of a TID moves to the cluster its majority belongs to
        tid_plurality <- df %>%
          dplyr::filter(cluster_tc1 != 0L) %>%
          dplyr::group_by(tid, cluster_tc1) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(tid) %>%
          dplyr::slice_max(n, n = 1L, with_ties = FALSE) %>%
          dplyr::ungroup() %>%
          dplyr::select(tid, cluster_tc1_new = cluster_tc1)
        df <- df %>%
          dplyr::left_join(tid_plurality, by = "tid") %>%
          dplyr::mutate(cluster_tc1 = dplyr::coalesce(cluster_tc1_new, 0L)) %>%
          dplyr::select(-cluster_tc1_new)

        if (isTRUE(input$exclude_unclustered_tc1)) {
          df <- df %>% dplyr::filter(!is.na(cluster_tc1) & cluster_tc1 != 0)
        }
      }
      
      # --- 2-color secondary clustering (tc1 cluster medians if tc1 on, else TID medians) ---
      if (isTRUE(input$enable_two_color_cluster2)) {

        if (isTRUE(input$enable_two_color_cluster1) && "cluster_tc1" %in% names(df)) {
          med <- df %>%
            dplyr::filter(!is.na(cluster_tc1) & cluster_tc1 != 0L) %>%
            dplyr::group_by(cluster_tc1) %>%
            dplyr::summarise(
              mx = stats::median(x, na.rm = TRUE),
              my = stats::median(y, na.rm = TRUE),
              mz = stats::median(z, na.rm = TRUE),
              two_color_group = {
                v <- two_color_group[!is.na(two_color_group)]
                if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
              },
              .groups = "drop"
            ) %>%
            dplyr::filter(two_color_group %in% c(1L, 2L), is.finite(mx), is.finite(my), is.finite(mz))
          join_key <- "cluster_tc1"
        } else {
          med <- df %>%
            dplyr::group_by(tid) %>%
            dplyr::summarise(
              mx = stats::median(x, na.rm = TRUE),
              my = stats::median(y, na.rm = TRUE),
              mz = stats::median(z, na.rm = TRUE),
              two_color_group = dplyr::first(two_color_group),
              .groups = "drop"
            ) %>%
            dplyr::filter(two_color_group %in% c(1L, 2L), is.finite(mx), is.finite(my), is.finite(mz))
          join_key <- "tid"
        }

        med$cluster_tc2 <- 0L

        for (g in c(1L, 2L)) {
          idx <- which(med$two_color_group == g)
          if (length(idx) < (input$minpts_tc2 %||% 5L)) next

          pts <- if (isTRUE(input$mode3d == "3D")) {
            as.matrix(med[idx, c("mx","my","mz")])
          } else {
            as.matrix(med[idx, c("mx","my")])
          }

          cid <- run_cluster2_grouping(
            pts,
            method        = input$cluster2_method_tc2 %||% "dbscan",
            eps           = input$eps_tc2 %||% 100,
            minPts        = input$minpts_tc2 %||% 5L,
            xi            = input$cluster2_xi_tc2 %||% 0.05,
            prefilter     = isTRUE(input$cluster2_prefilter_tc2),
            prefilter_r   = input$cluster2_prefilter_r_tc2,
            prefilter_min = input$cluster2_prefilter_min_tc2,
            prune         = isTRUE(input$cluster2_prune_tc2),
            prune_tol     = input$cluster2_prune_tol_tc2 %||% 2
          )
          if (g == 2L) cid[cid > 0] <- cid[cid > 0] + 10000L
          med$cluster_tc2[idx] <- cid
        }

        df <- df %>%
          dplyr::select(-dplyr::any_of("cluster_tc2")) %>%
          dplyr::left_join(
            med %>% dplyr::select(dplyr::all_of(join_key), cluster_tc2),
            by = join_key
          )
      }
      
      # --- apply two-color point colors ---
      col1 <- input$color_group1 %||% "#1f77b4"
      col2 <- input$color_group2 %||% "#d62728"
      
      df$color <- dplyr::case_when(
        df$two_color_group == 1L ~ col1,
        df$two_color_group == 2L ~ col2,
        TRUE ~ df$color
      )
      
    } else {
      df$two_color_group <- NA_integer_
      df$cluster_tc1 <- NA_integer_
      df$cluster_tc2 <- NA_integer_
    }
    
    
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

    tc1_lbl <- format_tc_id(df$cluster_tc1)
    tc2_lbl <- format_tc_id(df$cluster_tc2)
    
    df <- df %>% dplyr::mutate(
      hover_text = paste0(
        if ("tid" %in% names(.)) paste0("Trace ID: ", tid) else "",
        "<br>X: ", round(x, 2),
        "<br>Y: ", round(y, 2),
        if ("z" %in% names(.)) paste0("<br>Z: ", round(z, 2)) else "",
        if (!all(is.na(cluster)))  paste0("<br>Cluster: ", cluster) else "",
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
          paste0("<br>2-color cluster1: ", tc1_lbl),
          ""),
        ifelse(
          isTRUE(input$enable_two_color) & !is.na(tc2_lbl),
          paste0("<br>2-color cluster2: ", tc2_lbl),
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
    content  = function(file) openxlsx::write.xlsx(nn_table(), file, row.names = FALSE),
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
  
  # 1) Cross-DCR NN table ------------- 
  output$download_nn_cross_dcr <- downloadHandler(
    filename = function() paste0("nn_two_color_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      openxlsx::write.xlsx(two_color_nn_edges(), file, row.names = FALSE)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
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
      if (is.null(tbl) || nrow(tbl) == 0) stop("No Cluster2 metrics to download.")
      openxlsx::write.xlsx(
        x         = tbl,
        file      = file,
        sheetName = "Cluster2 metrics",
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
        sheetName = "NN cluster2",
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
    shiny::validate(shiny::need(nrow(df) > 0, "No points to display — please check filters or upload a valid file."))
    
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
    
    
    cluster2_medians <- NULL
    if (isTRUE(input$enable_cluster2) && "cluster2" %in% names(df)) {
      cluster2_medians <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          cluster2 = dplyr::first(cluster2),
          .groups = "drop"
        ) %>%
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
      
      alpha_use <- as.numeric(input$alpha_value_hull %||% 200)
      if (!is.finite(alpha_use) || alpha_use <= 0) alpha_use <- 200
      
      for (gid_key in as.character(group_ids)) {
        
        subset_df <- pts[as.character(pts[[id_col]]) == gid_key, , drop = FALSE]
        if (nrow(subset_df) < 4) next
        
        pts_mat <- as.matrix(subset_df[, c("x","y","z"), drop = FALSE])
        pts_mat <- pts_mat[stats::complete.cases(pts_mat), , drop = FALSE]
        if (nrow(pts_mat) < 4) next
        
        zsd <- stats::sd(pts_mat[, 3], na.rm = TRUE)
        if (!is.finite(zsd) || zsd < 1e-12) {
          pts_mat[, 3] <- pts_mat[, 3] + seq_len(nrow(pts_mat)) * 1e-6
        }
        
        ashape <- tryCatch(
          alphashape3d::ashape3d(pts_mat, alpha = alpha_use),
          error = function(e) NULL
        )
        if (is.null(ashape) || is.null(ashape$triang) || !nrow(ashape$triang)) next
        
        tris <- as.matrix(ashape$triang)
        
        # name-based lookup (NOW guaranteed)
        hull_col <- unname(gid_hex_map[[as.character(gid_key)]])
        
        if (is.null(hull_col) || !nzchar(hull_col)) hull_col <- "#bdbdbd"
        
        fig <- fig %>% plotly::add_trace(
          type = "mesh3d",
          x = pts_mat[, 1], y = pts_mat[, 2], z = pts_mat[, 3],
          i = tris[, 1] - 1, j = tris[, 2] - 1, k = tris[, 3] - 1,
          opacity = input$hull_opacity %||% 0.2,
          
          # Force solid mesh color (Plotly sometimes ignores `color` for mesh3d)
          color       = I(hull_col),
          facecolor   = rep(hull_col, nrow(tris)),        # per-triangle color
          vertexcolor = rep(hull_col, nrow(pts_mat)),     # per-vertex color
          
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
    if (isTRUE(input$show_nn_lines) && !isTRUE(input$draw_nn_two_color) && nrow(df) > 1) {
      
      # Build geometry table: cluster medians if clustering; else trace medians
      if (isTRUE(input$enable_clustering)) {
        df_geo <- df %>%
          dplyr::filter(!is.na(cluster) & cluster != 0) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
            color = dplyr::first(color),
            cluster2 = if ("cluster2" %in% names(df)) dplyr::first(cluster2) else NA_integer_,
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
            mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
            .groups = "drop"
          )
      }

      if (nrow(df_geo) >= 2) {
        coords <- as.matrix(df_geo[, c("x","y","z")])
        pm <- .nn_pairs_scoped(df_geo, coords, input)
        if (nrow(pm)) {
          seg_x <- seg_y <- seg_z <- lab_x <- lab_y <- lab_z <- labs <- c()
          seen  <- new.env(parent = emptyenv())  # dedupe (i,j) == (j,i)

          rng <- input$distance_range %||% c(0, Inf)

          for (r in seq_len(nrow(pm))) {
            i <- pm[r, 1]; ni <- pm[r, 2]
            {
              # dedupe undirected pair
              key <- paste0(min(i, ni), "_", max(i, ni))
              if (!is.null(seen[[key]])) next
              seen[[key]] <- TRUE

              dist <- pm[r, 3]
              if (dist < rng[1] || dist > rng[2]) next

              # within-cluster2 filter
              if (isTRUE(input$nn_within_cluster2) && "cluster2" %in% names(df_geo)) {
                c2i <- df_geo$cluster2[i]; c2ni <- df_geo$cluster2[ni]
                if (is.na(c2i) || is.na(c2ni) || c2i == 0 || c2ni == 0 || c2i != c2ni) next
              }

              seg_x <- c(seg_x, coords[i,1], coords[ni,1], NA_real_)
              seg_y <- c(seg_y, coords[i,2], coords[ni,2], NA_real_)
              seg_z <- c(seg_z, coords[i,3], coords[ni,3], NA_real_)

              nn_off <- suppressWarnings(as.numeric(input$nn_label_offset %||% 5))
              if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
              mid   <- perp_offset_xy(coords[i,], coords[ni,], nn_off)
              lab_x <- c(lab_x, mid[1]); lab_y <- c(lab_y, mid[2]); lab_z <- c(lab_z, mid[3])
              labs  <- c(labs, round(dist, 2))
            }
          }

          if (length(seg_x)) {
            ln_col <- if (isTRUE(input$dark_mode)) "white" else "black"
            fig <- fig %>% plotly::add_trace(
              x = seg_x, y = seg_y, z = seg_z,
              type = "scatter3d", mode = "lines",
              line = list(color = ln_col, width = 2),
              hoverinfo = "none", showlegend = FALSE, inherit = FALSE
            )
            if (length(labs) && isTRUE(input$show_nn_labels)) {
              lab_col <- if (isTRUE(input$dark_mode)) "orange" else "red"
              nn_lfs <- suppressWarnings(as.numeric(input$nn_label_font_size %||% 12))
              if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
              fig <- fig %>% plotly::add_trace(
                x = lab_x, y = lab_y, z = lab_z,
                type = "scatter3d", mode = "text",
                text = labs, textfont = list(color = lab_col, size = nn_lfs),
                hoverinfo = "none", showlegend = FALSE, inherit = FALSE
              )
            }
          }
        }
      }
    }
    if (isTRUE(input$draw_nn_two_color) && isTRUE(input$enable_two_color)) {
      edges <- two_color_nn_edges()
      if (nrow(edges)) {
        col_g1 <- input$color_group1 %||% "#1f77b4"
        col_g2 <- input$color_group2 %||% "#d62728"
        edge_col <- ifelse(edges$g1 == edges$g2,
                           ifelse(edges$g1 == 1L, col_g1, col_g2),
                           "magenta")
        for (cc in unique(edge_col)) {
          sel <- which(edge_col == cc)
          fig <- fig %>% plotly::add_trace(
            x = as.numeric(t(cbind(edges$x1[sel], edges$x2[sel], NA))),
            y = as.numeric(t(cbind(edges$y1[sel], edges$y2[sel], NA))),
            z = as.numeric(t(cbind(edges$z1[sel], edges$z2[sel], NA))),
            type = "scatter3d", mode = "lines",
            line = list(color = cc, width = 3),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
        if (isTRUE(input$show_nn_labels_two_color)) {
          nn_off <- suppressWarnings(as.numeric(input$tc_label_offset %||% 5))
          if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
          nn_lfs <- suppressWarnings(as.numeric(input$tc_label_font_size %||% 12))
          if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
          offsets <- t(mapply(function(x1, y1, z1, x2, y2, z2)
            perp_offset_xy(c(x1, y1, z1), c(x2, y2, z2), nn_off),
            edges$x1, edges$y1, edges$z1, edges$x2, edges$y2, edges$z2))
          mids <- data.frame(
            mx = offsets[, 1],
            my = offsets[, 2],
            mz = offsets[, 3],
            lab = sprintf("%.1f", edges$dist)
          )
          fig <- fig %>% plotly::add_trace(
            x = mids$mx, y = mids$my, z = mids$mz,
            type = "scatter3d", mode = "text",
            text = mids$lab,
            textfont = list(color = if (isTRUE(input$dark_mode)) "white" else "black", size = nn_lfs),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      }
    }
    # ── 2-color alpha shapes ────────────────────────────────────────────────
    if (isTRUE(input$enable_two_color) && (isTRUE(input$show_tc1_alpha) || isTRUE(input$show_tc2_alpha)) &&
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
      
      if (isTRUE(input$show_tc1_alpha) && "cluster_tc1" %in% names(df)) {
        fig <- add_alpha_for(
          cluster_col = "cluster_tc1",
          alpha   = input$tc1_alpha %||% 120,
          opacity = input$tc1_alpha_opacity %||% 0.15,
          min_pts = input$tc1_alpha_min_pts %||% 30
        )
      }
      
      if (isTRUE(input$show_tc2_alpha) && "cluster_tc2" %in% names(df)) {
        fig <- add_alpha_for(
          cluster_col = "cluster_tc2",
          alpha   = input$tc2_alpha %||% 150,
          opacity = input$tc2_alpha_opacity %||% 0.18,
          min_pts = input$tc2_alpha_min_pts %||% 15
        )
      }
    }
    # ── 2-color tc2 isosurface ──────────────────────────────────────────────
    if (isTRUE(input$show_tc2_isosurface) &&
        isTRUE(input$enable_two_color) &&
        "cluster_tc2" %in% names(df)) {
      
      grid_res <- input$tc2_iso_grid_res %||% 22
      pad      <- input$tc2_iso_pad %||% 100
      qlvl     <- input$tc2_iso_level_q %||% 0.85
      op       <- input$tc2_iso_opacity %||% 0.25
      eps0     <- input$tc2_iso_eps %||% 1e-4
      
      # choose which clusters to draw
      ids <- sort(unique(df$cluster_tc2))
      ids <- ids[is.finite(ids) & ids != 0]
      
      if (isTRUE(input$tc2_iso_selected_only)) {
        sel <- parse_tc_label(input$tc2_iso_cluster_label)
        ids <- ids[ids == sel]
      }
      for (cid in ids) {
        sub <- df[df$cluster_tc2 == cid, , drop = FALSE]
        if (nrow(sub) < 20) next
        if (!all(c("x","y","z") %in% names(sub))) next
        
        pts <- as.matrix(sub[, c("x","y","z")])
        pts <- pts[stats::complete.cases(pts), , drop = FALSE]
        if (nrow(pts) < 20) next
        
        # build grid bounding box + padding
        xr <- range(pts[,1]); yr <- range(pts[,2]); zr <- range(pts[,3])
        xs <- seq(xr[1]-pad, xr[2]+pad, length.out = grid_res)
        ys <- seq(yr[1]-pad, yr[2]+pad, length.out = grid_res)
        zs <- seq(zr[1]-pad, zr[2]+pad, length.out = grid_res)
        
        grid <- as.matrix(expand.grid(x = xs, y = ys, z = zs))
        # nearest distance from grid points to point cloud
        nn <- RANN::nn2(data = pts, query = grid, k = 1)
        d  <- as.numeric(nn$nn.dists[,1])
        
        sigma <- input$tc2_iso_sigma %||% 150
        sigma <- as.numeric(sigma)
        if (!is.finite(sigma) || sigma <= 0) sigma <- 150
        
        val <- exp(-(d*d) / (2 * sigma*sigma))
        iso_level <- stats::quantile(val, probs = qlvl, na.rm = TRUE, names = FALSE)
        
        # group A = ids < 10000, group B = ids >= 10000
        tc2_iso_col <- if (cid < 10000L) {
          input$color_group1 %||% "#1f77b4"
        } else {
          input$color_group2 %||% "#d62728"
        }

        fig <- fig %>% plotly::add_trace(
          type = "isosurface",
          x = grid[,1], y = grid[,2], z = grid[,3],
          value = val,
          isomin = iso_level,
          isomax = iso_level * 1.001,
          surface = list(count = 1),
          opacity = op,
          colorscale = list(list(0, tc2_iso_col), list(1, tc2_iso_col)),
          caps = list(x = list(show = FALSE),
                      y = list(show = FALSE),
                      z = list(show = FALSE)),
          showscale = FALSE,
          hoverinfo = "skip",
          inherit = FALSE
        )
      }
    }
    # ── cluster2 isosurface (3D) ────────────────────────────────────────────
    if (isTRUE(input$show_cluster2_isosurface) &&
        isTRUE(input$enable_cluster2) &&
        !is.null(cluster2_medians) &&
        nrow(cluster2_medians) > 0 &&
        isTRUE(input$mode3d == "3D")) {
      
      grid_res <- input$cluster2_iso_grid_res %||% 22
      pad      <- input$cluster2_iso_pad %||% 100
      qlvl     <- input$cluster2_iso_level_q %||% 0.85
      op       <- input$cluster2_iso_opacity %||% 0.25
      sigma    <- input$cluster2_iso_sigma %||% 150
      sigma    <- as.numeric(sigma)
      if (!is.finite(sigma) || sigma <= 0) sigma <- 150
      
      ids <- sort(unique(cluster2_medians$cluster2))
      ids <- ids[is.finite(ids) & ids != 0]

      if (isTRUE(input$cluster2_iso_selected_only)) {
        sel <- as.integer(input$cluster2_iso_id %||% 1L)
        ids <- ids[ids == sel]
      }

      if (length(ids)) {

      # build per-cluster color map using the same palette as cluster2 points
      qual_pals_iso <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
      base_colors_iso <- unlist(mapply(RColorBrewer::brewer.pal, qual_pals_iso$maxcolors, rownames(qual_pals_iso)))
      pal_iso <- rep(base_colors_iso, length.out = max(ids) + 1)

      for (cid in ids) {
        sub <- cluster2_medians[cluster2_medians$cluster2 == cid, , drop = FALSE]
        if (nrow(sub) < 4) next

        pts <- as.matrix(sub[, c("x","y","z")])
        pts <- pts[stats::complete.cases(pts), , drop = FALSE]
        if (nrow(pts) < 4) next

        xr <- range(pts[,1]); yr <- range(pts[,2]); zr <- range(pts[,3])
        xs <- seq(xr[1]-pad, xr[2]+pad, length.out = grid_res)
        ys <- seq(yr[1]-pad, yr[2]+pad, length.out = grid_res)
        zs <- seq(zr[1]-pad, zr[2]+pad, length.out = grid_res)

        grid <- as.matrix(expand.grid(x = xs, y = ys, z = zs))
        nn   <- RANN::nn2(data = pts, query = grid, k = 1)
        d    <- as.numeric(nn$nn.dists[,1])

        val <- exp(-(d*d) / (2 * sigma*sigma))
        iso_level <- stats::quantile(val, probs = qlvl, na.rm = TRUE, names = FALSE)

        iso_col <- if (isTRUE(input$cluster2_iso_unicolor)) {
          input$cluster2_iso_color %||% "#00CED1"
        } else {
          pal_iso[(cid %% length(pal_iso)) + 1]
        }

        fig <- fig %>% plotly::add_trace(
          type = "isosurface",
          x = grid[,1], y = grid[,2], z = grid[,3],
          value = val,
          isomin = iso_level,
          isomax = iso_level * 1.001,
          surface = list(count = 1),
          opacity = op,

          colorscale = list(list(0, iso_col), list(1, iso_col)),

          caps = list(
            x = list(show = FALSE),
            y = list(show = FALSE),
            z = list(show = FALSE)
          ),
          showscale = FALSE,
          hoverinfo = "skip",
          inherit = FALSE
        )
      }
      } # end else(!length(ids))
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
      #    - if primary clustering enabled -> cluster medians
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
        fit_log_rv("No cluster2 structures available for fitting. Check Cluster2 DBSCAN settings (eps2/minpts2).")
      }
      fit_alpha <- as.numeric(input$alpha_value_hull %||% 200)
      if (!is.finite(fit_alpha) || fit_alpha <= 0) fit_alpha <- 200
      fit_alpha <- min(300, max(0.1, fit_alpha))
      
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
      surface_triangles <- function(ashape_obj, pts) {
        tri <- NULL
        if (!is.null(ashape_obj) && !is.null(ashape_obj$triang) && nrow(ashape_obj$triang) > 0) {
          tt <- as.matrix(ashape_obj$triang)
          if (ncol(tt) >= 4) {
            tt <- tt[tt[, 4] >= 2L, 1:3, drop = FALSE]
          } else if (ncol(tt) >= 3) {
            tt <- tt[, 1:3, drop = FALSE]
          } else {
            tt <- NULL
          }
          if (!is.null(tt) && nrow(tt) > 0) tri <- tt
        }
        if (is.null(tri) || !nrow(tri)) {
          tri <- tryCatch({
            hh <- geometry::convhulln(pts, options = "Qt")
            as.matrix(hh[, 1:3, drop = FALSE])
          }, error = function(e) NULL)
        }
        tri
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
              alpha_value = fit_alpha,
              use_cluster1 = isTRUE(input$enable_clustering)
            )
            if (!is.null(fit)) {
              fit$radius <- cap_radius_for_cluster(fit$radius, X)
              if (!is.finite(fit$radius) || fit$radius <= 0) next
              pts <- fit$virtual_points
              ashape <- tryCatch(alphashape3d::ashape3d(pts, alpha = fit_alpha), error = function(e) NULL)
              surf_tris <- surface_triangles(ashape, pts)
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
              alpha_value = fit_alpha,
              use_cluster1 = isTRUE(input$enable_clustering)
            )
            if (!is.null(fit)) {
              fit$axes <- cap_axes_for_cluster(fit$axes, X)
              if (!all(is.finite(fit$axes)) || any(fit$axes <= 0)) next
              pts <- fit$virtual_points
              ashape <- fit$ashape
              surf_tris <- surface_triangles(ashape, pts)
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
              ashape <- tryCatch(alphashape3d::ashape3d(pts, alpha = fit_alpha), error = function(e) NULL)
              surf_tris <- surface_triangles(ashape, pts)
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
              ashape <- tryCatch(alphashape3d::ashape3d(pts, alpha = fit_alpha), error = function(e) NULL)
              surf_tris <- surface_triangles(ashape, pts)
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
          # Use cluster medians (primary DBSCAN)
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
      for (cl2 in unique(df$cluster2[df$cluster2 != 0])) {
        group_df <- df %>% dplyr::filter(cluster2 == cl2)
        if (nrow(group_df) >= 4) {
          pts <- as.matrix(group_df[, c("x", "y", "z")])
          alpha_shape <- tryCatch(alphashape3d::ashape3d(pts, alpha = input$alpha_value_hull), error = function(e) NULL)
          if (!is.null(alpha_shape)) {
            tris <- alpha_shape$triang
            if (nrow(tris) > 0) {
              
              # pick a single solid color: first non-NA cluster2 color (fallback to light gray)
              col_hex <- group_df$cluster2_color[which(!is.na(group_df$cluster2_color) & nzchar(group_df$cluster2_color))][1]
              if (is.na(col_hex) || !nzchar(col_hex)) col_hex <- "#bdbdbd"
              
              fig <- fig %>% plotly::add_trace(
                type = "mesh3d",
                x = pts[, 1], y = pts[, 2], z = pts[, 3],
                i = tris[, 1] - 1, j = tris[, 2] - 1, k = tris[, 3] - 1,
                opacity = input$hull_opacity,
                color = I(col_hex),           
                name = paste("Cluster2", cl2),
                hoverinfo = "skip",
                hovertemplate = "",
                showlegend = FALSE,
                inherit = FALSE
              )
            }
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
      txt3 <- sprintf("GT candidate<br>± %.0f nm%s<br>%d contributing unit(s) · %d TID / %d loc",
                      gpts$pos_sd_nm, ifelse(gpts$fused, " (fused: multiple GTs)", ""),
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

  # file -> group-name map (only the first n_groups slots, only assigned files)
  file_group_map <- reactive({
    ng <- max(1L, min(6L, as.integer(input$n_groups %||% 6L)))
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
    ng <- max(1L, min(6L, as.integer(input$n_groups %||% 6L)))
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
    n <- length(labels)
    cols <- if (n <= 8) RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)]
            else viridisLite::viridis(n)
    stats::setNames(cols, labels)
  }

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

  # ------ Batch NN Distance Histogram ----------------------------------------
  nn_batch_results <- reactiveVal(NULL)

  # Two-color group filter — only visible when two-color is enabled
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
    unit  <- input$nn_batch_unit  %||% "tid"
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
        df <- store[[fname]]
        if (!is.data.frame(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }

        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
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
          # two-color pipeline: assigns analysis_group AND aliases cluster cols
          df <- tryCatch(apply_two_color_clustering_pipeline(df, input), error = function(e) NULL)
          if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }
        } else if (needs_clustering) {
          df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
          if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }
        } else if (needs_grouping) {
          df <- tryCatch(assign_analysis_group(df, input), error = function(e) df)
        }

        d <- tryCatch(
          compute_batch_nn(df, unit = unit, k = k, group_filter = group_filter,
                           cross_groups = cross_groups),
          error = function(e) data.frame()
        )
        if (!nrow(d)) { n_skipped <- n_skipped + 1L; next }
        d$dataset <- fname
        d$tag     <- if (".dataset_tag" %in% names(df)) df$.dataset_tag[[1L]] else fname
        results[[fname]] <- d
      }
    })

    if (!length(results)) {
      hint <- switch(unit,
        pcid          = "Check that Primary point DBSCAN Clustering is enabled and eps/minPts produce clusters.",
        cluster2_tid  = "Check that Median Clustering is enabled.",
        cluster2_pcid = "Check that both clustering levels are enabled.",
        "Check filters — all rows may have been filtered out."
      )
      showNotification(paste0("No distances computed. ", hint), type = "warning")
      nn_batch_results(NULL)
      return()
    }
    nn_batch_results(do.call(rbind, results))
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
    apply_dataset_groups(res)   # relabel dataset -> group in group mode
  })

  output$nn_batch_hist <- renderPlotly({
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
      datasets <- unique(res$dataset)
      n_ds     <- length(datasets)
      grouped_mode <- (input$analysis_by %||% "file") == "group"
      # group mode: distinct group colour per series; file mode: shade one base
      # colour by dataset index (alpha) so bars stay distinct.
      pal      <- if (grouped_mode) series_palette(datasets) else NULL
      alphas   <- if (n_ds == 1L) "cc"
                  else sprintf("%02x", as.integer(seq(204, 80, length.out = n_ds)))

      for (i in seq_along(datasets)) {
        sub      <- res[res$dataset == datasets[i], , drop = FALSE]
        lbl      <- if (grouped_mode) datasets[i]
                    else if ("tag" %in% names(sub)) sub$tag[[1L]] else datasets[i]
        col_base <- if (grouped_mode) unname(pal[datasets[i]]) else group_col
        fill_col <- if (grouped_mode) paste0(col_base, "99") else paste0(group_col, alphas[i])
        fig <- fig %>% plotly::add_histogram(
          x        = sub$distance,
          name     = lbl,
          legendgroup = datasets[i],
          marker   = list(color = fill_col),
          xbins    = xbins,
          histnorm = if (density) "probability density" else "",
          opacity  = 0.75
        )
        fig <- add_spread(fig, sub$distance, col_base, datasets[i])
      }
    }

    indicator_note <- switch(error_type,
      iqr  = "— median  ··· Q25 / Q75",
      mad  = "— median  ··· ±MAD",
      "— median"
    )
    unit_lbl <- switch(input$nn_batch_unit %||% "tid",
      tid           = "TID centroid",
      pcid          = "Primary cluster centroid",
      cluster2_tid  = "TID-within-Cluster2",
      cluster2_pcid = "Primary cluster centroid-within-Cluster2",
      input$nn_batch_unit %||% "tid"
    )
    fig %>% plotly::layout(
      barmode = "overlay",
      xaxis   = list(title = sprintf("k=%d NN distance (nm)  [%s]",
                                     as.integer(input$nn_batch_k %||% 1L), unit_lbl)),
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
  })

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

  observeEvent(input$run_tid_count, {
    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }

    unit   <- input$tid_count_unit %||% "cluster2"
    cl_col <- if (unit == "pcid") "cluster" else "cluster2"

    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }
    if (isTRUE(input$tid_count_by_channel)) {
      cs <- mplex_split_active()
      if (length(cs)) store <- cs
    }

    all_counts <- list()
    withProgress(message = "Counting TIDs per cluster...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- store[[fname]]
        if (!is.data.frame(df) || !nrow(df)) next
        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) next
        df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) next
        if (!cl_col %in% names(df) || all(is.na(df[[cl_col]]) | df[[cl_col]] == 0L)) next
        cnt <- df %>%
          dplyr::filter(!is.na(.data[[cl_col]]) & .data[[cl_col]] > 0L) %>%
          dplyr::group_by(.data[[cl_col]]) %>%
          dplyr::summarise(n_tids = dplyr::n_distinct(tid), .groups = "drop") %>%
          dplyr::mutate(dataset = fname)
        all_counts[[fname]] <- cnt
      }
    })
    counts <- do.call(rbind, Filter(Negate(is.null), all_counts))

    if (is.null(counts) || !nrow(counts)) {
      hint <- if (cl_col == "cluster")
        "Enable Primary Clustering and press Plot for each file."
      else
        "Enable Median Clustering and press Plot for each file."
      showNotification(paste0("No clusters found. ", hint), type = "warning")
      tid_count_results(NULL)
      return()
    }
    tid_count_results(counts)
    showNotification(
      sprintf("Counted TIDs in %d clusters across %d file(s) (median = %g).",
              nrow(counts), length(unique(counts$dataset)),
              stats::median(counts$n_tids)),
      type = "message", duration = 4)
  })

  output$tid_count_hist <- renderPlotly({
    counts <- tid_count_results()
    shiny::validate(shiny::need(!is.null(counts) && nrow(counts) > 0,
                                "Press 'Run' to count traces per cluster."))
    counts <- apply_dataset_groups(counts)
    shiny::validate(shiny::need(nrow(counts) > 0, "No files assigned to groups (group mode)."))

    bins      <- max(5L, as.integer(input$tid_count_bins %||% 30L))
    density   <- identical(input$tid_count_yaxis %||% "count", "density")
    unit      <- input$tid_count_unit %||% "cluster2"
    lbl       <- if (unit == "pcid") "Primary point DBSCAN Cluster" else "Median Cluster"
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
        xaxis = list(title = paste0("TIDs per ", lbl)),
        yaxis = list(title = if (density) "Density" else "Number of clusters"),
        legend = list(orientation = "h", y = -0.18, groupclick = "togglegroup"),
        paper_bgcolor = "white",
        plot_bgcolor  = "white"
      )
  })

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
    if (!isTRUE(input$enable_cluster2)) {
      showNotification("Enable Median Clustering in the Display tab first.", type = "warning"); return()
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
        df <- store[[fname]]
        if (!is.data.frame(df) || !nrow(df)) next
        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df) || !"tim" %in% names(df)) next
        df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) next
        per_tid <- df %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(
            tim_med = stats::median(tim, na.rm = TRUE),
            in_c2   = any(!is.na(cluster2) & cluster2 > 0L),
            .groups = "drop"
          ) %>%
          dplyr::filter(is.finite(tim_med)) %>%
          dplyr::mutate(dataset = fname)
        if (nrow(per_tid)) all_tid[[fname]] <- per_tid
      }
    })
    res <- do.call(rbind, Filter(Negate(is.null), all_tid))
    if (is.null(res) || !nrow(res)) {
      showNotification("No traces with a tim column and Median cluster found. Enable Median Clustering and press Plot for each file.",
                       type = "warning")
      det_time_results(NULL); return()
    }
    det_time_results(res)
    showNotification(
      sprintf("Detection over time: %d traces across %d file(s); %.0f%% in Cluster2.",
              nrow(res), length(unique(res$dataset)), 100 * mean(res$in_c2)),
      type = "message", duration = 4)
  })

  output$det_time_plot <- renderPlotly({
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
                          line = list(color = "#2C7FB8"), marker = list(color = "#2C7FB8"), name = "In Cluster2") %>%
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
                            name = paste0(g, " — In Cluster2", kslab)) %>%
          plotly::add_trace(x = ctr, y = yy$bg, type = "scatter", mode = "lines",
                            line = list(color = pal[g], width = 2, dash = "dash"),
                            name = paste0(g, " — Unclustered"))
      }
    }
    fig %>% plotly::layout(
      title  = ttl,
      xaxis  = list(title = "Acquisition time (tim)"),
      yaxis  = list(title = ytitle),
      legend = list(orientation = "h")
    )
  })

  output$download_det_time <- downloadHandler(
    filename = function()
      paste0("cluster2_detection_over_time_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- det_time_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  ## ── Cluster2 Crowding vs Detection-Time Shift ─────────────────────────────
  # One row per Cluster2 structure (built from its member TID centroids):
  #   nn_median  = median nearest-neighbour distance among members (smaller = denser)
  #   time_shift = median(member tim_med) - median(unclustered background tim_med),
  #                computed per file (positive = structure detected later)
  # Tests whether denser (more crowded) structures are detected later.
  crowd_results <- reactiveVal(NULL)

  observeEvent(input$run_crowd, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    if (!isTRUE(input$enable_cluster2)) {
      showNotification("Enable Median Clustering in the Display tab first.", type = "warning"); return()
    }
    if (!isTRUE(input$crowd_all_files)) {
      af <- input$active_file
      if (!is.null(af) && af %in% names(store)) store <- store[af]
    }

    min_n  <- max(3L, as.integer(input$crowd_min_members %||% 4L))
    use_3d <- isTRUE(input$crowd_use_3d)
    cols   <- if (use_3d) c("x", "y", "z") else c("x", "y")

    rows <- list()
    withProgress(message = "Computing crowding vs detection-time shift...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- store[[fname]]
        if (!is.data.frame(df) || !nrow(df) || !"tim" %in% names(df)) next
        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) next
        df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) next
        if (use_3d && !"z" %in% names(df)) next

        per_tid <- df %>%
          dplyr::group_by(tid) %>%
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

        bg_med <- stats::median(per_tid$tim_med[per_tid$cl2 == 0L], na.rm = TRUE)
        if (!is.finite(bg_med)) next     # need an unclustered background to reference

        for (cid in unique(per_tid$cl2[per_tid$cl2 > 0L])) {
          sub <- per_tid[per_tid$cl2 == cid, , drop = FALSE]
          n   <- nrow(sub)
          if (n < min_n) next
          pts <- as.matrix(sub[, cols, drop = FALSE])
          pts <- pts[stats::complete.cases(pts), , drop = FALSE]
          if (nrow(pts) < 2L) next
          nnd <- tryCatch(dbscan::kNN(pts, k = 1)$dist, error = function(e) NULL)
          if (is.null(nnd)) next

          rows[[length(rows) + 1L]] <- data.frame(
            dataset    = fname,
            cluster2   = cid,
            n_tids     = n,
            nn_median  = stats::median(nnd),
            time_shift = stats::median(sub$tim_med) - bg_med,
            stringsAsFactors = FALSE
          )
        }
      }
    })

    res <- if (length(rows)) do.call(rbind, rows) else NULL
    if (is.null(res) || !nrow(res)) {
      showNotification(
        sprintf("No Cluster2 structures with >= %d traces and an unclustered background found.", min_n),
        type = "warning")
      crowd_results(NULL); return()
    }
    crowd_results(res)
    showNotification(
      sprintf("Crowding vs detection-time shift: %d structures across %d file(s).",
              nrow(res), length(unique(res$dataset))),
      type = "message", duration = 4)
  })

  output$crowd_plot <- renderPlotly({
    res <- crowd_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Press 'Run' to compute crowding vs detection-time shift."))
    res <- apply_dataset_groups(res)
    shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
    fig <- plotly::plot_ly()
    datasets <- unique(res$dataset)
    n_ds     <- length(datasets)
    pal <- series_palette(datasets)

    for (ds in datasets) {
      sub <- res[res$dataset == ds, , drop = FALSE]
      fig <- fig %>% plotly::add_markers(
        x = sub$nn_median, y = sub$time_shift,
        name = if (n_ds > 1) ds else "Structures",
        marker = list(color = pal[ds], size = 7, opacity = 0.7,
                      line = list(color = "white", width = 0.5)),
        text = paste0("cluster2 ", sub$cluster2, "<br>n = ", sub$n_tids,
                      "<br>NN median = ", round(sub$nn_median, 1), " nm",
                      "<br>time shift = ", round(sub$time_shift, 2)),
        hoverinfo = "text"
      )
    }

    # optional loess trend — per group when grouped, else over the pooled cloud
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

    fig %>% plotly::layout(
      margin = list(l = 70),
      xaxis  = list(title = "Median NN distance among members (nm)  -  smaller = more crowded"),
      yaxis  = list(title = "Detection-time shift (tim)",
                    zeroline = TRUE, zerolinecolor = "#888888", zerolinewidth = 1),
      legend = list(orientation = "h")
    )
  })

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
  # One row per Cluster2: size = # member units (trace IDs, or primary clusters
  # if primary clustering is on); nn_median = median 1-NN distance among the
  # member centroids. Tests whether spacing is constant with size or "balloons".
  size_nn_results <- reactiveVal(NULL)

  observeEvent(input$run_size_nn, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
    if (!isTRUE(input$enable_cluster2)) {
      showNotification("Enable Median Clustering in the Display tab first.", type = "warning"); return()
    }
    if (!isTRUE(input$size_nn_all_files)) {
      af <- input$active_file
      if (!is.null(af) && af %in% names(store)) store <- store[af]
    }
    min_n       <- max(2L, as.integer(input$size_nn_min_members %||% 3L))
    kk          <- max(1L, as.integer(input$size_nn_k %||% 1L))
    use_3d      <- isTRUE(input$size_nn_use_3d)
    cols        <- if (use_3d) c("x", "y", "z") else c("x", "y")
    use_primary <- isTRUE(input$enable_clustering)

    rows <- list()
    withProgress(message = "Computing within-cluster NN vs size...", value = 0, {
      nf <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]; incProgress(1 / nf, detail = fname)
        df <- store[[fname]]
        if (!is.data.frame(df) || !nrow(df)) next
        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) next
        df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) next
        if (use_3d && !"z" %in% names(df)) next

        id_col <- if (use_primary && "cluster" %in% names(df) &&
                      any(!is.na(df$cluster) & df$cluster > 0L)) "cluster" else "tid"
        cent <- df %>%
          dplyr::filter(!is.na(cluster2) & cluster2 > 0L) %>%
          dplyr::group_by(cluster2, .data[[id_col]]) %>%
          dplyr::summarise(
            x = stats::median(x, na.rm = TRUE),
            y = stats::median(y, na.rm = TRUE),
            z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
            .groups = "drop")
        if (!nrow(cent)) next

        for (cid in unique(cent$cluster2)) {
          sub <- cent[cent$cluster2 == cid, , drop = FALSE]
          n   <- nrow(sub)
          if (n < min_n) next
          pts <- as.matrix(sub[, cols, drop = FALSE])
          pts <- pts[stats::complete.cases(pts), , drop = FALSE]
          if (nrow(pts) < 2L) next
          # mean distance to the k nearest members (k clamped to cluster size - 1),
          # per member, then median over members. k = 1 -> nearest-neighbour spacing.
          k_eff <- min(kk, nrow(pts) - 1L)
          nn <- tryCatch(dbscan::kNN(pts, k = k_eff)$dist, error = function(e) NULL)
          if (is.null(nn)) next
          spacing <- if (k_eff == 1L) as.numeric(nn) else rowMeans(nn, na.rm = TRUE)
          # convex-hull size: area (2D) / volume (3D); NA if degenerate
          hull <- tryCatch(geometry::convhulln(pts, options = "FA")$vol, error = function(e) NA_real_)
          rows[[length(rows) + 1L]] <- data.frame(
            dataset   = fname, cluster2 = cid, size = n,
            nn_median = stats::median(spacing),
            hull_size = if (is.finite(hull) && hull > 0) hull else NA_real_,
            dim       = length(cols),
            unit      = if (id_col == "cluster") "primary_cluster" else "tid",
            stringsAsFactors = FALSE)
        }
      }
    })

    res <- if (length(rows)) do.call(rbind, rows) else NULL
    if (is.null(res) || !nrow(res)) {
      showNotification(sprintf("No Cluster2 structures with >= %d members found.", min_n), type = "warning")
      size_nn_results(NULL); return()
    }
    size_nn_results(res)
    rho <- suppressWarnings(stats::cor(res$size, res$nn_median, method = "spearman"))
    showNotification(
      sprintf("Within-cluster NN vs size: %d structures across %d file(s); Spearman rho = %.2f.",
              nrow(res), length(unique(res$dataset)), rho),
      type = "message", duration = 4)
  })

  output$size_nn_plot <- renderPlotly({
    res <- size_nn_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Press 'Run' to compute within-cluster NN vs size."))
    res <- apply_dataset_groups(res)
    shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
    unit_lbl <- if (all(res$unit == "primary_cluster")) "primary clusters" else "trace IDs"
    xmode    <- input$size_nn_xmode %||% "count"

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
      fig <- fig %>% plotly::add_markers(x = sub$xv, y = sub$nn_median,
        name = if (n_ds > 1) ds else "Structures",
        marker = list(color = pal[ds], size = 7, opacity = 0.7, line = list(color = "white", width = 0.5)),
        text = paste0("cluster2 ", sub$cluster2, "<br>", xname, " ", signif(sub$xv, 4),
                      "<br>median NN ", round(sub$nn_median, 1), " nm<br>", sub$dataset),
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
      ttl <- "Within-cluster NN vs size — per-group fits (slope & rho in legend)"
    } else {
      fig <- add_trends(fig, res, "#555555", "")
      ttl <- sprintf("Spearman rho = %.2f   |   linear slope = %s %s",
                     rho, if (is.finite(sl)) signif(sl, 3) else "NA", slope_u)
    }
    kk   <- max(1L, as.integer(input$size_nn_k %||% 1L))
    ylab <- if (kk == 1L) "Within-cluster median NN distance (nm)"
            else sprintf("Within-cluster median spacing (mean of %d nearest, nm)", kk)
    fig %>% plotly::layout(
      title  = list(text = ttl, font = list(size = 13)),
      xaxis  = list(title = xlab,
                    type = if (isTRUE(input$size_nn_logx)) "log" else "linear"),
      yaxis  = list(title = ylab),
      legend = list(orientation = "h"))
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
  # Inspect one Cluster2 structure at a time. A Run button snapshots the plotted
  # file's clusters; Prev/Next and view options then work on that snapshot, so it
  # does not refresh on filter/plot changes until Run is pressed again.
  # MLE best-fit via mle_fit_*; overlay from the fitted parameters.
  viewer_snapshot <- reactiveVal(NULL)   # list(df, ids) captured on Run
  viewer_idx      <- reactiveVal(1L)

  observeEvent(input$run_viewer, {
    df <- tryCatch(plottable_data(), error = function(e) NULL)
    if (is.null(df) || !"cluster2" %in% names(df)) {
      showNotification("No Median Clusters in the plotted data. Enable Median Clustering and press Plot.",
                       type = "warning")
      viewer_snapshot(NULL); return()
    }
    ids <- sort(unique(df$cluster2[!is.na(df$cluster2) & df$cluster2 > 0L]))
    if (!length(ids)) {
      showNotification("No Cluster2 structures found in the plotted data.", type = "warning")
      viewer_snapshot(NULL); return()
    }
    viewer_snapshot(list(df = df, ids = ids))
    viewer_idx(1L)
    showNotification(sprintf("Loaded %d Cluster2 structure(s).", length(ids)),
                     type = "message", duration = 3)
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
      "Press 'Run' to load the plotted file's Cluster2 structures."))
    ids <- s$ids; n <- length(ids)
    idx <- min(max(viewer_idx(), 1L), n)
    id  <- ids[idx]
    df  <- s$df
    sub <- df[!is.na(df$cluster2) & df$cluster2 == id, , drop = FALSE]
    cent <- sub %>% dplyr::group_by(tid) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = if ("z" %in% names(sub)) stats::median(z, na.rm = TRUE) else 0,
        .groups = "drop")
    zc    <- if ("z" %in% names(sub)) sub$z else rep(0, nrow(sub))
    cloud <- cbind(x = sub$x, y = sub$y, z = zc)
    cloud <- cloud[stats::complete.cases(cloud), , drop = FALSE]
    list(idx = idx, n = n, id = id, sub = sub,
         mat = as.matrix(cent[, c("x", "y", "z")]),
         cloud = cloud)
  })

  # Full MLE fit (optim orthogonal-distance) of all 5 shapes to the trace-median
  # centroids — same pathway the Display tab uses when fit_method = "mle", not
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

    fitS <- tryCatch(do.call(mle_fit_sphere,    c(list(X = X, w = NULL), cargs)), error = function(e) NULL)
    fitE <- tryCatch(do.call(mle_fit_ellipsoid, c(list(X = X, w = NULL), cargs)), error = function(e) NULL)
    fitH <- tryCatch(do.call(mle_fit_helix,     c(list(X = X, w = NULL, turns_hint = turns), cargs)), error = function(e) NULL)
    fitR <- tryCatch(do.call(mle_fit_ring,      c(list(X = X, w = NULL), cargs)), error = function(e) NULL)
    fitA <- tryCatch(do.call(mle_fit_arc,       c(list(X = X, w = NULL), cargs)), error = function(e) NULL)

    rmsd_of <- function(d) { d <- d[is.finite(d)]; if (length(d)) sqrt(mean(d^2)) else NA_real_ }
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

    # ── model selection: AIC / AICc / BIC under iid isotropic Gaussian residuals ──
    # sigma^2_hat = RMSD^2, so -2logL = n*log(2*pi*sigma2) + n. k = geometric DOF + 1
    # (the variance). AICc adds the small-sample correction and -> Inf when a model
    # has too few points to estimate (n <= k+1), which refuses over-parameterised fits.
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
    sprintf("Cluster %d / %d   |   id %d   |   %d traces%s",
            cur$idx, cur$n, cur$id, nrow(cur$mat),
            if (is.finite(rb)) sprintf("   |   best MLE fit (%s): %s (RMSD %.1f nm)", crit, best, rb) else "")
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
      shiny::validate(shiny::need(nrow(mat) >= 3, "Need >= 3 traces for a PCA projection."))
      # PCA basis from trace-median centroids so trace length (localization count)
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
          name = "Trace medians")
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
          marker = list(color = "#08306B", size = 4), name = "Trace medians")
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

  observeEvent(input$run_rog, {
    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }
    col   <- input$rog_cluster_col %||% "cluster"
    store <- smlm_file_store()
    if (!length(store)) {
      showNotification("No files loaded.", type = "warning"); return()
    }
    if (isTRUE(input$rog_by_channel)) {
      cs <- mplex_split_active()
      if (length(cs)) store <- cs
    }

    all_res <- list()
    withProgress(message = "Computing RoG...", value = 0, {
      n_files <- length(store)
      for (i in seq_along(store)) {
        fname <- names(store)[[i]]
        incProgress(1 / n_files, detail = fname)
        df <- store[[fname]]
        if (!is.data.frame(df) || !nrow(df)) next
        df <- tryCatch(apply_meta_filters_batch(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) next
        df <- tryCatch(apply_clustering_pipeline(df, input), error = function(e) NULL)
        if (is.null(df) || !nrow(df) || !col %in% names(df)) next
        r <- tryCatch(compute_cluster_rog(df, cluster_col = col), error = function(e) NULL)
        if (is.null(r) || !nrow(r)) next
        r$dataset <- fname
        all_res[[fname]] <- r
      }
    })
    res <- do.call(rbind, Filter(Negate(is.null), all_res))

    if (is.null(res) || !nrow(res)) {
      hint <- if (col == "cluster")
        "Enable Primary Clustering and press Plot for each file."
      else
        "Enable Median Clustering and press Plot for each file."
      showNotification(paste0("No RoG computed. ", hint), type = "warning")
      rog_results(NULL); return()
    }
    rog_results(res)
    showNotification(
      sprintf("RoG computed for %d cluster(s) across %d file(s).",
              nrow(res), length(unique(res$dataset))),
      type = "message", duration = 3)
  })

  output$rog_hist <- renderPlotly({
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
        xaxis = list(title = "Radius of Gyration (nm)"),
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
  ripley_results <- reactiveVal(NULL)

  observeEvent(input$run_ripley, {
    unit      <- input$ripley_unit %||% "tid"
    use_3d    <- isTRUE(input$ripley_use_3d)
    r_max_in  <- as.numeric(input$ripley_r_max  %||% 0)
    r_steps   <- max(10L, as.integer(input$ripley_r_steps %||% 60L))
    all_files <- isTRUE(input$ripley_all_files)

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

    compute_ripley_for_df <- function(df, tag, r_seq_use) {
      pts <- get_pts_flat(df)
      if (is.null(pts) || nrow(pts) < 2L) return(NULL)
      curve <- compute_ripley_L(pts$x, pts$y, if (use_3d) pts$z else NULL,
                                r_seq = r_seq_use, use_3d = use_3d)
      if (!nrow(curve)) return(NULL)
      curve$dataset <- tag
      curve
    }

    get_df_list <- function() {
      if (all_files) {
        store <- smlm_file_store()
        if (isTRUE(input$ripley_by_channel)) {
          cs <- mplex_split_active()
          if (length(cs)) store <- cs
        }
        shiny::validate(shiny::need(length(store) > 0, "No files loaded."))
        lapply(names(store), function(nm) {
          df <- apply_meta_filters_batch(store[[nm]], input)
          if (needs_cluster) {
            if (isTRUE(input$enable_two_color)) {
              df <- apply_two_color_clustering_pipeline(df, input)
            } else {
              df <- apply_clustering_pipeline(df, input)
            }
          }
          list(df = df, tag = nm)
        })
      } else {
        df <- plottable_data()
        shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "Press Plot first."))
        nm <- input$active_file %||% "Active file"
        list(list(df = df, tag = nm))
      }
    }

    df_list <- get_df_list()

    all_curves <- lapply(df_list, function(item) {
      r_seq_use <- if (r_max_in > 0) seq(0, r_max_in, length.out = r_steps + 1L)[-1] else NULL
      compute_ripley_for_df(item$df, item$tag, r_seq_use)
    })
    all_curves <- do.call(rbind, Filter(Negate(is.null), all_curves))
    shiny::validate(shiny::need(!is.null(all_curves) && nrow(all_curves) > 0, "Ripley computation returned no data."))
    ripley_results(all_curves)

  })

  output$ripley_plot <- renderPlotly({
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
    fig <- plotly::plot_ly()

    if (!grouped) {
      datasets <- unique(res$dataset)
      pal <- series_palette(datasets)
      for (nm in datasets) {
        sub <- res[res$dataset == nm, ]
        fig <- fig %>% plotly::add_trace(
          x = sub$r, y = smooth_vec(sub$r, sub$Lr_minus_r),
          type = "scatter", mode = "lines", name = nm,
          line = list(color = unname(pal[nm]), width = 2))
      }
    } else {
      gm <- file_group_map()
      shiny::validate(shiny::need(length(gm) > 0, "No files assigned to groups (group mode)."))
      res <- res[res$dataset %in% names(gm), , drop = FALSE]
      shiny::validate(shiny::need(nrow(res) > 0, "No group files in the Ripley results."))
      res$group <- unname(gm[res$dataset])
      pal  <- series_palette(unique(res$group))
      band <- input$ripley_group_band %||% "iqr"

      # common r-grid up to the smallest per-file r max (avoid extrapolation),
      # then average each group's per-file curves + between-file spread band.
      file_maxr <- tapply(res$r, res$dataset, function(v) max(v, na.rm = TRUE))
      rmax_c <- suppressWarnings(min(file_maxr, na.rm = TRUE))
      steps  <- max(10L, as.integer(input$ripley_r_steps %||% 60L))
      shiny::validate(shiny::need(is.finite(rmax_c) && rmax_c > 0, "Could not form a common r-grid."))
      rr <- seq(0, rmax_c, length.out = steps)

      for (g in unique(res$group)) {
        files <- unique(res$dataset[res$group == g])
        mat <- vapply(files, function(f) {
          d <- res[res$dataset == f, ]
          stats::approx(d$r, d$Lr_minus_r, xout = rr, rule = 2)$y
        }, numeric(length(rr)))
        if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
        meanc <- smooth_vec(rr, rowMeans(mat, na.rm = TRUE))
        col   <- unname(pal[g]); nfl <- length(files)

        if (band != "none" && nfl >= 2) {
          if (band == "iqr") {
            lo <- apply(mat, 1, stats::quantile, 0.25, na.rm = TRUE)
            hi <- apply(mat, 1, stats::quantile, 0.75, na.rm = TRUE)
          } else {
            sdv <- apply(mat, 1, stats::sd, na.rm = TRUE)
            lo <- rowMeans(mat, na.rm = TRUE) - sdv; hi <- rowMeans(mat, na.rm = TRUE) + sdv
          }
          fig <- fig %>%
            plotly::add_trace(x = rr, y = hi, type = "scatter", mode = "lines",
              line = list(width = 0), legendgroup = g, showlegend = FALSE, hoverinfo = "skip") %>%
            plotly::add_trace(x = rr, y = lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = hex2rgba(col, 0.18), line = list(width = 0),
              legendgroup = g, showlegend = FALSE, hoverinfo = "skip")
        }
        fig <- fig %>% plotly::add_trace(x = rr, y = meanc, type = "scatter", mode = "lines",
          name = sprintf("%s (n=%d)", g, nfl), legendgroup = g,
          line = list(color = col, width = 2.5))
      }
    }

    fig %>%
      plotly::add_trace(
        x = range(res$r), y = c(0, 0), type = "scatter", mode = "lines",
        name = "CSR (L−r = 0)",
        line = list(color = "black", width = 1, dash = "dash"), showlegend = TRUE) %>%
      plotly::layout(
        xaxis  = list(title = "r (nm)"),
        yaxis  = list(title = "L(r) − r (nm)"),
        hovermode = "x unified",
        legend = list(groupclick = "togglegroup"),
        paper_bgcolor = "white", plot_bgcolor  = "white")
  })

  output$download_ripley <- downloadHandler(
    filename = function() paste0("ripley_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- ripley_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ── Two-Color: Cross-Channel NN ──────────────────────────────────────────
  cross_nn_results <- reactiveVal(NULL)

  observeEvent(input$run_cross_nn, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

    unit <- input$cross_nn_unit %||% "tid"
    k    <- max(1L, as.integer(input$cross_nn_k %||% 1L))
    needs_cluster  <- unit %in% c("cluster", "cluster2_cluster")
    needs_cluster2 <- unit %in% c("cluster2_tid", "cluster2_cluster")
    base_unit <- if (unit %in% c("tid","cluster2_tid")) "tid" else "cluster"

    two_color_mode <- isTRUE(input$enable_two_color)
    mplex_pair     <- mplex_cross_pair()

    results <- list()
    withProgress(message = "Cross-channel NN...", value = 0, {
      for (i in seq_along(store)) {
        incProgress(1 / length(store))
        nm <- names(store)[i]
        df <- apply_meta_filters_batch(store[[nm]], input)
        if (!nrow(df)) next

        if (!is.null(mplex_pair)) {
          df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
          if (nrow(df)) df <- apply_two_color_clustering_pipeline(df, input)
        } else if (two_color_mode) {
          df <- apply_two_color_clustering_pipeline(df, input)
        } else {
          df <- assign_analysis_group(df, input)
          if (needs_cluster) df <- apply_clustering_pipeline(df, input)
        }
        if (!nrow(df)) next
        if (!all(c(1L, 2L) %in% df$analysis_group)) next

        if (needs_cluster2) {
          if (!"cluster2" %in% names(df) || !any(!is.na(df$cluster2) & df$cluster2 > 0L)) next
          grp_ids <- unique(df$cluster2[!is.na(df$cluster2) & df$cluster2 > 0L])
          d_list <- lapply(grp_ids, function(g2) {
            sub <- dplyr::filter(df, cluster2 == g2)
            if (!all(c(1L,2L) %in% sub$analysis_group)) return(NULL)
            compute_cross_nn(sub, unit = base_unit, k = k)
          })
          d <- do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0, d_list))
        } else {
          d <- compute_cross_nn(df, unit = base_unit, k = k)
        }
        if (!is.null(d) && nrow(d)) results[[nm]] <- cbind(dataset = nm, d)
      }
    })

    if (!length(results)) {
      showNotification("No cross-channel data found. Check that both channels are present and analysis_group is set correctly.", type = "warning")
      cross_nn_results(NULL)
      return()
    }
    cross_nn_results(do.call(rbind, results))
  })

  output$cross_nn_hist <- renderPlotly({
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

    fig %>% plotly::layout(
      barmode = "overlay",
      xaxis = list(title = "Distance (nm)"),
      yaxis = list(title = if (hist_type == "density") "Density" else "Count"),
      title = "Cross-Channel Nearest-Neighbour Distances",
      paper_bgcolor = "white", plot_bgcolor = "white"
    )
  })

  output$coloc_curve <- renderPlotly({
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
  })

  output$download_cross_nn <- downloadHandler(
    filename = function() paste0("cross_nn_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- cross_nn_results()
      shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ── Two-Color: Bivariate Ripley's L12 ────────────────────────────────────
  ripley12_results <- reactiveVal(NULL)

  observeEvent(input$run_ripley12, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

    unit      <- input$ripley12_unit %||% "tid"
    use_3d    <- isTRUE(input$ripley12_use_3d)
    r_max_raw <- as.numeric(input$ripley12_r_max %||% 0)
    r_steps   <- as.integer(input$ripley12_r_steps %||% 60L)
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

    results <- list()
    withProgress(message = "Bivariate Ripley L12...", value = 0, {
      for (i in seq_along(store)) {
        incProgress(1 / length(store))
        nm <- names(store)[i]
        df <- apply_meta_filters_batch(store[[nm]], input)
        if (!nrow(df)) next

        if (!is.null(mplex_pair)) {
          df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
          if (nrow(df)) df <- apply_two_color_clustering_pipeline(df, input)
        } else if (two_color_mode) {
          df <- apply_two_color_clustering_pipeline(df, input)
        } else {
          df <- assign_analysis_group(df, input)
          if (needs_cluster) df <- apply_clustering_pipeline(df, input)
        }
        if (!nrow(df)) next
        if (!all(c(1L, 2L) %in% df$analysis_group)) next

        r_seq <- if (r_max_raw > 0) seq(0, r_max_raw, length.out = r_steps) else NULL

        pts1 <- get_pts_for_grp(df, 1L)
        pts2 <- get_pts_for_grp(df, 2L)
        lr <- if (nrow(pts1) && nrow(pts2)) {
          compute_ripley_L12(
            x1=pts1$x, y1=pts1$y, z1=if(use_3d) pts1$z else NULL,
            x2=pts2$x, y2=pts2$y, z2=if(use_3d) pts2$z else NULL,
            r_seq=r_seq, use_3d=use_3d
          )
        } else NULL
        if (is.null(lr) || !nrow(lr)) next
        results[[nm]] <- cbind(dataset = nm, lr)
      }
    })

    if (!length(results)) {
      showNotification("No results. Ensure both channels are present.", type = "warning")
      ripley12_results(NULL)
      return()
    }
    ripley12_results(do.call(rbind, results))
  })

  output$ripley12_plot <- renderPlotly({
    res <- ripley12_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run Bivariate Ripley L12 first."))

    aggregate <- isTRUE(input$ripley12_aggregate)
    do_smooth <- isTRUE(input$ripley12_smooth)
    span_val  <- as.numeric(input$ripley12_smooth_span %||% 0.2)
    smooth_y  <- function(sub) {
      if (!do_smooth || nrow(sub) < 5L) return(sub$Lr_minus_r)
      tryCatch(predict(loess(Lr_minus_r ~ r, data = sub, span = span_val)),
               error = function(e) sub$Lr_minus_r)
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
        fig <- fig %>% plotly::add_trace(
          x = sub$r, y = smooth_y(sub),
          type = "scatter", mode = "lines", name = nm,
          line = list(color = pal[nm], width = 2)
        )
      }
    }

    fig %>%
      plotly::add_trace(
        x = range(res$r), y = c(0, 0),
        type = "scatter", mode = "lines",
        name = "Random (L−r = 0)",
        line = list(color = "black", width = 1, dash = "dash"),
        showlegend = TRUE
      ) %>%
      plotly::layout(
        xaxis = list(title = "r (nm)"),
        yaxis = list(title = "L₁₂(r) − r (nm)"),
        title = "Bivariate Ripley’s L₁₂(r)",
        hovermode = "x unified",
        paper_bgcolor = "white", plot_bgcolor = "white"
      )
  })

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
    if (!is.null(mplex_pair)) {
      df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
      if (nrow(df)) df <- apply_two_color_clustering_pipeline(df, input)
    } else if (isTRUE(input$enable_two_color)) {
      df <- apply_two_color_clustering_pipeline(df, input)
    } else {
      df <- assign_analysis_group(df, input)
      if (needs_cluster) df <- apply_clustering_pipeline(df, input)
    }
    if (!nrow(df) || !all(c(1L, 2L) %in% df$analysis_group)) return(NULL)
    df
  }

  # ── Pair Cross-Correlation g₁₂(r) ────────────────────────────────────────
  pcf12_results <- reactiveVal(NULL)

  observeEvent(input$run_pcf12, {
    store <- smlm_file_store()
    if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

    unit          <- input$pcf12_unit %||% "tid"
    use_3d        <- isTRUE(input$pcf12_use_3d)
    r_max_raw     <- as.numeric(input$pcf12_r_max %||% 0)
    r_steps       <- as.integer(input$pcf12_r_steps %||% 60L)
    needs_cluster <- unit == "cluster"

    results <- list()
    withProgress(message = "Pair cross-correlation g₁₂...", value = 0, {
      for (i in seq_along(store)) {
        incProgress(1 / length(store))
        nm <- names(store)[i]
        df <- prep_tc_df(store[[nm]], input, needs_cluster)
        if (is.null(df)) next
        pts1 <- get_tc_pts(df, 1L, unit)
        pts2 <- get_tc_pts(df, 2L, unit)
        if (!nrow(pts1) || !nrow(pts2)) next
        r_seq <- if (r_max_raw > 0) seq(0, r_max_raw, length.out = r_steps) else NULL
        g <- compute_pcf12(pts1$x, pts1$y, if (use_3d) pts1$z else NULL,
                           pts2$x, pts2$y, if (use_3d) pts2$z else NULL,
                           r_seq = r_seq, use_3d = use_3d)
        if (nrow(g)) results[[nm]] <- cbind(dataset = nm, g)
      }
    })

    if (!length(results)) {
      showNotification("No results. Check that both channels are present.", type = "warning")
      pcf12_results(NULL); return()
    }
    pcf12_results(do.call(rbind, results))
  })

  output$pcf12_plot <- renderPlotly({
    res <- pcf12_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run g₁₂ first."))

    aggregate <- isTRUE(input$pcf12_aggregate)
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
        fig <- fig %>% plotly::add_trace(x = sub$r, y = sub$g12,
          type = "scatter", mode = "lines", name = nm,
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
        title = "Pair Cross-Correlation Function g₁₂(r)",
        hovermode = "x unified",
        paper_bgcolor = "white", plot_bgcolor = "white"
      )
  })

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
    res <- compute_kde_coloc(pts1[[ax[1]]], pts1[[ax[2]]],
                             pts2[[ax[1]]], pts2[[ax[2]]], bw, n_grid)
    if (is.null(res)) {
      showNotification("KDE computation failed.", type = "warning")
      kde_coloc_results(NULL); return()
    }
    res$dataset <- nm
    res$ax      <- ax
    kde_coloc_results(res)
  })

  output$kde_coloc_plot <- renderPlotly({
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
        title = sprintf("KDE — %s  [%s]  (Ch1 / Ch2)",
                        r$dataset, paste(ax, collapse = "/")),
        paper_bgcolor = "white", plot_bgcolor = "white"
    )
  })

  output$kde_pearson_display <- renderUI({
    r <- kde_coloc_results()
    if (is.null(r)) return(NULL)
    tags$div(
      style = "margin-top:6px;",
      tags$b("Density-density Pearson r:"),
      tags$p(sprintf("%.3f  (bandwidth = %.1f nm)", r$pearson_r, r$bandwidth),
             style = "font-size:13px; margin:2px 0;")
    )
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





