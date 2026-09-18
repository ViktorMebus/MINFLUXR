## ─────────────────────────────────────────────────────────────────────────────
## ui_tracking.R  -  Tracking tab UI (assigned to tracking_ui)
## All input/output IDs prefixed with trk_ to avoid collisions with SMLM IDs.
## ─────────────────────────────────────────────────────────────────────────────

# NOTE: .method_box() (the collapsible "Method & MINFLUX adaptation" infobox used
# throughout this file) now lives in helpers_shared.R, so the SMLM tab can use the
# identical widget - helpers_shared.R is sourced before both ui_*.R files.

tracking_ui <- fluidPage(
  # Tracking's shinyjs calls (reset, delay, toggleState) previously worked only
  # because smlm_ui declares useShinyjs() into the same navbarPage document.
  # Declared here too so this tab does not depend on the SMLM tab existing;
  # Shiny de-duplicates the html dependency, so the double call costs nothing.
  useShinyjs(),
  tags$head(tags$style(HTML("
    details.sim-group { margin: 2px 0 6px 0; }
    details.sim-group > summary {
      cursor: pointer; font-size: 1.05em; font-weight: 600;
      margin: 6px 0; padding: 2px 0; list-style: none;
    }
    details.sim-group > summary::-webkit-details-marker { display: none; }
    details.sim-group > summary::before {
      content: \"\\25B8\"; color: #888; margin-right: 6px;
    }
    details.sim-group[open] > summary::before { content: \"\\25BE\"; }
    /* A panel that does not apply to the loaded data: dimmed but still visible,
       so it stays discoverable and can say why it is inactive. */
    details.sim-group.dimmed > summary { color: #9aa0a6; }
    details.sim-group.dimmed { opacity: 0.55; }
    /* A radio choice that does not apply to the loaded data. Shiny wraps each
       choice label in a <span> right after its <input>, so the disabled input
       can grey its own text without touching the other choices. */
    .radio input:disabled + span,
    .radio-inline input:disabled + span,
    .shiny-input-radiogroup input:disabled + span { color: #9aa0a6; cursor: not-allowed; }
  "))),
  titlePanel("MINFLUXR - Tracking"),

  tabsetPanel(
    id = "trk_tabs",

    # ── Display ───────────────────────────────────────────────────────────────
    tabPanel("Display",
      # Sticky Plot button bar
      tags$div(
        style = paste(
          "position: sticky; top: 0; z-index: 200;",
          "background: #f8f9fa; border-bottom: 1px solid #dee2e6;",
          "padding: 6px 15px; display: flex; align-items: center; gap: 12px;"
        ),
        actionButton("trk_plot_data", "Plot Data",
                     icon  = icon("play"),
                     class = "btn-primary btn-sm"),
        uiOutput("trk_bar_label")
      ),
      sidebarLayout(
        sidebarPanel(
          style = "height: 90vh; overflow-y: auto;",

          tags$details(class = "sim-group", open = NA,
          tags$summary(tags$b("Data Upload")),
          fileInput("trk_file", "Upload MINFLUX File(s) (.npy, .json or .csv)",
                    accept = c(".csv", ".csv.gz", ".rds", ".fst", ".qs", ".json", ".npy"),
                    multiple = TRUE),
          helpText(style = "font-size:0.82em; color:#666; margin-top:-6px;",
            "Multiple files merge into one dataset; track IDs are prefixed F1_, F2_, … so they cannot collide."),
          # Chunking splits the dataset into blocks of whole TIDs so a modest machine
          # renders one block at a time. Shown in both view modes: the chunk size
          # applies at upload, and "Limit Full tracks" below makes the selected chunk
          # scope the static view as well as playback.
          numericInput(
            "trk_tids_per_chunk",
            "TIDs per chunk (for big-file splitting)",
            value = TIDS_PER_CHUNK_DEFAULT,
            min = 5, step = 10
          ),
          helpText(style = "font-size:0.82em; color:#666; margin-top:-6px;",
                   "Fewer TIDs per chunk = smoother rendering on a slow machine. ",
                   "Applies on upload and on every Plot Data."),
          uiOutput("trk_chunk_ui"),
          checkboxInput("trk_chunk_scope_full",
                        "Limit Full tracks to the selected chunk", FALSE),
          helpText(style = "font-size:0.82em; color:#666; margin-top:-6px;",
                   "Off: Full tracks draws every track and the chunk affects playback only. ",
                   "Filtering and analysis always use all TIDs."),
          uiOutput("trk_chunk_stat"),
          ),  # /Data Upload
          tags$hr(),

          # ── Axial (z) calibration ───────────────────────────────────────────
          # Mirrors the SMLM tab: a property of the DATASET, set next to the
          # upload it describes, not buried in the display controls. Files are
          # read uncorrected (.npy/.json at z_scale = 1, CSVs stripped of
          # whatever factor they arrived with), so this is the only axial
          # correction applied anywhere in the tracking pipeline.
          tags$details(class = "sim-group", id = "trk_z_scale_group",
          tags$summary(tags$b("Axial (z) calibration")),
          numericInput("trk_z_scale", "Axial (z) scaling factor",
                       value = TRK_Z_SCALE_DEFAULT, min = 0.05, max = 3, step = 0.01),
          helpText(style = "font-size:0.82em; color:#666; margin-top:-6px;",
            "Refractive-index mismatch between the immersion and the mounting medium. ",
            tags$b("Rescales every z-dependent result"), " - check the mountant's index."),
          uiOutput("trk_z_scale_status"),
          helpText(style = "font-size:0.82em; color:#666;",
            "Applies on Plot Data. A file that records its own factor overrides this."),
          ),  # /Axial calibration
          tags$hr(),

          tags$hr(),

          tags$details(class = "sim-group", open = NA,
          tags$summary(tags$b("Track Cleaning & Filtering")),
          helpText(style = "font-size:0.8em; color:#666; margin:2px 0 6px 0;",
            "Order: Split → Trim ends → Quality gate → Metadata → Length → Duration. ",
            "Cleaning runs on the raw acquisition timing."),
          uiOutput("trk_lifetime_stat"),
          checkboxInput("trk_show_metadata_hist", "Show Filtering Controls", FALSE),
          conditionalPanel(
            condition = "input.trk_show_metadata_hist == true",


            checkboxInput("trk_split_gaps_enable",
                          "Split tracks at time gaps", FALSE),
            conditionalPanel(
              condition = "input.trk_split_gaps_enable == true",
              checkboxInput("trk_split_adaptive",
                            "Adaptive threshold", FALSE),
              conditionalPanel(
                condition = "input.trk_split_adaptive == false",
                numericInput("trk_split_gap_cycles",
                             "Gap threshold (cycles) - split if interval exceeds this",
                             value = 7, min = 2, step = 1)
              ),
              conditionalPanel(
                condition = "input.trk_split_adaptive == true",
                uiOutput("trk_split_adaptive_hint")
              ),
              # Applies in BOTH modes: segments shorter than this are discarded, so a
              # split never manufactures 2–3-localization stubs. Also keeps the plot
              # fast by not emitting a trace per tiny fragment.
              numericInput("trk_split_min_len",
                           "Min localizations per segment (shorter segments discarded)",
                           value = 50, min = 1, step = 1),
              uiOutput("trk_split_gap_ms_hint"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "An interval spanning more cycles than this is a gap, and the track splits there ",
                "(track 42 becomes 42_s1, 42_s2, …). ", tags$b("Adaptive"), " sets the threshold from ",
                "the data as median + SD + 1 of the per-step cycle counts, capped by maxOffTime.")
            ),

            tags$hr(),

            checkboxInput("trk_endtrim_enable",
                          "Trim degraded track ends (bleaching tails)", FALSE),
            conditionalPanel(
              condition = "input.trk_endtrim_enable == true",
              checkboxInput("trk_endtrim_dark", "Trim dark-gap-separated ends", TRUE),
              conditionalPanel(
                condition = "input.trk_endtrim_dark == true",
                numericInput("trk_endtrim_dark_cycles",
                             "Dark-gap threshold (cycles)", value = 3, min = 2, step = 1)
              ),
              checkboxInput("trk_endtrim_efo", "Trim dim (low-EFO) ends", TRUE),
              conditionalPanel(
                condition = "input.trk_endtrim_efo == true",
                numericInput("trk_endtrim_efo_frac",
                             "EFO threshold (× track median)", value = 0.5, min = 0.05, max = 1, step = 0.05)
              ),
              uiOutput("trk_endtrim_hint"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Peels localizations off each end while they are reached across a dark gap ≥ the ",
                "threshold or are dim, until a solid one is reached. Targets the photobleaching tail, ",
                "which otherwise inflates MSD tails and fakes exports at track ends.")
            ),

            tags$hr(),

            checkboxInput("trk_quality_enable",
                          "Track-continuity gate (drop dark-cycle-heavy tracks)", FALSE),
            conditionalPanel(
              condition = "input.trk_quality_enable == true",
              sliderInput("trk_quality_darktol",
                          "Dark-cycle tolerance (half-life, cycles)",
                          min = 0.5, max = 10, value = 3, step = 0.5),
              sliderInput("trk_quality_min_integ",
                          "Min continuity score",
                          min = 0, max = 1, value = 0.6, step = 0.05),
              uiOutput("trk_quality_hint"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Each step is weighted 2^(−dark/half-life), so short gaps barely hurt and long ones do. ",
                "The per-track mean is the continuity score; tracks below the threshold are dropped.")
            ),

            tags$hr(),

            tags$b("Metadata (per-localization) filters", style = "font-size:0.9em;"),
            selectInput(
              "trk_hist_param", "Variable for histogram:",
              choices = c(
                "efo"                    = "efo",
                "cfr"                    = "cfr",
                "dcr"                    = "dcr",
                "tim"                    = "tim",
                "eco"                    = "eco",
                "Dwell"                  = "dwell",
                "Fluo"                   = "fluo",
                "Trace ID length"        = "tid_length",
                "Trace ID duration (ms)" = "tid_duration"
              ),
              selected = "efo"
            ),
            plotlyOutput("trk_meta_hist", height = "260px"),
            sliderInput("trk_hist_bins", "Histogram bins", min = 10, max = 500, value = 80, step = 5),
            checkboxInput("trk_hist_show_raw", "Show unfiltered data (grey)", value = TRUE),
            uiOutput("trk_hist_range_ui"),

            tags$hr(),

            uiOutput("trk_efo_ui"),
            uiOutput("trk_cfr_ui"),
            uiOutput("trk_dcr_ui"),
            uiOutput("trk_tim_ui"),
            uiOutput("trk_eco_ui"),
            uiOutput("trk_dwell_ui"),
            uiOutput("trk_fluo_ui"),

            tags$hr(),
            uiOutput("trk_tid_len_ui"),

            tags$hr(),

            uiOutput("trk_dur_filter_ui"),

            actionButton("trk_reset_filters", "Reset ranges", style = "width:100%;"),
            helpText(style = "font-size:0.82em; color:#666;",
                     "Filters apply on ", tags$b("Plot Data"), ", which also re-chunks the surviving ",
                     "tracks. Reset ranges restores the full uploaded set."),

            tags$hr(style = "margin:8px 0;"),
            tags$b("Spatial crop", style = "font-size:0.9em;"),
            uiOutput("trk_spatial_info"),
            fluidRow(
              column(7, actionButton("trk_crop_to_view", "Crop to view",
                                     icon = icon("crop-alt"), class = "btn-info btn-sm",
                                     style = "width:100%")),
              column(5, actionButton("trk_clear_spatial", "Clear",
                                     class = "btn-default btn-sm", style = "width:100%"))
            ),
            numericInput("trk_crop_pct", "Min % of TID localizations inside box",
                         value = 50, min = 1, max = 100, step = 5),
            helpText(style = "font-size:0.82em; color:#666;",
                     "Zoom on the plot, then Crop to view. Reset ranges also clears this."),
            tags$hr(style = "margin:8px 0;"),

            downloadButton("trk_download_filtered_csv", "Export filtered data (.csv)"),
            helpText(style = "font-size:0.82em; color:#666;",
                     "Canonical MINFLUX .csv. Also the input for Simulator → Diffusion/Tracking → ",
                     tags$i("Calibrate from real MINFLUX data"), "."),
            tags$hr(),
            downloadButton("trk_download_filter_settings", "Save filter settings (.json)"),
            fileInput("trk_load_filter_settings", "Load filter settings (.json)", accept = ".json")
          ),
          ),  # /Metadata Filtering
          tags$hr(),

          radioButtons("trk_view_mode_main", NULL,
            choices = c("Playback" = "playback", "Full tracks" = "full"),
            selected = "full", inline = TRUE),

          conditionalPanel(
            condition = "input.trk_view_mode_main == 'full'",
            checkboxInput("trk_full_connect", "Connect each trace by time (lines)", FALSE),
            tags$hr()
          ),

          conditionalPanel(
            condition = "input.trk_view_mode_main == 'playback'",

            tags$details(class = "sim-group", open = NA,
            tags$summary(tags$b("Playback Controls")),

            fluidRow(
              column(3, actionButton("trk_btn_play",  "Play ▶",  class = "btn-success btn-sm")),
              column(3, actionButton("trk_btn_pause", "Pause ⏸", class = "btn-warning btn-sm")),
              column(3, actionButton("trk_btn_stop",  "Stop ⏹",  class = "btn-danger btn-sm")),
              column(3, actionButton("trk_auto_tune", "Auto-tune",    class = "btn-info btn-sm"))
            ),
            tags$div(style = "margin-top: 4px;"),
            uiOutput("trk_frame_slider_ui"),
            sliderInput("trk_play_speed", "Playback speed (×)",
                        min = 0.25, max = 8, value = 1, step = 0.25, width = "100%"),
            uiOutput("trk_frame_meta"),
            tags$hr(),
            checkboxInput("trk_show_playback_params", "Show Playback Parameters", FALSE),
            conditionalPanel(
              condition = "input.trk_show_playback_params == true",
              selectInput(
                "trk_play_mode", "Playback mode",
                choices = c("Simultaneous (normalized time)" = "simul",
                            "Sequential per TID"             = "seq"),
                selected = "simul"
              ),
              numericInput("trk_play_total",     "Total playback window (ms)", value = 4000, min = 0.001, step = 0.001),
              numericInput("trk_bin_dt",         "Bin width (ms)",             value = 5,   min = 0.0001, step = 0.001),
              helpText("Frame count is derived automatically: window ÷ bin width."),

              checkboxInput("trk_playback_ctrl", "Show/Hide Playback Controls", TRUE),
              conditionalPanel(
                condition = "input.trk_playback_ctrl == true",
                checkboxInput("trk_show_moving_trails", "Animate moving trails", FALSE),
                conditionalPanel(
                  condition = "input.trk_show_moving_trails == true",
      
                ),
                numericInput("trk_trail_len",     "Trail length (ms)",                  value = 20, min = 0,    step = 1),
                uiOutput("trk_trail_warn"),
                checkboxInput("trk_show_advanced", "Show advanced options", FALSE),
                conditionalPanel(
                  condition = "input.trk_show_advanced == true",
                  numericInput("trk_target_markers", "Target markers per frame",                value = 5000, min = 100, step = 500),
                  numericInput("trk_per_tid_cap",  "Per-TID cap per frame (0 = off)",          value = 0,   min = 0, step = 50),
                  numericInput("trk_trail_nmax",   "Max points per trail (per TID, per frame)", value = 50,  min = 5, step = 5),
                  numericInput("trk_max_pts_tid",  "Pre-decimate max points per TID",           value = 400, min = 50, step = 50)
                )
              ),

              downloadButton("trk_download_playback_settings", "Save playback settings (.json)"),
              fileInput("trk_load_playback_settings", "Load playback settings (.json)", accept = ".json"),
              tags$hr()
            )
            ),  # /Playback Controls
          ),
          tags$hr(),

          tags$details(class = "sim-group", open = NA,
          tags$summary(tags$b("Display Options")),
          # Static (not renderUIs) so the disabled state applied by the planar-data
          # observer in server_tracking.R cannot be wiped by a re-render. On 2D
          # data these are greyed rather than removed, matching the axial
          # calibration panel - the control stays visible and says why it is off.
          radioButtons("trk_view_mode", NULL,
                       choices = c("2D" = "2d", "3D" = "3d"),
                       selected = "2d", inline = TRUE),
          conditionalPanel(
            condition = "input.trk_view_mode == '2d'",
            radioButtons("trk_view_plane", "Viewing plane",
                         choices = c("XY", "XZ", "YZ"), selected = "XY", inline = TRUE)
          ),
          uiOutput("trk_z_axes_note"),
          numericInput("trk_pt_size", "Marker size", value = 6, min = 1),
          radioButtons("trk_color_mode", "Color localizations by",
                       choices = c("Track ID" = "tid", "Dark cycles (integration)" = "dark",
                                   "Speed (drift velocity)" = "speed"),
                       selected = "tid"),
          conditionalPanel(
            condition = "input.trk_color_mode == 'speed'",
            numericInput("trk_color_speed_win", "Drift window W (ms)",
                         value = 20, min = 0.1, step = 1),
            uiOutput("trk_speed_win_note"),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
              "Blue (slow) → red (fast) by ", tags$b("drift velocity"), " = |net displacement| / span ",
              "over a centred window W; grey = window too short (track ends). ",
              tags$b("Keep W well below a typical track duration"), ", or a whole track shares one ",
              "window and goes flat. Diffusion alone gives v ≈ √(4D/W) and noise a floor ≈ √2·σ/W, ",
              "so only sustained red well above those is directed transport.")
          ),
          conditionalPanel(
            condition = "input.trk_color_mode == 'dark'",
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Blue = integrated, reddening with the number of dark cycles preceding each ",
            "localization (cycle time from Acquisition Parameters).")
          ),
          radioButtons("trk_dark_show", "Show localizations",
                       choices = c("All" = "all", "Integrated only" = "integrated",
                                   "Dark-cycle only" = "dark"),
                       selected = "all"),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Markers only - trails keep the full path for context."),
          conditionalPanel(
            condition = "input.trk_view_mode_main == 'playback'",
            checkboxInput("trk_show_full_tracks_faint", "Show full tracks faintly", TRUE),
            helpText("Markers & (optional) static lines use one color per TID.")
          ),
          checkboxInput("trk_full_raw", "Show ALL raw tracks faintly (background)", FALSE),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Faint grey overlay of every uploaded track, so you can see what the filters removed."),
          checkboxInput("trk_dark_mode", "Dark mode", value = FALSE),
          checkboxInput("trk_hide_grid", "Hide gridlines and axis titles", FALSE),
          checkboxInput("trk_show_scale_bar", "Add scale bar", value = FALSE),
          conditionalPanel(
            condition = "input.trk_show_scale_bar == true",
            numericInput("trk_scale_bar_length",    "Scale bar length (nm)",    value = 100, min = 1, step = 10),
            numericInput("trk_scale_bar_thickness", "Scale bar thickness (nm)", value = 10,  min = 1, step = 1),
            numericInput("trk_scale_font_size",     "Font size",                value = 18,  min = 1, step = 1)
          ),
          ),  # /Display Options
          tags$hr(),

          tags$details(class = "sim-group",
          tags$summary(tags$b("Background Image Overlay")),
          checkboxInput("trk_show_confocal", "Show overlay settings", FALSE),
          conditionalPanel(
            condition = "input.trk_show_confocal == true",
            tags$div(id = "trk_overlay_form",
                     style = "border:0; padding:0; margin:0; background:none; box-shadow:none;",
                     checkboxInput("trk_show_bg", "Show background image", FALSE),
                     fileInput("trk_bg_image", "Upload image",
                               accept = c(".png",".jpg",".jpeg",".tif",".tiff")),
                     numericInput("trk_px_size",  "Pixel size (units/pixel)",    value = 10, min = 0.0001),
                     numericInput("trk_img_x0",   "Image origin X (data units)", value = 0),
                     tags$div(
                       style = "display:flex; gap:3px; margin:-4px 0 8px 0; flex-wrap:wrap;",
                       tags$small("X:", style = "line-height:26px; color:#555; margin-right:2px;"),
                       actionButton("trk_nx_m1000", "−1k",  class = "btn-xs btn-default"),
                       actionButton("trk_nx_m100",  "−100", class = "btn-xs btn-default"),
                       actionButton("trk_nx_m10",   "−10",  class = "btn-xs btn-default"),
                       actionButton("trk_nx_p10",   "+10",       class = "btn-xs btn-default"),
                       actionButton("trk_nx_p100",  "+100",      class = "btn-xs btn-default"),
                       actionButton("trk_nx_p1000", "+1k",       class = "btn-xs btn-default")
                     ),
                     numericInput("trk_img_y0",   "Image origin Y (data units)", value = 0),
                     tags$div(
                       style = "display:flex; gap:3px; margin:-4px 0 8px 0; flex-wrap:wrap;",
                       tags$small("Y:", style = "line-height:26px; color:#555; margin-right:2px;"),
                       actionButton("trk_ny_m1000", "−1k",  class = "btn-xs btn-default"),
                       actionButton("trk_ny_m100",  "−100", class = "btn-xs btn-default"),
                       actionButton("trk_ny_m10",   "−10",  class = "btn-xs btn-default"),
                       actionButton("trk_ny_p10",   "+10",       class = "btn-xs btn-default"),
                       actionButton("trk_ny_p100",  "+100",      class = "btn-xs btn-default"),
                       actionButton("trk_ny_p1000", "+1k",       class = "btn-xs btn-default")
                     ),
                     tags$div(
                       style = "display:flex; gap:5px; margin-bottom:8px;",
                       actionButton("trk_auto_foci",     "Auto: foci",
                                    icon = icon("bullseye"), class = "btn-sm btn-info"),
                       actionButton("trk_auto_boundary", "Auto: boundary",
                                    icon = icon("draw-polygon"), class = "btn-sm btn-info")
                     ),
                     checkboxInput("trk_show_diag", "Show detection diagnostics", FALSE),
                     tags$details(
                       tags$summary(style = "cursor:pointer; font-size:12px; color:#555; margin:3px 0 6px 0;",
                                    "Foci detection settings"),
                       tags$div(style = "padding:4px 0 4px 8px; border-left:2px solid #ddd;",
                         sliderInput("trk_foci_n_peaks",       "Max foci to find",                               min=1,  max=10,  value=5,   step=1),
                         sliderInput("trk_foci_bright_pct",    "Confocal: sensitivity (70=permissive → 99=strict)", min=70, max=99,  value=90,  step=1),
                         sliderInput("trk_foci_border_pct",    "Confocal: edge mask (%)",                        min=2,  max=25,  value=10,  step=1),
                         sliderInput("trk_foci_eps_pct",       "Confocal: focus radius (% of image)",            min=1,  max=10,  value=3,   step=0.5),
                         sliderInput("trk_foci_density_thr",    "MINFLUX: density threshold (% of peak)",        min=5,  max=60,  value=20,  step=5),
                         sliderInput("trk_foci_density_cell",   "MINFLUX: cell size (nm)",                       min=20, max=500, value=100, step=10),
                         sliderInput("trk_foci_density_smooth", "MINFLUX: smoothing radius (nm; 0=off)",         min=0,  max=500, value=100, step=25)
                       )
                     ),
                     tags$p("Pixel size must be set first. Results are a starting point - nudge to refine.",
                            style = "font-size:11px; color:#888; margin:0 0 6px 0;"),
                     checkboxInput("trk_img_invert_y", "Image origin at top-left (invert Y)", TRUE),
                     checkboxInput("trk_img_flip_x",   "Mirror horizontally (left↔right)", TRUE),
                     sliderInput("trk_img_opacity", "Image opacity", min = 0, max = 1, value = 0.6, step = 0.05),
                     selectInput("trk_img_lut", "Colour (LUT)",
                                 choices = c("Original (as uploaded)" = "original",
                                             "Grayscale"              = "gray",
                                             "Green"                  = "green",
                                             "Magenta"                = "magenta",
                                             "Yellow"                 = "yellow",
                                             "Cyan"                   = "cyan",
                                             "Red"                    = "red",
                                             "Blue"                   = "blue",
                                             "Orange"                 = "orange",
                                             "Custom colour…"         = "custom"),
                                 selected = "original"),
                     conditionalPanel(
                       condition = "input.trk_img_lut == 'custom'",
                       colourpicker::colourInput("trk_img_lut_custom", "Custom colour",
                                                 value = "#00FF00")
                     ),
                     tags$p("Recolours the image to match a fluorophore: converts to grayscale intensity, then tints black→colour.",
                            style = "font-size:11px; color:#888; margin:-4px 0 8px 0;"),
                     sliderInput("trk_img_brightness", "Image brightness (%)",
                                 min = 10, max = 300, value = 100, step = 5),
                     sliderInput("trk_img_contrast", "Image contrast",
                                 min = 0, max = 90, value = 0, step = 5),
                     selectInput("trk_img_rotation", "Rotate image",
                                 choices = c("0","90","180","270"), selected = "180"),
                     downloadButton("trk_download_overlay_settings", "Download Overlay Settings"),
                     fileInput("trk_load_overlay_settings", "Load Overlay Settings (.json)",
                               accept = ".json"),
                     actionButton("trk_reset_overlay", "Reset Overlay", icon = icon("undo"))
            )
          ),
          ),  # /Background Image Overlay
          tags$hr(),

          tags$details(class = "sim-group",
          tags$summary(tags$b("Dwell-time Heatmap")),
          checkboxInput("trk_show_heatmap", "Show dwell-time heatmap", FALSE),
          conditionalPanel(
            condition = "input.trk_show_heatmap == true",
            selectInput(
              "trk_hm_colors", "Heatmap colors",
              choices = c("Magenta glow" = "magenta",
                          "Fire"         = "fire",
                          NUMERIC_PALETTES),
              selected = "fire"
            ),
            sliderInput("trk_hm_clip_pct", "Clipping percentiles",      min = 0, max = 100, value = c(1, 99), step = 0.5),
            sliderInput("trk_hm_gamma",    "Gamma ( <1 boosts faint )", min = 0.30, max = 2.00, value = 1.00, step = 0.05),
            helpText("Click 'Localizations' in the plot legend to hide/show points."),
            numericInput("trk_hm_bins",   "Heatmap grid (N×N)",    value = 400, min = 20, max = 1000, step = 10),
            sliderInput("trk_hm_opacity", "Heatmap opacity",            min = 0, max = 1, value = 0.25, step = 0.05)
          ),
          ),  # /Dwell-time Heatmap
          tags$hr()
        ),
        mainPanel(
          uiOutput("trk_track2D_wrapper")
        )
      )
    ),  # end tabPanel("Visualization")

    # ── Analysis ─────────────────────────────────────────────────────────────
    tabPanel("Analysis",
      tags$br(),

      # ── Data Quality & Calibration (do this first) ───────────────────────────
      fluidRow(
        column(4,
          h4("Data Quality & Calibration"),
          .method_box("amber",
            tagList("Two estimators, from two different channels. ", tags$b("D"), " from the ",
                    "covariance-based estimator, CVE (Vestergaard, Blainey & Flyvbjerg 2014), read from consecutive steps with ",
                    "no MSD fit involved. ", sigloc(), " comes from the MINFLUX photon budget, ",
                    "L", tags$sub("TCP"), "/(a\u221AN) per localization, with a an instrument constant."),
            HTML("<b>D</b> comes from the CVE, restricted to consecutive <b>single-cycle (integrated)</b> pairs so intervals spanning dark cycles stay out of the covariance. D is unaffected by motion blur.<br><b>&sigma;<sub>loc</sub> comes from the photon budget</b>, not from the covariance: &sigma;<sub>i</sub> = L<sub>TCP</sub>/(a·√N<sub>i</sub>) per localization, with N = ECO and L from Acquisition Parameters, times a small term for the distance the emitter travels during one localization. It is a direct instrument model: one value per localization, and nothing in it depends on how the molecule moves, so confinement, drift and directed motion cannot bias it.<br>The prefactor a defaults to <b>4</b>, the ideal Cramér-Rao value for a doughnut TCP. Background, finite beam depth and estimator losses can only make &sigma; worse, so 4 is a ceiling on a and the default &sigma; is a <b>lower bound</b> on the localization error. An immobile sample measures the other side of the bracket, and is entered under Acquisition Parameters. The beam-size denominator of the doughnut CRLB is reached through the &nu; field and is off by default, since &sigma; &prop; L holds unless an L-series shows otherwise. The same pass recovers the cycle time, which gap-splitting, the track-quality gate and the D-map all read; TSA below gives the full Δt and integrated/dark breakdown."),
            HTML("<b>The rating is amber because of D, not &sigma;.</b> The reported position wanders even on a genuinely immobile target, so a measured D carries an instrumental floor. The wander is photon-driven rather than sample motion, so it shrinks as the photon budget grows and worsens at larger TCP diameters. Measure it on an immobile sample and enter it as the wander floor under Acquisition Parameters, and it is reported next to D. A D of the same order as that floor is not resolved above the instrument.<br>The <b>CVE &sigma;</b> shown beside the photon value is a sanity read, not a measurement. At realistic blur roughly 90% of it is the assumed R⟨Δx²⟩ term, so it tracks √R × rms step and inherits population heterogeneity and drift; on immobile samples the increments come out <i>positively</i> correlated and it fails outright. Under strongly directed motion the covariance is swamped, which the panel flags as degenerate.")),
          actionButton("trk_cve_run", HTML("Estimate &sigma;<sub>loc</sub> &amp; quality"),
                       icon = icon("ruler-combined"),
                       class = "btn-primary", style = "margin-bottom:6px; width:100%;"),
          helpText(style = "font-size:0.82em; color:#666; margin-top:0;",
            "Run this before the panels below, which use ", sigloc(),
            " as their noise floor and the cycle time as their clock."),
          # ── Acquisition parameters (optional cycle-time model) ──────────────
          tags$details(class = "sim-group",
          tags$summary(tags$b("Acquisition Parameters")),
          checkboxInput("trk_hw_enable",
                        "Use hardware model for cycle time (optional)", FALSE),
          checkboxInput("trk_cycle_use_regression",
                        "Fit cycle time from data (EOD/SW regression)", FALSE),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Line fit of the cluster centres (cycle = slope, SW = intercept). Off = cluster median. ",
            "Both are shown below; ▶ marks the active one."),
          numericInput("trk_msd_cycle_ms",
                       "Cycle time manual override (ms, 0 = use auto)",
                       value = 0, min = 0, step = 0.01),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "> 0 overrides both the hardware model and the data-driven estimates."),
          uiOutput("trk_cycle_time_status"),
          numericInput("trk_L_tcp",
                       HTML("TCP diameter L<sub>TCP</sub>, final iteration (nm)"),
                       value = DMAX_L_TCP_DEFAULT, min = 10, max = 2000, step = 10),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "The iteration repeated during tracking. Sets ", sigloc(), " = L/(a·√N) and the ",
            "D", tags$sub("max"), " ceiling. Uploading the sequence under Time Signature ",
            "Analysis fills it in exactly."),
          numericInput("trk_crlb_pref",
                       HTML("&sigma;<sub>loc</sub> prefactor a, in &sigma; = L/(a&middot;&radic;N)"),
                       value = TCP_CRLB_PREFACTOR, min = 0.5, max = 4, step = 0.1),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            tags$b("4 is the ideal bound, and the ceiling."), " Background, beam depth and estimator ",
            "losses can only lower it, so the default makes ", sigloc(), " a lower bound on the ",
            "localization error. A value measured on an immobile sample is an upper bound on ",
            sigloc(), "; entering it reports the two as a range. The porting note below gives the recipe."),
          numericInput("trk_geofactor_base",
                       HTML("patGeoFactor base (nm), L = patGeoFactor &times; base"),
                       value = TCP_L_PER_GEOFACTOR_NM, min = 50, max = 2000, step = 10),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Converts the sequence's patGeoFactor into a TCP diameter. 360 nm reproduces the ",
            "284 / 302 / 151 / 40 nm ladders on Abberior systems, but check it against a diameter ",
            "your own vendor quotes before trusting the auto-filled L."),
          numericInput("trk_beam_fwhm",
                       HTML("Excitation doughnut FWHM &nu; (nm)"),
                       value = TCP_BEAM_FWHM_NM, min = 0, max = 1000, step = 5),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "0 = off, the large-beam limit, in which ", sigloc(), " \u221d L. A finite \u03bd applies the ",
            "beam-size denominator of the doughnut CRLB, which makes ", sigloc(), " rise faster than ",
            "L as L approaches \u03bd. Worth setting only if an L-series on your system shows that rise."),
          numericInput("trk_dw_floor",
                       HTML("Tracker wander floor D<sub>w</sub> (nm²/ms, 0 = off)"),
                       value = TRK_DW_FLOOR_DEFAULT, min = 0, step = 10),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Apparent diffusion the tracking loop produces on an immobile target. Measure it on a ",
            "fixed sample at the same TCP diameter and photon budget you track with. Reported next ",
            "to D, never subtracted silently, because it scales with the photon budget."),
          numericInput("trk_sigma_loc_nm",
                       HTML("Localization precision &sigma;<sub>loc</sub> (nm, 0 = auto)"),
                       value = 0, min = 0, step = 0.5),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "0 = the app-wide value from ", tags$b("Data Quality & Calibration"),
            " (the photon budget), falling back to the CVE and then the MSD-fit intercept. ",
            "> 0 overrides it everywhere."),
          uiOutput("trk_sigma_loc_status"),
          conditionalPanel(
            condition = "input.trk_hw_enable == true",
            selectInput("trk_hw_pattern", "Scan pattern (vertices)",
              choices = c("Hexagonal (6)" = "6", "Triangular (3)" = "3",
                          "Square (4)" = "4", "Other" = "other"),
              selected = "6"),
            conditionalPanel(
              condition = "input.trk_hw_pattern == 'other'",
              numericInput("trk_hw_vertices", "Number of vertices", value = 6, min = 1, step = 1)
            ),
            fluidRow(
              column(6, numericInput("trk_hw_dwell_us", "Dwell (µs)",        value = 100, min = 1, step = 1)),
              column(6, numericInput("trk_hw_eod_us",   "EOD overhead/move (µs)", value = 5.3, min = 0, step = 0.1))
            ),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              tags$b("Cycle time = dwell + EOD overhead × vertices"), ". The EOD overhead is the ",
              "doughnut move between pattern vertices, typically ~5 µs. When enabled this drives ",
              "gap-splitting, the MSD and the Integrated/Dark classification; TSA estimates the ",
              "same quantity from the data as a cross-check.")
          ),
          ),  # /Acquisition Parameters
          tags$hr(style = "margin:8px 0;"),
          tags$hr(style = "margin:8px 0;"),
          tags$b(HTML("Per-track CVE &sigma; (diagnostic)"), style = "font-size:0.9em;"),
          fluidRow(
            column(6, numericInput("trk_sig_min_integ", "Min integrated steps / track",
                                   value = 50, min = 10, step = 10)),
            column(6, numericInput("trk_sig_skip", "Skip first N localizations",
                                   value = 5, min = 0, step = 1))
          ),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "One CVE ", sigloc(), " per track that is long enough. This does not feed the app-wide ",
            sigloc(), ": the blur term dominates it, so the spread shows how much the step size varies ",
            "between tracks rather than how much the precision does."),
          .detail(label = "About these two settings",
            HTML("A short track gives a covariance too noisy to fit, so tracks below the minimum are left out of the per-track figure. Fifty integrated steps is a practical floor; raising it gives cleaner individual estimates from fewer tracks."),
            tags$br(), tags$br(),
            HTML("The first localizations of a MINFLUX track are taken while the beam pattern is still centring on the molecule and scatter more widely than the rest, which would inflate &sigma;<sub>loc</sub>. Skipping the opening few removes them. Set it to 0 to see how large that transient is in your data.")),
        ),
        column(8,
          uiOutput("trk_cve_status"),
          plotlyOutput("trk_sigma_pertrack_plot", height = "220px"),
          tags$details(style = "margin-top:6px;",
            tags$summary(style = "cursor:pointer; font-size:0.86em; color:#555;",
                         "Data-selection funnel (what's excluded, for your methods)"),
            tags$div(style = "padding:4px 0 0 4px;", uiOutput("trk_quality_funnel"))),
          .detail(label = "ℹ Setup for your MINFLUX system (what to measure)",
            HTML("Three values in <b>Acquisition Parameters</b> belong to the individual microscopes. The current defaults are reasonable starting points, yet may be incorrect for your system. Measure them on the system that produced your data, or state in the methods that the defaults were used."),
            tags$br(), tags$br(),
            HTML("<b>1. The &sigma;<sub>loc</sub> prefactor a.</b> Acquire an immobile sample, i.e. fixed cells, labelled as usual, because the (a) absorbs background and estimator losses in moving environment. Set galvo damping off (as recommended by Vogler 2025), at two or more TCP diameters. Fit MSD(τ) = 4&sigma;² + 4D<sub>w</sub>τ on short lags; a = L·√N̄ / &sigma;. Read the result as (a) <b>lower bound on (a)</b>. Tracker wander is sub-diffusive, so a straight-line fit lands above the true intercept and inflates &sigma;. The true a lies between what you measure and 4. Two diameters also test whether &sigma; ∝ L holds on your system, which is what the <b>beam FWHM &nu;</b> field switches on if it does not."),
            tags$br(), tags$br(),
            HTML("<b>2. The patGeoFactor base.</b> Sequence files store a dimensionless pattern factor. 360 nm reproduces the published Abberior ladders, but confirm one diameter against your vendor before trusting the auto-filled L, since every &sigma; and D<sub>max</sub> scales with it."),
            tags$br(), tags$br(),
            HTML("<b>3. The wander floor D<sub>w</sub>.</b> Tracking on an immobile equivalent gives a slope of that MSD fit. This is the apparent diffusion the tracking loop produces on a stationary target, photon-driven, and it sets the smallest D your system can resolve."),
            tags$br(), tags$br(),
            HTML("<b>For a methods section</b>, quote the acquisition parameters from the sequence (TCP diameter of the reiterated iteration, photon limit, dwell time, damping, pattern, maxOffTime), the values measured from the data (cycle time, median ECO and EFO, time-to-localization, tracks and localizations after cleaning), the three constants above with how they were obtained, and both D and &sigma;<sub>loc</sub> with the route used for each. Vogler et&nbsp;al. 2025 Table&nbsp;2 is a good template."))
        )
      ),

      tags$hr(),

      # ── Time Signature Analysis (TSA) ────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Time Signature Analysis (TSA)"),
          .method_box("green",
            HTML("Time-Signature Analysis for MINFLUX (MFX-TSA), after <b>Vogler et&nbsp;al. 2025</b> (<i>Commun. Biol.</i>; reference implementation <i>Toolbox_for_MINFLUX_SPT / the_road_to_minflux_msd</i>). It separates integrated from dark cycles, recovers the cycle time the instrument actually ran at, and proposes thresholds for splitting tracks at gaps."),
            HTML("An R port of MFX-TSA, covering cycle-time and hardware-parameter recovery (EOD, software overhead, pattern vertices), the integrated/dark cycle balance, per-track cycle progression (cf. <code>print_cycle_stats</code> and <code>plot_cycle_progression_for_ID</code>), and adaptive gap removal from TOT-based thresholds with a 50-localization minimum (cf. <code>suggest_split_thresholds</code> and <code>remove_time_gaps</code>). Every estimator further down rests on this layer, because it builds the linearized cycle-index timeline that makes the irregular Δt of MINFLUX tractable. Each interval is labelled by the number of cycles it spans, round((Δt − SW) / cycle): <b>1 is integrated</b>, meaning the molecule was localized on every cycle, while <b>2 or more</b> means dark cycles came first. A high integrated fraction describes gap-free tracking. A long dark tail describes dropouts, from photophysics, weak signal, or the molecule moving out of range. The whole analysis is independent of the MSD, so it can tell you whether the data support diffusion analysis before any D is fitted."),
            "These are descriptive statistics and carry no significance test. The cycle time is an estimate from the data, so cross-check it against the sequence JSON or the hardware parameters when you have them."),
          checkboxInput("trk_show_tsa", "Enable TSA", FALSE),
          conditionalPanel(
            condition = "input.trk_show_tsa == true",
            tags$b("Sequence file (optional)", style = "font-size:0.9em;"),
            fileInput("trk_seq_json", "Upload MINFLUX sequence (.json)",
                      accept = ".json"),
            uiOutput("trk_seq_json_status"),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Fills in the pattern vertices, the dimensionality, the gap threshold (maxOffTime) and the ",
              "pattern diameter L", tags$sub("TCP"), " behind D", tags$sub("max"), ". ",
              "The cycle time is always measured from the data itself."),
            tags$hr(style = "margin:8px 0;"),
            checkboxInput("trk_dt_log_y", "Log count axis", value = FALSE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "The dark tail sits far below the 1-cycle peak and is hard to read on a linear axis."),
            checkboxInput("trk_dt_mark_cycles", "Mark cycle multiples (1×, 2×, …)", value = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Δt only takes the values SW + k × cycle. The verticals mark each k."),
            tags$hr(style = "margin:8px 0;"),
            numericInput("trk_cycle_dark_max", "Max cycles to display", value = 10, min = 4, step = 1),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Intervals spanning more cycles than this are pooled into a single tail bar."),
            tags$hr(style = "margin:8px 0;"),
            actionButton("trk_tsa_run", "Run TSA", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_tsa == true",
            uiOutput("trk_tsa_quality"),
            tags$hr(style = "margin:8px 0;"),
            h5("Inter-localization interval (Δt)"),
            plotlyOutput("trk_dt_hist_jdd", height = "280px"),
            uiOutput("trk_dt_stats_ui"),
            tags$hr(style = "margin:8px 0;"),
            h5("Integrated / Dark cycle distribution"),
            plotlyOutput("trk_cycle_class_plot", height = "300px"),
            uiOutput("trk_cycle_class_summary"),
            uiOutput("trk_cycle_class_beforeafter"),
            tags$hr(style = "margin:8px 0;"),
            h5("Per-track cycle progression"),
            uiOutput("trk_cycle_prog_tid_ui"),
            plotlyOutput("trk_cycle_prog_plot", height = "260px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "Blue marks localizations captured on every cycle. Red marks localizations that ",
              "follow one or more dark cycles, and a burst of red is the tracker repeatedly ",
              "losing and re-finding the molecule.")
          )
        )
      ),

      tags$hr(),

      # ── Sample geometry (shared) ──────────────────────────────────────────────
      # Defined once and read by every tool that needs it: the MSD's projection
      # factor and free-on-geometry overlay, and the displacement test's Step-2 null.
      # Keeping one definition means the shape can never be set two different ways.
      fluidRow(
        column(4,
          h4("Sample geometry"),
          .method_box("green",
            "The shape the molecules move in or on, and how many directions they are free to move in. Read by every analysis that needs to know what free motion would look like in this sample.",
            HTML("The shape is always <b>three-dimensional</b>, since a specimen has depth whether or not the tracking recorded it, and it is projected into the image plane automatically for 2D datasets. Two places read it: the MSD, both for the factor that converts the measured curve into D and for the optional free-diffusion overlay, and the geometry-aware band of the displacement test. Setting it here once keeps those two from disagreeing."),
            HTML("Both settings matter only when motion is restricted. For molecules moving freely in three dimensions the shape has no effect until they begin reaching its walls. Choose <b>None</b> when there is no bounding structure, such as molecules in solution or a container far larger than anything the molecule explores. \"How these settings enter the calculation\" opposite gives the exact factors and ceilings.")),
          selectInput("trk_geom_type", "Shape",
            choices = c("None (unbounded)"                          = "none",
                        "Sphere / round cell"                       = "sphere",
                        "Rod / spherocylinder"                      = "rod",
                        "Ellipsoid"                                 = "ellipsoid",
                        "Flat membrane / supported bilayer"         = "membrane",
                        "Filament / linear track"                   = "filament"),
            selected = "none"),
          conditionalPanel("input.trk_geom_type == 'sphere'",
            numericInput("trk_disp_R", "Sphere radius (nm)", value = 500, min = 50, step = 50)),
          conditionalPanel("input.trk_geom_type == 'rod'",
            fluidRow(
              column(6, numericInput("trk_disp_r", "Radius (nm)", value = 400, min = 50, step = 25)),
              column(6, numericInput("trk_disp_L", "Cyl. length (nm)", value = 2000, min = 100, step = 100)))),
          conditionalPanel("input.trk_geom_type == 'ellipsoid'",
            fluidRow(
              column(4, numericInput("trk_disp_a", "a (nm)", value = 500, min = 50, step = 50)),
              column(4, numericInput("trk_disp_b", "b (nm)", value = 500, min = 50, step = 50)),
              column(4, numericInput("trk_disp_c", "c (nm)", value = 900, min = 50, step = 50)))),
          conditionalPanel("input.trk_geom_type == 'membrane'",
            numericInput("trk_disp_mtilt", "Tilt from the image plane (°)",
                         value = 0, min = 0, max = 89, step = 5)),
          radioButtons("trk_motion", "Molecule moves",
            choices = c("Through a volume (3 degrees of freedom)"        = "volume",
                        "On a surface (2 degrees of freedom)"            = "surface",
                        "Along a filament or track (1 degree of freedom)" = "filament"),
            selected = "volume"),
          helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
            "This describes how the molecule moves, so a membrane protein has 2 degrees of freedom ",
            "whether the file holds 2 or 3 coordinates. The MSD and the displacement test both read it.")
        ),
        column(8,
          # Reference material, not a running commentary: collapsed by default so the
          # panel shows controls, and the d_eff / ceiling tables stay one click away.
          tags$details(style = "margin-top:6px;",
          tags$summary(style = "cursor:pointer; font-size:0.86em; color:#2c3e50; font-weight:600;",
                       "ℹ How these settings enter the calculation"),
          tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
            "The shape acts in two independent places, on different timescales.",
            tags$br(), tags$br(),

            tags$b("1 · Relating MSD to D"), tags$br(),
            HTML("Applies at every lag, as soon as motion is restricted, and does <b>not</b> depend on the shape's size:"),
            tags$div(style = "margin:3px 0 3px 10px;", tags$code("MSD = 2 · d_eff · D · τ")),
            HTML("With all three coordinates recorded, d<sub>eff</sub> is simply the number of directions the molecule can move in: 3 for a volume, 2 for a surface, 1 for a filament."),
            tags$br(),
            HTML("With only x and y recorded, part of every step points out of the image plane and is not measured. How much survives depends on how the surface is tilted where the molecule sits, averaged over the shape: d<sub>eff</sub> = 1 + ⟨n<sub>z</sub>²⟩, with n<sub>z</sub> the z-component of the surface normal."),
            tags$table(style = "margin:5px 0 5px 10px; border-collapse:collapse;",
              tags$tr(style = "color:#666;",
                tags$th(style = "text-align:left; padding-right:14px;", "surface, 2D data"),
                tags$th(style = "text-align:left; padding-right:14px;", "⟨n_z²⟩"),
                tags$th(style = "text-align:left; padding-right:14px;", "d_eff"),
                tags$th(style = "text-align:left;", "MSD")),
              tags$tr(tags$td("sphere"), tags$td("1/3"), tags$td("4/3"), tags$td(tags$code("8Dτ/3"))),
              tags$tr(tags$td("rod, axis in plane"), tags$td("1/2"), tags$td("3/2"), tags$td(tags$code("3Dτ"))),
              tags$tr(tags$td("flat membrane"), tags$td("1"), tags$td("2"), tags$td(tags$code("4Dτ"))),
              tags$tr(tags$td("filament in plane"), tags$td("-"), tags$td("1"), tags$td(tags$code("2Dτ")))),
            HTML("A membrane tilted by θ gives 1 + cos²θ, between the rod and flat values. Reading restricted motion as free motion under-reports D by a third on a sphere and a quarter on a rod."),
            tags$br(), tags$br(),

            tags$b("2 · The ceiling the MSD approaches"), tags$br(),
            HTML("Applies only once the molecule has had time to cover the shape and/or boundaries of the imaged molecule. e.g. cell shape."),
            tags$table(style = "margin:5px 0 5px 10px; border-collapse:collapse;",
              tags$tr(style = "color:#666;",
                tags$th(style = "text-align:left; padding-right:14px;", "shape"),
                tags$th(style = "text-align:left; padding-right:14px;", "3D"),
                tags$th(style = "text-align:left;", "2D projection")),
              tags$tr(tags$td("sphere, radius R"), tags$td(tags$code("2R²")), tags$td(tags$code("4R²/3"))),
              tags$tr(tags$td("rod, radius r, length L"), tags$td(tags$code("L²/6 + 2r²")), tags$td(tags$code("L²/6 + r²")))),
            HTML("These are squared lengths: size determines where the curve levels off. "),
            tags$br(), tags$br(),
            HTML("Enabling the overlay in the MSD panel fits whichever model is selected to the simulated free curve, so a fitted parameter can be compared against what the shape alone would produce. The effect grows with the fit window: over a short window the molecule covers far less of the shape and the correction is negligible.")))
        )
      ),

      tags$hr(),

      # ── MSD ─────────────────────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("MSD (diffusion)"),
          .method_box("green",
            "Mean squared displacement with optimal least-squares fitting (Michalet 2010; Michalet & Berglund 2012). Cycle-based binning and the D-ceiling follow Vogler et al. 2025.",
            HTML("Per-track displacements are binned by <b>integer cycle count</b>, which is the true elapsed time, so a skipped cycle cannot corrupt the curve. The per-track TAMSDs are then averaged into a track-weighted ensemble (eMSD) on the linearized lag×cycle timeline and fitted by iterative OLSF with Michalet optimal-point selection. Both eDIFF, the ensemble value, and tDIFF, the mean over tracks, are reported. The <b>fit window</b> defaults to 5–50 cycles. Below that, &sigma;<sub>loc</sub> and linker motion dominate; above it, few pairs contribute and the variance grows steeply (Qian 1991), so the fit sits in the diffusive regime in between. The lower bound depends on the probe, so read it off the kink at low lags in your own curve. The upper limit also sets how many lags are computed, capped at 5000, so raise it if you need to reach a confinement plateau."),
            HTML("Optimal-point theory holds for <b>free diffusion on a straight line</b>, so the Confined, Directed and Anomalous models are flagged exploratory. Above the trackable ceiling of roughly 2–3 µm²/s, D comes out too low, and the D<sub>max</sub> panel warns about this. Inside a closed cell the MSD of a freely diffusing molecule bends and saturates on its own, so a bend by itself says little about trapping. Turn on the geometry overlay to see how much of it the shape already accounts for.")),
          checkboxInput("trk_show_msd", "Enable MSD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_msd == true",
            numericInput("trk_msd_min_pts",  "Min localizations per TID (skip short)", value = 50, min = 5, step = 5),
            radioButtons("trk_msd_ndim", "Coordinates used",
              choices  = c("Auto-detect" = "auto", "2D  (x, y)" = "2", "3D  (x, y, z)" = "3"),
              selected = "auto", inline = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              tags$b("Auto"), " picks 3D when z carries real variance (SD above 1 nm) and 2D otherwise. ",
              "The degrees of freedom and the shape come from ", tags$b("Sample geometry"), " above."),
            tags$hr(style = "margin:6px 0;"),
            checkboxInput("trk_show_dmax",     HTML("Diffusion limit analysis (D<sub>max</sub>)"), FALSE),
            checkboxInput("trk_msd_show_hist", "Per-track D distribution", FALSE),
            conditionalPanel(
              condition = "input.trk_msd_show_hist == true || input.trk_msd_cluster == true",
              checkboxInput("trk_msd_hist_log", "Show D as log₁₀ scale", TRUE)
            ),
            checkboxInput("trk_msd_cluster",   "Cluster tracks by D & CI (k-means)", FALSE),
            conditionalPanel(
              condition = "input.trk_msd_cluster == true",
              fluidRow(
                column(7, numericInput("trk_msd_n_clusters", "k (clusters)", value = 3, min = 2, max = 6, step = 1)),
                column(5, checkboxInput("trk_msd_cluster_color_tracks", "Color tracks", FALSE))
              )
            ),
            tags$hr(style = "margin:6px 0;"),
            # The cycle-based MSD (Vogler's the_road_to_minflux_msd workflow) is now
            # THE MSD - always on. The checkbox is kept hidden (constant TRUE) only so the
            # server's isTRUE(input$trk_msd_unified) gates and conditionalPanels keep working.
            div(style = "display:none;", checkboxInput("trk_msd_unified", "MSD", TRUE)),
            div(
              fluidRow(
                column(6, numericInput("trk_uni_min_lag",
                             "Lower cutoff (cycles)", value = 5, min = 1, step = 1)),
                column(6, numericInput("trk_uni_max_lag",
                             "Upper fit limit (cycles)", value = 50, min = 10, max = 5000, step = 5))
              ),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "The lags the fit runs over, in cycles. The upper limit also sets how many lags are ",
                "computed, capped at 5000. Click Plot MSD again to apply a change."),
              checkboxInput("trk_uni_track_weight",
                            "Track-weighted ensemble (eMSD)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Averages the per-track TAMSDs so every track counts once. Turning it off weights ",
                "by pair count, which lets long tracks dominate the long lags."),
              checkboxInput("trk_uni_linear_tau",
                            "Linearized τ axis (lag × cycle)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Places each point at lag × cycle time. Turning it off places it at the true mean Δt ",
                "of the bin instead."),
              selectInput("trk_uni_model", "Fit model",
                          choices = c("Linear (Brownian)"             = "linear",
                                      "Confined (plateau → size L)"   = "confined",
                                      "Directed (drift → velocity v)" = "directed",
                                      "Anomalous (power-law α)"       = "anomalous"),
                          selected = "linear"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Linear returns D and is the only model the OLSF options below apply to. Confined adds ",
                "a corral size L, Directed adds a velocity v, and Anomalous fits Dτ^α. The three ",
                "non-linear models use the whole fit window."),
              conditionalPanel(
                condition = "input.trk_uni_model != 'linear'",
                div(style = "margin-top:4px;padding:5px 8px;border-left:3px solid #e67e22;background:rgba(230,126,34,0.06);font-size:11px;color:#7f8c8d;line-height:1.4;",
                    HTML(paste0(
                      "<b style='color:#e67e22;'>Exploratory.</b> Confined / directed / anomalous MSD ",
                      "combined with MINFLUX irregular cycle-time sampling is beyond the standard, ",
                      "literature-validated free-diffusion workflow (Vogler et&nbsp;al. 2025 build the ",
                      "irregular-Δt machinery but apply it only to free diffusion; their toolbox flags ",
                      "anomalous/confined as “implement on your own”). The model forms themselves are standard ",
                      "(Kusumi et&nbsp;al. 1993 for the corral, Metzler et&nbsp;al. 2014 for the power law), but ",
                      "nobody has validated their fit to short, cycle-skipped MINFLUX tracks, so treat L_c, ",
                      "v and α as exploratory. ",
                      "See docs/MINFLUX_analysis_methods.md §3a.")))),
              checkboxInput("trk_msd_geom",
                          HTML("Overlay <b>free diffusion on the cell geometry</b>"), FALSE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "The green line shows what free diffusion on the ", tags$b("Sample geometry"), " shape ",
              "would produce at the fitted D. A measured curve lying on it is explained by the shape ",
              "alone, and one falling below it carries confinement beyond the shape. Only long fit ",
              "windows show a visible difference."),
            checkboxInput("trk_uni_iterative",
                            "Iterative OLSF (converges on lowest lags)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Repeats the fit, updating the optimal number of points each time, until it settles. ",
                "Turning it off applies the Michalet point count once over the whole window."),
              checkboxInput("trk_uni_drop_negd",
                            "Drop D≤0 tracks from ensemble (QC)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "The eDIFF convention keeps every track, including the ones that fit to D≤0 through ",
                "noise. tDIFF always keeps only the tracks with D>0."),
              checkboxInput("trk_uni_bootstrap",
                            "Bootstrap 95% CI on eDIFF (slower)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Resamples the tracks with replacement about 300 times and refits each time. This ",
                "adds a few seconds."),
              checkboxInput("trk_uni_show_band", "Show IQR band", TRUE),
              checkboxInput("trk_uni_show_ntau", "Show N(τ)", TRUE),
              checkboxInput("trk_uni_show_spag", "Show individual track TAMSDs (spaghetti)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Draws every track as a faint curve, coloured by cluster when clustering is on."),
              checkboxInput("trk_msd_show_bindiag",
                            "Show cycle-binning diagnostic", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Shows the number of pairs at each cycle-lag, and the ratio of real elapsed time to ",
                "lag × cycle, which should come out flat. A wobble means the cycle assignment drifts ",
                "where long dark gaps occur.")
            ),
            tags$br(),
            actionButton("trk_msd_compute", "Plot MSD",
                         icon = icon("chart-line"),
                         class = "btn-primary btn-sm",
                         style = "margin-bottom:4px;"),
            fluidRow(
              column(6, downloadButton("trk_dl_msd_csv", "MSD (CSV)", class = "btn-xs btn-default")),
              column(6, downloadButton("trk_dl_msd_zip", "MSD + Params (ZIP)", class = "btn-xs btn-default"))
            ),
            downloadButton("trk_dl_msd_pertk_csv", "Per-track fits (CSV)", class = "btn-xs btn-default")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_msd == true",
            # MSD sub-panel (cycle-based, always shown when MSD is enabled)
            conditionalPanel(
              condition = "input.trk_msd_unified == true",
              plotlyOutput("trk_msd_unified_plot", height = paste0(TRK_MSD_PLOT_PX, "px")),
              uiOutput("trk_msd_unified_estimates")
            ),
            conditionalPanel(
              condition = "input.trk_msd_unified == true && input.trk_msd_show_bindiag == true",
              tags$hr(style = "margin:8px 0;"),
              h5("Cycle-binning diagnostic"),
              plotlyOutput("trk_msd_bindiag_plot", height = "260px"),
              uiOutput("trk_msd_bindiag_note")
            ),
            # Diffusion limit sub-panel
            conditionalPanel(
              condition = "input.trk_show_dmax == true",
              tags$hr(style = "margin:10px 0;"),
              numericInput("trk_dmax_t_hw", "Hardware overhead tₕₐ (ms)",
                           value = 0.05, min = 0, max = 1, step = 0.01),
              helpText(style = "font-size:0.82em; color:#666;",
                tags$b("L"), " is the diameter of the ", tags$b("final iteration"), ", the one repeated ",
                "continuously while tracking, so it sets how far the particle may move between fixes. ",
                "It is set once under ", tags$b("Acquisition Parameters"), " above, where ",
                sigloc(), " reads the same value; upload the sequence under ",
                tags$b("Time Signature Analysis"), " to fill it in exactly, or type it in yourself. ",
                "Tracking iterations usually land between 100 and 300 nm, and the coarse search ",
                "iterations before lock-on are a different number. ", tags$b("tₕₐ"), " is about ",
                "0.05 ms for a hexagonal TCP. ECO, EFO and the median Δt all come from the data."
              ),
              uiOutput("trk_dmax_panel")
            ),
            # Per-track D histogram sub-panel
            conditionalPanel(
              condition = "input.trk_msd_show_hist == true",
              tags$hr(style = "margin:10px 0;"),
              plotlyOutput("trk_msd_hist_plot", height = "300px"),
              helpText(style = "font-size:0.82em; color:#666;",
                "One bar per track. A narrow peak means the population is homogeneous, while a broad ",
                "or two-humped distribution points to sub-populations.")
            ),
            # Cluster sub-panel
            conditionalPanel(
              condition = "input.trk_msd_cluster == true",
              tags$hr(style = "margin:10px 0;"),
              plotlyOutput("trk_msd_cluster_plot", height = "320px"),
              uiOutput("trk_msd_cluster_model_ui"),
              plotlyOutput("trk_msd_cluster_curves", height = "320px"),
              uiOutput("trk_msd_cluster_models"),
              helpText(style = "font-size:0.82em; color:#666;",
                "log₁₀(D) against the confinement index on a log scale. The dashed line marks CI = 1, ",
                "free Brownian motion; tracks below it are confined and tracks above it are directed. ",
                "Cluster 1 is always the slowest and most confined. Both values are in the per-track CSV.")
            )
          )
        )
      ),

      tags$hr(),

      # ── Confinement (MSD confinement index) ───────────────────────────────────
      # The Rg/Rmax shape test that used to live here was retired 2026-08-03; its UI
      # and server code are in retired_ui.R / retired_tracking_server.R (gitignored,
      # not sourced).
      fluidRow(
        column(4,
          h4("Confinement (MSD index)"),
          .method_box("green",
            "Per-track confinement from the shape of each track's own MSD curve: CI = (MSD_late / MSD_early) / (τ_late / τ_early), compared against a free-diffusion null simulated at the dataset's measured conditions.",
            HTML("A freely diffusing molecule has an MSD that keeps growing in proportion to time, so its late lags exceed its early lags by the same factor as the lag times do. A confined molecule runs out of room, its MSD flattens towards a plateau, and the late rise falls short. A molecule under directed transport overshoots instead. CI captures all of that as one number per track.<br>The comparison is between a track's own short lags and its own long lags, so no diffusion coefficient has to be supplied from outside, and a molecule simply being fast or slow leaves CI where it was. Short-lag and long-lag behaviour carry independent information, which is what makes the ratio meaningful. The MSD it reads is the cycle-binned one, so irregular Δt is already handled.<br><b>The free range is simulated for each dataset.</b> Ideal noise-free free diffusion gives CI&nbsp;=&nbsp;1. Real data sits lower, because localization error adds a constant 2·d·&sigma;<sub>loc</sub>² to every MSD point and inflates the early lags proportionally more than the late ones. How much lower depends on the noise, the mobility, the cycle time, the fit window and the track length, so any fixed threshold would be wrong for most datasets. The range is therefore built by simulating free tracks under the measured &sigma;<sub>loc</sub>, D, cycle time and dark-gap pattern, matched to each track's localization count, and to its own D where the two interact. A track is called <b>confined</b> below the free 5th percentile and <b>directed</b> above the 95th."),
            HTML("<b>There is a hard resolution floor.</b> Confinement becomes visible only when the molecule would otherwise have travelled further than the boundary holding it. The distance a free molecule covers in one track, √(2·d·D·T), is reported alongside the result, and a domain much larger than that stays invisible. A low confined&nbsp;% therefore describes the scales the data can resolve, and says nothing about domains above them. As a rule of thumb the free extent has to exceed the domain size roughly twofold before per-track calls carry information.<br><b>Confined&nbsp;% cannot be compared between datasets with different track lengths or mobilities</b>, because the detection power differs. Longer tracks and faster molecules pick up more confinement for identical biology. Quote the resolution floor next to the percentage.<br><b>The directed tail is insensitive.</b> A modest drift barely lifts the late-versus-early ratio, so a directed&nbsp;% near zero is weak evidence either way. Directionality and Speed&nbsp;&amp;&nbsp;Transport answer that question better.<br><b>Other limits.</b> Each tail holds about 5% of purely free tracks by construction, so small percentages sit at the floor of the method. Tracks too short to fill the MSD fit window are dropped. The test also assumes one behaviour per track, so a molecule that is confined for part of its trajectory and free for the rest returns a blend of the two. Split those with Change-Point or the HMM panel first.")),
          checkboxInput("trk_show_conf", "Enable confinement index", FALSE),
          conditionalPanel(
            condition = "input.trk_show_conf == true",
            helpText(style = "font-size:0.82em; color:#666;",
              "Run Plot MSD first. ", sigloc(), " and D are taken from Data Quality & Calibration."),
            numericInput("trk_conf_ci_nsim", "Null simulations per length", value = 200,
                         min = 50, max = 500, step = 50),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "How many free tracks are simulated at each localization-count grid point. More gives a ",
              "smoother null and takes longer."),
            actionButton("trk_conf_ci_run", "Run confinement index", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$br(), tags$br(),
            downloadButton("trk_dl_conf_ci_csv", "Download CI table (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_conf == true",
            plotlyOutput("trk_conf_ci_plot", height = "380px"),
            uiOutput("trk_conf_ci_summary")
          )
        )
      ),

      tags$hr(),

      # ── Track Segmentation (CPD + GMM) ───────────────────────────────────────
      fluidRow(
        column(4,
          h4("Track Segmentation (Change-Point + GMM)"),
          .method_box("orange",
            "PELT change-point detection (Killick, Fearnhead & Eckley 2012) on the local-D signal, then a Gaussian mixture (BIC-selected) over the per-segment diffusivities.",
            HTML("The signal being segmented is <b>log(r²/Δt)</b>, the local diffusivity, scaled to unit variance (its variance is fixed by the dimensionality) so the BIC penalty is calibrated. The cost therefore accounts for the interval length and survives irregular sampling. Each segment's D is computed as mean(r²/Δt)/(2d), which stays unbiased when r² and Δt vary together, whereas dividing the mean r² by the mean Δt does not."),
            tagList("A diffusivity estimated from one step is noisy, and PELT still cuts a few percent of single-state tracks. The minimum segment length is the guard against that. ",
                    tags$b("The mixture treats every segment as equally precise."), " Short segments scatter more than long ones, so a single population split into segments of very different lengths can return two components. Treat the component count as exploratory and check it against the segment lengths. The ",
                    sigloc(), " noise floor σ²·⟨1/Δt⟩ is subtracted from every segment's D, using the toggle below, so an immobile segment reads close to zero instead of a small positive D.")),
          checkboxInput("trk_show_cpd", "Enable segmentation + GMM", FALSE),
          conditionalPanel(
            condition = "input.trk_show_cpd == true",
            fluidRow(
              column(6, numericInput("trk_cpd_min_seg", "Min segment length (steps)", value = 8L, min = 3L, max = 50L, step = 1L)),
              column(6, numericInput("trk_cpd_max_g",   "Max GMM components",         value = 4L, min = 2L, max = 6L,  step = 1L))
            ),
            checkboxInput("trk_cpd_sub_noise",
                          HTML("Subtract localization-noise floor (&sigma;<sub>loc</sub>) from segment D"), TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Removes the σ²·⟨1/Δt⟩ offset so an immobile segment reads close to zero. Requires ",
              sigloc(), " from Data Quality & Calibration."),
            actionButton("trk_cpd_compute", "Run Segmentation",
                         class = "btn-primary", style = "margin-bottom:8px;"),
            downloadButton("trk_dl_cpd_csv", "Download segments (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_cpd == true",
            plotlyOutput("trk_cpd_hist", height = "320px"),
            uiOutput("trk_cpd_summary")
          )
        )
      ),

      tags$hr(),

      # ── Jump Distance Distribution ────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Jump Distance Distribution (JDD)"),
          .method_box("green",
            "Jump-distance distribution fit with an N-component Rayleigh mixture (Schütz, Schindler & Schmidt 1997; benchmarked in Weimann et al. 2013).",
            HTML("Irregular Δt is handled in one of two ways. Normalizing each jump by its own <b>√Δt</b> lets steps spanning different numbers of cycles pool into one distribution. <b>Integrated steps only</b> keeps single-cycle steps, so every jump covers the same true Δt and the pooling assumption drops away entirely. The fitted spread is corrected for localization error, which widens every jump, and for motion blur, which narrows it because each localization averages the path over the cycle; both corrections use the mean of 1/Δt over the jumps in the fit. Because MINFLUX collects so many localizations per track, single-cycle steps are usually abundant, and the familiar worry about short tracks yielding too few jumps carries little weight here.<br><b>Setting the Δt window.</b> Read the modal Δt off the TSA Δt histogram and put <b>Min Δt</b> just below it. Steps shorter than one acquisition interval carry noise rather than motion and will manufacture a spurious fast component. Set <b>Max Δt</b> to three to five times the mode so untracked gaps stay out. Leaving the &sigma;<sub>loc</sub> field at 0 uses the app-wide value from Data Quality &amp; Calibration, and any other number overrides it.<br><b>Fast acquisition, below about 1 ms.</b> Individual steps become noise-dominated at that rate. <b>Bin-then-median</b> pools consecutive localizations into time bins two to five times the modal Δt wide and takes each bin's median position, which brings slow components back out. <b>Subsample</b> simply takes every N-th localization, which is simpler but averages nothing away.<br><b>Reading D.</b> D appears only when √Δt normalization is on. 1000 nm²/ms equals 1 µm²/s. Membrane proteins usually fall between 0.01 and 0.5 µm²/s and cytoplasmic proteins between 1 and 20 µm²/s."),
            HTML("The √Δt-normalized mode assumes free diffusion and distorts under confinement or anomalous motion, so use Integrated-only when the motion is in question. Keep the component count low, since a Rayleigh mixture treats every state as free Brownian motion. A component whose corrected D sits at or near zero reflects noise rather than a genuine slow population. When a slow population refuses to separate, the MSD handles it better, because it isolates diffusion from noise through the y-intercept rather than through the size of individual steps.")),
          checkboxInput("trk_show_jdd", "Enable JDD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_jdd == true",
            numericInput("trk_jdd_min_pts", "Min points per TID (skip short)",
                         value = 5, min = 2, step = 1),
            numericInput("trk_jdd_n_comp", "Number of diffusion components",
                         value = 2, min = 1, max = 4, step = 1),
            radioButtons("trk_jdd_step_mode", "How to compute steps",
              choices = c(
                "Consecutive localizations" = "consec",
                "Bin-then-median"           = "binmed",
                "Subsample (every N-th)"    = "subsamp"
              ), selected = "consec"),
            conditionalPanel(
              condition = "input.trk_jdd_step_mode == 'binmed'",
              numericInput("trk_jdd_bin_width", "Bin width (ms)",
                           value = 1, min = 0.1, step = 0.1)
            ),
            conditionalPanel(
              condition = "input.trk_jdd_step_mode == 'subsamp'",
              numericInput("trk_jdd_subsamp_n", "Use every N-th localization",
                           value = 4, min = 2, step = 1)
            ),
            checkboxInput("trk_jdd_norm_dt", "Normalize steps by √Δt", TRUE),
            checkboxInput("trk_jdd_integrated_only",
                          "Integrated steps only (single-cycle Δt)", FALSE),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
              "Keeps only the steps with no dark gap, so every jump covers the same true Δt and the ",
              "distribution holds for any motion model. Requires the cycle time from TSA."),
            tags$hr(style = "margin:6px 0;"),
            numericInput("trk_jdd_min_dt", "Min Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 0.1),
            numericInput("trk_jdd_max_dt", "Max Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 1),
            numericInput("trk_jdd_sigma_loc", HTML("localization error &sigma;<sub>loc</sub> (nm, 0 = auto: the app-wide value)"),
                         value = 0, min = 0, step = 0.5),
            actionButton("trk_jdd_run", "Run JDD", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$br(), tags$br(),
            downloadButton("trk_dl_jdd_csv", "Download JDD (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_jdd == true",
            plotlyOutput("trk_jdd_plot", height = "380px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "The bars are the observed step sizes pooled over all tracks, the dashed curves are the ",
              "fitted components and the red curve is their sum. A single peak means one population, ",
              "while a shoulder or a second peak means two."),
            uiOutput("trk_jdd_estimates")
          )
        )
      ),

      tags$hr(),

      # ── Displacement vs Diffusion (Step 1 → Step 2 geometry-aware null) ───────
      fluidRow(
        column(4,
          h4("Displacement vs Diffusion"),
          .method_box("orange",
            "A MINFLUXR original. It compares each track's net start-to-end displacement with the free-diffusion expectation √(4·D·T) evaluated at the ensemble D from CVE.",
            HTML("For 2D Brownian motion ⟨R<sub>net</sub>²⟩ = 4·D·T, and R<sub>net</sub> follows a Rayleigh distribution with scale √(2DT), so each track is judged against that whole spread rather than a single expected value. The displacement is measured in x and y, the image-plane projection, and the &sigma;<sub>loc</sub> noise on the two endpoints is subtracted. Path length appears in the output but takes no part in the classification, since it grows with the sampling rate in the same way path speed does."),
            HTML("The null used here is an <b>infinite plane</b>. A finite or curved cell caps the displacement of a freely diffusing molecule as well, and that is not in this null, so a track landing below the band may be confined or may simply be following the shape of its container. Step 2, free diffusion on a chosen geometry, tells the two apart. Per-track D is far too noisy for this comparison, so the ensemble D is used throughout, which means Data Quality &amp; Calibration has to be run first.<br><b>What the single ensemble D costs.</b> When per-track D varies widely across the population, molecules slower than average drop below the band while diffusing perfectly freely, and the sub-diffusive fraction comes out too high. The panel warns when the spread is large enough for this to bite. Each tail also holds about 5% of free tracks by construction.<br><b>How this relates to the other panels.</b> R<sub>net</sub> uses only the two endpoints and the true elapsed time, so irregular Δt leaves it alone, but it does depend on D. Confinement Analysis tests the shape of the curve and needs no D at all. The two will report different percentages, and this one should be read as an upper bound on confinement.")),
          checkboxInput("trk_show_disp", "Enable displacement vs diffusion", FALSE),
          conditionalPanel(
            condition = "input.trk_show_disp == true",
            numericInput("trk_disp_min_pts", "Min points per TID", value = 10, min = 3, step = 1),
            radioButtons("trk_disp_dims", "Displacement dimensions",
              choices = c("Auto (3D if z present)" = "auto", "2D (x, y)" = "2d", "3D (x, y, z)" = "3d"),
              selected = "auto", inline = TRUE),
            # ONE control for both steps: Step 1 reads it for the motion degrees of
            # freedom, Step 2 for whether to simulate on the boundary or inside it.
            # Shown whenever it can change something: 3D data (Step 1), or a closed
            # geometry in Step 2.
            helpText(style = "font-size:0.82em; color:#666;",
              "Shape and degrees of freedom come from ", tags$b("Sample geometry"), " above. ",
              "The 3D null is a Maxwell distribution, √6DT, and the 2D null a Rayleigh, √4DT. The ",
              "Step-2 geometry null is always simulated in 3D and then projected."),
            selectInput("trk_disp_dsrc", "Diffusion coefficient D (sets the null)",
              choices = c("CVE, short-lag and noise-corrected (recommended)" = "cve",
                          "MSD fit (current MSD model)" = "msd",
                          "Manual value" = "manual"),
              selected = "cve"),
            conditionalPanel("input.trk_disp_dsrc == 'manual'",
              numericInput("trk_disp_dman", "D (µm²/s)", value = 0.1, min = 0, step = 0.01)),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
              "The null is √(2·d·D·T), so ", tags$b("D sets the verdict"), ". Use the short-time ",
              "mobility. CVE measures it before confinement has bent the MSD, while a D taken from an ",
              "MSD fit already carries that bend and will under-call confinement."),
            actionButton("trk_disp_run", "Run displacement", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$hr(style = "margin:8px 0;"),
            checkboxInput("trk_disp_geom",
                          HTML("<b>Step 2:</b> compare to free diffusion on a geometry"), FALSE),
            conditionalPanel(
              condition = "input.trk_disp_geom == true",
              uiOutput("trk_disp_geom_echo"),
              helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
                "Adds a ", tags$b(style = "color:#1B9E5A;", "green band"), " for free diffusion on the ",
                tags$b("Sample geometry"), " shape at your D, matched to your track durations. Tracks ",
                "below the green band are confined beyond what the shape explains, and tracks between ",
                "green and grey are held back by the shape while diffusing freely. With the shape set to ",
                tags$b("None"), " the band collapses onto the Step-1 unbounded null.")
            ),
            tags$br(),
            downloadButton("trk_dl_disp_csv", "Download displacement table (CSV)"),
            helpText(style = "font-size:0.8em; color:#666; margin-top:4px;",
              "One row per track, holding R", tags$sub("net"), ", R", tags$sub("max"),
              ", path length, duration, ρ and class, together with the ensemble D and ", sigloc(), ".")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_disp == true",
            plotlyOutput("trk_disp_plot", height = "380px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "One dot per track, duration against R", tags$sub("net"), ". The grey band spans the 5th to ",
              "95th percentile of free diffusion at the measured D. ",
              tags$b(style = "color:#2c7fb8;", "Blue"), " tracks fall below it and read as sub-diffusive, ",
              tags$b(style = "color:#9e9e9e;", "grey"), " tracks lie within it, and ",
              tags$b(style = "color:#e31a1c;", "red"), " tracks rise above it and read as directed. ",
              "About 5% of free tracks land in each tail by construction."),
            uiOutput("trk_disp_summary")
          )
        )
      ),

      tags$hr(),

      # ── Directionality / autocorrelation ─────────────────────────────────────
      fluidRow(
        column(4,
          h4("Directionality / Autocorrelation"),
          .method_box("green",
            "Turn-angle distribution and velocity autocorrelation, both read against an isotropic-Brownian null (Burov et al. 2013, distribution of directional change).",
            HTML("The turn-angle null holds for any Δt, because the angle between two isotropic steps is uniform whatever their durations. The velocity autocorrelation is binned by <b>elapsed time</b>, cycle lag × cycle time, so a skipped cycle lands in the right bin. Both use steps that span at most the allowed number of dark cycles. The small &sigma;<sub>loc</sub> of MINFLUX removes most of the anti-persistence artefact that noise creates at the first lag in camera-based data."),
            tagList("A single short track carries very little power, so read these curves pooled across tracks. A small negative value at the first lag is what ",
                    sigloc(), " produces for an ordinary diffusing molecule and appears in free data too.")),
          checkboxInput("trk_show_dir", "Enable directionality analysis", FALSE),
          conditionalPanel(
            condition = "input.trk_show_dir == true",
            numericInput("trk_dir_min_pts", "Min points per TID", value = 6, min = 4, step = 1),
            numericInput("trk_dir_max_dark", "Max dark cycles per step", value = 2, min = 0, step = 1),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Steps spanning more dark cycles than this are left out. One or two is usually safe, ",
              "while longer gaps can hide a merge with another molecule."),
            actionButton("trk_dir_run", "Run directionality", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$br(), tags$br(),
            downloadButton("trk_dl_dir_csv", "Download angle data (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_dir == true",
            plotlyOutput("trk_dir_angle_plot", height = "300px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "A peak at 0° means directed motion, a peak at 180° means confined or tethered motion, ",
              "and a flat distribution means Brownian motion."),
            plotlyOutput("trk_dir_autocorr_plot", height = "300px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "C(τ) against real elapsed time. A curve that stays positive and decays slowly means ",
              "directed motion. A small dip at the first lag followed by roughly zero means free ",
              "Brownian motion, and that dip is the ", sigloc(), " signature, C ≈ −σ²/MSD₁."),
            uiOutput("trk_dir_summary")
          )
        )
      ),

      tags$hr(),

      # ── HMM ──────────────────────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Diffusion State Segmentation (HMM)"),
          .method_box("amber",
            "Continuous-time hidden Markov model fitted by maximum likelihood over all tracks at once. Hidden Markov models for diffusive states: Das, Cairo & Coombs 2009; pooling short tracks: Persson et al. 2013 (vbSPT); a Markov jump process observed at irregular times: Jackson 2011 (msm).",
            HTML("Each step's transition matrix is <b>exp(Q·Δt)</b>, so steps across dark cycles are handled exactly and the <b>switching rates are fitted parameters</b>. The step length is modelled as diffusion plus motion blur plus the noise of both localizations, using each localization's own photon-model ",
                 as.character(sigloc()), " when that model is active, so the D values are noise-corrected. The decoded state path colours the tracks and feeds the per-state MSD; rates never come from it."),
            tagList("The model gives each step one state, so a switch inside a step is not represented, and it treats consecutive steps as independent although they share a localization. On simulated MINFLUX acquisitions, D and dwell times are recovered within about 20% when mean dwells last 10 cycles or more; shorter dwells mix the states and read too long. States outside that range, or with fewer switches than any validated case, are flagged. ",
                    "The likelihood intervals assume one homogeneous population and come out too narrow on real molecules; the bootstrap over tracks gives the intervals to report. ",
                    "The emission reads x and y only. Choose the number of states yourself; the summary shows how BIC ranks one state fewer and one more.")),
          checkboxInput("trk_show_hmm", "Fit HMM to tracks", FALSE),
          conditionalPanel(
            condition = "input.trk_show_hmm == true",
            numericInput("trk_hmm_n_states",    "Number of states",                value = 2,   min = 2, max = 4,   step = 1),
            numericInput("trk_hmm_min_pts",     "Min localizations per TID",       value = 10,  min = 4, step = 1),
            numericInput("trk_hmm_n_ex",        "Example tracks to display",       value = 6,   min = 1, max = 30,  step = 1),
            numericInput("trk_hmm_boot",        "Bootstrap refits over tracks (0 = off)", value = 0, min = 0, max = 200, step = 10),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Each refit takes about as long as the main fit. Use 50 or more for intervals you intend to report."),
            actionButton("trk_hmm_run", "Run HMM", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$br(), tags$br(),
            downloadButton("trk_dl_hmm_csv", "Download state assignments (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_hmm == true",
            plotlyOutput("trk_hmm_params_plot",  height = "280px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "One bar per state, ordered slow to fast, with its 95% interval. The percentage is the fraction of time spent in that state."),
            plotlyOutput("trk_hmm_trans_plot",   height = "260px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "The probability of being in the column state one cycle after being in the row state, from the fitted rates. ",
              "A strong diagonal means long dwells."),
            plotlyOutput("trk_hmm_example_plot", height = "400px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "The longest tracks, each centred on its own centroid and coloured by its most likely state under the fitted model."),
            uiOutput("trk_hmm_summary")
          )
        )
      ),

      tags$hr(),

      # ── Per-segment MSD by HMM state ──────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Per-Segment MSD (by HMM state)"),
          .method_box("orange",
            "Pools the contiguous same-state runs of each HMM state into a per-state ensemble MSD and fits it with the chosen MSD model. This is where mixed-motion tracks pay off.",
            "The cycle-index binning and the fit window are inherited from the MSD panel. Fitting each state separately lets the Confined model work on the slow, wall-associated state and return a confinement size L, while the Directed model works on the fast, mobile state and returns a velocity v. Fitting the whole track at once would average the two phases together.",
            "The selection is circular: states are assigned from the step sizes, and the MSD is then computed from those same steps, cut into runs at every state switch. Short, truncated runs flatten the MSD, so the slow state can show confinement that the selection produced. Confined and anomalous MSD models on MINFLUX data are exploratory, as the MSD panel notes. Everything here rests on the state assignment, so a poorly fitted HMM leaves the per-state MSD without meaning."),
          checkboxInput("trk_show_seg_msd", "Enable per-segment MSD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_seg_msd == true",
            numericInput("trk_seg_msd_min_pts", "Min points per segment", value = 15, min = 5, step = 5),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Run the HMM above first. The model and fit window come from the MSD panel."),
            actionButton("trk_seg_msd_run", "Run per-segment MSD", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$br(), tags$br(),
            downloadButton("trk_dl_seg_msd_csv", "Per-state MSD (CSV)", class = "btn-xs btn-default")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_seg_msd == true",
            plotlyOutput("trk_seg_msd_plot", height = "340px"),
            uiOutput("trk_seg_msd_summary"),
            helpText(style = "font-size:0.82em; color:#666;",
              "One eMSD per HMM state, ordered slow to fast, with each state's model fit dashed. ",
              "The table sets D from the MSD beside D from the HMM so the two can be compared.")
          )
        )
      ),

      tags$hr(),

      # ── Speed & Transport (exploratory - kept late in the tab) ────────────────
      fluidRow(
        column(4,
          h4("Speed & Transport"),
          .method_box("amber",
            "A MINFLUXR original, combining per-step path speed, a drift velocity averaged over a time window, and a noise-corrected local anomalous exponent α(t). It finds sustained fast segments and scores each one against evidence drawn from outside the speed signal.",
            tagList("Everything is computed on the irregular-Δt data as recorded. Speeds use the real Δt of each step, the drift velocity averages net displacement over a window W and comes out in µm/s, and the local α is de-biased using ",
                    sigloc(), "."),
            HTML("<b>Path speed cannot separate diffusion from transport.</b> Diffusion inflates the instantaneous speed on its own, so a high path speed leaves the question open. Segments are flagged on drift velocity instead, and since the flag is itself a speed threshold, the velocity of a flagged segment can only repeat what the flag already said. A tracker that has lost its molecule and gone searching leaves the same coarse signature, fast and outward at the end of a track. The corroboration plots exist to separate the two: a real event holds a high integrated fraction, whereas a search run accumulates dark cycles. The <b>edge guard</b> drops the first and last window of every track, because lock-on jitter otherwise piles false events at position 0. <b>Only a matched negative control put through the same processing settles artifact against signal.</b>")),
          checkboxInput("trk_show_speed", "Enable speed & transport", FALSE),
          conditionalPanel(
            condition = "input.trk_show_speed == true",

            tags$b("Finding drift segments", style = "font-size:0.9em;"),
            numericInput("trk_speed_drift_win", "Drift window W (ms)",
                         value = 10, min = 1, step = 2),
            numericInput("trk_speed_fold", "Threshold (× diffusive baseline)",
                         value = 3, min = 1, step = 0.5),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "A segment is called when the drift velocity stays above the threshold times the diffusive ",
              "baseline for the whole window. A larger W rejects more noise and blurs the onset."),
            numericInput("trk_speed_min_pts", "Min points per TID", value = 20, min = 5, step = 5),

            tags$details(style = "margin:8px 0;",
              tags$summary(tags$b("Advanced", style = "font-size:0.88em;")),
              tags$div(style = "padding-left:4px;",
                numericInput("trk_speed_drift_minpts", "Min localizations in window",
                             value = 5, min = 2, step = 1),
                numericInput("trk_speed_min_run", "Min sustained windows to call a segment",
                             value = 3, min = 1, step = 1),
                numericInput("trk_speed_edge_guard", "Edge guard (× W), ignoring segments this close to a track start or end",
                             value = 1, min = 0, max = 5, step = 0.5),
                numericInput("trk_speed_onset_thr", "Absolute threshold (µm/s, 0 = use × baseline)",
                             value = 0, min = 0, step = 0.1),
                numericInput("trk_speed_sigma_nm",
                             HTML("Localization precision &sigma; (nm, 0 = auto)"),
                             value = 0, min = 0, step = 1),
                helpText(style = "font-size:0.78em; color:#666; margin-top:-4px;",
                  "Sets the √2σ/W noise floor and the α de-bias. Auto takes the app-wide ", sigloc(),
                  " from Data Quality & Calibration."),
                tags$hr(style = "margin:6px 0;"),
                tags$b("Local α(t)", style = "font-size:0.86em;"),
                numericInput("trk_speed_alpha_win", "α window (steps)", value = 21, min = 7, max = 201, step = 2),
                numericInput("trk_speed_alpha_maxlag", "α max lag (steps)",
                             value = 6, min = 2, max = 30, step = 1),
                numericInput("trk_speed_smooth", "Rolling median window (steps)",
                             value = 5, min = 1, max = 51, step = 2))),

            actionButton("trk_speed_run", "Run speed & transport", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$hr(style = "margin:8px 0;"),
            selectInput("trk_export_hist_metric", "Distribution metric (left plot)",
                        choices = c("Peak drift velocity"          = "v_drift_peak",
                                    "Segment velocity"             = "export_v",
                                    "Fraction of track directed"   = "directed_frac",
                                    "Fraction immobile/confined"   = "immobile_frac",
                                    "Transport episodes / track"   = "n_events",
                                    "Max local α"                  = "alpha_max",
                                    "Event integrated fraction"    = "event_integ_frac",
                                    "Event photon ratio"           = "event_eco_ratio",
                                    "Track tail integrated frac"   = "tail_integ_frac"),
                        selected = "v_drift_peak"),
            downloadButton("trk_dl_export_csv", "Per-track summary (CSV)", class = "btn-xs btn-primary")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_speed == true",
            tags$div(class = "well well-sm", style = "margin-bottom:6px;",
                     uiOutput("trk_export_summary_txt")),
            plotlyOutput("trk_export_corrob", height = "330px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "One dot per candidate segment. The x axis is where along the track the segment starts, ",
              "from 0 to 1, and the y axis is its integrated fraction. Segments high on the plot were ",
              "localized on nearly every cycle throughout, which is what a real event looks like. ",
              "Segments low on the plot accumulated dark cycles while the tracker searched. ",
              tags$b("Read the evidence off the height.")),
            fluidRow(
              column(6, plotlyOutput("trk_export_hist",       height = "290px")),
              column(6, plotlyOutput("trk_export_onset_hist", height = "290px"))
            ),
            helpText(style = "font-size:0.82em; color:#666;",
              tags$b("Left"), " shows the chosen per-track metric, grey for unflagged tracks and green ",
              "for flagged ones. A green mode standing clear of the grey one suggests a distinct ",
              "subpopulation. ", tags$b("Right"), " shows where along each track the motion is strongest. ",
              "A pile-up near 1 occurs in real and artifactual data alike."),
            tags$details(style = "margin-top:8px;",
              tags$summary(tags$b("Inspect a single track")),
              uiOutput("trk_speed_tid_ui"),
              tags$div(style = "margin-bottom:4px;",
                checkboxInput("trk_speed_show_drift", "Drift (transport) velocity", TRUE),
                checkboxInput("trk_speed_show_alpha", "Local exponent α(t)", TRUE),
                checkboxInput("trk_speed_show_path",  "Path speed (reference)", TRUE)),
              plotlyOutput("trk_speed_plot", height = "340px"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:6px;",
                tags$b(style = "color:#1B9E5A;", "Green"), " is the drift velocity, with the noise floor ",
                "dashed in orange and the onset threshold dashed in red. ",
                tags$b(style = "color:#7B4FBE;", "Purple"), " is the local α, which reaches 2 under ",
                "directed motion and sits near 1 for Brownian motion. ",
                tags$b(style = "color:#888;", "Grey dotted"), " is the path speed, which noise dominates."),
              downloadButton("trk_dl_speed_csv", "Per-timepoint curves (CSV)", class = "btn-xs btn-default")
            )
          )
        )
      ),

      tags$hr(),

      # ── Spatially Resolved D-map ──────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Spatially Resolved D-map"),
          .method_box("green",
            HTML("Single-lag diffusivity per spatial cell, D = [mean(r²/4Δt) − σ²·mean(1/Δt)] / [1 − 2R·t<sub>c</sub>·mean(1/Δt)], which is the short-lag limit of the CVE (Vestergaard 2014, Vogler 2025). The denominator corrects for motion blur."),
            tagList("Every step is placed using its real Δt. The 'integrated only' option drops steps that span dark cycles, since the midpoint of such a step says little about where the molecule was. The apparent-diffusion floor ",
                    sigloc(), "²/Δt is subtracted once per cell so that immobile regions read as immobile. The cell value is a mean because the relation between r² and D holds for the mean only; a median of single steps reads about 30% low."),
            "A single-lag D is an instantaneous estimate and carries no information about confinement or anomalous motion. A mean is sensitive to a mis-linked jump, so check the per-cell standard error in the hover and the CSV. Cells with too few steps are left blank rather than shown unreliably. Uneven sampling across the field of view can imprint its own pattern on the map."),
          checkboxInput("trk_show_dmap", "Compute D-map", FALSE),
          conditionalPanel(
            condition = "input.trk_show_dmap == true",
            numericInput("trk_dmap_bins",    "Grid resolution (N×N cells)", value = 10, min = 3, max = 500, step = 1),
            numericInput("trk_dmap_min_pts", "Min steps per cell",               value = 5,  min = 3, step = 1),
            checkboxInput("trk_dmap_show_n", "Overlay step count per cell",      FALSE),
            checkboxInput("trk_dmap_integrated_only",
                          "Integrated steps only (drop dark-cycle steps)", TRUE),
            checkboxInput("trk_dmap_sub_noise",
                          HTML("Subtract localization-noise floor (&sigma;<sub>loc</sub>)"), TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Each cell holds the mean of r²/4Δt over the steps whose midpoint falls inside it. ",
              "Noise subtraction requires ", sigloc(), " from Data Quality & Calibration."),
            actionButton("trk_dmap_run", "Run D-map", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$br(), tags$br(),
            downloadButton("trk_dl_dmap_csv", "Download D-map table (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_dmap == true",
            plotlyOutput("trk_dmap_plot", height = "520px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "Cell colour is the local D in nm²/ms, cool for slow and warm for fast. ",
              "Blank cells held too few steps to estimate."),
            uiOutput("trk_dmap_summary")
          )
        )
      ),

      tags$br(),

      # ── Methods & References (consolidated citation box) ──────────────────────
      # Provenance for every analysis tool. Full evaluation of appropriateness for
      # MINFLUX data (irregular Δt, short tracks) is in docs/MINFLUX_analysis_methods.md.
      tags$details(
        style = "margin:8px 0 16px; font-size:0.85em; color:#444;",
        tags$summary(style = "cursor:pointer; font-weight:bold; color:#2c3e50;",
                     "Methods & References"),
        tags$div(style = "padding:8px 4px 4px 12px;",
          tags$p(style = "color:#666; margin:0 0 6px;",
            "Primary reference for each tool. See ",
            tags$code("docs/MINFLUX_analysis_methods.md"),
            " for the full bibliography and a MINFLUX-appropriateness assessment."),
          tags$table(style = "border-collapse:collapse; line-height:1.5;",
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Time-Signature Analysis")),
                    tags$td("MFX-TSA, Vogler et al. 2025 (Commun. Biol.); reference implementation Toolbox_for_MINFLUX_SPT (the_road_to_minflux_msd). Integrated and dark cycles, cycle time, adaptive gap removal.")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("MSD / OLSF")),
                    tags$td("Michalet 2010; Michalet & Berglund 2012 (optimal least-squares fitting); cycle-based binning & workflow per Vogler et al. 2025 (Toolbox_for_MINFLUX_SPT)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b(sigloc(), " (photon budget)")),
                    tags$td(HTML("MINFLUX CRLB for a doughnut TCP, &sigma; = L<sub>TCP</sub>/(a&radic;N) with N = ECO and a an instrument constant; L from the sequence's patGeoFactor. Balzarotti et al. 2017; Vogler et al. 2025 for the tracking-iteration geometry"))),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("CVE (D)")),
                    tags$td("Vestergaard, Blainey & Flyvbjerg 2014, restricted here to consecutive integrated (single-cycle) pairs")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("JDD")),
                    tags$td("Schütz, Schindler & Schmidt 1997; benchmarked in Weimann et al. 2013")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("HMM states")),
                    tags$td("Das, Cairo & Coombs 2009; vbSPT (Persson et al. 2013); continuous-time rates at irregular sampling: Jackson 2011")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Confinement")),
                    tags$td("τ-based MSD confinement index against a simulated free null (cf. Simson, Sheets & Jacobson 1995)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Directionality")),
                    tags$td("Burov et al. 2013 (distribution of directional change)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Change-point")),
                    tags$td("Killick, Fearnhead & Eckley 2012 (PELT), on log(r²/Δt)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b(HTML("D<sub>max</sub> ceiling"))),
                    tags$td("Vogler et al. 2025 (Commun. Biol.); SimuFLUX (Marin & Ries 2025)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Anomalous α")),
                    tags$td("Metzler et al. 2014; short-track caveats per Kepten et al. 2015 (ensemble estimate)"))
          ),
          tags$p(style = "color:#888; margin:8px 0 0; font-size:0.95em;",
            "MINFLUX foundational: Balzarotti et al. 2017 (Science); Gwosch et al. 2020; ",
            "Schmidt et al. 2021; Wolff et al. 2023; Deguchi et al. 2023. ",
            "Motion blur: Michalet & Berglund 2012. The CVE bounds R from the cycle timing and the photon collection window; the MSD fits use R = 1/6.2, the MinSpt toolbox value (Vogler et al. 2025). No direct measurement of R for MINFLUX has been published.")
        )
      )
    )  # end tabPanel("Analysis")

  )  # end tabsetPanel
)
