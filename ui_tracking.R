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
              "Blue (slow) to red (fast) by ", tags$b("drift velocity"), " over a centred window W; grey: ",
              "window too short. Keep W well below a typical track duration. Only sustained red well above ",
              "\u221a(4D/W) is directed transport.")
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
            tagList(tags$b("D"), " from the covariance-based estimator, CVE (Vestergaard, Blainey & ",
                    "Flyvbjerg 2014). ", sigloc(), " from the MINFLUX photon budget, L", tags$sub("TCP"),
                    "/(a\u221AN) per localization (Vogler et\u00a0al. 2025)."),
            HTML("The CVE uses only consecutive single-cycle steps, so dark-cycle intervals stay out of the covariance, and D is unaffected by motion blur. &sigma;<sub>loc</sub> is computed per localization from its photon count, so the way the molecule moves cannot bias it. With a = 4, the ideal value, &sigma;<sub>loc</sub> is a lower bound."),
            HTML("The tracker wanders even on an immobile target, which puts an instrumental floor under D: measure it and enter it as the wander floor D<sub>w</sub>. The CVE &sigma; shown alongside is a cross-check, not a measurement.")),
          actionButton("trk_cve_run", HTML("Estimate &sigma;<sub>loc</sub> &amp; quality"),
                       icon = icon("ruler-combined"),
                       class = "btn-primary", style = "margin-bottom:6px; width:100%;"),
          helpText(style = "font-size:0.82em; color:#666; margin-top:0;",
            "Run this first: the panels below use its ", sigloc(), " and cycle time."),
          # ── Acquisition parameters (optional cycle-time model) ──────────────
          tags$details(class = "sim-group",
          tags$summary(tags$b("Acquisition Parameters")),
          checkboxInput("trk_hw_enable",
                        "Use hardware model for cycle time (optional)", FALSE),
          checkboxInput("trk_cycle_use_regression",
                        "Fit cycle time from data (EOD/SW regression)", FALSE),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Cycle time and SW as the slope and intercept through the \u0394t clusters. Off: the cluster spacing."),
          numericInput("trk_msd_cycle_ms",
                       "Cycle time manual override (ms, 0 = use auto)",
                       value = 0, min = 0, step = 0.01),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "> 0 overrides every automatic estimate."),
          uiOutput("trk_cycle_time_status"),
          numericInput("trk_L_tcp",
                       HTML("TCP diameter L<sub>TCP</sub>, final iteration (nm)"),
                       value = DMAX_L_TCP_DEFAULT, min = 10, max = 2000, step = 10),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "The iteration repeated during tracking; sets ", sigloc(), " and D", tags$sub("max"),
            ". A sequence file loaded under TSA fills it in."),
          numericInput("trk_crlb_pref",
                       HTML("&sigma;<sub>loc</sub> prefactor a, in &sigma; = L/(a&middot;&radic;N)"),
                       value = TCP_CRLB_PREFACTOR, min = 0.5, max = 4, step = 0.1),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "4 is the ideal value, so ", sigloc(), " is a lower bound. A value measured on an ",
            "immobile sample brackets it from the other side."),
          numericInput("trk_geofactor_base",
                       HTML("patGeoFactor base (nm), L = patGeoFactor &times; base"),
                       value = TCP_L_PER_GEOFACTOR_NM, min = 50, max = 2000, step = 10),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Converts the sequence's patGeoFactor into a diameter. Check one diameter against your ",
            "vendor's figure."),
          numericInput("trk_beam_fwhm",
                       HTML("Excitation doughnut FWHM &nu; (nm)"),
                       value = TCP_BEAM_FWHM_NM, min = 0, max = 1000, step = 5),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "0 = off (", sigloc(), " \u221d L). Set it only if ", sigloc(), " rises faster than L on ",
            "your system."),
          numericInput("trk_dw_floor",
                       HTML("Tracker wander floor D<sub>w</sub> (nm²/ms, 0 = off)"),
                       value = TRK_DW_FLOOR_DEFAULT, min = 0, step = 10),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Apparent D of an immobile target, measured at the same TCP diameter and photon budget. ",
            "Reported beside D, never subtracted."),
          numericInput("trk_sigma_loc_nm",
                       HTML("Localization precision &sigma;<sub>loc</sub> (nm, 0 = auto)"),
                       value = 0, min = 0, step = 0.5),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "0 = automatic (the photon budget). > 0 overrides it in every panel."),
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
              "Cycle time = dwell + EOD overhead \u00d7 vertices. When enabled it replaces the ",
              "data-driven estimate.")
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
            "One CVE ", sigloc(), " per long-enough track, showing how steps vary between tracks. ",
            "It does not feed the app-wide ", sigloc(), "."),
          .detail(label = "About these two settings",
            HTML("Shorter tracks give a covariance too noisy to use. The first localizations of a track are taken while the pattern centres on the molecule and scatter more; skipping them removes that transient.")),
        ),
        column(8,
          uiOutput("trk_cve_status"),
          plotlyOutput("trk_sigma_pertrack_plot", height = "220px"),
          tags$details(style = "margin-top:6px;",
            tags$summary(style = "cursor:pointer; font-size:0.86em; color:#555;",
                         "Data-selection funnel (what's excluded, for your methods)"),
            tags$div(style = "padding:4px 0 0 4px;", uiOutput("trk_quality_funnel"))),
          tags$details(style = "margin-top:6px;",
            tags$summary(style = "cursor:pointer; font-size:0.86em; color:#555;",
                         "Does the cleaning select on mobility? (D by track continuity)"),
            tags$div(style = "padding:4px 0 0 4px;",
              helpText(style = "font-size:0.82em; color:#666; margin:0 0 4px 0;",
                "A molecule that moves fast relative to the beam pattern is lost for a cycle more often, ",
                "so its track collects dark cycles. The table gives D for each quarter of the analysed ",
                "tracks, ranked by continuity."),
              uiOutput("trk_continuity_mobility"))),
          .detail(label = "\u2139 Setup for your MINFLUX system (what to measure)",
            HTML("Three values under <b>Acquisition Parameters</b> belong to your microscope, and the defaults may not fit it:"),
            tags$ul(style = "margin:4px 0; padding-left:18px;",
              tags$li(HTML("<b>&sigma;<sub>loc</sub> prefactor a</b>: track an immobile sample at two or more TCP diameters.")),
              tags$li(HTML("<b>patGeoFactor base</b>: check one diameter against your vendor.")),
              tags$li(HTML("<b>Wander floor D<sub>w</sub></b>: the MSD slope of the same immobile sample."))),
            HTML("The user guide gives the procedure and what to report in a methods section."))
        )
      ),

      tags$hr(),

      # ── Time Signature Analysis (TSA) ────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Time Signature Analysis (TSA)"),
          .method_box("green",
            HTML("Time-signature analysis after <b>Vogler et&nbsp;al. 2025</b> (MinSpt)."),
            HTML("Each interval is labelled by the cycles it spans, round((\u0394t \u2212 SW) / cycle): <b>1 is integrated</b>, 2 or more follows dark cycles. It recovers the cycle time the instrument ran at and builds the cycle-index timeline the other panels use."),
            "Descriptive, with no significance test. The cycle time is estimated from the data; cross-check it against the sequence file."),
          checkboxInput("trk_show_tsa", "Enable TSA", FALSE),
          conditionalPanel(
            condition = "input.trk_show_tsa == true",
            tags$b("Sequence file (optional)", style = "font-size:0.9em;"),
            fileInput("trk_seq_json", "Upload MINFLUX sequence (.json)",
                      accept = ".json"),
            uiOutput("trk_seq_json_status"),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Fills in the pattern, dimensionality, gap threshold and L", tags$sub("TCP"), ". ",
              "The cycle time is always measured from the data."),
            tags$hr(style = "margin:8px 0;"),
            checkboxInput("trk_dt_log_y", "Log count axis", value = FALSE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "The dark tail is small next to the 1-cycle peak."),
            checkboxInput("trk_dt_mark_cycles", "Mark cycle multiples (1×, 2×, …)", value = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "\u0394t only takes the values SW + k \u00d7 cycle."),
            tags$hr(style = "margin:8px 0;"),
            numericInput("trk_cycle_dark_max", "Max cycles to display", value = 10, min = 4, step = 1),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Longer intervals are pooled into one tail bar."),
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
              "Blue: localized on consecutive cycles. Red: after one or more dark cycles.")
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
            "The shape the molecules move in or on, and how many directions they can move in.",
            HTML("The shape is three-dimensional and is projected into the image plane for 2D data. The MSD (its conversion to D and its free-diffusion overlay) and the displacement test both read it."),
            HTML("It matters only when motion is restricted. Choose <b>None</b> when nothing bounds the molecule.")),
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
            "How the molecule moves, not how many coordinates the file holds: a membrane protein has 2.")
        ),
        column(8,
          # Reference material, not a running commentary: collapsed by default so the
          # panel shows controls, and the d_eff / ceiling tables stay one click away.
          tags$details(style = "margin-top:6px;",
          tags$summary(style = "cursor:pointer; font-size:0.86em; color:#2c3e50; font-weight:600;",
                       "ℹ How these settings enter the calculation"),
          tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
            "The shape acts in two places, on different timescales.",
            tags$br(), tags$br(),

            tags$b("1 · Relating MSD to D"), tags$br(),
            HTML("At every lag, independent of the shape's size:"),
            tags$div(style = "margin:3px 0 3px 10px;", tags$code("MSD = 2 · d_eff · D · τ")),
            HTML("With x, y and z recorded, d<sub>eff</sub> is the number of directions the molecule can move in. With x and y only, part of each step leaves the image plane: d<sub>eff</sub> = 1 + \u27e8n<sub>z</sub>\u00b2\u27e9, with n<sub>z</sub> the z-component of the surface normal."),
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
            HTML("A membrane tilted by \u03b8 gives 1 + cos\u00b2\u03b8. Treating surface motion as free under-reports D by a third on a sphere and a quarter on a rod."),
            tags$br(), tags$br(),

            tags$b("2 · The ceiling the MSD approaches"), tags$br(),
            HTML("Only once the molecule has covered the shape:"),
            tags$table(style = "margin:5px 0 5px 10px; border-collapse:collapse;",
              tags$tr(style = "color:#666;",
                tags$th(style = "text-align:left; padding-right:14px;", "shape"),
                tags$th(style = "text-align:left; padding-right:14px;", "3D"),
                tags$th(style = "text-align:left;", "2D projection")),
              tags$tr(tags$td("sphere, radius R"), tags$td(tags$code("2R²")), tags$td(tags$code("4R²/3"))),
              tags$tr(tags$td("rod, radius r, length L"), tags$td(tags$code("L²/6 + 2r²")), tags$td(tags$code("L²/6 + r²")))),
            HTML("Size sets where the curve levels off. The MSD panel's overlay fits the selected model to the simulated free curve for comparison; the effect grows with the fit window.")))
        )
      ),

      tags$hr(),

      # ── MSD ─────────────────────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("MSD (diffusion)"),
          .method_box("green",
            "Mean squared displacement with optimal least-squares fitting (Michalet 2010; Michalet & Berglund 2012), with cycle-based binning after Vogler et al. 2025.",
            HTML("Displacements are binned by <b>integer cycle count</b>, the true elapsed time, so a skipped cycle cannot corrupt the curve. The fit window defaults to 5\u201350 cycles, between the noise-dominated short lags and the sparse long ones; read the lower bound off the kink in your own curve."),
            HTML("Optimal-point theory holds for free diffusion, so the non-linear models are exploratory. Inside a closed cell a free molecule's MSD bends on its own: use the geometry overlay before reading a bend as trapping.")),
          checkboxInput("trk_show_msd", "Enable MSD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_msd == true",
            numericInput("trk_msd_min_pts",  "Min localizations per TID (skip short)", value = 50, min = 5, step = 5),
            radioButtons("trk_msd_ndim", "Coordinates used",
              choices  = c("Auto-detect" = "auto", "2D  (x, y)" = "2", "3D  (x, y, z)" = "3"),
              selected = "auto", inline = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              tags$b("Auto"), ": 3D when z varies (SD above 1 nm), else 2D."),
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
                "The fit window in cycles. The upper limit also sets how many lags are computed (max 5000)."),
              checkboxInput("trk_uni_track_weight",
                            "Track-weighted ensemble (eMSD)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Every track counts once. Off: weighted by pair count, so long tracks dominate the long lags."),
              checkboxInput("trk_uni_linear_tau",
                            "Linearized τ axis (lag × cycle)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Off: each point sits at the mean \u0394t of its bin."),
              selectInput("trk_uni_model", "Fit model",
                          choices = c("Linear (Brownian)"             = "linear",
                                      "Confined (plateau → size L)"   = "confined",
                                      "Directed (drift → velocity v)" = "directed",
                                      "Anomalous (power-law α)"       = "anomalous"),
                          selected = "linear"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "The OLSF options below apply to Linear only; the other models fit the whole window."),
              conditionalPanel(
                condition = "input.trk_uni_model != 'linear'",
                div(style = "margin-top:4px;padding:5px 8px;border-left:3px solid #e67e22;background:rgba(230,126,34,0.06);font-size:11px;color:#7f8c8d;line-height:1.4;",
                    HTML(paste0(
                      "<b style='color:#e67e22;'>Exploratory.</b> The corral and power-law forms are ",
                      "standard (Kusumi et&nbsp;al. 1993; Metzler et&nbsp;al. 2014), but their fit to short, ",
                      "cycle-skipped MINFLUX tracks is not established. Treat L_c, v and \u03b1 as exploratory.")))),
              checkboxInput("trk_msd_geom",
                          HTML("Overlay <b>free diffusion on the cell geometry</b>"), FALSE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Free diffusion on the ", tags$b("Sample geometry"), " shape at the fitted D. ",
              "A curve below it shows confinement beyond the shape."),
            checkboxInput("trk_uni_iterative",
                            "Iterative OLSF (converges on lowest lags)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Refits with an updated optimal point count until it settles. Off: one pass."),
              checkboxInput("trk_uni_drop_negd",
                            "Drop D≤0 tracks from ensemble (QC)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "eDIFF keeps every track by default; tDIFF always keeps D > 0 only."),
              checkboxInput("trk_uni_bootstrap",
                            "Bootstrap 95% CI on eDIFF (slower)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Resamples tracks about 300 times; adds a few seconds."),
              checkboxInput("trk_uni_show_band", "Show IQR band", TRUE),
              checkboxInput("trk_uni_show_ntau", "Show N(τ)", TRUE),
              checkboxInput("trk_uni_show_spag", "Show individual track TAMSDs (spaghetti)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Every track as a faint curve, coloured by cluster when clustering is on."),
              checkboxInput("trk_msd_show_bindiag",
                            "Show cycle-binning diagnostic", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Pairs per cycle-lag, and real elapsed time over lag \u00d7 cycle, which should be flat.")
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
                tags$b("L"), " comes from ", tags$b("Acquisition Parameters"), " (the repeated tracking ",
                "iteration, usually 100\u2013300 nm). ", tags$b("t\u2095\u2090"), " is about 0.05 ms for a ",
                "hexagonal pattern; ECO, EFO and \u0394t come from the data."
              ),
              uiOutput("trk_dmax_panel")
            ),
            # Per-track D histogram sub-panel
            conditionalPanel(
              condition = "input.trk_msd_show_hist == true",
              tags$hr(style = "margin:10px 0;"),
              plotlyOutput("trk_msd_hist_plot", height = "300px"),
              helpText(style = "font-size:0.82em; color:#666;",
                "One bar per track. A broad or two-humped distribution points to sub-populations.")
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
                "Below the dashed line (CI = 1): confined; above: directed. Cluster 1 is the slowest ",
                "and most confined.")
            )
          )
        )
      ),

      tags$hr(),

      # ── Confinement (MSD confinement index) ───────────────────────────────────
      fluidRow(
        column(4,
          h4("Confinement (MSD index)"),
          .method_box("green",
            "Per-track confinement index CI = (MSD_late / MSD_early) / (\u03c4_late / \u03c4_early) against a simulated free null (after Simson, Sheets & Jacobson 1995). The primary call on confinement in this tab.",
            HTML("CI compares each track's own short and long lags, so it needs no external D. The free range is simulated for each dataset from the measured &sigma;<sub>loc</sub>, D, cycle time and dark-gap pattern, matched to each track's length and D. <b>Confined</b> is below the free 5th percentile, <b>directed</b> above the 95th."),
            HTML("Confinement is visible only when a free molecule would travel further than the domain in one track; that distance, \u221a(2dDT), is reported with the result. <b>Confined&nbsp;% is not comparable between datasets with different track lengths or mobilities.</b> The directed tail is insensitive: use Directionality or Speed &amp; Transport for transport.")),
          checkboxInput("trk_show_conf", "Enable confinement index", FALSE),
          conditionalPanel(
            condition = "input.trk_show_conf == true",
            helpText(style = "font-size:0.82em; color:#666;",
              "Run Plot MSD first."),
            numericInput("trk_conf_ci_nsim", "Null simulations per length", value = 200,
                         min = 50, max = 500, step = 50),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Free tracks simulated per grid point. More is smoother and slower."),
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
            HTML("The segmented signal is <b>log(r\u00b2/\u0394t)</b>, scaled to unit variance, so irregular sampling is handled and the BIC penalty is calibrated. Segment D = mean(r\u00b2/\u0394t)/(2d), corrected for noise and motion blur."),
            HTML("The mixture treats every segment as equally precise, so one population split into segments of very different lengths can return two components. Treat the component count as exploratory.")),
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
              "An immobile segment then reads close to zero."),
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
            "Jump-distance distribution fitted with a mixture of chi distributions, Rayleigh in 2D and Maxwell in 3D (Sch\u00fctz, Schindler & Schmidt 1997).",
            HTML("Irregular \u0394t is handled either by normalizing each jump by <b>\u221a\u0394t</b>, or by keeping <b>integrated steps only</b>, which all cover the same \u0394t. D is corrected for localization error and motion blur."),
            HTML("The \u221a\u0394t mode assumes free diffusion; use integrated steps when the motion is in question. Each component is modelled as free Brownian motion, so keep the count low. A component near D = 0 is noise, not a slow population.")),
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
            .detail(label = "Fast acquisition (below about 1 ms)",
              HTML("Single steps become noise-dominated. <b>Bin-then-median</b> takes the median position in time bins 2\u20135\u00d7 the modal \u0394t wide, which brings slow components back out. <b>Subsample</b> takes every N-th localization and averages nothing away.")),
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
              "Steps with no dark gap only, so every jump covers the same \u0394t. Needs the cycle time."),
            tags$hr(style = "margin:6px 0;"),
            numericInput("trk_jdd_min_dt", "Min Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 0.1),
            numericInput("trk_jdd_max_dt", "Max Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 1),
            .detail(label = "Setting the \u0394t window",
              HTML("Put <b>Min \u0394t</b> just below the modal \u0394t from TSA, so sub-cycle steps cannot create a spurious fast component, and <b>Max \u0394t</b> at 3\u20135\u00d7 the mode, so untracked gaps stay out.")),
            numericInput("trk_jdd_sigma_loc", HTML("Localization error &sigma;<sub>loc</sub> (nm, 0 = automatic)"),
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
              "Bars: step sizes pooled over all tracks. Dashed: fitted components. Red: their sum."),
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
            "Tool-original. Each track's net start-to-end displacement R_net against the free-diffusion spread at the ensemble D.",
            HTML("R<sub>net</sub> uses only the two endpoints and the true elapsed time, so irregular \u0394t does not affect it. The endpoint noise is subtracted, and each track is judged against the whole free distribution (Rayleigh in 2D, Maxwell in 3D)."),
            HTML("The null is an <b>infinite plane</b> at one ensemble D: a curved cell, or a wide spread of per-track D, pushes free tracks below the band. Step 2 adds the cell shape. Read this panel as corroboration of the confinement index, which is the primary call.")),
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
              "Shape and degrees of freedom come from ", tags$b("Sample geometry"), " above."),
            selectInput("trk_disp_dsrc", "Diffusion coefficient D (sets the null)",
              choices = c("CVE, short-lag and noise-corrected (recommended)" = "cve",
                          "MSD fit (current MSD model)" = "msd",
                          "Manual value" = "manual"),
              selected = "cve"),
            conditionalPanel("input.trk_disp_dsrc == 'manual'",
              numericInput("trk_disp_dman", "D (µm²/s)", value = 0.1, min = 0, step = 0.01)),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
              tags$b("D sets the verdict."), " Use the short-time D (CVE): a D from an MSD fit ",
              "already carries the confinement bend and under-calls it."),
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
                tags$b("Sample geometry"), " shape. Below it: confined beyond what the shape explains.")
            ),
            tags$br(),
            downloadButton("trk_dl_disp_csv", "Download displacement table (CSV)"),
            helpText(style = "font-size:0.8em; color:#666; margin-top:4px;",
              "One row per track, with the ensemble D and ", sigloc(), ".")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_disp == true",
            plotlyOutput("trk_disp_plot", height = "380px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "Grey band: 5th\u201395th percentile of free diffusion. ",
              tags$b(style = "color:#2c7fb8;", "Blue"), ": sub-diffusive. ",
              tags$b(style = "color:#e31a1c;", "Red"), ": directed. About 5% of free tracks fall in each tail."),
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
            "Turn-angle distribution and velocity autocorrelation against an isotropic-Brownian null (Burov et al. 2013).",
            HTML("The autocorrelation is binned by <b>elapsed time</b>, so a skipped cycle lands in the right bin, and both use only steps with few dark cycles."),
            tagList("Read the curves pooled across tracks. A small negative value at the first lag is the ",
                    sigloc(), " signature and appears in free data too.")),
          checkboxInput("trk_show_dir", "Enable directionality analysis", FALSE),
          conditionalPanel(
            condition = "input.trk_show_dir == true",
            numericInput("trk_dir_min_pts", "Min points per TID", value = 6, min = 4, step = 1),
            numericInput("trk_dir_max_dark", "Max dark cycles per step", value = 2, min = 0, step = 1),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Longer gaps are left out, since they can hide a merge with another molecule."),
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
              "Peak at 0\u00b0: directed. Peak at 180\u00b0: confined or tethered. Flat: Brownian."),
            plotlyOutput("trk_dir_autocorr_plot", height = "300px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "Positive and slowly decaying: directed. A small first-lag dip, then about zero: Brownian."),
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
            "Continuous-time hidden Markov model fitted by maximum likelihood over all tracks (Das, Cairo & Coombs 2009; Persson et al. 2013; Jackson 2011).",
            HTML("Each step's transition matrix is <b>exp(Q\u00b7\u0394t)</b>, so dark cycles are handled exactly and the <b>switching rates are fitted parameters</b>. D is corrected for noise and motion blur."),
            tagList("Rates are reliable when a state lasts about ten cycles or more; states outside that range are flagged. ",
                    "Report the bootstrap intervals: the likelihood ones assume identical molecules. The model reads x and y only.")),
          checkboxInput("trk_show_hmm", "Fit HMM to tracks", FALSE),
          conditionalPanel(
            condition = "input.trk_show_hmm == true",
            numericInput("trk_hmm_n_states",    "Number of states",                value = 2,   min = 2, max = 4,   step = 1),
            numericInput("trk_hmm_min_pts",     "Min localizations per TID",       value = 10,  min = 4, step = 1),
            numericInput("trk_hmm_n_ex",        "Example tracks to display",       value = 6,   min = 1, max = 30,  step = 1),
            numericInput("trk_hmm_max_steps",   "Step budget for the fit (0 = all steps)",
                         value = 150000, min = 0, step = 25000),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Above it, whole tracks are drawn at random. More steps cost time without narrowing D."),
            checkboxInput("trk_hmm_bic", "Also fit one state fewer and one more (BIC)", FALSE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Two further fits, so about three times the run time."),
            numericInput("trk_hmm_boot",        "Bootstrap refits over tracks (0 = off)", value = 0, min = 0, max = 200, step = 10),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Each refit takes as long as the main fit. Use 50 or more for reported intervals."),
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
              "One bar per state, slow to fast, with its 95% interval and share of time."),
            plotlyOutput("trk_hmm_trans_plot",   height = "260px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "Probability of the column state one cycle after the row state. A strong diagonal means long dwells."),
            plotlyOutput("trk_hmm_example_plot", height = "400px"),
            helpText(style = "font-size:0.82em; color:#666;",
              "The longest tracks, coloured by most likely state."),
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
            "Pools the same-state runs of each HMM state into a per-state ensemble MSD and fits it with the chosen MSD model.",
            "Binning and fit window come from the MSD panel. Fitting each state separately keeps a slow and a fast phase from being averaged together.",
            "States are assigned from step sizes, and short runs flatten the MSD, so the slow state can show confinement that the selection produced. It is only as good as the HMM fit."),
          checkboxInput("trk_show_seg_msd", "Enable per-segment MSD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_seg_msd == true",
            numericInput("trk_seg_msd_min_pts", "Min points per segment", value = 15, min = 5, step = 5),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Run the HMM first."),
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
              "One eMSD per state, slow to fast; the table sets MSD D beside HMM D.")
          )
        )
      ),

      tags$hr(),

      # ── Speed & Transport (exploratory - kept late in the tab) ────────────────
      fluidRow(
        column(4,
          h4("Speed & Transport"),
          .method_box("amber",
            "Tool-original. Per-step path speed, a drift velocity averaged over a time window W, and a noise-corrected local exponent \u03b1(t). Sustained fast segments are flagged and scored against evidence from outside the speed signal.",
            tagList("Speeds use each step's real \u0394t, and the local \u03b1 is corrected with ", sigloc(), "."),
            HTML("Path speed cannot separate diffusion from transport, so segments are flagged on drift velocity. A tracker that has lost its molecule and is searching looks similar; a real event keeps a high integrated fraction, while a search run collects dark cycles. <b>Only a matched negative control settles artefact against signal.</b>")),
          checkboxInput("trk_show_speed", "Enable speed & transport", FALSE),
          conditionalPanel(
            condition = "input.trk_show_speed == true",

            tags$b("Finding drift segments", style = "font-size:0.9em;"),
            numericInput("trk_speed_drift_win", "Drift window W (ms)",
                         value = 10, min = 1, step = 2),
            numericInput("trk_speed_fold", "Threshold (× diffusive baseline)",
                         value = 3, min = 1, step = 0.5),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "Called when the drift velocity stays above this multiple of the diffusive baseline for ",
              "the whole window. A larger W rejects more noise and blurs the onset."),
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
                  "Sets the \u221a2\u03c3/W noise floor and the \u03b1 correction. 0 = the app-wide ", sigloc(), "."),
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
              "One dot per candidate segment: where it starts along the track, against its integrated ",
              "fraction. A real event sits high; a tracker search sits low."),
            fluidRow(
              column(6, plotlyOutput("trk_export_hist",       height = "290px")),
              column(6, plotlyOutput("trk_export_onset_hist", height = "290px"))
            ),
            helpText(style = "font-size:0.82em; color:#666;",
              tags$b("Left"), ": the chosen metric, grey unflagged and green flagged. ",
              tags$b("Right"), ": where along each track the motion is strongest; a pile-up near 1 ",
              "occurs in real and artefactual data alike."),
            tags$details(style = "margin-top:8px;",
              tags$summary(tags$b("Inspect a single track")),
              uiOutput("trk_speed_tid_ui"),
              tags$div(style = "margin-bottom:4px;",
                checkboxInput("trk_speed_show_drift", "Drift (transport) velocity", TRUE),
                checkboxInput("trk_speed_show_alpha", "Local exponent α(t)", TRUE),
                checkboxInput("trk_speed_show_path",  "Path speed (reference)", TRUE)),
              plotlyOutput("trk_speed_plot", height = "340px"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:6px;",
                tags$b(style = "color:#1B9E5A;", "Green"), ": drift velocity, noise floor dashed orange, ",
                "threshold dashed red. ", tags$b(style = "color:#7B4FBE;", "Purple"), ": local \u03b1 ",
                "(2 directed, 1 Brownian). ", tags$b(style = "color:#888;", "Grey dotted"), ": path speed."),
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
            HTML("Single-lag diffusivity per spatial cell, D = [mean(r\u00b2/4\u0394t) \u2212 \u03c3\u00b2\u00b7mean(1/\u0394t)] / [1 \u2212 2R\u00b7t<sub>c</sub>\u00b7mean(1/\u0394t)], the short-lag limit of the CVE (Vestergaard 2014; Vogler 2025)."),
            tagList("The noise floor ", sigloc(), "\u00b2/\u0394t is subtracted and motion blur corrected once per cell. ",
                    "'Integrated only' drops dark-cycle steps, whose midpoint says little about position."),
            "An instantaneous D, blind to confinement. A mis-linked jump inflates a cell's mean, so check the per-cell SE in the hover. Uneven sampling across the field can imprint its own pattern."),
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
              "Each step is assigned to the cell holding its midpoint."),
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
              "Local D in nm\u00b2/ms, cool for slow and warm for fast. Blank: too few steps."),
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
