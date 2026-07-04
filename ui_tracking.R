## ─────────────────────────────────────────────────────────────────────────────
## ui_tracking.R  —  Tracking tab UI (assigned to tracking_ui)
## All input/output IDs prefixed with trk_ to avoid collisions with SMLM IDs.
## ─────────────────────────────────────────────────────────────────────────────
tracking_ui <- fluidPage(
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
            "Select multiple files to merge them into one dataset. Track IDs are ",
            "rescaled (prefixed F1_, F2_, …) so IDs from different files never collide."),
          numericInput(
            "trk_tids_per_chunk",
            "TIDs per chunk (for big-file splitting)",
            value = TIDS_PER_CHUNK_DEFAULT,
            min = 50, step = 50
          ),
          helpText("Lower values = more, smaller chunks (easier on RAM); higher values = fewer, larger chunks."),
          uiOutput("trk_chunk_ui"),
          uiOutput("trk_chunk_stat"),
          ),  # /Data Upload
          tags$hr(),

          # ── Acquisition parameters (optional cycle-time model) ──────────────
          tags$details(class = "sim-group",
          tags$summary(tags$b("Acquisition Parameters")),
          checkboxInput("trk_hw_enable",
                        "Use hardware model for cycle time (optional)", FALSE),
          checkboxInput("trk_cycle_use_regression",
                        "Fit cycle time from data (EOD/SW regression)", FALSE),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Line fit of the cluster centres (cycle = slope, SW = intercept), more precise ",
            "than the cluster median. Off = cluster median (current default). ",
            "Both are shown below for comparison; the active one is marked ▶."),
          numericInput("trk_msd_cycle_ms",
                       "Cycle time manual override (ms, 0 = use auto)",
                       value = 0, min = 0, step = 0.01),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Set > 0 to force a specific cycle time (overrides hardware model and the ",
            "data-driven estimates). Drives gap-splitting and the MSD."),
          uiOutput("trk_cycle_time_status"),
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
              tags$b("Cycle time = dwell + EOD overhead × vertices"),
              " Time Signature Analysis (TSA) estimates the cycle time from your data ",
              " and predicts theoretical cycle times from your hardware settings. ",
              " The EOD overhead is the donut move between each pattern vertex ",
              " typically ~5 µs. The number of vertices depends on your scan pattern: ",
              " the cycle time. When enabled, drives gap splitting and MSD calculations, ",
              " and Integrated/Dark classification.")
          ),
          ),  # /Acquisition Parameters
          tags$hr(),

          tags$details(class = "sim-group", open = NA,
          tags$summary(tags$b("Metadata Filtering")),
          uiOutput("trk_lifetime_stat"),
          checkboxInput("trk_show_metadata_hist", "Show Filtering Controls", FALSE),
          conditionalPanel(
            condition = "input.trk_show_metadata_hist == true",

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

            checkboxInput("trk_split_gaps_enable",
                          "Split tracks at time gaps", FALSE),
            conditionalPanel(
              condition = "input.trk_split_gaps_enable == true",
              checkboxInput("trk_split_adaptive",
                            "Adaptive threshold", FALSE),
              conditionalPanel(
                condition = "input.trk_split_adaptive == false",
                numericInput("trk_split_gap_cycles",
                             "Gap threshold (cycles) — split if interval exceeds this",
                             value = 7, min = 2, step = 1)
              ),
              conditionalPanel(
                condition = "input.trk_split_adaptive == true",
                numericInput("trk_split_min_len",
                             "Min localizations per segment", value = 50, min = 1, step = 1),
                uiOutput("trk_split_adaptive_hint")
              ),
              uiOutput("trk_split_gap_ms_hint"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Threshold in acquisition cycles: an interval spanning more than this many ",
                "cycles is treated as a gap and the track is split there. ",
                tags$b("Adaptive"), " Threshold is set from the data as ",
                "median + SD + 1 of the per-step total-cycle counts (capped by the hardware maxOffTime), ",
                "and segments shorter than the min-length are discarded. ",
                "The cycle time is auto-detected. Segment 1 of track 42 becomes '42_s1', '42_s2', etc.; ",
                "length and duration filters then apply per segment.")
            ),

            tags$hr(),

            uiOutput("trk_tid_len_ui"),

            tags$hr(),

            uiOutput("trk_dur_filter_ui"),

            fluidRow(
              column(6, actionButton("trk_reset_filters",    "Reset ranges")),
              column(6, actionButton("trk_rechunk_filtered", "Filter & Rechunk", class = "btn-warning"))
            ),
            helpText("Filter & Rechunk loads all chunks, applies the current metadata/length/duration filters, and saves new chunks containing only the filtered traces."),

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
            helpText("Zoom into a region on the plot, then press Crop to view. ",
                     "TIDs with at least the specified fraction of localizations ",
                     "inside the zoomed area are kept. Reset Ranges also clears this filter."),
            tags$hr(style = "margin:8px 0;"),

            downloadButton("trk_download_filtered_csv", "Export filtered data (.csv)"),
            helpText("Exports filtered data from the current chunk."),
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
              numericInput("trk_bin_dt",         "Bin width (ms)",             value = 20,   min = 0.0001, step = 0.001),
              helpText("Frame count is derived automatically: window ÷ bin width."),

              checkboxInput("trk_playback_ctrl", "Show/Hide Playback Controls", TRUE),
              conditionalPanel(
                condition = "input.trk_playback_ctrl == true",
                numericInput("trk_trail_len",     "Trail length (ms)",                  value = 3, min = 0,    step = 1),
                numericInput("trk_ms_per_unit",   "ms per time unit (per-TID mode)",    value = 1,   min = 0.01, step = 0.01),
                checkboxInput("trk_show_advanced", "Show advanced options", FALSE),
                conditionalPanel(
                  condition = "input.trk_show_advanced == true",
                  numericInput("trk_target_markers", "Target markers per frame",                value = 5000, min = 100, step = 500),
                  numericInput("trk_per_tid_cap",  "Per-TID cap per frame (0 = off)",          value = 0,   min = 0, step = 50),
                  numericInput("trk_trail_nmax",   "Max points per trail (per TID, per frame)", value = 50,  min = 5, step = 5),
                  numericInput("trk_max_pts_tid",  "Pre-decimate max points per TID",           value = 400, min = 50, step = 50)
                )
              ),

              checkboxInput("trk_show_dt_hist", "Show Δt histogram (per TID)", FALSE),
              conditionalPanel(
                condition = "input.trk_show_dt_hist == true",
                plotlyOutput("trk_dt_hist", height = "220px")
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
          uiOutput("trk_view_mode_ui"),
          numericInput("trk_pt_size", "Marker size", value = 6, min = 1),
          conditionalPanel(
            condition = "input.trk_view_mode_main == 'playback'",
            checkboxInput("trk_show_full_tracks_faint", "Show full tracks faintly", TRUE),
            helpText("Markers & (optional) static lines use one color per TID."),
            checkboxInput("trk_show_moving_trails", "Animate moving trails", FALSE)
          ),
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
                     tags$p("Pixel size must be set first. Results are a starting point — nudge to refine.",
                            style = "font-size:11px; color:#888; margin:0 0 6px 0;"),
                     checkboxInput("trk_img_invert_y", "Image origin at top-left (invert Y)", TRUE),
                     checkboxInput("trk_img_flip_x",   "Mirror horizontally (left↔right)", FALSE),
                     sliderInput("trk_img_opacity", "Image opacity", min = 0, max = 1, value = 0.6, step = 0.05),
                     selectInput("trk_img_rotation", "Rotate image",
                                 choices = c("0","90","180","270"), selected = "0"),
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

      # ── MSD ─────────────────────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("MSD (diffusion)"),
          tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin-bottom:8px;",
            tags$b("How MSD is computed:"), tags$br(),
            "Per-track displacements are binned by acquisition cycle (cumulative cycle index), ",
            "then averaged into the ensemble MSD(τ) and fit with OLSF. Reports the ensemble ",
            "diffusion coefficient (eDIFF) and the per-track mean (tDIFF). ",
            "Method follows Vogler et al. (2025, ", tags$em("Commun. Biol."), ")."
          ),
          checkboxInput("trk_show_msd", "Enable MSD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_msd == true",
            numericInput("trk_msd_min_pts",  "Min points per TID (skip short)", value = 20, min = 5, step = 5),
            radioButtons("trk_msd_ndim", "Diffusion dimensionality",
              choices  = c("Auto-detect" = "auto", "2D  (MSD = 4Dτ)" = "2", "3D  (MSD = 6Dτ)" = "3"),
              selected = "auto", inline = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              tags$b("Auto:"), " uses 3D if z has real variance (SD > 1 nm), else 2D. ",
              tags$b("Choose 2D"), " for membrane-confined proteins. ",
              tags$b("Choose 3D"), " for freely diffusing cytoplasmic or nuclear proteins."),
            tags$hr(style = "margin:6px 0;"),
            checkboxInput("trk_show_dmax",     "Diffusion limit analysis (Dₘₐˣ)", FALSE),
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
            checkboxInput("trk_msd_unified",
                          "Unified MSD", TRUE),
            conditionalPanel(
              condition = "input.trk_msd_unified == true",
              tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin:4px 0;",
                tags$b("MINFLUX-optimised MSD."), " Cumulative cycle-index binning, ",
                "track-weighted ensemble (eMSD) on the linearized lag×cycle timeline, iterative ",
                "OLSF that converges on the lowest lags. Reports both eDIFF (ensemble) and tDIFF ",
                "(per-track mean). Untick the toggles below for the robust / real-time alternatives. ",
                "Method follows Vogler et al. (2025, ", tags$em("Commun. Biol."), ")."),
              fluidRow(
                column(6, numericInput("trk_uni_min_lag",
                             "Lower cutoff (cycles)", value = 5, min = 1, step = 1)),
                column(6, numericInput("trk_uni_max_lag",
                             "Upper fit limit (cycles)", value = 50, min = 10, max = 5000, step = 5))
              ),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Fit window = [lower, upper] cycle-lags (default 5–50); OLSF + Michalet pick the ",
                "optimal point count within it. Lower skips the linker/noise regime. ",
                tags$b("The upper limit also sets how many lags are computed"), " (capped 5000) — ",
                "raise it (e.g. 500–2000) to reach a confinement plateau, then ", tags$b("re-click Plot MSD"), "."),
              checkboxInput("trk_uni_track_weight",
                            "Track-weighted ensemble (eMSD)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "On (default) = mean of per-track TAMSDs (eMSD). ",
                "Off = pair-weighted (smoother — long tracks contribute more at long lags). ",
                "Re-click Plot MSD to apply."),
              checkboxInput("trk_uni_linear_tau",
                            "Linearized τ axis (lag × cycle)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "On (default) = τ = lag × cycle time (linearized timeline). ",
                "Off = each bin at its true mean Δt (real time, cycle-time-independent D). ",
                "Applies instantly — no recompute needed."),
              selectInput("trk_uni_model", "Fit model",
                          choices = c("Linear (Brownian)"             = "linear",
                                      "Confined (plateau → size L)"   = "confined",
                                      "Directed (drift → velocity v)" = "directed",
                                      "Anomalous (power-law α)"       = "anomalous"),
                          selected = "linear"),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                tags$b("Linear"), " = OLSF Brownian D (the OLSF options below apply to this). ",
                tags$b("Confined"), " fits a plateauing corral → local D + confinement size L. ",
                tags$b("Directed"), " fits Brownian + drift → D + velocity v. ",
                tags$b("Anomalous"), " fits MSD = Dτ^α (α<1 sub-diffusive/confined, α>1 directed). ",
                "Non-linear models fit the whole [lower, upper] window; apply on Plot MSD."),
              checkboxInput("trk_uni_iterative",
                            "Iterative OLSF (converges on lowest lags)", TRUE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "On (default) = the iterative loop starts at a few points and iterates the optimal ",
                "point count until it repeats, settling on the lowest, least-noisy lags. ",
                "Off = robust full-window-anchored Michalet. Applies instantly."),
              checkboxInput("trk_uni_drop_negd",
                            "Drop D≤0 tracks from ensemble (QC)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Off (default) = keep every track in the ensemble (standard eDIFF convention). ",
                "On = remove non-diffusing tracks from the ensemble too. The per-track ",
                "tDIFF always keeps only D>0, regardless. Re-click Plot MSD to apply."),
              checkboxInput("trk_uni_show_band", "Show IQR band", TRUE),
              checkboxInput("trk_uni_show_ntau", "Show N(τ)", TRUE),
              checkboxInput("trk_uni_show_spag", "Show individual track TAMSDs (spaghetti)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Faint per-track MSD curves. When 'Cluster tracks by D & CI' is on, they are ",
                "coloured by cluster and a per-cluster eMSD + tDIFF is added.")
            ),
            tags$br(),
            actionButton("trk_msd_compute", "Plot MSD",
                         icon = icon("chart-line"),
                         class = "btn-primary btn-sm",
                         style = "margin-bottom:4px;"),
            helpText(style = "margin-top:0; font-size:0.82em; color:#666;",
                     "Adjust settings, then click to compute."),
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
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:2px;",
              tags$b("How to read the MSD:"), tags$br(),
              "x-axis = time lag (τ); y-axis = mean squared displacement (µm²). ",
              tags$b("Straight line through the origin"), " = free Brownian diffusion; the slope gives D (MSD = 4Dτ in 2D, 6Dτ in 3D). ",
              tags$b("Bends down / flattens"), " = confined diffusion. ",
              tags$b("Curves upward"), " = directed motion."
            ),
            # Unified MSD sub-panel (primary — shown first)
            conditionalPanel(
              condition = "input.trk_msd_unified == true",
              plotlyOutput("trk_msd_unified_plot", height = "300px"),
              uiOutput("trk_msd_unified_estimates")
            ),
            # Diffusion limit sub-panel
            conditionalPanel(
              condition = "input.trk_show_dmax == true",
              tags$hr(style = "margin:10px 0;"),
              tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin-bottom:8px;",
                tags$b("What this does:"), tags$br(),
                "Estimates the maximum trackable diffusion coefficient Dₘₐˣ from acquisition parameters. ",
                "Based on Vogler et al. (2025, ", tags$em("Commun. Biol."), ")."
              ),
              fluidRow(
                column(6, numericInput("trk_dmax_L",   "TCP diameter L (nm)",         value = 150,  min = 20, max = 2000, step = 10)),
                column(6, numericInput("trk_dmax_t_hw","Hardware overhead tₕₐ (ms)", value = 0.05, min = 0,  max = 1,    step = 0.01))
              ),
              helpText(style = "font-size:0.82em; color:#666;",
                "L is the TCP scan diameter. tₕₐ is hardware overhead per cycle (~0.05–0.1 ms). Median Δt and EFO are read from loaded data."
              ),
              uiOutput("trk_dmax_panel")
            ),
            # Per-track D histogram sub-panel
            conditionalPanel(
              condition = "input.trk_msd_show_hist == true",
              tags$hr(style = "margin:10px 0;"),
              plotlyOutput("trk_msd_hist_plot", height = "300px"),
              tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
                tags$b("How to read:"), tags$br(),
                "Each bar is one track’s fitted D (log₁₀ scale). ",
                "A narrow peak means a homogeneous population; broad or bimodal means fast and slow sub-populations."
              )
            ),
            # Cluster sub-panel
            conditionalPanel(
              condition = "input.trk_msd_cluster == true",
              tags$hr(style = "margin:10px 0;"),
              plotlyOutput("trk_msd_cluster_plot", height = "320px"),
              tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
                tags$b("How to read:"), tags$br(),
                "X-axis: log₁₀(D). Y-axis: confinement index (CI, log scale).", tags$br(), tags$br(),
                tags$b("What CI measures:"), tags$br(),
                "For free Brownian diffusion, MSD(τ) grows linearly in time, so the ratio ",
                "MSD(τ₂)/MSD(τ₁) equals τ₂/τ₁ exactly. ",
                "CI is defined as:", tags$br(),
                tags$code("CI = [MSD(τ_late) / MSD(τ_early)] / [τ_late / τ_early]"),
                tags$br(),
                "where τ_early and τ_late are the median lag times of the first and last 15% of ",
                "lag bins for that track (minimum 2 bins each; medians used). ",
                tags$b("CI = 1"), ": MSD ∝ τ - free Brownian. ",
                tags$b("CI ≪ 1"), ": MSD growth slows or plateaus - confined or caged. ",
                tags$b("CI ≫ 1"), ": MSD grows faster than linear - directed or active transport. ",
                "The dashed line marks CI = 1.", tags$br(), tags$br(),
                tags$b("Cluster 1"), " = most confined / slowest; higher clusters = faster / more mobile. ",
                "D and CI values are available in the per-track CSV download."
              )
            )
          )
        )
      ),

      tags$hr(),

      # ── Instantaneous speed ───────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Speed & Transport"),
          tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin-bottom:8px;",
            tags$b("Three per-track curves:"), tags$br(),
            tags$b("Path speed"), " vᵢ = |Δrᵢ|/Δtᵢ — dominated by diffusion + localization ",
            "noise, so it barely responds to transport. ",
            tags$b("Drift velocity"), " = |net displacement over a time window W| / W — measures ",
            "coherent transport and is what flags an export event. ",
            tags$b("Local α(t)"), " = rolling log-log MSD slope (1 = Brownian, 2 = directed) — the ",
            "per-track analog of a TAMSD bending upward."),
          checkboxInput("trk_show_speed", "Enable speed & transport", FALSE),
          conditionalPanel(
            condition = "input.trk_show_speed == true",
            numericInput("trk_speed_min_pts", "Min points per TID", value = 20, min = 5, step = 5),

            tags$hr(style = "margin:6px 0;"),
            tags$b("Export detection", style = "font-size:0.9em;"),
            numericInput("trk_speed_drift_win", "Drift window W (ms) — larger = more noise rejection, blurrier onset",
                         value = 200, min = 1, step = 10),
            numericInput("trk_speed_drift_minpts", "Min localizations in window",
                         value = 5, min = 2, step = 1),
            numericInput("trk_speed_sigma_nm", "Localization precision σ (nm) — noise floor √2σ/W & de-biases α",
                         value = 5, min = 0, step = 1),
            numericInput("trk_speed_fold", "Export threshold (× diffusive baseline) — higher = stricter",
                         value = 3, min = 1, step = 0.5),
            numericInput("trk_speed_onset_thr", "…or absolute threshold (µm/s, 0 = use × baseline above)",
                         value = 0, min = 0, step = 0.1),
            numericInput("trk_speed_min_run", "Min sustained windows to call an export",
                         value = 3, min = 1, step = 1),
            numericInput("trk_speed_edge_guard", "Edge guard (× window W) — ignore exports this close to track start/end",
                         value = 1, min = 0, max = 5, step = 0.5),
            selectInput("trk_export_hist_metric", "Distribution metric",
                        choices = c("Peak drift velocity"        = "v_drift_peak",
                                    "Export-segment velocity"     = "export_v",
                                    "Fraction of track directed"  = "directed_frac",
                                    "Fraction immobile/confined"  = "immobile_frac",
                                    "Transport episodes / track"  = "n_events",
                                    "Max local α"                 = "alpha_max"),
                        selected = "v_drift_peak"),
            downloadButton("trk_dl_export_csv", "Per-track summary (CSV)", class = "btn-xs btn-primary"),

            tags$hr(style = "margin:6px 0;"),
            tags$b("Local α(t)", style = "font-size:0.9em;"),
            numericInput("trk_speed_alpha_win", "α window (steps)", value = 21, min = 7, max = 201, step = 2),
            numericInput("trk_speed_alpha_maxlag", "α max lag (steps) — raise if α under-reads on strong diffusion",
                         value = 6, min = 2, max = 30, step = 1),
            numericInput("trk_speed_smooth", "Rolling median window (steps)",
                         value = 5, min = 1, max = 51, step = 2)
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_speed == true",
            tags$div(class = "well well-sm", style = "margin-bottom:6px;",
                     uiOutput("trk_export_summary_txt")),
            fluidRow(
              column(6, plotlyOutput("trk_export_hist",       height = "300px")),
              column(6, plotlyOutput("trk_export_onset_hist", height = "300px"))
            ),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("How to read:"), " Left = distribution of the chosen per-track metric ",
              "(grey = non-exporters, green = exporters); a separated green mode = a real ",
              "transported subpopulation. Right = where transport is strongest along each track ",
              "(0 = start, 1 = end, from the peak drift velocity); a pile-up near 1 confirms ",
              "exports at track ends. ", tags$b("Edge guard"), " ignores the first/last window of ",
              "each track — track starts carry localization jitter that otherwise piles false ",
              "events at position 0. Tune W, σ, threshold, min-run and edge guard on the left; the ",
              tags$b("Per-track summary (CSV)"), " has one row per track for your own stats."
            ),
            tags$details(style = "margin-top:8px;",
              tags$summary(tags$b("Inspect a single track")),
              uiOutput("trk_speed_tid_ui"),
              tags$div(style = "margin-bottom:4px;",
                checkboxInput("trk_speed_show_drift", "Drift (transport) velocity", TRUE),
                checkboxInput("trk_speed_show_alpha", "Local exponent α(t)", TRUE),
                checkboxInput("trk_speed_show_path",  "Path speed (reference)", TRUE)),
              plotlyOutput("trk_speed_plot", height = "340px"),
              tags$div(class = "well well-sm", style = "font-size:0.82em; margin-top:6px;",
                tags$b(style = "color:#1B9E5A;", "Green = drift velocity"), " (transport); a sustained rise ",
                "above the noise floor (orange dashed) and onset threshold (red dashed, points red) = export. ",
                tags$b(style = "color:#7B4FBE;", "Purple = local α"), " → 2 directed, ~1 Brownian. ",
                tags$b(style = "color:#888;", "Grey dotted = path speed"), ", stays high/flat (noise-dominated). ",
                downloadButton("trk_dl_speed_csv", "Per-timepoint curves (CSV)", class = "btn-xs btn-default")
              )
            )
          )
        )
      ),

      tags$hr(),

      # ── Time Signature Analysis (TSA) ────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Time Signature Analysis (TSA)"),
          tags$div(style = "font-size:0.85em; color:#666; margin-bottom:6px;",
            "Recovers the acquisition cycle structure from the data: the atomic cycle ",
            "time, the Integrated/Dark cycle balance, and an overall data-quality verdict. ",
            "Independent of the MSD — use it to judge whether the data are good for ",
            "diffusion analysis before you trust D."),
          checkboxInput("trk_show_tsa", "Enable TSA", FALSE),
          conditionalPanel(
            condition = "input.trk_show_tsa == true",
            tags$b("Sequence file (optional)", style = "font-size:0.9em;"),
            fileInput("trk_seq_json", "Upload MINFLUX sequence (.json)",
                      accept = ".json"),
            uiOutput("trk_seq_json_status"),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Auto-fills the scan pattern (vertices), dimensionality, and the gap ",
              "threshold (from maxOffTime) into the Acquisition Parameters and gap-split. ",
              "The cycle time itself is recovered from the data (cluster method), so the ",
              "sequence is for cross-checking and convenience."),
            tags$hr(style = "margin:8px 0;"),
            sliderInput("trk_dt_hist_bins", "Δt histogram bins", min = 10, max = 300, value = 80, step = 5),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Display-only binning for the inter-localization Δt histogram. The median ",
              "and bin-free typical Δt are marked with vertical lines."),
            tags$hr(style = "margin:8px 0;"),
            numericInput("trk_cycle_dark_max", "Max cycles to display", value = 10, min = 4, step = 1),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Each consecutive interval is classified by how many acquisition cycles it ",
              "spans (round((Δt − offset) / cycle time)). 1 cycle = Integrated (clean); ",
              "≥2 cycles = Dark cycles preceded it (molecule briefly unseen).")
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
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("How to read this:"), tags$br(),
              "Each bar is the fraction of consecutive intervals that spanned that many ",
              "acquisition cycles. ", tags$b("1 cycle (Integrated)"), " = the molecule was ",
              "localized every cycle ", tags$b("≥2 cycles (Dark)"),
              " = background/failed cycles occurred before the next localization. ",
              "A high Integrated fraction means high-quality, gap-free tracking; a long ",
              "Dark tail indicates frequent dropouts (photophysics, low signal, or fast ",
              "motion out of range)."
            ),
            tags$hr(style = "margin:8px 0;"),
            h5("Per-track cycle progression"),
            uiOutput("trk_cycle_prog_tid_ui"),
            plotlyOutput("trk_cycle_prog_plot", height = "260px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              "Dark cycles before each localization along one track ",
              tags$b("Blue"), " = Integrated (0 dark); ",
              tags$b("red"), " = Dark cycles present. Bursts of red mark where the ",
              "molecule was repeatedly lost and re-found."
            )
          )
        )
      ),

      tags$hr(),

      # ── Track Segmentation (CPD + GMM) ───────────────────────────────────────
      fluidRow(
        column(4,
          h4("Track Segmentation — Change-Point + GMM"),
          tags$div(style = "font-size:0.85em; color:#666; margin-bottom:6px;",
            "Splits each track into diffusive-state segments using PELT change-point",
            " detection on log(r²/Δt), then fits a Gaussian mixture model (GMM) to the",
            " resulting per-segment D values."
          ),
          checkboxInput("trk_show_cpd", "Enable segmentation + GMM", FALSE),
          conditionalPanel(
            condition = "input.trk_show_cpd == true",
            fluidRow(
              column(6, numericInput("trk_cpd_min_seg", "Min segment length (steps)", value = 8L, min = 3L, max = 50L, step = 1L)),
              column(6, numericInput("trk_cpd_max_g",   "Max GMM components",         value = 4L, min = 2L, max = 6L,  step = 1L))
            ),
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
            tags$hr(style = "margin:6px 0;"),
            numericInput("trk_jdd_min_dt", "Min Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 0.1),
            numericInput("trk_jdd_max_dt", "Max Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 1),
            numericInput("trk_jdd_sigma_loc", "localization error σ_loc (nm, 0 = off)",
                         value = 0, min = 0, step = 0.5),
            downloadButton("trk_dl_jdd_csv", "Download JDD (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_jdd == true",
            plotlyOutput("trk_jdd_plot", height = "380px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("What this plot shows"), tags$br(),
              "Each bar represents how often molecules took a step of that size across all tracks. ",
              "A single peak means molecules move as one uniform population. A shoulder or two distinct peaks indicates two populations — one slow (small steps, left) and one fast (large steps, right). ",
              "The dashed coloured curves are the fitted populations; the red curve is their sum.",
              tags$br(), tags$br(),

              tags$b("Diffusion coefficient D"), tags$br(),
              "D describes how fast a molecule moves on average. It is only shown when 'Normalize steps by √Δt' is on. ",
              "Typical values: membrane proteins 0.01–0.5 µm²/s; cytoplasmic proteins 1–20 µm²/s. ",
              "Unit note: 1000 nm²/ms = 1 µm²/s.",
              tags$br(), tags$br(),

              tags$b("Recommended workflow"), tags$br(),
              tags$ol(style = "padding-left:1.2em; margin:0;",
                tags$li("Open the ", tags$b("Δt histogram"), " in the Time Signature Analysis section. The tallest bar (mode) is your instrument's typical acquisition interval."),
                tags$li("Set ", tags$b("Min Δt"), " to just below that modal value. Steps shorter than one acquisition interval are dominated by measurement noise rather than real molecular motion, and will create a spurious fast peak."),
                tags$li("Set ", tags$b("Max Δt"), " to 3–5× the modal value to exclude long gaps where the molecule was not being tracked."),
                tags$li("Run the ", tags$b("MSD analysis"), " — the reported σ_loc value (localization error) is shown directly below the MSD plot. This is the random positional error per measurement: a smaller number means more accurate positions. Copy that value into the ", tags$b("σ_loc correction"), " field above to remove its effect from the D estimates. The MSD estimates panel also shows this as a blue hint.")
              ),
              tags$br(),

              tags$b("When acquisition is fast (Δt < 1 ms)"), tags$br(),
              "At very short intervals, individual step sizes are often dominated by measurement noise rather than real motion. ",
              "Try the ", tags$b("Bin-then-median"), " step mode: this groups several consecutive localizations into time bins and uses the median position of each bin as a single cleaner estimate. ",
              "Set the bin width to 2–5× the modal Δt. This reduces noise and makes the fitted D values more reliable for slow molecules. ",
              "The 'Subsample' mode instead takes every N-th localization — simpler, but without the noise-reduction benefit of averaging.",
              tags$br(), tags$br(),

              tags$b("When to trust the result"), tags$br(),
              "The JDD reliably resolves populations when the typical step size is clearly larger than the measurement noise. ",
              "If the corrected D values for a fitted component are near zero (or are clamped to zero), that component is likely noise rather than real fast motion. ",
              "If a slow population is hard to see, use the ", tags$b("MSD analysis"), " instead — it separates diffusion from noise via the y-intercept rather than relying on individual step sizes."
            ),
            uiOutput("trk_jdd_estimates")
          )
        )
      ),

      tags$hr(),

      # ── Confinement analysis ──────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Confinement Analysis"),
          checkboxInput("trk_show_conf", "Enable confinement analysis", FALSE),
          conditionalPanel(
            condition = "input.trk_show_conf == true",
            numericInput("trk_conf_min_pts", "Min points per TID", value = 10, min = 4, step = 1),
            downloadButton("trk_dl_conf_csv", "Download confinement table (CSV)"),
            tags$hr(style = "margin:8px 0;"),
            checkboxInput("trk_conf_color_tracks", "Color 2D tracks by confinement ratio", FALSE)
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_conf == true",
            plotlyOutput("trk_conf_plot", height = "380px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("How to read this plot:"), tags$br(),
              "Each dot is one track. X = how far the molecule ever strayed from its average position (Rₘₐˣ); ",
              "Y = its typical spread around that average (Rₒ, radius of gyration). ",
              "The colour and the ratio Rₒ/Rₘₐˣ tell you the shape: ",
              tags$b("ratio ≈ 0.25 (red)"), " = confined / caged. ",
              tags$b("ratio ≈ 0.50 (blue)"), " = freely diffusing. ",
              "Hover a dot to see the exact TID and numbers.", tags$br(),
              tags$b("Track-length caveat:"), " short tracks (<~30 steps) have a downward-biased Rₘₐˣ ",
              "because they have not had time to explore their full range. ",
              "This makes short freely-diffusing tracks appear more confined than they truly are. ",
              "Filter by track length before interpreting the ratio distribution."
            ),
            uiOutput("trk_conf_summary")
          )
        )
      ),

      tags$hr(),

      # ── Directionality / autocorrelation ─────────────────────────────────────
      fluidRow(
        column(4,
          h4("Directionality / Autocorrelation"),
          checkboxInput("trk_show_dir", "Enable directionality analysis", FALSE),
          conditionalPanel(
            condition = "input.trk_show_dir == true",
            numericInput("trk_dir_min_pts", "Min points per TID", value = 6, min = 4, step = 1),
            downloadButton("trk_dl_dir_csv", "Download angle data (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_dir == true",
            plotlyOutput("trk_dir_angle_plot", height = "300px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("Turn-angle histogram — how to read it:"), tags$br(),
              "Each bar is the fraction of consecutive-step pairs that turned by that angle. ",
              tags$b("Peak near 0°"), " = directed / persistent motion. ",
              tags$b("Peak near 180°"), " = anti-correlated; typical of confinement or elastic tethering. ",
              tags$b("Flat / uniform"), " = pure Brownian: no preferred turning angle."
            ),
            plotlyOutput("trk_dir_autocorr_plot", height = "300px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("Velocity autocorrelation — how to read it:"), tags$br(),
              "C(lag): if the molecule moved in direction X, is it still moving in direction X after N steps? ",
              tags$b("Positive and slowly decaying"), " = persistent, directed motion. ",
              tags$b("Negative at lag 1 and near zero at higher lags"), " = consistent with free Brownian diffusion. ",
              tags$b("Important caveat:"), " localization noise introduces an ",
              "inherent anti-correlation at lag 1 even for purely Brownian motion ",
              "(C(1) ≈ −σ²_loc / MSD₁step). For MINFLUX data a small negative C(1) is expected ",
              "for diffusing molecules and does not by itself indicate confinement or active transport. ",
              "Only C(lag) substantially more negative than C(lag 2, 3, …) is informative."
            ),
            uiOutput("trk_dir_summary")
          )
        )
      ),

      tags$hr(),

      # ── HMM ──────────────────────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Diffusion State Segmentation (HMM)"),
          checkboxInput("trk_show_hmm", "Fit HMM to tracks", FALSE),
          conditionalPanel(
            condition = "input.trk_show_hmm == true",
            numericInput("trk_hmm_n_states",    "Number of states",                value = 2,   min = 2, max = 4,   step = 1),
            numericInput("trk_hmm_min_pts",     "Min localizations per TID",       value = 10,  min = 4, step = 1),
            numericInput("trk_hmm_n_ex",        "Example tracks to display",       value = 6,   min = 1, max = 30,  step = 1),
            numericInput("trk_hmm_max_fit_trk", "Max tracks used for fitting",     value = 100, min = 10, max = 500, step = 10),
            checkboxInput("trk_hmm_norm_dt", "Normalize steps by √Δt (recommended)", TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Reducing max tracks speeds up fitting significantly with little loss of accuracy. ",
              "100 longest tracks is usually sufficient for stable EM convergence."),
            downloadButton("trk_dl_hmm_csv", "Download state assignments (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_hmm == true",
            plotlyOutput("trk_hmm_params_plot",  height = "280px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("State parameters — how to read it:"), tags$br(),
              "Each bar is one inferred diffusion state, sorted slow to fast. The length shows D (nm²/step). ",
              "The percentage is the fraction of all steps assigned to that state across all tracks.", tags$br(),
              tags$b("Caveats for MINFLUX data:"), " (1) The HMM emission model uses 2D step magnitudes ",
              "(x, y only) regardless of whether z data are present. ",
              "(2) ", tags$b("Normalize by √Δt"), " (on by default): step sizes are divided by √Δt before training and Viterbi decoding, making states correspond to diffusion coefficients rather than raw step sizes. With this on, D is reported in nm²/ms directly. With it off, D is in nm²/step (converted via median Δt), but variable Δt within a track will create spurious state switches. ",
              "(3) The number of states K must be specified manually; no automatic BIC selection is performed."
            ),
            plotlyOutput("trk_hmm_trans_plot",   height = "260px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("Transition matrix — how to read it:"), tags$br(),
              "Each cell shows the probability of switching from the row state to the column state in one step. ",
              tags$b("High diagonal values"), " mean the molecule tends to stay in the same state for many consecutive steps. ",
              tags$b("Low diagonal / high off-diagonal"), " means rapid interconversion between states."
            ),
            plotlyOutput("trk_hmm_example_plot", height = "400px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("Example tracks — how to read it:"), tags$br(),
              "Each track is centered on its own centroid. Point colors show the Viterbi-decoded state at each localization."
            ),
            uiOutput("trk_hmm_summary")
          )
        )
      ),

      tags$hr(),

      # ── Per-segment MSD by HMM state ──────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Per-Segment MSD (by HMM state)"),
          tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin-bottom:8px;",
            "Splits each track into HMM diffusion-state segments and computes a separate ",
            "MSD per state, fit with the Unified MSD model + window. Resolves a mixed ",
            "population: fit ", tags$b("Confined"), " on the slow/wall state for a real ",
            "confinement size L, and ", tags$b("Directed"), " on the fast/exported state for a ",
            "velocity v — instead of one washed-out whole-track MSD. ",
            "Requires the Diffusion State HMM (above) enabled and computed."),
          checkboxInput("trk_show_seg_msd", "Enable per-segment MSD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_seg_msd == true",
            numericInput("trk_seg_msd_min_pts", "Min points per segment", value = 15, min = 5, step = 5),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Model, fit window, linearized-τ and iterative settings are inherited from the ",
              "Unified MSD panel — switch its 'Fit model' to Confined or Directed and re-view here."),
            downloadButton("trk_dl_seg_msd_csv", "Per-state MSD (CSV)", class = "btn-xs btn-default")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_seg_msd == true",
            plotlyOutput("trk_seg_msd_plot", height = "340px"),
            uiOutput("trk_seg_msd_summary"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("How to read:"), tags$br(),
              "One eMSD curve per HMM state (slow→fast). The dashed line is the model fit ",
              "for that state. A slow state that bends to a plateau (Confined model) gives a ",
              "confinement size L; a fast state curving upward (Directed model) gives a velocity v. ",
              "Compare D (MSD) with D (HMM) in the table as a cross-check."
            )
          )
        )
      ),

      tags$hr(),

      # ── Spatially Resolved D-map ──────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Spatially Resolved D-map"),
          checkboxInput("trk_show_dmap", "Compute D-map", FALSE),
          conditionalPanel(
            condition = "input.trk_show_dmap == true",
            numericInput("trk_dmap_bins",    "Grid resolution (N×N cells)", value = 10, min = 3, max = 500, step = 1),
            numericInput("trk_dmap_min_pts", "Min steps per cell",               value = 5,  min = 3, step = 1),
            checkboxInput("trk_dmap_show_n", "Overlay step count per cell",      FALSE),
            helpText("D per cell = median(r² / 4Δt) across all steps whose midpoint falls in that cell."),
            downloadButton("trk_dl_dmap_csv", "Download D-map table (CSV)")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_dmap == true",
            plotlyOutput("trk_dmap_plot", height = "520px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("How to read this map:"), tags$br(),
              "Each cell is coloured by its local D (nm²/ms). ",
              tags$b("Cool colours (blue/purple)"), " = slow / confined regions. ",
              tags$b("Warm colours (yellow/red)"), " = fast / freely diffusing regions. ",
              "Blank cells had too few steps to fit reliably."
            ),
            uiOutput("trk_dmap_summary")
          )
        )
      ),

      tags$br()
    )  # end tabPanel("Analysis")

  )  # end tabsetPanel
)
