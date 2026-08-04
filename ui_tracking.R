## ─────────────────────────────────────────────────────────────────────────────
## ui_tracking.R  —  Tracking tab UI (assigned to tracking_ui)
## All input/output IDs prefixed with trk_ to avoid collisions with SMLM IDs.
## ─────────────────────────────────────────────────────────────────────────────

# NOTE: .method_box() (the collapsible "Method & MINFLUX adaptation" infobox used
# throughout this file) now lives in helpers_shared.R, so the SMLM tab can use the
# identical widget — helpers_shared.R is sourced before both ui_*.R files.

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
          # Chunking exists only to keep the animation smooth on modest hardware, so
          # these controls appear in Playback mode only. Full tracks always renders the
          # complete dataset. The inputs stay registered while hidden (conditionalPanel
          # is CSS-only), so the chunk size still applies at upload and the selected
          # chunk is remembered when switching back to Playback.
          conditionalPanel(
            condition = "input.trk_view_mode_main == 'playback'",
            numericInput(
              "trk_tids_per_chunk",
              "TIDs per chunk (for big-file splitting)",
              value = TIDS_PER_CHUNK_DEFAULT,
              min = 50, step = 50
            ),
            helpText("Lower values = more, smaller chunks (easier on RAM); higher values = fewer, larger chunks."),
            uiOutput("trk_chunk_ui"),
            uiOutput("trk_chunk_stat")
          ),
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
          numericInput("trk_sigma_loc_nm",
                       HTML("Localization precision &sigma;<sub>loc</sub> (nm, 0 = auto)"),
                       value = 0, min = 0, step = 0.5),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "The short-lag noise floor. Leave at 0 to use the app-wide ", sigloc(), ", which comes from the ",
            tags$b("Data Quality & Calibration"), " panel (CVE — MSD-free) when you have run it, ",
            "falling back to the MSD-fit intercept otherwise; set > 0 to override. Used by the D-map, ",
            "JDD, confinement, HMM, directionality and speed to subtract the localization-error floor ",
            "so slow regions aren't reported as diffusing."),
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
          tags$summary(tags$b("Track Cleaning & Filtering")),
          helpText(style = "font-size:0.8em; color:#666; margin:2px 0 6px 0;",
            "Order: Split → Trim ends → Quality gate → Metadata → Length → Duration. ",
            "Cleaning runs on the raw acquisition timing; metadata and length/duration ",
            "then apply to the cleaned tracks."),
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
                             "Gap threshold (cycles) — split if interval exceeds this",
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
                "Threshold in acquisition cycles: an interval spanning more than this many ",
                "cycles is treated as a gap and the track is split there. ",
                tags$b("Adaptive"), " Threshold is set from the data as ",
                "median + SD + 1 of the per-step total-cycle counts (capped by the hardware maxOffTime), ",
                "and segments shorter than the min-length are discarded. ",
                "The cycle time is auto-detected. Segment 1 of track 42 becomes '42_s1', '42_s2', etc.; ",
                "length and duration filters then apply per segment.")
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
                "Peels localizations from each end of a track (or segment) while they are ",
                "reached across a dark gap ≥ the threshold, or are dim (EFO below the fraction ",
                "of the track's median) — until a solid, well-connected localization is reached. ",
                "Targets the photobleaching tail that inflates Rmax / MSD tails and fakes ",
                "'exports at track ends'. Idempotent; runs after gap-splitting, before quality.")
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
                "Each step is scored by how many acquisition cycles it spans, weighted ",
                "2^(−dark/half-life): 1 = integrated, and a few 1–2-cycle gaps barely hurt ",
                "while long dark gaps do. ", tags$b("Tolerance"), " sets the half-life — higher ",
                "forgives more dark cycles. The per-track mean is the continuity score; tracks ",
                "below the threshold are dropped. Runs after gap-splitting/trim, before metadata.")
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
            helpText("Exports the post-filtered tracks (canonical MINFLUX .csv). ",
                     "Load it in the Simulator → Diffusion/Tracking → ",
                     tags$i("Calibrate from real MINFLUX data"),
                     " to build a ground-truth simulation matched to this dataset."),
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
          uiOutput("trk_view_mode_ui"),
          uiOutput("trk_view_plane_ui"),
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
              "Colors each localization blue (slow/confined) → red (fast/directed) by ",
              tags$b("drift velocity"), " = |net displacement| / span over a centred ",
              "window W. Per-step speed (|Δr|/Δt) is deliberately NOT used: Brownian ",
              "motion has no instantaneous velocity, so |Δr|/Δt diverges as Δt→0 and ",
              "measures your sampling rate, not the molecule. Grey = window too short ",
              "to measure (track ends). Scale maxes at the 95th percentile. ",
              tags$b("W must be well below a typical track duration"), " — if W exceeds ",
              "a track's span, every point in it shares one window and the track goes ",
              "one flat color. ",
              tags$b("Caveat:"), " plain diffusion still yields an apparent ",
              "v ≈ sqrt(4D/W) and noise a floor ≈ sqrt(2)·σ/W — so only sustained red ",
              "well above those baselines indicates real directed transport. Larger W ",
              "suppresses both baselines.")
          ),
          conditionalPanel(
            condition = "input.trk_color_mode == 'dark'",
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Dark-cycle mode colors each localization blue when integrated (a fix every ",
            "cycle) through a red gradient by how many dark cycles preceded it ",
            "(cycle time from Acquisition Parameters). Full-tracks view shows a colorbar; ",
            "in moving-trails playback the colored dots overlay the trails.")
          ),
          radioButtons("trk_dark_show", "Show localizations",
                       choices = c("All" = "all", "Integrated only" = "integrated",
                                   "Dark-cycle only" = "dark"),
                       selected = "all"),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Filters which localizations are drawn by integration status. Integrated = a ",
            "fix every cycle (n_dark = 0); dark-cycle = one or more missed cycles before it. ",
            "Markers only — trails keep the full path for context."),
          conditionalPanel(
            condition = "input.trk_view_mode_main == 'playback'",
            checkboxInput("trk_show_full_tracks_faint", "Show full tracks faintly", TRUE),
            helpText("Markers & (optional) static lines use one color per TID.")
          ),
          checkboxInput("trk_full_raw", "Show ALL raw tracks faintly (background)", FALSE),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "Adds a faint grey overlay of every uploaded (unfiltered) track behind the ",
            "display — in both Full-tracks and Playback modes — so you can see what the ",
            "filters removed."),
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

      # ── Data Quality & Calibration (do this first) ───────────────────────────
      fluidRow(
        column(4,
          h4("Data Quality & Calibration"),
          .method_box("green",
            tagList("Covariance-based estimator, CVE (Vestergaard, Blainey & Flyvbjerg 2014). Reads ",
                    sigloc(), " from the lag-1 anti-correlation of consecutive steps — no MSD fit needed."),
            HTML("Motion blur R=0 for MINFLUX (no exposure integration), so &sigma;<sub>loc</sub>² = −⟨Δxᵢ·Δxᵢ₊₁⟩. Restricted to consecutive <b>single-cycle (integrated)</b> pairs so irregular dark-cycle Δt cannot bias the covariance. Optionally inverse-variance-weighted by the per-localization photon count (ECO), reported as &sigma;<sub>loc-wECO</sub>. This is the app-wide &sigma;<sub>loc</sub> source."),
            tagList("Needs enough integrated pairs; strongly directed motion (drift) makes ",
                    sigloc(), " unresolvable (flagged degenerate — enter it manually then).")),
          tags$div(class = "well well-sm", style = "font-size:0.84em; color:#555; margin-bottom:8px;",
            tags$b("Do this first."), " Estimates the localization precision ", sigloc(), " directly from ",
            "the short-lag step covariance (CVE; ", tags$b("no MSD / diffusion fit needed"),
            ") and scores the data by its integrated/dark-cycle balance. ", sigloc(), " becomes the ",
            "app-wide noise floor for every analysis; the cycle time drives gap-splitting, the ",
            "track-quality gate and the D-map."),
          actionButton("trk_cve_run", HTML("Estimate &sigma;<sub>loc</sub> &amp; quality"),
                       icon = icon("ruler-combined"),
                       class = "btn-primary", style = "margin-bottom:6px; width:100%;"),
          helpText(style = "font-size:0.82em; color:#666; margin-top:0;",
            "Uses consecutive INTEGRATED (single-cycle) steps only, so dark cycles don't bias it. ",
            "Set the cycle time / manual ", sigloc(), " override in the Display tab's Acquisition Parameters. ",
            "Strong directed motion (drift) makes ", sigloc(), " unresolvable — enter it manually then."),
          tags$hr(style = "margin:8px 0;"),
          tags$b(HTML("Per-track &sigma;<sub>loc</sub> heterogeneity"), style = "font-size:0.9em;"),
          fluidRow(
            column(6, numericInput("trk_sig_min_integ", "Min integrated steps / track",
                                   value = 50, min = 10, step = 10)),
            column(6, numericInput("trk_sig_skip", "Skip first N localizations",
                                   value = 5, min = 0, step = 1))
          ),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
            "One CVE ", sigloc(), " per long track; the spread shows whether a single global value is ",
            "safe. The first few localizations carry MINFLUX lock-on transient and are skipped. ",
            "The app-wide ", sigloc(), " uses the per-track median; pooled is shown for cross-check."),
        ),
        column(8,
          uiOutput("trk_cve_status"),
          plotlyOutput("trk_sigma_pertrack_plot", height = "220px"),
          tags$details(style = "margin-top:6px;",
            tags$summary(style = "cursor:pointer; font-size:0.86em; color:#555;",
                         "Data-selection funnel (what's excluded, for your methods)"),
            tags$div(style = "padding:4px 0 0 4px;", uiOutput("trk_quality_funnel"))),
          tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
            tags$b("How this feeds every tool:"), tags$br(),
            sigloc(), " → noise-floor subtraction in the D-map, the JDD noise line, the confinement Rg ",
            "correction, the HMM noise-state flag, and the directionality & speed floors. ",
            "Cycle time → gap-split threshold, Integrated/Dark classification, D-map integrated filter. ",
            "Integrated fraction → the track-quality gate (Display tab). ",
            "The full Δt / Integrated-Dark breakdown is in the TSA section below.")
        )
      ),

      tags$hr(),

      # ── Time Signature Analysis (TSA) ────────────────────────────────────────
      fluidRow(
        column(4,
          h4("Time Signature Analysis (TSA)"),
          .method_box("green",
            HTML("Time-Signature Analysis for MINFLUX (MFX-TSA), after <b>Vogler et&nbsp;al. 2025</b> (<i>Commun. Biol.</i>; reference implementation <i>Toolbox_for_MINFLUX_SPT / the_road_to_minflux_msd</i>). Disentangles integrated vs dark cycles, recovers the atomic cycle time, and suggests gap-split thresholds."),
            HTML("Faithful R port of MFX-TSA: cycle-time / hardware-parameter (EOD, software overhead, pattern vertices) recovery, the integrated/dark cycle balance and per-track cycle progression (cf. <code>print_cycle_stats</code> / <code>plot_cycle_progression_for_ID</code>), and adaptive gap removal via TOT-based thresholds with a 50-localization minimum (cf. <code>suggest_split_thresholds</code> / <code>remove_time_gaps</code>). This is the MINFLUX-native layer every estimator below depends on — it builds the linearized cycle-index timeline that makes irregular Δt tractable."),
            "Descriptive, not inferential. The cycle time is an estimate; cross-check against the sequence JSON / hardware parameters when available."),
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
            checkboxInput("trk_dt_log_y", "Log count axis", value = FALSE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "One bar per cycle count: bin edges sit on the classification boundaries, ",
              "SW + (k ± ½) × cycle, so each bar is exactly one cycle count centred on its ",
              "peak — the same round((Δt − SW)/cycle) the analysis uses. Counts fall by ",
              "roughly two orders of magnitude from the 1-cycle peak into the dark tail, ",
              "which is flat against a linear axis. The median and bin-free typical Δt are ",
              "marked with vertical lines."),
            checkboxInput("trk_dt_mark_cycles", "Mark cycle multiples (1×, 2×, …)", value = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Dotted verticals at Δt = SW + k × cycle, anchored on the observed 1-cycle peak.",
              tags$br(),
              tags$b("Δt is hard-quantized"), ", at Δt = k × cycle + SW"),
            tags$hr(style = "margin:8px 0;"),
            numericInput("trk_cycle_dark_max", "Max cycles to display", value = 10, min = 4, step = 1),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Each consecutive interval is classified by how many acquisition cycles it ",
              "spans (round((Δt − offset) / cycle time)). 1 cycle = Integrated (clean); ",
              "≥2 cycles = Dark cycles preceded it (molecule briefly unseen)."),
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

      # ── Sample geometry (shared) ──────────────────────────────────────────────
      # Defined once and read by every tool that needs it: the MSD's projection
      # factor and free-on-geometry overlay, and the displacement test's Step-2 null.
      # Keeping one definition means the shape can never be set two different ways.
      fluidRow(
        column(4,
          h4("Sample geometry"),
          .method_box("green",
            "The shape the molecules move in or on, and how many directions they are free to move in. Read by every analysis that needs to know what free motion would look like in this sample.",
            HTML("The shape is always <b>three-dimensional</b> — a specimen is 3D whether or not the tracking is — and is projected into the image plane automatically for 2D datasets. It is used by the MSD (the factor that converts the measured curve into D, and the optional free-diffusion-on-this-shape overlay) and by the displacement test's geometry-aware band, so it is set here once rather than separately in each panel."),
            HTML("Both settings matter only when motion is restricted. For molecules moving freely in three dimensions the shape has no effect until they begin reaching its walls. Choose <b>None</b> when there is no bounding structure — molecules in solution, or a container far larger than anything the molecule explores. The panel opposite sets out exactly how each choice enters the calculation.")),
          selectInput("trk_geom_type", "Shape",
            choices = c("None — unbounded"                          = "none",
                        "Sphere / round cell"                       = "sphere",
                        "Rod / spherocylinder (e.g. a bacterium)"   = "rod",
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
            choices = c("Through a volume — 3 degrees of freedom"       = "volume",
                        "On a surface — 2 degrees of freedom"            = "surface",
                        "Along a filament / track — 1 degree of freedom" = "filament"),
            selected = "volume"),
          helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
            "Degrees of freedom are a property of the motion, not of the file: a molecule ",
            "held in a membrane moves in two directions whether two or three coordinates ",
            "were recorded, and a motor on a filament moves in one. Both the MSD and the ",
            "displacement test read this setting.")
        ),
        column(8,
          tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
            tags$b("How these settings enter the calculation"), tags$br(),
            "The shape acts in two independent places, on different timescales.",
            tags$br(), tags$br(),

            tags$b("1 · The factor relating MSD to D"), tags$br(),
            HTML("Applies at every lag, as soon as motion is restricted, and does <b>not</b> depend on the shape's size:"),
            tags$div(style = "margin:3px 0 3px 10px;", tags$code("MSD = 2 · d_eff · D · τ")),
            HTML("With all three coordinates recorded, d<sub>eff</sub> is simply the number of directions the molecule can move in — 3 for a volume, 2 for a surface, 1 for a filament."),
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
              tags$tr(tags$td("filament in plane"), tags$td("—"), tags$td("1"), tags$td(tags$code("2Dτ")))),
            HTML("A membrane tilted by θ gives 1 + cos²θ, between the rod and flat values. Reading restricted motion as free motion under-reports D by a third on a sphere and a quarter on a rod."),
            tags$br(), tags$br(),

            tags$b("2 · The ceiling the MSD approaches"), tags$br(),
            HTML("Applies only once the molecule has had time to cover the shape, and this is where size enters. On a closed shape the MSD cannot exceed the mean squared separation of two random points on it:"),
            tags$table(style = "margin:5px 0 5px 10px; border-collapse:collapse;",
              tags$tr(style = "color:#666;",
                tags$th(style = "text-align:left; padding-right:14px;", "shape"),
                tags$th(style = "text-align:left; padding-right:14px;", "3D"),
                tags$th(style = "text-align:left;", "2D projection")),
              tags$tr(tags$td("sphere, radius R"), tags$td(tags$code("2R²")), tags$td(tags$code("4R²/3"))),
              tags$tr(tags$td("rod, radius r, length L"), tags$td(tags$code("L²/6 + 2r²")), tags$td(tags$code("L²/6 + r²")))),
            HTML("These are squared lengths: size determines where the curve levels off, not how steeply it rises. Both expressions agree with simulation to better than 1%."),
            tags$br(), tags$br(),

            tags$b("Why the ceiling matters for interpretation"), tags$br(),
            HTML("Because a closed shape flattens the MSD on its own, a molecule that is diffusing perfectly freely can still produce a curve that the non-linear models read as biology. Simulating free diffusion on a 500&nbsp;nm sphere over a long fit window and fitting each model returns a confined size of L<sub>c</sub>&nbsp;≈&nbsp;1000&nbsp;nm — a converged fit to a corral that does not exist — and an anomalous exponent α&nbsp;≈&nbsp;0.86, sub-diffusive with nothing anomalous present. The directed model correctly returns v&nbsp;≈&nbsp;0, so geometry does not manufacture drift."),
            tags$br(),
            HTML("Enabling the overlay in the MSD panel fits whichever model is selected to the simulated free curve, so a fitted parameter can be compared against what the shape alone would produce. The effect grows with the fit window: over a short window the molecule covers far less than the shape and the correction is negligible."))
        )
      ),

      tags$hr(),

      # ── MSD ─────────────────────────────────────────────────────────────────
      fluidRow(
        column(4,
          h4("MSD (diffusion)"),
          .method_box("green",
            "Mean-squared-displacement with optimal least-squares fitting (Michalet 2010; Michalet & Berglund 2012); cycle-based binning and the D-ceiling from Vogler et al. 2025",
            HTML("Binned by <b>integer cycle count = true elapsed time</b>, not step index, so cycle-skipping doesn't corrupt the curve. OLSF + Michalet optimal-point selection, fit window 5–50 cycles, non-diffusing (D≤0) tracks dropped, eDIFF (ensemble) and tDIFF (per-track) both reported."),
            HTML("Optimal-point theory is <b>linear/free-diffusion only</b> — the Confined/Directed/Anomalous models are flagged exploratory. R=0 makes D under-estimated above the ~2–3 µm²/s trackable ceiling (D<sub>max</sub> warning).")),
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
            numericInput("trk_msd_min_pts",  "Min localizations per TID (skip short)", value = 50, min = 5, step = 5),
            radioButtons("trk_msd_ndim", "Coordinates used",
              choices  = c("Auto-detect" = "auto", "2D  (x, y)" = "2", "3D  (x, y, z)" = "3"),
              selected = "auto", inline = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              tags$b("Auto:"), " uses 3D if z has real variance (SD > 1 nm), else 2D. ",
              "This is which coordinates go into the MSD — how many directions the ",
              "molecule can actually move in is the separate setting below."),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "How many directions the molecule can move in, and the shape it moves on, ",
              "are set in ", tags$b("Sample geometry"), " above and applied here ",
              "automatically. The MSD-to-D relation actually in use is printed with the ",
              "results whenever it differs from the plain one."),
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
            # THE MSD — always on. The checkbox is kept hidden (constant TRUE) only so the
            # server's isTRUE(input$trk_msd_unified) gates and conditionalPanels keep working.
            div(style = "display:none;", checkboxInput("trk_msd_unified", "MSD", TRUE)),
            div(
              tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin:4px 0;",
                tags$b("MINFLUX MSD (cycle-based)."), " Cumulative cycle-index binning, ",
                "track-weighted ensemble (eMSD) on the linearized lag×cycle timeline, iterative ",
                "OLSF that converges on the lowest lags. Reports both eDIFF (ensemble) and tDIFF ",
                "(per-track mean). Method follows Vogler et al. (2025, ", tags$em("Commun. Biol."), ")."),
              fluidRow(
                column(6, numericInput("trk_uni_min_lag",
                             "Lower cutoff (cycles)", value = 5, min = 1, step = 1)),
                column(6, numericInput("trk_uni_max_lag",
                             "Upper fit limit (cycles)", value = 50, min = 10, max = 5000, step = 5))
              ),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Fit window = [lower, upper] cycle-lags (default 5–50, per Vogler 2025). ",
                tags$b("Low lags"), " are dominated by localization error (", sigloc(), ") and secondary/linker ",
                "motion; ", tags$b("high lags"), " have few pairs and huge variance (Qian 1991) — so the ",
                "fit uses the middle 'diffusion regime'. OLSF + Michalet pick the optimal point count within it. ",
                tags$b("The exact lower bound is probe-specific"), " (not a universal constant) — inspect the ",
                "MSD curve for where the noisy low-lag kink ends. ",
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
              conditionalPanel(
                condition = "input.trk_uni_model != 'linear'",
                div(style = "margin-top:4px;padding:5px 8px;border-left:3px solid #e67e22;background:rgba(230,126,34,0.06);font-size:11px;color:#7f8c8d;line-height:1.4;",
                    HTML(paste0(
                      "<b style='color:#e67e22;'>Exploratory.</b> Confined / directed / anomalous MSD ",
                      "combined with MINFLUX irregular cycle-time sampling is beyond the standard, ",
                      "literature-validated free-diffusion workflow (Vogler et&nbsp;al. 2025 build the ",
                      "irregular-Δt machinery but apply it only to free diffusion; their toolbox flags ",
                      "anomalous/confined as “implement on your own”). The model forms are standard ",
                      "(Kusumi et&nbsp;al. 1993 corral; Metzler et&nbsp;al. 2014 power-law) but their fit to short, ",
                      "cycle-skipped MINFLUX tracks is not validated — treat L_c, v and α as exploratory. ",
                      "See docs/MINFLUX_analysis_methods.md §3a.")))),
              checkboxInput("trk_msd_geom",
                          HTML("Overlay <b>free diffusion on the cell geometry</b>"), FALSE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "On a closed cell the MSD of a ", tags$b("free"), " molecule bends and eventually ",
              "saturates on its own, because there is nowhere further to go. A bend is therefore ",
              "not by itself evidence of trapping. This overlays the curve free diffusion would ",
              "produce on the chosen shape at the fitted D, on the same lags and coordinates as ",
              "the data: measured curve on the green line = the cell explains the bend, measured ",
              "curve below it = confinement beyond the cell.",
              tags$br(),
              tags$b("Matters mainly for long fit windows."), " Over a short window the molecule ",
              "covers far less than the cell and the two are indistinguishable; the gap grows as ",
              "the upper fit limit is raised. The estimates line reports how much the geometry ",
              "alone accounts for over your current window.",
              tags$br(),
              "The shape and its dimensions are taken from ", tags$b("Displacement vs Diffusion"),
              " below, so the cell is defined once and used by both tests. Surface vs volume ",
              "follows the ", tags$b("Molecule moves"), " setting above."),
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
              checkboxInput("trk_uni_bootstrap",
                            "Bootstrap 95% CI on eDIFF (slower)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Resamples tracks with replacement (~300×) and refits — the honest CI given ",
                "within-track correlation. Adds a few seconds; recomputes on Plot MSD."),
              checkboxInput("trk_uni_show_band", "Show IQR band", TRUE),
              checkboxInput("trk_uni_show_ntau", "Show N(τ)", TRUE),
              checkboxInput("trk_uni_show_spag", "Show individual track TAMSDs (spaghetti)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Faint per-track MSD curves. When 'Cluster tracks by D & CI' is on, they are ",
                "coloured by cluster and a per-cluster eMSD + tDIFF is added."),
              checkboxInput("trk_msd_show_bindiag",
                            "Show cycle-binning diagnostic (why the tail wobbles)", FALSE),
              helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
                "Re-inspects the cycle binning: pairs per cycle-lag (where counts collapse the ",
                "MSD gets noisy / saw-toothed) and the ratio of real elapsed time to lag × cycle ",
                "(should be flat — a rise or wobble means the cycle assignment drifts, e.g. under ",
                "long dark gaps). Reads the current MSD result; no recompute.")
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
            # MSD sub-panel (cycle-based, always shown when MSD is enabled)
            conditionalPanel(
              condition = "input.trk_msd_unified == true",
              plotlyOutput("trk_msd_unified_plot", height = "300px"),
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
              tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin-bottom:8px;",
                tags$b("What this does:"), tags$br(),
                HTML("Estimates the maximum trackable diffusion coefficient D<sub>max</sub> from acquisition parameters. "),
                "Based on Vogler et al. (2025, ", tags$em("Commun. Biol."), ")."
              ),
              fluidRow(
                column(6, numericInput("trk_dmax_L",   "TCP diameter L_TCP (nm)",         value = 40,  min = 10, max = 2000, step = 10)),
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
              uiOutput("trk_msd_cluster_model_ui"),
              plotlyOutput("trk_msd_cluster_curves", height = "320px"),
              uiOutput("trk_msd_cluster_models"),
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

      # ── Confinement (MSD confinement index) ───────────────────────────────────
      # The Rg/Rmax shape test that used to live here was retired 2026-08-03; its UI
      # and server code are in retired_ui.R / retired_tracking_server.R (gitignored,
      # not sourced).
      fluidRow(
        column(4,
          h4("Confinement (MSD index)"),
          .method_box("green",
            "Per-track confinement from the shape of each track's own MSD curve: CI = (MSD_late / MSD_early) / (τ_late / τ_early), compared against a free-diffusion null simulated at the dataset's measured conditions.",
            HTML("A freely diffusing molecule has an MSD that keeps growing in proportion to time, so its late lags exceed its early lags by the same factor as the lag times do. A confined molecule runs out of room: its MSD flattens towards a plateau, so the late rise falls short. A molecule under directed transport does the opposite. CI captures that as one number per track.<br>Because it compares a track against <b>itself</b> — its own short lags versus its own long lags — it needs no externally supplied diffusion coefficient, and a molecule being intrinsically fast or slow does not by itself move CI. It is not circular: short-lag and long-lag behaviour are independent pieces of information. The MSD it reads is the cycle-binned one, so irregular Δt is already accounted for.<br><b>Why the free range is simulated rather than fixed.</b> For ideal noise-free motion free diffusion gives CI&nbsp;=&nbsp;1, but real data never does: localization error adds a constant 2·d·&sigma;<sub>loc</sub>² to every MSD point, which inflates the early lags proportionally more than the late ones and pushes free CI below 1. How far below depends on the noise, the mobility, the cycle time, the fit window and the track length, so no fixed threshold can be right for every dataset. The free range is therefore generated by simulating free tracks under the measured &sigma;<sub>loc</sub>, D, cycle time and dark-gap pattern, and is matched to each track's localization count — and to its own D when the two interact (see below). <b>Confined</b> = below the free 5th percentile, <b>directed</b> = above the 95th."),
            HTML("<b>There is a hard resolution floor.</b> Confinement can only be seen if the molecule would otherwise have travelled further than the boundary it is trapped by. The distance a free molecule covers in one track, √(2·d·D·T), is reported with the result: a domain much larger than that is invisible, so a low confined&nbsp;% means \"nothing confined below this scale\", not \"everything is free\". As a rule of thumb the free extent needs to exceed the domain size by a factor of about two before per-track calls carry information.<br><b>Confined&nbsp;% cannot be compared between datasets with different track lengths or mobilities</b>, because detection power differs. Longer tracks and faster molecules see more confinement for the same biology. Always quote the resolution floor next to the percentage.<br><b>Weak for directed motion.</b> A modest drift barely lifts the late-versus-early ratio, so the directed tail is insensitive; a near-zero directed&nbsp;% is not evidence against transport. Use Directionality or Speed&nbsp;&amp;&nbsp;Transport for that question.<br><b>Other limits.</b> Each tail contains about 5% of purely free tracks by construction, so small percentages are the floor rather than a signal. Tracks too short to fill the MSD fit window are dropped. The test assumes a single behaviour per track — a molecule that is confined for part of its trajectory and free for the rest returns a blend; use Change-Point or the HMM panel to split those first.")),
          checkboxInput("trk_show_conf", "Enable confinement index", FALSE),
          conditionalPanel(
            condition = "input.trk_show_conf == true",
            helpText(style = "font-size:0.82em; color:#666;",
              "Requires the MSD to have been run (Plot MSD) — the confinement index is ",
              "computed per track from the same cycle-binned MSD. ", sigloc(), " and D come ",
              "from the Data Quality & Calibration panel."),
            numericInput("trk_conf_ci_nsim", "Null simulations per length", value = 200,
                         min = 50, max = 500, step = 50),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Free tracks simulated per localization-count grid point. Higher = a smoother ",
              "null band at the cost of run time."),
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
          h4("Track Segmentation — Change-Point + GMM"),
          .method_box("green",
            "PELT change-point detection (Killick, Fearnhead & Eckley 2012) on the local-D signal, then a Gaussian mixture (BIC-selected) over the per-segment diffusivities.",
            HTML("Segments on <b>log(r²/Δt)</b> — i.e. local D — so the change-point cost is already Δt-aware and robust to irregular sampling. Per-segment D = mean(r²/Δt)/(2d) (not mean(r²)/mean(Δt)), which avoids the Jensen bias when r² and Δt covary."),
            tagList("Single-step D is noisy, so PELT can over-segment (the min-segment length guards this). The ",
                    sigloc(), " noise floor σ²·⟨1/Δt⟩ is subtracted from each segment's D (toggle below) so immobile segments read ~0 rather than a small non-zero D.")),
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
            checkboxInput("trk_cpd_sub_noise",
                          HTML("Subtract localization-noise floor (&sigma;<sub>loc</sub>) from segment D"), TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Removes the σ²·⟨1/Δt⟩ apparent-diffusion offset so immobile segments read ~0, not a ",
              "small non-zero D. Needs ", sigloc(), " set (Data Quality & Calibration panel, or Acquisition Parameters)."),
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
            HTML("Two ways to handle irregular Δt: (1) normalize each jump by its own <b>√Δt</b> so steps of different cycle counts pool together, or (2) <b>Integrated steps only</b> — keep single-cycle steps so every jump is at the same true Δt and no pooling assumption is needed (the rigorous mode). A &sigma;<sub>loc</sub> bias term removes localization-error broadening (small for MINFLUX). MINFLUX's high localization count means plenty of single-cycle jumps, so the classic 'short track → few jumps' concern is weak here."),
            "The √Δt-normalized mode assumes free diffusion (distorts under confinement/anomaly) — prefer Integrated-only when in doubt. Keep the component count low; a Rayleigh mixture still assumes each state is free Brownian."),
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
                          "Integrated steps only (single-cycle Δt — rigorous)", FALSE),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
              "Keeps only single-cycle (no dark gap) steps, so every jump is at the same true Δt. ",
              "The √Δt-pooling assumption then no longer applies — the distribution is exact regardless ",
              "of the motion model. MINFLUX tracks carry many localizations, so there are usually plenty ",
              "of single-cycle steps. Needs the cycle time set."),
            tags$hr(style = "margin:6px 0;"),
            numericInput("trk_jdd_min_dt", "Min Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 0.1),
            numericInput("trk_jdd_max_dt", "Max Δt (ms, 0 = off)",
                         value = 0, min = 0, step = 1),
            numericInput("trk_jdd_sigma_loc", HTML("localization error &sigma;<sub>loc</sub> (nm, 0 = auto: CVE calibration, MSD fallback)"),
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
                tags$li("Get ", sigloc(), " from the ", tags$b("Data Quality & Calibration"), " panel (the CVE estimate, which needs no MSD fit), or read it below the MSD plot. This is the random positional error per measurement: a smaller number means more accurate positions. Leaving the ", tags$b(sigloc(), " correction"), " field at 0 already uses that app-wide value; set it only to override. The MSD estimates panel also shows this as a blue hint.")
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

      # ── Displacement vs Diffusion (Step 1 → Step 2 geometry-aware null) ───────
      fluidRow(
        column(4,
          h4("Displacement vs Diffusion"),
          .method_box("orange",
            "MINFLUXR original: relates each track's NET start→end displacement to the free-diffusion expectation √(4·D·T) at the ensemble D (CVE).",
            HTML("For 2D Brownian motion ⟨R<sub>net</sub>²⟩ = 4·D·T and R<sub>net</sub> is Rayleigh-distributed (scale √(2DT)), so each track is judged against that spread, not a point. Measured in x,y (the 2D image projection). Endpoint &sigma;<sub>loc</sub> noise is subtracted. Path length is shown but never used to classify (it scales with sampling rate, like path speed)."),
            HTML("This is the INFINITE-PLANE null — it does not account for finite/curved geometry, which caps a FREE molecule's displacement too. So 'sub-diffusive' = confined OR geometrically constrained; Step 2 (free diffusion on a chosen geometry) separates them. Per-track D is too noisy, so the ENSEMBLE D is used — run Data Quality &amp; Calibration first.<br><b>Consequence of the single ensemble D:</b> in a population with a wide spread of per-track D, slower-than-average molecules fall below the band even when they diffuse freely, inflating the sub-diffusive fraction — the panel warns when that spread is large. Each tail also holds ~5% of free tracks by construction.<br><b>Overlap with other panels:</b> R<sub>net</sub> depends only on the two endpoints and the true elapsed time, so it is robust to irregular Δt, but it is D-dependent. <b>Confinement Analysis</b> is the D-independent shape test. Expect different percentages; treat this one as an upper bound on confinement.")),
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
              "Degrees of freedom and the sample shape come from ", tags$b("Sample geometry"),
              " at the top of this tab and are applied here automatically."),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
              "3D uses a Maxwell (√6DT) free null; 2D uses Rayleigh (√4DT). Even with 2D tracking the ",
              "underlying structure is 3D, so the Step-2 geometry null is simulated in 3D and projected ",
              "to the image plane — pick 2D here and a curved geometry below to test that."),
            selectInput("trk_disp_dsrc", "Diffusion coefficient D (sets the null)",
              choices = c("CVE — short-lag, noise-corrected (recommended)" = "cve",
                          "MSD fit (current MSD model)" = "msd",
                          "Manual value" = "manual"),
              selected = "cve"),
            conditionalPanel("input.trk_disp_dsrc == 'manual'",
              numericInput("trk_disp_dman", "D (µm²/s)", value = 0.1, min = 0, step = 0.01)),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
              "The null is √(2·d·D·T), so ", tags$b("D sets the verdict"), " — a larger D raises the ",
              "band and flags more tracks sub-diffusive. It must be the ", tags$b("short-time"),
              " mobility. ", tags$b("CVE"), " is recommended: short-lag, noise-corrected and model-free, ",
              "so it is measured before confinement bends the MSD. ", tags$b("MSD fit"),
              " is fitted over the lag window — if the population is confined the MSD has already bent ",
              "inside that window, so this D is deflated and the test ", tags$b("under-calls"),
              " confinement. ", tags$b("Manual"), " is for a known/literature value or sensitivity testing. ",
              "The R", tags$sub("net"), " noise-subtraction uses the app-wide ", sigloc(), "."),
            actionButton("trk_disp_run", "Run displacement", icon = icon("play"),
                         class = "btn-primary btn-sm", style = "margin:4px 0; width:100%;"),
            tags$hr(style = "margin:8px 0;"),
            checkboxInput("trk_disp_geom",
                          HTML("<b>Step 2:</b> compare to free diffusion on a geometry"), FALSE),
            conditionalPanel(
              condition = "input.trk_disp_geom == true",
              uiOutput("trk_disp_geom_echo"),
              helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
                "Shape and dimensions come from ", tags$b("Sample geometry"), " at the top of ",
                "this tab, and surface vs volume from ", tags$b("Molecule moves"), " in the MSD ",
                "panel, so the sample is described once for every analysis. A surface constrains ",
                "motion to a curved 2D sheet and foreshortens under projection; an interior only ",
                "caps displacement at the enclosure size. With the shape set to ",
                tags$b("None"), " this band is the same unbounded null as Step 1."),
              helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
                "Adds the ", tags$b(style = "color:#1B9E5A;", "green band"), ": free diffusion on the ",
                "chosen 2D-projected geometry at your D (same shapes and physics as the Diffusion ",
                "Simulator), matched to your track durations. Points are re-classified against it — ",
                "tracks below green are ", tags$b("truly confined"), " beyond geometry; between green and ",
                "grey is ", tags$b("geometrically constrained but free"), ". Closed surfaces (sphere / rod / ",
                "ellipsoid) impose a size limit; membrane / free are unbounded and match the plane. ",
                "Cycle time comes from Acquisition Parameters.")
            ),
            tags$br(),
            downloadButton("trk_dl_disp_csv", "Download displacement table (CSV)"),
            helpText(style = "font-size:0.8em; color:#666; margin-top:4px;",
              "CSV has one row per track (R", tags$sub("net"), ", R", tags$sub("max"),
              ", path length, duration, ρ, class) plus the ensemble D and ", sigloc(),
              " — also the exact input for the Step 2 geometry test.")
          )
        ),
        column(8,
          conditionalPanel(
            condition = "input.trk_show_disp == true",
            plotlyOutput("trk_disp_plot", height = "380px"),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("How to read this plot:"), tags$br(),
              "Each dot is one track: X = duration, Y = how far it got from start to end ",
              "(R", tags$sub("net"), "). The grey band is where FREE diffusion at the measured D ",
              "lands 90% of the time (5–95%), and the dashed line its median. ",
              tags$b(style = "color:#2c7fb8;", "Blue = below the band"), " (sub-diffusive: confined ",
              "or geometrically constrained), ", tags$b(style = "color:#9e9e9e;", "grey = within"),
              " (free-consistent), ", tags$b(style = "color:#e31a1c;", "red = above"),
              " (super-diffusive: directed). A free population puts ~5% in each tail — a large ",
              "blue or red excess is the signal. ",
              tags$b("The blue call is not final:"), " free diffusion on a curved/finite surface is itself ",
              "capped by geometry, so Step 2 is needed to tell real confinement from curvature. ",
              "With Step 2 on, points are re-classified against the ", tags$b(style = "color:#1B9E5A;", "green"),
              " (free-on-geometry) band instead."),
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
            "Turn-angle distribution and velocity autocorrelation vs. the isotropic-Brownian null (Burov et al. 2013, distribution of directional change).",
            HTML("The turn-angle null is Δt-independent (the angle between two isotropic steps is uniform regardless of their durations). The velocity autocorrelation is binned by <b>elapsed time</b> (cycle lag × cycle time), not step index, and both use only integrated-ish steps (≤ max dark cycles). MINFLUX's small &sigma;<sub>loc</sub> removes most of the classic noise anti-persistence artefact."),
            tagList("Single short tracks lack power — pool across tracks. A small negative first-lag value is the expected ",
                    sigloc(), " signature for a diffusing molecule, not evidence of confinement.")),
          checkboxInput("trk_show_dir", "Enable directionality analysis", FALSE),
          conditionalPanel(
            condition = "input.trk_show_dir == true",
            numericInput("trk_dir_min_pts", "Min points per TID", value = 6, min = 4, step = 1),
            numericInput("trk_dir_max_dark", "Max dark cycles per step", value = 2, min = 0, step = 1),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Steps spanning more than this many skipped (dark) cycles are excluded from the ",
              "turn-angle and autocorrelation. 1–2 dark cycles are tolerable; long gaps break the ",
              "irregular-Δt link and may hide a track merged from another molecule. ",
              "The autocorrelation lag axis is real elapsed time (cycles × cycle time), not step count, ",
              "so irregular sampling is not misrepresented as evenly spaced."),
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
              "C(τ): if the molecule moved in direction X, is it still moving in direction X after a time τ? ",
              "The lag axis is ", tags$b("real elapsed time"), " (cycle count × cycle time), not step ",
              "count, so irregularly-spaced localizations are placed at their true time separation. ",
              tags$b("Positive and slowly decaying"), " = persistent, directed motion. ",
              tags$b("Negative at the first lag and near zero after"), " = consistent with free Brownian diffusion. ",
              tags$b("Important caveat:"), " localization noise introduces an ",
              "inherent anti-correlation at the shortest lag even for purely Brownian motion ",
              "(C ≈ −σ²_loc / MSD₁cycle). For MINFLUX data a small negative first-lag value is expected ",
              "for diffusing molecules and does not by itself indicate confinement or active transport. ",
              "Only a first-lag value substantially more negative than the next lags is informative."
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
          .method_box("green",
            "Per-track Gaussian/Rayleigh-emission hidden Markov model, Baum–Welch + Viterbi (Das, Cairo & Coombs 2009; short-track pooling in the spirit of vbSPT, Persson et al. 2013).",
            HTML("Emission steps are <b>r/√Δt</b>-normalised so states correspond to diffusivities, not raw step sizes, despite variable Δt. Transition <b>rates are continuous-time</b>: a generator Q estimated from the decoded path using real elapsed time gives switching rates (s⁻¹) and mean dwell times that are invariant to dark cycles."),
            tagList("The Viterbi decode still uses the per-step transition prior (emission-dominated, so state assignments stay robust); states at or below the ",
                    sigloc(), " noise floor are flagged noise/immobile-limited.")),
          checkboxInput("trk_show_hmm", "Fit HMM to tracks", FALSE),
          conditionalPanel(
            condition = "input.trk_show_hmm == true",
            numericInput("trk_hmm_n_states",    "Number of states",                value = 2,   min = 2, max = 4,   step = 1),
            numericInput("trk_hmm_min_pts",     "Min localizations per TID",       value = 10,  min = 4, step = 1),
            numericInput("trk_hmm_n_ex",        "Example tracks to display",       value = 6,   min = 1, max = 30,  step = 1),
            numericInput("trk_hmm_max_fit_trk", "Max tracks used for fitting",     value = 100, min = 10, max = 500, step = 10),
            numericInput("trk_hmm_max_fit_steps", "Max steps per track for fitting (0 = all)",
                         value = 0, min = 0, step = 100),
            checkboxInput("trk_hmm_norm_dt", "Normalize steps by √Δt (recommended)", TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Fitting time is roughly proportional to (tracks × steps × EM iterations). ",
              tags$b("Max tracks"), " keeps the longest N tracks; ", tags$b("Max steps"),
              " truncates each of those to its first N steps for fitting only ",
              "(0 = use all). Both bound runtime with little loss of accuracy — 100 tracks ",
              "and ~300–500 steps are usually ample for stable EM. The per-localization ",
              "state decode always uses the full tracks."),
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
          .method_box("amber",
            "Pools each HMM state's contiguous same-state runs into a per-state ensemble MSD and fits it with the chosen MSD model — the payoff for mixed-motion tracks.",
            "Inherits the MSD's cycle-index binning and fit window. Lets the Confined model fit the slow/wall state (→ confinement size L) and the Directed model the fast/mobile state (→ velocity v), avoiding the whole-track wash-out that mixes phases.",
            "Confined/anomalous MSD on MINFLUX is exploratory (see the MSD model note). State assignment quality caps this — a poor HMM makes the per-state MSD meaningless."),
          tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555; margin-bottom:8px;",
            "Splits each track into HMM diffusion-state segments and computes a separate ",
            "MSD per state, fit with the MSD model + window. Resolves a mixed ",
            "population: fit ", tags$b("Confined"), " on the slow/wall state for a real ",
            "confinement size L, and ", tags$b("Directed"), " on the fast/mobile state for a ",
            "velocity v — instead of one washed-out whole-track MSD. ",
            "Requires the Diffusion State HMM (above) enabled and computed."),
          checkboxInput("trk_show_seg_msd", "Enable per-segment MSD", FALSE),
          conditionalPanel(
            condition = "input.trk_show_seg_msd == true",
            numericInput("trk_seg_msd_min_pts", "Min points per segment", value = 15, min = 5, step = 5),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Model, fit window, linearized-τ and iterative settings are inherited from the ",
              "MSD panel — switch its 'Fit model' to Confined or Directed and re-view here."),
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

      # ── Speed & Transport (exploratory — kept late in the tab) ────────────────
      fluidRow(
        column(4,
          h4("Speed & Transport"),
          .method_box("amber",
            "MINFLUXR original: per-step path speed, a drift velocity averaged over a time window, and a noise-corrected local anomalous exponent α(t). Finds sustained fast segments and scores each against evidence that does not come from the speed itself.",
            tagList("Built directly on the irregular-Δt data: speeds use real Δt; the drift velocity averages net displacement over a window W (µm/s); local α is de-biased with ",
                    sigloc(), "."),
            HTML("<b>Path speed cannot separate diffusion from transport</b> — diffusion inflates instantaneous speed, so a high path speed is not evidence of directed motion. The drift velocity is what flags a segment, and because the flag <i>is</i> a speed threshold the segment's own velocity is not independent evidence. A tracker that has lost its molecule and is searching produces the same gross signature: fast, outward, at the end of a track. That is what the corroboration below is for.")),
          checkboxInput("trk_show_speed", "Enable speed & transport", FALSE),
          conditionalPanel(
            condition = "input.trk_show_speed == true",

            tags$b("Finding drift segments", style = "font-size:0.9em;"),
            numericInput("trk_speed_drift_win", "Drift window W (ms)",
                         value = 10, min = 1, step = 2),
            numericInput("trk_speed_fold", "Threshold (× diffusive baseline)",
                         value = 3, min = 1, step = 0.5),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "A segment is called when the drift velocity stays above ", tags$b("threshold × "),
              "the speed diffusion alone would produce over the same window. Larger W rejects ",
              "more noise but blurs the onset; a higher threshold is stricter."),
            numericInput("trk_speed_min_pts", "Min points per TID", value = 20, min = 5, step = 5),

            tags$details(style = "margin:8px 0;",
              tags$summary(tags$b("Advanced", style = "font-size:0.88em;")),
              tags$div(style = "padding-left:4px;",
                numericInput("trk_speed_drift_minpts", "Min localizations in window",
                             value = 5, min = 2, step = 1),
                numericInput("trk_speed_min_run", "Min sustained windows to call a segment",
                             value = 3, min = 1, step = 1),
                numericInput("trk_speed_edge_guard", "Edge guard (× W) — ignore segments this close to track start/end",
                             value = 1, min = 0, max = 5, step = 0.5),
                numericInput("trk_speed_onset_thr", "Absolute threshold (µm/s, 0 = use × baseline)",
                             value = 0, min = 0, step = 0.1),
                numericInput("trk_speed_sigma_nm",
                             HTML("Localization precision &sigma; (nm, 0 = auto)"),
                             value = 0, min = 0, step = 1),
                helpText(style = "font-size:0.78em; color:#666; margin-top:-4px;",
                  "Sets the noise floor √2σ/W and de-biases α. Auto uses the CVE calibration, ",
                  "falling back to the MSD fit."),
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
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin:6px 0;",
              tags$b("Is a late-track excursion real, or the tracker searching?"), tags$br(),
              "Each dot is one candidate segment: ", tags$b("x"), " = where it starts along the ",
              "track (0 = start, 1 = end), ", tags$b("y"), " = the fraction of its localizations ",
              "that were integrated rather than following dark cycles. A molecule that is ",
              "genuinely present and moving keeps being found, so its event sits ",
              tags$b("high"), ". A tracker that has lost its molecule accumulates dark cycles ",
              "while it searches, so its event sits ", tags$b("low"), " — and typically at the ",
              "far right, because that is where tracks end. ",
              tags$b("Position alone cannot separate the two; the height is the evidence."),
              " Points are coloured by the same verdict as the summary above."),
            fluidRow(
              column(6, plotlyOutput("trk_export_hist",       height = "290px")),
              column(6, plotlyOutput("trk_export_onset_hist", height = "290px"))
            ),
            tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
              tags$b("Left"), " = distribution of the chosen per-track metric (grey = unflagged, ",
              "green = flagged); a separated green mode suggests a distinct subpopulation. ",
              tags$b("Right"), " = where motion is strongest along each track. A pile-up near 1 ",
              "is expected either way and is not confirmation on its own.", tags$br(),
              tags$b("Edge guard"), " ignores the first and last window of each track, since ",
              "track starts carry lock-on jitter that otherwise piles false events at 0. The ",
              tags$b("CSV"), " has one row per track including the event_* corroboration columns. ",
              tags$b("A matched negative control processed identically is the only thing that ",
              "fully settles artifact versus signal.")),
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
                "above the noise floor (orange dashed) and onset threshold (red dashed, points red) = a ",
                "candidate directed segment (not a proven export — check the corroboration fields). ",
                tags$b(style = "color:#7B4FBE;", "Purple = local α"), " → 2 directed, ~1 Brownian. ",
                tags$b(style = "color:#888;", "Grey dotted = path speed"), ", stays high/flat (noise-dominated). ",
                downloadButton("trk_dl_speed_csv", "Per-timepoint curves (CSV)", class = "btn-xs btn-default")
              )
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
            HTML("Single-lag diffusivity per spatial cell, D = median(r²/2dΔt). This is the short-lag limit of the CVE (Vestergaard 2014 / Vogler 2025) — <b>not</b> Ostersehlt 2022, which is a DNA-PAINT imaging paper with no diffusion mapping."),
            tagList("The 'integrated only' option drops dark-cycle steps (whose midpoint is meaningless), and the ",
                    sigloc(), "²/Δt apparent-diffusion floor is subtracted so immobile regions don't read as diffusing. Uses real Δt per step."),
            "Single-lag D ignores confinement/anomaly (it is an instantaneous estimate); cells need enough steps to be reliable (blank otherwise). Spatial heterogeneity of the sampling can imprint on the map."),
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
              "D per cell = median(r²/4Δt) over steps whose midpoint falls in the cell. ",
              tags$b("Integrated only"), " drops dark-cycle steps (large Δt, meaningless ",
              "midpoint). ", tags$b("Subtract noise floor"), " removes the ", sigloc(), "²/Δt ",
              "apparent-diffusion offset so immobile regions don't read as diffusing — ",
              "needs ", sigloc(), " set (run the Data Quality & Calibration panel, or enter it in Acquisition Parameters)."),
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
                    tags$td("MFX-TSA — Vogler et al. 2025 (Commun. Biol.); reference implementation Toolbox_for_MINFLUX_SPT (the_road_to_minflux_msd). Integrated/dark cycles, cycle time, adaptive gap removal.")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("MSD / OLSF")),
                    tags$td("Michalet 2010; Michalet & Berglund 2012 (optimal least-squares fitting); cycle-based binning & workflow per Vogler et al. 2025 (Toolbox_for_MINFLUX_SPT)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("CVE (", sigloc(), " / D)")),
                    tags$td("Vestergaard, Blainey & Flyvbjerg 2014 — restricted to consecutive integrated (single-cycle) pairs for MINFLUX")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("JDD")),
                    tags$td("Schütz, Schindler & Schmidt 1997; benchmarked in Weimann et al. 2013")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("HMM states")),
                    tags$td("Das, Cairo & Coombs 2009; vbSPT (Persson et al. 2013)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Confinement")),
                    tags$td("Rg/Rmax shape ratio + τ-based MSD confinement index (cf. Simson, Sheets & Jacobson 1995)")),
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
            "Motion blur R≈0 is a modelling assumption for MINFLUX (short cycle + active re-centring), not a published value.")
        )
      )
    )  # end tabPanel("Analysis")

  )  # end tabsetPanel
)
