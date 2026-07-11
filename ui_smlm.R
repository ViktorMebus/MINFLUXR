## ─────────────────────────────────────────────────────────────────────────────
## ui_smlm.R  —  SMLM analysis tab UI (assigned to smlm_ui)
## ─────────────────────────────────────────────────────────────────────────────
smlm_ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$script(HTML("
(function() {
  var _rotateTimer   = null;
  var _recordTimer   = null;
  var _mediaRecorder = null;
  var _chunks        = [];
  var _rotateAngle   = 0;

  /* ── helpers ── */
  function getPlotDiv() {
    return document.getElementById('scatter3D');
  }
  function currentCamera() {
    var d = getPlotDiv();
    if (!d || !d.layout || !d.layout.scene) return null;
    return d.layout.scene.camera || null;
  }

  /* ── rotation ── */
  window.startRotate3D = function(speedDeg, elevationDeg, distance) {
    if (_rotateTimer !== null) clearInterval(_rotateTimer);
    var elev = (elevationDeg || 30) * Math.PI / 180;
    var dist = distance || 2.5;
    _rotateAngle = 0;
    _rotateTimer = setInterval(function() {
      _rotateAngle += (speedDeg || 1);
      if (_rotateAngle >= 360) _rotateAngle -= 360;
      var theta = _rotateAngle * Math.PI / 180;
      var eye = {
        x: dist * Math.cos(elev) * Math.sin(theta),
        y: dist * Math.cos(elev) * Math.cos(theta),
        z: dist * Math.sin(elev)
      };
      var d = getPlotDiv();
      if (d) Plotly.relayout(d, {'scene.camera.eye': eye});
    }, 50);
  };

  window.stopRotate3D = function() {
    if (_rotateTimer !== null) { clearInterval(_rotateTimer); _rotateTimer = null; }
  };

  /* ── recording ── */
  var _recordActive = false;
  window.startRecord3D = function(width, height, fps, bitrateMbps) {
    if (_recordTimer !== null) return;   // already recording
    var w        = width  || 800;
    var h        = height || 600;
    var fps_use  = fps || 20;
    var bps      = (bitrateMbps || 10) * 1000000;
    var interval = Math.round(1000 / fps_use);

    var canvas = document.createElement('canvas');
    canvas.width  = w;
    canvas.height = h;
    var ctx = canvas.getContext('2d');

    var stream = canvas.captureStream(fps_use);
    _chunks = [];
    var mimeTypes = ['video/webm; codecs=vp9', 'video/webm; codecs=vp8', 'video/webm'];
    var mimeType = '';
    for (var i = 0; i < mimeTypes.length; i++) {
      if (MediaRecorder.isTypeSupported(mimeTypes[i])) { mimeType = mimeTypes[i]; break; }
    }
    var recOpts = { videoBitsPerSecond: bps };
    if (mimeType) recOpts.mimeType = mimeType;
    try {
      _mediaRecorder = new MediaRecorder(stream, recOpts);
    } catch(e) {
      _mediaRecorder = new MediaRecorder(stream);
    }
    _mediaRecorder.ondataavailable = function(e) {
      if (e.data && e.data.size > 0) _chunks.push(e.data);
    };
    _mediaRecorder.start(200);   // collect chunks every 200 ms

    var d = getPlotDiv();
    _recordActive = true;
    function captureFrame() {
      if (!_recordActive || !d) return;
      Plotly.toImage(d, { format: 'png', width: w, height: h }).then(function(dataUrl) {
        if (!_recordActive) return;
        var img = new Image();
        img.onload = function() {
          ctx.drawImage(img, 0, 0, w, h);
          _recordTimer = setTimeout(captureFrame, interval);
        };
        img.src = dataUrl;
      }).catch(function() {
        if (_recordActive) _recordTimer = setTimeout(captureFrame, interval);
      });
    }
    _recordTimer = setTimeout(captureFrame, 0);
  };

  window.stopRecord3D = function() {
    _recordActive = false;
    if (_recordTimer !== null) { clearTimeout(_recordTimer); _recordTimer = null; }
    if (_mediaRecorder && _mediaRecorder.state !== 'inactive') {
      _mediaRecorder.onstop = function() {
        var blob = new Blob(_chunks, { type: 'video/webm' });
        var url  = URL.createObjectURL(blob);
        var a    = document.createElement('a');
        a.href     = url;
        a.download = 'minfluxr_3d_rotation.webm';
        document.body.appendChild(a);
        a.click();
        setTimeout(function() { URL.revokeObjectURL(url); a.remove(); }, 1000);
        _chunks = [];
      };
      _mediaRecorder.stop();
    }
  };
})();
  "))),
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
  titlePanel("MINFLUXR - SMLM"),

  tabsetPanel(
    id = "smlm_tabs",

    # ── Display ───────────────────────────────────────────────────────────────
    tabPanel("Display",
      # Sticky Plot button bar — always visible regardless of sidebar scroll
      tags$div(
        style = paste(
          "position: sticky; top: 0; z-index: 200;",
          "background: #f8f9fa; border-bottom: 1px solid #dee2e6;",
          "padding: 6px 15px; display: flex; align-items: center; gap: 12px;"
        ),
        actionButton("plot_data", "Plot Data",
                     icon  = icon("play"),
                     class = "btn-primary btn-sm"),
        uiOutput("plot_bar_file_label")
      ),
      sidebarLayout(
        sidebarPanel(
          style = "height: 90vh; overflow-y: auto;",
      ### 1) Upload ===
      tags$details(class = "sim-group", open = NA,
      tags$summary(tags$b("Data Upload")),
      helpText("Single- or multi-file loader"),
      fileInput("file", "Upload MINFLUX File (.npy, .json or .csv)",
                accept = c(".csv", ".json", ".npy"), multiple = TRUE),
      checkboxInput("pool_dcr_json",
                    "Pool DCR (JSON & .npy)",
                    value = FALSE),
      helpText("Active file:"),
      selectInput("active_file", NULL,
                  choices = character(0), width = "100%"),
      actionButton("flush_file", "Remove Selected File", icon = icon("trash"),
                   class = "btn-warning btn-sm", style = "margin-bottom:8px;"),
      checkboxInput("lock_filters", "Lock filters", value = FALSE),
      radioButtons(
        "mode3d", "Display mode",
        choices = c("2D", "3D"),
        inline = TRUE,
        selected = "2D"
      ),
      conditionalPanel(
        condition = "input.mode3d == '2D'",
        uiOutput("view_plane_2d_ui")
      ),
      ),  # /Data Upload
      tags$hr(),

      ### 1b) Multiplexed MINFLUX ===
      tags$details(class = "sim-group",
      tags$summary(tags$b("Multiplexed MINFLUX")),
      helpText("Build an N-channel dataset for exchange/sequential MINFLUX: overlay several loaded files, or split one continuous acquisition by time (TIM). Each channel gets its own color and feeds the plot and the Analysis tab."),
      checkboxInput("mplex_enable", "Enable multiplexing", value = FALSE),
      conditionalPanel(
        condition = "input.mplex_enable == true",
        radioButtons("mplex_mode", "Channel source",
                     choices = c("Overlay loaded files"          = "files",
                                 "Split active file by time (TIM)" = "tim"),
                     selected = "files"),
        conditionalPanel(
          condition = "input.mplex_mode == 'files'",
          selectInput("mplex_files", "Files to overlay (one channel each)",
                      choices = character(0), multiple = TRUE, width = "100%"),
          uiOutput("mplex_files_colors_ui")
        ),
        conditionalPanel(
          condition = "input.mplex_mode == 'tim'",
          uiOutput("mplex_tim_info"),
          radioButtons("mplex_tim_method", "Time assignment",
                       choices = c("Contiguous breakpoints"          = "breaks",
                                   "Explicit windows (drop washes)"   = "windows"),
                       selected = "breaks"),
          conditionalPanel(
            condition = "input.mplex_tim_method == 'breaks'",
            textInput("mplex_tim_breaks", "Breakpoints (seconds, comma-separated)",
                      value = "", placeholder = "e.g. 600, 1200"),
            helpText("N breakpoints create N+1 contiguous time-window channels. Up to 6 splits (7 channels).")
          ),
          conditionalPanel(
            condition = "input.mplex_tim_method == 'windows'",
            textInput("mplex_tim_windows", "Included windows (s): start-end, start-end, …",
                      value = "", placeholder = "e.g. 0-600, 700-1300, 1400-2000"),
            helpText("Each window = one channel. Localizations OUTSIDE all windows (wash/exchange steps) are removed. Up to 7 windows.")
          ),
          plotlyOutput("mplex_tim_hist", height = "200px"),
          sliderInput("mplex_tim_bins", "Histogram bins",
                      min = 10, max = 500, value = 60, step = 5),
          helpText("Trace IDs over time; bar color = channel, grey = removed (wash), dashed lines = window/split edges."),
          uiOutput("mplex_tim_colors_ui")
        ),
        fluidRow(
          column(7, actionButton("mplex_build", "Build multiplex",
                                 icon = icon("layer-group"),
                                 class = "btn-primary btn-sm", style = "width:100%")),
          column(5, actionButton("mplex_clear", "Clear",
                                 class = "btn-default btn-sm", style = "width:100%"))
        ),
        verbatimTextOutput("mplex_summary")
      ),
      ),  # /Multiplexed MINFLUX
      tags$hr(),

      ### 2) Metadata Filtering ===
      tags$details(class = "sim-group", open = NA,
      tags$summary(tags$b("Metadata Filtering")),
      tags$div(id = "filters_form",
               checkboxInput("show_metadata_hist", "Show Filtering Controls", TRUE),
               conditionalPanel(
                 condition = "input.show_metadata_hist == true",
                 selectInput("hist_variable","Select Variable for Histogram",
                             choices = c("efo","tid","cfr","tim","eco","dcr","dwell","fluo"),
                             selected = "efo"),
                 plotlyOutput("metadata_histogram", height = "300px"),
                 sliderInput("hist_bins","Histogram Bins", min = 10, max = 1000, value = 200),
                 
                 uiOutput("hist_range_ui"),
                 uiOutput("efo_range_ui"),
                 uiOutput("tid_length_range_ui"),
                 uiOutput("cfr_range_ui"),
                 uiOutput("tim_range_ui"),
                 uiOutput("dcr_range_ui"),
                 uiOutput("eco_range_ui"),
                 uiOutput("dwell_range_ui"),
                 uiOutput("fluo_range_ui"),
                 
                 checkboxInput("hist_show_raw", "Show raw data (grey)", value = TRUE),
                 
                 downloadButton("download_ui_settings", "Download Filter Settings"),
                 
                 fileInput("load_settings", "Load Filter Settings (.json)", accept = ".json"),
                 
                 
               ),
               # Trace spread z control appears only in 3D
               checkboxInput("view_SD", "Inspect trace spread (\u03C3\u0303)", value = FALSE),
               conditionalPanel(
                 condition = "input.view_SD == true",
                 numericInput("sd_min_locs",
                              "Min localizations per trace for spread (σ̃)",
                              value = MIN_LOCS_FOR_SPREAD_DEFAULT, min = 2, step = 1),
                 helpText("Traces shorter than this are excluded from the spread histogram, ",
                          "the median-spread labels, the allowable-spread ranges, and the σ filters. ",
                          "Literature typically uses ≥ 5 (a SD needs ≥ 2 points to exist)."),
                 uiOutput("sd_summary_ui"),
                 uiOutput("sd_hist_ui"),
                 sliderInput("sd_bins", "Histogram bins", min = 5, max = 500, value = 100, step = 1),
                 checkboxInput("sd_overlay_all", "Show raw trace spreads (grey bars)", value = FALSE),
                 checkboxInput("sd_show_mean_sd", "Show median trace spread (post filtering)", value = FALSE),
                 checkboxInput("sd_show_overall_mean_sd", "Show median trace spread (raw data)", value = FALSE),
                 tags$span(
                   " (What is \u03C3\u0303 \u00b1 nMAD?)",
                   title = paste(
                     "\u03c3\u0303 = per-trace localization spread (SD of x/y/z within each trace).",
                     "Line label: median(\u03c3\u0303) \u00b1 nMAD across all traces.",
                     "nMAD = 1.4826 \u00d7 median(|\u03c3\u0303_i \u2212 median(\u03c3\u0303)|).",
                     "Factor 1.4826 \u2248 1/\u03a6\u207b\u00b9(3/4): makes nMAD a consistent estimator of \u03c3 for normal data (same as pyMINFLUX).",
                     "Traces shorter than the 'Min localizations per trace' setting are excluded."
                   )
                 ),
                 checkboxInput("sd_use_q95", "Set spread ranges to the 95th percentile", value = FALSE),
                 uiOutput("sd_ranges_ui"),
                 checkboxInput("remove_rogues", "Remove rogue traces by rule", value = FALSE),
                 conditionalPanel(
                   condition = "input.remove_rogues == true",
                   radioButtons(
                     "rogue_rule", "Rogue rule",
                     choices = c("Spread > k \u00d7 dataset median"         = "factor_median",
                                 "Spread > median + k \u00d7 MAD (robust)"  = "mad_rule",
                                 "Any axis spread > k \u00d7 median(axis spread)" = "axis_sd"),
                     selected = "factor_median"
                   ),
                   numericInput("rogue_factor", "Multiplier k", value = 2.5, min = 1, step = 0.1),
                   helpText("k applies to whichever rule is selected above. Dropped traces are removed before grouping, clustering, fits, and surfaces.")
                 )
               ),
               actionButton("reset_filters", "Reset Filters", icon = icon("refresh")),
               downloadButton("download_filtered_csv", "Export MINFLUX data (.csv)"),
               helpText("Exports data with current filters applied.")


      ),
      ),  # /Metadata Filtering
      tags$hr(),

      ### Visualization Options (shared + mode-specific bits) ===
      tags$details(class = "sim-group", open = NA,
      tags$summary(tags$b("Visualization Options")),
        numericInput("median_point_size", "Point size", value = 4, min = 1, max = 20),
        
        selectInput(
          "color_mode", "Color points by:",
          choices = c(
            "Trace/Cluster color (default)" = "default",
            "Multiplex channel"                      = "multiplex",
            "Unicolor (one for all)"                    = "unicolor",
            "Depth (z)"                              = "z",
            "Time (tim)"                             = "tim",
            "Trace-median spatial density"           = "locdens"
          ),
          selected = "default"
        ),
        
        conditionalPanel(
          condition = "input.color_mode == 'unicolor'",
          colourpicker::colourInput("unicolor_hex", "Unicolor", value = "#1f77b4")
        ),
        
        conditionalPanel(
          condition = "['z','tim','locdens'].includes(input.color_mode)",
          selectInput(
            "numeric_palette", "Numeric colorscale",
            choices = NUMERIC_PALETTES,
            selected = "Spectral"
          ),
          conditionalPanel(
            condition = "input.color_mode == 'locdens'",
            sliderInput("dens_radius_nm", "Spatial density radius (nm)",
                        min = 10, max = 5000, value = 500, step = 10),
            helpText("Counts how many other trace-median positions fall within this radius. Colors traces by that local density (not localization count).")
          ),
          checkboxInput("numeric_palette_reverse", "Reverse colorscale", value = FALSE)),
        
        checkboxInput("show_medians", "Show trace median positions (hide individual localizations)", value = FALSE),
        checkboxInput("dark_mode", "Dark mode", value = FALSE),
        checkboxInput("hide_grid", "Hide gridlines & axis titles", value = FALSE),
        conditionalPanel(
          condition = "input.mode3d == '2D'",
                checkboxInput("show_scale_bar", "Add scale bar", value = FALSE),
        conditionalPanel(
          condition = "input.show_scale_bar == true && input.mode3d == '2D'",
          numericInput("scale_bar_length", "Scale bar length (nm)", value = 100, min = 1, step = 10),
          numericInput("scale_bar_thickness", "Scale bar thickness (nm)", value = 10, min = 1, step = 1),
          numericInput("scale_font_size", "Font size", value = 18, min = 1, step = 1)
          )
        ),

        checkboxInput("make_positive", "Transform coordinates to positive space", value = FALSE),
        
        # 3D-only: alpha/mesh wrapper
        conditionalPanel(
          condition = "input.mode3d == '3D'",
          checkboxInput("show_convex_hull", "Show Alpha Shape (3D)", value = FALSE),
          helpText("Wraps traces (or clusters) with a 3D alpha surface."),
          conditionalPanel(
            condition = "input.show_convex_hull == true",
            sliderInput("alpha_value_hull", "Alpha radius (nm) \u2014 smaller = tighter mesh", min = 0.1, max = 1000, value = 200, step = 0.1),
            sliderInput("hull_opacity", "Alpha shape opacity", min = 0.1, max = 1.0, value = 0.2, step = 0.05)
          )
        ),
        # --- 2D-only: Crop to current zoom view ---
        conditionalPanel(
          condition = "input.mode3d == '2D'",
          tags$hr(style = "margin:6px 0;"),
          tags$b("Crop to view", style = "font-size:0.9em;"),
          uiOutput("smlm_spatial_info"),
          fluidRow(
            column(7, actionButton("smlm_crop_to_view", "Crop to view",
                                   icon = icon("crop-alt"), class = "btn-info btn-sm",
                                   style = "width:100%")),
            column(5, actionButton("smlm_clear_crop", "Clear",
                                   class = "btn-default btn-sm", style = "width:100%"))
          ),
          helpText("Zoom into a region, press Crop to view. The X/Y crop is retained when switching to 3D (visible in the Crop to region controls below).")
        ),

        # --- 3D-only: Cropping controls ---
        conditionalPanel(
          condition = "input.mode3d == '3D'",
          uiOutput("smlm_spatial_info_3d"),
          actionButton("smlm_clear_crop_3d", "Clear 2D crop",
                       class = "btn-default btn-sm", style = "width:100%; margin-bottom:6px;"),
          checkboxInput("enable_trimmer", "Crop to region (3D)", value = FALSE),
          conditionalPanel(
            condition = "input.enable_trimmer == true",
            radioButtons(
              "trim_mode", "Trimmer shape",
              choices = c("Axis-aligned box" = "box",
                          "Sphere"           = "sphere",
                          "Ellipsoid (PCA)"  = "ellip"),
              selected = "box"
            ),
            conditionalPanel(condition = "input.trim_mode == 'box'",   uiOutput("trim_box_ui")),
            conditionalPanel(condition = "input.trim_mode == 'sphere'", uiOutput("trim_sphere_ui")),
            conditionalPanel(condition = "input.trim_mode == 'ellip'",  uiOutput("trim_ellip_ui"))
          ),
          ### 3D Auto-rotate / Record ===
          checkboxInput("show_rotate_opt","Auto-Rotate & Record Options", value = FALSE),
          conditionalPanel(
            condition = "input.show_rotate_opt == true && input.mode3d == '3D'",
            checkboxInput("auto_rotate_3d", "Auto-rotate scene", value = FALSE),
            conditionalPanel(
              condition = "input.auto_rotate_3d == true",
              sliderInput("rotate_speed_3d", "Speed (°/frame)", min = 0.5, max = 10,  value = 2,   step = 0.5),
              sliderInput("rotate_elev_3d",  "Elevation (°)",   min = -89, max = 89,  value = 25,  step = 1),
              sliderInput("rotate_dist_3d",  "Camera distance", min = 0.3, max = 6,   value = 2.5, step = 0.1)
            ),
            tags$hr(),
            h5("Record WebM movie"),
            fluidRow(
              column(6, numericInput("record_width_3d",  "Width (px)",  value = 800, min = 200, max = 2560, step = 10)),
              column(6, numericInput("record_height_3d", "Height (px)", value = 600, min = 200, max = 1440, step = 10))
            ),
            numericInput("record_fps_3d", "FPS", value = 20, min = 5, max = 60, step = 1),
            sliderInput("record_bitrate_3d", "Video bitrate (Mbps)", min = 2, max = 25, value = 10, step = 1),
            fluidRow(
              column(6, actionButton("record_start_3d", "Start recording", class = "btn-danger btn-sm")),
              column(6, actionButton("record_stop_3d",  "Stop & download", class = "btn-success btn-sm"))
            ),
            helpText("Start auto-rotate before recording. Click 'Stop & download' to save the .webm file.")
            
          ))

      ),  # /Visualization Options
      tags$hr(),

      conditionalPanel(
        condition = "input.mode3d == '2D'",
        tags$details(class = "sim-group",
        tags$summary(tags$b("Background image overlay")),
        checkboxInput("show_confocal", "Show overlay settings", FALSE),
        conditionalPanel(
          condition = "input.show_confocal == true",
          tags$div(id = "overlay_form", style = "border:0; padding:0; margin:0; background:none; box-shadow:none;",
                   checkboxInput("show_bg", "Show background image", FALSE),
                   fileInput("bg_image", "Upload image", accept = c(".png",".jpg",".jpeg",".tif",".tiff")),
                   numericInput("px_size", "Pixel size (units/pixel)", value = 10, min = 0.0001),
                   numericInput("img_x0", "Image origin X (data units)", value = 0),
                   tags$div(
                     style = "display:flex; gap:3px; margin:-4px 0 8px 0; flex-wrap:wrap;",
                     tags$small("X:", style = "line-height:26px; color:#555; margin-right:2px;"),
                     actionButton("nx_m1000", "−1k",  class = "btn-xs btn-default"),
                     actionButton("nx_m100",  "−100", class = "btn-xs btn-default"),
                     actionButton("nx_m10",   "−10",  class = "btn-xs btn-default"),
                     actionButton("nx_p10",   "+10",       class = "btn-xs btn-default"),
                     actionButton("nx_p100",  "+100",      class = "btn-xs btn-default"),
                     actionButton("nx_p1000", "+1k",       class = "btn-xs btn-default")
                   ),
                   numericInput("img_y0", "Image origin Y (data units)", value = 0),
                   tags$div(
                     style = "display:flex; gap:3px; margin:-4px 0 8px 0; flex-wrap:wrap;",
                     tags$small("Y:", style = "line-height:26px; color:#555; margin-right:2px;"),
                     actionButton("ny_m1000", "−1k",  class = "btn-xs btn-default"),
                     actionButton("ny_m100",  "−100", class = "btn-xs btn-default"),
                     actionButton("ny_m10",   "−10",  class = "btn-xs btn-default"),
                     actionButton("ny_p10",   "+10",       class = "btn-xs btn-default"),
                     actionButton("ny_p100",  "+100",      class = "btn-xs btn-default"),
                     actionButton("ny_p1000", "+1k",       class = "btn-xs btn-default")
                   ),
                   tags$div(
                     style = "display:flex; gap:5px; margin-bottom:8px;",
                     actionButton("auto_foci",     "Auto: foci",
                                  icon = icon("bullseye"), class = "btn-sm btn-info"),
                     actionButton("auto_boundary", "Auto: boundary",
                                  icon = icon("draw-polygon"), class = "btn-sm btn-info")
                   ),
                   checkboxInput("show_diag", "Show detection diagnostics", FALSE),
                   conditionalPanel(
                     "input.show_diag",
                     plotOutput("foci_debug_plot", height = "220px")
                   ),
                   tags$details(
                     tags$summary(style = "cursor:pointer; font-size:12px; color:#555; margin:3px 0 6px 0;",
                                  "Foci detection settings"),
                     tags$div(style = "padding:4px 0 4px 8px; border-left:2px solid #ddd;",
                       sliderInput("foci_n_peaks",     "Max foci to find",             min=1, max=10, value=5,    step=1),
                       sliderInput("foci_bright_pct",  "Confocal: sensitivity (70=permissive → 99=strict)", min=70, max=99, value=90, step=1),
                       sliderInput("foci_border_pct",  "Confocal: edge mask (%)",      min=2,  max=25, value=10,  step=1),
                       sliderInput("foci_eps_pct",     "Confocal: min separation between foci (% of image)", min=1, max=15, value=5, step=0.5),
                       sliderInput("foci_density_thr",    "MINFLUX: density threshold (% of peak)", min=5,  max=60,  value=20,  step=5),
                       sliderInput("foci_density_cell",   "MINFLUX: cell size (nm)",                min=20, max=500, value=100, step=10),
                       sliderInput("foci_density_smooth", "MINFLUX: smoothing radius (nm; 0=off)",  min=0,  max=500, value=100, step=25)
                     )
                   ),
                   tags$p("Pixel size must be set first. Results are a starting point — nudge to refine.",
                          style = "font-size:11px; color:#888; margin:0 0 6px 0;"),
                   checkboxInput("img_invert_y", "Image origin at top-left (invert Y)", TRUE),
                   checkboxInput("img_flip_x", "Mirror horizontally (left↔right)", FALSE),
                   sliderInput("img_opacity", "Image opacity", min = 0, max = 1, value = 0.6, step = 0.05),
                   selectInput("img_lut", "Colour (LUT)",
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
                     condition = "input.img_lut == 'custom'",
                     colourpicker::colourInput("img_lut_custom", "Custom colour",
                                               value = "#00FF00")
                   ),
                   tags$p("Recolours the image to match a fluorophore: converts to grayscale intensity, then tints black→colour.",
                          style = "font-size:11px; color:#888; margin:-4px 0 8px 0;"),
                   sliderInput("img_brightness", "Image brightness (%)",
                               min = 10, max = 300, value = 100, step = 5),
                   sliderInput("img_contrast", "Image contrast",
                               min = 0, max = 90, value = 0, step = 5),
                   selectInput("img_rotation", "Rotate image", choices = c("0","90","180","270"), selected = "0"),
                   downloadButton("download_overlay_settings", "Download Overlay Settings"),
                   fileInput("load_overlay_settings", "Load Overlay Settings (.json)", accept = ".json"),
                   actionButton("reset_overlay", "Reset Overlay", icon = icon("undo"))
          )
        ),
        ),  # /Background image overlay
        tags$hr(),
        tags$details(class = "sim-group",
        tags$summary(tags$b("Heatmap overlay")),
        checkboxInput("show_heatmap", "Show heatmap (KDE)", FALSE),
        conditionalPanel(
          condition = "input.show_heatmap == true",
          numericInput("heatmap_n", "Grid resolution (n × n)", value = 200, min = 25, max = 1000),
          sliderInput("heatmap_opacity", "Opacity", min = 0, max = 1, value = 0.35, step = 0.05),
          selectInput("heatmap_palette", "Palette",
                      choices = NUMERIC_PALETTES,
                      selected = "YlOrRd"),
          sliderInput("heatmap_bw_factor", "Smoothing (× Scott rule)", min = 0.25, max = 4, value = 1, step = 0.05),
          numericInput("heatmap_bw_abs", "Bandwidth (nm, overrides \u00d7 Scott rule factor)", value = NA, min = 0),
          sliderInput("heatmap_floor", "Suppress bottom percentile", min = 0, max = 0.50, value = 0, step = 0.01),
          checkboxInput("heatmap_use_raw", "Use raw points (instead of trace medians)", FALSE),
          sliderInput("heatmap_clip", "Clip top percentile", min = 0.80, max = 1.00, value = 0.99, step = 0.005),
          sliderInput("heatmap_gamma", "Gamma (power)", min = 0.2, max = 4, value = 1, step = 0.1),
          checkboxInput("heatmap_log", "Log scale", FALSE),
          checkboxInput("heatmap_showlegend", "Show colorbar", FALSE)
        ),
        ),  # /Heatmap overlay
        tags$hr()
      ),


      ### 2-Color settings (shared) ===
      tags$details(class = "sim-group",
      tags$summary(tags$b("2-color MINFLUX")),
      checkboxInput("enable_two_color", "Enable 2-color (per-trace median)", value = FALSE),
      helpText("Note: Numeric color modes (z/time/density) override 2-color. Use 'Default' color mode to see 2-color groups."),
      conditionalPanel(
        condition = "input.enable_two_color == true",
        # DCR ranges (default mode when not using fluo)
        checkboxInput("color_by_fluo", "Override with 'fluo' groups (1/2)", value = FALSE),
        conditionalPanel(
          condition = "input.color_by_fluo == false",
          shinyWidgets::numericRangeInput("dcr_group1", "DCR range: Group 1",
                                          value = c(0, 0.369), min = 0, max = 2),
          shinyWidgets::numericRangeInput("dcr_group2", "DCR range: Group 2",
                                          value = c(0.3691, 1), min = 0, max = 2),
          actionButton("auto_detect_dcr", "Auto-detect split (GMM)",
                       icon = icon("wand-magic-sparkles"), class = "btn-sm btn-info"),
          textOutput("dcr_gmm_info"),
          uiOutput("use_pooled_dcr_ui")
        ),
        
        colourpicker::colourInput("color_group1", "Color: Group 1", value = "#1f77b4"),
        colourpicker::colourInput("color_group2", "Color: Group 2", value = "#d62728"),
        
        # Scope pickers for NN/cluster
        # Shared scope pickers for NN/cluster (works for 2D and 3D)
        checkboxInput("draw_nn_two_color", "Draw NN lines (2-color mode)", value = FALSE),
        conditionalPanel(
          condition = "input.draw_nn_two_color == true",
          numericInput("nn_k_two_color", "2-color NN: Neighbors (k)", value = 4, min = 1, step = 1),
          checkboxInput("show_nn_labels_two_color", "Show distance labels (2-color NN)", value =  TRUE),
          conditionalPanel(
            condition = "input.show_nn_labels_two_color == true",
            numericInput("tc_label_font_size", "2-color NN label font size", value = 12, min = 6, step = 1),
            numericInput("tc_label_offset", "2-color NN label offset (nm)", value = 5, min = 0, step = 1),
          ),
          sliderInput(
            "distance_range_two_color", "2-color NN Distance Range", min = 0, max = 2000, value = c(0, 100)),
          selectInput("two_color_scope", "2-color NN scope", choices = c("Between groups" = "between",
                                                                         "Within groups"  = "within"),
                      selected = "between"),
          
          downloadButton("download_nn_cross_dcr", "Download 2-color NN table"),
        ),
        checkboxInput("enable_two_color_cluster1", "2-color Primary clustering (raw points)", FALSE),
        conditionalPanel(
          condition = "input.enable_two_color_cluster1",
          numericInput("eps_tc1", "2-color DBSCAN eps, raw points (nm)", value = 14, min = 1),
          numericInput("minpts_tc1", "2-color DBSCAN minPts (raw)", value = 6, min = 1),
          conditionalPanel(
            condition = "input.mode3d == '3D'",
            checkboxInput("show_tc1_alpha", "Wrap 2-color primary clusters (tc1)", value = FALSE),
            conditionalPanel(
              condition = "input.show_tc1_alpha == true && input.mode3d == '3D'",
              numericInput("tc1_alpha", "Alpha (smaller = tighter)", value = 120, min = 1, step = 1),
              sliderInput("tc1_alpha_opacity", "Alpha shape opacity", min = 0, max = 1, value = 0.15, step = 0.01),
              numericInput("tc1_alpha_min_pts", "Min points per cluster to wrap", value = 5, min = 2, step = 1)
            ),
          ),
          checkboxInput("exclude_unclustered_tc1", "Exclude unclustered (2-color cluster1)", FALSE)
        ),
        checkboxInput("enable_two_color_cluster2", "2-color Median Clustering (TID medians)", FALSE),
        conditionalPanel(
          condition = "input.enable_two_color_cluster2",
          radioButtons("cluster2_method_tc2", "2-color Median Clustering grouping method",
                       choices = c("DBSCAN" = "dbscan", "HDBSCAN" = "hdbscan", "OPTICS (Xi)" = "optics"),
                       selected = "dbscan", inline = TRUE),
          conditionalPanel(
            condition = "input.cluster2_method_tc2 == 'dbscan'",
            numericInput("eps_tc2", "2-color eps, trace medians (nm)", value = 150, min = 1)
          ),
          conditionalPanel(
            condition = "input.cluster2_method_tc2 == 'optics'",
            numericInput("cluster2_xi_tc2", "2-color OPTICS Xi (0-1)", value = 0.05, min = 0.001, max = 0.5, step = 0.005)
          ),
          numericInput("minpts_tc2", "2-color MinPts (medians)", value = 8, min = 1),
          checkboxInput("cluster2_prefilter_tc2", "2-color density prefilter", value = FALSE),
          conditionalPanel(
            condition = "input.cluster2_prefilter_tc2 == true",
            numericInput("cluster2_prefilter_r_tc2", "2-color prefilter radius (nm)", value = 150, min = 1),
            numericInput("cluster2_prefilter_min_tc2", "2-color min neighbours", value = 3, min = 1)
          ),
          checkboxInput("cluster2_prune_tc2", "2-color prune sparse cluster members", value = FALSE),
          conditionalPanel(
            condition = "input.cluster2_prune_tc2 == true",
            numericInput("cluster2_prune_tol_tc2", "2-color prune tolerance (x median spacing)", value = 2, min = 1.1, step = 0.1)
          ),
          conditionalPanel(
            condition = "input.mode3d == '3D'",
            checkboxInput("show_tc2_alpha", "Wrap 2-color secondary clusters (tc2)", value = FALSE),
            conditionalPanel(
              condition = "input.show_tc2_alpha == true && input.mode3d == '3D'",
              numericInput("tc2_alpha", "tc2 alpha", value = 150, min = 1, step = 1),
              sliderInput("tc2_alpha_opacity", "tc2 opacity", min = 0, max = 1, value = 0.18, step = 0.01),
              numericInput("tc2_alpha_min_pts", "tc2 min points to wrap", value = 15, min = 4, step = 1)
            ),
            checkboxInput("show_tc2_isosurface", "Show isosurface for 2-color secondary clusters (tc2)", value = FALSE),
            conditionalPanel(
              condition = "input.show_tc2_isosurface",
              checkboxInput("tc2_iso_selected_only", "Only selected tc2 cluster", value = TRUE),
              textInput("tc2_iso_cluster_label", "Selected tc2 cluster (A-# or B-#)", value = "A-1"),
              sliderInput("tc2_iso_grid_res", "Grid resolution (higher = slower)", min = 12, max = 40, value = 22, step = 1),
              sliderInput("tc2_iso_pad", "Padding around cluster (nm)", min = 0, max = 500, value = 100, step = 10),
              sliderInput("tc2_iso_level_q", "Surface level (quantile)", min = 0.50, max = 0.98, value = 0.85, step = 0.01),
              sliderInput("tc2_iso_sigma", "Isosurface smoothness (sigma, nm)", min = 10, max = 500, value = 150, step = 10),
              sliderInput("tc2_iso_opacity", "Isosurface opacity", min = 0, max = 1, value = 0.25, step = 0.01),
              sliderInput("tc2_iso_eps", "Smoothness (epsilon)", min = 1e-6, max = 1e-2, value = 1e-4, step = 1e-6)
            )
          )
        ),
        uiOutput("two_color_settings_file_ui"),
        downloadButton("download_two_color_settings", "Save 2-color settings"),
        actionButton("reset_two_color_settings", "Reset 2-color settings", icon = icon("undo"))
      ),
      ),  # /2-color MINFLUX
      tags$hr(),

      ### 5) Nearest Neighbor (shared) ===
      tags$details(class = "sim-group",
      tags$summary(tags$b("Nearest Neighbor Analysis")),
      checkboxInput("show_nn_lines", "Show NN Lines", value = FALSE),
      helpText("NN computed from median positions (or cluster medians if clustering is enabled)."),
      conditionalPanel(
        condition = "input.show_nn_lines == true",
        numericInput("nn_k", "Nearest Neighbors (k)", min = 1, max = 10, value = 4),
        checkboxInput("show_nn_labels", "Show NN distance labels", value = TRUE),
        conditionalPanel(
          condition = "input.show_nn_labels == true",
          numericInput("nn_label_font_size", "NN label font size", value = 12, min = 6, step = 1),
          numericInput("nn_label_offset", "NN label offset (nm)", value = 5, min = 0, step = 1),
        ),
        
        sliderInput("distance_range", "NN Distance Range", min = 0, max = 500, value = c(0, 100)),
        conditionalPanel(
          condition = "output.mplex_active == true",
          radioButtons("nn_mplex_scope", "Multiplex NN scope",
                       choices = c("Within channel"   = "within",
                                   "Between channels" = "across"),
                       selected = "within", inline = TRUE),
          helpText("Within: nearest neighbour searched only inside each channel. Between: nearest neighbour searched only in OTHER channels (cross-channel pairs, never within the same channel).")
        ),
        conditionalPanel(
          condition = "input.enable_cluster2 == true",
          checkboxInput("nn_within_cluster2", "Within cluster2 only", value = FALSE)
        ),
        downloadButton("download_nn_table", "Download NN Distance Table")
      ),
      ),  # /Nearest Neighbor Analysis
      tags$hr(),

      ### DBSCAN Point Clustering (shared) ===
      tags$details(class = "sim-group",
      tags$summary(tags$b("DBSCAN Point Clustering")),
      checkboxInput("enable_clustering", "Enable Primary Clustering (Raw Points)", value = FALSE),
      helpText("Clusters individual localizations and overrides Trace ID grouping."),
      conditionalPanel(
        condition = "input.enable_clustering == true",
        numericInput("eps", "DBSCAN Epsilon (nm)", value = 10, min = 1),
        numericInput("minpts", "Minimum Points", value = 4, min = 1),
        checkboxInput("exclude_unclustered", "Exclude unclustered points", value = FALSE),
        helpText("TIDs are not split, all localizations of a Trace ID are assigned to the same cluster (majority vote). TIDs that overlap spatially are merged into one cluster.")
      ),
      ),  # /DBSCAN Point Clustering
      tags$hr(),

      ### DBSCAN Median Clustering (Cluster2 calls) ===
      tags$details(class = "sim-group",
      tags$summary(tags$b("Median Clustering")),
      checkboxInput("enable_cluster2", "Enable Median-based Clustering", value = FALSE),
      helpText("Clusters trace medians into larger structures."),
      conditionalPanel(
        condition = "input.enable_cluster2 == true",
        radioButtons("cluster2_method", "Median Clustering grouping method",
                     choices = c("DBSCAN" = "dbscan", "HDBSCAN" = "hdbscan", "OPTICS (Xi)" = "optics"),
                     selected = "dbscan", inline = TRUE),
        helpText("DBSCAN: fixed radius (eps). HDBSCAN: adaptive, MinPts only (= min cluster size). OPTICS (Xi): adaptive density, keeps small structures; tuned by Xi."),
        conditionalPanel(
          condition = "input.cluster2_method == 'dbscan'",
          numericInput("eps2", "Secondary DBSCAN Epsilon (nm)", value = 200, min = 3)
        ),
        conditionalPanel(
          condition = "input.cluster2_method == 'optics'",
          numericInput("cluster2_xi", "OPTICS Xi (steepness, 0-1; lower = more/smaller clusters)", value = 0.05, min = 0.001, max = 0.5, step = 0.005)
        ),
        numericInput("minpts2", "Secondary MinPts (all methods)", value = 15, min = 1),
        helpText("Tip: set MinPts near the smallest number of TIDs you'd call a real structure (often 3-6), not higher."),
        # ── local-density prefilter (drop scattered background TIDs) ──
        checkboxInput("cluster2_prefilter", "Density prefilter (drop sparse background TIDs first)", value = FALSE),
        conditionalPanel(
          condition = "input.cluster2_prefilter == true",
          numericInput("cluster2_prefilter_r", "Prefilter radius (nm)", value = 200, min = 1),
          numericInput("cluster2_prefilter_min", "Min neighbours within radius", value = 3, min = 1),
          helpText("Keeps only centroids with at least this many neighbours within the radius, then groups the survivors.")
        ),
        # ── sparse-member pruning (post-clustering trim of loose members) ──
        checkboxInput("cluster2_prune", "Prune sparse cluster members (trim loosely-attached traces)", value = FALSE),
        conditionalPanel(
          condition = "input.cluster2_prune == true",
          numericInput("cluster2_prune_tol", "Prune tolerance (x cluster-median spacing; lower = stricter)", value = 2, min = 1.1, step = 0.1),
          helpText("After clustering, demotes to noise any member whose local spacing exceeds this multiple of its cluster's typical spacing. ~2 trims clear outliers; lower is stricter.")
        ),
        # ── kNN-distance diagnostic to pick eps ──
        checkboxInput("show_cluster2_knn", "Show kNN-distance plot (tune eps / see density gap)", value = FALSE),
        conditionalPanel(
          condition = "input.show_cluster2_knn == true",
          plotOutput("cluster2_knn_plot", height = "220px"),
          helpText("Sorted distance to the (MinPts-1)-th neighbour of each centroid. The elbow marks the structure/background split = a good eps (red line = current eps).")
        ),
        conditionalPanel(
          condition = "input.mode3d =='3D' && input.enable_cluster2 == true",
          checkboxInput("show_cluster2_shape", "Wrap Cluster in Alpha Shape", value = FALSE),
          checkboxInput("show_cluster2_isosurface", "Show isosurface of Median Clusters", value = FALSE),
          conditionalPanel(
            condition = "input.show_cluster2_isosurface == true",
            checkboxInput("cluster2_iso_selected_only", "Only selected Median Cluster ID", value = FALSE),
            numericInput("cluster2_iso_id", "Selected Median Cluster ID", value = 1, min = 1, step = 1),
            checkboxInput("cluster2_iso_unicolor", "Single color for all clusters", value = FALSE),
            conditionalPanel(
              condition = "input.cluster2_iso_unicolor == true",
              colourpicker::colourInput(
                "cluster2_iso_color",
                "Isosurface color",
                value = "#00CED1"
              )
            ),
            sliderInput("cluster2_iso_grid_res", "Grid resolution (higher = slower)", min = 12, max = 40, value = 22, step = 1),
            sliderInput("cluster2_iso_pad", "Padding around cluster (nm)", min = 0, max = 500, value = 100, step = 10),
            sliderInput("cluster2_iso_level_q", "Surface level (quantile)", min = 0.50, max = 0.98, value = 0.85, step = 0.01),
            sliderInput("cluster2_iso_opacity", "Isosurface opacity", min = 0, max = 1, value = 0.25, step = 0.01),
            sliderInput("cluster2_iso_sigma", "Smoothness (sigma, nm)", min = 10, max = 500, value = 150, step = 10)
          )
        ),
        conditionalPanel(
          condition = "input.mode3d =='2D' && input.enable_cluster2 == true",
          checkboxInput("show_cluster2_shape_2D", "Wrap Median Cluster in Convex Hull", value = FALSE),
        ),
        conditionalPanel(
          condition = "input.mode3d =='3D' && input.enable_cluster2 == true",
          # 3D metrics/fit panel
          checkboxInput("show_cluster2_metrics", "Fitting & Metrics", value = FALSE)),
        
        # --- 3D fit controls ---
        conditionalPanel(
          condition = "input.mode3d == '3D' && input.show_cluster2_metrics == true",
          checkboxInput("compute_rmsd_metrics", "Compute all-model RMSD comparison (slower)", value = FALSE),
          checkboxInput("membrane_protein_options", "Enable Membrane Protein Options", value = FALSE),
          conditionalPanel(condition = "input.mode3d == '3D' && input.membrane_protein_options == true",
                           checkboxInput("show_cluster2_plane_angle", "Show Angle to Plane", value = FALSE),
                           checkboxInput("show_cluster2_plane_distance", "Show Distance to Plane", value = FALSE),
                           checkboxInput("show_bending_radius", "Show Bending Radius", value = FALSE),
                           helpText("ASSUMPTION: the cluster is a curved arc where PC1 is the chord half-length and PC2 is the sagitta (bow height). Formula: R = (chord\u00b2 + sagitta\u00b2)\u2009/\u2009(2\u00d7sagitta). Only valid for arc-like clusters (e.g. a partial ring viewed edge-on); not meaningful for complete rings, spheres, or helices.")),
          checkboxInput("show_cluster2_pca_axes", "Show PCA Axes", value = FALSE),
          checkboxInput("show_cluster2_labels", "Show Median Cluster summary text", value = FALSE),
          
          radioButtons(
            "fit_method", "Fit method",
            choices = c(
              "Heuristic (PCA / quantile extents)"    = "heuristic",
              "Weighted orthogonal distance" = "mle"
            ),
            selected = "heuristic"
          ),
          helpText("Heuristic: geometry derived analytically from PCA axes and quantile extents, uses no iteration, no centre optimization. MLE: minimizes the sum of squared Euclidean distances from each representative point to the nearest surface point. RMSD metrics are computed in the same units (nm) for both methods and are directly comparable."),
          conditionalPanel(
            condition = "input.fit_method == 'mle'",
            checkboxInput("fit_weighted", "Weight by per-trace precision (n/spread\u00b2, \u2248 1/SEM\u00b2)", FALSE),
            helpText("Each representative point is weighted by n\u2009/\u2009spread\u00b2, proportional to the inverse squared standard error of its median position. Traces with a single localization (spread undefined) receive the median spread of multi-localization traces as a conservative fallback"),
            radioButtons(
              "mle_center_mode", "MLE centre handling",
              choices = c(
                "Regularized: MAP prior toward robust centre"  = "reg",
                "Fixed to robust centre (analytical radius)"    = "fixed",
                "Free: bounded to data bounding box"           = "free"
              ),
              selected = "reg"
            ),
            helpText("Fixed: centre = convex-hull centroid; radius solved analytically (no iteration). Regularized: centre is free but penalised for moving away from the robust centre (Gaussian MAP prior, strength \u03bb). Free: centre fully optimised within an expanded data bounding box."),
            conditionalPanel(
              condition = "input.mle_center_mode == 'reg'",
              sliderInput("mle_center_lambda",
                          "Regularization strength \u03bb",
                          min = 0, max = 5, value = 0.1, step = 0.05),
              helpText("\u03bb = 0.1: a 5 nm centre shift adds \u22482.5 nm\u00b2 to the mean squared residual (\u224810% of a typical 5 nm residual cluster). \u03bb = 1: a 1 nm shift adds as much as 1 nm\u00b2. Scale-independent of cluster size.")
            )
          ),
          selectInput(
            "fit_type", "Visualize Cluster Fit",
            choices = c("None" = "none", "Spherical" = "sphere", "Ellipsoid" = "ellipsoid",
                        "Helix" = "helix", "Ring (3D circle)" = "ring", "Arc (partial ring)" = "arc"),
            selected = "none"
          ),
          conditionalPanel(
            condition = "input.fit_type == 'helix'",
            sliderInput("helix_turns", "Helix turns (initialisation hint)", min = 1, max = 10, value = 3),
            helpText("Used only to seed the pitch estimate. The fitted turn count is determined by the data extent and the optimised pitch \u2014 it will not be exactly this value.")
          ),
          conditionalPanel(
            condition = "input.fit_type == 'arc'",
            sliderInput("arc_span_deg", "Arc span hint (\u00b0, 0 = auto)", min = 0, max = 360, value = 270, step = 10),
            helpText("Sets the initial arc angular extent for the fit. 0 = auto-detected from the data\u2019s angular spread. The optimiser refines this value.")
          ),
          conditionalPanel(
            condition = "input.fit_type == 'helix' || input.fit_type == 'ring' || input.fit_type == 'arc'",
            checkboxInput("show_fit_tube", "Render as tube", value = FALSE),
            conditionalPanel(
              condition = "input.show_fit_tube == true",
              numericInput("fit_tube_radius", "Tube radius (nm)", value = 5, min = 1, step = 1)
            )
          ),
          checkboxInput("show_fit_distances", "Show Distances to Fit", value = FALSE),
          conditionalPanel(
            condition = "input.show_fit_distances == true",
            numericInput("max_dist_lines_nm", "Max line length (nm, 0 = all)", value = 0, min = 0, step = 5)
          ),
          numericInput("fit_opacity", "Fit Surface Opacity", min = 0.0, max = 1.0, value = 0.4, step = 0.05),
          downloadButton("download_cluster2_table", "Download Median Cluster Metrics")
        ),
        
        # NN within Median Clusters (shared)
        checkboxInput("show_cluster2_nn_lines", "Show NN Lines (Median Clusters)", value = FALSE),
        conditionalPanel(
          condition = "input.show_cluster2_nn_lines == true && input.enable_cluster2 == true",
          numericInput("nn_k_cluster2", "NN Count (Median Clusters)", value = 3, min = 1),
          checkboxInput("show_c2_nn_labels", "Show distance labels (Median Cluster NN)", value = TRUE),
          conditionalPanel(
            condition = "input.show_c2_nn_labels == true && input.show_cluster2_nn_lines == true",
            numericInput("c2_label_font_size", "Median Cluster NN label font size", value = 12, min = 6, step = 1),
            numericInput("c2_label_offset", "Median Cluster NN label offset (nm)", value = 5, min = 0, step = 1),
          ),
          sliderInput("distance_range_cluster2", "NN Distance Range (Median Clusters)", min = 0, max = 20000, value = c(0, 1250)),
          downloadButton("download_nn_table2", "Download NN Table (Median Clusters)")
        ),
        uiOutput("cluster2_settings_file_ui"),
        downloadButton("download_cluster2_settings", "Save Median Cluster settings"),
        actionButton("reset_cluster2_settings", "Reset Median Cluster settings", icon = icon("undo"))
      ),
      ),  # /Median Clustering (Cluster2)
      tags$hr(),

      ### Global reference plane / surface (3D only) ===
      conditionalPanel(
        condition = "input.mode3d == '3D'",
        tags$details(class = "sim-group",
        tags$summary(tags$b("Global Reference Plane and Surface")),
        checkboxInput("view_surface", "Compute Surface Mesh", value = FALSE),
        conditionalPanel(
          condition = "input.view_surface == true",
          textInput(
            "surface_exclude_tids",
            "Exclude TIDs (comma/space separated)",
            value = ""
          ),
          checkboxInput("show_xy_plane", "Overlay XY Reference Plane (PCA based)", value = FALSE),
          checkboxInput("plane_from_unclustered", "Use only unclustered points", value = FALSE),
          conditionalPanel(
            condition = "input.show_xy_plane == true",
            numericInput("plane_offset", "Z offset for XY plane", value = 0),
            downloadButton("download_plane_distance_table", "Download Plane Distances")
          ),
          
          selectInput("surface_method", "Surface method",
                      choices = c("None" = "blnk","LOESS heightmap" = "loess",
                                  "Polynomial Surface" = "poly",
                                  "Thin plate spline (TPS) (heavy)" = "tps"),
                      selected = "blnk"),
          conditionalPanel(
            condition = "input.surface_method != 'blnk'",
            sliderInput("surface_opacity", "Surface opacity", min = 0.05, max = 1, value = 0.35, step = 0.05),
            selectInput("surface_palette", "Surface colorscale",
                        choices = NUMERIC_PALETTES, selected = "Viridis"),
            checkboxInput("surface_palette_reverse", "Reverse colorscale", value = FALSE),
            selectInput("surface_color_group", "Use points from",
                        choices = c("All groups" = "all", "Group 1 only" = "1", "Group 2 only" = "2"),
                        selected = "all"),
            helpText("Requires 2-color mode or a DCR/fluo grouping mode to be active."),
            
            # LOESS-only controls
            conditionalPanel(
              condition = "input.surface_method == 'loess'",
              sliderInput("surface_span", "Smoothness (loess span)",
                          min = 0.05, max = 1.0, value = 0.25, step = 0.05)
            ),
            
            # Shared (LOESS + poly + tps): grid resolution
            conditionalPanel(
              condition = "input.surface_method == 'loess' || input.surface_method == 'poly' || input.surface_method == 'tps'",
              numericInput("surface_grid_n", "Grid resolution", value = 120, min = 40, max = 400)
            ),
            
            # Polynomial-only
            conditionalPanel(
              condition = "input.surface_method == 'poly'",
              sliderInput("poly_degree", "Polynomial degree", min = 1, max = 5, value = 2)
            ),
            conditionalPanel(
              condition = "input.surface_method == 'tps'",
              sliderInput("tps_lambda", "TPS smoothness (λ)", min = 0, max = 1, value = 0.01, step = 0.01),
              helpText("Higher λ = smoother surface. 0 = interpolating (can overfit).")
            )
          )
        ),
        ),  # /Global Reference Plane and Surface
        tags$hr()
      )
    ),

    ### Output Panel ===
    mainPanel(
      # 2D plot
      conditionalPanel(
        condition = "input.mode3d == '2D'",
        plotlyOutput("scatter2D", height = "1000px", width = "100%")
      ),
      # 3D plot
      conditionalPanel(
        condition = "input.mode3d == '3D'",
        plotlyOutput("scatter3D", height = "1000px", width = "100%")
      ),
      verbatimTextOutput("fit_log")
    )
  )
),  # end tabPanel("Display")

    # ── Analysis ──────────────────────────────────────────────────────────────
    tabPanel("Analysis",
      tags$br(),

      # ── Dataset Groups (shared across analyses) ───────────────────────────
      wellPanel(
        h4("Dataset Groups"),
        radioButtons("analysis_by", "Analyze by",
                     choices = c("Individual files" = "file", "Dataset groups" = "group"),
                     selected = "file", inline = TRUE),
        conditionalPanel(
          condition = "input.analysis_by == 'group'",
          helpText("Assign loaded files to named groups. In group mode the analyses below pool files",
                   "within each group (one series per group) so you can compare groups. Files left",
                   "unassigned are ignored. Spatial statistics (Ripley/PCF) will average per-file curves",
                   "within a group rather than pooling coordinates (not yet wired in this pass)."),
          numericInput("n_groups", "Number of groups", value = 2, min = 1, max = 6, step = 1),
          conditionalPanel("input.n_groups >= 1", fluidRow(
            column(3, textInput("grp_name_1", "Group 1 name", "Group 1")),
            column(2, colourpicker::colourInput("grp_color_1", "Colour", value = "#1f77b4")),
            column(7, selectInput("grp_files_1", "Files", choices = NULL, multiple = TRUE)))),
          conditionalPanel("input.n_groups >= 2", fluidRow(
            column(3, textInput("grp_name_2", "Group 2 name", "Group 2")),
            column(2, colourpicker::colourInput("grp_color_2", "Colour", value = "#d62728")),
            column(7, selectInput("grp_files_2", "Files", choices = NULL, multiple = TRUE)))),
          conditionalPanel("input.n_groups >= 3", fluidRow(
            column(3, textInput("grp_name_3", "Group 3 name", "Group 3")),
            column(2, colourpicker::colourInput("grp_color_3", "Colour", value = "#2ca02c")),
            column(7, selectInput("grp_files_3", "Files", choices = NULL, multiple = TRUE)))),
          conditionalPanel("input.n_groups >= 4", fluidRow(
            column(3, textInput("grp_name_4", "Group 4 name", "Group 4")),
            column(2, colourpicker::colourInput("grp_color_4", "Colour", value = "#9467bd")),
            column(7, selectInput("grp_files_4", "Files", choices = NULL, multiple = TRUE)))),
          conditionalPanel("input.n_groups >= 5", fluidRow(
            column(3, textInput("grp_name_5", "Group 5 name", "Group 5")),
            column(2, colourpicker::colourInput("grp_color_5", "Colour", value = "#ff7f0e")),
            column(7, selectInput("grp_files_5", "Files", choices = NULL, multiple = TRUE)))),
          conditionalPanel("input.n_groups >= 6", fluidRow(
            column(3, textInput("grp_name_6", "Group 6 name", "Group 6")),
            column(2, colourpicker::colourInput("grp_color_6", "Colour", value = "#17becf")),
            column(7, selectInput("grp_files_6", "Files", choices = NULL, multiple = TRUE))))
        )
      ),

      # ── Batch NN ──────────────────────────────────────────────────────────
      fluidRow(
        column(3,
          h4("Batch NN Distance Histogram"),
          selectInput("nn_batch_unit", "NN unit",
            choices = c(
              "TID centroids (all files)"             = "tid",
              "Primary cluster centroids"             = "pcid",
              "TID centroids within Median Clusters"         = "cluster2_tid",
              "Primary cluster centroids in Median Clusters" = "cluster2_pcid"
            ),
            selected = "tid"
          ),
          helpText(
            "All modes run on every loaded file using current Display tab filters and clustering settings.",
            "Median Cluster modes require Median Clustering to be enabled."
          ),
          numericInput("nn_batch_k", "Nearest neighbors (k)", value = 1, min = 1, max = 20),
          uiOutput("nn_batch_group_ui"),
          shinyWidgets::numericRangeInput(
            "nn_batch_dist_range", "Distance range filter (nm)",
            value = c(0, 2000), min = 0
          ),
          helpText("Distances outside this range are excluded from the histogram and statistics."),
          radioButtons("nn_batch_yaxis", "Y axis",
                       choices = c("Count" = "count", "Density" = "density"),
                       selected = "count", inline = TRUE),
          numericInput("nn_batch_bins", "Bins", value = 60, min = 5, max = 500),
          selectInput("nn_batch_error", "Spread indicator",
            choices = c(
              "IQR  (Q25 – Q75)"  = "iqr",
              "MAD  (median ± MAD)" = "mad",
              "None"              = "none"
            ),
            selected = "iqr"
          ),
          checkboxInput("nn_batch_aggregate", "Aggregate all datasets into one", value = FALSE),
          conditionalPanel(
            condition = "output.mplex_active == true",
            tags$div(
              style = "background:#eef6ff; border:1px solid #b8daff; border-radius:4px; padding:6px 10px; margin-bottom:8px;",
              tags$b("Multiplex NN", style = "font-size:0.9em;"),
              selectInput("nn_batch_mplex_mode", "Mode",
                choices = c("All channels (one curve each)" = "split",
                            "Within channel X"              = "within",
                            "Between channel X and Y"       = "between"),
                selected = "split"),
              conditionalPanel(
                condition = "input.nn_batch_mplex_mode != 'split'",
                selectInput("nn_batch_chx", "Channel X", choices = character(0))
              ),
              conditionalPanel(
                condition = "input.nn_batch_mplex_mode == 'between'",
                selectInput("nn_batch_chy", "Channel Y", choices = character(0))
              ),
              helpText("Within: NN among points of channel X. Between: NN from channel X to channel Y.")
            )
          ),
          actionButton("run_batch_nn", "Run Batch NN",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_batch_nn", "Download distances", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("nn_batch_hist", height = "420px"),
          tableOutput("nn_batch_stats")
        )
      ),

      tags$hr(),

      # ── Traces per Cluster ────────────────────────────────────────────────
      fluidRow(
        column(3,
          h4("Traces per Cluster"),
          helpText("Counts distinct TIDs (trace IDs) per cluster across all loaded files.",
                   "In MINFLUX, one fluorophore can produce multiple TIDs across separate blinking cycles,",
                   "so this reflects blinking-event count, not molecular stoichiometry."),
          radioButtons("tid_count_unit", "Cluster type",
            choices = c(
              "Primary (Raw point clusters)" = "pcid",
              "Secondary (Median Clusters)" = "cluster2"
            ),
            selected = "cluster2"
          ),
          radioButtons("tid_count_yaxis", "Y axis",
                       choices = c("Count" = "count", "Density" = "density"),
                       selected = "count", inline = TRUE),
          numericInput("tid_count_bins", "Bins", value = 30, min = 5, max = 200),
          checkboxInput("tid_count_aggregate", "Aggregate all datasets into one", value = FALSE),
          conditionalPanel(
            condition = "output.mplex_active == true",
            checkboxInput("tid_count_by_channel", "Split by multiplex channel", value = FALSE)
          ),
          actionButton("run_tid_count", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_tid_count", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("tid_count_hist", height = "380px")
        )
      ),

      tags$hr(),

      # ── Cluster2 Detection over Time ──────────────────────────────────────
      fluidRow(
        column(3,
          h4("Median Clusters Detection over Time"),
          helpText("Each trace (TID) is counted once at its median acquisition time (tim),",
                   "split into traces inside Median Cluster structures vs unclustered background.",
                   "Reveals how detection of structures vs background evolves over the acquisition",
                   "(e.g. photobleaching drop-off) and feeds the simulator's acquisition model.",
                   "Requires Secondary Clustering enabled in the Display tab."),
          numericInput("det_time_bins", "Time bins", value = 40, min = 5, max = 300),
          radioButtons("det_time_yaxis", "Y axis",
                       choices = c("Count per bin"         = "count",
                                   "Cumulative"            = "cumulative",
                                   "Normalized rate"       = "norm",
                                   "Normalized cumulative" = "norm_cumulative"),
                       selected = "count"),
          helpText("Normalized modes scale each population (clustered / unclustered) to its own total,",
                   "so a difference in shape reveals detection timing that depends on clustering",
                   "(a crowding signature). The KS gap = max distance between the two cumulative curves."),
          checkboxInput("det_time_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          actionButton("run_det_time", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_det_time", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("det_time_plot", height = "400px")
        )
      ),

      tags$hr(),

      # ── Cluster2 Crowding vs Detection-Time Shift ─────────────────────────
      fluidRow(
        column(3,
          h4("Median Clusters Crowding vs Detection-Time Shift"),
          helpText("One point per Median Cluster structure.",
                   "Crowding = median nearest-neighbour distance among the cluster's trace centroids",
                   "(smaller = more crowded). Detection-time shift = the cluster's median detection time (tim)",
                   "minus the unclustered background median in the same file (positive = detected later).",
                   "If crowding delays detection, denser clusters (small NN distance) show a positive shift.",
                   "Note: this is a density-correlated timing effect; crowding is the likely but not sole cause.",
                   "Requires Secondary Clustering and a tim column."),
          numericInput("crowd_min_members", "Min traces per cluster", value = 4, min = 3),
          checkboxInput("crowd_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          checkboxInput("crowd_trend", "Show trend (loess)", value = TRUE),
          checkboxInput("crowd_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          actionButton("run_crowd", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_crowd", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("crowd_plot", height = "420px")
        )
      ),

      tags$hr(),

      # ── Within-Cluster2 NN vs Cluster Size ────────────────────────────────
      fluidRow(
        column(3,
          h4("Within Median Clusters NN vs Cluster Size"),
          helpText("One point per Median Cluster structure across all files.",
                   "Y = median over members of each member's mean distance to its k nearest member",
                   "centroids (k = 1 is nearest-neighbour spacing; centroids are trace-ID medians,",
                   "or primary-cluster medians if Primary Clustering is on).",
                   "Flat = constant spacing (grows at fixed density); rising = the structure 'balloons'",
                   "(members spread apart with size); falling = it densifies.",
                   "Note: under DBSCAN within-cluster NN is <= eps2 (the link radius, default 200 nm);",
                   "HDBSCAN/OPTICS impose no such bound. Requires Secondary Clustering."),
          radioButtons("size_nn_xmode", "X axis (cluster size)",
                       choices = c("Member count" = "count",
                                   "Convex-hull size (nm2 / nm3)" = "hull"),
                       selected = "count"),
          numericInput("size_nn_k", "Neighbours (k) for spacing", value = 1, min = 1),
          numericInput("size_nn_min_members", "Min members per cluster", value = 3, min = 2),
          checkboxInput("size_nn_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          checkboxInput("size_nn_logx", "Log x-axis (size)", value = FALSE),
          checkboxInput("size_nn_trend", "Show trend (loess)", value = TRUE),
          checkboxInput("size_nn_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          actionButton("run_size_nn", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_size_nn", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("size_nn_plot", height = "420px")
        )
      ),

      tags$hr(),

      # ── Cluster Viewer (2D / 3D + shape fit) ──────────────────────────────
      fluidRow(
        column(3,
          h4("Cluster Viewer and Shape Fitting"),
          helpText("Inspect individual Median Cluster structures in the active (plotted) file.",
                   "Press Run to load the current clusters, then step through them; view raw 3D or a",
                   "PCA-normalized 2D projection (PC1 vs PC2 — most differential axes; PCA runs in 3D",
                   "and reduces to 2D when z exists). Shape fitting is optional — enable it below.",
                   "Requires Secondary Clustering; press Plot in the Display tab first."),
          actionButton("run_viewer", "Run / Load clusters", icon = icon("play"), class = "btn-primary btn-sm"),
          tags$br(), tags$br(),
          fluidRow(
            column(6, actionButton("viewer_prev", "< Prev", class = "btn-sm btn-block")),
            column(6, actionButton("viewer_next", "Next >", class = "btn-sm btn-block"))
          ),
          tags$br(),
          tags$b(textOutput("viewer_header", inline = TRUE)),
          tags$br(), tags$br(),
          radioButtons("viewer_mode", "View",
                       choices = c("3D" = "3d", "2D PCA (PC1 vs PC2)" = "pca"),
                       selected = "3d", inline = TRUE),
          checkboxInput("viewer_enable_fit", "Enable shape fitting (MLE + best-fit overlay)", value = FALSE),
          conditionalPanel(
            condition = "input.viewer_enable_fit == true",
            radioButtons("viewer_rank", "Rank best fit by",
                         choices = c("AICc" = "aicc", "BIC" = "bic", "RMSD" = "rmsd"),
                         selected = "aicc", inline = TRUE),
            helpText("Fits are MLE (orthogonal-distance). RMSD always favours more-flexible shapes;",
                     "AICc/BIC penalise parameters. AICc is recommended for small clusters and shows",
                     "blank when a cluster has too few traces to justify that model (n <= params + 2)."),
            selectInput("viewer_fit_shape", "Fit overlay",
                        choices = c("Best (by ranking)" = "best", "None" = "none",
                                    "Sphere" = "sphere", "Ellipsoid" = "ellipsoid",
                                    "Helix" = "helix", "Ring" = "ring", "Arc" = "arc"),
                        selected = "best"),
            tableOutput("viewer_rmsd_table")
          )
        ),
        column(9,
          plotlyOutput("viewer_plot", height = "520px")
        )
      ),

      tags$hr(),

      # ── Radius of Gyration ────────────────────────────────────────────────
      fluidRow(
        column(3,
          h4("Radius of Gyration (RoG)"),
          helpText("RMS displacement of per-TID median positions from their cluster centroid.",
                   "Each TID contributes one point (its spatial median across all its localizations).",
                   "Runs on all loaded files using current clustering settings."),
          radioButtons("rog_cluster_col", "Cluster type",
            choices = c(
              "Primary (Raw point clusters)"   = "cluster",
              "Secondary (Median Clusters)" = "cluster2"
            ),
            selected = "cluster"
          ),
          radioButtons("rog_yaxis", "Y axis",
                       choices = c("Count" = "count", "Density" = "density"),
                       selected = "count", inline = TRUE),
          numericInput("rog_bins", "Bins", value = 30, min = 5, max = 200),
          checkboxInput("rog_aggregate", "Aggregate all datasets into one", value = FALSE),
          conditionalPanel(
            condition = "output.mplex_active == true",
            checkboxInput("rog_by_channel", "Split by multiplex channel", value = FALSE)
          ),
          actionButton("run_rog", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_rog", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("rog_hist", height = "380px")
        )
      ),

      tags$hr(),

      # ── Ripley's L ────────────────────────────────────────────────────────
      fluidRow(
        column(3,
          h4("Ripley's L(r)"),
          helpText("Spatial point-process statistic. L(r) − r = 0 under complete spatial randomness (CSR);",
                   "> 0 indicates clustering, < 0 indicates regularity.",
                   "No confidence envelope is shown — deviations should be interpreted relative to",
                   "the point density and window size, not as statistically tested. Uses an uncorrected",
                   "K estimator; curves may dip below zero at large r due to border effects.",
                   "Large datasets are subsampled to 2000 points per file."),
          selectInput("ripley_unit", "Unit",
            choices = c(
              "TID centroids"             = "tid",
              "Primary cluster centroids" = "cluster"
            ),
            selected = "tid"
          ),
          helpText("Unit: which points the statistic runs on — one centroid per trace ID, or one per",
                   "primary (Cluster1) cluster if primary clustering is enabled."),
          checkboxInput("ripley_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          helpText("All files: pool every loaded dataset into one curve; uncheck to analyse only the active file."),
          conditionalPanel(
            condition = "output.mplex_active == true",
            checkboxInput("ripley_by_channel", "Split by multiplex channel (with 'all files')", value = FALSE)
          ),
          checkboxInput("ripley_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          helpText("3-D: measure neighbours in spheres of radius r (uses z) instead of 2-D circles.",
                   "Use only for genuinely 3-D data; the CSR baseline differs between 2-D and 3-D."),
          numericInput("ripley_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          helpText("r max: largest spatial scale (radius) evaluated. 0 auto-picks roughly a quarter of the",
                   "field size. Read clustering peaks at the r where L(r)-r is highest = the cluster scale."),
          numericInput("ripley_r_steps", "r steps", value = 60, min = 10, max = 300),
          helpText("r steps: how many radii between 0 and r max are computed — higher = smoother curve, slower."),
          checkboxInput("ripley_smooth", "Smooth curve (loess)", value = TRUE),
          numericInput("ripley_smooth_span", "Smoothing span", value = 0.2, min = 0.05, max = 0.8, step = 0.05),
          helpText("Smoothing: a loess pass over the L(r)-r curve to reduce sampling noise. Span = fraction",
                   "of points in each local window (higher = smoother but can flatten real peaks). Purely",
                   "cosmetic — it does not change the underlying statistic."),
          conditionalPanel(
            condition = "input.analysis_by == 'group'",
            radioButtons("ripley_group_band", "Group spread band",
                         choices = c("IQR (Q25-Q75)" = "iqr", "±SD" = "sd", "None" = "none"),
                         selected = "iqr", inline = TRUE),
            helpText("In group mode each file's L(r)-r is computed on its own window, then averaged",
                     "within the group (one mean curve per group); the band shows between-file spread.",
                     "Curves are interpolated onto a common r-grid (up to the smallest per-file r max).",
                     "Read a group difference against the band, not as a formal test (no CSR envelope).")
          ),
          actionButton("run_ripley", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_ripley", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("ripley_plot", height = "400px")
        )
      ),
      tags$hr(),

      # ── Two-Color Analysis ────────────────────────────────────────────────
      fluidRow(
        column(12, h3("Two-Color Analysis"))
      ),

      conditionalPanel(
        condition = "output.mplex_active == true",
        fluidRow(
          column(12,
            tags$div(
              style = "background:#eef6ff; border:1px solid #b8daff; border-radius:4px; padding:8px 12px; margin-bottom:10px;",
              tags$b("Multiplex channel pair"),
              helpText("A multiplex dataset is active. Pick two channels to compare; the cross-channel tools below (NN, Ripley L₁₂, PCF g₁₂, KDE) run on this pair instead of the DCR/fluo 2-color split. Group colors come from the 2-color settings in the Display tab."),
              fluidRow(
                column(6, selectInput("mplex_cross_a", "Channel A (group 1)", choices = character(0))),
                column(6, selectInput("mplex_cross_b", "Channel B (group 2)", choices = character(0)))
              )
            )
          )
        )
      ),

      fluidRow(
        column(3,
          h4("Cross-Channel NN & Co-localization"),
          helpText("Requires '2-color MINFLUX' to be enabled in the Display tab.",
                   "Channel assignment and clustering (tc1/tc2) are taken from the 2-color settings.",
                   "Median  modes restrict comparisons to molecules within the same Cluster2 group."),
          selectInput("cross_nn_unit", "Unit",
            choices = c(
              "TID centroids"                            = "tid",
              "Primary cluster centroids"                = "cluster",
              "TID centroids within Median Clusters"            = "cluster2_tid",
              "Primary cluster centroids within Median Clusters"= "cluster2_cluster"
            ),
            selected = "tid"),
          numericInput("cross_nn_k", "k-th nearest neighbour", value = 1, min = 1, max = 10, step = 1),
          radioButtons("cross_nn_hist_type", "Histogram type",
                       choices = c("Count" = "count", "Density" = "density"),
                       selected = "density", inline = TRUE),
          numericInput("cross_nn_bins", "Bins", value = 40, min = 5, max = 200),
          checkboxInput("cross_nn_aggregate", "Aggregate all datasets into one", value = FALSE),
          actionButton("run_cross_nn", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_cross_nn", "Download CSV", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("cross_nn_hist", height = "340px"),
          plotlyOutput("coloc_curve",   height = "340px")
        )
      ),

      tags$hr(),

      fluidRow(
        column(3,
          h4("Bivariate Ripley's L₁₂(r)"),
          helpText("Requires '2-color MINFLUX' to be enabled in the Display tab.",
                   "L₁₂(r) − r > 0 indicates co-localization; < 0 indicates spatial exclusion.",
                   "Median Cluster modes compute L₁₂ within each Median Cluster group and average the curves.",
                   "Large datasets are subsampled to 2000 points per channel per file/group."),
          selectInput("ripley12_unit", "Unit",
            choices = c(
              "TID centroids"             = "tid",
              "Primary cluster centroids" = "cluster"
            ),
            selected = "tid"),
          checkboxInput("ripley12_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          numericInput("ripley12_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          numericInput("ripley12_r_steps", "r steps", value = 60, min = 10, max = 300),
          checkboxInput("ripley12_smooth", "Smooth curve (loess)", value = TRUE),
          numericInput("ripley12_smooth_span", "Smoothing span", value = 0.2, min = 0.05, max = 0.8, step = 0.05),
          checkboxInput("ripley12_aggregate", "Aggregate all datasets into one", value = FALSE),
          actionButton("run_ripley12", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_ripley12", "Download CSV", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("ripley12_plot", height = "400px")
        )
      ),
      tags$hr(),

      # ── Pair Cross-Correlation g₁₂(r) ────────────────────────────────────
      fluidRow(
        column(3,
          h4("Pair Cross-Correlation g₁₂(r)"),
          helpText("g₁₂(r) = 1 under CSR; > 1 indicates attraction (co-localization) at scale r;",
                   "< 1 indicates repulsion. Complements Ripley L₁₂ — shows the",
                   "interaction scale as a peak rather than a cumulative curve.",
                   "Requires '2-color MINFLUX' to be enabled in the Display tab."),
          selectInput("pcf12_unit", "Unit",
            choices = c("TID centroids" = "tid", "Primary cluster centroids" = "cluster"),
            selected = "tid"),
          checkboxInput("pcf12_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          numericInput("pcf12_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          numericInput("pcf12_r_steps", "r steps", value = 60, min = 10, max = 300),
          checkboxInput("pcf12_aggregate", "Aggregate all datasets into one", value = FALSE),
          actionButton("run_pcf12", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_pcf12", "Download CSV", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("pcf12_plot", height = "400px")
        )
      ),

      tags$hr(),

      # ── 2D KDE Density Co-localization Map ───────────────────────────────
      fluidRow(
        column(3,
          h4("2D KDE Density Map"),
          helpText("Kernel density estimate for each channel on the active file's spatial grid.",
                   "Ch1 shown as red contours, Ch2 as blue contours. Pearson r of the two",
                   "density surfaces serves as a global co-localization coefficient,",
                   "but is strongly bandwidth-dependent — a small bandwidth gives near-zero r",
                   "even for perfectly co-localizing molecules, while a large bandwidth inflates r.",
                   "Treat the Pearson r as relative rather than absolute.",
                   "Select the dataset to inspect via 'Active file' in the Display tab.",
                   "Requires '2-color MINFLUX' to be enabled in the Display tab."),
          selectInput("kde_unit", "Unit",
            choices = c("TID centroids" = "tid", "Primary cluster centroids" = "cluster"),
            selected = "tid"),
          uiOutput("kde_plane_ui"),
          selectInput("kde_style", "Display style",
            choices = c("Lines" = "lines", "Filled contours" = "fill", "Heatmap" = "heatmap"),
            selected = "lines"),
          selectInput("kde_colorpair", "Color pair",
            choices = c(
              "Red / Blue"       = "redblue",
              "Orange / Teal"    = "orangeteal",
              "Magenta / Green"  = "magentagreen",
              "Yellow / Purple"  = "yellowpurple"
            ), selected = "redblue"),
          sliderInput("kde_opacity", "Opacity", min = 0.2, max = 1.0, value = 0.75, step = 0.05),
          numericInput("kde_bandwidth", "Bandwidth (nm, 0 = auto)", value = 0, min = 0, step = 5),
          numericInput("kde_n_grid", "Grid resolution", value = 100, min = 20, max = 300, step = 10),
          actionButton("run_kde_coloc", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_kde_coloc", "Download CSV", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("kde_coloc_plot", height = "450px"),
          uiOutput("kde_pearson_display")
        )
      ),

    )  # end tabPanel("Analysis")

  )  # end tabsetPanel
)  # end fluidPage
