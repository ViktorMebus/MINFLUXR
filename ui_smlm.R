## ─────────────────────────────────────────────────────────────────────────────
## ui_smlm.R  —  SMLM analysis tab UI (assigned to smlm_ui)
## ─────────────────────────────────────────────────────────────────────────────

# Small orange "EXPLORATORY" chip reused by bespoke / not-literature-validated
# analysis panels (same visual language as the GT-estimation panel).
.exploratory_badge <- function() {
  tags$span("EXPLORATORY",
            style = "margin-left:8px;font-size:10px;font-weight:700;color:#fff;background:#e67e22;border-radius:3px;padding:1px 5px;vertical-align:middle;")
}

# NOTE: the former .exploratory_note() footnote helper was removed — the per-panel
# caveats it carried now live in the shared .method_box() (helpers_shared.R), so the
# Analysis tab no longer stated the same warning twice. The badge above still marks
# exploratory panels in their headings.

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
                    "Pool DCR across iterations (.json / .npy)",
                    value = TRUE),
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
        # Static (not a renderUI) so the disabled state applied by the planar-data
        # observer in server_smlm.R cannot be wiped by a re-render. XZ/YZ are
        # greyed rather than removed, so it stays visible WHY they are unavailable.
        radioButtons("view_plane_2d", "Viewing plane",
                     choices = c("XY", "XZ", "YZ"), inline = TRUE, selected = "XY")
      ),
      uiOutput("z_axes_note"),
      ),  # /Data Upload
      tags$hr(),

      
      ### Axial calibration ===
      # Deliberately OUTSIDE tags$div(id = "filters_form"): that div is wiped by
      # shinyjs::reset("filters_form") whenever a dataset changes, which would
      # snap the calibration back to 0.7 — and since applying it updates the
      # dataset, the two would ping-pong. It is also not a filter.
      tags$details(class = "sim-group", id = "z_scale_group",
      tags$summary(tags$b("Axial (z) calibration")),
      numericInput("z_scale", "Axial (z) scaling factor",
                   value = 0.7, min = 0.05, max = 3, step = 0.01),
      uiOutput("z_scale_status"),
      helpText("Takes effect on Plot Data. Loading a file that records its own factor ",
               "adopts that value."),
      ),  # /Axial calibration
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
          helpText("Wraps each trace (or each de-blinked ID when De-blinking is enabled) in a 3D alpha ",
                   "surface."),
          conditionalPanel(
            condition = "input.show_convex_hull == true",
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

            ### De-blink correction (shared) ===
      tags$details(class = "sim-group",
      tags$summary(tags$b("De-blink correction")),
      checkboxInput("enable_clustering", "Enable de-blink correction", value = FALSE),
      helpText("Merges a fluorophore's repeated blinks (nearby localizations / trace IDs) into one",
               "de-blinked molecule, so downstream analyses count molecules rather than blinking events."),
      conditionalPanel(
        condition = "input.enable_clustering == true",
        radioButtons("deblink_method", "Merge criterion",
                     choices = c("Trace centroid distance" = "centroid",
                                 "Majority vote of localizations" = "vote"),
                     selected = "centroid"),
        helpText("Centroid: merge when two traces sit in the same place. Majority vote: merge when ",
                 "most of a trace's localizations are closer to a neighbouring trace than to their ",
                 "own. Centroid is more accurate; the vote needs no distance in nm."),

        conditionalPanel(
          condition = "input.deblink_method == 'centroid'",
          numericInput("eps_centroid", "Centroid distance threshold (nm)", value = 10, min = 0.1),
          helpText("Two traces merge when their median positions are closer than this. ",
                   "This is the value the merge actually uses; k below just rewrites it in ",
                   "units of the measured trace spread."),
          numericInput("centroid_k", "Threshold in units of trace spread (k)", value = 2.0,
                       min = 0.5, max = 6, step = 0.1),
          helpText("Changing k rewrites the distance above. Too low leaves re-blinks counted ",
                   "twice; too high fuses neighbouring molecules."),
          actionButton("suggest_eps_centroid", "Measure threshold from data", class = "btn-sm"),
          helpText("Sets k and the threshold from your own data."),
          uiOutput("centroid_suggest_summary"),
          checkboxInput("view_centroid_gap", "Show centroid threshold diagnostic", value = FALSE),
          conditionalPanel(
            condition = "input.view_centroid_gap == true",
            plotOutput("centroid_gap_plot", height = "280px"),
            helpText("The orange curve has two humps: re-blinks of one fluorophore (left) and ",
                     "separate molecules (right). Your threshold should sit in the gap between them.")
          ),
          checkboxInput("view_deblink_knee", "Molecule count vs threshold (exploratory)", value = FALSE),
          conditionalPanel(
            condition = "input.view_deblink_knee == true",
            plotOutput("deblink_knee_plot", height = "280px"),
            uiOutput("deblink_knee_summary"),
            helpText("How the molecule count falls as the threshold grows. Steep while re-blinks ",
                     "are being collected, then slow as separate molecules start fusing.")
          )
        ),
        conditionalPanel(
          condition = "input.deblink_method == 'vote'",
          numericInput("vote_frac", "Fraction of localizations that must agree (f)",
                       value = 0.2, min = 0.02, max = 0.9, step = 0.05),
          helpText("Two traces of one fluorophore score about 0.5, so useful values are well ",
                   "below that. Scores fall to ~0.23 at 10 nm apart and ~0.07 at 20 nm."),
          numericInput("vote_rho", "Distance tolerance (rho)", value = 1.2, min = 0.5, max = 3, step = 0.1),
          helpText("1.0 = the neighbouring trace must be strictly closer; higher is more willing to merge."),
          radioButtons("vote_link", "Agreement", inline = TRUE,
                       choices = c("Both traces must agree" = "mutual", "Either is enough" = "either"),
                       selected = "mutual")
        ),

        # ── anisotropic metric (see deblink_axis_weights in helpers_smlm.R) ──
        checkboxInput("deblink_anisotropic",
                      "Scale axes by their own spread (anisotropic radius)", value = FALSE),
        helpText("Weights x, y and z by their measured trace spread before merging, so one radius ",
                 "fits all three axes. Affects the merge decision only — reported distances, ",
                 "cluster sizes and Rg stay in true nm."),


        helpText("TIDs are not split: all localizations of a trace ID join the same de-blinked molecule (majority vote). ",
                 "Trace IDs that overlap spatially are merged into one de-blinked molecule. ",
                 "De-blinking only merges — a trace it cannot merge becomes its own single-blink molecule, ",
                 "so no molecule is removed here. Removal is the job of the metadata and spread filters.")
      ),
      ),  # /De-blink correction
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
          checkboxInput("nn_within_cluster2", "Within Median Cluster only", value = FALSE)
        ),
        downloadButton("download_nn_table", "Download NN Distance Table")
      ),
      ),  # /Nearest Neighbor Analysis
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
            helpText("Fixed: centre = convex-hull centroid; radius solved analytically (no iteration). Regularized: centre is free but penalized for moving away from the robust centre (Gaussian MAP prior, strength \u03bb). Free: centre fully optimized within an expanded data bounding box."),
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
            helpText("Used only to seed the pitch estimate. The fitted turn count is determined by the data extent and the optimized pitch \u2014 it will not be exactly this value.")
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
          uiOutput("dcr_gmm_info"),
          uiOutput("use_pooled_dcr_ui"),
          tags$details(style = "margin:6px 0 4px; font-size:0.85em;",
            tags$summary(style = "cursor:pointer; color:#2c3e50; font-weight:600;",
                         "Per-trace split confidence"),
            helpText(style = "font-size:0.95em; color:#666; margin:4px 0;",
              "Posterior probability that each trace belongs to the group it was assigned to, ",
              "from the fitted mixture. 1.0 = unambiguous; 0.5 = the molecule sits exactly on ",
              "the split and is a coin flip. ", tags$b("Orange = on the edge"),
              " (posterior below the cut) — these are the molecules that would change channel ",
              "under a small shift of the threshold, and the ones that manufacture false ",
              "short-range cross-channel pairs. Use the export to exclude or inspect them."),
            plotlyOutput("dcr_conf_plot", height = "260px"),
            downloadButton("download_dcr_conf", "Download per-trace confidence", class = "btn-sm"))
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
        checkboxInput("enable_two_color_cluster1", "2-color de-blink correction", FALSE),
        conditionalPanel(
          condition = "input.enable_two_color_cluster1",
          radioButtons("deblink_method_tc1", "2-color merge criterion",
                       choices = c("Trace centroid distance" = "centroid",
                                   "Majority vote of localizations" = "vote"),
                       selected = "centroid"),
          conditionalPanel(
            condition = "input.deblink_method_tc1 == 'centroid'",
            numericInput("eps_tc1", "2-color centroid distance threshold (nm)", value = 12, min = 0.1),
            helpText("About 2\u00d7 the trace spread.")
          ),
          conditionalPanel(
            condition = "input.deblink_method_tc1 == 'vote'",
            numericInput("vote_frac_tc1", "Fraction of localizations that must agree (f)",
                         value = 0.2, min = 0.02, max = 0.9, step = 0.05),
            numericInput("vote_rho_tc1", "Distance tolerance (rho)",
                         value = 1.2, min = 0.5, max = 3, step = 0.1)
          ),
          actionButton("suggest_eps_tc1", "Measure threshold from data", class = "btn-sm"),
          helpText("Measures the trace spread pooled over both channels and sets the threshold. ",
                   "Merging still happens strictly within a channel."),
          uiOutput("tc1_suggest_summary"),
          checkboxInput("view_tc1_gap", "Show 2-color threshold diagnostic", value = FALSE),
          conditionalPanel(
            condition = "input.view_tc1_gap == true",
            plotOutput("tc1_gap_plot", height = "280px"),
            helpText("The orange curve has two humps: re-blinks of one fluorophore (left) and ",
                     "separate molecules (right). Your threshold should sit in the gap between them.")
          ),
          helpText("De-blinked IDs carry their channel: A-1, A-2 for one fluorophore and ",
                   "B-1, B-2 for the other. Only members of the same channel can merge, so ",
                   "two colours can never be de-blinked into one molecule."),
          conditionalPanel(
            condition = "input.mode3d == '3D'",
            checkboxInput("show_tc1_alpha", "Wrap 2-color de-blinked molecules", value = FALSE),
            conditionalPanel(
              condition = "input.show_tc1_alpha == true && input.mode3d == '3D'",
              numericInput("tc1_alpha", "Alpha (smaller = tighter)", value = 120, min = 1, step = 1),
              sliderInput("tc1_alpha_opacity", "Alpha shape opacity", min = 0, max = 1, value = 0.15, step = 0.01),
              numericInput("tc1_alpha_min_pts", "Min points per cluster to wrap", value = 5, min = 2, step = 1)
            ),
          ),
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


      ### Estimate ground-truth positions (labeling-aware inverse) ===
      tags$details(class = "sim-group",
        tags$summary(tags$b("Estimate Ground-Truth Positions"),
                     tags$span("EXPLORATORY",
                               style = "margin-left:8px;font-size:10px;font-weight:700;color:#fff;background:#e67e22;border-radius:3px;padding:1px 5px;vertical-align:middle;")),
        helpText("Exploratory, model-based estimate (not a direct measurement). Reverse-engineers the labeling model to turn traces / median-clusters back into candidate ground-truth positions, an approximate count, and a crowding flag. Hollow circles are coloured green (low) → red (high) by position uncertainty. See the method/reference notes below the read-outs."),
        checkboxInput("show_gt_estimate", "Show estimated GT candidates", value = FALSE),
        conditionalPanel(
          condition = "input.show_gt_estimate == true",
          selectInput("gt_label_preset", "Labeling method (same presets as the Simulator)",
                      choices  = names(minflux_label_presets()),
                      selected = "Primary monoclonal IgG + secondary IgG"),
          conditionalPanel(
            condition = "input.gt_label_preset == 'Custom'",
            numericInput("gt_fluor_custom",     "Fluorophores / molecule (mean)",     value = 1,  min = 1, step = 1),
            numericInput("gt_link_custom",      "Linker length (nm)",                 value = 5,  min = 0, step = 0.5),
            numericInput("gt_spread_custom",    "Dye spread within complex (nm)",     value = 0,  min = 0, step = 0.5),
            numericInput("gt_footprint_custom", "Footprint / exclusion diam. (nm)",   value = 8,  min = 0, step = 1),
            numericInput("gt_blink_custom",     "Extra blinks / fluorophore (mean)",  value = 0.3, min = 0, step = 0.1)
          ),
          sliderInput("gt_label_eff", "Labeling / detection efficiency (%)",
                      min = 1, max = 100, value = 70, step = 1),
          radioButtons("gt_model_level", "Model from",
                       choices = c(
                         "Auto (de-blinked molecules if present, else raw traces)" = "auto",
                         "De-blinked molecules"                 = "primary",
                         "Raw trace IDs"                                 = "tid"),
                       selected = "auto"),
          checkboxInput("gt_scope_cluster2",
                        "Restrict to points inside Median Clusters", value = TRUE),
          helpText("Auto is recommended: de-blinking collapses re-blinks and multivalent labels to ~1 unit per molecule, turning an unstable divide-by-multiplicity count into a divide-by-~1 count."),
          numericInput("gt_complex_r",
                       "Label-complex merge radius (nm, blank = auto)",
                       value = NA, min = 0, step = 5),
          helpText("Traces / de-blinked molecules within this radius are pooled into ONE candidate GT (a multivalent label scatters many fluorophores around one molecule). Blank = auto, set per labeling preset (max of footprint and 2.5 × dye spread) — big for polyclonal+secondary, tiny for a nanobody; 0 = no merge (one unit = one GT)."),
          checkboxInput("gt_guided_merge",
                        "Stoichiometry-guided merge (model-driven)", value = FALSE),
          helpText("MODEL-DRIVEN: where the spatial merge leaves more GTs than the labeling stoichiometry predicts (TIDs ÷ expected TIDs-per-molecule), collapse further to that target. Only ever reduces the count toward stoichiometry, never splits. Use when the spatial and stoichiometric counts disagree — but treat the result as prior-driven, and always sanity-check against the two-count readout below."),
          checkboxInput("gt_displace_inward",
                        "Displace GT inward to structure (linker correction)", value = TRUE),
          numericInput("gt_displace_nm",
                       "Inward displacement (nm, blank = auto = linker length)",
                       value = NA, min = 0, step = 1),
          helpText("The antibody binds a surface epitope and the linker throws the dye OUTWARD, so the trace cloud sits ~linker-length outside the true structure. This pushes each GT back toward the local structure centre (perpendicular to a filament's axis; toward the centroid for a compact/ring shape). Needs Median Clusters defined to know the structure."),
          checkboxInput("gt_show_lines", "Draw contribution lines (unit → GT)", value = TRUE),
          sliderInput("gt_circle_size", "GT marker size", min = 4, max = 40, value = 15, step = 1),
          sliderInput("gt_support_strength",
                      "Support-penalty strength (position certainty)",
                      min = 0, max = 2, value = 1, step = 0.1),
          helpText("How hard a GT backed by few trace IDs is penalized in the uncertainty (redder). 0 = off (only the linker floor + sampling); 1 = default (a 1-TID GT ≈ √2× less certain than a well-supported one); lower it if legitimately sparse labelling reads too red."),
          uiOutput("gt_estimate_readout"),
          tags$hr(style = "margin:8px 0;"),
          tags$b("Structure shape & size"),
          checkboxInput("gt_show_shape",
                        "Classify shape + measure size (contour, extent, spacing)", value = FALSE),
          conditionalPanel(
            condition = "input.gt_show_shape == true",
            helpText("Infers each cluster's shape from the point distribution (filament / ring / punctum / sheet / diffuse) and reports neutral size descriptors: contour length, straight extent, width, radius of gyration, and the spacing between estimated points."),
            numericInput("gt_subunit_spacing",
                         "Assumed subunit spacing (nm, optional — blank = off)",
                         value = NA, min = 0.1, step = 0.5),
            helpText("Optional and generic: if your structure has a known periodic subunit spacing, a filament's contour length is divided by it to give an estimated subunit count. Leave blank for pure geometry."),
            uiOutput("gt_morphology_readout")
          )
        )
      ),
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


      # ── SECTION: Quality Control ──
      fluidRow(column(12, tags$h3('Quality Control',
        style = "border-bottom:2px solid #d7dde3; padding-bottom:4px; margin-top:18px; color:#2c3e50;"))),

      # ── QC: Localization precision cloud ──────────────────────────────────
      fluidRow(
        column(3,
          h4("Localization Precision Cloud"),
          .method_box("green",
            HTML("The direct measure of localization precision. Every localization is plotted relative to its own trace centroid and pooled at (0,0); the RMS radius is σ_loc, drawn as the 1σ and 2σ ellipses (per-axis σ<sub>x</sub>, σ<sub>y</sub>, σ<sub>z</sub>)."),
            HTML("Pooled within-trace residuals average out real motion and keep only the measurement scatter, so the cloud IS the precision. Uses the <b>filtered</b> localizations — the same set the Display tab plots (metadata / length / spread filters and any fiducial exclusion applied) — so it reports the precision of the data you actually keep."),
            HTML("The <b>shape</b> is the QC: a compact isotropic disc is good; <b>elongation</b> (anisotropy &gt; ~1.3) points to drift within a trace or an anisotropic PSF; <b>heavy tails</b> (kurtosis ≫ 3) point to mis-links / double events; a net <b>skew</b> points to a systematic bias.<br><br>The frame is zoomed to the <b>robust</b> (per-trace-median) precision so the real cloud fills it. With <b>Include z</b> on, three panels are shown: XY, XZ and YZ.")),
          helpText("Pooled scatter of every localization about its trace centroid. Reports σ_x/σ_y",
                   "(σ_z), anisotropy, skew and kurtosis. Uses the filtered (Display-tab) data."),
          helpText(style = "font-size:0.78em; color:#888;",
                   "Always assesses the active dataset only — precision is never pooled across files."),
          numericInput("prec_cloud_min_loc", "Min localizations / trace", value = 5, min = 3, step = 1),
          checkboxInput("prec_cloud_use_z", "Include z (report σ_z)", value = FALSE),
          actionButton("run_prec_cloud", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_prec_cloud", "Download σ", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("prec_cloud_plot", height = "340px"),
          uiOutput("prec_cloud_stats")
        )
      ),

      tags$hr(),


      # ── QC: Fiducial / bright-contaminant flagging (report only) ──────────
      fluidRow(
        column(3,
          h4("Fiducial / Contaminant Scan"),
          .method_box("amber",
            HTML("Flags traces that behave like <b>gold fiducials or bright aggregates</b> rather than single blinking molecules, using the two signatures that actually separate them (learned from a labelled fiducial-vs-clean dataset pair): <b>BRIGHT</b> — efo is an extreme outlier (a gold bead reflects ~100× a dye, efo ~10⁷ vs ~10⁵); and <b>PERSISTENT + RECURRENT</b> — a long trace that also <b>piles up at the same spot</b> as many others (fiducials recur in one place; real labels are spread out)."),
            HTML("Scores are <b>per-file and robust</b> (median + k·MAD on log, so a heavy fiducial section can't hide itself, and a dim file is judged against itself). <b>Spread is reported but not used to flag</b> — it does not distinguish fiducials from real structures. A trace is flagged if it is bright, OR long-and-recurrent, OR extremely long."),
            HTML("<b>Report by default</b> — it tabulates suspects; nothing is removed unless you tick the box below. On a clean file this flags essentially nothing; treat the list as candidates to inspect.")),
          helpText("Robust brightness / length outliers + spatial recurrence. Review the scatter",
                   "and table, then optionally exclude the flagged traces."),
          numericInput("fid_efo_z",    "Brightness outlier (robust z ≥)", value = 4, min = 2, max = 10, step = 0.5),
          numericInput("fid_len_z",    "Length outlier (robust z ≥)",     value = 4, min = 2, max = 10, step = 0.5),
          numericInput("fid_recur_r",  "Recurrence radius (nm)",          value = 100, min = 10, max = 1000, step = 10),
          numericInput("fid_recur_min","Min co-located traces",           value = 4, min = 1, max = 50, step = 1),
          checkboxInput("fid_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          checkboxInput("fid_exclude", "Exclude flagged traces from analyses", value = FALSE),
          helpText(style = "margin-top:-6px; font-size:0.78em; color:#888;",
                   "Non-destructive & reversible: drops the flagged traces (per file) from the",
                   "display and the Analysis tools while ticked. Precision Cloud and this scan",
                   "always see the raw data. Re-run the scan after changing thresholds."),
          actionButton("run_fiducial", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_fiducial", "Download per-trace table", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("fiducial_plot", height = "360px"),
          div(style = "max-height:220px; overflow-y:auto; font-size:0.82em; margin-top:6px;",
              tableOutput("fiducial_table"))
        )
      ),

      tags$hr(),

      # ── Trace Detection over Time (unique TIDs vs acquisition time) ────────
      fluidRow(
        column(3,
          h4("Trace Detection over Time"),
          .method_box("green",
            HTML("Counts <b>unique trace IDs</b> per acquisition-time bin. A trace ID is one continuous detection of a blinking fluorophore, so this is the <b>detection-event rate</b>."),
            HTML("Each trace is counted once, at its median <code>tim</code>, so a long trace does not dominate the early bins. Multiple files are drawn as <b>one series each</b> (or pooled per group in group mode). Each dataset has its own tim timeline, so a line is drawn only across that dataset's own window and <b>stops when its tim ends</b> — it does not drop to zero. The per-localization metadata filters apply; the trace-length filter is deliberately not applied, because short blinking events are real detections."),
            HTML("Descriptive QC, not an estimator. A declining rate over time reflects photobleaching, exhausted blinking, drift and label/sample depletion — read it as “is my acquisition stable?”, not as a single physical quantity.")),
          helpText("Unique traces (TIDs) counted once at their median acquisition time (tim),",
                   "as one series per file — or per Dataset Group when group mode is on above.",
                   "Useful for detection-rate and acquisition-stability QC."),
          tags$b("Trace-ID length prefilter", style = "font-size:0.9em;"),
          fluidRow(
            column(6, numericInput("tid_time_len_min", "Min locs/trace", value = 1, min = 1, step = 1)),
            column(6, numericInput("tid_time_len_max", "Max (0 = none)", value = 0, min = 0, step = 1))
          ),
          helpText("Keeps only trace IDs whose localization count falls in this range,",
                   "applied identically to every dataset before counting."),
          numericInput("tid_time_bins", "Time bins", value = 40, min = 5, max = 300),
          radioButtons("tid_time_yaxis", "Y axis",
                       choices = c("Count per bin"         = "count",
                                   "Cumulative"            = "cumulative",
                                   "Normalized rate"       = "norm",
                                   "Normalized cumulative" = "norm_cumulative"),
                       selected = "count"),
          helpText("Normalized modes scale each series to its own total, so files/groups with",
                   "different acquisition lengths or brightness are compared on detection timing (shape)",
                   "rather than absolute counts (amount)."),
          checkboxInput("tid_time_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          actionButton("run_tid_time", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_tid_time", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("tid_time_plot", height = "400px")
        )
      ),

      tags$hr(),


      # ── SECTION: Nearest-Neighbour & Spatial Statistics ──
      fluidRow(column(12, tags$h3('Nearest-Neighbour & Spatial Statistics',
        style = "border-bottom:2px solid #d7dde3; padding-bottom:4px; margin-top:18px; color:#2c3e50;"))),

      # ── Batch NN ──────────────────────────────────────────────────────────
      fluidRow(
        column(3,
          h4("Batch NN Distance Histogram"),
          .method_box("green",
            HTML("Distribution of the k-th nearest-neighbour distance, computed with a k-d tree (<code>dbscan::kNN</code>) — the classic first-order description of how tightly points sit together. Nearest-neighbour distance as a spatial statistic follows <b>Clark &amp; Evans 1954</b>; their observed-to-Poisson-expected ratio R is not computed here, as this panel reports the distance distribution itself with no null model attached (see the caveat)."),
            HTML("The choice of <b>NN unit</b> is the MINFLUX adaptation. Raw localizations are <i>not</i> molecules — one fluorophore re-appears as several traces and a multivalent label scatters dyes ~10–40&nbsp;nm around one epitope. Running NN on <b>TID centroids</b> (one point per trace) or <b>de-blinked-molecule / median-cluster centroids</b> de-blinks first, so the distribution approximates inter-<i>molecule</i> spacing. MINFLUX's nm precision means the histogram is limited by labelling geometry, not by localization error."),
            HTML("A sharp peak below ~40&nbsp;nm is usually the <b>label footprint, not biology</b> (Annibale et&nbsp;al. 2011; Sengupta et&nbsp;al. 2011). The distance-range filter trims the histogram but does not correct the bias. NN is a summary statistic with no CSR reference here — use Ripley&nbsp;L(r) with its envelope if you need significance.")),
          selectInput("nn_batch_unit", "NN unit",
            choices = c(
              "TID centroids (all files)"             = "tid",
              "De-blinked molecule centroids"             = "pcid",
              "TID centroids within Median Clusters"         = "cluster2_tid",
              "De-blinked molecule centroids within Median Clusters" = "cluster2_pcid"
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

      # ── Ripley's L ────────────────────────────────────────────────────────
      fluidRow(
        column(3,
          h4("Ripley's L(r)"),
          .method_box("green",
            HTML("Ripley's K-function (<b>Ripley 1977</b>) with the variance-stabilising L-transform (<b>Besag 1977</b>). For <i>n</i> points in a window of area <i>A</i>, the (uncorrected) estimator is<br><span style=\"font-family:serif; margin-left:6px;\">K̂(r) = (A / [n(n−1)]) · Σ<sub>i≠j</sub> 𝟙(d<sub>ij</sub> ≤ r)</span><br>i.e. area × (fraction of the n(n−1) ordered point pairs closer than r). The L-transform then linearises it: <span style=\"font-family:serif;\">L(r) = √(K(r)/π)</span> in 2D, <span style=\"font-family:serif;\">L(r) = (3·K(r)/4π)<sup>1/3</sup></span> in 3D. Plotted as L(r)&nbsp;−&nbsp;r: 0 under complete spatial randomness, &gt;0 clustering, &lt;0 dispersion. The SMLM lineage of this statistic runs through PC-PALM (<b>Sengupta et&nbsp;al. 2011</b>; <b>Veatch et&nbsp;al. 2012</b>)."),
            HTML("Runs on <b>de-blinked units</b> (TID or primary-cluster centroids) so that repeated appearances of one molecule do not manufacture a short-r clustering peak. Points are randomly subsampled to ≤2000 per file (pairwise distances are O(n²)); the window is the axis-aligned bounding box. The optional <b>Monte-Carlo CSR envelope</b> simulates n random points in the <i>identical frozen window</i>, so the band and the observed curve carry the same edge bias and it cancels."),
            HTML("There is <b>no analytic edge correction</b>: near the top of the r range points close to the box lose neighbours and L−r bends downward. Read significance as the observed curve leaving the envelope, <b>not</b> as departure from the flat zero line (which is the uncorrected theoretical CSR). Without the envelope, read shape and peak position, not magnitude. The dominant biological confound remains labelling: incomplete labelling flattens real clustering, multivalent labels fake it. In group mode per-file curves are averaged (never pooled coordinates) — pooling would mix windows and invent structure.")),
          helpText("Spatial point-process statistic. L(r) − r = 0 under complete spatial randomness (CSR);",
                   "> 0 indicates clustering, < 0 indicates regularity.",
                   "No confidence envelope is shown — deviations should be interpreted relative to",
                   "the point density and window size, not as statistically tested. Uses an uncorrected",
                   "K estimator; curves may dip below zero at large r due to border effects.",
                   "Large datasets are subsampled to 2000 points per file."),
          selectInput("ripley_unit", "Unit",
            choices = c(
              "TID centroids"             = "tid",
              "De-blinked molecule centroids" = "cluster"
            ),
            selected = "tid"
          ),
          helpText("Unit: which points the statistic runs on — one centroid per trace ID, or one per",
                   "the de-blinked molecule if de-blinking is enabled."),
          checkboxInput("ripley_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          helpText("All files: pool every loaded dataset into one curve; uncheck to analyze only the active file."),
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
          checkboxInput("ripley_envelope", "Show CSR envelope (Monte-Carlo)", value = FALSE),
          helpText("Simulates random (CSR) patterns in the same window as the data and shades the",
                   "2.5–97.5% band. The observed curve leaving the band = clustering (above) or",
                   "dispersion (below) beyond chance. The band captures the edge bias, so it — not",
                   "the flat zero line — is the correct reference. In 'all files' mode the envelope",
                   "is computed for the first few files only (cost = simulations × files)."),
          conditionalPanel(
            condition = "input.ripley_envelope == true",
            numericInput("ripley_nsim", "Envelope simulations", value = 99, min = 19, max = 199)
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


      # ── SECTION: Median-Cluster Analyses ──
      fluidRow(column(12, tags$h3('Median-Cluster Analyses',
        style = "border-bottom:2px solid #d7dde3; padding-bottom:4px; margin-top:18px; color:#2c3e50;"))),

      # ── Traces per Cluster ────────────────────────────────────────────────
      fluidRow(
        column(3,
          h4("Trace IDs / De-blinked Molecules per Median Cluster"),
          .method_box("amber",
            HTML("Cluster-occupancy counting, in the lineage of SMLM molecular-counting / blink-correction work (<b>Annibale et&nbsp;al. 2011</b>, <i>Nat. Methods</i> — quantitative photo-activated localization microscopy; <b>Sengupta et&nbsp;al. 2011</b>, PC-PALM). The number of units the median cluster contains."),
            HTML("The grouping is always the <b>median cluster</b> (secondary DBSCAN structure). You choose what to count inside it: raw <b>trace IDs</b>, or <b>de-blinked molecules</b> (which collapse one molecule's repeated blinks into one unit — needs de-blinking on). Comparing the two shows how much of the raw count is blinking."),
            HTML("<b>This is an occupancy count, not stoichiometry.</b> A trace-ID count over-counts molecules by an unknown, dye- and buffer-dependent blinking factor; de-blinked molecules partly correct for that but still depend on the de-blink eps. Annibale's blink-correction requires a dark-time distribution that is not calibrated here. Read it as a relative comparison between conditions acquired identically, never as an absolute subunit count.")),
          helpText("Counts, per median cluster, either distinct trace IDs or distinct primary",
                   "(de-blinked) clusters, across all loaded files. In MINFLUX one fluorophore",
                   "produces several trace IDs across blinking cycles, so a trace-ID count reflects",
                   "blinking events; de-blinked molecules are closer to molecule count."),
          radioButtons("tid_count_unit", "Count within each Median Cluster",
            choices = c(
              "Trace IDs (blinking events)"          = "tid",
              "De-blinked molecules (primary DBSCAN)" = "pcid"
            ),
            selected = "tid"
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
          h4("Median Clusters Detection over Time", .exploratory_badge()),
          .method_box("orange",
            HTML("<b>MINFLUXR original.</b> No published MINFLUX detection-timing estimator exists. The building blocks are standard: a detection-rate histogram over acquisition time, and a two-sample <b>Kolmogorov–Smirnov</b> comparison (max distance between two empirical cumulative distributions) of the clustered vs unclustered populations."),
            HTML("Each trace is counted <b>once</b>, at its median <code>tim</code>, so a long trace with many localizations does not dominate the early bins. The split is clustered (inside a median-cluster structure) vs unclustered background <i>within the same file</i>, so acquisition length, laser power and sample differences cancel in the comparison. The normalized modes scale each population to its own total, isolating a difference in <i>shape</i> (timing) from a difference in <i>amount</i> (abundance)."),
            HTML("<b>Descriptive, not causal.</b> The drop-off is a convolution of photobleaching, blinking kinetics, labelling density and the tracker's own re-targeting — crowding is one candidate cause among several. The KS gap is reported as a distance, <b>not as a p-value</b>: the two samples are not independent (they come from the same acquisition and the same molecules compete for the same detection budget), so the KS null does not hold. The empirical trace-detection-rate decay used to calibrate the simulator comes from in-house long acquisitions, not from a publication.")),
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
          h4("Median Clusters Crowding vs Detection-Time Shift", .exploratory_badge()),
          .method_box("orange",
            HTML("<b>MINFLUXR original.</b> A correlational diagnostic built from two standard quantities: a local-density proxy (median k-NN distance among a cluster's de-blinked molecule centroids) and a per-cluster timing offset, with an optional loess trend (<b>Cleveland 1979</b>)."),
            HTML("Both axes are measured <b>per cluster and referenced within the same file</b> — the timing shift subtracts that file's unclustered-background median <code>tim</code>, so between-file differences in acquisition length, bleaching and labelling cancel. Crowding uses <b>de-blinked molecule centroids</b> (or trace-ID centroids when de-blinking is off) — the <i>same</i> unit and k-NN definition as the Within-cluster-NN panel, so the two panels report the same spacing."),
            HTML("<b>A correlation, not a mechanism.</b> Local density co-varies with labelling efficiency, structure size and local photophysics, so a crowding–timing trend is a hypothesis to test with controls (e.g. a labelling-density titration), not a measured crowding effect. The sampling artefact — denser clusters have more members, and a larger sample's median lands later under a decaying detection rate, so a positive trend is partly expected with <i>no</i> crowding — is <b>now corrected for</b> by the permutation null band and its Spearman trend test. What the test cannot rule out is a third variable driving both axes; that still needs a control experiment.")),
          helpText("One point per Median Cluster structure.",
                   "Crowding = median k-NN distance among the cluster's de-blinked molecule centroids",
                   "(or trace-ID centroids if de-blinking is off; smaller = more crowded), ± MAD.",
                   "Detection-time shift = the cluster's median detection time (tim) minus the unclustered",
                   "background median in the same file (positive = detected later).",
                   "If crowding delays detection, denser clusters (small NN distance) show a positive shift.",
                   "Note: this is a density-correlated timing effect; crowding is the likely but not sole cause.",
                   "Requires Median Clustering and a tim column."),
          numericInput("crowd_min_members", "Min members per cluster", value = 4, min = 3),
          numericInput("crowd_k", "Nearest neighbours (k)", value = 1, min = 1),
          checkboxInput("crowd_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          checkboxInput("crowd_trend", "Show trend (loess)", value = TRUE),
          checkboxInput("crowd_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          checkboxInput("crowd_null", "Show no-crowding null (label permutation)", value = TRUE),
          conditionalPanel(
            condition = "input.crowd_null == true",
            numericInput("crowd_nsim", "Permutations", value = 200, min = 49, max = 2000, step = 50),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "Shuffles detection times among traces within each file, keeping every cluster's ",
              tags$b("membership, size and geometry fixed"), ". This reproduces the part of the ",
              "time-shift that comes purely from cluster size (a larger sample's median lands later ",
              "under a decaying detection rate) while destroying any real crowding-timing link — so ",
              "the shaded band is what you would see with ", tags$b("zero"), " crowding effect. ",
              "Read the points against the band, and the reported Spearman p for the ",
              tags$b("trend"), " (which is insensitive to a constant clustered-vs-background offset).")
          ),
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
          h4("Within Median Clusters NN vs Cluster Size", .exploratory_badge()),
          .method_box("orange",
            HTML("<b>MINFLUXR original.</b> A scaling diagnostic assembled from standard parts: k-nearest-neighbour spacing (<b>Clark &amp; Evans 1954</b>) inside each cluster, plotted against cluster size measured either as member count or as convex-hull area/volume (<b>quickhull</b>, Barber, Dobkin &amp; Huhdanpaa 1996)."),
            HTML("Spacing is computed on <b>de-blinked centroids</b> (trace medians, or de-blinked-molecule medians when de-blinking is on), so it measures molecule-to-molecule spacing inside a structure rather than the blink footprint — a measurement MINFLUX's nm precision actually supports. Both quantities are per-structure, so each cluster contributes one point regardless of how many localizations it collected."),
            HTML("<b>The slope is confounded by the clustering algorithm.</b> Under DBSCAN every within-cluster NN distance is bounded above by the link radius eps2 (default 200&nbsp;nm), so the curve is <i>censored</i> and cannot rise past it — a flat trend at high size may be the ceiling, not constant spacing. HDBSCAN/OPTICS impose no such bound; re-run under those before reading a slope. Reading the slope as 'ballooning' vs 'densifying' is this tool's interpretation, not a validated model, and convex-hull size is upward-biased by single outlying members.")),
          helpText("One point per Median Cluster structure across all files.",
                   "Y = median over members of each member's mean distance to its k nearest member",
                   "centroids (k = 1 is nearest-neighbour spacing; centroids are trace-ID medians,",
                   "or de-blinked-molecule medians if de-blinking is on).",
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
          h4("Cluster Viewer and Shape Fitting", .exploratory_badge()),
          .method_box("orange",
            HTML("Machinery is standard, dictionary is bespoke. Fits are <b>orthogonal-distance (total least squares) MLE</b> — minimizing the squared point-to-surface distance — ranked by <b>AIC</b> (Akaike 1974), its small-sample correction <b>AICc</b> (Hurvich &amp; Tsai 1989) or <b>BIC</b> (Schwarz 1978). The 2D view is a PCA projection onto PC1×PC2 (Pearson 1901). Under an iid isotropic-Gaussian residual model σ̂²&nbsp;=&nbsp;RMSD², so −2logL&nbsp;=&nbsp;n·log(2πσ̂²)+n and AICc&nbsp;=&nbsp;AIC&nbsp;+&nbsp;2k(k+1)/(n−k−1)."),
            HTML("Fits <b>de-blinked trace centroids</b>, not raw localizations — otherwise the many localizations of one long trace would dominate the residual sum and the fit would describe the brightest molecule rather than the structure. MINFLUX's nm precision is what makes a ~50–200&nbsp;nm shape distinguishable at all. AICc → ∞ when n&nbsp;≤&nbsp;k+1, so the ranking blanks out rather than over-fitting a 5-point cluster."),
            HTML("<b>The shape dictionary is imposed by this tool</b> (sphere / ellipsoid / helix / ring / arc) and 'best' is a model-selection outcome, not proof of biological geometry — the true shape may not be in the list, and the ranking will still name a winner. With few, blinking-inflated points several shapes are statistically indistinguishable; RMSD always favours the most flexible model (prefer AICc). The unequal-precision problem is <b>now addressed</b> by the optional photon (<code>eco</code>) weighting above — with it off, every trace is treated as equally precise, which MINFLUX traces are not.")),
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
            checkboxInput("viewer_weight_eco",
                          HTML("Weight fits by photon count (<code>eco</code>)"), value = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "MINFLUX localization variance scales as 1/photons, so a trace's centroid has variance ",
              "∝ 1/Σeco — its total photon count ", tags$b("is"), " its inverse-variance weight. ",
              "Unweighted (off) treats a 50-photon trace as being as trustworthy as a 5000-photon one. ",
              "The header reports ", tags$b("n_eff"), " (Kish effective sample size): if it falls far ",
              "below the trace count, a few very bright traces are carrying the fit. Needs an ",
              tags$code("eco"), " column; falls back to unweighted silently if absent. ",
              "Criterion values are not comparable between weighted and unweighted mode."),
            tags$div(style = "margin-top:2px; padding:5px 8px; border-left:3px solid #b26a00; background:rgba(178,106,0,0.06); font-size:0.8em; color:#7a5c00; line-height:1.45;",
              HTML(paste0(
                "<b>Filter double events first.</b> A high <code>eco</code> / <code>efo</code> is not only better precision,",
                "it is also the signature of <b>two fluorophores blinking in the same spot</b>, which is ",
                "implausible for a single dye. On unfiltered data, weighting therefore hands the most ",
                "influence to exactly the traces most likely to be double events, so the fit is pulled ",
                "toward artefacts. Apply proper <code>efo</code>/<code>cfr</code>/<code>eco</code> metadata filters in the ",
                "Display tab before enabling this. On a properly filtered dataset the weighted and ",
                "unweighted fits usually agree closely — a <b>large</b> difference between them is itself ",
                "a warning that unfiltered double events or a very uneven photon distribution are ",
                "driving the result."))),
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
          .method_box("green",
            HTML("Second moment of a point cloud about its centroid, Rg²&nbsp;=&nbsp;⟨‖p&nbsp;−&nbsp;p̄‖²⟩ — an RMS size. This tool answers <b>two</b> different questions depending on the <i>Measure</i> setting: the <b>size of your structures</b> (grouping by median cluster) or the <b>footprint of your labelling</b> (spread of one molecule's own repeat blinks). The gyration-tensor eigenvalues also give the shape descriptors used elsewhere (asphericity, <b>Rudnick &amp; Gaspari 1986</b>; linear/planar/spherical coefficients, <b>Westin et&nbsp;al. 1997/2002</b>)."),
            HTML("Points are <b>median positions</b>, never raw localizations, so localization noise is suppressed (a median of N localizations has SD ~ σ_loc/√N) and Rg reflects real geometry. For <b>structure size</b> with de-blinking on, each de-blinked molecule contributes <b>one</b> point, so a multivalent label's several traces do not inflate the size. For the <b>label footprint</b> measure, each of a molecule's traces is a point, so Rg is exactly the blink/label scatter (~linker length for a tag, ~10–40&nbsp;nm for an antibody)."),
            HTML("Rg is biased low at small member counts, because the centroid is estimated from the same points whose spread is measured. The bias grows as members are lost (29% at 2 members, 11% at 5, 1% at 50) and so is <b>member-count dependent</b>, which would create an artificial trend of Rg with cluster size. A Bessel <b>n/(n−1)</b> factor on the mean squared distance corrects it; at 2 members a residual underestimate remains, so treat those clusters as indicative. Rg is also a single number: a filament and a disc of the same extent share an Rg — that is what the anisotropy coefficients are for. The structure-size measure carries the label footprint as an additive floor (no footprint deconvolution is applied); compare against the label-footprint measure to gauge how much of a small structure's Rg is really labelling.")),
          helpText("Distribution of per-cluster Rg across all loaded files. Two views:",
                   "structure size (spread of members in the median cluster) or the per-molecule",
                   "labelling footprint."),
          radioButtons("rog_cluster_col", "Measure",
            choices = c(
              "Median-cluster size (structure)"                = "cluster2",
              "Label footprint — within-molecule blink spread" = "cluster"
            ),
            selected = "cluster2"
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
          .method_box("green",
            HTML("Cross-type nearest-neighbour distance distribution — for every point of channel A, the distance to its k-th nearest channel-B point (<b>Clark &amp; Evans 1954</b> extended to the bivariate case, cf. <b>Dixon 2002</b>). The second panel is the empirical cumulative distribution of those distances: the fraction of A molecules with a B partner within distance d — an object-based co-localization curve in the sense of <b>Lachmanovich et&nbsp;al. 2003</b> / <b>Dunn et&nbsp;al. 2011</b>, rather than a pixel-intensity coefficient."),
            HTML("Object-based co-localization is the <b>right</b> family for MINFLUX: there are no pixels and no intensity image, only coordinates, and nm precision makes the distance axis physically meaningful down to the label size. Cross-NN is computed strictly between channels (never within), on de-blinked centroids, and the Median-Cluster modes restrict pairing to molecules inside the same structure so a partner across the cell cannot be counted."),
            HTML("Cross-NN has <b>no null model built in</b> — every A point has a nearest B point, so the curve is never zero and a small median distance can arise from density alone. Judge it against Ripley&nbsp;L₁₂ or g₁₂ with their CSR envelopes. The dominant two-colour artefact here is <b>DCR crosstalk</b>: the split thresholds an overlapping DCR distribution, and every mis-assigned molecule contributes a near-zero distance — landing in the first histogram bin and the steep part of the co-localization curve. Read the co-localized fraction against the <b>per-group contamination</b> reported by the GMM auto-split, and use <b>pooled DCR</b> rather than raw DCR when available — it usually separates the dyes considerably better.")),
          helpText("Requires '2-color MINFLUX' to be enabled in the Display tab.",
                   "Channel assignment and clustering (tc1/tc2) are taken from the 2-color settings.",
                   "Median  modes restrict comparisons to molecules within the same Median Cluster group."),
          selectInput("cross_nn_unit", "Unit",
            choices = c(
              "TID centroids"                            = "tid",
              "De-blinked molecule centroids"                = "cluster",
              "TID centroids within Median Clusters"            = "cluster2_tid",
              "De-blinked molecule centroids within Median Clusters"= "cluster2_cluster"
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
          .method_box("green",
            HTML("Cross-type K-function for two point patterns (<b>Lotwick &amp; Silverman 1982</b>; <b>Dixon 2002</b>): K₁₂(r)&nbsp;=&nbsp;area·#{d&nbsp;≤&nbsp;r}/(n₁·n₂) over all A–B pairs, with the same Besag L-transform as the univariate case. L₁₂(r)&nbsp;−&nbsp;r&nbsp;&gt;&nbsp;0 = co-localization at scale r, &lt;&nbsp;0 = exclusion."),
            HTML("Runs on de-blinked centroids of the two channels (DCR/fluorophore split, or a chosen multiplex channel pair) with the ≤2000-point subsample per channel. The <b>CSR envelope randomizes both channels independently in the shared window</b>, which is the correct null for 'these two species are placed without regard to each other' and again cancels the uncorrected edge bias. Median-Cluster modes compute L₁₂ inside each structure and average the curves, so a co-localization signal cannot be produced by the two channels merely occupying the same region of the cell. <b>A MINFLUX advantage applies here:</b> the DCR split excites both dyes with the same laser through the same beam path, so both channels share one coordinate frame — there is <b>no chromatic offset and no channel registration to calibrate</b>, and L₁₂ is trustworthy down to the localization precision."),
            HTML("Same missing analytic edge correction as L(r) — read against the envelope, not the zero line. The real two-colour error is <b>DCR mis-assignment (crosstalk)</b>, not offset: the split thresholds a continuous, overlapping DCR distribution, and a wrongly-assigned molecule sits by construction right next to its true neighbours — so crosstalk biases L₁₂ <b>upward at short r, toward false co-localization</b>, exactly where a real interaction shows. <b>Quantify it before trusting a result:</b> the GMM auto-split in the Display tab now reports d′, the expected mis-assignment rate and the per-group contamination — a co-localization signal is only meaningful if it is large compared with those. <b>Pool DCR at upload</b> (photon-weighted across the iteration sequence) if you have not: it usually separates the two dyes far better than raw DCR, and a better split is worth more here than any statistical correction. (Drift can offset the multiplex-by-files / by-time-window modes, which overlay separate acquisitions — that does not apply to the DCR split.) Asymmetric labelling density also skews the estimator; and nm-scale co-localization is proximity, <b>not</b> evidence of interaction.")),
          helpText("Requires '2-color MINFLUX' to be enabled in the Display tab.",
                   "L₁₂(r) − r > 0 indicates co-localization; < 0 indicates spatial exclusion.",
                   "Median Cluster modes compute L₁₂ within each Median Cluster group and average the curves.",
                   "Large datasets are subsampled to 2000 points per channel per file/group."),
          selectInput("ripley12_unit", "Unit",
            choices = c(
              "TID centroids"             = "tid",
              "De-blinked molecule centroids" = "cluster"
            ),
            selected = "tid"),
          checkboxInput("ripley12_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          numericInput("ripley12_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          numericInput("ripley12_r_steps", "r steps", value = 60, min = 10, max = 300),
          checkboxInput("ripley12_smooth", "Smooth curve (loess)", value = TRUE),
          numericInput("ripley12_smooth_span", "Smoothing span", value = 0.2, min = 0.05, max = 0.8, step = 0.05),
          checkboxInput("ripley12_aggregate", "Aggregate all datasets into one", value = FALSE),
          checkboxInput("ripley12_envelope", "Show CSR envelope (Monte-Carlo)", value = FALSE),
          helpText("Randomizes both channels independently in the shared window and shades the",
                   "2.5–97.5% band. Observed L₁₂−r above the band = co-localization at that scale;",
                   "below = exclusion. Computed for the first few files only; not shown in aggregate mode."),
          conditionalPanel(
            condition = "input.ripley12_envelope == true",
            numericInput("ripley12_nsim", "Envelope simulations", value = 99, min = 19, max = 199)
          ),
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
          .method_box("green",
            HTML("The pair-correlation (radial distribution) function — the derivative form of Ripley's K, counting cross pairs in an <b>annulus/shell</b> (r−dr/2, r+dr/2] normalized by its area (2D) or volume (3D). g₁₂&nbsp;=&nbsp;1 under CSR, &gt;1 attraction, &lt;1 repulsion. In SMLM this is the PC-PALM statistic (<b>Sengupta et&nbsp;al. 2011</b>; <b>Veatch et&nbsp;al. 2012</b>), developed precisely to separate multiple appearances of one protein from genuine clustering."),
            HTML("Being <b>scale-local rather than cumulative</b>, g₁₂ is the more informative of the two cross statistics on MINFLUX data: it puts the interaction distance at a peak instead of smearing it into a rising curve, and a short-r blinking/labelling artefact stays confined to the first few bins instead of contaminating every larger r as it does in L₁₂. Same de-blinked units, shared window and independent two-channel CSR envelope as L₁₂."),
            HTML("<b>Noisier than L₁₂</b> — each shell holds few pairs, so short-r bins are unstable and sensitive to the r-step choice (an apparent peak in the first bin is often binning, not biology). Same missing analytic edge correction; read against the envelope, which brackets g&nbsp;=&nbsp;1. DCR crosstalk and asymmetric labelling affect g₁₂ exactly as they affect L₁₂ — and because g₁₂ is scale-local, crosstalk piles into the <i>first shells</i>, which is also where a genuine interaction peak would appear.")),
          helpText("g₁₂(r) = 1 under CSR; > 1 indicates attraction (co-localization) at scale r;",
                   "< 1 indicates repulsion. Complements Ripley L₁₂ — shows the",
                   "interaction scale as a peak rather than a cumulative curve.",
                   "Requires '2-color MINFLUX' to be enabled in the Display tab."),
          selectInput("pcf12_unit", "Unit",
            choices = c("TID centroids" = "tid", "De-blinked molecule centroids" = "cluster"),
            selected = "tid"),
          checkboxInput("pcf12_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          numericInput("pcf12_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          numericInput("pcf12_r_steps", "r steps", value = 60, min = 10, max = 300),
          checkboxInput("pcf12_aggregate", "Aggregate all datasets into one", value = FALSE),
          checkboxInput("pcf12_envelope", "Show CSR envelope (Monte-Carlo)", value = FALSE),
          helpText("Shades the 2.5–97.5% band of g₁₂(r) under independent randomization of both",
                   "channels in the shared window. g₁₂ above the band = attraction at that scale;",
                   "below = repulsion. Computed for the first few files only; not shown in aggregate mode."),
          conditionalPanel(
            condition = "input.pcf12_envelope == true",
            numericInput("pcf12_nsim", "Envelope simulations", value = 99, min = 19, max = 199)
          ),
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
          .method_box("amber",
            HTML("2D kernel density estimation (<b>Silverman 1986</b>; <code>MASS::kde2d</code>) of each channel on a shared grid, with the auto bandwidth from Silverman's rule bw&nbsp;=&nbsp;0.9·min(sd,&nbsp;IQR/1.34)·n<sup>−1/5</sup>. The co-localization index is the <b>Pearson correlation</b> of the two density surfaces (the coordinate-data analogue of Pearson co-localization on images, cf. <b>Dunn et&nbsp;al. 2011</b>), with an optional randomization test for significance."),
            HTML("MINFLUX has no pixels, so the density surface — and therefore the 'image' the Pearson r is computed on — is <b>manufactured by the bandwidth</b>, which is a free parameter rather than a property of the microscope. Two tools compensate: the <b>r-vs-bandwidth sweep</b> shows the whole family of answers instead of one, and the <b>randomization test</b> gives a p-value at the chosen bandwidth. Its default null is a <b>toroidal random shift</b> of channel 2 (<b>Lotwick &amp; Silverman 1982</b>), which preserves that channel's own clustering and point count and breaks only its registration with channel 1; a uniform-CSR null is selectable for comparison."),
            HTML("<b>The Pearson r is not an absolute quantity.</b> A small bandwidth drives r toward 0 even for perfectly co-localizing molecules (the two surfaces become disjoint spikes); a large one inflates r toward 1 because both channels simply fill the cell. Only compare r values computed at the same bandwidth on similarly dense data, and prefer the sweep to any single number. <b>The p-value is calibrated only when the signal fills the field:</b> both nulls randomize over the bounding box rather than the region the signal occupies, so for a structure confined to a small part of the view they reject too readily — the uniform-CSR null severely so, since p&nbsp;≈&nbsp;0 then mostly reports 'both channels are non-uniform'. For a bandwidth-free, scale-resolved answer use g₁₂(r) above.")),
          helpText("Kernel density estimate for each channel on the active file's spatial grid.",
                   "Ch1 shown as red contours, Ch2 as blue contours. Pearson r of the two",
                   "density surfaces serves as a global co-localization coefficient,",
                   "but is strongly bandwidth-dependent — a small bandwidth gives near-zero r",
                   "even for perfectly co-localizing molecules, while a large bandwidth inflates r.",
                   "Treat the Pearson r as relative rather than absolute.",
                   "Select the dataset to inspect via 'Active file' in the Display tab.",
                   "Requires '2-color MINFLUX' to be enabled in the Display tab."),
          selectInput("kde_unit", "Unit",
            choices = c("TID centroids" = "tid", "De-blinked molecule centroids" = "cluster"),
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
          checkboxInput("kde_bw_sweep", "Show r vs bandwidth sweep", value = FALSE),
          helpText("Recomputes Pearson r across a range of bandwidths to expose how strongly the",
                   "single-number coefficient depends on the (arbitrary) smoothing choice."),
          numericInput("kde_costes_nsim", "Randomizations (0 = off)", value = 100, min = 0, max = 500, step = 10),
          helpText("Randomizes Ch2 this many times and recomputes r, yielding a p-value",
                   "(fraction of random r ≥ observed). Low p = co-localization beyond chance."),
          radioButtons("kde_null_mode", "Randomization null",
            choices = c("Random shift of Ch2 (recommended)" = "shift",
                        "Uniform CSR resampling of Ch2"     = "csr"),
            selected = "shift"),
          helpText(HTML("<b>Random shift</b> translates Ch2 with wrap-around, preserving its own clustering and point count, so only the cross-channel registration is broken (Lotwick &amp; Silverman 1982). <b>Uniform CSR</b> replaces Ch2 with uniform random points, testing against complete spatial randomness instead.<br><b>Limitation:</b> both nulls randomize over the bounding box rather than the region the signal actually occupies. When the signal fills the field the p-value is well calibrated; when it is confined to a small part of the field (a single cell in a wider view) both nulls reject far too readily, CSR much more so. In that case treat a low p as weak evidence and confirm with g₁₂(r) against its envelope.")),
          actionButton("run_kde_coloc", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_kde_coloc", "Download CSV", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("kde_coloc_plot", height = "450px"),
          uiOutput("kde_pearson_display"),
          conditionalPanel(
            condition = "input.kde_bw_sweep == true",
            plotlyOutput("kde_bw_sweep_plot", height = "280px")
          )
        )
      ),

      tags$br(),

      # ── Methods & References (consolidated citation box) ──────────────────────
      # Provenance for the SMLM spatial-analysis tools. Panels flagged EXPLORATORY
      # (orange) are bespoke to this tool and not literature-validated for MINFLUX.
      tags$details(
        style = "margin:8px 0 16px; font-size:0.85em; color:#444;",
        tags$summary(style = "cursor:pointer; font-weight:bold; color:#2c3e50;",
                     "Methods & References"),
        tags$div(style = "padding:8px 4px 4px 12px;",
          tags$p(style = "color:#666; margin:0 0 6px;",
            "Primary reference for each spatial-analysis tool. MINFLUX is recent and has little",
            "method development of its own, so most tools are well-established SMLM / point-pattern",
            "methods applied to MINFLUX localizations. Panels marked ",
            tags$span("EXPLORATORY", style = "color:#e67e22; font-weight:700;"),
            " are bespoke to this tool and are not validated against the literature."),
          tags$table(style = "border-collapse:collapse; line-height:1.5;",
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Nearest-neighbour")),
                    tags$td("k-th NN distance distribution via k-d tree; NN distance as a spatial statistic after Clark & Evans 1954 — their observed/expected R index is not computed (no null model here)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("DBSCAN / HDBSCAN / OPTICS")),
                    tags$td("Ester, Kriegel, Sander & Xu 1996 (DBSCAN); Ankerst, Breunig, Kriegel & Sander 1999 (OPTICS); Campello, Moulavi & Sander 2013 (HDBSCAN)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Median (de-blinked) clustering")),
                    tags$td("Two-tier grouping (localizations → trace → de-blinked molecule → median cluster) collapses multi-blink / multivalent overcounting before clustering — an adaptation for MINFLUX/SMLM in the spirit of Annibale et al. 2011; Sengupta et al. 2011 (PC-PALM)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Trace IDs / de-blinked molecules per median cluster")),
                    tags$td("Median-cluster occupancy counting in the SMLM molecular-counting lineage (Annibale et al. 2011; Sengupta et al. 2011) — trace IDs (blinking events) or de-blinked molecules, not stoichiometry")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Cross-channel NN & co-loc curve")),
                    tags$td("Bivariate NN distances (Clark & Evans 1954; Dixon 2002) with an object-based co-localization read-out (Lachmanovich et al. 2003; Dunn et al. 2011)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Radius of gyration")),
                    tags$td("Standard second-moment cluster size (gyration tensor); anisotropy after Rudnick & Gaspari 1986")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Ripley's L(r)")),
                    tags$td("Ripley 1977; variance-stabilising L-transform, Besag 1977")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Bivariate Ripley L₁₂(r)")),
                    tags$td("Lotwick & Silverman 1982; Dixon 2002 (bivariate/cross point-pattern statistics)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Pair cross-correlation g₁₂(r)")),
                    tags$td("Sengupta et al. 2011 (PC-PALM); Veatch et al. 2012 (correlation functions / over-counting)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("KDE density + Pearson co-loc")),
                    tags$td("Kernel density estimate (Silverman 1986); Pearson coefficient as a bandwidth-dependent, relative co-localization measure (cf. Dunn et al. 2011); significance by structure-preserving random-shift randomization of Ch2 (Lotwick & Silverman 1982), with a uniform-CSR null selectable for comparison")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Shape fitting (viewer)")),
                    tags$td(tags$span("EXPLORATORY", style="color:#e67e22;font-weight:700;"),
                            " — orthogonal-distance MLE with AICc/BIC model selection (Akaike 1974; Schwarz 1978; Hurvich & Tsai 1989); the shape dictionary and biological reading are bespoke")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Ground-truth estimate & structure size")),
                    tags$td(tags$span("EXPLORATORY", style="color:#e67e22;font-weight:700;"),
                            " — labelling-model inverse (this tool), validated in silico against the simulator; contour via MST (Prim 1957), anisotropy coefficients (Westin et al. 1997/2002)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Detection-time / crowding / NN-vs-size")),
                    tags$td(tags$span("EXPLORATORY", style="color:#e67e22;font-weight:700;"),
                            " — descriptive diagnostics of this tool with no published MINFLUX equivalent; correlational, not causal"))
          ),
          tags$p(style = "color:#888; margin:8px 0 0; font-size:0.95em;",
            tags$b("MINFLUX caveat for all spatial statistics: "),
            "multi-blink and multivalent-antibody labelling scatter several localizations around one",
            "molecule, which inflates short-range clustering, NN peaks and pair-correlation signal as a",
            "label artefact rather than biology (Annibale et al. 2011; Sengupta et al. 2011). Working from",
            "median (de-blinked) clusters mitigates this; results on raw localizations should be read with",
            "that bias in mind. MINFLUX foundational: Balzarotti et al. 2017; Gwosch et al. 2020; ",
            "Ostersehlt et al. 2022 (DNA-PAINT MINFLUX).")
        )
      )
    )  # end tabPanel("Analysis")

  )  # end tabsetPanel
)  # end fluidPage
