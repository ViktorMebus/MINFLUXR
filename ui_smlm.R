## -----------------------------------------------------------------------------
## ui_smlm.R  -  SMLM analysis tab UI (assigned to smlm_ui)
## -----------------------------------------------------------------------------

# Small orange "EXPLORATORY" chip reused by bespoke / not-literature-validated
# analysis panels (same visual language as the GT-estimation panel).
.exploratory_badge <- function() {
  tags$span("EXPLORATORY",
            style = "margin-left:8px;font-size:10px;font-weight:700;color:#fff;background:#e67e22;border-radius:3px;padding:1px 5px;vertical-align:middle;")
}

# NOTE: the former .exploratory_note() footnote helper was removed - the per-panel
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

  /* -- helpers -- */
  function getPlotDiv() {
    return document.getElementById('scatter3D');
  }
  function currentCamera() {
    var d = getPlotDiv();
    if (!d || !d.layout || !d.layout.scene) return null;
    return d.layout.scene.camera || null;
  }

  /* -- rotation -- */
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

  /* -- recording -- */
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

    # -- Display ---------------------------------------------------------------
    tabPanel("Display",
      # Sticky Plot button bar - always visible regardless of sidebar scroll
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
        downloadButton("download_filtered_csv", "Export MINFLUX data (.csv)"),
        helpText("Exports data with current filters applied."),
        uiOutput("z_axes_note"),
      ),  # /Data Upload
      tags$hr(),

      
      ### Axial calibration ===
      # Deliberately OUTSIDE tags$div(id = "filters_form"): that div is wiped by
      # shinyjs::reset("filters_form") whenever a dataset changes, which would
      # snap the calibration back to 0.7 - and since applying it updates the
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

                  # -- Outlier localization trimming --------------------------
                 checkboxInput("trim_outliers",
                               "Trim outlier localizations within traces", value = FALSE),
                 conditionalPanel(
                   condition = "input.trim_outliers == true",
                   numericInput("trim_k", "Reject beyond k × σ of the trace",
                                value = 4, min = 2, max = 10, step = 0.5),
                   uiOutput("trim_summary_ui")
                 ),
                 helpText("Flags single localizations sitting far from the rest of their own trace ",
                          "and excludes them from the spread, the trace length and every downstream ",
                          "analysis. The distance is measured per axis in units of that trace's own ",
                          "spread."),
                 
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
                          "the median-spread labels, the allowable-spread ranges, and the σ filters. "),
                 uiOutput("sd_summary_ui"),
                 uiOutput("sd_hist_ui"),
                 sliderInput("sd_bins", "Histogram bins", min = 5, max = 500, value = 100, step = 1),
                 checkboxInput("sd_hist_shared",
                               "Shared axes across σ X / Y / Z", value = FALSE),
                 helpText("Puts the three histograms on one set of bin edges and one count ",
                          "axis, so they can be read against each other and used as a figure ",
                          "set. Off by default because axial precision differs from lateral, ",
                          "so per-axis ranges show each panel at full width."),
                 conditionalPanel(
                   condition = "input.sd_hist_shared == true",
                   shinyWidgets::numericRangeInput(
                     "sd_hist_range", "Shared σ axis range (nm)", value = c(0, 30)),
                   helpText("Display range only, clear both boxes ",
                            "to fall back to the full range of the data."),
                   uiOutput("sd_hist_range_note")
                 ),
                 checkboxInput("sd_overlay_all", "Show raw trace spreads (grey bars)", value = FALSE),
                 checkboxInput("sd_show_mean_sd", "Show median trace spread (post filtering)", value = FALSE),
                 checkboxInput("sd_show_overall_mean_sd", "Show median trace spread (raw data)", value = FALSE),
                 tags$span(
                   " (What is \u03c3\u0303 \u00b1 nMAD?)",
                   title = paste(
                     "\u03c3\u0303 = per-trace localization spread: the standard deviation of x/y/z within each trace, normalized by n\u22121.",
                     "Reported per axis as median(\u03c3\u0303) \u00b1 nMAD, with the interquartile range of \u03c3\u0303 alongside -",
                     "nMAD = 1.4826 \u00d7 median(|\u03c3\u0303_i \u2212 median(\u03c3\u0303)|).",
                     "Factor 1.4826 \u2248 1/\u03a6\u207b\u00b9(3/4): makes nMAD a consistent estimator of \u03c3 for normal data.",
                     "nMAD and IQR describe how \u03c3\u0303 varies BETWEEN traces."
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
                                 "Spread > median + k \u00d7 MAD"  = "mad_rule",
                                 "Any axis spread > k \u00d7 median(axis spread)" = "axis_sd"),
                     selected = "factor_median"
                   ),
                   numericInput("rogue_factor", "Multiplier k", value = 2.5, min = 1, step = 0.1),
                   helpText("k applies to whichever rule is selected above.")
                 )
               ),
               actionButton("reset_filters", "Reset Filters", icon = icon("refresh"))



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
            helpText("Counts how many other trace-median positions fall within this radius. Colors traces by that local density.")
          ),
          checkboxInput("numeric_palette_reverse", "Reverse colorscale", value = FALSE)),
        
        checkboxInput("show_medians", "Show trace median positions.", value = FALSE),
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
          helpText("Wraps each trace (or each de-blinked ID) in a 3D alpha ",
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
              "trim_mode", "Crop shape",
              choices = c("Axis-aligned box" = "box",
                          "Sphere"           = "sphere"),
              selected = "box"
            ),
            conditionalPanel(condition = "input.trim_mode == 'box'",   uiOutput("trim_box_ui")),
            conditionalPanel(condition = "input.trim_mode == 'sphere'", uiOutput("trim_sphere_ui")),
            checkboxInput("trim_preview", "Preview region in 3D", value = TRUE),
            helpText(style = "font-size:0.78em; color:#888;",
                     "The coordinates draw a low-opacity shape for you to check against the points.",
                     "Nothing is cropped until you press Apply."),
            fluidRow(
              column(7, actionButton("trim_apply", "Apply crop",
                                     class = "btn-primary btn-sm", style = "width:100%;")),
              column(5, actionButton("trim_clear", "Clear",
                                     class = "btn-default btn-sm", style = "width:100%;"))
            )
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
               "de-blinked ID, so downstream analyses don't overcount on the same molecule."),
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
          numericInput("centroid_k", "Threshold in units of trace spread (k)", value = 2.8,
                       min = 0.5, max = 6, step = 0.1),
          helpText("Changing k rewrites the distance above. Too low leaves re-blinks counted ",
                   "twice; too high fuses neighbouring molecules. Grouping accuracy is flat ",
                   "between about 2.6 and 3.1, and because a molecule can never be wider than ",
                   "the threshold, erring high costs less than erring low."),
          helpText("Every trace in a molecule is within the threshold of every other, so a ",
                   "molecule can never be wider than the threshold above."),
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
                       value = 0.3, min = 0.02, max = 0.9, step = 0.05),
          helpText("Two traces of one fluorophore score about 0.5, so useful values are well ",
                   "below that. Scores fall to ~0.23 at 10 nm apart and ~0.07 at 20 nm. ",
                   "Below ~0.25 a single localization can decide a short trace."),
          numericInput("vote_rho", "Distance tolerance (rho)", value = 3.0, min = 0.5, max = 8, step = 0.1),
          helpText("1.0 = the neighbouring trace must be strictly closer; higher is more willing to merge. ",
                   "Grouping accuracy rises with rho and then flattens, so erring high costs little."),
          radioButtons("vote_link", "Agreement", inline = TRUE,
                       choices = c("Both traces must agree" = "mutual", "Either is enough" = "either"),
                       selected = "mutual")
        ),

        # -- CHANNELS ---------------------------------------------------------
        # Only appears once the data actually has more than one channel (a
        # 2-colour split or a multiplex build). There is no "pool the channels"
        # option: merging across a channel is what the channel exists to prevent.
        conditionalPanel(
          condition = "output.n_channels > 1",
          tags$div(style = "margin:8px 0; padding:7px 10px; border-left:4px solid #cfe0ee; background:#f6f9fc;",
            uiOutput("channel_scope_ui"),
            helpText("Blinks are only ever merged inside one channel, so no de-blinked ",
                     "molecule can contain two colours or two multiplex rounds."),
            checkboxInput("deblink_shared_spread",
                          "Use one combined trace spread for all channels", value = FALSE),
            helpText("Off (default): each channel is merged at its own measured spread, ",
                     "because localization precision is a property of the dye. On: one ",
                     "pooled spread, which is faster on very large files but under-serves ",
                     "the tighter channel."),
            uiOutput("channel_spread_summary"),
            uiOutput("channel_integrity_summary")
          )
        ),

        # -- anisotropic metric (see deblink_axis_weights in helpers_smlm.R) --
        checkboxInput("deblink_anisotropic",
                      "Scale axes by their own spread (anisotropic radius)", value = FALSE),
        helpText("Weights x, y and z by their measured trace spread before merging, so one radius ",
                 "fits all three axes. Affects the merge decision only - reported distances, ",
                 "cluster sizes and Rg stay in true nm."),


        helpText("TIDs are not split: all localizations of a trace ID join the same de-blinked molecule (majority vote). ",
                 "Trace IDs that overlap spatially are merged into one de-blinked molecule. ",
                 "De-blinking only merges - a trace it cannot merge becomes its own single-blink molecule, ",
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
        # NN scope is THREE-valued, and unlike clustering "between" is legitimate
        # here: NN is a distance measure, not a merge. Comparing within against
        # between is one of the questions this panel exists to answer. One control
        # covers both a 2-colour split and a multiplex build - they used to be two
        # separate controls doing the same job.
        conditionalPanel(
          condition = "output.n_channels > 1",
          radioButtons("nn_mplex_scope", "Channel scope",
                       choices = c("Within channel"   = "within",
                                   "Between channels" = "across",
                                   "All pairs"        = "all"),
                       selected = "within", inline = TRUE),
          helpText("Within: the nearest neighbour is searched only inside each channel. ",
                   "Between: only in OTHER channels (never a same-channel pair). ",
                   "All: both, which is the distribution the other two partition.")
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
          numericInput("eps2", "DBSCAN Epsilon (nm)", value = 200, min = 3)
        ),
        conditionalPanel(
          condition = "input.cluster2_method == 'optics'",
          numericInput("cluster2_xi", "OPTICS Xi (steepness, 0-1; lower = more/smaller clusters)", value = 0.05, min = 0.001, max = 0.5, step = 0.005)
        ),
        numericInput("minpts2", "MinPts (all methods)", value = 15, min = 1),
        helpText("Tip: set MinPts near the smallest number of TIDs you'd call a real structure."),
        # -- local-density prefilter (drop scattered background TIDs) --
        checkboxInput("cluster2_prefilter", "Density prefilter (drop sparse background TIDs first)", value = FALSE),
        conditionalPanel(
          condition = "input.cluster2_prefilter == true",
          numericInput("cluster2_prefilter_r", "Prefilter radius (nm)", value = 200, min = 1),
          numericInput("cluster2_prefilter_min", "Min neighbours within radius", value = 3, min = 1),
          helpText("Keeps only centroids with at least this many neighbours within the radius, then groups the survivors.")
        ),
        # -- core-only membership (DBSCAN*) --
        conditionalPanel(
          condition = "input.cluster2_method == 'dbscan'",
          checkboxInput("cluster2_core_only", "Core points only (drop border points)", value = FALSE),
          helpText("DBSCAN also assigns border points: centroids within eps of a core point but sparse themselves.",
                   "This keeps only centroids that reach MinPts neighbours within eps.",
                   "HDBSCAN and OPTICS are already defined this way.")
        ),
        # -- minimum members per cluster --
        numericInput("cluster2_min_members", "Min members per Median Cluster", value = 0, min = 0, step = 1),
        helpText("Clusters left with fewer members than this are dropped whole. 0 = off. Counted after any core-only or pruning step."),
        # -- sparse-member pruning (post-clustering trim of loose members) --
        checkboxInput("cluster2_prune", "Prune sparse cluster members (trim loosely-attached traces)", value = FALSE),
        conditionalPanel(
          condition = "input.cluster2_prune == true",
          numericInput("cluster2_prune_tol", "Prune tolerance (x local core distance; lower = stricter)", value = 2, min = 1.1, step = 0.1),
          helpText("Demotes to noise any member whose core distance (radius needed to reach MinPts neighbours)",
                   "exceeds this multiple of the median core distance of the members around it.",
                   "The baseline is local, so a density gradient within a structure does not condemn its sparser side.")
        ),
        # -- kNN-distance diagnostic to pick eps --
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
          # Shared by the alpha shape and the isosurface. Only meaningful with more
          # than one channel, so it only appears then.
          conditionalPanel(
            condition = "output.n_channels > 1 && (input.show_cluster2_shape == true || input.show_cluster2_isosurface == true)",
            checkboxInput("cluster2_iso_by_channel", "Color wraps by channel", value = TRUE),
            helpText("Each Median Cluster is wrapped within its own channel, so a wrap never",
                     "joins two colours. This colours each wrap by that channel.")
          ),
          conditionalPanel(
            condition = "input.show_cluster2_isosurface == true",

            # Live readout of the geometry the surface was actually built with.
            # It sits directly above the shape controls because those adjust the
            # numbers it reports.
            uiOutput("cluster2_iso_status"),

            tags$b("Which clusters"),
            checkboxInput("cluster2_iso_selected_only", "Only selected Median Cluster ID", value = FALSE),
            conditionalPanel(
              condition = "input.cluster2_iso_selected_only == true",
              numericInput("cluster2_iso_id", "Selected Median Cluster ID", value = 1, min = 1, step = 1),
              helpText("IDs are unique across channels, so one ID is one structure in one colour.")
            ),

            tags$hr(style = "margin:10px 0;"),
            tags$b("Shape"),
            helpText(tags$em("Applied when you press Plot Data.")),
            sliderInput("cluster2_iso_pad", "Tighten or loosen the wrap (nm)", min = -40, max = 200, value = 0, step = 5),
            helpText("Shifts the wrap distance reported above. Negative tightens, positive loosens."),
            sliderInput("cluster2_iso_smooth_rel", "Smoothing", min = 0.1, max = 2, value = 0.4, step = 0.05),
            helpText("Scales the automatic smoothing. Lower keeps more of the measured detail,",
                     "higher gives a simpler surface that absorbs small protrusions."),
            sliderInput("cluster2_iso_grid_res", "Grid resolution (higher = slower)", min = 12, max = 40, value = 22, step = 1),
            helpText("Voxel budget for the surface. The grid is refined automatically when the",
                     "wrap distance is finer than the voxel."),

            tags$hr(style = "margin:10px 0;"),
            tags$b("Appearance"),
            sliderInput("cluster2_iso_opacity", "Isosurface opacity", min = 0, max = 1, value = 0.25, step = 0.01),
            checkboxInput("cluster2_iso_unicolor", "Single color for all clusters", value = FALSE),
            conditionalPanel(
              condition = "input.cluster2_iso_unicolor == true",
              colourpicker::colourInput(
                "cluster2_iso_color",
                "Isosurface color",
                value = "#00CED1"
              )
            ),

            tags$details(
              style = "margin-top:10px;",
              tags$summary(tags$b("Override the automatic geometry")),
              helpText("The wrap sits at 3 x the measured localization precision, held to at least",
                       "half the spacing between the molecules being wrapped. Smoothing sits at",
                       "0.85 x that spacing. Both are set per axis, so a structure measured more",
                       "precisely along z is wrapped more tightly along z."),
              helpText(tags$em("Applied when you press Plot Data.")),
              checkboxInput("cluster2_iso_manual_offset", "Set the wrap distance manually", value = FALSE),
              conditionalPanel(
                condition = "input.cluster2_iso_manual_offset == true",
                sliderInput("cluster2_iso_offset", "Wrap distance (nm)", min = 2, max = 200, value = 30, step = 2),
                helpText("Measured laterally. The axial wrap follows the axial precision.")
              ),
              checkboxInput("cluster2_iso_manual_sigma", "Set smoothing manually", value = FALSE),
              conditionalPanel(
                condition = "input.cluster2_iso_manual_sigma == true",
                sliderInput("cluster2_iso_sigma", "Smoothing (sigma, nm)", min = 2, max = 200, value = 30, step = 2),
                helpText("Absolute bandwidth, replacing the automatic value and the Smoothing dial.",
                         "The wrap is re-measured at whatever is set here.")
              )
            )
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
                "Regularized: MAP prior toward optimized centre"  = "reg",
                "Fixed to optimized centre (analytical radius)"    = "fixed",
                "Free: bounded to data bounding box"           = "free"
              ),
              selected = "reg"
            ),
            helpText("Fixed: centre = convex-hull centroid; radius solved analytically (no iteration). Regularized: centre is free but penalized for moving away from the optimized centre (Gaussian MAP prior, strength \u03bb). Free: centre fully optimized within an expanded data bounding box."),
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
            helpText("Used only to seed the pitch estimate. The fitted turn count is determined by the data extent and the optimized pitch - it will not be exactly this value.")
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
                   tags$p("Pixel size must be set first. Results are a starting point - nudge to refine.",
                          style = "font-size:11px; color:#888; margin:0 0 6px 0;"),
                   checkboxInput("img_invert_y", "Image origin at top-left (invert Y)", TRUE),
                   checkboxInput("img_flip_x", "Mirror horizontally (left↔right)", TRUE),
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
                   selectInput("img_rotation", "Rotate image", choices = c("0","90","180","270"), selected = "180"),
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

            ### Channels: 2-colour split ===
      # This block DEFINES channels. It sets no clustering parameters - those
      # live once, in De-blink correction and Median Clustering above, and apply
      # to every channel. The other source of channels is the Multiplexed MINFLUX
      # block below.
      tags$details(class = "sim-group",
      tags$summary(tags$b("Channels: 2-colour split")),
      checkboxInput("enable_two_color", "Split into 2 channels by DCR (per-trace median)", value = FALSE),
      helpText("Splits traces into two channels on their median DCR (or on 'fluo'). ",
               "Nothing is ever grouped across a channel: de-blinked molecules and ",
               "Median Clusters each belong to exactly one."),
      helpText("Note: numeric colour modes (z/time/density) override channel colours. ",
               "Use 'Default' colour mode to see the channels."),
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
              " (posterior below the cut) - these are the molecules that would change channel ",
              "under a small shift of the threshold, and the ones that manufacture false ",
              "short-range cross-channel pairs. Use the export to exclude or inspect them."),
            plotlyOutput("dcr_conf_plot", height = "260px"),
            downloadButton("download_dcr_conf", "Download per-trace confidence", class = "btn-sm"))
        ),
        
        colourpicker::colourInput("color_group1", "Color: Group 1", value = "#1f77b4"),
        colourpicker::colourInput("color_group2", "Color: Group 2", value = "#d62728"),
        
        # NN is NOT duplicated here either. The Nearest Neighbor Analysis panel
        # above ADAPTS when channels exist: it gains a Within / Between / All
        # channel scope, and each line is coloured by the channels it joins.
        # Two NN sections meant two sets of k, distance range and label controls
        # that could disagree about what a neighbour is.
        # De-blink and Median Clustering parameters are NOT duplicated here any
        # more. There is one clustering pipeline, so there is one set of controls
        # (De-blink correction / Median Clustering, above), and enabling 2-colour
        # only changes WHICH points may be grouped together - never how.
        conditionalPanel(
          condition = "input.enable_clustering == true",
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
        conditionalPanel(
          condition = "input.enable_cluster2 == true",
          helpText("Median Clusters belong to one channel: a cluster can never hold ",
                   "molecules from two colours. The eps and MinPts above apply to every ",
                   "channel - they are set by the structure, which the choice of dye ",
                   "does not change."),
          helpText("The Median Cluster wrap and isosurface are set under Median Clustering. ",
                   "With two colours they are drawn separately for each channel, so one ",
                   "surface never encloses both.")
        ),
        uiOutput("two_color_settings_file_ui"),
        downloadButton("download_two_color_settings", "Save channel settings"),
        actionButton("reset_two_color_settings", "Reset channel settings", icon = icon("undo"))
      ),
      ),  # /Channels: 2-colour split
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
        helpText("Exploratory, model-based estimate. Reverse-engineers the labeling model to turn traces / median-clusters back into candidate ground-truth positions, an approximate count, and a crowding flag. Hollow circles are coloured green (low) → red (high) by position uncertainty. See the method/reference notes below the read-outs."),
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
          tags$details(style = "margin:6px 0 10px; border:1px solid #dfe6ec; border-radius:4px; background:#fbfcfd;",
            tags$summary(style = "cursor:pointer; padding:6px 10px; font-weight:600; color:#2c3e50;",
              "Calibrate multiplicity from a standard",
              tags$span(style = "font-weight:400; color:#7b8794; margin-left:6px;",
                        "measured TIDs per molecule")),
            tags$div(style = "padding:4px 10px 8px;",
              helpText(HTML("The count divides trace IDs by the expected TIDs per molecule, which by default comes from the preset. With a <b>counting standard</b> (a sparse, isolated sample of known stoichiometry imaged with the same label and acquisition) you can measure that number directly. Each Median Cluster in the standard is one structure.")),
              selectInput("gt_calib_file", "Standard dataset", choices = NULL),
              numericInput("gt_calib_k", "True molecules per structure", value = 1, min = 1, step = 1),
              helpText("1 for isolated single molecules, such as sparse DNA docking strands with one site each."),
              actionButton("gt_calib_run", "Measure multiplicity", icon = icon("ruler"), class = "btn-sm"),
              uiOutput("gt_calib_readout"),
              checkboxInput("gt_use_calib", "Use the measured multiplicity for the count", value = FALSE))),
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
          helpText("Traces / de-blinked molecules within this radius are pooled into ONE candidate GT (a multivalent label scatters many fluorophores around one molecule). Blank = auto, set per labeling preset (max of footprint and 2.5 × dye spread) - big for polyclonal+secondary, tiny for a nanobody; 0 = no merge (one unit = one GT)."),
          checkboxInput("gt_guided_merge",
                        "Stoichiometry-guided merge (model-driven)", value = FALSE),
          helpText("MODEL-DRIVEN: where the spatial merge leaves more GTs than the labeling stoichiometry predicts (TIDs ÷ expected TIDs-per-molecule), collapse further to that target. Only ever reduces the count toward stoichiometry, never splits. Use when the spatial and stoichiometric counts disagree - but treat the result as prior-driven, and always sanity-check against the two-count readout below."),
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
                         "Assumed subunit spacing (nm, optional - blank = off)",
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

    # -- Analysis --------------------------------------------------------------
    tabPanel("Analysis",
      tags$br(),

      # How filtering relates to the Display tab. Worth stating once, up here: the
      # batch analyses read the SAME filter inputs as the Display tab and apply them
      # identically to every loaded file, but the sliders' BOUNDS are auto-detected
      # from the active file, which is the part that surprises people.
      tags$div(style = "margin:0 0 12px; padding:7px 12px; border-left:4px solid #cfe0ee; background:#f6f9fc; font-size:0.86em; color:#445;",
        tags$b("Filtering: "),
        "every analysis on this tab applies the ", tags$b("Display-tab filter ranges"),
        " (EFO, CFR, DCR, ECO, dwell, fluo, trace-ID length, trace spread, outlier trimming) ",
        "identically to every loaded file, including the ones you are not viewing. Note the sliders' ",
        tags$i("limits"), " are auto-detected from the active file, so if another dataset runs ",
        "higher on a metric, check the range still means what you intend before comparing groups. ",
        "The time (", tags$code("tim"), ") range is per-file by nature, since each acquisition has its own timeline."),

      # -- Density and area normalization (opt-in) --------------------------------------
      # Tucked into a closed <details> because these only matter for a specific
      # situation - comparing cells, where each cell got an
      # undefined amount of label load so labelling density is not a property of the
      # condition. 
      fluidRow(column(12,
        tags$details(style = "margin:0 0 14px; border:1px solid #dfe6ec; border-radius:4px; background:#fbfcfd;",
          tags$summary(style = "cursor:pointer; padding:7px 12px; font-weight:600; color:#2c3e50;",
                       "⚙ Density & area normalization",
                       tags$span(style = "font-weight:400; color:#7b8794; margin-left:8px;",
                                 "for samples with variable labelling density or images with empty regions")),
          tags$div(style = "padding:4px 12px 10px;",
            tags$div(style = "font-size:0.85em; color:#7a5c00; margin-bottom:8px;",
              tags$b("Note: "),
              "trace-ID and de-blinked-molecule counts depend on acquisition length, labelling ",
              "efficiency, expression level and buffer, so two acquisitions can differ in density ",
              "without differing in biology.",
              tags$br(), tags$br(),
              tags$b("Empty regions manufacture clustering. "),
              "Ripley's L, g(r) and the co-localization null all divide by the observed area and ",
              "assume a point could have landed anywhere in it. Counting empty field as sample ",
              "puts the expected density too low, so every real point appears to have more ",
              "neighbours than chance allows. K is inflated by the ratio of counted area to true ",
              "sample area, which turns L(r) − r into a straight rising ramp, the same signature ",
              "as real clustering but produced by framing alone. The ramp grows with that ratio ",
              "(as its square root in 2-D, its cube root in 3-D), so a cell filling half the frame ",
              "is hit far harder than one that fills it. The footprint window fixes this.",
              tags$br(), tags$br(),
              tags$b("Density moves the distance-based measures on its own. "),
              "Under complete randomness the mean nearest-neighbour distance falls as the inverse ",
              "square root of density in 2-D and the inverse cube root in 3-D. Doubling the label ",
              "density shortens it by about 29% and 21% respectively, with no change in ",
              "organisation. DBSCAN is more direct still: a fixed eps with minPts is itself a ",
              "density threshold, so a denser field crosses it more easily and the same random ",
              "pattern yields larger and fewer, then finally one merged, cluster. ",
              tags$b("K(r) and g(r) divide density out by construction"), ", so the correction ",
              "here is to thin the denser file down to a common density instead of rescaling a curve."),
          # -- Support-aware observation window (one switch, several panels) --------
          fluidRow(column(12,
            tags$div(style = "margin:0 0 12px; padding:7px 12px; border-left:4px solid #cfe0ee; background:#f6f9fc;",
              checkboxInput("spatial_edge_correct",
                            "Apply boundary (edge) correction (recommended)",
                            value = TRUE),
              checkboxInput("spatial_use_support",
                            "Also replace the bounding box with the sample footprint (if your cell(s) do not fill the frame)",
                            value = FALSE),
              tags$div(style = "font-size:0.85em; color:#445; line-height:1.5;",
                "Both apply to ", tags$b("Ripley's L, g(r), L\u2081\u2082, g\u2081\u2082 and the co-localization null"),
                ". They fix different things:",
                tags$br(), tags$br(),
                tags$b("Edge correction"), " A point near the edge of the ",
                "window has fewer neighbours available to it than one in the middle, simply because ",
                "part of its surroundings falls outside the imaged area. Uncorrected, every curve drifts as r grows ",
                "and a random pattern stops reading as random.",
                tags$br(), tags$br(),
                tags$b("Footprint window"), " changes what counts as the observed area. If your cell ",
                "fills the scan field the bounding box already is the cell, so leave this off. Turn ",
                "it on when the cell is concave or only partly in frame, because then the box ",
                "contains a large empty region and the pattern looks clustered purely from framing.",
                tags$br(), tags$br(),
                tags$span(style = "color:#7a5c00;",
                  "Two things to know. Correcting the edge of a bounding box that contains a lot of ",
                  "empty space is worse than leaving it uncorrected, so on a concave cell turn on ",
                  "both settings together. This window traces the ",
                  tags$i("outline of the cell"), ". A window drawn tightly around your structures ",
                  "would erase the clustering it is meant to measure, so if the footprint comes out ",
                  "in disconnected pieces the bounding box is used instead and you are told."),
                " The support map below shows the region in use."),
              uiOutput("spatial_window_readout"))
          )),

          # -- Labelling-density normalizsation (opt-in) ---------------------------
          fluidRow(
            column(3,
              h4("Labelling Density", .exploratory_badge()),
              .method_box("amber",
                HTML("<b>Independent random thinning</b>: keep each unit with a fixed probability. This leaves the pair correlation g(r) and Ripley's K unchanged in expectation (<b>Illian et&nbsp;al. 2008</b>; <b>Chiu et&nbsp;al. 2013</b>), so genuine clustering survives while a density-only difference is removed. Support area comes from a smoothed-density region (<b>Hyndman 1996</b>) or a morphological occupancy estimate (<b>Serra 1982</b>). The area-free density check uses nearest-neighbour distances (<b>Thompson 1956</b>)."),
                HTML("Whole <b>units</b> are dropped (a trace ID, or a de-blinked molecule with all of its localizations), never individual localizations. Dropping single localizations would shorten traces and corrupt de-blinking and per-trace spread. The default unit is the <b>de-blinked molecule</b>, since a trace ID is one blink and trace-ID density mixes expression level with photophysics. Density divides by the estimated support rather than the bounding box, because a MINFLUX field routinely contains the cell edge and interior dark regions. The draw is seeded per file and reproduces on re-run."),
                HTML("Use this only where density varies for reasons outside the biology, such as expression level, labelling efficiency or acquisition length. When density is itself the biology, for example a concentration-dependent process, equalising it removes the effect you are measuring. Files already at comparable density need no thinning, and Ripley's L and g(r) need none either since they divide density out. Thinning is one-way (nothing is up-sampled) and can push a small cluster below minPts. A difference that survives at matched density is evidence, and how clustering <i>scales</i> with density is the stronger test.")),
              helpText(HTML("For cells whose labelling density varies for reasons outside the biology. Leave off otherwise.")),
              .detail(HTML("Transient transfection, expression level and labelling efficiency all move density independently of the biology. Nearest-neighbour distances shrink as density rises, and DBSCAN at a fixed eps and minPts is itself a density threshold, so both read a high expresser as more clustered even when nothing has changed.")),
              checkboxInput("dens_norm_enable", "Measure density per file", value = FALSE),
              conditionalPanel(
                condition = "input.dens_norm_enable == true",
                selectInput("dens_norm_unit", "Unit",
                            choices = c("Auto (molecules if de-blinking on)" = "auto",
                                        "De-blinked molecules" = "cluster",
                                        "Trace IDs"            = "tid"),
                            selected = "auto"),
                checkboxInput("dens_norm_3d", "3-D (per µm³)", value = FALSE),
                numericInput("dens_nn_k", "k for local density / NN check", value = 3,
                             min = 1, max = 20, step = 1),

                tags$hr(style = "margin:8px 0;"),
                selectInput("dens_support_mode", "Support",
                            choices = c("Smoothed density, 99% region (recommended)" = "kde",
                                        "Occupancy grid (edge & dark-region aware)"  = "grid",
                                        "Convex hull"  = "hull",
                                        "Bounding box" = "bbox"),
                            selected = "kde"),
                helpText(HTML("The area density is divided by. A field usually contains the <b>cell edge</b> and interior <b>dark regions</b>; the bounding box counts those as sample.")),
                radioButtons("dens_support_from", "Outline drawn from",
                             choices = c("Trace / molecule centroids" = "units",
                                         "All localizations"          = "locs"),
                             selected = "units", inline = TRUE),
                helpText("Which points define the cell edge. Centroids are what the analyses count; localizations are a denser sample of the same region."),
                .detail(HTML(
                  "Neither choice is free of a bias, so this is a real setting rather than a default to leave alone.<br><br>",
                  "<b>Centroids</b> are one point per trace, or per de-blinked molecule, which is what every analysis on this tab counts. A sparsely labelled cell gives a raggeder, smaller footprint, so the area carries a little of the unit count: on synthetic pairs with identical geometry and a 4&times; labelling difference the supports came out 5.01e6 against 3.82e6&nbsp;nm&sup2;, reading a 4.0&times; density ratio as 3.05&times;.<br><br>",
                  "<b>All localizations</b> sample the same region more densely, so the outline is smoother. But the number of localizations in a trace is how long that fluorophore stayed on before going dark, which is photophysics, and it feeds the gap-filling sweep. On the membrane-protein set this inflated the area enough to loosen the randomisation-calibrated eps, and the largest-focus size stopped separating the two conditions (10&ndash;24 against 19&ndash;52 units, where centroids gave 6&ndash;12 against 19&ndash;51).")),
                numericInput("dens_grid_h", "Grid cell (nm)", value = 100, min = 10, step = 10),
                conditionalPanel(
                  condition = "input.dens_support_mode == 'kde'",
                  numericInput("dens_kde_bw", "Smoothing bandwidth (nm)", value = 150, min = 10, step = 25),
                  numericInput("dens_kde_frac", "Fraction of units enclosed", value = 0.99,
                               min = 0.5, max = 1, step = 0.01),
                  helpText(HTML("Keeps the smallest region containing this share of the units. <b>Keep the fraction and bandwidth identical across every sample you compare.</b>")),
              .detail(HTML("Smooths the localizations into a density map and takes the highest-density region (<b>Hyndman 1996</b>). Dark regions fall below the cut on their own, so there is no fill radius to tune. The fraction and the bandwidth together define what \"the sample area\" means, which is why changing them per file makes the areas incomparable."))
                ),
                conditionalPanel(
                  condition = "input.dens_support_mode == 'grid'",
                  checkboxInput("dens_close_auto", "Auto gap fill", value = TRUE),
                  helpText(HTML("<b>The most important setting here.</b> Set it too small and the measured area traces your structures instead of your cell.")),
              .detail(HTML("Gaps have two causes and this radius tells them apart: <b>sparse sampling</b> inside the cell should be filled in, because it is still sample, while a <b>dark region</b> or the space outside the cell should stay excluded. Auto picks a radius from each file and applies one common value to all of them.")),
                  conditionalPanel(
                    condition = "input.dens_close_auto == false",
                    numericInput("dens_close_nm", "Fill gaps up to (nm)", value = 200, min = 0, step = 50)
                  ),
                  helpText(HTML("<b>The gap fill is also what protects dark regions.</b> It can only close a hole narrower than roughly twice its radius, so anything wider than that stays excluded. If dark regions in your cell are being swallowed, lower the radius and check the result on the support map.")),
                  checkboxInput("dens_show_sweep", "Show gap-fill sweep", value = FALSE)
                ),
                checkboxInput("dens_show_map", "Show support map", value = FALSE),
                helpText(HTML("Counted area in <b style=\"color:#1a8a1a\">green</b>, disregarded in <b style=\"color:#c03030\">red</b>. Check the green region looks like your cell, not like a tracing of your structures.")),
              .detail(HTML("Works for every support mode. The bounding box is always the most generous and the occupancy grid the least.")),
                conditionalPanel(
                  condition = "input.dens_show_map == true",
                  selectInput("dens_map_file", "Dataset to render", choices = NULL)
                ),
                checkboxInput("dens_plot_log", "Log density axis", value = TRUE),
                helpText("Local density can span a wide range within a single cell, and a log axis keeps the whole distribution readable. Untick if your values are tightly grouped."),

                tags$hr(style = "margin:8px 0;"),
                checkboxInput("dens_thin_enable", "Thin all files to a common density", value = FALSE),
            tags$div(style = "font-size:0.82em; color:#7a5c00; margin:2px 0 6px; padding:4px 8px; border-left:3px solid #e0c890; background:#fffdf6;",
              tags$b("Only needed for count-based read-outs."), " Thinning is what makes ",
              "DBSCAN cluster counts, % in clusters, cluster sizes and nearest-neighbour ",
              "distances comparable between files, because all of those change with how ",
              "densely the sample is labelled.", tags$br(),
              tags$b("Ripley's L and g(r) do not need it"), ", since they already divide the ",
              "density out, so thinning them discards data without changing the answer. ",
              "Leave this off unless you are reporting counts."),
                conditionalPanel(
                  condition = "input.dens_thin_enable == true",
                  radioButtons("dens_target_mode", "Target",
                               choices = c("Lowest in the set"      = "min",
                                           "Control group's median" = "control",
                                           "Manual"                 = "manual"),
                               selected = "min"),
                  helpText("Lowest-in-set never needs up-sampling, so it is the default."),
                  conditionalPanel(
                    condition = "input.dens_target_mode == 'control' && input.analysis_by == 'group'",
                    selectInput("dens_target_group", "Control group", choices = NULL)
                  ),
                  conditionalPanel(
                    condition = "input.dens_target_mode == 'manual'",
                    numericInput("dens_target_manual", "Target (units/µm² or µm³)",
                                 value = NA, min = 0, step = 1)
                  )
                ),
                actionButton("run_density", "Measure density", icon = icon("play"),
                             class = "btn-primary btn-sm"),
                downloadButton("download_density", "Download", class = "btn-sm")
              )
            ),
            column(9,
              conditionalPanel(
                condition = "input.dens_norm_enable == true",
                uiOutput("density_status"),
                conditionalPanel(
                  condition = "output.dens_norm_active == true",
                  tags$div(style = "margin:6px 0; padding:5px 10px; border-left:4px solid #c03030; background:#fdf3f3; font-size:0.82em;",
                    tags$b("Thinning active."), " Analyses below run on thinned data (after filters, before clustering). Quality-Control panels are exempt.")
                ),
                plotlyOutput("density_plot", height = "290px"),
                conditionalPanel(
                  condition = "input.dens_show_map == true",
                  tags$div(style = "margin-top:10px;",
                           plotlyOutput("density_map", height = "460px"))
                ),
                conditionalPanel(
                  condition = "input.dens_show_sweep == true && input.dens_support_mode == 'grid'",
                  tags$div(style = "margin-top:10px;",
                           plotlyOutput("density_sweep_plot", height = "330px"))
                ),
                tags$div(style = "margin-top:10px;", DT::dataTableOutput("density_table"))
              )
            )
          ),
          tags$hr()
          )
        )
      )),

      # -- Dataset Groups (shared across analyses) ---------------------------
      wellPanel(
        h4("Dataset Groups"),
        radioButtons("analysis_by", "Analyze by",
                     choices = c("Individual files" = "file", "Dataset groups" = "group"),
                     selected = "file", inline = TRUE),
        conditionalPanel(
          condition = "input.analysis_by == 'group'",
          helpText("Assign loaded files to named groups. In group mode the analyses below pool files",
                   "within each group (one series per group) so you can compare groups.",
                   "Spatial statistics (Ripley/PCF) will average per-file curves."),
          helpText(style = "color:#7a5c00;",
                   "Distribution histograms pool units (distances, clusters) from every file in a group.",
                   "Units from one cell are correlated, so a test on pooled units overstates n.",
                   "For statistics between groups use one value per file from the per-file summary."),
          downloadButton("download_per_file_summary", "Download per-file summary (CSV)", class = "btn-sm"),
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

      # -- SECTION: Quality Control --
      fluidRow(column(12, tags$h3('Quality Control',
        style = "border-bottom:2px solid #d7dde3; padding-bottom:4px; margin-top:18px; color:#2c3e50;"))),

      # -- QC: Localization precision cloud ----------------------------------
      fluidRow(
        column(3,
          h4("Localization Precision & Sample drift estimation"),
          .method_box("green",
            HTML("Every localization is plotted relative to its own trace centroid, pooled at (0,0). The RMS radius is σ_loc, drawn as 1σ and 2σ ellipses with a separate σ for each axis (σ<sub>x</sub>, σ<sub>y</sub>, σ<sub>z</sub>)."),
            HTML("Runs on the post-filtered localizations set in the Display tab (metadata, length and spread filters, plus any fiducial exclusion). MINFLUX has <b>no fitted PSF</b>. The doughnut is the excitation intensity pattern, and a position comes from the ratio of photon counts collected at a few known beam positions around its central zero. Precision is set by the size of that targeting pattern and by the photon count."),
            HTML("σ is an <b>upper bound on precision</b>, since it also contains label and linker motion, molecule motion and any residual drift within the trace. Anisotropy means σ differs by direction, and the tight axis is genuinely tight. Axes are reported separately and the rings are drawn as ellipses so a single pooled σ never masks that. A tilted cloud has correlated x and y errors. A tilt orientation shared across many traces points to something fixed in the instrument or to residual, directional drift. A tilt that is radial with respect to the beam pattern is expected, since MINFLUX precision depends on where the molecule sits inside the targeting pattern. In 3D the axial and lateral directions are probed by different beam geometry, so σ<sub>z</sub> and σ<sub>xy</sub> differ in either direction.")),
          .method_box("amber",
            HTML("<b>Sample drift and jitter.</b> Drift is estimated from the displacement of trace-ID centroids between detected blinks, sorted by time. The de-blinked option assumes the same molecule blinked and was detected more than once."),
            HTML("Set your de-blinking to a satisfactory level first; otherwise the tool de-blinks by proximity on its own. Within each de-blinked group the first detected trace ID is the origin (0,0,0), and later trace IDs are regressed against time from there. <b>Drift</b> is the part that lines up with time and accumulates. <b>Jitter</b> is the scatter left once drift is removed. Two different molecules merged into one group scatter evenly about zero, so they widen the interval while leaving the slope flat."),
            HTML("Grouping assumes one molecule blinked and was detected several times. If the sample moved further than the grouping radius, blinks of one molecule stop pairing up and the largest displacements drop out, so a large drift is reported as a lower bound. Fiducials that lock the frame, optionally with beamline monitoring, already correct for sample drift.")),
          helpText("Pooled scatter of every localization about its trace centroid, for the active dataset.",
                   "Reports per-axis σ, anisotropy, skew and kurtosis from the filtered data."),
          numericInput("prec_cloud_min_loc", "Min localizations / trace", value = 5, min = 3, step = 1),
          checkboxInput("prec_cloud_use_z", "Include z (report σ_z)", value = FALSE),
          actionButton("run_prec_cloud", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_prec_cloud", "Download σ", class = "btn-sm"),
          checkboxInput("view_drift_plot", "Show sample drift plot", value = FALSE),
          conditionalPanel(
            condition = "input.view_drift_plot == true",
            helpText("Each molecule that blinked more than once starts at the origin; every later ",
                     "blink is one point, positioned relative to that molecule's first blink. ",
                     "One point = one trace ID that was merged into that de-blinked ID."),
            radioButtons("drift_plot_color", "Colour points by", inline = TRUE,
                         choices = c("Time since first blink" = "time",
                                     "Acquisition time" = "abs",
                                     "De-blinked ID" = "mol"),
                         selected = "time"),
            conditionalPanel(
              condition = "input.drift_plot_color == 'abs'",
              helpText(style = "font-size:0.78em; color:#888;",
                       "Blue early in the run through to red at the end. If displacements from one ",
                       "part of the acquisition favour one direction and those from another part ",
                       "favour another, the drift changed course during the run.")
            ),
            conditionalPanel(
              condition = "input.drift_plot_color == 'mol'",
              helpText(style = "font-size:0.78em; color:#888;",
                       "Blinks of one molecule share a colour and are joined from the origin in ",
                       "time order, so each molecule draws its own small track. Colours repeat and ",
                       "only group points; a colour does not identify a molecule.")
            ),
            helpText(style = "font-size:0.78em; color:#888;",
                     "A cloud centred on the origin means the sample held still. Points marching ",
                     "away in one direction, further the longer the gap, mean it drifted."),
                   .detail("The heavy markers are the geometric median displacement in each time ",
                     "bin. Jitter cancels there, so under a still sample every bin sits on the origin."),
            .detail("The geometric median is used because a merged molecule occasionally holds two ",
                     "different fluorophores, and it stays fixed if the x and y axes are relabelled.")
          )
        ),
        column(9,
          plotlyOutput("prec_cloud_plot", height = "340px"),
          uiOutput("prec_cloud_stats"),
          conditionalPanel(
            condition = "input.view_drift_plot == true",
            plotlyOutput("prec_drift_plot", height = "330px")
          )
        )
      ),

      tags$hr(),


      # -- QC: Fiducial / bright-contaminant flagging (report only) ----------
      fluidRow(
        column(3,
          h4("Fiducial / Contaminant Scan"),
          .method_box("amber",
            HTML("Flags traces that behave like gold fiducials or bright aggregates. <b>LIFETIME</b> is the primary signature: a fiducial does not bleach, so the sequence keeps re-finding it and its trace outlives every real fluorophore in the same acquisition. <b>BRIGHTNESS</b> is the second: gold reflects far more than a dye, though on some samples the bead is not the brightest thing in the file."),
            HTML("The lifetime reference is the file&#39;s own 99th percentile of the traces that live at all, so it adapts to the acquisition instead of assuming a number of seconds, and a file of mostly single-localization traces cannot collapse it. Brightness is a per-file robust z of log efo and additionally requires the trace to have persisted, since a single very bright localization is a blip. Length-outlier, recurrence and spread are tabulated as evidence. The rule is calibrated against labelled acquisitions with and without beads, from several samples and strains."),
            HTML("Reporting only by default. It tabulates suspects and removes nothing until you tick the box below.")),
          helpText("Flags traces that outlive the rest of the file, and bright traces that also persist.",
                   "Review the scatter and table, then optionally exclude them."),
          numericInput("fid_dur_x",    "Lifetime, as a multiple of the file's 99th percentile",
                       value = 5, min = 1.5, max = 50, step = 0.5),
          numericInput("fid_efo_z",    "Brightness outlier (robust z)", value = 4, min = 2, max = 10, step = 0.5),
          numericInput("fid_min_locs", "Min localizations for the brightness rule", value = 5, min = 1, max = 100, step = 1),
          .detail("The lifetime rule is the primary one: a fiducial does not bleach, so the tracker ",
                  "keeps re-finding it and its trace outlives every real fluorophore in the same file. ",
                  "The reference is that file's own 99th percentile of the traces that live at all, so it ",
                  "adapts to the acquisition rather than assuming a number of seconds. Lowering the ",
                  "multiple below about 5 starts to catch ordinary long traces. ",
                  tags$b("Brightness"), " is the second rule, and it only applies to a trace that ",
                  "persisted: a single very bright localization is a blip rather than a bead. ",
                  "Length-outlier and recurrence scores are tabulated as evidence."),
          numericInput("fid_len_z",    "Length outlier (reported only)",     value = 4, min = 2, max = 10, step = 0.5),
          numericInput("fid_recur_r",  "Recurrence radius (nm, reported only)", value = 100, min = 10, max = 1000, step = 10),
          numericInput("fid_recur_min","Min co-located traces (reported only)", value = 4, min = 1, max = 50, step = 1),
          checkboxInput("fid_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          checkboxInput("fid_exclude", "Exclude flagged traces from analyses", value = FALSE),
          helpText(style = "margin-top:-6px; font-size:0.78em; color:#888;",
                   "Reversible. While ticked, the flagged traces (per file) are dropped from the",
                   "display and the Analysis tools. Re-run the scan after changing thresholds."),
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

      # -- Trace Detection over Time (unique TIDs vs acquisition time) --------
      fluidRow(
        column(3,
          h4("Trace Detection over Time"),
          .method_box("green",
            HTML("Counts <b>unique trace IDs</b> per acquisition-time bin. A trace ID is one continuous detection of a blinking fluorophore, so this is the <b>detection-event rate</b> over the run."),
            HTML("Each trace is counted once at its median <code>tim</code>, so long traces do not crowd the early bins. Multiple files are drawn as <b>one series each</b>, or pooled per group in group mode. Each dataset has its own tim timeline, so its line spans only that window and <b>stops when the tim ends</b>. The per-localization metadata filters apply. The trace-length filter is left off, since short blinking events are real detections."),
            HTML("This is descriptive QC. A declining rate over time reflects photobleaching, exhausted blinking, drift and label or sample depletion. Read it as a check on acquisition stability.")),
          helpText("Unique traces counted once at their median acquisition time,",
                   "one series per file (or per Dataset Group in group mode)."),
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
          helpText("Normalized modes scale each series to its own total, comparing detection",
                   "timing (shape) across files or groups of different length or brightness."),
          checkboxInput("tid_time_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          actionButton("run_tid_time", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_tid_time", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("tid_time_plot", height = "400px")
        )
      ),

      tags$hr(),


      # -- SECTION: Nearest-Neighbour & Spatial Statistics --
      fluidRow(column(12, tags$h3('Nearest-Neighbour & Spatial Statistics',
        style = "border-bottom:2px solid #d7dde3; padding-bottom:4px; margin-top:18px; color:#2c3e50;"))),

      # -- Batch NN ----------------------------------------------------------
      fluidRow(
        column(3,
          h4("Batch NN Distance Histogram"),
          .method_box("green",
            HTML("Distribution of the k-th nearest-neighbour distance, computed with a k-d tree (<code>dbscan::kNN</code>). This is the classic first-order description of how tightly points sit together, following <b>Clark &amp; Evans 1954</b>. The panel reports the distance distribution itself, with no null model attached (see the caution)."),
            HTML("The choice of <b>NN unit</b> is the MINFLUX adaptation. One fluorophore re-appears as several traces, and a multivalent label scatters dyes ~10–40&nbsp;nm around one epitope, so raw localizations over-count molecules. Running NN on <b>TID centroids</b> (one point per trace) or <b>de-blinked-molecule / median-cluster centroids</b> de-blinks first, so the distribution approximates inter-<i>molecule</i> spacing. At MINFLUX's nm precision the histogram is set by labelling geometry."),
            HTML("A sharp peak below ~40&nbsp;nm usually reflects the <b>label footprint</b> (Annibale et&nbsp;al. 2011; Sengupta et&nbsp;al. 2011). The distance-range filter trims the histogram without correcting that bias. This is a summary statistic with no CSR reference. For significance, use Ripley&nbsp;L(r) with its envelope.")),
          selectInput("nn_batch_unit", "NN unit",
            choices = c(
              "Auto (molecules if de-blinking is on)"  = "auto",
              "TID centroids (all files)"             = "tid",
              "De-blinked molecule centroids"             = "pcid",
              "TID centroids within Median Clusters"         = "cluster2_tid",
              "De-blinked molecule centroids within Median Clusters" = "cluster2_pcid"
            ),
            selected = "auto"
          ),
          # Channel facet. Appears only once the data has channels; a cross-channel
          # run ignores it, because a statistic about the PAIR cannot be split by
          # the members of the pair.
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("analysis_by_channel",
                          "One distribution per channel", value = FALSE),
            helpText("One distribution per channel, labelled 'file · channel'."),
            .detail("Ignored for the cross-channel modes below. Ripley's L and g(r) have ",
                    "their own per-channel option, because they divide by the observed area ",
                    "and need a shared window.")
          ),
          helpText("Runs on every loaded file, using the Display tab's filters and clustering settings."),
          .detail("Auto follows the Display tab's de-blinking switch, so turning de-blinking on ",
                  "applies it here and to every dataset. Median Cluster modes require Median ",
                  "Clustering to be enabled."),
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
          uiOutput("nn_batch_stale"),
          plotlyOutput("nn_batch_hist", height = "540px"),
          tableOutput("nn_batch_stats")
        )
      ),

      tags$hr(),

      # -- Ripley's L --------------------------------------------------------
      fluidRow(
        column(3,
          h4("Ripley's L(r)"),
          .method_box("green",
            HTML("Ripley's K-function (<b>Ripley 1977</b>) with the variance-stabilising L-transform (<b>Besag 1977</b>). For <i>n</i> points in a window of area <i>A</i>, the (uncorrected) estimator is<br><span style=\"font-family:serif; margin-left:6px;\">K̂(r) = (A / [n(n−1)]) · Σ<sub>i≠j</sub> 𝟙(d<sub>ij</sub> ≤ r)</span><br>i.e. area × (fraction of the n(n−1) ordered point pairs closer than r). The L-transform then linearises it: <span style=\"font-family:serif;\">L(r) = √(K(r)/π)</span> in 2D, <span style=\"font-family:serif;\">L(r) = (3·K(r)/4π)<sup>1/3</sup></span> in 3D. Plotted as L(r)&nbsp;−&nbsp;r: 0 under complete spatial randomness, &gt;0 clustering, &lt;0 dispersion. The SMLM lineage of this statistic runs through PC-PALM (<b>Sengupta et&nbsp;al. 2011</b>; <b>Veatch et&nbsp;al. 2012</b>)."),
            HTML("Runs on <b>de-blinked units</b> (trace-ID or de-blinked molecule centroids) so that repeated appearances of one molecule do not manufacture a short-r clustering peak. Points are randomly subsampled to ≤2000 per file (pairwise distances are O(n²)); the window is the axis-aligned bounding box. The optional <b>Monte-Carlo CSR envelope</b> simulates n random points in the <i>identical frozen window</i>, so the band and the observed curve carry the same edge bias and it cancels."),
            HTML("Edge correction (top of the tab, on by default) weights each pair by the translation correction. With it off, points close to the window boundary lose neighbours and L−r bends downward near the top of the r range, and only the envelope is a valid reference. The shaded envelope is pointwise: it holds at one radius, so a curve that leaves it somewhere is not significant at 5%. Read significance from the dashed <b>global band</b> and the p-value above the plot, which cover the whole r range in one test. Without an envelope, read shape and peak position. The dominant biological confound is labelling: incomplete labelling flattens real clustering, multivalent labels fake it. In group mode per-file curves are averaged, never pooled coordinates, since pooling would mix windows and invent structure.")),
          helpText("L(r) − r = 0 under complete spatial randomness; > 0 is clustering, < 0 is regularity."),
          .detail("Without an envelope, read deviations relative to the point density and window ",
                  "size, and switch the CSR envelope on below for significance. The K estimator is ",
                  "uncorrected unless edge correction is on, so curves may dip below zero at large r ",
                  "from border effects. Large datasets are subsampled to 2,000 points per file."),
          selectInput("ripley_unit", "Unit",
            choices = c(
              "Auto (molecules if de-blinking is on)" = "auto",
              "TID centroids"             = "tid",
              "De-blinked molecule centroids" = "cluster"
            ),
            selected = "auto"
          ),
          helpText("Which points the statistic runs on. Auto follows the Display tab's de-blinking",
                   "switch and applies it to every dataset, so the unit never silently disagrees",
                   "with the other panels."),
          checkboxInput("ripley_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          helpText("All files: pool every loaded dataset into one curve; uncheck to analyze only the active file."),
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("ripley_by_channel", "One curve per channel", value = FALSE),
            helpText("One curve per channel, both measured in the same observation window."),
            .detail("Covers a 2-colour DCR split and a multiplex build alike. The window comes ",
                    "from the pooled file, so a difference between the curves reflects how the two ",
                    "targets are arranged in that field, independent of one channel having fewer ",
                    "localizations than the other.")
          ),
          checkboxInput("ripley_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          helpText("3-D: measure neighbours in spheres of radius r (uses z) instead of 2-D circles.",
                   "Use only for genuinely 3-D data; the CSR baseline differs between 2-D and 3-D."),
          numericInput("ripley_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          helpText("r max: largest spatial scale (radius) evaluated. 0 auto-picks roughly a quarter of the",
                   "field size. Read clustering peaks at the r where L(r)-r is highest = the cluster scale."),
          numericInput("ripley_r_steps", "r steps", value = 60, min = 10, max = 300),
          helpText("r steps: how many radii between 0 and r max are computed. Higher gives a smoother curve and runs slower."),
          checkboxInput("ripley_smooth", "Smooth curve (loess)", value = TRUE),
          numericInput("ripley_smooth_span", "Smoothing span", value = 0.2, min = 0.05, max = 0.8, step = 0.05),
          helpText("Smoothing: a loess pass over the L(r)-r curve to reduce sampling noise. Span = fraction",
                   "of points in each local window (higher = smoother but can flatten real peaks). Purely",
                   "cosmetic, and it does not change the underlying statistic."),
          conditionalPanel(
            condition = "input.analysis_by == 'group'",
            radioButtons("ripley_group_band", "Group spread band",
                         choices = c("IQR (Q25-Q75)" = "iqr", "±SD" = "sd", "None" = "none"),
                         selected = "iqr", inline = TRUE),
            checkboxInput("ripley_show_files", "Show individual files (shaded by group)", value = FALSE),
            helpText("Draws each file as a thin shade of its group's colour behind the group mean,",
                     "so a group difference can be judged against its own replicate spread."),
            helpText("One mean curve per group; the band is between-file spread, not a test."),
            .detail("Each file's L(r)−r is computed on its own window, then averaged within the ",
                    "group; coordinates are never pooled. Curves are interpolated onto a common ",
                    "r-grid, up to the smallest per-file r max.")
          ),
          checkboxInput("ripley_envelope", "Show CSR envelope (Monte-Carlo)", value = FALSE),
          helpText(HTML("The shaded band is <b>pointwise</b> (one radius at a time), so leaving it somewhere is not a 5% result. The dashed <b>global band</b> and the p-value above the plot test the whole r range at once: above it is clustering, below it is dispersion.")),
          .detail("Simulates CSR patterns in the same window as the data and shades the 2.5–97.5% ",
                  "band. The band captures the edge bias, so it is the correct reference. In 'all ",
                  "files' mode the envelope is computed for the first six files only, since the cost ",
                  "is simulations × files. ",
                  tags$b("In group mode"), " the envelope is drawn as a dotted outline in the ",
                  "group's colour (the mean of its files' envelopes), separate from the solid ",
                  "between-file spread band. The spread band shows how much your cells differ from ",
                  "each other; the envelope shows what randomness would look like. The legend reports ",
                  "how many files the envelope covers when that is fewer than the group."),
          conditionalPanel(
            condition = "input.ripley_envelope == true",
            numericInput("ripley_nsim", "Envelope simulations", value = 99, min = 19, max = 199)
          ),
          actionButton("run_ripley", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_ripley", "Download", class = "btn-sm")
        ),
        column(9,
          uiOutput("ripley_stale"),
          plotlyOutput("ripley_plot", height = "540px")
        )
      ),
      tags$hr(),

      # -- Univariate pair-correlation g(r) ----------------------------------
      fluidRow(
        column(3,
          h4("Pair Correlation g(r)"),
          .method_box("green",
            HTML("The pair-correlation (radial distribution) function, the <b>derivative</b> form of Ripley's K, counting point pairs in a shell (r−dr/2, r+dr/2] normalized by the shell area (2D) or volume (3D):<br><span style=\"font-family:serif; margin-left:6px;\">ĝ(r) = A · 2·N<sub>pairs</sub>(shell) / [n(n−1) · shell]</span><br>g&nbsp;=&nbsp;1 under complete spatial randomness, &gt;1 clustering at scale r, &lt;1 exclusion. In SMLM this is the PC-PALM statistic (<b>Sengupta et&nbsp;al. 2011</b>; <b>Veatch et&nbsp;al. 2012</b>), developed precisely to separate repeated appearances of one molecule from genuine clustering. Theory: <b>Illian et&nbsp;al. 2008</b>; <b>Stoyan &amp; Stoyan 1994</b>."),
            HTML("The one-colour version of the g₁₂(r) below, run on the same de-blinked units as Ripley's L. g divides out the intensity twice over, so it is insensitive to how densely the sample is labelled and suits comparing cells of different expression level. Each radius is measured independently, so you read the aggregate size off <b>where g comes back down to 1</b>, and a short-range blinking or label artefact stays in the first few bins."),
            HTML("Points near the window boundary lose some of their neighbours, which pulls g downward at larger radii. Switch on edge correction (top of the tab) to remove that. Without it, judge the curve against the Monte-Carlo envelope, which carries the same bias, in place of the flat g&nbsp;=&nbsp;1 line. Significance comes from the dashed global band and the p-value above the plot, never from the shaded pointwise band alone. g is noisier than L at short radii, because a thin ring holds far fewer pairs than a filled disc, so widen the bins or smooth before reading much into a single spike. In group mode the per-file curves are averaged; coordinates are never pooled.")),
          helpText(HTML("The <b>r at which g falls back to 1</b> is the characteristic aggregate size; the height of g is how many times denser than random that neighbourhood is.")),
          .detail(HTML("A solid aggregate makes g <i>decrease monotonically</i> from the smallest r. An interior peak (g low at short r, rising to a maximum further out) means <b>short-range exclusion</b>: molecules that cannot approach each other closer than some distance, so pairs pile up at that spacing. A steric footprint, a lattice, or a regular array all do this.<br><br>A hollow or ring-shaped structure does not produce it. Points spread around a ring still sit arbitrarily close to their neighbours along the arc, so a ring keeps its short-r excess and g decays monotonically (verified in simulation up to 300 points per ring). Only depletion at short r produces an interior peak.")),
          selectInput("pcf_unit", "Unit",
            choices = c("Auto (molecules if de-blinking is on)" = "auto",
                        "TID centroids" = "tid",
                        "De-blinked molecule centroids" = "cluster"),
            selected = "auto"),
          checkboxInput("pcf_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("pcf_by_channel", "One curve per channel", value = FALSE),
            helpText("Splits each file into its channels, both measured in the same ",
                     "observation window (the pooled file).")
          ),
          checkboxInput("pcf_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          helpText("3-D: count pairs in spherical shells instead of 2-D annuli. Use only for",
                   "genuinely 3-D data."),
          numericInput("pcf_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          numericInput("pcf_r_steps", "r steps", value = 50, min = 10, max = 300),
          helpText("Fewer steps = wider shells = more pairs per bin = a less noisy curve, at coarser",
                   "spatial resolution. 40-60 is usually the right trade-off."),
          checkboxInput("pcf_smooth", "Smooth curve (loess)", value = TRUE),
          numericInput("pcf_smooth_span", "Smoothing span", value = 0.25, min = 0.05, max = 0.8, step = 0.05),
          helpText("Cosmetic only. It does not change the statistic or the summary numbers."),
          fluidRow(
            column(6, numericInput("pcf_win_lo", "Summary window lo (nm)", value = 0, min = 0, step = 10)),
            column(6, numericInput("pcf_win_hi", "hi (nm)", value = 0, min = 0, step = 10))
          ),
          helpText("The window the summary mean and integrated excess are taken over. Leave both at",
                   "0 for an automatic 5-40% of r max. Set it to the scale you want to compare."),
          conditionalPanel(
            condition = "input.analysis_by == 'group'",
            radioButtons("pcf_group_band", "Group spread band",
                         choices = c("IQR (Q25-Q75)" = "iqr", "±SD" = "sd", "None" = "none"),
                         selected = "iqr", inline = TRUE),
            checkboxInput("pcf_show_files", "Show individual files (shaded by group)", value = FALSE)
          ),
          checkboxInput("pcf_envelope", "Show CSR envelope (Monte-Carlo)", value = FALSE),
          helpText(HTML("<b>Reading the band.</b> It is a spread around random, and <b>nothing is subtracted from your curve</b>. Significance is the curve leaving the band.")),
          .detail(HTML("g already sits at 1 for a random pattern. The band is widest at short radii because a thin ring contains very few pairs there, so a random pattern's g is genuinely noisier. A peak at short range has to clear a wider bar than a modest rise further out. Quote strength as g itself, or as the area between g and 1.")),
          conditionalPanel(
            condition = "input.pcf_envelope == true",
            numericInput("pcf_nsim", "Envelope simulations", value = 99, min = 19, max = 199)
          ),
          checkboxInput("pcf_fit_xi", "Fit correlation length ξ", value = FALSE),
          conditionalPanel(
            condition = "input.pcf_fit_xi == true",
            numericInput("pcf_fit_rmin", "Fit from r (nm)", value = 40, min = 0, step = 5),
            helpText("Fits g(r) - 1 = A·exp(-r/ξ), giving one number per curve. ",
                     "Start above the labelling floor: below a few tens of nm the curve is ",
                     "dominated by the label and by residual blinking, not by structure.")
          ),
          actionButton("run_pcf", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_pcf", "Download", class = "btn-sm")
        ),
        column(9,
          uiOutput("pcf_stale"),
          plotlyOutput("pcf_plot", height = "540px"),
          uiOutput("pcf_summary"),
          conditionalPanel(
            condition = "input.pcf_fit_xi == true",
            tags$div(style = "margin-top:10px;",
              tags$b("Correlation length ξ"),
              tableOutput("pcf_xi_table"),
              uiOutput("pcf_xi_note"))
          )
        )
      ),
      tags$hr(),


      # -- Local Density & Hotspots ------------------------------------------
      fluidRow(
        column(3,
          h4("Local Density & Hotspots"),
          .method_box("amber",
            HTML("Three established descriptors of a point pattern, each computed per cell and expressed relative to that cell's own mean density. <b>Local density</b> comes from the k-th nearest-neighbour distance, &lambda;<sub>i</sub>&nbsp;=&nbsp;k&nbsp;/&nbsp;(&omega;<sub>d</sub>&nbsp;r<sub>k</sub><sup>d</sup>) (<b>Loftsgaarden &amp; Quesenberry 1965</b>). The <b>index of dispersion</b> is the variance-to-mean ratio of quadrat counts, 1 under complete spatial randomness (<b>Fisher, Thornton &amp; Mackenzie 1922</b>). <b>Foci</b> are DBSCAN clusters (<b>Ester et&nbsp;al. 1996</b>) at an eps calibrated per cell against a randomised pattern in that cell's own footprint."),
            HTML("Every quantity is divided by the cell's own mean density, so a cell that expressed twice as much protein reads the same. The footprint from the density block supplies that mean, and DBSCAN's eps is re-derived per cell instead of being a fixed setting, so a sparse cell and a dense cell are held to the same standard."),
            HTML("These read the <b>upper tail</b> of the local-density distribution, so they answer whether a minority of units sit in dense spots. They say nothing about the shape or the identity of those spots, and the counts inside them are detected units, which is a lower bound on molecules. <b>Number of foci</b> and <b>fraction of units in foci</b> stay density-sensitive even after calibration; the size of the largest focus does not. The per-cell eps calibration is specific to this app. Density unevenness, the index of dispersion, the fraction above the threshold and the largest focus relative to chance hold their random baseline whatever the labelling density and cell shape; raw focus counts and raw largest-focus size grow with density, so compare cells with <b>Largest focus (x random 95th pct)</b>. The index of dispersion depends on the quadrat size, so compare it only at a fixed quadrat.")),
          helpText("Measures each unit's local density as a multiple of its own cell's mean, and finds the dense spots. Density unevenness reads 1.0 for a random pattern."),
          .detail(HTML(
            "<b>Local-density enhancement.</b> For each unit, the density of units around it from the k-th nearest neighbour, divided by that cell's units-per-area over the footprint. 1&times; is the cell average, 20&times; is twenty times denser than the cell average.<br><br>",
            "<b>Density unevenness.</b> How much the local density varies across the cell: the standard deviation of the per-unit local densities divided by their mean, then divided again by the value a random pattern gives. <b>1.0 means random.</b> 2 means twice as uneven as random, 5 means a few dense spots on an otherwise sparse background.<br><br>",
            "The random baseline depends only on k, not on the density, the area or the number of units, which is what makes the number comparable between cells. It is 1/&radic;(k&minus;2): 1.00 at k&nbsp;=&nbsp;3, 0.577 at k&nbsp;=&nbsp;5, 0.354 at k&nbsp;=&nbsp;10. The raw coefficient of variation is kept in the table and the metric list for anyone who wants it unnormalised.<br><br>",
            "<b>Index of dispersion.</b> Counts in fixed quadrats, variance over mean. It is 1 for a random pattern by construction, so it reads without a simulation. It is scale-dependent, so set the quadrat near the feature size you are after and quote the size with the number.<br><br>",
            "<b>Calibrated eps.</b> DBSCAN with a fixed eps is a fixed density threshold, which a dense cell crosses more easily than a sparse one. This picks, per cell, the largest eps at which a randomised pattern in that cell's own footprint still clusters no more than the target fraction.")),
          selectInput("hot_unit", "Unit",
            choices = c("Trace ID centroids (one point per blink)" = "tid",
                        "De-blinked molecule centroids"           = "cluster",
                        "Auto (molecules if de-blinking is on)"   = "auto"),
            selected = "tid"),
          helpText("One point per trace, or per de-blinked molecule. Raw localizations are not offered."),
          .detail(HTML(
            "Localizations inside a single trace sit a few nanometres apart, so their local density measures how long that fluorophore stayed on before going dark. On the membrane-protein pair the enhancement then reaches 3,000 to 17,000 times the cell mean and neither the CV nor the 90th percentile separates the two conditions. The same statistics on trace-ID centroids separate them cleanly.<br><br>",
            "<b>Trace IDs</b> count blinking events, <b>de-blinked molecules</b> merge repeated blinks of one fluorophore. Both work here. Trace IDs gave the wider margin on the membrane-protein pair (90th-percentile enhancement 2.7&ndash;5.0 against 13.1&ndash;72.4, versus 2.1&ndash;5.0 against 6.6&ndash;53.1 for molecules), because de-blinking removes a few points from the densest spots.")),
          checkboxInput("hot_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          helpText("3-D measures density in spheres and quadrats as cubes. Use it for genuinely 3-D data; a thin membrane sheet is 2-D."),
          numericInput("hot_k", "k for local density", value = 5, min = 1, max = 50, step = 1),
          helpText("Neighbour order for the density estimate. Small k is noisy, large k blurs a focus into its surroundings."),
          .detail(HTML("k is a smoothing parameter. A focus holding fewer than about 2k units cannot be resolved at that k. 5 is a reasonable default for foci of tens of units; raise it for large diffuse domains.")),
          numericInput("hot_enh_thresh", "Enhancement threshold (x cell mean)",
                       value = 5, min = 1.5, max = 100, step = 0.5),
          helpText("The multiple that the reported fraction counts, and the dashed line on the distribution plot."),
          numericInput("hot_quadrat", "Quadrat size (nm)", value = 250, min = 20, step = 10),
          helpText("Side of the quadrats for the index of dispersion. Set it near the size of the features you are looking for."),
          tags$hr(style = "margin:10px 0;"),
          tags$div(style = "font-weight:600; color:#2c3e50; font-size:0.92em;", "Focus detection"),
          numericInput("hot_minpts", "minPts", value = 5, min = 2, max = 50, step = 1),
          helpText("Smallest number of units DBSCAN will call a focus. Held fixed across cells."),
          numericInput("hot_target_frac", "False-cluster target", value = 0.05,
                       min = 0.005, max = 0.5, step = 0.01),
          helpText("The fraction of a randomised pattern that may cluster. eps is set to the largest value that stays at or below this."),
          numericInput("hot_nsim", "Randomisations per eps", value = 15, min = 3, max = 100, step = 1),
          helpText("Randomised patterns drawn per candidate eps. 15 places the crossing between two grid steps."),
          tags$hr(style = "margin:10px 0;"),
          tags$div(style = "font-weight:600; color:#2c3e50; font-size:0.92em;", "Footprint"),
          tags$div(style = "font-size:0.85em; color:#445; line-height:1.5; margin-bottom:8px;",
            "The estimator mode, grid pitch and gap-fill radius come from ",
            tags$b("Density & area normalization"), " at the top of this tab, so this panel, ",
            "Ripley's L and g(r) all measure the same cell area. The gap-fill radius is resolved ",
            "once and shared across every loaded file."),
          .detail(HTML("The gap-fill radius is what sets the area, and the area sets the density every statistic here divides by. It should be a few times the typical spacing between units and well below the size of any dark region you want kept out of the cell. It is shared across files on purpose: a per-file radius would put each cell on its own density scale, which is the comparison this panel exists to make.<br><br>This panel always uses the footprint, whatever the footprint switch at the top of the tab says, because a bounding box would not merely weaken these numbers, it would change them: a cell imaged half in frame reads at half its true density and every enhancement doubles. Setting the estimator mode to <b>bbox</b> or <b>hull</b> has the same effect, since neither produces the cell mask that the calibrated eps and the dispersion index need.")),
          selectInput("hot_palette", "Map colour scale",
                      choices = NUMERIC_PALETTES,
                      selected = "Turbo"),
          checkboxInput("hot_show_foci", "Box the detected foci on the map", value = TRUE),
          conditionalPanel(
            condition = "input.hot_show_foci == true",
            numericInput("hot_box_top", "Box the largest N foci (0 = all)",
                         value = 5, min = 0, step = 1),
            helpText("Each box is drawn at the size of the focus it encloses and labelled with its count.")),
          actionButton("run_hotspot", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_hotspot", "Per file", class = "btn-sm"),
          downloadButton("download_hotspot_foci", "Per focus", class = "btn-sm")
        ),
        column(9,
          uiOutput("hotspot_stale"),
          uiOutput("hotspot_status"),
          fluidRow(
            column(7,
              selectInput("hot_map_file", "Map dataset", choices = NULL, width = "100%"),
              plotlyOutput("hotspot_map", height = "470px")),
            column(5,
              radioButtons("hot_curve_mode", NULL, inline = TRUE,
                choices = c("Upper tail" = "tail", "Cumulative" = "cdf"),
                selected = "tail"),
              plotlyOutput("hotspot_ecdf", height = "430px"))
          ),

          # This selector drives BOTH plots on the row below it, so it sits on that
          # row rather than under the curve it has nothing to do with.
          tags$div(style = "margin:14px 0 4px; padding:8px 12px; background:#f4f6f8; border:1px solid #e3e8ee; border-radius:4px;",
            fluidRow(
              column(5, selectInput("hot_metric", "Metric shown in both plots below",
                choices = c("Density unevenness (x random)"  = "cv_ratio",
                            "CV of local density (raw)"      = "cv_local",
                            "Largest focus (n)"              = "max_focus_n",
                            "Enhancement, 90th percentile"   = "enh_q90",
                            "Enhancement, 99th percentile"   = "enh_q99",
                            "Fraction above the threshold"   = "frac_above",
                            "Index of dispersion (quadrats)" = "vmr",
                            "Number of foci"                 = "n_foci",
                            "Fraction in foci"               = "frac_in_foci",
                            "Largest focus Rg (nm)"          = "max_focus_rg",
                            "Calibrated eps (nm)"            = "eps_nm",
                            "Median NN distance (nm)"        = "nn_med",
                            "Clark-Evans R"                  = "clark_evans",
                            "Density (per um2)"              = "density"),
                selected = "cv_ratio", width = "100%")),
              column(7, tags$div(style = "padding-top:26px; font-size:0.85em; color:#556;",
                "Left: the metric per cell, grouped. Right: the same metric against that cell's labelling density.")))
          ),
          fluidRow(
            column(6, plotlyOutput("hotspot_metric",     height = "330px")),
            column(6, plotlyOutput("hotspot_vs_density", height = "330px"))
          ),
          tags$div(style = "margin:6px 0 10px; padding:7px 12px; border-left:4px solid #cfe0ee; background:#f6f9fc; font-size:0.85em; color:#445;",
            tags$b("Reading the pair of plots. "),
            "Groups that separate vertically at matched density are separated by organisation. Groups that lie along one shared line are separated by expression level. ",
            tags$b("Median NN distance"), " and ", tags$b("Clark-Evans R"), " are in the metric list as the worked example of the second case."),
          tags$div(style = "overflow-x:auto;", tableOutput("hotspot_table"))
        )
      ),
      tags$hr(),

      # -- SECTION: Multi-panel view --
      fluidRow(column(12, tags$h3('Multi-Panel View',
        style = "border-bottom:2px solid #d7dde3; padding-bottom:4px; margin-top:18px; color:#2c3e50;"))),

      fluidRow(column(12,
        tags$div(style = "margin:0 0 10px; padding:7px 12px; border-left:4px solid #cfe0ee; background:#f6f9fc; font-size:0.86em; color:#445;",
          tags$b("Six slots, one analysis each. "),
          "Each slot draws the same figure its full-size panel draws, using that panel's controls and its stored results. ",
          "Run a panel before selecting it here; an unrun panel shows its own placeholder message."),
        fluidRow(uiOutput("panel_slot_pickers")),
        fluidRow(
          column(4, numericInput("panel_row_height", "Panel height (px)",
                                 value = 340, min = 200, max = 900, step = 20)),
          column(8, helpText("Figures keep their own legends and axes, so a slot is readable on its own."))
        ),
        uiOutput("panel_grid")
      )),
      tags$hr(),

      # -- SECTION: Median-Cluster Analyses --
      fluidRow(column(12, tags$h3('Median Cluster Analyses',
        style = "border-bottom:2px solid #d7dde3; padding-bottom:4px; margin-top:18px; color:#2c3e50;"))),

      # -- Traces per Cluster ------------------------------------------------
      fluidRow(
        column(3,
          h4("Trace IDs / De-blinked Molecules per Median Cluster"),
          .method_box("amber",
            HTML("Cluster-occupancy counting, in the lineage of SMLM molecular-counting and blink-correction work (<b>Annibale et&nbsp;al. 2011</b>, <i>Nat. Methods</i>, quantitative photo-activated localization microscopy; <b>Sengupta et&nbsp;al. 2011</b>, PC-PALM). The number of units the median cluster contains."),
            HTML("The grouping is always the <b>median cluster</b>. You choose what to count inside it: raw <b>trace IDs</b>, or <b>de-blinked molecules</b> (repeated blinks of one molecule collapsed into one unit; needs de-blinking on). Comparing the two shows how much of the raw count is blinking."),
            HTML("<b>This is an occupancy count.</b> A trace-ID count over-counts molecules by an unknown, dye- and buffer-dependent blinking factor. De-blinked molecules partly correct for that but still depend on the de-blink eps. Annibale's blink-correction needs a dark-time distribution that is not calibrated here. Read it as a relative comparison between conditions acquired identically, never as an absolute subunit count.")),
          helpText("Counts trace IDs or de-blinked molecules per Median Cluster, across all loaded files."),
          .detail("One fluorophore produces several trace IDs across blinking cycles, so a ",
                  "trace-ID count reflects blinking events. De-blinked molecules are closer to a ",
                  "molecule count."),
          radioButtons("tid_count_unit", "Count within each Median Cluster",
            choices = c(
              "Trace IDs"          = "tid",
              "De-blinked molecules" = "pcid"
            ),
            selected = "tid"
          ),
          radioButtons("tid_count_yaxis", "Y axis",
                       choices = c("Count" = "count", "Density" = "density"),
                       selected = "count", inline = TRUE),
          numericInput("tid_count_bins", "Bins", value = 30, min = 5, max = 200),
          checkboxInput("tid_count_aggregate", "Aggregate all datasets into one", value = FALSE),
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("tid_count_by_channel", "One distribution per channel", value = FALSE)
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

      # -- Cluster2 Detection over Time --------------------------------------
      fluidRow(
        column(3,
          h4("Median Clusters Detection over Time", .exploratory_badge()),
          .method_box("orange",
            HTML("<b>MINFLUXR original.</b> No published MINFLUX detection-timing estimator exists. The building blocks are standard: a detection-rate histogram over acquisition time, and a two-sample <b>Kolmogorov–Smirnov</b> comparison (max distance between two empirical cumulative distributions) of the clustered vs unclustered populations."),
            HTML("Each trace is counted <b>once</b>, at its median <code>tim</code>, so a long trace with many localizations does not dominate the early bins. The split is clustered (inside a median-cluster structure) vs unclustered background <i>within the same file</i>, so acquisition length, laser power and sample differences cancel in the comparison. The normalized modes scale each population to its own total, isolating a difference in <i>shape</i> (timing) from a difference in <i>amount</i> (abundance)."),
            HTML("<b>This is descriptive.</b> The drop-off is a convolution of photobleaching, blinking kinetics, labelling density and the tracker's own re-targeting, so crowding is one candidate cause among several. The KS gap is reported as a distance. The two samples are not independent, since they come from the same acquisition and the same molecules compete for one detection budget, so a KS p-value would not hold. The trace-detection-rate decay used to calibrate the simulator comes from in-house long acquisitions.")),
          helpText("How detection of clustered vs background traces evolves over the acquisition."),
          .detail("Each trace is counted once at its median acquisition time (tim), split into ",
                  "traces inside Median Cluster structures vs unclustered background. Shows ",
                  "photobleaching drop-off, and feeds the simulator's acquisition model. ",
                  "Requires Median Clustering."),
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
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("det_time_by_channel", "One curve per channel", value = FALSE),
            helpText("With a TIM-split multiplex build, a channel losing traces over its own ",
                     "window becomes visible here, in place of averaging into the whole acquisition.")
          ),
          actionButton("run_det_time", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_det_time", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("det_time_plot", height = "400px")
        )
      ),

      tags$hr(),

      # -- Cluster2 Crowding vs Detection-Time Shift -------------------------
      fluidRow(
        column(3,
          h4("Median Clusters Crowding vs Detection-Time Shift", .exploratory_badge()),
          .method_box("orange",
            HTML("<b>MINFLUXR original.</b> A correlational diagnostic built from two standard quantities: a local-density proxy (median k-NN distance among a cluster's de-blinked molecule centroids) and a per-cluster timing offset, with an optional loess trend (<b>Cleveland 1979</b>)."),
            HTML("Both axes are measured <b>per cluster and referenced within the same file</b>. The timing shift subtracts that file's unclustered-background median <code>tim</code>, so between-file differences in acquisition length, bleaching and labelling cancel. Crowding uses <b>de-blinked molecule centroids</b> (or trace-ID centroids when de-blinking is off), the <i>same</i> unit and k-NN definition as the Within-cluster-NN panel, so the two panels report the same spacing."),
            HTML("<b>This is a correlation.</b> Local density co-varies with labelling efficiency, structure size and local photophysics, so a crowding–timing trend is a hypothesis to test with controls such as a labelling-density titration. The sampling artefact (denser clusters have more members, and a larger sample's median lands later under a decaying detection rate, so a positive trend is partly expected with <i>no</i> crowding) is <b>now corrected for</b> by the permutation null band and its Spearman trend test. The test cannot rule out a third variable driving both axes; that still needs a control experiment.")),
          helpText("One point per Median Cluster: does a denser structure get detected later?"),
          .detail("Crowding = median k-NN distance among the cluster's de-blinked molecule ",
                  "centroids (trace-ID centroids if de-blinking is off; smaller = more crowded), ",
                  "± MAD. Detection-time shift = the cluster's median tim minus the unclustered ",
                  "background median in the same file, so positive means detected later. If ",
                  "crowding delays detection, dense clusters show a positive shift. ",
                  tags$b("This is a density-correlated timing effect"), "; crowding is the likely ",
                  "but not the only cause. Requires Median Clustering and a tim column."),
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
              "The band is what you would see with ", tags$b("zero"), " crowding effect. Read the ",
              "points against it, and the Spearman p for the trend."),
            .detail("Shuffles detection times among traces within each file, keeping every ",
              "cluster's ", tags$b("membership, size and geometry fixed"), ". That reproduces the ",
              "part of the time-shift coming purely from cluster size (a larger sample's median ",
              "lands later under a decaying detection rate) while destroying any real ",
              "crowding-timing link. The Spearman p tests the trend, so it is insensitive to a ",
              "constant clustered-vs-background offset.")
          ),
          actionButton("run_crowd", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_crowd", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("crowd_plot", height = "420px")
        )
      ),

      tags$hr(),

      # -- Within-Cluster2 NN vs Cluster Size --------------------------------
      fluidRow(
        column(3,
          h4("Within Median Clusters NN vs Cluster Size", .exploratory_badge()),
          .method_box("orange",
            HTML("<b>MINFLUXR original.</b> A scaling diagnostic assembled from standard parts: k-nearest-neighbour spacing (<b>Clark &amp; Evans 1954</b>) inside each cluster, plotted against cluster size measured either as member count or as convex-hull area/volume (<b>quickhull</b>, Barber, Dobkin &amp; Huhdanpaa 1996)."),
            HTML("Spacing is computed on <b>de-blinked centroids</b> (trace medians, or de-blinked-molecule medians when de-blinking is on), so it measures molecule-to-molecule spacing inside a structure, which MINFLUX's nm precision supports. Both quantities are per-structure, so each cluster contributes one point regardless of how many localizations it collected."),
            HTML("<b>The slope is confounded by the clustering algorithm.</b> Under DBSCAN every within-cluster NN distance is bounded above by the link radius eps2 (default 200&nbsp;nm), so the curve is <i>censored</i> and cannot rise past it. A flat trend at high size may be that ceiling rather than constant spacing. HDBSCAN and OPTICS impose no such bound, so re-run under those before reading a slope. Reading the slope as 'ballooning' or 'densifying' is this tool's interpretation rather than a validated model, and convex-hull size is upward-biased by single outlying members.")),
          helpText("Does a bigger structure pack its members tighter, looser, or the same?"),
          .detail("One point per Median Cluster across all files. Y = median over members of each ",
                  "member's mean distance to its k nearest member centroids (k = 1 is ",
                  "nearest-neighbour spacing). Flat = constant spacing, growing at fixed density; ",
                  "rising = the structure balloons; falling = it densifies. ",
                  tags$b("Under DBSCAN within-cluster NN is bounded by eps2"), " (the link radius, ",
                  "default 200 nm); HDBSCAN and OPTICS impose no such bound. Requires Median Clustering."),
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
          plotlyOutput("size_nn_plot", height = "420px"),
          uiOutput("size_nn_small_n")
        )
      ),

      tags$hr(),

      # -- Cluster Viewer (2D / 3D + shape fit) ------------------------------
      fluidRow(
        column(3,
          h4("Cluster Viewer and Shape Fitting", .exploratory_badge()),
          .method_box("orange",
            HTML("Machinery is standard, dictionary is bespoke. Fits are <b>orthogonal-distance (total least squares) MLE</b> (minimizing the squared point-to-surface distance), ranked by <b>AIC</b> (Akaike 1974), its small-sample correction <b>AICc</b> (Hurvich &amp; Tsai 1989) or <b>BIC</b> (Schwarz 1978). The 2D view is a PCA projection onto PC1×PC2 (Pearson 1901). Under an iid isotropic-Gaussian residual model σ̂²&nbsp;=&nbsp;RMSD², so −2logL&nbsp;=&nbsp;n·log(2πσ̂²)+n and AICc&nbsp;=&nbsp;AIC&nbsp;+&nbsp;2k(k+1)/(n−k−1)."),
            HTML("Fits <b>de-blinked trace centroids</b>, so that the many localizations of one long trace do not dominate the residual sum and pull the fit toward the brightest molecule. MINFLUX's nm precision is what makes a ~50–200&nbsp;nm shape distinguishable at all. AICc → ∞ when n&nbsp;≤&nbsp;k+1, so the ranking blanks out instead of over-fitting a 5-point cluster."),
            HTML("<b>The shape dictionary is imposed by this tool</b> (sphere / ellipsoid / helix / ring / arc), and 'best' is a model-selection outcome rather than proof of biological geometry. The true shape may be absent from the list, and the ranking will still name a winner. With few, blinking-inflated points several shapes are statistically indistinguishable, and RMSD always favours the most flexible model (prefer AICc). The optional photon (<code>eco</code>) weighting above addresses the unequal-precision problem; with it off, every trace is treated as equally precise.")),
          helpText("Step through individual Median Cluster structures from every loaded file."),
          .detail("Run loads the clusters from all loaded datasets, using the same de-blinking ",
                  "and Median Clustering settings as everything else, so an active file with no ",
                  "clusters does not block the viewer. Structures are keyed ",
                  "by file and id, since cluster ids are only unique within a file.",
                  tags$br(), tags$br(),
                  tags$b("The plotted points, and the shape fit, follow de-blinking."), " With ",
                  "de-blinking on you get one point per de-blinked molecule; with it off, one per ",
                  "trace ID. This matters for the fit: one fluorophore reappears as several trace ",
                  "IDs, so on trace IDs a molecule that blinked six times pulls the surface six ",
                  "times as hard as one that blinked once, and the fitted shape tracks blinking as ",
                  "much as geometry. The header and the legend name whichever unit is in use.",
                  tags$br(), tags$br(),
                  "View raw 3D, or a PCA-normalized 2D projection (PC1 vs PC2, the most ",
                  "differential axes; PCA runs in 3D and reduces to 2D when z exists). Shape ",
                  "fitting is optional; enable it below. Requires Median Clustering."),
          checkboxInput("viewer_active_only", "Active file only", value = FALSE),
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
              "Weights each trace by its photon count. Off treats a 50-photon trace as being as ",
              "trustworthy as a 5,000-photon one."),
            .detail("MINFLUX localization variance scales as 1/photons, so a trace's centroid has ",
              "variance ∝ 1/Σeco, so its total photon count ", tags$b("is"), " its inverse-variance ",
              "weight. The header reports ", tags$b("n_eff"), " (Kish effective sample size): if it ",
              "falls far below the trace count, a few very bright traces are carrying the fit. ",
              "Needs an ", tags$code("eco"), " column and falls back to unweighted if absent. ",
              "Criterion values are not comparable between weighted and unweighted mode."),
            tags$div(style = "margin-top:2px; padding:5px 8px; border-left:3px solid #b26a00; background:rgba(178,106,0,0.06); font-size:0.8em; color:#7a5c00; line-height:1.45;",
              HTML(paste0(
                "<b>De-blink and filter double events first.</b> Fits use the median coordinates of each trace ID ",
                "(or de-blinked molecule), so trace length does not pull the fit. Several trace IDs stacked on the same ",
                "spot still pull it, so de-blink to avoid skewing the fit. Filtering double events with the efo threshold ",
                "is also recommended, since the fit assumes a cluster of distinct de-blinked molecules."))),
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

      # -- Radius of Gyration ------------------------------------------------
      fluidRow(
        column(3,
          h4("Radius of Gyration (RoG)"),
          .method_box("green",
            HTML("Second moment of a point cloud about its centroid, Rg²&nbsp;=&nbsp;⟨‖p&nbsp;−&nbsp;p̄‖²⟩, an RMS size. This tool answers <b>two</b> different questions depending on the <i>Measure</i> setting: the <b>size of your structures</b> (grouping by median cluster) or the <b>footprint of your labelling</b> (spread of one molecule's own repeat blinks). The gyration-tensor eigenvalues also give the shape descriptors used elsewhere (asphericity, <b>Rudnick &amp; Gaspari 1986</b>; linear/planar/spherical coefficients, <b>Westin et&nbsp;al. 1997/2002</b>)."),
            HTML("Points are <b>median positions</b>, never raw localizations, so localization noise is suppressed (a median of N localizations has SD ~ σ_loc/√N) and Rg reflects real geometry. For <b>structure size</b> with de-blinking on, each de-blinked molecule contributes <b>one</b> point, so a multivalent label's several traces do not inflate the size. For the <b>label footprint</b> measure, each of a molecule's traces is a point, so Rg is exactly the blink/label scatter (~linker length for a tag, ~10–40&nbsp;nm for an antibody)."),
            HTML("Rg is biased low at small member counts, because the centroid is estimated from the same points whose spread is measured. The bias grows as members are lost (29% at 2 members, 11% at 5, 1% at 50) and so is <b>member-count dependent</b>, which would create an artificial trend of Rg with cluster size. A Bessel <b>n/(n−1)</b> factor on the mean squared distance corrects it; at 2 members a residual underestimate remains, so treat those clusters as indicative. Rg is also a single number: a filament and a disc of the same extent share an Rg, which is what the anisotropy coefficients are for. The structure-size measure carries the label footprint as an additive floor (no footprint deconvolution is applied); compare against the label-footprint measure to gauge how much of a small structure's Rg is really labelling.")),
          helpText("Distribution of per-cluster Rg across all loaded files. Two views:",
                   "structure size (spread of members in the median cluster) or the per-molecule",
                   "labelling footprint."),
          radioButtons("rog_cluster_col", "Measure",
            choices = c(
              "Median-cluster size (structure)"                = "cluster2",
              "Label footprint (within-molecule blink spread)" = "cluster"
            ),
            selected = "cluster2"
          ),
          radioButtons("rog_yaxis", "Y axis",
                       choices = c("Count" = "count", "Density" = "density"),
                       selected = "count", inline = TRUE),
          numericInput("rog_bins", "Bins", value = 30, min = 5, max = 200),
          checkboxInput("rog_aggregate", "Aggregate all datasets into one", value = FALSE),
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("rog_by_channel", "One distribution per channel", value = FALSE)
          ),
          actionButton("run_rog", "Run",
                       icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_rog", "Download", class = "btn-sm")
        ),
        column(9,
          plotlyOutput("rog_hist", height = "380px"),
          uiOutput("rog_small_n")
        )
      ),

      tags$hr(),

      # -- Cross-Channel Analysis --------------------------------------------
      # Everything in this section measures a relationship BETWEEN two channels.
      # That is what separates it from the panels above, which each describe one
      # distribution and gained a per-channel facet instead. These tools are
      # pairwise by mathematics - L12, g12 and cross-NN each take exactly two
      # point patterns - so with more than two channels a pair has to be chosen.
      fluidRow(
        column(12, h3("Cross-Channel Analysis"))
      ),

      fluidRow(
        column(12,
          tags$div(
            style = "margin:0 0 10px; padding:7px 12px; border-left:4px solid #cfe0ee; background:#f6f9fc; font-size:0.86em; color:#445;",
            "These tools measure how two channels sit relative to each other, so ",
            "distances DO cross the channel boundary here. Clustering never does. ",
            "A 2-colour split has exactly two channels and needs no ",
            "choice; with a multiplex build, pick the pair below.")
        )
      ),

      conditionalPanel(
        condition = "output.mplex_active == true",
        fluidRow(
          column(12,
            tags$div(
              style = "background:#eef6ff; border:1px solid #b8daff; border-radius:4px; padding:8px 12px; margin-bottom:10px;",
              tags$b("Channel pair"),
              helpText("Pick two channels to compare; the cross-channel tools below (NN, Ripley L₁₂, PCF g₁₂, KDE) run on this pair instead of the DCR/fluo split. Channel colours come from the Display tab."),
              radioButtons("coloc_null_mode", "Null model for the L\u2081\u2082 and g\u2081\u2082 envelopes",
                           choices = c("Random labelling of the observed points" = "labels",
                                       "Independent random patterns in the window" = "csr"),
                           selected = "labels"),
              helpText(style = "font-size:0.85em; color:#445;",
                HTML("<b>Random labelling</b> keeps every point where it is and shuffles which channel it belongs to, so both channels keep their own shape and the cell outline cancels. It tests whether the two channels are a random labelling of one pattern, so two species that occupy genuinely different parts of the cell will separate on that alone. <b>Independent random patterns</b> scatters both channels uniformly in the window, so it inherits any error in the window outline and reports co-localization too readily for a cell that does not fill its frame.")),
              helpText(style = "color:#7a5c00;",
                HTML("<b>Registration.</b> Multiplex channels come from separate files or time windows, so drift or stage movement between rounds offsets one channel against the other. The cross-channel tools read such an offset as reduced co-localization or as exclusion. Check the overlay on a structure both channels share, or on fiducials, before quoting a cross-channel result. The DCR split shares one coordinate frame and needs no registration.")),
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
            HTML("Cross-type nearest-neighbour distance distribution: for every point of channel A, the distance to its k-th nearest channel-B point (<b>Clark &amp; Evans 1954</b> extended to the bivariate case, cf. <b>Dixon 2002</b>). The second panel is the empirical cumulative distribution of those distances, the fraction of A molecules with a B partner within distance d. This is an object-based co-localization curve in the sense of <b>Lachmanovich et&nbsp;al. 2003</b> / <b>Dunn et&nbsp;al. 2011</b>."),
            HTML("Object-based co-localization suits MINFLUX: there are no pixels and no intensity image, only coordinates, and nm precision makes the distance axis physically meaningful down to the label size. Cross-NN is computed strictly between channels (never within), on de-blinked centroids, and the Median-Cluster modes restrict pairing to molecules inside the same structure so a partner across the cell cannot be counted."),
            HTML("Cross-NN has <b>no null model built in</b>. Every A point has a nearest B point, so the curve is never zero and a small median distance can arise from density alone. Judge it against Ripley&nbsp;L₁₂ or g₁₂ with their CSR envelopes. The dominant two-colour artefact here is <b>DCR crosstalk</b>: the split thresholds an overlapping DCR distribution, and every mis-assigned molecule contributes a near-zero distance that lands in the first histogram bin and the steep part of the co-localization curve. Read the co-localized fraction against the <b>per-group contamination</b> reported by the GMM auto-split, and use <b>pooled DCR</b> when available, since it usually separates the dyes considerably better.")),
          helpText("Every distance is measured strictly between the two channels."),
          .detail("Channel assignment comes from the Display tab; clustering uses the single ",
                  "de-blink / Median Clustering setup, run separately within each channel, so the ",
                  "channels are never pooled or clustered together. The Median Cluster modes keep ",
                  "only molecules inside a Median Cluster of their own channel, dropping the ",
                  "unclustered background from both sides, so the question becomes how far channel-1 ",
                  "structures sit from channel-2 structures, ignoring stray background. ",
                  "'Median Cluster centroids' goes further and collapses each structure to one ",
                  "point, giving one distance per structure."),
          selectInput("cross_nn_unit", "Unit",
            choices = c(
              "TID centroids"                            = "tid",
              "De-blinked molecule centroids"                = "cluster",
              "TID centroids in Median Clusters (own channel)"            = "cluster2_tid",
              "De-blinked molecule centroids in Median Clusters (own channel)" = "cluster2_cluster",
              "Median Cluster centroids (own channel)"                    = "cluster2_centroid"
            ),
            selected = "tid"),
          numericInput("cross_nn_k", "k-th nearest neighbour", value = 1, min = 1, max = 10, step = 1),
          radioButtons("cross_nn_hist_type", "Histogram type",
                       choices = c("Count" = "count", "Density" = "density"),
                       selected = "density", inline = TRUE),
          numericInput("cross_nn_bins", "Bins", value = 40, min = 5, max = 200),
          checkboxInput("cross_nn_median_lines", "Show median lines", value = TRUE),
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
            HTML("Runs on de-blinked centroids of the two channels (DCR/fluorophore split, or a chosen multiplex channel pair) with the ≤2000-point subsample per channel. The <b>CSR envelope randomizes both channels independently in the shared window</b>, which is the correct null for 'these two species are placed without regard to each other' and carries any remaining edge bias. Median-Cluster modes compute L₁₂ inside each structure and average the curves, so a co-localization signal cannot be produced by the two channels merely occupying the same region of the cell. <b>A MINFLUX advantage applies here:</b> the DCR split excites both dyes with the same laser through the same beam path, so both channels share one coordinate frame. There is <b>no chromatic offset and no channel registration to calibrate</b>, and L₁₂ is trustworthy down to the localization precision."),
            HTML("Edge correction and the footprint window are shared with L(r); read significance against the envelope. The real two-colour error is <b>DCR mis-assignment (crosstalk)</b>. The split thresholds a continuous, overlapping DCR distribution, and a wrongly-assigned molecule sits by construction right next to its true neighbours, so crosstalk biases L₁₂ <b>upward at short r, toward false co-localization</b>, exactly where a real interaction shows. <b>Quantify it before trusting a result:</b> the GMM auto-split in the Display tab reports d′, the expected mis-assignment rate and the per-group contamination, and a co-localization signal is only meaningful if it is large compared with those. <b>Pool DCR at upload</b> (photon-weighted across the iteration sequence) if you have not, since it usually separates the two dyes far better than raw DCR, and a better split is worth more here than any statistical correction. Drift can offset the multiplex-by-files and by-time-window modes, which overlay separate acquisitions; this does not affect the DCR split. Asymmetric labelling density also skews the estimator. Finally, nm-scale co-localization measures proximity, and proximity is not interaction.")),
          helpText("L₁₂(r) − r > 0 is co-localization; < 0 is spatial exclusion."),
          .detail("Median Cluster modes compute L₁₂ within each Median Cluster group and average ",
                  "the curves. Large datasets are subsampled to 2,000 points per channel per ",
                  "file or group."),
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

      # -- Pair Cross-Correlation g₁₂(r) ------------------------------------
      fluidRow(
        column(3,
          h4("Pair Cross-Correlation g₁₂(r)"),
          .method_box("green",
            HTML("The pair-correlation (radial distribution) function, the derivative form of Ripley's K, counting cross pairs in an <b>annulus/shell</b> (r−dr/2, r+dr/2] normalized by its area (2D) or volume (3D). g₁₂&nbsp;=&nbsp;1 under CSR, &gt;1 attraction, &lt;1 repulsion. In SMLM this is the PC-PALM statistic (<b>Sengupta et&nbsp;al. 2011</b>; <b>Veatch et&nbsp;al. 2012</b>), developed to separate multiple appearances of one protein from genuine clustering."),
            HTML("g₁₂ is <b>scale-local</b>, so it is the more informative of the two cross statistics on MINFLUX data. It puts the interaction distance at a peak, and a short-r blinking or labelling artefact stays confined to the first few bins. Same de-blinked units, shared window and independent two-channel CSR envelope as L₁₂."),
            HTML("<b>Noisier than L₁₂</b>, since each shell holds few pairs, so short-r bins are unstable and sensitive to the r-step choice (an apparent peak in the first bin is often a binning effect). Edge correction is shared with L(r); read significance against the envelope, which brackets g&nbsp;=&nbsp;1. DCR crosstalk and asymmetric labelling affect g₁₂ exactly as they affect L₁₂, and because g₁₂ is scale-local the crosstalk piles into the <i>first shells</i>, which is also where a genuine interaction peak would appear.")),
          helpText("g₁₂(r) = 1 under CSR; > 1 is attraction at scale r, < 1 is repulsion."),
          .detail("Complements Ripley L₁₂ by showing the interaction scale as a peak."),
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

      # -- 2D KDE Density Co-localization Map -------------------------------
      fluidRow(
        column(3,
          h4("2D KDE Density Map"),
          .method_box("amber",
            HTML("2D kernel density estimation (<b>Silverman 1986</b>; <code>MASS::kde2d</code>) of each channel on a shared grid, with the auto bandwidth from Silverman's rule bw&nbsp;=&nbsp;0.9·min(sd,&nbsp;IQR/1.34)·n<sup>−1/5</sup>. The co-localization index is the <b>Pearson correlation</b> of the two density surfaces (the coordinate-data analogue of Pearson co-localization on images, cf. <b>Dunn et&nbsp;al. 2011</b>), with an optional randomization test for significance."),
            HTML("MINFLUX has no pixels, so the density surface, and therefore the 'image' the Pearson r is computed on, is <b>manufactured by the bandwidth</b>, a free parameter of the analysis. Two tools compensate: the <b>r-vs-bandwidth sweep</b> shows the whole family of answers, and the <b>randomization test</b> gives a p-value at the chosen bandwidth. Its default null is a <b>toroidal random shift</b> of channel 2 (<b>Lotwick &amp; Silverman 1982</b>), which preserves that channel's own clustering and point count and breaks only its registration with channel 1. A uniform-CSR null is selectable for comparison."),
            HTML("<b>The Pearson r is relative.</b> A small bandwidth drives r toward 0 even for perfectly co-localizing molecules (the two surfaces become disjoint spikes), and a large one inflates r toward 1 because both channels simply fill the cell. Only compare r values computed at the same bandwidth on similarly dense data, and prefer the sweep to any single number. <b>Where the null is drawn decides whether the p-value holds.</b> With the sample footprint switched on at the top of the tab, both nulls randomize inside the footprint and the test is conservative: a p-value above 0.05 can be quoted as not significant. With it off, they randomize over the bounding box, and for a structure occupying a small part of the view they reject too readily, the uniform-CSR null severely so. For a positive claim, use g₁₂(r) with its envelope.")),
          helpText("Density map per channel for the active file. ",
                   tags$b("Treat the Pearson r as a relative measure.")),
          .detail("Ch1 is drawn as red contours, Ch2 as blue. Pearson r of the two density ",
                  "surfaces is a global co-localization coefficient, but it is strongly ",
                  "bandwidth-dependent: a small bandwidth gives near-zero r even for perfectly ",
                  "co-localizing molecules, and a large one inflates it. Choose the dataset with ",
                  "'Active file' in the Display tab."),
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
          helpText(HTML("<b>Random shift</b> breaks only the cross-channel registration; <b>Uniform CSR</b> tests against complete randomness.")),
          .detail(HTML("Random shift translates Ch2 with wrap-around, preserving its own clustering and point count (Lotwick &amp; Silverman 1982). Uniform CSR replaces Ch2 with uniform random points.<br><b>Limitation:</b> both nulls randomize over the bounding box rather than the region the signal actually occupies. When the signal fills the field the p-value is well calibrated; when it is confined to part of the field, such as a single cell in a wider view, both reject far too readily, CSR much more so. Treat a low p there as weak evidence and confirm with g₁₂(r) against its envelope.")),
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

      # -- Methods & References (consolidated citation box) ----------------------
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
                    tags$td("k-th NN distance distribution via k-d tree; NN distance as a spatial statistic after Clark & Evans 1954; their observed/expected R index is not computed (no null model here)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("DBSCAN / HDBSCAN / OPTICS")),
                    tags$td("Ester, Kriegel, Sander & Xu 1996 (DBSCAN); Ankerst, Breunig, Kriegel & Sander 1999 (OPTICS); Campello, Moulavi & Sander 2013 (HDBSCAN)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Median (de-blinked) clustering")),
                    tags$td("Two-tier grouping (localizations → trace → de-blinked molecule → median cluster) collapses multi-blink / multivalent overcounting before clustering, an adaptation for MINFLUX/SMLM in the spirit of Annibale et al. 2011; Sengupta et al. 2011 (PC-PALM)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Trace IDs / de-blinked molecules per median cluster")),
                    tags$td("Median-cluster occupancy counting in the SMLM molecular-counting lineage (Annibale et al. 2011; Sengupta et al. 2011); trace IDs (blinking events) or de-blinked molecules, read as occupancy")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Cross-channel NN & co-loc curve")),
                    tags$td("Bivariate NN distances (Clark & Evans 1954; Dixon 2002) with an object-based co-localization read-out (Lachmanovich et al. 2003; Dunn et al. 2011)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Radius of gyration")),
                    tags$td("Standard second-moment cluster size (gyration tensor); anisotropy after Rudnick & Gaspari 1986")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Ripley's L(r)")),
                    tags$td("Ripley 1977; variance-stabilising L-transform, Besag 1977")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Pair correlation g(r)")),
                    tags$td("Pair-correlation / radial distribution function (Illian et al. 2008; Stoyan & Stoyan 1994); the SMLM lineage is PC-PALM, Sengupta et al. 2011 and Veatch et al. 2012; density-invariant by construction (normalized by intensity squared)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Labelling-density normalization")),
                    tags$td(tags$span("EXPLORATORY", style="color:#e67e22;font-weight:700;"),
                            ". Independent p-thinning, whose exact invariance of g(r) and K(r) is standard point-process theory (Illian et al. 2008; Chiu et al. 2013); morphological support estimate after Serra 1982; support-free k-NN intensity after Thompson 1956. Applying it as a between-condition correction, and the choice of thinning target, are this tool's design and validated in silico only")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Bivariate Ripley L₁₂(r)")),
                    tags$td("Lotwick & Silverman 1982; Dixon 2002 (bivariate/cross point-pattern statistics)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Pair cross-correlation g₁₂(r)")),
                    tags$td("Sengupta et al. 2011 (PC-PALM); Veatch et al. 2012 (correlation functions / over-counting)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("KDE density + Pearson co-loc")),
                    tags$td("Kernel density estimate (Silverman 1986); Pearson coefficient as a bandwidth-dependent, relative co-localization measure (cf. Dunn et al. 2011); significance by structure-preserving random-shift randomization of Ch2 (Lotwick & Silverman 1982), with a uniform-CSR null selectable for comparison")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Shape fitting (viewer)")),
                    tags$td(tags$span("EXPLORATORY", style="color:#e67e22;font-weight:700;"),
                            ". Orthogonal-distance MLE with AICc/BIC model selection (Akaike 1974; Schwarz 1978; Hurvich & Tsai 1989); the shape dictionary and biological reading are bespoke")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Ground-truth estimate & structure size")),
                    tags$td(tags$span("EXPLORATORY", style="color:#e67e22;font-weight:700;"),
                            ". Labelling-model inverse (this tool), validated in silico against the simulator; contour via MST (Prim 1957), anisotropy coefficients (Westin et al. 1997/2002)")),
            tags$tr(tags$td(style="padding:1px 14px 1px 0; vertical-align:top; color:#666;", tags$b("Detection-time / crowding / NN-vs-size")),
                    tags$td(tags$span("EXPLORATORY", style="color:#e67e22;font-weight:700;"),
                            ". Descriptive diagnostics of this tool with no published MINFLUX equivalent; correlational only"))
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
