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
                 helpText("One set of bins and one count axis for all three, for a figure. Off by default because ",
                          "the axes differ in precision."),
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


        helpText("Whole traces merge into molecules and are never split. De-blinking only merges; removal ",
                 "is done by the metadata and spread filters.")
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
              helpText("Wrap: 3\u00d7 the localization precision, at least half the molecule spacing; set per axis,",
                       "so a more precise z is wrapped more tightly."),
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
                           helpText("Assumes an arc: PC1 is the chord half-length and PC2 the bow height, R = (chord\u00b2 + sagitta\u00b2)/(2\u00d7sagitta). Not meaningful for full rings, spheres or helices.")),
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
          helpText("Heuristic: analytic geometry from PCA axes, no optimisation. MLE: minimises squared distances to the surface. RMSD is in nm for both."),
          conditionalPanel(
            condition = "input.fit_method == 'mle'",
            checkboxInput("fit_weighted", "Weight by per-trace precision (n/spread\u00b2, \u2248 1/SEM\u00b2)", FALSE),
            helpText("Each point is weighted by n\u2009/\u2009spread\u00b2, the inverse squared standard error of its median. Single-localization traces get the median spread."),
            radioButtons(
              "mle_center_mode", "MLE centre handling",
              choices = c(
                "Regularized: MAP prior toward optimized centre"  = "reg",
                "Fixed to optimized centre (analytical radius)"    = "fixed",
                "Free: bounded to data bounding box"           = "free"
              ),
              selected = "reg"
            ),
            helpText("Fixed: hull centroid, analytic radius. Regularized: free centre, penalised for moving (strength \u03bb). Free: centre fully optimised."),
            conditionalPanel(
              condition = "input.mle_center_mode == 'reg'",
              sliderInput("mle_center_lambda",
                          "Regularization strength \u03bb",
                          min = 0, max = 5, value = 0.1, step = 0.05),
              helpText("\u03bb = 0.1: a 5 nm centre shift adds about 2.5 nm\u00b2 to the mean squared residual. Independent of cluster size.")
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
            helpText("Seeds the pitch estimate only; the fitted turn count follows the data.")
          ),
          conditionalPanel(
            condition = "input.fit_type == 'arc'",
            sliderInput("arc_span_deg", "Arc span hint (\u00b0, 0 = auto)", min = 0, max = 360, value = 270, step = 10),
            helpText("Starting arc angle for the fit. 0 = from the data\u2019s angular spread.")
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
              "Posterior probability that each trace belongs to its assigned group: 1.0 unambiguous, 0.5 a ",
              "coin flip. ", tags$b("Orange = on the edge"), ": these would switch channel under a small shift ",
              "of the threshold and create false short-range cross-channel pairs. Export them to inspect."),
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
      helpText("Build an N-channel dataset: overlay loaded files, or split one acquisition by time (TIM). Each channel gets its own colour."),
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
        helpText("Exploratory. Inverts the labelling model to estimate ground-truth positions, an approximate count and a crowding flag. Circles run green (certain) to red (uncertain)."),
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
              helpText(HTML("The count divides trace IDs by the expected traces per molecule, from the preset or measured on a <b>counting standard</b>: a sparse sample of known stoichiometry, same label and acquisition.")),
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
          helpText("Auto is recommended: de-blinking brings each molecule close to one unit, so the count divides by about 1."),
          numericInput("gt_complex_r",
                       "Label-complex merge radius (nm, blank = auto)",
                       value = NA, min = 0, step = 5),
          helpText("Units within this radius merge into one candidate position. Blank: set by the labelling preset; 0: no merge."),
          checkboxInput("gt_guided_merge",
                        "Stoichiometry-guided merge (model-driven)", value = FALSE),
          helpText("Where the spatial merge leaves more positions than the stoichiometry predicts, merges further toward that target; never splits. The result is prior-driven: check it against the two counts below."),
          checkboxInput("gt_displace_inward",
                        "Displace GT inward to structure (linker correction)", value = TRUE),
          numericInput("gt_displace_nm",
                       "Inward displacement (nm, blank = auto = linker length)",
                       value = NA, min = 0, step = 1),
          helpText("A linker places the dye outside the structure; this moves each position back toward the local structure centre. Needs Median Clusters."),
          checkboxInput("gt_show_lines", "Draw contribution lines (unit → GT)", value = TRUE),
          sliderInput("gt_circle_size", "GT marker size", min = 4, max = 40, value = 15, step = 1),
          sliderInput("gt_support_strength",
                      "Support-penalty strength (position certainty)",
                      min = 0, max = 2, value = 1, step = 0.1),
          helpText("How much a position backed by few traces is penalised in the uncertainty. 0 = off; lower it if sparse labelling reads too red."),
          uiOutput("gt_estimate_readout"),
          tags$hr(style = "margin:8px 0;"),
          tags$b("Structure shape & size"),
          checkboxInput("gt_show_shape",
                        "Classify shape + measure size (contour, extent, spacing)", value = FALSE),
          conditionalPanel(
            condition = "input.gt_show_shape == true",
            helpText("Classifies each cluster's shape (filament, ring, punctum, sheet, diffuse) and reports contour length, extent, width, Rg and point spacing."),
            numericInput("gt_subunit_spacing",
                         "Assumed subunit spacing (nm, optional - blank = off)",
                         value = NA, min = 0.1, step = 0.5),
            helpText("If the structure has a known periodic spacing, a filament's contour length is divided by it to estimate a subunit count."),
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
                HTML("<b>Independent random thinning</b>: each unit is kept with a fixed probability, which leaves g(r) and Ripley's K unchanged in expectation (<b>Illian et&nbsp;al. 2008</b>; <b>Chiu et&nbsp;al. 2013</b>). Support from a smoothed-density region (<b>Hyndman 1996</b>) or morphological occupancy (<b>Serra 1982</b>)."),
                HTML("Whole <b>units</b> are dropped, never single localizations, so traces and de-blinking stay intact. Density divides by the estimated support, since a field includes the cell edge and dark regions. The draw is seeded per file."),
                HTML("Use it only when density varies for reasons outside the biology. When density is the biology, equalising it removes the effect. Thinning is one-way and can push small clusters below minPts.")),
              helpText("For cells whose labelling density varies for reasons outside the biology. Leave off otherwise."),
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
                helpText("The area density is divided by. The bounding box counts the cell edge and dark regions as sample."),
                radioButtons("dens_support_from", "Outline drawn from",
                             choices = c("Trace / molecule centroids" = "units",
                                         "All localizations"          = "locs"),
                             selected = "units", inline = TRUE),
                helpText("Which points define the cell edge."),
                .detail(HTML(
                  "<b>Centroids</b> are what every analysis on this tab counts, but a sparsely labelled cell gives a smaller, raggeder outline, so the area carries some of the unit count.<br><br>",
                  "<b>All localizations</b> give a smoother outline, but a brighter or longer-lived label enlarges it on its own, which loosens the calibrated eps.")),
                numericInput("dens_grid_h", "Grid cell (nm)", value = 100, min = 10, step = 10),
                conditionalPanel(
                  condition = "input.dens_support_mode == 'kde'",
                  numericInput("dens_kde_bw", "Smoothing bandwidth (nm)", value = 150, min = 10, step = 25),
                  numericInput("dens_kde_frac", "Fraction of units enclosed", value = 0.99,
                               min = 0.5, max = 1, step = 0.01),
                  helpText(HTML("The smallest region holding this share of units. <b>Keep fraction and bandwidth identical across compared samples.</b>"))
                ),
                conditionalPanel(
                  condition = "input.dens_support_mode == 'grid'",
                  checkboxInput("dens_close_auto", "Auto gap fill", value = TRUE),
                  helpText(HTML("<b>The key setting.</b> Too small and the area traces your structures instead of the cell.")),
              .detail(HTML("It fills gaps from sparse sampling inside the cell, and keeps out dark regions and the outside, since it only closes holes narrower than about twice its radius. Auto picks a radius per file and applies one common value to all.")),
                  conditionalPanel(
                    condition = "input.dens_close_auto == false",
                    numericInput("dens_close_nm", "Fill gaps up to (nm)", value = 200, min = 0, step = 50)
                  ),
                  checkboxInput("dens_show_sweep", "Show gap-fill sweep", value = FALSE)
                ),
                checkboxInput("dens_show_map", "Show support map", value = FALSE),
                helpText(HTML("Counted area <b style=\"color:#1a8a1a\">green</b>, excluded <b style=\"color:#c03030\">red</b>. The green region should look like your cell.")),
                conditionalPanel(
                  condition = "input.dens_show_map == true",
                  selectInput("dens_map_file", "Dataset to render", choices = NULL)
                ),
                checkboxInput("dens_plot_log", "Log density axis", value = TRUE),
                helpText("Local density spans a wide range within one cell."),
                tags$hr(style = "margin:8px 0;"),
                checkboxInput("dens_thin_enable", "Thin all files to a common density", value = FALSE),
            tags$div(style = "font-size:0.82em; color:#7a5c00; margin:2px 0 6px; padding:4px 8px; border-left:3px solid #e0c890; background:#fffdf6;",
              tags$b("Only for count-based read-outs"), " (DBSCAN counts, % in clusters, cluster sizes, NN ",
              "distances). Ripley's L and g(r) already divide density out."),
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
          helpText("Assign files to named groups; the analyses below then show one series per group."),
          helpText(style = "color:#7a5c00;",
                   "Pooled units from one cell are correlated, so compare groups on one value per file ",
                   "from the per-file summary."),
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
            HTML("Every localization relative to its own trace centroid, pooled at (0,0). The RMS radius is &sigma;<sub>loc</sub>, drawn as 1&sigma; and 2&sigma; ellipses with a &sigma; per axis."),
            HTML("Runs on the filtered localizations. The ellipses are rotated to the cloud's own axes, so a tilted cloud is not reported as round; a tilt shared across traces points to directional drift or the instrument."),
            HTML("&sigma; is an <b>upper bound on precision</b>: it also holds label, linker and molecule motion and any drift within a trace. In 3D, &sigma;<sub>z</sub> and &sigma;<sub>xy</sub> can differ in either direction.")),
          .method_box("amber",
            HTML("<b>Sample drift and jitter</b> from the displacement of re-blinks of one molecule, regressed on the time between them."),
            HTML("Within each de-blinked group the first trace is the origin. <b>Drift</b> is the part that grows with time; <b>jitter</b> is the scatter left over. Two molecules merged by mistake widen the interval without tilting the slope. De-blink first; otherwise the tool groups by proximity."),
            HTML("Blinks further apart than the grouping radius stop pairing, so a large drift is reported as a lower bound. Fiducial-locked acquisitions are already drift-corrected.")),
          helpText("Per-axis \u03c3, anisotropy, skew and kurtosis of the filtered data, for the active dataset."),
          numericInput("prec_cloud_min_loc", "Min localizations / trace", value = 5, min = 3, step = 1),
          checkboxInput("prec_cloud_use_z", "Include z (report σ_z)", value = FALSE),
          actionButton("run_prec_cloud", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          downloadButton("download_prec_cloud", "Download σ", class = "btn-sm"),
          checkboxInput("view_drift_plot", "Show sample drift plot", value = FALSE),
          conditionalPanel(
            condition = "input.view_drift_plot == true",
            helpText("Each re-blinking molecule starts at the origin; each later blink is one point."),
            radioButtons("drift_plot_color", "Colour points by", inline = TRUE,
                         choices = c("Time since first blink" = "time",
                                     "Acquisition time" = "abs",
                                     "De-blinked ID" = "mol"),
                         selected = "time"),
            conditionalPanel(
              condition = "input.drift_plot_color == 'abs'",
              helpText(style = "font-size:0.78em; color:#888;",
                       "Blue early, red late. Direction changing with colour means the drift changed course.")
            ),
            conditionalPanel(
              condition = "input.drift_plot_color == 'mol'",
              helpText(style = "font-size:0.78em; color:#888;",
                       "Blinks of one molecule share a colour and are joined in time order. Colours repeat.")
            ),
            helpText(style = "font-size:0.78em; color:#888;",
                     "Centred on the origin: the sample held still. Moving away with longer gaps: drift."),
            .detail("Heavy markers: the geometric median displacement per time bin, where jitter cancels. ",
                    "It resists the occasional merged molecule and does not depend on axis labels.")
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
            HTML("Tool-original. Flags traces that behave like fiducials or bright aggregates. <b>Lifetime</b> first: a fiducial does not bleach, so its trace outlives every real fluorophore in the file. <b>Brightness</b> second, for traces that also persisted."),
            HTML("The lifetime threshold is a multiple of the file's own 99th-percentile trace duration, so it adapts to the acquisition. Brightness is a per-file robust z of log efo. Length, recurrence and spread are reported as evidence."),
            HTML("Not validated by other groups. It reports only, and removes nothing until you tick the box below.")),
          helpText("Flags traces that outlive the rest of the file, and bright traces that persist."),
          numericInput("fid_dur_x",    "Lifetime, as a multiple of the file's 99th percentile",
                       value = 5, min = 1.5, max = 50, step = 0.5),
          numericInput("fid_efo_z",    "Brightness outlier (robust z)", value = 4, min = 2, max = 10, step = 0.5),
          numericInput("fid_min_locs", "Min localizations for the brightness rule", value = 5, min = 1, max = 100, step = 1),
          .detail("A multiple below about 5 starts to catch ordinary long traces. The brightness rule ",
                  "needs persistence because a single very bright localization is a blip, not a bead."),
          numericInput("fid_len_z",    "Length outlier (reported only)",     value = 4, min = 2, max = 10, step = 0.5),
          numericInput("fid_recur_r",  "Recurrence radius (nm, reported only)", value = 100, min = 10, max = 1000, step = 10),
          numericInput("fid_recur_min","Min co-located traces (reported only)", value = 4, min = 1, max = 50, step = 1),
          checkboxInput("fid_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          checkboxInput("fid_exclude", "Exclude flagged traces from analyses", value = FALSE),
          helpText(style = "margin-top:-6px; font-size:0.78em; color:#888;",
                   "Reversible. Re-run the scan after changing thresholds."),
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
            HTML("Unique <b>trace IDs</b> per acquisition-time bin: the detection-event rate over the run."),
            HTML("Each trace counts once, at its median time, so long traces do not crowd the early bins. One series per file, or per group; each line ends where its file's time ends."),
            HTML("Descriptive QC. A falling rate reflects bleaching, exhausted blinking, drift and depletion.")),
          tags$b("Trace-ID length prefilter", style = "font-size:0.9em;"),
          fluidRow(
            column(6, numericInput("tid_time_len_min", "Min locs/trace", value = 1, min = 1, step = 1)),
            column(6, numericInput("tid_time_len_max", "Max (0 = none)", value = 0, min = 0, step = 1))
          ),
          helpText("Applied identically to every dataset before counting."),
          numericInput("tid_time_bins", "Time bins", value = 40, min = 5, max = 300),
          radioButtons("tid_time_yaxis", "Y axis",
                       choices = c("Count per bin"         = "count",
                                   "Cumulative"            = "cumulative",
                                   "Normalized rate"       = "norm",
                                   "Normalized cumulative" = "norm_cumulative"),
                       selected = "count"),
          helpText("Normalized modes compare timing (shape) across files of different length."),
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
            HTML("Distribution of the k-th nearest-neighbour distance (<code>dbscan::kNN</code>), after <b>Clark &amp; Evans 1954</b>."),
            HTML("Runs on trace or de-blinked-molecule centroids rather than raw localizations, so re-blinks and multivalent labels do not over-count molecules and the distribution approximates the spacing between molecules."),
            HTML("A sharp peak below ~40&nbsp;nm usually reflects the label footprint (Annibale et&nbsp;al. 2011; Sengupta et&nbsp;al. 2011). There is no random reference: for significance use Ripley L(r) with its envelope.")),
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
          helpText("Every loaded file, with the Display tab's filters and clustering. Auto follows the de-blinking switch."),
          numericInput("nn_batch_k", "Nearest neighbors (k)", value = 1, min = 1, max = 20),
          uiOutput("nn_batch_group_ui"),
          shinyWidgets::numericRangeInput(
            "nn_batch_dist_range", "Distance range filter (nm)",
            value = c(0, 2000), min = 0
          ),
          helpText("Distances outside it are left out of the histogram and statistics."),
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
            HTML("Ripley's K (<b>Ripley 1977</b>) with the L-transform (<b>Besag 1977</b>), plotted as L(r)&nbsp;\u2212&nbsp;r: 0 under complete spatial randomness, &gt;0 clustering, &lt;0 dispersion. SMLM lineage: <b>Sengupta et&nbsp;al. 2011</b>; <b>Veatch et&nbsp;al. 2012</b>."),
            HTML("Runs on de-blinked units, so repeated appearances of one molecule do not create a short-range peak. Up to 2000 points per file. The Monte-Carlo envelope simulates random points in the same window, so edge effects appear in the band too."),
            HTML("The shaded band is pointwise; read significance from the dashed <b>global band</b> and the p-value above the plot. Labelling is the main confound: incomplete labelling flattens clustering, multivalent labels create it. In group mode per-file curves are averaged, never pooled.")),
          helpText("L(r) \u2212 r = 0 under complete spatial randomness; > 0 is clustering, < 0 is regularity."),
          selectInput("ripley_unit", "Unit",
            choices = c(
              "Auto (molecules if de-blinking is on)" = "auto",
              "TID centroids"             = "tid",
              "De-blinked molecule centroids" = "cluster"
            ),
            selected = "auto"
          ),
          helpText("Auto follows the Display tab's de-blinking switch."),
          checkboxInput("ripley_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          helpText("One curve per loaded file, each in its own window. Uncheck for the active file only."),
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("ripley_by_channel", "One curve per channel", value = FALSE),
            helpText("One curve per channel, both measured in the same observation window."),
            .detail("The window comes from the whole file, so a difference between the curves reflects ",
                    "arrangement rather than one channel having fewer localizations.")
          ),
          checkboxInput("ripley_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          helpText("Spheres instead of circles. Only for genuinely 3-D data."),
          numericInput("ripley_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          helpText("Largest radius evaluated. 0 = a quarter of the shortest window side."),
          numericInput("ripley_r_steps", "r steps", value = 60, min = 10, max = 300),
          checkboxInput("ripley_smooth", "Smooth curve (loess)", value = TRUE),
          numericInput("ripley_smooth_span", "Smoothing span", value = 0.2, min = 0.05, max = 0.8, step = 0.05),
          helpText("Display only; a large span can flatten real peaks."),
          conditionalPanel(
            condition = "input.analysis_by == 'group'",
            radioButtons("ripley_group_band", "Group spread band",
                         choices = c("IQR (Q25-Q75)" = "iqr", "±SD" = "sd", "None" = "none"),
                         selected = "iqr", inline = TRUE),
            checkboxInput("ripley_show_files", "Show individual files (shaded by group)", value = FALSE),
            helpText("One mean curve per group; the band is the spread between files, not a test."),
            .detail("Each file's L(r)−r is computed on its own window, then averaged within the ",
                    "group; coordinates are never pooled. Curves are interpolated onto a common ",
                    "r-grid, up to the smallest per-file r max.")
          ),
          checkboxInput("ripley_envelope", "Show CSR envelope (Monte-Carlo)", value = FALSE),
          helpText(HTML("Read significance from the dashed <b>global band</b> and its p-value; the shaded band is pointwise.")),
          .detail("In all-files mode the envelope covers the first six files. In group mode it is a ",
                  "dotted outline (the mean of the files' envelopes), separate from the solid ",
                  "between-file band. The smallest p is 1/(simulations + 1)."),
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
            HTML("The pair-correlation function, the derivative form of Ripley's K: pairs in a shell around r, divided by the shell's area (2D) or volume (3D). g = 1 for a random pattern (Illian et&nbsp;al. 2008)."),
            HTML("Runs on the same de-blinked units as Ripley's L. g divides out the intensity, so it compares cells of different labelling density, and it is local in r, so a short-range labelling artefact stays in the first bins."),
            HTML("Short-r bins hold few pairs and are noisy, so a first-bin peak is often binning. Read significance from the global band and its p-value.")),
          helpText(HTML("The <b>r where g returns to 1</b> is the aggregate size; the height is how many times denser than random.")),
          .detail(HTML("A solid aggregate makes g fall monotonically from the smallest r. A peak further out, with g low at short r, means short-range exclusion.")),
          selectInput("pcf_unit", "Unit",
            choices = c("Auto (molecules if de-blinking is on)" = "auto",
                        "TID centroids" = "tid",
                        "De-blinked molecule centroids" = "cluster"),
            selected = "auto"),
          checkboxInput("pcf_all_files", "Use all loaded files (uncheck for active file only)", value = TRUE),
          conditionalPanel(
            condition = "output.n_channels > 1",
            checkboxInput("pcf_by_channel", "One curve per channel", value = FALSE),
            helpText("Each channel in the same window, the whole file.")
          ),
          checkboxInput("pcf_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          helpText("Spherical shells instead of annuli. Only for genuinely 3-D data."),
          numericInput("pcf_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          numericInput("pcf_r_steps", "r steps", value = 50, min = 10, max = 300),
          helpText("Fewer steps: wider shells, less noise, coarser resolution. 40\u201360 usually suits."),
          checkboxInput("pcf_smooth", "Smooth curve (loess)", value = TRUE),
          numericInput("pcf_smooth_span", "Smoothing span", value = 0.25, min = 0.05, max = 0.8, step = 0.05),
          helpText("Display only."),
          fluidRow(
            column(6, numericInput("pcf_win_lo", "Summary window lo (nm)", value = 0, min = 0, step = 10)),
            column(6, numericInput("pcf_win_hi", "hi (nm)", value = 0, min = 0, step = 10))
          ),
          helpText("Range for the summary mean and integrated excess. 0 and 0: 5\u201340% of r max."),
          conditionalPanel(
            condition = "input.analysis_by == 'group'",
            radioButtons("pcf_group_band", "Group spread band",
                         choices = c("IQR (Q25-Q75)" = "iqr", "±SD" = "sd", "None" = "none"),
                         selected = "iqr", inline = TRUE),
            checkboxInput("pcf_show_files", "Show individual files (shaded by group)", value = FALSE)
          ),
          checkboxInput("pcf_envelope", "Show CSR envelope (Monte-Carlo)", value = FALSE),
          helpText(HTML("Nothing is subtracted from your curve: significance is the curve leaving the <b>global band</b>. The band is widest at short r, where rings hold few pairs.")),
          conditionalPanel(
            condition = "input.pcf_envelope == true",
            numericInput("pcf_nsim", "Envelope simulations", value = 99, min = 19, max = 199)
          ),
          checkboxInput("pcf_fit_xi", "Fit correlation length ξ", value = FALSE),
          conditionalPanel(
            condition = "input.pcf_fit_xi == true",
            numericInput("pcf_fit_rmin", "Fit from r (nm)", value = 40, min = 0, step = 5),
            helpText("Fits g(r) \u2212 1 = A\u00b7exp(\u2212r/\u03be), one number per curve. Start above the ",
                     "labelling floor, where the label no longer dominates.")
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
            HTML("k-NN local density (<b>Loftsgaarden &amp; Quesenberry 1965</b>), the quadrat index of dispersion (<b>Fisher, Thornton &amp; Mackenzie 1922</b>) and DBSCAN foci (<b>Ester et&nbsp;al. 1996</b>), each relative to the cell's own mean density. The per-cell eps calibration is tool-original."),
            HTML("Dividing by the cell's own density puts a high and a low expresser on the same scale, and eps is set per cell against a randomised pattern in that cell's footprint."),
            HTML("These read the <b>upper tail</b>: whether a minority of units sit in dense spots, not what those spots are. Focus count and fraction in foci stay density-sensitive; compare cells on <b>Largest focus (\u00d7 random 95th pct)</b>. Quote the quadrat size with the dispersion index.")),
          helpText("Each unit's local density as a multiple of its cell's mean, and the dense spots. Density unevenness is 1.0 for a random pattern."),
          .detail(HTML(
            "<b>Enhancement:</b> local density from the k-th neighbour over the cell mean; 20\u00d7 is twenty times denser than average.<br>",
            "<b>Density unevenness:</b> SD over mean of the local densities, divided by the random value 1/\u221a(k\u22122), so 1.0 is random at any density or area.<br>",
            "<b>Index of dispersion:</b> quadrat variance over mean, 1 for a random pattern.<br>",
            "<b>Calibrated eps:</b> the largest eps at which a randomised pattern in the cell's footprint clusters no more than the target fraction.")),
          selectInput("hot_unit", "Unit",
            choices = c("Trace ID centroids (one point per blink)" = "tid",
                        "De-blinked molecule centroids"           = "cluster",
                        "Auto (molecules if de-blinking is on)"   = "auto"),
            selected = "tid"),
          helpText("One point per trace, or per de-blinked molecule. Raw localizations are not offered."),
          .detail(HTML(
            "Localizations within one trace sit nanometres apart, so their density reflects how long a fluorophore stayed on, not how molecules are arranged.<br><br>",
            "<b>Trace IDs</b> count blinking events; <b>de-blinked molecules</b> merge repeat blinks and give a slightly more conservative result.")),
          checkboxInput("hot_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          helpText("Spheres and cubes instead of circles and squares. A thin membrane is 2-D."),
          numericInput("hot_k", "k for local density", value = 5, min = 1, max = 50, step = 1),
          helpText("Small k is noisy; large k blurs a focus. A focus needs about 2k units to be resolved."),
          numericInput("hot_enh_thresh", "Enhancement threshold (x cell mean)",
                       value = 5, min = 1.5, max = 100, step = 0.5),
          helpText("The multiple that the reported fraction counts, and the dashed line on the distribution plot."),
          numericInput("hot_quadrat", "Quadrat size (nm)", value = 250, min = 20, step = 10),
          helpText("Set it near the size of the features you are looking for."),
          tags$hr(style = "margin:10px 0;"),
          tags$div(style = "font-weight:600; color:#2c3e50; font-size:0.92em;", "Focus detection"),
          numericInput("hot_minpts", "minPts", value = 5, min = 2, max = 50, step = 1),
          helpText("Smallest number of units DBSCAN will call a focus. Held fixed across cells."),
          numericInput("hot_target_frac", "False-cluster target", value = 0.05,
                       min = 0.005, max = 0.5, step = 0.01),
          helpText("eps is the largest value at which a randomised pattern clusters no more than this."),
          numericInput("hot_nsim", "Randomisations per eps", value = 15, min = 3, max = 100, step = 1),
          helpText("Randomised patterns per candidate eps."),
          tags$hr(style = "margin:10px 0;"),
          tags$div(style = "font-weight:600; color:#2c3e50; font-size:0.92em;", "Footprint"),
          tags$div(style = "font-size:0.85em; color:#445; line-height:1.5; margin-bottom:8px;",
            "Taken from ", tags$b("Density & area normalization"), " at the top of the tab, so this panel, ",
            "Ripley's L and g(r) measure the same cell area. The gap-fill radius is shared across files."),
          .detail(HTML("This panel always uses the footprint, whatever the footprint switch says: a cell imaged half in frame would otherwise read at half its density, doubling every enhancement. The bbox and hull modes have the same effect, since neither gives a cell mask.")),
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
            "Groups apart vertically at matched density differ in organisation; groups along one line differ in expression level, as median NN distance and Clark-Evans R show."),
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
            HTML("Cluster-occupancy counting, in the SMLM molecular-counting lineage (<b>Annibale et&nbsp;al. 2011</b>; <b>Sengupta et&nbsp;al. 2011</b>)."),
            HTML("Counts raw <b>trace IDs</b> or <b>de-blinked molecules</b> inside each median cluster; comparing the two shows how much of the count is blinking."),
            HTML("An occupancy count, not stoichiometry: blinking inflates it by an unknown, dye- and buffer-dependent factor, and no dark-time calibration is applied. Compare only conditions acquired identically.")),
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
            HTML("Tool-original. Detection rate over acquisition time for clustered and background traces, compared with the two-sample <b>Kolmogorov\u2013Smirnov</b> distance."),
            HTML("Each trace counts once, at its median time. Clustered and background traces come from the same file, so acquisition length, laser power and sample cancel; normalized modes compare timing rather than abundance."),
            HTML("Descriptive: bleaching, blinking, labelling density and re-targeting all shape the curve. The KS gap is a distance, not a test, because the two samples share one detection budget.")),
          helpText("How detection of clustered vs background traces evolves over the acquisition."),
          .detail("Traces inside Median Cluster structures against unclustered background. Requires Median Clustering."),
          numericInput("det_time_bins", "Time bins", value = 40, min = 5, max = 300),
          radioButtons("det_time_yaxis", "Y axis",
                       choices = c("Count per bin"         = "count",
                                   "Cumulative"            = "cumulative",
                                   "Normalized rate"       = "norm",
                                   "Normalized cumulative" = "norm_cumulative"),
                       selected = "count"),
          helpText("Normalized modes compare timing; the KS gap is the largest distance between the two cumulative curves."),
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
            HTML("Tool-original. Per cluster, crowding (median k-NN distance among its members) against the shift of its median detection time from the file's background, with an optional loess trend (<b>Cleveland 1979</b>)."),
            HTML("Both axes are referenced within the same file, so differences in acquisition length, bleaching and labelling cancel. A permutation null corrects for larger clusters landing later by sampling alone."),
            HTML("A correlation: crowding co-varies with labelling efficiency, structure size and photophysics, so a trend is a hypothesis for a control experiment such as a labelling titration.")),
          helpText("One point per Median Cluster: does a denser structure get detected later?"),
          .detail("Crowding: median k-NN distance among the cluster's members (smaller is more crowded). ",
                  "Shift: the cluster's median time minus the background median in that file; positive is later."),
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
            .detail("Shuffles detection times within each file while keeping every cluster's membership, size ",
                    "and geometry, so the band holds the part of the shift due to cluster size alone. The ",
                    "Spearman p tests the trend and ignores a constant offset.")
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
            HTML("Tool-original. Nearest-neighbour spacing inside each cluster (<b>Clark &amp; Evans 1954</b>) against its size, as member count or convex-hull area/volume (Barber, Dobkin &amp; Huhdanpaa 1996)."),
            HTML("Spacing is measured between de-blinked centroids, and each cluster contributes one point regardless of how many localizations it has."),
            HTML("Under DBSCAN every within-cluster distance is capped by eps2, so a flat trend at large size may be that ceiling; re-run under HDBSCAN or OPTICS first. Hull size is inflated by single outlying members.")),
          helpText("Does a bigger structure pack its members tighter, looser, or the same?"),
          .detail("Y: median over members of the mean distance to their k nearest members. Flat: constant ",
                  "spacing; rising: the structure balloons; falling: it densifies. Requires Median Clustering."),
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
            HTML("Orthogonal-distance maximum-likelihood fits ranked by <b>AICc</b> (Hurvich &amp; Tsai 1989), <b>AIC</b> (Akaike 1974) or <b>BIC</b> (Schwarz 1978); the 2D view is a PCA projection (Pearson 1901). The shape dictionary is tool-original."),
            HTML("Fits run on de-blinked centroids, so a long trace cannot pull the fit, and optionally weight each trace by its photon count. AICc blanks out when a cluster has too few points for a shape."),
            HTML("The dictionary (sphere, ellipsoid, helix, ring, arc) is imposed, and the ranking always names a winner, even when the true shape is absent. With few points several shapes are indistinguishable. Each parameter is a point estimate without an error bar.")),
          helpText("Step through individual Median Cluster structures from every loaded file."),
          .detail("Loads clusters from every file with the current de-blinking and clustering. Points follow ",
                  "de-blinking: one per molecule when it is on, one per trace ID when off, where a molecule ",
                  "that blinked six times pulls the fit six times as hard. Requires Median Clustering."),
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
            helpText("RMSD favours flexible shapes; AICc suits small clusters."),
            checkboxInput("viewer_weight_eco",
                          HTML("Weight fits by photon count (<code>eco</code>)"), value = TRUE),
            helpText(style = "font-size:0.82em; color:#666; margin-top:-4px;",
              "A trace's total photon count is its inverse-variance weight."),
            .detail("The header reports n_eff (Kish effective sample size); far below the trace count means a few ",
              "bright traces carry the fit. Criterion values are not comparable between weighted and ",
              "unweighted mode. Needs an eco column."),
            tags$div(style = "margin-top:2px; padding:5px 8px; border-left:3px solid #b26a00; background:rgba(178,106,0,0.06); font-size:0.82em;",
              HTML("<b>De-blink and filter double events first:</b> the fit assumes distinct molecules, and a bright double event gets the most weight.")),
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
            HTML("Rg\u00b2 = \u27e8\u2016p \u2212 p\u0304\u2016\u00b2\u27e9, an RMS size, with gyration-tensor shape descriptors (asphericity, <b>Rudnick &amp; Gaspari 1986</b>; linear/planar/spherical, <b>Westin et&nbsp;al. 1997/2002</b>)."),
            HTML("Points are median positions, never raw localizations. <b>Structure size</b> uses one point per molecule; <b>label footprint</b> uses one per trace of a molecule, so it measures the blink and label scatter."),
            HTML("Small clusters read low; the n/(n\u22121) correction removes most of that, but 2-member clusters stay indicative. Structure size carries the label footprint as a floor, so compare the two measures.")),
          helpText("Per-cluster Rg across all loaded files: structure size, or the labelling footprint."),
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
            "These tools measure how two channels sit relative to each other. With a multiplex build, pick the pair below.")
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
                HTML("<b>Random labelling</b> keeps the points and shuffles channel membership, so the cell outline cancels; two species in different parts of the cell separate on that alone. <b>Independent random patterns</b> inherit any error in the window outline and over-call co-localization when the cell does not fill the frame.")),
              helpText(style = "color:#7a5c00;",
                HTML("<b>Registration.</b> Channels from separate files or time windows can be offset by drift, which reads as reduced co-localization. Check a shared structure or fiducials first; the DCR split shares one frame.")),
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
            HTML("For each channel-A point, the distance to its k-th nearest channel-B point (<b>Clark &amp; Evans 1954</b>; <b>Dixon 2002</b>), and its cumulative curve: the fraction of A with a B partner within d. Object-based co-localization after <b>Lachmanovich et&nbsp;al. 2003</b> and <b>Dunn et&nbsp;al. 2011</b>."),
            HTML("Suits MINFLUX, which has coordinates rather than an image. Distances are always between channels, on de-blinked centroids; the Median Cluster modes pair only molecules inside structures."),
            HTML("No null model: every A point has a nearest B, so judge against L\u2081\u2082 or g\u2081\u2082 with their envelopes. <b>DCR crosstalk</b> puts mis-assigned molecules at near-zero distance; read the co-localized fraction against the per-group contamination from the DCR split.")),
          helpText("Every distance is measured strictly between the two channels."),
          .detail("Each channel is de-blinked and clustered separately. The Median Cluster modes drop ",
                  "unclustered background from both sides, so the result is how far channel-1 structures ",
                  "sit from channel-2 structures."),
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
          h4("Bivariate Ripley's L\u2081\u2082(r)"),
          .method_box("green",
            HTML("Cross-type K-function (<b>Lotwick &amp; Silverman 1982</b>; <b>Dixon 2002</b>) with the L-transform: L\u2081\u2082(r)&nbsp;\u2212&nbsp;r &gt; 0 is co-localization at scale r, &lt; 0 exclusion."),
            HTML("Runs on de-blinked centroids of the two channels, with the edge correction and window of L(r). The default null is <b>random labelling</b> of the observed points, so the cell outline cancels. Median Cluster modes average curves computed inside each structure. The DCR split shares one coordinate frame, so there is no channel registration to calibrate."),
            HTML("<b>DCR crosstalk</b> biases L\u2081\u2082 upward at short r, toward false co-localization: compare a signal with the contamination the DCR split reports, and pool DCR at upload where possible. Asymmetric labelling skews it too, and proximity is not interaction.")),
          helpText("L\u2081\u2082(r) \u2212 r > 0 is co-localization; < 0 is spatial exclusion."),
          .detail("Median Cluster modes average the curves over structures. Up to 2,000 points per ",
                  "channel per file or group."),
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
          checkboxInput("ripley12_envelope", "Show envelope (Monte-Carlo)", value = FALSE),
          helpText("Read significance from the dashed global band and its p-value. First few files only; not in aggregate mode."),
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
          h4("Pair Cross-Correlation g\u2081\u2082(r)"),
          .method_box("green",
            HTML("Cross pair-correlation: A\u2013B pairs in a shell around r over the shell's area or volume; g\u2081\u2082 = 1 at random, &gt;1 attraction, &lt;1 repulsion. The PC-PALM statistic (<b>Sengupta et&nbsp;al. 2011</b>; <b>Veatch et&nbsp;al. 2012</b>)."),
            HTML("Local in r, so it puts the interaction distance at a peak and keeps a short-range labelling artefact in the first bins. Same units, window and random-labelling null as L\u2081\u2082."),
            HTML("Noisier than L\u2081\u2082, since each shell holds few pairs; a first-bin peak is often binning. DCR crosstalk piles into the first shells, where a real interaction would also appear.")),
          helpText("g\u2081\u2082(r) = 1 at random; > 1 is attraction at scale r, < 1 is repulsion."),
          selectInput("pcf12_unit", "Unit",
            choices = c("TID centroids" = "tid", "De-blinked molecule centroids" = "cluster"),
            selected = "tid"),
          checkboxInput("pcf12_use_3d", "Use 3-D coordinates (z)", value = FALSE),
          numericInput("pcf12_r_max", "r max (nm, 0 = auto)", value = 0, min = 0, step = 10),
          numericInput("pcf12_r_steps", "r steps", value = 60, min = 10, max = 300),
          checkboxInput("pcf12_aggregate", "Aggregate all datasets into one", value = FALSE),
          checkboxInput("pcf12_envelope", "Show envelope (Monte-Carlo)", value = FALSE),
          helpText("Read significance from the dashed global band and its p-value. First few files only; not in aggregate mode."),
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
            HTML("Kernel density of each channel on a shared grid (<b>Silverman 1986</b>; <code>MASS::kde2d</code>), and the Pearson correlation of the two surfaces as a co-localization index (cf. <b>Dunn et&nbsp;al. 2011</b>), with an optional randomization test."),
            HTML("MINFLUX has no pixels, so the density surface is made by the bandwidth. The r-vs-bandwidth sweep shows the whole family of answers, and the default null is a toroidal random shift of channel 2 (<b>Lotwick &amp; Silverman 1982</b>), which keeps its own clustering."),
            HTML("<b>r is relative:</b> a small bandwidth drives it toward 0 and a large one toward 1, so compare only at one bandwidth. With the footprint window on the test is conservative; with it off it rejects too readily for a structure filling part of the view. For a positive claim use g\u2081\u2082(r).")),
          helpText("Density map per channel for the active file: Ch1 as red contours, Ch2 as blue."),
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
          helpText("Pearson r across bandwidths, showing how much it depends on the smoothing."),
          numericInput("kde_costes_nsim", "Randomizations (0 = off)", value = 100, min = 0, max = 500, step = 10),
          helpText("p = fraction of randomized r at or above the observed."),
          radioButtons("kde_null_mode", "Randomization null",
            choices = c("Random shift of Ch2 (recommended)" = "shift",
                        "Uniform CSR resampling of Ch2"     = "csr"),
            selected = "shift"),
          helpText(HTML("<b>Random shift</b> breaks only the cross-channel registration; <b>Uniform CSR</b> tests against complete randomness.")),
          .detail(HTML("Random shift translates Ch2 with wrap-around; uniform CSR replaces it with uniform points. Both randomize over the analysis window, so switch the footprint window on for a cell that fills only part of the view.")),
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
            "Primary reference for each tool. Analysis methods specific to MINFLUX are thinly covered ",
            "in the literature, so most tools are established SMLM and point-pattern methods adapted ",
            "to MINFLUX data by this app, one considered reading rather than a community standard. Panels marked ",
            tags$span("EXPLORATORY", style = "color:#e67e22; font-weight:700;"),
            " were built for this app and have not been validated by other groups."),
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
                            ". Labelling-model inverse (this app), checked against the simulator; contour via MST (Prim 1957), anisotropy coefficients (Westin et al. 1997/2002)")),
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
