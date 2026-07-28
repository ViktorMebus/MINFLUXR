## ─────────────────────────────────────────────────────────────────────────────
## sim_ui.R  —  Simulation tab UI  (MINFLUX-only; style mirrors SMLM tab)
## Phase 1: Synthetic localizations  |  Phase 2: Diffusion tracks
## ─────────────────────────────────────────────────────────────────────────────

sim_ui <- fluidPage(
  titlePanel("MINFLUXR - MFX Simulator"),

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

  tabsetPanel(
    id = "sim_phase_tabs",

    # ══════════════════════════════════════════════════════════════════════════
    #  Phase 1 — Localizations
    # ══════════════════════════════════════════════════════════════════════════
    tabPanel("SMLM",
      # ── Sticky action bar ─────────────────────────────────────────────────
      tags$div(
        style = paste(
          "position: sticky; top: 0; z-index: 200;",
          "background: #f8f9fa; border-bottom: 1px solid #dee2e6;",
          "padding: 6px 15px; display: flex; align-items: center; gap: 12px;"
        ),
        actionButton("sim_run", "Generate",
                     icon  = icon("play"),
                     class = "btn-primary btn-sm"),
        radioButtons("sim_view_mode", NULL,
                     choices  = c("2D" = "2d", "3D" = "3d"),
                     selected = "2d", inline = TRUE,
                     width    = "auto"),
        tags$span(style = "color:#888; font-size:0.82em;",
                  "Fix ground truth seed to compare labeling strategies on the same structure.")
      ),
      sidebarLayout(
        sidebarPanel(
          style = "height:90vh; overflow-y:auto;",

          # ── Simulation mode ─────────────────────────────────────────────────
          radioButtons("sim_mode", "Simulation mode",
                       choices = c("Normal"      = "normal",
                                   "Acquisition" = "acquisition"),
                       selected = "normal", inline = TRUE),
          conditionalPanel("input.sim_mode == 'normal'",
            helpText("Visualize a structure with a labeling strategy.")),
          conditionalPanel("input.sim_mode == 'acquisition'",
            helpText("Lay out a field of view, place structures, and model the",
                     "acquisition (detection rate, duration, background). Adds",
                     "the Acquisition timeline tab.")),
          tags$hr(),

          # ── Ground-truth structure ──────────────────────────────────────────
          tags$details(
            open = NA, class = "sim-group",
            tags$summary(tags$b("Ground-truth Structure")),

          conditionalPanel("input.sim_shape != 'oligoring'",
            numericInput("sim_n_proteins", "Number of proteins / molecules",
                         value = 10, min = 1, max = 50000, step = 10)
          ),

          selectInput("sim_shape", "Shape",
                      choices = c(
                        "Cluster"               = "cluster",
                        "Ring"                  = "ring",
                        "Arc"                   = "arc",
                        "Double ring"           = "double_ring",
                        "Ellipse"               = "ellipse",
                        "Hollow sphere"         = "hollow_sphere",
                        "Filament"              = "filament",
                        "Oligomeric ring stack" = "oligoring",
                        "Image (confocal)"      = "image"
                      ),
                      selected = "filament"),

          # Cluster
          conditionalPanel("input.sim_shape == 'cluster'",
            numericInput("sim_cluster_radius", "Cluster radius (nm)",
                         value = 100, min = 1, max = 5000, step = 5),
            selectInput("sim_cluster_dist", "Point distribution",
                        choices = c("Uniform (disk / sphere)" = "uniform",
                                    "Gaussian"                = "gaussian"),
                        selected = "uniform")
          ),

          # Ring
          conditionalPanel("input.sim_shape == 'ring'",
            numericInput("sim_ring_diam",  "Ring diameter (nm)",
                         value = 200, min = 1, max = 5000, step = 10),
            numericInput("sim_ring_width", "Ring band width (nm)",
                         value = 20, min = 0, max = 500, step = 2),
            conditionalPanel("input.sim_view_mode == '3d'",
              numericInput("sim_ring_z_thick", "Z thickness (nm)",
                           value = 20, min = 0, max = 500, step = 2)
            )
          ),

          # Ellipse
          conditionalPanel("input.sim_shape == 'ellipse'",
            fluidRow(
              column(6, numericInput("sim_ellipse_a", "Semi-axis a (nm)",
                                     value = 150, min = 1, max = 5000, step = 5)),
              column(6, numericInput("sim_ellipse_b", "Semi-axis b (nm)",
                                     value = 80,  min = 1, max = 5000, step = 5))
            ),
            conditionalPanel("input.sim_view_mode == '3d'",
              numericInput("sim_ellipse_c", "Semi-axis c / Z (nm)",
                           value = 50, min = 1, max = 5000, step = 5)
            ),
            radioButtons("sim_ellipse_fill", NULL,
                         choices  = c("Filled" = "filled",
                                      "Surface only" = "perimeter"),
                         selected = "filled", inline = TRUE)
          ),

          # Filament
          conditionalPanel("input.sim_shape == 'filament'",
            helpText("N proteins are spaced at regular intervals along",
                     "a straight filament (actin ~7 nm, microtubule ~8 nm)."),
            numericInput("sim_fil_spacing",
                         "Molecule spacing (nm)",
                         value = 15, min = 0.1, max = 1000, step = 0.5),
            numericInput("sim_fil_angle",
                         "Orientation angle in XY (\u00B0)",
                         value = 0, min = 0, max = 360, step = 5),
            fluidRow(
              column(6, numericInput("sim_fil_jitter_lat",
                                     "Lateral curl (nm / step SD)",
                                     value = 4, min = 0, max = 200, step = 0.5)),
              column(6, numericInput("sim_fil_jitter_ax",
                                     "Axial jitter (nm SD)",
                                     value = 0, min = 0, max = 200, step = 0.5))
            ),
            helpText("Lateral curl is a random walk: each molecule inherits",
                     "the previous one\u2019s offset plus a small step, so the",
                     "filament curves and meanders rather than scatters."),
            checkboxInput("sim_fil_helix",
                          "Helical twist (e.g. DNA-bound protein)", FALSE),
            conditionalPanel("input.sim_fil_helix == true",
              fluidRow(
                column(6, numericInput("sim_fil_helix_r",
                                       "Helix radius (nm)",
                                       value = 1, min = 0.1, max = 50, step = 0.1)),
                column(6, numericInput("sim_fil_helix_pitch",
                                       "Pitch (nm / turn)",
                                       value = 36, min = 1, max = 500, step = 1))
              ),
              helpText("Default 36 nm pitch = B-DNA repeat.")
            ),
            conditionalPanel("input.sim_view_mode == '3d'",
              numericInput("sim_fil_z_angle",
                           "Z tilt (\u00B0, \u22120\u200990 to 90)",
                           value = 0, min = -90, max = 90, step = 5)
            )
          ),

          # Hollow sphere
          conditionalPanel("input.sim_shape == 'hollow_sphere'",
            helpText("Uniform surface sampling of a sphere shell.",
                     "Models vesicles, nuclei, plasma membrane proteins."),
            numericInput("sim_hsph_radius", "Sphere radius (nm)",
                         value = 150, min = 1, max = 5000, step = 5),
            numericInput("sim_hsph_thick", "Shell thickness (nm)",
                         value = 10, min = 0, max = 500, step = 1),
            helpText("In 2D mode rendered as a ring (XY cross-section).")
          ),

          # Double ring
          conditionalPanel("input.sim_shape == 'double_ring'",
            helpText("Two parallel rings separated along Z.",
                     "Models nuclear pore complex (~35 nm sep.),",
                     "centrioles, etc."),
            numericInput("sim_dring_diam", "Ring diameter (nm)",
                         value = 200, min = 1, max = 5000, step = 10),
            numericInput("sim_dring_width", "Ring band width (nm)",
                         value = 20, min = 0, max = 500, step = 2),
            numericInput("sim_dring_sep", "Ring separation (nm)",
                         value = 35, min = 1, max = 2000, step = 5),
            conditionalPanel("input.sim_view_mode == '3d'",
              numericInput("sim_dring_z_thick", "Z band thickness (nm)",
                           value = 10, min = 0, max = 500, step = 2)
            ),
            helpText("In 2D mode both rings are projected onto the XY plane.")
          ),

          # Arc
          conditionalPanel("input.sim_shape == 'arc'",
            helpText("A partial ring spanning a defined angular extent.",
                     "Models crescent domains, incomplete assemblies."),
            numericInput("sim_arc_diam", "Arc diameter (nm)",
                         value = 200, min = 1, max = 5000, step = 10),
            numericInput("sim_arc_width", "Band width (nm)",
                         value = 20, min = 0, max = 500, step = 2),
            fluidRow(
              column(6, numericInput("sim_arc_start",
                                     "Start angle (\u00B0)",
                                     value = 0, min = 0,
                                     max = 360, step = 5)),
              column(6, numericInput("sim_arc_extent",
                                     "Arc extent (\u00B0)",
                                     value = 270, min = 1,
                                     max = 360, step = 5))
            ),
            conditionalPanel("input.sim_view_mode == '3d'",
              numericInput("sim_arc_z_thick", "Z thickness (nm)",
                           value = 20, min = 0, max = 500, step = 2)
            )
          ),

          # Oligomeric ring stack
          conditionalPanel("input.sim_shape == 'oligoring'",
            helpText("Fixed subunit positions for X-mer rings stacked along Z.",
                     "subcomplexes, etc. GT count = subunits × stacks."),
            fluidRow(
              column(6, numericInput("sim_oligoring_nsub",
                                    "Subunits per ring",
                                    value = 6, min = 2, max = 48, step = 1)),
              column(6, numericInput("sim_oligoring_nstack",
                                    "Rings (stacks)",
                                    value = 1, min = 1, max = 20, step = 1))
            ),
            numericInput("sim_oligoring_diam", "Ring diameter (nm)",
                         value = 12, min = 1, max = 5000, step = 0.5),
            numericInput("sim_oligoring_base_angle",
                         "Base rotation (°)",
                         value = 0, min = -180, max = 180, step = 1),
            helpText("Fixed angular offset applied to the whole assembly.",
                     "Stacks on top of ‘Random XY orientation’ if that is also enabled."),
            conditionalPanel("input.sim_oligoring_nstack > 1",
              numericInput("sim_oligoring_sep",
                           "Axial ring separation (nm)",
                           value = 7, min = 0.5, max = 2000, step = 0.5),
              numericInput("sim_oligoring_angle_off",
                           "Angular offset per ring (°)",
                           value = 0, min = -180, max = 180, step = 1),
              helpText("Rotate each successive ring by this angle.",
                       "E.g., 30° for a hexamer/heptamer symmetry mismatch.")
            )
          ),

          # Image-guided
          conditionalPanel("input.sim_shape == 'image'",
            fileInput("sim_img_file",
                      "Upload confocal image (PNG / JPEG / TIFF)",
                      accept = c("image/png","image/jpeg","image/tiff",
                                 ".tif",".tiff")),
            numericInput("sim_img_px_size", "Pixel size (nm / pixel)",
                         value = 10, min = 0.1, max = 10000, step = 1),
            sliderInput("sim_img_threshold",
                        "Intensity threshold (ignore below)",
                        min = 0, max = 1, value = 0.05, step = 0.01),
            fluidRow(
              column(6, numericInput("sim_img_min_trace",
                                     "Min trace length (dark regions)",
                                     value = 2, min = 1, max = 100, step = 1)),
              column(6, numericInput("sim_img_max_trace",
                                     "Max trace length (bright regions)",
                                     value = 15, min = 1, max = 500, step = 1))
            ),
            helpText("Bright pixels attract more molecules and generate",
                     "longer traces. Source image shown as faint overlay."),
            conditionalPanel("input.sim_view_mode == '3d'",
              tags$hr(),
              h5("3D Depth Model"),
              helpText("Bright regions get more axial depth;",
                       "dim regions stay near the focal plane (Z\u22480)."),
              sliderInput("sim_img_z_spread",
                          "Z depth for bright spots (nm)",
                          min = 10, max = 1000, value = 200, step = 10),
              numericInput("sim_img_z_bg",
                           "Background Z scatter \u03C3 (nm, all regions)",
                           value = 20, min = 0, max = 300, step = 5),
              helpText("Z\u2009SD = max(bg\u2009\u03C3,\u2009depth\u2009\u00D7\u2009local\u2009intensity).")
            )
          ),

          # ── Random orientation & multi-copy ─────────────────────────────────
          conditionalPanel("input.sim_shape != 'image'",
            tags$hr(style = "margin-top:6px; margin-bottom:6px;"),
            checkboxInput("sim_rotate",
                          "Random XY orientation each Generate", FALSE),
            conditionalPanel("input.sim_rotate == true && !input.sim_multi_copy",
              sliderInput("sim_rotate_range",
                          "Max rotation (\u00B0, uniform \u00B1)",
                          min = 5, max = 180, value = 180, step = 5),
              helpText("Rotates the whole structure by a random angle on every",
                       "Generate. Disabled when Multi-copy is on (each copy",
                       "is then rotated independently).")
            ),
            # Multi-copy tiling \u2014 Normal mode only (Acquisition mode places
            # structures in the FOV via its own "Number of structures" control).
            conditionalPanel("input.sim_mode == 'normal'",
              tags$hr(style = "margin-top:6px; margin-bottom:6px;"),
              checkboxInput("sim_multi_copy",
                            "Multi-copy tiling", FALSE),
              helpText("Place N independent copies of the structure at random",
                       "positions within a square FOV."),
              conditionalPanel("input.sim_multi_copy == true",
                fluidRow(
                  column(6, numericInput("sim_copies_n",
                                         "Number of copies",
                                         value = 5, min = 2, max = 200, step = 1)),
                  column(6, numericInput("sim_copies_fov",
                                         "FOV side (nm)",
                                         value = 2000, min = 100,
                                         max = 100000, step = 100))
                ),
                checkboxInput("sim_copies_rotate",
                              "Random orientation per copy", TRUE),
                helpText("Each copy is independently rotated (full 360\u00B0).")
              )
            )
          ),

          # Acquisition setup (Acquisition mode). The FOV size itself lives in the
          # Display section (shared with Normal-mode ROI preview).
          conditionalPanel("input.sim_mode == 'acquisition'",
            conditionalPanel("input.sim_shape != 'image'",
              tags$hr(),
              numericInput("sim_acq_n_structures",
                           "Number of structures in FOV",
                           value = 4, min = 1, max = 500, step = 1),
              helpText("Independently rotated copies placed at random positions",
                       "across the X\u00D7Y field (e.g. 4 filaments).",
                       "Set the FOV size under Display.")
            )
          )
          ),  # /Ground-truth Structure (collapsible)

          tags$hr(),

          # ── Labeling ────────────────────────────────────────────────────────
          tags$details(
            open = NA, class = "sim-group",
            tags$summary(tags$b("Labeling Strategy")),
          helpText("Each labeled molecule generates N fluorophore positions.",
                   "Each fluorophore = one trace ID in the output."),

          selectInput("sim_label_type", "Label type",
                      choices = c(
                        "Direct tag (SNAP / Halo / Click, DOL=1)",
                        "Nanobody (DOL \u22481)",
                        "Primary monoclonal IgG (DOL \u22481)",
                        "Primary monoclonal IgG (NHS, DOL \u22482\u20134)",
                        "Primary polyclonal IgG (NHS, DOL \u22482\u20134)",
                        "Primary monoclonal IgG + secondary nanobody",
                        "Primary monoclonal IgG + secondary IgG",
                        "Primary polyclonal IgG + secondary IgG",
                        "Streptavidin (biotinylated Ab + SA-dye)",
                        "Custom"
                      ),
                      selected = "Direct tag (SNAP / Halo / Click, DOL=1)"),

          uiOutput("sim_label_info"),

          conditionalPanel("input.sim_label_type == 'Custom'",
            fluidRow(
              column(6, numericInput("sim_fluor_custom",
                                     "Fluorophores / molecule",
                                     value = 1, min = 1, max = 50, step = 1)),
              column(6, numericInput("sim_link_custom_nm",
                                     "Linker distance (nm)",
                                     value = 5, min = 0, max = 200, step = 0.5))
            ),
            numericInput("sim_dye_spread_custom",
                         "Dye spread within complex \u03C3 (nm)",
                         value = 0, min = 0, max = 100, step = 0.5),
            helpText("For multi-dye labels: Gaussian SD of individual dye",
                     "positions around the shared complex centre.",
                     "Set to 0 for a single dye."),
            numericInput("sim_footprint_custom",
                         "Label complex footprint (nm)",
                         value = 8, min = 0, max = 100, step = 0.5),
            helpText("Approximate diameter of excluded zone around each",
                     "labeled protein \u2014 used for steric exclusion check.",
                     "~2 nm (direct tag) to ~20 nm (2\u00B0 IgG).")
          ),

          sliderInput("sim_label_eff", "Labeling efficiency (%)",
                      min = 1, max = 100, value = 80, step = 1, post = "%"),

          checkboxInput("sim_ldir_aware",
                        "Structure-aware linker displacement",
                        value = TRUE),
          helpText("When enabled, linker error is constrained by structure geometry:",
                   "filaments show no axial displacement (antibody can only attach",
                   "from the side); rings/spheres displace radially outward."),

          numericInput("sim_protein_mw",
                       "Target protein MW (kDa)  — 0 = not modelled",
                       value = 0, min = 0, max = 5000, step = 1),
          helpText("Displaces the epitope from the protein centroid by the",
                   "protein surface radius (r ≈ 0.66 × MW¹³ nm).",
                   "Adds to the linker error budget. Larger proteins push",
                   "the antibody further from the true molecule centre.",
                   "Also tightens the steric exclusion check."),

          # TID over-counting belongs to NORMAL mode only. In Acquisition mode
          # over-counting instead emerges as detection rate × duration, capped by
          # the "Detections per fluorophore before bleaching" photon budget.
          conditionalPanel("input.sim_mode == 'normal'",
          tags$hr(style = "margin:6px 0;"),
          tags$b("Blinking / TID over-counting"),
          checkboxInput("sim_blink_model",
                        "Model dye re-activation (TID splitting)",
                        value = FALSE),
          helpText("Each fluorophore can re-activate after bleaching, generating",
                   "additional TIDs from the same physical dye — apparent over-counting."),
          conditionalPanel("input.sim_blink_model == true",
            sliderInput("sim_blink_mean",
                        "Mean extra TIDs per fluorophore (λ)",
                        min = 0, max = 3, value = 0.3, step = 0.05),
            helpText("Total apparent TIDs ≈ (1 + λ) × labeled fluorophores.",
                     "Presets use typical values for organic dyes in aqueous buffer.")
          )
          )  # /Blinking over-counting (Normal mode only)
          ),  # /Labeling (collapsible)

          tags$hr(),

          # ── Display ─────────────────────────────────────────────────────────
          tags$details(
            class = "sim-group",
            tags$summary(tags$b("Display")),

          numericInput("sim_point_size", "Point size",
                       value = 4, min = 1, max = 20),

          selectInput("sim_color_mode", "Color localizations by",
                      choices = c("Trace ID (default)" = "tid",
                                  "Depth (z)"          = "z",
                                  "Time"               = "tim",
                                  "Local density"      = "locdens"),
                      selected = "tid"),
          conditionalPanel(
            "['z','tim','locdens'].includes(input.sim_color_mode)",
            selectInput("sim_numeric_palette", "Numeric colorscale",
                        choices = NUMERIC_PALETTES, selected = "Spectral"),
            checkboxInput("sim_palette_reverse", "Reverse colorscale", FALSE),
            conditionalPanel("input.sim_color_mode == 'locdens'",
              sliderInput("sim_dens_radius_nm", "Density radius (nm)",
                          min = 10, max = 5000, value = 500, step = 10),
              helpText("Counts neighbouring localizations within this radius."))
          ),

          checkboxInput("sim_show_gt",
                        "Show molecule positions (ground truth)", TRUE),
          checkboxInput("sim_show_locs",
                        "Show localizations (coloured by trace ID)", TRUE),
          checkboxInput("sim_show_gt_lines",
                        "Show lines: GT molecule \u2192 each TID median",
                        FALSE),
          checkboxInput("sim_dark_mode", "Dark mode", FALSE),
          checkboxInput("sim_hide_grid", "Hide gridlines", FALSE),

          tags$hr(style = "margin:6px 0;"),
          # Field of view / ROI window — pins the plot axes to a real-camera-style
          # FOV. Always on in Acquisition mode (it also lays out the structures);
          # optional in Normal mode to preview the structure inside an ROI.
          conditionalPanel("input.sim_mode == 'normal'",
            checkboxInput("sim_view_fov",
                          "Pin axes to a field of view (ROI)", FALSE),
            helpText("Preview the structure inside a real-camera-style ROI window",
                     "instead of auto-zooming to the localizations.")
          ),
          conditionalPanel(
            "input.sim_view_fov == true || input.sim_mode == 'acquisition'",
            fluidRow(
              column(4, numericInput("sim_fov_x", "FOV X (nm)",
                                     value = 4000, min = 10, step = 100)),
              column(4, numericInput("sim_fov_y", "FOV Y (nm)",
                                     value = 4000, min = 10, step = 100)),
              column(4, numericInput("sim_fov_z", "FOV Z (nm)",
                                     value = 1000, min = 10, step = 100))
            ),
            helpText("Pins the plot to a real-camera-style FOV with true aspect",
                     "ratios (Z compressed in 3D). Localizations are not cropped.")
          ),

          tags$hr(style = "margin:6px 0;"),
          conditionalPanel("input.sim_view_mode == '2d'",
            checkboxInput("sim_show_scale_bar", "Add scale bar", FALSE),
            conditionalPanel("input.sim_show_scale_bar == true",
              numericInput("sim_scale_bar_length", "Scale bar length (nm)",
                           value = 500, min = 1, step = 10),
              numericInput("sim_scale_bar_thickness", "Scale bar thickness (nm)",
                           value = 10, min = 1, step = 1),
              numericInput("sim_scale_font_size", "Font size",
                           value = 18, min = 1, step = 1)
            )
          )
          ),  # /Display (collapsible)

          tags$hr(),
          # ── Localization precision ──────────────────────────────────────────
          tags$details(
            class = "sim-group",
            tags$summary(tags$b("Localization precision")),

          radioButtons("sim_prec_mode", "Precision input",
                       choices  = c("Manual \u03C3"             = "manual",
                                    "CRLB from photon budget" = "crlb"),
                       selected = "manual", inline = TRUE),

          # Manual mode
          conditionalPanel("input.sim_prec_mode == 'manual'",
            numericInput("sim_precision_xy",
                         "MINFLUX precision XY \u03C3 (nm)",
                         value = 6, min = 0.1, max = 500, step = 0.5),
            conditionalPanel("input.sim_view_mode == '3d'",
              numericInput("sim_precision_z",
                           "MINFLUX precision Z \u03C3 (nm)",
                           value = 4, min = 0.1, max = 500, step = 0.5),
              helpText("Z precision can differ from XY;",
                       "typical MINFLUX Z \u03C3 \u2248 3\u201310 nm.")
            )
          ),

          # CRLB mode
          conditionalPanel("input.sim_prec_mode == 'crlb'",
            helpText("\u03C3 \u2248 L / (2\u221AN_eff),",
                     "N_eff = N_photons \u00D7 (1 \u2212 bg fraction).",
                     "L = effective TCP scan-pattern diameter; it isn't stored",
                     "in the data, but \u201CCalibrate from a real run\u201D derives it",
                     "as L \u2248 2\u03C3\u221AN from your measured precision + photons."),
            numericInput("sim_crlb_L",
                         "L \u2014 effective TCP diameter (nm)",
                         value = 50, min = 10, max = 2000, step = 10),
            numericInput("sim_crlb_photons",
                         "Photons per localization",
                         value = 80, min = 1, max = 5000, step = 5),
            sliderInput("sim_crlb_bg",
                        "Background fraction",
                        min = 0, max = 0.5, value = 0.05, step = 0.01),
            conditionalPanel("input.sim_view_mode == '3d'",
              numericInput("sim_crlb_z_factor",
                           "\u03C3_z / \u03C3_xy ratio",
                           value = 0.7, min = 0.1, max = 5, step = 0.05),
              helpText("\u03C3_z = \u03C3_xy \u00D7 ratio.",
                       "\u03C3_xy: ratio \u2248 0.6\u20130.8 (e.g. 4\u20135 nm Z",
                       "vs 6\u20137 nm XY).")
            ),
            uiOutput("sim_crlb_preview")
          ),

          numericInput("sim_trace_length",
                       "Trace length (mean localizations per fluorophore)",
                       value = 10, min = 1, max = 1000, step = 1),

          checkboxInput("sim_tl_geometric",
                        "Geometric trace length (constant bleach rate)",
                        value = FALSE),
          helpText("Replaces Poisson with a geometric distribution:",
                   "same mean but heavy right tail",
                   "(most traces short, a few very long)."),

          checkboxInput("sim_tl_scale_fluor",
                        "Scale trace length with fluorophore count",
                        value = FALSE),
          helpText("Divides each TID\u2019s budget by the number of TIDs on",
                   "its molecule \u2014 keeps total locs per molecule constant",
                   "regardless of labeling strategy."),
          tags$hr(style = "margin:6px 0;"),
          tags$b("Background / noise TIDs"),
          checkboxInput("sim_bg_noise",
                        "Add background / non-specific TIDs",
                        value = FALSE),
          conditionalPanel("input.sim_bg_noise == true",
            # Normal mode → fixed count; Acquisition mode → areal rate over time.
            conditionalPanel("input.sim_mode == 'normal'",
              numericInput("sim_bg_n",
                           "Number of background TIDs",
                           value = 20, min = 1, max = 50000, step = 5)
            ),
            conditionalPanel("input.sim_mode == 'acquisition'",
              numericInput("sim_bg_rate",
                           "Background detection rate (traces / hour)",
                           value = 5, min = 0, max = 1000000, step = 1),
              helpText("Total background TIDs = rate × duration; longer runs",
                       "accumulate more non-specific traces across the FOV.")
            ),
            helpText("Spread uniformly across the field of view when an ROI is",
                     "pinned (Display → Pin axes to a FOV), otherwise across the",
                     "structure's area. Models non-specific binding or",
                     "autofluorescence.")
          )
          ),  # /Localization precision (collapsible)

          tags$hr(),

          # ── Acquisition ─────────────────────────────────────────────────────
          tags$details(
            class = "sim-group",
            tags$summary(tags$b("Acquisition")),

          conditionalPanel("input.sim_mode == 'acquisition'",
          tags$b("Calibrate from a real run"),
          fileInput("sim_real_file", NULL,
                    accept = c(".npy", ".json", ".csv"),
                    buttonLabel = "Load run", placeholder = "no file"),
          uiOutput("sim_real_status"),
          helpText("Loads a MINFLUX .npy / .json and fills the FOV, duration and",
                   "detection-rate profile below from the real data, then overlays",
                   "its measured rate on the Acquisition plot."),
          sliderInput("sim_real_struct_frac",
                      "Structure fraction of detections",
                      min = 0, max = 1, value = 0.5, step = 0.05),
          helpText("Auto-estimated from spatial density (dense = structure,",
                   "sparse = background); nudge it if the split looks off. The",
                   "rest of the rate becomes background, so a low value = a",
                   "densely labeled structure the donut resolves slowly while",
                   "background piles up fast."),
          checkboxInput("sim_real_show_split",
                        "Split rate into structure + background",
                        value = FALSE),
          helpText("On the Acquisition plot, break the predicted detection rate",
                   "into its structure (decaying + bleaching) and background",
                   "(flat) parts — solid lines. If a real run is loaded, also",
                   "overlays the measured spatial split (dotted)."),

          tags$hr(style = "margin:6px 0;"),
          tags$b("Acquisition time"),
          fluidRow(
            column(6, numericInput("sim_acq_hours", "Hours",
                                   value = 4, min = 0, max = 240, step = 1)),
            column(6, numericInput("sim_acq_minutes", "Minutes",
                                   value = 0, min = 0, max = 59, step = 1))
          ),
          numericInput("sim_struct_rate",
                       "Structure detection rate (traces / hour)",
                       value = 20, min = 0, max = 1000000, step = 5),
          numericInput("sim_acq_blinks",
                       "Detections per fluorophore before bleaching",
                       value = 5, min = 1, max = 1000, step = 1),
          helpText(tags$b("Acquisition-mode TID over-counting:"),
                   "the detection rate generates repeat TIDs over time, capped",
                   "by this photon budget — the time-resolved equivalent of",
                   "the Normal-mode blink split above."),
          helpText("Total structure traces = rate × duration, capped by the",
                   "photon budget (fluorophores × detections-before-bleaching).",
                   "A short run resolves part of the structure; a long run",
                   "resolves more then over-counts; once the budget is spent the",
                   "structure plateaus while background keeps accumulating."),
          sliderInput("sim_acq_resolve_q",
                      "\u201cResolved\u201d threshold (fraction of molecules)",
                      min = 0.5, max = 1, value = 0.9, step = 0.05),
          helpText("Threshold for when the structure counts as resolved.",
                   "Traces collected after that mostly re-detect counted",
                   "molecules. See the Acquisition tab for the curves."),

          radioButtons("sim_acq_det_model", "Detection-rate model",
                       choices = c("Decaying (blinking dyes bleach)" = "decay",
                                   "Constant"                        = "constant",
                                   "Rising to plateau (activation)"  = "rise"),
                       selected = "decay"),
          helpText(tags$b("Decaying"), "= bright early then bleaches toward a floor",
                   "(front-loaded, typical of blinking fluorophores).",
                   tags$b("Constant"), "= a flat rate (only Poisson variation).",
                   tags$b("Rising"), "= almost no detection early, then climbs to a",
                   "plateau (activation-limited, e.g. slow buffer / exchange onset)."),

          conditionalPanel(
            condition = "input.sim_acq_det_model != 'constant'",
            numericInput("sim_acq_decay_min",
                         "Rate-change time constant (min)",
                         value = 15, min = 0.5, max = 6000, step = 5),
            sliderInput("sim_acq_floor_frac",
                        "Level (× characteristic rate)",
                        min = 0, max = 1, value = 0.25, step = 0.05),
            conditionalPanel(
              condition = "input.sim_acq_det_model == 'decay'",
              helpText("The characteristic rate is the", tags$b("initial peak"), ".",
                       "The rate eases from the peak toward (this × peak) over the",
                       "time constant.", tags$b("0"), "= decays to zero;",
                       tags$b("1"), "= no decay (flat).")),
            conditionalPanel(
              condition = "input.sim_acq_det_model == 'rise'",
              helpText("The characteristic rate is the", tags$b("plateau (cap)"), ".",
                       "The rate starts at (this × plateau) and climbs to the plateau",
                       "over the time constant.", tags$b("0"), "= starts from zero;",
                       tags$b("1"), "= already at the plateau (flat)."))
          ),

          checkboxInput("sim_crowding",
                        "Model emitter crowding", value = FALSE),
          conditionalPanel("input.sim_crowding == true",
            numericInput("sim_crowd_radius", "Crowding radius (nm)",
                         value = 50, min = 1, max = 1000, step = 5),
            helpText("Where labeled molecules sit closer than this, the donut",
                     "sees several emitters at once and singles them out slowly,",
                     "so each molecule's discovery rate is scaled by",
                     "1/(1 + neighbours within the radius). Densely labeled",
                     "structures then resolve slower (lots of TIDs, few unique).")
          ),
          ),

          checkboxInput("sim_fix_seed",
                        "Fix ground truth seed (locks molecule positions)",
                        value = FALSE),
          conditionalPanel("input.sim_fix_seed == true",
            numericInput("sim_seed", "Ground truth seed",
                         value = 42, min = 1,
                         max = .Machine$integer.max, step = 1)
          ),
          checkboxInput("sim_fix_label_seed",
                        "Fix labeling seed (locks fluorophore placement)",
                        value = FALSE),
          conditionalPanel("input.sim_fix_label_seed == true",
            numericInput("sim_label_seed", "Labeling seed",
                         value = 1, min = 1,
                         max = .Machine$integer.max, step = 1)
          ),
          helpText("Fix the ground truth seed to compare labeling strategies",
                   "on the same underlying molecule positions.")

          ),  # /Acquisition (collapsible)

          tags$hr(),



          tags$br(), tags$br(),
          downloadButton("sim_dl_csv", "Download localizations (.csv)",
                         style = "width:100%;"),
          tags$br(),
          downloadButton("sim_dl_gt_csv", "Download GT positions (.csv)",
                         style = "width:100%;"),
          tags$br(),
          uiOutput("sim_tid_summary"),
          uiOutput("sim_summary"),
          uiOutput("sim_precision_budget")
        ),

        mainPanel(
          tabsetPanel(
            id = "sim_main_tabs",

            # ── MINFLUX localizations ───────────────────────────────────────
            tabPanel("MINFLUX",
              plotlyOutput("sim_plot", height = "82vh")
            ),

            # ── Confocal view ────────────────────────────────────────────────
            tabPanel("Confocal view",
              fluidRow(
                column(3,
                  tags$br(),
                  helpText("640 nm excitation (Alexa647 / Atto647N),",
                           "\u03BB_em \u2248 670 nm."),
                  numericInput("sim_diff_NA", "Objective NA",
                               value = 1.4, min = 0.5,
                               max = 1.5, step = 0.05),
                  sliderInput("sim_diff_AU", "Pinhole (Airy Units)",
                              min = 0.3, max = 3.0,
                              value = 0.83, step = 0.01),
                  uiOutput("sim_diff_psf_preview"),
                  tags$hr(),
                  sliderInput("sim_diff_px",
                              "Pixel size (nm)",
                              min = 10, max = 120,
                              value = 50, step = 5),
                  helpText("Nyquist: px < FWHM / 2 for proper sampling.",
                           "Image capped at 500 \u00D7 500 px."),
                  tags$hr(),
                  numericInput("sim_diff_photons",
                               "Mean photons / fluorophore",
                               value = 45, min = 1, max = 100000, step = 5),
                  numericInput("sim_diff_bg",
                               "Background (photons / pixel)",
                               value = 35, min = 0, max = 10000, step = 1),
                  helpText("Rendered from GT molecule positions.",
                           "Poisson shot noise included.",
                           "Capped at 5 000 molecules for speed.")
                ),
                column(9,
                  plotlyOutput("sim_diff_plot", height = "82vh")
                )
              )
            ),

            # ── Acquisition timeline ─────────────────────────────────────────
            tabPanel("Acquisition",
              tags$br(),
              fluidRow(
                column(8, plotlyOutput("sim_acq_plot", height = "72vh")),
                column(4, uiOutput("sim_acq_summary"))
              )
            ),

            # ── Compare labeling strategies ──────────────────────────────────
            tabPanel("Compare strategies",
              tags$br(),
              fluidRow(
                column(4,
                  tags$b("Strategy A", style = "color:#dc2626;"),
                  selectInput("sim_cmp_label_a", NULL,
                    choices = c(
                      "Direct tag (SNAP / Halo / Click, DOL=1)",
                      "Nanobody (DOL ≈1)",
                      "Primary monoclonal IgG (DOL ≈1)",
                      "Primary monoclonal IgG (NHS, DOL ≈2–4)",
                      "Primary polyclonal IgG (NHS, DOL ≈2–4)",
                      "Primary monoclonal IgG + secondary nanobody",
                      "Primary monoclonal IgG + secondary IgG",
                      "Primary polyclonal IgG + secondary IgG",
                      "Streptavidin (biotinylated Ab + SA-dye)",
                      "Custom"
                    ),
                    selected = "Direct tag (SNAP / Halo / Click, DOL=1)"),
                  sliderInput("sim_cmp_eff_a", "Labeling efficiency A (%)",
                              min = 1, max = 100, value = 80, step = 1,
                              post = "%"),
                  conditionalPanel("input.sim_cmp_label_a == 'Custom'",
                    tags$hr(style = "margin:6px 0;"),
                    tags$small(tags$b("Custom A parameters")),
                    numericInput("sim_cmp_a_fluor", "Fluorophores / molecule",
                                 value = 2, min = 1, max = 50, step = 1),
                    numericInput("sim_cmp_a_link", "Linker distance (nm)",
                                 value = 10, min = 0, max = 200, step = 0.5),
                    numericInput("sim_cmp_a_spread", "Dye spread σ (nm)",
                                 value = 0, min = 0, max = 100, step = 0.5),
                    numericInput("sim_cmp_a_foot", "Label footprint (nm)",
                                 value = 8, min = 0, max = 100, step = 0.5)
                  )
                ),
                column(4,
                  tags$b("Strategy B", style = "color:#2563eb;"),
                  selectInput("sim_cmp_label_b", NULL,
                    choices = c(
                      "Direct tag (SNAP / Halo / Click, DOL=1)",
                      "Nanobody (DOL ≈1)",
                      "Primary monoclonal IgG (DOL ≈1)",
                      "Primary monoclonal IgG (NHS, DOL ≈2–4)",
                      "Primary polyclonal IgG (NHS, DOL ≈2–4)",
                      "Primary monoclonal IgG + secondary nanobody",
                      "Primary monoclonal IgG + secondary IgG",
                      "Primary polyclonal IgG + secondary IgG",
                      "Streptavidin (biotinylated Ab + SA-dye)",
                      "Custom"
                    ),
                    selected = "Primary monoclonal IgG + secondary IgG"),
                  sliderInput("sim_cmp_eff_b", "Labeling efficiency B (%)",
                              min = 1, max = 100, value = 80, step = 1,
                              post = "%"),
                  conditionalPanel("input.sim_cmp_label_b == 'Custom'",
                    tags$hr(style = "margin:6px 0;"),
                    tags$small(tags$b("Custom B parameters")),
                    numericInput("sim_cmp_b_fluor", "Fluorophores / molecule",
                                 value = 2, min = 1, max = 50, step = 1),
                    numericInput("sim_cmp_b_link", "Linker distance (nm)",
                                 value = 10, min = 0, max = 200, step = 0.5),
                    numericInput("sim_cmp_b_spread", "Dye spread σ (nm)",
                                 value = 0, min = 0, max = 100, step = 0.5),
                    numericInput("sim_cmp_b_foot", "Label footprint (nm)",
                                 value = 8, min = 0, max = 100, step = 0.5)
                  )
                ),
                column(4,
                  tags$br(),
                  helpText("Uses the same ground truth from the last",
                           "Generate run.  Fix the ground truth seed",
                           "to keep the structure identical across runs."),
                  actionButton("sim_compare_run", "Compare",
                               icon = icon("columns"),
                               class = "btn-primary btn-lg",
                               style = "width:100%; margin-top:4px;")
                )
              ),
              tags$hr(style = "margin:6px 0;"),
              fluidRow(
                column(6,
                  plotlyOutput("sim_cmp_plot_a", height = "55vh")
                ),
                column(6,
                  plotlyOutput("sim_cmp_plot_b", height = "55vh")
                )
              ),
              fluidRow(
                column(12,
                  uiOutput("sim_cmp_table")
                )
              )
            )
          )
        )
      )
    ), # end SMLM (Localizations)

    # ══════════════════════════════════════════════════════════════════════════
    #  Phase 2 — Diffusion / Tracking (ground-truth simulator, own sidebar)
    #  Free surface Brownian on a cell shape, sampled the MINFLUX way (cycle clock +
    #  dark cycles + heteroscedastic σ_loc + photobleaching) → MINFLUX-format tracks
    #  with a ground-truth D, for validating CVE / MSD / confinement estimators.
    #  See docs/simulation_module.md.
    # ══════════════════════════════════════════════════════════════════════════
    tabPanel("Diffusion / Tracking",
      sidebarLayout(
        sidebarPanel(
          style = "height:90vh; overflow-y:auto;",
          tags$div(class = "well well-sm", style = "font-size:0.82em; color:#555;",
            tags$b("Ground-truth diffusion simulator."), " Free Brownian motion on a cell surface, ",
            "sampled with MINFLUX's irregular Δt (dark cycles) + per-localization σ_loc + ",
            "photobleaching. Download the tracks and load them in the Tracking tab to check that ",
            "CVE / MSD recover the known D. ", tags$b("Exploratory validation tool.")),
          # ── Calibrate the sampling model from a real MINFLUX file ──────────
          tags$details(class = "sim-group", open = NA,
            tags$summary("Calibrate from real MINFLUX data"),
            helpText(style = "font-size:0.8em; color:#666;",
              "Load a real .csv / .npy / .json to auto-fill cycle time, photons, dark-cycle ",
              "probability and median track length from that dataset. ",
              tags$b("Tip:"), " to calibrate on exactly the tracks you analyse, first apply your ",
              "filters in the Tracking tab and use its ", tags$i("Export filtered data (.csv)"),
              " button, then load that file here."),
            fileInput("sim_trk_calib_file", NULL, accept = c(".csv", ".npy", ".json"),
                      buttonLabel = "Browse…", placeholder = "MINFLUX .csv / .npy / .json"),
            numericInput("sim_trk_calib_minlen", "Min trace length", value = 100, min = 5, step = 10),
            helpText(style = "font-size:0.78em; color:#666; margin-top:-8px;",
              "Calibrate only on traces ≥ this length (up to the longest), so short junk traces ",
              "don't drag the statistics down. Uses the median length of the kept traces."),
            actionButton("sim_trk_calib", "Calibrate parameters", icon = icon("ruler-combined"),
                         class = "btn-default btn-sm btn-block"),
            uiOutput("sim_trk_calib_report")),
          tags$hr(style = "margin:6px 0;"),
          selectInput("sim_trk_geom", "Cell geometry",
            choices = c("Sphere (coccus)" = "sphere",
                        "Rod / spherocylinder (E. coli)" = "rod",
                        "Ellipsoid (ovoid / coccobacillus)" = "ellipsoid",
                        "Membrane plane (flat sheet in 3D)" = "membrane",
                        "Free diffusion (no boundaries)" = "free"),
            selected = "rod"),
          conditionalPanel("input.sim_trk_geom == 'sphere'",
            numericInput("sim_trk_R", "Sphere radius (nm)", value = 500, min = 50, step = 50)),
          conditionalPanel("input.sim_trk_geom == 'membrane'",
            fluidRow(
              column(6, numericInput("sim_trk_mthick", "Thickness (nm)", value = 6, min = 0, step = 1)),
              column(6, numericInput("sim_trk_mtilt", "Tilt angle (°)", value = 15, min = 0, max = 89, step = 5))),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "Sheet primarily in XY, tilted by this angle (0° = flat) so z varies across the plane. ",
              "2D lateral diffusion in the (tilted) plane; z also wobbles within the thickness. ",
              "No curvature → no projection bias, unlike the curved cell wall.")),
          conditionalPanel("input.sim_trk_geom == 'free'",
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "Unbounded 3D Brownian motion — the simplest null (no geometry, no walls).")),
          conditionalPanel("input.sim_trk_geom == 'rod'",
            fluidRow(
              column(6, numericInput("sim_trk_r", "Radius (nm)", value = 400, min = 50, step = 25)),
              column(6, numericInput("sim_trk_Lcyl", "Cyl. length (nm)", value = 2000, min = 100, step = 100))),
            selectInput("sim_trk_rod_axis", "Long axis (cell orientation)",
              choices = c("X — lying in image plane" = "x",
                          "Y — lying in image plane" = "y",
                          "Z — standing (vertical)"  = "z"),
              selected = "x"),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "E. coli on a coverslip usually lie flat → long axis in XY. Z = end-on / vertical.")),
          conditionalPanel("input.sim_trk_geom == 'ellipsoid'",
            fluidRow(
              column(4, numericInput("sim_trk_ea", "a (nm)", value = 500, min = 50, step = 50)),
              column(4, numericInput("sim_trk_eb", "b (nm)", value = 500, min = 50, step = 50)),
              column(4, numericInput("sim_trk_ec", "c (nm)", value = 900, min = 50, step = 50))),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "Semi-axes along x, y, z. a=b=c → sphere; a=b<c → prolate ovoid.")),
          conditionalPanel(
            "input.sim_trk_geom == 'sphere' || input.sim_trk_geom == 'rod' || input.sim_trk_geom == 'ellipsoid'",
            radioButtons("sim_trk_mode", "Diffusion",
              choices = c("Surface (membrane / wall protein)" = "surface",
                          "Volume (intracellular / cytoplasmic)" = "volume"),
              selected = "surface"),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              tags$b("Surface"), " = on the cell boundary (2 DOF). ",
              tags$b("Volume"), " = free 3D diffusion inside, reflecting off the cell wall (3 DOF).")),
          selectInput("sim_trk_motion", "Motion type",
            choices = c("Free diffusion" = "free",
                        "Confined (OU tether → MSD plateau)" = "confined",
                        "Directed (drift velocity)" = "directed",
                        "Two-state (mobile ↔ immobile)" = "twostate"),
            selected = "free"),
          conditionalPanel("input.sim_trk_motion == 'confined'",
            numericInput("sim_trk_Lc", "Confinement radius L_c (nm)", value = 150, min = 10, step = 10),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "Harmonic (OU) tether to a home point: MSD plateaus at ~2·L_c², τ = L_c²/D. ",
              "Local D at short lags is unchanged.")),
          conditionalPanel("input.sim_trk_motion == 'directed'",
            numericInput("sim_trk_vel", "Drift velocity v (µm/s)", value = 0.5, min = 0, step = 0.1),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "Along the rod long axis (MreB-like); a fixed random direction for sphere/ellipsoid.")),
          conditionalPanel("input.sim_trk_motion == 'twostate'",
            numericInput("sim_trk_kswitch", "State-switch prob. / cycle", value = 0.03, min = 0, max = 1, step = 0.01),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "Alternates mobile (D) ↔ immobile (~0). Ground truth for HMM state recovery.")),
          checkboxInput("sim_trk_surv", "Survivorship: lose the molecule on a large dark-gap jump", FALSE),
          conditionalPanel("input.sim_trk_surv == true",
            numericInput("sim_trk_rreacq", "Re-acquisition radius (nm)", value = 300, min = 20, step = 20),
            helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
              "If the TRUE displacement across a dark gap exceeds this, the track ends — the ",
              "dark-cycle bias that manufactures apparent confinement.")),
          tags$hr(style = "margin:6px 0;"),
          numericInput("sim_trk_D", "Diffusion D (µm²/s)", value = 0.1, min = 0, step = 0.01),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
            "Free surface Brownian. Fast D (≳2–3 µm²/s) is near the MINFLUX trackable ceiling."),
          fluidRow(
            column(6, numericInput("sim_trk_ntracks", "# tracks", value = 200, min = 1, step = 10)),
            column(6, numericInput("sim_trk_meanlen", "Median locs/track", value = 60, min = 5, step = 5))),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
            "Photobleaching-limited (geometric); this sets the MEDIAN recorded length per track."),
          tags$hr(style = "margin:6px 0;"),
          fluidRow(
            column(6, numericInput("sim_trk_cycle", "Cycle time (ms)", value = 0.2, min = 0.05, step = 0.05)),
            column(6, numericInput("sim_trk_pdark", "Dark-cycle prob.", value = 0.5, min = 0, max = 0.95, step = 0.05))),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
            "Each cycle is recorded with prob (1−dark); skipped cycles make Δt an integer multiple of the cycle time."),
          fluidRow(
            column(6, numericInput("sim_trk_photons", "Mean photons/loc", value = 40, min = 5, step = 5)),
            column(6, numericInput("sim_trk_Ltcp", "Eff. TCP L (nm) → σ_loc", value = 40, min = 10, step = 5))),
          helpText(style = "font-size:0.8em; color:#666; margin-top:-4px;",
            "σ_loc = L / (2·√N_photons) per localization. L is the EFFECTIVE (final-iteration) TCP ",
            "diameter — much smaller than the initial scan pattern (≈30–60 nm gives MINFLUX's ~2–5 nm σ_loc)."),
          tags$hr(style = "margin:6px 0;"),
          checkboxInput("sim_trk_project2d", "Project to 2D (drop z)", FALSE),
          numericInput("sim_trk_seed", "Random seed", value = 1, min = 1, step = 1),
          actionButton("sim_trk_run", "Simulate tracks",
                       icon = icon("dna"), class = "btn-primary btn-block",
                       style = "margin-top:4px;"),
          tags$br(),
          fluidRow(
            column(6, downloadButton("sim_trk_dl_csv", "Tracks (CSV)", class = "btn-xs btn-default")),
            column(6, downloadButton("sim_trk_dl_gt", "Ground truth (CSV)", class = "btn-xs btn-default")))
        ),
        mainPanel(
          plotlyOutput("sim_trk_plot", height = "74vh"),
          uiOutput("sim_trk_summary")
        )
      )
    )
  )
)
