## ─────────────────────────────────────────────────────────────────────────────
## packages.R  —  unified package install + load for MINFLUXR combined app
## ─────────────────────────────────────────────────────────────────────────────
options(repos = c(CRAN = "https://cloud.r-project.org"))

# rgl (required by alphashape3d) needs XQuartz on macOS, which many users don't
# have.  Setting useNULL = TRUE makes rgl run headlessly — safe here because all
# 3-D visualisation is done via Plotly; rgl/alphashape3d are used for geometry
# computation only.
options(rgl.useNULL = TRUE)

# Shiny upload limit — must be set before any session starts.
# Default is 5 MB; we raise it to 2 GB.  This is also set in global.R so it
# applies regardless of whether the app is launched via app.R or server.R/ui.R.
options(
  scipen              = 999,
  digits              = 6,
  shiny.maxRequestSize = 2 * 1024^3   # 2 GB
)

# SMLM analysis packages
pkgs_smlm <- c(
  "shiny", "shinyWidgets", "shinyjs",
  "colourpicker",
  "plotly", "RColorBrewer",
  "dplyr", "rlang",
  "dbscan", "RANN",
  "MASS", "mclust", "openxlsx", "viridisLite",
  "magick", "fields", "base64enc", "geometry", "alphashape3d", "pracma", "Gmedian", "minpack.lm",
  # misc3d: marching cubes for the Median-Cluster isosurface (computeContour3d).
  # Matrix: the sparse neighbour operator the mesh smoothing uses. Matrix ships
  # with R, but it is named here so the dependency is visible rather than assumed.
  # digest: keys the isosurface mesh cache on its geometry. Ships with Shiny's
  # dependency tree, but named here so the dependency is visible rather than assumed.
  "misc3d", "Matrix", "digest"
)

# Tracking-specific extra packages
pkgs_tracking <- c(
  "data.table", "readr", "fst", "qs", "zip", "later", "changepoint"
)

# jsonlite: must be INSTALLED but NOT attached via library().
# jsonlite exports validate(txt) which shadows shiny::validate(), breaking all
# shiny::validate(shiny::need(...)) calls with a cryptic "is.character(txt) is
# not TRUE" error.  All code uses jsonlite::fromJSON() with explicit namespace.
pkgs_install_only <- c("jsonlite")

pkgs_all <- unique(c(pkgs_smlm, pkgs_tracking, pkgs_install_only))

# ── Robust installer ───────────────────────────────────────────────────────────
# Installs each missing package individually so a single compile failure does not
# abort the rest.  On macOS, a second attempt using pre-built binaries is made
# for packages that fail from source (qs and fst use C++/OpenMP and often fail
# to compile on macOS without Xcode command-line tools or gfortran).
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))

  message("\n[MINFLUXR] Installing missing packages: ",
          paste(missing, collapse = ", "), "\n")

  is_mac <- identical(Sys.info()[["sysname"]], "Darwin")
  failed  <- character(0)

  for (pkg in missing) {
    installed <- tryCatch({
      install.packages(pkg, dependencies = TRUE, quiet = FALSE)
      requireNamespace(pkg, quietly = TRUE)
    }, error   = function(e) FALSE,
       warning = function(w) requireNamespace(pkg, quietly = TRUE))

    # macOS fallback: try a pre-built binary when source compilation fails
    if (!installed && is_mac) {
      message("[MINFLUXR] Source install of '", pkg,
              "' failed — retrying with binary build (macOS)...")
      installed <- tryCatch({
        install.packages(pkg, dependencies = TRUE, type = "binary", quiet = FALSE)
        requireNamespace(pkg, quietly = TRUE)
      }, error   = function(e) FALSE,
         warning = function(w) requireNamespace(pkg, quietly = TRUE))
    }

    if (!installed) {
      warning("[MINFLUXR] Could not install '", pkg,
              "'. Features that depend on it may be unavailable.", call. = FALSE)
      failed <- c(failed, pkg)
    } else {
      message("[MINFLUXR] '", pkg, "' ready.")
    }
  }

  if (length(failed)) {
    message("\n[MINFLUXR] WARNING — the following packages could not be",
            " installed and some features will be disabled:\n  ",
            paste(failed, collapse = ", "),
            "\nOn macOS, make sure Xcode Command Line Tools are installed",
            " (run: xcode-select --install) and restart R before trying again.\n")
  }

  invisible(!length(failed))
}

install_if_missing(pkgs_all)

# ── Load packages (skip those that failed to install) ─────────────────────────
pkgs_to_load <- setdiff(pkgs_all, pkgs_install_only)
suppressPackageStartupMessages({
  for (pkg in pkgs_to_load) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning("[MINFLUXR] Skipping library('", pkg,
              "') — package not available.", call. = FALSE)
      next
    }
    library(pkg, character.only = TRUE)
  }
})

# ── SMLM trace-spread (σ̃) defaults ─────────────────────────────────────────────
# Minimum number of localizations a trace must have before its localization
# spread (per-axis SD) is computed. Traces shorter than this are excluded from
# the σ̃ histogram, the median-spread labels, the allowable-spread ranges, and
# both spread filters. Literature (e.g. pyMINFLUX) typically uses ≥5 points.
MIN_LOCS_FOR_SPREAD_DEFAULT <- 5L

# ── Axial (z) calibration default ──────────────────────────────────────────────
# Refractive-index mismatch correction applied to z. 0.7 is the pyMINFLUX
# default and only a starting point: the physically correct value depends on the
# immersion/sample index pair and imaging depth, so it belongs to the DATASET,
# not the software. Both tabs expose it and neither hard-codes it any longer.
TRK_Z_SCALE_DEFAULT <- 0.7

# ── Trackable-diffusion ceiling (D_max) default ────────────────────────────────
# Diameter of the FINAL targeted-coordinate-pattern iteration - the one a tracking
# sequence reiterates continuously once locked on (headstart = -1). That is the
# pattern the particle must stay inside between consecutive fixes, so it is the L
# that sets D_max; the coarse search iterations before lock-on do not enter, and
# neither does the much tighter final iteration an IMAGING sequence uses. SPT stops
# at a larger L on purpose: precision improves only linearly as L shrinks
# (σ_CRLB ≈ L/(4√N)) while the escape risk grows quadratically (D_max ∝ L²).
#
# 100 nm is the round value of the 100.8 nm (patGeoFactor 0.28) final iteration
# shared by every stock tracking sequence to hand - Tracking_2D_Oct2022,
# Tracking_2D_Fast_Oct2022 and seqTrkOctahedron_3D alike, 2D and 3D. Vogler et al.
# 2025 quote 150 nm because their sequence stops one iteration earlier, at
# patGeoFactor 0.42; the imaging sequences run down to 39.6 nm, which is why an
# imaging-derived L reads far too low here.
#
# This is only the fallback. When a sequence .json is uploaded the exact L is read
# from its patGeoFactor and fills the input, which is always preferable: the ceiling
# scales as L², so L belongs to the SEQUENCE, not to the software.
# Read by ui_tracking.R's input default and by dmax_estimate()'s fallback, so the
# two cannot drift apart.
DMAX_L_TCP_DEFAULT <- 100

# Base TCP diameter (nm) at patGeoFactor = 1. Imspector sequence files do not store
# the pattern diameter directly; they store `patGeoFactor`, a dimensionless
# multiplier, and L = patGeoFactor x this constant. It is wavelength-independent -
# the same factor gives the same L at 642 and 680 nm excitation.
#
# Validated against three independent references:
#   patGeoFactor 0.79 -> 284.4 nm  and  0.84 -> 302.4 nm, the two search-iteration
#     diameters reported for MINFLUX SPT in Vogler et al. 2025 (284 and 302 nm)
#   patGeoFactor 0.11 -> 39.6 nm, the final imaging iteration, matching the ~40 nm
#     final-iteration diameter quoted in the imaging literature
#   patGeoFactor 0.42 -> 151.2 nm, the tracking iteration Vogler quote as 150 nm
TCP_L_PER_GEOFACTOR_NM <- 360

# ── MINFLUX localization precision from the photon budget ──────────────────────
# sigma_loc = L_TCP / (TCP_CRLB_PREFACTOR * sqrt(N)), N = detected photons (eco).
# The Cramer-Rao bound for a doughnut TCP of diameter L scales as L/sqrt(N); the
# prefactor depends on pattern geometry and is 4 for the hexagonal TCP used for
# tracking. Cross-checked on the GP67H files: L = 100.8 nm with a median eco of
# 36-38 photons gives 4.09 / 4.20 nm, and inverting the blur-aware CVE sigma the
# other way (L = 4*sigma*sqrt(N)) returns 101.9 / 98.0 nm against the sequence's
# 100.8 nm. Shared by the tracking sigma_loc estimator and the tracking simulator
# so the forward and inverse models cannot disagree.
#
# NOT the SMLM simulator's constant. That one uses an EFFECTIVE L in
# sigma = L_eff/(2*sqrt(N)), where L_eff absorbs drift, linker motion and imperfect
# centring and so runs several times the physical L. Different quantity, different
# prefactor - see docs/simulator.md.
TCP_CRLB_PREFACTOR <- 4

# ── MSD plot height (px) ───────────────────────────────────────────────────────
# Change this one number to make the MSD curve taller or shorter. The server reads
# it too: the legend sits BELOW the axes at a paper-coordinate offset, which is a
# FRACTION of the plotting area, so a taller plot would push the legend past the
# fixed bottom margin and clip it. Deriving the offset from this constant keeps the
# legend the same distance below the axes at any height.
TRK_MSD_PLOT_PX <- 460
TRK_MSD_MARGIN_T <- 10    # px, top margin used by the MSD layout
TRK_MSD_MARGIN_B <- 110   # px, bottom margin reserved for the legend
TRK_MSD_LEGEND_GAP <- 61  # px, gap from the x-axis down to the legend

# ── Tracking chunking defaults ─────────────────────────────────────────────────
# Fall back to fst if qs could not be installed (common on macOS without
# a working C++ toolchain / OpenMP support).
TIDS_PER_CHUNK_DEFAULT <- 100L
CHUNK_FORMAT <- if (requireNamespace("qs", quietly = TRUE)) "qs" else {
  message("[MINFLUXR] 'qs' unavailable — using 'fst' for chunk storage.")
  "fst"
}
