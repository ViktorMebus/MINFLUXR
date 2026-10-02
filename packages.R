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
  "misc3d", "Matrix", "digest",
  # DT: the sortable data tables in the SMLM tab (DT::datatable / DTOutput). Nothing
  # else in the tree pulls it, so a clean machine fails with
  # "there is no package called 'DT'" the first time a table renders.
  "DT",
  # scales: axis/label formatting helpers used directly as scales::.
  "scales",
  # htmltools: called explicitly as htmltools:: in ~160 places across the tracking
  # and simulator servers. It arrives with shiny in practice, but the code does not
  # go through shiny to reach it, so it is named rather than assumed.
  "htmltools",
  # stringi: string handling in helpers_shared.R. Comes with stringr, which is NOT
  # otherwise a dependency here, so it has to be named.
  "stringi"
)

# Tracking-specific extra packages
pkgs_tracking <- c(
  "data.table", "readr", "fst", "zip", "later", "changepoint"
)

# Optional: nice to have, never required. `qs` needs a C++ toolchain plus
# RcppParallel and stringfish and fails to install on many machines, so the app
# writes chunks with `fst` instead. qs is only used to READ a chunk set that was
# produced on a machine that had it. Install failures here are not an error.
pkgs_optional <- c("qs")

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
# DEPENDENCY SCOPE: `dependencies = TRUE` also pulls **Suggests**, which for
# plotly/dplyr is dozens of unrelated packages (nycflights13, lobstr, ...) and, worse,
# lets the installer UPGRADE packages that are already loaded in this session.
# Replacing a loaded package's files leaves it half-removed - observed here: an
# install of plotly upgraded `rlang` underneath a running session and broke every
# package that imports it. DEP_SCOPE keeps it to what the app actually needs.
DEP_SCOPE <- c("Depends", "Imports", "LinkingTo")

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))

  # Never try to (re)install something already loaded - on Windows especially, the
  # files are locked and a failed replace corrupts the installed copy.
  loaded <- intersect(missing, loadedNamespaces())
  if (length(loaded)) {
    message("[MINFLUXR] Already loaded, not reinstalling: ", paste(loaded, collapse = ", "),
            "
           Restart R first if one of these really needs updating.")
    missing <- setdiff(missing, loaded)
    if (!length(missing)) return(invisible(TRUE))
  }

  message("\n[MINFLUXR] Installing missing packages: ",
          paste(missing, collapse = ", "), "\n")

  is_mac <- identical(Sys.info()[["sysname"]], "Darwin")
  failed  <- character(0)

  for (pkg in missing) {
    installed <- tryCatch({
      install.packages(pkg, dependencies = DEP_SCOPE, quiet = FALSE)
      requireNamespace(pkg, quietly = TRUE)
    }, error   = function(e) FALSE,
       warning = function(w) requireNamespace(pkg, quietly = TRUE))

    # macOS fallback: try a pre-built binary when source compilation fails
    if (!installed && is_mac) {
      message("[MINFLUXR] Source install of '", pkg,
              "' failed — retrying with binary build (macOS)...")
      installed <- tryCatch({
        install.packages(pkg, dependencies = DEP_SCOPE, type = "binary", quiet = FALSE)
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
# prefactor depends on pattern geometry and on the instrument. Shared by the
# tracking sigma_loc estimator and the tracking simulator so the forward and
# inverse models cannot disagree.
#
# NOT the SMLM simulator's constant. That one models a 2D imaging pattern with
# SMLM_CRLB_PREFACTOR = sqrt(8) (see below); same physical L, different geometry.
#
# The prefactor a has no system-independent value. The ideal-estimator (unbiased)
# value follows from the pattern geometry: sqrt(8) = 2.83 for any symmetric planar
# pattern, sqrt(16/3) = 2.31 for an idealised octahedron (Balzarotti 2017;
# Scheiderer, Marin & Ries 2024 Eq. 8; TCP_PATTERN_GEOM in R/est_cve.R). Exported
# positions carry the instrument's own estimator, background and tracking loop, so
# the effective a differs per instrument and per sequence and is MEASURED on an
# immobile sample (Instrument Calibration, R/est_calib.R), kept as a named profile.
# This constant is only the uncalibrated fallback for a 2D pattern, and the value
# the tracking simulator draws its noise with; the app resolves the fallback per
# pattern, and labels sigma computed with it as uncalibrated.
TCP_CRLB_PREFACTOR <- sqrt(8)

# The TCP pattern geometry (axes spanned, ideal prefactor, off-centre coefficient)
# lives with the estimator, as TCP_PATTERN_GEOM in R/est_cve.R.

# ── SMLM simulator: 2D MINFLUX precision with a photon-independent floor ──────
# Ideal per-axis Cramer-Rao bound at the centre of a 2D doughnut pattern (three
# exposures on a circle of diameter L plus one at the centre), small-L limit:
#   sigma_CRLB = L / (sqrt(8) * sqrt(N)) * g(SBR)
#   g(SBR)     = sqrt( (1 + 1/SBR) * (1 + 3/(4*SBR)) )
# N = total detected photons (signal + background). Balzarotti et al. 2017
# Science Suppl. Eq. S27/S31; re-derived as Eq. 44/45 in Marin & Ries 2025
# (arXiv:2410.12427). With the background FRACTION f = B/N, SBR = (1-f)/f and
#   g = sqrt(4 - f) / (2 * (1 - f)).
# The small-L limit drops a factor 8*sq^2/(8*sq^2 - L^2) that is ~1.03 at the
# ~40 nm final imaging L and grows only as L approaches the beam size.
SMLM_CRLB_PREFACTOR <- sqrt(8)

# GOAL: localization precision for the SMLM simulator from a PHYSICAL pattern
# diameter. Two excess terms, because they scale differently with photons:
#   k      multiplies the bound: pattern geometry (hexagonal TCP), mLMS estimator
#          loss, detector effects. Scales as 1/sqrt(N) like the bound itself.
#   floor  adds in quadrature: drift, vibration, residual motion within a trace.
#          Independent of N, so more photons cannot beat it.
# Returns the ideal bound, the photon-limited term (k x bound) and the total.
smlm_sigma_loc <- function(L_nm, n_photons, bg_frac = 0, k = 1, floor_nm = 0) {
  f     <- max(0, min(0.99, bg_frac))                 # background fraction of N
  n     <- max(1, n_photons)                          # total detected photons
  g     <- sqrt(4 - f) / (2 * (1 - f))                # SBR penalty, 1 at f = 0
  crlb  <- L_nm / (SMLM_CRLB_PREFACTOR * sqrt(n)) * g # ideal per-axis bound
  phot  <- max(0, k) * crlb                           # photon-limited term
  list(crlb = crlb, photon = phot, g = g,
       total = sqrt(phot^2 + max(0, floor_nm)^2))     # quadrature sum
}

# ── Excitation doughnut FWHM (nm), the nu of Vogler Eq. 1 ─────────────────────
# sigma = L/(4*sqrt(N)) is only the nu -> infinity limit of the doughnut CRLB. The
# full form carries a denominator |1 + ln2*(L/nu)^2*((2x/L)^2 - 1)| that grows the
# estimate as L approaches nu, plus an x-dependence on the emitter's offset from
# the TCP centre.
#
# 195.5 nm is fitted from Vogler et al. 2025 Table 1/2: their listed sigma_CRLB of
# 20.03 / 14.35 / 9.08 nm at PL 10 / 20 / 50 with a 150 nm final iteration are
# reproduced by Eq. 1 at x = 0 to within 1.3%. It implies sigma = L/(2.35*sqrt(N))
# at that L, not L/(4*sqrt(N)).
#
# CAVEAT: this is transferred from THEIR microscope. nu is set by the objective and
# wavelength, so it belongs to the instrument, not the software. It is exposed as
# an Acquisition Parameter for exactly that reason. An immobile-sample run at two or
# three different L over-determines nu and tests the model - see
# docs/sigma_loc_model.md.
# 0 = OFF, i.e. the large-beam (nu -> infinity) limit, which is what this microscope
# measures. Vogler Eq. 1's denominator predicts sigma*sqrt(N)/L rising 39% from
# L = 100.8 to 151.2 nm; the fixed-cell series measured a ratio of **1.04**. The
# nu = 195.5 nm fitted from their Table 2 does not transfer to this instrument.
# Set a finite value here only if a future measurement shows an L-dependence.
TCP_BEAM_FWHM_NM <- 0

# ── Continuous-time HMM: validated range (est_hmm.R) ──────────────────────────
# From tests/validation/hmm_recovery_sweep.R (simulated tracking acquisition: t_c 0.192
# ms, SW 0.011 ms, dark-cycle mix, photon σ, blur R from the timing model, 200 tracks × 300 locs).
# A state is flagged in the HMM summary when its fitted mean dwell is below
# HMM_RESOLV_MIN_DWELL_CYC cycles, or when fewer than HMM_RESOLV_MIN_SWITCHES
# switches out of it are expected over the whole data set. See the sweep results
# in docs/MINFLUX_analysis_methods.md §I-B.5 for where these numbers come from.
# Sweep result (2026-09-17, 54 fits): dwell 1-3 cycles fails (dwell +34 to +180%, CIs never cover);
# 10 cycles within +1 to +19%; 30-300 cycles within ±10% on average. The smallest validated case
# had ~290 expected switches, so fewer than 300 is outside what has been checked.
HMM_RESOLV_MIN_DWELL_CYC <- 10
HMM_RESOLV_MIN_SWITCHES  <- 300

# ── Tracker wander floor, D_w (nm^2/ms) ───────────────────────────────────────
# The reported position wanders even on a genuinely immobile target, so a measured
# D contains an instrumental term. Measured 2026-09-02:
#   40 nm beads on glass, L = 100.8   D_w =  36 nm^2/ms   (ideal rigid target)
#   fixed E. coli on agarose, L=100.8 D_w = 160 nm^2/ms   (real photon budget)
#   fixed E. coli on agarose, L=151.2 D_w = 408 nm^2/ms
#
# It is PHOTON-DRIVEN, not sample motion: binning fixed-cell segments by their own
# ECO gives D_w = 444 -> 159 nm^2/ms as ECO goes 25 -> 34, i.e. falling FASTER than
# 1/N (0.36 observed against 0.74 for pure 1/N). Agarose motion would be flat in ECO.
# Bead wander is also incoherent between beads, so it is not stage drift either.
#
# 0 = off. It is NOT subtracted automatically: D_w depends on the photon budget, so a
# single constant would over-correct dim tracks and under-correct bright ones. Set it
# to a value measured on YOUR fixed sample at the same L and photon budget, and read
# the corrected D as an estimate with that caveat attached.
TRK_DW_FLOOR_DEFAULT <- 0

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
# Chunks are written with `fst`, always.
#
# It used to prefer `qs` and fall back to `fst`, which made the on-disk format
# depend on whether a given machine had managed to compile qs - so a chunk set was
# not reliably portable between people. fst is the better fit anyway: it is a hard
# dependency already, needs no special toolchain, and can read a COLUMN SUBSET
# without touching the rest of the file, which the playback path uses to pull just
# the tid column. Existing .qs chunk sets are still READ, via read_chunk_file().
TIDS_PER_CHUNK_DEFAULT <- 100L
CHUNK_FORMAT <- "fst"
