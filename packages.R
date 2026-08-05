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
  "magick", "fields", "base64enc", "geometry", "alphashape3d", "pracma", "Gmedian", "minpack.lm"
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

# ── Tracking chunking defaults ─────────────────────────────────────────────────
# Fall back to fst if qs could not be installed (common on macOS without
# a working C++ toolchain / OpenMP support).
TIDS_PER_CHUNK_DEFAULT <- 100L
CHUNK_FORMAT <- if (requireNamespace("qs", quietly = TRUE)) "qs" else {
  message("[MINFLUXR] 'qs' unavailable — using 'fst' for chunk storage.")
  "fst"
}
