## ─────────────────────────────────────────────────────────────────────────────
## global.R  —  Auto-sourced by Shiny before server.R / ui.R
##
## Shiny sets the working directory to the app folder before sourcing this file,
## so all relative source() calls below are safe on any OS / any user path.
##
## This ensures packages are installed and the upload limit is raised whether
## the app is launched via app.R or directly via server.R / ui.R (e.g. the
## RStudio "Run App" button or shiny::runApp()).
## ─────────────────────────────────────────────────────────────────────────────

source("packages.R")
source("helpers_shared.R")
source("helpers_smlm.R")
source("helpers_tracking.R")
# R/ holds the estimators as plain functions, outside the Shiny server closures, so
# they can be called (and tested) without a running app. The server panels gather
# the inputs and call these; see tests/ for the recovery tests that use them.
source("R/est_blur.R")       # motion-blur coefficient R from the sequence timing
source("R/est_cve.R")        # CVE: D and sigma_loc from step covariance, photon-budget sigma
source("R/est_msd.R")        # MSD fits: OLSF (linear) + confined / directed / anomalous
source("R/est_cpd.R")        # change-point segmentation of a track + per-segment D
source("R/est_jdd.R")        # jump-distance distribution mixture fit
source("R/est_dmap.R")       # spatially resolved D-map
source("R/est_hmm.R")        # continuous-time diffusive-state HMM
source("ui_smlm.R")
source("ui_tracking.R")
source("sim_ui.R")

## ── Plugin registry ───────────────────────────────────────────────────────────
## Plugins live in plugins/<name>/plugin.R.
## Each plugin.R calls register_plugin(tabs = list(...), server_fn = fn).
## Add a plugin: drop a folder in plugins/.  Remove: delete the folder.

.plugin_tabs    <- list()
.plugin_servers <- list()

register_plugin <- function(tabs = list(), server_fn = NULL) {
  .plugin_tabs    <<- c(.plugin_tabs, if (is.list(tabs)) tabs else list(tabs))
  if (is.function(server_fn))
    .plugin_servers[[length(.plugin_servers) + 1L]] <<- server_fn
}

if (dir.exists("plugins")) {
  for (.plugin_dir in list.dirs("plugins", recursive = FALSE, full.names = TRUE)) {
    .entry <- file.path(.plugin_dir, "plugin.R")
    if (file.exists(.entry)) source(.entry)
  }
  rm(list = intersect(c(".plugin_dir", ".entry"), ls()))
}

source("ui.R")
source("server_smlm.R")
source("server_tracking.R")
source("sim_server.R")
source("server.R")
