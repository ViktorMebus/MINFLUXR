## ─────────────────────────────────────────────────────────────────────────────
## helpers_smlm.R  -  SMLM helper functions (loader)
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   All the computation behind the SMLM (imaging) tab. Unlike the Tracking tab -
##   which asks how molecules MOVE - this tab asks where molecules ARE and how
##   they are ARRANGED: clustering, shape fitting, spatial statistics and
##   two-colour comparison. Nothing here is reactive; every function is a plain
##   input -> output transform testable in the console.
##
## THE DATA HIERARCHY (understand this and the helpers make sense)
##   localizations          raw points. One BLINK produces many of them.
##      -> tid              trace ID: one continuous blink of one fluorophore
##      -> cluster          de-blinked molecule: traces of one fluorophore merged by
##                          complete linkage on their centroids (merge-only, no
##                          noise). Its job is DE-BLINKING.
##      -> cluster2         "median"/secondary clustering, run over the medians
##                          of the level above. Its job is STRUCTURE: group
##                          molecules into the biological object they form.
##   Most functions take a level as input and summarise to the level above, so
##   "which level am I on?" is the first question to ask of any function here.
##
##   WHY THIS MATTERS SCIENTIFICALLY: a fluorophore blinking 10 times looks like
##   10 molecules, and an antibody carries several dyes ~10-40 nm from its
##   target. Feeding RAW localizations into a clustering or spatial statistic
##   therefore measures the LABEL, not the biology. The tiering exists to
##   collapse that artefact before any counting or spatial statistic is applied.
##
## WHERE THINGS ARE (R/, loaded below in this order)
##
##   smlm_plot.R          Axes, colour and plot points
##   smlm_channels.R      Channels: two-colour split and multiplex
##   smlm_filters.R       Filters, trace spread and localization precision
##   smlm_clustering.R    De-blinking and median clustering
##   smlm_nn.R            Nearest-neighbour distances
##   smlm_shapes.R        Shape fitting
##   smlm_structures.R    Structures: planes, surfaces, isosurfaces, curvature
##   smlm_gt.R            Ground-truth estimate and structure size (exploratory)
##   smlm_spatial.R       Spatial statistics and co-localization
##   smlm_support.R       Footprint (support) and labelling-density normalisation
##   smlm_hotspots.R      Local density and hotspots
##   smlm_summary.R       Per-file summaries for group comparisons
##
##   R/smlm_settings.R (settings and display options) is loaded by global.R.
##
## UNITS: nanometres throughout. `tim` is seconds in the SMLM tab (the Tracking
## tab converts to ms; this tab does not).
##
## Method provenance and the MINFLUX-appropriateness of each spatial statistic
## are documented in docs/MINFLUX_analysis_methods.md (Part II §6).
## ─────────────────────────────────────────────────────────────────────────────


# The parts are sourced into the environment this file is being sourced into, so
# source("helpers_smlm.R") behaves exactly as when everything lived in one file.
# The folder is found from this file's own path, so it works from any working
# directory (the app, the tests, analysis scripts).
.smlm_helpers_env <- environment()
local({
  here <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { here <- dirname(normalizePath(of)); break }
  }
  if (is.null(here)) here <- getwd()
  parts <- c("smlm_plot.R",
             "smlm_channels.R",
             "smlm_filters.R",
             "smlm_clustering.R",
             "smlm_nn.R",
             "smlm_shapes.R",
             "smlm_structures.R",
             "smlm_gt.R",
             "smlm_spatial.R",
             "smlm_support.R",
             "smlm_hotspots.R",
             "smlm_summary.R")
  for (p in parts) source(file.path(here, "R", p), local = .smlm_helpers_env, encoding = "UTF-8")
})
rm(.smlm_helpers_env)
