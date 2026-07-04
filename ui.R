## ─────────────────────────────────────────────────────────────────────────────
## ui.R  —  Combined navbarPage wrapping SMLM, Tracking and Simulation tabs
## Requires: smlm_ui, tracking_ui, sim_ui
## ─────────────────────────────────────────────────────────────────────────────

## Core tabs are listed explicitly; plugin tabs are appended automatically from
## .plugin_tabs, which is populated by register_plugin() calls in global.R.

ui <- do.call(navbarPage, c(
  list(
    title = "MINFLUXR",
    id    = "main_navbar",

    tabPanel(
      title = "SMLM Analysis",
      value = "tab_smlm",
      smlm_ui
    ),

    tabPanel(
      title = "Tracking Analysis",
      value = "tab_tracking",
      tracking_ui
    ),

    tabPanel(
      title = "MINFLUX Simulator",
      value = "tab_sim",
      sim_ui
    )
  ),
  .plugin_tabs   # zero or more tabPanels contributed by plugins
))
