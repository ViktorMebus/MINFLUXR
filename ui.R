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
    # GPL-3 "Appropriate Legal Notice" for an interactive program: who owns it,
    # under what licence, and where the source is. See LICENSE and THIRD-PARTY.md.
    footer = tags$div(
      style = paste0("margin:18px 12px 8px; padding-top:8px; border-top:1px solid #e3e3e3;",
                     "font-size:0.78em; color:#888;"),
      HTML(paste0("MINFLUXR &copy; 2026 Viktor Mebus and contributors. ",
                  "Free software under the <a href='https://www.gnu.org/licenses/gpl-3.0.html' ",
                  "target='_blank' rel='noopener'>GNU GPL v3</a>, with ABSOLUTELY NO WARRANTY. ",
                  "Source: <a href='https://github.com/ViktorMebus/MINFLUXR' target='_blank' ",
                  "rel='noopener'>github.com/ViktorMebus/MINFLUXR</a>. ",
                  "Time-signature analysis and the MSD/OLSF workflow are ported from MinSpt ",
                  "(Vogler et al. 2025, GPL-3.0); DCR pooling and per-trace spread are ",
                  "reimplemented from pyMINFLUX (Apache-2.0)."))),

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
