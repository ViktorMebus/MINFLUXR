## ─────────────────────────────────────────────────────────────────────────────
## plugin.R  —  MINFLUXR plugin template (starting point for your own plugin)
##
## WHAT THIS IS
##   A minimal, self-contained example plugin. Copy this whole folder into
##   `plugins/<yourname>/`, rename things, and you have a working extra tab in
##   the MINFLUXR navbar — without touching any of the core app files.
##
## HOW PLUGINS ARE LOADED
##   At startup, global.R scans the `plugins/` folder (note the plural) and
##   sources every `plugins/*/plugin.R` it finds. Before sourcing yours it sets
##   `.plugin_dir` to your folder's path, so relative source() calls are safe.
##   Each plugin.R must finish by calling register_plugin(...) (see bottom).
##
##   This template lives in `plugin/` (singular) so it is NOT auto-loaded — it
##   is just the copy-me source. Nothing here runs until you place a copy under
##   `plugins/`.
##
## ENABLE  : copy this folder to  plugins/<yourname>/
## DISABLE : delete or rename that plugins/<yourname>/ folder
##
## Full walk-through + the `smlm` data reference:  see  plugin_readme.md
## ─────────────────────────────────────────────────────────────────────────────

pd <- .plugin_dir   # path to THIS plugin's folder, set by the loader in global.R

## ── (Optional) split your code across files ──────────────────────────────────
## For anything beyond a small tab, keep plugin.R as the thin entry point and put
## the real code in sibling files, sourced here. Create the files, then uncomment:
#
#   source(file.path(pd, "packages.R"))   # install/load extra CRAN packages
#   source(file.path(pd, "helpers.R"))    # pure helper functions
#   source(file.path(pd, "ui.R"))         # defines your tabPanel(s)
#   source(file.path(pd, "server.R"))     # defines your server function
#
## The example below is inlined so this single file works on its own.


## ── 1. UI: one tab in the MINFLUXR navbar ────────────────────────────────────
## A plugin tab is just a Shiny tabPanel(). Give `value` a unique id so it does
## not clash with the core tabs (tab_smlm, tab_tracking, tab_sim) or other plugins.

my_ui_tab <- tabPanel(
  title = "My Plugin",              # label shown in the navbar
  value = "tab_myplugin",           # unique internal id

  fluidPage(
    titlePanel("My Plugin"),
    sidebarLayout(
      sidebarPanel(
        helpText("This is the starting template. Replace it with your own UI."),
        actionButton("myplugin_go", "Summarise loaded SMLM data")
      ),
      mainPanel(
        verbatimTextOutput("myplugin_summary")
      )
    )
  )
)


## ── 2. Server: reactive logic for the tab ────────────────────────────────────
## The loader calls your server function with the live Shiny session AND `smlm`,
## a list of reactive handles exported by the core SMLM tab. That is how a plugin
## reads whatever the user has currently loaded/filtered. See plugin_readme.md
## for the full list of what `smlm` gives you.
##
## Namespacing note: this template uses the SAME input/output namespace as the
## core app (there is no shiny::moduleServer wrapper). Prefix every input/output
## id you create with your plugin name (e.g. "myplugin_...") so ids stay unique.

my_server <- function(input, output, session, smlm) {

  output$myplugin_summary <- renderPrint({
    input$myplugin_go                       # re-run when the button is clicked
    df <- smlm$plottable_data()             # current, filtered localizations
    shiny::validate(
      shiny::need(is.data.frame(df) && nrow(df) > 0,
                  "Load a file in the SMLM Analysis tab first.")
    )
    cat("Localizations currently shown:", nrow(df), "\n")
    cat("Columns:", paste(names(df), collapse = ", "), "\n")
    if (all(c("x", "y") %in% names(df))) {
      cat(sprintf("X range: %.1f  to  %.1f nm\n", min(df$x), max(df$x)))
      cat(sprintf("Y range: %.1f  to  %.1f nm\n", min(df$y), max(df$y)))
    }
  })

}


## ── 3. Register with the app ─────────────────────────────────────────────────
## This is the one required call. `tabs` are appended to the navbar; `server_fn`
## is run inside the main Shiny session. Pass a list of tabs if you add several.

register_plugin(
  tabs      = list(my_ui_tab),
  server_fn = my_server
)
