## ─────────────────────────────────────────────────────────────────────────────
## app.R  —  MINFLUXR: combined SMLM, tracking and simulator app
##
## Copyright (C) 2026 Viktor Mebus and contributors
##
## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version. It is distributed WITHOUT ANY WARRANTY; without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License (LICENSE) for more details.
##
## GPL-3 is required rather than chosen: the time-signature analysis and MSD/OLSF
## workflow are ports of MinSpt (GPL-3.0). DCR pooling and the per-trace spread
## statistic are reimplemented from pyMINFLUX (Apache-2.0), which GPL-3 permits.
## See THIRD-PARTY.md.
##
## File layout:
##   global.R            Auto-sourced by Shiny; loads packages and all helpers
##   packages.R          Package install + load, upload limit, global constants
##   helpers_shared.R    Helper functions used by both SMLM and tracking
##   helpers_smlm.R      SMLM-only helpers
##   helpers_tracking.R  Tracking-only helper functions
##   ui_smlm.R           smlm_ui  — SMLM fluidPage interior
##   ui_tracking.R       tracking_ui — Tracking fluidPage interior
##   ui.R                ui  — navbarPage wrapping both UIs
##   server_smlm.R       smlm_server()   — SMLM reactive logic
##   server_tracking.R   tracking_server() — Tracking reactive logic
##   server.R            server()  — delegates to both sub-servers
## ─────────────────────────────────────────────────────────────────────────────

# Source everything through global.R (single code path shared with the
# server.R / ui.R entry point used by RStudio's "Run App" button).
# RStudio sets the working directory to the app folder before running app.R.

# Plugins are supported using following template: add a folder containing 
# desired plugin to the /plugins/ folder. Follow the guide from the plugin_readme.md
# to correctly assemble the plugin to the rest of the app.

source("global.R")

shinyApp(ui, server)
