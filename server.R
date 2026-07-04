## ─────────────────────────────────────────────────────────────────────────────
## server.R  —  Combined server: delegates to smlm_server() and tracking_server()
## Both sub-servers run within the same Shiny session so reactiveValues and
## observers are all live regardless of which tab is visible.
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Core modules ─────────────────────────────────────────────────────────────
  smlm_exports <- smlm_server(input, output, session)
  tracking_server(input, output, session)
  sim_server(input, output, session)

  # ── Plugin modules (populated by register_plugin() calls in global.R) ────────
  for (.srv in .plugin_servers) .srv(input, output, session, smlm = smlm_exports)

}
