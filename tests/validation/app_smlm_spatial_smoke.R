## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/app_smlm_spatial_smoke.R  -  SMLM spatial panels on real data
## ─────────────────────────────────────────────────────────────────────────────
## Loads two membrane-protein acquisitions into the SMLM tab in a test session, runs Ripley
## L(r) and g(r) with the CSR envelope on, and checks that the curves, the global
## band and the global p-value are produced and that the figures build.
## RUN  Rscript tests/validation/app_smlm_spatial_smoke.R   (needs Test data/)
## ─────────────────────────────────────────────────────────────────────────────
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings(source("global.R")))
fs <- list.files("Test data/collaborator membrane-protein project", "[.]npy$", full.names = TRUE)[1:2]

shiny::testServer(smlm_server, {
  store <- list()
  for (i in seq_along(fs)) {
    df <- .load_one_file(fs[i], basename(fs[i])); df$.dataset_tag <- sprintf("t%d", i)
    store[[basename(fs[i])]] <- df
  }
  smlm_file_store(store)
  session$setInputs(active_file = names(store)[1],
                    spatial_edge_correct = TRUE, spatial_use_support = TRUE,
                    dens_grid_h = 100, dens_support_mode = "grid", dens_support_from = "units",
                    ripley_all_files = TRUE, ripley_envelope = TRUE, ripley_nsim = 19,
                    ripley_rmax = 1000, ripley_steps = 40, ripley_unit = "tid",
                    pcf_all_files = TRUE, pcf_envelope = TRUE, pcf_nsim = 19,
                    pcf_rmax = 1000, pcf_steps = 40, pcf_unit = "tid")
  session$setInputs(run_ripley = 1)
  env <- ripley_env_results()
  cat("Ripley: curves", nrow(ripley_results() %||% data.frame()),
      "| envelope rows", nrow(env %||% data.frame()),
      "| global p:", if (!is.null(env)) paste(signif(unique(env$global_p), 3), collapse = ", ") else "-", "\n")
  cat("  global band present:", !is.null(env) && all(c("g_lo", "g_hi") %in% names(env)), "\n")
  cat("  plot:", tryCatch({ build_ripley_plot_fig(); "ok" }, error = function(e) paste("ERR", conditionMessage(e))), "\n")
  session$setInputs(run_pcf = 1)
  env2 <- pcf_env_results()
  cat("g(r):   curves", nrow(pcf_results() %||% data.frame()),
      "| envelope rows", nrow(env2 %||% data.frame()),
      "| global p:", if (!is.null(env2)) paste(signif(unique(env2$global_p), 3), collapse = ", ") else "-", "\n")
  cat("  plot:", tryCatch({ build_pcf_plot_fig(); "ok" }, error = function(e) paste("ERR", conditionMessage(e))), "\n")
  cat("  window readout:", tryCatch({ output$spatial_window_readout; "ok" }, error = function(e) paste("ERR", conditionMessage(e))), "\n")
})
