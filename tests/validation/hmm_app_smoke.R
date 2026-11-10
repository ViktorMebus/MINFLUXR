# The repository root, two folders above this script, so it runs from any machine.
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings(source("global.R")))
# The test data is local and not part of the repository. Point MINFLUXR_TRACK_PATTERN at a
# filename pattern under "Test data/" to run this on a different tracking acquisition.
f <- normalizePath(list.files("Test data", Sys.getenv("MINFLUXR_TRACK_PATTERN", "track.*_filtered.*csv$"), full.names = TRUE)[1])
# keep the app test quick: a subset of tracks
stopifnot(!is.na(f))
cat("file:", basename(f), "
")
d <- utils::read.csv(f); keep <- utils::head(unique(d$tid), 150); tmp <- file.path(tempdir(), "track_sub.csv")
utils::write.csv(d[d$tid %in% keep, ], tmp, row.names = FALSE)
shiny::testServer(tracking_server, {
  # A real session runs once at startup before any upload. The upload observer
  # uses ignoreInit = TRUE, so without this pass it would treat the upload as
  # its initial run and skip it.
  session$flushReact()
  session$setInputs(trk_file = data.frame(name = "track_sub.csv", size = file.size(tmp),
                                          type = "text/csv", datapath = tmp, stringsAsFactors = FALSE))
  cat("rows in trk():", tryCatch(nrow(trk()), error = function(e) paste("ERR", conditionMessage(e))), "\n")
  session$setInputs(trk_msd_cycle_ms = 0.1914); cat("cycle:", tryCatch(trk_cycle_time(), error = function(e) conditionMessage(e)), "\n")
  session$setInputs(trk_sigma_loc_nm = 4.2, trk_msd_cycle_ms = 0.1914)
  session$setInputs(trk_show_hmm = TRUE, trk_hmm_n_states = 2, trk_hmm_min_pts = 10, trk_hmm_n_ex = 6, trk_hmm_boot = 0,
                     trk_hmm_max_steps = 150000, trk_hmm_bic = FALSE)
  session$setInputs(trk_hmm_run = 1)
  bd <- tryCatch(hmm_bundle(), error = function(e) { cat("HMM ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(bd)) {
    cat("D um2/s:", signif(bd$D_vals / 1e3, 3), "| occupancy:", round(bd$fracs, 3), "| dwell ms:", signif(bd$dwell_ms, 3), "\n")
    cat("sigma mode:", bd$sigma_mode, "| BIC others:", bd$bic_other, "| tracks in seg_list:", length(bd$seg_list), "\n")
    print(bd$resolv)
    for (o in c("trk_hmm_params_plot", "trk_hmm_trans_plot", "trk_hmm_example_plot", "trk_hmm_summary"))
      cat(o, ":", tryCatch({ x <- output[[o]]; "ok" }, error = function(e) paste("ERR", conditionMessage(e))), "\n")
    session$setInputs(trk_show_seg_msd = TRUE, trk_seg_msd_run = 1)
    cat("per-state MSD:", tryCatch({ r <- seg_msd_rv(); "ok" }, error = function(e) paste("ERR", conditionMessage(e))), "\n")
  }
})
