## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/fiducial_labelled.R  -  fiducial detector on labelled files
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The Fiducial / contaminant scan flags traces that behave like gold beads.
##   This runs it, through the app's own observer at default settings, on every
##   labelled acquisition available:
##     Test data/Fiducial test/Samples w beads   (3 files, beads present)
##     Test data/Fiducial test/Samples wo beads  (4 files, no beads)
##     the original pair: ...targetP3-P3img...-001 (fiducial section) vs ...-002 (clean)
##
##   A file WITH beads should flag at least one trace; a file WITHOUT should flag
##   none. Reports which rule fired (lifetime / brightness) and the evidence.
##
## RUN  Rscript tests/validation/fiducial_labelled.R   (needs Test data/)
## ─────────────────────────────────────────────────────────────────────────────
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings(source("global.R")))

files <- rbind(
  do.call(rbind, lapply(c(beads = "Test data/Fiducial test/Samples w beads",
                          clean = "Test data/Fiducial test/Samples wo beads"), function(d) {
    f <- list.files(d, "[.](npy|csv|json)$", full.names = TRUE)
    if (!length(f)) return(NULL)
    data.frame(label = if (grepl("wo beads", d)) "clean" else "beads",
               path = f, name = basename(f), stringsAsFactors = FALSE) })),
  # the original labelled pair (001 contains a fiducial section, 002 is clean)
  do.call(rbind, lapply(list(c("beads", "P3img-2uL-12perc-001"), c("clean", "2uLP3-12perc-002")), function(p) {
    f <- list.files("Test data", paste0(p[2], "[.]npy$"), full.names = TRUE)
    if (!length(f)) return(NULL)
    data.frame(label = p[1], path = f, name = basename(f), stringsAsFactors = FALSE) })))
stopifnot(!is.null(files), nrow(files) > 0)

shiny::testServer(smlm_server, {
  store <- list()
  for (i in seq_len(nrow(files))) {
    df <- tryCatch(.load_one_file(files$path[i], files$name[i]), error = function(e) NULL)
    if (is.null(df)) { cat("could not load:", files$name[i], "\n"); next }
    df$.dataset_tag <- sprintf("t%02d", i)
    store[[files$name[i]]] <- df
  }
  smlm_file_store(store)
  session$setInputs(fid_all_files = TRUE, fid_efo_z = 4, fid_dur_x = 5, fid_min_locs = 5,
                    fid_len_z = 4, fid_recur_r = 100, fid_recur_min = 4, run_fiducial = 1)
  res <- fiducial_results()
  stopifnot(!is.null(res), nrow(res) > 0)
  res$label <- files$label[match(res$dataset, files$name)]

  out <- do.call(rbind, lapply(unique(res$dataset), function(nm) {
    r <- res[res$dataset == nm, ]
    data.frame(label = r$label[1], dataset = substr(nm, 1, 40), n_traces = nrow(r),
               flagged = sum(r$flagged),
               by_lifetime = sum(r$rule %in% c("long-lived", "long+bright")),
               by_brightness = sum(r$rule == "bright"),
               max_dur_x_q99 = round(max(r$dur_x_q99, na.rm = TRUE), 1),
               max_dur_s = round(max(r$duration, na.rm = TRUE), 1),
               stringsAsFactors = FALSE)
  }))
  print(out[order(out$label, out$dataset), ], row.names = FALSE)

  fl <- res[res$flagged, ]
  if (nrow(fl)) {
    fl <- fl[order(fl$label, fl$dataset, -fl$duration), ]
    cat("\nFlagged traces:\n")
    for (i in seq_len(nrow(fl)))
      cat(sprintf("  %-5s %-30s tid %-9s %5d locs %7.2f s (%5.1fx q99) spread %7.1f nm efo %5.0f kHz [%s]\n",
                  fl$label[i], substr(fl$dataset[i], 1, 30), as.character(fl$tid[i]), fl$n_loc[i],
                  fl$duration[i], fl$dur_x_q99[i], fl$spread[i], fl$efo_med[i] / 1000, fl$rule[i]))
  }
  cat(sprintf("\nDETECTED in %d/%d bead files | FALSE POSITIVES in %d/%d bead-free files (%d traces of %d)\n",
              sum(out$flagged[out$label == "beads"] > 0), sum(out$label == "beads"),
              sum(out$flagged[out$label == "clean"] > 0), sum(out$label == "clean"),
              sum(out$flagged[out$label == "clean"]), sum(out$n_traces[out$label == "clean"])))
  utils::write.csv(out, file.path(root, "tests", "validation", "results", "fiducial_labelled.csv"), row.names = FALSE)
})
