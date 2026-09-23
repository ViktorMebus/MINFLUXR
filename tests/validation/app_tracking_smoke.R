## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/app_tracking_smoke.R  -  Tracking-tab smoke test on real data
## ─────────────────────────────────────────────────────────────────────────────
## Loads a the tracking test data subset through the Tracking tab's own upload path in a test
## session, runs every panel whose estimator lives in R/est_*.R (CVE, MSD, JDD,
## change-point, D-map, HMM, per-state MSD) through its Run button, and renders
## each output. Prints the headline numbers so a refactor that changes results
## is visible. Needs Test data/ (not in the repo).
## RUN  Rscript tests/validation/app_tracking_smoke.R
## ─────────────────────────────────────────────────────────────────────────────
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings(source("global.R")))
# The test data is local and not part of the repository. Point MINFLUXR_TRACK_PATTERN at a
# filename pattern under "Test data/" to run this on a different tracking acquisition.
f <- normalizePath(list.files("Test data", Sys.getenv("MINFLUXR_TRACK_PATTERN", "track.*_filtered.*csv$"), full.names = TRUE)[1])
stopifnot(!is.na(f))
cat("file:", basename(f), "
")
d <- utils::read.csv(f); keep <- unique(d$tid)[1:150]; tmp <- file.path(tempdir(), "track_sub.csv")
utils::write.csv(d[d$tid %in% keep, ], tmp, row.names = FALSE)

try_out <- function(output, ids) for (o in ids)
  cat(sprintf("  %-28s %s\n", o, tryCatch({ output[[o]]; "ok" }, error = function(e) paste("ERR:", conditionMessage(e)))))

shiny::testServer(tracking_server, {
  session$setInputs(trk_file = data.frame(name = "track_sub.csv", size = file.size(tmp), type = "text/csv",
                                          datapath = tmp, stringsAsFactors = FALSE))
  session$setInputs(trk_msd_cycle_ms = 0.1914, trk_L_tcp = 100.8, trk_cve_run = 0)   # 0 = unclicked, as at app start
  cat("rows:", nrow(trk()), " cycle:", trk_cycle_time(), "\n")

  session$setInputs(trk_cve_run = 1)
  b <- tryCatch(trk_cve_bundle(), error = function(e) { cat("CVE ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(b)) cat(sprintf("CVE: D %.4g um2/s, sigma %.3g nm, R %.3f, photon sigma %.3g nm\n",
                               b$D_cve / 1e3, b$sigma, b$R_blur, b$sig_phot_med))
  session$flushReact()                                  # let the observer publish the photon σ
  try_out(output, "trk_cve_status")
  cat("sigma source:", trk_sigma_loc_src(), trk_sigma_loc_eff(), "\n")

  session$setInputs(trk_show_jdd = TRUE, trk_jdd_min_pts = 5, trk_jdd_n_comp = 1, trk_jdd_step_mode = "consec",
                    trk_jdd_bin_width = 1, trk_jdd_subsamp_n = 4, trk_jdd_norm_dt = TRUE,
                    trk_jdd_integrated_only = TRUE, trk_jdd_min_dt = 0, trk_jdd_max_dt = 0, trk_jdd_sigma_loc = 0)
  session$setInputs(trk_jdd_run = 1)
  j <- tryCatch(jdd_bundle(), error = function(e) { cat("JDD ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(j)) cat(sprintf("JDD (1 comp, single-cycle): D %.4g um2/s, jumps %d\n", j$D_vals / 1e3, j$n_jumps))
  try_out(output, c("trk_jdd_plot", "trk_jdd_estimates"))

  session$setInputs(trk_show_dmap = TRUE, trk_dmap_bins = 10, trk_dmap_min_pts = 5,
                    trk_dmap_integrated_only = TRUE, trk_dmap_sub_noise = TRUE)
  session$setInputs(trk_dmap_run = 1)
  m <- tryCatch(dmap_bundle(), error = function(e) { cat("D-map ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(m)) { w <- m$N_mat > 0
    cat(sprintf("D-map: step-weighted mean D %.4g um2/s over %d cells\n",
                weighted.mean(m$D_mat[w], m$N_mat[w], na.rm = TRUE) / 1e3, sum(w))) }
  try_out(output, c("trk_dmap_plot", "trk_dmap_summary"))

  session$setInputs(trk_show_cpd = TRUE, trk_cpd_min_seg = 8, trk_cpd_max_g = 4, trk_cpd_sub_noise = TRUE)
  session$setInputs(trk_cpd_compute = 1)
  s <- tryCatch(cpd_segments_rv(), error = function(e) { cat("CPD ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(s)) cat(sprintf("CPD: %d segments, step-weighted D %.4g um2/s\n", nrow(s),
                               weighted.mean(s$D_um2_s, s$n_steps)))
  try_out(output, c("trk_cpd_hist", "trk_cpd_summary"))

  session$setInputs(trk_show_msd = TRUE, trk_msd_unified = TRUE, trk_uni_min_lag = 4, trk_uni_max_lag = 50,
                    trk_uni_linear_tau = TRUE, trk_uni_iterative = FALSE, trk_uni_model = "linear")
  session$setInputs(trk_msd_compute = 1)
  ub <- tryCatch(msd_unified_bundle(), error = function(e) { cat("MSD ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ub)) { fr <- msd_unified_fit(ub)
    if (!is.null(fr)) cat(sprintf("MSD OLSF: D %.4g um2/s, sigma %.3g nm\n", fr$D / 1e3, fr$sigma_loc)) }
  try_out(output, c("trk_msd_unified_plot"))

  session$setInputs(trk_show_hmm = TRUE, trk_hmm_n_states = 2, trk_hmm_min_pts = 10, trk_hmm_n_ex = 6, trk_hmm_boot = 0,
                     trk_hmm_max_steps = 150000, trk_hmm_bic = FALSE)
  session$setInputs(trk_hmm_run = 1)
  h <- tryCatch(hmm_bundle(), error = function(e) { cat("HMM ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(h)) cat(sprintf("HMM: D %s um2/s, dwell %s ms, sigma mode: %s\n",
                               paste(signif(h$D_vals / 1e3, 3), collapse = "/"),
                               paste(signif(h$dwell_ms, 3), collapse = "/"), h$sigma_mode))
  try_out(output, c("trk_hmm_params_plot", "trk_hmm_trans_plot", "trk_hmm_example_plot", "trk_hmm_summary"))
  session$setInputs(trk_show_seg_msd = TRUE, trk_seg_msd_run = 1)
  cat("per-state MSD:", tryCatch({ seg_msd_rv(); "ok" }, error = function(e) paste("ERR", conditionMessage(e))), "\n")
})
