## ─────────────────────────────────────────────────────────────────────────────
## tests/run_tests.R  -  entry point for the MINFLUXR validation tests
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Recovery tests for the estimators that already live in top-level helper
##   files (helpers_tracking.R, helpers_smlm.R). Each test states a claim the
##   methods doc makes and checks it against simulated ground truth, so the claim
##   cannot silently stop being true. Estimators still nested inside the Shiny
##   server closure (CVE, MSD OLSF, JDD, HMM, change-point, D-map) cannot be
##   called from here yet; extracting them is Phase 2 of the publication plan.
##
## RUN
##   From the repo root:  Rscript tests/run_tests.R
##   Slow Monte-Carlo tests are skipped unless MINFLUXR_SLOW_TESTS=true.
## ─────────────────────────────────────────────────────────────────────────────
library(testthat)
root <- if (file.exists("global.R")) "." else ".."          # repo root or tests/
test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
