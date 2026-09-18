# Source the pure helper files once for every test file. Located relative to this
# file so the tests run from the repo root or from tests/.
.repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages({
  source(file.path(.repo_root, "helpers_shared.R"))
  source(file.path(.repo_root, "helpers_tracking.R"))
  source(file.path(.repo_root, "helpers_smlm.R"))
})
# Slow Monte-Carlo recovery tests run only when asked for.
skip_if_not_slow <- function()
  testthat::skip_if_not(identical(Sys.getenv("MINFLUXR_SLOW_TESTS"), "true"),
                        "slow test (set MINFLUXR_SLOW_TESTS=true)")
