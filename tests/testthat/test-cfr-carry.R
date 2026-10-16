# CFR / ECC / EFC carry-forward in .minflux_finalize_df (helpers_shared.R).
# CFR is measured in one iteration only, so each final localization must take the
# value of its OWN localization, never one value for the whole trace.

# One raw row per iteration; coordinates in metres as the readers deliver them.
.raw_rows <- function(tid, itr, cfr, fnl = NULL) {
  n <- length(itr)
  d <- data.frame(x = 1e-7 * seq_len(n), y = 1e-7 * seq_len(n), z = 0,
                  tid = tid, itr = itr, tim = seq_len(n) * 1e-3,
                  efo = 1e5, eco = 100, cfr = cfr, ecc = 0, efc = 0, dcr1 = 0.5)
  if (!is.null(fnl)) d$fnl <- fnl
  d
}

# CLAIM: 3D imaging writes itr 0..9 for the first localization and repeats
# 6,7,8,9 for each further one; each itr-9 row carries the itr-6 CFR of its own
# localization (itr 4's coarse CFR is overridden by the later itr 6).
test_that("3D imaging: each final localization gets its own itr-6 CFR", {
  itr <- c(0:9, 6:9, 6:9, 0:9, 6:9)
  cfr <- numeric(length(itr))
  cfr[c(5, 7)]  <- c(0.85, 0.55)   # trace 1, loc 1: itr 4 then itr 6
  cfr[11]       <- 0.75            # trace 1, loc 2: itr 6
  cfr[15]       <- 0.30            # trace 1, loc 3: itr 6
  cfr[c(23, 25)] <- c(0.90, 0.62)  # trace 2, loc 1
  cfr[29]       <- 0.41            # trace 2, loc 2
  tid <- rep(c(1L, 2L), c(18, 14))
  out <- .minflux_finalize_df(.raw_rows(tid, itr, cfr), use_fnl = FALSE)
  expect_equal(nrow(out), 5)
  expect_equal(out$cfr, c(0.55, 0.75, 0.30, 0.62, 0.41))
})

# CLAIM: the fnl flag (.npy) defines the final rows the same way as itr == max.
test_that("fnl flag gives the same per-localization CFR", {
  itr <- c(0:9, 6:9)
  cfr <- c(0, 0, 0, 0, 0.85, 0, 0.55, 0, 0, 0, 0.75, 0, 0, 0)
  out <- .minflux_finalize_df(.raw_rows(1L, itr, cfr, fnl = itr == 9), use_fnl = TRUE)
  expect_equal(out$cfr, c(0.55, 0.75))
})

# CLAIM: in 2D tracking the final iteration measures CFR itself, so every final
# row keeps its own value (the old per-trace carry overwrote all of them).
test_that("2D tracking: final rows keep their own CFR", {
  itr <- c(0:3, 3, 3, 3)
  cfr <- c(0, 0.31, 0, 1.34, 1.25, 0.20, 0.99)
  out <- .minflux_finalize_df(.raw_rows(1L, itr, cfr), use_fnl = FALSE)
  expect_equal(out$cfr, c(1.34, 1.25, 0.20, 0.99))
})

# CLAIM: a localization with no positive CFR in its own run keeps its own value;
# it never borrows from a neighbouring localization.
test_that("no borrowing across localizations", {
  itr <- c(0:3, 3, 3)
  cfr <- c(0, 0.31, 0, 1.34, 0, 0.99)
  out <- .minflux_finalize_df(.raw_rows(1L, itr, cfr), use_fnl = FALSE)
  expect_equal(out$cfr, c(1.34, 0, 0.99))
})

# CLAIM: tracking sequences without the CFR check (shorter cycle time) store no
# CFR; the column becomes 0 throughout so the CFR filter is a no-op.
test_that("tracking without CFR measurements gives cfr = 0", {
  itr <- c(0:3, 3, 3)
  r <- .raw_rows(1L, itr, rep(NA_real_, length(itr)))
  expect_equal(.minflux_finalize_df(r)$cfr, c(0, 0, 0))
  r$cfr <- NULL
  expect_equal(.minflux_finalize_df(r)$cfr, c(0, 0, 0))
})

# CLAIM: "itr6" takes the value from iteration 6 only, "none" leaves the final rows raw.
test_that("itr6 and none strategies", {
  itr <- c(0:9, 6:9)
  cfr <- c(0, 0, 0, 0, 0.85, 0, 0.55, 0, 0, 0, 0.75, 0, 0, 0)
  r <- .raw_rows(1L, itr, cfr)
  expect_equal(.minflux_finalize_df(r, cfr_strategy = "itr6")$cfr, c(0.55, 0.75))
  expect_equal(.minflux_finalize_df(r, cfr_strategy = "none")$cfr, c(0, 0))
})
