# The shared cleaning + filter pipeline (R/est_clean.R) that the live Tracking
# view, Plot Data's rechunk and Compare Cells all run.
source(file.path(.repo_root, "R/est_clean.R"))

# One track per row group; dt in ms, one long dark gap in track "a".
.trk <- function() {
  a <- data.frame(tid = "a", tim = c(0, 0.2, 0.4, 0.6, 5, 5.2, 5.4, 5.6, 5.8),
                  efo = c(10, 100, 100, 100, 100, 100, 100, 100, 10), x = 0, y = 0)
  b <- data.frame(tid = "b", tim = seq(0, 3, by = 0.2), efo = 100, x = 0, y = 0)   # continuous
  c <- data.frame(tid = "c", tim = c(0, 1, 2, 3, 4), efo = 100, x = 0, y = 0)       # every step 5 cycles
  rbind(a, b, c)
}

test_that("trim_track_ends peels dim ends but never below the floor", {
  d <- .trk()
  t <- trim_track_ends(d[d$tid == "a", ], t_c = 0.2, use_dark = FALSE, use_efo = TRUE, efo_frac = 0.5)
  expect_equal(range(t$tim), c(0.2, 5.6))                          # both dim ends gone
  f <- trim_track_ends(d[d$tid == "a", ], t_c = 0.2, use_dark = FALSE, use_efo = TRUE, floor_len = 9L)
  expect_equal(nrow(f), 9L)                                        # floor stops the trim
})

test_that("gate_continuity drops dark-cycle-heavy tracks only", {
  g <- gate_continuity(.trk(), t_c = 0.2, min_score = 0.6, half_life = 3)
  expect_setequal(unique(g$tid), c("a", "b"))                      # c: every step 4 dark cycles
})

test_that("range, length and duration filters, with their no-op cases", {
  d <- .trk()
  expect_equal(nrow(filter_ranges(d, list(efo = c(50, 200)))), nrow(d) - 2L)
  expect_equal(nrow(filter_ranges(d, list(efo = c(1e6, 2e6)))), nrow(d))   # no overlap: stale, passes
  expect_equal(nrow(filter_ranges(d, list(efo = NULL))), nrow(d))
  expect_setequal(unique(filter_tid_length(d, c(6, 100))$tid), c("a", "b"))
  expect_setequal(unique(filter_duration(d, c(0, 4.5))$tid), c("b", "c"))   # a lasts 5.8 ms
  expect_equal(nrow(filter_duration(d, c(0, 100))), nrow(d))                # covers all: no-op
})

test_that("filter_pipeline equals running the stages in order", {
  d <- .trk()
  p <- list(t_c = 0.2, sw = 0,
            split = list(enable = TRUE, adaptive = FALSE, gap_cycles = 7, min_len = 2L),
            trim  = list(enable = TRUE, dark = TRUE, dark_cycles = 3, efo = TRUE, efo_frac = 0.5),
            gate  = list(enable = TRUE, min_score = 0.6, half_life = 3),
            ranges = list(efo = c(50, 200)), tid_len = c(3, 100), dur = NULL)
  by_hand <- split_gaps(d, 7 * 0.2, 2L)
  by_hand <- trim_track_ends(by_hand, 0.2, TRUE, 3, TRUE, 0.5, floor_len = 2L)
  by_hand <- gate_continuity(by_hand, 0.2, 0.6, 3)
  by_hand <- filter_ranges(by_hand, list(efo = c(50, 200)))
  by_hand <- filter_tid_length(by_hand, c(3, 100))
  expect_identical(filter_pipeline(d, p), by_hand)
  expect_true(all(grepl("^a_s[12]$|^b$", unique(by_hand$tid))))  # a split at its 4.4 ms gap
})
