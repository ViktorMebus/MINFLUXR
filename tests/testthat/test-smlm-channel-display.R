# Channel-aware display: the channel count shared by the points and the Median
# Cluster wraps (n_frame_channels), and the DCR histogram stacked by channel
# (channel_hist_layers).

# CLAIM: the count covers every row, so a channel with no clusters still counts.
test_that("n_frame_channels counts channels over all rows", {
  expect_equal(n_frame_channels(data.frame(x = 1:3)), 1L)
  expect_equal(n_frame_channels(NULL), 1L)
  expect_equal(n_frame_channels(data.frame(channel = c(1L, 1L, NA))), 1L)
  # Only channel 2 holds Median Clusters; the frame still has two channels.
  df <- data.frame(channel = c(1L, 1L, 2L, 2L), cluster2 = c(0L, 0L, 5L, 5L))
  expect_equal(n_frame_channels(df), 2L)
  expect_equal(length(unique(wrap_groups(df, "cluster2")$channel)), 1L)
})

.hist_df <- function() {
  df <- data.frame(tid = 1:6, channel = c(1L, 1L, 1L, 2L, 2L, 2L),
                   dcr = c(0.1, 0.2, 0.6, 0.5, 0.8, 0.9))
  attr(df, "channel_table") <- data.frame(channel = 1:2, label = c("A", "B"),
                                          color = c("#112233", "#445566"))
  df
}

# CLAIM: the channels stack up to exactly the filtered histogram, each one
# starting where the previous ended, in the colors of the channel table.
test_that("channel layers stack to the filtered histogram in the channel colors", {
  df <- .hist_df()
  L <- channel_hist_layers(df, "dcr", df$dcr, n_bins = 4L)
  expect_length(L$layers, 2L)
  expect_identical(vapply(L$layers, `[[`, "", "color"), c("#112233", "#445566"))
  expect_identical(vapply(L$layers, `[[`, "", "name"), c("Channel A", "Channel B"))
  expect_equal(L$layers[[2]]$base, L$layers[[1]]$count)
  total <- L$layers[[2]]$base + L$layers[[2]]$count
  expect_equal(total, graphics::hist(df$dcr, breaks = L$breaks, plot = FALSE,
                                     right = FALSE, include.lowest = TRUE)$counts)
  expect_equal(sum(total), nrow(df))
  expect_null(L$raw)
})

# CLAIM: raw data uses the same breaks, and may extend the axis beyond the filtered data.
test_that("raw values are binned on the shared breaks", {
  df <- .hist_df()
  L <- channel_hist_layers(df, "dcr", df$dcr, vals_raw = c(df$dcr, 1.3), n_bins = 6L)
  expect_equal(max(L$breaks), 1.3)
  expect_equal(sum(L$raw), nrow(df) + 1L)
  expect_equal(sum(L$layers[[1]]$count) + sum(L$layers[[2]]$count), nrow(df))
})

# CLAIM: traces the split dropped appear as their own layer, so the top of the
# stack stays the filtered histogram.
test_that("localizations without a channel get a Not assigned layer", {
  df <- .hist_df()
  L <- channel_hist_layers(df, "dcr", c(df$dcr, 0.4, 0.4), n_bins = 5L)
  expect_length(L$layers, 3L)
  expect_identical(L$layers[[3]]$name, "Not assigned")
  expect_equal(sum(L$layers[[3]]$count), 2)
})

test_that("channel_hist_layers declines when there is nothing to split", {
  df <- .hist_df()
  expect_null(channel_hist_layers(NULL, "dcr", df$dcr))
  expect_null(channel_hist_layers(df, "missing", df$dcr))
  expect_null(channel_hist_layers(df, "dcr", numeric(0)))
  no_tab <- df; attr(no_tab, "channel_table") <- NULL
  expect_null(channel_hist_layers(no_tab, "dcr", df$dcr))
})
