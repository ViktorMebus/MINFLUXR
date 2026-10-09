# Display helpers (step 7 of the settings migration): the plot helpers take
# display options (SMLM_DISPLAY_DEFAULTS, R/smlm_settings.R) and plain arguments.

.plot_df <- function() {
  set.seed(2)
  df <- sim_smlm_traces(grid_molecules(4, 200), blink_mean = 1, sigma = 5)
  df$z <- stats::rnorm(nrow(df), 0, 30); df$tim <- seq_len(nrow(df)) / 10
  df$cluster <- df$mol; df$color <- "#123456"; df$hover_text <- paste("tid", df$tid)
  df
}

test_that("display options are completed from the defaults, whatever the input", {
  expect_identical(.smlm_display(NULL), SMLM_DISPLAY_DEFAULTS)
  d <- .smlm_display(list(color_mode = "z", not_a_display_option = 1))
  expect_identical(d$color_mode, "z")
  expect_identical(d$numeric_palette, SMLM_DISPLAY_DEFAULTS$numeric_palette)
  expect_false("not_a_display_option" %in% names(d))
  rv <- shiny::reactiveValues(color_mode = "tim", dens_radius_nm = 250)
  d2 <- shiny::isolate(.smlm_display(rv))
  expect_identical(d2$color_mode, "tim"); expect_equal(d2$dens_radius_nm, 250)
  expect_identical(smlm_channel_colors(list(color_group2 = "#00ff00")),
                   c(SMLM_DISPLAY_DEFAULTS$color_group1, "#00ff00"))
})

test_that("channel colours follow the channel order and keep multiplex colours", {
  df <- data.frame(channel = c(1L, 2L, 2L, 1L))
  out <- attach_channel_colors(df, c("#aa0000", "#0000bb"))
  expect_identical(out$channel_color, c("#aa0000", "#0000bb", "#0000bb", "#aa0000"))
  expect_identical(attach_channel_colors(df)$channel_color[1], SMLM_DISPLAY_DEFAULTS$color_group1)
  mp <- data.frame(channel = c(1L, 2L), mplex_color = c("#111111", "#222222"))
  expect_identical(attach_channel_colors(mp, c("#aa0000", "#0000bb"))$channel_color,
                   c("#111111", "#222222"))                 # data colour wins
})

# CLAIM: show_medians draws one point per trace, or per molecule when de-blinking
# ran; otherwise every localization, each carrying its trace's local density.
test_that("build_pts draws localizations, trace medians or molecule medians", {
  df <- .plot_df()
  raw <- build_pts_2d(df, list(show_medians = FALSE), deblinked = TRUE)
  expect_equal(nrow(raw), nrow(df))
  med <- build_pts_2d(df, list(show_medians = TRUE), deblinked = FALSE)
  expect_equal(nrow(med), length(unique(df$tid)))
  mol <- build_pts_2d(df, list(show_medians = TRUE), deblinked = TRUE)
  expect_equal(nrow(mol), length(unique(df$cluster)))
  expect_equal(nrow(build_pts_3d(df, list(show_medians = TRUE), deblinked = TRUE)),
               length(unique(df$cluster)))
  # the density radius is the display option, not a hidden constant
  small <- build_pts_2d(df, list(show_medians = TRUE, dens_radius_nm = 10))
  large <- build_pts_2d(df, list(show_medians = TRUE, dens_radius_nm = 1000))
  expect_lt(sum(small$dens), sum(large$dens))
})

test_that("point colours and the colour map follow the display options", {
  pts <- build_pts_2d(.plot_df(), list(show_medians = TRUE), deblinked = FALSE)
  expect_true(all(build_colors_for_pts(pts, list(color_mode = "default")) == "#123456"))
  expect_true(all(build_colors_for_pts(pts, list(color_mode = "unicolor",
                                                 unicolor_hex = "#abcdef")) == "#abcdef"))
  z  <- build_colors_for_pts(pts, list(color_mode = "z"))
  zr <- build_colors_for_pts(pts, list(color_mode = "z", numeric_palette_reverse = TRUE))
  lo <- which.min(pts$z_med); hi <- which.max(pts$z_med)
  expect_identical(z[lo], zr[hi])                         # reversed scale swaps the ends
  cm <- color_mapping_for_pts(pts, list(color_mode = "z"))
  expect_identical(cm$kind, "scale")
  expect_equal(c(cm$cmin, cm$cmax), range(pts$z_med))
  expect_identical(color_mapping_for_pts(pts, list(color_mode = "unicolor",
                                                   unicolor_hex = "#abcdef"))$hex[1], "#abcdef")
})

test_that("a plain list and the app's input give the same plot points", {
  df <- .plot_df()
  opts <- list(color_mode = "locdens", show_medians = TRUE, dens_radius_nm = 300)
  rv <- do.call(shiny::reactiveValues, opts)
  a <- build_pts_2d(df, opts, deblinked = TRUE)
  b <- shiny::isolate(build_pts_2d(df, rv, deblinked = TRUE))
  expect_identical(a, b)
  expect_identical(build_colors_for_pts(a, opts), shiny::isolate(build_colors_for_pts(a, rv)))
})

# GUARD: no helper reads Shiny's input directly any more, and the display widgets
# take their initial values from SMLM_DISPLAY_DEFAULTS.
test_that("helpers do not read input$ and display widgets use the shared defaults", {
  h <- paste(readLines(file.path(.repo_root, "helpers_smlm.R"), warn = FALSE, encoding = "UTF-8"),
             collapse = "\n")
  expect_false(grepl("input$", gsub("#[^\n]*", "", h), fixed = TRUE))
  ui <- paste(readLines(file.path(.repo_root, "ui_smlm.R"), warn = FALSE, encoding = "UTF-8"),
              collapse = "\n")
  for (id in names(SMLM_DISPLAY_DEFAULTS)) {
    # up to 14 lines: color_mode lists its choices before `selected =`
    m <- regexpr(sprintf('"%s"[^\\n]*?(\\n[^\\n]*?){0,14}?\\b(value|selected)\\s*=\\s*[^,)]+', id),
                 ui, perl = TRUE)
    expect_true(m > 0, info = id)
    expect_match(regmatches(ui, m), "SMLM_DISPLAY_DEFAULTS\\$", info = id)
  }
})
