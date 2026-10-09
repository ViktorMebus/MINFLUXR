## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_plot.R  -  Axes, colour and plot points
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Equal-aspect axis limits, palettes and per-point colours, and the point
##   tables the 2D/3D plots draw (localizations, trace medians or molecule
##   medians). Display only: nothing here changes a computed number.
##
## CONTENTS (in file order)
##   get_equal_limits_2d, NUMERIC_PALETTES, get_palette_hex,
##   build_colors_for_pts, compute_spatial_density, build_pts_2d,
##   build_pts_3d, color_mapping_for_pts, palette_colors, plotly_colorscale,
##   get_equal_limits, QUAL_COLORS, qual_palette, .adjust_lum, shade_colors,
##   density_axis_range, .draw_centroid_gap, expand_limits_with_image,
##   heatmap_input_df, compute_kde, sd_hist_clip, .int_input, .name_or
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

# ── get_equal_limits_2d() ────────────────────────────────────────────────────
# GOAL: return x/y axis limits that span the SAME distance, so a nanometre on
# screen is the same length on both axes.
# This is not cosmetic: with unequal limits plotly stretches one axis, and a
# round cluster renders as an ellipse - which would corrupt every visual shape
# judgement the user makes. Padding adds a small margin so points do not sit on
# the frame.
get_equal_limits_2d <- function(x, y, padding = 0.05) {
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  
  # Handle degenerate ranges (all points same)
  if (!is.finite(diff(xr)) || diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  if (!is.finite(diff(yr)) || diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
  
  span <- max(diff(xr), diff(yr))
  pad  <- span * padding
  
  x_mid <- mean(xr); y_mid <- mean(yr)
  x_lim <- c(x_mid - span/2 - pad, x_mid + span/2 + pad)
  y_lim <- c(y_mid - span/2 - pad, y_mid + span/2 + pad)
  
  list(x = x_lim, y = y_lim)
}

# Color palette utilities for heatmaps and numeric mapping
NUMERIC_PALETTES <- c(
  "Viridis", "Plasma", "Inferno", "Magma", "Cividis",
  "Turbo",
  "Spectral", "RdBu", "RdYlBu", "RdYlGn",
  "YlOrRd", "Greys"
)

# Return a vector of n hex colors for the given palette name
get_palette_hex <- function(name, n = 256, reverse = FALSE) {
  name <- if (is.null(name) || !nzchar(name)) "Viridis" else as.character(name)
  
  cols <- switch(
    name,
    # viridisLite family
    "Viridis" = viridisLite::viridis(n),
    "Plasma"  = viridisLite::plasma(n),
    "Inferno" = viridisLite::inferno(n),
    "Magma"   = viridisLite::magma(n),
    "Cividis" = viridisLite::cividis(n),
    
    "Turbo" = tryCatch(
      viridisLite::turbo(n),
      error = function(e1) tryCatch(
        grDevices::hcl.colors(n, "turbo"),
        error = function(e2) viridisLite::viridis(n)
      )
    ),
    
    # Brewer (diverging/sequential)
    "Spectral" = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n),
    "RdBu"     = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(n),
    "RdYlBu"   = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n),
    "RdYlGn"   = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(n),
    "YlOrRd"   = colorRampPalette(RColorBrewer::brewer.pal(9,  "YlOrRd"))(n),
    "Greys"    = colorRampPalette(RColorBrewer::brewer.pal(9,  "Greys"))(n),
    
    # default
    viridisLite::viridis(n)
  )
  
  if (isTRUE(reverse)) cols <- rev(cols)
  cols
}

# GOAL: one colour per point for the current colouring mode.
#   display : display options (SMLM_DISPLAY_DEFAULTS; smlm_display_from_input)
# Always returns a color vector length nrow(df_pts)
build_colors_for_pts <- function(df_pts, display = NULL, default_palette_col = "color") {
  d    <- .smlm_display(display)
  mode <- d$color_mode
  map_to_gradient <- function(v, pal_name = NULL, reverse = FALSE) {
    v <- as.numeric(v)
    if (!length(v) || all(!is.finite(v))) return(rep("#bdbdbd", length(v)))
    lo <- suppressWarnings(min(v[is.finite(v)], na.rm = TRUE))
    hi <- suppressWarnings(max(v[is.finite(v)], na.rm = TRUE))
    # avoid flat (lo==hi) -> single color
    if (!is.finite(lo) || !is.finite(hi)) return(rep("#bdbdbd", length(v)))
    if (lo == hi) {
      # add tiny epsilon so a scale can still render different bins if needed
      eps <- max(1e-9, abs(lo) * 1e-6)
      v <- v + seq_len(length(v)) * eps
      lo <- min(v, na.rm = TRUE); hi <- max(v, na.rm = TRUE)
    }
    pal <- get_palette_hex(pal_name %||% d$numeric_palette, n = 256, reverse = isTRUE(reverse))
    idx <- round(1 + 255 * (v - lo) / (hi - lo))
    idx[!is.finite(idx)] <- 128
    pal[pmax(1, pmin(256, idx))]
  }
  # -----------------------------
  
  if (mode == "default") {
    # Respect existing categorical colors (trace/cluster) you already computed
    col <- df_pts[[default_palette_col]]
    if (is.null(col)) col <- rep(d$unicolor_hex, nrow(df_pts))
    col[!nzchar(col)] <- "#bdbdbd"
    col[is.na(col)]   <- "#bdbdbd"
    return(col)
  }
  
  if (mode == "unicolor") {
    return(rep(d$unicolor_hex, nrow(df_pts)))
  }

  if (mode == "multiplex") {
    col <- df_pts[["mplex_color"]]
    if (is.null(col)) return(rep("#bdbdbd", nrow(df_pts)))
    col[is.na(col) | !nzchar(col)] <- "#bdbdbd"
    return(col)
  }

  # numeric modes: pick the right vector depending on raw vs median
  if (mode == "z") {
    # raw points: 'z' exists; medians: we put 'z_med' in builders
    base <- if ("z" %in% names(df_pts)) df_pts$z else df_pts$z_med
    return(map_to_gradient(base, pal_name = d$numeric_palette, reverse = d$numeric_palette_reverse))
    
  } else if (mode == "tim") {
    # medians: 'tim_med'; raw: fallback to per-point 'tim'
    base <- if ("tim_med" %in% names(df_pts)) df_pts$tim_med else df_pts$tim
    return(map_to_gradient(base, pal_name = d$numeric_palette, reverse = d$numeric_palette_reverse))
    
  }  else if (mode %in% c("density","locdens")) {
    # precomputed in builders: 'dens' (counts per TID/cluster for medians; per-trace for raw)
    base <- if ("dens" %in% names(df_pts)) df_pts$dens else rep(1, nrow(df_pts))
    # if constant, nudge to avoid flat coloring
    if (length(unique(base[is.finite(base)])) <= 1L) {
      eps <- max(1e-9, abs(as.numeric(base[which(is.finite(base))[1]]) %||% 1) * 1e-6)
      base <- as.numeric(base) + seq_len(length(base)) * eps
    }
    return(map_to_gradient(base, pal_name = d$numeric_palette, reverse = d$numeric_palette_reverse))
  }
  
  # fallback
  rep("#bdbdbd", nrow(df_pts))
}

# --- compute_spatial_density ---
compute_spatial_density <- function(x, y, z = NULL, r = 150, k_max = 1000) {
  # Returns integer vector dens, same length as x
  stopifnot(length(x) == length(y))
  n <- length(x)
  if (n == 0) return(integer(0))
  if (n == 1) return(0L)
  
  pts <- if (is.null(z)) cbind(x, y) else cbind(x, y, z)
  
  # frNN returns neighbor lists within eps=r (includes self)
  nn <- dbscan::frNN(pts, eps = r, sort = FALSE)
  
  # density = neighbors within r EXCLUDING self
  dens <- lengths(nn$id) - 1L
  dens[dens < 0] <- 0L
  as.integer(dens)
}

# --- build_pts_2d ---
#   display   : display options (dens_radius_nm, show_medians are read)
#   deblinked : TRUE when de-blinking ran, so medians are per molecule (`cluster`)
build_pts_2d <- function(df, display = NULL, deblinked = FALSE) {
  d <- .smlm_display(display)
  stopifnot(all(c("x","y","tid") %in% names(df)))
  
  # --- 1) compute per-tid medians once ---
  tid_meds <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z_med   = if ("z"   %in% names(df)) stats::median(z,   na.rm = TRUE) else NA_real_,
      tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
      color = dplyr::first(color),
      hover_text = dplyr::first(hover_text),
      mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
      mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
      mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(x), is.finite(y))

  # --- 2) spatial density of trace medians (2D) ---
  r_nm <- d$dens_radius_nm
  tid_meds$dens <- compute_spatial_density(tid_meds$x, tid_meds$y, r = r_nm)
  
  if (isTRUE(d$show_medians)) {
    
    # If clustering AND you are showing cluster medians:
    # density should probably be computed on *cluster medians* (not tid medians),
    # because that's what you're plotting. We'll keep that consistent.
    if (isTRUE(deblinked)) {
      pts <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z_med   = if ("z"   %in% names(df)) stats::median(z,   na.rm = TRUE) else NA_real_,
          tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
          color  = dplyr::first(color),
          hover_text = dplyr::first(hover_text),
          mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
          .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(x), is.finite(y))

      pts$dens <- compute_spatial_density(pts$x, pts$y, r = r_nm)
      pts <- pts %>% dplyr::mutate(group_id = cluster)
      
    } else {
      # Plot the tid medians directly (already has dens)
      pts <- tid_meds %>% dplyr::mutate(group_id = tid)
    }
    
  } else {
    # RAW points: join tid-level spatial density onto every localization
    pts <- df %>%
      dplyr::left_join(tid_meds %>% dplyr::select(tid, dens), by = "tid") %>%
      dplyr::mutate(
        z_med   = if ("z"   %in% names(.)) z   else NA_real_,
        tim_med = if ("tim" %in% names(.)) tim else NA_real_,
        group_id = if (isTRUE(deblinked) && "cluster" %in% names(.))
          ifelse(is.na(cluster) | cluster == 0, -1L, cluster) else tid
      )
  }
  
  # safety columns
  for (nm in c("x","y","z_med","tim_med","dens","color","hover_text","group_id",
               "mplex_color","mplex_channel","mplex_label")) {
    if (!nm %in% names(pts)) pts[[nm]] <- NA
  }
  pts %>% dplyr::filter(is.finite(x), is.finite(y))
}

# --- build_pts_3d ---
#   display   : display options (dens_radius_nm, show_medians are read)
#   deblinked : TRUE when de-blinking ran, so medians are per molecule (`cluster`)
build_pts_3d <- function(df, display = NULL, deblinked = FALSE) {
  d <- .smlm_display(display)
  stopifnot(all(c("x","y","z","tid") %in% names(df)))
  
  tid_meds <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z = stats::median(z, na.rm = TRUE),
      tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
      color = dplyr::first(color),
      hover_text = dplyr::first(hover_text),
      mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
      mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
      mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(x), is.finite(y), is.finite(z))

  r_nm <- d$dens_radius_nm
  tid_meds$dens <- compute_spatial_density(tid_meds$x, tid_meds$y, tid_meds$z, r = r_nm)
  
  if (isTRUE(d$show_medians)) {
    if (isTRUE(deblinked)) {
      pts <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          tim_med = if ("tim" %in% names(df)) stats::median(tim, na.rm = TRUE) else NA_real_,
          color  = dplyr::first(color),
          hover_text = dplyr::first(hover_text),
          mplex_color   = if ("mplex_color"   %in% names(df)) dplyr::first(mplex_color)   else NA_character_,
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          mplex_label   = if ("mplex_label"   %in% names(df)) dplyr::first(mplex_label)   else NA_character_,
          .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z))

      pts$dens <- compute_spatial_density(pts$x, pts$y, pts$z, r = r_nm)
      pts <- pts %>% dplyr::mutate(group_id = cluster)
      
    } else {
      pts <- tid_meds %>% dplyr::mutate(group_id = tid)
    }
    
  } else {
    pts <- df %>%
      dplyr::left_join(tid_meds %>% dplyr::select(tid, dens), by = "tid") %>%
      dplyr::mutate(
        tim_med = if ("tim" %in% names(.)) tim else NA_real_,
        group_id = if (isTRUE(deblinked) && "cluster" %in% names(.))
          ifelse(is.na(cluster) | cluster == 0, -1L, cluster) else tid
      )
  }
  
  for (nm in c("x","y","z","tim_med","dens","color","hover_text","group_id",
               "mplex_color","mplex_channel","mplex_label")) {
    if (!nm %in% names(pts)) pts[[nm]] <- NA
  }
  pts %>% dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
}

# --- color_mapping_for_pts ---
# GOAL: the plotly colour specification for the points (a hex vector, or a
# numeric colour scale with its range) in the current colouring mode.
#   display : display options (SMLM_DISPLAY_DEFAULTS; smlm_display_from_input)
color_mapping_for_pts <- function(pts, display = NULL) {
  d        <- .smlm_display(display)
  mode     <- d$color_mode
  pal_name <- d$numeric_palette
  pal_rev  <- isTRUE(d$numeric_palette_reverse)
  
  # Helper: build a numeric colorscale mapping
  build_num_map <- function(v, title) {
    v <- as.numeric(v)
    # All-NA: still return a colorscale shell so Plotly doesn’t choke
    if (!any(is.finite(v))) {
      return(list(
        showscale  = FALSE,
        kind       = "hex",                # fall back to hex (no scale)
        hex        = rep("#bdbdbd", nrow(pts))
      ))
    }
    v_ok <- v[is.finite(v)]
    cmin <- min(v_ok); cmax <- max(v_ok)
    if (!is.finite(cmin) || !is.finite(cmax) || cmin == cmax) {
      # No spread → solid color (pick mid color from scale to be deterministic)
      mid_cols <- plotly_colorscale(pal_name, reverse = pal_rev)
      mid_hex  <- mid_cols[[ceiling(length(mid_cols)/2)]][[2]]
      return(list(
        showscale  = FALSE,
        kind       = "hex",
        hex        = rep(mid_hex, nrow(pts))
      ))
    }
    list(
      showscale  = TRUE,
      kind       = "scale",
      colorvar   = v,
      cmin       = cmin,
      cmax       = cmax,
      colorscale = plotly_colorscale(pal_name, reverse = pal_rev),
      title      = title
    )
  }
  
  if (mode == "default") {
    col <- pts$color
    if (is.null(col)) col <- rep(d$unicolor_hex, nrow(pts))
    col[!nzchar(col)] <- "#bdbdbd"; col[is.na(col)] <- "#bdbdbd"
    return(list(showscale = FALSE, kind = "hex", hex = col))
  }
  
  if (mode == "unicolor") {
    return(list(
      showscale = FALSE, kind = "hex",
      hex = rep(d$unicolor_hex, nrow(pts))
    ))
  }

  if (mode == "multiplex") {
    col <- if ("mplex_color" %in% names(pts)) pts$mplex_color else pts$color
    if (is.null(col)) col <- rep("#bdbdbd", nrow(pts))
    col[is.na(col) | !nzchar(col)] <- "#bdbdbd"
    return(list(showscale = FALSE, kind = "hex", hex = col))
  }

  if (mode == "z") {
    zvec <- if ("z" %in% names(pts)) pts$z else if ("z_med" %in% names(pts)) pts$z_med else NA_real_
    return(build_num_map(zvec, "Depth (z)"))
  }
  if (mode == "tim")     return(build_num_map(if ("tim_med" %in% names(pts)) pts$tim_med else pts$tim, "Time"))
  if (mode %in% c("density","locdens")) {
    return(build_num_map(
      if ("dens" %in% names(pts)) pts$dens else rep(1, nrow(pts)),
      "Local trace-median density"
    ))
  }
  
  # Fallback
  list(showscale = FALSE, kind = "hex", hex = rep("#bdbdbd", nrow(pts)))
}

# palette helper (uses RColorBrewer that you already load)
palette_colors <- function(name, n = 256) {
  get_palette_hex(name, n = n, reverse = FALSE)
}

plotly_colorscale <- function(name = "Viridis", reverse = FALSE, n = 11) {
  cols <- get_palette_hex(name, n = n, reverse = reverse)
  pos  <- seq(0, 1, length.out = length(cols))
  Map(function(p, c) list(p, c), pos, cols)
}

#3D Helpers


get_equal_limits <- function(x, y, z = NULL, padding = 0.05) {
  xr <- range(x, na.rm = TRUE); if (!is.finite(diff(xr)) || diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  yr <- range(y, na.rm = TRUE); if (!is.finite(diff(yr)) || diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
  
  if (is.null(z)) {
    span <- max(diff(xr), diff(yr))
    pad  <- span * padding
    xm <- mean(xr); ym <- mean(yr)
    return(list(
      x = c(xm - span/2 - pad, xm + span/2 + pad),
      y = c(ym - span/2 - pad, ym + span/2 + pad)
    ))
  } else {
    zr <- range(z, na.rm = TRUE); if (!is.finite(diff(zr)) || diff(zr) == 0) zr <- zr + c(-0.5, 0.5)
    span <- max(diff(xr), diff(yr), diff(zr))
    pad  <- span * padding
    xm <- mean(xr); ym <- mean(yr); zm <- mean(zr)
    return(list(
      x = c(xm - span/2 - pad, xm + span/2 + pad),
      y = c(ym - span/2 - pad, ym + span/2 + pad),
      z = c(zm - span/2 - pad, zm + span/2 + pad)
    ))
  }
}

## ── Series colours for the batch/spatial panels ──────────────────────────────
## GOAL: two palettes with different jobs.
##
## `qual_palette(n)` - CATEGORICAL colours for one-series-per-file plots. Files are
## nominal, not ordered, so a sequential ramp is the wrong encoding: viridis on 10
## files reads as a blue-green gradient in which neighbouring files are nearly
## indistinguishable. This is a hand-ordered set of maximally separated hues
## (Dark2 then Set1 then darker Paired entries), recycled with a darkening pass if
## more series than colours. DETERMINISTIC, not random - a random draw would
## re-colour every series on every redraw and can collide by chance.
##
## `shade_colors(base, n)` - ORDERED shades of one base colour, for individual
## files inside a group. Here the encoding is right: the group identity is the hue
## and the file is a shade of it, so a red group and a blue group stay readable at a
## glance while individual files remain separable.
QUAL_COLORS <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
  "#a6761d", "#666666", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
  "#ff7f00", "#a65628", "#f781bf", "#1f78b4", "#33a02c", "#6a3d9a",
  "#b15928", "#b2182b"
)

qual_palette <- function(n) {
  n <- max(1L, as.integer(n))
  base <- QUAL_COLORS
  if (n <= length(base)) return(base[seq_len(n)])
  # More series than distinct hues: recycle, darkening each successive lap so a
  # repeat is still visually distinct from its first appearance.
  reps <- ceiling(n / length(base))
  out  <- unlist(lapply(seq_len(reps) - 1L, function(k)
    if (k == 0L) base else .adjust_lum(base, -0.18 * k)))
  out[seq_len(n)]
}

# Shift luminance of hex colours: amt > 0 lightens toward white, < 0 darkens.
.adjust_lum <- function(hex, amt) {
  m <- grDevices::col2rgb(hex) / 255
  m <- if (amt >= 0) m + (1 - m) * amt else m * (1 + amt)
  m <- pmin(pmax(m, 0), 1)
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

shade_colors <- function(base, n) {
  n <- max(1L, as.integer(n))
  if (n == 1L) return(base)
  # Light tint -> base -> dark shade. Kept away from the extremes so no file ends
  # up near-white (invisible) or near-black (indistinguishable from another group).
  grDevices::colorRampPalette(c(.adjust_lum(base, 0.55), base,
                                .adjust_lum(base, -0.30)))(n)
}

## GOAL: the x-axis range for the density-distribution plot, and the clipping span
## that goes with it.
##
## WHY NOT AUTORANGE: a violin's kernel tails run well past the data, and on a log
## axis plotly then pads out to whole decade boundaries - which leaves the panel
## long and mostly empty. So the span is taken from a robust quantile of the local
## densities (`probs`), widened to include the per-file global densities so those
## markers can never fall off frame, and set explicitly.
##
## Values outside `lim` are meant to be CLIPPED to the edge by the caller rather
## than dropped, so each violin ends at the frame instead of trailing off it.
##
## Returns list(lim, range, clip) where `range` is what plotly's xaxis wants
## (log10 of the padded limits when use_log, the padded limits otherwise) and
## `clip` is a function mapping values into `lim`.
density_axis_range <- function(local_vals, global_vals, use_log = TRUE,
                               probs = c(0.005, 0.995)) {
  lv <- local_vals[is.finite(local_vals) & local_vals > 0]
  gv <- global_vals[is.finite(global_vals)]
  # Both empty is reachable (a panel rendering before any measurement); range() on
  # nothing warns and returns c(Inf, -Inf), so short-circuit to a harmless frame.
  if (!length(lv) && !length(gv))
    return(list(lim = c(1, 10), range = if (use_log) c(0, 1) else c(0, 10),
                clip = function(v) v))
  lim <- if (length(lv) >= 10L)
    as.numeric(stats::quantile(lv, probs, na.rm = TRUE)) else range(c(lv, gv))
  lim <- range(c(lim, gv))
  if (!all(is.finite(lim)) || lim[2] <= lim[1])
    lim <- c(max(1e-9, min(lim, na.rm = TRUE)), max(1, max(lim, na.rm = TRUE)) * 10)
  if (use_log) {
    pos <- c(lv, gv[gv > 0])
    if (length(pos)) lim[1] <- max(lim[1], min(pos) / 1.6)
    lim[1] <- max(lim[1], 1e-9)
    rng <- log10(pmax(c(lim[1] / 1.25, lim[2] * 1.25), 1e-12))
  } else {
    sp  <- diff(lim)
    rng <- c(max(0, lim[1] - 0.04 * sp), lim[2] + 0.04 * sp)
  }
  list(lim = lim, range = rng,
       clip = function(v) pmin(pmax(v, lim[1]), lim[2]))
}


# =============================================================================
# PANEL HELPERS (moved from server_smlm.R; plain functions of their arguments)
# =============================================================================

# Drawn by BOTH the single-colour and two-colour panels. Factored so the two
# can never drift apart visually or in what the lines mean.
.draw_centroid_gap <- function(r, cur_thr) {
  wi <- r$within[is.finite(r$within) & r$within > 0]
  bw <- r$between[is.finite(r$between) & r$between > 0]
  validate(need(length(wi) >= 2, "Not enough within-trace distances to plot."))

  COL_W <- "#1f78b4"; COL_B <- "#e08214"; COL_S <- "#1a9850"; COL_E <- "#d7191c"
  # LOG x-axis. The re-blink hump sits at a few nm and the distinct-molecule
  # hump at hundreds; on a linear axis clipped to a quantile of the pooled
  # distances the second hump fell off the plot entirely, which defeats the
  # point - the threshold is supposed to be placed in the GAP between them.
  xhi <- max(stats::quantile(bw, 0.99, names = FALSE), stats::quantile(wi, 0.99, names = FALSE),
             3 * r$suggested, na.rm = TRUE)
  xlo <- max(0.5, min(stats::quantile(c(wi, bw), 0.005, names = FALSE), r$sigma / 4))
  scal <- function(v) {
    v <- v[v > 0]
    d <- stats::density(log10(v), n = 512L, from = log10(xlo), to = log10(xhi))
    d$x <- 10^d$x; d$y <- d$y / max(d$y); d
  }
  # Legend goes in the BOTTOM margin. Note xpd is NOT set in par(): doing so
  # let abline() draw its vertical markers out through the margins too, which
  # is what made the lines run past the plot box.
  op <- graphics::par(mar = c(8.2, 4.5, 3.0, 1.2))
  on.exit(graphics::par(op), add = TRUE)
  dw <- scal(wi)
  graphics::plot(dw$x, dw$y, type = "l", col = COL_W, lwd = 2.5, ylim = c(0, 1.05),
                 xlim = c(xlo, xhi), log = "x",
                 xlab = "distance (nm, log scale)", ylab = "density (scaled to peak)",
                 main = "Within-trace scatter vs between-trace centroid spacing")
  if (length(bw) >= 2) { db <- scal(bw); graphics::lines(db$x, db$y, col = COL_B, lwd = 2.5) }
  graphics::abline(v = r$sigma,     col = COL_W, lwd = 2, lty = 3)
  graphics::abline(v = r$suggested, col = COL_S, lwd = 2, lty = 4)
  if (length(cur_thr) == 1 && is.finite(cur_thr))
    graphics::abline(v = cur_thr, col = COL_E, lwd = 2.5, lty = 2)
  # Placed in device coords so the log x-axis cannot distort it; xpd = NA is
  # scoped to this call alone.
  graphics::legend(x = graphics::grconvertX(0.02, "ndc", "user"),
                   y = graphics::grconvertY(0.14, "ndc", "user"),
    bty = "n", cex = 0.85, ncol = 2, xpd = NA,
    legend = c("within trace (loc to its centroid)", "between traces (centroid NN)",
               "sigma", "suggested", "current threshold"),
    col = c(COL_W, COL_B, COL_W, COL_S, COL_E),
    lty = c(1, 1, 3, 4, 2), lwd = c(2.5, 2.5, 2, 2, 2.5))
  invisible(NULL)
}

# helper: expand axes to include the image bbox; safe if overlay is NULL
expand_limits_with_image <- function(lim, overlay, invert_y = TRUE) {
  if (is.null(overlay)) return(lim)
  x0 <- overlay$x;  w <- overlay$sizex
  y  <- overlay$y;  h <- overlay$sizey
  xr <- c(x0, x0 + w)
  yr <- if (invert_y) c(y - h, y) else c(y, y + h)
  list(x = range(c(lim$x, xr), na.rm = TRUE),
       y = range(c(lim$y, yr), na.rm = TRUE))
}

# Use trace medians (or de-blinked molecule medians when de-blinking is on),
# unless the user explicitly requests raw points.
heatmap_input_df <- function(df, use_raw = FALSE, clustering = FALSE) {
  if (isTRUE(use_raw)) {
    return(df %>% dplyr::select(x, y))
  }
  if (isTRUE(clustering) && "cluster" %in% names(df) && any(df$cluster > 0, na.rm = TRUE)) {
    return(
      df %>%
        dplyr::filter(is.finite(cluster) & cluster > 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          .groups = "drop"
        )
    )
  }
  if ("tid" %in% names(df)) {
    return(
      df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          .groups = "drop"
        )
    )
  }
  df %>% dplyr::select(x, y)
}

compute_kde <- function(df, lim, n = 200, bw = NA_real_, bw_factor = 1) {
  if (!all(c("x","y") %in% names(df)) || nrow(df) < 2) return(NULL)
  n <- max(25, min(1000, as.integer(n)))
  if (any(!is.finite(c(lim$x, lim$y)))) return(NULL)
  
  if (is.finite(bw) && bw > 0) {
    # ABSOLUTE bandwidth (same for x and y)
    h <- c(bw, bw)
  } else {
    # Scott rule (isotropic), then scale by bw_factor
    vx  <- stats::var(df$x, na.rm = TRUE)
    vy  <- stats::var(df$y, na.rm = TRUE)
    sdx <- sqrt(max(vx, 0)); sdy <- sqrt(max(vy, 0))
    sx  <- min(sdx, stats::IQR(df$x, na.rm = TRUE) / 1.34)
    sy  <- min(sdy, stats::IQR(df$y, na.rm = TRUE) / 1.34)
    s   <- sqrt((sx^2 + sy^2) / 2)
    npts <- nrow(df)
    hiso <- 1.06 * s * (npts ^ (-1/6))
    if (!is.finite(hiso) || hiso <= 0) {
      rngx <- diff(range(df$x, na.rm = TRUE))
      rngy <- diff(range(df$y, na.rm = TRUE))
      hiso <- max((rngx + rngy)/2 * 0.02, 1e-6)
    }
    hiso <- hiso * if (is.finite(bw_factor) && bw_factor > 0) bw_factor else 1
    h <- c(hiso, hiso)
  }
  
  MASS::kde2d(
    df$x, df$y, n = n,
    lims = c(lim$x[1], lim$x[2], lim$y[1], lim$y[2]),
    h = h
  )
}

# Values outside a MANUAL range must be dropped, not drawn. sd_counts() clamps
# anything beyond the outermost break into the edge bin, which is right for an
# auto range (nothing is outside it) but would put every broad trace into one
# tall bar at the right edge of a user-set range - a spike that is pure
# artifact and would be read as a real population in a figure.
sd_hist_clip <- function(v, sh) {
  if (is.null(sh) || !isTRUE(sh$manual)) return(v)
  v <- as.numeric(v)
  v[is.finite(v) & v >= sh$xrange[1] & v <= sh$xrange[2]]
}

## GOAL: read a numericInput as a whole number inside a range, and survive the
## user emptying the box.
##
## WHY `%||%` IS NOT ENOUGH. An emptied numericInput does not send NULL, it sends
## **NA**, and `%||%` only substitutes for NULL. NA then survives as.integer(),
## and min()/max() propagate it rather than clamping it, so
## `max(1L, min(6L, NA))` is NA and the loop below died on `seq_len(NA)` the
## moment the group-count box was cleared. Anything that SIZES something from a
## user-typed number needs this, not `%||%`.
.int_input <- function(x, default, lo = -Inf, hi = Inf) {
  v <- suppressWarnings(as.integer(x))
  if (!length(v) || is.na(v[1])) v <- as.integer(default)
  as.integer(max(lo, min(hi, v[1])))
}

## Look a name up in a named vector WITHOUT the `[[` trap. `v[["absent"]]` raises
## "subscript out of bounds" rather than returning NULL, so `%||%` never sees it
## and cannot supply the default - which is how an ungrouped file took the whole
## enhancement curve down. Single-bracket indexing yields NA for a missing name,
## which this turns into the default.
.name_or <- function(v, key, default) {
  if (is.null(v) || !length(v) || is.null(key) || !length(key)) return(default)
  out <- unname(v[as.character(key)[1]])
  if (length(out) != 1L || is.na(out)) default else out
}

