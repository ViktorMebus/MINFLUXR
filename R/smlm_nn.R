## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_nn.R  -  Nearest-neighbour distances
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   k-NN distances within and between channels and per Median Cluster, for
##   the NN and cross-channel NN panels.
##
## CONTENTS (in file order)
##   nn_index_pairs, compute_batch_nn, safe_nn2, compute_cross_nn,
##   .nn_pairs_scoped, build_nn_segments, perp_offset_xy, member_nn_stats
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

# Nearest-neighbour index pairs from a coordinate matrix. Returns a 3-col matrix
# (i, ni, dist); indices are into the original coords rows.
#   channel = NULL          -> global k-NN (self excluded).
#   channel + mode="within" -> per-channel k-NN (neighbours only from same channel).
#   channel + mode="between"-> cross-channel k-NN: each point's nearest neighbours
#                              are drawn ONLY from OTHER channels (never its own).
# These are true searches (not post-hoc filters of a global search).
nn_index_pairs <- function(coords, k, channel = NULL, mode = "within") {
  n <- nrow(coords)
  if (is.null(n) || n < 2L) return(matrix(numeric(0), ncol = 3))
  kk <- as.integer(k)

  # within one index set (self excluded)
  build <- function(idx) {
    if (length(idx) < 2L) return(NULL)
    sub   <- coords[idx, , drop = FALSE]
    k_eff <- min(kk + 1L, length(idx))
    if (k_eff < 2L) return(NULL)
    nn <- RANN::nn2(sub, k = k_eff)
    out <- vector("list", length(idx))
    for (a in seq_along(idx)) {
      pr <- lapply(2:k_eff, function(b) {
        ni <- nn$nn.idx[a, b]
        if (is.na(ni) || ni < 1L) return(NULL)
        c(idx[a], idx[ni], nn$nn.dists[a, b])
      })
      out[[a]] <- do.call(rbind, pr)
    }
    do.call(rbind, out)
  }

  if (is.null(channel)) {
    m <- build(seq_len(n))
    return(if (is.null(m)) matrix(numeric(0), ncol = 3) else m)
  }

  if (identical(mode, "between")) {
    chs <- sort(unique(channel[!is.na(channel)]))
    res <- list()
    for (chv in chs) {
      qi <- which(!is.na(channel) & channel == chv)   # query points (this channel)
      di <- which(!is.na(channel) & channel != chv)   # data points (other channels)
      if (!length(qi) || !length(di)) next
      k_eff <- min(kk, length(di))
      if (k_eff < 1L) next
      nn <- RANN::nn2(coords[di, , drop = FALSE],
                      query = coords[qi, , drop = FALSE], k = k_eff)
      for (a in seq_along(qi)) {
        for (b in seq_len(k_eff)) {
          ni <- nn$nn.idx[a, b]
          if (is.na(ni) || ni < 1L) next
          res[[length(res) + 1L]] <- c(qi[a], di[ni], nn$nn.dists[a, b])
        }
      }
    }
    m <- do.call(rbind, res)
    return(if (is.null(m)) matrix(numeric(0), ncol = 3) else m)
  }

  # mode == "within"
  m <- do.call(rbind, lapply(split(seq_len(n), channel), build))
  if (is.null(m)) matrix(numeric(0), ncol = 3) else m
}

# Compute k-NN distances between unit centroids.
# unit:         "tid"           - one centroid per TID (global)
#               "pcid"          - one centroid per de-blinked molecule (global)
#               "cluster2_tid"  - TID centroids within each cluster2
#               "cluster2_pcid" - de-blinked molecule centroids within each cluster2
# group_filter: NA (all), 1L, or 2L to restrict to an analysis_group
# Returns a data.frame with column `distance`.
compute_batch_nn <- function(df, unit = "tid", k = 1L, group_filter = NA_integer_,
                             cross_groups = NULL) {
  if (!all(c("x", "y", "tid") %in% names(df))) return(data.frame())
  if (!"z" %in% names(df)) df$z <- 0

  # Between-group (cross-channel) NN: query from cross_groups[1] to cross_groups[2]
  if (!is.null(cross_groups) && length(cross_groups) == 2L &&
      "analysis_group" %in% names(df)) {
    get_pts_grp <- function(sub) {
      if (unit == "pcid" && "cluster" %in% names(sub) &&
          any(!is.na(sub$cluster) & sub$cluster > 0L)) {
        pts <- sub %>%
          dplyr::filter(!is.na(cluster), cluster > 0L) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                           y = stats::median(y, na.rm = TRUE),
                           z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
          dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
      } else {
        pts <- sub %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                           y = stats::median(y, na.rm = TRUE),
                           z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
          dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
      }
      pts
    }
    df_a <- df[!is.na(df$analysis_group) & df$analysis_group == cross_groups[1L], ]
    df_b <- df[!is.na(df$analysis_group) & df$analysis_group == cross_groups[2L], ]
    pts_a <- get_pts_grp(df_a); pts_b <- get_pts_grp(df_b)
    if (nrow(pts_a) < 1L || nrow(pts_b) < 1L) return(data.frame())
    k_eff <- min(k, nrow(pts_b))
    nn <- RANN::nn2(as.matrix(pts_b[, c("x","y","z")]),
                    as.matrix(pts_a[, c("x","y","z")]), k = k_eff)
    d <- as.vector(nn$nn.dists[, k_eff])
    return(data.frame(distance = d[is.finite(d) & d > 1e-9]))
  }

  if (!is.na(group_filter) && "analysis_group" %in% names(df))
    df <- df[!is.na(df$analysis_group) & df$analysis_group == as.integer(group_filter), , drop = FALSE]

  run_nn <- function(coords, k_eff) {
    k_eff <- min(k_eff, nrow(coords) - 1L)
    if (k_eff < 1L) return(numeric(0))
    nn <- tryCatch(RANN::nn2(coords, coords, k = k_eff + 1L), error = function(e) NULL)
    if (is.null(nn)) return(numeric(0))
    d <- as.vector(nn$nn.dists[, seq_len(k_eff) + 1L, drop = FALSE])
    d[is.finite(d) & d > 1e-9]
  }

  if (unit == "tid") {
    pts <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(pts) < 2L) return(data.frame())
    return(data.frame(distance = run_nn(as.matrix(pts[, c("x","y","z")]), k)))

  } else if (unit == "pcid") {
    if (!"cluster" %in% names(df) || all(is.na(df$cluster) | df$cluster == 0L))
      return(data.frame())
    pts <- df %>%
      dplyr::filter(!is.na(cluster), cluster > 0L) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(pts) < 2L) return(data.frame())
    return(data.frame(distance = run_nn(as.matrix(pts[, c("x","y","z")]), k)))

  } else if (unit == "cluster2_tid") {
    if (!"cluster2" %in% names(df) || all(is.na(df$cluster2) | df$cluster2 == 0L))
      return(data.frame())
    tid_pts <- df %>%
      dplyr::filter(!is.na(cluster2), cluster2 > 0L) %>%
      dplyr::group_by(cluster2, tid) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(tid_pts) < 2L) return(data.frame())
    rows <- lapply(split(tid_pts, tid_pts$cluster2), function(g) {
      d <- run_nn(as.matrix(g[, c("x","y","z")]), k)
      if (!length(d)) return(NULL)
      data.frame(distance = d)
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    return(if (is.null(out)) data.frame() else out)

  } else if (unit == "cluster2_pcid") {
    if (!all(c("cluster", "cluster2") %in% names(df))) return(data.frame())
    if (all(is.na(df$cluster2) | df$cluster2 == 0L)) return(data.frame())
    if (all(is.na(df$cluster)  | df$cluster  == 0L)) return(data.frame())
    pcid_pts <- df %>%
      dplyr::filter(!is.na(cluster2), cluster2 > 0L,
                    !is.na(cluster),  cluster  > 0L) %>%
      dplyr::group_by(cluster2, cluster) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    if (nrow(pcid_pts) < 2L) return(data.frame())
    rows <- lapply(split(pcid_pts, pcid_pts$cluster2), function(g) {
      d <- run_nn(as.matrix(g[, c("x","y","z")]), k)
      if (!length(d)) return(NULL)
      data.frame(distance = d)
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    return(if (is.null(out)) data.frame() else out)
  }

  data.frame()
}

safe_nn2 <- function(X, k = 4) {
  X <- as.matrix(X)
  k <- max(1, min(k, nrow(X)))
  out <- tryCatch(RANN::nn2(X, k = k), error = function(e) NULL)
  if (is.null(out)) return(list(nn = NULL, k_eff = 0))
  # Drop all-NA rows (shouldn’t happen, but be safe)
  if (any(!is.finite(out$nn.dists))) {
    bad <- !is.finite(rowSums(out$nn.dists))
    if (any(bad)) {
      out$nn.dists <- out$nn.dists[!bad, , drop = FALSE]
      out$nn.idx   <- out$nn.idx  [!bad, , drop = FALSE]
    }
  }
  list(nn = out, k_eff = ncol(out$nn.idx))
}

## Two-colour clustering used to live here as apply_two_color_clustering_pipeline(),
## a near line-by-line copy of apply_clustering_pipeline() with the group loop
## hardcoded to c(1L, 2L). It has been folded into that function: channels come
## from resolve_channels() and the loop covers 1..N, so 2-colour and multiplex
## are the same code path and neither can merge across a channel boundary.


## ── Two-color cross-channel NN ──────────────────────────────────────────────
compute_cross_nn <- function(df, unit = "tid", k = 1L) {
  if (!"analysis_group" %in% names(df)) return(data.frame())
  df1 <- dplyr::filter(df, analysis_group == 1L)
  df2 <- dplyr::filter(df, analysis_group == 2L)
  if (!nrow(df1) || !nrow(df2)) return(data.frame())

  get_pts <- function(d, u) {
    if (u == "tid") {
      d %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(x = median(x, na.rm = TRUE),
                         y = median(y, na.rm = TRUE),
                         z = median(z, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    } else {
      # "cluster2" collapses each Median Cluster to ONE point, so the distance is
      # structure-to-structure ("is my group-1 cluster near a group-2 cluster")
      # rather than molecule-to-molecule. Both channels keep their own clusters -
      # nothing here mixes the two groups.
      want   <- if (u == "cluster2") "cluster2" else "cluster"
      cl_col <- if (want %in% names(d) && any(!is.na(d[[want]]) & d[[want]] > 0L)) want else NULL
      if (is.null(cl_col)) return(data.frame(x = numeric(), y = numeric(), z = numeric()))
      d %>%
        dplyr::filter(!is.na(.data[[cl_col]]) & .data[[cl_col]] > 0L) %>%
        dplyr::group_by(.data[[cl_col]]) %>%
        dplyr::summarise(x = median(x, na.rm = TRUE),
                         y = median(y, na.rm = TRUE),
                         z = median(z, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    }
  }

  pts1 <- get_pts(df1, unit)
  pts2 <- get_pts(df2, unit)
  if (!nrow(pts1) || !nrow(pts2)) return(data.frame())

  k1 <- min(k, nrow(pts2))
  k2 <- min(k, nrow(pts1))

  # Use 2-D (x, y) distances so that the CSR null model (which is 2-D) is consistent.
  run_nn <- function(query, ref, kk) {
    nn <- RANN::nn2(as.matrix(ref[, c("x","y")]),
                    as.matrix(query[, c("x","y")]),
                    k = kk)
    as.vector(nn$nn.dists[, kk])
  }

  bbox_area <- diff(range(c(pts1$x, pts2$x), na.rm = TRUE)) *
               diff(range(c(pts1$y, pts2$y), na.rm = TRUE))
  rbind(
    data.frame(direction = "Group 1 → Group 2", distance = run_nn(pts1, pts2, k1),
               n_ref = nrow(pts2), area = bbox_area),
    data.frame(direction = "Group 2 → Group 1", distance = run_nn(pts2, pts1, k2),
               n_ref = nrow(pts1), area = bbox_area)
  )
}


# =============================================================================
# PANEL HELPERS (moved from server_smlm.R; plain functions of their arguments)
# =============================================================================

# GOAL: the options of the NN-line overlay and NN table, read once from the app.
#   nn_k            neighbours per point
#   nn_mplex_scope  "within" / "across" / "all" channel pairs
#   deblinked       de-blinking ran, so lines join molecule medians
#   distance_range  keep lines whose length lies in this range (nm)
#   within_cluster2 keep only lines inside one Median Cluster
#   label_offset    distance-label offset (nm)
smlm_nn_options <- function(input) {
  list(nn_k            = input[["nn_k"]] %||% 4,
       nn_mplex_scope  = input[["nn_mplex_scope"]] %||% "within",
       deblinked       = isTRUE(input[["enable_clustering"]]),
       distance_range  = input[["distance_range"]] %||% c(0, Inf),
       within_cluster2 = isTRUE(input[["nn_within_cluster2"]]) && isTRUE(input[["enable_cluster2"]]),
       label_offset    = suppressWarnings(as.numeric(input[["nn_label_offset"]] %||% 5)))
}

# NN index pairs honouring the channel scope radio. Reads `channel`, so one
# control covers a 2-colour split and a multiplex build alike; single-channel
# data falls back to a global k-NN. "within" = same-channel only, "across" =
# strictly between different channels (never within), "all" = both.
# Used by every NN-line / NN-table path.
#   opts : smlm_nn_options()
.nn_pairs_scoped <- function(df_geo, coords, opts) {
  k <- opts$nn_k %||% 4
  ch <- if ("channel" %in% names(df_geo)) df_geo$channel
        else if ("mplex_channel" %in% names(df_geo)) df_geo$mplex_channel else NULL
  if (is.null(ch) || length(unique(ch[!is.na(ch)])) < 2L)
    return(nn_index_pairs(coords, k))
  scope <- opts$nn_mplex_scope %||% "within"
  if (identical(scope, "all")) return(nn_index_pairs(coords, k))
  nn_index_pairs(coords, k, channel = ch,
                 mode = if (identical(scope, "across")) "between" else "within")
}

# core builder used by both 2D/3D
#   opts : smlm_nn_options()
build_nn_segments <- function(df, opts, is3d = TRUE) {
  if (!nrow(df)) return(NULL)
  
  # If analysis_group isn't present, create NA column (scope=="all" will still work)
  if (!"analysis_group" %in% names(df)) df$analysis_group <- NA_integer_
  
  scope <- "all"
  
  # 1) Choose medians: by cluster if enabled, else by tid
  if (isTRUE(opts$deblinked) && "cluster" %in% names(df)) {
    df_geo <- df %>%
      dplyr::filter(!is.na(cluster) & cluster != 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
        color = dplyr::first(color),
        cluster2 = if ("cluster2" %in% names(df)) dplyr::first(cluster2) else NA_integer_,
        analysis_group = suppressWarnings(dplyr::first(analysis_group)),
        channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
        channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
        channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
        mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
        .groups = "drop"
      )
  } else {
    df_geo <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
        color = dplyr::first(color),
        cluster2 = if ("cluster2" %in% names(df)) dplyr::first(cluster2) else NA_integer_,
        analysis_group = suppressWarnings(dplyr::first(analysis_group)),
        channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
        channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
        channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
        mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
        .groups = "drop"
      )
  }

  if (nrow(df_geo) < 2) return(NULL)

  coords <- as.matrix(df_geo[, c("x","y","z")])
  pm <- .nn_pairs_scoped(df_geo, coords, opts)
  if (!nrow(pm)) return(NULL)

  seg_x <- c(); seg_y <- c(); seg_z <- c()
  lab_x <- c(); lab_y <- c(); lab_z <- c(); labs <- c()
  e_i <- c(); e_j <- c(); e_d <- c()      # per-edge record, for channel colouring
  seen  <- new.env(parent = emptyenv())

  dmin <- (opts$distance_range %||% c(0, Inf))[1]
  dmax <- (opts$distance_range %||% c(0, Inf))[2]

  for (r in seq_len(nrow(pm))) {
    i <- pm[r, 1]; ni <- pm[r, 2]
    {
      # dedupe undirected
      key <- paste0(min(i, ni), "_", max(i, ni))
      if (!is.null(seen[[key]])) next
      seen[[key]] <- TRUE

      # within/between/all
      g1 <- df_geo$analysis_group[i]
      g2 <- df_geo$analysis_group[ni]
      if (!pair_passes_scope(g1, g2, scope)) next

      dist <- pm[r, 3]
      if (dist < dmin || dist > dmax) next

      # within-cluster2 filter
      if (isTRUE(opts$within_cluster2) && "cluster2" %in% names(df_geo)) {
        c2i <- df_geo$cluster2[i]; c2ni <- df_geo$cluster2[ni]
        if (is.na(c2i) || is.na(c2ni) || c2i == 0 || c2ni == 0 || c2i != c2ni) next
      }

      seg_x <- c(seg_x, coords[i,1], coords[ni,1], NA_real_)
      seg_y <- c(seg_y, coords[i,2], coords[ni,2], NA_real_)
      seg_z <- c(seg_z, coords[i,3], coords[ni,3], NA_real_)

      nn_off <- opts$label_offset %||% 5
      if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
      mid <- perp_offset_xy(coords[i,], coords[ni,], nn_off)
      lab_x <- c(lab_x, mid[1]); lab_y <- c(lab_y, mid[2]); lab_z <- c(lab_z, mid[3])
      labs  <- c(labs, round(dist, 2))
      e_i <- c(e_i, i); e_j <- c(e_j, ni); e_d <- c(e_d, dist)
    }
  }

  if (!length(seg_x)) return(NULL)

  # Per-edge table, so the plot can colour a line by the channels it joins.
  # This is what the separate 2-colour NN panel used to exist for: a
  # within-channel line drawn in that channel's colour, a cross-channel line in
  # magenta. It is the one thing that panel did which the shared one could not,
  # so it moved here rather than being dropped with it.
  edges <- data.frame(
    x1 = coords[e_i, 1], y1 = coords[e_i, 2], z1 = coords[e_i, 3],
    x2 = coords[e_j, 1], y2 = coords[e_j, 2], z2 = coords[e_j, 3],
    dist = e_d,
    ch1 = df_geo$channel[e_i], ch2 = df_geo$channel[e_j],
    lab1 = df_geo$channel_label[e_i], lab2 = df_geo$channel_label[e_j],
    stringsAsFactors = FALSE
  )
  same <- !is.na(edges$ch1) & !is.na(edges$ch2) & edges$ch1 == edges$ch2
  col_same <- df_geo$channel_color[e_i]
  edges$color <- ifelse(same,
                        ifelse(is.na(col_same), "#666666", col_same),
                        "magenta")

  list(
    seg_x = seg_x, seg_y = seg_y, seg_z = seg_z,
    lab_x = lab_x, lab_y = lab_y, lab_z = lab_z, labs = labs,
    edges = edges
  )
}

# Shift midpoint of segment (p1→p2) perpendicularly in XY plane by offset_nm.
# Works for both 2-element (x,y) and 3-element (x,y,z) vectors.
perp_offset_xy <- function(p1, p2, offset_nm) {
  p1 <- as.numeric(p1); p2 <- as.numeric(p2)
  mid <- (p1 + p2) / 2
  if (!is.finite(offset_nm) || offset_nm <= 0) return(mid)
  dx <- p2[1] - p1[1]; dy <- p2[2] - p1[2]
  len <- sqrt(dx^2 + dy^2)
  if (len < 1e-10) return(mid)
  perp <- c(-dy / len, dx / len)
  if (length(mid) == 3) perp <- c(perp, 0)
  mid + perp * offset_nm
}

# SINGLE SOURCE OF TRUTH for "how close are the members of a median cluster",
# used by BOTH the Crowding panel and the Within-cluster-NN-vs-size panel so the
# two report the SAME number. Both panels build their member centroids on the
# same unit (de-blinked molecule when de-blinking is on, else trace ID; see
# smlm_nn_unit below) and pass the member coordinates here.
#
# For each member, mean distance to its k nearest members (k clamped to n-1;
# k=1 => nearest-neighbour spacing), then the MEDIAN over members with a robust
# spread (MAD and the inter-quartile range) as the ± error.
# pts = numeric matrix of member coordinates (one row per member).
member_nn_stats <- function(pts, k = 1L) {
  n <- nrow(pts)
  if (is.null(n) || n < 2L) return(NULL)
  k_eff <- min(max(1L, as.integer(k)), n - 1L)
  nn <- tryCatch(dbscan::kNN(pts, k = k_eff)$dist, error = function(e) NULL)
  if (is.null(nn)) return(NULL)
  spacing <- if (k_eff == 1L) as.numeric(nn) else rowMeans(nn, na.rm = TRUE)
  qs <- as.numeric(stats::quantile(spacing, c(0.25, 0.75), na.rm = TRUE))
  list(median = stats::median(spacing, na.rm = TRUE),
       mad    = stats::mad(spacing, na.rm = TRUE),
       q1 = qs[1], q3 = qs[2], k_used = k_eff, n = n)
}

