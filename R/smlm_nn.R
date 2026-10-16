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


# =============================================================================
# batch_nn_for_file() - ONE FILE'S NEAREST-NEIGHBOUR DISTANCES FOR BATCH NN
# =============================================================================
# GOAL: the Batch NN panel's work for one file that has already been filtered
# (and thinned): cluster or group it as the chosen comparison needs, facet it by
# channel when asked, and compute the k-th NN distances per facet.
#
# WHICH PREPARATION
#   multiplex pair  cluster if the unit needs it, then channel X -> group 1 and
#                   channel Y -> group 2, for a between-channel distance
#   two-colour with a group filter or cross-group comparison: the pipeline,
#                   whose channel split sets analysis_group
#   molecule unit   the pipeline (de-blinking)
#   group filter only (trace unit): assign_analysis_group
#   A cross-channel run is never faceted: it is a statistic about the pair.
#
# ARGUMENTS
#   df            filtered localization table of one file
#   fname         dataset name (facet names are built from it)
#   settings      smlm_settings or input (the run's snapshot)
#   unit          "tid", "pcid", "cluster2_tid" or "cluster2_pcid"
#   k             neighbour order
#   group_filter  1L / 2L to stay within one group, NA for all
#   cross_groups  c(from, to) for between-group distances, or NULL
#   mplex_pair    c(channel X, channel Y) labels for a multiplex between run, or NULL
#   by_channel    facet by channel (the "analyse each channel" switch)
#
# RETURNS  named list of data.frames (distance, dataset, tag, ...), one per facet
#          that produced distances; empty list when none did.
# =============================================================================
batch_nn_for_file <- function(df, fname, settings, unit = "tid", k = 1L,
                              group_filter = NA_integer_, cross_groups = NULL,
                              mplex_pair = NULL, by_channel = FALSE) {
  s <- as_smlm_settings(settings)
  needs_clustering <- unit != "tid"
  needs_grouping   <- !is.na(group_filter) || !is.null(cross_groups)
  if (!is.null(mplex_pair)) {
    if (needs_clustering)
      df <- tryCatch(apply_clustering_pipeline(df, settings), error = function(e) df)
    if (is.null(df) || !nrow(df)) return(list())
    df$analysis_group <- dplyr::case_when(
      df$mplex_label == mplex_pair[1] ~ 1L,
      df$mplex_label == mplex_pair[2] ~ 2L,
      TRUE ~ NA_integer_)
  } else if (isTRUE(s$enable_two_color) && needs_grouping) {
    df <- tryCatch(apply_clustering_pipeline(df, settings), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(list())
  } else if (needs_clustering) {
    df <- tryCatch(apply_clustering_pipeline(df, settings), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(list())
  } else if (needs_grouping) {
    df <- tryCatch(assign_analysis_group(df, settings), error = function(e) df)
  }
  facets <- if (!is.null(mplex_pair) || needs_grouping) list(list(label = "", df = df))
            else channel_facets(df, by_channel)
  out <- list()
  for (fc in facets) {
    d <- tryCatch(compute_batch_nn(fc$df, unit = unit, k = k, group_filter = group_filter,
                                   cross_groups = cross_groups),
                  error = function(e) data.frame())
    if (!nrow(d)) next
    nm2 <- .facet_name(fname, fc$label)
    d$dataset <- nm2
    d$tag     <- if (".dataset_tag" %in% names(df)) df$.dataset_tag[[1L]] else fname
    out[[nm2]] <- d
  }
  out
}

# =============================================================================
# crowding_for_file() - CROWDING VS DETECTION-TIME SHIFT, ONE FILE
# =============================================================================
# GOAL: for every Median Cluster (cluster2) structure in one prepared file, the
# member spacing (median k-NN distance among member centroids) and the
# detection-time shift (median member detection time minus the median of the
# unclustered background). Optionally the no-crowding null: detection times
# shuffled among the units of this file B times, geometry and membership kept.
#
# Members are de-blinked molecules when id_col = "cluster" (primary noise dropped),
# else trace IDs: the same unit as the Within-cluster-NN panel. A member's
# detection time is the median tim of its localizations.
#
# ARGUMENTS
#   df      prepared localization table (filter -> thin -> cluster chain), with
#           tim and cluster2
#   fname   dataset name written into the rows
#   id_col  member unit: "cluster" (de-blinked molecules, the server's
#           smlm_nn_unit()) or "tid"
#   min_n   members a structure needs
#   k       neighbour order for the spacing
#   use_3d  use z in the spacing
#   B       permutations for the null (0 = no null)
#
# RETURNS list(rows = data.frame or NULL, null = structures x B matrix or NULL,
#              diag = counts of why structures or the file were skipped).
#          Draws random numbers when B > 0.
# =============================================================================
crowding_for_file <- function(df, fname, id_col = "tid", min_n = 4L, k = 1L,
                              use_3d = FALSE, B = 0L) {
  diag <- c(no_c2 = 0L, no_z = 0L, no_clusters = 0L, no_background = 0L, too_small = 0L)
  out  <- list(rows = NULL, null = NULL, diag = diag)
  # NOTE: the column is `cluster2`. "Median cluster" is only the display name.
  if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) { out$diag["no_c2"] <- 1L; return(out) }
  if (use_3d && !"z" %in% names(df)) { out$diag["no_z"] <- 1L; return(out) }
  cols <- if (use_3d) c("x", "y", "z") else c("x", "y")

  # one row per member unit: detection time, centroid, structure id
  base   <- if (id_col == "cluster") dplyr::filter(df, !is.na(cluster) & cluster > 0L) else df
  per_tid <- base %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      tim_med = stats::median(tim, na.rm = TRUE),
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
      cl2 = { v <- cluster2[!is.na(cluster2) & cluster2 > 0L]; if (length(v)) v[1] else 0L },
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(tim_med))
  if (!nrow(per_tid)) return(out)

  cids <- unique(per_tid$cl2[per_tid$cl2 > 0L])
  if (!length(cids)) { out$diag["no_clusters"] <- 1L; return(out) }

  bg_med <- stats::median(per_tid$tim_med[per_tid$cl2 == 0L], na.rm = TRUE)
  if (!is.finite(bg_med)) {        # need an unclustered background to reference
    out$diag["no_background"] <- 1L; return(out)
  }

  rows <- list()
  keep_idx <- list()   # row indices of each qualifying cluster, for the permutation
  for (cid in cids) {
    which_c <- which(per_tid$cl2 == cid)
    sub <- per_tid[which_c, , drop = FALSE]
    n   <- nrow(sub)
    if (n < min_n) { out$diag["too_small"] <- out$diag["too_small"] + 1L; next }
    pts <- as.matrix(sub[, cols, drop = FALSE])
    pts <- pts[stats::complete.cases(pts), , drop = FALSE]
    if (nrow(pts) < 2L) next
    st <- member_nn_stats(pts, k = k)
    if (is.null(st)) next
    rows[[length(rows) + 1L]] <- data.frame(
      dataset    = fname,
      cluster2   = cid,
      n_tids     = n,
      nn_median  = st$median,
      nn_mad     = st$mad,
      unit       = if (id_col == "cluster") "molecule" else "tid",
      time_shift = stats::median(sub$tim_med) - bg_med,
      stringsAsFactors = FALSE
    )
    keep_idx[[length(keep_idx) + 1L]] <- which_c
  }
  if (length(rows)) out$rows <- do.call(rbind, rows)

  # Permute detection times WITHIN this file; geometry and membership untouched.
  if (B > 0L && length(keep_idx)) {
    tv     <- per_tid$tim_med
    is_bg  <- per_tid$cl2 == 0L
    nullm  <- matrix(NA_real_, nrow = length(keep_idx), ncol = B)
    for (b in seq_len(B)) {
      tv_b   <- tv[sample.int(length(tv))]          # shuffle times across all units
      bg_b   <- stats::median(tv_b[is_bg], na.rm = TRUE)
      nullm[, b] <- vapply(keep_idx,
                           function(ix) stats::median(tv_b[ix], na.rm = TRUE) - bg_b,
                           numeric(1))
    }
    out$null <- nullm
  }
  out
}

# =============================================================================
# crowding_null_test() - NULL BAND AND TREND TEST FOR THE CROWDING PANEL
# =============================================================================
# GOAL: attach each structure's 95% no-crowding interval (null_lo, null_med,
# null_hi) and test the TREND: Spearman rho between spacing and time shift,
# observed against its permutation distribution. Rank-based, so a constant
# clustered-vs-background offset and a skewed shift distribution do not move it.
#
# ARGUMENTS
#   res  structure rows (crowding_for_file rows, all files bound)
#   NM   structures x B matrix of permuted time shifts, rows matching res; NULL
#        when no null was computed
#
# RETURNS list(res = res with the null columns (NA without a null),
#              stats = list(rho_obs, p, B, outside) or NULL)
# =============================================================================
crowding_null_test <- function(res, NM = NULL) {
  if (is.null(NM) || nrow(NM) != nrow(res)) {
    res$null_lo <- NA_real_; res$null_med <- NA_real_; res$null_hi <- NA_real_
    return(list(res = res, stats = NULL))
  }
  qs <- t(apply(NM, 1, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  res$null_lo  <- qs[, 1]; res$null_med <- qs[, 2]; res$null_hi <- qs[, 3]
  rho_of <- function(y) suppressWarnings(
    stats::cor(res$nn_median, y, method = "spearman", use = "complete.obs"))
  rho_obs  <- rho_of(res$time_shift)
  rho_null <- apply(NM, 2, rho_of)
  rho_null <- rho_null[is.finite(rho_null)]
  p_two <- if (is.finite(rho_obs) && length(rho_null))
    (1 + sum(abs(rho_null) >= abs(rho_obs))) / (length(rho_null) + 1) else NA_real_
  outside <- mean(res$time_shift < res$null_lo | res$time_shift > res$null_hi, na.rm = TRUE)
  list(res = res, stats = list(rho_obs = rho_obs, p = p_two, B = length(rho_null),
                               outside = outside))
}

# =============================================================================
# size_nn_for_file() - WITHIN-STRUCTURE NN SPACING VS STRUCTURE SIZE, ONE FILE
# =============================================================================
# GOAL: one row per Median Cluster (cluster2) structure: size (member count),
# median k-NN spacing among member centroids (+- MAD, the same definition as the
# Crowding panel) and convex-hull area (2D) or volume (3D). Tests whether the
# spacing stays constant as structures grow or "balloons".
#
# ARGUMENTS
#   df      prepared localization table with cluster2
#   fname   dataset name written into the rows
#   id_col  member unit: "cluster" (de-blinked molecules) or "tid"
#   min_n   members a structure needs
#   k       neighbour order for the spacing
#   use_3d  use z for the spacing and the hull
#
# RETURNS list(rows = data.frame or NULL, diag = counts of why structures or the
#              file were skipped)
# =============================================================================
size_nn_for_file <- function(df, fname, id_col = "tid", min_n = 3L, k = 1L, use_3d = FALSE) {
  diag <- c(no_c2 = 0L, no_z = 0L, no_clusters = 0L, too_small = 0L, degenerate = 0L)
  out  <- list(rows = NULL, diag = diag)
  # NOTE: the column is `cluster2`. "Median Cluster" is only the display name.
  if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) { out$diag["no_c2"] <- 1L; return(out) }
  if (use_3d && !"z" %in% names(df)) { out$diag["no_z"] <- 1L; return(out) }
  cols <- if (use_3d) c("x", "y", "z") else c("x", "y")

  # member centroids per structure
  cent <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 > 0L) %>%
    { if (id_col == "cluster") dplyr::filter(., !is.na(cluster) & cluster > 0L) else . } %>%
    dplyr::group_by(cluster2, .data[[id_col]]) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z = if ("z" %in% names(df)) stats::median(z, na.rm = TRUE) else 0,
      .groups = "drop")
  if (!nrow(cent)) { out$diag["no_clusters"] <- 1L; return(out) }

  rows <- list()
  for (cid in unique(cent$cluster2)) {
    sub <- cent[cent$cluster2 == cid, , drop = FALSE]
    n   <- nrow(sub)
    if (n < min_n) { out$diag["too_small"] <- out$diag["too_small"] + 1L; next }
    pts <- as.matrix(sub[, cols, drop = FALSE])
    pts <- pts[stats::complete.cases(pts), , drop = FALSE]
    if (nrow(pts) < 2L) { out$diag["degenerate"] <- out$diag["degenerate"] + 1L; next }
    st <- member_nn_stats(pts, k = k)
    if (is.null(st)) next
    # convex-hull size: area (2D) / volume (3D); NA if degenerate
    hull <- tryCatch(geometry::convhulln(pts, options = "FA")$vol, error = function(e) NA_real_)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset   = fname, cluster2 = cid, size = n,
      nn_median = st$median,
      nn_mad    = st$mad,
      hull_size = if (is.finite(hull) && hull > 0) hull else NA_real_,
      dim       = length(cols),
      unit      = if (id_col == "cluster") "molecule" else "tid",
      stringsAsFactors = FALSE)
  }
  if (length(rows)) out$rows <- do.call(rbind, rows)
  out
}
