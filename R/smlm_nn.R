## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_nn.R  -  Nearest-neighbour distances
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   k-NN distances within and between channels and per Median Cluster, for
##   the NN and cross-channel NN panels.
##
## CONTENTS (in file order)
##   nn_index_pairs, compute_batch_nn, safe_nn2, compute_cross_nn
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

