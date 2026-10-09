## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_clustering.R  -  De-blinking and median clustering
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   apply_clustering_pipeline (the entry point) with its two tiers:
##   de-blinking (complete linkage on trace centroids, merge-only) and median
##   clustering (DBSCAN / HDBSCAN / OPTICS on unit medians), both strictly per
##   channel; plus the merge-distance and drift diagnostics.
##
## CONTENTS (in file order)
##   run_cluster2_grouping, .prune_sparse_members, .cluster2_core,
##   deblink_axis_weights, apply_deblink_weights, deblink_coord_cols,
##   matrixStats_colMedian, link_complete_threshold, deblink_merge_centroids,
##   suggest_centroid_threshold, deblink_threshold_for_channel,
##   sample_drift_estimate, deblink_merge_majority_vote, deblink_minpts,
##   deblink_assign_singletons, deblink_count_curve, apply_clustering_pipeline
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

## ── Secondary (cluster2) grouping engine ──────────────────────────────────────
# Runs the secondary grouping on a matrix of centroid coordinates (trace or
# de-blinked molecule medians) and returns an integer vector of cluster ids using
# the dbscan convention: 0 = noise/unclustered, 1..k = clusters. This keeps the
# output identical in shape regardless of method, so every downstream consumer
# (cluster2 columns, isosurfaces, fitting pipelines, NN, exports) is unaffected.
#
#   method = "dbscan"  -> density clustering with fixed radius `eps` + `minPts`
#   method = "hdbscan" -> hierarchical density clustering, `minPts` only (no eps)
#   method = "optics"  -> OPTICS reachability + Xi extraction; `minPts` + `xi`,
#                         adaptive density, keeps small structures (no big eps)
#
# Optional local-density prefilter (prefilter = TRUE): before the chosen method
# runs, any centroid with fewer than `prefilter_min` neighbours within radius
# `prefilter_r` is forced to noise. This encodes "structures have higher local
# TID density than scattered background" as an explicit, tunable gate, and works
# in front of any method.
#
# Three membership refinements run AFTER the method, in this order:
#
#   core_only    (DBSCAN only) keep only points that are themselves core, i.e.
#                have >= minPts neighbours within eps. Classic DBSCAN also
#                absorbs BORDER points - anything within eps of a core point,
#                however isolated it is itself - which is what attaches lone
#                background centroids to the rim of a real structure. Dropping
#                them is the DBSCAN* cluster definition of Campello et al.
#                (2013), the same definition HDBSCAN uses, so HDBSCAN/OPTICS
#                need no equivalent and ignore this flag.
#
#   prune        graded version of the same idea - see .prune_sparse_members().
#
#   min_members  clusters left with fewer than this many members are dropped
#                whole. Applied LAST so it counts what actually survived, and
#                separate from minPts so the density scale and the "how big is
#                a real structure" decision are not the same knob.
#
# `eps` is used only by dbscan; `xi` only by optics. Works in 2D or 3D depending
# on ncol(pts).
run_cluster2_grouping <- function(pts, method = "dbscan", eps = 200, minPts = 15,
                                  xi = 0.05,
                                  prefilter = FALSE, prefilter_r = NULL,
                                  prefilter_min = NULL,
                                  prune = FALSE, prune_tol = 2,
                                  core_only = FALSE, min_members = 0) {
  if (is.null(pts)) return(integer(0))
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  n <- nrow(pts)
  if (n < 1L) return(integer(0))

  method <- tolower(as.character(method %||% "dbscan"))

  # ── optional local-density prefilter ──
  keep <- rep(TRUE, n)
  if (isTRUE(prefilter) && !is.null(prefilter_r) && isTRUE(prefilter_r > 0)) {
    pf_min <- max(1L, as.integer(round(prefilter_min %||% 1L)))
    nn <- tryCatch(dbscan::frNN(pts, eps = prefilter_r, sort = FALSE),
                   error = function(e) NULL)
    if (!is.null(nn)) keep <- lengths(nn$id) >= pf_min   # neighbours exclude self
  }

  out <- rep(0L, n)
  if (!any(keep)) return(out)
  sub <- pts[keep, , drop = FALSE]
  lab <- .cluster2_core(sub, method, eps, minPts, xi)

  # ── core-only membership (DBSCAN*) ──
  # Only DBSCAN has border points to drop: HDBSCAN/OPTICS are already defined on
  # core points, and neither has a fixed eps to test against.
  if (isTRUE(core_only) && identical(method, "dbscan")) {
    mp <- max(2L, as.integer(round(minPts)))
    nn <- tryCatch(dbscan::frNN(sub, eps = eps, sort = FALSE), error = function(e) NULL)
    if (!is.null(nn)) {
      is_core <- (lengths(nn$id) + 1L) >= mp        # +1: frNN excludes self, minPts includes it
      lab[!is_core] <- 0L
    }
  }

  # ── optional sparse-member pruning (post-clustering) ──
  if (isTRUE(prune)) lab <- .prune_sparse_members(sub, lab, tol = prune_tol, minPts = minPts)

  # ── minimum members per cluster, on what survived ──
  # as.numeric() first: a restored settings file can deliver this as a string,
  # and round() errors on character rather than coercing.
  mm <- suppressWarnings(as.numeric(min_members %||% 0))
  mm <- if (is.finite(mm)) as.integer(round(mm)) else 0L
  if (is.finite(mm) && mm > 1L) {
    tb <- table(lab[lab > 0])
    small <- as.integer(names(tb)[tb < mm])
    if (length(small)) lab[lab %in% small] <- 0L
  }

  # Renumber survivors to a contiguous 1..k so the documented output contract
  # holds after members and whole clusters have been removed.
  pos <- lab > 0
  if (any(pos)) lab[pos] <- as.integer(factor(lab[pos]))

  out[keep] <- lab
  out
}

# Post-clustering refinement: within each cluster, demote to noise (0) any
# member sitting in a much sparser spot than the members around it.
#
# STATISTIC: each member's CORE DISTANCE - the distance to its (minPts-1)'th
# nearest fellow member, i.e. the radius it needs to reach minPts neighbours.
# Small = deep inside a density maximum, large = on the outskirts. This judges
# density at the SAME scale the clustering used (minPts), and it is a density
# statistic, not a Euclidean distance to a centre or to a mode - so it behaves
# on rings, arcs, rods and helices, where "far from the middle" is a property of
# legitimate members rather than of stragglers.
#
# BASELINE: the median core distance of the member's OWN neighbourhood, not of
# the whole cluster. A single per-cluster median is set by the densest region,
# so on any structure with an internal density gradient it condemns the sparser
# half wholesale even where those members are interior. A local baseline moves
# with the gradient and only flags members that are loose RELATIVE TO THEIR
# SURROUNDINGS.
#
# `tol` is the only knob and the cutoff stays relative (no absolute nm): a
# member is dropped when its core distance exceeds `tol` x its neighbourhood's
# median (tol ~ 2 = "twice as loose as the members around me").
.prune_sparse_members <- function(pts, labels, tol = 2, minPts = 15) {
  if (!length(labels)) return(labels)
  tol <- suppressWarnings(as.numeric(tol %||% 2))
  if (!is.finite(tol) || tol <= 1) return(labels)   # tol<=1 would trim the core too
  mp  <- max(2L, as.integer(round(minPts %||% 15)))
  out <- labels
  for (cl in unique(labels[labels > 0])) {
    idx <- which(labels == cl)
    nm  <- length(idx)
    if (nm < 4L) next                                # no baseline from <4 members
    sub <- pts[idx, , drop = FALSE]

    k  <- max(1L, min(mp - 1L, nm - 1L))             # core distance at the clustering scale
    kd <- tryCatch(dbscan::kNN(sub, k = k)$dist, error = function(e) NULL)
    if (is.null(kd)) next
    cd <- kd[, k]

    m   <- max(4L, min(mp, nm - 1L))                 # baseline neighbourhood
    nbr <- tryCatch(dbscan::kNN(sub, k = m)$id, error = function(e) NULL)
    if (is.null(nbr)) next
    ref <- apply(nbr, 1, function(j) stats::median(cd[j]))   # neighbours', not own

    drop <- is.finite(ref) & ref > 0 & cd > tol * ref
    out[idx[drop]] <- 0L
  }
  out
}

# Single-method core; assumes points already prefiltered. Returns 0=noise/1..k.
.cluster2_core <- function(pts, method, eps, minPts, xi) {
  n <- nrow(pts)
  if (n < 1L) return(integer(0))

  if (identical(method, "hdbscan")) {
    # hdbscan needs minPts >= 2 and at least minPts points to form any cluster.
    mp <- max(2L, as.integer(round(minPts)))
    if (n < mp) return(rep(0L, n))
    res <- tryCatch(dbscan::hdbscan(pts, minPts = mp), error = function(e) NULL)
    if (is.null(res) || is.null(res$cluster)) return(rep(0L, n))
    return(as.integer(res$cluster))
  }

  if (identical(method, "optics")) {
    mp <- max(2L, as.integer(round(minPts)))
    if (n < mp) return(rep(0L, n))
    oo <- tryCatch(dbscan::optics(pts, minPts = mp), error = function(e) NULL)
    if (is.null(oo)) return(rep(0L, n))
    xv <- min(max(as.numeric(xi %||% 0.05), 1e-3), 0.999)
    cx <- tryCatch(dbscan::extractXi(oo, xi = xv), error = function(e) NULL)
    if (is.null(cx) || is.null(cx$cluster)) return(rep(0L, n))
    return(as.integer(cx$cluster))
  }

  # default: classic DBSCAN
  as.integer(dbscan::dbscan(pts, eps = eps, minPts = minPts)$cluster)
}

# =============================================================================
# deblink_axis_weights() - PUT X, Y AND Z ON A COMMON PRECISION SCALE
# =============================================================================
#
# GOAL
#   DBSCAN uses a spherical eps, but MINFLUX trace spread is not isotropic - in
#   this app z is additionally multiplied by the axial calibration (default 0.7,
#   a refractive-index mismatch correction, user-set per dataset). A single
#   radius therefore cannot fit all three axes: set it for xy and traces
#   fragment along the tight axis, set it for z and xy over-merges.
#
# WHY RESCALING IS CORRECT *FOR DE-BLINKING* SPECIFICALLY
#   "Are these two localizations the same fluorophore?" is a statistical
#   question, not a geometric one. If the axial scatter is smaller, the same
#   displacement in z is MORE evidence of two distinct molecules than in xy.
#   Weighting each axis by its own spread is the Mahalanobis form of that test.
#
#   This is emphatically NOT valid for measuring physical size - cluster radii,
#   Rg, nearest-neighbour distances must stay in true nm. The weights are for
#   the merge decision only.
#
# UNITS ARE PRESERVED
#   Weights are normalised by the geometric mean of the per-axis spreads, so the
#   transform is volume-preserving and eps keeps its nm meaning. Perfectly
#   isotropic data yields w = (1,1,1) and behaviour is identical to before,
#   which keeps saved settings comparable.
#
# ARGUMENTS
#   df    : data.frame with tid, x, y [, z]
#   min_n : minimum localizations for a trace to contribute a spread estimate
#   use_z : FALSE forces wz = 1 (2D data carries z = 0 and must not be rescaled)
#
# RETURNS
#   list(w = c(wx, wy, wz), sigma = c(sx, sy, sz), has_z, ok)
#   `w` is ALWAYS length 3 so callers can apply it to an (x,y,z) matrix
#   unconditionally. ok = FALSE means the spreads were unusable and w is
#   c(1,1,1) - the caller then behaves exactly as it did before.
# =============================================================================
deblink_axis_weights <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, use_z = TRUE) {
  flat <- list(w = c(1, 1, 1), sigma = c(NA_real_, NA_real_, NA_real_),
               has_z = FALSE, ok = FALSE)
  if (!is.data.frame(df) || !nrow(df)) return(flat)

  has_z <- isTRUE(use_z) && "z" %in% names(df) &&
    any(is.finite(df$z)) && stats::sd(df$z, na.rm = TRUE) > 1e-6

  sds <- tryCatch(sd_by_tid(df, min_n = min_n), error = function(e) NULL)
  if (!is.data.frame(sds) || !nrow(sds)) return(flat)

  med <- function(v) {
    v <- v[is.finite(v) & v > 0]
    if (length(v) < 3L) NA_real_ else stats::median(v)
  }
  s <- c(med(sds$sd_x), med(sds$sd_y), if (has_z) med(sds$sd_z) else NA_real_)

  # xy must both be measurable; z only when the data is 3D.
  if (!all(is.finite(s[1:2]))) return(flat)
  if (has_z && !is.finite(s[3])) has_z <- FALSE

  use <- if (has_z) 1:3 else 1:2
  gm  <- exp(mean(log(s[use])))                    # geometric mean -> volume-preserving
  if (!is.finite(gm) || gm <= 0) return(flat)

  w <- c(1, 1, 1)
  w[use] <- gm / s[use]
  list(w = w, sigma = s, has_z = has_z, ok = TRUE)
}

# GOAL: apply per-axis de-blink weights to an (x,y,z) coordinate matrix.
# Kept as one function so every de-blink call site transforms identically
# - the suggestion and the clustering must share a metric or eps is meaningless.
apply_deblink_weights <- function(pts, w) {
  if (is.null(w)) return(pts)
  # A weight vector shorter than the coordinate matrix used to pad with NA
  # (w[seq_len(3)] on a length-2 w), which poisons every centroid downstream and
  # silently switched de-blinking off. Anything missing or non-finite means
  # "leave that axis alone", which is the isotropic default.
  ww <- rep(1, ncol(pts))
  n  <- min(length(w), ncol(pts))
  if (n) ww[seq_len(n)] <- as.numeric(w)[seq_len(n)]
  ww[!is.finite(ww)] <- 1
  if (all(ww == 1)) return(pts)
  sweep(pts, 2L, ww, "*")
}

# GOAL: the coordinate columns actually present, in canonical order.
# The de-blink family must agree on this: suggest_centroid_threshold() and
# deblink_count_curve() already tolerated 2D input while the two merge functions
# hard-coded c("x","y","z") and errored with "undefined columns selected".
deblink_coord_cols <- function(df) intersect(c("x", "y", "z"), names(df))

# GOAL: per-column median of a coordinate matrix -> one numeric vector.
# Base R only (no matrixStats dependency); the matrices here are 2-3 columns and
# a handful of rows per trace, so a column loop costs nothing. Kept as its own
# function because every de-blink call site must collapse a trace to a centroid
# the same way - median, NA-tolerant, always length ncol even for empty input.
matrixStats_colMedian <- function(m) {
  if (!is.matrix(m)) m <- as.matrix(m)
  storage.mode(m) <- "double"
  if (!nrow(m)) return(rep(NA_real_, ncol(m)))
  vapply(seq_len(ncol(m)), function(j) stats::median(m[, j], na.rm = TRUE),
         numeric(1))
}

# =============================================================================
# link_complete_threshold() - COMPLETE LINKAGE AT A DISTANCE, AT SCALE
# =============================================================================
#
# GOAL
#   Group points so that EVERY pair inside a group is within `thr`, not merely
#   connected through a chain of near neighbours.
#
# WHY NOT JUST hclust()
#   `hclust(dist(M))` needs the full n^2 distance matrix - at 60k traces that is
#   ~1.8e9 pairs and will not run. The way round it is exact, not an
#   approximation: two points in different connected components of the
#   threshold graph are more than `thr` apart, so complete linkage can never
#   merge them. Cutting the global dendrogram at `thr` therefore gives exactly
#   the same partition as cutting each component's dendrogram separately. The
#   components are tiny at de-blinking thresholds, so the n^2 cost is paid only
#   within each of them.
#
# THE GUARD
#   A pathological threshold (or genuinely dense data) can still produce one
#   huge component. Above `max_comp` members the component falls back to its
#   single-linkage answer rather than attempting an allocation that would take
#   the session down - the result is then locally chain-prone, which is the old
#   behaviour, and never an error.
#
# ARGUMENTS
#   M        : numeric matrix, one row per object
#   thr      : distance threshold, same units as M
#   max_comp : component size above which single linkage is used instead
#
# RETURNS  integer group id per ROW of M, all > 0.
# =============================================================================
link_complete_threshold <- function(M, thr, max_comp = 2000L) {
  n <- nrow(M)
  if (!n) return(integer(0))
  if (n == 1L) return(1L)

  # Connected components of the threshold graph, by union-find over frNN.
  fr <- dbscan::frNN(M, eps = thr)
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  for (i in seq_len(n)) {
    nb <- fr$id[[i]]
    if (!length(nb)) next
    for (j in nb) { ra <- find(i); rb <- find(j); if (ra != rb) parent[rb] <- ra }
  }
  comp <- vapply(seq_len(n), find, integer(1))

  out <- integer(n); nxt <- 0L
  for (ix in split(seq_len(n), comp)) {
    if (length(ix) == 1L) { nxt <- nxt + 1L; out[ix] <- nxt; next }
    if (length(ix) > max_comp) {          # too big to resolve - keep the chain
      nxt <- nxt + 1L; out[ix] <- nxt; next
    }
    h <- stats::cutree(stats::hclust(stats::dist(M[ix, , drop = FALSE]),
                                     method = "complete"), h = thr)
    out[ix] <- nxt + as.integer(h); nxt <- nxt + max(h)
  }
  out
}

# =============================================================================
# deblink_merge_centroids() - MERGE TRACES ON CENTROID AGREEMENT
# =============================================================================
#
# GOAL
#   Decide "same fluorophore?" by comparing where the traces ARE, not whether
#   their point clouds touch. Each trace collapses to its median position and
#   single-linkage joins centroids closer than `thr`.
#
# WHY THIS BEATS RAW-POINT DBSCAN (measured, not assumed)
#   Raw-point single-linkage merges when ANY localization of one cloud comes
#   within eps of ANY localization of another, so the outlying tails do the
#   chaining and the effective resolution degrades as ~eps + 2*sigma*sqrt(2 ln n)
#   - it gets WORSE as traces get longer. Centroids have no tails.
#
#   On ground-truth simulation (360 molecules, 923 traces, sigma 6.5 nm):
#     raw-point DBSCAN   best at 0.62*sigma -> 35 errors (splits + merges)
#     centroid linkage   best at 1.54*sigma -> 10 errors
#   The reason is visible in the distance distributions: two re-blinks of one
#   fluorophore have centroids a median 3.1 nm apart (95th pct 7.1 nm ~ 1.1
#   sigma), while the nearest DIFFERENT molecule's centroid sits ~54 nm away.
#   Those two populations barely overlap; raw-point linkage discards that.
#
# THRESHOLD SCALE
#   ~1.5x the median per-axis trace spread. Note this is NOT the same number as
#   the raw-point eps (~0.6 sigma) - the two methods measure different distances,
#   which is why they have separate inputs rather than sharing one.
#
# NO NOISE BY CONSTRUCTION
#   minPts = 1, so every centroid is a core point and every trace lands in a
#   cluster. That matches the merge-only rule: de-blinking never removes a
#   molecule (see deblink_assign_singletons).
#
# COMPLETE LINKAGE IS THE DEFAULT (measured on a real imaging file, 2026-08-09)
#   Single linkage merges through chains, so a group's diameter is not bounded
#   by `thr`. That is visible in the real data without any ground truth - on
#   250218 imaging file A at k = 2.5 (thr = 15.2 nm) single linkage produces a
#   group spanning 27.5 nm, 1.8x the threshold and 4.5 sigma across; at k = 5 it
#   spans 73.2 nm. Complete linkage caps the diameter at `thr` by construction
#   (15.1 nm and 30.1 nm for the same two settings).
#
#   Semi-synthetic ground truth built ON that file's own point cloud (real
#   per-trace residuals, real molecule positions, 50% of molecules given a
#   second blink), best F1 over k:
#
#     injected motion   1 nm   2 nm   4 nm   6 nm   9 nm
#     single            0.954  0.944  0.914  0.872  0.786
#     complete          0.955  0.953  0.946  0.920  0.868
#
#   On a fixed, fiducial-corrected sample the motion term is small and the two
#   are nearly tied - the gain at 1-2 nm is within noise. What complete linkage
#   buys at EVERY motion level is tolerance to an over-large threshold, which is
#   the error a user makes when correcting a threshold by hand, which is to
#   overshoot. At 2 nm motion, F1 across k = 3.0 to
#   6.0 falls 0.928 -> 0.715 under single linkage and only 0.952 -> 0.843 under
#   complete. It never scored WORSE at any motion level tested.
#
# ARGUMENTS
#   sub_df  : rows to cluster (needs tid, x, y, z)
#   thr     : centroid-distance threshold in nm (weighted units if w is given)
#   w       : optional per-axis de-blink weights, same metric as everything else
#   linkage : "complete" - every member within thr of every other. The only
#             behaviour the app offers; there is no UI control for this.
#             "single" is the historical chain-through-neighbours rule, kept
#             REACHABLE FROM CODE ONLY so analyses run before 2026-08-09 can be
#             reproduced exactly. It lost to complete linkage at every motion
#             level tested (table above) - do not re-expose it.
#
# RETURNS  integer cluster id per ROW of sub_df, all > 0.
# =============================================================================
deblink_merge_centroids <- function(sub_df, thr, w = NULL,
                                    linkage = c("complete", "single")) {
  linkage <- match.arg(linkage)
  n <- nrow(sub_df)
  if (!n) return(integer(0))
  thr <- suppressWarnings(as.numeric(thr))
  if (!length(thr) || !is.finite(thr) || thr <= 0) thr <- 1e-9

  cols <- deblink_coord_cols(sub_df)
  if (length(cols) < 2L) return(rep(1L, n))
  tid <- as.character(sub_df$tid)
  P   <- apply_deblink_weights(as.matrix(sub_df[, cols, drop = FALSE]), w)

  g <- split(seq_len(n), tid)
  M <- t(vapply(g, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]),
                numeric(ncol(P))))

  # A trace whose centroid is non-finite (e.g. its z is entirely NA) makes
  # dbscan throw. That used to hit the blanket tryCatch below and return
  # "every trace is its own molecule" for the WHOLE dataset - de-blinking
  # silently disabled, with an inflated molecule count and no warning. Instead,
  # cluster the traces that do have a usable centroid and give each unusable one
  # its own id, which is the merge-only rule applied to a trace we cannot place.
  # An axis that is unmeasured across the WHOLE dataset (a z column that is all
  # NA rather than the synthesised 0) is dropped rather than allowed to
  # invalidate every centroid - de-blinking then proceeds in the axes that exist.
  keep_ax <- apply(is.finite(M), 2L, any)
  if (any(keep_ax) && !all(keep_ax)) M <- M[, keep_ax, drop = FALSE]

  ok <- is.finite(rowSums(M))
  cl <- integer(nrow(M))
  if (any(ok)) {
    cl[ok] <- if (identical(linkage, "complete")) {
      tryCatch(link_complete_threshold(M[ok, , drop = FALSE], thr),
               error = function(e) tryCatch(
                 dbscan::dbscan(M[ok, , drop = FALSE], eps = thr, minPts = 1L)$cluster,
                 error = function(e2) seq_len(sum(ok))))
    } else {
      # Single linkage at `thr`: DBSCAN with minPts = 1 is exactly that, and it
      # cannot emit noise.
      tryCatch(dbscan::dbscan(M[ok, , drop = FALSE], eps = thr, minPts = 1L)$cluster,
               error = function(e) seq_len(sum(ok)))
    }
  }
  if (any(!ok)) cl[!ok] <- max(0L, max(cl[ok], 0L)) + seq_len(sum(!ok))

  lab <- stats::setNames(as.integer(cl), names(g))
  unname(lab[tid])
}

# =============================================================================
# suggest_centroid_threshold() - THE MERGE DISTANCE FROM THE TRACE SPREAD
# =============================================================================
#
# GOAL
#   Set the centroid merge distance from the per-axis trace spread the data
#   actually shows, and hand back the two distance distributions needed to see
#   whether that choice is sane.
#
# THE SCALE, AND WHERE 1.5 COMES FROM
#   Two re-blinks of ONE fluorophore have centroids separated by roughly the
#   distance the label wanders between blinks. Measured on ground-truth
#   simulation, that separation is median 0.48*sigma and 95th percentile
#   1.09*sigma, while the nearest DIFFERENT molecule's centroid sits ~8*sigma
#   away. A threshold of ~1.5*sigma therefore collects essentially every
#   re-blink while still sitting far below the distinct-molecule population, and
#   the simulated error count bottoms out there (10 errors at k = 1.54, 11 at
#   the k = 2.0 default this function ships with).
#
#   THE DEFAULT IS k = 2.8, MEASURED AGAINST GROUND TRUTH
#   Ground truth was built on a real file's own point cloud: real per-trace
#   residuals, real molecule positions, half the molecules given a second blink.
#   Scoring which trace PAIRS were correctly grouped, over a range of injected
#   between-blink displacement:
#
#     k        1 nm     2 nm     4 nm      mean
#     2.00    0.928    0.916    0.864     0.903   <- the old default
#     2.355   0.947    0.944    0.897     0.929   <- sigma -> FWHM, still tight
#     2.80    0.954    0.952    0.931     0.946
#     3.00    0.954    0.953    0.942     0.950
#
#   The optimum sits at 2.7-3.1 and the curve is flat across it, so 2.8 is taken
#   as the default rather than the exact peak. Complete linkage bounds a
#   molecule's width by the threshold, so overshooting is far cheaper than
#   undershooting - erring high is the safe direction.
#
#   FWHM (2.355 sigma) was tested as a principled-looking alternative and beats
#   the old 2.0, but is consistently too tight. It also does not describe the
#   right quantity: 2.355 is the width of the single-localization density
#   profile, whereas merging asks whether two CENTROIDS agree.
#
#   Normalising by each trace's own sigma, or by the centroid standard error
#   sigma/sqrt(n), was tested too. Per-trace sigma shrunk toward the global
#   median wins by ~0.004 F1 - real but negligible against the ~0.05 that
#   choosing k buys. sigma/sqrt(n) is markedly WORSE, because it tightens the
#   threshold for long traces while the physical displacement between blinks
#   does not shrink with n.
#
#   The old default assumed the label moves appreciably between blinks. Measured
#   on fiducial-locked SMLM (35 and 105 min runs) there is no coherent sample
#   drift and no dependence of pair separation on the time gap between blinks,
#   so the point-cloud spread - not motion - sets the scale. See
#   sample_drift_estimate(), which flags the rare file where that fails.
#
# ANISOTROPY
#   When per-axis weights are supplied the coordinates are whitened first, so
#   every axis contributes at the same precision and sigma collapses to the
#   geometric mean of the per-axis spreads. That is why the returned `sigma` can
#   differ from a plain RMS of sigma_x/y/z.
#
# ARGUMENTS
#   df    : filtered data.frame (tid, x, y [, z])
#   min_n : minimum localizations for a trace to contribute a spread estimate
#   use_z : include z when the data is genuinely 3D
#   w     : optional per-axis de-blink weights (same metric as the clustering)
#   k     : multiplier on sigma (default 2.0, see above)
#
# RETURNS
#   list(sigma_ax, sigma, suggested, within, between, frac_merged, n_traces,
#        dims) - `within` and `between` are the raw distances the plot draws.
#   NULL when the spread cannot be estimated.
# =============================================================================
suggest_centroid_threshold <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT,
                                       use_z = TRUE, w = NULL, k = 2.8) {
  if (!is.data.frame(df) || !nrow(df) || !all(c("tid","x","y") %in% names(df))) return(NULL)
  has_z <- isTRUE(use_z) && has_axial_data(df)
  cols  <- if (has_z) c("x","y","z") else c("x","y")

  keep <- stats::complete.cases(df[, cols, drop = FALSE])
  d    <- df[keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  P   <- as.matrix(d[, cols, drop = FALSE])
  if (!is.null(w)) P <- sweep(P, 2L, w[seq_len(ncol(P))], "*")
  tid <- as.character(d$tid)

  g <- split(seq_len(nrow(P)), tid)
  g <- g[lengths(g) >= max(2L, as.integer(min_n))]
  if (length(g) < 3L) return(NULL)

  # per-axis spread in the metric that will actually be clustered in
  sd_ax <- vapply(seq_len(ncol(P)), function(k2) {
    v <- vapply(g, function(ix) stats::sd(P[ix, k2]), numeric(1))
    stats::median(v[is.finite(v)])
  }, numeric(1))
  if (!all(is.finite(sd_ax)) || any(sd_ax <= 0)) return(NULL)
  sigma <- sqrt(mean(sd_ax^2))

  # within-trace: how far a localization sits from its own trace centroid
  within <- unlist(lapply(g, function(ix) {
    c0 <- matrixStats_colMedian(P[ix, , drop = FALSE])
    sqrt(rowSums(sweep(P[ix, , drop = FALSE], 2L, c0, "-")^2))
  }), use.names = FALSE)

  # between-trace: nearest-neighbour distance between trace centroids. The near
  # population is re-blinks of one fluorophore, the far population is distinct
  # molecules - the threshold has to land between them.
  M  <- t(vapply(g, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]), numeric(ncol(P))))
  nn <- tryCatch(dbscan::kNN(M, k = 1L)$dist[, 1L], error = function(e) numeric(0))
  nn <- nn[is.finite(nn)]

  sug <- k * sigma
  list(sigma_ax   = sd_ax,
       sigma      = sigma,
       suggested  = round(sug, 1),
       k          = k,
       within     = within,
       between    = nn,
       # what fraction of traces would find a partner at this threshold
       frac_merged = if (length(nn)) mean(nn <= sug) else NA_real_,
       n_traces   = length(g),
       dims       = ncol(P))
}

# =============================================================================
# deblink_threshold_for_channel() - ONE MERGE DISTANCE PER CHANNEL
# =============================================================================
#
# GOAL
#   Return a named list, channel id -> centroid merge distance in nm, so that
#   de-blinking uses a threshold matched to the localization precision of the
#   dye that produced those traces.
#
# WHY NOT ONE SHARED DISTANCE
#   The merge distance is k x the measured trace spread, and trace spread IS the
#   localization precision. Two fluorophores rarely share a photon yield, so they
#   rarely share a precision. A single pooled sigma sits between the two: it
#   over-merges the broader channel and under-merges the tighter one, and the
#   error is invisible because both channels report the same number.
#
#   `k` stays SHARED. It is a modelling choice about how many spreads apart two
#   blinks can be and still be one molecule - a statement about blinking, not
#   about the dye. Only sigma is per channel.
#
# THE FALLBACK IS EXPLICIT, NEVER SILENT
#   A channel with too few long traces to measure a spread (`suggest_centroid_
#   threshold` returns NULL) falls back to the shared `eps_centroid`. The
#   returned attribute `measured` records which channels were measured and which
#   fell back, so the UI can say so rather than showing a number that looks
#   measured but is not.
#
# ARGUMENTS
#   df       : data.frame carrying `channel` (from resolve_channels)
#   settings : smlm_settings or Shiny input. Reads eps_centroid, centroid_k,
#              sd_min_locs, deblink_anisotropic, deblink_shared_spread.
#
# RETURNS
#   named list (names = as.character(channel)) of distances in nm, with
#   attr "measured" (named logical) and attr "sigma" (named numeric).
# =============================================================================
deblink_threshold_for_channel <- function(df, settings) {
  s <- as_smlm_settings(settings)
  shared <- s$eps_centroid
  if (!is.finite(shared) || shared <= 0) shared <- SMLM_DEFAULTS$eps_centroid

  ids <- if ("channel" %in% names(df)) sort(unique(df$channel[!is.na(df$channel)])) else 1L
  out <- stats::setNames(as.list(rep(shared, length(ids))), as.character(ids))
  meas <- stats::setNames(rep(FALSE, length(ids)), as.character(ids))
  sig  <- stats::setNames(rep(NA_real_, length(ids)), as.character(ids))

  # One channel, or the user has asked for a single combined spread: nothing to
  # split, so the shared value stands.
  if (length(ids) < 2L || isTRUE(s$deblink_shared_spread))
    return(structure(out, measured = meas, sigma = sig))

  kk <- s$centroid_k
  if (!is.finite(kk) || kk <= 0) kk <- SMLM_DEFAULTS$centroid_k
  min_n <- s$sd_min_locs

  for (ch in ids) {
    sub <- df[which(df$channel == ch), , drop = FALSE]
    w <- if (isTRUE(s$deblink_anisotropic))
      tryCatch(deblink_axis_weights(sub, min_n = min_n)$w, error = function(e) NULL)
      else NULL
    r <- tryCatch(suggest_centroid_threshold(sub, min_n = min_n, use_z = TRUE,
                                            w = w, k = kk),
                  error = function(e) NULL)
    if (!is.null(r) && is.finite(r$sigma) && r$sigma > 0) {
      out[[as.character(ch)]]  <- kk * r$sigma
      meas[[as.character(ch)]] <- TRUE
      sig[[as.character(ch)]]  <- r$sigma
    }
  }
  structure(out, measured = meas, sigma = sig)
}

# =============================================================================
# sample_drift_estimate() - HOW FAR DID THE SAMPLE MOVE DURING THE RUN?
# =============================================================================
#
# GOAL
#   Put a number, with an interval, on how much the sample moved across the
#   acquisition - so the de-blink radius can be widened on the rare file that
#   drifted, and left on the trace spread everywhere else.
#
# WHY NOT WITHIN-TRACE SLOPES (tried, and it cannot work)
#   The obvious test - regress each axis on time inside a trace and look for a
#   common direction - has no power here. A trace lasts ~0.15 s, so drift of
#   even 0.3 nm/s moves the sample 0.05 nm while it is being recorded, against a
#   localization sigma of ~7 nm. Measured on synthetic data with drift injected
#   from 0 to 0.3 nm/s, the coherence statistic never left its null value.
#   Within-trace slopes see fast jitter only, never slow drift, and scaling such
#   a slope up to the length of the run produced absurd numbers (31 um).
#
# THE LEVER ARM IS THE GAP BETWEEN BLINKS
#   Take pairs of traces whose centroids are close enough to be re-blinks of one
#   molecule. If the sample drifted at velocity v, the later blink of a pair sits
#   at displacement v * dt from the earlier one. Regressing the SIGNED
#   displacement on the SIGNED time gap, through the origin, estimates v with a
#   lever arm of the whole acquisition instead of 0.15 s.
#
#   Pairs that are actually two DIFFERENT molecules do not bias this. Their
#   displacement is isotropic with zero mean, so they widen the interval but
#   cannot tilt the slope. That is why the radius can be generous.
#
# THE ONE REAL BIAS
#   Selecting pairs by separation truncates exactly the large displacements fast
#   drift produces, pulling the estimate toward zero. The radius is therefore set
#   well above the merge threshold, and drift large enough to push same-molecule
#   pairs outside it is UNDERSTATED rather than missed - on a badly drifting file
#   this is a lower bound.
#
# ARGUMENTS
#   df     : localization table (tid, tim, x, y[, z])
#   min_n  : traces shorter than this are ignored
#   radius : pair-selection radius in nm; defaults to 10 x the trace spread
#   conf   : confidence level for the interval on the drift
#
# RETURNS
#   list(rate, rate_lo, rate_hi, drift_run, drift_run_lo, drift_run_hi, dt_typ,
#   drift_typ, jitter, frac_directed, direction, n_units, span, n_pairs, radius,
#   dims, significant) - rates in nm/s, drift_run in nm across the acquisition,
#   drift_typ in nm across a typical gap between blinks, jitter in nm - or NULL
#   when it cannot be measured.
#
#   DRIFT vs JITTER: drift is the displacement that lines up with time and
#   accumulates; jitter is the scatter left once it is removed, which does not.
#   frac_directed says which of the two dominates (0 = pure jitter, 1 = a
#   coherent march). `group` supplies molecule ids (de-blinked IDs) so that a
#   fluorophore seen at several times gives a trajectory instead of one step.
# =============================================================================
sample_drift_estimate <- function(df, min_n = 5L, radius = NULL, conf = 0.95,
                                  group = NULL) {
  if (!is.data.frame(df) || !nrow(df) || !all(c("tid", "tim", "x", "y") %in% names(df)))
    return(NULL)
  cols <- intersect(c("x", "y", "z"), names(df))
  keep <- stats::complete.cases(df[, c("tim", cols), drop = FALSE])
  d <- df[keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  # A synthesised constant z would add a zero-variance axis to every distance.
  cols <- cols[vapply(cols, function(a) stats::sd(d[[a]], na.rm = TRUE) > 0, logical(1))]
  if (length(cols) < 2L) return(NULL)

  gid_row <- if (!is.null(group)) as.character(group)[keep] else NULL

  ord <- split(seq_len(nrow(d)), as.character(d$tid))
  ord <- ord[lengths(ord) >= max(3L, as.integer(min_n))]
  if (length(ord) < 20L) return(NULL)

  M  <- t(vapply(ord, function(ix) matrixStats_colMedian(as.matrix(d[ix, cols, drop = FALSE])),
                 numeric(length(cols))))
  tm <- vapply(ord, function(ix) stats::median(d$tim[ix], na.rm = TRUE), numeric(1))
  # One molecule id per TRACE. Traces never straddle molecules, so the first row
  # of a trace names its molecule.
  tg <- if (!is.null(gid_row)) vapply(ord, function(ix) gid_row[ix[1]], character(1)) else NULL
  ok <- is.finite(rowSums(M)) & is.finite(tm) & (if (is.null(tg)) TRUE else !is.na(tg))
  M  <- M[ok, , drop = FALSE]; tm <- tm[ok]
  if (!is.null(tg)) tg <- tg[ok]
  if (nrow(M) < 20L) return(NULL)

  if (is.null(radius) || !is.finite(radius) || radius <= 0) {
    sg <- stats::median(vapply(ord[ok], function(ix)
      sqrt(mean(apply(as.matrix(d[ix, cols, drop = FALSE]), 2L, stats::var), na.rm = TRUE)),
      numeric(1)), na.rm = TRUE)
    radius <- if (is.finite(sg) && sg > 0) 10 * sg else 60
  }

  # Pair traces that belong to the same molecule. Two routes:
  #   GROUPED - an explicit molecule id (the de-blinked ID). Highest purity, and
  #             a fluorophore seen at three separate times contributes three
  #             pairs, so one long-lived molecule gives a whole mini-trajectory
  #             rather than a single displacement.
  #   SPATIAL - no grouping supplied: pair traces sitting within `radius`.
  # The grouped route is TRUNCATED by whatever threshold produced the grouping:
  # a molecule whose blinks drifted further than that radius was never merged,
  # so it cannot appear here. That is exactly why this is a QC readout and must
  # not feed back into the merge threshold - it would be measuring its own input.
  n_units <- NA_integer_
  if (!is.null(tg)) {
    gs <- split(seq_along(tg), tg)
    gs <- gs[lengths(gs) >= 2L]
    if (!length(gs)) return(NULL)
    cb <- lapply(gs, function(v) utils::combn(v, 2L))
    ii <- unlist(lapply(cb, function(m) m[1, ]), use.names = FALSE)
    jj <- unlist(lapply(cb, function(m) m[2, ]), use.names = FALSE)
    n_units <- length(gs)
  } else {
    nb <- tryCatch(dbscan::frNN(M, eps = radius), error = function(e) NULL)
    if (is.null(nb)) return(NULL)
    ii <- rep(seq_along(nb$id), lengths(nb$id))
    jj <- unlist(nb$id, use.names = FALSE)
    sel <- jj > ii                                 # keep each pair once
    ii <- ii[sel]; jj <- jj[sel]
  }
  if (length(ii) < 20L) return(NULL)

  # Orient every pair the same way in time, so the displacement reads "later
  # minus earlier" and a common drift adds instead of cancelling.
  late  <- ifelse(tm[jj] >= tm[ii], jj, ii)
  early <- ifelse(tm[jj] >= tm[ii], ii, jj)
  dt <- tm[late] - tm[early]
  dP <- M[late, , drop = FALSE] - M[early, , drop = FALSE]
  use <- is.finite(dt) & dt > 0 & is.finite(rowSums(dP))
  dt <- dt[use]; dP <- dP[use, , drop = FALSE]
  if (length(dt) < 20L || sum(dt^2) <= 0) return(NULL)

  # Regression through the origin, per axis: v = sum(dt * dP) / sum(dt^2).
  sdt2  <- sum(dt^2)
  v     <- as.numeric(crossprod(dt, dP)) / sdt2
  resid <- dP - outer(dt, v)
  dfree <- max(length(dt) - 1L, 1L)
  se    <- sqrt(colSums(resid^2) / dfree / sdt2)
  tcrit <- stats::qt(1 - (1 - conf) / 2, dfree)

  span <- diff(range(d$tim, na.rm = TRUE))
  rate <- sqrt(sum(v^2))
  # The displacement that matters is NOT the drift across the whole run - it is
  # the drift across the gap between two blinks of one molecule, which is what
  # the merge radius has to cover. The median gap among candidate pairs is the
  # available proxy for it, and it is typically a fraction of the run.
  dt_typ <- stats::median(dt)
  # Interval on the SPEED: propagate the per-axis errors along the fitted
  # direction, which is the quantity the caller acts on.
  se_r <- if (rate > 0) sqrt(sum((v / rate)^2 * se^2)) else sqrt(sum(se^2))
  lo <- max(0, rate - tcrit * se_r); hi <- rate + tcrit * se_r

  # DRIFT vs JITTER. Drift is the part of the displacement that lines up with
  # time - one direction, growing with the gap. Jitter is what is left once that
  # is removed: how far two blinks of one molecule land apart for reasons that do
  # NOT accumulate. Both are reported because they mean different things about a
  # sample, and because jitter is the term that actually sets how wide a de-blink
  # radius has to be.
  #   jitter here still contains the localization error of the two centroids, so
  #   it is an upper bound on any real positional wobble.
  jitter <- sqrt(mean(rowSums(resid^2)))

  # Displacements FROM EACH MOLECULE'S FIRST BLINK, for the drift plot: every
  # molecule starts at the origin and each later blink is one point. This is a
  # subset of the pairs used for the fit (the fit also uses second-to-third
  # steps, which carry information but have no common origin to plot against),
  # so the plot shows the same phenomenon on a readable frame rather than the
  # exact rows the regression saw.
  if (!is.null(tg)) {
    firsts <- vapply(split(seq_along(tg), tg), function(v) v[which.min(tm[v])], integer(1))
    is_first <- logical(length(tm)); is_first[firsts] <- TRUE
    ssel <- is_first[early[use]]
  } else {
    ssel <- rep(TRUE, length(dt))
  }
  steps <- if (any(ssel)) {
    st <- as.data.frame(dP[ssel, , drop = FALSE])
    names(st) <- paste0("d", cols)
    st$dt <- dt[ssel]
    # Which molecule each step belongs to, so the caller can colour by it and
    # join one molecule's blinks into a trajectory. NA when there is no grouping
    # (spatial pairing has no molecule to name).
    st$mol <- if (!is.null(tg)) tg[early[use]][ssel] else NA_character_
    # WHEN in the acquisition this blink happened, as distinct from dt (how long
    # AFTER its molecule's first blink it happened). Colouring by this shows
    # whether displacements point one way early in the run and another way
    # later - a reversal or a step that a gap-based colour cannot reveal.
    st$t_abs <- tm[late[use]][ssel]
    st[order(st$mol, st$dt), , drop = FALSE]
  } else NULL
  # Directional purity: how much of the total displacement the fitted drift
  # explains. Near 0 = scattered (jitter only), near 1 = a coherent march.
  tot <- sum(dP^2)
  frac_dir <- if (tot > 0) max(0, 1 - sum(resid^2) / tot) else NA_real_

  list(rate = rate, rate_lo = lo, rate_hi = hi,
       drift_run = rate * span, drift_run_lo = lo * span, drift_run_hi = hi * span,
       dt_typ = dt_typ, drift_typ = rate * dt_typ,
       jitter = jitter, frac_directed = frac_dir, n_units = n_units,
       direction = if (rate > 0) v / rate else rep(NA_real_, length(cols)),
       velocity = stats::setNames(v, cols), axes = cols, steps = steps,
       span = span, n_pairs = length(dt), radius = radius, dims = length(cols),
       significant = is.finite(lo) && lo > 0)
}

# =============================================================================
# deblink_merge_majority_vote() - LET THE LOCALIZATIONS VOTE
# =============================================================================
#
# GOAL
#   Replace "do the point clouds touch anywhere?" with "do MOST of this trace's
#   localizations think they belong to that other trace?".
#
# THE CRITERION
#   For every localization p in trace X, compare two distances:
#     d_own  |p - centroid of X, computed WITHOUT p|
#     d_alt  |p - centroid of the nearest OTHER trace Y|
#   p votes to merge X into Y when d_alt <= rho * d_own. X and Y merge when at
#   least a fraction `f` of X's localizations vote the same way.
#
# WHY LEAVE-ONE-OUT IS NOT OPTIONAL
#   The centroid of X is fitted TO X's own points, so d_own is systematically
#   too small - every point is unfairly close to its own centre. The effect goes
#   as ~1/n, so it hits the short traces hardest, which are exactly the ones
#   most in need of merging. Excluding p from its own centroid removes the bias
#   exactly, and is cheap for a mean: c^(-p) = (n*c - p) / (n - 1). That is why
#   this uses the MEAN centroid rather than the median used elsewhere.
#
# WHY THIS IS BETTER THAN THE RAW-POINT DBSCAN IT REPLACES
#   Single linkage merged when ANY point of one cloud came within eps of ANY
#   point of another, so one stray localization could fuse two molecules, and
#   the effective resolution degraded as traces got longer. Here a stray point
#   casts one vote and is outvoted.
#
#   It is also SCALE-FREE: no nanometre threshold. A diffuse trace automatically
#   demands proportionally more separation than a tight one, so heterogeneous
#   data (mixed labels, varying photon counts) needs no global sigma.
#
# SYMMETRY AND TRANSITIVITY
#   Votes are asymmetric - a short diffuse trace may vote for a long tight one
#   without the reverse. `mutual = TRUE` requires both directions, which is the
#   conservative reading. Accepted pairs then form a graph and connected
#   components become molecules; chaining is still possible but through a far
#   stricter edge test than single linkage.
#
# ARGUMENTS
#   sub_df : rows to cluster (tid, x, y, z)
#   f      : fraction of a trace's localizations that must vote to merge (0.5)
#   rho    : distance tolerance; 1.0 = the other centroid must be strictly
#            closer, >1 loosens
#   mutual : require both traces to vote for each other
#   w      : optional per-axis de-blink weights
#   k_cand : how many candidate centroids to consider per localization
#
# RETURNS  integer cluster id per ROW of sub_df, all > 0.
# =============================================================================
deblink_merge_majority_vote <- function(sub_df, f = 0.3, rho = 3.0, mutual = TRUE,
                                        w = NULL, k_cand = 4L) {
  n <- nrow(sub_df)
  if (!n) return(integer(0))
  f   <- min(max(as.numeric(f   %||% 0.3), 0.01), 1)
  rho <- max(as.numeric(rho %||% 3.0), 0.01)

  cols <- deblink_coord_cols(sub_df)
  if (length(cols) < 2L) return(rep(1L, n))
  tid  <- as.character(sub_df$tid)
  P    <- apply_deblink_weights(as.matrix(sub_df[, cols, drop = FALSE]), w)
  g    <- split(seq_len(n), tid)
  tn   <- names(g); nT <- length(g)
  if (nT < 2L) return(rep(1L, n))

  nt  <- lengths(g)
  tix <- match(tid, tn)                                   # trace index per row
  C   <- t(vapply(g, function(ix) colMeans(P[ix, , drop = FALSE]), numeric(ncol(P))))

  # leave-one-out own-centroid distance
  Cl    <- (C[tix, , drop = FALSE] * nt[tix] - P) / pmax(nt[tix] - 1L, 1L)
  d_own <- sqrt(rowSums((P - Cl)^2))

  # nearest OTHER trace centroid: ask for a few candidates and skip self
  kq <- min(as.integer(k_cand) + 1L, nT)
  kn <- tryCatch(dbscan::kNN(C, k = kq, query = P), error = function(e) NULL)
  if (is.null(kn)) return(as.integer(tix))
  best_id <- integer(n); best_d <- rep(Inf, n)
  for (j in seq_len(kq)) {
    cand <- kn$id[, j]; dd <- kn$dist[, j]
    take <- !is.finite(best_d) | (best_d == Inf)
    take <- take & cand != tix
    best_id[take] <- cand[take]; best_d[take] <- dd[take]
  }

  vote <- is.finite(best_d) & best_id > 0L & (best_d <= rho * pmax(d_own, 1e-9))
  if (!any(vote)) return(as.integer(tix))

  # Decide on the TOTAL vote, then wire up every sibling that shares it.
  #
  # Each point names only its single nearest other centroid, so a fluorophore
  # that blinked three times has its votes SPLIT across the two siblings and no
  # individual candidate clears f - the trace stays unmerged even though nearly
  # every point voted to merge with something. So the decision uses the total
  # vote fraction, and edges then go to each candidate holding a real share.
  vt <- data.frame(x = tix[vote], y = best_id[vote])
  tot <- tabulate(vt$x, nbins = nT) / nt                   # share of X's points voting
  keep_x <- which(tot >= f)
  if (!length(keep_x)) return(as.integer(tix))

  key   <- paste(vt$x, vt$y, sep = "_")
  cnt   <- table(key)
  parts <- do.call(rbind, strsplit(names(cnt), "_", fixed = TRUE))
  ed <- data.frame(x = as.integer(parts[, 1]), y = as.integer(parts[, 2]),
                   nv = as.integer(cnt))
  ed <- ed[ed$x %in% keep_x, , drop = FALSE]
  # a candidate counts if it holds a real share of the votes cast, not of all
  # points - otherwise siblings cancel each other out again
  ed$share <- ed$nv / pmax(tabulate(vt$x, nbins = nT)[ed$x], 1L)
  ed <- ed[ed$share >= 0.25, , drop = FALSE]

  if (isTRUE(mutual) && nrow(ed)) {
    fwd <- paste(ed$x, ed$y, sep = "_")
    rev <- paste(ed$y, ed$x, sep = "_")
    ed  <- ed[rev %in% fwd, , drop = FALSE]
  }
  # NOTE on f: the vote fraction for two traces of the SAME molecule is ~0.5 by
  # construction (leave-one-out makes d_own and d_alt identically distributed,
  # since both centroids estimate the same position). It falls to ~0.23 at 10 nm
  # and ~0.07 at 20 nm separation. So f must sit well BELOW 0.5 - thresholding
  # at 0.5 tests against the null value and rejects half of all true merges.

  if (!nrow(ed)) return(as.integer(tix))

  # connected components by union-find
  parent <- seq_len(nT)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  for (r in seq_len(nrow(ed))) {
    ra <- find(ed$x[r]); rb <- find(ed$y[r])
    if (ra != rb) parent[rb] <- ra
  }
  comp <- vapply(seq_len(nT), find, integer(1))
  as.integer(factor(comp))[tix]
}

# =============================================================================
# deblink_minpts() - minPts IS THE TRACE-LENGTH FILTER, NOT A FREE KNOB
# =============================================================================
# GOAL
#   Derive this diagnostic's minPts from the Trace ID Length filter's lower
#   bound instead of exposing it as a separate control.
#
# WHY IT IS NOT A FREE PARAMETER
#   minPts set independently of the trace-length filter silently contradicts it:
#   a trace long enough to pass filtering could still be called noise purely
#   because minPts was set higher, which is a length-dependent bias applied
#   outside the filter pipeline and invisible in the filter UI.
#
#   Tying them makes the rule self-consistent: any trace that survives the
#   length filter carries enough localizations to seed a molecule on its own,
#   so the only thing that can still leave it unmerged is its own spread -
#   which is exactly what the spread filter governs.
#
# ARGUMENTS  tid_length_range : the trace-length filter, c(lo, hi), or NULL
# RETURNS    integer >= 2 (DBSCAN needs at least 2)
# =============================================================================
deblink_minpts <- function(tid_length_range) {
  lo  <- suppressWarnings(as.numeric(tid_length_range[1]))
  if (!length(lo) || !is.finite(lo)) lo <- 4       # filter UI not rendered yet
  max(2L, as.integer(round(lo)))
}

# =============================================================================
# deblink_assign_singletons() - DE-BLINKING MERGES, IT NEVER REMOVES
# =============================================================================
# GOAL
#   Give every trace DBSCAN labelled noise (cluster 0) its own cluster ID, so
#   cluster 0 no longer exists after de-blinking.
#
# WHY
#   A fluorophore that blinked once, or whose localizations were too scattered
#   to link, is still a molecule. Dropping it here would delete molecules on a
#   spatial criterion applied outside the metadata and spread filters, and would
#   bias every count toward multi-blink species. Downstream code universally
#   filters `cluster != 0`, so leaving traces at 0 would remove them silently -
#   promoting them to singletons is what actually keeps them.
#
#   The cost is explicit and preferable: with a too-small radius, counts inflate
#   toward the trace count rather than quietly shrinking.
#
# ARGUMENTS  df : data.frame with `tid` and an integer `cluster` column
# RETURNS    df with cluster > 0 everywhere
# =============================================================================
deblink_assign_singletons <- function(df) {
  if (!nrow(df) || !"cluster" %in% names(df)) return(df)
  noise <- which(df$cluster == 0L | is.na(df$cluster))
  if (!length(noise)) return(df)

  nxt <- suppressWarnings(max(df$cluster, na.rm = TRUE))
  if (!is.finite(nxt) || nxt < 0) nxt <- 0L

  # One new ID per TRACE (not per point) - a noise trace stays one molecule.
  # Keyed on (tid, channel) for the same reason the plurality step is: under a
  # per-localization channel split, one trace can sit in two channels, and a
  # single shared ID would be a molecule spanning both.
  key <- if ("channel" %in% names(df))
    paste0(as.character(df$tid[noise]), "
", as.character(df$channel[noise]))
    else as.character(df$tid[noise])
  df$cluster[noise] <- as.integer(nxt) + as.integer(factor(key))
  df
}

# =============================================================================
# deblink_count_curve() - MOLECULE COUNT vs RADIUS, AND ITS KNEE   [EXPLORATORY]
# =============================================================================
#
# GOAL
#   Pick the de-blink radius from the data alone, with no model of precision and
#   no ground truth: sweep eps, run the REAL merge at each value, and plot how
#   many molecules come out.
#
# THE SHAPE OF THE CURVE AND WHY THE KNEE IS THE ANSWER
#   small eps  count is inflated - re-blinks of one fluorophore stay separate
#              and long traces fragment. The count falls STEEPLY as eps grows
#              and genuine re-blinks start to merge.
#   the knee   the point where that steep fall stops: the re-blinks have been
#              collected and nothing much is left to merge correctly.
#   large eps  a slow, steady decline - this is DISTINCT molecules being fused.
#              It never plateaus perfectly, because there is always another
#              neighbour to swallow.
#   So the knee is the largest radius that is still merging mostly re-blinks.
#
# WHY THIS NEEDS NO LABEL / TAG INFORMATION
#   Unlike a prior-based estimate (linker length, dye size), the knee is read off
#   the data's own behaviour. The label geometry is already baked into where the
#   knee falls, so asking the user to declare a tag would add an input that
#   changes nothing.
#
# WHY IT IS MARKED EXPLORATORY
#   The knee is a curvature maximum on a noisy, subsampled curve. It is stable
#   when the two regimes are well separated and unreliable when they are not -
#   which is the same crowding condition that makes de-blinking ambiguous in the
#   first place. Treat it as evidence, not as an answer, and look at the curve.
#
# ARGUMENTS
#   df       : filtered data.frame (x, y, z, tid)
#   minpts   : accepted for signature compatibility and currently UNUSED. The
#              sweep counts connected components at each radius (minPts = 1), so
#              there is no density criterion to match. De-blinking itself has no
#              minPts either: it merges on a distance threshold.
#   eps_seq  : radii to try (nm)
#   w        : optional per-axis de-blink weights, so the curve is measured in
#              the same metric the clustering will use
#   max_pts  : subsample cap - the sweep runs DBSCAN once per eps
#
# RETURNS  list(curve = data.frame(eps, n_molecules), knee, n_locs, n_traces)
# =============================================================================
deblink_count_curve <- function(df, minpts = 4L, eps_seq = seq(1, 40, by = 1),
                                w = NULL, max_pts = 60000L, method = "raw") {
  if (!is.data.frame(df) || !nrow(df) || !all(c("x", "y", "tid") %in% names(df)))
    return(NULL)
  cols <- intersect(c("x", "y", "z"), names(df))
  keep <- stats::complete.cases(df[, cols, drop = FALSE])
  d    <- df[keep, , drop = FALSE]
  if (!"z" %in% names(d)) d$z <- 0
  if (nrow(d) < 10L) return(NULL)
  if (nrow(d) > max_pts) {
    set.seed(1L)                                   # reproducible curve
    d <- d[sort(sample.int(nrow(d), max_pts)), , drop = FALSE]
  }
  P   <- apply_deblink_weights(as.matrix(d[, cols, drop = FALSE]), w)
  tid <- as.character(d$tid)
  n_tr <- length(unique(tid))

  # One eps -> the molecule count the pipeline would actually report: DBSCAN,
  # then TID-plurality, then unmerged traces become their own molecule.
  # Both criteria put every trace in a cluster, so the count is just the number
  # of distinct labels. For "vote" the swept quantity is the vote FRACTION, not
  # a distance - the caller supplies the right sequence.
  Mc <- NULL
  if (identical(method, "centroid")) {
    gg <- split(seq_len(nrow(P)), tid)
    Mc <- t(vapply(gg, function(ix) matrixStats_colMedian(P[ix, , drop = FALSE]),
                   numeric(ncol(P))))
  }
  # CENTROID ONLY. A vote sweep was tried and withdrawn: it did not return in a
  # usable time, and an unresponsive diagnostic is worse than an absent one.
  # The vote criterion has no distance to sweep anyway - its parameter is a
  # dimensionless fraction whose useful range is already documented on the input.
  if (is.null(Mc)) return(NULL)
  count_at <- function(e) {
    cl <- tryCatch(dbscan::dbscan(Mc, eps = e, minPts = 1L)$cluster,
                   error = function(err) NULL)
    if (is.null(cl)) NA_integer_ else length(unique(cl))
  }
  cnt <- vapply(eps_seq, count_at, numeric(1))
  ok  <- is.finite(cnt)
  if (sum(ok) < 5L) return(NULL)
  es <- eps_seq[ok]; cn <- cnt[ok]

  # Knee = maximum curvature of log(count) against eps. Logs because the fall is
  # multiplicative; curvature because the knee is where the SLOPE changes, and
  # the slope itself never reaches zero.
  knee <- NA_real_
  if (length(es) >= 5L) {
    lc <- log(pmax(cn, 1))
    d1 <- diff(lc) / diff(es)
    d2 <- diff(d1) / diff(es[-1])
    if (any(is.finite(d2))) knee <- es[which.max(replace(d2, !is.finite(d2), -Inf)) + 1L]
  }
  list(curve = data.frame(eps = es, n_molecules = cn),
       knee = knee, n_locs = nrow(d), n_traces = n_tr)
}

# =============================================================================
# apply_clustering_pipeline() - THE TWO-TIER CLUSTERING ENTRY POINT
# =============================================================================
#
# GOAL
#   Turn raw localizations into the two grouping levels the whole SMLM tab is
#   built on. This mirrors exactly what the interactive plot does
#   (`plottable_data()` in server_smlm.R), so batch tools and the on-screen plot
#   can never disagree - that duplication is the reason this function exists.
#
# THE TWO TIERS AND WHY THERE ARE TWO
#   DE-BLINKING (`cluster`).
#     Merges TRACE CENTROIDS - one point per trace, not the raw localizations -
#     when they sit within a threshold of a few times the MEASURED TRACE SPREAD
#     (~10-30 nm). A single fluorophore blinking on and off produces several
#     traces in the same place; this collapses them into one object. Without it,
#     one molecule is counted many times and every downstream count and spatial
#     statistic is inflated.
#
#     ⚠ THIS IS NOT DBSCAN, and it has not been for some time. The two criteria
#     are complete-linkage on the centroids (`deblink_merge_centroids`, the
#     default) and a majority vote of each trace's localizations
#     (`deblink_merge_majority_vote`). Neither has a minPts and neither has a
#     density criterion, so neither can emit noise: every trace comes out with
#     an id and `cluster` has no 0. The single-linkage option calls
#     `dbscan::dbscan(minPts = 1)`, but that is exactly single-linkage connected
#     components at the threshold, used as a components engine rather than as
#     density clustering. Older comments called this stage "primary DBSCAN";
#     they were wrong and have been corrected. The input VALUES "primary" and
#     "pcid" survive in saved settings and still mean the de-blinked level.
#
#   MEDIAN CLUSTERING (`cluster2`) - STRUCTURE.
#     Runs over the MEDIAN position of each unit from the level above, with a
#     LARGE eps (~100-300 nm). Its job is biological: group molecules into the
#     structure they form (a filament, a ring, a complex). THIS is where DBSCAN
#     actually runs, and it is the only density-clustering stage in the pipeline.
#
#   The two distances are doing completely different jobs - the de-blink
#   threshold is set by the MICROSCOPE (precision + label size), cluster2's eps
#   by the BIOLOGY (how far apart molecules in one structure sit). Tuning them
#   together is the most common source of confusion.
#
# WHY DBSCAN FOR THE MEDIAN LEVEL (and its alternatives here)
#   DBSCAN needs no pre-set cluster count, finds arbitrary shapes, and labels
#   sparse points as noise (cluster 0) rather than forcing them into a group -
#   all three matter for SMLM data. It does assume roughly uniform density, so
#   the median level also offers:
#     HDBSCAN ... adaptive density; only a minimum cluster size, no eps
#     OPTICS .... adaptive; keeps small structures a single eps would merge
#
# THE CHANNEL RULE (enforced, not optional) - BOTH TIERS
#   Both tiers run SEPARATELY PER CHANNEL, with IDs offset so they can never
#   collide. Channels come from resolve_channels(): a 2-colour DCR/fluo split, a
#   multiplex build, or a single channel when neither is active.
#
#   Molecules of different colours must never be de-blinked into one molecule,
#   and molecules imaged in different exchange rounds must never be merged into
#   one cluster - the second would invent a structure that was never
#   simultaneously present. These are one rule, so there is one loop.
#
#   There is no pooled option. Pooling across channels is the thing the channel
#   column exists to prevent; `cluster_group_scope` chooses between "within"
#   (every channel, separately) and "only:<k>" (one channel, for inspection).
#
#   PARAMETERS SPLIT WHERE THE PHYSICS SPLITS. The de-blink merge distance is per
#   channel, because it is k x the trace spread and spread is a property of the
#   dye (see deblink_threshold_for_channel). eps2 is shared, because it is set by
#   the biology - how far apart molecules in one structure sit - which no choice
#   of fluorophore changes.
#
# DE-BLINKING MERGES, IT NEVER REMOVES
#   De-blinking has no "drop unclustered" behaviour and no minPts at all: it
#   merges trace centroids on a distance threshold (complete linkage or majority
#   vote - not DBSCAN), and a trace nothing merges with becomes its own
#   single-blink molecule (see deblink_assign_singletons), so `cluster` is > 0
#   everywhere.
#   Removing molecules is the job of the metadata and spread filters, where it
#   is visible and reversible - not a side effect of a clustering radius.
#
# ARGUMENTS
#   df    : an ALREADY-FILTERED data.frame (x, y, z, tid [, mplex_channel])
#   settings : smlm_settings() for scripts, or the Shiny input (read through
#              as_smlm_settings, so batch tools use exactly the settings the user
#              sees). Uses the "channels", "deblink" and "cluster2" settings plus
#              sd_min_locs; see R/smlm_settings.R.
#
# RETURNS
#   the same df with `channel`/`channel_label`/`channel_color` (from
#   resolve_channels) and `cluster` and/or `cluster2` integer columns appended.
#   For `cluster2`, ID 0 still means "noise / unassigned" and consumers must
#   exclude it before counting. For `cluster`, 0 no longer occurs - the existing
#   `cluster != 0` filters downstream are harmless no-ops that now keep every
#   trace which survived filtering.
# =============================================================================
apply_clustering_pipeline <- function(df, settings) {
  if (!nrow(df)) return(df)
  s <- as_smlm_settings(settings)

  # One channel definition for the whole run. Everything below loops over it.
  df <- resolve_channels(df, s)
  if (!nrow(df)) return(df)

  chans <- channels_in_scope(df, s$cluster_group_scope)

  # Run `fn` once per channel, offsetting the returned IDs so they can never
  # collide across channels. Running-max, not a fixed constant: the 2-colour
  # path used +10000L, which breaks silently at 3 channels or 10k clusters.
  per_channel <- function(n_rows, subset_idx, fn) {
    out <- integer(n_rows)
    next_id <- 1L
    for (ch in chans) {
      idx <- subset_idx(ch)
      if (!length(idx)) next
      cid <- fn(idx)
      if (!length(cid)) next
      cid[cid > 0] <- cid[cid > 0] + (next_id - 1L)
      out[idx] <- cid
      if (any(cid > 0)) next_id <- max(next_id, max(cid, na.rm = TRUE) + 1L)
    }
    out
  }

  # ── De-blinking ─────────────────────────────────────────────────────
  if (isTRUE(s$enable_clustering)) {
    df$cluster <- 0L

    # Per-axis metric, computed ONCE over the whole df so every channel subset
    # below is clustered in the same space (a per-subset metric would make
    # cluster IDs incomparable across channels).
    aw1 <- if (isTRUE(s$deblink_anisotropic))
      deblink_axis_weights(df, min_n = s$sd_min_locs)$w
      else NULL

    # Two merge criteria, both trace-level. "centroid" compares where the traces
    # ARE; "vote" asks whether most of a trace's localizations think they belong
    # to a neighbouring trace. Raw-point DBSCAN was removed: one stray
    # localization could bridge two molecules. Ground-truth errors at each
    # method's own best threshold: raw-point 35, vote 29, centroid 10.
    # (Centroid scores 11 at the k = 2.0 default and 10 at its optimum 1.54 -
    #  quote the basis, not a bare number.)
    method1 <- s$deblink_method

    # THE MERGE DISTANCE IS PER CHANNEL (see deblink_threshold_for_channel).
    # It is dictated by the trace spread, and trace spread is a property of the
    # dye: two fluorophores with different photon yields have different
    # localization precision, so one shared distance under-serves the tighter
    # channel. The shared value remains available as an explicit fallback.
    thr_by_ch <- deblink_threshold_for_channel(df, s)

    run_primary <- function(idx) {
      sub <- df[idx, , drop = FALSE]
      if (identical(method1, "vote"))
        return(deblink_merge_majority_vote(sub,
                 f = s$vote_frac, rho = s$vote_rho,
                 mutual = !identical(s$vote_link, "either"), w = aw1))
      ch  <- sub$channel[1]
      thr <- thr_by_ch[[as.character(ch)]] %||% s$eps_centroid
      deblink_merge_centroids(sub, thr, aw1)
    }

    df$cluster <- per_channel(nrow(df),
                              function(ch) which(df$channel == ch),
                              run_primary)

    # TID-plurality reassignment, keyed on (tid, CHANNEL) rather than tid alone.
    #
    # THE KEY MATTERS, and getting it wrong silently undoes the channel loop.
    # A 2-colour split is assigned per TRACE, so there (tid, channel) and tid are
    # the same key and this changes nothing. A MULTIPLEX split is assigned per
    # LOCALIZATION, so a trace whose blinks straddle a TIM breakpoint belongs to
    # two channels at once. Keyed on tid alone, plurality then pulled all of that
    # trace into whichever channel held the majority - re-merging across the
    # boundary the loop above had just respected. Measured on a 3-window TIM split
    # of a real file: 2 de-blinked molecules spanned two channels, after the
    # per-channel clustering had produced none.
    #
    # Splitting the trace instead is the honest reading: the localizations really
    # were recorded in different exchange rounds, so they are not evidence of one
    # simultaneously-present molecule. It does mean a straddling trace is counted
    # twice - the boundary case is split so each channel counts it once.
    tid_plurality <- df %>%
      dplyr::filter(cluster != 0L) %>%
      dplyr::group_by(tid, channel, cluster) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(tid, channel) %>%
      dplyr::slice_max(n, n = 1L, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(tid, channel, cluster_new = cluster)
    df <- df %>%
      dplyr::left_join(tid_plurality, by = c("tid", "channel")) %>%
      dplyr::mutate(cluster = dplyr::coalesce(cluster_new, 0L)) %>%
      dplyr::select(-cluster_new)

    # Merge-only: unmerged traces become their own molecule rather than being
    # dropped. See deblink_assign_singletons().
    df <- deblink_assign_singletons(df)
  } else {
    df$cluster <- NA_integer_
  }

  # ── Secondary clustering (structure) ──────────────────────────────────────
  if (isTRUE(s$enable_cluster2)) {
    eps2    <- s$eps2
    minpts2 <- s$minpts2
    use_primary <- isTRUE(s$enable_clustering) && "cluster" %in% names(df)

    # eps2 is SHARED across channels: it is set by the biology (how far apart
    # molecules in one structure sit), not by the dye. The clustering still runs
    # strictly per channel, so a Median Cluster can never hold molecules from
    # two channels.
    if (use_primary) {
      medians <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0L) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          id = dplyr::first(cluster),
          mx = stats::median(x, na.rm = TRUE),
          my = stats::median(y, na.rm = TRUE),
          mz = stats::median(z, na.rm = TRUE),
          channel = dplyr::first(channel),
          analysis_group = {
            v <- analysis_group[!is.na(analysis_group)]
            if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_
          },
          .groups = "drop"
        )
    } else {
      medians <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          id = dplyr::first(tid),
          mx = stats::median(x, na.rm = TRUE),
          my = stats::median(y, na.rm = TRUE),
          mz = stats::median(z, na.rm = TRUE),
          channel = dplyr::first(channel),
          analysis_group = dplyr::first(analysis_group),
          .groups = "drop"
        )
    }
    medians <- medians[is.finite(medians$mx) & is.finite(medians$my) &
                       is.finite(medians$mz), , drop = FALSE]

    run_secondary <- function(idx) {
      if (length(idx) < minpts2) return(integer(length(idx)))
      # DIMENSIONALITY IS A PROPERTY OF THE DATA, NOT OF THE VIEW. Median
      # clustering always measures distances in 3D. Following the `mode3d`
      # DISPLAY toggle here meant flipping the view to 2D silently re-clustered
      # the data. A placeholder z (all zeros) contributes nothing to a Euclidean
      # distance, so always-3D is safe for 2D data; on real z it is simply
      # correct. Measured on real MINFLUX files (z sd 112-138 nm), 3D
      # nearest-neighbour distances are 1.46-1.89x the 2D projection, so the
      # choice moved eps2 by nearly a factor of two.
      pts <- as.matrix(medians[idx, c("mx", "my", "mz")])
      run_cluster2_grouping(
        pts,
        method        = s$cluster2_method,
        eps           = eps2,
        minPts        = minpts2,
        xi            = s$cluster2_xi,
        prefilter     = isTRUE(s$cluster2_prefilter),
        prefilter_r   = s$cluster2_prefilter_r,
        prefilter_min = s$cluster2_prefilter_min,
        prune         = isTRUE(s$cluster2_prune),
        prune_tol     = s$cluster2_prune_tol,
        core_only     = isTRUE(s$cluster2_core_only),
        min_members   = s$cluster2_min_members
      )
    }

    medians$cluster2 <- per_channel(nrow(medians),
                                    function(ch) which(medians$channel == ch),
                                    run_secondary)

    join_key <- if (use_primary) c("cluster" = "id") else c("tid" = "id")
    df <- df %>%
      dplyr::select(-dplyr::any_of("cluster2")) %>%
      dplyr::left_join(medians[, c("id", "cluster2")], by = join_key)
  } else {
    df$cluster2 <- NA_integer_
  }

  df
}

