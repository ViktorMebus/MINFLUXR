## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_hotspots.R  -  Local density and hotspots
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Local enhancement (k-NN density over the cell mean), the quadrat
##   dispersion index, the CSR-calibrated DBSCAN eps, the focus table and
##   the per-file hotspot summary.
##
## CONTENTS (in file order)
##   compute_local_enhancement, compute_dispersion_index,
##   calibrate_dbscan_eps, compute_foci_table, csr_cv_local_density,
##   hotspot_summary_row, hotspot_unit_label
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

## ═══════════════════════════════════════════════════════════════════════════════
## 8c. LOCAL DENSITY & HOTSPOTS
## ═══════════════════════════════════════════════════════════════════════════════
## WHY THIS BLOCK EXISTS
## Block 8b removes the labelling-density confound by THINNING every file down to a
## common density. That is the right move when you want to compare a curve. It is
## the wrong move when the thing you are looking for is a small number of very
## dense spots, because thinning throws away exactly the points that make them
## dense.
##
## This block takes the other route: leave the data alone and build statistics that
## are DIMENSIONLESS BY CONSTRUCTION, each unit measured against its own cell's
## mean density rather than against an absolute scale. A cell that expressed twice
## as much protein has twice the density everywhere, so a ratio to the cell mean is
## unchanged.
##
## THE FAILURE THIS BLOCK IS ANSWERING
## Mean and median nearest-neighbour distance are CENTRAL statistics. When a
## minority of molecules sit in foci and the rest are diffuse, the diffuse majority
## sets the median and the foci move it hardly at all. Measured on a real
## membrane-protein pair (a self-associating variant against a non-associating one,
## 10 acquisitions): median NN correlated with density at Spearman rho = -1.00
## within one condition and did not separate the conditions at all, while every
## statistic in this block separated them with non-overlapping ranges. NN was not
## noisy - it was measuring density.
##
## SO: read the TAIL of the local-density distribution, not its centre.
##
## WHAT IS HERE
##   compute_local_enhancement .. per-unit local density / that cell's mean density
##   compute_dispersion_index ... variance-to-mean of quadrat counts (1 under CSR)
##   calibrate_dbscan_eps ....... eps that a CSR pattern in THIS cell's own mask
##                                would barely cluster - makes eps comparable
##                                across cells of different density
##   compute_foci_table ......... one row per detected focus (n, Rg, centroid)
##   hotspot_summary_row ........ one row per file, everything above collapsed
##
## EVERY FUNCTION HERE TAKES A `support` (block 8b's footprint) AND USES IT. The
## cell mean density is n / support area, so with the bounding box in place of the
## footprint a cell imaged half in frame reads at half its true density and every
## enhancement doubles. The footprint is not optional for these statistics.

# =============================================================================
# compute_local_enhancement() - LOCAL DENSITY, IN UNITS OF THIS CELL'S MEAN
# =============================================================================
#
# GOAL
#   For every unit, estimate the density of units in its immediate neighbourhood
#   and express it as a MULTIPLE of the density the same cell would have if its
#   units were spread evenly over the footprint. 1.0 means "as dense as this cell
#   is on average"; 20 means "twenty times denser than this cell's own average".
#
# HOW THE LOCAL DENSITY IS ESTIMATED
#   The k-th nearest-neighbour distance r_k. In d dimensions, k units inside a
#   ball of radius r_k gives lambda_i = k / (omega_d * r_k^d), omega_2 = pi,
#   omega_3 = 4pi/3. This is the standard k-NN density estimator.
#
#   WHY k-NN AND NOT A FIXED-RADIUS COUNT: a fixed radius has to be chosen, and the
#   right choice depends on the density you are trying to measure - which is the
#   unknown. k-NN adapts its own scale, so one k works across a 7x density range.
#   WHY NOT A KDE: a KDE needs a bandwidth, same problem, and it smooths away the
#   isolated units that carry the "sparser background" half of the signal.
#
# CHOOSING k
#   k is a smoothing parameter, not a cluster size. Small k is noisy (k = 1 has a
#   ~52% coefficient of variation on a random pattern); large k blurs a focus into
#   its surroundings and cannot resolve a focus holding fewer than about 2k units.
#   k = 5 is the default: stable, and still resolves the ~17-unit foci that the
#   membrane-protein data set showed. Raise it for large diffuse domains, lower it
#   only if you are chasing pairs.
#
# THE BIAS THAT IS LEFT, AND WHICH WAY IT POINTS
#   The footprint estimator shrinks on a sparse cell, which raises lambda_cell and
#   therefore LOWERS enhancement. That acts AGAINST finding hotspots in sparse
#   cells, so a positive result on a sparse cell is conservative. A negative one on
#   a sparse cell is not evidence of absence.
#
# ARGUMENTS
#   pts     : data.frame with x, y (and z when use_3d) in nm, one row per unit
#   k       : neighbour order for the density estimate
#   use_3d  : measure in 3D (spheres) instead of 2D (discs)
#   support : an estimate_support() object; its $support field is the area/volume.
#             If NULL, falls back to the bounding box and flags it in used_bbox.
#
# RETURNS
#   list(enh, lambda_i, lambda_cell, r_k, k, d, n, area, used_bbox, keep)
#   `enh` follows the order of the input rows that had finite coordinates, which
#   `keep` identifies.
# =============================================================================
compute_local_enhancement <- function(pts, k = 5L, use_3d = FALSE, support = NULL) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  keep <- stats::complete.cases(m)
  m <- m[keep, , drop = FALSE]
  n <- nrow(m); d <- ncol(m)
  if (n < 3L) return(NULL)

  k_eff <- min(max(1L, as.integer(k)), n - 1L)
  nn <- tryCatch(RANN::nn2(m, m, k = k_eff + 1L), error = function(e) NULL)
  if (is.null(nn)) return(NULL)
  r_k <- nn$nn.dists[, k_eff + 1L]        # column 1 is the unit itself (d = 0)

  omega <- if (d == 3L) 4 * pi / 3 else pi
  # Guard the degenerate r_k = 0 (two units at identical coordinates, which
  # de-blinking can leave behind): it would divide by zero and produce Inf, which
  # then poisons every quantile. Floor it at half the smallest positive radius.
  pos <- r_k[is.finite(r_k) & r_k > 0]
  floor_r <- if (length(pos)) max(min(pos) / 2, 1e-6) else 1e-6
  r_k <- pmax(r_k, floor_r)
  lambda_i <- k_eff / (omega * r_k^d)

  # The cell's own mean density: units per unit area/volume of the FOOTPRINT.
  used_bbox <- FALSE
  area <- support$support %||% NA_real_
  if (!is.finite(area) || area <= 0) {
    used_bbox <- TRUE
    rng <- apply(m, 2, range)
    area <- prod(pmax(rng[2, ] - rng[1, ], 1))
  }
  lambda_cell <- n / area

  list(enh = lambda_i / lambda_cell, lambda_i = lambda_i, lambda_cell = lambda_cell,
       r_k = r_k, k = k_eff, d = d, n = n, area = area, used_bbox = used_bbox,
       keep = keep)
}

# =============================================================================
# compute_dispersion_index() - VARIANCE-TO-MEAN OF QUADRAT COUNTS
# =============================================================================
#
# GOAL
#   One dimensionless number per cell that is 1.0 under complete spatial
#   randomness BY CONSTRUCTION, so it needs no simulated null to be read.
#
# THE STATISTIC
#   Lay a grid of square (2D) or cubic (3D) quadrats of side `quadrat_nm` over the
#   window, count units in each, and take var(count)/mean(count). A Poisson process
#   has var = mean, so VMR = 1. Clustering pushes it above 1; a regular or inhibited
#   pattern below. This is the index of dispersion (Fisher, Thornton & Mackenzie
#   1922), and the reason to keep it alongside the k-NN work is precisely that it is
#   a DIFFERENT KIND of estimator - counts in fixed cells rather than distances to
#   neighbours - so agreement between the two is real evidence rather than one
#   statistic restated.
#
# THE PART THAT MUST BE RIGHT: THE EMPTY QUADRATS
#   Only quadrats inside the window may be counted, and the empty ones among them
#   MUST be included. Counting only quadrats that contain a unit throws away every
#   zero and drives VMR down towards 0; counting quadrats outside the cell invents
#   zeros and drives it up without limit. Both errors are silent. Hence the window
#   area comes from `support` when there is one, and `used_bbox` records when it
#   did not.
#
# CHOOSING THE QUADRAT SIZE
#   VMR is scale-dependent and that is a feature, not a nuisance - it is asking
#   "is there structure AT THIS SCALE". Set it near the size of the feature you are
#   looking for. On the membrane-protein data set 250 nm quadrats separated the
#   conditions with non-overlapping ranges (0.89-2.07 against 2.46-5.17) while
#   500 nm quadrats overlapped, because the foci are ~200-400 nm across and a
#   500 nm quadrat averages a focus together with its surrounding empty membrane.
#   ALWAYS report the quadrat size with the number.
#
# ARGUMENTS
#   pts        : data.frame with x, y (and z when use_3d), nm
#   quadrat_nm : quadrat side in nm
#   use_3d     : cubic quadrats instead of square
#   support    : estimate_support() object; only its area is used here
#
# RETURNS
#   list(vmr, mean_count, var_count, n_quadrats, quadrat_nm, used_bbox, edge_quadrats)
#   edge_quadrats = TRUE when the area-based fallback (partial edge quadrats) was used.
# =============================================================================
compute_dispersion_index <- function(pts, quadrat_nm = 250, use_3d = FALSE,
                                     support = NULL) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < 3L) return(NULL)
  q <- max(1, as.numeric(quadrat_nm))
  d <- ncol(m)

  # WHOLE QUADRATS ONLY (fixed 2026-09-17). With a support mask, count only the
  # quadrats that lie entirely inside the footprint, empty ones included. Edge
  # quadrats that are partly outside have a smaller expected count, which adds
  # variance that is not clustering; the excess grows with density and with the
  # footprint's edge-to-area ratio. Measured on CSR (tests/validation/
  # hotspot_csr_null.R): a 2 x 8 µm rod read VMR 2.56 at 5000 units and 0.79 at
  # 200, where whole quadrats read ~1. Needs >= 10 whole quadrats; otherwise the
  # area-based count below is used and `edge_quadrats` is TRUE.
  if (!is.null(support) && !is.null(support$mask) && is.finite(support$h %||% NA_real_)) {
    mk <- support$mask; dm <- dim(mk); h <- support$h; o_m <- support$origin
    if (length(dm) == d) {
      # quadrat lattice anchored at the mask origin; quadrat a spans [o + a q, o + (a+1) q)
      n_ax <- floor(dm * h / q)
      if (all(n_ax >= 1)) {
        qgrid <- as.matrix(expand.grid(lapply(n_ax, function(k) seq_len(k) - 1L)))
        whole <- vapply(seq_len(nrow(qgrid)), function(r) {
          rng_cells <- lapply(seq_len(d), function(k) {
            lo <- floor(qgrid[r, k] * q / h) + 1L
            hi <- min(dm[k], ceiling((qgrid[r, k] + 1) * q / h))
            lo:hi
          })
          all(do.call(`[`, c(list(mk), rng_cells)))           # every overlapping cell inside
        }, logical(1))
        if (sum(whole) >= 10L) {
          qi   <- floor(sweep(m, 2, o_m, "-") / q)               # quadrat index of each unit
          ukey <- apply(qi, 1, paste, collapse = ",")
          wkey <- apply(qgrid[whole, , drop = FALSE], 1, paste, collapse = ",")
          cnt  <- as.integer(table(factor(ukey[ukey %in% wkey], levels = wkey)))
          mu   <- mean(cnt)
          if (is.finite(mu) && mu > 0)
            return(list(vmr = stats::var(cnt) / mu, mean_count = mu, var_count = stats::var(cnt),
                        n_quadrats = length(cnt), quadrat_nm = q, used_bbox = FALSE,
                        edge_quadrats = FALSE))
        }
      }
    }
  }

  # Fallback without a usable mask: index every unit into a quadrat on a lattice
  # anchored at the data minimum and infer the empty quadrats from the area.
  o   <- apply(m, 2, min)
  idx <- floor(sweep(m, 2, o, "-") / q)
  key <- apply(idx, 1, paste, collapse = ",")
  occupied <- table(key)

  # How many quadrats does the window hold? Empty ones must be counted, so this is
  # an AREA divided by a quadrat area, not a count of occupied cells.
  used_bbox <- FALSE
  area <- support$support %||% NA_real_
  if (!is.finite(area) || area <= 0) {
    used_bbox <- TRUE
    rng  <- apply(m, 2, range)
    area <- prod(pmax(rng[2, ] - rng[1, ], q))
  }
  n_q <- max(length(occupied), round(area / q^d))

  cnt <- c(as.integer(occupied), rep(0L, max(0L, n_q - length(occupied))))
  mu  <- mean(cnt)
  if (!is.finite(mu) || mu <= 0) return(NULL)
  list(vmr = stats::var(cnt) / mu, mean_count = mu, var_count = stats::var(cnt),
       n_quadrats = length(cnt), quadrat_nm = q, used_bbox = used_bbox,
       edge_quadrats = TRUE)
}

# =============================================================================
# calibrate_dbscan_eps() - ONE eps THAT MEANS THE SAME THING IN EVERY CELL
# =============================================================================
#
# GOAL
#   Choose DBSCAN's eps per cell so that "clustered" means the same thing in a
#   sparse cell and a dense one.
#
# WHY A FIXED eps CANNOT WORK ACROSS CELLS
#   DBSCAN with (eps, minPts) is a DENSITY THRESHOLD: a unit is core when at least
#   minPts units lie within eps, i.e. when the local density exceeds
#   minPts / (omega_d * eps^d). That threshold is absolute, but the cells are not -
#   across a 7x range of labelling density the same eps is crossed trivially in the
#   dense cells and never in the sparse ones. The number of "clusters" then tracks
#   transfection load, which is the confound the whole Analysis tab exists to avoid.
#
# THE CALIBRATION
#   For a candidate eps, scatter n units UNIFORMLY AT RANDOM INSIDE THIS CELL'S OWN
#   FOOTPRINT, run the same DBSCAN, and record what fraction of that random pattern
#   is assigned to a cluster. Repeat `nsim` times and average. Take the LARGEST eps
#   whose random-clustering fraction is still at or below `target_frac`. That eps is
#   the most permissive setting at which chance alone would not produce clusters, so
#   a cluster found at it is a real density excess for THIS cell and the resulting
#   focus statistics are comparable between cells.
#
#   Randomising inside the FOOTPRINT rather than the bounding box is what makes this
#   correct for a cell that does not fill the frame; scattering into empty field
#   would spread the random units too thin and hand back an eps far too permissive.
#
# WHAT IT DOES NOT FIX
#   Cluster COUNT is still not comparable and neither is the fraction of units in a
#   cluster - on real data that fraction reached 31% in a control cell, because at
#   low density the calibrated eps grows large enough to sweep up diffuse
#   background. The statistic that survived was the SIZE OF THE LARGEST focus.
#
# ARGUMENTS
#   pts         : data.frame with x, y (and z when use_3d)
#   support     : estimate_support() object WITH $mask - the whole point is to
#                 randomise inside the cell. Without a mask it falls back to the
#                 bounding box, which weakens the calibration.
#   minPts      : DBSCAN minPts, held fixed across cells
#   eps_grid    : candidate eps values in nm
#   target_frac : the false-cluster rate to calibrate to (0.05 = 5%)
#   nsim        : CSR patterns per candidate. 15 is enough to place the crossing
#                 between two grid rungs; more only sharpens a number the grid
#                 spacing has already quantised.
#
# RETURNS
#   list(eps, curve = data.frame(eps, frac), target_frac, minPts, hit_floor,
#        grid_extended, null_max_focus_mean, null_max_focus_q95)
#   The null_max_focus_* entries are the largest focus found in random patterns in
#   this cell's footprint at the calibrated eps (at least 20 draws).
#   `grid_extended` is TRUE when the grid had to be extended below its floor.
#   `hit_floor` is TRUE when even the extended grid (down to ~1/12 of the floor)
#   clustered too much - the returned eps is then the smallest tried and the
#   calibration did NOT converge.
# =============================================================================
calibrate_dbscan_eps <- function(pts, support = NULL, minPts = 5L,
                                 eps_grid = c(40, 60, 80, 100, 125, 150, 175,
                                              200, 250, 300, 400, 500),
                                 target_frac = 0.05, nsim = 15L, use_3d = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m)
  if (n < max(3L, minPts)) return(NULL)
  eps_grid <- sort(unique(eps_grid[is.finite(eps_grid) & eps_grid > 0]))
  if (!length(eps_grid)) return(NULL)
  nsim <- max(3L, as.integer(nsim))

  # Randomiser: inside the footprint when there is a mask, else inside the bbox.
  rng <- apply(m, 2, range)
  draw <- function() {
    p <- if (!is.null(support) && !is.null(support$mask))
      tryCatch(sample_in_support(support, n, use_3d = use_3d), error = function(e) NULL)
    else NULL
    if (is.null(p)) {
      out <- vapply(seq_len(ncol(m)),
                    function(j) stats::runif(n, rng[1, j], rng[2, j]), numeric(n))
      return(matrix(out, nrow = n))
    }
    if (use_3d) cbind(p$x, p$y, p$z %||% rep(0, n)) else cbind(p$x, p$y)
  }

  frac <- vapply(eps_grid, function(e) {
    mean(vapply(seq_len(nsim), function(i) {
      cl <- tryCatch(dbscan::dbscan(draw(), eps = e, minPts = minPts)$cluster,
                     error = function(err) NULL)
      if (is.null(cl)) NA_real_ else mean(cl > 0L)
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))

  ok <- is.finite(frac) & frac <= target_frac
  # BELOW THE GRID FLOOR (fixed 2026-09-17). In a dense cell even the smallest grid
  # eps clusters random units above the target: a CSR rod with 5000 units failed
  # to calibrate in 20/20 runs and then put 31% of random units "in foci". Step
  # the eps down by 20% at a time (to ~1/12 of the floor, 11 more rungs) until the
  # random pattern stays under the target. Cells that calibrated before are
  # unchanged, because the extension only runs when the whole grid failed.
  grid_extended <- FALSE
  if (!any(ok)) {
    grid_extended <- TRUE
    for (e in eps_grid[1] * 0.8^(1:11)) {
      fe <- mean(vapply(seq_len(nsim), function(i) {
        cl <- tryCatch(dbscan::dbscan(draw(), eps = e, minPts = minPts)$cluster,
                       error = function(err) NULL)
        if (is.null(cl)) NA_real_ else mean(cl > 0L)
      }, numeric(1)), na.rm = TRUE)
      eps_grid <- c(e, eps_grid); frac <- c(fe, frac)
      if (is.finite(fe) && fe <= target_frac) break
    }
    ok <- is.finite(frac) & frac <= target_frac
  }
  hit_floor <- !any(ok)
  eps <- if (hit_floor) min(eps_grid) else max(eps_grid[ok])

  # NULL FOR THE HEADLINE METRIC (added 2026-09-17). The largest focus found in a
  # RANDOM pattern at the calibrated eps grows with density (CSR in a filled frame:
  # mean 2.7 units at 200 units per cell, 9.4 at 5000), so a raw max_focus_n
  # compared across cells of different density is confounded. Record this cell's
  # own random distribution of the largest focus, at the eps actually used, so the
  # summary row can express the observed focus relative to it.
  null_max <- vapply(seq_len(max(nsim, 20L)), function(i) {
    cl <- tryCatch(dbscan::dbscan(draw(), eps = eps, minPts = minPts)$cluster,
                   error = function(err) NULL)
    if (is.null(cl)) NA_real_ else if (any(cl > 0L)) max(tabulate(cl[cl > 0L])) else 0
  }, numeric(1))
  null_max <- null_max[is.finite(null_max)]

  list(eps = eps, curve = data.frame(eps = eps_grid, frac = frac),
       target_frac = target_frac, minPts = minPts, hit_floor = hit_floor,
       grid_extended = grid_extended,
       null_max_focus_mean = if (length(null_max)) mean(null_max) else NA_real_,
       null_max_focus_q95  = if (length(null_max)) unname(stats::quantile(null_max, 0.95)) else NA_real_)
}

# =============================================================================
# compute_foci_table() - ONE ROW PER DETECTED HOTSPOT
# =============================================================================
#
# GOAL
#   Run DBSCAN at a given eps and describe each cluster it finds, so the headline
#   number can be a count of UNITS IN THE BIGGEST FOCUS rather than a curve.
#
# WHY Rg IS BIAS-CORRECTED
#   The plain radius of gyration is biased low at small n by the factor (n-1)/n,
#   because the centroid is estimated from the same points. Uncorrected, that
#   manufactures a "bigger clusters are bigger" trend out of nothing (29% low at
#   n = 2, 11% at n = 5, 1% at n = 50). The Bessel correction is applied here, the
#   same as compute_cluster_rog(). A residual Jensen bias survives the square root
#   at very small n and has no closed form.
#
# WHAT Rg DOES NOT TELL YOU
#   It is a single scalar for a shape that need not be round. Real foci in the
#   membrane-protein data came out as two lobes and as elongated chains, and Rg
#   reports both as "about 150 nm". `aniso` (the ratio of the largest to smallest
#   principal standard deviation) is returned alongside for exactly that reason:
#   read it before treating a focus as a disc.
#
# ARGUMENTS
#   pts    : data.frame with x, y (and z when use_3d)
#   eps    : DBSCAN eps in nm, normally from calibrate_dbscan_eps()
#   minPts : DBSCAN minPts
#
# RETURNS
#   list(table, cluster, keep, eps). `table` is one row per focus (focus, n, rg,
#   aniso, x, y, z) ordered largest first; `cluster` is the per-unit assignment
#   (0 = not in any focus) for the rows `keep` marks.
# =============================================================================
compute_foci_table <- function(pts, eps, minPts = 5L, use_3d = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  keep <- stats::complete.cases(m)
  m <- m[keep, , drop = FALSE]
  if (nrow(m) < max(3L, minPts) || !is.finite(eps) || eps <= 0) return(NULL)

  cl <- tryCatch(dbscan::dbscan(m, eps = eps, minPts = minPts)$cluster,
                 error = function(e) NULL)
  if (is.null(cl)) return(NULL)

  empty <- data.frame(focus = integer(0), n = integer(0), rg = numeric(0),
                      aniso = numeric(0), x = numeric(0), y = numeric(0),
                      z = numeric(0))
  ids <- sort(unique(cl[cl > 0L]))
  if (!length(ids))
    return(list(table = empty, cluster = cl, keep = keep, eps = eps))

  rows <- lapply(ids, function(id) {
    P   <- m[cl == id, , drop = FALSE]
    ctr <- colMeans(P)
    Q   <- sweep(P, 2, ctr, "-")
    nn  <- nrow(Q)
    # Bessel-corrected mean squared distance from the centroid.
    rg <- sqrt(sum(Q^2) / nn * (if (nn > 1L) nn / (nn - 1L) else 1))
    # Principal spread ratio: above ~1.6 the focus is not a disc, so look at it.
    sdv <- tryCatch(stats::prcomp(Q, center = FALSE)$sdev, error = function(e) NA_real_)
    an  <- if (all(is.finite(sdv)) && min(sdv) > 0) max(sdv) / min(sdv) else NA_real_
    data.frame(focus = id, n = nn, rg = rg, aniso = an,
               x = ctr[1], y = ctr[2],
               z = if (length(ctr) >= 3L) ctr[3] else NA_real_)
  })
  tb <- do.call(rbind, rows)
  tb <- tb[order(-tb$n), , drop = FALSE]
  list(table = tb, cluster = cl, keep = keep, eps = eps)
}

# =============================================================================
# csr_cv_local_density() - WHAT THE CV READS WHEN NOTHING IS CLUSTERED
# =============================================================================
#
# GOAL
#   Give the CV of local density a zero point, so it can be read as a multiple of
#   random rather than as a bare number.
#
# THE DERIVATION
#   Under a Poisson process the k-th neighbour distance satisfies
#   lambda*omega_d*r_k^d ~ Gamma(k, 1), and the estimator is
#   lambda_i = k / (omega_d r_k^d), so lambda_i is proportional to 1/Gamma(k, 1).
#   For X ~ Gamma(k, 1): E[1/X] = 1/(k-1) and Var(1/X) = 1/((k-1)^2 (k-2)), so
#
#       CV = sd/mean = 1 / sqrt(k - 2)
#
#   It depends ONLY on k. Not on the density, not on the area, not on the number
#   of units, and not on the dimension. That is what makes the raw CV comparable
#   between cells in the first place.
#
#   k = 3 -> 1.000, k = 5 -> 0.577, k = 10 -> 0.354, k = 20 -> 0.236.
#
# UNDEFINED BELOW k = 3
#   The variance of 1/Gamma(k,1) exists only for k > 2, and the mean only for
#   k > 1. At k = 1 or 2 the sampling distribution has no finite spread and the
#   ratio is meaningless, so this returns NA and the ratio column goes NA with it.
#
# VALIDATED AGAINST SIMULATION, NOT ASSUMED
#   The formula ignores two things the real estimator has: neighbourhoods overlap,
#   so the per-unit estimates are spatially correlated, and units near the edge of
#   the footprint lose neighbours. Both could shift the baseline. Measured on the
#   10 two-condition footprints at k = 5, 60 randomisations each, drawing into each file's
#   own footprint at its own n: the simulated median CV came out 0.53-0.59 against
#   the predicted 0.577, and it did NOT track the unit count (Spearman -0.07 over
#   an 86-546 range). So the closed form is usable as-is and no simulation is
#   needed at run time.
# =============================================================================
csr_cv_local_density <- function(k) {
  k <- suppressWarnings(as.numeric(k))
  if (!length(k) || !is.finite(k) || k <= 2) return(NA_real_)
  1 / sqrt(k - 2)
}

# =============================================================================
# hotspot_summary_row() - EVERYTHING ABOVE, COLLAPSED TO ONE ROW PER FILE
# =============================================================================
#
# GOAL
#   The row that the comparison table and the metric-versus-density plot both read.
#
# WHICH COLUMNS TO TRUST, MEASURED
#   On the membrane-protein pair (4 control against 6 self-associating cells, whose
#   densities were NOT significantly different, rank-sum p = 0.92) the columns with
#   completely non-overlapping ranges, ranked by how significant condition stayed
#   after regressing out log density, were:
#     max_focus_n (p 0.0005) > cv_local (0.0007) > enh_q90 (0.0025) >
#     vmr (0.0030) > frac_above (0.0057).
#   `cv_ratio` is `cv_local` divided by a constant that depends only on k, so it
#   ranks identically; it exists to be readable, not to add information. On that
#   data cond01 ran 1.6-2.4x random and cond02 3.6-6.8x.
#   `nn_med` and `clark_evans` are carried in the same row deliberately, so the
#   table shows the two that FAILED next to the ones that worked rather than
#   quietly omitting them.
#
# ARGUMENTS
#   pts        : one row per unit (from unit_centroids())
#   support    : estimate_support() object with $mask
#   k          : k for the local density estimate
#   quadrat_nm : quadrat side for the dispersion index
#   enh_thresh : the enhancement multiple that `frac_above` counts
#   eps_cal    : a calibrate_dbscan_eps() result, or NULL to skip the focus columns
#   minPts     : DBSCAN minPts
#
# RETURNS
#   a one-row data.frame, or NULL when there are too few units.
# =============================================================================
hotspot_summary_row <- function(pts, support = NULL, k = 5L, quadrat_nm = 250,
                                enh_thresh = 5, eps_cal = NULL, minPts = 5L,
                                use_3d = FALSE, dataset = NA_character_) {
  le <- compute_local_enhancement(pts, k = k, use_3d = use_3d, support = support)
  if (is.null(le)) return(NULL)
  di <- compute_dispersion_index(pts, quadrat_nm = quadrat_nm, use_3d = use_3d,
                                 support = support)

  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m); d <- ncol(m)

  # 1-NN distance and the Clark-Evans ratio, carried as the negative controls.
  nn1 <- tryCatch(RANN::nn2(m, m, k = 2L)$nn.dists[, 2], error = function(e) NA_real_)
  lam <- le$lambda_cell
  # Expected mean 1-NN under CSR: 0.5 * lambda^(-1/2) in 2D,
  # gamma(4/3) * (3 / (4 pi lambda))^(1/3) in 3D.
  exp_nn <- if (d == 3L) gamma(4/3) * (3 / (4 * pi * lam))^(1/3) else 0.5 / sqrt(lam)
  ce <- if (is.finite(exp_nn) && exp_nn > 0) mean(nn1, na.rm = TRUE) / exp_nn else NA_real_

  ft <- if (!is.null(eps_cal))
    compute_foci_table(pts, eps = eps_cal$eps, minPts = minPts, use_3d = use_3d) else NULL
  tb <- ft$table
  # Volume in um^3 for 3D, area in um^2 for 2D - the scale factor differs.
  scale_unit <- if (d == 3L) 1e9 else 1e6

  data.frame(
    dataset      = as.character(dataset),
    n_units      = n,
    area_um2     = le$area / scale_unit,
    density      = lam * scale_unit,
    # --- the readouts that worked ---
    cv_local     = stats::sd(le$lambda_i) / mean(le$lambda_i),
    # The same number as a multiple of what a random pattern gives at this k, so
    # 1.0 reads as "random" and 3.0 as "three times as uneven as random".
    cv_ratio     = (stats::sd(le$lambda_i) / mean(le$lambda_i)) /
                     csr_cv_local_density(le$k),
    enh_q90      = unname(stats::quantile(le$enh, 0.90, na.rm = TRUE)),
    enh_q99      = unname(stats::quantile(le$enh, 0.99, na.rm = TRUE)),
    frac_above   = mean(le$enh > enh_thresh, na.rm = TRUE),
    vmr          = di$vmr %||% NA_real_,
    eps_nm       = eps_cal$eps %||% NA_real_,
    n_foci       = if (is.null(tb)) NA_integer_ else nrow(tb),
    max_focus_n  = if (is.null(tb) || !nrow(tb)) 0L else as.integer(tb$n[1]),
    # Largest focus relative to this cell's random 95th percentile (1 = what chance
    # alone reaches); comparable across densities where max_focus_n is not.
    max_focus_null95 = eps_cal$null_max_focus_q95 %||% NA_real_,
    max_focus_ratio  = if (is.null(tb)) NA_real_ else
                         (if (!nrow(tb)) 0 else tb$n[1]) / max(1, eps_cal$null_max_focus_q95 %||% NA_real_),
    max_focus_rg = if (is.null(tb) || !nrow(tb)) NA_real_ else tb$rg[1],
    max_focus_aniso = if (is.null(tb) || !nrow(tb)) NA_real_ else tb$aniso[1],
    frac_in_foci = if (is.null(ft)) NA_real_ else mean(ft$cluster > 0L),
    # --- the two that FAILED, kept visible on purpose ---
    nn_med       = stats::median(nn1, na.rm = TRUE),
    clark_evans  = ce,
    used_bbox    = isTRUE(le$used_bbox),
    stringsAsFactors = FALSE)
}


# =============================================================================
# PANEL HELPERS (moved from server_smlm.R; plain functions of their arguments)
# =============================================================================

## THE UNIT IS THE WHOLE MEASUREMENT, so it is named in words on every axis and
## title rather than left to the control panel.
##
## RAW LOCALIZATIONS ARE NOT OFFERED, and this is not a convenience choice.
## Localizations inside one trace sit a few nm apart, so the k-th neighbour
## distance collapses to the localization precision and the local density
## measures HOW LONG THAT FLUOROPHORE STAYED ON, not how the protein is
## arranged. Measured on the two-condition pair: on localizations the enhancement runs
## to 3,000-17,000x the cell mean and NEITHER the CV nor the 90th percentile
## separates the conditions (CV 1.46-1.74 against 1.23-2.54, overlapping), while
## the same statistics on trace-ID centroids are cleanly disjoint. Trace length
## is photophysics; the spatial arrangement of separate traces is the biology.
hotspot_unit_label <- function(u) switch(u %||% "tid",
  tid     = "trace ID centroids",
  cluster = "de-blinked molecules",
  "units")

