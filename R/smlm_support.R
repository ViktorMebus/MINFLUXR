## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_support.R  -  Footprint (support) and labelling-density normalisation
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The region the sample occupies (occupancy grid, KDE, hull, box), its
##   set covariance for edge correction and its debiasing, sampling inside it,
##   per-file density and unit-level p-thinning.
##
## CONTENTS (in file order)
##   unit_centroids, estimate_support, .array_shift, kde_hdr_support,
##   support_area, sample_in_support, in_support, attach_gamma_cache,
##   .gamma_at, bbox_support, mask_components, set_covariance_iso,
##   mask_outline_segments, mask_project_2d, support_closing_sweep,
##   plateau_close_nm, nn_local_intensity, unit_local_density,
##   compute_unit_density, thin_units_to_p, dilate_mask_cells,
##   debias_support_window
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

## ═══════════════════════════════════════════════════════════════════════════════
## 8b. LABELLING-DENSITY NORMALISATION  (opt-in)
## ═══════════════════════════════════════════════════════════════════════════════
## WHY THIS BLOCK EXISTS
## Transient transfection gives each cell an undefined plasmid load, so labelling
## density varies cell to cell and run to run even within one condition. Almost
## every clustering read-out is density-dependent in the direction that MIMICS
## aggregation:
##   * k-NN distance scales as rho^(-1/d) - under pure CSR in 2D the mean 1-NN
##     distance is 0.5 * rho^(-1/2), so a high expresser "looks clustered".
##   * DBSCAN with fixed (eps, minPts) IS a density threshold: rho_crit =
##     minPts / (pi eps^2) in 2D. Raise rho and more clusters cross it, and they
##     chain into larger ones.
## So two transfections cannot be compared on those statistics as they stand.
##
## THE FIX AND WHY IT IS SOUND
## Independent ("Poisson") p-thinning - keep each unit with probability p,
## independently - sends the intensity lambda -> p*lambda and the second-order
## product density rho2 -> p^2 * rho2. The pair-correlation function
## g = rho2 / lambda^2 is therefore EXACTLY invariant, as is Ripley's K. Real
## clustering survives thinning; a pure density artefact does not. Thinning every
## file down to one common density therefore makes even the fixed-threshold
## DBSCAN read-outs comparable, without erasing the signal being measured.
##
## WHAT THE THINNING TARGET MUST BE - THE ONE REAL TRAP
## The target must be a density that is NOT itself inflated by clustering:
##   * GLOBAL density = n_units / support  is rearrangement-invariant. Clustering
##     moves points around inside the support but does not change n or A, so this
##     is the CORRECT thinning target.
##   * LOCAL NN-based intensity (nn_local_intensity) is CONFOUNDED with
##     clustering: aggregates shorten NN distances, inflating lambda_hat, so a
##     genuinely aggregating cell would be thinned HARDER and part of the signal
##     removed. It is offered as a DIAGNOSTIC and as a support-free cross-check,
##     and is flagged in the UI - it is not the recommended target.
## The high quantiles of the NN distribution are dominated by the sparse
## background rather than the aggregates, which is why nn_local_intensity takes a
## quantile range rather than the median (see `probs`).
## ─────────────────────────────────────────────────────────────────────────────

## One point per de-blinked unit - the point pattern every density statistic and
## every spatial statistic in this app should be measured on. A trace ID is ONE
## BLINK, so trace-ID density confounds expression level with photophysics and
## laser power and is NOT comparable between runs; prefer unit = "cluster"
## (de-blinked molecules) whenever de-blinking is enabled.
##   unit: "tid"     - one median position per trace ID
##         "cluster" - one median position per de-blinked molecule (noise dropped)
## Returns a data.frame(x, y, z) of finite centroids, or NULL if unavailable.
unit_centroids <- function(df, unit = "tid", use_3d = FALSE) {
  if (is.null(df) || !nrow(df) || !all(c("x", "y") %in% names(df))) return(NULL)
  if (!"z" %in% names(df)) df$z <- 0
  if (unit == "cluster") {
    if (!"cluster" %in% names(df) || !any(!is.na(df$cluster) & df$cluster > 0L)) return(NULL)
    pts <- df %>%
      dplyr::filter(!is.na(cluster), cluster > 0L) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop")
  } else {
    if (!"tid" %in% names(df)) return(NULL)
    pts <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                       y = stats::median(y, na.rm = TRUE),
                       z = stats::median(z, na.rm = TRUE), .groups = "drop")
  }
  pts <- pts[is.finite(pts$x) & is.finite(pts$y) & (!use_3d | is.finite(pts$z)), , drop = FALSE]
  if (!nrow(pts)) return(NULL)
  as.data.frame(pts[, c("x", "y", "z"), drop = FALSE])
}

## GOAL: estimate the AREA (2D) or VOLUME (3D) actually occupied by the sample,
## so that n_units / support is a labelling density rather than an artefact of how
## much empty frame was scanned.
##
## WHY NOT THE BOUNDING BOX: MINFLUX scan fields routinely include the edge of a
## cell and interior dark regions (organelle exclusions, unlabelled compartments).
## The bounding box counts all of that as "sample", so a cell imaged half-in-frame
## reads at half its true density purely from framing. That is exactly the
## confound this whole block exists to remove.
##
## modes:
##   "grid"  (DEFAULT, edge- and dark-region-aware) - occupancy grid. Bin the unit
##           centroids into cells of side `h_nm`, then MORPHOLOGICALLY CLOSE the
##           occupied set by `close_nm` (dilate then erode) and count. The closing
##           is the knob that separates the two things that both look like holes:
##             - sparse SAMPLING inside the cell -> should be filled in (it is
##               sample), which closing does;
##             - a genuine DARK REGION or the outside of the cell -> should stay
##               excluded, which closing leaves alone as long as close_nm is
##               smaller than the dark region.
##           So close_nm should be a few times the typical unit spacing and well
##           below the size of the dark features you want excluded.
##   "hull"  - convex hull (geometry::convhulln). Follows the outer edge but
##           bridges concavities and dark regions, so it sits between bbox and grid.
##   "bbox"  - axis-aligned bounding box. The implicit choice everywhere else in
##           the app; kept for comparison and as the fallback.
##
## Returns list(support, mode, h, close_nm, n_occupied, bbox, frac_of_bbox, d).
## `support` is in nm^2 (2D) or nm^3 (3D). NULL only if there is nothing to measure.
##
## HONEST LIMITS: the grid estimate is h-dependent (it can only ever be a multiple
## of h^d) and, like any union-of-cells support, is biased upward at the boundary
## by roughly perimeter x h/2. Neither matters much for the intended use, which is
## a RATIO of densities between files measured with the SAME h and close_nm - the
## boundary bias largely cancels between comparable fields. It does matter if you
## quote an absolute molecules/um^2, so both the value and the settings are
## reported, and the fraction of the bounding box is shown so an implausible
## support is visible at a glance.
## `return_mask = TRUE` additionally returns the occupancy array itself (`mask`)
## plus `origin`/`pad`/`dims`/`win`, so the support can be DRAWN over the data -
## which is the only way to see whether the gap-fill radius did something sane.
##
## DARK REGIONS AND THE GAP FILL - the relationship worth knowing, because it means
## no extra control is needed. A morphological closing of radius r can only fill an
## interior hole narrower than about 2r, so the gap fill IS the dark-region control:
##     a dark region wider than ~2 * close_nm stays excluded.
## MEASURED on a 5 um field with a 1500 nm dark region and sparse sampling: at the
## auto plateau (200 nm) the dark region is correctly excluded and the support lands
## within -0.1% of truth. A separate "keep holes bigger than X" knob was built and
## REMOVED - re-excluding interior holes of the raw occupancy made the estimate
## WORSE (-9.0% at 400 nm, -2.0% at 800 nm vs -0.1% without it), because at the raw
## scale sparse-sampling voids are face-connected to genuine dark regions and get
## re-excluded with them. If dark regions are being swallowed, lower close_nm - the
## support map shows immediately whether they are.
estimate_support <- function(pts, use_3d = FALSE, mode = c("grid", "hull", "bbox", "kde"),
                             h_nm = 100, close_nm = 200, return_mask = FALSE,
                             bw_nm = 150, frac = 0.99) {
  mode <- match.arg(mode)
  if (is.null(pts) || !nrow(pts)) return(NULL)
  # KDE / highest-density region: a separate estimator, always 2D (the field is a
  # plan-view density; a 3D HDR would need a volumetric kernel and is not offered).
  if (mode == "kde") {
    k <- kde_hdr_support(pts, res_nm = h_nm, bw_nm = bw_nm, frac = frac,
                         return_mask = return_mask)
    if (is.null(k)) return(NULL)
    rng2 <- apply(as.matrix(pts[, c("x", "y"), drop = FALSE]), 2, range)
    bb2  <- prod(pmax(rng2[2, ] - rng2[1, ], h_nm))
    k$bbox <- bb2
    k$frac_of_bbox <- k$support / bb2
    k$d <- 2L
    # Normalise the returned shape to match the grid/hull/bbox branches, so every
    # caller (density row, map renderer, status readout) can read the same fields.
    k$h        <- k$res
    k$close_nm <- NA_real_
    return(k)
  }
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < 1L) return(NULL)
  d <- ncol(m)

  rng  <- apply(m, 2, range)
  ext  <- rng[2, ] - rng[1, ]
  # A degenerate axis (all points equal) would zero the product; floor it at h so
  # a genuinely flat slab still reports a finite support instead of 0.
  ext_f <- pmax(ext, h_nm)
  bbox  <- prod(ext_f)

  out <- list(mode = mode, h = h_nm, close_nm = close_nm, d = d,
              bbox = bbox, n_occupied = NA_integer_)

  if (mode == "bbox" || nrow(m) < (d + 2L)) {
    out$support <- bbox
    if (mode != "bbox") out$mode <- "bbox (too few units for the chosen mode)"
    out$frac_of_bbox <- 1
    return(out)
  }

  if (mode == "hull") {
    hv <- tryCatch({
      h <- geometry::convhulln(m, options = "FA")
      if (d == 2L) h$vol else h$vol   # convhulln$vol is area in 2D, volume in 3D
    }, error = function(e) NA_real_)
    out$support      <- if (is.finite(hv) && hv > 0) hv else bbox
    out$frac_of_bbox <- out$support / bbox
    return(out)
  }

  # ── occupancy grid + morphological closing ────────────────────────────────
  h <- max(1e-6, as.numeric(h_nm))
  # Integer voxel index per point, 1-based, per axis.
  idx <- sweep(m, 2, rng[1, ], "-") / h
  idx <- floor(idx) + 1L
  dims <- as.integer(pmax(1L, apply(idx, 2, max)))
  # Pad by the closing radius in voxels so a dilation near the border has room.
  rad  <- max(0L, as.integer(floor(close_nm / h)))
  pad  <- rad + 1L
  dimp <- dims + 2L * pad
  idxp <- matrix(as.integer(sweep(idx, 2, rep(pad, d), "+")), ncol = d)

  occ <- array(FALSE, dim = dimp)
  occ[idxp] <- TRUE                      # matrix-index assignment: one row = one cell
  n_raw <- sum(occ)

  if (rad > 0L) {
    # Separable Chebyshev (box) dilation then erosion = morphological closing with
    # a cubic structuring element. Separable => O(rad) passes per axis, not O(rad^d).
    box_max <- function(a, axis, r) {
      res <- a
      for (s in seq_len(r)) {
        res <- res | .array_shift(res, axis, 1L) | .array_shift(res, axis, -1L)
      }
      res
    }
    dil <- occ
    for (ax in seq_len(d)) dil <- box_max(dil, ax, rad)
    # Erosion of X = complement of the dilation of the complement.
    ero <- !dil
    for (ax in seq_len(d)) ero <- box_max(ero, ax, rad)
    closed <- !ero
  } else {
    closed <- occ
  }

  n_occ <- sum(closed)
  out$n_occupied   <- as.integer(n_occ)
  out$n_raw_cells  <- as.integer(n_raw)
  out$support      <- n_occ * h^d
  out$frac_of_bbox <- out$support / bbox
  if (isTRUE(return_mask)) {
    # The mask is PADDED by `pad` cells on every side so the dilation had room, so
    # padded cell index i on an axis starts at world coordinate
    # rng[1, axis] + (i - 1 - pad) * h. `origin` is that world coordinate of the
    # first padded cell, which is what a renderer needs to place the grid.
    out$mask   <- closed
    out$origin <- rng[1, ] - pad * h
    out$pad    <- pad
    out$dims   <- dimp
    out$win    <- rng
  }
  out
}

# Shift a logical array by one cell along `axis` (vacated cells become FALSE).
# Helper for the separable morphological closing in estimate_support; internal.
# Uses `TRUE` as a whole-axis index so the same code works for 2D and 3D.
.array_shift <- function(a, axis, by) {
  dm <- dim(a)
  n  <- dm[axis]
  if (by == 0L) return(a)
  out <- array(FALSE, dim = dm)
  if (n < 2L) return(out)
  src <- if (by > 0L) seq_len(n - 1L) else 2:n
  dst <- if (by > 0L) 2:n             else seq_len(n - 1L)
  i_s <- rep(list(TRUE), length(dm)); i_s[[axis]] <- src
  i_d <- rep(list(TRUE), length(dm)); i_d[[axis]] <- dst
  slice <- do.call("[", c(list(a), i_s, list(drop = FALSE)))
  do.call("[<-", c(list(out), i_d, list(value = slice)))
}

## GOAL: support from a SMOOTHED DENSITY FIELD instead of a binary occupancy grid -
## a kernel density estimate of the localizations, thresholded at the level that
## encloses a given FRACTION of the units (a highest-density region, Hyndman 1996).
##
## WHY THIS FORM AND NOT AN ABSOLUTE DENSITY LEVEL. Thresholding at a fixed
## density level is the obvious version and it is CIRCULAR: the support becomes a
## function of the density it is the denominator for, so a sparsely labelled cell
## falls below the level over most of its area and its measured density is dragged
## toward the threshold. Measured on two scenes with identical geometry and a known
## 4.0x labelling difference, the recovered ratio was **0.77 and 0.55** - not merely
## compressed but INVERTED - and the support of the sparse scene collapsed to zero
## by n = 1000. The HDR form asks instead "where do 99% of the units live", which is
## invariant to how many there are: same scenes, recovered ratio **3.90 and 3.90**
## (occupancy grid: 3.85 / 3.78).
##
## Measured support error vs known truth at n = 3000, bandwidth 150 nm:
##   plain square   grid -0%, KDE-HDR -0.6%
##   L-shape        grid -0%, KDE-HDR +2.3%
##   square + dark  grid -0%, KDE-HDR -0.6%
## Both degrade at very low n (n = 500: grid -15%, KDE-HDR -18%). KDE-HDR holds its
## support better under thinning (to 25% of units: -3.7% vs the grid's -6.2%; to 10%:
## -18.3% vs -26.3%), which matters because thinning is a core feature here.
##
## WHY IT IS THE DEFAULT - measured on CLUSTERED data (60 structures in a 5 um cell,
## 4.0x labelling difference), recovered density ratio:
##   KDE-HDR  bw 150 -> 4.03,  250 -> 3.95,  400 -> 4.03,  600 -> 4.02
##   grid     fill 200 -> 3.20, 400 -> 3.53, 800 -> 3.83, 1200 -> 3.83
## The ratio is nearly INDEPENDENT of bandwidth, and right, even where the support is
## only 30% of the cell - because the level is calibrated per file to enclose the same
## FRACTION of units, so the two supports agree with each other (7.46e6 vs 7.51e6 nm2)
## whether or not either equals the cell. The grid's threshold is an absolute count
## (">=1 unit per cell"), so a sparsely labelled file gets a systematically smaller
## mask (8.05e6 vs 6.45e6, 25% apart) and the ratio is dragged down.
## Consequence for use: bandwidth sets the ABSOLUTE area (over-smoothing spills past
## the cell - 107% of it at bw 600, 123% at 1400) but barely affects the RATIO, which
## is what a between-condition comparison depends on.
##
## `frac` is the fraction of UNITS enclosed; `bw_nm` MUST be held common across the
## samples being compared (it is the smoothing, and a per-file bandwidth would make
## each file's gaps close differently). Dark regions are excluded whenever they are
## wider than roughly 2 * 3 * bw_nm, since no kernel reaches the middle.
##
## Scale-invariance note: the level is a quantile of the SAME field it thresholds, so
## multiplying the field (e.g. more localizations per unit) scales the level equally
## and the mask does not move. That is what makes it thinning-invariant.
##
## The field is built by BINNING then convolving with a separable Gaussian, so cost
## is set by the grid, not by the number of localizations.
kde_hdr_support <- function(pts, res_nm = 50, bw_nm = 150, frac = 0.99,
                           return_mask = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  if (!all(c("x", "y") %in% names(pts))) return(NULL)
  p <- pts[stats::complete.cases(pts[, c("x", "y")]), , drop = FALSE]
  n <- nrow(p)
  if (n < 3L) return(NULL)
  res <- max(1e-6, as.numeric(res_nm)); bw <- max(res / 2, as.numeric(bw_nm))
  pad <- 3 * bw

  xr <- range(p$x); yr <- range(p$y)
  gx <- seq(xr[1] - pad, xr[2] + pad, by = res)
  gy <- seq(yr[1] - pad, yr[2] + pad, by = res)
  nx <- length(gx); ny <- length(gy)
  if (nx < 3L || ny < 3L) return(NULL)

  ix <- pmin(pmax(findInterval(p$x, gx), 1L), nx)
  iy <- pmin(pmax(findInterval(p$y, gy), 1L), ny)
  cnt <- matrix(0, nx, ny)
  tb  <- table(ix, iy)
  cnt[cbind(as.integer(rownames(tb)[row(tb)]), as.integer(colnames(tb)[col(tb)]))] <-
    as.vector(tb)

  # Separable Gaussian smoothing via banded kernel matrices (exact, and independent
  # of n). Kernels are normalised so the field is units per nm^2.
  kern <- function(m) {
    d <- outer(seq_len(m), seq_len(m), function(a, b) (a - b) * res)
    exp(-0.5 * (d / bw)^2)
  }
  z <- kern(nx) %*% cnt %*% kern(ny)
  z <- z / (2 * pi * bw^2)

  # Level = the (1 - frac) quantile of the field AT the unit positions.
  vals <- z[cbind(ix, iy)]
  lvl  <- as.numeric(stats::quantile(vals, 1 - frac, na.rm = TRUE))
  mask <- z >= lvl
  out <- list(support = sum(mask) * res^2, level = lvl, bw = bw, res = res,
              frac = frac, mode = "kde")
  if (isTRUE(return_mask)) {
    out$mask   <- mask
    out$origin <- c(gx[1], gy[1])
    out$win    <- rbind(c(xr[1], yr[1]), c(xr[2], yr[2]))
    out$h      <- res
  }
  out
}

## ═══════════════════════════════════════════════════════════════════════════════
## SUPPORT-AWARE OBSERVATION WINDOW for the spatial statistics
## ═══════════════════════════════════════════════════════════════════════════════
## WHY: Ripley's K/L, the pair correlation g and the co-localization null all need
## an observation window, and all of them used the axis-aligned BOUNDING BOX. For a
## cell that does not fill its box - a concave shape, an interior dark region, a
## scan field larger than the cell - the box overstates the area, so
##     lambda_hat = n / A  is UNDERSTATED,
## and since K = A * (pair fraction), an understated lambda inflates K at every r.
## The pattern then looks clustered because of the framing. The same error makes a
## bounding-box CSR null anti-conservative: the 2026-07-30 audit measured 53%
## (random-shift) and 93% (CSR) false-positive rates for co-localization when both
## channels shared one cell inside a larger window, against a nominal 5%, and
## identified mask-aware randomization as the outstanding fix. This is it.
##
## A `support` here is the list returned by estimate_support(..., return_mask = TRUE):
## $mask (logical array), $origin (world coords of the first cell's low corner),
## $h (cell side) and $support (the measure). `support_area()` reads the measure and
## `sample_in_support()` draws uniformly on it, so the observed statistic and its
## envelope share one geometry.
##
## MEASURED (scratchpad/test_support_window.R). A CSR pattern generated INSIDE an
## L-shaped cell whose support is 63% of its bounding box - truth is "not clustered":
##   Ripley  mean L(r)-r    bbox window +68.9 nm  ->  support window -13.5 nm
##   g(r)    mean g, r>150  bbox window 1.422     ->  support window 0.871   (want 1.0)
## Co-localization false-positive rate, two channels placed INDEPENDENTLY in one cell
## inside a larger window (truth = not co-localized, nominal 5%):
##   random-shift null   bbox 17%  ->  support 0%
##   CSR null            bbox 37%  ->  support 0%
## with power retained (10/10 genuinely co-localized pairs detected, median p 0.025).
##
## ⚠ RESIDUAL, do not oversell: the correction is only as good as the support
## estimate. Here the estimated support came out 6% larger than the true L-shape, and
## that alone leaves the observed curve outside a narrow envelope most of the time
## even though the bias in L-r is down ~80%. And 0/30 rejections cannot distinguish
## "calibrated" from "conservative" - the honest statement is that the support-aware
## null is CONSERVATIVE rather than anti-conservative, which is the safe direction.
support_area <- function(support, fallback_ranges = NULL) {
  if (!is.null(support) && is.finite(support$support %||% NA_real_) && support$support > 0)
    return(as.numeric(support$support))
  if (!is.null(fallback_ranges)) return(prod(fallback_ranges[2, ] - fallback_ranges[1, ]))
  NA_real_
}

## Uniform sample on the support: pick occupied cells with replacement, then a
## uniform position inside each. Exact uniform-on-support up to the cell size, and
## it works for any shape including one with holes - which a rejection sampler on
## the bounding box would also do but far more slowly for a thin or sparse support.
## Returns a list(x, y[, z]) of length n, or NULL if the support has no cells.
sample_in_support <- function(support, n, use_3d = FALSE) {
  if (is.null(support) || is.null(support$mask)) return(NULL)
  m <- support$mask
  d <- length(dim(m))
  occ <- which(m)
  if (!length(occ)) return(NULL)
  h  <- support$h
  o  <- support$origin
  pick <- occ[sample.int(length(occ), n, replace = TRUE)]
  ai <- arrayInd(pick, dim(m))
  # Cell (i,..) spans [origin + (i-1)h, origin + i*h] on each axis.
  jitter_axis <- function(k) o[k] + (ai[, k] - 1L) * h + stats::runif(n, 0, h)
  out <- list(x = jitter_axis(1L), y = jitter_axis(2L))
  out$z <- if (use_3d && d >= 3L) jitter_axis(3L) else NULL
  out
}

## Is a set of points inside the support? Used to restrict a randomized pattern to
## the cell (the co-localization shift null) rather than to the box.
in_support <- function(support, x, y, z = NULL) {
  if (is.null(support) || is.null(support$mask)) return(rep(TRUE, length(x)))
  m <- support$mask; dm <- dim(m); d <- length(dm); h <- support$h; o <- support$origin
  idx <- function(v, k) pmin(pmax(floor((v - o[k]) / h) + 1L, 1L), dm[k])
  ix <- idx(x, 1L); iy <- idx(y, 2L)
  if (d >= 3L && !is.null(z)) m[cbind(ix, iy, idx(z, 3L))] else {
    if (d >= 3L) {
      # 3D mask but 2D query: a column counts if any slice in it is occupied.
      p2 <- apply(m, c(1L, 2L), any); p2[cbind(ix, iy)]
    } else m[cbind(ix, iy)]
  }
}

## GOAL: compute the window's set covariance ONCE and cache it on the support object.
##
## WHY THIS EXISTS: gamma-bar depends only on the WINDOW, never on the point pattern,
## but the estimators are called once per Monte-Carlo simulation. Recomputing an FFT
## plus a million-shift radial binning inside a 99-simulation envelope loop made
## Ripley's L appear to hang - measured 0.78 s per simulation on a 210x210 mask, so
## ~77 s for one file, with the progress bar only advancing per FILE and therefore
## sitting at 0 the whole time. Attach once, reuse for the observed curve and every
## simulation. `gcache_r` is the grid gamma-bar was evaluated on; the estimators
## interpolate off it rather than rebuilding.
attach_gamma_cache <- function(support, r_max, use_3d = FALSE, n_grid = 256L) {
  if (is.null(support) || is.null(support$mask)) return(support)
  if (!is.null(support$gcache) && isTRUE(support$gcache_max >= r_max)) return(support)
  rg <- seq(0, max(r_max, support$h), length.out = as.integer(n_grid))
  gb <- tryCatch(set_covariance_iso(
    if (use_3d) support$mask else mask_project_2d(support$mask), support$h, rg),
    error = function(e) NULL)
  if (is.null(gb)) return(support)
  support$gcache     <- gb
  support$gcache_r   <- rg
  support$gcache_max <- max(rg)
  support
}

# gamma-bar at arbitrary distances, from the cache when present, else computed.
.gamma_at <- function(support, d, use_3d = FALSE) {
  if (is.null(support) || is.null(support$mask) || !length(d)) return(NULL)
  if (!is.null(support$gcache))
    return(stats::approx(support$gcache_r, support$gcache, xout = d, rule = 2)$y)
  rg <- seq(0, max(d, na.rm = TRUE), length.out = 256L)
  gb <- set_covariance_iso(
    if (use_3d) support$mask else mask_project_2d(support$mask), support$h, rg)
  if (is.null(gb)) return(NULL)
  stats::approx(rg, gb, xout = d, rule = 2)$y
}

## GOAL: a rectangular "support" over the bounding box, so the edge correction can be
## applied WITHOUT any support estimate.
##
## This matters because the edge correction and the support window fix two different
## things and only one of them is usually needed. For a cell that FILLS its scan
## field the bounding box already IS the correct window - measured on clustered data
## in a 5 um cell, bbox gave peak L-r +196 against the true window's +202 - so all
## that is missing there is the boundary correction. Requiring a support estimate to
## get it would have made the safe, always-valid fix conditional on the risky one.
bbox_support <- function(pts, use_3d = FALSE, h_nm = 100) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (!all(axes %in% names(pts))) return(NULL)
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < 3L) return(NULL)
  rng <- apply(m, 2, range)
  ext <- pmax(rng[2, ] - rng[1, ], h_nm)
  dims <- pmax(2L, as.integer(ceiling(ext / h_nm)))
  mask <- array(TRUE, dim = dims)
  list(mask = mask, origin = rng[1, ], h = h_nm,
       support = prod(dims) * h_nm^length(dims),
       win = rng, mode = "bbox", d = length(dims))
}

## GOAL: how many face-connected pieces the window mask has - the diagnostic that
## separates "a cell" from "a pile of structures".
##
## WHY IT IS NEEDED. A density-derived support answers "where is the sample dense",
## which on clustered data is the CLUSTERS, not the cell. Using that as the
## observation window conditions on the very thing a clustering statistic tests -
## inside the clusters the points look locally uniform, so the signal collapses.
## Measured on clustered data in a 5 um cell: the KDE/HDR support came to 29% of the
## cell and crushed peak L(r)-r from +202 to +49 while dragging the tail to -34,
## which is the "straight line into the negative" symptom. A footprint estimator
## (occupancy grid at its closing plateau) gave 79% and peak +163.
## A genuine cell footprint is ONE piece; a structures-hugging mask is many, and that
## is checkable without knowing the true cell.
## Returns list(n, largest_frac) - the piece count AND the share of the window that
## lives in its biggest piece.
##
## ⚠ USE `largest_frac`, NOT `n`. A raw count is the wrong criterion and rejected
## perfectly good real data: measured on four real MINFLUX acquisitions the footprint
## came out in 2, 3, 1 and 1 pieces, because sparse labelling at a cell edge leaves a
## few detached slivers even though 90%+ of the window is one cell. A structures-
## hugging mask is different in kind - many pieces of COMPARABLE size, so no single
## piece dominates. `largest_frac` separates those two cases; a count cannot.
mask_components <- function(mask, h) {
  if (is.null(mask) || !any(mask)) return(list(n = 0L, largest_frac = 0))
  idx <- which(mask)
  if (length(idx) == 1L) return(list(n = 1L, largest_frac = 1))
  sub <- arrayInd(idx, dim(mask))
  cl <- tryCatch(dbscan::dbscan(sub * h, eps = 1.05 * h, minPts = 1L)$cluster,
                 error = function(e) NULL)
  if (is.null(cl)) return(list(n = NA_integer_, largest_frac = NA_real_))
  tb <- table(cl[cl > 0])
  if (!length(tb)) return(list(n = 0L, largest_frac = 0))
  list(n = length(tb), largest_frac = max(tb) / sum(tb))
}

## GOAL: the isotropised SET COVARIANCE of the support, which is what turns Ripley's
## K and the pair correlation g into unbiased estimators in a window of ANY shape.
##
## WHY THIS IS THE REAL FIX, AND WHY THE AREA WAS NOT. Both estimators divide by an
## area, so a wrong area is a MULTIPLICATIVE scale error. But even with the exactly
## correct area they are still biased, because a point near the boundary has fewer
## neighbours available at separation r than one in the middle - nothing about the
## area fixes that. Measured on CSR (truth g = 1) inside an L-shaped cell:
##     bounding-box area, no correction   1.398
##     support area, no correction        0.857
##     EXACT true area, no correction     0.827   <- the area was never the problem
##     exact area + this correction       1.003
##     support area + this correction     1.030
## and on a plain square with the exact area, uncorrected g decays 0.949 -> 0.783
## across r (mean 0.882), corrected 0.984. The bounding box only ever looked closer
## to 1 because its area error and the edge deficit had opposite signs and partly
## cancelled; enabling the support removed one and exposed the other, which is why
## the curves started looking wrong.
##
## THE MATHS. For a window W, gamma(h) = |W ∩ W_h| is the set covariance. The number
## of ordered pairs with separation in a shell S has expectation
##     lambda^2 * integral_S g(h) gamma(h) dh  ~=  lambda^2 g(r) gammabar(r) |S|,
## so the estimator must divide by gammabar(r), not by |W|. Substituting gamma = |W|
## recovers the uncorrected form, i.e. the uncorrected estimator is the special case
## that assumes NO area is lost to the boundary. This is the translation (Ohser)
## correction; see Illian et al. 2008 §4.3.3, Baddeley et al. 2015 §7.4.
##
## COMPUTATION. gamma is the autocorrelation of the mask indicator, so one FFT gives
## it for every shift at once. The array is ZERO-PADDED to twice its size per axis -
## without that the FFT autocorrelation is circular and wraps opposite edges
## together. gamma(0) must come back equal to the support measure; that is asserted
## by the callers' tests.
##
## RADIAL BINNING: shifts live on a lattice of spacing h, so a shell thinner than h
## can contain NO lattice point. Binning at the requested r-spacing therefore leaves
## most shells empty (this bit me in the prototype: NAs at scattered radii, which
## silently dropped those radii from the curve). gammabar is therefore accumulated in
## bins of width h and then INTERPOLATED onto the requested radii.
set_covariance_iso <- function(mask, h, r_out) {
  if (is.null(mask) || !any(mask) || !length(r_out)) return(NULL)
  dm <- dim(mask); nd <- length(dm)
  if (nd < 2L || nd > 3L) return(NULL)
  pad <- vapply(dm, function(k) 2^ceiling(log2(2 * k)), numeric(1))
  # Guard against an absurd FFT (a very fine 3D grid); fall back to no correction.
  if (prod(pad) > 3e7) return(NULL)
  m <- array(0, dim = pad)
  m[slice.index(m, 1) <= dm[1]] <- 0   # (initialised zero; explicit for clarity)
  idx <- lapply(dm, seq_len)
  m <- do.call("[<-", c(list(m), idx, list(value = as.numeric(mask))))

  F  <- stats::fft(m)
  Cc <- Re(stats::fft(F * Conj(F), inverse = TRUE)) / prod(pad)   # overlap in CELLS

  # Shift lattice per axis: index 1 is shift 0, upper half is negative shifts.
  shifts <- lapply(pad, function(k) c(0:(k / 2 - 1), -(k / 2):(-1)) * h)
  d2 <- if (nd == 2L) outer(shifts[[1]]^2, shifts[[2]]^2, "+") else
    outer(outer(shifts[[1]]^2, shifts[[2]]^2, "+"), shifts[[3]]^2, "+")
  dist <- sqrt(d2)

  # Bin at the lattice spacing, then interpolate - see RADIAL BINNING above.
  rmax <- max(r_out) + h
  brk  <- seq(0, rmax + h, by = h)
  nb   <- length(brk) - 1L
  bin  <- findInterval(as.vector(dist), brk, rightmost.closed = TRUE)
  keep <- bin >= 1L & bin <= nb
  cv   <- as.vector(Cc)[keep]; bin <- bin[keep]
  if (!length(bin)) return(NULL)
  # rowsum is C-level; tapply over ~1e6 lattice shifts was a measurable cost.
  rs_  <- rowsum(cv, bin, reorder = TRUE)
  cnt  <- tabulate(bin, nbins = nb)
  ok   <- as.integer(rownames(rs_))
  sums <- as.numeric(rs_[, 1]) / cnt[ok]
  ctr  <- (brk[-1] + brk[-length(brk)]) / 2
  if (length(ok) < 2L) return(NULL)
  g_at <- stats::approx(ctr[ok], as.numeric(sums), xout = r_out, rule = 2)$y
  # (binning above uses rowsum, which is C-level; tapply over ~1e6 lattice shifts
  #  was a measurable share of the cost when this ran inside an envelope loop.)
  cell <- h^nd
  gb   <- g_at * cell                                   # cells -> nm^2 / nm^3
  # gamma is non-increasing in |h| for a convex window and bounded by |W| always;
  # clamp tiny numerical excursions so the correction can never shrink the estimate.
  gb <- pmin(pmax(gb, .Machine$double.eps), sum(mask) * cell)
  gb
}

## GOAL: turn an occupancy mask into drawable OUTLINE segments, so the support can
## be shown as a boundary over the localizations rather than as a number.
##
## An edge is on the boundary when an occupied cell's neighbour is not occupied (or
## is off the grid). Emitting exactly those edges outlines the support including any
## interior holes, which is what makes a swallowed dark region visible at a glance.
## Segments are returned as one x/y pair with NA separators, so the whole outline is
## a SINGLE plotly trace however many disjoint pieces it has.
##
## `mask` is 2D (project a 3D mask first - see mask_project_2d). `origin` is the
## world coordinate of the first cell's low corner, `h` the cell side.
mask_outline_segments <- function(mask, origin, h) {
  if (is.null(mask) || !any(mask)) return(list(x = numeric(0), y = numeric(0)))
  d <- dim(mask)
  occ <- which(mask, arr.ind = TRUE)          # columns: row = x index, col = y index
  ix <- occ[, 1L]; iy <- occ[, 2L]
  # Cell (i,j) spans [x0 + (i-1)h, x0 + i*h] x [y0 + (j-1)h, y0 + j*h].
  x_lo <- origin[1] + (ix - 1L) * h; x_hi <- x_lo + h
  y_lo <- origin[2] + (iy - 1L) * h; y_hi <- y_lo + h

  # Is the neighbour occupied? Off-grid counts as not occupied.
  nb <- function(dx, dy) {
    jx <- ix + dx; jy <- iy + dy
    ok <- jx >= 1L & jx <= d[1] & jy >= 1L & jy <= d[2]
    out <- rep(FALSE, length(ix))
    if (any(ok)) out[ok] <- mask[cbind(jx[ok], jy[ok])]
    out
  }
  xs <- ys <- numeric(0)
  push <- function(x1, y1, x2, y2, sel) {
    if (!any(sel)) return(invisible(NULL))
    xs <<- c(xs, as.numeric(rbind(x1[sel], x2[sel], NA)))
    ys <<- c(ys, as.numeric(rbind(y1[sel], y2[sel], NA)))
  }
  push(x_lo, y_lo, x_hi, y_lo, !nb(0L, -1L))   # bottom edge
  push(x_lo, y_hi, x_hi, y_hi, !nb(0L,  1L))   # top edge
  push(x_lo, y_lo, x_lo, y_hi, !nb(-1L, 0L))   # left edge
  push(x_hi, y_lo, x_hi, y_hi, !nb( 1L, 0L))   # right edge
  list(x = xs, y = ys)
}

# Collapse a 3D occupancy mask onto XY: a column counts as sample if ANY z slice in
# it is occupied. The 2D map is a projection, so a 3D support looks larger in plan
# view than it is - stated in the UI rather than silently.
mask_project_2d <- function(mask) {
  if (is.null(mask)) return(NULL)
  if (length(dim(mask)) < 3L) return(mask)
  apply(mask, c(1L, 2L), any)
}

## GOAL: make the gap-fill (closing) radius CHOOSABLE instead of guessed, by
## sweeping it and reporting the resulting support.
##
## The curve has a characteristic shape and the shape is the whole point:
##   * at close_nm = 0 the support is just the occupied cells, so it tracks the
##     STRUCTURES rather than the cell - strongly n-dependent and not a density
##     denominator you want;
##   * it then RISES steeply as sampling gaps inside the sample get filled;
##   * it PLATEAUS once the sample region is solid - this plateau is the footprint,
##     and its left edge is the radius to use;
##   * it rises again only if the radius grows big enough to bridge genuine dark
##     regions or reach across the frame, which is what you must NOT do.
## Pick the radius at the start of the plateau. A support that is still climbing
## steeply at your chosen radius means the estimate is not yet a footprint; one
## that has climbed to near the bounding box means real dark regions are being
## swallowed.
##
## Returns data.frame(close_nm, support, frac_of_bbox, rel_step) where rel_step is
## the fractional increase from the previous radius - the plateau is where
## rel_step falls near 0. `plateau_close_nm` reads the plateau start off this.
##
## MEASURED, do not re-derive: on two synthetic cells with identical cluster
## geometry and a 4.0x labelling-density difference, a 200 nm radius (the old
## default) gave supports of 30% and 24% of frame and a density ratio of 3.13,
## while the plateau at 600-800 nm gave 71% and 70% - agreeing to 1 percentage
## point - and a ratio of 3.94. Choosing this radius badly is the single largest
## error source in the whole density normalisation, which is why it is swept.
support_closing_sweep <- function(pts, use_3d = FALSE, h_nm = 100,
                                  close_seq = NULL) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  if (is.null(close_seq)) close_seq <- h_nm * c(0, 1, 2, 3, 4, 6, 8, 12, 16, 24)
  close_seq <- sort(unique(pmax(0, as.numeric(close_seq))))
  out <- lapply(close_seq, function(cn) {
    s <- tryCatch(estimate_support(pts, use_3d = use_3d, mode = "grid",
                                   h_nm = h_nm, close_nm = cn),
                  error = function(e) NULL)
    if (is.null(s)) return(NULL)
    data.frame(close_nm = cn, support = s$support, frac_of_bbox = s$frac_of_bbox)
  })
  res <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(res) || !nrow(res)) return(NULL)
  res$rel_step <- c(NA_real_, diff(res$support) / utils::head(res$support, -1))
  res
}

## GOAL: read the plateau start off a closing sweep - the smallest radius at which
## the support has stopped growing appreciably, i.e. sampling gaps are filled but
## dark regions are not yet bridged.
##
## `tol` is the fractional support increase below which the curve counts as flat.
## Guards: never returns a radius whose support already fills more than
## `max_frac` of the bounding box (that is the "swallowed the dark regions" end of
## the curve), and falls back to the largest swept radius below that cap.
plateau_close_nm <- function(sweep, tol = 0.05, max_frac = 0.97) {
  if (is.null(sweep) || !nrow(sweep)) return(NA_real_)
  s <- sweep[order(sweep$close_nm), , drop = FALSE]
  ok <- is.finite(s$frac_of_bbox) & s$frac_of_bbox <= max_frac
  if (!any(ok)) return(NA_real_)
  s <- s[ok, , drop = FALSE]
  hit <- which(is.finite(s$rel_step) & s$rel_step < tol)
  if (length(hit)) return(s$close_nm[hit[1]])
  s$close_nm[nrow(s)]
}

## GOAL: a SUPPORT-FREE intensity estimate, from the k-th nearest-neighbour
## distances of the units themselves.
##
## For a Poisson process of intensity lambda in R^d, the k-th NN distance R_k
## satisfies lambda * omega_d * R_k^d ~ Gamma(k, 1), where omega_d is the volume
## of the unit ball (omega_2 = pi, omega_3 = 4pi/3). Hence
##     E[R_k] = Gamma(k + 1/d) / (Gamma(k) * (lambda * omega_d)^(1/d))
## which inverts to the moment estimator
##     lambda_hat = ( Gamma(k + 1/d) / (Gamma(k) * mean(R_k)) )^d / omega_d.
##
## WHY IT IS USEFUL: it needs no window and no support estimate, and because it
## only measures where points actually ARE, it is automatically blind to the frame
## edge, to the cell boundary and to interior dark regions. That makes it the
## natural cross-check on estimate_support: if the grid support is sane, the two
## density numbers agree to within a factor of order one.
##
## ⚠ WHY IT IS NOT THE THINNING TARGET: it is a LOCAL intensity, so clustering
## inflates it. An aggregating cell has short NN distances and reads as denser
## than it is, and thinning to match it would remove part of the aggregation
## signal. `probs` therefore defaults to the UPPER-MIDDLE of the NN distribution
## (0.5-0.9), which is dominated by the sparse background rather than by the
## aggregates - a partial, not a complete, mitigation. Use the global density for
## thinning; use this to diagnose.
##
## CALIBRATION OF THE QUANTILE BAND (do not remove - it is worth 30%):
## the Gamma relation above gives E[R_k] over ALL points. Averaging R_k over a
## quantile band instead raises the mean, so feeding a band mean into the
## full-distribution formula under-estimates lambda by ~30% at probs = (0.5, 0.9).
## Measured before the fix: ratio to truth 0.64 / 0.70 / 0.74 across three decades
## of lambda in 2D and 0.63 in 3D - a systematic bias, not noise. The band is
## therefore calibrated against its OWN conditional expectation,
##     c(k, d, a, b) = E[ G^(1/d) | G in (qgamma(a, k), qgamma(b, k)) ],  G ~ Gamma(k, 1)
## evaluated numerically, so lambda_hat = (c / mean(R_k | band))^d / omega_d is
## unbiased on CSR for any band. Verified back to ratio ~1.0 (see test 3).
##
## Returns list(lambda, mean_r, median_r, k, d, n) with lambda in units per nm^d.
nn_local_intensity <- function(pts, k = 3L, use_3d = FALSE, probs = c(0.5, 0.9)) {
  if (is.null(pts) || !nrow(pts)) return(NULL)
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m); d <- ncol(m)
  if (n < 3L) return(NULL)
  k_eff <- min(max(1L, as.integer(k)), n - 1L)
  nn <- tryCatch(RANN::nn2(m, m, k = k_eff + 1L), error = function(e) NULL)
  if (is.null(nn)) return(NULL)
  rk <- nn$nn.dists[, k_eff + 1L]        # column 1 is the point itself (d = 0)
  rk <- rk[is.finite(rk) & rk > 0]
  if (length(rk) < 3L) return(NULL)

  # Restrict to the quantile band that the sparse background dominates.
  pr  <- sort(probs)
  qq  <- stats::quantile(rk, probs = pr, na.rm = TRUE)
  sel <- rk >= qq[1] & rk <= qq[2]
  used_band <- sum(sel) >= 3L
  rb  <- if (used_band) rk[sel] else rk

  omega <- if (d == 3L) 4 * pi / 3 else pi
  mr    <- mean(rb, na.rm = TRUE)
  if (!is.finite(mr) || mr <= 0) return(NULL)

  # Calibration constant for the band actually used (see the note above): the
  # conditional mean of G^(1/d) over the matching Gamma quantile band. Falls back
  # to the closed-form full-distribution constant if the integration fails.
  cst <- if (!used_band || (pr[1] <= 0 && pr[2] >= 1)) {
    gamma(k_eff + 1 / d) / gamma(k_eff)
  } else {
    ga <- stats::qgamma(pr[1], shape = k_eff)
    gb <- stats::qgamma(pr[2], shape = k_eff)
    val <- tryCatch(
      stats::integrate(function(g) g^(1 / d) * stats::dgamma(g, shape = k_eff),
                       lower = ga, upper = gb, rel.tol = 1e-8)$value / (pr[2] - pr[1]),
      error = function(e) NA_real_)
    if (is.finite(val) && val > 0) val else gamma(k_eff + 1 / d) / gamma(k_eff)
  }

  lam <- (cst / mr)^d / omega
  list(lambda = lam, mean_r = mr, median_r = stats::median(rk, na.rm = TRUE),
       k = k_eff, d = d, n = n, probs = pr, band_const = cst, banded = used_band)
}

## GOAL: the PER-UNIT local density, one value per unit, for the distribution plot.
##
## From the same Gamma law as nn_local_intensity but read per point rather than
## averaged: lambda_i = k / (omega_d * R_k,i^d) is the unbiased per-point estimator
## of the local intensity (E[R_k^d] = k / (lambda * omega_d)).
##
## WHY A DISTRIBUTION AND NOT JUST THE MEAN: the single global number n/support says
## how much label there is; the SHAPE says how it is arranged. A cell with genuine
## aggregates is heterogeneous - a sparse background mode plus a high-density tail -
## whereas a uniformly-labelled cell of the same mean density is unimodal. So the
## spread of this distribution is itself clustering information, which is exactly
## why its MEAN must not be used as the thinning target (see nn_local_intensity).
## Returned in units per um^d to match compute_unit_density's `density`.
unit_local_density <- function(pts, k = 3L, use_3d = FALSE) {
  if (is.null(pts) || !nrow(pts)) return(numeric(0))
  axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
  m <- as.matrix(pts[, axes, drop = FALSE])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m); d <- ncol(m)
  if (n < 3L) return(numeric(0))
  k_eff <- min(max(1L, as.integer(k)), n - 1L)
  nn <- tryCatch(RANN::nn2(m, m, k = k_eff + 1L), error = function(e) NULL)
  if (is.null(nn)) return(numeric(0))
  rk <- nn$nn.dists[, k_eff + 1L]          # column 1 is the point itself
  ok <- is.finite(rk) & rk > 0
  if (!any(ok)) return(numeric(0))
  omega <- if (d == 3L) 4 * pi / 3 else pi
  (k_eff / (omega * rk[ok]^d)) * 1000^d
}

## GOAL: the per-file row of the density panel - everything needed to decide
## whether two acquisitions are comparable, and to thin them if they are not.
##
## Returns a one-row data.frame:
##   n_units      units in the file (trace IDs or de-blinked molecules)
##   support_nm   estimated support (nm^2 or nm^3), by `support_mode`
##   density      units per um^2 (2D) or per um^3 (3D)   <- the THINNING TARGET
##   density_nn   the same quantity from nn_local_intensity (support-free, but
##                clustering-confounded - a cross-check, not a target)
##   nn_median    median k-th NN distance (nm), the raw observable
##   frac_of_bbox support / bounding box - a framing/edge diagnostic. Well below
##                1 means a lot of empty frame or dark region was excluded, which
##                is the case where the bounding-box density would have been wrong.
##
## ⚠ `support_override`: THE SUPPORT ESTIMATE IS ITSELF DENSITY-DEPENDENT at low n.
## An occupancy grid built from fewer points finds fewer occupied cells, so
## re-estimating the support after thinning shrinks A and inflates n/A - which
## would make the thinning look like it had failed to equalise anything (measured:
## a 4x labelling difference read as 3.4x before thinning and still 1.16x after,
## purely from this). A cell's footprint does not change when you thin it, so pass
## the PRE-THINNING support here whenever the pattern has been thinned. The app's
## thinning gate does exactly that.
##
## ⚠ `support_from`: WHICH POINTS DEFINE THE FOOTPRINT. The support is a property
## of WHERE SIGNAL WAS DETECTED, so it is best estimated from EVERY LOCALIZATION
## ("locs", the default) rather than from the unit centroids. Using centroids makes
## the support depend on the unit count, which is the very thing being normalised:
## measured on two synthetic cells with identical geometry and a 4x labelling
## difference, centroid-based supports were 5.01e6 vs 3.82e6 nm^2 for the same
## region, so the density ratio read 3.05 instead of 4.0.
compute_unit_density <- function(df, unit = "tid", use_3d = FALSE,
                                support_mode = "grid", h_nm = 100, close_nm = 200,
                                k = 3L, support_override = NULL,
                                support_from = c("locs", "units"),
                                bw_nm = 150, frac = 0.99) {
  support_from <- match.arg(support_from)
  pts <- unit_centroids(df, unit = unit, use_3d = use_3d)
  if (is.null(pts) || nrow(pts) < 3L) return(NULL)
  sup_pts <- if (support_from == "locs") {
    keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
    if (!all(keep %in% names(df))) pts else {
      d2 <- df[, keep, drop = FALSE]
      if (!use_3d) d2$z <- 0
      d2[stats::complete.cases(d2), , drop = FALSE]
    }
  } else pts
  sup <- if (!is.null(support_override) && is.finite(support_override) && support_override > 0) {
    list(support = as.numeric(support_override), mode = "inherited (pre-thinning)",
         h = h_nm, close_nm = close_nm, frac_of_bbox = NA_real_)
  } else {
    estimate_support(sup_pts, use_3d = use_3d, mode = support_mode,
                     h_nm = h_nm, close_nm = close_nm,
                     bw_nm = bw_nm, frac = frac)
  }
  if (is.null(sup) || !is.finite(sup$support) || sup$support <= 0) return(NULL)
  nni <- nn_local_intensity(pts, k = k, use_3d = use_3d)

  # Dimensionality follows the ESTIMATOR, not the 3-D checkbox: the KDE/HDR support
  # is a plan-view area even on 3-D data, so treating it as a volume would divide by
  # nm^3 and report a density three orders of magnitude wrong.
  d      <- if (!is.null(sup$d)) as.integer(sup$d) else if (use_3d) 3L else 2L
  per_um <- 1000^d                       # nm^d -> um^d
  data.frame(
    n_units      = nrow(pts),
    support_nm   = sup$support,
    support_mode = sup$mode,
    h_nm         = sup$h %||% NA_real_,
    close_nm     = sup$close_nm %||% NA_real_,
    frac_of_bbox = sup$frac_of_bbox,
    density      = nrow(pts) / sup$support * per_um,
    density_nn   = if (!is.null(nni)) nni$lambda * per_um else NA_real_,
    nn_median    = if (!is.null(nni)) nni$median_r else NA_real_,
    nn_k         = if (!is.null(nni)) nni$k else NA_integer_,
    dims         = d,
    stringsAsFactors = FALSE
  )
}

## GOAL: independent p-thinning of a localization table, AT THE UNIT LEVEL.
##
## Each unit (trace ID, or de-blinked molecule) is kept or dropped as a WHOLE,
## with all of its localizations. Thinning individual localizations instead would
## corrupt the thing being normalised: it would shorten traces, change their
## spread and precision, and break de-blinking - and it would not thin the point
## pattern at all, since a partially-sampled trace still contributes one centroid.
##
## `p` >= 1 is a no-op (never up-samples - you cannot invent molecules).
## The draw is deterministic given `seed`, and the GLOBAL RNG stream is SAVED AND
## RESTORED around it, so enabling thinning cannot perturb any other randomised
## result in the app (the permutation nulls, the CSR envelopes, the alpha-shape
## subsampling all draw from the same stream).
thin_units_to_p <- function(df, unit = "tid", p = 1, seed = 1L) {
  if (is.null(df) || !nrow(df)) return(df)
  if (!is.finite(p) || p >= 1) return(df)
  if (p <= 0) return(df[0, , drop = FALSE])
  key <- if (unit == "cluster" && "cluster" %in% names(df)) "cluster" else "tid"
  if (!key %in% names(df)) return(df)

  ids <- unique(df[[key]])
  # De-blinking noise (id 0 / NA) is not a unit; it is dropped from the thinned
  # pattern rather than kept wholesale, which would leave a density floor that no
  # thinning could reach.
  if (key == "cluster") ids <- ids[!is.na(ids) & ids > 0L]
  if (!length(ids)) return(df[0, , drop = FALSE])

  had  <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_backup <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", seed_backup, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  set.seed(seed)
  # Fixed-count draw rather than n independent coin flips: same expectation, no
  # binomial noise in the achieved density, which matters because the whole point
  # is to land every file on ONE common density.
  n_keep <- max(1L, as.integer(round(length(ids) * p)))
  keep   <- sample(ids, n_keep)
  df[df[[key]] %in% keep, , drop = FALSE]
}

# =============================================================================
# dilate_mask_cells() - grow a boolean mask by r whole cells (Chebyshev/box)
# =============================================================================
# Separable box dilation, the same structuring element estimate_support() closes
# with. The array is padded by r on every side first so the growth has room, so
# the returned mask is larger than the input and `origin` shifts by -r*h.
# RETURNS list(mask, pad) with pad = r (cells added on each side).
dilate_mask_cells <- function(mask, r) {
  r <- max(0L, as.integer(r))
  if (r == 0L) return(list(mask = mask, pad = 0L))
  dm <- dim(mask); d <- length(dm)
  if (!d %in% c(2L, 3L)) return(list(mask = mask, pad = 0L))
  # pad by r on every side so the growth has room
  big <- array(FALSE, dim = dm + 2L * r)
  if (d == 2L) big[seq_len(dm[1]) + r, seq_len(dm[2]) + r] <- mask
  else         big[seq_len(dm[1]) + r, seq_len(dm[2]) + r, seq_len(dm[3]) + r] <- mask
  db <- dim(big)
  # separable: OR the array with itself shifted by one cell, r times per axis
  for (ax in seq_len(d)) {
    n_ax <- db[ax]
    if (n_ax < 2L) next
    for (step in seq_len(r)) {
      lo <- seq_len(n_ax - 1L); hi <- lo + 1L
      if (d == 2L) {
        if (ax == 1L) big[lo, ] <- big[lo, ] | big[hi, , drop = FALSE]
        else          big[, lo] <- big[, lo] | big[, hi, drop = FALSE]
        if (ax == 1L) big[hi, ] <- big[hi, ] | big[lo, , drop = FALSE]
        else          big[, hi] <- big[, hi] | big[, lo, drop = FALSE]
      } else {
        if (ax == 1L) { big[lo, , ] <- big[lo, , ] | big[hi, , , drop = FALSE]
                        big[hi, , ] <- big[hi, , ] | big[lo, , , drop = FALSE] }
        else if (ax == 2L) { big[, lo, ] <- big[, lo, ] | big[, hi, , drop = FALSE]
                             big[, hi, ] <- big[, hi, ] | big[, lo, , drop = FALSE] }
        else { big[, , lo] <- big[, , lo] | big[, , hi, drop = FALSE]
               big[, , hi] <- big[, , hi] | big[, , lo, drop = FALSE] }
      }
    }
  }
  list(mask = big, pad = r)
}

# =============================================================================
# debias_support_window() - grow the footprint to the region the points really fill
# =============================================================================
#
# GOAL
#   Remove the area bias of the occupancy-grid footprint before it is used as the
#   observation window of a spatial statistic.
#
# WHY IT IS NEEDED
#   The footprint is built from occupied grid cells, so at sparse labelling its
#   boundary sits INSIDE the true cell edge: measured on CSR in an L-shaped cell,
#   estimated/true area = 0.87 / 0.93 / 0.99 / 1.00 at 5.6 / 14 / 42 / 139 units per
#   um^2. Ripley L, g(r), L12 and g12 compare the observed pattern with simulations
#   drawn INSIDE that window. The observed points fill the real cell, the simulated
#   ones only the smaller window, so random data look more spread out than the null
#   and read as regular: L-r at r = 200 nm came out -6.9 nm against a null mean of
#   -0.7 nm, with 33% of random replicates below the pointwise band.
#
# HOW
#   The shrink is measured on the window itself: draw `n_pts` uniform points inside
#   the mask, re-estimate the footprint from them with the same pitch and gap fill,
#   and average the ratio rho = A_re-estimated / A_mask over `n_draw` draws. The cell
#   the observed points came from is then about A_mask / rho, so the mask is dilated
#   by whole cells until its area is closest to that target (capped at the bounding
#   box). Dilation, not rescaling: the simulations have to be drawn from the larger
#   REGION, and only a mask can do that.
#
# support : estimate_support() object WITH $mask
# n_pts   : number of points the window was estimated from
# RETURNS the support with mask, origin, support (area), frac_of_bbox and
#   debias_rho / debias_cells updated; unchanged when there is no mask or rho ~ 1.
#   Any gamma cache is dropped, so re-attach attach_gamma_cache() afterwards.
debias_support_window <- function(support, n_pts, n_draw = 6L, max_cells = 8L, tol = 0.02) {
  if (is.null(support) || is.null(support$mask) || !is.finite(support$h %||% NA_real_)) return(support)
  if (!is.finite(n_pts) || n_pts < 10L) return(support)
  h <- support$h; d <- length(dim(support$mask))
  A0 <- support$support
  if (!is.finite(A0) || A0 <= 0) return(support)
  rhos <- vapply(seq_len(max(1L, n_draw)), function(i) {
    p <- tryCatch(sample_in_support(support, n_pts, use_3d = d >= 3L), error = function(e) NULL)
    if (is.null(p)) return(NA_real_)
    pts <- data.frame(x = p$x, y = p$y, z = if (!is.null(p$z)) p$z else 0)
    e <- tryCatch(estimate_support(pts, use_3d = d >= 3L, mode = support$mode %||% "grid",
                                   h_nm = h, close_nm = support$close_nm %||% 200),
                  error = function(err) NULL)
    a <- e$support %||% NA_real_
    if (is.finite(a) && a > 0) a / A0 else NA_real_
  }, numeric(1))
  rho <- mean(rhos, na.rm = TRUE)
  if (!is.finite(rho) || rho <= 0 || rho >= 1 - tol) return(support)   # already unbiased
  target <- min(A0 / rho, support$bbox %||% Inf)
  cell_v <- h^d
  best <- support; best_diff <- abs(A0 - target); best_r <- 0L
  for (r in seq_len(max_cells)) {
    dl <- dilate_mask_cells(support$mask, r)
    a  <- sum(dl$mask) * cell_v
    if (abs(a - target) < best_diff) {
      best_diff <- abs(a - target); best_r <- r
      best <- support
      best$mask   <- dl$mask
      best$origin <- support$origin - r * h
      best$support <- a
      best$frac_of_bbox <- if (is.finite(support$bbox %||% NA_real_)) a / support$bbox else NA_real_
    }
    if (a >= target) break
  }
  best$gcache <- NULL; best$gcache_r <- NULL      # the covariance belongs to the old mask
  best$debias_rho <- rho; best$debias_cells <- best_r
  best
}


# =============================================================================
# density_profile_files() - PER-FILE LABELLING DENSITY ON ONE COMMON SUPPORT
# =============================================================================
# GOAL: the numbers behind the labelling-density panel. For every prepared file:
# the closing-radius sweep (optional), then ONE gap-fill radius shared by all
# files (the largest plateau start when auto), then each file's unit density on
# that radius and a per-unit local-density sample for the distribution plot.
#
# The radius must be the same for every file: different denominators measure
# different kinds of region, and that is exactly what makes a density ratio
# wrong. The largest plateau is taken because a radius below a file's plateau
# leaves its support tracking its structures.
#
# ARGUMENTS
#   prepped       named list of filtered (and de-blinked when unit = "cluster")
#                 tables, NOT thinned
#   unit          "cluster" (de-blinked molecules) or "tid"
#   use_3d        volume density instead of area density
#   support_mode  "grid" or a KDE mode (see compute_unit_density)
#   h_nm, close_nm, k, support_from, bw_nm, frac
#                 passed to compute_unit_density; close_nm is overridden when
#                 auto_close finds plateaus
#   auto_close    pick close_nm from the sweeps (grid support only)
#   want_sweep    compute the sweeps (auto_close needs them)
#   progress      optional function(frac, detail), called once per density file
#
# RETURNS list(profile = one row per file or NULL, sweeps = bound sweep table or
#              NULL, locals = data.frame(dataset, local) or NULL,
#              close_nm = the radius used)
# =============================================================================
density_profile_files <- function(prepped, unit = "tid", use_3d = FALSE,
                                  support_mode = "grid", h_nm = 100, close_nm = 200,
                                  auto_close = FALSE, want_sweep = auto_close, k = 3L,
                                  support_from = "units", bw_nm = 150, frac = 0.99,
                                  progress = NULL) {
  # sweep on the same points the support is built from (all localizations), so
  # the plateau the user reads is the plateau the density uses
  sweeps <- list()
  if (want_sweep) {
    keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
    for (fname in names(prepped)) {
      df <- prepped[[fname]]
      if (!all(keep %in% names(df))) next
      sp <- df[, keep, drop = FALSE]
      if (!use_3d) sp$z <- 0
      sw <- tryCatch(support_closing_sweep(sp[stats::complete.cases(sp), , drop = FALSE],
                                           use_3d = use_3d, h_nm = h_nm),
                     error = function(e) NULL)
      if (!is.null(sw)) { sw$dataset <- fname; sweeps[[fname]] <- sw }
    }
  }

  # one common gap-fill radius
  if (auto_close && length(sweeps)) {
    pls <- vapply(sweeps, function(s) plateau_close_nm(s), numeric(1))
    pls <- pls[is.finite(pls)]
    if (length(pls)) close_nm <- max(pls)
  }

  # densities, all on that one radius
  rows <- list(); locals <- list()
  for (fname in names(prepped)) {
    if (is.function(progress)) progress(0.4 / max(1L, length(prepped)), fname)
    df <- prepped[[fname]]
    r <- tryCatch(compute_unit_density(df, unit = unit, use_3d = use_3d,
                                      support_mode = support_mode, h_nm = h_nm,
                                      close_nm = close_nm, k = k,
                                      support_from = support_from,
                                      bw_nm = bw_nm, frac = frac),
                  error = function(e) NULL)
    if (is.null(r)) next
    r$dataset    <- fname
    r$unit       <- unit
    r$close_auto <- auto_close
    # the settings this row was measured with, so the status line and the plot
    # describe the run rather than the live sliders
    r$run_mode   <- support_mode
    r$run_bw     <- bw_nm
    r$run_frac   <- frac
    r$run_3d     <- use_3d
    rows[[fname]] <- r
    # per-unit local density, subsampled to 4000: the plot needs only the shape
    lv <- tryCatch(unit_local_density(unit_centroids(df, unit, use_3d),
                                     k = k, use_3d = use_3d),
                   error = function(e) numeric(0))
    lv <- lv[is.finite(lv) & lv > 0]
    if (length(lv) > 4000L) lv <- lv[seq(1, length(lv), length.out = 4000L)]
    if (length(lv)) locals[[fname]] <- data.frame(dataset = fname, local = lv)
  }
  list(profile  = do.call(rbind, rows),
       sweeps   = if (length(sweeps)) do.call(rbind, sweeps) else NULL,
       locals   = if (length(locals)) do.call(rbind, locals) else NULL,
       close_nm = close_nm)
}
