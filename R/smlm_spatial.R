## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_spatial.R  -  Spatial statistics and co-localization
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Ripley's L and L12, the pair correlation g(r) and g12(r), their
##   randomisation envelopes and global test (CSR, random labelling, toroidal
##   shift), the correlation length, radius of gyration and the KDE
##   co-localization test.
##
## CONTENTS (in file order)
##   compute_cluster_rog, .stable_subsample, compute_ripley_L,
##   compute_ripley_L12, compute_pcf12, fit_correlation_length, compute_pcf,
##   .shift_channel, compute_ripley_envelope, envelope_global_p,
##   compute_kde_coloc, .add_global_band, .global_p_caption, .global_p_annot,
##   .rmax_axis_title
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

## ── Radius of Gyration ───────────────────────────────────────────────────────
# Returns a data.frame with one row per cluster: cluster_id, n_traces, rog (nm)
# ── compute_cluster_rog() ─────────────────────────────────────────────────────
# GOAL: radius of gyration  Rg = sqrt(<|p - centroid|^2>)  per cluster - the RMS
# size of a point cloud about its centre.
#
# TWO PARAMETERS, TWO QUESTIONS:
#   cluster_col = what a "cluster" IS (the grouping):
#       "cluster2" -> a median cluster -> Rg measures STRUCTURE SIZE.
#       "cluster"  -> one de-blinked molecule -> Rg measures that molecule's own
#                     blink/label FOOTPRINT (spread of its repeat detections).
#   member_col  = the UNIT reduced to one point before measuring spread:
#       "tid"      -> raw traces (one median point per trace).
#       "cluster"  -> de-blinked molecules (one median point per molecule). Use
#                     this for structure size when de-blinking is ON, so a
#                     multivalent label's several traces count as ONE molecule
#                     rather than inflating Rg with several nearby points.
#
# MINFLUX adaptation: points are MEDIAN member positions, never raw localizations.
# A median of N localizations has SD ~ sigma_loc/sqrt(N), so localization noise is
# already suppressed and Rg reflects real geometry, not measurement scatter. (What
# remains is the label footprint - real between-trace spread - which is exactly the
# thing the "cluster" grouping is meant to report and the "cluster2" grouping is
# meant to look past by de-blinking to one point per molecule.)
#
# SMALL-N BIAS (bias_correct): the centroid is estimated from the SAME points whose
# spread is being measured, so it sits inside the cloud and E[Rg²_obs] = Rg²·(N−1)/N.
# That is a 29% underestimate at N=2, 11% at N=5, 1% at N=50 - i.e. an N-DEPENDENT
# bias, which manufactures a spurious "bigger clusters have bigger Rg" trend whenever
# clusters of different member counts are pooled or plotted against size. Multiplying
# the mean squared distance by N/(N−1) (Bessel, exactly as for a sample variance)
# removes it, giving an unbiased estimate of the underlying structure's Rg. Set
# bias_correct = FALSE for the raw configurational Rg of the observed point set.
#
# Requires >= 2 members per cluster (Rg is undefined for one point). Noise ids (0)
# are dropped for both the grouping and, when members are a cluster level, the
# member. Returns data.frame(cluster_id, n_units, rog).  nm throughout.
compute_cluster_rog <- function(df, cluster_col = "cluster2", member_col = "tid",
                                bias_correct = TRUE) {
  if (!all(c(cluster_col, member_col) %in% names(df))) return(data.frame())
  if (!all(c("x", "y") %in% names(df))) return(data.frame())
  use_z <- "z" %in% names(df)

  # 1) one median point per MEMBER, carrying its cluster id.
  mem <- df %>%
    dplyr::group_by(.data[[member_col]]) %>%
    dplyr::summarise(
      mx = median(x, na.rm = TRUE),
      my = median(y, na.rm = TRUE),
      mz = if (use_z) median(z, na.rm = TRUE) else 0,
      cl = dplyr::first(.data[[cluster_col]]),
      .groups = "drop"
    )
  # Drop noise: a member id of 0 (only meaningful when members are a cluster
  # level), and any point whose cluster id is 0 / NA.
  if (member_col != "tid")
    mem <- mem[is.finite(mem[[member_col]]) & mem[[member_col]] != 0L, , drop = FALSE]
  mem <- mem[is.finite(mem$cl) & mem$cl != 0L, , drop = FALSE]
  if (!nrow(mem)) return(data.frame())

  # 2) per cluster: centroid = mean of member points; Rg = RMS distance to it.
  cents <- mem %>%
    dplyr::group_by(cl) %>%
    dplyr::summarise(cx = mean(mx), cy = mean(my), cz = mean(mz),
                     n_units = dplyr::n(), .groups = "drop")
  # Bessel factor n/(n−1) on the MEAN SQUARED distance removes the N-dependent
  # small-N bias described in the header (no-op as n grows).
  bc <- isTRUE(bias_correct)
  out <- dplyr::left_join(mem, cents, by = "cl") %>%
    dplyr::mutate(r2 = (mx - cx)^2 + (my - cy)^2 + (mz - cz)^2) %>%
    dplyr::group_by(cl, n_units) %>%
    dplyr::summarise(
      rog = sqrt(mean(r2, na.rm = TRUE) *
                   (if (bc && dplyr::first(n_units) > 1L)
                      dplyr::first(n_units) / (dplyr::first(n_units) - 1) else 1)),
      .groups = "drop") %>%
    dplyr::filter(n_units >= 2L, is.finite(rog)) %>%
    dplyr::rename(cluster_id = cl)

  as.data.frame(out)
}

## ── Ripley's L(r) ────────────────────────────────────────────────────────────
# Computes L(r) - r for a set of 2-D or 3-D points.
# Returns a data.frame: r, Kr, Lr, Lr_minus_r
# K(r) = A * 2*pairs(r) / (n*(n-1)). With no `support`, A is the bounding box and no
# edge correction is applied, so K(r) is underestimated at large r; read it against the
# CSR envelope. With a `support` carrying a mask and edge_correct = TRUE, each pair is
# weighted by the translation correction 1/gammabar (see below) and A is the footprint.
## GOAL: a REPRODUCIBLE subsample for the O(n^2) spatial statistics.
##
## These estimators cap the point count (pairwise distances are quadratic), and the
## cap used a bare `sample.int`, so pressing Run twice on the same data returned
## DIFFERENT curves whenever a file exceeded the cap. That is not acceptable for a
## figure, and it made cached and uncached results compare unequal.
##
## A seeded draw rather than a sorted head: taking the first k rows would be biased
## whenever the table is written in acquisition or spatial order, which MINFLUX
## exports are. The global RNG stream is saved and restored so seeding here cannot
## perturb the CSR envelopes, permutation nulls or thinning draws.
.stable_subsample <- function(n, cap, seed = 20260813L) {
  if (!is.finite(n) || n <= cap) return(seq_len(n))
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  bk  <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", bk, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  sort(sample.int(n, cap))
}

## `support`: optional support object (estimate_support(..., return_mask = TRUE)).
## When given, its measure replaces the bounding-box area in lambda_hat = n / A. See
## the SUPPORT-AWARE OBSERVATION WINDOW banner for why the box is wrong.
compute_ripley_L <- function(x, y, z = NULL, r_seq = NULL, use_3d = FALSE,
                             max_subsample = 2000L, window_ranges = NULL,
                             support = NULL, edge_correct = TRUE) {
  coords <- cbind(x, y)
  if (use_3d && !is.null(z)) coords <- cbind(x, y, z)

  n <- nrow(coords)
  if (n < 3L) return(data.frame(r = numeric(0), Lr_minus_r = numeric(0)))

  # subsample for large datasets
  idx <- .stable_subsample(n, max_subsample)
  if (length(idx) < n) { coords <- coords[idx, , drop = FALSE]; n <- nrow(coords) }

  # observation window (bounding box). A pre-supplied window_ranges (a 2×d matrix
  # of per-axis min/max, as from apply(coords, 2, range)) lets the observed curve
  # and the Monte-Carlo CSR envelope share one identical window so the edge bias
  # cancels between them.
  rng <- if (!is.null(window_ranges)) window_ranges else apply(coords, 2, range)
  # Support measure when supplied, else the bounding box (the historical behaviour).
  area_vol <- if (!is.null(support)) support_area(support, rng) else prod(rng[2, ] - rng[1, ])
  if (!is.finite(area_vol) || area_vol <= 0) {
    return(data.frame(r = numeric(0), Lr_minus_r = numeric(0)))
  }

  if (is.null(r_seq)) {
    # r_max is capped by the SHORTEST window axis (edge bias grows once r approaches
    # half the smallest extent). In 3D the z extent is usually far smaller than x/y
    # (thin optical section), so it MUST be included in the cap.
    r_max <- min(diff(rng[, 1]), diff(rng[, 2])) / 4
    if (use_3d && ncol(rng) > 2) r_max <- min(r_max, diff(rng[, 3]) / 4)
    r_seq <- seq(0, r_max, length.out = 60)[-1]
  }

  # pairwise distances (only upper triangle)
  dm <- as.matrix(stats::dist(coords))

  d_vals <- dm[lower.tri(dm)]

  # Translation edge correction: each PAIR is weighted by 1/gammabar at its own
  # separation, since K(r) = integral over |h| <= r of g(h) dh accumulates pairs of
  # every separation below r. gammabar = area_vol for all pairs recovers the
  # uncorrected estimator exactly.
  # SINGLE-PASS accumulation. The old form swept every pair distance once PER RADIUS
  # (60 radii x ~2e6 pairs), which dominated the cost; binning each distance once and
  # cumulating is one pass. `left.open` makes bin k hold the pairs whose distance is
  # the smallest r_seq value that is >= d, so cumsum(k) counts exactly d <= r_seq[k].
  nb    <- length(r_seq)
  cnt_b <- tabulate(findInterval(d_vals, r_seq, left.open = TRUE) + 1L, nbins = nb)

  # Edge weight applied PER BIN rather than per pair. gamma-bar is already the
  # direction-averaged (isotropised) set covariance, so resolving it finer than the
  # r bins buys nothing - and interpolating it at ~2e6 individual pair distances per
  # simulation was the remaining bottleneck inside the envelope loop.
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d && !is.null(z)) else NULL

  Kr <- if (is.null(gb)) {
    # unbiased K: multiply by 2 for both directions, divide by n^2/area
    (area_vol * 2 * cumsum(cnt_b)) / (n * (n - 1))
  } else {
    w <- ifelse(is.finite(gb) & gb > 0, area_vol / gb, 1)  # 1 = no area lost
    (area_vol * 2 * cumsum(cnt_b * w)) / (n * (n - 1))
  }

  Lr <- if (use_3d && !is.null(z)) {
    (3 * Kr / (4 * pi))^(1/3)
  } else {
    sqrt(Kr / pi)
  }

  data.frame(r = r_seq, Kr = Kr, Lr = Lr, Lr_minus_r = Lr - r_seq)
}

## ── Bivariate Ripley's L12(r) ───────────────────────────────────────────────
compute_ripley_L12 <- function(x1, y1, z1 = NULL, x2, y2, z2 = NULL,
                               r_seq = NULL, use_3d = FALSE,
                               max_subsample = 2000L, window_ranges = NULL,
                               support = NULL, edge_correct = TRUE) {
  if (use_3d && (!is.null(z1)) && (!is.null(z2))) {
    pts1 <- cbind(x1, y1, z1)
    pts2 <- cbind(x2, y2, z2)
  } else {
    pts1 <- cbind(x1, y1)
    pts2 <- cbind(x2, y2)
    use_3d <- FALSE
  }

  n1 <- nrow(pts1); n2 <- nrow(pts2)
  if (n1 < 2L || n2 < 2L) return(data.frame())

  pts1 <- pts1[.stable_subsample(n1, max_subsample), , drop = FALSE]
  pts2 <- pts2[.stable_subsample(n2, max_subsample, seed = 20260814L), , drop = FALSE]
  n1 <- nrow(pts1); n2 <- nrow(pts2)

  all_pts <- rbind(pts1, pts2)
  # per-axis extents; a pre-supplied window_ranges (2×d min/max matrix) makes the
  # observed curve and its CSR envelope share one identical window.
  rng <- if (!is.null(window_ranges)) (window_ranges[2, ] - window_ranges[1, ])
         else apply(all_pts, 2, function(v) diff(range(v, na.rm = TRUE)))
  if (!is.null(support)) {
    area_vol <- support_area(support)
    if (!is.finite(area_vol) || area_vol <= 0) return(data.frame())
  } else if (use_3d) {
    area_vol <- prod(rng)
    if (!is.finite(area_vol) || area_vol == 0) return(data.frame())
  } else {
    area_vol <- prod(rng[1:2])
    if (!is.finite(area_vol) || area_vol == 0) return(data.frame())
  }

  dmat <- as.matrix(dist(rbind(pts1, pts2)))[seq_len(n1), (n1 + 1):(n1 + n2)]

  if (is.null(r_seq) || length(r_seq) == 0) {
    # Cap by the SHORTEST window axis, not the longest: edge bias is set by the
    # dimension that runs out first (in 3D that is almost always z). Matches
    # compute_ripley_L's cap.
    r_max <- min(if (use_3d) rng else rng[1:2]) / 4
    r_seq <- seq(0, r_max, length.out = 60L)[-1]  # drop r=0 (trivially L12(0)-0=0)
  }

  # Translation edge correction, identical in form to the univariate case: divide by
  # the window's isotropised set covariance instead of its area. Applied PER SHELL
  # (gamma-bar is already direction-averaged, so finer resolution buys nothing).
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d) else NULL
  nb <- length(r_seq)
  cnt_b <- tabulate(findInterval(as.vector(dmat), r_seq, left.open = TRUE) + 1L,
                    nbins = nb)
  Kr <- if (is.null(gb)) {
    (area_vol * cumsum(cnt_b)) / (n1 * n2)
  } else {
    w <- ifelse(is.finite(gb) & gb > 0, area_vol / gb, 1)
    (area_vol * cumsum(cnt_b * w)) / (n1 * n2)
  }

  Lr <- if (use_3d) (3 * Kr / (4 * pi))^(1/3) else sqrt(Kr / pi)
  data.frame(r = r_seq, Kr = Kr, Lr = Lr, Lr_minus_r = Lr - r_seq)
}

## ── Pair Cross-Correlation Function g₁₂(r) ──────────────────────────────────
## g12(r) = (area × N_pairs_in_annulus) / (n1 × n2 × annulus_area)
## g12 = 1 under CSR; > 1 = attraction; < 1 = repulsion.
compute_pcf12 <- function(x1, y1, z1 = NULL, x2, y2, z2 = NULL,
                           r_seq = NULL, use_3d = FALSE, max_subsample = 2000L,
                           window_ranges = NULL, support = NULL,
                           edge_correct = TRUE) {
  if (use_3d && !is.null(z1) && !is.null(z2)) {
    pts1 <- cbind(x1, y1, z1); pts2 <- cbind(x2, y2, z2)
  } else {
    pts1 <- cbind(x1, y1);     pts2 <- cbind(x2, y2); use_3d <- FALSE
  }
  n1 <- nrow(pts1); n2 <- nrow(pts2)
  if (n1 < 2L || n2 < 2L) return(data.frame())

  pts1 <- pts1[.stable_subsample(n1, max_subsample), , drop = FALSE]
  pts2 <- pts2[.stable_subsample(n2, max_subsample, seed = 20260814L), , drop = FALSE]
  n1 <- nrow(pts1); n2 <- nrow(pts2)

  all_pts <- rbind(pts1, pts2)
  # per-axis extents; a pre-supplied window_ranges (2×d min/max matrix) makes the
  # observed curve and its CSR envelope share one identical window.
  rng     <- if (!is.null(window_ranges)) (window_ranges[2, ] - window_ranges[1, ])
             else apply(all_pts, 2, function(v) diff(range(v, na.rm = TRUE)))
  area_vol <- if (!is.null(support)) support_area(support)
              else if (use_3d) prod(rng) else prod(rng[1:2])
  if (!is.finite(area_vol) || area_vol <= 0) return(data.frame())

  dmat <- as.matrix(dist(rbind(pts1, pts2)))[seq_len(n1), (n1 + 1):(n1 + n2)]

  if (is.null(r_seq) || !length(r_seq)) {
    # Cap by the SHORTEST window axis (see compute_ripley_L / compute_ripley_L12).
    r_max <- min(if (use_3d) rng else rng[1:2]) / 4
    r_seq <- seq(0, r_max, length.out = 60L)
  }
  dr <- if (length(r_seq) > 1L) mean(diff(r_seq)) else r_seq[1L]

  # See compute_pcf: same translation correction, same single-pass shell binning.
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d) else NULL
  nb    <- length(r_seq)
  edges <- c(r_seq - dr / 2, r_seq[nb] + dr / 2)
  cnt_b <- tabulate(findInterval(as.vector(dmat), edges, left.open = TRUE), nbins = nb)

  g12 <- vapply(seq_along(r_seq), function(k) {
    r <- r_seq[k]
    lo <- max(0, r - dr / 2); hi <- r + dr / 2
    shell <- if (use_3d) (4 / 3) * pi * (hi^3 - lo^3) else pi * (hi^2 - lo^2)
    if (!is.finite(shell) || shell <= 0) return(NA_real_)
    denom <- if (is.null(gb)) area_vol else gb[k]
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    area_vol^2 * cnt_b[k] / (n1 * n2 * denom * shell)
  }, numeric(1L))

  data.frame(r = r_seq, g12 = g12)
}

## ── Univariate Pair-Correlation Function g(r) ────────────────────────────────
## GOAL: the one-channel counterpart of compute_pcf12 - the scale-LOCAL companion
## to Ripley's L(r). Counts ORDERED pairs of the SAME point set falling in the
## annulus/shell (r-dr/2, r+dr/2], normalised by the shell volume and by the
## intensity squared:
##
##     g(r) = A * N_pairs(annulus) / (n(n-1) * shell)      [pairs counted twice]
##
## g = 1 under CSR, > 1 clustering at scale r, < 1 exclusion at scale r.
##
## WHY THIS AND NOT ONLY L(r): K/L are CUMULATIVE, so a short-r label or blinking
## artefact is carried into every larger r and a real cluster shows up as a broad
## rising ramp whose peak position is not the cluster size. g(r) is the derivative
## form, so an aggregate of radius R gives a PEAK near R and a short-r artefact
## stays confined to the first bins. This is the PC-PALM statistic (Sengupta et
## al. 2011; Veatch et al. 2012).
##
## DENSITY NORMALISATION: g is divided by the intensity SQUARED, so - like K/L -
## it is invariant to how densely the sample is labelled, which is what makes it
## the endpoint of choice when comparing cells of differing expression level.
## That invariance is exact only for INDEPENDENT (Poisson) thinning; see the
## LABELLING-DENSITY block below.
##
## Self-pairs (i == j) are excluded via the n(n-1) denominator and by never
## counting the zero diagonal (dr/2 > 0, and the diagonal is dropped explicitly).
## No analytic edge correction - use compute_ripley_envelope(estimator = "g"),
## which carries the identical bias in the same window.
# =============================================================================
# fit_correlation_length() - ONE NUMBER FOR "HOW BIG ARE THE CLUSTERS"
# =============================================================================
#
# GOAL
#   Reduce a g(r) curve to a single comparable number, the correlation length xi,
#   by fitting the excess pair correlation
#
#       g(r) - 1  =  A * exp(-r / xi)
#
#   so that two conditions can be compared with a test on one value per cell
#   instead of by eye across two curve shapes.
#
# WHY AN EXPONENTIAL, AND WHAT xi MEANS
#   g(r) - 1 is the excess probability of finding a second molecule at
#   separation r. For clusters of a characteristic size it decays roughly
#   exponentially, and xi is the distance over which that excess falls by 1/e.
#   Read it as the scale over which molecules "know about" each other - close to
#   a cluster radius for compact clusters, longer for diffuse ones. A is the
#   amplitude at r -> 0, which is NOT a clean measure of clustering strength (see
#   the floor below), so it is reported but should not be the headline.
#
# HOW IT IS FITTED, AND WHY NOT THE OBVIOUS WAY
#   The obvious fit is log-linear: log(g - 1) = log(A) - r/xi is a straight line,
#   so ordinary least squares solves it and cannot fail to converge. It is also
#   BIASED, badly, and the bias is not obvious. Points with g <= 1 have no
#   logarithm and must be dropped - but in a noisy tail those are exactly the
#   downward noise excursions, so dropping them keeps only the upward ones, lifts
#   the tail, flattens the slope and OVERESTIMATES xi. Measured on synthetic
#   curves with a known answer (20 repeats, noise sd 0.05 on g):
#
#       true xi     40      60     120     250
#       log-linear 121      94     115     248      <- unusable below ~100 nm
#       nls         36      56     116     244
#
#   So the fit is nonlinear least squares on g itself, g = 1 + A*exp(-r/xi),
#   which uses every point including those below 1 and needs no logarithm. The
#   log-linear solution is still computed - it seeds nls, and it is the fallback
#   if nls fails to converge, with `method` recording which one produced the
#   number.
#
# THE FLOOR (r_min) IS THE IMPORTANT ARGUMENT
#   Short-r g(r) in SMLM is dominated by LABELLING, not biology: residual blinks
#   of one fluorophore and the several dyes of a multivalent label all pile up
#   within a few tens of nm. Fitting from r = 0 therefore measures the label,
#   not the structure. Start the fit above that floor - a sensible default is a
#   few times the measured trace spread, or the r where a non-clustering control
#   returns to 1.
#
#   Note the floor cuts both ways: if xi turns out to be SMALLER than r_min, the
#   fit began past the decay it was meant to measure and the estimate is not
#   trustworthy. That case is flagged rather than silently returned.
#
# ARGUMENTS
#   r, g   : the curve, as returned by compute_pcf()
#   r_min  : lower bound of the fit (nm). Points below it are ignored.
#   r_max  : upper bound (nm); NULL = use all of the curve above r_min.
#
# RETURNS
#   one-row data.frame: xi, xi_se, amplitude, r2, method, below_floor, n_used,
#   n_avail, r_min, r_max. xi is NA when fewer than 4 usable points remain or the
#   excess does not decay (no decay = no correlation length; reporting one there
#   would be fitting noise).
# =============================================================================
fit_correlation_length <- function(r, g, r_min = 0, r_max = NULL) {
  na_row <- data.frame(xi = NA_real_, xi_se = NA_real_, amplitude = NA_real_,
                       r2 = NA_real_, method = NA_character_, below_floor = FALSE,
                       n_used = 0L, n_avail = 0L,
                       r_min = r_min, r_max = r_max %||% NA_real_,
                       stringsAsFactors = FALSE)
  if (!length(r) || length(r) != length(g)) return(na_row)
  keep <- is.finite(r) & is.finite(g) & r >= r_min
  if (!is.null(r_max) && is.finite(r_max)) keep <- keep & r <= r_max
  na_row$n_avail <- sum(keep)
  if (na_row$n_avail < 4L) return(na_row)

  rr <- r[keep]; gg <- g[keep]
  na_row$n_used <- length(rr)
  na_row$r_max  <- max(rr, na.rm = TRUE)

  # Seed from the log-linear solution, on the excess points only.
  pos <- gg > 1
  if (sum(pos) < 4L) return(na_row)
  lf <- tryCatch(stats::lm(log(gg[pos] - 1) ~ rr[pos]), error = function(e) NULL)
  if (is.null(lf)) return(na_row)
  co <- stats::coef(lf); slope <- unname(co[2])
  if (!is.finite(slope) || slope >= 0) return(na_row)   # no decay to measure
  seed <- list(A = unname(exp(co[1])), xi = -1 / slope)

  fit <- tryCatch(stats::nls(gg ~ 1 + A * exp(-rr / xi), start = seed,
                             control = stats::nls.control(warnOnly = TRUE, maxiter = 200)),
                  error = function(e) NULL)

  if (!is.null(fit) && all(is.finite(stats::coef(fit))) &&
      is.finite(stats::coef(fit)[["xi"]]) && stats::coef(fit)[["xi"]] > 0) {
    cf <- stats::coef(fit)
    se <- tryCatch(summary(fit)$coefficients["xi", "Std. Error"], error = function(e) NA_real_)
    xi <- unname(cf[["xi"]]); amp <- unname(cf[["A"]])
    # Pseudo-R2 against the mean of g: nls has no r.squared, and 1 - RSS/TSS is
    # the comparable quantity.
    rss <- sum(stats::residuals(fit)^2); tss <- sum((gg - mean(gg))^2)
    r2  <- if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_
    meth <- "nls"
  } else {
    # Fallback: the seed itself. Recorded as such - it is the biased estimator,
    # and a reader should know which one produced the number.
    xi  <- seed$xi; amp <- seed$A
    se  <- tryCatch(summary(lf)$coefficients[2, 2] / slope^2, error = function(e) NA_real_)
    r2  <- tryCatch(summary(lf)$r.squared, error = function(e) NA_real_)
    meth <- "log-linear (nls did not converge)"
  }

  data.frame(xi = xi, xi_se = se, amplitude = amp, r2 = r2, method = meth,
             below_floor = is.finite(xi) && xi < r_min,
             n_used = length(rr), n_avail = na_row$n_avail,
             r_min = r_min, r_max = max(rr, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

compute_pcf <- function(x, y, z = NULL, r_seq = NULL, use_3d = FALSE,
                        max_subsample = 2000L, window_ranges = NULL,
                        support = NULL, edge_correct = TRUE) {
  pts <- if (use_3d && !is.null(z)) cbind(x, y, z) else { use_3d <- FALSE; cbind(x, y) }
  n <- nrow(pts)
  if (is.null(n) || n < 3L) return(data.frame())

  idx <- .stable_subsample(n, max_subsample)
  if (length(idx) < n) { pts <- pts[idx, , drop = FALSE]; n <- nrow(pts) }

  # Per-axis extents. A pre-supplied window_ranges (2×d min/max matrix) lets the
  # observed curve and its CSR envelope share one identical window.
  rng      <- if (!is.null(window_ranges)) window_ranges else apply(pts, 2, range)
  ext      <- rng[2, ] - rng[1, ]
  area_vol <- if (!is.null(support)) support_area(support, rng)
              else if (use_3d) prod(ext[1:3]) else prod(ext[1:2])
  if (!is.finite(area_vol) || area_vol <= 0) return(data.frame())

  if (is.null(r_seq) || !length(r_seq)) {
    # Cap by the SHORTEST window axis (matches compute_ripley_L / _L12 / _pcf12):
    # edge bias is set by the dimension that runs out first, in 3D almost always z.
    r_max <- min(if (use_3d) ext[1:3] else ext[1:2]) / 4
    r_seq <- seq(0, r_max, length.out = 60L)
  }
  dr <- if (length(r_seq) > 1L) mean(diff(r_seq)) else r_seq[1L]

  # Upper-triangle distances only: each unordered pair once, so the annulus count
  # is doubled below to get ordered pairs (matching the n(n-1) denominator).
  dm     <- as.matrix(stats::dist(pts))
  d_vals <- dm[lower.tri(dm)]

  # Translation edge correction from the support's set covariance (see
  # set_covariance_iso). Substituting gammabar = area_vol recovers the uncorrected
  # estimator, so `gb = NULL` is exactly the old behaviour.
  gb <- if (!is.null(support) && !is.null(support$mask) && isTRUE(edge_correct))
    .gamma_at(support, r_seq, use_3d = use_3d) else NULL

  # SINGLE-PASS shell counts. The shells (r - dr/2, r + dr/2] tile because r_seq is
  # equally spaced, so one findInterval over the shared edges replaces one full sweep
  # of the pair distances per radius.
  nb    <- length(r_seq)
  edges <- c(r_seq - dr / 2, r_seq[nb] + dr / 2)
  cnt_b <- tabulate(findInterval(d_vals, edges, left.open = TRUE), nbins = nb)

  gr <- vapply(seq_along(r_seq), function(k) {
    r <- r_seq[k]
    lo <- max(0, r - dr / 2); hi <- r + dr / 2
    shell <- if (use_3d) (4 / 3) * pi * (hi^3 - lo^3) else pi * (hi^2 - lo^2)
    if (!is.finite(shell) || shell <= 0) return(NA_real_)
    denom <- if (is.null(gb)) area_vol else gb[k]
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    area_vol^2 * 2 * cnt_b[k] / (n * (n - 1) * denom * shell)
  }, numeric(1L))

  out <- data.frame(r = r_seq, g = gr)
  attr(out, "edge_corrected") <- !is.null(gb)
  out
}

# =============================================================================
# .shift_channel() - ONE RANDOM RIGID SHIFT OF A CHANNEL (the independence null)
# =============================================================================
#
# GOAL
#   Produce one draw of channel B under "A and B are placed independently": B keeps
#   its own spatial structure exactly and only its registration with A is broken.
#
# HOW
#   Translate every B point by the same random vector in x and y, wrapping around
#   inside the analysis box (a toroidal shift; Lotwick & Silverman 1982). z is left
#   alone: MINFLUX volumes are thin slabs, and an in-plane shift already breaks the
#   registration.
#
# TWO RULES THAT MAKE IT WORK ON A REAL CELL
#   1. Minimum shift = the largest radius analysed. A smaller shift leaves B close
#      to where it was, so the "null" still contains part of the association and the
#      test loses power without saying so.
#   2. With a footprint, a shift is accepted only if >= 90% of the shifted points
#      still land inside the cell; the rest are dropped. A shift that throws most of
#      B into empty space is not a plausible arrangement of this cell, and accepting
#      it makes the null anti-conservative (the KDE panel measured this). If no
#      shift within 100 tries passes, the best candidate is used and flagged, since
#      falling back to CSR would destroy B's clustering and bias the null the
#      other way.
#
# ARGUMENTS
#   pts           : list/data.frame with x, y (and z) of channel B
#   window_ranges : 2 x d matrix of the box min/max, the same window as the estimator
#   support       : estimate_support() object with $mask, or NULL
#   min_shift     : smallest accepted toroidal shift, nm
#   use_3d        : carry z through unchanged
#
# RETURNS list(x, y, z) with attr "fallback" TRUE when rule 2 could not be met.
# =============================================================================
.shift_channel <- function(pts, window_ranges, support = NULL, min_shift = 0,
                           use_3d = FALSE, accept_frac = 0.9, max_tries = 100L) {
  lo <- window_ranges[1, ]; wx <- window_ranges[2, 1] - lo[1]; wy <- window_ranges[2, 2] - lo[2]
  x <- pts$x; y <- pts$y; z <- if (use_3d) pts$z else NULL
  # The largest toroidal shift possible is half the box on each axis; never ask for
  # more than half of that, or no shift can qualify in a small box.
  min_shift <- min(min_shift, 0.5 * sqrt((wx / 2)^2 + (wy / 2)^2))
  use_sup <- !is.null(support) && !is.null(support$mask)
  best <- NULL; best_frac <- -1
  for (i in seq_len(max_tries)) {
    sx <- stats::runif(1, 0, wx); sy <- stats::runif(1, 0, wy)
    # toroidal length of the shift: the short way round on each axis
    if (sqrt(min(sx, wx - sx)^2 + min(sy, wy - sy)^2) < min_shift) next
    cx <- lo[1] + ((x - lo[1] + sx) %% wx)
    cy <- lo[2] + ((y - lo[2] + sy) %% wy)
    if (!use_sup) return(list(x = cx, y = cy, z = z))
    inside <- in_support(support, cx, cy, if (use_3d) z else NULL)
    frac <- mean(inside)
    if (frac >= accept_frac)
      return(list(x = cx[inside], y = cy[inside], z = if (use_3d) z[inside] else NULL))
    if (frac > best_frac) { best_frac <- frac; best <- list(cx = cx, cy = cy, inside = inside) }
  }
  if (is.null(best)) {
    # No shift met even the minimum length (tiny box): use an unconstrained one.
    sx <- stats::runif(1, 0, wx); sy <- stats::runif(1, 0, wy)
    best <- list(cx = lo[1] + ((x - lo[1] + sx) %% wx), cy = lo[2] + ((y - lo[2] + sy) %% wy))
    best$inside <- if (use_sup) in_support(support, best$cx, best$cy, if (use_3d) z else NULL)
                   else rep(TRUE, length(x))
  }
  keep <- best$inside
  out <- list(x = best$cx[keep], y = best$cy[keep], z = if (use_3d) z[keep] else NULL)
  attr(out, "fallback") <- TRUE
  out
}

## ── Monte-Carlo CSR envelope for L(r) / L₁₂(r) / g₁₂(r) / g(r) ───────────────
## Generates `nsim` complete-spatial-randomness (CSR) patterns of `n` points
## (and `n2` for the cross statistics) uniformly inside the SAME window as the
## observed pattern (`window_ranges`, a 2×d matrix of per-axis min/max), and runs
## the identical estimator on each. Returns pointwise lo/hi quantiles (default
## 2.5% / 97.5% → a 95% band; a two-sided α≈0.05 test needs nsim ≥ 39) and the
## mean null curve, on the supplied `r_seq`.
##
## Because the sims share the observed window, the (uncorrected) edge bias that
## drags L(r)−r negative at large r is present in the envelope too and cancels
## when the observed curve is read against the band - this envelope is the
## interpretive tool, superseding an explicit edge correction.
## `support`: when supplied, the CSR patterns are drawn uniformly ON THE SUPPORT
## instead of in the bounding box, and the estimators use the support measure. This
## is the mask-aware randomization the 2026-07-30 audit flagged as outstanding: a
## box-CSR null puts simulated points where the sample never was, so it under-states
## how clustered a random pattern INSIDE the cell already looks, and the observed
## curve leaves the band too easily.
## `progress`: optional function(frac) called after each simulation. The envelope is
## the slow part by a wide margin, so a bar that only advances per FILE looks frozen -
## pass a callback and it advances per simulation instead.
compute_ripley_envelope <- function(n, window_ranges, r_seq, use_3d = FALSE,
                                     nsim = 99L, estimator = c("L", "L12", "pcf", "g"),
                                     n2 = NULL, probs = c(0.025, 0.975),
                                     support = NULL, edge_correct = TRUE,
                                     null_mode = c("csr", "labels", "shift"), pts1 = NULL, pts2 = NULL,
                                     progress = NULL) {
  estimator <- match.arg(estimator)
  if (is.null(r_seq) || !length(r_seq)) return(NULL)
  if (is.null(window_ranges) || any(!is.finite(window_ranges))) return(NULL)
  nsim <- max(1L, as.integer(nsim))
  n    <- as.integer(n)
  if (n < 3L) return(NULL)
  # The one-sample estimators ("L", "g") need no second point count; the cross
  # statistics do.
  if (!estimator %in% c("L", "g")) {
    if (is.null(n2)) return(NULL)
    n2 <- as.integer(n2)
    if (n2 < 2L) return(NULL)
  }

  d    <- ncol(window_ranges)
  lo_v <- window_ranges[1, ]; hi_v <- window_ranges[2, ]
  is3d <- isTRUE(use_3d) && d >= 3L

  # Draw on the support when one is given, else uniformly in the box.
  use_sup <- !is.null(support) && !is.null(support$mask)
  gen <- function(m) {
    if (use_sup) {
      s <- sample_in_support(support, m, use_3d = is3d)
      if (!is.null(s)) return(s)
    }
    list(x = stats::runif(m, lo_v[1], hi_v[1]),
         y = stats::runif(m, lo_v[2], hi_v[2]),
         z = if (is3d) stats::runif(m, lo_v[3], hi_v[3]) else NULL)
  }

  nr  <- length(r_seq)
  mat <- matrix(NA_real_, nrow = nr, ncol = nsim)
  cap <- if (estimator %in% c("L", "g")) n else max(n, n2)

  # ── NULL MODE FOR THE TWO-CHANNEL ESTIMATORS ───────────────────────────────
  # "csr"    simulates two independent uniform patterns INSIDE the window. It is
  #          only as good as the window: the observed points fill the real cell and
  #          the simulated ones only its estimate, and the leftover SHAPE error is
  #          not removable by any area correction. Measured on independent CSR
  #          channels in an L-shaped cell (global test, 60 reps): L12 rejected 22%
  #          of the time and g12 8%, against a nominal 5%.
  # "labels" keeps the OBSERVED points exactly where they are and permutes which
  #          channel each one belongs to, so both marginal patterns and the window
  #          cancel out of the comparison. This is the standard random-labelling
  #          null for two types in one pattern. What it tests is therefore "the two
  #          channels are a random labelling of one pattern": if the two species
  #          genuinely occupy different regions of the cell, that difference alone
  #          makes the test reject, which is why it is offered rather than forced.
  #          It CANNOT answer "are A and B placed independently": two species that
  #          share the same clusters ARE a random labelling of one pattern (measured
  #          0% detection), and two unrelated clustered species read as segregated
  #          (100% "exclusion"), because their pooled pattern is more clustered
  #          than either cross pair.
  # "shift"  keeps channel A where it is and translates channel B as a rigid
  #          whole by a random vector, wrapping inside the box (Lotwick & Silverman
  #          1982). Each channel keeps its own clustering exactly; only their
  #          registration is broken. This is the null for "are A and B placed
  #          independently of each other", i.e. the co-localization question.
  #          Measured on clustered channels filling the frame: 5% false positives on
  #          independent species, 100% detection of shared clusters and heterodimers.
  #          See .shift_channel() for the minimum shift and the footprint rule.
  null_mode <- match.arg(null_mode)
  use_labels <- identical(null_mode, "labels") && !is.null(pts1) && !is.null(pts2) &&
                estimator %in% c("L12", "pcf")
  use_shift  <- identical(null_mode, "shift") && !is.null(pts1) && !is.null(pts2) &&
                estimator %in% c("L12", "pcf")
  if (use_labels) {
    pool <- rbind(data.frame(x = pts1$x, y = pts1$y, z = if (!is.null(pts1$z)) pts1$z else 0),
                  data.frame(x = pts2$x, y = pts2$y, z = if (!is.null(pts2$z)) pts2$z else 0))
    n1_lab <- nrow(pts1); n_pool <- nrow(pool)
  }
  shift_fallbacks <- 0L
  if (use_shift) {
    # Channel A stays fixed (as observed); channel B is re-drawn per simulation.
    fixed1 <- list(x = pts1$x, y = pts1$y, z = if (is3d) pts1$z else NULL)
    min_shift <- max(r_seq, na.rm = TRUE)
  }

  # Set covariance ONCE for the whole loop - it is a property of the window, and
  # rebuilding it per simulation is what made this look like a hang.
  if (!is.null(support) && isTRUE(edge_correct)) {
    d_max <- sqrt(sum((hi_v - lo_v)^2))
    support <- attach_gamma_cache(support, r_max = max(max(r_seq), d_max),
                                 use_3d = is3d)
  }

  # WHY THE NULL IS NOT RE-ESTIMATED PER SIMULATION (tried and rejected 2026-09-17).
  # The observed points occupy the CELL, while a simulated pattern is confined to the
  # estimated footprint, which is smaller (see debias_support_window). Re-estimating
  # each simulation's area from its own points makes the mismatch worse, not better:
  # random data then read ABOVE the band in 100% of replicates instead of below it.
  # The region, not the normalising area, is what has to match, so the fix lives in
  # the window (debias_support_window) and this loop uses the window it is given.
  for (s in seq_len(nsim)) {
    p1 <- gen(n)
    sup_s <- support
    col <- if (estimator == "L") {
      cv <- compute_ripley_L(p1$x, p1$y, p1$z, r_seq = r_seq, use_3d = is3d,
                             max_subsample = cap, window_ranges = window_ranges,
                             support = sup_s, edge_correct = edge_correct)
      cv$Lr_minus_r
    } else if (estimator == "g") {
      cv <- compute_pcf(p1$x, p1$y, p1$z, r_seq = r_seq, use_3d = is3d,
                        max_subsample = cap, window_ranges = window_ranges,
                        support = sup_s, edge_correct = edge_correct)
      cv$g
    } else if (estimator == "L12") {
      if (use_labels) {
        i1 <- sample.int(n_pool, n1_lab)
        p1 <- as.list(pool[i1, ]); p2 <- as.list(pool[-i1, ])
        if (!is3d) { p1$z <- NULL; p2$z <- NULL }
      } else if (use_shift) {
        p1 <- fixed1
        p2 <- .shift_channel(pts2, window_ranges, support, min_shift, is3d)
        shift_fallbacks <- shift_fallbacks + isTRUE(attr(p2, "fallback"))
      } else p2 <- gen(n2)
      # edge_correct passed through (fixed 2026-09-17): the default TRUE here made
      # the envelope edge-corrected even when the observed curve was not.
      cv <- compute_ripley_L12(p1$x, p1$y, p1$z, p2$x, p2$y, p2$z,
                               r_seq = r_seq, use_3d = is3d,
                               max_subsample = cap, window_ranges = window_ranges,
                               support = sup_s, edge_correct = edge_correct)
      cv$Lr_minus_r
    } else {
      if (use_labels) {
        i1 <- sample.int(n_pool, n1_lab)
        p1 <- as.list(pool[i1, ]); p2 <- as.list(pool[-i1, ])
        if (!is3d) { p1$z <- NULL; p2$z <- NULL }
      } else if (use_shift) {
        p1 <- fixed1
        p2 <- .shift_channel(pts2, window_ranges, support, min_shift, is3d)
        shift_fallbacks <- shift_fallbacks + isTRUE(attr(p2, "fallback"))
      } else p2 <- gen(n2)
      cv <- compute_pcf12(p1$x, p1$y, p1$z, p2$x, p2$y, p2$z,
                          r_seq = r_seq, use_3d = is3d,
                          max_subsample = cap, window_ranges = window_ranges,
                          support = sup_s, edge_correct = edge_correct)
      cv$g12
    }
    if (length(col) == nr) {
      mat[, s] <- col
    } else if (length(col) > 0L) {
      k <- min(length(col), nr)
      mat[seq_len(k), s] <- col[seq_len(k)]
    }
    if (is.function(progress)) try(progress(1 / nsim), silent = TRUE)
  }

  lo <- apply(mat, 1, stats::quantile, probs = probs[1], na.rm = TRUE)
  hi <- apply(mat, 1, stats::quantile, probs = probs[2], na.rm = TRUE)

  # ── GLOBAL ENVELOPE (added 2026-09-17) ──────────────────────────────────────
  # The band above is POINTWISE: it holds at each r separately, so reading
  # "the curve leaves the band ANYWHERE" over 40 radii is not a 5% test. Measured
  # on CSR in an L-shaped cell with the true window, random data left the pointwise
  # band somewhere in 73% of replicates. The global band is the studentized
  # maximum-deviation (MAD) test: for each simulated curve
  #     T = max_r |f(r) - mean(r)| / sd(r)
  # and the band is mean +/- T_crit * sd with T_crit the (1 - alpha) quantile of the
  # simulated T. Crossing it ANYWHERE is then a single test at level alpha
  # (Loosmore & Ford 2006; Myllymaki et al. 2017). mean and sd for each simulated
  # curve exclude that curve (leave-one-out), so its T is on the same footing as the
  # observed one, which was not part of the simulations at all.
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, stats::sd, na.rm = TRUE)
  alpha <- max(1e-3, min(0.5, 2 * probs[1]))
  ok_r <- is.finite(mu) & is.finite(sdv) & sdv > 0
  T_sims <- rep(NA_real_, nsim)
  if (any(ok_r) && nsim >= 3L) {
    m_ok <- mat[ok_r, , drop = FALSE]
    n_j  <- rowSums(!is.na(m_ok))
    s1   <- rowSums(m_ok, na.rm = TRUE)
    s2   <- rowSums(m_ok^2, na.rm = TRUE)
    for (j in seq_len(nsim)) {
      v  <- m_ok[, j]
      nj <- n_j - !is.na(v)
      mj <- (s1 - ifelse(is.na(v), 0, v)) / pmax(nj, 1)
      vj <- (s2 - ifelse(is.na(v), 0, v^2) - nj * mj^2) / pmax(nj - 1, 1)
      sj <- sqrt(pmax(vj, 0))
      d  <- abs(v - mj) / sj
      T_sims[j] <- suppressWarnings(max(d[is.finite(d)], na.rm = TRUE))
    }
  }
  T_sims <- T_sims[is.finite(T_sims)]
  T_crit <- if (length(T_sims) >= 3L) unname(stats::quantile(T_sims, 1 - alpha)) else NA_real_
  out <- data.frame(r = r_seq, lo = as.numeric(lo), hi = as.numeric(hi),
                    null_mean = mu, null_sd = sdv,
                    g_lo = mu - T_crit * sdv, g_hi = mu + T_crit * sdv)
  attr(out, "null_mode") <- if (use_labels) "labels" else if (use_shift) "shift" else "csr"
  # Shifts that could not meet the footprint rule and used the best candidate instead.
  attr(out, "shift_fallbacks") <- if (use_shift) shift_fallbacks else NA_integer_
  attr(out, "T_sims") <- T_sims
  attr(out, "T_crit") <- T_crit
  attr(out, "alpha")  <- alpha
  attr(out, "nsim")   <- nsim
  out
}

# =============================================================================
# envelope_global_p() - one p-value for a curve against its CSR envelope
# =============================================================================
# GOAL: turn "does the curve leave the band" into a single test over the whole r
# range, using the studentized maximum deviation stored on the envelope.
#
# env : a compute_ripley_envelope() result (needs null_mean, null_sd, T_sims)
# obs : the observed curve, aligned with env$r
#
# RETURNS list(p, T_obs, T_crit, alpha, nsim, r_at_max) where p is the Monte-Carlo
# p-value (1 + #{T_sim >= T_obs}) / (nsim + 1), i.e. never below 1/(nsim+1), and
# r_at_max is the radius where the largest standardized deviation sits.
envelope_global_p <- function(env, obs) {
  if (is.null(env) || is.null(obs) || !all(c("null_mean", "null_sd") %in% names(env)))
    return(NULL)
  T_sims <- attr(env, "T_sims")
  n_ok <- min(length(obs), nrow(env))
  d <- abs(obs[seq_len(n_ok)] - env$null_mean[seq_len(n_ok)]) / env$null_sd[seq_len(n_ok)]
  d[!is.finite(d)] <- NA_real_
  if (all(is.na(d))) return(NULL)
  T_obs <- max(d, na.rm = TRUE)
  p <- if (length(T_sims) >= 3L)
         (1 + sum(T_sims >= T_obs, na.rm = TRUE)) / (length(T_sims) + 1) else NA_real_
  list(p = p, T_obs = T_obs, T_crit = attr(env, "T_crit"), alpha = attr(env, "alpha"),
       nsim = attr(env, "nsim"), r_at_max = env$r[which.max(replace(d, is.na(d), -Inf))])
}

## ── 2D KDE density co-localization map ──────────────────────────────────────
## Returns KDE grids for Ch1 and Ch2 on a shared spatial grid,
## plus Pearson r of the two density surfaces.
## Optional honesty add-ons:
##   bw_sweep    – vector of bandwidths; returns $bw_sweep = data.frame(bw,
##                 pearson_r) so the bandwidth-dependence of Pearson r is visible.
##   costes_nsim – >0 builds a randomization null for Pearson r this many times and
##                 returns a one-sided p-value (fraction of null r ≥ observed) plus
##                 the null mean/SD.
##   null_mode   – WHICH null. Measured false-positive rates are in the note below;
##                 read it before quoting any p-value from this function.
##     "shift" (default) – TOROIDAL RANDOM SHIFT of channel 2 (Lotwick & Silverman
##                 1982). Channel 2 is translated by a random offset with wrap-around,
##                 so its OWN clustering, density and point count are preserved exactly
##                 and only its spatial REGISTRATION with channel 1 is destroyed. This
##                 is the null that isolates cross-channel dependence, and it is the
##                 same random-shift logic already used for the bivariate statistics.
##     "csr"     – channel 2 replaced by uniform random points, i.e. a test against
##                 complete spatial randomness rather than against independent
##                 placement. It is NOT the Costes et al. 2004 test (theirs preserves
##                 the randomized channel's autocorrelation by block randomization).
##
## ⚠ CALIBRATION - measured, not assumed (2026-07-30). Simulated 30× per scenario,
## independent placement (truth = NOT co-localized), so a calibrated test should
## reject at p<0.05 about 5% of the time:
##   • Channels clustered but spread over the WHOLE window:
##         shift 3% · csr 3%   → both fine.
##   • Channels sharing ONE CELL that fills a small part of the window (the realistic
##     SMLM case), placed independently within it:
##         shift 53% · csr 93% → BOTH anti-conservative; csr severely so.
## So "csr" really does mostly report "both channels are non-uniform", but the
## random-shift default is only a PARTIAL fix: wrapping Ch2 moves it outside the cell,
## so the null still contrasts "Ch2 inside the cell" with "Ch2 partly outside it". The
## root cause is that BOTH nulls randomize over the bounding box rather than over the
## actual support (the cell mask). Until randomization is mask-aware, treat a low
## p-value here as WEAK evidence whenever the signal does not fill the field, and get
## significance from g12(r) against its CSR envelope instead. Power is not the problem:
## on a genuinely co-localized positive control the shift null fired 10/10, median
## p = 0.02.
## `support`: when supplied, BOTH nulls randomize over the sample support instead of
## the bounding box. This is the fix the 2026-07-30 audit named as outstanding - it
## measured 53% (shift) and 93% (CSR) false-positive rates when both channels shared
## one cell inside a larger window, because a box-randomized Ch2 gets scattered into
## regions where no sample exists, decorrelating it far more than a real
## rearrangement inside the cell would. See the measured post-fix rates in
## scratchpad/test_support_window.R and the README.
compute_kde_coloc <- function(x1, y1, x2, y2, bandwidth = NULL, n_grid = 100L,
                              bw_sweep = NULL, costes_nsim = 0L,
                              null_mode = c("shift", "csr"), support = NULL) {
  null_mode <- match.arg(null_mode)
  all_x <- c(x1, x2); all_y <- c(y1, y2)
  xlim  <- range(all_x, na.rm = TRUE); ylim <- range(all_y, na.rm = TRUE)
  if (diff(xlim) == 0 || diff(ylim) == 0) return(NULL)

  if (is.null(bandwidth)) {
    n   <- length(all_x)
    bw  <- 0.9 * min(stats::sd(all_x, na.rm = TRUE),
                     stats::IQR(all_x, na.rm = TRUE) / 1.34,
                     na.rm = TRUE) * n^(-0.2)
    bandwidth <- if (is.finite(bw) && bw > 0) bw else diff(xlim) / 20
  }

  k1 <- MASS::kde2d(x1, y1, h = bandwidth, n = n_grid, lims = c(xlim, ylim))
  k2 <- MASS::kde2d(x2, y2, h = bandwidth, n = n_grid, lims = c(xlim, ylim))

  pearson_r <- stats::cor(as.vector(k1$z), as.vector(k2$z))

  out <- list(x = k1$x, y = k1$y, z1 = k1$z, z2 = k2$z,
              pearson_r = pearson_r, bandwidth = bandwidth)

  # ── Bandwidth sweep: r(bw) ─────────────────────────────────────────────────
  if (!is.null(bw_sweep) && length(bw_sweep)) {
    bw_sweep <- sort(unique(bw_sweep[is.finite(bw_sweep) & bw_sweep > 0]))
    if (length(bw_sweep)) {
      sweep_r <- vapply(bw_sweep, function(h) {
        ka <- MASS::kde2d(x1, y1, h = h, n = n_grid, lims = c(xlim, ylim))
        kb <- MASS::kde2d(x2, y2, h = h, n = n_grid, lims = c(xlim, ylim))
        stats::cor(as.vector(ka$z), as.vector(kb$z))
      }, numeric(1L))
      out$bw_sweep <- data.frame(bw = bw_sweep, pearson_r = sweep_r)
    }
  }

  # ── Spatial randomization test (see null_mode in the header) ──────────────
  costes_nsim <- as.integer(costes_nsim %||% 0L)
  if (isTRUE(costes_nsim > 0L)) {
    n2 <- length(x2)
    wx <- diff(xlim); wy <- diff(ylim)
    use_sup <- !is.null(support) && !is.null(support$mask)
    null_r <- vapply(seq_len(costes_nsim), function(i) {
      if (identical(null_mode, "shift")) {
        # Toroidal shift: wrap channel 2 inside the SAME window. Preserves its own
        # clustering and n exactly; only the cross-channel registration is broken.
        # With a support, shifts are REJECTION-SAMPLED so the shifted channel still
        # lands inside the cell: a shift that throws most of Ch2 into empty space is
        # not a plausible rearrangement of the sample, and accepting it is exactly
        # what made this null anti-conservative.
        accept <- function(rx, ry) mean(in_support(support, rx, ry)) >= 0.9
        rx <- ry <- NULL
        for (tries in seq_len(if (use_sup) 60L else 1L)) {
          sx <- stats::runif(1, 0, wx); sy <- stats::runif(1, 0, wy)
          cx <- xlim[1] + ((x2 - xlim[1] + sx) %% wx)
          cy <- ylim[1] + ((y2 - ylim[1] + sy) %% wy)
          if (!use_sup || accept(cx, cy)) { rx <- cx; ry <- cy; break }
        }
        if (is.null(rx)) {
          # No acceptable shift found (a cell filling most of its box leaves little
          # room). Fall back to a support-restricted CSR for this replicate rather
          # than silently using an out-of-cell shift.
          s <- sample_in_support(support, n2)
          rx <- s$x; ry <- s$y
        }
      } else if (use_sup) {
        s <- sample_in_support(support, n2)
        rx <- s$x; ry <- s$y
      } else {
        rx <- stats::runif(n2, xlim[1], xlim[2])
        ry <- stats::runif(n2, ylim[1], ylim[2])
      }
      kb <- MASS::kde2d(rx, ry, h = bandwidth, n = n_grid, lims = c(xlim, ylim))
      stats::cor(as.vector(k1$z), as.vector(kb$z))
    }, numeric(1L))
    # one-sided p: how often randomized co-localization matches/exceeds observed
    out$costes_p         <- (sum(null_r >= pearson_r, na.rm = TRUE) + 1) /
                            (sum(is.finite(null_r)) + 1)
    out$costes_null_mean <- mean(null_r, na.rm = TRUE)
    out$costes_null_sd   <- stats::sd(null_r, na.rm = TRUE)
    out$costes_nsim      <- costes_nsim
    out$null_mode        <- null_mode
    out$null_support     <- use_sup
    out$null_label       <- paste0(
      if (identical(null_mode, "shift"))
        "random toroidal shift of Ch2 (structure-preserving; Lotwick & Silverman 1982)"
      else "uniform CSR resampling of Ch2 (tests against complete spatial randomness)",
      if (use_sup) ", randomized over the sample support" else ", randomized over the bounding box")
  }

  out
}


# =============================================================================
# spatial_curve_for_points() - ONE FILE'S L(r) OR g(r), WITH ITS ENVELOPE
# =============================================================================
#
# GOAL
#   The computation behind the Ripley's L and g(r) panels for one point set: the
#   curve, and optionally its randomisation envelope with the global
#   (studentized maximum-deviation) test. The panels call this once per file or
#   channel; resolving the footprint, the progress bar and storing results stay
#   in the server.
#
# THE WINDOW AND THE FOOTPRINT ARE INPUTS
#   `win` and `support` are passed in, not derived from `pts`: when a file is
#   split into channels, both channels must be measured in the pooled file's
#   window, or the sparser channel's window shrinks around it and its curve is
#   deflated (see split_df_by_channel()). The curve and its envelope always share
#   the one window given here, so the edge bias cancels between them.
#
# ARGUMENTS
#   pts          data.frame x, y (and z), one row per unit (unit_centroids())
#   win          2 x d matrix of the window's per-axis min/max
#   estimator    "L" (Ripley's L, reported as L(r) - r) or "g" (pair correlation)
#   r_seq        radii (nm); NULL lets the estimator choose from the window
#   use_3d       measure in 3D
#   support      footprint (estimate_support(..., return_mask = TRUE)) or NULL
#   edge_correct translation edge correction from the footprint's set covariance
#   envelope     also compute the envelope and global p-value
#   nsim         simulations for the envelope
#   progress     function(frac) called per simulation, or NULL
#   tag          dataset label written into the result rows
#
# RETURNS  list(curve, env) - env NULL without an envelope - or NULL when the
#          estimator returns no curve.
# =============================================================================
spatial_curve_for_points <- function(pts, win, estimator = c("L", "g"), r_seq = NULL,
                                     use_3d = FALSE, support = NULL, edge_correct = TRUE,
                                     envelope = FALSE, nsim = 99L, progress = NULL,
                                     tag = NULL) {
  estimator <- match.arg(estimator)
  # One set-covariance build per file, reused by the curve AND every simulation.
  if (!is.null(support) && isTRUE(edge_correct))
    support <- attach_gamma_cache(support, r_max = sqrt(sum((win[2, ] - win[1, ])^2)),
                                  use_3d = use_3d)
  z <- if (use_3d) pts$z else NULL
  curve <- if (estimator == "L")
    compute_ripley_L(pts$x, pts$y, z, r_seq = r_seq, use_3d = use_3d, window_ranges = win,
                     support = support, edge_correct = isTRUE(edge_correct))
  else
    compute_pcf(pts$x, pts$y, z, r_seq = r_seq, use_3d = use_3d, window_ranges = win,
                support = support, edge_correct = isTRUE(edge_correct))
  if (is.null(curve) || !nrow(curve)) return(NULL)
  curve$dataset <- tag
  if (estimator == "g") curve$n_units <- nrow(pts)
  env <- NULL
  if (isTRUE(envelope)) {
    env <- compute_ripley_envelope(n = min(nrow(pts), 2000L), window_ranges = win,
                                   r_seq = curve$r, use_3d = use_3d, nsim = nsim,
                                   estimator = estimator, support = support,
                                   edge_correct = isTRUE(edge_correct), progress = progress)
    # Global test, computed here because rbind drops the attributes the envelope
    # carries. One p per file for the whole r range.
    if (!is.null(env)) {
      gp <- envelope_global_p(env, if (estimator == "L") curve$Lr_minus_r else curve$g)
      env$global_p <- gp$p %||% NA_real_
      env$dataset <- tag
    }
  }
  list(curve = curve, env = env)
}

# =============================================================================
# PANEL HELPERS (moved from server_smlm.R; plain functions of their arguments)
# =============================================================================

## The shaded envelope is POINTWISE: it holds at each r on its own, so a curve
## that leaves it somewhere is not significant at 5% - on simulated random data
## in a cell footprint that happens in 96-100% of runs. `g_lo`/`g_hi` (from
## compute_ripley_envelope) are the studentized maximum-deviation band, which the
## curve crosses anywhere with probability 5% under the null, and `global_p` is
## the matching p-value for the whole r range. Both are drawn/reported per file.
## `smooth_fn` is the builder's own smoother (each panel defines its own from its
## smoothing controls), so the band gets exactly the same treatment as the curve.
.add_global_band <- function(fig, e, col, nm, smooth_fn = NULL) {
  if (is.null(e$g_lo) || all(!is.finite(e$g_lo))) return(fig)
  sm <- if (is.function(smooth_fn)) smooth_fn else function(r, y) y
  gl <- sm(e$r, e$g_lo); gh <- sm(e$r, e$g_hi)
  fig %>%
    plotly::add_trace(x = e$r, y = gh, type = "scatter", mode = "lines",
      line = list(color = col, width = 1, dash = "dash"), legendgroup = nm,
      name = paste0(nm, " global 95% band"), showlegend = TRUE, hoverinfo = "skip") %>%
    plotly::add_trace(x = e$r, y = gl, type = "scatter", mode = "lines",
      line = list(color = col, width = 1, dash = "dash"), legendgroup = nm,
      showlegend = FALSE, hoverinfo = "skip")
}

## One line naming each file's global p, for the plot caption. Files beyond the
## fourth are summarised so the caption cannot grow without bound.
.global_p_caption <- function(env_res) {
  if (is.null(env_res) || !"global_p" %in% names(env_res)) return(NULL)
  d <- unique(env_res[, c("dataset", "global_p")])
  d <- d[is.finite(d$global_p), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  fmt <- function(p) if (p <= 1 / 100) sprintf("p \u2264 %.2g", p) else sprintf("p = %.2f", p)
  txt <- paste(sprintf("%s %s", d$dataset, vapply(d$global_p, fmt, character(1))), collapse = "  \u00b7  ")
  if (nrow(d) > 4) txt <- paste0(paste(sprintf("%s %s", d$dataset[1:4], vapply(d$global_p[1:4], fmt, character(1))),
                                       collapse = "  \u00b7  "),
                                 sprintf("  \u00b7  (+%d more, %d with p < 0.05)", nrow(d) - 4,
                                         sum(d$global_p < 0.05, na.rm = TRUE)))
  # Two-channel panels: say which null the p-value is against, since the three
  # nulls answer different questions (see compute_ripley_envelope).
  nl <- if ("null_mode" %in% names(env_res)) unique(stats::na.omit(env_res$null_mode)) else character(0)
  nl_txt <- c(shift = "null: channel B shifted at random (independent placement)",
              labels = "null: random labelling (do the channels differ in arrangement)",
              csr = "null: independent random patterns in the window")[nl]
  nl_txt <- nl_txt[!is.na(nl_txt)]
  paste0("Global test over the whole r range: ", txt,
         if (length(nl_txt)) paste0("  |  ", paste(nl_txt, collapse = "; ")) else "")
}

## Caption as a plot annotation (NULL when no envelope was computed).
.global_p_annot <- function(env_res) {
  cap <- .global_p_caption(env_res)
  if (is.null(cap)) return(NULL)
  list(list(text = cap, xref = "paper", yref = "paper", x = 0, y = 1.06,
            xanchor = "left", showarrow = FALSE, font = list(size = 11, color = "#555")))
}

## GOAL: say on the axis why a group curve stops where it does.
##
## In group mode the curves are averaged on a grid capped at the SMALLEST
## per-file r max, so if one cell was measured to 400 nm and another to 1500 the
## group curve ends at 400. That is deliberate - averaging past it would compare
## a mean of five files against a mean of one - but with the axis now ending
## there too, nothing on screen says the longer files were truncated. This puts
## it in the axis title, where it costs no space and cannot be missed.
.rmax_axis_title <- function(grouped, res, rmax_c, show_files) {
  base <- "r (nm)"
  if (!isTRUE(grouped) || isTRUE(show_files)) return(base)
  rmax_all <- suppressWarnings(max(res$r, na.rm = TRUE))
  if (!is.finite(rmax_c) || !is.finite(rmax_all) || rmax_all <= rmax_c * 1.02)
    return(base)
  sprintf("r (nm) - group curves stop at %.0f nm, the smallest r max in the set (longest file reaches %.0f)",
          rmax_c, rmax_all)
}

