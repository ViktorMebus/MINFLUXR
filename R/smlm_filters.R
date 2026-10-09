## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_filters.R  -  Filters, trace spread and localization precision
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The metadata / outlier / length / spread filters applied before
##   clustering, the per-trace spread statistics and their plot overlays,
##   rogue-trace removal and the localization-precision estimate.
##
## CONTENTS (in file order)
##   compute_sd_medians, trace_spread_table, sd_breaks, sd_counts, grp_median,
##   grp_mad, .MAD_SMALL_N_FACTOR, mad_finite_sample_factor, sd_by_tid,
##   trace_point_outliers, trim_summary, cov_ellipse_params, ellipse_xy,
##   density_hdr_2d, robust_median_cov, sigma_label_med_pm_sd,
##   add_sigma_iqr_lines, add_sigma_vline, add_median_vlines,
##   apply_meta_filters_batch, remove_rogue_traces, geometric_median,
##   smlm_loc_precision_nm
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

# =============================================================================
# compute_sd_medians() - THE REPORTED SPREAD SUMMARY
# =============================================================================
#
# GOAL
#   Collapse the per-trace spreads into the numbers the UI reports: a centre, a
#   robust width, and the interquartile range, per axis.
#
# REPORT vs DECIDE - WHY robust DEFAULTS TO FALSE HERE
#   The per-trace sigma fed in is the plain n-1 SD (see sd_by_tid), because that
#   is the one number with an outside reference. Robustness enters ACROSS traces
#   (median and MAD), which is where pyMINFLUX puts it too, and where it costs
#   no comparability. The corrected-MAD branch stays reachable for the internal
#   callers that need a breakdown-resistant per-trace sigma.
#
# RETURNS
#   list(x, y[, z]) medians as before, plus mad_x/mad_y[/mad_z] and
#   q1_x/q3_x/... quartiles. NA where an axis cannot be measured.
# =============================================================================
compute_sd_medians <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, robust = FALSE) {
  tbl <- trace_spread_table(df, min_n = min_n, robust = robust)
  has_z <- "sd_z" %in% names(tbl) && any(is.finite(tbl$sd_z))
  if (!nrow(tbl)) {
    # Consistent empty result so the callers don't crash
    return(list(x = NA_real_, y = NA_real_, z = if ("sd_z" %in% names(tbl)) NA_real_ else NULL,
                mad_x = NA_real_, mad_y = NA_real_, mad_z = NA_real_,
                q1_x = NA_real_, q3_x = NA_real_, q1_y = NA_real_, q3_y = NA_real_,
                q1_z = NA_real_, q3_z = NA_real_, n_traces = 0L))
  }

  # One axis -> centre, robust width, quartiles. Every statistic drops NA first
  # so a trace that is too short on one axis does not blank the whole summary.
  ax_stats <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(list(med = NA_real_, mad = NA_real_, q1 = NA_real_, q3 = NA_real_))
    q <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    list(med = stats::median(v),
         # constant = 1.4826 matches scipy's scale = 0.67449 used by pyMINFLUX
         mad = stats::mad(v, constant = 1.4826, na.rm = TRUE),
         q1  = q[1], q3 = q[2])
  }

  sx <- ax_stats(tbl$sd_x); sy <- ax_stats(tbl$sd_y)
  sz <- if (has_z) ax_stats(tbl$sd_z) else NULL

  res <- list(
    x = sx$med, y = sy$med,
    mad_x = sx$mad, mad_y = sy$mad,
    q1_x = sx$q1, q3_x = sx$q3,
    q1_y = sy$q1, q3_y = sy$q3,
    n_traces = nrow(tbl)
  )
  if (!is.null(sz)) {
    res$z <- sz$med; res$mad_z <- sz$mad; res$q1_z <- sz$q1; res$q3_z <- sz$q3
  }
  res
}

trace_spread_table <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, robust = FALSE) {
  stopifnot(is.data.frame(df))
  if (!nrow(df)) {
    return(data.frame(
      tid = character(),
      sd_x = numeric(), sd_y = numeric(), sd_z = numeric(),
      spread_rmsd = numeric(), sd_max = numeric()
    ))
  }

  min_n <- max(2L, as.integer(min_n %||% MIN_LOCS_FOR_SPREAD_DEFAULT))

  has_z <- "z" %in% names(df) && any(is.finite(df$z))

  if (isTRUE(robust)) {
    # Robust branch: per-axis spread comes from sd_by_tid (corrected MAD), and
    # the overall spread is rebuilt from those axes rather than measured as an
    # RMSD. RMSD squares the residuals, so one stray localization would dominate
    # spread_rmsd even though every sd_ column had shrugged it off - the table
    # would then contradict itself. For Gaussian residuals sqrt(sum of per-axis
    # sigma^2) has the same expectation as the RMSD it replaces (sqrt(d)*sigma),
    # so the column keeps its meaning and its scale.
    out <- sd_by_tid(df, min_n = min_n, robust = TRUE)
    if (!nrow(out)) return(data.frame(
      tid = character(), sd_x = numeric(), sd_y = numeric(), sd_z = numeric(),
      spread_rmsd = numeric(), sd_max = numeric()))
    if (!has_z) out$sd_z <- NA_real_
    sq <- out$sd_x^2 + out$sd_y^2 + if (has_z) out$sd_z^2 else 0
    out$spread_rmsd <- sqrt(sq)
    return(out %>%
      dplyr::filter(.data$n_locs >= min_n) %>%
      dplyr::mutate(sd_max = pmax(sd_x, sd_y, sd_z, na.rm = TRUE)) %>%
      dplyr::select(-dplyr::any_of("n_locs")))
  }

  base <- df %>% dplyr::group_by(tid)

  if (has_z) {
    out <- base %>% dplyr::summarise(
      n_locs = dplyr::n(),
      mx = stats::median(x, na.rm = TRUE),
      my = stats::median(y, na.rm = TRUE),
      mz = stats::median(z, na.rm = TRUE),
      sd_x = stats::sd(x, na.rm = TRUE),
      sd_y = stats::sd(y, na.rm = TRUE),
      sd_z = stats::sd(z, na.rm = TRUE),
      spread_rmsd = sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  } else {
    out <- base %>% dplyr::summarise(
      n_locs = dplyr::n(),
      mx = stats::median(x, na.rm = TRUE),
      my = stats::median(y, na.rm = TRUE),
      sd_x = stats::sd(x, na.rm = TRUE),
      sd_y = stats::sd(y, na.rm = TRUE),
      sd_z = NA_real_,
      spread_rmsd = sqrt(mean((x - mx)^2 + (y - my)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  }
  out %>%
    dplyr::filter(.data$n_locs >= min_n) %>%
    dplyr::mutate(sd_max = pmax(sd_x, sd_y, sd_z, na.rm = TRUE)) %>%
    dplyr::select(-dplyr::any_of(c("mx", "my", "mz", "n_locs")))
}

# Generate histogram breaks for a numeric vector over given number of bins
sd_breaks <- function(x, bins) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    rng <- c(rng[1] - 0.5, rng[2] + 0.5)
  }
  seq(rng[1], rng[2], length.out = max(2, as.integer(bins) + 1L))
}

# Turn a numeric vector into (centers, counts) on given breaks
sd_counts <- function(x, breaks) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (is.null(breaks)) return(list(centers = numeric(0), counts = numeric(0), width = 1))
  
  bw  <- diff(breaks)
  eps <- if (length(bw)) max(1e-9, 1e-6 * max(bw, na.rm = TRUE)) else 1e-6
  brk <- breaks
  brk[1] <- brk[1] - eps
  brk[length(brk)] <- brk[length(brk)] + eps
  
  if (length(x)) {
    lo <- brk[1]; hi <- brk[length(brk)]
    x  <- pmin(hi, pmax(lo, x))
  }
  
  h <- hist(x, breaks = brk, plot = FALSE, include.lowest = TRUE, right = TRUE)
  
  w <- suppressWarnings(as.numeric(diff(h$breaks)[1]))
  if (!is.finite(w) || w <= 0) w <- 1
  
  list(centers = h$mids, counts = h$counts, width = w)
}

# =============================================================================
# grp_median() / grp_mad() - PER-GROUP MEDIAN AND MAD IN ONE RADIX SORT
# =============================================================================
#
# GOAL
#   Per-trace median and MAD fast enough to sit in a reactive that refires on
#   every filter change.
#
# WHY NOT stats::mad() UNDER group_by()
#   Measured on 600k localizations / 60k traces, three axes:
#     dplyr + stats::sd    1.5 s
#     dplyr + stats::mad   8.6 s     <- 5.9x, unusable in a reactive
#     this                 0.7 s     <- 0.5x, i.e. FASTER than the sd it replaces
#   stats::mad is a generic that re-sorts for each of its two medians and pays
#   full argument dispatch per group. Sorting ONCE by (group, value) makes every
#   group's values both contiguous and ordered, so the median is an index
#   lookup, not a search. Same trick as the cleaning-pipeline radix fix.
#
#   Verified identical to stats::mad on 60k traces (max abs difference 0).
#
# ARGUMENTS
#   g  : integer group code, 1..ng (NOT character - that is the slow path)
#   v  : values, NAs allowed (dropped per group, like na.rm = TRUE)
#   ng : number of groups
#
# RETURNS  numeric(ng), NA for groups with no non-NA value.
# =============================================================================
grp_median <- function(g, v, ng) {
  ok <- !is.na(v)
  if (!all(ok)) { g <- g[ok]; v <- v[ok] }
  o  <- order(g, v, method = "radix")
  gs <- g[o]; vs <- v[o]

  cnt   <- tabulate(gs, nbins = ng)
  end   <- cumsum(cnt)
  start <- end - cnt + 1L

  res <- rep(NA_real_, ng)
  hit <- cnt > 0L
  # even n -> mean of the two central values; odd n -> lo and hi coincide
  lo  <- start[hit] + (cnt[hit] - 1L) %/% 2L
  hi  <- start[hit] + cnt[hit] %/% 2L
  res[hit] <- (vs[lo] + vs[hi]) / 2
  res
}

grp_mad <- function(g, v, ng) {
  m <- grp_median(g, v, ng)
  1.4826 * grp_median(g, abs(v - m[g]), ng)
}

# =============================================================================
# mad_finite_sample_factor() - MAKE MAD AN UNBIASED sigma AT SHORT TRACE LENGTHS
# =============================================================================
#
# GOAL
#   Remove MAD's small-sample downward bias so a spread threshold in nm means
#   the same thing for a 6-localization trace as for a 40-localization one.
#
# WHY THIS IS NOT OPTIONAL HERE
#   The 1.4826 constant makes MAD consistent for sigma only ASYMPTOTICALLY. At
#   the trace lengths this app actually sees, E[MAD]/sigma is:
#     n = 4  0.735      n = 8  0.887      n = 14 0.940
#     n = 6  0.841      n = 12 0.931      n = 20 0.959
#   Uncorrected, the same nm threshold is therefore far more permissive for
#   SHORT traces - a length-dependent bias applied invisibly inside the filter,
#   which is exactly the failure mode deblink_minpts() documents for minPts.
#
# THE FACTORS
#   1 / E[MAD] under a Gaussian, simulated at 2e5 replicates per n. These
#   reproduce the published Croux & Rousseeuw (1992) correction table to within
#   0.5% (their n = 8: 1.129, here 1.1269), which is the standard reference for
#   this correction. Above n = 20 the usual n / (n - 0.8) form takes over; it is
#   continuous with the table at the join (1.0417 vs 1.0427).
#
# ARGUMENTS  n : integer vector of per-trace localization counts
# RETURNS    numeric multiplier, 1 where n is not usable
# =============================================================================
.MAD_SMALL_N_FACTOR <- c(
  1.1948, 1.4862, 1.3596, 1.2181, 1.1889, 1.1373, 1.1269, 1.0992, 1.0956,
  1.0798, 1.0769, 1.0662, 1.0633, 1.0558, 1.0552, 1.0486, 1.0482, 1.0426,
  1.0427)                                            # n = 2 .. 20

mad_finite_sample_factor <- function(n) {
  n   <- as.integer(n)
  out <- rep(1, length(n))
  small <- !is.na(n) & n >= 2L & n <= 20L
  out[small] <- .MAD_SMALL_N_FACTOR[n[small] - 1L]
  big <- !is.na(n) & n > 20L
  out[big] <- n[big] / (n[big] - 0.8)
  out
}

# =============================================================================
# sd_by_tid() - PER-TRACE, PER-AXIS SPREAD
# =============================================================================
#
# GOAL
#   One spread number per trace per axis, for the Stage-2 spread filter, the
#   spread histograms and the de-blink axis weights.
#
# WHY THE DEFAULT IS THE PLAIN SD
#   This is the number the app REPORTS, and it is the one quantity here that has
#   an outside reference: pyMINFLUX computes exactly this per trace ID
#   (processor/_processor.py, calculate_statistics_on):
#       sx = df_grouped["x"].std()   # pandas, ddof = 1
#   Keeping it identical is what makes MINFLUXR's precision numbers comparable
#   with pyMINFLUX's. pyMINFLUX has NO spread-based trace filter and NO
#   per-localization outlier removal, so nothing downstream of this number is
#   constrained by that comparison - only the number itself.
#
# WHY sd's ZERO BREAKDOWN POINT IS NOT A PROBLEM HERE ANY MORE
#   One stray localization sets an sd. That used to delete whole traces over two
#   bad points. The fix is trace_point_outliers(), which removes the stray
#   localizations BEFORE this runs, rather than swapping the estimator: the
#   contamination is dealt with where it lives instead of being averaged over.
#   With trimming on there is one spread statistic in the whole app.
#
# ARGUMENTS
#   df     : localization table (tid, x, y, optionally z)
#   min_n  : below this, spread is left NA (undefined) rather than guessed
#   robust : FALSE = pyMINFLUX-matching stats::sd (default and user-facing).
#            TRUE  = finite-sample-corrected MAD. INTERNAL machinery only -
#            trace_point_outliers() uses it to set the trimming scale, where a
#            breakdown-resistant sigma is required and comparability is
#            irrelevant. Not exposed as a user choice: two competing readings of
#            "trace spread" in the UI is exactly the confusion to avoid.
#
# RETURNS  data.frame(tid, n_locs, sd_x, sd_y, sd_z), tid in its input type.
# =============================================================================
# The plain-SD form is reimplemented in R from pyMINFLUX (Apache-2.0); the robust
# variant and the spread-based trace filter are additions. See THIRD-PARTY.md section 2.
sd_by_tid <- function(df, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT, robust = FALSE) {
  need_cols <- intersect(c("tid","x","y","z"), names(df))
  if (!all(c("tid","x","y") %in% need_cols)) return(data.frame())

  min_n <- max(2L, as.integer(min_n %||% MIN_LOCS_FOR_SPREAD_DEFAULT))
  has_z <- "z" %in% need_cols

  if (!isTRUE(robust)) {
    d <- df[, need_cols, drop = FALSE]
    d$x <- as.numeric(d$x); d$y <- as.numeric(d$y)
    if (has_z) d$z <- as.numeric(d$z)
    return(d %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        n_locs = dplyr::n(),
        sd_x = if (dplyr::n() >= min_n) stats::sd(x, na.rm = TRUE) else NA_real_,
        sd_y = if (dplyr::n() >= min_n) stats::sd(y, na.rm = TRUE) else NA_real_,
        sd_z = if (has_z && dplyr::n() >= min_n) stats::sd(z, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ))
  }

  # Integer trace codes: grp_median sorts on them, and a character sort here is
  # the ~25x penalty the cleaning pipeline already had to remove once.
  f  <- factor(as.character(df$tid))
  gi <- as.integer(f)
  ng <- nlevels(f)
  if (!ng) return(data.frame())

  n_locs <- tabulate(gi, nbins = ng)
  # The correction is per trace, because it depends on that trace's own n.
  fac <- mad_finite_sample_factor(n_locs)

  ax <- function(nm) {
    if (!nm %in% need_cols) return(rep(NA_real_, ng))
    grp_mad(gi, as.numeric(df[[nm]]), ng) * fac
  }

  # tid must come back in the type it went in as. factor() stringifies it, and
  # returning character where the classic branch returns integer breaks every
  # downstream join on tid ("can't join <integer> with <character>"). The
  # level -> value mapping survives the conversion because as.integer(f) indexes
  # into levels(f), so the pairing is positional, not alphabetical.
  tid_out <- levels(f)
  proto   <- df$tid
  if (is.factor(proto))       tid_out <- factor(tid_out, levels = levels(proto))
  else if (is.integer(proto)) tid_out <- as.integer(tid_out)
  else if (is.numeric(proto)) tid_out <- as.numeric(tid_out)

  out <- data.frame(tid = tid_out, n_locs = n_locs,
                    sd_x = ax("x"), sd_y = ax("y"), sd_z = ax("z"),
                    stringsAsFactors = FALSE)

  # Short traces get NA, not a guess - an unreliable spread must not drive the
  # filter. Matches the non-robust branch above.
  short <- out$n_locs < min_n
  out$sd_x[short] <- NA_real_
  out$sd_y[short] <- NA_real_
  out$sd_z[short] <- NA_real_
  out
}

# =============================================================================
# trace_point_outliers() - FLAG SINGLE LOCALIZATIONS THAT DO NOT FIT THEIR TRACE
# =============================================================================
#
# GOAL
#   Find the individual localizations that sit far from the rest of their OWN
#   trace, so a 12-point trace is not deleted over 2 stray points. Flags rows;
#   it never deletes and never touches whole traces.
#
# THE QUESTION IT ASKS, AND WHY THAT CHOICE MATTERS
#   "Is this point unlike the REST OF THIS TRACE" - not "unlike the dataset".
#   A uniformly broad trace (drift, two molecules, a fiducial) then loses
#   NOTHING here, because all of its points agree with each other; deciding its
#   fate stays with the spread filter and the rogue rules. Only a trace that is
#   tight apart from a few strays loses those strays. Measured on real data
#   (imaging file A, k = 4): a dataset-wide sigma flags 7.6% of all
#   localizations, this rule flags 3.4%, and the difference is almost entirely
#   points inside legitimately broad traces.
#
# THE SCALE: PER-TRACE SIGMA SHRUNK TOWARD THE DATASET
#   A trace's own robust sigma is the right scale but is noisy at n ~ 10: used
#   raw it flags 2.4% of points in CLEAN simulated data, because a trace that
#   happens to draw a small MAD then rejects its own good points. The dataset
#   sigma is stable but wrong for broad traces. So they are blended,
#       sigma_eff = w * sigma_trace + (1 - w) * sigma_pooled,   w = n / (n + 6)
#   which is the usual shrinkage form: short traces borrow the dataset's
#   stability, long traces trust themselves. Measured on clean simulation this
#   costs 0.88% false flags (vs 2.43% per-trace, 0.22% pooled) while still
#   catching 100% of injected contaminants at k = 4.
#
#   sigma is robust (corrected MAD) at BOTH levels on purpose. An sd-based scale
#   would be inflated by the very outliers it is meant to find, so the test
#   would hide them - the estimator is internal machinery, not the reported
#   number (see sd_by_tid).
#
# THE DISTANCE
#   Per-axis whitened, r^2 = sum over axes of ((v - median) / sigma_eff)^2, so
#   the anisotropy of MINFLUX precision (z is not x) does not leak into the
#   threshold. Same metric as apply_deblink_weights.
#
# WHY THE GAUSSIAN k DOES NOT PREDICT THE REJECTION RATE
#   Under a Gaussian null, k = 4 rejects essentially nothing. On real data it
#   flags ~3%. The residuals are genuinely heavier-tailed than Gaussian, so the
#   only honest rejection estimate is the one MEASURED on the loaded data -
#   which is what trim_summary() reports rather than a theoretical number.
#
# PROTECTING SHORT TRACES
#   Below min_trim_n a trace is never trimmed: at n = 5 a single point is 20% of
#   the evidence and the trace median that defines "far" is itself barely
#   determined. Trimming also never takes a trace below keep_n retained points.
#
# ARGUMENTS
#   df         : localization table (tid, x, y, optionally z)
#   k          : reject beyond this many sigma (user-set, default 4)
#   min_trim_n : traces shorter than this are left untouched
#   keep_n     : never trim a trace below this many retained localizations
#   w0         : shrinkage half-weight; n = w0 splits trace and pooled evenly
#
# RETURNS  logical vector, one per ROW of df. TRUE = outlier.
# =============================================================================
trace_point_outliers <- function(df, k = 4, min_trim_n = 6L,
                                 keep_n = MIN_LOCS_FOR_SPREAD_DEFAULT, w0 = 6) {
  n <- nrow(df)
  if (!is.data.frame(df) || !n || !all(c("tid", "x", "y") %in% names(df)))
    return(rep(FALSE, max(0L, n)))
  k <- suppressWarnings(as.numeric(k))
  if (!length(k) || !is.finite(k) || k <= 0) return(rep(FALSE, n))

  axes <- intersect(c("x", "y", "z"), names(df))
  f  <- factor(as.character(df$tid))
  gi <- as.integer(f); ng <- nlevels(f)
  n_tr <- tabulate(gi, nbins = ng)

  # Robust scale at both levels: per trace, and pooled over traces long enough
  # for their own estimate to mean anything.
  sg  <- sd_by_tid(df, min_n = 2L, robust = TRUE)
  sg  <- sg[match(levels(f), as.character(sg$tid)), , drop = FALSE]
  fac <- n_tr / (n_tr + w0)

  r2 <- numeric(n)
  for (ax in axes) {
    v   <- as.numeric(df[[ax]])
    col <- paste0("sd_", ax)
    per <- if (col %in% names(sg)) sg[[col]] else rep(NA_real_, ng)
    pooled <- stats::median(per[n_tr >= max(4L, min_trim_n)], na.rm = TRUE)
    if (!is.finite(pooled) || pooled <= 0) pooled <- stats::median(per, na.rm = TRUE)
    if (!is.finite(pooled) || pooled <= 0) next          # axis carries no spread
    use <- ifelse(is.finite(per) & per > 0, per, pooled)
    sig <- fac * use + (1 - fac) * pooled
    sig[!is.finite(sig) | sig <= 0] <- pooled
    ctr <- grp_median(gi, v, ng)
    r2  <- r2 + ((v - ctr[gi]) / sig[gi])^2
  }
  flag <- sqrt(r2) > k
  flag[is.na(flag)] <- FALSE

  # Short traces are exempt entirely.
  flag[n_tr[gi] < as.integer(min_trim_n)] <- FALSE

  # Never trim a trace below keep_n retained points. A trace that would breach
  # that is not "a good trace with a few strays" - it is a trace this rule
  # cannot describe, so it is left whole for the trace-level filters to judge.
  n_flag <- tabulate(gi[flag], nbins = ng)
  breach <- which(n_tr - n_flag < as.integer(keep_n))
  if (length(breach)) flag[gi %in% breach] <- FALSE

  flag
}

# =============================================================================
# trim_summary() - WHAT THE TRIM ACTUALLY DID, FOR THE INFO BOX
# =============================================================================
# GOAL
#   Report the trim as a measured outcome on THIS dataset, not as a theoretical
#   rejection rate - see trace_point_outliers() on why the Gaussian number is
#   wrong here by a factor of ten.
#
#   The per-trace trimmed FRACTION is the composite-trace diagnostic: a few
#   percent means stray localizations, a large fraction means the trace was
#   never one molecule (simulated, a second molecule contributing 1/12 of the
#   points trims 9%, 2/12 trims 18%, 4/12 trims 36%, 6/12 trims 86%). Reported
#   as a count of traces past a fraction, not applied as a rule, because the
#   fraction is n-dependent and the threshold would be arbitrary.
#
# ARGUMENTS  df : the table the flags were computed on;  flag : that logical
# RETURNS    list of scalars for renderUI
# =============================================================================
trim_summary <- function(df, flag) {
  n <- length(flag)
  if (!n) return(list(n_pts = 0L, n_flag = 0L, pct = 0, n_tr_touched = 0L,
                      pct_tr = 0, n_tr_heavy = 0L, n_tr = 0L))
  tid  <- as.character(df$tid)
  per  <- tapply(flag, tid, mean)
  list(
    n_pts        = n,
    n_flag       = sum(flag),
    pct          = 100 * mean(flag),
    n_tr         = length(per),
    n_tr_touched = sum(per > 0),
    pct_tr       = 100 * mean(per > 0),
    n_tr_heavy   = sum(per >= 0.30)
  )
}

# =============================================================================
# cov_ellipse_params() - THE ERROR ELLIPSE OF A 2x2 COVARIANCE MATRIX
# =============================================================================
#
# GOAL
#   Turn a covariance matrix into a drawable, ROTATED ellipse plus the two
#   numbers that describe the cloud's shape independently of the axes.
#
# WHY ROTATION MATTERS (the thing σx / σy alone cannot express)
#   An axis-aligned ellipse has two parameters and cannot represent a TILTED
#   cloud. Elongation along a diagonal inflates σx and σy by the SAME amount, so
#   σx ≈ σy and the ellipse renders as a circle over a visibly stretched cloud -
#   all the elongation lives in cov(x,y), which an axis-aligned drawing discards.
#
#   The same blind spot afflicts max(σx,σy)/min(σx,σy) as an anisotropy metric:
#   a 45° cloud scores ~1.0, i.e. "isotropic". Since drift within a trace is
#   DIRECTIONAL, a tilted cloud is the expected signature of drift - precisely
#   the artifact the precision panel exists to catch.
#
#   The eigen-decomposition is the fix: eigenvectors give the true major/minor
#   directions, sqrt(eigenvalues) the semi-axes. Their ratio is rotation-
#   invariant, so it detects elongation at ANY angle.
#
# ARGUMENTS  S : symmetric 2x2 covariance matrix;  k : how many σ to scale to
# RETURNS
#   list(ra, rb, theta, ratio, rho, angle_deg) - ra/rb are semi-major/minor in
#   the data's units, theta the major-axis angle in radians. NULL if unusable.
# =============================================================================
cov_ellipse_params <- function(S, k = 1) {
  if (!is.matrix(S) || any(!is.finite(S))) return(NULL)
  ev <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev)) return(NULL)
  vals <- pmax(ev$values, 0)                      # clamp tiny negative round-off
  if (!all(is.finite(vals)) || vals[1] <= 0) return(NULL)
  th <- atan2(ev$vectors[2, 1], ev$vectors[1, 1]) # direction of the major axis
  den <- sqrt(max(S[1, 1] * S[2, 2], 1e-12))
  list(
    ra    = k * sqrt(vals[1]),
    rb    = k * sqrt(vals[2]),
    theta = th,
    ratio = sqrt(vals[1] / max(vals[2], 1e-12)),  # rotation-invariant anisotropy
    rho   = S[1, 2] / den,                        # correlation, for reference
    angle_deg = ((th * 180 / pi) + 180) %% 180    # fold to [0,180): an axis has no sign
  )
}

# GOAL: the outline of a rotated ellipse as plain x/y vectors.
# Drawn as a TRACE rather than a layout shape: a shape cannot be switched off
# from the legend, and these rings are exactly what a reader wants to toggle.
ellipse_xy <- function(ra, rb, theta, cx = 0, cy = 0, n = 181L) {
  if (!all(is.finite(c(ra, rb, theta, cx, cy))) || ra <= 0 || rb <= 0) return(NULL)
  t  <- seq(0, 2 * pi, length.out = n)
  ct <- cos(theta); st <- sin(theta)
  list(x = cx + ra * cos(t) * ct - rb * sin(t) * st,
       y = cy + ra * cos(t) * st + rb * sin(t) * ct)
}

# =============================================================================
# density_hdr_2d() - SMOOTH DENSITY + HIGHEST-DENSITY-REGION CONTOURS
# =============================================================================
#
# GOAL
#   Estimate the 2D density of a point cloud on a fixed frame and return both
#   (a) a grid smooth enough to draw as a continuous heatmap, and (b) contour
#   lines at chosen PROBABILITY MASSES rather than at arbitrary density levels.
#
# WHY NOT A BINNED HISTOGRAM-CONTOUR (what the precision cloud used to draw)
#   A 2D histogram estimates each bin independently, so every bin carries full
#   Poisson noise. For the precision cloud (~40k points, 60x60 bins) that is
#   ~9% noise at the peak but ~26% at the 2-sigma radius and ~90% at 3 sigma -
#   i.e. exactly the outer contours the eye reads as "the shape of the cloud"
#   are almost pure shot noise, which is what makes the plot look ragged and
#   under-smoothed. Re-binning cannot fix it: coarser bins blur the core,
#   finer bins make the tails worse. A kernel estimate is the fix, because each
#   point contributes to a neighbourhood and the estimate is continuous.
#
# WHY HDR LEVELS
#   A contour at "density = 0.004" means nothing on its own. A contour holding
#   half the points means something a reader can act on. Drawn together with a
#   fitted ellipse enclosing the SAME mass, the empirical contour and the
#   theoretical ring coincide when the cloud really is Gaussian and separate
#   visibly where it is not. That comparison is the whole diagnostic, and it is
#   unavailable with count-based levels.
#
#   The masses are round numbers on purpose. Any pair works, so they are chosen
#   to be readable (50%, 90%) and the sigma multiple is derived from them - for
#   a 2D Gaussian the k-sigma ellipse holds 1 - exp(-k^2/2), so 50% sits at
#   1.177 sigma and 90% at 2.146 sigma. Choosing the masses the other way round
#   (0.393, 0.865, i.e. exactly 1 and 2 sigma) makes the model ring easier to
#   label and the plot harder to read; the reader is the one to optimise for.
#
# BANDWIDTH
#   Scott's rule for 2D (sigma * n^(-1/6)) with a floor at the grid pitch -
#   below the pitch a kernel cannot smooth anything, it just re-creates the
#   binning artefact. sigma is passed in (use a ROBUST sigma) rather than taken
#   from the data, so a fiducial tail cannot inflate the bandwidth.
#
# ARGUMENTS
#   a, b    : coordinates (same units)
#   lim     : frame half-width; the grid spans [-lim, lim] on both axes
#   sigma   : scale for the bandwidth rule (robust sigma of the cloud)
#   n_grid  : grid points per axis
#   probs   : probability masses to contour
#   max_pts : subsample cap for the kernel sum (cost is O(n * n_grid))
# RETURNS
#   list(x, y, z, bw_sd, mass_in_frame, paths) where paths is a list of
#   list(prob=, d=, x=, y=) contour outlines, or NULL if the cloud is unusable.
#   A requested mass that extends beyond the frame is skipped rather than drawn
#   wrong. Draw the outlines with plotly::add_paths(), NOT add_lines() -
#   add_lines() sorts by x and turns a closed contour into a filled zigzag.
# =============================================================================
density_hdr_2d <- function(a, b, lim, sigma, n_grid = 160L,
                           probs = c(0.50, 0.90), max_pts = 20000L) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 20L || !is.finite(lim) || lim <= 0) return(NULL)
  if (length(a) > max_pts) {                       # kernel sum is O(n * n_grid)
    s <- sort(sample.int(length(a), max_pts)); a <- a[s]; b <- b[s]
  }
  pitch <- 2 * lim / (n_grid - 1)                  # grid spacing
  if (!is.finite(sigma) || sigma <= 0) sigma <- stats::mad(sqrt(a^2 + b^2))
  if (!is.finite(sigma) || sigma <= 0) return(NULL)
  # Scott's rule, floored at the grid pitch so the kernel always spans >1 cell.
  bw_sd <- max(sigma * length(a)^(-1 / 6), pitch)
  # MASS::kde2d's `h` is 4x the kernel SD (it uses dnorm(sd = h/4) internally).
  k <- tryCatch(MASS::kde2d(a, b, h = c(4 * bw_sd, 4 * bw_sd), n = n_grid,
                            lims = c(-lim, lim, -lim, lim)),
                error = function(e) NULL)
  if (is.null(k) || !all(is.finite(k$z))) return(NULL)

  # HDR level for mass p: walk the grid cells from densest down, accumulating
  # probability (density x cell area) until p is reached. The density there is
  # the level whose contour encloses exactly that mass.
  cell <- (k$x[2] - k$x[1]) * (k$y[2] - k$y[1])
  srt  <- sort(as.vector(k$z), decreasing = TRUE)
  cm   <- cumsum(srt) * cell                       # mass enclosed at each level
  in_frame <- cm[length(cm)]                       # total mass inside the frame
  paths <- list()
  for (p in probs) {
    # If the frame does not even contain p of the cloud, no honest contour for p
    # exists here - skip it rather than draw the frame edge and imply otherwise.
    if (!is.finite(in_frame) || in_frame < p) next
    idx <- which(cm >= p)[1]
    lev <- srt[idx]
    if (!is.finite(lev) || lev <= 0) next
    cls <- tryCatch(grDevices::contourLines(k$x, k$y, k$z, levels = lev),
                    error = function(e) NULL)
    for (cl in cls) {
      if (length(cl$x) < 3L) next
      # Plain coordinates: the caller draws these as a TRACE (so the legend can
      # toggle them), which needs x/y, not an SVG path string.
      paths[[length(paths) + 1L]] <- list(prob = p, d = lev, x = cl$x, y = cl$y)
    }
  }
  list(x = k$x, y = k$y, z = k$z, bw_sd = bw_sd, mass_in_frame = in_frame, paths = paths)
}

# =============================================================================
# robust_median_cov() - ONE VOTE PER TRACE, BUT STILL TILT-AWARE
# =============================================================================
# GOAL
#   Combine many per-trace covariance matrices into one robust matrix by taking
#   the element-wise median, so a single 2000-localization fiducial cannot
#   dominate the shape the way it dominates a pooled estimate.
#
#   Element-wise medians are not guaranteed to stay positive semi-definite, so
#   the off-diagonal is shrunk if it would exceed what the diagonal supports -
#   without that, eigen() can return a negative eigenvalue and the ellipse
#   collapses or inverts.
#
# ARGUMENTS  vxx, vyy : per-trace variances;  vxy : per-trace covariances
# RETURNS    2x2 matrix, or NULL when there is nothing usable.
# =============================================================================
robust_median_cov <- function(vxx, vyy, vxy) {
  m <- function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) }
  a <- m(vxx); b <- m(vyy); c0 <- m(vxy)
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(NULL)
  if (!is.finite(c0)) c0 <- 0
  lim <- 0.99 * sqrt(a * b)                        # keep |rho| < 1
  if (abs(c0) > lim) c0 <- sign(c0) * lim
  matrix(c(a, c0, c0, b), 2L, 2L)
}

# =============================================================================

sigma_label_med_pm_sd <- function(v, digits = 2) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)

  med <- stats::median(v, na.rm = TRUE)
  mad <- stats::mad(v, constant = 1.4826, na.rm = TRUE)  # normalized MAD (= MAD / 0.67449), matches pyMINFLUX

  if (!is.finite(med)) return(NULL)
  if (!is.finite(mad)) mad <- NA_real_

  # σ̃ is "sigma with tilde" meaning median-ish notation
  if (is.finite(mad)) {
    sprintf("\u03C3\u0303 %.2f \u00B1 %.2f nm", round(med, digits), round(mad, digits))
  } else {
    sprintf("\u03C3\u0303 %.2f nm", round(med, digits))
  }
}

# =============================================================================
# add_sigma_iqr_lines() - THE QUARTILE PAIR THAT FLANKS THE MEDIAN LINE
# =============================================================================
#
# GOAL
#   Draw thin lines at the 25th and 75th percentile of a spread distribution,
#   flanking the median line that add_sigma_vline() already draws.
#
# WHY
#   Matches pyMINFLUX's add_median_line(..., add_iqr = True), which puts dashed
#   quartile lines either side of its median marker on the sx/sy/sz histograms.
#   It also earns its place independently: median +/- nMAD is symmetric, while
#   the distribution of per-trace spread is right-skewed (a spread cannot go
#   below zero but has no ceiling). The quartiles show that asymmetry, which the
#   +/- label cannot.
#
# ARGUMENTS
#   fig    : plotly figure to add to
#   v      : the per-trace spread values the median line was drawn from
#   y_top  : same panel ceiling used for the median line
#   color  : match the median line this pair belongs to
#
# RETURNS  the figure, unchanged when the quartiles are not finite.
# =============================================================================
add_sigma_iqr_lines <- function(fig, v, y_top, color) {
  v <- as.numeric(v); v <- v[is.finite(v)]
  if (length(v) < 4L || !is.finite(y_top) || y_top <= 0) return(fig)
  q <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(q)) || q[1] == q[2]) return(fig)

  for (xq in q) {
    fig <- fig %>% plotly::add_segments(
      x = xq, xend = xq,
      y = 0, yend = y_top,
      line = list(color = color, width = 1, dash = "dot"),
      inherit = FALSE, hoverinfo = "skip", showlegend = FALSE
    )
  }
  fig
}

add_sigma_vline <- function(fig, x, y_top, color, label, dash = NULL) {
  if (!is.finite(x) || !is.finite(y_top) || y_top <= 0) return(fig)
  
  # Add the line
  fig <- fig %>% plotly::add_segments(
    x = x, xend = x,
    y = 0, yend = y_top,
    line = c(list(color = color, width = 2), if (!is.null(dash)) list(dash = dash) else list()),
    inherit = FALSE,
    hoverinfo = "skip",
    showlegend = FALSE
  )
  
  # Add the sigma label slightly above the line
  fig %>% plotly::layout(
    annotations = c(
      fig$x$layout$annotations %||% list(),
      list(
        x = x, y = y_top,
        xref = "x", yref = "y",
        text = label,
        showarrow = FALSE,
        xanchor = "left",
        yanchor = "bottom",
        font = list(color = color, size = 12),
        bgcolor = "rgba(255,255,255,0.0)"
      )
    )
  )
}

# 2D convex-hull centroid (polygon centroid of the hull)
# Input: coords_matrix (N x 2) with columns [x, y]
# Output: numeric length-2 c(cx, cy)
# --- per-trace SD table from currently filtered data ---
# sd_by_tid is defined above (canonical helper).

# --- median marker lines for any plotly histogram -------------------------
# GOAL: put the summary statistic ON the distribution it summarises, so a shifted
# or bimodal histogram can be read against a number instead of by eye.
#
# Draws one dashed vertical line per series, spanning the full panel height
# (yref = "paper", so it does not need the histogram's y scale) with the median
# value as a label at the top. Lines are drawn as layout SHAPES rather than
# traces, so they add nothing to the legend and cannot be toggled off with a
# series - which also means they do NOT follow legend clicks; hide the panel's
# median lines with the panel's own checkbox instead.
#
# Input : fig     - a plotly object
#         meds    - named numeric vector, name = series label, value = median
#         colors  - named character vector of line colors (same names as `meds`)
#         digits  - significant digits in the label
#         unit    - appended to the label ("nm" etc.)
# Output: the plotly object with shapes + annotations appended
add_median_vlines <- function(fig, meds, colors = NULL, digits = 3, unit = "") {
  meds <- meds[is.finite(meds)]
  if (!length(meds)) return(fig)
  nms  <- names(meds) %||% rep("", length(meds))
  # Stagger the labels so two nearby medians do not print on top of each other.
  ytop <- 1 - 0.07 * (seq_along(meds) - 1L)

  shapes <- lapply(seq_along(meds), function(i) {
    col <- if (!is.null(colors) && nms[i] %in% names(colors)) unname(colors[nms[i]]) else "#333333"
    list(type = "line", xref = "x", yref = "paper",
         x0 = meds[[i]], x1 = meds[[i]], y0 = 0, y1 = 1,
         line = list(color = col, width = 2, dash = "dash"))
  })
  anns <- lapply(seq_along(meds), function(i) {
    col <- if (!is.null(colors) && nms[i] %in% names(colors)) unname(colors[nms[i]]) else "#333333"
    list(x = meds[[i]], y = ytop[i], xref = "x", yref = "paper",
         yanchor = "bottom", xanchor = "left",
         text = sprintf("med %s%s", signif(meds[[i]], digits),
                        if (nzchar(unit)) paste0(" ", unit) else ""),
         showarrow = FALSE, font = list(color = col, size = 11))
  })

  fig %>% plotly::layout(
    shapes      = c(fig$x$layout$shapes %||% list(), shapes),
    annotations = c(fig$x$layout$annotations %||% list(), anns))
}

# ---- Batch NN helpers -------------------------------------------------------

# Apply the metadata, outlier, length and spread filters to a raw stored df (no
# bypass logic - the ranges are applied as-is so that the same thresholds work
# across all loaded files in batch mode).
#   settings : smlm_settings or Shiny input (the "filter" settings)
apply_meta_filters_batch <- function(df, settings) {
  s <- as_smlm_settings(settings)
  hard_filter <- function(dat, nm, rng_val) {
    if (!nm %in% names(dat)) return(dat)
    if (is.null(rng_val) || length(rng_val) != 2 || any(!is.finite(as.numeric(rng_val)))) return(dat)
    lo <- rng_val[1]; hi <- rng_val[2]
    dplyr::filter(dat, is.na(.data[[nm]]) | (.data[[nm]] >= lo & .data[[nm]] <= hi))
  }
  df <- df %>%
    hard_filter("efo",   s$efo_range)   %>%
    hard_filter("cfr",   s$cfr_range)   %>%
    hard_filter("dcr",   s$dcr_range)   %>%
    hard_filter("eco",   s$eco_range)   %>%
    hard_filter("tim",   s$tim_range)   %>%
    hard_filter("dwell", s$dwell_range) %>%
    hard_filter("fluo",  s$fluo_range)
  # Outlier localizations - Stage 1b, and it MUST come before the length filter
  # below, exactly as in the interactive pipeline: a trimmed point is not a
  # localization of that trace, so it must not count toward its length either.
  if (isTRUE(s$trim_outliers) && "tid" %in% names(df) && nrow(df)) {
    fl <- tryCatch(trace_point_outliers(df, k = s$trim_k, keep_n = s$sd_min_locs),
          error = function(e) NULL)
    if (!is.null(fl) && length(fl) == nrow(df)) df <- df[!fl, , drop = FALSE]
  }
  # TID length filter
  tlr <- s$tid_length_range
  if ("tid" %in% names(df) && !is.null(tlr) && length(tlr) == 2 && all(is.finite(tlr))) {
    tl   <- table(df$tid)
    keep <- names(tl)[tl >= tlr[1] & tl <= tlr[2]]
    df   <- dplyr::filter(df, tid %in% keep)
  }
  # Trace spread (σ) filter, the batch form of the Display tab's spread filter
  if (isTRUE(s$view_SD) && "tid" %in% names(df)) {
    tx <- s$target_sd_x; ty <- s$target_sd_y
    if (!is.null(tx) && !is.null(ty) && is.finite(tx) && is.finite(ty) && tx > 0 && ty > 0) {
      min_n <- max(2L, s$sd_min_locs)
      # Same estimator as the interactive Stage-2 filter, on the same trimmed
      # point set, or batch and interactive report different trace counts from
      # one setting. Trimming already happened above.
      spread <- sd_by_tid(df, min_n = min_n)
      # Traces with undefined spread (NA, too short) are kept, not filtered.
      valid <- spread$tid[dplyr::coalesce(spread$sd_x, 0) <= tx &
                          dplyr::coalesce(spread$sd_y, 0) <= ty]
      df <- dplyr::filter(df, tid %in% valid)
    }
  }
  df
}

remove_rogue_traces <- function(df, rule = c("factor_median","mad_rule","axis_sd"),
                                factor = 2.5, min_n = MIN_LOCS_FOR_SPREAD_DEFAULT,
                                robust = TRUE) {
  rule <- match.arg(rule)

  # The per-trace spread MUST use the same estimator as the Stage-2 filter.
  # These rules were robust ACROSS traces (median/MAD of the per-trace values)
  # but the per-trace value itself was a plain sd, so a trace with 2 stray
  # localizations that the spread filter deliberately kept would be deleted here
  # instead - the two stages would silently disagree about what a rogue is.
  spread <- sd_by_tid(df, min_n = 2L, robust = robust)
  if (!nrow(spread)) return(list(df = df, removed_tids = character(0),
                                 removed_n = 0L, note = "no traces"))
  names(spread)[names(spread) == "n_locs"] <- "n"

  keep <- rep(TRUE, nrow(spread))
  note <- ""
  
  if (rule == "factor_median") {
    med_all <- c(stats::median(spread$sd_x, na.rm = TRUE),
                 stats::median(spread$sd_y, na.rm = TRUE),
                 if ("z" %in% names(df)) stats::median(spread$sd_z, na.rm = TRUE) else NA)
    thr_x <- factor * med_all[1]
    thr_y <- factor * med_all[2]
    thr_z <- if (!is.na(med_all[3])) factor * med_all[3] else Inf
    keep <- (spread$sd_x <= thr_x) & (spread$sd_y <= thr_y) & (spread$sd_z <= thr_z | is.na(spread$sd_z))
    note <- sprintf("factor_median (×%.2f)", factor)
    
  } else if (rule == "mad_rule") {
    mx <- stats::median(spread$sd_x, na.rm = TRUE); madx <- stats::mad(spread$sd_x, na.rm = TRUE)
    my <- stats::median(spread$sd_y, na.rm = TRUE); mady <- stats::mad(spread$sd_y, na.rm = TRUE)
    if ("z" %in% names(df)) { mz <- stats::median(spread$sd_z, na.rm = TRUE); madz <- stats::mad(spread$sd_z, na.rm = TRUE) } else { mz <- NA; madz <- NA }
    thr_x <- mx + factor * madx; thr_y <- my + factor * mady; thr_z <- if (!is.na(mz)) mz + factor * madz else Inf
    keep <- (spread$sd_x <= thr_x) & (spread$sd_y <= thr_y) & (spread$sd_z <= thr_z | is.na(spread$sd_z))
    note <- sprintf("mad_rule (median + %.2f × MAD)", factor)
    
  } else if (rule == "axis_sd") {
    medx <- stats::median(spread$sd_x, na.rm = TRUE)
    medy <- stats::median(spread$sd_y, na.rm = TRUE)
    medz <- if ("z" %in% names(df)) stats::median(spread$sd_z, na.rm = TRUE) else NA
    thr_x <- factor * medx; thr_y <- factor * medy; thr_z <- if (!is.na(medz)) factor * medz else Inf
    keep <- (spread$sd_x <= thr_x) & (spread$sd_y <= thr_y) & (spread$sd_z <= thr_z | is.na(spread$sd_z))
    note <- sprintf("axis_sd (any axis > %.2f × median(axis SD))", factor)
  }
  
  keep <- keep & (spread$n >= min_n)
  kept_tids <- spread$tid[keep]
  removed_tids <- spread$tid[!keep]
  
  list(
    df = df %>% dplyr::filter(tid %in% kept_tids),
    removed_tids = removed_tids,
    removed_n = length(removed_tids),
    note = note
  )
}

# =============================================================================
# geometric_median() - A ROBUST CENTRE THAT DOES NOT DEPEND ON THE AXES
# =============================================================================
#
# GOAL
#   The point minimising the sum of DISTANCES to a set of points, found by
#   Weiszfeld iteration. Used to summarise a cloud of displacement vectors where
#   the answer is a direction, so the estimator must not care how x and y happen
#   to be oriented.
#
# WHY NOT A PER-AXIS MEDIAN
#   Taking median(x) and median(y) separately is NOT rotation-equivariant:
#   rotate the frame, take the medians, rotate back, and you get a different
#   vector. Measured on displacement clouds like the drift plot's, that
#   disagreement reaches 6.8 nm - which for a plot whose entire content is a
#   DIRECTION means the arrow would depend on the microscope's axis convention.
#   The mean and the geometric median are both exactly equivariant (0.000 nm).
#
# WHY NOT THE MEAN
#   The mean has no breakdown point. On these clouds the contamination is real:
#   a "molecule" can hold two genuinely different fluorophores, which puts a
#   displacement of tens of nm into an otherwise ~10 nm cloud. Measured RMSE
#   against a known 10 nm shift, jitter sd 12 nm, n = 20 per bin:
#
#     contamination        mean    per-axis median   geometric median
#     none                 4.55        5.69              4.95
#     10% wrong molecule   5.70        6.01              5.54
#     25% wrong molecule   7.59        6.81              6.07
#     10% gross outliers  15.65        6.08              5.60
#
#   The geometric median is best or tied everywhere contamination exists and
#   costs ~9% efficiency when there is none, which is the right trade for a
#   panel whose whole purpose is judging a real sample.
#
# The implementation (Weiszfeld from the coordinate-wise median) follows.
# =============================================================================

# =============================================================================
# geometric_median() - ROBUST MULTIVARIATE CENTRE (WEISZFELD)
# =============================================================================
#
# GOAL
#   Return the point m minimising the SUM OF EUCLIDEAN DISTANCES to every row of
#   `pts` (the L1 / spatial median). This is the centre used by the heuristic
#   shape fits and as the MLE starting centre.
#
# WHY NOT THE CONVEX-HULL CENTROID (what these fits used before)
#   compute_convex_hull_centroid() is an area-weighted centroid of the hull
#   SURFACE, so it is defined entirely by the extreme points: one stray member
#   drags it, i.e. its breakdown point is 0. The geometric median has breakdown
#   0.5 - half the members must be corrupted before it can be moved without
#   bound - which is what we want on MINFLUX clusters that may still carry a
#   loosely-attached trace or two.
#
# WHY NOT A COORDINATE-WISE median(x), median(y), median(z)
#   That is not rotation-equivariant: rotate the structure and the centre moves
#   relative to it. Since these centres are handed to prcomp() to define ring
#   normals and ellipsoid axes, the centre must rotate WITH the data. The
#   geometric median does; the per-axis median does not.
#
# WHAT THIS IS NOT
#   It is a centre OF THE POINTS, not a centre of curvature. When a structure is
#   only partially labelled (a cap of a vesicle, an arc of a ring) no
#   centre-of-the-points statistic recovers the sphere/ring centre - the
#   geometric median then sits inside the observed patch. Use the MLE pathway
#   (mle_fit_*, which optimises the centre itself) when the centre is the
#   measurement; this is a robust starting point and a heuristic answer.
#
# ALGORITHM
#   Weiszfeld (1937) fixed-point iteration
#       m <- SUM_i (p_i / d_i) / SUM_i (1 / d_i),    d_i = ||p_i - m||
#   started from the coordinate-wise median (already robust, so it lands close
#   and converges in a handful of steps). The classic degeneracy is d_i = 0 when
#   an iterate comes to rest exactly on a data point; d_i is floored at `eps_d`
#   so the weight stays finite.
#
#   Deterministic and exact for the small member counts here (tens to hundreds),
#   which is why this is used rather than Gmedian::Gmedian - that is a stochastic
#   averaged-gradient estimator built for large n.
geometric_median <- function(pts, max_iter = 200L, tol = 1e-7, eps_d = 1e-9) {
  if (is.null(pts)) return(NULL)
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]   # a partial row has no position
  n <- nrow(pts)
  if (n == 0L) return(rep(NA_real_, ncol(pts)))
  if (n <= 2L) return(colMeans(pts))       # 1 point is its own median; for 2 the
                                           # whole segment is optimal - take the midpoint
  m <- apply(pts, 2, stats::median)        # robust starting point
  for (i in seq_len(max_iter)) {
    d <- sqrt(rowSums(sweep(pts, 2, m)^2))
    d[d < eps_d] <- eps_d                  # iterate sitting on a data point
    w <- 1 / d
    m_new <- colSums(pts * w) / sum(w)
    if (max(abs(m_new - m)) < tol) { m <- m_new; break }
    m <- m_new
  }
  unname(m)
}

# =============================================================================
# smlm_loc_precision_nm() - PER-LOCALIZATION PRECISION OF THE CURRENT DATASET
# =============================================================================
#
# GOAL
#   Give the isosurface a PHYSICAL length scale to hang its geometry on, so the
#   drawn size of a structure stops depending on how well it happened to be
#   labelled.
#
# HOW
#   Each trace is one fluorophore sitting still, so the scatter of its own
#   localizations around their centre IS the localization precision. Take the
#   per-trace SD on each axis and report the median trace, matching the
#   convention already used by estimate_gt_positions() (sigma_loc_nm) and by the
#   de-blinking radius.
#
#   The MEDIAN over traces, not the mean: a trace that caught two molecules in
#   sequence has a spread of tens of nm and would drag a mean far off the real
#   single-molecule precision.
#
# WHY NOT THE REPORTED CRLB
#   Not every export carries one, and where it does it is a lower bound on a
#   single localization rather than a measured scatter. The measured spread
#   already contains whatever residual drift and stage vibration survived, which
#   is what actually blurs the rendered structure.
#
# THE AXES ARE NOT EQUAL, AND NOT IN THE DIRECTION CAMERA SMLM TRAINS YOU TO
# EXPECT. MINFLUX z is typically TIGHTER than x/y, the opposite of astigmatism
# or biplane imaging. Measured on the a real imaging file files: sx 6.84 / sy 6.81 /
# sz 4.11 nm on B2-010 and 7.07 / 7.16 / 4.17 nm on B3-008, a z:lateral ratio of
# 0.60 and 0.59. Anything that wraps a structure isotropically therefore
# over-wraps it axially by about two thirds, so `per_axis = TRUE` returns the
# three numbers separately and the isosurface uses them.
#
# RETURNS
#   `per_axis = FALSE`: one number in nm (the RMS over the axes in use).
#   `per_axis = TRUE` : c(x, y, z) in nm; z is the lateral mean on 2-D data.
#   `fallback` (recycled) when there are too few multi-localization traces.
smlm_loc_precision_nm <- function(df, min_n = 3L, fallback = 2, per_axis = FALSE) {
  bail <- function() if (isTRUE(per_axis)) rep(fallback, 3L) else fallback
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(bail())
  if (!all(c("tid", "x", "y") %in% names(df))) return(bail())

  # z only counts when the dataset really is 3-D; a 2-D file carries a constant
  # (or all-zero) z whose SD is 0 and would halve the estimate.
  use_z <- "z" %in% names(df) && any(is.finite(df$z) & df$z != 0)

  per_tid <- tryCatch(
    df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        sx = stats::sd(x, na.rm = TRUE),
        sy = stats::sd(y, na.rm = TRUE),
        sz = if (use_z) stats::sd(z, na.rm = TRUE) else NA_real_,
        n  = dplyr::n(), .groups = "drop"
      ) %>%
      dplyr::filter(n >= min_n, is.finite(sx), is.finite(sy)),
    error = function(e) NULL
  )
  if (is.null(per_tid) || !nrow(per_tid)) return(bail())

  # median trace on each axis independently, so a real anisotropy survives
  ax <- c(
    x = stats::median(per_tid$sx, na.rm = TRUE),
    y = stats::median(per_tid$sy, na.rm = TRUE),
    z = if (use_z) stats::median(per_tid$sz, na.rm = TRUE) else NA_real_
  )
  if (!is.finite(ax[["z"]]) || ax[["z"]] <= 0) ax[["z"]] <- mean(ax[c("x", "y")])
  ax[!is.finite(ax) | ax <= 0] <- fallback

  if (isTRUE(per_axis)) return(ax)

  out <- sqrt(mean(ax^2))                       # RMS over the axes
  if (!is.finite(out) || out <= 0) fallback else out
}

