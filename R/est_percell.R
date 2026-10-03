## ============================================================================
## est_percell.R - one row of results per cell (file), for comparing conditions
## ----------------------------------------------------------------------------
## Each file is one field of view, i.e. one cell, and is the unit of replication.
## Every file is analysed on its own: its own cycle time, its own gap split, and
## the TCP diameter, pattern, prefactor and wander floor of the sequence it was
## acquired with. Tracks from one cell are not independent measurements, so
## conditions are compared on these per-cell values, never on pooled tracks.
##
##   percell_summary() ...... the per-cell pipeline for one file -> one-row data.frame
##   percell_condition_stats() mean, SD and t-based 95% CI across cells per condition
##   PERCELL_METRICS ........ the reported columns, with labels and units
## ============================================================================

# Reported per-cell columns: name -> axis label. Shared by the server's metric
# picker and plot so the two cannot disagree.
PERCELL_METRICS <- c(
  D_msd_um2_s  = "D (upper-regime MSD, µm²/s)",
  D_um2_s      = "D (CVE, lag 1, µm²/s)",
  D_corr_um2_s = "D − D_w (µm²/s)",
  Dw_frac      = "D_w / D",
  sigma_nm     = "σ_loc, photon budget (nm)",
  eco_med      = "median ECO (photons)",
  t_c_ms       = "cycle time (ms)",
  frac_single  = "fraction of single-cycle steps",
  n_tracks     = "tracks after filtering",
  len_med      = "median track length (localizations)",
  alpha_short  = "short-lag α (MSD lags 1-8)",
  r1_excess    = "r₁ above the motion bound (>0: correlated error)"
)

# ── percell_summary() ─────────────────────────────────────────────────────────
# GOAL: the per-cell numbers for ONE file, computed exactly as the main tab
# computes them for that file when it is the active dataset: cycle time from the
# data (cycle_cluster), the whole cleaning and filter pipeline (filter_pipeline:
# gap split, end-trim, continuity gate, range, length and duration filters), then
# the CVE, the photon-budget sigma and the smoothing check on what survives.
#
# df      : localizations of one file (tid, x, y, tim in ms, optional z, eco)
# file    : label for the row
# L       : TCP diameter of the final iteration (nm)
# pattern : TCP pattern (hexagon, square, triangle, octahedron)
# a       : photon-budget prefactor for this sequence
# D_w     : wander floor for this sequence (nm^2/ms), 0 = not measured
# nu      : doughnut FWHM (nm), 0 = beam term off
# msd     : the MSD panel's fit settings (window, track_weight, iterative, linear_tau,
#           R, coef = c(d2, d3), min_pts) for the headline upper-regime D; NULL skips it
# filters : the main tab's filter settings, as filter_pipeline() takes them without
#           t_c / sw (list(split, trim, gate, ranges, tid_len, dur)). The spatial crop
#           is not applied: it belongs to one field of view.
# RETURNS one-row data.frame with the PERCELL_METRICS columns plus file, n_locs,
#   r1_obs, r1_bound and a `note` for anything that could not be computed.
percell_summary <- function(df, file, L, pattern = TCP_PATTERN_DEFAULT, a = TCP_CRLB_PREFACTOR,
                            D_w = 0, nu = 0, msd = NULL,
                            filters = list(split = list(enable = TRUE, adaptive = FALSE,
                                                        gap_cycles = 7, min_len = 50L, maxoff_ms = 3.0))) {
  row <- data.frame(file = file, n_locs = nrow(df), n_tracks = NA_integer_, len_med = NA_real_,
                    t_c_ms = NA_real_, frac_single = NA_real_, eco_med = NA_real_,
                    sigma_nm = NA_real_, D_msd_um2_s = NA_real_, D_um2_s = NA_real_, D_corr_um2_s = NA_real_,
                    Dw_frac = NA_real_, alpha_short = NA_real_, r1_obs = NA_real_,
                    r1_bound = NA_real_, r1_excess = NA_real_, note = "",
                    stringsAsFactors = FALSE)
  if (!nrow(df) || !all(c("tid", "x", "y", "tim") %in% names(df))) {
    row$note <- "no localizations"; return(row)
  }
  # 1. this cell's own cycle time: sequences with different dwell times differ
  cc <- cycle_cluster(df$tid, df$tim)
  if (is.null(cc) || !is.finite(cc$cycle) || cc$cycle <= 0) {
    row$note <- "too few intervals for a cycle time"; return(row)
  }
  t_c <- cc$cycle; row$t_c_ms <- t_c

  # 2. the whole filter pipeline with this cell's cycle time (the main tab's settings)
  fp <- filters; fp$t_c <- t_c; fp$sw <- cc$sw
  df <- filter_pipeline(df, fp)
  if (!nrow(df)) { row$note <- "nothing left after the filters"; return(row) }
  tl <- table(df$tid)
  row$n_tracks <- length(tl); row$len_med <- stats::median(as.numeric(tl))

  # 3. CVE on the filtered data, as the Data Quality panel does
  has_z <- has_axial_data(df, min_sd = 1, min_n = 11L)
  st <- cve_steps(df, t_c, has_z)
  row$frac_single <- st$frac_integ
  cv <- cve_core(st$dx, st$dy, st$dz, st$dt, st$integ_step, has_z, t_c)
  if (is.null(cv)) { row$note <- "fewer than 20 single-cycle step pairs"; return(row) }
  D_cve <- cv$D_cve                                      # nm^2/ms
  row$D_um2_s <- D_cve / 1e3                             # 1 nm^2/ms = 1e-3 um^2/s

  # headline D: the upper lag regime of the ensemble MSD (msd_headline, est_msd.R)
  if (!is.null(msd)) {
    hd <- tryCatch(msd_headline(df, t_c, cc$sw, msd$window, msd$track_weight %||% TRUE,
                                msd$iterative %||% TRUE, msd$linear_tau %||% TRUE, msd$R %||% 0,
                                if (has_z) msd$coef[["d3"]] else msd$coef[["d2"]],
                                msd$min_pts %||% 50L, has_z), error = function(e) NULL)
    if (!is.null(hd)) row$D_msd_um2_s <- hd$D / 1e3
  }
  D_head <- if (is.finite(row$D_msd_um2_s)) row$D_msd_um2_s * 1e3 else D_cve
  if (is.finite(D_w) && D_w > 0) {                       # the floor against the headline D
    row$D_corr_um2_s <- (D_head - D_w) / 1e3
    row$Dw_frac      <- if (D_head > 0) D_w / D_head else NA_real_
  }

  # 4. photon-budget sigma with THIS cell's sequence settings
  if ("eco" %in% names(df)) {
    eco <- df$eco[is.finite(df$eco) & df$eco > 0]
    if (length(eco)) {
      row$eco_med  <- stats::median(eco)
      row$sigma_nm <- stats::median(.crlb_doughnut(L, nu, eco, D = max(0, D_cve),
                                                   t_loc = cv$dt_int, pref = a, pattern = pattern))
    }
  }

  # 5. smoothing check: r1 above the fGn bound means damped / smoothed positions
  sm <- smoothing_check(st$df, t_c, cv$dt_int)
  row$alpha_short <- sm$alpha_short
  row$r1_obs      <- cv$r1_norm
  row$r1_bound    <- sm$r1_fgn
  row$r1_excess   <- cv$r1_norm - sm$r1_fgn
  row
}

# ── percell_condition_stats() ─────────────────────────────────────────────────
# GOAL: summarise one metric per condition with the CELL as the unit: number of
# cells, mean, SD and a t-based 95% CI of the mean across cells.
# tab    : per-cell table with a `condition` column (blank = unassigned, skipped)
# metric : column name
# RETURNS data.frame(condition, n_cells, mean, sd, ci_lo, ci_hi)
percell_condition_stats <- function(tab, metric) {
  tab <- tab[nzchar(trimws(tab$condition)) & is.finite(tab[[metric]]), , drop = FALSE]
  if (!nrow(tab)) return(data.frame(condition = character(0), n_cells = integer(0), mean = numeric(0),
                                    sd = numeric(0), ci_lo = numeric(0), ci_hi = numeric(0)))
  out <- do.call(rbind, lapply(split(tab[[metric]], tab$condition), function(v) {
    n <- length(v); m <- mean(v); s <- if (n > 1) stats::sd(v) else NA_real_
    h <- if (n > 1) stats::qt(0.975, n - 1) * s / sqrt(n) else NA_real_   # CI half-width
    data.frame(n_cells = n, mean = m, sd = s, ci_lo = m - h, ci_hi = m + h)
  }))
  cbind(condition = rownames(out), out, row.names = NULL)
}
