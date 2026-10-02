## ============================================================================
## est_clean.R - acquisition timing and track cleaning, as plain functions
## ----------------------------------------------------------------------------
## The Tracking tab's cycle-time estimate and gap splitting, outside the Shiny
## server, so the live view, Plot Data's rechunk and the per-cell comparison
## (est_percell.R) run the SAME code. The server reactives are thin wrappers.
##
##   half_sample_mode() ..... bin-free mode of a sample (Robertson-Cryer / Bickel)
##   cycle_cluster() ........ cycle time and SW from the integer-cycle clusters of dt
##   total_cycles() ......... cycle count of every interval (adaptive split input)
##   adaptive_gap_cycles() .. Vogler's median + SD + 1 split threshold, capped
##   split_gaps() ........... split tracks at long gaps, drop short fragments
##   trim_track_ends() ...... peel dark-gap / dim localizations off track ends
##   continuity_scores() .... dark-cycle-weighted continuity per track
##   gate_continuity() ...... drop tracks below a continuity score
##   filter_ranges() ........ per-localization range filters (EFO, ECO, CFR, ...)
##   filter_tid_length() .... track length filter
##   filter_duration() ...... track duration filter
##   filter_pipeline() ...... all of the above, in the live order
## ============================================================================

# ── half_sample_mode() ────────────────────────────────────────────────────────
# GOAL: the mode of a continuous sample without choosing a bin width. Keeps the
# densest half of the sorted sample until three or fewer points remain.
half_sample_mode <- function(x) {
  x <- sort(x[is.finite(x)]); n <- length(x)
  if (n == 0L) return(NA_real_)
  while (n > 3L) {
    h <- (n + 1L) %/% 2L
    w <- x[h:n] - x[1:(n - h + 1L)]                   # width of every half-sample window
    j <- which.min(w)
    x <- x[j:(j + h - 1L)]; n <- length(x)
  }
  if (n == 1L) return(x[1])
  if (n == 2L) return(mean(x))
  d1 <- x[2] - x[1]; d2 <- x[3] - x[2]
  if (d1 <= d2) mean(x[1:2]) else mean(x[2:3])
}

# ── cycle_cluster() ───────────────────────────────────────────────────────────
# GOAL: the acquisition cycle time from the data alone. Consecutive intervals
# fall into clusters at integer numbers of cycles; the median SPACING between
# adjacent cluster centres is the cycle, whichever cluster is most populated, so
# it holds even when most intervals span 2-3 cycles.
# tid, tim : track id and time (ms) of every localization (raw, before cleaning)
# RETURNS list(cycle, sw, n_clusters, centers, modal, frac_one_cycle), or NULL
#   when there are fewer than 50 positive intervals.
cycle_cluster <- function(tid, tim) {
  df  <- data.frame(tid = tid, tim = tim)
  df  <- df[order(df$tid, df$tim), ]
  dts <- unlist(tapply(df$tim, df$tid, function(t) diff(sort(t))), use.names = FALSE)
  dts <- dts[is.finite(dts) & dts > 0]
  if (length(dts) < 50L) return(NULL)
  dts <- sort(dts)
  modal <- half_sample_mode(dts)
  rough <- as.numeric(quantile(dts, 0.02))            # rough smallest-cycle scale
  if (!is.finite(rough) || rough <= 0) return(NULL)
  thr  <- 0.4 * rough                                  # gap threshold for clustering
  grp  <- cumsum(c(0, diff(dts) > thr))
  cl   <- split(dts, grp)
  centers <- vapply(cl, mean, numeric(1))
  counts  <- vapply(cl, length, integer(1))
  keep    <- counts >= max(5L, as.integer(0.005 * length(dts)))
  centers <- sort(centers[keep])
  if (length(centers) < 2L) {
    cyc <- if (length(centers) == 1L) centers[1] else rough
    return(list(cycle = cyc, sw = 0, n_clusters = length(centers),
                centers = centers, modal = modal, frac_one_cycle = NA_real_))
  }
  cyc <- stats::median(diff(centers))                  # spacing = atomic cycle
  sw  <- max(0, centers[1] - cyc)                      # intercept (one-time offset)
  # fraction of intervals that are a single (atomic) cycle
  n_cyc_each <- pmax(1L, as.integer(round((dts - sw) / cyc)))
  list(cycle = cyc, sw = sw, n_clusters = length(centers), centers = centers,
       modal = modal, frac_one_cycle = mean(n_cyc_each == 1L))
}

# ── total_cycles() ────────────────────────────────────────────────────────────
# GOAL: the number of acquisition cycles each interval spans, round((dt - SW)/cycle).
# min_len : when given, only tracks with at least this many localizations count,
#           as Vogler estimates his split threshold after filter_by_length().
# RETURNS integer vector of cycle counts >= 1.
total_cycles <- function(tid, tim, cyc, sw = 0, min_len = NULL) {
  if (!length(tid) || !is.finite(cyc) || cyc <= 0) return(integer(0))
  df <- data.frame(tid = tid, tim = tim)
  df <- df[order(df$tid, df$tim), ]
  if (!is.null(min_len)) {
    tl <- table(df$tid); keep_tid <- names(tl)[tl >= min_len]
    df <- df[df$tid %in% keep_tid, , drop = FALSE]
    if (nrow(df) == 0) return(integer(0))
  }
  dts <- unlist(lapply(split(df$tim, df$tid),
                       function(t) if (length(t) > 1L) diff(t) else numeric(0)),
                use.names = FALSE)
  tot <- round((dts - sw) / cyc)
  tot <- tot[is.finite(tot) & tot >= 1]
  as.integer(tot)
}

# ── adaptive_gap_cycles() ─────────────────────────────────────────────────────
# GOAL: Vogler's adaptive split threshold in cycles, floor(median + SD + 1) of the
# cycle counts, capped by the hardware "particle lost" timeout (maxOffTime).
# RETURNS integer threshold, or NA when there are fewer than 2 intervals.
adaptive_gap_cycles <- function(tot, cyc, maxoff_ms = 3.0) {
  if (length(tot) < 2L) return(NA_integer_)
  med <- as.integer(stats::median(tot))
  sdv <- stats::sd(tot); if (!is.finite(sdv)) sdv <- 0
  med_th <- as.integer(floor(med + sdv + 1))
  if (!is.finite(maxoff_ms)) maxoff_ms <- 3.0
  high_th <- if (is.finite(cyc) && cyc > 0) as.integer(round(maxoff_ms / cyc)) else NA_integer_
  if (is.finite(high_th) && high_th >= 2L) min(med_th, high_th) else med_th
}

# ── split_gaps() ──────────────────────────────────────────────────────────────
# GOAL: split every track wherever consecutive localizations are further apart
# than gap_thresh_ms, naming the pieces "<tid>_s1", "<tid>_s2", ..., and drop
# pieces shorter than min_len localizations. Idempotent: an existing "_s<N>"
# suffix is stripped first, so splitting twice equals splitting once.
# Fully vectorised with an integer-coded radix sort (sorting character track ids
# directly was ~25x slower on full raw data).
split_gaps <- function(df, gap_thresh_ms, min_len = 1L) {
  if (!nrow(df) || !is.finite(gap_thresh_ms) || gap_thresh_ms <= 0) return(df)
  if (any(grepl("_s", df$tid, fixed = TRUE)))
    df$tid <- sub("_s\\d+$", "", as.character(df$tid), perl = TRUE)
  tf  <- match(df$tid, unique(df$tid))                # integer code per track
  ord <- order(tf, df$tim, method = "radix")
  df  <- df[ord, , drop = FALSE]; tf <- tf[ord]
  n   <- nrow(df); tid <- as.character(df$tid)
  newtrk <- c(TRUE, tf[-1] != tf[-n])                 # first row of each track
  dt     <- c(NA_real_, diff(df$tim)); dt[newtrk] <- NA_real_
  is_gap <- is.finite(dt) & (dt > gap_thresh_ms)      # split point at this row
  grp    <- cumsum(newtrk)                            # track id 1..K
  cg     <- cumsum(is_gap)
  seg    <- cg - cg[newtrk][grp] + 1L                 # 1-based segment # within track
  nseg   <- (as.integer(rowsum(as.integer(is_gap), grp)) + 1L)[grp]  # #segments in track
  new_tid <- tid; spl <- nseg > 1L
  if (any(spl)) new_tid[spl] <- paste0(tid[spl], "_s", seg[spl])
  if (min_len > 1L) {                                 # drop segments shorter than min_len
    u    <- unique(new_tid); mi <- match(new_tid, u)
    keep <- tabulate(mi, length(u))[mi] >= min_len
    df   <- df[keep, , drop = FALSE]; new_tid <- new_tid[keep]
  }
  if (!nrow(df)) return(df[0, , drop = FALSE])
  df$tid <- new_tid
  df
}

## ----------------------------------------------------------------------------
## The rest of the cleaning and filter pipeline. The live Tracking pipeline, the
## Plot Data rechunk and the per-cell comparison all call these, in this order:
##   split_gaps -> trim_track_ends -> gate_continuity -> filter_ranges ->
##   filter_tid_length -> filter_duration        (filter_pipeline() runs them all)
## ----------------------------------------------------------------------------

# ── trim_track_ends() ─────────────────────────────────────────────────────────
# GOAL: peel degraded localizations off each end of every track: those reached
# across a dark gap of at least dark_cycles cycles, or dimmer than efo_frac times
# the track's median EFO, working inward until a solid localization is reached.
# Stops at floor_len so a track never drops below the gap split's minimum length.
trim_track_ends <- function(df, t_c, use_dark = TRUE, dark_cycles = 3L, use_efo = FALSE,
                            efo_frac = 0.5, floor_len = 1L) {
  if (!nrow(df) || !"tim" %in% names(df)) return(df)
  use_efo  <- isTRUE(use_efo) && "efo" %in% names(df)
  use_dark <- isTRUE(use_dark)
  if (!use_dark && !use_efo) return(df)
  dark_cyc <- max(2L, as.integer(dark_cycles))
  dark_gap <- if (use_dark && is.finite(t_c) && t_c > 0) dark_cyc * t_c else Inf
  tf  <- match(df$tid, unique(df$tid))                  # integer track codes (fast sort/group)
  ord <- order(tf, df$tim, method = "radix")
  df  <- df[ord, , drop = FALSE]
  sp  <- split(seq_len(nrow(df)), tf[ord])              # group by integer, not character
  keep <- logical(nrow(df))
  for (idx in sp) {
    m <- length(idx)
    if (m < 2L) { keep[idx] <- TRUE; next }
    d <- diff(df$tim[idx])                              # d[j] = gap idx[j] -> idx[j+1]
    efo_thr <- if (use_efo) efo_frac * stats::median(df$efo[idx], na.rm = TRUE) else -Inf
    efo_v   <- if (use_efo) df$efo[idx] else rep(Inf, m)
    lo <- 1L
    while (lo < m && (m - lo + 1L) > floor_len &&
           ((use_efo && is.finite(efo_v[lo]) && efo_v[lo] < efo_thr) ||
            (use_dark && d[lo]      >= dark_gap))) lo <- lo + 1L
    hi <- m
    while (hi > lo && (hi - lo + 1L) > floor_len &&
           ((use_efo && is.finite(efo_v[hi]) && efo_v[hi] < efo_thr) ||
            (use_dark && d[hi - 1L] >= dark_gap))) hi <- hi - 1L
    if (hi >= lo) keep[idx[lo:hi]] <- TRUE
  }
  df[keep, , drop = FALSE]
}

# ── continuity_scores() ───────────────────────────────────────────────────────
# GOAL: one continuity score per track in [0, 1]. Each step spanning k cycles
# (k = 1 integrated, k >= 2 = k - 1 dark cycles) is weighted 2^(-(k-1)/half_life);
# the score is the mean weight. half_life 0 or NA gives the plain integrated fraction.
# RETURNS a named vector (tid -> score), or NULL.
continuity_scores <- function(df, t_c, half_life = 3) {
  if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df)) ||
      !is.finite(t_c) || t_c <= 0) return(NULL)
  utid <- unique(df$tid); tf <- match(df$tid, utid)     # integer track codes
  ord  <- order(tf, df$tim, method = "radix")           # radix sort (fast on integers)
  tf   <- tf[ord]; tim <- df$tim[ord]; n <- length(tf)
  step <- c(FALSE, tf[-1] == tf[-n])
  dt   <- c(NA_real_, diff(tim)); dt[!step] <- NA_real_
  k    <- pmax(1L, as.integer(round(dt / t_c)))
  w    <- if (is.finite(half_life) && half_life > 0) 2^(-(k - 1) / half_life)
          else as.numeric(k == 1L)                      # half-life 0 -> integrated fraction
  sc   <- tapply(w[step], tf[step], mean)               # by integer code
  stats::setNames(as.numeric(sc), utid[as.integer(names(sc))])   # map codes -> tids
}

# ── gate_continuity() ─────────────────────────────────────────────────────────
# GOAL: keep tracks whose continuity score is at least min_score. Single-
# localization tracks have no score and are kept.
gate_continuity <- function(df, t_c, min_score = 0.6, half_life = 3) {
  if (!nrow(df) || !"tim" %in% names(df) || !is.finite(t_c) || t_c <= 0) return(df)
  thr <- max(0, min(1, as.numeric(min_score)))
  sc  <- continuity_scores(df, t_c, half_life)
  if (is.null(sc)) return(df)
  keep <- names(sc)[is.finite(sc) & sc >= thr]
  keep <- union(keep, setdiff(unique(df$tid), names(sc)))   # keep unjudged single-loc tracks
  df[df$tid %in% keep, , drop = FALSE]
}

# ── filter_ranges() ───────────────────────────────────────────────────────────
# GOAL: per-localization range filters. ranges: named list, column -> c(lo, hi)
# or NULL (no filter). Missing values pass. A range with no overlap at all with
# the column's values is treated as stale and passes everything.
filter_ranges <- function(df, ranges) {
  for (nm in names(ranges)) {
    rng <- ranges[[nm]]
    if (is.null(rng) || length(rng) != 2L || !nm %in% names(df) || !nrow(df)) next
    v <- df[[nm]]; v <- v[is.finite(v)]
    if (length(v) && (rng[1] > max(v) || rng[2] < min(v))) next   # stale guard
    x <- df[[nm]]
    df <- df[is.na(x) | (x >= rng[1] & x <= rng[2]), , drop = FALSE]
  }
  df
}

# ── filter_tid_length() ───────────────────────────────────────────────────────
# GOAL: keep tracks whose localization count lies in rng. A range with no
# overlap with the data's lengths is treated as stale and passes everything.
filter_tid_length <- function(df, rng) {
  if (is.null(rng) || length(rng) != 2L || !nrow(df) || !"tid" %in% names(df)) return(df)
  tl <- table(df$tid); lv <- as.numeric(tl)
  if (rng[1] > max(lv) || rng[2] < min(lv)) return(df)            # stale guard
  keep <- names(tl)[lv >= rng[1] & lv <= rng[2]]
  df[df$tid %in% keep, , drop = FALSE]
}

# ── filter_duration() ─────────────────────────────────────────────────────────
# GOAL: keep tracks whose duration max(tim) - min(tim) (ms) lies in rng. A range
# covering every duration is a no-op.
filter_duration <- function(df, rng) {
  if (is.null(rng) || length(rng) != 2L || !nrow(df) || !"tim" %in% names(df)) return(df)
  dur <- tapply(df$tim, df$tid, function(t) max(t) - min(t))
  dv  <- as.numeric(dur)
  if (!length(dv) || (rng[1] <= min(dv) && rng[2] >= max(dv))) return(df)
  keep <- names(dur)[dv >= rng[1] & dv <= rng[2]]
  df[df$tid %in% keep, , drop = FALSE]
}

# ── filter_pipeline() ─────────────────────────────────────────────────────────
# GOAL: the whole cleaning + filter pipeline on one dataset, in the live order.
# df : localizations with tim in ms
# p  : list(t_c, sw,
#           split   = list(enable, adaptive, gap_cycles, min_len, maxoff_ms),
#           trim    = list(enable, dark, dark_cycles, efo, efo_frac),
#           gate    = list(enable, min_score, half_life),
#           ranges  = named list column -> c(lo, hi) or NULL,
#           tid_len = c(lo, hi) or NULL, dur = c(lo, hi) or NULL)
# The spatial crop is not part of it: it belongs to one field of view.
filter_pipeline <- function(df, p) {
  t_c <- p$t_c; s <- p$split %||% list(); tr <- p$trim %||% list(); g <- p$gate %||% list()
  min_len <- max(1L, as.integer(s$min_len %||% 50L))
  if (isTRUE(s$enable) && is.finite(t_c) && t_c > 0) {
    gap_cyc <- if (isTRUE(s$adaptive))
                 adaptive_gap_cycles(total_cycles(df$tid, df$tim, t_c, p$sw %||% 0, min_len),
                                     t_c, s$maxoff_ms %||% 3.0)
               else suppressWarnings(as.numeric(s$gap_cycles))
    if (length(gap_cyc) && is.finite(gap_cyc) && gap_cyc > 0) df <- split_gaps(df, gap_cyc * t_c, min_len)
  }
  if (isTRUE(tr$enable))
    df <- trim_track_ends(df, t_c, use_dark = isTRUE(tr$dark), dark_cycles = tr$dark_cycles %||% 3L,
                          use_efo = isTRUE(tr$efo), efo_frac = tr$efo_frac %||% 0.5,
                          floor_len = if (isTRUE(s$enable)) min_len else 1L)
  if (isTRUE(g$enable))
    df <- gate_continuity(df, t_c, g$min_score %||% 0.6, g$half_life %||% 3)
  df <- filter_ranges(df, p$ranges %||% list())
  df <- filter_tid_length(df, p$tid_len)
  filter_duration(df, p$dur)
}
