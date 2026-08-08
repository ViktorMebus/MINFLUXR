## ─────────────────────────────────────────────────────────────────────────────
## helpers_tracking.R  -  tracking-mode-only helper functions
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Small, self-contained utilities used ONLY by the Tracking tab
##   (ui_tracking.R + server_tracking.R). Nothing in this file is reactive:
##   every function is a plain input -> output transform, so each one can be
##   tested on its own in the console without launching the app.
##
## WHAT LIVES HERE (grouped by job)
##   1. Colour          tid_color_map()          stable colour per trace ID
##   2. On-disk chunking split_by_tid_to_dir()   split a big dataset into
##                                               .qs/.fst chunks + manifest.csv
##   3. Data prep       prep_tracking()          validate / sort / zero the time axis
##   4. Plot helpers    decimate_per_tid()       thin tracks so plotly stays fast
##                      nearest_idx()            nearest time index for scrubbing
##
## WHY THE CHUNKING EXISTS
##   MINFLUX tracking files are large (often ~1e6 localizations). The app does
##   not keep the raw upload in memory: it writes per-TID chunks to disk once,
##   records them in manifest.csv, and the server reads chunks back on demand.
##
## CONVENTIONS USED THROUGHOUT
##   - jsonlite is NOT library()'d - always call jsonlite:: explicitly, because
##     its `validate` would otherwise shadow shiny::validate.
##   - Every input/output ID in the tracking tab is prefixed `trk_`, to avoid
##     collisions with the SMLM tab's IDs.
##   - Coordinates are in nanometres. `tim` is in MILLISECONDS by the time it
##     reaches these helpers (server_tracking.R multiplies the on-disk seconds
##     by 1000 right after loading).
##   - "TID" = trace ID = one tracked molecule = one trajectory.
##
## REMOVED CODE
##   convert_to_qs(), split_by_tid(), na_lines() and even_idx() used to live here
##   but had no caller anywhere in the app. They were moved verbatim to
##   deadcode_tracking.R (2026-07-20) - see that file to restore them.
## ─────────────────────────────────────────────────────────────────────────────


# ── tid_color_map() ──────────────────────────────────────────────────────────
# GOAL: give every trace ID (TID) its own visually distinct, STABLE colour.
#
# Why stability matters: the same TID must keep the same colour across
# re-renders, filter changes and playback frames. Otherwise tracks appear to
# swap identity mid-animation. That is the whole reason the shuffle is seeded.
#
# tids    : vector of trace IDs (any type - coerced to character)
# palette : an RColorBrewer palette name; silently falls back to "Set3"
# seed    : fixes the colour shuffle so the mapping is reproducible
# returns : named character vector, tid -> "#RRGGBB"
#
# SIDE EFFECT: set.seed() below resets R's GLOBAL random state. Harmless here
# (colours only), but do not call this in the middle of a simulation that
# depends on the RNG stream.
tid_color_map <- function(tids, palette = "Set3", seed = 123) {
  lev <- sort(unique(as.character(tids)))   # unique TIDs in a deterministic order
  n   <- length(lev)                        # how many distinct colours we need
  info    <- RColorBrewer::brewer.pal.info  # lookup table: palette -> max colours
  if (!palette %in% rownames(info)) palette <- "Set3"   # unknown name -> safe default
  maxcols <- info[palette, "maxcolors"]     # largest size this palette supports
  # brewer.pal() errors below 3 colours, and errors above the palette's max -
  # so clamp the request into [3, maxcols].
  base    <- RColorBrewer::brewer.pal(max(3, min(maxcols, n)), palette)
  # If we need MORE colours than the palette holds, interpolate a continuous
  # ramp through them; otherwise just take the first n.
  cols    <- if (n > length(base)) grDevices::colorRampPalette(base)(n) else base[seq_len(n)]
  # Shuffle so that ADJACENT TIDs get contrasting hues. In palette order,
  # neighbouring entries are often near-identical, which makes overlapping
  # tracks hard to tell apart. Seeded => same shuffle every time.
  set.seed(seed); cols <- sample(cols, n)
  stats::setNames(as.character(cols), lev)  # -> lookup table keyed by TID
}


# ── prep_tracking() ──────────────────────────────────────────────────────────
# GOAL: put a raw tracking data.frame into the canonical shape that the rest of
# the tracking code assumes. Everything downstream (MSD, playback, filtering)
# relies on these guarantees:
#   - the four required columns exist          (hard stop if not)
#   - no non-finite coordinates or times       (rows dropped)
#   - rows sorted by (tid, time)               (step differences assume this!)
#   - a per-track time axis starting at zero   (column `tim0`)
#
# Why `tim0`: different tracks start at different wall-clock times. To draw or
# animate them on a shared axis you need "time since THIS track began", not
# absolute time.
#
# df                : raw data.frame, must contain x, y, tid, tim
# normalize_per_tid : TRUE -> tim0 is relative to each track's first
#                     localization; FALSE -> tim0 is just a copy of tim
# returns           : same rows (minus non-finite), reordered, with tim0 added
prep_tracking <- function(df, normalize_per_tid = TRUE) {
  stopifnot(all(c("x","y","tid","tim") %in% names(df)))   # fail loudly, early
  df %>%
    # Put the key columns first for readability; keep every other column as-is.
    dplyr::select(tid, x, y, tim, dplyr::everything()) %>%
    dplyr::filter(is.finite(x), is.finite(y), is.finite(tim)) %>%   # drop NA/Inf rows
    dplyr::arrange(tid, tim) %>%                          # REQUIRED: diff() assumes time order
    dplyr::group_by(tid) %>%                              # zero the clock per track...
    dplyr::mutate(tim0 = if (normalize_per_tid) tim - dplyr::first(tim) else tim) %>%
    dplyr::ungroup()                                      # ...and always ungroup before returning
}


# ── .geom_ops() ───────────────────────────────────────────────────────────────
# Shared geometry engine for the free-diffusion nulls. Returns the operators for one
# shape so that every null (net displacement, MSD, ...) steps molecules through the
# SAME physics - if these ever diverged the two nulls would disagree for reasons that
# had nothing to do with the data.
#
# geom  : "sphere" | "rod" | "ellipsoid" | "membrane" | "free"
# params: list(R, r, Lcyl, a, b, c, mtilt) - only the relevant entries are read
# mode  : "surface" = confined to the boundary (2 DOF on a curved sheet)
#         "volume"  = free 3D motion inside, reflecting off the boundary
# returns list(one_step(P, sd), rand_start(m), is_flat, is_vol)
.geom_ops <- function(geom, params = list(), mode = c("surface", "volume")) {
  mode <- match.arg(mode)
  `%|NA|%` <- function(a, b) if (!is.null(a) && is.finite(a)) a else b
  R    <- params$R    %|NA|% 500;  rr   <- params$r %|NA|% 400;  Lcyl <- params$Lcyl %|NA|% 2000
  abc  <- if (identical(geom, "sphere")) c(R, R, R)
          else c(params$a %|NA|% 500, params$b %|NA|% 500, params$c %|NA|% 900)
  mtilt <- (params$mtilt %|NA|% 0) * pi / 180
  is_ell  <- geom %in% c("sphere", "ellipsoid")
  is_flat <- geom %in% c("free", "membrane")
  n_hat <- c(sin(mtilt), 0, cos(mtilt)); e1 <- c(cos(mtilt), 0, -sin(mtilt)); e2 <- c(0, 1, 0)

  dotrows <- function(A, B) rowSums(A * B)
  geo_normal <- function(P) {
    if (is_ell) { N <- sweep(P, 2, abc^2, "/"); N / sqrt(dotrows(N, N)) }
    else {                                                     # rod: body vs hemispherical caps
      a <- P[, 1]; body <- a >= 0 & a <= Lcyl
      N <- matrix(0, nrow(P), 3); N[body, 2:3] <- P[body, 2:3]
      cl <- a < 0; cr <- a > Lcyl
      if (any(cl)) { N[cl, ] <- P[cl, , drop = FALSE] }
      if (any(cr)) { N[cr, ] <- P[cr, , drop = FALSE]; N[cr, 1] <- P[cr, 1] - Lcyl }
      m <- sqrt(dotrows(N, N)); m[m == 0] <- 1; N / m
    }
  }
  geo_project <- function(P) {
    if (is_ell) { U <- sweep(P, 2, abc, "/"); sweep(U / sqrt(dotrows(U, U)), 2, abc, "*") }
    else {
      a <- P[, 1]; body <- a >= 0 & a <= Lcyl; Q <- P
      rho <- sqrt(P[, 2]^2 + P[, 3]^2); rho[rho == 0] <- 1e-9
      Q[body, 2] <- P[body, 2] * rr / rho[body]; Q[body, 3] <- P[body, 3] * rr / rho[body]
      for (side in list(c(m = 0), c(m = Lcyl))) {
        cx <- side[["m"]]; msk <- if (cx == 0) a < 0 else a > Lcyl
        if (any(msk)) { d <- P[msk, , drop = FALSE]; d[, 1] <- d[, 1] - cx
          dd <- sqrt(dotrows(d, d)); dd[dd == 0] <- 1e-9; d <- d * rr / dd
          d[, 1] <- d[, 1] + cx; Q[msk, ] <- d }
      }
      Q
    }
  }
  geo_rand <- function(m) {                                    # uniform ON the boundary
    if (is_ell) geo_project(sweep(matrix(stats::rnorm(3 * m), m, 3), 2, abc, "*"))
    else { s <- stats::runif(m, 0, Lcyl); phi <- stats::runif(m, 0, 2 * pi)
           cbind(s, rr * cos(phi), rr * sin(phi)) }
  }
  geo_reflect <- function(P) {                                 # back INSIDE the shape
    if (is_ell) {
      U <- sweep(P, 2, abc, "/"); nu <- sqrt(dotrows(U, U))
      f <- ifelse(nu > 1, (2 - nu) / nu, 1); f[!is.finite(f)] <- 1; f <- pmax(f, 0)
      sweep(U * f, 2, abc, "*")
    } else {
      a <- P[, 1]; Q <- P; body <- a >= 0 & a <= Lcyl
      rho <- sqrt(P[, 2]^2 + P[, 3]^2)
      fb <- ifelse(rho > rr, (2 * rr - rho) / rho, 1); fb[!is.finite(fb)] <- 1; fb <- pmax(fb, 0)
      Q[body, 2] <- P[body, 2] * fb[body]; Q[body, 3] <- P[body, 3] * fb[body]
      for (cx in c(0, Lcyl)) {
        msk <- if (cx == 0) a < 0 else a > Lcyl
        if (any(msk)) {
          d <- P[msk, , drop = FALSE]; d[, 1] <- d[, 1] - cx
          dd <- sqrt(dotrows(d, d))
          f2 <- ifelse(dd > rr, (2 * rr - dd) / dd, 1); f2[!is.finite(f2)] <- 1; f2 <- pmax(f2, 0)
          d <- d * f2; d[, 1] <- d[, 1] + cx; Q[msk, ] <- d
        }
      }
      Q
    }
  }
  geo_rand_in <- function(m) {                                 # uniform WITHIN the volume
    if (is_ell) {
      d <- matrix(stats::rnorm(3 * m), m, 3); d <- d / sqrt(dotrows(d, d))
      sweep(d * stats::runif(m)^(1 / 3), 2, abc, "*")
    } else {
      out <- matrix(NA_real_, m, 3); need <- seq_len(m)
      while (length(need)) {
        k <- length(need)
        p <- cbind(stats::runif(k, -rr, Lcyl + rr), stats::runif(k, -rr, rr),
                   stats::runif(k, -rr, rr))
        a <- p[, 1]
        ins <- ifelse(a >= 0 & a <= Lcyl, p[, 2]^2 + p[, 3]^2 <= rr^2,
               ifelse(a < 0, dotrows(p, p) <= rr^2,
                      (p[, 1] - Lcyl)^2 + p[, 2]^2 + p[, 3]^2 <= rr^2))
        if (any(ins)) { out[need[ins], ] <- p[ins, , drop = FALSE]; need <- need[!ins] }
      }
      out
    }
  }
  is_vol <- identical(mode, "volume") && !is_flat
  one_step <- function(P, sd) {
    G <- matrix(stats::rnorm(nrow(P) * 3, 0, sd), nrow(P), 3)
    if (is_flat) {
      if (identical(geom, "free")) P + G
      else P + outer(G[, 1], e1) + outer(G[, 2], e2)        # membrane: 2 DOF in the tilted plane
    } else if (is_vol) geo_reflect(P + G)
    else { N <- geo_normal(P); geo_project(P + (G - N * dotrows(G, N))) }
  }
  rand_start <- function(m) {
    if (is_flat) matrix(0, m, 3) else if (is_vol) geo_rand_in(m) else geo_rand(m)
  }
  list(one_step = one_step, rand_start = rand_start, is_flat = is_flat, is_vol = is_vol)
}

# ── sim_free_geom_rnet() ──────────────────────────────────────────────────────
# GOAL: the GEOMETRY-AWARE free-diffusion null for the displacement test (Step 2), for
# ANY of the simulator geometries. Simulate FREE Brownian diffusion on the chosen shape
# and return each track's NET start→end displacement in the 2D image projection - so
# "how far a molecule gets" is judged against free diffusion ON THAT GEOMETRY, not an
# infinite plane. Closed surfaces (sphere/rod/ellipsoid) cap displacement at ~their size;
# an unbounded plane (membrane/free) does not, so those reproduce the infinite-plane null.
#
# PHYSICS: tangent Brownian step on the surface then project back to it (the same operators
# as the Diffusion Simulator), vectorised across all M tracks (one matrix op per cycle).
# 2D image projection drops z (the long axis lies in-plane, the usual mounting).
#
# PROJECTION MATCHING (match_projection = TRUE): a measured 2D D from a curved surface is
# the PROJECTED D, which is < the true surface D (the azimuthal component is foreshortened;
# e.g. exactly 3/4 for a cylinder). So we EMPIRICALLY calibrate that factor for the chosen
# geometry - one short single-step pass - and drive the null at D_surf = D_meas / factor,
# so both nulls agree at short T and the closed-surface band diverges only where finite
# size + curvature bite. Set FALSE to feed D straight through.
#
# dur_ms  : track durations (ms) to match (resample from real data)
# D_nm2ms : measured per-axis diffusion coefficient (nm²/ms)
# geom    : "sphere" | "rod" | "ellipsoid" | "membrane" | "free"
# params  : list(R, r, Lcyl, a, b, c, mtilt) - only the relevant ones are read
# returns : numeric vector (length M) of net displacements R_net (nm), + attr "D_surf"
sim_free_geom_rnet <- function(dur_ms, D_nm2ms, geom = "rod", params = list(),
                               cyc_ms = 0.5, proj2d = TRUE, match_projection = TRUE,
                               mode = c("surface", "volume")) {
  mode <- match.arg(mode)
  M <- length(dur_ms)
  if (M == 0L || !is.finite(D_nm2ms) || D_nm2ms <= 0) return(numeric(0))
  n_cyc <- pmax(2L, as.integer(round(dur_ms / cyc_ms))); nmax <- max(n_cyc)
  ops <- .geom_ops(geom, params, mode)

  # Observed axes: 2 when the data is a 2D image projection (drop z), 3 when tracked in 3D.
  nobs   <- if (proj2d) 2L else 3L
  to_obs <- function(P) P[, seq_len(nobs), drop = FALSE]

  # ── projection-factor calibration → surface D that reproduces the MEASURED per-axis D ──
  # f is the per-axis MSD ratio of the observed motion to the input surface motion. For a
  # 2D projection of a cylinder f = 3/4; for a curved surface tracked in 3D f = 2/3 (2 DOF
  # spread over 3 axes). Dividing by nobs makes it the per-axis factor in either case.
  D_surf <- D_nm2ms
  if (match_projection && !identical(geom, "free")) {
    ncal <- 4000L; sdc <- sqrt(2 * D_nm2ms * cyc_ms)
    P0 <- ops$rand_start(ncal)
    d2 <- to_obs(ops$one_step(P0, sdc)) - to_obs(P0)
    f2 <- mean(rowSums(d2^2)) / (nobs * sdc^2)
    if (is.finite(f2) && f2 > 0.05) D_surf <- D_nm2ms / f2
  }
  sd <- sqrt(2 * D_surf * cyc_ms)

  # ── main pass ─────────────────────────────────────────────────────────────────
  P  <- ops$rand_start(M)
  S2 <- to_obs(P); E2 <- S2
  for (i in 2:nmax) {
    alive <- i <= n_cyc
    Pn <- ops$one_step(P, sd); P[alive, ] <- Pn[alive, ]
    ended <- i == n_cyc
    if (any(ended)) E2[ended, ] <- to_obs(P)[ended, ]
  }
  out <- sqrt(rowSums((E2 - S2)^2))
  attr(out, "D_surf") <- D_surf
  out
}


# ── nearest_idx() ────────────────────────────────────────────────────────────
# GOAL: given a SORTED vector of times `tt`, return the index of the entry
# closest to the query time `t`. Used to map a playback scrub position onto an
# actual localization.
#
# tt      : sorted numeric vector of times (ASSUMED sorted - not checked)
# t       : the query time
# returns : integer index into tt, or NA_integer_ if tt is empty
nearest_idx <- function(tt, t) {
  if (!length(tt)) return(NA_integer_)
  # findInterval gives the index of the last element <= t. That neighbour is
  # not necessarily the CLOSEST one, so we compare it against the next entry.
  i <- base::findInterval(t, tt)
  if (i <= 0) return(1L)                     # t is before the first time -> clamp to start
  if (i >= length(tt)) return(length(tt))    # t is after the last time  -> clamp to end
  # Genuinely between two samples: pick whichever is nearer in time.
  if (abs(tt[i+1] - t) < abs(tt[i] - t)) i + 1L else i
}


# ── decimate_per_tid() ───────────────────────────────────────────────────────
# GOAL: cap how many points each track contributes to a plot, so plotly stays
# responsive on datasets with ~1e6 localizations.
#
# KEY DESIGN CHOICE: thinning happens PER TRACK, not globally. A global sample
# would erase short tracks entirely whenever a few long tracks dominate the
# dataset. Per-track thinning means every track keeps its overall shape.
#
# The points kept are EVENLY SPACED along the track (not random), so the
# thinned path still traces the original trajectory.
#
# NOTE: this is for DISPLAY only. Never feed a decimated frame into MSD/CVE or
# any other quantitative estimator - dropping points changes the step statistics.
#
# df         : data.frame with tid, x, y, [z], and the time column `tcol`
# max_pts    : maximum points to keep per track
# tcol       : name of the time column to sort/keep (default "tim0")
# extra_cols : extra columns to carry through untouched (e.g. per-point colour)
# returns    : data.frame(tid, x, y, t, [z], <extra_cols>) - note the time
#              column is RENAMED to `t` in the output
decimate_per_tid <- function(df, max_pts = 400, tcol = "tim0", extra_cols = character()) {
  has_z <- "z" %in% names(df)                        # 3D data is optional
  extra_cols <- intersect(extra_cols, names(df))     # carry these through unchanged (e.g. per-point colour)
  # Accumulators. Built as plain vectors and combined once at the end.
  outx <- outy <- outt <- outid <- NULL
  outz <- NULL
  out_extra <- stats::setNames(vector("list", length(extra_cols)), extra_cols)
  # Keep only the columns we actually need, then sort by (tid, time) so the
  # even-spacing selection below walks each track in chronological order.
  sel_cols <- c("tid","x","y", if (has_z) "z", tcol, extra_cols)
  df_ord <- df[order(df$tid, df[[tcol]]), sel_cols, drop = FALSE]
  sp <- split(df_ord, df_ord$tid)                    # one data.frame per track
  for (nm in names(sp)) {
    tb <- sp[[nm]]; n <- nrow(tb)
    # Short track -> keep everything. Long track -> keep max_pts evenly spaced
    # indices. unique() guards against round() producing duplicates when the
    # track is only slightly longer than max_pts.
    keep  <- if (n <= max_pts) seq_len(n) else unique(round(seq(1, n, length.out = max_pts)))
    outid <- c(outid, rep(nm, length(keep)))         # repeat the TID for each kept row
    outx  <- c(outx,  tb$x[keep])
    outy  <- c(outy,  tb$y[keep])
    if (has_z) outz <- c(outz, tb$z[keep])
    outt  <- c(outt,  tb[[tcol]][keep])
    for (ec in extra_cols) out_extra[[ec]] <- c(out_extra[[ec]], tb[[ec]][keep])
  }
  # Assemble the result. The time column is deliberately renamed to `t` so
  # callers get one predictable name regardless of what `tcol` was.
  out <- data.frame(tid = outid, x = outx, y = outy, t = outt,
                    row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  if (has_z) out$z <- outz
  for (ec in extra_cols) out[[ec]] <- out_extra[[ec]]
  out
}


# ── split_by_tid_to_dir() ────────────────────────────────────────────────────
# GOAL: the app's main "make a large upload manageable" step. Takes a dataset
# ALREADY IN MEMORY, splits it into chunk files on disk (grouped by whole
# TIDs), and writes a manifest.csv index next to them.
#
# This is what lets the app handle ~1e6-localization files: after chunking, the
# server only ever reads back the chunks it needs, and `manifest.csv` is the
# lookup that says which file holds which TIDs.
#
# INVARIANT: a track is NEVER split across two chunk files. Chunk boundaries
# always fall between whole TIDs, so any per-track computation can be done by
# loading a single chunk.
#
# df             : data.frame with at least tid and tim
# out_dir        : directory to create and write into
# tids_per_chunk : how many whole TIDs per chunk file (global default)
# fmt            : "qs" (smaller on disk) or "fst" (faster partial reads)
# returns        : invisibly, the path to the written manifest.csv
split_by_tid_to_dir <- function(df, out_dir,
                                tids_per_chunk = TIDS_PER_CHUNK_DEFAULT,
                                fmt = CHUNK_FORMAT) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  df$tid <- as.character(df$tid)                       # TIDs are character everywhere in-app
  df <- df[order(df$tid, df$tim), , drop = FALSE]      # chronological within each track
  tids <- sort(unique(df$tid))
  if (!length(tids)) stop("No TIDs to split.")         # nothing to write -> fail loudly
  # Group whole TIDs into blocks. max(1L, ...) guards against a 0/NA setting
  # that would otherwise make ceiling() produce Inf groups.
  chunks <- split(tids, ceiling(seq_along(tids) / max(1L, as.integer(tids_per_chunk))))

  # The manifest is the index the server reads later to find its data.
  manifest <- data.frame(chunk  = seq_along(chunks),
                         n_tids = lengths(chunks),
                         file   = NA_character_)

  # Chunking a big file is slow, so report progress in the Shiny UI.
  withProgress(message = "Chunking large dataset...", value = 0, {
    for (i in seq_along(chunks)) {
      incProgress(1/length(chunks), detail = paste("Chunk", i, "of", length(chunks)))
      sub <- df[df$tid %in% chunks[[i]], , drop = FALSE]   # all rows for this chunk's TIDs
      f <- file.path(out_dir, sprintf("chunk_%03d.%s", i, fmt))  # zero-padded => sorts correctly
      if (identical(fmt, "qs")) {
        qs::qsave(sub, f, preset = "high")             # "high" = smaller file, slower write
      } else if (identical(fmt, "fst")) {
        fst::write_fst(sub, f, compress = 50)          # 50 = balanced size vs speed
      } else {
        stop("Unsupported chunk format: ", fmt)        # caught early rather than half-writing
      }
      manifest$file[i] <- f                            # record where this chunk landed
    }
  })
  write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)
  invisible(file.path(out_dir, "manifest.csv"))
}




# ── sim_free_ci_null() ────────────────────────────────────────────────────────
# GOAL: the free-diffusion null for the per-track MSD CONFINEMENT INDEX (CI).
#
# CI = (MSD_late/MSD_early) / (tau_late/tau_early) over a track's own MSD window.
# For ideal free Brownian motion CI = 1, but the OBSERVED free value is NOT 1: the
# localization-noise term 2*d*sigma^2 sits on every MSD point and inflates the EARLY
# window more than the late one, pushing CI below 1. How far below depends on
# sigma_loc, D, the cycle time, the fit window and the track length - so there is no
# universal threshold and the null has to be simulated at the data's own conditions.
#
# Free tracks are generated with the SAME irregular sampling as the data: each step
# spans k cycles drawn from the measured dark-gap distribution `gap_p`, so the null
# inherits the real integrated/dark mix rather than assuming an even clock.
#
# WHY THE NULL IS ALSO MATCHED ON D. The free CI is governed by the dimensionless
# ratio  x = sigma^2 / (D * t_c) - how large the localization noise is relative to the
# motion in one cycle. A SLOW track has a bigger x, so the noise intercept lifts its
# early lags more and its free CI sits lower. Measured at N=200, sigma = 10 nm: the
# free median runs 0.325 at D = 20 nm^2/ms and 0.922 at D = 2000, i.e. a factor of
# ~2.8 across the D range. Calibrating every track against a single ensemble-D null
# would therefore call the slow half of a heterogeneous population "confined". At
# small sigma the dependence nearly vanishes (at sigma = 1.5 nm the same sweep moves
# the median only 0.917 -> 0.967), which is why the caller may collapse the D grid to
# one point when the noise is negligible.
#
# n_grid   : localization counts to simulate the null at (nearest-N lookup by caller)
# D_grid   : per-axis D values (nm^2/ms) to simulate at; length 1 = the old behaviour
# n_sim    : free tracks per (n, D) grid point
# t_c      : cycle time (ms)         sigma : localization sigma (nm)
# ndim     : 2 or 3
# gap_p    : probabilities of a step spanning 1,2,3,... cycles (measured from the data)
# min_lag/max_lag : the SAME cycle-lag fit window the real CI uses
# returns  : data.frame(n, D, lo05, med, hi95) - CI quantiles at each (n, D) pair
sim_free_ci_null <- function(n_grid, D_grid = 100, n_sim = 200L, t_c = 0.2, sigma = 10,
                             ndim = 2L, gap_p = c(0.55, 0.22, 0.09, 0.05, 0.04, 0.05),
                             min_lag = 4L, max_lag = 50L, progress = NULL) {
  gap_p <- gap_p[is.finite(gap_p) & gap_p > 0]
  if (!length(gap_p)) gap_p <- 1
  gap_p <- gap_p / sum(gap_p)
  kmax  <- length(gap_p)
  D_grid <- D_grid[is.finite(D_grid) & D_grid > 0]
  if (!length(D_grid)) D_grid <- 100

  # CI of one simulated free track, using cycle-index binning exactly as the app does.
  ci_one <- function(N, D) {
    k   <- sample.int(kmax, N - 1L, TRUE, prob = gap_p)
    cyc <- c(0L, cumsum(k))
    pos <- matrix(0, N, ndim)
    sd_step <- sqrt(2 * D * t_c * k)
    for (j in seq_len(ndim))
      pos[, j] <- cumsum(c(0, stats::rnorm(N - 1L, 0, sd_step)))
    pos <- pos + matrix(stats::rnorm(N * ndim, 0, sigma), N, ndim)
    # per-track MSD by integer cycle lag, capped at max_lag (bounds the inner loop)
    s <- numeric(max_lag); cnt <- integer(max_lag)
    for (i in seq_len(N - 1L)) {
      # Every step spans >= 1 cycle, so lag >= (j - i); j beyond i + max_lag can never
      # land in the window. That bounds the inner loop without dropping any pair.
      jj  <- seq.int(i + 1L, min(N, i + max_lag))
      lag <- cyc[jj] - cyc[i]
      ok  <- lag >= 1L & lag <= max_lag
      if (!any(ok)) next
      jj <- jj[ok]; b <- lag[ok]
      dr2 <- rowSums((pos[jj, , drop = FALSE] -
                      matrix(pos[i, ], length(jj), ndim, byrow = TRUE))^2)
      s[b]   <- s[b] + dr2
      cnt[b] <- cnt[b] + 1L
    }
    m   <- ifelse(cnt > 0L, s / cnt, NA_real_)
    lags <- seq_len(max_lag)
    keep <- which(is.finite(m) & lags >= min_lag & lags <= max_lag)
    if (length(keep) < 5L) return(NA_real_)
    tt <- lags[keep] * t_c; yy <- m[keep]
    ne <- max(2L, floor(0.15 * length(keep)))
    tf <- stats::median(tt[seq_len(ne)]); tl <- stats::median(tt[seq(length(tt) - ne + 1L, length(tt))])
    mf <- stats::median(yy[seq_len(ne)]); ml <- stats::median(yy[seq(length(yy) - ne + 1L, length(yy))])
    if (!is.finite(mf) || mf <= 0 || !is.finite(tl) || tl <= tf) return(NA_real_)
    (ml / mf) / (tl / tf)
  }

  cells <- expand.grid(n = pmax(5L, as.integer(n_grid)), D = D_grid)
  out <- lapply(seq_len(nrow(cells)), function(g) {
    if (is.function(progress)) progress(g / nrow(cells))
    N <- cells$n[g]; Dg <- cells$D[g]
    v <- vapply(seq_len(n_sim), function(i) ci_one(N, Dg), numeric(1))
    v <- v[is.finite(v)]
    if (length(v) < 10L)
      return(data.frame(n = N, D = Dg, lo05 = NA_real_, med = NA_real_, hi95 = NA_real_))
    q <- stats::quantile(v, c(0.05, 0.50, 0.95), names = FALSE)
    data.frame(n = N, D = Dg, lo05 = q[1], med = q[2], hi95 = q[3])
  })
  do.call(rbind, out)
}


# ── sim_free_geom_msd() ───────────────────────────────────────────────────────
# GOAL: the free-diffusion MSD curve ON a given cell geometry - the reference the
# measured MSD should be read against when the cell is small enough to matter.
#
# WHY IT IS NEEDED. On a closed surface the MSD cannot grow forever: once a molecule
# has explored the shape its displacement saturates, so FREE diffusion already
# produces a bending, plateauing MSD. Fit a straight line across that bend and D comes
# out low; fit the confined model and it reports a "confinement size" that is really
# the cell. The distortion is negligible while the fit window is short and grows with
# it - on a rod (r 400, L 2000 nm, D 200 nm^2/ms) the fitted D runs 0.995 of truth over
# a 50-cycle window, 0.965 at 500 cycles and 0.893 at 2000. This curve makes that
# visible instead of leaving it to be discovered as a fake plateau.
#
# The molecules are FREE throughout: any bend is geometry alone. Comparing the measured
# curve against it separates "the cell is small" from "the molecule is trapped".
#
# lag_cycles : integer cycle lags to evaluate (the MSD panel's own lag grid)
# D_surf     : the TRUE motion D (nm^2/ms) - i.e. the d_eff-corrected D the MSD panel
#              now reports, NOT a raw per-axis value; no projection re-calibration is
#              applied here, unlike sim_free_geom_rnet which is fed a measured per-axis D
# geom/params/mode : as .geom_ops()
# cyc_ms     : simulation time step (the acquisition cycle, or a coarser step)
# proj2d     : TRUE = observe x,y only (2D tracking); FALSE = all three coordinates
# M          : number of simulated molecules
# returns    : data.frame(lag, tau, msd) in cycles / ms / nm^2
sim_free_geom_msd <- function(lag_cycles, D_surf, geom = "rod", params = list(),
                              cyc_ms = 0.5, proj2d = TRUE,
                              mode = c("surface", "volume"), M = 1500L) {
  mode <- match.arg(mode)
  lag_cycles <- sort(unique(as.integer(lag_cycles[is.finite(lag_cycles) & lag_cycles >= 1])))
  if (!length(lag_cycles) || !is.finite(D_surf) || D_surf <= 0) return(NULL)
  ops   <- .geom_ops(geom, params, mode)
  nobs  <- if (proj2d) 2L else 3L
  nmax  <- max(lag_cycles)
  sd    <- sqrt(2 * D_surf * cyc_ms)

  P  <- ops$rand_start(M)
  # Keep only the positions the requested lags actually need, so memory stays flat
  # even for very long windows.
  want <- c(0L, lag_cycles)
  keep <- new.env(parent = emptyenv())
  assign("0", P[, seq_len(nobs), drop = FALSE], envir = keep)
  for (i in seq_len(nmax)) {
    P <- ops$one_step(P, sd)
    if (i %in% lag_cycles) assign(as.character(i), P[, seq_len(nobs), drop = FALSE], envir = keep)
  }
  P0 <- get("0", envir = keep)
  msd <- vapply(lag_cycles, function(L) {
    Pl <- get(as.character(L), envir = keep)
    mean(rowSums((Pl - P0)^2))
  }, numeric(1))
  data.frame(lag = lag_cycles, tau = lag_cycles * cyc_ms, msd = msd)
}
