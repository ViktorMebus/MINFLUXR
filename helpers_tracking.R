## ─────────────────────────────────────────────────────────────────────────────
## helpers_tracking.R  —  tracking-mode-only helper functions
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
##   - jsonlite is NOT library()'d — always call jsonlite:: explicitly, because
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
##   deadcode_tracking.R (2026-07-20) — see that file to restore them.
## ─────────────────────────────────────────────────────────────────────────────


# ── tid_color_map() ──────────────────────────────────────────────────────────
# GOAL: give every trace ID (TID) its own visually distinct, STABLE colour.
#
# Why stability matters: the same TID must keep the same colour across
# re-renders, filter changes and playback frames. Otherwise tracks appear to
# swap identity mid-animation. That is the whole reason the shuffle is seeded.
#
# tids    : vector of trace IDs (any type — coerced to character)
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
  # brewer.pal() errors below 3 colours, and errors above the palette's max —
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


# ── sim_free_geom_rnet() ──────────────────────────────────────────────────────
# GOAL: the GEOMETRY-AWARE free-diffusion null for the displacement test (Step 2), for
# ANY of the simulator geometries. Simulate FREE Brownian diffusion on the chosen shape
# and return each track's NET start→end displacement in the 2D image projection — so
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
# geometry — one short single-step pass — and drive the null at D_surf = D_meas / factor,
# so both nulls agree at short T and the closed-surface band diverges only where finite
# size + curvature bite. Set FALSE to feed D straight through.
#
# dur_ms  : track durations (ms) to match (resample from real data)
# D_nm2ms : measured per-axis diffusion coefficient (nm²/ms)
# geom    : "sphere" | "rod" | "ellipsoid" | "membrane" | "free"
# params  : list(R, r, Lcyl, a, b, c, mtilt) — only the relevant ones are read
# returns : numeric vector (length M) of net displacements R_net (nm), + attr "D_surf"
sim_free_geom_rnet <- function(dur_ms, D_nm2ms, geom = "rod", params = list(),
                               cyc_ms = 0.5, proj2d = TRUE, match_projection = TRUE,
                               mode = c("surface", "volume")) {
  mode <- match.arg(mode)
  `%|NA|%` <- function(a, b) if (!is.null(a) && is.finite(a)) a else b
  M <- length(dur_ms)
  if (M == 0L || !is.finite(D_nm2ms) || D_nm2ms <= 0) return(numeric(0))
  n_cyc <- pmax(2L, as.integer(round(dur_ms / cyc_ms))); nmax <- max(n_cyc)

  R    <- params$R    %|NA|% 500;  rr   <- params$r %|NA|% 400;  Lcyl <- params$Lcyl %|NA|% 2000
  abc  <- if (identical(geom, "sphere")) c(R, R, R)
          else c(params$a %|NA|% 500, params$b %|NA|% 500, params$c %|NA|% 900)   # ellipsoid
  mtilt <- (params$mtilt %|NA|% 0) * pi / 180
  is_rod <- identical(geom, "rod"); is_ell <- geom %in% c("sphere", "ellipsoid")
  is_flat <- geom %in% c("free", "membrane")

  # Tilted-plane basis for the membrane (flat sheet inclined by mtilt); free = untilted 3D.
  n_hat <- c(sin(mtilt), 0, cos(mtilt)); e1 <- c(cos(mtilt), 0, -sin(mtilt)); e2 <- c(0, 1, 0)

  # ── vectorised geometry operators (P is M×3) ──────────────────────────────────
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
    if (is_ell) geo_project(sweep(matrix(rnorm(3 * m), m, 3), 2, abc, "*"))
    else { s <- runif(m, 0, Lcyl); phi <- runif(m, 0, 2 * pi); cbind(s, rr * cos(phi), rr * sin(phi)) }
  }
  # Reflect a point back INSIDE the shape (used for interior/volume diffusion).
  geo_reflect <- function(P) {
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
      d <- matrix(rnorm(3 * m), m, 3); d <- d / sqrt(dotrows(d, d))
      sweep(d * runif(m)^(1 / 3), 2, abc, "*")
    } else {                                                   # rod: rejection sample the capsule
      out <- matrix(NA_real_, m, 3); need <- seq_len(m)
      while (length(need)) {
        k <- length(need)
        p <- cbind(runif(k, -rr, Lcyl + rr), runif(k, -rr, rr), runif(k, -rr, rr))
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
  # One free step for all M points at step SD `sd`, returning the new M×3 positions.
  # surface → tangent step then re-project onto the boundary;
  # volume  → isotropic 3D step then reflect back inside.
  one_step <- function(P, sd) {
    G <- matrix(rnorm(nrow(P) * 3, 0, sd), nrow(P), 3)
    if (is_flat) {
      if (identical(geom, "free")) P + G
      else P + outer(G[, 1], e1) + outer(G[, 2], e2)          # membrane: 2 DOF in the tilted plane
    } else if (is_vol) geo_reflect(P + G)
    else { N <- geo_normal(P); geo_project(P + (G - N * dotrows(G, N))) }  # tangent step + reproject
  }
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
    P0 <- if (is_flat) matrix(0, ncal, 3) else if (is_vol) geo_rand_in(ncal) else geo_rand(ncal)
    d2 <- to_obs(one_step(P0, sdc)) - to_obs(P0)
    f2 <- mean(rowSums(d2^2)) / (nobs * sdc^2)
    if (is.finite(f2) && f2 > 0.05) D_surf <- D_nm2ms / f2
  }
  sd <- sqrt(2 * D_surf * cyc_ms)

  # ── main pass ─────────────────────────────────────────────────────────────────
  P  <- if (is_flat) matrix(0, M, 3) else if (is_vol) geo_rand_in(M) else geo_rand(M)
  S2 <- to_obs(P); E2 <- S2
  for (i in 2:nmax) {
    alive <- i <= n_cyc
    Pn <- one_step(P, sd); P[alive, ] <- Pn[alive, ]
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
# tt      : sorted numeric vector of times (ASSUMED sorted — not checked)
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
# any other quantitative estimator — dropping points changes the step statistics.
#
# df         : data.frame with tid, x, y, [z], and the time column `tcol`
# max_pts    : maximum points to keep per track
# tcol       : name of the time column to sort/keep (default "tim0")
# extra_cols : extra columns to carry through untouched (e.g. per-point colour)
# returns    : data.frame(tid, x, y, t, [z], <extra_cols>) — note the time
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


