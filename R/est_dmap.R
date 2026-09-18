## ─────────────────────────────────────────────────────────────────────────────
## est_dmap.R  -  spatially resolved diffusion map (D per spatial cell)
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The estimator behind the Tracking tab's D-map panel, as a pure function.
##   Moved out of dmap_bundle() on 2026-09-17 (publication audit C1).
##
## METHOD
##   Every step is assigned to the grid cell containing its midpoint. Per cell,
##   with first moments only (a median of single-step r² is ln2 × the mean in 2D):
##       D = max(0, mean(r²/(4Δt)) − σ²·mean(1/Δt)) / (1 − 2·R_abs·mean(1/Δt))
##   i.e. the 2D Brownian step variance 4D(Δt − 2R_abs) + 4σ² inverted for D.
##   SE = SD of the per-step contributions / √n, scaled by the same blur factor.
##
## BLUR TERM (added 2026-09-17). The 2026-09-16 fix was verified on simulated
## tracks without motion blur. With blur (R = 1/6.2, GP67H timings) the unblurred
## form reads D about 30% low, the same bias found and fixed in the JDD.
##
## UNITS: nm, ms, D nm²/ms.
## ─────────────────────────────────────────────────────────────────────────────


# ── dmap_compute() ────────────────────────────────────────────────────────────
# df              : data.frame(tid, x, y, tim), tim in ms
# n_bins          : cells per axis (≥ 5)
# min_pts         : minimum steps per cell (≥ 3)
# t_c             : cycle time (ms), for integrated_only
# integrated_only : keep only single-cycle steps (their midpoints are meaningful)
# sigma_loc       : σ_loc (nm) to subtract; NA or ≤ 0 = no subtraction
# R_abs           : blur offset in ms (MSD_BLUR_R · t_c); 0 = no blur correction
#
# RETURNS list(D_mat, SE_mat, N_mat, x_centers, y_centers, cell_w_nm, cell_h_nm,
#   n_bins, n_steps, dt_med_ms, integrated_only, noise_subtracted, sigma_loc_nm),
#   or list(error = "<message>") when there is too little data.
dmap_compute <- function(df, n_bins = 20L, min_pts = 10L, t_c = NA_real_,
                         integrated_only = FALSE, sigma_loc = NA_real_, R_abs = 0) {
  n_bins  <- max(5L, as.integer(n_bins))
  min_pts <- max(3L, as.integer(min_pts))

  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
  df <- df[order(df$tid, df$tim), ]
  sp <- split(df, df$tid)

  # Build step table: midpoint position, dt, r²
  step_rows <- lapply(sp, function(tb) {
    if (nrow(tb) < 2L) return(NULL)
    n  <- nrow(tb)
    i1 <- seq_len(n - 1L); i2 <- i1 + 1L
    dt <- tb$tim[i2] - tb$tim[i1]
    data.frame(
      mx = 0.5 * (tb$x[i1] + tb$x[i2]),
      my = 0.5 * (tb$y[i1] + tb$y[i2]),
      dt = dt,
      r2 = (tb$x[i2] - tb$x[i1])^2 + (tb$y[i2] - tb$y[i1])^2,
      stringsAsFactors = FALSE)
  })
  steps <- do.call(rbind, step_rows[!sapply(step_rows, is.null)])
  if (is.null(steps)) return(list(error = "Too few steps - relax filters or lower 'Min steps per cell'."))
  steps <- steps[is.finite(steps$dt) & steps$dt > 0 & is.finite(steps$r2), ]
  if (nrow(steps) < min_pts * 2)
    return(list(error = "Too few steps - relax filters or lower 'Min steps per cell'."))

  # Integrated (single-cycle) steps only: a step spanning dark cycles has a large Δt
  # and a spatially meaningless midpoint, so mixing them in smears the map.
  integrated_only <- isTRUE(integrated_only) && is.finite(t_c) && t_c > 0
  if (integrated_only) {
    keep_step <- pmax(1L, as.integer(round(steps$dt / t_c))) == 1L
    steps <- steps[keep_step, , drop = FALSE]
    if (nrow(steps) < min_pts * 2)
      return(list(error = "Too few integrated (single-cycle) steps - untick 'Integrated steps only' or relax filters."))
  }

  dt_med <- median(steps$dt, na.rm = TRUE)
  if (!is.finite(dt_med) || dt_med <= 0) return(list(error = "Cannot determine median Δt."))

  # Per-step apparent 2D diffusion D_raw = r²/(4Δt); the static noise floor σ²/Δt
  # is subtracted PER CELL (never per step, which would clamp and rectify).
  steps$D_raw <- steps$r2 / (4 * steps$dt)
  noise_subtracted <- is.finite(sigma_loc) && sigma_loc > 0
  sig2 <- if (noise_subtracted) sigma_loc^2 else 0
  R_abs <- if (is.finite(R_abs) && R_abs > 0) R_abs else 0

  x_rng    <- range(steps$mx); y_rng <- range(steps$my)
  x_breaks <- seq(x_rng[1], x_rng[2], length.out = n_bins + 1L)
  y_breaks <- seq(y_rng[1], y_rng[2], length.out = n_bins + 1L)
  cell_w   <- diff(x_rng) / n_bins
  cell_h   <- diff(y_rng) / n_bins

  steps$ix <- pmin(findInterval(steps$mx, x_breaks), n_bins)
  steps$iy <- pmin(findInterval(steps$my, y_breaks), n_bins)

  D_mat  <- matrix(NA_real_, n_bins, n_bins)
  SE_mat <- matrix(NA_real_, n_bins, n_bins)
  N_mat  <- matrix(0L,       n_bins, n_bins)

  cell_idx <- steps$ix + (steps$iy - 1L) * n_bins
  by_cell  <- split(seq_len(nrow(steps)), cell_idx)

  for (nm in names(by_cell)) {
    rows <- by_cell[[nm]]
    rows <- rows[is.finite(steps$D_raw[rows])]
    if (length(rows) < min_pts) next
    idx <- as.integer(nm)
    ix  <- ((idx - 1L) %% n_bins) + 1L
    iy  <- ((idx - 1L) %/% n_bins) + 1L
    # per-step contribution to the estimate, so its SD gives the SE directly
    contrib  <- steps$D_raw[rows] - sig2 / steps$dt[rows]
    # blur shortens each step's effective diffusion time to Δt − 2·R_abs
    blur_fac <- if (R_abs > 0) max(0.05, 1 - 2 * R_abs * mean(1 / steps$dt[rows])) else 1
    D_mat[ix, iy]  <- max(0, mean(contrib)) / blur_fac
    SE_mat[ix, iy] <- stats::sd(contrib) / sqrt(length(rows)) / blur_fac
    N_mat[ix, iy]  <- length(rows)
  }

  x_centers <- 0.5 * (x_breaks[-1] + x_breaks[-length(x_breaks)])
  y_centers <- 0.5 * (y_breaks[-1] + y_breaks[-length(y_breaks)])

  list(D_mat = D_mat, SE_mat = SE_mat, N_mat = N_mat,
       x_centers = x_centers, y_centers = y_centers,
       cell_w_nm = cell_w, cell_h_nm = cell_h,
       n_bins = n_bins, n_steps = nrow(steps),
       dt_med_ms = dt_med,
       integrated_only = integrated_only,
       noise_subtracted = noise_subtracted,
       sigma_loc_nm = if (noise_subtracted) sigma_loc else NA_real_,
       R_abs = R_abs)
}
