## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_shapes.R  -  Shape fitting
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Fits of sphere, ellipsoid, helix, ring and arc to a cluster (geometric
##   start, then maximum likelihood), their projections, surface samples and
##   meshes, residual distances, and the alpha-shape surface.
##
## CONTENTS (in file order)
##   calculate_distance_to_fit, build_center_args, robust_center,
##   make_alpha_mesh3d, mle_fit_sphere, mle_project_sphere, mle_sample_sphere,
##   .eul2R, mle_fit_ellipsoid, mle_project_ellipsoid, mle_sample_ellipsoid,
##   .ang_to_u, .rot_to_u, .helix_point, .helix_newton, .helix_d2,
##   .project_t_helix, .grid_search_orient, mle_fit_helix, mle_sample_helix,
##   .ring_basis, mle_project_ring, mle_fit_ring, mle_sample_ring,
##   mle_project_arc, mle_fit_arc, mle_sample_arc, .tube_mesh_from_frames,
##   make_ring_tube_mesh, make_arc_tube_mesh, make_helix_tube_mesh,
##   add_distance_lines_to_plot, compute_convex_hull_centroid,
##   medianize_coords, spread_from_group, fit_sphere_geometry,
##   fit_ellipsoid_geometry, fit_helix_geometry, fit_ring_geometry,
##   fit_arc_geometry, generate_ring_from_cluster2,
##   generate_arc_from_cluster2, alpha_surface_3d,
##   generate_sphere_from_cluster2, generate_ellipsoid_from_cluster2,
##   generate_helix_from_cluster2
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

calculate_distance_to_fit <- function(df           = NULL,
                                      coords_matrix = NULL,
                                      fit_type     = c("sphere", "ellipsoid", "helix", "ring", "arc", "best"),
                                      cluster2_id  = NULL,
                                      use_cluster1 = FALSE,
                                      turns        = 3,
                                      helix_n      = 5000,
                                      arc_span_hint = 0) {

  fit_type <- match.arg(fit_type)
  turns    <- max(1, as.numeric(turns))

  ## ── which models to compute (lazy: skip models not required) ─────────────
  need_s <- fit_type %in% c("sphere",    "best")
  need_e <- fit_type %in% c("ellipsoid", "best")
  need_h <- fit_type %in% c("helix",     "best")
  need_r <- fit_type %in% c("ring",      "best")
  need_a <- fit_type %in% c("arc",       "best")
  .empty_dd <- list(rmsd = NA_real_, meandist = NA_real_, sd = NA_real_,
                    dist = numeric(0), nearest = NULL)

  ## ── resolve coordinate matrix (skip re-medianing if already provided) ----
  if (is.null(coords_matrix)) {
    if (is.null(df)) return(NULL)
    if (!is.null(cluster2_id)) df <- dplyr::filter(df, cluster2 == cluster2_id)
    if (nrow(df) == 0) return(NULL)
    medians <- if (use_cluster1 && "cluster" %in% names(df)) {
      df %>% dplyr::group_by(cluster) %>%
        dplyr::summarise(dplyr::across(c(x, y, z), median), .groups = "drop")
    } else {
      df %>% dplyr::group_by(tid) %>%
        dplyr::summarise(dplyr::across(c(x, y, z), median), .groups = "drop")
    }
    if (nrow(medians) == 0) return(NULL)
    coords_matrix <- as.matrix(medians[, c("x", "y", "z")])
  } else {
    medians <- as.data.frame(coords_matrix)
    colnames(medians) <- c("x", "y", "z")
  }

  med_pts <- coords_matrix

  ## ── inner helper: distances → RMSD + mean distance + SD ------------------
  rmsd_sd <- function(d) {
    d <- as.numeric(d); d <- d[is.finite(d) & d >= 0]
    list(
      rmsd     = if (length(d)) sqrt(mean(d^2)) else NA_real_,
      meandist = if (length(d)) mean(d)          else NA_real_,
      sd       = if (length(d) > 1) stats::sd(d) else NA_real_
    )
  }

  ## ── Sphere ───────────────────────────────────────────────────────────────
  dd_sphere <- if (need_s) {
    sph_geo <- fit_sphere_geometry(med_pts)
    pr_s    <- mle_project_sphere(med_pts, sph_geo$center, sph_geo$radius)
    d_s     <- as.numeric(pr_s$dist); rs_s <- rmsd_sd(d_s)
    list(rmsd=rs_s$rmsd, meandist=rs_s$meandist, sd=rs_s$sd, dist=d_s, nearest=pr_s$nearest)
  } else .empty_dd

  ## ── Ellipsoid ─────────────────────────────────────────────────────────────
  dd_ellipsoid <- if (need_e) {
    ell_geo <- fit_ellipsoid_geometry(med_pts)
    if (!is.null(ell_geo)) {
      pr_e <- mle_project_ellipsoid(med_pts, ell_geo$center, ell_geo$axes, ell_geo$R)
      d_e  <- as.numeric(pr_e$dist); rs_e <- rmsd_sd(d_e)
      list(rmsd=rs_e$rmsd, meandist=rs_e$meandist, sd=rs_e$sd, dist=d_e, nearest=pr_e$nearest)
    } else .empty_dd
  } else .empty_dd

  ## ── Helix ─────────────────────────────────────────────────────────────────
  dd_helix <- if (need_h) {
    hel_geo <- fit_helix_geometry(med_pts, turns = turns, n_points = helix_n)
    t_h <- .project_t_helix(med_pts, hel_geo$proj_center, hel_geo$axis,
                             hel_geo$radius, hel_geo$pitch, hel_geo$phi0,
                             check_neighbors = TRUE)
    H_h <- t(vapply(t_h, .helix_point, numeric(3),
                    c=hel_geo$proj_center, u=hel_geo$axis,
                    r=hel_geo$radius, p=hel_geo$pitch, phi0=hel_geo$phi0))
    d_h <- sqrt(rowSums((med_pts - H_h)^2)); rs_h <- rmsd_sd(d_h)
    list(rmsd=rs_h$rmsd, meandist=rs_h$meandist, sd=rs_h$sd, dist=d_h, nearest=H_h)
  } else .empty_dd

  ## ── Ring ──────────────────────────────────────────────────────────────────
  dd_ring <- if (need_r) {
    ring_geo <- fit_ring_geometry(med_pts)
    pr_r <- mle_project_ring(med_pts, ring_geo$center, ring_geo$normal, ring_geo$radius)
    d_r  <- as.numeric(pr_r$dist); rs_r <- rmsd_sd(d_r)
    list(rmsd=rs_r$rmsd, meandist=rs_r$meandist, sd=rs_r$sd, dist=d_r, nearest=pr_r$nearest)
  } else .empty_dd

  ## ── Arc ───────────────────────────────────────────────────────────────────
  dd_arc <- if (need_a) {
    arc_geo <- fit_arc_geometry(med_pts, arc_span_hint = arc_span_hint)
    pr_a <- mle_project_arc(med_pts, arc_geo$center, arc_geo$normal, arc_geo$radius,
                             arc_geo$t_center, arc_geo$arc_half)
    d_a  <- as.numeric(pr_a$dist); rs_a <- rmsd_sd(d_a)
    list(rmsd=rs_a$rmsd, meandist=rs_a$meandist, sd=rs_a$sd, dist=d_a, nearest=pr_a$nearest)
  } else .empty_dd

  ## ── choose which fit's distances to embed in the return dataframe ---------
  if (fit_type == "best") {
    rmsd_vals <- c(sphere = dd_sphere$rmsd, ellipsoid = dd_ellipsoid$rmsd,
                   helix  = dd_helix$rmsd,  ring = dd_ring$rmsd, arc = dd_arc$rmsd)
    rmsd_vals <- rmsd_vals[is.finite(rmsd_vals)]
    if (!length(rmsd_vals)) return(NULL)
    best_name <- names(which.min(rmsd_vals))
  } else {
    best_name <- fit_type
  }

  pick <- switch(best_name,
                 sphere    = dd_sphere,
                 ellipsoid = dd_ellipsoid,
                 helix     = dd_helix,
                 ring      = dd_ring,
                 arc       = dd_arc)
  if (is.null(pick) || !length(pick$dist)) return(NULL)

  ## ── build return dataframe (same columns as before for line-drawing compat)
  n <- nrow(med_pts)
  medians$distance_to_fit <- pick$dist[seq_len(n)]
  if (!is.null(pick$nearest) && nrow(pick$nearest) == n) {
    medians$nearest_x <- pick$nearest[, 1]
    medians$nearest_y <- pick$nearest[, 2]
    medians$nearest_z <- pick$nearest[, 3]
  } else {
    medians$nearest_x <- NA_real_
    medians$nearest_y <- NA_real_
    medians$nearest_z <- NA_real_
  }

  attr(medians, "rmsd_sphere")        <- dd_sphere$rmsd
  attr(medians, "meandist_sphere")    <- dd_sphere$meandist
  attr(medians, "sd_sphere")          <- dd_sphere$sd
  attr(medians, "rmsd_ellipsoid")     <- dd_ellipsoid$rmsd
  attr(medians, "meandist_ellipsoid") <- dd_ellipsoid$meandist
  attr(medians, "sd_ellipsoid")       <- dd_ellipsoid$sd
  attr(medians, "rmsd_helix")         <- dd_helix$rmsd
  attr(medians, "meandist_helix")     <- dd_helix$meandist
  attr(medians, "sd_helix")           <- dd_helix$sd
  attr(medians, "rmsd_ring")          <- dd_ring$rmsd
  attr(medians, "meandist_ring")      <- dd_ring$meandist
  attr(medians, "sd_ring")            <- dd_ring$sd
  attr(medians, "rmsd_arc")           <- dd_arc$rmsd
  attr(medians, "meandist_arc")       <- dd_arc$meandist
  attr(medians, "sd_arc")             <- dd_arc$sd
  attr(medians, "best_fit")           <- best_name

  return(medians)
}

# GOAL: the centre arguments for an MLE shape fit.
#   mode   : "fixed" (centre held at the robust centre), "reg" (free, pulled toward
#            it with strength `lambda`) or "free" (free, started there)
#   lambda : regularisation strength, used by "reg" only
build_center_args <- function(X, mode, lambda) {
  c_robust <- robust_center(X)
  switch(mode,
         "fixed" = list(center_fixed = c_robust),
         "reg"   = list(center_init = c_robust,
                        center_prior = c_robust,
                        lambda = lambda),
         "free"  = list(center_init = c_robust)
  )
}

# ---------- MLE helpers (orthogonal-distance) ----------
# (A) Sphere
# --- Robust center guess (geometric median -> mean) ---
# Starting centre for the MLE fits, via build_center_args(). Was the convex-hull
# centroid, which is set by the extreme points (breakdown 0); geometric_median()
# has breakdown 0.5, so a loosely-attached member no longer drags the start.
# The MLE optimises the centre from here, so this changes where the search
# begins, not what it converges to - except where the old start was pulled far
# enough out to find a different local optimum.
robust_center <- function(X) {
  tryCatch(geometric_median(X), error = function(e) colMeans(X))
}

# =============================================================================
# SHAPE FITTING - fit a geometric primitive to a cluster of localizations
# =============================================================================
#
# GOAL OF THIS WHOLE BLOCK
#   Given the points belonging to one cluster, ask "what SHAPE is this?" and fit
#   the best sphere / ellipsoid / helix / ring / arc to it. That turns a cloud of
#   localizations into interpretable numbers: a radius, a pitch, an axis length.
#
# THE FITTING PRINCIPLE - ORTHOGONAL-DISTANCE (total least squares)
#   For a shape S with parameters θ, minimise the sum of squared PERPENDICULAR
#   distances from each point to the SURFACE:
#
#        θ* = argmin_θ  Σ_i  || x_i − P_S(x_i; θ) ||²
#
#   where P_S(x) is the nearest point on the shape. This is the right criterion
#   here because localization error is isotropic in space - there is no
#   privileged "vertical" axis to measure residuals along, as ordinary least
#   squares would assume. It is the same reason these are called MLE fits: for
#   isotropic Gaussian error, minimising orthogonal distance IS the maximum
#   likelihood solution.
#
#   Because P_S has no closed form for the curved shapes, several fits iterate:
#   project -> re-fit parameters -> re-project, until the residual stops falling.
#
# THE NAMING CONVENTION (this is the key to navigating ~30 functions)
#   For each shape S in {sphere, ellipsoid, helix, ring, arc} there is a family:
#
#     fit_S_geometry(X) ...... convenience wrapper: fit and return the geometry
#                              in one call (what most callers actually use)
#     mle_fit_S(X, w) ........ the real optimiser. Returns the fitted parameters
#                              plus the residual/RMSD used for model ranking.
#     mle_project_S(X, ...) .. nearest point ON the shape for each input point.
#                              This is what residuals and the "distance to fit"
#                              lines are computed from.
#     mle_sample_S(...) ...... sample points along/over the fitted shape, used
#                              to DRAW it (not to fit it).
#     make_S_mesh() /
#       make_S_tube_mesh() ... turn the sampled shape into plotly 3D geometry.
#                              "tube" variants sweep a small circle along a
#                              curve, because a 1D curve is invisible in 3D.
#     generate_S_from_cluster2(...)  per-cluster wrapper used by the server.
#
#   So: to change how a shape is FITTED edit mle_fit_S; to change how it LOOKS
#   edit mle_sample_S / make_S_*_mesh. They are deliberately independent.
#
# MODEL SELECTION - which shape wins?
#   RMSD always favours the most flexible shape (more parameters always fit
#   better), so ranking uses AICc/BIC, which penalise parameter count. AICc is
#   preferred for the small point counts typical of one cluster, and is left
#   blank when a cluster has too few points to justify that model (n <= params+2).
#
# HONEST CAVEAT (mirrored by the method box in the UI)
#   The shape dictionary is IMPOSED by this tool - a cluster is scored only
#   against the five shapes offered, so "best" means best OF THESE, never proof
#   of biological geometry. With few, blinking-inflated localizations several
#   shapes are often statistically indistinguishable. Treat the ranking as a
#   hypothesis, and check the residuals rather than the label.
# =============================================================================

make_alpha_mesh3d <- function(pts_xyz, alpha) {
  sh <- alphashape3d::ashape3d(pts_xyz, alpha = alpha)
  tri <- sh$triang
  if (is.null(tri) || !nrow(tri)) return(NULL)
  
  v <- sh$x  # vertices (n x 3)
  
  # plotly mesh3d uses 0-based indices
  list(
    x = v[, 1], y = v[, 2], z = v[, 3],
    i = tri[, 1] - 1L,
    j = tri[, 2] - 1L,
    k = tri[, 3] - 1L
  )
}

# ========= Sphere (free / fixed / regularized) =========
mle_fit_sphere <- function(X, w = NULL,
                           center_init  = NULL,
                           center_fixed = NULL,
                           center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 4) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)
  
  if (!is.null(center_fixed)) {
    c <- as.numeric(center_fixed)
    L <- sqrt(rowSums((X - matrix(c, n, 3, TRUE))^2))
    r <- sum(w * L) / sum(w)
    return(list(center = c, radius = as.numeric(abs(r))))
  }
  
  c0 <- if (!is.null(center_init)) as.numeric(center_init) else colMeans(X)
  d0 <- sqrt(rowSums((X - matrix(c0, n, 3, TRUE))^2))
  r0 <- median(d0)
  # --- NEW: robust upper bound for radius ---
  r_max <- max(stats::quantile(d0, 0.97, na.rm = TRUE), r0 * 2)
  
  par0 <- c(c0, r0)
  # Bound center to an expanded data bounding box and r to [0, r_max].
  # Without centre bounds, the joint (c, r) objective is non-convex: the
  # optimizer can drift c outside the data cloud and inflate r to compensate
  # while keeping all residuals small - a geometrically meaningless solution.
  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)  # ≥10 nm margin
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand
  obj <- function(par) {
    c <- par[1:3]; r <- abs(par[4])
    L <- sqrt(rowSums((X - matrix(c, n, 3, TRUE))^2))
    # Divide by n so the objective is mean squared residual (nm²), independent of
    # cluster size.  This makes lambda directly interpretable across clusters:
    # lambda ≈ 1 means a 1 nm centre displacement costs as much as 1 nm² of mean residual.
    val <- sum(w * (L - r)^2) / n
    if (!is.null(center_prior) && lambda > 0) val <- val + lambda * sum((c - center_prior)^2)
    val
  }
  fit <- optim(
    par0, obj, method = "L-BFGS-B", control = list(maxit = 500),
    lower = c(c_lo, 0),       # center bounded to expanded data box; r >= 0
    upper = c(c_hi, r_max)    # r <= r_max (hard bound replaces soft penalty)
  )
  list(center = fit$par[1:3], radius = abs(fit$par[4]))
}

mle_project_sphere <- function(X, center, radius) {
  V <- sweep(X, 2, center); L <- sqrt(rowSums(V^2))
  scale <- ifelse(L > 0, radius / L, 0)
  Q <- sweep(V, 1, scale, "*")
  nearest <- sweep(Q, 2, center, "+")
  dist <- abs(L - radius)
  list(nearest = nearest, dist = dist)
}

mle_sample_sphere <- function(center, radius, n = 4000) {
  phi <- runif(n, 0, 2*pi)
  u   <- runif(n, -1, 1); theta <- acos(u)
  x <- radius * sin(theta) * cos(phi)
  y <- radius * sin(theta) * sin(phi)
  z <- radius * cos(theta)
  cbind(x, y, z) + matrix(center, n, 3, TRUE)
}

# (B) Ellipsoid  (rotation via yaw/pitch/roll)
.eul2R <- function(yaw, pitch, roll) {
  Ry <- matrix(c(cos(yaw),-sin(yaw),0, sin(yaw),cos(yaw),0, 0,0,1),3,3,TRUE)
  Rp <- matrix(c(cos(pitch),0,sin(pitch), 0,1,0, -sin(pitch),0,cos(pitch)),3,3,TRUE)
  Rr <- matrix(c(1,0,0, 0,cos(roll),-sin(roll), 0,sin(roll),cos(roll)),3,3,TRUE)
  Ry %*% Rp %*% Rr
}

# ========= Ellipsoid (free / fixed / regularized center) =========
# Objective: sum of squared Euclidean distances from each localization to the
# nearest point on the ellipsoid surface - identical units (nm²) to the sphere
# fitter and to the reported RMSD, making fits directly comparable.
#
# Previous version minimized the dimensionless radial ratio sum((s-1)^2) where
# s = ||Xp/axes||.  That objective has a scale degeneracy: the optimizer could
# translate the centre and simultaneously scale the axes while keeping s ≈ 1
# for all points, producing centres far outside the data cloud.  The orthogonal-
# distance objective eliminates this degeneracy because any displacement of the
# centre away from the data increases the distances unambiguously.
#
# Axis bounds (log-space) prevent runaway growth; centre bounds (expanded data
# bounding box) catch any remaining drift for the free-centre case.
mle_fit_ellipsoid <- function(X, w = NULL,
                              center_init  = NULL,
                              center_fixed = NULL,
                              center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 6) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)

  # Shared initialisation: PCA at a given centre + axis/rotation starting values
  .ell_init <- function(cc) {
    Xc <- sweep(X, 2, cc)
    pc <- prcomp(Xc, center = FALSE); R0 <- pc$rotation
    Q  <- Xc %*% R0
    a0 <- pmax(diff(stats::quantile(Q[,1], c(0.05, 0.95))) / 2, 1e-6)
    b0 <- pmax(diff(stats::quantile(Q[,2], c(0.05, 0.95))) / 2, 1e-6)
    c0 <- pmax(diff(stats::quantile(Q[,3], c(0.05, 0.95))) / 2, 1e-6)
    yaw0   <- atan2(R0[2,1], R0[1,1])
    pitch0 <- asin(max(-1, min(1, -R0[3,1])))  # clamp avoids NaN at ±1
    roll0  <- atan2(R0[3,2], R0[3,3])
    # Log-axis bounds: axes must stay in [1% of smallest initial half-extent, 5× RMS spread]
    d_rms  <- sqrt(mean(rowSums(Xc^2)))
    ax_lo  <- log(pmax(min(a0, b0, c0) * 0.01, 1e-3))
    ax_hi  <- log(max(d_rms * 5, 1))
    list(a0 = a0, b0 = b0, c0 = c0, R0 = R0,
         yaw0 = yaw0, pitch0 = pitch0, roll0 = roll0,
         ax_lo = ax_lo, ax_hi = ax_hi)
  }

  if (!is.null(center_fixed)) {
    cc  <- as.numeric(center_fixed)
    ini <- .ell_init(cc)
    par0 <- c(log(c(ini$a0, ini$b0, ini$c0)), ini$yaw0, ini$pitch0, ini$roll0)
    obj <- function(par) {
      axes_now <- exp(par[1:3])
      R  <- .eul2R(par[4], par[5], par[6])
      pr <- mle_project_ellipsoid(X, cc, axes_now, R)
      sum(w * pr$dist^2) / n   # mean squared distance - scale-independent of n
    }
    fit <- optim(par0, obj, method = "L-BFGS-B",
                 lower = c(rep(ini$ax_lo, 3), -Inf, -Inf, -Inf),
                 upper = c(rep(ini$ax_hi, 3),  Inf,  Inf,  Inf),
                 control = list(maxit = 500))
    p <- fit$par
    return(list(center = cc, axes = exp(p[1:3]), R = .eul2R(p[4], p[5], p[6])))
  }

  # Free centre (optionally regularized toward center_prior)
  cc0 <- if (!is.null(center_init)) as.numeric(center_init) else colMeans(X)
  ini <- .ell_init(cc0)
  par0 <- c(cc0, log(c(ini$a0, ini$b0, ini$c0)), ini$yaw0, ini$pitch0, ini$roll0)

  obj <- function(par) {
    cc       <- par[1:3]
    axes_now <- exp(par[4:6])
    R        <- .eul2R(par[7], par[8], par[9])
    pr  <- mle_project_ellipsoid(X, cc, axes_now, R)
    val <- sum(w * pr$dist^2) / n   # mean squared distance - scale-independent of n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((cc - center_prior)^2)
    val
  }

  # Centre bounded to expanded data bounding box to prevent drift
  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand

  fit <- optim(par0, obj, method = "L-BFGS-B",
               lower = c(c_lo, rep(ini$ax_lo, 3), -Inf, -Inf, -Inf),
               upper = c(c_hi, rep(ini$ax_hi, 3),  Inf,  Inf,  Inf),
               control = list(maxit = 500))
  p <- fit$par
  list(center = p[1:3], axes = exp(p[4:6]), R = .eul2R(p[7], p[8], p[9]))
}

mle_project_ellipsoid <- function(X, center, axes, R) {
  # Exact orthogonal nearest-point via vectorised Newton on the KKT condition.
  # Nearest point p to query q (in PCA frame) satisfies:
  #   p_j = q_j * a_j^2 / (a_j^2 + lambda)
  # where lambda solves  F(lambda) = sum_j (q_j*a_j / (a_j^2+lambda))^2 - 1 = 0
  # F'(lambda) = -2 * sum_j q_j^2 * a_j^2 / (a_j^2+lambda)^3
  n    <- nrow(X)
  Xp   <- (X - matrix(center, n, 3, TRUE)) %*% R   # in PCA frame  (n x 3)
  a2   <- axes^2                                     # [a^2, b^2, c^2]
  Xp2  <- Xp^2                                       # precompute  (n x 3)
  r2   <- rowSums(Xp2)                               # squared distance from center

  min_a2 <- min(a2)
  # Detect interior points: sum_j (q_j/a_j)^2 < 1
  inside_ratio <- rowSums(Xp2 / matrix(a2, n, 3, TRUE))
  lam_outside  <- pmax(max(axes), 1) * sqrt(pmax(r2, 1e-20)) - max(a2)
  # Interior points: lambda root lies in (-min_a2, 0); start at midpoint to avoid clamp divergence
  is_interior <- inside_ratio < 1
  lam <- ifelse(is_interior, -min_a2 * 0.5, lam_outside)
  lam <- pmax(lam, -min_a2 + 1e-6)

  for (iter in 1:25) {
    denom <- outer(lam, a2, "+")                      # n x 3 : a_j^2 + lambda_i
    denom <- pmax(denom, 1e-12)
    A2X2  <- matrix(a2, n, 3, TRUE) * Xp2            # a_j^2 * q_j^2
    fval  <- rowSums(A2X2 / denom^2) - 1             # F(lambda)
    fpval <- -2 * rowSums(A2X2 / denom^3)            # F'(lambda)
    step  <- ifelse(abs(fpval) > 1e-30, fval / fpval, 0)
    lam   <- pmax(lam - step, -min_a2 + 1e-6)        # keep denominators positive
    lam[is_interior] <- pmin(lam[is_interior], -1e-12) # interior root < 0; prevent overshoot above 0
    if (max(abs(step)) < 1e-8) break
  }

  denom_f       <- pmax(outer(lam, a2, "+"), 1e-12)  # n x 3
  nearest_local <- matrix(a2, n, 3, TRUE) * Xp / denom_f  # p_j = a_j^2*q_j/(a_j^2+lambda)
  nearest       <- nearest_local %*% t(R) + matrix(center, n, 3, TRUE)
  dist          <- sqrt(rowSums((X - nearest)^2))

  # Guard: verify the projected point actually lies on the ellipsoid surface.
  # A correct KKT solution satisfies sum_j (p_j/a_j)^2 == 1 exactly.
  # If Newton diverged (lambda stuck near -min_a2), nearest_local is astronomically
  # large → surface_check >> 1 → fall back to radial projection.
  surface_check <- rowSums((nearest_local / matrix(axes, n, 3, TRUE))^2)
  runaway <- !is.finite(surface_check) | abs(surface_check - 1) > 0.1 | !is.finite(dist)

  if (any(runaway)) {
    # Radial fallback: project Xp along its own direction to the ellipsoid surface.
    # Point p = Xp * t lands on surface when sum_j (t*q_j/a_j)^2 = 1
    # → t = 1/sqrt(inside_ratio).  Always finite; always on-surface.
    ir_run  <- pmax(inside_ratio[runaway], 1e-20)
    t_rad   <- 1 / sqrt(ir_run)
    Xp_run  <- Xp[runaway, , drop = FALSE]
    xp_norm <- sqrt(rowSums(Xp_run^2))
    nn_local_rad <- Xp_run * t_rad

    # For points at the center (Xp ≈ 0), radial projection is undefined;
    # place the nearest point on the tip of the shortest semi-axis.
    near_ctr <- xp_norm < 1e-12
    if (any(near_ctr)) {
      j_min <- which.min(axes)
      p_short <- numeric(3); p_short[j_min] <- axes[j_min]
      nn_local_rad[near_ctr, ] <- matrix(p_short, sum(near_ctr), 3, TRUE)
    }

    nearest[runaway, ] <- nn_local_rad %*% t(R) + matrix(center, sum(runaway), 3, TRUE)
    dist[runaway]      <- sqrt(rowSums((X[runaway, , drop = FALSE] -
                                        nearest[runaway, , drop = FALSE])^2))
  }

  # Edge case: point at ellipsoid center -> nearest on shortest axis
  at_ctr <- r2 < 1e-20
  if (any(at_ctr)) {
    min_ax            <- which.min(axes)
    p_deg             <- numeric(3); p_deg[min_ax] <- axes[min_ax]
    nearest[at_ctr, ] <- matrix(p_deg %*% t(R) + center, sum(at_ctr), 3, TRUE)
    dist[at_ctr]      <- axes[min_ax]
  }
  list(nearest = nearest, dist = dist)
}

mle_sample_ellipsoid <- function(center, axes, R, n = 4000) {
  phi <- runif(n, 0, 2*pi)
  u   <- runif(n, -1, 1); theta <- acos(u)
  xl <- axes[1]*sin(theta)*cos(phi)
  yl <- axes[2]*sin(theta)*sin(phi)
  zl <- axes[3]*cos(theta)
  (cbind(xl,yl,zl) %*% t(R)) + matrix(center, n, 3, TRUE)
}

# (C) Helix
.ang_to_u <- function(th, ph) c(sin(th)*cos(ph), sin(th)*sin(ph), cos(th))

.rot_to_u <- function(u) {
  ez <- c(0,0,1); v <- pracma::cross(ez, u); s <- sqrt(sum(v*v)); c <- sum(ez*u)
  if (s < 1e-12) return(diag(3))
  Vx <- matrix(c( 0,-v[3], v[2], v[3], 0,-v[1], -v[2], v[1], 0), 3,3, TRUE)
  diag(3) + Vx + Vx%*%Vx*((1-c)/(s^2))
}

.helix_point <- function(t, c, u, r, p, phi0) {
  R <- .rot_to_u(u)
  H <- cbind(r*cos(t+phi0), r*sin(t+phi0), (p/(2*pi))*t)
  drop(H %*% t(R) + c)
}

.helix_newton <- function(t, Xloc, r, phi0, a, iters) {
  for (k in seq_len(iters)) {
    cx  <- r * cos(t + phi0); sx <- r * sin(t + phi0)
    fx  <- cbind(cx, sx, a * t) - Xloc
    d   <- cbind(-r * sin(t + phi0), r * cos(t + phi0), a)
    t   <- t - rowSums(d * fx) / (rowSums(d * d) + 1e-12)
  }
  t
}

.helix_d2 <- function(t, Xloc, r, phi0, a) {
  H <- cbind(r * cos(t + phi0), r * sin(t + phi0), a * t)
  rowSums((H - Xloc)^2)
}

.project_t_helix <- function(X, c, u, r, p, phi0, t0 = NULL, iters = 15,
                              check_neighbors = FALSE) {
  n <- nrow(X); R <- .rot_to_u(u); Xloc <- (X - matrix(c, n, 3, TRUE)) %*% R
  a <- p/(2*pi)
  if (is.null(t0)) {
    # Angular estimate: which angle on any turn?
    t_xy <- atan2(Xloc[,2], Xloc[,1]) - phi0
    if (abs(a) > 1e-12) {
      # Axial estimate: helix z-coordinate = a*t, so t ≈ z/a.
      # Use this to pick which turn (n_turns) and angular position for fine-tuning.
      # This prevents Newton from converging to a local minimum on the wrong turn.
      t_z     <- Xloc[,3] / a
      n_turns <- round((t_z - t_xy) / (2*pi))
      t0 <- t_xy + n_turns * 2*pi
    } else {
      t0 <- t_xy
    }
  }
  t <- .helix_newton(t0, Xloc, r, phi0, a, iters)

  # Check adjacent turns (±2π) and keep whichever gives shortest distance.
  # Handles cases where the axial/angular initial estimate lands on the wrong turn.
  if (check_neighbors && abs(a) > 1e-12) {
    t_m <- .helix_newton(t - 2*pi, Xloc, r, phi0, a, iters)
    t_p <- .helix_newton(t + 2*pi, Xloc, r, phi0, a, iters)
    d0  <- .helix_d2(t,   Xloc, r, phi0, a)
    dm  <- .helix_d2(t_m, Xloc, r, phi0, a)
    dp  <- .helix_d2(t_p, Xloc, r, phi0, a)
    t   <- ifelse(dm < pmin(d0, dp), t_m,
                  ifelse(dp < pmin(d0, dm), t_p, t))
  }
  t
}

# Coarse half-sphere orientation grid search.
# Sweeps n_theta × n_phi candidate (theta, phi) pairs covering the half-sphere
# (theta 0→π/2; the normal/axis and its negation define the same geometry).
# obj_fn(theta, phi) must return a scalar objective (lower = better).
# Returns the (theta, phi) pair with the lowest objective, or the PCA fallback
# supplied via th_init/ph_init if no grid candidate beats it.
.grid_search_orient <- function(obj_fn, th_init, ph_init, n_theta = 12L, n_phi = 16L) {
  # Strip any R names so the return vector has clean "theta"/"phi" keys.
  th_init <- unname(as.numeric(th_init))
  ph_init <- unname(as.numeric(ph_init))

  # cos-spaced theta rows → roughly equal solid-angle coverage
  th_vals <- acos(seq(1, 0, length.out = n_theta))   # 0 → π/2
  ph_vals <- seq(0, 2*pi, length.out = n_phi + 1L)[-(n_phi + 1L)]

  pca_val  <- tryCatch(obj_fn(th_init, ph_init), error = function(e) Inf)
  best_val <- pca_val
  best_th  <- th_init
  best_ph  <- ph_init

  for (th in th_vals) {
    for (ph in ph_vals) {
      v <- tryCatch(obj_fn(th, ph), error = function(e) Inf)
      if (is.finite(v) && v < best_val) {
        best_val <- v; best_th <- th; best_ph <- ph
      }
    }
  }
  c(theta = best_th, phi = best_ph)
}

mle_fit_helix <- function(X, w = NULL, turns_hint = 3, max_outer = 20,
                          center_init  = NULL,
                          center_fixed = NULL,
                          center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 6) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)
  
  # initiation
  c0 <- if (!is.null(center_init)) as.numeric(center_init) else colMeans(X)
  Xc <- sweep(X, 2, c0)
  pc <- prcomp(Xc, center = FALSE); u0 <- pc$rotation[,1]
  v1 <- pc$rotation[,2]; v2 <- pracma::cross(u0, v1)
  XY <- Xc %*% cbind(v1, v2, u0)
  r0 <- median(sqrt(XY[,1]^2 + XY[,2]^2))
  h  <- diff(range(XY[,3])); p0 <- if (turns_hint > 0) h/turns_hint else h
  par <- c(c0, th = acos(u0[3]), ph = atan2(u0[2], u0[1]), r = r0, p = p0, phi0 = 0)
  t   <- .project_t_helix(X, par[1:3], .ang_to_u(par[4],par[5]), par[6], par[7], par[8])

  # ── Orientation grid search: finds best axis direction before L-BFGS-B ────
  # Cheap radial proxy: distance from candidate axis should equal r0.
  # Avoids Newton projection; captures the right basin for the optimizer.
  c_gs <- if (!is.null(center_fixed)) as.numeric(center_fixed) else par[1:3]
  gs_h <- .grid_search_orient(
    obj_fn = function(th, ph) {
      ax  <- .ang_to_u(th, ph)
      Xc  <- sweep(X, 2, c_gs)
      d   <- as.numeric(Xc %*% ax)
      Xp  <- Xc - outer(d, ax)
      sum(w * (sqrt(rowSums(Xp^2)) - r0)^2) / n
    },
    th_init = par[4], ph_init = par[5]
  )
  if (gs_h[1] != par[4] || gs_h[2] != par[5]) {
    par[4] <- gs_h[1]; par[5] <- gs_h[2]
    t <- .project_t_helix(X, par[1:3], .ang_to_u(par[4], par[5]),
                           par[6], par[7], par[8])
  }

  for (outer in 1:max_outer) {
    obj <- function(pv) {
      # choose center according to mode
      c <- if (!is.null(center_fixed)) as.numeric(center_fixed) else pv[1:3]
      u <- .ang_to_u(pv[4],pv[5]); r <- abs(pv[6]); p <- pv[7]; phi0 <- pv[8]
      H <- t(vapply(t, .helix_point, numeric(3), c=c, u=u, r=r, p=p, phi0=phi0))
      res2 <- rowSums((X - H)^2)
      val <- sum(w * res2) / n   # mean squared distance - scale-independent of n
      if (is.null(center_fixed) && !is.null(center_prior) && lambda > 0) {
        val <- val + lambda * sum((c - center_prior)^2)
      }
      val
    }
    # Optimization dimension depends on whether c is fixed
    if (!is.null(center_fixed)) {
      par_opt <- par[-(1:3)]
      fit <- optim(par_opt, function(q) obj(c(par[1:3], q)),
                   method = "L-BFGS-B", control = list(maxit = 200))
      par[-(1:3)] <- fit$par
    } else {
      fit <- optim(par, obj, method = "L-BFGS-B", control = list(maxit = 200))
      par <- fit$par
    }
    c_now <- if (is.null(center_fixed)) par[1:3] else as.numeric(center_fixed)
    t_new <- .project_t_helix(X, c_now, .ang_to_u(par[4],par[5]), abs(par[6]), par[7], par[8], t0 = t)
    if (mean(abs(t_new - t)) < 1e-6) { t <- t_new; break } else t <- t_new
  }
  c <- if (is.null(center_fixed)) par[1:3] else as.numeric(center_fixed)
  u <- .ang_to_u(par[4],par[5]); r <- abs(par[6]); p <- par[7]; phi0 <- par[8]
  H <- t(vapply(t, .helix_point, numeric(3), c=c, u=u, r=r, p=p, phi0=phi0))
  list(center = c, axis = u, radius = r, pitch = p, phi0 = phi0, t = t, H = H)
}

mle_sample_helix <- function(fit, n = 3000, pad_turns = 0.25) {
  tmin <- min(fit$t) - 2*pi*pad_turns
  tmax <- max(fit$t) + 2*pi*pad_turns
  tt <- seq(tmin, tmax, length.out = n)
  H  <- t(vapply(tt, .helix_point, numeric(3),
                 c = fit$center, u = fit$axis, r = fit$radius, p = fit$pitch, phi0 = fit$phi0))
  as.data.frame(H)
}

# ── (D) Ring ─────────────────────────────────────────────────────────────────
# A 3D ring (circle) lives in a plane whose unit normal is n.
# Orthogonal distance from any point to the ring = sqrt((r_plane - R)^2 + d_ax^2)
# where r_plane = in-plane distance from the ring axis, d_ax = axial offset.

.ring_basis <- function(normal) {
  # Return [u | v] - two orthonormal vectors spanning the plane ⊥ to normal.
  n   <- as.numeric(normal) / sqrt(sum(normal^2))
  ref <- if (abs(n[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0)
  u   <- pracma::cross(ref, n); u <- u / sqrt(sum(u^2))
  v   <- pracma::cross(n, u)           # already unit (n ⊥ u, both unit)
  cbind(u, v)
}

mle_project_ring <- function(X, center, normal, radius) {
  n   <- as.numeric(normal) / sqrt(sum(normal^2))
  Xc  <- sweep(X, 2, center)
  d_ax <- as.numeric(Xc %*% n)
  Xp   <- Xc - outer(d_ax, n)          # in-plane component
  r_pl <- sqrt(rowSums(Xp^2))           # distance from ring axis in-plane

  # Unit in-plane direction; fall back to ring basis x-axis when on the axis
  basis  <- .ring_basis(n)
  safe   <- r_pl > 1e-12
  dir    <- matrix(basis[, 1], nrow(X), 3, TRUE)
  dir[safe, ] <- Xp[safe, , drop = FALSE] / r_pl[safe]

  nearest <- matrix(center, nrow(X), 3, TRUE) + radius * dir
  dist    <- sqrt((r_pl - radius)^2 + d_ax^2)
  list(nearest = nearest, dist = dist)
}

mle_fit_ring <- function(X, w = NULL,
                         center_init  = NULL,
                         center_fixed = NULL,
                         center_prior = NULL, lambda = 0) {
  n <- nrow(X); if (n < 4) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)

  geo  <- fit_ring_geometry(X)
  c0   <- if (!is.null(center_init)) as.numeric(center_init) else geo$center
  nrm0 <- geo$normal
  th0  <- acos(max(-1, min(1, nrm0[3])))
  ph0  <- atan2(nrm0[2], nrm0[1])
  r0   <- geo$radius
  r_max <- max(sqrt(rowSums(sweep(X, 2, c0)^2))) * 2

  # ── Orientation grid search: best normal direction before L-BFGS-B ────────
  c_gs <- if (!is.null(center_fixed)) as.numeric(center_fixed) else c0
  gs_r <- .grid_search_orient(
    obj_fn  = function(th, ph) sum(w * mle_project_ring(X, c_gs, .ang_to_u(th, ph), r0)$dist^2) / n,
    th_init = th0, ph_init = ph0
  )
  th0 <- gs_r[1]; ph0 <- gs_r[2]

  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand

  if (!is.null(center_fixed)) {
    cc   <- as.numeric(center_fixed)
    par0 <- c(th0, ph0, log(max(r0, 1e-3)))
    obj  <- function(par) {
      pr <- mle_project_ring(X, cc, .ang_to_u(par[1], par[2]), exp(par[3]))
      sum(w * pr$dist^2) / n
    }
    fit <- optim(par0, obj, method = "L-BFGS-B",
                 lower = c(-Inf, -Inf, log(1e-3)),
                 upper = c( Inf,  Inf, log(r_max)),
                 control = list(maxit = 500))
    p <- fit$par
    return(list(center = cc, normal = .ang_to_u(p[1], p[2]), radius = exp(p[3])))
  }

  par0 <- c(c0, th0, ph0, log(max(r0, 1e-3)))
  obj  <- function(par) {
    pr  <- mle_project_ring(X, par[1:3], .ang_to_u(par[4], par[5]), exp(par[6]))
    val <- sum(w * pr$dist^2) / n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((par[1:3] - center_prior)^2)
    val
  }
  fit <- optim(par0, obj, method = "L-BFGS-B",
               lower = c(c_lo, -Inf, -Inf, log(1e-3)),
               upper = c(c_hi,  Inf,  Inf, log(r_max)),
               control = list(maxit = 500))
  p <- fit$par
  list(center = p[1:3], normal = .ang_to_u(p[4], p[5]), radius = exp(p[6]))
}

mle_sample_ring <- function(fit, n = 500) {
  basis <- .ring_basis(fit$normal)
  tt    <- seq(0, 2*pi, length.out = n + 1L)[-(n + 1L)]
  pts   <- fit$radius * (outer(cos(tt), basis[, 1]) + outer(sin(tt), basis[, 2])) +
           matrix(fit$center, n, 3, TRUE)
  colnames(pts) <- c("x", "y", "z")
  as.data.frame(pts)
}

# ── (E) Arc ───────────────────────────────────────────────────────────────────
# An arc is a partial ring: same center/normal/radius, plus t_center (midpoint
# angle in the ring-plane frame) and arc_half (half-span, 0 < arc_half <= π).
# The arc runs from t_center − arc_half to t_center + arc_half.

mle_project_arc <- function(X, center, normal, radius, t_center, arc_half) {
  n     <- as.numeric(normal) / sqrt(sum(normal^2))
  basis <- .ring_basis(n)
  u <- basis[, 1]; vv <- basis[, 2]

  Xc   <- sweep(X, 2, center)
  d_ax <- as.numeric(Xc %*% n)
  Xp   <- Xc - outer(d_ax, n)

  # Angle of each in-plane projection in ring frame
  t_ang <- atan2(as.numeric(Xp %*% vv), as.numeric(Xp %*% u))

  # Wrap relative angle to [−π, π] then clamp to the arc half-span
  dt      <- ((t_ang - t_center + pi) %% (2*pi)) - pi
  dt_cl   <- pmax(-arc_half, pmin(arc_half, dt))
  t_near  <- t_center + dt_cl

  nearest <- radius * (outer(cos(t_near), u) + outer(sin(t_near), vv)) +
             matrix(center, nrow(X), 3, TRUE)
  dist    <- sqrt(rowSums((X - nearest)^2))
  list(nearest = nearest, dist = dist)
}

mle_fit_arc <- function(X, w = NULL,
                        center_init  = NULL,
                        center_fixed = NULL,
                        center_prior = NULL, lambda = 0,
                        arc_span_hint = 0) {
  # The arc_half parameter has a monotone incentive to grow (more arc coverage
  # = lower sum-of-squared-distances).  To prevent collapse to a full ring we
  # use a two-stage strategy:
  #   Stage 1 – fix arc_half at the heuristic estimate; optimise center /
  #              normal / radius / t_center.  This anchors the arc geometry.
  #   Stage 2 – release arc_half but bound it to ±40 % of the heuristic value;
  #              fine-tune all parameters together.
  n <- nrow(X); if (n < 4) return(NULL)
  if (is.null(w)) w <- rep(1, n); w <- w / mean(w)

  geo  <- fit_arc_geometry(X, arc_span_hint = arc_span_hint)
  c0   <- if (!is.null(center_init)) as.numeric(center_init) else geo$center
  nrm0 <- geo$normal
  th0  <- acos(max(-1, min(1, nrm0[3])))
  ph0  <- atan2(nrm0[2], nrm0[1])
  r0   <- geo$radius
  tc0  <- geo$t_center
  ah0  <- geo$arc_half          # heuristic half-span; used as anchor
  r_max <- max(sqrt(rowSums(sweep(X, 2, c0)^2))) * 2

  # ── Orientation grid search: best normal direction before L-BFGS-B ────────
  # Uses the arc projection with tc0/ah0 fixed - only orientation is swept.
  c_gs <- if (!is.null(center_fixed)) as.numeric(center_fixed) else c0
  gs_a <- .grid_search_orient(
    obj_fn  = function(th, ph)
      sum(w * mle_project_arc(X, c_gs, .ang_to_u(th, ph), r0, tc0, ah0)$dist^2) / n,
    th_init = th0, ph_init = ph0
  )
  th0 <- gs_a[1]; ph0 <- gs_a[2]

  # arc_half bounds: ±40 % of heuristic, clamped to [1°, 179°]
  ah_lo_log <- log(pmax(ah0 * 0.60, 0.02))
  ah_hi_log <- log(pmin(ah0 * 1.40, pi - 1e-6))

  X_range    <- apply(X, 2, range)
  box_expand <- pmax((X_range[2,] - X_range[1,]) * 0.5, 10)
  c_lo <- X_range[1,] - box_expand
  c_hi <- X_range[2,] + box_expand

  # Shared objective: fixed arc_half supplied externally
  .obj_fixed_ah <- function(par, cc, ah) {
    pr  <- mle_project_arc(X, cc, .ang_to_u(par[1], par[2]), exp(par[3]), par[4], ah)
    sum(w * pr$dist^2) / n
  }
  # Full objective with arc_half as par[5]
  .obj_full_cc <- function(par, cc) {
    ah  <- pmin(exp(par[5]), pi - 1e-6)
    pr  <- mle_project_arc(X, cc, .ang_to_u(par[1], par[2]), exp(par[3]), par[4], ah)
    sum(w * pr$dist^2) / n
  }

  if (!is.null(center_fixed)) {
    cc <- as.numeric(center_fixed)
    # Stage 1
    p1 <- optim(c(th0, ph0, log(max(r0, 1e-3)), tc0),
                .obj_fixed_ah, cc = cc, ah = ah0, method = "L-BFGS-B",
                lower = c(-Inf, -Inf, log(1e-3), -Inf),
                upper = c( Inf,  Inf, log(r_max),  Inf),
                control = list(maxit = 300))$par
    # Stage 2
    p2 <- optim(c(p1, log(ah0)), .obj_full_cc, cc = cc, method = "L-BFGS-B",
                lower = c(-Inf, -Inf, log(1e-3), -Inf, ah_lo_log),
                upper = c( Inf,  Inf, log(r_max),  Inf, ah_hi_log),
                control = list(maxit = 300))$par
    return(list(center = cc, normal = .ang_to_u(p2[1], p2[2]), radius = exp(p2[3]),
                t_center = p2[4], arc_half = pmin(exp(p2[5]), pi - 1e-6)))
  }

  # Free centre ── Stage 1: fix arc_half
  obj1 <- function(par) {
    pr  <- mle_project_arc(X, par[1:3], .ang_to_u(par[4], par[5]),
                           exp(par[6]), par[7], ah0)
    val <- sum(w * pr$dist^2) / n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((par[1:3] - center_prior)^2)
    val
  }
  f1 <- optim(c(c0, th0, ph0, log(max(r0, 1e-3)), tc0), obj1,
              method = "L-BFGS-B",
              lower = c(c_lo, -Inf, -Inf, log(1e-3), -Inf),
              upper = c(c_hi,  Inf,  Inf, log(r_max),  Inf),
              control = list(maxit = 300))$par

  # Stage 2: release arc_half within ±40 % of heuristic
  obj2 <- function(par) {
    ah  <- pmin(exp(par[8]), pi - 1e-6)
    pr  <- mle_project_arc(X, par[1:3], .ang_to_u(par[4], par[5]),
                           exp(par[6]), par[7], ah)
    val <- sum(w * pr$dist^2) / n
    if (!is.null(center_prior) && lambda > 0)
      val <- val + lambda * sum((par[1:3] - center_prior)^2)
    val
  }
  f2 <- optim(c(f1, log(ah0)), obj2, method = "L-BFGS-B",
              lower = c(c_lo, -Inf, -Inf, log(1e-3), -Inf, ah_lo_log),
              upper = c(c_hi,  Inf,  Inf, log(r_max),  Inf, ah_hi_log),
              control = list(maxit = 300))$par

  list(center = f2[1:3], normal = .ang_to_u(f2[4], f2[5]), radius = exp(f2[6]),
       t_center = f2[7], arc_half = pmin(exp(f2[8]), pi - 1e-6))
}

mle_sample_arc <- function(fit, n = 300) {
  basis   <- .ring_basis(fit$normal)
  t_start <- fit$t_center - fit$arc_half
  t_end   <- fit$t_center + fit$arc_half
  tt      <- seq(t_start, t_end, length.out = n)
  pts     <- fit$radius * (outer(cos(tt), basis[, 1]) + outer(sin(tt), basis[, 2])) +
             matrix(fit$center, n, 3, TRUE)
  colnames(pts) <- c("x", "y", "z")
  as.data.frame(pts)
}

## ── Tube mesh helpers for helix, ring, and arc ─────────────────────────────
# Internal: given centerline positions and two cross-section frame vectors
# (frame1, frame2, each n_along × 3), build an n_along × n_tube vertex grid
# and return a mesh3d-compatible list (x, y, z, i, j, k - 0-indexed for plotly).
.tube_mesh_from_frames <- function(positions, frame1, frame2, tube_r,
                                   n_tube = 16, wrap_along = FALSE) {
  n_along <- nrow(positions)
  phi <- seq(0, 2*pi, length.out = n_tube + 1L)[-(n_tube + 1L)]

  vert <- matrix(NA_real_, n_along * n_tube, 3L)
  for (a in seq_len(n_along)) {
    base <- (a - 1L) * n_tube
    for (t in seq_len(n_tube)) {
      vert[base + t, ] <- positions[a, ] +
        tube_r * (cos(phi[t]) * frame1[a, ] + sin(phi[t]) * frame2[a, ])
    }
  }

  vidx <- function(a, t) (a - 1L) * n_tube + (t - 1L)   # 0-based for plotly
  n_a_max <- if (wrap_along) n_along else n_along - 1L
  ii <- integer(n_a_max * n_tube * 2L)
  jj <- integer(n_a_max * n_tube * 2L)
  kk <- integer(n_a_max * n_tube * 2L)
  k_out <- 1L
  for (a in seq_len(n_a_max)) {
    a2 <- if (a == n_along) 1L else a + 1L
    for (t in seq_len(n_tube)) {
      t2 <- if (t == n_tube) 1L else t + 1L
      ii[k_out]     <- vidx(a,  t);  jj[k_out]     <- vidx(a2, t);  kk[k_out]     <- vidx(a,  t2)
      ii[k_out + 1] <- vidx(a2, t);  jj[k_out + 1] <- vidx(a2, t2); kk[k_out + 1] <- vidx(a,  t2)
      k_out <- k_out + 2L
    }
  }

  list(x = vert[, 1], y = vert[, 2], z = vert[, 3], i = ii, j = jj, k = kk)
}

# Torus mesh for a ring fit (closed tube around a 3D circle).
make_ring_tube_mesh <- function(center, normal, radius, tube_r,
                                n_ring = 64, n_tube = 16) {
  basis <- .ring_basis(normal)
  u <- basis[, 1]; v <- basis[, 2]
  n_vec <- as.numeric(normal) / sqrt(sum(normal^2))
  theta  <- seq(0, 2*pi, length.out = n_ring + 1L)[-(n_ring + 1L)]
  cos_t  <- cos(theta); sin_t <- sin(theta)
  pos    <- matrix(center, n_ring, 3L, TRUE) +
            radius * (outer(cos_t, u) + outer(sin_t, v))
  radial <- outer(cos_t, u) + outer(sin_t, v)       # outward from ring axis
  axial  <- matrix(n_vec, n_ring, 3L, TRUE)
  .tube_mesh_from_frames(pos, radial, axial, tube_r, n_tube = n_tube, wrap_along = TRUE)
}

# Open-ended tube for an arc fit.
make_arc_tube_mesh <- function(center, normal, radius, t_center, arc_half, tube_r,
                               n_arc = 64, n_tube = 16) {
  basis <- .ring_basis(normal)
  u <- basis[, 1]; v <- basis[, 2]
  n_vec  <- as.numeric(normal) / sqrt(sum(normal^2))
  theta  <- seq(t_center - arc_half, t_center + arc_half, length.out = n_arc)
  cos_t  <- cos(theta); sin_t <- sin(theta)
  pos    <- matrix(center, n_arc, 3L, TRUE) +
            radius * (outer(cos_t, u) + outer(sin_t, v))
  radial <- outer(cos_t, u) + outer(sin_t, v)
  axial  <- matrix(n_vec, n_arc, 3L, TRUE)
  .tube_mesh_from_frames(pos, radial, axial, tube_r, n_tube = n_tube, wrap_along = FALSE)
}

# Open-ended tube for a helix fit.
# center/axis/radius/pitch/phi0 are the helix parameters as used by .helix_point.
# tmin/tmax are the parameter range to render.
make_helix_tube_mesh <- function(center, axis, radius, pitch, phi0, tmin, tmax,
                                 tube_r, n_helix = 200, n_tube = 12) {
  tt  <- seq(tmin, tmax, length.out = n_helix)
  pos <- t(vapply(tt, .helix_point, numeric(3),
                  c = center, u = axis, r = radius, p = pitch, phi0 = phi0))
  ax  <- axis / sqrt(sum(axis^2))

  frame1 <- matrix(NA_real_, n_helix, 3L)
  frame2 <- matrix(NA_real_, n_helix, 3L)
  for (i in seq_len(n_helix)) {
    # Outward radial: component of (pos - center) perpendicular to axis
    diff_i   <- pos[i, ] - center
    p_axial  <- sum(diff_i * ax) * ax
    p_radial <- diff_i - p_axial
    r_len    <- sqrt(sum(p_radial^2))
    f1 <- if (r_len > 1e-12) p_radial / r_len else {
      ref <- if (abs(ax[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0)
      tmp <- pracma::cross(ax, ref); tmp / sqrt(sum(tmp^2))
    }
    frame1[i, ] <- f1
    frame2[i, ] <- pracma::cross(ax, f1)
  }
  .tube_mesh_from_frames(pos, frame1, frame2, tube_r, n_tube = n_tube, wrap_along = FALSE)
}

add_distance_lines_to_plot <- function(fig, medians, max_dist = Inf) {
  required_cols <- c("x", "y", "z", "nearest_x", "nearest_y", "nearest_z", "distance_to_fit")
  if (!all(required_cols %in% colnames(medians))) return(fig)

  # Drop rows where either endpoint is missing/non-finite; apply length threshold
  valid <- is.finite(medians$x) & is.finite(medians$y) & is.finite(medians$z) &
           is.finite(medians$nearest_x) & is.finite(medians$nearest_y) & is.finite(medians$nearest_z) &
           is.finite(medians$distance_to_fit)
  if (is.finite(max_dist) && max_dist > 0) valid <- valid & (medians$distance_to_fit <= max_dist)
  medians <- medians[valid, , drop = FALSE]
  if (nrow(medians) == 0) return(fig)

  # Build NA-separated x/y/z for a single line trace (each triplet = one segment)
  xs <- as.vector(rbind(medians$x, medians$nearest_x, NA_real_))
  ys <- as.vector(rbind(medians$y, medians$nearest_y, NA_real_))
  zs <- as.vector(rbind(medians$z, medians$nearest_z, NA_real_))

  fig <- fig %>% add_trace(
    x = xs, y = ys, z = zs,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "blue", width = 2),
    showlegend = FALSE,
    inherit = FALSE
  )

  # Text labels at midpoints - one trace for all labels
  fig <- fig %>% add_trace(
    x = (medians$x + medians$nearest_x) / 2,
    y = (medians$y + medians$nearest_y) / 2,
    z = (medians$z + medians$nearest_z) / 2,
    type = "scatter3d",
    mode = "text",
    text = as.character(round(medians$distance_to_fit, 1)),
    textfont = list(color = "blue", size = 10),
    showlegend = FALSE,
    inherit = FALSE
  )

  return(fig)
}

compute_convex_hull_centroid <- function(coords_matrix) {
  if (nrow(coords_matrix) < 4 || ncol(coords_matrix) != 3) {
    warning("Invalid input to convex hull centroid computation.")
    return(colMeans(coords_matrix))
  }
  
  result <- tryCatch({
    hull <- geometry::convhulln(coords_matrix, options = "FA")
    if (is.null(hull$hull)) stop("Hull computation failed")
    
    triangles <- hull$hull
    pts <- coords_matrix
    
    tri_centers <- matrix(0, nrow = nrow(triangles), ncol = 3)
    tri_areas <- numeric(nrow(triangles))
    
    for (i in seq_len(nrow(triangles))) {
      tri_idx <- triangles[i, ]
      if (length(tri_idx) != 3) next  # skip invalid
      tri_pts <- pts[tri_idx, , drop = FALSE]
      if (nrow(tri_pts) != 3 || ncol(tri_pts) != 3) next  # skip if malformed
      
      v1 <- tri_pts[2, ] - tri_pts[1, ]
      v2 <- tri_pts[3, ] - tri_pts[1, ]
      if (length(v1) != 3 || length(v2) != 3) next  # safety
      
      centroid <- colMeans(tri_pts)
      area <- 0.5 * sqrt(sum(crossprod(pracma::cross(v1, v2))))
      
      tri_centers[i, ] <- centroid
      tri_areas[i] <- area
    }
    
    if (sum(tri_areas) == 0) stop("Degenerate triangles")
    
    center <- colSums(tri_centers * tri_areas) / sum(tri_areas)
    center
  }, error = function(e) {
    warning("Convex hull fallback: ", conditionMessage(e))
    return(colMeans(coords_matrix))  # fallback
  })
  
  return(result)
}

medianize_coords <- function(df, use_cluster1 = FALSE) {
  if (use_cluster1 &&
      "cluster" %in% names(df) &&
      any(is.finite(df$cluster) & df$cluster != 0, na.rm = TRUE)) {
    df %>%
      dplyr::filter(is.finite(cluster) & cluster != 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(x = median(x), y = median(y), z = median(z), .groups = "drop") %>%
      dplyr::select(x, y, z)
  } else {
    df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = median(x), y = median(y), z = median(z), .groups = "drop") %>%
      dplyr::select(x, y, z)
  }
}

# ---- lightweight geometry helpers (parameters only, no random sampling) ----
spread_from_group <- function(x, y, z) {
  if (length(x) <= 1) return(0)
  mx <- stats::median(x, na.rm = TRUE)
  my <- stats::median(y, na.rm = TRUE)
  mz <- stats::median(z, na.rm = TRUE)
  sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2, na.rm = TRUE))
}

fit_sphere_geometry <- function(coords_matrix) {
  center <- tryCatch(geometric_median(coords_matrix), error = function(e) colMeans(coords_matrix))
  r      <- stats::median(sqrt(rowSums(sweep(coords_matrix, 2, center)^2)))
  list(center = center, radius = r)
}

fit_ellipsoid_geometry <- function(coords_matrix) {
  center <- tryCatch(geometric_median(coords_matrix), error = function(e) colMeans(coords_matrix))
  Xc     <- sweep(coords_matrix, 2, center)   # center before PCA so axes are aligned to the ellipsoid center
  pca    <- prcomp(Xc, center = FALSE, scale. = FALSE)
  rot    <- Xc %*% pca$rotation
  d <- diff(stats::quantile(rot[, 1], c(0.05, 0.95))) / 2
  e <- diff(stats::quantile(rot[, 2], c(0.05, 0.95))) / 2
  f <- diff(stats::quantile(rot[, 3], c(0.05, 0.95))) / 2
  sf <- stats::quantile(sqrt(rowSums(rot^2)), 0.90) / max(d, e, f)
  d <- d * sf; e <- e * sf; f <- f * sf
  if (min(d, e, f) < 1e-6) return(NULL)
  list(center = center, axes = c(d, e, f), R = pca$rotation)
}

fit_helix_geometry <- function(coords_matrix, turns = 3, n_points = 5000) {
  center <- tryCatch(geometric_median(coords_matrix), error = function(e) colMeans(coords_matrix))
  Xc     <- sweep(coords_matrix, 2, center)
  pca    <- prcomp(Xc, center = FALSE)
  a1 <- pca$rotation[, 1]; a2 <- pca$rotation[, 2]
  a3 <- pracma::cross(a1, a2)
  ht <- diff(range(Xc %*% a1))
  br <- stats::quantile(abs(Xc %*% a2), 0.9, na.rm = TRUE)
  t  <- seq(0, 2 * pi * max(1, turns), length.out = n_points)
  Hloc <- cbind(br * cos(t), br * sin(t),
                seq(-ht / 2, ht / 2, length.out = n_points))
  H <- sweep(Hloc %*% t(cbind(a2, a3, a1)), 2, center, "+")
  colnames(H) <- c("x", "y", "z")

  # Parameters compatible with .project_t_helix / .helix_point so that
  # Newton projection operates on the exact same 3D curve as H above.
  # .rot_to_u(a1) uses a different x/y-axis convention than cbind(a2,a3,a1),
  # so phi0 aligns the phase; proj_center shifts z-origin from -ht/2 to 0.
  #
  # Derivation: .helix_point(t=0) = proj_center + r*(cos(phi0)*R[,1]+sin(phi0)*R[,2])
  # must equal H[1,] = proj_center + r*a2.
  # Projecting a2 onto the {R[,1], R[,2]} Newton frame gives:
  #   cos(phi0) = a2 . R[,1]  →  sin(phi0) = a2 . R[,2]
  R_nwt       <- .rot_to_u(a1)
  phi0        <- atan2(sum(a2 * R_nwt[, 2]), sum(a2 * R_nwt[, 1]))
  proj_center <- center - a1 * (ht / 2)

  list(points = H, center = center,
       axis = a1, radius = br, pitch = ht / max(1, turns),
       proj_center = proj_center, phi0 = phi0)
}

fit_ring_geometry <- function(coords_matrix) {
  center <- tryCatch(geometric_median(coords_matrix),
                     error = function(e) colMeans(coords_matrix))
  Xc     <- sweep(coords_matrix, 2, center)
  pca    <- prcomp(Xc, center = FALSE)
  # PC3 = direction of smallest variance = normal to the ring plane
  normal <- pca$rotation[, 3]

  d_ax   <- as.numeric(Xc %*% normal)
  Xp     <- Xc - outer(d_ax, normal)      # in-plane projections
  radius <- stats::median(sqrt(rowSums(Xp^2)))
  list(center = center, normal = normal, radius = radius)
}

fit_arc_geometry <- function(coords_matrix, arc_span_hint = 0) {
  geo   <- fit_ring_geometry(coords_matrix)
  basis <- .ring_basis(geo$normal)
  u <- basis[, 1]; vv <- basis[, 2]

  Xc    <- sweep(coords_matrix, 2, geo$center)
  d_ax  <- as.numeric(Xc %*% geo$normal)
  Xp    <- Xc - outer(d_ax, geo$normal)
  t_ang <- atan2(as.numeric(Xp %*% vv), as.numeric(Xp %*% u))

  # Circular mean → arc midpoint angle
  t_center <- atan2(mean(sin(t_ang)), mean(cos(t_ang)))

  arc_half <- if (arc_span_hint > 1) {
    pmin((arc_span_hint / 2) * pi / 180, pi - 1e-6)
  } else {
    # Auto: 97th-percentile angular deviation from the arc midpoint + 5% buffer
    dt <- ((t_ang - t_center + pi) %% (2*pi)) - pi
    pmin(stats::quantile(abs(dt), 0.97, names = FALSE) * 1.05, pi - 1e-6)
  }
  list(center = geo$center, normal = geo$normal, radius = geo$radius,
       t_center = t_center, arc_half = arc_half)
}

generate_ring_from_cluster2 <- function(cluster_data, n_points = 500,
                                        use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  geo <- fit_ring_geometry(as.matrix(coords))
  pts <- mle_sample_ring(geo, n = n_points)
  list(points = pts, center = geo$center, normal = geo$normal, radius = geo$radius)
}

generate_arc_from_cluster2 <- function(cluster_data, n_points = 300,
                                       arc_span_hint = 0, use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  geo <- fit_arc_geometry(as.matrix(coords), arc_span_hint = arc_span_hint)
  pts <- mle_sample_arc(geo, n = n_points)
  list(points = pts, center = geo$center, normal = geo$normal, radius = geo$radius,
       t_center = geo$t_center, arc_half = geo$arc_half)
}

# =============================================================================
# alpha_surface_3d() - WRAP ONE GROUP OF LOCALIZATIONS IN AN ALPHA SHAPE,
#                      RADIUS CHOSEN FROM THE DATA
# =============================================================================
#
# GOAL
#   Given the localizations of ONE group (a trace, a de-blinked molecule, a
#   Median Cluster), return the closed triangular surface that wraps it, ready
#   to hand straight to plotly's mesh3d. Callers only supply points; there is
#   deliberately no radius argument and no UI control for one.
#
# WHY THERE IS NO ALPHA RADIUS CONTROL ANY MORE (this was a real bug)
#   `alphashape3d::ashape3d()` computes ONE Delaunay tetrahedralization and then,
#   for each requested alpha, adds a COLUMN classifying every triangle. The
#   triangle list in `$triang[, 1:3]` is therefore identical for every alpha -
#   alpha lives only in the classification columns. The old code drew
#   `$triang[, 1:3]` wholesale, i.e. the entire tetrahedralization, so the mesh
#   was byte-for-byte the same at every slider position. That is why moving the
#   "Alpha radius" slider visibly did nothing.
#
#   Column layout of `$triang` (verified against alphashape3d):
#     1-3  vertex indices        4  on.ch (on convex hull)   5  attached
#     6-8  rhoT / muT / MuT      9..  one class per alpha:
#          0 = exterior, 1 = singular, 2 = REGULAR, 3 = interior
#   Only REGULAR triangles form the closed boundary surface. (Check: as alpha
#   grows the regular set converges exactly onto the convex-hull faces, while
#   the singular count explodes - singular triangles are dangling slivers and
#   drawing them is what produced the old opaque blob.)
#
# HOW THE RADIUS IS PICKED
#   A ladder of candidate radii is built from the group's OWN spacing (multiples
#   of its median nearest-neighbour distance) and classified in a SINGLE
#   ashape3d call - extra alphas cost one column each, not another
#   tetrahedralization. We then take the SMALLEST radius that passes two tests:
#     (1) every point sits in ONE connected component (nothing shed), and
#     (2) the regular triangles form a CLOSED surface - each edge used by
#         exactly two of them. Test (2) matters because connectivity alone still
#         admits radii whose boundary has holes, and plotly renders a holed mesh
#         as a torn, see-through shell.
#   That is the tightest honest wrap: smaller and the surface fragments, larger
#   and detail is lost for nothing. It is scale-free, so a 30 nm de-blinked
#   molecule and a 500 nm structure are each wrapped correctly without anyone
#   touching a control. Falls back to the convex hull (largest rung, then
#   geometry::convhulln) if no candidate ever qualifies.
#
# ARGUMENTS
#   pts_mat : n x 3 numeric matrix of coordinates (nm)
#   max_pts : subsample cap. Meshing cost grows steeply with n and a surface is
#             not made more truthful by tens of thousands of vertices.
# RETURNS
#   list(x, y, z, i, j, k, alpha, n_tri, n_pts) with i/j/k already 0-BASED for
#   plotly mesh3d and x/y/z the (possibly subsampled) vertices they index, or
#   NULL when the group cannot be wrapped (<4 usable points, or no triangle
#   survives). Callers MUST use the returned x/y/z, not their own copy of the
#   points, because of the subsampling.
# =============================================================================
alpha_surface_3d <- function(pts_mat, max_pts = 2000L) {
  if (!is.matrix(pts_mat)) pts_mat <- as.matrix(pts_mat)     # accept data.frame
  if (ncol(pts_mat) < 3L) return(NULL)
  pts_mat <- pts_mat[, 1:3, drop = FALSE]
  storage.mode(pts_mat) <- "double"
  pts_mat <- pts_mat[stats::complete.cases(pts_mat), , drop = FALSE]
  if (nrow(pts_mat) < 4L) return(NULL)                        # a tetrahedron is the minimum
  # Deterministic subsample (sorted) so the mesh does not flicker between renders.
  if (nrow(pts_mat) > max_pts)
    pts_mat <- pts_mat[sort(sample.int(nrow(pts_mat), max_pts)), , drop = FALSE]

  # A group that is flat in any axis (2D data promoted to 3D, or a trace whose z
  # never moved) has zero volume and makes the tetrahedralization degenerate.
  # Nudge by sub-picometre amounts: invisible at nm scale, enough for the solver.
  for (k in 1:3) {
    sdk <- stats::sd(pts_mat[, k], na.rm = TRUE)
    if (!is.finite(sdk) || sdk < 1e-9)
      pts_mat[, k] <- pts_mat[, k] + seq_len(nrow(pts_mat)) * 1e-6
  }

  # --- candidate radii, expressed in the group's own units of spacing ---------
  nnd <- tryCatch({
    d <- RANN::nn2(pts_mat, pts_mat, k = 2L)$nn.dists[, 2]    # distance to nearest OTHER point
    d[is.finite(d) & d > 0]
  }, error = function(e) numeric(0))
  # Bounding-box diagonal bounds the useful range: an alpha far above it can only
  # ever give the convex hull, so it is the natural top rung.
  diag_nm <- sqrt(sum((apply(pts_mat, 2, max) - apply(pts_mat, 2, min))^2))
  base <- if (length(nnd)) stats::median(nnd) else diag_nm / 10
  if (!is.finite(base) || base <= 0) base <- diag_nm / 10
  if (!is.finite(base) || base <= 0) return(NULL)             # all points coincide
  alphas <- base * c(1, 1.5, 2, 3, 4, 6, 8, 12, 20, 40)
  alphas <- c(alphas[alphas < diag_nm], 10 * diag_nm)         # last rung = convex hull
  alphas <- sort(unique(alphas[is.finite(alphas) & alphas > 0]))
  if (!length(alphas)) return(NULL)

  # Degenerate inputs (n == 4, all-coincident points) make ashape3d emit internal
  # warnings on its way to returning nothing useful; they are not actionable.
  as3d <- suppressWarnings(tryCatch(alphashape3d::ashape3d(pts_mat, alpha = alphas),
                                    error = function(e) NULL))

  # Is this triangle set a closed surface? Every edge of a closed, watertight
  # mesh belongs to exactly two faces; an edge used once is a hole rim.
  is_closed <- function(tt) {
    if (is.null(tt) || nrow(tt) < 4L) return(FALSE)
    e <- rbind(tt[, c(1, 2)], tt[, c(2, 3)], tt[, c(1, 3)])
    e <- cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2]))    # undirected edge key
    all(tabulate(match(paste(e[, 1], e[, 2]), unique(paste(e[, 1], e[, 2])))) == 2L)
  }

  tt <- NULL; alpha_used <- NA_real_
  if (!is.null(as3d) && !is.null(as3d$triang) && nrow(as3d$triang)) {
    tri  <- as.matrix(as3d$triang)
    # The per-alpha classification columns are the LAST length(alphas) columns, in
    # the order the alphas were passed. Indexing from the end keeps this correct
    # if the package ever adds a descriptive column at the front.
    cls0 <- ncol(tri) - length(alphas)
    # components_ashape3d() labels each point with its connected component
    # (1 = largest, < 0 for points isolated at this radius). One component
    # covering every point == "wraps the whole group in one piece".
    comp <- suppressWarnings(tryCatch(
      alphashape3d::components_ashape3d(as3d, indexAlpha = "all"),
      error = function(e) NULL))
    if (!is.null(comp) && !is.list(comp)) comp <- list(comp)
    if (cls0 >= 3L) {
      for (ii in seq_along(alphas)) {
        cc <- if (!is.null(comp) && length(comp) == length(alphas)) comp[[ii]] else NULL
        # Test (1): nothing shed. Skipped only if components could not be computed.
        if (!is.null(cc) && !(length(cc) && all(is.finite(cc)) && all(cc == 1L))) next
        cand <- tri[tri[, cls0 + ii] == 2L, 1:3, drop = FALSE]  # REGULAR = the boundary
        if (!nrow(cand)) next
        # Test (2): watertight. The last rung is the convex hull, always accepted.
        if (is_closed(cand) || ii == length(alphas)) {
          tt <- cand; alpha_used <- alphas[ii]; break
        }
      }
    }
  }
  # Last resort: the convex hull. Correct (if loose) for any point set, and the
  # only thing that works for the n == 4 / near-coincident cases ashape3d rejects.
  if (is.null(tt) || !nrow(tt)) {
    tt <- tryCatch(as.matrix(geometry::convhulln(pts_mat, options = "Qt")[, 1:3, drop = FALSE]),
                   error = function(e) NULL)
    alpha_used <- Inf
  }
  if (is.null(tt) || !nrow(tt)) return(NULL)

  list(x = pts_mat[, 1], y = pts_mat[, 2], z = pts_mat[, 3],
       i = tt[, 1] - 1L, j = tt[, 2] - 1L, k = tt[, 3] - 1L,  # plotly mesh3d is 0-based
       alpha = alpha_used, n_tri = nrow(tt), n_pts = nrow(pts_mat))
}

# NOTE: `virtual_points` is a dense uniform sampling of the fitted SPHERE, not
# data. Its triangulation is the convex hull by construction, so no alpha shape
# is built here - the caller hulls it. (An `ashape` element used to be returned
# alongside; it was never actually drawn.)
generate_sphere_from_cluster2 <- function(cluster_data, n_points = 4000,
                                          use_cluster1 = FALSE) {

  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  coords_matrix <- as.matrix(coords)

  geo <- fit_sphere_geometry(coords_matrix)

  phi <- runif(n_points, 0, 2*pi)
  theta <- acos(runif(n_points, -1, 1))
  virtual_points <- cbind(
    geo$radius * sin(theta) * cos(phi),
    geo$radius * sin(theta) * sin(phi),
    geo$radius * cos(theta)
  ) + matrix(geo$center, n_points, 3, TRUE)

  list(center = geo$center, radius = geo$radius, virtual_points = virtual_points)
}

# Same as generate_sphere_from_cluster2: a dense sampling of the fitted (convex)
# ELLIPSOID, hulled by the caller. No alpha shape involved.
generate_ellipsoid_from_cluster2 <- function(cluster_data, n_points = 4000,
                                             use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 4) return(NULL)
  coords_matrix <- as.matrix(coords)

  geo <- fit_ellipsoid_geometry(coords_matrix)
  if (is.null(geo)) return(NULL)

  phi <- runif(n_points, 0, 2*pi)
  theta <- acos(runif(n_points, -1, 1))
  ellipsoid_points <- cbind(
    geo$axes[1] * sin(theta) * cos(phi),
    geo$axes[2] * sin(theta) * sin(phi),
    geo$axes[3] * cos(theta)
  ) %*% t(geo$R) + matrix(geo$center, n_points, 3, TRUE)

  list(center = geo$center, axes = geo$axes, R = geo$R,
       virtual_points = ellipsoid_points)
}

generate_helix_from_cluster2 <- function(cluster_data, n_points = 10000, turns = 3,
                                         use_cluster1 = FALSE) {
  coords <- medianize_coords(cluster_data, use_cluster1)
  if (nrow(coords) < 5) return(NULL)
  coords_matrix <- as.matrix(coords)

  geo <- fit_helix_geometry(coords_matrix, turns = turns, n_points = n_points)
  list(points = as.data.frame(geo$points), center = geo$center)
}

