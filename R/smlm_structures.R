## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_structures.R  -  Structures: planes, surfaces, isosurfaces, curvature
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Reference planes and surfaces, the Median Cluster lines/angles overlay,
##   the cluster density field and its isosurface, alpha volumes, and
##   curvature per Median Cluster.
##
## CONTENTS (in file order)
##   fit_poly_surface, add_cluster2_lines_and_angles,
##   compute_cluster2_plane_distances, generate_xy_plane,
##   create_global_surface, cluster2_member_col, SIGMA_GAP_FACTOR,
##   cluster_density_field, density_field_mesh, iso_channel_color,
##   add_cluster_isosurface, iso_mesh_trace, iso_status_ui, isTRUE_vec,
##   cluster_alpha_volume, estimate_cluster2_curvature, build_curvature_table,
##   .viewer_curve_points, .viewer_surface_mesh
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

fit_poly_surface <- function(df, degree = 2, grid_n = 50) {
  stopifnot(all(c("x","y","z") %in% names(df)))
  df <- df %>% dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
  if (nrow(df) < (degree + 1L) * (degree + 2L) / 2L) return(NULL) # not enough points
  
  # build formula like: z ~ poly(x,deg,raw=TRUE) + poly(y,deg,raw=TRUE) + cross terms
  # easiest in base R via manual term expansion:
  terms <- c()
  for (i in 0:degree) {
    for (j in 0:(degree - i)) {
      if (i == 0 && j == 0) next
      ti <- if (i == 0) "1" else paste0("I(x^", i, ")")
      tj <- if (j == 0) "1" else paste0("I(y^", j, ")")
      t  <- if (i == 0) tj else if (j == 0) ti else paste0(ti, ":", tj)
      terms <- c(terms, t)
    }
  }
  fml <- stats::as.formula(paste("z ~", paste(terms, collapse = " + ")))
  
  fit <- stats::lm(fml, data = df)
  
  xr <- range(df$x, na.rm = TRUE)
  yr <- range(df$y, na.rm = TRUE)
  
  gx <- seq(xr[1], xr[2], length.out = grid_n)
  gy <- seq(yr[1], yr[2], length.out = grid_n)
  grid <- expand.grid(x = gx, y = gy)
  grid$z <- as.numeric(stats::predict(fit, newdata = grid))
  
  list(
    fit = fit,
    x = matrix(grid$x, nrow = grid_n, ncol = grid_n, byrow = FALSE),
    y = matrix(grid$y, nrow = grid_n, ncol = grid_n, byrow = FALSE),
    z = matrix(grid$z, nrow = grid_n, ncol = grid_n, byrow = FALSE)
  )
}

add_cluster2_lines_and_angles <- function(fig, df, plane,
                                          use_pca_normal = TRUE,
                                          cluster2_table_store = NULL,
                                          helix_turns = 3,
                                          show_plane_angles = FALSE,
                                          show_plane_distance = FALSE,
                                          show_pca_axes = FALSE,
                                          enable_clustering = FALSE,
                                          dark_mode = FALSE,
                                          bending_radius = FALSE,
                                          show_info_labels = FALSE,
                                          compute_rmsd_metrics = FALSE,
                                          selected_fit_type = "sphere",
                                          selected_fit_method = "heuristic",
                                          fit_weighted = FALSE,
                                          mle_center_mode = "free",
                                          mle_center_lambda = 0.1,
                                          arc_span_hint = 0) {
  
  df <- df %>% dplyr::filter(!is.na(cluster2) & cluster2 != 0)
  if (nrow(df) == 0) return(list(fig = fig, log = "No Median Cluster points found."))
  
  log_msgs <- character(0)
  result_table <- data.frame()
  clusters <- df %>% dplyr::group_by(cluster2) %>% dplyr::group_split()
  
  for (cluster_data in clusters) {

  cluster_id <- unique(cluster_data$cluster2)[1]
  
  fmt <- function(x, digits = 2) {
    ifelse(is.na(x) || !is.finite(x), "NA",
           format(round(x, digits), nsmall = digits, trim = TRUE))
  }

  # --- representative points for THIS cluster2 (used for ALL metrics/fit geometry) ---
  has_valid_clusters <- "cluster" %in% names(cluster_data) &&
    any(!is.na(cluster_data$cluster), na.rm = TRUE)

  if (isTRUE(enable_clustering) && has_valid_clusters) {
    medians_df <- cluster_data %>%
      dplyr::filter(!is.na(cluster) & cluster != 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        x        = stats::median(x, na.rm = TRUE),
        y        = stats::median(y, na.rm = TRUE),
        z        = stats::median(z, na.rm = TRUE),
        spread   = spread_from_group(x, y, z),
        n_locs   = dplyr::n(),
        group_id = dplyr::first(cluster),
        .groups  = "drop"
      )
  } else {
    medians_df <- cluster_data %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        x        = stats::median(x, na.rm = TRUE),
        y        = stats::median(y, na.rm = TRUE),
        z        = stats::median(z, na.rm = TRUE),
        spread   = spread_from_group(x, y, z),
        n_locs   = dplyr::n(),
        group_id = dplyr::first(tid),
        .groups  = "drop"
      )
  }

  medians_df <- medians_df %>% dplyr::filter(is.finite(x), is.finite(y), is.finite(z))

  if (nrow(medians_df) < 4) {
    log_msgs <- c(log_msgs, paste("Median Cluster", cluster_id, "- Not enough representative medians (<4)."))
    next
  }
  need_rmsd <- isTRUE(compute_rmsd_metrics)
  
  # ---- representative point matrix for this cluster2 ----
  coords_matrix <- as.matrix(medians_df[, c("x","y","z"), drop = FALSE])
  
  # Weights for MLE fits: proportional to n_locs / spread² (≈ 1/SEM² of the
  # median position).  Traces with a single localization have spread = 0 by
  # construction (variance is undefined for n = 1); they receive the median
  # spread of multi-localization traces as a conservative fallback so they
  # participate at roughly average weight rather than near-infinite weight.
  # RMSD metrics are always reported unweighted.
  w <- if (isTRUE(fit_weighted) && "spread" %in% names(medians_df)) {
    s   <- medians_df$spread
    n_l <- if ("n_locs" %in% names(medians_df)) medians_df$n_locs else rep(1L, nrow(medians_df))
    s_ref <- stats::median(s[n_l > 1 & s > 0], na.rm = TRUE)
    if (!is.finite(s_ref) || s_ref <= 0) s_ref <- max(s[s > 0], na.rm = TRUE)
    if (!is.finite(s_ref) || s_ref <= 0) s_ref <- 1
    s_eff <- ifelse(n_l > 1 & s > 0, s, s_ref)
    n_l / s_eff^2
  } else NULL
  
  # defaults (fast path)
  rmsd_h_sphere <- meandist_h_sphere <- sd_h_sphere <- NA_real_
  rmsd_h_ellipsoid <- meandist_h_ellipsoid <- sd_h_ellipsoid <- NA_real_
  rmsd_h_helix <- meandist_h_helix <- sd_h_helix <- NA_real_
  best_fit_name <- NA_character_

  rmsd_m_sphere <- meandist_m_sphere <- sd_m_sphere <- NA_real_
  rmsd_m_ellipsoid <- meandist_m_ellipsoid <- sd_m_ellipsoid <- NA_real_
  rmsd_m_helix <- meandist_m_helix <- sd_m_helix <- NA_real_
  rmsd_h_ring  <- meandist_h_ring  <- sd_h_ring  <- NA_real_
  rmsd_h_arc   <- meandist_h_arc   <- sd_h_arc   <- NA_real_
  rmsd_m_ring  <- meandist_m_ring  <- sd_m_ring  <- NA_real_
  rmsd_m_arc   <- meandist_m_arc   <- sd_m_arc   <- NA_real_

  # ── Active model RMSD: always computed (regardless of compute_rmsd_metrics) ─
  active_rmsd <- active_meandist <- active_sd <- NA_real_
  if (!is.null(selected_fit_type) && selected_fit_type != "none") {
    .active_ft <- tryCatch(
      calculate_distance_to_fit(
        coords_matrix = coords_matrix,
        fit_type      = selected_fit_type,
        turns         = helix_turns,
        arc_span_hint = arc_span_hint
      ),
      error = function(e) NULL
    )
    if (!is.null(.active_ft)) {
      active_rmsd     <- attr(.active_ft, paste0("rmsd_",     selected_fit_type))
      active_meandist <- attr(.active_ft, paste0("meandist_", selected_fit_type))
      active_sd       <- attr(.active_ft, paste0("sd_",       selected_fit_type))
    }
    # MLE path: override with the MLE projection for the selected model only
    if (identical(selected_fit_method, "mle")) {
      center_args_act <- build_center_args(
        X = coords_matrix, mode = mle_center_mode, lambda = mle_center_lambda
      )
      mle_act <- tryCatch(
        switch(selected_fit_type,
          sphere    = do.call(mle_fit_sphere,
                              c(list(X = coords_matrix, w = w), center_args_act)),
          ellipsoid = do.call(mle_fit_ellipsoid,
                              c(list(X = coords_matrix, w = w), center_args_act)),
          helix     = do.call(mle_fit_helix,
                              c(list(X = coords_matrix, w = w, turns_hint = helix_turns), center_args_act)),
          ring      = do.call(mle_fit_ring,
                              c(list(X = coords_matrix, w = w), center_args_act)),
          arc       = do.call(mle_fit_arc,
                              c(list(X = coords_matrix, w = w, arc_span_hint = arc_span_hint), center_args_act)),
          NULL
        ),
        error = function(e) NULL
      )
      if (!is.null(mle_act)) {
        d_act <- switch(selected_fit_type,
          sphere    = as.numeric(mle_project_sphere(coords_matrix, mle_act$center, mle_act$radius)$dist),
          ellipsoid = as.numeric(mle_project_ellipsoid(coords_matrix, mle_act$center, mle_act$axes, mle_act$R)$dist),
          helix = {
            if (!is.null(mle_act$t)) {
              t_a <- .project_t_helix(coords_matrix, mle_act$center, mle_act$axis,
                                       mle_act$radius, mle_act$pitch, mle_act$phi0,
                                       t0 = mle_act$t, check_neighbors = TRUE)
              H_a <- t(vapply(t_a, .helix_point, numeric(3),
                               c = mle_act$center, u = mle_act$axis,
                               r = mle_act$radius, p = mle_act$pitch, phi0 = mle_act$phi0))
              sqrt(rowSums((coords_matrix - H_a)^2))
            } else numeric(0)
          },
          ring = as.numeric(mle_project_ring(coords_matrix, mle_act$center, mle_act$normal, mle_act$radius)$dist),
          arc  = as.numeric(mle_project_arc(coords_matrix, mle_act$center, mle_act$normal,
                                             mle_act$radius, mle_act$t_center, mle_act$arc_half)$dist),
          numeric(0)
        )
        d_act <- d_act[is.finite(d_act)]
        if (length(d_act)) {
          active_rmsd     <- sqrt(mean(d_act^2))
          active_meandist <- mean(d_act)
          active_sd       <- if (length(d_act) > 1) stats::sd(d_act) else NA_real_
        }
      }
    }
  }

  if (need_rmsd) {
    # ---- Heuristic RMSD/SD (passes pre-computed coords_matrix; no re-medianing) ----
    fit_result <- calculate_distance_to_fit(
      coords_matrix = coords_matrix,
      fit_type      = "best",
      turns         = helix_turns
    )

    if (!is.null(fit_result)) {
      rmsd_h_sphere    <- attr(fit_result, "rmsd_sphere");        sd_h_sphere    <- attr(fit_result, "sd_sphere")
      meandist_h_sphere    <- attr(fit_result, "meandist_sphere")
      rmsd_h_ellipsoid <- attr(fit_result, "rmsd_ellipsoid");    sd_h_ellipsoid <- attr(fit_result, "sd_ellipsoid")
      meandist_h_ellipsoid <- attr(fit_result, "meandist_ellipsoid")
      rmsd_h_helix     <- attr(fit_result, "rmsd_helix");        sd_h_helix     <- attr(fit_result, "sd_helix")
      meandist_h_helix     <- attr(fit_result, "meandist_helix")
      rmsd_h_ring      <- attr(fit_result, "rmsd_ring");         sd_h_ring      <- attr(fit_result, "sd_ring")
      meandist_h_ring      <- attr(fit_result, "meandist_ring")
      rmsd_h_arc       <- attr(fit_result, "rmsd_arc");          sd_h_arc       <- attr(fit_result, "sd_arc")
      meandist_h_arc       <- attr(fit_result, "meandist_arc")
    }
    
    # ---- MLE RMSD/SD (only when MLE is the selected method) ----
    if (identical(selected_fit_method, "mle")) {
      center_args <- build_center_args(
        X = coords_matrix, mode = mle_center_mode, lambda = mle_center_lambda
      )

      fS <- tryCatch(do.call(mle_fit_sphere, c(list(X = coords_matrix, w = w), center_args)), error = function(e) NULL)
      if (!is.null(fS)) {
        pr <- mle_project_sphere(coords_matrix, fS$center, fS$radius)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_sphere     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_sphere <- if (length(d)) mean(d)         else NA_real_
        sd_m_sphere       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fE <- tryCatch(do.call(mle_fit_ellipsoid, c(list(X = coords_matrix, w = w), center_args)), error = function(e) NULL)
      if (!is.null(fE)) {
        pr <- mle_project_ellipsoid(coords_matrix, fE$center, fE$axes, fE$R)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_ellipsoid     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_ellipsoid <- if (length(d)) mean(d)         else NA_real_
        sd_m_ellipsoid       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fH <- tryCatch(do.call(mle_fit_helix, c(list(X = coords_matrix, w = w, turns_hint = helix_turns), center_args)), error = function(e) NULL)
      if (!is.null(fH) && !is.null(fH$t)) {
        t_rmsd <- .project_t_helix(coords_matrix, fH$center, fH$axis, fH$radius, fH$pitch, fH$phi0,
                                   t0 = fH$t, check_neighbors = TRUE)
        H_rmsd <- t(vapply(t_rmsd, .helix_point, numeric(3),
                           c = fH$center, u = fH$axis, r = fH$radius, p = fH$pitch, phi0 = fH$phi0))
        d <- sqrt(rowSums((coords_matrix - H_rmsd)^2)); d <- d[is.finite(d)]
        rmsd_m_helix     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_helix <- if (length(d)) mean(d)         else NA_real_
        sd_m_helix       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fR <- tryCatch(do.call(mle_fit_ring, c(list(X = coords_matrix, w = w), center_args)), error = function(e) NULL)
      if (!is.null(fR)) {
        pr <- mle_project_ring(coords_matrix, fR$center, fR$normal, fR$radius)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_ring     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_ring <- if (length(d)) mean(d)         else NA_real_
        sd_m_ring       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }

      fA <- tryCatch(do.call(mle_fit_arc, c(list(X = coords_matrix, w = w, arc_span_hint = arc_span_hint), center_args)), error = function(e) NULL)
      if (!is.null(fA)) {
        pr <- mle_project_arc(coords_matrix, fA$center, fA$normal, fA$radius, fA$t_center, fA$arc_half)
        d  <- as.numeric(pr$dist); d <- d[is.finite(d)]
        rmsd_m_arc     <- if (length(d)) sqrt(mean(d^2)) else NA_real_
        meandist_m_arc <- if (length(d)) mean(d)         else NA_real_
        sd_m_arc       <- if (length(d) > 1) stats::sd(d) else NA_real_
      }
    }
  }
    
  # ---- best fit across all 6 models (heuristic + MLE), by RMSD -----------
  if (need_rmsd) {
    all_rmsd <- c(
      sphere_h    = rmsd_h_sphere,    ellipsoid_h    = rmsd_h_ellipsoid,    helix_h    = rmsd_h_helix,
      ring_h      = rmsd_h_ring,      arc_h          = rmsd_h_arc,
      sphere_mle  = rmsd_m_sphere,    ellipsoid_mle  = rmsd_m_ellipsoid,    helix_mle  = rmsd_m_helix,
      ring_mle    = rmsd_m_ring,      arc_mle        = rmsd_m_arc
    )
    finite_rmsd <- all_rmsd[is.finite(all_rmsd)]
    best_fit_name <- if (length(finite_rmsd)) names(which.min(finite_rmsd)) else NA_character_
  }

  # ---- convex hull (compute ONCE) -> center + area + volume ----
  hull <- NULL
  surface_area  <- NA_real_
  convex_volume <- NA_real_
  
  if (nrow(coords_matrix) >= 4) {
    hull <- tryCatch(
      geometry::convhulln(coords_matrix, options = "FA"),
      error = function(e) {
        log_msgs <<- c(log_msgs, paste("Median Cluster", cluster_id, "- Convex hull error:", e$message))
        NULL
      }
    )
  }
  
  if (!is.null(hull)) {
    center <- compute_convex_hull_centroid(coords_matrix)
    surface_area  <- hull$area
    convex_volume <- hull$vol
  } else {
    center <- colMeans(coords_matrix)
  }
    
    distances <- sqrt(rowSums((coords_matrix - matrix(center, nrow = nrow(coords_matrix), ncol = 3, byrow = TRUE))^2))
    radius <- median(distances)
    # ONE curvature convention app-wide: MEAN curvature H = 1/R (nm^-1), the
    # quantity membrane-bending work is expressed in. This used to be 1/R^2 here
    # (Gaussian curvature K for a sphere) while estimate_cluster2_curvature()
    # reported 1/R, so the same word meant two different things in one codebase.
    # For a fitted sphere both are just reparametrisations of `radius` and carry
    # no information beyond it.
    mean_curvature <- if (radius > 0) 1 / radius else NA
    sphere_area <- if (!is.na(radius) && radius > 0) 4 * pi * radius^2 else NA
    sphere_volume <- if (!is.na(radius) && radius > 0) (4/3) * pi * radius^3 else NA
    
    
    pca <- prcomp(coords_matrix, scale. = FALSE)
    rotated <- (coords_matrix - matrix(center, nrow = nrow(coords_matrix), ncol = 3, byrow = TRUE)) %*% pca$rotation
    d <- diff(stats::quantile(rotated[, 1], c(0.05, 0.95))) / 2
    e <- diff(stats::quantile(rotated[, 2], c(0.05, 0.95))) / 2
    f <- diff(stats::quantile(rotated[, 3], c(0.05, 0.95))) / 2
    scale_factor <- stats::quantile(sqrt(rowSums(rotated^2)), 0.90) / max(c(d, e, f))
    d <- d * scale_factor; e <- e * scale_factor; f <- f * scale_factor
    
    if (min(d, e, f) < 1e-6) {
      ellipsoid_area <- NA; ellipsoid_volume <- NA
    } else {
      # Knud Thomsen approximation (g = 1.6075); exact only for a sphere.
      # Error < ~1.2% for moderate axis ratios; increases for highly elongated shapes.
      g <- 1.6075
      ellipsoid_area <- 4 * pi * ((d^g * e^g + d^g * f^g + e^g * f^g) / 3)^(1/g)
      ellipsoid_volume <- (4/3) * pi * d * e * f
    }
    
    ellipticity <- round((max(d,e,f) - min(d,e,f))/max(d,e,f), 2)
    axes <- sort(c(d, e, f), decreasing = TRUE); a <- axes[1]; b <- axes[2]; c <- axes[3]
    # Eccentricities of the three principal cross-sections (2-D ellipse formula).
    # ecc1: major (PC1) vs intermediate (PC2) plane.
    # ecc2: major (PC1) vs minor (PC3) plane.
    # ecc3: intermediate (PC2) vs minor (PC3) plane.
    # All three are needed to fully characterise a general 3-D ellipsoid.
    #
    # ecc3 == 0 IS A RESULT, NOT A FAILURE: it is zero exactly when the two
    # SMALLER axes are equal, i.e. the structure is rotationally symmetric about
    # its long axis (rod / cigar / spherocylinder). Verified against built
    # shapes: prolate 400:150:150 -> ecc1 = ecc2 = 0.930, ecc3 = 0.057;
    # triaxial 400:250:100 -> ecc3 = 0.915 against a truth of 0.917.
    #
    # CAVEAT (see docs I-A.6b): PCA orders axes by SAMPLE variance, so the
    # largest sampled direction always becomes PC1 - the ordering manufactures
    # apparent elongation. On a TRUE SPHERE the median reported ecc2 is 0.81 at
    # 10 members, 0.53 at 40, 0.34 at 200, 0.18 at 2000. These numbers are only
    # evidence of elongation when they clearly exceed that floor for the
    # cluster's own member count.
    ecc1 <- if (is.finite(d) && d > 1e-6) sqrt(max(0, 1 - (e/d)^2)) else NA_real_
    ecc2 <- if (is.finite(d) && d > 1e-6) sqrt(max(0, 1 - (f/d)^2)) else NA_real_
    ecc3 <- if (is.finite(e) && e > 1e-6) sqrt(max(0, 1 - (f/e)^2)) else NA_real_
    # Bending (perceived) radius - ASSUMPTION: the cluster represents a curved arc
    # where PC1 spans the chord half-length and PC2 is the sagitta (bow height).
    # Formula: R = (chord²/2 + sagitta²) / (2 * sagitta).
    # Valid for ring segments (e.g. S. aureus division ring viewed on edge) but
    # not meaningful for spheres, arbitrary ellipsoids, or helices.
    perceived_radius_1 <- ifelse(e > 0, (d^2 + e^2) / (2 * e), NA)
    perceived_radius_2 <- ifelse(f > 0, (d^2 + f^2) / (2 * f), NA)
    perceived_curv_1 <- 1/(perceived_radius_1^2)
    perceived_curv_2 <- 1/(perceived_radius_2^2)
    perceived_area_1 <- if (!is.na(perceived_radius_1) && perceived_radius_1 > 0) 4 * pi * perceived_radius_1^2 else NA
    prog_area <- surface_area/perceived_area_1
    
    coords <- coords_matrix
    tmp <- safe_nn2(coords, k = 4)
    if (is.null(tmp$nn)) {
      nn_idx  <- matrix(integer(0), nrow = 0, ncol = 0)
      nn_dist <- matrix(numeric(0),  nrow = 0, ncol = 0)
    } else {
      nn_res <- tmp$nn
      usedK  <- max(0, tmp$k_eff - 1)
      if (usedK == 0) {
        nn_idx  <- matrix(integer(0), nrow = nrow(coords_matrix), ncol = 0)
        nn_dist <- matrix(numeric(0),  nrow = nrow(coords_matrix), ncol = 0)
      } else {
        nn_idx  <- nn_res$nn.idx[, 2:(usedK + 1), drop = FALSE]
        nn_dist <- nn_res$nn.dists[, 2:(usedK + 1), drop = FALSE]
      }
    }
    
    nn_df <- data.frame()
    for (i in seq_len(nrow(coords_matrix))) {
      origin_id <- medians_df$group_id[i]
      spread_i <- medians_df$spread[i]
      for (j in seq_len(min(3, ncol(nn_idx)))) {
        nbr_idx <- nn_idx[i, j]
        if (is.na(nbr_idx) || nbr_idx < 1) next
        nbr_id <- medians_df$group_id[nbr_idx]
        spread_j <- medians_df$spread[nbr_idx]
        distance <- nn_dist[i, j]
        spread_cutoff <- max(1.5 * max(spread_i, spread_j, na.rm = TRUE), 6)
        if (is.finite(spread_cutoff) && distance > spread_cutoff) {
          nn_df <- rbind(nn_df, data.frame(
            from = origin_id, to = nbr_id, distance = distance,
            spread_i = spread_i, spread_j = spread_j, spread_cutoff = spread_cutoff
          ))
        }
      }
    }
    nn_stats <- list(
      nn_n = nrow(nn_df),
      nn_mean = if (nrow(nn_df)) mean(nn_df$distance) else NA_real_,
      nn_sd   = if (nrow(nn_df))  sd(nn_df$distance)  else NA_real_,
      nn_min  = if (nrow(nn_df)) min(nn_df$distance) else NA_real_,
      nn_max  = if (nrow(nn_df)) max(nn_df$distance) else NA_real_,
      nn_med  = if (nrow(nn_df)) median(nn_df$distance) else NA_real_
    )
    
    
    symmetry_plane_normal <- pca$rotation[,3]
    symmetry_plane_normal <- symmetry_plane_normal / sqrt(sum(symmetry_plane_normal^2))
    vectors_to_center <- t(t(coords) - center)
    distances_to_plane <- vectors_to_center %*% symmetry_plane_normal
    reflected_coords <- coords - 2 * distances_to_plane %*% t(symmetry_plane_normal)
    nn_reflect <- RANN::nn2(data = coords, query = reflected_coords, k = 1)
    mean_nn_distance <- mean(nn_reflect$nn.dists)
    bounding_radius <- max(distances)
    semi_score <- mean_nn_distance / bounding_radius
    
    if (!is.null(plane) && is.list(plane) && !is.null(plane$normal) && !is.null(plane$center)) {
      n_vec <- plane$normal / sqrt(sum(plane$normal^2))
      p0 <- plane$center
      p <- center
      vector_p0p <- p - p0
      distance_to_plane <- sum(n_vec * vector_p0p)
      foot_point <- p - distance_to_plane * n_vec
    } else {
      distance_to_plane <- NA
      foot_point <- NA
    }
    
    fit_center <- center
    if (!is.null(plane) && use_pca_normal && all(c("normal", "center") %in% names(plane))) {
      base_vector <- plane$normal / sqrt(sum(plane$normal^2))
      orange_vector <- if (distance_to_plane < 0) -base_vector else base_vector
    } else {
      orange_vector <- c(0, 0, 1)
    }
    
    angle_deg <- NA
    end_point_orange <- rep(NA_real_, 3)
    end_point_blue   <- rep(NA_real_, 3)
    if (!is.na(distance_to_plane) && abs(distance_to_plane) <= radius) {
      n_vec <- plane$normal / sqrt(sum(plane$normal^2))
      v_perp <- -distance_to_plane * n_vec
      l_parallel <- sqrt(radius^2 - distance_to_plane^2)
      z_up <- c(0, 0, 1)
      if (abs(sum(n_vec * z_up)) > 0.999) z_up <- c(0, 1, 0)
      in_plane_dir <- pracma::cross(n_vec, z_up)
      in_plane_dir <- in_plane_dir / sqrt(sum(in_plane_dir^2))
      blue_vector <- v_perp + l_parallel * in_plane_dir
      blue_vector <- blue_vector / sqrt(sum(blue_vector^2))
      if (sign(sum(blue_vector * n_vec)) == sign(distance_to_plane)) blue_vector <- -blue_vector
      end_point_orange <- center + orange_vector * radius
      end_point_blue   <- center + radius * blue_vector
      dot_product <- sum(orange_vector * blue_vector)
      angle_rad <- acos(min(max(dot_product, -1), 1))
      angle_deg <- round(angle_rad * 180 / pi, 2)
      label_offset <- 0.1 * radius
    }
    
    if (show_plane_angles) {
      fig <- fig %>% plotly::add_trace(
        x = c(center[1], end_point_orange[1]),
        y = c(center[2], end_point_orange[2]),
        z = c(center[3], end_point_orange[3]),
        type = "scatter3d", mode = "lines",
        line = list(color = "orange", width = 4),
        inherit = FALSE, showlegend = FALSE
      ) %>% plotly::add_trace(
        x = c(center[1], end_point_blue[1]),
        y = c(center[2], end_point_blue[2]),
        z = c(center[3], end_point_blue[3]),
        type = "scatter3d", mode = "lines",
        line = list(color = "blue", width = 4),
        inherit = FALSE, showlegend = FALSE
      ) %>% plotly::add_trace(
        x = c(center[1] + label_offset),
        y = c(center[2] + label_offset),
        z = c(center[3] + label_offset),
        type = "scatter3d", mode = "text",
        text = paste0(angle_deg, "°"),
        textfont = list(color = if (dark_mode) "lightblue" else "darkblue", size = 11),
        inherit = FALSE, showlegend = FALSE
      )
    }
    
    if (isTRUE(bending_radius) && !is.na(perceived_radius_1)) {
      end_point_orange <- center + orange_vector * radius
      bend_dir <- center - end_point_orange
      bend_dir <- bend_dir / sqrt(sum(bend_dir^2))
      end_point_perceived <- end_point_orange + bend_dir * perceived_radius_1
      fig <- fig %>% plotly::add_trace(
        x = c(end_point_orange[1], end_point_perceived[1]),
        y = c(end_point_orange[2], end_point_perceived[2]),
        z = c(end_point_orange[3], end_point_perceived[3]),
        type  = "scatter3d", mode  = "lines",
        line  = list(color = "black", width = 4, dash = "dash"),
        name  = "Perceived radius",
        inherit = FALSE, showlegend = TRUE
      )
    }
    
    if (show_plane_distance && !any(is.na(foot_point))) {
      fig <- fig %>% plotly::add_trace(
        x = c(center[1], foot_point[1]),
        y = c(center[2], foot_point[2]),
        z = c(center[3], foot_point[3]),
        type = "scatter3d", mode = "lines",
        line = list(color = if (dark_mode) "#F6CEFC" else "purple", width = 6, dash = 'dot'),
        inherit = FALSE, showlegend = FALSE
      ) %>% plotly::add_trace(
        x = (center[1] + foot_point[1]) / 2,
        y = (center[2] + foot_point[2]) / 2,
        z = (center[3] + foot_point[3]) / 2,
        type = "scatter3d", mode = "text",
        text = paste0(round(abs(distance_to_plane), 1)),
        textfont = list(color = if (dark_mode) "#F6CEFC" else "purple", size = 10),
        inherit = FALSE, showlegend = FALSE
      )
    }
    
    if (show_pca_axes) {
      axis_colors <- c("red", "green", "blue")
      axis_labels <- c("PC1", "PC2", "PC3")
      for (i in 1:3) {
        axis_vector <- pca$rotation[, i] * c(d, e, f)[i]
        fig <- fig %>% plotly::add_trace(
          x = c(center[1], center[1] + axis_vector[1]),
          y = c(center[2], center[2] + axis_vector[2]),
          z = c(center[3], center[3] + axis_vector[3]),
          type = "scatter3d", mode = "lines+text",
          line = list(color = axis_colors[i], width = 6, dash = "dash"),
          text = axis_labels[i], textposition = "top right",
          textfont = list(color = if (dark_mode) "white" else "black", size = 10),
          showlegend = FALSE, inherit = FALSE
        )
      }
    }
    
    if (isTRUE(show_info_labels)) {
      
      
      hover_txt <- paste0(
        "Median cluster: ", cluster_id,
        "<br>CH Surface Area=", fmt(surface_area, 1)
      )
      
      if (identical(selected_fit_type, "sphere")) {
        hover_txt <- paste0(hover_txt,
                            "<br>Radius=", fmt(radius, 1),
                            "<br>Sphere Surface Area=", fmt(sphere_area, 1))
      } else if (identical(selected_fit_type, "ellipsoid")) {
        hover_txt <- paste0(hover_txt,
                            "<br>Ellipticity=", fmt(ellipticity, 2),
                            "<br>Ecc1(PC1-PC2)=", fmt(ecc1, 2),
                            "<br>Ecc2(PC1-PC3)=", fmt(ecc2, 2),
                            "<br>Ecc3(PC2-PC3)=", fmt(ecc3, 2),
                            "<br>Ellipsoid Surface Area (approx)=", fmt(ellipsoid_area, 1))
      } else if (identical(selected_fit_type, "helix")) {
        hover_txt <- paste0(hover_txt, "<br>Helix turns=", helix_turns)
      }
      
      if (!is.na(active_rmsd)) {
        hover_txt <- paste0(
          hover_txt,
          "<br>RMSD(", selected_fit_method, "/", selected_fit_type, ")=",
          fmt(active_rmsd, 2),
          ", mean=", fmt(active_meandist, 2),
          " ± ", fmt(active_sd, 2)
        )
      }
      
      fig <- fig %>% plotly::add_trace(
        type = "scatter3d",
        mode = "markers",
        x = center[1], y = center[2], z = center[3],
        marker = list(size = 6, opacity = 0.001),  # basically invisible but hoverable
        hovertemplate = paste0(hover_txt, "<extra></extra>"),
        hoverlabel = list(font = list(size = 13)),
        showlegend = FALSE,
        inherit = FALSE,
        name = paste0("Median cluster ", cluster_id)
      )
    }

    result_table <- rbind(result_table, data.frame(
      cluster2 = cluster_id,
      center_x = center[1], center_y = center[2], center_z = center[3],
      radius = radius, mean_curvature_inv_nm = mean_curvature,
      surface_area_convex_hull = surface_area,
      surface_area_sphere = sphere_area,
      surface_area_ellipsoid = ellipsoid_area,
      ellipticity = ellipticity,
      eccentricity_1 = ecc1, eccentricity_2 = ecc2, eccentricity_3 = ecc3,
      volume_convex_hull = convex_volume,
      volume_sphere = sphere_volume,
      volume_ellipsoid = ellipsoid_volume,
      best_fit = best_fit_name,
      rmsd_sphere_h = rmsd_h_sphere, meandist_sphere_h = meandist_h_sphere, sd_sphere_h = sd_h_sphere,
      rmsd_ellipsoid_h = rmsd_h_ellipsoid, meandist_ellipsoid_h = meandist_h_ellipsoid, sd_ellipsoid_h = sd_h_ellipsoid,
      rmsd_helix_h = rmsd_h_helix, meandist_helix_h = meandist_h_helix, sd_helix_h = sd_h_helix,
      rmsd_sphere_mle = rmsd_m_sphere, meandist_sphere_mle = meandist_m_sphere, sd_sphere_mle = sd_m_sphere,
      rmsd_ellipsoid_mle = rmsd_m_ellipsoid, meandist_ellipsoid_mle = meandist_m_ellipsoid, sd_ellipsoid_mle = sd_m_ellipsoid,
      rmsd_helix_mle = rmsd_m_helix, meandist_helix_mle = meandist_m_helix, sd_helix_mle = sd_m_helix,
      rmsd_ring_h = rmsd_h_ring, meandist_ring_h = meandist_h_ring, sd_ring_h = sd_h_ring,
      rmsd_arc_h  = rmsd_h_arc,  meandist_arc_h  = meandist_h_arc,  sd_arc_h  = sd_h_arc,
      rmsd_ring_mle = rmsd_m_ring, meandist_ring_mle = meandist_m_ring, sd_ring_mle = sd_m_ring,
      rmsd_arc_mle  = rmsd_m_arc,  meandist_arc_mle  = meandist_m_arc,  sd_arc_mle  = sd_m_arc,
      pc1_length = 2*d, pc2_length = 2*e, pc3_length = 2*f,
      pc2_pc1 = (2*e)/(2*d), pc3_pc1 = (2*f)/(2*d),
      nn_n = nn_stats$nn_n, nn_mean = nn_stats$nn_mean, nn_sd = nn_stats$nn_sd,
      nn_min = nn_stats$nn_min, nn_max = nn_stats$nn_max, nn_med = nn_stats$nn_med,
      bending_radius_1 = perceived_radius_1, bending_radius_2 = perceived_radius_2,
      perceived_curv_1 = perceived_curv_1, perceived_curv_2 = perceived_curv_2,
      perceived_area_1 = perceived_area_1, prog_area = prog_area,
      distance_to_plane = distance_to_plane, closing_angle_deg = angle_deg
    ))
    
    msg <- paste0(
      "Median cluster ", cluster_id,
      ": CH Area=", fmt(surface_area, 2),
      ", CH Vol=",  fmt(convex_volume, 2)
    )
    
    shape_msg <- switch(
      as.character(selected_fit_type),
      "sphere" = paste0(
        ", Radius=", fmt(radius, 2),
        ", Sphere Area=", fmt(sphere_area, 2),
        ", Sphere Vol=", fmt(sphere_volume, 2),
        ", Mean curv=", ifelse(is.na(mean_curvature), "NA", paste0(signif(mean_curvature, 3), " nm^-1"))
      ),
      "ellipsoid" = paste0(
        ", Ellipticity=", fmt(ellipticity, 2),
        ", Ecc1(PC1-PC2)=", fmt(ecc1, 2),
        ", Ecc2(PC1-PC3)=", fmt(ecc2, 2),
        ", Ecc3(PC2-PC3)=", fmt(ecc3, 2),
        ", Ellipsoid Area(approx)=", fmt(ellipsoid_area, 2),
        ", Ellipsoid Vol=", fmt(ellipsoid_volume, 2)
      ),
      "helix" = paste0(", Helix turns=", helix_turns),
      ""
    )
    
    msg <- paste0(msg, shape_msg)
    
    if (!is.na(active_rmsd)) {
      msg <- paste0(
        msg,
        ", RMSD(", selected_fit_method, "/", selected_fit_type, ")=",
        fmt(active_rmsd, 2),
        ", mean=", fmt(active_meandist, 2),
        " ± ", fmt(active_sd, 2)
      )
    }

    log_msgs <- c(log_msgs, msg)
  }
  if (!is.null(cluster2_table_store)) cluster2_table_store(result_table)
  list(fig = fig, log = paste(log_msgs, collapse = "\n"), table = result_table)
}

compute_cluster2_plane_distances <- function(df) {
  df_medians <- df %>% group_by(tid) %>% summarize(x = median(x), y = median(y), z = median(z), .groups = "drop")
  pca <- prcomp(df_medians[, c("x", "y", "z")], scale. = FALSE)
  normal <- pca$rotation[, 3]
  center <- colMeans(df_medians[, c("x", "y", "z")])
  
  df_cluster2 <- df %>%
    filter(!is.na(cluster2) & cluster2 != 0) %>%
    group_by(cluster2) %>%
    summarize(
      x = median(x), y = median(y), z = median(z),
      trace_ids = paste(unique(tid), collapse = ","),
      .groups = "drop"
    )
  
  df_cluster2$distance_to_plane <- apply(df_cluster2[, c("x", "y", "z")], 1, function(row) {
    sum(normal * (row - center)) / sqrt(sum(normal^2))
  })
  
  df_cluster2
}

generate_xy_plane <- function(df, limits, z_offset = 0) {
  medians <- df %>%
    group_by(tid) %>%
    summarize(x = median(x), y = median(y), z = median(z), .groups = "drop")
  
  pca <- prcomp(medians[, c("x", "y", "z")], scale. = FALSE)
  center <- colMeans(medians[, c("x", "y", "z")])
  center[3] <- center[3] + z_offset  # Apply offset
  
  normal <- pca$rotation[, 3]
  grid_size <- max(diff(limits$x), diff(limits$y)) / 2
  grid_x <- seq(center[1] - grid_size, center[1] + grid_size, length.out = 20)
  grid_y <- seq(center[2] - grid_size, center[2] + grid_size, length.out = 20)
  grid <- expand.grid(x = grid_x, y = grid_y)
  
  a <- normal[1]
  b <- normal[2]
  c <- normal[3]
  d <- sum(normal * center)
  grid$z <- (d - a * grid$x - b * grid$y) / c
  
  list(
    x = matrix(grid$x, nrow = length(grid_y)),
    y = matrix(grid$y, nrow = length(grid_y)),
    z = matrix(grid$z, nrow = length(grid_y)),
    normal = normal,  
    center = center          
  )
}

#Helper: surface plot
create_global_surface <- function(df, grid_n = 120, span = 0.25) {
  grid_x <- seq(min(df$x, na.rm = TRUE), max(df$x, na.rm = TRUE), length.out = grid_n)
  grid_y <- seq(min(df$y, na.rm = TRUE), max(df$y, na.rm = TRUE), length.out = grid_n)
  grid   <- expand.grid(x = grid_x, y = grid_y)
  
  fit    <- loess(z ~ x + y, data = df, span = span, control = loess.control(surface = "direct"))
  z_pred <- predict(fit, newdata = grid)
  
  stopifnot(length(z_pred) == nrow(grid))
  list(
    x = matrix(grid$x, nrow = length(grid_y), byrow = TRUE),
    y = matrix(grid$y, nrow = length(grid_y), byrow = TRUE),
    z = matrix(z_pred,  nrow = length(grid_y), byrow = TRUE)
  )
}

# =============================================================================
# cluster2_member_col() - WHICH UNIT IS A "MEMBER" OF A MEDIAN CLUSTER
# =============================================================================
#
# GOAL
#   One rule, shared by every per-structure statistic: a member is a DE-BLINKED
#   MOLECULE when de-blinking has run, otherwise a TRACE ID.
#
# WHY IT MATTERS
#   A trace ID is one blinking event. One fluorophore commonly emits several, so
#   counting trace IDs counts blinks, not molecules, and a structure's measured
#   composition then moves with photophysics rather than with biology. This
#   mirrors the server-side smlm_nn_unit() used by the Crowding and
#   NN-vs-size panels, so all per-structure numbers agree on the unit.
cluster2_member_col <- function(df, member_col = NULL) {
  if (!is.null(member_col)) return(member_col)
  if ("cluster" %in% names(df) && any(!is.na(df$cluster) & df$cluster > 0L)) "cluster" else "tid"
}

# Bandwidth of the isosurface KDE, as a multiple of the member gap. See the
# `sigma` block inside cluster_density_field() for how this was chosen: it is a
# compromise between a surface that fragments (too low) and one that has smoothed
# the structure away (too high), and the mesh smoothing downstream means less
# bandwidth is needed here than when the voxel grid was rendered directly.
SIGMA_GAP_FACTOR <- 0.85

# =============================================================================
# cluster_density_field() - GAUSSIAN KDE OVER A VOXEL GRID + A DENSITY LEVEL
#                           THAT ENCLOSES A STATED FRACTION OF THE LABELLING
# =============================================================================
#
# GOAL
#   Give the 3-D viewer a smooth envelope to wrap a Median Cluster in, built
#   from a field that actually responds to how dense the labelling is.
#
# WHAT IT REPLACED, AND WHY THAT MATTERED
#   The isosurface used to score each voxel by its distance to the SINGLE
#   NEAREST member: val = exp(-d^2 / 2 sigma^2). That is a distance transform,
#   not a density - measured, a probe 60 nm from a lone isolated centroid and a
#   probe 60 nm from a centroid with 30 neighbours within 200 nm both scored
#   0.9231, identical. Visually that means every straggler grows its own
#   detached bubble of exactly the same size as the real structure's rim, which
#   is both misleading and ugly. Summing a kernel over ALL members instead makes
#   the field proportional to local labelling density, so the body of the
#   structure dominates and isolated members fade out.
#
# THE LEVEL IS CALIBRATED TO A SURFACE OFFSET IN NANOMETRES
#   `offset` = "how far outside the member centroids should the skin sit". The
#   level is read off the field at that distance: voxels whose nearest-member
#   distance falls in a band around `offset` vote, and their MEDIAN value is the
#   contour. The surface then follows the offset while density still modulates
#   it - sparse patches pull in slightly, dense ones push out slightly, which is
#   what makes it look like a skin rather than a union of balls.
#
# THE OFFSET IS PHYSICAL, THE BANDWIDTH IS COSMETIC
#   `offset` (AUTO = 3 x localization precision, floored at half the member gap)
#   sets HOW TIGHT the skin is, in nanometres that mean something on the
#   microscope. `sigma` (AUTO = SIGMA_GAP_FACTOR x gap) sets HOW SMOOTH it is.
#   Both are applied PER AXIS, in units of that axis's own precision. They are
#   independent by construction: the level is re-measured at `offset` using the
#   current sigma, so smoothing more does not inflate the structure.
#
#   This is the property that makes the surface comparable between clusters. The
#   earlier rule keyed both dials to the member spacing, so a poorly labelled
#   copy of a structure drew larger than a well labelled one - on an R = 150 nm
#   shell the realised radius ran from 304 nm at 30 members down to 192 nm at
#   400, all of it labelling artifact. Anchored to precision the same four cases
#   land at 166 / 163 / 173 / 169 nm.
#
#   `sampling_limited` in the returned list says which regime you are in: TRUE
#   means the gap floor set the offset, so the surface is as tight as this
#   labelling density permits and the remaining bulk is missing molecules rather
#   than real structure.
#
#   WHY NOT AN ENCLOSED-DENSITY FRACTION (what this replaced). A Gaussian
#   spreads its mass over ~3 sigma in every direction, so "enclose 50% of the
#   labelling density" is already a thick band: measured on an R=150 shell with
#   40 members, frac = 0.50 - the lowest the slider allowed - put the skin
#   56.4 nm out (max 128.5), and frac = 0.85 put it 108.7 nm out. A mass
#   fraction cannot express "hug the points"; a distance can.
#
#   WHY NOT A FRACTION OF THE PEAK either: it collapses too fast. On the same
#   shell, peak fraction 0.4 -> 681 voxels, 0.6 -> 122, 0.8 -> 16, i.e. the
#   surface shatters into fragments before it gets tight.
#
# THE DIALS
#   `offset` sets the wrap and `sigma` the smoothing, both defaulting as above;
#   `pad` adds to or subtracts from the offset and is the coarser tighter/looser
#   control. Grid extent stays derived, so no knob can silently fight another.
#
# BOX
#   `grid_res` is a COST BUDGET of grid_res^3 voxels, not a per-axis count. The
#   voxels are cubic and the box is sized to members + offset + 2 voxels, which
#   is just enough for marching cubes to close the surface without clipping it.
#   Spending the budget per-axis instead gave anisotropic voxels on elongated
#   structures (28 nm along a rod, 10 nm across it, so it broke up lengthwise)
#   and the old offset + 2*sigma margin left 32-70% of the box empty. The box is
#   NOT a user control: when padding was a separate box setting it both overrode
#   a request for a tight box (forced to 3*sigma) and leaked into the contour.
#
# NOT A CALIBRATED DENSITY ESTIMATE
#   The kernel is unnormalised and sampled on a coarse voxel grid, so `value` is
#   not a density in localizations/nm^3. It is a visualization envelope; use the
#   Ripley/PCF panels for density statistics.
cluster_density_field <- function(pts, grid_res = 22L, pad = 0, sigma = NULL,
                                  offset = NULL, precision = NULL,
                                  close_surface = TRUE, sigma_scale = 1) {
  if (is.null(pts)) return(NULL)
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  if (ncol(pts) < 3L) return(NULL)
  pts <- pts[stats::complete.cases(pts), 1:3, drop = FALSE]
  storage.mode(pts) <- "double"
  n <- nrow(pts)
  if (n < 4L) return(NULL)

  # ── work in PRECISION UNITS, not nanometres ──
  # `precision` may be one number or c(x, y, z). Dividing each axis by its own
  # precision makes the space isotropic in the only currency that matters here:
  # how well a molecule's position is known along that axis. Everything below is
  # then a plain isotropic problem, and the mesh is scaled back to nm on the way
  # out. The skin that comes back is an ellipsoidal shell in real space, tighter
  # on whichever axis is measured better.
  #
  # For MINFLUX that axis is usually z (measured z:lateral 0.60 on a real imaging file),
  # so an isotropic wrap inflates the structure axially by about two thirds. Note
  # the direction: this is the opposite of astigmatic camera SMLM, and assuming
  # the camera convention here would make the error worse rather than better.
  prec_ax <- suppressWarnings(as.numeric(precision %||% NA_real_))
  if (length(prec_ax) == 1L) prec_ax <- rep(prec_ax, 3L)
  if (length(prec_ax) != 3L) prec_ax <- rep(prec_ax[1], 3L)
  prec_ax[!is.finite(prec_ax) | prec_ax <= 0] <- NA_real_
  if (all(is.na(prec_ax))) prec_ax <- rep(2, 3L)
  prec_ax[is.na(prec_ax)] <- mean(prec_ax, na.rm = TRUE)

  # nm <-> precision-unit conversion for the dials the user types in nm. The
  # LATERAL precision is the reference: it is the axis people quote, and it keeps
  # a manual "30 nm wrap" meaning 30 nm laterally.
  prec_ref <- mean(prec_ax[1:2])
  pts_nm   <- pts
  pts      <- sweep(pts, 2, prec_ax, "/")

  # ── two length scales, and they measure different things ──
  # `spacing`  = the TYPICAL member separation (median NN). Reported for
  #              context; it is deliberately NOT what the geometry keys off.
  # `gap`      = the LARGE member separations (90th percentile NN). A closed
  #              surface has to bridge the big gaps, not the typical one, so
  #              this is the scale that decides whether the skin closes.
  #
  # WHY THE DISTINCTION IS NOT PEDANTRY. For a quasi-regular point set the two
  # coincide, but real labelling is Poisson-like, and there the median NN sits
  # far below the gaps that actually need bridging: measured on random points on
  # a 600 x 120 nm rod, median NN was 19.8 nm while the 90th percentile was
  # 39.7 nm, a factor of 2.0. That factor is exactly what the old
  # `offset = 2 x median spacing` rule was silently paying for. Naming the gap
  # directly lets the offset shrink to something physical without the surface
  # shattering.
  nnd <- tryCatch({
    d <- RANN::nn2(pts, pts, k = 2L)$nn.dists[, 2]
    d[is.finite(d) & d > 0]
  }, error = function(e) numeric(0))
  spacing <- if (length(nnd)) stats::median(nnd) else NA_real_
  gap     <- if (length(nnd)) as.numeric(stats::quantile(nnd, 0.90, na.rm = TRUE)) else NA_real_
  if (!is.finite(spacing) || spacing <= 0 || !is.finite(gap) || gap <= 0) {
    diag_nm <- sqrt(sum((apply(pts, 2, max) - apply(pts, 2, min))^2))
    spacing <- if (is.finite(diag_nm) && diag_nm > 0) diag_nm / 10 else 20
    gap     <- spacing
  }

  # In precision units one unit IS one localization sigma, on every axis.
  precision <- 1

  # ── surface offset: how far outside the members the skin sits ──
  # AUTO = 3 x localization precision, i.e. a 3-sigma skin around where each
  # molecule actually is, floored at half the gap.
  #
  # WHY 3 x PRECISION AND NOT 2 x SPACING (what this replaced). The old rule read
  # its scale off the labelling, so the SAME structure drew bigger when it was
  # sampled worse. Measured on an R = 150 nm shell, the realised outer radius was
  # 304 nm at 30 members (+103%), 263 nm at 60 (+76%), 216 nm at 150 (+44%) and
  # 192 nm at 400 (+28%): the drawn size tracked labelling density, and a sparse
  # shell even rendered hollow where a dense one filled in solid. Anchoring to
  # precision instead put the same four cases at 166 / 163 / 173 / 169 nm, which
  # is the truth (150) plus the shell's own thickness plus the offset.
  #
  # WHY THE GAP FLOOR IS NOT OPTIONAL. A continuous surface passing outside a set
  # of points cannot come nearer than half their separation to all of them, so an
  # offset below gap/2 cannot close whatever the bandwidth is. The floor states
  # that limit instead of letting the surface shatter into one bead per molecule.
  # On well-labelled data the precision term wins and the skin is genuinely tight;
  # on sparse data the floor wins and the offset reports the real sampling limit.
  offset <- suppressWarnings(as.numeric(offset %||% NA_real_)) / prec_ref  # nm in
  if (!is.finite(offset) || offset <= 0) offset <- max(3 * precision, 0.5 * gap)

  # `pad` adds to (or, negative, subtracts from) that - the coarser/tighter dial.
  pad <- suppressWarnings(as.numeric(pad %||% 0)) / prec_ref                # nm in
  if (!is.finite(pad)) pad <- 0
  offset <- max(offset + pad, 0.5 * precision)     # floor: never collapse to nothing

  # ── bandwidth: the SMOOTHNESS dial, largely independent of the offset ──
  # AUTO = SIGMA_GAP_FACTOR x gap. Because the contour level is re-measured at
  # `offset` with the current sigma, raising sigma does not inflate the surface
  # the way it did when the two were tied together. Measured on a real imaging file, the
  # MEDIAN vertex-to-molecule distance moved 25.6 -> 29.0 nm against a 25.3 nm
  # target across a 5x sigma sweep.
  #
  # It is not a perfect invariance and should not be sold as one: the same sweep
  # took the MAXIMUM distance from 39.8 to 65.7 nm, because a wider kernel
  # absorbs small protrusions into the body rather than following them. Sigma
  # holds the typical wrap and trades away fine relief.
  #
  # Measured at fixed offset, sweeping sigma over 0.75 to 2.0 x gap on shells, a
  # solid ball and rods: components fell to 1 by ~1.25x, shells stayed HOLLOW up
  # to 1.5x and collapsed into a solid ball at 2.0x. The factor now sits BELOW
  # that closure point on purpose, because two things downstream did not exist
  # when it was measured: the level search closes the surface without spending
  # bandwidth on it (see `close_surface`), and the mesh is smoothed as geometry
  # rather than as density. Carrying enough bandwidth to close the surface on its
  # own is what made real structures look oversmoothed - it blurs the measurement
  # to fix what is really a topology problem. Smoothing the mesh instead leaves
  # the density where it was measured.
  #
  # It was previously tied to the offset (max(0.5 * offset, spacing)), so
  # tightening the wrap also sharpened the field and re-fragmented the surface -
  # the two dials fought each other.
  # `sigma_scale` tunes the AUTOMATIC bandwidth only; an explicit `sigma` in nm
  # is a full override and is left alone. Keeping them separate means the
  # relative dial stays meaningful no matter what the data does to the auto value.
  sigma_scale <- suppressWarnings(as.numeric(sigma_scale %||% 1))
  if (!is.finite(sigma_scale) || sigma_scale <= 0) sigma_scale <- 1
  sigma <- suppressWarnings(as.numeric(sigma %||% NA_real_)) / prec_ref      # nm in
  if (!is.finite(sigma) || sigma <= 0)
    sigma <- max(SIGMA_GAP_FACTOR * sigma_scale * gap, precision)

  # ── box: cubic voxels, spent where the surface actually is ──
  # `grid_res` is a COST budget (grid_res^3 voxels), not a per-axis count. Two
  # things were wrong with using it per-axis:
  #
  #   1. margin = offset + 2 * sigma padded the box far past anything the surface
  #      could reach - the skin stops at ~offset from the outermost member - so
  #      32-70% of the voxels held empty space and coarsened the rest. Measured
  #      box occupancy on 300 nm structures: 30% at 30 members, 68% at 400.
  #   2. an equal count per axis gives ANISOTROPIC voxels on anything elongated.
  #      A 600 x 120 nm rod got 28.2 nm voxels along x and 10 nm across it, so it
  #      broke up lengthwise while being oversampled across.
  #
  # Cubic voxels under a fixed budget fix both at identical cost: the same rod
  # now gets 11.7 nm voxels on a 59x20x20 grid, and a 300 nm shell goes from 25.6
  # to 18.3 nm. Two voxels of headroom past the offset keeps marching cubes from
  # clipping the surface against the wall of the box.
  grid_res <- max(8L, as.integer(round(grid_res)))
  budget   <- as.numeric(grid_res)^3

  rng <- apply(pts, 2, range)
  ext <- pmax(rng[2, ] - rng[1, ], 1e-9)           # guard a flat/degenerate axis

  # margin depends on the voxel size and vice versa; one fixed-point pass from
  # margin = offset is enough, the correction is a few percent
  vox_for <- function(budget) {                    # solve vox against its own margin
    v <- (prod(ext + 2 * offset) / budget)^(1 / 3)
    for (.i in 1:2) v <- (prod(ext + 2 * (offset + 2 * v)) / budget)^(1 / 3)
    v
  }
  vox <- vox_for(budget)

  # ── resolution floor: the voxel must be able to RESOLVE the skin ──
  # Marching cubes needs the enclosed shell to be a couple of voxels thick or it
  # drops pieces of it, and the shell is only ~`offset` thick. Now that the
  # offset is precision-anchored it is often 10-13 nm, well under the ~19 nm
  # voxel that grid_res = 22 buys - measured, that alone broke a closed shell
  # into 3 fragments purely as a sampling artifact, with nothing wrong with the
  # field. Refine until the skin is ~1.5 voxels thick.
  #
  # This raises the user's grid_res, so it is capped: VOXEL_CAP bounds the redraw
  # at roughly two seconds for the largest clusters. Beyond that the surface is
  # under-resolved and `resolved` comes back FALSE rather than the app quietly
  # spending ten seconds on one cluster.
  VOXEL_CAP <- 2e5
  vox_max   <- offset / 1.5
  resolved  <- TRUE
  if (vox > vox_max) {
    need <- prod(ext + 2 * (offset + 2 * vox_max)) / (vox_max^3)
    if (need <= VOXEL_CAP) {
      vox <- vox_max
    } else {
      vox      <- vox_for(VOXEL_CAP)
      resolved <- vox <= vox_max
    }
  }
  margin <- offset + 2 * vox

  n_ax <- pmax(6L, as.integer(ceiling((ext + 2 * margin) / vox)) + 1L)
  xs <- seq(rng[1, 1] - margin, rng[2, 1] + margin, length.out = n_ax[1])
  ys <- seq(rng[1, 2] - margin, rng[2, 2] + margin, length.out = n_ax[2])
  zs <- seq(rng[1, 3] - margin, rng[2, 3] + margin, length.out = n_ax[3])
  grid <- as.matrix(expand.grid(x = xs, y = ys, z = zs))

  # ── truncated KDE: a Gaussian is numerically dead past 3 sigma, so only the
  #    nearest members can contribute.
  #
  # k IS MEASURED, NOT ASSUMED. It used to be min(n, 512), which sets the cost of
  # every grid point by the SIZE OF THE CLUSTER rather than by how many members
  # can actually reach a probe. The nn2 result is a (grid points x k) matrix, so
  # at k = 512 a 280k-voxel grid allocates 145M distances (~1.2 GB) to multiply
  # almost all of them by a weight of zero: measured 22.8 s for one 1500-member
  # cluster where the real neighbour count within 3 sigma was a few dozen.
  #
  # Instead, count how many members genuinely fall within 3 sigma of a member and
  # size k from that, with headroom. The probe subset is a deterministic stride,
  # never sample() - drawing at random here would advance the GLOBAL RNG stream
  # and silently perturb every downstream random draw (the permutation nulls, the
  # alpha-surface subsample), and a redraw could jitter the surface.
  k_probe <- as.integer(round(seq(1, n, length.out = min(n, 256L))))
  k_probe <- unique(k_probe)
  k_need <- tryCatch({
    kn <- RANN::nn2(pts, pts[k_probe, , drop = FALSE], k = min(n, 512L))
    max(rowSums(kn$nn.dists <= 3 * sigma))
  }, error = function(e) min(n, 512L))
  # 1.5x headroom covers members denser than anywhere in the probe subset; the
  # floor keeps small clusters exact. Truncation past this is not an
  # approximation being traded away - those weights are already numerically 0.
  kmax <- min(n, max(8L, as.integer(ceiling(1.5 * k_need)) + 4L))
  nn <- tryCatch(RANN::nn2(pts, grid, k = kmax), error = function(e) NULL)
  if (is.null(nn)) return(NULL)
  w <- exp(-(nn$nn.dists^2) / (2 * sigma * sigma))
  w[nn$nn.dists > 3 * sigma] <- 0
  val <- as.numeric(rowSums(w))
  if (!any(is.finite(val)) || max(val, na.rm = TRUE) <= 0) return(NULL)

  # ── level calibrated to the offset, measured OFF-GRID ──
  # Probe points are scattered on spheres of radius `offset` around randomly
  # chosen members; the median field value over the probes that really do sit
  # ~offset from the point set is the contour.
  #
  # This is deliberately independent of the voxel grid. Reading the level from
  # grid voxels instead lets the BOX size leak into the surface: one far-flung
  # member stretches the bounding box, the voxels coarsen, the "within one
  # voxel of the offset" band widens to tens of nm, and the median is taken
  # over voxels far from any member - which dropped the level low enough to
  # re-draw the isolated straggler that the density field exists to suppress.
  # Directions come from a Fibonacci sphere and members are cycled, so the probe
  # set is deterministic (a redraw cannot jitter the surface) and evenly spread.
  # Drawing them at random would also have meant calling set.seed() here, which
  # would reset the GLOBAL RNG stream and silently perturb every downstream
  # random draw - the permutation nulls, the alpha-surface subsample.
  n_probe <- 1024L
  kk    <- seq_len(n_probe) - 0.5
  phi   <- acos(1 - 2 * kk / n_probe)
  theta <- pi * (1 + sqrt(5)) * kk
  dir   <- cbind(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi))
  mi    <- rep_len(seq_len(n), n_probe)
  probes <- pts[mi, , drop = FALSE] + dir * offset

  pn <- tryCatch(RANN::nn2(pts, probes, k = kmax), error = function(e) NULL)
  if (is.null(pn)) return(NULL)
  # keep probes that are genuinely at the offset - the rest landed back inside
  # the structure, where the field says nothing about where the skin belongs
  keep_p <- pn$nn.dists[, 1] >= 0.8 * offset
  if (!any(keep_p)) keep_p <- rep(TRUE, nrow(probes))
  pw <- exp(-(pn$nn.dists^2) / (2 * sigma * sigma))
  pw[pn$nn.dists > 3 * sigma] <- 0
  level <- stats::median(rowSums(pw)[keep_p])

  # plotly draws nothing if the level sits at or above the peak
  vmax <- max(val, na.rm = TRUE)
  if (!is.finite(level) || level <= 0) level <- vmax * 1e-3
  if (level >= vmax) level <- vmax * 0.999

  # ── solve the level down until the surface is ONE body ──
  # The offset fixes where the skin belongs; this asks the far cheaper question
  # of whether that level actually yields a connected surface, and if not, how
  # far it has to give. Lowering the level grows the surface, which is the exact
  # inflation the offset rewrite removed, so the search is BOUNDED and the cost
  # is reported rather than absorbed: `closed` says whether it succeeded and
  # `level_relaxed` how many times the level had to step down.
  #
  # The field is computed once; only the threshold moves. Connectivity is tested
  # on the voxel mask by label propagation, which is a few vectorised array
  # shifts over a grid this size.
  lvl_components <- function(thr) {
    mask <- array(val >= thr, n_ax)
    if (!any(mask)) return(.Machine$integer.max)
    lab <- array(0L, n_ax); idx <- which(mask); lab[idx] <- idx
    repeat {
      prev <- lab
      for (ax in 1:3) for (sh in c(-1L, 1L)) {
        sl <- switch(ax,
          lab[if (sh == 1L) c(1L, seq_len(n_ax[1] - 1L)) else c(2:n_ax[1], n_ax[1]), , , drop = FALSE],
          lab[, if (sh == 1L) c(1L, seq_len(n_ax[2] - 1L)) else c(2:n_ax[2], n_ax[2]), , drop = FALSE],
          lab[, , if (sh == 1L) c(1L, seq_len(n_ax[3] - 1L)) else c(2:n_ax[3], n_ax[3]), drop = FALSE])
        hit <- mask & sl > 0L & (lab == 0L | sl < lab)
        lab[hit] <- sl[hit]
      }
      lab[!mask] <- 0L
      if (identical(lab, prev)) break
    }
    length(unique(lab[mask]))
  }

  # NA, not TRUE: when the search is off nothing has checked connectivity, and
  # reporting an untested surface as closed is worse than reporting nothing
  closed <- NA; level_relaxed <- 0L
  if (isTRUE(close_surface)) {
    ncomp <- lvl_components(level)
    # 0.85 per step, at most 6 steps: the level can fall to ~0.38 of the
    # offset-calibrated value and no further. Past that the surface has stopped
    # describing the offset at all and a fragmented honest answer beats a
    # connected invented one.
    while (ncomp > 1L && level_relaxed < 6L) {
      level <- level * 0.85
      level_relaxed <- level_relaxed + 1L
      ncomp <- lvl_components(level)
    }
    closed <- ncomp == 1L
  }

  # ── back to nanometres ──
  # Everything above ran in precision units. `grid` is scaled per axis (so the
  # caller gets real coordinates); the scalar diagnostics are reported against
  # the lateral precision, which is the axis they are quoted in.
  grid_nm <- sweep(grid, 2, prec_ax, "*")

  list(grid = grid_nm, value = val, level = level,
       sigma = sigma * prec_ref, offset = offset * prec_ref,
       spacing = spacing * prec_ref, gap = gap * prec_ref,
       precision = prec_ax, precision_ref = prec_ref,
       margin = margin * prec_ref, voxel = vox * prec_ref,
       # axis scale factors, so a caller can rebuild the anisotropic geometry
       axis_scale = prec_ax, dims = n_ax, n_members = n,
       # grid/level in the isotropic working space, for the mesh builder
       grid_units = grid, breaks_units = list(x = xs, y = ys, z = zs),
       # FALSE when the voxel cap stopped the grid short of resolving the skin;
       # the surface is then under-sampled and may show gaps that are not real
       resolved = resolved,
       # did the level search reach a single connected body, and at what cost
       closed = closed, level_relaxed = level_relaxed,
       # TRUE when the sampling, not the microscope, set the wrap distance -
       # i.e. the surface is as tight as this labelling density allows, no tighter
       sampling_limited = (0.5 * gap) > (3 * precision))
}

# =============================================================================
# density_field_mesh() - MARCHING CUBES OVER A DENSITY FIELD, WELDED + SMOOTHED
# =============================================================================
#
# GOAL
#   Turn the voxel field from cluster_density_field() into a genuine triangle
#   mesh that renders as a smooth skin, instead of handing the raw voxels to
#   plotly's `isosurface` trace.
#
# WHY NOT plotly's OWN isosurface TRACE (what this replaced)
#   plotly runs marching cubes per cell and shades the result FLAT, so every
#   voxel face is visible as its own facet. That is the blocky look: at the
#   default grid the voxels were tens of nm across, so the structure was drawn
#   with tens-of-nm flat plates regardless of how good the underlying data was.
#   It also had to be driven with a hairline isomin/isomax band around the level,
#   which renders a thin SLAB with an inner and an outer wall - under a
#   translucent opacity those two walls shade separately and read as bulk.
#
# THE THREE STEPS THAT MAKE IT LOOK LIKE A SURFACE
#   1. MARCHING CUBES (misc3d::computeContour3d) on the field. This is the same
#      algorithm plotly uses; the point is what happens next.
#   2. WELD. computeContour3d emits independent triangles, so a shared corner
#      appears once per triangle touching it - measured 83% duplicate vertex rows
#      on a test sphere. plotly can only average normals across triangles that
#      literally share a vertex INDEX, so without welding `flatshading = FALSE`
#      has nothing to work with and the mesh shades flat no matter what. Welding
#      is what buys smooth shading; it also cuts the payload about six-fold.
#   3. TAUBIN SMOOTH. Marching cubes leaves a staircase at the voxel scale.
#      Taubin's lambda/mu pass (a positive smoothing step followed by a slightly
#      larger negative one) removes it WITHOUT the shrinkage that plain Laplacian
#      smoothing causes - a Laplacian-only mesh creeps inward every iteration and
#      would quietly undo the offset calibration the field just established.
#
# WHY SMOOTH THE MESH RATHER THAN THE DENSITY
#   Both remove the staircase, but they are not the same operation. Raising the
#   KDE bandwidth blurs the MEASUREMENT: it merges neighbouring molecules and
#   pulls real detail out of the structure. Smoothing the mesh only regularises
#   the ISOSURFACE EXTRACTION, which is a rendering artifact to begin with. That
#   split is why SIGMA_GAP_FACTOR could be lowered from 1.25 to 0.85 here.
#
# NOT A VOLUME ESTIMATOR
#   The mesh is for looking at. Use cluster_alpha_volume() for enclosed volume;
#   the smoothing step alone moves this surface by a few percent.
#
# RETURNS
#   list(x, y, z, i, j, k) with 0-based i/j/k ready for plotly::add_trace(type =
#   "mesh3d"), plus n_vertices / n_triangles, or NULL when the level yields no
#   surface at all.
density_field_mesh <- function(fld, smooth_iter = 8L, lambda = 0.5, mu = -0.53,
                               weld_tol = 1e-6) {
  if (is.null(fld) || is.null(fld$value) || is.null(fld$dims)) return(NULL)
  if (!requireNamespace("misc3d", quietly = TRUE)) return(NULL)

  n_ax <- fld$dims
  vol  <- array(fld$value, n_ax)
  br   <- fld$breaks_units          # isotropic working space (precision units)
  if (is.null(br)) return(NULL)

  # ── 1. marching cubes ──
  tri <- tryCatch(
    misc3d::computeContour3d(vol, level = fld$level, x = br$x, y = br$y, z = br$z),
    error = function(e) NULL
  )
  if (is.null(tri) || !is.matrix(tri) || nrow(tri) < 3L) return(NULL)

  # ── 2. weld ──
  # Quantise to a tolerance relative to the voxel, so vertices that marching
  # cubes computed independently but that ARE the same point collapse together.
  # match() on the packed key is the vectorised way to do this; paste() of three
  # rounded numbers is stable because they came from the same interpolation.
  vox_u <- if (length(br$x) > 1L) abs(br$x[2] - br$x[1]) else 1
  q     <- round(tri / (vox_u * weld_tol))
  key   <- paste(q[, 1], q[, 2], q[, 3], sep = ",")
  first <- !duplicated(key)
  verts <- tri[first, , drop = FALSE]
  idx   <- match(key, key[first])            # every row -> its welded vertex

  faces <- matrix(idx, ncol = 3L, byrow = TRUE)
  # marching cubes can emit slivers where two corners welded together; a triangle
  # with a repeated vertex has no area and no normal, and would poison the
  # averaged normal of every vertex it touches
  keep_f <- faces[, 1] != faces[, 2] & faces[, 2] != faces[, 3] & faces[, 1] != faces[, 3]
  faces  <- faces[keep_f, , drop = FALSE]
  if (!nrow(faces)) return(NULL)

  # ── 3. Taubin smoothing ──
  if (smooth_iter > 0L && nrow(verts) > 3L) {
    # Undirected edge list, each edge once per direction so the neighbour average
    # is symmetric. Built once and reused for every iteration.
    e_from <- c(faces[, 1], faces[, 2], faces[, 3], faces[, 2], faces[, 3], faces[, 1])
    e_to   <- c(faces[, 2], faces[, 3], faces[, 1], faces[, 1], faces[, 2], faces[, 3])
    ord    <- !duplicated(paste(e_from, e_to, sep = ","))
    e_from <- e_from[ord]; e_to <- e_to[ord]
    deg    <- tabulate(e_from, nbins = nrow(verts))
    deg[deg == 0L] <- 1L                     # isolated vertex: leave it alone

    # Neighbour-averaging operator as a sparse matrix, built once. Each
    # iteration is then one sparse matvec instead of three tapply() calls over
    # tens of thousands of edges, which is what made this the slow step.
    A <- Matrix::sparseMatrix(i = e_from, j = e_to, x = 1 / deg[e_from],
                              dims = c(nrow(verts), nrow(verts)))
    for (it in seq_len(smooth_iter)) {
      verts <- verts + lambda * (as.matrix(A %*% verts) - verts)   # shrink
      verts <- verts + mu     * (as.matrix(A %*% verts) - verts)   # push back out
    }
  }

  # ── back to nanometres: undo the per-axis precision scaling ──
  sc <- fld$axis_scale %||% c(1, 1, 1)
  verts <- sweep(verts, 2, sc, "*")

  list(x = verts[, 1], y = verts[, 2], z = verts[, 3],
       i = faces[, 1] - 1L, j = faces[, 2] - 1L, k = faces[, 3] - 1L,   # plotly is 0-based
       n_vertices = nrow(verts), n_triangles = nrow(faces))
}

# GOAL: one hex colour for a set of members, or NA when they carry none.
# Used by the isosurface so a surface is coloured by the SAME channel_color the
# localizations inside it use. Returns NA (not a guess) when the members span
# more than one channel colour, because a two-channel structure has no single
# honest colour and the caller should fall back to its own palette.
iso_channel_color <- function(members) {
  if (!is.data.frame(members) || !nrow(members)) return(NA_character_)
  if (!"channel_color" %in% names(members)) return(NA_character_)
  cc <- unique(members$channel_color[!is.na(members$channel_color)])
  if (length(cc) == 1L) as.character(cc) else NA_character_
}

# =============================================================================
# add_cluster_isosurface() - ONE ISOSURFACE TRACE, USED BY EVERY VIEWING MODE
# =============================================================================
#
# GOAL
#   Wrap one structure in one smooth skin and add it to a plotly 3-D figure.
#   Normal mode, the 2-colour split and multiplex channels all call this; they
#   differ only in which points they hand over and which colour they ask for.
#
# WHY IT EXISTS
#   The cluster2 and tc2 render blocks were near-duplicates that had already
#   drifted apart: only one of them had been switched to de-blinked molecule
#   centroids, and only one honoured channel colours. A wrapping rule that lives
#   in two places is a wrapping rule that will disagree with itself.
#
# COLOUR IS THE CALLER'S DECISION
#   This takes a single hex colour and does not know what it means. The caller
#   picks it from the cluster palette, the channel table or a fixed user colour,
#   which is what keeps "colour by channel" from having to be re-implemented per
#   viewing mode.
#
# RENDERING NOTES
#   The trace is a `mesh3d`, not plotly's `isosurface`. flatshading = FALSE plus
#   the welded vertices from density_field_mesh() is what produces smooth shading;
#   either alone gives flat facets. Lighting is deliberately low-specular: a
#   translucent skin with a hot highlight reads as a plastic bead rather than as
#   a density envelope.
#
# RETURNS
#   The figure, with one trace added. On any failure (too few members, no
#   surface at that level) the figure comes back UNCHANGED rather than erroring,
#   because one bad cluster must not blank the whole plot.
add_cluster_isosurface <- function(fig, pts, color = "#00CED1", opacity = 0.25,
                                   name = NULL, showlegend = FALSE,
                                   grid_res = 22L, pad = 0, sigma = NULL,
                                   offset = NULL, precision = NULL,
                                   smooth_iter = 8L, close_surface = TRUE,
                                   sigma_scale = 1, diag_env = NULL, cache = NULL) {
  # `diag_env` collects one row per surface so the UI can report what the
  # automatic geometry actually chose. An environment rather than a return value
  # because the caller wants the figure back in a pipe; nothing here reads it.
  cache_row <- NULL
  note <- function(...) {
    row <- data.frame(..., stringsAsFactors = FALSE)
    cache_row <<- row
    if (is.null(diag_env)) return(invisible(NULL))
    diag_env$rows <- c(diag_env$rows, list(row))
  }
  if (is.null(pts)) return(fig)
  pts <- as.matrix(pts)
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]
  if (nrow(pts) < 4L || ncol(pts) < 3L) {
    note(name = name %||% NA_character_, n_members = nrow(pts), wrap = NA_real_,
         smoothing = NA_real_, voxel = NA_real_, sampling_limited = NA,
         closed = NA, resolved = NA, drawn = FALSE)
    return(fig)
  }

  # ── geometry cache ──
  # Colour and opacity do not enter the key, because they do not change the mesh.
  # A Shiny render re-runs on ANY of its reactive dependencies, so without this a
  # change of opacity or colour mode rebuilt every field and every mesh from
  # scratch. The caller clears the cache when the committed geometry changes.
  ck <- NULL
  if (!is.null(cache)) {
    ck <- digest::digest(list(pts, grid_res, pad, sigma, offset, precision,
                              smooth_iter, close_surface, sigma_scale))
    hit <- cache[[ck]]
    if (!is.null(hit)) {
      if (!is.null(diag_env) && !is.null(hit$diag))
        diag_env$rows <- c(diag_env$rows, list(hit$diag))
      if (is.null(hit$msh)) return(fig)
      return(iso_mesh_trace(fig, hit$msh, color, opacity, name, showlegend))
    }
  }

  fld <- tryCatch(
    cluster_density_field(pts, grid_res = grid_res, pad = pad, sigma = sigma,
                          offset = offset, precision = precision,
                          close_surface = close_surface, sigma_scale = sigma_scale),
    error = function(e) NULL
  )
  if (is.null(fld)) {
    note(name = name %||% NA_character_, n_members = nrow(pts), wrap = NA_real_,
         smoothing = NA_real_, voxel = NA_real_, sampling_limited = NA,
         closed = NA, resolved = NA, drawn = FALSE)
    return(fig)
  }

  msh <- tryCatch(density_field_mesh(fld, smooth_iter = smooth_iter),
                  error = function(e) NULL)
  note(name = name %||% NA_character_, n_members = nrow(pts),
       wrap = fld$offset, smoothing = fld$sigma, voxel = fld$voxel,
       sampling_limited = isTRUE(fld$sampling_limited),
       closed = fld$closed, resolved = isTRUE(fld$resolved),
       drawn = !is.null(msh))
  if (!is.null(ck)) cache[[ck]] <- list(msh = msh, diag = cache_row)
  if (is.null(msh)) return(fig)

  iso_mesh_trace(fig, msh, color, opacity, name, showlegend)
}

# GOAL: add one finished mesh to a figure. Split out so a cache hit and a fresh
# build draw through exactly the same code.
#
# COLOUR MUST NOT BE PASSED AS `color`. plotly's R wrapper intercepts a `color`
# argument and routes it through its own scale-mapping machinery, which has no
# output slot for mesh3d - so the attribute is silently DROPPED and every surface
# renders in plotly's default. Verified against plotly_build(): both
# `color = "#FF0000"` and `color = I("#FF0000")` vanish from the built trace,
# while intensity/colorscale, facecolor and vertexcolor all survive.
#
# A flat intensity against a two-stop colorscale is the lightest of the three:
# one number per vertex rather than a hex STRING per vertex (vertexcolor) or per
# triangle (facecolor), which matters when a view draws tens of surfaces.
# cmin/cmax are pinned because an all-equal intensity would otherwise autoscale
# to a degenerate cmin == cmax range.
iso_mesh_trace <- function(fig, msh, color, opacity, name, showlegend) {
  nv <- length(msh$x)
  fig %>% plotly::add_trace(
    type = "mesh3d",
    x = msh$x, y = msh$y, z = msh$z,
    i = msh$i, j = msh$j, k = msh$k,
    intensity = rep(0, nv),
    colorscale = list(list(0, color), list(1, color)),
    cmin = 0, cmax = 1, showscale = FALSE,
    opacity = opacity,
    flatshading = FALSE,
    lighting = list(ambient = 0.55, diffuse = 0.8, specular = 0.12,
                    roughness = 0.9, fresnel = 0.1),
    lightposition = list(x = 1e4, y = 1e4, z = 1e4),
    name = name %||% "Isosurface",
    showlegend = isTRUE(showlegend),
    hoverinfo = "skip",
    inherit = FALSE
  )
}

# =============================================================================
# iso_status_ui() - WHAT THE AUTOMATIC ISOSURFACE GEOMETRY ACTUALLY CHOSE
# =============================================================================
#
# GOAL
#   Report, under the isosurface controls, the numbers the surface was built
#   with and anything that stopped it being drawn properly.
#
# WHY THIS IS NOT OPTIONAL POLISH
#   The wrap distance, the smoothing and the grid are all derived from the data,
#   so "Tighten / loosen" adjusts a value the user cannot otherwise see. Without
#   this readout, -20 nm is a blind edit. It also separates the two ways a
#   surface can look wrong that would otherwise be indistinguishable on screen:
#   a structure that is genuinely fragmented, and one that is merely
#   under-resolved by the voxel grid.
#
# THE TWO REGIMES
#   The wrap is 3 x localization precision unless that falls below half the
#   molecule spacing, in which case spacing wins. Which one bound it is the
#   single most useful thing to know about a surface: a spacing-bound surface is
#   as tight as this labelling density permits and its remaining bulk is missing
#   molecules, not structure.
#
# RETURNS
#   A shiny tag, or NULL when nothing has been drawn yet.
iso_status_ui <- function(rows) {
  if (is.null(rows) || !length(rows)) return(NULL)
  d <- if (is.data.frame(rows)) rows else do.call(rbind, rows)
  if (is.null(d) || !nrow(d)) return(NULL)

  ok  <- d[isTRUE_vec(d$drawn), , drop = FALSE]
  n_f <- sum(!isTRUE_vec(d$drawn))
  if (!nrow(ok)) {
    return(shiny::div(class = "iso-status",
      shiny::tags$span(style = "color:#b45309;",
        sprintf("No surface drawn (%d cluster%s too small or too sparse).",
                n_f, if (n_f == 1L) "" else "s"))))
  }

  rng <- function(v, digits = 0) {
    v <- v[is.finite(v)]
    if (!length(v)) return("n/a")
    if (diff(range(v)) < 0.5) sprintf("%.*f nm", digits, v[1])
    else sprintf("%.*f to %.*f nm", digits, min(v), digits, max(v))
  }

  n_samp <- sum(isTRUE_vec(ok$sampling_limited))
  n_prec <- nrow(ok) - n_samp
  regime <- if (n_samp == 0L) "Wrap set by localization precision."
    else if (n_prec == 0L)   "Wrap set by molecule spacing."
    else sprintf("Wrap set by localization precision on %d, by molecule spacing on %d.",
                 n_prec, n_samp)

  warn <- character(0)
  n_open <- sum(!is.na(ok$closed) & !ok$closed)
  if (n_open) warn <- c(warn, sprintf(
    "%d surface%s did not close into one body. The structure is sampled too sparsely to wrap this tightly; loosen the wrap or raise smoothing.",
    n_open, if (n_open == 1L) "" else "s"))
  n_unres <- sum(!isTRUE_vec(ok$resolved))
  if (n_unres) warn <- c(warn, sprintf(
    "%d surface%s hit the voxel cap and %s under-resolved. Gaps in %s may not be real.",
    n_unres, if (n_unres == 1L) "" else "s",
    if (n_unres == 1L) "is" else "are", if (n_unres == 1L) "it" else "them"))
  if (n_f) warn <- c(warn, sprintf(
    "%d cluster%s had too few molecules to wrap.", n_f, if (n_f == 1L) "" else "s"))

  shiny::div(
    class = "iso-status",
    style = paste("margin:6px 0 10px 0; padding:8px 10px; border-radius:4px;",
                  "background:rgba(127,127,127,0.10); font-size:12px; line-height:1.45;"),
    shiny::tags$div(shiny::tags$b(sprintf("%d surface%s drawn.", nrow(ok),
                                          if (nrow(ok) == 1L) "" else "s")),
                    sprintf(" Wrap %s, smoothing %s, voxel %s.",
                            rng(ok$wrap), rng(ok$smoothing), rng(ok$voxel, 1))),
    shiny::tags$div(style = "opacity:0.8;", regime),
    if (length(warn)) shiny::tags$div(
      style = "margin-top:5px; color:#b45309;",
      lapply(warn, function(w) shiny::tags$div(w))
    )
  )
}

# isTRUE() is not vectorised; these columns can carry NA (a diagnostic that was
# never tested), and NA must read as FALSE rather than propagate.
isTRUE_vec <- function(v) !is.na(v) & v

# =============================================================================
# cluster_alpha_volume() - ENCLOSED VOLUME OF ONE STRUCTURE, ALPHA CHOSEN FROM
#                          THAT STRUCTURE'S OWN SPACING
# =============================================================================
#
# GOAL
#   Return list(volume, alpha): the alpha-shape volume enclosed by a structure's
#   member centroids, plus the alpha radius actually used, so the number is
#   reportable rather than dependent on a hidden constant.
#
# WHY THE ALPHA CANNOT BE A CONSTANT
#   Alpha is a probe radius: gaps wider than it are left open, so a value below
#   the members' own spacing returns an empty shape and a volume of ZERO. The
#   previous fixed `alpha = 10` did exactly that - on a 150 nm shell whose
#   members sat 18.5 nm apart it returned 0, as it did for every closed
#   structure tested; that shell needs alpha ~ 200 nm to close. A 30 nm molecule
#   and a 500 nm vesicle cannot share one probe radius, which is the same
#   argument alpha_surface_3d() makes for the rendered surfaces.
#
# HOW ALPHA IS CHOSEN
#   Walk a ladder of radii built from the group's own median nearest-neighbour
#   spacing and take the SMALLEST that both (1) wraps every member in ONE
#   connected component and (2) yields a WATERTIGHT boundary - every edge shared
#   by exactly two faces. Same two tests alpha_surface_3d() uses for rendering.
#
#   Connectivity alone is not enough for SOLID structures. A filled cloud is
#   connected at a very small alpha, but the boundary is then pitted with
#   concavities that eat the volume: measured on a uniformly filled R=150 ball,
#   connectivity-only recovered 21.5% of the true volume (14.3% at n=1500 - it
#   gets WORSE with more points, since denser sampling connects at a still
#   smaller alpha), while adding the watertight test gave 75.1% and 83.8%.
#   Hollow shells are unaffected (96.0% either way), because they only connect
#   once alpha spans the gaps, which is also when they close.
#
# WHAT THE NUMBER MEANS
#   For a hollow structure this is the volume ENCLOSED BY the labels once alpha
#   bridges the gaps between them - the volume of the compartment they bound,
#   NOT the volume of labelled material. It also grows with alpha, so `alpha`
#   is returned alongside and both belong in any report.
cluster_alpha_volume <- function(pts) {
  none <- list(volume = NA_real_, alpha = NA_real_)
  if (is.null(pts)) return(none)
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  if (ncol(pts) < 3L) return(none)
  pts <- pts[stats::complete.cases(pts), 1:3, drop = FALSE]
  storage.mode(pts) <- "double"
  if (nrow(pts) < 4L) return(none)                 # a tetrahedron is the minimum

  # A group flat in any axis makes the tetrahedralization degenerate; nudge by
  # sub-picometre amounts (invisible at nm scale). Same guard as alpha_surface_3d().
  for (k in 1:3) {
    sdk <- stats::sd(pts[, k], na.rm = TRUE)
    if (!is.finite(sdk) || sdk < 1e-9) pts[, k] <- pts[, k] + seq_len(nrow(pts)) * 1e-6
  }

  nnd <- tryCatch({
    d <- RANN::nn2(pts, pts, k = 2L)$nn.dists[, 2]  # distance to nearest OTHER member
    d[is.finite(d) & d > 0]
  }, error = function(e) numeric(0))
  diag_nm <- sqrt(sum((apply(pts, 2, max) - apply(pts, 2, min))^2))
  base <- if (length(nnd)) stats::median(nnd) else diag_nm / 10
  if (!is.finite(base) || base <= 0) base <- diag_nm / 10
  if (!is.finite(base) || base <= 0) return(none)   # all members coincide

  alphas <- c(base * c(1, 1.5, 2, 3, 4, 6, 8, 12, 20, 40), diag_nm)
  alphas <- sort(unique(alphas[is.finite(alphas) & alphas > 0]))
  if (!length(alphas)) return(none)

  as3d <- suppressWarnings(tryCatch(alphashape3d::ashape3d(pts, alpha = alphas),
                                    error = function(e) NULL))
  if (is.null(as3d)) return(none)

  comp <- suppressWarnings(tryCatch(
    alphashape3d::components_ashape3d(as3d, indexAlpha = "all"),
    error = function(e) NULL))
  if (!is.null(comp) && !is.list(comp)) comp <- list(comp)

  # Watertight test: every edge of a closed mesh belongs to exactly two faces;
  # an edge used once is a hole rim. Same criterion as alpha_surface_3d().
  is_closed <- function(tt) {
    if (is.null(tt) || nrow(tt) < 4L) return(FALSE)
    ed  <- rbind(tt[, c(1, 2)], tt[, c(2, 3)], tt[, c(1, 3)])
    ed  <- cbind(pmin(ed[, 1], ed[, 2]), pmax(ed[, 1], ed[, 2]))   # undirected key
    key <- paste(ed[, 1], ed[, 2])
    all(tabulate(match(key, unique(key))) == 2L)
  }
  tri  <- if (!is.null(as3d$triang) && nrow(as3d$triang)) as.matrix(as3d$triang) else NULL
  # per-alpha classification columns are the LAST length(alphas) columns, in the
  # order the alphas were passed
  cls0 <- if (!is.null(tri)) ncol(tri) - length(alphas) else 0L

  for (ii in seq_along(alphas)) {
    cc <- if (!is.null(comp) && length(comp) == length(alphas)) comp[[ii]] else NULL
    # (1) one component covering every member == nothing shed at this radius
    if (!is.null(cc) && !(length(cc) && all(is.finite(cc)) && all(cc == 1L))) next
    if (ii < length(alphas)) {
      # (2) watertight boundary
      if (!is.null(tri) && cls0 >= 3L) {
        cand <- tri[tri[, cls0 + ii] == 2L, 1:3, drop = FALSE]      # REGULAR = boundary
        if (!nrow(cand) || !is_closed(cand)) next
      }
      # (3) the shape must actually CONTAIN the members. Tests (1) and (2) alone
      # admit a degenerate sliver: components_ashape3d() reports connectivity of
      # POINTS, and is_closed() accepts any closed sub-mesh, so a single
      # tetrahedron spanning 4 of the members passes both. Observed on a real
      # 31-member shell - alpha 100.4 gave a 4-triangle shape and a volume of
      # 3,155 nm^3 where the correct rung (alpha 150.5, 58 triangles) gave
      # 10,629,935. Requiring containment is the direct test and works for
      # hollow and solid structures alike, where a "fraction of points on the
      # boundary" rule would not.
      inside <- suppressWarnings(tryCatch(
        alphashape3d::inashape3d(as3d, indexAlpha = ii, points = pts),
        error = function(e) NULL))
      if (!is.null(inside) && !(length(inside) == nrow(pts) && all(inside))) next
    }
    v <- suppressWarnings(tryCatch(
      alphashape3d::volume_ashape3d(as3d, indexAlpha = ii),
      error = function(e) NA_real_))
    if (length(v) && is.finite(v[1]) && v[1] > 0)
      return(list(volume = as.numeric(v[1]), alpha = alphas[ii]))
  }
  none
}

# =============================================================================
# estimate_cluster2_curvature() - PER-STRUCTURE SIZE, CURVATURE AND VOLUME
# =============================================================================
#
# GOAL
#   One row per Median Cluster: how many members, how big (fitted sphere
#   radius), how curved (mean curvature) and what volume its members enclose.
#
# UNIT
#   Members are de-blinked molecules when de-blinking has run, else trace IDs -
#   see cluster2_member_col(). Each member is reduced to its MEDIAN position
#   first, so localization precision enters as sigma/sqrt(n) rather than as
#   the raw point cloud, and a long trace does not outvote a short one.
#
# CURVATURE
#   MEAN curvature H = 1/R in nm^-1, one convention app-wide (see
#   add_cluster2_lines_and_angles). For a fitted sphere this is a
#   reparametrisation of `sphere_radius` and carries no extra information; it is
#   reported because bending is conventionally quoted this way.
#
# WHAT WAS REMOVED AND WHY
#   The old `curvature` / `curvature_deg` columns were lambda3/sum(lambda), the
#   variance fraction on the smallest principal axis. That is a PLANARITY index
#   (the Westin anisotropy family), not curvature: measured, R = 50, 150 and
#   500 nm shells all scored ~0.30, so it could not distinguish a small vesicle
#   from a large one, and a flat disc (0.0002) scored the same as a rod (0.0010),
#   so it could not distinguish a plane from a line. `asin(sqrt(.))` in degrees
#   had no geometric meaning, and `sphere_curvature_deg` was the identical
#   formula under a name implying it came from the sphere fit.
estimate_cluster2_curvature <- function(df, member_col = NULL) {
  mcol <- cluster2_member_col(df, member_col)
  unit_lbl <- if (identical(mcol, "cluster")) "de-blinked molecule" else "trace ID"

  base <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0)
  if (identical(mcol, "cluster"))
    base <- base %>% dplyr::filter(!is.na(cluster) & cluster > 0L)

  base %>%
    dplyr::group_by(cluster2) %>%
    dplyr::group_modify(~{
      med <- .x %>%
        dplyr::group_by(.data[[mcol]]) %>%
        dplyr::summarise(x = stats::median(x, na.rm = TRUE),
                         y = stats::median(y, na.rm = TRUE),
                         z = stats::median(z, na.rm = TRUE), .groups = "drop")
      coords <- as.matrix(med[, c("x", "y", "z")])
      coords <- coords[stats::complete.cases(coords), , drop = FALSE]
      n_mem  <- nrow(coords)

      if (n_mem >= 5) {
        radius <- fit_sphere_geometry(coords)$radius     # geometric-median centre
        av     <- cluster_alpha_volume(coords)
        data.frame(
          member_unit   = unit_lbl,
          n_members     = n_mem,
          sphere_radius = radius,
          mean_curvature_inv_nm = if (is.finite(radius) && radius > 0) 1 / radius else NA_real_,
          volume        = av$volume,
          volume_alpha  = av$alpha
        )
      } else {
        data.frame(
          member_unit   = unit_lbl,
          n_members     = n_mem,
          sphere_radius = NA_real_,
          mean_curvature_inv_nm = NA_real_,
          volume        = NA_real_,
          volume_alpha  = NA_real_
        )
      }
    }) %>%
    dplyr::ungroup()
}

# =============================================================================
# build_curvature_table() - THE EXPORTED PER-STRUCTURE TABLE
# =============================================================================
#
# GOAL
#   One row per Median Cluster, joining membership/identity columns to the
#   geometry from estimate_cluster2_curvature().
#
# COUNTS ARE REPORTED SEPARATELY AND LABELLED
#   n_traces      = blinking events (trace IDs) - photophysics-dependent
#   n_molecules   = de-blinked molecules, NA when de-blinking has not run
#   n_members     = the unit the geometry was actually measured on (see
#                   `member_unit`), so the table never implies the geometry used
#                   a unit it did not.
#   n_points      = raw localizations
build_curvature_table <- function(df) {
  curvature_df <- estimate_cluster2_curvature(df)
  has_mol <- "cluster" %in% names(df) && any(!is.na(df$cluster) & df$cluster > 0L)

  merged <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
    dplyr::group_by(cluster2) %>%
    dplyr::summarize(
      trace_ids   = paste(unique(tid), collapse = ", "),
      molecules   = if (has_mol) paste(unique(cluster[!is.na(cluster) & cluster > 0L]), collapse = ", ") else NA_character_,
      n_traces    = dplyr::n_distinct(tid),
      n_molecules = if (has_mol) dplyr::n_distinct(cluster[!is.na(cluster) & cluster > 0L]) else NA_integer_,
      n_points    = dplyr::n(),
      # Structure position on the same footing as the geometry: the geometric
      # median of the MEMBER centroids, not a coordinate-wise median of every
      # raw localization (which weights long traces and is not rotation-equivariant).
      .groups = "drop"
    )

  mcol <- cluster2_member_col(df)
  centres <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
    { if (identical(mcol, "cluster")) dplyr::filter(., !is.na(cluster) & cluster > 0L) else . } %>%
    dplyr::group_by(cluster2, .data[[mcol]]) %>%
    dplyr::summarise(mx = stats::median(x, na.rm = TRUE),
                     my = stats::median(y, na.rm = TRUE),
                     mz = stats::median(z, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(cluster2) %>%
    dplyr::group_modify(~{
      g <- geometric_median(as.matrix(.x[, c("mx", "my", "mz")]))
      data.frame(x = g[1], y = g[2], z = g[3])
    }) %>%
    dplyr::ungroup()

  out <- merged %>%
    dplyr::left_join(centres, by = "cluster2") %>%
    dplyr::left_join(curvature_df, by = "cluster2")

  # ── SMALL-N RELIABILITY FLAG ────────────────────────────────────────────────
  # Every shape metric here is estimated from the cluster's own members, and all of
  # them are biased at low member counts - in DIFFERENT directions, which is why a
  # single "n" column is not enough for a reader to judge them. Measured against
  # known ground truth (see docs and the code notes on estimate_cluster2_curvature):
  #
  #   eccentricity - PCA orders axes by SAMPLE variance, so the ordering itself
  #     manufactures elongation. On a TRUE SPHERE the median ecc2 is 0.81 at n = 10,
  #     0.53 at n = 40, 0.34 at n = 200 and 0.18 at n = 2000. Never read an
  #     eccentricity without comparing it to that floor for the member count.
  #   alpha-shape volume - flat facets under-fill a curved surface: 29.9% of truth
  #     at n = 10, 54.6% at n = 20, 67.6% at n = 31, 86.2% at n = 80, 96.0% at
  #     n = 300. Below ~20 members it is not usable.
  #   sphere_radius - the centre is fitted to the same points, so the radius is
  #     biased low: 147.4 nm at n = 40 against a truth of 150, 149.9 at n = 1000.
  #
  # The flag is advisory and deliberately does not blank anything: which metric the
  # reader cares about decides what counts as "enough", and silently dropping values
  # would hide the size dependence rather than expose it.
  nm_col <- if ("n_members" %in% names(out)) out$n_members else
            if ("n_molecules" %in% names(out) && any(!is.na(out$n_molecules))) out$n_molecules else
            out$n_traces
  out$n_members_used <- nm_col
  out$smallN_flag <- dplyr::case_when(
    is.na(nm_col)  ~ NA_character_,
    nm_col <  20L  ~ "severe (volume <55% of truth; eccentricity floor >0.5)",
    nm_col <  40L  ~ "strong (volume ~68-86% of truth; eccentricity floor ~0.5)",
    nm_col < 200L  ~ "moderate (eccentricity floor ~0.3-0.5)",
    TRUE           ~ "ok"
  )
  # Expected eccentricity of a TRUE SPHERE at this member count, log-interpolated
  # from the measured floor - the number to compare an observed eccentricity against.
  kn <- c(10, 40, 200, 2000); kv <- c(0.81, 0.53, 0.34, 0.18)
  out$ecc_sphere_floor <- ifelse(is.na(nm_col), NA_real_,
    stats::approx(log(kn), kv, xout = log(pmax(nm_col, 2)), rule = 2)$y)
  out
}


# =============================================================================
# PANEL HELPERS (moved from server_smlm.R; plain functions of their arguments)
# =============================================================================

# Curve overlays (ring/arc/helix): ordered points to draw as a connected line.
.viewer_curve_points <- function(f, shape) {
  pts <- tryCatch(switch(shape,
    ring = { p <- as.matrix(mle_sample_ring(f, n = 240)); rbind(p, p[1, ]) },  # close loop
    arc  = as.matrix(mle_sample_arc(f, n = 240)),
    helix = {
      if (!is.null(f$t)) {
        tr <- range(f$t, na.rm = TRUE)
        ts <- seq(tr[1], tr[2], length.out = 600L)
        t(vapply(ts, .helix_point, numeric(3),
                 c = f$center, u = f$axis, r = f$radius, p = f$pitch, phi0 = f$phi0))
      } else if (!is.null(f$H)) as.matrix(f$H) else NULL
    },
    NULL), error = function(e) NULL)
  if (is.null(pts) || !nrow(pts)) return(NULL)
  pts <- as.matrix(pts)[, 1:3, drop = FALSE]; colnames(pts) <- c("x", "y", "z"); pts
}

# Surface overlays (sphere/ellipsoid): a (theta,phi) vertex grid + triangle
# faces for a translucent mesh3d in 3D; the projected convex hull is the exact
# orthographic silhouette outline for these convex bodies in 2D.
.viewer_surface_mesh <- function(f, shape, nu = 48L, nv = 28L) {
  res <- tryCatch({
    uu <- rep(seq(0, 2*pi, length.out = nu), times = nv)
    vv <- rep(seq(0, pi,  length.out = nv), each  = nu)
    base <- cbind(sin(vv) * cos(uu), sin(vv) * sin(uu), cos(vv))
    verts <- switch(shape,
      sphere    = base * f$radius,
      ellipsoid = sweep(base, 2, f$axes, "*") %*% t(f$R),
      NULL)
    if (is.null(verts)) stop("unsupported")
    verts <- sweep(verts, 2, f$center, "+")
    iu <- rep(seq_len(nu - 1L), times = nv - 1L)
    iv <- rep(seq_len(nv - 1L), each  = nu - 1L)
    a <- (iv - 1L) * nu + (iu - 1L); b <- (iv - 1L) * nu + iu     # 0-based vertex ids
    cc <- iv * nu + iu;              d <- iv * nu + (iu - 1L)
    list(verts = verts, i = c(a, a), j = c(b, cc), k = c(cc, d))
  }, error = function(e) NULL)
  res
}

