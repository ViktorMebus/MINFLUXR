## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_gt.R  -  Ground-truth estimate and structure size (exploratory)
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The inverse of the labelling model (candidate molecule positions, two
##   molecule counts and their consistency test), the label multiplicity
##   from a counting standard, and structure shape and size descriptors.
##
## CONTENTS (in file order)
##   measure_label_multiplicity, gt_tid_model, estimate_gt_from_labeling,
##   gt_uncertainty_colors, mst_contour_length, .grid_thin,
##   analyze_cluster2_morphology
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

## ── Empirical label-multiplicity calibration (counting standard) ──────────────
# Measures the observed labelling multiplicity from a CALIBRATION STANDARD of
# known stoichiometry, so the count estimator can divide by a measured value
# instead of the preset guess m_tid = n_fluor * eff * (1 + blink_lambda).
#
# The standard is a sparse sample of isolated structures each holding a KNOWN
# number of true molecules (`k_per_structure`; 1 for single, well-separated
# molecules such as sparse DNA docking strands). It MUST be imaged with the same
# label and acquisition as the sample being counted, and clustered with the same
# settings, so each secondary structure (cluster2 > 0) is one standard object.
#
# Observed TIDs per molecule = median over structures of (TIDs in the structure)
# / k. That IS the m_tid the count divides by, measured rather than assumed. The
# per-structure spread (MAD) reports how reproducible the label is.
#
# Returns a list; ok = FALSE (never an error) when nothing is measurable.
measure_label_multiplicity <- function(df, k_per_structure = 1L, minpts = 3L) {
  out <- list(ok = FALSE, n_struct = 0L, k = as.integer(k_per_structure),
              m_tid = NA_real_, m_tid_mad = NA_real_, m_tid_iqr = c(NA_real_, NA_real_),
              m_prim = NA_real_, m_loc = NA_real_, per_struct = data.frame())
  if (is.null(df) || !nrow(df) || !"tid" %in% names(df) || !"cluster2" %in% names(df))
    return(out)
  k <- max(1L, as.integer(k_per_structure))
  has_prim <- "cluster" %in% names(df)
  d <- df[!is.na(df$cluster2) & df$cluster2 > 0L, , drop = FALSE]
  if (!nrow(d)) return(out)
  per <- d %>%
    dplyr::group_by(cluster2) %>%
    dplyr::summarise(
      n_loc  = dplyr::n(),
      n_tid  = dplyr::n_distinct(tid),
      n_prim = if (has_prim) dplyr::n_distinct(cluster[!is.na(cluster) & cluster > 0L]) else NA_integer_,
      .groups = "drop") %>%
    dplyr::filter(n_loc >= minpts)
  if (!nrow(per)) return(out)
  out$ok         <- TRUE
  out$n_struct   <- nrow(per)
  out$k          <- k
  out$m_tid      <- stats::median(per$n_tid, na.rm = TRUE) / k
  out$m_tid_mad  <- stats::mad(per$n_tid / k, na.rm = TRUE)
  out$m_tid_iqr  <- as.numeric(stats::quantile(per$n_tid / k, c(0.25, 0.75), na.rm = TRUE))
  out$m_prim     <- if (has_prim) stats::median(per$n_prim, na.rm = TRUE) / k else NA_real_
  out$m_loc      <- stats::median(per$n_loc, na.rm = TRUE) / k
  out$per_struct <- as.data.frame(per)
  out
}

# =============================================================================
# gt_tid_model() - HOW MANY TRACE IDS ONE MOLECULE YIELDS UNDER A LABEL PRESET
# =============================================================================
#
# GOAL
#   The exact distribution of the trace-ID count T of one molecule under the same
#   forward model the Simulator draws from (docs/simulator.md 1.2-1.4):
#       K  ~ max(1, Poisson(n_fluor))     dyes on the molecule (K = 1 if n_fluor <= 1)
#       K' ~ Binomial(K, eff)             dyes detected
#       T  = sum over detected dyes of (1 + Poisson(blink_lambda))
#   so the count estimator, its uncertainty and the crowding test all rest on one
#   model instead of separate approximations.
#
# WHY TWO SETS OF MOMENTS
#   "unconditional" moments count undetected molecules as T = 0. Dividing a
#   structure's total trace count by mean_u estimates ALL molecules, detected or
#   not (the stoichiometric count). "detected" moments condition on K' >= 1; they
#   describe one molecule that produced at least one trace, which is what a spatial
#   candidate is.
#
# ARGUMENTS  n_fluor, eff (fraction 0-1), blink_lambda : the preset values
# RETURNS    list(eta, mean_u, var_u, mean_d, var_d, mean_dyes, p_ge)
#   eta       P(K' >= 1), the molecule detection probability
#   mean_dyes E[K'], detected dyes per molecule (unconditional)
#   p_ge(t)   P(T >= t | detected), vectorised over t
# =============================================================================
gt_tid_model <- function(n_fluor, eff, blink_lambda) {
  p   <- min(1, max(1e-6, eff))
  lam <- max(0, blink_lambda)
  if (!is.finite(n_fluor) || n_fluor <= 1) {
    k <- 1L; pk <- 1
  } else {
    kmax <- max(2L, as.integer(stats::qpois(1 - 1e-12, n_fluor)) + 5L)
    k  <- seq_len(kmax)
    pk <- stats::dpois(k, n_fluor)
    pk[1] <- pk[1] + stats::dpois(0, n_fluor)           # max(1, .) folds 0 into 1
  }
  # Detected-dye count: a mixture of binomials over the dye count.
  kd  <- 0:max(k)
  pkd <- vapply(kd, function(j) sum(pk * stats::dbinom(j, k, p)), numeric(1))
  eta <- 1 - pkd[1]
  # Given K' = j, T = j + Poisson(j * lam): mean j(1+lam), variance j*lam.
  m_j  <- kd * (1 + lam)
  e2_j <- kd * lam + m_j^2
  mean_u <- sum(pkd * m_j)
  var_u  <- sum(pkd * e2_j) - mean_u^2
  pc <- pkd[-1] / max(eta, 1e-12); kc <- kd[-1]
  mean_d <- sum(pc * m_j[-1])
  var_d  <- sum(pc * e2_j[-1]) - mean_d^2
  p_ge <- function(t) vapply(t, function(tt)
    sum(pc * stats::ppois(tt - kc - 1, kc * lam, lower.tail = FALSE)), numeric(1))
  list(eta = eta, mean_u = mean_u, var_u = max(0, var_u),
       mean_d = mean_d, var_d = max(0, var_d),
       mean_dyes = sum(pkd * kd), p_ge = p_ge)
}

## ── Inverse GT-position estimator (labeling-aware) ────────────────────────────
# Reverse-engineers the Simulator's forward labeling model to turn observed
# traces/clusters back into CANDIDATE ground-truth molecule positions, an
# approximate molecule count, and a crowding validity flag.
#
# Scientific model (inverts helpers_shared::minflux_label_presets() +
# sim_server's expand_to_fluorophores / split_for_blinking / gen_localizations):
#
#   one GT molecule
#     --(link_nm)-->    ONE common-mode complex offset, shared by all its dyes
#     --(n_fluor)-->    Poisson(n_fluor) fluorophores, each scattered dye_spread_nm
#     --(blink)-->      each fluorophore -> (1 + Poisson(blink_lambda)) TIDs
#     --(sigma_loc)-->  each TID -> many localizations, precision sigma_loc
#
# So the labeling multiplicity at the RAW-TID level is
#     m_tid = E[detected dyes] * (1 + blink_lambda)        (TIDs per molecule)
# computed exactly by gt_tid_model(), together with its variance; and the DE-BLINKING step (merge threshold ~ dye_spread) collapses re-blinks and
# the dyes of one molecule back into ~one object => ~one unit per molecule,
# turning an unstable /m_tid calibration into a /~1 count. De-blinked molecules
# are therefore the preferred unit when present ("auto").
#
# KEY LIMIT: the link_nm offset is COMMON-MODE within a molecule (all its
# localizations share the one draw), so it never averages out - it is the
# irreducible position-uncertainty floor. Reported as `pos_sd_nm`.
#
# Args:
#   df       clustered SMLM frame (needs x,y[,z],tid; uses `cluster`/`cluster2`)
#   preset   a preset list from minflux_label_presets(), OR its name, OR a list
#            with n_fluor/link_nm/dye_spread_nm/footprint_nm/blink_lambda
#   eff_pct  labeling/detection efficiency (%) - folds into the multiplicity
#   level    "auto" (de-blinked molecules when available else TIDs), "primary", "tid"
#            ("primary" is the stored VALUE for the de-blinked level; kept for
#            settings compatibility, it does not mean a primary DBSCAN)
#   fuse_k   a unit whose per-axis spread exceeds fuse_k x the single-molecule
#            expectation is flagged as fused (multiple GTs merged)
#   minpts   localizations below which a unit is ignored
#   sigma_loc_nm  per-localization precision; estimated from data when NULL
#   scope_cluster2 restrict to points inside a secondary structure (cluster2>0)
#
# Returns list(points, links, summary, meta) - see field construction below. Empty
# data.frames (never NULL) when nothing is estimable, so callers stay simple.
#   points$count_fused      candidate holds more TIDs than one molecule gives (p < 0.01)
#   summary$count_verdict   "consistent" / "under-resolved" / "over-split", from the
#                           TID total given the spatial count (count_z); under-resolved
#                           also sets `crowded`, any non-consistent verdict sets
#                           `stoich_mismatch`
#   summary$n_gt_lo / _hi   range spanning both count estimators with their errors
estimate_gt_from_labeling <- function(df, preset, eff_pct = 70,
                                      level = c("auto", "primary", "tid"),
                                      fuse_k = 2.5, minpts = 3L,
                                      sigma_loc_nm = NULL,
                                      scope_cluster2 = TRUE,
                                      complex_r = NULL,
                                      displace_inward = TRUE,
                                      displace_nm = NULL,
                                      guided_merge = FALSE,
                                      support_strength = 1,
                                      m_tid_override = NULL) {
  level <- match.arg(level)
  empty_pts <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0),
                          pos_sd_nm = numeric(0), n_loc = integer(0),
                          n_tid = integer(0), n_units = integer(0),
                          fused = logical(0), cluster2 = integer(0),
                          mol = integer(0), level = character(0))
  empty_links <- data.frame(mol = integer(0), ux = numeric(0), uy = numeric(0),
                            uz = numeric(0), gx = numeric(0), gy = numeric(0),
                            gz = numeric(0))
  empty <- list(points = empty_pts, links = empty_links,
                summary = data.frame(), meta = list())
  if (is.null(df) || !nrow(df) || !all(c("x", "y", "tid") %in% names(df)))
    return(empty)

  # ── resolve preset parameters ──
  if (is.character(preset)) preset <- minflux_label_presets()[[preset]]
  gp <- function(k, d) { v <- suppressWarnings(as.numeric(preset[[k]])); if (length(v) && is.finite(v)) v else d }
  n_fluor      <- max(1,  gp("n_fluor", 1))
  link_nm      <- max(0,  gp("link_nm", 5))
  dye_spread   <- max(0,  gp("dye_spread_nm", 0))
  footprint    <- max(0,  gp("footprint_nm", link_nm * 1.5))
  blink_lambda <- max(0,  gp("blink_lambda", 0.3))
  eff          <- min(1, max(0.01, eff_pct / 100))
  preset_note  <- if (!is.null(preset$note)) as.character(preset$note) else ""

  # TIDs expected per molecule (m_tid) and the detection probability eta, both
  # from the exact forward model (gt_tid_model). A measured multiplicity from a
  # counting standard replaces the preset's mean further down; the model's blink
  # rate is then re-solved so its mean matches, keeping the shape of the
  # distribution for the uncertainty and crowding tests.
  tid_mdl <- gt_tid_model(n_fluor, eff, blink_lambda)
  m_tid   <- tid_mdl$mean_u
  eta     <- min(1, max(1e-3, tid_mdl$eta))

  use_z <- "z" %in% names(df) && any(is.finite(df$z) & df$z != 0)

  # ── estimate per-localization precision from data when not supplied ──
  # Each TID's within-trace coordinate SD ~ sigma_loc; take the median over TIDs.
  if (is.null(sigma_loc_nm) || !is.finite(sigma_loc_nm)) {
    per_tid_sd <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(sx = stats::sd(x, na.rm = TRUE),
                       sy = stats::sd(y, na.rm = TRUE),
                       n  = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n >= 3L, is.finite(sx), is.finite(sy))
    sigma_loc_nm <- if (nrow(per_tid_sd))
      stats::median(sqrt((per_tid_sd$sx^2 + per_tid_sd$sy^2) / 2), na.rm = TRUE) else 2
    if (!is.finite(sigma_loc_nm) || sigma_loc_nm <= 0) sigma_loc_nm <- 2
  }
  # expected per-axis spread of a SINGLE molecule's localizations (dye scatter + noise)
  r_single <- sqrt(dye_spread^2 + sigma_loc_nm^2)

  # ── choose the unit level ──
  has_primary <- "cluster" %in% names(df) && any(!is.na(df$cluster) & df$cluster > 0L)
  unit_level  <- switch(level,
                        auto    = if (has_primary) "primary" else "tid",
                        primary = "primary",
                        tid     = "tid")
  if (unit_level == "primary" && !has_primary) {
    if (level == "primary") return(empty)      # explicitly asked for primary, none exist
    unit_level <- "tid"
  }

  # ── optionally scope to points inside a secondary structure ──
  has_c2 <- "cluster2" %in% names(df) && any(!is.na(df$cluster2) & df$cluster2 > 0L)
  work <- df
  if (isTRUE(scope_cluster2) && has_c2)
    work <- work[!is.na(work$cluster2) & work$cluster2 > 0L, , drop = FALSE]
  if (!nrow(work)) return(empty)
  if (!"cluster2" %in% names(work)) work$cluster2 <- 0L
  work$cluster2[is.na(work$cluster2)] <- 0L
  if (unit_level == "primary")
    work <- work[!is.na(work$cluster) & work$cluster > 0L, , drop = FALSE]
  if (!nrow(work)) return(empty)

  unit_col <- if (unit_level == "primary") "cluster" else "tid"

  # ── build candidate units ──
  if (!"z" %in% names(work)) work$z <- 0    # ensure column exists; zeroed later if 2D
  work$z[!is.finite(work$z)] <- 0
  # NOTE: compute the SDs BEFORE the median columns overwrite x/y/z - dplyr
  # summarise evaluates sequentially, so `sd_x` must reference the raw x, not the
  # just-created scalar median (which would give NA).
  units <- work %>%
    dplyr::group_by(.uid = .data[[unit_col]]) %>%
    dplyr::summarise(
      sd_x     = stats::sd(x, na.rm = TRUE),
      sd_y     = stats::sd(y, na.rm = TRUE),
      sd_z     = if (use_z) stats::sd(z, na.rm = TRUE) else 0,
      n_loc    = dplyr::n(),
      n_tid    = dplyr::n_distinct(tid),
      cluster2 = { v <- cluster2[cluster2 > 0L]; if (length(v)) as.integer(names(which.max(table(v)))) else 0L },
      x        = stats::median(x, na.rm = TRUE),
      y        = stats::median(y, na.rm = TRUE),
      z        = if (use_z) stats::median(z, na.rm = TRUE) else 0,
      .groups  = "drop"
    ) %>%
    dplyr::filter(n_loc >= minpts, is.finite(x), is.finite(y))
  if (!nrow(units)) return(empty)

  units$sd_x[!is.finite(units$sd_x)] <- 0
  units$sd_y[!is.finite(units$sd_y)] <- 0
  units$sd_z[!is.finite(units$sd_z)] <- 0

  # observed per-axis spread of each unit (RMS across axes)
  ndim        <- if (use_z) 3 else 2
  units$unit_spread <- sqrt((units$sd_x^2 + units$sd_y^2 +
                             if (use_z) units$sd_z^2 else 0) / ndim)

  # ── MERGE units into candidate molecules ──────────────────────────────────
  # A multivalent label (n_fluor>1, dye_spread>0) scatters several fluorophores,
  # hence several de-blinked molecules / TIDs, around ONE molecule. Pool unit
  # centroids that fall within `complex_r` (the label-complex radius) so those
  # contributors collapse into a single candidate GT - this is what makes
  # "many traces -> one molecule" happen, and it produces the contribution
  # fan-out lines. Merging is confined within each cluster2 so distinct
  # structures never pool together.
  #
  # AUTO radius is PRESET-SPECIFIC: max(footprint, 2.5 x dye_spread). Both terms
  # come from the labeling preset, so a big polyclonal+secondary complex pools
  # over ~45 nm while a nanobody pools over ~3 nm - i.e. the labeling strategy
  # itself sets how many traces collapse into one molecule, as expected.
  if (is.null(complex_r) || !is.finite(complex_r) || complex_r < 0)
    complex_r <- max(footprint, 2.5 * dye_spread)
  coord_cols <- if (use_z) c("x", "y", "z") else c("x", "y")

  # `m_tid` = expected TIDs per molecule from the labeling strategy; the
  # stoichiometric molecule count of a structure is (its TID count) / m_tid.
  # A measured multiplicity from a counting standard (m_tid_override) replaces the
  # preset guess when supplied, so the count divides by a calibrated value.
  m_tid_measured <- !is.null(m_tid_override) && is.finite(m_tid_override) && m_tid_override > 0
  tids_per_mol <- if (m_tid_measured) as.numeric(m_tid_override) else max(1e-6, m_tid)
  if (m_tid_measured && tid_mdl$mean_dyes > 0) {
    # mean_u = E[K'] * (1 + lambda), so the blink rate that reproduces the measured
    # multiplicity is m / E[K'] - 1 (floored at 0: a standard below one trace per
    # detected dye means the efficiency, not the blinking, is off).
    tid_mdl <- gt_tid_model(n_fluor, eff, max(0, tids_per_mol / tid_mdl$mean_dyes - 1))
  }

  units$mol <- NA_integer_
  next_mol  <- 1L
  for (sid in unique(units$cluster2)) {
    ix <- which(units$cluster2 == sid)
    if (length(ix) == 1L || complex_r <= 0) {
      units$mol[ix] <- next_mol + seq_along(ix) - 1L
      next_mol <- next_mol + length(ix); next
    }
    cc <- as.matrix(units[ix, coord_cols, drop = FALSE])
    g  <- tryCatch(dbscan::dbscan(cc, eps = complex_r, minPts = 1L)$cluster,
                   error = function(e) seq_along(ix))     # minPts=1 -> no noise
    g  <- as.integer(g)

    # ── stoichiometry-guided merge (opt-in, MODEL-DRIVEN) ──
    # If the spatial merge left MORE groups than the labeling stoichiometry
    # predicts (N_TID / m_tid), collapse further to that target via hierarchical
    # clustering. Only ever REDUCES the count toward stoichiometry - it never
    # splits (that would fabricate molecules). Confined to real structures.
    if (isTRUE(guided_merge) && sid > 0L) {
      n_grp    <- length(unique(g))
      k_stoich <- max(1L, round(sum(units$n_tid[ix]) / tids_per_mol))
      if (k_stoich < n_grp) {
        hc <- tryCatch(stats::hclust(stats::dist(cc), method = "average"),
                       error = function(e) NULL)
        if (!is.null(hc)) g <- as.integer(stats::cutree(hc, k = k_stoich))
      }
    }

    g <- g + next_mol - 1L
    units$mol[ix] <- g
    next_mol <- max(g) + 1L
  }

  # ── aggregate contributing units into molecules ───────────────────────────
  # Member-centroid SDs are computed BEFORE the weighted-mean position columns
  # (sequential-summarise trap), then the n_loc-weighted molecule position.
  # Order matters (dplyr summarise is sequential): SD/spread and the n_loc-
  # weighted positions must be computed while x/y/z and n_loc still hold the raw
  # per-unit vectors; the summed n_tid/n_loc come last so they don't shadow them.
  mols <- units %>%
    dplyr::group_by(mol) %>%
    dplyr::summarise(
      cluster2 = cluster2[1],
      n_units  = dplyr::n(),
      msd_x    = if (dplyr::n() > 1L) stats::sd(x) else 0,
      msd_y    = if (dplyr::n() > 1L) stats::sd(y) else 0,
      msd_z    = if (use_z && dplyr::n() > 1L) stats::sd(z) else 0,
      inner_sp = stats::median(unit_spread, na.rm = TRUE),
      x        = stats::weighted.mean(x, pmax(1, n_loc)),
      y        = stats::weighted.mean(y, pmax(1, n_loc)),
      z        = if (use_z) stats::weighted.mean(z, pmax(1, n_loc)) else 0,
      n_tid    = sum(n_tid),
      n_loc    = sum(n_loc),
      .groups  = "drop"
    )
  mols$msd_x[!is.finite(mols$msd_x)] <- 0
  mols$msd_y[!is.finite(mols$msd_y)] <- 0
  mols$msd_z[!is.finite(mols$msd_z)] <- 0
  mols$inner_sp[!is.finite(mols$inner_sp)] <- 0

  # spread of the fluorophore cloud that formed this molecule (per-axis RMS)
  member_spread <- sqrt((mols$msd_x^2 + mols$msd_y^2 +
                         if (use_z) mols$msd_z^2 else 0) / ndim)

  # Total spatial extent of a molecule = its members' fan-out (member_spread)
  # combined with the internal spread of those members (inner_sp). This catches
  # BOTH failure modes: several units spread wide (crowding merge) AND a single
  # over-wide de-blinked molecule that already chained many traces together.
  mol_extent <- sqrt(member_spread^2 + mols$inner_sp^2)

  # ── fusion: a molecule far wider than the expected single-molecule label
  # cloud is really several merged molecules (crowding / under-resolution).
  # Expected cloud = the dye-spread scatter plus localization noise; a tightly
  # labeled single molecule stays well under fuse_k x that.
  expected_cloud <- max(dye_spread, r_single, 1e-6)
  mols$fused <- mol_extent > (fuse_k * expected_cloud)

  # ── displace GT inward toward the local structure (linker correction) ─────
  # In the forward model the antibody binds a surface epitope and the linker
  # throws the fluorophore OUTWARD, so a single molecule's trace-cloud centroid
  # sits ~link_nm OUTSIDE the true structure. Push each (non-fused) molecule back
  # toward its structure by `displace_nm` (auto = link_nm). Geometry-aware: for
  # an elongated structure (a filament) "inward" is PERPENDICULAR toward the axis
  # (the axial component is dropped - there is no inside along the length); for a
  # compact / ring structure it points straight at the centroid. FUSED molecules
  # are skipped: their centroid geometry is already corrupted by merging, so the
  # clean outward-offset assumption doesn't hold. Applied only within a defined
  # structure (cluster2 > 0), and it never overshoots the centre.
  disp <- if (is.null(displace_nm) || !is.finite(displace_nm)) link_nm else max(0, displace_nm)
  meta_disp <- 0
  if (isTRUE(displace_inward) && disp > 0 && nrow(mols) >= 2L) {
    mols$disp_step <- 0                          # how far each GT was pushed (nm)
    mols$dir_conf  <- 0.4                         # confidence in the inward direction
    for (sid in unique(mols$cluster2)) {
      if (sid <= 0L) next                       # need a defined structure to aim at
      ix <- which(mols$cluster2 == sid)
      if (length(ix) < 3L) next
      M <- as.matrix(mols[ix, coord_cols, drop = FALSE])
      w <- pmax(1, mols$n_loc[ix])
      C <- colSums(M * w) / sum(w)              # n_loc-weighted structure centre
      axis1 <- NULL; elongated <- FALSE
      cv <- tryCatch(stats::cov(M), error = function(e) NULL)
      if (!is.null(cv)) {
        ev   <- eigen(cv, symmetric = TRUE)
        lam  <- pmax(ev$values, 0)
        rest <- if (length(lam) > 1) mean(lam[-1]) else 0
        if (lam[1] > 0 && rest > 0 && lam[1] / rest > 3) {
          elongated <- TRUE; axis1 <- ev$vectors[, 1]
        }
      }
      # inward direction is well-defined for an elongated structure (perp to a
      # clear axis) but ambiguous for a compact blob -> lower directional conf.
      mols$dir_conf[ix] <- if (elongated) 0.8 else 0.4
      for (k in ix) {
        if (isTRUE(mols$fused[k])) next          # don't displace ambiguous fused GTs
        p  <- as.numeric(mols[k, coord_cols])
        d  <- C - p
        if (elongated && !is.null(axis1)) d <- d - sum(d * axis1) * axis1  # perp to axis
        nd <- sqrt(sum(d^2))
        if (nd < 1e-6) next
        step <- min(disp, nd)                   # cap: never cross the centre
        newp <- p + step * d / nd
        mols$x[k] <- newp[1]; mols$y[k] <- newp[2]
        if (use_z) mols$z[k] <- newp[3]
        mols$disp_step[k] <- step
      }
    }
    meta_disp <- disp
  }
  if (!"disp_step" %in% names(mols)) { mols$disp_step <- 0; mols$dir_conf <- 1 }

  # ── position uncertainty (this drives the green->red colour) ──────────────
  # Independent contributions added in quadrature (all in nm):
  #  s_link    linker floor - irreducible, common-mode per molecule (dominant)
  #  s_sample  how well the label centroid is pinned = fan-out spread averaged
  #            over the contributing units, plus localization precision / #locs
  #  s_support RELIABILITY penalty for few independent labels: a GT from 1 TID
  #            is far less trustworthy than one from many (~link_nm/sqrt(n_tid))
  #  s_disp    the inward displacement guessed a direction; its error scales with
  #            how far we pushed and how ambiguous the direction was (1-dir_conf)
  #  s_fused   full extent when several molecules were merged into one GT
  s_link    <- link_nm
  s_sample  <- sqrt(member_spread^2 / pmax(1, mols$n_units) +
                    sigma_loc_nm^2   / pmax(1, mols$n_loc))
  s_support <- max(0, support_strength) * link_nm / sqrt(pmax(1, mols$n_tid))
  s_disp    <- mols$disp_step * (1 - mols$dir_conf)
  s_fused   <- ifelse(mols$fused, mol_extent, 0)
  mols$pos_sd_nm <- sqrt(s_link^2 + s_sample^2 + s_support^2 + s_disp^2 + s_fused^2)

  # ── per-candidate count check ──
  # P(one detected molecule yields >= this many TIDs). A small value means the
  # candidate carries more traces than one molecule plausibly gives, i.e. several
  # unresolved molecules. This catches merges that the spread test cannot see:
  # de-blinking caps a unit's width at its merge threshold, so a unit holding two
  # molecules a few nm apart is no wider than a single one.
  mols$count_p <- tid_mdl$p_ge(mols$n_tid)

  points <- data.frame(
    x = mols$x, y = mols$y, z = if (use_z) mols$z else 0,
    pos_sd_nm = mols$pos_sd_nm, n_loc = mols$n_loc, n_tid = mols$n_tid,
    n_units = mols$n_units, fused = mols$fused,
    count_p = mols$count_p, count_fused = mols$count_p < 0.01,
    cluster2 = mols$cluster2, mol = mols$mol, level = unit_level,
    stringsAsFactors = FALSE
  )

  # ── contribution links: one row per member unit -> its molecule centroid ──
  links <- data.frame(mol = units$mol, ux = units$x, uy = units$y,
                      uz = if (use_z) units$z else 0, stringsAsFactors = FALSE)
  links <- merge(links, mols[, c("mol", "x", "y", "z")], by = "mol", all.x = TRUE)
  names(links)[match(c("x", "y", "z"), names(links))] <- c("gx", "gy", "gz")
  # keep only real fan-outs (molecules with >=2 contributing units)
  multi <- mols$mol[mols$n_units >= 2L]
  links <- links[links$mol %in% multi, , drop = FALSE]

  # ── per-molecule count contribution (fused = several TIDs-per-molecule) ────
  mol_est <- ifelse(mols$fused, pmax(1, round(mols$n_tid / tids_per_mol)), 1)

  # ── per-structure summary (grouped by cluster2; 0 = ungrouped/global) ─────
  # Two INDEPENDENT molecule-count estimators per structure:
  #   spatial        = merged GT groups, efficiency-corrected (/eta)
  #   stoichiometric = total TIDs / m_tid  (already efficiency-aware via eff)
  #
  # CONSISTENCY TEST. If the spatial candidates are the true detected molecules,
  # the structure's total TID count is a sum of that many draws from the
  # detected-molecule distribution: mean n_det * mean_d, variance n_det * var_d.
  #   z >  z_crit  more traces than those molecules can carry -> under-resolved
  #                (candidates hold several molecules; crowding)
  #   z < -z_crit  fewer traces than expected -> over-split (re-blinks or dyes of
  #                one molecule left as separate candidates)
  # The variance is floored at the mean (Poisson): real blink counts are more
  # variable than the model, and a test at the bare model variance would flag
  # ordinary data. A ratio rule (the previous "outside 0.5-2x") ignores sample
  # size: it misses a large structure that is off by 40% and flags a 3-molecule
  # structure on noise.
  #
  # COUNT RANGE. Both estimators carry their own error and the range spans both:
  # spatial +- 1.96 SE (binomial detection), stoichiometric +- 1.96 SE with
  # SE = sqrt(n * var_u) / m_tid.
  # Not included: the uncertainty of the preset itself (or of a measured standard).
  z_crit <- 2.58                                    # two-sided 1%
  v_d    <- max(tid_mdl$var_d, tid_mdl$mean_d)
  v_u    <- max(tid_mdl$var_u, tid_mdl$mean_u)
  # Stoichiometric SE in the units of the reported multiplicity (a measured
  # standard rescales the mean; keep the model's variance-to-mean ratio).
  disp_u <- v_u / max(tid_mdl$mean_u, 1e-9)
  struct_ids <- sort(unique(mols$cluster2))
  summ <- lapply(struct_ids, function(sid) {
    idx <- which(mols$cluster2 == sid)
    nu  <- length(idx)
    mean_nn <- NA_real_
    if (nu >= 2L) {
      cc <- as.matrix(mols[idx, coord_cols, drop = FALSE])
      nn <- tryCatch(RANN::nn2(cc, cc, k = 2L)$nn.dists[, 2], error = function(e) NULL)
      if (!is.null(nn)) mean_nn <- stats::median(nn[is.finite(nn)])
    }
    sum_mol    <- sum(mol_est[idx])
    n_spatial  <- sum_mol / eta
    n_tid_tot  <- sum(mols$n_tid[idx])
    n_stoich   <- n_tid_tot / tids_per_mol
    tids_per_gt <- n_tid_tot / max(1L, nu)
    count_z    <- (n_tid_tot - sum_mol * tid_mdl$mean_d) / sqrt(max(1, sum_mol) * v_d)
    verdict    <- if (!is.finite(count_z)) "consistent"
                  else if (count_z >  z_crit) "under-resolved"
                  else if (count_z < -z_crit) "over-split"
                  else "consistent"
    mismatch   <- verdict != "consistent"
    fused_frac <- mean(mols$fused[idx])
    # Close or wide candidates suggest crowding, unless the trace count says the
    # structure is over-split: then the close neighbours are copies of one molecule.
    geom_crowd <- (is.finite(mean_nn) && mean_nn < max(footprint, link_nm + dye_spread)) ||
                  fused_frac > 0.3
    crowded <- verdict == "under-resolved" || (geom_crowd && verdict != "over-split")
    # Spatial count error = detection error: the detected count is
    # Binomial(N, eta), so SE(N_hat) = sqrt(n_det (1 - eta)) / eta. Resolution
    # error is not in this SE; it is what the consistency test reports.
    cnt_se  <- sqrt(max(1, sum_mol) * (1 - eta)) / eta
    if (crowded) cnt_se <- cnt_se * 1.5
    st_se   <- sqrt(max(1, n_stoich) * disp_u * tids_per_mol) / tids_per_mol
    # The range brackets BOTH estimators, each with its own sampling error.
    lo <- max(1, min(n_spatial - 1.96 * cnt_se, n_stoich - 1.96 * st_se))
    hi <- max(n_spatial + 1.96 * cnt_se, n_stoich + 1.96 * st_se)
    data.frame(
      cluster2      = sid,
      n_units       = nu,                # rendered GT groups (spatial)
      n_fused       = sum(mols$fused[idx]),
      fused_frac    = fused_frac,
      n_gt_est      = n_spatial,
      n_stoich      = n_stoich,
      n_tid_total   = n_tid_tot,
      tids_per_gt   = tids_per_gt,
      stoich_mismatch = mismatch,
      count_z       = count_z,
      count_verdict = verdict,
      n_count_fused = sum(mols$count_p[idx] < 0.01),
      n_gt_lo       = lo,
      n_gt_hi       = hi,
      eta           = eta,
      tids_per_mol  = tids_per_mol,
      mean_nn_nm    = mean_nn,
      crowded       = crowded,
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summ)

  meta <- list(
    level        = unit_level,
    n_fluor      = n_fluor, link_nm = link_nm, dye_spread_nm = dye_spread,
    footprint_nm = footprint, blink_lambda = blink_lambda, eff_pct = eff_pct,
    m_tid        = m_tid, eta = eta, sigma_loc_nm = sigma_loc_nm,
    tids_per_mol = tids_per_mol,
    m_tid_source = if (m_tid_measured) "measured" else "preset",
    r_single_nm  = r_single, complex_r_nm = complex_r,
    displace_nm  = meta_disp, guided_merge = isTRUE(guided_merge), use_z = use_z,
    preset_note = preset_note, scoped_to_cluster2 = isTRUE(scope_cluster2) && has_c2
  )
  list(points = points, links = links, summary = summary_df, meta = meta)
}

# Map per-candidate position uncertainty (nm) to a green->amber->red hex ramp.
# lo (best achievable, ~link_nm) maps to green; hi (poor) maps to red.
gt_uncertainty_colors <- function(pos_sd, lo, hi) {
  lo <- suppressWarnings(as.numeric(lo)); hi <- suppressWarnings(as.numeric(hi))
  if (!is.finite(lo)) lo <- 0
  if (!is.finite(hi) || hi <= lo) hi <- lo + 1
  tt   <- pmin(1, pmax(0, (pos_sd - lo) / (hi - lo)))
  tt[!is.finite(tt)] <- 1
  ramp <- grDevices::colorRamp(c("#2ecc71", "#f1c40f", "#e74c3c"))
  m    <- ramp(tt)
  grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
}

## ── Phase B: per-structure morphology / shape classification ──────────────────
# Contour length of a set of points via the DIAMETER of their Euclidean minimum
# spanning tree (longest path in the MST). Robust for CURVED filaments where a
# straight principal-axis extent would underestimate. O(n^2) Prim + tree BFS;
# fine for the tens-to-hundreds of backbone points per structure.
mst_contour_length <- function(P) {
  P <- as.matrix(P); n <- nrow(P)
  if (n < 2L) return(0)
  D <- as.matrix(stats::dist(P))
  inX <- logical(n); inX[1] <- TRUE
  best <- D[1, ]; bestfrom <- rep(1L, n); best[1] <- Inf
  from <- integer(n - 1L); to <- integer(n - 1L); w <- numeric(n - 1L)
  for (it in seq_len(n - 1L)) {
    j <- which.min(ifelse(inX, Inf, best))
    inX[j] <- TRUE
    from[it] <- bestfrom[j]; to[it] <- j; w[it] <- best[j]
    upd <- !inX & (D[j, ] < best)
    best[upd] <- D[j, upd]; bestfrom[upd] <- j
    best[inX] <- Inf
  }
  adj <- vector("list", n)
  for (e in seq_along(from)) {
    adj[[from[e]]] <- c(adj[[from[e]]], e)
    adj[[to[e]]]   <- c(adj[[to[e]]], e)
  }
  far <- function(src) {
    dist <- rep(Inf, n); dist[src] <- 0
    visited <- logical(n); visited[src] <- TRUE; queue <- src
    while (length(queue)) {
      u <- queue[1]; queue <- queue[-1]
      for (e in adj[[u]]) {
        v <- if (from[e] == u) to[e] else from[e]
        if (!visited[v]) { visited[v] <- TRUE; dist[v] <- dist[u] + w[e]; queue <- c(queue, v) }
      }
    }
    dist
  }
  d1 <- far(1L); a <- which.max(d1)
  d2 <- far(a); max(d2[is.finite(d2)])
}

# Collapse points onto a coarse grid (cell = h nm) and return cell centroids.
# Used to thin a backbone before the MST so lateral localization noise doesn't
# inflate the contour length with zig-zag, while curvature is preserved.
.grid_thin <- function(P, h) {
  if (h <= 0 || nrow(P) < 3L) return(P)
  key <- apply(round(P / h), 1L, paste, collapse = "_")
  cen <- lapply(split(seq_len(nrow(P)), key),
                function(ix) colMeans(P[ix, , drop = FALSE]))
  do.call(rbind, cen)
}

# Classify each cluster2's shape and report neutral geometric descriptors.
# EXPLORATORY. Shape from the covariance (gyration-tensor) eigenvalues; contour
# length from the backbone's minimum spanning tree. Not tied to any biology -
# the caller supplies an optional generic subunit spacing if a count is wanted.
#
#   df              clustered SMLM frame (x,y[,z], tid, cluster2)
#   backbone_pts    optional data.frame(x,y[,z],cluster2) of representative points
#                   (e.g. estimated GT candidates) used for the contour length and
#                   inter-point spacing; falls back to per-TID medians when NULL.
#   min_pts         structures with fewer representative points are skipped
#   subunit_spacing_nm  optional generic repeat spacing (nm); when finite and >0
#                   a filament's contour length is divided by it to give an
#                   estimated subunit count (no biological meaning imposed).
#
# Shape from the covariance eigenvalues l1>=l2(>=l3):
#   linearity  = (l1-l2)/l1   planarity = (l2-l3)/l1   sphericity = l3/l1
#   filament (linearity high), ring (hollow, ~constant radius), sheet (planar,
#   3D only), punctum (compact), else diffuse.
# Returns one row per cluster2; empty data.frame when nothing is analysable.
analyze_cluster2_morphology <- function(df, backbone_pts = NULL,
                                        min_pts = 4L, contour_smooth_nm = 10,
                                        subunit_spacing_nm = NA_real_) {
  if (is.null(df) || !nrow(df) || !all(c("x", "y", "tid") %in% names(df)))
    return(data.frame())
  if (!"cluster2" %in% names(df)) return(data.frame())
  use_z <- "z" %in% names(df) && any(is.finite(df$z) & df$z != 0)
  cc    <- if (use_z) c("x", "y", "z") else c("x", "y")

  work <- df[!is.na(df$cluster2) & df$cluster2 > 0L, , drop = FALSE]
  if (!nrow(work)) return(data.frame())
  if (!"z" %in% names(work)) work$z <- 0

  sids <- sort(unique(work$cluster2))
  out <- lapply(sids, function(sid) {
    sub <- work[work$cluster2 == sid, , drop = FALSE]
    # representative points for shape: per-TID medians (robust to loc noise)
    reps <- sub %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(x = stats::median(x), y = stats::median(y),
                       z = if (use_z) stats::median(z) else 0, .groups = "drop")
    P <- as.matrix(reps[, cc, drop = FALSE])
    if (nrow(P) < min_pts) return(NULL)

    ev  <- eigen(stats::cov(P), symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    l1  <- lam[1]; l2 <- if (length(lam) > 1) lam[2] else 0
    l3  <- if (length(lam) > 2) lam[3] else 0
    if (l1 <= 0) return(NULL)
    linearity  <- (l1 - l2) / l1
    planarity  <- (l2 - l3) / l1
    sphericity <- l3 / l1
    rg <- sqrt(sum(lam))

    # ring test: hollow centre + points at ~constant radius in the best-fit plane
    ctr <- colMeans(P)
    rad <- sqrt(rowSums(sweep(P, 2, ctr)^2))
    r_cv    <- if (mean(rad) > 0) stats::sd(rad) / mean(rad) else Inf
    hollow  <- mean(rad < 0.4 * stats::quantile(rad, 0.95)) < 0.12   # few points near centre
    is_ring <- hollow && r_cv < 0.35 && linearity < 0.5

    # contour + straight extent along principal axis
    axis1  <- ev$vectors[, 1]
    proj   <- as.numeric(P %*% axis1)
    extent <- diff(range(proj))                    # straight-line span
    width  <- 2 * sqrt(max(l2, 0))                 # transverse spread
    bpts <- if (!is.null(backbone_pts) && "cluster2" %in% names(backbone_pts)) {
      bp <- backbone_pts[backbone_pts$cluster2 == sid,
                         intersect(cc, names(backbone_pts)), drop = FALSE]
      if (nrow(bp) >= 2L) as.matrix(bp) else P
    } else P
    # thin lateral noise before measuring contour so it doesn't zig-zag-inflate
    thinned <- .grid_thin(bpts, contour_smooth_nm)
    contour <- mst_contour_length(thinned)
    contour <- max(contour, extent)                # never below the straight span

    # Shape class - dimension-aware. In 2D there is no sheet/sphericity axis, so
    # sphericity (l3/l1) is always 0; classify by linearity + ring only.
    shape <- if (is_ring) "ring"
             else if (linearity > 0.6) "filament"
             else if (use_z && planarity > 0.5 && sphericity < 0.25) "sheet"
             else if (use_z && sphericity > 0.45) "punctum"
             else if (!use_z && linearity < 0.35) "punctum"
             else "diffuse"

    # inter-point spacing among the backbone points (generic "distance between
    # GTs") and, for filaments, mean spacing along the contour + linear density.
    n_bb    <- nrow(bpts)
    nn_med  <- NA_real_
    if (n_bb >= 2L) {
      nn <- tryCatch(RANN::nn2(bpts, bpts, k = 2L)$nn.dists[, 2], error = function(e) NULL)
      if (!is.null(nn)) nn_med <- stats::median(nn[is.finite(nn)])
    }
    spacing_along <- if (shape == "filament" && n_bb >= 2L) contour / (n_bb - 1L) else NA_real_
    lin_density   <- if (shape == "filament" && contour > 0) n_bb / contour * 1000 else NA_real_  # per µm
    n_subunit <- if (shape == "filament" && is.finite(subunit_spacing_nm) &&
                     subunit_spacing_nm > 0) contour / subunit_spacing_nm else NA_real_
    data.frame(
      cluster2      = sid,
      shape         = shape,
      n_rep         = nrow(P),
      n_backbone    = n_bb,
      rg_nm         = rg,
      linearity     = linearity,
      planarity     = planarity,
      sphericity    = sphericity,
      extent_nm     = extent,
      contour_nm    = contour,
      width_nm      = width,
      ring_r_nm     = if (is_ring) stats::median(rad) else NA_real_,
      spacing_nn_nm = nn_med,          # median nearest-neighbour distance between points
      spacing_along_nm = spacing_along, # mean spacing along a filament backbone
      lin_density_um   = lin_density,   # points per µm of contour (filament)
      n_subunit_est    = n_subunit,     # only if a generic subunit spacing was given
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out)) data.frame() else out
}

