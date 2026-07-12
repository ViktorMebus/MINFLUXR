## ─────────────────────────────────────────────────────────────────────────────
## sim_server.R  —  Simulation tab server logic
## Synthetic MINFLUX localizations (shapes OR image-guided)
## ─────────────────────────────────────────────────────────────────────────────

sim_server <- function(input, output, session) {

  # ── Zoom-aware scale bar (mirrors the SMLM tab) ──────────────────────────────
  # The bar lives in paper coordinates with width = length_nm / visible_x_span,
  # so it must be re-pushed whenever the user zooms (the visible span changes).
  sim_axis_x_span <- reactiveVal(NULL)

  # TRUE once each plot has actually rendered WITH event_register(), so the
  # relayout observers never call event_data() before the source is registered.
  # That call would otherwise print a benign "source not registered" warning,
  # deferred via session$onFlushed() so suppressWarnings() can't catch it. They
  # stay FALSE until the user clicks Generate (the render short-circuits on a
  # validate() until then, so the source is never registered before that).
  sim_plot_ready <- reactiveVal(FALSE)
  sim_diff_ready <- reactiveVal(FALSE)

  # Single source of truth: is the FOV / ROI window active?  TRUE in Acquisition
  # mode, or in Normal mode when the user pins an ROI ("Pin axes to a FOV").
  # Used everywhere the FOV matters — axis pinning, background spread, and the
  # confocal window — so they can never disagree.
  sim_use_fov <- function() {
    identical(input$sim_mode %||% "normal", "acquisition") ||
      isTRUE(input$sim_view_fov)
  }

  # Full-view x-span: the FOV X when the FOV is active, else the data extent.
  sim_full_x_span <- function() {
    if (isolate(sim_use_fov()))
      return(max(1, as.numeric(isolate(input$sim_fov_x) %||% 4000)))
    res <- tryCatch(isolate(sim_result()), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    ref  <- if (nrow(res$locs) > 0) res$locs else res$gt
    span <- suppressWarnings(diff(range(ref$x, na.rm = TRUE)))
    if (length(span) == 1 && is.finite(span) && span > 0) span else NULL
  }

  sim_push_scale_bar <- function() {
    proxy <- plotlyProxy("sim_plot", session)
    if (!isTRUE(input$sim_show_scale_bar) ||
        !identical(input$sim_view_mode, "2d")) {
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = list(), annotations = list()))
      return()
    }
    span <- sim_axis_x_span()
    if (length(span) != 1 || !is.finite(span) || span <= 0) {
      span <- sim_full_x_span()
      if (length(span) == 1 && is.finite(span) && span > 0) sim_axis_x_span(span)
    }
    if (length(span) != 1 || !is.finite(span) || span <= 0) return()
    sb_font <- suppressWarnings(as.numeric(input$sim_scale_font_size))
    sb <- make_scale_bar(
      suppressWarnings(as.numeric(input$sim_scale_bar_length)),
      span,
      suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
      dark = isTRUE(input$sim_dark_mode),
      scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
    if (!is.null(sb))
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = sb$shapes, annotations = sb$annotations))
  }

  # Same zoom-aware scale bar for the Confocal view (always white on black).
  sim_diff_x_span <- reactiveVal(NULL)   # zoomed visible x-span
  sim_diff_full   <- reactiveVal(NULL)   # full rendered-image x-span

  sim_push_diff_scale_bar <- function() {
    proxy <- plotlyProxy("sim_diff_plot", session)
    if (!isTRUE(input$sim_show_scale_bar)) {
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = list(), annotations = list()))
      return()
    }
    span <- sim_diff_x_span()
    if (length(span) != 1 || !is.finite(span) || span <= 0) span <- sim_diff_full()
    if (length(span) != 1 || !is.finite(span) || span <= 0) return()
    sb_font <- suppressWarnings(as.numeric(input$sim_scale_font_size))
    sb <- make_scale_bar(
      suppressWarnings(as.numeric(input$sim_scale_bar_length)),
      span,
      suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
      dark = TRUE,
      scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
    if (!is.null(sb))
      plotlyProxyInvoke(proxy, "relayout",
                        list(shapes = sb$shapes, annotations = sb$annotations))
  }

  # Front-loaded detection times: the rate decays from a peak toward a constant
  # floor, so most trace IDs are detected early and the rate flattens later.
  # Inverse-CDF sampling of n times on [0, T] from density
  #   f(t) = floor_frac + (1 - floor_frac) · exp(-t / tau).
  sim_sample_detection_times <- function(n, T_sec, tau_sec, floor_frac) {
    if (n <= 0) return(numeric(0))
    if (!is.finite(T_sec) || T_sec <= 0) return(rep(0, n))
    ng   <- 2000L
    tg   <- seq(0, T_sec, length.out = ng)
    dens <- pmax(1e-9,
                 floor_frac + (1 - floor_frac) * exp(-tg / max(1e-6, tau_sec)))
    cdf  <- cumsum(dens); cdf <- cdf / cdf[ng]
    approx(cdf, tg, runif(n), rule = 2, ties = "ordered")$y
  }

  # Expected number of structure detections over [0, T] for an ABSOLUTE rate
  # r(t) = struct_rate · (floor + (1-floor)·e^{-t/τ})  [traces/min]:
  #   N_T = struct_rate · ( floor·T + (1-floor)·τ·(1 - e^{-T/τ}) )   (T, τ in min)
  # The early (exp) term is fixed, so the per-elapsed-time count is independent
  # of the total duration — only the floor·T term keeps growing with T.
  acq_expected_count <- function(struct_rate, T_sec, decay_min, floor_frac) {
    T_min   <- T_sec / 60
    tau_min <- max(1e-6, decay_min)
    struct_rate * (floor_frac * T_min +
                   (1 - floor_frac) * tau_min * (1 - exp(-T_min / tau_min)))
  }

  # Acquisition timeline from a fixed fluorophore pool + live acquisition
  # parameters.  Operates only on molecule IDs (no positions/localizations), so
  # it is cheap enough to recompute live as rate / duration / background change
  # without re-running the whole structure.  pool_mol = mol_id of every unique
  # fluorophore (length = pool size); n_mol_total = molecules placed (for
  # display).  Returns the same shape of list the acquisition plot/summary use.
  compute_acq_timeline <- function(pool_mol, n_mol_total, struct_rate, T_sec,
                                   blinks, decay_min, floor_frac, bg_rate,
                                   resolve_q, crowd_m = NULL) {
    n_pool <- length(pool_mol)
    if (n_pool == 0 || !is.finite(T_sec) || T_sec <= 0 || struct_rate <= 0)
      return(NULL)
    n_mol_resolvable <- length(unique(pool_mol))

    # Deterministic EXPECTED timeline (no per-realization jitter), so the curves
    # and the time-to-resolve are stable as the user drags rate / duration.
    # struct_rate / bg_rate are TRACES PER HOUR (easier to match real settings);
    # converted to per-minute internally.  Detections follow an ABSOLUTE rate
    # r(t) = X·(floor + (1-floor)·e^{-t/τ}), so how fast the structure resolves
    # depends on the rate and structure — NOT on the chosen duration.
    T_min     <- T_sec / 60
    tau_min   <- max(1e-6, decay_min)
    rate_min  <- struct_rate / 60                    # per hour → per minute
    bg_min    <- bg_rate / 60

    # Photon budget: each fluorophore yields ~blinks detections before bleaching.
    # This caps the TOTAL TIDs (over-counting) only — it does NOT gate resolution
    # (a molecule is resolved by its FIRST detection; the extra blinks just pile
    # more TIDs on top).
    struct_cap <- n_pool * max(1, blinks)
    f_m        <- as.numeric(tapply(rep(1, n_pool), pool_mol, sum))  # dyes / molecule

    ng     <- 240L
    tg_min <- seq(0, T_min, length.out = ng)
    tg_sec <- tg_min * 60
    # Cumulative fluorophore-detection attempts from the absolute rate.
    cum_R <- rate_min * (floor_frac * tg_min +
               (1 - floor_frac) * tau_min * (1 - exp(-tg_min / tau_min)))
    # Photobleaching consumes the photon budget GRADUALLY (fluorophores bleach at
    # different times), so the total saturates smoothly toward struct_cap and the
    # structure rate tapers off — no hard cliff when the budget is reached.
    bleach <- exp(-cum_R / struct_cap)          # fraction of budget still emitting
    N_t    <- struct_cap * (1 - bleach)         # total recorded TIDs (over-counting)
    n_s    <- N_t[ng]

    # Expected distinct molecules resolved by time t: a molecule is found once
    # ANY of its f_m fluorophores is first detected.  Driven by cum_R (rate ×
    # time × pool) — INDEPENDENT of the per-fluorophore blink budget, so
    # "detections before bleaching" never shifts the road to X% resolved.
    # crowd_m (≤ 1, per molecule, aligned to the sorted unique mol_ids) slows
    # discovery where emitters crowd the donut: a crowded molecule's effective
    # detection rate is scaled down, so dense structures resolve slower.
    cw <- if (is.null(crowd_m) || length(crowd_m) != length(f_m))
            rep(1, length(f_m)) else crowd_m
    cum_unique <- vapply(cum_R,
                         function(c) sum(1 - exp(-f_m * cw * c / n_pool)),
                         numeric(1))
    n_mol_det  <- cum_unique[ng]
    denom      <- max(1, n_mol_resolvable)

    # Instantaneous detection rate (traces/HOUR): structure (decaying with τ and
    # tapering as the dyes bleach) + constant background floor.
    rate_struct <- struct_rate *
                   (floor_frac + (1 - floor_frac) * exp(-tg_min / tau_min)) * bleach
    rate_v <- rate_struct + bg_rate
    n_bg   <- max(0L, min(200000L, as.integer(round(bg_min * T_min))))

    q        <- max(0.5, min(1, resolve_q))
    res_time <- if (n_mol_det >= q * denom)
                  approx(cum_unique, tg_sec, q * denom, ties = "ordered")$y
                else NA_real_

    list(
      T_sec = T_sec, struct_rate = struct_rate, bg_rate = bg_rate,
      n_pool = n_pool, n_mol_total = n_mol_total,
      n_mol_resolvable = n_mol_resolvable, struct_cap = struct_cap, n_bg = n_bg,
      t = tg_sec, cum_unique = cum_unique, rate_t = tg_sec, rate_v = rate_v,
      rate_struct = rate_struct,          # structure-only component of rate_v
      rate_peak  = max(rate_v),
      rate_floor = stats::median(tail(rate_v, max(1L, length(rate_v) %/% 4L))),
      rate_mean  = n_s / max(1e-9, T_min) * 60,      # traces/hour
      n_mol_det = round(n_mol_det), resolved_frac = n_mol_det / denom,
      n_real_tid = round(n_s), det_per_mol = n_s / max(1, n_mol_det),
      resolve_q = q, res_time = res_time,
      overcount_frac = 1 - n_mol_det / max(1, n_s))
  }

  # Calibrate the acquisition model from a REAL MINFLUX run: measure the
  # detection-rate profile (peak / mean / floor / decay τ), the run duration and
  # the FOV extent from the data.  Works on the first detection time of each
  # unique trace (tid).  df = canonical MINFLUX frame (tid, tim [s], x, y[, z] nm).
  analyze_real_run <- function(df) {
    if (!all(c("tid", "tim", "x", "y") %in% names(df))) return(NULL)
    df <- df[is.finite(df$tim) & is.finite(df$x) & is.finite(df$y), , drop = FALSE]
    if (nrow(df) < 10) return(NULL)

    # First detection time + median position of each trace (same tid order).
    t0_tid <- as.numeric(tapply(df$tim, df$tid, min))
    tx     <- as.numeric(tapply(df$x,   df$tid, median))
    ty     <- as.numeric(tapply(df$y,   df$tid, median))
    t0_tid <- t0_tid - min(t0_tid)                    # seconds from start
    n_tr   <- length(t0_tid)
    T_sec  <- max(t0_tid)
    if (!is.finite(T_sec) || T_sec <= 0) return(NULL)

    # Binned detection rate (traces / hour).
    nb     <- max(8L, min(60L, n_tr %/% 10L))
    br     <- seq(0, T_sec, length.out = nb + 1L)
    bin_hr <- max(1e-9, (br[2] - br[1]) / 3600)
    h      <- hist(t0_tid, breaks = br, plot = FALSE)
    rate_t <- h$mids                                  # seconds
    rate_v <- h$counts / bin_hr

    R0     <- mean(head(rate_v, max(1L, nb %/% 10L)))           # early rate
    Rf     <- stats::median(tail(rate_v, max(1L, nb %/% 4L)))   # late rate

    R_mean <- n_tr / (T_sec / 3600)

    # τ from a log-linear fit of |rate − late_rate| toward 0 — works whether the
    # rate falls (dye bleaching) OR rises (e.g. buffer activating over the run).
    tau_min <- 15
    exc     <- abs(rate_v - Rf)
    ok      <- exc > 0.05 * max(exc, 1e-9) & rate_t < T_sec * 0.8
    if (sum(ok) >= 3) {
      fit <- tryCatch(stats::lm(log(exc[ok]) ~ rate_t[ok]), error = function(e) NULL)
      if (!is.null(fit)) {
        slope <- stats::coef(fit)[[2]]
        if (is.finite(slope) && slope < 0) tau_min <- max(0.5, (-1 / slope) / 60)
      }
    }
    # If the rate visibly changed over the run, τ can't exceed it — a noisy
    # (especially rising) profile can fit an absurd τ ≫ run, which would make the
    # simulated change look flat.  Cap at the run length.
    tau_min <- min(tau_min, T_sec / 60)

    fov_x <- diff(range(df$x, na.rm = TRUE))
    fov_y <- diff(range(df$y, na.rm = TRUE))
    fov_z <- if ("z" %in% names(df)) diff(range(df$z, na.rm = TRUE)) else 0

    # Spatial split: a trace sitting in a locally dense neighbourhood is ON the
    # structure; an isolated trace is background.  Uses the DBSCAN core-point
    # criterion (≥ minPts neighbours within eps) — absolute density, so it works
    # even when the structure itself is densely labeled.  Approximate; the user
    # can override the resulting fraction.  Also bins the two classes separately
    # so the measured structure / background rate curves can be overlaid.
    struct_frac   <- 0.5
    rate_v_struct <- NULL
    rate_v_bg     <- NULL
    if (n_tr >= 20) {
      eps    <- 60          # nm neighbourhood radius
      minPts <- 3L          # core-point threshold
      K      <- min(n_tr, 64L)
      nn     <- RANN::nn2(cbind(tx, ty), cbind(tx, ty), k = K,
                          searchtype = "radius", radius = eps)
      cnt    <- rowSums(nn$nn.idx > 0) - 1L          # neighbours within eps
      is_struct   <- cnt >= minPts
      struct_frac <- max(0, min(1, mean(is_struct)))
      rate_v_struct <- hist(t0_tid[is_struct],  breaks = br, plot = FALSE)$counts / bin_hr
      rate_v_bg     <- hist(t0_tid[!is_struct], breaks = br, plot = FALSE)$counts / bin_hr
    }

    # Localization precision = SD of the repeated final localizations WITHIN a
    # trace (the molecule re-localized while ON); median over traces with ≥3.
    # Photons/loc from the effective count (eco).  Effective TCP diameter that
    # reproduces this precision at this photon count: L_eff = 2·σ·√N.
    nl  <- as.numeric(tapply(df$x, df$tid, length))
    sx  <- as.numeric(tapply(df$x, df$tid, stats::sd))
    sy  <- as.numeric(tapply(df$y, df$tid, stats::sd))
    okp <- nl >= 3 & is.finite(sx) & is.finite(sy)
    sigma_x <- if (any(okp)) stats::median(sx[okp]) else NA_real_
    sigma_y <- if (any(okp)) stats::median(sy[okp]) else NA_real_
    sigma_z <- NA_real_
    if ("z" %in% names(df)) {
      sz <- as.numeric(tapply(df$z, df$tid, stats::sd))
      if (any(okp & is.finite(sz))) sigma_z <- stats::median(sz[okp & is.finite(sz)])
    }
    photons  <- if ("eco" %in% names(df)) stats::median(df$eco, na.rm = TRUE) else NA_real_
    sigma_xy <- mean(c(sigma_x, sigma_y), na.rm = TRUE)
    L_eff    <- if (is.finite(sigma_xy) && is.finite(photons) && photons > 0)
                  2 * sigma_xy * sqrt(photons) else NA_real_

    list(n_traces = n_tr, duration_hr = T_sec / 3600,
         rate_peak = R0, rate_mean = R_mean, rate_floor = Rf,
         tau_min = tau_min, struct_frac = struct_frac,
         fov_x = fov_x, fov_y = fov_y, fov_z = fov_z,
         rate_t = rate_t, rate_v = rate_v,
         rate_v_struct = rate_v_struct, rate_v_bg = rate_v_bg,
         sigma_x = sigma_x, sigma_y = sigma_y, sigma_z = sigma_z,
         photons = photons, L_eff = L_eff)
  }

  # ╔══════════════════════════════════════════════════════════════════════════╗
  # ║  PHASE 1 — Synthetic localizations                                      ║
  # ╚══════════════════════════════════════════════════════════════════════════╝

  # ── Shape generators ────────────────────────────────────────────────────────

  gen_cluster <- function(n, radius, distribution = "uniform", use_3d = FALSE) {
    sigma <- radius / 2
    if (use_3d) {
      if (distribution == "gaussian") {
        df <- data.frame(x = rnorm(n, 0, sigma),
                         y = rnorm(n, 0, sigma),
                         z = rnorm(n, 0, sigma))
      } else {
        dirs <- matrix(rnorm(n * 3), ncol = 3)
        lens <- sqrt(rowSums(dirs^2)); lens[lens == 0] <- 1
        dirs <- dirs / lens
        r3   <- radius * runif(n)^(1/3)
        df <- data.frame(x = r3 * dirs[,1], y = r3 * dirs[,2], z = r3 * dirs[,3])
      }
    } else {
      r     <- if (distribution == "gaussian") abs(rnorm(n, 0, sigma))
               else radius * sqrt(runif(n))
      theta <- runif(n, 0, 2 * pi)
      df <- data.frame(x = r * cos(theta), y = r * sin(theta))
    }
    df$ldir_type <- NA_character_
    df
  }

  gen_ring <- function(n, diameter, ring_width, z_thickness = 0,
                        use_3d = FALSE) {
    radius <- diameter / 2
    theta  <- runif(n, 0, 2 * pi)
    r      <- rnorm(n, mean = radius, sd = max(0.1, ring_width / 4))
    df     <- data.frame(x = r * cos(theta), y = r * sin(theta))
    if (use_3d) df$z <- rnorm(n, 0, max(0.1, z_thickness / 4))
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  gen_ellipse <- function(n, semi_a, semi_b, semi_c = 1,
                           filled = TRUE, use_3d = FALSE) {
    if (use_3d) {
      phi   <- acos(2 * runif(n) - 1)
      theta <- runif(n, 0, 2 * pi)
      r     <- if (filled) runif(n)^(1/3) else rep(1, n)
      df <- data.frame(x = semi_a * r * sin(phi) * cos(theta),
                       y = semi_b * r * sin(phi) * sin(theta),
                       z = semi_c * r * cos(phi))
    } else {
      theta <- runif(n, 0, 2 * pi)
      if (filled) {
        r  <- sqrt(runif(n))
        df <- data.frame(x = semi_a * r * cos(theta), y = semi_b * r * sin(theta))
      } else {
        df <- data.frame(x = semi_a * cos(theta), y = semi_b * sin(theta))
      }
    }
    df$ldir_type <- NA_character_
    df
  }

  # n molecules spaced at spacing_nm along a straight filament.
  # jitter_lat_nm: noise perpendicular to the filament axis (isotropic in perp plane).
  # jitter_ax_nm:  noise along the filament axis (arc-length perturbation).
  # helix_r / helix_pitch: optional helical arrangement (0 = disabled);
  #   default pitch 36 nm mirrors B-DNA major-groove repeat.
  gen_filament <- function(n, spacing_nm, angle_deg = 0,
                            jitter_lat_nm = 0, jitter_ax_nm = 0,
                            helix_r = 0, helix_pitch = 36,
                            use_3d = FALSE, z_angle_deg = 0) {
    n      <- max(2L, n)
    s      <- (seq_len(n) - (n + 1) / 2) * spacing_nm   # centred arc-lengths
    if (jitter_ax_nm > 0) s <- s + rnorm(n, 0, jitter_ax_nm)

    xy_rad <- angle_deg  * pi / 180
    z_rad  <- z_angle_deg * pi / 180

    # Filament axis unit vector
    ax <- cos(xy_rad) * cos(z_rad)
    ay <- sin(xy_rad) * cos(z_rad)
    az <- sin(z_rad)

    # Two perpendicular unit vectors (analytical, both already unit length):
    # p1 = lateral perpendicular in XY plane; p2 = cross(axis, p1) → toward Z
    p1x <- -sin(xy_rad); p1y <- cos(xy_rad); p1z <- 0
    p2x <- -sin(z_rad) * cos(xy_rad)
    p2y <- -sin(z_rad) * sin(xy_rad)
    p2z <-  cos(z_rad)

    # Lateral jitter as a random walk along the filament: each molecule's
    # perpendicular offset = previous offset + small increment.  This makes
    # the filament curve and meander rather than scatter independently.
    # jitter_lat_nm is the per-step SD (nm); larger values = more curling.
    if (jitter_lat_nm > 0) {
      jit1 <- cumsum(rnorm(n, 0, jitter_lat_nm))
      jit2 <- if (use_3d) cumsum(rnorm(n, 0, jitter_lat_nm)) else rep(0, n)
    } else {
      jit1 <- rep(0, n)
      jit2 <- rep(0, n)
    }

    # Helical displacement: cosine in p1, sine in p2 → right-handed helix
    if (helix_r > 0 && helix_pitch > 0) {
      phi  <- 2 * pi * s / helix_pitch
      hel1 <- helix_r * cos(phi)
      hel2 <- helix_r * sin(phi)
    } else {
      hel1 <- rep(0, n)
      hel2 <- rep(0, n)
    }

    d1 <- jit1 + hel1   # total displacement in p1 direction
    d2 <- jit2 + hel2   # total displacement in p2 direction (0 in 2D XY)

    df <- data.frame(
      x = s * ax + d1 * p1x + d2 * p2x,
      y = s * ay + d1 * p1y + d2 * p2y
    )
    if (use_3d) df$z <- s * az + d1 * p1z + d2 * p2z
    # Store filament axis tangent for direction-aware linker displacement.
    # All molecules share the same global axis direction (straight filament).
    df$ldir_x    <- ax
    df$ldir_y    <- ay
    df$ldir_z    <- if (use_3d) az else 0
    df$ldir_type <- "perp"
    df
  }

  # Hollow sphere — uniform surface sampling with adjustable shell thickness.
  # In 2D collapses to a ring (circle cross-section).
  gen_hollow_sphere <- function(n, radius, shell_thick, use_3d = FALSE) {
    if (!use_3d) {
      r     <- rnorm(n, mean = radius, sd = max(0.1, shell_thick / 4))
      theta <- runif(n, 0, 2 * pi)
      df <- data.frame(x = r * cos(theta), y = r * sin(theta))
    } else {
      phi   <- acos(2 * runif(n) - 1)
      theta <- runif(n, 0, 2 * pi)
      r     <- radius + runif(n, -shell_thick / 2, shell_thick / 2)
      df <- data.frame(x = r * sin(phi) * cos(theta),
                       y = r * sin(phi) * sin(theta),
                       z = r * cos(phi))
    }
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  # Double ring — two parallel rings separated along Z (or offset along Y in 2D).
  # Useful for nuclear pore complex (~35 nm separation), centrioles, etc.
  gen_double_ring <- function(n, diameter, ring_width, separation_nm,
                               z_thickness = 0, use_3d = FALSE) {
    radius <- diameter / 2
    n1 <- round(n / 2); n2 <- n - n1
    make_one <- function(n_pts, z_pos) {
      theta <- runif(n_pts, 0, 2 * pi)
      r     <- rnorm(n_pts, mean = radius, sd = max(0.1, ring_width / 4))
      df    <- data.frame(x = r * cos(theta), y = r * sin(theta))
      if (use_3d)
        df$z <- rnorm(n_pts, z_pos, max(0.1, z_thickness / 4))
      df
    }
    df <- rbind(make_one(n1, -separation_nm / 2),
                make_one(n2,  separation_nm / 2))
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  # Arc — a partial ring spanning arc_extent_deg degrees.
  # arc_start_deg sets the angular position of the start of the arc.
  gen_arc <- function(n, diameter, ring_width, arc_start_deg, arc_extent_deg,
                       z_thickness = 0, use_3d = FALSE) {
    radius     <- diameter / 2
    start_rad  <- arc_start_deg  * pi / 180
    extent_rad <- max(1, arc_extent_deg) * pi / 180
    theta <- runif(n, start_rad, start_rad + extent_rad)
    r     <- rnorm(n, mean = radius, sd = max(0.1, ring_width / 4))
    df    <- data.frame(x = r * cos(theta), y = r * sin(theta))
    if (use_3d)
      df$z <- rnorm(n, 0, max(0.1, z_thickness / 4))
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    df
  }

  # Oligomeric ring stack — fixed, evenly-spaced subunit positions (no randomness).
  # Models protein complexes: hexamers, heptamers,
  # nuclear pore subcomplexes, etc.
  # GT positions are rigidly determined: exactly n_subunits × n_stacks molecules.
  # angle_offset_deg: each successive ring is rotated by this amount (e.g., for
  # symmetry mismatch between stacked assemblies).
  # In 3D, rings are separated along Z by separation_nm.
  # In 2D, all rings project onto XY (stacks overlap); use angle_offset_deg to
  # make the rotational offset visible.
  gen_oligoring_stack <- function(n_subunits, diameter, n_stacks,
                                   separation_nm, angle_offset_deg = 0,
                                   base_angle_deg = 0,
                                   use_3d = FALSE) {
    n_subunits <- max(2L, as.integer(n_subunits))
    n_stacks   <- max(1L, as.integer(n_stacks))
    radius     <- diameter / 2
    off_rad    <- angle_offset_deg * pi / 180
    base_rad   <- base_angle_deg   * pi / 180

    layers <- lapply(seq_len(n_stacks), function(k) {
      z_pos  <- (k - (n_stacks + 1) / 2) * separation_nm
      rot    <- base_rad + (k - 1L) * off_rad
      angles <- rot + 2 * pi * (seq_len(n_subunits) - 1L) / n_subunits
      df     <- data.frame(x = radius * cos(angles),
                           y = radius * sin(angles))
      if (use_3d) df$z <- z_pos
      df
    })
    df <- do.call(rbind, layers)
    df$ldir_type <- "radial"
    df$struct_cx <- 0; df$struct_cy <- 0
    if (use_3d) df$struct_cz <- 0
    rownames(df) <- NULL
    df
  }

  # ── Image-guided helpers ─────────────────────────────────────────────────────

  read_image_gray <- function(path) {
    img   <- magick::image_read(path)
    img_g <- magick::image_convert(img, colorspace = "Gray")
    arr   <- magick::image_data(img_g, "gray")
    # dim(arr) = [1, width, height]  →  [H x W] matrix
    W <- dim(arr)[2]; H <- dim(arr)[3]
    mat <- t(matrix(as.integer(arr[1,,]), nrow = W, ncol = H)) / 255
    rng <- range(mat)
    if (diff(rng) > .Machine$double.eps)
      mat <- (mat - rng[1]) / diff(rng)
    list(mat = mat, img = img_g)
  }

  img_to_b64 <- function(img_magick, max_px = 800L) {
    info <- magick::image_info(img_magick)
    if (max(info$width, info$height) > max_px)
      img_magick <- magick::image_scale(
        img_magick, paste0(max_px, "x", max_px))
    raw_bytes <- magick::image_write(img_magick, format = "png")
    paste0("data:image/png;base64,", base64enc::base64encode(raw_bytes))
  }

  # Sample proteins from confocal image intensity.
  # XY is always drawn from pixel intensity (proportional to brightness).
  # In 3D mode, Z spread scales with local pixel intensity: bright regions
  # receive more axial depth (z_spread_nm), dim regions stay near Z = 0.
  gen_from_image <- function(img_gray, n_proteins, px_size_nm,
                              threshold, min_trace, max_trace,
                              use_3d = FALSE,
                              z_spread_nm = 200, z_bg_nm = 20) {
    m     <- img_gray; m[m < threshold] <- 0
    probs <- as.vector(m)
    if (sum(probs) == 0)
      stop("All pixels are below the threshold — lower the threshold value.")
    probs <- probs / sum(probs)
    H <- nrow(m); W <- ncol(m)

    idx     <- sample(H * W, n_proteins, replace = TRUE, prob = probs)
    row_idx <- ((idx - 1L) %% H) + 1L
    col_idx <- ((idx - 1L) %/% H) + 1L
    x_nm    <- (col_idx - 0.5 + runif(n_proteins, -0.5, 0.5)) * px_size_nm -
               (W / 2) * px_size_nm
    y_nm    <- ((H - row_idx + 0.5) + runif(n_proteins, -0.5, 0.5)) * px_size_nm -
               (H / 2) * px_size_nm
    local_int  <- probs[idx] / max(probs, na.rm = TRUE)   # 0–1
    n_locs_exp <- min_trace + (max_trace - min_trace) * local_int

    if (!use_3d) {
      df <- data.frame(x = x_nm, y = y_nm, n_locs_expected = n_locs_exp)
    } else {
      # 3D: Z spread = pmax(z_bg_nm, z_spread_nm × local_intensity)
      z_sd  <- pmax(max(0.5, z_bg_nm), z_spread_nm * local_int)
      z_nm  <- rnorm(n_proteins, 0, z_sd)
      z_nm  <- pmax(-1000, pmin(1000, z_nm))
      df <- data.frame(x = x_nm, y = y_nm, z = z_nm, n_locs_expected = n_locs_exp)
    }
    df$ldir_type <- NA_character_
    df
  }

  # ── Labeling presets ─────────────────────────────────────────────────────────
  # n_fluor: mean fluorophores per labeled molecule (Poisson-drawn, min 1).
  # link_nm: mean distance fluorophore→true molecule centre (nm).
  # In MINFLUX each blinking fluorophore = one trace ID.
  # n_fluor is the *theoretical* count; actual detected TIDs are ~60–80 % of
  # this due to incomplete blinking/detection efficiency.
  #
  # Key assumptions:
  #   IgG arm (Fab tip to Fc base) ≈ 11 nm; secondary adds another ≈ 11 nm.
  #   Nanobody body ≈ 2.5 nm → effective linker ~3 nm.
  #   NHS-ester DOL for commercial IgG kits: 2–4 (mean ~3).
  #   Monoclonal: 1 antibody per epitope (bivalent, but typically 1 per target).
  #   Polyclonal: 3–5 antibodies per antigen (multiple epitopes accessible).
  #   Secondary binding stoichiometry: 2–4 secondaries per primary IgG.
  # dye_spread_nm: Gaussian sigma of individual dye positions within the label
  # complex (around the shared complex-centre displaced by link_nm).
  # n_fluor > 1 on the SAME complex → dyes cluster within dye_spread_nm;
  # each molecule's complex is displaced independently by link_nm.
  LABEL_PRESETS <- list(

    # ── Covalent / minimal-linkage tags ─────────────────────────────────────
    # footprint_nm:  approximate excluded-zone diameter around each labeled protein.
    #                If mean protein spacing < footprint_nm, steric exclusion is likely.
    # blink_lambda:  mean extra TID appearances per fluorophore from dye re-activation.
    "Direct tag (SNAP / Halo / Click, DOL=1)" = list(
      n_fluor = 1L, link_nm = 3.5, dye_spread_nm = 0,
      footprint_nm = 2,  blink_lambda = 0.2,
      note    = "1 fluorophore / molecule \u2014 ~3.5 nm linkage (SNAP/Halo tag domain)"
    ),
    "Nanobody (DOL \u22481)" = list(
      n_fluor = 1L, link_nm = 2.5, dye_spread_nm = 0,
      footprint_nm = 3,  blink_lambda = 0.2,
      note    = "1 fluorophore / molecule \u2014 ~2.5 nm (nanobody VHH domain)"
    ),

    # ── Primary monoclonal (1 antibody / epitope) ────────────────────────────
    "Primary monoclonal IgG (DOL \u22481)" = list(
      n_fluor = 1L, link_nm = 8, dye_spread_nm = 0,
      footprint_nm = 10, blink_lambda = 0.3,
      note    = "1 fluorophore / molecule \u2014 ~8 nm (IgG Fab arm to epitope)"
    ),
    "Primary monoclonal IgG (NHS, DOL \u22482\u20134)" = list(
      n_fluor = 3L, link_nm = 8, dye_spread_nm = 4,
      footprint_nm = 10, blink_lambda = 0.5,
      note    = "~3 fluorophores on one IgG body (~4 nm spread) \u2014 ~8 nm from target"
    ),

    # ── Primary polyclonal (3–5 antibodies / antigen) ────────────────────────
    "Primary polyclonal IgG (NHS, DOL \u22482\u20134)" = list(
      n_fluor = 9L, link_nm = 8, dye_spread_nm = 9,
      footprint_nm = 12, blink_lambda = 0.5,
      note    = "~9 fluorophores across 3\u20135 Ab bodies (~9 nm spread) \u2014 ~8 nm from target"
    ),

    # ── Monoclonal + secondary ───────────────────────────────────────────────
    "Primary monoclonal IgG + secondary nanobody" = list(
      n_fluor = 2L, link_nm = 10.5, dye_spread_nm = 3,
      footprint_nm = 10, blink_lambda = 0.4,
      note    = "~2 fluorophores, nanobody on primary IgG \u2014 ~10.5 nm from target"
    ),
    "Primary monoclonal IgG + secondary IgG" = list(
      n_fluor = 8L, link_nm = 20, dye_spread_nm = 12,
      footprint_nm = 16, blink_lambda = 0.6,
      note    = "~8 fluorophores across ~3 secondary IgG bodies (~12 nm spread) \u2014 ~20 nm from target"
    ),

    # ── Polyclonal + secondary ───────────────────────────────────────────────
    "Primary polyclonal IgG + secondary IgG" = list(
      n_fluor = 18L, link_nm = 20, dye_spread_nm = 18,
      footprint_nm = 20, blink_lambda = 0.8,
      note    = "~18 fluorophores, large complex (~18 nm spread) \u2014 ~20 nm from target"
    ),

    # ── Biotin–streptavidin ──────────────────────────────────────────────────
    "Streptavidin (biotinylated Ab + SA-dye)" = list(
      n_fluor = 3L, link_nm = 15, dye_spread_nm = 3,
      footprint_nm = 8,  blink_lambda = 0.4,
      note    = "~3 fluorophores on SA tetramer (~3 nm spread) \u2014 ~15 nm from target"
    ),

    "Custom" = list(
      n_fluor = NA_integer_, link_nm = NA_real_,
      dye_spread_nm = NA_real_, footprint_nm = NA_real_,
      blink_lambda = NA_real_, note = ""
    )
  )

  # Expand labeled molecules to individual fluorophore positions.
  # eff_pct:       labeling efficiency (% of molecules that get labeled)
  # n_fluor:       mean fluorophores per labeled molecule (Poisson-drawn, min 1)
  # link_nm:       mean distance from label-complex centre to molecule centre (nm)
  # dye_spread_nm: Gaussian SD of individual dye positions within the complex
  #
  # Two-level displacement model:
  #   Level 1 (per molecule): one complex-centre offset of ~link_nm, shared by
  #     all fluorophores on that molecule.
  #   Level 2 (per fluorophore): each dye scatters by dye_spread_nm around the
  #     complex centre. For single-dye labels dye_spread_nm = 0.
  # Each returned row = one fluorophore = one future trace ID.
  expand_to_fluorophores <- function(gt_df, eff_pct, n_fluor,
                                      link_nm, dye_spread_nm = 0,
                                      ldir_aware = TRUE,
                                      protein_r_nm = 0) {
    use_3d <- "z" %in% names(gt_df)

    # Step 1: all GT molecules can be labeled; orig_rows tracks indices into
    # the full gt_df so mol_id can look up res$gt in the plotting code.
    orig_rows <- seq_len(nrow(gt_df))
    n_mol     <- nrow(gt_df)
    if (n_mol == 0) {
      base <- gt_df[, intersect(c("x","y","z","n_locs_expected"), names(gt_df)),
                    drop = FALSE]
      base$mol_id <- integer(0)
      return(base)
    }

    # Step 2: draw fluorophore count per molecule (Poisson, min 1)
    if (n_fluor <= 1L) {
      n_per <- rep(1L, n_mol)
    } else {
      n_per <- pmax(1L, rpois(n_mol, n_fluor))
    }

    # Step 3: expand molecule rows to fluorophore rows.
    # mol_id stores the original GT row index for use in the plotting code.
    mol_idx <- rep(seq_len(n_mol), n_per)
    keep    <- intersect(c("x","y","z","n_locs_expected"), names(gt_df))
    f_df    <- gt_df[mol_idx, keep, drop = FALSE]
    f_df$mol_id <- orig_rows[mol_idx]
    n_f     <- nrow(f_df)

    # Step 3b: per-fluorophore detection filter.
    # Labeling efficiency applies to each individual fluorophore independently —
    # this is the biologically relevant quantity for MINFLUX: what fraction of
    # the theoretical TIDs are actually detected?  Molecules where every
    # fluorophore is missed simply produce no TIDs (equivalent to being
    # unlabeled), which emerges naturally without a separate molecule-level step.
    if (eff_pct < 100 && n_f > 0) {
      detected <- runif(n_f) < (eff_pct / 100)
      f_df  <- f_df[detected, , drop = FALSE]
      n_f   <- nrow(f_df)
      if (n_f == 0) return(f_df)
      mol_idx <- mol_idx[detected]
    }

    # Local helper: draw n_mol structured displacement vectors of given magnitude.
    # Respects ldir_type: "perp" (filament) → lateral only; "radial" (ring/sphere)
    # → outward from structure centre; NA → isotropic.
    # Returns a list with elements dx, dy, and (when use_3d) dz.
    struct_displace <- function(mag_vec) {
      ldt <- if (ldir_aware && "ldir_type" %in% names(gt_df))
               gt_df$ldir_type[1L] else NA_character_

      if (!is.na(ldt) && ldt == "perp" &&
          all(c("ldir_x","ldir_y") %in% names(gt_df))) {
        if (use_3d) {
          tang <- cbind(gt_df$ldir_x, gt_df$ldir_y, gt_df$ldir_z)
          u    <- matrix(rnorm(n_mol * 3L), ncol = 3L)
          u    <- u - rowSums(u * tang) * tang
          u_n  <- sqrt(rowSums(u^2)); u_n[u_n == 0] <- 1
          u    <- u / u_n
          list(dx = mag_vec * u[,1L], dy = mag_vec * u[,2L], dz = mag_vec * u[,3L])
        } else {
          sgn <- sample(c(-1L, 1L), n_mol, replace = TRUE)
          list(dx = mag_vec * sgn * (-gt_df$ldir_y),
               dy = mag_vec * sgn * ( gt_df$ldir_x))
        }
      } else if (!is.na(ldt) && ldt == "radial") {
        cx <- if ("struct_cx" %in% names(gt_df)) gt_df$struct_cx else rep(0, n_mol)
        cy <- if ("struct_cy" %in% names(gt_df)) gt_df$struct_cy else rep(0, n_mol)
        rx <- gt_df$x - cx; ry <- gt_df$y - cy
        if (use_3d) {
          cz <- if ("struct_cz" %in% names(gt_df)) gt_df$struct_cz else rep(0, n_mol)
          rz <- gt_df$z - cz
          rr <- sqrt(rx^2 + ry^2 + rz^2); rr[rr == 0] <- 1
          list(dx = mag_vec * rx/rr, dy = mag_vec * ry/rr, dz = mag_vec * rz/rr)
        } else {
          rr <- sqrt(rx^2 + ry^2); rr[rr == 0] <- 1
          list(dx = mag_vec * rx/rr, dy = mag_vec * ry/rr)
        }
      } else {
        if (use_3d) {
          dirs <- matrix(rnorm(n_mol * 3L), ncol = 3L)
          lens <- sqrt(rowSums(dirs^2)); lens[lens == 0] <- 1
          dirs <- dirs / lens
          list(dx = mag_vec * dirs[,1L], dy = mag_vec * dirs[,2L],
               dz = mag_vec * dirs[,3L])
        } else {
          theta <- runif(n_mol, 0, 2 * pi)
          list(dx = mag_vec * cos(theta), dy = mag_vec * sin(theta))
        }
      }
    }

    # Step 3c (Level 0): epitope offset from protein centroid.
    # The antibody binds a surface epitope, not the geometric centre of the protein.
    # Displacement magnitude is exactly the protein radius (r = 0.66 × MW^(1/3) nm);
    # direction follows the same structure-aware rules as the linker arm.
    if (protein_r_nm > 0 && n_mol > 0) {
      d0 <- struct_displace(rep(protein_r_nm, n_mol))
      f_df$x <- f_df$x + d0$dx[mol_idx]
      f_df$y <- f_df$y + d0$dy[mol_idx]
      if (use_3d) f_df$z <- f_df$z + d0$dz[mol_idx]
    }

    # Step 4 (Level 1): one complex-centre displacement per MOLECULE.
    # All fluorophores on the same molecule share this offset.
    # Direction is constrained by the local structure geometry when ldir_aware=TRUE.
    if (link_nm > 0 && n_mol > 0) {
      sigma_link <- link_nm * sqrt(pi / 2)   # half-normal → mean ≈ link_nm
      d1 <- struct_displace(abs(rnorm(n_mol, 0, sigma_link)))
      f_df$x <- f_df$x + d1$dx[mol_idx]
      f_df$y <- f_df$y + d1$dy[mol_idx]
      if (use_3d) f_df$z <- f_df$z + d1$dz[mol_idx]
    }

    # Step 5 (Level 2): per-fluorophore scatter within the complex body.
    # Skipped when dye_spread_nm = 0 (single-dye labels, custom dye_spread = 0).
    if (dye_spread_nm > 0 && n_f > 0) {
      f_df$x <- f_df$x + rnorm(n_f, 0, dye_spread_nm)
      f_df$y <- f_df$y + rnorm(n_f, 0, dye_spread_nm)
      if (use_3d) f_df$z <- f_df$z + rnorm(n_f, 0, dye_spread_nm)
    }

    rownames(f_df) <- NULL
    f_df
  }

  # ── Blinking / TID splitting ────────────────────────────────────────────────
  # Models dye re-activation: each fluorophore may appear as K additional TIDs,
  # K ~ Poisson(lambda). Extra rows carry the same mol_id and position; each TID
  # generates its own localizations independently via gen_localizations, giving
  # apparent over-counting proportional to (1 + lambda) × n_fluor_per_molecule.
  split_for_blinking <- function(fluoro_df, lambda) {
    if (lambda <= 0 || nrow(fluoro_df) == 0L) return(fluoro_df)
    extra <- rpois(nrow(fluoro_df), lambda)
    if (sum(extra) == 0L) return(fluoro_df)
    rbind(fluoro_df, fluoro_df[rep(seq_len(nrow(fluoro_df)), extra), , drop = FALSE])
  }

  # ── Localization generator ──────────────────────────────────────────────────
  # Uses n_locs_expected column in fluoro_df when present (image-guided mode);
  # otherwise draws from Poisson(trace_length) for all fluorophores.
  # sigma_loc_xy / sigma_loc_z are the per-localization precision in nm.
  #
  # use_geometric:  replace Poisson with Geometric(p = 1/lambda).  Same mean
  #   but heavy right tail: constant per-step bleaching gives many short traces
  #   and a few long ones, matching real MINFLUX trace-length distributions.
  #
  # scale_by_fluor: divide each fluorophore's lambda by the number of
  #   fluorophores on its molecule (requires mol_id column in fluoro_df).
  #   Keeps total locs per molecule ≈ trace_length regardless of labeling
  #   strategy — models a fixed observation budget per molecule.

  gen_localizations <- function(fluoro_df, trace_length,
                                 sigma_loc_xy, sigma_loc_z = sigma_loc_xy,
                                 use_geometric  = FALSE,
                                 scale_by_fluor = FALSE) {
    n      <- nrow(fluoro_df)
    use_3d <- "z" %in% names(fluoro_df)
    if (n == 0) {
      cols  <- if (use_3d) c("x","y","z","tid","tim")
               else        c("x","y","tid","tim")
      empty <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
      names(empty) <- cols
      return(empty)
    }

    lambda <- if ("n_locs_expected" %in% names(fluoro_df))
                fluoro_df$n_locs_expected
              else
                rep(max(1, trace_length), n)

    # Divide each fluorophore's budget by how many siblings it has on its molecule
    if (scale_by_fluor && "mol_id" %in% names(fluoro_df)) {
      mol_counts <- ave(rep(1, n), fluoro_df$mol_id, FUN = sum)
      lambda     <- lambda / pmax(1, mol_counts)
    }
    lambda <- pmax(1, lambda)

    if (use_geometric) {
      # rgeom(n, p) = failures before first success (min 0); +1 → min trace = 1
      # mean = 1/p = lambda, heavy right tail from constant bleach probability
      p     <- pmin(0.999, 1 / lambda)
      n_loc <- pmax(1L, rgeom(n, p) + 1L)
    } else {
      n_loc <- pmax(1L, rpois(n, lambda))
    }

    idx    <- rep(seq_len(n), n_loc)
    total  <- length(idx)
    rep_df <- fluoro_df[idx, intersect(c("x","y","z"), names(fluoro_df)),
                        drop = FALSE]
    rep_df$x <- rep_df$x + rnorm(total, 0, sigma_loc_xy)
    rep_df$y <- rep_df$y + rnorm(total, 0, sigma_loc_xy)
    if (use_3d) {
      rep_df$z <- rep_df$z + rnorm(total, 0, sigma_loc_z)
      rep_df$z <- pmax(-1000, pmin(1000, rep_df$z))
    }
    rep_df$tid <- as.character(idx)
    rep_df$tim <- seq(0, by = 10, length.out = total)
    rownames(rep_df) <- NULL
    rep_df
  }

  # ── Confocal PSF computation from optical parameters ────────────────────────
  # Gaussian approximation of the lateral confocal PSF sigma.
  #
  #   σ_ex   = 0.21 × λ_ex / NA               excitation PSF lateral sigma
  #   σ_em   = 0.21 × λ_em / NA               emission PSF lateral sigma
  #   1 AU   = 1.22 × λ_em / NA  (nm)         one Airy Unit diameter (sample plane)
  #   r_ph   = AU × (1 AU_nm) / 2             pinhole radius in nm
  #   σ_det  = √(σ_em² + (r_ph / 2.355)²)    detection sigma incl. pinhole size
  #   σ_conf = σ_ex × σ_det / √(σ_ex² + σ_det²)   product of two Gaussians
  #
  # For 640 nm / NA 1.4 / 0.83 AU:  σ ≈ 80 nm, FWHM ≈ 188 nm
  # For 560 nm / NA 1.4 / 1.0 AU:   σ ≈ 72 nm, FWHM ≈ 170 nm
  # For 488 nm / NA 1.4 / 1.0 AU:   σ ≈ 63 nm, FWHM ≈ 148 nm
  compute_confocal_sigma <- function(lambda_ex, lambda_em, na, au) {
    sigma_ex  <- 0.21 * lambda_ex / na
    sigma_em  <- 0.21 * lambda_em / na
    au_nm     <- 1.22 * lambda_em / na        # 1 AU diameter in nm (sample plane)
    ph_radius <- au * au_nm / 2               # pinhole radius in nm
    sigma_det <- sqrt(sigma_em^2 + (ph_radius / 2.355)^2)
    max(1, sigma_ex * sigma_det / sqrt(sigma_ex^2 + sigma_det^2))
  }


  # ── Diffraction-limited image renderer ──────────────────────────────────────
  # Renders each GT molecule as a 2D Gaussian PSF on a pixel grid, sums them,
  # adds uniform background, then draws Poisson shot noise.
  # Only computes the PSF within a 3-sigma radius per molecule (sparse update)
  # so performance stays acceptable for up to ~5 000 molecules.
  # Grid is capped to 500 × 500 px.
  # px_nm_user: when supplied, that pixel size is used directly (honouring a
  # user-chosen camera/scanner pixel size in nm).  If the structure extent
  # would push the grid over 500 px at that pixel size, px_nm is silently
  # increased just enough to stay within the cap.
  render_diffraction_image <- function(x_nm, y_nm, psf_sigma_nm,
                                        photons_per_mol, bg_photons_per_px,
                                        px_nm_user = NULL,
                                        photons_vec = NULL,
                                        extent = NULL) {
    n <- length(x_nm)
    if (n == 0) return(NULL)
    pad  <- 4 * psf_sigma_nm
    if (!is.null(extent) && length(extent) == 4 && all(is.finite(extent))) {
      # Explicit window (pinned to the FOV) so the confocal view spans the same
      # ROI as the MINFLUX view instead of cropping to the data bounding box.
      x0 <- extent[1]; x1 <- extent[2]; y0 <- extent[3]; y1 <- extent[4]
    } else {
      x0 <- min(x_nm) - pad; x1 <- max(x_nm) + pad
      y0 <- min(y_nm) - pad; y1 <- max(y_nm) + pad
    }
    span <- max(x1 - x0, y1 - y0, 1)
    # Use user pixel size if given; auto-size otherwise.
    # Either way, enforce the 500 px cap.
    px_nm <- if (!is.null(px_nm_user) && px_nm_user > 0)
               max(px_nm_user, span / 500)
             else
               max(psf_sigma_nm / 2, span / 500)
    nx    <- max(10L, min(500L, ceiling((x1 - x0) / px_nm)))
    ny    <- max(10L, min(500L, ceiling((y1 - y0) / px_nm)))
    # Keep pixels square at the requested size (don't rescale to fit)
    py_nm <- px_nm
    xs <- x0 + (seq_len(nx) - 0.5) * px_nm
    ys <- y0 + (seq_len(ny) - 0.5) * py_nm
    img  <- matrix(0, nrow = ny, ncol = nx)
    r_px <- ceiling(3 * psf_sigma_nm / min(px_nm, py_nm))
    for (i in seq_len(n)) {
      xi <- floor((x_nm[i] - x0) / px_nm) + 1L
      yi <- floor((y_nm[i] - y0) / py_nm) + 1L
      ix <- max(1L, xi - r_px):min(nx, xi + r_px)
      iy <- max(1L, yi - r_px):min(ny, yi + r_px)
      if (length(ix) == 0 || length(iy) == 0) next
      gx   <- exp(-0.5 * ((xs[ix] - x_nm[i]) / psf_sigma_nm)^2)
      gy   <- exp(-0.5 * ((ys[iy] - y_nm[i]) / psf_sigma_nm)^2)
      ph_i <- if (!is.null(photons_vec)) photons_vec[i] else photons_per_mol
      img[iy, ix] <- img[iy, ix] + ph_i * outer(gy, gx)
    }
    img <- img + bg_photons_per_px
    img <- matrix(rpois(length(img), as.vector(img)), nrow = ny)
    list(img = img, xs = xs, ys = ys, px_nm = px_nm)
  }

  # ── Linkage presets ──────────────────────────────────────────────────────────
  LINK_PRESETS <- c(
    "Direct tag (HaloTag / SNAP / click-EdU)" = 2,
    "Nanobody"                                 = 4,
    "Primary antibody (IgG)"                  = 10,
    "Primary + secondary antibody"            = 20,
    "Custom"                                   = NA_real_
  )

  # ── Phase 1 simulation (button-gated) ───────────────────────────────────────
  sim_result <- eventReactive(input$sim_run, {
    if (isTRUE(input$sim_fix_seed)) {
      seed <- suppressWarnings(as.integer(input$sim_seed))
      if (!is.finite(seed)) seed <- 42L
    } else {
      seed <- sample.int(.Machine$integer.max, 1L)
    }
    set.seed(seed)

    if (isTRUE(input$sim_fix_label_seed)) {
      label_seed <- suppressWarnings(as.integer(input$sim_label_seed))
      if (!is.finite(label_seed)) label_seed <- 1L
    } else {
      label_seed <- sample.int(.Machine$integer.max, 1L)
    }

    use_3d <- identical(input$sim_view_mode, "3d")
    n_prot <- max(1L, as.integer(input$sim_n_proteins))

    # Acquisition mode adds the FOV, multi-structure placement, timing model and
    # rate-based background; Normal mode is the simple structure + labeling view.
    acq_mode <- identical(input$sim_mode %||% "normal", "acquisition")
    # FOV active → background spreads across the whole FOV (see sim_use_fov()).
    use_fov  <- sim_use_fov()
    fov_x    <- max(1, as.numeric(input$sim_fov_x %||% 4000))
    fov_y    <- max(1, as.numeric(input$sim_fov_y %||% 4000))
    fov_z    <- max(1, as.numeric(input$sim_fov_z %||% 1000))

    img_b64  <- NULL
    img_w_nm <- NULL
    img_h_nm <- NULL

    if (identical(input$sim_shape, "image")) {
      shiny::validate(shiny::need(
        !is.null(input$sim_img_file), "Please upload a confocal image."))
      img_read <- read_image_gray(input$sim_img_file$datapath)
      img_gray <- img_read$mat
      img_b64  <- img_to_b64(img_read$img)
      px_nm    <- max(1, as.numeric(input$sim_img_px_size))
      img_w_nm <- ncol(img_gray) * px_nm
      img_h_nm <- nrow(img_gray) * px_nm
      gt <- gen_from_image(
        img_gray    = img_gray,
        n_proteins  = n_prot,
        px_size_nm  = px_nm,
        threshold   = suppressWarnings(as.numeric(input$sim_img_threshold)),
        min_trace   = max(1L, as.integer(input$sim_img_min_trace)),
        max_trace   = max(1L, as.integer(input$sim_img_max_trace)),
        use_3d      = use_3d,
        z_spread_nm = max(0, as.numeric(input$sim_img_z_spread %||% 200)),
        z_bg_nm     = max(0, as.numeric(input$sim_img_z_bg     %||% 20))
      )
    } else {
      gt <- switch(input$sim_shape,
        cluster = gen_cluster(
          n = n_prot,
          radius       = max(1, as.numeric(input$sim_cluster_radius)),
          distribution = input$sim_cluster_dist,
          use_3d       = use_3d
        ),
        ring = gen_ring(
          n = n_prot,
          diameter    = max(1, as.numeric(input$sim_ring_diam)),
          ring_width  = max(0, as.numeric(input$sim_ring_width)),
          z_thickness = if (use_3d)
                          max(0, as.numeric(input$sim_ring_z_thick))
                        else 0,
          use_3d = use_3d
        ),
        ellipse = gen_ellipse(
          n = n_prot,
          semi_a = max(1, as.numeric(input$sim_ellipse_a)),
          semi_b = max(1, as.numeric(input$sim_ellipse_b)),
          semi_c = if (use_3d) max(1, as.numeric(input$sim_ellipse_c)) else 1,
          filled = identical(input$sim_ellipse_fill, "filled"),
          use_3d = use_3d
        ),
        filament = gen_filament(
          n             = n_prot,
          spacing_nm    = max(0.1, as.numeric(input$sim_fil_spacing)),
          angle_deg     = suppressWarnings(as.numeric(input$sim_fil_angle)),
          jitter_lat_nm = max(0, as.numeric(input$sim_fil_jitter_lat)),
          jitter_ax_nm  = max(0, as.numeric(input$sim_fil_jitter_ax)),
          helix_r       = if (isTRUE(input$sim_fil_helix))
                            max(0, as.numeric(input$sim_fil_helix_r))
                          else 0,
          helix_pitch   = max(0.1, as.numeric(input$sim_fil_helix_pitch)),
          use_3d        = use_3d,
          z_angle_deg   = if (use_3d)
                            suppressWarnings(as.numeric(input$sim_fil_z_angle))
                          else 0
        ),
        hollow_sphere = gen_hollow_sphere(
          n           = n_prot,
          radius      = max(1, as.numeric(input$sim_hsph_radius)),
          shell_thick = max(0, as.numeric(input$sim_hsph_thick)),
          use_3d      = use_3d
        ),
        double_ring = gen_double_ring(
          n              = n_prot,
          diameter       = max(1, as.numeric(input$sim_dring_diam)),
          ring_width     = max(0, as.numeric(input$sim_dring_width)),
          separation_nm  = max(0, as.numeric(input$sim_dring_sep)),
          z_thickness    = if (use_3d)
                             max(0, as.numeric(input$sim_dring_z_thick))
                           else 0,
          use_3d         = use_3d
        ),
        arc = gen_arc(
          n              = n_prot,
          diameter       = max(1, as.numeric(input$sim_arc_diam)),
          ring_width     = max(0, as.numeric(input$sim_arc_width)),
          arc_start_deg  = suppressWarnings(as.numeric(input$sim_arc_start)),
          arc_extent_deg = max(1, as.numeric(input$sim_arc_extent)),
          z_thickness    = if (use_3d)
                             max(0, as.numeric(input$sim_arc_z_thick))
                           else 0,
          use_3d         = use_3d
        ),
        oligoring = gen_oligoring_stack(
          n_subunits       = max(2L, as.integer(input$sim_oligoring_nsub       %||% 6L)),
          diameter         = max(1,  as.numeric(input$sim_oligoring_diam       %||% 12)),
          n_stacks         = max(1L, as.integer(input$sim_oligoring_nstack     %||% 1L)),
          separation_nm    = max(0.5, as.numeric(input$sim_oligoring_sep       %||% 7)),
          angle_offset_deg = suppressWarnings(as.numeric(input$sim_oligoring_angle_off  %||% 0)),
          base_angle_deg   = suppressWarnings(as.numeric(input$sim_oligoring_base_angle %||% 0)),
          use_3d           = use_3d
        )
      )
      # Oligoring GT count is fixed by structure; update n_prot so downstream
      # labeling code (labeling efficiency, steric check) uses the correct N.
      if (identical(input$sim_shape, "oligoring"))
        n_prot <- nrow(gt)
    }

    # Clamp GT Z to ±1000 nm regardless of shape parameters
    if (use_3d && "z" %in% names(gt))
      gt$z <- pmax(-1000, pmin(1000, gt$z))

    # ── Random XY orientation ────────────────────────────────────────────────
    # Rotates the entire GT cloud by a uniformly sampled angle.  Skipped in
    # multi-copy mode (each copy is rotated independently there) and in image
    # mode (the image defines orientation implicitly).
    if (isTRUE(input$sim_rotate) && !identical(input$sim_shape, "image") &&
        !isTRUE(input$sim_multi_copy)) {
      rot_range <- max(0, as.numeric(input$sim_rotate_range %||% 180))
      angle_rad <- runif(1, -rot_range * pi / 180, rot_range * pi / 180)
      cos_a <- cos(angle_rad); sin_a <- sin(angle_rad)
      x_new <- gt$x * cos_a - gt$y * sin_a
      y_new <- gt$x * sin_a + gt$y * cos_a
      gt$x  <- x_new
      gt$y  <- y_new
      # Rotate stored filament tangent with the same angle
      if ("ldir_type" %in% names(gt) && identical(gt$ldir_type[1L], "perp")) {
        lx_new <- gt$ldir_x * cos_a - gt$ldir_y * sin_a
        ly_new <- gt$ldir_x * sin_a + gt$ldir_y * cos_a
        gt$ldir_x <- lx_new; gt$ldir_y <- ly_new
      }
    }

    # ── Multi-structure placement ────────────────────────────────────────────
    # Place N copies of the base structure at random positions in a rectangular
    # FOV (centred at the origin), each independently rotated.  Two drivers:
    #   • Normal mode  → optional "Multi-copy tiling" (square copies FOV).
    #   • Acquisition mode → "Number of structures" placed in the X×Y FOV.
    # Disabled for image mode because the image already fills the full FOV.
    if (acq_mode && !identical(input$sim_shape, "image")) {
      do_tiling <- TRUE
      n_copies  <- max(1L, as.integer(input$sim_acq_n_structures %||% 1L))
      tile_w    <- fov_x; tile_h <- fov_y
      rot_ea    <- TRUE
    } else if (!acq_mode && isTRUE(input$sim_multi_copy) &&
               !identical(input$sim_shape, "image")) {
      do_tiling <- TRUE
      n_copies  <- max(2L, as.integer(input$sim_copies_n   %||% 5L))
      tile_w    <- tile_h <- max(1, as.numeric(input$sim_copies_fov %||% 2000))
      rot_ea    <- isTRUE(input$sim_copies_rotate)
    } else {
      do_tiling <- FALSE
    }

    if (do_tiling && n_copies > 1L) {
      copies <- lapply(seq_len(n_copies), function(i) {
        cp <- gt
        if (rot_ea) {
          a     <- runif(1, 0, 2 * pi)
          ca    <- cos(a); sa <- sin(a)
          x_rot <- cp$x * ca - cp$y * sa
          y_rot <- cp$x * sa + cp$y * ca
          cp$x  <- x_rot; cp$y <- y_rot
          # Rotate filament tangent with each per-copy rotation
          if ("ldir_type" %in% names(cp) && identical(cp$ldir_type[1L], "perp")) {
            lx_r <- cp$ldir_x * ca - cp$ldir_y * sa
            ly_r <- cp$ldir_x * sa + cp$ldir_y * ca
            cp$ldir_x <- lx_r; cp$ldir_y <- ly_r
          }
        }
        # Keep the first structure near the centre (~0-100 nm), like a real
        # MINFLUX ROI is usually centred on a feature; the rest tile randomly.
        if (i == 1L) {
          off_x <- runif(1, -50, 50)
          off_y <- runif(1, -50, 50)
        } else {
          off_x <- runif(1, -tile_w / 2, tile_w / 2)
          off_y <- runif(1, -tile_h / 2, tile_h / 2)
        }
        cp$x <- cp$x + off_x
        cp$y <- cp$y + off_y
        # Track per-copy structure centre for radial displacement
        if ("ldir_type" %in% names(cp) && identical(cp$ldir_type[1L], "radial")) {
          if ("struct_cx" %in% names(cp)) cp$struct_cx <- cp$struct_cx + off_x
          if ("struct_cy" %in% names(cp)) cp$struct_cy <- cp$struct_cy + off_y
        }
        if (use_3d && "z" %in% names(cp))
          cp$z <- pmax(-1000, pmin(1000, cp$z))
        cp
      })
      gt <- do.call(rbind, copies)
      rownames(gt) <- NULL
    }

    # ── Labeling: preset or custom ───────────────────────────────────────────
    set.seed(label_seed)
    label_key <- input$sim_label_type %||% "Custom"
    if (identical(label_key, "Custom")) {
      n_fluor       <- max(1L, as.integer(input$sim_fluor_custom     %||% 1L))
      link_nm       <- max(0,  as.numeric(input$sim_link_custom_nm   %||% 5))
      dye_spread_nm <- max(0,  as.numeric(input$sim_dye_spread_custom %||% 0))
      footprint_nm  <- max(0,  as.numeric(input$sim_footprint_custom  %||% (link_nm * 1.5)))
      blink_lambda  <- max(0,  as.numeric(input$sim_blink_mean        %||% 0.3))
      lbl_note      <- sprintf("%d fluorophore(s)/molecule, %.0f nm linker, %.0f nm dye spread",
                               n_fluor, link_nm, dye_spread_nm)
    } else {
      p             <- LABEL_PRESETS[[label_key]]
      n_fluor       <- p$n_fluor
      link_nm       <- p$link_nm
      dye_spread_nm <- p$dye_spread_nm
      footprint_nm  <- p$footprint_nm
      blink_lambda  <- if (!is.null(p$blink_lambda)) p$blink_lambda else 0.3
      lbl_note      <- p$note
    }
    eff_pct        <- max(1, min(100, as.numeric(input$sim_label_eff)))
    protein_mw_kda <- max(0, as.numeric(input$sim_protein_mw %||% 0))
    protein_r_nm   <- if (protein_mw_kda > 0) 0.66 * protein_mw_kda^(1/3) else 0
    fluoro  <- expand_to_fluorophores(gt, eff_pct, n_fluor, link_nm, dye_spread_nm,
                                       ldir_aware   = isTRUE(input$sim_ldir_aware),
                                       protein_r_nm = protein_r_nm)
    # Unique physical fluorophores (one dye = one position), before any blinking
    # or rate-based over-counting.  The confocal view renders from THIS so a
    # re-detected molecule is not rendered brighter than it physically is.
    fluoro_pool <- fluoro

    # ── Blinking / TID splitting (Normal mode only) ──────────────────────────
    # In Acquisition mode over-counting comes from the detection rate × duration
    # below, so the λ-based blink split is skipped there.
    n_fluoro_base <- nrow(fluoro)   # unique fluorophores (detectable pool)
    if (!acq_mode && isTRUE(input$sim_blink_model)) {
      bl_lambda <- max(0, as.numeric(input$sim_blink_mean %||% blink_lambda))
      fluoro    <- split_for_blinking(fluoro, bl_lambda)
    }

    # ── Localization precision: manual or CRLB-derived ──────────────────────
    # CRLB mode: σ ≈ L / (2√N_eff),  N_eff = N_photons × (1 − bg_fraction).
    # Simplified 2D MINFLUX CRLB (equal pattern arms, ideal PSF, on-axis).
    if (identical(input$sim_prec_mode %||% "manual", "crlb")) {
      crlb_L   <- max(1,    as.numeric(input$sim_crlb_L       %||% 100))
      crlb_ph  <- max(1,    as.numeric(input$sim_crlb_photons %||% 50))
      crlb_bg  <- max(0, min(0.99, as.numeric(input$sim_crlb_bg %||% 0.05)))
      n_eff    <- max(1, crlb_ph * (1 - crlb_bg))
      sigma_loc_xy <- max(0.1, crlb_L / (2 * sqrt(n_eff)))
      sigma_loc_z  <- if (use_3d) {
        z_fac <- max(0.1, as.numeric(input$sim_crlb_z_factor %||% 0.7))
        max(0.1, sigma_loc_xy * z_fac)
      } else sigma_loc_xy
    } else {
      sigma_loc_xy <- max(0.1, as.numeric(input$sim_precision_xy))
      sigma_loc_z  <- if (use_3d)
                        max(0.1, as.numeric(input$sim_precision_z %||% sigma_loc_xy))
                      else sigma_loc_xy
    }
    trace_length <- max(1, as.numeric(input$sim_trace_length))

    # ── Acquisition window + over-counting (Acquisition mode only) ────────────
    # The user sets the duration (hours + minutes) AND the structure detection
    # rate (traces/hour).  Requested structure traces = rate × duration.  Each
    # fluorophore yields a finite number of detections before PERMANENT bleaching
    # (mean = sim_acq_blinks), so the total is capped at that photon budget: a
    # short run under-samples the pool (low % resolved), a longer run resolves
    # more and over-counts, and beyond the budget the structure plateaus (dead)
    # while background keeps accumulating.  Detections are drawn WITHOUT
    # replacement from the per-fluorophore capacity, so molecules with more
    # fluorophores / more blinks are found sooner.
    # (The resolution/saturation curves are computed live in sim_acq_live();
    # here we only need the over-counted detection list for the MINFLUX view.)
    fov_area_nm2 <- fov_x * fov_y
    acq_T_sec    <- 0
    struct_rate  <- max(0, as.numeric(input$sim_struct_rate %||% 20))   # traces/hour
    n_pool       <- nrow(fluoro)        # unique fluorophores (detectable pool)
    if (acq_mode) {
      hrs       <- max(0, as.numeric(input$sim_acq_hours   %||% 0))
      mins      <- max(0, as.numeric(input$sim_acq_minutes %||% 0))
      acq_T_sec <- max(1, (hrs * 60 + mins) * 60)
      blinks    <- max(1, as.numeric(input$sim_acq_blinks %||% 5))   # detections/fluorophore
      dec_min   <- max(0.1, as.numeric(input$sim_acq_decay_min  %||% 15))
      flr_frac  <- max(0, min(4, as.numeric(input$sim_acq_floor_frac %||% 0.25)))
      n_req     <- as.integer(round(                                # absolute-rate count
                     acq_expected_count(struct_rate / 60, acq_T_sec,  # per hour → per min
                                        dec_min, flr_frac)))
      if (n_pool > 0 && n_req > 0) {
        k_i        <- 1L + rpois(n_pool, max(0, blinks - 1))         # per-fluorophore budget
        struct_cap <- sum(k_i)
        n_s        <- min(n_req, struct_cap, 200000L)
        det_idx    <- if (struct_cap <= 5e6)
                        sample(rep(seq_len(n_pool), k_i), n_s)        # exact capacity
                      else
                        sample(seq_len(n_pool), n_s, replace = TRUE, prob = k_i)
        fluoro     <- fluoro[det_idx, , drop = FALSE]
        rownames(fluoro) <- NULL
      } else {
        fluoro     <- fluoro[0, , drop = FALSE]
      }
    }

    locs <- gen_localizations(fluoro, trace_length, sigma_loc_xy, sigma_loc_z,
                               use_geometric  = isTRUE(input$sim_tl_geometric),
                               scale_by_fluor = isTRUE(input$sim_tl_scale_fluor))

    # ── Background / noise TIDs ──────────────────────────────────────────────
    # Acquisition mode: count = background detection rate (traces/hour) × duration,
    # mirroring the structure — longer runs accumulate more non-specific TIDs.
    # Normal mode: a plain fixed count.  Background is spread across the FOV when
    # set, else the multi-copy FOV, else the GT bounding box + 10 % padding.
    # Each background detection is a fresh random position. TID numbering
    # continues from the last real TID.
    n_bg   <- 0L
    bg_pos <- NULL
    bg_rate <- max(0, as.numeric(input$sim_bg_rate %||% 5))   # traces/hour
    if (isTRUE(input$sim_bg_noise)) {
      copies_mode <- !acq_mode && isTRUE(input$sim_multi_copy) &&
                     !identical(input$sim_shape, "image")
      if (use_fov) {
        bg_xc <- 0; bg_yc <- 0
        bg_hx <- fov_x / 2; bg_hy <- fov_y / 2
      } else if (copies_mode) {
        fov_half <- max(1, as.numeric(input$sim_copies_fov %||% 2000)) / 2
        bg_xc <- 0; bg_yc <- 0
        bg_hx <- fov_half; bg_hy <- fov_half
      } else {
        pad   <- 0.1
        x_rng <- range(gt$x, na.rm = TRUE); y_rng <- range(gt$y, na.rm = TRUE)
        x_pad <- diff(x_rng) * pad; y_pad <- diff(y_rng) * pad
        bg_xc <- mean(x_rng); bg_yc <- mean(y_rng)
        bg_hx <- diff(x_rng) / 2 + x_pad; bg_hy <- diff(y_rng) / 2 + y_pad
      }

      if (acq_mode) {
        n_bg <- as.integer(round(bg_rate * acq_T_sec / 3600))   # per hour × hours
      } else {
        n_bg <- as.integer(input$sim_bg_n %||% 20L)
      }
      n_bg <- max(0L, min(200000L, n_bg))

      if (n_bg > 0L) {
        bg_df <- data.frame(x = runif(n_bg, bg_xc - bg_hx, bg_xc + bg_hx),
                            y = runif(n_bg, bg_yc - bg_hy, bg_yc + bg_hy))
        if (use_3d && "z" %in% names(gt)) {
          z_rng  <- range(gt$z, na.rm = TRUE)
          # FOV mode → fill the FOV Z slab; otherwise track the GT extent, with
          # a precision-based fallback for degenerate (planar) structures.
          z_half <- if (use_fov)
                      fov_z / 2
                    else if (diff(z_rng) < 1) max(50, sigma_loc_z * 10)
                    else diff(z_rng) / 2 * 1.1
          z_mid  <- if (use_fov) 0 else mean(z_rng)
          bg_df$z <- runif(n_bg, z_mid - z_half, z_mid + z_half)
        }
        bg_df$mol_id <- NA_integer_
        # Save XY(Z) positions so the confocal view can render bg spots
        bg_pos <- bg_df[, intersect(c("x","y","z"), names(bg_df)), drop = FALSE]
        # Each background TID gets its own locs via gen_localizations
        bg_locs <- gen_localizations(bg_df, trace_length, sigma_loc_xy,
                                      sigma_loc_z,
                                      use_geometric  = isTRUE(input$sim_tl_geometric),
                                      scale_by_fluor = FALSE)
        # Offset TID numbers so they don't clash with real TIDs
        bg_locs$tid <- as.character(as.integer(bg_locs$tid) + nrow(fluoro))
        locs <- rbind(locs, bg_locs)
      }
    }

    # ── RMSE: per-TID displacement error (GT position vs trace median) ──────────
    # Only computed for real TIDs (1 … nrow(fluoro)); background TIDs excluded.
    # Pre-split locs by TID for O(1) lookup per trace — essential when multi-
    # copy produces thousands of TIDs.  Capped at 1 000 sampled TIDs so the
    # computation stays fast even with very large multi-copy scenes.
    rmse_vals <- NULL
    if ("mol_id" %in% names(fluoro) && nrow(gt) > 0 && nrow(locs) > 0) {
      n_real    <- nrow(fluoro)
      tid_samp  <- if (n_real > 1000L) sample(n_real, 1000L) else seq_len(n_real)
      real_rows <- which(suppressWarnings(as.integer(locs$tid)) <= n_real &
                         !is.na(suppressWarnings(as.integer(locs$tid))))
      locs_by_tid <- split(real_rows, locs$tid[real_rows])
      errs <- vapply(tid_samp, function(k) {
        rows <- locs_by_tid[[as.character(k)]]
        if (is.null(rows) || length(rows) == 0) return(NA_real_)
        lk <- locs[rows, , drop = FALSE]
        gi <- fluoro$mol_id[k]
        if (is.na(gi) || gi < 1L || gi > nrow(gt)) return(NA_real_)
        dx <- median(lk$x) - gt$x[gi]
        dy <- median(lk$y) - gt$y[gi]
        dz <- if (use_3d && "z" %in% names(lk) && "z" %in% names(gt))
                median(lk$z) - gt$z[gi]
              else 0
        sqrt(dx^2 + dy^2 + dz^2)
      }, numeric(1L))
      rmse_vals <- errs[!is.na(errs)]
    }

    # ── Steric exclusion check ───────────────────────────────────────────────
    # Compute mean nearest-neighbour distance among GT molecule positions and
    # compare it to the label complex footprint.  A footprint larger than the
    # protein spacing means antibodies physically cannot label every molecule.
    steric_warn <- NULL
    gt_nnd      <- NA_real_
    if (nrow(gt) >= 2L) {
      gt_pts <- as.matrix(gt[, intersect(c("x","y","z"), names(gt)), drop = FALSE])
      samp_i <- if (nrow(gt_pts) > 300L) sample(nrow(gt_pts), 300L) else seq_len(nrow(gt_pts))
      d_mat  <- as.matrix(dist(gt_pts[samp_i, , drop = FALSE]))
      diag(d_mat) <- Inf
      gt_nnd <- mean(apply(d_mat, 1, min), na.rm = TRUE)
    }
    # Effective exclusion zone = label footprint + protein diameter.
    # Proteins physically can't overlap, so even before the antibody is considered
    # the centre-to-centre spacing must be ≥ 2 × protein radius.
    eff_exclusion_nm <- max(footprint_nm, 2 * protein_r_nm, na.rm = TRUE)
    if (!is.na(eff_exclusion_nm) && !is.na(gt_nnd)) {
      if (gt_nnd < eff_exclusion_nm * 0.5) {
        steric_warn <- list(
          level = "danger",
          msg   = sprintf(
            paste0("Severe steric exclusion: combined exclusion zone (~%.0f nm: ",
                   "%.0f nm label footprint + %.0f nm protein diameter) exceeds ",
                   "mean protein spacing (~%.0f nm). Full stoichiometric labeling is ",
                   "physically impossible at this density."),
            eff_exclusion_nm, footprint_nm, 2 * protein_r_nm, gt_nnd)
        )
      } else if (gt_nnd < eff_exclusion_nm) {
        steric_warn <- list(
          level = "warning",
          msg   = sprintf(
            paste0("Steric exclusion likely: combined exclusion zone (~%.0f nm) ",
                   "approaches mean protein spacing (~%.0f nm). Some molecules may be ",
                   "inaccessible; effective labeling efficiency could be lower than expected."),
            eff_exclusion_nm, gt_nnd)
        )
      }
    }

    # ── Snapshot detection times (Acquisition mode) ──────────────────────────
    # Assign front-loaded times to THIS realization's TIDs so the MINFLUX Time
    # colour mode and CSV export have realistic timestamps.  The Acquisition-tab
    # saturation/rate curves are computed live in sim_acq_live() (decoupled so
    # they react to rate/duration/background without re-Generating).
    if (acq_mode && nrow(locs) > 0 && acq_T_sec > 0) {
      all_tids   <- unique(locs$tid)
      n_real     <- nrow(fluoro)
      decay_min  <- max(0.1, as.numeric(input$sim_acq_decay_min  %||% 15))
      floor_frac <- max(0, min(4, as.numeric(input$sim_acq_floor_frac %||% 0.25)))
      tid_times  <- setNames(rep(NA_real_, length(all_tids)), all_tids)
      real_ids   <- intersect(as.character(seq_len(n_real)), all_tids)
      if (length(real_ids) > 0)
        tid_times[real_ids] <- sim_sample_detection_times(
          length(real_ids), acq_T_sec, decay_min * 60, floor_frac)
      bg_ids <- setdiff(all_tids, real_ids)
      if (length(bg_ids) > 0)
        tid_times[bg_ids] <- runif(length(bg_ids), 0, acq_T_sec)
      within   <- ave(seq_len(nrow(locs)), locs$tid,
                      FUN = function(i) (seq_along(i) - 1L) * 0.01)
      locs$tim <- unname(tid_times[locs$tid]) + within
    }

    list(gt = gt, fluoro = fluoro, fluoro_pool = fluoro_pool,
         locs = locs,
         link_nm = link_nm, n_fluor = n_fluor, lbl_note = lbl_note,
         protein_r_nm = protein_r_nm,
         sigma_loc_xy = sigma_loc_xy, sigma_loc_z = sigma_loc_z,
         n_prot = n_prot, n_fluoro = nrow(fluoro),
         n_fluoro_base = n_fluoro_base,
         n_locs = nrow(locs), n_bg = n_bg, use_3d = use_3d,
         img_b64 = img_b64,
         img_w_nm = img_w_nm, img_h_nm = img_h_nm,
         rmse_vals = rmse_vals, bg_pos = bg_pos,
         steric_warn = steric_warn,
         footprint_nm = footprint_nm,
         gt_nnd = gt_nnd, acq_mode = acq_mode,
         blink_active = (!acq_mode) && isTRUE(input$sim_blink_model),
         blink_lambda_used = if ((!acq_mode) && isTRUE(input$sim_blink_model))
                               max(0, as.numeric(input$sim_blink_mean %||% blink_lambda))
                             else 0)
  }, ignoreNULL = TRUE)

  # ── Live acquisition timeline ────────────────────────────────────────────────
  # Recomputes the saturation / detection-rate curves from the last-generated
  # fluorophore pool whenever the acquisition parameters change — so the
  # Acquisition tab reacts to rate / duration / background / bleaching live, with
  # no need to re-Generate the structure.  A fixed seed keeps the curve stable
  # for a given parameter set (it only moves when an input actually changes).
  sim_acq_live <- reactive({
    if (!identical(input$sim_mode %||% "normal", "acquisition")) return(NULL)
    res <- sim_result()
    if (is.null(res) || is.null(res$fluoro_pool)) return(NULL)
    pool_mol <- res$fluoro_pool$mol_id
    if (length(pool_mol) == 0) return(NULL)
    hrs   <- max(0, as.numeric(input$sim_acq_hours   %||% 0))
    mins  <- max(0, as.numeric(input$sim_acq_minutes %||% 0))

    # Emitter-crowding weight per molecule: a molecule with many other labeled
    # molecules inside the donut radius is harder to single out, so its effective
    # discovery rate is scaled by 1/(1 + neighbours).  Aligned to the sorted
    # unique mol_ids (= the order compute_acq_timeline builds f_m in).
    crowd_m <- NULL
    if (isTRUE(input$sim_crowding)) {
      u  <- sort(unique(pool_mol))
      gp <- res$gt[u, , drop = FALSE]
      cols <- intersect(c("x", "y", "z"), names(gp))
      if (length(u) >= 2 && length(cols) >= 2) {
        R   <- max(1, as.numeric(input$sim_crowd_radius %||% 50))
        pts <- as.matrix(gp[, cols, drop = FALSE])
        nn  <- RANN::nn2(pts, pts, k = min(nrow(pts), 64L),
                         searchtype = "radius", radius = R)
        neigh   <- rowSums(nn$nn.idx > 0) - 1L
        crowd_m <- 1 / (1 + neigh)
      }
    }

    set.seed(1L)                          # stable curve per parameter set
    compute_acq_timeline(
      pool_mol    = pool_mol,
      n_mol_total = nrow(res$gt),
      struct_rate = max(0, as.numeric(input$sim_struct_rate %||% 20)),
      T_sec       = max(1, (hrs * 60 + mins) * 60),
      blinks      = max(1, as.numeric(input$sim_acq_blinks %||% 5)),
      decay_min   = max(0.1, as.numeric(input$sim_acq_decay_min %||% 15)),
      floor_frac  = max(0, min(4, as.numeric(input$sim_acq_floor_frac %||% 0.25))),
      bg_rate     = max(0, as.numeric(input$sim_bg_rate %||% 5)),
      resolve_q   = as.numeric(input$sim_acq_resolve_q %||% 0.9),
      crowd_m     = crowd_m)
  })

  # ── Calibrate from a real MINFLUX run ────────────────────────────────────────
  # Load a .npy / .json / .csv run, measure its detection-rate profile, duration
  # and FOV, pre-fill the Acquisition inputs, and keep the measured rate curve to
  # overlay on the Acquisition plot.
  sim_real_rv <- reactiveVal(NULL)

  # Split the measured total rate into structure + constant background given the
  # fraction of traces that are on-structure, and push struct_rate / bg_rate /
  # late-rate ratio.  Background is the constant offset; the structure is the
  # time-varying part on top — its rate can FALL (bleaching, ratio < 1) or RISE
  # (buffer activating, ratio > 1) toward its late value.  Background is capped
  # at the minimum of the run so the early/late structure rate stays positive.
  apply_real_rates <- function(a, frac) {
    frac  <- max(0, min(1, frac))
    bg    <- min((1 - frac) * a$rate_mean,
                 0.9 * min(a$rate_peak, a$rate_floor))      # ≤ the run minimum
    bg    <- max(0, bg)
    s0    <- max(0.1, a$rate_peak - bg)                     # initial structure rate
    ratio <- max(0, min(4, (a$rate_floor - bg) / s0))       # late / early
    updateNumericInput(session, "sim_struct_rate", value = round(s0))
    updateNumericInput(session, "sim_bg_rate",     value = round(bg))
    updateSliderInput (session, "sim_acq_floor_frac", value = round(ratio, 2))
  }

  observeEvent(input$sim_real_file, {
    f <- input$sim_real_file
    req(f)
    ext <- get_ext(f$name)
    df <- tryCatch(
      withProgress(message = "Reading run…", value = 0.2, {
        switch(ext,
          npy  = read_minflux_npy_as_df(f$datapath),
          json = read_minflux_json_as_df(f$datapath),
          csv  = read.csv(f$datapath, check.names = FALSE),
          stop("Unsupported file type (use .npy, .json or .csv)."))
      }),
      error = function(e) {
        showNotification(paste("Load failed:", conditionMessage(e)),
                         type = "error", duration = 8); NULL
      })
    if (is.null(df)) return()
    a <- analyze_real_run(df)
    if (is.null(a)) {
      showNotification("Could not analyze run (need tid, tim, x, y).",
                       type = "error", duration = 6); return()
    }
    sim_real_rv(a)

    # Pre-fill the Acquisition inputs from the measured values.
    hrs  <- floor(a$duration_hr)
    mins <- min(59L, as.integer(round((a$duration_hr - hrs) * 60)))
    updateRadioButtons(session, "sim_mode", selected = "acquisition")
    updateNumericInput(session, "sim_acq_hours",     value = hrs)
    updateNumericInput(session, "sim_acq_minutes",   value = mins)
    updateNumericInput(session, "sim_acq_decay_min", value = round(a$tau_min, 1))
    updateNumericInput(session, "sim_fov_x", value = ceiling(a$fov_x * 1.1 / 100) * 100)
    updateNumericInput(session, "sim_fov_y", value = ceiling(a$fov_y * 1.1 / 100) * 100)
    if (a$fov_z > 1)
      updateNumericInput(session, "sim_fov_z",
                         value = ceiling(a$fov_z * 1.1 / 100) * 100)
    # Structure/background split from the spatial classification (overridable).
    updateSliderInput(session, "sim_real_struct_frac",
                      value = round(a$struct_frac, 2))
    apply_real_rates(a, a$struct_frac)

    # Measured localization precision → Manual σ; photons + effective L → CRLB.
    if (is.finite(a$sigma_x) && is.finite(a$sigma_y))
      updateNumericInput(session, "sim_precision_xy",
                         value = round(mean(c(a$sigma_x, a$sigma_y)), 1))
    if (is.finite(a$sigma_z))
      updateNumericInput(session, "sim_precision_z", value = round(a$sigma_z, 1))
    if (is.finite(a$photons))
      updateNumericInput(session, "sim_crlb_photons", value = round(a$photons))
    if (is.finite(a$L_eff))
      updateNumericInput(session, "sim_crlb_L", value = round(a$L_eff))

    showNotification("Acquisition + precision set from the real run.",
                     type = "message", duration = 5)
  })

  # Manual override of the structure/background split (re-derives the rates).
  observeEvent(input$sim_real_struct_frac, {
    a <- sim_real_rv()
    if (!is.null(a)) apply_real_rates(a, input$sim_real_struct_frac)
  }, ignoreInit = TRUE)

  output$sim_real_status <- renderUI({
    a <- sim_real_rv()
    if (is.null(a))
      return(tags$small(style = "color:#888;",
        "Load a .npy / .json MINFLUX run to calibrate the model from real data."))
    tags$div(class = "alert alert-success",
             style = "padding:6px 10px; font-size:0.82em; margin:4px 0;",
      tags$b(sprintf("%d traces over %.1f h", a$n_traces, a$duration_hr)),
      tags$br(),
      sprintf("Total rate: peak %.0f, mean %.0f, floor %.0f /hr  (τ≈%.0f min)",
              a$rate_peak, a$rate_mean, a$rate_floor, a$tau_min),
      tags$br(),
      sprintf("Spatial split: %.0f%% on-structure / %.0f%% background",
              100 * a$struct_frac, 100 * (1 - a$struct_frac)),
      tags$br(),
      sprintf("FOV: %.0f × %.0f%s nm", a$fov_x, a$fov_y,
              if (a$fov_z > 1) sprintf(" × %.0f", a$fov_z) else ""),
      if (is.finite(a$sigma_x)) tagList(tags$br(),
        sprintf("Precision σ: x %.1f, y %.1f%s nm%s",
                a$sigma_x, a$sigma_y,
                if (is.finite(a$sigma_z)) sprintf(", z %.1f", a$sigma_z) else "",
                if (is.finite(a$photons))
                  sprintf("  |  %.0f photons/loc → L≈%.0f nm",
                          a$photons, a$L_eff)
                else "")))
  })

  # ── Compare-strategies reactive ─────────────────────────────────────────────
  # Applies two labeling strategies to the same GT from sim_result().
  # Geometry is NOT re-run; seeds come from current label_seed + offsets.
  compare_result <- eventReactive(input$sim_compare_run, {
    res <- sim_result()
    shiny::validate(shiny::need(
      !is.null(res),
      "Run the main simulation first to generate a ground truth."))

    gt   <- res$gt
    use_3d <- res$use_3d
    sigma_loc_xy <- res$sigma_loc_xy
    sigma_loc_z  <- res$sigma_loc_z
    trace_length <- max(1, as.numeric(input$sim_trace_length))

    # Base seed: use the fixed label seed if set, else a fresh one per compare run
    base_seed <- if (isTRUE(input$sim_fix_label_seed)) {
      suppressWarnings(as.integer(input$sim_label_seed))
    } else {
      sample.int(.Machine$integer.max, 1L)
    }
    if (!is.finite(base_seed)) base_seed <- 1L

    cmp_protein_mw_kda <- max(0, as.numeric(input$sim_protein_mw %||% 0))
    cmp_protein_r_nm   <- if (cmp_protein_mw_kda > 0) 0.66 * cmp_protein_mw_kda^(1/3) else 0

    run_strategy <- function(label_key, eff_pct, seed_val, suffix) {
      set.seed(seed_val)
      if (identical(label_key, "Custom")) {
        n_fluor       <- max(1L, as.integer(input[[paste0("sim_cmp_", suffix, "_fluor")]] %||% 2L))
        link_nm       <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_link")]]  %||% 10))
        dye_spread_nm <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_spread")]] %||% 0))
        footprint_nm  <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_foot")]]  %||% 8))
        blink_lambda  <- max(0,  as.numeric(input$sim_blink_mean %||% 0.3))
        lbl_note      <- sprintf("Custom: %d fluor / %.0f nm linker",
                                 n_fluor, link_nm)
      } else {
        p             <- LABEL_PRESETS[[label_key]]
        n_fluor       <- p$n_fluor
        link_nm       <- p$link_nm
        dye_spread_nm <- p$dye_spread_nm
        footprint_nm  <- p$footprint_nm
        blink_lambda  <- p$blink_lambda %||% 0.3
        lbl_note      <- p$note
      }

      fluoro <- expand_to_fluorophores(
        gt, eff_pct, n_fluor, link_nm, dye_spread_nm,
        ldir_aware   = isTRUE(input$sim_ldir_aware),
        protein_r_nm = cmp_protein_r_nm)

      if (isTRUE(input$sim_blink_model))
        fluoro <- split_for_blinking(
          fluoro, max(0, as.numeric(input$sim_blink_mean %||% blink_lambda)))

      locs <- gen_localizations(fluoro, trace_length, sigma_loc_xy, sigma_loc_z,
                                 use_geometric  = isTRUE(input$sim_tl_geometric),
                                 scale_by_fluor = isTRUE(input$sim_tl_scale_fluor))

      # RMSE (capped at 500 TIDs for speed)
      rmse_vals <- NULL
      if ("mol_id" %in% names(fluoro) && nrow(gt) > 0 && nrow(locs) > 0) {
        n_samp <- min(nrow(fluoro), 500L)
        real_rows <- which(suppressWarnings(as.integer(locs$tid)) <= nrow(fluoro) &
                           !is.na(suppressWarnings(as.integer(locs$tid))))
        locs_by_tid <- split(real_rows, locs$tid[real_rows])
        errs <- vapply(seq_len(n_samp), function(k) {
          rows <- locs_by_tid[[as.character(k)]]
          if (is.null(rows) || length(rows) == 0L) return(NA_real_)
          lk <- locs[rows, , drop = FALSE]
          gi <- fluoro$mol_id[k]
          if (is.na(gi) || gi < 1L || gi > nrow(gt)) return(NA_real_)
          dx <- median(lk$x) - gt$x[gi]
          dy <- median(lk$y) - gt$y[gi]
          sqrt(dx^2 + dy^2)
        }, numeric(1L))
        rmse_vals <- errs[!is.na(errs)]
      }

      list(fluoro = fluoro, locs = locs, link_nm = link_nm,
           n_fluor = n_fluor, n_fluoro = nrow(fluoro),
           lbl_note = lbl_note, footprint_nm = footprint_nm,
           protein_r_nm = cmp_protein_r_nm,
           rmse_vals = rmse_vals)
    }

    eff_a <- max(1, min(100, as.numeric(input$sim_cmp_eff_a %||% 80)))
    eff_b <- max(1, min(100, as.numeric(input$sim_cmp_eff_b %||% 80)))
    a <- run_strategy(input$sim_cmp_label_a %||%
                        "Direct tag (SNAP / Halo / Click, DOL=1)",
                      eff_a, base_seed, "a")
    b <- run_strategy(input$sim_cmp_label_b %||%
                        "Primary monoclonal IgG + secondary IgG",
                      eff_b, base_seed + 1L, "b")

    list(gt = gt, use_3d = use_3d, a = a, b = b,
         sigma_loc_xy = sigma_loc_xy,
         lbl_a = input$sim_cmp_label_a %||%
                   "Direct tag (SNAP / Halo / Click, DOL=1)",
         lbl_b = input$sim_cmp_label_b %||%
                   "Primary monoclonal IgG + secondary IgG")
  }, ignoreNULL = TRUE)

  # Per-localization local density: number of neighbours within radius_nm
  # (self excluded). Used by the "Local density" colour mode. k is capped so
  # very dense scenes stay fast; counts saturate at k-1 in that regime.
  sim_local_density <- function(locs, radius_nm) {
    cols <- intersect(c("x", "y", "z"), names(locs))
    pts  <- as.matrix(locs[, cols, drop = FALSE])
    n    <- nrow(pts)
    if (n < 2) return(rep(0, n))
    k    <- min(n, 256L)
    nn   <- RANN::nn2(pts, pts, k = k)
    rowSums(nn$nn.dists <= radius_nm) - 1L
  }

  # ── Phase 1 plot ─────────────────────────────────────────────────────────────
  output$sim_plot <- renderPlotly({
    res <- sim_result()
    shiny::validate(
      shiny::need(!is.null(res), "Click 'Generate' to run the simulation."))

    dark    <- isTRUE(input$sim_dark_mode)
    bg_col  <- if (dark) "#111111" else "#FFFFFF"
    fg_col  <- if (dark) "#DDDDDD" else "#444444"
    use_3d  <- isTRUE(res$use_3d)
    tr_type <- if (use_3d) "scatter3d" else "scattergl"
    pt_size <- max(1L, as.integer(input$sim_point_size))
    gc      <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"
    show_grid <- !isTRUE(input$sim_hide_grid)

    # Field of view — display window only, nothing is cropped. Axes are pinned to
    # a real-camera-style FOV; the 3D aspect ratio reflects the true X:Y:Z
    # proportions so Z appears compressed (4000×4000×1000). Always on in
    # Acquisition mode; in Normal mode it's opt-in via "Pin axes to a FOV (ROI)".
    use_fov <- sim_use_fov()
    fov_x   <- max(1, as.numeric(input$sim_fov_x %||% 4000))
    fov_y   <- max(1, as.numeric(input$sim_fov_y %||% 4000))
    fov_z   <- max(1, as.numeric(input$sim_fov_z %||% 1000))

    plt <- plot_ly(source = "sim2d")

    # Helper: add a plain marker layer
    add_layer <- function(p, df, symbol, size, color, name, htxt = NULL) {
      if (nrow(df) == 0) return(p)
      z_val <- if (use_3d && "z" %in% names(df)) df$z else NULL
      args  <- list(
        x = df$x, y = df$y,
        type = tr_type, mode = "markers",
        marker = list(symbol = symbol, size = size, color = color,
                      line = list(color = color, width = 1)),
        name = name,
        hoverinfo = if (is.null(htxt)) "skip" else "text",
        text = htxt
      )
      if (!is.null(z_val)) args$z <- z_val
      do.call(add_trace, c(list(p), args))
    }

    # Build hover text helpers
    hover_xy <- function(df, prefix = "") {
      z_part <- if (use_3d && "z" %in% names(df))
                  paste0("<br>z: ", round(df$z, 1), " nm")
                else rep("", nrow(df))
      paste0(prefix,
             "x: ", round(df$x, 1), " nm",
             "<br>y: ", round(df$y, 1), " nm",
             z_part)
    }

    # GT layer — bright distinct markers so they stand out in both 2D and 3D.
    # In image mode these are the intensity-weighted sampled molecule positions.
    if (isTRUE(input$sim_show_gt) && nrow(res$gt) > 0) {
      gt_clean <- res$gt[, intersect(c("x","y","z"), names(res$gt)), drop = FALSE]
      gt_label <- if (identical(input$sim_shape, "image"))
                    "GT (sampled positions)" else "Ground truth"
      ht <- hover_xy(gt_clean, "GT<br>")
      plt <- add_layer(plt, gt_clean,
                       if (use_3d) "circle-open" else "cross",
                       round(pt_size * 2),
                       if (dark) "#FFEE44" else "#CC2200",
                       gt_label, ht)
    }


    # Localization layer — coloured by TID, Z, time, or local density.
    if (isTRUE(input$sim_show_locs) && nrow(res$locs) > 0) {
      cmode <- input$sim_color_mode %||% "tid"
      if (identical(cmode, "z") && !(use_3d && "z" %in% names(res$locs)))
        cmode <- "tid"                                  # no depth in 2D \u2192 TID
      z_val <- if (use_3d && "z" %in% names(res$locs)) res$locs$z else NULL
      z_ht  <- if (!is.null(z_val))
                  paste0("<br>z: ", round(z_val, 1), " nm")
               else rep("", nrow(res$locs))
      ht <- paste0("TID: ", res$locs$tid,
                   "<br>x: ", round(res$locs$x, 1), " nm",
                   "<br>y: ", round(res$locs$y, 1), " nm",
                   z_ht)

      if (identical(cmode, "tid")) {
        pal    <- tid_color_map(unique(res$locs$tid))
        marker <- list(size = pt_size,
                       color = unname(pal[res$locs$tid]), opacity = 0.85)
      } else {
        cv <- switch(cmode,
          z       = list(v = res$locs$z,   title = "Depth z (nm)"),
          tim     = list(v = res$locs$tim, title = "Time (s)"),
          locdens = list(v = sim_local_density(res$locs,
                              max(1, as.numeric(input$sim_dens_radius_nm %||% 500))),
                         title = "Local density"))
        v   <- as.numeric(cv$v)
        rng <- range(v[is.finite(v)])
        if (!all(is.finite(rng))) rng <- c(0, 1)
        marker <- list(
          size = pt_size, color = v, opacity = 0.85,
          colorscale = plotly_colorscale(
            input$sim_numeric_palette %||% "Spectral",
            reverse = isTRUE(input$sim_palette_reverse)),
          cmin = rng[1], cmax = rng[2], showscale = TRUE,
          colorbar = list(title = cv$title,
                          titlefont = list(color = fg_col),
                          tickfont  = list(color = fg_col)))
      }
      args <- list(
        x = res$locs$x, y = res$locs$y,
        type = tr_type, mode = "markers",
        marker = marker,
        name = sprintf("Localizations (n\u202F=\u202F%d)", nrow(res$locs)),
        hoverinfo = "text", text = ht
      )
      if (!is.null(z_val)) args$z <- z_val
      plt <- do.call(add_trace, c(list(plt), args))
    }

    # Lines: GT molecule position → median of each TID.
    #
    # One line per TID, anchored at the true GT molecule the TID came from.
    # For single-dye labels one short line per molecule shows the small
    # link offset; for polyclonal + secondary many longer lines fan out
    # from the same GT point, reflecting the greater linker distance and
    # the number of fluorophores.  Works in all shape modes including image.
    if (isTRUE(input$sim_show_gt_lines) &&
        nrow(res$gt) > 0 && nrow(res$fluoro) > 0 && nrow(res$locs) > 0 &&
        "mol_id" %in% names(res$fluoro)) {
      mol_ids  <- res$fluoro$mol_id
      has_gz   <- use_3d && "z" %in% names(res$gt)
      has_lz   <- use_3d && "z" %in% names(res$locs)
      seg_col  <- if (dark) "rgba(255,120,80,0.75)" else "rgba(190,50,30,0.75)"
      n_fluoro <- nrow(res$fluoro)

      # One segment per TID: GT[mol_id] → median(locs for that TID)
      segs <- do.call(rbind, lapply(seq_len(n_fluoro), function(k) {
        lk <- res$locs[res$locs$tid == as.character(k), , drop = FALSE]
        if (nrow(lk) == 0) return(NULL)
        gi <- mol_ids[k]
        if (is.na(gi) || gi < 1L || gi > nrow(res$gt)) return(NULL)
        data.frame(
          x0 = res$gt$x[gi],  y0 = res$gt$y[gi],
          x1 = median(lk$x),  y1 = median(lk$y),
          z0 = if (has_gz) res$gt$z[gi]  else NA_real_,
          z1 = if (has_lz) median(lk$z)  else NA_real_
        )
      }))

      if (!is.null(segs) && nrow(segs) > 0) {
        if (use_3d) {
          xs <- c(rbind(segs$x0, segs$x1, NA))
          ys <- c(rbind(segs$y0, segs$y1, NA))
          zs <- c(rbind(segs$z0, segs$z1, NA))
          plt <- do.call(add_trace,
                         list(plt, x = xs, y = ys, z = zs,
                              type = "scatter3d", mode = "lines",
                              line = list(color = seg_col, width = 2),
                              name = "GT \u2192 TID median",
                              hoverinfo = "skip", showlegend = TRUE))
        } else {
          plt <- plt %>%
            add_segments(x = ~x0, xend = ~x1, y = ~y0, yend = ~y1,
                         data = segs, inherit = FALSE,
                         line = list(color = seg_col, width = 1.5),
                         name = "GT \u2192 TID median",
                         hoverinfo = "skip", showlegend = TRUE)
        }
      }
    }

    title_txt <- sprintf("Simulated localizations%s",
                         if (use_3d) " (3D)" else "")

    if (use_3d) {
      ref_pts <- if (nrow(res$locs) > 0) res$locs else res$gt
      if (use_fov) {
        # FOV pins each axis; aspect ratio = true X:Y:Z so Z is compressed.
        xr <- c(-fov_x / 2, fov_x / 2)
        yr <- c(-fov_y / 2, fov_y / 2)
        zr <- c(-fov_z / 2, fov_z / 2)
        aspect <- list(x = 1, y = fov_y / fov_x, z = fov_z / fov_x)
      } else {
        # Equal-scale axes: Z gets the same span as XY so the plot looks cubic.
        # Localizations remain at their actual Z positions near the centre.
        x_rng <- range(ref_pts$x, na.rm = TRUE)
        y_rng <- range(ref_pts$y, na.rm = TRUE)
        half  <- max(diff(x_rng), diff(y_rng)) / 2 * 1.05   # 5% padding
        if (half < 1) half <- 1
        x_ctr <- mean(x_rng); y_ctr <- mean(y_rng)
        z_ctr <- if ("z" %in% names(ref_pts))
                    median(ref_pts$z, na.rm = TRUE) else 0
        xr <- c(x_ctr - half, x_ctr + half)
        yr <- c(y_ctr - half, y_ctr + half)
        zr <- c(z_ctr - half, z_ctr + half)
        aspect <- list(x = 1, y = 1, z = 1)
      }

      plt <- plt %>% layout(
        paper_bgcolor = bg_col,
        scene = list(
          bgcolor = bg_col,
          xaxis = list(title = "X (nm)", color = fg_col, range = xr,
                       showgrid = show_grid, zeroline = show_grid),
          yaxis = list(title = "Y (nm)", color = fg_col, range = yr,
                       showgrid = show_grid, zeroline = show_grid),
          zaxis = list(title = "Z (nm)", color = fg_col, range = zr,
                       showgrid = show_grid, zeroline = show_grid),
          aspectmode  = "manual",
          aspectratio = aspect
        ),
        font = list(color = fg_col), showlegend = TRUE,
        legend = list(x = 0.01, y = 0.99, xanchor = "left", yanchor = "top",
                      bgcolor = if (dark) "rgba(17,17,17,0.6)"
                                else "rgba(255,255,255,0.6)",
                      font = list(size = 11)),
        title = list(text = title_txt, font = list(color = fg_col))
      )
    } else {
      img_layout <- if (!is.null(res$img_b64)) {
        list(list(source  = res$img_b64,
                  xref    = "x", yref = "y",
                  x       = -res$img_w_nm / 2,
                  y       = res$img_h_nm / 2,
                  sizex   = res$img_w_nm, sizey = res$img_h_nm,
                  sizing  = "stretch", opacity = 0.25, layer = "below"))
      } else NULL

      xaxis <- list(title = "X (nm)", color = fg_col,
                    scaleanchor = "y", scaleratio = 1, gridcolor = gc,
                    showgrid = show_grid, zeroline = show_grid)
      yaxis <- list(title = "Y (nm)", color = fg_col, gridcolor = gc,
                    showgrid = show_grid, zeroline = show_grid)
      if (use_fov) {
        xaxis$range <- c(-fov_x / 2, fov_x / 2)
        yaxis$range <- c(-fov_y / 2, fov_y / 2)
      }

      # Scale bar (2D only) — reuses the shared SMLM make_scale_bar helper.
      sb <- NULL
      if (isTRUE(input$sim_show_scale_bar)) {
        x_span <- if (use_fov) fov_x else {
          ref <- if (nrow(res$locs) > 0) res$locs else res$gt
          diff(range(ref$x, na.rm = TRUE))
        }
        sb <- make_scale_bar(
          suppressWarnings(as.numeric(input$sim_scale_bar_length)),
          x_span,
          suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
          dark = dark,
          scale_font_size = suppressWarnings(
            as.numeric(input$sim_scale_font_size %||% 18)))
      }

      lay_args <- list(plt,
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font   = list(color = fg_col),
        xaxis  = xaxis, yaxis = yaxis,
        showlegend = TRUE,
        legend = list(x = 0.01, y = 0.99, xanchor = "left", yanchor = "top",
                      bgcolor = if (dark) "rgba(17,17,17,0.6)"
                                else "rgba(255,255,255,0.6)",
                      font = list(size = 11)),
        title  = list(text = title_txt, font = list(color = fg_col)),
        images = img_layout,
        uirevision = "sim"
      )
      if (!is.null(sb)) {
        lay_args$shapes      <- sb$shapes
        lay_args$annotations <- sb$annotations
      }
      plt <- do.call(layout, lay_args)
    }
    plt <- plt %>% config(displaylogo = FALSE,
                   toImageButtonOptions = list(format = "svg",
                                               filename = "sim_localisations"))
    # Register the relayout event so the "sim2d" source is recognised in both
    # 2D and 3D (avoids the "source not registered" console warning). Flag it as
    # rendered so the relayout observer only reads the event after registration.
    sim_plot_ready(TRUE)
    plt %>% plotly::event_register("plotly_relayout")
  })

  # Capture the visible x-span on zoom/pan and re-push the scale bar so it
  # tracks the current view (matches the SMLM tab's behaviour).
  observe({
    # Only read the relayout event when the plot has rendered + registered AND is
    # on screen (sim tab + MINFLUX subtab + 2D); otherwise plotly warns "source
    # not registered". sim_plot_ready stays FALSE until the first Generate render.
    if (!isTRUE(sim_plot_ready()) ||
        !identical(input$main_navbar, "tab_sim") ||
        !identical(input$sim_main_tabs, "MINFLUX") ||
        !identical(input$sim_view_mode, "2d")) return()
    ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "sim2d"))
    if (!is.list(ed)) return()
    x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
    x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
    if (length(x0) == 1 && length(x1) == 1 &&
        is.finite(x0) && is.finite(x1) && x1 > x0) {
      sim_axis_x_span(x1 - x0)
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      span <- sim_full_x_span()
      if (length(span) == 1 && is.finite(span) && span > 0) sim_axis_x_span(span)
    } else return()
    if (isTRUE(input$sim_show_scale_bar)) sim_push_scale_bar()
  })

  # Confocal view: capture its own zoom and re-push its (white) scale bar.
  observe({
    if (!isTRUE(sim_diff_ready()) ||
        !identical(input$main_navbar, "tab_sim") ||
        !identical(input$sim_main_tabs, "Confocal view")) return()
    ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "sim_diff"))
    if (!is.list(ed)) return()
    x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
    x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
    if (length(x0) == 1 && length(x1) == 1 &&
        is.finite(x0) && is.finite(x1) && x1 > x0) {
      sim_diff_x_span(x1 - x0)
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      sim_diff_x_span(NULL)
    } else return()
    if (isTRUE(input$sim_show_scale_bar)) sim_push_diff_scale_bar()
  })

  observeEvent(list(input$sim_show_scale_bar, input$sim_scale_bar_length,
                    input$sim_scale_bar_thickness, input$sim_scale_font_size,
                    input$sim_dark_mode, input$sim_view_mode,
                    input$sim_mode, input$sim_view_fov, input$sim_fov_x), {
    sim_push_scale_bar()
    sim_push_diff_scale_bar()
  }, ignoreInit = TRUE)

  # ── Phase 1 summary ──────────────────────────────────────────────────────────
  output$sim_summary <- renderUI({
    res <- sim_result()
    if (is.null(res)) return(NULL)
    nn_dist <- if (nrow(res$locs) >= 2) {
      idx  <- if (nrow(res$locs) > 500) sample(nrow(res$locs), 500)
              else seq_len(nrow(res$locs))
      cols <- intersect(c("x","y","z"), names(res$locs))
      pts  <- as.matrix(res$locs[idx, cols])
      d    <- as.matrix(dist(pts)); diag(d) <- Inf
      mean(apply(d, 1, min), na.rm = TRUE)
    } else NA_real_

    row <- function(lbl, val)
      tags$tr(tags$td(style = "padding:2px 8px;", lbl),
              tags$td(style = "padding:2px 8px;", tags$b(val)))

    # Steric warning block
    steric_ui <- if (!is.null(res$steric_warn)) {
      alert_cls <- if (identical(res$steric_warn$level, "danger"))
                     "alert alert-danger" else "alert alert-warning"
      tags$div(class = alert_cls,
               style = "padding:6px 10px; font-size:0.84em; margin-top:6px;",
               tags$b(if (identical(res$steric_warn$level, "danger"))
                        "\u26A0 Severe steric exclusion" else "\u26A0 Steric exclusion warning"),
               tags$br(), res$steric_warn$msg)
    } else NULL

    # Blinking TID inflation row (shown only when blinking is active)
    blink_row <- if (isTRUE(res$blink_active) && !is.null(res$n_fluoro_base)) {
      row("Blinking TID inflation",
          sprintf("%d base \u2192 %d apparent  (\u00D7%.1f, \u03BB\u202F=\u202F%.2f)",
                  res$n_fluoro_base, res$n_fluoro,
                  res$n_fluoro / max(1, res$n_fluoro_base),
                  res$blink_lambda_used))
    } else NULL

    tagList(
      steric_ui,
      tags$div(
        class = "well well-sm", style = "font-size:0.88em; margin-top:10px;",
        tags$b("Simulation summary"), tags$br(),
        tags$table(style = "width:100%; border-collapse:collapse;",
          row("Molecules placed",  res$n_prot),
          row("Labeling strategy", res$lbl_note),
          row("Fluorophores (trace IDs)",
              sprintf("%d  (~%.1f/molecule)", res$n_fluoro,
                      res$n_fluoro / max(1, res$n_prot))),
          blink_row,
          row("Linker distance",   sprintf("~%.0f nm", res$link_nm)),
          if (!is.na(res$footprint_nm))
            row("Label footprint (steric check)", sprintf("~%.0f nm", res$footprint_nm)),
          row("Total localizations", res$n_locs),
          row("Avg locs / trace ID",
              sprintf("%.1f", res$n_locs / max(1, res$n_fluoro))),
          row("Precision XY \u03C3",    sprintf("%.1f nm", res$sigma_loc_xy)),
          if (isTRUE(res$use_3d))
            row("Precision Z \u03C3",   sprintf("%.1f nm", res$sigma_loc_z)),
          if (is.finite(nn_dist))
            row("Mean nearest-neighbour", sprintf("%.1f nm", nn_dist)),
          if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
            row("Mean TID displacement error",
                sprintf("%.1f nm (SD %.1f, n\u202F=\u202F%d)",
                        mean(res$rmse_vals), sd(res$rmse_vals),
                        length(res$rmse_vals))),
          if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
            row("Theoretical combined \u03C3  (\u221A(link\u00B2 + \u03C3\u00B2_xy))",
                sprintf("%.1f nm", sqrt(res$link_nm^2 + res$sigma_loc_xy^2)))
        )
      )
    )
  })

  # ── Labeling strategy info card (live, no Generate needed) ─────────────────
  output$sim_label_info <- renderUI({
    key <- input$sim_label_type %||% "Custom"
    if (identical(key, "Custom")) return(NULL)
    p <- LABEL_PRESETS[[key]]
    if (is.null(p) || is.na(p$n_fluor)) return(NULL)
    tags$div(
      class = "alert alert-info",
      style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
      tags$b(sprintf("%d fluorophore(s) per molecule", p$n_fluor)),
      sprintf(" \u2014 linker ~%d nm", round(p$link_nm)),
      tags$br(), p$note
    )
  })

  # ── CRLB precision preview (live — no Generate needed) ──────────────────────
  output$sim_crlb_preview <- renderUI({
    req(identical(input$sim_prec_mode %||% "manual", "crlb"))
    L_nm    <- max(1,    as.numeric(input$sim_crlb_L       %||% 100))
    n_ph    <- max(1,    as.numeric(input$sim_crlb_photons %||% 50))
    bg_frac <- max(0, min(0.99, as.numeric(input$sim_crlb_bg %||% 0.05)))
    n_eff   <- max(1, n_ph * (1 - bg_frac))
    sig_xy  <- max(0.1, L_nm / (2 * sqrt(n_eff)))
    use_3d  <- identical(input$sim_view_mode, "3d")
    z_fac   <- max(0.1, as.numeric(input$sim_crlb_z_factor %||% 0.7))
    sig_z   <- sig_xy * z_fac
    tags$div(
      class = "alert alert-success",
      style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
      tags$b(sprintf("\u03C3_xy \u2248 %.2f nm", sig_xy)),
      if (use_3d) sprintf("   |   \u03C3_z \u2248 %.2f nm", sig_z),
      tags$br(),
      tags$small(sprintf(
        "L\u202F=\u202F%.0f nm,  N_eff\u202F=\u202F%.0f ph",
        L_nm, n_eff))
    )
  })

  # ── Phase 1 per-molecule TID count summary ───────────────────────────────────
  output$sim_tid_summary <- renderUI({
    res <- sim_result()
    if (is.null(res) || !("mol_id" %in% names(res$fluoro)) ||
        nrow(res$fluoro) == 0) return(NULL)

    # Count detected TIDs per molecule (background TIDs have mol_id = NA)
    real_tids  <- res$fluoro[!is.na(res$fluoro$mol_id), , drop = FALSE]
    if (nrow(real_tids) == 0) return(NULL)
    tids_per_mol <- as.integer(table(real_tids$mol_id))
    # Molecules that lost all fluorophores to detection efficiency = 0 TIDs
    n_zero  <- res$n_prot - length(tids_per_mol)
    counts  <- c(tids_per_mol, rep(0L, max(0L, n_zero)))

    row <- function(lbl, val)
      tags$tr(tags$td(style = "padding:2px 8px;", lbl),
              tags$td(style = "padding:2px 8px;", tags$b(val)))

    tagList(tags$div(
      class = "well well-sm",
      style = "font-size:0.88em; margin-top:6px;",
      tags$b("TIDs per molecule"),
      tags$table(style = "width:100%; border-collapse:collapse;",
        row("Mean",    sprintf("%.2f", mean(counts))),
        row("SD",      sprintf("%.2f", sd(counts))),
        row("Median",  sprintf("%d",   as.integer(median(counts)))),
        row("Range",   sprintf("%d \u2013 %d", min(counts), max(counts))),
        row("Molecules with 0 TIDs",
            sprintf("%d (%.0f%%)", n_zero,
                    100 * n_zero / max(1, res$n_prot))),
        if (res$n_bg > 0)
          row("Background TIDs added", res$n_bg)
      )
    ))
  })

  # ── Phase 1 precision budget ─────────────────────────────────────────────────
  # Bar chart breaking down: linker error | localization precision | combined
  # (+ observed RMSE when GT data are available).
  output$sim_precision_budget <- renderUI({
    res <- sim_result()
    if (is.null(res)) return(NULL)

    link       <- res$link_nm
    sig_xy     <- res$sigma_loc_xy
    prot_r     <- res$protein_r_nm %||% 0
    theory     <- sqrt(prot_r^2 + link^2 + sig_xy^2)
    obs        <- if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
                    mean(res$rmse_vals) else NA_real_

    # Resolution limit for separating two molecules: ~2× the combined
    # localization spread σ (Gaussian clouds separate at ≳2σ). MINFLUX
    # resolution is the measured spread (σ / MAD per TID), not an optical FWHM.
    resolve_limit <- 2 * theory
    max_val       <- max(theory, obs, 1, na.rm = TRUE) * 1.05

    make_bar <- function(label, value, color) {
      if (is.na(value) || value < 0) return(NULL)
      w_pct <- sprintf("%.1f%%", 100 * min(value, max_val) / max_val)
      tags$div(style = "margin-bottom:7px;",
        tags$div(
          style = paste0("font-size:0.80em; margin-bottom:2px;",
                         " display:flex; justify-content:space-between;"),
          tags$span(label),
          tags$b(sprintf("%.1f nm", value))
        ),
        tags$div(
          style = "background:#e8e8e8; border-radius:3px; height:11px;",
          tags$div(style = sprintf(
            "background:%s; width:%s; height:11px; border-radius:3px;",
            color, w_pct))
        )
      )
    }

    tags$div(
      class = "well well-sm",
      style = "font-size:0.88em; margin-top:8px;",
      tags$b("Precision Budget"),
      tags$br(),
      tags$small("Error components (\u03C3, nm):"),
      tags$br(), tags$br(),
      if (prot_r > 0) make_bar("Protein size (epitope offset)", prot_r, "#be185d"),
      make_bar("Linker / label error",           link,   "#e07b39"),
      make_bar("Localization precision \u03C3_xy", sig_xy, "#3b82c4"),
      make_bar("\u221A(link\u00B2 + \u03C3\u00B2_xy)  — combined theoretical",
               theory, "#7c3aed"),
      if (!is.na(obs))
        make_bar(sprintf("Observed mean TID displacement  (n\u202F=\u202F%d)",
                         length(res$rmse_vals)),
                 obs, "#16a34a"),
      tags$small(style = "color:#888;",
        "Observed should approach combined theoretical when traces are long."),
      if (!is.na(res$gt_nnd)) {
        resolvable <- res$gt_nnd >= resolve_limit
        tags$div(
          style = sprintf(
            paste0("margin-top:8px; padding:5px 8px; border-radius:4px;",
                   " font-size:0.82em; background:%s; color:%s;"),
            if (resolvable) "#dcfce7" else "#fef9c3",
            if (resolvable) "#166534" else "#854d0e"),
          tags$b(if (resolvable) "\u2714 Structure resolvable"
                 else "\u26A0 Below resolution limit"),
          tags$br(),
          sprintf("GT molecule spacing \u2248 %.0f nm vs 2\u03c3 limit = %.0f nm.  %s",
            res$gt_nnd, resolve_limit,
            if (resolvable) "Individual molecules can be resolved."
            else paste0("Proteins are closer than ~2\u03c3 of the localization",
                        " spread and will appear as unresolved clusters."))
        )
      }
    )
  })

  # ── Compare strategies: helper to build a stripped-down scatter plot ─────────
  make_cmp_plot <- function(gt, locs, use_3d, title,
                             loc_color = "rgba(220,38,38,0.65)",
                             dark = FALSE) {
    bg_col <- if (dark) "#111111" else "#FFFFFF"
    fg_col <- if (dark) "#DDDDDD" else "#444444"
    tr_type <- if (use_3d) "scatter3d" else "scattergl"
    show_grid <- !isTRUE(input$sim_hide_grid)

    # Compute isotropic axis ranges from GT so both panels share identical scales
    # and 1 nm maps to the same pixel distance on every axis.
    has_z  <- use_3d && "z" %in% names(gt)
    cx     <- mean(range(gt$x, na.rm = TRUE))
    cy     <- mean(range(gt$y, na.rm = TRUE))
    cz     <- if (has_z) mean(range(gt$z, na.rm = TRUE)) else 0
    spans  <- c(diff(range(gt$x, na.rm = TRUE)),
                diff(range(gt$y, na.rm = TRUE)),
                if (has_z) diff(range(gt$z, na.rm = TRUE)) else 0)
    half   <- max(spans, 1) / 2 * 1.15   # 15 % padding, same on all axes
    xrng   <- cx + c(-half, half)
    yrng   <- cy + c(-half, half)
    zrng   <- cz + c(-half, half)

    gt_col <- if (dark) "#FFFFFF" else "#000000"

    # Draw localizations first so GT markers render on top of them.
    plt <- plot_ly() %>%
      add_trace(type = tr_type,
                x = locs$x, y = locs$y,
                z = if (use_3d && "z" %in% names(locs)) locs$z else NULL,
                mode = "markers",
                marker = list(size = if (use_3d) 4 else 7,
                              color = loc_color,
                              opacity = 0.75),
                name = "Localizations",
                hoverinfo = "skip")

    # GT layer drawn on top so it is never hidden by localizations.
    # All shapes use open circle markers (no connecting lines).
    plt <- add_trace(plt,
              type = tr_type,
              x = gt$x, y = gt$y,
              z = if (has_z) gt$z else NULL,
              mode = "markers",
              marker = list(size   = if (use_3d) 10 else 14,
                            color  = gt_col,
                            symbol = "circle-open"),
              name = "GT molecules",
              hoverinfo = "skip")

    plt <- plt %>%
      layout(
        title = list(text = title, font = list(size = 12, color = fg_col)),
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        xaxis = list(color = fg_col, title = "x (nm)",
                     scaleanchor = "y", scaleratio = 1,
                     range = xrng, constrain = "domain", showgrid = show_grid, zeroline = show_grid),
        yaxis = list(color = fg_col, title = "y (nm)",
                     range = yrng, constrain = "domain", showgrid = show_grid, zeroline = show_grid),
        showlegend = TRUE,
        legend = list(font = list(size = 10, color = fg_col))
      )
    if (use_3d) plt <- plt %>%
      layout(scene = list(
        aspectmode  = "manual",
        aspectratio = list(x = 1, y = 1, z = 1),
        xaxis = list(title = "x (nm)", color = fg_col, range = xrng,
                     showgrid = show_grid, zeroline = show_grid),
        yaxis = list(title = "y (nm)", color = fg_col, range = yrng,
                     showgrid = show_grid, zeroline = show_grid),
        zaxis = list(title = "z (nm)", color = fg_col, range = zrng,
                     showgrid = show_grid, zeroline = show_grid),
        bgcolor = bg_col))
    plt
  }

  output$sim_cmp_plot_a <- renderPlotly({
    cmp <- compare_result()
    shiny::validate(shiny::need(!is.null(cmp), "Click Compare."))
    make_cmp_plot(cmp$gt, cmp$a$locs, cmp$use_3d,
                  "Strategy A", "rgba(220,38,38,0.75)",
                  dark = isTRUE(input$sim_dark_mode))
  })

  output$sim_cmp_plot_b <- renderPlotly({
    cmp <- compare_result()
    shiny::validate(shiny::need(!is.null(cmp), "Click Compare."))
    make_cmp_plot(cmp$gt, cmp$b$locs, cmp$use_3d,
                  "Strategy B", "rgba(37,99,235,0.75)",
                  dark = isTRUE(input$sim_dark_mode))
  })

  output$sim_cmp_table <- renderUI({
    cmp <- compare_result()
    if (is.null(cmp)) return(NULL)

    sig <- cmp$sigma_loc_xy
    a   <- cmp$a;  b <- cmp$b

    pr_a     <- a$protein_r_nm %||% 0
    pr_b     <- b$protein_r_nm %||% 0
    theory_a <- sqrt(pr_a^2 + a$link_nm^2 + sig^2)
    theory_b <- sqrt(pr_b^2 + b$link_nm^2 + sig^2)
    obs_a    <- if (!is.null(a$rmse_vals) && length(a$rmse_vals) >= 1)
                  sprintf("%.1f nm (n=%d)", mean(a$rmse_vals), length(a$rmse_vals))
                else "—"
    obs_b    <- if (!is.null(b$rmse_vals) && length(b$rmse_vals) >= 1)
                  sprintf("%.1f nm (n=%d)", mean(b$rmse_vals), length(b$rmse_vals))
                else "—"

    cell_style <- "padding:4px 8px; white-space:normal; word-break:break-word;"
    hdr  <- function(...) tags$th(style = paste0(cell_style, " text-align:left;"), ...)
    cel  <- function(...) tags$td(style = cell_style, ...)
    bcel <- function(...) tags$td(style = cell_style, tags$b(...))
    row2 <- function(lbl, va, vb)
      tags$tr(cel(lbl), bcel(va), bcel(vb))

    tags$div(
      class = "well well-sm",
      style = "font-size:0.86em; margin-top:10px; overflow-x:auto;",
      tags$b("Strategy comparison"),
      tags$table(
        style = "width:100%; border-collapse:collapse; margin-top:6px; table-layout:fixed;",
        tags$colgroup(
          tags$col(style = "width:36%;"),
          tags$col(style = "width:32%;"),
          tags$col(style = "width:32%;")
        ),
        tags$thead(tags$tr(
          hdr("Metric"),
          hdr(tags$span(style = "color:#dc2626;", "Strategy A")),
          hdr(tags$span(style = "color:#2563eb;", "Strategy B"))
        )),
        tags$tbody(
          row2("Label",
               a$lbl_note,
               b$lbl_note),
          row2("Fluorophores / molecule (mean)",
               sprintf("%.1f", a$n_fluoro / max(1, nrow(cmp$gt))),
               sprintf("%.1f", b$n_fluoro / max(1, nrow(cmp$gt)))),
          row2("Total TIDs",
               a$n_fluoro, b$n_fluoro),
          row2("Linker distance",
               sprintf("%.0f nm", a$link_nm),
               sprintf("%.0f nm", b$link_nm)),
          row2("Label footprint",
               if (!is.na(a$footprint_nm)) sprintf("~%.0f nm", a$footprint_nm) else "—",
               if (!is.na(b$footprint_nm)) sprintf("~%.0f nm", b$footprint_nm) else "—"),
          if (pr_a > 0 || pr_b > 0)
            row2("Protein radius",
                 sprintf("%.1f nm", pr_a),
                 sprintf("%.1f nm", pr_b)),
          row2("Combined σ",
               sprintf("%.1f nm", theory_a),
               sprintf("%.1f nm", theory_b)),
          row2("Observed TID displacement",
               obs_a, obs_b)
        )
      )
    )
  })

  # ── Confocal view: live PSF preview  ────────────────────
  output$sim_diff_psf_preview <- renderUI({
    na      <- max(0.5, min(1.5, as.numeric(input$sim_diff_NA %||% 1.4)))
    au      <- max(0.1, as.numeric(input$sim_diff_AU %||% 0.83))
    sig     <- compute_confocal_sigma(640, 670, na, au)
    fwhm    <- round(sig * 2.355)
    au_nm   <- round(1.22 * 670 / na)
    tags$div(
      class = "alert alert-info",
      style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
      tags$b(sprintf("\u03C3 \u2248 %.1f nm  (FWHM \u2248 %d nm)", sig, fwhm)),
      tags$br(),
      tags$small(sprintf(
        "1 AU = %d nm  \u2014  pinhole radius = %.0f nm",
        au_nm, au * au_nm / 2))
    )
  })

  # ── Confocal view plot ───────────────────────────────────────────────────────
  output$sim_diff_plot <- renderPlotly({
    res <- sim_result()
    shiny::validate(
      shiny::need(!is.null(res) && nrow(res$gt) > 0,
                  "Click 'Generate' first to build the GT positions."))

    dark   <- isTRUE(input$sim_dark_mode)
    bg_col <- if (dark) "#111111" else "#FFFFFF"
    fg_col <- if (dark) "#DDDDDD" else "#444444"
    gc     <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"

    na      <- max(0.5, min(1.5, as.numeric(input$sim_diff_NA %||% 1.4)))
    au      <- max(0.1, as.numeric(input$sim_diff_AU %||% 0.83))
    psf_sig <- compute_confocal_sigma(640, 670, na, au)
    n_ph    <- max(1, as.numeric(input$sim_diff_photons %||% 200))
    bg_ph   <- max(0, as.numeric(input$sim_diff_bg      %||% 5))

    # Render from the UNIQUE fluorophore pool (one physical dye = one PSF), not
    # the over-counted MINFLUX detection list — re-detecting a molecule over the
    # acquisition does not make its confocal spot brighter.
    #   • Each fluorophore = one PSF at its displaced position, n_ph photons.
    #     Multi-dye labels (e.g. polyclonal + secondary) still stack into a
    #     brighter spot because they are distinct dyes.
    #   • Background TIDs (non-specific binding) render as single n_ph spots.
    fp     <- if (!is.null(res$fluoro_pool)) res$fluoro_pool else res$fluoro
    src_x  <- if (nrow(fp) > 0) fp$x else numeric(0)
    src_y  <- if (nrow(fp) > 0) fp$y else numeric(0)
    src_ph <- rep(n_ph, length(src_x))

    if (!is.null(res$bg_pos) && nrow(res$bg_pos) > 0) {
      src_x  <- c(src_x,  res$bg_pos$x)
      src_y  <- c(src_y,  res$bg_pos$y)
      src_ph <- c(src_ph, rep(n_ph, nrow(res$bg_pos)))
    }

    # Cap at 5 000 sources total to keep rendering fast (sample uniformly)
    if (length(src_x) > 5000L) {
      idx    <- sample(length(src_x), 5000L)
      src_x  <- src_x[idx]; src_y <- src_y[idx]; src_ph <- src_ph[idx]
    }

    shiny::validate(shiny::need(length(src_x) > 0,
      "No fluorophores to render. Adjust labeling efficiency or generate first."))

    px_user  <- max(1, as.numeric(input$sim_diff_px %||% 50))
    # Pin the confocal window to the FOV when it's active so the confocal spans
    # the same ROI as the MINFLUX view rather than cropping to the sources.
    use_fov  <- sim_use_fov()
    fov_ext  <- if (use_fov) {
                  fx <- max(1, as.numeric(input$sim_fov_x %||% 4000))
                  fy <- max(1, as.numeric(input$sim_fov_y %||% 4000))
                  c(-fx / 2, fx / 2, -fy / 2, fy / 2)
                } else NULL
    img_data <- render_diffraction_image(src_x, src_y, psf_sig, n_ph, bg_ph,
                                          px_nm_user = px_user,
                                          photons_vec = src_ph,
                                          extent = fov_ext)
    shiny::validate(shiny::need(!is.null(img_data), "Rendering failed."))

    # FWHM = 2.355 × σ: optical resolution — apparent width of a point source
    # and the minimum resolvable separation between two identical emitters.
    fwhm_nm   <- round(psf_sig * 2.355)
    title_txt <- sprintf(
      "Confocal view  (\u03C3\u202F=\u202F%.0f nm,  opt. resolution FWHM\u202F\u2248\u202F%d nm)",
      psf_sig, fwhm_nm)

    # Hot Red LUT (ImageJ "Red Hot") — fixed regardless of dark mode; the
    # confocal image is always rendered on black. Low → high:
    # black → red → orange → yellow → white.
    cscale <- list(
      list(0.00, "rgb(0,0,0)"),
      list(0.25, "rgb(128,0,0)"),
      list(0.45, "rgb(255,0,0)"),
      list(0.65, "rgb(255,140,0)"),
      list(0.85, "rgb(255,255,0)"),
      list(1.00, "rgb(255,255,255)"))

    # The confocal image is inherently a black-background picture, so render it
    # on black with no gridlines and white axis text regardless of dark mode.
    cf_bg <- "#000000"; cf_fg <- "#DDDDDD"

    # Scale bar — same controls as the MINFLUX view, forced white so it stands
    # out on the black image, measured against the rendered image width.  The
    # full image span is stored so the zoom observer can re-scale the bar.
    x_span <- diff(range(img_data$xs, na.rm = TRUE))
    sim_diff_full(x_span)
    sb <- NULL
    if (isTRUE(input$sim_show_scale_bar)) {
      sb <- make_scale_bar(
        suppressWarnings(as.numeric(input$sim_scale_bar_length)),
        x_span,
        suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
        dark = TRUE,
        scale_font_size = suppressWarnings(
          as.numeric(input$sim_scale_font_size %||% 18)))
    }

    lay_args <- list(
      paper_bgcolor = cf_bg, plot_bgcolor = cf_bg,
      font  = list(color = cf_fg),
      xaxis = list(title = "X (nm)", color = cf_fg,
                   scaleanchor = "y", scaleratio = 1,
                   showgrid = FALSE, zeroline = FALSE),
      yaxis = list(title = "Y (nm)", color = cf_fg,
                   showgrid = FALSE, zeroline = FALSE),
      title = list(text = title_txt, font = list(color = cf_fg)))
    if (!is.null(sb)) {
      lay_args$shapes      <- sb$shapes
      lay_args$annotations <- sb$annotations
    }

    fig <- plot_ly(source = "sim_diff") %>%
      add_trace(
        type       = "heatmap",
        x          = img_data$xs, y = img_data$ys, z = img_data$img,
        colorscale = cscale,
        showscale  = FALSE, hoverinfo = "skip",
        name       = "Diffraction-limited"
      )
    sim_diff_ready(TRUE)
    do.call(layout, c(list(fig), lay_args)) %>%
      plotly::event_register("plotly_relayout") %>%
      config(displaylogo = FALSE,
             toImageButtonOptions = list(format = "svg",
                                         filename = "diffraction_limited"))
  })

  # ── Phase 1 CSV export ───────────────────────────────────────────────────────
  # Each row is one localisation.  When GT data exist, mol_id and the GT
  # coordinates (gt_x, gt_y, gt_z) are appended so every loc can be traced
  # back to its source molecule.  Background TIDs get NA in these columns.
  output$sim_dl_csv <- downloadHandler(
    filename = function()
      paste0("sim_localisations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- sim_result()
      shiny::validate(shiny::need(!is.null(res) && nrow(res$locs) > 0,
                                  "No localizations to export."))
      out <- res$locs[, intersect(c("x", "y", "z", "tid", "tim"),
                                  names(res$locs)), drop = FALSE]

      if ("mol_id" %in% names(res$fluoro) && nrow(res$gt) > 0) {
        n_real  <- nrow(res$fluoro)
        tid_idx <- suppressWarnings(as.integer(out$tid))
        valid   <- !is.na(tid_idx) & tid_idx >= 1L & tid_idx <= n_real
        mol_ids <- rep(NA_integer_, nrow(out))
        mol_ids[valid] <- res$fluoro$mol_id[tid_idx[valid]]
        out$mol_id <- mol_ids
        out$gt_x   <- ifelse(is.na(mol_ids), NA_real_, res$gt$x[mol_ids])
        out$gt_y   <- ifelse(is.na(mol_ids), NA_real_, res$gt$y[mol_ids])
        if ("z" %in% names(res$gt))
          out$gt_z <- ifelse(is.na(mol_ids), NA_real_, res$gt$z[mol_ids])
      }

      write.csv(out, file, row.names = FALSE)
    }
  )

  # ── Acquisition timeline: resolution-over-time plot ──────────────────────────
  output$sim_acq_plot <- renderPlotly({
    shiny::validate(shiny::need(
      identical(input$sim_mode %||% "normal", "acquisition"),
      "Switch to Acquisition mode (top of the sidebar) to model an acquisition timeline."))
    a <- sim_acq_live()
    shiny::validate(shiny::need(
      !is.null(a),
      "Click 'Generate' to build the acquisition timeline."))
    dark   <- isTRUE(input$sim_dark_mode)
    bg_col <- if (dark) "#111111" else "#FFFFFF"
    fg_col <- if (dark) "#DDDDDD" else "#444444"
    gc     <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"

    # Work in minutes; lead with the detection rate (traces/hour) since that is
    # the curve directly comparable to a real run (early peak → constant floor).
    rt_min  <- a$rate_t / 60
    t_min   <- a$t / 60
    # % resolved is relative to the LABELED molecules, so the curve plateaus at
    # the fraction actually reachable in this run (reaches 100 % when all
    # labeled molecules are found).
    cum_pct <- 100 * a$cum_unique / max(1, a$n_mol_resolvable)
    rmax    <- max(a$rate_v, 1, na.rm = TRUE)

    p <- plot_ly()
    if (!is.null(a$rate_t))
      p <- p %>% add_trace(x = rt_min, y = a$rate_v, type = "scatter",
                           mode = "lines", fill = "tozeroy",
                           fillcolor = "rgba(224,123,57,0.22)",
                           line = list(color = "#e07b39", width = 1.5),
                           name = "Detection rate", hoverinfo = "x+y")
    # Optional decomposition of the PREDICTED rate into its structure (decaying +
    # bleaching) and constant background components — makes explicit that the
    # orange total = structure + background, and the floor is background.
    # Convention: solid = model, dotted = measured (real run, below).
    if (isTRUE(input$sim_real_show_split) && !is.null(a$rate_struct)) {
      p <- p %>%
        add_trace(x = rt_min, y = a$rate_struct, type = "scatter",
                  mode = "lines",
                  line = list(color = "#2a9d8f", width = 1.5),
                  name = "Structure (model)", hoverinfo = "x+y") %>%
        add_segments(x = 0, xend = max(rt_min),
                     y = a$bg_rate, yend = a$bg_rate,
                     line = list(color = "#b07aa1", width = 1.5),
                     name = "Background (model)", inherit = FALSE)
    }
    # Late-phase floor line (median of the last quarter of the run).
    if (!is.null(a$rate_floor) && is.finite(a$rate_floor))
      p <- p %>% add_segments(x = 0, xend = max(rt_min),
                              y = a$rate_floor, yend = a$rate_floor,
                              line = list(color = "#e07b39", width = 1,
                                          dash = "dot"),
                              name = sprintf("floor ~%.0f/hr", a$rate_floor),
                              inherit = FALSE)
    # Cumulative molecules resolved (%) on the secondary axis.
    p <- p %>% add_trace(x = t_min, y = cum_pct, type = "scatter",
                         mode = "lines", yaxis = "y2",
                         line = list(color = "#3b82c4", width = 2),
                         name = "Molecules resolved", hoverinfo = "x+y")
    # Measured detection rate from a loaded real run (calibration overlay).
    # Optionally decomposed into the spatially-split structure / background.
    real <- sim_real_rv()
    if (!is.null(real)) {
      if (isTRUE(input$sim_real_show_split) && !is.null(real$rate_v_struct)) {
        p <- p %>%
          add_trace(x = real$rate_t / 60, y = real$rate_v_struct,
                    type = "scatter", mode = "lines",
                    line = list(color = "#2a9d8f", width = 1.5, dash = "dot"),
                    name = "Measured structure", hoverinfo = "x+y") %>%
          add_trace(x = real$rate_t / 60, y = real$rate_v_bg,
                    type = "scatter", mode = "lines",
                    line = list(color = "#b07aa1", width = 1.5, dash = "dot"),
                    name = "Measured background", hoverinfo = "x+y")
      } else {
        p <- p %>% add_trace(x = real$rate_t / 60, y = real$rate_v,
                             type = "scatter", mode = "lines",
                             line = list(color = "#777777", width = 1.5,
                                         dash = "dot"),
                             name = "Measured (real run)", hoverinfo = "x+y")
      }
    }
    # Resolution marker only when the threshold is actually reached.
    if (!is.null(a$res_time) && is.finite(a$res_time))
      p <- p %>% add_segments(x = a$res_time / 60, xend = a$res_time / 60,
                              y = 0, yend = rmax,
                              line = list(color = "#16a34a", width = 1.5,
                                          dash = "dash"),
                              name = sprintf("%.0f%% resolved", 100 * a$resolve_q),
                              inherit = FALSE)
    p %>%
      layout(
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        font   = list(color = fg_col),
        margin = list(l = 60, r = 95, t = 50, b = 110),
        xaxis  = list(title = "Acquisition time (min)", color = fg_col,
                      gridcolor = gc),
        yaxis  = list(title = "Detection rate (traces/hour)", color = "#e07b39",
                      gridcolor = gc, rangemode = "tozero", automargin = TRUE),
        yaxis2 = list(title = list(text = "Molecules resolved (%)",
                                   standoff = 22),
                      color = "#3b82c4", overlaying = "y", side = "right",
                      range = c(0, 105), showgrid = FALSE, automargin = TRUE),
        title  = list(text = "Detection rate & resolution over time",
                      font = list(color = fg_col), x = 0.02, xanchor = "left",
                      y = 0.97, yanchor = "top"),
        showlegend = TRUE,
        legend = list(orientation = "h", x = 1, xanchor = "right",
                      y = -0.28, yanchor = "top", font = list(size = 10))
      ) %>%
      config(displaylogo = FALSE)
  })

  # ── Acquisition timeline summary ─────────────────────────────────────────────
  output$sim_acq_summary <- renderUI({
    a <- sim_acq_live()
    if (is.null(a)) return(NULL)

    fmt_t <- function(s) if (s >= 60) sprintf("%.1f min", s / 60)
                         else          sprintf("%.1f s", s)
    reached <- !is.null(a$res_time) && is.finite(a$res_time)
    capped  <- !is.null(a$struct_cap) && a$struct_cap > 0 &&
               a$n_real_tid >= a$struct_cap

    row <- function(lbl, val)
      tags$tr(tags$td(style = "padding:2px 8px;", lbl),
              tags$td(style = "padding:2px 8px;", tags$b(val)))

    tags$div(
      class = "well well-sm", style = "font-size:0.88em; margin-top:10px;",
      tags$b("Acquisition timeline"), tags$br(),
      tags$table(style = "width:100%; border-collapse:collapse;",
        row("Total acquisition time", fmt_t(a$T_sec)),
        row("Structure detection rate",
            sprintf("%.0f traces/hour  (mean %.0f)", a$struct_rate, a$rate_mean)),
        row("Structure trace IDs",
            if (capped) sprintf("%d  (photon budget reached)", a$n_real_tid)
            else        sprintf("%d", a$n_real_tid)),
        row("Molecules resolved",
            sprintf("%d / %d labeled  (%.0f%%; %d of %d total)",
                    a$n_mol_det, a$n_mol_resolvable, 100 * a$resolved_frac,
                    a$n_mol_det, a$n_mol_total)),
        row(sprintf("Time to %.0f%% resolved", 100 * a$resolve_q),
            if (reached) fmt_t(a$res_time)
            else sprintf("not reached (%.0f%% max)", 100 * a$resolved_frac)),
        if (reached)
          row("Time after resolved",
              sprintf("%s  (%.0f%% of run, mostly re-detecting)",
                      fmt_t(max(0, a$T_sec - a$res_time)),
                      100 * max(0, a$T_sec - a$res_time) / max(1e-9, a$T_sec))),
        row("Detections per molecule",
            sprintf("%.1f× on average  (%.0f%% of traces are repeats)",
                    a$det_per_mol, 100 * a$overcount_frac)),
        if (!is.null(a$n_bg) && a$n_bg > 0)
          row("Background TIDs",
              sprintf("%d  (@ %.0f/hr)", a$n_bg, a$bg_rate))
      ),
      tags$small(style = "color:#888;",
        if (reached)
          paste("Resolved early, then later traces re-detect the same",
                "molecules — the over-counting you group as repeat blinks.",
                "Increase the rate or duration to resolve more; decrease them",
                "(or shorten the run) to see partial resolution.")
        else
          paste(sprintf("Only %.0f%% of the structure resolved.",
                        100 * a$resolved_frac),
                "Increase the detection rate or duration (or improve labeling)",
                "to reach the resolve threshold."))
    )
  })

  # ── Phase 1 GT positions export ──────────────────────────────────────────────
  output$sim_dl_gt_csv <- downloadHandler(
    filename = function()
      paste0("sim_gt_positions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- sim_result()
      shiny::validate(shiny::need(!is.null(res) && nrow(res$gt) > 0,
                                  "No ground-truth positions to export."))
      cols <- intersect(c("x","y","z"), names(res$gt))
      gt_out <- res$gt[, cols, drop = FALSE]
      gt_out$mol_id <- seq_len(nrow(gt_out))
      write.csv(gt_out, file, row.names = FALSE)
    }
  )

  # ══════════════════════════════════════════════════════════════════════════
  # Diffusion / Tracking ground-truth simulator (Phase 1)
  # See docs/simulation_module.md. Free surface Brownian on a sphere or rod,
  # sampled the MINFLUX way (cycle clock + dark cycles + heteroscedastic σ_loc +
  # photobleaching) → MINFLUX-format tracks + a ground-truth D. The CVE round-trip
  # in the summary demonstrates the validation loop in-place.
  # ══════════════════════════════════════════════════════════════════════════
  # --- Calibrate the sampling model from a real MINFLUX file -----------------
  # Reads a real .csv/.npy/.json and auto-fills cycle time, photons, dark-cycle
  # probability and mean track length so the simulator matches a real experiment.
  sim_trk_calib_rv <- reactiveVal(NULL)
  observeEvent(input$sim_trk_calib, {
    fi <- input$sim_trk_calib_file
    shiny::validate(shiny::need(!is.null(fi), "Choose a MINFLUX file first."))
    p <- fi$datapath; nm <- tolower(fi$name)
    df <- withProgress(message = "Reading MINFLUX file…", value = 0, {
      tryCatch({
        if (grepl("\\.npy$",  nm)) read_minflux_npy_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
        else if (grepl("\\.json$", nm)) read_minflux_json_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
        else if (grepl("\\.csv$",  nm)) {
          x <- as.data.frame(data.table::fread(p, showProgress = FALSE))
          if (all(c("tid","tim") %in% names(x))) x else canonicalize_minflux_df(x)
        } else NULL
      }, error = function(e) NULL)
    })
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0 && all(c("tid","tim") %in% names(df)),
      "Could not read tid / tim from that file."))
    df <- df[is.finite(df$tim), ]
    # Calibrate ONLY on analysable traces (≥ MIN_LEN localizations, no upper bound), so
    # the short-junk traces don't drag the statistics down. The representative length is
    # the MEDIAN of the kept traces.
    MIN_LEN  <- max(2L, as.integer(input$sim_trk_calib_minlen %||% 100L))
    len_all  <- table(df$tid)
    keep_tid <- names(len_all)[len_all >= MIN_LEN]
    shiny::validate(shiny::need(length(keep_tid) >= 5, sprintf("Fewer than 5 traces have ≥%d localizations.", MIN_LEN)))
    n_total  <- length(len_all)
    df <- df[as.character(df$tid) %in% keep_tid, ]        # keep only the long traces
    len     <- as.numeric(table(df$tid)); med_len <- stats::median(len)
    photons <- if ("eco" %in% names(df)) stats::median(df$eco[is.finite(df$eco) & df$eco > 0]) else NA_real_
    o    <- order(df$tid, df$tim); tid <- df$tid[o]; tim <- df$tim[o]; nrw <- length(tid)
    same <- c(FALSE, tid[-1] == tid[-nrw]); dt <- c(NA_real_, diff(tim)) * 1000; dt[!same] <- NA_real_  # ms
    dtv  <- dt[is.finite(dt) & dt > 0]
    shiny::validate(shiny::need(length(dtv) > 20, "Too few within-track steps to calibrate."))
    hb   <- hist(dtv[dtv < stats::quantile(dtv, 0.98)], breaks = 400, plot = FALSE)
    cyc  <- hb$mids[which.max(hb$counts)]                          # modal Δt = cycle time
    ncyc <- pmax(1L, as.integer(round(dtv / cyc)))
    p_dark <- sum(ncyc - 1L) / sum(ncyc); integ <- mean(ncyc == 1L)
    if (is.finite(cyc))     updateNumericInput(session, "sim_trk_cycle",   value = signif(cyc, 3))
    if (is.finite(p_dark))  updateNumericInput(session, "sim_trk_pdark",   value = round(p_dark, 2))
    if (is.finite(photons)) updateNumericInput(session, "sim_trk_photons", value = round(photons))
    if (is.finite(med_len)) updateNumericInput(session, "sim_trk_meanlen", value = round(med_len))
    sim_trk_calib_rv(list(file = fi$name, n_tracks = length(keep_tid), n_total = n_total,
      min_len = MIN_LEN, n_loc = nrow(df), med_len = med_len, mean_len = mean(len),
      photons = photons, cyc = cyc, p_dark = p_dark, integ = integ, med_dt = stats::median(dtv)))
    showNotification(sprintf("Calibrated on %d traces (≥%d locs) of %d total.",
                             length(keep_tid), MIN_LEN, n_total), type = "message", duration = 6)
  })

  output$sim_trk_calib_report <- renderUI({
    cb <- sim_trk_calib_rv(); if (is.null(cb)) return(NULL)
    rowf <- function(k, v) htmltools::tags$tr(
      htmltools::tags$td(style = "padding:1px 10px 1px 0; color:#666;", k),
      htmltools::tags$td(htmltools::tags$b(v)))
    htmltools::div(style = "font-size:0.8em; color:#333; margin-top:4px;",
      htmltools::tags$div(style = "color:#228833;", sprintf("✓ Calibrated from %s", cb$file)),
      htmltools::tags$table(style = "border-collapse:collapse;",
        rowf(sprintf("Traces used (≥%d locs)", cb$min_len %||% 50L),
             sprintf("%s of %s  ·  %s locs", format(cb$n_tracks, big.mark = ","),
                     format(cb$n_total %||% cb$n_tracks, big.mark = ","), format(cb$n_loc, big.mark = ","))),
        rowf("Cycle time", sprintf("%.3g ms", cb$cyc)),
        rowf("Photons/loc (median eco)", if (is.finite(cb$photons)) sprintf("%.0f", cb$photons) else "—"),
        rowf("Dark-cycle prob.", sprintf("%.2f  (%.0f%% integrated)", cb$p_dark, 100 * cb$integ)),
        rowf("Track length (median)", sprintf("%.0f · mean %.0f", cb$med_len, cb$mean_len))),
      htmltools::tags$p(style = "color:#888; margin-top:3px;",
        "Auto-filled the sampling fields (mean locs/track uses the robust median). Set geometry, D and ",
        "the effective TCP L, then Simulate."))
  })

  sim_trk_result <- eventReactive(input$sim_trk_run, {
    set.seed(as.integer(input$sim_trk_seed %||% 1L))
    geom     <- input$sim_trk_geom %||% "rod"
    mode     <- input$sim_trk_mode %||% "surface"      # surface | volume (intracellular)
    d_mot    <- if (identical(mode, "volume")) 3L else 2L   # motion DOF (used by the CVE check)
    D_um2_s  <- max(0, as.numeric(input$sim_trk_D %||% 0.05))
    D_nm2ms  <- D_um2_s * 1000                        # µm²/s → nm²/ms
    n_trk    <- max(1L, as.integer(input$sim_trk_ntracks %||% 200L))
    mean_len <- max(5L, as.integer(input$sim_trk_meanlen %||% 80L))
    cyc      <- max(0.01, as.numeric(input$sim_trk_cycle %||% 0.5))     # ms
    p_dark   <- min(0.95, max(0, as.numeric(input$sim_trk_pdark %||% 0.2)))
    photons  <- max(1, as.numeric(input$sim_trk_photons %||% 100))
    L_tcp    <- max(1, as.numeric(input$sim_trk_Ltcp %||% 150))
    proj2d   <- isTRUE(input$sim_trk_project2d)
    sd_step  <- sqrt(2 * D_nm2ms * cyc)               # per-axis step SD per cycle (nm)
    R    <- max(1, as.numeric(input$sim_trk_R %||% 500))
    rr   <- max(1, as.numeric(input$sim_trk_r %||% 400))
    Lcyl <- max(1, as.numeric(input$sim_trk_Lcyl %||% 2000))
    circ <- 2 * pi * rr
    ea   <- max(1, as.numeric(input$sim_trk_ea %||% 500))
    eb   <- max(1, as.numeric(input$sim_trk_eb %||% 500))
    ec   <- max(1, as.numeric(input$sim_trk_ec %||% 900))
    abc  <- if (identical(geom, "sphere")) c(R, R, R) else c(ea, eb, ec)   # ellipsoid semi-axes
    # Long-axis orientation for the rod: E. coli on a coverslip lie flat → axis in the
    # image plane. Generated in canonical coords (axis ∥ X) then permuted to the choice.
    rod_axis <- input$sim_trk_rod_axis %||% "x"
    rod_perm <- switch(rod_axis, x = c(1L,2L,3L), y = c(2L,1L,3L), z = c(2L,3L,1L))

    # ── unified geometry operators (rod = spherocylinder body + hemispherical caps;
    #    else ellipsoid, sphere = a=b=c). project = nearest surface point; reflect =
    #    back inside the volume; normal = unit outward normal (for the tangent step). ─
    is_rod <- identical(geom, "rod")
    geo_normal <- function(p) {
      if (is_rod) { a <- p[1]
        n <- if (a >= 0 && a <= Lcyl) c(0, p[2], p[3]) else { cx <- if (a < 0) 0 else Lcyl; p - c(cx,0,0) }
        m <- sqrt(sum(n^2)); if (m == 0) c(0,1,0) else n / m
      } else { n <- p / abc^2; n / sqrt(sum(n^2)) }
    }
    geo_project <- function(p) {
      if (is_rod) { a <- p[1]
        if (a >= 0 && a <= Lcyl) { rad <- c(p[2],p[3]); rho <- sqrt(sum(rad^2)); rad <- if (rho == 0) c(rr,0) else rad*rr/rho; c(a, rad) }
        else { cx <- if (a < 0) 0 else Lcyl; d <- p - c(cx,0,0); dd <- sqrt(sum(d^2)); d <- if (dd == 0) c(if (a<0) -rr else rr,0,0) else d*rr/dd; c(cx,0,0)+d }
      } else { u <- p/abc; abc * u / sqrt(sum(u^2)) }
    }
    geo_reflect <- function(p) {
      if (is_rod) { a <- p[1]
        if (a >= 0 && a <= Lcyl) { rad <- c(p[2],p[3]); rho <- sqrt(sum(rad^2)); if (rho > rr) rad <- rad*(2*rr - rho)/rho; c(a, rad) }
        else { cx <- if (a < 0) 0 else Lcyl; d <- p - c(cx,0,0); dd <- sqrt(sum(d^2)); if (dd > rr) d <- d*(2*rr - dd)/dd; c(cx,0,0)+d }
      } else { u <- p/abc; nu <- sqrt(sum(u^2)); if (nu > 1) u <- u*(2-nu)/nu; abc*u }
    }
    geo_rand_surf <- function() if (is_rod) { s <- runif(1,0,Lcyl); phi <- runif(1,0,2*pi); c(s, rr*cos(phi), rr*sin(phi)) } else geo_project(rnorm(3) * abc)
    geo_rand_in   <- function() if (is_rod) {
        repeat { p <- c(runif(1,-rr,Lcyl+rr), runif(1,-rr,rr), runif(1,-rr,rr)); a <- p[1]
          if ((if (a >= 0 && a <= Lcyl) (p[2]^2+p[3]^2) <= rr^2 else if (a < 0) sum(p^2) <= rr^2 else sum((p-c(Lcyl,0,0))^2) <= rr^2)) return(p) }
      } else { d <- rnorm(3); d <- d/sqrt(sum(d^2)); abc * d * runif(1)^(1/3) }

    # ── motion modifier (Phase 3): free | confined (OU tether) | directed (drift) |
    #    twostate (mobile↔immobile CTMC). Applied to the proposed 3D point, then the
    #    geometry constraint re-projects/reflects it. ──────────────────────────────
    motion  <- input$sim_trk_motion  %||% "free"
    L_c     <- max(1, as.numeric(input$sim_trk_Lc %||% 150))         # confinement radius (nm)
    v_nmms  <- max(0, as.numeric(input$sim_trk_vel %||% 0))          # directed speed (µm/s = nm/ms)
    kswitch <- min(1, max(0, as.numeric(input$sim_trk_kswitch %||% 0.03)))  # per-cycle state-switch prob
    tau_c   <- if (D_nm2ms > 0) L_c^2 / D_nm2ms else Inf             # OU relaxation time (ms)
    phi_c   <- exp(-cyc / tau_c)                                     # per-cycle mean-reversion factor
    drift0  <- if (is_rod) c(1,0,0) else { v <- rnorm(3); v/sqrt(sum(v^2)) }   # rod drift ∥ long axis

    gen_path <- function(n_cyc) {
      surface <- identical(mode, "surface")
      p <- if (surface) geo_rand_surf() else geo_rand_in(); home <- p; st <- 1L
      pos <- matrix(NA_real_, n_cyc, 3); pos[1,] <- p; states <- integer(n_cyc); states[1] <- 1L
      for (i in 2:n_cyc) {
        Dstep <- if (identical(motion,"twostate") && st == 0L) D_nm2ms * 0.01 else D_nm2ms
        sdi   <- sqrt(2 * Dstep * cyc)
        if (surface) { g <- rnorm(3,0,sdi); nn <- geo_normal(p); p2 <- p + (g - sum(g*nn)*nn) }  # tangent step
        else         { p2 <- p + rnorm(3,0,sdi) }
        if (identical(motion,"directed")) p2 <- p2 + drift0 * v_nmms * cyc
        if (identical(motion,"confined")) p2 <- home + (p2 - home) * phi_c                       # OU pull to home
        p  <- if (surface) geo_project(p2) else geo_reflect(p2)
        if (identical(motion,"twostate") && runif(1) < kswitch) st <- 1L - st
        pos[i,] <- p; states[i] <- st
      }
      if (is_rod) pos <- pos[, rod_perm, drop = FALSE]
      attr(pos, "states") <- states
      pos
    }
    # Survivorship: a molecule that moves > R_reacq during a dark gap is lost — the
    # track ends there. Reproduces the dark-cycle bias (long jumps culled → apparent
    # confinement). Off unless enabled.
    surv    <- isTRUE(input$sim_trk_surv)
    R_reacq <- max(1, as.numeric(input$sim_trk_rreacq %||% 300))

    rows <- vector("list", n_trk); gt <- vector("list", n_trk); k <- 0L
    withProgress(message = "Simulating tracks…", value = 0, {
    for (ti in seq_len(n_trk)) {
      if (ti %% 20L == 0L) incProgress(20 / n_trk)
      # mean_len = target MEDIAN recorded localizations/track. Recorded count is ≈
      # geometric, whose median ≈ ln(2)·mean, and only (1−p_dark) of cycles are recorded,
      # so the bleach lifetime in CYCLES has mean = median / (ln2 · (1−p_dark)) to land
      # the recorded MEDIAN on mean_len despite dark-cycling.
      n_cyc <- max(3L, stats::rgeom(1, 1 / (mean_len / (log(2) * max(0.05, 1 - p_dark)))) + 1L)
      pos   <- gen_path(n_cyc)                          # ground-truth path at cycle resolution
      integ <- runif(n_cyc) > p_dark; integ[1] <- TRUE  # recorded (integrated) vs dark cycles
      idx   <- which(integ)
      if (surv && length(idx) > 1L) {                   # survivorship: cut at first over-range dark-gap jump
        jmp  <- sqrt(rowSums((pos[idx[-1], , drop = FALSE] - pos[idx[-length(idx)], , drop = FALSE])^2))
        lost <- which(jmp > R_reacq)
        if (length(lost)) idx <- idx[seq_len(lost[1])]  # keep up to (and incl.) the pre-jump loc
      }
      if (length(idx) < 5L) next
      N     <- pmax(1L, stats::rpois(length(idx), photons))
      sig   <- L_tcp / (2 * sqrt(N))                   # per-localization σ_loc (nm)
      obs   <- pos[idx, , drop = FALSE] +
               matrix(rnorm(length(idx) * 3), ncol = 3) * sig
      k <- k + 1L
      rows[[k]] <- data.frame(
        tid = k, tim = ((idx - 1L) * cyc) / 1000,      # SECONDS (importer ×1000 → ms)
        x = obs[, 1], y = obs[, 2], z = obs[, 3],
        efo = N / (cyc / 1000), eco = N, stringsAsFactors = FALSE)
      gt[[k]] <- data.frame(tid = k, D_true_um2_s = D_um2_s, geometry = geom, mode = mode,
                            motion = motion, L_c_nm = if (identical(motion,"confined")) L_c else NA_real_,
                            v_um_s = if (identical(motion,"directed")) v_nmms else NA_real_,
                            n_locs = length(idx), sigma_loc_med = stats::median(sig),
                            stringsAsFactors = FALSE)
    }
    })
    if (k == 0L) return(NULL)
    df <- dplyr::bind_rows(rows[seq_len(k)])
    if (proj2d) df$z <- NULL
    list(df = df, gt = dplyr::bind_rows(gt[seq_len(k)]),
         D_um2_s = D_um2_s, cyc = cyc, p_dark = p_dark, proj2d = proj2d,
         mode = mode, d_mot = d_mot, geom = geom,
         motion = motion, L_c = L_c, v_nmms = v_nmms, surv = surv, R_reacq = R_reacq)
  })

  # CVE round-trip: recover D from consecutive single-cycle (Δt = cyc) steps, using the
  # TOTAL displacement with the MOTION dimensionality. D = (⟨|Δr|²⟩ + 2⟨Δr_i·Δr_{i+1}⟩)
  # / (2·d·Δt). SURFACE diffusion has d=2 DOF but spreads over 3 Cartesian axes (a
  # per-axis CVE would under-estimate D by ⅓); VOLUME diffusion has d=3. When projected
  # to 2D only n_cart=2 axes are seen, so d = min(motion-DOF, n_cart). n_cart Cartesian
  # axes carry the noise (cov = −n_cart·σ²).
  sim_trk_cve <- reactive({
    res <- sim_trk_result(); if (is.null(res)) return(NULL)
    df <- res$df; cyc <- res$cyc
    axes  <- c("x", "y", if (!is.null(df$z)) "z")
    ncart <- length(axes)
    d_mot <- min(res$d_mot %||% 2L, ncart)              # 2 surface / 3 volume, capped by observed axes
    r2_sum <- 0; r2_n <- 0L; cov_sum <- 0; cov_n <- 0L
    for (tb in split(df, df$tid)) {
      tb  <- tb[order(tb$tim), ]
      dt  <- c(NA_real_, diff(tb$tim) * 1000)           # ms
      one <- is.finite(dt) & abs(dt - cyc) < cyc * 0.25 # single-cycle step (ending at this row)
      # per-step displacement vector (rows = steps ending at row i)
      dmat <- sapply(axes, function(ax) c(NA_real_, diff(tb[[ax]])))
      if (is.null(dim(dmat))) dmat <- matrix(dmat, ncol = ncart)
      r2  <- rowSums(dmat^2)                             # |Δr|² total
      s1  <- one & is.finite(r2)
      r2_sum <- r2_sum + sum(r2[s1]); r2_n <- r2_n + sum(s1)
      nxt  <- c(s1[-1], FALSE); pair <- s1 & nxt         # consecutive single-cycle pair
      dot  <- rowSums(dmat * rbind(dmat[-1, , drop = FALSE], NA_real_))  # Δr_i · Δr_{i+1}
      cov_sum <- cov_sum + sum(dot[pair], na.rm = TRUE); cov_n <- cov_n + sum(pair)
    }
    if (r2_n < 1L) return(NULL)
    r2m  <- r2_sum / r2_n                                # ⟨|Δr|²⟩ (nm²)
    cov1 <- if (cov_n > 0L) cov_sum / cov_n else 0       # ⟨Δr_i·Δr_{i+1}⟩ = −n_cart·σ²
    list(D_um2_s  = ((r2m + 2 * cov1) / (2 * d_mot * cyc)) / 1000,   # CVE D (µm²/s)
         sigma_nm = sqrt(max(0, -cov1 / ncart)), n_pairs = cov_n)
  })

  output$sim_trk_plot <- renderPlotly({
    res <- sim_trk_result()
    shiny::validate(shiny::need(!is.null(res), "Set parameters and click 'Simulate tracks'."))
    df <- res$df; is3d <- !is.null(df$z)
    show <- head(unique(df$tid), 60L); d <- df[df$tid %in% show, ]
    if (is3d)
      plotly::plot_ly(d, x = ~x, y = ~y, z = ~z, split = ~factor(tid),
                      type = "scatter3d", mode = "lines",
                      line = list(width = 2), showlegend = FALSE) %>%
        plotly::layout(title = "Simulated MINFLUX tracks (first 60)",
          scene = list(xaxis = list(title = "x (nm)"), yaxis = list(title = "y (nm)"),
                       zaxis = list(title = "z (nm)")))
    else
      plotly::plot_ly(d, x = ~x, y = ~y, split = ~factor(tid),
                      type = "scatter", mode = "lines",
                      line = list(width = 1), showlegend = FALSE) %>%
        plotly::layout(title = "Simulated MINFLUX tracks — 2D projection (first 60)",
          xaxis = list(title = "x (nm)"),
          yaxis = list(title = "y (nm)", scaleanchor = "x"))
  })

  output$sim_trk_summary <- renderUI({
    res <- sim_trk_result(); if (is.null(res)) return(NULL)
    df <- res$df; cve <- sim_trk_cve()
    dt_all <- unlist(lapply(split(df$tim, df$tid), function(t) diff(sort(t)) * 1000), use.names = FALSE)
    frac_int <- if (length(dt_all)) mean(abs(dt_all - res$cyc) < res$cyc * 0.25) else NA_real_
    recov <- if (!is.null(cve)) {
      rr  <- cve$D_um2_s / res$D_um2_s
      col <- if (is.finite(rr) && abs(rr - 1) < 0.25) "#228833" else "#BB7700"
      htmltools::tags$span(style = sprintf("color:%s; font-weight:bold;", col),
        sprintf("%.4g µm²/s  (%.0f%% of truth)", cve$D_um2_s, 100 * rr))
    } else "—"
    motion_lbl <- switch(res$motion %||% "free",
      confined = sprintf("confined (OU tether, L_c = %.0f nm)", res$L_c),
      directed = sprintf("directed (v = %.3g µm/s)", res$v_nmms),
      twostate = "two-state (mobile↔immobile)",
      "free")
    geom_lbl <- sprintf("%s · %s · %s%s", res$geom %||% "?",
                        if (identical(res$mode, "volume")) "volume (intracellular)" else "surface",
                        motion_lbl, if (isTRUE(res$surv)) sprintf(" · survivorship R=%.0fnm", res$R_reacq) else "")
    htmltools::div(style = "font-size:0.9em; color:#333; margin-top:8px;",
      htmltools::tags$table(style = "border-collapse:collapse;",
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Model"),
          htmltools::tags$td(htmltools::tags$b(geom_lbl))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Tracks / localizations"),
          htmltools::tags$td(sprintf("%d / %s", length(unique(df$tid)), format(nrow(df), big.mark=",")))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Ground-truth D"),
          htmltools::tags$td(htmltools::tags$b(sprintf("%.4g µm²/s", res$D_um2_s)))),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Recovered D (CVE round-trip)"),
          htmltools::tags$td(recov)),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Recovered σ_loc"),
          htmltools::tags$td(if (!is.null(cve)) sprintf("%.2f nm", cve$sigma_nm) else "—")),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Integrated (single-cycle) fraction"),
          htmltools::tags$td(if (is.finite(frac_int)) sprintf("%.0f%%", 100 * frac_int) else "—")),
        htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Median Δt / median track length"),
          htmltools::tags$td(sprintf("%.3g ms / %.0f locs",
            if (length(dt_all)) stats::median(dt_all) else NA_real_, stats::median(as.numeric(table(df$tid))))))),
      htmltools::tags$p(style = "color:#888; font-size:0.9em; margin-top:6px;",
        if (identical(res$motion, "confined")) {
          tau_c   <- if (res$D_um2_s > 0) (res$L_c^2 / (res$D_um2_s * 1000)) else Inf   # ms
          med_dur <- stats::median(unlist(lapply(split(df$tim, df$tid), function(t) (max(t) - min(t)) * 1000)))
          verdict <- if (is.finite(tau_c) && tau_c < 0.5 * med_dur) "✓ RESOLVABLE" else "⚠ NOT resolvable — looks free"
          sprintf("Confined (OU, L_c=%.0f nm): CVE gives the LOCAL D at single-cycle lags (≈ truth); confinement is an MSD plateau ≈ 2·d·L_c² in the Tracking tab. Resolvability: τ_c = L_c²/D ≈ %.3g ms vs median track duration ≈ %.3g ms → %s. (If τ_c ≳ track duration, free diffusion over a short track mimics the same extent — the key artifact to test.)",
                  res$L_c, tau_c, med_dur, verdict)
        }
        else if (identical(res$motion, "directed"))
          "Directed: the CVE local D reflects the diffusive part; the drift shows as an upward-curving MSD / a sustained drift-velocity in the Tracking tab."
        else if (identical(res$motion, "twostate"))
          "Two-state: mixed mobile/immobile steps — the CVE gives a population-average D; load in the Tracking tab and fit the HMM to recover the two states."
        else
          "Free diffusion: the CVE round-trip should land near 100 % of truth — that is the validation loop.",
        " Download the tracks (CSV) and load them in the Tracking tab to test the full MSD / confinement pipeline against this known ground truth."),
      # Projection artifact: 2D imaging of a CURVED-surface protein under-measures the
      # azimuthal motion (its out-of-plane component is dropped), so the recovered D is a
      # fixed fraction of truth (~75% for a side-on cylinder) — NOT a bug, and the same
      # bias applies to real 2D data. The factor is itself a useful geometric calibration.
      if (identical(res$mode, "surface") && isTRUE(res$proj2d))
        htmltools::tags$p(style = "color:#b26a00; font-size:0.9em; margin-top:4px;",
          htmltools::HTML("&#9888; "),
          htmltools::tags$b("2D projection of a curved surface:"),
          " recovered D is expected BELOW truth (≈75% for a side-on cylinder) because the ",
          "out-of-plane part of the azimuthal motion is lost — this is a real geometric bias that ",
          "also affects your real 2D data. Use this ratio as the correction factor, or acquire in 3D. ",
          "It applies to free AND confined alike, so it is not the confinement signal."))
  })

  output$sim_trk_dl_csv <- downloadHandler(
    filename = function() paste0("sim_minflux_tracks_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      res <- sim_trk_result()
      shiny::validate(shiny::need(!is.null(res), "Simulate tracks first."))
      write.csv(res$df, file, row.names = FALSE)
    })
  output$sim_trk_dl_gt <- downloadHandler(
    filename = function() paste0("sim_ground_truth_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      res <- sim_trk_result()
      shiny::validate(shiny::need(!is.null(res), "Simulate tracks first."))
      write.csv(res$gt, file, row.names = FALSE)
    })

}
