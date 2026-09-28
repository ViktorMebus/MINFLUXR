## ─────────────────────────────────────────────────────────────────────────────
## est_hmm.R  -  continuous-time diffusive-state HMM for irregularly sampled tracks
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Estimate K diffusive states (D_k) and the SWITCHING RATES between them
##   (generator Q, per ms) directly by maximum likelihood, from tracks whose
##   sampling interval Δt varies from step to step (MINFLUX dark cycles).
##   Pure functions, no Shiny: the Tracking tab's hmm_bundle() calls them, and
##   tests/ can call them against simulated ground truth.
##
## WHY THIS REPLACED THE OLD HMM
##   The previous version fitted a per-STEP transition matrix by Baum-Welch and
##   then counted switches along the Viterbi path to get rates. Counting on a
##   decoded path misses dwells shorter than a few steps (they are absorbed into
##   the neighbouring state), so dwell times read long and rates read low, and the
##   per-step matrix silently assumed every step spans the same time. Here the
##   rates are PARAMETERS of the likelihood, and each step's transition matrix is
##   P(Δt) = expm(Q·Δt), which is exact for any sampling pattern.
##
## THE MODEL
##   Hidden state s_j ∈ 1..K for step j (the displacement from localization j to
##   j+1), a continuous-time Markov chain with generator Q:
##       P(s_{j+1} = b | s_j = a) = [expm(Q·Δt_j)]_{ab}
##   Emission: the squared step length r² over d observed axes. Per axis the
##   displacement is Gaussian with variance
##       v_j = 2·D_s·(Δt_j − 2·R_abs) + σ_j² + σ_{j+1}²
##   (Brownian motion + motion blur + the localization noise of both ends;
##   R_abs = R·t_c is the blur offset in ms, from blur_offsets() in est_blur.R), so
##       r² ~ Gamma(shape = d/2, scale = 2·v_j).
##   Initial state distribution: the stationary distribution of Q.
##
## STATED APPROXIMATIONS (measured by the recovery sweep, tests/validation/)
##   1. A step is assigned ONE state, so a switch inside a step is not modelled;
##      this biases rates when dwells approach Δt.
##   2. Consecutive steps share a localization, so their noise is correlated
##      (−σ² covariance); the HMM treats steps as conditionally independent. Small
##      when D·Δt ≫ σ², larger for slow states.
##   3. Tracks are cut into chunks of `chunk_len` steps, each starting from the
##      stationary distribution (a composite likelihood). Only transitions across
##      chunk boundaries are lost (1/chunk_len of them); the marginal is exact
##      for a stationary chain.
##
## DATA FLOW
##   hmm_ct_prepare(df)  -> chunked step matrices (r², Δt, noise, Δt index)
##   hmm_ct_fit(prep)    -> D, Q, rates, dwell times, CIs, log-likelihood, BIC
##   hmm_ct_viterbi(prep, fit) -> per-localization state labels (display only)
##   hmm_ct_bootstrap()  -> track-resampled CIs (the honest ones; see below)
##   sim_ct_state_tracks() -> ground-truth tracks for validation
##   hmm_path_rates()    -> the OLD decoded-path rate estimate, kept only so the
##                          validation can show its bias next to the MLE
##
## UNITS: positions nm, time ms, D nm²/ms (÷1000 → µm²/s), rates per ms (×1000 → s⁻¹).
## ─────────────────────────────────────────────────────────────────────────────


# ── .expm_small() ─────────────────────────────────────────────────────────────
# GOAL: matrix exponential of a small (K ≤ 4) generator times a time.
# Scaling and squaring with a 12-term Taylor series: scale A so its 1-norm is
# ≤ 0.5, where the truncation error is < 1e-16, then square back. No package
# dependency (expm is not installed with the app).
# Returns a row-stochastic matrix (tiny negative round-off clamped, rows renormalised).
.expm_small <- function(A) {
  n1 <- max(colSums(abs(A)))                              # 1-norm
  s  <- if (is.finite(n1) && n1 > 0.5) ceiling(log2(n1 / 0.5)) else 0L
  A  <- A / 2^s                                           # scaled so the series converges fast
  E  <- diag(nrow(A)); term <- E
  for (k in 1:12) { term <- term %*% A / k; E <- E + term }
  for (i in seq_len(s)) E <- E %*% E                      # undo the scaling
  E[E < 0] <- 0                                           # round-off only
  E / rowSums(E)                                          # keep rows exactly stochastic
}

# ── .ctmc_stationary() ────────────────────────────────────────────────────────
# GOAL: stationary distribution π of a generator Q (π Q = 0, Σπ = 1).
.ctmc_stationary <- function(Q) {
  K <- nrow(Q)
  A <- rbind(t(Q), rep(1, K))                             # K equations + normalisation
  p <- tryCatch(qr.solve(A, c(rep(0, K), 1)), error = function(e) rep(1 / K, K))
  p[!is.finite(p) | p < 0] <- 0
  if (sum(p) <= 0) rep(1 / K, K) else p / sum(p)
}

# ── .hmm_unpack() ─────────────────────────────────────────────────────────────
# GOAL: map the unconstrained parameter vector to (D, Q).
# D_1 = exp(θ_1), D_k = D_{k−1} + exp(θ_k): positive AND ordered slow → fast, so
# state labels cannot swap between optimiser iterations. Off-diagonal rates are
# exp(θ), filled row-wise; the diagonal makes rows sum to zero.
.hmm_unpack <- function(theta, K) {
  D <- cumsum(exp(theta[seq_len(K)]))
  Q <- matrix(0, K, K)
  # fill row-wise: order off-diagonal cells by (row, col)
  idx <- which(row(Q) != col(Q), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  Q[idx] <- exp(theta[K + seq_len(nrow(idx))])
  diag(Q) <- -rowSums(Q)
  list(D = D, Q = Q, idx = idx)
}

# ── .hmm_pack() ───────────────────────────────────────────────────────────────
# GOAL: inverse of .hmm_unpack for starting values (D must be increasing).
.hmm_pack <- function(D, Q) {
  K  <- length(D)
  D  <- cummax(pmax(D, 1e-6))
  dD <- c(D[1], pmax(diff(D), 1e-6 * max(D)))             # strictly positive increments
  idx <- which(row(Q) != col(Q), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  c(log(dD), log(pmax(Q[idx], 1e-12)))
}


# ── hmm_ct_prepare() ──────────────────────────────────────────────────────────
# GOAL: turn localizations into the chunked step matrices the likelihood reads.
#
# df        : data.frame(tid, x, y[, z], tim) with tim in ms; optional `sig`
#             column = per-localization σ (nm). Rows need not be sorted.
# ndim      : 2 (x, y) or 3 (x, y, z)
# sigma     : fallback σ_loc (nm) when there is no usable `sig` column
# chunk_len : steps per chunk (see approximation 3 in the file header)
# min_pts   : minimum localizations per track
# max_dt    : steps longer than this (ms) are treated as breaks, like an invalid Δt.
#             Set just above one cycle to fit single-cycle steps only: on real
#             files multi-cycle steps are compressed (⟨Δx²⟩ barely grows with the
#             number of cycles), which the Δt-linear emission reads as a slow state.
#
# A non-positive or non-finite Δt breaks a track into two segments (no step is
# invented across it). Chunks are sorted longest first, so at column j the
# active chunks are always rows 1..n_act[j]; the loops then subset a prefix.
#
# RETURNS list(R2, DT, NOISE, UIDX [n_chunk × chunk_len], len (steps per chunk),
#   n_act (active rows per column), dt_u (unique Δt values), ndim,
#   steps = data.frame(tid, loc, chunk, col) mapping every step back to its
#   track and localization index (loc = position of the step's START within the
#   sorted track), tids, n_steps, total_time_ms)
hmm_ct_prepare <- function(df, ndim = 2L, sigma = 0, chunk_len = 256L, min_pts = 10L,
                           max_dt = Inf) {
  axes <- if (ndim == 3L) c("x", "y", "z") else c("x", "y")
  need <- c("tid", "tim", axes)
  if (is.null(df) || !nrow(df) || !all(need %in% names(df))) return(NULL)
  df <- df[stats::complete.cases(df[, need]), , drop = FALSE]
  df$tid <- as.character(df$tid)
  df <- df[order(df$tid, df$tim), , drop = FALSE]
  cnt <- table(df$tid)
  df  <- df[df$tid %in% names(cnt)[cnt >= min_pts], , drop = FALSE]
  n <- nrow(df); if (n < 3L) return(NULL)

  sig <- if ("sig" %in% names(df)) suppressWarnings(as.numeric(df$sig)) else rep(NA_real_, n)
  sig_fallback <- if (is.finite(sigma) && sigma >= 0) sigma else 0
  sig[!is.finite(sig) | sig < 0] <- sig_fallback

  same <- df$tid[-1] == df$tid[-n]                        # step j joins row j and j+1
  dt   <- diff(df$tim)
  r2   <- Reduce(`+`, lapply(axes, function(a) diff(df[[a]])^2))
  noise <- sig[-n]^2 + sig[-1]^2                          # per-axis noise variance of the step
  ok   <- same & is.finite(dt) & dt > 0 & dt <= max_dt & is.finite(r2)
  # position of each row inside its track (1-based), for mapping states back
  loc_in_track <- stats::ave(seq_len(n), df$tid, FUN = seq_along)

  # Segments: a new segment starts at every invalid step or new track.
  j   <- which(ok)                                        # valid step indices (row of step start)
  if (length(j) < 2L) return(NULL)
  brk <- c(TRUE, diff(j) != 1L | df$tid[j[-1]] != df$tid[j[-length(j)]])
  seg <- cumsum(brk)
  pos <- stats::ave(seq_along(j), seg, FUN = seq_along)   # step position inside its segment
  chunk_key <- paste(seg, (pos - 1L) %/% chunk_len)
  chunk_id  <- match(chunk_key, unique(chunk_key))
  col       <- (pos - 1L) %% chunk_len + 1L
  len       <- tabulate(chunk_id)
  # longest chunks first -> active rows are always a prefix
  ord_chunk <- order(-len, seq_along(len))
  new_id    <- match(chunk_id, ord_chunk)
  len       <- len[ord_chunk]
  nc <- length(len); L <- max(len)

  # Unique Δt values (MINFLUX Δt is quantised at k·t_c + SW, so there are few).
  dt_key <- round(dt[j], 6)
  dt_u   <- sort(unique(dt_key))
  if (length(dt_u) > 400L) {                              # jittered clock: bin coarsely
    dt_key <- signif(dt[j], 3); dt_u <- sort(unique(dt_key))
  }
  mk <- function(v, fill) { m <- matrix(fill, nc, L); m[cbind(new_id, col)] <- v; m }
  R2    <- mk(r2[j], NA_real_)
  DT    <- mk(dt[j], NA_real_)
  NOISE <- mk(noise[j], NA_real_)
  UIDX  <- mk(match(dt_key, dt_u), NA_integer_)
  n_act <- vapply(seq_len(L), function(cc) sum(len >= cc), integer(1))

  # r2 transforms the 3D emission density needs; constant across the fit, so they
  # are computed once here rather than inside the optimiser's inner loop.
  R2_SQ  <- if (ndim == 3L) sqrt(R2) else NULL
  R2_LOG <- if (ndim == 3L) log(pmax(R2, 1e-300)) else NULL

  list(R2 = R2, R2_SQ = R2_SQ, R2_LOG = R2_LOG,
       DT = DT, NOISE = NOISE, UIDX = UIDX, len = len, n_act = n_act,
       dt_u = dt_u, ndim = as.integer(ndim), chunk_len = as.integer(chunk_len),
       steps = data.frame(tid = df$tid[j], loc = loc_in_track[j], row = j,
                          chunk = new_id, col = col, stringsAsFactors = FALSE),
       n_steps = length(j), total_time_ms = sum(dt[j]),
       tids = unique(df$tid[j]), df = df)
}


# ── .hmm_emissions() ──────────────────────────────────────────────────────────
# GOAL: emission densities for every step and state, n_chunk × L × K.
# Padding cells get density 1 (they are never read by the prefix loops, but a
# finite value keeps array arithmetic clean).
#
# The gamma density is written out rather than called through stats::dgamma().
# This runs once per state on every optimiser iteration, over every step in the
# data, so it dominates the fit time. The shape is fixed by the dimensionality
# (1 in 2D, 3/2 in 3D), and at those two values the density is elementary:
#   shape 1   : e = exp(-r2/c) / c
#   shape 3/2 : e = sqrt(r2) * exp(-r2/c) / (Gamma(3/2) * c^(3/2))
# with c = 2v. The sqrt(r2) and log(r2) that the 3D form needs do not change
# between iterations, so hmm_ct_prepare() computes them once.
.hmm_emissions <- function(prep, D, R_abs, log = FALSE) {
  K <- length(D); shape <- prep$ndim / 2
  eff <- pmax(prep$DT - 2 * R_abs, 1e-9)                  # blur-shortened diffusion time
  E <- array(if (log) 0 else 1, dim = c(dim(prep$R2), K))
  pad <- is.na(prep$R2)
  r2  <- prep$R2
  # Cached r2 transforms for the 3D form (recomputed if an older prep lacks them).
  r2_sq  <- if (shape != 1) prep$R2_SQ  %||% sqrt(r2)         else NULL
  r2_log <- if (shape != 1) prep$R2_LOG %||% log(pmax(r2, 1e-300)) else NULL
  for (s in seq_len(K)) {
    v <- 2 * D[s] * eff + prep$NOISE                      # per-axis variance of the step
    cc <- 2 * pmax(v, 1e-9)                               # gamma scale
    e <- if (shape == 1) {
      if (log) -r2 / cc - log(cc) else exp(-r2 / cc) / cc
    } else if (shape == 1.5) {
      if (log) 0.5 * r2_log - r2 / cc - lgamma(1.5) - 1.5 * log(cc)
      else     r2_sq * exp(-r2 / cc) / (gamma(1.5) * cc^1.5)
    } else stats::dgamma(r2, shape = shape, scale = cc, log = log)
    if (!log) e <- pmax(e, 1e-300)                        # avoid exact zeros in products
    e[pad] <- if (log) 0 else 1
    E[, , s] <- e
  }
  E
}

# ── .hmm_transitions() ────────────────────────────────────────────────────────
# GOAL: P(Δt) for every unique Δt, as a U × K × K array.
.hmm_transitions <- function(Q, dt_u) {
  K <- nrow(Q); U <- length(dt_u)
  P <- array(0, dim = c(U, K, K))
  for (u in seq_len(U)) P[u, , ] <- .expm_small(Q * dt_u[u])
  P
}


# ── hmm_ct_loglik() ───────────────────────────────────────────────────────────
# GOAL: composite log-likelihood of all chunks (scaled forward algorithm),
# vectorised over chunks. theta as in .hmm_unpack. Returns a finite number
# (−1e300 on numerical failure, so the optimiser steps away).
hmm_ct_loglik <- function(theta, prep, K, R_abs = 0) {
  par <- .hmm_unpack(theta, K)
  if (any(!is.finite(par$D)) || any(!is.finite(par$Q))) return(-1e300)
  E  <- .hmm_emissions(prep, par$D, R_abs)
  P  <- .hmm_transitions(par$Q, prep$dt_u)
  pi0 <- .ctmc_stationary(par$Q)
  nc <- nrow(prep$R2); L <- ncol(prep$R2)

  alpha <- sweep(matrix(E[, 1, ], nc, K), 2, pi0, `*`)   # first step: stationary start
  cs <- rowSums(alpha); ll <- sum(log(cs)); alpha <- alpha / cs
  if (L >= 2L) for (j in 2:L) {
    na <- prep$n_act[j]; if (na == 0L) break
    rows <- seq_len(na)
    u    <- prep$UIDX[rows, j - 1L]                       # Δt of the step that just ended
    a    <- alpha[rows, , drop = FALSE]
    new  <- matrix(0, na, K)
    for (b in seq_len(K)) {
      acc <- 0
      for (s in seq_len(K)) acc <- acc + a[, s] * P[u, s, b]
      new[, b] <- acc * E[rows, j, b]
    }
    cs <- rowSums(new)
    ll <- ll + sum(log(cs))
    alpha[rows, ] <- new / cs
  }
  if (!is.finite(ll)) -1e300 else ll
}


# ── hmm_ct_fit() ──────────────────────────────────────────────────────────────
# GOAL: maximum-likelihood D_k and Q, with Hessian-based CIs.
#
# prep     : from hmm_ct_prepare()
# K        : number of states (2-4)
# R_abs    : blur offset in ms (blur_offsets()$R_abs, est_blur.R); 0 = no blur
# t_c      : cycle time (ms), used only to spread the starting rates
# n_starts : starting points (rate scales ~5, 30 and 300 cycles dwell); the best
#            log-likelihood wins. The likelihood is multimodal in the rates.
# theta0   : optional start (e.g. the point estimate, for bootstrap refits)
#
# CIs: the Hessian is the observed information of the composite likelihood and
# treats every step as coming from ONE homogeneous population. Real molecules
# differ, so these intervals are too narrow on real data; hmm_ct_bootstrap()
# (resampling whole tracks) gives the honest ones.
#
# RETURNS list(K, D (nm²/ms), Q (per ms), rates_per_s (K×K, off-diagonal),
#   dwell_ms, pi (stationary occupancy), loglik, n_par, BIC, converged,
#   theta, vcov, ci = list(D, dwell_ms, rates_per_s) with lo/hi matrices,
#   expected_switches (K×K, total observed time × occupancy × rate))
hmm_ct_fit <- function(prep, K = 2L, R_abs = 0, t_c = NA_real_, n_starts = 3L,
                       theta0 = NULL, hessian = TRUE, maxit = 300L, seed = 1L) {
  stopifnot(!is.null(prep), K >= 2L)
  obs <- prep$R2[!is.na(prep$R2)]
  eff <- pmax(prep$DT[!is.na(prep$R2)] - 2 * R_abs, 1e-9)
  nz  <- prep$NOISE[!is.na(prep$R2)]
  # per-step D estimate for starting values: (r²/d − noise) / (2·eff)
  Dhat <- (obs / prep$ndim - nz) / (2 * eff)
  pos  <- Dhat[is.finite(Dhat) & Dhat > 0]
  base <- if (length(pos)) stats::median(pos) else 100
  D0   <- if (length(pos) > 20)
            stats::quantile(pos, seq(0.15, 0.85, length.out = K), names = FALSE)
          else base * 2^(seq_len(K) - 1)
  D0   <- sort(pmax(D0, base * 1e-3))
  dt_ref <- if (is.finite(t_c) && t_c > 0) t_c else stats::median(prep$DT, na.rm = TRUE)

  starts <- if (!is.null(theta0)) list(theta0) else
    lapply(c(5, 30, 300)[seq_len(max(1L, min(3L, n_starts)))], function(cyc) {
      Q0 <- matrix(1 / (cyc * dt_ref * (K - 1)), K, K); diag(Q0) <- 0
      diag(Q0) <- -rowSums(Q0)
      .hmm_pack(D0, Q0)
    })
  nll <- function(th) -hmm_ct_loglik(th, prep, K, R_abs)
  best <- NULL
  for (th0 in starts) {
    o <- tryCatch(stats::optim(th0, nll, method = "BFGS",
                               control = list(maxit = maxit, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(NULL)
  par <- .hmm_unpack(best$par, K)
  ll  <- -best$value
  n_par <- K + K * (K - 1)

  V <- NULL
  if (isTRUE(hessian)) {
    H <- tryCatch(stats::optimHess(best$par, nll), error = function(e) NULL)
    V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
    if (!is.null(V) && any(!is.finite(diag(V)) | diag(V) < 0)) V <- NULL
  }
  out <- .hmm_summarise(best$par, K, prep, V, seed)
  c(out, list(loglik = ll, n_par = n_par,
              BIC = -2 * ll + n_par * log(prep$n_steps),
              converged = identical(best$convergence, 0L), R_abs = R_abs))
}

# ── .hmm_summarise() ──────────────────────────────────────────────────────────
# GOAL: derived quantities + CIs from θ and its covariance (parametric draws, so
# the nonlinear transforms to D, rates and dwell times need no delta method).
.hmm_summarise <- function(theta, K, prep, V = NULL, seed = 1L, n_draw = 2000L) {
  par  <- .hmm_unpack(theta, K)
  Q    <- par$Q
  pi_s <- .ctmc_stationary(Q)
  rates_s <- Q * 1000; diag(rates_s) <- NA
  dwell <- -1 / diag(Q)
  T_tot <- prep$total_time_ms
  exp_sw <- T_tot * pi_s * Q; diag(exp_sw) <- NA          # expected observed switches
  ci <- NULL
  if (!is.null(V)) {
    set.seed(seed)
    Lc <- tryCatch(chol((V + t(V)) / 2), error = function(e) NULL)
    if (!is.null(Lc)) {
      Z  <- matrix(stats::rnorm(n_draw * length(theta)), n_draw) %*% Lc
      TH <- sweep(Z, 2, theta, `+`)
      draws <- lapply(seq_len(n_draw), function(i) .hmm_unpack(TH[i, ], K))
      qfun <- function(M) apply(M, 2, stats::quantile, c(0.025, 0.975), na.rm = TRUE)
      Dm <- qfun(t(vapply(draws, function(p) p$D, numeric(K))))
      Wm <- qfun(t(vapply(draws, function(p) -1 / diag(p$Q), numeric(K))))
      Rm <- qfun(t(vapply(draws, function(p) as.numeric(p$Q * 1000), numeric(K * K))))
      ci <- list(D = Dm, dwell_ms = Wm,
                 rates_lo = matrix(Rm[1, ], K), rates_hi = matrix(Rm[2, ], K))
    }
  }
  list(K = K, D = par$D, Q = Q, rates_per_s = rates_s, dwell_ms = dwell, pi = pi_s,
       theta = theta, vcov = V, ci = ci, expected_switches = exp_sw,
       n_steps = prep$n_steps, n_tracks = length(prep$tids))
}


# ── hmm_ct_viterbi() ──────────────────────────────────────────────────────────
# GOAL: most likely state per step under the FITTED continuous-time model, for
# colouring tracks and for the per-state MSD. Rates are never estimated from it.
# Vectorised over chunks like the likelihood. Returns an integer vector of
# states aligned with prep$steps.
hmm_ct_viterbi <- function(prep, fit) {
  K  <- fit$K
  lE <- .hmm_emissions(prep, fit$D, fit$R_abs %||% 0, log = TRUE)
  lP <- log(pmax(.hmm_transitions(fit$Q, prep$dt_u), 1e-300))
  nc <- nrow(prep$R2); L <- ncol(prep$R2)
  delta <- sweep(matrix(lE[, 1, ], nc, K), 2, log(pmax(fit$pi, 1e-300)), `+`)
  psi   <- array(0L, dim = c(nc, L, K))
  if (L >= 2L) for (j in 2:L) {
    na <- prep$n_act[j]; if (na == 0L) break
    rows <- seq_len(na); u <- prep$UIDX[rows, j - 1L]
    d  <- delta[rows, , drop = FALSE]
    nd <- matrix(-Inf, na, K); arg <- matrix(1L, na, K)
    for (b in seq_len(K)) {
      best <- rep(-Inf, na); ba <- rep(1L, na)
      for (s in seq_len(K)) {
        cand <- d[, s] + lP[u, s, b]
        better <- cand > best
        best[better] <- cand[better]; ba[better] <- s
      }
      nd[, b] <- best + lE[rows, j, b]; arg[, b] <- ba
    }
    delta[rows, ] <- nd; psi[rows, j, ] <- arg
  }
  # backtrack: each chunk ends at its own length
  S <- matrix(NA_integer_, nc, L)
  S[cbind(seq_len(nc), prep$len)] <- max.col(delta, ties.method = "first")
  if (L >= 2L) for (j in L:2) {
    rows <- which(prep$len >= j)
    if (!length(rows)) next
    S[cbind(rows, j - 1L)] <- psi[cbind(rows, j, S[cbind(rows, j)])]
  }
  S[cbind(prep$steps$chunk, prep$steps$col)]
}

# ── hmm_ct_states_by_loc() ────────────────────────────────────────────────────
# GOAL: per-localization state labels (the step leaving a localization carries
# its state; a track's last localization inherits its last step). Returns the
# prepared data.frame with a `state` column; localizations with no valid step
# get NA.
hmm_ct_states_by_loc <- function(prep, step_states) {
  df <- prep$df
  st <- rep(NA_integer_, nrow(df))
  st[prep$steps$row] <- step_states
  # last localization of each run of valid steps: copy the previous state
  nxt <- prep$steps$row + 1L
  fill <- nxt[is.na(st[nxt]) & df$tid[nxt] == df$tid[prep$steps$row]]
  st[fill] <- st[fill - 1L]
  df$state <- st
  df
}


# ── hmm_ct_bootstrap() ────────────────────────────────────────────────────────
# GOAL: CIs that include between-molecule heterogeneity, by resampling whole
# TRACKS with replacement and refitting from the point estimate.
# Returns list(D, dwell_ms, rates_lo, rates_hi, n_ok) in the same layout as
# fit$ci, or NULL if too few refits succeed.
hmm_ct_bootstrap <- function(df, fit, B = 50L, ndim = 2L, sigma = 0, t_c = NA_real_,
                             chunk_len = 256L, min_pts = 10L, seed = 1L, progress = NULL) {
  set.seed(seed)
  tids <- unique(as.character(df$tid))
  sp   <- split(df, as.character(df$tid))
  K    <- fit$K
  res  <- vector("list", B)
  for (b in seq_len(B)) {
    if (is.function(progress)) progress(b / B)
    pick <- sample(tids, length(tids), replace = TRUE)
    bd   <- do.call(rbind, lapply(seq_along(pick), function(i) {
      d <- sp[[pick[i]]]; d$tid <- paste0(pick[i], "#", i); d }))  # duplicates stay distinct tracks
    pr <- hmm_ct_prepare(bd, ndim = ndim, sigma = sigma, chunk_len = chunk_len, min_pts = min_pts)
    if (is.null(pr)) next
    f  <- tryCatch(hmm_ct_fit(pr, K = K, R_abs = fit$R_abs, t_c = t_c,
                              theta0 = fit$theta, hessian = FALSE), error = function(e) NULL)
    if (!is.null(f)) res[[b]] <- f
  }
  res <- Filter(Negate(is.null), res)
  if (length(res) < 10L) return(NULL)
  q <- function(M) apply(M, 2, stats::quantile, c(0.025, 0.975), na.rm = TRUE)
  Rm <- q(t(vapply(res, function(f) as.numeric(f$Q * 1000), numeric(K * K))))
  list(D = q(t(vapply(res, function(f) f$D, numeric(K)))),
       dwell_ms = q(t(vapply(res, function(f) f$dwell_ms, numeric(K)))),
       rates_lo = matrix(Rm[1, ], K), rates_hi = matrix(Rm[2, ], K),
       n_ok = length(res))
}


# ── hmm_path_rates() ──────────────────────────────────────────────────────────
# GOAL: the OLD estimator, rates from counting switches along a decoded path:
# q[a,b] = N[a→b] / T[a], steps spanning more than `gap_tol_cyc` cycles skipped.
# Kept ONLY so the validation can report its bias next to the MLE; the app no
# longer shows it.
hmm_path_rates <- function(prep, step_states, K, t_c, gap_tol_cyc = 4L) {
  st <- prep$steps; st$state <- step_states
  st$dt <- prep$DT[cbind(st$chunk, st$col)]
  st <- st[order(st$row), ]
  Tst <- numeric(K); N <- matrix(0, K, K)
  nxt_same <- c(st$row[-1] == st$row[-nrow(st)] + 1L, FALSE)
  for (i in seq_len(nrow(st))) {
    if (is.finite(t_c) && round(st$dt[i] / t_c) > gap_tol_cyc) next
    a <- st$state[i]; Tst[a] <- Tst[a] + st$dt[i]
    if (nxt_same[i]) { b <- st$state[i + 1L]; if (b != a) N[a, b] <- N[a, b] + 1 }
  }
  Q <- sweep(N, 1, pmax(Tst, 1e-12), `/`); diag(Q) <- -rowSums(Q)
  list(Q = Q, dwell_ms = -1 / diag(Q))
}


# ── sim_ct_state_tracks() ─────────────────────────────────────────────────────
# GOAL: ground-truth tracks for validating the HMM, with the physics the model
# SIMPLIFIES kept in: states switch in continuous time (also inside a step),
# every localization is the average of the path over its collection window
# (motion blur), steps span 1..k cycles (dark cycles) plus the SW overhead, and
# each localization has its own σ.
#
# n_tracks, n_locs : tracks and localizations per track (n_locs may be a vector)
# D                : per-state D (nm²/ms), length K
# Q                : generator (per ms), K × K
# t_c, sw          : cycle time and per-step overhead (ms)
# gap_p            : probabilities of a step spanning 1, 2, 3, ... cycles
# sigma            : σ per localization (nm): a scalar or a vector to resample
# R_blur           : blur coefficient; collection window τ = 6·R·t_c
# ndim             : 2 or 3 recorded axes
# n_sub            : points averaged per window (midpoint rule)
#
# RETURNS data.frame(tid, tim, x, y[, z], sig, state) - `state` is the true state
# at the start of each localization's collection window.
sim_ct_state_tracks <- function(n_tracks, n_locs, D, Q, t_c = 0.192, sw = 0.011,
                                gap_p = c(0.72, 0.20, 0.05, 0.03), sigma = 4,
                                R_blur = 0.106, ndim = 2L, n_sub = 8L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- length(D); pi0 <- .ctmc_stationary(Q)
  tau <- min(t_c, 6 * R_blur * t_c)
  n_locs <- rep_len(n_locs, n_tracks)
  out <- vector("list", n_tracks)
  for (tr in seq_len(n_tracks)) {
    N  <- n_locs[tr]
    k  <- sample.int(length(gap_p), N - 1L, TRUE, prob = gap_p)
    t_loc <- c(0, cumsum(k * t_c + sw))                   # localization start times
    T_end <- t_loc[N] + tau
    # continuous-time state path (Gillespie)
    s <- sample.int(K, 1L, prob = pi0); tj <- 0; states <- s; times <- 0
    while (tj < T_end) {
      rate <- -Q[s, s]
      if (rate <= 0) break
      tj <- tj - log(stats::runif(1)) / rate
      p  <- Q[s, ]; p[s] <- 0
      s  <- sample.int(K, 1L, prob = p)
      times <- c(times, tj); states <- c(states, s)
    }
    # G(t) = ∫ D(s) ds, piecewise linear between jump times
    seg_end <- c(times[-1], max(T_end, times[length(times)]) + 1)
    G_knots <- c(0, cumsum(D[states] * (seg_end - times)))
    t_knots <- c(times, seg_end[length(seg_end)])
    Gfun <- function(t) stats::approx(t_knots, G_knots, xout = t, rule = 2)$y
    # sample times: n_sub midpoints inside every collection window
    offs <- (seq_len(n_sub) - 0.5) / n_sub * tau
    ts   <- as.vector(outer(t_loc, offs, `+`))
    ord  <- order(ts); ts_s <- ts[ord]
    dG   <- diff(c(0, Gfun(ts_s)))
    pos  <- matrix(0, N, ndim)
    for (a in seq_len(ndim)) {
      path <- numeric(length(ts))
      path[ord] <- cumsum(stats::rnorm(length(ts), 0, sqrt(2 * pmax(dG, 0))))
      pos[, a] <- rowMeans(matrix(path, N, n_sub))        # window average = blur
    }
    sig_i <- if (length(sigma) > 1L) sample(sigma, N, replace = TRUE) else rep(sigma, N)
    pos <- pos + matrix(stats::rnorm(N * ndim), N, ndim) * sig_i
    st_true <- states[findInterval(t_loc, times)]
    d <- data.frame(tid = paste0("t", tr), tim = t_loc, x = pos[, 1], y = pos[, 2],
                    sig = sig_i, state = st_true, stringsAsFactors = FALSE)
    if (ndim == 3L) d$z <- pos[, 3]
    out[[tr]] <- d
  }
  do.call(rbind, out)
}
