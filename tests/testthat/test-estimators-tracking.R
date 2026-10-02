# Recovery tests for the tracking estimators extracted from the server
# (est_cve.R, est_msd.R, est_jdd.R, est_dmap.R, est_cpd.R), against ground truth
# from sim_ct_state_tracks() (est_hmm.R) at tracking-like acquisition settings:
# t_c 0.192 ms, SW 0.011 ms, dark-cycle mix, per-localization photon σ, blur R from the
# timing model of the standard 2D sequence (MSD_R_FALLBACK, est_blur.R). These tests check the
# algebra at a known R; test-blur.R checks that R against a simulated MINFLUX estimator.
TCP_CRLB_PREFACTOR <- if (exists("TCP_CRLB_PREFACTOR")) TCP_CRLB_PREFACTOR else 3.6
for (f in c("R/est_blur.R", "R/est_hmm.R", "R/est_cve.R", "R/est_msd.R", "R/est_jdd.R", "R/est_dmap.R", "R/est_cpd.R"))
  source(file.path(.repo_root, f))

.t_c <- 0.192; .sw <- 0.011; .R <- MSD_R_FALLBACK
.free_tracks <- function(D, n_tracks = 200, n_locs = 300, seed = 1, gap_p = c(0.72, 0.2, 0.05, 0.03)) {
  set.seed(seed)
  eco <- pmax(20, rlnorm(1e5, log(36), 0.6))
  sig <- 100.8 / (4 * sqrt(eco))
  Q0  <- matrix(c(-1e-9, 1e-9, 1e-9, -1e-9), 2)          # one effective state
  df  <- sim_ct_state_tracks(n_tracks, n_locs, D = c(D, D * 1.000001), Q = Q0, t_c = .t_c, sw = .sw,
                             gap_p = gap_p, sigma = sig, R_blur = .R, seed = seed)
  attr(df, "sig_rms") <- sqrt(mean(sig^2))
  df
}

# CLAIM (methods doc §I-B.1): the CVE D is blur-independent and unbiased.
test_that("CVE recovers D on blurred, noisy, irregularly sampled tracks", {
  skip_if_not_slow()
  for (D in c(50, 215, 1000)) {
    df <- .free_tracks(D, seed = 21)
    ce <- cve_estimate(df, .t_c)
    expect_lt(abs(ce$D_cve / D - 1), 0.05)
  }
})

# CLAIM (§I-B MSD): OLSF with the blur term recovers σ; R = 0 under-reads it.
test_that("OLSF blur term recovers sigma from the MSD intercept", {
  set.seed(1)
  D <- 215; sig <- 5; dt <- 0.192; nsub <- 24; N <- 1e5
  # Each localization averages the path over a window of 6*R*dt at the start of the
  # cycle, which is exactly blur R (Berglund 2010), so the fit is tested at a known R.
  n_win <- round(6 * .R * nsub)
  sim1 <- function() { h <- dt / nsub
    path <- cumsum(rnorm(N * nsub, 0, sqrt(2 * D * h)))
    rowMeans(matrix(path, ncol = nsub, byrow = TRUE)[, seq_len(n_win), drop = FALSE]) + rnorm(N, 0, sig) }
  x <- sim1(); y <- sim1(); lags <- 1:50
  msd <- sapply(lags, function(k) mean((x[-(1:k)] - x[1:(N - k)])^2 + (y[-(1:k)] - y[1:(N - k)])^2))
  tau <- lags * dt
  f0 <- msd_fit_olsf(tau, msd, rep(1L, 50), k = 50, auto_pmin = TRUE, R = 0, Npts = 1e4, dt_cycle = dt)
  fb <- msd_fit_olsf(tau, msd, rep(1L, 50), k = 50, auto_pmin = TRUE, R = .R, Npts = 1e4, dt_cycle = dt)
  fi <- msd_fit_olsf_iterative(tau, msd, R = .R, dt_cycle = dt, N = 50)
  expect_lt(abs(fb$D / D - 1), 0.03)
  expect_lt(abs(fb$sigma_loc / sig - 1), 0.05)
  expect_lt(abs(fi$sigma_loc / sig - 1), 0.05)
  expect_lt(f0$sigma_loc / sig, 0.9)                        # R = 0 reads ~4.0 nm at the physical R
})

# CLAIM (§I-B.4, fixed 2026-09-17): JDD D is noise- and blur-corrected in both modes.
test_that("JDD recovers D with blur, single-cycle and pooled", {
  skip_if_not_slow()
  for (D in c(50, 1000)) {
    df <- .free_tracks(D, n_tracks = 300, seed = 3)
    for (io in c(TRUE, FALSE)) {
      cj <- jdd_collect_jumps(df, min_pts = 10, integ_only = io, t_c = .t_c, norm_dt = TRUE)
      jf <- jdd_fit(cj$jumps, cj$dt_all, n_comp = 1, sigma_loc = attr(df, "sig_rms"), R_abs = .R * .t_c)
      expect_lt(abs(jf$D_vals / D - 1), 0.05)
    }
  }
})

# CLAIM (§I-B.9, fixed 2026-09-17): the D-map is unbiased on a uniform blurred field.
test_that("D-map recovers a uniform D with blur", {
  skip_if_not_slow()
  df <- .free_tracks(215, n_tracks = 300, seed = 6)
  m  <- dmap_compute(df, 5L, 50L, .t_c, TRUE, attr(df, "sig_rms"), R_abs = .R * .t_c)
  w  <- m$N_mat > 0
  expect_lt(abs(weighted.mean(m$D_mat[w], m$N_mat[w]) / 215 - 1), 0.05)
  m0 <- dmap_compute(df, 5L, 50L, .t_c, TRUE, attr(df, "sig_rms"), R_abs = 0)
  expect_lt(weighted.mean(m0$D_mat[w], m0$N_mat[w]) / 215, 0.88)  # unblurred form ~ −20% at the physical R
})

# CLAIM (§I-B.3): scaled PELT cuts few single-state tracks, and segment D is unbiased.
test_that("change-point: few spurious segments on one free population, D unbiased", {
  skip_if_not_slow()
  df   <- .free_tracks(215, n_tracks = 150, n_locs = 300, seed = 7, gap_p = 1)
  segs <- cpd_segment_tracks(df, min_seg = 8L, noise_sig2 = attr(df, "sig_rms")^2,
                             R_abs = .R * .t_c)
  frac_cut <- mean(tapply(segs$n_states, segs$tid, max) > 1)
  expect_lt(frac_cut, 0.15)
  expect_lt(abs(stats::weighted.mean(segs$D_nm2_ms, segs$n_steps) / 215 - 1), 0.05)
})

test_that("photon CRLB reduces to L/(a sqrt N) for an immobile emitter with nu off", {
  expect_equal(.crlb_doughnut(100.8, 0, c(36, 100), pref = 3.6), 100.8 / (3.6 * sqrt(c(36, 100))))
  # motion adds the pattern's off-centre term, <r^2> = dim * D * t_loc
  g <- TCP_PATTERN_GEOM$hexagon
  expect_equal(.crlb_doughnut(100.8, 0, 36, D = 200, t_loc = 0.2, pref = 3.6, pattern = "hexagon"),
               100.8 / (3.6 * 6) * sqrt(1 + g$c_off * g$dim * 200 * 0.2 / 100.8^2))
  # the triangle loses precision off-centre faster than the hexagon
  expect_gt(.crlb_doughnut(100.8, 0, 36, D = 200, t_loc = 0.2, pref = 3.6, pattern = "triangle"),
            .crlb_doughnut(100.8, 0, 36, D = 200, t_loc = 0.2, pref = 3.6, pattern = "hexagon"))
  # unknown pattern names fall back to the default geometry
  expect_equal(.crlb_doughnut(100.8, 0, 36, D = 200, t_loc = 0.2, pref = 3.6, pattern = "zline"),
               .crlb_doughnut(100.8, 0, 36, D = 200, t_loc = 0.2, pref = 3.6, pattern = TCP_PATTERN_DEFAULT))
})

# CLAIM (packages.R, TCP_PATTERN_GEOM): the ideal prefactor and the off-centre
# coefficient of each pattern follow from the Fisher matrix of the multinomial photon
# split over an ideal doughnut zero. Recomputed here independently, so the table
# cannot drift from the physics it states.
test_that("TCP pattern constants match the Fisher-matrix CRLB", {
  verts <- function(pat, L) {
    r <- L / 2
    switch(pat,
      triangle   = { a <- 2 * pi * (0:2) / 3; cbind(r * cos(a), r * sin(a)) },
      square     = { a <- 2 * pi * (0:3) / 4; cbind(r * cos(a), r * sin(a)) },
      hexagon    = { a <- 2 * pi * (0:5) / 6; cbind(r * cos(a), r * sin(a)) },
      octahedron = r * rbind(diag(3), -diag(3)))
  }
  fisher <- function(b, pos, h = 1e-6) {
    d <- ncol(b)
    p <- function(q) { lam <- rowSums((matrix(q, nrow(b), d, byrow = TRUE) - b)^2) + 1e-14; lam / sum(lam) }
    p0 <- p(pos)
    G <- sapply(seq_len(d), function(j) { e <- numeric(d); e[j] <- h; (p(pos + e) - p(pos - e)) / (2 * h) })
    t(G) %*% (G / p0)                                   # per photon
  }
  dirs <- function(d, n = 200) {
    if (d == 2) { a <- 2 * pi * (0:(n - 1)) / n; return(cbind(cos(a), sin(a))) }
    i <- 0:(n - 1) + 0.5; ph <- acos(1 - 2 * i / n); th <- pi * (1 + sqrt(5)) * i
    cbind(cos(th) * sin(ph), sin(th) * sin(ph), cos(ph))
  }
  L <- 100; us <- c(0.05, 0.1, 0.15)
  for (pat in names(TCP_PATTERN_GEOM)) {
    g <- TCP_PATTERN_GEOM[[pat]]; b <- verts(pat, L)
    expect_equal(ncol(b), g$dim, info = pat)
    v0 <- mean(diag(solve(fisher(b, rep(0, g$dim)))))
    expect_equal(L / sqrt(v0), g$a_ideal, tolerance = 1e-4, info = pat)
    f <- sapply(us, function(u) mean(apply(dirs(g$dim), 1, function(dv)
           mean(diag(solve(fisher(b, u * L * dv)))) / v0)))
    c_fit <- unname(coef(lm(I(f - 1) ~ 0 + I(us^2))))
    expect_equal(c_fit, g$c_off, tolerance = 0.05, info = pat)   # within 5%
  }
})

test_that("pattern names accept short forms and 2D vertex counts", {
  expect_equal(.tcp_pattern_name("6"), "hexagon")
  expect_equal(.tcp_pattern_name(" Hexagonal "), "hexagon")
  expect_equal(.tcp_pattern_name("3"), "triangle")
  expect_equal(.tcp_pattern_name("oct"), "octahedron")
  expect_equal(.tcp_pattern_name(""), "")
  expect_true(is.na(.tcp_pattern_name("7")))
  expect_identical(.tcp_geom("3"), TCP_PATTERN_GEOM$triangle)    # was a silent hexagon fallback
})
