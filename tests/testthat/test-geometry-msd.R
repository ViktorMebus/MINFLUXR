# CLAIM (docs/MINFLUX_analysis_methods.md, MSD geometry reference): the free-diffusion
# MSD simulated on a closed surface is the analytic surface-diffusion curve, so the
# green "free on geometry" overlay is a correct null.
# Sphere surface (3D observation): MSD(t) = 2 R^2 (1 - exp(-2 D t / R^2)).
# The x,y projection is 2/3 of that by isotropy.
test_that("sphere-surface MSD matches the analytic limit (3D and 2D projection)", {
  skip_if_not_slow()
  set.seed(11)
  R <- 500; D <- 200; dt <- 0.2
  lags <- c(1, 10, 50, 200, 500, 1000, 2500, 5000)
  for (proj in c(FALSE, TRUE)) {
    reps <- replicate(8, sim_free_geom_msd(lags, D_surf = D, geom = "sphere",
                                           params = list(R = R), cyc_ms = dt,
                                           proj2d = proj, mode = "surface", M = 1500)$msd)
    theo <- 2 * R^2 * (1 - exp(-2 * D * lags * dt / R^2)) * (if (proj) 2/3 else 1)
    rel  <- rowMeans(reps) / theo - 1
    expect_true(all(abs(rel) < 0.03), info = paste("proj2d =", proj))   # measured < 1.4%
  }
})
