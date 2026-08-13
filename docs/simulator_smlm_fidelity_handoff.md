# SMLM Simulator — MINFLUX Fidelity Handoff

*Findings from a camera-SMLM-assumption audit of the SMLM app and simulator (2026-08-11). The app
side is done; **this file is the un-started simulator work**, split out so it can be picked up in a
fresh session.* Companion to [simulation_module.md](simulation_module.md) (the diffusion/tracking
simulator, which is **already MINFLUX-faithful** and is the reference implementation here) and
[MINFLUX_analysis_methods.md](MINFLUX_analysis_methods.md).

---

## 0. One-paragraph brief for a cold start

The SMLM simulator (`sim_server.R`, the `gen_localizations` / CRLB path — **not** the `sim_trk_*`
tracking simulator) draws every localization with a **single constant σ_loc**, and its CSV export
carries **no `eco` / `efo` / `cfr`**. Real MINFLUX σ varies per localization with photon count, and
those three columns are the MINFLUX-native metadata the SMLM tab's own tools key off. The result is
that simulated SMLM data cannot exercise the MINFLUX-specific parts of the app it was built to
validate. The fix is small, and a correct implementation of the same idea already exists ~1800 lines
away in the same file.

---

## 1. What is wrong (ranked by value)

### 1.1 σ_loc is constant across every localization — **fix this first**

`gen_localizations()` (`sim_server.R:848`) takes scalar `sigma_loc_xy` / `sigma_loc_z` and applies:

```r
rep_df$x <- rep_df$x + rnorm(total, 0, sigma_loc_xy)
```

Real MINFLUX precision is set per localization by that localization's photon count. Consequences:

- Simulated precision clouds are **homogeneous Gaussians**. The app's own precision-cloud QC states
  that a pooled cloud of Gaussians with *different* widths is leptokurtic by construction, and that
  MINFLUX σ varies with photon count — it explicitly predicts kurtosis from the σ spread rather than
  testing against 3. Simulated data has no σ spread, so **it can never exercise that check
  realistically**, and anything tuned against it is tuned against the wrong null.
- Downstream, every inverse-variance weighting in the app becomes a no-op on simulated data, because
  all points really are equally precise.

**There is already a correct implementation in the same file.** The tracking simulator
(`sim_server.R:3097`) does exactly the right thing:

```r
N   <- pmax(1L, stats::rpois(length(idx), photons))   # per-localization photon draw
sig <- L_tcp / (2 * sqrt(N))                          # per-localization sigma
obs <- pos[idx, ] + matrix(rnorm(length(idx) * 3), ncol = 3) * sig
```

Port that pattern into `gen_localizations`. `simulation_module.md` §1 lists "per-localization σ_loc"
as a design goal of the tracking simulator; the SMLM simulator was simply never held to it.

### 1.2 The CSV export drops the MINFLUX-native columns

`output$sim_dl_csv` (`sim_server.R:2673`) exports only `x, y, z, tid, tim` plus ground truth. The
tracking export emits `efo` and `eco` (`sim_server.R:3105`). What silently breaks when simulated SMLM
data is loaded back into the SMLM tab:

| Feature | Default | Behaviour without `eco`/`efo` |
|---|---|---|
| Shape-fit photon weighting (`viewer_weight_eco`) | **ON** | Silently falls back to unweighted |
| Fiducial / contaminant scan | — | Its primary BRIGHT signature is `efo`; cannot run |
| efo / cfr metadata filters | — | Inert |

Fixing 1.1 produces the per-localization photon count `N`, so `eco = N` and
`efo = N / dwell_seconds` fall out for free. `cfr` has no counterpart in the current model — see the
open question in §3.

### 1.3 Background is modelled as a plain photon loss

CRLB mode (`sim_server.R:1263`, preview at `sim_server.R:2212`) uses:

```
N_eff = N_photons × (1 − bg_fraction)
sigma = L / (2 * sqrt(N_eff))
```

This is the one place a **camera-SMLM intuition genuinely does not transfer**. On a camera, background
adds noise and effectively costs you photons. In MINFLUX the position information lives in the
*contrast at the doughnut zero*, and background **fills that zero in** — so the CRB degrades faster
than the photon loss alone implies. Treating background as `N × (1 − bg)` is therefore optimistic,
increasingly so at high background.

Impact is limited to CRLB mode used *predictively*. If L is obtained via "Calibrate from a real run",
the error is largely absorbed into the calibrated L (see §2).

---

## 2. Do NOT change these — they were checked and are correct

Recorded so a future session does not re-litigate them.

- **The `L / (2√N)` prefactor.** The literature convention varies with pattern geometry (1D vs 2D,
  number of beam positions, on- vs off-axis). It does not matter here, because
  **"Calibrate from a real run" inverts the same formula** (`L ≈ 2σ√N`, `sim_ui.R:532`). The
  convention cancels as long as both directions use it. Do not "correct" the 2 to a 4 in isolation —
  that would silently rescale every previously calibrated L.
- **`sim_crlb_z_factor` default 0.7** (σ_z = 0.7 × σ_xy, i.e. **z tighter than xy**). This is correct
  for 3D MINFLUX and is *not* the camera-SMLM astigmatism reflex. See
  `reference_minflux_axial_precision` in memory. It is also correctly labelled as a precision ratio
  and is not confused with the separate 0.7 refractive-index axial calibration.
- **The confocal PSF model** (`compute_confocal_sigma`, `sim_server.R:899`): `0.21·λ/NA`, Airy
  pinhole, Bessel encircled-energy. This is a genuine *detection* PSF for the confocal comparison
  view, so "PSF" is the right word there and the physics is right.
- **Variance ∝ 1/photons.** The MINFLUX CRB is σ ∝ L/√N; only the prefactor differs from camera SMLM
  (targeting-pattern size instead of spot width). The 1/√N scaling is shared, so photon weighting is
  sound.

---

## 3. Open questions for the user (decide before implementing)

1. **Should the simulator emit `cfr`?** `cfr` reports how well centred the molecule is in the
   doughnut. Emitting it honestly requires modelling the molecule's position *within* the TCP and
   the resulting position-dependent precision degradation — a bigger change than 1.1. Options:
   omit it (current), emit a constant placeholder (risks being read as real), or model it.
   **Recommendation: omit, and state that simulated data cannot exercise cfr-based tools.** A
   placeholder is the worst option.
2. **Should σ vary per localization only, or also per trace?** Photon count varies both within a
   trace and between molecules (labelling, depth, local background). 1.1 as written gives
   within-trace variation from the Poisson draw alone.
3. **Should background degrade the CRB non-linearly (1.3), or stay a documented simplification?**
   If CRLB mode is only ever used with "Calibrate from a real run", leaving it is defensible —
   but then the UI should say so.

---

## 4. How to verify a change

Patterns that worked in this session:

- Source the helpers into an env and call the functions directly, rather than launching the app:
  `e <- new.env(); sys.source("helpers_smlm.R", envir = e)`. For UI checks, `sys.source("global.R")`
  first (constants like `MIN_LOCS_FOR_SPREAD_DEFAULT`) and load `shinyjs`.
- **The acceptance test for 1.1** is the app's own precision-cloud QC: simulate, export, load into the
  SMLM tab, run the precision cloud. Before the fix, the pooled cloud is a clean Gaussian. After it,
  the σ spread should make it **leptokurtic**, and the panel's predicted-kurtosis check should sit
  near its expectation rather than firing spuriously. That closes the loop the simulator exists for.
- Round-trip check for 1.2: export CSV, confirm `eco`/`efo` are present, load it, confirm
  `viewer_weight_eco` reports a non-trivial `n_eff` rather than falling back.
- Parse check every touched file: `Rscript -e "invisible(parse('sim_server.R'))"`.

---

## 5. Scope boundary

This file covers **only** the SMLM simulator's MINFLUX fidelity. Not in scope, and being handled
separately: Median Clustering, the localization-spread / precision-cloud panel, and sample-drift
estimation. The tracking/diffusion simulator is already MINFLUX-faithful — see
[simulation_module.md](simulation_module.md) — and is the reference to copy from, not to change.
