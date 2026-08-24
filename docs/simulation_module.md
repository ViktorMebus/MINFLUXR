# MINFLUX Diffusion-Simulation Module — Design

*Design spec for a ground-truth diffusion simulator in the MINFLUXR simulator tab.*
Status: **implemented** (Phases 1-3; see §8 for what remains pending) · last reviewed 2026-08-06 · Companion to [MINFLUX_analysis_methods.md](MINFLUX_analysis_methods.md).

## 1. Purpose — close the validation loop

Every non-standard estimator we want (the OU-Kalman confined model, bacterial/geometry-aware
fits, continuous-time HMM) is **exploratory** and cannot be trusted on real data until it recovers a
**known** answer. This module simulates a *known* diffusion process on a *known* geometry, samples it
with **MINFLUX's own pathologies** (irregular Δt, dark cycles, per-localization σ_loc, photobleaching,
survivorship), and emits data in the **same format as real MINFLUX** so it flows through the **exact
same analysis pipeline**. The loop:

```
   ground truth (D, L_c, v, state, geometry)
        │  simulate continuous path
        ▼
   MINFLUX sampling overlay  ──►  MINFLUX-format table (tid,tim,x,y,z,efo,eco)
        │                              │  (+ ground-truth sidecar)
        │                              ▼
        └──────── score ◄──── run through the normal Tracking analysis
```

**Design rule:** the simulator must reproduce the constraints that break estimators — irregular Δt,
mixed dark cycles, heteroscedastic σ_loc, short tracks, and re-acquisition survivorship — otherwise it
validates nothing useful.

## 1a. Workflow — "is my (2D wall) protein confined, or just free diffusion?"

The comparison lives in a **round-trip**: the simulator makes ground-truth data in the *same format*
as real MINFLUX, so you push **both** real and simulated data through the **same Tracking-tab
analysis** and overlay the results. The simulator does not plot the comparison — the Tracking tab does.

1. **Analyse the real data.** Tracking tab → upload → apply your filters (split / trim / quality gate /
   length / duration). Run **MSD** (note the local D and whether the MSD bends over) and **Confinement**
   (per-track D, confinement index / Rg-Rmax, apparent extent).
2. **Export the filtered tracks.** Tracking tab → *Export filtered data (.csv)* — this is the
   post-filter `trk_spatial_filtered()` in canonical MINFLUX format (tim in seconds).
3. **Calibrate the simulator** on it. Simulator → Diffusion/Tracking → *Calibrate from real MINFLUX
   data* → load that CSV (min length 100). It fills cycle time, dark-cycle prob, photons, median length.
4. **Simulate the FREE null.** Geometry = **Rod**, long axis in the image plane; **Surface** diffusion;
   tick **Project to 2D** (your data is 2D); Motion = **Free**; set **D** to the real local D from step 1.
   Simulate, then *Download tracks (CSV)*.
5. **Analyse the simulated-free data** through the **same** Tracking-tab MSD + Confinement.
6. **Compare.** Overlay real vs simulated-free:
   - If the simulated **free** diffusion **reproduces** the real MSD shape / confinement index / extent →
     the apparent confinement is explained by geometry + 2D projection + finite short tracks. **You
     cannot claim confinement.**
   - If the **real data is tighter** (lower plateau, smaller extent, CI further below 1) than
     simulated-free → there is genuine confinement beyond geometry.
7. **Optional — quantify.** If real is tighter, simulate Motion = **Confined** and tune `L_c` until the
   simulated confinement statistics match the real ones; that L_c is your estimate. The summary's
   **resolvability** read-out warns when τ_c = L_c²/D exceeds the track duration (below which free and
   confined are indistinguishable — see §Phase 3).

**Key caveat:** at typical wall-protein numbers (D ≈ 0.1–0.2 µm²/s, ~20–100 ms tracks) a ~100 nm extent
is exactly what free diffusion produces, so step 6 is the whole point — do not read confinement off the
raw extent/MSD alone.

## 2. Architecture & integration

- New tab **"Diffusion / Tracking"** in `sim_ui.R` (alongside MINFLUX / Confocal / Acquisition /
  Compare strategies), driven from `sim_server.R`.
- **Reuse** the existing MINFLUX precision model: σ_loc = L_TCP / (2·√N_photons) (the sim already uses
  L_eff = 2·σ·√N). Photon count N per localization ties to EFO·t_dwell; σ_loc is therefore
  **heteroscedastic** (varies per localization), like reality.
- **Output** the canonical frame `data.frame(tid, tim[s], x, y[, z] nm, efo, eco[, itr])` — identical
  to the importer's expectation — plus a **ground-truth sidecar** `data.frame(tid, D_true_um2_s,
  geometry, mode, state, L_c_true, v_true, …)`. Downloadable, and loadable straight into the Tracking
  tab (or piped in-memory).

## 3. Geometries

### Whole-cell (ROI contains the cell) — microbiology
One parametric **ellipsoid** + one **spherocylinder** cover all the requested shapes:

| Shape | Parameters | Notes |
|---|---|---|
| **Sphere** (coccus) | radius R | ellipsoid with a=b=c=R |
| **Ovoid / coccobacillus** | semi-axes a=b, c | prolate ellipsoid |
| **Rod / bacillus** (*E. coli*) | cylinder radius r, cylinder length ℓ; caps radius r | spherocylinder; total length ℓ+2r |

### Partial (ROI < cell) — other systems
| Shape | Parameters | Notes |
|---|---|---|
| **Membrane plane** | Lx × Ly patch, orientation (nx,ny,nz), z-offset | a wall/membrane protein imaged over a patch |
| **Cytoplasmic slab / half-space** ⚠ *not implemented (see §8)* | box or single boundary plane | free 3D diffusion, one membrane, no far wall |

## 4. What diffuses, and how

- **Surface diffusion** — on the ellipsoid/rod surface or the membrane plane. *The cell-wall-protein
  case.* Overdamped Brownian in the local tangent plane, retracted onto the surface each step.
- **Volume diffusion** — 3D Brownian inside the shape, **reflecting** at the boundary. Cytoplasmic.
- **Modifiers** (composable; each becomes a ground-truth label):
  - **Corrals** ⚠ *not implemented.* The confinement modifier that WAS built is an
    **Ornstein-Uhlenbeck harmonic tether** to a per-track home position
    (φ = e^{−t_c/τ_c}, τ_c = L_c²/D) — see §8 Phase 3 and
    `MINFLUX_analysis_methods.md` §I-C.2. No corral tessellation exists.
  - **Directed** — drift velocity v along a direction (e.g. cell axis, MreB-like), optionally
    run-and-tumble state switching.
  - **Immobile fraction** — a subpopulation with D≈0.
  - **Multi-state** — a continuous-time Markov chain over {states}, each with its own (D, v); switches
    at a known rate matrix → ground truth for HMM / OU-HMM validation.

### 4.1 Surface-diffusion maths
Generate the continuous path at a fine internal step δt_sim (≪ cycle time), variance per axis = 2D·δt_sim.
- **Sphere**: tangent-plane step at point **p**, then retract `p ← R·(p+step)/|p+step|` (small-step
  accurate) or apply a random small rotation (exact).
- **Spherocylinder body**: **superseded.** The rod surface now uses **tangent-step +
  cylinder-project**, unified with the ellipsoid path (§8 Phase 3). The unrolled isometric chart
  below is kept as a historical note on why the flat treatment is valid, not as a description of the
  code: the cylinder is **developable** (zero Gaussian curvature), so unrolled coordinates
  (axial s ∈ [0, ℓ], arclength u = r·φ) are a **true isometric flat chart** —
  Brownian in flat (s, u), **periodic** in φ, and at s∈{0,ℓ} **transition onto the
  hemispherical cap**. Map (s,φ) → 3D: (r·cosφ, r·sinφ, s).
- **Ellipsoid**: nonzero curvature → no flat chart. Use tangent-step + **project back** onto the
  ellipsoid (normalize in scaled coordinates, or one Newton step). Accurate for small δt_sim.
- **Membrane plane**: flat 2D Brownian, embedded in 3D at the chosen orientation/offset.

### 4.2 Volume diffusion
3D Brownian; if a step exits the boundary, **specularly reflect** across the surface normal (analytic
inside-test and normal for ellipsoid & spherocylinder). Repeat until the step lands inside.

## 5. The MINFLUX sampling overlay (the load-bearing part)

Subsample the fine ground-truth path the MINFLUX way:

1. **Cycle clock + SW offset** — base cycle time t_cyc, plus a fixed per-localization overhead SW
   (settling / window). Real files obey **Δt = k·t_cyc + SW exactly**, so the *modal* Δt is
   t_cyc + SW, not t_cyc. The generator gives every cycle a duration of t_cyc and adds SW once for
   each *recorded* cycle, which reproduces that quantisation. The molecule diffuses through the SW
   too — putting SW on the clock but not in the step duration under-recovers D by t_cyc/(t_cyc+SW).
2. **Dark cycles → irregular Δt** — each attempted cycle is **integrated** (recorded) with probability
   p_int, else **dark** (skipped). Consecutive darks accumulate, so Δt between recorded localizations =
   k·t_cyc + SW, with k geometric (or empirically-calibrated). This is
   what produces the mixed dark-cycle structure across tracks.
3. **σ_loc (heteroscedastic)** — per recorded localization draw N_photons (∼ Poisson around the mean
   `eco`), σ_loc = **L_TCP/(4·√N)**; add Gaussian noise N(0, σ_loc²) per axis. Record `eco=N`, `efo`.
   L_TCP is the final-iteration TCP diameter, read off a sequence .json as `patGeoFactor × 360 nm`
   (see `TCP_L_PER_GEOFACTOR_NM`); the stock tracking sequences end at 0.28 → **100.8 nm**.
4. **Motion blur** — implemented, on by default. A localization is the *time-average* of the path over
   the photon-collection window τ (an input, as a fraction of the cycle). Conditioned on both end
   points of the step, that average is

   > mean = x_i + (τ/2Δt)·(x_{i+1} − x_i),  var = (2Dτ/3)·(1 − 3τ/4Δt) per axis

   which reduces to the midpoint plus D·Δt/6 when τ = Δt. The corresponding blur coefficient is
   **R = τ/(6·Δt)**: a full-cycle window gives 1/6, a short burst tends to 0. Exact for free / volume
   diffusion; on a curved surface the bridge noise is kept in the tangent plane (the average of two
   surface points sits marginally inside the surface). Switching blur off gives R = 0 snapshot sampling.
5. **Photobleaching → finite length** — end the track after a photon-budget-limited number of
   localizations (geometric/exponential), so track lengths are short and realistic.
5b. **Heterogeneous D** — each track may draw its own D from a lognormal whose median is the D input
   and whose spread is set by a CV input. 0 (the default) gives one D for the whole population. With
   spread on, the pooled CVE cannot recover the median: it recovers the step-weighted mean of the
   drawn values, and the summary scores the round-trip against that.
5c. **Speed-dependent dark cycles** — optionally p_dark = 1 − (1 − p_dark₀)·exp(−d/scale), with d the
   true displacement over that cycle, decided inside the path loop so it uses the realized step. This
   is the *population selection* mechanism: fast molecules dark-cycle more, so the surviving gapped
   steps oversample slow ones.
5d. **Re-acquisition shrinkage** — optionally the first localization after a dark gap is reported a
   fixed fraction of the way back toward the previous *reported* position, the estimator-bias
   mechanism. See §10 for what this predicts and why the real data rejects it in this form.
6. **Survivorship / re-acquisition radius R_reacq** — if the molecule moves > R_reacq during a dark
   gap, the track **terminates** (lost). *This reproduces the dark-cycle survivorship bias* (long-Δt
   steps preferentially retain small displacements → manufactured apparent confinement). Off by default,
   toggle-on to measure its size.
7. **Projection** — output (x,y) for 2D or (x,y,z) for 3D. For a surface protein on a curved cell, the
   2D projection distorts the transverse coordinate — the artifact we want to *study*, not hide.

## 6. Output

- **Localizations**: `data.frame(tid, tim[s], x, y[, z] nm, efo, eco[, itr])` — importer-ready.
- **Ground truth sidecar**: `data.frame(tid, D_true_um2_s, geometry, mode, state_seq?, L_c_true?,
  v_true?, n_dark_frac, sigma_loc_med)` for scoring.
- Buttons: download both; "Load into Tracking tab" (in-memory hand-off).

### 6a. Built-in CVE round-trip (`sim_trk_cve`, `sim_server.R:3128`)

The simulator scores itself: it runs the same CVE estimator the Tracking tab uses over the tracks it
just generated, so the recovered D can be compared against the D that was dialled in. That is the
validation loop of §1 closed in one number.

**What it computes.** Only **single-cycle steps** are used (`|Δt − t_c| < 0.25·t_c`), which
is the same integrated-pair restriction the real CVE applies — mixing multi-cycle steps would inflate
D. Then, per Vestergaard et al. 2014:

> D = ( ⟨|Δr|²⟩ + 2·⟨Δrᵢ·Δrᵢ₊₁⟩ ) / (2·d_mot·Δt),  Δt = t_c + SW
> σ_loc² = ( R·⟨|Δr|²⟩ + (2R−1)·⟨Δrᵢ·Δrᵢ₊₁⟩ ) / n_cart,  R = 1/6 with blur, 0 without

D is blur-independent; σ_loc is not. Reading a blurred acquisition with R = 0 (σ² = −cov₁) is what
makes a real 4 nm σ_loc read as 1.4 nm.

The consecutive-step covariance is negative for localization noise (an error displaces step *i* one way
and step *i+1* the other), which is what lets the same two moments separate D from σ_loc without an
MSD fit.

**The load-bearing subtlety: d_mot is not n_cart.** The denominator of D uses the **motion**
dimensionality — 2 for surface diffusion, 3 for volume — while σ_loc is normalised by the
number of **Cartesian axes observed**. For a protein diffusing on a membrane or a rod surface these
differ: the motion is 2D but you observe x, y and z. Using n_cart = 3 in the D denominator would
under-report D by a factor of 3/2 on exactly the geometry most MINFLUX membrane work targets. The app
already corrects for this (`server_tracking.R:8449-8470`, driven by the **Sample geometry** setting);
the round-trip is what *verifies* the correction rather than what discovers the need for it — the
"simulated: 0.675" in that source comment is this loop reporting 2/3, as predicted.

**How to read it.** A round-trip ratio near 1.0 means the estimator is unbiased for that geometry,
sampling rate and σ_loc. A systematic departure is the *geometric bias* §7 asks about — and
it is a property of the estimator-plus-geometry, not of the data, so it transfers to real acquisitions
with the same shape and cycle time.

**Caveat.** The round-trip validates the estimator against the *simulator's* motion model. It cannot
detect a mismatch between that model and real biology — if a real protein hops between domains and
the simulation assumes free surface diffusion, a perfect round-trip still says nothing about the real
D. It bounds estimator error, not model error.

## 7. Validation targets (what each simulation checks)

| Simulate | Checks |
|---|---|
| Free surface/volume Brownian, known D | CVE, cycle-MSD, OLSF recover D under irregular Δt + noise |
| Confined (small ellipsoid / corral), known L_c | plateau detection, Rg/CI, **OU-Kalman L_c** |
| Directed (drift v) | drift-velocity, directionality, local α(t) |
| Multi-state CTMC, known rates | HMM / continuous-time OU-HMM state & rate recovery |
| Dark-cycle survivorship on/off | size of the manufactured-confinement bias |
| 2D projection of a rod-surface protein | geometric confinement artifact vs true corralling |

## 8. Build phases

- **Phase 1 (first):** rod (spherocylinder) + sphere **surface** diffusion (free Brownian), the MINFLUX
  sampling overlay (cycle + dark cycles + heteroscedastic σ_loc + photobleaching), MINFLUX-format output
  + ground-truth D, and a 3D/2D preview. Validates CVE / cycle-MSD immediately.
- **Phase 2 (DONE 2026-07-12):** ellipsoid (ovoid/coccobacillus) geometry, and **volume /
  intracellular diffusion** for sphere/rod/ellipsoid — free 3D Brownian inside the cell,
  reflecting off the wall (radial reflect in scaled coords for ellipsoid; cylinder-wall +
  hemispherical-cap reflect for the rod). A **Surface / Volume** selector. The CVE round-trip uses
  d = min(motion-DOF, observed axes) — 2 for surface, 3 for volume, capped to 2 under 2D
  projection. Validated: D recovered 97–101% of truth across all 6 geometry×mode combos.
- **Phase 2b (DONE 2026-07-12):** **Membrane plane** (flat XY sheet with a thickness; 2D lateral
  diffusion, z wobbles within the thickness) and **Free diffusion (no boundaries)** (unbounded 3D).
  Both FLAT → recover D at 99–101% in *both* 3D and 2D (no projection bias) — the clean contrast to
  the curved cylinder's 75%. A **geometric bias factor** read-out (recovered/truth) is shown in the
  summary so the projection correction can be lifted directly. *Still pending:* cytoplasmic slab.
- **Phase 3 (DONE 2026-07-12):** motion modifiers on top of any geometry×mode via a unified
  geometry framework (project / reflect / normal). **Confined** = OU harmonic tether to a per-track
  home (validated: stationary SD per DOF = L_c exactly; total MSD plateau = 2·d·L_c²). **Directed**
  = drift v along the rod long axis / random for ellipsoid (validated: recovers v). **Two-state** =
  mobile↔immobile CTMC at a per-cycle switch prob (HMM ground truth). **Survivorship** = track ends
  if the true dark-gap displacement exceeds a re-acquisition radius (manufactures apparent
  confinement). Rod surface now uses tangent-step + cylinder-project (unified with ellipsoid), not
  the unrolled chart. KEY INSIGHT surfaced: confinement is only visible when τ_c = L_c²/D is shorter
  than the track — else free diffusion over a short track mimics the same ~extent (the artifact to
  test for the wall protein). *Still pending:* membrane plane / cyto slab; per-localization state
  export for full HMM scoring.

Everything is labelled a **simulation / validation tool**; it does not alter the standard analysis path.

## 9. Default parameters — where they come from

The Diffusion / Tracking defaults are the measured values of the GP67H tracking test data
(`251027-200445_GP67H-track-013_filtered…csv`, `251027-202407_GP67H-track-014_filtered…csv`) and of the
stock 2D tracking sequence (`Test data/Tracking Imaging Sequences/Tracking_2D_Oct2022.json`).

| Field | Default | Source |
|---|---|---|
| TCP diameter L | 100.8 nm | sequence .json, final iteration `patGeoFactor` 0.28 × 360 nm |
| Mean photons/loc | 38 | median `eco` (013: 38, 014: 36) |
| Cycle time | 0.192 ms | spacing of the Δt levels (013: 0.1939, 014: 0.1923) |
| SW offset | 0.011 ms | intercept of Δt = k·t_c + SW (013: 11.4 µs, 014: 10.4 µs) |
| Dark-cycle prob. | 0.28 | Σ(k−1)/Σk on filtered data (013: 0.26, 014: 0.29; unfiltered 014: 0.61) |
| Median locs/track | 100 | median track length of 014 filtered (013 is a 101-track subset with median 666) |
| # tracks | 300 | 679 tracks in 014 filtered, 342 of them ≥100 locs |
| D | 0.3 µm²/s | per-track CVE median (013: 0.28, 014: 0.33); the pooled CVE is 0.21–0.22, weighted toward the long, slower tracks |

σ_loc is *not* an input — it follows from L and N. The cross-check: L/(4·√38) = 4.09 nm, and the
blur-aware CVE (R = 1/6) on both real files gives 4.13 and 4.08 nm. The R = 0 form reads 1.5 and 1.3 nm
on the same data, which is the blur term, not a precision difference.

Sequence-file numbers not yet used by the simulator: `patDwellTime` 100 µs, `phtLimit` 20 (final
iteration), 4 iterations, hexagon pattern, `ctrDwellFactor` 0.16, `stickiness` 4, `damping` 2.
The dwell and photon limit together imply the cycle time and the photon budget, so they are a route to
predicting t_cyc and `eco` rather than calibrating them from data.

### Still missing from the sampling overlay

- **Photon budget as the cycle terminator.** N_photons is drawn independently of the cycle; on the
  instrument the cycle ends when `phtLimit` is reached *or* `patDwellTime` times out, which couples
  N, t_cyc and brightness. Dim molecules should get longer cycles, not just fewer photons.
- **`efo` / `cfr` / `dcr`.** `efo` is written as N/t_cyc and the other two are absent, so simulated
  tracks cannot exercise the quality filters (real medians: efo ≈ 3.0–3.2e5, cfr ≈ 0.40–0.62,
  dcr ≈ 0.64–1.0).
- **Track-to-track heterogeneity.** Every simulated molecule gets the same D. Real per-track CVE D
  spans 0.23–0.43 µm²/s (IQR) on these files, so an estimator that looks unbiased on the simulator can
  still be biased on a heterogeneous population.
- **Background / mislocalization.** No off-target or background localizations, so no test of what a
  stray point does to MSD, JDD or confinement.
- **Trace ID re-use and re-acquisition.** A molecule lost and re-found gets a new tid on the
  instrument; the simulator has no re-acquisition model, only the survivorship cut.

## 10. Which mechanism flattens the k-curve?

Real MINFLUX tracking files show per-axis ⟨Δx²⟩ growing far more slowly with the dark-gap length k
than free diffusion allows (GP67H: 100 / 75 / 53 / 42 / 36 / 30 % of expectation for k = 1…6, see
`tracking_sigma_and_cve_plan.md` §2.3). The simulator now carries the three candidate mechanisms, and
the summary panel prints the same k-curve, so they can be told apart. Sweeps at D = 0.3 µm²/s,
free geometry, 250 tracks:

| model | k = 1…6, % of free expectation | verdict |
|---|---|---|
| baseline | 99 100 104 100 107 100 | the null: flat, as it must be |
| heterogeneous D (CV 0.7) | 100 100 105 107 96 144 | **no effect on its own** |
| speed-dependent dark (scale 20 nm) | 81 81 88 103 115 114 | depresses *all* k alike, no gradient |
| gap shrinkage 0.5 | **118** 27 27 25 28 24 | crushes k ≥ 2 but *inflates k = 1* |
| survivorship R_reacq = 50 nm | 100 96 95 78 66 67 | **matches the real shape** |

Two conclusions:

- **Shrinkage in this form is falsified by the data.** Pulling a post-gap localization back toward
  the previous one displaces the reported track from the true one, and the following single-cycle step
  has to repay the difference — so k = 1 goes *up*. Real files show k = 1 fully consistent with the
  CVE D, so whatever compresses the gapped steps does not leave that footprint. The knob is kept
  because it is the honest implementation of the mechanism, and because a falsifiable prediction that
  fails is a result.
- **Censoring is the surviving candidate.** It is the only mechanism tested that leaves k = 1 intact
  and produces a *gradient* over k. The earlier objection — "tracks do not preferentially end after a
  large gap" — was not a valid test: a censored step is **absent** from the file, it does not appear
  as a track-ending step, so this signature is invisible by construction.

Censoring also reproduces the second real-data signature. With R_reacq = 50 nm and a heterogeneous
population, the step-weighted mean D drops to 0.231 against a median of 0.428, because fast molecules
have their tracks cut short and contribute fewer steps. That is the same direction and rough size as
the real pooled-vs-per-track gap (0.22 pooled against a 0.30 per-track median).

The caveat worth carrying: the R_reacq that reproduces the real k-curve is ~50 nm, well below the
323 nm search field of the sequence's first iteration. Either the effective re-acquisition range is
much tighter than the field, or censoring is not the whole story.
