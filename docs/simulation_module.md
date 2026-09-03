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
4b. **Anomalous motion** — optional fractional Brownian motion (`.sim_fgn`, Davies-Harte circulant
   embedding) with MSD ∝ t^α. Increments are correlated, so they are drawn up-front for the whole
   track rather than step by step. D keeps its meaning at the single-cycle lag.
4c. **Brightness** — per-molecule lognormal offset, exponential fade over the molecule's lifetime, and
   an AR(1) wobble within the track. Drives eco, σ_loc, efo and (under the photon budget) the gaps.
4d. **Photon budget** — the pattern repeats until phtLimit photons are collected, each repeat costing
   a cycle. Dark cycles, the eco floor, the τ quantisation and the rising gap fraction along a track
   all follow from this one rule. See §11.
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

D is blur-independent. σ_loc is not, and reading a blurred acquisition with R = 0 (σ² = −cov₁)
under-reports it roughly threefold.

> **Scope of the σ_loc half (2026-09-03).** The Tracking tab no longer takes σ_loc from the CVE
> (see `sigma_loc_model.md`), which narrows what this part of the round-trip can tell you.
> Simulated positions are perturbed with the same photon-budget expression the Tracking tab applies,
> so recovering σ from them checks that the forward and inverse arithmetic agree, and says nothing
> about whether that expression fits a real microscope. It is also gentler than the instrument: each
> simulated position is drawn independently around the truth, with no tracking feedback loop, whereas
> the real TCP is recentred on its own noisy estimate and turns part of the error into a correlated
> wander (log-log MSD slope 0.80 on immobile cells). That correlation is what defeats covariance-based
> σ estimation on real data, so a σ round-trip can succeed here and still fail on the microscope.
> Treat the σ number as a unit test under ideal noise. **The D round-trip is unaffected and now
> carries more weight, since D is the one quantity the CVE still supplies.** Adding a recentring loop
> to the simulator is the change that would make the σ half informative again.

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

## 11. Audit against GP67H-track-014 (2026-08-26)

Every statistic below is measured the same way on the real file and on simulator output.
GP67H is a cell-wall protein on the E. coli surface, so the reference model is surface diffusion on a
rod with little freedom.

| statistic | REAL filtered | REAL raw | SIM before | SIM after |
|---|---|---|---|---|
| track length, mean/median | 5.79 | 21.9 | 1.56 | 4.73 |
| longest 5% hold | 66% of locs | 75% | 24% | 49% |
| eco var/mean | 18.9 | 16.0 | 1.00 | 15.2 |
| eco floor | 20 (= phtLimit) | 20 | none | 20 |
| σ_loc q10–q90 | 2.87–5.25 nm | 3.13–5.50 | 3.72–4.60 | 2.85–5.63 |
| per-track D | 0.329 µm²/s | 0.335 | 0.219 | 0.31 (GP67H preset) |
| MSD α (3–30 cycles) | 0.78 | 0.67 | 0.97 | 0.88 |
| ⟨Δx²⟩ decline k=1→6 | ×3.0 | ×3.0 | ×1.05 | ×1.44 |
| length-vs-D (Spearman) | −0.42 | −0.36 | 0.00 | −0.18 |
| field of view | 9.2 × 3.7 µm, 2 cells | 10.8 × 5.6 | 2.4 × 0.8, 1 cell | multi-cell |
| cell size | ~4.0–4.9 µm long, ~0.96 µm wide | | 2.8 × 0.8 µm | 3.96 × 0.96 µm |

### Fixed in this pass

- **Photon model.** Real `eco` has a hard floor at exactly 20 — the sequence's `phtLimit` — with a long
  right tail from brightness differences (var/mean 19, not 1). The old Poisson draw could not produce
  either. Now the rate is lognormal and the count Poisson around it, floored at the photon limit, which
  reproduces the measured σ_loc spread almost exactly. This matters because σ_loc heterogeneity was
  understated roughly two-fold.
- **Track-length heterogeneity.** Real lengths are extremely heavy-tailed: the longest 5% of tracks
  carry two thirds of all localizations, so every pooled estimator is dominated by a handful of
  molecules. A single geometric lifetime cannot do this. Per-molecule lifetime spread (CV ≈ 3) gets
  most of the way there.
- **Collection window τ, measured not assumed.** Real files report τ directly as **eco/efo**: 0.1061 ms
  of a 0.192 ms cycle, i.e. **55%**, quantised in units of 0.1061 ms. So R = τ/(6Δt) ≈ **0.087**, not
  1/6. τ is now an input, defaulted to 0.55.
- **Gapped localizations integrate for longer.** τ before a single-cycle step is 1 repeat; before a
  two-cycle step it is 2 (98% of cases); it saturates near 3. **A dark cycle is not a skipped cycle —
  it is a longer integration.** The simulator now extends the averaging window across a gap, up to a
  repeat cap.
- **Multi-cell field of view** with random placement and in-plane orientation, and rod defaults moved
  to the measured cell size (radius 480 nm, cylinder 3000 nm).

### What the k-curve now looks like

Longer integration over gaps explains part of the flat k-curve but not all of it. A joint blur-aware
fit using the *measured* per-step τ over k = 1…6 gives D = 0.05–0.06 µm²/s and σ = 6.2 nm at R² = 0.78
— a poor fit, and four times below the single-cycle CVE D. Blur is a contributor, not the cause. τ
itself saturates near 3 repeats while the free-diffusion growth keeps going, so the deficit at k ≥ 4
survives. Censoring remains the leading explanation for the remainder (§10).

### Second pass: brightness, the photon budget, and anomalous diffusion

**The dark-cycle probability is not a free parameter.** `efo` falls monotonically along a track — to
44–56% of its starting rate by the last decile — and the gap fraction climbs with it, from 18% at the
start of a track to 43% at the end. Both follow from one mechanism: **the tracker repeats the pattern
until phtLimit photons are in, and each extra repeat costs a cycle.** A bright molecule finishes in one
repeat; a faded one needs two or three, and those extra cycles are exactly what appears in the file as
a "dark cycle". This is now the default (`Dark cycles come from the photon budget`), and p_dark becomes
an *output* rather than an input.

| | REAL filtered | SIM photon budget |
|---|---|---|
| eco median / mean | 36 / 45.3 | 36 / 45.5 |
| eco var/mean | 18.9 | 18.7 |
| eco floor | 20 | 20 |
| τ = 1 / 2 / 3 repeats | 0.72 / 0.23 / 0.05 | 0.74 / 0.23 / 0.04 |
| dark-cycle fraction | 0.28 | **0.26, emergent** |
| efo fade over a track | to 56% | to 52% |
| gap fraction, first → last decile | 0.19 → 0.39 | 0.11 → 0.46 |

The brightness model behind it: a per-molecule lognormal offset (sdlog 0.45), an exponential fade over
the molecule's lifetime, and an AR(1) wobble within the track (sdlog 0.45, lag-1 correlation 0.4) —
all three measured off the real files. Because eco drives σ_loc, σ_loc now also grows along a track,
as it does in the real data.

`efo` is written as N/τ rather than N/cycle, so simulated files carry the same τ = eco/efo
quantisation the real ones do and can be read back with the same measurement.

**Anomalous diffusion** is implemented as fractional Brownian motion (Davies–Harte circulant
embedding, exact), parameterised by α with Hurst = α/2. D still sets the MSD at the single-cycle lag,
so only the growth with lag changes. Round-trip through the whole sampling overlay, 2D projection
included: α = 0.78 in → 0.75 out.

This exposes an estimator failure worth carrying to the tracking tab: **on subdiffusive tracks the CVE
mis-attributes the anti-correlated fBm increments to localization noise.** Across all geometries at
α = 0.78 it returns D about 65% of truth while inflating σ_loc from 4.4 to ~4.9 nm. Both errors point
the same way the real data does, which is a reason to be careful about reading a low CVE D on GP67H as
a slow molecule rather than a subdiffusive one.

### Hop diffusion, and the length-vs-D coupling

**Hop diffusion** is now the second route to α < 1: the molecule diffuses freely inside a corral,
reflects off its barrier, and crosses into the neighbour with a set probability per encounter. Unlike
fBm it has no long-range memory — the subdiffusion comes from a crossover between free motion inside
a corral and slower motion between them. Measured α over 3–40 cycles:

| corral radius / permeability | 40 / 0.05 | 60 / 0.05 | 60 / 0.12 | 60 / 0.20 | 100 / 0.10 |
|---|---|---|---|---|---|
| α | 0.61 | 0.69 | 0.81 | 0.88 | 0.86 |

So a 60 nm corral at ~12% permeability sits on the measured GP67H α of 0.78, and is at least as
plausible for a cell-wall protein in a peptidoglycan mesh as fBm is. **The two are not distinguished
by α alone** — telling them apart is a question for the tracking tab, not the simulator.

**The length-vs-D anticorrelation** is now an explicit knob: D and lifetime are drawn as a correlated
pair. It is a phenomenological coupling, not a mechanism, and that is deliberate — the mechanisms that
could produce it (censoring, or a genuinely bound-versus-mobile population) are separately available
and neither reproduces the observed size on its own.

| D–lifetime correlation | 0 | −0.4 | −0.7 | −0.9 | −0.7 + censoring |
|---|---|---|---|---|---|
| recovered Spearman(length, D) | +0.01 | −0.23 | −0.24 | −0.33 | −0.05 |

Real is −0.42. The recovered value is attenuated because both quantities are noisy — short tracks give
noisy per-track D — so the latent coupling has to be much stronger than the number it produces. Worth
noting from the last column: **censoring does not add to the anticorrelation, it destroys it**, because
truncating the fast tracks compresses the length range that the correlation lives on.

Why it matters at all: this coupling is what decides whether a pooled D and a per-track median D
agree. Pooling weights by step count, so when the slow molecules are the long ones they dominate the
pooled number. On GP67H that gap is 0.22 against 0.33, a 40% disagreement between two defensible ways
of reporting one D.

### Remaining gaps

- **The k-curve is still only half reproduced.** Real declines ×3.0 from k=1 to k=6; the closest
  preset manages ×1.5. Blur over gaps, censoring, heterogeneity, hop diffusion and the photon budget
  have all been tried; each contributes and none closes it.
- **Full anticorrelation strength.** −0.33 against −0.42, and only with a near-maximal coupling.
- `cfr` and `dcr` remain absent — deliberately out of scope.

### A GP67H-like parameter set

Rod surface, radius 480 nm, cylinder 3000 nm, 2 cells; D 0.45 µm²/s with CV 0.35; confined,
L_c 60 nm; cycle 0.192 ms, SW 0.011, p_dark 0.28; photons 38, limit 20, spread 0.55; τ 0.55 with
3 repeats over gaps; lifetime CV 3; survivorship on at R_reacq 60 nm; project to 2D. This lands
per-track D at 0.31 against 0.329, Rg at 53 nm against 46, and α at 0.88 against 0.78.

## 12. Control-surface review (2026-08-26)

The panel had grown to 47 controls. Each was put to one test: **does the user both have a reason to
change it and a way to know what to set it to?** A control that fails the second half is not a
control, it is a constant with a text box around it.

That splits them three ways.

**The question you are asking** — geometry and cell size, number of cells, surface vs volume, motion
model and its parameter, D and its spread, track count and length, the D–lifetime coupling, 2D
projection, seed. These stay in the open. They are the reason to run the simulator at all.

**The microscope run** — cycle time, overhead, photons, photon limit, collecting fraction, cycles
before giving up, the four brightness parameters, dark-cycle probability. Every one is measurable from
a real file, so **Calibrate now fills all of them** and they live in a collapsed group. A user should
not be typing a number the data already knows.

**Instrument artifacts** — motion blur, survivorship, post-gap shrinkage. Deliberate distortions with
a known size, used to measure what an estimator does when the data misbehaves. Collapsed and labelled
as such, rather than mixed in with settings that describe reality.

### Removed

- **Brightness loss from motion.** A dimensionless coupling with no measurable value. The real
  evidence for it is a −0.14 correlation between per-track D and efo, it defaulted to 0, and nobody
  could know what to set it to. Deleted rather than left as decoration.
- **Dark probability rises with displacement.** Superseded by the photon budget, which is the actual
  mechanism, is measurable, and is on by default. Keeping a second, phenomenological route to the same
  effect only invites setting both.

### What Calibrate now reads from a file

| field | measured as | GP67H |
|---|---|---|
| cycle time, SW | regression of Δt on cycle count | 0.192 ms, 0.0103 ms |
| photons | median eco | 37 |
| photon limit | 0.1st percentile of eco | 20 |
| collecting fraction | (eco/efo) ÷ cycle time | 0.55 |
| cycles before giving up | 99.9th percentile of τ/τ₀ | 3 |
| brightness spread, between / within | SD of per-track median log efo / median within-track SD | 0.39 / 0.42 |
| brightness carry-over | median within-track lag-1 autocorrelation of log efo | 0.29 |
| fading over a track | pooled log drop, first decile to last | 0.59 |
| length spread | mean/median of track lengths, inverted through the generator | CV 3.0 |
| D–lifetime correlation | Spearman(log length, log per-track D), de-attenuated ×2 | −0.45 → −0.90 |
| D | CVE on single-cycle steps | 0.22 µm²/s |

Two of these needed care:

- **Fading must be pooled over localizations, not averaged over tracks.** The per-track median fade is
  0.12 where the pooled figure is 0.59, a factor of seven, because short tracks mostly end by being
  lost rather than by bleaching and so barely fade at all. The pooled figure is the one that
  reproduces the observable.
- **The D–lifetime correlation is filled de-attenuated.** Both quantities are noisy, so a simulation
  set to ρ reproduces a visibly weaker correlation; the report shows the measured and the set value
  side by side rather than hiding the adjustment.

Fading and eco dispersion trade off against each other in the current model: the fade that reproduces
the real efo decline (0.59–0.85) compresses eco var/mean below what the file shows, because faint
localizations get merged into multi-repeat windows that accumulate more photons. The model cannot hit
both at once, and that is a genuine limitation rather than a tuning failure.

**The one acquisition setting Calibrate cannot fill is the TCP diameter L**, because it lives in the
sequence .json rather than the data. It stays in the open with its provenance in the help text.

