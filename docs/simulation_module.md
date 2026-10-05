# MINFLUX Diffusion Simulator

*The ground-truth diffusion simulator in the Simulator tab: what it models and how to use it to test a
tracking result.* Companion to [MINFLUX_analysis_methods.md](MINFLUX_analysis_methods.md).

## 1. Purpose

An estimator can only be trusted on real data once it recovers a **known** answer. This module simulates a *known* diffusion process on a *known* geometry, samples it
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

## 1a. Workflow — is a molecule confined, or diffusing freely on a curved cell?

The comparison lives in a **round-trip**: the simulator makes ground-truth data in the *same format*
as real MINFLUX, so you push **both** real and simulated data through the **same Tracking-tab
analysis** and overlay the results. The simulator does not plot the comparison — the Tracking tab does.

1. **Analyse the real data.** Tracking tab → upload → apply your filters (split / trim / quality gate /
   length / duration). Run **MSD** (note the local D and whether the MSD bends over) and **Confinement**
   (per-track D, confinement index, apparent extent).
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
   confined are indistinguishable).

**Key caveat:** at D ≈ 0.1–0.2 µm²/s and ~20–100 ms tracks, a ~100 nm extent
is exactly what free diffusion produces, so step 6 is the whole point — do not read confinement off the
raw extent/MSD alone.

## 2. Architecture & integration

- New tab **"Diffusion / Tracking"** in `sim_ui.R` (alongside MINFLUX / Confocal / Acquisition /
  Compare strategies), driven from `sim_server.R`.
- **Reuse** the existing MINFLUX precision model: σ_loc = L_TCP / (a·√N_photons), a = √8 (`TCP_CRLB_PREFACTOR`,
  shared with the Tracking tab). Photon count N per localization ties to EFO·t_dwell; σ_loc is therefore
  **heteroscedastic** (varies per localization), like reality.
- **Output** the canonical frame `data.frame(tid, tim[s], x, y[, z] nm, efo, eco[, itr])` — identical
  to the importer's expectation — plus a **ground-truth sidecar** `data.frame(tid, D_true_um2_s,
  geometry, mode, state, L_c_true, v_true, …)`. Both download as CSV, and the tracks load into the
  Tracking tab like any other file.

## 3. Geometries

### Whole-cell (ROI contains the cell) — microbiology
One parametric **ellipsoid** + one **spherocylinder** cover all the requested shapes:

| Shape | Parameters | Notes |
|---|---|---|
| **Sphere** | radius R | ellipsoid with a=b=c=R |
| **Ovoid** | semi-axes a=b, c | prolate ellipsoid |
| **Rod** | cylinder radius r, cylinder length ℓ; caps radius r | spherocylinder; total length ℓ+2r |

### Partial (ROI < cell) — other systems
| Shape | Parameters | Notes |
|---|---|---|
| **Membrane plane** | Lx × Ly patch, orientation (nx,ny,nz), z-offset | a flat membrane imaged over a patch |
| **Free** | none | unbounded 3D diffusion |

## 4. What diffuses, and how

- **Surface diffusion** — on the ellipsoid/rod surface or the membrane plane. Overdamped Brownian
  motion in the local tangent plane, projected back onto the surface each step.
- **Volume diffusion** — 3D Brownian motion inside the shape, **reflecting** at the boundary.
- **Motion type** (applies on top of any geometry and mode; each is a ground-truth label):
  - **Free** — Brownian.
  - **Anomalous** — fractional Brownian motion, MSD ∝ t^α with correlated increments (Hurst = α/2).
    D sets the MSD at the single-cycle lag, so the CVE estimate is unchanged and only the growth with
    lag differs.
  - **Hop diffusion** — free diffusion inside corrals with permeable barriers.
  - **Confined** — an Ornstein-Uhlenbeck harmonic tether to a per-track home position
    (φ = e^{−t_c/τ_c}, τ_c = L_c²/D), giving an MSD plateau of 2·d·L_c².
  - **Directed** — drift velocity v along the rod axis, or a random direction on an ellipsoid.
  - **Two-state** — continuous-time switching between a mobile and a slow state at known rates, the
    ground truth for the HMM.
- **Survivorship** (optional) — a track ends when the true displacement across a dark gap exceeds a
  re-acquisition radius. This manufactures apparent confinement, which is what it is there to measure.

### 4.1 Surface-diffusion maths
Generate the continuous path at a fine internal step δt_sim (≪ cycle time), variance per axis = 2D·δt_sim.
- **Sphere**: tangent-plane step at point **p**, then retract `p ← R·(p+step)/|p+step|` (small-step
  accurate) or apply a random small rotation (exact).
- **Spherocylinder**: tangent step, then project back onto the cylinder body or the hemispherical cap.
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
   t_cyc + SW, one SW above the cycle time. The generator gives every cycle a duration of t_cyc and adds SW once for
   each *recorded* cycle, which reproduces that quantisation. The molecule diffuses through the SW
   too — putting SW on the clock but not in the step duration under-recovers D by t_cyc/(t_cyc+SW).
2. **Dark cycles → irregular Δt** — each attempted cycle is **integrated** (recorded) with probability
   p_int, else **dark** (skipped). Consecutive darks accumulate, so Δt between recorded localizations =
   k·t_cyc + SW, with k geometric (or empirically-calibrated). This is
   what produces the mixed dark-cycle structure across tracks.
3. **σ_loc (heteroscedastic)** — per recorded localization draw N_photons (∼ Poisson around the mean
   `eco`), σ_loc = **L_TCP/(√8·√N)**; add Gaussian noise N(0, σ_loc²) per axis. Record `eco=N`, `efo`.
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
   embedding) with MSD ∝ t^α. Increments are correlated, so they are drawn for the whole track at
   once. D keeps its meaning at the single-cycle lag.
4c. **Brightness** — per-molecule lognormal offset, exponential fade over the molecule's lifetime, and
   an AR(1) wobble within the track. Drives eco, σ_loc, efo and (under the photon budget) the gaps.
4d. **Photon budget** — the pattern repeats until phtLimit photons are collected, each repeat costing
   a cycle. Dark cycles, the eco floor, the τ quantisation and the rising gap fraction along a track
   all follow from this one rule.
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
   mechanism.
6. **Survivorship / re-acquisition radius R_reacq** — if the molecule moves > R_reacq during a dark
   gap, the track **terminates** (lost). *This reproduces the dark-cycle survivorship bias* (long-Δt
   steps preferentially retain small displacements → manufactured apparent confinement). Off by default,
   toggle-on to measure its size.
7. **Projection** — output (x,y) for 2D or (x,y,z) for 3D. For a surface protein on a curved cell, the
   2D projection distorts the transverse coordinate — the artifact the simulator exists to study.

## 6. Output

- **Localizations**: `data.frame(tid, tim[s], x, y[, z] nm, efo, eco[, itr])` — importer-ready.
- **Ground truth sidecar**: `data.frame(tid, D_true_um2_s, geometry, mode, state_seq?, L_c_true?,
  v_true?, n_dark_frac, sigma_loc_med)` for scoring.
- Buttons: *Tracks (CSV)* and *Ground truth (CSV)*. Load the tracks CSV into the Tracking tab to analyse them.

### 6a. Built-in CVE round-trip (`sim_trk_cve`, `sim_server.R:3128`)

The simulator scores itself: it runs the same CVE estimator the Tracking tab uses over the tracks it
just generated, so the recovered D can be compared against the D that was dialled in. That is the
validation loop of §1 closed in one number.

**What it computes.** Only **single-cycle steps** are used (`|Δt − t_c| < 0.25·t_c`), which
is the same integrated-pair restriction the real CVE applies — mixing multi-cycle steps would inflate
D. Then, per Vestergaard et al. 2014:

> D = ( ⟨|Δr|²⟩ + 2·⟨Δrᵢ·Δrᵢ₊₁⟩ ) / (2·d_mot·Δt),  Δt = t_c + SW
> σ_loc² = ( R·⟨|Δr|²⟩ + (2R−1)·⟨Δrᵢ·Δrᵢ₊₁⟩ ) / n_cart,  R = f·t_c/(6·Δt) with blur (f the collecting
> fraction of the cycle), 0 without

D is blur-independent. σ_loc is not, and reading a blurred acquisition with R = 0 (σ² = −cov₁)
under-reports it roughly threefold.

> **The σ half is a unit test.** Simulated positions use the same photon-budget σ the Tracking tab
> applies (see `sigma_loc_model.md`), so recovering σ only confirms the forward and inverse arithmetic
> agree. Each simulated position is drawn independently around the truth, while the real TCP recentres
> on its own noisy estimate and turns part of the error into a correlated wander. A σ round-trip can
> therefore pass here and still fail on the microscope. The D round-trip is the informative one.

The consecutive-step covariance is negative for localization noise (an error displaces step *i* one way
and step *i+1* the other), which is what lets the same two moments separate D from σ_loc without an
MSD fit.

**The load-bearing subtlety: d_mot is not n_cart.** The denominator of D uses the **motion**
dimensionality — 2 for surface diffusion, 3 for volume — while σ_loc is normalised by the
number of **Cartesian axes observed**. For a protein diffusing on a membrane or a rod surface these
differ: the motion is 2D but you observe x, y and z. Using n_cart = 3 in the D denominator would
under-report D by a factor of 3/2 on exactly the geometry most MINFLUX membrane work targets. The
Tracking tab corrects for this through its **Sample geometry** setting, and the round-trip verifies
the correction: without it the recovered ratio is 2/3.

**How to read it.** A round-trip ratio near 1.0 means the estimator is unbiased for that geometry,
sampling rate and σ_loc. A systematic departure is the *geometric bias* of §7. It belongs to the
estimator and the geometry, so it transfers to real acquisitions with the same shape and cycle time.

**Caveat.** The round-trip validates the estimator against the *simulator's* motion model. It cannot
detect a mismatch between that model and real biology — if a real protein hops between domains and
the simulation assumes free surface diffusion, a perfect round-trip still says nothing about the real
D. It bounds estimator error, not model error.

## 7. What each simulation tests

| Simulate | Tests |
|---|---|
| Free surface or volume Brownian, known D | CVE, cycle-binned MSD and OLSF recover D under irregular Δt and noise |
| Confined (OU tether) or hop diffusion, known L_c | plateau detection and the confinement index |
| Directed (drift v) | drift velocity, directionality, local α(t) |
| Two-state switching, known rates | HMM state and rate recovery |
| Survivorship on or off | size of the confinement bias that dark-cycle losses manufacture |
| Rod-surface diffusion projected to 2D | geometric confinement against true corralling |

The simulator is a validation tool: it never alters the standard analysis path.

## 8. Default parameters

The Diffusion / Tracking defaults describe a 2D tracking acquisition with the stock 2D tracking sequence.

| Field | Default | Derived from |
|---|---|---|
| TCP diameter L | 100.8 nm | sequence .json, final iteration `patGeoFactor` 0.28 × 360 nm |
| Mean photons/loc | 38 | median `eco` |
| Cycle time | 0.192 ms | spacing of the Δt levels |
| SW offset | 0.011 ms | intercept of Δt = k·t_c + SW |
| Dark-cycle prob. | 0.28 | Σ(k−1)/Σk on filtered data |
| Median locs/track | 100 | median track length |
| # tracks | 300 | |
| D | 0.3 µm²/s | per-track CVE median |

σ_loc follows from L and N as L/(√8·√N), the ideal value for a planar pattern.

## 9. Calibrate from a file

*Calibrate from real MINFLUX data* fills the acquisition fields from a loaded tracking file:

| Field | Measured as |
|---|---|
| cycle time, SW | regression of Δt on cycle count |
| photons | median `eco` |
| photon limit | 0.1st percentile of `eco` |
| collecting fraction | (eco/efo) ÷ cycle time |
| cycles before giving up | 99.9th percentile of τ/τ₀ |
| brightness spread, between / within | SD of per-track median log efo / median within-track SD |
| brightness carry-over | median within-track lag-1 autocorrelation of log efo |
| fading over a track | pooled log drop, first decile to last |
| length spread | mean/median of track lengths, inverted through the generator |
| D–lifetime correlation | Spearman(log length, log per-track D), de-attenuated ×2 |
| D | CVE on single-cycle steps |

Fading is pooled over localizations: short tracks mostly end by being lost, so averaging over tracks
would dilute it. The D–lifetime correlation is
filled de-attenuated, because both quantities are noisy and a simulation set to ρ reproduces a weaker
correlation; the report shows the measured and the set value side by side.

Calibrate cannot fill the **TCP diameter L**, which is stored in the sequence .json.

## 10. Limitations

- **Fading and `eco` dispersion trade off.** A fade that reproduces the observed `efo` decline
  compresses the `eco` variance, so the model cannot match both at once.
- **No `cfr` or `dcr`**, and `efo` is written as N/t_cyc, so simulated tracks cannot exercise the
  quality filters.
- **No background or mislocalized points**, so no test of what a stray localization does to MSD, JDD
  or confinement.
- **No re-acquisition model.** A molecule lost and re-found keeps its track here, while on the
  instrument it gets a new trace ID; only the survivorship cut is modelled.
- **No tracking feedback loop**, so the correlated wander of a real TCP (§6a) is absent.
