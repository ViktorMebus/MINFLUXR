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

1. **Cycle clock** — base cycle time t_cyc (e.g. 0.5–1 ms). A localization is *attempted* each cycle.
2. **Dark cycles → irregular Δt** — each attempted cycle is **integrated** (recorded) with probability
   p_int, else **dark** (skipped). Consecutive darks accumulate, so Δt between recorded localizations =
   (1 + n_dark)·t_cyc, with n_dark from a geometric (or empirically-calibrated) distribution. This is
   what produces the mixed dark-cycle structure across tracks.
3. **σ_loc (heteroscedastic)** — per recorded localization draw N_photons (∼ Poisson around EFO·t_dwell),
   σ_loc = L_TCP/(2·√N); add Gaussian noise N(0, σ_loc²) per axis. Record `eco=N`, `efo`.
4. **Motion blur** — R≈0 (sample the instantaneous position at the cycle time). ⚠ The optional nonzero-R path below is **not implemented**; the sampling overlay has no blur step. Optional:
   average the path over the acquisition sub-window for a small nonzero R.
5. **Photobleaching → finite length** — end the track after a photon-budget-limited number of
   localizations (geometric/exponential), so track lengths are short and realistic.
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
