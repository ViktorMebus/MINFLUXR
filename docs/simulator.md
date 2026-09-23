# MINFLUX Simulator — Reference

*Reference for the Simulator tab: what it models, how each stage works, and where it currently
departs from real MINFLUX.* Last reviewed 2026-08-11; §7 intro and §8 status checked against the code 2026-09-17 · Code: [`sim_server.R`](../sim_server.R),
[`sim_ui.R`](../sim_ui.R) · Companions: [simulation_module.md](simulation_module.md) (Phase 2 design
spec), [MINFLUX_analysis_methods.md](MINFLUX_analysis_methods.md) (the analysis side),
the 2026-08-11 fidelity audit, now superseded by §7 and §8 here.

---

## 0. What the simulator is for

MINFLUX localizes with an **excitation doughnut**: the molecule is found by minimising fluorescence at
the doughnut zero, and each *trace* (`tid`) is one molecule held and re-localized many times. A trace
therefore comes out as a **point cloud**, not a single point, and the cloud width — not any optical
PSF — is what sets the achievable resolution.

The simulator exists to answer, forward:

> *"If my structure looks like **X** and I label it with **Y**, how much acquisition time do I need
> before the structure is resolved?"*

It does that by building a known ground truth, pushing it through a labeling model, a detection
model and a localization-precision model, and emitting MINFLUX-shaped output that can be fed back
through the app's own analysis. Two independent simulators live in the tab:

| Sub-tab | What it simulates | Status |
|---|---|---|
| **SMLM** (Phase 1) | Static structures → labels → traces → localization point clouds; plus an acquisition-timeline model | Production tool; fidelity gaps in §7 |
| **Diffusion / Tracking** (Phase 2) | Brownian motion on a cell geometry, sampled the MINFLUX way | MINFLUX-faithful; see [simulation_module.md](simulation_module.md) |

Phase 2 is the **reference implementation** for MINFLUX sampling fidelity. Where Phase 1 and Phase 2
disagree about how MINFLUX works, Phase 2 is right.

---

## 1. Phase 1 — the SMLM generative pipeline

One pass of `sim_result()` ([`sim_server.R:1000`](../sim_server.R#L1000)) runs seven stages. Each is a
separate, inspectable model; understanding which stage a parameter acts on is most of understanding
the tab.

```
  (1) ground-truth molecules          gen_* / gen_from_image        → gt
        │  positions of the real protein, nm
  (2) labeling                        expand_to_fluorophores()      → fluoro
        │  efficiency, N dyes/molecule, epitope offset, linker, dye spread
  (3) TID multiplication              split_for_blinking()  [Normal]
        │  or  rate × duration sampling of the pool  [Acquisition]
  (4) precision model                 manual σ  |  CRLB σ = L/(2√N_eff)
        │
  (5) localizations per trace         gen_localizations()           → locs
        │  Poisson or geometric trace length, σ added per localization
  (6) background TIDs                 uniform over the FOV
        │
  (7) scoring                         RMSE vs GT, steric check, precision budget
```

### 1.1 Ground truth (stage 1)

Nine generators, all centred on the origin, 2D or 3D
([`sim_server.R:374-657`](../sim_server.R#L374)): cluster, ring, arc, double ring, ellipse/ellipsoid,
hollow sphere, filament (with optional helical twist and random-walk lateral curl), oligomeric ring
stack, and **image-guided** (`gen_from_image`) which samples molecule positions with probability
proportional to pixel intensity of a loaded confocal image.

Each generator may attach a **label-direction hint** used later by the linker model:

- `ldir_type = "perp"` (filament) with a per-molecule tangent `ldir_x/y/z` → the label is displaced
  *perpendicular* to the filament.
- `ldir_type = "radial"` (ring, sphere, arc) with a structure centre → the label is displaced
  *outward*.
- `NA` → isotropic displacement.

This is the "Structure-aware label direction" toggle (`sim_ldir_aware`). It matters: an isotropic
linker on a 10 nm-radius ring blurs the ring both inward and outward and fills the lumen; a radial one
inflates the apparent radius instead. Both are real, for different labels.

**Multi-copy** places N independently rotated copies at random positions in a tile
([`sim_server.R:1150-1219`](../sim_server.R#L1150)). Copy 1 is deliberately kept near the origin
(±50 nm) because a real ROI is usually centred on a feature. Per-copy rotation also rotates the
`ldir` vectors and the per-copy structure centre, so the label-direction model survives the copy.

### 1.2 Labeling (stage 2)

`expand_to_fluorophores()` ([`sim_server.R:693`](../sim_server.R#L693)) turns molecules into
**fluorophores**, and *one fluorophore = one future trace ID*. Presets come from
`minflux_label_presets()` in `helpers_shared.R` — a single source of truth shared with the SMLM tab's
inverse GT estimator, so the estimator inverts exactly what this generates.

#### 1.2.1 How many dyes — and how many of them are seen

For molecule *m*, the dye count is a **zero-truncated Poisson**, and each dye is then detected
independently at the labeling efficiency *p* = `eff_pct`/100:

```
K_m  ~ max(1, Poisson(n_fluor))          dyes on the molecule   (K_m = 1 when n_fluor ≤ 1)
K'_m ~ Binomial(K_m, p)                  dyes actually detected
```

`E[max(1, Poisson(λ))] = λ + e^{−λ}`, so the truncation inflates the mean slightly at small `n_fluor`
(λ = 3 → 3.05) and is negligible above ~5.

The consequence people usually get wrong is what *p* does at the **molecule** level. A molecule is
seen if at least one of its dyes is:

```
P(molecule seen) = 1 − (1 − p)^{K_m}
```

Because `K_m` is itself random, the honest quantity is `E[1 − (1−p)^{K_m}]`, **not**
`1 − (1−p)^{n_fluor}` — Jensen's inequality makes the naive form optimistic (at `n_fluor` = 3, 80 %
rather than 88 %). The table below reports the correctly averaged values.

At p = 0.5 a single-dye tag finds half your molecules, but a 9-dye polyclonal finds **98.9 %**. This
is why multi-dye labels are so forgiving of poor efficiency — and simultaneously why they wreck
counting, since the same property multiplies TIDs (below). The simulator gets this for free by
applying efficiency per *dye*, never per molecule.

#### 1.2.2 Where each dye lands

Three displacement levels compose, in order. **Direction** for levels 0 and 1 comes from
`struct_displace()`, which honours the structure's `ldir_type`: `perp` (⊥ the filament tangent),
`radial` (outward from the structure centre), or isotropic.

| Level | Term | Magnitude | Direction | Shared by |
|---|---|---|---|---|
| 0 | Epitope offset | fixed `r = 0.66·MW^(1/3)` nm | structure-aware | all dyes on the molecule |
| 1 | Linker / complex centre | `\|N(0, σ_link)\|`, `σ_link = link_nm·√(π/2)` | structure-aware | all dyes on the molecule |
| 2 | Dye spread in complex | `N(0, dye_spread²)` **per axis** | isotropic | per dye |

The level-1 scale is set so the **mean** magnitude is exactly `link_nm`: a half-normal with scale
`σ_link` has mean `σ_link·√(2/π) = link_nm`.

> ⚠️ **`link_nm` is a mean, not an RMS.** For quadrature you need the second moment:
> `E[|d₁|²] = σ_link² = (π/2)·link_nm²`, so the RMS displacement is **1.253 × `link_nm`**. The
> Precision Budget used to add `link_nm` itself, understating the term by 20 % — corrected, along
> with three other quadrature errors, in §7.10.

#### 1.2.3 Dyes → trace IDs

One detected dye becomes **one or more TIDs**. In Normal mode each dye is split by blinking, `B_i ~
Poisson(λ_blink)` extra appearances, all at the **identical position** (§7.5):

```
T_m = Σ_{i ∈ detected}  (1 + B_i)                 TIDs from molecule m
E[T_m] = n_fluor · p · (1 + λ_blink)
```

In Acquisition mode the blink split is skipped and multiplication comes from the detection rate
instead (§1.3). Either way, the quantity that matters is **TIDs per molecule**, because that is the
over-counting factor the SMLM tab's de-blinking and clustering then have to undo.

#### 1.2.4 The presets, and what they imply

Derived from `minflux_label_presets()` (`helpers_shared.R:105`), all values **measured by simulating
the actual draw** (4×10⁵ molecules) rather than read off the nominal parameters. `E[K]` is the true
mean dye count after the `max(1, ·)` truncation; "TIDs/mol" is `E[K]·(1+λ)` at full efficiency;
"floor" is the irreducible label error `√(π/2)·link_nm`; "min spacing" is the closest two molecules
can sit and still satisfy the existing 2× criterion.

| Preset | n_fluor | E[K] | link nm | spread nm | λ | **TIDs/mol** | seen @ p=0.5 | floor nm | min spacing nm |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Nanobody (DOL≈1) | 1 | 1.00 | 2.5 | 0 | 0.2 | 1.2 | 50 % | 3.1 | 6 |
| Direct tag (SNAP/Halo/Click) | 1 | 1.00 | 3.5 | 0 | 0.2 | 1.2 | 50 % | 4.4 | 9 |
| Primary mAb IgG (DOL≈1) | 1 | 1.00 | 8 | 0 | 0.3 | 1.3 | 50 % | 10.0 | 20 |
| Primary mAb IgG (NHS) | 3 | 3.05 | 8 | 4 | 0.5 | 4.6 | 80 % | 10.0 | 20 |
| mAb + secondary nanobody | 2 | 2.14 | 10.5 | 3 | 0.4 | 3.0 | 70 % | 13.2 | 26 |
| Streptavidin | 3 | 3.05 | 15 | 3 | 0.4 | 4.3 | 80 % | 18.8 | 38 |
| Primary polyclonal IgG (NHS) | 9 | 9.00 | 8 | 9 | 0.5 | 13.5 | 98.9 % | 10.0 | 20 |
| mAb + secondary IgG | 8 | 8.00 | 20 | 12 | 0.6 | 12.8 | 98.2 % | 25.1 | 50 |
| Polyclonal + secondary IgG | 18 | 17.99 | 20 | 18 | 0.8 | **32.4** | ~100 % | 25.1 | 50 |

Note `n_fluor` = 1 takes a separate code branch (`K` = 1 exactly, no Poisson draw), which is why
`E[K]` is 1.00 there and 2.14 for `n_fluor` = 2.

Read the two right-hand columns together with the TID column — they are the whole labeling
trade-off in one place:

- **A polyclonal + secondary sandwich turns one protein into ~32 trace IDs**, none of which can be
  told apart from 32 real proteins without de-blinking or clustering. Counting experiments should not
  use it.
- **The floor is set by the linker alone** and no acquisition time touches it. An IgG sandwich cannot
  resolve anything closer than ~50 nm — so for the 8-molecule filament question, a 15 nm spacing needs
  a nanobody or a direct tag, and *any* antibody route fails before the run even starts.
- Multi-dye labels buy *detection* (99.8 % vs 50 %) and pay in *over-counting* and *spread*. They do
  **not** buy precision: `dye_spread` averages down over dyes, but the shared linker never does.

#### 1.2.5 What survives averaging

Only the terms drawn **per localization** or **per dye** shrink with sampling:

```
ε_m(t) = √( r² + (π/2)·link²  +  dye_spread²/n_dye(t)  +  σ_loc²/n_loc(t) )
          └──── fixed per molecule ────┘   └──── averages away ────┘
```

Levels 0 and 1 are drawn **once per molecule** and shared by every dye and every localization on it,
so they are a bias, not a noise. This is exactly why time-to-resolve has a hard floor (§7.1, §8.2).

A **steric check** compares the label footprint + protein diameter against the mean GT
nearest-neighbour distance and warns when stoichiometric labeling is physically impossible at that
density ([`sim_server.R:1419`](../sim_server.R#L1419)).

### 1.3 TID multiplication (stage 3)

Real MINFLUX over-counts: one physical dye can be re-acquired many times and appears as several
`tid`s. Two mutually exclusive models, selected by the mode:

- **Normal mode — `split_for_blinking()`** ([`sim_server.R:827`](../sim_server.R#L827)). Each
  fluorophore gets `K ~ Poisson(λ)` extra TIDs at the **same position**. Apparent TIDs ≈ (1+λ) × dyes.
- **Acquisition mode — rate × duration sampling** ([`sim_server.R:1298`](../sim_server.R#L1298)).
  Each pool fluorophore gets a photon budget `k_i = 1 + Poisson(blinks−1)`; the requested number of
  detections `n_req` is drawn **without replacement** from that capacity, so dyes with a larger budget
  are found sooner and the total saturates at `Σk_i`.

### 1.4 Precision model (stage 4)

Two modes:

- **Manual** — you type σ_xy and σ_z directly.
- **CRLB from photon budget** ([`sim_server.R:1264`](../sim_server.R#L1264)):

  ```
  N_eff = N_photons × (1 − bg_fraction)
  σ_xy  = L / (2 · √N_eff)
  σ_z   = σ_xy × z_factor            (default 0.7 — z TIGHTER than xy)
  ```

  `L` is the effective targeting-pattern (TCP) diameter of the final iteration. The **σ_z < σ_xy**
  default is correct for MINFLUX and is *not* the camera-SMLM astigmatism reflex — see
  `reference_minflux_axial_precision`.

  **Defaults measured on real data.** Running the calibrator over the fixed-sample test runs
  (ungated, same code path as the app):

  | Run | traces | locs/trace | σ_x | σ_y | σ_z | σ_z/σ_x | photons | L_eff |
  |---|--:|--:|--:|--:|--:|--:|--:|--:|
  | collaborator B2 cond02 | 546 | 8 | 7.06 | 7.16 | 6.17 | 0.87 | 73 | **121** |
  | collaborator B3 cond02 | 472 | 7 | 6.72 | 7.31 | 5.63 | 0.84 | 71 | **118** |
  | JE2 001 (fiducials) | 3679 | 2 | 10.56 | 9.21 | 5.50 | 0.49 | 137 | **231** |
  | JE2 002 (clean) | 4583 | 2 | 11.01 | 10.17 | 6.84 | 0.61 | 105 | **217** |

  The default `L` is **120 nm** (was 50), from the two collaborator runs, which agree to 2 %. At the default
  80 photons that gives σ ≈ 6.9 nm — consistent with the manual-mode defaults (σ_xy = 6, σ_z = 4)
  so the two precision modes no longer disagree by 2×. The JE2 runs sit at ~220 nm: a harder sample
  (median 2 locs/trace, so σ comes from the longer-lived minority), which is the spread a single
  default cannot cover — calibrate rather than trust it.

  **L_eff is not the instrument's L.** Published final-iteration TCP diameters are ~40–75 nm, so
  these effective values run 2–3× larger. That is expected: `L_eff = 2σ√N` absorbs everything that
  inflates σ beyond the ideal CRB — drift, molecular motion within the trace, background, imperfect
  centring. It is the right quantity to simulate with, and the wrong one to quote as an instrument
  setting.

  The measured σ_z/σ_x of 0.49–0.87 brackets the `sim_crlb_z_factor` default of 0.7, and every run
  is below 1 — independent confirmation that MINFLUX z is *tighter* than xy.

  The `2` in the prefactor is a convention that varies in the literature with pattern geometry. It is
  **self-consistent here** because "Calibrate from a real run" inverts the same formula
  (`L ≈ 2·σ·√N`, [`sim_ui.R:532`](../sim_ui.R#L532)). Do not "fix" it in isolation — that silently
  rescales every previously calibrated L.

### 1.5 Localizations (stage 5)

`gen_localizations()` ([`sim_server.R:848`](../sim_server.R#L848)) draws a trace length per
fluorophore and scatters that many points around it:

- **Poisson(λ)** by default; **Geometric(1/λ)** with the "heavy-tailed trace length" toggle — same
  mean, heavy right tail, which is what constant per-step bleaching actually produces and what real
  MINFLUX trace-length histograms look like.
- **Scale by fluorophore count** divides λ by the number of dyes on the molecule, modelling a fixed
  observation budget per molecule rather than per dye.
- Image mode overrides λ per molecule via `n_locs_expected` (brighter pixel → longer trace).

### 1.6 Background (stage 6)

Uniform random positions over the FOV (or the multi-copy tile, or the GT bounding box + 10 %), each
generating its own trace. Count is a fixed number in Normal mode, and `bg_rate × duration` in
Acquisition mode so longer runs accumulate more non-specific TIDs. Background TIDs are numbered above
all real TIDs and carry `mol_id = NA`, so scoring excludes them.

### 1.7 Scoring (stage 7)

- **RMSE** — per-TID distance between the trace median and its source GT molecule (real TIDs only,
  capped at 1000 sampled TIDs).
- **Precision budget** ([`sim_server.R:2594`](../sim_server.R#L2594)) — a bar chart of the error terms
  as **radial RMS**, their quadrature sum, the observed RMS overlaid, the **label floor**, and a
  **resolvability verdict** comparing GT nearest-neighbour spacing against 2× the per-molecule error.
  The formula and its validation are in §7.10.

---

## 2. Normal vs Acquisition mode

| | Normal | Acquisition |
|---|---|---|
| Question | "what does this structure look like with this label?" | "how long do I need to run?" |
| Over-counting | blink λ | detection rate × duration, capped by photon budget |
| FOV | optional pinned ROI | always on |
| Timestamps | synthetic index (`0, 10, 20 …`) | sampled from the detection-rate profile |
| Extra output | — | Acquisition timeline tab |

## 3. The acquisition model

All three rate profiles are one family, `r(t) = struct_rate · g(t)`, with `g` normalised to the
**characteristic** rate ([`sim_server.R:111`](../sim_server.R#L111)):

| Model | `g(t)` | Characteristic rate | Physical picture |
|---|---|---|---|
| `decay` | `level + (1−level)·e^{−t/τ}` | initial **peak** | dyes bright early, bleach to a floor |
| `constant` | `1` | the flat rate | variation is Poisson only |
| `rise` | `1 − (1−level)·e^{−t/τ}` | the **plateau** | activation-limited (buffer, photoactivation) |

`compute_acq_timeline()` ([`sim_server.R:158`](../sim_server.R#L158)) then integrates deterministically
(no per-realization jitter, so the curve is stable while you drag sliders):

```
cum_R(t)   = rate · ∫₀ᵗ g            cumulative detection attempts
bleach(t)  = exp(−cum_R / struct_cap)   fraction of the photon budget still emitting
N(t)       = struct_cap · (1 − bleach)  total recorded TIDs (over-counting included)
unique(t)  = Σ_m [ 1 − exp(−f_m · c_m · cum_R / n_pool) ]   distinct molecules found
```

`f_m` = dyes on molecule *m*; `c_m` = the optional **crowding** weight `1/(1+neighbours within R)`,
which slows discovery where emitters sit inside one doughnut. **Time-to-resolve is driven by
`cum_R` only** — it is deliberately independent of the per-dye blink budget, so "detections before
bleaching" changes over-counting but never the road to X %.

This timeline is recomputed **live** (`sim_acq_live()`) from the last-generated fluorophore pool, so
rate/duration/background respond without re-Generating. See §7.5 for the consequence.

## 4. Calibrate from a real run

`analyze_real_run()` ([`sim_server.R:244`](../sim_server.R#L244)) accepts a real `.npy` / `.json` /
`.csv` and measures, from first-detection time per `tid`, everything below.

**Reading conventions.** The loader mirrors the SMLM tab: `.npy` / `.json` are read with
`z_scale = 1.0` (raw, uncorrected z), a `.csv` is returned to raw by `strip_csv_z_calibration()`,
and every source then goes through `canonicalize_minflux_df()` — so alternative headers
(`trace_id`, `x_nm`, `photons`, …) are accepted here as they are elsewhere. **One** axial factor is
then applied, under user control:

| Source | Detected | Input pre-filled with |
|---|---|---|
| `.npy` / `.json` | raw instrument data | 0.7 (the correction still to apply) |
| `.csv` from this app | factor recorded in the file (`z_scale_applied`) | that exact factor, marked *certain* |
| `.csv` from pyMINFLUX / unknown | assumed | 0.7, marked *assumed* |
| any 2D run | no axial structure | nothing — the control is hidden |

The factor is applied in one place (`sim_real_gated()`), so σ_z and the FOV depth can never disagree
with the number on screen, and changing it re-runs the calibration. Reading with the readers' own
0.7 default instead would bake the correction into `.npy`/`.json` while leaving a `.csv` at whatever
its producer applied — the same run would then yield a different σ_z per format.

> **Two different 0.7s.** This axial factor is a refractive-index correction on the *coordinates*.
> `sim_crlb_z_factor` (also ≈0.7) is the σ_z/σ_xy *precision ratio*. They are unrelated and must not
> be conflated — see `reference_minflux_axial_precision`.

**The calibration gate** (`gate_real_run()`, [`sim_server.R:1559`](../sim_server.R#L1559)) restricts
*which* localizations and traces those statistics come from, in two stages:

1. **EFO** and **CFR** ranges gate individual localizations.
2. Trace length is then **recounted on the survivors**, and the length range drops whole traces.

So "localizations per trace" means *usable* locs per trace — gating hard on EFO and then asking for
≥ 20 locs is far stricter than either gate alone, which is why the info box reports both stages.
Ranges default to each file's own full extent (a true no-op); absent columns disable their gate
rather than silently passing everything; non-finite values are kept and counted separately, so a
full-range gate never drops anything. Changing any range re-runs the whole calibration.

Gating is **not neutral**: EFO, CFR and length all select the brighter, better-centred, longer-lived
molecules, which are also the most precise ones — so a narrowed gate yields a σ tighter than the run
as a whole. The info box says so whenever a range is actually narrowed.

Measured, after gating:

| Measured | How | Fills |
|---|---|---|
| Duration | `max(t0) − min(t0)` | hours / minutes |
| Rate profile | binned traces/hour, log-linear fit of `|rate − late|` → τ | τ, and the model choice |
| Model | exponential vs flat residual comparison; sign of (early − late) picks decay/rise | `sim_acq_det_model`, `level` |
| FOV | coordinate ranges × 1.1 | `sim_fov_x/y/z` |
| Structure fraction | DBSCAN core-point test (≥3 neighbours within 60 nm) on trace medians | `sim_real_struct_frac`, overridable |
| σ_loc, per axis | median over traces of the within-trace SD (≥3 locs), **small-n debiased**, reported ± nMAD | manual σ_xy / σ_z |
| Axis ratios | median of the **per-trace** `σ_y/σ_x`, `σ_z/σ_x` | cloud anisotropy (8.1c) |
| Photon distribution | the whole `eco` vector (capped at 200 k) | per-localization σ (8.1c) |
| Trace-length distribution | localizations per `tid`, the whole vector | trace lengths (8.1c) |
| Photons | median `eco` | CRLB photons |
| L | `L_eff = 2·σ·√N` | CRLB L |

**Small-n debiasing.** The sample SD under-estimates σ at small *n*, and the length gate lets the
user push short traces into the median, so each trace's SD is divided by its own correction factor.
The textbook constant `c₄(n) = E[s]/σ` is **the wrong one here**: it un-biases the *mean* of the
sampling distribution, while this aggregates with a *median*, and `median(s) < E[s]` at small *n*.
The matching constant is the median of the same chi distribution,
`m(n) = √(qchisq(0.5, n−1)/(n−1))` — at *n* = 3, `c₄ = 0.886` but `m = 0.833`. Measured on a
short-trace-heavy population (truth σ = 10): raw median **−7.6 %**, c₄-corrected **−2.5 %**,
m-corrected **−0.05 %** (200 realizations, MC se 0.05 %). Either constant **cancels exactly** in the
axis ratios, since both axes share the same *n*.

The ratios are taken as the **median of per-trace ratios**, not the ratio of medians: within one
trace both axes share a photon count, so the per-trace ratio divides the photon dependence out and
measures cloud shape alone. `median(log(s_y/s_x))` is exactly `log(σ_y/σ_x)` because the two noise
terms are i.i.d. and their difference is symmetric — verified unbiased to +0.11 % (MC se 0.08 %).

`apply_real_rates()` splits the measured total rate into a **constant background offset** and a
**time-varying structure component**, using the structure fraction, and pushes both plus the level.

**What it does *not* transfer:** the structure geometry, the point-cloud shape, the trace-length
distribution, or the labeling. It calibrates the *acquisition clock and the precision*, and you still
have to pick a synthetic shape by hand. That gap is the subject of §8.

## 5. The Confocal comparison view

Genuinely a *detection* PSF model, unlike the MINFLUX view — the word "PSF" is correct here:

```
σ_ex  = 0.21·λ_ex/NA ;  σ_em = 0.21·λ_em/NA
1 AU  = 1.22·λ_em/NA ;  r_ph = AU·(1 AU)/2
σ_det = √(σ_em² + (r_ph/2.355)²)
σ_conf = σ_ex·σ_det / √(σ_ex² + σ_det²)
```

with Airy encircled-energy collection efficiency `1 − J₀(v)² − J₁(v)²`, `v = 3.8317·AU`, which
reproduces the real pinhole trade-off (closing it sharpens and sections but costs photons). Rendered
from the **unique fluorophore pool**, not the over-counted detection list, so a re-detected molecule
is not drawn brighter than it physically is. Poisson shot noise on top.

## 6. Phase 2 — Diffusion / Tracking, in one paragraph

Brownian motion on a sphere / rod / ellipsoid / membrane plane / unbounded volume, surface (2 DOF) or
volume (3 DOF), with free / OU-confined / directed / two-state motion, sampled through a MINFLUX
overlay: cycle clock, geometric dark cycles → irregular Δt, **per-localization** `N ~ Poisson(photons)`
→ `σ = L/(2√N)`, photobleaching-limited length, optional re-acquisition survivorship. Emits
`tid, tim, x, y, z, efo, eco` plus a ground-truth sidecar, and scores itself with a built-in CVE
round-trip. Full spec: [simulation_module.md](simulation_module.md).

---

## 7. Known fidelity gaps and inconsistencies

Ranked by how much they distort the tab's headline question. Items 1, 5, 6, 7 were found in this
review; 2, 3, 4 come from the 2026-08-11 audit
(from the 2026-08-11 audit). *Status 2026-09-17:* 7.10 is fixed. 7.1 is partly
addressed: the Precision Budget now headlines the label floor, the per-molecule error ε and a
resolvable / not-resolvable verdict (§8.2), but the acquisition timeline's "% resolved" and
time-to-resolve still count first detection. The others are open.

### 7.1 "% resolved" measures *discovery*, not *resolution* — the biggest one

The Acquisition tab reports "time to 90 % resolved", and a molecule counts as resolved the moment
**any one of its dyes is detected once** ([`sim_server.R:208`](../sim_server.R#L208)). Meanwhile the
Precision Budget calls the structure "resolvable" from a purely static criterion — GT spacing vs
`2·√(r² + link² + σ²)` — with no time in it at all. Two panels, two different meanings of the same
word, and neither is the quantity the user actually wants.

The physical quantity is time-dependent and neither panel computes it. Per molecule, the centroid
error after `n` localizations is

```
ε(n) = √( r_protein² + link²  +  σ_loc²/n )
        └── fixed label bias ──┘   └ shrinks with sampling ┘
```

So more time buys you `σ_loc/√n` and **nothing else** — the label bias is a hard floor. A structure
with 12 nm spacing and an IgG sandwich (link ≈ 20 nm) never resolves, no matter how long you run,
and the tab will still cheerfully report "90 % resolved in 8 minutes". This is the single most
misleading output in the simulator.

**Fix direction:** make resolution a time-dependent, per-molecule quantity — accumulate localizations
per molecule over the timeline, compute ε(n(t)), and report time-to-`ε(n) < spacing/2` alongside (not
instead of) the discovery curve. Report the floor `√(r²+link²)` explicitly and say "never" when it
exceeds the criterion.

### 7.2 σ_loc is constant across every localization

`gen_localizations()` applies one scalar σ to every point. Real MINFLUX σ is set **per localization**
by that localization's photon count. Consequences:

- Simulated precision clouds are homogeneous Gaussians. The app's own precision-cloud QC predicts
  **leptokurtosis from the σ spread** rather than testing against 3 — so simulated data can never
  exercise that check, and anything tuned against it is tuned against the wrong null.
- Every inverse-variance weighting downstream becomes a no-op on simulated data.

The correct implementation already exists ~1800 lines away, in the Phase 2 overlay
([`sim_server.R:3097`](../sim_server.R#L3097)): draw `N ~ Poisson(photons)` per localization, set
`σ = L/(2√N)`, add noise with that σ. Port it.

### 7.3 The CSV export drops the MINFLUX-native columns

`sim_dl_csv` ([`sim_server.R:2673`](../sim_server.R#L2673)) writes only `x, y, z, tid, tim` + GT
columns. No `eco` / `efo` / `cfr`. Loading simulated SMLM data back into the SMLM tab therefore:

| Feature | Default | Behaviour without `eco`/`efo` |
|---|---|---|
| Shape-fit photon weighting (`viewer_weight_eco`) | **ON** | silently falls back to unweighted |
| Fiducial / contaminant scan | — | its primary BRIGHT signature is `efo` — cannot run |
| efo / cfr metadata filters | — | inert |

Fixing 7.2 produces `N` per localization, so `eco = N` and `efo = N/dwell` fall out for free. `cfr`
has no counterpart in the model; **recommended: omit it and state that simulated data cannot exercise
cfr-based tools.** A constant placeholder is the worst option — it reads as real.

### 7.4 Background is modelled as a plain photon loss

CRLB mode uses `N_eff = N × (1 − bg)`. This is the one place a camera-SMLM intuition genuinely does
not transfer: MINFLUX position information lives in the **contrast at the doughnut zero**, and
background **fills that zero in**, so the CRB degrades faster than the photon loss alone implies.
The current model is optimistic, increasingly so at high background. Impact is confined to CRLB mode
used *predictively* — if L came from "Calibrate from a real run", the error is largely absorbed into
the calibrated L.

### 7.5 Blink copies sit at *exactly* the same position

Both over-counting models (`split_for_blinking`, and acquisition re-sampling of the pool) duplicate a
fluorophore row **without re-drawing its position**. Real re-detections of the same dye are separated
by the between-blink motion of the molecule plus any drift, which is precisely the quantity the SMLM
tab's de-blink radius has to be set against — and which we already know `k_clear` is blind to
(`reference_deblink_k_motion_blind`). Simulated data has zero between-blink displacement, so it
**cannot be used to calibrate or validate the de-blink threshold**. Anything that looks like a clean
de-blink result on simulated data is an artifact of this simplification.

### 7.6 Timestamps in Normal mode are not physical

`tim = seq(0, by = 10, length.out = total)` ([`sim_server.R:894`](../sim_server.R#L894)) — a single
global ramp across *all* localizations of *all* traces. So every trace occupies a disjoint,
non-overlapping time window, and the within-trace Δt is 10 (units unstated) instead of the ~ms cycle
time. Any time-domain tool (TSA, cycle time, detection-over-time, dark-cycle colouring) reads garbage
from a Normal-mode export. Acquisition mode overwrites this with detection-profile times + 0.01 s
within-trace spacing, which is sane but still a fixed cadence with no dark cycles.

### 7.7 Crowding affects the timeline but not the point cloud

`sim_crowding` scales the per-molecule *discovery rate* in `sim_acq_live()`, but the generated
localizations are untouched — crowded molecules get the same σ and the same trace lengths as isolated
ones. In reality, a second emitter inside the doughnut degrades the localization itself (biased
position, elevated `cfr`). The Acquisition summary and the MINFLUX view therefore tell different
stories about the same crowding.

### 7.8 Live timeline vs frozen point cloud can disagree

`sim_result()` is an `eventReactive` on **Generate**; `sim_acq_live()` recomputes on every parameter
change. Change the duration or rate and the timeline updates immediately while the displayed
localizations, RMSE, and precision budget still reflect the parameters at the last Generate. This is
intentional (responsiveness) but undocumented in the UI.

### 7.9 Documentation drift

The main [README](../README.md) claims "Simulated data can be sent directly to the SMLM or Tracking
tabs for immediate analysis." There is no in-memory hand-off — the only route is
*Download localizations (.csv)* and re-upload. (Corrected in the README as part of this review.)

### 7.10 Precision Budget quadrature — **FIXED 2026-08-11**

Found while documenting the labeling maths (§1.2.2). The panel used to compute

```r
theory <- sqrt(prot_r^2 + link^2 + sig_xy^2)      # every term wrong in a different way
```

Four errors, all optimistic:

1. **`link` is a mean, not an RMS.** The level-1 displacement is half-normal with
   `E[|d|²] = (π/2)·link²`, so the term is `√(π/2)·link` — using `link` understated it by **20 %**.
2. **Units mixed.** `prot_r`/`link` are *radial* magnitudes, `sig_xy` a *per-axis* σ. The radial
   variance of one localization is `2σ_xy² + σ_z²` (3D) or `2σ_xy²` (2D).
3. **`dye_spread` omitted entirely** — for a polyclonal + secondary that is an 18 nm term.
4. **No trace averaging.** A TID's position is the *median of its localizations*, so the
   localization term must be divided down by the trace length. Without it, "observed vs theoretical"
   compared two different quantities — and the panel's own note ("observed should approach
   theoretical when traces are long") had it backwards.

Now:

```
ε_TID = √( r² + (π/2)·link² + D·dye_spread² + (2σ_xy² [+ σ_z²])·f(n) )
f(n) = (π/2) / (n + π/2 − 1)          variance of the median of n Gaussian draws
```

`f` is the finite-*n* form: the textbook asymptotic `π/(2n)` runs **+17 % at n = 3** and +9 % at
n = 5, exactly the short traces MINFLUX produces. It is evaluated as the **mean over the realized
trace lengths**, not at the median length — `f` is convex and Poisson lengths otherwise cost another
5–7 % by Jensen.

**Validated against the generator itself**: the formula is compared with the observed radial RMS
`|TID median − GT|` produced by the real `expand_to_fluorophores` → `gen_localizations` chain, across
12 regimes (2D/3D, label- and localization-dominated, trace lengths 3–60, all preset shapes). Worst
error **2.4 %**; the old formula was off by −42 % to +76 % over the same set.

The panel now also reports the **label floor** `√(r² + (π/2)·link²)` — the part no acquisition time
removes — and the resolvability verdict distinguishes "⚠ below the limit *for now*" from
"✖ not resolvable with this label", which is a different problem with a different fix.

### What was checked and is **correct** — do not re-litigate

- The `L/(2√N)` prefactor convention (self-consistent with the calibration inverse).
- `σ_z = 0.7·σ_xy` — z tighter than xy is right for 3D MINFLUX.
- The confocal PSF model (`0.21·λ/NA`, Airy pinhole, Bessel encircled energy).
- Variance ∝ 1/photons — the 1/√N scaling is shared with camera SMLM; only the prefactor differs
  (targeting-pattern size instead of spot width), so photon weighting is sound.
- Phase 2's sampling overlay in full.

### How to verify a change

- Source helpers into an env rather than launching the app:
  `e <- new.env(); sys.source("helpers_smlm.R", envir = e)` (`sys.source("global.R")` first for
  constants; load `shinyjs` for UI checks).
- Parse check every touched file: `Rscript -e "invisible(parse('sim_server.R'))"`.
- **Acceptance test for 7.2:** simulate → export → load into the SMLM tab → run the precision cloud.
  Before the fix the pooled cloud is a clean Gaussian; after it the σ spread should make it
  leptokurtic and the panel's predicted-kurtosis check should sit near its expectation instead of
  firing spuriously.
- **Round-trip for 7.3:** export, confirm `eco`/`efo` present, reload, confirm `viewer_weight_eco`
  reports a non-trivial `n_eff` instead of falling back.

---

## 8. PLANNED — real-run calibration + honest time-to-resolve

*Agreed scope, 2026-08-11. Two work items only.* **Build status (checked against `sim_server.R`,
2026-09-17):** 8.1a gate **built** (`sim_real_gate_ui`: EFO / CFR / trace-length ranges with before and
after counts). 8.1b profile **built** (small-n debiased per-axis σ, per-trace axis ratios, `eco_vec` and
`len_vec` stored). 8.1c **not built**: the stored vectors are displayed but `gen_localizations` does not
resample them, so the simulated cloud still uses one σ and Poisson/geometric trace lengths, and the
bias-3 guard in 8.1d has nothing to guard yet. 8.2 **partly built**: the label floor, ε per molecule and
the 2ε criterion are in the Precision Budget; ε_m(t) and a time-to-resolve based on it are not, so the
timeline still reports discovery (§7.1). Explicitly
**deferred** (considered and dropped as over-engineering for now): empirical cluster-geometry fitting,
the simulated-vs-real round-trip validation, blink-displacement modelling, physical Normal-mode
timestamps, the non-linear background CRB, and `cfr` emission. §7 keeps them on record.

### 8.0 The workflow this has to serve

> *"I think this structure may be an 8-molecule filament, labeled with Y. Given a **real** MINFLUX run,
> how long do I need to image before I resolve it?"*

The structure stays **hypothesised and hand-picked** — that is the point of the question, not a
limitation. What the real run supplies is the *instrument and sample behaviour*: how fast traces
arrive, how many localizations each trace yields, and how tight the resulting point cloud is.

```
  real run  ──gate──►  point-cloud profile  ──┐
  (efo, cfr, trace length)                     │
                                               ├──►  simulate  ──►  ε_m(t)  ──►  time to resolve
  hypothesised structure (8-mol filament)  ────┤                                 (or "never")
  labeling preset Y                        ────┘
```

---

### 8.1 Item 1 — gated point-cloud calibration

**Goal:** replace the single scalar σ with a point-cloud model measured from a *user-selected subset*
of a real run.

#### 8.1a New gating controls

Applied **before** any statistic is computed, in the "Calibrate from a real run" block:

| Control | Acts on | Notes |
|---|---|---|
| `efo` range (min / max, kHz) | localizations | drop dim and abnormally bright locs |
| `cfr` range (min / max) | localizations | drop poorly-centred locs |
| Trace length range (min / max locs per `tid`) | whole traces | drop stubs and fiducial-like runaways |

Columns may be absent (a plain CSV often has no `efo`/`cfr`): each control disables itself with a
visible note rather than silently passing everything. Report **n before / n after** for both
localizations and traces, so the gate's severity is never invisible.

#### 8.1b What is measured from the surviving traces

| Quantity | How | Used for |
|---|---|---|
| σ_x, σ_y, σ_z | within-trace SD per `tid`, **small-n debiased**, median over traces | per-axis cloud shape |
| Axis ratios | `σ_y/σ_x`, `σ_z/σ_x` | anisotropy of the simulated cloud |
| Photon distribution | the `eco` **vector** (not just its median) | per-localization σ |
| Trace-length distribution | locs per `tid`, the **vector** | how many points a trace gets |
| `L_eff` | `2·σ̄_xy·√N̄` (existing) | ties the two together |

The small-n correction matters once the user controls the length gate — and the constant to use is
**not** c₄; see §4, where the choice is derived and measured. Uncorrected, a short-trace-heavy gate
biases the calibrated σ low by ~7.6 %.

#### 8.1c How the simulator consumes the profile

Two new options, each a *third* choice alongside what exists today — the current behaviour stays the
default so nothing silently changes:

1. **Trace length → "Empirical (from real run)"**, alongside Poisson and Geometric. Draw each
   simulated trace's localization count by **resampling the measured vector** with replacement. No
   distributional assumption, which is the honest choice given we already know real trace lengths are
   neither Poisson nor cleanly geometric.

2. **Precision → per-localization, anisotropic.** Per localization *i*, resample `N_i` from the
   measured `eco` vector and set

   ```
   σ_x,i = L_eff / (2·√N_i)
   σ_y,i = σ_x,i · (σ_y/σ_x)          ← measured ratio
   σ_z,i = σ_x,i · (σ_z/σ_x)          ← measured ratio, NOT the 0.7 default
   ```

   This is 7.2 and 7.3 delivered together: `N_i` is exactly what `eco` should be, and
   `efo_i = N_i / dwell` follows, so the CSV export gains both columns for free. The construction is
   self-consistent — because `L_eff` was derived by inverting the same formula, the median simulated
   σ reproduces the measured σ by construction, while the *spread* around it now comes from the real
   photon distribution.

   `cfr` stays **omitted**, per §7.3.

#### 8.1d Three biases the UI must state, not hide

These are properties of gating itself, not bugs, but they all push the same direction — **optimistic**
— and they compound:

1. **Gating on `efo`/`eco` selects bright localizations**, which are the precise ones. The calibrated
   σ is therefore tighter than the run as a whole.
2. **Gating on `cfr` selects well-centred molecules**, again the precise ones.
3. **Gating on trace length is circular for the length model.** If you filter traces to ≥ 20 locs and
   then resample that filtered distribution to predict how many points a trace yields, you have
   hard-coded your own filter into the prediction and will over-estimate the photons available per
   molecule. This is the one that most directly corrupts the time-to-resolve answer.

Bias 3 deserves a hard guard rather than prose: when the length gate is active, either resample the
length distribution from the **ungated** traces (recommended default — the gate then only cleans the
σ estimate, which is what it is actually for), or require an explicit "use the gated length
distribution" tick with the consequence spelled out next to it.

There is also a physical conflation worth stating once in the doc and once in the UI: **within-trace
SD is not pure localization precision.** It contains any real molecular motion, stage drift and
sample wobble during the trace. On a static SMLM target that is small; on anything mobile it is not,
and the calibrated σ is then an upper bound on precision rather than an estimate of it.

---

### 8.2 Item 2 — resolution as a function of time

**Goal:** replace "% resolved = % detected at least once" (§7.1) with the real quantity.

#### 8.2a The error model, derived from the labeling model

The four displacement levels of §1.2 average differently, and that is the whole result. For molecule
*m* with `n_dye` of its dyes detected and `n_loc` localizations accumulated in total:

```
                 ┌───── never averages ─────┐  ┌ averages over dyes ┐ ┌ averages over locs ┐
ε_m(t) = √(  r_protein²  +  (π/2)·link²      +  D·dye_spread²/n_dye(t) + D·σ̄²/n_loc(t)     )
```

with `D` = 2 or 3 (the observed dimensionality). The `π/2` and the `D` are the two corrections §7.10
records as missing from the current Precision Budget — `link_nm` is a *mean* magnitude while
quadrature needs the RMS, and `dye_spread`/`σ̄` are *per-axis* while the others are radial.

`r_protein` (epitope offset) and `link` (complex centre) are drawn **once per molecule** and shared by
every dye and every localization on it — so no amount of acquisition time removes them.

**The floor is `√(r_protein² + link²)`.** This single number decides whether the question has an
answer at all, and it is a property of the *label*, not the run.

#### 8.2b The criterion

A molecule is resolved from its neighbour when `ε_m(t) < d_NN(m) / k`, with `d_NN` the ground-truth
nearest-neighbour distance and `k` a user-set separation criterion (default 2, matching the existing
Precision Budget). Flag in the UI that `k = 2` compares a spacing against a *single* cloud's error;
two independent clouds each of error ε separate at `d ≳ 2√2·ε`, so `k = 2√2 ≈ 2.83` is the stricter,
arguably more correct reading. Offering both and defaulting to the existing convention keeps old
numbers comparable while making the choice explicit.

#### 8.2c What gets reported

Replace the headline number; **keep the discovery curve as a second series** — it is still the right
answer to "how much of the structure have I even seen", and `n_loc = 0 ⟹ ε = ∞` means the new
criterion strictly subsumes it.

| Output | Why |
|---|---|
| **Label floor** `√(r² + link²)` vs the criterion `d_NN/k` | the go / no-go, shown first |
| **"Never resolvable" fraction** — molecules whose floor already exceeds the criterion | with the verdict *"change the label, not the run time"* |
| Time-to-X %-resolved curve, for the molecules that can resolve | the actual answer |
| Discovery curve (current "% resolved"), relabelled **"% detected"** | honest naming of what it always was |
| ε(t) for the median molecule | shows how fast the σ/√n term is still buying anything |

The relabelling is not cosmetic — with both curves on one plot, the gap between "detected" and
"resolved" *is* the thing the tool exists to show.

#### 8.2d One modelling assumption to state

`dye_spread²/n_dye` assumes the analysis **merges all TIDs belonging to one molecule** into a single
position — which is what clustering and de-blinking actually do downstream. If the dyes are *not*
merged, each TID is its own apparent molecule and a polyclonal label makes the structure look denser,
not more precise. The merged case is the right default because it matches the pipeline, but it is an
assumption about the *analysis*, not about the physics, and the summary should say so.

---

### 8.3 Build order

1. ~~**8.1a** gating controls + before/after counts (self-contained, no model change).~~
   **DONE 2026-08-11** — `gate_real_run()` + `sim_real_gate_ui` / `sim_real_gate_info`; calibration
   now re-runs from one path on both upload and gate change (see §4).
1b. **DONE 2026-08-11** — read path aligned with the SMLM tab (raw z + canonicalize) and a
   user-controlled axial factor pre-filled from the file's own provenance. Prerequisite for the
   σ_z/σ_x axis ratio in 8.1b, which would otherwise be biased by 0.7 and differ per file format.
2. ~~**8.1b** profile measurement incl. the c₄ correction and the axis ratios.~~
   **DONE 2026-08-11** — per-axis debiased σ ± nMAD, per-trace axis ratios, and the whole `eco` /
   trace-length vectors kept for 8.1c. The c₄ constant in the original plan was wrong for a median
   aggregator and was replaced (§4). Recovery verified against known σ and anisotropy.
3. **8.1c** per-localization anisotropic σ + empirical trace lengths + `eco`/`efo` in the CSV export.
   *This is where §7.2 and §7.3 close.*
4. **8.2** the ε(t) model, the floor, the "never resolvable" verdict, and the relabelled curves.

Steps 1–3 are additive and low-risk. Step 4 changes a headline number the tab already reports, so it
lands last and keeps the old curve visible beside the new one.

### 8.4 Open question

**Whose spacing sets `d_NN`** — the ground-truth molecules (what you want to resolve) or the detected
clusters (what you can actually see)? They differ by exactly the labeling efficiency: at 50 %
efficiency the detected spacing is ~1.4× the true spacing, so using it reports a *shorter* time to
resolve a structure you are only half seeing. Ground-truth spacing is the honest default for the
"can I resolve my hypothesised structure" question; the detected spacing answers "can I resolve what
I got", which is a different question worth exposing separately.
