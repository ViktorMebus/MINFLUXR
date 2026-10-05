# MINFLUX Simulator — Reference

*Reference for the Simulator tab: what it models, how each stage works, and where it departs from
real MINFLUX (§7).* Code: [`sim_server.R`](../sim_server.R), [`sim_ui.R`](../sim_ui.R) ·
Companions: [simulation_module.md](simulation_module.md) (the diffusion simulator),
[MINFLUX_analysis_methods.md](MINFLUX_analysis_methods.md) (the analysis side).

---

## 0. What the simulator is for

MINFLUX localizes with an **excitation doughnut**: the molecule is found by minimising fluorescence at
the doughnut zero, and each *trace* (`tid`) is one molecule held and re-localized many times. A trace
therefore comes out as a **point cloud**, and the cloud width sets the achievable resolution.

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
  (4) precision model                 manual σ  |  CRLB σ = √[(k·L/√(8N)·g)² + σ₀²]
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

Because `K_m` is itself random, the detection probability is `E[1 − (1−p)^{K_m}]`. The plug-in form
`1 − (1−p)^{n_fluor}` is optimistic by Jensen's inequality: at `n_fluor` = 3 it gives 88 %
where the true value is 80 %. The table below reports the averaged values.

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
> `E[|d₁|²] = σ_link² = (π/2)·link_nm²`, so the RMS displacement is **1.253 × `link_nm`**, and the
> Precision Budget uses that second moment (§1.4).

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

Derived from `minflux_label_presets()` (`helpers_shared.R:105`), all values **computed by simulating
the actual draw** (4×10⁵ molecules). `E[K]` is the true
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
so averaging over localizations cannot remove them, and time-to-resolve has a hard floor (§7).

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
- **CRLB from photon budget** (`smlm_sigma_loc()` in [`packages.R`](../packages.R)):

  ```
  σ_CRLB = L / (√8 · √N) · g,   g = √(4 − f) / (2 · (1 − f))
  σ_xy   = √( (k · σ_CRLB)² + σ₀² )
  σ_z    = σ_xy × z_factor            (default 0.7 — z TIGHTER than xy)
  ```

  `L` is the physical diameter of the final targeting pattern (TCP), the instrument setting (imaging
  sequences end near 40 nm). `N` is the total detected photons per localization and `f` the background
  fraction of them. `σ_z` defaults to 0.7 × `σ_xy`: in MINFLUX the axial precision is usually tighter
  than the lateral one.

  **The bound.** `σ_CRLB` is the ideal per-axis Cramér-Rao bound at the centre of a 2D doughnut
  pattern (three exposures on a circle of diameter L plus one at the centre), in the small-L limit.
  It is Balzarotti et al. 2017, Suppl. Eq. S27/S31, re-derived as Eq. 44/45 in Marin & Ries 2025
  (arXiv:2410.12427). With SBR = (1 − f)/f, their background penalty
  `√((1 + 1/SBR)(1 + 3/(4·SBR)))` reduces to `g` above; g = 1 at f = 0 and 1.05 at f = 0.05. The
  small-L limit drops a beam-size factor `8σq²/(8σq² − L²)` that is about 1.03 at L = 40 nm.

  **Two excess terms, because they scale differently with photons.**
  - `k` multiplies the bound. It absorbs everything that degrades the photon-limited precision:
    the hexagonal TCP of commercial systems (the bound above is for the triangular pattern), the
    iterative estimator, detector effects. k = 1 is the ideal; k < 1 is possible only if the real
    pattern beats the reference geometry.
  - `σ₀` adds in quadrature. It absorbs everything that does not shrink with photons: residual drift,
    vibration, motion of the dye within a trace. More photons cannot beat it, so a photon sweep with
    σ₀ > 0 levels off where real data does.

  A single multiplier cannot do both jobs. Folding drift into L (or k) makes the drift term shrink as
  1/√N and over-promises what brighter dyes or longer dwell buy. Both terms have literature support:
  camera SMLM finds a photon-independent multiplicative excess over the CRLB (Brenlla et al. 2025,
  Biophys. Rep.), and MINFLUX reviews name setup instability as the leading error at the nanometre
  scale (Scheiderer, Marin & Ries 2025, Nat. Photonics).

  **Calibration.** *Calibrate from a real run* fits both terms. Since σ² = A/N + B, the within-trace
  variance is linear in 1/N: traces are binned by their mean 1/eco (3–8 quantile bins, ≥ ~20 traces
  each), each bin's variance is pooled (Σ(n−1)s² / Σ(n−1), unbiased at any n), and a weighted line
  (weights dof/σ⁴) gives A (slope) and B = σ₀² (intercept). k = √A / (L·g/√8) at the L and background
  set when calibrating, so change L first if yours is not 40 nm. A negative term is set to zero and the
  other refitted. If 1/N spans less than 1.5× across bins, slope and intercept trade off freely: the
  fit is refused, k is set to 1 and the rest of the measured variance goes to σ₀. On synthetic traces
  (800 traces, lognormal photons) the fit recovers k within ~3% and σ₀ within ~0.3 nm.

  The fit assumes the excess variance does not itself depend on brightness. Bright and dim traces
  can differ in other ways (focal position, dye state), which the fit would read as a photon effect.
  The within-trace variance also contains real motion, so on a mobile target σ₀ is an upper bound.

**Precision Budget.** The panel predicts the error of each molecule's position, a trace median,
against ground truth:

```
ε_TID = √( r² + (π/2)·link² + D·dye_spread² + (2σ_xy² [+ σ_z²])·f(n) )
f(n)  = (π/2) / (n + π/2 − 1)        variance of the median of n Gaussian draws
```

`link` is the mean of a half-normal displacement, hence the `π/2`. `f(n)` is the finite-n form of the
median's variance, averaged over the realized trace lengths. The panel reports the **label floor**
`√(r² + (π/2)·link²)` separately, since no acquisition time removes it, and distinguishes "below the
limit for now" from "not resolvable with this label".

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
> be conflated.

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
| Rate profile | binned traces/hour, log-linear fit of `\|rate − late\|` → τ | τ, and the model choice |
| Model | exponential vs flat residual comparison; sign of (early − late) picks decay/rise | `sim_acq_det_model`, `level` |
| FOV | coordinate ranges × 1.1 | `sim_fov_x/y/z` |
| Structure fraction | DBSCAN core-point test (≥3 neighbours within 60 nm) on trace medians | `sim_real_struct_frac`, overridable |
| σ_loc, per axis | median over traces of the within-trace SD (≥3 locs), **small-n debiased**, reported ± nMAD | manual σ_xy / σ_z |
| Axis ratios | median of the **per-trace** `σ_y/σ_x`, `σ_z/σ_x` | cloud anisotropy (8.1c) |
| Photon distribution | the whole `eco` vector (capped at 200 k) | per-localization σ (8.1c) |
| Trace-length distribution | localizations per `tid`, the whole vector | trace lengths (8.1c) |
| Photons | median `eco` | CRLB photons |
| k, σ₀ | weighted fit of pooled within-trace variance vs 1/eco across photon bins (§1.4) | CRLB k, floor |

**Small-n debiasing.** The sample SD under-estimates σ at small *n*, and the length gate lets the
user push short traces into the median, so each trace's SD is divided by its own correction factor.
The textbook constant `c₄(n) = E[s]/σ` is **the wrong one here**: it un-biases the *mean* of the
sampling distribution, while this aggregates with a *median*, and `median(s) < E[s]` at small *n*.
The matching constant is the median of the same chi distribution,
`m(n) = √(qchisq(0.5, n−1)/(n−1))`; at *n* = 3, `c₄ = 0.886` and `m = 0.833`. On simulated short traces
the m-corrected median recovers σ, where the raw median reads several percent low. Either constant **cancels exactly** in the
axis ratios, since both axes share the same *n*.

The ratios are taken as the **median of per-trace ratios**, not the ratio of medians: within one
trace both axes share a photon count, so the per-trace ratio divides the photon dependence out and
measures cloud shape alone. `median(log(s_y/s_x))` is exactly `log(σ_y/σ_x)` because the two noise
terms are i.i.d. and their difference is symmetric — verified unbiased to +0.11 % (MC se 0.08 %).

`apply_real_rates()` splits the measured total rate into a **constant background offset** and a
**time-varying structure component**, using the structure fraction, and pushes both plus the level.

**What it does *not* transfer:** the structure geometry, the point-cloud shape, the trace-length
distribution, or the labeling. It calibrates the *acquisition clock and the precision*, and you still
have to pick a synthetic shape by hand.

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
from the **unique fluorophore pool**, so a re-detected molecule is drawn at its physical brightness. Poisson shot noise on top.

## 6. Phase 2 — Diffusion / Tracking, in one paragraph

Brownian motion on a sphere / rod / ellipsoid / membrane plane / unbounded volume, surface (2 DOF) or
volume (3 DOF), with free / OU-confined / directed / two-state motion, sampled through a MINFLUX
overlay: cycle clock, geometric dark cycles → irregular Δt, **per-localization** `N ~ Poisson(photons)`
→ `σ = L/(√8·√N)`, photobleaching-limited length, optional re-acquisition survivorship. Emits
`tid, tim, x, y, z, efo, eco` plus a ground-truth sidecar, and scores itself with a built-in CVE
round-trip. Full spec: [simulation_module.md](simulation_module.md).

---

## 7. Limitations

What the simulator does not model, and what that means for results built on its output.

**"% resolved" counts discovery.** On the acquisition timeline a molecule counts as resolved when it is
first detected. The resolution of a molecule is bounded by `ε(n) = √(r² + link² + σ_loc²/n)`: more
localizations shrink only the last term, and the label terms are a floor that no acquisition time
removes. Read the timeline as a detection curve and use the Precision Budget verdict for resolvability.

**One σ for every localization.** Real MINFLUX precision is set per localization by its photon count.
Simulated precision clouds are therefore homogeneous Gaussians, so they cannot exercise the precision-
cloud QC, and inverse-variance weighting has no effect on simulated data.

**The CSV export has no `eco`, `efo` or `cfr`.** Reloaded simulated data cannot run photon-weighted
shape fits (they fall back to unweighted), the fiducial scan, or `efo`/`cfr` filters.

**Background uses the ideal-pattern penalty.** The SBR factor `g` is exact only for the triangular
2D pattern at its centre. It ignores the emitter's offset from the pattern centre and the hexagonal
TCP of commercial systems; a calibrated k absorbs the average of both.

**Re-blinks sit at exactly the same position.** Repeat detections of one fluorophore carry no
between-blink motion or drift, so simulated data cannot calibrate or validate the de-blinking radius.

**Normal-mode timestamps are not physical.** They are a single ramp across all localizations, so
time-domain tools (TSA, cycle time, detection over time, dark-cycle colouring) should only be run on
Acquisition-mode output, which uses a detection profile but still a fixed cadence without dark cycles.

**Crowding changes the timeline only.** It lowers the discovery rate, but crowded molecules get the
same σ and trace lengths as isolated ones, whereas a second emitter in the doughnut biases the real
localization.

**The timeline updates live, the point cloud on Generate.** After a parameter change, the localizations,
RMSE and Precision Budget still show the last Generate until it is pressed again.
