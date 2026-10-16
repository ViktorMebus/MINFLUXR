# MINFLUX Analysis — Methods, Mathematics & Appropriateness

*Reference document for the analysis and simulation modules of MINFLUXR.* Compiled from the source
code and from primary literature (see [Bibliography](#bibliography)).

> **In-app counterpart.** Every analysis panel in the Tracking tab and the SMLM **Analysis** tab carries a
> collapsible *"Method & MINFLUX adaptation"* box with an Origin / MINFLUX adaptation / Caveat summary
> and a 🟢/🟡/🟠 chip. Those boxes summarise the entries here.

This document has **two layers**, both organised around the app's three method families —
**SMLM (imaging)**, **Tracking (SPT)**, and **Simulation**:

- **Part I — Implementation & Mathematics ([§I](#part-i--implementation--mathematics)).** What each tool
  computes: the estimator, its formula as coded, the modelling assumptions, the MINFLUX-specific
  handling, and how to read the output.
- **Part II — Literature & Appropriateness ([§1](#1-why-minflux-data-needs-its-own-treatment) onward).**
  The published source of each method, and a rating of how well it suits the statistical character of
  MINFLUX data.

Every Part I entry names its Part II rating, and vice versa.

Part II ratings:

- 🟢 **GREEN** - well-suited; MINFLUX's properties are neutral or an advantage.
- 🟡 **AMBER** - usable, but an assumption is strained; needs a MINFLUX-aware adaptation.
- 🔴 **RED** - a core assumption is violated; results are biased unless mitigated.

> A rating describes the *textbook method*. Where MINFLUXR adapts the method for MINFLUX, that is
> described under **MINFLUXR's handling**, and several tools are more defensible than their textbook
> rating suggests.

> **Scope.** Analysis methods specific to MINFLUX are thinly covered in the literature. Most tools here
> are MINFLUXR's own adaptation of an established method, and some (marked *tool-original*) were built
> for MINFLUX data. They describe one considered way to handle such data, not a community standard, and
> none has been validated by other groups.

---

# Part I - Implementation

*What each tool computes, as coded. Formulas are transcribed from the source (`server_tracking.R`,
`server_smlm.R`, `R/smlm_*.R`, `sim_server.R`); each entry names its Part II appropriateness rating.
Units in code: positions in **nm**, time internally in **ms** (imported `tim` is seconds → ×1000). D is
carried internally in **nm²/ms** and displayed in **µm²/s**, with `D_um2_s = D_nm2ms / 1e3`
(1 µm²/s = 10³ nm²/ms).*

**Notation.** Symbols used in more than one section; symbols that belong to one formula are defined
where they appear.

| Symbol | Meaning | Unit |
|---|---|---|
| d | number of spatial dimensions (2 or 3) | – |
| D | diffusion coefficient | nm²/ms (shown as µm²/s) |
| σ, σ_loc | localization precision, per axis | nm |
| L, L_TCP | diameter of the targeted coordinate pattern (TCP) of the final iteration | nm |
| N | photon count of one localization (ECO) | photons |
| a | prefactor in σ = L/(a·√N) | – |
| t_c | cycle time: the repeat unit of the localization sequence | ms |
| SW | fixed overhead added once per recorded localization | ms |
| Δt | time between two consecutive localizations of a track, Δt = k·t_c + SW | ms |
| k, n_cyc | number of cycles a step spans (1 = integrated, ≥ 2 = after dark cycles) | – |
| Δt₁ | single-cycle interval, t_c + SW | ms |
| R | motion-blur coefficient of a single-cycle step, typically 0.07–0.13 (§I-B.1a) | – |
| R_abs | blur offset R·Δt₁: every step's effective diffusion time is Δt − 2·R_abs | ms |
| τ | lag time of an MSD point | ms |
| n | number of points (localizations, units or members) in the set analysed | – |

Free Brownian motion with localization error and motion blur gives, per step and for the MSD:

**E[|Δr|²] = 2d·D·(Δt − 2·R_abs) + 2d·σ²**,  **MSD(τ) = 2d·D·τ + 2d·σ² − 4d·D·R_abs**.

## I-A. SMLM (imaging) analysis

The SMLM tab treats each acquisition as a spatial point pattern of localizations, optionally reduced
through a **de-blinking hierarchy** before any spatial statistic:

> **localizations → trace (`tid`) median → de-blinked molecule (`cluster`) median → median cluster
> (`cluster2`)**

The two tiers answer different questions and should not be conflated. **De-blinking** is an instrument
artifact correction: one fluorophore that blinks several times yields several trace IDs, and this tier
merges them back into one molecule. **Median clustering** is biology: it groups distinct molecules into
the structure they form. The word *cluster* in this document refers to the second tier; the first is
called a *de-blinked molecule* (its data column is named `cluster`).

Collapsing to trace/molecule medians is the mitigation for the labelling/blinking over-counting bias
(Part II §6.0). Most estimators can run on raw localizations or on medians (`enable_clustering`).

### I-A.1 Nearest-neighbour (NN) distances — [Part II §6.1]
- **Compute:** `nn_index_pairs(coords, k)` returns each point's `k` nearest neighbours (via `dbscan::kNN`).
  Multiplex/two-colour mode supports `mode = "within"` (same channel) or `"between"` (cross-channel only,
  never same channel) so a bivariate NN distance distribution can be built.
- **Output:** the distribution of the k-th NN distance; median/quartiles per group. On de-blinked medians
  this approximates inter-molecule spacing; on raw localizations it is dominated by the intra-molecule
  blink/label footprint.
- **Cross-channel co-localization curve:** the two-colour panel additionally plots the **empirical
  cumulative distribution** of the cross-NN distances — the fraction of channel-A molecules with a
  channel-B partner within distance d. This is an **object-based** co-localization measure (Lachmanovich
  et al. 2003; Dunn et al. 2011), the right family for MINFLUX, which has coordinates and no image.
- **Caveat (no null):** cross-NN has **no null model built in** — every A point has a nearest B point, so
  the curve is never zero and a small median distance can arise from density alone. Significance must
  come from L₁₂/g₁₂ with their CSR envelopes (§I-A.4).
- **Attribution.** NN distance as a spatial statistic originates with **Clark & Evans 1954**, whose
  index R = d̄_obs / d̄_Poisson is tested against a Poisson null. That index is not computed here,
  because the panel has no null.

### I-A.2a De-blinking — merging re-blinks into molecules — [Part II §6.1]

Operates at **trace level**: a whole trace joins a molecule or it does not, so a trace can never be
split across two molecules. Two criteria, selected by `deblink_method`.

- **Centroid distance** (default, `deblink_merge_centroids()`). Each trace collapses to its median
  position, and centroids are grouped by **complete linkage** at `eps_centroid`: every member of a
  molecule lies within the threshold of every other, so a group's diameter is capped at the threshold.
  (Single linkage would merge through chains and let a group grow well beyond it.) `centroid_k`
  expresses the threshold in units of the measured per-axis trace spread σ and drives it.
- **Majority vote** (`deblink_merge_majority_vote()`). For each localization *p* in trace X, compare
  d_own = |p − centroid of X computed **without** p| against d_alt = |p − nearest other trace centroid|.
  *p* votes to merge when d_alt ≤ ρ·d_own; X and Y merge when a fraction ≥ *f* of X's localizations
  agree (`vote_frac`, `vote_rho`, `vote_link` = mutual/either), with connected components over the
  accepted edges. Leave-one-out is not optional: the centroid is fitted to X's own points, so without
  it d_own is systematically too small and short traces under-merge.

  The vote fraction for two traces of the **same** molecule is **≈ 0.5 by construction**: leave-one-out
  makes d_own and d_alt identically distributed, since both centroids estimate the same position. *f*
  must therefore sit well below 0.5. A fluorophore that blinked three or more times splits its votes
  across siblings, so the decision uses the **total** vote and then wires edges to every candidate
  holding a real share.

- **Merge-only, always.** De-blinking never removes a molecule. A trace that cannot be merged becomes
  its own single-blink molecule (`deblink_assign_singletons()`; the centroid rule guarantees this by
  construction), so `cluster > 0` everywhere and **ID 0 does not occur for this tier**. Removal is the
  job of the metadata and spread filters, where it is visible and reversible. With a too-tight threshold,
  counts visibly inflate toward the trace count.
- **minPts follows the length filter.** Where it applies it is derived from the Trace ID Length filter's
  lower bound (`deblink_minpts()`), so it cannot contradict the filter the user set.
- **Anisotropic option** (`deblink_anisotropic`, `deblink_axis_weights()` / `apply_deblink_weights()`).
  Per-axis whitening by the measured spread, normalised by the geometric mean so the transform is
  volume-preserving and the threshold keeps its nm meaning. A sphere in whitened space is an
  **ellipsoid in real nm**, tighter along the better-resolved axis. Valid for the *merge decision only*
  — reported cluster radii, R_g and neighbour distances stay in true nm.
- **Several channels.** De-blinking runs once per channel, so cross-channel merging is impossible by
  construction. Molecule IDs are offset by a running maximum (channel B starts at max(A) + 1), and the
  channel letter shown (`A-1`, `B-1`) comes from the channel column, never from the ID.

**Choosing the threshold.** Three data-driven aids, in increasing order of trust:
- `suggest_centroid_threshold()` — threshold = *k* × measured per-axis spread. The *k* is the weak
  point: it depends on how far a label moves *between* blinks, which the app cannot measure without
  already knowing which traces are the same molecule. The right *k* grows with that motion, from ≈1.4
  for a static label to ≈3 at 6 nm of between-blink displacement. Default 2.0.
- `deblink_reblink_excess()` — measures *k* **from the data**. Leftover re-blinks sit a few nm apart, so
  after merging they show as an excess of molecule pairs closer than 2σ against a 30–60 nm reference
  shell. 🟠 It is a **lower bound only**: once the excess clears the test goes silent, so it detects
  under-merging and is blind to over-merging.
- `deblink_count_curve()` — 🟠 EXPLORATORY. Molecule count vs threshold; the knee separates "still
  collecting re-blinks" from "now fusing distinct molecules". **Centroid criterion only**: the knee
  detector keys off the steep fragmentation drop that only that curve has.

**Why centroids.** Linking whole traces by their centroids means one stray localization cannot fuse two
molecules, which it can when raw localizations are linked point to point. The thresholds are therefore
not comparable across the two approaches: the centroid threshold is ≈1.5–2σ, while a raw-point radius
sits near 0.6σ.

### I-A.2b Median clustering (DBSCAN / HDBSCAN / OPTICS) — [Part II §6.1]
- **Median clustering** (`cluster2`): density clustering run on de-blinked molecule medians (or trace
  medians when de-blinking is off) → groups molecules into structures. Unlike the de-blinking tier,
  **ID 0 means noise/unassigned here and must be excluded before counting.** With several channels,
  each channel is clustered separately and IDs are offset by a running maximum.
- **eps tuning:** an optional **kNN-distance plot** (sorted distance to the k-th neighbour) exposes the
  density "knee" that `eps2` should sit at. This diagnostic belongs to *this* tier only.
- **Density field:** `compute_spatial_density()` uses `dbscan::frNN(pts, eps = r)` to count neighbours
  within radius `r` per point (fixed-radius local density used for colouring and the acquisition
  structure/background split).

#### Membership refinements (`run_cluster2_grouping`)

Three optional steps run after the chosen method, in this order: **core-only → prune → min members**.
All are off by default; ordering matters, because the member-count filter must count what survived.

- **Core points only (DBSCAN\*)** — classic DBSCAN assigns *border* points: centroids within `eps` of a
  core point but sparse themselves, so a lone background centroid near the rim of a structure is drawn
  as a member. Keeping only core points is the cluster definition of Campello et al. 2013, which is why
  HDBSCAN and OPTICS need no equivalent (the flag is ignored for them).
- **Sparse-member prune** — graded form of the same criterion. The statistic is each member's **core
  distance** (the radius needed to reach `minPts` fellow members), compared with the **median core
  distance of that member's own neighbourhood**. The cutoff is relative (`tol` × baseline; no absolute
  nm).
  - *Why core distance.* Distance to the centroid or to a density maximum assumes
    a convex, roughly isotropic blob, and would amputate the ends of a rod or the far side of a ring,
    which are legitimately far from any centre. Core distance asks about local density instead.
  - *A local baseline.* A cluster with an internal density gradient is not judged by its densest region,
    so a genuinely looser periphery is kept.
- **Min members per Median Cluster** — drops whole clusters below a member count. Kept separate from
  `minPts`, so density and structure size are set independently. The prune
  cannot do this job: it skips clusters with <4 members outright, and because it is scale-free a tight
  spurious 5–15 member cluster is internally consistent and loses nothing.

🟡 **These change measurements.** Demoted members become `cluster2 = 0` and leave
curvature, sphere radius, R_g, alpha-shape volume and the NN-vs-size panels. Tightening membership
shrinks measured structure size systematically — report which refinements were on.

⚠️ **Identifiability limit.** No density method separates a genuinely diffuse periphery of a structure
from background centroids; in density terms they are the same object. Every setting here is a stated
choice about what counts as belonging.

### I-A.2c Localization precision cloud (QC)

Pools every localization minus **its own trace centroid** and recentres at (0,0): the empirical
point-spread of the localization precision, with each molecule's true position subtracted out. Its RMS
radius is σ_loc; the *shape* is the diagnostic.

- **Rings are covariance ellipses, rotated.** Semi-axes and orientation come from the eigen-decomposition
  of the 2×2 covariance (`cov_ellipse_params()`), drawn as an SVG path because plotly's `circle` shape is
  defined by a bounding box and can only ever be axis-aligned.
- **Anisotropy is rotation-invariant.** An axis-aligned ellipse cannot represent a *tilted* cloud:
  diagonal elongation inflates σx and σy equally, so a max(σx,σy)/min(σx,σy) ratio reads a tilted 3:1
  cloud as isotropic. Drift within a trace is directional, so a tilted cloud is its expected signature.
  Anisotropy is therefore the **eigenvalue ratio** plus the **tilt angle**, which says *which way* the
  drift runs.
- **Two populations, both drawn.** Dotted rings are the *pooled* per-localization covariance (matching the
  heatmap, so long traces weigh more); the dashed ring is a *robust* per-trace covariance
  (`robust_median_cov()`, element-wise median of per-trace matrices with a positive-definiteness guard,
  one vote per trace). When the two separate, that gap **is** the "a few long/wide traces dominate"
  diagnostic — usually fiducials.
- **Caveat:** the spread measured here is not the microscope's precision. It is
  σ²_photon + σ²_linker + σ²_motion, so a σ of 5–7 nm on antibody-labelled biology is expected and does
  *not* contradict a 1–4 nm instrument specification.

### I-A.2d Axial (z) calibration

MINFLUX z is multiplied by a refractive-index mismatch correction. This is a **property of the dataset**
(immersion and sample indices, imaging depth), so it is user-set in both tabs
(`z_scale`, `trk_z_scale`; shared default `TRK_Z_SCALE_DEFAULT = 0.7`).

- **One factor, applied once.** Readers are called with `z_scale = 1.0` so `.npy`/`.json` arrive raw;
  CSVs are returned to raw by `strip_csv_z_calibration()` first. The user's factor is then applied at a
  single point per tab. Without this, formats drift apart and a CSV gets the correction applied twice.
- **CSV provenance** (`detect_csv_z_calibration()`): a `z_scale_applied` column (written by this app's
  exports) is **exact**; a `dataset_tag` marks an older MINFLUXR export; anything else is assumed to be a
  pyMINFLUX export. 🟡 **Cases 2 and 3 are inferred, not measured** — nothing in those files records the
  factor, and a CSV written with z-scaling disabled would be misread. The app says which branch it took.
- **It multiplies every axial quantity** — z trace spread, cluster extent in z, 3D neighbour distances,
  R_g. Applying it twice compresses z by 30%; the guard above is what prevents that.
- **Applied on Plot Data**, so a saved analysis reproduces. Greyed out when the dataset has no axial
  variation (`has_axial_data()` tests for variation, since 2D MINFLUX exports carry z = 0).

### I-A.2e Sample drift & jitter (QC) — [Part II §6.8]
`sample_drift_estimate(df, min_n, radius, conf, group)`. Estimates how far the sample moved across the
run, with an interval.
- **Lever arm = the gap between blinks.** A trace lasts ~0.15 s, so within-trace slopes
  have no power against slow drift. Instead, pairs of traces
  close enough to be re-blinks of one molecule are taken, and the **signed displacement is regressed on
  the signed time gap through the origin**, giving a velocity with the whole acquisition as lever arm.
  With de-blinked IDs as `group`, a molecule seen several times contributes a trajectory.
- **Different molecules do not tilt the slope.** A pair that is really two molecules has isotropic,
  zero-mean displacement, so it widens the interval without biasing the rate; the pair radius (default
  10× the trace spread) can therefore be generous.
- **Drift vs jitter.** Drift is the time-aligned, accumulating part; jitter is the residual scatter;
  `frac_directed` reports which dominates.
- **The one real bias:** selecting pairs by separation truncates the large displacements that fast drift
  produces, so on a badly drifting file the estimate is a **lower bound**.

### I-A.2f Fiducial / contaminant scan (QC) — [Part II §6.8]
`observeEvent(input$run_fiducial)` in `server_smlm.R`. Per trace: localization count, duration, median
`efo`/`cfr`, centroid, spread.
- **Lifetime (primary).** A fiducial does not bleach, so the sequence keeps re-finding it and its trace
  outlives every real fluorophore in the same acquisition. Flag when
  **duration ≥ `fid_dur_x` × q99(duration)**, default 5×, with the 99th percentile taken over the traces
  of that file that live at all (duration > 0). The threshold is relative to the file, so it adapts to
  the acquisition, and a file of mostly single-localization traces cannot collapse the reference.
- **Brightness (secondary).** Robust z of log10 `efo` ≥ `fid_efo_z` (default 4) **and** at least
  `fid_min_locs` localizations (default 5). The localization guard excludes 1–3 localization blips,
  which can be very bright without being fiducials.
- Length-outlier z, recurrence and spread are tabulated as evidence only.

**Report-only** unless *Exclude flagged* is ticked, which then removes the traces upstream of every
analysis and is saved in the settings JSON.

### I-A.3 Ripley's K / L(r) — [Part II §6.2]
`compute_ripley_L(x,y,z, r_seq, use_3d)`:
- Random **subsample to ≤2000 points** for tractable O(n²) pairwise distances (`stats::dist`).
- Observation window = **axis-aligned bounding box**; `area_vol` = box area (2D) or volume (3D).
- For each radius r, with |A| the window area (volume in 3D) and P(r) the number of point pairs closer
  than r, **K(r) = |A| · 2·P(r) / (n(n−1))**; the 2 counts each pair in both directions.
- **L-transform:** `L(r) = √(K/π)` (2D) or `L(r) = (3K/4π)^{1/3}` (3D). Plotted as **L(r) − r**: 0 under
  complete spatial randomness (CSR), **>0 clustering**, **<0 dispersion**.
- **Window and edge correction (two switches at the top of the tab).** By default the window is the
  bounding box **with edge correction on**. Without it L(r) − r is biased downward near large r. With
  *Apply boundary (edge) correction* on, each pair is weighted by the **translation correction**
  1/γ̄(d) from the window's set covariance (Illian et al. 2008; Baddeley et al. 2015). With *Replace the
  bounding box with the sample footprint* on, the window and its area are the estimated support
  (§I-A.3a), so a cell that fills only part of its frame is not read as clustered. The two interact:
  edge-correcting a box that contains a lot of empty space is worse than not correcting, so on a concave
  or half-in-frame cell use both together. A bounding-box run reports how much of the box the footprint
  covers, and the window readout prompts for the footprint window below 70%.
- **Auto r_max is capped by the shortest window axis** (`min(extent)/4`, z included in 3D), since edge
  bias is set by whichever dimension runs out first. The same cap applies to L₁₂ and g₁₂. An explicit
  r_max from the UI is always honoured.
- **Monte-Carlo CSR envelope (`compute_ripley_envelope`, optional):** generates `nsim` CSR patterns of
  the *same n* in the **same frozen window** (`window_ranges`, a pre-supplied 2×d bounding box shared by
  the observed curve and the sims) and returns the pointwise 2.5–97.5% band + mean null curve. Because
  the sims inherit the identical window, the edge bias is present in the band too and **cancels when the
  observed curve is read against it**, so the envelope is the reference (the flat line is relabelled "theoretical CSR (uncorrected)"). A two-sided α≈0.05 test needs
  `nsim ≥ 39`; default 99, capped 19–199. In all-files mode the envelope is computed for the first few
  files only (cost = `nsim × files × O(n²)`).

#### CSR envelopes: what they test
`compute_ripley_envelope` + `envelope_global_p` (`R/smlm_spatial.R`), used by L(r), g(r), L₁₂ and g₁₂.
- **Pointwise band (2.5–97.5%)** at each r. It is descriptive: it holds at ONE radius, so a curve
  leaving it *somewhere* over many radii is not a 5% test, and random data does so in most replicates.
- **Global band + p**, the studentized maximum-deviation test (Loosmore & Ford 2006; Myllymäki et al.
  2017): T = max_r |f(r) − mean(r)| / sd(r), band = mean ± T_(1−α) · sd,
  p = (1 + #{T_sim ≥ T_obs})/(nsim + 1). Crossing it anywhere IS a 5% test. Simulated T uses
  leave-one-out mean/sd so it is on the same footing as the observed curve. Drawn dashed per file, with
  the p-values captioned above each plot. The smallest reachable p is 1/(nsim + 1): at nsim = 19,
  p = 0.05 means no simulation deviated as far as the data.
- **Window debias.** The observed curve is normalised by an area estimated from the observed points,
  and the simulations live inside that same estimate. At sparse labelling the occupancy footprint stops
  inside the true cell edge, so the null would be denser than the data and random points would read as
  regular. `debias_support_window()` measures that shrink on the window itself (draw n points in the
  mask, re-estimate, average the ratio) and dilates the mask by whole cells to the implied area. It is
  applied to the Ripley, g, L₁₂, g₁₂ and co-localization window, and not to the hotspot panel, whose eps
  calibration randomises inside the same mask and is self-consistent.
- **Two-channel null — three choices, three different questions.**
  - `"shift"` (default): channel A stays where it is and channel B is translated as a rigid whole by a
    random vector, wrapping inside the box (toroidal shift; Lotwick & Silverman 1982; `.shift_channel`).
    Each channel keeps its own clustering exactly and only their registration is broken, so this is the
    null for **"are A and B placed independently"**, the co-localization question. Shifts shorter than
    the largest radius analysed are rejected (a near-zero shift keeps part of the association and costs
    power silently); with a footprint, a shift is accepted when ≥ 90% of B lands in the cell, otherwise
    the best of 100 candidates is used, and points outside the cell are dropped. In a concave cell the
    90% rule is rarely met, so the best-of-100 path is the usual one; it is what was validated.
  - `"labels"` keeps the observed points and permutes which channel each belongs to. It tests **"do the
    two channels differ in arrangement"**. Two species that share the same clusters ARE a random
    labelling of one pattern, so it cannot detect them; two unrelated clustered species read as
    exclusion, because the pooled pattern is more clustered than either cross pair.
  - `"csr"` draws independent uniform patterns in the window. It ignores each channel's own clustering
    and calls almost any pair of clustered channels co-localized.
- **Calibration** (`tests/validation/coloc_null_clustered.R`; clustered channels in an L-shaped cell,
  footprint window and edge correction on, 60 replicates × 99 simulations, global test at 5%):

  | truth | null | L₁₂ rejects | g₁₂ rejects |
  |---|---|---|---|
  | independent clustered species | shift | **6.7%** | **3.3%** |
  | | labels | 100% (all "exclusion") | 100% (all "exclusion") |
  | | csr | 97% | 100% |
  | same clusters | shift | 100% | 100% |
  | | labels | **1.7%** | **3.3%** |
  | B within ~8 nm of A | shift | 100% | 100% |
  | | labels | 100% | **5%** |

  On independent **unclustered** (CSR) channels, the only case where random labelling is exact, it holds
  1.7% (L₁₂) and 3.3% (g₁₂), and the CSR-in-window null 21.7% and 8.3%
  (`spatial_envelope_coverage.R`).

### I-A.3a Sample footprint & labelling-density thinning (opt-in) — [Part II §6.2b]
`R/smlm_support.R`.
- **Unit:** `unit_centroids()`, one median position per trace ID or, preferably, per de-blinked molecule
  (a trace ID is one blink, so trace-ID density mixes expression with photophysics).
- **Footprint:** `estimate_support(mode = grid | hull | bbox | kde)`. The default is a morphological
  occupancy estimate (grid occupancy closed over a gap-fill radius); `kde` takes a highest-density
  region. For multi-file panels the gap-fill radius is **shared across the file set**
  (`spatial_resolve_close`): each file's closing sweep is read with a bridging-jump rule (the largest
  relative jump is the gaps between structures closing) and the maximum across files is used, so files
  with identical geometry are given the same window. `frac_of_bbox < 40%` is
  flagged as a possibly structure-hugging window.
- **Thinning:** `thin_units_to_p()` keeps each whole unit independently with probability p (never single
  localizations), seeded per file. Independent p-thinning sends λ → pλ and ρ₂ → p²ρ₂, so **g(r) and K(r)
  are exactly invariant** (Illian et al. 2008) while fixed-threshold statistics (NN distance, DBSCAN at
  fixed eps/minPts) are brought onto a common density.
- **Target = global density n / support**, which clustering cannot change. `nn_local_intensity()` is a
  diagnostic only: aggregates shorten NN distances, so an NN-based target thins a genuinely clustered cell
  harder and removes part of the signal.
- **Error budget:** the support estimate, above all its gap-fill radius, dominates. Thinning is one-way.

### I-A.3b Local density & hotspots — [Part II §6.2c]
`R/smlm_hotspots.R`. Dimensionless per-cell statistics that read the **tail** of the
local-density distribution, for samples where a minority of molecules sit in foci. A median NN distance
mostly tracks overall density, so it cannot separate a clustered cell from a denser uniform one; these
statistics can.
- **Local enhancement** `compute_local_enhancement`: k-NN density **λᵢ = k/(ω_d·r_k^d)** (ω₂ = π,
  ω₃ = 4π/3; Loftsgaarden & Quesenberry 1965) divided by the cell mean **n/support**. Default k = 5.
  The footprint shrinks on a sparse cell, which *lowers* enhancement, so a positive result there is
  conservative.
- **Index of dispersion** `compute_dispersion_index`: var/mean of quadrat counts inside the window,
  **including empty quadrats** (1 under CSR; Fisher, Thornton & Mackenzie 1922). Only quadrats lying
  **entirely inside the footprint mask** are counted (an area-based fallback is flagged
  `edge_quadrats`), because a partial edge quadrat has a smaller expected count and inflates the ratio.
  Scale-dependent by design: always report the quadrat size.
- **Calibrated DBSCAN eps** `calibrate_dbscan_eps`: for each candidate eps, scatter n units uniformly
  **inside this cell's footprint** `nsim` times and take the largest eps whose random-clustering fraction
  stays ≤ `target_frac` (5%). `hit_floor` flags non-convergence. When even the grid floor (40 nm)
  over-clusters, as in a very dense cell, the grid is extended downward in 20% steps (to ~3.4 nm,
  `grid_extended`).
- **Foci** `compute_foci_table`: DBSCAN at that eps; per focus n, Bessel-corrected Rg and anisotropy.
- **What survives calibration:** the **size of the largest focus**. Focus count and fraction of units in
  foci remain density-sensitive, because the calibrated eps grows in sparse cells.
- **Largest focus relative to chance:** the raw `max_focus_n` of a random pattern grows with density,
  so it is confounded across cells. The calibration records the largest focus in ≥ 20 random patterns in
  the cell's own footprint at the calibrated eps; `max_focus_ratio = max_focus_n / that 95th percentile`
  (1 = what chance reaches). This is the metric to compare between cells.
- **Calibration.** On random patterns across density and footprint shape, cv_ratio, VMR, the fraction
  above the enhancement cut and `max_focus_ratio` hold their chance baseline, while raw `max_focus_n` and
  focus count move with density (`tests/validation/results/`). Power is weak in sparse cells.

### I-A.4 Bivariate Ripley L₁₂(r) & pair cross-correlation g₁₂(r) — [Part II §6.2]
- **L₁₂:** cross-type K from the n₁×n₂ distance matrix. Write |A| for the window area (volume in 3D)
  and P(r) = #{(i,j) : d_ij ≤ r} for the number of cross pairs within r — every pair is counted once,
  because the two channels are distinct point sets and the matrix is already full, so there is no
  ×2 as in the univariate case. Without edge correction

  **K₁₂(r) = |A| · P(r) / (n₁·n₂)**

  With edge correction on (the default), each pair at distance d contributes 1/γ̄(d) in place of
  1/|A|, γ̄ being the window's isotropised set covariance, so |A| in the numerator is replaced by
  |A|²/γ̄(d) and the weight is applied per shell. Same L-transform as the univariate case,
  `L₁₂ = √(K₁₂/π)` (2D) or `(3K₁₂/4π)^{1/3}` (3D). Plotted as L₁₂(r) − r: >0 attraction,
  <0 segregation.
- **g₁₂ (PCF):** the *derivative* form. Count only the cross pairs whose distance falls in the shell
  (lo, hi] = (r−dr/2, r+dr/2], call that P_shell(r), and normalise by the shell's own area or volume
  V_shell = π(hi²−lo²) in 2D, (4/3)·π(hi³−lo³) in 3D:

  **g₁₂(r) = |A| · P_shell(r) / (n₁·n₂·V_shell(r))**

  |A| is the window and V_shell the shell: the two are different areas, and only the ratio of the
  expected to the observed pair density is dimensionless. The same per-shell edge weighting applies.
  **g₁₂ = 1 under CSR**, >1 attraction, <1 repulsion.
  PCF is scale-local (not cumulative like K), so it localises the interaction distance more sharply.
- Same ≤2000 subsample and the same two window switches as L(r) (bounding box, edge-corrected, by
  default).
- **Envelope:** the same `compute_ripley_envelope` (with `estimator = "L12"` / `"pcf"` and both `n₁`,
  `n₂`), with the global band and p-value of §I-A.3. The default null is a random shift of channel B
  (§I-A.3, two-channel null); random labelling and CSR-in-window are the alternatives, and the plot
  caption names the null used. For L₁₂ the observed curve above the band
  = co-localization at that scale, below = exclusion; for g₁₂ the band **brackets g = 1** and the
  observed curve rising above it at short r = attraction.

### I-A.5 KDE density map + Pearson co-localization — [Part II §6.4]
- `MASS::kde2d` on each channel over a **shared grid** with a common bandwidth.
- **Bandwidth (auto):** Silverman's rule `bw = 0.9·min(sd, IQR/1.34)·n^{-1/5}`.
- **Co-localization index:** Pearson r of the two flattened density surfaces. **Explicitly bandwidth-
  dependent** — small bw → r≈0 even for real co-localization, large bw → inflated; treat as a *relative*
  measure at a fixed bandwidth (stated in-UI).
- **Bandwidth checks:**
  - **`bw_sweep`** recomputes r over a log-spaced range of bandwidths (0.4–20% of the field span) and
    returns `data.frame(bw, pearson_r)`, plotted as r(bw) so the smoothing-dependence is explicit.
  - **`costes_nsim`** randomizes Ch2 `nsim` times, recomputes r at the chosen bandwidth, and returns a
    one-sided **p-value** = (#{random r ≥ observed}+1)/(nsim+1) plus the null mean/SD.
  - **`null_mode` selects WHICH null** — this is what decides whether the p-value means anything:
    - **`"shift"` (default): toroidal random shift** of Ch2 — translate it by a random offset with
      wrap-around inside the same window. Ch2's own clustering, density and point count are preserved
      **exactly**; only its spatial registration with Ch1 is destroyed. This isolates cross-channel
      dependence, and is the same random-shift logic (Lotwick & Silverman 1982) already behind the
      bivariate statistics.
    - **`"csr"`: uniform resampling** of Ch2, a test against *complete spatial randomness*. Kept for
      comparison and labelled in-UI.
  - The UI points users to **g₁₂(r)** for a bandwidth-free, scale-resolved alternative.
- **Window matters for the p-value.** Both nulls randomise over the analysis window. With the bounding
  box (the default), channels that share one cell inside a larger frame are contrasted with a Ch2 that
  partly lands outside the cell, so both nulls are anti-conservative there; the random shift less so
  than CSR. With *Replace the bounding box with the sample footprint* on, the null is drawn inside the
  estimated support (`sample_in_support()`) and becomes conservative: "not significant" can be quoted,
  while a positive claim should come from g₁₂(r) against its envelope. The random-shift null follows
  Lotwick & Silverman 1982; it is not the block randomisation of Costes et al. 2004.

### I-A.6 Radius of gyration & gyration-tensor anisotropy — [Part II §6.3]
- **Rg:** `Rg² = mean(‖p − p̄‖²) · n/(n−1)` per cluster (second moment about the centroid, Bessel-
  corrected). Points are **median member positions**, never raw localizations.
- **Small-N bias correction.** The centroid is estimated from the same points whose spread is
  measured, so it sits inside the cloud and **E[Rg²_obs] = Rg²·(n−1)/n**: a 29% underestimate at n = 2,
  11% at n = 5, 1% at n = 50. That bias depends on n, so uncorrected it manufactures a "larger clusters
  have larger Rg" trend whenever clusters of differing member count are compared. The `n/(n−1)` factor
  (a sample-variance Bessel correction) removes it. At n = 2 the corrected mean still sits ~12% low,
  because Rg = √(Rg²) and E[√x] < √E[x], so treat 2-member clusters as indicative only.
  `bias_correct = FALSE` returns the raw configurational Rg of the observed point set.
- **Anisotropy:** eigenvalues λ₁≥λ₂≥λ₃ of the gyration tensor give **Westin linear / planar / spherical**
  coefficients (c_l = (λ₁−λ₂)/Σλ, c_p = 2(λ₂−λ₃)/Σλ, c_s = 3λ₃/Σλ) and the **Rudnick–Gaspari asphericity**
  — used to *infer* cluster shape (linear filament vs planar vs isotropic) without the user declaring it.

### I-A.6b Per-structure geometry table (`estimate_cluster2_curvature` / export)

One row per Median Cluster. **Members are de-blinked molecules when de-blinking has run, else trace
IDs** (`cluster2_member_col()`, the same rule as `smlm_nn_unit()`), each reduced to its median
position first. The exported table reports `n_traces` (blinking events), `n_molecules` (NA when
de-blinking is off) and `n_members` (the unit the geometry actually used) separately, plus
`member_unit`, so a count can never be mistaken for the measurement basis.

- **Curvature is MEAN curvature H = 1/r_s in nm⁻¹**, with r_s the fitted sphere radius. 🟡 For a fitted sphere it
  is a one-to-one function of `sphere_radius` and carries **no information beyond it**; the column exists
  because bending is conventionally quoted this way.
- **`volume` + `volume_alpha` (`cluster_alpha_volume`).** Alpha is a probe radius, so a value below the
  members' own spacing returns an *empty* shape. Alpha is walked up a ladder built from the group's own
  median NN spacing, taking the smallest radius whose shape passes three tests: **one connected
  component**, **watertight**, and **containing every member** (`inashape3d`). Connectivity alone fails
  on filled structures, because denser sampling connects at a smaller α; connectivity plus
  watertightness still admits a degenerate sliver such as one tetrahedron spanning four members, which
  the containment test rejects. The radius used is reported alongside. Hollow shells of known radius are
  recovered to ~96%.
  - 🟡 *Strong downward bias at low member counts.* Alpha shapes interpolate between sampled points with
    flat facets, so they under-fill a curved surface and the deficit grows as sampling thins. On a
    synthetic R = 150 nm shell (median of 25 draws), the recovered fraction of the true volume is:

    | members | 10 | 15 | 20 | 31 | 50 | 80 | 150 | 300 |
    |---|---|---|---|---|---|---|---|---|
    | volume / truth | 30% | 44% | 55% | 68% | 79% | 86% | 92% | 96% |

    Below ~20 members the number is not usable as a volume. Read `volume` against `n_members`, always.
  - 🟡 *What the number is:* for a hollow structure this is the volume **enclosed by** the labels once
    α bridges them — the compartment they bound, not the volume of labelled material. It grows with α,
    so `volume` and `volume_alpha` must be reported together.
- 🟡 **Small-N bias in `sphere_radius`.** The geometric-median centre is fit to the same points the
  radius is measured from, so it is biased low on few members (about −2% at 40 members, −0.1% at 1000,
  on a synthetic shell). Same family
  as the R_g small-N bias (§I-A.6); left uncorrected, since the magnitude is well under the labelling
  uncertainty.

#### Eccentricity (Ecc1/Ecc2/Ecc3 in the metrics overlay)

Half-ranges `d ≥ e ≥ f` are taken along PC1/PC2/PC3, and each Ecc is the 2-D ellipse eccentricity of one
principal cross-section: **Ecc1 = PC1 vs PC2, Ecc2 = PC1 vs PC3, Ecc3 = PC2 vs PC3**.

- **Ecc3 = 0 is expected for a prolate structure.** It is `√(1−(f/e)²)`, which
  is zero exactly when the two *smaller* axes are equal — i.e. the structure is rotationally symmetric
  about its long axis (a cigar, a rod, a spherocylinder). Ecc3 ≈ 0 while Ecc1, Ecc2 > 0 therefore
  *reports* that only one axis is distinct.
- ⚠️ **Eccentricity has a large small-N floor and is biased away from zero.** PCA assigns axes by
  *sample* variance, so the largest sampled direction always becomes PC1 and the smallest PC3 — the
  ordering itself manufactures apparent elongation. On a synthetic **sphere** (all axes equal, true
  eccentricity 0), median over 400 draws:

    | members | Ecc1 | Ecc2 | Ecc3 |
    |---|---|---|---|
    | 10 | 0.538 | 0.806 | 0.691 |
    | 20 | 0.424 | 0.666 | 0.563 |
    | 40 | 0.335 | 0.526 | 0.417 |
    | 80 | 0.267 | 0.414 | 0.339 |
    | 200 | 0.227 | 0.336 | 0.247 |
    | 500 | 0.179 | 0.267 | 0.203 |
    | 2000 | 0.122 | 0.184 | 0.135 |

  A sphere with 40 members reports Ecc2 ≈ 0.53. **An eccentricity is evidence of elongation only if it
  clearly exceeds the row for its own member count.** No null is computed in-app; compare by hand.

### I-A.6c Median-cluster isosurface (`cluster_density_field`)

A **Gaussian KDE over the member centroids**, contoured to wrap a structure so it reads as an object.
It is a visualization.

- **Built from members.** The field uses member centroids (de-blinked molecules, else
  trace medians), per channel, so a long trace cannot dominate the density.
- **The level is a surface offset in nanometres**: how far outside the member centroids the skin sits.
  Probe points on spheres of radius `offset` around the members are evaluated in the field, and the
  median over probes that sit at that distance is the contour, so the skin reaches about the requested
  offset. The level is read off-grid, so the box size cannot leak into the surface, and the probe set is
  a deterministic Fibonacci sphere, so a redraw cannot jitter it.
- **Default wrap = 2 × the structure's own median member spacing**, and **σ = max(0.5 × wrap, spacing)**.
  At one spacing the kernels around neighbouring members barely touch and the surface breaks into
  fragments; at two spacings it closes into one connected surface.
- **One dial: padding (nm)**, added to (or, negative, subtracted from) the wrap. Bandwidth and grid
  extent are derived from it, and the box is sized to members + (wrap + 2σ), so no second control can
  fight the first.
- 🟡 **Sparse structures cannot be wrapped tightly.** The wrap follows member spacing, so a sparsely
  labelled structure wraps loosely, and negative padding beyond about −10 nm fragments the surface.

- 🟡 **An isolated member still draws a small speck**, much smaller than the structure's rim. Remove
  stragglers upstream with core-only / prune / min members (§I-A.2b) if they must not appear at all.
- 🟡 **Not a density estimate.** The kernel is unnormalised and the mass integral is a voxel sum on a
  coarse grid, so the field is not in localizations/nm³. Use the Ripley/PCF panels for density
  statistics.

### I-A.7 Cluster shape fitting (3D viewer) — [Part II §6.5, 🟠 EXPLORATORY]
- **Dictionary:** sphere, ellipsoid, helix, ring, arc. Each is fit by **orthogonal-distance MLE**
  (minimise Σ of squared point-to-surface distances). The ellipsoid projection solves the exact
  nearest-point KKT condition by vectorised Newton iteration in the PCA frame
  (`mle_project_ellipsoid`), with a radial fallback if Newton diverges.
- **Model selection** (`viewer_rank` = AICc | BIC | RMSD): under an iid isotropic-Gaussian residual model,
  σ̂² = RMSD², so **−2logL = n·log(2πσ̂²) + n**; **AIC = −2logL + 2k**, **AICc = AIC + 2k(k+1)/(n−k−1)**,
  **BIC = −2logL + k·log n**, with k = (geometric DOF + 1). AICc is the default (penalises complex shapes
  hard at small n); it → ∞ when n ≤ k+1 (too few points), which correctly refuses to over-fit tiny clusters.
- **Input units:** fits run on **de-blinked trace centroids**, so one long trace cannot dominate the
  residual sum. The 2D view is a PCA projection onto PC1×PC2 (Pearson 1901).
- **Why exploratory:** the shape is imposed by the dictionary and "best" is a model-selection outcome, not
  proof of biology; with few blink-inflated points several shapes are statistically indistinguishable.
- **Centre estimate — geometric median (Weiszfeld 1937).** `geometric_median()` returns the point
  minimising Σ‖pᵢ − m‖ (the L1 / spatial median), used by every `fit_*_geometry` heuristic and as the
  MLE starting centre (`robust_center`). It has breakdown 0.5, so one loosely attached member cannot
  move it, and unlike a coordinate-wise median it is **rotation-equivariant**, which matters because
  these centres define ring normals and ellipsoid axes through `prcomp()`. On a clean, well-sampled
  hollow shell a surface centroid is slightly more accurate, but that error enters the radius only at
  second order.
  - ⚠️ *Not a centre of curvature.* For a partially labelled structure (a cap of a vesicle, an arc of a
    ring) no centre-of-the-points statistic recovers the geometric centre: it sits inside the observed
    cap. Use the MLE pathway, which optimises the centre against the shape model, when the centre is
    itself the measurement.
- **Photon (`eco`) weighting (`viewer_weight_eco`, default on).** The unweighted fit assumes **iid
  isotropic Gaussian** residuals with a single σ̂, but MINFLUX precision scales as 1/photons. A trace
  centroid averages that trace's localizations, so its variance goes as 1/Σ`eco`, and the trace's
  **total photon count is its inverse-variance weight**. The `mle_fit_*` optimisers minimise
  `sum(w·(L−r)²)`, and the reported RMSD is the weighted √(Σw·d²/Σw).
  **Effect on the criteria:** with per-point variance σ²/wᵢ the log-likelihood gains a −Σlog(wᵢ) term
  that is *identical for every shape*, so it cancels in all AIC/AICc/BIC **differences** and model
  ranking stays valid. Only the absolute values shift, so **never compare criterion values across
  weighting modes**. `n_eff` (Kish effective sample size, `(Σw)²/Σw²`) is reported in the viewer header,
  because weighting trades bias for variance.
  ⚠ **Filter double events before weighting.** A high `eco` means better precision, *and* it is the
  signature of two fluorophores emitting in the same spot. Inverse-variance weighting hands the most
  influence to the brightest traces, so on unfiltered data it amplifies double-emitter artefacts; it is
  only sound downstream of the `eco`/`cfr` metadata filters. **A large weighted-vs-unweighted difference
  is a data-quality alarm.** Quote the weighting mode and `n_eff` with any shape
  parameter, and treat AICc differences of a few units as indistinguishable either way.
  The weight assumes `eco` is comparable across traces, which holds within one sequence and laser power
  but not across files acquired at different powers: check the `eco` distributions before pooling
  weighted fits across such files.
- 🟡 **No parameter uncertainty.** Each shape parameter is a single point estimate with no error bar, so
  quote it as a description of the structure.

### I-A.8 Ground-truth position/count estimator — [Part II §6.6, 🟠 EXPLORATORY]
`estimate_gt_from_labeling(df, preset, eff)` — a **labelling-model inverse** (`R/smlm_gt.R`):
- **Label preset** parameters: `n_fluor` (mean dyes/molecule), `link_nm` (linker length), `dye_spread_nm`
  (dye scatter σ within the complex), `footprint_nm`, `blink_lambda` (extra re-appearances), `eff`
  (labelling efficiency).
- **Forward model it inverts** (the Simulator's, `docs/simulator.md` §1.2–1.4): molecule → K =
  max(1, Poisson(`n_fluor`)) fluorophores, K′ ~ Binomial(K, eff) of them detected, each displaced by a
  half-normal linker (mean ≈ `link_nm`), scattered by `dye_spread_nm`, and yielding 1 + Poisson(`blink_lambda`)
  traces. `gt_tid_model` gives the exact distribution of the per-molecule trace count T: the multiplicity
  **m_tid = E[K′]·(1 + blink_lambda)**, its variance, the detection probability **η = P(K′ ≥ 1)**, and the
  tail P(T ≥ t | detected).
- **Inverse:** build candidate units (trace or de-blinked molecule medians, ≥`minpts` locs), then **merge**
  unit centroids within a **complex radius** `complex_r` (auto = max(footprint, 2.5·dye_spread)) via
  `dbscan(eps = complex_r, minPts = 1)` so several dyes/blinks of one molecule collapse to **one candidate
  GT position** (this produces the contribution fan-out lines). Optional **stoichiometry-guided merge**
  further collapses to k = round(ΣN_tid / m_tid) via hierarchical clustering — only ever *reduces* toward
  the labelling-predicted count, never splits.
- **Two counts per structure:** spatial (candidates / η) and stoichiometric (total traces / m_tid).
- **Consistency test:** given n spatial candidates, the trace total should be a sum of n draws from the
  detected-molecule distribution, so **z = (N_tid − n·E[T|det]) / √(n·Var[T|det])**, with the variance
  floored at the mean (real blink counts are more variable than Poisson). z > 2.58 = *under-resolved*
  (candidates hold several molecules; sets the crowding flag), z < −2.58 = *over-split* (re-blinks or
  dyes of one molecule kept apart). Each candidate also gets P(T ≥ its trace count | one molecule);
  below 0.01 it is marked as holding more traces than one molecule gives. This catches merges the
  spread test cannot: de-blinking caps a unit's width at its threshold, so two molecules a few nm apart
  make a unit no wider than one.
- **Output:** candidate GT positions coloured **green→red by positional uncertainty** (`pos_sd_nm`), and a
  **count range** [lo, hi] spanning both counts, each ± 1.96 SE (spatial: binomial detection error;
  stoichiometric: √(n·Var T)/m_tid). The range does not include uncertainty in the preset itself.
  Checked against the simulator (`tests/testthat/test-smlm-gt-estimator.R`); the counts are only as good
  as the labelling preset or the measured multiplicity.

### I-A.9 Detection-over-time / crowding / NN-vs-size — [Part II §6.7, 🟠 EXPLORATORY]
Bespoke descriptive diagnostics (detection-rate vs time, crowding vs timing, NN distance vs cluster size).
**Correlational, not causal** — detection timing and crowding co-vary with photobleaching/blinking/labelling,
so a trend is a hypothesis to test with controls.

Each trace is counted **once**, at its median `tim`, so a long trace with many localizations cannot
dominate the early bins; the clustered-vs-background split is taken **within the same file**, so
acquisition length, laser power and sample differences cancel in the comparison. Normalized modes scale
each population to its own total, separating a difference in *shape* (timing) from one in *amount*
(abundance). Beyond the generic "correlational" caveat, three specific pitfalls:

- **The KS gap is a distance, not a p-value — and cannot become one.** The code reports only
  `max|F_clustered(t) − F_background(t)|` (no `ks.test`), which is correct: the two samples are **not
  independent**. They come from the same acquisition, from molecules competing for the same detection
  budget, so the KS null distribution does not hold and any p-value computed from it would be
  anti-conservative. Read the gap as an effect size only.
- **Crowding-vs-timing has a built-in positive trend, corrected by a permutation null.** Denser
  clusters have **more members**, and the median of a larger sample drawn from a decaying detection-rate
  distribution lands **later** purely by sampling, before any crowding mechanism. **Implemented null
  (`crowd_null`, default on):** detection times are shuffled among the traces **within each file**, while
  cluster membership, size and geometry stay fixed — so every point keeps its exact x coordinate and its
  exact n, and only the time↔geometry association is destroyed. Each structure therefore gets its own
  2.5–97.5% null interval (wide for small clusters, narrow for large ones), drawn as a shaded band.
  *Why a permutation and not the simulator:* the permutation keeps the real sizes, the real geometry and
  the real detection-time distribution of that acquisition and is assumption-free, whereas a simulator
  null would impose a bleaching/detection model that is itself calibrated from data of this kind — partly
  circular, and strictly more assumptions for a weaker answer.
  **Global test:** Spearman ρ between `nn_median` and `time_shift`, observed vs its permutation
  distribution (two-sided p = (1 + #{|ρ_b| ≥ |ρ_obs|})/(B+1), default B = 200). Rank-based, so it tests
  the **trend** and is unaffected by a constant clustered-vs-background offset — which the band alone is
  *not*: a global offset pushes every point off the band without any crowding dependence.
  *What it still cannot do:* rule out a third variable driving both axes (labelling density, structure
  identity). That needs a control experiment.
- **NN-vs-size is censored by the clustering algorithm.** Under DBSCAN every within-cluster NN distance is
  bounded above by the link radius `eps2` (default 200 nm), so the y-axis is **right-censored**: a curve
  that flattens at large cluster size may be hitting that ceiling.
  HDBSCAN/OPTICS impose no such bound — re-run under those before reading a slope. Convex-hull size is
  additionally upward-biased by single outlying members.

### I-A.10 Traces per cluster (blinking-event count) — [Part II §6.1a, 🟡 AMBER]
Counts **distinct TIDs** per de-blinked molecule or median cluster. Counting traces already
collapses the many localizations MINFLUX collects within one trace, so the number is a *blinking-event*
count. It is **not stoichiometry**: one fluorophore re-enters the on-state repeatedly and is re-assigned a
new TID each time, and a multivalent antibody carries several dyes — so the count over-estimates molecules
by an unknown, dye- and buffer-dependent factor. The Annibale et al. (2011) blink-correction that would
divide this out requires a **dark-time distribution that is not calibrated in this app**. Valid as a
relative comparison between conditions acquired identically; never as an absolute subunit count.

## I-B. Tracking (SPT) analysis

### I-B.0 Cycle time & Time-Signature Analysis (the foundation for everything Δt-aware)
The **Time Signature Analysis (TSA)** tab is an R port of the time-signature analysis in MinSpt
(*Toolbox_for_MINFLUX_SPT*; Vogler et al. 2025, *Commun. Biol.*). The correspondence:
- **cycle time / hardware parameters** — MinSpt `initialize_tsa` + `estimate_time_signal_parameter`
  (EOD/SW overhead, pattern vertices, dwell) ↔ `trk_cycle_time`, the EOD/SW regression and the hardware
  pattern (hexagon 6 / square 4 / triangle 3 vertices).
- **integrated vs dark cycles + per-track progression** — `print_cycle_stats` /
  `plot_cycle_progression_for_ID` ↔ the Integrated/Dark distribution and per-track cycle-progression plot.
- **adaptive gap removal** — `suggest_split_thresholds` (TOT-based `med_th = floor(median+std+1)`, capped
  by `high_th = maxOffTime/cycle`) + `remove_time_gaps` (`min_len = 50`) ↔ the adaptive gap-split.
- **linearized cycle-index timeline** — `prepare_for_msd` ↔ cumulative cycle-index binning.

MINFLUX has no fixed clock; the base **cycle time t_c** is inferred (`trk_cycle_time`), in priority order:
(1) the hardware model `t_c = N_iter · t_dwell` if enabled; (2) a manual entry; (3) the regression cycle,
if enabled; (4) the **cluster cycle** — cluster the inter-localization Δt values and take the spacing
between cluster centres; (5) the modal Δt as a last resort. The resulting t_c drives gap-splitting, the
quality gate, cycle-index MSD binning, the blur correction and the D-map. A step's **cycle count** is
`n_cyc = round((Δt − SW) / t_c)`; **integrated (single-cycle)** steps have n_cyc = 1.

**Δt is quantised.** MINFLUX intervals follow **Δt = k·t_c + SW** with no meaningful jitter within a
cycle: each peak of the Δt distribution is a near-delta spike. Two consequences follow.
- **The modal Δt is not the cycle time.** It is `1·t_c + SW`, and treating it as t_c biases every D
  low by the ratio SW/t_c. The cluster cycle therefore takes t_c as the **spacing between interval
  clusters** (`median(diff(centers))`) and SW as the leftover offset (`centers[1] − t_c`).
- **The TSA histogram has one bar per cycle count**, with edges on the classification boundaries
  `SW + (k ± ½)·t_c`, so the plot shows exactly the `round((Δt − SW)/t_c)` the analysis performs.
  Fixed-width bins would alias against the comb and misplace the peaks. No density curve is overlaid,
  because the width of a smoothed peak would be kernel bandwidth. A log count axis is
  offered, since counts fall steeply into the dark tail.

Because intervals are exact integer cycles plus a constant offset, cumulative cycle-index binning is not
an approximation to noisy timing: it recovers the true elapsed time.

### I-B.1 σ_loc & data-quality calibration (photon model + CVE) — [Part II §2.2]
`trk_cve_bundle`. Two estimators from two independent channels (full derivation in
[sigma_loc_model.md](sigma_loc_model.md)).

**App-wide σ_loc priority** (`trk_sigma_loc_eff`): manual entry > **photon-budget model** > blur-aware
CVE σ > MSD-fit intercept. The photon model comes first because it assumes nothing about the motion,
while the CVE σ is ~90% blur term near the full-cycle ceiling R = 1/6: it mostly reports √R × RMS step and inherits population
heterogeneity and drift.

- **Photon-budget σ_loc** (`.crlb_doughnut`, est_cve.R): per localization
  **σᵢ = L_TCP/(a·√Nᵢ) · √(1 + c·⟨r²⟩/L²)**, N = `eco`, ⟨r²⟩ = d·D_cve·t_loc the
  exposure-averaged offset of the moving emitter, reported as the median. c and d come from the **TCP
  pattern** (`TCP_PATTERN_GEOM`: hexagon c = 8.6, square 9.0, triangle 13.7, octahedron 9.0 with d = 3),
  derived from the Fisher matrix of an ideal doughnut zero and recomputed in the tests. **a** is a measured, effective value, taken from a calibration profile (`R/est_calib.R`; without one,
  the pattern's ideal value, labelled uncalibrated): the ideal unbiased per-axis bound is
  √8 = 2.83 for any planar pattern, but the instrument's estimator is biased toward the pattern centre
  and the loop correlates successive errors, so exported positions beat it and a must be calibrated per
  sequence on an immobile sample ([sigma_loc_model.md](sigma_loc_model.md)). The beam-size factor
  1/|1 − ln2·(L/ν)²| (ν, `TCP_BEAM_FWHM_NM`) is **off (0)** by default; enable it if an L-series on
  your system shows that σ does not scale with L.
- **CVE with motion blur** (Vestergaard 2014): **σ² = R·⟨Δx²⟩ + (2R − 1)·cov₁**, averaged over x, y[, z],
  using **only consecutive integrated (single-cycle) step pairs** so variable dark-cycle Δt cannot bias
  the covariance. R is bounded: `R_hi = min(1/6, t_c/(6·Δt))` (collection spans the
  cycle) and R from the photon-collection timeline of the sequence (§I-B.1a); `R_hi` is the upper bound. The instantaneous-exposure form (R = 0,
  σ² = −cov₁) under-reports σ severalfold on blurred data.
- **Free bias-corrected D:** `D_cve = ⟨Δ²⟩/(2Δt) + cov₁/Δt` per axis. **Blur-independent**, so D is
  unaffected by the choice of R.
- **Tracker wander floor D_w** (`trk_dw_floor`, manual, default off): the apparent D the tracking loop
  produces on an immobile target. Reported beside D as `D − D_w` with a warning when D_w > 20% / 50% of
  D; **never subtracted silently**. Measure it on an immobile sample at the same TCP diameter and photon
  budget; it scales steeply with L.
- **Per-track σ heterogeneity:** one CVE σ per track (≥ `n_min` integrated steps, skipping the first
  ~5 locs of lock-on transient); the median, IQR, CV and a **σ–EFO Spearman ρ** (precision should improve
  with photons ⇒ negative ρ validates the model) plus a bootstrap 95% CI on the median σ.
- **No photon weighting of D or σ.** Photon counts differ between tracks for several reasons besides
  precision (mobility, bleaching, emitter brightness, background), so weighting steps by them would
  weight tracks by whatever drives those differences. Every step counts equally.
- **No increment-correlation checks.** The lag-1 covariance of the steps is
  cov₁ = 2·D·R_abs − σ² (Vestergaard 2014): motion blur adds a positive term and noise a negative one. Its
  sign therefore cannot show smoothing (galvo damping), and −2·r₁ is no longer the noise share of the
  step variance once R > 0, so neither is reported; galvo damping is read from the sequence file.
- **Real z only.** A z column counts as real only with more than 10 finite values **and** sd > 1 nm
  (`.has_real_z(df, min_sd = 1)`), because 2D MINFLUX exports often carry an all-zero z column. The same
  predicate decides the dimensionality of the 3D toggle, CVE, MSD, JDD, confinement and
  displacement-vs-D, so the tools cannot disagree about it.
- **Per-axis σ: lateral vs axial.** The pooled σ_loc averages x, y **and z**, and in 3D MINFLUX the axial
  and lateral precision differ (the axial spread is often the tighter one). Several consumers are purely
  lateral — the **D-map** is an x–y map and the **HMM** never reads z — so they need the lateral value.
  `sigma_xy`, `sigma_z` and their ratio `sig_aniso` are reported separately; when the two differ by more than 1.3× either way the
  panel warns and gives the lateral value to enter manually for the 2D tools.
- **QC flags:** *implausible jumps* where `|Δr|² > 25·E[|Δr|²]` (i.e. |Δr| > 5× the expected RMS for that
  gap, using E[|Δr|²] = 2d·D_cve·Δt + 2d·σ² so long dark gaps aren't falsely flagged); *degenerate* tracks
  with σ² ≤ 0 (drift-dominated) or σ far in the outlier tail.
- **Continuity score:** tolerance-weighted `w = 2^{−(n_cyc−1)/H}` averaged over steps (H = dark-cycle
  tolerance), so a dataset of gentle 1–2-cycle gaps scores higher than its bare integrated fraction.
  A fast molecule collects dark cycles, so continuity can follow mobility; the *D by track continuity*
  table reports D for each continuity quarter.
- The app-wide σ_loc (priority above) is the **noise floor** consumed by the D-map, JDD, confinement,
  HMM, directionality, change-point and speed tools.

### I-B.1a Motion blur R — one value for every estimator (`R/est_blur.R`, `trk_blur()`)
A localization averages the molecule's path over the time its photons are collected, which lowers every
step's squared displacement by a constant: per axis **E[Δx²] = 2D(Δt − 2·R_abs) + 2σ²** (Berglund 2010;
Vestergaard, Blainey & Flyvbjerg 2014). Every estimator that reads single steps or an MSD intercept needs R.

- **What sets R.** R = (1/Δt)∫ S(t)(1 − S(t)) dt, with S(t) the cumulative weight the estimate gives the
  path. So R is set by **when** the photons are collected inside the step: dwell per vertex, number and
  order of vertices, pattern repeats, passes to reach the photon limit (η), beam moves between vertices
  (t_EOD ≈ 5 µs), the CFR centre dwell and re-centring (t_REP ≈ 20 µs) (Vogler et al. 2025, Eq. 4–7).
  **The TCP diameter L does not enter**: simulated with a MINFLUX estimator at L = 50, 101 and 302 nm and
  D = 0.2 and 2 µm²/s, R stayed at 0.084–0.095, inside the Monte-Carlo scatter.
- **The MINFLUX-specific part.** A vertex senses only the displacement along its own direction, so the
  linearised estimator weights vertex i by u²ᵢ,ₐ on axis a (hexagon x: 1, ¼, ¼, 1, ¼, ¼; y: 0, ¾, ¾, 0,
  ¾, ¾). R differs per axis and depends on the visiting order, taken as angular. The pooled estimators use
  the mean over x and y.
- **Inputs, best first:** the sequence JSON (pattern, `patDwellTime`, `patRepeat`, CFR from
  `ccrLimit` > 0 with `ctrDwellFactor`), with η = ECO/(EFO·dwell) from the data; otherwise the data alone
  (collection window ECO/EFO, the hardware-pattern vertex count). Δt is the mean single-cycle interval.
  Fallback without either: R = 0.10 per step (the standard 2D sequence). A sequence file is **not used**
  when its collection timeline alone is longer than the measured step (for example an imaging sequence
  loaded with tracking data: 2175 µs against a 205 µs step) or when its final pattern is axial-only
  (`zline*`); the panel says which, and the data route is used. An unknown pattern name is treated as a
  regular polygon with the parsed vertex count.
- **Conversions** (`blur_offsets`): R_abs = R·Δt₁ (ms; each step's effective diffusion time is Δt − 2R_abs,
  JDD / D-map / change-point / HMM), R_msd = R_abs/t_c (MSD / OLSF and the confinement null).
- **Validation** (`tests/testthat/test-blur.R`; standard 2D tracking sequence: hexagon, 100 µs dwell, CFR, Δt 203 µs):

  | R per step | value | single-step D on simulated MINFLUX localizations |
  |---|---|---|
  | full simulation (photon counts + linearised estimator) | x 0.096, y 0.097 | truth |
  | timing model, projection weights (used) | x 0.103, y 0.097 → 0.100 | within 2% |
  | timing model, equal vertex weights | 0.106 | +2–3% |
  | MinSpt constant 1/6.2 of t_c | 0.152 | **+15–17%** |

  R moves with the sequence (fast triangle 50 µs: 0.067; two pattern repeats: 0.129; two passes: 0.120),
  so no constant stands in for it.
- **Deliberate deviation from MinSpt.** MinSpt passes R = 1/6.2, cited in its code to Schmidt et al. 2020
  with a question mark; Schmidt et al. 2021 (Nat. Commun. 12:1478) mention blur only qualitatively and give
  no coefficient, and no sequence field sets one. Until 2026-09-28 the app used it for the MSD, JDD, D-map,
  change-point and HMM.
- **Damping is not blur.** A sequence with `damping` > 0 smooths the exported positions **across** steps
  (the increment-consistency check in §I-B.1 detects it). No R describes that; the Data Quality panel
  warns when the sequence says damping > 0, and D from single steps is then approximate.

### I-B.2 Unified MSD, optimal-fit OLSF & fit models — [Part II §2.1, §2.7, §3a]
`msd_unified_bundle` + `msd_unified_fit` + `.msd_model_fit`.
- **Binning:** per-track **cumulative cycle-index** lag (the MinSpt timeline) — `dcyc = round((Δt − sw)/t_c)`,
  `cycidx = cumsum(dcyc)` — so pairs are grouped by *integer elapsed cycles* (true elapsed time), removing
  the `round(total/t_c)` sawtooth. Per-track TAMSD = Σdr²/Σpairs; x-axis τ = mean real Δt per bin.
- **Ensemble:** **pair-weighted** (default, smoother) or **track-weighted** eMSD = mean of per-track TAMSDs
  (Vogler et al. 2025). Optional **D>0 filter** drops tracks whose TAMSD has a negative OLS slope; IQR
  band from per-track quantiles; bins kept only where ≥10% of tracks contribute.
- **Linear fit — OLSF (Michalet & Berglund optimal point count):** unweighted linear fit
  `MSD = β₀ + β₁·τ`, `D = β₁/(2d)`, with **separate optimal point counts** for the slope (→D) and the
  intercept (→σ²). With the reduced localization error `x = σ²/(D·t_c) − 2R_msd` and M the number of
  positions in the track, the optimal number of points is `p = f·ℓ/(f³+ℓ³)^{1/3}`, using
  `f_b = 2 + 1.35·x^{0.6}, ℓ_b = 0.8 + 0.564·M` for the slope and `f_a = 2 + 1.6·x^{0.51},
  ℓ_a = 3 + (4.5·M^{0.4} − 8.5)^{1.2}` for the intercept (Michalet & Berglund 2012). Two variants: **anchored** (default — anchors x
  on a stable whole-window fit so p is reproducible) and **iterative** (the MinSpt loop — converges onto
  the lowest, least-noisy lags). Fit window default **5–50 cycles** (Vogler `fit_limits`).
- **Motion blur in the OLSF.** The intercept of a blurred MSD is `β₀ = 2d(σ² − 2D·R_abs)` (Michalet &
  Berglund 2012), so the fits report `σ² = β₀/(2d) + 2D·R_abs`, and R_msd = R_abs/t_c enters the reduced
  error x (§I-B.1a). The slope, and so D, is unaffected; reading the intercept
  with R = 0 under-reports σ (by ~20% at the physical R). The drawn fit line uses
  the raw fitted intercept, because the data sit on the blurred line. The non-linear models (confined,
  directed, anomalous) read σ from their offset C without the blur term; their offset is rarely
  separable anyway.
- **Non-linear models** (whole window — they need curvature):
  - **Confined:** `MSD = A(1 − e^{−τ/τ_c}) + C`, with plateau A, relaxation time τ_c and offset C →
    local `D = A/(2d·τ_c)` and confinement size **L_c = √(6A/d)**; flagged *degenerate* if
    τ_c > 2·max(τ) (no plateau resolved).
    **L_c is fitted from the data, never entered** — it is read off the plateau A, so it is only as
    meaningful as the plateau is resolved (hence the *degenerate* flag) and is additionally warned
    about when L_c < 2σ_loc.
    **Why √(6A/d):** a particle uniformly exploring a d-dimensional box of side L_c has
    Var(Δx) = L_c²/6 per axis, so A = d·L_c²/6. This reduces to √(3A) in 2D.
    **Symbol note:** the confinement size is written **L_c** and the MINFLUX targeted-coordinate-pattern
    diameter **L_TCP** (D_max panel).
  - **Directed:** `MSD = 2d·D·τ + (v·τ)² + C` → D and drift speed **v**.
  - **Anomalous:** `MSD = 2d·D·τ^α + C` → exponent **α** (fit on the *ensemble* MSD, with a parametric
    α SE from the fit covariance, labelled a lower bound since MSD lags are correlated).
- **Per-track outputs** (`trk_unified_tdiff`): per-track OLSF D (→ tDIFF mean + 95% CI), and a
  **confinement index CI = (MSD_late/MSD_early)/(τ_late/τ_early)** from the robust median of the first/last
  15% of windowed bins: **CI < 1 confined, ≈1 Brownian, >1 directed.** These feed the D/CI **k-means
  clustering** (feature = z-scored α doubled + log₁₀CI; clusters sorted slow→fast).

### I-B.3 Change-point / motion-state segmentation (PELT + GMM) — [Part II §2.8]
`cpd_segments_rv`. Per track, segment the series **log(r²/Δt)** (= local log-D, already Δt-aware) with
**PELT** (`changepoint::cpt.mean`, BIC penalty, `minseglen`). Per-segment **D = mean(r²/Δt)/(2d)** — the
mean-of-ratios form avoids the Jensen bias of `mean(r²)/(2d·mean(Δt))` when r² and Δt covary. An optional
**σ_loc noise-floor subtraction** (default on) removes the immobile offset **σ_loc²·⟨1/Δt⟩** at the D
level (clamped ≥0 once per segment, since clamping each step would bias it upward), the same floor the D-map uses, so
truly immobile segments read ~0 instead of a small non-zero D. Motion blur is corrected in the segment D
in the same way as in the JDD (§I-B.4). The pooled
per-segment log₁₀D distribution is then fit with a **BIC-selected Gaussian mixture** (`mclust`, up to
`max_g` states).
- **Variance scaling.** `cpt.mean`'s Normal cost assumes unit variance, while the log of a scaled χ²_d
  has fixed variance trigamma(d/2) (π²/6 ≈ 1.64 in 2D, 0.93 in 3D). The signal is divided by that SD
  before PELT; unscaled, PELT cuts a large share of single-population 2D tracks.
- **The GMM component count is confounded by segment length.** The mixture treats every per-segment
  log₁₀D as equally precise, so one population split into segments of different lengths can return
  G = 2. A length-aware model (known per-segment variance ≈ trigamma(d/2)/n in log space) would remove
  this; until then the state count is 🟠 exploratory.

### I-B.4 Jump-distance distribution (JDD) — [Part II §2.3]
`R/est_jdd.R` (`jdd_collect_jumps`, `jdd_fit`), called by `jdd_bundle`. Fit the step-length histogram with an **N-component scaled-chi mixture** — the law of
|Δr| for a d-dimensional isotropic Gaussian displacement:

**p(r) = Σ_k w_k · C_d · r^(d−1)/σ_k^d · exp(−r²/2σ_k²)**,  **C_d = 1/(2^(d/2−1)·Γ(d/2))**

with w_k the weight and σ_k the per-axis spread of component k, and Γ the gamma function, fitted by
maximum likelihood (Nelder–Mead). d = 2 gives the **Rayleigh** form (C_d = 1, Schütz 1997) and d = 3
the **Maxwell** form (C_d = √(2/π)), used when z is included. MINFLUX handling:
- **√Δt normalization ON by default** — pool `r/√Δt` so steps of different cycle counts are comparable
  (so states map to D).
- **D per component:** `D = σ²/2`, **dimension-independent** — because E[r²] = d·σ² and MSD = 2·d·D·Δt,
  the d cancels. This also makes the σ_loc term below dimensionally consistent with it.
- **Noise and blur correction:** with √Δt normalisation a component's per-axis variance is
  **σ² = 2D·(1 − 2R_abs·⟨1/Δt⟩) + 2σ_loc²·⟨1/Δt⟩**, so **D = (σ² − 2σ_loc²⟨1/Δt⟩) / (2·(1 − 2R_abs·⟨1/Δt⟩))**,
  clamped ≥ 0, R from §I-B.1a. Both corrections use the mean of 1/Δt over the jumps in the fit, so
  pooled steps of mixed cycle count are handled exactly. Without the blur term D reads about 20% low at
  MINFLUX tracking timings; with the old constant 1/6.2 it read 15–17% high. With one component the MLE is closed form (σ² = ⟨r²⟩/d). A **noise-floor jump
  distance** `√(2d)·σ_loc·√⟨1/Δt⟩` marks where peaks are noise-dominated.
- **2D and 3D.** The d-dimensional law is used when z is included, so the fitted σ is per axis in both
  cases and the σ_loc term subtracts the right noise floor. The **HMM** and **D-map** read x and y only,
  so the 2D form is correct there, and since the 2D projection of isotropic 3D diffusion has the same D,
  their D is unbiased on 3D data.

### I-B.5 Diffusive-state HMM (continuous-time, maximum likelihood) — [Part II §2.4]
`R/est_hmm.R` (pure functions), called by `hmm_bundle`.

- **Model.** Hidden state per step, a continuous-time Markov chain with generator Q. A step of any Δt has
  transition matrix **P(Δt) = expm(Q·Δt)** (scaling-and-squaring Taylor, `.expm_small`), so dark cycles are
  exact and the **rates are parameters of the likelihood**. Emission: squared step length over x, y,
  **r² ~ Gamma(d/2, 2v)** with per-axis variance **v = 2·D_s·(Δt − 2R_abs) + σ_j² + σ_{j+1}²**
  (diffusion, motion blur with R from §I-B.1a, and the noise of both localizations). With the
  photon model active each localization carries its own σ from its photon count, which the panel derives
  itself if the calibration has not been run. D is therefore noise- and blur-corrected by construction.
  Initial distribution = stationary distribution of Q.
- **Fitting.** All tracks at once. Tracks are cut into 256-step chunks, each starting from the stationary
  distribution (composite likelihood; only transitions across chunk edges are lost), sorted longest first so
  the scaled forward recursion runs vectorised over chunks. P(Δt) is computed once per unique Δt, which is
  cheap because MINFLUX Δt is quantised. Parameters on log scale (D ordered slow → fast by cumulative
  increments), BFGS from three rate starts (dwell ≈ 5, 30, 300 cycles), best log-likelihood kept.
- **Step budget.** Fit time scales with the number of steps. Above the budget (default 150,000), whole
  tracks are drawn at random until it is met, and the summary says so. The precision of D is limited by how
  much molecules differ, so the extra steps would cost time without narrowing the interval.
- **Uncertainty.** Parametric draws from the inverse Hessian give CIs for D, rates and dwell times; these
  assume one homogeneous population and are too narrow on real data. `hmm_ct_bootstrap` refits on tracks
  resampled with replacement (the intervals to report).
- **Model size.** Optionally refits K − 1 and K + 1 and reports BIC differences (two further fits).
- **Decoding.** Viterbi under the fitted continuous-time model, used only to colour tracks and feed the
  per-state MSD panel. Rates are never counted from the decoded path, which reads dwell times long.
- **Stated approximations.** (1) one state per step, so a switch inside a step is not modelled; (2) steps
  sharing a localization have correlated noise (−σ²), treated as independent; (3) the chunking above.
- **Usable range.** On simulated two-state tracks at MINFLUX tracking timings, D and dwell times are
  recovered when the mean dwell is **≥ 10 cycles**. At 1–3 cycles the states mix, dwell times read long
  and the likelihood CIs do not cover the truth (`tests/validation/results/`). States whose fitted dwell is
  under `HMM_RESOLV_MIN_DWELL_CYC` (10) cycles, or with fewer than `HMM_RESOLV_MIN_SWITCHES` (300) expected
  switches, are flagged, as are two states with the same D. Not covered by that simulation: asymmetric
  dwells, K > 2, heterogeneous D within a state, subdiffusion and galvo damping.
- **Reading a flagged fit.** A K-state HMM always returns K states. When every state is flagged, the fit
  is describing a spread of step sizes rather than persistent switching. An estimate that moves a lot
  when the step budget is lowered points the same way.
- **Ground truth.** `sim_ct_state_tracks` switches states in continuous time (including inside a step),
  averages the path over each collection window (blur), draws dark cycles and per-localization σ. The
  in-app simulator's two-state mode has a slow-state D fraction, separate switch probabilities in each
  direction, and exports `state_true`.

### I-B.6 Confinement (MSD confinement index + simulated null) — [Part II §2.5]
`conf_ci_bundle` (server) + `sim_free_ci_null` (`helpers_tracking.R`). Per track, from that track's own
cycle-binned MSD (`trk_unified_tdiff`):

**CI = (MSD_late / MSD_early) / (τ_late / τ_early)**, early/late = robust medians of the first/last 15% of
the windowed lag bins. Ideal free motion gives CI = 1, confinement < 1, directed motion > 1.

- **Why this statistic.** It compares a track's own short lags with its own long lags, so it needs no
  external D and is not biased by D varying between molecules (the confound of the displacement test).
- **Why the null must be simulated.** The observed free CI is **not** 1: the 2d·σ² intercept lifts the
  early lags more than the late ones. `sim_free_ci_null` simulates free tracks with the data's own
  t_c, σ_loc, measured cycle-gap distribution and fit window, binned exactly as the app bins. Each track
  is compared with the null **at its own localization count**, because the free range moves strongly
  with track length, and, when noise matters (x = σ²/(D·t_c) > 0.3 somewhere and the D spread is ≥ 2×),
  **at its own D**, so the slow half of a heterogeneous population is not judged against a null built for
  the fast half. Null tracks are blur-averaged over a collection window τ = 6R·t_c with R from §I-B.1a, and
  when σ_loc comes from the photon model each null localization draws its σ from the measured
  per-localization values (the MSD offset is set by mean σ², which a long-tailed photon count pushes
  above median σ²).
- **Classes:** confined = CI below the matched free 5th percentile; directed = above the 95th; else
  free-like. This is the primary call on confinement in the app.
- **Sensitivity.** CI detects confinement well within the free extent of a track, and is weak for
  directed transport by design; use Directionality or Speed & Transport for transport.
- **Resolution floor:** the free extent √(2·d·D·T) is reported. Confinement is detectable only when that
  extent exceeds the domain roughly twofold, so a low confined % means "none below this scale", and
  **confined % is not comparable between datasets with different track lengths or mobilities**.

#### Effective dimensionality of the MSD
`MSD = 2·d_eff·D·τ`, where **d_eff follows the MOTION and its projection, not the number of
recorded coordinates** (`.msd_d_eff`, threaded into `msd_fit_olsf`, `msd_fit_olsf_iterative` and
`.msd_model_fit` via `coef_override`; `ndim` still carries the geometric dimension because Michalet's
optimal point count is derived for real spatial dimensions).

| motion | tracked | d_eff | note |
|---|---|---|---|
| volume | 3D | 3 | the plain case |
| volume | 2D | 2 | dropping z leaves two independent axes |
| surface | 3D | 2 | 2 DOF however many coordinates are recorded |
| surface | 2D | 1 + ⟨n_z²⟩ | shape-dependent — see below |

**The cell is three-dimensional even when tracking is not.** A molecule on a curved surface steps in its
local tangent plane, and how much of that step survives projection into the image plane depends on how
the surface is tilted there. Averaging over the shape gives ⟨n_z²⟩ = 1/3 for a sphere, 1/2 for a
cylinder whose axis lies in the image plane, and 1 for a flat membrane parallel to it — hence
**d_eff = 4/3 (sphere), 3/2 (rod), 2 (flat membrane)**, i.e. MSD = 8Dτ/3, 3Dτ and 4Dτ. A membrane
tilted by θ gives 1 + cos²θ, between the rod and flat values; not offered separately. Dividing by the
number of recorded coordinates instead reads D 25–30% low for surface motion.

Exposed as **Molecule moves** (volume / surface) plus, when surface motion is recorded in 2D, a
**Cell shape** selector. Defaults to volume. The estimates line states the relation in force whenever
d_eff ≠ ndim. Applied to the ensemble fit, the per-track fits, the per-cluster models, the bootstrap and
the per-segment (HMM state) MSD.

#### Free-diffusion-on-a-geometry MSD reference
`sim_free_geom_msd` + `msd_geom_null`. On a **closed** cell the MSD of a **free** molecule bends and
eventually saturates on its own, so a bend is not by itself evidence of trapping. This simulates free
diffusion on the chosen shape at the fitted D — same lag grid, same observed coordinates — and overlays
the curve it produces. Measured curve *on* the reference ⇒ the cell explains the bend; *below* it ⇒
confinement beyond the cell.

**The result depends on the fit window.** Over a short window the molecule covers far less than the cell
and the two are indistinguishable; the gap opens as the upper fit limit is raised. The estimates line
reports how much of D the geometry alone costs over the window in use, rounded to 5% because it is a
Monte-Carlo estimate, and warns below 95%.

**The shape enters the MSD in two separate places:**
1. **The divisor**, at every lag, independent of size — the d_eff table above.
2. **The ceiling**, only once the molecule has explored the shape. The MSD saturates at the mean
   squared separation of two random points on it: **sphere 2R² (3D), 4R²/3 (projected); rod
   L²/6 + 2r² (3D), L²/6 + r² (projected)**. Every term is a **length squared** — the volume
   (4/3·πR³) never appears. Size sets where the curve bends over.

**The ceiling fakes biology in every non-linear model.** Free diffusion on a 500 nm sphere (2D, surface)
fitted over a 4–2000 cycle window returns:

| model | what pure geometry gives |
|---|---|
| linear | ~95% of the true D |
| confined | **L_c ≈ 1000 nm** — a converged, entirely spurious corral |
| anomalous | **α ≈ 0.86** — sub-diffusive with nothing anomalous present |
| directed | v ≈ 0 — geometry does not fake drift |

With the overlay enabled, the estimates line fits **whichever model is selected** to the simulated free
curve and reports what the shape alone would return, next to the measured value, warning when the shape
already accounts for most of it.

**Where the geometry is set.** One shared **Sample geometry** block at the top of the Analysis tab
(`trk_geom_type` plus the dimension inputs) is read by the MSD projection factor, the MSD overlay and the
displacement Step-2 null, so the sample cannot be described two different ways. Options are None /
sphere / rod / ellipsoid / flat membrane / **filament**, and *Molecule moves* carries the matching
1/2/3-DOF choice (a filament covers in-vitro motor assays). `None` and `filament` have no closed surface
to simulate on and map to the unbounded case, so Step 2 reduces to the Step-1 band. The D fed in is the
d_eff-corrected D, a true motion D. Cost scales with the largest lag.

**Curvature vs the DOF term.** With d_eff handled, curvature costs only a few per cent of D at the
default short window; the DOF/projection factor is the larger error. Geometry matters mainly when the fit
window is long, which is when this overlay should be enabled.

### I-B.6a Displacement vs diffusion (net displacement test)
`trk_disp_bundle`. Per track, net start→end displacement `R_net = |r(T) − r(0)|` against the free
expectation. Depends only on the two endpoints and the true elapsed time T, so it is **inherently robust
to irregular Δt** — no cycle binning needed. Endpoint noise `2·d·σ_loc²` is subtracted.

**Motion vs observation degrees of freedom.** The band is `R_net² = 2·D_mot·T·χ²(d_mot)`, where
**d_mot is the degrees of freedom of the MOTION, not the number of recorded coordinates** (`d_obs`).
They differ for a molecule confined to a *surface* (membrane, cell envelope) but tracked in 3D: 2 DOF
observed on 3 axes. Two consequences:
1. A per-axis estimator (CVE, or MSD with `ndim = 3`) spreads 2-DOF motion over 3 axes and reports
   `D_meas = D_true · d_mot/d_obs`, exactly **2/3** for a surface in 3D. The panel rescales it:
   `D_mot = D_meas · d_obs/d_mot`. A **manual** D is taken as the true motion D and passed through
   unchanged.
2. The χ² degrees of freedom follow the motion, so the band uses `χ²(d_mot)`.

Correcting (1) without (2) leaves the *mean* right and the *band shape* wrong, so free surface-bound
molecules are called sub-diffusive at about twice the nominal rate. When `d_mot == d_obs` (all 2D data,
and true volume diffusion in 3D) both corrections vanish. The choice is exposed as
**Motion type: volume (3 DOF) / surface (2 DOF)**, shown only for 3D data, defaulting to volume.

**Relation to the confinement index.** The confinement index (§I-B.6) is the primary call on
confinement: its null is matched to each track's own length and D, while this band uses an infinite
plane and one ensemble D. This panel is read as corroboration and as an upper bound, and reports how far
the two agree track by track.

**Step 2 (geometry null)** already handled this independently: `sim_free_geom_rnet` has a
`mode = surface | volume` switch and empirically calibrates the projection factor (`f = 2/3` for a curved
surface tracked in 3D, `3/4` for a 2D projection of a cylinder), so it is fed the **raw per-axis D**, not
the rescaled `D_mot` — passing the rescaled value would double-correct.

### I-B.7 Directionality / turn-angle autocorrelation — [Part II §2.6]
`dir_bundle`. **Turn angle** between consecutive displacement vectors θ = acos( (Δr_i·Δr_{i+1}) /
(|Δr_i||Δr_{i+1}|) ), clamped to [−1,1]. **Velocity autocorrelation** of unit step vectors,
`AC(lag) = ⟨v̂_i · v̂_{i+lag}⟩`. Summaries: median angle, %forward (<60°), %backward (>120°), AC(first
lag). A **noise floor** RMS pure-noise step = 2·σ_loc; the fraction of steps below it quantifies how
much of the backward-peak / first-lag anti-correlation is localization noise rather than real
anti-persistence (MINFLUX's small σ_loc largely removes this classic artefact).
**Irregular-Δt handling:** both quantities use only **integrated-ish steps** — those spanning
≤ `input$trk_dir_max_dark` dark cycles (default 2; `n_dark = round(Δt/t_cycle) − 1`); longer-gap steps
are dropped (unreliable direction, possible cross-molecule merge). The autocorrelation lag is
**elapsed time**, binned by integer cumulative-cycle index (as in the MSD), so
the x-axis is real time (cycle lag × t_cycle) and unevenly-spaced points are not misrepresented as
evenly spaced. Step-index fallback only when t_cycle is unknown.

### I-B.8 Speed & transport (path speed, drift velocity, local α) — [Part II — MINFLUXR original]
`trk_export_summary`. Three complementary per-track curves:
- **Path speed:** `dr/dt` per step (µm/s), lightly median-smoothed. Cannot distinguish diffusion from
  transport (diffusion inflates instantaneous speed).
- **Drift (transport) velocity:** `v_drift(t) = |net displacement over a time window W| / W`. Diffusion's
  net displacement grows only as √t, so its apparent v ~ √(4D/W) *shrinks* with W, while directed transport
  gives a sustained v_drift — this is what actually detects export. Windows are defined by **time (ms)**
  (robust to irregular Δt) with a min-points guard and an **end-guard** (`min_span_frac`) that drops
  end-of-track windows whose duration collapses (they otherwise fake a terminal velocity spike).
- **Local α(t):** rolling log-log slope of a short-window MSD (τ = *real* mean elapsed time per step-lag,
  not step index), with a static-noise offset `2d·σ²` subtracted per lag. α≈1 Brownian, α→2 ballistic.
- **Export call:** a window is "directed" if v_drift exceeds `max(fold × diffusive baseline, √2·σ/W)` —
  fold over the track's own baseline (≈5% false-positive on free tracks at fold ≈ 3), floored at the
  localization-noise velocity.

### I-B.9 Spatially-resolved D-map — [single-lag = short-lag CVE limit]
`R/est_dmap.R` (`dmap_compute`), called by `dmap_bundle`. Single-lag estimator per step: **D_raw = r²/(4·Δt)**
(2D), midpoint-binned on an n×n grid; per-cell
**D = max(0, mean(D_raw) − σ_loc²·mean(1/Δt)) / (1 − 2R_abs·mean(1/Δt))** with its standard error (same
scaling), cells with
< `min_pts` steps left NA. Options: **integrated-only** (drop dark-gap steps whose midpoint is spatially
meaningless) and **noise subtraction** (removes the spurious immobile-region diffusion floor), applied once
per cell. The cell value is the **mean** of per-step D, not the median: a single 2D step's r² is
exponential, so its median is ln 2 ≈ 0.69 of the mean. A mis-linked jump therefore shows up as a large
per-cell SE. Without the blur factor the map reads about 30% low at MINFLUX
cycle times. Far more robust than a per-cell MSD fit in small cells, where no τ-bin has enough steps.
This is the short-lag limit of the CVE (Vestergaard 2014; Vogler et al. 2025).

### I-B.10 Fast-diffusion ceiling (D_max)
`D_max = (L²/16)·ψ / (N_PL + ψ·t_hw)` (Vogler 2025), with ψ the photon count rate (median EFO), N_PL the
photons per localization (median ECO) and t_hw the hardware overhead per localization (default 0.05 ms):
the trackable-D limit before the localization loop loses lock; surfaced as a green/amber/red flag. Theoretical form, so it may under-warn near the
demonstrated ~2.5–3 µm²/s practical ceiling (Part II §3a item 8).

## I-C. Simulation

The two simulators use **different** precision constants, deliberately:
- **Diffusion / tracking simulator:** **σ = L_TCP / (a·√N)** with `TCP_CRLB_PREFACTOR` (the 2D ideal-estimator value, √8)
  as the Tracking tab's σ_loc model, so the forward and inverse models cannot disagree. σ varies per
  localization with the drawn photon count N.
- **SMLM simulator (CRLB mode):** **σ = √[(k·σ_CRLB)² + σ₀²]**, σ_CRLB = L/(√8·√N)·g(SBR), the ideal
  per-axis bound for a 2D doughnut pattern (Balzarotti 2017 Suppl. S27/S31; Marin & Ries 2025 Eq. 44).
  L is the physical final-iteration diameter (~40 nm), k scales the photon-limited term and σ₀ is a
  photon-independent floor (drift, vibration); both are fitted from a real run as the slope and intercept
  of within-trace variance against 1/N. One σ is used for every localization. See [simulator.md](simulator.md) §1.4.

### I-C.1 SMLM structure + labelling simulator (`sim_server.R`, SMLM tab)
Generates a *known* super-resolution scene to validate the SMLM estimators and the GT inverse.
- **Structures:** sphere/shell, ring, double-ring (Z- or Y-separated), stacked rings (helix-like),
  filament (random-walk backbone), each with a chosen point count and geometry.
- **Three-level labelling model** (shared with the GT inverse via `helpers_shared.R`):
  1. *molecule* at the structure point;
  2. *fluorophores* — count `K ~ Poisson(n_fluor)` (min 1), the shared complex displaced from the epitope
     by a **half-normal linker** (σ = link_nm·√(π/2) ⇒ mean ≈ link_nm), optionally **direction-aware**
     (outward from the local structure normal, since an antibody binds a surface epitope not the centre);
  3. *dyes* — each fluorophore scattered by `dye_spread_nm` (Gaussian).
- **Blinking / trace length:** each fluorophore emits `n_loc` localizations, drawn **Poisson(trace_length)**
  or **Geometric(1/λ)** (same mean, heavy right tail — constant per-step bleaching gives many short traces);
  extra blink re-appearances via `blink_lambda`. Each localization gets Gaussian noise σ_loc.
- **Acquisition model** (Acquisition mode): detection **rate(t) = structure·(floor + (1−floor)·e^{−t/τ})·
  bleach(t) + background_floor**, with photobleaching `bleach = e^{−cumulative/cap}` consuming a photon
  budget; τ and floor can be calibrated from a real run. A DBSCAN
  core-point test (≥minPts within eps=60 nm) splits recorded traces into **structure vs background** so the
  measured rate curves can be overlaid on the model.
- **Output:** MINFLUX-format `data.frame(tid, tim[s], x, y[,z], efo, eco)` + a molecule/fluorophore
  ground-truth sidecar — loads straight into the SMLM tab so the GT inverse can be scored against truth.

### I-C.2 Diffusion / tracking ground-truth simulator (`sim_server.R`, Diffusion tab)
Simulates a *known* diffusion process on a *known* geometry, samples it with **MINFLUX's own pathologies**,
and emits MINFLUX-format tracks so they flow through the *exact same* Tracking analysis — closing the
validation loop (see [simulation_module.md](simulation_module.md) for the full description).

**Geometry operators** (unified `project` / `reflect` / `normal`):
- **Rod** = spherocylinder (cylinder body radius r, length L_cyl + two hemispherical caps): `project` =
  nearest surface point, `reflect` = fold back inside, `normal` = outward unit normal.
- **Sphere / ellipsoid** (a=b=c for sphere): radial project/reflect in scaled coordinates.
- **Membrane** = flat sheet with thickness (2D lateral + z wobble within thickness); **Free** = unbounded 3D.

**Continuous path** (`gen_path`, per cycle): per-axis step SD **sdi = √(2·D·t_c)**.
- *Surface* diffusion: tangent-plane step `g − (g·n̂)n̂` then `project` back onto the surface.
- *Volume* diffusion: 3D step then `reflect` at the wall. *Free/membrane*: unbounded / 2D-lateral.

**Motion type** (each a ground-truth label):
- **Free** Brownian.
- **Anomalous:** fractional Brownian motion, MSD ∝ t^α, with D set at the single-cycle lag.
- **Hop diffusion:** free diffusion inside corrals with permeable barriers.
- **Confined = OU tether** to a per-track home: `p ← home + (p − home)·φ`, **φ = e^{−t_c/τ_c}**,
  **τ_c = L_c²/D**, giving a stationary SD of L_c per DOF and an MSD plateau of 2d·L_c².
- **Directed:** add drift `v·t_c` along the rod long axis (or a random direction otherwise).
- **Two-state:** continuous-time switching between a mobile and a slow state, with a slow-state D
  fraction and separate switch probabilities in each direction; the ground truth for the HMM.

**MINFLUX sampling overlay:**
1. **Cycle clock** t_c; a localization is attempted each cycle.
2. **Dark cycles → irregular Δt:** each cycle recorded with prob `1 − p_dark`, else skipped, giving
   variable-cycle gaps.
3. **Photobleaching → short tracks:** track length in cycles is **Geometric**, with the mean scaled so the
   *recorded median* (after dark-cycling) lands on the target `mean_len` (median ≈ ln2·mean, only
   (1−p_dark) recorded ⇒ mean_cycles = median/(ln2·(1−p_dark))).
4. **Heteroscedastic σ_loc:** per recorded loc draw **N ~ Poisson(photons)**, σ = L_TCP/(a√N) with a = `TCP_CRLB_PREFACTOR`, add Gaussian
   noise per axis; record `eco = N`, `efo = N/(t_c in s)`.
5. **Survivorship (optional):** if the true displacement across a dark gap exceeds a re-acquisition radius
   `R_reacq`, the track terminates there — reproduces the dark-cycle survivorship bias (long jumps culled →
   *manufactured* apparent confinement).
6. **Projection:** output 2D (drop z) or 3D — a surface protein on a curved cell shows the 2D-projection
   distortion the simulator exists to study.

**Built-in CVE round-trip** (`sim_trk_cve`) scores the recovery immediately:
**D = (⟨|Δr|²⟩ + 2⟨Δr_i·Δr_{i+1}⟩) / (2·d·Δt₁)** over consecutive single-cycle steps, using the **total**
displacement with **d = min(motion-DOF, observed axes)** — surface diffusion has 2 DOF but spreads over 3
Cartesian axes (a per-axis CVE would under-estimate D by ⅓); volume has 3; 2D projection caps d at 2.
Flat geometries recover D; a curved surface projected to 2D does not, and the summary reports that
ratio as a geometric-bias factor.

**Key resolvability caveat:** confinement is only visible when **τ_c = L_c²/D is shorter than the track
duration**; otherwise free diffusion over a short track mimics the same ~extent — the whole reason the
round-trip (not the raw MSD/extent) is the test for a candidate confined protein.

---

# Part II — Literature & Appropriateness

## 1. Why MINFLUX data needs its own treatment

MINFLUX localizes an emitter by probing it with a doughnut-shaped excitation minimum at several
positions (a targeted coordinate pattern) and estimating position from the relative photon counts
near the zero-intensity point. This gives a data structure that differs from camera-based SMLM in
ways that matter for every downstream estimator:

| Property | Value / behaviour | Consequence for analysis | Source |
|---|---|---|---|
| **Localization precision σ_loc** | ~1–3 nm in cells; <1 nm (focal plane) / ~2 nm (3D) on origami | Localization-error terms shrink; precision co-estimation is less valuable | Balzarotti 2017; Gwosch 2020; Ostersehlt 2022 |
| **Cycle time** (one localization) | ~0.4 ms–a few ms; set by the photon budget | No single nominal Δt | Balzarotti 2017; Wolff 2023; Deguchi 2023 |
| **Irregular Δt** | consecutive locs separated by integer multiples of a base cycle (skipped/"dark" cycles) | **Estimators assuming constant Δt are biased** | Scheiderer/Marin/Ries 2025; Marin & Ries 2025 |
| **Gaps / dark cycles** | photobleaching budget + µs–ms blinking abort/skip cycles — *structural* | Must be modelled | Marin & Ries 2025 (SimuFLUX) |
| **Motion blur (R factor)** | ≈ 0.07–0.13 per step, set by the sequence timing (§I-B.1a) | σ from a covariance or MSD intercept is under-read unless R enters; D is unaffected | Michalet & Berglund 2012; Vestergaard 2014; MinSpt uses R = 1/6.2 |
| **Track length N** | photobleaching-limited; ~100× more locs/track than camera-SPT, but still short (tens–low hundreds) | Short-track estimator bias dominates | Balzarotti 2017; Deguchi 2023 |
| **Trackable D ceiling** | ≈ 3 µm²/s (standard settings) before the loop loses lock | Fast diffusion is systematically under-estimated | Marin & Ries 2025; Vogler 2025 |

**The two properties that drive almost every rating below are (a) irregular Δt and (b) short tracks.**

---

## 2. Methods implemented — sources & appropriateness

### 2.1 MSD with optimal point-count fitting — 🟡 AMBER
- **Implemented as:** unweighted OLSF `MSD = β₀ + β₁·τ` with Michalet–Berglund optimal number of fit
  points (`p_a`/`p_b`), R from the sequence timing (§I-B.1a; MinSpt uses 1/6.2), plus a **cycle-based MSD** that rounds each lag to the
  nearest integer number of acquisition cycles.
- **Source:** Michalet 2010; Michalet & Berglund 2012 (optimal least-squares fitting / optimal `p*`).
  Cycle-based binning and the trackable-D ceiling: **Vogler et al. 2025**.
- **MINFLUX caveat:** the optimal-`p*` formalism is *built for* short, noisy tracks (a match), and
  the blur term R is part of the same theory — but the framework assumes a single
  nominal Δt for lag binning. Binning by *step index* on irregular-Δt data corrupts the MSD curve.
- **MINFLUXR's handling:** the cycle-based MSD (Vogler 2025) bins by **integer cycle count** = true elapsed
  time, which is exactly the recommended fix. ✅ This lifts the practical rating toward GREEN *for the
  cycle-based path*; the plain step-index MSD retains the AMBER caveat.

### 2.2 CVE — covariance-based D & σ_loc estimator — 🔴 RED (textbook) → 🟡 mitigated
- **Implemented as:** the blur-aware form `σ² = R·⟨Δx²⟩ + (2R − 1)·cov₁` with R bounded from the cycle
  timing and the measured collection window (§I-B.1), over **consecutive integrated (single-cycle) step
  pairs** only; also returns the blur-independent `D_cve`. **The app-wide σ_loc comes from the
  photon-budget model** (§I-B.1); the CVE σ is a cross-check and D_cve is the D.
- **Source:** Vestergaard, Blainey & Flyvbjerg 2014.
- **MINFLUX caveat:** the CVE algebra is derived for a **constant** Δt. Applied naively across
  variable-Δt MINFLUX steps it biases D — hence the raw RED rating.
- **MINFLUXR's handling:** the covariance is restricted to **consecutive single-cycle pairs** that share
  a midpoint within one track (`n_cyc = round(Δt/t_cycle) == 1` on both steps), which pre-segments into
  constant-Δt
  sub-runs — precisely the mitigation the literature prescribes. ✅ Statistically valid on MINFLUX at
  the cost of discarding dark-gap-adjacent steps.
- **No photon weighting.** Photon counts differ between tracks for several reasons besides precision
  (mobility, bleaching, emitter brightness, background), so photon weighting would weight tracks by
  whatever drives those differences. Every step counts equally; `noise_frac` is reported.
- **The dominant caveat is the D floor.** The tracking loop's own wander adds an apparent D on an immobile
  target, and it grows steeply with the TCP diameter. The panel reports this floor beside D when it is
  entered; a published D needs the matching per-sequence measurement. This is why the in-app rating is 🟡.
- **Lateral vs axial σ:** σ is pooled over x, y and z, and in 3D MINFLUX the axial and lateral precision
  differ, so the purely lateral consumers (D-map, HMM) need the lateral value. `sigma_xy` / `sigma_z` /
  their ratio are reported separately (§I-B.1).

### 2.3 JDD — jump-distance / chi-mixture — 🟢 GREEN (with MINFLUXR's handling)
- **Implemented as:** jump-distance distribution with an N-component scaled-chi mixture (**Rayleigh in
  2D, Maxwell in 3D** — see I-B.4); σ_loc bias correction, √Δt-normalization (default on), and an
  **integrated-steps-only** mode.
- **Source:** Schütz, Schindler & Schmidt 1997; benchmarking in Weimann et al. 2013.
- **MINFLUX caveat:** small σ_loc reduces error-broadening of the Rayleigh peaks (an advantage), but
  jumps must be pooled at a **fixed real Δt**; mixing steps of different cycle counts pollutes the
  distribution.
- **MINFLUXR's handling:** two options — (1) **√Δt-normalise** each jump (default on) so mixed cycle
  counts pool correctly for free diffusion; or (2) **Integrated steps only** (`round(Δt/t_c) == 1`),
  which keeps a single true Δt so *no* pooling assumption is needed at all, valid regardless of the
  motion model. Both correct for localization error and motion blur. The "short track → few jumps"
  concern is weak for MINFLUX (100× more locs/track than camera SPT). Residual: the √Δt mode still
  assumes free diffusion; prefer integrated-only in doubt, and keep the component count low (each
  component still assumes free Brownian). The fitted law is the d-dimensional one, so 3D data is fitted
  with the Maxwell form.

### 2.4 HMM diffusive-state segmentation — 🟡 AMBER
- **Implemented as:** continuous-time HMM fitted by maximum likelihood across all tracks (§I-B.5):
  P(Δt) = expm(Q·Δt), emissions with blur and per-localization noise, rates as parameters, bootstrap CIs
  over tracks.
- **Source:** Das, Cairo & Coombs 2009 (HMM for diffusive states); Persson et al. 2013 (vbSPT, pooling
  short tracks); Jackson 2011 (msm: Markov jump processes observed at irregular times).
- **MINFLUX caveat:** a discrete per-step transition matrix assumes uniform Δt; MINFLUX Δt is irregular.
  The continuous-time generator removes that assumption exactly. What remains: one state per step, and
  shared-localization noise correlation.
- **MINFLUXR's handling:** rates are parameters of the likelihood rather than counts along the decoded
  path, which read dwell times long. Simulation sets the usable range (mean dwell ≥ 10 cycles), and
  states outside it are flagged in the app.
- **Why 🟡:** below that range the fit describes the step-size distribution rather than switching, which
  is common at MINFLUX cycle times; the flags say when that is the case.

### 2.5 Confinement analysis — 🟢 GREEN (with the N- and D-matched simulated null)
- **Implemented as:** the per-track **MSD confinement index** CI = (MSD_late/MSD_early)/(τ_late/τ_early)
  classified against a **simulated free-diffusion null** matched on localization count and, where noise
  matters, on per-track D (§I-B.6).
- **Source:** confinement-detection SPT lineage (Simson, Sheets & Jacobson 1995; Saxton; Kusumi 1993).
  The CI with a simulated, matched null is this tool's own construction.
- **MINFLUX fit:** CI inherits the cycle-index binning, so irregular Δt is handled; it needs no external D,
  so D heterogeneity does not inflate the confined fraction once the null is D-matched.
- **Residual caveats:** detection is limited by the free extent √(2dDT) relative to the domain, so
  **confined % is not comparable across datasets with different track lengths or mobilities**; CI is weak
  for directed motion. The null is blur-averaged and resamples photon σ; its false-positive rate on
  simulated free tracks is covered by the test suite.

### 2.6 Directionality / turn-angle autocorrelation — 🟢 GREEN
- **Implemented as:** turn-angle distribution + directional autocorrelation vs. the isotropic-Brownian null.
- **Source:** Burov et al. 2013 (distribution of directional change).
- **MINFLUX fit:** low σ_loc largely removes the classic **localization-noise anti-persistence
  artifact** that plagues turn-angle analysis of noisy camera steps. The Brownian null is
  Δt-independent, so irregular sampling is not a major confound. Still needs pooling across tracks for
  statistical power (short single tracks).

### 2.7 Anomalous-diffusion exponent α — 🔴 RED
- **Implemented as:** power-law fit to (time-averaged) MSD to extract α, per track and pooled.
- **Source:** Metzler et al. 2014 (models/pitfalls); Kepten et al. 2015 (fitting guidelines & bias maps).
- **MINFLUX caveat:** Kepten's own precision maps show α is **least reliable exactly in the short-track
  regime** MINFLUX produces, and irregular Δt cuts the number of valid lag bins for the power-law fit.
- **Recommendation:** report **pooled / ensemble** α with confidence intervals and short-track
  exclusion; treat **single-track α as qualitative**.

### 2.8 Change-point detection of motion states — 🟠 EXPLORATORY
- **Implemented as:** PELT segmentation of trajectories (via `changepoint`), min-segment / penalty controls.
- **Source:** Killick, Fearnhead & Eckley 2012 (PELT).
- **MINFLUX caveat:** PELT is sampling-agnostic *if* the per-segment cost is written in **true elapsed
  time**; off-the-shelf costs assume evenly spaced points. Short tracks limit resolvable segments.
- **MINFLUXR's handling:** the segmented signal is log(r²/Δt), a local log-D that is Δt-aware by
  construction, scaled to unit variance before PELT. The mixture's state count is confounded by segment
  length (§I-B.3), hence 🟠.

### Summary

| Method | Source | Textbook | With MINFLUXR's adaptation |
|---|---|---|---|
| MSD optimal-fit | Michalet 2010/2012; Vogler 2025 | 🟡 | 🟢 cycle-based path / 🟡 step-index |
| CVE (D) + photon-model σ_loc | Vestergaard 2014; Vogler 2025 | 🔴 | 🟡 (integrated pairs, blur-aware; wander floor on D) |
| JDD | Schütz 1997; Weimann 2013 | 🟡 | 🟢 (√Δt norm + integrated-only mode) |
| HMM (continuous-time MLE) | Das 2009; Persson 2013; Jackson 2011 | 🟡 | 🟡 rates by likelihood, validated range flagged (§2.4) |
| Confinement (MSD CI) | Simson 1995; tool-original null | 🟢 | 🟢 (N- and D-matched simulated null) |
| Turn-angle | Burov 2013 | 🟢 | 🟢 (elapsed-time autocorrelation) |
| Anomalous α | Metzler 2014; Kepten 2015 | 🔴 | 🔴 (pool + exclude short) |
| Change-point / PELT | Killick 2012 | 🟡 | 🟠 (variance-scaled cost; GMM state count confounded by segment length) |
| D-map (single-lag) | Vestergaard 2014 (short-lag limit) | 🟢 | 🟢 (mean estimator, noise and blur corrected) |

---

## 3. MINFLUX-specific parameter tools

- **D_max ceiling / trackable-D limit** — `D_max = (L²/16)·ψ / (N_PL + ψ·t_hw)`, cited in-code to
  **Vogler et al. 2025** (*Commun. Biol.*). Consistent with the independent SimuFLUX ceiling of
  ~3 µm²/s (Marin & Ries 2025). 🟢 Correctly attributed.
- **σ_loc noise floor / dark-cycle handling / continuity gate** — MINFLUXR-original engineering
  built on the data properties in §1.

---

## 3a. MSD on MINFLUX — computation, confinement & what is non-standard

The MSD tool follows the MINFLUX MSD workflow of **Vogler et al. 2025** (MinSpt), the only published MSD
treatment specific to MINFLUX. The research landscape here is sparse;
this section records what is *standard* (safe) versus *non-standard/unsolved* (deliberately not done).

### How MSD is (and should be) computed on MINFLUX — checklist
1. **Bin by integer cycle-count lag, not wall-clock time.** MINFLUX cycle-skipping makes Δt
   non-uniform; group displacement pairs by cycle-index η (Vogler 2025). The cycle-index binning does this.
2. **Segment traces at abnormal cycle gaps before MSD** (Vogler split when η > ⟨η⟩+SD). The gap-split
   + continuity gate is the equivalent; prefer integrated (single-cycle) steps with tolerance.
3. **Never fit a fixed/arbitrary lag range.** Michalet & Berglund 2012: a fixed choice ("first 4
   points", "10 %") has much larger uncertainty than the computed optimum p*. MINFLUXR uses OLSF + Michalet
   optimal-point selection for the **linear** model. *(NB: this optimal-p theory is derived for the
   linear/free-diffusion case ONLY — it is NOT valid for a plateauing confined curve; the
   confined/anomalous fits deliberately use the whole window instead.)*
4. **Both ends of the MSD are unreliable, for different reasons:** short lags — variance inflates when
   the reduced error x = σ²/(DΔt) > 1 (Michalet & Berglund 2012), plus probe/linker wobble; long lags —
   fewer, statistically *correlated* pairs (Qian, Sheetz & Elson 1991). Hence the fit window [lower,
   upper] cutoff. The default 5–50 cycles matches Vogler's `fit_limits=(5,50)`.
5. **The exact 5-cycle lower bound is probe- and dataset-specific.** Vogler's
   two-regime "kink" (5–9 lags) reflects their QD/PEG-linker probe's non-diffusive wobble; port the
   *diagnostic logic* (inspect the curve for where the low-lag regime ends) not the number.
6. **Model σ_loc and motion blur in the linear fit:** ρ(iΔt)=2d(σ²−2RDΔt)+2dD·iΔt. MINFLUX integrates
   photons over part of the step, so R ≈ 0.07–0.13 depending on the sequence (§I-B.1a; MinSpt's 1/6.2 overstates it); ignoring it under-reads σ from the intercept. Small σ_loc → small x → favours *fewer* MSD points, not more.
7. **Minimum track length is a JUDGEMENT CALL, not a literature constant.** There is **no
   MINFLUX-calibrated minimum-localizations rule**. The default of **50** follows MinSpt's practical
   choice (`filter_by_length(min_length=50)`). Kepten et al. 2015's
   L≤300 guidance is for camera SPT with far larger σ_loc and is not directly transferable.
8. **D_max: theoretical ≠ demonstrated.** Vogler's Eq. 3 gives D_max ≈ 10–11 µm²/s theoretically, but
   they only *demonstrated* tracking to D ≈ 2.5 µm²/s (lipids). The D-ceiling warning uses the
   theoretical form, so it may under-warn near the ~2.5–3 µm²/s practical ceiling — interpret the
   green/amber/red flag conservatively.

### Confined diffusion in closed cells
MinSpt (Vogler et al. 2025) covers free diffusion only and leaves confined and anomalous analysis to the
user, so confined MSD on MINFLUX is **non-standard**.
- **Standard confined-MSD model** (when free diffusion is ruled out): the plateauing corral form,
  MSD(t)=MSD∞·[1−A₁·exp(−4A₂D t/MSD∞)] introduced by **Kusumi, Sako & Yamamoto 1993**; the "generic"
  square-domain version MSD(t)=(L²/3)·[1−(8/π²)·exp(−t/τ)], τ=L²/(π²D). Companion corral/escape
  framework: **Saxton 1993/1994/1995**. Destainville & Salomé 2006 is a motion-blur correction to this
  form.
- **Detect confinement:** a resolvable MSD plateau vs. a straight line, or α<1 in MSD=4D·t^α (Saxton
  1994; Metzler & Klafter 2000; Metzler et al. 2014). α<1 has several possible causes besides
  confinement.
- **A rod has no closed-form MSD.** English et al. 2011 and Bakshi et al. 2011 handle it empirically
  (long- vs short-axis decomposition; Monte-Carlo of a spherocylinder). vbSPT (Persson et al. 2013)
  avoids geometry with an HMM over diffusive states, which answers a different question.

### Deliberately NOT implemented (non-standard / unsolved)
- Confined-MSD **combined with MINFLUX cycle-based irregular Δt**: not in the literature; Vogler et al.
  2025 apply the irregular-Δt machinery to free diffusion only.
- Michalet optimal-point selection on a **plateauing** curve (an unsupported extrapolation).
- Rod-cell **geometry correction** (no closed form).
- Non-Gaussian, viscoelastic and superstatistical cytoplasm models (e.g. Sadoon & Wang 2018): active
  research without consensus.

In the app these appear as the **Confined** / **Anomalous** fit models, labelled *beyond-standard for
MINFLUX* and exploratory, with a free-diffusion note on the Linear result.

---

## 6. SMLM (imaging) analysis tools — sources & appropriateness

The SMLM tab applies established single-molecule-localization-microscopy (SMLM) and spatial
point-pattern methods to MINFLUX localizations. MINFLUX-specific method development for *imaging*
statistics barely exists: **the estimators are standard; the MINFLUX
caveat is the labelling and blinking geometry**.

### 6.0 The pervasive caveat — labelling & blinking inflate every spatial statistic — 🟡 AMBER
A single molecule is reported as several localizations: fluorophores blink and re-appear, and a
multivalent antibody label scatters several dyes ~10–40 nm around one epitope. This **over-counts and
spatially spreads** each molecule, which inflates short-range nearest-neighbour peaks, Ripley/PCF
signal, and apparent cluster counts as a *label artefact*, not biology (Annibale et al. 2011;
Sengupta et al. 2011; Veatch et al. 2012). **MINFLUXR's handling:** the two-tier hierarchy
collapses localizations → trace median → **de-blinked molecule** → **median cluster** before spatial
statistics (de-blinking is the artifact correction; median clustering is the structural tier),
which removes most multi-blink over-counting; results on *raw* localizations should be read with this
bias in mind. It mitigates the bias but cannot recover the true molecule count below the label
footprint.

### 6.1 Nearest-neighbour and median clustering (DBSCAN/HDBSCAN/OPTICS) — 🟢 GREEN (method)
- **Source:** Clark & Evans 1954 (NN statistic); Ester, Kriegel, Sander & Xu 1996 (DBSCAN);
  Ankerst et al. 1999 (OPTICS); Campello, Moulavi & Sander 2013 (HDBSCAN).
- **Applies to the *structural* tier only** (`cluster2`). De-blinking is a different method with a
  different provenance — see §6.1b.
- **Fit:** density-based clustering is well-adapted to SMLM and widely used; DBSCAN in particular is a
  community standard. MINFLUX's fine precision *helps*. The main caveat is §6.0 (feed de-blinked
  medians, tune `eps2` to the structure scale, not the label footprint).
- **Membership refinements** (core-only / prune / min members, §I-A.2b) — the core-only option is
  literature-backed (Campello et al. 2013 define clusters on core points precisely to avoid DBSCAN's
  border-point assignment, which is also order-dependent when two clusters reach the same border
  point). The **core-distance prune with a local baseline is tool-original** 🟡 and checked against
  simulated structures of six shapes. The min-members filter is a plain size threshold.
- ⚠️ **Not resolvable by any of these:** separating a genuinely diffuse periphery of a structure from
  background centroids. Density-based clustering has no information distinguishing them. Refinement
  settings are declared choices; state them alongside any cluster count or structure size.

### 6.1b De-blinking (merging re-blinks) — 🟡 AMBER (tool-original)
- **Source:** none directly. The *problem* is standard (Annibale et al. 2011 established blink-merging
  and the dark-time calibration for camera SMLM), but Annibale's temporal grouping does not transfer:
  a MINFLUX trace ID *is* already one continuous emission event, so the instrument has done the
  temporal part and what remains is merging **separate traces** that are re-blinks of one fluorophore.
  The centroid and majority-vote criteria used here (§I-A.2a) are tool-original.
- **Validation:** against simulated ground truth, where centroid
  linkage makes fewer split and merge errors than linking raw localizations. The SMLM simulator places a
  fluorophore's re-blinks at the **same true position** (docs/simulator.md §7), so that test covers the
  zero-motion case only; movement of the label between blinks is not simulated. Threshold *k* can be
  measured per dataset from residual re-blink excess (§I-A.2a).
- **Caveat (why AMBER):** the threshold depends on how far a label moves *between* blinks,
  and that quantity cannot be measured without already knowing which traces are the same molecule —
  the question being asked. The residual-excess test bounds it from below only. Report the threshold
  and the de-blinked count together; a molecule count from this tier still depends on photophysics.
- **Not applicable:** Nyquist/Shannon sampling arguments (Shroff et al. 2008) constrain *labelling
  density* and do not set this threshold.

### 6.1a Traces per cluster — 🟡 AMBER
- **Source:** cluster-occupancy counting in the SMLM molecular-counting lineage — Annibale et al. 2011
  (blinking artefacts in PALM); Sengupta et al. 2011 (PC-PALM).
- **Fit:** counting *traces* is the right unit for MINFLUX and removes the
  within-trace over-count for free.
- **Caveat:** it remains a **blinking-event count, not stoichiometry** — see §I-A.10. The Annibale
  correction is not applicable here because no dark-time distribution is calibrated. Relative comparisons
  between identically-acquired conditions only.

### 6.2 Ripley's L(r), bivariate L₁₂(r), pair cross-correlation g₁₂(r) — 🟢 GREEN (method) / 🟡 (labelling)
- **Source:** Ripley 1977; Besag 1977 (L-transform); Lotwick & Silverman 1982 and Dixon 2002
  (bivariate/cross statistics); Sengupta et al. 2011 (PC-PALM) and Veatch et al. 2012 (correlation
  functions, over-counting) for the SMLM pair-correlation lineage.
- **Fit:** correct point-pattern estimators. The §6.0 labelling bias is the dominant confound — PC-PALM
  was *designed* to separate multiple appearances of one protein from real clusters, so citing it also
  signals the correct mitigation.
- **Edge bias.** A translation edge correction (on by default) and an opt-in sample-footprint window
  handle it (§I-A.3), and the Monte-Carlo envelope simulates in the identical window, so any remaining
  edge bias is present in the band too. The envelope does **not** cancel a window that is mostly empty
  space: for a cell that does not fill its frame, the footprint window is what removes spurious
  "clustering at every scale".
- **✅ No inter-channel registration error — a genuine MINFLUX advantage.** Unlike camera-based two-colour
  SMLM, the DCR (ratiometric) two-colour mode excites both fluorophores with the **same laser** and
  localizes them through the **same beam path in the same iterative sequence**, so both channels are
  reported in one common coordinate frame. There is no chromatic offset and no affine channel
  registration to calibrate — the usual dominant error term of two-colour co-localization simply does not
  exist here, and L₁₂/g₁₂ are trustworthy right down to the localization precision. *(Exception: the
  multiplex-by-files and multiplex-by-time-window modes overlay **separate acquisitions or separate time
  windows**, so sample or stage drift between rounds can offset them. This applies only to those modes.)*
- **⚠ The real two-colour error is DCR mis-assignment (crosstalk), not offset.** The channel split is a
  threshold on a **continuous DCR distribution whose two populations overlap**, so a fraction of
  molecules is assigned to the wrong channel. This does not shift distances — it moves *points between
  the two sets*, which **manufactures short-range cross-channel pairs**: a molecule wrongly labelled
  channel B sits, by construction, right next to its true channel-A neighbours. Crosstalk therefore
  biases L₁₂ and g₁₂ **upward at short r — toward apparent co-localization**, exactly where a real
  interaction would show. Asymmetric labelling density between the two channels also skews the
  estimators. Finally, co-localization at nm scale is **proximity, not interaction**.
- **Quantified (`auto_detect_dcr_split`).** A two-component Gaussian mixture is fitted to the DCR values;
  component j has mean μ_j, SD σ_j and weight π_j, and X_j is a DCR value drawn from it. Alongside the
  threshold s it returns:
  - **d′ = |μ₂−μ₁| / √((σ₁²+σ₂²)/2)** — standardised separation of the two modes. Verdicts used in-UI:
    **d′ ≥ 3 clean, ≥ 2 workable, < 2 poor** (channels overlap heavily).
  - **Expected mis-assignment rate** `π₁·P(X₁>s) + π₂·P(X₂≤s)`, and the more directly usable
    **per-group contamination** — what fraction of everything *called* group k actually came from the
    other dye: `contam₁ = π₂P(X₂≤s) / [π₁P(X₁≤s) + π₂P(X₂≤s)]`, symmetrically for group 2.
  - **Ambiguous fraction** — traces whose best GMM posterior is < 0.95, i.e. the ones that would switch
    channel under a small threshold change. These are the natural candidates for a strict gate.

  **Read any cross-channel result against these numbers:** a co-localized fraction of the same order as
  the contamination is not evidence of interaction.

  On synthetic data the predicted mis-assignment rate matches the true one for a good split and drifts
  for a poor one, so **"poor" means "do not use"**.
  Contamination is correctly asymmetric under unequal populations.
- **Scope of d′: it grades a split, it does not detect one.** The GMM is a sound estimator of the split
  and is the right tool here — the point is only that requesting G = 2 always yields two components,
  because a single-dye DCR distribution is **skewed** and a two-Gaussian mixture fits skewed unimodal
  data better than one Gaussian, so a single-colour file can return a "workable" split. A BIC G = 1 vs
  G = 2 comparison is not used as a gate, because it detects skewness rather than bimodality. Use d′ and
  the crosstalk numbers to judge **how good** a two-colour split is, once the experiment establishes
  that there are two dyes.
- **✅ Per-molecule split confidence (`dcr_split_confidence`).** The aggregate crosstalk numbers cannot
  say *which* molecules are uncertain, and those are exactly the ones that manufacture false short-range
  cross-channel pairs. For each trace the fitted mixture gives the posterior of its assigned component,
  `post_k(v) = π_k·N(v|μ_k,σ_k) / Σ_j π_j·N(v|μ_j,σ_j)` — 1.0 unambiguous, 0.5 sitting exactly on the
  split — plus a log-odds form and an `on_edge` flag (posterior < 0.95). Shown in the 2-colour panel as
  posterior-vs-DCR with the edge band highlighted, and exportable per trace. Edge traces are flagged but
  not excluded from the cross-channel tools.
  On synthetic data the posterior is close to calibrated, and dropping the edge traces removes most of
  the crosstalk for a small loss of data, which is what makes a strict gate worth offering.
- **✅ Pool DCR at upload — the most effective mitigation, where it is available.** `dcr_pooled`
  (`helpers_shared.R`) is the photon-weighted DCR across the whole iteration sequence, `Σ(dcr·eco)/Σeco`.
  Averaging over the sequence with photon weighting **sharpens the
  two modes substantially**, so the split quality — and therefore the crosstalk floor on every
  cross-channel statistic — improves at source. This beats any downstream statistical correction: the
  split is the error, so fix the split.
  ⚠ **Only .npy and .json uploads can be pooled** — pooling happens in `read_minflux_npy_as_df()` /
  `read_minflux_json_as_df()` via their `pool_dcr` argument. A **.csv is a pyMINFLUX export of the final
  iteration**, already through that pipeline, so there is nothing left to pool and the UI must not nag
  about it. The in-app hint is therefore conditional: it points at the `dcr_pooled` column when one
  exists and is unused, and otherwise explains the format restriction.
  The split-quality numbers are only comparable between runs that used the **same** DCR column.
- **g₁₂ vs L₁₂ on MINFLUX data.** g₁₂ is the more informative of the two here because it is scale-local:
  a short-r blinking/labelling artefact stays confined to the first few bins, whereas in the cumulative
  L₁₂ it contaminates every larger r. The price is noise — each annulus holds few pairs, so short-r bins
  are unstable and an apparent first-bin peak is often binning.

### 6.2b Labelling-density thinning (opt-in) — 🟡 AMBER
- **Source:** independent p-thinning leaves g(r) and K(r) invariant (Illian et al. 2008; Chiu et al.
  2013); support from a smoothed-density region (Hyndman 1996) or morphological occupancy (Serra 1982);
  area-free NN density cross-check (Thompson 1956).
- **Fit:** the invariance is exact, and whole de-blinked units are dropped, so traces and per-trace
  spread are untouched. Transient transfection makes this correction relevant for most MINFLUX cell data.
- **Caveat (why AMBER):** the thinning target is n/support, so the **support estimate, above all its
  gap-fill radius, is the dominant error**. Thinning is one-way, can push small clusters below minPts,
  and removes the effect being measured when density is itself the biology.

### 6.2c Local density & hotspots — 🟡 AMBER (tool-original calibration)
- **Source:** k-NN density (Loftsgaarden & Quesenberry 1965); index of dispersion (Fisher, Thornton &
  Mackenzie 1922); DBSCAN (Ester et al. 1996). The per-cell eps calibration against CSR inside the cell's
  own footprint, and the choice of the largest focus as the headline statistic, are this tool's.
- **Fit:** each statistic is dimensionless relative to the cell's own mean density, which answers the
  failure of central NN statistics on transfected cells (§I-A.3b).
- **Caveat (why AMBER):** the calibration is tool-original. On random patterns cv_ratio, VMR,
  frac_above and max_focus_ratio hold their chance baseline across density and footprint shape, while raw
  focus count and max_focus_n move with density; compare cells on the former.

### 6.2a Cross-channel NN & object-based co-localization curve — 🟢 GREEN (method)
- **Source:** bivariate nearest-neighbour distances (Clark & Evans 1954; Dixon 2002) with an
  **object-based** co-localization read-out (Lachmanovich et al. 2003; Dunn et al. 2011).
- **Fit:** object-based co-localization is the right family for MINFLUX — there are no pixels and no
  intensity image, only coordinates, and nm precision makes the distance axis physically meaningful down
  to the label size. This is strictly better founded than an intensity-correlation coefficient here.
- **Caveat:** no null model (§I-A.1) — every A point has a nearest B point, so judge against L₁₂/g₁₂
  envelopes. DCR crosstalk applies exactly as in §6.2, and hits cross-NN hardest of all the two-colour
  tools: a mis-assigned molecule contributes a near-zero nearest-neighbour distance, which is precisely
  the first bin of the histogram and the steep part of the co-localization curve.

### 6.3 Radius of gyration / gyration-tensor anisotropy — 🟢 GREEN
- **Source:** standard second-moment cluster size; anisotropy eigenvalue coefficients Rudnick & Gaspari
  1986; linear/planar/spherical coefficients Westin et al. 1997/2002.
- **MINFLUXR's handling:** computed on **per-TID median positions**, one point per trace. This matters — on raw
  localizations each molecule is weighted by how long it was tracked, so a single bright fluorophore
  contributing thousands of localizations pulls Rg toward itself and shrinks it. De-blinking first makes
  Rg a property of the *structure*.
- **Caveat:** Rg is inflated by the **label footprint** (a multivalent antibody scatters dyes ~10–40 nm
  from the epitope), so structures near that scale are systematically over-sized; **no footprint
  deconvolution is applied**. Rg is also biased low at small member counts (the centroid is estimated
  from the same points), and it is a single number — a filament and a disc of equal extent share an Rg,
  which is what the anisotropy coefficients are for.

### 6.4 KDE density map + Pearson co-localization — 🟡 AMBER
- **Source:** kernel density estimation (Silverman 1986); Pearson correlation of two density surfaces as
  a co-localization index (cf. Dunn et al. 2011).
- **Caveat:** the Pearson r is strongly bandwidth-dependent (small bandwidth → r≈0 even for true
  co-localization; large → inflated). Now made explicit by the optional **r(bandwidth) sweep** plot, and
  the single-bandwidth r is accompanied by a **random-shift spatial-randomization p-value** (§I-A.5) that
  tests it against independent placement. For a bandwidth-free, scale-resolved alternative the UI points
  to the pair cross-correlation g₁₂(r) (§6.2). Still treat the raw r as a *relative* measure.
- **Deeper caveat — the "image" is manufactured.** MINFLUX has no pixels, so the density surface the
  Pearson r is computed on does not exist until the bandwidth creates it: the bandwidth is a free
  parameter. Only compare r values computed at the **same bandwidth on
  similarly dense data**; prefer the sweep to any single number.
- **The p-value depends on the window** (§I-A.5). Both nulls randomise over the analysis window. With
  the bounding box (the default) they are anti-conservative for channels that share one cell inside a
  larger frame; with the footprint switch on the test is conservative. Either way, take a positive claim
  from g₁₂(r) with its envelope.

### 6.5 Cluster shape fitting (viewer) — 🟠 EXPLORATORY
- **Method (standard):** orthogonal-distance MLE + AICc/BIC model selection (Akaike 1974; Schwarz 1978;
  Hurvich & Tsai 1989).
- **Why exploratory:** the shape dictionary (sphere/ellipsoid/helix/ring/arc) is imposed by this tool and
  the "best" shape is a model-selection outcome; with few,
  blinking-inflated points several shapes are statistically indistinguishable. **The true shape may not
  be in the dictionary at all — the ranking will still name a winner.**
- **Caveat (residual model):** the information criteria assume iid isotropic-Gaussian residuals, but
  MINFLUX precision varies with photon count, so points are not equally precise (§I-A.7). The likelihood
  is approximate; AICc differences of a few units are not meaningful. RMSD always favours the most
  flexible model — prefer AICc.

### 6.6 Ground-truth position/count estimate & structure size — 🟠 EXPLORATORY
- Labelling-model inverse developed in this tool, checked against the simulator. Contour via
  minimum-spanning-tree (Prim 1957); anisotropy coefficients as §6.3. Reported as a range below the labelling resolution.

### 6.7 Detection-over-time / crowding-vs-timing / NN-vs-size — 🟠 EXPLORATORY
- Bespoke descriptive diagnostics of this tool with **no published MINFLUX equivalent**. They are
  correlational: detection timing and crowding co-vary with photobleaching, blinking and labelling, so a
  trend is a hypothesis to test with controls. The detection-rate decay used by the simulator is
  empirical.
- **Three specific pitfalls beyond "correlational"** (detail and rationale in §I-A.9):
  1. **The KS gap is an effect size, not a test.** The clustered and background samples are not
     independent (same acquisition, same detection budget), so the KS null does not hold. The code
     correctly reports only the gap — do not convert it to a p-value.
  2. **Crowding-vs-timing has a positive trend built in, corrected by a permutation null** (§I-A.9). A within-file
     label-permutation null gives each structure its own 95% interval and a Spearman trend test, so the
     slope is read against the no-crowding expectation rather than against zero. The band alone still
     reflects offset + trend; the ρ test isolates the trend.
  3. **NN-vs-size is right-censored under DBSCAN** at the link radius `eps2` (default 200 nm), so a
     flattening curve may be that ceiling. Re-run under HDBSCAN/OPTICS, which
     impose no such bound, before interpreting a slope.

### 6.8 QC tools: sample drift & fiducial scan — 🟡 AMBER
- **Sample drift & jitter (§I-A.2e).** Regression of signed re-blink displacement on signed time gap
  through the origin, i.e. a standard through-origin least-squares velocity with the whole acquisition as
  lever arm. Sound in construction: pairs of different molecules widen the interval without biasing the
  slope. Caveat: selecting pairs by separation truncates large displacements, so fast drift is
  **understated** (a lower bound).
- **Fiducial / contaminant scan (§I-A.2f).** Per-file lifetime tail (duration ≥ 5× the file's own q99)
  with a brightness rule that requires persistence. The physical basis is that a fiducial does not
  bleach, so the tracker stays locked on it, and the threshold is relative to the file, so it adapts to
  the acquisition. Report-only by default. A fiducial that is neither long-lived nor bright, for example
  one the sequence never locks onto, would be missed.

### SMLM summary

| Tool | Source | Rating |
|---|---|---|
| NN / median clustering (DBSCAN / HDBSCAN / OPTICS) | Clark & Evans 1954; Ester 1996; Ankerst 1999; Campello 2013 (also the core-only option) | 🟢 (feed de-blinked medians) |
| Median-cluster membership prune (core distance, local baseline) | this tool | 🟡 tool-original |
| De-blinking (centroid / majority vote) | tool-original; problem framing after Annibale 2011 | 🟡 threshold assumption (§6.1b) |
| Traces per cluster | Annibale 2011; Sengupta 2011 | 🟡 blinking (not stoichiometry) |
| Ripley L / L₁₂ / PCF g₁₂ | Ripley 1977; Besag 1977; Lotwick & Silverman 1982; Dixon 2002; Sengupta 2011; Veatch 2012; Illian 2008 (edge correction) | 🟢 / 🟡 labelling, window choice |
| Cross-channel NN & co-loc curve | Clark & Evans 1954; Dixon 2002; Lachmanovich 2003; Dunn 2011 | 🟢 (no null of its own) |
| Radius of gyration / anisotropy | Rudnick & Gaspari 1986; Westin 1997/2002 | 🟢 |
| Labelling-density thinning | Illian 2008; Chiu 2013; Hyndman 1996; Serra 1982 | 🟡 support estimate (§6.2b) |
| Local density & hotspots | Loftsgaarden & Quesenberry 1965; Fisher 1922; Ester 1996; tool-original calibration | 🟡 tool-original calibration (§6.2c) |
| Sample drift & jitter | through-origin regression on re-blink pairs | 🟡 lower bound (§6.8) |
| Fiducial / contaminant scan | per-file lifetime tail + brightness; tool-original | 🟡 not validated by other groups; report-only by default (§6.8) |
| KDE + Pearson co-loc | Silverman 1986; Dunn 2011 (*not* the Costes 2004 algorithm — see §6.4) | 🟡 bandwidth; p-value needs footprint on |
| Shape fitting | Akaike 1974; Schwarz 1978; Hurvich & Tsai 1989; Pearson 1901; Weiszfeld 1937 (centre) | 🟠 exploratory |
| GT estimate & structure size | this tool; Prim 1957 | 🟠 exploratory |
| Detection-time / crowding / NN-vs-size | this tool; Cleveland 1979; Barber 1996 | 🟠 exploratory |

---

## 4. Bibliography

**MINFLUX microscopy & tracking**
1. Balzarotti F, Eilers Y, Gwosch KC, Gynnå AH, Westphal V, Stefani FD, Elf J, Hell SW (2017). *Nanometer resolution imaging and tracking of fluorescent molecules with minimal photon fluxes.* Science 355(6325):606–612. doi:10.1126/science.aak9913
2. Gwosch KC, Pape JK, Balzarotti F, Hoess P, Ellenberg J, Ries J, Hell SW (2020). *MINFLUX nanoscopy delivers 3D multicolor nanometer resolution in cells.* Nature Methods 17:217–224. doi:10.1038/s41592-019-0688-0
3. Schmidt R, Weihs T, Wurm CA, Jansen I, Rehman J, Sahl SJ, Hell SW (2021). *MINFLUX nanometer-scale 3D imaging and microsecond-range tracking on a common fluorescence microscope.* Nature Communications 12:1478. doi:10.1038/s41467-021-21652-z
4. Wolff JO, Scheiderer L, Engelhardt T, Engelhardt J, Matthias J, Hell SW (2023). *MINFLUX dissects the unimpeded walking of kinesin-1.* Science 379(6636):1004–1010. doi:10.1126/science.ade2650
5. Deguchi T, Iwanski MK, Schentarra EM, et al. (2023). *Direct observation of motor protein stepping in living cells using MINFLUX.* Science 379(6636):1010–1015. doi:10.1126/science.ade2676
6. Ostersehlt LM, Jans DC, Wittek A, et al. (2022). *DNA-PAINT MINFLUX nanoscopy.* Nature Methods 19:1072–1075. doi:10.1038/s41592-022-01577-1
7. Vogler BTL, De Angelis G, Zhao Z, Eggeling C, Reina F (2025). *Parameter optimization for MINFLUX microscopy enabled single-particle tracking.* Communications Biology 8:1573. doi:10.1038/s42003-025-09060-1. Reference implementation: *Toolbox_for_MINFLUX_SPT* (MinSpt / `the_road_to_minflux_msd`) — source of the MINFLUX Time-Signature Analysis (MFX-TSA) and the cycle-based MSD/OLSF workflow ported here.
8. Scheiderer L, Marin Z, Ries J (2025). *MINFLUX achieves molecular resolution with minimal photons.* Nature Photonics 19(3):238–247. doi:10.1038/s41566-025-01625-0
9. Marin Z, Ries J (2025). *Evaluating MINFLUX experimental performance in silico* (SimuFLUX). Nature Communications 17:246. doi:10.1038/s41467-025-66952-w

**SPT analysis methods**
10. Michalet X (2010). *Mean square displacement analysis of single-particle trajectories with localization error.* Phys. Rev. E 82:041914. doi:10.1103/PhysRevE.82.041914
11. Michalet X, Berglund AJ (2012). *Optimal diffusion coefficient estimation in single-particle tracking.* Phys. Rev. E 85:061916. doi:10.1103/PhysRevE.85.061916
12. Vestergaard CL, Blainey PC, Flyvbjerg H (2014). *Optimal estimation of diffusion coefficients from single-particle trajectories.* Phys. Rev. E 89:022726. doi:10.1103/PhysRevE.89.022726
13. Schütz GJ, Schindler H, Schmidt T (1997). *Single-molecule microscopy on model membranes reveals anomalous diffusion.* Biophys. J. 73:1073–1080.
14. Weimann L, Ganzinger KA, McColl J, et al. (2013). *A Quantitative Comparison of Single-Dye Tracking Analysis Tools Using Monte Carlo Simulations.* PLOS ONE 8:e64287. doi:10.1371/journal.pone.0064287
15. Das R, Cairo CW, Coombs D (2009). *A Hidden Markov Model for Single Particle Tracks…* PLoS Comput. Biol. 5:e1000556.
16. Persson F, Lindén M, Unoson C, Elf J (2013). *Extracting intracellular diffusive states and transition rates from single-molecule tracking data* (vbSPT). Nature Methods 10:265–269.
17. Simson R, Sheets ED, Jacobson K (1995). *Detection of temporary lateral confinement of membrane proteins using single-particle tracking analysis.* Biophys. J. 69:989–993.
18. Burov S, Tabei SMA, Huynh T, et al. (2013). *Distribution of directional change as a signature of complex dynamics.* PNAS 110(49):19689–19694. doi:10.1073/pnas.1319473110
19. Metzler R, Jeon JH, Cherstvy AG, Barkai E (2014). *Anomalous diffusion models and their properties…* Phys. Chem. Chem. Phys. 16:24128–24164. doi:10.1039/c4cp03465a
20. Kepten E, Weron A, Sikora G, Burnecki K, Garini Y (2015). *Guidelines for the Fitting of Anomalous Diffusion Mean Square Displacement Graphs…* PLOS ONE 10:e0117722. doi:10.1371/journal.pone.0117722
21. Killick R, Fearnhead P, Eckley IA (2012). *Optimal Detection of Changepoints With a Linear Computational Cost* (PELT). J. Am. Stat. Assoc. 107(500):1590–1598. doi:10.1080/01621459.2012.737745

**MSD on MINFLUX & confined diffusion (§3a)**
25. Qian H, Sheetz MP, Elson EL (1991). *Single particle tracking. Analysis of diffusion and flow in two-dimensional systems.* Biophys. J. 60:910–921. doi:10.1016/S0006-3495(91)82125-7
26. Kusumi A, Sako Y, Yamamoto M (1993). *Confined lateral diffusion of membrane receptors as studied by single particle tracking (nanovid microscopy)…* Biophys. J. 65:2021–2040. doi:10.1016/S0006-3495(93)81253-0
27. Saxton MJ (1993). *Lateral diffusion in an archipelago. Single-particle diffusion.* Biophys. J. 64:1766–1780.
28. Saxton MJ (1994). *Anomalous diffusion due to obstacles: a Monte Carlo study.* Biophys. J. 66:394–401. doi:10.1016/S0006-3495(94)80789-1
29. Saxton MJ (1995). *Single-particle tracking: effects of corrals.* Biophys. J. 69:389–398. doi:10.1016/S0006-3495(95)79911-8
30. Destainville N, Salomé L (2006). *Quantification and correction of systematic errors due to detector time-averaging in single-molecule tracking experiments.* Biophys. J. 90:L17–L19. doi:10.1529/biophysj.105.075176 *(motion-blur correction to the confined-MSD form, not the confinement derivation)*
31. Metzler R, Klafter J (2000). *The random walk's guide to anomalous diffusion: a fractional dynamics approach.* Phys. Rep. 339:1–77. doi:10.1016/S0370-1573(00)00070-3
32. English BP, Hauryliuk V, Sanamrad A, et al. (2011). *Single-molecule investigations of the stringent response machinery in living bacterial cells.* PNAS 108:E365–E373. doi:10.1073/pnas.1102255108
33. Bakshi S, Bratton BP, Weisshaar JC (2011). *Subdiffraction-limit study of Kaede diffusion and spatial distribution in live Escherichia coli.* Biophys. J. 101:2535–2544. doi:10.1016/j.bpj.2011.10.013
34. Sadoon AA, Wang Y (2018). *Anomalous, non-Gaussian, viscoelastic, and age-dependent dynamics of histone-like nucleoid-structuring proteins in live E. coli.* Phys. Rev. E 98:042411. doi:10.1103/PhysRevE.98.042411 *(experimental/non-standard, noted for contrast)*

**SMLM / imaging spatial-analysis methods (§6)**
35. Clark PJ, Evans FC (1954). *Distance to nearest neighbor as a measure of spatial relationships in populations.* Ecology 35:445–453. doi:10.2307/1931034
36. Ester M, Kriegel HP, Sander J, Xu X (1996). *A density-based algorithm for discovering clusters in large spatial databases with noise* (DBSCAN). KDD-96:226–231.
37. Ankerst M, Breunig MM, Kriegel HP, Sander J (1999). *OPTICS: Ordering points to identify the clustering structure.* ACM SIGMOD Record 28(2):49–60. doi:10.1145/304181.304187
38. Campello RJGB, Moulavi D, Sander J (2013). *Density-based clustering based on hierarchical density estimates* (HDBSCAN). PAKDD 2013, LNCS 7819:160–172. doi:10.1007/978-3-642-37456-2_14
39. Ripley BD (1977). *Modelling spatial patterns.* J. R. Stat. Soc. B 39(2):172–212.
40. Besag J (1977). *Contribution to the discussion of Dr Ripley's paper* (the L-transform of K). J. R. Stat. Soc. B 39:193–195.
41. Lotwick HW, Silverman BW (1982). *Methods for analysing spatial processes of several types of points.* J. R. Stat. Soc. B 44(3):406–413.
42. Dixon PM (2002). *Ripley's K function.* In: Encyclopedia of Environmetrics, Wiley. doi:10.1002/9780470057339.var046
43. Sengupta P, Jovanovic-Talisman T, Skoko D, et al. (2011). *Probing protein heterogeneity in the plasma membrane using PALM and pair correlation analysis* (PC-PALM). Nature Methods 8:969–975. doi:10.1038/nmeth.1704
44. Veatch SL, Machta BB, Shelby SA, et al. (2012). *Correlation functions quantify super-resolution images and estimate apparent clustering due to over-counting.* PLoS ONE 7:e31457. doi:10.1371/journal.pone.0031457
45. Annibale P, Vanni S, Scarselli M, Rothlisberger U, Radenovic A (2011). *Identification of clustering artifacts in photoactivated localization microscopy.* Nature Methods 8:527–528. doi:10.1038/nmeth.1627
46. Silverman BW (1986). *Density Estimation for Statistics and Data Analysis.* Chapman & Hall.
47. Dunn KW, Kamocka MM, McDonald JH (2011). *A practical guide to evaluating colocalization in biological microscopy.* Am. J. Physiol. Cell Physiol. 300:C723–C742. doi:10.1152/ajpcell.00462.2010
48. Akaike H (1974). *A new look at the statistical model identification.* IEEE Trans. Autom. Control 19(6):716–723. doi:10.1109/TAC.1974.1100705
49. Schwarz G (1978). *Estimating the dimension of a model* (BIC). Ann. Stat. 6(2):461–464. doi:10.1214/aos/1176344136
50. Hurvich CM, Tsai CL (1989). *Regression and time series model selection in small samples* (AICc). Biometrika 76(2):297–307. doi:10.1093/biomet/76.2.297
51. Rudnick J, Gaspari G (1986). *The asphericity of random walks.* J. Phys. A: Math. Gen. 19:L191–L193. doi:10.1088/0305-4470/19/4/004
52. Westin CF, Peled S, Gudbjartsson H, Kikinis R, Jolesz FA (1997). *Geometrical diffusion measures for MRI from tensor basis analysis.* ISMRM 1997. (planar/linear/spherical anisotropy coefficients; see also Westin et al. 2002, Med. Image Anal. 6:93–108.)
53. Prim RC (1957). *Shortest connection networks and some generalizations* (minimum spanning tree). Bell Syst. Tech. J. 36(6):1389–1401. doi:10.1002/j.1538-7305.1957.tb01515.x
54. Costes SV, Daelemans D, Cho EH, Dobbin Z, Pavlakis G, Lockett S (2004). *Automatic and quantitative measurement of protein-protein colocalization in live cells.* Biophys. J. 86(6):3993–4003. doi:10.1529/biophysj.103.038422 *(the randomization-based significance test behind the KDE panel's p-value)*
55. Lachmanovich E, Shvartsman DE, Malka Y, Botvin C, Henis YI, Weiss AM (2003). *Co-localization analysis of complex formation among membrane proteins by computerized fluorescence microscopy: application to immunofluorescence co-patching studies.* J. Microsc. 212(Pt 2):122–131. doi:10.1046/j.1365-2818.2003.01239.x *(object-based, distance-threshold co-localization — the family the cross-NN ECDF curve belongs to)*
56. Cleveland WS (1979). *Robust locally weighted regression and smoothing scatterplots* (loess). J. Am. Stat. Assoc. 74(368):829–836. doi:10.1080/01621459.1979.10481038
57. Barber CB, Dobkin DP, Huhdanpaa H (1996). *The quickhull algorithm for convex hulls.* ACM Trans. Math. Softw. 22(4):469–483. doi:10.1145/235815.235821 *(Qhull, via `geometry::convhulln` — the convex-hull cluster-size option)*
58. Pearson K (1901). *On lines and planes of closest fit to systems of points in space* (PCA). Philosophical Magazine 2(11):559–572. doi:10.1080/14786440109462720
59. Weiszfeld E (1937). *Sur le point pour lequel la somme des distances de n points donnés est minimum.* Tôhoku Math. J. 43:355–386. *(the geometric-median iteration used by `geometric_median()`; see Vardi & Zhang 2000, PNAS 97(4):1423–1426, for the modern convergence treatment)*
60. Illian J, Penttinen A, Stoyan H, Stoyan D (2008). *Statistical Analysis and Modelling of Spatial Point Patterns.* Wiley. *(independent thinning invariance of g and K; translation edge correction)*
61. Chiu SN, Stoyan D, Kendall WS, Mecke J (2013). *Stochastic Geometry and Its Applications*, 3rd ed. Wiley.
62. Baddeley A, Rubak E, Turner R (2015). *Spatial Point Patterns: Methodology and Applications with R.* CRC Press. *(translation edge correction, set covariance)*
63. Hyndman RJ (1996). *Computing and graphing highest density regions.* Am. Stat. 50(2):120–126. doi:10.1080/00031305.1996.10474359 *(KDE footprint)*
64. Serra J (1982). *Image Analysis and Mathematical Morphology.* Academic Press. *(morphological-closing footprint)*
65. Thompson HR (1956). *Distribution of distance to nth neighbour in a population of randomly distributed individuals.* Ecology 37(2):391–394. doi:10.2307/1933159
66. Loftsgaarden DO, Quesenberry CP (1965). *A nonparametric estimate of a multivariate density function.* Ann. Math. Stat. 36(3):1049–1051. doi:10.1214/aoms/1177700079 *(k-NN density, hotspot enhancement)*
67. Fisher RA, Thornton HG, Mackenzie WA (1922). *The accuracy of the plating method of estimating the density of bacterial populations.* Ann. Appl. Biol. 9(3–4):325–359. doi:10.1111/j.1744-7348.1922.tb05962.x *(index of dispersion)*

**Reviews**
22. Saxton MJ, Jacobson K (1997). *Single-Particle Tracking: Applications to Membrane Dynamics.* Annu. Rev. Biophys. Biomol. Struct. 26:373–399. doi:10.1146/annurev.biophys.26.1.373
23. Manzo C, Garcia-Parajo MF (2015). *A review of progress in single particle tracking…* Rep. Prog. Phys. 78:124601. doi:10.1088/0034-4885/78/12/124601
24. Shen H, Tauzin LJ, Baiyasi R, et al. (2017). *Single Particle Tracking: From Theory to Biophysical Applications.* Chem. Rev. 117(11):7331–7376. doi:10.1021/acs.chemrev.6b00815

---

## 4a. A note on the citations

Citations were checked against the primary literature; journal references with a DOI resolve on Crossref
with authors, title, journal, volume, pages and year as cited. Weiszfeld 1937 has no DOI and is cited in
its standard form.
