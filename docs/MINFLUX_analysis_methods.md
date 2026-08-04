# MINFLUX Analysis — Methods, Mathematics & Appropriateness

*Reference document for the analysis and simulation modules of MINFLUXR.*
Last reviewed: 2026-07-24 · Compiled from the source code and from verified primary literature (see [Bibliography](#bibliography)).

> **In-app counterpart.** Every analysis panel in the Tracking tab and the SMLM **Analysis** tab carries a
> collapsible *"ℹ Method & MINFLUX adaptation"* box (`.method_box()` in `helpers_shared.R`) with the
> Origin / MINFLUX adaptation / Caveat triple and a 🟢/🟡/🟠 chip. Those boxes are the user-facing summary
> of §I-A and §6 here — **keep the two in sync when either changes.**

This document has **two layers**, both organised around the app's three method families —
**SMLM (imaging)**, **Tracking (SPT)**, and **Simulation**:

- **Part I — Implementation & Mathematics ([§I](#part-i--implementation--mathematics)).** What each tool
  actually computes: the estimator, its formula as coded, the modelling assumptions, the MINFLUX-specific
  handling, and how to read the output. This is the "how it works" layer, extracted from the source.
- **Part II — Literature & Appropriateness ([§1](#1-why-minflux-data-needs-its-own-treatment) onward).**
  Which researcher each method is ported from (so the work is cited), and an honest **appropriateness
  rating** of each method against the statistical character of MINFLUX data.

The two layers cross-reference each other: every Part I entry names its Part II rating, and vice versa.

Part II ratings:

- 🟢 **GREEN** — well-suited; MINFLUX's properties are neutral or an advantage.
- 🟡 **AMBER** — usable, but an assumption is strained; needs a MINFLUX-aware adaptation.
- 🔴 **RED** — a core assumption is genuinely violated; results are biased unless mitigated.

> A rating describes the *textbook method*. Where our implementation already adapts the method for
> MINFLUX, that is called out under **Our handling** — several tools are more defensible than their
> raw textbook rating suggests.

---

# Part I — Implementation & Mathematics

*What each tool computes, as coded. Formulas are transcribed from the source (`server_tracking.R`,
`server_smlm.R`, `helpers_smlm.R`, `sim_server.R`); each entry names its Part II appropriateness rating.
Units in code: positions in **nm**, time internally in **ms** (imported `tim` is seconds → ×1000). D is
carried internally in **nm²/ms** and displayed in **µm²/s**, with `D_um2_s = D_nm2ms / 1e3`
(1 µm²/s = 10³ nm²/ms).*

Conventions used throughout:
- `d` = spatial dimensionality (2 or 3); the MSD prefactor is `2d` (4 in 2D, 6 in 3D).
- Free Brownian MSD with static localization error σ: **MSD(τ) = 2d·D·τ + 2d·σ²** (blur term R≈0 for
  MINFLUX, see §5 action item 2).
- Expected squared step over a gap Δt: **E[|Δr|²] = 2d·D·Δt + 2d·σ²**.

## I-A. SMLM (imaging) analysis

The SMLM tab treats each acquisition as a spatial point pattern of localizations, optionally reduced
through a **de-blinking hierarchy** before any spatial statistic:

> **localizations → trace (`tid`) median → primary cluster (DBSCAN `cluster`) median → secondary
> structure (`cluster2`)**

Collapsing to trace/cluster medians is the mitigation for the labelling/blinking over-counting bias
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
  et al. 2003; Dunn et al. 2011) rather than a pixel-intensity coefficient, which is the correct family
  for MINFLUX: there is no intensity image, only coordinates.
- **Caveat (no null):** cross-NN has **no null model built in** — every A point has a nearest B point, so
  the curve is never zero and a small median distance can arise from density alone. Significance must
  come from L₁₂/g₁₂ with their CSR envelopes (§I-A.4).
- **Attribution precision (tightened 2026-07-30):** this panel was credited to **Clark & Evans 1954**
  as though it implemented their statistic. NN distance *as a spatial statistic* does originate with them,
  but their contribution is specifically the index **R = d̄_obs / d̄_Poisson** tested against a Poisson
  null — and this panel deliberately computes no null (see above). The UI now says so explicitly rather
  than implying a Clark–Evans test is being run.

### I-A.2 Density-based clustering (DBSCAN / HDBSCAN / OPTICS) — [Part II §6.1]
- **Primary clustering** (`cluster`): `dbscan::dbscan(pts, eps, minPts)` on raw localizations → labels
  blinks/dyes of one region together. **Median clustering** (`cluster2`): the same family run on
  trace/cluster medians → groups molecules into structures. Two-colour variants cluster each channel
  separately (channel-2 ids offset by +10000) then lock every TID to its majority cluster.
- **eps tuning:** an optional **kNN-distance plot** (sorted distance to the k-th neighbour) exposes the
  density "knee" that DBSCAN's `eps` should sit at.
- **Density field:** `compute_spatial_density()` uses `dbscan::frNN(pts, eps = r)` to count neighbours
  within radius `r` per point (fixed-radius local density used for colouring and the acquisition
  structure/background split).

### I-A.3 Ripley's K / L(r) — [Part II §6.2]
`compute_ripley_L(x,y,z, r_seq, use_3d)`:
- Random **subsample to ≤2000 points** for tractable O(n²) pairwise distances (`stats::dist`).
- Observation window = **axis-aligned bounding box**; `area_vol` = box area (2D) or volume (3D).
- For each radius r, count ordered pairs within r and form the unbiased estimator
  **K(r) = area · 2·#{d_ij ≤ r} / (n(n−1))** (the ×2 restores both directions from the lower triangle).
- **L-transform:** `L(r) = √(K/π)` (2D) or `L(r) = (3K/4π)^{1/3}` (3D). Plotted as **L(r) − r**: 0 under
  complete spatial randomness (CSR), **>0 clustering**, **<0 dispersion**.
- **Caveat (implementation):** there is **no analytic edge correction** and the window is the bounding
  box, so L(r) − r is biased downward near large r (points near the box edge lose neighbours). A
  `# TODO` remains for a translation/isotropic correction; instead this is handled by the envelope
  below.
- **Auto r_max is capped by the SHORTEST window axis** (`min(extent)/4`), since edge bias is set by
  whichever dimension runs out first. **Fixed 2026-07-30:** an inverted condition (`!use_3d` where
  `use_3d` was meant) meant the **z** extent never entered the cap, so on a typical thin 3D slab
  (5000×5000×300 nm) auto r_max came out at 1250 nm — over 4× the entire z depth — instead of 75 nm,
  making the whole large-r tail of the 3D curve edge-artefact. The same "cap by the longest axis rather
  than the shortest" error was present in `compute_ripley_L12` and `compute_pcf12` (`max(rng)/4`) and is
  likewise fixed. Only the **auto** path was affected: an explicit r_max from the UI was always honoured.
- **Monte-Carlo CSR envelope (`compute_ripley_envelope`, optional):** generates `nsim` CSR patterns of
  the *same n* in the **same frozen window** (`window_ranges`, a pre-supplied 2×d bounding box shared by
  the observed curve and the sims) and returns the pointwise 2.5–97.5% band + mean null curve. Because
  the sims inherit the identical window, the edge bias is present in the band too and **cancels when the
  observed curve is read against it** — so the envelope, not the flat zero line, is the interpretive
  reference (the flat line is relabelled "theoretical CSR (uncorrected)"). A two-sided α≈0.05 test needs
  `nsim ≥ 39`; default 99, capped 19–199. In all-files mode the envelope is computed for the first few
  files only (cost = `nsim × files × O(n²)`).

### I-A.4 Bivariate Ripley L₁₂(r) & pair cross-correlation g₁₂(r) — [Part II §6.2]
- **L₁₂:** cross-type K from the n₁×n₂ distance matrix, **K₁₂(r) = area · #{d ≤ r} / (n₁·n₂)**, same
  L-transform. Measures co-clustering of two channels (>0 attraction, <0 segregation).
- **g₁₂ (PCF):** the *derivative* form — count cross pairs in the **annulus/shell** (r−dr/2, r+dr/2],
  normalise by shell area (2D `π(hi²−lo²)`) or volume (3D `4/3·π(hi³−lo³)`):
  **g₁₂(r) = area · #{annulus} / (n₁·n₂·shell)**. **g₁₂ = 1 under CSR**, >1 attraction, <1 repulsion.
  PCF is scale-local (not cumulative like K), so it localises the interaction distance more sharply.
- Same ≤2000 subsample + bounding-box + no-edge-correction caveat as L(r).
- **CSR envelope:** the same `compute_ripley_envelope` (with `estimator = "L12"` / `"pcf"` and both `n₁`,
  `n₂`) randomizes **both channels independently** in the shared window and shades the 2.5–97.5% band.
  For L₁₂ the observed curve above the band = co-localization at that scale, below = exclusion; for g₁₂
  the band **brackets g = 1** under CSR and the observed curve rising above it at short r = attraction.
  The bias cancels in the band exactly as for the univariate L(r).

### I-A.5 KDE density map + Pearson co-localization — [Part II §6.4]
- `MASS::kde2d` on each channel over a **shared grid** with a common bandwidth.
- **Bandwidth (auto):** Silverman's rule `bw = 0.9·min(sd, IQR/1.34)·n^{-1/5}`.
- **Co-localization index:** Pearson r of the two flattened density surfaces. **Explicitly bandwidth-
  dependent** — small bw → r≈0 even for real co-localization, large bw → inflated; treat as a *relative*
  measure at a fixed bandwidth (stated in-UI).
- **Bandwidth-honesty add-ons:**
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
    - **`"csr"`: uniform resampling** of Ch2 — a test against *complete spatial randomness* rather than
      against independent placement. Kept for comparison and labelled in-UI.
  - The UI points users to **g₁₂(r)** for a bandwidth-free, scale-resolved alternative.
- **Null changed + CALIBRATION MEASURED (2026-07-30).** The old uniform null tests against CSR, not
  against "Ch2 has this cell's structure but is placed independently of Ch1" — the failure mode Costes
  et al. 2004 designed *against* (their test preserves the randomized channel's autocorrelation by block
  randomization), so the field is also no longer attributed to Costes as if it implemented their
  algorithm. The default is now the structure-preserving random shift. **But the fix is partial, and the
  measured numbers matter more than the argument.** Simulating independent placement 30× per scenario
  (truth = NOT co-localized, so a calibrated test rejects at p<0.05 ~5% of the time):

  | scenario | random shift | uniform CSR |
  |---|---|---|
  | channels clustered across the **whole window** | 3% ✅ | 3% ✅ |
  | channels sharing **one cell** inside a larger window (realistic SMLM) | **53%** ❌ | **93%** ❌ |

  So the CSR null genuinely is severely anti-conservative in the realistic case — but the random shift is
  only *better*, not calibrated. Root cause: **both nulls randomize over the bounding box rather than over
  the actual support** (the cell mask), so a shifted Ch2 wraps outside the cell and the null still
  contrasts "Ch2 inside the cell" with "Ch2 partly outside it". **Mask-aware randomization is the real
  fix and is not implemented.** Until then, read a low p here as *weak* evidence whenever the signal does
  not fill the field, and take significance from g₁₂(r) against its CSR envelope. Power is not the
  limitation: on a genuinely co-localized positive control the shift null fired 10/10 at median p = 0.02.

### I-A.6 Radius of gyration & gyration-tensor anisotropy — [Part II §6.3]
- **Rg:** `Rg² = mean(‖p − p̄‖²) · n/(n−1)` per cluster (second moment about the centroid, Bessel-
  corrected). Points are **median member positions**, never raw localizations.
- **Small-N bias correction (added 2026-07-30).** The centroid is estimated from the same points whose
  spread is measured, so it sits inside the cloud and **E[Rg²_obs] = Rg²·(n−1)/n** — a 29% underestimate
  at n=2, 11% at n=5, 1% at n=50. That bias is **n-dependent**, so uncorrected it manufactures a spurious
  "larger clusters have larger Rg" trend whenever clusters of differing member count are pooled or plotted
  against size. The `n/(n−1)` factor (exactly a sample-variance Bessel correction) removes it; verified in
  silico, recovering 0.98–1.00 of the true Rg for n≥5 versus 0.87–0.99 uncorrected. **Residual:** at n=2
  the corrected mean still sits ~12% low, because Rg = √(Rg²) and E[√x] < √E[x] (Jensen) — no closed-form
  fix, so treat 2-member clusters as indicative only. `bias_correct = FALSE` returns the raw
  configurational Rg of the observed point set.
- **Anisotropy:** eigenvalues λ₁≥λ₂≥λ₃ of the gyration tensor give **Westin linear / planar / spherical**
  coefficients (c_l = (λ₁−λ₂)/Σλ, c_p = 2(λ₂−λ₃)/Σλ, c_s = 3λ₃/Σλ) and the **Rudnick–Gaspari asphericity**
  — used to *infer* cluster shape (linear filament vs planar vs isotropic) without the user declaring it.

### I-A.7 Cluster shape fitting (3D viewer) — [Part II §6.5, 🟠 EXPLORATORY]
- **Dictionary:** sphere, ellipsoid, helix, ring, arc. Each is fit by **orthogonal-distance MLE**
  (minimise Σ of squared point-to-surface distances). The ellipsoid projection solves the exact
  nearest-point KKT condition by vectorised Newton iteration in the PCA frame
  (`mle_project_ellipsoid`), with a radial fallback if Newton diverges.
- **Model selection** (`viewer_rank` = AICc | BIC | RMSD): under an iid isotropic-Gaussian residual model,
  σ̂² = RMSD², so **−2logL = n·log(2πσ̂²) + n**; **AIC = −2logL + 2k**, **AICc = AIC + 2k(k+1)/(n−k−1)**,
  **BIC = −2logL + k·log n**, with k = (geometric DOF + 1). AICc is the default (penalises complex shapes
  hard at small n); it → ∞ when n ≤ k+1 (too few points), which correctly refuses to over-fit tiny clusters.
- **Input units:** fits run on **de-blinked trace centroids**, not raw localizations — otherwise the many
  localizations of one long trace dominate the residual sum and the fit describes the brightest molecule
  rather than the structure. The 2D view is a PCA projection onto PC1×PC2 (Pearson 1901).
- **Why exploratory:** the shape is imposed by the dictionary and "best" is a model-selection outcome, not
  proof of biology; with few blink-inflated points several shapes are statistically indistinguishable.
- **Photon (`eco`) weighting — implemented, `viewer_weight_eco`, default on.** The unweighted fit assumes
  **iid isotropic Gaussian** residuals with a single σ̂, but MINFLUX precision scales as 1/photons. A
  trace centroid averages that trace's localizations, so its variance goes as 1/Σ`eco` — meaning the
  trace's **total photon count is exactly its inverse-variance weight**. The five `mle_fit_*` optimisers
  already accepted a `w` argument (`sum(w·(L−r)²)`); the viewer now supplies it, and the reported RMSD
  becomes the weighted √(Σw·d²/Σw).
  **Effect on the criteria:** with per-point variance σ²/wᵢ the log-likelihood gains a −Σlog(wᵢ) term that
  is *identical for every shape*, so it cancels in all AIC/AICc/BIC **differences** — model ranking stays
  valid. Only the absolute values shift, so **never compare criterion values across weighting modes**.
  **`n_eff` (Kish effective sample size, `(Σw)²/Σw²`) is reported in the viewer header** because
  weighting is not free: it trades bias for variance.
  ⚠ **Filter double events before weighting — `eco` is dual-signed.** A high `eco` means better
  precision, *but it is equally the signature of two fluorophores blinking in the same spot*, which is
  implausible for a single dye. Inverse-variance weighting hands the most influence to the brightest
  traces, i.e. **on unfiltered data it preferentially amplifies double-emitter artefacts** — the
  weighting is only sound downstream of the `eco`/`cfr` metadata filters that remove implausibly bright
  events. This is a real limit on the method, not a caveat about the implementation: the same statistic
  that makes a trace informative also makes it suspect.
  ⚠ *Measured on the collaborator test file (unfiltered):* per-trace photon budgets span **50 → 13 174 (263×)**,
  and on the largest structure weighting moved the fitted sphere radius from **359 nm to 149 nm** while
  `n_eff` fell from 170 to 83. On a **filtered** dataset the two fits agree closely (user-observed),
  which is the expected behaviour and supports the reading above: most of that 2.4× swing was
  double-event contamination being up-weighted, not a genuine precision correction. Practical rule —
  **a large weighted-vs-unweighted difference is a data-quality alarm, not a better answer.** Quote the
  weighting mode and `n_eff` with any shape parameter; treat AICc differences of a few units as
  indistinguishable either way.

### I-A.8 Ground-truth position/count estimator — [Part II §6.6, 🟠 EXPLORATORY]
`estimate_gt_from_labeling(df, preset, eff)` — a **labelling-model inverse** (`helpers_smlm.R`):
- **Label preset** parameters: `n_fluor` (mean dyes/molecule), `link_nm` (linker length), `dye_spread_nm`
  (dye scatter σ within the complex), `footprint_nm`, `blink_lambda` (extra re-appearances), `eff`
  (labelling efficiency).
- **Forward model it inverts:** molecule → Poisson(`n_fluor`) fluorophores, each displaced by a
  half-normal linker (mean ≈ `link_nm`) and scattered by `dye_spread_nm`, each blinking into
  ~(1+`blink_lambda`) traces. Expected traces per molecule: **m_tid = n_fluor · eff · (1+blink_lambda)**;
  detection efficiency **η = eff** (single dye) or **η = 1 − e^{−n_fluor·eff}** (multivalent).
- **Inverse:** build candidate units (trace or primary-cluster medians, ≥`minpts` locs), then **merge**
  unit centroids within a **complex radius** `complex_r` (auto = max(footprint, 2.5·dye_spread)) via
  `dbscan(eps = complex_r, minPts = 1)` so several dyes/blinks of one molecule collapse to **one candidate
  GT position** (this produces the contribution fan-out lines). Optional **stoichiometry-guided merge**
  further collapses to k = round(ΣN_tid / m_tid) via hierarchical clustering — only ever *reduces* toward
  the labelling-predicted count, never splits.
- **Output:** candidate GT positions coloured **green→red by positional uncertainty** (`pos_sd_nm`), a GT
  **count range** [lo, hi] rather than a single number, and a **crowding flag** when spacing < label
  footprint (under-resolved). Validated *in silico* against the simulator; quantitatively a heuristic.

### I-A.9 Detection-over-time / crowding / NN-vs-size — [Part II §6.7, 🟠 EXPLORATORY]
Bespoke descriptive diagnostics (detection-rate vs time, crowding vs timing, NN distance vs cluster size).
**Correlational, not causal** — detection timing and crowding co-vary with photobleaching/blinking/labelling,
so a trend is a hypothesis for controls, not a measurement.

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
- **Crowding-vs-timing has a built-in positive trend — now corrected by a permutation null.** Denser
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
  *Validated on the collaborator test file* (472 traces, 20 qualifying structures): null-band coverage 0.052
  against a nominal 0.05, null ρ centred on 0, and band width vs cluster size ρ = −0.92 (small clusters
  correctly get wide intervals). Observed on that file: ρ = −0.23, p ≈ 0.40 — **no significant
  crowding–timing trend**.
  *What it still cannot do:* rule out a third variable driving both axes (labelling density, structure
  identity). That needs a control experiment, not a null.
- **NN-vs-size is censored by the clustering algorithm.** Under DBSCAN every within-cluster NN distance is
  bounded above by the link radius `eps2` (default 200 nm), so the y-axis is **right-censored**: a curve
  that flattens at large cluster size may be hitting the ceiling rather than showing constant spacing.
  HDBSCAN/OPTICS impose no such bound — re-run under those before reading a slope. Convex-hull size is
  additionally upward-biased by single outlying members.

### I-A.10 Traces per cluster (blinking-event count) — [Part II §6.1a, 🟡 AMBER]
Counts **distinct TIDs** per primary or median cluster. Counting traces rather than localizations already
collapses the many localizations MINFLUX collects within one trace, so the number is a *blinking-event*
count. It is **not stoichiometry**: one fluorophore re-enters the on-state repeatedly and is re-assigned a
new TID each time, and a multivalent antibody carries several dyes — so the count over-estimates molecules
by an unknown, dye- and buffer-dependent factor. The Annibale et al. (2011) blink-correction that would
divide this out requires a **dark-time distribution that is not calibrated in this app**. Valid as a
relative comparison between conditions acquired identically; never as an absolute subunit count.

## I-B. Tracking (SPT) analysis

### I-B.0 Cycle time & Time-Signature Analysis (the foundation for everything Δt-aware)
The **Time Signature Analysis (TSA)** tab is a faithful R port of **Bela Vogler's MINFLUX Time-Signature
Analysis (MFX-TSA)** — §2 of his `the_road_to_minflux_msd` tutorial (Vogler et al. 2025, *Commun. Biol.*;
reference implementation *Toolbox_for_MINFLUX_SPT / MinSpt*). The correspondence is direct:
- **cycle time / hardware parameters** — MinSpt `initialize_tsa` + `estimate_time_signal_parameter`
  (EOD/SW overhead, pattern vertices, dwell) ↔ our `trk_cycle_time` + the EOD/SW regression and the
  hardware pattern (hexagon 6 / square 4 / triangle 3 vertices).
- **integrated vs dark cycles + per-track progression** — `print_cycle_stats` /
  `plot_cycle_progression_for_ID` ↔ our Integrated/Dark distribution and per-track cycle-progression plot.
- **adaptive gap removal** — `suggest_split_thresholds` (TOT-based `med_th = floor(median+std+1)`, capped
  by `high_th = maxOffTime/cycle`) + `remove_time_gaps` (`min_len = 50`) ↔ our adaptive gap-split.
- **linearized cycle-index timeline** — `prepare_for_msd` ↔ our cumulative cycle-index binning.

MINFLUX has no fixed clock; the base **cycle time t_c** is inferred (`trk_cycle_time`), in priority order:
(1) hardware model `t_c = N_iter · t_dwell` if enabled; (2) manual override; (3) **data-driven cluster
cycle** — cluster the inter-localization Δt values and take the smallest robust cluster centre (the modal
Δt is *misleading* because dark cycles inflate it); (4) a regression cycle from Δt-vs-index slope. The
resulting t_c drives gap-splitting, the quality gate, cycle-index MSD binning and the D-map. A step's
**cycle count** is `n_cyc = round(Δt / t_c)`; **integrated (single-cycle)** steps have n_cyc = 1.

#### Measured on real data (2026-07-30) — Δt is QUANTISED, not distributed
Checked directly against `Test data/251027-202407_GP67H-track-014.npy` (650,042 within-track steps,
`tim` stored as float64 so this is real structure, not a storage artefact):

- **Δt = k·t_c + SW, essentially exactly.** Only 4,432 distinct Δt values across 650k steps. Regressing
  the interval clusters gives **t_c = 0.19161 ms, SW = 0.01153 ms (R² = 0.9999978)**. Testing whether
  intervals are integer multiples: dividing by the *mode* succeeds for only **52%** of steps, while
  `(Δt − SW)/t_c` is within 1% of an integer for **97.6%** and within 5% for **99.98%**.
- **⚠ The modal Δt is NOT the cycle time** — it is `1·t_c + SW` (0.2027 ms here). Mistaking the mode for
  t_c would bias every D by **−5.5%**. This is exactly why `trk_cluster_cycle` takes the cycle as the
  **spacing between interval clusters** (`median(diff(centers))`) and SW as the leftover offset
  (`centers[1] − cyc`) rather than using the mode. Verified end-to-end: the app's own estimator returns
  **t_c = 0.191483, SW = 0.011165**, agreeing with the independent regression to **0.065%** on t_c.
- **There is no meaningful jitter within a cycle.** The integrated (n_cyc = 1) population has
  **CV = 0.0425%**, a total spread of 0.72% of one cycle, ~19 distinct values separated by 25–200 ns
  timestamp ticks, and **46% of steps on a single identical value**. Each "peak" is a near-delta spike.
- **Consequence for the TSA Δt histogram binning (2026-07-30).** Fixed-width bins alias against the
  comb. At the old default of 80 bins the grid ran **1.63 bins per cycle**, so each peak landed at a
  different offset inside its bin and the bar centres sat **5–25% of a cycle away from the true peaks**,
  drifting systematically across the range (+25%, −23%, −19%, −14%, −10%, −5% for k = 1…6) with empty
  bins interleaved — which is what made the plot look misaligned. Binning is now always
  **one bar per cycle count**: edges on the classification boundaries `SW + (k ± ½)·t_c`, so each bar is
  exactly one cycle count centred on its peak (verified: bar centres match the peaks to 0.0e+00 ms) and
  the histogram is literally a picture of the `round((Δt − SW)/t_c)` the analysis performs. Fixed-width
  binning survives only as an internal fallback when no cycle time is available. A **log count axis**
  option is offered, since counts span ~100× from the 1-cycle peak into the dark tail.
- **Consequence for the TSA Δt plot: the KDE overlay was REMOVED (2026-07-30).** It had been described
  in-UI as revealing a "bell-shaped cluster around each cycle count" — but no such bell exists. Since
  each peak is a near-delta spike, 100% of the curve's width was kernel bandwidth rather than timing
  jitter, and at higher smoothing settings it merged the comb into a single invented unimodal blob. The
  `trk_dt_show_density` and `trk_dt_density_bw` controls are gone; the histogram plus the
  **cycle-multiple markers** (`trk_dt_mark_cycles`, now a standalone option rather than gated behind the
  density checkbox) convey the structure without fabricating a continuous density. The real distribution
  is the **discrete** one over cycle counts (here 52.4% integrated, 21.6% 2-cycle, 8.2% 3-cycle, with a
  2.8% tail beyond 12), shown by the cycle-classification plot.
- **This validates the cycle-index design.** Because intervals are exact integer cycles plus a constant
  offset, the cumulative-cycle-index binning is not an approximation to noisy timing — it recovers the
  true elapsed time, and the √Δt / integrated-only handling elsewhere rests on solid ground.

### I-B.1 σ_loc & data-quality calibration (CVE, MSD-free) — [Part II §2.2]
`trk_cve_bundle`. Consecutive per-axis steps share one localization's error and are anti-correlated by −σ²
(Vestergaard 2014). With MINFLUX blur R≈0:
- **σ_loc² = −⟨Δx_i · Δx_{i+1}⟩**, averaged over x, y[, z], using **only consecutive integrated
  (single-cycle) step pairs** so variable dark-cycle Δt cannot bias the covariance (the constant-Δt
  restriction the literature prescribes).
- **Free bias-corrected D:** `D_cve = ⟨Δ²⟩/(2Δt) + cov₁/Δt` per axis (the +cov₁ term removes the
  localization-error inflation).
- **Per-track σ heterogeneity:** one CVE σ per track (≥ `n_min` integrated steps, skipping the first
  ~5 locs of lock-on transient); the median, IQR, CV and a **σ–EFO Spearman ρ** (precision should improve
  with photons ⇒ negative ρ validates the model) plus a bootstrap 95% CI on the median σ.
- **Photon-quality weighting: REMOVED (was added 2026-07-19, withdrawn since).** The doc previously
  described an inverse-variance CVE weighted by a fitted `log σ ~ log EFO` model as "default on"; that no
  longer exists in the code and this entry is corrected. It was withdrawn because on MINFLUX **tracking**
  data the photon channel is structurally coupled to **mobility**, not just brightness: a fast molecule
  leaves the beam pattern sooner and so collects *fewer* photons per localization. Weighting by photons
  therefore up-weighted slow molecules and biased D downward — the estimator was measuring the very thing
  it was meant to correct for. The pooled σ_loc / D_cve never used the photon model and are unaffected.
  The two diagnostics that exposed it are retained and reported:
  - **`noise_frac` = −2·r₁** (r₁ = cov₁/⟨Δ²⟩), the share of step variance that is localization noise.
    Small ⇒ per-track σ differences are motion, not photons.
  - **`phot_mob_rho`** = Spearman ρ[per-track photons, per-track step variance]. Strongly negative ⇒ the
    photon channel is a mobility proxy, which also means **any ECO/EFO filter is silently a mobility
    filter** — surfaced as a UI warning at ρ ≤ −0.35.
- **⚠ Placeholder z columns were being treated as 3D (fixed 2026-07-30).** 2D MINFLUX exports routinely
  carry an all-zero `z` column, and `0` *is* finite — so the widespread presence test
  `sum(is.finite(df$z)) > 0` returned TRUE for 2D data and the dimension-aware estimators switched to 3D,
  pooling a **zero-variance axis** into the maths. Measured effect on 2D data with such a column: pooled
  **σ_loc −18.4%** (a factor √(2/3)) and **D_cve −33.3%** (a factor 2/3). It would also have made the
  now-d-correct JDD fit a Maxwell law to 2D distances, and built a 3D free-walk null for 2D confinement
  data. There is now a single shared predicate **`.has_real_z(df, min_sd = 1)`** (>10 finite z values
  *and* sd > 1 nm) used by the 3D toggle, CVE, MSD, JDD, confinement and displacement-vs-D. Previously
  MSD and displacement-vs-D tested variance correctly while CVE, JDD and confinement did not — i.e. the
  dimension decision had silently drifted between tools, which is exactly what the shared predicate
  prevents. **Re-running CVE on such a dataset raises σ_loc by ~22% and D_cve by ~50%.**
- **Per-axis σ: lateral vs axial (added 2026-07-30).** The pooled σ_loc above averages x, y **and z**.
  In 3D MINFLUX the axial precision is typically worse than the lateral, so one pooled number over-states
  σ_xy and under-states σ_z. Several consumers are purely lateral — the **D-map** is an x–y map and the
  **HMM** never reads z — so feeding them a z-inflated floor makes them **over-subtract** noise. `sigma_xy`,
  `sigma_z` and their ratio `sig_aniso` are therefore computed separately and reported; at
  `sig_aniso > 1.3` the panel warns and gives the lateral value to enter manually for the 2D tools. The
  app-wide σ_loc default is left pooled (changing it silently would move every existing result).
- **QC flags:** *implausible jumps* where `|Δr|² > 25·E[|Δr|²]` (i.e. |Δr| > 5× the expected RMS for that
  gap, using E[|Δr|²] = 2d·D_cve·Δt + 2d·σ² so long dark gaps aren't falsely flagged); *degenerate* tracks
  with σ² ≤ 0 (drift-dominated) or σ far in the outlier tail.
- **Continuity score:** tolerance-weighted `w = 2^{−(n_cyc−1)/H}` averaged over steps (H = dark-cycle
  tolerance), so a dataset of gentle 1–2-cycle gaps scores higher than its bare integrated fraction.
- σ_loc from here becomes the **app-wide noise floor** consumed by the D-map, JDD, confinement, HMM,
  directionality and speed tools.

### I-B.2 Unified MSD, optimal-fit OLSF & fit models — [Part II §2.1, §2.7, §3a]
`msd_unified_bundle` + `msd_unified_fit` + `.msd_model_fit`.
- **Binning:** per-track **cumulative cycle-index** lag (Bela's TOT) — `dcyc = round((Δt − sw)/t_c)`,
  `cycidx = cumsum(dcyc)` — so pairs are grouped by *integer elapsed cycles* (true elapsed time), removing
  the `round(total/t_c)` sawtooth. Per-track TAMSD = Σdr²/Σpairs; x-axis τ = mean real Δt per bin.
- **Ensemble:** **pair-weighted** (default, smoother) or **track-weighted** eMSD = mean of per-track TAMSDs
  (Vogler). Optional **D>0 filter** drops tracks whose TAMSD has a negative OLS slope (Bela's rule); IQR
  band from per-track quantiles; bins kept only where ≥10% of tracks contribute.
- **Linear fit — OLSF (Michalet & Berglund optimal point count):** unweighted linear fit
  `MSD = a + b·τ`, `D = b/(2d)`, `σ² = a/(2d)`, with **separate optimal point counts** for the slope
  (→D) and intercept (→σ²). Reduced localization error `x = σ²/(D·Δt) − 2R`; optimal p = `f·L/(f³+L³)^{1/3}`
  with slope `f_b = 2 + 1.35·x^{0.6}, L_b = 0.8 + 0.564·N` and intercept `f_a = 2 + 1.6·x^{0.51},
  L_a = 3 + (4.5·N^{0.4} − 8.5)^{1.2}` (R = 0 for MINFLUX). Two variants: **anchored** (default — anchors x
  on a stable whole-window fit so p is reproducible) and **iterative** (Bela's exact loop — converges onto
  the lowest, least-noisy lags). Fit window default **5–50 cycles** (Vogler `fit_limits`).
- **Non-linear models** (whole window — they need curvature):
  - **Confined:** `MSD = A(1 − e^{−τ/t_c}) + C` → local `D = A/(2d·t_c)`, confinement size
    **L_c = √(6A/d)**, plateau A; flagged *degenerate* if t_c > 2·max(τ) (no plateau resolved).
    **L_c is fitted from the data, never entered** — it is read off the plateau A, so it is only as
    meaningful as the plateau is resolved (hence the *degenerate* flag) and is additionally warned
    about when L_c < 2σ_loc.
    **Dimensionality fix (2026-07-30):** this was `√(3A)` for all d. A particle uniformly exploring a
    d-dimensional box of side L_c has Var(Δx) = L_c²/6 per axis, so A = d·L_c²/6 and L_c = √(6A/d) —
    which reduces to √(3A) in 2D but **over-estimated L_c by √(3/2) = 1.22× in 3D**. Verified against
    reflecting-boundary simulations and by round-tripping ideal curves through `.msd_model_fit`
    (L_c recovered to <0.01% at d = 2 and 3, for L_c = 140/300/800 nm).
    **Symbol note:** the confinement size is written **L_c** and the MINFLUX targeted-coordinate-pattern
    diameter **L_TCP** (D_max panel). Both used to be plain "L", which read as one quantity.
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
level (clamped ≥0, not per-step — that would rectify and bias upward), the same floor the D-map uses, so
truly immobile segments read ~0 instead of a small non-zero D (validated in silico). The pooled
per-segment log₁₀D distribution is then fit with a **BIC-selected Gaussian mixture** (`mclust`, up to
`max_g` states).

### I-B.4 Jump-distance distribution (JDD) — [Part II §2.3]
`jdd_bundle`. Fit the step-length histogram with an **N-component scaled-chi mixture** — the law of
|Δr| for a d-dimensional isotropic Gaussian displacement:

**p(r) = Σ_k w_k · C_d · r^(d−1)/σ_k^d · exp(−r²/2σ_k²)**,  **C_d = 1/(2^(d/2−1)·Γ(d/2))**

by MLE (Nelder–Mead on the NLL; weights via softmax, σ via log-parametrisation). d = 2 gives the
classic **Rayleigh** form (C=1, Schütz 1997); d = 3 gives **Maxwell** (C=√(2/π)), used when z is
included. σ_k is the **per-axis** spread in both cases. MINFLUX handling:
- **√Δt normalization ON by default** — pool `r/√Δt` so steps of different cycle counts are comparable
  (states then map to D, not raw step size).
- **D per component:** `D = σ²/2`, **dimension-independent** — because E[r²] = d·σ² and MSD = 2·d·D·Δt,
  the d cancels. This also makes the σ_loc term below dimensionally consistent with it.
- **σ_loc bias correction:** subtract `2·σ_loc²/dt0` from σ² before D (localization noise adds 2σ_loc² per
  coordinate), clamped ≥0. A **noise-floor jump distance** `√(2d)·σ_loc/√dt0` marks where peaks are
  noise-dominated.

> **Fixed 2026-07-30.** This previously fitted the 2D **Rayleigh** shape even when z was included, then
> divided by 3. Because a Rayleigh MLE is a moment estimator on Σr², the *total* D happened to come out
> right for a single noiseless population — but the fitted σ was inflated by √1.5, so the σ_loc term
> (defined on the per-axis variance) subtracted only ⅔ of the true noise floor, and the PDF *shape* was
> wrong. Verified in silico at realistic MINFLUX values (D = 200 nm²/ms, t_c = 0.5 ms, σ_loc = 10 nm):
> single 3D population D biased **+32%**; a 100/900 nm²/ms two-population 3D mixture recovered as
> **281/840** with weights 0.63/0.37. After the fix: −1.0% and **102/903** at 0.60/0.40. The 2D path is
> bit-identical (d = 2 reproduces the old expression exactly). Note the **HMM** and **D-map** are
> deliberately 2D-only (they never read z), so the Rayleigh emission there is correct — and since the 2D
> projection of isotropic 3D diffusion has the same D, their D is unbiased on 3D data.

### I-B.5 Diffusive-state HMM (Baum–Welch + Viterbi) — [Part II §2.4]
`hmm_bundle`. Per-step displacement `r` (√Δt-normalised by default) with a **Rayleigh emission**
`log p(r|state) = log r − 2log σ_k − r²/2σ_k²`, K = 2–4 states. **Baum–Welch EM** entirely in log-space
(vectorised forward/backward, log-sum-exp), ≤30 iterations with early-stop; **σ_k update =
√(Σγ·r²/(2·Σγ))** (Rayleigh MLE), transition matrix A from the accumulated log-ξ. Fit on the longest
`max_fit_trk` tracks (optional contiguous-prefix step cap for speed); **Viterbi** decodes *all* tracks at
full length. States sorted slow→fast by σ.
**Continuous-time transition rates (2026-07):** the fitted per-step matrix A conflates switching with
each step's irregular Δt (a 3-dark-cycle step counts as one transition opportunity, like a 1-cycle
step), so per-step probabilities are *not* rates. We therefore also estimate the **generator Q** of the
underlying Markov jump process from the decoded state path using *real elapsed time*: the
continuously-observed-path MLE `q[a,b] = N[a,b]/T[a]` (b≠a; N = decoded a→b jumps, T = time spent in a),
`q[a,a] = −Σq[a,b]`, giving rates in s⁻¹ and **mean dwell times** τ_a = −1/q[a,a] that are invariant to
dark cycles. Intervals spanning > `gap_tol_cyc` (=4) cycles are excluded from the counts so one big gap
can't masquerade as a long dwell. Recovers injected rates within ~15% in silico. A is retained (heatmap)
for reference; state σ/D remain the √Δt-robust emission quantities.

### I-B.6 Confinement (Rg / Rmax) — [Part II §2.5]
`conf_bundle`. Per track: **Rg² = mean(‖p−p̄‖²)**, **Rmax = max‖p−p̄‖**, ratio **Rg/Rmax**. Optional
noise correction removes localization inflation: **Rg²_true = Rg²_meas − d·σ_loc²·(1 − 1/n)** (the 1−1/n
accounts for estimating the centroid from the same noisy points), clamped ≥0. **Rmax is now given the same
first-order subtraction** (**Rmax²_true = Rmax²_meas − d·σ_loc²·(1 − 1/n)**, clamped ≥0): this removes the
*average* noise inflation of the extreme point. The residual **extreme-value selection bias** (the max is
picked partly because its noise was positive) has no closed form and is left in place, so Rmax stays
mildly conservative on short tracks (validated in silico: an immobile track's Rmax shrinks substantially
but not to zero).

**Confinement call — N-matched free-diffusion null (2026-07, corrected direction).** The raw ratio is
track-length-dependent, so a single fixed threshold is unreliable. `.free_walk_ratio_dist()` Monte-Carlos
free Gaussian walks of the SAME localization count (the ratio is scale-free, so D/Δt drop out) on an
N-grid, and each track is classified against that null: **confined = ratio above the free 95th pct**,
**directed/escape = below the 5th pct**, else free-like. **Direction correction:** diffusion simulations
show a bounded (reflected-boundary) walk has a *higher* Rg/Rmax than free (≈0.60 vs 0.53) — it fills its
domain while free diffusion wanders to extreme excursions. The old fixed thresholds *"confined <0.35 /
free >0.45"* were therefore both mis-calibrated (below every regime, ~0.53–0.61) and pointed the wrong
way; they are removed in favour of the null. Validated in silico: free tracks classify ~5% confined (=
the null false-positive rate by construction), confined tracks are ~3× enriched, directed tracks are not
called confined. Residuals (both conservative): the Rmax extreme-value bias above; and a genuinely
confined but very slow molecule can be missed if it barely explores. Corroborate with the τ-based MSD
confinement index / plateau.

#### Effective dimensionality of the MSD (added 2026-08-03)
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

The last row is the subtle one, and it is the case the user raised: **the cell is three-dimensional
even when tracking is not.** A molecule on a curved surface steps in its local tangent plane, and how
much of that step survives projection into the image plane depends on how the surface is tilted there.
Averaging over the shape gives ⟨n_z²⟩ = 1/3 for a sphere, 1/2 for a cylinder whose axis lies in the
image plane, and 1 for a flat membrane parallel to it — hence **d_eff = 4/3 (sphere), 3/2 (rod),
2 (flat membrane)**, i.e. MSD = 8Dτ/3, 3Dτ and 4Dτ. A membrane tilted by θ gives 1 + cos²θ, between the
rod and flat values; not offered separately. All four factors verified by simulation to three decimals.

**This is a first-order error on D, not a refinement.** Round-tripped against free diffusion on a
sphere (R = 800 nm) and a rod (r = 400, L = 2000 nm), true D = 200 nm²/ms:

| case | old | fixed |
|---|---|---|
| 3D coords, surface motion | 0.697 | 1.045 |
| 2D coords, surface on a sphere | 0.675 | 1.012 |
| 2D coords, surface on a rod | 0.759 | 1.012 |

Exposed as **Molecule moves** (volume / surface) plus, when surface motion is recorded in 2D, a
**Cell shape** selector. Defaults to volume, so existing results are unchanged. The estimates line
states the relation in force whenever d_eff ≠ ndim. Applied to the ensemble fit, the per-track fits,
the per-cluster models, the bootstrap and the per-segment (HMM state) MSD.

#### Free-diffusion-on-a-geometry MSD reference (added 2026-08-03)
`sim_free_geom_msd` + `msd_geom_null`. On a **closed** cell the MSD of a **free** molecule bends and
eventually saturates on its own, so a bend is not by itself evidence of trapping. This simulates free
diffusion on the chosen shape at the fitted D — same lag grid, same observed coordinates — and overlays
the curve it produces. Measured curve *on* the green line ⇒ the cell explains the bend; *below* it ⇒
confinement beyond the cell.

**This is window-dependent, which is the whole point.** Over a short window the molecule covers far
less than the cell and the two are indistinguishable; the gap opens as the upper fit limit is raised.
Free diffusion on a rod (r = 400, L = 2000 nm, D = 200 nm²/ms, 2D projection), fitted over [4, upper]:

| upper (cycles) | D recovered |
|---|---|
| 50 | 100% |
| 250 | 97% |
| 500 | 97% |
| 1000 | 94% |
| 2000 | 90% |

so a 2000-cycle window costs ~10% of D to geometry alone. The estimates line reports this figure for
the window in use and warns below 95%. It is a Monte-Carlo estimate with ~2% run-to-run scatter
(measured over 10 seeds at M = 4000), so it is rounded to 5% and labelled as approximate rather than
quoted to the percent.

**The shape enters the MSD in two separate places** (verified against simulation to <1%):
1. **The divisor**, at every lag, independent of size — the d_eff table above.
2. **The ceiling**, only once the molecule has explored the shape. The MSD saturates at the mean
   squared separation of two random points on it: **sphere 2R² (3D), 4R²/3 (projected); rod
   L²/6 + 2r² (3D), L²/6 + r² (projected)**. Every term is a **length squared** — the volume
   (4/3·πR³) never appears. Size sets where the curve bends over, not how steeply it rises.

**The ceiling is what fakes biology, and it fakes it in every non-linear model.** Simulating free
diffusion on a 500 nm sphere (2D, surface) over a 4–2000 cycle window and fitting each model returns:

| model | what pure geometry gives |
|---|---|
| linear | ~95% of the true D |
| confined | **L_c ≈ 1000 nm** — a converged, entirely spurious corral |
| anomalous | **α ≈ 0.86** — sub-diffusive with nothing anomalous present |
| directed | v ≈ 0 ✅ — geometry does not fake drift |

So the overlay does not only serve the free/linear case: with it enabled, the estimates line fits
**whichever model is selected** to the simulated free curve and reports what the shape alone would
return, next to the measured value, warning when the shape already accounts for most of it. Note the
d_eff correction itself applies to all four models regardless (it flows through `coef_override` into
`.msd_model_fit`).

**Design notes.** The geometry (shape + dimensions) lives in one shared **Sample geometry** block at the
top of the Analysis tab — `trk_geom_type` plus the dimension inputs — read by the MSD projection factor,
the MSD overlay and the displacement Step-2 null, so the sample cannot be described two different ways.
Options are None / sphere / rod / ellipsoid / flat membrane / **filament**, and *Molecule moves* carries
the matching 1/2/3-DOF choice (a filament covers in-vitro motor assays). `None` and `filament` have no
closed surface to simulate on and map to the unbounded case via `.geom_for_null`, so Step 2 correctly
reduces to the Step-1 band rather than inventing a shape. Earlier drafts had the geometry read from the
**Displacement vs Diffusion** panel, so the cell is defined once and both tests use it; surface/volume follows
the MSD panel's *Molecule moves*. D fed in is the **d_eff-corrected** D this panel now reports, i.e. a
true motion D, so no projection re-calibration is applied — unlike `sim_free_geom_rnet`, which receives
a raw per-axis D and calibrates internally. The geometry operators were factored out into a shared
**`.geom_ops()`** used by both nulls, so they cannot drift apart; the refactor was verified
bit-identical across all five shapes × both modes. Cost scales with the largest lag, not the number of
grid points: ~0.1 s at 50 cycles, 1.6 s at 2000, 4.2 s at the 5000 cap.

**Curvature vs the DOF term.** With d_eff handled, curvature costs only a few per cent of D at the
default short window (~3% at R = 400 nm, ~15% only by R = 100 nm) — the DOF/projection factor above is
much the larger error. Geometry matters mainly when the fit window is long, which is exactly when this
overlay should be enabled.

### I-B.6a Displacement vs diffusion (net displacement test)
`trk_disp_bundle`. Per track, net start→end displacement `R_net = |r(T) − r(0)|` against the free
expectation. Depends only on the two endpoints and the true elapsed time T, so it is **inherently robust
to irregular Δt** — no cycle binning needed. Endpoint noise `2·d·σ_loc²` is subtracted.

**Motion vs observation degrees of freedom (added 2026-07-30).** The band is
`R_net² = 2·D_mot·T·χ²(d_mot)`, where **d_mot is the degrees of freedom of the MOTION, not the number of
recorded coordinates** (`d_obs`). They differ for a molecule confined to a *surface* (membrane, cell
envelope) but tracked in 3D: 2 DOF observed on 3 axes. Two consequences, both now handled:
1. A per-axis estimator (CVE, or MSD with `ndim = 3`) spreads 2-DOF motion over 3 axes and reports
   `D_meas = D_true · d_mot/d_obs` — exactly **2/3** for a surface in 3D (simulated: 0.675). The panel
   rescales it: `D_mot = D_meas · d_obs/d_mot`. A **manual** D is taken as the true motion D and passed
   through unchanged.
2. The χ² degrees of freedom follow the motion, so the band uses `χ²(d_mot)`.

Getting (1) right but not (2) leaves the *mean* correct while the *band shape* is wrong — the two errors
cancel in `⟨R_net²⟩` but not in the quantiles. Measured on simulated free-on-surface tracks
(D = 200 nm²/ms, T = 100 ms, all genuinely free):

| band | sub-diffusive | free | super-diffusive |
|---|---|---|---|
| old (`d_obs` only) | **11.4%** | 81.4% | 7.2% |
| corrected (`d_mot = 2`) | 5.2% | 90.0% | 4.8% |

so surface-bound molecules were being called sub-diffusive at ~2.3× the nominal rate. When
`d_mot == d_obs` (all 2D data, and true volume diffusion in 3D) both terms cancel and the expression
reduces exactly to the previous one — verified unchanged on simulated 3D volume diffusion (5.0/90.1/4.9).
The choice is exposed as **Motion type: volume (3 DOF) / surface (2 DOF)**, shown only for 3D data,
defaulting to volume so existing results are unaffected.

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
**Irregular-Δt handling (2026-07):** both quantities use only **integrated-ish steps** — those spanning
≤ `input$trk_dir_max_dark` dark cycles (default 2; `n_dark = round(Δt/t_cycle) − 1`); longer-gap steps
are dropped (unreliable direction, possible cross-molecule merge). The autocorrelation lag is
**elapsed time**, binned by integer cumulative-cycle index (as in the MSD) rather than step index, so
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
  fold-over-the-track's-own-baseline (validated ≈5% false-positive at fold≈3), floored at the
  localization-noise velocity.

### I-B.9 Spatially-resolved D-map — [Part II §5 note; single-lag = short-lag CVE limit]
`dmap_bundle`. Single-lag estimator per step: **D_raw = r²/(4·Δt)** (2D), midpoint-binned on an
n×n grid; per-cell **D = median(D_raw)** (robust to outlier jumps), cells with < `min_pts` steps left NA.
Options: **integrated-only** (drop dark-gap steps whose midpoint is spatially meaningless) and **noise
subtraction** `D_raw − σ_loc²/Δt` (removes the spurious immobile-region diffusion floor). Far more robust
than a per-cell MSD fit in small cells (e.g. *E. coli*) where no τ-bin has enough steps. *(Citation
corrected: this is the CVE short-lag limit, Vestergaard 2014 / Vogler 2025 — not Ostersehlt 2022.)*

### I-B.10 Fast-diffusion ceiling (D_max)
`D_max = (L²/16)·ψ / (N_PL + ψ·t_hw)` (Vogler 2025) — the trackable-D limit before the localization loop
loses lock; surfaced as a green/amber/red flag. Theoretical form, so it may under-warn near the
demonstrated ~2.5–3 µm²/s practical ceiling (Part II §3a item 8).

## I-C. Simulation

Two independent simulators share one precision model: **σ_loc = L_TCP / (2·√N_photons)** (so σ_loc is
*heteroscedastic* — varies per localization with the drawn photon count, as in reality; equivalently the
effective TCP diameter L_eff = 2·σ·√N).

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
  budget; τ and floor calibrated from in-house real runs (`reference_minflux_detection_rate`). A DBSCAN
  core-point test (≥minPts within eps=60 nm) splits recorded traces into **structure vs background** so the
  measured rate curves can be overlaid on the model.
- **Output:** MINFLUX-format `data.frame(tid, tim[s], x, y[,z], efo, eco)` + a molecule/fluorophore
  ground-truth sidecar — loads straight into the SMLM tab so the GT inverse can be scored against truth.

### I-C.2 Diffusion / tracking ground-truth simulator (`sim_server.R`, Diffusion tab)
Simulates a *known* diffusion process on a *known* geometry, samples it with **MINFLUX's own pathologies**,
and emits MINFLUX-format tracks so they flow through the *exact same* Tracking analysis — closing the
validation loop (see [simulation_module.md](simulation_module.md) for the full design + build phases).

**Geometry operators** (unified `project` / `reflect` / `normal`):
- **Rod** = spherocylinder (cylinder body radius r, length L_cyl + two hemispherical caps): `project` =
  nearest surface point, `reflect` = fold back inside, `normal` = outward unit normal.
- **Sphere / ellipsoid** (a=b=c for sphere): radial project/reflect in scaled coordinates.
- **Membrane** = flat sheet with thickness (2D lateral + z wobble within thickness); **Free** = unbounded 3D.

**Continuous path** (`gen_path`, per cycle): per-axis step SD **sdi = √(2·D·t_c)**.
- *Surface* diffusion: tangent-plane step `g − (g·n̂)n̂` then `project` back onto the surface.
- *Volume* diffusion: 3D step then `reflect` at the wall. *Free/membrane*: unbounded / 2D-lateral.

**Motion modifiers** (each a ground-truth label):
- **Free** Brownian.
- **Confined = OU tether** to a per-track home: `p ← home + (p − home)·φ`, **φ = e^{−t_c/τ_c}**,
  **τ_c = L_c²/D**. Validated: stationary SD per DOF = L_c exactly; total MSD plateau = 2d·L_c².
- **Directed:** add drift `v·t_c` along the rod long axis (or random direction otherwise).
- **Two-state:** mobile↔immobile CTMC at per-cycle switch prob `kswitch`; immobile D = 0.01·D. HMM ground truth.

**MINFLUX sampling overlay** (the load-bearing part):
1. **Cycle clock** t_c; a localization is attempted each cycle.
2. **Dark cycles → irregular Δt:** each cycle recorded with prob `1 − p_dark`, else skipped, giving
   variable-cycle gaps.
3. **Photobleaching → short tracks:** track length in cycles is **Geometric**, with the mean scaled so the
   *recorded median* (after dark-cycling) lands on the target `mean_len` (median ≈ ln2·mean, only
   (1−p_dark) recorded ⇒ mean_cycles = median/(ln2·(1−p_dark))).
4. **Heteroscedastic σ_loc:** per recorded loc draw **N ~ Poisson(photons)**, σ = L_TCP/(2√N), add Gaussian
   noise per axis; record `eco = N`, `efo = N/(t_c in s)`.
5. **Survivorship (optional):** if the true displacement across a dark gap exceeds a re-acquisition radius
   `R_reacq`, the track terminates there — reproduces the dark-cycle survivorship bias (long jumps culled →
   *manufactured* apparent confinement).
6. **Projection:** output 2D (drop z) or 3D — a surface protein on a curved cell shows the 2D-projection
   distortion we want to *study*.

**Built-in CVE round-trip** (`sim_trk_cve`) scores the recovery immediately:
**D = (⟨|Δr|²⟩ + 2⟨Δr_i·Δr_{i+1}⟩) / (2·d·t_c)** over consecutive single-cycle steps, using the **total**
displacement with **d = min(motion-DOF, observed axes)** — surface diffusion has 2 DOF but spreads over 3
Cartesian axes (a per-axis CVE would under-estimate D by ⅓); volume has 3; 2D projection caps d at 2.
Validated D recovery 97–101% across geometry×mode combos (flat geometries; the curved cylinder's 2D
projection recovers ~75%, the geometric-bias factor the projection correction can be lifted from).

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
| **Cycle time** (one localization) | ~0.4 ms–a few ms; set by photon budget, **not** a fixed clock | No single nominal Δt | Balzarotti 2017; Wolff 2023; Deguchi 2023 |
| **Irregular Δt** | consecutive locs separated by integer multiples of a base cycle (skipped/"dark" cycles) | **Estimators assuming constant Δt are biased** | Scheiderer/Marin/Ries 2025; Marin & Ries 2025 |
| **Gaps / dark cycles** | photobleaching budget + µs–ms blinking abort/skip cycles — *structural*, not noise | Must be modelled, not smoothed away | Marin & Ries 2025 (SimuFLUX) |
| **Motion blur (R factor)** | effectively ≈ 0 (short cycle + active re-centering) | Blur-correction terms ≈ 0 (see caveat §5) | Scheiderer/Marin/Ries 2025 (inferred) |
| **Track length N** | photobleaching-limited; ~100× more locs/track than camera-SPT, but still short (tens–low hundreds) | Short-track estimator bias dominates | Balzarotti 2017; Deguchi 2023 |
| **Trackable D ceiling** | ≈ 3 µm²/s (standard settings) before the loop loses lock | Fast diffusion is systematically under-estimated | Marin & Ries 2025; Vogler 2025 |

**The two properties that drive almost every rating below are (a) irregular Δt and (b) short tracks.**

---

## 2. Methods implemented — sources & appropriateness

### 2.1 MSD with optimal point-count fitting — 🟡 AMBER
- **Implemented as:** unweighted OLSF `MSD = a + b·τ` with Michalet–Berglund optimal number of fit
  points (`p_a`/`p_b`), `R = 0` (MINFLUX), plus a **cycle-based MSD** that rounds each lag to the
  nearest integer number of acquisition cycles.
- **Source:** Michalet 2010; Michalet & Berglund 2012 (optimal least-squares fitting / optimal `p*`).
  Cycle-based binning and the trackable-D ceiling: **Vogler et al. 2025**.
- **MINFLUX caveat:** the optimal-`p*` formalism is *built for* short, noisy tracks (a match), and
  `R≈0` is a clean limiting case the theory already covers — but the framework assumes a single
  nominal Δt for lag binning. Binning by *step index* on irregular-Δt data corrupts the MSD curve.
- **Our handling:** the cycle-based MSD (Vogler 2025) bins by **integer cycle count** = true elapsed
  time, which is exactly the recommended fix. ✅ This lifts the practical rating toward GREEN *for the
  cycle-based path*; the plain step-index MSD retains the AMBER caveat.

### 2.2 CVE — covariance-based D & σ_loc estimator — 🔴 RED (textbook) → 🟡 mitigated
- **Implemented as:** `σ_loc² = −⟨Δxᵢ·Δxᵢ₊₁⟩` (Vestergaard's R=0 / no-motion-blur case) over
  **consecutive integrated (single-cycle) step pairs** only; also returns bias-corrected `D_cve`.
- **Source:** Vestergaard, Blainey & Flyvbjerg 2014.
- **MINFLUX caveat:** the CVE algebra is derived for a **constant** Δt. Applied naively across
  variable-Δt MINFLUX steps it biases D — hence the raw RED rating.
- **Our handling:** we **restrict to consecutive single-cycle pairs** that share a midpoint within one
  track (`n_cyc = round(Δt/t_cycle) == 1` on both steps), i.e. we pre-segment into constant-Δt
  sub-runs — precisely the mitigation the literature prescribes. ✅ Statistically valid on MINFLUX at
  the cost of discarding dark-gap-adjacent steps.
- **Photon-quality weighting: tried, then REMOVED** (added 2026-07-19, withdrawn; this entry corrected
  2026-07-30 — it previously described the weighting as live). The idea was that precision scales with
  photons, so fitting `σ_loc ∝ photons^p` and inverse-variance-weighting the pooled σ_loc and D would
  exploit the quality channel instead of a hard EFO cut. It was withdrawn on discovering that in
  **tracking** the photon channel is a **mobility** proxy: a fast molecule leaves the beam pattern sooner
  and collects fewer photons, so photon-weighting up-weights slow molecules and biases D **downward**.
  The diagnostics that revealed this (`noise_frac`, `phot_mob_rho`) are retained and reported — and the
  same coupling means **any ECO/EFO quality filter is also a mobility filter**, which the panel now warns
  about. σ_loc co-estimation is inherently less critical for MINFLUX (σ is already ~nm), so CVE serves
  mainly as an MSD-free precision/plausibility check.
- **Lateral vs axial σ (2026-07-30):** σ is pooled over x, y and z, but the axial precision is usually
  worse; the purely lateral consumers (D-map, HMM) then over-subtract their noise floor. `sigma_xy` /
  `sigma_z` / their ratio are now reported separately, with a warning above 1.3× (see §I-B.1).

### 2.3 JDD — jump-distance / chi-mixture — 🟢 GREEN (with our handling)
- **Implemented as:** jump-distance distribution with an N-component scaled-chi mixture (**Rayleigh in
  2D, Maxwell in 3D** — see I-B.4); σ_loc bias correction, √Δt-normalization (default on), and an
  **integrated-steps-only** mode.
- **Source:** Schütz, Schindler & Schmidt 1997; benchmarking in Weimann et al. 2013.
- **MINFLUX caveat:** small σ_loc reduces error-broadening of the Rayleigh peaks (an advantage), but
  jumps must be pooled at a **fixed real Δt**; mixing steps of different cycle counts pollutes the
  distribution.
- **Our handling (2026-07):** two rigorous options — (1) **√Δt-normalise** each jump (default on) so
  mixed cycle counts pool correctly for free diffusion; or (2) **Integrated steps only** (`round(Δt/t_c)
  ==1`), which keeps a single true Δt so *no* pooling assumption is needed at all — valid regardless of
  the motion model. The "short track → few jumps" concern is weak for MINFLUX (100× more locs/track than
  camera SPT). Residual: the √Δt mode still assumes free diffusion; prefer integrated-only in doubt, and
  keep the component count low (each component still assumes free Brownian).
- **Dimensionality fix (2026-07-30):** the fitted PDF is now the **d-correct** distance law rather than
  Rayleigh-regardless. On 3D data the old form biased single-population D by +32% and badly distorted
  multi-population splits (a true 100/900 mixture read as 281/840). 2D results are unchanged. See the
  boxed note in I-B.4 for the validation numbers.

### 2.4 HMM diffusive-state segmentation — 🟡 AMBER
- **Implemented as:** per-track Gaussian-emission HMM over displacement steps, N states, with
  Δt-normalization option.
- **Source:** Das, Cairo & Coombs 2009; variational-Bayes short-track pooling (vbSPT) Persson et al. 2013.
- **MINFLUX caveat:** vbSPT's strength (pooling power across many *short* tracks) and MINFLUX's fine
  precision (better state separability) are a good match — but the discrete-time transition matrix
  assumes a **uniform Δt per step**. Irregular cycles need either resampling to a common grid or a
  continuous-time-rate transition model.
- **Our handling (updated 2026-07):** the `norm_dt` option addresses emission scaling (state σ/D robust),
  and the reported **transition rates are now continuous-time**: a generator Q estimated from the decoded
  path using real elapsed time (rates in s⁻¹, mean dwell τ = −1/qₐₐ), invariant to dark cycles. The
  per-step matrix A is kept for reference. Residual: the *decode* still uses the Δt-agnostic A as its
  transition prior (emission-dominated, so state assignments are robust); a fully continuous-time decode
  (A^Δt propagation) is the remaining refinement. Rating lifts toward 🟢 for the rate outputs.

### 2.5 Confinement analysis — 🟢 GREEN (with the N-matched null)
- **Implemented as:** a **per-track spatial shape ratio Rg/Rmax** classified against an **N-matched
  free-diffusion null**, with σ_loc² subtraction from both Rg and Rmax. *(A static spatial descriptor, not
  Simson's temporal sliding-window index; the app also exposes a τ-based MSD confinement index separately.)*
- **Source:** confinement-detection SPT lineage — Simson, Sheets & Jacobson 1995; Saxton, Kusumi.
- **MINFLUX fit:** being purely spatial it is **Δt-independent** (irregular sampling does not bias it), and
  nm precision resolves domains far below the camera-SMLM ~20–30 nm floor — a genuine strength.
- **Fixed a real defect (2026-07):** the old fixed thresholds ("confined <0.35 / free >0.45") were both
  mis-calibrated and **inverted** — diffusion simulation shows confined motion gives a *higher* Rg/Rmax
  than free (≈0.60 vs 0.53). Replaced by the N-matched null (confined = above free 95th pct at that N;
  directed/escape = below the 5th), which also removes the track-length confound. See I-B.6 for the
  validation.
- **Residual caveats (both conservative):** the Rmax extreme-value selection bias (no closed form); and a
  genuinely confined but very slow molecule can be missed if it barely explores. Corroborate with the
  τ-based MSD confinement index / plateau.

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
- **Our handling / recommendation:** report **pooled / ensemble** α with confidence intervals and
  short-track exclusion; treat **single-track α as qualitative**. Consider surfacing Kepten's
  track-length caveat in the UI info box.

### 2.8 Change-point detection of motion states — 🟡 AMBER
- **Implemented as:** PELT segmentation of trajectories (via `changepoint`), min-segment / penalty controls.
- **Source:** Killick, Fearnhead & Eckley 2012 (PELT).
- **MINFLUX caveat:** PELT is sampling-agnostic *if* the per-segment cost is written in **true elapsed
  time**; off-the-shelf costs assume evenly spaced points. Short tracks limit resolvable segments.
- **Our handling:** verify the segment cost is Δt-aware (see §5).

### Summary

| Method | Source | Textbook | With our adaptation |
|---|---|---|---|
| MSD optimal-fit | Michalet 2010/2012; Vogler 2025 | 🟡 | 🟢 cycle-based path / 🟡 step-index |
| CVE | Vestergaard 2014 | 🔴 | 🟡 (integrated-pair restriction + EFO weighting) |
| JDD | Schütz 1997; Weimann 2013 | 🟡 | 🟢 (√Δt norm + integrated-only mode) |
| HMM / vbSPT | Das 2009; Persson 2013 | 🟡 | 🟢 rates (continuous-time generator) |
| Confinement | Simson 1995 | 🟢 | 🟢 (N-matched free-diffusion null) |
| Turn-angle | Burov 2013 | 🟢 | 🟢 (elapsed-time autocorrelation) |
| Anomalous α | Metzler 2014; Kepten 2015 | 🔴 | 🔴 (pool + exclude short) |
| Change-point / PELT | Killick 2012 | 🟡 | 🟢 (Δt-aware cost + σ_loc floor) |
| D-map (single-lag) | see §5 citation note | 🟢 | 🟢 |

---

## 3. MINFLUX-specific parameter tools

- **D_max ceiling / trackable-D limit** — `D_max = (L²/16)·ψ / (N_PL + ψ·t_hw)`, cited in-code to
  **Vogler et al. 2025** (*Commun. Biol.*). Consistent with the independent SimuFLUX ceiling of
  ~3 µm²/s (Marin & Ries 2025). 🟢 Correctly attributed.
- **σ_loc noise floor / dark-cycle handling / continuity gate** — MINFLUXR-original engineering
  built on the data properties in §1; not a port of a specific paper.

---

## 3a. MSD on MINFLUX — computation, confinement & what is non-standard

The MSD tool follows **Bela Vogler's `the_road_to_minflux_msd` workflow** (Vogler et al. 2025), the
only published MSD treatment specific to MINFLUX. The research landscape here is genuinely sparse;
this section records what is *standard* (safe) versus *non-standard/unsolved* (deliberately not done).

### How MSD is (and should be) computed on MINFLUX — checklist
1. **Bin by integer cycle-count lag, not wall-clock time.** MINFLUX cycle-skipping makes Δt
   non-uniform; group displacement pairs by cycle-index η (Vogler 2025). Our cycle-index binning does this.
2. **Segment traces at abnormal cycle gaps before MSD** (Vogler split when η > ⟨η⟩+SD). Our gap-split
   + continuity gate is the equivalent; prefer integrated (single-cycle) steps with tolerance.
3. **Never fit a fixed/arbitrary lag range.** Michalet & Berglund 2012: a fixed choice ("first 4
   points", "10 %") has much larger uncertainty than the computed optimum p*. We use OLSF + Michalet
   optimal-point selection for the **linear** model. *(NB: this optimal-p theory is derived for the
   linear/free-diffusion case ONLY — it is NOT valid for a plateauing confined curve; our
   confined/anomalous fits deliberately use the whole window instead.)*
4. **Both ends of the MSD are unreliable, for different reasons:** short lags — variance inflates when
   the reduced error x = σ²/(DΔt) > 1 (Michalet & Berglund 2012), plus probe/linker wobble; long lags —
   fewer, statistically *correlated* pairs (Qian, Sheetz & Elson 1991). Hence the fit window [lower,
   upper] cutoff. Our default 5–50 cycles matches Vogler's `fit_limits=(5,50)`.
5. **The exact 5-cycle lower bound is probe/dataset-specific, NOT a universal constant.** Vogler's
   two-regime "kink" (5–9 lags) reflects their QD/PEG-linker probe's non-diffusive wobble; port the
   *diagnostic logic* (inspect the curve for where the low-lag regime ends) not the number.
6. **Model σ_loc and motion blur in the linear fit:** ρ(iΔt)=2d(σ²−2RDΔt)+2dD·iΔt, with **R≈0** for
   MINFLUX (no integration). Small σ_loc → small x → favours *fewer* MSD points, not more.
7. **Minimum track length is a JUDGEMENT CALL, not a literature constant.** There is **no
   MINFLUX-calibrated minimum-localizations rule** (verified). Our default of **50** is *Bela's
   practical choice* (`filter_by_length(min_length=50)`), not a derived threshold. Kepten et al. 2015's
   L≤300 guidance is for camera SPT with far larger σ_loc and is not directly transferable.
8. **D_max: theoretical ≠ demonstrated.** Vogler's Eq. 3 gives D_max ≈ 10–11 µm²/s theoretically, but
   they only *demonstrated* tracking to D ≈ 2.5 µm²/s (lipids). Our D-ceiling warning uses the
   theoretical form, so it may under-warn near the ~2.5–3 µm²/s practical ceiling — interpret the
   green/amber/red flag conservatively.

### Confined diffusion (the E. coli case)
Vogler's toolbox is **free-diffusion only**; his README states anomalous/confined analysis is "under
development / implement on your own." So confined MSD on MINFLUX is genuinely **non-standard**.
- **Standard confined-MSD model** (when free diffusion is ruled out): the plateauing corral form,
  MSD(t)=MSD∞·[1−A₁·exp(−4A₂D t/MSD∞)] introduced by **Kusumi, Sako & Yamamoto 1993**; the "generic"
  square-domain version MSD(t)=(L²/3)·[1−(8/π²)·exp(−t/τ)], τ=L²/(π²D). Companion corral/escape
  framework: **Saxton 1993/1994/1995**. (Destainville & Salomé 2006 is a motion-blur *correction* to
  this form, **not** the confinement derivation — do not mis-cite it as such.)
- **Detect confinement:** a resolvable MSD plateau vs. a straight line, or α<1 in MSD=4D·t^α (Saxton
  1994; Metzler & Klafter 2000; Metzler et al. 2014). α<1 is a *symptom with several possible causes*,
  not proof of confinement.
- **E. coli / rod geometry has NO closed-form MSD.** English et al. 2011 and Bakshi et al. 2011 handle
  it empirically (long- vs short-axis decomposition; Monte-Carlo of a spherocylinder), not with a
  portable formula. vbSPT (Persson et al. 2013) sidesteps geometry with an HMM over diffusive states —
  a different question than an MSD curve, not a drop-in replacement.

### Deliberately NOT implemented (non-standard / unsolved)
- Confined-MSD **combined with MINFLUX cycle-based irregular Δt** — does not exist in the literature;
  Vogler built the irregular-Δt machinery but applied it only to free diffusion and explicitly *ruled
  out* confinement rather than modelling it.
- Michalet optimal-point selection on a **plateauing** curve (an unsupported extrapolation).
- Rod-cell **geometry correction** for E. coli (no closed form).
- Non-Gaussian / viscoelastic / superstatistical bacterial-cytoplasm models (e.g. Sadoon & Wang 2018)
  — active research, not consensus.

In the app these appear as the **Confined** / **Anomalous** fit models, labelled *beyond-standard for
MINFLUX* and exploratory, with a free-diffusion honesty note on the Linear result. The user decides.

---

## 6. SMLM (imaging) analysis tools — sources & appropriateness

The SMLM tab applies established single-molecule-localization-microscopy (SMLM) and spatial
point-pattern methods to MINFLUX localizations. MINFLUX-specific method development for *imaging*
statistics barely exists, so the honest position is: **the estimators are standard; the MINFLUX
caveat is the labelling/blinking geometry**, not the estimator maths.

### 6.0 The pervasive caveat — labelling & blinking inflate every spatial statistic — 🟡 AMBER
A single molecule is reported as several localizations: fluorophores blink and re-appear, and a
multivalent antibody label scatters several dyes ~10–40 nm around one epitope. This **over-counts and
spatially spreads** each molecule, which inflates short-range nearest-neighbour peaks, Ripley/PCF
signal, and apparent cluster counts as a *label artefact*, not biology (Annibale et al. 2011;
Sengupta et al. 2011; Veatch et al. 2012). **Our handling:** the "median (de-blinked) clustering"
tier collapses localizations → trace → primary cluster → median cluster before spatial statistics,
which removes most multi-blink over-counting; results on *raw* localizations should be read with this
bias in mind. This is a mitigation, not a cure (it cannot recover true molecule count below the label
footprint).

### 6.1 Nearest-neighbour, DBSCAN/HDBSCAN/OPTICS clustering — 🟢 GREEN (method)
- **Source:** Clark & Evans 1954 (NN statistic); Ester, Kriegel, Sander & Xu 1996 (DBSCAN);
  Ankerst et al. 1999 (OPTICS); Campello, Moulavi & Sander 2013 (HDBSCAN).
- **Fit:** density-based clustering is well-adapted to SMLM and widely used; DBSCAN in particular is a
  community standard. MINFLUX's fine precision *helps*. The only caveat is §6.0 (feed de-blinked
  medians, tune eps to the label footprint).

### 6.1a Traces per cluster — 🟡 AMBER
- **Source:** cluster-occupancy counting in the SMLM molecular-counting lineage — Annibale et al. 2011
  (blinking artefacts in PALM); Sengupta et al. 2011 (PC-PALM).
- **Fit:** counting *traces* rather than localizations is the right unit for MINFLUX and removes the
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
- **Edge-bias caveat — now mitigated by the same-window CSR envelope.** The lack of an analytic edge
  correction previously meant L−r was read only for shape/peak, not magnitude. The optional Monte-Carlo
  envelope (§I-A.3/§I-A.4) simulates CSR in the **identical frozen window**, so the edge bias is present
  in the band and cancels; significance is read as the observed curve leaving the 2.5–97.5% band rather
  than departing from the (biased) flat zero. This supersedes an explicit translation/isotropic
  correction for interpretation (a `# TODO` for the analytic correction remains).
- **✅ No inter-channel registration error — a genuine MINFLUX advantage.** Unlike camera-based two-colour
  SMLM, the DCR (ratiometric) two-colour mode excites both fluorophores with the **same laser** and
  localizes them through the **same beam path in the same iterative sequence**, so both channels are
  reported in one common coordinate frame. There is no chromatic offset and no affine channel
  registration to calibrate — the usual dominant error term of two-colour co-localization simply does not
  exist here, and L₁₂/g₁₂ are trustworthy right down to the localization precision. *(Exception: the
  multiplex-by-files and multiplex-by-time-window modes overlay **separate acquisitions or separate time
  windows**, so sample/stage drift between rounds can offset them — that is drift, not chromatic
  registration, and it applies only to those modes, not to the DCR split.)*
- **⚠ The real two-colour error is DCR mis-assignment (crosstalk), not offset.** The channel split is a
  threshold on a **continuous DCR distribution whose two populations overlap**, so a fraction of
  molecules is assigned to the wrong channel. This does not shift distances — it moves *points between
  the two sets*, which **manufactures short-range cross-channel pairs**: a molecule wrongly labelled
  channel B sits, by construction, right next to its true channel-A neighbours. Crosstalk therefore
  biases L₁₂ and g₁₂ **upward at short r — toward apparent co-localization**, exactly where a real
  interaction would show. Asymmetric labelling density between the two channels also skews the
  estimators. Finally, co-localization at nm scale is **proximity, not interaction**.
- **Now quantified (`auto_detect_dcr_split`).** The GMM auto-split returns, alongside the threshold s:
  - **d′ = |μ₂−μ₁| / √((σ₁²+σ₂²)/2)** — standardised separation of the two modes. Verdicts used in-UI:
    **d′ ≥ 3 clean, ≥ 2 workable, < 2 poor** (channels overlap heavily).
  - **Expected mis-assignment rate** `π₁·P(X₁>s) + π₂·P(X₂≤s)`, and the more directly usable
    **per-group contamination** — what fraction of everything *called* group k actually came from the
    other dye: `contam₁ = π₂P(X₂≤s) / [π₁P(X₁≤s) + π₂P(X₂≤s)]`, symmetrically for group 2.
  - **Ambiguous fraction** — traces whose best GMM posterior is < 0.95, i.e. the ones that would switch
    channel under a small threshold change. These are the natural candidates for a strict gate.

  **Read any cross-channel result against these numbers:** a co-localized fraction of the same order as
  the contamination is not evidence of interaction.

  *Validation (synthetic, known truth):* the predicted mis-assignment rate matches the empirical one
  **exactly when the split is good** (d′ = 3.9 → predicted 2.56% vs actual 2.53%; d′ = 8 → 0% vs 0%),
  but drifts ~3 pp when it is not (d′ ≈ 1 → predicted 27.7% vs actual 30.7%). The crosstalk estimate is
  itself least reliable exactly where it matters most, so **"poor" should be read as "do not use", not
  as "use with this correction"**. Contamination is correctly asymmetric under unequal populations
  (80/20 mix → minority group 4.1% contaminated vs 1.3% for the majority).
- **Scope of d′: it grades a split, it does not detect one.** The GMM is a sound estimator of the split
  and is the right tool here — the point is only that requesting G = 2 always yields two components,
  because a single-dye DCR distribution is **skewed** and a two-Gaussian mixture fits skewed unimodal
  data better than one Gaussian. On a **confirmed single-colour file** (collaborator, no DCR pooling applied)
  the routine still returned a *"workable"* split at d′ = 2.65. **A BIC G=1 vs G=2 comparison does not
  help and was deliberately not added as a gate:** it also favoured two components on that file
  (ΔBIC = +16.9) and on a synthetic skewed-unimodal gamma control (+21.1), while *missing* a genuine
  two-colour pair at d′ = 1 (ΔBIC = −9.6) — it detects skewness, not bimodality. So: use d′ and the
  crosstalk numbers to judge **how good** a two-colour split is, once the experiment already establishes
  that there are two dyes.
- **✅ Per-molecule split confidence (`dcr_split_confidence`).** The aggregate crosstalk numbers cannot
  say *which* molecules are uncertain, and those are exactly the ones that manufacture false short-range
  cross-channel pairs. For each trace the fitted mixture gives the posterior of its assigned component,
  `post_k(v) = π_k·N(v|μ_k,σ_k) / Σ_j π_j·N(v|μ_j,σ_j)` — 1.0 unambiguous, 0.5 sitting exactly on the
  split — plus a log-odds form and an `on_edge` flag (posterior < 0.95). Shown in the 2-colour panel as
  posterior-vs-DCR with the edge band highlighted, and exportable per trace.
  *Validated against known truth (synthetic, d′ = 4.1):* the posterior is **calibrated** — traces at
  posterior 0.5–0.7 were 55% correct, 0.7–0.9 → 79%, 0.9–0.99 → 95.5%, ≥0.99 → 99.9%, each within a few
  points of nominal (slightly over-confident in the middle bins, as expected when the mixture parameters
  are estimated rather than known). Practically: **dropping the edge traces cut the mis-assignment rate
  from 2.43% to 0.40% while keeping 90.7% of the data** — a ~6× reduction in crosstalk for ~9% data loss,
  which is what makes a strict gate worth offering.
  *Implementation note:* the confidence uses **unrounded** mixture parameters (`gmm$params`); the
  top-level `mean/sd/pi` fields are rounded for display, and using those made the per-trace edge count
  disagree with `ambiguous_frac` and cost ~0.03% of assignments.
- **✅ Pool DCR at upload — the most effective mitigation, where it is available.** `dcr_pooled`
  (`helpers_shared.R`) is the photon-weighted DCR across the whole iteration sequence, `Σ(dcr·eco)/Σeco`,
  rather than the final-iteration value. Averaging over the sequence with photon weighting **sharpens the
  two modes substantially**, so the split quality — and therefore the crosstalk floor on every
  cross-channel statistic — improves at source. This beats any downstream statistical correction: the
  split is the error, so fix the split.
  ⚠ **Only .npy and .json uploads can be pooled** — pooling happens in `read_minflux_npy_as_df()` /
  `read_minflux_json_as_df()` via their `pool_dcr` argument. A **.csv is a pyMINFLUX export of the final
  iteration**, already through that pipeline, so there is nothing left to pool and the UI must not nag
  about it. The in-app hint is therefore conditional: it points at the `dcr_pooled` column when one
  exists and is unused, and otherwise explains the format restriction.
  The split-quality numbers are only comparable between runs that used the **same** DCR column.
  *Quantified on a simulation of the mechanism* (noisy per-iteration DCR, photon-weighted pooling vs
  final-iteration value): pooling took **d′ from 1.56 → 3.59** and expected crosstalk from
  **17.9% → 3.6%** — moving the split from "poor" to "clean" and cutting mis-assignment ~5×. That is a
  simulation of the averaging effect rather than a measurement on real two-colour data, but it is
  consistent with the improvement reported in practice.
- **g₁₂ vs L₁₂ on MINFLUX data.** g₁₂ is the more informative of the two here because it is scale-local:
  a short-r blinking/labelling artefact stays confined to the first few bins, whereas in the cumulative
  L₁₂ it contaminates every larger r. The price is noise — each annulus holds few pairs, so short-r bins
  are unstable and an apparent first-bin peak is often binning, not biology.

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
- **Our handling:** computed on **per-TID median positions**, one point per trace. This matters — on raw
  localizations each molecule is weighted by how long it was tracked, so a single bright fluorophore
  contributing thousands of localizations pulls Rg toward itself and shrinks it. De-blinking first makes
  Rg a property of the *structure*, not of the acquisition.
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
  parameter, not a property of the microscope. Only compare r values computed at the **same bandwidth on
  similarly dense data**; prefer the sweep to any single number.
- **The randomization null was anti-conservative — improved but NOT solved, 2026-07-30** (§I-A.5). It
  scattered Ch2 *uniformly*, testing against complete spatial randomness rather than independent
  placement. Measured on simulated independent channels confined to one cell inside a larger window, that
  null falsely called co-localization **93%** of the time (calibrated = 5%). The default is now a
  **toroidal random shift** of Ch2, which preserves Ch2's clustering and count — but it still rejected
  **53%** in the same scenario, because both nulls randomize over the **bounding box, not the cell mask**.
  Mask-aware randomization is the outstanding fix. Until then this p-value is only trustworthy when the
  signal fills the field; otherwise use g₁₂(r) with its envelope. Also a **citation fix**: this is not the
  Costes et al. 2004 block-randomization algorithm and is no longer credited as such. Rating stays 🟡.

### 6.5 Cluster shape fitting (viewer) — 🟠 EXPLORATORY
- **Method (standard):** orthogonal-distance MLE + AICc/BIC model selection (Akaike 1974; Schwarz 1978;
  Hurvich & Tsai 1989).
- **Why exploratory:** the shape dictionary (sphere/ellipsoid/helix/ring/arc) is imposed by this tool and
  the "best" shape is a model-selection outcome, not proof of biological geometry; with few,
  blinking-inflated points several shapes are statistically indistinguishable. **The true shape may not
  be in the dictionary at all — the ranking will still name a winner.**
- **Caveat (residual model):** the information criteria assume iid isotropic-Gaussian residuals, but
  MINFLUX precision varies with photon count, so points are not equally precise (§I-A.7). The likelihood
  is approximate; AICc differences of a few units are not meaningful. RMSD always favours the most
  flexible model — prefer AICc.

### 6.6 Ground-truth position/count estimate & structure size — 🟠 EXPLORATORY
- Labelling-model inverse developed in this tool, validated *in silico* against the simulator. Contour
  via minimum-spanning-tree (Prim 1957); anisotropy coefficients as §6.3. Reported as a range, not a
  single number, below the labelling resolution. See memory `project_gt_estimator`.

### 6.7 Detection-over-time / crowding-vs-timing / NN-vs-size — 🟠 EXPLORATORY
- Bespoke descriptive diagnostics of this tool with **no published MINFLUX equivalent**. They are
  correlational, not causal: detection timing and crowding co-vary with photobleaching, blinking and
  labelling, so a trend is a hypothesis for controls, not a measurement. The detection-rate decay used
  to calibrate the simulator is characterised empirically from in-house acquisitions (memory
  `reference_minflux_detection_rate`), not from a publication.
- **Three specific pitfalls beyond "correlational"** (detail and rationale in §I-A.9):
  1. **The KS gap is an effect size, not a test.** The clustered and background samples are not
     independent (same acquisition, same detection budget), so the KS null does not hold. The code
     correctly reports only the gap — do not convert it to a p-value.
  2. **Crowding-vs-timing has a positive trend built in — now handled** (§I-A.9). A within-file
     label-permutation null gives each structure its own 95% interval and a Spearman trend test, so the
     slope is read against the no-crowding expectation rather than against zero. The band alone still
     reflects offset + trend; the ρ test isolates the trend.
  3. **NN-vs-size is right-censored under DBSCAN** at the link radius `eps2` (default 200 nm), so a
     flattening curve may be the ceiling rather than constant spacing. Re-run under HDBSCAN/OPTICS, which
     impose no such bound, before interpreting a slope.

### SMLM summary

| Tool | Source | Rating |
|---|---|---|
| NN / DBSCAN / HDBSCAN / OPTICS | Clark & Evans 1954; Ester 1996; Ankerst 1999; Campello 2013 | 🟢 (feed de-blinked) |
| Traces per cluster | Annibale 2011; Sengupta 2011 | 🟡 blinking (not stoichiometry) |
| Ripley L / L₁₂ / PCF g₁₂ | Ripley 1977; Besag 1977; Lotwick & Silverman 1982; Dixon 2002; Sengupta 2011; Veatch 2012 | 🟢 / 🟡 labelling + registration |
| Cross-channel NN & co-loc curve | Clark & Evans 1954; Dixon 2002; Lachmanovich 2003; Dunn 2011 | 🟢 (no null of its own) |
| Radius of gyration / anisotropy | Rudnick & Gaspari 1986; Westin 1997/2002 | 🟢 |
| KDE + Pearson co-loc | Silverman 1986; Dunn 2011; Costes 2004 | 🟡 bandwidth |
| Shape fitting | Akaike 1974; Schwarz 1978; Hurvich & Tsai 1989; Pearson 1901 | 🟠 exploratory |
| GT estimate & structure size | this tool (in-silico validated); Prim 1957 | 🟠 exploratory |
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

**Reviews**
22. Saxton MJ, Jacobson K (1997). *Single-Particle Tracking: Applications to Membrane Dynamics.* Annu. Rev. Biophys. Biomol. Struct. 26:373–399. doi:10.1146/annurev.biophys.26.1.373
23. Manzo C, Garcia-Parajo MF (2015). *A review of progress in single particle tracking…* Rep. Prog. Phys. 78:124601. doi:10.1088/0034-4885/78/12/124601
24. Shen H, Tauzin LJ, Baiyasi R, et al. (2017). *Single Particle Tracking: From Theory to Biophysical Applications.* Chem. Rev. 117(11):7331–7376. doi:10.1021/acs.chemrev.6b00815

---

## 4a. Implementation audit (2026-07) — most Δt handling was already in place

A code audit against the ratings above found the implementation is **more MINFLUX-appropriate
than the textbook ratings suggest** — several mitigations were already present:

- **CVE** — restricts to consecutive *integrated single-cycle* step pairs (constant-Δt sub-runs). ✅
- **JDD** — `Normalize steps by √Δt` defaults **ON**, plus a σ_loc bias correction. ✅
- **HMM** — emission steps `r/√Δt`-normalised by default (state σ/D handle variable Δt). ✅
  *Residual:* the transition matrix is discrete-time (per-step) → transition **rates** approximate
  under heavy dark-gapping. **Now disclosed** via a caveat in the HMM summary panel.
- **PELT/change-point** — segments on `log(r²/Δt)` (= local D), so the cost is already Δt-aware;
  per-segment D uses `mean(r²/Δt)` (avoids Jensen bias). ✅ Effectively GREEN as implemented.
- **Anomalous α** — fit on the **ensemble** (averaged) MSD, not naive per-track; `min points/TID`
  default 20 excludes very short tracks. **Now adds** a 95% CI (from the fit covariance, labelled a
  lower bound since MSD lags are correlated) + track count + a Kepten-2015 short-track caveat.

**Tier-1 (correctness) fixes applied:** α CI + caveat; HMM transition-rate caveat; D-map citation
corrected (was Ostersehlt 2022 → Vestergaard 2014 / Vogler 2025).

## 5. Action items & caveats surfaced by this review

1. **⚠ Citation fix — D-map.** The single-lag D-map estimator (`D = ⟨r²⟩/(4·median_Δt)`) is cited
   in-code to *Ostersehlt et al. 2022*, but that paper is a **DNA-PAINT imaging** paper with **no
   diffusion-mapping** component (verified). The single-lag method itself is sound (it is the CVE
   short-lag limit, Vestergaard 2014). **Recommend** re-attributing to Vestergaard 2014 / Vogler 2025,
   or describing it as a standard single-lag estimator, and removing the Ostersehlt citation.
2. **"R ≈ 0" is an inference, not a cited number.** No MINFLUX paper reports a numeric Berglund blur
   coefficient; R≈0 follows from short cycle + active re-centering. Fine to keep as a modelling
   assumption, but label it as such rather than as a published value. **Also note** SimuFLUX shows this
   breaks down at high D / coarse sampling, where MINFLUX under-estimates D — worth an info-box note
   near the fast-diffusion (D_max) tools.
3. **Δt-awareness audit.** Confirm JDD pooling, HMM transitions, and PELT segment costs operate on
   **true elapsed time / integrated steps**, not raw step index. CVE (integrated pairs) and cycle-based
   MSD already do this; the others should be checked.
4. **Short-track honesty for α.** Single-track anomalous-α is unreliable at MINFLUX track lengths
   (Kepten 2015). Prefer pooled/ensemble α with CIs and short-track exclusion; consider a UI caveat.
5. **✅ DONE — DCR crosstalk quantified** (superseded an earlier, mistaken "registration error" item —
   ratiometric MINFLUX has no inter-channel offset, see §6.2). `auto_detect_dcr_split` now returns d′,
   the expected mis-assignment rate, per-group contamination and the ambiguous-trace fraction, reported
   with a clean/workable/poor verdict and a format-aware hint about pooled DCR. Per-molecule confidence
   (`dcr_split_confidence`) identifies which individual traces sit on the edge of the split, with a plot
   and a per-trace export. **Follow-up not implemented — the strict gate:** actually *excluding*
   edge traces and re-running the cross-channel tools, so a co-localization signal can be shown to
   survive removal of the ambiguous band. The case for it is now quantified (edge-drop took
   mis-assignment 2.43% → 0.40% for 9% data loss), but it would change channel assignment, so it was
   left as an explicit choice rather than done silently.
9. **Split confidence is not exposed in the Display-tab colouring.** `dcr_split_confidence` is currently
   surfaced only in the 2-colour panel (plot + CSV). Colouring localizations by split confidence would
   make edge molecules visible in situ, but the point builders (`build_pts_2d` / 3D) use explicit
   column whitelists in six `summarise()` blocks, so it needs a deliberate pass rather than a patch.
6. **✅ DONE — crowding-vs-timing null.** Implemented as a within-file label permutation with a per-structure
   95% band and a Spearman trend test (§I-A.9), rather than the simulator null originally proposed: the
   permutation is exact and assumption-free, while a simulator null would import a bleaching model
   calibrated from this same kind of data. Validated on the collaborator file (coverage 0.052 vs nominal 0.05).
7. **✅ DONE — photon-weighted shape fitting.** `viewer_weight_eco` supplies w = Σ`eco` per trace to the
   `mle_fit_*` optimisers, with weighted RMSD and `n_eff` reported (§6.5). **Follow-up:** the resulting
   sensitivity (radius 359 → 149 nm on one test structure) suggests the viewer should report a
   **parameter uncertainty**, e.g. a bootstrap over traces, before any shape parameter is quoted as a
   measurement. Currently a single point estimate is shown with no error bar.
8. **Weighting is only as good as `eco`.** The photon weight assumes `eco` is proportional to collected
   photons and comparable across traces within an acquisition. That holds within one sequence/power
   setting; it is **not** safe across files acquired at different laser powers, so photon-weighted fits
   should not be pooled across such files without checking the `eco` distributions.

---

*Compiled with two literature-research passes; every citation's title/venue/year was verified via
search (DOIs cross-checked where shown). Items that could not be double-confirmed are flagged in §5.*

*Addendum 2026-07-24 — per-panel method boxes added to the SMLM Analysis tab.* This pass added §6.1a,
§6.2a, §I-A.10, the pitfall lists in §I-A.9 / §6.7, and action items 5–7. **Refs 54–58 (Costes 2004,
Lachmanovich 2003, Cleveland 1979, Barber/Dobkin/Huhdanpaa 1996, Pearson 1901) were added from standing
knowledge and have NOT been through the search-verification pass** that refs 1–53 received — they are
long-established, widely-cited works, but re-confirm the DOIs before quoting them in a manuscript.
