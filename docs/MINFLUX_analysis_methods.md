# MINFLUX Analysis — Methods, Mathematics & Appropriateness

*Reference document for the analysis and simulation modules of Mini-MINFLUXR.*
Last reviewed: 2026-07 · Compiled from the source code and from verified primary literature (see [Bibliography](#bibliography)).

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
  - **`costes_nsim`** spatially randomizes Ch2's points uniformly in the shared window `nsim` times,
    recomputes r at the chosen bandwidth, and returns a one-sided **p-value** = (#{random r ≥ observed}+1)
    /(nsim+1) plus the null mean/SD — a significance test against spatial independence (p≈0.5 for
    independent channels, p→0 for co-localized). A Costes-style label-randomization null.
  - The UI points users to **g₁₂(r)** for a bandwidth-free, scale-resolved alternative.

### I-A.6 Radius of gyration & gyration-tensor anisotropy — [Part II §6.3]
- **Rg:** `Rg² = mean(‖p − p̄‖²)` per cluster (second moment about the centroid).
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
- **Why exploratory:** the shape is imposed by the dictionary and "best" is a model-selection outcome, not
  proof of biology; with few blink-inflated points several shapes are statistically indistinguishable.

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

## I-B. Tracking (SPT) analysis

### I-B.0 Cycle time & time-signature (the foundation for everything Δt-aware)
MINFLUX has no fixed clock; the base **cycle time t_c** is inferred (`trk_cycle_time`), in priority order:
(1) hardware model `t_c = N_iter · t_dwell` if enabled; (2) manual override; (3) **data-driven cluster
cycle** — cluster the inter-localization Δt values and take the smallest robust cluster centre (the modal
Δt is *misleading* because dark cycles inflate it); (4) a regression cycle from Δt-vs-index slope. The
resulting t_c drives gap-splitting, the quality gate, cycle-index MSD binning and the D-map. A step's
**cycle count** is `n_cyc = round(Δt / t_c)`; **integrated (single-cycle)** steps have n_cyc = 1.

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
  - **Confined:** `MSD = A(1 − e^{−τ/t_c}) + C` → local `D = A/(2d·t_c)`, confinement size **L = √(3A)**,
    plateau A; flagged *degenerate* if t_c > 2·max(τ) (no plateau resolved).
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
mean-of-ratios form avoids the Jensen bias of `mean(r²)/(2d·mean(Δt))` when r² and Δt covary. The pooled
per-segment log₁₀D distribution is then fit with a **BIC-selected Gaussian mixture** (`mclust`, up to
`max_g` states).

### I-B.4 Jump-distance distribution (JDD) — [Part II §2.3]
`jdd_bundle`. Fit the step-length histogram with an **N-component Rayleigh mixture**
**p(r) = Σ_k w_k · (r/σ_k²)·exp(−r²/2σ_k²)** by MLE (Nelder–Mead on the NLL; weights via softmax, σ via
log-parametrisation). MINFLUX handling:
- **√Δt normalization ON by default** — pool `r/√Δt` so steps of different cycle counts are comparable
  (states then map to D, not raw step size).
- **D per component:** `D = σ²/divisor` (divisor 2 in 2D / 3 in 3D).
- **σ_loc bias correction:** subtract `2·σ_loc²/dt0` from σ² before D (localization noise adds 2σ_loc² per
  coordinate), clamped ≥0. A **noise-floor jump distance** `√(2d)·σ_loc/√dt0` marks where peaks are
  noise-dominated.

### I-B.5 Diffusive-state HMM (Baum–Welch + Viterbi) — [Part II §2.4]
`hmm_bundle`. Per-step displacement `r` (√Δt-normalised by default) with a **Rayleigh emission**
`log p(r|state) = log r − 2log σ_k − r²/2σ_k²`, K = 2–4 states. **Baum–Welch EM** entirely in log-space
(vectorised forward/backward, log-sum-exp), ≤30 iterations with early-stop; **σ_k update =
√(Σγ·r²/(2·Σγ))** (Rayleigh MLE), transition matrix A from the accumulated log-ξ. Fit on the longest
`max_fit_trk` tracks (optional contiguous-prefix step cap for speed); **Viterbi** decodes *all* tracks at
full length. States sorted slow→fast by σ. **Caveat:** emissions are Δt-aware but the transition matrix is
**discrete-time (per-step)**, so transition *rates* are approximate under heavy dark-gapping (disclosed in
the summary panel).

### I-B.6 Confinement (Rg / Rmax) — [Part II §2.5]
`conf_bundle`. Per track: **Rg² = mean(‖p−p̄‖²)**, **Rmax = max‖p−p̄‖**, ratio **Rg/Rmax**. Optional
noise correction removes localization inflation: **Rg²_true = Rg²_meas − d·σ_loc²·(1 − 1/n)** (the 1−1/n
accounts for estimating the centroid from the same noisy points), clamped ≥0; Rmax (an extreme value) is
left uncorrected. Low Rg/Rmax + small extent ⇒ confined — but see the simulator round-trip caveat (free
diffusion over a short track can mimic the same extent).

### I-B.7 Directionality / turn-angle autocorrelation — [Part II §2.6]
`dir_bundle`. **Turn angle** between consecutive displacement vectors θ = acos( (Δr_i·Δr_{i+1}) /
(|Δr_i||Δr_{i+1}|) ), clamped to [−1,1]. **Velocity autocorrelation** of unit step vectors over lags
1…20: `AC(lag) = ⟨v̂_i · v̂_{i+lag}⟩`. Summaries: median angle, %forward (<60°), %backward (>120°),
AC(1). A **noise floor** RMS pure-noise step = 2·σ_loc; the fraction of steps below it quantifies how much
of the backward-peak / lag-1 anti-correlation is localization noise rather than real anti-persistence
(MINFLUX's small σ_loc largely removes this classic artefact).

### I-B.8 Speed & transport (path speed, drift velocity, local α) — [Part II — Mini-MINFLUXR original]
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
- **Implemented as:** `σ_loc² = −⟨Δxᵢ·Δxᵢ₊₁⟩` over **consecutive integrated (single-cycle) step
  pairs** only; also returns bias-corrected `D_cve`.
- **Source:** Vestergaard, Blainey & Flyvbjerg 2014.
- **MINFLUX caveat:** the CVE algebra is derived for a **constant** Δt. Applied naively across
  variable-Δt MINFLUX steps it biases D — hence the raw RED rating.
- **Our handling:** we **restrict to consecutive single-cycle pairs**, i.e. we pre-segment into
  constant-Δt sub-runs — precisely the mitigation the literature prescribes. ✅ This makes our CVE
  statistically valid on MINFLUX data at the cost of discarding dark-gap-adjacent steps. Note the
  σ_loc co-estimation is inherently *less* valuable for MINFLUX (σ_loc is already ~nm), so we use CVE
  primarily as an MSD-free precision/plausibility check, which is the appropriate role.

### 2.3 JDD — jump-distance / Rayleigh-mixture — 🟡 AMBER
- **Implemented as:** jump-distance distribution with an N-component Rayleigh mixture; user σ_loc and
  Δt-normalization controls.
- **Source:** Schütz, Schindler & Schmidt 1997; benchmarking in Weimann et al. 2013.
- **MINFLUX caveat:** small σ_loc reduces error-broadening of the Rayleigh peaks (an advantage), but
  jumps must be pooled at a **fixed real Δt**; mixing steps of different cycle counts pollutes the
  distribution, and short tracks give few jumps per state → unstable multi-component fits.
- **Our handling:** exposes Δt-normalization / step-mode controls — **confirm these pool by true Δt
  (integrated steps), not raw step index** (see §5 action items).

### 2.4 HMM diffusive-state segmentation — 🟡 AMBER
- **Implemented as:** per-track Gaussian-emission HMM over displacement steps, N states, with
  Δt-normalization option.
- **Source:** Das, Cairo & Coombs 2009; variational-Bayes short-track pooling (vbSPT) Persson et al. 2013.
- **MINFLUX caveat:** vbSPT's strength (pooling power across many *short* tracks) and MINFLUX's fine
  precision (better state separability) are a good match — but the discrete-time transition matrix
  assumes a **uniform Δt per step**. Irregular cycles need either resampling to a common grid or a
  continuous-time-rate transition model.
- **Our handling:** the `norm_dt` option addresses emission scaling; the **transition model remains
  discrete-time** — treat transition-rate outputs as approximate under heavy dark-gapping.

### 2.5 Confinement analysis — 🟢 GREEN
- **Implemented as:** sliding-window local-variance vs. free-diffusion expectation, noise-subtraction option.
- **Source:** Simson, Sheets & Jacobson 1995 (confinement index); lineage Saxton, Kusumi.
- **MINFLUX fit:** nm precision resolves confinement domains far below the camera-SMLM ~20–30 nm floor
  — a genuine strength. Only caveat: the window should be **time-weighted** (elapsed time), not
  point-count-weighted, for irregular Δt — a minor implementation detail, not a fundamental mismatch.

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
| CVE | Vestergaard 2014 | 🔴 | 🟡 (integrated-pair restriction) |
| JDD | Schütz 1997; Weimann 2013 | 🟡 | 🟡 (verify Δt pooling) |
| HMM / vbSPT | Das 2009; Persson 2013 | 🟡 | 🟡 (discrete-time transitions) |
| Confinement | Simson 1995 | 🟢 | 🟢 (time-weight window) |
| Turn-angle | Burov 2013 | 🟢 | 🟢 |
| Anomalous α | Metzler 2014; Kepten 2015 | 🔴 | 🔴 (pool + exclude short) |
| Change-point / PELT | Killick 2012 | 🟡 | 🟡 (Δt-aware cost) |
| D-map (single-lag) | see §5 citation note | 🟢 | 🟢 |

---

## 3. MINFLUX-specific parameter tools

- **D_max ceiling / trackable-D limit** — `D_max = (L²/16)·ψ / (N_PL + ψ·t_hw)`, cited in-code to
  **Vogler et al. 2025** (*Commun. Biol.*). Consistent with the independent SimuFLUX ceiling of
  ~3 µm²/s (Marin & Ries 2025). 🟢 Correctly attributed.
- **σ_loc noise floor / dark-cycle handling / continuity gate** — Mini-MINFLUXR-original engineering
  built on the data properties in §1; not a port of a specific paper.

---

## 3a. MSD on MINFLUX — computation, confinement & what is non-standard

The MSD tool follows **Béla Vogler's `the_road_to_minflux_msd` workflow** (Vogler et al. 2025), the
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

### 6.3 Radius of gyration / gyration-tensor anisotropy — 🟢 GREEN
- **Source:** standard second-moment cluster size; anisotropy eigenvalue coefficients Rudnick & Gaspari
  1986; linear/planar/spherical coefficients Westin et al. 1997/2002.

### 6.4 KDE density map + Pearson co-localization — 🟡 AMBER
- **Source:** kernel density estimation (Silverman 1986); Pearson correlation of two density surfaces as
  a co-localization index (cf. Dunn et al. 2011).
- **Caveat:** the Pearson r is strongly bandwidth-dependent (small bandwidth → r≈0 even for true
  co-localization; large → inflated). Now made explicit by the optional **r(bandwidth) sweep** plot, and
  the single-bandwidth r is accompanied by a **Costes-style spatial-randomization p-value** (§I-A.5) that
  tests it against spatial independence. For a bandwidth-free, scale-resolved alternative the UI points
  to the pair cross-correlation g₁₂(r) (§6.2). Still treat the raw r as a *relative* measure.

### 6.5 Cluster shape fitting (viewer) — 🟠 EXPLORATORY
- **Method (standard):** orthogonal-distance MLE + AICc/BIC model selection (Akaike 1974; Schwarz 1978;
  Hurvich & Tsai 1989).
- **Why exploratory:** the shape dictionary (sphere/ellipsoid/helix/ring/arc) is imposed by this tool and
  the "best" shape is a model-selection outcome, not proof of biological geometry; with few,
  blinking-inflated points several shapes are statistically indistinguishable.

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

### SMLM summary

| Tool | Source | Rating |
|---|---|---|
| NN / DBSCAN / HDBSCAN / OPTICS | Clark & Evans 1954; Ester 1996; Ankerst 1999; Campello 2013 | 🟢 (feed de-blinked) |
| Ripley L / L₁₂ / PCF g₁₂ | Ripley 1977; Besag 1977; Lotwick & Silverman 1982; Dixon 2002; Sengupta 2011; Veatch 2012 | 🟢 / 🟡 labelling |
| Radius of gyration / anisotropy | Rudnick & Gaspari 1986; Westin 1997/2002 | 🟢 |
| KDE + Pearson co-loc | Silverman 1986; Dunn 2011 | 🟡 bandwidth |
| Shape fitting | Akaike 1974; Schwarz 1978; Hurvich & Tsai 1989 | 🟠 exploratory |
| GT estimate & structure size | this tool (in-silico validated); Prim 1957 | 🟠 exploratory |
| Detection-time / crowding / NN-vs-size | this tool | 🟠 exploratory |

---

## 4. Bibliography

**MINFLUX microscopy & tracking**
1. Balzarotti F, Eilers Y, Gwosch KC, Gynnå AH, Westphal V, Stefani FD, Elf J, Hell SW (2017). *Nanometer resolution imaging and tracking of fluorescent molecules with minimal photon fluxes.* Science 355(6325):606–612. doi:10.1126/science.aak9913
2. Gwosch KC, Pape JK, Balzarotti F, Hoess P, Ellenberg J, Ries J, Hell SW (2020). *MINFLUX nanoscopy delivers 3D multicolor nanometer resolution in cells.* Nature Methods 17:217–224. doi:10.1038/s41592-019-0688-0
3. Schmidt R, Weihs T, Wurm CA, Jansen I, Rehman J, Sahl SJ, Hell SW (2021). *MINFLUX nanometer-scale 3D imaging and microsecond-range tracking on a common fluorescence microscope.* Nature Communications 12:1478. doi:10.1038/s41467-021-21652-z
4. Wolff JO, Scheiderer L, Engelhardt T, Engelhardt J, Matthias J, Hell SW (2023). *MINFLUX dissects the unimpeded walking of kinesin-1.* Science 379(6636):1004–1010. doi:10.1126/science.ade2650
5. Deguchi T, Iwanski MK, Schentarra EM, et al. (2023). *Direct observation of motor protein stepping in living cells using MINFLUX.* Science 379(6636):1010–1015. doi:10.1126/science.ade2676
6. Ostersehlt LM, Jans DC, Wittek A, et al. (2022). *DNA-PAINT MINFLUX nanoscopy.* Nature Methods 19:1072–1075. doi:10.1038/s41592-022-01577-1
7. Vogler BTL, De Angelis G, Zhao Z, Eggeling C, Reina F (2025). *Parameter optimization for MINFLUX microscopy enabled single-particle tracking.* Communications Biology 8:1573. doi:10.1038/s42003-025-09060-1
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

---

*Compiled with two literature-research passes; every citation's title/venue/year was verified via
search (DOIs cross-checked where shown). Items that could not be double-confirmed are flagged in §5.*
