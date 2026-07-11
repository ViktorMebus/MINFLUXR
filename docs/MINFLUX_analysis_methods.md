# MINFLUX Tracking Analysis — Literature & Method-Appropriateness Review

*Reference document for the tracking-analysis module of Mini-MINFLUXR.*
Last reviewed: 2026-07 · Compiled from verified primary literature (see [Bibliography](#bibliography)).

This document (1) records **which analysis tools are ported from which researchers** so the work is
properly cited, and (2) gives an honest **appropriateness rating** of each method against the specific
statistical character of MINFLUX single-particle-tracking (SPT) data. Ratings:

- 🟢 **GREEN** — well-suited; MINFLUX's properties are neutral or an advantage.
- 🟡 **AMBER** — usable, but an assumption is strained; needs a MINFLUX-aware adaptation.
- 🔴 **RED** — a core assumption is genuinely violated; results are biased unless mitigated.

> A rating describes the *textbook method*. Where our implementation already adapts the method for
> MINFLUX, that is called out under **Our handling** — several tools are more defensible than their
> raw textbook rating suggests.

---

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
