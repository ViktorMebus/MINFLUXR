# CVE and σ_loc in MINFLUX tracking — verdict and plan

Status: **plan only. No `server_tracking.R` change has been made.** Written 2026-08-24 from the
GP67H tracking test data and the simulator round-trip built the same day
(see `simulation_module.md` §5 and §9).

---

## 1. Verdict in one paragraph

CVE is a sound estimator for **D** on MINFLUX tracking data, provided it is restricted to
single-cycle steps, and the simulator now confirms it is unbiased to ~1% across every geometry and
motion model the app offers. CVE is **not** a usable estimator for **σ_loc** on fast tracking data,
and the app's current form of it is additionally wrong by a factor of ~3. σ_loc should come from the
photon budget instead, with the covariance form kept only as a cross-check.

---

## 2. Evidence

### 2.1 D — CVE is fine

Simulator round-trip, 100.8 nm TCP, 38 photons, t_c 0.192 ms, SW 0.011 ms, D = 0.30 µm²/s:

| geometry | free | confined | directed | two-state |
|---|---|---|---|---|
| sphere | 0.305 | 0.302 | 0.305 | 0.156 |
| rod | 0.298 | 0.295 | 0.298 | 0.161 |
| ellipsoid | 0.304 | 0.301 | 0.304 | 0.156 |
| membrane | 0.311 | 0.308 | 0.311 | 0.164 |

(Two-state is *supposed* to halve — it is a population average over mobile and immobile.)

The conditions that make this hold are load-bearing:

- **Single-cycle steps only.** See §2.3 — multi-cycle steps are not just noisier, they are biased.
- **Δt = k·t_c + SW.** Using t_c alone over-reports D by SW/t_c ≈ 5.5%.
- **Free motion.** On confined tracks CVE returns the *local* D at single-cycle lags, which is the
  right answer to a different question.
- **Survivorship off.** Not yet swept; see §5.

### 2.2 σ_loc — CVE is structurally unable to answer

Vestergaard's σ estimator is σ² = R·⟨Δx²⟩ + (2R−1)·cov₁ per axis. On file 013:

```
⟨Δx²⟩ = 92.95 nm²          cov₁ = −2.38 nm²
R = 0      (app today)     σ² =        2.38   → σ = 1.54 nm
R = 0.087  (measured τ)    σ² =  8.09 + 1.97  → σ = 3.17 nm
R = 1/6    (τ = full cycle) σ² = 15.49 + 1.59  → σ = 4.13 nm
```

**Correction to an earlier draft of this document.** It reported σ = 4.13 nm and treated its agreement
with L/(4√N) as a confirmation of the prefactor 4. That agreement was contingent on assuming the
collection window fills the whole cycle. It does not: eco/efo says τ is 55% of the cycle, which moves
the covariance estimate to 3.1 nm and the implied prefactor to ≈ 5.3. **The prefactor is not pinned by
this data** — it sits somewhere between 4 and 5.3 depending on the blur model, which is the argument
for §6's immobile-sample calibration, not against it.

Three separate problems, in increasing order of severity:

1. **The app assumes the wrong R.** A MINFLUX localization integrates photons over its collection
   window, so R = τ/(6·Δt) > 0. **τ is measurable, not assumable: τ = eco/efo**, which on both GP67H
   files is 0.1061 ms of a 0.192 ms cycle — 55%, giving **R ≈ 0.087**. Using R = 0 under-reports
   σ_loc by ~2.4× (1.5 and 1.3 nm against **3.17 and 3.09 nm**). This is blur, not precision.
2. **The estimate is dominated by the assumption, not by the data.** At R = 1/6, 91% of σ² comes
   from the R·⟨Δx²⟩ term. σ̂ ≈ √R · rms_step is therefore close to a deterministic function of the
   step size and the assumed R. The actual noise signature — the part cov₁ carries — is a 4%
   residual (r₁ = −0.026). Whatever R is chosen, this number is not really measuring σ_loc.
3. **It inherits every motion-model error.** ⟨Δx²⟩ is inflated by population heterogeneity, by drift,
   and by any confinement structure, and all of that flows straight into σ̂.

The photon route has none of these properties: **σ_i = L/(a·√N_i)** with L = 100.8 nm from the
sequence `patGeoFactor` and N = `eco` is per-localization rather than pooled, and assumes nothing about
the motion. Its one free parameter is the prefactor a, and the covariance route brackets it at 4–5.3
rather than fixing it. `TCP_CRLB_PREFACTOR` in `packages.R` currently holds 4; settling it needs the
experiment in §6.

### 2.3 A separate finding: multi-cycle steps are badly biased

Mean per-axis ⟨Δx²⟩ against cycle count k (file 014, and 013 behaves the same):

| k | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| ⟨Δx²⟩ nm² | 93.2 | 131.3 | 138.1 | 142.5 | 152.0 | 153.0 |
| free-diffusion expectation | 93 | 176 | 259 | 342 | 425 | 508 |

Free Brownian motion demands linear growth. The observed growth is ~35% of expectation at k=2 and
essentially flat beyond. A naive slope fit over k returns D = 0.046–0.056 µm²/s, four times below the
single-cycle CVE value of 0.22.

Candidate mechanisms, all now implemented in the simulator and swept in `simulation_module.md` §10:

| mechanism | simulated k-curve (k = 1…6) | verdict |
|---|---|---|
| heterogeneous D alone | 100 100 105 107 96 144 | no effect |
| speed-dependent dark cycles | 81 81 88 103 115 114 | depresses all k alike, no gradient |
| post-gap estimator shrinkage | **118** 27 27 25 28 24 | inflates k = 1 — **rejected** |
| survivorship / censoring (R_reacq 50 nm) | 100 96 95 78 66 67 | matches the real shape |
| longer integration over gaps (τ×k, cap 3) | 92 80 75 79 81 94 | real, measured — explains part |
| *real data* | 100 75 53 42 36 30 | |

A fifth mechanism turned up while auditing the simulator against the raw file and is now implemented:
**a dark cycle is not a skipped cycle, it is a longer integration.** τ = eco/efo before a single-cycle
step is one repeat; before a two-cycle step it is two, in 98% of cases, saturating near three. The
longer window blurs away part of the displacement. It is a genuine contributor, but a joint blur-aware
fit using the measured per-step τ across k = 1…6 returns D = 0.05–0.06 µm²/s at R² = 0.78 — four times
below the single-cycle CVE D and a poor fit — so it is not the whole story either. τ saturates at
three repeats while free diffusion keeps growing, which is why the k ≥ 4 deficit survives it.

**Censoring is the surviving candidate.** It is the only mechanism that leaves k = 1 intact while
producing a gradient over k, and with a heterogeneous population it also reproduces the second
signature: the step-weighted mean D falls well below the per-track median, because fast molecules have
their tracks cut short and contribute fewer steps.

An earlier draft of this document argued censoring was excluded because tracks do not preferentially
end after a large gap (0.07% of k=1 and 0.6% of k=2 steps are a track's last step). **That test is
invalid** — a censored step is *absent* from the file, so it never appears as a track-ending step.

Whatever the cause, the practical consequence is immediate: **any analysis that admits multi-cycle
steps is biased downward on this data**, which includes the dark-tolerance continuity path and any MSD
lag built from gapped steps. This should be checked before the next round of D numbers is quoted.

---

## 3. Proposed tracking-tab changes

Ordered by value; each is independently shippable.

### Phase 0 — read what the file already tells you

Three quantities the app currently infers or assumes are written in the data:

| quantity | where it is | value on GP67H |
|---|---|---|
| collection window τ | **eco / efo**, per localization | 0.1061 ms, quantised in units of 0.1061 |
| blur coefficient R | τ / (6·Δt) | 0.087 |
| repeats used | τ / τ₀ | 1 before a k=1 step, 2 before k=2, saturating at 3 |

τ = eco/efo is exact and per-localization. Any estimator that needs R should take it from there
rather than assuming a value. This is a small change with a large effect: R = 0 versus R = 0.087
versus R = 1/6 moves the σ_loc estimate over a factor of 2.7.

### Phase 1 — make σ_loc photon-based (the actual fix)

- Add a **photon-CRLB σ_loc** as the app's primary estimate: σ_i = L/(4√N_i) per localization, with
  L pre-filled from the uploaded sequence .json (`trk_seq_params()$L_tcp_nm` already parses it) and
  N from `eco`. Falls back to `DMAX_L_TCP_DEFAULT` when no sequence is loaded.
- Re-point `trk_sigma_loc()`'s priority chain at it. Today the chain prefers the CVE value
  (`trk_sigma_loc_cve`), which is the number §2.2 says not to trust.
- Because it is per-localization, it can feed the heteroscedastic consumers directly rather than
  through a pooled median.

### Phase 2 — fix, demote and label the CVE σ

- Add the R term: σ² = R⟨Δx²⟩ + (2R−1)cov₁ with R = τ/(6Δt) and Δt = t_c + SW. **Take τ from the data
  as eco/efo per localization** — it is quantised and unambiguous — rather than from the sequence or a
  default. Fall back to τ = t_c only when efo is missing.
- Because τ grows over a gap, the single-cycle steps adjacent to a gap have mismatched windows, which
  the CVE assumes away. In the simulator this biases D high by ~7% when the repeat behaviour is on;
  expect the same on real data.
- Present it as a **cross-check on L**, not as the σ_loc: report the implied L = 4·σ̂·√N̄ next to the
  sequence's L. Agreement validates the photon model; a gap is a real finding about the acquisition.
- State the r₁ noise fraction next to it so its (lack of) resolving power is visible.

### Phase 3 — audit the multi-cycle steps

- Ship the k-vs-⟨Δx²⟩ table of §2.3 as a diagnostic panel. It is three lines of arithmetic and it
  exposes the single largest bias found in this data.
- Audit every consumer for whether it admits k>1 steps (MSD lags, JDD, D-map, dark-tolerance
  continuity), and gate the ones that do behind an explicit, defaulted-off toggle until the
  mechanism is known.

### Phase 4 — D itself

- Keep CVE as the D estimator. Add the Δt = t_c + SW correction if it is not already using the
  measured spacing rather than the modal Δt.
- Report a per-track D distribution alongside the pooled value. On these files the pooled CVE
  (0.21–0.22) and the per-track median (0.28–0.33) differ by 40%, because pooling weights by step
  count and the long tracks are the slow ones. Both are correct answers to different questions and
  only one of them is currently shown.

---

## 4. Simulator support — **done**

All four are implemented in the Diffusion / Tracking tab, and the summary panel now prints the
k-curve of §2.3 so any of them can be scored against the real files:

- **Heterogeneous D** — lognormal per track, CV input, median = the D field.
- **Re-acquisition shrinkage** — post-gap localization pulled back toward the previous reported one.
- **Speed-dependent dark cycles** — p_dark = 1 − (1 − p_dark₀)·exp(−d/scale).
- **Exposure fraction τ** — sets R = τ/(6Δt) instead of a hard-coded 1/6.

Verified: the CVE round-trip holds (98% of the step-weighted mean D) and σ_loc recovers 4.12 against a
truth of 4.09 nm across a τ sweep of 1.0 → 0.5 → 0.1 (R = 0.158 → 0.079 → 0.016), so the Phase-2
estimator of §3 does track the true R.

## 5. Sweeps still worth running

| Sweep | Question |
|---|---|
| R_reacq 20 → 300 nm × D | which R_reacq reproduces the real k-curve, and is it physically credible against the 323 nm search field? |
| censoring on/off at fixed D | how much apparent confinement does it manufacture in the Tracking tab's confinement tools? |
| heterogeneous D × censoring | does it reproduce the real 0.22 pooled vs 0.30 per-track gap quantitatively? |
| D sweep at fixed σ | at what D does CVE σ become unusable (r₁ → 0)? Predicted: already, here. |
| α sweep 0.6 → 1.0 at fixed D, σ | how much of the app's D deficit on GP67H is subdiffusion rather than slowness? |
| hop corral size / permeability | can an α-matched hop model be told apart from fBm by any statistic the app computes? |

## 6b. Estimator failure modes, measured in the simulator

Each of these was produced by simulating data with a known ground truth and running the app's own
estimator over it. They are properties of the estimator, so they transfer to real acquisitions with
the same character — and GP67H has all three characters.

| condition | effect on D | effect on σ_loc | why |
|---|---|---|---|
| **subdiffusion** (fBm, α = 0.78) | **≈ 65% of truth** | **inflated, 4.4 → 4.9 nm** | CVE reads the anti-correlated increments of subdiffusion as localization noise. The lag-1 covariance is meant to carry only σ², but α < 1 adds its own negative term |
| **wrong R** (τ assumed, not measured) | unaffected | 1.5 / 3.1 / 4.1 nm at R = 0 / 0.087 / 1/6 | σ² = R⟨Δx²⟩ + (2R−1)cov₁ is dominated by the R term whenever D·Δt ≫ σ² |
| **τ growing over a gap** | **≈ +7%** | small | single-cycle steps next to a gap have mismatched averaging windows; the CVE assumes they are equal |
| population heterogeneity | pooled D ≠ per-track median | inflated | ⟨Δx²⟩ pools molecules of different D |

**The three σ_loc errors stack, and they all point the same way.** On GP67H the app reports 1.3–1.5 nm.
Correcting R to the measured value takes it to ≈ 3.1 nm. Correcting for α = 0.78 subdiffusion would
take it lower again, ≈ 2.7 nm. Each correction is itself model-dependent, which is the point: the
covariance route cannot deliver a trustworthy σ_loc on data like this, at any level of care.

### So what *should* estimate σ_loc

1. **Photon CRLB, σ_i = L/(a·√N_i).** L from the sequence, N from `eco`. Per-localization, no motion
   model, immune to every row of the table above. This is the one to ship as primary.
2. **Pin `a` with an immobile acquisition** (§6). This is the only measurement in the chain that is
   motion-free. Until it exists, `a` is bracketed at 4 (theory, current `TCP_CRLB_PREFACTOR`) to about
   5.3–6 (covariance route with measured τ, corrected for subdiffusion), and the app should say so
   rather than print a bare number.
3. **Keep the covariance σ as a labelled cross-check only**, with R from eco/efo, and print the r₁
   noise fraction next to it so its resolving power is visible. Never let it feed `trk_sigma_loc()`.

Because eco varies within a track (brightness fades to ~50% by a molecule's last decile), σ_loc is
**not constant along a track**. The photon route gives that for free; a single pooled σ cannot.

## 6. Experiment worth requesting

None of the above measures σ_loc **directly**. The clean measurement is an **immobile sample acquired
with the tracking sequence** — fixed cells or fiducials under `Tracking_2D_Oct2022`, where the
within-trace spread *is* σ_loc with no motion model in the way. Fitting σ against `eco` over that
acquisition gives the constant in σ = a/√N empirically, which is the one number the entire photon
route rests on. The existing fiducial pair is imaged with an SMLM sequence (L = 39.6 nm), so it
calibrates the imaging iteration, not the tracking one.

---

## 7. Mode-invariance, measured on GP67H (added 2026-08-27)

**Question:** does the photon route work for directed, free *and* confined motion, given these tracks
are confined?

**They are confined.** Per-track TAMSD(k=100)/TAMSD(k=10) against the free-diffusion expectation of
10x: median **0.35** (013) and **0.37** (014), 10th-90th percentile 0.13-1.09. Strongly sub-linear,
consistent with the alpha < 1 of section 6b.

**Answer: yes, and it is the only one of the options that does.** Per track, on the 84 (013) and 200
(014) tracks with >= 120 localizations and >= 40 single-cycle steps, split into mobility terciles:

| 013 | D (nm2/ms) | confinement | sigma_CVE | sigma_phot |
|---|---|---|---|---|
| slow | 201.5 | 0.37 | 3.33 | 4.20 |
| mid  | 282.7 | 0.28 | 4.47 | 4.68 |
| fast | 402.7 | 0.39 | 4.78 | 4.76 |
| **swing slow->fast** | 2.0x | | **+44%** | **+13%** |

014 behaves the same: sigma_CVE +34%, sigma_phot +7%. Spread across tracks: **CV 37% for sigma_CVE
against 9-14% for sigma_phot**.

**The CVE sigma is tracking sqrt(D), not precision.** The blur term gives sigma^2 ~ R<dx^2> ~ 2R*D*dt,
so sigma ~ sqrt(D). D doubles across the terciles (201 -> 403), predicting a 1.41x rise; the observed
sigma_CVE rise is **1.44x** on 013 and 1.34x on 014. It is a mobility readout wearing a precision
label, and on confined tracks it is contaminated a second time because anti-correlated increments
mimic localization noise (section 6b).

**sigma_phot is not perfectly flat either, and should not be.** Its residual 7-13% swing is the real
mobility-photon coupling: a faster molecule leaves the TCP sooner, collects fewer photons, and *is*
localized less precisely. Spearman rho[sigma_phot, ECO] = -1.00 by construction; the D and confinement
correlations (+0.28, -0.26) are that coupling seen indirectly, not a motion-model dependence.

**Caveat if Vogler's full Eq. 1 is adopted.** The x-term makes precision depend on the emitter's offset
from the TCP centre, x = sqrt(2*D*t_loc), which *is* mode-dependent: smaller under confinement, larger
under directed motion. That dependence is physically correct rather than a bias, and on this data it is
worth +0.6%. The simple L/(a*sqrt(N)) form has no motion dependence at all.

**Consequence.** The prefactor a is still unpinned (section 6), but a is a common scale factor - it
shifts free, confined and directed alike and cannot break mode-invariance. So the immobile-sample
calibration and the mode-invariance property are independent: the route is already the right one for
confined data, the experiment only fixes its scale.

---

## 8. The increments are smoothed - independent confirmation of damping (2026-08-27)

A subdiffusive walk has anti-correlated increments *before* any localization noise. For fBm with
Hurst H (alpha = 2H), fractional Gaussian noise gives cov1/var = 2^(2H-1) - 1. Localization noise and
blur can only make r1 **more** negative. So for any combination of subdiffusion + noise + blur:

    r1_observed  <=  r1_motion-only  <  0

Measured on single-cycle steps, with alpha fitted on cycle-lags 1-10 (dark-gap pairs excluded):

| file | alpha | H | r1 predicted by motion alone | r1 observed | |
|---|---|---|---|---|---|
| 013 | 0.799 | 0.400 | **-0.130** | **-0.026** | 5.0x too weak |
| 014 | 0.774 | 0.387 | **-0.145** | **-0.019** | 7.6x too weak |

(The alpha values independently reproduce section 6b's 0.78.)

**The inequality is violated, and not marginally.** The observed anti-correlation is 5-7x weaker than
the motion alone demands. No choice of sigma_loc, R, or blur model can produce that - all of them push
r1 the other way. Something is injecting *positive* increment correlation, which is what a smoothing
filter does. `damping: 2`.

This is a second, independent line of evidence for the damping hypothesis of
`sigma_loc_estimation.md` section 1b, reached without touching the sigma ceiling argument.

### What it rules out

It is not only that the CVE sigma is imprecise on this data. **The increments have been filtered
before the file was written**, so the information the covariance route needs is not present to be
recovered. That kills the obvious "more correct model" upgrades:

- an **alpha-aware CVE** (correcting cov1 for the fGn term) would fit a motion model to filtered
  increments and return a confident wrong answer;
- likewise **MLE**, for the same reason - a better estimator of a quantity the data no longer carries.

The fix is at acquisition (DMP = 0), not in the estimator.

### Cheap QC that would have caught this

r1 and alpha are both already computed. Flagging `r1 > 2^(alpha-1) - 1` (i.e. less anti-correlated
than the fitted alpha demands) detects instrument smoothing directly, and costs two lines. Worth
shipping regardless of what happens to the sigma model.

---

## 9. Implemented 2026-08-27

**Vogler Eq. 1 is now the sigma_loc model.** `.crlb_doughnut(L, nu, N, x)` in `server_tracking.R`;
verified against their Table 2 at L = 150 nm, x = 0: PL 10/20/50 give 20.03 / 14.17 / 8.96 nm against
the paper's 20.03 / 14.35 / 9.08 (within 1.3%).

- `TCP_BEAM_FWHM_NM <- 195.5` in `packages.R`, exposed as **Excitation doughnut FWHM nu** under
  Acquisition Parameters. Fitted from Vogler, **not measured on this microscope**.
- x = sqrt(2 D t_loc) from the CVE D, so the estimate is mode-aware in the physically correct
  direction (smaller under confinement, larger under directed motion).
- Degrades to L/(4 sqrt(N)) if nu is cleared - the nu -> infinity limit.

**Effect on GP67H: sigma_loc 4.09 -> 5.04 nm (013) and 4.20 -> 5.18 nm (014), +23%.** Every downstream
noise floor moves with it: D-map, JDD bias term, confinement correction, HMM noise state,
directionality and speed. Results quoted before this change are not comparable.

**The app now prints a value its own covariance check calls impossible** (section 1b of
`sigma_loc_estimation.md`: the data caps sigma at 4.13 nm). That is intended - the disagreement is the
finding, and section 8 says why the covariance side is the untrustworthy one here.

**Smoothing check shipped.** The status line reports alpha (from cycle-lags 1 and 8), the fGn bound
r1 <= 2^(alpha-1) - 1, and the observed r1, flagging a violation as probable galvo damping.

Still open: nu is transferred, not measured. The immobile run should be done **at two or three
different L** so nu is over-determined and Eq. 1's L-dependence is tested rather than assumed.

---

## 10. The fixed-cell calibration, measured (2026-09-02)

24 files, `Test data/Fixed Tracking datasets/`, fixed E. coli + 4 nM JF646, three strains
(GP67H, GP170H, HaloCtrl) across TrkSeq2D (DMP 2) and cal_L101 / L151 / L302 (DMP 0).
Naming verified against the data: iteration count identifies the sequence (4/3/2 -> L101/L151/L302)
and cfr ceiling separates CFR-on from CFR-off. **All 24 filenames are consistent with their contents.**

### 10.1 The sample is not stationary, and it is the tracker

MSD on a fixed sample should be FLAT at 4 sigma^2. It is not:

| L | sigma from MSD intercept | wander D_w (nm^2/ms) |
|---|---|---|
| 100.8 | 5.1 nm | 83-186 |
| 151.2 | 7.2 nm | 330-508 |
| 302.4 | fit fails | 10300-12300 |

r1 of the increments is **positive** (+0.05 to +0.18). For a stationary emitter with independent
noise r1 = -0.5 exactly, and any damped feedback loop gives r1 in [-0.5, 0]. Positive r1 is
impossible for a stationary process - so the reported position genuinely wanders.

It is **not biology**: the three strains agree at the same L (within ~10%), and D_w scales
super-quadratically with L (80 -> 400 -> 11000 across L = 101 -> 151 -> 302). No biological
process tracks the TCP diameter. It is the loop re-centring on its own estimation noise.

### 10.2 sigma_loc: the prefactor is ~3.7, and there is NO L-dependence

Fitting MSD(tau) = 4 sigma^2 + 4 D_w tau on lags 1-5 (the curve saturates, so a 1-20 fit
extrapolates a concave curve and returns an intercept above MSD(1) - impossible):

| L | a = sigma*sqrt(N)/L, intercept | a, from MSD(1) (upper bound) |
|---|---|---|
| 100.8 | **0.265** | 0.336 |
| 151.2 | **0.272** | 0.360 |

**measured a(151)/a(101) = 1.03.** Vogler Eq. 1 at nu = 195.5 predicts **1.39**; the simple
`L/(4 sqrt(N))` predicts **1.00**. The flatness holds on the upper-bound column too (ratio 1.07),
so it is not an artefact of how the wander was handled.

**Conclusions:**

1. **The simple form wins on this instrument.** Eq. 1's L-dependence is not present. nu = 195.5 nm,
   transferred from Vogler's microscope, does not describe this one - effectively nu is large and the
   `nu -> infinity` limit is the right model here.
2. **The prefactor is 1/0.27 = 3.7**, against the app's `TCP_CRLB_PREFACTOR = 4`. The app is ~8% low.
   That is within the spread across strains and files, so 4 is defensible; 3.7 is the point estimate.
3. **L302 is unusable** - the loop cannot hold at that diameter (D_w ~ 11000 nm^2/ms, intercept fit
   collapses). As predicted for a diameter past Eq. 1's singularity.

### 10.3 The damping cost, measured

Same L = 100.8, same strains, DMP 2 vs DMP 0:

| | sigma (nm) | D_w (nm^2/ms) |
|---|---|---|
| TrkSeq2D, DMP 2 | 4.3-4.8 | 70-125 |
| cal_L101, DMP 0 | 4.7-6.6 | 83-186 |

Damping suppresses the apparent sigma by ~13% and the wander by roughly 2x - exactly the direction
section 8 inferred from the r1 violation, now measured directly.

### 10.4 The finding that matters most for the biology

**D_w on FIXED cells with the live-cell sequence (DMP 2, L = 100.8) is 70-125 nm^2/ms
= 0.07-0.125 um^2/s.** The live GP67H D_CVE is 0.215 um^2/s. So an instrumental floor accounts for
**roughly a third to a half** of the diffusion coefficient measured on the live data.

This is not a sigma_loc issue - it is a D issue, and it is larger than any of the sigma corrections
argued about in sections 1-9. Any quoted D from this instrument needs this floor subtracted or at
minimum reported alongside. It should be re-measured per sequence, since D_w scales steeply with L.

---

## 11. Through the cleaning pipeline, and what shipped (2026-09-02)

Section 10 analysed raw final-iteration data. Re-run through the app's cleaning - split at time gaps
(>7 cycles), trim dark-gap-separated ends (>=3 cycles), continuity gate (H=3, >=0.6), and long
segments only (>=200 localizations, 60-90% of localizations retained) - the numbers barely move:

| | raw | cleaned |
|---|---|---|
| a(151)/a(101) | 1.03 | **1.04** |
| prefactor 1/a | 3.7 | **3.62** |

So the conclusion is not an artefact of dirty traces. Per L, calibration sequences, excluding
HaloCtrl (whose larger wander makes its intercept extrapolation unreliable): a = 0.276 (L=100.8),
0.304 (L=151.2).

L302 survives cleaning with only **3% of localizations** - confirmed unusable.

### Shipped

- `TCP_CRLB_PREFACTOR` 4 -> **3.6** (measured, not assumed). The app ran ~10% low.
- `TCP_BEAM_FWHM_NM` 195.5 -> **0 = off**. Eq. 1's beam-size denominator predicts a 39% rise in
  a from L=101 to L=151; the measurement says 4%. nu does not transfer from Vogler's instrument.
- `.crlb_doughnut()` now keeps the **moving-emitter term when nu is off**: the large-beam limit is
  `sigma = L/(a sqrt(N)) * [1 + 2(x/L)^2]`, x = sqrt(2 D t_loc). Verified: a = 0.278 at both L, and
  the x-term adds +2% at D = 0.21 um2/s (GP67H), +8% at 1.0, **+39% at 5.0** (a fast cytosolic
  protein). That is the mode-awareness - smaller under confinement, larger under fast/directed motion.

### The mounting caveat on D_w

The samples are on a 2% agarose pad in a Gene Frame, not adsorbed to a coverslip, so part of the
wander may be genuine sample movement rather than tracker jitter. The L-scaling (D_w 160 -> 408 ->
6075 nm2/ms across L = 101 -> 151 -> 302) is far too steep for agarose creep, so an instrumental
component certainly exists at large L. But at L = 100.8 an L-independent mounting term could dominate,
and these data cannot separate the two.

**Either way it is a floor on D, not on sigma**, and it applies to the live-cell data because the
mounting is the same. A rigid target (surface-adsorbed dye or beads) run with the same sequences
would separate mounting from instrument.

---

## 12. Rigid beads on a coverslip - the wander is the tracker (2026-09-02)

8 files, `Test data/Fixed beads dataset/`, 40 nm red beads dried on a coverslip (no agarose), same
cleaning pipeline and short-lag MSD intercept as sections 10-11.

### 12.1 L302 is dead

Nothing survives cleaning on either L302 bead file. It failed on cells (3% retained) and it fails on
a bright, rigid, non-bleaching target - so the 2-iteration 302 nm sequence is broken in itself, not
starved of photons. **Drop it.** The earlier prediction that a bead would rescue it was wrong.

### 12.2 The stage is excellent, and the wander is NOT drift

Per-bead drift velocity from a linear fit over traces of >=2000 localizations, spans up to 10 min:

| file | beads | median drift | direction concentration R | random R |
|---|---|---|---|---|
| L101 CalL101B_002 | 33 | **0.2 nm/s** | 0.71 | 0.17 |
| L151 CAL151B_002 | 45 | 0.3 nm/s | 0.07 | 0.15 |
| L151 L151B_006 | 64 | 0.7 nm/s | 0.05 | 0.12 |
| L151 L151B_007 | 43 | 0.3 nm/s | 0.30 | 0.15 |

**Drift is 0.2-0.7 nm/s** - over a 9 s trace that is 2-6 nm total, below sigma. Negligible.

And the directions are **incoherent** (R at or near the random expectation in three of four files; the
one coherent case, R = 0.71, has a mean vector of only 0.3 nm/s, so it is coherent and irrelevant).
Beads do not move together. **The wander is generated per target, by the tracking loop - not by the
stage and not by the sample mount.**

### 12.3 The instrument floor

| target / mounting | L = 100.8 | L = 151.2 |
|---|---|---|
| beads on coverslip (rigid) | **D_w = 36 nm^2/ms** | 395 |
| fixed cells on agarose | 160 | 408 |

At L = 100.8 the tracker's own floor is 36 nm^2/ms = **0.036 um^2/s**; the fixed cells show 160, so
~124 is contributed by the agarose mount, incomplete fixation, single-molecule photophysics, or some
mixture. These data cannot split those further.

**Against a live GP67H D of 0.215 um^2/s: the instrument floor alone is ~17%, and the full fixed-cell
floor would be ~74%.** The latter is an upper bound - it includes anything the fixation failed to
immobilise.

### 12.4 Do NOT take the prefactor from beads

Bead sigma comes out *higher* than cell sigma at the same L and comparable ECO: a = 0.389 / 0.555
(L101 / L151) against 0.276 / 0.304 from the cells. Two reasons not to trust it:

- a **40 nm bead is not a point emitter** against a 100.8 nm TCP - 40% of the pattern diameter.
- at L151 the wander is so large (D_w ~ 400) that extrapolating lags 1-5 back to tau = 0 is a long
  extrapolation, and it inflates the intercept.

**The cell-derived `TCP_CRLB_PREFACTOR = 3.6` stands.** Beads answered the wander question, not the
sigma question, and that is the right division of labour: a rigid target tells you about the
instrument, a real sample tells you about the measurement.
