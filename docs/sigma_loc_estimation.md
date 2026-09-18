# Estimating σ_loc for MINFLUX tracking — options, and what MINFLUX sampling does to each

> **Superseded (2026-09-17) by [sigma_loc_model.md](sigma_loc_model.md)**, which documents the model the
> app ships and its immobile-sample calibration. This file is the earlier survey of estimation options
> (2026-08-25) and its status line below is historical; for the fixed-cell and bead measurements see
> [tracking_sigma_and_cve_plan.md](tracking_sigma_and_cve_plan.md) §10–12.

**Status (2026-08-25):** the app uses σ = L_TCP/(4√N). The literature says that form is only a
limiting case and should read larger; the GP67H data says it cannot be larger. §1b resolves the
contradiction — those files were acquired with galvo damping on, which suppresses exactly the step
variance every motion-based estimator reads. **Conclusion: change nothing on this data.** No
motion-based estimator can settle σ, and one immobile-sample run with a DMP = 0 sequence settles both
the prefactor and the size of the damping bias.

## Resume here — where this was left

**Nothing in the code was changed on the strength of this analysis, deliberately.** The app still uses
σ = L_TCP/(4√N), fed by `TCP_CRLB_PREFACTOR` in `packages.R`.

What to do, in order:

1. **Acquire a fixed/immobile sample with `Tracking_2D_L150_fastcycle.json`** (or the 3D one). Those
   set `damping: 0`, which the existing GP67H files do not. This single run answers both open
   questions — see §1b for why damping is the blocker.
2. **Fit σ against 1/√ECO** on that run. The slope *is* the prefactor, measured rather than assumed,
   and it makes ν irrelevant. Use short windows or extrapolate spread vs. window length so stage drift
   does not inflate it.
3. **Compare D on the same sample against a damped run.** That quantifies how much `damping: 2` has
   been costing the existing tracking results — D = 0.215 / 0.221 µm²/s on GP67H is probably low.
4. Only then decide whether to adopt the full Eq. 1 (§1) or keep the simple form.

Open, and not resolvable from the data in hand:

- The prefactor: 4 (ours), 2.35 (Vogler's own table), or ~1.64 (Imspector's `_lemd`). The data caps σ
  at 4.13 nm, which excludes the latter two — but the data is damped.
- ν for this microscope. Fitted at 195.5 nm from Vogler's table and transferred; step 2 makes it moot.
- What `_lemd` / `_lemn` in the sequence JSON actually mean. Worth one question to Abberior.

---

Sources: Vogler et al. 2025, *Commun. Biol.* 8:1573 (`Test data/Litterature/s42003-025-09060-1.pdf`);
Fitzpatrick et al. 2025, *Nanoscale* 17:16609 (`d5nr00948k.pdf`); Vestergaard 2014; Michalet 2010;
Berglund 2010.

---

## 0. The one thing to understand first

Every **motion-based** estimator (CVE, MSD intercept, MLE) sees σ_loc only through the same channel:
the negative lag-1 correlation that localization noise imprints on consecutive displacements. On this
data that channel carries **r₁ ≈ −0.02**, i.e. about 4% of the step variance. The other ~90%+ of the
fitted σ² comes from the motion-blur term R⟨Δx²⟩, where R is *assumed*, not measured.

So the choice between CVE, MSD and MLE is a choice of how to process a 4% residual. It is not a choice
that creates information. **The estimator is not the problem; identifiability is.**

The **photon-based** route (CRLB from ECO) uses a completely different channel — the photon count —
and is therefore the only genuinely independent option short of measuring precision directly.

---

## 1. What Vogler actually report, and what it implies about our prefactor

Their Eq. 1, the doughnut CRLB:

```
σ_CRLB(x) = L/(4√N) · [1 + 2(x/L)²] · |cosh(a) − (x/L)·sinh(a)| / |1 + ln2·(L/ν)²·((2x/L)² − 1)|
                                       a = 4·ln2·x·L/ν²
```

`L` = TCP diameter, `ν` = **FWHM of the excitation doughnut**, `N` = photons, `x` = emitter offset from
the TCP centre.

`σ = L/(4√N)` is the **ν → ∞, x = 0 limit** — it drops the denominator, which is always < 1 in the
practical regime. Their Table 2 lists σ_CRLB against PL at their final iteration L = 150 nm:

| PL | σ_CRLB listed | implied a in σ = a/√N | implied prefactor L/a |
|---|---|---|---|
| 10 | 20.03 nm | 63.3 | 2.37 |
| 20 | 14.35 nm | 64.2 | 2.34 |
| 50 | 9.08 nm | 64.2 | 2.34 |

**Their own numbers correspond to σ = L/(2.35·√N), not L/(4·√N).** Solving Eq. 1 at x = 0 for ν
reproduces the table with **ν ≈ 195 nm** — a plausible doughnut FWHM at 642 nm:

```
ν = 195.5 nm  ->  PL=10: 20.03 (listed 20.03),  PL=20: 14.16 (14.35),  PL=50: 8.96 (9.08)
ν = 360   nm  ->  PL=10: 13.48  -- does not fit
```

(So the 360 nm in `patGeoFactor × 360` is a base diameter, **not** ν. Both facts stand independently.)

### Consequence for our data

GP67H, L = 100.8 nm, ν = 195.5 nm transferred from the fit above:

| estimate | 013 (N=38) | 014 (N=36) |
|---|---|---|
| blur-aware CVE, R range | 3.84–4.03 nm | 3.58–3.98 nm |
| **app now, L/(4√N)** | **4.09 nm** | **4.20 nm** |
| Eq. 1 at x = 0 (immobile) | 5.01 nm | 5.15 nm |
| Eq. 1 at x = √(2D·t_loc) (moving) | 5.04 nm | 5.18 nm |

**The app is ~20% low**, and the CVE is lower still. The moving-emitter term adds only ~0.6% here
(x ≈ 9.4 nm against L = 100.8 nm) — it matters for faster D or smaller L, not for this dataset.

⚠ ν = 195.5 nm is *fitted from their table assuming L = 150 nm exactly*, then transferred to a
different microscope. That transfer is the weakest link in the chain above.

---

## 1b. The data refuses both literature values — and says why

Given the measured ⟨Δx²⟩ and cov₁, σ is fixed once R is chosen:
`R = (σ² + cov₁)/(⟨Δx²⟩ + 2·cov₁)`. R is physically capped at **1/6** (collection cannot outlast the
cycle). So each candidate σ implies a required R, and that is testable:

| candidate σ (013) | value | required R | verdict |
|---|---|---|---|
| app now, L/(4√N) | 4.09 nm | 0.162 | possible |
| Vogler Eq. 1, ν = 195.5 | 5.01 nm | 0.258 | **impossible** |
| Imspector `_lemd` rescaled | 9.76 nm | 1.054 | **impossible** |

**The maximum σ this data can support is 4.13 nm (013) / 4.08 nm (014)** — and `L/(4√N)` lands
exactly there. D is unaffected throughout: it is blur-independent.

### The likely reason, and it is actionable

`Tracking_2D_Oct2022.json` sets **`damping: 2`**. Vogler report that DMP > 0 "led to an
underestimation of particle motility" and kept it at 0 throughout their work. Damping smooths
consecutive position estimates, which suppresses the step variance — and every motion-based estimator
reads σ and D out of exactly that variance.

So the GP67H files plausibly under-report **both** σ_loc and D, and the ceiling above is an artefact of
the acquisition rather than a property of the instrument. This is not resolvable by choosing a better
estimator.

> The derived `*_L150_fastcycle` sequences already set `damping: 0`. A fixed-sample run with one of
> those is therefore the clean measurement for *both* questions at once: it calibrates the prefactor,
> and it quantifies how much damping was costing the existing data.

---

## 2. The options

### A. CVE (Vestergaard 2014) — current cross-check

`σ² = R⟨Δx²⟩ + (2R−1)·cov₁`, `D = ⟨Δx²⟩/(2Δt) + cov₁/Δt`

| | |
|---|---|
| **Pro (MINFLUX)** | Δt is hard-quantised, so integrated steps give *exactly* constant Δt — cleaner than camera SPT. ~10⁵ step pairs, so negligible statistical error. D is blur-independent and provably good (simulator round-trip recovers 0.30 from 0.30 to ~1%). |
| **Con** | σ² is ~90% the assumed-R term. R is bounded, not known: 0.137 (duty-cycle staircase) to 0.157 (full window). Population heterogeneity and drift inflate ⟨Δx²⟩ and flow straight into σ. |
| **Verdict** | **Keep for D. Never for σ.** |

### B. MSD intercept via OLSF — what Vogler use for σ_eMSD

`MSD(τ) = 2d·D·τ + 2d·σ² − 4d·D·R·Δt` — σ from the intercept.

The same R appears, so it inherits the identical confound, plus sensitivity to curve shape. Their own
Table 2 lets us test it: σ_eMSD/σ_CRLB across their 13 conditions runs **0.51 to 1.64**, mean ≈ 1.23.

| | |
|---|---|
| **Pro** | Independent of the CVE's covariance path; already implemented (`trk_sigma_loc_rv`); ensemble averaging over many lags. |
| **Con** | A 3× spread against the CRLB across acquisition settings. Confinement or anomaly bends the curve and biases the intercept. It is a small difference of large numbers. Damping (DMP > 0) visibly suppresses it — their two DMP > 0 rows give the lowest ratios. |
| **Verdict** | **Third cross-check at best.** Its disagreement with the CRLB is diagnostic, not a measurement. |

### C. MLE (Berglund 2010; Vestergaard's MLE variant)

Full likelihood for a blurred, noisy Brownian trajectory; jointly estimates D and σ².

| | |
|---|---|
| **Pro** | Statistically efficient — attains the CRB where CVE does not, which matters at low SNR and short tracks. Handles the correlation structure exactly. |
| **Con** | **Does not fix identifiability.** The likelihood is parameterised by (D, σ², R); R is not identifiable from the data and must be fixed by assumption, exactly as in A and B. Assumes constant Δt (integrated steps only) and free diffusion. Materially more code and compute for a σ that remains as assumption-bound as the CVE's. |
| **Verdict** | **Not worth implementing for σ.** Would be a modest gain for D, and D is already fine. |

### D. Photon-budget CRLB — current primary

`σ_i = L/(4√N_i)` today; should be the full Eq. 1.

| | |
|---|---|
| **Pro (MINFLUX)** | Uses a channel the motion estimators cannot see. Per-localization, not pooled. Assumes nothing about the motion, so heterogeneity, drift and confinement do not touch it. MINFLUX reports ECO per localization — camera SPT has no comparable per-frame photon count. |
| **Con** | It is a **bound**: real σ ≥ CRLB. Ignores background (`bgcThreshold`), imperfect TCP centring, and estimator sub-optimality. Needs ν, an instrument constant we do not have. Our current form additionally drops the denominator and is ~20% below even the bound. |
| **Verdict** | **Best default — but relabel it as a lower bound and add the denominator.** |

### E. Direct measurement on an immobile sample, tracking sequence

Within-trace spread with D = 0. No motion model anywhere in the chain.

| | |
|---|---|
| **Pro** | Measures the achieved precision, including background, centring and estimator losses. Calibrates the prefactor (equivalently ν) in one number. Removes the transfer assumption in §1. |
| **Con** | Needs an acquisition. Stage drift inflates within-trace SD over long traces — use short windows, or fit spread vs. window length and extrapolate to zero. Must use the **tracking** sequence: the existing fiducial pair is an imaging sequence (L = 39.6 nm), a different iteration. |
| **Verdict** | **The thing that actually settles it.** |

---

## 3. What MINFLUX sampling does to all of this

| Feature | Effect on motion-based σ (A/B/C) | Effect on photon-based σ (D) |
|---|---|---|
| **Dark cycles** | Break constant Δt, so only integrated steps qualify. That is a *selection*: dark cycles correlate with dim and fast molecules, so the surviving subset is biased bright and slow. | None — per-localization. |
| **Hard Δt quantisation** (Δt = k·t_c + SW) | Helps: lags are exact, no jitter. But t_c ≠ Δt, and R = t_c/(6Δt) depends on that offset. | None. |
| **Duty cycle** (collection only during dwells) | The crux. Sets R below the 1/6 ideal — 0.137 here — and R is ~90% of the fitted σ². | None. |
| **EOD movement** | Spreads collection across the cycle, keeping S(t) near-linear (R near 1/6), but each move is dead time. Vertex count is a weak lever: triangle 0.131 → hexagon 0.137. | None. |
| **Terminal localization step** | The strong lever: 0.165 without it, 0.137 with. ~17% on R, ~9% on σ. | None. |
| **PL termination** | Cycle ends at the photon limit, so N clusters near PL — but measured ECO is the honest count and runs above or below PL depending on settings. | Use **ECO**, not PL. We do. |
| **Re-centring each cycle** | The emitter offset x stays small while the loop keeps lock, so the CRLB stays near its x = 0 value. It degrades as D rises toward D_max — Vogler's central point. | Enters through the `[1+2(x/L)²]` term; +0.6% here. |

---

## 4. Recommendation

1. **Keep D from the CVE.** Verified, blur-independent, nothing to change.
2. **Keep σ from the photon budget, but fix and relabel it.** Add the Eq. 1 denominator with ν as an
   instrument parameter, and call it what it is: **a lower bound (CRLB)**, not the achieved precision.
3. **Do not implement MLE for σ.** It is a better estimator of an unidentifiable quantity.
4. **Report the three routes side by side** — CVE range, MSD intercept, CRLB — and treat disagreement
   as the honest error bar. On this data they span 3.6 to 5.2 nm.
5. **Acquire the immobile sample with `Tracking_2D_Oct2022`.** Fitting σ against 1/√ECO gives the
   prefactor directly; that collapses the whole question to one measured constant.

### What is still open

- **ν for this microscope.** Fitted as 195.5 nm from Vogler's table and transferred. A PSF measurement
  or the vendor's doughnut FWHM would replace the transfer.
- Whether the achieved precision sits at the CRLB or above it. Only the immobile sample answers this.
