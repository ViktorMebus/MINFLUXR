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
