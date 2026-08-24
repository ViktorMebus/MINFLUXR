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
⟨Δx²⟩ = 92.95 nm²      cov₁ = −2.38 nm²
R = 0    (app today)   σ² =        2.38   → σ = 1.54 nm
R = 1/6  (correct)     σ² = 15.49 + 1.59  → σ = 4.13 nm
```

Three separate problems, in increasing order of severity:

1. **The app assumes the wrong R.** A MINFLUX localization integrates photons across the cycle, so
   R = τ/(6·Δt) ≈ 0.158, not 0. Using R = 0 under-reports σ_loc by ~3× on both GP67H files
   (1.5 and 1.3 nm against 4.13 and 4.08 nm). This is blur, not precision.
2. **The estimate is dominated by the assumption, not by the data.** At R = 1/6, 91% of σ² comes
   from the R·⟨Δx²⟩ term. σ̂ ≈ √R · rms_step is therefore close to a deterministic function of the
   step size and the assumed R. The actual noise signature — the part cov₁ carries — is a 4%
   residual (r₁ = −0.026). Whatever R is chosen, this number is not really measuring σ_loc.
3. **It inherits every motion-model error.** ⟨Δx²⟩ is inflated by population heterogeneity, by drift,
   and by any confinement structure, and all of that flows straight into σ̂.

The photon route has none of these properties: **σ_i = L/(4·√N_i)** with L = 100.8 nm from the
sequence `patGeoFactor` and N = `eco` gives **4.09 / 4.20 nm** median on the two files — agreeing with
the blur-aware CVE to within 3%, per-localization rather than pooled, and with no motion model at all.

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
| *real data* | 100 75 53 42 36 30 | |

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

### Phase 1 — make σ_loc photon-based (the actual fix)

- Add a **photon-CRLB σ_loc** as the app's primary estimate: σ_i = L/(4√N_i) per localization, with
  L pre-filled from the uploaded sequence .json (`trk_seq_params()$L_tcp_nm` already parses it) and
  N from `eco`. Falls back to `DMAX_L_TCP_DEFAULT` when no sequence is loaded.
- Re-point `trk_sigma_loc()`'s priority chain at it. Today the chain prefers the CVE value
  (`trk_sigma_loc_cve`), which is the number §2.2 says not to trust.
- Because it is per-localization, it can feed the heteroscedastic consumers directly rather than
  through a pooled median.

### Phase 2 — fix, demote and label the CVE σ

- Add the R term: σ² = R⟨Δx²⟩ + (2R−1)cov₁ with R = τ/(6Δt), τ from the sequence
  (`patDwellTime` × pattern positions) and Δt = t_c + SW. Default R = 1/6 when no sequence is loaded.
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

## 6. Experiment worth requesting

None of the above measures σ_loc **directly**. The clean measurement is an **immobile sample acquired
with the tracking sequence** — fixed cells or fiducials under `Tracking_2D_Oct2022`, where the
within-trace spread *is* σ_loc with no motion model in the way. Fitting σ against `eco` over that
acquisition gives the constant in σ = a/√N empirically, which is the one number the entire photon
route rests on. The existing fiducial pair is imaged with an SMLM sequence (L = 39.6 nm), so it
calibrates the imaging iteration, not the tracking one.
