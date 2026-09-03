# Localization precision for MINFLUX tracking: the model, its derivation, and its calibration

This document derives the localization-error model the Tracking tab uses, states what each term is
for, and shows what the immobile-sample data does and does not support. It supersedes the
options survey in `sigma_loc_estimation.md`.

---

## 1. What σ_loc is, and what it is not

For a trace, write each reported position as

```math
x_i = x_i^\ast + \varepsilon_i
```

with $x^\ast$ the true emitter position and $\varepsilon$ the estimation error.
**σ_loc is the per-axis standard deviation of ε.** Three things are commonly folded into it and
should not be:

| | what it is | where it belongs |
|---|---|---|
| **σ_loc** | shot-noise-limited error of one position estimate | this document |
| **motion blur** | the estimate is a time-average of the path over the cycle, so it is a *shrunken* position, not a noisier one | the D estimator (CVE, `R = t_c/6Δt`) |
| **tracker wander** | the feedback loop repositions the TCP using its own noisy estimate, so the reported coordinate drifts even for a stationary emitter | §5, reported separately as $D_w$ |

Conflating the third with the first is the specific failure mode of every trajectory-based
estimator, and §5 shows it measurably.

---

## 2. Derivation of the MINFLUX bound

### 2.1 The measurement

An ideal doughnut has an exact zero at its center, so near the minimum the intensity is quadratic,
$I(r) \propto r^2$. Place two beams at $b_\pm = \pm L/2$ and put the emitter at $x$. The expected
detection rates are

```math
\lambda_\pm \;\propto\; \left(x \mp \tfrac{L}{2}\right)^2
```

Collect $N$ photons in total. Which beam a photon arrives under is a Bernoulli draw, so the count in
the $+$ beam is $n_+ \sim \text{Binomial}(N, p)$ with

```math
p(x) \;=\; \frac{\left(x - \tfrac{L}{2}\right)^2}{\left(x - \tfrac{L}{2}\right)^2 + \left(x + \tfrac{L}{2}\right)^2}
\;=\; \frac{x^2 - xL + L^2/4}{2x^2 + L^2/2}
```

**All the position information is in this ratio.** That is the whole of MINFLUX: the emitter's
position is read off how the photons split between beam positions, and $N$ enters only through the
counting statistics of that split.

### 2.2 Fisher information

For a binomial with $N$ trials and success probability $p(x)$,

```math
I_F(x) \;=\; N\,\frac{\left(p'\right)^2}{p\,(1-p)}
```

Write $u = (x - L/2)^2$, $v = (x + L/2)^2$, $D = u + v = 2x^2 + L^2/2$. Then

```math
p' = \frac{u'v - uv'}{D^2}
   = \frac{(2x-L)\left(2x^2 + \tfrac{L^2}{2}\right) - 4x\left(x^2 - xL + \tfrac{L^2}{4}\right)}{D^2}
   = \frac{L\left(2x^2 - \tfrac{L^2}{2}\right)}{D^2}
   = \frac{L\left(4x^2 - L^2\right)}{2D^2}
```

and, using $\left(x^2 + \tfrac{L^2}{4}\right)^2 - x^2L^2 = \left(x - \tfrac{L}{2}\right)^2\left(x + \tfrac{L}{2}\right)^2$,

```math
p(1-p) = \frac{\left(x-\tfrac{L}{2}\right)^2\left(x+\tfrac{L}{2}\right)^2}{D^2}
```

Since $\left(4x^2 - L^2\right)^2 = 16\left(x-\tfrac{L}{2}\right)^2\left(x+\tfrac{L}{2}\right)^2$, the
emitter-position factors cancel exactly:

```math
I_F(x) = N \cdot \frac{L^2\left(4x^2-L^2\right)^2}{4D^4} \cdot \frac{D^2}{\left(x-\tfrac{L}{2}\right)^2\left(x+\tfrac{L}{2}\right)^2}
       = \frac{4NL^2}{D^2}
       = \frac{16N}{L^2\left(1 + 4x^2/L^2\right)^2}
```

### 2.3 The bound

```math
\boxed{\;\sigma_\text{CRLB}(x) \;=\; I_F^{-1/2} \;=\; \frac{L}{4\sqrt{N}}\left[1 + 4\left(\frac{x}{L}\right)^2\right]\;}
```

Two things to read off this:

- **At the center, $\sigma = L/(4\sqrt{N})$.** The 4 is not a fitted number. It is
  $\sqrt{16}$, and the 16 is $\left(p'\right)^2/[p(1-p)] \cdot L^2 = (2/L)^2 \cdot 4 \cdot L^2$ evaluated
  at $x = 0$, where $p' = -2/L$ and $p(1-p) = 1/4$. Tighter beams (smaller $L$) make $p$ swing faster
  per nanometre of displacement, which is the entire reason MINFLUX beats camera SMLM at equal $N$.
- **Precision degrades away from the center**, quadratically. During tracking the emitter is
  displaced from the TCP center by whatever it moved since the last update, $x_\text{rms} \approx \sqrt{4Dt_\text{loc}}$.

> The app currently uses $\left[1 + 2(x/L)^2\right]$. The exact 1D two-beam result derived above is
> $\left[1 + 4(x/L)^2\right]$. At $D = 0.2\ \mu\text{m}^2/\text{s}$, $t_c = 0.2$ ms and $L = 100.8$ nm
> this is the difference between a 3.1% and a 6.2% correction, so it does not change any conclusion,
> but the coefficient should be checked against the 2D hexagonal-TCP geometry rather than assumed.

### 2.4 What the ideal 4 does not include

| effect | direction | why it is not derivable |
|---|---|---|
| finite doughnut depth / beam size ν | σ worse, and **faster than ∝ L** as $L \to \nu$ | needs the real PSF, not the parabola |
| background photons | σ worse | carries no position information, flattens $p(x)$ |
| 2D/3D TCP geometry (hexagon + center vs. two beams) | modest factor | different beam layout, different $p$ |
| estimator inefficiency (mLMS is not maximum-likelihood) | σ worse | estimator-specific |
| detector and electronics | σ worse | instrument-specific |

Every one of these multiplies σ by a constant **at fixed L**. So they collapse into one number:

```math
\boxed{\;\sigma_i \;=\; \frac{L_\text{TCP}}{a\sqrt{N_i}}\left[1 + c\left(\frac{x}{L}\right)^2\right],
\qquad N_i = \text{ECO}_i \;}
```

with $a \le 4$, and $a = 4$ the unreachable ideal. **The model has exactly one free parameter**, and
it has a hard theoretical ceiling — which is what makes it falsifiable.

---

## 3. Why σ cannot be read off the trajectory

This is the part worth being explicit about, because it is the standard SPT approach and it does not
work here.

### 3.1 The identifiability argument

CVE, MSD-intercept and MLE all separate σ from motion using the same fact: localization noise is
white, motion is correlated, so the lag-1 covariance of the increments is negative and equal to
$-\sigma^2$. That identity holds **only** for instantaneous exposure. With a finite exposure the
general result is

```math
\sigma^2 = R\,\langle \Delta x^2\rangle + (2R-1)\,\text{cov}_1, \qquad R = \frac{t_c}{6\Delta t}
```

At $R \approx 1/6$ roughly 90% of the answer is the first term, which is just
$\sqrt{R} \times$ rms step. Blur and precision enter the covariance in the same place, so no estimator
built on the increment covariance can separate them from a single trajectory. This is an
identifiability limit, not a choice of estimator. MLE inherits it.

### 3.2 The measurement that confirms it

On an immobile sample all apparent motion is instrumental, so the MSD should be flat at $4\sigma^2$
if the only error were white localization noise. Pooled over the fixed-cell datasets at $L = 100.8$ nm:

| lag k | τ (ms) | MSD (nm²) | local log-log slope |
|---|---|---|---|
| 1 | 0.20 | 140 | – |
| 2 | 0.41 | 247 | 0.80 |
| 4 | 0.82 | 425 | 0.80 |
| 7 | 1.43 | 596 | 0.40 |
| 15 | 3.07 | 732 | 0.21 |
| 45 | 9.22 | 1019 | 0.40 |
| 174 | 35.7 | 1623 | 0.32 |

Three regimes: a steep sub-diffusive rise at the shortest lags, a shoulder near 700–800 nm² around
2–4 ms, then a slow secondary creep. Neither a plateau (pure white noise) nor a straight line (white
noise plus free diffusion).

**This is a feedback signature.** The TCP center is updated from the previous, noisy estimate with a
gain $g < 1$, so the error obeys $\varepsilon_i \approx (1-g)\varepsilon_{i-1} + g\,\eta_i$ — an AR(1)
process, which relaxes on a timescale $t_c/g$ and looks sub-diffusive over exactly the lags we can
observe. The residual motion is not sample motion and not free diffusion; it is the localization
noise itself, filtered through the tracking loop.

Two consequences:

1. **The immobile MSD has power at every lag we can reach**, so there is no lag short enough for the
   noise to stand alone. The wander and the noise are the same physical quantity, not two independent
   processes to be separated.
2. **Fitting `MSD(k) = 4σ² + 4D_w·k·t_c` over short lags overestimates σ.** The curve is concave,
   so a straight-line fit extrapolated back to zero returns an intercept above the true value.

### 3.3 The internal-replicate channel, and why it fails too

Each trace opens with an iteration ramp: the same molecule localized at $L$ = 302.4, 151.2 and 100.8 nm
within ~1.7 ms, with **disjoint photons**. If those three errors were independent, a three-cornered-hat
solve would give σ at each L with no motion model at all:

```math
\sigma_2^2 = \tfrac{1}{2}\left(V_{12} + V_{23} - V_{13}\right), \qquad V_{jk} = \text{Var}(x_j - x_k)
```

and any common-mode drift cancels identically, because $\Delta t_{13} = \Delta t_{12} + \Delta t_{23}$.

Run on 4067 fixed-cell triples: $V_{12} = 16.56^2$, $V_{23} = 8.54^2$, $V_{13} = 21.02^2$, giving
$\sigma_2^2 = -47$ nm². **Negative**, so the premise is false: consecutive iteration errors are
positively correlated, again through the feedback loop that centers each iteration on the previous
one's estimate. The replicates are not independent, and the channel is unusable.

---

## 4. What the immobile data does establish

Since σ cannot be isolated, test the *model* instead. Bin segments by their own median ECO, fit the
short-lag MSD within each bin, and regress the intercept on $\langle 1/N\rangle$. If
$\sigma^2 = L^2/(a^2 N)$ the relation is linear **through the origin**; a photon-independent error
floor (TCP placement, background, stage) would appear as a positive intercept.

| sample | L (nm) | bins | N range | fit | R² | implied a |
|---|---|---|---|---|---|---|
| fixed cells | 100.8 | 39 | 24–41 | $\sigma^2 = 2247/N - 46$ | 0.66 | 3.80 |
| fixed cells | 151.2 | 25 | 25–48 | $\sigma^2 = 3767/N - 51$ | 0.43 | 3.50 |
| beads | 100.8 | 5 | 23–27 | $\sigma^2 = 2260/N - 25$ | 0.08 | 2.52 |
| beads | 151.2 | 12 | **25–144** | $\sigma^2 = 7489/N - 19$ | **0.94** | 1.84 |

**1. The functional form is right.** Every intercept is negative, i.e. consistent with zero: there is
no photon-independent error floor. The strongest single line is the bead L151 set, where the photon
count spans 6× and $\sigma^2 \propto 1/N$ holds at R² = 0.94. The model is not fitted to survive; it
had a free way to fail and did not.

**2. The prefactor is bracketed, not measured.** Because §3.2 shows the MSD intercept is biased
*upward*, the fitted $a$ is a **lower bound on $a$** — equivalently an **upper bound on σ**. Theory
gives the other side, $a \le 4$. So for the fixed cells:

```math
3.80 \le a \le 4.00 \quad (L = 100.8\ \text{nm}) \qquad
3.50 \le a \le 4.00 \quad (L = 151.2\ \text{nm})
```

At $L = 100.8$ that is a **5% window on σ**. This is a much stronger statement than a point estimate:
whatever the true prefactor is, σ on this instrument is within 5% of the theoretical MINFLUX bound.

**3. Beads are the wrong calibration standard.** They bracket only $1.84 \le a \le 4$. A 40 nm bead is
an extended, bright, coverslip-mounted object: it is not a point emitter, it sits at a different axial
position relative to the doughnut, and its brightness drives the tracker harder. Use fixed cells
labelled exactly as the live sample.

**4. There is mild evidence for the beam-size term.** $a$ falls from 3.80 to 3.50 between $L$ = 100.8
and 151.2, i.e. σ rises ~8.6% faster than $\propto L$. Solving the Vogler denominator
$\left|1 + \ln 2\,(L/\nu)^2\left((2x/L)^2 - 1\right)\right|$ for that excess gives $\nu \approx 255$ nm,
a plausible doughnut FWHM. But two L values cannot distinguish a beam term from an L-dependent
background difference, so this stays **off by default**. The practical escape is that you do not need
the L-scaling to be right — **calibrate at the L you actually track at**, and the question is moot.

---

## 5. The wander floor $D_w$

The same fits give the wander as the short-lag slope. It is large and it scales with $L$:

| sample | L (nm) | $D_w$ (nm²/ms) |
|---|---|---|
| fixed cells | 100.8 | 70 – 380 |
| fixed cells | 151.2 | 297 – 1132 |
| beads | 100.8 | 22 – 141 |

Because it is the localization noise fed back through the loop, it scales the way $\sigma^2$ does:
it falls with photon count and worsens with $L$. It is therefore **not a constant of the instrument**
and must never be subtracted silently — report $D - D_w$ beside $D$ and flag when the floor is a
large fraction of the measurement.

Note also that $D_w$ from a straight-line short-lag fit is a *rate over those lags*, not a diffusion
coefficient: the underlying process relaxes rather than diffuses. Treat it as "apparent diffusivity
of a stationary target over the same lags you fit your data on", which is exactly the comparison that
matters, and quote the lag range with it.

---

## 6. Which recorded quantities are needed

| field | verdict | reason |
|---|---|---|
| **ECO** | **required** | $N$ in the CRLB. The only per-localization noise driver. |
| **L_TCP** | **required** | from the sequence's `patGeoFactor`; sets the scale. |
| `tim` | required | cycle time, gap structure, and the lag axis. |
| `itr` | required | select the final (reiterated) iteration. |
| `vld` | required | drop invalid records. |
| **EFO** | **not needed for σ** | $= \text{ECO}/t$; adds no information the pair (ECO, $t_c$) lacks. Useful for QC (fiducial detection), not for precision. |
| **cfr** | **not needed** | tested below. |
| `dcr` | not needed for σ | axial/color discrimination. |

**The cfr test.** `cfr` is the center-to-outer frequency ratio, the natural background/SBR indicator,
so if background mattered beyond the photon count, high-`cfr` localizations should be noisier at
matched ECO. Over 1.70 million single-cycle steps from the fixed cells, MSD(1) by cfr quartile within
ECO tertiles:

| ECO | cfr Q1 | Q2 | Q3 | Q4 | ratio Q4/Q1 |
|---|---|---|---|---|---|
| ~27 | 181 | 180 | 173 | 168 | 0.93 |
| ~39 | 116 | 110 | 105 | 102 | 0.88 |
| ~59 | 60 | 47 | 40 | 43 | 0.71 |

The effect is 7–12% in MSD, so ~4–6% in σ, and it runs **the wrong way** for a background
explanation: higher cfr gives *tighter* steps, because cfr also tracks how well-centered the emitter
is. There is no case for adding it. (The third row is confounded — ECO itself varies 59–79 inside
that tertile.)

**Conclusion: ECO and L_TCP are sufficient.** One photon count, one pattern diameter, one prefactor.

---

## 7. Is Vogler's approach the right one for MINFLUX?

**Yes for the form, with two qualifications.**

The photon-budget CRLB is the correct channel, and it is the only one available: §3 shows the
trajectory cannot supply σ independently on a feedback-tracked instrument, no matter which estimator
is used. Vogler's Eq. 1 is that CRLB with the two refinements that matter — the emitter-offset term
(essential for tracking, where the molecule is never at the TCP center) and the beam-size denominator.

The qualifications:

1. **The beam-size term should be opt-in, not default.** It requires ν, an instrument constant most
   users cannot state, and its effect is confined to $L \to \nu$. Our own data supports it only
   weakly. Defaulting it on imports an unverifiable number into everyone's σ.
2. **The prefactor is not universal.** Vogler's parameterization implicitly fixes the collection
   efficiency, background and estimator of one microscope. Exposing $a$ with the ideal 4 as its
   ceiling is what makes the model portable, and gives every user a self-check: if their calibration
   returns $a > 4$, something is wrong (usually σ underestimated by an over-cleaned dataset); if it
   returns $a \ll 4$, background or alignment dominates and the photon model alone will understate σ.

---

## 8. Recommendation

**Use $\sigma_i = \dfrac{L_\text{TCP}}{a\sqrt{N_i}}\left[1 + c\left(x/L\right)^2\right]$, per localization,
with $a$ defaulting to 4.**

Defaulting to the theoretical bound rather than to a fitted number is the more honest choice, because:

- the fitted value is provably biased in a known direction (§3.2), so it is a bound, not an estimate;
- $a = 4$ makes σ a stated *lower bound* on the localization error, which is the conservative direction
  for every downstream test that subtracts σ (JDD, confinement, D-map, HMM) — those tests are then
  anti-conservative in D and conservative in "is there structure", which is the right way round;
- it is portable with no calibration at all, and any user can improve on it.

Then:

1. **Report the bracket, not a point value.** The calibration panel should show
   "σ = 4.85 nm at N = 27 (a = 4, ideal); immobile-sample upper bound 5.11 nm (a = 3.80)".
2. **Keep $a$ user-adjustable, with 4 as a hard upper limit** and a warning if a calibration returns
   more.
3. **Calibrate on fixed cells, not beads**, at the L you will track at, and treat the result as an
   upper bound on σ.
4. **Report $D_w$ alongside D**, never subtract it, and quote the lag range it was fitted over.
5. **Fix the offset coefficient** $c$: the exact 1D two-beam value is 4, the app uses 2, and the
   correct value for the 2D hexagonal TCP should be derived rather than inherited.

### What to report in a methods section

> Localization precision was computed per localization from the MINFLUX Cramér-Rao bound,
> σ = L/(a√N), with N the photon count at offset (ECO), L the TCP diameter taken from the acquisition
> sequence (patGeoFactor × 360 nm = XXX nm) and a = 4 the ideal doughnut prefactor. A fixed,
> identically labelled sample gave an upper bound of a ≥ X.XX from the short-lag MSD intercept,
> bounding σ within X% of the ideal. The apparent diffusion coefficient of that immobile sample was
> D_w = XXX nm²/ms over lags 1–5, quoted alongside all reported D.

---

## Appendix: reproducing the numbers

The analyses in §3–§6 were one-off scripts against local `Test data/`, which is gitignored, so they
are not part of the repo — they could not run in a fresh clone. What they did, in enough detail to
rewrite:

Common preparation: read the `.npy`, keep `vld == TRUE` and the final iteration, split traces at gaps
above 7 cycles, trim leading/trailing dark cycles of 3+, keep segments with continuity ≥ 0.6 and
≥ 150 localizations. `L` follows from the iteration count (2 → 302.4 nm, 3 → 151.2, 4 → 100.8);
`t_c` is the modal single-cycle Δt.

| § | analysis |
|---|---|
| 3.2 | pooled MSD over log-spaced lags 1–300, all fixed-cell files at L = 100.8; report the local log-log slope |
| 3.3 | three-cornered hat on the opening iteration ramp: first record of itr 1/2/3 per trace, robust variance of each pairwise coordinate difference, then $\sigma_2^2 = (V_{12}+V_{23}-V_{13})/2$ |
| 4, 5 | bin segments into ECO quartiles, fit `MSD(k) = 4σ² + 4D_w k t_c` on lags 1–5 within each bin, regress the intercept on $\langle 1/Nangle$ weighted by pair count |
| 6 | single-cycle squared steps only; cross-tabulate mean MSD(1) by cfr quartile within ECO tertiles |

Datasets: `Test data/Fixed Tracking datasets/` (fixed E. coli, JF646, three strains),
`Test data/Fixed beads dataset/` (40 nm red beads on coverslip).
