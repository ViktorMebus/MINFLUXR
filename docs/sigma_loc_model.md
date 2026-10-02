# Localization precision for MINFLUX tracking: the model, its derivation, and its calibration

This document derives the localization-error model the Tracking tab uses, states what each term is
for, explains why σ cannot be read off the trajectory, and describes how to calibrate the model on an
immobile sample.

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

### 2.3 The bound in 1D

```math
\sigma_\text{1D}(x) \;=\; I_F^{-1/2} \;=\; \frac{L}{4\sqrt{N}}\left[1 + 4\left(\frac{x}{L}\right)^2\right]
```

At the center, $\sigma = L/(4\sqrt{N})$. The 4 is $\sqrt{16}$, and the 16 is $(p')^2/[p(1-p)] \cdot L^2$
evaluated at $x = 0$, where $p' = -2/L$ and $p(1-p) = 1/4$. Tighter beams (smaller $L$) make $p$ swing
faster per nanometre of displacement, which is why MINFLUX beats camera SMLM at equal $N$.

**This is the two-beam, one-axis case. Tracking data is 2D or 3D, and the prefactor changes.**

### 2.4 The bound for a real pattern

A 2D pattern places $K$ beams on a circle of diameter $L$ (3 for a triangle, 4 for a square, 6 for a
hexagon). The same argument with a multinomial split instead of a binomial one:

1. Beam $k$ sits at $\mathbf{b}_k$ with $|\mathbf{b}_k| = L/2$, and $\sum_k \mathbf{b}_k = 0$ by symmetry.
   Its rate is $\lambda_k \propto |\mathbf{r} - \mathbf{b}_k|^2$, and $p_k = \lambda_k / \sum_j \lambda_j$.
2. At the center every beam is equally bright: $\lambda_k = L^2/4$, the total is $KL^2/4$, so $p_k = 1/K$.
3. Moving the emitter along $x$ changes $\lambda_k$ at rate $-2b_{k,x}$ and leaves the total unchanged
   (its derivative is $2Kx = 0$ at the center), so
   ```math
   \frac{\partial p_k}{\partial x} = \frac{-2\,b_{k,x}}{KL^2/4} = -\frac{8\,b_{k,x}}{KL^2}
   ```
4. Fisher information for a multinomial is $I_x = N \sum_k (\partial p_k/\partial x)^2 / p_k$:
   ```math
   I_x = N \sum_k \frac{64\, b_{k,x}^2}{K^2 L^4} \cdot K = \frac{64N}{KL^4} \sum_k b_{k,x}^2
   ```
5. For $K \ge 3$ evenly spaced beams, $\sum_k b_{k,x}^2 = K (L/2)^2 / 2 = KL^2/8$ (each beam contributes its
   squared radius times $\cos^2$, which averages to $1/2$). So
   ```math
   I_x = \frac{8N}{L^2}, \qquad \boxed{\;\sigma_\text{2D} = \frac{L}{\sqrt{8N}}\;}
   ```

That is Balzarotti et al. 2017 and Eq. 8 of Scheiderer, Marin & Ries 2024. It does not depend on $K$:
triangle, square and hexagon share $a = \sqrt{8} \approx 2.83$. The two-beam case gave 4 because both of
its beams lie on the measured axis; in 2D, each beam's displacement is split between $x$ and $y$, so
each photon carries less information about either one. An idealised octahedron (six beams on $\pm x$,
$\pm y$, $\pm z$) gives $\sum_k b_{k,x}^2 = L^2/2$ with $K = 6$, so $I_x = 16N/(3L^2)$ and
$a = \sqrt{16/3} \approx 2.31$.

**Off the center, the patterns differ.** Solving the full Fisher matrix numerically (off-center the
axes couple, so the per-axis variance is the diagonal of $I^{-1}$) and averaging over axes and over
the direction of the offset gives, for small offsets,

```math
\sigma^2(r) \;\approx\; \sigma_0^2\left(1 + c\,\frac{r^2}{L^2}\right)
```

| pattern | axes | ideal $a$ | off-center coefficient $c$ |
|---|---|---|---|
| triangle | 2 | $\sqrt{8} = 2.83$ | 13.7 |
| square | 2 | $\sqrt{8} = 2.83$ | 9.0 |
| hexagon | 2 | $\sqrt{8} = 2.83$ | 8.6 |
| octahedron (isotropic idealisation) | 3 | $\sqrt{16/3} = 2.31$ | 9.0 |

The triangle loses precision fastest off-center because it samples the fewest directions. Real 3D
MINFLUX uses an axial beam whose profile differs from the lateral one, so the octahedron row is a
guide rather than a prediction.

During tracking the emitter drifts off the center while photons arrive. Its squared offset grows as
$2\,d\,D\,s$ over time $s$ ($d$ = axes the pattern spans), and photons arrive throughout the exposure,
so the exposure-averaged value is $\langle r^2 \rangle = d\,D\,t_\text{loc}$. At $D = 0.2\ \mu\text{m}^2/\text{s}$,
$t_\text{loc} = 0.2$ ms and $L = 100.8$ nm this raises σ by about 3.5% for the hexagon and 5.5% for the
triangle. The table is recomputed from the Fisher matrix in `tests/testthat/test-estimators-tracking.R`.

### 2.5 Why the prefactor has to be measured

The ideal bound leaves out effects that make σ worse:

| effect | direction |
|---|---|
| finite doughnut depth / beam size ν | σ worse, and **faster than ∝ L** as $L \to \nu$ |
| background photons | σ worse: they carry no position information and flatten $p$ |
| detector and electronics | σ worse |

and one that can make σ look *better*: **the instrument's estimator is biased.** The Cramér-Rao bound
applies to unbiased estimators only. The Abberior estimator corrects a least-mean-squares estimate with
pre-computed coefficients and reports positions that are pulled toward the TCP center, and the
tracking loop centers each pattern on the previous estimate, which correlates successive errors. Both
let the exported positions scatter less than an unbiased estimate would.

Immobile-sample calibrations bear this out: fixed cells tracked with a 2D hexagonal sequence at
$L = 100.8$ nm give $a$ of about 3.5 to 4.2, above the ideal $\sqrt{8}$. So **no Cramér-Rao value bounds $a$ in either direction for exported data**,
and $a$ is an instrument constant to be measured, per sequence, rather than a derived number:

```math
\boxed{\;\sigma_i \;=\; \frac{L_\text{TCP}}{a\sqrt{N_i}}\,\sqrt{1 + c\,\frac{\langle r^2\rangle}{L^2}}\,,
\qquad N_i = \text{ECO}_i \;}
```

What survives from the theory is the functional form: $\sigma^2 \propto L^2/N$, which §4 tests directly,
and the pattern-specific coefficient $c$.

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

### 3.2 What an immobile sample shows

On an immobile sample all apparent motion is instrumental, so the MSD would be flat at $4\sigma^2$ if
the only error were white localization noise. On a feedback-tracked instrument it is not: it rises
sub-diffusively at short lags, flattens into a shoulder and then creeps upward. It is neither a plateau
(white noise alone) nor a straight line (white noise plus free diffusion).

**This is a feedback signature.** The TCP center is updated from the previous, noisy estimate with a
gain $g < 1$, so the error obeys $\varepsilon_i \approx (1-g)\varepsilon_{i-1} + g\,\eta_i$, an AR(1)
process that relaxes on a timescale $t_c/g$ and looks sub-diffusive over exactly the lags that can be
observed. The residual motion is the localization noise itself, filtered through the tracking loop.

Two consequences:

1. **The immobile MSD has power at every reachable lag**, so there is no lag short enough for the
   noise to stand alone. The wander and the noise are one physical quantity, not two independent
   processes to be separated.
2. **Fitting `MSD(k) = 4σ² + 4D_w·k·t_c` over short lags overestimates σ.** The curve is concave, so a
   straight-line fit extrapolated back to zero returns an intercept above the true value.

### 3.3 The internal-replicate channel, and why it fails too

Each trace opens with an iteration ramp: the same molecule localized at successively smaller $L$ within
a few milliseconds, with **disjoint photons**. If those errors were independent, a three-cornered-hat
solve would give σ at each L with no motion model:

```math
\sigma_2^2 = \tfrac{1}{2}\left(V_{12} + V_{23} - V_{13}\right), \qquad V_{jk} = \text{Var}(x_j - x_k)
```

and any common-mode drift would cancel, because $\Delta t_{13} = \Delta t_{12} + \Delta t_{23}$. The
premise fails for the same reason as §3.2: each iteration is centered on the previous one's estimate,
so consecutive iteration errors are positively correlated, and the solve can return a negative
variance. The channel cannot isolate σ.

---

## 4. Testing the model on an immobile sample

Since σ cannot be isolated, test the *model* instead. Bin segments by their own median ECO, fit the
short-lag MSD within each bin, and regress the intercept on $\langle 1/N\rangle$. If
$\sigma^2 = L^2/(a^2 N)$ the relation is linear **through the origin**; a photon-independent error
floor (TCP placement, background, stage) would appear as a positive intercept.

**What the fit returns is an effective prefactor.** Two biases act on the short-lag intercept in
opposite directions. The concave wander of §3.2 pushes a straight-line intercept up, and the positive
correlation that the tracking loop puts between successive errors pulls the lag-1 MSD down. Neither
can be bounded from the data alone, so the fitted $a$ is neither an upper nor a lower bound. It is the
value that reproduces the noise as it enters short-lag displacements, which is the quantity the
downstream estimators subtract, and that is the sense in which it is the right calibration. Record it
per sequence: pattern, $L$ and photon limit all change it.

**Calibrate on fixed cells labelled as the live sample, not on beads.** A bead is an extended, bright,
coverslip-mounted object: it is not a point emitter, it sits at a different axial position relative to
the doughnut, and its brightness drives the tracker harder, so it brackets $a$ only loosely.

**The beam-size term needs more than two diameters.** A prefactor that falls as $L$ grows suggests the
beam-size denominator, but two $L$ values cannot distinguish it from an L-dependent background
difference, so it is off by default. **Calibrate at the $L$ you track at**, and the L-scaling question
does not arise.

---

## 5. The wander floor $D_w$

The same fits give the wander as the short-lag slope. Because it is the localization noise fed back
through the loop, it scales the way $\sigma^2$ does: it falls with photon count and worsens steeply
with $L$. It is therefore **not a constant of the instrument** and must never be subtracted silently:
report $D - D_w$ beside $D$, and flag when the floor is a large fraction of the measurement.

$D_w$ from a straight-line short-lag fit is a *rate over those lags*, not a diffusion coefficient: the
underlying process relaxes rather than diffuses. Treat it as the apparent diffusivity of a stationary
target over the lags you fit your data on, which is the comparison that matters, and quote the lag
range with it.

---

## 6. Which recorded quantities are needed

| field | verdict | reason |
|---|---|---|
| **ECO** | **required** | $N$ in the CRLB. The only per-localization noise driver. |
| **L_TCP** | **required** | from the sequence's `patGeoFactor`; sets the scale. |
| `tim` | required | cycle time, gap structure, and the lag axis. |
| `itr` | required | select the final (reiterated) iteration. |
| `vld` | required | drop invalid records. |
| **EFO** | **not needed for σ** | $= \text{ECO}/t$; adds nothing the pair (ECO, $t_c$) lacks. Useful for QC (fiducial detection). |
| **cfr** | **not needed for σ** | it also reflects how well the emitter is centered, so it is not a clean background measure. |
| `dcr` | not needed for σ | axial/color discrimination. |

**ECO and L_TCP are sufficient:** one photon count, one pattern diameter, one prefactor.

---

## 7. The photon-budget model for MINFLUX

**The form is right, with two qualifications.**

The photon-budget CRLB is the correct channel, and the only one available: §3 shows the trajectory
cannot supply σ independently on a feedback-tracked instrument, whatever the estimator. Eq. 1 of Vogler
et al. 2025 is that CRLB with the two refinements that matter for tracking: the emitter-offset term,
since the molecule is never exactly at the TCP center, and the beam-size denominator.

1. **The beam-size term is opt-in.** It requires ν, an instrument constant most users cannot state,
   and its effect is confined to $L \to \nu$. Defaulting it on would import an unverifiable number into
   every σ.
2. **The prefactor is not universal, and theory does not fix it.** The ideal unbiased value depends
   only on the pattern ($\sqrt{8}$ in 2D), but exported positions come from a biased estimator inside a
   feedback loop, and an immobile-sample calibration lands above that value (§2.5). Exposing $a$ as a
   measured instrument constant is what makes the model portable. A calibration far *below*
   $\sqrt{8}$ is the informative warning sign: background or alignment then dominates, and $\sigma$ grows
   faster with falling $N$ than the photon model alone predicts.

Eq. 1's offset term is written for one axis. The app replaces it with the coefficient of the actual
pattern (§2.4), and keeps Eq. 1's beam factor evaluated at the center.

---

## 8. How the app uses it

```math
\sigma_i = \frac{L_\text{TCP}}{a\sqrt{N_i}}\,B\,\sqrt{1 + c\,\frac{d\,D\,t_\text{loc}}{L^2}}
```

per localization, where $B = 1/|1 - \ln 2\,(L/\nu)^2|$ is the beam factor (1 when ν is off), $c$ and $d$
come from the **TCP pattern** setting, $D$ from the CVE and $t_\text{loc}$ is the single-cycle interval.

- **$a$ comes from a calibration profile**, measured on an immobile sample with the sequence in
  question (*Instrument Calibration*, `R/est_calib.R`) and kept with the sequence it belongs to.
  Without one, $a$ is the ideal-estimator value of the pattern ($\sqrt{8}$ in 2D) and σ is labelled
  uncalibrated. There is no system-independent value to ship: the effective $a$ depends on the
  instrument's estimator, background, photon limit, dwell, L and damping.
- **The TCP pattern is read from the sequence file** (final iteration) when one is loaded, and can be
  set by hand. Unrecognised patterns fall back to the hexagon.
- **Calibrate on fixed cells labelled as the live sample**, with the sequence used on the live sample.
- **Report $D_w$ alongside D**, never subtract it, and quote the lag range it was fitted over.

### What to report in a methods section

> Localization precision was computed per localization from the MINFLUX photon budget,
> σ = L/(a√N)·√(1 + c·⟨r²⟩/L²), with N the photon count at offset (ECO), L the TCP diameter taken
> from the acquisition sequence (patGeoFactor × 360 nm = XXX nm), c = X.X the off-center coefficient
> of the [hexagonal] pattern and ⟨r²⟩ the exposure-averaged offset of the moving emitter. The prefactor
> a = X.XX was measured as the short-lag MSD intercept of a fixed, identically labelled sample tracked
> with the same sequence. The apparent diffusion coefficient of that immobile sample was
> D_w = XXX nm²/ms over lags 1–5, quoted alongside all reported D.
