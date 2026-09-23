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

**The prefactor is bracketed, not measured.** Because §3.2 shows the MSD intercept is biased *upward*,
the fitted $a$ is a **lower bound on $a$**, equivalently an **upper bound on σ**. Theory gives the other
side, $a \le 4$. The result is a window on σ, which is a stronger statement than a point estimate.

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
et al. 2025 is that CRLB with the two refinements that matter: the emitter-offset term, essential for
tracking, where the molecule is never at the TCP center, and the beam-size denominator.

1. **The beam-size term is opt-in.** It requires ν, an instrument constant most users cannot state,
   and its effect is confined to $L \to \nu$. Defaulting it on would import an unverifiable number into
   every σ.
2. **The prefactor is not universal.** A fixed prefactor implicitly fixes the collection efficiency,
   background and estimator of one microscope. Exposing $a$, with the ideal 4 as its ceiling, makes the
   model portable and gives a self-check: a calibration returning $a > 4$ means something is wrong
   (usually σ underestimated by an over-cleaned dataset), and $a \ll 4$ means background or alignment
   dominates, so the photon model alone will understate σ.

---

## 8. How the app uses it

$\sigma_i = \dfrac{L_\text{TCP}}{a\sqrt{N_i}}\left[1 + c\left(x/L\right)^2\right]$ per localization,
with $a$ defaulting to 4. The default is the theoretical bound rather than a fitted number because:

- a fitted value is biased in a known direction (§3.2), so it is a bound, not an estimate;
- $a = 4$ makes σ a stated *lower bound* on the localization error, the conservative direction for
  every downstream test that subtracts σ (JDD, confinement, D-map, HMM);
- it is portable with no calibration, and any user can improve on it.

In practice:

1. **Report the bracket**: σ at $a = 4$, and the immobile-sample bound when one has been measured.
2. **$a$ is user-adjustable, with 4 as its upper limit.**
3. **Calibrate on fixed cells, not beads**, at the $L$ you track at, and treat the result as an upper
   bound on σ.
4. **Report $D_w$ alongside D**, never subtract it, and quote the lag range it was fitted over.

**Limitation: the offset coefficient $c$.** The app uses $c = 2$. The exact value for a 1D two-beam
pattern is 4, and the value for the 2D hexagonal TCP has not been derived. The term matters only when
the emitter moves a sizeable fraction of $L$ within one localization, i.e. for fast diffusion.

### What to report in a methods section

> Localization precision was computed per localization from the MINFLUX Cramér-Rao bound,
> σ = L/(a√N), with N the photon count at offset (ECO), L the TCP diameter taken from the acquisition
> sequence (patGeoFactor × 360 nm = XXX nm) and a = 4 the ideal doughnut prefactor. A fixed,
> identically labelled sample gave an upper bound of a ≥ X.XX from the short-lag MSD intercept,
> bounding σ within X% of the ideal. The apparent diffusion coefficient of that immobile sample was
> D_w = XXX nm²/ms over lags 1–5, quoted alongside all reported D.
