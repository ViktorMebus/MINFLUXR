# Third-party code and attribution

MINFLUXR is released under the **GNU General Public License v3.0** (see `LICENSE`).
This file records the external work it builds on, what was taken, and what was changed.

---

## 1. MinSpt / Toolbox_for_MINFLUX_SPT — GPL-3.0

**Vogler et al. (2025)**, *Commun. Biol.*, `Toolbox_for_MINFLUX_SPT` (MinSpt), licensed **GPL-3.0**.

**Ported to R in MINFLUXR:**
- **Time Signature Analysis (TSA):** cycle-time recovery and the cumulative cycle-index timeline
  (`TOT`-style classification), `server_tracking.R`.
- **MSD workflow:** cycle-index lag binning, the fit window convention (`fit_limits`, default 5–50), and
  the **OLSF** optimal-point-count fit after Michalet & Berglund, including the exact iterative variant
  (`est_msd.R`: `msd_fit_olsf`, `msd_fit_olsf_iterative`, `.michalet_p`).

**Changes made:**
- Rewritten in R; the app infers the cycle time from the Δt distribution, where MinSpt simulates it from
  the declared sequence geometry and fits the overheads.
- The linearized MSD axis is `k·t_c + SW`, matching MinSpt's `T_LIN`.
- An "anchored" OLSF variant was added alongside the faithful iterative one, because the iterative loop is
  initialisation-dependent on short tracks.
- Motion blur: both variants take `R` (default `MSD_BLUR_R = 1/6.2`, the MinSpt value) and add the blur
  offset back when converting the fitted intercept to σ².

**Because this code is a port of GPL-3.0 software, MINFLUXR as a whole is GPL-3.0.**

Citation: Vogler et al. 2025, *Communications Biology*; see `docs/MINFLUX_analysis_methods.md`.

---

## 2. pyMINFLUX — Apache License 2.0

**pyMINFLUX** (Aaron Ponti, Single Cell Facility, D-BSSE, ETH Zurich),
<https://github.com/bsse-scf/pyMINFLUX>, licensed **Apache-2.0**.

Apache-2.0 code may be used inside a GPL-3.0 work (the combination is GPL-3.0), provided attribution and
notice of modification are kept. That is the purpose of this section.

**Reimplemented in R in MINFLUXR:**
- **DCR pooling across iterations** (`helpers_shared.R`, `.minflux_finalize_df` / `pool_dcr`): the
  ECO-weighted mean of `dcr1` over iterations with non-zero DCR and ECO, per trace ID.
- **Per-trace spread** (`helpers_smlm.R`, `sd_by_tid`): the per-axis standard deviation of the
  localizations within one trace ID, used as a quality statistic.

**Changes made:**
- Rewritten in R from the documented behaviour; no pyMINFLUX source is copied or bundled.
- The spread statistic is offered in two forms: the pyMINFLUX-matching `stats::sd` (default) **and** a
  robust MAD-based variant (`robust = TRUE`, scale constant 1.4826). pyMINFLUX uses the plain SD within a
  trace and applies MAD only across traces.
- **MINFLUXR adds a spread-based trace filter, which pyMINFLUX does not have**, and it removes whole
  traces rather than individual localizations. This deviation is documented in the app's UI and in
  `docs/MINFLUX_analysis_methods.md`.
- The z-scaling default of 0.7 for CSV exports follows pyMINFLUX's convention and is flagged in the app as
  *assumed* rather than recorded.

Citation: Ponti A., pyMINFLUX, ETH Zurich D-BSSE Single Cell Facility.

---

## 3. R packages

MINFLUXR imports CRAN packages (shiny, plotly, dplyr, dbscan, RANN, mclust, changepoint, alphashape3d,
geometry, fst and others; see `packages.R`). They are **dependencies, not bundled code**: they are
installed from CRAN at run time and are not redistributed with this repository. Their own licences apply
to them. Two are worth naming because they are more restrictive than GPL-3:

| Package | Licence | Note |
|---|---|---|
| `alphashape3d` | GPL-2 only | Used for 3D alpha shapes (visualisation + enclosed volume). Not bundled. |
| `fst` | AGPL-3 (or its own terms) | Used for on-disk chunking of large tracking files. Not bundled. |

If MINFLUXR is ever distributed as a bundle that *includes* these packages, revisit this table: a GPL-2-only
component cannot be combined with GPL-3 code in one distributed work.

---

## 4. Methods implemented from published descriptions

These are implementations from the literature, not ports of anyone's code. They carry no licence
obligation, only a citation one, and every one of them is cited in the app's method boxes and in
`docs/MINFLUX_analysis_methods.md`: Vestergaard et al. 2014 (CVE), Michalet & Berglund 2012 (OLSF, motion
blur), Schütz et al. 1997 (jump-distance distribution), Killick et al. 2012 (PELT), Das et al. 2009 and
Persson et al. 2013 (diffusive-state HMM), Jackson 2011 (Markov jump processes at irregular times),
Ripley 1977 / Illian et al. 2008 / Baddeley et al. 2015 (spatial statistics and edge correction),
Loosmore & Ford 2006 and Myllymäki et al. 2017 (global envelope tests), Ester et al. 1996 (DBSCAN),
Loftsgaarden & Quesenberry 1965 (k-NN density), Fisher et al. 1922 (index of dispersion),
Lotwick & Silverman 1982 and Dixon 2002 (bivariate point patterns), Costes et al. 2004 (randomisation
co-localization test), Simson et al. 1995 (confinement index).
