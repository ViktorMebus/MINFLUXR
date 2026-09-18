# MINFLUXR

An interactive R/Shiny application for visualization and analysis of **Abberior MINFLUX** data.

---

## Overview

MINFLUXR provides three main analysis environments:

| Tab | Purpose |
|-----|---------|
| **SMLM** | Static super-resolution data: filtering, de-blinking, clustering, co-localization, spatial statistics |
| **Tracking** | Single-molecule tracking: trajectory visualization, MSD, diffusion analysis |
| **MINFLUX Simulator** | Generate synthetic MINFLUX datasets - labelled structures, whole acquisitions, and ground-truth diffusion tracks |

Each tab is split into a **Display** sub-tab (sidebar controls + plot) and an **Analysis** sub-tab
(per-panel plots and exportable tables).

The app is designed around Abberior MINFLUX `.npy` and `.json` exports, but also accepts `.csv` files
in pyMINFLUX format or MINFLUXR's own export. In principle any `.csv` with x/y coordinates works,
though this has not been tested.

---

## Installation

MINFLUXR requires **R ≥ 4.1**. All dependencies are installed automatically on first launch via `packages.R`.

```r
# From the app directory:
source("packages.R")   # install dependencies (only needed once)
```

Key dependencies include: `shiny`, `shinyWidgets`, `shinyjs`, `plotly`, `DT`, `dplyr`, `dbscan`, `RANN`, `mclust`, `magick`, `alphashape3d`, `data.table`, `fst`, `readr`, `changepoint`, `minpack.lm`.

To check that every package the code calls is actually declared:

```sh
Rscript tools/check_dependencies.R
```

It parses the sourced files and reports anything referenced as `pkg::` or `library(pkg)` that
`packages.R` does not install, and exits non-zero if there is a gap.

`qs` is **optional**. It needs a C++ toolchain plus `RcppParallel` and `stringfish` and often fails to
install; the app therefore writes its tracking chunks with `fst`, which is a hard dependency and needs
no special toolchain. `qs` is only used to *read* a chunk set produced on a machine that had it, and
its absence is reported rather than thrown.

> **macOS note:** `fst` uses C++/OpenMP. If it fails to compile, install Xcode Command Line Tools
> (`xcode-select --install`) and retry.

> **If a package fails to install:** restart R before retrying. The installer skips anything already
> loaded in the session, because replacing a loaded package's files can leave it half-removed.

---

## Launching the app

```r
shiny::runApp("path/to/MINFLUXR")
```

Or open `app.R` in RStudio and click **Run App**. The maximum upload size is set to **2 GB** by default.

---

## A first run

1. **Launch** the app and pick a tab: **SMLM Analysis** for fixed-cell localization data,
   **Tracking Analysis** for single-particle tracks, **MINFLUX Simulator** to generate synthetic data.
2. **Load** an Abberior `.npy` or `.json` export (or a `.csv` from pyMINFLUX or from MINFLUXR itself).
   The Display sub-tab shows the localizations; the metadata filters, trace-length and trace-spread
   controls live in its sidebar.
3. **Press Plot** after changing anything in the clustering or filter blocks: the display and every
   analysis read the same filtered data, so they cannot disagree with each other.
4. **Open the Analysis sub-tab** for the quantitative panels. Each has a collapsed
   *Method & MINFLUX adaptation* box giving the origin of the method, what was adapted for MINFLUX and
   what it cannot tell you.
5. **Save your settings** (Download settings, SMLM tab) to reproduce a figure later: every input that
   changes a number is written to the JSON.

Tracking data: run **Data Quality & Calibration** first. It recovers the cycle time and the
localization precision that the MSD, JDD, confinement, HMM and D-map panels all read.

---

## How far to trust a panel

Every analysis panel carries a rating chip in its method box:

| Chip | Meaning |
|---|---|
| 🟢 well-adapted | Standard method, or a faithful port, correctly adapted to MINFLUX sampling |
| 🟡 adapted, with caution | Sound, but with a stated limitation you need to read before quoting it |
| 🟠 exploratory for MINFLUX | Research tool. Interpretation is not settled; do not quote it as a result |

The estimators behind the tracking panels live in `R/est_*.R` as plain functions, and
`tests/` checks them against ground truth from the simulator (`Rscript tests/run_tests.R`;
slow Monte-Carlo tests with `MINFLUXR_SLOW_TESTS=true`). The longer validation runs and their
results are in `tests/validation/`.

---

## Documentation

| File | Contents |
|------|----------|
| [`docs/user_guide.md`](docs/user_guide.md) | Panel-by-panel reference: every control, every output, and the calculations behind them |
| [`docs/MINFLUX_analysis_methods.md`](docs/MINFLUX_analysis_methods.md) | Part I: the implementation and maths of every analysis. Part II: a literature review of whether each method suits MINFLUX data, with GREEN / AMBER / RED verdicts |
| [`docs/simulator.md`](docs/simulator.md) | The Simulator tab: what each stage models and where it departs from real MINFLUX |
| [`docs/simulation_module.md`](docs/simulation_module.md) | The ground-truth diffusion simulator: model, validation and the CVE round-trip on D |
| [`docs/sigma_loc_model.md`](docs/sigma_loc_model.md) | The tracking localization-precision model: derivation and immobile-sample calibration |
| [`docs/tracking_sigma_and_cve_plan.md`](docs/tracking_sigma_and_cve_plan.md) | Record of the σ_loc / CVE decisions, including the fixed-cell and bead measurements of tracker wander |
| [`THIRD-PARTY.md`](THIRD-PARTY.md) | External work this builds on, what was taken and what was changed |
| File headers in `helpers_*.R`, `server_*.R`, `R/est_*.R` | Code architecture: each file opens with a map of its contents, data flow and conventions |

---

## Citing MINFLUXR

If you use MINFLUXR in published work, please cite it using the metadata in [`CITATION.cff`](CITATION.cff)
(GitHub's "Cite this repository" button reads it). The time-signature analysis and MSD workflow are R ports
of Bela Vogler's [Toolbox_for_MINFLUX_SPT](https://github.com/Zatyrus/Toolbox_for_MINFLUX_SPT) (GPL-3.0);
please also cite Vogler et al. 2025 when you use them.

---

## Licence

MINFLUXR is licensed under the **GNU General Public License v3.0 or later** ([`LICENSE`](LICENSE)).

GPL-3 is required rather than chosen: the time-signature analysis and the MSD/OLSF workflow are ports of
MinSpt, which is GPL-3.0, and a port of GPL-3 code must itself be GPL-3.

DCR pooling and the per-trace spread statistic are reimplemented in R from
[pyMINFLUX](https://github.com/bsse-scf/pyMINFLUX) (Apache-2.0), which may be combined into a GPL-3 work.
No pyMINFLUX or MinSpt source is copied or bundled here.

[`THIRD-PARTY.md`](THIRD-PARTY.md) records what was taken from each project, what was changed, and the
licences of the R packages MINFLUXR depends on (installed from CRAN, not redistributed).

---

## AI declaration

This app was in parts coded by the author, in collaboration with Claude Code (Opus 4.6, 4.8 and 5). Claude Code was used to optimize, research and design experiments to uncover various implementation possibilities. Claude Code has been implemented in our workflow, to create or revise code, where various functionalities were tested by the authors, before being refined. I started this as a simple MINFLUX viewer, yet AI enabled expansion of the tool to be more profound. We hope this app can be used as a way to explore MINFLUX data, without having to have extensive coding experience.

