# MINFLUXR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24194196.svg)](https://doi.org/10.5281/zenodo.24194196)

An interactive R/Shiny application for visualization and analysis of **Abberior MINFLUX** data.

---

## Statement of need

MINFLUX localizes single fluorophores with nanometre precision, but its data do not behave like
camera-based single-molecule data. Each molecule produces a trace of repeated localizations, and in
tracking data the time steps come in fixed multiples of the cycle time. Localization precision differs
between axes and depends on the photon count at each step. Methods written for camera data ignore
these properties, and the outputs can look right while being biased.

Open tools for MINFLUX data already exist. [pyMINFLUX](https://github.com/bsse-scf/pyMINFLUX) offers a
graphical interface in Python for reading, filtering and visualizing Abberior exports, and
[Toolbox_for_MINFLUX_SPT](https://github.com/Zatyrus/Toolbox_for_MINFLUX_SPT) (Vogler et al. 2025)
provides a tracking workflow. MINFLUXR builds on ideas from both (see [`THIRD-PARTY.md`](THIRD-PARTY.md))
and brings fixed-cell imaging and single-particle tracking analysis together in one R application.

It is written for researchers who have MINFLUX data and want to analyse it without writing code. It
reads Abberior exports directly and covers imaging (filtering, de-blinking, clustering,
colocalization, spatial statistics) and tracking (MSD, diffusion estimators, confinement, state
switching). Each analysis explains what it is suited for, what to be aware of and how its output can
be taken further. A built-in simulator generates data with a known answer, so you can check an
analysis before applying it to real data.

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
the app reports when it is missing.

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
   *About this method* box giving the origin of the method, what was adapted for MINFLUX and
   what to be aware of before using the result.
5. **Save your settings** (Download settings, SMLM tab) to reproduce a figure later: every input that
   changes a number is written to the JSON.

Tracking data: run **Data Quality & Calibration** first. It recovers the cycle time and the
localization precision that the MSD, JDD, confinement, HMM and D-map panels all read.

---

## Result panels and when to trust them

Analysis methods specific to MINFLUX are thinly covered in the literature. Most panels here are
MINFLUXR's own adaptation and interpretation of an established method to MINFLUX data, and some are 
tools built for it. Infoboxes and method sections are considered reading for how MINFLUXR approaches
the data. Always make sure that the data fits a given analysis method.

Whether a method suits your data depends largely on the dataset (track length,
labelling density, photon counts), so each method box states its assumptions and limits, and
leaves the judgement to you. Tools built for this app say so in their method box.

The estimators behind the tracking panels live in `R/est_*.R` and the SMLM computations in
`R/smlm_*.R`, all as plain functions, and `tests/` checks them against ground truth from the simulator (`Rscript tests/run_tests.R`;
slow Monte-Carlo tests with `MINFLUXR_SLOW_TESTS=true`). The longer validation runs and their
results are in `tests/validation/`.

---

## Documentation

| File | Contents |
|------|----------|
| [`docs/user_guide.md`](docs/user_guide.md) | Panel-by-panel reference: every control, every output, and the calculations behind them |
| [`docs/MINFLUX_analysis_methods.md`](docs/MINFLUX_analysis_methods.md) | Part I: the implementation and maths of every analysis. Part II: for each method, the published sources it rests on, what MINFLUXR adapts for MINFLUX data, and the assumptions to check on your own data |
| [`docs/simulator.md`](docs/simulator.md) | The Simulator tab: what each stage models and where it departs from real MINFLUX |
| [`docs/simulation_module.md`](docs/simulation_module.md) | The ground-truth diffusion simulator: model, validation and the CVE round-trip on D |
| [`docs/sigma_loc_model.md`](docs/sigma_loc_model.md) | The tracking localization-precision model: derivation and immobile-sample calibration |
| [`THIRD-PARTY.md`](THIRD-PARTY.md) | External work this builds on, what was taken and what was changed |
| File headers in `helpers_*.R`, `server_*.R`, `R/est_*.R`, `R/smlm_*.R` | Code architecture: each file opens with a map of its contents, data flow and conventions |

---

## Getting help and contributing

- **Questions and bug reports:** open an issue at https://github.com/ViktorMebus/MINFLUXR/issues. Include the MINFLUXR version
  (in `CITATION.cff`), your R version and operating system, the tab and panel involved, and the steps
  that led to the problem. If your data cannot be shared, the Simulator tab can usually produce a file
  that shows the same behaviour.
- **Feature requests:** open an issue describing the analysis you need and, if possible, the paper or
  method it comes from.
- **Code contributions:** see [`CONTRIBUTING.md`](CONTRIBUTING.md).

Everyone taking part is expected to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

---

## Citing MINFLUXR

If you use MINFLUXR in published work, please cite it using the metadata in [`CITATION.cff`](CITATION.cff)
(GitHub's "Cite this repository" button reads it). Each release is archived on Zenodo; the DOI
[10.5281/zenodo.24194196](https://doi.org/10.5281/zenodo.24194196) always points to the latest version, and
the Zenodo page lists a DOI for each individual version. The time-signature analysis and MSD workflow are R ports
of Bela Vogler's [Toolbox_for_MINFLUX_SPT](https://github.com/Zatyrus/Toolbox_for_MINFLUX_SPT) (GPL-3.0);
please also cite Vogler et al. 2025 when you use them.

---

## Licence

MINFLUXR is licensed under the **GNU General Public License v3.0 or later** ([`LICENSE`](LICENSE)).

The MINFLUXR time-signature analysis and the MSD/OLSF workflow are R-adapted ports of Bela Vogler's Toolbox_for_MINFLUX_SPT (GPL-3.0)

DCR pooling and the per-trace spread statistic are R-adapted ports from
[pyMINFLUX](https://github.com/bsse-scf/pyMINFLUX) (Apache-2.0).
No pyMINFLUX or MinSpt source is copied or bundled here.

[`THIRD-PARTY.md`](THIRD-PARTY.md) records what was taken from each project, what was changed, and the
licences of the R packages MINFLUXR depends on (installed from CRAN, not redistributed).

---

## AI declaration

This app was in parts coded by the author, in collaboration with Claude Code (Opus 4.6, 4.8 and 5.5). Claude Code was used to optimize, research and design experiments to uncover various implementation possibilities. Claude Code has been implemented in our workflow, to create or revise code, where various functionalities were tested by the authors, before being refined. Code created by AI, was reviewed by the authors and other AI models. AI was used to create markdowns in the code to clarify its inputs, outputs and functionalities. 

I started this as a simple MINFLUX viewer, yet AI enabled expansion of the tool to be more profound. We hope this app can be used as a way to explore MINFLUX data, without having to have extensive coding experience.

