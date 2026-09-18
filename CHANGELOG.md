# Changelog

Notable changes to MINFLUXR. Versions follow [semantic versioning](https://semver.org);
while the app is below 1.0 the panel layout and the settings JSON may still change between releases.

## 0.6.0 — 2026-09-18

First public release (Zenodo).

### Tracking estimators

- The diffusion estimators now live in `R/est_*.R` as plain functions (CVE, MSD/OLSF,
  jump-distance, D-map, change-point, diffusive-state HMM). The server panels gather inputs and
  call them, so each one can be run and tested without a running app.
- **Motion blur** is now applied in the jump-distance fit, the D-map and the change-point
  segmentation, which previously read D about 30% low at MINFLUX cycle times. Recovery on
  simulated tracks is now within a few percent.
- The MSD intercept converts to σ_loc with the blur offset added back; without it σ was
  under-read by roughly a third.
- **The diffusive-state HMM was rebuilt** as a continuous-time model fitted by maximum
  likelihood on the irregular Δt (`expm(QΔt)`), replacing the previous path-based rate estimate.
  It reports the range over which it recovers D and dwell time, and flags dwell times too short
  to resolve and tracks with too few switches to constrain the rate.
- Change-point detection uses a variance-scaled penalty, so the segment count no longer depends
  on the magnitude of D.

### Spatial statistics (SMLM)

- Ripley L, g(r), L₁₂ and g₁₂ now report a **global (studentized maximum-deviation) p-value**
  alongside the envelope band. The pointwise band alone was never a 5% test: on random data it
  was crossed almost always, where the global test now sits near its nominal level.
- The bivariate envelopes can use a **random-labelling** null in addition to CSR-in-window, which
  is the correct null when both channels share the same support.
- The analysis window is estimated from the localizations by occupancy and morphological closing,
  and its area is debiased, so a sparse pattern is no longer measured against a window that hugs
  its own points.
- Hotspot statistics: whole-quadrat variance-to-mean ratio (edge quadrats were biasing it), an
  eps floor in the DBSCAN calibration, and a density-aware focus null. Hotspot numbers from
  before this release are not comparable with the current ones.

### Quality control

- **Fiducial detection** is now driven by trace lifetime relative to the file's own distribution,
  with brightness as a secondary test, matching how a bead behaves: the tracker holds onto it long
  after real fluorophores have bleached.
- The spread filter offers a robust MAD variant alongside the pyMINFLUX-matching SD.
- Per-file summary export for batch runs.

### Documentation, licensing and tests

- Licensed **GPL-3.0-or-later** (`LICENSE`), required by the MinSpt port. `THIRD-PARTY.md`
  records what was taken from MinSpt and pyMINFLUX and what was changed.
- `CITATION.cff` with full author metadata and the works to cite alongside it.
- The README is now installation and orientation; the panel-by-panel reference moved to
  `docs/user_guide.md`.
- `tests/` holds a testthat suite (`Rscript tests/run_tests.R`, slow Monte-Carlo tests behind
  `MINFLUXR_SLOW_TESTS=true`) and `tests/validation/` the longer runs with their result CSVs.
