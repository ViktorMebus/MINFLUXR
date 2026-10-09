# Changelog

Notable changes to MINFLUXR. Versions follow [semantic versioning](https://semver.org);
while the app is below 1.0 the panel layout and the settings JSON may still change between releases.

## Unreleased

### Tracking: headline D and the immobile MSD

- **The headline D is now the upper-regime MSD** (the MSD panel's fit window and settings) in Data
  Quality and Compare Cells. The CVE stays beside it as the lag-1 check. D − D_w is taken against the
  headline D.
- **Instrument Calibration plots the immobile MSD** per cell over many cycles, with the noise level 4σ²
  and the fit window. The wander floor in a profile is now the immobile slope **over the MSD fit
  window**, matching the headline D; the short-lag rise is reported beside it.
- **The MSD is computed about 25 times faster** (vectorised over index offsets); results are identical.
- **The filter pipeline always starts from the original upload.** After Plot Data it previously re-ran
  on the rechunked, already-filtered data, so results depended on whether Plot Data had been pressed.
  The rechunked set now only serves the playback chunk selector. Data Quality reads the plotted data
  like every other panel, and the cycle time comes from the original upload.
- **Acquisition Parameters regrouped**: the sequence file (moved here from TSA), calibration, and a
  collapsed *Advanced* group.
- The TSA track selector loads its track IDs on demand (server-side selectize).

### Tracking: datasets

- **Uploads add to a list of loaded datasets** instead of replacing them. *Active dataset* shows one
  file at a time (the default) or all loaded files merged, and *Remove selected dataset* drops a file.
  Switching rebuilds the view as an upload does (chunks, filter ranges, the Filter & Rechunk reset
  point).
- **Filters apply to every loaded dataset** and survive switching between them. Each range filter is
  stored as entered and shown clamped to the active dataset's range; *Reset Ranges* clears them. The
  spatial crop stays per dataset. Settings files save and restore the stored filters.
- **One filter pipeline** (`R/est_clean.R::filter_pipeline`) serves the live view, Plot Data's rechunk
  and Compare Cells. The rechunk previously had its own copy of every stage.
- **The Data Quality CVE reads the filtered data**, like every other analysis. ECO / EFO range filters
  now move the CVE D as they move the other panels.
- **Removed the photon–mobility warning**, which read a correlation between per-track photons and step
  size as "an ECO/EFO filter removes fast molecules".
- **Fixed: 3D CSV re-imports had their axial factor applied twice.**

### Tracking: compare cells

- **New Compare Cells panel** (Analysis tab) for experiments where each file is one cell. Each file is
  analysed on its own: cycle time from the data, gap split, CVE D, D − D_w, and photon-budget σ_loc with
  that file's TCP diameter, pattern, prefactor and wander floor. Conditions are compared on one value
  per cell (mean and t-based 95% CI across cells). Settings are entered in an editable table, filled by
  file name, and saved or loaded as CSV; blank settings follow Acquisition Parameters.
- Compare Cells reads every loaded file directly, whichever dataset is active.
- Cycle-time estimation, the adaptive split threshold and gap splitting moved into `R/est_clean.R`, so
  the main pipeline and Compare Cells share one implementation (`R/est_percell.R`).

### Tracking: localization precision

- **New Instrument Calibration panel and calibration profiles.** It measures the photon-budget prefactor
  a and the wander floor D_w on an immobile sample, per cell, from the short-lag MSD after the Display
  filters, and keeps the median and quartiles as a named profile recording its sequence (pattern, L,
  dwell, photon limit, damping) and instrument. Profiles save and load as JSON, are applied under
  Acquisition Parameters, and per file in Compare Cells (`calibration` column).
- **Without a calibration, a is the ideal value for the TCP pattern** (√8 in 2D) and Data Quality says
  "uncalibrated". The previous default was 3.6, so uncalibrated σ is ×1.27 larger until a profile is
  loaded.
- **The prefactor input no longer has an upper limit of 4.** 4 is the one-axis two-beam value; the
  ideal per-axis value for a planar pattern is √8.
- **The off-center term depends on the TCP pattern** (new *TCP pattern* setting, read from the sequence
  file): σ² grows as 1 + c·⟨r²⟩/L² with c = 8.6 (hexagon), 9.0 (square), 13.7 (triangle), 9.0
  (octahedron). The previous term was 1 + 2(x/L)².
- σ_loc rises at typical tracking settings. Re-run the jump-distance fit, D-map, HMM and confinement
  null; CVE D and the MSD slope are unaffected.

### Tracking: calibration accuracy

- **The calibration's σ includes the motion-blur term**, which lowers the short-lag MSD intercept by
  8R·D·Δt (Berglund 2010). R comes from each file's photon collection window; the summary shows how far
  a moves for R ±30%.
- **Removed the suggested lower cutoff** for the MSD fit window; the immobile-MSD plot remains.
- **Removed the GOOD / MODERATE / POOR verdict** on the single-cycle fraction from Data Quality, the
  Time Signature Analysis summary and the cycle status under Acquisition Parameters; the fractions
  (raw and filtered) are still reported.

### Tracking: increment-correlation checks removed

- **Removed the "increments … too weak: correlated error" check** from Data Quality and Compare Cells.
  With motion blur the lag-1 step covariance is cov₁ = 2R·D·Δt − σ² (Vestergaard 2014), so a positive
  r₁ is expected and does not indicate smoothing. Galvo damping is read from the sequence file instead.
- **Removed the "noise share of step variance" line and its warning** (−2·r₁ assumes no blur).
- The "uncalibrated" status and the diagnostic CVE σ are no longer shown in warning colours.

### Co-localization (SMLM, two channels)

- **New default null for L₁₂ and g₁₂: channel B shifted at random** as a whole (a toroidal shift,
  Lotwick & Silverman 1982), restricted to shifts that keep it in the cell footprint and at least as long
  as the largest radius analysed. Each channel keeps its own clustering and only their registration is
  broken. Random labelling stays available for asking whether the two channels differ in arrangement.
- The plot caption names the null each p-value was computed against.
- L₁₂ / g₁₂ results computed with the old default should be re-run.

### Tracking: motion blur

- **Motion blur is computed from the acquisition timing** (`R/est_blur.R`), and one value is used by the
  CVE, MSD fits, jump-distance fit, D-map, change-point segmentation, HMM and confinement null. R is set
  by when a localization's photons are collected (pattern, vertex order, dwell, pattern repeats, passes
  to the photon limit, beam moves, CFR centre dwell), read from the sequence file when one is loaded and
  from the collection window eco/efo otherwise.
- This replaces the constant R = 1/6.2 from MinSpt. Single-step D values (jump distance, D-map,
  change-point, HMM) and the σ read from the MSD intercept change on re-analysis; the CVE D and the MSD
  slope do not.
- Data Quality reports the blur value and its source, and warns when the sequence file has damping
  switched on.
- A sequence file that does not fit the data, or whose final pattern is axial-only, is not used for the
  blur; the panel says why and takes the value from the data.

### SMLM: ground-truth estimator

- **The trace count per molecule follows the Simulator's labelling model exactly** (`gt_tid_model`):
  multiplicity, variance and detection probability. The multiplicity and η change slightly for
  multivalent presets.
- **New consistency test** between the spatial count and the trace total, reported as *under-resolved*
  or *over-split*. Under-resolved now sets the crowding flag; the previous flags missed molecules that
  de-blinking had already merged. Replaces the fixed "counts differ by more than 2×" rule.
- **The count range includes each count's sampling error** (binomial detection for the spatial count,
  the trace-count variance for the stoichiometric one).
- Candidates that hold more traces than one molecule plausibly gives are marked in the hover text.

### SMLM: one set of default settings

- **Every setting that changes SMLM results now has one default** (`R/smlm_settings.R`). The widgets,
  the helper functions and the reset buttons all read it. They previously disagreed:
  - Median clustering MinPts was 15 in the UI, 3 in the helper functions, 8 after *Reset* and 5 in
    the kNN-distance plot.
  - The majority-vote de-blink settings were 0.3 / 3.0 in the UI and 0.2 / 1.2 in the helpers.
  - *Reset* on the Median Cluster settings set eps to 150 nm; a new session shows 200 nm.
- *Reset* now restores the values a new session shows: eps 200 nm, MinPts 15.
- `smlm_settings()` builds the same settings from a script, with checked types and limits:
  `apply_clustering_pipeline(df, smlm_settings(enable_clustering = TRUE))`.
- The two-colour reset restores channel A's colour to the UI default (#1f77b4; it was #1f77a4).
- The filter, channel, de-blink and median-clustering functions take a `settings` argument
  (`smlm_settings()` from a script, or the app's input). Small helpers take plain arguments:
  `channels_in_scope(df, scope)`, `deblink_minpts(tid_length_range)`, `build_center_args(X, mode, lambda)`.
  Results are unchanged.
- The per-file cache of filtered and clustered data is keyed on the settings object that the
  clustering actually runs with, plus the data version, thinning and fiducial exclusions. A setting
  can no longer be left out of the key.
- **Fixed: loading a settings file did not restore de-blinking, median clustering, the two-colour
  split or the axial calibration.** The file held most of them, but the loader restored only the
  filters, crop and QC options, so a reloaded analysis clustered with whatever the session had.
  Ten settings were also never saved (the Median Cluster method and its options, and the density
  unit). Settings files now carry every setting that changes results, and loading restores all of
  them; an older file restores what it holds.
- The plot helpers take display options (colour mode, palette, medians, density radius, channel
  colours) as plain arguments with one set of defaults (`SMLM_DISPLAY_DEFAULTS`), kept apart from
  the settings that change results. No helper reads the app's input directly any more.
- Removed three unused helper functions.

### Multiplex: plain dataset name

- The built multiplex dataset is called "Multiplex" (it was "★ Multiplex"). The name is part of
  downloaded file names (settings, CSV export), where the star was rejected or garbled on some
  systems.

### Tests

- SMLM test files for de-blinking, median clustering, density normalisation, Ripley's L and g(r),
  local enhancement, two-channel co-localization, colour separation and the ground-truth estimator,
  each against simulated data with a known answer.

### Documentation

- README and docs/ shortened and made sample-neutral. A notation table opens the methods document, and
  every formula's symbols are defined where it appears.
- Formulas brought in line with the code: the motion-blur offset R_abs in the MSD, JDD, D-map and HMM
  expressions; the √8 prefactor in the simulators; the user guide's HMM section describes the
  continuous-time maximum-likelihood model.
- The JDD panel's formula line shows R_abs, which is what the fit uses.

### Multiplex

- Channels of a multiplex can be aligned by nudging one channel in x, y (and z for 3D data). Offsets
  are stored per source file or time window, re-applied when the multiplex is rebuilt, and saved in the
  settings JSON. The panel warns against aligning two different targets onto each other.
- The multiplex recipe (source files or time windows, channel labels and colours) is saved with the
  settings and rebuilds on load. A CSV exported from a multiplex keeps its channels and loads back as
  one.

## 0.6.0 — 2026-09-18

First public release (Zenodo).

### Tracking estimators

- The diffusion estimators live in `R/est_*.R` as plain functions (CVE, MSD/OLSF, jump-distance, D-map,
  change-point, diffusive-state HMM). The server panels gather inputs and call them.
- **Motion blur** is applied in the jump-distance fit, the D-map and the change-point segmentation.
- The MSD intercept converts to σ_loc with the blur offset added back.
- **The diffusive-state HMM was rebuilt** as a continuous-time model fitted by maximum likelihood on the
  irregular Δt (`expm(QΔt)`), replacing the path-based rate estimate. It flags dwell times too short to
  resolve and tracks with too few switches to constrain the rate.
- Change-point detection uses a variance-scaled penalty, so the segment count does not depend on the
  magnitude of D.
- The HMM fits at most a set number of steps (default 150,000), drawing whole tracks at random when the
  data hold more, and says so in the summary. Fitting one state fewer and one more for BIC is opt-in.
- The HMM derives the photon-budget σ_loc itself when the calibration panel has not published one.
- The jump-distance, D-map and change-point panels require the cycle time instead of falling back to no
  blur correction.

### Spatial statistics (SMLM)

- Ripley L, g(r), L₁₂ and g₁₂ report a **global (studentized maximum-deviation) p-value** alongside the
  envelope band.
- The bivariate envelopes can use a **random-labelling** null in addition to CSR-in-window. (Superseded
  as the default; see Unreleased.)
- The analysis window is estimated from the localizations by occupancy and morphological closing, and
  its area is debiased.
- Hotspot statistics: whole-quadrat variance-to-mean ratio, an eps floor in the DBSCAN calibration, and
  a density-aware focus null. Hotspot numbers from earlier versions are not comparable.

### Reading the tracking panels

- The **confinement index is the primary call on confinement**. The displacement-vs-diffusion panel
  reports how far it agrees with it track by track instead of a second percentage.

### Quality control

- **D by track continuity** (Data Quality & Calibration): splits the analysed tracks into quarters by
  continuity score and gives each its own D, with the gradient across quarters and the effect of the
  most continuous quarter on the pooled D.
- Trimming degraded track ends no longer shortens a segment below the gap-splitting minimum length.
- **Fiducial detection** is driven by trace lifetime relative to the file's own distribution, with
  brightness as a secondary test.
- The spread filter offers a robust MAD variant alongside the pyMINFLUX-matching SD.
- Per-file summary export for batch runs.

### Documentation, licensing and tests

- Documentation and in-app text describe what each tool does and where it stops being reliable. Help
  text is one or two sentences per control, with longer explanations in `docs/user_guide.md` and
  `docs/MINFLUX_analysis_methods.md`.
- The lateral/axial precision warning fires in both directions.
- The documentation states that most tools are MINFLUXR's own adaptation of established methods.
- The fiducial scan is rated amber.
- The KDE co-localization caveat follows the window the randomisation used.
- Help text for the L₁₂ and g₁₂ envelopes describes the null they use; Ripley's "all files" mode draws
  one curve per file.
- Licensed **GPL-3.0-or-later** (`LICENSE`), required by the MinSpt port. `THIRD-PARTY.md` records what
  was taken from MinSpt and pyMINFLUX and what was changed.
- `CITATION.cff` with full author metadata and the works to cite alongside it.
- The README covers installation and orientation; the panel-by-panel reference is in
  `docs/user_guide.md`.
- `tests/` holds a testthat suite (`Rscript tests/run_tests.R`, slow Monte-Carlo tests behind
  `MINFLUXR_SLOW_TESTS=true`) and `tests/validation/` the longer runs with their result CSVs.
