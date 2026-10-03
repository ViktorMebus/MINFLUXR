# Changelog

Notable changes to MINFLUXR. Versions follow [semantic versioning](https://semver.org);
while the app is below 1.0 the panel layout and the settings JSON may still change between releases.

## Unreleased

### Tracking: headline D and the immobile MSD

- **The headline D is now the upper-regime MSD** (the MSD panel's fit window and settings) in Data
  Quality and Compare Cells, the estimator Vogler et al. 2025 validated against TIRF. The CVE stays
  beside it as the lag-1 check. D − D_w is taken against the headline D.
- **Instrument Calibration plots the immobile MSD** per cell over many cycles, with the noise level 4σ²,
  the fit window and a suggested lower cutoff where the wander has flattened. The wander floor in a
  profile is now the immobile slope **over the MSD fit window**, which matches the headline D; the
  steeper short-lag rise is reported beside it. On fixed cells at L 100.8 nm the two differ about
  tenfold, because the wander saturates within ~20 cycles.
- **The MSD is computed about 25 times faster** (vectorised over index offsets); results are identical,
  checked against the previous code and against a pair-by-pair count.
- **The filter pipeline always starts from the original upload.** After Plot Data it read the rechunked,
  already-filtered data and ran a second time: localizations removed by a range filter left holes that
  the gap split then cut at, so results depended on whether Plot Data had been pressed. The rechunked
  set now only serves the playback chunk selector. Data Quality reads the plotted data like every other
  panel, and the cycle time comes from the original upload.
- **Acquisition Parameters regrouped**: the sequence file (moved here from TSA), calibration, and a
  collapsed *Advanced* group for values that need measurements a first analysis does not have.
- The TSA track selector loads its track IDs on demand (server-side selectize).

### Tracking: datasets

- **Uploads now add to a list of loaded datasets** instead of replacing them. *Active dataset* shows
  one file at a time (the default) or all loaded files merged, and *Remove selected dataset* drops a
  file, so loading the wrong file no longer needs an app restart. Switching rebuilds the view exactly
  as an upload does (chunks, filter ranges, the Filter & Rechunk reset point).
- **Filters now apply to every loaded dataset** and survive switching between them. Each range
  filter is stored as entered and shown clamped to the active dataset's range, so viewing a dataset
  with a narrower range no longer changes it; *Reset Ranges* clears them. The spatial crop stays per
  dataset. Settings files save and restore the stored filters.
- **One filter pipeline** (`R/est_clean.R::filter_pipeline`) now serves the live view, Plot Data's
  rechunk and Compare Cells. The rechunk had its own copy of every stage, which differed from the live
  one (end-trim without the minimum-length floor, a different length-filter rule); results of the live
  view are unchanged, verified stage by stage on real data.
- **The Data Quality CVE now reads the filtered data**, like every other analysis. It read the data
  after the gap split only, so its D and the MSD's D were computed on different data. With the default
  settings (no cleaning) nothing changes; ECO / EFO range filters now move the CVE D as they already
  moved the other panels.
- **Removed the photon–mobility warning.** It read a negative correlation between per-track photons and
  step size as "an ECO/EFO filter removes fast molecules". Low photon counts have several causes
  (bleaching, dim emitters, background, mobility), and one correlation does not identify which.
- **Fixed: 3D CSV re-imports had their axial factor applied twice.** The upload path never returned a
  CSV's z to raw before the session's calibration was applied; only an unused fallback reader did.

### Tracking: compare cells

- **New Compare Cells panel** (Analysis tab) for experiments where each file is one cell. Every file is
  analysed on its own: its own cycle time from the data, the gap split, CVE D, D − D_w, photon-budget
  σ_loc with that file's TCP diameter, pattern, prefactor and wander floor, and the correlated-error
  check. Conditions are compared on one value per cell (mean and t-based 95% CI across cells), so
  correlated tracks from one cell are not counted as independent. Settings are entered in an editable
  table, filled by file name, and saved or loaded as CSV; blank settings follow Acquisition Parameters.
- Compare Cells reads every loaded file directly, whichever dataset is active.
- The per-file pipeline reproduces the Data Quality panel exactly for a single file (`R/est_percell.R`,
  checked on real data and in `tests/testthat/test-percell.R`).
- Cycle-time estimation, the adaptive split threshold and gap splitting moved out of the server into
  `R/est_clean.R`, so the main pipeline and the per-cell comparison share one implementation. Results
  are unchanged (verified value by value on a raw tracking file).

### Tracking: localization precision

- **New Instrument Calibration panel and calibration profiles.** It measures the photon-budget prefactor
  a and the wander floor D_w on an immobile sample, per cell, from the short-lag MSD after the Display
  filters, and keeps the median and quartiles as a named profile recording its sequence (pattern, L,
  dwell, photon limit, damping) and instrument. Profiles save and load as JSON, are applied under
  Acquisition Parameters, and per file in Compare Cells (`calibration` column).
- **Without a calibration, a is the ideal-estimator value of the TCP pattern** (√8 in 2D) and the
  Data Quality readout says "uncalibrated". There is no system-independent value; the previous fixed
  default of 3.6 was one instrument's measurement. Uncalibrated σ is larger than before (×1.27 against
  3.6) until a profile is loaded.
- **The photon-budget prefactor is documented as a measured, effective value.** The previous
  default of 4 was presented as the ideal 2D Cramér-Rao bound; 4 is the one-axis two-beam value, and the
  ideal per-axis bound for any planar pattern is √8 = 2.83. Exported positions come out tighter than
  that (biased estimator, loop-correlated errors), so no Cramér-Rao value bounds the prefactor and it
  has to be calibrated per sequence. The input no longer has an upper limit of 4.
- **The off-center term now depends on the TCP pattern** (new *TCP pattern* setting, read from the
  sequence file): σ² grows as 1 + c·⟨r²⟩/L² with c = 8.6 (hexagon), 9.0 (square), 13.7 (triangle),
  9.0 (octahedron), from the Fisher matrix of an ideal doughnut zero. The previous term,
  1 + 2(x/L)², understated the motion penalty about threefold.
- σ_loc at typical 2D tracking settings rises by about 13% (σ² by about 28%). Results from the
  jump-distance fit, D-map, HMM and confinement null computed before this change should be re-run;
  CVE D and the MSD slope are unaffected.

### Tracking: correlated-error check

- **The smoothing check now uses the single-cycle step duration** (t_c + SW) for its lag-8 pairs. It
  compared against 8 · t_c, and at fast cycles (CFR off, short dwell) 8 · SW exceeds half a cycle,
  so the check silently found no pairs. Results on sequences where it did report are unchanged.
- Its message no longer attributes the signal to galvo damping alone: immobile samples exceed the
  motion bound with damping off, because the tracking loop correlates successive errors.

### Co-localization (SMLM, two channels)

- **New default null for L₁₂ and g₁₂: channel B shifted at random** as a whole (a toroidal shift,
  Lotwick & Silverman 1982), restricted to shifts that keep it in the cell footprint and at least as long
  as the largest radius analysed. Each channel keeps its own clustering and only their registration is
  broken, which is the question "are the two placed independently". The previous default, random
  labelling, answers a different question: on clustered channels it never detected two species sharing
  the same clusters and read two unrelated clustered species as exclusion. Random labelling stays
  available for the question it does answer (do the two channels differ in arrangement). Measured in
  `tests/validation/coloc_null_clustered.R`.
- The plot caption now names the null each p-value was computed against.
- L₁₂ / g₁₂ results computed with the old default should be re-run before they are reported.

### Tracking: motion blur

- **Motion blur is now computed from the acquisition timing** (`R/est_blur.R`) and one value is used
  by the CVE, MSD fits, jump-distance fit, D-map, change-point segmentation, HMM and confinement null.
  R is set by when a localization's photons are collected: pattern, vertex order, dwell, pattern repeats,
  passes to the photon limit, beam moves and the CFR centre dwell, read from the sequence file when one is
  loaded and from the collection window eco/efo otherwise. The TCP diameter L does not enter.
- The previous constant, MinSpt's R = 1/6.2, is about 50% above the value the timing gives on the
  standard 2D tracking sequence (0.152 vs 0.100 per step). On simulated MINFLUX localizations it made the
  single-step D estimates (jump distance, D-map, change-point, HMM) read 15–17% high; the timing model
  reads within 2%. **Expect those D values to drop by roughly 12% on re-analysis.** The CVE D does not
  depend on R and does not change; the MSD slope does not change, the σ read from its intercept does.
- The Data Quality panel reports the blur value and where it came from, and warns when the sequence file
  has damping switched on, which smooths positions across steps and is not a blur that any R describes.
- A sequence file that does not fit the data (its collection alone outlasts the measured step, as an
  imaging sequence does against tracking data) or whose final pattern is axial-only is not used for the
  blur; the panel says why and takes the value from the data instead.

### Multiplex

- Channels of a multiplex can be aligned by nudging one channel in x, y (and z for 3D data). Offsets
  are stored per source file or time window, re-applied when the multiplex is rebuilt, and saved in the
  settings JSON. The panel warns against aligning two different targets onto each other, which would
  build the co-localization into the data.
- The multiplex recipe (source files or time windows, channel labels and colours) is saved with the
  settings and rebuilds on load. A CSV exported from a multiplex keeps its channels and loads back as
  one.

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
- The HMM fits at most a set number of steps (default 150,000), drawing whole tracks at random
  when the data hold more, and says so in the summary. Recovery of D is unchanged at that budget
  (`tests/validation/hmm_step_budget.csv`). Fitting one state fewer and one more for BIC is now
  opt-in, since it is two further fits.
- The HMM derives the photon-budget σ_loc itself when the calibration panel has not published one,
  instead of refusing to run.
- The jump-distance, D-map and change-point panels now require the cycle time rather than falling
  back to no blur correction, which would have read D about a third low without saying so.

### Spatial statistics (SMLM)

- Ripley L, g(r), L₁₂ and g₁₂ now report a **global (studentized maximum-deviation) p-value**
  alongside the envelope band. The pointwise band alone was never a 5% test: on random data it
  was crossed almost always, where the global test now sits near its nominal level.
- The bivariate envelopes can use a **random-labelling** null in addition to CSR-in-window. (Superseded
  as the default before release; see Unreleased.)
- The analysis window is estimated from the localizations by occupancy and morphological closing,
  and its area is debiased, so a sparse pattern is no longer measured against a window that hugs
  its own points.
- Hotspot statistics: whole-quadrat variance-to-mean ratio (edge quadrats were biasing it), an
  eps floor in the DBSCAN calibration, and a density-aware focus null. Hotspot numbers from
  before this release are not comparable with the current ones.

### Reading the tracking panels

- The **confinement index is the primary call on confinement**. The displacement-vs-diffusion panel
  now reports how far it agrees with it track by track, instead of a second percentage that invited
  comparison: the two ask different questions (the shape of the MSD curve against a null matched on
  each track's length and D, versus net displacement against an infinite-plane null at one ensemble D).

### Quality control

- **D by track continuity** (Data Quality & Calibration): splits the analysed tracks into quarters by
  the continuity score and gives each its own D. It reports the gradient across molecules and,
  separately, what the most continuous quarter would do to the pooled D. On MINFLUX tracking data a
  fast molecule collects dark cycles, so a continuity threshold can select on speed.
- Trimming degraded track ends no longer shortens a segment below the gap-splitting minimum length.

- **Fiducial detection** is now driven by trace lifetime relative to the file's own distribution,
  with brightness as a secondary test, matching how a bead behaves: the tracker holds onto it long
  after real fluorophores have bleached.
- The spread filter offers a robust MAD variant alongside the pyMINFLUX-matching SD.
- Per-file summary export for batch runs.

### Documentation, licensing and tests

- Documentation and in-app text describe what each tool does and where it stops being reliable.
  Help text in the Tracking tab is one or two sentences per control, with the longer explanations
  in `docs/user_guide.md` and `docs/MINFLUX_analysis_methods.md`.
- The lateral/axial precision warning fires in both directions. It previously warned only when the
  axial precision was worse, although in MINFLUX it is usually the tighter of the two.
- The displacement panel's resolution floor now states the right direction: confinement on a scale
  larger than the free extent is what the test cannot see.
- The SMLM and simulator tabs got the same pass. The documentation now says plainly that most tools
  are MINFLUXR's own adaptation of established methods, not a community standard.
- The fiducial scan is rated amber: it is tool-original and has not been validated by other groups.
- The KDE co-localization caveat now follows the window the randomisation used. Inside the sample
  footprint the test is conservative; over the bounding box it rejects too readily.
- Help text for the bivariate L₁₂ and g₁₂ envelopes describes the null they use. Help text for Ripley's "all files" mode now says it draws one curve per file.

- Licensed **GPL-3.0-or-later** (`LICENSE`), required by the MinSpt port. `THIRD-PARTY.md`
  records what was taken from MinSpt and pyMINFLUX and what was changed.
- `CITATION.cff` with full author metadata and the works to cite alongside it.
- The README is now installation and orientation; the panel-by-panel reference moved to
  `docs/user_guide.md`.
- `tests/` holds a testthat suite (`Rscript tests/run_tests.R`, slow Monte-Carlo tests behind
  `MINFLUXR_SLOW_TESTS=true`) and `tests/validation/` the longer runs with their result CSVs.
