# Changelog

Notable changes to MINFLUXR. Versions follow [semantic versioning](https://semver.org);
while the app is below 1.0 the panel layout and the settings JSON may still change between releases.

## Unreleased

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
