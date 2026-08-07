# MINFLUXR

An interactive R/Shiny application for visualization and analysis of **Abberior MINFLUX** data.

---

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Launching the App](#launching-the-app)
- [Loading Data](#loading-data)
- [MINFLUX Data Columns Reference](#minflux-data-columns-reference)
- [SMLM Display Tab](#smlm-display-tab)
  - [Data Upload](#data-upload)
  - [Axial (z) calibration](#axial-z-calibration)
  - [Metadata Filtering](#metadata-filtering)
  - [Visualization Options](#visualization-options)
  - [De-blink correction](#de-blink-correction)
  - [Nearest Neighbor Analysis](#nearest-neighbor-analysis)
  - [Median Clustering](#median-clustering)
  - [Background image overlay](#background-image-overlay)
  - [Heatmap overlay](#heatmap-overlay)
  - [2-color MINFLUX](#2-color-minflux)
  - [Multiplexed MINFLUX](#multiplexed-minflux)
  - [Estimate Ground-Truth Positions](#estimate-ground-truth-positions)
  - [Global Reference Plane and Surface](#global-reference-plane-and-surface)
- [SMLM Analysis Tab](#smlm-analysis-tab)
  - [Panels](#panels)
  - [Two-colour method details](#two-colour-method-details)
  - [Export](#export)
- [Tracking Display Tab](#tracking-display-tab)
- [Tracking Analysis Tab](#tracking-analysis-tab)
- [MINFLUX Simulator Tab](#minflux-simulator-tab)
- [Key Calculations](#key-calculations)
- [Plugin System](#plugin-system)
- [Further Documentation](#further-documentation)
- [Notes & Tips](#notes--tips)

## Overview

MINFLUXR provides three main analysis environments:

| Tab | Purpose |
|-----|---------|
| **SMLM** | Static super-resolution data: filtering, de-blinking, clustering, co-localization, spatial statistics |
| **Tracking** | Single-molecule tracking: trajectory visualization, MSD, diffusion analysis |
| **MINFLUX Simulator** | Generate synthetic MINFLUX datasets — labelled structures, whole acquisitions, and ground-truth diffusion tracks |

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

Key dependencies include: `shiny`, `shinyWidgets`, `shinyjs`, `plotly`, `dplyr`, `dbscan`, `RANN`, `mclust`, `magick`, `alphashape3d`, `data.table`, `qs`, `fst`, `readr`, `changepoint`, `minpack.lm`.

> **macOS note:** `qs` and `fst` use C++/OpenMP. If they fail to compile, install Xcode Command Line Tools (`xcode-select --install`) and retry, or the app will fall back to `fst`-only chunk storage automatically.

---

## Launching the App

```r
shiny::runApp("path/to/MINFLUXR")
```

Or open `app.R` in RStudio and click **Run App**. The maximum upload size is set to **2 GB** by default.

---

## Loading Data

The importer handles four sources and normalises them to one internal form:

| Source | Notes |
|--------|-------|
| **`.npy`** | Abberior's **tabular** binary export. The primary format: 5–10× smaller than the equivalent JSON and much faster to read. |
| **`.json`** | Imspector export. Same content, much larger. |
| `.csv` (pyMINFLUX) | Already axially corrected by pyMINFLUX; the app detects and removes that factor on read. |
| `.csv` (this app) | Our own export. Carries `z_scale_applied` and `dataset_tag`, so its axial factor is known **exactly** rather than assumed. |

### About the `.npy` reader

`.npy` is read by a **hand-written pure-R NumPy decoder** (`read_npy()` in `helpers_shared.R`) — there
is no Python or `reticulate` dependency. It decodes the header, the structured dtype descriptor and
the raw data block directly.

Two constraints follow from that:

- **Only the newer, flat "tabular" MINFLUX format is supported.** The **legacy nested `m2205` dtype is
  not read** — the app stops with a clear error asking you to export the tabular format instead.
  Nested structured dtypes have no flat column layout to decode.
- Fortran-ordered arrays are rejected (re-save C-ordered). A plain unnamed array is rejected too: the
  reader needs the named fields (`loc`/`lnc`, `tid`, `tim`, `efo`, `eco`, `cfr`, `dcr`, `vld`, `fnl`, …).

`.npy` and `.json` converge on the same post-processing: `vld = FALSE` records are dropped, iterations
are collapsed to one row per localization, and CFR is carried forward. The one difference is how the
final iteration is chosen — `.npy` uses the `fnl` flag when present, `.json` uses `itr == max(itr)`.

Multiple files can be loaded and switched between with the *Active file* selector. `.npy` and `.json`
are read with **raw z** and get the axial factor applied afterwards; CSVs are first returned to raw so
every source lands on the same footing (see *Axial (z) calibration*). **DCR pooling** is on by
default for `.npy`/`.json`, adding an ECO-weighted `dcr_pooled` column.

---

## MINFLUX Data Columns Reference

The following columns are read from Abberior exports. Understanding them is useful for interpreting the filter controls.

| Column | Full name | Typical units | Description |
|--------|-----------|---------------|-------------|
| `x`, `y`, `z` | Position | nm | Localization coordinates |
| `tid` | Trace ID | integer | Groups localizations belonging to the same fluorophore event |
| `itr` | Iteration | integer | MINFLUX iteration step; only the final valid iteration is kept by default |
| `tim` | Time | s (on disk), ms (tracking tab) | Absolute acquisition timestamp |
| `efo` | Effective Frequency at Offset | Hz | Photon emission rate at circumference of the beampattern (like eco, but in frequency Hz). EFO = ECO/(n·t_dwell), with n the number of TCP cycles |
| `eco` | Effective Counts at Offset | photons | Summarized photon count at circumference of the beampattern of all TCP cycles |
| `cfr` | Center Frequency Ratio | 0–1 | Ratio of Effective Frequency at Center (EFC) versus Effective Frequency at Offset (EFO). The ratio of emitters from center to periphery of the beampattern |
| `dcr` | Detector channel ratio | 0–1 | Fraction of emitters detected by 2 open detectors (i.e. Cy5 near and Cy5 far); used for 2-color separation |
| `dwell` | Dwell time | steps | Derived: number of MINFLUX steps per localization (see [Dwell Time](#dwell-time)) |
| `fluo` | Fluorophore channel | 1 or 2 | Groups assigned by the pyMINFLUX pipeline for 2-color MINFLUX |
| `gri` | Grid index | integer | Identifies which scan grid the localization was acquired on |
| `bot` / `eot` / `fnl` | Begin/End/Final of trace flags | logical | Trajectory boundary markers |

---

## SMLM Display Tab

Sidebar groups appear in this order. **Plot Data** (sticky button at the top) applies filters,
calibration and clustering and re-renders — settings do not live-update the plot.

### Data Upload

Single- or multi-file loader accepting **`.npy`, `.json` and `.csv`**. Multiple files can be loaded
simultaneously; use the **Active file** selector to switch between them, and **Remove Selected File**
to drop one. **Lock filters** preserves current filter values when switching files.

- **Pool DCR across iterations** (on by default for `.npy` / `.json`) adds an ECO-weighted
  `dcr_pooled` column, a less noisy DCR summary than the final iteration alone.
- **Display mode** toggles 2D / 3D. In 2D, **Viewing plane** selects XY, XZ or YZ; XZ and YZ are
  greyed out (not hidden) when the dataset has no axial variation, so it stays visible *why* they are
  unavailable.

### Axial (z) calibration

A per-dataset refractive-index correction applied to z. The default of **0.7 is the pyMINFLUX
standard**; change it if your refractive index differs. Files are read with raw z and
the factor is applied once, **on Plot Data** — so a saved analysis reproduces. Loading one of our own
exports adopts the factor recorded in the file. It multiplies every axial quantity (z spread, cluster
extent in z, 3D neighbour distances, Rg), so it is a calibration to measure on a known standard, not
a display preference. Greyed out when the dataset has no axial variation.

### Metadata Filtering

Filters are applied sequentially in the order listed in the sidebar. All filters operate on the currently active dataset.

```
Raw data
 └── [1] Metadata filters  (EFO, CFR, DCR, TIM, ECO, dwell, fluo)
      └── [2] TID length filter  (number of points per trace)
           └── [3] TID spread filter  (per-trace spread of localizations)
                └── [4] Rogue-trace rule (optional)
```

**Metadata filters.** Pick a variable, see its histogram (with the unfiltered data behind it in grey
if *Show raw data* is on), and set a numeric range:

- `EFO` — Effective Frequency at Offset. The averaged emission frequency over all points of the scan
  pattern except the centre; exceedingly high values may indicate multi-emitters.
- `CFR` — Center Frequency Ratio, EFC/EFO. Commonly used to exclude multi-emitter events; **no CFR
  cutoff is applied automatically anywhere in the app** — set it yourself if you want one.
- `TIM` — acquisition timestamp; useful for selecting a time window
- `DCR` — Detector Channel Ratio; photons at detector 1 vs detector 2, used to separate two labels in
  2-colour MINFLUX
- `ECO` — Effective photon Count at Offset: the raw photon count per localization at the periphery of
  the donut
- `Dwell` — dwell time of the localization, in MINFLUX steps
- `FLUO` — fluorophore channel index; a pyMINFLUX-exported column, usable here when present

**TID length filter**: keeps only traces with a number of localizations within the specified range.
Useful for excluding very short (potentially noise) or very long traces (potentially fiducials or
moving emitters).

**Inspect trace spread (σ̃)** opens the spread panel:

- *Min localizations per trace* — traces shorter than this are excluded from the spread histogram,
  the median-spread labels, the allowable ranges and the σ filters. Literature typically uses ≥ 5
  (an SD needs ≥ 2 points to exist at all).
- The histogram shows per-TID spreads, optionally overlaying the raw distribution and the median
  spread before/after filtering. Labels report median(σ̃) ± nMAD (nMAD = 1.4826 × MAD, matching
  pyMINFLUX).
- *Set spread ranges to the 95th percentile* seeds the σ filter ranges from the data.
- **TID spread filter**: per-axis x/y/z bounds. Setting the upper bound low retains only spatially
  tight traces — useful for removing gold beads or freely diffusing molecules. See
  [Trace Spread](#trace-spread-sd-filter).
- **Remove rogue traces by rule** — an alternative to hand-set bounds. Three rules, each with a
  multiplier *k*: spread > k × dataset median, spread > median + k × MAD (robust), or any single axis
  spread > k × that axis's median. Dropped traces are removed **before** grouping, clustering, fits
  and surfaces.

**Download / Load Filter Settings** saves and restores all filter parameters as JSON.
**Export MINFLUX data (.csv)** writes the currently filtered localizations.

### Visualization Options

- **Point size**
- **Color points by** — Trace/Cluster colour (default), Multiplex channel, Unicolor (one colour for
  everything, with a picker), Depth (z), Time (tim), or Trace-median spatial density. The three
  numeric modes get a colorscale picker and a reverse toggle; density mode adds a radius (nm) that
  counts how many other *trace-median* positions fall within it.
- **Show trace median positions** — one point per TID (or per de-blinked molecule / cluster when
  those are active) instead of individual localizations
- **Dark mode**, **Hide gridlines & axis titles**
- **Add scale bar** (2D) — length, thickness and font size in nm
- **Transform coordinates to positive space** — shifts the dataset so all coordinates are positive
- **Show Alpha Shape (3D)** — wraps each trace (or each de-blinked ID) in a 3D alpha surface, with an
  opacity slider
- **Crop to view** (2D) — zoom into a region and press the button; the X/Y crop is retained when you
  switch to 3D
- **Crop to region (3D)** — an explicit trimmer with three shapes: axis-aligned box, sphere, or
  PCA-aligned ellipsoid
- **Auto-Rotate & Record (3D)** — a spinning animation (speed, elevation, camera distance) and a WebM
  recorder with width/height, FPS and bitrate controls. Start auto-rotate first, then record.

### De-blink correction

Combats the instrument artifact where one fluorophore that blinks several times produces several trace
IDs stacked on top of each other, overcounting the number of molecules. De-blinking merges them back
into one **de-blinked ID** with a shared centroid, so downstream analyses count molecules rather than
blinking events.

Works at **trace level** — a whole trace joins a molecule or it does not, so a trace is never split.

Two criteria:

- **Trace centroid distance** (default) — merge when two traces' median positions are closer than the
  threshold. Set it in nm (`Centroid distance threshold`) or as **k × the measured trace spread**;
  changing k rewrites the distance.
- **Majority vote of localizations** — merge when a fraction *f* of a trace's localizations sit closer
  to a neighbouring trace's centroid than to their own. No distance in nm. Useful values of *f* are
  well below 0.5 (two traces of the *same* molecule score ≈0.5 by construction). A *distance
  tolerance* ρ loosens the "strictly closer" test, and **Agreement** decides whether both traces must
  agree (mutual) or either is enough.

**De-blinking only merges — it never deletes.** A trace that cannot be merged becomes its own
single-blink molecule, so every trace that survived filtering is still counted. Removing data is the
job of the metadata and spread filters, where it is visible and reversible.

**Choosing the threshold.** *Measure threshold from data* sets it from your own file: it measures the
per-axis trace spread and finds the smallest k at which no unmerged re-blinks remain. Two diagnostics
sit beside it — a density plot whose orange curve should show two humps (re-blinks left, separate
molecules right) with your threshold in the gap, and an exploratory molecule-count-vs-threshold curve.

**Anisotropic option** — weights x, y and z by their own measured spread, so one threshold fits all
three axes (the merge volume becomes an ellipsoid in nm, tighter along the better-resolved axis).
Affects the merge decision only; reported distances and sizes stay in true nm.

De-blinked molecules are labelled **De-blinked ID** in hover text, and their medians (not individual
TID medians) become the representative points for downstream analysis including Median Clustering.

> **Why not DBSCAN here?** Raw-point DBSCAN was the previous method and has been retired for merging:
> it merged whenever any two localizations touched, so a single stray point could fuse two molecules.
> Against ground-truth simulation its error count was 35 (splits + merges) versus 10 for centroid
> distance and 29 for majority vote.
>
> **Naming in the code.** The build originally had two DBSCAN levels, named `cluster` and `cluster2`
> to tell them apart. Those names are still used in the source even though the first level is now
> de-blinking and the second is Median Clustering — so `cluster` in the code means de-blinking, and
> `cluster2` means Median Clustering.

### Nearest Neighbor Analysis

Nearest-neighbour distances computed from median positions — de-blinked molecule medians when
de-blink correction is on, trace medians otherwise. Controls: *k*, optional distance labels (font
size and offset), a distance-range filter, and a downloadable table. NN lines are drawn on the plot.

When multiplexing is active a **scope** picker appears (within channel / between channels), and when
Median Clustering is active a **Within Median Cluster only** option restricts the search.

### Median Clustering

The app's one real clustering step. It runs on **representative points** — one point per trace median,
or per de-blinked molecule median when de-blinking is on — and groups them into larger spatial
structures, for example a collection of receptors arranged in a ring, or a row of nuclear pores along
the nuclear envelope.

Three grouping methods:

| Method | Parameters | Behaviour |
|--------|-----------|-----------|
| **DBSCAN** | `eps` (nm, default 200), `MinPts` | Fixed neighbourhood radius |
| **HDBSCAN** | `MinPts` only | Adaptive density; MinPts doubles as minimum cluster size |
| **OPTICS (Xi)** | `Xi` (0–1, default 0.05), `MinPts` | Adaptive density that keeps small structures; lower Xi = more, smaller clusters |

`MinPts` (default 15) applies to all three. Set it near the smallest number of traces you would call a
real structure — often 3–6, not higher.

Two optional cleanup steps:

- **Density prefilter** — before clustering, drop centroids with fewer than *n* neighbours within a
  radius, removing scattered background.
- **Prune sparse cluster members** — after clustering, demote to noise any member whose local spacing
  exceeds a multiple of its cluster's typical spacing (~2 trims clear outliers; lower is stricter).

**kNN-distance plot** — sorted distance to the (MinPts−1)-th neighbour of each centroid. The elbow
marks the structure/background split and therefore a good `eps`; the current `eps` is drawn as a red
line.

Settings can be saved and reset with the buttons at the bottom of the group.

#### Visualization

In **2D**: *Wrap Median Cluster in Convex Hull*.

In **3D**: *Wrap Cluster in Alpha Shape* (concave hull via `alphashape3d`), or *Show isosurface* — a
Gaussian-smoothed density surface with controls for grid resolution (12–40), padding (0–500 nm),
surface level (density quantile 0.50–0.98), opacity, and smoothness sigma (10–500 nm). The isosurface
can be restricted to one cluster ID and drawn in a single colour or per-cluster colours.

#### Fitting & Metrics (3D only)

Each cluster's representative points are fitted to a geometric model independently. Two fitting
methods:

| Method | Description |
|--------|-------------|
| **Heuristic** | Analytical: geometry derived directly from PCA axes and 5th–95th percentile extents. No iteration, no centre optimization. Fast and always converges. |
| **Weighted orthogonal distance** | Iterative: minimizes the sum of squared Euclidean distances from each representative point to the nearest point on the model surface. Equivalent to an orthogonal regression on the geometry. |

RMSD is reported in nm for both, so the two are directly comparable.

The **weighted orthogonal distance** method has additional options:

- *Weight by per-trace precision* — each representative point is weighted by $n / \sigma_\text{spread}^2$ (proportional to $1/\text{SEM}^2$ of the median position). Traces with a single localization (undefined spread) receive the dataset-median spread as a fallback.
- *Centre handling*:
  - **Fixed** — centre = convex-hull centroid; radius/semi-axes solved analytically (no iteration)
  - **Regularized (MAP)** — centre is free but penalized by a Gaussian prior anchored to the robust centre; strength $\lambda$ controls the trade-off. At $\lambda = 0.1$, a 5 nm centre shift adds ~2.5 nm² to the mean squared residual
  - **Free** — centre fully optimized within the data bounding box; use when the centre is expected to deviate significantly from the data centroid

**Fit shapes**:

| Shape | What is fitted | Key outputs |
|-------|---------------|-------------|
| **Sphere** | Sphere of radius R centred at the hull centroid | Centre (x, y, z), radius R, RMSD |
| **Ellipsoid** | Triaxial ellipsoid aligned to PCA axes | Centre, semi-axes a/b/c, axis orientations, RMSD |
| **Ring** | Flat circle (3D); PC3 = ring-plane normal | Centre, normal vector, radius, RMSD |
| **Arc** | Partial ring; angular extent auto-detected or set via hint | Centre, normal, radius, arc span (°), RMSD |
| **Helix** | Right-handed helix; PCA axis = helix axis | Centre, axis direction, helix radius, pitch (nm/turn), RMSD |

For **Ring**, **Arc** and **Helix**, *Render as tube* draws the fitted curve as a 3D cylinder of
user-defined radius (nm) for visual comparison against the data.

- *Arc span hint* (°) — initial angular extent for arc fitting (0 = auto-detected). The optimiser refines it.
- *Helix turns hint* — seeds the pitch estimate only; the fitted turn count follows from the data extent and the optimised pitch.

**Compute all-model RMSD comparison** (slower) fits every shape model to each cluster and reports all
RMSDs side-by-side, useful for model selection.

**Additional 3D overlays**: PCA axes as coloured arrows; per-cluster summary text (centre, radius,
RMSD) drawn on the plot; *Show distances to fit* (line segments from each representative point to the
nearest surface point, with a maximum line-length filter); and a fit-surface opacity control. The
metrics table is downloadable.

**Membrane protein options** (requires the global reference plane to be active) relate each cluster to
a global PCA-fitted reference surface, e.g. the plasma membrane:

- *Angle to plane* — angle between the cluster's primary PCA axis and the reference plane normal;
  0° = perpendicular to membrane, 90° = lying flat
- *Distance to plane* — signed distance from the cluster centroid to the reference plane
- *Bending radius* — curvature radius of an arc-like cluster:
  $$R_\text{bend} = \frac{\ell^2 + s^2}{2s}$$
  where $\ell$ is the PC1 half-extent (chord half-length) and $s$ is the PC2 extent (sagitta).
  **Only valid for arc-like clusters viewed roughly edge-on** — not meaningful for complete rings,
  spheres, or helices.

#### NN distances within Median Clusters

Nearest-neighbour distances between the representative points inside each cluster, drawn as lines:
*k* (default 3), a distance-range filter, optional numeric labels (font size and offset), and a CSV
export.

### Background image overlay

*(2D only.)* A confocal or widefield image (PNG / JPG / TIFF) placed underneath the localizations.
You set the **pixel size** and the image origin in data units, with ±10/±100/±1000 nudge buttons for
X and Y. Display controls: invert Y (origin at top-left), mirror horizontally, opacity, rotation
(0/90/180/270), brightness, contrast, and a **LUT** that recolours the image to match a fluorophore
(grayscale, green, magenta, yellow, cyan, red, blue, orange, or a custom colour). Overlay settings
save and load as JSON.

Two **auto-alignment** modes are available once an image and pixel size are provided:

**Foci mode**:
1. Bright blobs are detected in the confocal image by Gaussian smoothing followed by iterative peak finding with a minimum peak-separation constraint.
2. Dense clusters in the MINFLUX data are detected via a 2D density histogram with Gaussian smoothing.
3. Each possible assignment of MINFLUX clusters to confocal peaks is scored by the sum of squared residuals of the implied offsets (geometric consistency). The assignment with the lowest score and a physically valid bounding box is used.
4. If no assignment satisfies the bounds constraint, a rank-weighted Hough voting fallback is used.

**Boundary mode**:
1. The bright region of the confocal image is thresholded at the 50th percentile.
2. The intensity-weighted centroid of the bright region is aligned to the centroid of the MINFLUX data.
3. The convex hull of bright pixels is extracted as a cell boundary polygon in data coordinates.
4. Optionally, MINFLUX traces whose median position falls outside this polygon are removed.

*Show detection diagnostics* plots what the detector found, and a collapsible **Foci detection
settings** block exposes the sensitivity, edge mask, minimum separation, and the MINFLUX-side density
threshold / cell size / smoothing radius. Results are a starting point — nudge the origin to refine.

### Heatmap overlay

*(2D only.)* A KDE density heatmap layered over the localizations: grid resolution, opacity, palette,
smoothing (as a multiple of the Scott rule, or an absolute bandwidth in nm that overrides it),
bottom-percentile suppression, top-percentile clipping, gamma, log scale, and a colorbar toggle. By
default it uses trace medians; *Use raw points* switches to individual localizations.

### 2-color MINFLUX

Splits one acquisition into two channels using the per-trace **median DCR** (or the `fluo` column if
present, via *Override with 'fluo' groups*).

**Manual mode**: set two non-overlapping DCR ranges. Each trace is assigned to group 1 or 2 by its
median DCR; traces outside both ranges are unassigned.

**Auto-detect split (GMM)**: fits a two-component Gaussian mixture to the per-TID median DCR
distribution and solves for the crossing point of the two weighted densities. See
[DCR-Based 2-Color Separation](#dcr-based-2-color-separation).

**Per-trace split confidence** (collapsible) plots the posterior probability that each trace belongs
to the group it was assigned to. 1.0 = unambiguous; 0.5 = the molecule sits exactly on the split.
Traces below the cut are marked orange — these are the ones that would change channel under a small
threshold shift, and the ones that manufacture false short-range cross-channel pairs. Exportable so
you can exclude or inspect them.

Each group gets its own colour. Additional options:

- **Draw NN lines (2-color mode)** — cross-channel or within-channel nearest neighbours with *k*, a
  distance range, optional labels, and a CSV export.
- **2-color de-blink correction** — the same centroid/vote criteria as the mono-colour de-blinker,
  with its own threshold, *Measure threshold from data* button and gap diagnostic. De-blinked IDs
  carry their channel (`A-1`, `A-2` … / `B-1`, `B-2` …), and **only members of the same channel can
  merge**, so two colours can never be de-blinked into one molecule. In 3D the de-blinked molecules
  can be wrapped in alpha shapes.
- **2-color Median Clustering** — Median Clustering run per channel, with the same DBSCAN / HDBSCAN /
  OPTICS choice, prefilter and pruning, plus per-channel alpha shapes and isosurfaces in 3D.

Settings save and reset with the buttons at the bottom of the group.

> Numeric colour modes (z / time / density) override the 2-color colours. Use the default colour mode
> to see the groups.

### Multiplexed MINFLUX

Builds an **N-channel** dataset for exchange / sequential MINFLUX. Two channel sources:

- **Overlay loaded files** — pick several loaded files; each becomes one channel with its own colour.
- **Split active file by time (TIM)** — split one continuous acquisition. Either *contiguous
  breakpoints* (N breakpoints → N+1 channels, up to 7) or *explicit windows* (`start-end, start-end,
  …`), where localizations falling **outside** all windows — the wash / exchange steps — are removed.
  A trace-IDs-over-time histogram shows the split live: bar colour = channel, grey = removed, dashed
  lines = window edges.

Clustering runs strictly **per channel** — molecules imaged in different exchange rounds must never
merge. Channels feed both the plot and the Analysis tab.

### Estimate Ground-Truth Positions

> **Exploratory.** This inverts an assumed labelling model — read it as a hypothesis about the
> labelling, not a measurement.

Reverses the labelling model to propose where the underlying epitopes were. You choose a **labelling
preset** (the same presets as the Simulator: linker length, valency, dye spread, footprint, extra
blinks) or supply custom values, plus a labelling/detection efficiency.

- **Model from** — auto (de-blinked molecules if present, else raw traces), de-blinked molecules, or
  raw trace IDs. Auto is recommended: de-blinking collapses re-blinks and multivalent labels to ~1
  unit per molecule, turning an unstable divide-by-multiplicity count into a divide-by-~1 count.
- **Restrict to points inside Median Clusters** (on by default).
- **Label-complex merge radius** — units within this radius are pooled into one candidate GT, since a
  multivalent label scatters many fluorophores around one molecule. Blank = auto per preset (max of
  footprint and 2.5 × dye spread); 0 = no merge.
- **Stoichiometry-guided merge** *(model-driven)* — where the spatial merge leaves more GTs than the
  labelling stoichiometry predicts, collapse further toward that target. Only ever reduces the count;
  never splits. Sanity-check it against the two-count readout.
- **Displace GT inward** — the antibody binds a surface epitope and the linker throws the dye
  *outward*, so the trace cloud sits ~one linker-length outside the true structure. This pushes each
  GT back toward the local structure centre (perpendicular to a filament axis, toward the centroid
  for a compact or ring shape). Needs Median Clusters defined to know the structure.
- **Support-penalty strength** — how hard a GT backed by few trace IDs is penalized in the
  uncertainty. 0 = off; 1 = default (a 1-TID GT ≈ √2× less certain than a well-supported one).

Candidates are drawn as hollow circles coloured **green (low) → red (high) by position uncertainty**,
with optional contribution lines from each unit to its GT, an approximate count, and a crowding flag.

**Structure shape & size** classifies each cluster's shape from the point distribution (filament /
ring / punctum / sheet / diffuse) and reports neutral size descriptors: contour length, straight
extent, width, radius of gyration, and the spacing between estimated points. An optional assumed
subunit spacing divides a filament's contour length into an estimated subunit count.

### Global Reference Plane and Surface

*(3D only.)* Fits a reference plane or curved surface through the dataset, giving heights relative to
that surface rather than to the instrument axes. Useful when the sample is tilted or curved and raw z
is not the biologically meaningful coordinate.

- **Overlay XY Reference Plane (PCA based)** with a z offset and a downloadable per-point
  plane-distance table. It can be fitted from unclustered points only.
- **Surface method** — LOESS heightmap (with a span control), polynomial surface (degree 1–5), or
  thin-plate spline (with a smoothness λ; heavy). All share a grid resolution, opacity, colorscale and
  a *use points from* group selector.
- Specific TIDs can be excluded from the fit by ID.

---

## SMLM Analysis Tab

### Panels

The Analysis tab is organised into four sections. Each panel produces an interactive plot and, where
a table makes sense, an exportable CSV. Panels marked ⚠ carry an *exploratory* badge in the app —
they are research tools whose interpretation is not yet settled, not validated readouts.

**Dataset Groups** sits above everything: analyse **by individual file** (default) or **by named
group**. In group mode you define up to 6 groups (name, colour, member files) and the distribution
analyses relabel and pool by group so you can compare conditions. Files left unassigned are ignored.
Spatial statistics (Ripley / PCF) average per-file curves within a group rather than pooling
coordinates — pooling would mix observation windows and invent structure.

**Quality Control**

| Panel | What it answers |
|-------|-----------------|
| Localization Precision Cloud | Every localization minus its own trace centroid — the empirical precision σ_loc, with 1σ/2σ covariance ellipses whose tilt exposes drift |
| Fiducial / Contaminant Scan | Flags traces behaving like gold beads or aggregates rather than single molecules (report only) |
| Trace Detection over Time | Detection rate across the acquisition — photobleaching and drop-off |

**Nearest-Neighbour & Spatial Statistics**

| Panel | What it answers |
|-------|-----------------|
| Batch NN Distance Histogram | k-th nearest-neighbour distances across files or channels |
| Ripley's L(r) | Clustering vs randomness across length scales, with an optional Monte-Carlo CSR envelope |

**Median-Cluster Analyses**

| Panel | What it answers |
|-------|-----------------|
| Trace IDs / De-blinked Molecules per Median Cluster | Occupancy — how many traces, or how many de-blinked molecules, per structure |
| Median Clusters Detection over Time ⚠ | When structures appear during acquisition |
| Median Clusters Crowding vs Detection-Time Shift ⚠ | Whether crowding biases when a structure is detected |
| Within Median Clusters NN vs Cluster Size ⚠ | Internal spacing as a function of structure size |
| Cluster Viewer and Shape Fitting ⚠ | Per-structure 3D inspection and geometric fits |
| Radius of Gyration (RoG) | Structure size from the second moment, with anisotropy |

**Two-Color Analysis**

| Panel | What it answers |
|-------|-----------------|
| Cross-Channel NN & Co-localization | Two-colour proximity, as distances and as a co-localization curve |
| Bivariate Ripley's L₁₂(r) | Cross-channel clustering across scales |
| Pair Cross-Correlation g₁₂(r) | Cross-channel correlation as a function of separation |
| 2D KDE Density Map | Smoothed density field per channel, with a pixel-level overlap metric |

All analyses operate on the filtered, de-blinked and clustered data visible in the Display tab. In
2-color mode, group membership comes from the DCR split (or the `fluo` column).

**NN unit.** The Batch NN and Cross-Channel panels share a *unit* selector that decides what counts as
one point: TID centroids, de-blinked molecule centroids, or either of those restricted to within each
Median Cluster. Raw localizations are deliberately not an option — they are not molecules.

### Two-colour method details

The mono-colour statistics (NN, Ripley's L, RoG) are described under
[Key Calculations](#key-calculations). The four cross-channel panels:

#### Cross-Channel Nearest Neighbor & Co-localization Curve

Directional nearest-neighbour distance from each point in one group to its nearest point in the other,
computed in **2D (x, y)** with `RANN::nn2`:

$$d_i^{G_1 \to G_2} = \min_{j \in G_2} \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}, \quad i \in G_1$$

Both directions are computed independently. The **co-localization curve** is the empirical CDF:

$$F(r) = \frac{1}{|G_1|}\; \#\!\left\{i \in G_1 : d_i^{G_1 \to G_2} \leq r\right\}$$

— the fraction of group-1 points with a group-2 neighbour within *r*. A steep rise at small *r*
indicates strong co-localization; a gradual rise suggests coincidental proximity.

#### Bivariate Ripley's L₁₂(r)

Tests whether two populations co-localize beyond chance given their individual densities.

$$\hat{K}_{12}(r) = \frac{A}{n_1 n_2} \sum_{i \in G_1}\sum_{j \in G_2} \mathbf{1}\!\left(\|\mathbf{p}_i - \mathbf{p}_j\| \leq r\right)$$

with the same L-transform as the univariate case. The plotted $\hat{L}_{12}(r) - r$: > 0 = co-localize
more than two independent CSR processes, ≈ 0 = spatial independence, < 0 = exclusion.

When Median Clustering is active, $\hat{L}_{12}$ is computed **per cluster** and averaged, giving a
within-structure measure rather than a global one. The 2000-point subsampling limit applies per group.

#### Pair Cross-Correlation g₁₂(r)

The annulus-normalized derivative of $K_{12}$ — more sensitive to the *scale* of co-localization than
the cumulative $L_{12}$:

$$g_{12}(r) = \frac{A}{n_1 n_2} \cdot \frac{\#\text{cross-pairs in annulus}\bigl(r-\tfrac{\delta}{2},\, r+\tfrac{\delta}{2}\bigr)}{\pi\!\left[(r+\tfrac{\delta}{2})^2 - (r-\tfrac{\delta}{2})^2\right]}$$

with the 3D form using a spherical shell and $\tfrac{4}{3}\pi[(r+\tfrac{\delta}{2})^3 - (r-\tfrac{\delta}{2})^3]$.
Under CSR with independent groups, $g_{12}(r) = 1$; > 1 indicates co-localization at scale *r*, < 1
exclusion. $L_{12}(r)$ accumulates all pairs up to *r* and can mask scale-specific effects — the PCF
pinpoints them.

#### 2D KDE Density Map & Co-localization

KDE maps are computed for each group on a shared grid (`MASS::kde2d`) and rendered as heatmaps.
Bandwidth follows Silverman's rule of thumb:

$$h = 0.9 \cdot \min\!\left(\hat{\sigma},\; \frac{\text{IQR}}{1.34}\right) \cdot n^{-1/5}$$

computed on the combined x-coordinates of both groups and applied identically in x and y, so the two
channels share a grid. Co-localization is the Pearson *r* between the two density surfaces.

> **Caveat**: the Pearson *r* between KDE surfaces is strongly bandwidth-dependent. A large bandwidth
> smooths both channels and inflates the correlation even for spatially distinct populations. Use it
> for relative comparisons across conditions, not as an absolute co-localization measure.

### Export

- **Export MINFLUX data (.csv)** — the currently filtered localizations
- **Download Filter Settings** — all filter parameters as JSON, reloadable
- **NN distance tables** — CSV export of nearest-neighbour results (mono-colour, 2-colour, and within-cluster)
- **Median Cluster metrics table** — per-cluster morphology and fit results
- **Per-trace split confidence** — 2-colour posterior probabilities

---

## Tracking Display Tab

Sidebar groups appear in this order.

### Data Upload

Accepts `.npy`, `.json`, `.csv` (pyMINFLUX or this app's export), and the fast cache formats
`.qs` / `.fst` / `.rds`. Multiple files merge into one dataset with track IDs prefixed `F1_`, `F2_`
so IDs from different files never collide.

Large files are **chunked** on upload — split into per-TID groups written to a manifest directory —
so playback stays responsive. Chunk size is set by *TIDs per chunk* and only affects Playback mode;
Full-tracks mode always renders everything. The filter pipeline always reads **every** chunk, so
filtering is never restricted to the chunk you are viewing.

**No filter is applied automatically on upload** — CFR in particular is left alone, since it is not a
strong artefact discriminator for tracking data. Filter it yourself via the CFR range below if you
want to.

### Axial (z) calibration

Same control as the SMLM tab: an absolute refractive-index factor (default 0.7), applied once after
the file is read raw. Set both tabs to the same value when comparing axial distances between them.
Greyed out when the dataset has no axial variation.

### Acquisition Parameters

Optional hardware model for cycle time, a manual cycle-time override, and a manual σ_loc override.
Cycle time drives gap-splitting and the MSD lag axis, so getting it right matters more here than
almost anywhere else — the status line shows the value in use and where it came from.

### Track Cleaning & Filtering

Cleaning runs on the raw acquisition timing; the metadata, length and duration filters then apply per
segment:

```
Raw data
 └── [1] Split at time gaps → Trim degraded ends → Continuity gate
      └── [2] Metadata filters  (EFO, CFR, DCR, TIM, ECO, dwell, fluo)
           └── [3] TID length filter  (number of points per trace)
                └── [4] Trace duration filter  (max(tim) − min(tim), ms)
```

The three cleaning stages run **before** any filter, in this order. Each is optional and each reports
how many localizations it removed.

| Stage | Controls | What it does |
|-------|----------|--------------|
| **Split at time gaps** | `Split tracks at time gaps`, `Adaptive threshold`, `Gap threshold (cycles)`, `Min localizations per segment` | A MINFLUX track ID can span a long dark interval and rejoin — which is not one continuous trajectory. Splits at intervals longer than the threshold and discards segments too short to analyse. Adaptive mode derives the threshold from the track's own interval distribution. |
| **Trim degraded ends** | `Trim dark-gap-separated ends`, `Dark-gap threshold (cycles)`, `Trim dim (low-EFO) ends`, `EFO threshold (× track median)` | Bleaching tails: the end of a track is often sparse and dim before it dies. Trims trailing localizations separated by a long dark gap, or whose EFO has fallen below a fraction of the track median. |
| **Continuity gate** | `Track-continuity gate`, `Dark-cycle tolerance (half-life, cycles)`, `Min continuity score` | Drops whole tracks that are mostly dark cycles. The score weights each localization by how many cycles were skipped before it, with an exponential half-life, so a track that only nominally exists across a long window is rejected. |

Then the metadata filters: pick a variable (EFO, CFR, DCR, ECO, track length, duration…), see its
histogram with the unfiltered data behind it in grey, and set a range. **Crop to view** restricts to
a drawn box, keeping only tracks with at least a set percentage of their localizations inside.

Two things to know about how filters apply:

- **Filters do not live-update the plot.** They apply on **Plot Data** (or Filter & Rechunk). The
  histogram and the read-outs update live so you can see what you are selecting; the rendered tracks
  stay put until you ask. This is deliberate — re-rendering a WebGL track view on every slider move is
  slow and makes the filters feel laggy.
- **Filter & Rechunk** makes the current filtered set the working dataset and re-chunks it, so
  subsequent operations are faster. *Reset ranges* restores the original manifest.

Filter settings save and load as `.json`, and the filtered set exports as `.csv` (carrying
`z_scale_applied` for exact re-import).

### Playback Controls, Display Options, Overlays

**Playback** animates the tracks through time with speed, trail length and marker controls, plus
*Auto-tune* to fit the frame budget. **Display Options** carries the 2D/3D toggle, the XY/XZ/YZ
viewing plane, marker size and the colour mode (track ID, dark cycles, or speed). A **Background
Image Overlay** (same controls and foci detection as the SMLM tab) and a **Dwell-time Heatmap** can be
layered underneath.

---

## Tracking Analysis Tab

Panels in UI order. Several depend on σ_loc and the cycle time, so run **Data Quality & Calibration**
first.

### Data Quality & Calibration

The foundation panel. Estimates **σ_loc** (localization precision) and **D** by the CVE method
(Vestergaard 2014), using only single-cycle steps. Both come from two moments of the displacement
series — ⟨|Δr|²⟩ and the covariance of consecutive steps — so no MSD fit is needed and the two are
separated cleanly. The covariance is negative for localization noise, which is what makes the
separation possible.

Reports cycle time, free D, σ_loc and the step count. The D shown is corrected for the **motion
degrees of freedom** set below, and says so when the correction is active.

### Time Signature Analysis (TSA)

Examines the temporal structure of the acquisition: the inter-localization interval (Δt) distribution,
integrated / dark-cycle statistics, and per-track cycle progression. Useful for confirming the cycle
time and for seeing whether gap-splitting is doing the right thing before you trust anything
downstream.

### Sample geometry

Sets the **motion degrees of freedom**: volume (3), surface (2), or filament (1). This is not
cosmetic. A per-axis estimator spreads 2-DOF motion over 3 recorded axes, so a membrane protein
tracked in 3D reports D at 2/3 of its true value unless the geometry says otherwise. It also sets the
χ² degrees of freedom for the confinement band. **Defaults to volume** — set it to surface for
membrane work or every D and confinement call below inherits the bias.

### MSD (diffusion)

Mean-squared displacement, binned by **integer cycle-count lag** rather than wall-clock time
(MINFLUX cycle-skipping makes Δt non-uniform). Fits free, confined and anomalous models. The linear
fit uses Michalet & Berglund optimal-point selection rather than a fixed lag range; the confined and
anomalous fits deliberately use the whole window, because that optimality theory is derived for free
diffusion only. A **cycle-binning diagnostic** shows how steps distribute across cycle counts.

### Confinement (MSD index)

Confinement ratio from the MSD curve's departure from linearity, with a simulated null for
comparison. Confined percentages are **not comparable across different track lengths** — a longer
track has more chance to reveal a plateau.

### Track Segmentation — Change-Point + GMM

Splits tracks at statistically detected change points, then clusters the resulting segments by their
diffusive behaviour into a Gaussian mixture. Use when a single D per track is clearly the wrong model.

### Jump Distance Distribution (JDD)

Fits the distribution of single-step displacements as an N-component mixture — **Rayleigh in 2D,
Maxwell in 3D** — with σ_loc bias correction and an integrated-steps-only mode. Complementary to MSD:
it uses the full displacement distribution rather than its second moment, so it separates populations
that share a mean. Steps must be pooled at a fixed real Δt; mixing cycle counts pollutes the fit.

### Displacement vs Diffusion

Compares each track's net displacement against the free-diffusion expectation √(4DT) for its own
duration, with a Rayleigh null. Tracks far above the band are candidates for directed transport;
far below, for confinement. Uses the motion-DOF-corrected D and the χ² degrees of freedom from
**Sample geometry**.

### Directionality / Autocorrelation

Velocity autocorrelation and turning-angle statistics — persistent directional motion shows positive
correlation at short lags where free diffusion shows none.

### Diffusion State Segmentation (HMM) and Per-Segment MSD

A hidden Markov model over the displacement series assigns each step to a diffusive state, giving
per-state D and dwell times, and the transition structure between them. **Per-Segment MSD** then
re-runs the MSD independently within each HMM state, so a two-state track yields two curves instead
of one average that describes neither.

### Speed & Transport

Drift velocity and a noise-corrected local α(t). Path speed alone cannot distinguish directed
transport from fast diffusion — a random walk sampled finely has a high path speed — so the drift
term and the local anomalous exponent are the informative quantities here. An **Advanced** block
exposes candidate directed segments with per-event corroboration.

### Spatially Resolved D-map

Bins the field of view and estimates D per bin, showing where in the sample diffusion differs.
Subtracts the σ_loc noise floor so slow regions are not reported as diffusing.

### Inspect a single track

Per-track inspection: raw coordinates, intervals, per-step displacements and the fits for one
selected track. The first place to look when an aggregate result seems wrong.

---

## MINFLUX Simulator Tab

Two phases, each its own sub-tab.

### Phase 1 — SMLM (synthetic localizations)

Generates labelled synthetic structures in 2D or 3D. Two **simulation modes**:

- **Normal** — visualize one structure with one labelling strategy.
- **Acquisition** — lay out a field of view, place structures in it, and model the whole acquisition
  (detection rate, duration, background). Adds the *Acquisition timeline* tab.

**Ground-truth Structure**:

| Shape | Description |
|-------|-------------|
| Cluster | Round cluster with uniform or Gaussian radial distribution |
| Ring | Circular arrangement with controllable band width and Z thickness |
| Arc | Partial ring spanning a user-defined angular range |
| Double ring | Two parallel rings separated along Z (e.g. nuclear pore complex) |
| Ellipse | Filled or hollow ellipse/ellipsoid (semi-axes a, b, c) |
| Hollow sphere | Shell with adjustable thickness |
| Filament | Linear arrangement with optional lateral jitter and helical twist |
| Oligomeric ring stack | Stacked oligomeric rings |
| Image (confocal) | Places molecules proportionally to local pixel intensity of a loaded confocal image |

A fixed ground-truth seed lets you compare labelling strategies on the *same* structure.

**Labeling Strategy** — each labelled molecule generates N fluorophore positions, and each fluorophore
becomes one trace ID in the output. Presets cover direct tags (SNAP / Halo / Click, DOL = 1),
nanobodies, primary monoclonal and polyclonal IgG (with and without NHS multi-labelling), primary +
secondary (nanobody or IgG) sandwiches, streptavidin, and a Custom option. These are the same presets
the SMLM tab's *Estimate Ground-Truth Positions* inverts.

**Localization precision** — either a manual σ (separate XY and Z) or a **CRLB from photon budget**.

**Acquisition** *(acquisition mode)* — acquisition time in hours/minutes, structure and background
detection rates, and a **Calibrate from a real run** loader: point it at a real `.npy` / `.json` /
`.csv` and it fills the FOV, duration and detection-rate profile from that data, then overlays the
measured rate on the Acquisition plot. A *structure fraction of detections* slider (auto-estimated
from spatial density) splits the rate into structure and background components.

**Main-panel tabs**: *MINFLUX* (the synthetic localizations), *Confocal view* (a diffraction-limited
render of the same ground truth — objective NA, pinhole in Airy units, pixel size, photons per
fluorophore, background, Poisson shot noise), *Acquisition* (the timeline plot and summary), and
*Compare strategies* (two labelling strategies A and B side by side on identical ground truth).

Simulated data can be sent directly to the SMLM or Tracking tabs for immediate analysis.

### Phase 2 — Diffusion / Tracking (ground-truth simulator)

> **Exploratory validation tool.** Its purpose is to check that the Tracking tab's estimators recover
> a *known* D.

Free Brownian motion on a cell shape, sampled the MINFLUX way — cycle clock, dark cycles,
heteroscedastic σ_loc and photobleaching — producing MINFLUX-format tracks with a ground-truth D.
Download the tracks and load them in the Tracking tab to check CVE / MSD / confinement.

- **Cell geometry**: sphere (coccus), rod / spherocylinder (E. coli, with a choice of long axis),
  ellipsoid (ovoid / coccobacillus), membrane plane (a tilted flat sheet in 3D, no curvature → no
  projection bias), or free (unbounded 3D — the simplest null).
- **Diffusion**: surface (2 DOF, membrane/wall protein) or volume (3 DOF, reflecting off the wall).
- **Motion type**: free, confined (OU tether → MSD plateau at ~2·L_c²), directed (drift velocity along
  the rod long axis), or two-state (mobile ↔ immobile, ground truth for HMM state recovery).
- **Survivorship** — optionally lose the molecule when the *true* displacement across a dark gap
  exceeds a re-acquisition radius. This is the dark-cycle bias that manufactures apparent confinement.
- **Sampling model**: cycle time, dark-cycle probability, mean photons per localization, and the
  effective (final-iteration) TCP diameter L, from which σ_loc = L / (2·√N_photons).
- **Calibrate from real MINFLUX data** — load a real file to auto-fill cycle time, photons, dark-cycle
  probability and median track length. Tip: apply your filters in the Tracking tab first, export the
  filtered CSV, and calibrate on that, so the model matches exactly the tracks you analyse.

Both the tracks and the ground truth export as CSV.

---

## Key Calculations

### Unit Conversions & Import

Abberior files store coordinates in **meters**. On load, the app converts to nanometers:

$$x_\text{nm} = x_\text{m} \times 10^9$$

**Z-axis scaling** is a user-set calibration, not a fixed constant. Files are read with **raw z**;
the factor you set (SMLM sidebar *Axial (z) calibration*, Tracking sidebar *Axial (z) scaling factor*,
both defaulting to 0.7) is applied once afterwards, on Plot Data. CSVs are first returned to raw z —
the app detects what factor they already carry, exactly when the file was exported from here (a
`z_scale_applied` column) and by assumption otherwise. Exports write that column so a round-trip is
exact. Because the factor multiplies every axial quantity — z spread, cluster extent, 3D neighbour
distances, Rg — applying it twice would compress z by 30%, which the detection prevents.

Scaling is applied after unit conversion, to correct for the refractive-index mismatch between the
sample medium and the objective immersion medium:

$$z_\text{corrected} = z_\text{nm} \times s_z \quad (s_z = 0.7 \text{ by default})$$

**Multi-iteration MINFLUX**: each TID produces localizations at several iteration steps (`itr`). By
default only the final iteration is retained — chosen by the `fnl` flag for `.npy`, by
`itr == max(itr)` for `.json`. The CFR value used is the last non-zero CFR recorded across all
iterations for that TID (*last_nonzero* carry-forward), because CFR is often only populated in early
iterations.

**Valid-only filter**: records with `vld = FALSE` are discarded on load.

---

### Dwell Time

The dwell time (number of MINFLUX measurement steps per localization) is derived from the photon-count columns:

$$\text{dwell} = \mathrm{round}\!\left(\frac{\text{eco}}{\text{efo} / 1000} \cdot \frac{1}{T_\text{step}}\right)$$

where `eco` is the total estimated photon count, `efo` is the emission rate in kHz, and $T_\text{step}$ is the user-specified dwell time per step in ms (default 1 ms). The ratio eco / (efo / 1000) gives the effective dwell duration in ms; dividing by $T_\text{step}$ converts to step counts.

---

### Trace Spread (SD Filter)

For each TID, the per-axis standard deviation of all localizations is computed:

$$\sigma_x(\text{tid}) = \mathrm{SD}(\{x_i : \text{tid}_i = \text{tid}\})$$

and analogously for $\sigma_y$ and $\sigma_z$.

The SD filter ranges set bounds $[\sigma_\text{min},\, \sigma_\text{max}]$; a TID is kept only if its
$\sigma$ falls within those bounds on every active axis. Setting the upper bound to a small value
retains only spatially tight traces (useful for removing freely diffusing molecules), while a broad
range passes all traces. Traces below the *Min localizations per trace* setting are excluded from the
statistic entirely.

The **spread panel** histogram shows the distribution of per-TID spreads and overlays the current
filter range. The median spread is reported as median(σ̃) ± nMAD, with nMAD = 1.4826 × MAD — the
1.4826 makes nMAD a consistent estimator of σ for normal data (same convention as pyMINFLUX).

---

### DCR-Based 2-Color Separation

**Per-TID median DCR** is used as a robust summary to classify each trace into one of two spectral groups:

$$\widetilde{\text{dcr}}_\text{tid} = \mathrm{median}(\{\text{dcr}_i : \text{tid}_i = \text{tid}\})$$

**Manual mode**: the user defines two ranges $[d_{1,\text{lo}},\, d_{1,\text{hi}})$ and $[d_{2,\text{lo}},\, d_{2,\text{hi}}]$. A TID is assigned to group 1 or 2 if its median DCR falls in the respective range; TIDs outside both ranges are unassigned.

**Auto-detect mode (GMM)**: fits a two-component Gaussian mixture model to the distribution of $\{\widetilde{\text{dcr}}_\text{tid}\}$ using the `mclust` package:

$$p(\text{dcr}) = \pi_1 \,\mathcal{N}(\mu_1, \sigma_1^2) + \pi_2 \,\mathcal{N}(\mu_2, \sigma_2^2)$$

Components are sorted so $\mu_1 < \mu_2$ (group 1 = lower DCR). The split threshold $d^*$ is the crossing point of the two weighted component densities:

$$\pi_1 \,\mathcal{N}(d^*;\mu_1, \sigma_1^2) = \pi_2 \,\mathcal{N}(d^*;\mu_2, \sigma_2^2)$$

solved numerically in $(\mu_1, \mu_2)$ with `uniroot`. If the two component means differ by less than
0.05 DCR units, the distribution is considered unimodal and auto-detection is aborted with a warning.
The same fit supplies the **per-trace posterior** shown in the split-confidence panel.

---

### Spatial Density

For each point (localization or TID median), the local density is the count of other points within a sphere (or disk in 2D) of radius $r$:

$$\rho_i = \#\{j \neq i : \|p_i - p_j\| \leq r\}$$

Computed via `dbscan::frNN` (fixed-radius nearest-neighbor search), which uses a kd-tree for $O(n \log n)$ performance. The radius $r$ is user-controlled (default 500 nm in the colour-by-density mode, 150 nm elsewhere).

---

### Nearest-Neighbour Distances

The $k$-th NN distance for point $i$ is:

$$d_i^{(k)} = k\text{-th smallest} \;\|\mathbf{p}_i - \mathbf{p}_j\|, \quad j \neq i$$

computed with a kd-tree (`RANN::nn2` / `dbscan::kNN`). Histograms are summarized by the median and by
either the IQR or the MAD as a spread indicator.

> The first-NN distance distribution reveals the typical inter-molecule or inter-structure spacing. A
> narrow peak indicates a regular, lattice-like arrangement; a broad tail suggests heterogeneity.
> **A sharp peak below ~40 nm is usually the label footprint, not biology** — the distance-range
> filter trims the histogram but does not correct that bias. NN carries no CSR reference; use
> Ripley's L(r) with its envelope if you need significance.

---

### Ripley's L(r)

A spatial clustering statistic that tests whether points are more aggregated than expected under
complete spatial randomness (CSR). MINFLUXR uses a **naive (boundary-uncorrected) estimator**:

$$\hat{K}(r) = \frac{A \cdot 2}{n(n-1)} \sum_{i < j} \mathbf{1}(\|\mathbf{p}_i - \mathbf{p}_j\| \leq r)$$

where $A$ is the bounding area (or volume $V$ in 3D) and $n$ is the number of points. The
L-transformation:

$$\hat{L}(r) = \sqrt{\frac{\hat{K}(r)}{\pi}} \quad \text{(2D)}, \qquad \hat{L}(r) = \left(\frac{3\,\hat{K}(r)}{4\pi}\right)^{1/3} \quad \text{(3D)}$$

The plotted quantity is $\hat{L}(r) - r$:

- $> 0$: clustering at scale $r$ (more pairs than expected under CSR)
- $\approx 0$: complete spatial randomness
- $< 0$: regularity / inhibition

**Implementation notes**:

- Runs on **de-blinked units** (trace or de-blinked-molecule centroids), so repeated appearances of one
  molecule cannot manufacture a short-*r* clustering peak.
- $r_\text{max}$ defaults to one quarter of the shorter data span, keeping the estimator in the
  statistically reliable range.
- If $n > 2000$, a random subsample of 2000 points is used and a warning is shown (pairwise distances
  are $O(n^2)$).
- An optional **Monte-Carlo CSR envelope** simulates *n* random points in the *identical frozen
  window*, so the band and the observed curve carry the same edge bias and it cancels.

> There is no analytic edge correction: near the top of the *r* range, points close to the boundary
> lose neighbours and L−r bends downward. Read significance as the observed curve **leaving the
> envelope**, not as departure from the flat zero line. Without the envelope, read shape and peak
> position rather than magnitude.

---

### Radius of Gyration

The radius of gyration quantifies the spatial compactness of a Median Cluster. **RoG is computed from
representative point positions** (per-TID medians, or de-blinked molecule medians when de-blinking is
on) rather than from individual localizations, so each molecule contributes one point and localization
noise is suppressed.

For a cluster of $m$ points at $(\tilde{x}_k, \tilde{y}_k, \tilde{z}_k)$, the centroid is the mean of
those positions and:

$$R_g = \sqrt{\frac{1}{m}\sum_{k=1}^{m}\left[(\tilde{x}_k - \bar{x})^2 + (\tilde{y}_k - \bar{y})^2 + (\tilde{z}_k - \bar{z})^2\right]}$$

Only clusters with $m \geq 2$ are included (a single point has $R_g = 0$ by definition). The
exportable table reports cluster ID, point count, centroid coordinates and $R_g$. RoG is
small-*N*-biased — it is reported with its anisotropy and should not be used as a classifier on its
own.

---

### MSD Computation

MINFLUX tracking data is irregularly sampled in time, because cycles are skipped when the molecule is
dark. The Tracking tab therefore bins MSD by **integer cycle-count lag**, not wall-clock time.

For a single track $\{(x_i, y_i, z_i, t_i)\}$ sorted by time, all pairs $(i, j)$ with $i < j$ are
considered, and each pair contributes its squared displacement to the bin whose centre $\tau$
satisfies $\tau_\text{lo} \leq t_j - t_i < \tau_\text{hi}$:

$$\text{MSD}(\tau) = \frac{1}{N_\tau} \sum_{\substack{(i,j)\,:\,\\ \tau_\text{lo} \leq t_j - t_i < \tau_\text{hi}}} \left[(x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2\right]$$

where $N_\tau$ is the number of pairs in that bin.

The **ensemble MSD** across all tracks is computed bin-by-bin as either:
- the **pair-weighted mean**, $\overline{\text{MSD}}(\tau) = \sum_k N_{\tau,k} \cdot \text{MSD}_k(\tau) \,/\, \sum_k N_{\tau,k}$ (downweights bins with few observations), or
- the **track-unweighted mean**, a simple average of per-track MSD values.

Bins where fewer than a user-specified fraction of tracks contribute are masked as NA.

---

### Diffusion Coefficient Fitting

**Brownian diffusion** model:

$$\text{MSD}(\tau) = 2n D \tau + 2n \sigma_\text{loc}^2$$

where $n$ is the number of dimensions, $D$ the diffusion coefficient (nm²/ms = µm²/s × 10⁻³), and
$\sigma_\text{loc}$ the localization precision (nm). The intercept $2n\sigma_\text{loc}^2$ arises from
localization noise and does not reflect true motion. Note that $n$ here is the **motion degrees of
freedom** set in *Sample geometry*, not the number of recorded axes.

Fitting is two-pass weighted least squares:
1. Unweighted linear regression for initial estimates $D_0$, $\sigma_0$
2. Weighted regression with $w_k \propto N_{\tau,k} / \widehat{\text{MSD}}(\tau_k)^2$, which minimizes
   variance for Poisson-like pair-count noise

The linear fit range comes from **Michalet & Berglund optimal-point selection** rather than a fixed
percentage of the curve.

**Anomalous diffusion** model:

$$\text{MSD}(\tau) = 2n D_\alpha \tau^\alpha + 2n \sigma_\text{loc}^2$$

fitted by non-linear least squares (Levenberg–Marquardt via `minpack.lm`). The anomalous exponent:
- $\alpha < 1$: subdiffusion (confined or obstructed)
- $\alpha = 1$: normal Brownian diffusion
- $\alpha > 1$: superdiffusion (directed transport)

The confined and anomalous fits use the whole lag window, because the optimal-point theory is derived
for free diffusion only. Reported units: $D$ in µm²/s, $\sigma_\text{loc}$ in nm.

---

### CVE (Covariance-Based Estimator)

The Data Quality panel estimates $D$ and $\sigma_\text{loc}$ without fitting an MSD curve, following
Vestergaard et al. 2014. Using only single-cycle steps, both follow from two moments of the
displacement series — the mean squared step $\langle|\Delta r|^2\rangle$ and the covariance of
*consecutive* steps. Localization noise makes that covariance **negative**, which is exactly what
separates it from real motion. This is the app's primary σ_loc and free-D readout; the MSD panel is
the cross-check, not the other way around.

---

### Confinement Index

The confinement index (CI) provides a model-free measure of diffusion character:

$$\text{CI} = \frac{\text{MSD}(\tau_\text{late}) / \text{MSD}(\tau_\text{early})}{\tau_\text{late} / \tau_\text{early}}$$

where $\tau_\text{early}$ and $\tau_\text{late}$ are the medians of the first and last 15% of lag-time
bins (minimum 2 bins each). For free Brownian diffusion MSD ∝ τ, so CI = 1:

- **CI ≪ 1**: MSD plateaus → confined diffusion
- **CI ≈ 1**: linear growth → Brownian
- **CI ≫ 1**: superlinear growth → directed motion

Per-track fits produce a table of ($D$, $\alpha$, $\sigma_\text{loc}$, CI). **MSD clustering** applies
k-means in the (log₁₀ $D$, log₁₀ CI) feature space (z-scored, with log₁₀ $D$ upweighted ×2 as the
primary discriminator) to identify mobility sub-populations, ordered from most confined (cluster 1) to
most mobile.

> **Track-length dependence**: CI, and therefore the "confined %" readout, is **not comparable across
> datasets with different track lengths** — a longer track has more opportunity to reveal a plateau.
> Compare only at matched track lengths, and read the simulated null band shown alongside.

---

## Plugin System

Additional analysis modules can be added as plugins without modifying the base app:

```
MINFLUXR/
└── plugins/
    └── MyPlugin/
        └── plugin.R        ← entry point; calls register_plugin(tabs, server_fn)
```

`global.R` auto-discovers and loads every `plugins/*/plugin.R` at startup. Each plugin receives `.plugin_dir` (the path to its own folder) for loading relative resources. Plugin tabs are appended to the navbar after the three built-in tabs.

---

## Further Documentation

| File | Contents |
|------|----------|
| `docs/MINFLUX_analysis_methods.md` | Part I: the implementation and maths of every analysis, taken from the source. Part II: a literature review of whether each method is appropriate for MINFLUX data, with GREEN / AMBER / RED verdicts |
| `docs/simulation_module.md` | The ground-truth diffusion simulator: model, validation and CVE round-trip |
| `docs/arch/` | Architecture overview |

---

## Notes & Tips

- **No automatic quality filter**: nothing is removed on load except invalid records (`vld = FALSE`).
  In particular there is **no automatic CFR cutoff** in either tab — CFR is not a reliable artefact
  discriminator for tracking data, and in SMLM it is your call. Set the CFR range yourself if you want
  one.
- **Z-axis scaling**: the default of 0.7 corrects for the refractive-index mismatch between water
  (n ≈ 1.33) and glass (n ≈ 1.52) in a typical sample. Adjust it if your mounting medium differs, and
  keep the SMLM and Tracking tabs on the same value when comparing between them.
- **CFR carry-forward**: because CFR is measured during early MINFLUX iterations and may be NaN in the
  final one, the app takes the last non-zero CFR across all iterations. This keeps CFR filtering
  working even though only the final iteration is retained.
- **De-blink before you count**: a single fluorophore produces several trace IDs across sequential
  blinking events, so a raw TID count is a count of **blinking events**, not molecules. Enable de-blink
  correction before interpreting per-cluster occupancy. Experiments with irreversible photoactivation
  (e.g. PALM) are less affected.
- **Multiple datasets**: in the SMLM tab several files can be loaded simultaneously and compared with
  the Active file selector. Lock filters when switching to prevent accidental resets, and use Dataset
  Groups in the Analysis tab to pool and compare conditions.
- **Filter & Rechunk is non-destructive**: the original chunk set is frozen at upload time. Every
  Filter & Rechunk reads from that frozen original, so resetting all filters and rechunking recovers
  the full dataset.
- **MSD fit range**: the linear fit uses Michalet & Berglund optimal-point selection for you. If you
  override it, keep to the early lags — later points rest on fewer pairs, have higher variance, and
  deviate from the linear model through confinement or undersampling.
