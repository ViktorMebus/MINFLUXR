# MINFLUXR user guide

Panel-by-panel reference for the app: what each control does, what each number means, and the
calculations behind them. Start at [the README](../README.md) for installation and a first run.

For the methods themselves (the maths, the literature and how appropriate each method is for MINFLUX
data) see [MINFLUX_analysis_methods.md](MINFLUX_analysis_methods.md).

- [Loading Data](#loading-data)
- [MINFLUX Data Columns Reference](#minflux-data-columns-reference)
- [SMLM Display Tab](#smlm-display-tab)
- [SMLM Analysis Tab](#smlm-analysis-tab)
- [Tracking Display Tab](#tracking-display-tab)
- [Tracking Analysis Tab](#tracking-analysis-tab)
- [MINFLUX Simulator Tab](#minflux-simulator-tab)
- [Key Calculations](#key-calculations)
- [Plugin System](#plugin-system)
- [Notes & Tips](#notes--tips)
- [References](#references)

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

`.npy` is read by a **hand-written pure-R NumPy decoder** (`read_npy()` in `helpers_shared.R`) - there
is no Python or `reticulate` dependency. It decodes the header, the structured dtype descriptor and
the raw data block directly.

Two constraints follow from that:

- **Only the newer, flat "tabular" MINFLUX format is supported.** The **legacy nested `m2205` dtype is
  not read** - the app stops with a clear error asking you to export the tabular format instead.
  Nested structured dtypes have no flat column layout to decode.
- Fortran-ordered arrays are rejected (re-save C-ordered). A plain unnamed array is rejected too: the
  reader needs the named fields (`loc`/`lnc`, `tid`, `tim`, `efo`, `eco`, `cfr`, `dcr`, `vld`, `fnl`, …).

`.npy` and `.json` converge on the same post-processing: `vld = FALSE` records are dropped, iterations
are collapsed to one row per localization, and CFR is carried forward. The one difference is how the
final iteration is chosen - `.npy` uses the `fnl` flag when present, `.json` uses `itr == max(itr)`.

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
calibration and clustering and re-renders - settings do not live-update the plot.

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
the factor is applied once, **on Plot Data** - so a saved analysis reproduces. Loading one of our own
exports adopts the factor recorded in the file. It multiplies every axial quantity (z spread, cluster
extent in z, 3D neighbour distances, Rg), so it is a calibration to measure on a known standard, not
a display preference. Greyed out when the dataset has no axial variation.

### Metadata Filtering

Filters are applied sequentially in the order listed in the sidebar. All filters operate on the currently active dataset.

```
Raw data
 └── [1] Metadata filters      (EFO, CFR, DCR, TIM, ECO, dwell, fluo)  → removes POINTS
      └── [2] Outlier trimming (localizations far from their own trace) → removes POINTS
           └── [3] TID length filter  (number of points per trace)      → removes TRACES
                └── [4] TID spread filter  (per-trace spread)           → removes TRACES
```

Stages [1] and [2] remove individual localizations; [3] and [4] remove whole traces. Because the two
trace-level stages are both evaluated on the same trimmed point set, their relative order does not
change the result. Trimming, when enabled, runs before both - so a trimmed localization counts
neither toward a trace's length nor toward its spread.

**Metadata filters.** Pick a variable, see its histogram (with the unfiltered data behind it in grey
if *Show raw data* is on), and set a numeric range:

- `EFO` - Effective Frequency at Offset. The averaged emission frequency over all points of the scan
  pattern except the centre; exceedingly high values may indicate multi-emitters.
- `CFR` - Center Frequency Ratio, EFC/EFO. Commonly used to exclude multi-emitter events; **no CFR
  cutoff is applied automatically anywhere in the app** - set it yourself if you want one.
- `TIM` - acquisition timestamp; useful for selecting a time window
- `DCR` - Detector Channel Ratio; photons at detector 1 vs detector 2, used to separate two labels in
  2-colour MINFLUX
- `ECO` - Effective photon Count at Offset: the raw photon count per localization at the periphery of
  the donut
- `Dwell` - dwell time of the localization, in MINFLUX steps
- `FLUO` - fluorophore channel index; a pyMINFLUX-exported column, usable here when present

**TID length filter**: keeps only traces with a number of localizations within the specified range.
Useful for excluding very short (potentially noise) or very long traces (potentially fiducials or
moving emitters).

**Inspect trace spread (σ̃)** opens the spread panel:

- *Min localizations per trace* - traces shorter than this are excluded from the spread histogram,
  the median-spread labels, the allowable ranges and the σ filters. Literature typically uses ≥ 5
  (an SD needs ≥ 2 points to exist at all).
- **Trim outlier localizations within traces** (off by default) - flags and excludes single
  localizations sitting far from the rest of *their own* trace, so a long trace is not discarded
  over two stray points. *Reject beyond k × σ of the trace* sets the cutoff (default `k = 4`). An
  info box reports the outcome measured on the loaded data: how many localizations and traces were
  affected, how many traces lost ≥ 30 % of their points, and how many fell below the trace-length
  minimum as a result. See [Outlier trimming](#outlier-localization-trimming).
- The histogram shows per-TID spreads, optionally overlaying the raw distribution and the median
  spread before/after filtering. Labels report median(σ̃) ± nMAD (nMAD = 1.4826 × MAD), with dotted
  lines at the quartiles. The per-trace σ̃, the median ± nMAD label and the interquartile lines all
  follow pyMINFLUX. The info box reports the same three numbers per axis, alongside the trace count.
- **Shared axes across σ X / Y / Z** (off by default) - puts the three histograms on one set of bin
  edges and one count axis so they can be read against each other and used as a figure set. When
  enabled, *Shared σ axis range (nm)* sets the displayed span (default 0–30); clearing it falls back
  to the full data range. This is a display range only and filters nothing - values outside it are
  omitted from the plot rather than piled into the edge bin, and a note reports how many.
- *Set spread ranges to the 95th percentile* seeds the σ filter ranges from the data.
- **TID spread filter**: per-axis x/y/z bounds. Setting the upper bound low retains only spatially
  tight traces - useful for removing gold beads or freely diffusing molecules. An empty box means no
  bound on that side (blank lower reads as 0, blank upper as no upper limit), so clearing both leaves
  that axis unfiltered. See [Trace Spread](#trace-spread-sd-filter).
- **Remove rogue traces by rule** - an alternative to hand-set bounds. Three rules, each with a
  multiplier *k*: spread > k × dataset median, spread > median + k × MAD (robust), or any single axis
  spread > k × that axis's median. Dropped traces are removed **before** grouping, clustering, fits
  and surfaces.

**Download / Load Filter Settings** saves and restores all filter parameters as JSON.
**Export MINFLUX data (.csv)** writes the currently filtered localizations.

### Visualization Options

- **Point size**
- **Color points by** - Trace/Cluster colour (default), Multiplex channel, Unicolor (one colour for
  everything, with a picker), Depth (z), Time (tim), or Trace-median spatial density. The three
  numeric modes get a colorscale picker and a reverse toggle; density mode adds a radius (nm) that
  counts how many other *trace-median* positions fall within it.
- **Show trace median positions** - one point per TID (or per de-blinked molecule / cluster when
  those are active) instead of individual localizations
- **Dark mode**, **Hide gridlines & axis titles**
- **Add scale bar** (2D) - length, thickness and font size in nm
- **Transform coordinates to positive space** - shifts the dataset so all coordinates are positive
- **Show Alpha Shape (3D)** - wraps each trace (or each de-blinked ID) in a 3D alpha surface, with an
  opacity slider
- **Crop to view** (2D) - zoom into a region and press the button; the X/Y crop is retained when you
  switch to 3D
- **Crop to region (3D)** - an explicit trimmer with three shapes: axis-aligned box, sphere, or
  PCA-aligned ellipsoid
- **Auto-Rotate & Record (3D)** - a spinning animation (speed, elevation, camera distance) and a WebM
  recorder with width/height, FPS and bitrate controls. Start auto-rotate first, then record.

### De-blink correction

Combats the instrument artifact where one fluorophore that blinks several times produces several trace
IDs stacked on top of each other, overcounting the number of molecules. De-blinking merges them back
into one **de-blinked ID** with a shared centroid, so downstream analyses count molecules rather than
blinking events.

Works at **trace level** - a whole trace joins a molecule or it does not, so a trace is never split.

#### The two criteria

**Trace centroid distance** (default) - each trace collapses to its median position, and traces are
grouped when those positions agree. Grouping is *complete*: every trace in a molecule lies within the
threshold of every other trace in it, so a molecule can never be wider than the threshold. Set the
threshold in nm (`Centroid distance threshold`) or as **k × the measured trace spread**; changing k
rewrites the distance.

- *Advantages* - the criterion is a single distance in nanometres, which is easy to state, to compare
  between datasets and to report in a methods section. Because grouping is complete rather than
  chained, a stray pair of nearby traces cannot pull a long string of molecules into one group.
- *Disadvantages* - it needs a distance, so it depends on the measured trace spread being
  representative of the whole file. Heterogeneous data (mixed labels, widely varying photon counts)
  is served less well by one global number. The result is sensitive to the value chosen: too low
  leaves re-blinks counted separately, too high fuses genuinely distinct neighbours.

**Majority vote of localizations** - merge when a fraction *f* of a trace's localizations sit closer
to a neighbouring trace's centroid than to their own. Each localization is compared against its own
trace's centroid computed *without* that localization, so a point is not unfairly close to its own
centre. **Agreement** decides whether both traces must vote for each other (mutual) or either is
enough.

- *Parameters* - *f*, the fraction that must agree (default 0.3). Two traces of the *same* molecule
  score ≈0.5 by construction, so useful values sit well below 0.5; below about 0.25 a single
  localization can decide a short trace. *Distance tolerance* ρ (default 3.0) loosens the "strictly
  closer" test; grouping accuracy rises with ρ and then flattens, so erring high costs little.
- *Advantages* - no distance in nm and no dependence on the measured spread, so it needs neither a
  threshold nor a well-behaved σ. It is self-scaling: a diffuse trace automatically demands
  proportionally more separation than a tight one, which suits heterogeneous data. A stray
  localization casts one vote and is outvoted.
- *Disadvantages* - three controls instead of one, and a criterion that is harder to state concisely
  than "merged within X nm". Because it is relative, there is no single nanometre figure to quote or
  to compare against another dataset.

**De-blinking only merges - it never deletes.** A trace that cannot be merged becomes its own
single-blink molecule, so every trace that survived filtering is still counted. Removing data is the
job of the metadata filters, the optional outlier trimming and the spread filters, where it is
visible and reversible.

#### Choosing the threshold

*Measure threshold from data* sets it from your own file: it measures the per-axis trace spread and
returns **k × that spread**, with *k* defaulting to 2.8. Grouping accuracy is flat between about
k = 2.6 and 3.1, and because a molecule can never be wider than the threshold, overshooting costs less
than undershooting - so the default sits in that plateau rather than at its edge.

The threshold is the trace spread times *k*, and nothing else. Sample drift is measured, but it is
reported in the [Localization Precision Cloud](#panels) QC panel and deliberately kept out of this
number: the only molecules whose motion can be observed are the ones that were already merged - that
is, the ones that moved *less* than the very radius being set - so feeding it back would be circular.

Two diagnostics sit beside it - a density plot whose orange curve should show two humps (re-blinks
left, separate molecules right) with your threshold in the gap, and an exploratory
molecule-count-vs-threshold curve. If the density plot shows no clear gap, no threshold will cleanly
separate re-blinks from neighbours in that dataset, and the molecule count should be treated as
uncertain whichever value is used.

**Anisotropic option** - weights x, y and z by their own measured spread, so one threshold fits all
three axes (the merge volume becomes an ellipsoid in nm, tighter along the better-resolved axis).
Affects the merge decision only; reported distances and sizes stay in true nm.

De-blinked molecules are labelled **De-blinked ID** in hover text, and their medians (not individual
TID medians) become the representative points for downstream analysis including Median Clustering.


> **Naming in the code.** The source calls these two stages `cluster` and `cluster2`: `cluster` is
> de-blinking, `cluster2` is Median Clustering.

### Nearest Neighbor Analysis

Nearest-neighbour distances computed from median positions - de-blinked molecule medians when
de-blink correction is on, trace medians otherwise. Controls: *k*, optional distance labels (font
size and offset), a distance-range filter, and a downloadable table. NN lines are drawn on the plot.

When multiplexing is active a **scope** picker appears (within channel / between channels), and when
Median Clustering is active a **Within Median Cluster only** option restricts the search.

### Median Clustering

The app's one real clustering step. It runs on **representative points** - one point per trace median,
or per de-blinked molecule median when de-blinking is on - and groups them into larger spatial
structures, for example a collection of receptors arranged in a ring, or a row of nuclear pores along
the nuclear envelope.

Three grouping methods:

| Method | Parameters | Behaviour |
|--------|-----------|-----------|
| **DBSCAN** | `eps` (nm, default 200), `MinPts` | Fixed neighbourhood radius |
| **HDBSCAN** | `MinPts` only | Adaptive density; MinPts doubles as minimum cluster size |
| **OPTICS (Xi)** | `Xi` (0–1, default 0.05), `MinPts` | Adaptive density that keeps small structures; lower Xi = more, smaller clusters |

`MinPts` (default 15) applies to all three. Set it near the smallest number of traces you would call a
real structure - often 3–6, not higher.

Optional cleanup steps, applied in this order:

- **Density prefilter** - before clustering, drop centroids with fewer than *n* neighbours within a
  radius, removing scattered background.
- **Core points only** (DBSCAN only) - keep only centroids that reach `MinPts` neighbours within
  `eps`. Plain DBSCAN also absorbs *border* points, which is what attaches lone background centroids
  to the rim of a real structure. HDBSCAN and OPTICS are already defined this way.
- **Prune sparse cluster members** - demote to noise any member whose core distance (the radius it
  needs to reach `MinPts` fellow members) exceeds a multiple of the median core distance of the
  members around it. The baseline is local, so a density gradient inside a structure does not condemn
  its sparser side. ~2 trims clear outliers; lower is stricter.
- **Min members per Median Cluster** - drop whole clusters below a member count. Separate from
  `MinPts` so the density scale and the "how big is a real structure" decision stay independent.

**kNN-distance plot** - sorted distance to the (MinPts−1)-th neighbour of each centroid. The elbow
marks the structure/background split and therefore a good `eps`; the current `eps` is drawn as a red
line.

Settings can be saved and reset with the buttons at the bottom of the group.

#### Visualization

In **2D**: *Wrap Median Cluster in Convex Hull*.

In **3D**: *Wrap Cluster in Alpha Shape* (concave hull via `alphashape3d`), or *Show isosurface* - a
surface render of the cluster, drawn as a contour of a Gaussian kernel density summed over the same
representative points the clustering used (de-blinked molecule medians when de-blinking is on, else
trace medians).

The wrap sits at **2 × the structure's own median centroid spacing** by default, so each structure is
wrapped at its own scale. One slider adjusts it: **padding (−40 to +200 nm)**, negative tightens and
positive loosens. Tightening much below the default stops the surface closing and breaks it into
separate lumps. Smoothing (sigma) and the sampling box follow from the wrap; both can be set manually
behind checkboxes if needed. Remaining controls are grid resolution (12–40) and opacity. The
isosurface can be restricted to one cluster ID and drawn in a single colour or per-cluster colours.

#### Fitting & Metrics (3D only)

Each cluster's representative points are fitted to a geometric model independently. Two fitting
methods:

| Method | Description |
|--------|-------------|
| **Heuristic** | Analytical: geometry derived directly from PCA axes and 5th–95th percentile extents. No iteration, no centre optimization. Fast and always converges. |
| **Weighted orthogonal distance** | Iterative: minimizes the sum of squared Euclidean distances from each representative point to the nearest point on the model surface. Equivalent to an orthogonal regression on the geometry. |

RMSD is reported in nm for both, so the two are directly comparable.

The **weighted orthogonal distance** method has additional options:

- *Weight by per-trace precision* - each representative point is weighted by $n / \sigma_\text{spread}^2$ (proportional to $1/\text{SEM}^2$ of the median position). Traces with a single localization (undefined spread) receive the dataset-median spread as a fallback.
- *Centre handling*:
  - **Fixed** - centre = convex-hull centroid; radius/semi-axes solved analytically (no iteration)
  - **Regularized (MAP)** - centre is free but penalized by a Gaussian prior anchored to the robust centre; strength $\lambda$ controls the trade-off. At $\lambda = 0.1$, a 5 nm centre shift adds ~2.5 nm² to the mean squared residual
  - **Free** - centre fully optimized within the data bounding box; use when the centre is expected to deviate significantly from the data centroid

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

- *Arc span hint* (°) - initial angular extent for arc fitting (0 = auto-detected). The optimiser refines it.
- *Helix turns hint* - seeds the pitch estimate only; the fitted turn count follows from the data extent and the optimised pitch.

**Compute all-model RMSD comparison** (slower) fits every shape model to each cluster and reports all
RMSDs side-by-side, useful for model selection.

**Additional 3D overlays**: PCA axes as coloured arrows; per-cluster summary text (centre, radius,
RMSD) drawn on the plot; *Show distances to fit* (line segments from each representative point to the
nearest surface point, with a maximum line-length filter); and a fit-surface opacity control. The
metrics table is downloadable.

**Membrane protein options** (requires the global reference plane to be active) relate each cluster to
a global PCA-fitted reference surface, e.g. the plasma membrane:

- *Angle to plane* - angle between the cluster's primary PCA axis and the reference plane normal;
  0° = perpendicular to membrane, 90° = lying flat
- *Distance to plane* - signed distance from the cluster centroid to the reference plane
- *Bending radius* - curvature radius of an arc-like cluster

The bending radius is:

```math
R_\text{bend} = \frac{\ell^2 + s^2}{2s}
```

where $\ell$ is the PC1 half-extent (chord half-length) and $s$ is the PC2 extent (sagitta).
**Only valid for arc-like clusters viewed roughly edge-on** - not meaningful for complete rings,
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
threshold / cell size / smoothing radius. Results are a starting point - nudge the origin to refine.

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
Traces below the cut are marked orange - these are the ones that would change channel under a small
threshold shift, and the ones that manufacture false short-range cross-channel pairs. Exportable so
you can exclude or inspect them.

Each group gets its own colour. Additional options:

- **Draw NN lines (2-color mode)** - cross-channel or within-channel nearest neighbours with *k*, a
  distance range, optional labels, and a CSV export.
- **2-color de-blink correction** - the same centroid/vote criteria as the mono-colour de-blinker,
  with its own threshold, *Measure threshold from data* button and gap diagnostic. De-blinked IDs
  carry their channel (`A-1`, `A-2` … / `B-1`, `B-2` …), and **only members of the same channel can
  merge**, so two colours can never be de-blinked into one molecule. In 3D the de-blinked molecules
  can be wrapped in alpha shapes.
- **2-color Median Clustering** - Median Clustering run per channel, with the same DBSCAN / HDBSCAN /
  OPTICS choice and the same cleanup steps (prefilter, core points only, prune, min members), plus
  per-channel alpha shapes and isosurfaces in 3D. The isosurface uses the same density field and the
  same padding dial as the mono-colour one, built from that channel's de-blinked molecule centroids
  (trace medians when 2-color de-blinking is off).

Settings save and reset with the buttons at the bottom of the group.

> Numeric colour modes (z / time / density) override the 2-color colours. Use the default colour mode
> to see the groups.

### Multiplexed MINFLUX

Builds an **N-channel** dataset for exchange / sequential MINFLUX. Two channel sources:

- **Overlay loaded files** - pick several loaded files; each becomes one channel with its own colour.
- **Split active file by time (TIM)** - split one continuous acquisition. Either *contiguous
  breakpoints* (N breakpoints → N+1 channels, up to 7) or *explicit windows* (`start-end, start-end,
  …`), where localizations falling **outside** all windows - the wash / exchange steps - are removed.
  A trace-IDs-over-time histogram shows the split live: bar colour = channel, grey = removed, dashed
  lines = window edges.

Clustering runs strictly **per channel** - molecules imaged in different exchange rounds must never
merge. Channels feed both the plot and the Analysis tab.

### Estimate Ground-Truth Positions

> **Exploratory.** This inverts an assumed labelling model - read it as a hypothesis about the
> labelling, not a measurement.

Reverses the labelling model to propose where the underlying epitopes were. You choose a **labelling
preset** (the same presets as the Simulator: linker length, valency, dye spread, footprint, extra
blinks) or supply custom values, plus a labelling/detection efficiency.

- **Model from** - auto (de-blinked molecules if present, else raw traces), de-blinked molecules, or
  raw trace IDs. Auto is recommended: de-blinking collapses re-blinks and multivalent labels to ~1
  unit per molecule, turning an unstable divide-by-multiplicity count into a divide-by-~1 count.
- **Restrict to points inside Median Clusters** (on by default).
- **Label-complex merge radius** - units within this radius are pooled into one candidate GT, since a
  multivalent label scatters many fluorophores around one molecule. Blank = auto per preset (max of
  footprint and 2.5 × dye spread); 0 = no merge.
- **Stoichiometry-guided merge** *(model-driven)* - where the spatial merge leaves more GTs than the
  labelling stoichiometry predicts, collapse further toward that target. Only ever reduces the count;
  never splits. Sanity-check it against the two-count readout.
- **Displace GT inward** - the antibody binds a surface epitope and the linker throws the dye
  *outward*, so the trace cloud sits ~one linker-length outside the true structure. This pushes each
  GT back toward the local structure centre (perpendicular to a filament axis, toward the centroid
  for a compact or ring shape). Needs Median Clusters defined to know the structure.
- **Support-penalty strength** - how hard a GT backed by few trace IDs is penalized in the
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
- **Surface method** - LOESS heightmap (with a span control), polynomial surface (degree 1–5), or
  thin-plate spline (with a smoothness λ; heavy). All share a grid resolution, opacity, colorscale and
  a *use points from* group selector.
- Specific TIDs can be excluded from the fit by ID.

---

## SMLM Analysis Tab

### Panels

The Analysis tab is organised into sections: Labelling Density and the observation window, Dataset
Groups, Quality Control, Nearest-Neighbour & Spatial Statistics (including Local Density & Hotspots),
Multi-Panel View, Median Cluster Analyses and Cross-Channel Analysis. Each panel produces an interactive plot and, where
a table makes sense, an exportable CSV. Panels marked ⚠ carry an *exploratory* badge in the app -
they are research tools whose interpretation is not yet settled, not validated readouts.

**Dataset Groups** sits above everything: analyse **by individual file** (default) or **by named
group**. In group mode you define up to 6 groups (name, colour, member files) and the distribution
analyses relabel and pool by group so you can compare conditions. Files left unassigned are ignored.
Spatial statistics (Ripley / PCF) average per-file curves within a group rather than pooling
coordinates - pooling would mix observation windows and invent structure.

**The cell (file) is the unit of replication.** A group-mode histogram pools distances or clusters from
every file, and units from one cell are correlated, so a statistical test on those pooled units overstates
the sample size. **Download per-file summary (CSV)** in the Dataset Groups panel writes one row per file
(per channel, per metric) for every distribution analysis you have run: n, mean, SD, median and quartiles,
with the file's group. Run between-group statistics on those rows.

### De-blinking applies to every dataset

De-blinking is a **global** switch (Display tab), and every batch analysis re-derives it per file, so
enabling it once applies it to all loaded datasets - there is nothing to set per file.

The panels with a **Unit** selector (Batch NN, Ripley's L, g(r)) default to **Auto**, which resolves to
de-blinked molecules when de-blinking is on and trace IDs when it is off, so every panel follows the
same de-blinking state. Choosing a molecule unit explicitly while de-blinking is off falls back to trace
IDs with a warning rather
than returning nothing.

Note the **Crowding** panel needs unclustered background traces to reference its timing shift against.
If every molecule ends up inside a median cluster it will report nothing - that is the background rule,
not a de-blinking problem.

### Observation window and edge correction

Two checkboxes at the top of the tab, inside the collapsed labelling-density section. They apply to
Ripley's L, g(r), L12, g12 and the co-localization null, and they fix **different** things - usually only
the first is wanted.

**Edge correction.** A point near the edge of the window has fewer neighbours available to it than one in
the middle, simply because part of its surroundings falls outside the window. Left uncorrected, that
deficit grows with r, so every curve drifts downward at larger radii and a genuinely random pattern stops
reading as random. The correction accounts for how much of each neighbourhood actually fell inside the
window (the translation correction; Illian et al. 2008; Baddeley et al. 2015). It is on by default.

**Footprint window.** Changes what counts as the observed area. If your cell fills the scan field, the
bounding box already *is* the cell and you should leave this off. Turn it on when the cell is concave or
only partly in frame - then the box contains a large empty region, the density is under-stated, and the
pattern reads as clustered purely because of how it was framed.

Two things to know:

- The settings are **not independent**. Correcting the edge of a bounding box that contains a lot of
  empty space is worse than not correcting it at all, because it amplifies a deficit that was partly
  cancelling the area error. On a concave cell, use both together rather than the first alone. When the
  footprint window is off, the window readout reports how much of the bounding box the sample covers and
  asks you to turn the footprint window on when that is below 70%.
- The window is the **outline of the cell**, not the region where the sample happens to be densest. A
  window drawn tightly around your structures would erase the clustering it is meant to measure, so the
  footprint is estimated independently of the Labelling Density panel, and if it comes out in
  disconnected pieces the bounding box is used instead and you are told.

The co-localization null also randomises within this window rather than within the bounding box, which is
what makes it trustworthy for a cell that does not fill its frame. It is now conservative rather than
exactly calibrated: safe to quote as "not significant", but for a positive claim prefer g12(r) read
against its envelope. The two-colour statistics get the window and the edge correction; the support map
in the Labelling Density panel shows which region is in use.

### What a Ripley L(r)-r curve should look like

- **Random:** flat at 0 across the whole range.
- **Clustered:** rises from 0, peaks somewhere near the characteristic cluster scale, then decays slowly
  back toward 0 as larger radii average over the whole field.
- **Regular / dispersed:** dips below 0 at short range.

A curve that heads **straight down into the negative in a near-straight line** is not a biological result.
K scales with the window area, so if the area being divided by is much too small the whole curve is
pulled down proportionally, which looks like a line. The usual cause is a window drawn around the
structures instead of around the cell - see the decision rule above.

Note also that L(r)-r **does not return to zero** once clustering stops, and that is expected rather than
a sign of structure at larger scales. K counts every pair closer than r, so an excess at short range stays
counted at every larger radius; the curve decays only slowly from there. That is also why the position of
the L(r)-r peak is *not* the cluster size - read the size from g(r) instead, at the radius where it comes
back down to 1.

### How filtering works on this tab

Every analysis on the Analysis tab applies the **Display-tab filter ranges** - EFO, CFR, DCR, ECO,
dwell, fluo, trace-ID length, trace spread and outlier trimming - identically to *all* loaded files,
not only the active one. There is no separate per-file filtering to set up.

The one thing to watch: the sliders' *limits* are auto-detected from the **active file**, so if another
dataset runs higher on a metric, check that a range which looks like "no filter" still means that
before comparing groups. The time (`tim`) range is per-file by nature, since each acquisition has its
own timeline. This is stated as a note at the top of the tab.

### Labelling density normalisation (opt-in) ⚠

Collapsed by default, and off by default - most experiments do not need it.

It exists for the case where labelling density varies between acquisitions for reasons that are not
biology: transient transfection giving each cell an undefined plasmid load, differing expression level,
run-to-run labelling efficiency. That variation is not cosmetic, and it pushes in the direction that
*mimics* aggregation. Nearest-neighbour distances shrink as density rises, and DBSCAN at a fixed eps and
minPts is itself a density threshold, so a high expresser reads as more clustered even when nothing about
its organisation has changed.

**Measure** reports, per file, how many units there are and the density over an estimated sample area:

- **Unit** is a de-blinked molecule when de-blinking is on, and a trace ID otherwise. This matters: a
  trace ID is a single blink, so trace-ID density mixes expression level together with photophysics and
  laser power, and is not comparable between runs.
- **Area** is the region the sample actually occupies, not the bounding box. A MINFLUX field routinely
  contains the edge of a cell and interior dark regions, and a bounding box counts all of that as sample,
  so a cell imaged half in frame reads at half its true density purely from framing. Four estimators are
  offered; the default smooths the localizations into a density map and keeps the smallest region holding
  a given share of the units (Hyndman 1996). Dark regions fall below the cut on their own and there is no
  fill radius to tune. The alternatives are an occupancy grid with a morphological fill, the convex hull,
  and the bounding box.
- **Density (NN)** is a second estimate from nearest-neighbour distances that needs no area at all, so it
  is a useful check on the area-based number. It is deliberately not used as the thinning target, because
  it measures density locally and clustering inflates it.

The plot shows one **horizontal violin per dataset** over the per-unit local density, with the overall
density marked as a diamond. The single number tells you how much label there is; the shape tells you how
it is arranged. A cell with real aggregates is heterogeneous - a sparse background plus a dense tail -
where a uniformly labelled cell at the same average is narrow. Keep the smoothing settings identical
across the samples you compare: together they define what "the sample area" means.

**Show support map** draws one dataset with the counted area outlined in green and the disregarded part
shaded red, over the localizations. This is what makes the area estimate checkable - the number says
"30% of the frame", but only the picture says whether that 30% is the *cell* or a tracing of your own
*structures*. Two honest limits: the frame is the bounding box of the data, because a localization table
carries no record of the real scan field, so red only ever marks interior dark regions and concave
margins; and in 3-D the map is a plan view, so the area looks larger in projection than it is.

**Thin all files to a common density** is the second, separate opt-in. Keeping each unit with a fixed
probability leaves the pair correlation g(r) and Ripley's K unchanged, so real clustering survives
thinning while a difference caused only by labelling density does not. Thinning every file to one density
therefore makes even the fixed-threshold DBSCAN read-outs comparable. It drops whole units - a trace or a
molecule with all its localizations - after the metadata filters and before clustering. Quality-Control
panels are exempt: they must see the data as acquired.

**Scope: what thinning is and is not for.** It matters only for **count-based** read-outs - DBSCAN cluster
counts, % in clusters, cluster sizes, nearest-neighbour distributions - which change with labelling
density by construction, and are not comparable between files of differing density unless you match first.
**Ripley's L and g(r) do not need it at all**, because they already divide the density out; thinning those
discards data without changing the answer. The useful summary: **for K, L and g the observation window is
everything and the labelling density is nothing** - they divide density out but multiply by the area,
which is why the area estimate and the edge correction matter so much and this section usually does not.

With thinning off, the density measurement is still worth having: it tells you how comparable your
acquisitions are, whether you have enough molecules for a stable g(r) estimate, and whether density
happens to confound with condition.

Two limits the panel will not exceed. It **never up-samples**, so a file already below the target is left
untouched; the default target is the lowest density in the set for exactly that reason. And thinning
**cannot separate density from biology when the two are genuinely linked** - for a self-aggregating
protein they are, since aggregation depends on concentration. A difference that survives at matched
density is evidence; a difference in how clustering *scales* with density is stronger evidence, and needs
a spread of expression levels per condition rather than one matched pair.

### Two-colour method details

The mono-colour statistics (NN, Ripley's L, RoG) are described under
[Key Calculations](#key-calculations). The four cross-channel panels:

#### Cross-Channel Nearest Neighbor & Co-localization Curve

Directional nearest-neighbour distance from each point in one group to its nearest point in the other,
computed in **2D (x, y)** with `RANN::nn2`:

```math
d_i^{G_1 \to G_2} = \min_{j \in G_2} \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}, \quad i \in G_1
```

Both directions are computed independently. The **co-localization curve** is the empirical CDF:

```math
F(r) = \frac{1}{|G_1|}\; \#\!\left\{i \in G_1 : d_i^{G_1 \to G_2} \leq r\right\}
```

- the fraction of group-1 points with a group-2 neighbour within *r*. A steep rise at small *r*
indicates strong co-localization; a gradual rise suggests coincidental proximity.

#### Bivariate Ripley's L₁₂(r)

Tests whether two populations co-localize beyond chance given their individual densities.

```math
\hat{K}_{12}(r) = \frac{A}{n_1 n_2} \sum_{i \in G_1}\sum_{j \in G_2} \mathbf{1}\!\left(\|\mathbf{p}_i - \mathbf{p}_j\| \leq r\right)
```

with the same L-transform as the univariate case. The plotted L̂₁₂(r) − r: > 0 = co-localize
more than two independent CSR processes, ≈ 0 = spatial independence, < 0 = exclusion.

When Median Clustering is active, L̂₁₂ is computed **per cluster** and averaged, giving a
within-structure measure rather than a global one. The 2000-point subsampling limit applies per group.

#### Pair Cross-Correlation g₁₂(r)

The annulus-normalized derivative of $K_{12}$ - more sensitive to the *scale* of co-localization than
the cumulative $L_{12}$:

```math
g_{12}(r) = \frac{A}{n_1 n_2} \cdot \frac{\#\text{cross-pairs in annulus}\bigl(r-\tfrac{\delta}{2},\, r+\tfrac{\delta}{2}\bigr)}{\pi\!\left[(r+\tfrac{\delta}{2})^2 - (r-\tfrac{\delta}{2})^2\right]}
```

with the 3D form using a spherical shell and $\tfrac{4}{3}\pi[(r+\tfrac{\delta}{2})^3 - (r-\tfrac{\delta}{2})^3]$.
Under CSR with independent groups, $g_{12}(r) = 1$; > 1 indicates co-localization at scale *r*, < 1
exclusion. $L_{12}(r)$ accumulates all pairs up to *r* and can mask scale-specific effects - the PCF
pinpoints them.

#### 2D KDE Density Map & Co-localization

KDE maps are computed for each group on a shared grid (`MASS::kde2d`) and rendered as heatmaps.
Bandwidth follows Silverman's rule of thumb:

```math
h = 0.9 \cdot \min\!\left(\hat{\sigma},\; \frac{\text{IQR}}{1.34}\right) \cdot n^{-1/5}
```

computed on the combined x-coordinates of both groups and applied identically in x and y, so the two
channels share a grid. Co-localization is the Pearson *r* between the two density surfaces.

> **Caveat**: the Pearson *r* between KDE surfaces is strongly bandwidth-dependent. A large bandwidth
> smooths both channels and inflates the correlation even for spatially distinct populations. Use it
> for relative comparisons across conditions, not as an absolute co-localization measure.

### Export

- **Export MINFLUX data (.csv)** - the currently filtered localizations
- **Download Filter Settings** - all filter parameters as JSON, reloadable
- **NN distance tables** - CSV export of nearest-neighbour results (mono-colour, 2-colour, and within-cluster)
- **Median Cluster metrics table** - per-cluster morphology and fit results
- **Per-trace split confidence** - 2-colour posterior probabilities

---

## Tracking Display Tab

Sidebar groups appear in this order.

### Data Upload

Accepts `.npy`, `.json`, `.csv` (pyMINFLUX or this app's export), and the fast cache formats
`.qs` / `.fst` / `.rds`. Multiple files merge into one dataset with track IDs prefixed `F1_`, `F2_`
so IDs from different files never collide.

Large files are **chunked** on upload - split into blocks of whole TIDs written to a manifest
directory - so rendering stays responsive on modest hardware. Chunk size is set by *TIDs per chunk*.
The selected chunk scopes the **playback animation**; tick *Limit Full tracks to the selected chunk*
to scope the static view as well (off by default, where Full tracks draws every track). The filter
pipeline and every analysis always read **all** TIDs, so nothing is ever restricted to the chunk you
are viewing.

**No filter is applied automatically on upload** - CFR in particular is left alone, since it is not a
strong artefact discriminator for tracking data. Filter it yourself via the CFR range below if you
want to.

### Axial (z) calibration

Same control as the SMLM tab: an absolute refractive-index factor (default 0.7), applied once after
the file is read raw. Set both tabs to the same value when comparing axial distances between them.
Greyed out when the dataset has no axial variation.

### Acquisition Parameters

Optional hardware model for cycle time, a manual cycle-time override, and a manual σ_loc override.
Cycle time drives gap-splitting and the MSD lag axis, so getting it right matters more here than
almost anywhere else - the status line shows the value in use and where it came from.

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
| **Split at time gaps** | `Split tracks at time gaps`, `Adaptive threshold`, `Gap threshold (cycles)`, `Min localizations per segment` | A MINFLUX track ID can span a long dark interval and rejoin - which is not one continuous trajectory. Splits at intervals longer than the threshold and discards segments too short to analyse. Adaptive mode derives the threshold from the track's own interval distribution. |
| **Trim degraded ends** | `Trim dark-gap-separated ends`, `Dark-gap threshold (cycles)`, `Trim dim (low-EFO) ends`, `EFO threshold (× track median)` | Bleaching tails: the end of a track is often sparse and dim before it dies. Trims trailing localizations separated by a long dark gap, or whose EFO has fallen below a fraction of the track median. |
| **Continuity gate** | `Track-continuity gate`, `Dark-cycle tolerance (half-life, cycles)`, `Min continuity score` | Drops whole tracks that are mostly dark cycles. The score weights each localization by how many cycles were skipped before it, with an exponential half-life, so a track that only nominally exists across a long window is rejected. |

Then the metadata filters: pick a variable (EFO, CFR, DCR, ECO, track length, duration…), see its
histogram with the unfiltered data behind it in grey, and set a range. **Crop to view** restricts to
a drawn box, keeping only tracks with at least a set percentage of their localizations inside.

Two things to know about how filters apply:

- **Filters do not live-update the plot.** They apply on **Plot Data**. The histogram and the
  read-outs update live so you can see what you are selecting; the rendered tracks stay put until you
  ask. This is deliberate - re-rendering a WebGL track view on every slider move is slow and makes the
  filters feel laggy.
- **Plot Data is the single apply action.** It runs the cleaning stages, applies the filters,
  re-chunks the surviving tracks at the current *TIDs per chunk*, and re-renders. There is no separate
  rechunk step. *Reset ranges* restores the original uploaded manifest.

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

The foundation panel. Estimates **D** by the CVE method (Vestergaard 2014), using only single-cycle
steps, so no MSD fit is needed. It reads two moments of the displacement series, ⟨|Δr|²⟩ and the
covariance of consecutive steps, and the covariance is negative for localization noise, which is what
separates motion from noise. The panel also sets the app-wide **σ_loc**, which comes from the photon
budget rather than from the covariance, and holds the instrument constants everything else reads.

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
χ² degrees of freedom for the confinement band. **Defaults to volume** - set it to surface for
membrane work or every D and confinement call below inherits the bias.

### MSD (diffusion)

Mean-squared displacement, binned by **integer cycle-count lag** rather than wall-clock time
(MINFLUX cycle-skipping makes Δt non-uniform), following the cycle-based workflow of Vogler et al.
2025. Fits free, confined and anomalous models. The linear fit uses Michalet & Berglund optimal-point
selection rather than a fixed lag range; the confined and anomalous fits deliberately use the whole
window, because that optimality theory is derived for free diffusion only. A **cycle-binning
diagnostic** shows how steps distribute across cycle counts.

### Confinement (MSD index)

Confinement ratio from the MSD curve's departure from linearity, with a simulated null for
comparison. Confined percentages are **not comparable across different track lengths** - a longer
track has more chance to reveal a plateau.

### Track Segmentation - Change-Point + GMM

Splits tracks at statistically detected change points, then clusters the resulting segments by their
diffusive behaviour into a Gaussian mixture. Use when a single D per track is clearly the wrong model.

### Jump Distance Distribution (JDD)

Fits the distribution of single-step displacements as an N-component mixture - **Rayleigh in 2D,
Maxwell in 3D** - with σ_loc bias correction and an integrated-steps-only mode. Complementary to MSD:
it uses the full displacement distribution rather than its second moment, so it separates populations
that share a mean. Steps must be pooled at a fixed real Δt; mixing cycle counts pollutes the fit.

### Displacement vs Diffusion

Compares each track's net displacement against the free-diffusion expectation √(4DT) for its own
duration, with a Rayleigh null. Tracks far above the band are candidates for directed transport;
far below, for confinement. Uses the motion-DOF-corrected D and the χ² degrees of freedom from
**Sample geometry**.

### Directionality / Autocorrelation

Velocity autocorrelation and turning-angle statistics - persistent directional motion shows positive
correlation at short lags where free diffusion shows none.

### Diffusion State Segmentation (HMM) and Per-Segment MSD

A hidden Markov model over the displacement series assigns each step to a diffusive state, giving
per-state D and dwell times, and the transition structure between them. **Per-Segment MSD** then
re-runs the MSD independently within each HMM state, so a two-state track yields two curves instead
of one average that describes neither.

### Speed & Transport

Drift velocity and a noise-corrected local α(t). Path speed alone cannot distinguish directed
transport from fast diffusion - a random walk sampled finely has a high path speed - so the drift
term and the local anomalous exponent are the informative quantities here. An **Advanced** block
exposes candidate directed segments with per-event corroboration.

### Spatially Resolved D-map

Bins the field of view and estimates D per bin, showing where in the sample diffusion differs.
Subtracts the σ_loc noise floor so slow regions are not reported as diffusing.

### Inspect a single track

Per-track inspection: raw coordinates, intervals, per-step displacements and the fits for one
selected track. The first place to look when an aggregate result seems wrong.

### Where the explanations live

Panels on this tab print **numbers and short verdicts only**. Every method's rationale - what it
assumes, how it is adapted to MINFLUX, and where it fails - is in the collapsed
**ℹ Method & MINFLUX adaptation** box at the head of each panel, with the maths in *Key Calculations*
below. Warnings still appear inline, but as one line each rather than a paragraph, so a result that
needs attention is visible instead of buried.

---

## MINFLUX Simulator Tab

Two phases, each its own sub-tab.

### Phase 1 - SMLM (synthetic localizations)

Generates labelled synthetic structures in 2D or 3D. Two **simulation modes**:

- **Normal** - visualize one structure with one labelling strategy.
- **Acquisition** - lay out a field of view, place structures in it, and model the whole acquisition
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

**Labeling Strategy** - each labelled molecule draws `K ~ max(1, Poisson(n_fluor))` fluorophores,
each detected independently at the labelling efficiency, and each detected fluorophore becomes one or
more trace IDs (blinking splits it further). Presets cover direct tags (SNAP / Halo / Click, DOL = 1),
nanobodies, primary monoclonal and polyclonal IgG (with and without NHS multi-labelling), primary +
secondary (nanobody or IgG) sandwiches, streptavidin, and a Custom option. These are the same presets
the SMLM tab's *Estimate Ground-Truth Positions* inverts.

Two consequences worth knowing before choosing one: a polyclonal + secondary IgG sandwich produces
**~32 trace IDs per protein**, and the shared linker sets a resolution **floor** that no amount of
acquisition time removes (~50 nm for that sandwich, ~6 nm for a nanobody). The full generative model -
dye counts, the three displacement levels, TID multiplication, and a table of every preset's
TIDs/molecule, detection probability and resolution floor - is derived in
[docs/simulator.md §1.2](docs/simulator.md).

**Localization precision** - either a manual σ (separate XY and Z) or a **CRLB from photon budget**.

**Acquisition** *(acquisition mode)* - acquisition time in hours/minutes, structure and background
detection rates, and a **Calibrate from a real run** loader: point it at a real `.npy` / `.json` /
`.csv` and it fills the FOV, duration and detection-rate profile from that data, then overlays the
measured rate on the Acquisition plot. A *structure fraction of detections* slider (auto-estimated
from spatial density) splits the rate into structure and background components.

**Main-panel tabs**: *MINFLUX* (the synthetic localizations), *Confocal view* (a diffraction-limited
render of the same ground truth - objective NA, pinhole in Airy units, pixel size, photons per
fluorophore, background, Poisson shot noise), *Acquisition* (the timeline plot and summary), and
*Compare strategies* (two labelling strategies A and B side by side on identical ground truth).

Simulated localizations and the ground-truth positions export as CSV; re-upload them in the SMLM tab
to analyse them through the normal pipeline.

> Full simulator reference - the generative pipeline stage by stage, the acquisition maths, and the
> known fidelity gaps - is in [docs/simulator.md](docs/simulator.md).

### Phase 2 - Diffusion / Tracking (ground-truth simulator)

> **Exploratory validation tool.** Its purpose is to check that the Tracking tab's estimators recover
> a *known* D.

Free Brownian motion on a cell shape, sampled the MINFLUX way - cycle clock, dark cycles,
heteroscedastic σ_loc and photobleaching - producing MINFLUX-format tracks with a ground-truth D.
Download the tracks and load them in the Tracking tab to check CVE / MSD / confinement.

- **Cell geometry**: sphere (coccus), rod / spherocylinder (E. coli, with a choice of long axis),
  ellipsoid (ovoid / coccobacillus), membrane plane (a tilted flat sheet in 3D, no curvature → no
  projection bias), or free (unbounded 3D - the simplest null).
- **Diffusion**: surface (2 DOF, membrane/wall protein) or volume (3 DOF, reflecting off the wall).
- **Motion type**: free, confined (OU tether → MSD plateau at ~2·L_c²), directed (drift velocity along
  the rod long axis), or two-state (mobile ↔ immobile, ground truth for HMM state recovery).
- **Survivorship** - optionally lose the molecule when the *true* displacement across a dark gap
  exceeds a re-acquisition radius. This is the dark-cycle bias that manufactures apparent confinement.
- **Sampling model**: cycle time, dark-cycle probability, mean photons per localization, and the
  effective (final-iteration) TCP diameter L, from which σ_loc = L / (2·√N_photons).
- **Calibrate from real MINFLUX data** - load a real file to auto-fill cycle time, photons, dark-cycle
  probability and median track length. Tip: apply your filters in the Tracking tab first, export the
  filtered CSV, and calibrate on that, so the model matches exactly the tracks you analyse.

Both the tracks and the ground truth export as CSV.

---

## Key Calculations

### Unit Conversions & Import

Abberior files store coordinates in **meters**. On load, the app converts to nanometers:

```math
x_\text{nm} = x_\text{m} \times 10^9
```

**Z-axis scaling** is a user-set calibration, not a fixed constant. Files are read with **raw z**;
the factor you set (SMLM sidebar *Axial (z) calibration*, Tracking sidebar *Axial (z) scaling factor*,
both defaulting to 0.7) is applied once afterwards, on Plot Data. CSVs are first returned to raw z -
the app detects what factor they already carry, exactly when the file was exported from here (a
`z_scale_applied` column) and by assumption otherwise. Exports write that column so a round-trip is
exact. Because the factor multiplies every axial quantity - z spread, cluster extent, 3D neighbour
distances, Rg - applying it twice would compress z by 30%, which the detection prevents.

Scaling is applied after unit conversion, to correct for the refractive-index mismatch between the
sample medium and the objective immersion medium:

```math
z_\text{corrected} = z_\text{nm} \times s_z \quad (s_z = 0.7 \text{ by default})
```

**Multi-iteration MINFLUX**: each TID produces localizations at several iteration steps (`itr`). By
default only the final iteration is retained - chosen by the `fnl` flag for `.npy`, by
`itr == max(itr)` for `.json`. The CFR value used is the last non-zero CFR recorded across all
iterations for that TID (*last_nonzero* carry-forward), because CFR is often only populated in early
iterations.

**Valid-only filter**: records with `vld = FALSE` are discarded on load.

---

### Dwell Time

The dwell time (number of MINFLUX measurement steps per localization) is derived from the photon-count columns:

```math
\text{dwell} = \mathrm{round}\!\left(\frac{\text{eco}}{\text{efo} / 1000} \cdot \frac{1}{T_\text{step}}\right)
```

where `eco` is the total estimated photon count, `efo` is the emission rate in kHz, and $T_\text{step}$ is the user-specified dwell time per step in ms (default 1 ms). The ratio eco / (efo / 1000) gives the effective dwell duration in ms; dividing by $T_\text{step}$ converts to step counts.

---

### Trace Spread (SD Filter)

For each TID, the per-axis standard deviation of all localizations is computed:

```math
\sigma_x(\text{tid}) = \mathrm{SD}(\{x_i : \text{tid}_i = \text{tid}\})
```

and analogously for $\sigma_y$ and $\sigma_z$.

This is the plain standard deviation, normalized by $n-1$, and it is set up to match pyMINFLUX conventional computation of per trace ID (`processor/_processor.py`, `calculate_statistics_on`:
`sx = df_grouped["x"].std()`), so precision numbers should be comparable between the two tools. MINFLUXR allows the user to filter out traces with extreme trace spreads.

The SD filter ranges set bounds $[\sigma_\text{min}, \sigma_\text{max}]$; a TID is kept only if its
$\sigma$ falls within those bounds on every active axis. Setting the upper bound to a small value
retains only spatially tight traces (useful for removing freely diffusing molecules), while a broad
range passes all traces. Traces below the *Min localizations per trace* setting are excluded from the
statistic entirely. When outlier trimming is enabled, $\sigma$ is computed on the retained
localizations.

The **spread panel** histogram shows the distribution of per-TID spreads and overlays the current
filter range. The median spread is reported as median(σ̃) ± nMAD, with nMAD = 1.4826 × MAD - the
1.4826 makes nMAD a consistent estimator of σ for normal data - and the quartiles are drawn as dotted
lines beside it. This matches pyMINFLUX convention summarizing the same distribution, which labels
its median line with median ± MAD and marks the interquartile range alongside it.

nMAD and the interquartile range both describe how σ̃ varies **between traces**. Axes are reported
separately rather than pooled into one isotropic figure, because axial and lateral precision differ.
σ̃ is a per-localization spread and an upper bound on the instrument's precision.

---

### Outlier Localization Trimming

Optional, off by default. Flags individual localizations that do not fit their own trace, so a trace
is not discarded whole because of a few stray points. It removes points only - deciding the fate of a
whole trace stays with the spread filter and the rogue-trace rules.

For trace $t$ with $n_t$ localizations, the per-axis centre is the median and the distance is
whitened per axis, so the anisotropy of MINFLUX precision does not leak into one threshold:

```math
r_i^2 \;=\; \sum_{a \in \{x,y,z\}} \left( \frac{a_i - \mathrm{median}_t(a)}{\hat{\sigma}_{a,t}} \right)^{2},
\qquad \text{flag } i \text{ if } r_i > k
```

**The scale.** $\hat{\sigma}_{a,t}$ is a robust (MAD-based) spread, shrunk from the trace's own
estimate toward the dataset-wide one:

```math
\hat{\sigma}_{a,t} \;=\; w_t\,\sigma^{\text{trace}}_{a,t} \;+\; (1-w_t)\,\sigma^{\text{pooled}}_{a},
\qquad w_t = \frac{n_t}{n_t + 6}
```

Short traces borrow the dataset's stability; long traces trust their own estimate. Both levels use a
robust estimator because an SD-based scale would be inflated by the very outliers the test is meant
to find. This is internal machinery - the reported σ̃ remains the plain SD above.

Asking whether a point is unlike *its own trace*, rather than unlike the dataset, is what makes a
uniformly broad trace (drift, a fiducial, two merged molecules) lose nothing here: all of its points
agree with each other, and judging it is the trace-level filters' job.

**Guards.** Traces shorter than 6 localizations are never trimmed - at that length a single point is
a large share of the evidence and the median defining "far" is itself poorly determined. Trimming
also never reduces a trace below the *Min localizations per trace* setting.

**Choosing $k$.** Gaussian theory is not a useful guide here: real MINFLUX within-trace residuals are
heavier-tailed than Gaussian, so a $k$ that rejects almost nothing in theory rejects a few percent in
practice. The info box therefore reports the rejection rate **measured on the loaded data** rather
than a theoretical figure. Note that trimming lowers the reported σ̃ partly by genuinely removing
outliers and partly by selection, since the retained points are those closest to the median; a
conservative $k$ keeps the selection component small.

**Interaction with the length filter.** Trimming shortens traces, so a trace can drop below the TID
length minimum and be removed by that filter even though it passed before. This is the correct
reading - a trace with four valid localizations is a four-point trace - and it is reported rather
than prevented, because tying the trim guard to the length minimum would retain known-bad
localizations purely to satisfy a count, corrupting the median position and spread that the length
filter exists to protect.

---

### DCR-Based 2-Color Separation

**Per-TID median DCR** is used as a robust summary to classify each trace into one of two spectral groups:

```math
\widetilde{\text{dcr}}_\text{tid} = \mathrm{median}(\{\text{dcr}_i : \text{tid}_i = \text{tid}\})
```

**Manual mode**: the user defines two ranges $[d_{1,\text{lo}}, d_{1,\text{hi}})$ and $[d_{2,\text{lo}}, d_{2,\text{hi}}]$. A TID is assigned to group 1 or 2 if its median DCR falls in the respective range; TIDs outside both ranges are unassigned.

**Auto-detect mode (GMM)**: fits a two-component Gaussian mixture model to the distribution of those
per-TID median DCR values using the `mclust` package:

```math
p(\text{dcr}) = \pi_1 \,\mathcal{N}(\mu_1, \sigma_1^2) + \pi_2 \,\mathcal{N}(\mu_2, \sigma_2^2)
```

Components are sorted so $\mu_1 < \mu_2$ (group 1 = lower DCR). The split threshold $d^*$ is the crossing point of the two weighted component densities:

```math
\pi_1 \,\mathcal{N}(d^*;\mu_1, \sigma_1^2) = \pi_2 \,\mathcal{N}(d^*;\mu_2, \sigma_2^2)
```

solved numerically in $(\mu_1, \mu_2)$ with `uniroot`. If the two component means differ by less than
0.05 DCR units, the distribution is considered unimodal and auto-detection is aborted with a warning.
The same fit supplies the **per-trace posterior** shown in the split-confidence panel.

---

### Spatial Density

For each point (localization or TID median), the local density is the count of other points within a sphere (or disk in 2D) of radius $r$:

```math
\rho_i = \#\{j \neq i : \|p_i - p_j\| \leq r\}
```

Computed via `dbscan::frNN` (fixed-radius nearest-neighbor search), which uses a kd-tree for $O(n \log n)$ performance. The radius $r$ is user-controlled (default 500 nm in the colour-by-density mode, 150 nm elsewhere).

---

### Nearest-Neighbour Distances

The $k$-th NN distance for point $i$ is:

```math
d_i^{(k)} = k\text{-th smallest} \;\|\mathbf{p}_i - \mathbf{p}_j\|, \quad j \neq i
```

computed with a kd-tree (`RANN::nn2` / `dbscan::kNN`). Histograms are summarized by the median and by
either the IQR or the MAD as a spread indicator.

> The first-NN distance distribution reveals the typical inter-molecule or inter-structure spacing. A
> narrow peak indicates a regular, lattice-like arrangement; a broad tail suggests heterogeneity.
> **A sharp peak below ~40 nm is usually the label footprint* - the distance-range
> filter trims the histogram but does not correct that bias. NN carries no CSR reference; use
> Ripley's L(r) with its envelope if you need significance.

---

### Ripley's L(r)

A spatial clustering statistic that tests whether points are more aggregated than expected under
complete spatial randomness (CSR). MINFLUXR uses a **naive (boundary-uncorrected) estimator**:

```math
\hat{K}(r) = \frac{A \cdot 2}{n(n-1)} \sum_{i < j} \mathbf{1}(\|\mathbf{p}_i - \mathbf{p}_j\| \leq r)
```

where $A$ is the bounding area (or volume $V$ in 3D) and $n$ is the number of points. The
L-transformation:

```math
\hat{L}(r) = \sqrt{\frac{\hat{K}(r)}{\pi}} \quad \text{(2D)}, \qquad \hat{L}(r) = \left(\frac{3\,\hat{K}(r)}{4\pi}\right)^{1/3} \quad \text{(3D)}
```

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

```math
R_g = \sqrt{\frac{1}{m}\sum_{k=1}^{m}\left[(\tilde{x}_k - \bar{x})^2 + (\tilde{y}_k - \bar{y})^2 + (\tilde{z}_k - \bar{z})^2\right]}
```

Only clusters with $m \geq 2$ are included (a single point has $R_g = 0$ by definition). The
exportable table reports cluster ID, point count, centroid coordinates and $R_g$. RoG is
small-*N*-biased - it is reported with its anisotropy and should not be used as a classifier on its
own.

---

### Time Signature Analysis (TSA) & Cycle Time

Everything Δt-aware in the Tracking tab rests on the base **cycle time** $t_c$ - it drives
gap-splitting, the quality gate, cycle-index MSD binning, the D-map and the CVE integrated-step
selection. Get it wrong and every downstream number moves. The TSA panel is a faithful R port of Bela
Vogler's MINFLUX Time-Signature Analysis (Vogler et al. 2025, *Commun. Biol.*).

#### The interval law

MINFLUX has no fixed clock. The molecule is localized once per cycle *if it emits or is detected*; during a dark
cycle nothing is recorded, so the next localization simply lands one or more cycles later. Intervals
are therefore **hard quantized by a TCP cycle time**:

```math
\Delta t = k \cdot t_c + \text{SW}, \qquad k \in \lbrace 1, 2, 3, \ldots \rbrace
```

where SW is a constant per-step overhead. Inverting it classifies each step:

```math
n_\text{cyc} = \mathrm{round}\!\left(\frac{\Delta t - \text{SW}}{t_c}\right), \qquad n_\text{dark} = n_\text{cyc} - 1
```

Steps with $n_\text{cyc} = 1$ are **integrated** (single-cycle) steps - the constant-Δt subset that CVE,
JDD and the directionality tools restrict themselves to.

#### Estimating $t_c$

In priority order: (1) the hardware model $t_c = N_\text{iter} \cdot t_\text{dwell}$ if enabled;
(2) a manual override; (3) a **data-driven cluster cycle** - cluster the observed Δt values, take
$t_c = \mathrm{median}(\mathrm{diff}(\text{cluster centres}))$ and
$\text{SW} = \text{centre}_1 - t_c$; (4) a regression from the Δt-vs-index slope.

Method (3) deliberately uses the **spacing between clusters**, not the mode. This is the trap:

> ⚠ **The modal Δt is not the cycle time** - it is $1 \cdot t_c + \text{SW}$. Mistaking one for the
> other biases every D by about **−5.5%**.

#### What the interval structure looks like

The comb is essentially exact: each cycle-count peak is a near-delta spike rather than a broadened
distribution, with a within-cycle spread of a few hundredths of a percent. Two properties of the panel follow from that:

- **Binning is one bar per cycle count**, with edges on the classification boundaries
  $\text{SW} + (k \pm \tfrac{1}{2}) t_c$. Fixed-width bins alias against the comb and put bar centres
  off the true peaks; one-bar-per-cycle makes the histogram literally a picture of the
  $n_\text{cyc}$ the analysis computes. Fixed-width binning is used only as a fallback when no cycle
  time is available. A log count axis is offered, since counts span ~100× from the 1-cycle peak into
  the dark tail.
- **The Δt distribution is shown as a discrete distribution over cycle counts**, via the histogram plus
  optional cycle-multiple markers. There is no kernel-density overlay: with near-delta peaks its width
  would be pure bandwidth rather than timing jitter.

Because intervals are exact integer cycles plus a constant offset, cumulative-cycle-index binning is
not an approximation to noisy timing - it recovers the true elapsed time. This is what the cycle-index
design elsewhere in the tab rests on.

> Full measurement, on 650k within-track steps: `docs/MINFLUX_analysis_methods.md` §I-B.0.

#### Continuity score

A tolerance-weighted measure of how much of a track is genuinely observed rather than inferred across
darkness, using the same weighting as the continuity gate:

```math
w_i = 2^{-(n_{\text{cyc},i} - 1)/H}, \qquad \text{score} = \langle w_i \rangle
```

$H$ is the dark-cycle tolerance (half-life, in cycles). An integrated step scores 1 and the weight
halves every $H$ dark cycles, so a dataset of gentle 1–2-cycle gaps scores higher than its bare
integrated fraction would suggest. Reported for raw (frozen at upload) and analysed data side by side,
so the effect of filtering is visible.

---

### MSD Computation

MINFLUX tracking data is irregularly sampled in time, because cycles are skipped when the molecule is
dark. The Tracking tab therefore bins MSD by **integer cycle-count lag**, not wall-clock time.

For a single track $\lbrace (x_i, y_i, z_i, t_i) \rbrace$ sorted by time, all pairs $(i, j)$ with $i < j$ are
considered, and each pair contributes its squared displacement to the bin whose centre $\tau$
satisfies $\tau_\text{lo} \leq t_j - t_i < \tau_\text{hi}$:

```math
\text{MSD}(\tau) = \frac{1}{N_\tau} \sum_{\substack{(i,j)\,:\,\\ \tau_\text{lo} \leq t_j - t_i < \tau_\text{hi}}} \left[(x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2\right]
```

where $N_\tau$ is the number of pairs in that bin.

The **ensemble MSD** across all tracks is computed bin-by-bin as either the **pair-weighted mean**,
which downweights bins with few observations:

```math
\overline{\text{MSD}}(\tau) = \frac{\sum_k N_{\tau,k} \cdot \text{MSD}_k(\tau)}{\sum_k N_{\tau,k}}
```

or the **track-unweighted mean**, a simple average of per-track MSD values that treats each track
equally regardless of how many pairs it contributes.

Bins where fewer than a user-specified fraction of tracks contribute are masked as NA. An optional
D > 0 filter drops tracks whose MSD has a negative OLS slope.

The MINFLUX adaptation here follows **Vogler 2025** [7] and its MinSpt reference implementation: pairs
are grouped by **cumulative cycle index** rather than wall-clock lag,

```math
d_\text{cyc} = \mathrm{round}\!\left(\frac{\Delta t - \text{SW}}{t_c}\right), \qquad \text{cycidx} = \mathrm{cumsum}(d_\text{cyc})
```

so a bin holds pairs separated by the same number of *elapsed cycles*. Binning on
$\mathrm{round}(\text{total}/t_c)$ instead introduces a sawtooth; the cumulative index removes it. The
track-unweighted ensemble and the D > 0 filter are that work's conventions too.

---

### Diffusion Coefficient Fitting

**Brownian diffusion** model:

```math
\text{MSD}(\tau) = 2n D \tau + 2n \sigma_\text{loc}^2
```

where $n$ is the number of dimensions, $D$ the diffusion coefficient (nm²/ms = µm²/s × 10⁻³), and
$\sigma_\text{loc}$ the localization precision (nm). The intercept $2n\sigma_\text{loc}^2$ arises from
localization noise and does not reflect true motion. Note that $n$ here is the **motion degrees of
freedom** set in *Sample geometry*, not the number of recorded axes.

Fitting is two-pass weighted least squares:
1. Unweighted linear regression for initial estimates $D_0$, $\sigma_0$
2. Weighted regression with $w_k \propto N_{\tau,k} / \widehat{\text{MSD}}(\tau_k)^2$, which minimizes
   variance for Poisson-like pair-count noise

The linear fit range comes from **Michalet & Berglund optimal-point selection** (OLSF) rather than a
fixed percentage of the curve, with separate optimal point counts for the slope (→ D) and the intercept
(→ σ²). Its adaptation to MINFLUX cycle lags - the iterative variant that converges onto the lowest,
least-noisy lags, and the default 5–50 cycle fit window - follows **Vogler 2025** [7].

**Anomalous diffusion** model:

```math
\text{MSD}(\tau) = 2n D_\alpha \tau^\alpha + 2n \sigma_\text{loc}^2
```

fitted by non-linear least squares (Levenberg–Marquardt via `minpack.lm`). The anomalous exponent:
- $\alpha < 1$: subdiffusion (confined or obstructed)
- $\alpha = 1$: normal Brownian diffusion
- $\alpha > 1$: superdiffusion (directed transport)

The confined and anomalous fits use the whole lag window, because the optimal-point theory is derived
for free diffusion only. Reported units: $D$ in µm²/s, $\sigma_\text{loc}$ in nm.

---

### CVE (Covariance-Based Estimator)

The Data Quality panel estimates $D$ **without fitting an MSD curve** (Vestergaard et al. 2014), and
this is the app's primary free-D readout, with the MSD panel serving as its cross-check.

$\sigma_\text{loc}$ is **not** taken from the CVE. It comes from the photon budget,
$\sigma_i = L_\text{TCP}/(a\sqrt{N_i})$ with $N$ = ECO and $a$ an adjustable
instrument constant, and that is what becomes the app-wide noise floor consumed by the D-map, JDD, confinement, HMM, directionality and speed tools. The reasons are
below, under *Why σ does not come from the CVE*.

#### Why consecutive steps carry the precision

Write each recorded position as truth plus independent localization error, $x_i = x_i^\ast + \varepsilon_i$
with $\text{Var}(\varepsilon) = \sigma_\text{loc}^2$. Then a step is

```math
\Delta x_i = x_{i+1} - x_i = \underbrace{(x_{i+1}^\ast - x_i^\ast)}_{\text{motion}} + \varepsilon_{i+1} - \varepsilon_i
```

Two **consecutive** steps share the localization at index $i+1$ - it enters $\Delta x_i$ with a $+$ and
$\Delta x_{i+1}$ with a $-$. For Brownian motion the true increments are independent, so every term
vanishes in the product except that shared error, which survives with a negative sign:

```math
\langle \Delta x_i \, \Delta x_{i+1} \rangle = -\langle \varepsilon_{i+1}^2 \rangle = -\sigma_\text{loc}^2
```

That sign is the whole trick. Real motion cannot produce a negative lag-1 covariance, so the two
contributions separate cleanly with no model fitting:

```math
D_\text{CVE} = \frac{\langle \Delta x^2 \rangle}{2\Delta t} + \frac{\text{cov}_1}{\Delta t}
```

The added $\text{cov}_1$ term (negative in value) is what removes the localization-error inflation from the
naive $\langle \Delta x^2 \rangle / 2\Delta t$. It is computed per axis and averaged over x, y and,
when real z is present, z. **D is independent of motion blur** (the identity above holds for any exposure model), so D
survives the correction below unchanged.

**Only consecutive integrated (single-cycle) step pairs are used.** The derivation assumes constant
$\Delta t$; admitting dark-gap steps of varying length would bias the covariance. See TSA above for
what "integrated" means.

#### Why σ does not come from the CVE

The identity $\sigma^2 = -\text{cov}_1$ is the **instantaneous-exposure** case, $R = 0$ in
Vestergaard's notation. A MINFLUX localization is not a snapshot: photons are collected across the TCP
cycle, so the reported position is a time-average of the path over that window and the blur term does
not vanish. The general form is

```math
\sigma_\text{loc}^2 = R\,\langle \Delta x^2 \rangle + (2R - 1)\,\text{cov}_1,
\qquad R = \frac{t_c}{6\,\Delta t}
```

with $R \to 1/6$ for a full-cycle window. Reading MINFLUX data at $R = 0$ under-reports σ roughly
threefold.

Correcting R is necessary but not sufficient. At $R \approx 1/6$ roughly **90% of $\sigma^2$ comes
from the $R\langle\Delta x^2\rangle$ term**, so the estimate is close to
$\sqrt{R} \times \text{rms step}$, a function of the assumed R and of the step-size distribution,
which population heterogeneity and any residual drift both inflate. What the data contributes
independently, $\text{cov}_1$, is a few-percent residual. On an immobile sample the increments can
even come out *positively* correlated, at which point the estimator returns no solution at all.

The limit is identifiability. Within one trajectory, blur and precision enter the covariance in the
same place, so MSD-intercept and MLE fits inherit the same degeneracy.

So σ comes from the photon budget instead:

```math
\sigma_i = \frac{L_\text{TCP}}{a\sqrt{N_i}}\left[1 + 2\left(\frac{x}{L}\right)^2\right],
\qquad x = \sqrt{2 D\, t_\text{loc}}
```

per localization, with $N$ = ECO and $L_\text{TCP}$ from *Acquisition Parameters* (auto-filled from the
sequence's `patGeoFactor`). This is a direct instrument model, evaluated per localization and free of any assumption about the
motion. The bracketed term is the one place motion enters. It is the
distance the emitter covers during one localization, so precision degrades for fast molecules and
tightens under confinement.

The blur-aware CVE σ is still shown beside it as a sanity read on the covariance. It is not the
app-wide value and nothing downstream consumes it.

> **σ ∝ L, or steeper?** The full doughnut CRLB carries a denominator
> $\left|1 + \ln 2\,(L/\nu)^2\left((2x/L)^2 - 1\right)\right|$, with ν the excitation-doughnut FWHM,
> which makes σ rise faster than $L$ as $L$ approaches ν. It is reached through the **ν** field and is
> **off by default**, since $\sigma \propto L$ holds unless an L-series on the instrument shows
> otherwise.

#### Instrument constants, and porting to another microscope

Three values in *Acquisition Parameters* describe the microscope itself. Treat the defaults as
starting points: measure them on the system that produced the data, or state in the methods that the
defaults were used.

| field | what it is | how to measure |
|---|---|---|
| **σ prefactor a** | absorbs background and estimator losses on top of the ideal Cramér-Rao bound | immobile sample, damping off, two or more TCP diameters; fit $\text{MSD}(\tau) = 4\sigma^2 + 4D_w\tau$ on short lags, then $a = L\sqrt{\bar{N}}\,/\,\sigma$ |
| **patGeoFactor base** | converts the sequence's pattern factor to a diameter, $L = \text{patGeoFactor} \times \text{base}$ | check one diameter against a figure the vendor quotes |
| **wander floor $D_w$** | apparent diffusion the tracking loop produces on a stationary target | the slope of the same MSD fit |

Fixed cells labelled as usual are a better calibration sample than dye on glass, because $a$ should
absorb the background of the real optical environment. Tracker wander is sub-diffusive, so a straight
line fitted to short lags lands above the true intercept and inflates σ: the measured $a$ is a lower
bound, and the true value lies between it and 4. Two or more diameters also test whether σ ∝ L holds,
which is what the beam FWHM ν field is for if it does not.

#### The wander floor

The reported position wanders even on a genuinely immobile target, so any measured D carries an
instrumental floor. It is **photon-driven**, shrinking as the photon budget grows, and it worsens at
larger TCP diameters.

Entering it as $D_w$ makes the calibration panel report $D - D_w$ beside D, and warn when the floor
exceeds 20% or 50% of the measured value. It is never subtracted silently: $D_w$ scales with the
photon budget, so one constant would over-correct dim tracks and under-correct bright ones.

A D of the same order as the floor is not resolved above the instrument. Because the floor falls
faster than the photon count rises, the effective lever is to **raise the photon limit**, trading
$D_\text{max}$ headroom for a floor you can get under.

#### What to report in a methods section

- The acquisition parameters from the sequence: TCP diameter of the repeated iteration, photon limit,
  dwell time, damping, pattern and maxOffTime.
- The values measured from the data: cycle time, median ECO and EFO, time to localization, and the
  number of tracks and localizations after cleaning, with the data-selection funnel.
- The three instrument constants above, and how each was obtained.
- D and σ_loc, with the route used for each.

Table 2 of Vogler et al. 2025 is a good template.

#### What the simulator round-trip still covers

The simulator generates tracks with a dialled-in D and then runs the same CVE over them
(`docs/simulation_module.md` §6a). Demoting the CVE for σ_loc leaves that loop in place, and if
anything it now carries more weight, because D is the one quantity the CVE still supplies to the rest
of the app. It remains the check on the parts of the D path that are easy to get wrong and hard to
see: the single-cycle step restriction, the Δt = t_c + SW convention, and the motion-dimensionality
denominator (2 for surface diffusion, 3 for volume) that would otherwise under-report D by a third on
membrane and rod geometries.

The σ_loc half of the same scorecard has narrowed in scope, for two reasons.

Its target moved. Simulated positions are now perturbed with the same photon-budget expression the
Tracking tab applies ($\sigma = L/(a\sqrt{N})$), so recovering σ from those tracks confirms the
arithmetic is consistent between the forward and inverse paths, and confirms nothing about whether
the expression describes a real microscope. Only immobile-sample data can speak to that.

Its conditions are gentler than reality. The simulator draws each position independently around the
true one, with no model of the tracking feedback loop. On the instrument the TCP is recentred on its
own noisy estimate, which turns part of that error into a correlated wander (§ *The wander floor*),
and it is precisely that correlation which defeats covariance-based σ estimation on real data. A
covariance σ estimator can therefore round-trip well in simulation while failing on the microscope,
which is roughly what happened here.

Read a σ_loc round-trip as a unit test on the estimator arithmetic under ideal noise. Read the D
round-trip as the validation of the D path it has always been. Adding a recentring loop to the
simulator would close the gap, and is the change that would make the σ half informative again.

#### How well is σ actually resolved?

The normalized lag-1 correlation $r_1 = \text{cov}_1 / \langle \Delta^2 \rangle$ has, for Brownian
motion plus noise, the closed form

```math
r_1 = \frac{-\sigma_\text{loc}^2}{2D\Delta t + 2\sigma_\text{loc}^2}, \qquad \text{noise fraction} = -2 r_1
```

- $r_1 \to -0.5$: noise dominates the step variance, σ is cleanly resolved
- $r_1 \to 0$: σ is a marginal residual on top of motion, which is the normal regime for MINFLUX
  tracking - per-track σ differences then mostly track mobility, not precision
- $r_1 > 0$: correlated or directed motion dominates and σ is not resolvable there

The panel reports this as a percentage and warns below 15%.

#### Lateral vs axial precision

The pooled σ averages x, y **and** z, and in MINFLUX the axial precision differs from the lateral one,
usually by being tighter. One pooled number then misstates $\sigma_{xy}$. Both are computed separately,
with the ratio $\sigma_z / \sigma_{xy}$ (reported as `sig_aniso`). When the two differ by more than 1.3×
in either direction, the panel warns and gives the lateral value to enter manually: the **D-map is an
x–y map and the HMM never reads z**, so a pooled floor makes them over-subtract noise when z is wider
and under-subtract it when z is tighter. The app-wide default stays pooled.

#### Diagnostics and QC flags

- **Photon–mobility coupling** $\rho$ = Spearman correlation between per-track photons and per-track
  step variance. On MINFLUX *tracking* data a fast molecule leaves the beam pattern sooner and collects
  *fewer* photons, so the photon channel is a mobility proxy rather than a pure brightness channel.
  Strongly negative ρ therefore means **any ECO/EFO filter is silently a mobility filter** - surfaced
  as a warning at $\rho \leq -0.35$.
- **Implausible jumps**: flagged where $|\Delta r|^2 > 25 \cdot E[|\Delta r|^2]$, using
  $E[|\Delta r|^2] = 2dD_\text{CVE}\Delta t + 2d\sigma_\text{loc}^2$ so that long dark gaps are not
  falsely flagged.
- **Degenerate tracks**: $\sigma^2 \leq 0$ (drift-dominated) or σ far in the outlier tail.
- **Per-track heterogeneity**: one CVE σ per track with at least $n_\text{min}$ integrated steps,
  skipping the first ~5 localizations of lock-on transient; reported as median, IQR, CV and a bootstrap
  95% CI, plus a σ–EFO Spearman ρ (precision should improve with photons, so negative ρ validates the
  model).

The estimate is deliberately **unweighted by photons**. Weighting by brightness would up-weight slow
molecules - for the coupling reason above - and bias D downward, so the pooled σ_loc and D_CVE treat
every qualifying step equally.

> **Detecting real z.** 2D exports routinely carry an all-zero `z` column, and `0` is finite, so a bare
> presence test would treat 2D data as 3D and pool a zero-variance axis into the maths. A shared
> predicate - more than 10 finite z values **and** sd > 1 nm - gates the 3D path for CVE, MSD, JDD,
> confinement and displacement-vs-D alike, so the dimension decision cannot drift between tools.

---

### Confinement Index

The confinement index (CI) provides a model-free measure of diffusion character:

```math
\text{CI} = \frac{\text{MSD}(\tau_\text{late}) / \text{MSD}(\tau_\text{early})}{\tau_\text{late} / \tau_\text{early}}
```

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
> datasets with different track lengths** - a longer track has more opportunity to reveal a plateau.
> Compare only at matched track lengths, and read the simulated null band shown alongside.

---

### Jump Distance Distribution (JDD)

The step-length histogram is fitted with an **N-component scaled-chi mixture** - the law of $|\Delta r|$
for a $d$-dimensional isotropic Gaussian displacement:

```math
p(r) = \sum_k w_k \, C_d \, \frac{r^{\,d-1}}{\sigma_k^{\,d}} \exp\!\left(-\frac{r^2}{2\sigma_k^2}\right), \qquad C_d = \frac{1}{2^{\,d/2-1}\,\Gamma(d/2)}
```

$d = 2$ gives the classic **Rayleigh** form, $d = 3$ the **Maxwell** form used when z is included.
$\sigma_k$ is the **per-axis** spread in both cases. Fitted by MLE (Nelder–Mead on the negative
log-likelihood). Per component:

```math
D = \frac{\sigma^2}{2} \quad \text{(dimension-independent)}, \qquad \sigma^2 \leftarrow \sigma^2 - \frac{2\sigma_\text{loc}^2}{\Delta t_0} \ \ (\text{clamped} \geq 0)
```

$D = \sigma^2/2$ carries no $d$ because $E[r^2] = d\sigma^2$ while $\text{MSD} = 2dD\Delta t$ - the $d$
cancels. Steps are pooled as $r/\sqrt{\Delta t}$ by default so different cycle counts are comparable.
A noise-floor jump distance $\sqrt{2d} \sigma_\text{loc}/\sqrt{\Delta t_0}$ marks where peaks are
noise-dominated.

> The dimension $d$ must match the data: fitting the 2D Rayleigh shape to 3D distances inflates the
> fitted σ, so the σ_loc term subtracts only part of the true noise floor and the PDF shape is wrong.
> $d$ follows the same real-z predicate used by CVE and MSD.

---

### Diffusive-State HMM

Per-step displacement $r$ (√Δt-normalised) with a **Rayleigh emission** over $K = 2$–4 states:

```math
\log p(r \mid k) = \log r - 2\log\sigma_k - \frac{r^2}{2\sigma_k^2}, \qquad \sigma_k = \sqrt{\frac{\sum \gamma \, r^2}{2\sum \gamma}}
```

fitted by **Baum–Welch EM** entirely in log-space, then decoded for all tracks with **Viterbi**.

The fitted per-step transition matrix $A$ is *not* a rate: a 3-dark-cycle step counts as one transition
opportunity just like a 1-cycle step. So the **generator** $Q$ is estimated separately from the decoded
path using real elapsed time, giving dark-cycle-invariant dwell times:

```math
q_{ab} = \frac{N_{ab}}{T_a} \ (b \neq a), \qquad q_{aa} = -\sum_{b \neq a} q_{ab}, \qquad \tau_a = -\frac{1}{q_{aa}}
```

where $N_{ab}$ counts decoded $a \to b$ jumps and $T_a$ is time spent in state $a$.

---

### Displacement vs Diffusion

Per track, the net start-to-end displacement $R_\text{net} = |r(T) - r(0)|$ is compared against the
free-diffusion expectation. It depends only on the two endpoints and the true elapsed time $T$, so it
is **inherently robust to irregular Δt** - no cycle binning needed:

```math
R_\text{net}^2 = 2 D_\text{mot} \, T \cdot \chi^2(d_\text{mot}), \qquad R_\text{net}^2 \leftarrow R_\text{net}^2 - 2 d \sigma_\text{loc}^2
```

$d_\text{mot}$ is the degrees of freedom of the **motion**, not the number of recorded coordinates - a
membrane protein tracked in 3D has $d_\text{mot} = 2$ on $d_\text{obs} = 3$ axes. This is what *Sample
geometry* sets.

---

### Directionality / Autocorrelation

**Turn angle** between consecutive displacement vectors, and the **velocity autocorrelation** of unit
step vectors:

```math
\theta = \arccos\!\left(\frac{\Delta r_i \cdot \Delta r_{i+1}}{|\Delta r_i|\,|\Delta r_{i+1}|}\right), \qquad AC(\text{lag}) = \langle \hat{v}_i \cdot \hat{v}_{i+\text{lag}} \rangle
```

Both use only steps spanning at most a few dark cycles - a long-gap step has an unreliable direction.
The autocorrelation lag is **elapsed time** (cycle lag × $t_c$), not step index. A noise floor at RMS
step $= 2\sigma_\text{loc}$ quantifies how much of any backward peak is localization noise rather than
real anti-persistence.

---

### Speed & Transport

Path speed $dr/dt$ cannot separate transport from diffusion - a finely sampled random walk has a high
path speed. Two quantities that can:

```math
v_\text{drift}(t) = \frac{|\text{net displacement over window } W|}{W}, \qquad \alpha(t) = \frac{d \log \text{MSD}}{d \log \tau}
```

Diffusion's net displacement grows only as $\sqrt{t}$, so its apparent $v \sim \sqrt{4D/W}$ *shrinks*
as the window $W$ grows, while directed transport sustains $v_\text{drift}$. Windows are defined in
**time (ms)**, not step count, for robustness to irregular Δt. The local anomalous exponent $\alpha$ is
a rolling log-log slope with the static-noise offset $2d\sigma^2$ subtracted per lag: $\alpha \approx 1$
Brownian, $\alpha \to 2$ ballistic.

A window is called directed when $v_\text{drift}$ exceeds the larger of two thresholds:

```math
v_\text{drift} > \max\!\left(\text{fold} \times \text{diffusive baseline},\ \frac{\sqrt{2}\,\sigma}{W}\right)
```

- measured against the track's own baseline and floored at the localization-noise velocity.

---

### Spatially Resolved D-map

A single-lag estimator per step, midpoint-binned on an $n \times n$ grid:

```math
D_\text{raw} = \frac{r^2}{4\,\Delta t} \quad \text{(2D)}, \qquad D_\text{cell} = \max\!\Big(0,\ \langle D_\text{raw} \rangle - \sigma_\text{loc}^2 \big\langle \tfrac{1}{\Delta t} \big\rangle\Big)
```

Each cell reports the mean with its standard error, because $\langle r^2 \rangle = 4D\Delta t + 4\sigma^2$
holds for the mean only; the median of single-step $r^2$ sits at ln 2 of the mean. Cells with too few
steps are left NA. The noise subtraction removes the spurious diffusion floor that would otherwise make
immobile regions look mobile. This is the **short-lag limit of CVE**, far more robust than a per-cell
MSD fit in small cells such as a bacterium.

---

### Fast-Diffusion Ceiling

The trackable-D limit, before the localization loop loses lock, in the form given by **Vogler 2025**
[7]:

```math
D_\text{max} = \frac{L_\text{TCP}^2}{16} \cdot \frac{\psi}{N_\text{PL} + \psi \, t_\text{hw}}
```

where $L_\text{TCP}$ is the targeted-coordinate-pattern diameter, $\psi$ the photon rate, $N_\text{PL}$
the photons per localization and $t_\text{hw}$ the hardware overhead. In the app $\psi$ is the median
**EFO**, $N_\text{PL}$ the median **ECO** (which *is* the photon count per localization), and
$t_\text{loc} = \text{ECO}/\text{EFO}$ feeds the simpler $L^2/(16\,t_\text{loc})$ fallback used when
either column is missing. Surfaced as a green/amber/red flag.

**Which $L$ to enter (and how to avoid entering it).** Upload the acquisition sequence `.json` under *Time
Signature Analysis* and $L_\text{TCP}$ is read from it directly. Sequence files do not store the pattern
diameter; they store `patGeoFactor`, a dimensionless multiplier, and

```math
L_\text{TCP} = \text{patGeoFactor} \times 360\ \text{nm}
```

wavelength-independently. The panel takes the **final iteration**, the one the sequence reiterates continuously once locked on
(`headstart = -1`), because that is the pattern the particle has to stay inside between consecutive
fixes. Coarse search iterations do not enter, and neither does the much
tighter final iteration an *imaging* sequence uses. SPT stops at a larger $L$ on purpose: precision
improves only linearly as $L$ shrinks ($\sigma_\text{CRLB} \propto L/4\sqrt{N}$) while the risk of the
particle escaping the pattern grows quadratically.

Typical ladders, for orientation:

| Sequence | Iteration diameters (nm) | Final |
|---|---|---|
| Tracking, 2D | 284 · 302 · 151 · 101 | **100.8** |
| Tracking, 3D octahedron | 288 · 288 · 288 · 151 · 101 | **100.8** |
| Imaging, 2D | 288 · 288 · 151 · 76 · 40 | 39.6 |

Without a sequence file the app falls back to 100 nm. Vogler quote 150 nm because their sequence stops
one iteration earlier; an imaging sequence's 39.6 nm is a different regime entirely. Since the ceiling
scales as $L^2$, this is the input that most rewards being read from the file rather than assumed.
$t_\text{hw}$ defaults to 0.05 ms, Vogler's $t_{\text{hw},\eta=1} \approx 52\ \mu\text{s}$ for a
hexagonal TCP.

It is a theoretical form, so it may under-warn near the demonstrated practical ceiling of roughly
2.5–3 µm²/s: Vogler compute $D_\text{max} \approx 10\text{–}11$ µm²/s under optimized parameters but
demonstrate tracking only to $\approx 2.5$ µm²/s on lipid analogues.

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

## Notes & Tips

- **No automatic quality filter**: nothing is removed on load except invalid records (`vld = FALSE`).
  In particular there is **no automatic CFR cutoff** in either tab - CFR is not a reliable artefact
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
- **Outlier trimming shortens traces**: it removes localizations, so a trace can drop below the TID
  length minimum and be removed by that filter. The info box in the spread panel counts this when it
  happens. Trimming also lowers the reported σ̃ partly by selection, so prefer a conservative *k* and
  quote the setting alongside any precision figure.
- **Multiple datasets**: in the SMLM tab several files can be loaded simultaneously and compared with
  the Active file selector. Lock filters when switching to prevent accidental resets, and use Dataset
  Groups in the Analysis tab to pool and compare conditions.
- **Re-chunking is non-destructive**: the original chunk set is frozen at upload time. Every
  Plot Data reads from that frozen original, so resetting all filters and pressing it again recovers
  the full dataset.
- **MSD fit range**: the linear fit uses Michalet & Berglund optimal-point selection for you. If you
  override it, keep to the early lags - later points rest on fewer pairs, have higher variance, and
  deviate from the linear model through confinement or undersampling.

---

## References

The published literature each implemented method follows. Numbers link to the bibliography below.
`docs/MINFLUX_analysis_methods.md` Part II additionally rates every method for **appropriateness on
MINFLUX data specifically** (green / amber / red) - the more important read before trusting a number.

### Method → literature

| Method | Follows |
|---|---|
| [CVE - D](#cve-covariance-based-estimator) | Vestergaard 2014 [12] |
| [σ_loc - photon budget](#cve-covariance-based-estimator) | MINFLUX CRLB for a doughnut TCP; Balzarotti 2017 [1], Vogler 2025 [7] |
| [TSA & cycle time](#time-signature-analysis-tsa--cycle-time) | Vogler 2025 [7] (MFX-TSA / MinSpt) |
| [MSD computation](#msd-computation) | Qian 1991 [22]; Michalet 2010 [10]; **Vogler 2025 [7]** - cycle-index binning and ensemble convention |
| [MSD fit-range selection](#diffusion-coefficient-fitting) | Michalet & Berglund 2012 [11]; **Vogler 2025 [7]** - OLSF as adapted to MINFLUX cycle lags |
| [Anomalous exponent α](#diffusion-coefficient-fitting) | Metzler & Klafter 2000 [26]; Metzler 2014 [19]; Kepten 2015 [20] |
| [Confinement](#confinement-index) | Simson 1995 [17]; Kusumi 1993 [23]; Saxton 1993/1995 [24][25]; Destainville & Salomé 2006 [27] |
| [Jump distance distribution](#jump-distance-distribution-jdd) | Schütz 1997 [13]; Weimann 2013 [14] |
| [Diffusive-state HMM](#diffusive-state-hmm) | Das, Cairo & Coombs 2009 [15]; Persson 2013 [16] (vbSPT) |
| [Directionality / turn angles](#directionality--autocorrelation) | Burov 2013 [18] |
| [Spatially resolved D-map](#spatially-resolved-d-map) | Vestergaard 2014 [12] (short-lag CVE limit) |
| [Fast-diffusion ceiling](#fast-diffusion-ceiling) | Vogler 2025 [7] |
| [Track segmentation (change points)](#track-segmentation--change-point--gmm) | Killick, Fearnhead & Eckley 2012 [21] (PELT) |
| [Nearest-neighbour distances](#nearest-neighbour-distances) | Clark & Evans 1954 [28] |
| [Median Clustering](#median-clustering) | Ester 1996 [29] (DBSCAN); Ankerst 1999 [30] (OPTICS); Campello 2013 [31] (HDBSCAN) |
| [Ripley's L(r)](#ripleys-lr) | Ripley 1977 [32]; Besag 1977 [33]; Dixon 2002 [35] |
| [Bivariate Ripley's L₁₂(r)](#bivariate-ripleys-l₁₂r) | Lotwick & Silverman 1982 [34] |
| [Pair cross-correlation g₁₂(r)](#pair-cross-correlation-g₁₂r) | Sengupta 2011 [36]; Veatch 2012 [37] |
| [Cross-channel NN & co-localization](#cross-channel-nearest-neighbor--co-localization-curve) | Lachmanovich 2003 [40] |
| [KDE density map & Pearson r](#2d-kde-density-map--co-localization) | Silverman 1986 [38] (bandwidth); Dunn 2011 [41]; Costes 2004 [42] |
| [Radius of gyration & anisotropy](#radius-of-gyration) | Rudnick & Gaspari 1986 [43]; Westin 1997 [44] |
| [Cluster shape fitting / PCA axes](#fitting--metrics-3d-only) | Pearson 1901 [45]; Barber 1996 [47] (Qhull) |
| [Reference surface (LOESS)](#global-reference-plane-and-surface) | Cleveland 1979 [46] |
| Blinking / over-counting caveat | Annibale 2011 [39] |
| MINFLUX platform & tracking context | Balzarotti 2017 [1]; Gwosch 2020 [2]; Schmidt 2021 [3]; Ostersehlt 2022 [6]; Scheiderer 2025 [8] |
| MINFLUX tracking applications (what the data can resolve) | Wolff 2023 [4]; Deguchi 2023 [5] |
| Simulator validation context | Marin & Ries 2025 [9] (SimuFLUX) |
| General SPT background | Saxton & Jacobson 1997 [48]; Manzo & Garcia-Parajo 2015 [49]; Shen 2017 [50] |

### Bibliography

**MINFLUX microscopy & tracking**

1. Balzarotti F, Eilers Y, Gwosch KC, Gynnå AH, Westphal V, Stefani FD, Elf J, Hell SW (2017). *Nanometer resolution imaging and tracking of fluorescent molecules with minimal photon fluxes.* Science 355(6325):606–612. doi:10.1126/science.aak9913
2. Gwosch KC, Pape JK, Balzarotti F, Hoess P, Ellenberg J, Ries J, Hell SW (2020). *MINFLUX nanoscopy delivers 3D multicolor nanometer resolution in cells.* Nature Methods 17:217–224. doi:10.1038/s41592-019-0688-0
3. Schmidt R, Weihs T, Wurm CA, Jansen I, Rehman J, Sahl SJ, Hell SW (2021). *MINFLUX nanometer-scale 3D imaging and microsecond-range tracking on a common fluorescence microscope.* Nature Communications 12:1478. doi:10.1038/s41467-021-21652-z
4. Wolff JO, Scheiderer L, Engelhardt T, Engelhardt J, Matthias J, Hell SW (2023). *MINFLUX dissects the unimpeded walking of kinesin-1.* Science 379(6636):1004–1010. doi:10.1126/science.ade2650
5. Deguchi T, Iwanski MK, Schentarra EM, et al. (2023). *Direct observation of motor protein stepping in living cells using MINFLUX.* Science 379(6636):1010–1015. doi:10.1126/science.ade2676
6. Ostersehlt LM, Jans DC, Wittek A, et al. (2022). *DNA-PAINT MINFLUX nanoscopy.* Nature Methods 19:1072–1075. doi:10.1038/s41592-022-01577-1
7. Vogler BTL, De Angelis G, Zhao Z, Eggeling C, Reina F (2025). *Parameter optimization for MINFLUX microscopy enabled single-particle tracking.* Communications Biology 8:1573. doi:10.1038/s42003-025-09060-1 - reference implementation *Toolbox_for_MINFLUX_SPT* (MinSpt), the source of the Time-Signature Analysis and the cycle-based MSD workflow ported here.
8. Scheiderer L, Marin Z, Ries J (2025). *MINFLUX achieves molecular resolution with minimal photons.* Nature Photonics 19(3):238–247. doi:10.1038/s41566-025-01625-0
9. Marin Z, Ries J (2025). *Evaluating MINFLUX experimental performance in silico* (SimuFLUX). Nature Communications 17:246. doi:10.1038/s41467-025-66952-w

**Single-particle-tracking methods**

10. Michalet X (2010). *Mean square displacement analysis of single-particle trajectories with localization error.* Phys. Rev. E 82:041914. doi:10.1103/PhysRevE.82.041914
11. Michalet X, Berglund AJ (2012). *Optimal diffusion coefficient estimation in single-particle tracking.* Phys. Rev. E 85:061916. doi:10.1103/PhysRevE.85.061916
12. Vestergaard CL, Blainey PC, Flyvbjerg H (2014). *Optimal estimation of diffusion coefficients from single-particle trajectories.* Phys. Rev. E 89:022726. doi:10.1103/PhysRevE.89.022726
13. Schütz GJ, Schindler H, Schmidt T (1997). *Single-molecule microscopy on model membranes reveals anomalous diffusion.* Biophys. J. 73:1073–1080.
14. Weimann L, Ganzinger KA, McColl J, et al. (2013). *A quantitative comparison of single-dye tracking analysis tools using Monte Carlo simulations.* PLOS ONE 8:e64287. doi:10.1371/journal.pone.0064287
15. Das R, Cairo CW, Coombs D (2009). *A hidden Markov model for single particle tracks quantifies dynamic interactions between LFA-1 and the actin cytoskeleton.* PLoS Comput. Biol. 5:e1000556.
16. Persson F, Lindén M, Unoson C, Elf J (2013). *Extracting intracellular diffusive states and transition rates from single-molecule tracking data* (vbSPT). Nature Methods 10:265–269.
17. Simson R, Sheets ED, Jacobson K (1995). *Detection of temporary lateral confinement of membrane proteins using single-particle tracking analysis.* Biophys. J. 69:989–993.
18. Burov S, Tabei SMA, Huynh T, et al. (2013). *Distribution of directional change as a signature of complex dynamics.* PNAS 110(49):19689–19694. doi:10.1073/pnas.1319473110
19. Metzler R, Jeon JH, Cherstvy AG, Barkai E (2014). *Anomalous diffusion models and their properties.* Phys. Chem. Chem. Phys. 16:24128–24164. doi:10.1039/c4cp03465a
20. Kepten E, Weron A, Sikora G, Burnecki K, Garini Y (2015). *Guidelines for the fitting of anomalous diffusion mean square displacement graphs.* PLOS ONE 10:e0117722. doi:10.1371/journal.pone.0117722
21. Killick R, Fearnhead P, Eckley IA (2012). *Optimal detection of changepoints with a linear computational cost* (PELT). J. Am. Stat. Assoc. 107(500):1590–1598. doi:10.1080/01621459.2012.737745
22. Qian H, Sheetz MP, Elson EL (1991). *Single particle tracking. Analysis of diffusion and flow in two-dimensional systems.* Biophys. J. 60:910–921. doi:10.1016/S0006-3495(91)82125-7
23. Kusumi A, Sako Y, Yamamoto M (1993). *Confined lateral diffusion of membrane receptors as studied by single particle tracking.* Biophys. J. 65:2021–2040. doi:10.1016/S0006-3495(93)81253-0
24. Saxton MJ (1993). *Lateral diffusion in an archipelago. Single-particle diffusion.* Biophys. J. 64:1766–1780.
25. Saxton MJ (1995). *Single-particle tracking: effects of corrals.* Biophys. J. 69:389–398. doi:10.1016/S0006-3495(95)79911-8
26. Metzler R, Klafter J (2000). *The random walk's guide to anomalous diffusion: a fractional dynamics approach.* Phys. Rep. 339:1–77. doi:10.1016/S0370-1573(00)00070-3
27. Destainville N, Salomé L (2006). *Quantification and correction of systematic errors due to detector time-averaging in single-molecule tracking experiments.* Biophys. J. 90:L17–L19. doi:10.1529/biophysj.105.075176

**SMLM & spatial statistics**

28. Clark PJ, Evans FC (1954). *Distance to nearest neighbor as a measure of spatial relationships in populations.* Ecology 35:445–453. doi:10.2307/1931034
29. Ester M, Kriegel HP, Sander J, Xu X (1996). *A density-based algorithm for discovering clusters in large spatial databases with noise* (DBSCAN). KDD-96:226–231.
30. Ankerst M, Breunig MM, Kriegel HP, Sander J (1999). *OPTICS: ordering points to identify the clustering structure.* ACM SIGMOD Record 28(2):49–60. doi:10.1145/304181.304187
31. Campello RJGB, Moulavi D, Sander J (2013). *Density-based clustering based on hierarchical density estimates* (HDBSCAN). PAKDD 2013, LNCS 7819:160–172. doi:10.1007/978-3-642-37456-2_14
32. Ripley BD (1977). *Modelling spatial patterns.* J. R. Stat. Soc. B 39(2):172–212.
33. Besag J (1977). *Contribution to the discussion of Dr Ripley's paper* (the L-transform of K). J. R. Stat. Soc. B 39:193–195.
34. Lotwick HW, Silverman BW (1982). *Methods for analysing spatial processes of several types of points.* J. R. Stat. Soc. B 44(3):406–413.
35. Dixon PM (2002). *Ripley's K function.* In: Encyclopedia of Environmetrics, Wiley. doi:10.1002/9780470057339.var046
36. Sengupta P, Jovanovic-Talisman T, Skoko D, et al. (2011). *Probing protein heterogeneity in the plasma membrane using PALM and pair correlation analysis.* Nature Methods 8:969–975. doi:10.1038/nmeth.1704
37. Veatch SL, Machta BB, Shelby SA, et al. (2012). *Correlation functions quantify super-resolution images and estimate apparent clustering due to over-counting.* PLoS ONE 7:e31457. doi:10.1371/journal.pone.0031457
38. Silverman BW (1986). *Density Estimation for Statistics and Data Analysis.* Chapman & Hall.
39. Annibale P, Vanni S, Scarselli M, Rothlisberger U, Radenovic A (2011). *Identification of clustering artifacts in photoactivated localization microscopy.* Nature Methods 8:527–528. doi:10.1038/nmeth.1627
40. Lachmanovich E, Shvartsman DE, Malka Y, Botvin C, Henis YI, Weiss AM (2003). *Co-localization analysis of complex formation among membrane proteins by computerized fluorescence microscopy.* J. Microsc. 212(Pt 2):122–131. doi:10.1046/j.1365-2818.2003.01239.x
41. Dunn KW, Kamocka MM, McDonald JH (2011). *A practical guide to evaluating colocalization in biological microscopy.* Am. J. Physiol. Cell Physiol. 300:C723–C742. doi:10.1152/ajpcell.00462.2010
42. Costes SV, Daelemans D, Cho EH, Dobbin Z, Pavlakis G, Lockett S (2004). *Automatic and quantitative measurement of protein-protein colocalization in live cells.* Biophys. J. 86(6):3993–4003. doi:10.1529/biophysj.103.038422
43. Rudnick J, Gaspari G (1986). *The asphericity of random walks.* J. Phys. A: Math. Gen. 19:L191–L193. doi:10.1088/0305-4470/19/4/004
44. Westin CF, Peled S, Gudbjartsson H, Kikinis R, Jolesz FA (1997). *Geometrical diffusion measures for MRI from tensor basis analysis.* ISMRM 1997. See also Westin et al. (2002), Med. Image Anal. 6:93–108.
45. Pearson K (1901). *On lines and planes of closest fit to systems of points in space* (PCA). Philosophical Magazine 2(11):559–572. doi:10.1080/14786440109462720
46. Cleveland WS (1979). *Robust locally weighted regression and smoothing scatterplots* (loess). J. Am. Stat. Assoc. 74(368):829–836. doi:10.1080/01621459.1979.10481038
47. Barber CB, Dobkin DP, Huhdanpaa H (1996). *The quickhull algorithm for convex hulls.* ACM Trans. Math. Softw. 22(4):469–483. doi:10.1145/235815.235821

**Reviews**

48. Saxton MJ, Jacobson K (1997). *Single-particle tracking: applications to membrane dynamics.* Annu. Rev. Biophys. Biomol. Struct. 26:373–399. doi:10.1146/annurev.biophys.26.1.373
49. Manzo C, Garcia-Parajo MF (2015). *A review of progress in single particle tracking: from methods to biophysical insights.* Rep. Prog. Phys. 78:124601. doi:10.1088/0034-4885/78/12/124601
50. Shen H, Tauzin LJ, Baiyasi R, et al. (2017). *Single particle tracking: from theory to biophysical applications.* Chem. Rev. 117(11):7331–7376. doi:10.1021/acs.chemrev.6b00815
