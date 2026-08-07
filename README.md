# MINFLUXR

An interactive R/Shiny application for visualization and analysis of **Abberior MINFLUX** data.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Launching the App](#launching-the-app)
4. [MINFLUX Data Columns Reference](#minflux-data-columns-reference)
5. [SMLM Analysis Tab](#smlm-analysis-tab)
   - [Loading Data](#loading-data)
   - [Filtering](#filtering)
   - [Visualization](#visualization)
   - [2-Color Analysis](#2-color-analysis)
   - [DBSCAN Clustering](#dbscan-clustering)
   - [Confocal Overlay & Auto-Alignment](#confocal-overlay--auto-alignment)
   - [Analysis Outputs](#analysis-outputs)
   - [Export](#export)
6. [Tracking Analysis Tab](#tracking-analysis-tab)
   - [Loading & Chunking](#loading--chunking)
   - [Filter Pipeline](#filter-pipeline)
   - [Playback](#playback)
   - [MSD Analysis](#msd-analysis)
   - [Export](#export-1)
7. [MINFLUX Simulator Tab](#minflux-simulator-tab)
8. [Key Calculations](#key-calculations)
   - [Unit Conversions & JSON Import](#unit-conversions--json-import)
   - [Dwell Time](#dwell-time)
   - [Trace Spread (SD Filter)](#trace-spread-sd-filter)
   - [DCR-Based 2-Color Separation](#dcr-based-2-color-separation)
   - [Spatial Density](#spatial-density)
   - [Ripley's L(r)](#ripleys-lr)
   - [Radius of Gyration](#radius-of-gyration)
   - [MSD Computation](#msd-computation)
   - [Diffusion Coefficient Fitting](#diffusion-coefficient-fitting)
   - [Confinement Index](#confinement-index)
9. [Plugin System](#plugin-system)
10. [Notes & Tips](#notes--tips)

---

## Overview

MINFLUXR provides three main analysis environments:

| Tab | Purpose |
|-----|---------|
| **SMLM Analysis** | Static super-resolution data: filtering, clustering, co-localization, spatial statistics |
| **Tracking Analysis** | Single-molecule tracking: trajectory visualization, MSD, diffusion analysis |
| **MINFLUX Simulator** | Generate synthetic MINFLUX datasets from geometric shapes or confocal images |

The app is designed around Abberior MINFLUX JSON exports and also accepts CSV files in pyMINFLUX format or any column-based CSV with x/y coordinates.

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

## MINFLUX Data Columns Reference

The following columns are read from Abberior JSON exports. Understanding them is useful for interpreting the filter controls.

| Column | Full name | Typical units | Description |
|--------|-----------|---------------|-------------|
| `x`, `y`, `z` | Position | nm | Localization coordinates |
| `tid` | Trace ID | integer | Groups localizations belonging to the same fluorophore event |
| `itr` | Iteration | integer | MINFLUX iteration step; only the final valid iteration is kept by default |
| `tim` | Time | s (on disk), ms (tracking tab) | Absolute acquisition timestamp |
| `efo` | Effective Frequency at Offset | Hz | Photon emission rate at circumference of the beampattern (like eco, but in frequency Hz) EFO = ECO/(n*t_dwell), with n being the number of TCP cycles |
| `eco` | Effective Counts at Offset | photons | Summarized photon count at circumference of the beampattern of all TCP cycles |
| `cfr` | Center Frequency Ratio | 0–1 | Ratio of Effective Frequency at Center (EFC) versus Effective Frequency at Offset (EFO). The ratio of emitters from center to periphery of the beampattern. |
| `dcr` | Detector channel ratio | 0–1 | Fraction of emitters being detected by 2 open detectors (i.e. Cy5 near and Cy5 far; used for 2-color separation in 2-color-MINFLUX |
| `dwell` | Dwell time | steps | Derived: number of MINFLUX steps per localization (see [Dwell Time](#dwell-time)) |
| `fluo` | Fluorophore channel | 1 or 2 | fluo groups assigned using the pyMINFLUX pipeline for 2-color MINFLUX |
| `gri` | Grid index | integer | Identifies which scan grid the localization was acquired on |
| `bot` / `eot` / `fnl` | Begin/End/Final of trace flags | logical | Trajectory boundary markers |

---

## SMLM Analysis Tab

### Loading Data

An upload widget browser supporting single, or multiple datasets. Supported formats:

- **`.json`** — Abberior MINFLUX native export. The app reads the full multi-iteration dataset, keeps only the final (highest `itr`) localization per trace, and carries the last non-zero CFR value forward per TID (see [JSON Import](#unit-conversions--json-import)).
- **`.csv`** — pyMINFLUX-format CSV or any flat table with x/y columns. Column names are canonicalized automatically (e.g. `x_nm`, `pos_x`, `xPos` all map to `x`).

Multiple files can be loaded simultaneously; use the **Active file** selector to switch between them. The **Lock filters** checkbox preserves current filter values when switching files.

### Filtering

Filters are applied sequentially in the order listed in the sidebar. All filters operate on the currently active dataset.

Raw data
   └── [1] Metadata filters  (EFO, CFR, DCR, TIM, ECO, dwell, fluo)
        └── [2] TID length filter  (number of points per trace)
             └── [3] TID spread filter  (per trace spread of localizations)


**Metadata filters** (Numeric range inputs):
- `EFO` — Effective Frequency at Offset. The averaged emission frequency over all points in the scan pattern except for the center; exceedingly high EFO values may indicate multi-emitters. 
- `CFR` — Center Frequency Ratio; the ratio of EFC (effective frequency at center) to EFO (effective frequency at offset), i.e. EFC/EFO. A common hard cutoff is CFR ≤ 0.5 to exclude multi-emitter events. For tracking data, a CFR ≤ 0.5 cutoff is applied automatically during chunking in the Tracking tab.

- `TIM` — acquisition timestamp; useful for selecting a time window
- `DCR` — Detector Channel Ratio; ratio of photons arriving at detector 1 vs detector 2, used to separate two fluorescent labels in 2-colour MINFLUX
- `ECO` — Total summarized photon count
- `Dwell` — number of MINFLUX steps
- `FLUO` — fluorophore channel index

**TID length filter**: keeps only traces with a number of localizations within the specified range. Useful for excluding very short (potentially noise) or very long traces (potentially fiducials or moving emitters).

**TID spread filter**: per trace x-, y-, or z-axis spread. See [Trace Spread](#trace-spread-sd-filter) for the math. Setting the upper bound filters out traces with large spatial scatter — useful for excluding gold bead detections, or freely diffusing molecules in a SMLM context.

### Visualization

The **Display** panel controls the scatter plot:

- **2D / 3D** toggle - switches between a 2D top-down view and an interactive 3D Plotly scene. 2D supports XY, XZ and YZ viewmodes (if z exists) to inspect 3D data in 2D modes.
- **Color by** - color localizations by unique trace ID or clusters (default), or with a colorscale for trace ID density, Z (depth), or time (tim) column.
- **Show medians** - plots one point per TID (or per cluster if clustering is active) at the median x/y/z position
- **Scale bar** - adds a physical scale bar to the 2D view; length and thickness are user-controlled
- **Dark mode** - inverts the plot background from white to black
- **Hide axis and gridlines** - hides gridlines and axis lines for a cleaner view.  
- **3D auto-rotate & recording** - enables a spinning animation and records it as a `.webm` video file

### 2-Color Analysis

When data contains two spectrally distinct labels recorded on two detectors, the DCR column (or the `fluo` column if present) can be used to separate them.

**Manual mode**: the user sets two non-overlapping DCR ranges. Localizations are assigned to group 1 or group 2 based on the median DCR of their parent TID.

**Auto-detect DCR split**: fits a two-component Gaussian mixture model (GMM) to the distribution of per-TID median DCR values and finds the optimal threshold automatically. See [DCR-Based 2-Color Separation](#dcr-based-2-color-separation).

Once groups are defined, the **nearest-neighbor (NN) distance** panel computes cross-group, within group, or all, distances. 

<<<<<<< HEAD
### De-blink correction and Median Clustering

Two independent grouping steps, answering different questions. **De-blink correction** is an instrument
artifact fix: one fluorophore that blinks several times produces several trace IDs, and this merges them
back into one molecule. **Median clustering** is biology: it groups distinct molecules into the structure
they form. Only the second is DBSCAN.

#### De-blink correction (Round 1)

Works at **trace level** — a whole trace joins a molecule or it does not, so a trace is never split.
Two criteria:

- **Trace centroid distance** (default) — merge when two traces' median positions are closer than the
  threshold. Set it in nm (`Centroid distance threshold`) or as **k × the measured trace spread**;
  changing k rewrites the distance.
- **Majority vote of localizations** — merge when a fraction *f* of a trace's localizations sit closer to
  a neighbouring trace's centroid than to their own. No distance in nm. Useful values of *f* are well
  below 0.5 (two traces of the *same* molecule score ≈0.5 by construction).

**De-blinking only merges — it never deletes.** A trace that cannot be merged becomes its own
single-blink molecule, so every trace that survived filtering is still counted. Removing data is the job
of the metadata and spread filters, where it is visible and reversible. (There is consequently no
"exclude unclustered" option; cluster 0 applies only to Median Clustering below.)

**Choosing the threshold.** *Measure threshold from data* sets it from your own file: it measures the
per-axis trace spread and finds the smallest k at which no unmerged re-blinks remain. Two diagnostics sit
beside it — a density plot whose orange curve should show two humps (re-blinks left, separate molecules
right) with your threshold in the gap, and an exploratory molecule-count-vs-threshold curve.

**Anisotropic option** — weights x, y and z by their own measured spread, so one threshold fits all three
axes (the merge volume becomes an ellipsoid in nm, tighter along the better-resolved axis). Affects the
merge decision only; reported distances and sizes stay in true nm.

De-blinked molecules are labelled **De-blinked ID** in hover text, and their medians (not individual TID
medians) become the representative points for downstream analysis including Median Clustering.

*Note:* raw-point DBSCAN was the previous method here and has been retired. It merged whenever any two
localizations touched, so a single stray point could fuse two molecules. Against ground-truth simulation
its error count was 35 (splits + merges) versus 10 for centroid distance and 29 for majority vote.
=======
### DBSCAN Clustering

Two independent rounds of DBSCAN clustering are available. Both follow the standard DBSCAN definition: a point is a *core point* if it has at least `minPts` neighbors within distance `eps`; all density-reachable points from a core point form a cluster; points with fewer than `minPts` neighbors are labeled noise (cluster 0).

#### Primary DBSCAN (Round 1)

Runs on **raw localizations** (all points after filtering). Its purpose is to group localizations from the same or closely co-located molecules, particularly useful when one molecule produces multiple nearby trace IDs due to repeated blinking events.

- `eps` — neighborhood radius in nm
- `minPts` — minimum core-point occupancy
- Once active, the app replaces per-TID coloring with per-cluster coloring and uses **primary cluster centroids** (rather than individual TID medians) as the representative points for all downstream analysis including Cluster2.
- **Exclude unclustered** — optionally removes all noise points (cluster = 0) from the analysis and plots.
>>>>>>> c4383f3519c4268414a76a4243a2691a25c5543b

#### Secondary DBSCAN — Cluster2

Runs on **representative points** rather than raw localizations: one point per TID (median x/y/z), or one point per primary cluster (centroid) if primary DBSCAN is active. This second pass groups those representative points into larger spatial structures — for example, a collection of receptor clusters arranged in a ring, or a row of nuclear pores along the nuclear envelope.

- `eps2` - neighborhood radius for the second pass (nm; default 200)
- `minPts2` - minimum representative-point count for a Cluster2 core point (default 15)

**Group scope** (in the 2-color MINFLUX branch): controls which analysis groups are clustered.
- *All* - clusters all representative points together regardless of group
- *Group 1 only* / *Group 2 only* - restricts Cluster2 to one color channel
- *Split* - runs Cluster2 independently on each group with globally unique IDs, so group-1 clusters and group-2 clusters never merge

---

**Visualization**

In **2D mode**:
- *Convex hull* - draws a convex-hull outline around each Cluster2 in the scatter plot

In **3D mode**:
- *Alpha shape* - wraps each Cluster2 in a 3D alpha shape (concave hull via `alphashape3d`)
- *Isosurface* - renders a Gaussian-smoothed density isosurface around each Cluster2, with controls for:
  - Grid resolution (12–40; higher = more detail, slower)
  - Padding around the cluster boundary (0–500 nm)
  - Surface level - quantile of the smoothed density at which the surface is drawn (0.50–0.98; higher = tighter fit)
  - Opacity (0–1)
  - Smoothness sigma (10–500 nm; larger sigma = smoother, more featureless surface)
  - Option to show a selected single Cluster2 ID, or all; single-color or per-cluster color

---

**Fitting & Metrics** (3D only; enabled via "Fitting & Metrics" checkbox)

Each Cluster2's representative points are fitted to a geometric model independently. The **representative point coordinates** used for fitting are per-TID medians (or per-primary-cluster centroids). Two fitting methods are available:

| Method | Description |
|--------|-------------|
| **Heuristic** | Analytical: geometry derived directly from PCA axes and 5th–95th percentile extents. No iteration, no centre optimisation. Fast and always converges. |
| **Weighted orthogonal distance** | Iterative: minimises the sum of squared Euclidean distances from each representative point to the nearest point on the model surface. Equivalent to an orthogonal regression on the geometry. |

The **weighted orthogonal distance** method has additional options:
- *Weight by per-trace precision* - each representative point is weighted by $n / \sigma_\text{spread}^2$ (proportional to $1/\text{SEM}^2$ of the median position). Traces with a single localization (undefined spread) receive the dataset-median spread as a fallback.
- *Centre handling*:
  - **Fixed** - centre = convex-hull centroid; radius/semi-axes solved analytically (no iteration)
  - **Regularized (MAP)** - centre is free but penalized by a Gaussian prior anchored to the robust centre; strength $\lambda$ controls the trade-off. At $\lambda = 0.1$, a 5 nm centre shift adds ~2.5 nm² to the mean squared residual
  - **Free** - centre fully optimized within the data bounding box; use when the centre is expected to deviate significantly from the data centroid

**Fit shapes** available:

| Shape | What is fitted | Key outputs |
|-------|---------------|-------------|
| **Sphere** | Sphere of radius R centred at the hull centroid | Centre (x, y, z), radius R, RMSD |
| **Ellipsoid** | Triaxial ellipsoid aligned to PCA axes | Centre, semi-axes a/b/c, axis orientations, RMSD |
| **Ring** | Flat circle (3D); PC3 = ring-plane normal | Centre, normal vector, radius, RMSD |
| **Arc** | Partial ring; angular extent auto-detected or set via hint | Centre, normal, radius, arc span (°), RMSD |
| **Helix** | Right-handed helix; PCA axis = helix axis | Centre, axis direction, helix radius, pitch (nm/turn), RMSD |

For **Ring**, **Arc**, and **Helix**, a *Render as tube* option draws the fitted curve as a 3D cylinder of user-defined radius (nm), which is useful for visual comparison against the data.

- *Arc span hint* (°) — sets the initial angular extent for arc fitting (0 = auto-detected from the angular spread of the data). The optimiser refines this value.
- *Helix turns hint* — seeds the pitch estimate only; the fitted turn count is determined by the data extent and the optimised pitch.

**RMSD comparison** (*Compute all-model RMSD*, slower): fits all available shape models to each Cluster2 and reports all RMSDs side-by-side in the metrics table, useful for model selection.

---

**Additional 3D overlays**

- *PCA axes* - draws the three principal-component vectors of each Cluster2 as coloured arrows, scaled by their extent
- *Cluster2 summary labels* - text overlay reporting key metrics (centre coordinates, fitted radius, RMSD) per cluster directly on the 3D plot
- *Show distances to fit* - draws line segments from each representative point to its nearest point on the fitted surface; a maximum line length filter (nm) can suppress long outlier lines
- *Fit surface opacity* - controls transparency of the fitted mesh

---

**Membrane protein options** (requires global reference plane to be active)

These options relate each Cluster2 to a global PCA-fitted reference surface (e.g. the plasma membrane):

- *Angle to plane* - angle between the primary PCA axis of the Cluster2 and the reference plane normal; 0° = perpendicular to membrane, 90° = lying flat
- *Distance to plane* - signed distance from the Cluster2 centroid to the reference plane (positive = above, negative = below)
- *Bending radius* - estimates the curvature radius of an arc-like cluster using:
  $$R_\text{bend} = \frac{\ell^2 + s^2}{2s}$$
  where $\ell$ is the PC1 half-extent (chord half-length) and $s$ is the PC2 extent (sagitta / bow height). **This formula is only valid for arc-like clusters viewed roughly edge-on**; it is not meaningful for complete rings, spheres, or helices.

---

**NN distances within Cluster2's**

Nearest-neighbor distances between the representative points inside each Cluster2 are computed and visualized as lines:
- `k` - number of nearest neighbors to find per point (default 3)
- Distance range filter - limits displayed lines to a user-defined nm range
- Optional numeric labels on each line showing the distance in nm
- Exportable as a CSV table

**NN distances between Cluster2's**

Nearest-neighbor distances between the cluster2 centroid between other cluster2 centroids. Is computed and visualized as lines:

- `k` - number of nearest neighbors to find per point (default 3)
- Distance range filter - limits displayed lines to a user-defined nm range
- Optional numeric labels on each line showing the distance in nm
- Exportable as a CSV table


### Confocal Overlay & Auto-Alignment

A confocal image (PNG, TIFF, etc.) can be overlaid on the MINFLUX scatter plot as a reference. Two auto-alignment modes are available once an image and pixel size are provided:

**Foci mode**:
1. Bright blobs are detected in the confocal image by Gaussian smoothing followed by iterative peak finding with a minimum peak-separation constraint.
2. Dense clusters in the MINFLUX data are detected via a 2D density histogram with Gaussian smoothing.
3. Each possible assignment of MINFLUX clusters to confocal peaks is scored by the sum of squared residuals of the implied offsets (geometric consistency). The assignment with the lowest score and a physically valid bounding box is used as the alignment.
4. If no assignment satisfies the bounds constraint, a rank-weighted Hough voting fallback is used.

**Boundary mode**:
1. The bright region of the confocal image is thresholded at the 50th percentile.
2. The intensity-weighted centroid of the bright region is aligned to the centroid of the MINFLUX data.
3. The convex hull of bright pixels is extracted as a cell boundary polygon in data coordinates.
4. Optionally, MINFLUX traces whose median position falls outside this polygon are removed.

### Analysis Outputs

<<<<<<< HEAD
The **Analysis** tab provides sixteen panels. Each produces an interactive plot and, where a table
makes sense, an exportable CSV. Panels marked ⚠ carry an *exploratory* badge in the app — they are
research tools whose interpretation is not yet settled, not validated readouts.

| Panel | What it answers |
|-------|-----------------|
| Dataset Groups | Define named file groups so distribution analyses relabel and pool by group |
| Localization Precision Cloud | Every localization minus its own trace centroid — the empirical precision σ_loc, with covariance ellipses whose tilt exposes drift |
| Fiducial / Contaminant Scan | Flags traces behaving like gold beads or aggregates rather than single molecules (report only) |
| Trace Detection over Time | Detection rate across the acquisition — photobleaching and drop-off |
| Batch NN Distance Histogram | Nearest-neighbour distances across files or channels |
| Ripley's L(r) | Clustering vs randomness across length scales |
| Trace IDs / De-blinked Molecules per Median Cluster | Occupancy — how many traces, or how many de-blinked molecules, per structure |
| Median Clusters Detection over Time ⚠ | When structures appear during acquisition |
| Median Clusters Crowding vs Detection-Time Shift ⚠ | Whether crowding biases when a structure is detected |
| Within Median Clusters NN vs Cluster Size ⚠ | Internal spacing as a function of structure size |
| Cluster Viewer and Shape Fitting ⚠ | Per-structure 3D inspection and geometric fits |
| Radius of Gyration (RoG) | Structure size from the second moment, with anisotropy |
| Cross-Channel NN & Co-localization | Two-colour proximity |
| Bivariate Ripley's L₁₂(r) | Cross-channel clustering across scales |
| Pair Cross-Correlation g₁₂(r) | Cross-channel correlation as a function of separation |
| 2D KDE Density Map | Smoothed density field |

 All analyses operate on the filtered, clustered data visible in the Visualization tab. In 2-color mode, group membership is determined by DCR thresholding (or the `fluo` column if present).
=======
The **Analysis** tab provides eight quantitative tools. Each produces an interactive plot and an exportable CSV table. All analyses operate on the filtered, clustered data visible in the Visualization tab. In 2-color mode, group membership is determined by DCR thresholding (or the `fluo` column if present).
>>>>>>> c4383f3519c4268414a76a4243a2691a25c5543b

---

#### 1. Batch Nearest-Neighbor (NN) Distance Histogram

For each representative point, computes the distance to its $k$-th nearest neighbor within the same scope. The **Unit** selector determines what counts as a representative point:

| Unit | Point definition |
|------|-----------------|
| `tid` | Per-TID median x/y/z position |
| `pcid` | Per-primary-cluster centroid |
| `cluster2_tid` | Per-TID median, restricted to within each Cluster2 |
| `cluster2_pcid` | Per-primary-cluster centroid, restricted to within each Cluster2 |

The $k$-th NN distance for point $i$ is:

$$d_i^{(k)} = k\text{-th smallest} \;\|\mathbf{p}_i - \mathbf{p}_j\|, \quad j \neq i$$

Computed with a kd-tree (`RANN::nn2`). The histogram is summarized by the median and by either the **IQR** or **MAD** as a spread indicator. In 2-color mode, cross-group distances (G1→G2, G2→G1, and within-group) can be selected independently.

> The first-NN distance distribution reveals the typical inter-molecule or inter-cluster spacing. A narrow peak indicates a regular lattice-like arrangement; a broad tail suggests spatial heterogeneity.

---

#### 2. Traces per Cluster (TID Count Histogram)

Counts the number of unique TIDs per Cluster2 (or per primary cluster when Cluster2 is not active) and displays the distribution as a histogram.

> **Caution**: in SMLM experiments a single fluorophore can produce multiple distinct trace IDs across sequential blinking events. The TID count therefore reflects the number of **blinking events**, not necessarily the number of molecules. Experiments with irreversible photoactivation (e.g. PALM) are less affected by this ambiguity.

---

#### 3. Radius of Gyration (RoG)

Quantifies the spatial compactness of each Cluster2. **MINFLUXR computes RoG from per-TID median positions, not from individual localizations**, reducing the influence of localization noise and treating each blinking event as one representative point.

For a Cluster2 containing $m$ TIDs with median positions $(\tilde{x}_k,\tilde{y}_k,\tilde{z}_k)$:

**Step 1** — cluster centroid (mean of TID medians):

$$\bar{x} = \frac{1}{m}\sum_{k=1}^m \tilde{x}_k, \qquad \bar{y} = \frac{1}{m}\sum_{k=1}^m \tilde{y}_k, \qquad \bar{z} = \frac{1}{m}\sum_{k=1}^m \tilde{z}_k$$

**Step 2** — RoG is the RMS displacement of TID medians from the centroid:

$$R_g = \sqrt{\frac{1}{m}\sum_{k=1}^m \!\left[(\tilde{x}_k - \bar{x})^2 + (\tilde{y}_k - \bar{y})^2 + (\tilde{z}_k - \bar{z})^2\right]}$$

Only Cluster2s with $m \geq 2$ TIDs are included (a single point has $R_g = 0$ by definition). The exportable table reports cluster ID, TID count, centroid coordinates, and $R_g$.

---

#### 4. Ripley's L(r)

A spatial clustering statistic that tests whether points are more aggregated than expected under complete spatial randomness (CSR). MINFLUXR uses a **naive (uncorrected) estimator** — no edge correction is applied, which is fast but introduces boundary bias for points near the data extent.

**Naive K estimator**:

$$\hat{K}(r) = \frac{A \cdot 2}{n(n-1)} \sum_{i < j} \mathbf{1}\!\left(\|\mathbf{p}_i - \mathbf{p}_j\| \leq r\right)$$

where $A$ is the bounding area of the convex hull (or bounding volume $V$ in 3D) and $n$ is the number of points. The factor of 2 accounts for the symmetric count of unordered pairs.

**L-transformation**:

$$\hat{L}(r) = \sqrt{\frac{\hat{K}(r)}{\pi}} \;\text{(2D)}, \qquad \hat{L}(r) = \left(\frac{3\,\hat{K}(r)}{4\pi}\right)^{1/3} \;\text{(3D)}$$

The plotted quantity is $\hat{L}(r) - r$:
- $> 0$: clustering at scale $r$ (more pairs than CSR predicts)
- $\approx 0$: CSR
- $< 0$: regularity / spatial inhibition

**Implementation notes**:
- $r_\text{max}$ defaults to $\min(\Delta x, \Delta y)/4$ (one quarter of the shorter data span), keeping the estimator in the statistically reliable range
- If $n > 2000$, a random subsample of 2000 points is used and a warning is shown
- The analysis operates on the current 2D or 3D point set (localizations or TID medians, depending on the "Show medians" toggle)

No Monte Carlo CSR envelope is provided. Deviations from zero should be interpreted in the context of data set size and the observation window geometry.

---

#### 5. Cross-Channel Nearest Neighbor & Co-localization Curve

For 2-color data, computes the **directional nearest-neighbor distance** from each point in one group to its nearest point in the other group. Distances are computed in **2D (x, y only)** using `RANN::nn2`:

$$d_i^{G_1 \to G_2} = \min_{j \in G_2} \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}, \quad i \in G_1$$

Both directions ($G_1 \to G_2$ and $G_2 \to G_1$) are computed independently.

**Co-localization curve** (empirical CDF of NN distances):

$$F(r) = \frac{1}{|G_1|}\; \#\!\left\{i \in G_1 : d_i^{G_1 \to G_2} \leq r\right\}$$

This shows the fraction of group-1 representative points that have a group-2 neighbor within distance $r$. A steep rise at small $r$ indicates strong co-localization; a gradual rise suggests coincidental proximity. The same Unit selector as Batch NN controls whether per-TID medians or per-cluster centroids are used.

---

#### 6. Bivariate Ripley's L₁₂(r)

Tests whether two point populations co-localize beyond chance given their individual densities.

**Cross K function** (naive estimator):

$$\hat{K}_{12}(r) = \frac{A}{n_1 n_2} \sum_{i \in G_1}\sum_{j \in G_2} \mathbf{1}\!\left(\|\mathbf{p}_i - \mathbf{p}_j\| \leq r\right)$$

where $n_1 = |G_1|$, $n_2 = |G_2|$, and $A$ is the joint bounding area.

**L-transformation** (same as univariate):

$$\hat{L}_{12}(r) = \sqrt{\frac{\hat{K}_{12}(r)}{\pi}} \;\text{(2D)}, \qquad \hat{L}_{12}(r) = \left(\frac{3\,\hat{K}_{12}(r)}{4\pi}\right)^{1/3} \;\text{(3D)}$$

The plotted $\hat{L}_{12}(r) - r$:
- $> 0$: groups co-localize more than expected for two independent CSR processes
- $\approx 0$: spatial independence
- $< 0$: spatial exclusion

When Cluster2 is active, $\hat{L}_{12}$ is computed **per Cluster2** and the results are averaged, giving a within-cluster measure of inter-channel co-localization rather than a global one. The same 2000-point subsampling limit as univariate $L(r)$ is applied per group.

---

#### 7. Pair Cross-Correlation g₁₂(r)

The pair cross-correlation function (PCF) is the annulus-normalized derivative of $K_{12}$ and is more sensitive to the *scale* of co-localization than the cumulative $L_{12}$:

**2D**:

$$g_{12}(r) = \frac{A}{n_1 n_2} \cdot \frac{\#\text{cross-pairs in annulus}\bigl(r-\tfrac{\delta}{2},\, r+\tfrac{\delta}{2}\bigr)}{\pi\!\left[(r+\tfrac{\delta}{2})^2 - (r-\tfrac{\delta}{2})^2\right]}$$

**3D**:

$$g_{12}(r) = \frac{V}{n_1 n_2} \cdot \frac{\#\text{cross-pairs in shell}}{\tfrac{4}{3}\pi\!\left[(r+\tfrac{\delta}{2})^3 - (r-\tfrac{\delta}{2})^3\right]}$$

where $\delta$ is the bin width (the same $r$-step used for L₁₂). Under CSR with spatial independence between groups, $g_{12}(r) = 1$ for all $r$. Values $> 1$ indicate co-localization at scale $r$; values $< 1$ indicate exclusion.

The PCF pinpoints the scale of maximum co-localization, whereas $L_{12}(r)$ accumulates all pairs up to $r$ and can mask scale-specific effects.

---

#### 8. 2D KDE Density Map & Co-localization

Kernel density estimation (KDE) maps are computed for each color group on a shared 2D grid and rendered as heatmaps. Pixel-level overlap between the two maps quantifies co-localization without requiring explicit clustering.

**Bandwidth (Silverman's rule of thumb)**:

$$h = 0.9 \cdot \min\!\left(\hat{\sigma},\; \frac{\text{IQR}}{1.34}\right) \cdot n^{-1/5}$$

where $\hat{\sigma}$ is the sample standard deviation and IQR/1.34 is a robust scale estimate ($\approx$ normal equivalent SD). The bandwidth is computed on the combined x-coordinates of both groups and applied identically in x and y, ensuring a consistent grid for both channels. KDE is computed via `MASS::kde2d`.

**Co-localization metric (Pearson r of density surfaces)**:

$$r = \frac{\sum_{i,j}(z_{1,ij} - \bar{z}_1)(z_{2,ij} - \bar{z}_2)}{\sqrt{\sum_{i,j}(z_{1,ij}-\bar{z}_1)^2 \cdot \sum_{i,j}(z_{2,ij}-\bar{z}_2)^2}}$$

where $z_{1,ij}$ and $z_{2,ij}$ are the density values at pixel $(i,j)$ for groups 1 and 2.

> **Caveat**: the Pearson $r$ between KDE surfaces is strongly bandwidth-dependent. A large bandwidth smooths both channels, artificially inflating correlation even for spatially distinct populations. Use this metric for relative comparisons across conditions rather than as an absolute co-localization measure.

### Export

- **Export MINFLUX data (.csv)** - the currently filtered localizations
- **Download Filter Settings** - saves all filter parameters as a JSON file that can be reloaded
- **NN distance tables** - CSV export of nearest-neighbor analysis results
- **Cluster metrics table** - per-cluster morphology (position, size, localization count)

---

## Tracking Analysis Tab

### Loading & Chunking

On upload, every file is automatically split into chunks (`.qs` or `.fst` binary files, 100 TIDs per chunk by default) and indexed by a `manifest.csv`. This keeps RAM usage bounded regardless of dataset size.

Supported input formats: `.json`, `.csv` (pyMINFLUX), `.qs`, `.fst`, `.rds`.

A hard pre-filter of **CFR ≤ 0.5** is applied during chunking to remove multi-emitter events before any analysis.

Time is stored in **seconds** on disk; the app converts to **milliseconds** internally for all displays and calculations.

### Filter Pipeline

Filters are applied sequentially in three stages. Each stage operates on the output of the previous one:

```
Raw data
   └── [1] Metadata filters  (EFO, CFR, DCR, TIM, ECO, dwell, fluo)
         └── [2] TID length filter  (number of points per trace)
               └── [3] Trace duration filter  (max(tim) − min(tim) in ms)
```

This sequential ordering matters: the TID length and duration distributions shown in the filter histograms reflect data that has already passed the upstream filters.

**Filter & Rechunk**: applies the current filter state to the *original* uploaded data (frozen at upload time, never overwritten by repeated rechunking) and writes a new set of chunk files. This ensures repeated Filter & Rechunk cycles are consistent where resetting all filters and rechunking always restores the full dataset.

### Playback

Traces can be animated in two modes:
- **Sequential**: plays through each TID one at a time, showing a configurable trailing window of past positions
- **Simultaneous**: plays all TIDs in real time, synchronized by their absolute timestamps

### MSD Analysis

Mean Squared Displacement analysis is available via the **MSD** panel. Click **Plot MSD** to compute.

See [MSD Computation](#msd-computation) and [Diffusion Coefficient Fitting](#diffusion-coefficient-fitting) for the mathematical details.

Key controls:

| Parameter | Description |
|-----------|-------------|
| **τ max (ms)** | Maximum lag time to compute |
| **τ bin (ms)** | Width of each lag-time bin |
| **Min points per track** | Tracks with fewer localizations are excluded |
| **Fit % (first N%)** | Fraction of the MSD curve used for linear fitting (early time points are most reliable) |
| **Anomalous diffusion** | If checked, fits MSD ∝ τ^α instead of the linear Brownian model |
| **Weight by TID** | If checked, ensemble MSD is the unweighted mean across tracks; if unchecked, weights by number of pairs per bin |
| **Min track fraction (%)** | Masks ensemble MSD bins where fewer than this fraction of tracks contribute |

**Per-track D and α histograms**: after clicking Plot MSD, individual diffusion coefficients (D) and anomalous exponents (α) are fitted per track. Their distributions are shown as histograms.

**MSD clustering**: k-means clustering of tracks in the (log₁₀ D, log₁₀ CI) feature space identifies mobility sub-populations. See [Confinement Index](#confinement-index).

**dt histogram**: distribution of consecutive inter-localization time intervals within each track; useful for assessing temporal sampling regularity.

### Export

- **Filtered data (.csv or .qs)** — the rechunked/filtered dataset
- **MSD table (.csv)** — tau, ensemble MSD, and pair counts

---

## MINFLUX Simulator Tab

Generates synthetic MINFLUX datasets for testing analysis pipelines or illustrating spatial organization patterns.

**Geometric shapes** (2D and 3D):

| Shape | Description |
|-------|-------------|
| Cluster | Round cluster with uniform or Gaussian radial distribution |
| Ring | Circular arrangement with controllable width |
| Ellipse | Filled or hollow ellipse/ellipsoid |
| Filament | Linear arrangement with optional lateral jitter and helical twist |
| Double ring | Two parallel rings separated along Z (e.g. nuclear pore complex) |
| Arc | Partial ring spanning a user-defined angular range |
| Hollow sphere | Shell with adjustable thickness |

**Image-guided mode**: loads a confocal image and places synthetic molecules proportionally to local pixel intensity, mimicking the underlying fluorophore distribution.

Localizations are placed according to the chosen geometry, then Gaussian noise at a user-specified precision (nm) is added to simulate MINFLUX localization uncertainty.

Simulated data can be sent directly to the **SMLM Analysis** or **Tracking Analysis** tabs for immediate visualization and analysis.

---

## Key Calculations

### Unit Conversions & JSON Import

Abberior JSON files store coordinates in **meters**. On load, the app converts to nanometers:

$$x_\text{nm} = x_\text{m} \times 10^9$$

<<<<<<< HEAD
**Z-axis scaling** is a user-set calibration, not a fixed constant. Files are read with **raw z**;
the factor you set (SMLM sidebar *Axial (z) calibration*, Tracking sidebar *Axial (z) scaling factor*,
both defaulting to 0.7) is applied once afterwards, on Plot Data. CSVs are first returned to raw z — the
app detects what factor they already carry, exactly when the file was exported from here (a
`z_scale_applied` column) and by assumption otherwise. Exports write that column so a round-trip is
exact. Because the factor multiplies every axial quantity — z spread, cluster extent, 3D neighbour
distances, Rg — applying it twice would compress z by 30%, which the detection prevents.

The historical note below describes the mechanics of the correction itself; scaling is applied after
unit conversion to correct for refractive index mismatch between the sample medium and the objective immersion medium:
=======
**Z-axis scaling** is applied after unit conversion to correct for refractive index mismatch between the sample medium and the objective immersion medium:
>>>>>>> c4383f3519c4268414a76a4243a2691a25c5543b

$$z_\text{corrected} = z_\text{nm} \times s_z \quad (s_z = 0.7 \text{ by default})$$

**Multi-iteration MINFLUX**: each TID produces localizations at several iteration steps (`itr`). By default only the final iteration (highest `itr` value, typically the most precise) is retained. The CFR value used is the last non-zero CFR recorded across all iterations for that TID (*last_nonzero* carry-forward strategy), because CFR is often only populated in early iterations.

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

The SD filter slider ranges set bounds $[\sigma_\text{min},\, \sigma_\text{max}]$; a TID is kept only if its $\sigma$ falls within those bounds on every active axis. Setting the upper bound to a small value retains only spatially tight traces (useful for removing freely diffusing molecules), while a broad range passes all traces.

The **SD panel** histogram shows the distribution of per-TID spreads and overlays the current filter range. The median SD across all traces is reported numerically.

---

### DCR-Based 2-Color Separation

**Per-TID median DCR** is used as a robust summary to classify each trace into one of two spectral groups:

$$\widetilde{\text{dcr}}_\text{tid} = \mathrm{median}(\{\text{dcr}_i : \text{tid}_i = \text{tid}\})$$

**Manual mode**: the user defines two ranges $[d_{1,\text{lo}},\, d_{1,\text{hi}})$ and $[d_{2,\text{lo}},\, d_{2,\text{hi}}]$. A TID is assigned to group 1 or 2 if its median DCR falls in the respective range; TIDs outside both ranges are unassigned.

**Auto-detect mode (GMM)**: fits a two-component Gaussian mixture model to the distribution of $\{\widetilde{\text{dcr}}_\text{tid}\}$ using the `mclust` package:

$$p(\text{dcr}) = \pi_1 \,\mathcal{N}(\mu_1, \sigma_1^2) + \pi_2 \,\mathcal{N}(\mu_2, \sigma_2^2)$$

Components are sorted so $\mu_1 < \mu_2$ (group 1 = lower DCR = donor-channel-dominant fluorophore). The split threshold $d^*$ is the crossing point of the two weighted component densities:

$$\pi_1 \,\mathcal{N}(d^*;\mu_1, \sigma_1^2) = \pi_2 \,\mathcal{N}(d^*;\mu_2, \sigma_2^2)$$

This equation is solved numerically in $(\mu_1, \mu_2)$ using `uniroot`. If the two component means differ by less than 0.05 DCR units, the distribution is considered unimodal and auto-detection is aborted with a warning.

---

### Spatial Density

For each point (localization or TID median), the local density is defined as the count of other points within a sphere (or disk in 2D) of radius $r$:

$$\rho_i = \#\{j \neq i : \|p_i - p_j\| \leq r\}$$

Computed via `dbscan::frNN` (fixed-radius nearest-neighbor search), which uses a kd-tree for $O(n \log n)$ performance. The radius $r$ is user-controlled (default 150 nm).

---

### Ripley's L(r)

The $L(r)$ function is a variance-stabilized version of Ripley's $K(r)$. MINFLUXR uses a **naive (boundary-uncorrected) estimator**:

$$\hat{K}(r) = \frac{A \cdot 2}{n(n-1)} \sum_{i < j} \mathbf{1}(\|\mathbf{p}_i - \mathbf{p}_j\| \leq r)$$

where $A$ is the convex-hull area (or volume in 3D), $n$ is the number of points, and the sum counts all unique pairs within distance $r$. The L-transformation:

$$\hat{L}(r) = \sqrt{\frac{\hat{K}(r)}{\pi}} \quad \text{(2D)}, \qquad \hat{L}(r) = \left(\frac{3\,\hat{K}(r)}{4\pi}\right)^{1/3} \quad \text{(3D)}$$

The plotted quantity is $\hat{L}(r) - r$:

- $> 0$: clustering at scale $r$ (more pairs than expected under CSR)
- $\approx 0$: complete spatial randomness (CSR)
- $< 0$: regularity / inhibition

No edge correction is applied; this introduces upward bias for points near the data boundary. $r_\text{max}$ is automatically capped at $\min(\Delta x, \Delta y)/4$ to keep the estimator in the reliable central region. If $n > 2000$, a random subsample of 2000 points is used. See the [Ripley's L(r)](#4-ripleys-lr) analysis section for full details.

---

### Radius of Gyration

The radius of gyration quantifies the spatial compactness of a Cluster2. In MINFLUXR, **RoG is computed from per-TID median positions** rather than from individual localizations, so each blinking event contributes one representative point and localization noise is suppressed.

For a cluster of $m$ TIDs with median positions $(\tilde{x}_k, \tilde{y}_k, \tilde{z}_k)$, the centroid is $(\bar{x}, \bar{y}, \bar{z})$ (mean of TID medians) and:

$$R_g = \sqrt{\frac{1}{m}\sum_{k=1}^{m}\left[(\tilde{x}_k - \bar{x})^2 + (\tilde{y}_k - \bar{y})^2 + (\tilde{z}_k - \bar{z})^2\right]}$$

Only clusters with $m \geq 2$ TIDs are included. See the [Radius of Gyration](#3-radius-of-gyration-rog) analysis section for full details.

---

### MSD Computation

MINFLUX tracking data is irregularly sampled in time (the inter-localization interval varies within and between traces). Standard equal-lag MSD analysis is therefore not applicable. Instead, the app uses **time-binned MSD**:

For a single track $\{(x_i, y_i, z_i, t_i)\}$ sorted by time, all pairs $(i, j)$ with $i < j$ are considered. Each pair contributes its squared displacement to the bin whose center $\tau$ satisfies $\tau_\text{lo} \leq t_j - t_i < \tau_\text{hi}$:

$$\text{MSD}(\tau) = \frac{1}{N_\tau} \sum_{\substack{(i,j)\,:\,\\ \tau_\text{lo} \leq t_j - t_i < \tau_\text{hi}}} \left[(x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2\right]$$

where $N_\tau$ is the number of pairs in that bin.

The **ensemble MSD** across all tracks is computed bin-by-bin as either:
- The **pair-weighted mean**: $\overline{\text{MSD}}(\tau) = \sum_k N_{\tau,k} \cdot \text{MSD}_k(\tau) \,/\, \sum_k N_{\tau,k}$ (downweights bins with few observations)
- The **track-unweighted mean**: simple average of per-track MSD values, treating each track equally regardless of how many pairs it contributes

Bins where fewer than a user-specified fraction of tracks contribute are masked as NA to avoid noise from statistically unreliable late time points.

---

### Diffusion Coefficient Fitting

**Brownian diffusion** model:

$$\text{MSD}(\tau) = 2n D \tau + 2n \sigma_\text{loc}^2$$

where $n$ is the number of dimensions (2 or 3), $D$ is the diffusion coefficient (nm²/ms = µm²/s × 10⁻³), and $\sigma_\text{loc}$ is the localization precision (nm). The intercept $2n\sigma_\text{loc}^2$ arises from localization noise and does not reflect true motion.

The fit uses the first $k$% of lag-time bins (default 30%) where the linear approximation holds best. Fitting is two-pass weighted least squares:
1. Unweighted linear regression to obtain initial estimates $D_0$, $\sigma_0$
2. Weighted linear regression with optimal weights $w_k \propto N_{\tau,k} / \widehat{\text{MSD}}(\tau_k)^2$, which minimizes variance for Poisson-like pair-count noise

**Anomalous diffusion** model (optional):

$$\text{MSD}(\tau) = 2n D_\alpha \tau^\alpha + 2n \sigma_\text{loc}^2$$

Fitted by non-linear least squares (Levenberg–Marquardt via `minpack.lm`). The anomalous exponent $\alpha$ characterizes the diffusion mode:
- $\alpha < 1$: subdiffusion (e.g. confined or obstructed)
- $\alpha = 1$: normal Brownian diffusion
- $\alpha > 1$: superdiffusion (e.g. directed transport)

Reported units: $D$ in µm²/s, $\sigma_\text{loc}$ in nm.

---

### Confinement Index

The confinement index (CI) provides a model-free measure of diffusion character:

$$\text{CI} = \frac{\text{MSD}(\tau_\text{late}) / \text{MSD}(\tau_\text{early})}{\tau_\text{late} / \tau_\text{early}}$$

where $\tau_\text{early}$ and $\tau_\text{late}$ are the medians of the first and last 15% of lag-time bins (minimum 2 bins each). For free Brownian diffusion MSD ∝ τ, so CI = 1. Deviations:

- **CI ≪ 1**: MSD plateaus -> confined diffusion
- **CI ≈ 1**: linear growth -> Brownian
- **CI ≫ 1**: superlinear growth -> directed motion

Per-track fits produce a table of ($D$, $\alpha$, $\sigma_\text{loc}$, CI). **MSD clustering** applies k-means in the (log₁₀ $D$, log₁₀ CI) feature space (z-scored, with log₁₀ $D$ upweighted ×2 as the primary discriminator) to identify mobility sub-populations. Clusters are ordered from most confined (cluster 1) to most mobile (cluster $k$).

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

- **CFR pre-filter**: In the Tracking tab, events with CFR > 0.5 are removed automatically during chunking. This hard cutoff is standard for excluding multi-emitter contamination and is applied before all user-defined filters.
- **Z-axis scaling**: The default Z scale of 0.7 corrects for the refractive index mismatch between water (n ≈ 1.33) and glass (n ≈ 1.52) in a typical sample. Adjust this if your mounting medium differs.
- **JSON CFR carry-forward**: Because CFR is measured during early MINFLUX iterations and may be NaN in the final iteration, the app takes the last non-zero CFR value across all iterations. This ensures CFR filtering works correctly on JSON data even when only the final iteration is retained.
- **Multiple datasets**: In the SMLM tab, up to several files can be loaded simultaneously and compared using the Active file selector. Filter settings can be locked when switching to prevent accidental resets.
- **Filter & Rechunk is non-destructive**: the original chunk set is frozen at upload time. Every Filter & Rechunk operation reads from that frozen original, so you can always recover the full dataset by resetting all filters and rechunking again.
- **MSD fit range**: fitting only the first 20–30% of the MSD curve is strongly recommended. Later time points are based on fewer pairs, have higher variance, and deviate from the simple linear model due to confinement or insufficient sampling.
