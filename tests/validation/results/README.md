# Validation results

Each script in `tests/validation/` checks one method against a known answer. The runs below use
**simulated data** with a known ground truth, so their results are committed here; the method
documentation cites them.

| File | Script | What it measures |
|---|---|---|
| `spatial_envelope_coverage.csv` | `spatial_envelope_coverage.R` | How often random data are called significant by the Ripley L, g(r), L₁₂ and g₁₂ envelopes |
| `coloc_null_clustered.csv` | `coloc_null_clustered.R` | Which two-channel null answers the co-localization question on clustered channels |
| `hotspot_csr_null.csv` | `hotspot_csr_null.R` | The baseline each hotspot statistic takes on random patterns, across density and cell shape |
| `hmm_recovery_sweep.csv`, `_summary.csv` | `hmm_recovery_sweep.R` | Where the diffusive-state HMM recovers D and dwell times, over dwell time and D contrast |
| `hmm_step_budget.csv` | `hmm_step_budget.R` | What the HMM step budget costs: D and dwell recovery at each budget |

## What the simulations show

**Spatial envelopes** (CSR, or independent CSR channels, in an L-shaped cell; footprint window and edge
correction on; 60 replicates × 99 simulations; 40 radii; nominal 5%):

| Statistic | Global test, false-positive rate | Pointwise band left somewhere |
|---|---|---|
| L(r) | 1.7% | 45% |
| g(r) | 6.7% | 92% |
| L₁₂(r), random-labelling null | 1.7% | 58% |
| g₁₂(r), random-labelling null | 3.3% | 93% |
| L₁₂(r), CSR-in-window null | 21.7% | 83% |
| g₁₂(r), CSR-in-window null | 8.3% | 100% |

The global test holds near its nominal level with random labelling; the pointwise band is not a test
over many radii. The CSR-in-window null is anti-conservative for L₁₂ because the estimated window's
shape differs from the cell.

**Hotspots** (random patterns in a filled frame, a 2 × 8 µm rod and an L-shape, at 200 / 1000 / 5000
units, 20 replicates each; positive control 20% of units in six Gaussian foci): `cv_ratio`, VMR, the
fraction above the enhancement cut and `max_focus_ratio` hold their chance baseline across density and
shape. Raw `max_focus_n` and focus count move with density. Power is weak in sparse cells.

**Diffusive-state HMM** (two states at MINFLUX tracking timings, 200 tracks × 300 localizations, three
replicates per cell): D and dwell times are recovered for mean dwells of **10 cycles or more**. At 1–3
cycles the states mix, dwell times read 35–180% long and the likelihood intervals do not cover the
truth. Counting dwells along the decoded path reads them long at every dwell time.

**HMM step budget** (dwells 20 and 10 cycles): D stays within 3% of the truth at every budget down to
20,000 steps; the dwell estimate drifts only at the smallest budget.

## Scripts that run on real acquisitions

`fiducial_labelled.R`, `hotspot_realdata_recheck.R` and `hmm_real_tracking*.R` apply the methods to
real MINFLUX files. The acquisitions are not part of the repository, and neither are these scripts'
results. They take their input from an environment variable:

| Variable | Used by | Default |
|---|---|---|
| `MINFLUXR_TRACK_PATTERN` | the tracking scripts | `track.*_filtered.*csv$`, first match under `Test data/` |
| `MINFLUXR_SMLM_TESTDIR` | the SMLM and hotspot scripts | `Test data/`, searched recursively |

Each script prints the file it resolved to, because the pattern can match more than one acquisition.

Run logs (`*.log`) are not tracked; re-running a script recreates its log.

**Two-channel nulls on clustered channels** (Thomas clusters, L-shaped cell, footprint window and edge
correction on, 60 replicates × 99 simulations, global test at 5%). Rejection rate:

| Truth | Shift (default) L₁₂ / g₁₂ | Random labelling L₁₂ / g₁₂ | CSR in window L₁₂ / g₁₂ |
|---|---|---|---|
| Independent clustered species | 6.7% / 3.3% | 100% / 100%, all as exclusion | 97% / 100% |
| Both species in the same clusters | 100% / 100% | 1.7% / 3.3% | 100% / 100% |
| B within ~8 nm of A | 100% / 100% | 100% / 5% | 100% / 100% |

Only the shift null is near nominal on independent species while detecting both kinds of
co-localization. Random labelling is exact only for unclustered channels (table above).
