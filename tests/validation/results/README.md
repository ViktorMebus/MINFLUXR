# Validation results

One CSV per validation run, each produced by the script of the same name in `tests/validation/`.
These are the measurements the method documentation cites; the scripts regenerate them from scratch.

| File | Script | What it measures |
|---|---|---|
| `hmm_recovery_sweep.csv`, `_summary.csv` | `hmm_recovery_sweep.R` | Where the diffusive-state HMM recovers D and dwell times, over dwell time and D contrast |
| `spatial_envelope_coverage.csv` | `spatial_envelope_coverage.R` | How often random data are called significant by the Ripley L, g(r), L₁₂ and g₁₂ envelopes |
| `hotspot_csr_null.csv` | `hotspot_csr_null.R` | The baseline every hotspot statistic takes on random patterns, across density and cell shape |
| `hotspot_realdata_recheck.csv` | `hotspot_realdata_recheck.R` | The same statistics on the labelled two-condition acquisitions |
| `fiducial_labelled.csv` | `fiducial_labelled.R` | Fiducial detection on labelled acquisitions with and without beads |
| `hmm_step_budget.csv` | `hmm_step_budget.R` | What the HMM step budget costs: D and dwell recovery at each budget, against known values |

The `dataset` column is anonymised (`file_01`, `file_02`, …): the acquisitions are local and not part
of the repository. The scripts that read real data take their input from an environment variable, so
they run against whatever you point them at:

| Variable | Used by | Default |
|---|---|---|
| `MINFLUXR_TRACK_PATTERN` | the tracking scripts | `track.*_filtered.*csv$`, first match under `Test data/` |
| `MINFLUXR_SMLM_TESTDIR` | the SMLM and hotspot scripts | `Test data/`, searched recursively |

Each script prints the file it resolved to, because the pattern can match more than one acquisition and
the numbers only mean something against a known input.

Run logs (`*.log`) are not tracked: they are the console output of the run that produced the CSV next to
it, and they are recreated by re-running the script. The scripts that need `Test data/` say so in their
header; that folder is not part of the repository.
