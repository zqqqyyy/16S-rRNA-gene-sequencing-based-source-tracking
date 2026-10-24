# SeawaterTracker research core

Standalone, readable research code recovered from **SeawaterTracker_portable.rar**. The numerical modules use descriptive variables, explicit configuration, and structured results. This folder contains the sequencing source-tracking workflow and runs independently of the combined antibiotic/source-tracking release.

[中文使用说明](README_zh-CN.md) · [Numerical methods](docs/methods.md) · [Data formats](docs/data_format.md)

## Files and responsibilities

| File | Task |
| --- | --- |
| `source_tracking.py` | Source-mixture sampling, independent runs, mean/std/RSD summaries |
| `primer_trimming.py` | Paired-read primer removal through Cutadapt |
| `data_preparation.py` | ASV count-table loading, feature alignment, environment aggregation |
| `dada2_denoise.R` | DADA2 filtering, denoising, merging, chimera removal and ASV export |
| `run_source_tracking.py` | Standalone entry point for analysis of count tables |
| `configs/source_tracking.json` | Parameters from the supplied data-folder configuration |

The first two Python modules and the R script are byte-identical to their validated versions in the combined research release. Count-table functions were retained without numerical changes. The standalone runner is a new interface around the recovered sampler.

## Quick start

Use Python 3.11. From this folder, create a virtual environment and install the recorded numerical dependencies. Windows PowerShell:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_source_tracking.py --sources examples/source_counts.tsv --sinks examples/sink_counts.tsv --config examples/demo_parameters.json --output results/demo
```

On Linux/macOS, create and activate a Python 3.11 virtual environment, then use `python -m pip install -r requirements.txt` and the same command beginning with `python run_source_tracking.py`.

The example uses synthetic counts and a short configuration. It writes `mean_proportions.tsv`, `std_proportions.tsv`, `relative_std_percent.tsv`, and `run_metadata.json`. The last file records parameters and input-file hashes. Means and standard deviations are fractions; relative standard deviations are percentages.

Analyze the archived study inputs with:

```bash
python run_source_tracking.py --sources data/sequencing/source_environment_counts.tsv.gz --sinks data/sequencing/aligned_sink_counts.tsv.gz --config configs/source_tracking.json --output results/study
```

These inputs contain 21 source environments, 3 sinks and 78,705 aligned features. This analysis takes longer than the synthetic example and was not rerun in full during release preparation. See `data/sequencing/input_derivation.json` for their origin.

## Calling the functions

```python
from data_preparation import read_count_table
from source_tracking import SourceTrackingParameters, estimate_source_contributions

sources = read_count_table("examples/source_counts.tsv")
sinks = read_count_table("examples/sink_counts.tsv")
parameters = SourceTrackingParameters(
    burn_in_iterations=20,
    posterior_draws=20,
    thinning_interval=2,
    source_mixture_prior=0.1,
    unknown_feature_prior=0.01,
    independent_runs=3,
)
result = estimate_source_contributions(sinks, sources, parameters)
print(result.mean_proportions)
```

For importable installation, use `python -m pip install . --no-deps` after installing the requirements. It also exposes a `seawater-source` command with the same arguments as `run_source_tracking.py`.

## Read preprocessing

Install `python -m pip install ".[preprocessing]"` to enable the Cutadapt wrapper. `primer_trimming.trim_paired_reads` receives forward/reverse read paths, an output directory, and optional `PrimerTrimmingParameters`. Its archived ambiguous-base complement behavior is documented in [Methods](docs/methods.md#primer-trimming-and-dada2).

The R script's recorded environment is R 4.4.3, DADA2 1.34.0 and digest 0.6.37:

```bash
Rscript dada2_denoise.R trimmed_reads results/dada2 2 2 0 0 12 4
```

Arguments after the directories are forward/reverse maximum expected errors, forward/reverse truncation lengths, minimum overlap and thread count. Full FASTQ processing was not run for this release.

## Validation and provenance

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
```

The frozen sampler fixtures were produced from the original Python 3.11 executable code objects. Direct original-code comparisons and standalone execution results are recorded in `validation/`. R syntax trees were verified to match after the documented local-variable renaming. The source archive SHA-256 and module hashes are recorded in `validation/archive_provenance.json`.

This is a recovery and reorganization of the supplied software's implementation. The original resampling rule, unknown-source update and primer mapping are retained and explained in Methods. These checks do not demonstrate reproduction of all manuscript results. Authorship, licence and paper citation metadata should reflect the authors' publication decisions.
