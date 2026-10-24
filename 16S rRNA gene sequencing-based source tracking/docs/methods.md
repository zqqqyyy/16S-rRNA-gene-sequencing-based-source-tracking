# SeawaterTracker numerical methods

This source release retains the numerical behavior recovered from the supplied software. It does not establish equivalence to another SourceTracker implementation or reproduce all manuscript results.

## Source-mixture estimation

`source_tracking.py` is the custom sampler recovered from `app.core.sourcetracker`. Its behavior has not been established as equivalent to the official SourceTracker algorithm.

| Public parameter | Original parameter | Mathematical/operational role | Archived class default | Supplied data configuration |
| --- | --- | --- | --- | --- |
| `burn_in_iterations` | `burnin` | Iterations before retained draws | 1000 | 100 |
| `posterior_draws` | `n_draws` | Number of retained draws per run | 100 | 100 |
| `thinning_interval` | `delay` | Iteration spacing between retained draws | 10 | 10 |
| `source_feature_prior` | `alpha1` | Additive smoothing of source feature counts | 0.001 | 0.001 |
| `source_mixture_prior` | `alpha2` | Mixture-assignment and Dirichlet prior | 0.001 | 0.1 |
| `unknown_feature_prior` | `beta` | Unknown-source feature prior | 10 | 0.01 |
| `resampling_depth` | `rarefaction_depth` | Multinomial sink resampling depth; zero disables it | 0 | 1000 |
| `random_seed` | `seed` | Base MT19937 seed | 42 | 42 |
| `independent_runs` | `n_independent_runs` | Independent run means used for summary statistics | 5 | 5 |

`configs/source_tracking.json` records the supplied data-folder configuration. The class defaults remain explicit in Python. The manuscript's “restart 10” cannot be mapped to `thinning_interval`: the latter is thinning, and there is no separate restart parameter in the recovered class.

For known source `k` and feature `j`, source probabilities are normalized from `source_counts[k, j] + source_feature_prior`. Sink feature counts are initially assigned among known sources and an extra unknown component by a uniform multinomial draw. At each update, assignments for the current feature are removed and reassigned with probabilities based on remaining source totals and feature probabilities. Each iteration also draws a Dirichlet mixture from source totals plus `source_mixture_prior`.

Two implementation details are retained explicitly:

1. The operation formerly called rarefaction uses multinomial sampling **with replacement**. Sinks below the requested depth are omitted. It is not conventional sampling without replacement.
2. The current feature's unknown-source count is read after its assignment column has been reset to zero. Consequently that count is zero in the unknown-feature probability numerator at the update. This rule has not been silently corrected.

Run `i`, numbered from zero, uses `random_seed + 1000 * i`; the same sequence of run seeds is reused for each sink. A local `numpy.random.RandomState` preserves the original MT19937 draws without resetting NumPy's global RNG.

Per-run proportions are the mean retained draws. Final proportions are the mean of those run means. Final standard deviations are population standard deviations (`ddof=0`) across run means, not the pooled posterior-draw uncertainty. Relative standard deviation is `100 * std / mean`, with zero assigned for means at or below `1e-10`.

The archived output tables contain RP contributions of approximately 72.5570% (JST), 92.8748% (SLR), and 78.4914% (YY). These do not match the manuscript range 73.12–93.24%. The exported inputs and recovered configuration therefore should not be represented as a demonstrated reconstruction of all paper figures. This release verifies implementation preservation on controlled source-tracking inputs, not the full scientific result chain.

## Primer trimming and DADA2

Cutadapt adapter orientations, filenames, read-length thresholds, and discard settings preserve the archived wrapper. The reverse-complement mapping in the software interchanges `S` and `W`; standard IUPAC `S` and `W` are each self-complementary. The compatibility function is therefore named `archived_reverse_complement`. The default reverse primer contains `W`, so this issue can affect the derived adapter sequence. Its downstream scientific effect has not been quantified. The standard mapping can be inspected in [Biopython's IUPAC data definitions](https://github.com/biopython/biopython/blob/master/Bio/Data/IUPACData.py).

The public wrapper invokes `python -m cutadapt` in the active Python environment. The original GUI-oriented fallback that modified process-global CLI/logging state was removed. Command arguments were compared against the executable with a mocked backend; raw reads were not processed for this release.

The DADA2 script performs filtering, error learning, denoising, paired-read merging, sequence-table construction, consensus chimera removal, MD5-based ASV identifiers, and tracking-table output. Its local variables and formatting have been normalized. R 4.4.3 parsed both scripts, and the syntax trees were identical after applying the recorded variable mapping. DADA2 calls, argument names, literals, mathematical expressions, and output schemas remain unchanged. That check does not validate raw-read pairing, biological parameters, or full DADA2 execution.

