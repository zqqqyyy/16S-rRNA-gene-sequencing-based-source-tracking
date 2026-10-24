# Sequencing input and output formats

Count tables are tab-separated text, optionally gzip-compressed. The first column contains row identifiers and the remaining columns contain nonnegative integer-valued ASV counts. Sources have one row per already-aggregated environment; sinks have one row per query sample.

Source and sink feature columns must have identical names and order. Identifiers must be unique. `align_feature_tables` aligns the feature union, fills absent counts with zero and removes features with zero total counts in both inputs. `aggregate_source_environments` accepts sample counts plus metadata indexed by sample identifier with `role` and `environment` columns. It sums the available samples whose role is `source` by environment.

`configs/source_tracking.json` records the supplied study data-folder settings. `examples/demo_parameters.json` uses fewer iterations and synthetic data for a quick run. Detailed parameter names and their original counterparts are listed in Methods.

Output rows use `sample_id`; output columns are source labels followed by `Unknown`. `Unknown` is reserved for the inferred extra component. Sinks below the selected resampling depth are omitted; depth zero disables resampling. A zero-count sink with resampling disabled produces zero summaries.

`mean_proportions.tsv` and `std_proportions.tsv` use fractions, while `relative_std_percent.tsv` uses percentages. Standard deviations are population standard deviations across independent run means. The JSON run record stores parameters, basenames and SHA-256 hashes; rerunning into the same output directory replaces the corresponding output files.

DADA2 files retain their original schemas, including `SampleID`, `ASV_ID` and `Sequence`. The count-table reader uses the first column as the row identifier regardless of its header. Source-table derivation and data hashes are recorded under `data/`.
