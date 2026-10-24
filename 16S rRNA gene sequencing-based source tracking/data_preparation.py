"""ASV count-table reading, feature alignment, and source-environment aggregation."""

from pathlib import Path

import pandas as pd


def read_count_table(path: str | Path) -> pd.DataFrame:
    """Read a TSV/TSV.GZ with row identifiers in the first column."""
    table = pd.read_csv(path, sep="\t", index_col=0)
    if table.index.has_duplicates or table.columns.has_duplicates:
        raise ValueError("Count-table row and feature identifiers must be unique.")
    return table


def align_feature_tables(
    source_table: pd.DataFrame, sink_table: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Align the union of ASVs, fill absent counts with zero and drop empty ASVs.

    The union and zero-column selection match the recovered analysis pipeline.
    Rows and source labels are retained in their existing order.
    """
    if source_table.columns.has_duplicates or sink_table.columns.has_duplicates:
        raise ValueError("Feature identifiers must be unique before alignment.")
    combined_features = source_table.columns.union(sink_table.columns)
    aligned_sources = source_table.reindex(columns=combined_features, fill_value=0)
    aligned_sinks = sink_table.reindex(columns=combined_features, fill_value=0)
    nonzero_features = aligned_sources.sum(axis=0) + aligned_sinks.sum(axis=0) > 0
    return aligned_sources.loc[:, nonzero_features], aligned_sinks.loc[:, nonzero_features]


def aggregate_source_environments(
    sample_count_table: pd.DataFrame, sample_metadata: pd.DataFrame
) -> pd.DataFrame:
    """Sum reference samples by environment in sorted environment order.

    Metadata is indexed by sample identifier and has ``role`` and ``environment``
    columns. Only rows whose role is ``source`` are used. Source samples absent
    from the feature table are ignored, matching the recovered aggregation rule.
    """
    if not {"role", "environment"}.issubset(sample_metadata.columns):
        raise ValueError("Metadata requires 'role' and 'environment' columns.")
    if sample_count_table.index.has_duplicates or sample_metadata.index.has_duplicates:
        raise ValueError("Sample identifiers must be unique.")
    source_metadata = sample_metadata.loc[sample_metadata["role"] == "source"]
    if source_metadata["environment"].isna().any():
        raise ValueError("Every source sample must have an environment label.")
    environment_counts = {}
    for environment in sorted(source_metadata["environment"].unique()):
        sample_ids = source_metadata.index[source_metadata["environment"] == environment]
        available_sample_ids = sample_ids.intersection(sample_count_table.index)
        environment_counts[environment] = sample_count_table.loc[available_sample_ids].sum(axis=0)
    return pd.DataFrame(environment_counts, index=sample_count_table.columns).T.rename_axis(
        "environment"
    )
