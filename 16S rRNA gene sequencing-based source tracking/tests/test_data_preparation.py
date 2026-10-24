"""Retained count-table alignment and source-aggregation regression check."""

import numpy as np
import pandas as pd

from data_preparation import aggregate_source_environments, align_feature_tables


def test_feature_alignment_and_source_aggregation():
    counts = pd.DataFrame(
        [[3, 0], [5, 0], [100, 0]], index=["A", "B", "C"], columns=["f1", "empty"]
    )
    metadata = pd.DataFrame(
        {"role": ["source", "source", "sink"], "environment": ["river", "river", "sea"]},
        index=counts.index,
    )
    sources = aggregate_source_environments(counts, metadata)
    assert sources.loc["river", "f1"] == 8
    sinks = pd.DataFrame([[2, 7]], index=["query"], columns=["f1", "new"])
    aligned_sources, aligned_sinks = align_feature_tables(sources, sinks)
    assert list(aligned_sources.columns) == ["f1", "new"]
    np.testing.assert_array_equal(aligned_sources.values, [[8, 0]])
    np.testing.assert_array_equal(aligned_sinks.values, [[2, 7]])
