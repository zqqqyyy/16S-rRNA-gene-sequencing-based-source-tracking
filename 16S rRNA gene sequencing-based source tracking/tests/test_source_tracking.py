"""Numerical regression, reproducibility and input-contract checks."""

import numpy as np
import pandas as pd
import pytest

from source_tracking import (
    SourceTrackingParameters,
    estimate_source_contributions,
    resample_sink_counts,
    sample_source_proportions,
)


@pytest.fixture
def parameters():
    return SourceTrackingParameters(
        burn_in_iterations=4,
        posterior_draws=6,
        thinning_interval=2,
        source_feature_prior=0.001,
        source_mixture_prior=0.1,
        unknown_feature_prior=0.01,
        random_seed=42,
        independent_runs=3,
    )


def test_all_posterior_draws_match_archived_code(archived_results, parameters):
    result = sample_source_proportions(
        archived_results["sink_counts"], archived_results["source_counts"], parameters
    )
    expected = archived_results["single_run"]
    np.testing.assert_array_equal(result.posterior_samples, expected["draws"])
    np.testing.assert_array_equal(result.mean_proportions, expected["mean"])
    np.testing.assert_array_equal(result.std_proportions, expected["std"])


def test_independent_runs_and_rsd_match_archived_code(archived_results, parameters):
    expected = archived_results["multiple_sinks"]
    sources = pd.DataFrame(archived_results["source_counts"], index=["source_A", "source_B"])
    sinks = pd.DataFrame(expected["counts"], index=expected["sample_ids"])
    result = estimate_source_contributions(sinks, sources, parameters)
    np.testing.assert_array_equal(result.mean_proportions.to_numpy(), expected["mean"])
    np.testing.assert_array_equal(result.std_proportions.to_numpy(), expected["std"])
    np.testing.assert_array_equal(result.relative_std_percent.to_numpy(), expected["rsd"])
    assert list(result.mean_proportions.columns) == ["source_A", "source_B", "Unknown"]
    assert result.mean_proportions.index.name == "sample_id"
    for sample_id, expected_runs in expected["runs"].items():
        np.testing.assert_array_equal(result.run_proportions[sample_id], expected_runs)


def test_resampling_preserves_draws_and_drops_low_depth_samples(archived_results):
    data = archived_results["multiple_sinks"]
    sinks = pd.DataFrame(data["counts"], index=data["sample_ids"])
    result = resample_sink_counts(sinks, 20, random_seed=42)
    expected = archived_results["resampled_depth_20"]
    assert list(result.index) == expected["sample_ids"]
    np.testing.assert_array_equal(result.values, expected["counts"])
    assert (result.sum(axis=1) == 20).all()


def test_sampler_does_not_change_process_wide_random_state(parameters):
    before = np.random.get_state()
    sample_source_proportions([2, 4], [[20, 1], [1, 20]], parameters)
    after = np.random.get_state()
    assert before[0] == after[0]
    np.testing.assert_array_equal(before[1], after[1])
    assert before[2:] == after[2:]


def test_zero_sink_has_zero_summary_without_division_warning(parameters):
    result = estimate_source_contributions(
        pd.DataFrame([[0, 0]], index=["empty"]),
        pd.DataFrame([[2, 3]], index=["source_A"]),
        parameters,
    )
    assert not result.relative_std_percent.isna().any().any()
    assert result.mean_proportions.to_numpy().sum() == 0
    assert result.relative_std_percent.to_numpy().sum() == 0


@pytest.mark.parametrize(
    "changes",
    [
        {"posterior_draws": 0},
        {"independent_runs": 0},
        {"thinning_interval": 0},
        {"burn_in_iterations": -1},
        {"source_feature_prior": 0},
        {"unknown_feature_prior": float("nan")},
        {"random_seed": 2**32 - 1},
        {"posterior_draws": True},
    ],
)
def test_invalid_parameters_are_rejected(changes):
    with pytest.raises(ValueError):
        SourceTrackingParameters(**changes)


@pytest.mark.parametrize("sink", [[1, -1], [1, 0.5], [1, float("nan")], [1, 2, 3]])
def test_invalid_count_vectors_are_rejected(sink, parameters):
    with pytest.raises(ValueError):
        sample_source_proportions(sink, [[2, 3]], parameters)


def test_feature_order_is_checked(parameters):
    with pytest.raises(ValueError, match="match exactly"):
        estimate_source_contributions(
            pd.DataFrame([[1, 2]], columns=["b", "a"]),
            pd.DataFrame([[3, 4]], columns=["a", "b"]),
            parameters,
        )
