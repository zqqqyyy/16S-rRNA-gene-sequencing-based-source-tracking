"""The study's source-mixture sampler with explicit parameters and results.

This module preserves the recovered numerical update rules. It is a custom
implementation, not a claim of equivalence to the official SourceTracker package.
See docs/methods.md for resampling and unknown-source limitations.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class SourceTrackingParameters:
    """Sampling parameters; priors retain their original mathematical roles.

    ``source_feature_prior`` corresponds to alpha1, ``source_mixture_prior``
    to alpha2 and ``unknown_feature_prior`` to beta. ``thinning_interval``
    is an iteration interval, not a restart count. The defaults below match
    the archived class; configs/source_tracking.json records the data-folder
    configuration supplied with the application.
    """

    burn_in_iterations: int = 1000
    posterior_draws: int = 100
    thinning_interval: int = 10
    source_feature_prior: float = 0.001
    source_mixture_prior: float = 0.001
    unknown_feature_prior: float = 10.0
    resampling_depth: int = 0
    random_seed: int = 42
    independent_runs: int = 5

    def __post_init__(self) -> None:
        for name in ("burn_in_iterations", "resampling_depth", "random_seed"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer.")
        for name in ("posterior_draws", "thinning_interval", "independent_runs"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
                raise ValueError(f"{name} must be a positive integer.")
        for name in ("source_feature_prior", "source_mixture_prior", "unknown_feature_prior"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.random_seed + 1000 * (self.independent_runs - 1) > 2**32 - 1:
            raise ValueError("Seeds for all independent runs must fit into the MT19937 seed range.")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> "SourceTrackingParameters":
        """Construct from a JSON object, rejecting misspelled/unknown keys."""
        return cls(**dict(values))


@dataclass(frozen=True)
class SamplingResult:
    """One run's mean, population standard deviation and retained draws."""

    mean_proportions: FloatArray
    std_proportions: FloatArray
    posterior_samples: FloatArray


@dataclass(frozen=True)
class SourceTrackingResult:
    """Per-sample summaries across independent run means.

    ``relative_std_percent`` is 100 * population_std / mean. Zero/near-zero
    means receive zero, as in the recovered code. ``run_proportions`` stores
    each sample's array of shape (independent_runs, number_of_sources).
    """

    mean_proportions: pd.DataFrame
    std_proportions: pd.DataFrame
    relative_std_percent: pd.DataFrame
    run_proportions: dict[str, FloatArray]


def _validate_counts(values: ArrayLike, *, name: str, dimensions: int) -> NDArray:
    counts = np.asarray(values)
    if counts.ndim != dimensions:
        raise ValueError(f"{name} must have {dimensions} dimension(s).")
    if not np.issubdtype(counts.dtype, np.number):
        raise ValueError(f"{name} must contain numeric counts.")
    if not np.isfinite(counts).all() or (counts < 0).any() or (counts != np.floor(counts)).any():
        raise ValueError(f"{name} must contain finite nonnegative integer-valued counts.")
    if (counts > np.iinfo(np.int64).max).any():
        raise ValueError(f"{name} contains counts outside the supported integer range.")
    return counts


def resample_sink_counts(
    sink_table: pd.DataFrame, depth: int, *, random_seed: int = 42
) -> pd.DataFrame:
    """Multinomial resampling with replacement, preserving the study behavior.

    Samples below ``depth`` reads are dropped. This is not conventional
    without-replacement rarefaction, despite the name used in the old software.
    """
    if isinstance(depth, bool) or not isinstance(depth, Integral) or depth < 1:
        raise ValueError("depth must be a positive integer.")
    if sink_table.index.has_duplicates:
        raise ValueError("Sink sample identifiers must be unique.")
    _validate_counts(sink_table.to_numpy(), name="sink_table", dimensions=2)
    random_generator = np.random.RandomState(random_seed)
    resampled_table = pd.DataFrame(0, index=sink_table.index, columns=sink_table.columns)
    retained_samples = []
    for sample_id in sink_table.index:
        sample_counts = sink_table.loc[sample_id].to_numpy().astype(int)
        total_reads = sample_counts.sum()
        if total_reads < depth:
            continue
        feature_probabilities = sample_counts / total_reads
        resampled_table.loc[sample_id] = random_generator.multinomial(depth, feature_probabilities)
        retained_samples.append(sample_id)
    return resampled_table.loc[retained_samples]


def sample_source_proportions(
    sink_counts: ArrayLike,
    source_counts: ArrayLike,
    parameters: SourceTrackingParameters,
    *,
    random_seed: int | None = None,
) -> SamplingResult:
    """Sample a single sink using the recovered feature-block update rule.

    ``source_counts`` has shape (known sources, features), and ``sink_counts``
    has shape (features,). The final output component is the unknown source.
    A private MT19937 generator gives the same seeded draws without changing
    NumPy's process-wide random state.
    """
    source_array = _validate_counts(source_counts, name="source_counts", dimensions=2)
    sink_array = _validate_counts(sink_counts, name="sink_counts", dimensions=1).astype(int)
    known_source_count, feature_count = source_array.shape
    if known_source_count < 1 or feature_count < 1 or sink_array.size != feature_count:
        raise ValueError("Provide at least one source and aligned, nonempty feature vectors.")
    total_source_count = known_source_count + 1
    seed = parameters.random_seed if random_seed is None else random_seed
    random_generator = np.random.RandomState(seed)
    if sink_array.sum() == 0:
        return SamplingResult(
            np.zeros(total_source_count),
            np.zeros(total_source_count),
            np.empty((0, total_source_count)),
        )
    source_feature_probabilities = source_array.astype(float) + parameters.source_feature_prior
    source_feature_probabilities /= source_feature_probabilities.sum(axis=1, keepdims=True)
    assigned_feature_counts = np.zeros((total_source_count, feature_count), dtype=int)
    for feature_index in range(feature_count):
        if sink_array[feature_index] > 0:
            initial_probabilities = np.ones(total_source_count) / total_source_count
            assigned_feature_counts[:, feature_index] = random_generator.multinomial(
                sink_array[feature_index], initial_probabilities
            )
    source_read_totals = assigned_feature_counts.sum(axis=1)
    total_iterations = parameters.burn_in_iterations + (
        parameters.posterior_draws * parameters.thinning_interval
    )
    retained_draws = []
    for iteration in range(total_iterations):
        for feature_index in range(feature_count):
            if sink_array[feature_index] == 0:
                continue
            previous_assignments = assigned_feature_counts[:, feature_index].copy()
            assigned_feature_counts[:, feature_index] = 0
            source_read_totals -= previous_assignments
            assignment_probabilities = np.zeros(total_source_count)
            for source_index in range(known_source_count):
                assignment_probabilities[source_index] = (
                    source_read_totals[source_index] + parameters.source_mixture_prior
                ) * source_feature_probabilities[source_index, feature_index]
            unknown_read_total = source_read_totals[known_source_count]
            # This value is zero after the column reset, as in the archived rule.
            unknown_feature_count = assigned_feature_counts[known_source_count, feature_index]
            unknown_feature_probability = (
                unknown_feature_count + parameters.unknown_feature_prior
            ) / (unknown_read_total + feature_count * parameters.unknown_feature_prior)
            assignment_probabilities[known_source_count] = (
                source_read_totals[known_source_count] + parameters.source_mixture_prior
            ) * unknown_feature_probability
            total_probability = assignment_probabilities.sum()
            if total_probability > 0:
                assignment_probabilities /= total_probability
            else:
                assignment_probabilities = np.ones(total_source_count) / total_source_count
            new_assignments = random_generator.multinomial(
                sink_array[feature_index], assignment_probabilities
            ).astype(int)
            assigned_feature_counts[:, feature_index] = new_assignments
            source_read_totals += new_assignments
        source_proportions = random_generator.dirichlet(
            source_read_totals + parameters.source_mixture_prior
        )
        if (
            iteration >= parameters.burn_in_iterations
            and (iteration - parameters.burn_in_iterations) % parameters.thinning_interval == 0
        ):
            retained_draws.append(source_proportions.copy())
    posterior_samples = np.array(retained_draws)
    return SamplingResult(
        posterior_samples.mean(axis=0),
        posterior_samples.std(axis=0, ddof=0),
        posterior_samples,
    )


def estimate_source_contributions(
    sink_table: pd.DataFrame,
    source_table: pd.DataFrame,
    parameters: SourceTrackingParameters | None = None,
) -> SourceTrackingResult:
    """Estimate source contributions for aligned feature-count tables.

    Rows are sink samples/source environments, columns are feature identifiers.
    Row order is retained, and feature columns must match exactly. Each sample
    uses seeds ``random_seed + 1000 * run_index``, preserving the archived
    sequence of independent runs. Output tables append an ``Unknown`` column.
    """
    parameters = parameters or SourceTrackingParameters()
    if not sink_table.columns.equals(source_table.columns):
        raise ValueError("Source and sink feature columns must match exactly and in order.")
    if sink_table.index.has_duplicates or source_table.index.has_duplicates:
        raise ValueError("Sample and source identifiers must be unique.")
    if source_table.columns.has_duplicates:
        raise ValueError("Feature identifiers must be unique.")
    if "Unknown" in source_table.index:
        raise ValueError("'Unknown' is reserved for the additional inferred source.")
    _validate_counts(source_table.to_numpy(), name="source_table", dimensions=2)
    _validate_counts(sink_table.to_numpy(), name="sink_table", dimensions=2)
    if source_table.shape[0] < 1 or source_table.shape[1] < 1:
        raise ValueError("source_table must contain at least one source and one feature.")
    if parameters.resampling_depth > 0:
        sink_table = resample_sink_counts(
            sink_table, parameters.resampling_depth, random_seed=parameters.random_seed
        )
    source_labels = list(source_table.index) + ["Unknown"]
    mean_by_sample = {}
    std_by_sample = {}
    relative_std_by_sample = {}
    run_proportions = {}
    for sample_id in sink_table.index:
        independent_means = []
        for run_index in range(parameters.independent_runs):
            run_result = sample_source_proportions(
                sink_table.loc[sample_id].to_numpy(),
                source_table.to_numpy(),
                parameters,
                random_seed=parameters.random_seed + run_index * 1000,
            )
            independent_means.append(run_result.mean_proportions)
        independent_means = np.array(independent_means)
        mean_proportions = independent_means.mean(axis=0)
        std_proportions = independent_means.std(axis=0, ddof=0)
        relative_std = np.zeros_like(mean_proportions)
        np.divide(
            std_proportions, mean_proportions, out=relative_std, where=mean_proportions > 1e-10
        )
        mean_by_sample[sample_id] = mean_proportions
        std_by_sample[sample_id] = std_proportions
        relative_std_by_sample[sample_id] = relative_std * 100
        run_proportions[sample_id] = independent_means

    def result_table(values: dict) -> pd.DataFrame:
        return pd.DataFrame(values, index=source_labels).T.rename_axis("sample_id")

    return SourceTrackingResult(
        result_table(mean_by_sample),
        result_table(std_by_sample),
        result_table(relative_std_by_sample),
        run_proportions,
    )
