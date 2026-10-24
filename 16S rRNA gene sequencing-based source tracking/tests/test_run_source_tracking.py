"""Standalone entry-point checks for result tables and configuration errors."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from run_source_tracking import main

PROJECT = Path(__file__).resolve().parents[1]


def test_command_preserves_proportion_units_and_input_hashes(tmp_path):
    source_path = PROJECT / "examples/source_counts.tsv"
    assert (
        main(
            [
                "--sources",
                str(source_path),
                "--sinks",
                str(PROJECT / "examples/sink_counts.tsv"),
                "--config",
                str(PROJECT / "examples/demo_parameters.json"),
                "--output",
                str(tmp_path),
            ]
        )
        == 0
    )
    proportions = pd.read_csv(tmp_path / "mean_proportions.tsv", sep="\t", index_col=0)
    assert proportions.shape == (2, 3)
    assert proportions.columns[-1] == "Unknown"
    np.testing.assert_allclose(proportions.sum(axis=1), 1, rtol=0, atol=1e-12)
    metadata = json.loads((tmp_path / "run_metadata.json").read_text(encoding="utf-8"))
    assert (
        metadata["input_files"][0]["sha256"] == hashlib.sha256(source_path.read_bytes()).hexdigest()
    )


def test_misspelled_configuration_is_reported(tmp_path, capsys):
    configuration = tmp_path / "invalid.json"
    configuration.write_text('{"random_sead": 42}', encoding="utf-8")
    with pytest.raises(SystemExit) as error:
        main(
            [
                "--sources",
                str(PROJECT / "examples/source_counts.tsv"),
                "--sinks",
                str(PROJECT / "examples/sink_counts.tsv"),
                "--config",
                str(configuration),
                "--output",
                str(tmp_path / "results"),
            ]
        )
    assert error.value.code == 2
    assert "random_sead" in capsys.readouterr().err
    assert not (tmp_path / "results").exists()
