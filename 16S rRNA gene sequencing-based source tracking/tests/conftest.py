"""Shared fixtures generated independently from the archived implementation."""

import json
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def archived_results():
    return json.loads(
        (Path(__file__).parent / "fixtures" / "archived_expected_results.json").read_text(
            encoding="utf-8"
        )
    )
