"""Run source-mixture estimation from aligned ASV count tables."""

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from data_preparation import read_count_table
from source_tracking import SourceTrackingParameters, estimate_source_contributions

RELEASE_VERSION = "0.1.0"


def main(argv: list[str] | None = None) -> int:
    """Load explicit inputs, estimate source fractions, and save reproducible outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=RELEASE_VERSION)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--sinks", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args(argv)
    try:
        parameters = SourceTrackingParameters.from_mapping(
            json.loads(arguments.config.read_text(encoding="utf-8"))
        )
        result = estimate_source_contributions(
            read_count_table(arguments.sinks), read_count_table(arguments.sources), parameters
        )
        arguments.output.mkdir(parents=True, exist_ok=True)
        result.mean_proportions.to_csv(arguments.output / "mean_proportions.tsv", sep="\t")
        result.std_proportions.to_csv(arguments.output / "std_proportions.tsv", sep="\t")
        result.relative_std_percent.to_csv(arguments.output / "relative_std_percent.tsv", sep="\t")
        metadata = {
            "release": "SeawaterTracker research core",
            "version": RELEASE_VERSION,
            "parameters": asdict(parameters),
            "input_files": [
                {"name": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                for path in [arguments.sources, arguments.sinks, arguments.config]
            ],
        }
        (arguments.output / "run_metadata.json").write_text(
            json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
        )
    except (ValueError, OSError, TypeError, KeyError) as error:
        parser.error(str(error))
    print(f"Results saved to {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
