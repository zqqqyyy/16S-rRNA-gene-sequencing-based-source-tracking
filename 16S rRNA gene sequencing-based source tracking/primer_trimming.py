"""Paired-read primer trimming with explicit Cutadapt arguments and results."""

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PrimerTrimmingParameters:
    """Primer sequences and Cutadapt options from the archived data configuration."""

    forward_primer: str = "ACTCCTACGGGAGGCAGCA"
    reverse_primer: str = "GGACTACHVGGGTWTCTAAT"
    discard_untrimmed: bool = True
    minimum_read_length: int = 100
    threads: int = 4

    def __post_init__(self) -> None:
        if not self.forward_primer or not self.reverse_primer:
            raise ValueError("Both primer sequences must be nonempty.")
        for name in ("minimum_read_length", "threads"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer.")


@dataclass(frozen=True)
class PrimerTrimmingResult:
    """Output read paths, parsed counts and the complete Cutadapt log."""

    forward_reads: Path
    reverse_reads: Path
    processed_pairs: int | None
    retained_pairs: int | None
    log_text: str


def archived_reverse_complement(sequence: str) -> str:
    """Reproduce the archived primer mapping, including its S/W interchange.

    This explicitly named compatibility function is not a general IUPAC utility:
    the archived implementation maps S to W and W to S, whereas both standard
    symbols are self-complementary. Retaining it avoids changing preprocessing
    during a source-code refactor; its scientific impact remains unassessed.
    """
    complement_by_base = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "R": "Y",
        "Y": "R",
        "M": "K",
        "K": "M",
        "S": "W",
        "W": "S",
        "H": "D",
        "D": "H",
        "B": "V",
        "V": "B",
        "N": "N",
    }
    return "".join(complement_by_base.get(base, base) for base in reversed(sequence.upper()))


def build_cutadapt_arguments(
    forward_reads: str | Path,
    reverse_reads: str | Path,
    output_directory: str | Path,
    parameters: PrimerTrimmingParameters | None = None,
) -> tuple[list[str], Path, Path]:
    """Build the archived adapter options without launching a subprocess."""
    parameters = parameters or PrimerTrimmingParameters()
    sample_stem = Path(forward_reads).stem
    sample_stem = (
        sample_stem.replace(".R1", "").replace("_R1", "").replace(".fq", "").replace(".fastq", "")
    )
    output_directory = Path(output_directory)
    forward_output = output_directory / f"{sample_stem}_trimmed_R1.fastq.gz"
    reverse_output = output_directory / f"{sample_stem}_trimmed_R2.fastq.gz"
    arguments = [
        "-g",
        parameters.forward_primer,
        "-G",
        parameters.reverse_primer,
        "-a",
        archived_reverse_complement(parameters.reverse_primer),
        "-A",
        archived_reverse_complement(parameters.forward_primer),
        "-j",
        str(parameters.threads),
        "-m",
        str(parameters.minimum_read_length),
        "-o",
        str(forward_output),
        "-p",
        str(reverse_output),
    ]
    if parameters.discard_untrimmed:
        arguments.append("--discard-untrimmed")
    arguments.extend([str(forward_reads), str(reverse_reads)])
    return arguments, forward_output, reverse_output


def parse_cutadapt_statistics(log_text: str) -> dict[str, int]:
    """Read processed/retained pair counts from a Cutadapt text report."""
    statistics = {}
    for line in log_text.splitlines():
        if "Total read pairs processed" in line:
            statistics["processed_pairs"] = int(line.split(":")[-1].strip().replace(",", ""))
        elif "Pairs written (passing filters)" in line:
            count_text = line.split(":")[-1].strip().split("(")[0].strip().replace(",", "")
            statistics["retained_pairs"] = int(count_text)
    return statistics


def trim_paired_reads(
    forward_reads: str | Path,
    reverse_reads: str | Path,
    output_directory: str | Path,
    parameters: PrimerTrimmingParameters | None = None,
) -> PrimerTrimmingResult:
    """Run Cutadapt in this Python environment using archived adapter options.

    Install the ``preprocessing`` extra to make ``python -m cutadapt`` available.
    Process execution is separate from numerical analysis and does not change
    global argv, stdin or logging handlers.
    """
    for reads_path in (forward_reads, reverse_reads):
        if not Path(reads_path).is_file():
            raise FileNotFoundError(reads_path)
    arguments, forward_output, reverse_output = build_cutadapt_arguments(
        forward_reads, reverse_reads, output_directory, parameters
    )
    Path(output_directory).mkdir(parents=True, exist_ok=True)
    process = subprocess.run(
        [sys.executable, "-m", "cutadapt", *arguments],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdin=subprocess.DEVNULL,
        check=False,
    )
    log_text = process.stderr + process.stdout
    if process.returncode != 0:
        raise RuntimeError(f"Cutadapt failed with exit code {process.returncode}:\n{log_text}")
    statistics = parse_cutadapt_statistics(log_text)
    return PrimerTrimmingResult(
        forward_output,
        reverse_output,
        statistics.get("processed_pairs"),
        statistics.get("retained_pairs"),
        log_text,
    )
