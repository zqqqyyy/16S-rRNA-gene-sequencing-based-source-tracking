"""Command-construction checks do not require Cutadapt or sequencing reads."""

from primer_trimming import (
    PrimerTrimmingParameters,
    archived_reverse_complement,
    build_cutadapt_arguments,
    parse_cutadapt_statistics,
)


def test_archived_complement_is_explicitly_preserved():
    assert archived_reverse_complement("ACGT") == "ACGT"
    assert archived_reverse_complement("S") == "W"
    assert archived_reverse_complement("W") == "S"


def test_cutadapt_arguments_preserve_adapter_orientation(tmp_path):
    parameters = PrimerTrimmingParameters()
    arguments, forward_path, reverse_path = build_cutadapt_arguments(
        "sample_R1.fastq.gz", "sample_R2.fastq.gz", tmp_path, parameters
    )
    assert arguments[arguments.index("-g") + 1] == parameters.forward_primer
    assert arguments[arguments.index("-a") + 1] == archived_reverse_complement(
        parameters.reverse_primer
    )
    assert "--discard-untrimmed" in arguments
    assert forward_path.name == "sample_trimmed_R1.fastq.gz"
    assert reverse_path.name == "sample_trimmed_R2.fastq.gz"


def test_pair_statistics_are_machine_readable():
    log_text = "Total read pairs processed: 1,200\nPairs written (passing filters): 1,000 (83.3%)"
    assert parse_cutadapt_statistics(log_text) == {"processed_pairs": 1200, "retained_pairs": 1000}
