"""Tests for filter_mature_reads.py — must produce identical output to legacy removePrecursor.pl.

Test fixture logic (flank=50, flankEnd=44, tRNA block length=171, mature end=127):
  - readA: start=55 > 50, end=89 <= 127 → KEEP (mature region)
  - readB: start=10 <= 50 → SKIP (in 5' flanking region)
  - readC: start=60 > 50, end=94 <= 127 → KEEP (mature region, different tRNA)
  - readD: start=55, CIGAR=33M2I, end=87 <= 127 → KEEP (insertion adjusts end)
  - readE: start=55 > 50, end=89 <= 127 → KEEP (mature region)
"""

import subprocess
import tempfile
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_filter(bed12, sam, fastq):
    result = subprocess.run(
        ["python3", str(BIN / "filter_mature_reads.py"),
         str(bed12), str(sam), str(fastq)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def test_output_matches_expected():
    output = run_filter(
        FIXTURES / "filter_mature_pretrnas.bed12",
        FIXTURES / "filter_mature_input.sam",
        FIXTURES / "filter_mature_input.fastq",
    )
    expected = (FIXTURES / "filter_mature_expected.fastq").read_text()
    assert output == expected, (
        f"Output mismatch.\n\nExpected:\n{expected}\n\nActual:\n{output}"
    )


def test_flanking_read_excluded():
    """readB maps to the 5' flanking region (start=10 <= 50) and must be excluded."""
    output = run_filter(
        FIXTURES / "filter_mature_pretrnas.bed12",
        FIXTURES / "filter_mature_input.sam",
        FIXTURES / "filter_mature_input.fastq",
    )
    assert "@readB" not in output


def test_mature_reads_kept():
    """Reads mapping within the mature tRNA body must be present."""
    output = run_filter(
        FIXTURES / "filter_mature_pretrnas.bed12",
        FIXTURES / "filter_mature_input.sam",
        FIXTURES / "filter_mature_input.fastq",
    )
    assert "@readA" in output
    assert "@readC" in output
    assert "@readE" in output


def test_insertion_cigar_handled():
    """readD has a 2I insertion — CIGAR must be parsed to compute correct end position."""
    output = run_filter(
        FIXTURES / "filter_mature_pretrnas.bed12",
        FIXTURES / "filter_mature_input.sam",
        FIXTURES / "filter_mature_input.fastq",
    )
    assert "@readD" in output


def test_full_fasta_header_ref_names():
    """SAM ref_names with full FASTA headers (name::chrom:start-end(strand))
    must be correctly matched to BED12 entries.
    Regression test: real data uses headers like
      tRNA-AlaAGC-1-1::chr6:28868691-28868862(+)
    but BED12 keys are tRNA-AlaAGC-1-1(+)."""
    output = run_filter(
        FIXTURES / "filter_mature_pretrnas.bed12",
        FIXTURES / "filter_mature_fullheader_input.sam",
        FIXTURES / "filter_mature_input.fastq",
    )
    # readA: start=55 > 50, end=89 <= 127 → KEEP
    assert "@readA" in output
    # readB: start=10 <= 50 → SKIP
    assert "@readB" not in output
    # readC: start=60 > 50, end=94 <= 127 → KEEP
    assert "@readC" in output


def test_output_is_valid_fastq():
    """Output must have lines in groups of 4: @header, seq, +, qual."""
    output = run_filter(
        FIXTURES / "filter_mature_pretrnas.bed12",
        FIXTURES / "filter_mature_input.sam",
        FIXTURES / "filter_mature_input.fastq",
    )
    lines = output.strip().split("\n")
    assert len(lines) % 4 == 0, f"Expected multiple of 4 lines, got {len(lines)}"
    for i in range(0, len(lines), 4):
        assert lines[i].startswith("@"), f"Line {i} should be a FASTQ header: {lines[i]}"
        assert lines[i + 2] == "+", f"Line {i+2} should be '+': {lines[i+2]}"
