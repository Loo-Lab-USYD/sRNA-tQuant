"""Tests for filter_genome_reads.py — must produce identical output to legacy removeGenomeMapper.pl.

Test fixture logic:
  - read1: maps to pre-tRNA contig, length 35 → KEEP
  - read2: maps to chr6 (genome), length 35 → DISCARD
  - read3: maps to pre-tRNA contig, length 35 → KEEP
  - read4: two alignments — one to pre-tRNA, one to chr1 (genome) → DISCARD both
  - read5: maps to pre-tRNA contig, length 20 (< 30, exempt from check) → KEEP
"""

import subprocess
import tempfile
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_filter(pretrna_fa, sam_in):
    with tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as f:
        out_path = f.name

    subprocess.run(
        ["python3", str(BIN / "filter_genome_reads.py"),
         str(pretrna_fa), str(sam_in), out_path],
        check=True,
    )
    output = Path(out_path).read_text()
    Path(out_path).unlink()
    return output


def test_output_matches_expected():
    output = run_filter(
        FIXTURES / "filter_genome_pretrnas.fa",
        FIXTURES / "filter_genome_input.sam",
    )
    expected = (FIXTURES / "filter_genome_expected.sam").read_text()
    assert output == expected, (
        f"Output mismatch.\n\nExpected:\n{expected}\n\nActual:\n{output}"
    )


def test_genome_mapped_reads_removed():
    """Reads mapping to genome contigs must not appear in output."""
    output = run_filter(
        FIXTURES / "filter_genome_pretrnas.fa",
        FIXTURES / "filter_genome_input.sam",
    )
    assert "read2" not in output
    assert "read4" not in output


def test_pretrna_reads_kept():
    """Reads mapping only to pre-tRNA contigs must be kept."""
    output = run_filter(
        FIXTURES / "filter_genome_pretrnas.fa",
        FIXTURES / "filter_genome_input.sam",
    )
    assert "read1" in output
    assert "read3" in output


def test_short_read_kept():
    """Reads shorter than 30 nt are exempt from genome-mapping check."""
    output = run_filter(
        FIXTURES / "filter_genome_pretrnas.fa",
        FIXTURES / "filter_genome_input.sam",
    )
    assert "read5" in output


def test_multimapper_with_genome_hit_fully_removed():
    """If any alignment of a read hits the genome, ALL alignments are removed."""
    output = run_filter(
        FIXTURES / "filter_genome_pretrnas.fa",
        FIXTURES / "filter_genome_input.sam",
    )
    read4_lines = [l for l in output.strip().split("\n") if l.startswith("read4")]
    assert len(read4_lines) == 0
