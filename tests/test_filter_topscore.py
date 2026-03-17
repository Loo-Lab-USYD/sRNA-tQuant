"""Tests for filter_topscore.py — keep only equal top-scoring alignments per read.

Test fixture logic:
  - read1: 3 alignments with AS:i scores -5, -10, -5 → keep the two with -5
  - read2: 1 alignment with AS:i:0 → keep (single alignment)
  - read3: 2 alignments both AS:i:-8 → keep both (tied top score)
"""

import subprocess
import tempfile
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_filter(sam_in):
    with tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as f:
        out_path = f.name

    subprocess.run(
        ["python3", str(BIN / "filter_topscore.py"), str(sam_in), out_path],
        check=True,
    )
    output = Path(out_path).read_text()
    Path(out_path).unlink()
    return output


def test_headers_preserved():
    output = run_filter(FIXTURES / "filter_topscore_input.sam")
    assert "@HD" in output
    assert "@SQ" in output


def test_suboptimal_alignment_removed():
    """read1's cluster2 alignment (AS:-10) should be removed; cluster1 and cluster3 (AS:-5) kept."""
    output = run_filter(FIXTURES / "filter_topscore_input.sam")
    read1_lines = [l for l in output.strip().split("\n") if l.startswith("read1")]
    assert len(read1_lines) == 2
    for line in read1_lines:
        assert "AS:i:-5" in line
    # Verify the suboptimal one is gone
    assert not any("cluster2" in l for l in read1_lines)


def test_single_alignment_kept():
    """read2 has only one alignment — should be kept."""
    output = run_filter(FIXTURES / "filter_topscore_input.sam")
    read2_lines = [l for l in output.strip().split("\n") if l.startswith("read2")]
    assert len(read2_lines) == 1


def test_tied_scores_all_kept():
    """read3 has two alignments with equal scores — both should be kept."""
    output = run_filter(FIXTURES / "filter_topscore_input.sam")
    read3_lines = [l for l in output.strip().split("\n") if l.startswith("read3")]
    assert len(read3_lines) == 2


def test_total_output_lines():
    """4 headers + 2 (read1) + 1 (read2) + 2 (read3) = 9 lines."""
    output = run_filter(FIXTURES / "filter_topscore_input.sam")
    lines = [l for l in output.strip().split("\n") if l]
    assert len(lines) == 9, f"Expected 9 lines, got {len(lines)}"
