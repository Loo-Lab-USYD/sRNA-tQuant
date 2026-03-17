"""Tests for resolve_multimappers.py — anticodon-aware multimapper resolution.

Test fixture logic:
  cluster1 → AlaAGC, cluster2 → ValCAC, cluster3 → GlyGCC

  - unique1: NH:i:1 → output directly
  - multi_same1: 2 hits, both cluster1 (AlaAGC) → same anticodon → output one alignment
  - multi_diff1: 2 hits, cluster1 (AlaAGC) + cluster2 (ValCAC) → different → multimappers
  - multi_last: 2 hits, both cluster3 (GlyGCC) → same anticodon → output (tests final flush)
"""

import subprocess
import tempfile
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_resolve(sam_in, cluster_info):
    with tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as f1, \
         tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as f2:
        out_path = f1.name
        multi_path = f2.name

    subprocess.run(
        ["python3", str(BIN / "resolve_multimappers.py"),
         str(sam_in), str(cluster_info), out_path, multi_path],
        check=True,
    )
    out = Path(out_path).read_text()
    multi = Path(multi_path).read_text()
    Path(out_path).unlink()
    Path(multi_path).unlink()
    return out, multi


def test_sam_headers_preserved():
    out, _ = run_resolve(
        FIXTURES / "resolve_mm_input.sam",
        FIXTURES / "resolve_mm_clusterinfo.fa",
    )
    assert "@HD" in out
    assert "@SQ" in out


def test_unique_reads_kept():
    out, _ = run_resolve(
        FIXTURES / "resolve_mm_input.sam",
        FIXTURES / "resolve_mm_clusterinfo.fa",
    )
    assert "unique1" in out


def test_same_anticodon_multimapper_kept():
    """Multimapper where all hits share AlaAGC → one alignment in output."""
    out, multi = run_resolve(
        FIXTURES / "resolve_mm_input.sam",
        FIXTURES / "resolve_mm_clusterinfo.fa",
    )
    assert "multi_same1" in out
    # Should be exactly one line for this read in output
    out_lines = [l for l in out.strip().split("\n") if l.startswith("multi_same1")]
    assert len(out_lines) == 1
    # Should NOT appear in multimappers file
    assert "multi_same1" not in multi


def test_different_anticodon_multimapper_discarded():
    """Multimapper spanning AlaAGC + ValCAC → goes to multimappers file."""
    out, multi = run_resolve(
        FIXTURES / "resolve_mm_input.sam",
        FIXTURES / "resolve_mm_clusterinfo.fa",
    )
    assert "multi_diff1" not in out
    assert "multi_diff1" in multi
    # Both alignments should be in multimappers
    multi_lines = [l for l in multi.strip().split("\n") if l.startswith("multi_diff1")]
    assert len(multi_lines) == 2


def test_final_group_flushed():
    """The last group in the file must be processed (legacy bug fix)."""
    out, _ = run_resolve(
        FIXTURES / "resolve_mm_input.sam",
        FIXTURES / "resolve_mm_clusterinfo.fa",
    )
    assert "multi_last" in out


def test_output_read_count():
    """Output should have: 4 headers + unique1 + multi_same1 + multi_last = 7 lines."""
    out, _ = run_resolve(
        FIXTURES / "resolve_mm_input.sam",
        FIXTURES / "resolve_mm_clusterinfo.fa",
    )
    lines = [l for l in out.strip().split("\n") if l]
    assert len(lines) == 7, f"Expected 7 lines, got {len(lines)}: {lines}"
