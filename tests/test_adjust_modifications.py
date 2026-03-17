"""Tests for adjust_modifications.py — post-quantification anticodon mod redistribution.

Fixtures use the same 3-cluster setup as compute_expression tests.
cluster1 = AlaAGC, anticodon triplet at 1-based positions 23-25 in sequence.
"""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_adjust(rpm_csv, ase_csv, cluster_info, trnas):
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        out_path = f.name

    subprocess.run(
        ["python3", str(BIN / "adjust_modifications.py"),
         str(rpm_csv), str(ase_csv), str(cluster_info), str(trnas), out_path],
        check=True,
    )
    df = pd.read_csv(out_path, index_col=0)
    Path(out_path).unlink()
    return df


def test_no_variants_unchanged():
    """When ASE file has no variants, RPM should be unchanged."""
    df = run_adjust(
        FIXTURES / "adjust_rpm_input.csv",
        FIXTURES / "expression_ase.csv",  # empty ASE
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    assert df.loc["AlaAGC", "rpm"] == 500000.0
    assert df.loc["ValCAC", "rpm"] == 300000.0
    assert df.loc["GlyGCC", "rpm"] == 200000.0


def test_variant_outside_anticodon_unchanged():
    """When variant is outside the anticodon triplet, RPM should be unchanged."""
    df = run_adjust(
        FIXTURES / "adjust_rpm_input.csv",
        FIXTURES / "adjust_ase_outside.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    assert df.loc["AlaAGC", "rpm"] == 500000.0


def test_variant_in_anticodon_redistributes():
    """When variant falls in the anticodon triplet, RPM should be redistributed."""
    df = run_adjust(
        FIXTURES / "adjust_rpm_input.csv",
        FIXTURES / "adjust_ase_with_mod.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    # cluster1 has 20% alt allele at anticodon position
    # 20% of 500000 = 100000 should be redistributed
    assert df.loc["AlaAGC", "rpm"] < 500000.0
    assert df.loc["AlaAGC", "rpm"] == 500000.0 - 100000.0


def test_rpm_nomod_unchanged_after_adjustment():
    """rpm_nomod column should never change, regardless of modifications."""
    df = run_adjust(
        FIXTURES / "adjust_rpm_input.csv",
        FIXTURES / "adjust_ase_with_mod.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    assert df.loc["AlaAGC", "rpm_nomod"] == 500000.0
    assert df.loc["ValCAC", "rpm_nomod"] == 300000.0
    assert df.loc["GlyGCC", "rpm_nomod"] == 200000.0


def test_reads_row_preserved():
    """READS row should be preserved through adjustment."""
    df = run_adjust(
        FIXTURES / "adjust_rpm_input.csv",
        FIXTURES / "adjust_ase_with_mod.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    assert "READS" in df.index
    assert df.loc["READS", "rpm"] == 1000.0
