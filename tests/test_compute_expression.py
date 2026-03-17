"""Tests for compute_expression.py — RPM calculation with optional modification adjustment.

Test fixture: 3 clusters, 1000 total reads, no ASE variants (simple case).
  - cluster1 (AlaAGC): 500 reads → 500000 RPM
  - cluster2 (ValCAC): 300 reads → 300000 RPM
  - cluster3 (GlyGCC): 200 reads → 200000 RPM
  - rpm and rpm_nomod should be equal (no modifications)
  - READS row should contain total reads (1000)
"""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_compute(expr, ase, info, trnas):
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        out_path = f.name

    subprocess.run(
        ["python3", str(BIN / "compute_expression.py"),
         str(expr), str(ase), str(info), str(trnas), out_path],
        check=True,
    )
    df = pd.read_csv(out_path, index_col=0)
    Path(out_path).unlink()
    return df


def test_rpm_values_correct():
    """RPM = reads * 1e6 / total_reads."""
    df = run_compute(
        FIXTURES / "expression_input.txt",
        FIXTURES / "expression_ase.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    assert df.loc["AlaAGC", "rpm"] == 500000.0
    assert df.loc["ValCAC", "rpm"] == 300000.0
    assert df.loc["GlyGCC", "rpm"] == 200000.0


def test_rpm_nomod_equals_rpm_without_modifications():
    """When no ASE variants exist, rpm and rpm_nomod should be equal."""
    df = run_compute(
        FIXTURES / "expression_input.txt",
        FIXTURES / "expression_ase.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    for trna in ["AlaAGC", "ValCAC", "GlyGCC"]:
        assert df.loc[trna, "rpm"] == df.loc[trna, "rpm_nomod"]


def test_reads_row_present():
    """A READS row should contain the total mapped read count."""
    df = run_compute(
        FIXTURES / "expression_input.txt",
        FIXTURES / "expression_ase.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    assert "READS" in df.index
    assert df.loc["READS", "rpm"] == 1000
    assert df.loc["READS", "rpm_nomod"] == 1000


def test_rpms_sum_to_million():
    """RPM values (excluding READS row) should sum to 1,000,000."""
    df = run_compute(
        FIXTURES / "expression_input.txt",
        FIXTURES / "expression_ase.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    trna_rpm = df.drop("READS")
    assert abs(trna_rpm["rpm"].sum() - 1_000_000) < 0.01
    assert abs(trna_rpm["rpm_nomod"].sum() - 1_000_000) < 0.01


def test_all_trnas_present():
    """All tRNAs from the list should appear in the output index."""
    df = run_compute(
        FIXTURES / "expression_input.txt",
        FIXTURES / "expression_ase.csv",
        FIXTURES / "expression_clusterinfo.fa",
        FIXTURES / "expression_trnas.txt",
    )
    for trna in ["AlaAGC", "ValCAC", "GlyGCC", "UndetNNN"]:
        assert trna in df.index
