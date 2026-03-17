"""Tests for salmon_to_rpm.py — convert Salmon quant.sf to anticodon-level RPM.

Test fixture: 3 clusters, 1000 total reads (same as compute_expression tests).
  - cluster1 (AlaAGC): 500 reads → 500000 RPM
  - cluster2 (ValCAC): 300 reads → 300000 RPM
  - cluster3 (GlyGCC): 200 reads → 200000 RPM
"""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_salmon_to_rpm(quant_sf, cluster_info):
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        out_path = f.name

    subprocess.run(
        ["python3", str(BIN / "salmon_to_rpm.py"),
         str(quant_sf), str(cluster_info), out_path],
        check=True,
    )
    df = pd.read_csv(out_path, index_col=0)
    Path(out_path).unlink()
    return df


def test_rpm_values_correct():
    """RPM = NumReads * 1e6 / total_reads."""
    df = run_salmon_to_rpm(
        FIXTURES / "salmon_quant.sf",
        FIXTURES / "expression_clusterinfo.fa",
    )
    assert df.loc["AlaAGC", "rpm"] == 500000.0
    assert df.loc["ValCAC", "rpm"] == 300000.0
    assert df.loc["GlyGCC", "rpm"] == 200000.0


def test_rpm_equals_rpm_nomod():
    """salmon_to_rpm always sets rpm == rpm_nomod (no mod adjustment)."""
    df = run_salmon_to_rpm(
        FIXTURES / "salmon_quant.sf",
        FIXTURES / "expression_clusterinfo.fa",
    )
    for trna in ["AlaAGC", "ValCAC", "GlyGCC"]:
        assert df.loc[trna, "rpm"] == df.loc[trna, "rpm_nomod"]


def test_reads_row_present():
    """A READS row should contain the total mapped read count."""
    df = run_salmon_to_rpm(
        FIXTURES / "salmon_quant.sf",
        FIXTURES / "expression_clusterinfo.fa",
    )
    assert "READS" in df.index
    assert df.loc["READS", "rpm"] == 1000.0


def test_rpms_sum_to_million():
    """RPM values (excluding READS row) should sum to 1,000,000."""
    df = run_salmon_to_rpm(
        FIXTURES / "salmon_quant.sf",
        FIXTURES / "expression_clusterinfo.fa",
    )
    trna_rpm = df.drop("READS")
    assert abs(trna_rpm["rpm"].sum() - 1_000_000) < 0.01


def test_zero_reads_handled():
    """When no reads are mapped, output should have zero RPM and zero READS."""
    df = run_salmon_to_rpm(
        FIXTURES / "salmon_quant_empty.sf",
        FIXTURES / "expression_clusterinfo.fa",
    )
    assert "READS" in df.index
    assert df.loc["READS", "rpm"] == 0
    trna_rpm = df.drop("READS")
    assert trna_rpm["rpm"].sum() == 0.0
