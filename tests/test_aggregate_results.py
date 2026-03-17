"""Tests for aggregate_results.py — merge per-sample RPM into a matrix."""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_aggregate(rpm_files, ase_files=None):
    tmpdir = tempfile.mkdtemp()
    cmd = ["python3", str(BIN / "aggregate_results.py"), tmpdir]
    cmd.extend([str(f) for f in rpm_files])
    if ase_files:
        cmd.append("--ase")
        cmd.extend([str(f) for f in ase_files])

    subprocess.run(cmd, check=True)
    return tmpdir


def test_rpm_matrix_shape():
    """Output matrix should have 4 rows (3 tRNAs + READS) × 2 samples."""
    tmpdir = run_aggregate([
        FIXTURES / "sampleA.RPM.csv",
        FIXTURES / "sampleB.RPM.csv",
    ])
    df = pd.read_csv(Path(tmpdir) / "results_rpm.csv", index_col=0)
    assert df.shape == (4, 2)
    assert "sampleA" in df.columns
    assert "sampleB" in df.columns


def test_rpm_values_correct():
    tmpdir = run_aggregate([
        FIXTURES / "sampleA.RPM.csv",
        FIXTURES / "sampleB.RPM.csv",
    ])
    df = pd.read_csv(Path(tmpdir) / "results_rpm.csv", index_col=0)
    assert df.loc["AlaAGC", "sampleA"] == 500000.0
    assert df.loc["AlaAGC", "sampleB"] == 400000.0
    assert df.loc["ValCAC", "sampleB"] == 350000.0


def test_nomod_matrix_produced():
    tmpdir = run_aggregate([
        FIXTURES / "sampleA.RPM.csv",
        FIXTURES / "sampleB.RPM.csv",
    ])
    nomod = pd.read_csv(Path(tmpdir) / "results_rpm_nomod.csv", index_col=0)
    assert nomod.shape == (4, 2)


def test_reads_row_preserved():
    tmpdir = run_aggregate([
        FIXTURES / "sampleA.RPM.csv",
        FIXTURES / "sampleB.RPM.csv",
    ])
    df = pd.read_csv(Path(tmpdir) / "results_rpm.csv", index_col=0)
    assert "READS" in df.index
    assert df.loc["READS", "sampleA"] == 1000
    assert df.loc["READS", "sampleB"] == 2000


def test_single_sample():
    tmpdir = run_aggregate([FIXTURES / "sampleA.RPM.csv"])
    df = pd.read_csv(Path(tmpdir) / "results_rpm.csv", index_col=0)
    assert df.shape == (4, 1)
    assert "sampleA" in df.columns
