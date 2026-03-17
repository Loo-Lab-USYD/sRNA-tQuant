"""Tests for cluster_trnas.py — must produce identical output to legacy clustering.pl."""

import subprocess
import tempfile
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def run_clustering(input_fa):
    """Run cluster_trnas.py, return (cluster_output, info_output) as strings."""
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f1, \
         tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f2:
        cluster_path = f1.name
        info_path = f2.name

    subprocess.run(
        ["python3", str(BIN / "cluster_trnas.py"), str(input_fa), cluster_path, info_path],
        check=True,
    )
    cluster_out = Path(cluster_path).read_text()
    info_out = Path(info_path).read_text()
    Path(cluster_path).unlink()
    Path(info_path).unlink()
    return cluster_out, info_out


def test_cluster_fa_matches_expected():
    """cluster.fa should have one entry per unique sequence."""
    cluster_out, _ = run_clustering(FIXTURES / "clustering_input.fa")
    expected = (FIXTURES / "clustering_expected_cluster.fa").read_text()
    assert cluster_out == expected


def test_cluster_info_matches_expected():
    """clusterInfo.fa should map every original entry to its cluster."""
    _, info_out = run_clustering(FIXTURES / "clustering_input.fa")
    expected = (FIXTURES / "clustering_expected_info.fa").read_text()
    assert info_out == expected


def test_cluster_count():
    """6 input entries with 3 unique sequences → 3 clusters."""
    cluster_out, _ = run_clustering(FIXTURES / "clustering_input.fa")
    headers = [l for l in cluster_out.strip().split("\n") if l.startswith(">")]
    assert len(headers) == 3


def test_all_members_present_in_info():
    """All 6 original entries must appear in clusterInfo."""
    _, info_out = run_clustering(FIXTURES / "clustering_input.fa")
    headers = [l for l in info_out.strip().split("\n") if l.startswith(">")]
    assert len(headers) == 6


def test_identical_sequences_share_cluster():
    """Entries with the same sequence must have the same cluster ID."""
    _, info_out = run_clustering(FIXTURES / "clustering_input.fa")
    # Parse cluster assignments
    cluster_for = {}
    for line in info_out.strip().split("\n"):
        if line.startswith(">"):
            cluster_id = line.split(":")[0]  # e.g. ">cluster1"
            original = line.split(":", 1)[1]  # e.g. "tRNA-AlaAGC-1-1(...)"
            cluster_for[original] = cluster_id

    # AlaAGC-1-1 and AlaAGC-2-1 must share a cluster
    assert cluster_for["tRNA-AlaAGC-1-1(chr6:28868741-28868812)"] == \
           cluster_for["tRNA-AlaAGC-2-1(chr1:12345-12416)"]

    # All three GlyGCC entries must share a cluster
    gly_clusters = {cluster_for[k] for k in cluster_for if "GlyGCC" in k}
    assert len(gly_clusters) == 1
