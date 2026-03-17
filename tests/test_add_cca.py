"""Tests for add_cca.py — must produce identical output to legacy addCCA.pl."""

import subprocess
import tempfile
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
BIN = Path(__file__).parent.parent / "pipeline" / "bin"


def test_add_cca_matches_expected():
    """Run add_cca.py on test input and compare to expected output."""
    input_fa = FIXTURES / "add_cca_input.fa"
    expected_fa = FIXTURES / "add_cca_expected.fa"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name

    subprocess.run(
        ["python3", str(BIN / "add_cca.py"), str(input_fa), tmp_path],
        check=True,
    )

    actual = Path(tmp_path).read_text()
    expected = expected_fa.read_text()
    assert actual == expected, (
        f"Output mismatch.\n\nExpected:\n{expected}\n\nActual:\n{actual}"
    )
    Path(tmp_path).unlink()


def test_pseudogene_filtered():
    """Entries with 'pseudo' in header (case-insensitive) must be excluded."""
    input_fa = FIXTURES / "add_cca_input.fa"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name

    subprocess.run(
        ["python3", str(BIN / "add_cca.py"), str(input_fa), tmp_path],
        check=True,
    )

    output = Path(tmp_path).read_text()
    assert "pseudo" not in output.lower()
    Path(tmp_path).unlink()


def test_cca_appended():
    """Every sequence line must end with 'cca'."""
    input_fa = FIXTURES / "add_cca_input.fa"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name

    subprocess.run(
        ["python3", str(BIN / "add_cca.py"), str(input_fa), tmp_path],
        check=True,
    )

    for line in Path(tmp_path).read_text().strip().split("\n"):
        if not line.startswith(">"):
            assert line.endswith("cca"), f"Sequence does not end with 'cca': {line}"
    Path(tmp_path).unlink()


def test_sequences_lowercased():
    """All sequence characters must be lowercase."""
    input_fa = FIXTURES / "add_cca_input.fa"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name

    subprocess.run(
        ["python3", str(BIN / "add_cca.py"), str(input_fa), tmp_path],
        check=True,
    )

    for line in Path(tmp_path).read_text().strip().split("\n"):
        if not line.startswith(">"):
            assert line == line.lower(), f"Sequence not lowercase: {line}"
    Path(tmp_path).unlink()


def test_entry_count():
    """Should output 3 entries (4 input minus 1 pseudogene)."""
    input_fa = FIXTURES / "add_cca_input.fa"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name

    subprocess.run(
        ["python3", str(BIN / "add_cca.py"), str(input_fa), tmp_path],
        check=True,
    )

    headers = [l for l in Path(tmp_path).read_text().strip().split("\n") if l.startswith(">")]
    assert len(headers) == 3, f"Expected 3 entries, got {len(headers)}"
    Path(tmp_path).unlink()
