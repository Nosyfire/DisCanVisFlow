"""
tests/test_create_blast_table_worker.py

Tests for bin/create_blast_table_worker.py

Uses the 5-record XML fixtures extracted from real BLAST runs:
    tests/dummy_data/dummy_uniprotdb_gencode_query.xml
    tests/dummy_data/dummy_gencodedb_uniprot_query.xml

Run from the project root:
    pytest tests/test_create_blast_table_worker.py -v
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).parent.parent
BIN          = PROJECT_ROOT / "bin" / "create_blast_table_worker.py"
DUMMY_DIR    = PROJECT_ROOT / "tests" / "dummy_data"
BLAST1_XML   = DUMMY_DIR / "dummy_uniprotdb_gencode_query.xml"  # query=GENCODE, db=UniProt
BLAST2_XML   = DUMMY_DIR / "dummy_gencodedb_uniprot_query.xml"  # query=UniProt, db=GENCODE

EXPECTED_BEST_COLS = [
    "Gencode",
    "Uniprot",
    "alignmentpuntcuality_x",
    "coverage_x",
    "alignmentpuntcuality_y",
    "coverage_y",
]
EXPECTED_ISO_COLS = [
    "Gencode", "Uniprot", "alignmentpuntcuality", "coverage",
    "identity", "identity/aln_len", "aln_len",
    "seq_len_x", "seq_len_y",
    "region_len_x", "region_len_y",
    "Starx", "Endx", "Stary", "Endy",
]


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def run_worker(tmp_path: Path, coverage: float = 70.0) -> dict[str, pd.DataFrame]:
    """Invoke the worker and return a dict of output DataFrames."""
    cmd = [
        sys.executable, str(BIN),
        "--blast1_xml", str(BLAST1_XML),
        "--blast2_xml", str(BLAST2_XML),
        "--output_dir", str(tmp_path),
        "--coverage",   str(coverage),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, (
        f"Worker failed (code {r.returncode})\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    )
    return {
        "best": pd.read_csv(tmp_path / "bestsequences.tsv",    sep="\t"),
        "all":  pd.read_csv(tmp_path / "allsequences.tsv",     sep="\t"),
        "iso":  pd.read_csv(tmp_path / "isoformssequences.tsv",sep="\t"),
    }


# ---------------------------------------------------------------------------
# Tests: output schema
# ---------------------------------------------------------------------------

class TestOutputSchema:

    def test_bestsequences_columns(self, tmp_path):
        dfs = run_worker(tmp_path)
        assert list(dfs["best"].columns) == EXPECTED_BEST_COLS

    def test_allsequences_columns(self, tmp_path):
        dfs = run_worker(tmp_path)
        assert list(dfs["all"].columns) == EXPECTED_BEST_COLS

    def test_isoforms_columns(self, tmp_path):
        dfs = run_worker(tmp_path)
        assert list(dfs["iso"].columns) == EXPECTED_ISO_COLS

    def test_three_files_created(self, tmp_path):
        run_worker(tmp_path)
        for fname in ("bestsequences.tsv", "allsequences.tsv", "isoformssequences.tsv"):
            assert (tmp_path / fname).exists(), f"Missing: {fname}"


# ---------------------------------------------------------------------------
# Tests: column content correctness
# ---------------------------------------------------------------------------

class TestColumnContent:

    def test_gencode_column_contains_ensp(self, tmp_path):
        dfs = run_worker(tmp_path)
        assert dfs["best"]["Gencode"].str.startswith("ENSP").all(), (
            "Gencode column should contain ENSP headers"
        )

    def test_uniprot_column_starts_with_sp_or_tr(self, tmp_path):
        dfs = run_worker(tmp_path)
        starts_ok = dfs["best"]["Uniprot"].str.startswith(("sp|", "tr|"))
        assert starts_ok.all(), (
            "Uniprot column should start with 'sp|' or 'tr|'"
        )

    def test_alignment_quality_values(self, tmp_path):
        dfs = run_worker(tmp_path)
        valid = {"identical", "aligned"}
        for col in ("alignmentpuntcuality_x", "alignmentpuntcuality_y"):
            bad = set(dfs["best"][col].unique()) - valid
            assert not bad, f"Unexpected quality values in {col}: {bad}"

    def test_coverage_is_float_0_to_100(self, tmp_path):
        dfs = run_worker(tmp_path)
        for col in ("coverage_x", "coverage_y"):
            assert dfs["best"][col].between(0, 100).all(), (
                f"Coverage values out of [0,100] in {col}"
            )


# ---------------------------------------------------------------------------
# Tests: reciprocal filtering
# ---------------------------------------------------------------------------

class TestReciprocalFilter:

    def test_best_subset_of_all(self, tmp_path):
        """Every (Gencode, Uniprot) pair in bestsequences is also in allsequences."""
        dfs = run_worker(tmp_path, coverage=0)
        best_pairs = set(zip(dfs["best"]["Gencode"], dfs["best"]["Uniprot"]))
        all_pairs  = set(zip(dfs["all"]["Gencode"],  dfs["all"]["Uniprot"]))
        assert best_pairs.issubset(all_pairs), "bestsequences contains rows not in allsequences"

    def test_coverage_filter_reduces_rows(self, tmp_path):
        """bestsequences (high cov) ≤ allsequences (no cov filter) in row count."""
        dfs = run_worker(tmp_path, coverage=70)
        assert len(dfs["best"]) <= len(dfs["all"])

    def test_all_best_rows_meet_coverage(self, tmp_path):
        """All rows in bestsequences have coverage ≥ the threshold on BOTH sides."""
        threshold = 70.0
        dfs = run_worker(tmp_path, coverage=threshold)
        mask = (dfs["best"]["coverage_x"] >= threshold) & (dfs["best"]["coverage_y"] >= threshold)
        assert mask.all(), "Some rows in bestsequences are below the coverage threshold"


# ---------------------------------------------------------------------------
# Tests: isoforms table
# ---------------------------------------------------------------------------

class TestIsoformsTable:

    def test_isoforms_has_more_rows_than_best(self, tmp_path):
        """Isoforms TSV (no filter) should have ≥ as many rows as bestsequences."""
        dfs = run_worker(tmp_path, coverage=70)
        # isoforms are from blast1 only (no merge), so they're independent
        assert len(dfs["iso"]) > 0, "isoformssequences.tsv is unexpectedly empty"

    def test_isoforms_coverage_formula(self, tmp_path):
        """iso.coverage == iso.identity / iso.seq_len_x * 100 (rounded to 3 dp)."""
        dfs = run_worker(tmp_path)
        expected = (dfs["iso"]["identity"] / dfs["iso"]["seq_len_x"] * 100).round(3)
        pd.testing.assert_series_equal(
            dfs["iso"]["coverage"].round(3),
            expected.round(3),
            check_names=False,
        )

    def test_region_lengths_consistent(self, tmp_path):
        """region_len_x == Endx - Starx + 1."""
        dfs = run_worker(tmp_path)
        expected_x = dfs["iso"]["Endx"] - dfs["iso"]["Starx"] + 1
        pd.testing.assert_series_equal(
            dfs["iso"]["region_len_x"],
            expected_x,
            check_names=False,
        )


# ---------------------------------------------------------------------------
# Tests: edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:

    def test_help_exits_cleanly(self):
        r = subprocess.run(
            [sys.executable, str(BIN), "--help"],
            capture_output=True, text=True
        )
        assert r.returncode == 0

    def test_missing_file_fails(self, tmp_path):
        r = subprocess.run(
            [sys.executable, str(BIN),
             "--blast1_xml", str(tmp_path / "nonexistent.xml"),
             "--blast2_xml", str(BLAST2_XML),
             "--output_dir", str(tmp_path)],
            capture_output=True, text=True
        )
        assert r.returncode != 0
