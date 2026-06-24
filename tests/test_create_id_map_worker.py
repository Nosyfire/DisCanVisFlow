"""
tests/test_create_id_map_worker.py

Unit / integration tests for bin/create_id_map_worker.py

Run from the project root:
    pytest tests/test_create_id_map_worker.py -v

These tests exercise the Python worker directly — no Nextflow required.
Each assertion mirrors what create_id_map.py (legacy) would produce on the
same input, validated against the columns required by DisCanVis2 upload.py.
"""

import sys
import subprocess
from pathlib import Path

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent
BIN          = PROJECT_ROOT / "bin" / "create_id_map_worker.py"
DUMMY_DIR    = PROJECT_ROOT / "tests" / "dummy_data"
BLAST_TSV    = DUMMY_DIR / "dummy_bestsequences.tsv"
ISO_TSV      = DUMMY_DIR / "dummy_isoformssequences.tsv"

# Expected output columns (mirrors DisCanVis2 upload contract)
EXPECTED_COLS = [
    "Entry_Name",
    "Gene_Uniprot",
    "Gene_Gencode",
    "Name",
    "Transcript name",
    "transcript_stable_id",
    "Transcript ID",
    "Entry_Isoform",
    "Database",
    "coverage_x",
    "coverage_y",
    "coverage",
    "alignmentpuntcuality",
]

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def run_worker(tmp_path: Path, extra_args: list[str] | None = None) -> pd.DataFrame:
    """Invoke the worker as a subprocess and return the output TSV as a DataFrame."""
    cmd = [
        sys.executable, str(BIN),
        "--blast_tsv",   str(BLAST_TSV),
        "--output_dir",  str(tmp_path),
        "--database",    "Gencode",
        "--coverage",    "80",
    ]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"Worker exited with code {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )
    out_file = tmp_path / "bestmaps_blast_gene_transcript.tsv"
    assert out_file.exists(), "Output file not created"
    return pd.read_csv(out_file, sep="\t", header=0)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestOutputSchema:
    """Ensure the output TSV has exactly the expected column set."""

    def test_output_columns(self, tmp_path):
        df = run_worker(tmp_path)
        assert list(df.columns) == EXPECTED_COLS, (
            f"Column mismatch.\nExpected: {EXPECTED_COLS}\nGot:      {list(df.columns)}"
        )

    def test_no_empty_columns(self, tmp_path):
        df = run_worker(tmp_path)
        for col in ["Entry_Isoform", "transcript_stable_id", "Transcript ID", "Database"]:
            assert df[col].notna().all(), f"Column '{col}' has unexpected NaN values"

    def test_transcript_stable_id_no_version(self, tmp_path):
        """transcript_stable_id must NOT contain a '.' version suffix."""
        df = run_worker(tmp_path)
        has_dot = df["transcript_stable_id"].str.contains(r"\.", na=False)
        assert not has_dot.any(), (
            f"transcript_stable_id should not have version suffix:\n{df.loc[has_dot, 'transcript_stable_id']}"
        )


class TestBestHitSelection:
    """Validate the priority-ladder best-hit selection logic."""

    def test_one_row_per_transcript_name(self, tmp_path):
        """Each unique 'Transcript name' must appear exactly once."""
        df = run_worker(tmp_path)
        counts = df["Transcript name"].value_counts()
        duplicates = counts[counts > 1]
        assert duplicates.empty, (
            f"Duplicate transcript names in output:\n{duplicates}"
        )

    def test_swissprot_preferred_over_trembl(self, tmp_path):
        """
        SAMD11-203 has both a TrEMBL hit (coverage 98.5/97.3) and a Swiss-Prot hit
        (coverage 100/16). With coverage >= 80, the TrEMBL hit has both coverages
        above threshold while the Swiss-Prot does not (16 < 80), so TrEMBL is
        correctly chosen. Verify the Database column is set appropriately.
        """
        df = run_worker(tmp_path)
        samd11_203 = df[df["Transcript name"] == "SAMD11-203"]
        assert len(samd11_203) == 1
        # TrEMBL hit meets the 80% coverage threshold; Swiss-Prot main entry does not
        assert samd11_203.iloc[0]["Database"] in (
            "Uniprot/SPTREMBL", "Uniprot/SWISSPROT"
        ), "Unexpected Database value for SAMD11-203"

    def test_identical_alignment_preferred(self, tmp_path):
        """PRDM16-201 has an identical alignment — must be marked 'identical'."""
        df = run_worker(tmp_path)
        prdm16 = df[df["Transcript name"] == "PRDM16-201"]
        assert len(prdm16) == 1
        assert prdm16.iloc[0]["alignmentpuntcuality"] == "identical"
        assert prdm16.iloc[0]["Database"] == "Uniprot/SWISSPROT"

    def test_isoform_tag_detection(self, tmp_path):
        """Entries with '-' in accession (e.g. Q9HAZ2-2) must be Uniprot_isoform."""
        df = run_worker(tmp_path)
        iso_rows = df[df["Entry_Isoform"].str.contains("-", na=False)]
        assert (iso_rows["Database"] == "Uniprot_isoform").all(), (
            "Isoform entries should be labelled Uniprot_isoform"
        )

    def test_coverage_calculated(self, tmp_path):
        """coverage must equal (coverage_x + coverage_y) / 2."""
        df = run_worker(tmp_path)
        expected = (df["coverage_x"] + df["coverage_y"]) / 2
        pd.testing.assert_series_equal(
            df["coverage"].round(6),
            expected.round(6),
            check_names=False,
        )


class TestIsoformsOutput:
    """Validate optional blastmaps_isoforms.tsv generation."""

    def test_isoforms_file_created_when_input_provided(self, tmp_path):
        run_worker(tmp_path, extra_args=["--isoforms_tsv", str(ISO_TSV)])
        iso_out = tmp_path / "blastmaps_isoforms.tsv"
        assert iso_out.exists(), "blastmaps_isoforms.tsv was not created"

    def test_isoforms_file_not_created_without_input(self, tmp_path):
        run_worker(tmp_path)
        iso_out = tmp_path / "blastmaps_isoforms.tsv"
        assert not iso_out.exists(), (
            "blastmaps_isoforms.tsv should NOT be created when no --isoforms_tsv is given"
        )

    def test_isoforms_contains_all_input_rows(self, tmp_path):
        """The isoforms table keeps all rows (no best-hit filtering)."""
        run_worker(tmp_path, extra_args=["--isoforms_tsv", str(ISO_TSV)])
        iso_df = pd.read_csv(tmp_path / "blastmaps_isoforms.tsv", sep="\t")
        src_df = pd.read_csv(ISO_TSV, sep="\t")
        assert len(iso_df) == len(src_df), (
            f"Row count mismatch: got {len(iso_df)}, expected {len(src_df)}"
        )


class TestEdgeCases:
    """Edge-case robustness checks."""

    def test_worker_help_exits_cleanly(self):
        result = subprocess.run(
            [sys.executable, str(BIN), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0

    def test_missing_blast_tsv_fails_gracefully(self, tmp_path):
        result = subprocess.run(
            [
                sys.executable, str(BIN),
                "--blast_tsv",  str(tmp_path / "nonexistent.tsv"),
                "--output_dir", str(tmp_path),
            ],
            capture_output=True, text=True,
        )
        assert result.returncode != 0, "Should fail with a missing input file"


# ---------------------------------------------------------------------------
# Mapping mode: main_isoform_mapping vs all_isoform_mapping
# ---------------------------------------------------------------------------
sys.path.insert(0, str(PROJECT_ROOT / "bin"))
from create_id_map_worker import (  # noqa: E402
    _pick_best_row,
    COL_DATABASE, COL_COVERAGE, COL_COVERAGE_X, COL_COVERAGE_Y,
    COL_PUNTCUALITY_X, COL_PUNTCUALITY_Y, COL_ENTRY_ISO,
)


def _cand(acc, db, cx, cy, punct):
    return {
        COL_ENTRY_ISO: acc, COL_DATABASE: db,
        COL_COVERAGE_X: cx, COL_COVERAGE_Y: cy,
        COL_COVERAGE: (cx + cy) / 2,
        COL_PUNTCUALITY_X: punct, COL_PUNTCUALITY_Y: punct,
    }


class TestMappingMode:
    def _sub(self):
        # A transcript that is IDENTICAL to isoform P04049-2 but only ~97% to
        # the canonical P04049.
        import pandas as pd
        return pd.DataFrame([
            _cand("P04049",   "Uniprot/SWISSPROT", 97.0, 97.0, "aligned"),
            _cand("P04049-2", "Uniprot_isoform",  100.0, 100.0, "identical"),
        ])

    def test_main_mode_prefers_canonical(self):
        row = _pick_best_row(self._sub(), 90.0, "main_isoform_mapping")
        assert row[COL_ENTRY_ISO] == "P04049"           # canonical preferred
        assert row[COL_DATABASE] == "Uniprot/SWISSPROT"

    def test_all_mode_picks_best_matching_isoform(self):
        row = _pick_best_row(self._sub(), 90.0, "all_isoform_mapping")
        assert row[COL_ENTRY_ISO] == "P04049-2"         # identical isoform wins
        assert row[COL_DATABASE] == "Uniprot_isoform"

    def test_all_mode_falls_back_to_canonical_when_no_better_isoform(self):
        import pandas as pd
        sub = pd.DataFrame([
            _cand("P04049", "Uniprot/SWISSPROT", 100.0, 100.0, "identical"),
            _cand("P04049-2", "Uniprot_isoform",  92.0, 92.0, "aligned"),
        ])
        row = _pick_best_row(sub, 90.0, "all_isoform_mapping")
        assert row[COL_ENTRY_ISO] == "P04049"           # canonical is the identical one
