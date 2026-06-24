"""Tests for create_alphamissense_worker.py (raw GENCODE gz and plain TSV formats)."""
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_alphamissense_worker.py"

_HEADER = "# AlphaMissense_isoforms_hg38\n#CHROM\tPOS\tREF\tALT\tgenome\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n"


def _run(seq_table, am_gz, outdir):
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--seq_table", str(seq_table),
         "--alphamissense_gz", str(am_gz),
         "--outdir", str(outdir)],
        capture_output=True, text=True
    )


def _make_seq(tmp, rows):
    """rows: list of dicts with at least Protein_ID and transcript_stable_id."""
    p = tmp / "seq.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _make_gz(tmp, data_rows):
    """data_rows: list of dicts with raw AM columns (transcript_id has .version suffix)."""
    p = tmp / "am.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_HEADER)
        for r in data_rows:
            fh.write("\t".join([
                r.get("CHROM", "chr3"),
                r.get("POS", "12625100"),
                r.get("REF", "A"),
                r.get("ALT", "T"),
                r.get("genome", "hg38"),
                r["transcript_id"],
                r["protein_variant"],
                r["am_pathogenicity"],
                r["am_class"],
            ]) + "\n")
    return p


def _make_tsv(tmp, data_rows):
    """Plain (decompressed) TSV — the format produced by DECOMPRESS_ALPHAMISSENSE."""
    p = tmp / "am.tsv"
    with open(p, "wt") as fh:
        fh.write(_HEADER)
        for r in data_rows:
            fh.write("\t".join([
                r.get("CHROM", "chr3"),
                r.get("POS", "12625100"),
                r.get("REF", "A"),
                r.get("ALT", "T"),
                r.get("genome", "hg38"),
                r["transcript_id"],
                r["protein_variant"],
                r["am_pathogenicity"],
                r["am_class"],
            ]) + "\n")
    return p


class TestBasicOutput:
    def test_output_created(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        gz = _make_gz(tmp_path, [
            {"transcript_id": "ENST00000251849.2",
             "protein_variant": "V1A", "am_pathogenicity": "0.3", "am_class": "likely_benign"}
        ])
        r = _run(seq, gz, tmp_path / "out")
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "out" / "alphamissense.tsv").exists()

    def test_required_columns_present(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        gz = _make_gz(tmp_path, [
            {"transcript_id": "ENST00000251849.2",
             "protein_variant": "V1A", "am_pathogenicity": "0.3", "am_class": "likely_benign"}
        ])
        _run(seq, gz, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        for col in ["Protein_ID", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"]:
            assert col in df.columns


class TestFiltering:
    def test_only_target_proteins(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        gz = _make_gz(tmp_path, [
            {"transcript_id": "ENST00000251849.2",
             "protein_variant": "V1A", "am_pathogenicity": "0.3", "am_class": "likely_benign"},
            {"transcript_id": "ENST00000288602.5",   # BRAF, not in seq
             "protein_variant": "V1A", "am_pathogenicity": "0.8", "am_class": "likely_pathogenic"},
        ])
        _run(seq, gz, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "RAF1-201"

    def test_missing_file_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        r = _run(seq, tmp_path / "nonexistent.tsv.gz", tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert len(df) == 0

    def test_version_suffix_stripped(self, tmp_path):
        """transcript_id 'ENST00000251849.12' should match 'ENST00000251849'."""
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        gz = _make_gz(tmp_path, [
            {"transcript_id": "ENST00000251849.12",
             "protein_variant": "A100T", "am_pathogenicity": "0.9", "am_class": "likely_pathogenic"}
        ])
        _run(seq, gz, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "RAF1-201"

    def test_multiple_variants_per_protein(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        rows = [
            {"transcript_id": "ENST00000251849.2",
             "protein_variant": f"V{i}A", "am_pathogenicity": "0.3", "am_class": "likely_benign"}
            for i in range(1, 6)
        ]
        gz = _make_gz(tmp_path, rows)
        _run(seq, gz, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert len(df) == 5

    def test_multiple_isoforms_mapped(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"},
            {"Protein_ID": "RAF1-202", "transcript_stable_id": "ENST00000399990"},
        ])
        gz = _make_gz(tmp_path, [
            {"transcript_id": "ENST00000251849.2",
             "protein_variant": "V1A", "am_pathogenicity": "0.3", "am_class": "likely_benign"},
            {"transcript_id": "ENST00000399990.4",
             "protein_variant": "V1A", "am_pathogenicity": "0.4", "am_class": "likely_benign"},
        ])
        _run(seq, gz, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert set(df["Protein_ID"]) == {"RAF1-201", "RAF1-202"}

    def test_no_enst_column_raises(self, tmp_path):
        """seq_table without transcript_stable_id or 'Transcript ID' should fail."""
        seq = tmp_path / "seq.tsv"
        pd.DataFrame({"Protein_ID": ["RAF1-201"]}).to_csv(seq, sep="\t", index=False)
        gz = _make_gz(tmp_path, [
            {"transcript_id": "ENST00000251849.2",
             "protein_variant": "V1A", "am_pathogenicity": "0.3", "am_class": "likely_benign"}
        ])
        r = _run(seq, gz, tmp_path / "out")
        assert r.returncode != 0


class TestPlainTsv:
    """Verify the grep pre-filter path (DECOMPRESS_ALPHAMISSENSE output)."""

    def test_plain_tsv_matches_gz_output(self, tmp_path):
        """Plain TSV and gzip versions must produce identical output."""
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"},
            {"Protein_ID": "RAF1-202", "transcript_stable_id": "ENST00000399990"},
        ])
        rows = [
            {"transcript_id": "ENST00000251849.2",
             "protein_variant": "V1A", "am_pathogenicity": "0.3", "am_class": "likely_benign"},
            {"transcript_id": "ENST00000399990.4",
             "protein_variant": "K5R", "am_pathogenicity": "0.7", "am_class": "ambiguous"},
            {"transcript_id": "ENST00000288602.5",  # not in seq_table
             "protein_variant": "A1T", "am_pathogenicity": "0.9", "am_class": "likely_pathogenic"},
        ]
        gz  = _make_gz(tmp_path, rows)
        tsv = _make_tsv(tmp_path, rows)

        _run(seq, gz,  tmp_path / "out_gz")
        _run(seq, tsv, tmp_path / "out_tsv")

        df_gz  = pd.read_csv(tmp_path / "out_gz"  / "alphamissense.tsv", sep="\t")
        df_tsv = pd.read_csv(tmp_path / "out_tsv" / "alphamissense.tsv", sep="\t")

        pd.testing.assert_frame_equal(
            df_gz.sort_values(list(df_gz.columns)).reset_index(drop=True),
            df_tsv.sort_values(list(df_tsv.columns)).reset_index(drop=True),
        )

    def test_plain_tsv_filters_correctly(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "TP53-201", "transcript_stable_id": "ENST00000269305"}
        ])
        tsv = _make_tsv(tmp_path, [
            {"transcript_id": "ENST00000269305.9",
             "protein_variant": "R175H", "am_pathogenicity": "0.99", "am_class": "likely_pathogenic"},
            {"transcript_id": "ENST00000288602.5",
             "protein_variant": "V1A",  "am_pathogenicity": "0.1",  "am_class": "likely_benign"},
        ])
        _run(seq, tsv, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "TP53-201"
        assert df["protein_variant"].iloc[0] == "R175H"

    def test_plain_tsv_version_suffix_stripped(self, tmp_path):
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        tsv = _make_tsv(tmp_path, [
            {"transcript_id": "ENST00000251849.12",
             "protein_variant": "A100T", "am_pathogenicity": "0.9", "am_class": "likely_pathogenic"}
        ])
        _run(seq, tsv, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert len(df) == 1

    def test_plain_tsv_no_matches_returns_empty(self, tmp_path):
        """grep finds nothing → empty output, no crash."""
        seq = _make_seq(tmp_path, [
            {"Protein_ID": "RAF1-201", "transcript_stable_id": "ENST00000251849"}
        ])
        tsv = _make_tsv(tmp_path, [
            {"transcript_id": "ENST00000288602.5",   # not RAF1
             "protein_variant": "V1A", "am_pathogenicity": "0.1", "am_class": "likely_benign"},
        ])
        r = _run(seq, tsv, tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "alphamissense.tsv", sep="\t")
        assert len(df) == 0
