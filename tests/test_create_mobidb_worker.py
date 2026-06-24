"""Tests for create_mobidb_worker.py (Module 5o — MobiDBDisorder)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_mobidb_worker.py"


def _run(args: list[str], tmpdir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(WORKER)] + args,
        capture_output=True, text=True, cwd=tmpdir,
    )


def _seq(tmpdir: Path) -> Path:
    p = tmpdir / "seq.tsv"
    p.write_text(
        "Protein_ID\tEntry_Isoform\tSequence\n"
        "GENE1-201\tP11111\t" + "A" * 100 + "\n"
        "GENE1-204\tP11111-2\t" + "A" * 80 + "\n"
        "GENE2-201\tP22222\t" + "A" * 50 + "\n",
        encoding="utf-8",
    )
    return p


def _mobidb(tmpdir: Path) -> Path:
    p = tmpdir / "mobidb.tsv"
    p.write_text(
        "acc\tfeature\tsource\tstart..end\n"
        "P11111\tcurated-disorder-merge\tmobidb_curated\t10-50\n"
        "P11111\tcurated-disorder-merge\tmobidb_curated\t70-90\n"
        "P22222\thomology-disorder-merge\tmobidb_homol\t1-20\n",
        encoding="utf-8",
    )
    return p


class TestBasicOutput:
    def test_output_file_created(self, tmp_path):
        r = _run(["--seq_table", str(_seq(tmp_path)),
                  "--mobidb_tsv", str(_mobidb(tmp_path)),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "mobidb_disorder.tsv").exists()

    def test_required_columns(self, tmp_path):
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        for col in ["Protein_ID", "Entry_Isoform", "feature",
                    "start_end", "content_fraction", "content_count", "length"]:
            assert col in df.columns, f"Missing: {col}"

    def test_protein_mapped(self, tmp_path):
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        assert "GENE1-201" in df["Protein_ID"].values
        assert "GENE2-201" in df["Protein_ID"].values


class TestRegionAggregation:
    def test_start_end_string_format(self, tmp_path):
        """Two regions for same protein → combined 'start-end,start-end' string."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["feature"] == "curated-disorder-merge")].iloc[0]
        assert "10-50" in str(row["start_end"])
        assert "70-90" in str(row["start_end"])

    def test_content_count(self, tmp_path):
        """Two regions → content_count = 2."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["feature"] == "curated-disorder-merge")].iloc[0]
        assert int(row["content_count"]) == 2

    def test_length_sum(self, tmp_path):
        """10-50 = 41 residues, 70-90 = 21 residues → total 62."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["feature"] == "curated-disorder-merge")].iloc[0]
        assert int(row["length"]) == 62

    def test_content_fraction(self, tmp_path):
        """62 disordered / 100 total = 0.62."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["feature"] == "curated-disorder-merge")].iloc[0]
        assert abs(float(row["content_fraction"]) - 0.62) < 0.01


class TestIsoformExpansion:
    def test_canonical_acc_maps_all_isoforms(self, tmp_path):
        """P11111 in MobiDB maps to both GENE1-201 (P11111) and GENE1-204 (P11111-2)."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        g1_pids = set(df[df["feature"] == "curated-disorder-merge"]["Protein_ID"].values)
        assert "GENE1-201" in g1_pids
        assert "GENE1-204" in g1_pids

    def test_content_fraction_per_isoform_length(self, tmp_path):
        """GENE1-204 has length 80, so content_fraction = 62/80."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(_mobidb(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-204") & (df["feature"] == "curated-disorder-merge")].iloc[0]
        expected = 62 / 80
        assert abs(float(row["content_fraction"]) - expected) < 0.01


class TestEdgeCases:
    def test_no_file_returns_empty(self, tmp_path):
        r = _run(["--seq_table", str(_seq(tmp_path)),
                  "--mobidb_tsv", str(tmp_path / "NO_FILE"),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        assert len(df) == 0

    def test_protein_not_in_mobidb_excluded(self, tmp_path):
        """Proteins not in MobiDB produce no rows (not empty rows)."""
        mob = tmp_path / "mob_partial.tsv"
        mob.write_text("acc\tfeature\tsource\tstart..end\nP22222\tcurated-disorder-merge\tx\t1-20\n",
                       encoding="utf-8")
        _run(["--seq_table", str(_seq(tmp_path)),
              "--mobidb_tsv", str(mob),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "mobidb_disorder.tsv", sep="\t")
        assert "GENE1-201" not in df["Protein_ID"].values
        assert "GENE2-201" in df["Protein_ID"].values
