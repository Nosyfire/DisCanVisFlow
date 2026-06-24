"""Tests for create_snp_worker.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_snp_worker.py"


def _run(seq_table, snp_tsv, outdir):
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--seq_table",   str(seq_table),
         "--snp_pos_tsv", str(snp_tsv),
         "--outdir",      str(outdir)],
        capture_output=True, text=True
    )


def _make_seq(tmp, proteins):
    p = tmp / "seq.tsv"
    pd.DataFrame({"Protein_ID": proteins}).to_csv(p, sep="\t", index=False)
    return p


def _make_snp_pos(tmp, rows):
    p = tmp / "snp.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


class TestBasicOutput:
    def test_output_created(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        snp = _make_snp_pos(tmp_path, [
            {"AccessionPosition": "RAF1-201|270", "Polymorphism": "All Polymorphisms"}
        ])
        r = _run(seq, snp, tmp_path / "out")
        assert r.returncode == 0
        assert (tmp_path / "out" / "snp_polymorphisms.tsv").exists()

    def test_accession_position_split(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        snp = _make_snp_pos(tmp_path, [
            {"AccessionPosition": "RAF1-201|270", "Polymorphism": "Common Polymorphisms"}
        ])
        _run(seq, snp, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "snp_polymorphisms.tsv", sep="\t")
        assert "Protein_ID" in df.columns
        assert "Position" in df.columns
        assert df["Protein_ID"].iloc[0] == "RAF1-201"
        assert str(df["Position"].iloc[0]) == "270"


class TestFiltering:
    def test_only_target_proteins(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        snp = _make_snp_pos(tmp_path, [
            {"AccessionPosition": "RAF1-201|270", "Polymorphism": "All Polymorphisms"},
            {"AccessionPosition": "BRAF-201|100", "Polymorphism": "All Polymorphisms"},
        ])
        _run(seq, snp, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "snp_polymorphisms.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "RAF1-201"

    def test_common_and_all_kept(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        snp = _make_snp_pos(tmp_path, [
            {"AccessionPosition": "RAF1-201|100", "Polymorphism": "Common Polymorphisms"},
            {"AccessionPosition": "RAF1-201|200", "Polymorphism": "All Polymorphisms"},
        ])
        _run(seq, snp, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "snp_polymorphisms.tsv", sep="\t")
        assert len(df) == 2

    def test_missing_file_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        r = _run(seq, tmp_path / "nonexistent.tsv", tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "snp_polymorphisms.tsv", sep="\t")
        assert len(df) == 0
