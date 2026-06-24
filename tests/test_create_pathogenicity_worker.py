"""Tests for create_pathogenicity_worker.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_pathogenicity_worker.py"


def _run(seq_table, dbnsfp_tsv, outdir):
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--seq_table", str(seq_table),
         "--dbnsfp_tsv", str(dbnsfp_tsv),
         "--outdir", str(outdir)],
        capture_output=True, text=True
    )


def _make_seq(tmp, proteins):
    p = tmp / "seq.tsv"
    pd.DataFrame({"Protein_ID": proteins}).to_csv(p, sep="\t", index=False)
    return p


def _make_dbnsfp(tmp, rows):
    p = tmp / "dbnsfp.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


_SAMPLE_ROW = {
    "Protein_ID": "RAF1-201",
    "chr": "3",
    "Start_Position": "12590959",
    "End_Position": "12590960",
    "Protein_position": "403",
    "aaref": "V",
    "aaalt": "A",
    "aapos": "403",
    "ref": "T",
    "alt": "C",
    "rs_dbSNP": "rs12345",
    "AlphaMissense_score": "0.42",
    "CADD_phred": "22.3",
    "CADD_raw": "2.1",
    "ClinPred_score": "0.7",
    "ESM1b_score": "-1.2",
    "EVE_score": "0.5",
    "Polyphen2_HDIV_score": "0.9",
    "Polyphen2_HVAR_score": "0.85",
    "PrimateAI_score": "0.6",
    "SIFT_score": "0.01",
    "VARITY_ER_LOO_score": "0.8",
    "VARITY_R_LOO_score": "0.75",
    "REVEL_score": "0.65",
    "gMVP_score": "0.55",
}


class TestBasicOutput:
    def test_output_created(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        db = _make_dbnsfp(tmp_path, [_SAMPLE_ROW])
        r = _run(seq, db, tmp_path / "out")
        assert r.returncode == 0
        assert (tmp_path / "out" / "pathogenicity_scores.tsv").exists()

    def test_predictor_columns_present(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        db = _make_dbnsfp(tmp_path, [_SAMPLE_ROW])
        _run(seq, db, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "pathogenicity_scores.tsv", sep="\t")
        for col in ["AlphaMissense_score", "CADD_phred", "SIFT_score",
                    "Polyphen2_HDIV_score", "REVEL_score", "gMVP_score"]:
            assert col in df.columns


class TestFiltering:
    def test_only_target_proteins(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        other = dict(_SAMPLE_ROW, Protein_ID="BRAF-201")
        db = _make_dbnsfp(tmp_path, [_SAMPLE_ROW, other])
        _run(seq, db, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "pathogenicity_scores.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "RAF1-201"

    def test_missing_file_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        r = _run(seq, tmp_path / "nonexistent.tsv", tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "pathogenicity_scores.tsv", sep="\t")
        assert len(df) == 0

    def test_present_file_without_protein_id_column(self, tmp_path):
        # A non-empty dbNSFP file lacking Protein_ID must degrade to an empty
        # output, not crash with KeyError mid-stream.
        seq = _make_seq(tmp_path, ["RAF1-201"])
        db = tmp_path / "dbnsfp.tsv"
        pd.DataFrame([
            {"chr": "3", "aaref": "V", "AlphaMissense_score": "0.4"},
        ]).to_csv(db, sep="\t", index=False)
        r = _run(seq, db, tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "pathogenicity_scores.tsv", sep="\t")
        assert len(df) == 0


class TestMultipleVariants:
    def test_multiple_variants_same_position(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        rows = [dict(_SAMPLE_ROW, aaalt=aa) for aa in ["A", "G", "L", "P"]]
        db = _make_dbnsfp(tmp_path, rows)
        _run(seq, db, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "pathogenicity_scores.tsv", sep="\t")
        assert len(df) == 4
