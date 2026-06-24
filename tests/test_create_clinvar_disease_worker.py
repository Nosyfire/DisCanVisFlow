"""Tests for create_clinvar_disease_worker.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_clinvar_disease_worker.py"


def _run(seq_table, clinvar_disease, outdir, cat_tsv=None):
    no_file = Path(__file__).parent.parent / "assets" / "NO_FILE"
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--seq_table",           str(seq_table),
         "--clinvar_disease",     str(clinvar_disease),
         "--clinvar_category_tsv", str(cat_tsv) if cat_tsv else str(no_file),
         "--outdir",              str(outdir)],
        capture_output=True, text=True
    )


def _make_seq(tmp, proteins):
    p = tmp / "seq.tsv"
    pd.DataFrame({"Protein_ID": proteins, "seq": ["ACDE"] * len(proteins)}).to_csv(p, sep="\t", index=False)
    return p


def _make_clinvar(tmp, rows):
    p = tmp / "clinvar_disease.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _make_cat(tmp, rows):
    p = tmp / "cat.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


class TestBasicOutput:
    def test_output_created(self, tmp_path):
        seq = _make_seq(tmp_path, ["P1-201"])
        cv = _make_clinvar(tmp_path, [
            {"Accession": "P1-201", "Disease": "CancerX", "DOID": "DOID:123"}
        ])
        r = _run(seq, cv, tmp_path / "out")
        assert r.returncode == 0
        assert (tmp_path / "out" / "clinvar_disease.tsv").exists()

    def test_protein_id_renamed(self, tmp_path):
        seq = _make_seq(tmp_path, ["P1-201"])
        cv = _make_clinvar(tmp_path, [
            {"Accession": "P1-201", "Disease": "CancerX", "DOID": "DOID:1"}
        ])
        _run(seq, cv, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "clinvar_disease.tsv", sep="\t")
        assert "Protein_ID" in df.columns
        assert "Accession" not in df.columns


class TestFiltering:
    def test_only_target_proteins(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        cv = _make_clinvar(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "NSML", "DOID": "DOID:1"},
            {"Accession": "BRAF-201", "Disease": "Other", "DOID": "DOID:2"},
        ])
        _run(seq, cv, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "clinvar_disease.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "RAF1-201"

    def test_no_match_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, ["UNKN-201"])
        cv = _make_clinvar(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "NSML", "DOID": "DOID:1"},
        ])
        _run(seq, cv, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "clinvar_disease.tsv", sep="\t")
        assert len(df) == 0

    def test_missing_clinvar_file_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        r = _run(seq, tmp_path / "nonexistent.tsv", tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "clinvar_disease.tsv", sep="\t")
        assert len(df) == 0


class TestCategoryJoin:
    def test_final_category_added(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        cv = _make_clinvar(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "Noonan syndrome", "DOID": "DOID:3490"},
        ])
        cat = _make_cat(tmp_path, [
            {"Disease": "Noonan syndrome", "Final_Category": "Cardiovascular/Hematopoietic",
             "disease_category": "Cardiovascular/Hematopoietic"},
        ])
        _run(seq, cv, tmp_path / "out", cat_tsv=cat)
        df = pd.read_csv(tmp_path / "out" / "clinvar_disease.tsv", sep="\t")
        assert "Final_Category" in df.columns
        assert df["Final_Category"].iloc[0] == "Cardiovascular/Hematopoietic"

    def test_unknown_when_no_category_match(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        cv = _make_clinvar(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "UnknownDisease", "DOID": "DOID:999"},
        ])
        cat = _make_cat(tmp_path, [
            {"Disease": "Noonan syndrome", "Final_Category": "Cardiovascular/Hematopoietic",
             "disease_category": "Cardiovascular/Hematopoietic"},
        ])
        _run(seq, cv, tmp_path / "out", cat_tsv=cat)
        df = pd.read_csv(tmp_path / "out" / "clinvar_disease.tsv", sep="\t")
        assert df["Final_Category"].iloc[0] == "Unknown"

    def test_no_category_file_still_works(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        cv = _make_clinvar(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "NSML", "DOID": "DOID:1"},
        ])
        r = _run(seq, cv, tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "clinvar_disease.tsv", sep="\t")
        assert df["Final_Category"].iloc[0] == "Unknown"
