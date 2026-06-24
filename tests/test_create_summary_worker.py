"""Tests for create_summary_worker.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_summary_worker.py"


def _run(gene, results_dir, outdir, extra=None):
    cmd = [sys.executable, str(WORKER),
           "--gene_name",   gene,
           "--results_dir", str(results_dir),
           "--outdir",      str(outdir)]
    if extra:
        cmd += extra
    return subprocess.run(cmd, capture_output=True, text=True)


def _make_tsv(path: Path, n_rows: int, cols=("Protein_ID",)):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({c: ["RAF1-201"] * n_rows for c in cols}).to_csv(path, sep="\t", index=False)


class TestBasicOutput:
    def test_output_created(self, tmp_path):
        r = _run("RAF1", tmp_path, tmp_path)
        assert r.returncode == 0
        assert (tmp_path / "annotation_summary.tsv").exists()

    def test_columns(self, tmp_path):
        _run("RAF1", tmp_path, tmp_path)
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        for col in ["gene", "annotation_type", "count", "note"]:
            assert col in df.columns

    def test_gene_name_stored(self, tmp_path):
        _run("RAF1", tmp_path, tmp_path)
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        assert (df["gene"] == "RAF1").all()


class TestCounting:
    def test_elm_count(self, tmp_path):
        # final/annotations/elm.tsv is the Protein_ID-keyed output (label "ELM motifs")
        elm = tmp_path / "final" / "annotations" / "elm.tsv"
        _make_tsv(elm, 4, ("Protein_ID", "Start", "End"))
        _run("RAF1", tmp_path, tmp_path)
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        row = df[df["annotation_type"] == "ELM motifs"]
        assert row["count"].iloc[0] == 4

    def test_raw_elm_count(self, tmp_path):
        # unfinal/annotations/elm.tsv is the Entry_Isoform-keyed raw input (label "ELM motifs (raw)")
        elm = tmp_path / "intermediate" / "annotations" / "elm.tsv"
        _make_tsv(elm, 6, ("Protein_ID", "Start", "End"))
        _run("RAF1", tmp_path, tmp_path)
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        row = df[df["annotation_type"] == "ELM motifs (raw)"]
        assert row["count"].iloc[0] == 6

    def test_zero_when_file_missing(self, tmp_path):
        _run("RAF1", tmp_path, tmp_path)
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        row = df[df["annotation_type"] == "ELM motifs"]
        assert row["count"].iloc[0] == 0

    def test_file_override(self, tmp_path):
        custom_elm = tmp_path / "my_elm.tsv"
        _make_tsv(custom_elm, 7, ("Protein_ID",))
        _run("RAF1", tmp_path, tmp_path,
             extra=["--file", "ELM motifs=" + str(custom_elm)])
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        row = df[df["annotation_type"] == "ELM motifs"]
        assert row["count"].iloc[0] == 7

    def test_mutations_counted(self, tmp_path):
        mut = tmp_path / "final" / "mutations" / "ClinVar" / "Missense_filter_mutations_mapped.tsv"
        _make_tsv(mut, 5, ("Protein_ID", "HGVSp_Short"))
        _run("RAF1", tmp_path, tmp_path)
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        row = df[df["annotation_type"] == "Missense mutations"]
        assert row["count"].iloc[0] == 5

    def test_conservation_counted(self, tmp_path):
        cons = tmp_path / "final" / "conservation" / "conservation_multiple_level.tsv"
        _make_tsv(cons, 3, ("Protein_ID", "Entry_Isoform", "level", "conservationscores"))
        _run("RAF1", tmp_path, tmp_path)
        df = pd.read_csv(tmp_path / "annotation_summary.tsv", sep="\t")
        row = df[df["annotation_type"] == "GOPHER conservation entries"]
        assert row["count"].iloc[0] == 3
