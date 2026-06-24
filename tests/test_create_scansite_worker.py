"""Tests for create_scansite_worker.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_scansite_worker.py"


def _run(seq_table, scansite_tsv, outdir, extra=None):
    cmd = [sys.executable, str(WORKER),
           "--seq_table",    str(seq_table),
           "--scansite_tsv", str(scansite_tsv),
           "--outdir",       str(outdir)]
    if extra:
        cmd += extra
    return subprocess.run(cmd, capture_output=True, text=True)


def _make_seq(tmp, rows):
    p = tmp / "seq.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _make_scansite(tmp, rows):
    p = tmp / "scansite.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


_SS_ROW = {
    "Protein_ID": "RAF1-201", "motifName": "PKC delta",
    "motifShortName": "PKC_d", "score": "0.301",
    "site": "S259", "siteSequence": "RLSSSSVGSSEDASStT",
    "Start": "255", "End": "271",
}


class TestPrecomputed:
    def test_output_created(self, tmp_path):
        seq = _make_seq(tmp_path, [{"Protein_ID": "RAF1-201", "Sequence": "MAAAA"}])
        ss  = _make_scansite(tmp_path, [_SS_ROW])
        r   = _run(seq, ss, tmp_path / "out")
        assert r.returncode == 0
        assert (tmp_path / "out" / "scansite.tsv").exists()

    def test_columns_present(self, tmp_path):
        seq = _make_seq(tmp_path, [{"Protein_ID": "RAF1-201", "Sequence": "MAAAA"}])
        ss  = _make_scansite(tmp_path, [_SS_ROW])
        _run(seq, ss, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "scansite.tsv", sep="\t")
        for col in ["Protein_ID", "motifName", "motifShortName", "score", "site", "Start", "End"]:
            assert col in df.columns

    def test_filtered_to_target_proteins(self, tmp_path):
        seq = _make_seq(tmp_path, [{"Protein_ID": "RAF1-201", "Sequence": "MAAAA"}])
        other = dict(_SS_ROW, Protein_ID="BRAF-201")
        ss  = _make_scansite(tmp_path, [_SS_ROW, other])
        _run(seq, ss, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "scansite.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "RAF1-201"

    def test_no_file_sentinel_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, [{"Protein_ID": "RAF1-201", "Sequence": "MAAAA"}])
        no_file = tmp_path / "NO_FILE"
        no_file.write_text("")
        r = _run(seq, no_file, tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "scansite.tsv", sep="\t")
        assert len(df) == 0

    def test_missing_file_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, [{"Protein_ID": "RAF1-201", "Sequence": "MAAAA"}])
        r = _run(seq, tmp_path / "nonexistent.tsv", tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "scansite.tsv", sep="\t")
        assert len(df) == 0
