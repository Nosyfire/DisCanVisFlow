# tests/test_create_plaac_worker.py
"""Tests for bin/create_plaac_worker.py (PLAAC prion-like domain HMM)."""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin" / "create_plaac_worker.py"
JAR = ROOT / "External_Programs" / "plaac" / "plaac.jar"
OUT_COLS = ["Protein_ID", "Position", "plaac_score", "in_PRD"]

pytestmark = pytest.mark.skipif(
    not JAR.exists() or shutil.which("java") is None,
    reason="PLAAC jar or Java not available")


def _run(args):
    return subprocess.run([sys.executable, str(BIN)] + args,
                          capture_output=True, text=True)


def test_prd_called_in_qn_rich(tmp_path):
    seq = tmp_path / "seq.tsv"
    prd = "MSDSNQGNNQQNYQQYSQNGNQQQGNNRYQGYQAYNAQAQPAGGYYQNYQGYSGYQQGGYQ"
    glob = "MKVLAAGDEFRHIKPWYTNCQ"
    pd.DataFrame(
        [["PRD-201", "yes", prd], ["GLOB-201", "yes", glob]],
        columns=["Protein_ID", "main_isoform", "Sequence"]
    ).to_csv(seq, sep="\t", index=False)

    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--plaac_jar", str(JAR)])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "plaac.tsv", sep="\t", dtype=str)
    assert list(out.columns) == OUT_COLS

    prd_rows = out[out["Protein_ID"] == "PRD-201"]
    glob_rows = out[out["Protein_ID"] == "GLOB-201"]
    # one row per residue, 1-based contiguous
    assert prd_rows["Position"].astype(int).tolist() == list(range(1, len(prd) + 1))
    assert glob_rows["Position"].astype(int).tolist() == list(range(1, len(glob) + 1))
    # the Q/N-rich prion protein has in-PRD residues; the globular control does not
    assert (prd_rows["in_PRD"].astype(int) == 1).any()
    assert (glob_rows["in_PRD"].astype(int) == 0).all()


def test_only_main_isoforms_filter(tmp_path):
    seq = tmp_path / "seq.tsv"
    pd.DataFrame(
        [["A-201", "yes", "MKVLAAGDEFRHIKPWYTNCQ"],
         ["A-202", "no",  "MKVLAAGDEFRHIKPWYTNCQ"]],
        columns=["Protein_ID", "main_isoform", "Sequence"]
    ).to_csv(seq, sep="\t", index=False)
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--plaac_jar", str(JAR), "--only_main_isoforms"])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "plaac.tsv", sep="\t", dtype=str)
    assert set(out["Protein_ID"]) <= {"A-201"}
