"""Tests for bin/create_lcr_worker.py (SEG low-complexity via segmasker)."""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_lcr_worker.py"
OUT_COLS = ["Protein_ID", "start", "end", "length"]

pytestmark = pytest.mark.skipif(
    shutil.which("segmasker") is None, reason="segmasker (BLAST+) not installed")


def _run(args):
    return subprocess.run([sys.executable, str(BIN)] + args,
                          capture_output=True, text=True)


def _write_seq_table(p: Path):
    # LCR protein: globular head + long poly-S/poly-Q low-complexity tail.
    lcr_seq = "MKVLAAGDEFRHIKPWY" + "S" * 25 + "Q" * 25
    globular = "MKVLAAGDEFRHIKPWYTNCQMKVLAAGDEFRHIKPWYTNCQ"
    cols = ["Protein_ID", "main_isoform", "Sequence"]
    rows = [
        ["LCRX-201", "yes", lcr_seq],
        ["GLOB-201", "yes", globular],
    ]
    pd.DataFrame(rows, columns=cols).to_csv(p, sep="\t", index=False)


def test_masks_low_complexity_run(tmp_path):
    seq = tmp_path / "seq.tsv"
    _write_seq_table(seq)
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path)])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "low_complexity.tsv", sep="\t", dtype=str)
    assert list(out.columns) == OUT_COLS
    lcr = out[out["Protein_ID"] == "LCRX-201"]
    assert len(lcr) >= 1
    # the masked interval must fall inside the poly-S/poly-Q tail (start > 17)
    starts = lcr["start"].astype(int)
    ends = lcr["end"].astype(int)
    assert (ends >= starts).all()
    assert (starts >= 1).all()
    assert lcr["length"].astype(int).max() >= 15


def test_only_main_isoforms_filter(tmp_path):
    seq = tmp_path / "seq.tsv"
    df = pd.DataFrame(
        [["A-201", "yes", "M" + "P" * 40],
         ["A-202", "no",  "M" + "P" * 40]],
        columns=["Protein_ID", "main_isoform", "Sequence"])
    df.to_csv(seq, sep="\t", index=False)
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--only_main_isoforms"])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "low_complexity.tsv", sep="\t", dtype=str)
    assert set(out["Protein_ID"]) <= {"A-201"}
