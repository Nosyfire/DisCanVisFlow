# tests/test_create_catgranule_worker.py
"""Tests for bin/create_catgranule_worker.py (catGRANULE 2.0 LLMS)."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_catgranule_worker.py"
OUT_COLS = ["Protein_ID", "Position", "catgranule_score", "catgranule_total"]
CATG_PY = Path("/opt/anaconda3/envs/catgranule/bin/python")


def _run(args):
    return subprocess.run([sys.executable, str(BIN)] + args,
                          capture_output=True, text=True)


def _seq_table(p: Path):
    pd.DataFrame(
        [["FUS-201", "yes", "MKKGGYYQNQGGYSGYQQGGYQSNYGQQSYGGGGQQGNNRYGGYQ"],
         ["GLOB-201", "yes", "MKVLAAGDEFRHIKPWYTNCQ"]],
        columns=["Protein_ID", "main_isoform", "Sequence"]
    ).to_csv(p, sep="\t", index=False)


def test_graceful_skip_when_unavailable(tmp_path):
    # No --catgranule_python and catGRANULE not importable in the test interpreter
    # (discanvis/test env lacks matplotlib) → header-only output, exit 0.
    seq = tmp_path / "seq.tsv"
    _seq_table(seq)
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--catgranule_lib", "/nonexistent/catgranule"])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "catgranule.tsv", sep="\t", dtype=str)
    assert list(out.columns) == OUT_COLS
    assert len(out) == 0


@pytest.mark.skipif(not CATG_PY.exists(),
                    reason="catgranule env not built")
def test_real_compute_via_env_python(tmp_path):
    seq = tmp_path / "seq.tsv"
    _seq_table(seq)
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--catgranule_python", str(CATG_PY), "--only_main_isoforms"])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "catgranule.tsv", sep="\t", dtype=str)
    assert list(out.columns) == OUT_COLS
    assert len(out) > 0
    for pid, grp in out.groupby("Protein_ID"):
        # one row per residue, 1-based contiguous
        positions = grp["Position"].astype(int).tolist()
        assert positions[0] == 1
        assert positions == list(range(1, len(positions) + 1))
        # total is constant within a protein and numeric
        assert grp["catgranule_total"].nunique() == 1
        grp["catgranule_score"].astype(float)      # must parse
    # FUS-like (45 aa) should have more residues than the short globular control
    lens = out.groupby("Protein_ID").size()
    assert lens["FUS-201"] == 45
