"""Tests for the dbNSFP pack + slice tooling (bin/dbnsfp_pack.sh + slice_dbnsfp.py).

Skipped when bgzip (htslib) is not on PATH.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).parent.parent / "bin"
PACK = BIN / "dbnsfp_pack.sh"
SLICE = BIN / "slice_dbnsfp.py"

pytestmark = pytest.mark.skipif(shutil.which("bgzip") is None,
                                reason="bgzip (htslib) not installed")


def _make_bundle(tmp_path: Path) -> Path:
    src = tmp_path / "in.tsv"
    src.write_text(
        "Protein_ID\tProtein_position\tREVEL_score\tgnomAD4.1_joint_AF\n"
        "BRAF-201\t5\t0.9\t0.01\n"
        "RAF1-201\t10\t0.5\t0.02\n"
        "BRAF-201\t2\t0.8\t0.03\n"       # unsorted + interleaved on purpose
        "RAF1-201\t3\t0.4\t0.05\n"
        "RAF1-204\t7\t0.6\t0.09\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    env = dict(os.environ, SORT_BUF="1G", SORT_PAR="2", TMPDIR=str(tmp_path / "st"))
    r = subprocess.run(["bash", str(PACK), str(src), str(out)],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    return out / "dbnsfp_scores.tsv.gz"


def _slice(bgz: Path, ids: str, extra=None):
    args = [sys.executable, str(SLICE), "--bgz", str(bgz), "--id", ids] + (extra or [])
    r = subprocess.run(args, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout, r.stderr


def test_pidx_offsets(tmp_path):
    bgz = _make_bundle(tmp_path)
    pidx = (bgz.parent / "dbnsfp_scores.pidx").read_text().splitlines()
    rows = dict((p.split("\t")[0], p.split("\t")[1:]) for p in pidx)
    assert set(rows) == {"BRAF-201", "RAF1-201", "RAF1-204"}
    # BRAF sorts first → offset 0
    assert rows["BRAF-201"][0] == "0"


def test_slice_single_sorted_by_position(tmp_path):
    bgz = _make_bundle(tmp_path)
    out, _ = _slice(bgz, "RAF1-201")
    lines = out.strip().split("\n")
    assert lines[0].startswith("Protein_ID")          # header included
    assert lines[1].split("\t")[:2] == ["RAF1-201", "3"]   # residue 3 before 10
    assert lines[2].split("\t")[:2] == ["RAF1-201", "10"]
    assert len(lines) == 3


def test_slice_multiple_no_header(tmp_path):
    bgz = _make_bundle(tmp_path)
    out, _ = _slice(bgz, "BRAF-201,RAF1-204", extra=["--no_header"])
    pids = {ln.split("\t")[0] for ln in out.strip().split("\n")}
    assert pids == {"BRAF-201", "RAF1-204"}
    assert "Protein_ID" not in out


def test_slice_does_not_bleed_into_next_protein(tmp_path):
    """A slice returns only the requested protein's rows (offset+length exact)."""
    bgz = _make_bundle(tmp_path)
    out, _ = _slice(bgz, "BRAF-201", extra=["--no_header"])
    assert all(ln.startswith("BRAF-201\t") for ln in out.strip().split("\n"))


def test_missing_id_warns_not_fatal(tmp_path):
    bgz = _make_bundle(tmp_path)
    out, err = _slice(bgz, "NOPE-201", extra=["--no_header"])
    assert out.strip() == ""
    assert "not in index" in err
