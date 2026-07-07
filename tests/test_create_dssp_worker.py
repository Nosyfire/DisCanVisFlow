# tests/test_create_dssp_worker.py
"""Tests for bin/create_dssp_worker.py (DSSP SS + true RSA from local mmCIF)."""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_dssp_worker.py"
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "dssp"
OUT_COLS = ["Protein_ID", "Position", "aa", "ss8", "ss3", "rsa"]

pytestmark = pytest.mark.skipif(
    shutil.which("mkdssp") is None, reason="mkdssp not installed")


def _run(args):
    return subprocess.run([sys.executable, str(BIN)] + args,
                          capture_output=True, text=True)


def test_parses_local_cif(tmp_path):
    # Entry_Isoform P0TEST -> fixture file AF-P0TEST-F1.cif in cif_dir.
    # Fixture is the first 8 residues (MVLSPADK) of the real AlphaFold model
    # for human hemoglobin subunit alpha (UniProt P69905, AF-P69905-F1),
    # trimmed offline — see tests/fixtures/dssp/_make_fixture.py.
    seq = tmp_path / "seq.tsv"
    pd.DataFrame(
        [["PT-201", "P0TEST", "yes", "MVLSPADK"]],
        columns=["Protein_ID", "Entry_Isoform", "main_isoform", "Sequence"]
    ).to_csv(seq, sep="\t", index=False)

    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--cif_dir", str(FIXTURE_DIR), "--no_download"])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "dssp.tsv", sep="\t", dtype=str)
    assert list(out.columns) == OUT_COLS
    pt = out[out["Protein_ID"] == "PT-201"]
    assert len(pt) == 8                       # one row per residue
    assert list(pt["Position"].astype(int)) == [1, 2, 3, 4, 5, 6, 7, 8]
    assert "".join(pt.sort_values("Position")["aa"]) == "MVLSPADK"
    assert set(pt["ss3"]) <= {"H", "E", "C"}
    rsa = pt["rsa"].astype(float)
    assert (rsa >= 0).all() and (rsa <= 1.6).all()


def test_sanitize_cif_strips_ma_categories():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "dsspw", str(Path(__file__).resolve().parents[1] / "bin" / "create_dssp_worker.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    cif = (
        "data_AFTEST\n#\n"
        "_atom_site.group_PDB\nATOM\n#\n"
        "_ma_qa_metric.id\n1\n#\n"
        "_ma_qa_metric_local.label_asym_id\nA\n#\n"
        "_entity_poly_seq.mon_id\nMET\n#\n"
    )
    out = mod._sanitize_cif(cif)
    assert "_ma_qa_metric" not in out
    assert "_ma_qa_metric_local" not in out
    assert "_atom_site" in out
    assert "_entity_poly_seq" in out
    assert out.startswith("data_AFTEST")


def test_missing_structure_skipped(tmp_path):
    seq = tmp_path / "seq.tsv"
    pd.DataFrame(
        [["NOPE-201", "Q0NONE", "yes", "MKAAQR"]],
        columns=["Protein_ID", "Entry_Isoform", "main_isoform", "Sequence"]
    ).to_csv(seq, sep="\t", index=False)
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--cif_dir", str(FIXTURE_DIR), "--no_download"])
    assert res.returncode == 0, res.stderr
    out = pd.read_csv(tmp_path / "dssp.tsv", sep="\t", dtype=str)
    assert list(out.columns) == OUT_COLS
    assert len(out) == 0                       # no structure -> skipped, header only
