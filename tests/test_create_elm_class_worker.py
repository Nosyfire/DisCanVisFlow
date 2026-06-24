"""Tests for create_elm_class_worker.py (Module 5n — ElmProteomeClassMatch)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_elm_class_worker.py"
ELM_CLASSES = Path(__file__).parent.parent / "legacy_data" / "elm" / "elm_classes-2025.tsv"


def _run(args: list[str], tmpdir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(WORKER)] + args,
        capture_output=True, text=True, cwd=tmpdir,
    )


def _make_classes(tmpdir: Path) -> Path:
    p = tmpdir / "elm_classes.tsv"
    p.write_text(
        '#ELM_Classes_Download_Version: 1.4\n'
        '#ELM_Classes_Download_Date: 2025-01-01\n'
        '"Accession"\t"ELMIdentifier"\t"FunctionalSiteName"\t"Description"\t"Regex"\t"Probability"\t"#Instances"\t"#Instances_in_PDB"\n'
        '"ELME000321"\t"CLV_C14_Caspase3-7"\t"Caspase cleavage motif"\t"Caspase-3 cleavage site."\t"[DSTE][^P]D[GSAN]"\t"0.0031"\t"41"\t"0"\n'
        '"ELME000001"\t"DEG_APCC_DBOX_1"\t"APC/C D-box"\t"D-box for APC/C."\t"R.L"\t"0.0012"\t"120"\t"20"\n'
        '"ELME000002"\t"LIG_SH2_GRB2"\t"GRB2 SH2 domain ligand"\t"GRB2 SH2 ligand."\t"Y.N"\t"0.0008"\t"55"\t"10"\n',
        encoding="utf-8",
    )
    return p


class TestBasicParsing:
    def test_output_file_created(self, tmp_path):
        elm = _make_classes(tmp_path)
        r = _run(["--elm_classes", str(elm), "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "elm_classes.tsv").exists()

    def test_row_count(self, tmp_path):
        elm = _make_classes(tmp_path)
        _run(["--elm_classes", str(elm), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        assert len(df) == 3

    def test_required_columns(self, tmp_path):
        elm = _make_classes(tmp_path)
        _run(["--elm_classes", str(elm), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        for col in ["elm_accession", "elm_identifier", "regex", "elm_type",
                    "probability", "n_instances", "n_instances_in_pdb"]:
            assert col in df.columns, f"Missing: {col}"

    def test_accession_values(self, tmp_path):
        elm = _make_classes(tmp_path)
        _run(["--elm_classes", str(elm), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        assert "ELME000321" in df["elm_accession"].values

    def test_elm_type_extracted(self, tmp_path):
        elm = _make_classes(tmp_path)
        _run(["--elm_classes", str(elm), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        row = df[df["elm_accession"] == "ELME000321"].iloc[0]
        assert row["elm_type"] == "CLV"
        row2 = df[df["elm_accession"] == "ELME000001"].iloc[0]
        assert row2["elm_type"] == "DEG"

    def test_probability_numeric(self, tmp_path):
        elm = _make_classes(tmp_path)
        _run(["--elm_classes", str(elm), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        row = df[df["elm_accession"] == "ELME000321"].iloc[0]
        assert abs(float(row["probability"]) - 0.0031) < 1e-5

    def test_n_instances_integer(self, tmp_path):
        elm = _make_classes(tmp_path)
        _run(["--elm_classes", str(elm), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        row = df[df["elm_accession"] == "ELME000321"].iloc[0]
        assert int(row["n_instances"]) == 41


class TestEdgeCases:
    def test_missing_file_returns_empty(self, tmp_path):
        r = _run(["--elm_classes", str(tmp_path / "no_file.tsv"),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        assert len(df) == 0

    def test_duplicate_accessions_deduplicated(self, tmp_path):
        p = tmp_path / "dup.tsv"
        p.write_text(
            '"Accession"\t"ELMIdentifier"\t"FunctionalSiteName"\t"Description"\t"Regex"\t"Probability"\t"#Instances"\t"#Instances_in_PDB"\n'
            '"ELME000001"\t"CLV_X_1"\t"X"\t"X"\t"X.X"\t"0.001"\t"1"\t"0"\n'
            '"ELME000001"\t"CLV_X_1"\t"X"\t"X"\t"X.X"\t"0.001"\t"1"\t"0"\n',
            encoding="utf-8",
        )
        _run(["--elm_classes", str(p), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        assert len(df) == 1


@pytest.mark.skipif(not ELM_CLASSES.exists(),
                    reason="elm_classes-2025.tsv not in legacy_data")
class TestRealFile:
    def test_real_file_produces_rows(self, tmp_path):
        r = _run(["--elm_classes", str(ELM_CLASSES), "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        assert len(df) > 300, f"Expected >300 ELM classes, got {len(df)}"

    def test_real_file_has_regex_column(self, tmp_path):
        _run(["--elm_classes", str(ELM_CLASSES), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elm_classes.tsv", sep="\t")
        assert df["regex"].notna().sum() > 300
