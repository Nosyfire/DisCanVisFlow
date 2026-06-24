"""Tests for create_elm_switches_worker.py (Module 5p — Elm_Switches)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_elm_switches_worker.py"


def _run(args, tmpdir):
    return subprocess.run(
        [sys.executable, str(WORKER)] + args,
        capture_output=True, text=True, cwd=tmpdir,
    )


def _seq(tmpdir):
    p = tmpdir / "seq.tsv"
    p.write_text(
        "Protein_ID\tEntry_Isoform\tGene_Gencode\tSequence\n"
        "LAT-201\tO43561\tLAT\t" + "MDEFGHIKLMNPQRSTVWY" * 10 + "\n"
        "LAT-210\tO43561-2\tLAT\t" + "MDEFGHIKLMNPQRSTVWY" * 9 + "\n"
        "AKT1-206\tP31749\tAKT1\t" + "ACDEFGHIKLMNPQRSTVWY" * 25 + "\n",
        encoding="utf-8",
    )
    return p


def _switches(tmpdir):
    p = tmpdir / "elmswitches.tsv"
    p.write_text(
        "Switch ID\tStatus\tInteraction ID\tIntramolecular\t"
        "ID A\tBindingsite A ID\tBindingsite A Start\tBindingsite A End\t"
        "ID B\tBindingsite B ID\tBindingsite B Start\tBindingsite B End\t"
        "Affected interactor\tSwitch type\tSwitch subtype\tSwitch mechanism\t"
        "Switch direction\tSwitch outcome direction\tSwitch outcome\t"
        "Modification\tModification sites\tModifying enzymes\tEffector\t"
        "Cell cycle phase\tLocalisation\tPathway\tPMID\n"
        "SWTI000001\tActive\tINTI000001\t\t"
        "UNIPROT:O43561\tELM:LIG_SH2_STAT5\t5\t19\t"
        "UNIPROT:P19174\tPFAM:PF00017\t550\t639\t"
        "ID A\tBinary\tPhysicochemical compatibility\tPTM-dependent\t"
        "Reversible\tPositive\tInduction\tMOD:00696\tY161\t\t\t"
        "\tGO:0009898\tKEGG:hsa04660\tPMID: 20610546\n",
        encoding="utf-8",
    )
    return p


class TestBasicOutput:
    def test_output_file_created(self, tmp_path):
        r = _run(["--seq_table", str(_seq(tmp_path)),
                  "--switches", str(_switches(tmp_path)),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "elmswitches_mapped.tsv").exists()

    def test_required_columns(self, tmp_path):
        _run(["--seq_table", str(_seq(tmp_path)),
              "--switches", str(_switches(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        for col in ["Protein_ID", "Entry_Isoform", "homology_transfer",
                    "Switch ID", "Bindingsite A Start", "Bindingsite A End"]:
            assert col in df.columns, f"Missing: {col}"

    def test_canonical_isoform_mapped(self, tmp_path):
        _run(["--seq_table", str(_seq(tmp_path)),
              "--switches", str(_switches(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        assert "LAT-201" in df["Protein_ID"].values


class TestHomologyTransfer:
    def test_canonical_no_transfer(self, tmp_path):
        """Direct match for canonical isoform → homology_transfer=False."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--switches", str(_switches(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        row = df[df["Protein_ID"] == "LAT-201"].iloc[0]
        assert row["homology_transfer"] == False

    def test_other_isoform_mapped(self, tmp_path):
        """LAT-210 (O43561-2) should be mapped via homology transfer."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--switches", str(_switches(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        lat_pids = set(df["Protein_ID"].values)
        assert "LAT-210" in lat_pids

    def test_homology_transfer_true_for_isoform(self, tmp_path):
        """Isoform row should have homology_transfer=True."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--switches", str(_switches(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        isoform_rows = df[df["Entry_Isoform"] == "O43561-2"]
        if not isoform_rows.empty:
            assert isoform_rows.iloc[0]["homology_transfer"] == True

    def test_coordinates_preserved_for_canonical(self, tmp_path):
        """Canonical isoform keeps original Start=5, End=19."""
        _run(["--seq_table", str(_seq(tmp_path)),
              "--switches", str(_switches(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        row = df[df["Protein_ID"] == "LAT-201"].iloc[0]
        assert int(row["Bindingsite A Start"]) == 5
        assert int(row["Bindingsite A End"]) == 19


class TestEdgeCases:
    def test_no_file_returns_empty(self, tmp_path):
        r = _run(["--seq_table", str(_seq(tmp_path)),
                  "--switches", str(tmp_path / "NO_FILE"),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        assert len(df) == 0

    def test_protein_not_in_seq_skipped(self, tmp_path):
        """Switch for unknown accession produces no output rows."""
        sw = tmp_path / "unknown.tsv"
        sw.write_text(
            "Switch ID\tStatus\tInteraction ID\tIntramolecular\t"
            "ID A\tBindingsite A ID\tBindingsite A Start\tBindingsite A End\t"
            "ID B\tBindingsite B ID\tBindingsite B Start\tBindingsite B End\t"
            "Affected interactor\tSwitch type\tSwitch subtype\tSwitch mechanism\t"
            "Switch direction\tSwitch outcome direction\tSwitch outcome\t"
            "Modification\tModification sites\tModifying enzymes\tEffector\t"
            "Cell cycle phase\tLocalisation\tPathway\tPMID\n"
            "SWTI999\tActive\tINTI999\t\tUNIPROT:ZZZZZZ\t\t1\t10\t"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n",
            encoding="utf-8",
        )
        _run(["--seq_table", str(_seq(tmp_path)),
              "--switches", str(sw),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        assert len(df) == 0

    def test_real_file_if_available(self, tmp_path):
        """Smoke test against the real ELM switches file."""
        switches_path = (Path(__file__).parent.parent /
                         "legacy_data" / "elm" / "elmswitches-2023.tsv")
        if not switches_path.exists():
            pytest.skip("Real ELM switches file not available")
        seq_path = switches_path.parent.parent.parent / "legacy_data" / "elm" / "elmswitches-2023.tsv"
        # Just test parsing doesn't crash; use a minimal seq_table
        r = _run(["--seq_table", str(_seq(tmp_path)),
                  "--switches", str(switches_path),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "elmswitches_mapped.tsv", sep="\t")
        assert "Switch ID" in df.columns
