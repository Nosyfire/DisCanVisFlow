"""Tests for create_depmap_worker.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_depmap_worker.py"


def _run(seq_table, depmap_tsv, outdir):
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--seq_table", str(seq_table),
         "--depmap_tsv", str(depmap_tsv),
         "--outdir", str(outdir)],
        capture_output=True, text=True
    )


def _make_seq(tmp, proteins):
    p = tmp / "seq.tsv"
    pd.DataFrame({"Protein_ID": proteins}).to_csv(p, sep="\t", index=False)
    return p


def _make_depmap(tmp, rows):
    p = tmp / "depmap.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


_SAMPLE_ROW = {
    "Protein_ID": "RAF1-201",
    "Chrom": "chr3",
    "Start_Position": "12590959",
    "End_Position": "12590960",
    "HugoSymbol": "RAF1",
    "Protein_position": "403",
    "HGVSp_Short": "p.V403A",
    "VariantType": "SNV",
    "VariantInfo": "missense_variant",
    "DNAChange": "c.1208T>C",
    "ModelID": "ACH-000936",
    "Hotspot": "False",
    "EntrezGeneID": "5894",
    "Rescue": "False",
    "RescueReason": "",
}


class TestBasicOutput:
    def test_output_created(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        dm = _make_depmap(tmp_path, [_SAMPLE_ROW])
        r = _run(seq, dm, tmp_path / "out")
        assert r.returncode == 0
        assert (tmp_path / "out" / "depmap_mutations.tsv").exists()

    def test_required_columns(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        dm = _make_depmap(tmp_path, [_SAMPLE_ROW])
        _run(seq, dm, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        for col in ["Protein_ID", "HGVSp_Short", "ModelID", "VariantType"]:
            assert col in df.columns


def _make_seq_full(tmp, rows):
    """seq_table with the full schema the gene-isoform lookup needs."""
    p = tmp / "seq_full.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


class TestRawGeneKeyedMode:
    """fetch_depmap_worker.py emits HugoSymbol/HGVSp_Short with no Protein_ID;
    the worker must map gene+variant onto every run isoform."""

    def _seq(self, tmp):
        # two isoforms of TP53 sharing the region around residue 5
        return _make_seq_full(tmp, [
            {"Protein_ID": "TP53-201", "Entry_Isoform": "P04637",
             "Gene": "TP53", "Gene_Gencode": "TP53",
             "Sequence": "MEEPQSDPSVEPPLSQ"},
            {"Protein_ID": "TP53-204", "Entry_Isoform": "P04637-2",
             "Gene": "TP53", "Gene_Gencode": "TP53",
             "Sequence": "XXMEEPQSDPSVEPPLSQ"},  # +2 N-terminal shift
        ])

    def test_raw_maps_to_all_isoforms(self, tmp_path):
        seq = self._seq(tmp_path)
        # p.Q5K → residue 5 on TP53-201 is Q (MEEP*Q*), should map to both
        dm = _make_depmap(tmp_path, [{
            "HugoSymbol": "TP53", "Protein_position": "5",
            "HGVSp_Short": "p.Q5K", "ModelID": "ACH-000001",
            "Start_Position": "7676000", "EntrezGeneID": "7157", "Hotspot": "False",
        }])
        r = _run(seq, dm, tmp_path / "out")
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        assert "Protein_ID" in df.columns
        pids = set(df["Protein_ID"])
        assert "TP53-201" in pids and "TP53-204" in pids
        # direct hit on -201 at pos 5, context-transferred to -204 at pos 7
        d201 = df[df["Protein_ID"] == "TP53-201"].iloc[0]
        d204 = df[df["Protein_ID"] == "TP53-204"].iloc[0]
        assert str(d201["Protein_position"]) == "5"
        assert str(d201["isoform_mapped"]) == "False"
        assert str(d204["Protein_position"]) == "7"
        assert str(d204["isoform_mapped"]) == "True"

    def test_raw_wt_mismatch_not_forced(self, tmp_path):
        seq = self._seq(tmp_path)
        # claim p.A5K but residue 5 is Q on both → no direct hit, context (A at
        # centre) won't match either → variant dropped
        dm = _make_depmap(tmp_path, [{
            "HugoSymbol": "TP53", "Protein_position": "5",
            "HGVSp_Short": "p.A5K", "ModelID": "ACH-000002",
            "Start_Position": "7676000", "EntrezGeneID": "7157", "Hotspot": "False",
        }])
        r = _run(seq, dm, tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        assert len(df) == 0

    def test_raw_other_gene_filtered(self, tmp_path):
        seq = self._seq(tmp_path)
        dm = _make_depmap(tmp_path, [{
            "HugoSymbol": "BRAF", "Protein_position": "600",
            "HGVSp_Short": "p.V600E", "ModelID": "ACH-000003",
            "Start_Position": "140753336", "EntrezGeneID": "673", "Hotspot": "True",
        }])
        r = _run(seq, dm, tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        assert len(df) == 0

    def test_raw_three_letter_hgvsp(self, tmp_path):
        seq = self._seq(tmp_path)
        dm = _make_depmap(tmp_path, [{
            "HugoSymbol": "TP53", "Protein_position": "5",
            "HGVSp_Short": "p.Gln5Lys", "ModelID": "ACH-000004",
            "Start_Position": "7676000", "EntrezGeneID": "7157", "Hotspot": "False",
        }])
        r = _run(seq, dm, tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        assert "TP53-201" in set(df["Protein_ID"])


class TestFiltering:
    def test_only_target_proteins(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        other = dict(_SAMPLE_ROW, Protein_ID="BRAF-201")
        dm = _make_depmap(tmp_path, [_SAMPLE_ROW, other])
        _run(seq, dm, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        assert (df["Protein_ID"] == "RAF1-201").all()
        assert len(df) == 1

    def test_multiple_cell_lines(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        rows = [dict(_SAMPLE_ROW, ModelID=f"ACH-{i:06d}") for i in range(5)]
        dm = _make_depmap(tmp_path, rows)
        _run(seq, dm, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        assert len(df) == 5

    def test_missing_file_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        r = _run(seq, tmp_path / "nonexistent.tsv", tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "depmap_mutations.tsv", sep="\t")
        assert len(df) == 0
