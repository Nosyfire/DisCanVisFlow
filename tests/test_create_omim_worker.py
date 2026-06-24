"""Tests for create_omim_worker.py"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_omim_worker.py"


def _run(seq_table, omim_table, outdir):
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--seq_table", str(seq_table),
         "--omim_table", str(omim_table),
         "--outdir", str(outdir)],
        capture_output=True, text=True
    )


def _make_seq(tmp, proteins):
    p = tmp / "seq.tsv"
    pd.DataFrame({"Protein_ID": proteins}).to_csv(p, sep="\t", index=False)
    return p


def _make_omim(tmp, rows):
    p = tmp / "omim.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


class TestBasicOutput:
    def test_output_created(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        om = _make_omim(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "NSML", "MIMID": "611553"}
        ])
        r = _run(seq, om, tmp_path / "out")
        assert r.returncode == 0
        assert (tmp_path / "out" / "omim_disease.tsv").exists()

    def test_protein_id_renamed(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        om = _make_omim(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "NSML", "MIMID": "611553"}
        ])
        _run(seq, om, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "omim_disease.tsv", sep="\t")
        assert "Protein_ID" in df.columns
        assert "Accession" not in df.columns


class TestFiltering:
    def test_only_target_proteins(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        om = _make_omim(tmp_path, [
            {"Accession": "RAF1-201", "Disease": "NSML", "MIMID": "611553"},
            {"Accession": "BRAF-201", "Disease": "Other", "MIMID": "164757"},
        ])
        _run(seq, om, tmp_path / "out")
        df = pd.read_csv(tmp_path / "out" / "omim_disease.tsv", sep="\t")
        assert len(df) == 1
        assert df["Protein_ID"].iloc[0] == "RAF1-201"

    def test_missing_file_produces_empty(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        r = _run(seq, tmp_path / "nonexistent.tsv", tmp_path / "out")
        assert r.returncode == 0
        df = pd.read_csv(tmp_path / "out" / "omim_disease.tsv", sep="\t")
        assert len(df) == 0
"""Tests for create_omim_worker.py raw genemap2 parsing."""
import importlib.util, io, subprocess, sys
from pathlib import Path
import pandas as pd

BIN = Path(__file__).resolve().parent.parent / "bin" / "create_omim_worker.py"
spec = importlib.util.spec_from_file_location("create_omim_worker", BIN)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def _genemap2_line(symbols, approved, phenos):
    cols = [""] * 14
    cols[6] = symbols; cols[8] = approved; cols[12] = phenos
    return "\t".join(cols)


def test_parse_genemap2_extracts_name_and_mim():
    txt = io.StringIO(
        "# comment\n" +
        _genemap2_line("KRAS2, KRAS", "KRAS",
                       "Noonan syndrome 3, 609942 (3), Autosomal dominant; "
                       "{Leukemia}, 613065 (3)") + "\n"
    )
    rows = m.parse_genemap2(txt)
    genes = {g for g, _, _ in rows}
    assert genes == {"KRAS"}
    by_mim = {mim: name for _, name, mim in rows}
    assert by_mim["609942"] == "Noonan syndrome 3"
    assert by_mim["613065"] == "Leukemia"


def test_raw_mode_maps_gene_to_protein_id(tmp_path):
    seq = tmp_path / "seq.tsv"
    pd.DataFrame([{"Protein_ID": "KRAS-201", "Gene": "KRAS", "Entry_Isoform": "P01116"},
                  {"Protein_ID": "KRAS-202", "Gene": "KRAS", "Entry_Isoform": "P01116-2"}]
                 ).to_csv(seq, sep="\t", index=False)
    raw = tmp_path / "omim_raw"; raw.mkdir()
    (raw / "genemap2.txt").write_text(
        _genemap2_line("KRAS", "KRAS", "Noonan syndrome 3, 609942 (3)") + "\n")
    subprocess.run([sys.executable, str(BIN), "--mapping_mode", "raw",
                    "--seq_table", str(seq), "--omim_raw_dir", str(raw),
                    "--outdir", str(tmp_path)], check=True, capture_output=True)
    df = pd.read_csv(tmp_path / "omim_disease.tsv", sep="\t", dtype=str)
    assert set(df["Protein_ID"]) == {"KRAS-201", "KRAS-202"}
    assert df["Disease"].iloc[0] == "Noonan syndrome 3"
    assert set(df["MIMID"]) == {"609942"}


_HUMSAVAR = """\
        UniProt header junk
Main        Swiss-Prot             AA             Variant
gene name   AC         FTId        change         category dbSNP          Disease name
_________   __________ ___________ ______________ ________ ______________ _____________________
KRAS        P01116     VAR_006840  p.Gly12Asp     LP/P     rs121913529    Noonan syndrome 3 [MIM:609942]
KRAS        P01116     VAR_018347  p.His52Arg     LB/B     rs893184       -
TP53        P04637     VAR_044503  p.Arg175His    LP/P     rs28934578     Li-Fraumeni syndrome
"""


def test_parse_humsavar_block():
    rows = m.parse_humsavar(io.StringIO(_HUMSAVAR))
    # all three data rows parsed
    assert len(rows) == 3
    kras = rows[0]
    assert kras["acc"] == "P01116"
    assert kras["pos"] == "12"
    assert kras["wt"] == "G"
    assert kras["category"] == "LP/P"
    assert kras["dbSNP"] == "rs121913529"
    assert kras["disease"].startswith("Noonan syndrome 3")
    assert kras["MIMID"] == "609942"
    # benign row has no disease
    assert rows[1]["disease"] == ""


def test_raw_mode_humsavar_mutations_and_fanout(tmp_path):
    # canonical KRAS-201 (P01116) with G at pos12; alt KRAS-205 (P01116-2) with a
    # 2-residue N-term extension so the context window shifts the variant.
    seq = tmp_path / "seq.tsv"
    pd.DataFrame([
        {"Protein_ID": "KRAS-201", "Entry_Isoform": "P01116", "Gene": "KRAS",
         "Sequence": "MTEYKLVVVG", "main_isoform": "yes"},
        {"Protein_ID": "KRAS-205", "Entry_Isoform": "P01116-2", "Gene": "KRAS",
         "Sequence": "XXMTEYKLVVVG", "main_isoform": "no"},
    ]).to_csv(seq, sep="\t", index=False)
    # variant at pos 5 = K (MTEYK...), disease-named
    hum = tmp_path / "humsavar.txt"
    hum.write_text(
        "_________\n"
        "KRAS        P01116     VAR_000001  p.Lys5Asn      LP/P     rs1   Noonan syndrome [MIM:163950]\n"
        "KRAS        P01116     VAR_000002  p.Thr2Ala      LB/B     rs2   -\n"
    )
    subprocess.run([sys.executable, str(BIN), "--mapping_mode", "raw",
                    "--seq_table", str(seq), "--humsavar", str(hum),
                    "--outdir", str(tmp_path)], check=True, capture_output=True)
    mut = pd.read_csv(tmp_path / "omim_mutations.tsv", sep="\t", dtype=str)
    by_pid = {r["Protein_ID"]: r for _, r in mut.iterrows()}
    # benign (no disease) variant excluded; only the disease variant maps
    assert set(mut["aa_change"]) == {"p.Lys5Asn"}
    # canonical: direct at pos 5
    assert by_pid["KRAS-201"]["Protein_position"] == "5"
    assert by_pid["KRAS-201"]["MIMID"] == "163950"
    assert by_pid["KRAS-201"]["dbSNP"] == "rs1"
    # alt isoform: context window shifts K to pos 7
    assert by_pid["KRAS-205"]["Protein_position"] == "7"
    # disease table derived too
    dis = pd.read_csv(tmp_path / "omim_disease.tsv", sep="\t", dtype=str)
    assert set(dis["Protein_ID"]) == {"KRAS-201", "KRAS-205"}
    assert dis["Disease"].str.startswith("Noonan syndrome").all()
