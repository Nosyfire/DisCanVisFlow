"""Tests for bin/create_exon_worker.py — exon boundary extraction from combined_map.map."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_exon_worker.py"


def _write_map(path: Path):
    # Two "exons": residues 0-2 contiguous, then a >20bp genomic jump at residue 3.
    # Header is the GENCODE pipe-delimited format; Protein_ID is pipe field 4.
    lines = [
        "# ENST00000000001.1|ENSG00000000001.1|OTT1|-|GENE-201|GENE|999|CDS:1-9| chr1 - 100-200",
        "0 M 1,2,3 ATG M 1000,999,998, ATG M",
        "1 E 4,5,6 GAG E 997,996,995, GAG E",
        "2 H 7,8,9 CAC H 994,993,992, CAC H",
        "3 I 10,11,12 ATA I 900,899,898, ATA I",   # genomic jump 992 -> 900 (>20)
        "4 Q 13,14,15 CAG Q 897,896,895, CAG Q",
    ]
    path.write_text("\n".join(lines) + "\n")


def _write_loc(path: Path):
    pd.DataFrame({"Protein_ID": ["GENE-201"]}).to_csv(path, sep="\t", index=False)


def test_exon_extraction(tmp_path):
    mp = tmp_path / "combined_map.map"
    loc = tmp_path / "loc.tsv"
    _write_map(mp)
    _write_loc(loc)

    subprocess.run(
        [sys.executable, str(BIN), "--combined_map", str(mp),
         "--loc_chrom", str(loc), "--output_dir", str(tmp_path)],
        check=True,
    )

    df = pd.read_csv(tmp_path / "exon.tsv", sep="\t", dtype=str)
    assert list(df.columns) == ["Protein_ID", "exon_number", "total_exons",
                                "aa_start", "aa_end", "aa_length",
                                "genomic_start", "genomic_end"]
    # Two exons: aa 1-3 and aa 4-5 (1-based), all for GENE-201
    assert len(df) == 2
    assert (df["Protein_ID"] == "GENE-201").all()
    assert (df["total_exons"] == "2").all()
    e1, e2 = df.iloc[0], df.iloc[1]
    assert e1["exon_number"] == "1"
    assert (e1["aa_start"], e1["aa_end"], e1["aa_length"]) == ("1", "3", "3")
    assert e2["exon_number"] == "2"
    assert (e2["aa_start"], e2["aa_end"], e2["aa_length"]) == ("4", "5", "2")
    # genomic boundary coordinate of exon 1 start = first codon's first nt
    assert e1["genomic_start"] == "1000"


def test_protein_id_keyed_lookup(tmp_path):
    """Block must be findable by Protein_ID (pipe field 4), not just the ENST id."""
    mp = tmp_path / "combined_map.map"
    loc = tmp_path / "loc.tsv"
    _write_map(mp)
    _write_loc(loc)
    subprocess.run(
        [sys.executable, str(BIN), "--combined_map", str(mp),
         "--loc_chrom", str(loc), "--output_dir", str(tmp_path)],
        check=True,
    )
    df = pd.read_csv(tmp_path / "exon.tsv", sep="\t", dtype=str)
    assert not df.empty, "Protein_ID lookup failed — exon.tsv is empty"
