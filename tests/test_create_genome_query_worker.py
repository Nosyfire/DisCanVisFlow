"""Tests for bin/create_genome_query_worker.py — genome↔protein reference tables."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_genome_query_worker.py"


def _write_minus_map(path: Path):
    # Minus strand: codon in coding order, genomic positions decreasing.
    # ATG(M) at 12618720,19,18  →  genome plus-strand ref = complement(A,T,G)=T,A,C
    lines = [
        "# ENST00000000001.1|ENSG1|OTT1|-|GENE-201|GENE|999|CDS:1-6| chr3 - 100-200",
        "0 M 1,2,3 ATG M 12618720,12618719,12618718, ATG M",
        "1 F 4,5,6 TTC F 12584518,12584517,12584516, TTC F",
    ]
    path.write_text("\n".join(lines) + "\n")


def _write_plus_map(path: Path):
    # Plus strand: genomic ref base == coding base.
    lines = [
        "# ENST00000000002.1|ENSG2|OTT2|+|GENE-301|GENE|999|CDS:1-6| chr1 + 100-200",
        "0 M 1,2,3 ATG M 100,101,102, ATG M",
        "1 F 4,5,6 TTC F 103,104,105, TTC F",
    ]
    path.write_text("\n".join(lines) + "\n")


def _run(mp, outdir):
    subprocess.run([sys.executable, str(BIN), "--combined_map", str(mp),
                    "--outdir", str(outdir)], check=True)


def test_index_is_per_nucleotide(tmp_path):
    mp = tmp_path / "combined_map.map"
    _write_minus_map(mp)
    _run(mp, tmp_path)
    idx = pd.read_csv(tmp_path / "genome_protein_index.tsv", sep="\t", dtype=str)
    assert list(idx.columns) == ["chrom", "gpos", "strand", "Protein_ID",
                                 "prot_pos", "codon_offset", "aa", "codon"]
    assert len(idx) == 6                      # 2 residues × 3 nucleotides
    first = idx.iloc[0]
    assert first["prot_pos"] == "1" and first["aa"] == "M"
    assert first["chrom"] == "chr3" and first["strand"] == "-"


def test_all_snvs_enumerated(tmp_path):
    mp = tmp_path / "combined_map.map"
    _write_minus_map(mp)
    _run(mp, tmp_path)
    m = pd.read_csv(tmp_path / "genome_protein_mutations.tsv", sep="\t", dtype=str)
    # 2 residues × 3 nucleotides × 3 alternative bases = 18 reference SNVs
    assert len(m) == 18
    assert list(m.columns) == ["chrom", "gpos", "strand", "ref", "alt",
                               "Protein_ID", "prot_pos", "codon_offset",
                               "ref_codon", "alt_codon", "ref_aa", "alt_aa",
                               "consequence", "hgvs_g", "hgvs_p"]


def test_minus_strand_ref_is_complement(tmp_path):
    mp = tmp_path / "combined_map.map"
    _write_minus_map(mp)
    _run(mp, tmp_path)
    m = pd.read_csv(tmp_path / "genome_protein_mutations.tsv", sep="\t", dtype=str)
    # Codon offset 1 of residue 1: coding base A, genome ref = complement(A) = T
    first_nuc = m[(m["prot_pos"] == "1") & (m["codon_offset"] == "1")]
    assert set(first_nuc["ref"]) == {"T"}
    assert set(first_nuc["alt"]) == {"A", "C", "G"}


def test_plus_strand_ref_equals_coding(tmp_path):
    mp = tmp_path / "combined_map.map"
    _write_plus_map(mp)
    _run(mp, tmp_path)
    m = pd.read_csv(tmp_path / "genome_protein_mutations.tsv", sep="\t", dtype=str)
    first_nuc = m[(m["prot_pos"] == "1") & (m["codon_offset"] == "1")]
    assert set(first_nuc["ref"]) == {"A"}    # coding base A == genome ref on + strand


def test_consequence_and_hgvs(tmp_path):
    mp = tmp_path / "combined_map.map"
    _write_plus_map(mp)
    _run(mp, tmp_path)
    m = pd.read_csv(tmp_path / "genome_protein_mutations.tsv", sep="\t", dtype=str)
    # Residue 1 = M (ATG). ATG->ACG = Thr (missense), ATG->ATA/ATC = Ile (missense)
    row = m[(m["prot_pos"] == "1") & (m["alt_codon"] == "ACG")].iloc[0]
    assert row["ref_aa"] == "M" and row["alt_aa"] == "T"
    assert row["consequence"] == "missense"
    assert row["hgvs_p"] == "p.M1T"
    assert row["hgvs_g"].startswith("g.101")   # codon offset 2, gpos 101
    # synonymous example: TTC (F) -> TTT (F)
    syn = m[(m["prot_pos"] == "2") & (m["alt_codon"] == "TTT")].iloc[0]
    assert syn["consequence"] == "synonymous"


def test_nonsense_detected(tmp_path):
    # Codon TAC (Y) -> TAA / TAG = stop  => nonsense
    mp = tmp_path / "combined_map.map"
    lines = [
        "# ENST3|ENSG3|OTT3|+|GENE-401|GENE|9|CDS:1-3| chr1 + 1-9",
        "0 Y 1,2,3 TAC Y 10,11,12, TAC Y",
    ]
    mp.write_text("\n".join(lines) + "\n")
    _run(mp, tmp_path)
    m = pd.read_csv(tmp_path / "genome_protein_mutations.tsv", sep="\t", dtype=str)
    nonsense = m[m["consequence"] == "nonsense"]
    assert len(nonsense) >= 1
    assert set(nonsense["alt_aa"]) == {"*"}
