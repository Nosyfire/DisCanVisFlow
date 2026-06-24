"""Tests for run_gopher_worker.py — conservation scoring + taxonomic levels."""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parent.parent / "bin" / "run_gopher_worker.py"
spec = importlib.util.spec_from_file_location("run_gopher_worker", BIN)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_species_of():
    assert m.species_of("P04637_HUMAN") == "HUMAN"
    assert m.species_of("sp|P04637|P53_MOUSE extra") == "MOUSE"
    assert m.species_of("Q9XYZ1_DROME") == "DROME"


def test_column_conservation_bounds():
    assert m.column_conservation(["A", "A", "A"]) == 1.0          # identical
    assert m.column_conservation(["-", "-"]) == 0.0               # all gap
    # more diversity → strictly lower conservation than a mostly-conserved column
    diverse = m.column_conservation(["A", "C", "D", "E"])
    mostly = m.column_conservation(["A", "A", "A", "C"])
    assert 0.0 <= diverse < mostly < 1.0
    # a gap fraction lowers the score relative to the same residues with no gap
    assert m.column_conservation(["A", "A", "-"]) < m.column_conservation(["A", "A"])


def test_score_alignment_query_residues_only():
    seqs = [("Q_HUMAN", "MA-K"), ("o1_MOUSE", "MACK"), ("o2_CHICK", "MA-R")]
    scores = m.score_alignment(seqs, 0, [0, 1, 2])
    # query 'MA-K' has 3 non-gap residues (M,A,K) → 3 scores
    assert len(scores) == 3
    # column0 (M,M,M) fully conserved
    assert scores[0] == 1.0


def test_members_for_level_partitions_by_taxon():
    seqs = [("Q_HUMAN", "MK"), ("a_MOUSE", "MK"), ("b_DROME", "MR"), ("c_ARATH", "MN")]
    taxon = m.DEFAULT_TAXON
    mam = m.members_for_level(seqs, 0, "Mammalia", taxon)
    assert set(mam) == {0, 1}                       # human + mouse
    euk = m.members_for_level(seqs, 0, "Eukaryota", taxon)
    assert set(euk) == {0, 1, 2, 3}                 # all are eukaryotes
    plant = m.members_for_level(seqs, 0, "Viridiplantae", taxon)
    assert set(plant) == {0, 3}                     # query + arabidopsis
    glob = m.members_for_level(seqs, 0, "global", taxon)
    assert set(glob) == {0, 1, 2, 3}


def test_end_to_end_table(tmp_path):
    seq = tmp_path / "seq.tsv"
    pd.DataFrame([{"Protein_ID": "TP53-201", "Entry_Isoform": "P04637",
                   "Gene": "TP53", "Sequence": "MEK"}]).to_csv(seq, sep="\t", index=False)
    aln = tmp_path / "aln"; aln.mkdir()
    (aln / "P04637.orthaln.fas").write_text(
        ">P04637_HUMAN\nMEK\n>o1_MOUSE\nMEK\n>o2_DROME\nMDK\n")
    subprocess.run([sys.executable, str(BIN), "--seq_table", str(seq),
                    "--aln_dir", str(aln), "--out", str(tmp_path / "ct.tsv")],
                   check=True, capture_output=True)
    df = pd.read_csv(tmp_path / "ct.tsv", sep="\t", dtype=str)
    assert list(df.columns) == ["uniprot_acc", "level", "conservation_score"]
    assert set(df["uniprot_acc"]) == {"P04637"}
    # one row per level
    assert "global" in set(df["level"]) and "Mammalia" in set(df["level"])
    glob = df[df["level"] == "global"].iloc[0]["conservation_score"]
    # query MEK = 3 residues → 3 comma-separated scores
    assert len(glob.split(",")) == 3
