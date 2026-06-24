#!/usr/bin/env python3
"""Tests for bin/create_mavedb_worker.py.

Covers the original premapped mode (unchanged behaviour) and the new
uniprot mode (maps UniProt accession + protein position onto the run's
Gencode isoforms with direct / homology_similarity fan-out).
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).resolve().parent.parent / "bin" / "create_mavedb_worker.py"

OUT_COLS = ["Protein_ID", "Protein_position", "prot_expr", "score",
            "mavedb_id", "urn", "gene_name", "uniprot", "Transcript_ID",
            "is_double_mutant", "mapping_type"]


def _run(args):
    res = subprocess.run([sys.executable, str(BIN)] + args,
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    return res


# ---------------------------------------------------------------------------
# premapped mode (existing behaviour — must keep passing)
# ---------------------------------------------------------------------------

def _write_seq_table(path, rows):
    """rows: list of dicts → tsv"""
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_premapped_filters_by_protein_id(tmp_path):
    seq = tmp_path / "seq.tsv"
    _write_seq_table(seq, [{"Protein_ID": "RAF1-201"}, {"Protein_ID": "RAF1-205"}])

    mave = tmp_path / "mave.tsv"
    pd.DataFrame({
        "Protein_ID": ["RAF1-201", "OTHER-201"],
        "prot_expr": ["p.Arg89Tyr", "p.Lys10Glu"],
        "protein_start": ["89", "10"],
        "score": ["1.2", "0.4"],
        "mavedb_id": ["urn:x#1", "urn:y#2"],
        "urn": ["urn:x", "urn:y"],
        "gene_name": ["RAF1", "OTH"],
        "uniprot": ["P04049", "Q9"],
        "is_double_mutant": ["False", "False"],
    }).to_csv(mave, sep="\t", index=False)

    _run(["--seq_table", str(seq), "--mavedb", str(mave),
          "--outdir", str(tmp_path)])

    df = pd.read_csv(tmp_path / "mavedb.tsv", sep="\t", dtype=str)
    assert list(df.columns) == OUT_COLS
    assert set(df["Protein_ID"]) == {"RAF1-201"}
    assert df.iloc[0]["Protein_position"] == "89"
    assert df.iloc[0]["mapping_type"] == "direct"


def test_premapped_default_mode(tmp_path):
    """Without --mapping_mode the worker still does premapped filtering."""
    seq = tmp_path / "seq.tsv"
    _write_seq_table(seq, [{"Protein_ID": "G-201"}])
    mave = tmp_path / "mave.tsv"
    pd.DataFrame({
        "Protein_ID": ["G-201"], "prot_expr": ["p.Arg1Cys"],
        "protein_start": ["1"], "score": ["9"], "mavedb_id": ["u#1"],
        "urn": ["u"], "gene_name": ["G"], "uniprot": ["P1"],
        "is_double_mutant": ["False"],
    }).to_csv(mave, sep="\t", index=False)
    _run(["--seq_table", str(seq), "--mavedb", str(mave), "--outdir", str(tmp_path)])
    df = pd.read_csv(tmp_path / "mavedb.tsv", sep="\t", dtype=str)
    assert df.iloc[0]["Protein_ID"] == "G-201"


# ---------------------------------------------------------------------------
# uniprot mode (new behaviour)
# ---------------------------------------------------------------------------

def _seq_table_uniprot(path):
    """Canonical isoform RAF1-201 (P04049) and alternative RAF1-205 (P04049-2).

    Canonical seq:  M A R G L D K (pos 1..7) ; WT at pos 5 = L
    Context window (pos4,5,6) = G L D = 'GLD'.
    Alt seq has a 2-residue N-terminal extension: X X M A R G L D K
       → 'GLD' appears at alt pos 6,7,8, so the WT residue (canonical pos 5, L)
         lands at alt pos 7 → homology_similarity, Protein_position=7.
    no-match isoform RAF1-301 (P04049) sequence has none of the context.
    """
    canon = "MARGLDK"            # pos5 = L, context GLD
    alt = "XXMARGLDK"            # GLD at alt pos 6-8 → pos5 L lands at alt pos 7
    nomatch = "QQQQQQQ"
    pd.DataFrame([
        {"Protein_ID": "RAF1-201", "Entry_Isoform": "P04049", "Gene": "RAF1",
         "Sequence": canon, "main_isoform": "yes"},
        {"Protein_ID": "RAF1-205", "Entry_Isoform": "P04049-2", "Gene": "RAF1",
         "Sequence": alt, "main_isoform": "no"},
        {"Protein_ID": "RAF1-301", "Entry_Isoform": "P04049-3", "Gene": "RAF1",
         "Sequence": nomatch, "main_isoform": "no"},
    ]).to_csv(path, sep="\t", index=False)


def _raw_table(path):
    """One raw MaveDB row: P04049 p.Leu5Pro at protein_start 5."""
    pd.DataFrame({
        "uniprot": ["P04049"],
        "gene_name": ["RAF1"],
        "urn": ["urn:mavedb:001"],
        "mavedb_id": ["urn:mavedb:001#0"],
        "prot_expr": ["p.Leu5Pro"],
        "protein_start": ["5"],
        "score": ["2.5"],
        "is_double_mutant": ["False"],
    }).to_csv(path, sep="\t", index=False)


def test_uniprot_mode_direct_and_homology(tmp_path):
    seq = tmp_path / "seq.tsv"
    raw = tmp_path / "mavedb_raw.tsv"
    _seq_table_uniprot(seq)
    _raw_table(raw)

    _run(["--mapping_mode", "uniprot", "--seq_table", str(seq),
          "--mavedb_raw", str(raw), "--outdir", str(tmp_path)])

    df = pd.read_csv(tmp_path / "mavedb.tsv", sep="\t", dtype=str)
    assert list(df.columns) == OUT_COLS

    by_pid = {r["Protein_ID"]: r for _, r in df.iterrows()}

    # canonical: WT L is already at pos 5 → direct
    assert "RAF1-201" in by_pid
    assert by_pid["RAF1-201"]["mapping_type"] == "direct"
    assert by_pid["RAF1-201"]["Protein_position"] == "5"

    # alt: context window R G L shifts L to alt pos 7 → homology_similarity
    assert "RAF1-205" in by_pid
    assert by_pid["RAF1-205"]["mapping_type"] == "homology_similarity"
    assert by_pid["RAF1-205"]["Protein_position"] == "7"

    # no-match isoform must be skipped entirely
    assert "RAF1-301" not in by_pid

    # provenance columns preserved
    assert by_pid["RAF1-201"]["uniprot"] == "P04049"
    assert by_pid["RAF1-201"]["urn"] == "urn:mavedb:001"
    assert by_pid["RAF1-201"]["score"] == "2.5"
    assert by_pid["RAF1-201"]["is_double_mutant"] == "False"


def test_uniprot_mode_empty_raw(tmp_path):
    seq = tmp_path / "seq.tsv"
    _seq_table_uniprot(seq)
    raw = tmp_path / "mavedb_raw.tsv"
    pd.DataFrame(columns=["uniprot", "gene_name", "urn", "mavedb_id",
                          "prot_expr", "protein_start", "score",
                          "is_double_mutant"]).to_csv(raw, sep="\t", index=False)
    _run(["--mapping_mode", "uniprot", "--seq_table", str(seq),
          "--mavedb_raw", str(raw), "--outdir", str(tmp_path)])
    df = pd.read_csv(tmp_path / "mavedb.tsv", sep="\t", dtype=str)
    assert list(df.columns) == OUT_COLS
    assert df.empty


def test_uniprot_mode_unknown_accession_skipped(tmp_path):
    """A raw row whose UniProt accession has no matching gene/isoform yields nothing."""
    seq = tmp_path / "seq.tsv"
    _seq_table_uniprot(seq)
    raw = tmp_path / "mavedb_raw.tsv"
    pd.DataFrame({
        "uniprot": ["Q99999"], "gene_name": ["XYZ"], "urn": ["urn:z"],
        "mavedb_id": ["urn:z#0"], "prot_expr": ["p.Ala2Gly"],
        "protein_start": ["2"], "score": ["1"], "is_double_mutant": ["False"],
    }).to_csv(raw, sep="\t", index=False)
    _run(["--mapping_mode", "uniprot", "--seq_table", str(seq),
          "--mavedb_raw", str(raw), "--outdir", str(tmp_path)])
    df = pd.read_csv(tmp_path / "mavedb.tsv", sep="\t", dtype=str)
    assert df.empty


def test_uniprot_mode_wt_mismatch_direct_skips_but_context_may_recover(tmp_path):
    """If the WT aa from prot_expr doesn't match the canonical residue at pos and
    there's no recoverable context window, that isoform is skipped."""
    seq = tmp_path / "seq.tsv"
    pd.DataFrame([
        {"Protein_ID": "RAF1-201", "Entry_Isoform": "P04049", "Gene": "RAF1",
         "Sequence": "MARGLDK", "main_isoform": "yes"},
    ]).to_csv(seq, sep="\t", index=False)
    raw = tmp_path / "mavedb_raw.tsv"
    # claim WT = Trp at pos 5 but canonical pos 5 is L → context 'GWD' not present
    pd.DataFrame({
        "uniprot": ["P04049"], "gene_name": ["RAF1"], "urn": ["urn:w"],
        "mavedb_id": ["urn:w#0"], "prot_expr": ["p.Trp5Pro"],
        "protein_start": ["5"], "score": ["1"], "is_double_mutant": ["False"],
    }).to_csv(raw, sep="\t", index=False)
    _run(["--mapping_mode", "uniprot", "--seq_table", str(seq),
          "--mavedb_raw", str(raw), "--outdir", str(tmp_path)])
    df = pd.read_csv(tmp_path / "mavedb.tsv", sep="\t", dtype=str)
    assert df.empty


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
