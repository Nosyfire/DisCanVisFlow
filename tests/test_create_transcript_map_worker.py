"""
Tests for create_transcript_map_worker.py (Module 5e)

Covers:
- Isoform bounds check (same UniProt acc, different Protein_ID)
- Homology transfer (different UniProt acc, same gene, matching sequence)
- Out-of-bounds annotations are dropped
- Pfam is accepted via --pfam argument
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).parent.parent / "bin" / "create_transcript_map_worker.py"
PYTHON = sys.executable

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from create_transcript_map_worker import (
    build_lookup,
    map_annotations_to_transcripts,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_loc(rows):
    """rows: list of (acc, pid, seq, gene, is_main)"""
    data = {
        "Entry_Isoform": [r[0] for r in rows],
        "Protein_ID":    [r[1] for r in rows],
        "Sequence":      [r[2] for r in rows],
        "Gene_Gencode":  [r[3] for r in rows],
        "main_isoform":  [r[4] for r in rows],
    }
    return pd.DataFrame(data)


def _make_annot(acc, start, end):
    return pd.DataFrame([{"Entry_Isoform": acc, "Name": "test", "Start": start, "End": end}])


# ---------------------------------------------------------------------------
# Bounds check for same-acc different-pid isoforms
# ---------------------------------------------------------------------------

class TestBoundsCheck:
    def test_in_bounds_copied(self):
        seq_long  = "MAEAKLLPKL"  # len 10
        seq_short = "MAEAKLLL"    # len 8
        loc_df = _make_loc([
            ("P001", "G-201", seq_long,  "G", "yes"),
            ("P001", "G-202", seq_short, "G", "no"),
        ])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 5)  # End=5 ≤ len(seq_short)=8

        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid
        )
        pids = result["Protein_ID"].tolist()
        assert "G-201" in pids
        assert "G-202" in pids

    def test_out_of_bounds_dropped(self):
        seq_long  = "MAEAKLLPKL"  # len 10
        seq_short = "MAE"          # len 3
        loc_df = _make_loc([
            ("P001", "G-201", seq_long,  "G", "yes"),
            ("P001", "G-202", seq_short, "G", "no"),
        ])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 8)  # End=8 > len(seq_short)=3

        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid
        )
        pids = result["Protein_ID"].tolist()
        assert "G-201" in pids   # main isoform is kept
        assert "G-202" not in pids  # short isoform is dropped

    def test_same_pid_always_copied(self):
        """When Protein_ID matches source, copy directly without bounds check."""
        seq = "MAEAKLLPKL"
        loc_df = _make_loc([("P001", "G-201", seq, "G", "yes")])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 10)

        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid
        )
        assert len(result) == 1
        assert result.iloc[0]["homology_transfer"] is False or result.iloc[0]["homology_transfer"] == False


# ---------------------------------------------------------------------------
# Homology transfer tests
# ---------------------------------------------------------------------------

class TestHomologyTransfer:
    def test_exact_region_match_transfers(self):
        # seq1 has MAEAK at start; seq2 also has MAEAK but after different prefix
        seq1 = "MAEAKLLPKL"
        seq2 = "XXXMAEAKXXX"  # MAEAK present at offset 3
        loc_df = _make_loc([
            ("P001", "G-201", seq1, "G", "yes"),
            ("P002", "G-202", seq2, "G", "no"),
        ])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 5)  # region = seq1[0:5] = "MAEAK"

        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid
        )
        pids = result["Protein_ID"].tolist()
        assert "G-201" in pids
        assert "G-202" in pids
        row_202 = result[result["Protein_ID"] == "G-202"].iloc[0]
        assert row_202["homology_transfer"] == True
        assert row_202["Start"] == 4   # 0-indexed find=3 → 1-based=4
        assert row_202["End"]   == 8   # 4 + len("MAEAK") - 1

    def test_no_match_not_transferred(self):
        seq1 = "MAEAKLLPKL"
        seq2 = "ZZZZZZZZZZZ"  # no match
        loc_df = _make_loc([
            ("P001", "G-201", seq1, "G", "yes"),
            ("P002", "G-202", seq2, "G", "no"),
        ])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 5)

        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid
        )
        pids = result["Protein_ID"].tolist()
        assert "G-202" not in pids

    def test_exact_transfer_has_identity_1(self):
        seq1 = "MAEAKLLPKL"
        seq2 = "XXXMAEAKXXX"
        loc_df = _make_loc([
            ("P001", "G-201", seq1, "G", "yes"),
            ("P002", "G-202", seq2, "G", "no"),
        ])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 5)
        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid
        )
        row = result[result["Protein_ID"] == "G-202"].iloc[0]
        assert row["mapping_type"] == "homology_similarity"
        assert float(row["homology_identity"]) == 1.0

    def test_near_identical_region_transfers_above_threshold(self):
        # 10-residue region with a single substitution in the target = 90% identity
        region = "MAEAKLLPKL"
        seq1 = region + "GG"
        seq2 = "MAEAKLLPRL"  # position 9 K->R  → 9/10 = 0.9 identity
        loc_df = _make_loc([
            ("P001", "G-201", seq1, "G", "yes"),
            ("P002", "G-202", seq2, "G", "no"),
        ])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 10)
        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid,
            min_identity=0.9,
        )
        row = result[result["Protein_ID"] == "G-202"]
        assert len(row) == 1
        assert row.iloc[0]["homology_transfer"] == True
        assert abs(float(row.iloc[0]["homology_identity"]) - 0.9) < 1e-6

    def test_below_threshold_not_transferred(self):
        region = "MAEAKLLPKL"            # 10 aa
        seq1 = region + "GG"
        seq2 = "MAEZZZZZKL"              # only 5/10 identical = 0.5
        loc_df = _make_loc([
            ("P001", "G-201", seq1, "G", "yes"),
            ("P002", "G-202", seq2, "G", "no"),
        ])
        acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid = build_lookup(loc_df)
        annot_df = _make_annot("P001", 1, 10)
        result = map_annotations_to_transcripts(
            annot_df, "Entry_Isoform", acc_to_seq, acc_to_gene, gene_to_rows, acc_to_pid,
            min_identity=0.9,
        )
        assert "G-202" not in result["Protein_ID"].tolist()


# ---------------------------------------------------------------------------
# CLI test with pfam
# ---------------------------------------------------------------------------

class TestCLIWithPfam:
    def _minimal_files(self, tmp_path, with_pfam=True):
        loc_content = (
            "Entry_Isoform\tProtein_ID\tSequence\tGene_Gencode\tmain_isoform\n"
            "P001\tG-201\tMEAKLLP\tGENE\tyes\n"
        )
        (tmp_path / "loc.tsv").write_text(loc_content)

        for name in ["elm", "dibs", "mfib", "phasepro", "roi", "bind", "ptm", "pfam",
                     "disorder", "disorder_pos"]:
            (tmp_path / f"{name}.tsv").write_text("Entry_Isoform\tStart\tEnd\n")

        return tmp_path

    def test_cli_with_pfam_produces_output(self, tmp_path):
        self._minimal_files(tmp_path)
        result = subprocess.run(
            [PYTHON, str(BIN),
             "--loc_chrom",   str(tmp_path / "loc.tsv"),
             "--elm",         str(tmp_path / "elm.tsv"),
             "--dibs",        str(tmp_path / "dibs.tsv"),
             "--mfib",        str(tmp_path / "mfib.tsv"),
             "--phasepro",    str(tmp_path / "phasepro.tsv"),
             "--uniprot_roi", str(tmp_path / "roi.tsv"),
             "--uniprot_bind",str(tmp_path / "bind.tsv"),
             "--ptm",         str(tmp_path / "ptm.tsv"),
             "--pfam",        str(tmp_path / "pfam.tsv"),
             "--disorder",    str(tmp_path / "disorder.tsv"),
             "--disorder_pos",str(tmp_path / "disorder_pos.tsv"),
             "--output_dir",  str(tmp_path)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "pfam_domains.tsv").exists()
        assert (tmp_path / "transcript_map_stats.tsv").exists()

    def test_stats_contains_pfam(self, tmp_path):
        self._minimal_files(tmp_path)
        subprocess.run(
            [PYTHON, str(BIN),
             "--loc_chrom",   str(tmp_path / "loc.tsv"),
             "--elm",         str(tmp_path / "elm.tsv"),
             "--dibs",        str(tmp_path / "dibs.tsv"),
             "--mfib",        str(tmp_path / "mfib.tsv"),
             "--phasepro",    str(tmp_path / "phasepro.tsv"),
             "--uniprot_roi", str(tmp_path / "roi.tsv"),
             "--uniprot_bind",str(tmp_path / "bind.tsv"),
             "--ptm",         str(tmp_path / "ptm.tsv"),
             "--pfam",        str(tmp_path / "pfam.tsv"),
             "--disorder",    str(tmp_path / "disorder.tsv"),
             "--disorder_pos",str(tmp_path / "disorder_pos.tsv"),
             "--output_dir",  str(tmp_path)],
            capture_output=True, text=True,
        )
        stats = pd.read_csv(tmp_path / "transcript_map_stats.tsv", sep="\t")
        assert "pfam" in stats.columns
