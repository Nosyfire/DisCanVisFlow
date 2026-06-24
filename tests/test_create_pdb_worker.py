"""Tests for bin/create_pdb_worker.py — structure→transcript mapping + missing residues."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from create_pdb_worker import (  # noqa: E402
    missing_unp_ranges,
    best_window,
    build_gene_lookup,
)


class TestMissingResidues:
    def test_fully_observed_has_no_missing(self):
        seg = {"unp_start": 10, "unp_end": 20, "label_start": 1, "label_end": 11}
        # observed label 1..11 → unp 10..20 fully observed
        assert missing_unp_ranges(seg, [(1, 11)]) == []

    def test_internal_gap_is_missing(self):
        # unp 10..20, offset = 10-1 = 9. Observed label 1..3 (unp 10..12) and
        # label 8..11 (unp 17..20) → unp 13..16 missing.
        seg = {"unp_start": 10, "unp_end": 20, "label_start": 1, "label_end": 11}
        assert missing_unp_ranges(seg, [(1, 3), (8, 11)]) == [(13, 16)]

    def test_terminal_missing(self):
        seg = {"unp_start": 1, "unp_end": 10, "label_start": 1, "label_end": 10}
        # only label 3..8 observed → unp 1..2 and 9..10 missing
        assert missing_unp_ranges(seg, [(3, 8)]) == [(1, 2), (9, 10)]

    def test_no_coverage_all_missing(self):
        seg = {"unp_start": 5, "unp_end": 8, "label_start": 1, "label_end": 4}
        assert missing_unp_ranges(seg, []) == [(5, 8)]


class TestBestWindow:
    def test_exact(self):
        idx, ident = best_window("MAEAK", "XXMAEAKXX", 0.9)
        assert idx == 2 and ident == 1.0

    def test_near_identical(self):
        idx, ident = best_window("MAEAKLLPKL", "MAEAKLLPRL", 0.9)  # 9/10
        assert idx == 0 and abs(ident - 0.9) < 1e-9

    def test_below_threshold(self):
        idx, ident = best_window("MAEAKLLPKL", "ZZZZZZZZZZ", 0.9)
        assert idx == -1


class TestGeneLookup:
    def test_groups_by_gene(self):
        df = pd.DataFrame({
            "Entry_Isoform": ["P1", "P1", "P2"],
            "Protein_ID": ["G-201", "G-202", "G-301"],
            "Sequence": ["MAEAK", "MAEAK", "MMMMM"],
            "Gene_Gencode": ["G", "G", "H"],
            "main_isoform": ["yes", "no", "yes"],
        })
        gene_to_rows, acc_to_seq, acc_main, _ = build_gene_lookup(df)
        assert len(gene_to_rows["G"]) == 2
        assert acc_main["P1"] is True
        assert acc_to_seq["P2"] == "MMMMM"
