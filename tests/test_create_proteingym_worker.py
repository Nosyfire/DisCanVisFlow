#!/usr/bin/env python3
"""Tests for bin/create_proteingym_worker.py.

Covers:
  - premapped mode (existing behaviour, Protein_ID-keyed filter)
  - uniprot mapping mode with isoform fan-out (direct vs homology_similarity)
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_proteingym_worker.py"

OUT_COLS = ["Protein_ID", "Protein_position", "protein_variant", "DMS_score",
            "DMS_score_bin", "DMS_id", "uniprot_id", "mapping_type"]


def _run(args):
    res = subprocess.run([sys.executable, str(BIN)] + args,
                         capture_output=True, text=True)
    return res


# ---------------------------------------------------------------------------
# Fixtures: a two-isoform gene.
#   - GENEX-201 / P00001          (main isoform)  seq:  MKAAGAAR  (len 8)
#   - GENEX-202 / P00001-2        (alt isoform)   seq:  MKAAXGAAR (insertion at pos5)
# pos 7 in the main is 'A'; the 3-aa context 'GAA' (pos6-8) re-locates to the alt
# isoform shifted by +1.
# ---------------------------------------------------------------------------
SEQ_MAIN = "MKAAGAAR"
SEQ_ALT  = "MKAAXGAAR"


def _write_seq_table(p: Path):
    cols = ["Entry_Name", "Gene_Gencode", "Gene", "Protein_ID",
            "Entry_Isoform", "main_isoform", "Sequence"]
    rows = [
        ["GENEX_HUMAN", "GENEX", "GENEX", "GENEX-201", "P00001",   "yes", SEQ_MAIN],
        ["GENEX_HUMAN", "GENEX", "GENEX", "GENEX-202", "P00001-2", "no",  SEQ_ALT],
    ]
    pd.DataFrame(rows, columns=cols).to_csv(p, sep="\t", index=False)


# ===========================================================================
# premapped mode (unchanged behaviour)
# ===========================================================================
class TestPremapped:
    def test_filters_to_run_proteins(self, tmp_path):
        seq = tmp_path / "seq.tsv"
        _write_seq_table(seq)
        pg = tmp_path / "premapped.tsv"
        pd.DataFrame([
            ["GENEX-201", "5", "G5A", "1.1", "1", "ASSAY_A", "P00001"],
            ["OTHER-201", "9", "K9R", "0.2", "0", "ASSAY_B", "P99999"],
        ], columns=["Protein_ID", "pos", "protein_variant", "DMS_score",
                    "DMS_score_bin", "DMS_id", "uniprot_id"]
        ).to_csv(pg, sep="\t", index=False)

        res = _run(["--seq_table", str(seq), "--proteingym", str(pg),
                    "--outdir", str(tmp_path)])
        assert res.returncode == 0, res.stderr
        out = pd.read_csv(tmp_path / "proteingym.tsv", sep="\t", dtype=str)
        assert list(out.columns) == OUT_COLS
        assert set(out["Protein_ID"]) == {"GENEX-201"}
        assert (out["mapping_type"] == "direct").all()

    def test_default_mode_is_premapped(self, tmp_path):
        seq = tmp_path / "seq.tsv"
        _write_seq_table(seq)
        pg = tmp_path / "premapped.tsv"
        pd.DataFrame(columns=["Protein_ID", "pos", "protein_variant", "DMS_score",
                              "DMS_score_bin", "DMS_id", "uniprot_id"]
                     ).to_csv(pg, sep="\t", index=False)
        res = _run(["--seq_table", str(seq), "--proteingym", str(pg),
                    "--outdir", str(tmp_path)])
        assert res.returncode == 0, res.stderr
        out = pd.read_csv(tmp_path / "proteingym.tsv", sep="\t", dtype=str)
        assert list(out.columns) == OUT_COLS


# ===========================================================================
# uniprot mode — isoform fan-out
# ===========================================================================
def _write_raw(p: Path, rows):
    cols = ["uniprot", "gene_name", "DMS_id", "protein_variant", "pos",
            "DMS_score", "DMS_score_bin"]
    pd.DataFrame(rows, columns=cols).to_csv(p, sep="\t", index=False)


class TestUniprotMode:
    def test_direct_match_on_native_isoform(self, tmp_path):
        seq = tmp_path / "seq.tsv"
        _write_seq_table(seq)
        raw = tmp_path / "raw.tsv"
        # G5A : main isoform residue 5 is 'G' (MKAAG...) → direct on GENEX-201
        _write_raw(raw, [["P00001", "GENEX", "ASSAY_A", "G5A", "5", "1.1", "1"]])

        res = _run(["--seq_table", str(seq),
                    "--mapping_mode", "uniprot",
                    "--proteingym_raw", str(raw),
                    "--outdir", str(tmp_path)])
        assert res.returncode == 0, res.stderr
        out = pd.read_csv(tmp_path / "proteingym.tsv", sep="\t", dtype=str)
        assert list(out.columns) == OUT_COLS

        direct = out[out["Protein_ID"] == "GENEX-201"]
        assert len(direct) == 1
        assert direct.iloc[0]["mapping_type"] == "direct"
        assert direct.iloc[0]["Protein_position"] == "5"

    def test_homology_fanout_to_alt_isoform(self, tmp_path):
        seq = tmp_path / "seq.tsv"
        _write_seq_table(seq)
        raw = tmp_path / "raw.tsv"
        # G6A on main: main seq MKAAGAAR, residue 6 = 'A'. Wait: build a variant at
        # pos 7 where main residue is 'A' (MKAAGAAR -> 1M2K3A4A5G6A7A8R), context
        # 'GAA' (pos5-7) relocates to alt seq MKAAXGAAR (pos6-8) → shifted Protein_position 8.
        _write_raw(raw, [["P00001", "GENEX", "ASSAY_A", "A6V", "6", "0.5", "0"]])

        res = _run(["--seq_table", str(seq),
                    "--mapping_mode", "uniprot",
                    "--proteingym_raw", str(raw),
                    "--outdir", str(tmp_path)])
        assert res.returncode == 0, res.stderr
        out = pd.read_csv(tmp_path / "proteingym.tsv", sep="\t", dtype=str)

        # direct hit on the native isoform (residue 6 of MKAAGAAR is 'A')
        direct = out[out["Protein_ID"] == "GENEX-201"]
        assert len(direct) == 1
        assert direct.iloc[0]["mapping_type"] == "direct"

        # homology fan-out onto the alt isoform with shifted position
        homo = out[out["Protein_ID"] == "GENEX-202"]
        assert len(homo) == 1
        assert homo.iloc[0]["mapping_type"] == "homology_similarity"
        # context around pos6 in main = pos5..7 'GAA' → in alt MKAAXGAAR at index5..7,
        # remapped variant position = 7 (1-based)
        assert homo.iloc[0]["Protein_position"] == "7"

    def test_wt_mismatch_drops_row(self, tmp_path):
        seq = tmp_path / "seq.tsv"
        _write_seq_table(seq)
        raw = tmp_path / "raw.tsv"
        # WT aa says 'W' at pos 5 but main residue is 'G' → no direct, no context.
        _write_raw(raw, [["P00001", "GENEX", "ASSAY_A", "W5A", "5", "0.1", "0"]])

        res = _run(["--seq_table", str(seq),
                    "--mapping_mode", "uniprot",
                    "--proteingym_raw", str(raw),
                    "--outdir", str(tmp_path)])
        assert res.returncode == 0, res.stderr
        out = pd.read_csv(tmp_path / "proteingym.tsv", sep="\t", dtype=str)
        # the native isoform should not get a direct hit (WT mismatch)
        assert "GENEX-201" not in set(out["Protein_ID"])

    def test_unknown_uniprot_yields_empty(self, tmp_path):
        seq = tmp_path / "seq.tsv"
        _write_seq_table(seq)
        raw = tmp_path / "raw.tsv"
        _write_raw(raw, [["P12345", "NOGENE", "ASSAY_X", "G5A", "5", "1.1", "1"]])
        res = _run(["--seq_table", str(seq),
                    "--mapping_mode", "uniprot",
                    "--proteingym_raw", str(raw),
                    "--outdir", str(tmp_path)])
        assert res.returncode == 0, res.stderr
        out = pd.read_csv(tmp_path / "proteingym.tsv", sep="\t", dtype=str)
        assert len(out) == 0
        assert list(out.columns) == OUT_COLS
