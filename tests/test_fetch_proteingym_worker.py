#!/usr/bin/env python3
"""Tests for bin/fetch_proteingym_worker.py — exercise the pure parse/join
helpers with tiny inline fixtures. No network access.
"""

import importlib.util
import io
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "fetch_proteingym_worker.py"

spec = importlib.util.spec_from_file_location("fetch_proteingym_worker", BIN)
fpg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fpg)


# ---------------------------------------------------------------------------
# position parsing from the `mutant` HGVS-lite string
# ---------------------------------------------------------------------------
class TestPosFromMutant:
    def test_simple(self):
        assert fpg.pos_from_mutant("G145R") == (145, "G")

    def test_multi_digit(self):
        assert fpg.pos_from_mutant("M1234V") == (1234, "M")

    def test_first_of_multi_mutant(self):
        # multi-mutant separated by ':' — first single is used
        assert fpg.pos_from_mutant("G145R:A200T") == (145, "G")

    def test_garbage(self):
        assert fpg.pos_from_mutant("") == (None, "")
        assert fpg.pos_from_mutant("wt") == (None, "")
        assert fpg.pos_from_mutant(None) == (None, "")


# ---------------------------------------------------------------------------
# UniProt accession extraction from the reference's UniProt_ID field
# ---------------------------------------------------------------------------
class TestAccession:
    def test_strip_organism(self):
        assert fpg.uniprot_accession("A0A140D2T1_ZIKV") == "A0A140D2T1"

    def test_plain_accession(self):
        assert fpg.uniprot_accession("P04637") == "P04637"

    def test_blank(self):
        assert fpg.uniprot_accession("") == ""
        assert fpg.uniprot_accession(None) == ""


# ---------------------------------------------------------------------------
# reference parsing: DMS_id → (uniprot, gene, filename)
# ---------------------------------------------------------------------------
REF_CSV = (
    "DMS_index,DMS_id,DMS_filename,UniProt_ID,molecule_name\n"
    "1,ASSAY_A,ASSAY_A.csv,P04637_HUMAN,p53\n"
    "2,ASSAY_B,ASSAY_B.csv,Q9Y6K9,IKBKG\n"
)


class TestParseReference:
    def test_parse(self):
        ref = fpg.parse_reference(io.StringIO(REF_CSV))
        assert ref["ASSAY_A"]["uniprot"] == "P04637"
        assert ref["ASSAY_A"]["gene_name"] == "p53"
        assert ref["ASSAY_A"]["filename"] == "ASSAY_A.csv"
        assert ref["ASSAY_B"]["uniprot"] == "Q9Y6K9"


# ---------------------------------------------------------------------------
# per-assay CSV parsing → long rows
# ---------------------------------------------------------------------------
ASSAY_CSV = (
    "mutant,mutated_sequence,DMS_score,DMS_score_bin\n"
    "G145R,MAAA,1.23,1\n"
    "A200T,MAAB,-0.5,0\n"
    "wt,MAAA,0.0,1\n"          # unparseable mutant → skipped
)


class TestParseAssay:
    def test_rows(self):
        rows = fpg.parse_assay(
            io.StringIO(ASSAY_CSV), "ASSAY_A", "P04637", "p53")
        assert len(rows) == 2          # wt row skipped
        r0 = rows[0]
        assert r0["uniprot"] == "P04637"
        assert r0["gene_name"] == "p53"
        assert r0["DMS_id"] == "ASSAY_A"
        assert r0["protein_variant"] == "G145R"
        assert r0["pos"] == "145"
        assert r0["DMS_score"] == "1.23"
        assert r0["DMS_score_bin"] == "1"

    def test_limit(self):
        rows = fpg.parse_assay(
            io.StringIO(ASSAY_CSV), "ASSAY_A", "P04637", "p53", limit=1)
        assert len(rows) == 1

    def test_out_cols_present(self):
        rows = fpg.parse_assay(
            io.StringIO(ASSAY_CSV), "ASSAY_A", "P04637", "p53")
        for r in rows:
            assert set(r.keys()) == set(fpg.OUT_COLS)
