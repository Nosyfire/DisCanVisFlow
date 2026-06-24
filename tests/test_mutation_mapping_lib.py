"""Tests for mutation_mapping_lib.py"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from mutation_mapping_lib import (
    filter_hypermutated_samples,
    normalize_tcga_sample,
    parse_hgvsp_ref_pos,
    validate_hgvsp_aa,
)


class TestHgvspValidation:
    def test_parse_missense(self):
        ref, pos = parse_hgvsp_ref_pos("p.V403A")
        assert ref == "V"
        assert pos == 403

    def test_validate_match(self):
        assert validate_hgvsp_aa("p.V403A", "V", 403) is True

    def test_validate_mismatch(self):
        assert validate_hgvsp_aa("p.V403A", "L", 403) is False

    def test_empty_hgvs_passes(self):
        assert validate_hgvsp_aa("", "V", 403) is True


class TestTcgaSample:
    def test_truncates_barcode(self):
        s = "TCGA-AB-1234-01A-11D-1234-56"
        assert normalize_tcga_sample(s, "TCGA") == "TCGA-AB-1234"


class TestHypermutation:
    def test_drops_hot_sample(self):
        variants = [{"sample": "S1"}] * 1600 + [{"sample": "S2"}]
        out = filter_hypermutated_samples(variants, threshold=1500)
        assert all(v["sample"] == "S2" for v in out)
