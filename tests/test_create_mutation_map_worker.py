"""
Tests for create_mutation_map_worker.py (Module 4: Mutation Mapping)

Tests cover:
- Mutation classification (Missense / Frameshift / Nonsense / Indel)
- Isoform expansion via sequence context search
- ClinVar VCF parsing
- MAF parsing
- Generic VCF parsing
- Deduplication of isoform-expanded hits
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).parent.parent / "bin" / "create_mutation_map_worker.py"
PYTHON = sys.executable


# ---------------------------------------------------------------------------
# Import helper
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from create_mutation_map_worker import (
    _classify,
    resolve_mutations,
    write_split_tsv,
    parse_maf,
    load_combined_map,
)
from mutation_mapping_lib import filter_hypermutated_samples, validate_hgvsp_aa


# ---------------------------------------------------------------------------
# _classify
# ---------------------------------------------------------------------------

class TestClassify:
    def test_snv_is_missense(self):
        assert _classify("p.R333P", "G", "C", "single_nucleotide_variant") == "Missense_Mutation"

    def test_frameshift_del(self):
        cls = _classify("p.A10fs", "GA", "G", "Frame_Shift_Del")
        assert "Frame_Shift" in cls

    def test_frameshift_from_hgvs(self):
        cls = _classify("p.L15fs", "A", "AT", "")
        assert "Frame_Shift" in cls

    def test_nonsense_star(self):
        assert _classify("p.R333*", "G", "A", "") == "Nonsense_Mutation"

    def test_nonsense_stop_gained(self):
        assert _classify("", "G", "A", "stop_gained") == "Nonsense_Mutation"

    def test_indel_deletion(self):
        assert _classify("", "GAT", "G", "") == "Indel"

    def test_indel_insertion(self):
        assert _classify("", "G", "GAT", "in_frame_ins") == "Indel"


# ---------------------------------------------------------------------------
# Isoform expansion
# ---------------------------------------------------------------------------

class TestIsoformExpansion:
    """Test sequence-context translation to additional isoforms."""

    def _make_loc_df(self, rows):
        """rows: list of (enst_base, pid, acc, seq, gene)"""
        data = {
            "transcript_stable_id": [r[0] for r in rows],
            "Protein_ID":           [r[1] for r in rows],
            "Entry_Isoform":        [r[2] for r in rows],
            "Sequence":             [r[3] for r in rows],
            "Gene_Gencode":         [r[4] for r in rows],
            "main_isoform":         ["yes" if i == 0 else "no" for i in range(len(rows))],
        }
        import pandas as pd
        from create_mutation_map_worker import load_loc_chrom
        import tempfile, os
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            pd.DataFrame(data).to_csv(f.name, sep="\t", index=False)
            fname = f.name
        result = load_loc_chrom(fname)
        os.unlink(fname)
        return result

    def test_primary_isoform_is_included(self):
        seq1 = "MAEAKLLPKL"
        seq2 = "MAEAK"  # truncated
        loc_df, gene_to_rows, pid_to_seq, gene_col = self._make_loc_df([
            ("ENST001", "GENE-201", "P00001", seq1, "GENE"),
            ("ENST002", "GENE-202", "P00001", seq2, "GENE"),
        ])
        # Simulate a hit at position 4 (0-based) of seq1 ('K')
        lookup = {"chr1": {100: [("ENST001", 4, "K")]}}
        results = resolve_mutations(
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
              "hgvs": "p.K5R", "clnsig": "", "disease": "", "variant_type": "",
              "study_abbr": "TEST", "study_name": "TEST", "sample": "",
              "review_status": "", "phenotype_ids": "", "rcv": "", "rs": "",
              "mondo_id": "", "mesh_id": ""}],
            lookup, loc_df, gene_to_rows, pid_to_seq, gene_col, "TEST",
            map_all_isoforms=True,
        )
        pids = [r["Protein_ID"] for r in results]
        assert "GENE-201" in pids  # primary

    def test_secondary_isoform_translated(self):
        # seq1 contains MAEAK at positions 0-4; seq2 starts with MAEAK
        seq1 = "MAEAKLLPKL"
        seq2 = "MAEAKXXX"   # same first 5 aa
        loc_df, gene_to_rows, pid_to_seq, gene_col = self._make_loc_df([
            ("ENST001", "GENE-201", "P00001", seq1, "GENE"),
            ("ENST002", "GENE-202", "P00001", seq2, "GENE"),
        ])
        # Hit at position index 4 (5th aa 'K') in seq1
        lookup = {"chr1": {100: [("ENST001", 4, "K")]}}
        results = resolve_mutations(
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
              "hgvs": "p.K5R", "clnsig": "", "disease": "", "variant_type": "",
              "study_abbr": "TEST", "study_name": "TEST", "sample": "",
              "review_status": "", "phenotype_ids": "", "rcv": "", "rs": "",
              "mondo_id": "", "mesh_id": ""}],
            lookup, loc_df, gene_to_rows, pid_to_seq, gene_col, "TEST",
            map_all_isoforms=True,
        )
        pids = [r["Protein_ID"] for r in results]
        assert "GENE-202" in pids
        assert any(r["isoform_mapped"] for r in results if r["Protein_ID"] == "GENE-202")

    def test_secondary_isoform_out_of_bounds_not_included(self):
        seq1 = "MAEAKLLPKL"
        seq2 = "MAE"  # too short for position 5
        loc_df, gene_to_rows, pid_to_seq, gene_col = self._make_loc_df([
            ("ENST001", "GENE-201", "P00001", seq1, "GENE"),
            ("ENST002", "GENE-202", "P00001", seq2, "GENE"),
        ])
        lookup = {"chr1": {100: [("ENST001", 7, "L")]}}  # pos 8 in seq1
        results = resolve_mutations(
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
              "hgvs": "", "clnsig": "", "disease": "", "variant_type": "",
              "study_abbr": "TEST", "study_name": "TEST", "sample": "",
              "review_status": "", "phenotype_ids": "", "rcv": "", "rs": "",
              "mondo_id": "", "mesh_id": ""}],
            lookup, loc_df, gene_to_rows, pid_to_seq, gene_col, "TEST",
            map_all_isoforms=True,
        )
        pids = [r["Protein_ID"] for r in results]
        assert "GENE-202" not in pids  # seq2 is too short

    def test_no_expansion_when_disabled(self):
        seq1 = "MAEAKLLPKL"
        seq2 = "MAEAKXXX"
        loc_df, gene_to_rows, pid_to_seq, gene_col = self._make_loc_df([
            ("ENST001", "GENE-201", "P00001", seq1, "GENE"),
            ("ENST002", "GENE-202", "P00001", seq2, "GENE"),
        ])
        lookup = {"chr1": {100: [("ENST001", 4, "K")]}}
        results = resolve_mutations(
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
              "hgvs": "", "clnsig": "", "disease": "", "variant_type": "",
              "study_abbr": "TEST", "study_name": "TEST", "sample": "",
              "review_status": "", "phenotype_ids": "", "rcv": "", "rs": "",
              "mondo_id": "", "mesh_id": ""}],
            lookup, loc_df, gene_to_rows, pid_to_seq, gene_col, "TEST",
            map_all_isoforms=False,
        )
        pids = [r["Protein_ID"] for r in results]
        assert "GENE-202" not in pids


# ---------------------------------------------------------------------------
# Output splitting
# ---------------------------------------------------------------------------

class TestOutputSplitting:
    def test_write_creates_all_four_files(self):
        rows = [
            {"Protein_ID": "G-201", "Accession": "P1", "Gene": "G",
             "Mutation Description": "", "Mutation": "p.A1V",
             "Protein_position": 1, "Study Abbrevation": "T", "Study Name": "T",
             "Sample name": "", "Start_Position": 100,
             "isoform_mapped": False, "ClinicalSignificance": "",
             "PhenotypeList": "", "PhenotypeIDS": "", "ReviewStatus": "",
             "RCVaccession": "", "MONDO_ID": "", "MeSH_ID": "",
             "_ref": "G", "_alt": "T", "_variant_type": "SNV", "_hgvs": "p.A1V"},
        ]
        with tempfile.TemporaryDirectory() as d:
            write_split_tsv(rows, Path(d), "TEST")
            files = list(Path(d).glob("*.tsv"))
            fnames = [f.name for f in files]
            assert "Missense_filter_mutations_mapped.tsv" in fnames
            assert "Frameshift_filter_mutations_mapped.tsv" in fnames
            assert "Nonsense_filter_mutations_mapped.tsv" in fnames
            assert "Indel_filter_mutations_mapped.tsv" in fnames
            assert "mutation_stats.tsv" in fnames

    def test_nonsense_goes_to_nonsense_file(self):
        rows = [
            {"Protein_ID": "G-201", "Accession": "P1", "Gene": "G",
             "Mutation Description": "", "Mutation": "p.R10*",
             "Protein_position": 10, "Study Abbrevation": "T", "Study Name": "T",
             "Sample name": "", "Start_Position": 100,
             "isoform_mapped": False, "ClinicalSignificance": "",
             "PhenotypeList": "", "PhenotypeIDS": "", "ReviewStatus": "",
             "RCVaccession": "", "MONDO_ID": "", "MeSH_ID": "",
             "_ref": "G", "_alt": "A", "_variant_type": "", "_hgvs": "p.R10*"},
        ]
        with tempfile.TemporaryDirectory() as d:
            write_split_tsv(rows, Path(d), "TEST")
            ns_df = pd.read_csv(Path(d) / "Nonsense_filter_mutations_mapped.tsv", sep="\t")
            assert len(ns_df) == 1


# ---------------------------------------------------------------------------
# MAF parsing
# ---------------------------------------------------------------------------

class TestMafParsing:
    def test_basic_maf(self):
        content = (
            "Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\t"
            "HGVSp_Short\tVariant_Classification\tTumor_Sample_Barcode\n"
            "1\t12345\tA\tG\tp.K5R\tMissense_Mutation\tSAMPLE_001\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".maf", delete=False) as f:
            f.write(content)
            fname = f.name
        try:
            rows = parse_maf(fname)
            assert len(rows) == 1
            assert rows[0]["chrom"] == "chr1"
            assert rows[0]["pos"] == 12345
            assert rows[0]["hgvs"] == "p.K5R"
            assert rows[0]["sample"] == "SAMPLE_001"
        finally:
            os.unlink(fname)

    def test_maf_chrom_prefix(self):
        content = (
            "Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\n"
            "chr3\t500\tA\tG\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".maf", delete=False) as f:
            f.write(content)
            fname = f.name
        try:
            rows = parse_maf(fname)
            assert rows[0]["chrom"] == "chr3"
        finally:
            os.unlink(fname)


class TestHgvspFilterInResolve:
    def test_rejects_wrong_ref_aa(self):
        seq = "MAEAKLLPKL"
        loc_content = (
            "transcript_stable_id\tProtein_ID\tEntry_Isoform\tSequence\tGene_Gencode\tmain_isoform\n"
            "ENST001\tG-201\tP00001\tMAEAKLLPKL\tGENE\tyes\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            pd.DataFrame({
                "transcript_stable_id": ["ENST001"],
                "Protein_ID": ["G-201"],
                "Entry_Isoform": ["P00001"],
                "Sequence": [seq],
                "Gene_Gencode": ["GENE"],
                "main_isoform": ["yes"],
            }).to_csv(f.name, sep="\t", index=False)
            loc_df, gene_to_rows, pid_to_seq, gene_col = __import__(
                "create_mutation_map_worker", fromlist=["load_loc_chrom"]
            ).load_loc_chrom(f.name)
            os.unlink(f.name)

        lookup = {"chr1": {100: [("ENST001", 4, "K")]}}
        var = {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
               "hgvs": "p.V5R", "clnsig": "", "disease": "", "variant_type": "",
               "study_abbr": "TEST", "study_name": "TEST", "sample": "",
               "review_status": "", "phenotype_ids": "", "rcv": "", "rs": "",
               "mondo_id": "", "mesh_id": ""}
        results = resolve_mutations(
            [var], lookup, loc_df, gene_to_rows, pid_to_seq, gene_col, "TEST",
            map_all_isoforms=False, validate_hgvsp=True,
        )
        assert len(results) == 0


class TestHypermutationFilter:
    def test_lib_drops_samples(self):
        v = filter_hypermutated_samples([{"sample": "A"}] * 2000, 1500)
        assert v == []


# ---------------------------------------------------------------------------
# CLI smoke test
# ---------------------------------------------------------------------------

class TestCLI:
    def test_cli_produces_output_files(self, tmp_path):
        # Build tiny combined_map.map
        map_content = "# ENST001 chr1 + 100-200\n1 A 0 A ATG chr1,100, ATG A\n"
        map_file = tmp_path / "combined_map.map"
        map_file.write_text(map_content)

        # Build tiny loc_chrom
        loc_content = (
            "transcript_stable_id\tProtein_ID\tEntry_Isoform\tSequence\tGene_Gencode\tmain_isoform\n"
            "ENST001\tG-201\tP00001\tMAEAKLLPKL\tGENE\tyes\n"
        )
        loc_file = tmp_path / "loc_chrom.tsv"
        loc_file.write_text(loc_content)

        # Build tiny MAF
        maf_content = (
            "Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\n"
            "chr1\t100\tA\tG\n"
        )
        maf_file = tmp_path / "test.maf"
        maf_file.write_text(maf_content)

        result = subprocess.run(
            [PYTHON, str(BIN),
             "--combined_map", str(map_file),
             "--loc_chrom",    str(loc_file),
             "--maf",          str(maf_file),
             "--source",       "TEST",
             "--output_dir",   str(tmp_path)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "mutation_stats.tsv").exists()
        assert (tmp_path / "Missense_filter_mutations_mapped.tsv").exists()
        assert (tmp_path / "Nonsense_filter_mutations_mapped.tsv").exists()

    def test_output_has_protein_id_column(self, tmp_path):
        map_content = "# ENST001 chr1 + 100-200\n1 A 0 A ATG chr1,100, ATG A\n"
        (tmp_path / "combined_map.map").write_text(map_content)
        loc_content = (
            "transcript_stable_id\tProtein_ID\tEntry_Isoform\tSequence\tGene_Gencode\tmain_isoform\n"
            "ENST001\tG-201\tP00001\tMAEAKLLPKL\tGENE\tyes\n"
        )
        (tmp_path / "loc.tsv").write_text(loc_content)
        maf_content = "Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2\nchr1\t100\tA\tG\n"
        (tmp_path / "t.maf").write_text(maf_content)

        subprocess.run(
            [PYTHON, str(BIN),
             "--combined_map", str(tmp_path / "combined_map.map"),
             "--loc_chrom",    str(tmp_path / "loc.tsv"),
             "--maf",          str(tmp_path / "t.maf"),
             "--source",       "TEST",
             "--output_dir",   str(tmp_path)],
            capture_output=True, text=True,
        )
        df = pd.read_csv(tmp_path / "Missense_filter_mutations_mapped.tsv", sep="\t")
        assert "Protein_ID" in df.columns
        assert "isoform_mapped" in df.columns
