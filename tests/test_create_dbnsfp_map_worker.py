"""Tests for create_dbnsfp_map_worker.py (Module 8f — raw dbNSFP mapping).

Focus: the single merged dbNSFP 5.x file path (inverted (chr,pos) index +
pattern-based column selection: scores + rankscores + CADD + conservation +
gnomAD 4.1 joint AF).
"""

import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_dbnsfp_map_worker.py"
sys.path.insert(0, str(WORKER.parent))

import create_dbnsfp_map_worker as w  # noqa: E402


# --------------------------------------------------------------------------- #
# Unit tests — column selection                                               #
# --------------------------------------------------------------------------- #

_HEADER = [
    "#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt", "rs_dbSNP",
    "aapos", "genename", "SIFT_score", "SIFT_converted_rankscore", "SIFT_pred",
    "REVEL_score", "REVEL_rankscore", "AlphaMissense_score", "AlphaMissense_pred",
    "CADD_raw", "CADD_phred", "GERP++_RS", "phyloP100way_vertebrate",
    "phastCons17way_primate", "gnomAD4.1_joint_AF", "gnomAD4.1_joint_POPMAX_AF",
    "gnomAD4.1_joint_AFR_AF", "gnomAD2.1.1_exomes_AF",
]


class TestSelectKeepColumns:
    def test_chr_pos_indices(self):
        chr_idx, pos_idx, aaref_idx, file_cols = w.select_keep_columns(_HEADER)
        assert chr_idx == 0
        assert pos_idx == 1
        assert aaref_idx == 4

    def test_scores_and_rankscores_kept(self):
        _, _, _, file_cols = w.select_keep_columns(_HEADER)
        names = [n for n, _ in file_cols]
        assert "SIFT_score" in names
        assert "SIFT_converted_rankscore" in names
        assert "REVEL_score" in names
        assert "REVEL_rankscore" in names
        assert "AlphaMissense_score" in names

    def test_pred_columns_excluded(self):
        _, _, _, file_cols = w.select_keep_columns(_HEADER)
        names = [n for n, _ in file_cols]
        assert "SIFT_pred" not in names
        assert "AlphaMissense_pred" not in names

    def test_cadd_and_conservation_kept(self):
        _, _, _, file_cols = w.select_keep_columns(_HEADER)
        names = [n for n, _ in file_cols]
        assert "CADD_raw" in names
        assert "CADD_phred" in names
        assert "GERP++_RS" in names
        assert "phyloP100way_vertebrate" in names
        assert "phastCons17way_primate" in names

    def test_gnomad_joint_af_kept_others_dropped(self):
        _, _, _, file_cols = w.select_keep_columns(_HEADER)
        names = [n for n, _ in file_cols]
        assert "gnomAD4.1_joint_AF" in names
        assert "gnomAD4.1_joint_POPMAX_AF" in names
        # population-specific and old gnomAD blocks are NOT kept
        assert "gnomAD4.1_joint_AFR_AF" not in names
        assert "gnomAD2.1.1_exomes_AF" not in names

    def test_identity_columns_kept(self):
        _, _, _, file_cols = w.select_keep_columns(_HEADER)
        names = [n for n, _ in file_cols]
        for c in ["ref", "alt", "aaref", "aaalt", "aapos", "rs_dbSNP"]:
            assert c in names


class TestBuildGposIndex:
    def test_inverts_pid_map(self):
        pid_map = {
            "GENE1-201": {1000: (5, "A"), 1003: (6, "L")},
            "GENE2-201": {2000: (10, "M")},
        }
        pid_to_chrom = {"GENE1-201": "chr1", "GENE2-201": "chr2"}
        idx = w.build_gpos_index(pid_map, pid_to_chrom)
        assert idx[("chr1", 1000)] == [("GENE1-201", 5, "A")]
        assert idx[("chr2", 2000)] == [("GENE2-201", 10, "M")]

    def test_shared_position_groups_isoforms(self):
        pid_map = {
            "GENE1-201": {1000: (5, "A")},
            "GENE1-204": {1000: (5, "A")},
        }
        pid_to_chrom = {"GENE1-201": "chr1", "GENE1-204": "chr1"}
        idx = w.build_gpos_index(pid_map, pid_to_chrom)
        pids = {e[0] for e in idx[("chr1", 1000)]}
        assert pids == {"GENE1-201", "GENE1-204"}


# --------------------------------------------------------------------------- #
# End-to-end — single merged .gz file                                         #
# --------------------------------------------------------------------------- #

def _run(args, tmpdir):
    return subprocess.run([sys.executable, str(WORKER)] + args,
                          capture_output=True, text=True, cwd=tmpdir)


def _seq(tmpdir: Path) -> Path:
    p = tmpdir / "seq.tsv"
    p.write_text(
        "Protein_ID\tEntry_Isoform\tGene_Gencode\tChromosome\tSequence\n"
        "GENE1-201\tP11111\tGENE1\t1\t" + "MKALLV" + "\n",
        encoding="utf-8",
    )
    return p


def _combined_map(tmpdir: Path) -> Path:
    """combined_map.map: protein residue 3 (A) at genomic pos 1000, chr1.
    8 cols: protein_pos aa nuc_pos codon aa gpos_csv codon aa.
    load_combined_map_by_protein uses col0 (0-based prot pos), col1 aa, col5 gpos.
    """
    p = tmpdir / "combined_map.map"
    # header line names the protein via | pipe fields; field index 4 = Protein_ID
    p.write_text(
        "# GENE1|x|x|x|GENE1-201\n"
        "2\tA\t7\tGCT\tA\t1000,\tGCT\tA\n",
        encoding="utf-8",
    )
    return p


def _dbnsfp_gz(tmpdir: Path) -> Path:
    p = tmpdir / "dbNSFP_merged.gz"
    header = "\t".join([
        "#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt", "rs_dbSNP",
        "aapos", "SIFT_score", "REVEL_score", "REVEL_rankscore",
        "AlphaMissense_score", "SIFT_pred", "CADD_phred",
        "gnomAD4.1_joint_AF", "gnomAD4.1_joint_POPMAX_AF",
    ])
    # A at protein pos 3 (1-based) → aaref A ; genomic pos 1000
    row = "\t".join([
        "1", "1000", "G", "A", "A", "V", "rs777",
        "3", "0.02", "0.81", "0.90", "0.95", "D", "24.1", "0.001", "0.004",
    ])
    with gzip.open(p, "wt") as fh:
        fh.write(header + "\n")
        fh.write(row + "\n")
    return p


class TestEndToEndMergedFile:
    def test_maps_and_keeps_rankscore_and_gnomad(self, tmp_path):
        r = _run(["--seq_table", str(_seq(tmp_path)),
                  "--combined_map", str(_combined_map(tmp_path)),
                  "--dbnsfp_raw_dir", str(_dbnsfp_gz(tmp_path)),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "dbnsfp_scores.tsv", sep="\t", dtype=str)
        assert len(df) == 1, df.to_dict("records")
        row = df.iloc[0]
        assert row["Protein_ID"] == "GENE1-201"
        assert str(row["Protein_position"]) == "3"
        assert row["REVEL_score"] == "0.81"
        assert row["REVEL_rankscore"] == "0.90"
        assert row["AlphaMissense_score"] == "0.95"
        assert row["gnomAD4.1_joint_AF"] == "0.001"
        assert row["gnomAD4.1_joint_POPMAX_AF"] == "0.004"

    def test_pred_column_not_emitted(self, tmp_path):
        _run(["--seq_table", str(_seq(tmp_path)),
              "--combined_map", str(_combined_map(tmp_path)),
              "--dbnsfp_raw_dir", str(_dbnsfp_gz(tmp_path)),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "dbnsfp_scores.tsv", sep="\t", dtype=str)
        assert "SIFT_pred" not in df.columns

    def test_aaref_mismatch_filtered(self, tmp_path):
        """A dbNSFP row whose aaref disagrees with the mapped residue is dropped."""
        gz = tmp_path / "bad.gz"
        header = "\t".join(["#chr", "pos(1-based)", "ref", "alt", "aaref",
                            "aaalt", "rs_dbSNP", "aapos", "REVEL_score",
                            "gnomAD4.1_joint_AF"])
        row = "\t".join(["1", "1000", "G", "A", "W", "V", "rs1", "3", "0.5", "0.01"])
        with gzip.open(gz, "wt") as fh:
            fh.write(header + "\n" + row + "\n")
        _run(["--seq_table", str(_seq(tmp_path)),
              "--combined_map", str(_combined_map(tmp_path)),
              "--dbnsfp_raw_dir", str(gz),
              "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "dbnsfp_scores.tsv", sep="\t", dtype=str)
        assert len(df) == 0
