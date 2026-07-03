"""Tests for the provenance / base-statistics enrichment in
bin/create_mapping_report_worker.py (data-source versions, input-scale counts,
genome-mapped fallback)."""
import gzip
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_mapping_report_worker.py"


def _load_mod():
    spec = importlib.util.spec_from_file_location("mapping_report_worker", BIN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load_mod()


# ── count_fasta_entries ──────────────────────────────────────────────────────
def test_count_fasta_plain(tmp_path):
    fa = tmp_path / "x.fasta"
    fa.write_text(">a\nMKV\n>b\nAAA\n>c\nGGG\n")
    assert M.count_fasta_entries(str(fa)) == 3


def test_count_fasta_gz(tmp_path):
    fa = tmp_path / "x.fasta.gz"
    with gzip.open(fa, "wt") as fh:
        fh.write(">a\nMKV\n>b\nAAA\n")
    assert M.count_fasta_entries(str(fa)) == 2


def test_count_fasta_missing_returns_none():
    assert M.count_fasta_entries("/no/such/file.fa") is None
    assert M.count_fasta_entries(None) is None


# ── parse_gencode_version ────────────────────────────────────────────────────
def test_parse_gencode_version():
    assert M.parse_gencode_version(
        "/data/gencode_process/gencode.v44.pc_translations.fa") == "v44"
    assert M.parse_gencode_version(
        "/x/gencode.v46lift37.annotation.gtf.gz") == "v46lift37"
    assert M.parse_gencode_version("/x/random.fasta") is None
    assert M.parse_gencode_version(None) is None


# ── genome_mapped_pids: fallback to genome_protein_index.tsv ──────────────────
def test_genome_mapped_prefers_regions(tmp_path):
    (tmp_path / "genome").mkdir()
    regions = {"GENE-201": ("chr1", "+", "1", "9")}
    assert M.genome_mapped_pids(tmp_path, regions) == {"GENE-201"}


def test_genome_mapped_fallback_to_index(tmp_path):
    g = tmp_path / "genome"
    g.mkdir()
    pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "gpos": [100, 103, 200],
        "strand": ["+", "+", "-"],
        "Protein_ID": ["GENE-201", "GENE-201", "GENE-202"],
        "prot_pos": [1, 2, 1],
        "codon_offset": [0, 0, 0],
        "aa": ["M", "K", "A"],
        "codon": ["ATG", "AAA", "GCA"],
    }).to_csv(g / "genome_protein_index.tsv", sep="\t", index=False)
    # empty regions → must fall back to the index file
    assert M.genome_mapped_pids(tmp_path, {}) == {"GENE-201", "GENE-202"}


# ── compute_input_scale ──────────────────────────────────────────────────────
def _seq_df():
    return pd.DataFrame({
        "_pid":  ["G1-201", "G1-202", "G2-201"],
        "_gene": ["G1", "G1", "G2"],
        "_main": ["yes", "no", "yes"],
        "Database": ["Uniprot/SWISSPROT", "Uniprot_isoform", "Uniprot/SWISSPROT"],
    })


def test_compute_input_scale_counts(tmp_path):
    uni = tmp_path / "uni.fasta"
    uni.write_text("".join(f">sp|P{i}|X\nMK\n" for i in range(5)))
    gen = tmp_path / "gencode.v44.pc_translations.fa"
    gen.write_text("".join(f">e{i}\nMK\n" for i in range(9)))

    scale = M.compute_input_scale(
        _seq_df(), tmp_path, {"G1-201": 1, "G2-201": 1},
        gencode_fasta=str(gen), uniprot_fasta=str(uni))
    d = dict(scale)
    # reference sizes
    assert d["UniProt SwissProt entries (reference)"] == 5
    assert d["GENCODE protein-coding entries (reference)"] == 9
    # run-derived
    assert d["Genes mapped in run"] == 2
    assert d["Transcripts assigned in run"] == 3
    assert d["  via SwissProt canonical (direct)"] == 2
    assert d["  via curated isoform"] == 1
    assert d["Genome-mapped isoforms"] == 2


# ── end-to-end: summary must contain the new sections ────────────────────────
def test_summary_has_provenance_and_scale(tmp_path):
    final = tmp_path / "final"
    (final / "sequence").mkdir(parents=True)
    (final / "genome").mkdir(parents=True)
    (final / "annotations").mkdir(parents=True)
    seq = final / "sequence" / "loc_chrom_with_names.tsv"
    pd.DataFrame({
        "Protein_ID": ["G1-201", "G1-202"],
        "Entry_Isoform": ["P1", "P1-2"],
        "Gene": ["G1", "G1"],
        "main_isoform": ["yes", "no"],
        "coverage": ["100", "98"],
        "alignmentpuntcuality": ["identical", "high"],
        "Chromosome": ["chr1", "chr1"],
        "Database": ["Uniprot/SWISSPROT", "Uniprot_isoform"],
    }).to_csv(seq, sep="\t", index=False)
    pd.DataFrame({
        "chrom": ["chr1"], "gpos": [1], "strand": ["+"],
        "Protein_ID": ["G1-201"], "prot_pos": [1], "codon_offset": [0],
        "aa": ["M"], "codon": ["ATG"],
    }).to_csv(final / "genome" / "genome_protein_index.tsv", sep="\t", index=False)
    pd.DataFrame({"Protein_ID": ["G1-201"], "ELMType": ["DOC"]}).to_csv(
        final / "annotations" / "elm.tsv", sep="\t", index=False)

    gen = tmp_path / "gencode.v44.pc_translations.fa"
    gen.write_text(">e\nMK\n")
    uni = tmp_path / "uni.fasta"
    uni.write_text(">sp|P1|X\nMK\n")

    out = tmp_path / "reports"
    subprocess.run([sys.executable, str(BIN),
                    "--seq_table", str(seq),
                    "--final_dir", str(final),
                    "--outdir", str(out),
                    "--gencode_fasta", str(gen),
                    "--uniprot_fasta", str(uni),
                    "--mapping_mode", "all_isoform_mapping"],
                   check=True)
    summary = (out / "mapping_summary.md").read_text()
    assert "## Input scale" in summary
    assert "GENCODE" in summary and "v44" in summary
    assert "Genome-mapped isoforms" in summary
    # genome-mapped must be 1 (fallback via index), NOT 0
    assert "Data source versions" in summary

    # per-gene report must agree with the summary: G1-201 is genome-mapped via
    # the index fallback → ✅ and "1 / 2", not ❌ / "0 / 2"
    pg = (out / "G1_mapping_report.md").read_text()
    assert "Isoforms with a genomic location:** 1 / 2" in pg
    assert "❌" in pg  # G1-202 is NOT mapped
    assert "| G1-201 |" in pg and "✅" in pg

    # public-safe: no absolute filesystem path leaks into the summary
    assert str(tmp_path) not in summary
    assert "Launch dir" not in summary          # replaced by relative Output directory
    assert "Output locations" in summary
    # output-locations line is relative: "<run-dir-name>/final/", not absolute
    assert f"All outputs are under `{tmp_path.name}/final/`" in summary
    assert f"| Output directory | `{tmp_path.name}/` |" in summary
