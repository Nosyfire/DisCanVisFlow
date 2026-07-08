"""Tests for bin/create_mapping_report_worker.py — comprehensive mapping audit."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_mapping_report_worker.py"


def _build_run(tmp_path):
    """Create a minimal final/ + intermediate/ tree for gene GENE (2 isoforms)."""
    base = tmp_path
    final = base / "final"
    inter = base / "intermediate" / "annotations"
    (final / "annotations").mkdir(parents=True)
    (final / "disorder").mkdir(parents=True)
    (final / "structure").mkdir(parents=True)
    (final / "pathogenicity").mkdir(parents=True)
    (final / "genome").mkdir(parents=True)
    (final / "mutations" / "ClinVar").mkdir(parents=True)
    (final / "sequence").mkdir(parents=True)
    inter.mkdir(parents=True)

    # sequence table: 2 isoforms of GENE
    seq = final / "sequence" / "loc_chrom_with_names_isoforms_with_seq.tsv"
    pd.DataFrame({
        "Protein_ID":    ["GENE-201", "GENE-202"],
        "Entry_Isoform": ["P00001",   "P00001-2"],
        "Gene":          ["GENE",     "GENE"],
        "main_isoform":  ["yes",      "no"],
        "coverage":      ["100.0",    "98.0"],
        "alignmentpuntcuality": ["identical", "high"],
        "Chromosome":    ["chr1",     "chr1"],
    }).to_csv(seq, sep="\t", index=False)

    # combined_map: only GENE-201 gets a genomic location
    (final / "genome" / "combined_map.map").write_text(
        "# ENST1.1|ENSG1|OTT1|OTT2|GENE-201|GENE|999|CDS:1-6| chr1 + 100-200\n"
        "0 M 1,2,3 ATG M 100,101,102, ATG M\n")

    # raw ELM (2 motifs on P00001)
    pd.DataFrame({
        "Protein_ID": ["X", "X"], "Entry_Isoform": ["P00001", "P00001"],
        "ELMType": ["DOC", "LIG"], "ELMIdentifier": ["DOC_AAA", "LIG_BBB"],
        "Start": ["10", "20"], "End": ["15", "25"],
    }).to_csv(inter / "elm.tsv", sep="\t", index=False)

    # mapped ELM: only DOC_AAA mapped to GENE-201 (direct); LIG_BBB dropped
    pd.DataFrame({
        "Protein_ID": ["GENE-201"], "mapping_type": ["direct"],
        "homology_transfer": ["False"], "homology_identity": ["1.000"],
        "Entry_Isoform": ["P00001"], "ELMType": ["DOC"],
        "ELMIdentifier": ["DOC_AAA"], "Start": ["10"], "End": ["15"],
    }).to_csv(final / "annotations" / "elm.tsv", sep="\t", index=False)

    # disorder: AIUPred-Binding present only for GENE-201 (coverage gap)
    pd.DataFrame({"Protein_ID": ["GENE-201"],
                  "AIUPredBinding": ["0.1,0.2,0.3"]}).to_csv(
        final / "disorder" / "AIUPredBinding.tsv", sep="\t", index=False)
    # structure/pdb: both isoforms
    pd.DataFrame({"Protein_ID": ["GENE-201", "GENE-202"],
                  "pdb_id": ["1abc", "1abc"]}).to_csv(
        final / "structure" / "pdb_structures.tsv", sep="\t", index=False)
    # pathogenicity (mavedb lives here)
    pd.DataFrame({"Protein_ID": ["GENE-201"], "score": ["1.5"]}).to_csv(
        final / "pathogenicity" / "mavedb.tsv", sep="\t", index=False)
    # mutations under a source subdir
    pd.DataFrame({"Protein_ID": ["GENE-201", "GENE-202"],
                  "Mutation": ["p.A1V", "p.A1V"]}).to_csv(
        final / "mutations" / "ClinVar" / "Missense_filter_mutations_mapped.tsv",
        sep="\t", index=False)

    out = base / "out"
    subprocess.run([sys.executable, str(BIN),
                    "--seq_table", str(seq),
                    "--final_dir", str(final),
                    "--intermediate_dir", str(base / "intermediate"),
                    "--mapping_mode", "all_isoform_mapping",
                    "--command", "nextflow run main.nf -profile test_one_protein,conda",
                    "--pipeline_version", "0.5.0", "--nextflow_version", "26.04.3",
                    "--source", "pathogenicity/mavedb.tsv=local|data/mavedb.tsv",
                    "--outdir", str(out)], check=True)
    return out


def test_per_gene_report_created(tmp_path):
    out = _build_run(tmp_path)
    rep = out / "GENE_mapping_report.md"
    assert rep.exists()
    md = rep.read_text()
    assert "# Mapping report — GENE" in md
    assert "Isoforms selected:** 2" in md


def test_coverage_lists_all_categories_and_gaps(tmp_path):
    out = _build_run(tmp_path)
    md = (out / "GENE_mapping_report.md").read_text()
    # AIUPred-Binding present only for GENE-201 → 1/2, GENE-202 flagged missing
    assert "AIUPred-Binding" in md
    assert "GENE-202" in md
    # PDB structures cover both
    assert "PDB structures" in md
    # mutation source folded into label
    assert "ClinVar: Missense mutations" in md
    # pathogenicity (mavedb in pathogenicity folder)
    assert "MaveDB functional scores" in md
    # coverage table now carries a Source column with provenance
    assert "AIUPred-Binding" in md and "computed" in md
    assert "PDBe API" in md and "downloaded" in md


def test_before_after_unmapped(tmp_path):
    out = _build_run(tmp_path)
    md = (out / "GENE_mapping_report.md").read_text()
    assert "1 of 2" in md           # only DOC_AAA mapped
    assert "NOT mapped" in md
    assert "LIG_BBB" in md


def test_summary_reproducibility_and_locations(tmp_path):
    out = _build_run(tmp_path)
    s = (out / "mapping_summary.md").read_text()
    assert "# Mapping summary" in s
    assert "nextflow run main.nf -profile test_one_protein,conda" in s
    assert "0.5.0" in s and "26.04.3" in s
    # output locations listed with RELATIVE paths
    assert "## Output locations" in s
    assert "`final/pathogenicity/`" in s


def test_summary_overview_main_nonmain(tmp_path):
    out = _build_run(tmp_path)
    s = (out / "mapping_summary.md").read_text()
    assert "## Mapping overview (all annotations)" in s
    # scope header splits main / non-main
    assert "main: 1, non-main: 1" in s
    # every annotation is an independent row with a Source column
    assert "AIUPred-Binding" in s
    assert "PDB structures" in s and "PDBe API" in s
    # AIUPred-Binding only on the main isoform (GENE-201) → 1/1 main, 0/1 non-main
    # Search within the mapping overview section (the annotation sources section also
    # has an AIUPred-Binding row with a date, which would cause a false match on [0]).
    overview_section = s[s.find("## Mapping overview"):]
    line = [ln for ln in overview_section.splitlines() if "AIUPred-Binding" in ln][0]
    assert "1 / 1" in line and "0 / 1" in line
    # --source override (rel-path keyed) is honoured for MaveDB
    assert "data/mavedb.tsv" in s
    # two source columns now: origin + type
    assert "Source Type" in s


def test_genomic_location_gap_flagged(tmp_path):
    out = _build_run(tmp_path)
    md = (out / "GENE_mapping_report.md").read_text()
    assert "100-200" in md                       # GENE-201 located
    assert "without a genomic location" in md    # GENE-202 not in combined_map


def test_parse_map_regions_unit(tmp_path):
    sys.path.insert(0, str(BIN.parent))
    import importlib
    mod = importlib.import_module("create_mapping_report_worker")
    cmap = tmp_path / "cm.map"
    cmap.write_text(
        "# ENST1.1|ENSG1|OTT1|OTT2|GENE-201|GENE|999|CDS:1-6| chr1 + 100-200\n"
        "0 M 1,2,3 ATG M 100,101,102, ATG M\n")
    regions = mod.parse_map_regions(cmap)
    assert regions["GENE-201"] == ("chr1", "+", "100", "200")
