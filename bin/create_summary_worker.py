#!/usr/bin/env python3
"""
Summary — annotation count report per gene.

Walks a results directory (or individual --file arguments) and
produces a two-column table: annotation_type, count.

Usage:
  create_summary_worker.py
      --gene_name  RAF1
      --results_dir  results/chr3/raf1/
      --outdir       results/chr3/raf1/
      [--file KEY=path ...]   # override / add individual files

Output:
  annotation_summary.tsv  (gene, annotation_type, count, note)
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

# Map (relative path glob patterns) → (annotation_type label, count_col)
# count_col = None → count rows
_FILE_MAP = [
    # ── Sequence ──────────────────────────────────────────────────────────────
    ("sequence/loc_chrom_with_names_isoforms_with_seq.tsv",
     "Isoforms processed", None),

    # ── Genome / Exon ─────────────────────────────────────────────────────────
    ("genome/combined_map.map", "Genome mapped residues", None),
    ("genome/exon.tsv",         "Exon boundaries",        None),

    # ── Mutations (mapped/mutations/<source>/) ────────────────────────────────
    ("final/mutations/*/Missense_filter_mutations_mapped.tsv",   "Missense mutations",   None),
    ("final/mutations/*/Frameshift_filter_mutations_mapped.tsv", "Frameshift mutations", None),
    ("final/mutations/*/Nonsense_filter_mutations_mapped.tsv",   "Nonsense mutations",   None),
    ("final/mutations/*/Indel_filter_mutations_mapped.tsv",      "Indel mutations",      None),
    ("final/mutations/DepMap/depmap_mutations.tsv",              "DepMap mutations",     None),

    # ── Unmapped annotations (Entry_Isoform-keyed — inputs to TRANSCRIPT_MAP) ─
    ("intermediate/annotations/elm.tsv",           "ELM motifs (raw)",      None),
    ("intermediate/annotations/dibs.tsv",          "DIBS sites (raw)",      None),
    ("intermediate/annotations/mfib.tsv",          "MFIB sites (raw)",      None),
    ("intermediate/annotations/phasepro.tsv",      "PhasePro entries (raw)",None),
    ("intermediate/annotations/ptm_merged.tsv",    "PTM sites (raw)",       None),
    ("intermediate/annotations/pfam_domains.tsv",  "Pfam domains (raw)",    None),

    # ── Mapped annotations (Protein_ID-keyed) ────────────────────────────────
    ("final/annotations/elm.tsv",           "ELM motifs",           None),
    ("final/annotations/dibs.tsv",          "DIBS sites",           None),
    ("final/annotations/mfib.tsv",          "MFIB sites",           None),
    ("final/annotations/phasepro.tsv",      "PhasePro entries",     None),
    ("final/annotations/ptm_merged.tsv",    "PTM sites",            None),
    ("final/annotations/pfam_domains.tsv",  "Pfam domains",         None),
    ("final/annotations/go_terms.tsv",      "GO terms",             None),
    ("final/annotations/polymorphism.tsv",  "Polymorphisms (all + allele freq)", None),
    ("final/annotations/pem_core_motifs.tsv","PEM core motifs",     None),
    ("final/annotations/coiled_coils.tsv",  "Coiled-coil regions",  None),
    ("final/annotations/interactions.tsv",  "PPI interactions",     None),
    ("final/annotations/scansite.tsv",      "ScanSite motifs",      None),

    # ── Pathogenicity / functional scores ─────────────────────────────────────
    ("final/pathogenicity/mavedb.tsv",       "MaveDB functional scores", None),
    ("final/pathogenicity/proteingym.tsv",   "ProteinGym DMS scores", None),

    # ── Disorder ──────────────────────────────────────────────────────────────
    ("final/disorder/IUPredscores.tsv",        "IUPred3 residues",          None),
    ("final/disorder/CombinedDisorderNew.tsv", "Combined disorder regions", None),

    # ── Structure (AlphaFold pLDDT + RSA + PDB) ───────────────────────────────
    ("final/structure/AlphaFoldTable.tsv",  "AlphaFold pLDDT residues",  None),
    ("final/structure/rsa_scores.tsv",      "RSA residues (from pLDDT)", None),
    ("final/structure/pdb_structures.tsv",  "PDB structures",   None),
    ("final/structure/pdb_missing.tsv",     "PDB missing residues (disorder)", None),

    # ── Conservation (mapped/) ────────────────────────────────────────────────
    ("final/conservation/conservation_multiple_level.tsv",
     "GOPHER conservation entries", None),
    ("final/conservation/conservation_phastcons.tsv",
     "phastCons residues", None),

    # ── Disease (mapped/) ─────────────────────────────────────────────────────
    ("final/disease/clinvar_disease.tsv",      "ClinVar disease entries",  None),
    ("final/disease/omim_disease.tsv",         "OMIM disease entries",     None),

    # ── Drivers (mapped/) ────────────────────────────────────────────────────
    ("final/drivers/census_driver.tsv",        "CGC census driver entries",None),
    ("final/drivers/compendium_driver.tsv",    "Compendium driver entries",None),

    # ── Pathogenicity (mapped/) ───────────────────────────────────────────────
    ("final/pathogenicity/alphamissense.tsv",        "AlphaMissense variants",   None),
    ("final/pathogenicity/pathogenicity_scores.tsv", "dbNSFP scored variants",   None),
]


def _count_file(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        if path.suffix == ".map":
            return sum(1 for line in path.open()
                       if line.strip() and not line.startswith("#"))
        df = pd.read_csv(path, sep="\t", dtype=str, nrows=None)
        return len(df)
    except Exception:
        return 0


def _resolve(base: Path, pattern: str) -> list:
    """Expand glob patterns, return all matching Paths."""
    paths = list(base.glob(pattern))
    if not paths:
        direct = base / pattern
        if direct.exists():
            return [direct]
    return paths


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gene_name",   required=True)
    p.add_argument("--results_dir", required=True)
    p.add_argument("--outdir",      required=True)
    p.add_argument("--file",        nargs="*", default=[],
                   metavar="KEY=PATH",
                   help="Override: 'ELM motifs=/path/to/elm.tsv'")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base = Path(args.results_dir)

    overrides = {}
    for kv in (args.file or []):
        k, _, v = kv.partition("=")
        overrides[k.strip()] = Path(v.strip())

    rows = []
    for pattern, label, _ in _FILE_MAP:
        if label in overrides:
            paths = [overrides[label]]
        else:
            paths = _resolve(base, pattern)

        if not paths:
            rows.append({"gene": args.gene_name, "annotation_type": label,
                         "count": 0, "note": "file not found"})
            continue

        total = sum(_count_file(fp) for fp in paths)
        rows.append({"gene": args.gene_name, "annotation_type": label,
                     "count": total, "note": ""})

    df = pd.DataFrame(rows, columns=["gene", "annotation_type", "count", "note"])
    df.to_csv(outdir / "annotation_summary.tsv", sep="\t", index=False)
    log.info("Summary: %d annotation types for %s", len(df), args.gene_name)
    log.info("\n%s", df.to_string(index=False))


if __name__ == "__main__":
    main()
