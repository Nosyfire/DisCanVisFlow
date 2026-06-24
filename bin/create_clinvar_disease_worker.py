#!/usr/bin/env python3
"""
Module 8a — ClinVar Disease annotation with disease category classification.

Filters pre-processed ClinVar disease ontology table to proteins in this run,
then joins with the paper's disease category mapping to add Final_Category
(Cancer, Neurodegenerative, Cardiovascular/Hematopoietic, etc.).

Usage:
  create_clinvar_disease_worker.py
      --seq_table             <loc_chrom_with_names_isoforms_with_seq.tsv>
      --clinvar_disease        <clinvar_table_with_ontology_only_disease_pathogen_annotate.tsv>
      --clinvar_category_tsv  <clinvar_diseases.tsv  or  NO_FILE>
      --outdir                <output directory>

Output:
  clinvar_disease.tsv — Protein_ID + disease ontology columns + Final_Category
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

_KEEP_COLS = [
    "Protein_ID", "Disease", "DOID", "DO Subset", "synonyms",
    "PhenotypeList", "PhenotypeIDS", "xref_source", "reference_id",
    "level1", "level2", "level3", "level4", "level5", "level6",
    "level7", "level8", "level9", "level10", "level11", "level12", "level13",
    "Disordered", "Ordered", "Total Mutations", "Disordered Percent",
    "Final_Category",
]


def _build_category_map(cat_path: Path) -> dict:
    """Disease name → Final_Category from clinvar_diseases.tsv"""
    if not cat_path.exists() or cat_path.stat().st_size == 0 or cat_path.name == "NO_FILE":
        return {}
    try:
        df = pd.read_csv(cat_path, sep="\t", dtype=str,
                         usecols=["Disease", "Final_Category"])
        df = df.dropna(subset=["Disease", "Final_Category"])
        return dict(zip(df["Disease"].str.strip(), df["Final_Category"].str.strip()))
    except Exception as exc:
        log.warning("Could not load category map: %s", exc)
        return {}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table",            required=True)
    p.add_argument("--clinvar_disease",       required=True)
    p.add_argument("--clinvar_category_tsv",  required=True)
    p.add_argument("--outdir",               required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "clinvar_disease.tsv"

    empty = pd.DataFrame(columns=["Protein_ID", "Disease", "DOID", "Final_Category"])

    src = Path(args.clinvar_disease)
    if not src.exists() or src.stat().st_size == 0:
        log.info("ClinVar disease table not found — writing empty output")
        empty.to_csv(out, sep="\t", index=False)
        return

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
    protein_ids = set(seq_df["Protein_ID"].dropna())

    df = pd.read_csv(src, sep="\t", dtype=str)
    if "Accession" not in df.columns:
        log.warning("No 'Accession' column in ClinVar disease table")
        empty.to_csv(out, sep="\t", index=False)
        return

    df = df[df["Accession"].isin(protein_ids)].copy()
    df.rename(columns={"Accession": "Protein_ID"}, inplace=True)

    # Join disease category
    cat_map = _build_category_map(Path(args.clinvar_category_tsv))
    if cat_map:
        df["Final_Category"] = df["Disease"].str.strip().map(cat_map).fillna("Unknown")
        log.info("Category map applied: %d unique categories", df["Final_Category"].nunique())
    else:
        df["Final_Category"] = "Unknown"

    keep = [c for c in _KEEP_COLS if c in df.columns]
    df = df[keep].drop_duplicates()
    df.to_csv(out, sep="\t", index=False)
    log.info("ClinVar disease: %d rows for %d proteins", len(df), df["Protein_ID"].nunique())


if __name__ == "__main__":
    main()
