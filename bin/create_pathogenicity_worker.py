#!/usr/bin/env python3
"""
Module 8f — Pathogenicity predictor scores (dbNSFP custom).

Filters the pre-processed dbNSFP variant table (Protein_ID-keyed)
to proteins in this run.

Predictor columns included:
  AlphaMissense_score, CADD_phred, CADD_raw, ClinPred_score,
  ESM1b_score, EVE_score, Polyphen2_HDIV_score, Polyphen2_HVAR_score,
  PrimateAI_score, SIFT_score, VARITY_ER_LOO_score, VARITY_R_LOO_score,
  REVEL_score, gMVP_score

Usage:
  create_pathogenicity_worker.py
      --seq_table       <loc_chrom_with_names_isoforms_with_seq.tsv>
      --dbnsfp_tsv      <dbNSFP_custom/mapped_filtered_mutations.tsv>
      --outdir          <output directory>

Output:
  pathogenicity_scores.tsv
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
    "Protein_ID", "chr", "Start_Position", "End_Position",
    "Protein_position", "aaref", "aaalt", "aapos",
    "ref", "alt", "rs_dbSNP",
    "AlphaMissense_score", "CADD_phred", "CADD_raw",
    "ClinPred_score", "ESM1b_score", "EVE_score",
    "Polyphen2_HDIV_score", "Polyphen2_HVAR_score",
    "PrimateAI_score", "SIFT_score",
    "VARITY_ER_LOO_score", "VARITY_ER_score",
    "VARITY_R_LOO_score", "VARITY_R_score",
    "REVEL_score", "REVEL_rankscore",
    "gMVP_score",
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table", required=True)
    p.add_argument("--dbnsfp_tsv", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "pathogenicity_scores.tsv"

    empty = pd.DataFrame(columns=["Protein_ID", "Protein_position", "AlphaMissense_score"])

    src = Path(args.dbnsfp_tsv)
    if not src.exists() or src.stat().st_size == 0:
        log.info("dbNSFP table not found — writing empty output")
        empty.to_csv(out, sep="\t", index=False)
        return

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
    protein_ids = set(seq_df["Protein_ID"].dropna())

    kept = []
    total_in = 0
    for chunk in pd.read_csv(src, sep="\t", dtype=str, chunksize=200_000):
        total_in += len(chunk)
        if "Protein_ID" not in chunk.columns:
            log.warning("dbNSFP table has no 'Protein_ID' column — writing empty output")
            empty.to_csv(out, sep="\t", index=False)
            return
        c = chunk[chunk["Protein_ID"].isin(protein_ids)].copy()
        if c.empty:
            continue
        keep = [col for col in _KEEP_COLS if col in c.columns]
        kept.append(c[keep])

    if kept:
        df = pd.concat(kept, ignore_index=True)
    else:
        df = pd.DataFrame(columns=["Protein_ID", "Protein_position", "AlphaMissense_score"])

    df.to_csv(out, sep="\t", index=False)
    log.info("Pathogenicity: %d rows for %d proteins (scanned %d lines)",
             len(df), df["Protein_ID"].nunique() if len(df) else 0, total_in)


if __name__ == "__main__":
    main()
