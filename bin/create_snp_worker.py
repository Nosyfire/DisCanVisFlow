#!/usr/bin/env python3
"""
Module 5l — Population-level SNP / common polymorphism annotation.

Reads the pre-processed positional polymorphism table
(format: AccessionPosition = Protein_ID|pos, Polymorphism = Common/All)
and filters to proteins in this run.

Data source: DisCanVis_Data_Process positional_data_process/polymorphism_pos.tsv

Usage:
  create_snp_worker.py
      --seq_table    <loc_chrom_with_names_isoforms_with_seq.tsv>
      --snp_pos_tsv  <polymorphism_pos.tsv  or  NO_FILE>
      --outdir       <output directory>

Output:
  snp_polymorphisms.tsv  (Protein_ID, Position, Polymorphism)
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

_CHUNKSIZE = 500_000


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table",   required=True)
    p.add_argument("--snp_pos_tsv", required=True)
    p.add_argument("--outdir",      required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "snp_polymorphisms.tsv"

    empty = pd.DataFrame(columns=["Protein_ID", "Position", "Polymorphism"])

    src = Path(args.snp_pos_tsv)
    if not src.exists() or src.stat().st_size == 0 or src.name == "NO_FILE":
        log.info("SNP polymorphism table not found — writing empty output")
        empty.to_csv(out, sep="\t", index=False)
        return

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
    protein_ids = set(seq_df["Protein_ID"].dropna())

    kept = []
    for chunk in pd.read_csv(src, sep="\t", dtype=str, chunksize=_CHUNKSIZE):
        # AccessionPosition format: "Protein_ID|pos"
        if "AccessionPosition" in chunk.columns:
            split = chunk["AccessionPosition"].str.split("|", n=1, expand=True)
            chunk = chunk.copy()
            chunk["Protein_ID"] = split[0]
            chunk["Position"]   = split[1]
            chunk = chunk.drop(columns=["AccessionPosition"])
        mask = chunk["Protein_ID"].isin(protein_ids)
        if mask.any():
            kept.append(chunk.loc[mask, ["Protein_ID", "Position", "Polymorphism"]])

    df = pd.concat(kept, ignore_index=True) if kept else empty
    df.to_csv(out, sep="\t", index=False)
    log.info("SNP polymorphisms: %d rows for %d proteins",
             len(df), df["Protein_ID"].nunique() if len(df) else 0)


if __name__ == "__main__":
    main()
