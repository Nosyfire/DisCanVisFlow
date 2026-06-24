#!/usr/bin/env python3
"""
Module 5h — Transfer PEM core motifs to alternative isoforms via sequence homology.

Takes filtered pem_core_motifs.tsv (main-isoform PEM from supplementary dataset)
and maps regional motifs to other GENCODE isoforms of the same gene using exact
substring search (same logic as TRANSCRIPT_MAP regional transfer).

Usage:
  create_pem_transfer_worker.py
      --loc_chrom   <loc_chrom_with_names_isoforms_with_seq.tsv>
      --pem_tsv     <pem_core_motifs.tsv>
      --outdir      <output directory>

Output:
  pem_core_motifs_mapped.tsv
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


def _build_lookup(loc_df: pd.DataFrame):
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"] if c in loc_df.columns), None)
    pid_to_seq = {}
    pid_to_gene = {}
    gene_to_rows = {}

    for _, row in loc_df.iterrows():
        pid = str(row.get("Protein_ID", "") or "")
        seq = str(row.get("Sequence", "") or "")
        gene = str(row.get(gene_col, "") or "") if gene_col else ""
        if pid and seq and seq not in ("nan", ""):
            pid_to_seq[pid] = seq
            pid_to_gene[pid] = gene
            gene_to_rows.setdefault(gene, []).append((pid, seq))
    return pid_to_seq, pid_to_gene, gene_to_rows


def transfer_pem(pem_df: pd.DataFrame, loc_df: pd.DataFrame) -> pd.DataFrame:
    pid_to_seq, pid_to_gene, gene_to_rows = _build_lookup(loc_df)
    if pem_df.empty:
        return pem_df

    pem_df = pem_df.copy()
    for col in ("Start", "End"):
        if col in pem_df.columns:
            pem_df[col] = pd.to_numeric(pem_df[col], errors="coerce")

    out_rows = []
    for _, row in pem_df.iterrows():
        src_pid = str(row.get("Protein_ID", ""))
        src_seq = pid_to_seq.get(src_pid, "")
        gene = pid_to_gene.get(src_pid, "")
        if not src_seq or not gene:
            continue

        try:
            start = int(row["Start"])
            end = int(row["End"])
        except (ValueError, TypeError, KeyError):
            out_rows.append(row.to_dict())
            continue

        if start < 1 or end > len(src_seq) or start > end:
            continue
        region = src_seq[start - 1:end]

        for tgt_pid, tgt_seq in gene_to_rows.get(gene, []):
            if not tgt_seq:
                continue
            if tgt_pid == src_pid:
                r = row.to_dict()
                r["homology_transfer"] = "False"
                out_rows.append(r)
                continue
            idx = tgt_seq.find(region)
            if idx == -1:
                continue
            new_start = idx + 1
            new_end = idx + len(region)
            if new_end > len(tgt_seq):
                continue
            r = row.to_dict()
            r["Protein_ID"] = tgt_pid
            r["Start"] = new_start
            r["End"] = new_end
            r["homology_transfer"] = "True"
            out_rows.append(r)

    if not out_rows:
        return pem_df
    return pd.DataFrame(out_rows).drop_duplicates()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--loc_chrom", required=True)
    p.add_argument("--pem_tsv", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)
    pem_path = Path(args.pem_tsv)
    if not pem_path.exists() or pem_path.stat().st_size == 0:
        pd.DataFrame(columns=["Protein_ID", "Start", "End"]).to_csv(
            outdir / "pem_core_motifs_mapped.tsv", sep="\t", index=False)
        return

    pem_df = pd.read_csv(pem_path, sep="\t", dtype=str)
    out_df = transfer_pem(pem_df, loc_df)
    out_df.to_csv(outdir / "pem_core_motifs_mapped.tsv", sep="\t", index=False)
    log.info("PEM transfer: %d rows (%d proteins)",
             len(out_df), out_df["Protein_ID"].nunique() if len(out_df) else 0)


if __name__ == "__main__":
    main()
