#!/usr/bin/env python3
"""
split_seq_table.py — split the sequence table (loc_chrom_with_names_isoforms_with_seq.tsv)
into N balanced chunks for scatter-parallel per-isoform steps (DISORDER_MAP,
COILEDCOILS_MAP). At 20 k genes these steps are the wall (serial per-protein
loops); chunking lets Nextflow run K tasks concurrently (maxForks / SLURM).

Splitting is **by gene**: every isoform of a gene stays in the same chunk, so
per-gene/per-accession work (e.g. the AlphaFold pLDDT call, keyed by canonical
accession) is not duplicated across chunks and any cross-isoform aggregation
inside a chunk is correct. Genes are distributed greedily (largest first) across
chunks to balance isoform counts. Each chunk keeps the header. Empty chunks are
not written, so the number of files emitted is min(n_chunks, n_genes).

Usage:
  split_seq_table.py --loc_chrom <seq.tsv> --n_chunks 200 [--prefix chunk_] [--outdir .]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

_GENE_COLS = ["Gene_Gencode", "Gene_Uniprot", "Gene"]


def split_by_gene(df: pd.DataFrame, n_chunks: int) -> list[pd.DataFrame]:
    """Return up to n_chunks DataFrames, partitioning rows by gene with balanced
    isoform counts (greedy largest-first bin packing)."""
    gene_col = next((c for c in _GENE_COLS if c in df.columns), None)
    if gene_col is None:
        # no gene column → fall back to contiguous row chunks
        groups = [(str(i), g) for i, g in enumerate(_contiguous(df, n_chunks))]
    else:
        groups = [(str(gene), sub) for gene, sub in df.groupby(gene_col, sort=False)]

    n = max(1, min(n_chunks, len(groups)))
    # greedy: assign each gene (largest first) to the currently-smallest bin
    order = sorted(groups, key=lambda kv: len(kv[1]), reverse=True)
    bins: list[list[pd.DataFrame]] = [[] for _ in range(n)]
    sizes = [0] * n
    for _gene, sub in order:
        j = sizes.index(min(sizes))
        bins[j].append(sub)
        sizes[j] += len(sub)
    out = []
    for parts in bins:
        if parts:
            out.append(pd.concat(parts, ignore_index=True))
    return out


def _contiguous(df: pd.DataFrame, n_chunks: int):
    n = max(1, min(n_chunks, len(df) or 1))
    size = -(-len(df) // n)  # ceil
    for i in range(0, len(df), size):
        yield df.iloc[i:i + size]


def main():
    ap = argparse.ArgumentParser(description="Split seq table into N gene-balanced chunks")
    ap.add_argument("--loc_chrom", required=True)
    ap.add_argument("--n_chunks", type=int, required=True)
    ap.add_argument("--prefix", default="chunk_")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str).fillna("")
    if df.empty:
        # still emit one (header-only) chunk so downstream has a file
        p = outdir / f"{args.prefix}001.tsv"
        df.to_csv(p, sep="\t", index=False)
        print(f"wrote 1 empty chunk → {p}", file=sys.stderr)
        return

    chunks = split_by_gene(df, max(1, args.n_chunks))
    width = max(3, len(str(len(chunks))))
    for i, sub in enumerate(chunks, 1):
        p = outdir / f"{args.prefix}{str(i).zfill(width)}.tsv"
        sub.to_csv(p, sep="\t", index=False)
    print(f"split {len(df)} rows / {df.shape} into {len(chunks)} chunk(s)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
