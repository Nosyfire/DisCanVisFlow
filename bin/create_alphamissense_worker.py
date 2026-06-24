#!/usr/bin/env python3
"""
Module 8d — AlphaMissense per-variant pathogenicity scores (GENCODE isoform version).

Computes from the raw AlphaMissense_isoforms_hg38.tsv.gz file (Google DeepMind).
Joins on ENST transcript ID (stripping .version suffix) to map transcript_id → Protein_ID.

Raw file columns (tab-separated, comment lines start with #):
  #CHROM  POS  REF  ALT  genome  transcript_id  protein_variant  am_pathogenicity  am_class

Output: alphamissense.tsv
  Protein_ID | transcript_id | protein_variant | am_pathogenicity | am_class

Usage:
  create_alphamissense_worker.py
      --seq_table          <loc_chrom_with_names_isoforms_with_seq.tsv>
      --alphamissense_gz   <AlphaMissense_isoforms_hg38.tsv.gz>
      --outdir             <output directory>
"""

import argparse
import gzip
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_CHUNKSIZE = 500_000
_OUT_COLS  = ["Protein_ID", "transcript_id", "protein_variant",
              "am_pathogenicity", "am_class"]
# Raw file column order (from #CHROM header line).
_KNOWN_COLS = ["CHROM", "POS", "REF", "ALT", "genome",
               "transcript_id", "protein_variant", "am_pathogenicity", "am_class"]
# Columns we actually need downstream.
_USE_COLS   = ["transcript_id", "protein_variant", "am_pathogenicity", "am_class"]


def _build_enst_map(seq_df: pd.DataFrame) -> dict:
    """Return dict: base_enst_id → Protein_ID (e.g. 'ENST00000251849' → 'RAF1-201')."""
    if "transcript_stable_id" in seq_df.columns:
        col = "transcript_stable_id"
    elif "Transcript ID" in seq_df.columns:
        col = "Transcript ID"
    else:
        raise ValueError("seq_table has no 'transcript_stable_id' or 'Transcript ID' column")

    out = {}
    for _, row in seq_df[["Protein_ID", col]].dropna().iterrows():
        enst_base = str(row[col]).split(".")[0].strip()
        if enst_base:
            out[enst_base] = row["Protein_ID"]
    return out


def _scan_plain_tsv(src: Path, enst_map: dict) -> tuple[list, int]:
    """Fast path for decompressed TSV: grep pre-filter (C-level) → parse small result.

    grep -F -f scans the 9.4 GB file at I/O speed and returns only matching lines.
    For a single gene (~9 ENSTs) this yields ~15k rows out of 144M → pandas parse
    is instant.  For the full proteome all rows match anyway, so grep is still the
    fastest way to stream the data into pandas without the Python-loop overhead.
    """
    # Read column names from the #CHROM comment line (fast — stops at first data line).
    header_cols = _KNOWN_COLS[:]
    with open(src, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM") or line.startswith("# CHROM"):
                header_cols = line.lstrip("#").strip().split("\t")
                break
            elif not line.startswith("#"):
                break

    use_cols = [c for c in _USE_COLS if c in header_cols]

    # Write ENST base IDs to a temp file for grep -F -f.
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("\n".join(enst_map.keys()))
        ids_path = Path(f.name)

    kept = []
    total_in = 0
    try:
        proc = subprocess.Popen(
            ["grep", "-F", "-f", str(ids_path), str(src)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        try:
            for chunk in pd.read_csv(
                proc.stdout, sep="\t",
                header=None, names=header_cols,
                usecols=use_cols,
                dtype=str, chunksize=_CHUNKSIZE, engine="c", na_filter=False,
            ):
                total_in += len(chunk)
                enst_base = chunk["transcript_id"].str.split(".", n=1).str[0]
                mask = enst_base.isin(enst_map)
                if mask.any():
                    c = chunk.loc[mask].copy()
                    c["Protein_ID"] = enst_base[mask].map(enst_map)
                    kept.append(c)
        except pd.errors.EmptyDataError:
            pass  # grep found no matches
        proc.wait()
    finally:
        ids_path.unlink(missing_ok=True)

    return kept, total_in


def _scan_gzip(src: Path, enst_map: dict) -> tuple[list, int]:
    """Fallback: line-by-line gzip scan (used if plain TSV not available)."""
    kept = []
    total_in = 0
    opener = gzip.open if src.suffix == ".gz" else open

    with opener(src, "rt") as fh:
        header = None
        chunk_rows = []

        for line in fh:
            if line.startswith("#"):
                if line.startswith("#CHROM") or line.startswith("# CHROM"):
                    header = line.lstrip("#").strip().split("\t")
                continue
            if header is None:
                continue

            chunk_rows.append(line.rstrip("\n").split("\t"))
            total_in += 1

            if len(chunk_rows) >= _CHUNKSIZE:
                chunk = pd.DataFrame(chunk_rows, columns=header)
                chunk_rows = []
                enst_base = chunk["transcript_id"].str.split(".").str[0]
                mask = enst_base.isin(enst_map)
                if mask.any():
                    c = chunk.loc[mask].copy()
                    c["Protein_ID"] = enst_base[mask].map(enst_map)
                    kept.append(c)

        if chunk_rows:
            chunk = pd.DataFrame(chunk_rows, columns=header)
            enst_base = chunk["transcript_id"].str.split(".").str[0]
            mask = enst_base.isin(enst_map)
            if mask.any():
                c = chunk.loc[mask].copy()
                c["Protein_ID"] = enst_base[mask].map(enst_map)
                kept.append(c)
            total_in += len(chunk_rows)

    return kept, total_in


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table",       required=True)
    p.add_argument("--alphamissense_gz", required=True)   # accepts .tsv or .gz
    p.add_argument("--outdir",          required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "alphamissense.tsv"

    empty = pd.DataFrame(columns=_OUT_COLS)

    src = Path(args.alphamissense_gz)
    if not src.exists() or src.stat().st_size == 0 or src.name == "NO_FILE":
        log.info("AlphaMissense raw file not found — writing empty output")
        empty.to_csv(out, sep="\t", index=False)
        return

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)
    enst_map = _build_enst_map(seq_df)

    if not enst_map:
        log.warning("No ENST IDs found in seq_table — writing empty output")
        empty.to_csv(out, sep="\t", index=False)
        return

    log.info("Matching %d ENST IDs from seq_table (file: %s)", len(enst_map), src.name)

    if src.suffix != ".gz":
        kept, total_in = _scan_plain_tsv(src, enst_map)
    else:
        kept, total_in = _scan_gzip(src, enst_map)

    if kept:
        df = pd.concat(kept, ignore_index=True)
        keep_cols = [c for c in _OUT_COLS if c in df.columns]
        df = df[keep_cols]
    else:
        df = empty

    df.to_csv(out, sep="\t", index=False)
    log.info("AlphaMissense: %d rows for %d proteins (scanned %d lines)",
             len(df), df["Protein_ID"].nunique() if len(df) else 0, total_in)


if __name__ == "__main__":
    main()
