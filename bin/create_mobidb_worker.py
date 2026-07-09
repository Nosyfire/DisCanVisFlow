#!/usr/bin/env python3
"""
Module 5o — MobiDBDisorder standalone TSV output.

Converts the bulk MobiDB TSV (from FETCH_MOBIDB) to a per-protein feature
summary for the MobiDBDisorder Django model.

Input MobiDB bulk TSV format (from FETCH_MOBIDB / mobidb.bio.unipd.it):
    Headerless, tab-separated, one aggregated row per (acc, feature):
      acc          feature                  start..end                content_fraction  content_count  length
      P04637       curated-disorder-merge   1..96,288..312,361..393   0.392             154            393
    The `start..end` field may hold several comma-separated regions.
    FETCH_MOBIDB concatenates two per-feature downloads with `sort -u`, which
    leaves the source header row buried among the data (col 0 == "acc"); such
    rows are skipped. A legacy 4-column `acc/feature/source/start..end` layout
    is also accepted.

Output (mobidb_disorder.tsv):
    Protein_ID | Entry_Isoform | feature | start_end | content_fraction | content_count | length

Usage:
    create_mobidb_worker.py
        --seq_table   <loc_chrom_with_names_isoforms_with_seq.tsv>
        --mobidb_tsv  <mobidb_human.tsv>  (or NO_FILE)
        --outdir      <output directory>
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_OUT_COLS = [
    "Protein_ID", "Entry_Isoform", "feature",
    "start_end", "content_fraction", "content_count", "length",
]


def _parse_region(s: str) -> tuple[int, int] | None:
    """Parse a single '35-120' or '35..120' into (35, 120)."""
    s = str(s).strip()
    for sep in ("..", "-"):
        if sep in s:
            parts = s.split(sep, 1)
            try:
                return int(parts[0]), int(parts[1])
            except (ValueError, IndexError):
                pass
    return None


def _parse_regions(s: str) -> list[tuple[int, int]]:
    """Parse a possibly comma-separated region field.

    '1..96,288..312,361..393' → [(1, 96), (288, 312), (361, 393)].
    """
    out = []
    for chunk in str(s).split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        reg = _parse_region(chunk)
        if reg:
            out.append(reg)
    return out


def _read_mobidb(path: Path) -> pd.DataFrame:
    """Read the headerless MobiDB bulk TSV into acc/feature/region columns.

    Real FETCH_MOBIDB output is 6 columns with no header:
        acc, feature, start..end, content_fraction, content_count, length
    A legacy 4-column layout (acc, feature, source, start..end) is also read.
    A stray header row buried by `sort -u` (col 0 == "acc") is dropped later.
    """
    raw = pd.read_csv(path, sep="\t", dtype=str, header=None)
    ncol = raw.shape[1]
    if ncol >= 6:
        names = ["acc", "feature", "region", "content_fraction",
                 "content_count", "length"]
    elif ncol == 4:
        names = ["acc", "feature", "source", "region"]
    elif ncol == 3:
        names = ["acc", "feature", "region"]
    else:
        return pd.DataFrame(columns=["acc", "feature", "region"])
    # Keep only the leading named columns (ignore any extra trailing columns).
    raw = raw.iloc[:, :len(names)]
    raw.columns = names
    return raw


def build_mobidb_table(
    seq_df: pd.DataFrame,
    mob_df: pd.DataFrame,
) -> pd.DataFrame:
    """Return per-protein-feature rows mapped to Protein_ID."""
    # Build acc → [(Protein_ID, seq_length), ...]
    pid_col = next((c for c in ["Protein_ID", "Entry_Name"] if c in seq_df.columns), None)
    acc_col = next((c for c in ["Entry_Isoform", "Accession"] if c in seq_df.columns), None)
    if pid_col is None or acc_col is None:
        return pd.DataFrame(columns=_OUT_COLS)

    acc_to_pids: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for _, row in seq_df.iterrows():
        pid = str(row.get(pid_col, "")).strip()
        acc = str(row.get(acc_col, "")).strip()
        seq = str(row.get("Sequence", "")).strip()
        if not pid or pid == "nan" or not acc or acc == "nan":
            continue
        # Canonical base accession (P04049-2 → P04049)
        base_acc = acc.split("-")[0] if "-" in acc else acc
        seq_len  = len(seq) if seq and seq != "nan" else 0
        entry = (pid, acc, seq_len)
        if entry not in acc_to_pids[base_acc]:
            acc_to_pids[base_acc].append(entry)

    if mob_df.empty or "acc" not in mob_df.columns:
        return pd.DataFrame(columns=_OUT_COLS)

    region_col = "region" if "region" in mob_df.columns else None
    if region_col is None:
        return pd.DataFrame(columns=_OUT_COLS)

    # Group by (acc, feature) → list of (start, end) regions. Each cell may hold
    # several comma-separated regions, and an (acc, feature) pair may recur.
    regions_by_af: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for _, row in mob_df.iterrows():
        acc  = str(row.get("acc", "")).strip().split("-")[0]
        feat = str(row.get("feature", "")).strip()
        if not acc or acc.lower() == "acc" or not feat:
            continue
        regions_by_af[(acc, feat)].extend(_parse_regions(row.get(region_col, "")))

    rows = []
    for (acc, feat), regions in regions_by_af.items():
        if acc not in acc_to_pids or not regions:
            continue
        # De-duplicate and order regions; build "1-50,100-150" string.
        sorted_regions = sorted(set(regions), key=lambda x: x[0])
        start_end_str  = ",".join(f"{s}-{e}" for s, e in sorted_regions)
        covered        = sum(e - s + 1 for s, e in sorted_regions)

        for pid, full_acc, seq_len in acc_to_pids[acc]:
            cf = round(covered / seq_len, 4) if seq_len > 0 else 0.0
            rows.append({
                "Protein_ID":       pid,
                "Entry_Isoform":    full_acc,
                "feature":          feat,
                "start_end":        start_end_str,
                "content_fraction": cf,          # covered residues / isoform length
                "content_count":    covered,     # covered (disordered) residues
                "length":           seq_len,     # isoform sequence length
            })

    return pd.DataFrame(rows, columns=_OUT_COLS) if rows \
        else pd.DataFrame(columns=_OUT_COLS)


def main():
    p = argparse.ArgumentParser(
        description="Module 5o: MobiDB disorder → MobiDBDisorder TSV"
    )
    p.add_argument("--seq_table",  required=True)
    p.add_argument("--mobidb_tsv", required=True)
    p.add_argument("--outdir",     default=".")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out    = outdir / "mobidb_disorder.tsv"

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)

    mob_path = Path(args.mobidb_tsv)
    if mob_path.name == "NO_FILE" or not mob_path.exists() or mob_path.stat().st_size < 10:
        log.info("MobiDB TSV not available — writing empty output")
        pd.DataFrame(columns=_OUT_COLS).to_csv(out, sep="\t", index=False)
        return

    mob_df = _read_mobidb(mob_path)
    log.info("Loaded %d MobiDB region rows", len(mob_df))

    result = build_mobidb_table(seq_df, mob_df)
    result.to_csv(out, sep="\t", index=False)
    log.info("Done — %d MobiDB disorder rows written (%d proteins)",
             len(result),
             result["Protein_ID"].nunique() if not result.empty else 0)


if __name__ == "__main__":
    main()
