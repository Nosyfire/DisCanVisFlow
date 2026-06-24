#!/usr/bin/env python3
"""
Module 5o — MobiDBDisorder standalone TSV output.

Converts the bulk MobiDB TSV (from FETCH_MOBIDB) to a per-protein feature
summary for the MobiDBDisorder Django model.

Input MobiDB bulk TSV format (from mobidb.bio.unipd.it):
    acc   feature   source   start..end
    P12345  curated-disorder-merge  mobidb_curated  35-120
    ...

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
    """Parse '35-120' or '35..120' into (35, 120)."""
    for sep in ("-", ".."):
        if sep in str(s):
            parts = str(s).split(sep, 1)
            try:
                return int(parts[0]), int(parts[1])
            except (ValueError, IndexError):
                pass
    return None


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

    # Normalise start..end column name (may have dots in header)
    region_col = next((c for c in mob_df.columns if "end" in c.lower()), None)
    if region_col is None:
        return pd.DataFrame(columns=_OUT_COLS)

    # Group by (acc, feature) → list of (start, end) regions
    feature_col = "feature" if "feature" in mob_df.columns else mob_df.columns[1]
    regions_by_af: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for _, row in mob_df.iterrows():
        acc  = str(row.get("acc", "")).strip().split("-")[0]
        feat = str(row.get(feature_col, "")).strip()
        reg  = _parse_region(row.get(region_col, ""))
        if acc and feat and reg:
            regions_by_af[(acc, feat)].append(reg)

    rows = []
    for (acc, feat), regions in regions_by_af.items():
        if acc not in acc_to_pids:
            continue
        # Build start_end string: "1-50,100-150"
        sorted_regions = sorted(regions, key=lambda x: x[0])
        start_end_str  = ",".join(f"{s}-{e}" for s, e in sorted_regions)
        total_len      = sum(e - s + 1 for s, e in sorted_regions)

        for pid, full_acc, seq_len in acc_to_pids[acc]:
            cf = round(total_len / seq_len, 4) if seq_len > 0 else 0.0
            rows.append({
                "Protein_ID":       pid,
                "Entry_Isoform":    full_acc,
                "feature":          feat,
                "start_end":        start_end_str,
                "content_fraction": cf,
                "content_count":    len(sorted_regions),
                "length":           total_len,
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

    mob_df = pd.read_csv(mob_path, sep="\t", dtype=str)
    log.info("Loaded %d MobiDB region rows", len(mob_df))

    result = build_mobidb_table(seq_df, mob_df)
    result.to_csv(out, sep="\t", index=False)
    log.info("Done — %d MobiDB disorder rows written (%d proteins)",
             len(result),
             result["Protein_ID"].nunique() if not result.empty else 0)


if __name__ == "__main__":
    main()
