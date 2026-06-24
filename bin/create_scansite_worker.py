#!/usr/bin/env python3
"""
Module 5k — ScanSite 4.0 phosphorylation / kinase motif predictions.

Two modes:
  1. Pre-processed (default): filter the pre-computed proteome-wide scansite.tsv
  2. API mode (--use_api): call https://scansite4.mit.edu for each protein

Usage:
  create_scansite_worker.py
      --seq_table    <loc_chrom_with_names_isoforms_with_seq.tsv>
      --scansite_tsv <pre-processed scansite.tsv  or  NO_FILE>
      --outdir       <output directory>
      [--use_api]    # fall back to live API when --scansite_tsv is NO_FILE
      [--stringency  High|Medium|Low]   default: High
      [--workers     N]                 default: 10  (API mode only)

Output:
  scansite.tsv  (Protein_ID, motifName, motifShortName, score, site,
                 siteSequence, Start, End)
"""

import argparse
import logging
import re
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_COLS = ["Protein_ID", "motifName", "motifShortName", "score",
         "site", "siteSequence", "Start", "End"]


def _filter_precomputed(src: Path, protein_ids: set, out: Path):
    df = pd.read_csv(src, sep="\t", dtype=str)
    df = df[df["Protein_ID"].isin(protein_ids)].copy()
    keep = [c for c in _COLS if c in df.columns]
    df[keep].to_csv(out, sep="\t", index=False)
    log.info("ScanSite (pre-computed): %d rows for %d proteins",
             len(df), df["Protein_ID"].nunique() if len(df) else 0)


def _call_api(acc, sequence, protein_id, stringency):
    try:
        import requests
        url = (f"https://scansite4.mit.edu/webservice/proteinscan/"
               f"identifier={acc}/sequence={sequence}/"
               f"motifclass=MAMMALIAN/stringency={stringency}")
        r = requests.get(url, timeout=60)
        if r.status_code != 200:
            return []
        root = ET.fromstring(r.content)
        rows = []
        interest = {"motifName", "motifShortName", "score", "site", "siteSequence"}
        for child in root:
            if child.tag != "predictedSite":
                continue
            row = {}
            site_seq = None
            for node in child:
                if node.tag in interest:
                    if node.tag == "siteSequence":
                        site_seq = node.text.replace("*", "") if node.text else None
                        row[node.tag] = site_seq
                    else:
                        row[node.tag] = node.text
            if site_seq:
                m = re.search(re.escape(site_seq), sequence, re.IGNORECASE)
                row["Start"] = m.start() + 1 if m else None
                row["End"]   = m.end()       if m else None
            row["Protein_ID"] = protein_id
            rows.append(row)
        return rows
    except Exception as exc:
        log.warning("ScanSite API error for %s: %s", acc, exc)
        return []


def _run_api(seq_df, out: Path, stringency: str, max_workers: int):
    all_rows = []
    tasks = [
        (row["Entry_Isoform"].split("-")[0], row["Sequence"], row["Protein_ID"])
        for _, row in seq_df.iterrows()
        if pd.notna(row.get("Sequence")) and pd.notna(row.get("Entry_Isoform"))
    ]
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(_call_api, acc, seq, pid, stringency): pid
                for acc, seq, pid in tasks}
        for fut in as_completed(futs):
            all_rows.extend(fut.result())

    df = pd.DataFrame(all_rows) if all_rows else pd.DataFrame(columns=_COLS)
    keep = [c for c in _COLS if c in df.columns]
    df[keep].to_csv(out, sep="\t", index=False)
    log.info("ScanSite (API): %d rows for %d proteins",
             len(df), df["Protein_ID"].nunique() if len(df) else 0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table",    required=True)
    p.add_argument("--scansite_tsv", required=True)
    p.add_argument("--outdir",       required=True)
    p.add_argument("--use_api",      action="store_true")
    p.add_argument("--stringency",   default="High",
                   choices=["High", "Medium", "Low"])
    p.add_argument("--workers",      type=int, default=10)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "scansite.tsv"

    src = Path(args.scansite_tsv)
    use_precomputed = src.exists() and src.stat().st_size > 0 and src.name != "NO_FILE"

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)
    protein_ids = set(seq_df["Protein_ID"].dropna())

    if args.use_api:
        log.info("ScanSite: using live API (stringency=%s, workers=%d)",
                 args.stringency, args.workers)
        _run_api(seq_df, out, args.stringency, args.workers)
    elif use_precomputed:
        _filter_precomputed(src, protein_ids, out)
    else:
        log.info("ScanSite: no pre-computed file and --use_api not set — writing empty output")
        pd.DataFrame(columns=_COLS).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
