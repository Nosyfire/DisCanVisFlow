#!/usr/bin/env python3
"""create_lcr_worker.py — low-complexity regions (LCR) via NCBI SEG (segmasker).

Runs `segmasker` (BLAST+) over every protein sequence and emits the 1-based
inclusive masked intervals. Backbone parity with the legacy pipeline's LCR track.

Output
------
low_complexity.tsv — columns: Protein_ID  start  end  length
"""
import argparse
import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

OUT_COLS = ["Protein_ID", "start", "end", "length"]


def _write_fasta(seq_map: dict, path: Path):
    with open(path, "w") as fh:
        for pid, seq in seq_map.items():
            fh.write(f">{pid}\n{seq}\n")


def _parse_segmasker(text: str) -> list[dict]:
    """segmasker -outfmt interval:  header line '>id' then 'START - END' (0-based)."""
    rows: list[dict] = []
    current = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            continue
        if current is None or "-" not in line:
            continue
        a, _, b = line.partition("-")
        try:
            s0 = int(a.strip())
            e0 = int(b.strip())
        except ValueError:
            continue
        start = s0 + 1          # segmasker intervals are 0-based inclusive
        end = e0 + 1
        rows.append({"Protein_ID": current, "start": start, "end": end,
                     "length": end - start + 1})
    return rows


def main():
    ap = argparse.ArgumentParser(description="SEG low-complexity regions via segmasker")
    ap.add_argument("--seq_table", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--only_main_isoforms", action="store_true", default=False)
    ap.add_argument("--segmasker", default="segmasker")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "low_complexity.tsv"

    if shutil.which(args.segmasker) is None:
        log.warning("segmasker not found on PATH — skipping LCR track (empty output).")
        pd.DataFrame(columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
        sys.exit(0)

    df = pd.read_csv(args.seq_table, sep="\t", dtype=str).dropna(subset=["Sequence"])
    if args.only_main_isoforms and "main_isoform" in df.columns:
        df = df[df["main_isoform"] == "yes"]
    seq_map = dict(zip(df["Protein_ID"], df["Sequence"]))
    log.info("LCR: %d sequences", len(seq_map))

    with tempfile.TemporaryDirectory() as td:
        fasta = Path(td) / "in.fasta"
        _write_fasta(seq_map, fasta)
        proc = subprocess.run(
            [args.segmasker, "-in", str(fasta), "-infmt", "fasta",
             "-outfmt", "interval"],
            capture_output=True, text=True)
        if proc.returncode != 0:
            log.error("segmasker failed: %s", proc.stderr)
            sys.exit(1)
        rows = _parse_segmasker(proc.stdout)

    pd.DataFrame(rows, columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
    log.info("Done — %d masked intervals over %d proteins",
             len(rows), len({r["Protein_ID"] for r in rows}))


if __name__ == "__main__":
    main()
