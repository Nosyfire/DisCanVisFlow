#!/usr/bin/env python3
"""create_plaac_worker.py — PLAAC prion-like domain HMM (Lancaster et al. 2014).

Runs the PLAAC jar in per-residue mode over each sequence and emits the
per-residue PLAAC LLR score and the Viterbi in-PrD flag. Missing jar/Java →
empty output, exit 0 (never crashes a run).

Output
------
plaac.tsv — columns: Protein_ID  Position  plaac_score  in_PRD
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

OUT_COLS = ["Protein_ID", "Position", "plaac_score", "in_PRD"]
DEFAULT_JAR = "External_Programs/plaac/plaac.jar"


def _write_fasta(seq_map: dict, path: Path):
    with open(path, "w") as fh:
        for pid, seq in seq_map.items():
            fh.write(f">{pid}\n{seq}\n")


def _parse_plaac_perres(text: str) -> list[dict]:
    """Parse PLAAC `-p all` per-residue output (16 tab-separated columns)."""
    rows: list[dict] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        toks = line.split("\t")
        if toks[0] == "ORDER":          # header line
            continue
        if len(toks) < 16:
            continue
        try:
            pos = int(toks[2])
        except ValueError:
            continue
        pid = toks[1]
        raw_score = toks[9]
        score = "" if raw_score in ("NaN", "nan", "") else raw_score
        in_prd = 1 if toks[4] == "1" else 0
        rows.append({"Protein_ID": pid, "Position": pos,
                     "plaac_score": score, "in_PRD": in_prd})
    return rows


def main():
    ap = argparse.ArgumentParser(description="PLAAC prion-like domain HMM")
    ap.add_argument("--seq_table", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--only_main_isoforms", action="store_true", default=False)
    ap.add_argument("--plaac_jar", default=DEFAULT_JAR)
    ap.add_argument("--java", default="java")
    ap.add_argument("--bg_freqs", default=None,
                    help="Background AA freqs (-B). Defaults to bg_freqs_HUMAN.txt beside the jar.")
    ap.add_argument("--alpha", default="0", help="PLAAC -a (default 0 = use -B bg fully)")
    ap.add_argument("--core_length", default="60", help="PLAAC -c core PrD length")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "plaac.tsv"

    jar = Path(args.plaac_jar)
    if not jar.exists() or shutil.which(args.java) is None:
        log.warning("PLAAC jar or Java missing — skipping PLAAC track (empty output).")
        pd.DataFrame(columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
        sys.exit(0)

    bg = args.bg_freqs
    if bg is None:
        sib = jar.parent / "bg_freqs_HUMAN.txt"
        bg = str(sib) if sib.exists() else None

    df = pd.read_csv(args.seq_table, sep="\t", dtype=str).dropna(subset=["Sequence"])
    if args.only_main_isoforms and "main_isoform" in df.columns:
        df = df[df["main_isoform"] == "yes"]
    seq_map = dict(zip(df["Protein_ID"], df["Sequence"]))
    log.info("PLAAC: %d sequences", len(seq_map))

    with tempfile.TemporaryDirectory() as td:
        fasta = Path(td) / "in.fa"
        _write_fasta(seq_map, fasta)
        cmd = [args.java, "-jar", str(jar), "-i", str(fasta),
               "-a", str(args.alpha), "-c", str(args.core_length), "-p", "all"]
        if bg:
            cmd += ["-B", bg]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            log.error("PLAAC failed (rc=%s): %s", proc.returncode, proc.stderr[:300])
            sys.exit(1)
        rows = _parse_plaac_perres(proc.stdout)

    pd.DataFrame(rows, columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
    log.info("Done — %d residues over %d proteins",
             len(rows), len({r["Protein_ID"] for r in rows}))


if __name__ == "__main__":
    main()
