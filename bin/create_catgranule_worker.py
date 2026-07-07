#!/usr/bin/env python3
"""create_catgranule_worker.py — catGRANULE 2.0 LLPS propensity (Monti et al.).

Per-residue LLPS propensity profile + per-protein RandomForest LLPS score, from
the local catGRANULE 2.0 install. Used as the sequence-based LLPS track (the
specced fallback for PScore). catGRANULE's deps live in a separate env, so this
worker re-invokes itself under --catgranule_python (FINCHES/aiupred pattern),
falling back to direct import. Missing lib/env -> empty output, exit 0.

Output
------
catgranule.tsv — columns: Protein_ID  Position  catgranule_score  catgranule_total
"""
import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

OUT_COLS = ["Protein_ID", "Position", "catgranule_score", "catgranule_total"]
DEFAULT_LIB = "/dlab/home/norbi/PycharmProjects/catGRANULE2.0"


def _load_seqs(seq_table: str, only_main: bool) -> dict:
    df = pd.read_csv(seq_table, sep="\t", dtype=str).dropna(subset=["Sequence"])
    if only_main and "main_isoform" in df.columns:
        df = df[df["main_isoform"] == "yes"]
    seq_map = {}
    for pid, seq in zip(df["Protein_ID"], df["Sequence"]):
        if "U" in seq:                      # catGRANULE cannot score selenocysteine
            continue
        seq_map[pid] = seq
    return seq_map


def _compute_rows(seq_map: dict, lib: str) -> list[dict]:
    """Import catGRANULE from `lib` and compute per-residue rows. Raises on import failure."""
    prev = os.getcwd()
    os.chdir(lib)
    sys.path.insert(0, lib)
    try:
        from compute_profiles_and_predictions import compute_score_and_profile_from_text
        fasta = "".join(f">{pid}\n{seq}\n" for pid, seq in seq_map.items())
        prof, scores = compute_score_and_profile_from_text(fasta)
    finally:
        os.chdir(prev)

    rows = []
    for pid, seq in seq_map.items():
        if pid not in prof or pid not in scores.index:
            log.warning("catGRANULE produced no result for %s — skipping", pid)
            continue
        profile = prof[pid]                                  # numpy ndarray
        total = float(scores.loc[pid, "RandomForest"])       # per-protein LLPS score
        n = min(len(profile), len(seq))
        if len(profile) != len(seq):
            log.warning("profile/seq length mismatch for %s (%d vs %d)",
                        pid, len(profile), len(seq))
        for j in range(n):
            rows.append({"Protein_ID": pid, "Position": j + 1,
                         "catgranule_score": round(float(profile[j]), 6),
                         "catgranule_total": round(total, 6)})
    return rows


def main():
    ap = argparse.ArgumentParser(description="catGRANULE 2.0 LLPS propensity")
    ap.add_argument("--seq_table", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--only_main_isoforms", action="store_true", default=False)
    ap.add_argument("--catgranule_lib", default=DEFAULT_LIB)
    ap.add_argument("--catgranule_python", default=None)
    ap.add_argument("--_inner", action="store_true", default=False,
                    help=argparse.SUPPRESS)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "catgranule.tsv"

    # Route heavy compute into the correct env when a python is given.
    if args.catgranule_python and not args._inner:
        cmd = [args.catgranule_python, os.path.abspath(__file__),
               "--seq_table", args.seq_table,
               "--outdir", str(outdir),
               "--catgranule_lib", args.catgranule_lib,
               "--_inner"]
        if args.only_main_isoforms:
            cmd.append("--only_main_isoforms")
        log.info("Delegating catGRANULE compute to %s", args.catgranule_python)
        sys.exit(subprocess.call(cmd))

    seq_map = _load_seqs(args.seq_table, args.only_main_isoforms)
    log.info("catGRANULE: %d sequences", len(seq_map))

    try:
        rows = _compute_rows(seq_map, args.catgranule_lib)
    except Exception as exc:
        log.warning("catGRANULE unavailable (%s) — skipping track (empty output).", exc)
        pd.DataFrame(columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
        sys.exit(0)

    pd.DataFrame(rows, columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
    log.info("Done — %d residues over %d proteins",
             len(rows), len({r["Protein_ID"] for r in rows}))


if __name__ == "__main__":
    main()
