#!/usr/bin/env python3
"""
create_finches_worker.py — Site-wise FINCHES Δε saturation mutagenesis worker.

Computes the change in self-interaction (homotypic) epsilon for every possible
single-amino-acid substitution at every position in a protein sequence using the
FINCHES Mpipi_GGv1 force field.  The Δε value (mut_epsilon − wt_epsilon) serves
as a site-wise LLPS-propensity change score, analogous to AlphaMissense but for
condensate / phase-separation potential.

Output
------
finches_saturation.tsv  — columns:
  Protein_ID  Position  WT_AA  Mut_AA  WT_Epsilon  Mut_Epsilon  Delta_Epsilon

Positive Δε = mutation increases LLPS tendency.
Negative Δε = mutation decreases LLPS tendency.

Citation
--------
Ginell GM, Emenecker RJ, Lotthammer JM, Usher ET & Holehouse AS (2024)
FINCHES. bioRxiv 2024.06.03.597104. https://doi.org/10.1101/2024.06.03.597104

License: CC BY-NC 4.0 — non-commercial use only.

Usage
-----
  create_finches_worker.py \\
      --loc_chrom  loc_chrom_with_names.tsv \\
      --output_dir . \\
      [--n_cpu 4] \\
      [--finches_lib /path/to/finches-main] \\
      [--only_main_isoforms] \\
      [--max_seq_len 3000]
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

_WORKER_MODEL = None


# ---------------------------------------------------------------------------
# Pool-level functions (top-level so pickle works)
# ---------------------------------------------------------------------------

def _init_worker(finches_lib: str | None):
    global _WORKER_MODEL
    if finches_lib and finches_lib not in sys.path:
        sys.path.insert(0, finches_lib)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from finches.epsilon_calculation import InteractionMatrixConstructor
            from finches.forcefields.mpipi import Mpipi_model
        _WORKER_MODEL = InteractionMatrixConstructor(
            parameters=Mpipi_model(version='Mpipi_GGv1'))
    except Exception as exc:
        log.error("FINCHES model init failed: %s", exc)
        sys.exit(1)


def _compute_wt_epsilon(seq: str) -> float:
    """Return self-interaction epsilon for the wild-type sequence."""
    return float(_WORKER_MODEL.calculate_epsilon_value(seq, seq))


def _compute_mut_epsilon(args: tuple) -> tuple:
    """Return (pos, wt_aa, mut_aa, mut_epsilon) for one variant."""
    mut_seq, pos, wt_aa, mut_aa = args
    try:
        eps = float(_WORKER_MODEL.calculate_epsilon_value(mut_seq, mut_seq))
    except Exception:
        eps = float("nan")
    return (pos, wt_aa, mut_aa, eps)


# ---------------------------------------------------------------------------
# Per-protein saturation scan
# ---------------------------------------------------------------------------

def _saturation_scan(pool: Pool, protein_id: str, sequence: str) -> list[dict]:
    wt_epsilon = pool.apply(_compute_wt_epsilon, (sequence,))
    if wt_epsilon is None or np.isnan(wt_epsilon):
        log.warning("WT epsilon failed for %s — skipping", protein_id)
        return []

    tasks = []
    for i, wt_aa in enumerate(sequence):
        if wt_aa not in AMINO_ACIDS:
            continue
        pos    = i + 1
        prefix = sequence[:i]
        suffix = sequence[i + 1:]
        for mut_aa in AMINO_ACIDS:
            if mut_aa != wt_aa:
                tasks.append((prefix + mut_aa + suffix, pos, wt_aa, mut_aa))

    rows = []
    for pos, wt_aa, mut_aa, mut_eps in pool.imap(
            _compute_mut_epsilon, tasks, chunksize=50):
        delta = (mut_eps - wt_epsilon) if not np.isnan(mut_eps) else float("nan")
        rows.append({
            "Protein_ID":    protein_id,
            "Position":      pos,
            "WT_AA":         wt_aa,
            "Mut_AA":        mut_aa,
            "WT_Epsilon":    round(wt_epsilon, 6),
            "Mut_Epsilon":   round(mut_eps, 6)  if not np.isnan(mut_eps)  else "",
            "Delta_Epsilon": round(delta, 6)    if not np.isnan(delta)    else "",
        })
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="FINCHES site-wise saturation mutagenesis (Δε LLPS-change score)")
    p.add_argument("--loc_chrom",          required=True,
                   help="Sequence table TSV (Protein_ID + Sequence columns)")
    p.add_argument("--output_dir",         default=".")
    p.add_argument("--n_cpu",              type=int,
                   default=max(1, cpu_count() - 1))
    p.add_argument("--finches_lib",        default=None,
                   help="Path to finches repo root (added to sys.path)")
    p.add_argument("--only_main_isoforms", action="store_true", default=False)
    p.add_argument("--max_seq_len",        type=int, default=3000)
    p.add_argument("--batch_size",         type=int, default=100,
                   help="Flush output every N proteins")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "finches_saturation.tsv"

    log.info("Loading sequences from %s", args.loc_chrom)
    df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str).dropna(subset=["Sequence"])

    if args.only_main_isoforms and "main_isoform" in df.columns:
        df = df[df["main_isoform"] == "yes"]

    df = df[df["Sequence"].str.len() <= args.max_seq_len]
    seq_map: dict[str, str] = dict(zip(df["Protein_ID"], df["Sequence"]))
    log.info("FINCHES: %d proteins (≤%d aa), %d workers",
             len(seq_map), args.max_seq_len, args.n_cpu)

    cols = ["Protein_ID", "Position", "WT_AA", "Mut_AA",
            "WT_Epsilon", "Mut_Epsilon", "Delta_Epsilon"]
    pd.DataFrame(columns=cols).to_csv(out_tsv, sep="\t", index=False)

    buffer: list[dict] = []
    ok = skipped = 0

    with Pool(processes=args.n_cpu,
              initializer=_init_worker,
              initargs=(args.finches_lib,)) as pool:

        for idx, (pid, seq) in enumerate(seq_map.items()):
            if (idx + 1) % 50 == 0:
                log.info("  %d / %d processed", idx + 1, len(seq_map))

            rows = _saturation_scan(pool, pid, seq)
            if rows:
                buffer.extend(rows)
                ok += 1
            else:
                skipped += 1

            if ok > 0 and ok % args.batch_size == 0:
                pd.DataFrame(buffer, columns=cols).to_csv(
                    out_tsv, mode="a", header=False, sep="\t", index=False)
                buffer = []

        if buffer:
            pd.DataFrame(buffer, columns=cols).to_csv(
                out_tsv, mode="a", header=False, sep="\t", index=False)

    n_rows = sum(1 for _ in open(out_tsv)) - 1
    log.info("Done — %d proteins, %d rows written, %d skipped", ok, n_rows, skipped)


if __name__ == "__main__":
    main()
