#!/usr/bin/env python3
"""create_apr_worker.py — aggregation-prone regions (APR) from the AGGRESCAN a3v scale.

Slides a window over every protein sequence, averages the per-residue AGGRESCAN a3v
aggregation-propensity value, and emits maximal runs of consecutive positions whose
window mean clears a threshold. Sequence-derived like the LCR track: no external
reference, no external binary, so it has no "missing dependency" degradation path.

Scale: Conchillo-Sole et al. (2007) BMC Bioinformatics 8:65 — a3v intrinsic
aggregation propensity. Positive = aggregation-prone.

Threshold, in precedence order:
  1. --threshold VALUE                       -> used verbatim
  2. >= --min_proteins_for_threshold inputs  -> q(--quantile) of all window means
  3. otherwise                               -> DEFAULT_A3V_THRESHOLD, logged loudly

Output
------
aggregation_prone_regions.tsv — columns:
    Protein_ID  start  end  length  mean_a3v  peak_a3v
`start`/`end` are 1-based inclusive (matching low_complexity.tsv).
`mean_a3v` is the mean RAW a3v of the segment's residues; `peak_a3v` is the maximum
WINDOW mean inside the segment (the quantity the threshold acts on).
"""
import argparse
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

# AGGRESCAN a3v (Conchillo-Sole et al. 2007). Positive = aggregation-prone.
A3V = {
    "A": -0.036, "R": -1.240, "N": -1.302, "D": -1.836, "C": 0.604, "Q": -1.231,
    "E": -1.412, "G": -0.535, "H": -1.033, "I": 1.822, "L": 1.380, "K": -0.931,
    "M": 0.910, "F": 1.754, "P": -1.407, "S": -0.294, "T": -0.159, "W": 1.037,
    "Y": 1.159, "V": 1.594,
}

# PROVISIONAL — replaced by the measured proteome q90 in Task 5 of the plan.
# Only used when the input is too small to calibrate (single-gene runs).
DEFAULT_A3V_THRESHOLD = 0.5
MIN_PROTEINS_FOR_THRESHOLD = 1000

OUT_COLS = ["Protein_ID", "start", "end", "length", "mean_a3v", "peak_a3v"]


def window_means(seq: str, w: int) -> list:
    """Centered, end-truncated window means of the a3v scale, one per residue.

    A window containing any non-standard residue yields NaN, which excludes that
    position from every run (a NaN can never clear the threshold).
    """
    n = len(seq)
    vals = [A3V.get(a, math.nan) for a in seq]
    half = w // 2
    out = []
    for i in range(n):
        win = vals[max(0, i - half):min(n, i + half + 1)]
        if any(math.isnan(v) for v in win):
            out.append(math.nan)
        else:
            out.append(sum(win) / len(win))
    return out


def runs_above(means: list, threshold: float, min_len: int) -> list:
    """Maximal runs of >= min_len consecutive positions with means[i] >= threshold.

    Returns 0-based inclusive (start, end) tuples.
    """
    runs = []
    start = None
    for i, m in enumerate(means):
        ok = (not math.isnan(m)) and (m >= threshold)
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            if i - start >= min_len:
                runs.append((start, i - 1))
            start = None
    if start is not None and len(means) - start >= min_len:
        runs.append((start, len(means) - 1))
    return runs


def resolve_threshold(means_map: dict, args) -> float:
    """Apply the three-way threshold precedence and log which branch was taken."""
    if args.threshold is not None:
        log.info("APR threshold: %.4f (supplied via --threshold)", args.threshold)
        return float(args.threshold)
    if len(means_map) >= args.min_proteins_for_threshold:
        pool = np.fromiter((m for ms in means_map.values() for m in ms), dtype=float)
        thr = float(np.nanquantile(pool, args.quantile))
        log.info("APR threshold: %.4f (q%.0f over %d window means from %d proteins)",
                 thr, args.quantile * 100, pool.size, len(means_map))
        return thr
    log.warning("APR threshold: %.4f (DEFAULT — only %d proteins, need >= %d to "
                "calibrate; pass --threshold to override)",
                DEFAULT_A3V_THRESHOLD, len(means_map), args.min_proteins_for_threshold)
    return DEFAULT_A3V_THRESHOLD


def main():
    ap = argparse.ArgumentParser(
        description="Aggregation-prone regions from the AGGRESCAN a3v scale")
    ap.add_argument("--seq_table", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--only_main_isoforms", action="store_true", default=False)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--min_len", type=int, default=5)
    ap.add_argument("--quantile", type=float, default=0.90)
    ap.add_argument("--threshold", type=float, default=None)
    ap.add_argument("--min_proteins_for_threshold", type=int,
                    default=MIN_PROTEINS_FOR_THRESHOLD)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "aggregation_prone_regions.tsv"

    df = pd.read_csv(args.seq_table, sep="\t", dtype=str).dropna(subset=["Sequence"])
    if args.only_main_isoforms and "main_isoform" in df.columns:
        df = df[df["main_isoform"] == "yes"]
    seq_map = dict(zip(df["Protein_ID"], df["Sequence"]))
    log.info("APR: %d sequences, window=%d min_len=%d",
             len(seq_map), args.window, args.min_len)

    means_map = {pid: window_means(seq, args.window) for pid, seq in seq_map.items()}
    threshold = resolve_threshold(means_map, args)

    rows = []
    for pid, seq in seq_map.items():
        ms = means_map[pid]
        for s0, e0 in runs_above(ms, threshold, args.min_len):
            seg = [A3V[a] for a in seq[s0:e0 + 1]]
            rows.append({
                "Protein_ID": pid,
                "start": s0 + 1,
                "end": e0 + 1,
                "length": e0 - s0 + 1,
                "mean_a3v": round(sum(seg) / len(seg), 4),
                "peak_a3v": round(max(ms[s0:e0 + 1]), 4),
            })

    pd.DataFrame(rows, columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
    log.info("Done — %d APRs over %d proteins (threshold %.4f)",
             len(rows), len({r["Protein_ID"] for r in rows}), threshold)


if __name__ == "__main__":
    main()
