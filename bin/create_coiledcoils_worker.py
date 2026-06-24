#!/usr/bin/env python3
"""
create_coiledcoils_worker.py — Module 5i: Coiled-Coil Prediction (DeepCoil)

Runs DeepCoil on protein sequences to predict coiled-coil regions.
DeepCoil requires its own conda environment (deepcoil_env, TF 1.x + PyTorch).
This script dispatches predictions via subprocess to that environment.

Inputs
------
--loc_chrom        loc_chrom_with_names_isoforms_with_seq.tsv
--deepcoil_python  Python binary with DeepCoil installed
                   (default: /home/nosyfire/miniconda3/envs/deepcoil_env/bin/python)
--threshold        Per-position CC probability threshold (default: 0.5)
--batch_size       Sequences per DeepCoil batch (default: 32)
--output_dir       output directory (default: .)

Outputs
-------
coiled_coils.tsv  — Protein_ID | Prob_scores  (DeepCoil per-residue probabilities)
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

DEFAULT_DEEPCOIL_PYTHON = "python"  # override via --deepcoil_python or create discanvis_deepcoil env


# ---------------------------------------------------------------------------
# Subprocess dispatcher — runs DeepCoil in its own environment
# ---------------------------------------------------------------------------

_DEEPCOIL_SCRIPT = """
import sys, json
import warnings
warnings.filterwarnings('ignore')
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from deepcoil import DeepCoil

data = json.loads(sys.stdin.read())
dc = DeepCoil(use_gpu=False)
results = dc.predict(data)
out = {pid: scores['cc'].tolist() for pid, scores in results.items()}
print(json.dumps(out))
"""


def _run_deepcoil_batch(
    batch: dict,
    python_bin: str,
    timeout: int = 1800,
) -> dict:
    """Run DeepCoil on {Protein_ID: sequence} dict via subprocess."""
    try:
        result = subprocess.run(
            [python_bin, "-c", _DEEPCOIL_SCRIPT],
            input=json.dumps(batch),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode != 0:
            log.debug("DeepCoil stderr: %s", result.stderr[:500])
            return {}
        return json.loads(result.stdout.strip())
    except (subprocess.TimeoutExpired, json.JSONDecodeError, Exception) as e:
        log.warning("DeepCoil batch failed: %s", e)
        return {}


# ---------------------------------------------------------------------------
# Region extractor
# ---------------------------------------------------------------------------

def scores_to_regions(scores: list[float], threshold: float = 0.5) -> list[tuple]:
    """Convert per-position probability list to (start, end) 1-based regions."""
    regions = []
    in_region = False
    start = 0
    for i, s in enumerate(scores):
        if s > threshold:
            if not in_region:
                in_region = True
                start = i + 1  # 1-based
        else:
            if in_region:
                regions.append((start, i))  # i is already end (exclusive → last was i)
                in_region = False
    if in_region:
        regions.append((start, len(scores)))
    return regions


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Module 5i: CoiledCoil prediction via DeepCoil")
    p.add_argument("--loc_chrom",       required=True)
    p.add_argument("--deepcoil_python", default=DEFAULT_DEEPCOIL_PYTHON)
    p.add_argument("--threshold",       type=float, default=0.5)
    p.add_argument("--batch_size",      type=int,   default=32)
    p.add_argument("--output_dir",      default=".")
    args = p.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading loc_chrom …")
    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)

    proteins: dict[str, str] = {}
    for _, row in loc_df.iterrows():
        pid = str(row.get("Protein_ID", "") or "")
        seq = str(row.get("Sequence", "") or "")
        if pid and pid not in ("nan", "") and seq and seq not in ("nan", ""):
            proteins[pid] = seq

    log.info("Running DeepCoil on %d sequences (batch_size=%d) …",
             len(proteins), args.batch_size)

    all_scores: dict[str, list[float]] = {}
    pids = list(proteins.keys())

    for i in range(0, len(pids), args.batch_size):
        batch_pids = pids[i: i + args.batch_size]
        batch = {pid: proteins[pid] for pid in batch_pids}
        log.info("  batch %d/%d …", i // args.batch_size + 1,
                 (len(pids) + args.batch_size - 1) // args.batch_size)
        scores = _run_deepcoil_batch(batch, args.deepcoil_python)
        all_scores.update(scores)

    log.info("DeepCoil completed: %d/%d proteins scored", len(all_scores), len(proteins))

    # Single coiled-coil output: per-protein DeepCoil per-residue probabilities.
    score_rows = [
        {"Protein_ID": pid,
         "Prob_scores": ",".join(str(round(s, 4)) for s in scores)}
        for pid, scores in all_scores.items()
    ]
    out_df = (pd.DataFrame(score_rows)
              if score_rows
              else pd.DataFrame(columns=["Protein_ID", "Prob_scores"]))
    out_df.to_csv(outdir / "coiled_coils.tsv", sep="\t", index=False)
    log.info("CoiledCoils (DeepCoil): %d proteins scored → coiled_coils.tsv", len(out_df))


if __name__ == "__main__":
    main()
