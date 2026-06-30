#!/usr/bin/env python3
"""
create_coiledcoils_worker.py — Module 5i: Coiled-Coil Prediction (DeepCoil)

Runs DeepCoil on protein sequences to predict coiled-coil regions.
DeepCoil requires its own conda environment (discanvis_deepcoil, TF 2.x + PyTorch).
This script dispatches ALL predictions in ONE subprocess call to that environment
(loading the model once), writing per-protein progress to stderr in real time.

Inputs
------
--loc_chrom        loc_chrom_with_names_isoforms_with_seq.tsv
--deepcoil_python  Python binary with DeepCoil installed
--threshold        Per-position CC probability threshold (default: 0.5)
--batch_size       Sequences per DeepCoil internal batch (default: 32)
--n_cpu            CPU threads for DeepCoil (default: -1 = all). Set equal to
                   task.cpus when running with maxForks to avoid oversubscription.
--output_dir       output directory (default: .)

Outputs
-------
coiled_coils.tsv  — Protein_ID | Prob_scores  (DeepCoil per-residue probabilities)
"""

import argparse
import json
import logging
import os
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

DEFAULT_DEEPCOIL_PYTHON = "python"


# ---------------------------------------------------------------------------
# Single-subprocess dispatcher — loads model once for all sequences
# ---------------------------------------------------------------------------

# The subprocess: loads DeepCoil once, processes all sequences in batches,
# writes "PROGRESS done/total batch b/n_batches" to stderr after each batch,
# prints JSON result to stdout when done.
_DEEPCOIL_SCRIPT = r"""
import sys, json, os, warnings
warnings.filterwarnings('ignore')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import torch
n_cpu = int(os.environ.get('DEEPCOIL_N_CPU', '-1'))
if n_cpu > 0:
    torch.set_num_threads(n_cpu)

from deepcoil import DeepCoil

data      = json.loads(sys.stdin.read())
batch_sz  = int(os.environ.get('DEEPCOIL_BATCH_SIZE', '32'))
pids      = list(data.keys())
total     = len(pids)
n_batches = (total + batch_sz - 1) // batch_sz

print(f"DEEPCOIL_INIT n_cpu={n_cpu} total={total} batches={n_batches}", file=sys.stderr, flush=True)

dc = DeepCoil(use_gpu=False, n_cpu=n_cpu)

results = {}
for i in range(0, total, batch_sz):
    chunk = {pid: data[pid] for pid in pids[i:i+batch_sz]}
    b_num = i // batch_sz + 1
    pred  = dc.predict(chunk)
    results.update({pid: scores['cc'].tolist() for pid, scores in pred.items()})
    done = min(i + batch_sz, total)
    print(f"DEEPCOIL_PROGRESS {done}/{total} batch {b_num}/{n_batches}", file=sys.stderr, flush=True)

print(json.dumps(results))
"""


def _run_deepcoil_all(
    proteins: dict,
    python_bin: str,
    n_cpu: int = -1,
    batch_size: int = 32,
    timeout: int = 86400,
) -> dict:
    """
    Run DeepCoil on ALL proteins in one subprocess (model loaded once).
    Progress lines appear in real time in .command.err (stderr flows through).
    """
    env = os.environ.copy()
    env["DEEPCOIL_N_CPU"]       = str(n_cpu)
    env["DEEPCOIL_BATCH_SIZE"]  = str(batch_size)

    log.info("Launching DeepCoil subprocess for %d sequences (n_cpu=%d, batch=%d) …",
             len(proteins), n_cpu, batch_size)
    try:
        result = subprocess.run(
            [python_bin, "-c", _DEEPCOIL_SCRIPT],
            input=json.dumps(proteins),
            stdout=subprocess.PIPE,   # capture JSON result
            stderr=None,              # let progress lines flow to .command.err in real time
            text=True,
            timeout=timeout,
            env=env,
        )
        if result.returncode != 0:
            log.error("DeepCoil subprocess exited with code %d", result.returncode)
            return {}
        return json.loads(result.stdout.strip())
    except subprocess.TimeoutExpired:
        log.error("DeepCoil timed out after %ds — partial results lost", timeout)
        return {}
    except json.JSONDecodeError as e:
        log.error("DeepCoil JSON parse error: %s", e)
        return {}
    except Exception as e:
        log.error("DeepCoil error: %s", e)
        return {}


# ---------------------------------------------------------------------------
# Region extractor
# ---------------------------------------------------------------------------

def scores_to_regions(scores: list, threshold: float = 0.5) -> list:
    regions = []
    in_region = False
    start = 0
    for i, s in enumerate(scores):
        if s > threshold:
            if not in_region:
                in_region = True
                start = i + 1
        else:
            if in_region:
                regions.append((start, i))
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
    p.add_argument("--n_cpu",           type=int,   default=-1,
                   help="CPU threads for DeepCoil. Set to task.cpus to avoid oversubscription.")
    p.add_argument("--output_dir",      default=".")
    args = p.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading loc_chrom …")
    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)

    proteins: dict = {}
    for _, row in loc_df.iterrows():
        pid = str(row.get("Protein_ID", "") or "")
        seq = str(row.get("Sequence", "") or "")
        if pid and pid not in ("nan", "") and seq and seq not in ("nan", ""):
            proteins[pid] = seq

    log.info("Running DeepCoil on %d sequences …", len(proteins))

    all_scores = _run_deepcoil_all(
        proteins,
        args.deepcoil_python,
        n_cpu=args.n_cpu,
        batch_size=args.batch_size,
    )

    log.info("DeepCoil completed: %d/%d proteins scored", len(all_scores), len(proteins))

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
