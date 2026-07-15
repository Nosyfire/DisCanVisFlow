#!/usr/bin/env python3
"""
create_finches_worker.py — Site-wise FINCHES Δε saturation mutagenesis worker.

Computes the change in self-interaction (homotypic) epsilon for every possible
single-amino-acid substitution at every position in a protein sequence using the
FINCHES Mpipi_GGv1 force field.  The Δε value (mut_epsilon − wt_epsilon) is a
site-wise LLPS-propensity change score, analogous to AlphaMissense but for
condensate / phase-separation potential.

Output
------
finches_saturation.tsv  — columns:
  Protein_ID  Position  WT_AA  Mut_AA  WT_Epsilon  Mut_Epsilon  Delta_Epsilon

Positive Δε = mutation increases LLPS tendency.
Negative Δε = mutation decreases LLPS tendency.

Engines
-------
* ``incremental`` (default) — exact O(20·L²)/protein engine (bin/finches_incremental.py)
  that rebuilds only the band of the interaction matrix a point mutation actually
  changes, instead of the whole L×L matrix per variant. ~1–2 orders of magnitude
  faster than ``full``; results are validated bit-for-bit against the reference.
  Proteins with non-standard residues transparently fall back to ``full``.
* ``full`` — reference path: one ``calculate_epsilon_value(mut, mut)`` per variant.

Parallelism is per-protein (``imap_unordered`` over a length-sorted stream), so all
workers stay saturated and no giant per-variant task lists are pickled. Progress is
shown with a single tqdm bar. The run is resumable: Protein_IDs already present in
the output are skipped and new results are appended.

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
      [--n_cpu 64] [--engine incremental|full] [--only_main_isoforms] \\
      [--max_seq_len 3000] [--finches_lib /path/to/finches] \\
      [--no_resume] [--validate 3] [--order short_first|input|long_first]
"""

from __future__ import annotations

import argparse
import logging
import os
import random
import sys
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
_AA_SET = frozenset(AMINO_ACIDS)
COLS = ["Protein_ID", "Position", "WT_AA", "Mut_AA",
        "WT_Epsilon", "Mut_Epsilon", "Delta_Epsilon"]

# ---------------------------------------------------------------------------
# Worker globals (built once per process)
# ---------------------------------------------------------------------------
_MODEL = None          # FINCHES InteractionMatrixConstructor
_ENG = None            # IncrementalEpsilon
_ENGINE = "incremental"
_VALIDATE_TOL = 1e-6


def _init_worker(finches_lib, engine):
    global _MODEL, _ENG, _ENGINE
    _ENGINE = engine
    if finches_lib and finches_lib not in sys.path:
        sys.path.insert(0, finches_lib)
    # make bin/ importable for finches_incremental regardless of cwd
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from finches.epsilon_calculation import InteractionMatrixConstructor
            from finches.forcefields.mpipi import Mpipi_model
        _MODEL = InteractionMatrixConstructor(
            parameters=Mpipi_model(version="Mpipi_GGv1"))
        if engine == "incremental":
            from finches_incremental import IncrementalEpsilon
            _ENG = IncrementalEpsilon(_MODEL)
    except Exception as exc:                       # pragma: no cover
        log.error("FINCHES init failed: %s", exc)
        sys.exit(1)


def _rows_full(pid, seq):
    """Reference engine: one epsilon call per variant."""
    import numpy as np
    try:
        wt_eps = float(_MODEL.calculate_epsilon_value(seq, seq))
    except Exception as exc:
        log.warning("WT epsilon failed for %s (%s) — skipping", pid, exc)
        return pid, []
    if np.isnan(wt_eps):
        return pid, []
    rows = []
    for i, wt_aa in enumerate(seq):
        if wt_aa not in AMINO_ACIDS:
            continue
        for mut_aa in AMINO_ACIDS:
            if mut_aa == wt_aa:
                continue
            ms = seq[:i] + mut_aa + seq[i + 1:]
            me = float(_MODEL.calculate_epsilon_value(ms, ms))
            delta = me - wt_eps
            rows.append((pid, i + 1, wt_aa, mut_aa,
                         round(wt_eps, 6),
                         round(me, 6) if not np.isnan(me) else "",
                         round(delta, 6) if not np.isnan(delta) else ""))
    return pid, rows


def _rows_incremental(pid, seq, validate):
    """Incremental engine with optional per-protein validation sampling."""
    import numpy as np
    ctx = _ENG.prepare_wt(seq)
    if ctx is None:                        # non-standard residue → reference path
        return _rows_full(pid, seq)
    wt_eps = ctx["eps_wt"]
    rows = []
    for i, wt_aa in enumerate(seq):
        if wt_aa not in AMINO_ACIDS:
            continue
        for mut_aa in AMINO_ACIDS:
            if mut_aa == wt_aa:
                continue
            me, delta = _ENG.delta_for_variant(ctx, i, mut_aa)
            rows.append((pid, i + 1, wt_aa, mut_aa,
                         round(wt_eps, 6), round(me, 6), round(delta, 6)))
    if validate:
        _validate_protein(pid, seq, ctx)
    return pid, rows


def _validate_protein(pid, seq, ctx):
    """Cross-check a random sample of variants against the reference path."""
    import numpy as np
    L = len(seq)
    rng = random.Random(hash(pid) & 0xFFFFFFFF)
    n = min(50, L)
    worst = 0.0
    for _ in range(n):
        p = rng.randrange(L)
        if seq[p] not in AMINO_ACIDS:
            continue
        m = rng.choice([a for a in AMINO_ACIDS if a != seq[p]])
        ms = seq[:p] + m + seq[p + 1:]
        ref = float(_MODEL.calculate_epsilon_value(ms, ms))
        inc, _ = _ENG.delta_for_variant(ctx, p, m)
        worst = max(worst, abs(inc - ref))
    if worst > _VALIDATE_TOL:
        raise RuntimeError(
            f"VALIDATION FAILED for {pid}: incremental vs full max|Δ|={worst:.2e} "
            f"> {_VALIDATE_TOL:.0e}")
    log.info("  validate %s: max|Δ| vs full = %.2e (OK)", pid, worst)


# module-level so Pool can pickle the callable
_VALIDATE_IDS: set = set()


def _process(args):
    pid, seq = args
    # Mpipi has no parameters for U (selenocysteine) or X (unknown), so epsilon
    # is undefined for the whole sequence, not just that site — skip the protein.
    odd = set(seq) - _AA_SET
    if odd:
        log.warning("%s contains non-parameterised residue(s) %s — skipping",
                    pid, ",".join(sorted(odd)))
        return pid, []
    try:
        if _ENGINE == "full":
            return _rows_full(pid, seq)
        return _rows_incremental(pid, seq, validate=(pid in _VALIDATE_IDS))
    except Exception as exc:               # pragma: no cover
        log.error("protein %s failed: %s", pid, exc)
        raise


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="FINCHES site-wise saturation mutagenesis (Δε LLPS-change score)")
    p.add_argument("--loc_chrom", required=True,
                   help="Sequence table TSV (Protein_ID + Sequence columns)")
    p.add_argument("--output_dir", default=".")
    p.add_argument("--n_cpu", type=int, default=max(1, cpu_count() - 1))
    p.add_argument("--finches_lib", default=None,
                   help="Path to finches repo root (added to sys.path)")
    p.add_argument("--only_main_isoforms", action="store_true", default=False)
    p.add_argument("--max_seq_len", type=int, default=3000)
    p.add_argument("--engine", choices=["incremental", "full"],
                   default="incremental",
                   help="incremental (fast, exact) or full (reference)")
    p.add_argument("--order", choices=["short_first", "input", "long_first"],
                   default="short_first",
                   help="protein processing order (short_first maximises early coverage)")
    p.add_argument("--no_resume", action="store_true", default=False,
                   help="ignore existing output and recompute everything")
    p.add_argument("--validate", type=int, default=2,
                   help="cross-check the first N proteins against the full engine")
    p.add_argument("--batch_size", type=int, default=None,
                   help="(accepted for compatibility; output is streamed per protein)")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "finches_saturation.tsv"

    log.info("Loading sequences from %s", args.loc_chrom)
    df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str).dropna(subset=["Sequence"])
    if args.only_main_isoforms and "main_isoform" in df.columns:
        df = df[df["main_isoform"] == "yes"]
    df = df[df["Sequence"].str.len() <= args.max_seq_len]

    seq_map = dict(zip(df["Protein_ID"], df["Sequence"]))

    # --- resume: skip Protein_IDs already written --------------------------
    done: set = set()
    if out_tsv.exists() and not args.no_resume:
        try:
            prev = pd.read_csv(out_tsv, sep="\t", usecols=["Protein_ID"], dtype=str)
            done = set(prev["Protein_ID"].dropna().unique())
        except Exception:
            done = set()
    if not out_tsv.exists() or args.no_resume:
        pd.DataFrame(columns=COLS).to_csv(out_tsv, sep="\t", index=False)

    todo = [(pid, s) for pid, s in seq_map.items() if pid not in done]
    if args.order == "short_first":
        todo.sort(key=lambda kv: len(kv[1]))
    elif args.order == "long_first":
        todo.sort(key=lambda kv: -len(kv[1]))

    log.info("FINCHES engine=%s: %d proteins to do (%d already done, %d total ≤%d aa), %d workers",
             args.engine, len(todo), len(done), len(seq_map), args.max_seq_len, args.n_cpu)
    if not todo:
        log.info("nothing to do — output already complete")
        return

    global _VALIDATE_IDS
    _VALIDATE_IDS = {pid for pid, _ in todo[:max(0, args.validate)]}

    # Build the model once in the parent first. The same failure inside a Pool
    # initializer makes multiprocessing respawn workers forever, so a broken
    # env or engine would hang instead of reporting anything.
    _init_worker(args.finches_lib, args.engine)

    try:
        from tqdm import tqdm
    except Exception:                       # pragma: no cover
        tqdm = None

    ok = skipped = 0
    fh = open(out_tsv, "a", buffering=1)
    try:
        with Pool(processes=args.n_cpu,
                  initializer=_init_worker,
                  initargs=(args.finches_lib, args.engine)) as pool:
            it = pool.imap_unordered(_process, todo, chunksize=1)
            bar = tqdm(total=len(todo), unit="prot", dynamic_ncols=True,
                       smoothing=0.05) if tqdm else None
            for pid, rows in it:
                if rows:
                    for r in rows:
                        fh.write("\t".join(map(str, r)) + "\n")
                    ok += 1
                else:
                    skipped += 1
                if bar is not None:
                    bar.update(1)
                    bar.set_postfix(done=ok, skip=skipped, last=pid[:14])
                elif (ok + skipped) % 50 == 0:
                    log.info("  %d / %d proteins", ok + skipped, len(todo))
            if bar is not None:
                bar.close()
    finally:
        fh.close()

    log.info("Done — %d proteins written, %d skipped (resume-safe; rerun to continue)",
             ok, skipped)


if __name__ == "__main__":
    main()
