#!/usr/bin/env python3
"""
Module 7 — Conservation scores.

Outputs:
  conservation_multiple_level.tsv   — GOPHER trident scores per taxonomic level
  conservation_phastcons.tsv        — phastCons per-residue scores

Usage:
  create_conservation_worker.py
      --seq_table              <loc_chrom_with_names_isoforms_with_seq.tsv>
      --conservation_table     <GOPHER conservation_table.tsv>
      --combined_map           <combined_map.map>
      --outdir                 <output directory>
      [--phastcons_bedgraph    <pre-converted BedGraph file>]
      [--phastcons_dir         <directory of per-chrom .bw files>]
      [--bigwigtobedgraph      <path to bigWigToBedGraph binary>]
      [--skip_gopher]          skip GOPHER conservation
      [--skip_phastcons]       skip phastCons conservation
"""

import argparse
import logging
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# combined_map.map parsing
# ---------------------------------------------------------------------------

def parse_combined_map(map_path: Path) -> dict[str, dict]:
    """Return {Protein_ID: {'chrom': str, 'strand': str, 'residues': [(gpos1,gpos2,gpos3), ...]}}"""
    result: dict[str, dict] = {}
    current_pid: str | None = None
    current_chrom = ""
    current_strand = ""
    current_residues: list = []

    def _flush():
        if current_pid:
            result[current_pid] = {
                "chrom": current_chrom,
                "strand": current_strand,
                "residues": list(current_residues),
            }

    with open(map_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                _flush()
                current_residues = []
                parts = line.split()
                if len(parts) < 4:
                    current_pid = None
                    continue
                # Extract Protein_ID: field index 4 (0-based) in |-split of parts[1]
                fasta_id = parts[1]
                pid_parts = fasta_id.split("|")
                current_pid = pid_parts[4] if len(pid_parts) > 4 else None
                current_chrom = parts[2] if len(parts) > 2 else ""
                current_strand = parts[3] if len(parts) > 3 else "+"
            else:
                if current_pid is None:
                    continue
                cols = line.split()
                if len(cols) < 6:
                    continue
                gpos_field = cols[5]  # e.g. "12618720,12618719,12618718,"
                gpos_parts = [g.strip() for g in gpos_field.split(",") if g.strip()]
                g1 = gpos_parts[0] if len(gpos_parts) > 0 else "-"
                g2 = gpos_parts[1] if len(gpos_parts) > 1 else "-"
                g3 = gpos_parts[2] if len(gpos_parts) > 2 else "-"
                current_residues.append((g1, g2, g3))

    _flush()
    return result


# ---------------------------------------------------------------------------
# phastCons scoring from BedGraph
# ---------------------------------------------------------------------------

def _bedgraph_to_pos_score(bg_path: Path) -> dict[int, float]:
    """Parse a BedGraph file into {genomic_position: score}."""
    pos_score: dict[int, float] = {}
    with open(bg_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
                score = float(parts[3])
            except ValueError:
                continue
            for p in range(start, end):
                pos_score[p] = score
    return pos_score


def _run_bigwig_to_bedgraph(bw_path: Path, chrom: str, beg: int, end: int,
                             bigwig_bin: str) -> Path | None:
    """Run bigWigToBedGraph and return the temp file path, or None on failure."""
    tmp = tempfile.NamedTemporaryFile(suffix=".bedgraph", delete=False)
    tmp.close()
    cmd = [bigwig_bin, f"-chrom={chrom}", f"-start={beg}", f"-end={end}",
           str(bw_path), tmp.name]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return Path(tmp.name)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        log.warning("bigWigToBedGraph failed for %s %s:%d-%d: %s", bw_path.name, chrom, beg, end, exc)
        Path(tmp.name).unlink(missing_ok=True)
        return None


def compute_phastcons_scores(
    pid: str,
    prot_info: dict,
    phastcons_dir: Path | None,
    bigwig_bin: str,
    preloaded_pos_score: dict[int, float] | None = None,
) -> list[float] | None:
    """Return per-residue phastCons score list or None if not available."""
    residues = prot_info["residues"]
    chrom = prot_info["chrom"]

    if chrom == "chrM":
        return None

    # Collect all valid genomic positions
    all_positions: list[int] = []
    for g1, g2, g3 in residues:
        for g in (g1, g2, g3):
            if g != "-":
                try:
                    all_positions.append(int(g))
                except ValueError:
                    pass

    if not all_positions:
        return [0.0] * len(residues)

    # Get position→score map
    if preloaded_pos_score is not None:
        pos_score = preloaded_pos_score
    elif phastcons_dir is not None:
        bw_path = phastcons_dir / f"{chrom}.bw"
        if not bw_path.exists():
            log.warning("phastCons bigWig not found: %s", bw_path)
            return None
        beg = min(all_positions)
        end = max(all_positions) + 1
        bg_file = _run_bigwig_to_bedgraph(bw_path, chrom, beg, end, bigwig_bin)
        if bg_file is None:
            return None
        pos_score = _bedgraph_to_pos_score(bg_file)
        bg_file.unlink(missing_ok=True)
    else:
        return None

    # Map residues to mean phastCons score
    scores: list[float] = []
    for g1, g2, g3 in residues:
        vals: list[float] = []
        for g in (g1, g2, g3):
            if g != "-":
                try:
                    vals.append(pos_score.get(int(g), 0.0))
                except ValueError:
                    pass
        scores.append(sum(vals) / len(vals) if vals else 0.0)

    return scores


# ---------------------------------------------------------------------------
# GOPHER conservation
# ---------------------------------------------------------------------------

def load_gopher_conservation(cons_table_path: Path) -> dict[str, dict[str, str]]:
    """Return {canonical_acc: {level: scores_string}}"""
    result: dict[str, dict[str, str]] = defaultdict(dict)
    if not cons_table_path.exists():
        return result
    df = pd.read_csv(cons_table_path, sep="\t", dtype=str)
    if df.empty or "uniprot_acc" not in df.columns:
        return result
    for _, row in df.iterrows():
        acc = str(row["uniprot_acc"]).strip()
        level = str(row["level"]).strip()
        scores = str(row["conservation_score"]).strip()
        result[acc][level] = scores
    return result


def _canonical_acc(acc: str) -> str:
    """Strip isoform suffix: P04049-2 → P04049."""
    return acc.split("-")[0] if "-" in acc else acc


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Module 7 — Conservation scores")
    p.add_argument("--seq_table", required=True)
    p.add_argument("--conservation_table", required=True)
    p.add_argument("--combined_map", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--phastcons_bedgraph",
                   help="Pre-converted BedGraph file (for testing / single-gene runs)")
    p.add_argument("--phastcons_dir",
                   help="Directory containing per-chromosome .bw files")
    p.add_argument("--bigwigtobedgraph", default="bigWigToBedGraph",
                   help="Path to bigWigToBedGraph binary")
    p.add_argument("--skip_gopher", action="store_true")
    p.add_argument("--skip_phastcons", action="store_true")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)
    # Keep only needed columns; tolerate missing Chromosome
    needed = ["Protein_ID", "Entry_Isoform"]
    if "Chromosome" in seq_df.columns:
        needed.append("Chromosome")
    seq_df = seq_df[needed].drop_duplicates()

    # ── GOPHER multiple-level conservation ───────────────────────────────────
    gopher_rows: list[dict] = []
    if not args.skip_gopher:
        gopher_data = load_gopher_conservation(Path(args.conservation_table))
        for _, row in seq_df.iterrows():
            pid = str(row["Protein_ID"])
            acc = str(row["Entry_Isoform"])
            canon = _canonical_acc(acc)
            if canon not in gopher_data:
                continue
            for level, scores_str in gopher_data[canon].items():
                gopher_rows.append({
                    "Protein_ID": pid,
                    "Entry_Isoform": acc,
                    "level": level,
                    "conservationscores": scores_str,
                })

    gopher_df = pd.DataFrame(gopher_rows) if gopher_rows else pd.DataFrame(
        columns=["Protein_ID", "Entry_Isoform", "level", "conservationscores"])
    gopher_df.to_csv(outdir / "conservation_multiple_level.tsv", sep="\t", index=False)
    log.info("GOPHER conservation: %d rows (Protein_ID×level)", len(gopher_df))

    # ── phastCons conservation ────────────────────────────────────────────────
    phastcons_rows: list[dict] = []
    if not args.skip_phastcons:
        map_path = Path(args.combined_map)
        if map_path.exists() and map_path.stat().st_size > 0:
            prot_map = parse_combined_map(map_path)
        else:
            prot_map = {}

        # Preloaded BedGraph (single-gene / test mode)
        preloaded: dict[int, float] | None = None
        if args.phastcons_bedgraph:
            preloaded = _bedgraph_to_pos_score(Path(args.phastcons_bedgraph))

        phastcons_dir = Path(args.phastcons_dir) if args.phastcons_dir else None

        for _, row in seq_df.iterrows():
            pid = str(row["Protein_ID"])
            acc = str(row["Entry_Isoform"])
            if pid not in prot_map:
                continue
            prot_info = prot_map[pid]
            scores = compute_phastcons_scores(
                pid, prot_info,
                phastcons_dir=phastcons_dir,
                bigwig_bin=args.bigwigtobedgraph,
                preloaded_pos_score=preloaded,
            )
            if scores is None:
                continue
            phastcons_rows.append({
                "Protein_ID": pid,
                "Entry_Isoform": acc,
                "conservationscores": ", ".join(f"{v:.4f}" for v in scores),
            })

    phastcons_df = pd.DataFrame(phastcons_rows) if phastcons_rows else pd.DataFrame(
        columns=["Protein_ID", "Entry_Isoform", "conservationscores"])
    phastcons_df.to_csv(outdir / "conservation_phastcons.tsv", sep="\t", index=False)
    log.info("phastCons: %d proteins", len(phastcons_df))


if __name__ == "__main__":
    main()
