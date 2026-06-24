#!/usr/bin/env python3
"""
Module 5m — PositionBasedAnnotations + RSAscores aggregation.

Aggregates per-residue data from disorder, conservation, and annotation modules
into a single per-position table for the PositionBasedAnnotations Django model,
and derives RSA from pLDDT (rsa = (100 - plddt) / 100).

Usage:
    create_position_based_worker.py
        --seq_table          <loc_chrom_with_names_isoforms_with_seq.tsv>
        --iupred_tsv         <IUPredscores.tsv>
        --plddt_tsv          <AlphaFoldTable.tsv>
        --combined_pos_tsv   <CombinedDisorderNew_Pos.tsv>
        --phastcons_tsv      <conservation_phastcons.tsv>   (or NO_FILE)
        --conservation_tsv   <conservation_multiple_level.tsv>  (or NO_FILE)
        --pfam_tsv           <pfam_domains.tsv>   (or NO_FILE)
        --outdir             <output directory>

Output:
    rsa_scores.tsv                 — Protein_ID | rsascores (comma-sep)
    position_based_annotations.tsv — one row per (Protein_ID, position)
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

# GOPHER level → PositionBasedAnnotations column mapping
_LEVEL_TO_COL = {
    "global":        "conservationGlobal",
    "Mammalia":      "conservationMammal",
    "Vertebrata":    "conservationVertebrate",
    "Eukaryota":     "conservationEukaryota",
    "Eumetazoa":     "conservationEumetazoa",
    "Opisthokonta":  "conservationOpisthokonta",
    "Viridiplantae": "conservationViridiplantae",
}

_CONS_COLS = [
    "conservationGlobal", "conservationMammal", "conservationVertebrate",
    "conservationEukaryota", "conservationEumetazoa", "conservationOpisthokonta",
    "conservationViridiplantae",
]

_OUT_COLS = [
    "Protein_ID", "position",
    "plddt", "rsa", "iupred",
    "edisorder", "combineddisorder",
    "phastCons",
] + _CONS_COLS + ["pfam"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_no_file(path: Path) -> bool:
    return path.name == "NO_FILE" or not path.exists()


def _parse_scores(s: str) -> list[float]:
    if not s or str(s) in ("nan", "", "None"):
        return []
    try:
        return [float(x.strip()) for x in str(s).split(",") if x.strip()]
    except ValueError:
        return []


def _load_score_map(path: Path, id_col: str, score_col: str) -> dict[str, list[float]]:
    """Return {Protein_ID: [float, ...]} from a tab-sep file with a comma-sep score column."""
    if _is_no_file(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    if id_col not in df.columns or score_col not in df.columns:
        return {}
    result = {}
    for _, row in df.iterrows():
        pid = str(row[id_col]).strip()
        if pid and pid != "nan":
            result[pid] = _parse_scores(row.get(score_col, ""))
    return result


def _load_combined_disorder(path: Path) -> dict[str, dict[int, int]]:
    """Return {Protein_ID: {1-based_pos: 0_or_1}}."""
    if _is_no_file(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    pid_col = next((c for c in ["Protein_ID", "Entry_Name"] if c in df.columns), None)
    if pid_col is None or "Position" not in df.columns or "CombinedDisorder" not in df.columns:
        return {}
    result: dict[str, dict[int, int]] = {}
    for _, row in df.iterrows():
        pid = str(row[pid_col]).strip()
        if not pid or pid == "nan":
            continue
        try:
            pos = int(float(row["Position"]))
            val = int(float(row.get("CombinedDisorder", 0)))
        except (ValueError, TypeError):
            continue
        result.setdefault(pid, {})[pos] = val
    return result


def _load_conservation(path: Path) -> dict[str, dict[str, list[float]]]:
    """Return {Protein_ID: {level: [float, ...]}} from conservation_multiple_level.tsv."""
    if _is_no_file(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "Protein_ID" not in df.columns or "level" not in df.columns or "conservationscores" not in df.columns:
        return {}
    result: dict[str, dict[str, list[float]]] = {}
    for _, row in df.iterrows():
        pid = str(row["Protein_ID"]).strip()
        level = str(row["level"]).strip()
        if not pid or pid == "nan" or not level or level == "nan":
            continue
        scores = _parse_scores(row.get("conservationscores", ""))
        result.setdefault(pid, {})[level] = scores
    return result


def _load_phastcons(path: Path) -> dict[str, list[float]]:
    """Return {Protein_ID: [float, ...]} from conservation_phastcons.tsv."""
    if _is_no_file(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "Protein_ID" not in df.columns or "conservationscores" not in df.columns:
        return {}
    result = {}
    for _, row in df.iterrows():
        pid = str(row["Protein_ID"]).strip()
        if pid and pid != "nan":
            result[pid] = _parse_scores(row.get("conservationscores", ""))
    return result


def _load_pfam(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Return {Protein_ID: [(start, end, hmm_name), ...]} from pfam_domains.tsv.
    Uses envelope_start/envelope_end columns (InterPro API output).
    Only includes Domain-type entries.
    """
    if _is_no_file(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    pid_col = next((c for c in ["Protein_ID", "Accession"] if c in df.columns), None)
    if pid_col is None:
        return {}
    start_col = next((c for c in ["envelope_start", "alignment_start", "Start"]
                      if c in df.columns), None)
    end_col   = next((c for c in ["envelope_end",   "alignment_end",   "End"]
                      if c in df.columns), None)
    if start_col is None or end_col is None:
        return {}

    result: dict[str, list[tuple[int, int, str]]] = {}
    for _, row in df.iterrows():
        # Skip non-Domain entries (Repeat, Homologous_superfamily, etc.)
        if "Domain" not in str(row.get("type", "Domain")):
            continue
        pid = str(row[pid_col]).strip()
        if not pid or pid == "nan":
            continue
        try:
            s = int(float(row[start_col]))
            e = int(float(row[end_col]))
        except (ValueError, TypeError):
            continue
        name = str(row.get("hmm_name", row.get("hmm_acc", ""))).strip()
        result.setdefault(pid, []).append((s, e, name))
    return result


# ---------------------------------------------------------------------------
# Core builder
# ---------------------------------------------------------------------------

def build_position_based(
    seq_df:          pd.DataFrame,
    iupred_map:      dict[str, list[float]],
    plddt_map:       dict[str, list[float]],
    combined_dis:    dict[str, dict[int, int]],
    phastcons_map:   dict[str, list[float]],
    conservation_map: dict[str, dict[str, list[float]]],
    pfam_map:        dict[str, list[tuple[int, int, str]]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (position_df, rsa_df).

    position_df has one row per (Protein_ID, position).
    rsa_df has one row per Protein_ID with rsascores as a comma-separated string.
    """
    pid_col = next((c for c in ["Protein_ID", "Entry_Name"] if c in seq_df.columns), None)
    if pid_col is None:
        return pd.DataFrame(columns=_OUT_COLS), pd.DataFrame(columns=["Protein_ID", "rsascores"])

    seen_pids: set[str] = set()
    pos_rows = []
    rsa_rows = []

    for _, row in seq_df.iterrows():
        pid = str(row.get(pid_col, "")).strip()
        seq = str(row.get("Sequence", "")).strip()
        if not pid or pid == "nan" or not seq or seq in ("nan", ""):
            continue
        if pid in seen_pids:
            continue
        seen_pids.add(pid)

        n = len(seq)
        iupred_lst  = iupred_map.get(pid, [])
        plddt_lst   = plddt_map.get(pid, [])
        cdis_map    = combined_dis.get(pid, {})
        pcons_lst   = phastcons_map.get(pid, [])
        cons_levels = conservation_map.get(pid, {})
        pfam_list   = pfam_map.get(pid, [])

        # RSA derived from pLDDT: rsa = (100 - plddt) / 100
        rsa_lst = [(100.0 - p) / 100.0 if not np.isnan(p) else np.nan
                   for p in plddt_lst]

        # RSA scores output (comma-separated, same format as IUPredscores)
        rsa_str = ", ".join(f"{v:.4f}" for v in rsa_lst) if rsa_lst else ""
        if rsa_str:
            rsa_rows.append({"Protein_ID": pid, "rsascores": rsa_str})

        for i in range(n):
            pos   = i + 1
            plddt = plddt_lst[i]  if i < len(plddt_lst) else None
            rsa   = rsa_lst[i]    if i < len(rsa_lst)   else None
            iupred = iupred_lst[i] if i < len(iupred_lst) else None
            cdis  = cdis_map.get(pos, 0)
            pcons = pcons_lst[i]  if i < len(pcons_lst) else None

            # Pfam names covering this position
            pfam_names = [name for (s, e, name) in pfam_list if s <= pos <= e]
            pfam_str   = "|".join(pfam_names) if pfam_names else "-"

            r: dict = {
                "Protein_ID":       pid,
                "position":         pos,
                "plddt":            plddt,
                "rsa":              rsa,
                "iupred":           iupred,
                "edisorder":        bool(cdis),
                "combineddisorder": float(cdis),
                "phastCons":        pcons,
                "pfam":             pfam_str,
            }
            for level, col in _LEVEL_TO_COL.items():
                lvl_lst = cons_levels.get(level, [])
                r[col] = lvl_lst[i] if i < len(lvl_lst) else None

            pos_rows.append(r)

    pos_df = pd.DataFrame(pos_rows, columns=_OUT_COLS) if pos_rows \
        else pd.DataFrame(columns=_OUT_COLS)
    rsa_df = pd.DataFrame(rsa_rows, columns=["Protein_ID", "rsascores"]) if rsa_rows \
        else pd.DataFrame(columns=["Protein_ID", "rsascores"])
    return pos_df, rsa_df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Module 5m: per-residue PositionBasedAnnotations + RSAscores"
    )
    p.add_argument("--seq_table",        required=True)
    p.add_argument("--iupred_tsv",       required=True)
    p.add_argument("--plddt_tsv",        required=True)
    p.add_argument("--combined_pos_tsv", required=True)
    p.add_argument("--phastcons_tsv",    default="NO_FILE")
    p.add_argument("--conservation_tsv", default="NO_FILE")
    p.add_argument("--pfam_tsv",         default="NO_FILE")
    p.add_argument("--outdir",           default=".")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading seq_table …")
    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)

    pid_col = next((c for c in ["Protein_ID", "Entry_Name"] if c in seq_df.columns), "Protein_ID")
    n_pids = seq_df[pid_col].nunique() if pid_col in seq_df.columns else 0
    log.info("  %d unique Protein_IDs", n_pids)

    log.info("Loading disorder scores …")
    iupred_map   = _load_score_map(Path(args.iupred_tsv),       pid_col, "IUPredscores")
    plddt_map    = _load_score_map(Path(args.plddt_tsv),        pid_col, "Plldtscores")
    combined_dis = _load_combined_disorder(Path(args.combined_pos_tsv))
    log.info("  IUPred: %d proteins  pLDDT: %d proteins  CombinedDis: %d proteins",
             len(iupred_map), len(plddt_map), len(combined_dis))

    log.info("Loading conservation …")
    phastcons_map    = _load_phastcons(Path(args.phastcons_tsv))
    conservation_map = _load_conservation(Path(args.conservation_tsv))
    log.info("  phastCons: %d proteins  GOPHER: %d proteins",
             len(phastcons_map), len(conservation_map))

    log.info("Loading Pfam domains …")
    pfam_map = _load_pfam(Path(args.pfam_tsv))
    log.info("  Pfam: %d proteins with domain entries", len(pfam_map))

    log.info("Building per-position table …")
    pos_df, rsa_df = build_position_based(
        seq_df, iupred_map, plddt_map, combined_dis,
        phastcons_map, conservation_map, pfam_map,
    )

    pos_df.to_csv(outdir / "position_based_annotations.tsv", sep="\t", index=False)
    rsa_df.to_csv(outdir / "rsa_scores.tsv",                 sep="\t", index=False)

    log.info("Done — %d position rows across %d proteins; %d RSA rows",
             len(pos_df),
             pos_df["Protein_ID"].nunique() if not pos_df.empty else 0,
             len(rsa_df))


if __name__ == "__main__":
    main()
