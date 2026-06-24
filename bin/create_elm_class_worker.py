#!/usr/bin/env python3
"""
Module 5n — ElmProteomeClassMatch lookup table.

Parses ELM class definitions (elm_classes-*.tsv) and produces a flat TSV for
the ElmProteomeClassMatch Django model. This is a per-pipeline-run lookup
table (not protein-specific) that stores ELM regex patterns for client-side
WT vs mutant SLiM scanning in DisCanVis2.

Usage:
    create_elm_class_worker.py
        --elm_classes  <elm_classes-YYYY.tsv>
        --outdir       <output directory>

Output:
    elm_classes.tsv  — one row per ELM class
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_OUT_COLS = [
    "elm_accession", "elm_identifier", "functional_site_name",
    "description", "regex", "probability",
    "n_instances", "n_instances_in_pdb", "elm_type",
]

_COL_MAP = {
    "Accession":          "elm_accession",
    "ELMIdentifier":      "elm_identifier",
    "FunctionalSiteName": "functional_site_name",
    "Description":        "description",
    "Regex":              "regex",
    "Probability":        "probability",
    "#Instances":         "n_instances",
    "#Instances_in_PDB":  "n_instances_in_pdb",
}


def _elm_type(identifier: str) -> str:
    """Extract ELM type prefix: CLV_... → CLV, DEG_... → DEG, etc."""
    m = re.match(r"^([A-Z]+)_", str(identifier))
    return m.group(1) if m else ""


def parse_elm_classes(path: Path) -> pd.DataFrame:
    if not path.exists():
        log.warning("ELM classes file not found: %s — returning empty table", path)
        return pd.DataFrame(columns=_OUT_COLS)

    # File has comment lines starting with '#' above the header
    rows = []
    header = None
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            parts = [p.strip('"') for p in parts]
            if header is None:
                header = parts
                continue
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))

    if not rows:
        return pd.DataFrame(columns=_OUT_COLS)

    df = pd.DataFrame(rows)
    df = df.rename(columns=_COL_MAP)
    df["elm_type"] = df["elm_identifier"].apply(_elm_type)

    for col in ["probability"]:
        df[col] = pd.to_numeric(df.get(col, pd.Series(dtype=float)), errors="coerce")
    for col in ["n_instances", "n_instances_in_pdb"]:
        df[col] = pd.to_numeric(df.get(col, pd.Series(dtype=float)), errors="coerce").astype("Int64")

    # Keep only mapped columns (some source files may lack extras)
    keep = [c for c in _OUT_COLS if c in df.columns]
    df = df[keep].drop_duplicates(subset=["elm_accession"])

    # Ensure all output columns exist
    for c in _OUT_COLS:
        if c not in df.columns:
            df[c] = None
    return df[_OUT_COLS]


def main():
    p = argparse.ArgumentParser(
        description="Module 5n: ELM class definitions → ElmProteomeClassMatch lookup table"
    )
    p.add_argument("--elm_classes", required=True,
                   help="ELM classes TSV (e.g. legacy_data/elm/elm_classes-2025.tsv)")
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = parse_elm_classes(Path(args.elm_classes))
    df.to_csv(outdir / "elm_classes.tsv", sep="\t", index=False)
    log.info("Done — %d ELM class entries written to elm_classes.tsv", len(df))


if __name__ == "__main__":
    main()
