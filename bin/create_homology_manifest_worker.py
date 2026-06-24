#!/usr/bin/env python3
"""
create_homology_manifest_worker.py — Homology-similarity mapping manifest.

Scans mapped annotation TSVs and records every row that was transferred from a
main isoform onto an alternative isoform by sequence homology (mapping_type ==
'homology_similarity', i.e. legacy homology_transfer == True). The result is a
single audit table so it is traceable which annotations on alternative isoforms
are homology-derived rather than directly observed.

Output
------
  homology_similarity_manifest.tsv with columns:
    annotation  Protein_ID  source_accession  identifier  start  end
    position  mapping_type
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

OUT_COLS = ["annotation", "Protein_ID", "source_accession", "identifier",
            "start", "end", "position", "homology_identity", "mapping_type"]

# Candidate columns (first match wins) used to populate the generic manifest fields
_ACC_COLS = ["Accession", "Entry_Isoform", "uniprot", "Primary_Acc", "acc"]
_ID_COLS = ["ELMIdentifier", "Identifier", "Name", "name", "feature",
            "ID A", "Switch ID", "prot_expr", "Mutation", "Mutation Description"]
_START_COLS = ["Start", "start", "Bindingsite A Start"]
_END_COLS = ["End", "end", "Bindingsite A End"]
_POS_COLS = ["Position", "Protein_position", "position"]


def _first(row, cols):
    for c in cols:
        if c in row and pd.notna(row[c]) and str(row[c]) != "":
            return row[c]
    return ""


def _is_homology(row):
    mt = str(row.get("mapping_type", "")).strip().lower()
    if mt:
        return mt == "homology_similarity"
    ht = str(row.get("homology_transfer", "")).strip().lower()
    return ht in ("true", "1", "yes")


def main():
    ap = argparse.ArgumentParser(description="Build homology-similarity manifest")
    ap.add_argument("--inputs", nargs="*", default=[],
                    help="mapped annotation TSVs to scan")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "homology_similarity_manifest.tsv"

    records = []
    for path in args.inputs:
        p = Path(path)
        if not p.exists() or p.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(p, sep="\t", dtype=str)
        except Exception as exc:  # noqa: BLE001
            log.warning("Could not read %s: %s", p.name, exc)
            continue
        if df.empty or ("mapping_type" not in df.columns
                        and "homology_transfer" not in df.columns):
            continue
        annotation = p.stem
        for _, row in df.iterrows():
            if not _is_homology(row):
                continue
            records.append({
                "annotation": annotation,
                "Protein_ID": row.get("Protein_ID", ""),
                "source_accession": _first(row, _ACC_COLS),
                "identifier": _first(row, _ID_COLS),
                "start": _first(row, _START_COLS),
                "end": _first(row, _END_COLS),
                "position": _first(row, _POS_COLS),
                "homology_identity": row.get("homology_identity", ""),
                "mapping_type": "homology_similarity",
            })

    df_out = pd.DataFrame(records, columns=OUT_COLS)
    df_out.to_csv(out_path, sep="\t", index=False)
    log.info("homology_similarity_manifest.tsv: %d transferred rows across %d files",
             len(df_out), len(args.inputs))


if __name__ == "__main__":
    main()
