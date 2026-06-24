#!/usr/bin/env python3
"""
create_pem_worker.py — Module 5h: PEM (Predicted ELM Motifs) annotation

Filters the predicted_elm_dataset.tsv (from HotspotPEM / IDP Pathogenic
Mutations supplementary) to proteins present in the current run's loc_chrom.

The input TSV is already keyed by Protein_ID (GENCODE transcript name,
e.g. PEX5-221) — no coordinate remapping needed.

Inputs
------
--loc_chrom     loc_chrom_with_names_isoforms_with_seq.tsv
--pem_dataset   predicted_elm_dataset.tsv
--output_dir    output directory (default: .)

Outputs
-------
pem_core_motifs.tsv  — Protein_ID | ELM_Accession | ELMIdentifier | ELMType |
                        Start | End | InstanceLogic | References | Methods |
                        PDB | Organism
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

OUT_COLS = [
    "Protein_ID",
    "ELM_Accession",
    "ELMIdentifier",
    "ELMType",
    "Start",
    "End",
    "InstanceLogic",
    "References",
    "Methods",
    "PDB",
    "Organism",
    "Found_Known",
]


def main():
    p = argparse.ArgumentParser(description="Module 5h: PEM Core Motifs filter")
    p.add_argument("--loc_chrom",   required=True)
    p.add_argument("--pem_dataset", required=True,
                   help="predicted_elm_dataset.tsv (Protein_ID-keyed)")
    p.add_argument("--output_dir",  default=".")
    args = p.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading loc_chrom …")
    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)
    pid_set = set(loc_df["Protein_ID"].dropna().unique())
    log.info("Proteins in run: %d", len(pid_set))

    log.info("Loading PEM dataset …")
    try:
        pem_df = pd.read_csv(args.pem_dataset, sep="\t", dtype=str, low_memory=False)
    except Exception as e:
        log.warning("Cannot read PEM dataset: %s — writing empty output", e)
        pd.DataFrame(columns=OUT_COLS).to_csv(outdir / "pem_core_motifs.tsv", sep="\t", index=False)
        return

    pid_col = next((c for c in ["Protein_ID", "protein_id"] if c in pem_df.columns), None)
    if pid_col is None:
        log.warning("No Protein_ID column in PEM dataset — writing empty output")
        pd.DataFrame(columns=OUT_COLS).to_csv(outdir / "pem_core_motifs.tsv", sep="\t", index=False)
        return

    if pid_col != "Protein_ID":
        pem_df = pem_df.rename(columns={pid_col: "Protein_ID"})

    # Filter to proteins in this run
    pem_df = pem_df[pem_df["Protein_ID"].isin(pid_set)].copy()

    # Map InstanceLogic from Found_Known flag (legacy server contract)
    if "Found_Known" in pem_df.columns:
        def _logic(v):
            v = str(v).strip().lower()
            if v == "true":
                return "known_ELM"
            if v == "false":
                return "predicted_ELM"
            return "-"
        if "InstanceLogic" not in pem_df.columns:
            pem_df["InstanceLogic"] = pem_df["Found_Known"].apply(_logic)

    # Ensure all output columns exist
    for col in OUT_COLS:
        if col not in pem_df.columns:
            pem_df[col] = "-"

    out = pem_df[OUT_COLS].drop_duplicates()
    out.to_csv(outdir / "pem_core_motifs.tsv", sep="\t", index=False)
    log.info("PEM: %d motif rows for %d proteins", len(out), out["Protein_ID"].nunique())


if __name__ == "__main__":
    main()
