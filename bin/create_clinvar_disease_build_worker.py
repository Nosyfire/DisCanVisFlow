#!/usr/bin/env python3
"""
Module 8a — Build ClinVar disease ontology table from mapped mutation TSVs + MONDO OBO.

Replaces filter-only mode when MUTATION_MAP outputs are available. Derives
Final_Category from MONDO hierarchy (obonet) and IDP finalize rules instead of
lookup from a static clinvar_diseases.tsv.

Usage:
  create_clinvar_disease_build_worker.py
      --seq_table     <loc_chrom_with_names_isoforms_with_seq.tsv>
      --mondo_obo     <mondo.obo>
      --mutation_dir  <dir with *_filter_mutations_mapped.tsv>
      --outdir        <output directory>

Output:
  clinvar_disease.tsv
  clinvar_disease_mutations.tsv  — simplified mutation ↔ disease merge
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from clinvar_disease_lib import (
    categorize_mondo_id,
    extract_mondo_ids,
    finalize_disease_row,
    load_mondo_graph,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_OUT_COLS = [
    "Protein_ID", "Disease", "DOID", "DO Subset", "synonyms",
    "PhenotypeList", "PhenotypeIDS", "xref_source", "reference_id",
    "level1", "level2", "level3", "level4", "level5", "level6",
    "level7", "level8", "level9", "level10", "level11", "level12", "level13",
    "Disordered", "Ordered", "Total Mutations", "Disordered Percent",
    "Final_Category",
]

_MUTATION_GLOB = "*_filter_mutations_mapped.tsv"


def _load_mutations(mutation_dir: Path) -> pd.DataFrame:
    frames = []
    for p in sorted(mutation_dir.glob(_MUTATION_GLOB)):
        try:
            df = pd.read_csv(p, sep="\t", dtype=str)
        except Exception as exc:
            log.warning("Skip %s: %s", p.name, exc)
            continue
        if df.empty:
            continue
        if "Study Abbrevation" in df.columns:
            df = df[df["Study Abbrevation"].astype(str).str.strip().eq("ClinVar")]
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_disease_table(mut_df: pd.DataFrame, graph) -> pd.DataFrame:
    """Aggregate mutation rows to per-Protein_ID disease associations."""
    if mut_df.empty:
        return pd.DataFrame(columns=_OUT_COLS)

    needed = ["Protein_ID", "PhenotypeList", "PhenotypeIDS"]
    for c in needed:
        if c not in mut_df.columns:
            log.warning("Missing column %s in mutation table", c)
            return pd.DataFrame(columns=_OUT_COLS)

    mut_df = mut_df.fillna("")
    rows = []
    grouped = mut_df.groupby(["Protein_ID", "PhenotypeList", "PhenotypeIDS"], dropna=False)

    for (pid, pheno_list, pheno_ids), grp in grouped:
        diseases = [d.strip() for d in str(pheno_list).split("|") if d.strip()]
        if not diseases:
            diseases = [str(pheno_list).strip()] if str(pheno_list).strip() else ["Unknown"]

        mondo_ids = extract_mondo_ids(pheno_ids)
        mondo_id = mondo_ids[0] if mondo_ids else ""
        fc = categorize_mondo_id(graph, mondo_id) if mondo_id else "Other"

        for disease in diseases:
            row = finalize_disease_row({
                "Protein_ID": pid,
                "Disease": disease.replace("_", " "),
                "DOID": mondo_id.replace("MONDO:", "DOID:") if mondo_id else "",
                "PhenotypeList": pheno_list,
                "PhenotypeIDS": pheno_ids,
                "disease_group": disease.replace("_", " "),
                "Final_Category": fc,
                "Total Mutations": str(len(grp)),
                "Disordered": "0",
                "Ordered": str(len(grp)),
                "Disordered Percent": "0.0",
            })
            rows.append(row)

    if not rows:
        return pd.DataFrame(columns=_OUT_COLS)

    out = pd.DataFrame(rows)
    for col in _OUT_COLS:
        if col not in out.columns:
            out[col] = ""
    return out[_OUT_COLS].drop_duplicates()


def build_disease_mutation_table(mut_df: pd.DataFrame, graph) -> pd.DataFrame:
    """Per-mutation rows with simplified disease name + Final_Category."""
    if mut_df.empty:
        return pd.DataFrame(columns=[
            "Protein_ID", "Mutation", "Protein_position", "Disease",
            "Final_Category", "synonyms", "ClinicalSignificance",
        ])

    mut_df = mut_df.fillna("")
    rows = []
    for _, row in mut_df.iterrows():
        pheno_list = str(row.get("PhenotypeList", "")).strip()
        pheno_ids = str(row.get("PhenotypeIDS", "")).strip()
        diseases = [d.strip().replace("_", " ") for d in pheno_list.split("|") if d.strip()]
        if not diseases:
            diseases = [pheno_list.replace("_", " ")] if pheno_list else ["Unknown"]

        mondo_ids = extract_mondo_ids(pheno_ids)
        mondo_id = mondo_ids[0] if mondo_ids else ""
        fc = categorize_mondo_id(graph, mondo_id) if mondo_id else "Other"

        for disease in diseases:
            rows.append({
                "Protein_ID": row.get("Protein_ID", ""),
                "Mutation": row.get("Mutation", ""),
                "Protein_position": row.get("Protein_position", ""),
                "Disease": disease,
                "Final_Category": fc,
                "synonyms": "",
                "ClinicalSignificance": row.get("ClinicalSignificance", row.get("Mutation Description", "")),
            })

    out = pd.DataFrame(rows).drop_duplicates()
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table", required=True)
    p.add_argument("--mondo_obo", required=True)
    p.add_argument("--mutation_dir", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "clinvar_disease.tsv"
    mut_out_path = outdir / "clinvar_disease_mutations.tsv"
    empty = pd.DataFrame(columns=["Protein_ID", "Disease", "DOID", "Final_Category"])
    empty_mut = pd.DataFrame(columns=[
        "Protein_ID", "Mutation", "Protein_position", "Disease", "Final_Category",
    ])

    mondo = Path(args.mondo_obo)
    mut_dir = Path(args.mutation_dir)
    if not mondo.exists() or mondo.stat().st_size == 0:
        log.info("MONDO OBO not found — writing empty output")
        empty.to_csv(out_path, sep="\t", index=False)
        empty_mut.to_csv(mut_out_path, sep="\t", index=False)
        return
    if not mut_dir.is_dir():
        log.info("Mutation dir missing — writing empty output")
        empty.to_csv(out_path, sep="\t", index=False)
        empty_mut.to_csv(mut_out_path, sep="\t", index=False)
        return

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
    protein_ids = set(seq_df["Protein_ID"].dropna())

    log.info("Loading MONDO OBO …")
    graph = load_mondo_graph(str(mondo))

    log.info("Loading mutation tables from %s", mut_dir)
    mut_df = _load_mutations(mut_dir)
    log.info("Mutation rows loaded: %d", len(mut_df))

    df = build_disease_table(mut_df, graph)
    if not df.empty:
        df = df[df["Protein_ID"].isin(protein_ids)]

    df.to_csv(out_path, sep="\t", index=False)
    log.info("ClinVar disease (built): %d rows for %d proteins",
             len(df), df["Protein_ID"].nunique() if len(df) else 0)

    mut_merge = build_disease_mutation_table(mut_df, graph)
    if not mut_merge.empty:
        mut_merge = mut_merge[mut_merge["Protein_ID"].isin(protein_ids)]
    mut_merge.to_csv(mut_out_path, sep="\t", index=False)
    log.info("ClinVar disease mutations (merged): %d rows",
             len(mut_merge))


if __name__ == "__main__":
    main()
