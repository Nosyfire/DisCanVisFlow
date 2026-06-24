#!/usr/bin/env python3
"""
Module 5j — Protein-Protein Interactions (PPI).

Merges BioGRID + IntAct + HIPPIE processed interaction tables and filters to
proteins in the current run. Partner proteins (Protein_ID_B) are normalized to
the canonical main isoform (GENE-201); same-gene pairs are excluded.

Usage:
  create_ppi_worker.py
      --seq_table  <loc_chrom_with_names_isoforms_with_seq.tsv>
      --intact     <Interaction_intact.tsv  or NO_FILE>
      --biogrid    <Interaction_biogrid.tsv or NO_FILE>
      --hippie     <Interaction_hippie.tsv  or NO_FILE>
      --outdir     <output directory>

Output:
  interactions.tsv  — Protein_ID_A, Protein_ID_B, database, number_of_pubmed
"""

import argparse
import logging
import re
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

_ACC_A = "Accession A"
_ACC_B = "Accession B"
_PUBS = "Publication Identifiers"


def _count_pubmed(pub_str: str) -> int:
    if pd.isna(pub_str) or not str(pub_str).strip():
        return 0
    ids = re.findall(r"pubmed:(\d+)", str(pub_str), re.IGNORECASE)
    return len(set(ids))


def _gene_from_protein_id(protein_id: str) -> str:
    if not protein_id or "-" not in protein_id:
        return protein_id
    return protein_id.rsplit("-", 1)[0]


def load_interactions(path: Path, db_label: str) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as exc:
        log.warning("Could not read %s: %s", path, exc)
        return None
    if _ACC_A not in df.columns or _ACC_B not in df.columns:
        log.warning("Skipping %s — missing Accession columns", path.name)
        return None
    df = df[[_ACC_A, _ACC_B, _PUBS]].copy()
    df["_db"] = db_label
    return df


def main():
    p = argparse.ArgumentParser(description="Module 5j — PPI")
    p.add_argument("--seq_table", required=True)
    p.add_argument("--intact", required=True)
    p.add_argument("--biogrid", required=True)
    p.add_argument("--hippie", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_cols = ["Protein_ID_A", "Protein_ID_B", "database", "number_of_pubmed"]

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)
    run_pids = set(seq_df["Protein_ID"].dropna())
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"] if c in seq_df.columns), None)
    run_genes = set()
    gene_to_run_pids: dict[str, list[str]] = {}
    # UniProt base accession → list of Protein_IDs in this run (for raw-format input files)
    uniprot_to_pids: dict[str, list[str]] = {}
    for _, row in seq_df.iterrows():
        pid = str(row.get("Protein_ID", ""))
        gene = str(row.get(gene_col, "")) if gene_col else _gene_from_protein_id(pid)
        if not pid or not gene:
            continue
        run_genes.add(gene)
        gene_to_run_pids.setdefault(gene, [])
        if pid not in gene_to_run_pids[gene]:
            gene_to_run_pids[gene].append(pid)
        # Build UniProt→pids map from Entry_Isoform column if present
        eif = str(row.get("Entry_Isoform", ""))
        if eif and eif != "nan":
            base = eif.split("-")[0]
            uniprot_to_pids.setdefault(base, [])
            if pid not in uniprot_to_pids[base]:
                uniprot_to_pids[base].append(pid)

    frames: list[pd.DataFrame] = []
    for arg_name, label in [("intact", "IntAct"), ("biogrid", "BioGRID"), ("hippie", "HIPPIE")]:
        df = load_interactions(Path(getattr(args, arg_name)), label)
        if df is not None:
            frames.append(df)

    if not frames:
        log.info("No interaction databases available — writing empty output")
        pd.DataFrame(columns=out_cols).to_csv(outdir / "interactions.tsv", sep="\t", index=False)
        return

    combined = pd.concat(frames, ignore_index=True)

    rows: list[dict] = []
    seen: set[tuple] = set()

    for _, r in combined.iterrows():
        acc_a = str(r[_ACC_A]).strip()
        acc_b = str(r[_ACC_B]).strip()
        if not acc_a or not acc_b:
            continue

        gene_a = _gene_from_protein_id(acc_a)
        gene_b = _gene_from_protein_id(acc_b)
        if gene_a == gene_b:
            continue

        partner_a = acc_a
        partner_b = acc_b
        main_b = f"{gene_b}-201"

        # Support both GENE-201 format and raw UniProt accessions (from preprocessed raw files)
        pids_a = (gene_to_run_pids.get(gene_a) or
                  (uniprot_to_pids.get(acc_a) if acc_a in uniprot_to_pids else None))
        pids_b = (gene_to_run_pids.get(gene_b) or
                  (uniprot_to_pids.get(acc_b) if acc_b in uniprot_to_pids else None))

        a_in_run = bool(pids_a) or acc_a in run_pids
        b_in_run = bool(pids_b) or acc_b in run_pids

        if not a_in_run and not b_in_run:
            continue

        pubmed_count = _count_pubmed(r[_PUBS])
        db_label     = r["_db"]

        a_pid_list = pids_a or ([acc_a] if acc_a in run_pids else [])
        b_pid_list = pids_b or ([acc_b] if acc_b in run_pids else [])

        if a_in_run and b_in_run:
            # Both proteins in this run — emit A→B using their actual Protein_IDs
            for pid_a in a_pid_list:
                for pid_b in b_pid_list:
                    key = (pid_a, pid_b, db_label)
                    if key not in seen:
                        seen.add(key)
                        rows.append({"Protein_ID_A": pid_a, "Protein_ID_B": pid_b,
                                     "database": db_label, "number_of_pubmed": pubmed_count})
        elif a_in_run:
            # A is ours; B is an external partner (normalized to gene-201)
            partner = f"{gene_b}-201"
            for pid_a in a_pid_list:
                key = (pid_a, partner, db_label)
                if key not in seen:
                    seen.add(key)
                    rows.append({"Protein_ID_A": pid_a, "Protein_ID_B": partner,
                                 "database": db_label, "number_of_pubmed": pubmed_count})
        else:
            # B is ours; A is an external partner (normalized to gene-201)
            partner = f"{gene_a}-201"
            for pid_b in b_pid_list:
                key = (partner, pid_b, db_label)
                if key not in seen:
                    seen.add(key)
                    rows.append({"Protein_ID_A": partner, "Protein_ID_B": pid_b,
                                 "database": db_label, "number_of_pubmed": pubmed_count})

    if not rows:
        log.info("No interactions found for proteins in this run")
        pd.DataFrame(columns=out_cols).to_csv(outdir / "interactions.tsv", sep="\t", index=False)
        return

    out_df = pd.DataFrame(rows)
    agg = (
        out_df.groupby(["Protein_ID_A", "Protein_ID_B"], sort=False)
        .agg(
            database=("database", lambda s: "|".join(sorted(set(s.dropna())))),
            number_of_pubmed=("number_of_pubmed", "sum"),
        )
        .reset_index()
    )
    agg.to_csv(outdir / "interactions.tsv", sep="\t", index=False)
    log.info(
        "PPI: %d unique interactions for %d run isoforms",
        len(agg),
        agg["Protein_ID_A"].str.extract(r"^([^-]+)", expand=False).nunique()
        + agg["Protein_ID_B"].str.extract(r"^([^-]+)", expand=False).nunique(),
    )


if __name__ == "__main__":
    main()
