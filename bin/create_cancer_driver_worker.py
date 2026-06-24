#!/usr/bin/env python3
"""
Module 8c — Cancer Gene Census + Compendium driver gene annotation.

The combined legacy ``cancer_driver.tsv`` records *membership* only
(Protein_ID + 'Cancer Driver' = a comma list of {Census, Compendium}).  The
actual driver *role* lives in two small gene-keyed reference tables vendored
in ``legacy_data/drivers/``:

  census_roles.tsv      Gene | Tier | Role in Cancer | Tumour Types(Somatic) |
                        Tumour Types(Germline)            (Cancer Gene Census)
  compendium_roles.tsv  Gene | ROLE | CANCER_TYPE        (IntOGen Compendium)

Each Protein_ID (e.g. ``RAF1-244``) is reduced to its gene symbol (``RAF1``)
and joined to those tables so the outputs carry the real role, matching the
Django ``DriverGenesCensus`` (role_in_cancer, tumour types) and
``DriverGenesCompendium`` (cancer_type) models.

Usage:
  create_cancer_driver_worker.py
      --seq_table          <loc_chrom_with_names_isoforms_with_seq.tsv>
      --cancer_driver      <combined cancer_driver.tsv (membership)>
      --census_roles       <census_roles.tsv  or NO_FILE>
      --compendium_roles   <compendium_roles.tsv or NO_FILE>
      --outdir             <output directory>

Outputs:
  cancer_driver.tsv     Protein_ID | Cancer Driver | Role in Cancer | Compendium Role
  census_driver.tsv     Protein_ID | Gene | Tier | Role in Cancer |
                        Tumour Types(Somatic) | Tumour Types(Germline)
  compendium_driver.tsv Protein_ID | Gene | ROLE | CANCER_TYPE
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

CENSUS_COLS = ["Protein_ID", "Gene", "Tier", "Role in Cancer",
               "Tumour Types(Somatic)", "Tumour Types(Germline)"]
COMPENDIUM_COLS = ["Protein_ID", "Gene", "ROLE", "CANCER_TYPE"]
COMBINED_COLS = ["Protein_ID", "Cancer Driver", "Role in Cancer", "Compendium Role"]


def _gene_of(protein_id: str) -> str:
    """RAF1-244 -> RAF1  (strip a trailing -<number> transcript suffix)."""
    return re.sub(r"-\d+$", "", str(protein_id))


def _load_roles(path: str | None) -> pd.DataFrame | None:
    if not path or path == "NO_FILE":
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
    if "Gene" not in df.columns:
        log.warning("%s: no 'Gene' column — ignoring", p.name)
        return None
    return df.drop_duplicates(subset=["Gene"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table", required=True)
    p.add_argument("--cancer_driver", required=True,
                   help="Combined cancer_driver.tsv (Protein_ID + 'Cancer Driver')")
    p.add_argument("--census_roles", default=None,
                   help="Gene-keyed CGC role table (census_roles.tsv)")
    p.add_argument("--compendium_roles", default=None,
                   help="Gene-keyed Compendium role table (compendium_roles.tsv)")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    protein_ids = set(
        pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
        ["Protein_ID"].dropna()
    )

    src = Path(args.cancer_driver)
    if not src.exists() or src.stat().st_size == 0 or not protein_ids:
        pd.DataFrame(columns=COMBINED_COLS).to_csv(outdir / "cancer_driver.tsv", sep="\t", index=False)
        pd.DataFrame(columns=CENSUS_COLS).to_csv(outdir / "census_driver.tsv", sep="\t", index=False)
        pd.DataFrame(columns=COMPENDIUM_COLS).to_csv(outdir / "compendium_driver.tsv", sep="\t", index=False)
        log.warning("No combined driver file / proteins — wrote empty outputs")
        return

    df = pd.read_csv(src, sep="\t", dtype=str).fillna("")
    drv_col = "Cancer Driver" if "Cancer Driver" in df.columns else df.columns[1]
    df = df[df["Protein_ID"].isin(protein_ids)].copy()
    df["Gene"] = df["Protein_ID"].map(_gene_of)
    df["__drv"] = df[drv_col].fillna("")

    census_roles = _load_roles(args.census_roles)
    comp_roles = _load_roles(args.compendium_roles)

    # ----- census_driver.tsv -----
    cen = df[df["__drv"].str.contains("Census", case=False)][["Protein_ID", "Gene"]].copy()
    if census_roles is not None:
        cen = cen.merge(census_roles, on="Gene", how="left")
    for c in CENSUS_COLS:
        if c not in cen.columns:
            cen[c] = ""
    cen = cen[CENSUS_COLS].fillna("")
    cen.to_csv(outdir / "census_driver.tsv", sep="\t", index=False)

    # ----- compendium_driver.tsv -----
    comp = df[df["__drv"].str.contains("Compendium", case=False)][["Protein_ID", "Gene"]].copy()
    if comp_roles is not None:
        comp = comp.merge(comp_roles, on="Gene", how="left")
    for c in COMPENDIUM_COLS:
        if c not in comp.columns:
            comp[c] = ""
    comp = comp[COMPENDIUM_COLS].fillna("")
    comp.to_csv(outdir / "compendium_driver.tsv", sep="\t", index=False)

    # ----- combined cancer_driver.tsv (membership + roles) -----
    role_lookup = (census_roles.set_index("Gene")["Role in Cancer"].to_dict()
                   if census_roles is not None else {})
    comp_lookup = (comp_roles.set_index("Gene")["ROLE"].to_dict()
                   if comp_roles is not None else {})
    out = df[["Protein_ID", drv_col, "Gene"]].copy()
    out["Role in Cancer"] = out["Gene"].map(role_lookup).fillna("")
    out["Compendium Role"] = out["Gene"].map(comp_lookup).fillna("")
    out = out.rename(columns={drv_col: "Cancer Driver"})[COMBINED_COLS]
    out.to_csv(outdir / "cancer_driver.tsv", sep="\t", index=False)

    log.info("Cancer drivers — Census: %d rows, Compendium: %d rows, combined: %d rows",
             len(cen), len(comp), len(out))


if __name__ == "__main__":
    main()
