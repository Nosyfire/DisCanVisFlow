#!/usr/bin/env python3
"""Normalise raw LLPS-database downloads into three flat TSVs.

All database-specific format knowledge lives here so the mapping workers
(`create_llps_regions_worker.py`, `create_llps_variants_worker.py`) stay generic
and testable. Reads whichever raw sources are provided (each optional) and emits:

  <outdir>/llps_regions.tsv    acc  start  end  region_type  source
  <outdir>/llps_proteins.tsv   acc  gene  organelle_mlo  material_state  klass  source
  <outdir>/llps_variants.tsv   acc  pos  wt_aa  mut_aa  variant_class  disease_term  source

Sources
-------
PhaSepDB  (phasepdbv2_1_llps.xlsx): regions from the ``region`` column, protein
          localisation/MLO/material-state, and clean HGVS ``mutation`` entries
          (``p.V460P``) as *experimental* variants.
LLPSDB    (protein.xls, actually xlsx): IDR/LCR ranges as regions plus
          localisation. Its construct ``Mutation`` column uses opaque internal
          codes (``M14``) with no substitution, so it yields no variants.
DisPhaseDB: a disease-variant table (UniProt/ClinVar/COSMIC/DisGeNET merge).
          Parsed generically by column name when a dump is supplied; the public
          server is frequently offline, so this source is optional.

Empty / missing sources are silently skipped. Always exits 0 with header rows
present, so a run with no LLPS data still produces well-formed (empty) tables.
"""
import argparse
import logging
import re
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("parse_llps")

REGION_COLS  = ["acc", "start", "end", "region_type", "source"]
PROTEIN_COLS = ["acc", "gene", "organelle_mlo", "material_state", "klass", "source"]
VARIANT_COLS = ["acc", "pos", "wt_aa", "mut_aa", "variant_class", "disease_term", "source"]

_RANGE_RE = re.compile(r"(\d+)\s*-\s*(\d+)")
_HGVS_RE  = re.compile(r"p\.?\(?([A-Z])(\d+)([A-Z])\)?")
_EMPTY    = {"", "_", "-", "nan", "none", "na"}


def _clean(v) -> str:
    s = str(v).strip()
    return "" if s.lower() in _EMPTY else s


def _acc_base(v) -> str:
    """UniProt accession without isoform suffix; first token of a list."""
    s = _clean(v)
    if not s:
        return ""
    s = re.split(r"[;,|\s]", s)[0]
    return s.split("-")[0].strip()


def _ranges(cell) -> list[tuple[int, int]]:
    out = []
    for a, b in _RANGE_RE.findall(_clean(cell)):
        a, b = int(a), int(b)
        if 0 < a <= b:
            out.append((a, b))
    return out


def _hgvs(cell) -> list[tuple[str, int, str]]:
    return [(wt, int(pos), mut) for wt, pos, mut in _HGVS_RE.findall(_clean(cell))]


# --------------------------------------------------------------------------- #
# PhaSepDB
# --------------------------------------------------------------------------- #
def parse_phasepdb(path, regions, proteins, variants):
    df = pd.read_excel(path, engine="openpyxl")
    log.info("PhaSepDB: %d rows", len(df))
    for _, r in df.iterrows():
        acc = _acc_base(r.get("uniprot_entry"))
        if not acc:
            continue
        for a, b in _ranges(r.get("region")):
            regions.append([acc, a, b, "PS_region", "phasepdb"])
        proteins.append([acc, _clean(r.get("gene_name")), _clean(r.get("MLO")),
                         _clean(r.get("material_state")), _clean(r.get("class")),
                         "phasepdb"])
        for wt, pos, mut in _hgvs(r.get("mutation")):
            variants.append([acc, pos, wt, mut, "experimental", "", "phasepdb"])


# --------------------------------------------------------------------------- #
# LLPSDB  (protein.xls / LLPS.xls are really xlsx)
# --------------------------------------------------------------------------- #
def parse_llpsdb(protein_path, regions, proteins, variants):
    df = pd.read_excel(protein_path, engine="openpyxl")
    log.info("LLPSDB: %d proteins", len(df))
    for _, r in df.iterrows():
        acc = _acc_base(r.get("Uniprot ID"))
        if not acc:
            continue
        for a, b in _ranges(r.get("IDR")):
            regions.append([acc, a, b, "IDR", "llpsdb"])
        for a, b in _ranges(r.get("LCR")):
            regions.append([acc, a, b, "LCR", "llpsdb"])
        proteins.append([acc, _clean(r.get("Gene name")), _clean(r.get("Localization")),
                         "", "", "llpsdb"])
    # LLPS.xls Mutation column is opaque internal codes (M14) → no usable variants.


# --------------------------------------------------------------------------- #
# DisPhaseDB  (generic column-name parsing; server often offline)
# --------------------------------------------------------------------------- #
def parse_disphasedb(path, regions, proteins, variants):
    p = Path(path)
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(p, sep=sep, dtype=str, low_memory=False)
    cols = {c.lower(): c for c in df.columns}

    def col(*names):
        for n in names:
            if n in cols:
                return cols[n]
        return None

    c_acc = col("uniprot", "uniprot_acc", "accession", "acc", "protein")
    c_var = col("variant", "mutation", "protein_variant", "hgvs")
    c_pos = col("position", "pos", "protein_position")
    c_wt, c_mut = col("wt", "wt_aa", "ref"), col("mut", "mut_aa", "alt")
    c_dis = col("disease", "disease_term", "phenotype", "condition")
    if not c_acc:
        log.warning("DisPhaseDB: no accession column found in %s — skipping", cols)
        return
    log.info("DisPhaseDB: %d rows", len(df))
    for _, r in df.iterrows():
        acc = _acc_base(r.get(c_acc))
        if not acc:
            continue
        disease = _clean(r.get(c_dis)) if c_dis else ""
        parsed = _hgvs(r.get(c_var)) if c_var else []
        if parsed:
            for wt, pos, mut in parsed:
                variants.append([acc, pos, wt, mut, "disease", disease, "disphasedb"])
        elif c_pos and c_wt and c_mut:
            pos = _clean(r.get(c_pos))
            wt, mut = _clean(r.get(c_wt)), _clean(r.get(c_mut))
            if pos.isdigit() and len(wt) == 1 and len(mut) == 1:
                variants.append([acc, int(pos), wt, mut, "disease", disease, "disphasedb"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phasepdb", help="phasepdbv2_1_llps.xlsx")
    ap.add_argument("--llpsdb_protein", help="LLPSDB protein.xls (xlsx)")
    ap.add_argument("--disphasedb", help="DisPhaseDB variant table (csv/tsv)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    regions, proteins, variants = [], [], []
    if args.phasepdb and Path(args.phasepdb).exists():
        parse_phasepdb(args.phasepdb, regions, proteins, variants)
    if args.llpsdb_protein and Path(args.llpsdb_protein).exists():
        parse_llpsdb(args.llpsdb_protein, regions, proteins, variants)
    if args.disphasedb and Path(args.disphasedb).exists():
        parse_disphasedb(args.disphasedb, regions, proteins, variants)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(regions, columns=REGION_COLS).to_csv(
        outdir / "llps_regions.tsv", sep="\t", index=False)
    pd.DataFrame(proteins, columns=PROTEIN_COLS).drop_duplicates().to_csv(
        outdir / "llps_proteins.tsv", sep="\t", index=False)
    pd.DataFrame(variants, columns=VARIANT_COLS).drop_duplicates().to_csv(
        outdir / "llps_variants.tsv", sep="\t", index=False)
    log.info("Wrote regions=%d proteins=%d variants=%d → %s",
             len(regions), len(pd.DataFrame(proteins).drop_duplicates()),
             len(pd.DataFrame(variants).drop_duplicates()), outdir)


if __name__ == "__main__":
    main()
