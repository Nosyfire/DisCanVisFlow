#!/usr/bin/env python3
"""Map normalised LLPS regions + protein-level info onto every curated isoform.

Consumes the flat tables emitted by ``parse_llps_sources.py`` (source-agnostic)
plus the run's isoform sequence table, and writes two Protein_ID-keyed outputs:

  <outdir>/llps_regions.tsv   Protein_ID  Gene  start  end  region_type  source  mapping_type
  <outdir>/llps_proteins.tsv  Protein_ID  Gene  organelle_mlo  material_state  klass  source  mapping_type

Mapping is by UniProt accession: a region/protein annotated on accession A is
attached to every isoform in the proteome whose UniProt accession (isoform suffix
stripped) is A. Regions are residue-validated against each isoform's own sequence
length (start..end must fit) so coordinates that fall off a shorter isoform are
dropped rather than silently mis-placed. ``mapping_type`` is ``direct`` — every
placement here is within the same UniProt protein.

Empty / missing inputs yield header-only outputs and exit 0.
"""
import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("llps_regions")

REGION_OUT  = ["Protein_ID", "Gene", "start", "end", "region_type", "source", "mapping_type"]
PROTEIN_OUT = ["Protein_ID", "Gene", "organelle_mlo", "material_state", "klass",
               "source", "mapping_type"]
# normalised-input headers (must match parse_llps_sources.py, kept local to avoid
# a cross-module import that Nextflow's per-script staging would not satisfy)
NORM_REGION_COLS  = ["acc", "start", "end", "region_type", "source"]
NORM_PROTEIN_COLS = ["acc", "gene", "organelle_mlo", "material_state", "klass", "source"]


def _load_loc(path):
    df = pd.read_csv(path, sep="\t", dtype=str,
                     usecols=lambda c: c in {"Protein_ID", "Entry_Isoform", "Gene", "Sequence"})
    df["acc"] = df["Entry_Isoform"].astype(str).str.split("-").str[0]
    df["seqlen"] = df["Sequence"].astype(str).str.len()
    return df


def _read_norm(path, cols):
    if not path or not Path(path).exists():
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return df if not df.empty else pd.DataFrame(columns=cols)


def map_regions(regions, loc, outpath):
    if regions.empty:
        pd.DataFrame(columns=REGION_OUT).to_csv(outpath, sep="\t", index=False)
        return 0
    regions = regions.copy()
    regions["start"] = pd.to_numeric(regions["start"], errors="coerce")
    regions["end"] = pd.to_numeric(regions["end"], errors="coerce")
    regions = regions.dropna(subset=["start", "end"])
    regions[["start", "end"]] = regions[["start", "end"]].astype(int)

    merged = regions.merge(loc[["acc", "Protein_ID", "Gene", "seqlen"]], on="acc", how="inner")
    # residue-validate: region must fit inside this isoform's sequence
    merged = merged[(merged["start"] >= 1) & (merged["end"] <= merged["seqlen"])]
    merged["mapping_type"] = "direct"
    out = merged[REGION_OUT].drop_duplicates().sort_values(["Protein_ID", "start", "end"])
    out.to_csv(outpath, sep="\t", index=False)
    log.info("regions: %d placements across %d isoforms → %s",
             len(out), out["Protein_ID"].nunique(), outpath.name)
    return len(out)


def map_proteins(proteins, loc, outpath):
    if proteins.empty:
        pd.DataFrame(columns=PROTEIN_OUT).to_csv(outpath, sep="\t", index=False)
        return 0
    merged = proteins.merge(loc[["acc", "Protein_ID", "Gene"]], on="acc",
                            how="inner", suffixes=("_src", ""))
    merged["mapping_type"] = "direct"
    keep = [c for c in PROTEIN_OUT if c in merged.columns]
    out = merged[keep].drop_duplicates().sort_values(["Protein_ID", "source"])
    out.to_csv(outpath, sep="\t", index=False)
    log.info("proteins: %d rows across %d isoforms → %s",
             len(out), out["Protein_ID"].nunique(), outpath.name)
    return len(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", help="normalised llps_regions.tsv")
    ap.add_argument("--proteins", help="normalised llps_proteins.tsv")
    ap.add_argument("--loc_chrom", required=True, help="isoform sequence table")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    loc = _load_loc(args.loc_chrom)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    map_regions(_read_norm(args.regions, NORM_REGION_COLS), loc, outdir / "llps_regions.tsv")
    map_proteins(_read_norm(args.proteins, NORM_PROTEIN_COLS), loc, outdir / "llps_proteins.tsv")


if __name__ == "__main__":
    main()
