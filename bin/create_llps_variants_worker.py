#!/usr/bin/env python3
"""Map normalised LLPS variants onto every curated isoform (residue-validated).

Consumes the flat variant table from ``parse_llps_sources.py`` plus the run's
isoform sequence table and writes one Protein_ID-keyed output:

  <outdir>/llps_variants.tsv
      Protein_ID  Protein_position  WT_AA  Mut_AA  variant_class
      disease_term  source_db  mapping_type

Mapping is by UniProt accession *and* residue identity: a variant (acc, pos,
WT→Mut) is attached to an isoform only when that isoform's residue at ``pos``
actually equals WT (the ProteinGym convention). This prevents mis-placing a
variant onto an isoform whose coordinates have drifted. ``variant_class`` keeps
DisPhaseDB disease variants separate from PhaSepDB/LLPSDB experimental
construct mutations; ``source_db`` records provenance.

Empty / missing input yields a header-only output and exits 0.
"""
import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("llps_variants")

OUT_COLS = ["Protein_ID", "Protein_position", "WT_AA", "Mut_AA", "variant_class",
            "disease_term", "source_db", "mapping_type"]


def _load_loc(path):
    df = pd.read_csv(path, sep="\t", dtype=str,
                     usecols=lambda c: c in {"Protein_ID", "Entry_Isoform", "Sequence"})
    df["acc"] = df["Entry_Isoform"].astype(str).str.split("-").str[0]
    return df[["acc", "Protein_ID", "Sequence"]]


def map_variants(variants, loc, outpath):
    if variants.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(outpath, sep="\t", index=False)
        return 0, 0
    variants = variants.copy()
    variants["pos"] = pd.to_numeric(variants["pos"], errors="coerce")
    variants = variants.dropna(subset=["pos"])
    variants["pos"] = variants["pos"].astype(int)

    merged = variants.merge(loc, on="acc", how="inner")

    def residue_ok(row):
        seq = row["Sequence"]
        p = row["pos"]
        return isinstance(seq, str) and 1 <= p <= len(seq) and seq[p - 1] == row["wt_aa"]

    n_pairs = len(merged)
    merged = merged[merged.apply(residue_ok, axis=1)]
    out = pd.DataFrame({
        "Protein_ID": merged["Protein_ID"],
        "Protein_position": merged["pos"],
        "WT_AA": merged["wt_aa"],
        "Mut_AA": merged["mut_aa"],
        "variant_class": merged["variant_class"],
        "disease_term": merged["disease_term"],
        "source_db": merged["source"],
        "mapping_type": "direct",
    }).drop_duplicates().sort_values(["Protein_ID", "Protein_position", "Mut_AA"])
    out.to_csv(outpath, sep="\t", index=False)
    log.info("variants: %d placements across %d isoforms (%d acc-pairs, %d dropped on residue check) → %s",
             len(out), out["Protein_ID"].nunique(), n_pairs, n_pairs - len(merged), outpath.name)
    return len(out), out["Protein_ID"].nunique()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", help="normalised llps_variants.tsv")
    ap.add_argument("--loc_chrom", required=True, help="isoform sequence table")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    norm_cols = ["acc", "pos", "wt_aa", "mut_aa", "variant_class", "disease_term", "source"]
    if args.variants and Path(args.variants).exists():
        v = pd.read_csv(args.variants, sep="\t", dtype=str).fillna("")
    else:
        v = pd.DataFrame(columns=norm_cols)
    loc = _load_loc(args.loc_chrom)
    map_variants(v if not v.empty else pd.DataFrame(columns=norm_cols),
                 loc, outdir / "llps_variants.tsv")


if __name__ == "__main__":
    main()
