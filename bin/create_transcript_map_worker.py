#!/usr/bin/env python3
"""
create_transcript_map_worker.py — Module 5e: Transcript Annotation Mapping

Maps UniProt-keyed annotations onto every Gencode transcript (Protein_ID,
e.g. "RAF1-201") for the same gene.

Logic
-----
For each annotation (keyed by UniProt Entry_Isoform):
  1. Find source isoform sequence in loc_chrom.
  2. For each transcript (Protein_ID) of the same gene:
     a. If Entry_Isoform matches → copy directly (homology_transfer=False).
     b. Else extract region from source sequence, substring-search in target
        sequence → copy with adjusted Start/End + homology_transfer=True.
  3. For positional annotations (PTM, Position column): 3-aa context window.

Disorder files (CombinedDisorderNew*.tsv) are already keyed by Protein_ID from
DISORDER_MAP, so they are passed through unchanged into mapped/disorder/.

Inputs
------
  --loc_chrom      loc_chrom_with_names_isoforms_with_seq.tsv
  --elm / --dibs / --mfib / --phasepro
  --uniprot_roi / --uniprot_bind
  --ptm
  --disorder       CombinedDisorderNew.tsv (passed through, already Protein_ID-keyed)
  --disorder_pos   CombinedDisorderNew_Pos.tsv (pass-through)
  --only_main_isoforms   (flag; if set, only map to main isoform)
  --output_dir

Outputs (in mapped/ folder, same file names as unmapped/)
---------
  elm.tsv                 — regional annotations mapped to Protein_ID
  dibs.tsv
  mfib.tsv
  phasepro.tsv
  uniprot_roi.tsv
  uniprot_binding.tsv
  ptm_merged.tsv
  CombinedDisorderNew.tsv     — pass-through (already Protein_ID-keyed)
  CombinedDisorderNew_Pos.tsv — pass-through
  transcript_map_stats.tsv
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Load loc_chrom and build lookup structures
# ---------------------------------------------------------------------------

def load_loc(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["Sequence"] = df["Sequence"].fillna("")
    return df


def build_lookup(loc_df: pd.DataFrame, only_main: bool = False):
    """
    Returns:
      acc_to_seq       : Entry_Isoform → sequence string
      acc_to_gene      : Entry_Isoform → gene (Gene_Gencode / Gene_Uniprot)
      gene_to_rows     : gene → list of (Entry_Isoform, Protein_ID, sequence)
      acc_to_protein_id: Entry_Isoform → Protein_ID ("RAF1-201")
    """
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"]
                     if c in loc_df.columns), None)

    acc_to_seq:        dict[str, str]       = {}
    acc_to_gene:       dict[str, str]       = {}
    acc_to_protein_id: dict[str, str]       = {}
    gene_to_rows:      dict[str, list]      = {}

    for _, row in loc_df.iterrows():
        acc  = str(row.get("Entry_Isoform", ""))
        pid  = str(row.get("Protein_ID",    ""))
        seq  = str(row.get("Sequence",      ""))
        gene = str(row.get(gene_col, "")) if gene_col else ""
        is_main = str(row.get("main_isoform", "")).lower() == "yes"

        if not acc or acc == "nan":
            continue

        if only_main and not is_main:
            continue

        if seq and seq != "nan":
            acc_to_seq[acc] = seq
        acc_to_gene[acc] = gene
        if pid and pid != "nan":
            acc_to_protein_id[acc] = pid

        gene_to_rows.setdefault(gene, [])
        entry = (acc, pid, seq if (seq and seq != "nan") else "")
        if entry not in gene_to_rows[gene]:
            gene_to_rows[gene].append(entry)

    return acc_to_seq, acc_to_gene, gene_to_rows, acc_to_protein_id


# ---------------------------------------------------------------------------
# Core mapping function
# ---------------------------------------------------------------------------

def _detect_acc_col(df: pd.DataFrame) -> str | None:
    """Find the column that holds the UniProt accession."""
    for c in ["Entry_Isoform", "Accession", "acc"]:
        if c in df.columns:
            return c
    return None


def _best_similar_window(region: str, tgt_seq: str, min_identity: float):
    """Find where `region` aligns (ungapped, same length) in `tgt_seq` with the
    highest per-residue identity. Returns (start_index_0based, identity) or
    (-1, best_identity_seen) if nothing reaches `min_identity`.

    An exact substring is the identity==1.0 case (fast-pathed). A ≥min_identity
    window lets a motif (e.g. a SLiM) transfer onto an isoform whose sequence has
    drifted slightly, which is precisely what 'homology similarity' captures."""
    L = len(region)
    if L == 0 or len(tgt_seq) < L:
        return -1, 0.0
    idx = tgt_seq.find(region)
    if idx != -1:
        return idx, 1.0
    best_idx, best_id = -1, 0.0
    for i in range(0, len(tgt_seq) - L + 1):
        win = tgt_seq[i:i + L]
        matches = sum(1 for a, b in zip(region, win) if a == b)
        ident = matches / L
        if ident > best_id:
            best_id, best_idx = ident, i
    if best_id >= min_identity:
        return best_idx, best_id
    return -1, best_id


def map_annotations_to_transcripts(
    annot_df:          pd.DataFrame,
    acc_col:           str,
    acc_to_seq:        dict[str, str],
    acc_to_gene:       dict[str, str],
    gene_to_rows:      dict[str, list],
    acc_to_protein_id: dict[str, str],
    has_position_range: bool = True,
    min_identity:      float = 0.9,
) -> pd.DataFrame:
    """
    Map each annotation row to all transcripts (Protein_IDs) of the same gene.
    Adds 'Protein_ID' column identifying the target transcript.

    Same UniProt accession  → direct mapping (identical protein).
    Different accession      → homology-similarity transfer, accepted only when
                               the region aligns at ≥ `min_identity` (default 0.9).
    """
    if annot_df.empty:
        out = annot_df.copy()
        out["Protein_ID"]        = ""
        out["homology_transfer"] = False
        out["homology_identity"] = ""
        out["mapping_type"]      = "direct"
        return out

    out_rows = []

    try:
        from tqdm import tqdm as _tqdm
        _annot_iter = _tqdm(annot_df.iterrows(), total=len(annot_df),
                            desc='Transcript map', unit='annot', leave=False)
    except ImportError:
        _annot_iter = annot_df.iterrows()

    for _, row in _annot_iter:
        src_acc = str(row.get(acc_col, ""))
        if not src_acc or src_acc == "nan":
            continue

        src_seq  = acc_to_seq.get(src_acc, "")
        src_gene = acc_to_gene.get(src_acc, "")
        src_pid  = acc_to_protein_id.get(src_acc, "")

        # All transcripts of the same gene
        siblings = gene_to_rows.get(src_gene, [])

        for tgt_acc, tgt_pid, tgt_seq in siblings:

            # --- same UniProt accession as source ---
            if tgt_acc == src_acc:
                # Always verify positions fit within target sequence length.
                # Multiple GENCODE transcripts can share one UniProt acc with
                # different lengths; drop the annotation if it overruns the target.
                if tgt_seq:
                    if has_position_range:
                        try:
                            e_check = int(float(row.get("End", 0)))
                        except (ValueError, TypeError):
                            e_check = 0
                        if e_check > len(tgt_seq):
                            continue
                    else:
                        try:
                            pos_check = int(float(row.get("Position", 0)))
                        except (ValueError, TypeError):
                            pos_check = 0
                        if pos_check > len(tgt_seq):
                            continue
                nr = row.to_dict()
                nr["Protein_ID"]        = tgt_pid
                nr["homology_transfer"] = False
                nr["homology_identity"] = "1.000"
                out_rows.append(nr)
                continue

            if not tgt_seq or not src_seq:
                continue

            if has_position_range:
                try:
                    s = int(float(row.get("Start", 0)))
                    e = int(float(row.get("End",   0)))
                except (ValueError, TypeError):
                    continue
                if s < 1 or e < s or s > len(src_seq):
                    continue

                region = src_seq[s - 1: e]
                if not region:
                    continue

                idx, ident = _best_similar_window(region, tgt_seq, min_identity)
                if idx == -1:
                    continue

                nr = row.to_dict()
                nr[acc_col]             = tgt_acc
                nr["Protein_ID"]        = tgt_pid
                nr["Start"]             = idx + 1
                nr["End"]               = idx + len(region)
                nr["homology_transfer"] = True
                nr["homology_identity"] = f"{ident:.3f}"
                out_rows.append(nr)

            else:
                # Positional (PTM-like): 3-aa context window
                try:
                    pos = int(float(row.get("Position", 0)))
                except (ValueError, TypeError):
                    continue
                if pos < 1 or pos > len(src_seq):
                    continue

                ctx_s   = max(0, pos - 2)
                ctx_e   = min(len(src_seq), pos + 1)
                context = src_seq[ctx_s: ctx_e]

                if not context:
                    continue
                c_idx = tgt_seq.find(context)
                if c_idx == -1:
                    continue

                offset  = pos - 1 - ctx_s
                new_pos = c_idx + offset + 1

                if new_pos < 1 or new_pos > len(tgt_seq):
                    continue

                nr = row.to_dict()
                nr[acc_col]             = tgt_acc
                nr["Protein_ID"]        = tgt_pid
                nr["Position"]          = new_pos
                nr["homology_transfer"] = True
                nr["homology_identity"] = "1.000"   # exact 3-aa context window
                out_rows.append(nr)

    if not out_rows:
        out = annot_df.copy()
        out["Protein_ID"]        = ""
        out["homology_transfer"] = False
        out["homology_identity"] = ""
        out["mapping_type"]      = "direct"
        return out

    df = pd.DataFrame(out_rows)
    # mapping_type categorises provenance: 'direct' (same UniProt accession /
    # native isoform) vs 'homology_similarity' (transferred from the main
    # isoform onto an alternative isoform by sequence homology ≥ threshold).
    df["mapping_type"] = df["homology_transfer"].map(
        lambda v: "homology_similarity" if bool(v) else "direct")

    # Ensure Protein_ID + mapping_type + homology columns appear first
    lead = ["Protein_ID", "mapping_type", "homology_transfer", "homology_identity", acc_col]
    col_order = (lead +
                 [c for c in annot_df.columns if c not in lead])
    return df.reindex(columns=col_order, fill_value="")


# ---------------------------------------------------------------------------
# Map a single annotation file
# ---------------------------------------------------------------------------

def map_file(
    path:              str,
    acc_to_seq:        dict,
    acc_to_gene:       dict,
    gene_to_rows:      dict,
    acc_to_protein_id: dict,
    output_path:       Path,
    label:             str,
    min_identity:      float = 0.9,
) -> int:
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as e:
        log.warning("%s: cannot read %s — %s", label, path, e)
        pd.DataFrame().to_csv(output_path, sep="\t", index=False)
        return 0

    if df.empty:
        df.to_csv(output_path, sep="\t", index=False)
        return 0

    acc_col = _detect_acc_col(df)
    if acc_col is None:
        log.warning("%s: no accession column — copying unchanged", label)
        df.to_csv(output_path, sep="\t", index=False)
        return len(df)

    has_range = ("Start" in df.columns and "End" in df.columns)
    mapped = map_annotations_to_transcripts(
        df, acc_col, acc_to_seq, acc_to_gene, gene_to_rows, acc_to_protein_id,
        has_position_range=has_range, min_identity=min_identity,
    )
    mapped.to_csv(output_path, sep="\t", index=False)
    n_mapped = len(mapped)
    log.info("%s: %d input rows → %d transcript-mapped rows", label, len(df), n_mapped)
    return n_mapped


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Module 5e: Map UniProt annotations → Gencode transcripts")
    p.add_argument("--loc_chrom",         required=True)
    p.add_argument("--elm",               required=True)
    p.add_argument("--dibs",              required=True)
    p.add_argument("--mfib",              required=True)
    p.add_argument("--phasepro",          required=True)
    p.add_argument("--uniprot_roi",       required=True)
    p.add_argument("--uniprot_bind",      required=True)
    p.add_argument("--ptm",               required=True)
    p.add_argument("--pfam",              required=True)
    p.add_argument("--disorder",          required=True,
                   help="CombinedDisorderNew.tsv — pass-through, already Protein_ID-keyed")
    p.add_argument("--disorder_pos",      default=None,
                   help="CombinedDisorderNew_Pos.tsv — pass-through")
    p.add_argument("--only_main_isoforms", action="store_true", default=False,
                   help="If set, only map annotations to main isoforms")
    p.add_argument("--homology_min_identity", type=float, default=0.9,
                   help="Min per-residue identity (0-1) to transfer a region onto "
                        "a different isoform as homology_similarity (default 0.9)")
    p.add_argument("--output_dir",        default=".")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading loc_chrom for transcript mapping…")
    loc_df = load_loc(args.loc_chrom)
    acc_to_seq, acc_to_gene, gene_to_rows, acc_to_protein_id = build_lookup(
        loc_df, only_main=args.only_main_isoforms)

    n_seqs = sum(1 for s in acc_to_seq.values() if s)
    log.info("Loaded %d isoforms with sequences across %d genes",
             n_seqs, len(gene_to_rows))

    stats = {}

    # Annotation files: map Entry_Isoform → Protein_ID (Pfam already keyed by Protein_ID)
    ann_files = [
        ("elm",          args.elm,          "elm.tsv"),
        ("dibs",         args.dibs,         "dibs.tsv"),
        ("mfib",         args.mfib,         "mfib.tsv"),
        ("phasepro",     args.phasepro,     "phasepro.tsv"),
        ("uniprot_roi",  args.uniprot_roi,  "uniprot_roi.tsv"),
        ("uniprot_bind", args.uniprot_bind, "uniprot_binding.tsv"),
        ("ptm",          args.ptm,          "ptm_merged.tsv"),
    ]

    for label, src, out_name in ann_files:
        n = map_file(src, acc_to_seq, acc_to_gene, gene_to_rows,
                     acc_to_protein_id, outdir / out_name, label,
                     min_identity=args.homology_min_identity)
        stats[label] = n

    # Pfam: computed per Protein_ID in ANNOTATION_MAP — pass through unchanged
    if args.pfam and Path(args.pfam).exists():
        pfam_df = pd.read_csv(args.pfam, sep="\t", dtype=str, nrows=1)
        shutil.copy(args.pfam, outdir / "pfam_domains.tsv")
        n_pfam = sum(1 for _ in open(args.pfam)) - 1
        stats["pfam"] = max(0, n_pfam)
        log.info("pfam: passed through (%d rows, Protein_ID-keyed=%s)",
                 max(0, n_pfam), "Protein_ID" in pfam_df.columns)
    else:
        pd.DataFrame().to_csv(outdir / "pfam_domains.tsv", sep="\t", index=False)
        stats["pfam"] = 0

    # Disorder files: already keyed by Protein_ID — pass through unchanged
    for src_path, out_name, key in [
        (args.disorder,     "CombinedDisorderNew.tsv",     "disorder"),
        (args.disorder_pos, "CombinedDisorderNew_Pos.tsv", "disorder_pos"),
    ]:
        if src_path and Path(src_path).exists():
            shutil.copy(src_path, outdir / out_name)
            n = sum(1 for _ in open(src_path)) - 1
            stats[key] = max(0, n)
            log.info("%s: passed through (%d rows)", key, max(0, n))
        else:
            pd.DataFrame().to_csv(outdir / out_name, sep="\t", index=False)
            stats[key] = 0

    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(outdir / "transcript_map_stats.tsv", sep="\t", index=False)
    log.info("Transcript mapping done: %s", stats)


if __name__ == "__main__":
    main()
