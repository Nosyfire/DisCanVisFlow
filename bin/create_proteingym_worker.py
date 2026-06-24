#!/usr/bin/env python3
"""
create_proteingym_worker.py — Map ProteinGym DMS (deep mutational scanning)
substitution scores onto the run's protein isoforms.

ProteinGym is the standard variant-effect-prediction benchmark: each assay
provides per-variant experimental fitness/function (DMS) scores plus a binarised
label (DMS_score_bin: 1 = functional / benign-like, 0 = deleterious-like) used as
the pathogenicity proxy.

Two mapping modes
-----------------
premapped (default):
  The input table is already keyed on Protein_ID (Gencode transcript), so this is
  a chunked filter to the run's proteins — every kept row is mapping_type=direct.
  Input columns: Protein_ID  uniprot_id  protein_variant  pos  DMS_score
                 DMS_score_bin  DMS_id

uniprot:
  The input (--proteingym_raw, from fetch_proteingym_worker.py) is UniProt-keyed.
  Each row is fanned out onto EVERY run isoform of that accession's gene:
    * if the isoform's residue at `pos` equals the variant's WT aa
      → mapping_type=direct (the native isoform), Protein_position = pos.
    * else, if the 3-aa context window (residue pos±1 in the row's reference
      sequence — taken from the native isoform) matches in the isoform's sequence
      at a (possibly shifted) location → mapping_type=homology_similarity with the
      remapped Protein_position. This lets DMS scores transfer onto drifted /
      alternative isoforms (mirrors create_transcript_map_worker.py's PTM window).
  Input columns: uniprot  gene_name  DMS_id  protein_variant  pos  DMS_score
                 DMS_score_bin

Output
------
  proteingym.tsv with columns:
    Protein_ID  Protein_position  protein_variant  DMS_score  DMS_score_bin
    DMS_id  uniprot_id  mapping_type
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

OUT_COLS = ["Protein_ID", "Protein_position", "protein_variant", "DMS_score",
            "DMS_score_bin", "DMS_id", "uniprot_id", "mapping_type"]

# leading WT residue + position, e.g. "G145R" → ("G", 145)
_MUT_RE = re.compile(r"^([A-Za-z])(\d+)")


def _wt_aa(variant):
    if not isinstance(variant, str):
        return ""
    m = _MUT_RE.match(variant.split(":")[0].strip())
    return m.group(1).upper() if m else ""


# ===========================================================================
# premapped mode (Protein_ID-keyed filter; original behaviour)
# ===========================================================================
def run_premapped(args, protein_ids, out_path):
    src = Path(args.proteingym)
    if not src.exists() or src.stat().st_size == 0 or not protein_ids:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
        log.warning("ProteinGym source missing/empty or no proteins — wrote empty %s", out_path)
        return

    kept = []
    n_scanned = 0
    for chunk in pd.read_csv(src, sep="\t", dtype=str, chunksize=args.chunksize):
        n_scanned += len(chunk)
        if "Protein_ID" not in chunk.columns:
            log.error("ProteinGym file has no Protein_ID column — aborting filter")
            break
        sub = chunk[chunk["Protein_ID"].isin(protein_ids)]
        if not sub.empty:
            kept.append(sub)

    if not kept:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
        log.info("ProteinGym: scanned %d rows, 0 matched run proteins", n_scanned)
        return

    df = pd.concat(kept, ignore_index=True)
    out = pd.DataFrame()
    out["Protein_ID"] = df["Protein_ID"]
    out["Protein_position"] = df.get("pos", "")
    out["protein_variant"] = df.get("protein_variant", "")
    out["DMS_score"] = df.get("DMS_score", "")
    out["DMS_score_bin"] = df.get("DMS_score_bin", "")
    out["DMS_id"] = df.get("DMS_id", "")
    out["uniprot_id"] = df.get("uniprot_id", "")
    out["mapping_type"] = "direct"

    out[OUT_COLS].to_csv(out_path, sep="\t", index=False)
    log.info("ProteinGym: scanned %d rows, kept %d for %d proteins → %s",
             n_scanned, len(out), len(protein_ids), out_path)


# ===========================================================================
# uniprot mode — fan out a UniProt-keyed table onto all run isoforms
# ===========================================================================
def build_isoform_index(seq_df):
    """Group run isoforms by gene, keyed off the UniProt base accession.

    Returns acc_to_isoforms : base_accession (e.g. 'P00001') →
              list of (Protein_ID, Entry_Isoform, sequence).
    Every isoform of the gene that owns the accession is included, so a row keyed
    by P00001 fans out to P00001 and P00001-2 alike."""
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"]
                     if c in seq_df.columns), None)

    gene_to_isoforms = {}
    gene_accessions = {}      # gene → set of base accessions present
    for _, row in seq_df.iterrows():
        pid = str(row.get("Protein_ID", "")).strip()
        iso = str(row.get("Entry_Isoform", "")).strip()
        seq = str(row.get("Sequence", "")).strip()
        gene = str(row.get(gene_col, "")).strip() if gene_col else ""
        if not pid or pid == "nan" or not iso or iso == "nan":
            continue
        if seq == "nan":
            seq = ""
        gene_to_isoforms.setdefault(gene, []).append((pid, iso, seq))
        base = iso.split("-")[0]
        gene_accessions.setdefault(gene, set()).add(base)

    acc_to_isoforms = {}
    for gene, isoforms in gene_to_isoforms.items():
        for base in gene_accessions.get(gene, set()):
            acc_to_isoforms.setdefault(base, [])
            for entry in isoforms:
                if entry not in acc_to_isoforms[base]:
                    acc_to_isoforms[base].append(entry)
    return acc_to_isoforms


def _map_position(pos, wt_aa, native_seq, iso_seq):
    """Map a 1-based protein position from a reference (native) isoform onto a
    target isoform sequence.

    Returns (new_pos_1based, mapping_type) or None if it cannot be placed.
      - direct: target residue at `pos` equals the WT aa (same coordinate frame).
      - homology_similarity: a 3-aa context window (pos±1 in native_seq) is found
        in iso_seq at a possibly-shifted location; the variant residue is remapped.
    """
    # Direct: target isoform carries the WT residue at the exact position.
    if 1 <= pos <= len(iso_seq) and wt_aa and iso_seq[pos - 1] == wt_aa:
        return (pos, "direct")

    # Homology: locate the native 3-aa context (pos-1, pos, pos+1) in iso_seq.
    if not native_seq or pos < 1 or pos > len(native_seq):
        return None
    if wt_aa and native_seq[pos - 1] != wt_aa:
        return None      # row's WT aa disagrees with its own reference → drop
    ctx_s = max(0, pos - 2)              # 0-based start of context window
    ctx_e = min(len(native_seq), pos + 1)
    context = native_seq[ctx_s:ctx_e]
    if not context:
        return None
    c_idx = iso_seq.find(context)
    if c_idx == -1:
        return None
    offset = (pos - 1) - ctx_s           # position of the variant within context
    new_pos = c_idx + offset + 1         # 1-based remapped position
    if new_pos < 1 or new_pos > len(iso_seq):
        return None
    return (new_pos, "homology_similarity")


def run_uniprot(args, protein_ids, out_path):
    raw = Path(args.proteingym_raw)
    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str).fillna("")
    if not raw.exists() or raw.stat().st_size == 0 or seq_df.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
        log.warning("ProteinGym raw missing/empty or no proteins — wrote empty %s", out_path)
        return

    acc_to_isoforms = build_isoform_index(seq_df)
    # ProteinGym keys assays by the UniProt entry MNEMONIC (e.g. 'P53_HUMAN'),
    # not the accession. Build a mnemonic → base-accession map from the seq table
    # so those rows resolve onto the run's isoforms (Entry_Name → Entry_Isoform).
    name_to_acc = {}
    if "Entry_Name" in seq_df.columns:
        for _, sr in seq_df.iterrows():
            nm = str(sr.get("Entry_Name", "")).strip()
            base = str(sr.get("Entry_Isoform", "")).split("-")[0].strip()
            if nm and base:
                name_to_acc.setdefault(nm, base)
    # native (canonical) sequence per accession: the isoform whose Entry_Isoform
    # has no '-suffix' (or the first listed) is the reference frame for that acc.
    native_seq = {}
    for base, isoforms in acc_to_isoforms.items():
        for pid, iso, seq in isoforms:
            if iso == base and seq:
                native_seq[base] = seq
                break
        if base not in native_seq:
            # fall back to any isoform whose base matches and that has sequence
            for pid, iso, seq in isoforms:
                if iso.split("-")[0] == base and seq:
                    native_seq[base] = seq
                    break

    out_rows = []
    n_scanned = 0
    for chunk in pd.read_csv(raw, sep="\t", dtype=str, chunksize=args.chunksize):
        n_scanned += len(chunk)
        for _, r in chunk.iterrows():
            raw_up = str(r.get("uniprot", "")).strip()
            acc = raw_up.split("-")[0].strip()
            isoforms = acc_to_isoforms.get(acc)
            if not isoforms:
                # raw_up may be a UniProt mnemonic (e.g. 'P53_HUMAN') → resolve
                resolved = name_to_acc.get(raw_up)
                if resolved:
                    acc = resolved
                    isoforms = acc_to_isoforms.get(acc)
            if not isoforms:
                continue
            variant = str(r.get("protein_variant", ""))
            try:
                pos = int(float(str(r.get("pos", "")).strip()))
            except (ValueError, TypeError):
                continue
            wt_aa = _wt_aa(variant)
            nat = native_seq.get(acc, "")
            for pid, iso, seq in isoforms:
                if not seq:
                    continue
                placed = _map_position(pos, wt_aa, nat, seq)
                if placed is None:
                    continue
                new_pos, mtype = placed
                out_rows.append({
                    "Protein_ID": pid,
                    "Protein_position": str(new_pos),
                    "protein_variant": variant,
                    "DMS_score": r.get("DMS_score", ""),
                    "DMS_score_bin": r.get("DMS_score_bin", ""),
                    "DMS_id": r.get("DMS_id", ""),
                    "uniprot_id": r.get("uniprot", ""),
                    "mapping_type": mtype,
                })

    df = pd.DataFrame(out_rows, columns=OUT_COLS).drop_duplicates()
    df.to_csv(out_path, sep="\t", index=False)
    n_direct = (df["mapping_type"] == "direct").sum() if not df.empty else 0
    n_homo = (df["mapping_type"] == "homology_similarity").sum() if not df.empty else 0
    log.info("ProteinGym(uniprot): scanned %d raw rows → %d mapped (%d direct, "
             "%d homology) across %d proteins → %s",
             n_scanned, len(df), n_direct, n_homo,
             df["Protein_ID"].nunique() if not df.empty else 0, out_path)


def main():
    ap = argparse.ArgumentParser(description="Map ProteinGym DMS scores to run proteins")
    ap.add_argument("--seq_table", required=True,
                    help="loc_chrom_with_names_isoforms_with_seq.tsv "
                         "(Protein_ID / Entry_Isoform / Sequence columns)")
    ap.add_argument("--mapping_mode", choices=["premapped", "uniprot"],
                    default="premapped",
                    help="premapped: filter a Protein_ID-keyed table (default). "
                         "uniprot: fan out a UniProt-keyed raw table onto isoforms.")
    ap.add_argument("--proteingym", help="Pre-mapped ProteinGym TSV (premapped mode)")
    ap.add_argument("--proteingym_raw",
                    help="UniProt-keyed raw TSV from fetch_proteingym_worker.py "
                         "(uniprot mode)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--chunksize", type=int, default=500_000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "proteingym.tsv"

    protein_ids = set(
        pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
        ["Protein_ID"].dropna()
    )

    if args.mapping_mode == "uniprot":
        if not args.proteingym_raw:
            ap.error("--proteingym_raw is required for --mapping_mode uniprot")
        run_uniprot(args, protein_ids, out_path)
    else:
        if not args.proteingym:
            ap.error("--proteingym is required for --mapping_mode premapped")
        run_premapped(args, protein_ids, out_path)


if __name__ == "__main__":
    main()
