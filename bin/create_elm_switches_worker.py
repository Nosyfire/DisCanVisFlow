#!/usr/bin/env python3
"""
Module 5p — Elm_Switches annotation mapping.

Parses the ELM molecular switches TSV and maps each switch to all matching
Gencode transcripts (Protein_IDs) via sequence homology transfer.

Input ELM switches TSV columns (from elm.eu.org/switches.tsv):
    Switch ID | Status | Interaction ID | Intramolecular |
    ID A | Bindingsite A ID | Bindingsite A Start | Bindingsite A End |
    ID B | Bindingsite B ID | Bindingsite B Start | Bindingsite B End |
    Affected interactor | Switch type | Switch subtype | Switch mechanism |
    Switch direction | Switch outcome direction | Switch outcome |
    Modification | Modification sites | Modifying enzymes | Effector |
    Cell cycle phase | Localisation | Pathway | PMID

Mapping rule:
  - Extract canonical UniProt accession from "ID A" (format "UNIPROT:P12345")
  - For the matching canonical isoform: keep Bindingsite A Start/End unchanged
  - For other isoforms of the same gene: remap coordinates via 3-AA context
    substring search (homology_transfer=True)

Output (elmswitches_mapped.tsv):
    Protein_ID | Entry_Isoform | homology_transfer | <all original columns>

Published to: mapped/annotations/

Usage:
    create_elm_switches_worker.py
        --seq_table  <loc_chrom_with_names_isoforms_with_seq.tsv>
        --switches   <elmswitches.tsv>  (or NO_FILE)
        --outdir     <output directory>
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_OUT_PREFIX = ["Protein_ID", "Entry_Isoform", "homology_transfer"]

_RAW_COLS = [
    "Switch ID", "Status", "Interaction ID", "Intramolecular",
    "ID A", "Bindingsite A ID", "Bindingsite A Start", "Bindingsite A End",
    "ID B", "Bindingsite B ID", "Bindingsite B Start", "Bindingsite B End",
    "Affected interactor", "Switch type", "Switch subtype", "Switch mechanism",
    "Switch direction", "Switch outcome direction", "Switch outcome",
    "Modification", "Modification sites", "Modifying enzymes", "Effector",
    "Cell cycle phase", "Localisation", "Pathway", "PMID",
]


def _extract_acc(id_a: str) -> str:
    """'UNIPROT:P12345' or 'UNIPROT:P12345-2' → 'P12345' or 'P12345-2'."""
    if ":" in str(id_a):
        return str(id_a).split(":", 1)[1].strip()
    return str(id_a).strip()


def _canonical(acc: str) -> str:
    """'P12345-2' → 'P12345'."""
    return acc.split("-")[0] if "-" in acc else acc


def _build_lookup(seq_df: pd.DataFrame):
    """
    Returns:
      canon_to_isoforms: canonical_acc → [(full_acc, Protein_ID, sequence)]
      gene_to_isoforms:  gene → [(full_acc, Protein_ID, sequence)]
      acc_to_gene:       full_acc → gene
    """
    pid_col = next((c for c in ["Protein_ID"] if c in seq_df.columns), None)
    acc_col = next((c for c in ["Entry_Isoform", "Accession"] if c in seq_df.columns), None)
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"] if c in seq_df.columns), None)
    if pid_col is None or acc_col is None:
        return {}, {}, {}

    canon_to_isoforms: dict[str, list] = defaultdict(list)
    gene_to_isoforms:  dict[str, list] = defaultdict(list)
    acc_to_gene:       dict[str, str]  = {}

    for _, row in seq_df.iterrows():
        pid  = str(row.get(pid_col, "")).strip()
        acc  = str(row.get(acc_col, "")).strip()
        seq  = str(row.get("Sequence", "")).strip()
        gene = str(row.get(gene_col, "")).strip() if gene_col else ""

        if not pid or pid == "nan" or not acc or acc == "nan":
            continue

        seq = seq if seq and seq != "nan" else ""
        entry = (acc, pid, seq)
        canon = _canonical(acc)

        if entry not in canon_to_isoforms[canon]:
            canon_to_isoforms[canon].append(entry)
        if gene and entry not in gene_to_isoforms[gene]:
            gene_to_isoforms[gene].append(entry)
        acc_to_gene[acc] = gene

    return dict(canon_to_isoforms), dict(gene_to_isoforms), acc_to_gene


def _remap_position(src_seq: str, tgt_seq: str, pos: int) -> int:
    """3-AA context window position remapping. Returns 0 if not found."""
    if pos < 1 or pos > len(src_seq):
        return 0
    ctx_s = max(0, pos - 2)
    ctx_e = min(len(src_seq), pos + 1)
    context = src_seq[ctx_s:ctx_e]
    if not context:
        return 0
    idx = tgt_seq.find(context)
    if idx == -1:
        return 0
    return idx + (pos - 1 - ctx_s) + 1


def _remap_region(src_seq: str, tgt_seq: str, start: int, end: int):
    """Substring search remapping for a region. Returns (new_start, new_end) or None."""
    if start < 1 or end < start or start > len(src_seq):
        return None
    region = src_seq[start - 1: end]
    if not region:
        return None
    idx = tgt_seq.find(region)
    if idx == -1:
        return None
    return idx + 1, idx + len(region)


def map_switches(sw_df: pd.DataFrame, seq_df: pd.DataFrame) -> pd.DataFrame:
    """Map ELM switches to all matching Protein_IDs."""
    if sw_df.empty:
        out_cols = _OUT_PREFIX + [c for c in _RAW_COLS if c in sw_df.columns]
        return pd.DataFrame(columns=out_cols)

    canon_to_isoforms, gene_to_isoforms, acc_to_gene = _build_lookup(seq_df)

    out_rows = []
    for _, row in sw_df.iterrows():
        id_a = str(row.get("ID A", "")).strip()
        if not id_a or id_a == "nan":
            continue

        full_acc = _extract_acc(id_a)
        canon    = _canonical(full_acc)

        # Find all isoforms that share this canonical accession
        isoform_entries = canon_to_isoforms.get(canon, [])
        if not isoform_entries:
            continue

        # Parse binding site A coordinates
        try:
            bs_start_raw = str(row.get("Bindingsite A Start", "")).strip()
            bs_end_raw   = str(row.get("Bindingsite A End",   "")).strip()
            # Handle semicolon-separated multiple sites
            starts = [int(x.strip()) for x in bs_start_raw.split(";") if x.strip()]
            ends   = [int(x.strip()) for x in bs_end_raw.split(";") if x.strip()]
        except (ValueError, AttributeError):
            starts, ends = [], []

        # Find source isoform sequence (prefer exact match, fall back to canonical)
        src_acc = full_acc
        src_entry = next(((a, p, s) for a, p, s in isoform_entries if a == full_acc), None)
        if src_entry is None:
            # Fallback: use the main (first canonical) isoform as source
            src_entry = isoform_entries[0]
            src_acc = src_entry[0]
        src_seq = src_entry[2]

        for tgt_acc, tgt_pid, tgt_seq in isoform_entries:
            nr = row.to_dict()

            if tgt_acc == src_acc:
                # Same accession → direct copy
                nr["Protein_ID"]        = tgt_pid
                nr["Entry_Isoform"]     = tgt_acc
                nr["homology_transfer"] = False
                out_rows.append(nr)
                continue

            # Different isoform → try coordinate remapping
            if not src_seq or not tgt_seq or not starts:
                continue

            new_starts, new_ends = [], []
            ok = True
            for s, e in zip(starts, ends):
                mapped = _remap_region(src_seq, tgt_seq, s, e)
                if mapped is None:
                    ok = False
                    break
                new_starts.append(mapped[0])
                new_ends.append(mapped[1])

            if not ok:
                continue

            nr["Protein_ID"]        = tgt_pid
            nr["Entry_Isoform"]     = tgt_acc
            nr["homology_transfer"] = True
            # Write back possibly-remapped coordinates
            if new_starts:
                nr["Bindingsite A Start"] = ";".join(str(x) for x in new_starts)
                nr["Bindingsite A End"]   = ";".join(str(x) for x in new_ends)
            out_rows.append(nr)

    if not out_rows:
        out_cols = _OUT_PREFIX + [c for c in _RAW_COLS if c in sw_df.columns]
        return pd.DataFrame(columns=out_cols)

    result = pd.DataFrame(out_rows)
    col_order = _OUT_PREFIX + [c for c in sw_df.columns if c not in _OUT_PREFIX]
    return result.reindex(columns=col_order, fill_value="")


def main():
    p = argparse.ArgumentParser(
        description="Module 5p: ELM switches → Elm_Switches TSV"
    )
    p.add_argument("--seq_table", required=True,
                   help="loc_chrom_with_names_isoforms_with_seq.tsv")
    p.add_argument("--switches",  required=True,
                   help="elmswitches.tsv (or NO_FILE)")
    p.add_argument("--outdir",    default=".")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "elmswitches_mapped.tsv"

    sw_path = Path(args.switches)
    out_cols = _OUT_PREFIX + _RAW_COLS

    if sw_path.name == "NO_FILE" or not sw_path.exists() or sw_path.stat().st_size < 10:
        log.info("ELM switches file not available — writing empty output")
        pd.DataFrame(columns=out_cols).to_csv(out, sep="\t", index=False)
        return

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)
    sw_df  = pd.read_csv(sw_path, sep="\t", dtype=str)
    log.info("Loaded %d ELM switch entries", len(sw_df))

    result = map_switches(sw_df, seq_df)
    result.to_csv(out, sep="\t", index=False)
    log.info("Done — %d mapped switch rows (%d proteins)",
             len(result),
             result["Protein_ID"].nunique() if not result.empty else 0)


if __name__ == "__main__":
    main()
