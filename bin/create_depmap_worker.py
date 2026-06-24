#!/usr/bin/env python3
"""
Module 8e — DepMap cancer cell line somatic mutation annotation.

Filters DepMap mutations to proteins in this run and expands each hit to all
GENCODE isoforms of the same gene via sequence homology (Module 4 parity).

Usage:
  create_depmap_worker.py
      --seq_table   <loc_chrom_with_names_isoforms_with_seq.tsv>
      --depmap_tsv  <mapped_filtered_mutations.tsv>
      --outdir      <output directory>

Output:
  depmap_mutations.tsv
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import pandas as pd

from mutation_mapping_lib import expand_protein_position_to_isoforms, load_gene_isoform_lookup

# p.A119D / p.Ala119Asp (1- or 3-letter WT) → (wt_residue_1letter, position)
_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}
_HGVSP_RE = re.compile(r"p\.\(?([A-Za-z]{1,3}?)(\d+)")


def _parse_hgvsp_short(hgvsp: str):
    """Return (wt_1letter, position) from an HGVSp string, or ("", None)."""
    m = _HGVSP_RE.match(str(hgvsp or ""))
    if not m:
        return "", None
    wt, pos = m.group(1), int(m.group(2))
    if len(wt) == 3:
        wt = _AA3TO1.get(wt.title(), "")
    return wt.upper(), pos


def _gene_canon_seq(gene_to_rows):
    """gene → first non-empty isoform sequence (anchor for the context window)."""
    canon = {}
    for gene, rows in gene_to_rows.items():
        for _acc, _pid, seq in rows:
            if seq:
                canon[gene] = seq
                break
    return canon


def _fan_out_gene_variant(gene, pos, wt, gene_to_rows, gene_canon_seq):
    """Map a gene-keyed protein variant (gene + 1-based pos + claimed WT residue)
    onto every run isoform of that gene, mirroring the OMIM/MaveDB fan-out: an
    exact residue match at the position is 'direct' (isoform_mapped=False); a
    match of the 3-AA context window elsewhere is a coordinate transfer
    (isoform_mapped=True). Yields (Protein_ID, new_pos, is_transfer)."""
    canon = gene_canon_seq.get(gene, "")
    ctx = ""
    if canon and 1 <= pos <= len(canon):
        lo, hi = max(0, pos - 2), min(len(canon), pos + 1)
        ctx = canon[lo:hi]
    ctx_center = min(pos - 1, 1) if pos >= 1 else 0
    if ctx and wt and 0 <= ctx_center < len(ctx):
        ctx = ctx[:ctx_center] + wt + ctx[ctx_center + 1:]
    for _acc, pid, seq in gene_to_rows.get(gene, []):
        if not seq:
            continue
        if 1 <= pos <= len(seq) and (not wt or seq[pos - 1] == wt):
            yield (pid, pos, False)
            continue
        if ctx:
            idx = seq.find(ctx)
            if idx != -1:
                yield (pid, idx + ctx_center + 1, True)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_KEEP_COLS = [
    "Protein_ID", "Chrom", "Start_Position", "End_Position",
    "HugoSymbol", "Protein_position", "HGVSp_Short",
    "VariantType", "VariantInfo", "DNAChange",
    "ModelID", "Hotspot", "EntrezGeneID",
    "Rescue", "RescueReason", "isoform_mapped",
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table", required=True)
    p.add_argument("--depmap_tsv", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "depmap_mutations.tsv"

    empty = pd.DataFrame(columns=["Protein_ID", "HGVSp_Short", "ModelID", "isoform_mapped"])

    src = Path(args.depmap_tsv)
    if not src.exists() or src.stat().st_size == 0:
        log.info("DepMap table not found — writing empty output")
        empty.to_csv(out, sep="\t", index=False)
        return

    gene_to_rows, pid_to_seq, pid_to_gene, _ = load_gene_isoform_lookup(args.seq_table)
    run_pids = set(pid_to_seq) | set(pid_to_gene)
    if not run_pids:
        # seq_table may have only Protein_ID column (e.g., in tests or stub runs)
        _st = pd.read_csv(args.seq_table, sep="\t", dtype=str)
        if "Protein_ID" in _st.columns:
            run_pids = set(_st["Protein_ID"].dropna())

    df = pd.read_csv(src, sep="\t", dtype=str)

    # ── Raw gene-keyed mode (fetch_depmap_worker.py output: HugoSymbol/HGVSp_Short,
    #    no Protein_ID) — map each variant onto the run's isoforms by gene + WT. ──
    if "Protein_ID" not in df.columns:
        if "HugoSymbol" not in df.columns:
            log.warning("DepMap raw table lacks both Protein_ID and HugoSymbol "
                        "columns — writing empty output")
            empty.to_csv(out, sep="\t", index=False)
            return
        run_genes = {pid_to_gene[p] for p in run_pids if p in pid_to_gene}
        canon = _gene_canon_seq(gene_to_rows)
        raw_keep = ["HugoSymbol", "HGVSp_Short", "ModelID", "Start_Position",
                    "EntrezGeneID", "Hotspot"]
        rows, seen = [], set()
        for _, row in df.iterrows():
            gene = str(row.get("HugoSymbol", ""))
            if run_genes and gene not in run_genes:
                continue
            hgvsp = row.get("HGVSp_Short", "")
            wt, pos = _parse_hgvsp_short(hgvsp)
            if pos is None:
                try:
                    pos = int(float(row.get("Protein_position", 0)))
                except (ValueError, TypeError):
                    continue
            for tgt_pid, tgt_pos, is_transfer in _fan_out_gene_variant(
                gene, pos, wt, gene_to_rows, canon
            ):
                key = (tgt_pid, str(row.get("ModelID", "")), str(hgvsp))
                if key in seen:
                    continue
                seen.add(key)
                orow = {c: row.get(c, "") for c in raw_keep if c in row.index}
                orow["Protein_ID"] = tgt_pid
                orow["Protein_position"] = str(tgt_pos)
                orow["isoform_mapped"] = "True" if is_transfer else "False"
                rows.append(orow)
        if not rows:
            empty.to_csv(out, sep="\t", index=False)
            return
        out_df = pd.DataFrame(rows)
        lead = ["Protein_ID", "Protein_position", "HugoSymbol", "HGVSp_Short",
                "ModelID", "Start_Position", "EntrezGeneID", "Hotspot",
                "isoform_mapped"]
        out_df = out_df[[c for c in lead if c in out_df.columns]]
        out_df.to_csv(out, sep="\t", index=False)
        log.info("DepMap (raw): %d rows for %d proteins", len(out_df),
                 out_df["Protein_ID"].nunique())
        return

    df = df[df["Protein_ID"].isin(run_pids)].copy()
    if df.empty:
        empty.to_csv(out, sep="\t", index=False)
        return

    expanded_rows: list[dict] = []
    seen: set[tuple] = set()

    for _, row in df.iterrows():
        primary_pid = str(row.get("Protein_ID", ""))
        try:
            primary_pos = int(float(row.get("Protein_position", 0)))
        except (ValueError, TypeError):
            continue
        gene = pid_to_gene.get(primary_pid, str(row.get("HugoSymbol", "")))
        for tgt_pid, tgt_pos, is_transfer in expand_protein_position_to_isoforms(
            primary_pid, primary_pos, gene, gene_to_rows, pid_to_seq
        ):
            key = (
                tgt_pid,
                str(row.get("ModelID", "")),
                str(row.get("HGVSp_Short", "")),
                str(row.get("Start_Position", "")),
            )
            if key in seen:
                continue
            seen.add(key)
            out_row = {c: row.get(c, "") for c in _KEEP_COLS if c in row.index}
            out_row["Protein_ID"] = tgt_pid
            out_row["Protein_position"] = str(tgt_pos)
            out_row["isoform_mapped"] = "True" if is_transfer else "False"
            expanded_rows.append(out_row)

    if not expanded_rows:
        empty.to_csv(out, sep="\t", index=False)
        return

    out_df = pd.DataFrame(expanded_rows)
    keep = [c for c in _KEEP_COLS if c in out_df.columns]
    out_df = out_df[keep]
    out_df.to_csv(out, sep="\t", index=False)
    log.info(
        "DepMap: %d rows for %d proteins (%d primary hits expanded)",
        len(out_df),
        out_df["Protein_ID"].nunique(),
        len(df),
    )


if __name__ == "__main__":
    main()
