#!/usr/bin/env python3
"""
create_mavedb_worker.py — Map MaveDB single-mutant functional scores to the
run's protein isoforms.

Two selectable modes (--mapping_mode):

premapped (default — original behaviour, unchanged):
  Reads the large pre-mapped MaveDB table (already Protein_ID / Gencode-transcript
  keyed, e.g. Benchmark_Pathogenicity_Predictors/data/mavedb/
  mave_single_mutant_protein.tsv) in chunks and keeps only rows whose Protein_ID
  is in this run. Because the source is already keyed on the Gencode transcript,
  this is a direct mapping (mapping_type=direct), not a homology transfer.

uniprot (reproducible pipeline — fresh fetch_mavedb_worker.py output):
  Reads the UniProt-keyed raw table (mavedb_raw.tsv: uniprot, gene_name, urn,
  mavedb_id, prot_expr, protein_start, score, is_double_mutant) and maps each
  raw row's (UniProt accession + 1-based protein position + WT aa parsed from
  prot_expr) onto every run isoform of that accession's gene:
    - residue at that position equals the WT aa → mapping_type=direct,
      Protein_position = pos.
    - otherwise the 3-AA context window (pos-1,pos,pos+1) taken from the SOURCE
      canonical/main isoform sequence is located in the target isoform; if found
      at shifted position p' → mapping_type=homology_similarity, Protein_position=p'.
    - no match → that isoform is skipped.

Output
------
  mavedb.tsv with columns:
    Protein_ID  Protein_position  prot_expr  score  mavedb_id  urn
    gene_name  uniprot  Transcript_ID  is_double_mutant  mapping_type
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

# p.Arg89Tyr / p.(Arg89Tyr) / p.R89Y → position 89
_POS_RE = re.compile(r"p\.\(?[A-Za-z]{1,3}(\d+)")
# Capture the WT residue (1- or 3-letter) immediately after the 'p.' prefix.
_WT_RE = re.compile(r"p\.\(?([A-Za-z]{1,3})\d+")

OUT_COLS = ["Protein_ID", "Protein_position", "prot_expr", "score",
            "mavedb_id", "urn", "gene_name", "uniprot", "Transcript_ID",
            "is_double_mutant", "mapping_type"]

# 3-letter → 1-letter amino-acid map (incl. Ter/* for nonsense).
_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Xaa": "X",
}


def _pos_from_hgvs(expr):
    if not isinstance(expr, str):
        return ""
    m = _POS_RE.search(expr)
    return m.group(1) if m else ""


def _wt_aa_from_hgvs(expr):
    """Parse the wild-type residue (1-letter) from a protein HGVS like
    'p.Arg72Pro' (→ 'R') or 'p.R72P' (→ 'R'). Returns '' if unparseable."""
    if not isinstance(expr, str):
        return ""
    m = _WT_RE.search(expr)
    if not m:
        return ""
    tok = m.group(1)
    if len(tok) == 1:
        return tok.upper()
    return _AA3TO1.get(tok.capitalize(), "")


# ---------------------------------------------------------------------------
# uniprot mode helpers
# ---------------------------------------------------------------------------

def _base_acc(entry_isoform: str) -> str:
    """Base UniProt accession of an Entry_Isoform ('P04049-2' → 'P04049')."""
    return str(entry_isoform).split("-")[0] if entry_isoform else ""


def load_isoforms(seq_table: str):
    """Parse loc_chrom_with_names_isoforms_with_seq.tsv into per-gene isoform
    records. Returns (gene_rows, acc_to_gene, gene_canon_seq):
      gene_rows      : gene -> list of dicts (Protein_ID, Entry_Isoform, base_acc,
                       Sequence, is_main)
      acc_to_gene    : base UniProt accession -> gene
      gene_canon_seq : gene -> canonical/main isoform sequence (source for the
                       3-AA context window)
    """
    df = pd.read_csv(seq_table, sep="\t", dtype=str)
    df = df.fillna("")
    gene_col = next((c for c in ["Gene", "Gene_Gencode", "Gene_Uniprot"]
                     if c in df.columns), None)

    gene_rows: dict = {}
    acc_to_gene: dict = {}
    gene_canon_seq: dict = {}

    for _, row in df.iterrows():
        acc = str(row.get("Entry_Isoform", ""))
        pid = str(row.get("Protein_ID", ""))
        seq = str(row.get("Sequence", ""))
        gene = str(row.get(gene_col, "")) if gene_col else ""
        is_main = str(row.get("main_isoform", "")).lower() == "yes"
        if not acc or acc == "nan" or not pid or pid == "nan":
            continue
        base = _base_acc(acc)
        rec = {"Protein_ID": pid, "Entry_Isoform": acc, "base_acc": base,
               "Sequence": seq if seq != "nan" else "", "is_main": is_main}
        gene_rows.setdefault(gene, []).append(rec)
        if base:
            acc_to_gene[base] = gene
        # canonical = main isoform, or the accession without an isoform suffix
        if seq and seq != "nan":
            if is_main or "-" not in acc:
                gene_canon_seq.setdefault(gene, seq)
            gene_canon_seq.setdefault(gene, gene_canon_seq.get(gene, seq))

    return gene_rows, acc_to_gene, gene_canon_seq


def map_uniprot_row(acc, pos, wt_aa, gene_rows, acc_to_gene, gene_canon_seq):
    """Fan one raw MaveDB row (UniProt accession + 1-based pos + WT aa) out to
    every run isoform of that accession's gene. Yields per-isoform dicts with
    Protein_ID, Protein_position, mapping_type.

    direct             : the isoform residue at `pos` equals `wt_aa`.
    homology_similarity: the 3-AA context window (pos-1,pos,pos+1) from the gene's
                         canonical sequence is found elsewhere in the isoform.
    """
    base = _base_acc(acc)
    gene = acc_to_gene.get(base)
    if gene is None:
        return
    canon = gene_canon_seq.get(gene, "")
    # 3-AA context window centred on `pos` (1-based) from the source canonical seq.
    ctx = ""
    if canon and 1 <= pos <= len(canon):
        lo = max(0, pos - 2)               # pos-1 (0-based)
        hi = min(len(canon), pos + 1)      # pos+1 inclusive (0-based exclusive)
        ctx = canon[lo:hi]
    # offset of the WT residue within the context window (0,1, or 2)
    ctx_center = min(pos - 1, 1) if pos >= 1 else 0
    # Embed the CLAIMED WT residue at the window centre (canonical flanks + claimed
    # WT). If the claim disagrees with the canonical residue, this window exists
    # nowhere, so a mismatched variant is correctly skipped rather than spuriously
    # "recovered" against the canonical's own residue.
    if ctx and wt_aa and 0 <= ctx_center < len(ctx):
        ctx = ctx[:ctx_center] + wt_aa + ctx[ctx_center + 1:]

    for rec in gene_rows.get(gene, []):
        seq = rec["Sequence"]
        if not seq:
            continue
        # direct: residue at pos matches the WT aa
        if 1 <= pos <= len(seq) and (not wt_aa or seq[pos - 1] == wt_aa):
            yield {"Protein_ID": rec["Protein_ID"],
                   "Protein_position": pos, "mapping_type": "direct"}
            continue
        # homology: locate the context window in this isoform
        if len(ctx) >= 1 and ctx:
            idx = seq.find(ctx)
            if idx != -1:
                p_prime = idx + ctx_center + 1   # 1-based position of WT residue
                yield {"Protein_ID": rec["Protein_ID"],
                       "Protein_position": p_prime,
                       "mapping_type": "homology_similarity"}
                continue
        # no match → skip this isoform


def run_uniprot_mode(args, out_path):
    raw = Path(args.mavedb_raw)
    if not raw.exists() or raw.stat().st_size == 0:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
        log.warning("mavedb_raw missing/empty — wrote empty %s", out_path)
        return

    gene_rows, acc_to_gene, gene_canon_seq = load_isoforms(args.seq_table)

    df = pd.read_csv(raw, sep="\t", dtype=str).fillna("")
    out_rows = []
    n_in = len(df)
    for _, r in df.iterrows():
        acc = str(r.get("uniprot", ""))
        expr = str(r.get("prot_expr", ""))
        wt_aa = _wt_aa_from_hgvs(expr)
        pos_str = _pos_from_hgvs(expr) or str(r.get("protein_start", ""))
        try:
            pos = int(float(pos_str))
        except (ValueError, TypeError):
            continue
        if pos < 1:
            continue
        for hit in map_uniprot_row(acc, pos, wt_aa, gene_rows,
                                   acc_to_gene, gene_canon_seq):
            out_rows.append({
                "Protein_ID": hit["Protein_ID"],
                "Protein_position": hit["Protein_position"],
                "prot_expr": expr,
                "score": r.get("score", ""),
                "mavedb_id": r.get("mavedb_id", ""),
                "urn": r.get("urn", ""),
                "gene_name": r.get("gene_name", ""),
                "uniprot": acc,
                "Transcript_ID": "",
                "is_double_mutant": r.get("is_double_mutant", ""),
                "mapping_type": hit["mapping_type"],
            })

    out = pd.DataFrame(out_rows, columns=OUT_COLS)
    out.to_csv(out_path, sep="\t", index=False)
    n_direct = (out["mapping_type"] == "direct").sum() if not out.empty else 0
    n_hom = (out["mapping_type"] == "homology_similarity").sum() if not out.empty else 0
    log.info("MaveDB (uniprot): %d raw rows → %d isoform rows (%d direct, "
             "%d homology_similarity) → %s",
             n_in, len(out), n_direct, n_hom, out_path)


# ---------------------------------------------------------------------------
# premapped mode
# ---------------------------------------------------------------------------

def run_premapped_mode(args, out_path):
    protein_ids = set(
        pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
        ["Protein_ID"].dropna()
    )

    src = Path(args.mavedb)
    if not src.exists() or src.stat().st_size == 0 or not protein_ids:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
        log.warning("MaveDB source missing/empty or no proteins — wrote empty %s", out_path)
        return

    kept = []
    n_scanned = 0
    for chunk in pd.read_csv(src, sep="\t", dtype=str, chunksize=args.chunksize):
        n_scanned += len(chunk)
        if "Protein_ID" not in chunk.columns:
            log.error("MaveDB file has no Protein_ID column — aborting filter")
            break
        sub = chunk[chunk["Protein_ID"].isin(protein_ids)]
        if not sub.empty:
            kept.append(sub)

    if not kept:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
        log.info("MaveDB: scanned %d rows, 0 matched run proteins", n_scanned)
        return

    df = pd.concat(kept, ignore_index=True)
    out = pd.DataFrame()
    out["Protein_ID"] = df["Protein_ID"]
    # Authoritative 1-based position from the protein HGVS; fall back to protein_start
    out["Protein_position"] = df["prot_expr"].map(_pos_from_hgvs)
    if "protein_start" in df.columns:
        out.loc[out["Protein_position"] == "", "Protein_position"] = \
            df.loc[out["Protein_position"] == "", "protein_start"]
    out["prot_expr"] = df.get("prot_expr", "")
    out["score"] = df.get("score", "")
    out["mavedb_id"] = df.get("mavedb_id", "")
    out["urn"] = df.get("urn", "")
    out["gene_name"] = df.get("gene_name", "")
    out["uniprot"] = df.get("uniprot", "")
    out["Transcript_ID"] = df.get("Transcript ID", "")
    out["is_double_mutant"] = df.get("is_double_mutant", "")
    out["mapping_type"] = "direct"

    out[OUT_COLS].to_csv(out_path, sep="\t", index=False)
    log.info("MaveDB: scanned %d rows, kept %d for %d proteins → %s",
             n_scanned, len(out), len(protein_ids), out_path)


def main():
    ap = argparse.ArgumentParser(
        description="Map MaveDB single-mutant functional scores to run isoforms")
    ap.add_argument("--mapping_mode", choices=["premapped", "uniprot"],
                    default="premapped",
                    help="premapped: filter a Protein_ID-keyed table (default). "
                         "uniprot: fan out a UniProt-keyed fetch_mavedb_worker.py "
                         "table onto all run isoforms.")
    ap.add_argument("--seq_table", required=True,
                    help="loc_chrom_with_names_isoforms_with_seq.tsv")
    ap.add_argument("--outdir", required=True,
                    help="output directory (writes mavedb.tsv)")
    # premapped mode
    ap.add_argument("--mavedb",
                    help="pre-mapped Protein_ID-keyed MaveDB table (premapped mode)")
    ap.add_argument("--chunksize", type=int, default=500_000,
                    help="chunk size for the premapped table scan")
    # uniprot mode
    ap.add_argument("--mavedb_raw",
                    help="UniProt-keyed raw table from fetch_mavedb_worker.py "
                         "(uniprot mode)")
    args = ap.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "mavedb.tsv"

    if args.mapping_mode == "uniprot":
        if not args.mavedb_raw:
            ap.error("--mavedb_raw is required for --mapping_mode uniprot")
        run_uniprot_mode(args, out_path)
    else:
        if not args.mavedb:
            ap.error("--mavedb is required for --mapping_mode premapped")
        run_premapped_mode(args, out_path)


if __name__ == "__main__":
    main()
