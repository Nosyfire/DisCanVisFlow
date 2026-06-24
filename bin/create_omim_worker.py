#!/usr/bin/env python3
"""
Module 8b — OMIM disease + protein-level mutation annotation.

Outputs:
  omim_disease.tsv   — disease ontology rows (filter to run proteins)
  omim_mutations.tsv — protein-level OMIM variant rows (position + aa_change)
"""

import argparse
import csv
import logging
import re
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

# a phenotype entry in genemap2 "Phenotypes": "Name, 612345 (3), Inheritance"
_PHENO_RE = re.compile(r"^(.*?),\s*(\d{6})\s*(?:\(\d\))?", re.S)

# humsavar AA change "p.His52Arg" → (His, 52, Arg)
_HUMSAVAR_AA = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")
# OMIM id embedded in a disease string "[MIM:123456]" or "MIM:123456"
_MIM_RE = re.compile(r"MIM:?\s*(\d{6})")

_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Sec": "U", "Ter": "*",
}


def _base_acc(entry_isoform):
    return str(entry_isoform).split("-")[0] if entry_isoform else ""


def _load_isoforms_omim(seq_table):
    """Parse the seq table for isoform fan-out: per-gene isoform records
    (Protein_ID, base UniProt acc, sequence) + acc→gene + gene→canonical seq."""
    df = pd.read_csv(seq_table, sep="\t", dtype=str).fillna("")
    gene_col = next((c for c in ["Gene", "Gene_Gencode", "Gene_Uniprot"]
                     if c in df.columns), None)
    gene_rows, acc_to_gene, gene_canon_seq = {}, {}, {}
    for _, r in df.iterrows():
        acc = str(r.get("Entry_Isoform", ""))
        pid = str(r.get("Protein_ID", ""))
        seq = str(r.get("Sequence", ""))
        gene = str(r.get(gene_col, "")) if gene_col else ""
        if not acc or not pid or acc == "nan" or pid == "nan":
            continue
        base = _base_acc(acc)
        gene_rows.setdefault(gene, []).append(
            {"Protein_ID": pid, "base_acc": base,
             "Sequence": seq if seq != "nan" else ""})
        if base:
            acc_to_gene[base] = gene
        if seq and seq != "nan" and (str(r.get("main_isoform", "")).lower() == "yes"
                                     or "-" not in acc):
            gene_canon_seq.setdefault(gene, seq)
        if seq and seq != "nan":
            gene_canon_seq.setdefault(gene, seq)
    return gene_rows, acc_to_gene, gene_canon_seq


def parse_humsavar(handle):
    """Parse UniProt humsavar.txt → list of variant dicts. Only the fixed
    7-column data block (gene, AC, FTId, AA change, category, dbSNP, disease).
    Position/WT/ALT parsed from the p.XxxNNNYyy AA change."""
    out = []
    started = False
    for line in handle:
        # data block starts after the '____' ruler under the column header
        if not started:
            if line.startswith("___") or line.strip().startswith("_______"):
                started = True
            continue
        s = line.rstrip("\n")
        if not s.strip():
            # a blank line after data → footer/legend begins; stop
            if out:
                break
            continue
        parts = s.split(None, 6)
        if len(parts) < 6 or not parts[1] or not parts[2].startswith("VAR"):
            continue
        gene, acc, ftid, aa_change, category = parts[0], parts[1], parts[2], parts[3], parts[4]
        dbsnp = parts[5] if len(parts) > 5 else "-"
        disease = parts[6].strip() if len(parts) > 6 else "-"
        m = _HUMSAVAR_AA.search(aa_change)
        wt = _AA3TO1.get(m.group(1), "") if m else ""
        pos = m.group(2) if m else ""
        mim = _MIM_RE.search(disease)
        out.append({
            "gene": gene, "acc": acc, "ftid": ftid, "aa_change": aa_change,
            "category": category, "dbSNP": "" if dbsnp == "-" else dbsnp,
            "disease": "" if disease == "-" else disease,
            "wt": wt, "pos": pos, "MIMID": mim.group(1) if mim else "",
        })
    return out


def _gene_to_pids(seq_table):
    """Map every gene symbol in the run to its Protein_IDs (via the seq table)."""
    df = pd.read_csv(seq_table, sep="\t", dtype=str).fillna("")
    gene_col = next((c for c in ["Gene", "Gene_Gencode", "Gene_Uniprot"]
                     if c in df.columns), None)
    g2p = {}
    if gene_col:
        for _, r in df.iterrows():
            g = str(r.get(gene_col, "")).strip().upper()
            pid = str(r.get("Protein_ID", "")).strip()
            if g and pid:
                g2p.setdefault(g, []).append(pid)
    return g2p


def parse_genemap2(handle):
    """Parse OMIM genemap2.txt → list of (gene_symbol, disease_name, mim) tuples.

    genemap2 is tab-delimited with '#'-comment lines. The phenotype column
    holds ';'-separated entries 'Name, MIM (key), Inheritance'."""
    out = []
    for line in handle:
        if not line or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 13:
            continue
        # Approved Gene Symbol (col 8, 0-based) else first of Gene Symbols (col 6)
        approved = cols[8].strip() if len(cols) > 8 else ""
        symbols = [s.strip() for s in cols[6].split(",") if s.strip()] if len(cols) > 6 else []
        gene = (approved or (symbols[0] if symbols else "")).upper()
        if not gene:
            continue
        phenos = cols[12] if len(cols) > 12 else ""
        for entry in phenos.split(";"):
            entry = entry.strip()
            if not entry:
                continue
            m = _PHENO_RE.match(entry)
            if m:
                name, mim = m.group(1).strip(), m.group(2)
            else:
                name, mim = entry, ""
            name = name.lstrip("{[?").rstrip("}]").strip()
            if name:
                out.append((gene, name, mim))
    return out


def fan_out_variant(acc, pos, wt_aa, gene_rows, acc_to_gene, gene_canon_seq):
    """Map a UniProt-keyed variant (acc + 1-based pos + WT aa) onto every run
    isoform of its gene. Yields (Protein_ID, Protein_position, mapping_type)."""
    gene = acc_to_gene.get(_base_acc(acc))
    if gene is None:
        return
    canon = gene_canon_seq.get(gene, "")
    ctx = ""
    if canon and 1 <= pos <= len(canon):
        lo, hi = max(0, pos - 2), min(len(canon), pos + 1)
        ctx = canon[lo:hi]
    ctx_center = min(pos - 1, 1) if pos >= 1 else 0
    if ctx and wt_aa and 0 <= ctx_center < len(ctx):
        ctx = ctx[:ctx_center] + wt_aa + ctx[ctx_center + 1:]
    for rec in gene_rows.get(gene, []):
        seq = rec["Sequence"]
        if not seq:
            continue
        if 1 <= pos <= len(seq) and (not wt_aa or seq[pos - 1] == wt_aa):
            yield (rec["Protein_ID"], pos, "direct")
            continue
        if ctx:
            idx = seq.find(ctx)
            if idx != -1:
                yield (rec["Protein_ID"], idx + ctx_center + 1, "homology_similarity")


def run_raw_mode(args, disease_out, mutation_out):
    """Build OMIM disease + protein-level mutation tables from raw source files.

    Primary source = UniProt humsavar.txt (open): per-variant disease + position
    + dbSNP, fanned out onto all run isoforms. genemap2.txt (OMIM key-gated) adds
    gene→phenotype disease rows when present."""
    raw_dir = Path(args.omim_raw_dir) if args.omim_raw_dir else None
    humsavar = Path(args.humsavar) if args.humsavar \
        else (raw_dir / "humsavar.txt" if raw_dir else None)
    genemap2 = (raw_dir / "genemap2.txt") if raw_dir else None

    g2p = _gene_to_pids(args.seq_table)
    gene_rows, acc_to_gene, gene_canon_seq = _load_isoforms_omim(args.seq_table)

    disease_rows, mut_rows = [], []

    # ── humsavar: protein-level disease variants (the main reproducible source) ──
    if humsavar and humsavar.exists() and humsavar.stat().st_size > 0:
        with open(humsavar, encoding="utf-8", errors="replace") as fh:
            variants = parse_humsavar(fh)
        n_disease = 0
        for v in variants:
            if not v["disease"]:          # only disease/cancer-named rows
                continue
            try:
                pos = int(v["pos"])
            except (ValueError, TypeError):
                continue
            n_disease += 1
            for pid, p_prime, _mtype in fan_out_variant(
                    v["acc"], pos, v["wt"], gene_rows, acc_to_gene, gene_canon_seq):
                mut_rows.append({
                    "Protein_ID": pid, "Entry_Isoform": v["acc"],
                    "Protein_position": p_prime, "aa_change": v["aa_change"],
                    "Disease": v["disease"], "MIMID": v["MIMID"],
                    "dbSNP": v["dbSNP"], "FTId": v["ftid"]})
                disease_rows.append({"Protein_ID": pid, "Disease": v["disease"],
                                     "MIMID": v["MIMID"]})
        log.info("OMIM (humsavar): %d disease variants → %d isoform mutation rows",
                 n_disease, len(mut_rows))

    # ── genemap2: gene→phenotype disease ontology (optional enrichment) ──────────
    if genemap2 and genemap2.exists() and genemap2.stat().st_size > 0:
        with open(genemap2, encoding="utf-8", errors="replace") as fh:
            for gene, disease, mim in parse_genemap2(fh):
                for pid in g2p.get(gene, []):
                    disease_rows.append({"Protein_ID": pid, "Disease": disease,
                                         "MIMID": mim})

    disease_df = pd.DataFrame(disease_rows,
                              columns=["Protein_ID", "Disease", "MIMID"]).drop_duplicates()
    disease_df.to_csv(disease_out, sep="\t", index=False)
    mut_df = pd.DataFrame(mut_rows, columns=_MUTATION_COLS).drop_duplicates()
    mut_df.to_csv(mutation_out, sep="\t", index=False)
    log.info("OMIM (raw): %d disease rows, %d mutation rows for %d proteins",
             len(disease_df), len(mut_df),
             disease_df["Protein_ID"].nunique() if len(disease_df) else 0)

_DISEASE_COLS = [
    "Protein_ID", "Disease", "DOID", "DO Subset", "synonyms", "MIMID",
    "level1", "level2", "level3", "level4", "level5", "level6",
    "level7", "level8", "level9", "level10", "level11", "level12",
    "Disordered", "Ordered", "Total Mutations", "Disordered Percent", "Name",
]

_MUTATION_COLS = [
    "Protein_ID", "Entry_Isoform", "Protein_position", "aa_change",
    "Disease", "MIMID", "dbSNP", "FTId",
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table", required=True)
    p.add_argument("--mapping_mode", choices=["processed", "raw"], default="processed",
                   help="processed: pre-built OMIM disease/variant tables (default). "
                        "raw: parse FETCH_OMIM's genemap2.txt directly.")
    p.add_argument("--omim_table", help="processed OMIM disease table (processed mode)")
    p.add_argument("--omim_mutations", default=None,
                   help="Protein-level OMIM variants (omim_mapped.tsv)")
    p.add_argument("--omim_raw_dir", help="dir with raw humsavar.txt / genemap2.txt (raw mode)")
    p.add_argument("--humsavar", help="explicit path to UniProt humsavar.txt (raw mode)")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    disease_out = outdir / "omim_disease.tsv"
    mutation_out = outdir / "omim_mutations.tsv"

    if args.mapping_mode == "raw":
        if not args.omim_raw_dir and not args.humsavar:
            p.error("--mapping_mode raw needs --omim_raw_dir or --humsavar")
        run_raw_mode(args, disease_out, mutation_out)
        return
    if not args.omim_table:
        p.error("--omim_table is required for --mapping_mode processed")

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str, usecols=["Protein_ID"])
    protein_ids = set(seq_df["Protein_ID"].dropna())

    # ── Disease ontology table ───────────────────────────────────────────────
    empty_disease = pd.DataFrame(columns=["Protein_ID", "Disease", "MIMID"])
    src = Path(args.omim_table)
    if not src.exists() or src.stat().st_size == 0:
        log.info("OMIM disease table not found — writing empty output")
        empty_disease.to_csv(disease_out, sep="\t", index=False)
    else:
        df = pd.read_csv(src, sep="\t", dtype=str)
        if "Accession" not in df.columns:
            log.warning("No 'Accession' column in OMIM disease table")
            empty_disease.to_csv(disease_out, sep="\t", index=False)
        else:
            df = df[df["Accession"].isin(protein_ids)].copy()
            df.rename(columns={"Accession": "Protein_ID"}, inplace=True)
            keep = [c for c in _DISEASE_COLS if c in df.columns]
            df = df[keep].drop_duplicates()
            df.to_csv(disease_out, sep="\t", index=False)
            log.info("OMIM disease: %d rows for %d proteins",
                     len(df), df["Protein_ID"].nunique() if len(df) else 0)

    # ── Protein-level OMIM mutations ─────────────────────────────────────────
    empty_mut = pd.DataFrame(columns=["Protein_ID", "Protein_position", "aa_change", "Disease"])
    mut_src = Path(args.omim_mutations) if args.omim_mutations else None
    if not mut_src or not mut_src.exists() or mut_src.stat().st_size == 0:
        log.info("OMIM mutations table not found — writing empty output")
        empty_mut.to_csv(mutation_out, sep="\t", index=False)
        return

    mut_df = pd.read_csv(mut_src, sep="\t", dtype=str)
    if "Protein_ID" not in mut_df.columns:
        log.warning("No Protein_ID in OMIM mutations table")
        empty_mut.to_csv(mutation_out, sep="\t", index=False)
        return

    mut_df = mut_df[mut_df["Protein_ID"].isin(protein_ids)].copy()
    if "position" in mut_df.columns and "Protein_position" not in mut_df.columns:
        mut_df.rename(columns={"position": "Protein_position"}, inplace=True)
    keep_mut = [c for c in _MUTATION_COLS if c in mut_df.columns]
    mut_df = mut_df[keep_mut].drop_duplicates()
    mut_df.to_csv(mutation_out, sep="\t", index=False)
    log.info("OMIM mutations: %d rows for %d proteins",
             len(mut_df), mut_df["Protein_ID"].nunique() if len(mut_df) else 0)


if __name__ == "__main__":
    main()
