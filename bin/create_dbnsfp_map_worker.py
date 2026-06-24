#!/usr/bin/env python3
"""
Module 8f — Map raw dbNSFP chr*.gz variants to Protein_ID via combined_map.map.

Reads per-chromosome gzipped dbNSFP files, validates reference AA against
combined_map.map (legacy dbNSFP_custom parity), and filters to proteins in run.

Usage:
  create_dbnsfp_map_worker.py
      --seq_table       <loc_chrom_with_names_isoforms_with_seq.tsv>
      --combined_map    <combined_map.map>
      --dbnsfp_raw_dir  <directory with dbNSFP chr*.gz files>
      --outdir          <output directory>

Output:
  pathogenicity_scores.tsv
"""

import argparse
import gzip
import logging
import sys
from pathlib import Path

import pandas as pd

from mutation_mapping_lib import (
    expand_protein_position_to_isoforms,
    load_combined_map_by_protein,
    load_gene_isoform_lookup,
    validate_hgvsp_aa,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_KEEP_COLS = [
    "Protein_ID", "chr", "Start_Position", "End_Position",
    "Protein_position", "aaref", "aaalt", "aapos",
    "ref", "alt", "rs_dbSNP",
    "AlphaMissense_score", "CADD_phred", "CADD_raw",
    "ClinPred_score", "ESM1b_score", "EVE_score",
    "Polyphen2_HDIV_score", "Polyphen2_HVAR_score",
    "PrimateAI_score", "SIFT_score",
    "VARITY_ER_LOO_score", "VARITY_ER_score",
    "VARITY_R_LOO_score", "VARITY_R_score",
    "REVEL_score", "REVEL_rankscore",
    "gMVP_score",
]

# dbNSFP header aliases (with or without # prefix)
_COL_MAP = {
    "#chr": "chr",
    "chr": "chr",
    "pos(1-based)": "Start_Position",
    "ref": "ref",
    "alt": "alt",
    "aaref": "aaref",
    "aaalt": "aaalt",
    "aapos": "aapos",
    "rs_dbSNP": "rs_dbSNP",
}


def _chromosomes_from_seq(seq_df: pd.DataFrame) -> set[str]:
    if "Chromosome" not in seq_df.columns:
        return set()
    chrs = set()
    for c in seq_df["Chromosome"].dropna().unique():
        c = str(c).strip()
        if not c:
            continue
        chrs.add(c if c.startswith("chr") else f"chr{c}")
    return chrs


def _find_chr_files(raw_dir: Path, chromosomes: set[str]) -> list[Path]:
    files = []
    for chrom in chromosomes:
        num = chrom.replace("chr", "")
        for pat in (f"chr{num}.gz", f"*chr{num}.gz", f"chr{num}.bed", f"*chr{num}.bed"):
            hits = sorted(raw_dir.glob(pat))
            if hits:
                files.append(hits[0])
                break
    return files


def _load_bed_header(raw_dir: Path) -> list[str] | None:
    for name in ("dbNSFP_custom.bed.header", "dbNSFP.bed.header", "header.txt"):
        hp = raw_dir / name
        if hp.is_file():
            return hp.read_text().strip().split("\t")
        hp2 = raw_dir.parent / name
        if hp2.is_file():
            return hp2.read_text().strip().split("\t")
    return None


def _normalize_chrom(raw: str) -> str:
    c = str(raw).strip()
    if not c:
        return c
    return c if c.startswith("chr") else f"chr{c}"


def map_dbnsfp_bed(
    bed_path: Path,
    protein_ids: set[str],
    pid_map: dict,
    chromosomes: set[str],
    header: list[str] | None,
    gene_to_rows: dict,
    pid_to_seq: dict,
    pid_to_gene: dict,
) -> list[dict]:
    """Stream legacy dbNSFP BED/TSV (genomic coords + scores) and map via combined_map."""
    kept: list[dict] = []
    col_idx = {c: i for i, c in enumerate(header)} if header else {}
    seen: set[tuple] = set()

    def _col(parts: list[str], name: str, default: str = "") -> str:
        if name in col_idx and col_idx[name] < len(parts):
            return parts[col_idx[name]]
        return default

    with open(bed_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = _normalize_chrom(_col(parts, "chr", parts[0]))
            if chromosomes and chrom not in chromosomes:
                continue
            pos_raw = _col(parts, "Start_Position", parts[1])
            try:
                pos = int(pos_raw)
            except ValueError:
                continue
            aaref = _col(parts, "aaref", "").strip()

            score_fields = {
                c: _col(parts, c, "")
                for c in _KEEP_COLS
                if c in col_idx and c not in {
                    "Protein_ID", "chr", "Start_Position", "End_Position", "Protein_position",
                }
            }

            for pid in protein_ids:
                hit = pid_map.get(pid, {}).get(pos)
                if not hit:
                    continue
                prot_pos, map_aa = hit
                if aaref and aaref != "." and not validate_hgvsp_aa(
                    f"p.{aaref}{prot_pos}", map_aa, prot_pos
                ):
                    continue
                gene = pid_to_gene.get(pid, "")
                for tgt_pid, tgt_pos, _transfer in expand_protein_position_to_isoforms(
                    pid, prot_pos, gene, gene_to_rows, pid_to_seq
                ):
                    key = (tgt_pid, pos, tgt_pos)
                    if key in seen:
                        continue
                    seen.add(key)
                    out = {
                        "Protein_ID": tgt_pid,
                        "Protein_position": str(tgt_pos),
                        "chr": chrom,
                        "Start_Position": str(pos),
                        "End_Position": str(pos + 1),
                        **score_fields,
                    }
                    kept.append(out)
    return kept


def _normalize_header(cols: list[str]) -> list[str]:
    out = []
    for c in cols:
        c = c.lstrip("#").strip()
        out.append(_COL_MAP.get(c, c))
    return out


def map_dbnsfp_file(
    gz_path: Path,
    protein_ids: set[str],
    pid_map: dict,
    gene_to_rows: dict,
    pid_to_seq: dict,
    pid_to_gene: dict,
) -> list[dict]:
    kept = []
    seen: set[tuple] = set()
    with gzip.open(gz_path, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                if "chr" in line.lower() or "CHROM" in line:
                    header = _normalize_header(line.lstrip("#").strip().split("\t"))
                continue
            if header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            row = dict(zip(header, parts))
            if row.get("Start_Position", ".") == ".":
                continue
            try:
                pos = int(row["Start_Position"])
            except ValueError:
                continue
            chrom = row.get("chr", "")
            if chrom and not str(chrom).startswith("chr"):
                chrom = f"chr{chrom}"

            for pid in protein_ids:
                hit = pid_map.get(pid, {}).get(pos)
                if not hit:
                    continue
                prot_pos, map_aa = hit
                aaref = str(row.get("aaref", "")).strip()
                if aaref and aaref != "." and not validate_hgvsp_aa(
                    f"p.{aaref}{prot_pos}", map_aa, prot_pos
                ):
                    continue
                gene = pid_to_gene.get(pid, "")
                for tgt_pid, tgt_pos, _transfer in expand_protein_position_to_isoforms(
                    pid, prot_pos, gene, gene_to_rows, pid_to_seq
                ):
                    key = (tgt_pid, pos, tgt_pos)
                    if key in seen:
                        continue
                    seen.add(key)
                    out = {c: row.get(c, "") for c in _KEEP_COLS if c in row or c == "Protein_ID"}
                    out["Protein_ID"] = tgt_pid
                    out["Protein_position"] = str(tgt_pos)
                    out["chr"] = chrom
                    out["Start_Position"] = str(pos)
                    out["End_Position"] = str(pos + 1)
                    kept.append(out)
    return kept


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seq_table", required=True)
    p.add_argument("--combined_map", required=True)
    p.add_argument("--dbnsfp_raw_dir", required=True)
    p.add_argument("--dbnsfp_bed_header", default=None,
                   help="Path to dbNSFP_custom.bed.header (score column names)")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "pathogenicity_scores.tsv"
    empty_cols = ["Protein_ID", "Protein_position", "AlphaMissense_score"]
    empty = pd.DataFrame(columns=empty_cols)

    raw_dir = Path(args.dbnsfp_raw_dir)
    if not raw_dir.exists():
        log.info("dbNSFP raw path not found — empty output")
        empty.to_csv(out, sep="\t", index=False)
        return

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)
    protein_ids = set(seq_df["Protein_ID"].dropna())
    chromosomes = _chromosomes_from_seq(seq_df)
    gene_to_rows, pid_to_seq, pid_to_gene, _ = load_gene_isoform_lookup(args.seq_table)

    log.info("Loading combined_map …")
    pid_map = load_combined_map_by_protein(args.combined_map)

    all_rows: list[dict] = []
    header = None
    if args.dbnsfp_bed_header and Path(args.dbnsfp_bed_header).is_file():
        header = Path(args.dbnsfp_bed_header).read_text().strip().split("\t")
        log.info("Using dbNSFP header: %d columns", len(header))

    if raw_dir.is_file():
        inputs = [raw_dir]
        if header is None:
            header = _load_bed_header(raw_dir.parent)
    else:
        if header is None:
            header = _load_bed_header(raw_dir)
        mapped_bed = raw_dir / "mapped_mutations.bed"
        if mapped_bed.is_file():
            log.info("Streaming raw dbNSFP BED: %s", mapped_bed.name)
            all_rows.extend(map_dbnsfp_bed(
                mapped_bed, protein_ids, pid_map, chromosomes, header,
                gene_to_rows, pid_to_seq, pid_to_gene))
            inputs = []
        else:
            inputs = _find_chr_files(raw_dir, chromosomes) if chromosomes else sorted(raw_dir.glob("*.gz"))

    for gz in inputs:
        if gz.suffix == ".bed":
            log.info("Scanning BED %s …", gz.name)
            all_rows.extend(map_dbnsfp_bed(
                gz, protein_ids, pid_map, chromosomes, header,
                gene_to_rows, pid_to_seq, pid_to_gene))
        else:
            log.info("Scanning %s …", gz.name)
            all_rows.extend(map_dbnsfp_file(
                gz, protein_ids, pid_map, gene_to_rows, pid_to_seq, pid_to_gene))

    if not all_rows and not inputs and raw_dir.is_dir():
        log.info("No dbNSFP chr files for run chromosomes — empty output")

    if all_rows:
        df = pd.DataFrame(all_rows)
        keep = [c for c in _KEEP_COLS if c in df.columns]
        df = df[keep]
    else:
        df = empty

    df.to_csv(out, sep="\t", index=False)
    log.info("Pathogenicity (raw map): %d rows for %d proteins",
             len(df), df["Protein_ID"].nunique() if len(df) else 0)


if __name__ == "__main__":
    main()
