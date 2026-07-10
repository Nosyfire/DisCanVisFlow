#!/usr/bin/env python3
"""
Module 8f — Map raw dbNSFP variants to Protein_ID via combined_map.map.

Two raw input shapes are supported:
  * a single merged dbNSFP 5.x gzip (e.g. dbNSFP5.3.1a_grch38.gz, ~50 GB, 505
    cols) — auto-detected; streamed once via an inverted (chr, genomic_pos)
    index built from combined_map.map, keeping all predictor scores + rankscores
    + CADD + conservation + gnomAD 4.1 joint AF (see select_keep_columns).
  * a directory of legacy per-chromosome chr*.gz files (dbNSFP 4.x).
Validates reference AA against combined_map.map and filters to proteins in run.

Usage:
  create_dbnsfp_map_worker.py
      --seq_table       <loc_chrom_with_names_isoforms_with_seq.tsv>
      --combined_map    <combined_map.map>
      --dbnsfp_raw_dir  <merged dbNSFP 5.x .gz file OR directory of chr*.gz>
      --outdir          <output directory>
      --n_cpu           <number of parallel chromosome workers, per-chr mode (default: 1)>

Output:
  dbnsfp_scores.tsv
"""

import argparse
import gzip
import logging
import multiprocessing as mp
import os
import re
import shutil
import subprocess
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

# Module-level shared state — set in main() before Pool creation.
# Fork-inherited by worker processes via copy-on-write; no pickling needed.
_g_pid_map: dict = {}
_g_gene_to_rows: dict = {}
_g_pid_to_seq: dict = {}
_g_pid_to_gene: dict = {}
_g_chr_to_pids: dict = {}   # {chrN: set of protein_ids on that chromosome}
_g_protein_ids: set = set() # fallback if chr-level lookup unavailable
_g_chromosomes: set = set()
_g_bed_header: list | None = None


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


def _chr_to_pids_from_seq(seq_df: pd.DataFrame) -> dict[str, set[str]]:
    """Map each chromosome to the set of Protein_IDs located on it."""
    if "Chromosome" not in seq_df.columns or "Protein_ID" not in seq_df.columns:
        return {}
    tmp = seq_df[["Protein_ID", "Chromosome"]].dropna()
    tmp = tmp[tmp["Protein_ID"].str.strip() != ""]
    tmp = tmp[tmp["Chromosome"].str.strip() != ""]
    tmp = tmp.copy()
    tmp["Chromosome"] = tmp["Chromosome"].apply(
        lambda c: str(c) if str(c).startswith("chr") else f"chr{c}"
    )
    return tmp.groupby("Chromosome")["Protein_ID"].apply(set).to_dict()


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


def _chrom_from_path(gz_path: Path) -> str | None:
    """Extract chromosome name from filename, e.g. chr5 from dbNSFP4.8a_variant.chr5.gz."""
    m = re.search(r"\b(chr(?:\d+|X|Y|M|MT))\b", gz_path.name, re.IGNORECASE)
    return m.group(1) if m else None


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


def _process_one_file(gz_path_str: str) -> list[dict]:
    """Worker function for multiprocessing.Pool — inherits globals via fork (COW).

    Filters protein_ids to only those on this chromosome before scanning,
    reducing the inner loop from ~110k proteins to ~2-5k per chromosome.
    """
    gz_path = Path(gz_path_str)
    chrom = _chrom_from_path(gz_path)
    if chrom and _g_chr_to_pids:
        pids = _g_chr_to_pids.get(chrom, set())
        if not pids:
            log.info("pid=%d no proteins on %s — skip", os.getpid(), chrom)
            return []
    else:
        pids = _g_protein_ids

    log.info("pid=%d scanning %s (%d proteins) …", os.getpid(), gz_path.name, len(pids))
    if gz_path.suffix == ".bed":
        return map_dbnsfp_bed(
            gz_path, pids, _g_pid_map, _g_chromosomes, _g_bed_header,
            _g_gene_to_rows, _g_pid_to_seq, _g_pid_to_gene,
        )
    return map_dbnsfp_file(
        gz_path, pids, _g_pid_map,
        _g_gene_to_rows, _g_pid_to_seq, _g_pid_to_gene,
    )


# ─────────────────────────────────────────────────────────────────────────
# Single merged dbNSFP file path (dbNSFP 5.x academic release: one 50 GB gzip,
# 505 columns, no per-chr split, no tabix). We stream it once against an
# inverted (chr, genomic_pos) index built from combined_map.map — O(1) per line
# instead of O(#proteins) — and keep all predictor scores + rankscores + CADD +
# conservation + gnomAD 4.1 joint AF.
# ─────────────────────────────────────────────────────────────────────────

# Identity columns carried through from the dbNSFP row (in output-friendly names).
_ID_KEEP = {
    "ref": "ref", "alt": "alt", "aaref": "aaref", "aaalt": "aaalt",
    "aapos": "aapos", "rs_dbSNP": "rs_dbSNP",
}
# Synthesized columns (computed from the mapping, not read from the file).
_SYNTH_COLS = ["Protein_ID", "Protein_position", "chr", "Start_Position", "End_Position"]
# gnomAD allele-frequency columns kept (modern joint release only).
_GNOMAD_KEEP = {"gnomAD4.1_joint_AF", "gnomAD4.1_joint_POPMAX_AF"}


def _keep_column(name: str) -> bool:
    """Pattern rule: predictor scores + rankscores + CADD + conservation + gnomAD AF."""
    return (
        name.endswith("_score")
        or name.endswith("_rankscore")
        or name in ("CADD_raw", "CADD_phred")
        or name.startswith(("GERP", "phyloP", "phastCons"))
        or name in _GNOMAD_KEEP
    )


def select_keep_columns(raw_header: list[str]):
    """Given the dbNSFP header, return (chr_idx, pos_idx, aaref_idx, file_cols).

    file_cols is an ordered list of (out_name, col_index) for every column kept
    from the file (identity cols + all scores/rankscores/CADD/conservation/gnomAD
    AF), excluding the chr/pos columns which are emitted via the synthesized block.
    """
    cols = [c.lstrip("#").strip() for c in raw_header]
    chr_idx = pos_idx = aaref_idx = None
    for i, c in enumerate(cols):
        if c in ("chr", "#chr") and chr_idx is None:
            chr_idx = i
        elif c == "pos(1-based)":
            pos_idx = i
        if c == "aaref":
            aaref_idx = i

    file_cols: list[tuple[str, int]] = []
    for i, c in enumerate(cols):
        if i == chr_idx or i == pos_idx:
            continue
        if c in _ID_KEEP:
            file_cols.append((_ID_KEEP[c], i))
        elif _keep_column(c):
            file_cols.append((c, i))
    return chr_idx, pos_idx, aaref_idx, file_cols


def _pid_to_chrom_from_seq(seq_df: pd.DataFrame) -> dict[str, str]:
    if "Chromosome" not in seq_df.columns or "Protein_ID" not in seq_df.columns:
        return {}
    out: dict[str, str] = {}
    for pid, chrom in zip(seq_df["Protein_ID"], seq_df["Chromosome"]):
        pid = str(pid).strip()
        chrom = str(chrom).strip()
        if not pid or not chrom or chrom.lower() == "nan":
            continue
        out[pid] = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return out


def build_gpos_index(pid_map: dict, pid_to_chrom: dict) -> dict:
    """Invert combined_map: {(chrom, genomic_pos): [(pid, protein_pos, aa), ...]}."""
    index: dict = {}
    for pid, posmap in pid_map.items():
        chrom = pid_to_chrom.get(pid)
        if not chrom:
            continue
        for gpos, (ppos, aa) in posmap.items():
            index.setdefault((chrom, gpos), []).append((pid, ppos, aa))
    return index


def _open_dbnsfp(path: Path):
    """Yield decompressed text lines, using pigz when available (much faster)."""
    if shutil.which("pigz"):
        proc = subprocess.Popen(
            ["pigz", "-dc", str(path)],
            stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
        )
        return proc.stdout, proc
    return gzip.open(path, "rt"), None


def stream_merged_dbnsfp(
    gz_path: Path,
    gpos_index: dict,
    out_path: Path,
) -> tuple[int, int]:
    """Stream a single merged dbNSFP gzip once, writing mapped rows directly to
    out_path. Returns (rows_written, distinct_proteins). Memory stays flat — the
    only accumulator is the small set of Protein_IDs seen (~20 k), so the caller
    never has to re-read the (potentially >100 GB) output to count proteins."""
    fh, proc = _open_dbnsfp(gz_path)
    n_written = 0
    proteins_seen: set = set()
    try:
        # First non-comment-consumed line is the header (starts with '#chr').
        header_line = fh.readline()
        while header_line and not header_line.lstrip().startswith("#"):
            header_line = fh.readline()
        if not header_line:
            _write_empty(out_path)
            return 0, 0
        raw_header = header_line.rstrip("\n").split("\t")
        chr_idx, pos_idx, aaref_idx, file_cols = select_keep_columns(raw_header)
        if chr_idx is None or pos_idx is None:
            log.error("dbNSFP header missing #chr / pos(1-based) — aborting")
            _write_empty(out_path)
            return 0, 0
        fast = (chr_idx == 0 and pos_idx == 1)
        col_idxs = [idx for _, idx in file_cols]
        out_header = _SYNTH_COLS + [name for name, _ in file_cols]

        with open(out_path, "w", encoding="utf-8") as out:
            out.write("\t".join(out_header) + "\n")
            for line in fh:
                if not line or line[0] == "#":
                    continue
                if fast:
                    pre = line.split("\t", 2)
                    if len(pre) < 3:
                        continue
                    chrom_raw, pos_raw = pre[0], pre[1]
                else:
                    parts0 = line.split("\t")
                    if len(parts0) <= max(chr_idx, pos_idx):
                        continue
                    chrom_raw, pos_raw = parts0[chr_idx], parts0[pos_idx]
                try:
                    pos = int(pos_raw)
                except ValueError:
                    continue
                chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"
                entries = gpos_index.get((chrom, pos))
                if not entries:
                    continue

                full = line.rstrip("\n").split("\t")
                nfull = len(full)
                file_vals = [full[j] if j < nfull else "" for j in col_idxs]
                aaref = (full[aaref_idx].strip()
                         if aaref_idx is not None and aaref_idx < nfull else "")
                pos_s = str(pos)
                end_s = str(pos + 1)
                seen: set = set()
                # Emit index entries DIRECTLY. combined_map already contains every
                # curated isoform independently genome-mapped, so each isoform gets
                # its own true codon coordinates here. We deliberately do NOT run
                # expand_protein_position_to_isoforms (homology transfer): it is
                # redundant when all isoforms are already in combined_map, and in
                # low-complexity/repeat regions its naive context find() collapses
                # many source residues onto one target residue — attributing dozens
                # of foreign genomic positions to a single residue (25-100x row
                # inflation). Direct mapping caps a residue at 3 codon positions x
                # alt alleles.
                for pid, prot_pos, map_aa in entries:
                    if (aaref and aaref != "." and not validate_hgvsp_aa(
                            f"p.{aaref}{prot_pos}", map_aa, prot_pos)):
                        continue
                    key = (pid, prot_pos)
                    if key in seen:
                        continue
                    seen.add(key)
                    out.write("\t".join(
                        [pid, str(prot_pos), chrom, pos_s, end_s] + file_vals
                    ) + "\n")
                    n_written += 1
                    proteins_seen.add(pid)
    finally:
        try:
            fh.close()
        except Exception:
            pass
        if proc is not None:
            proc.wait()
    return n_written, len(proteins_seen)


def _write_empty(out_path: Path) -> None:
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("\t".join(_SYNTH_COLS) + "\n")


def main():
    global _g_pid_map, _g_gene_to_rows, _g_pid_to_seq, _g_pid_to_gene
    global _g_chr_to_pids, _g_protein_ids, _g_chromosomes, _g_bed_header

    p = argparse.ArgumentParser()
    p.add_argument("--seq_table", required=True)
    p.add_argument("--combined_map", required=True)
    p.add_argument("--dbnsfp_raw_dir", required=True)
    p.add_argument("--dbnsfp_bed_header", default=None,
                   help="Path to dbNSFP_custom.bed.header (score column names)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--n_cpu", type=int, default=1,
                   help="Number of parallel chromosome workers (default: 1)")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "dbnsfp_scores.tsv"
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
    chr_to_pids = _chr_to_pids_from_seq(seq_df)

    if chr_to_pids:
        total_pids = sum(len(v) for v in chr_to_pids.values())
        log.info("Per-chromosome protein filter built: %d chromosomes, %d protein-chr assignments",
                 len(chr_to_pids), total_pids)

    gene_to_rows, pid_to_seq, pid_to_gene, _ = load_gene_isoform_lookup(args.seq_table)

    log.info("Loading combined_map …")
    pid_map = load_combined_map_by_protein(args.combined_map)

    # Fast path — single self-describing merged dbNSFP gzip (dbNSFP 5.x academic
    # release). Stream once against an inverted (chr, pos) index; keep all
    # scores + rankscores + CADD + conservation + gnomAD 4.1 joint AF.
    is_merged_gz = (
        raw_dir.is_file()
        and str(raw_dir).endswith(".gz")
        and not args.dbnsfp_bed_header
    )
    if is_merged_gz:
        log.info("Merged dbNSFP file detected — building inverted (chr,pos) index …")
        pid_to_chrom = _pid_to_chrom_from_seq(seq_df)
        gpos_index = build_gpos_index(pid_map, pid_to_chrom)
        log.info("Index built: %d (chr,pos) keys — streaming %s …",
                 len(gpos_index), raw_dir.name)
        n, n_prot = stream_merged_dbnsfp(raw_dir, gpos_index, out)
        log.info("dbNSFP (merged map): %d rows for %d proteins", n, n_prot)
        return

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

    if inputs:
        n_workers = min(max(1, args.n_cpu), len(inputs))
        log.info("Processing %d chromosome files with %d workers …", len(inputs), n_workers)

        if n_workers > 1:
            # Set module globals before forking — child processes inherit via COW,
            # so large dicts are shared without pickling.
            _g_pid_map = pid_map
            _g_gene_to_rows = gene_to_rows
            _g_pid_to_seq = pid_to_seq
            _g_pid_to_gene = pid_to_gene
            _g_chr_to_pids = chr_to_pids
            _g_protein_ids = protein_ids
            _g_chromosomes = chromosomes
            _g_bed_header = header

            ctx = mp.get_context("fork")
            with ctx.Pool(processes=n_workers) as pool:
                chunks = pool.map(_process_one_file, [str(gz) for gz in inputs])
            for chunk in chunks:
                all_rows.extend(chunk)
        else:
            for gz in inputs:
                chrom = _chrom_from_path(gz)
                pids = chr_to_pids.get(chrom, protein_ids) if (chrom and chr_to_pids) else protein_ids
                if gz.suffix == ".bed":
                    log.info("Scanning BED %s (%d proteins) …", gz.name, len(pids))
                    all_rows.extend(map_dbnsfp_bed(
                        gz, pids, pid_map, chromosomes, header,
                        gene_to_rows, pid_to_seq, pid_to_gene))
                else:
                    log.info("Scanning %s (%d proteins) …", gz.name, len(pids))
                    all_rows.extend(map_dbnsfp_file(
                        gz, pids, pid_map, gene_to_rows, pid_to_seq, pid_to_gene))

    if not all_rows and not inputs and raw_dir.is_dir():
        log.info("No dbNSFP chr files for run chromosomes — empty output")

    if all_rows:
        df = pd.DataFrame(all_rows)
        keep = [c for c in _KEEP_COLS if c in df.columns]
        df = df[keep]
    else:
        df = empty

    df.to_csv(out, sep="\t", index=False)
    log.info("dbNSFP (raw map): %d rows for %d proteins",
             len(df), df["Protein_ID"].nunique() if len(df) else 0)


if __name__ == "__main__":
    main()
