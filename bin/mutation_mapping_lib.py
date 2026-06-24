"""
Shared mutation / variant mapping utilities for Module 4 and Module 8f.

Used by create_mutation_map_worker.py, create_dbnsfp_map_worker.py, and tests.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path

HGVSP_RE = re.compile(r"p\.([A-Z*])(\d+)", re.IGNORECASE)


def open_text(path: str | Path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, encoding="utf-8")


def parse_hgvsp_ref_pos(hgvs: str) -> tuple[str | None, int | None]:
    """Return (ref_aa, 1-based position) from HGVSp short notation, or (None, None)."""
    if not hgvs:
        return None, None
    m = HGVSP_RE.search(hgvs.replace("Ter", "*"))
    if not m:
        return None, None
    return m.group(1).upper(), int(m.group(2))


def validate_hgvsp_aa(hgvs: str, mapped_aa: str, protein_pos_1based: int) -> bool:
    """
    Return True if hgvs is empty/frameshift OR ref AA and position match the map.
    Legacy parity: drop missense rows where HGVSp ref != combined_map AA.
    """
    ref_aa, pos = parse_hgvsp_ref_pos(hgvs)
    if ref_aa is None:
        return True
    if pos != protein_pos_1based:
        return False
    map_aa = (mapped_aa or "").upper()
    if ref_aa == "*":
        return True
    return ref_aa == map_aa


def normalize_tcga_sample(sample: str, source: str) -> str:
    """Truncate TCGA barcodes to 12 characters (legacy parity)."""
    if not sample:
        return sample
    src = (source or "").upper()
    if src == "TCGA" and sample.upper().startswith("TCGA-"):
        return sample[:12]
    return sample


def filter_hypermutated_samples(
    variants: list[dict],
    threshold: int = 1500,
) -> list[dict]:
    """Drop variants from samples with more than *threshold* mutations."""
    if threshold <= 0:
        return variants
    counts: dict[str, int] = defaultdict(int)
    for v in variants:
        s = v.get("sample") or ""
        if s:
            counts[s] += 1
    bad = {s for s, n in counts.items() if n > threshold}
    if not bad:
        return variants
    return [v for v in variants if (v.get("sample") or "") not in bad]


def load_combined_map_by_chrom(map_path: str | Path) -> dict:
    """
    Parse combined_map.map → chrom → {genomic_pos → [(enst_id, protein_pos_0based, aa)]}.
    """
    lookup: dict = defaultdict(lambda: defaultdict(list))
    transcript_id = None
    chrom = None

    with open_text(map_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.split()
                transcript_id = parts[1].split("|")[0] if len(parts) > 2 else None
                chrom = parts[2] if len(parts) > 2 else None
                continue
            if transcript_id is None or chrom is None:
                continue
            cols = line.split()
            if len(cols) < 8:
                continue
            try:
                protein_pos = int(cols[0])
                aa = cols[1]
                gpos_str = cols[5]
                for gp in gpos_str.rstrip(",").split(","):
                    if gp and gp != "-":
                        lookup[chrom][int(gp)].append((transcript_id, protein_pos, aa))
            except (ValueError, IndexError):
                pass
    return lookup


NONCODING_VC = frozenset({
    "intron", "5'utr", "3'utr", "silent", "rna", "igr", "3'flank", "5'flank",
    "splice_site", "splice_region", "translation_start_site",
})


def is_protein_coding_variant(variant_type: str, description: str = "") -> bool:
    """Return False for intronic / UTR / silent variants that should not map to protein tables."""
    text = f"{variant_type or ''} {description or ''}".lower().replace("_", " ")
    return not any(k in text for k in NONCODING_VC)


def load_gene_isoform_lookup(loc_path: str | Path):
    """
    Returns (gene_to_rows, pid_to_seq, pid_to_gene, gene_col).
    gene_to_rows: gene → [(Entry_Isoform, Protein_ID, sequence)]
    """
    import pandas as pd

    df = pd.read_csv(loc_path, sep="\t", dtype=str).fillna("")
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"] if c in df.columns), None)
    gene_to_rows: dict[str, list] = {}
    pid_to_seq: dict[str, str] = {}
    pid_to_gene: dict[str, str] = {}

    for _, row in df.iterrows():
        acc = str(row.get("Entry_Isoform", ""))
        pid = str(row.get("Protein_ID", ""))
        seq = str(row.get("Sequence", ""))
        gene = str(row.get(gene_col, "")) if gene_col else ""
        if pid and seq and seq not in ("nan", ""):
            pid_to_seq[pid] = seq
        if pid and gene:
            pid_to_gene[pid] = gene
        entry = (acc, pid, seq if seq not in ("nan", "") else "")
        if gene:
            gene_to_rows.setdefault(gene, [])
            if entry not in gene_to_rows[gene]:
                gene_to_rows[gene].append(entry)
    return gene_to_rows, pid_to_seq, pid_to_gene, gene_col


def main_isoform_protein_id(protein_id: str) -> str:
    """DisCanVis convention: canonical transcript is usually GENE-201."""
    if not protein_id or "-" not in protein_id:
        return protein_id
    gene = protein_id.rsplit("-", 1)[0]
    return f"{gene}-201"


def expand_protein_position_to_isoforms(
    primary_pid: str,
    protein_pos_1based: int,
    gene: str,
    gene_to_rows: dict,
    pid_to_seq: dict,
) -> list[tuple[str, int, bool]]:
    """Translate a primary isoform position to all isoforms of the same gene."""
    out: list[tuple[str, int, bool]] = [(primary_pid, protein_pos_1based, False)]
    seq = pid_to_seq.get(primary_pid, "")
    if not seq or not gene:
        return out

    ctx_s = max(0, protein_pos_1based - 1 - 2)
    ctx_e = min(len(seq), protein_pos_1based - 1 + 1)
    context = seq[ctx_s:ctx_e]
    offset = (protein_pos_1based - 1) - ctx_s
    if not context:
        return out

    for _acc, tgt_pid, tgt_seq in gene_to_rows.get(gene, []):
        if not tgt_seq or tgt_pid == primary_pid:
            continue
        c_idx = tgt_seq.find(context)
        if c_idx == -1:
            continue
        new_pos = c_idx + offset + 1
        if 1 <= new_pos <= len(tgt_seq):
            out.append((tgt_pid, new_pos, True))
    return out


def load_combined_map_by_protein(map_path: str | Path) -> dict[str, dict[int, tuple[int, str]]]:
    """
    Parse combined_map.map → Protein_ID → {genomic_pos → (protein_pos_1based, aa)}.
    """
    result: dict[str, dict[int, tuple[int, str]]] = defaultdict(dict)
    current_pid: str | None = None

    with open_text(map_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.split()
                if len(parts) < 2:
                    current_pid = None
                    continue
                fasta_id = parts[1]
                pid_parts = fasta_id.split("|")
                current_pid = pid_parts[4] if len(pid_parts) > 4 else None
                continue
            if not current_pid:
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            try:
                protein_pos_0 = int(cols[0])
                aa = cols[1]
                gpos_str = cols[5]
                for gp in gpos_str.rstrip(",").split(","):
                    if gp and gp != "-":
                        result[current_pid][int(gp)] = (protein_pos_0 + 1, aa)
            except (ValueError, IndexError):
                pass
    return result
