#!/usr/bin/env python3
"""
create_mutation_map_worker.py — Module 4: Mutation Mapping

Maps genomic mutations to protein positions on ALL isoforms of each gene.

Approach
--------
1. Parse combined_map.map  → {chrom: {genomic_pos: [(enst_id, protein_pos, aa)]}}
2. Parse loc_chrom TSV     → {enst_base: metadata}, sequence + gene lookup
3. Parse mutation input    → list of variant dicts
4. Primary mapping         → genome position → protein position via combined_map.map
5. Isoform expansion       → translate each primary hit to all isoforms of the
                             same gene using a 3-AA context substring search
6. Classify & write        → Missense / Frameshift / Nonsense / Indel TSVs

Inputs
------
--combined_map  combined_map.map                         (Module 3 output)
--loc_chrom     loc_chrom_with_names_isoforms_with_seq.tsv (Module 2 output)
--clinvar_vcf   clinvar.vcf[.gz]   (mutually exclusive with --maf / --vcf)
--maf           mutations.maf[.gz] / mutations.tsv
--vcf           generic VCF (non-ClinVar)
--source        label: ClinVar | TCGA | CBioportal | custom
--output_dir    output directory (default: .)

Outputs
-------
Missense_filter_mutations_mapped.tsv
Frameshift_filter_mutations_mapped.tsv
Nonsense_filter_mutations_mapped.tsv
Indel_filter_mutations_mapped.tsv
mutation_stats.tsv
"""

import argparse
import gzip
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from mutation_mapping_lib import (
    filter_hypermutated_samples,
    is_protein_coding_variant,
    load_combined_map_by_chrom,
    normalize_tcga_sample,
    validate_hgvsp_aa,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Output columns
# ---------------------------------------------------------------------------
OUT_COLS = [
    "Protein_ID",           # GENCODE transcript name (e.g. RAF1-201)
    "Accession",            # UniProt Entry_Isoform (e.g. P04049)
    "Gene",
    "Mutation Description", # ClinicalSignificance / Variant_Classification
    "Mutation",             # HGVSp_Short
    "Protein_position",
    "Study Abbrevation",
    "Study Name",
    "Sample name",
    "Start_Position",
    "isoform_mapped",       # True when the same genomic variant was mapped onto
                            # a sibling isoform by exact sequence-context match.
                            # NB: this is a coordinate *mapping*, not annotation
                            # homology transfer (which is region-similarity based).
    # ClinVar-specific (empty for non-ClinVar sources)
    "ClinicalSignificance",
    "PhenotypeList",
    "PhenotypeIDS",
    "ReviewStatus",
    "RCVaccession",
    "MONDO_ID",
    "MeSH_ID",
]


# ---------------------------------------------------------------------------
# File helper
# ---------------------------------------------------------------------------

def _open(path: str):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, encoding="utf-8")


# ---------------------------------------------------------------------------
# combined_map.map parser
# ---------------------------------------------------------------------------

def load_combined_map(map_path: str) -> dict:
    """Wrapper around shared parser (tests import this name)."""
    lookup = load_combined_map_by_chrom(map_path)
    log.info("Loaded combined_map: %d chromosomes", len(lookup))
    return lookup


# ---------------------------------------------------------------------------
# loc_chrom loader + isoform lookup builder
# ---------------------------------------------------------------------------

def load_loc_chrom(loc_path: str):
    """
    Returns (df_indexed_by_enst_base, gene_to_rows, pid_to_seq).

    gene_to_rows: gene → [(Entry_Isoform, Protein_ID, sequence)]
    pid_to_seq:   Protein_ID → sequence
    """
    df = pd.read_csv(loc_path, sep="\t", dtype=str).fillna("")

    ts_col = "transcript_stable_id" if "transcript_stable_id" in df.columns else "Transcript ID"
    df["_ts_base"] = df[ts_col].str.split(".").str[0]

    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"] if c in df.columns), None)

    gene_to_rows: dict[str, list] = {}
    pid_to_seq: dict[str, str] = {}

    for _, row in df.iterrows():
        acc  = row.get("Entry_Isoform", "")
        pid  = row.get("Protein_ID", "")
        seq  = row.get("Sequence", "")
        gene = row.get(gene_col, "") if gene_col else ""

        if pid and seq and seq not in ("nan", ""):
            pid_to_seq[pid] = seq

        entry = (acc, pid, seq if seq not in ("nan", "") else "")
        gene_to_rows.setdefault(gene, [])
        if entry not in gene_to_rows[gene]:
            gene_to_rows[gene].append(entry)

    df_idx = df.set_index("_ts_base", drop=False)
    log.info("Loaded loc_chrom: %d isoforms across %d genes",
             len(df_idx), len(gene_to_rows))
    return df_idx, gene_to_rows, pid_to_seq, gene_col


# ---------------------------------------------------------------------------
# ClinVar VCF parser
# ---------------------------------------------------------------------------

CLNHGVS_RE  = re.compile(r"CLNHGVS=([^;]+)")
CLNSIG_RE   = re.compile(r"CLNSIG=([^;]+)")
CLNDN_RE    = re.compile(r"CLNDN=([^;]+)")
CLNVC_RE    = re.compile(r"CLNVC=([^;]+)")
GENEINFO_RE = re.compile(r"GENEINFO=([^;]+)")
CLNREVSTAT_RE = re.compile(r"CLNREVSTAT=([^;]+)")
CLNVI_RE    = re.compile(r"CLNVI=([^;]+)")
RS_RE       = re.compile(r"RS=([^;]+)")
MONDO_RE    = re.compile(r"MONDO:(\d+)")
MESH_RE     = re.compile(r"MeSH:([A-Z]\d+)")


def parse_clinvar_vcf(vcf_path: str) -> list[dict]:
    rows = []
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref  = parts[3]
            alt  = parts[4]
            info = parts[7]

            m_hgvs   = CLNHGVS_RE.search(info)
            m_sig    = CLNSIG_RE.search(info)
            m_dn     = CLNDN_RE.search(info)
            m_vc     = CLNVC_RE.search(info)
            m_gi     = GENEINFO_RE.search(info)
            m_rev    = CLNREVSTAT_RE.search(info)
            m_vi     = CLNVI_RE.search(info)
            m_rs     = RS_RE.search(info)
            mondo_ids = ",".join(MONDO_RE.findall(info))
            mesh_ids  = ",".join(MESH_RE.findall(info))

            clndn = m_dn.group(1).replace("_", " ") if m_dn else ""
            phenotype_ids = m_vi.group(1) if m_vi else ""

            rows.append({
                "chrom":           chrom,
                "pos":             pos,
                "ref":             ref,
                "alt":             alt,
                "hgvs":            m_hgvs.group(1) if m_hgvs else "",
                "clnsig":          m_sig.group(1).replace("_", " ") if m_sig else "",
                "disease":         clndn,
                "variant_type":    m_vc.group(1) if m_vc else "",
                "geneinfo":        m_gi.group(1) if m_gi else "",
                "review_status":   m_rev.group(1).replace("_", " ") if m_rev else "",
                "phenotype_ids":   phenotype_ids,
                "rcv":             m_vi.group(1).split(":")[0] if m_vi else "",
                "rs":              m_rs.group(1) if m_rs else "",
                "mondo_id":        mondo_ids,
                "mesh_id":         mesh_ids,
                "study_abbr":      "ClinVar",
                "study_name":      "ClinVar",
                "sample":          "",
            })
    log.info("Parsed ClinVar VCF: %d variants", len(rows))
    return rows


# ---------------------------------------------------------------------------
# Generic MAF parser
# ---------------------------------------------------------------------------

def parse_maf(
    maf_path: str,
    allowed_chroms: set[str] | None = None,
    allowed_genes: set[str] | None = None,
    chunksize: int = 100_000,
) -> list[dict]:
    """Stream MAF in chunks; restrict to chromosomes / genes present in the run."""
    rows = []
    for chunk in pd.read_csv(
        maf_path, sep="\t", comment="#", dtype=str,
        low_memory=False, chunksize=chunksize,
    ):
        chunk.columns = chunk.columns.str.strip()
        if allowed_chroms:
            chrom_col = chunk.get("Chromosome", pd.Series(dtype=str)).astype(str).str.strip()
            norm = chrom_col.apply(
                lambda c: c if c.startswith("chr") else (f"chr{c}" if c else c)
            )
            chunk = chunk[norm.isin(allowed_chroms)]
        if allowed_genes and "Hugo_Symbol" in chunk.columns:
            chunk = chunk[chunk["Hugo_Symbol"].astype(str).str.strip().isin(allowed_genes)]
        if chunk.empty:
            continue

        for _, r in chunk.iterrows():
            chrom = str(r.get("Chromosome", "")).strip()
            if chrom and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                pos = int(r.get("Start_Position", 0))
            except (ValueError, TypeError):
                pos = 0

            rows.append({
                "chrom":         chrom,
                "pos":           pos,
                "ref":           str(r.get("Reference_Allele", "")),
                "alt":           str(r.get("Tumor_Seq_Allele2", r.get("Tumor_Seq_Allele1", ""))),
                "hgvs":          str(r.get("HGVSp_Short", r.get("HGVSp", ""))),
                "clnsig":        str(r.get("Variant_Classification", "")),
                "disease":       str(r.get("DISEASE", r.get("oncotree_code", ""))),
                "variant_type":  str(r.get("Variant_Type", "")),
                "study_abbr":    str(r.get("Study_Abbr", r.get("STUDY_ABBR", ""))),
                "study_name":    str(r.get("Study_Name", r.get("STUDY_NAME", ""))),
                "sample":        str(r.get("Tumor_Sample_Barcode", r.get("sample_id", ""))).strip(),
                "review_status": "",
                "phenotype_ids": "",
                "rcv":           "",
                "rs":            "",
                "mondo_id":      "",
                "mesh_id":       "",
            })
    log.info("Parsed MAF: %d variants (chrom filter=%s, gene filter=%s)",
             len(rows),
             sorted(allowed_chroms) if allowed_chroms else "none",
             sorted(allowed_genes) if allowed_genes else "none")
    return rows


# ---------------------------------------------------------------------------
# Generic VCF parser (non-ClinVar)
# ---------------------------------------------------------------------------

def parse_generic_vcf(vcf_path: str, source: str) -> list[dict]:
    rows = []
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref = parts[3]
            alt = parts[4].split(",")[0]  # first alt allele

            rows.append({
                "chrom":         chrom,
                "pos":           pos,
                "ref":           ref,
                "alt":           alt,
                "hgvs":          "",
                "clnsig":        "",
                "disease":       "",
                "variant_type":  "",
                "study_abbr":    source,
                "study_name":    source,
                "sample":        "",
                "review_status": "",
                "phenotype_ids": "",
                "rcv":           "",
                "rs":            "",
                "mondo_id":      "",
                "mesh_id":       "",
            })
    log.info("Parsed generic VCF: %d variants", len(rows))
    return rows


# ---------------------------------------------------------------------------
# Mutation classifier
# ---------------------------------------------------------------------------

def _classify(hgvs: str, ref: str, alt: str, variant_type: str, description: str = "") -> str | None:
    """Return mutation class or None if the variant should be dropped (non-coding)."""
    if not is_protein_coding_variant(variant_type, description):
        return None

    vc = (variant_type or "").lower()
    h  = (hgvs or "").lower()
    desc = (description or "").lower()

    if "nonsense" in vc or "stop_gained" in vc or (
        (h.endswith("*") or h.endswith("ter")) and "intron" not in desc and "intron" not in vc
    ):
        return "Nonsense_Mutation"
    if "frame_shift" in vc or "frameshift" in h or "fs" in h:
        kind = "Del" if (ref and len(ref) > len(alt)) else "Ins"
        return f"Frame_Shift_{kind}"
    if any(k in vc for k in ("del", "ins", "indel", "in_frame")):
        return "Indel"
    if len(ref) != len(alt) or len(ref) != 1:
        return "Indel"
    return "Missense_Mutation"


# ---------------------------------------------------------------------------
# Resolve + isoform expansion
# ---------------------------------------------------------------------------

def resolve_mutations(
    variants:      list[dict],
    lookup:        dict,
    loc_df,
    gene_to_rows:  dict,
    pid_to_seq:    dict,
    gene_col:      str | None,
    source:        str,
    map_all_isoforms: bool = True,
    validate_hgvsp: bool = True,
) -> list[dict]:
    results = []
    unmapped = 0
    seen: set = set()  # deduplicate (Protein_ID, pos, hgvs) per variant

    for var in variants:
        chrom = var["chrom"]
        pos   = var["pos"]

        hits = lookup.get(chrom, {}).get(pos, [])
        if not hits:
            unmapped += 1
            continue

        # Collect primary hits
        primary_pids: list[tuple] = []  # (Protein_ID, protein_pos_0based, acc, gene)
        for (ts_id, protein_pos, aa) in hits:
            ts_base = ts_id.split(".")[0]
            if ts_base not in loc_df.index:
                continue
            meta = loc_df.loc[ts_base]
            if isinstance(meta, pd.DataFrame):
                meta = meta.iloc[0]

            pid  = str(meta.get("Protein_ID", "") or "")
            acc  = str(meta.get("Entry_Isoform", meta.get("Entry_Name", "")) or "")
            gene = str(meta.get(gene_col, "") if gene_col else "") or ""
            if not pid:
                continue
            primary_pids.append((pid, protein_pos, acc, gene))

        for (pid, protein_pos, acc, gene) in primary_pids:
            seq = pid_to_seq.get(pid, "")
            prot_pos_1based = protein_pos + 1
            if validate_hgvsp and not validate_hgvsp_aa(
                var.get("hgvs", ""), aa, prot_pos_1based
            ):
                continue

            def _make_row(tgt_pid, tgt_acc, tgt_pos, is_isoform_mapped):
                return {
                    "Protein_ID":           tgt_pid,
                    "Accession":            tgt_acc,
                    "Gene":                 gene,
                    "Mutation Description": var.get("clnsig", var.get("variant_type", "")),
                    "Mutation":             var.get("hgvs", ""),
                    "Protein_position":     tgt_pos,
                    "Study Abbrevation":    var.get("study_abbr", source),
                    "Study Name":           var.get("study_name", source),
                    "Sample name":          var.get("sample", ""),
                    "Start_Position":       pos,
                    "isoform_mapped":       is_isoform_mapped,
                    "ClinicalSignificance": var.get("clnsig", ""),
                    "PhenotypeList":        var.get("disease", ""),
                    "PhenotypeIDS":         var.get("phenotype_ids", ""),
                    "ReviewStatus":         var.get("review_status", ""),
                    "RCVaccession":         var.get("rcv", ""),
                    "MONDO_ID":             var.get("mondo_id", ""),
                    "MeSH_ID":              var.get("mesh_id", ""),
                    "_ref":                 var.get("ref", ""),
                    "_alt":                 var.get("alt", ""),
                    "_variant_type":        var.get("variant_type", ""),
                    "_hgvs":                var.get("hgvs", ""),
                }

            # Primary isoform
            key = (pid, prot_pos_1based, var.get("hgvs", ""), pos)
            if key not in seen:
                seen.add(key)
                results.append(_make_row(pid, acc, prot_pos_1based, False))

            if not map_all_isoforms or not seq:
                continue

            # Sequence-context translation to all isoforms of the same gene
            ctx_s   = max(0, protein_pos - 2)
            ctx_e   = min(len(seq), protein_pos + 1)
            context = seq[ctx_s:ctx_e]
            offset  = protein_pos - ctx_s  # 0-based offset in context

            if not context:
                continue

            for (tgt_acc, tgt_pid, tgt_seq) in gene_to_rows.get(gene, []):
                if not tgt_seq or tgt_pid == pid:
                    continue
                c_idx = tgt_seq.find(context)
                if c_idx == -1:
                    continue
                new_pos = c_idx + offset + 1  # 1-based
                if new_pos < 1 or new_pos > len(tgt_seq):
                    continue
                key2 = (tgt_pid, new_pos, var.get("hgvs", ""), pos)
                if key2 not in seen:
                    seen.add(key2)
                    results.append(_make_row(tgt_pid, tgt_acc, new_pos, True))

    log.info("Resolved: %d rows from %d variants (unmapped: %d)",
             len(results), len(variants), unmapped)
    return results


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_split_tsv(rows: list[dict], outdir: Path, source: str) -> dict:
    missense, frameshifts, nonsense, indels = [], [], [], []

    for r in rows:
        cls = _classify(
            r["_hgvs"], r["_ref"], r["_alt"], r["_variant_type"],
            r.get("Mutation Description", ""),
        )
        if cls is None:
            continue
        out_row = {c: r.get(c, "") for c in OUT_COLS}
        if cls == "Missense_Mutation":
            missense.append(out_row)
        elif cls == "Nonsense_Mutation":
            nonsense.append(out_row)
        elif cls.startswith("Frame_Shift"):
            frameshifts.append(out_row)
        else:
            indels.append(out_row)

    def _write(data, fname):
        df = pd.DataFrame(data, columns=OUT_COLS) if data else pd.DataFrame(columns=OUT_COLS)
        df.to_csv(outdir / fname, sep="\t", index=False)
        log.info("Written: %s (%d rows)", fname, len(df))
        return len(df)

    n_mis = _write(missense,    "Missense_filter_mutations_mapped.tsv")
    n_fs  = _write(frameshifts, "Frameshift_filter_mutations_mapped.tsv")
    n_non = _write(nonsense,    "Nonsense_filter_mutations_mapped.tsv")
    n_ind = _write(indels,      "Indel_filter_mutations_mapped.tsv")

    pd.DataFrame([{
        "source":          source,
        "total_resolved":  len(rows),
        "missense":        n_mis,
        "frameshift":      n_fs,
        "nonsense":        n_non,
        "indel":           n_ind,
    }]).to_csv(outdir / "mutation_stats.tsv", sep="\t", index=False)

    return {"missense": n_mis, "frameshift": n_fs, "nonsense": n_non, "indel": n_ind}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Module 4: map genomic mutations to protein positions (all isoforms)")
    p.add_argument("--combined_map",       required=True)
    p.add_argument("--loc_chrom",          required=True)
    p.add_argument("--clinvar_vcf",        default=None)
    p.add_argument("--maf",                default=None)
    p.add_argument("--vcf",                default=None, help="Generic VCF (non-ClinVar)")
    p.add_argument("--source",             default="ClinVar")
    p.add_argument("--no_isoform_expand",  action="store_true", default=False,
                   help="Skip sequence-based isoform expansion (primary transcript only)")
    p.add_argument("--hypermutation_threshold", type=int, default=1500,
                   help="Drop MAF samples with more than N mutations (0=disable)")
    p.add_argument("--no_hgvsp_validation", action="store_true",
                   help="Skip HGVSp ref-AA validation against combined_map.map")
    p.add_argument("--output_dir",         default=".")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not (args.clinvar_vcf or args.maf or args.vcf):
        log.error("Provide --clinvar_vcf, --maf, or --vcf")
        sys.exit(1)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading combined_map …")
    lookup = load_combined_map(args.combined_map)
    allowed_chroms = set(lookup.keys())

    log.info("Loading loc_chrom …")
    loc_df, gene_to_rows, pid_to_seq, gene_col = load_loc_chrom(args.loc_chrom)
    gene_col_name = gene_col or "Gene_Gencode"
    allowed_genes = set(
        loc_df[gene_col_name].dropna().astype(str).str.strip()
    ) - {"", "nan"}

    log.info("Parsing mutation input …")
    if args.clinvar_vcf:
        variants = parse_clinvar_vcf(args.clinvar_vcf)
        if allowed_chroms:
            variants = [v for v in variants if v["chrom"] in allowed_chroms]
    elif args.maf:
        variants = parse_maf(
            args.maf,
            allowed_chroms=allowed_chroms,
            allowed_genes=allowed_genes if allowed_genes else None,
        )
    else:
        variants = parse_generic_vcf(args.vcf, args.source)
        if allowed_chroms:
            variants = [v for v in variants if v["chrom"] in allowed_chroms]

    if args.maf and args.hypermutation_threshold > 0:
        before = len(variants)
        variants = filter_hypermutated_samples(variants, args.hypermutation_threshold)
        log.info("Hypermutation filter: %d → %d variants", before, len(variants))

    if args.maf and args.source.upper() == "TCGA":
        for v in variants:
            v["sample"] = normalize_tcga_sample(v.get("sample", ""), args.source)

    log.info("Resolving mutations → protein positions (all isoforms=%s) …",
             not args.no_isoform_expand)
    resolved = resolve_mutations(
        variants, lookup, loc_df, gene_to_rows, pid_to_seq,
        gene_col, args.source,
        map_all_isoforms=not args.no_isoform_expand,
        validate_hgvsp=not args.no_hgvsp_validation,
    )

    log.info("Writing output TSVs …")
    counts = write_split_tsv(resolved, outdir, args.source)

    log.info("Done — missense=%d  frameshift=%d  nonsense=%d  indel=%d",
             counts["missense"], counts["frameshift"],
             counts["nonsense"], counts["indel"])


if __name__ == "__main__":
    main()
