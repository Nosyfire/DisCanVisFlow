#!/usr/bin/env python3
"""
Mapping report — per-run annotation coverage audit.

Outputs:
  * ``mapping_summary.md``  — reproducibility header + run-wide coverage table
  * ``mapping_coverage.tsv`` — flat TSV: one row per (Gene, annotation),
    with isoform counts + annotation row counts (replaces per-gene MD files
    for large runs; machine-readable and avoids writing ~20k files).
  * Per-gene ``<GENE>_mapping_report.md`` — only when gene count ≤
    ``--per_gene_md_threshold`` (default 50) so single/small runs still get
    the rich Markdown, while full-proteome runs skip it.

Usage (driven by modules/annotation_mapping.nf :: MAPPING_REPORT):
  create_mapping_report_worker.py
      --seq_table loc_chrom_with_names_isoforms_with_seq.tsv
      --final_dir /abs/.../final
      --intermediate_dir /abs/.../intermediate
      --outdir .
      [--mapping_mode all_isoform_mapping]
      [--command "nextflow run ..."] [--pipeline_version 0.5.0]
      [--nextflow_version 26.04.3] [--profile test_one_protein,conda]
      [--run_name foo] [--start_time ...] [--work_dir ...] [--launch_dir ...]
      [--versions_file versions.txt]
      [--source "ELM motifs=local|/abs/path" ...]
      [--per_gene_md_threshold 50]
"""

import argparse
import logging
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

# ── Transcript-mapped annotation config (for the before→after detail) ─────────
ANNOTATION_CONFIG = [
    {"key": "elm",  "label": "ELM motifs", "raw": "elm.tsv", "mapped": "annotations/elm.tsv",
     "source_col": "Entry_Isoform", "id_cols": ["ELMIdentifier", "Start", "End"]},
    {"key": "dibs", "label": "DIBS sites", "raw": "dibs.tsv", "mapped": "annotations/dibs.tsv",
     "source_col": "Entry_Isoform", "id_cols": ["Name", "Start", "End"]},
    {"key": "mfib", "label": "MFIB sites", "raw": "mfib.tsv", "mapped": "annotations/mfib.tsv",
     "source_col": "Entry_Isoform", "id_cols": ["Name", "Start", "End"]},
    {"key": "ptm",  "label": "PTM sites", "raw": "ptm_merged.tsv", "mapped": "annotations/ptm_merged.tsv",
     "source_col": "Entry_Isoform", "id_cols": ["Type", "Position", "Database"]},
    {"key": "uniprot_roi", "label": "UniProt regions of interest",
     "raw": "uniprot_roi.tsv", "mapped": "annotations/uniprot_roi.tsv",
     "source_col": "Accession", "id_cols": ["Type", "Start", "End", "Note"]},
    {"key": "uniprot_binding", "label": "UniProt binding sites",
     "raw": "uniprot_binding.tsv", "mapped": "annotations/uniprot_binding.tsv",
     "source_col": "Accession", "id_cols": ["Type", "Start", "End"]},
    {"key": "pfam", "label": "Pfam domains", "raw": "pfam_domains.tsv", "mapped": "annotations/pfam_domains.tsv",
     "source_col": "Protein_ID", "id_cols": ["hmm_name", "envelope_start", "envelope_end"]},
]

# ── Friendly labels + grouping for the coverage section. Anything not listed is
#    auto-discovered and labelled by its filename. ``unit`` is just display text.
COVERAGE_LABELS = {
    # annotations
    "annotations/elm.tsv": ("ELM motifs", "rows"),
    "annotations/dibs.tsv": ("DIBS sites", "rows"),
    "annotations/mfib.tsv": ("MFIB sites", "rows"),
    "annotations/phasepro.tsv": ("PhasePro entries", "rows"),
    "annotations/ptm_merged.tsv": ("PTM sites", "rows"),
    "annotations/uniprot_roi.tsv": ("UniProt regions of interest", "rows"),
    "annotations/uniprot_binding.tsv": ("UniProt binding sites", "rows"),
    "annotations/pfam_domains.tsv": ("Pfam domains", "rows"),
    "annotations/go_terms.tsv": ("GO terms", "rows"),
    "annotations/polymorphism.tsv": ("Polymorphisms (rsid + allele freq)", "rows"),
    "annotations/pem_core_motifs.tsv": ("PEM core motifs", "rows"),
    "annotations/pem_core_motifs_mapped.tsv": ("PEM core motifs (isoform transfer)", "rows"),
    "annotations/coiled_coils.tsv": ("Coiled coils (DeepCoil)", "isoforms scored"),
    "annotations/interactions.tsv": ("PPI interactions", "pairs"),
    "annotations/scansite.tsv": ("ScanSite motifs", "rows"),
    "annotations/elmswitches_mapped.tsv": ("ELM molecular switches", "rows"),
    "annotations/homology_similarity_manifest.tsv": ("Homology-transfer manifest", "rows"),
    # disorder
    "disorder/IUPredscores.tsv": ("IUPred3 disorder", "isoforms scored"),
    "disorder/Anchorscores.tsv": ("ANCHOR2 binding", "isoforms scored"),
    "disorder/AIUPredscores.tsv": ("AIUPred disorder", "isoforms scored"),
    "disorder/AIUPredBinding.tsv": ("AIUPred-Binding", "isoforms scored"),
    "disorder/AlphaFoldTable.tsv": ("AlphaFold pLDDT", "isoforms scored"),
    "disorder/rsa_scores.tsv": ("RSA (from pLDDT)", "isoforms scored"),
    "disorder/CombinedDisorderNew.tsv": ("Combined disorder regions", "regions"),
    "disorder/CombinedDisorderNew_Pos.tsv": ("Combined disorder (per-residue)", "positions"),
    # pdb
    "pdb/pdb_structures.tsv": ("PDB structures", "rows"),
    "pdb/pdb_missing.tsv": ("PDB missing residues (disorder)", "rows"),
    # conservation
    "conservation/conservation_multiple_level.tsv": ("GOPHER conservation (multi-level)", "rows"),
    "conservation/conservation_phastcons.tsv": ("phastCons conservation", "isoforms scored"),
    # position
    "position/position_based_annotations.tsv": ("Position-based annotations", "positions"),
    # genome
    "genome/exon.tsv": ("Exon boundaries", "exons"),
    "genome/genome_protein_index.tsv": ("Genome↔protein index", "nucleotides"),
    "genome/genome_protein_mutations.tsv": ("All-SNV reference table", "variants"),
    "genome/combined_map.map": ("Genome map (per-residue)", "residues"),
    # disease
    "disease/clinvar_disease.tsv": ("ClinVar disease ontology", "rows"),
    "disease/clinvar_disease_mutations.tsv": ("ClinVar disease mutations", "rows"),
    "disease/omim_disease.tsv": ("OMIM disease ontology", "rows"),
    "disease/omim_mutations.tsv": ("OMIM mutations", "rows"),
    # drivers
    "drivers/cancer_driver.tsv": ("Cancer driver (combined)", "rows"),
    "drivers/census_driver.tsv": ("CGC census driver", "rows"),
    "drivers/compendium_driver.tsv": ("Compendium driver", "rows"),
    # pathogenicity
    "pathogenicity/pathogenicity_scores.tsv": ("dbNSFP pathogenicity scores", "variants"),
    "pathogenicity/alphamissense.tsv": ("AlphaMissense", "variants"),
    "pathogenicity/mavedb.tsv": ("MaveDB functional scores", "variants"),
    "pathogenicity/proteingym.tsv": ("ProteinGym DMS scores", "variants"),
}

# Files that are shared lookups / not per-protein — reported run-wide, not per gene.
SHARED_FILES = {
    "annotations/elm_classes.tsv": ("ELM class definitions (lookup)", "classes"),
    "annotations/transcript_map_stats.tsv": None,   # internal stats, skip
}

# Provenance of EVERY output — (kind, origin). kind ∈ local / downloaded /
# computed / derived. ``origin`` is a relative path (when the data comes from
# another directory in the repo) or a short description (API / tool / derivation).
# Param-driven external file paths are overridden at runtime via --source.
SOURCE_REGISTRY = {
    "annotations/elm.tsv": ("local", "legacy_data/elm/elm_instances-2023.tsv"),
    "annotations/dibs.tsv": ("local", "legacy_data/dibs/dibs_parsed.tsv"),
    "annotations/mfib.tsv": ("local", "legacy_data/mfib/mfib_parsed.tsv"),
    "annotations/phasepro.tsv": ("local", "legacy_data/phasepro/phasepro_parsed.tsv"),
    "annotations/ptm_merged.tsv": ("local", "legacy_data/ptm/{ptmdb,ptmphs}"),
    "annotations/pfam_domains.tsv": ("downloaded", "InterPro/Pfam REST API"),
    "annotations/uniprot_roi.tsv": ("downloaded", "UniProt SwissProt"),
    "annotations/uniprot_binding.tsv": ("downloaded", "UniProt SwissProt"),
    "annotations/go_terms.tsv": ("downloaded", "GOA goa_human.gaf + go.obo"),
    "annotations/polymorphism.tsv": ("downloaded", "dbSNP dbSnp155Common.bb"),
    "annotations/pem_core_motifs.tsv": ("local", "PEM predicted_elm_dataset"),
    "annotations/pem_core_motifs_mapped.tsv": ("derived", "isoform homology transfer"),
    "annotations/coiled_coils.tsv": ("computed", "DeepCoil"),
    "annotations/interactions.tsv": ("local", "IntAct + BioGRID + HIPPIE"),
    "annotations/scansite.tsv": ("downloaded", "MIT ScanSite 4.0"),
    "annotations/elmswitches_mapped.tsv": ("local", "legacy_data/elm/elmswitches-2023.tsv"),
    "annotations/elm_classes.tsv": ("local", "legacy_data/elm/elm_classes-2025.tsv"),
    "annotations/homology_similarity_manifest.tsv": ("derived", "transfer audit"),
    "disorder/IUPredscores.tsv": ("computed", "IUPred3"),
    "disorder/Anchorscores.tsv": ("computed", "ANCHOR2"),
    "disorder/AIUPredscores.tsv": ("computed", "AIUPred"),
    "disorder/AIUPredBinding.tsv": ("computed", "AIUPred-Binding"),
    "disorder/AlphaFoldTable.tsv": ("downloaded", "AlphaFold EBI API (pLDDT)"),
    "disorder/rsa_scores.tsv": ("derived", "from pLDDT"),
    "disorder/CombinedDisorderNew.tsv": ("derived", "MobiDB + RSA + IUPred + Pfam"),
    "disorder/CombinedDisorderNew_Pos.tsv": ("derived", "MobiDB + RSA + IUPred + Pfam"),
    "pdb/pdb_structures.tsv": ("downloaded", "PDBe API"),
    "pdb/pdb_missing.tsv": ("downloaded", "PDBe API"),
    "conservation/conservation_multiple_level.tsv": ("local", "GOPHER conservation_table"),
    "conservation/conservation_phastcons.tsv": ("local", "phastCons bigWig"),
    "position/position_based_annotations.tsv": ("derived", "aggregated per-residue"),
    "genome/combined_map.map": ("computed", "BLAT vs hg38.2bit"),
    "genome/exon.tsv": ("derived", "from combined_map"),
    "genome/genome_protein_index.tsv": ("derived", "from combined_map"),
    "genome/genome_protein_mutations.tsv": ("derived", "all-SNV reference"),
    "disease/clinvar_disease.tsv": ("derived", "MONDO OBO + ClinVar mutations"),
    "disease/clinvar_disease_mutations.tsv": ("derived", "MONDO OBO + ClinVar mutations"),
    "disease/omim_disease.tsv": ("local", "OMIM"),
    "disease/omim_mutations.tsv": ("local", "OMIM"),
    "drivers/cancer_driver.tsv": ("local", "legacy_data/drivers/cancer_driver.tsv"),
    "drivers/census_driver.tsv": ("local", "legacy_data/drivers/census_roles.tsv"),
    "drivers/compendium_driver.tsv": ("local", "legacy_data/drivers/compendium_roles.tsv"),
    "pathogenicity/pathogenicity_scores.tsv": ("local", "dbNSFP"),
    "pathogenicity/alphamissense.tsv": ("downloaded", "AlphaMissense (Google)"),
    "pathogenicity/mavedb.tsv": ("local", "MaveDB"),
    "pathogenicity/proteingym.tsv": ("local", "ProteinGym"),
}
MUTATION_SOURCE = {
    "ClinVar": ("downloaded", "NCBI ClinVar VCF"),
    "TCGA": ("local", "TCGA MAF"),
    "CBioportal": ("local", "cBioPortal MAF"),
    "DepMap": ("local", "DepMap TSV"),
}


def source_of(rel: str, overrides: dict):
    """Return (kind, origin) for an output rel-path, applying --source overrides."""
    if rel in overrides:
        return overrides[rel]
    parts = rel.split("/")
    if parts[0] == "mutations" and len(parts) == 3:
        return MUTATION_SOURCE.get(parts[1], ("local", parts[1]))
    return SOURCE_REGISTRY.get(rel, ("derived", "pipeline output"))


def source_disp(rel: str, overrides: dict) -> str:
    kind, origin = source_of(rel, overrides)
    return f"{kind}: {origin}" if origin else kind


def source_cols(rel: str, overrides: dict):
    """Return (source, source_type) for the two-column rendering."""
    kind, origin = source_of(rel, overrides)
    return (origin or "—", kind)


# Key column resolution: file rel-path → column holding the Protein_ID.
KEY_OVERRIDES = {
    "annotations/interactions.tsv": "Protein_ID_A",
}
CATEGORY_ORDER = ["annotations", "disorder", "pdb", "conservation", "position",
                  "genome", "mutations", "disease", "drivers", "pathogenicity",
                  "sequence"]


# ──────────────────────────────────────────────────────────────────────────
def _read(path: Path, **kw) -> pd.DataFrame:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, **kw).fillna("")
    except Exception as exc:                                    # pragma: no cover
        log.warning("Could not read %s: %s", path, exc)
        return pd.DataFrame()


def _id_of(row, id_cols):
    parts = [str(row[c]) for c in id_cols if c in row.index and str(row[c]) != ""]
    return " | ".join(parts) if parts else "(unnamed)"


def parse_sources(source_args, launch_dir=None):
    """Parse --source overrides ('rel=kind|location'). An absolute ``location``
    is made relative to ``launch_dir`` (we don't need full paths)."""
    import os
    out = {}
    for s in source_args or []:
        key, _, rest = s.partition("=")
        kind, _, loc = rest.partition("|")
        loc = loc.strip()
        if launch_dir and os.path.isabs(loc):
            try:
                loc = os.path.relpath(loc, launch_dir)
            except ValueError:                                  # pragma: no cover
                pass
        out[key.strip()] = (kind.strip() or "unknown", loc)
    return out


def parse_map_regions(combined_map: Path) -> dict:
    """combined_map.map header → {Protein_ID: (chrom, strand, start, end)}."""
    regions = {}
    if combined_map is None or not combined_map.exists():
        return regions
    for line in combined_map.open():
        if not line.startswith("#"):
            continue
        m = re.search(r"(chr[\w.]+)\s+([+-])\s+(\d+)-(\d+)", line)
        fields = line.lstrip("# ").split("|")
        if m and len(fields) >= 5:
            regions[fields[4]] = (m.group(1), m.group(2), m.group(3), m.group(4))
    return regions


def _detect_key_col(rel: str, columns) -> str:
    if rel in KEY_OVERRIDES:
        return KEY_OVERRIDES[rel] if KEY_OVERRIDES[rel] in columns else None
    for cand in ("Protein_ID", "Protein ID", "Transcript name"):
        if cand in columns:
            return cand
    return None


def _category(rel: str) -> str:
    return rel.split("/")[0] if "/" in rel else "other"


def build_coverage(final_dir: Path, pid2gene: dict):
    """Scan every file under final_dir. Return:
       coverage[gene][rel] = {'rows': n, 'isoforms': set(pids), 'counts': {pid:n}}
       shared[rel] = total_rows ; meta[rel] = (label, unit, category, keyed)
    """
    coverage, shared, meta, combined_regions = {}, {}, {}, {}
    files = sorted([p for p in final_dir.rglob("*")
                    if p.is_file() and p.suffix in (".tsv", ".map")])
    for fp in files:
        rel = str(fp.relative_to(final_dir))
        if rel in SHARED_FILES and SHARED_FILES[rel] is None:
            continue
        cat = _category(rel)
        # Internal / non-annotation outputs: sequence plumbing + *_stats tables.
        if cat == "sequence" or fp.name.endswith("_stats.tsv"):
            continue
        # Mutation files share basenames across sources → fold source into label.
        parts = rel.split("/")
        if cat == "mutations" and len(parts) == 3:
            mtype = parts[2].split("_")[0]
            label, unit = (f"{parts[1]}: {mtype} mutations", "rows")
        else:
            label, unit = COVERAGE_LABELS.get(rel, (parts[-1], "rows"))

        if rel == "genome/combined_map.map":
            regs = parse_map_regions(fp)
            combined_regions.update(regs)
            meta[rel] = (label, "isoforms", cat, True)
            for pid in regs:
                g = pid2gene.get(pid, "UNKNOWN")
                cov = coverage.setdefault(g, {}).setdefault(
                    rel, {"rows": 0, "isoforms": set(), "counts": {}})
                cov["rows"] += 1
                cov["isoforms"].add(pid)
                cov["counts"][pid] = cov["counts"].get(pid, 0) + 1
            continue

        # peek header
        try:
            head = pd.read_csv(fp, sep="\t", nrows=0)
        except Exception:
            continue
        key = _detect_key_col(rel, head.columns)
        if key is None:
            # shared/lookup table — count rows run-wide
            n = sum(1 for _ in fp.open()) - 1
            shared[rel] = max(n, 0)
            meta[rel] = (label, unit, cat, False)
            continue

        meta[rel] = (label, unit, cat, True)
        ser = _read(fp, usecols=[key])
        if ser.empty:
            continue
        col = ser[key].astype(str)
        genes = col.map(lambda v: pid2gene.get(v, "UNKNOWN"))
        for g, idx in col.groupby(genes).groups.items():
            sub = col.loc[idx]
            vc = sub.value_counts()
            cov = coverage.setdefault(g, {}).setdefault(
                rel, {"rows": 0, "isoforms": set(), "counts": {}})
            cov["rows"] += len(sub)
            for pid, n in vc.items():
                cov["isoforms"].add(pid)
                cov["counts"][pid] = cov["counts"].get(pid, 0) + int(n)
    return coverage, shared, meta, combined_regions


# ── per-gene report sections ─────────────────────────────────────────────────
def section_sources(intermediate_dir: Path, overrides: dict, gene_pids: set) -> str:
    raw_anno = intermediate_dir / "annotations"
    lines = ["## 2. Annotation sources (before mapping)\n",
             "_Raw transcript-mapping inputs, run-wide. Source from the registry._\n",
             "| Annotation | Source | Source Type | Raw rows | Source isoforms | Distinct identifiers |",
             "|---|---|---|---:|---:|---:|"]
    for cfg in ANNOTATION_CONFIG:
        raw = _read(raw_anno / cfg["raw"])
        if raw.empty:
            n_rows = n_src = n_ids = 0
        else:
            n_rows = len(raw)
            sc = cfg["source_col"] if cfg["source_col"] in raw.columns else (
                "Protein_ID" if "Protein_ID" in raw.columns else None)
            n_src = raw[sc].nunique() if sc else 0
            n_ids = raw.apply(lambda r: _id_of(r, cfg["id_cols"]), axis=1).nunique()
        src, stype = source_cols(cfg["mapped"], overrides)
        lines.append(f"| {cfg['label']} | {src} | {stype} "
                     f"| {n_rows} | {n_src} | {n_ids} |")
    return "\n".join(lines) + "\n"


def section_isoforms(seq_df: pd.DataFrame, regions: dict, gene: str) -> str:
    lines = ["## 3. Isoforms, alignment & genomic locations\n",
             "| Protein_ID | UniProt isoform | Main | Coverage | Alignment | Chr | Genomic region | Genome-mapped |",
             "|---|---|:---:|---:|---|---|---|:---:|"]
    sub = seq_df[seq_df["_gene"] == gene]
    no_genome = []
    for _, r in sub.sort_values("_main", ascending=False).iterrows():
        pid = r["_pid"]
        reg = regions.get(pid)
        if reg:
            chrom, strand, start, end = reg
            region, mapped, chr_disp = f"{strand}{start}-{end}", "✅", chrom
        else:
            region, mapped, chr_disp = "—", "❌", r["_chr"]
            no_genome.append(pid)
        main_flag = "**yes**" if str(r["_main"]).lower() in ("yes", "true", "1") else "no"
        lines.append(f"| {pid} | {r['_iso']} | {main_flag} | {r['_cov']} | {r['_aln']} "
                     f"| {chr_disp} | {region} | {mapped} |")
    lines.append("")
    if no_genome:
        lines.append(f"> ⚠️ **{len(no_genome)} isoform(s) without a genomic location**: "
                     + ", ".join(no_genome))
    else:
        lines.append("> All selected isoforms received a genomic location.")
    return "\n".join(lines) + "\n"


def section_coverage(coverage: dict, meta: dict, gene: str, gene_isos: list,
                     overrides: dict) -> str:
    n_iso = len(gene_isos)
    iso_set = set(gene_isos)
    out = ["## 4. Annotation coverage (all final outputs)\n",
           f"For every output: its **source**, rows belonging to **{gene}**, how "
           f"many of its {n_iso} isoforms carry data, and which isoforms are "
           "missing it (so you can spot any process that produced no data for an "
           "isoform).\n"]
    cov = coverage.get(gene, {})
    by_cat = {}
    for rel, (label, unit, cat, keyed) in meta.items():
        by_cat.setdefault(cat, []).append((rel, label, unit, keyed))
    cats = [c for c in CATEGORY_ORDER if c in by_cat] + \
           [c for c in by_cat if c not in CATEGORY_ORDER]
    for cat in cats:
        out.append(f"### {cat}\n")
        out.append("| Output | Source | Source Type | Rows | Isoforms with data | Missing isoforms |")
        out.append("|---|---|---|---:|:---:|---|")
        for rel, label, unit, keyed in sorted(by_cat[cat]):
            src, stype = source_cols(rel, overrides)
            c = cov.get(rel)
            if not keyed:
                out.append(f"| {label} | {src} | {stype} | _run-wide lookup_ | — | — |")
                continue
            if not c or not c["isoforms"]:
                out.append(f"| {label} | {src} | {stype} | 0 | 0/{n_iso} | ⚠️ all ({n_iso}) |")
                continue
            present = c["isoforms"] & iso_set
            missing = sorted(iso_set - present)
            miss_disp = "—" if not missing else (
                ", ".join(missing) if len(missing) <= 6 else f"{len(missing)} isoforms")
            out.append(f"| {label} | {src} | {stype} | {c['rows']} ({unit}) | "
                       f"{len(present)}/{n_iso} | {miss_disp} |")
        out.append("")
    return "\n".join(out) + "\n"


def section_before_after(intermediate_dir: Path, final_dir: Path, gene_pids: set) -> str:
    raw_dir = intermediate_dir / "annotations"
    out = ["## 5. Per-isoform mapping detail (before → after)\n",
           "For transcript-mapped annotations: **N of M** source annotations "
           "landed on each transcript, by source isoform and `direct` vs "
           "`homology_similarity`, listing which identifiers mapped / did not.\n"]
    for cfg in ANNOTATION_CONFIG:
        raw = _read(raw_dir / cfg["raw"])
        mapped = _read(final_dir / cfg["mapped"])
        out.append(f"### {cfg['label']}\n")
        if mapped.empty or "Protein_ID" not in mapped.columns:
            out.append("_No mapped output for this gene._\n")
            continue
        mapped = mapped[mapped["Protein_ID"].isin(gene_pids)]
        if mapped.empty:
            out.append("_No mapped output for this gene._\n")
            continue
        raw_ids_by_src = {}
        if not raw.empty:
            sc = cfg["source_col"] if cfg["source_col"] in raw.columns else "Protein_ID"
            rsrc = raw[sc].astype(str) if sc in raw.columns else pd.Series([""] * len(raw))
            for src_val, grp in raw.groupby(rsrc):
                raw_ids_by_src[src_val] = set(grp.apply(lambda r: _id_of(r, cfg["id_cols"]), axis=1))
        sc = cfg["source_col"] if cfg["source_col"] in mapped.columns else "Protein_ID"
        mapped = mapped.assign(_src=mapped[sc].astype(str) if sc in mapped.columns else "")
        has_mt = "mapping_type" in mapped.columns
        for pid, grp in mapped.groupby("Protein_ID"):
            grp = grp.assign(_id=grp.apply(lambda r: _id_of(r, cfg["id_cols"]), axis=1))
            n_distinct = grp["_id"].nunique()
            if has_mt:
                direct_ids = set(grp.loc[grp["mapping_type"] == "direct", "_id"])
                homol_only = set(grp.loc[grp["mapping_type"] == "homology_similarity", "_id"]) - direct_ids
                tag = []
                if direct_ids:
                    tag.append(f"{len(direct_ids)} direct")
                if homol_only:
                    tag.append(f"{len(homol_only)} via homology")
                tag_s = f" ({', '.join(tag)})" if tag else ""
            else:
                tag_s = ""
            dup = len(grp) - n_distinct
            dup_s = f" — {len(grp)} rows incl. {dup} redundant homology copies" if dup > 0 else ""
            out.append(f"- **{pid}** — {n_distinct} distinct annotations mapped{tag_s}{dup_s}")
            for src_val, sgrp in grp.groupby("_src"):
                mapped_ids = set(sgrp["_id"])
                raw_ids = raw_ids_by_src.get(src_val, set())
                m_total = len(raw_ids) if raw_ids else len(mapped_ids)
                out.append(f"    - from `{src_val}`: **{len(mapped_ids)} of {m_total}**")
                shown = sorted(mapped_ids)
                out.append(f"        - mapped: {', '.join(shown[:20])}"
                           + (" …" if len(shown) > 20 else ""))
                missing = sorted(raw_ids - mapped_ids)
                if missing:
                    out.append(f"        - NOT mapped ({len(missing)}): "
                               + ", ".join(missing[:20]) + (" …" if len(missing) > 20 else ""))
        out.append("")
    return "\n".join(out) + "\n"


def _write_coverage_tsv(outdir: Path, genes, seq_df, coverage, meta, regions,
                         main_pids, nonmain_pids, overrides):
    """Write mapping_coverage.tsv — flat per-(Gene, annotation) coverage table."""
    rows = []
    # Gene-level metadata
    gene_meta = {}
    for _, r in seq_df.iterrows():
        g = r.get("_gene", "")
        if g not in gene_meta:
            gene_meta[g] = {"n_iso": 0, "n_main": 0, "n_genome": 0, "main_pid": ""}
        gene_meta[g]["n_iso"] += 1
        pid = r.get("_pid", "")
        is_main = str(r.get("_main", "")).lower() in ("yes", "true", "1")
        if is_main:
            gene_meta[g]["n_main"] += 1
            gene_meta[g]["main_pid"] = pid
        if pid in regions:
            gene_meta[g]["n_genome"] += 1

    cats = [c for c in CATEGORY_ORDER if any(
        meta.get(rel, (None, None, None, None))[2] == c for rel in meta)] + \
           [c for c in {meta[r][2] for r in meta if meta[r][3]} if c not in CATEGORY_ORDER]

    for gene in genes:
        gm = gene_meta.get(gene, {"n_iso": 0, "n_main": 0, "n_genome": 0, "main_pid": ""})
        gene_pids = set(seq_df.loc[seq_df["_gene"] == gene, "_pid"])
        gene_main = gene_pids & main_pids
        gene_nonmain = gene_pids - main_pids
        for cat in cats:
            for rel, (label, unit, rcat, keyed) in meta.items():
                if not keyed or rcat != cat:
                    continue
                cov = coverage.get(gene, {}).get(rel, {})
                counts = cov.get("counts", {})
                src, stype = source_cols(rel, overrides)
                m_iso = sum(1 for p in gene_main if counts.get(p, 0) > 0)
                n_iso = sum(1 for p in gene_nonmain if counts.get(p, 0) > 0)
                m_ann = sum(v for p, v in counts.items() if p in gene_main)
                n_ann = sum(v for p, v in counts.items() if p in gene_nonmain)
                rows.append({
                    "Gene": gene,
                    "N_isoforms": gm["n_iso"],
                    "N_main_isoforms": gm["n_main"],
                    "N_genome_mapped": gm["n_genome"],
                    "Main_isoform": gm["main_pid"],
                    "Category": cat,
                    "Annotation": label,
                    "Source": src,
                    "Source_type": stype,
                    "Main_iso_with_data": m_iso,
                    "NonMain_iso_with_data": n_iso,
                    "Ann_main": m_ann,
                    "Ann_nonmain": n_ann,
                })

    df = pd.DataFrame(rows)
    out = outdir / "mapping_coverage.tsv"
    df.to_csv(out, sep="\t", index=False)
    log.info("Wrote mapping_coverage.tsv (%d rows, %d genes × %d annotation tracks)",
             len(df), len(genes), len(df) // max(len(genes), 1))


# ── summary report ───────────────────────────────────────────────────────────
def build_summary(args, seq_df, coverage, meta, regions, genes,
                  overrides, versions, main_pids, nonmain_pids, *, wrote_per_gene=True):
    final_dir = Path(args.final_dir).resolve()
    base = final_dir.name  # 'final'
    L = []
    L.append("# Mapping summary\n")
    L.append(f"_Generated: {datetime.now().isoformat(timespec='seconds')}_\n")

    # ── reproducibility ──
    L.append("## Reproducibility\n")
    L.append("| Field | Value |")
    L.append("|---|---|")
    L.append(f"| Run command | `{args.command or 'n/a'}` |")
    L.append(f"| Pipeline version | {args.pipeline_version or 'n/a'} |")
    L.append(f"| Nextflow version | {args.nextflow_version or 'n/a'} |")
    L.append(f"| Profile | {args.profile or 'n/a'} |")
    L.append(f"| Run name | {args.run_name or 'n/a'} |")
    L.append(f"| Started | {args.start_time or 'n/a'} |")
    L.append(f"| Mapping mode | {args.mapping_mode or 'n/a'} |")
    L.append(f"| Launch dir | `{args.launch_dir or 'n/a'}` |")
    L.append("")
    if versions:
        L.append("### Tool / data versions\n")
        L.append("| Component | Version |")
        L.append("|---|---|")
        for k, v in versions:
            L.append(f"| {k} | {v} |")
        L.append("")

    # ── output file locations (relative to the run output dir) ──
    L.append("## Output locations\n")
    L.append(f"All outputs are under `{final_dir}`. Paths below are relative to it.\n")
    L.append("| Category | Path | Files |")
    L.append("|---|---|---|")
    for d in sorted([p for p in final_dir.iterdir() if p.is_dir()]):
        fnames = sorted(p.name for p in d.iterdir() if p.is_file())
        subs = sorted(p.name for p in d.iterdir() if p.is_dir())
        extra = f" (+ {', '.join(subs)})" if subs else ""
        L.append(f"| {d.name} | `{base}/{d.name}/` | {', '.join(fnames) or '—'}{extra} |")
    L.append("")

    # ── run-wide overview: every annotation, its source, and main / non-main
    #    isoform coverage + annotation counts ──
    n_main, n_non = len(main_pids), len(nonmain_pids)
    gm_main = sum(1 for p in main_pids if p in regions)
    gm_non = sum(1 for p in nonmain_pids if p in regions)
    L.append("## Mapping overview (all annotations)\n")
    L.append(f"- **Genes / proteins:** {len(genes)}")
    L.append(f"- **Isoforms:** {n_main + n_non}  (main: {n_main}, non-main: {n_non})")
    L.append(f"- **Genome-mapped isoforms:** {gm_main + gm_non}  "
             f"(main: {gm_main}, non-main: {gm_non})")
    L.append("")
    L.append("Each output is an independent row. *Main isoforms w/ data* and "
             "*Non-main isoforms w/ data* count how many isoforms carry the "
             "annotation; *Ann. (main)* / *Ann. (non-main)* are the actual "
             "annotation row counts on those isoforms (run-wide).\n")
    L.append("| Category | Annotation | Source | Source Type | Main iso w/ data | Non-main iso w/ data | Ann. (main) | Ann. (non-main) |")
    L.append("|---|---|---|---|:---:|:---:|---:|---:|")

    # aggregate per-pid counts across all genes for each rel
    agg = {}  # rel -> {pid: count}
    for g, relmap in coverage.items():
        for rel, c in relmap.items():
            d = agg.setdefault(rel, {})
            for pid, n in c["counts"].items():
                d[pid] = d.get(pid, 0) + n

    by_cat = {}
    for rel, (label, unit, cat, keyed) in meta.items():
        if not keyed:
            continue
        by_cat.setdefault(cat, []).append((rel, label))
    cats = [c for c in CATEGORY_ORDER if c in by_cat] + \
           [c for c in by_cat if c not in CATEGORY_ORDER]
    for cat in cats:
        for rel, label in sorted(by_cat[cat]):
            counts = agg.get(rel, {})
            m_iso = sum(1 for p in main_pids if counts.get(p, 0) > 0)
            n_iso = sum(1 for p in nonmain_pids if counts.get(p, 0) > 0)
            m_ann = sum(v for p, v in counts.items() if p in main_pids)
            n_ann = sum(v for p, v in counts.items() if p in nonmain_pids)
            src, stype = source_cols(rel, overrides)
            L.append(f"| {cat} | {label} | {src} | {stype} "
                     f"| {m_iso} / {n_main} | {n_iso} / {n_non} | {m_ann} | {n_ann} |")
    L.append("")
    if wrote_per_gene:
        L.append("_Per-protein detailed reports are written alongside this summary as "
                 "`<GENE>_mapping_report.md` (one per gene)._")
    else:
        L.append("_Run contained too many genes for per-gene Markdown reports. "
                 "Machine-readable per-gene annotation coverage is in "
                 "`mapping_coverage.tsv` (columns: Gene, Category, Annotation, "
                 "Source, Source_type, Main_iso_with_data, NonMain_iso_with_data, "
                 "Ann_main, Ann_nonmain)._")
    return "\n".join(L) + "\n"


def load_versions(path):
    if not path:
        return []
    p = Path(path)
    if not p.exists():
        return []
    out = []
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        k, _, v = line.partition(":")
        out.append((k.strip(), v.strip()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq_table", required=True)
    ap.add_argument("--final_dir", required=True)
    ap.add_argument("--intermediate_dir", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mapping_mode", default="")
    ap.add_argument("--command", default="")
    ap.add_argument("--pipeline_version", default="")
    ap.add_argument("--nextflow_version", default="")
    ap.add_argument("--profile", default="")
    ap.add_argument("--run_name", default="")
    ap.add_argument("--start_time", default="")
    ap.add_argument("--work_dir", default="")
    ap.add_argument("--launch_dir", default="")
    ap.add_argument("--versions_file", default=None)
    ap.add_argument("--source", action="append", default=[])
    ap.add_argument("--per_gene_md_threshold", type=int, default=50,
                    help="Write per-gene MD only when gene count <= this (default 50). "
                         "Full-proteome runs use mapping_coverage.tsv instead.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    final_dir = Path(args.final_dir)
    intermediate_dir = Path(args.intermediate_dir) if args.intermediate_dir else final_dir.parent / "intermediate"

    seq_df = _read(Path(args.seq_table))
    if seq_df.empty:
        log.error("Empty sequence table — cannot build report")
        (outdir / "mapping_summary.md").write_text("# Mapping summary\n\n_No sequence data._\n")
        return

    def col(df, *names, default=""):
        for n in names:
            if n in df.columns:
                return df[n]
        return pd.Series([default] * len(df), index=df.index)

    seq_df = seq_df.copy()
    seq_df["_pid"] = col(seq_df, "Protein_ID", "Transcript name")
    seq_df["_iso"] = col(seq_df, "Entry_Isoform")
    seq_df["_gene"] = col(seq_df, "Gene", "Gene_Gencode", "Gene_Uniprot")
    seq_df["_main"] = col(seq_df, "main_isoform")
    seq_df["_cov"] = col(seq_df, "coverage")
    seq_df["_aln"] = col(seq_df, "alignmentpuntcuality")
    seq_df["_chr"] = col(seq_df, "Chromosome")
    pid2gene = dict(zip(seq_df["_pid"], seq_df["_gene"]))

    overrides = parse_sources(args.source, args.launch_dir)
    versions = load_versions(args.versions_file)
    coverage, shared, meta, regions = build_coverage(final_dir, pid2gene)

    is_main = seq_df["_main"].astype(str).str.lower().isin(["yes", "true", "1"])
    main_pids = set(seq_df.loc[is_main, "_pid"])
    nonmain_pids = set(seq_df.loc[~is_main, "_pid"])

    genes = sorted([g for g in seq_df["_gene"].unique() if g])
    log.info("Building reports for %d gene(s): %s", len(genes), ", ".join(genes[:20]))

    write_per_gene = len(genes) <= args.per_gene_md_threshold

    if write_per_gene:
        for gene in genes:
            sub = seq_df[seq_df["_gene"] == gene]
            gene_pids = set(sub["_pid"])
            gene_isos = list(sub["_pid"])
            n_iso = len(gene_isos)
            main_pid = ""
            mains = sub[sub["_main"].astype(str).str.lower().isin(["yes", "true", "1"])]
            if not mains.empty:
                main_pid = mains.iloc[0]["_pid"]
            n_gen = sum(1 for p in gene_pids if p in regions)

            head = [
                f"# Mapping report — {gene}\n",
                f"_Generated: {datetime.now().isoformat(timespec='seconds')}_\n",
                "## 1. Overview\n",
                f"- **Gene:** {gene}",
                f"- **Mapping mode:** {args.mapping_mode or 'n/a'}",
                f"- **Isoforms selected:** {n_iso}",
                f"- **Main isoform (transcript):** {main_pid or 'n/a'}",
                f"- **Isoforms with a genomic location:** {n_gen} / {n_iso}",
                f"- **Pipeline / Nextflow:** {args.pipeline_version or '?'} / {args.nextflow_version or '?'}",
                "",
            ]
            report = "\n".join(head) + "\n"
            if intermediate_dir.exists():
                report += section_sources(intermediate_dir, overrides, gene_pids) + "\n"
            report += section_isoforms(seq_df, regions, gene) + "\n"
            report += section_coverage(coverage, meta, gene, gene_isos, overrides) + "\n"
            if intermediate_dir.exists():
                report += section_before_after(intermediate_dir, final_dir, gene_pids)

            (outdir / f"{gene}_mapping_report.md").write_text(report)
        log.info("Wrote %d per-gene report(s)", len(genes))
    else:
        log.info("Gene count %d > threshold %d — writing mapping_coverage.tsv instead of per-gene MD files",
                 len(genes), args.per_gene_md_threshold)
        _write_coverage_tsv(outdir, genes, seq_df, coverage, meta, regions,
                            main_pids, nonmain_pids, overrides)

    summary = build_summary(args, seq_df, coverage, meta, regions,
                            genes, overrides, versions, main_pids, nonmain_pids,
                            wrote_per_gene=write_per_gene)
    (outdir / "mapping_summary.md").write_text(summary)
    log.info("Wrote mapping_summary.md")


if __name__ == "__main__":
    main()
