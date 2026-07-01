#!/usr/bin/env python3
"""
bin/extract_gene_from_results.py — Extract one gene's rows from a completed
full-proteome results directory without re-running the pipeline.

Every annotation TSV in results/<project>/final/ uses Protein_ID as its
primary key (GENCODE transcript names, e.g. RAF1-201).  This script filters
all TSVs to only rows whose Protein_ID starts with <GENE>- and writes the
filtered copies to <out_dir>/final/<category>/<file>.

Usage:
  python bin/extract_gene_from_results.py \
      --source results/discanvis \
      --gene   RAF1 \
      --out    results/discanvis_raf1

  # Multiple genes (comma-separated)
  python bin/extract_gene_from_results.py \
      --source results/discanvis \
      --gene   RAF1,BRAF,KRAS \
      --out    results/discanvis_kinases

  # Gene list from file (one HGNC name per line, # comments OK)
  python bin/extract_gene_from_results.py \
      --source           results/discanvis \
      --gene_list_file   config/gene_lists/cellular_vulnerability.txt \
      --out              results/cellular_vulnerability
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

# Some annotation TSVs have fields >128KB (e.g. long coiled-coil score arrays)
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

# TSV columns that carry a Protein_ID value (checked in order; first match wins)
_PID_COLS = ["Protein_ID", "protein_id", "TranscriptID", "transcript_id"]

# Files in final/ that use gene-name columns instead of Protein_ID — skip them
_SKIP_FILENAMES = {
    "cancer_driver.tsv", "census_driver.tsv", "compendium_driver.tsv",
    "elm_classes.tsv", "homology_similarity_manifest.tsv",
}


def _find_pid_col(header: list[str]) -> str | None:
    for col in _PID_COLS:
        if col in header:
            return col
    return None


def _gene_prefix(gene: str) -> str:
    return f"{gene}-"


def filter_tsv(src: Path, dst: Path, prefixes: list[str]) -> int:
    """Read src TSV, filter rows by Protein_ID prefix, write to dst.
    Returns number of rows written (excluding header)."""
    with open(src, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_text("")
            return 0

        pid_col = _find_pid_col(header)
        if pid_col is None:
            # No Protein_ID column — copy whole file unchanged
            dst.parent.mkdir(parents=True, exist_ok=True)
            import shutil
            shutil.copy2(src, dst)
            return -1

        pid_idx = header.index(pid_col)
        kept = []
        for row in reader:
            if len(row) <= pid_idx:
                continue
            pid = row[pid_idx]
            if any(pid.startswith(p) for p in prefixes):
                kept.append(row)

    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(kept)
    return len(kept)


def extract(source_dir: Path, genes: list[str], out_dir: Path) -> None:
    prefixes = [_gene_prefix(g) for g in genes]
    gene_label = ",".join(genes)
    log.info("Extracting %s from %s → %s", gene_label, source_dir, out_dir)

    final_dir = source_dir / "final"
    if not final_dir.exists():
        log.error("No final/ directory found in %s", source_dir)
        sys.exit(1)

    tsv_files = sorted(final_dir.rglob("*.tsv"))
    if not tsv_files:
        log.error("No TSV files found under %s", final_dir)
        sys.exit(1)

    total_rows = 0
    copied = 0
    skipped = 0
    for src in tsv_files:
        if src.name in _SKIP_FILENAMES:
            skipped += 1
            continue
        rel = src.relative_to(final_dir)
        dst = out_dir / "final" / rel
        n = filter_tsv(src, dst, prefixes)
        if n == -1:
            log.debug("  (no Protein_ID col) copied %s", rel)
            copied += 1
        else:
            log.info("  %s: %d rows", rel, n)
            total_rows += n

    # Copy mapping_reports if present
    report_dir = source_dir / "mapping_reports"
    if report_dir.exists():
        import shutil
        dst_rep = out_dir / "mapping_reports"
        if dst_rep.exists():
            shutil.rmtree(dst_rep)
        shutil.copytree(report_dir, dst_rep)
        log.info("  mapping_reports/ copied")

    log.info("Done: %d TSV rows extracted (%d files copied, %d skipped)",
             total_rows, copied, skipped)


def main():
    ap = argparse.ArgumentParser(description="Extract gene rows from completed results")
    ap.add_argument("--source", required=True,
                    help="source project results dir (e.g. results/discanvis)")
    ap.add_argument("--gene", default="",
                    help="HGNC gene name(s), comma-separated (e.g. RAF1 or RAF1,BRAF)")
    ap.add_argument("--gene_list_file", default="",
                    help="plain-text file: one HGNC name per line, # comments OK")
    ap.add_argument("--out", required=True,
                    help="output dir (e.g. results/discanvis_raf1)")
    args = ap.parse_args()

    genes: list[str] = []
    if args.gene:
        genes += [g.strip().upper() for g in args.gene.split(",") if g.strip()]
    if args.gene_list_file:
        gene_list_path = Path(args.gene_list_file)
        if not gene_list_path.exists():
            log.error("--gene_list_file not found: %s", gene_list_path)
            sys.exit(1)
        for line in gene_list_path.read_text(encoding="utf-8").splitlines():
            line = line.split("#")[0].strip()
            if line:
                genes.append(line.upper())

    if not genes:
        ap.error("Specify at least one gene via --gene or --gene_list_file")

    genes = list(dict.fromkeys(genes))  # deduplicate, preserve order

    source = Path(args.source)
    out = Path(args.out)
    if not source.exists():
        log.error("Source dir does not exist: %s", source)
        sys.exit(1)

    extract(source, genes, out)


if __name__ == "__main__":
    main()
