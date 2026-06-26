#!/usr/bin/env python3
"""
bin/generate_manifest.py — Scan references/ and produce a MANIFEST.tsv
summarising every cached reference file: source, filename, size, modification
date, and SHA-256 checksum.

The manifest is a plain TSV that can be included in mapping reports or checked
into version control to document exactly which data version was used for a run.

Usage:
  python bin/generate_manifest.py                          # write references/MANIFEST.tsv
  python bin/generate_manifest.py --ref_dir /path/to/refs  # custom ref dir
  python bin/generate_manifest.py --out report.tsv          # custom output path
  python bin/generate_manifest.py --no_checksum             # skip slow SHA-256
"""

import argparse
import csv
import hashlib
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

# Files/dirs to skip (work products, stubs, internal Nextflow metadata)
_SKIP_DIRS = {"_stub"}
_SKIP_SUFFIXES = {".idx", ".fai", ".nhr", ".nin", ".nsq", ".nsi", ".nsd", ".nog",
                  ".nto", ".not", ".ntf", ".nal"}
_MAX_CHECKSUM_MB = 500      # files larger than this get checksum skipped by default

# Known source labels (prefix of path relative to ref_dir)
_SOURCE_MAP = {
    "uniprot": "UniProt SwissProt",
    "uniprot_parsed": "UniProt SwissProt (parsed)",
    "gencode": "GENCODE v44",
    "clinvar": "ClinVar",
    "mobidb": "MobiDB",
    "go": "Gene Ontology (GOA)",
    "mondo": "MONDO disease ontology",
    "alphamissense": "AlphaMissense",
    "alphamissense_parsed": "AlphaMissense (decompressed)",
    "ppi": "PPI (IntAct / BioGRID / HIPPIE)",
    "sifts": "SIFTS (EBI)",
    "mavedb": "MaveDB",
    "proteingym": "ProteinGym",
    "depmap": "DepMap",
    "omim": "OMIM (humsavar)",
    "interpro": "InterPro / Pfam",
    "hg38": "hg38 2bit genome",
    "dbsnp": "dbSNP 155 Common bigBed",
    "alphafold": "AlphaFold pLDDT (EBI)",
    "elm": "ELM instances (local copy)",
    "dbnsfp": "dbNSFP",
    "cbioportal": "cBioPortal (local)",
}


def sha256_file(path: Path, max_mb: int) -> str:
    size_mb = path.stat().st_size / (1024 * 1024)
    if size_mb > max_mb:
        return f"(skipped >{max_mb} MB)"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect(ref_dir: Path, no_checksum: bool, max_checksum_mb: int) -> list[dict]:
    rows = []
    for root, dirs, files in os.walk(ref_dir):
        # Skip stub and hidden dirs in-place
        dirs[:] = [d for d in sorted(dirs) if d not in _SKIP_DIRS and not d.startswith(".")]
        root_path = Path(root)
        rel_root = root_path.relative_to(ref_dir)
        source_key = str(rel_root).split(os.sep)[0] if str(rel_root) != "." else ""
        source_label = _SOURCE_MAP.get(source_key, source_key or "other")

        for fname in sorted(files):
            fpath = root_path / fname
            if fpath.suffix in _SKIP_SUFFIXES:
                continue
            try:
                stat = fpath.stat()
            except OSError:
                continue
            size_bytes = stat.st_size
            mod_time = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M")
            checksum = ""
            if not no_checksum and size_bytes > 0:
                try:
                    checksum = sha256_file(fpath, max_checksum_mb)
                except Exception as exc:
                    checksum = f"(error: {exc})"

            rows.append({
                "source":      source_label,
                "file":        str(fpath.relative_to(ref_dir)),
                "size_bytes":  size_bytes,
                "size_human":  _human(size_bytes),
                "modified":    mod_time,
                "sha256":      checksum,
            })
    return rows


def _human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def main():
    ap = argparse.ArgumentParser(description="Generate references/MANIFEST.tsv")
    ap.add_argument("--ref_dir", default="references",
                    help="references directory (default: references/)")
    ap.add_argument("--out", default="",
                    help="output TSV path (default: <ref_dir>/MANIFEST.tsv)")
    ap.add_argument("--no_checksum", action="store_true",
                    help="skip SHA-256 checksums (faster)")
    ap.add_argument("--max_checksum_mb", type=int, default=_MAX_CHECKSUM_MB,
                    help=f"skip checksum for files > N MB (default: {_MAX_CHECKSUM_MB})")
    args = ap.parse_args()

    ref_dir = Path(args.ref_dir)
    if not ref_dir.exists():
        log.error("references dir not found: %s", ref_dir)
        sys.exit(1)

    out_path = Path(args.out) if args.out else ref_dir / "MANIFEST.tsv"

    log.info("Scanning %s ...", ref_dir)
    rows = collect(ref_dir, args.no_checksum, args.max_checksum_mb)

    cols = ["source", "file", "size_bytes", "size_human", "modified", "sha256"]
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    total_bytes = sum(r["size_bytes"] for r in rows)
    log.info("Manifest: %d files, total %s → %s", len(rows), _human(total_bytes), out_path)


if __name__ == "__main__":
    main()
