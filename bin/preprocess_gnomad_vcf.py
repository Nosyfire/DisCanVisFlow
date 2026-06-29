#!/usr/bin/env python3
"""
preprocess_gnomad_vcf.py — Extract compact AF table from one gnomAD exome VCF.

Called once per chromosome by the FETCH_GNOMAD_VCF Nextflow process.
Streams through a single bgzipped chromosome VCF, emits TSV lines to stdout
(no header — FETCH_GNOMAD_VCF writes the header and merges all chromosomes).

Output columns (written to stdout, tab-separated):
  chrom | pos | rsid | ref | alt | af | af_popmax | is_common

- chrom:      chr1..chr22, chrX, chrY (passed explicitly via --chrom)
- pos:        1-based (VCF convention)
- rsid:       rs ID from VCF ID field, or '.' if absent
- ref/alt:    single-base biallelic SNVs only
- af:         gnomAD global exome AF (AF INFO tag), empty if absent
- af_popmax:  gnomAD max population AF (AF_popmax INFO tag), empty if absent
- is_common:  1 if af >= 0.01 OR af_popmax >= 0.01, else 0

Only PASS variants are emitted.

Usage (called by FETCH_GNOMAD_VCF shell loop, output redirected to a file):
  preprocess_gnomad_vcf.py --vcf gnomad.exomes.v4.1.sites.chr1.vcf.bgz --chrom chr1
"""

import argparse
import gzip
import logging
import re
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_AF_RE    = re.compile(r"(?:^|;)AF=([^;,]+)")
_AFPOP_RE = re.compile(r"(?:^|;)AF_popmax=([^;,]+)")


def _extract_af(info: str, pattern: re.Pattern) -> str:
    m = pattern.search(info)
    if not m:
        return ""
    val = m.group(1).split(",")[0]   # take first value (single-allele for biallelic)
    try:
        return f"{float(val):.6g}"
    except ValueError:
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf",   required=True, help="gnomAD chromosome VCF (.bgz or .gz)")
    ap.add_argument("--chrom", required=True, help="chromosome label (e.g. chr1)")
    args = ap.parse_args()

    vcf_path = Path(args.vcf)
    chrom    = args.chrom

    n_written = n_skip = 0
    log.info("Processing %s → %s …", vcf_path.name, chrom)

    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t", 8)
            if len(cols) < 8:
                continue

            _chrom, pos, rsid, ref, alt, _qual, filt, info = (
                cols[0], cols[1], cols[2], cols[3], cols[4],
                cols[5], cols[6], cols[7])

            # Only PASS variants
            if filt not in ("PASS", "."):
                n_skip += 1
                continue

            # Only biallelic SNVs
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                n_skip += 1
                continue

            af        = _extract_af(info, _AF_RE)
            af_popmax = _extract_af(info, _AFPOP_RE)

            try:
                af_f    = float(af)    if af        else 0.0
                af_pop_f = float(af_popmax) if af_popmax else 0.0
            except ValueError:
                af_f = af_pop_f = 0.0

            is_common = 1 if (af_f >= 0.01 or af_pop_f >= 0.01) else 0

            sys.stdout.write(
                f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t{af}\t{af_popmax}\t{is_common}\n")
            n_written += 1

    log.info("%s: %d variants written, %d skipped", chrom, n_written, n_skip)


if __name__ == "__main__":
    main()
