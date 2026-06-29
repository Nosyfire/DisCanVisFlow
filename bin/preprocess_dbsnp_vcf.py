#!/usr/bin/env python3
"""
preprocess_dbsnp_vcf.py — Extract compact MAF table from NCBI dbSNP VCF.

Input:  GCF_000001405.40.gz  (NCBI dbSNP latest release, hg38, ~28 GiB)
Output: dbsnp_maf.tsv.gz     (bgzipped TSV, ready for tabix indexing)

Output columns (tab-separated, sorted by chrom+pos):
  chrom | pos | rsid | ref | alt | maf | is_common

- chrom:     chr1..chr22, chrX, chrY, chrM
- pos:       1-based (VCF convention)
- rsid:      rs123456
- ref/alt:   single-base SNVs only (multi-allelic and indels skipped)
- maf:       max per-population minor allele frequency from FREQ INFO tag
- is_common: 1 if COMMON flag is set in INFO (MAF ≥ 1% in ≥ 1 population)

Only variants with COMMON flag are written to keep the output compact.

Usage:
  preprocess_dbsnp_vcf.py --vcf GCF_000001405.40.gz --out dbsnp_maf.tsv.gz

After this script, tabix-index the output:
  tabix -s 1 -b 2 -e 2 dbsnp_maf.tsv.gz
"""

import argparse
import gzip
import logging
import re
import subprocess
import sys
import tempfile
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

# RefSeq accession → UCSC chromosome name (GRCh38 / hg38)
_NC_TO_CHR = {
    "NC_000001.11": "chr1",  "NC_000002.12": "chr2",  "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",  "NC_000005.10": "chr5",  "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",  "NC_000008.11": "chr8",  "NC_000009.12": "chr9",
    "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9":  "chr14", "NC_000015.10": "chr15",
    "NC_000016.10": "chr16", "NC_000017.11": "chr17", "NC_000018.10": "chr18",
    "NC_000019.10": "chr19", "NC_000020.11": "chr20", "NC_000021.9":  "chr21",
    "NC_000022.11": "chr22", "NC_000023.11": "chrX",  "NC_000024.10": "chrY",
    "NC_012920.1":  "chrM",
}

_FREQ_RE = re.compile(r"FREQ=([^;]+)")


def _parse_maf(info: str) -> str:
    """Extract max per-population minor allele frequency from FREQ= INFO tag.
    FREQ format: Pop1:ref_freq,alt_freq|Pop2:ref_freq,alt_freq...
    Returns empty string if no frequency data is available."""
    m = _FREQ_RE.search(info)
    if not m:
        return ""
    best = None
    for pop_entry in m.group(1).split("|"):
        _, _, freqs = pop_entry.partition(":")
        parts = freqs.split(",")
        if len(parts) < 2:
            continue
        try:
            ref_f = float(parts[0])
            alt_f = float(parts[1])
        except ValueError:
            continue
        maf = min(ref_f, alt_f)
        if best is None or maf > best:
            best = maf
    return "" if best is None else f"{best:.6g}"


def _is_common(info: str) -> bool:
    """Return True if the COMMON flag is present in INFO."""
    for field in info.split(";"):
        if field == "COMMON":
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="NCBI dbSNP VCF (.gz)")
    ap.add_argument("--out", required=True, help="Output bgzipped TSV path")
    args = ap.parse_args()

    vcf_path = Path(args.vcf)
    out_path = Path(args.out)

    # Write to a plain temp TSV first (sorted output required for tabix)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False,
                                      dir=out_path.parent)
    tmp_path = Path(tmp.name)

    n_written = n_skip_chrom = n_skip_snv = n_skip_nocommon = 0
    log.info("Streaming %s …", vcf_path)

    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if line.startswith("#"):
                continue

            # Fast pre-check: skip lines without COMMON before splitting
            if "COMMON" not in line:
                n_skip_nocommon += 1
                continue

            cols = line.rstrip("\n").split("\t", 8)
            if len(cols) < 8:
                continue

            chrom_nc, pos, rsid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            info = cols[7]

            chrom = _NC_TO_CHR.get(chrom_nc)
            if chrom is None:
                n_skip_chrom += 1
                continue

            # Only biallelic single-base SNVs
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                n_skip_snv += 1
                continue

            if not _is_common(info):
                n_skip_nocommon += 1
                continue

            maf = _parse_maf(info)
            tmp.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t{maf}\t1\n")
            n_written += 1

            if (i + 1) % 10_000_000 == 0:
                log.info("  … %d M lines scanned, %d written", (i + 1) // 1_000_000, n_written)

    tmp.close()
    log.info("Scan complete: %d common SNVs written, %d skipped (non-chr: %d, non-SNV: %d, "
             "non-common: %d)", n_written, n_skip_chrom + n_skip_snv + n_skip_nocommon,
             n_skip_chrom, n_skip_snv, n_skip_nocommon)

    # Sort by chrom+pos then gzip-compress
    log.info("Sorting by chrom+pos …")
    sorted_tmp = tmp_path.with_suffix(".sorted.tsv")
    sort_cmd = ["sort", "-k1,1", "-k2,2n", "-T", str(out_path.parent),
                "--output", str(sorted_tmp), str(tmp_path)]
    subprocess.run(sort_cmd, check=True)
    tmp_path.unlink()

    log.info("Compressing → %s …", out_path)
    header = "chrom\tpos\trsid\tref\talt\tmaf\tis_common\n"
    with open(sorted_tmp, "rb") as src, gzip.open(out_path, "wb") as dst:
        dst.write(header.encode())
        while chunk := src.read(4 * 1024 * 1024):
            dst.write(chunk)
    sorted_tmp.unlink()

    log.info("Done. Output: %s (%.1f MiB)", out_path, out_path.stat().st_size / 1024**2)


if __name__ == "__main__":
    main()
