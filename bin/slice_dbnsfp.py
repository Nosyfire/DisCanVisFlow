#!/usr/bin/env python3
"""
Slice a bgzipped dbNSFP map by Protein_ID without reading the whole file.

The companion packer (`bin/dbnsfp_pack.sh`) produces three files next to each
other, sharing the stem `dbnsfp_scores`:
  * `dbnsfp_scores.tsv.gz`   — BGZF (block-gzip, seekable) body, sorted by
                               Protein_ID then Protein_position
  * `dbnsfp_scores.tsv.gz.gzi` — bgzip random-access index (from `bgzip -r`)
  * `dbnsfp_scores.pidx`     — Protein_ID <tab> uncompressed_offset <tab> length
  * `dbnsfp_scores.header`   — the single header line

This tool looks up each requested Protein_ID in the `.pidx`, then asks `bgzip`
to decompress only that byte range (`bgzip -b <offset> -s <length>`), so a slice
reads a few KB instead of ~170 GB.

Usage:
  slice_dbnsfp.py --bgz results/discanvis/final/pathogenicity/dbnsfp_scores.tsv.gz \
      --id RAF1-201                 # one isoform
  slice_dbnsfp.py --bgz <...> --id RAF1-201,BRAF-201   # several
  slice_dbnsfp.py --bgz <...> --id_file ids.txt --out raf.tsv
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def _stem_paths(bgz: Path):
    """dbnsfp_scores.tsv.gz -> (pidx, header) sharing the 'dbnsfp_scores' stem."""
    name = bgz.name
    if name.endswith(".tsv.gz"):
        stem = name[: -len(".tsv.gz")]
    elif name.endswith(".gz"):
        stem = name[: -len(".gz")]
    else:
        stem = bgz.stem
    return bgz.parent / f"{stem}.pidx", bgz.parent / f"{stem}.header"


def load_pidx(path: Path) -> dict:
    idx = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            idx[parts[0]] = (int(parts[1]), int(parts[2]))
    return idx


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bgz", required=True, help="dbnsfp_scores.tsv.gz (BGZF)")
    ap.add_argument("--id", default="", help="Protein_ID(s), comma-separated")
    ap.add_argument("--id_file", default=None, help="file with one Protein_ID per line")
    ap.add_argument("--out", default=None, help="output TSV (default: stdout)")
    ap.add_argument("--no_header", action="store_true", help="omit the header line")
    args = ap.parse_args()

    if not shutil.which("bgzip"):
        sys.stderr.write("ERROR: bgzip not found on PATH (install htslib)\n")
        return 2

    bgz = Path(args.bgz)
    pidx_path, header_path = _stem_paths(bgz)
    for p in (bgz, pidx_path):
        if not p.is_file():
            sys.stderr.write(f"ERROR: missing {p}\n")
            return 2

    ids: list[str] = []
    if args.id:
        ids += [x.strip() for x in args.id.split(",") if x.strip()]
    if args.id_file:
        ids += [ln.strip() for ln in Path(args.id_file).read_text().splitlines()
                if ln.strip() and not ln.startswith("#")]
    if not ids:
        sys.stderr.write("ERROR: no --id / --id_file given\n")
        return 2

    idx = load_pidx(pidx_path)
    out = open(args.out, "w", encoding="utf-8") if args.out else sys.stdout
    try:
        if not args.no_header and header_path.is_file():
            out.write(header_path.read_text())
        missing = []
        for pid in ids:
            if pid not in idx:
                missing.append(pid)
                continue
            off, length = idx[pid]
            # bgzip random access by uncompressed offset (needs the .gzi index)
            proc = subprocess.run(
                ["bgzip", "-b", str(off), "-s", str(length), str(bgz)],
                stdout=subprocess.PIPE, check=True,
            )
            out.write(proc.stdout.decode("utf-8", "replace"))
        if missing:
            sys.stderr.write(f"WARN: {len(missing)} id(s) not in index: "
                             f"{', '.join(missing[:10])}"
                             f"{' …' if len(missing) > 10 else ''}\n")
    finally:
        if out is not sys.stdout:
            out.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
