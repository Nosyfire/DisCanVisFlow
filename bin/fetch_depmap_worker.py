#!/usr/bin/env python3
"""
fetch_depmap_worker.py — download DepMap somatic mutations (OPEN data) and
normalise them into the TSV that create_depmap_worker.py consumes.

DepMap public data is freely downloadable. The portal exposes a CSV catalogue
at https://depmap.org/portal/api/download/files whose rows are

    release,release_date,filename,url,md5_hash

and the `url` is a *time-limited presigned* Google-Storage link — so we resolve
it fresh at run time rather than hard-coding it. We pick the
`OmicsSomaticMutations.csv` of the newest (or a requested) public release,
download it, and emit a tab-separated table with the columns the mapper wants:

    HugoSymbol  Protein_position  HGVSp_Short  ModelID  Start_Position
    EntrezGeneID  Hotspot

Raw DepMap calls the protein-change column `ProteinChange` (e.g. "p.G12D") —
we split it into HGVSp_Short (the p. string) + Protein_position (the integer).

Usage:
  fetch_depmap_worker.py --out depmap_mutations.tsv
      [--catalogue_url URL] [--file_url URL] [--release "DepMap Public 26Q1"]
      [--filename OmicsSomaticMutations.csv] [--cache_dir DIR] [--limit N]
"""

import argparse
import csv
import io
import logging
import re
import sys
import time
from pathlib import Path

import requests

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

CATALOGUE_URL = "https://depmap.org/portal/api/download/files"

OUT_COLS = ["HugoSymbol", "Protein_position", "HGVSp_Short", "ModelID",
            "Start_Position", "EntrezGeneID", "Hotspot"]

# "p.G12D" / "p.(Gly12Asp)" / "p.Gly12Aspfs" → first integer = 12
_POS_RE = re.compile(r"p\.\(?[A-Za-z*]{1,3}(\d+)")


def pos_from_protein_change(expr):
    """Return the 1-based protein position from a ProteinChange string, or ''."""
    if not expr:
        return ""
    m = _POS_RE.search(str(expr))
    return m.group(1) if m else ""


def _session():
    s = requests.Session()
    s.headers.update({"User-Agent": "discanvis-pipeline/1.0"})
    return s


def _req(sess, url, timeout=300, retries=4, **kw):
    last = None
    for attempt in range(1, retries + 1):
        try:
            r = sess.get(url, timeout=timeout, **kw)
            if r.status_code == 200:
                return r
            last = f"HTTP {r.status_code}"
            if r.status_code in (429, 500, 502, 503, 504):
                wait = min(60, 2 ** attempt)
                log.warning("%s on %s (try %d/%d) — backoff %ss",
                            last, url, attempt, retries, wait)
                time.sleep(wait)
                continue
            log.warning("%s on %s — not retrying", last, url)
            return None
        except requests.RequestException as exc:
            last = str(exc)
            wait = min(60, 2 ** attempt)
            log.warning("%s on %s (try %d/%d) — backoff %ss",
                        exc, url, attempt, retries, wait)
            time.sleep(wait)
    log.error("giving up on %s (%s)", url, last)
    return None


def resolve_file_url(sess, catalogue_url, filename, release):
    """Query the DepMap catalogue CSV → presigned URL for `filename`.

    If `release` is given, match it exactly; otherwise take the lexicographically
    greatest release string (newest YYQn sorts correctly)."""
    r = _req(sess, catalogue_url)
    if r is None:
        return None
    rows = list(csv.DictReader(io.StringIO(r.text)))
    cand = [row for row in rows
            if (row.get("filename") or "").strip() == filename]
    if not cand:
        log.error("no '%s' in DepMap catalogue (%d rows)", filename, len(rows))
        return None
    if release:
        cand = [row for row in cand
                if (row.get("release") or "").strip() == release] or cand
    cand.sort(key=lambda row: (row.get("release") or ""))
    chosen = cand[-1]
    log.info("DepMap: %s from release '%s'", filename, chosen.get("release"))
    return (chosen.get("url") or "").strip()


def normalise(handle, limit=0):
    """Raw DepMap mutations CSV → list of OUT_COLS dicts (missense-style rows
    with a parseable protein change)."""
    rdr = csv.DictReader(handle)
    out = []
    for rec in rdr:
        pchange = (rec.get("ProteinChange") or rec.get("HGVSp_Short")
                   or rec.get("HGVSp") or "").strip()
        pos = pos_from_protein_change(pchange)
        if not pos:
            continue
        out.append({
            "HugoSymbol": (rec.get("HugoSymbol") or rec.get("Hugo_Symbol")
                           or rec.get("Gene") or "").strip(),
            "Protein_position": pos,
            "HGVSp_Short": pchange,
            "ModelID": (rec.get("ModelID") or rec.get("DepMap_ID")
                        or rec.get("model_id") or "").strip(),
            "Start_Position": (rec.get("Pos") or rec.get("Start_Position")
                               or rec.get("Position") or "").strip(),
            "EntrezGeneID": (rec.get("EntrezGeneID") or rec.get("Entrez_Gene_Id")
                             or "").strip(),
            "Hotspot": (rec.get("Hotspot") or rec.get("HessDriver") or "").strip(),
        })
        if limit and len(out) >= limit:
            break
    return out


def main():
    ap = argparse.ArgumentParser(description="Fetch + normalise DepMap mutations")
    ap.add_argument("--out", required=True, help="output normalised TSV")
    ap.add_argument("--catalogue_url", default=CATALOGUE_URL,
                    help="DepMap download catalogue CSV endpoint")
    ap.add_argument("--file_url", default="",
                    help="direct URL/local path to the mutations CSV (skip catalogue)")
    ap.add_argument("--release", default="",
                    help="exact DepMap release name (default: newest)")
    ap.add_argument("--filename", default="OmicsSomaticMutations.csv",
                    help="catalogue filename to resolve")
    ap.add_argument("--cache_dir", default="",
                    help="cache dir for the raw CSV (default: alongside --out)")
    ap.add_argument("--limit", type=int, default=0, help="cap rows (testing)")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir) if args.cache_dir else out_path.parent / "depmap_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    raw_csv = cache_dir / "OmicsSomaticMutations.csv"

    sess = _session()

    # 1. obtain the raw mutations CSV (local path, direct URL, or via catalogue)
    if args.file_url and Path(args.file_url).exists():
        raw_csv = Path(args.file_url)
    elif not raw_csv.exists() or raw_csv.stat().st_size == 0:
        url = args.file_url or resolve_file_url(sess, args.catalogue_url,
                                                args.filename, args.release)
        if not url:
            log.error("could not resolve a DepMap mutations URL — aborting")
            sys.exit(1)
        log.info("downloading DepMap mutations CSV...")
        r = _req(sess, url, stream=True)
        if r is None:
            log.error("download failed — aborting")
            sys.exit(1)
        with open(raw_csv, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        log.info("cached %s (%d bytes)", raw_csv, raw_csv.stat().st_size)

    # 2. normalise columns
    with open(raw_csv, newline="", encoding="utf-8", errors="replace") as fh:
        rows = normalise(fh, limit=args.limit)

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log.info("DepMap: wrote %d normalised rows → %s", len(rows), out_path)


if __name__ == "__main__":
    main()
