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

# Known filenames across DepMap releases (tried in order; first match wins)
_MUTATION_FILENAMES = [
    "OmicsSomaticMutations.csv",             # 22Q2 – 24Q1
    "OmicsSomaticMutationsProfile.csv",      # 24Q2+
    "OmicsSomaticMutationsMatrixDamaging.csv",
    "CCLE_mutations.csv",                    # pre-22Q2 legacy
]

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
    """Query the DepMap catalogue CSV → presigned URL for a somatic-mutations file.

    Tries `filename` first, then each entry in _MUTATION_FILENAMES, then falls
    back to any catalogue row whose filename contains both 'somatic' and
    'mutation' (case-insensitive). Takes the newest (lexicographically greatest)
    release unless `release` is specified."""
    r = _req(sess, catalogue_url)
    if r is None:
        return None
    rows = list(csv.DictReader(io.StringIO(r.text)))

    def _pick(candidates):
        if not candidates:
            return None
        if release:
            candidates = ([row for row in candidates
                           if (row.get("release") or "").strip() == release]
                          or candidates)
        candidates.sort(key=lambda row: (row.get("release") or ""))
        chosen = candidates[-1]
        log.info("DepMap: using '%s' from release '%s'",
                 chosen.get("filename", "?"), chosen.get("release", "?"))
        return (chosen.get("url") or "").strip()

    # 1. exact filename requested by caller
    cand = [row for row in rows if (row.get("filename") or "").strip() == filename]
    if cand:
        return _pick(cand)

    # 2. try known aliases in priority order
    for known in _MUTATION_FILENAMES:
        if known == filename:
            continue
        cand = [row for row in rows if (row.get("filename") or "").strip() == known]
        if cand:
            log.warning("'%s' not found — using '%s' instead", filename, known)
            return _pick(cand)

    # 3. fuzzy: any file whose name contains 'somatic' and 'mutation'
    cand = [row for row in rows
            if all(k in (row.get("filename") or "").lower()
                   for k in ("somatic", "mutation"))]
    if cand:
        names = sorted(set(r.get("filename", "") for r in cand))
        log.warning("'%s' not found — fuzzy match: %s", filename, names)
        return _pick(cand)

    # nothing matched — log everything to help diagnose future renames
    all_files = sorted(set((row.get("filename") or "").strip() for row in rows))
    omics = [f for f in all_files if "omics" in f.lower() or "mutation" in f.lower()]
    log.error("no somatic-mutation file found in DepMap catalogue (%d rows)", len(rows))
    log.error("mutation-related files available: %s", omics or all_files[:30])
    return None


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

    sess = _session()

    # 1. obtain the raw mutations CSV (local path, direct URL, or via catalogue)
    if args.file_url and Path(args.file_url).exists():
        raw_csv = Path(args.file_url)
    else:
        # Check whether any known cached file already exists (filename may differ from release to release)
        raw_csv = next(
            (cache_dir / fn for fn in _MUTATION_FILENAMES
             if (cache_dir / fn).exists() and (cache_dir / fn).stat().st_size > 0),
            cache_dir / args.filename,
        )
        if not raw_csv.exists() or raw_csv.stat().st_size == 0:
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
