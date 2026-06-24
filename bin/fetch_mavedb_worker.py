#!/usr/bin/env python3
"""
fetch_mavedb_worker.py — download MaveDB functional scores from the OFFICIAL
MaveDB API (https://api.mavedb.org, docs at https://api.mavedb.org/docs) and
emit a UniProt-keyed raw table that the mapping step turns into Protein_ID rows.

Why the official API (not the Broad g2p mirror): it is the authoritative,
versioned source, exposes every published score set, and carries the UniProt
cross-reference + protein HGVS we need for coordinate mapping.

Flow:
  1. enumerate score sets   POST /api/v1/score-sets/search        {"text": ""}
  2. keep protein targets that have a UniProt external identifier
  3. per score set          GET  /api/v1/score-sets/{urn}/scores  (CSV)
  4. parse hgvs_pro → 1-based position, classify single vs multi-mutant
  5. write mavedb_raw.tsv  (uniprot-keyed; downstream maps to Protein_ID)

Output columns (mavedb_raw.tsv):
  uniprot  gene_name  urn  mavedb_id  prot_expr  protein_start  score  is_double_mutant

Usage:
  fetch_mavedb_worker.py --outdir DIR [--uniprot P04637,P38398]
                         [--urn urn:mavedb:00000059-a-1] [--max_sets N]
                         [--timeout 120] [--retries 4] [--cache_dir DIR]
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

API = "https://api.mavedb.org/api/v1"
OUT_COLS = ["uniprot", "gene_name", "urn", "mavedb_id", "prot_expr",
            "protein_start", "score", "is_double_mutant"]

# p.Lys291Trp / p.(Lys291Trp) / p.Lys291= → capture the residue number
_POS_RE = re.compile(r"p\.\(?[A-Za-z*]{1,3}(\d+)")


def _session():
    s = requests.Session()
    s.headers.update({"User-Agent": "discanvis-pipeline/1.0",
                      "accept": "application/json"})
    return s


def _req(sess, url, timeout, retries, method="GET", **kw):
    """GET/POST with exponential backoff (MaveDB occasionally 502/504s)."""
    last = None
    for attempt in range(1, retries + 1):
        try:
            r = sess.request(method, url, timeout=timeout, **kw)
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


def _uniprot_of(target_gene: dict):
    """Pull the UniProt accession from a targetGene.

    Tries the structured externalIdentifiers first, then falls back to the
    newer ``uniprotIdFromMappedMetadata`` field the API added.
    """
    for ext in target_gene.get("externalIdentifiers", []) or []:
        ident = ext.get("identifier", {}) if isinstance(ext, dict) else {}
        db = (ident.get("dbName") or ident.get("db_name") or "").lower()
        if db == "uniprot":
            acc = ident.get("identifier") or ident.get("id")
            if acc:
                return acc
    return target_gene.get("uniprotIdFromMappedMetadata")


def _search_page(sess, args, text, limit, offset):
    """One page of POST /score-sets/search. Returns a list (possibly empty),
    or None on transport failure. Handles both the wrapped
    ``{"scoreSets": [...], "numScoreSets": N}`` response and a legacy bare list."""
    r = _req(sess, f"{API}/score-sets/search", args.timeout, args.retries,
             method="POST", json={"text": text, "limit": limit, "offset": offset},
             headers={"Content-Type": "application/json"})
    if r is None:
        return None
    j = r.json()
    if isinstance(j, dict):
        return j.get("scoreSets", []) or []
    return j or []


def _iter_score_sets(sess, args, text=""):
    """Paginate the search endpoint via body limit/offset until exhausted."""
    limit, offset = 100, 0
    while True:
        page = _search_page(sess, args, text, limit, offset)
        if not page:                      # None (error) or empty page → stop
            return
        for s in page:
            yield s
        if len(page) < limit:
            return
        offset += limit


def _row_from_score_set(s, uniprot_filter):
    """Return (urn, gene_name, uniprot) for a protein score set with a UniProt
    xref passing the filter, else None. Uses the targetGenes embedded in the
    search result (no extra per-URN GET needed)."""
    urn = s.get("urn")
    if not urn:
        return None
    for tg in s.get("targetGenes", []) or []:
        cat = (tg.get("category") or "").lower()
        if cat not in ("protein_coding", "protein coding", ""):
            continue
        up = _uniprot_of(tg)
        if not up:
            continue
        if uniprot_filter and up.split("-")[0] not in uniprot_filter:
            continue
        return (urn, tg.get("name", ""), up)
    return None


def enumerate_score_sets(sess, args, uniprot_filter):
    """Return [(urn, gene_name, uniprot)] for protein score sets with UniProt.

    Full-scan + in-memory filter (the search endpoint paginates 100/page; there
    are ~2.8k sets, so this is ~28 cached requests and avoids text-match misses).
    """
    out, seen = [], set()

    if args.urn:
        meta = _req(sess, f"{API}/score-sets/{args.urn}", args.timeout, args.retries)
        if meta is not None:
            row = _row_from_score_set(meta.json(), uniprot_filter)
            if row:
                out.append(row)
        return out

    n_scanned = 0
    for s in _iter_score_sets(sess, args, ""):
        n_scanned += 1
        row = _row_from_score_set(s, uniprot_filter)
        if not row or row[0] in seen:
            continue
        seen.add(row[0])
        out.append(row)
        if args.max_sets and len(out) >= args.max_sets:
            break
    log.info("MaveDB search scanned %d score sets, kept %d with a UniProt xref",
             n_scanned, len(out))
    return out


def fetch_scores(sess, urn, args, cache_dir):
    cached = cache_dir / f"{urn.replace(':', '_').replace('/', '_')}.csv"
    if cached.exists() and cached.stat().st_size > 0:
        return cached.read_text()
    r = _req(sess, f"{API}/score-sets/{urn}/scores", args.timeout, args.retries,
             headers={"accept": "text/csv"})
    if r is None:
        return None
    cached.write_text(r.text)
    return r.text


def parse_scores(csv_text, urn, gene, uniprot):
    rows = []
    rdr = csv.DictReader(io.StringIO(csv_text))
    for i, rec in enumerate(rdr):
        hp = rec.get("hgvs_pro") or rec.get("hgvs_p") or ""
        score = rec.get("score", "")
        if not hp or score in ("", "NA", "nan", None):
            continue
        is_double = ";" in hp or hp.count("p.") > 1
        m = _POS_RE.search(hp)
        if not m:
            continue
        rows.append({
            "uniprot": uniprot,
            "gene_name": gene,
            "urn": urn,
            "mavedb_id": f"{urn}#{i}",
            "prot_expr": hp,
            "protein_start": m.group(1),
            "score": score,
            "is_double_mutant": str(bool(is_double)),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description="Fetch MaveDB scores (official API)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--uniprot", default="",
                    help="comma-separated accessions to keep (default: all)")
    ap.add_argument("--urn", default="", help="fetch a single score set (test)")
    ap.add_argument("--max_sets", type=int, default=0, help="0 = no limit")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--retries", type=int, default=4)
    ap.add_argument("--cache_dir", default="")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "mavedb_raw.tsv"
    cache_dir = Path(args.cache_dir) if args.cache_dir else (outdir / "raw_scores")
    cache_dir.mkdir(parents=True, exist_ok=True)

    uniprot_filter = {u.strip() for u in args.uniprot.split(",") if u.strip()} or None
    sess = _session()

    sets = enumerate_score_sets(sess, args, uniprot_filter)
    if not sets:
        with out_path.open("w") as fh:        # always emit a valid (empty) file
            fh.write("\t".join(OUT_COLS) + "\n")
        log.warning("no score sets retrieved (API down or no matches) — wrote empty %s", out_path)
        return

    all_rows = []
    for urn, gene, up in sets:
        txt = fetch_scores(sess, urn, args, cache_dir)
        if not txt:
            continue
        n = len(all_rows)
        all_rows.extend(parse_scores(txt, urn, gene, up))
        log.info("  %s (%s/%s): +%d rows", urn, gene, up, len(all_rows) - n)

    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)
    log.info("MaveDB: %d score sets → %d rows → %s",
             len(sets), len(all_rows), out_path)


if __name__ == "__main__":
    main()
