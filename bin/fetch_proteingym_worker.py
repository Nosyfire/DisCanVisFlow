#!/usr/bin/env python3
"""
fetch_proteingym_worker.py — download ProteinGym DMS *substitution* assays from
the OFFICIAL ProteinGym release (https://proteingym.org,
https://github.com/OATML-Markslab/ProteinGym) and emit a UniProt-keyed long
table that the mapping step (create_proteingym_worker.py --mapping_mode uniprot)
turns into Protein_ID rows.

Why this source: ProteinGym is the standard variant-effect-prediction benchmark.
The substitution benchmark ships as a single zip of per-assay CSVs
(DMS_ProteinGym_substitutions.zip), and a reference CSV (DMS_substitutions.csv)
maps each DMS_id → UniProt_ID / target gene. We join the two so every variant
carries its UniProt accession for coordinate mapping onto Gencode isoforms.

Default hosts (override with --zip_url / --ref_csv):
  zip : https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip
  ref : reference_files/DMS_substitutions.csv in the OATML-Markslab/ProteinGym repo
        (raw.githubusercontent.com)

Per-assay CSV columns:  mutant  mutated_sequence  DMS_score  DMS_score_bin
Reference CSV columns :  DMS_id  DMS_filename  UniProt_ID  molecule_name ...
  (UniProt_ID is e.g. "P04637_HUMAN"; we keep the accession part "P04637".)

Output columns (proteingym_raw.tsv):
  uniprot  gene_name  DMS_id  protein_variant  pos  DMS_score  DMS_score_bin

Usage:
  fetch_proteingym_worker.py --out proteingym_raw.tsv
                             [--ref_csv FILE_OR_URL] [--zip_url URL]
                             [--cache_dir DIR] [--uniprot P04637,...]
                             [--limit N] [--max_assays N]
                             [--timeout 300] [--retries 4]
"""

import argparse
import csv
import io
import logging
import re
import sys
import time
import zipfile
from pathlib import Path

import requests

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

DEFAULT_VERSION = "v1.3"
DEFAULT_ZIP_URL = (
    f"https://marks.hms.harvard.edu/proteingym/ProteinGym_{DEFAULT_VERSION}/"
    "DMS_ProteinGym_substitutions.zip"
)
DEFAULT_REF_URL = (
    "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/"
    "reference_files/DMS_substitutions.csv"
)

OUT_COLS = ["uniprot", "gene_name", "DMS_id", "protein_variant", "pos",
            "DMS_score", "DMS_score_bin"]

# leading WT residue (1 letter) + position, e.g. "G145R" → (145, "G").
# Multi-mutant strings ("G145R:A200T") — take the first single substitution.
_MUT_RE = re.compile(r"^([A-Za-z])(\d+)")


# ---------------------------------------------------------------------------
# pure helpers (unit-tested without network)
# ---------------------------------------------------------------------------
def pos_from_mutant(mutant):
    """Return (int position, WT aa) from a ProteinGym mutant string, or (None, '')."""
    if not isinstance(mutant, str) or not mutant:
        return (None, "")
    first = mutant.split(":")[0].strip()
    m = _MUT_RE.match(first)
    if not m:
        return (None, "")
    return (int(m.group(2)), m.group(1).upper())


# A real UniProt accession (6 or 10 char form). ProteinGym's UniProt_ID column is
# usually the entry MNEMONIC (e.g. 'P53_HUMAN', 'BRCA1_HUMAN'), NOT
# '<accession>_HUMAN' — so 'P53_HUMAN'.split('_')[0] = 'P53' is a wrong accession.
_ACC_RE = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")


def uniprot_accession(uniprot_id):
    """Return a usable join key from a ProteinGym UniProt_ID.

    'P04637_HUMAN' → 'P04637' (real accession head);
    'P04637'       → 'P04637';
    'P53_HUMAN'    → 'P53_HUMAN' (mnemonic — kept whole; the mapper resolves it
                    to an accession via the seq table's Entry_Name).
    """
    if not uniprot_id:
        return ""
    s = str(uniprot_id).strip()
    head = s.split("_")[0].strip()
    if _ACC_RE.fullmatch(head):
        return head
    return s


def parse_reference(handle):
    """Parse the reference CSV → {DMS_id: {uniprot, gene_name, filename}}."""
    ref = {}
    rdr = csv.DictReader(handle)
    for rec in rdr:
        dms_id = (rec.get("DMS_id") or "").strip()
        if not dms_id:
            continue
        filename = (rec.get("DMS_filename") or f"{dms_id}.csv").strip()
        ref[dms_id] = {
            "uniprot": uniprot_accession(rec.get("UniProt_ID", "")),
            "gene_name": (rec.get("molecule_name") or rec.get("gene_name") or "").strip(),
            "filename": filename,
        }
    return ref


def parse_assay(handle, dms_id, uniprot, gene_name, limit=0):
    """Parse one per-assay CSV → list of long rows (OUT_COLS).

    Unparseable mutants (e.g. 'wt', multi-mutant with no leading single) are
    skipped. `limit` caps the number of emitted rows (for testing)."""
    rows = []
    rdr = csv.DictReader(handle)
    for rec in rdr:
        mutant = (rec.get("mutant") or "").strip()
        pos, _wt = pos_from_mutant(mutant)
        if pos is None:
            continue
        rows.append({
            "uniprot": uniprot,
            "gene_name": gene_name,
            "DMS_id": dms_id,
            "protein_variant": mutant,
            "pos": str(pos),
            "DMS_score": (rec.get("DMS_score") or "").strip(),
            "DMS_score_bin": (rec.get("DMS_score_bin") or "").strip(),
        })
        if limit and len(rows) >= limit:
            break
    return rows


# ---------------------------------------------------------------------------
# download (network) helpers
# ---------------------------------------------------------------------------
def _session():
    s = requests.Session()
    s.headers.update({"User-Agent": "discanvis-pipeline/1.0"})
    return s


def _req(sess, url, timeout, retries, **kw):
    """GET with exponential backoff (mirrors fetch_mavedb_worker._req).

    The official ProteinGym mirror (marks.hms.harvard.edu) serves an incomplete
    TLS chain (missing intermediate cert), so verified GETs raise SSLError even
    with an up-to-date CA bundle. On the first SSLError we fall back to an
    unverified download with a loud warning — the payload is a public benchmark
    zip whose own CRC covers integrity. Override with --no_insecure_ssl."""
    last = None
    insecure_done = "verify" in kw and kw["verify"] is False
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
        except requests.exceptions.SSLError as exc:
            last = str(exc)
            if not insecure_done and getattr(sess, "allow_insecure", True):
                log.warning("SSL verification failed for %s (%s) — retrying "
                            "WITHOUT certificate verification (public data; "
                            "integrity covered by the zip CRC).",
                            url, type(exc).__name__)
                kw["verify"] = False
                insecure_done = True
                try:
                    import urllib3
                    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
                except Exception:
                    pass
                continue
            wait = min(60, 2 ** attempt)
            log.warning("SSLError on %s (try %d/%d) — backoff %ss",
                        url, attempt, retries, wait)
            time.sleep(wait)
        except requests.RequestException as exc:
            last = str(exc)
            wait = min(60, 2 ** attempt)
            log.warning("%s on %s (try %d/%d) — backoff %ss",
                        exc, url, attempt, retries, wait)
            time.sleep(wait)
    log.error("giving up on %s (%s)", url, last)
    return None


def _load_reference(ref_arg, sess, cache_dir, args):
    """Read the reference CSV from a local path or URL (cached)."""
    if ref_arg and Path(ref_arg).exists():
        with open(ref_arg, newline="", encoding="utf-8") as fh:
            return parse_reference(fh)
    url = ref_arg or DEFAULT_REF_URL
    cached = cache_dir / "DMS_substitutions.csv"
    if cached.exists() and cached.stat().st_size > 0:
        with cached.open(newline="", encoding="utf-8") as fh:
            return parse_reference(fh)
    r = _req(sess, url, args.timeout, args.retries)
    if r is None:
        return None
    cached.write_text(r.text, encoding="utf-8")
    return parse_reference(io.StringIO(r.text))


def _load_zip(zip_arg, sess, cache_dir, args):
    """Return a zipfile.ZipFile for the substitutions archive (cached on disk)."""
    if zip_arg and Path(zip_arg).exists():
        return zipfile.ZipFile(zip_arg)
    url = zip_arg or DEFAULT_ZIP_URL
    cached = cache_dir / "DMS_ProteinGym_substitutions.zip"
    if not (cached.exists() and cached.stat().st_size > 0):
        log.info("downloading %s → %s", url, cached)
        r = _req(sess, url, args.timeout, args.retries, stream=True)
        if r is None:
            return None
        with cached.open("wb") as fh:
            for ch in r.iter_content(chunk_size=1 << 20):
                if ch:
                    fh.write(ch)
    try:
        return zipfile.ZipFile(cached)
    except zipfile.BadZipFile:
        log.error("downloaded archive is not a valid zip: %s", cached)
        return None


def _zip_member(zf, filename, dms_id):
    """Find the per-assay member in the zip by filename or DMS_id."""
    names = zf.namelist()
    targets = {filename, f"{dms_id}.csv"}
    for n in names:
        base = n.rsplit("/", 1)[-1]
        if base in targets:
            return n
    return None


def main():
    ap = argparse.ArgumentParser(description="Fetch ProteinGym DMS substitution scores")
    ap.add_argument("--out", required=True, help="output proteingym_raw.tsv")
    ap.add_argument("--ref_csv", default="",
                    help="reference CSV (local path or URL); default OATML repo")
    ap.add_argument("--zip_url", default="",
                    help="substitutions zip (local path or URL); default marks.hms.harvard.edu")
    ap.add_argument("--cache_dir", default="",
                    help="download cache dir (default: alongside --out)")
    ap.add_argument("--uniprot", default="",
                    help="comma-separated accessions to keep (default: all)")
    ap.add_argument("--limit", type=int, default=0,
                    help="cap rows PER ASSAY (testing); 0 = no limit")
    ap.add_argument("--max_assays", type=int, default=0,
                    help="cap number of assays processed (testing); 0 = no limit")
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--retries", type=int, default=4)
    ap.add_argument("--no_insecure_ssl", action="store_true",
                    help="do NOT fall back to an unverified download when the "
                         "ProteinGym mirror's TLS chain fails (default: fall back)")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir) if args.cache_dir else out_path.parent / "proteingym_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    uniprot_filter = {u.strip() for u in args.uniprot.split(",") if u.strip()} or None
    sess = _session()
    sess.allow_insecure = not args.no_insecure_ssl

    ref = _load_reference(args.ref_csv, sess, cache_dir, args)
    if ref is None:
        log.error("could not load ProteinGym reference CSV — aborting")
        sys.exit(1)
    log.info("reference: %d DMS assays", len(ref))

    zf = _load_zip(args.zip_url, sess, cache_dir, args)
    if zf is None:
        log.error("could not load ProteinGym substitutions zip — aborting")
        sys.exit(1)

    all_rows = []
    n_assays = 0
    for dms_id, meta in ref.items():
        up = meta["uniprot"]
        if not up:
            continue
        if uniprot_filter and up not in uniprot_filter:
            continue
        member = _zip_member(zf, meta["filename"], dms_id)
        if member is None:
            log.warning("assay %s (%s) not found in zip — skipped", dms_id, meta["filename"])
            continue
        try:
            with zf.open(member) as raw:
                txt = io.TextIOWrapper(raw, encoding="utf-8")
                rows = parse_assay(txt, dms_id, up, meta["gene_name"], limit=args.limit)
        except (KeyError, zipfile.BadZipFile, UnicodeDecodeError) as exc:
            log.warning("failed to read %s: %s", member, exc)
            continue
        all_rows.extend(rows)
        n_assays += 1
        if args.max_assays and n_assays >= args.max_assays:
            break

    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)
    log.info("ProteinGym: %d assays → %d variant rows → %s",
             n_assays, len(all_rows), out_path)


if __name__ == "__main__":
    main()
