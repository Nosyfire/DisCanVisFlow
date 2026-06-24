#!/usr/bin/env python3
"""
create_pdb_worker.py — Module 5c: PDB Structure Mapping (structure → transcript).

For every protein this maps experimental PDB structures onto each Gencode
transcript (Protein_ID) of the gene and records, for the DisCanVis2 database:

  1. WHICH region of the transcript is resolved in WHICH structure and chain.
  2. WHICH residues inside a structure's mapped range are *missing* (unobserved
     in the coordinates) — for an X-ray model these unobserved residues are the
     genuine disordered segments.

Data sources (PDBe Graph / REST API, SIFTS-based):
  • mappings/{acc}            UniProt↔structure residue mapping per PDB chain
  • pdb summary/{pdb_id}      resolution + experimental method (cached)
  • polymer_coverage/{pdb_id} observed residue ranges per chain (cached)

UniProt residue numbers are translated to each transcript's coordinates by
locating the covered region in the transcript sequence (exact, else ≥90 %
identity window — same convention as the transcript annotation mapping).

Outputs
-------
  pdb_structures.tsv  Protein_ID | Accession | pdb_id | chain_id | struct_asym_id |
                      entity_id | prot_start | prot_end | unp_start | unp_end |
                      resolution | experimental_method
  pdb_missing.tsv     Protein_ID | Accession | pdb_id | chain_id |
                      prot_start | prot_end | unp_start | unp_end | length
                      (one row per contiguous run of missing/unobserved residues)
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd

try:
    import requests
    _HAS_REQUESTS = True
except ImportError:
    _HAS_REQUESTS = False

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

MAPPINGS_URL = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/{acc}"
SUMMARY_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pid}"
EXPERIMENT_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/{pid}"
COVERAGE_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/{pid}"

STRUCT_COLS = ["Protein_ID", "Accession", "pdb_id", "chain_id", "struct_asym_id",
               "entity_id", "prot_start", "prot_end", "unp_start", "unp_end",
               "resolution", "experimental_method"]
MISSING_COLS = ["Protein_ID", "Accession", "pdb_id", "chain_id",
                "prot_start", "prot_end", "unp_start", "unp_end", "length"]


# ---------------------------------------------------------------------------
# HTTP + simple per-process caches
# ---------------------------------------------------------------------------
_SUMMARY_CACHE: dict = {}
_COVERAGE_CACHE: dict = {}


def _get_json(url: str, delay: float = 0.3, timeout: int = 30):
    if not _HAS_REQUESTS:
        return None
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
        except Exception as exc:  # noqa: BLE001
            log.debug("PDBe request failed (%s): %s", url, exc)
        time.sleep(delay * (attempt + 1))
    return None


def fetch_mappings(acc: str, delay: float) -> dict:
    """pdb_id -> list of chain segments {chain_id, struct_asym_id, entity_id,
    unp_start, unp_end, label_start, label_end}."""
    data = _get_json(MAPPINGS_URL.format(acc=acc), delay)
    out: dict[str, list] = {}
    if not data:
        return out
    pdbs = data.get(acc, {}).get("PDB", {})
    for pdb_id, segs in pdbs.items():
        rows = []
        for s in segs:
            try:
                rows.append({
                    "chain_id": s.get("chain_id", ""),
                    "struct_asym_id": s.get("struct_asym_id", ""),
                    "entity_id": s.get("entity_id", ""),
                    "unp_start": int(s["unp_start"]),
                    "unp_end": int(s["unp_end"]),
                    "label_start": int(s["start"]["residue_number"]),
                    "label_end": int(s["end"]["residue_number"]),
                })
            except (KeyError, ValueError, TypeError):
                continue
        if rows:
            out[pdb_id] = rows
    return out


def fetch_summary(pdb_id: str, delay: float):
    """(resolution, experimental_method) for a PDB entry, both cached."""
    if pdb_id in _SUMMARY_CACHE:
        return _SUMMARY_CACHE[pdb_id]
    method = ""
    data = _get_json(SUMMARY_URL.format(pid=pdb_id), delay)
    if data and data.get(pdb_id):
        methods = data[pdb_id][0].get("experimental_method", []) or []
        method = ", ".join(methods) if isinstance(methods, list) else str(methods)
    resolution = ""
    exp = _get_json(EXPERIMENT_URL.format(pid=pdb_id), delay)
    if exp and exp.get(pdb_id):
        res = exp[pdb_id][0].get("resolution")
        resolution = "" if res is None else str(res)
    _SUMMARY_CACHE[pdb_id] = (resolution, method)
    return resolution, method


def fetch_coverage(pdb_id: str, delay: float) -> dict:
    """struct_asym_id -> list of observed (label_start, label_end)."""
    if pdb_id in _COVERAGE_CACHE:
        return _COVERAGE_CACHE[pdb_id]
    data = _get_json(COVERAGE_URL.format(pid=pdb_id), delay)
    cov: dict[str, list] = {}
    if data:
        for mol in data.get(pdb_id, {}).get("molecules", []):
            for ch in mol.get("chains", []):
                asym = ch.get("struct_asym_id", "")
                ranges = []
                for obs in ch.get("observed", []):
                    try:
                        ranges.append((int(obs["start"]["residue_number"]),
                                       int(obs["end"]["residue_number"])))
                    except (KeyError, ValueError, TypeError):
                        continue
                cov.setdefault(asym, []).extend(ranges)
    _COVERAGE_CACHE[pdb_id] = cov
    return cov


# ---------------------------------------------------------------------------
# Missing-residue computation (unobserved = disorder)
# ---------------------------------------------------------------------------

def missing_unp_ranges(seg: dict, observed_label: list) -> list:
    """Given one chain segment and its observed label ranges, return the
    contiguous UniProt ranges that are mapped but NOT observed."""
    unp_s, unp_e = seg["unp_start"], seg["unp_end"]
    offset = seg["unp_start"] - seg["label_start"]   # unp = label + offset
    observed_unp = set()
    for ls, le in observed_label:
        a = max(unp_s, ls + offset)
        b = min(unp_e, le + offset)
        for u in range(a, b + 1):
            observed_unp.add(u)
    missing = [u for u in range(unp_s, unp_e + 1) if u not in observed_unp]
    # collapse into contiguous ranges
    ranges = []
    run_start = None
    prev = None
    for u in missing:
        if run_start is None:
            run_start = prev = u
        elif u == prev + 1:
            prev = u
        else:
            ranges.append((run_start, prev))
            run_start = prev = u
    if run_start is not None:
        ranges.append((run_start, prev))
    return ranges


# ---------------------------------------------------------------------------
# UniProt → transcript region mapping
# ---------------------------------------------------------------------------

def best_window(region: str, tgt_seq: str, min_identity: float):
    """Return (start_0based, identity) of best ungapped placement, or (-1, id)."""
    L = len(region)
    if L == 0 or len(tgt_seq) < L:
        return -1, 0.0
    idx = tgt_seq.find(region)
    if idx != -1:
        return idx, 1.0
    best_idx, best_id = -1, 0.0
    for i in range(0, len(tgt_seq) - L + 1):
        m = sum(1 for a, b in zip(region, tgt_seq[i:i + L]) if a == b)
        ident = m / L
        if ident > best_id:
            best_id, best_idx = ident, i
    return (best_idx, best_id) if best_id >= min_identity else (-1, best_id)


def build_gene_lookup(loc_df: pd.DataFrame):
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"]
                     if c in loc_df.columns), None)
    gene_to_rows: dict[str, list] = {}
    acc_to_seq: dict[str, str] = {}
    acc_main: dict[str, bool] = {}
    for _, row in loc_df.iterrows():
        acc = str(row.get("Entry_Isoform", ""))
        pid = str(row.get("Protein_ID", ""))
        seq = str(row.get("Sequence", ""))
        gene = str(row.get(gene_col, "")) if gene_col else ""
        if not acc or acc == "nan":
            continue
        seq = "" if seq == "nan" else seq
        if seq:
            acc_to_seq[acc] = seq
        is_main = str(row.get("main_isoform", "")).lower() == "yes"
        acc_main[acc] = acc_main.get(acc, False) or is_main
        gene_to_rows.setdefault(gene, [])
        entry = (acc, pid, seq, gene)
        if entry not in gene_to_rows[gene]:
            gene_to_rows[gene].append(entry)
    return gene_to_rows, acc_to_seq, acc_main, gene_col


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Module 5c: PDB structure → transcript mapping")
    p.add_argument("--loc_chrom", required=True)
    p.add_argument("--output_dir", default=".")
    p.add_argument("--request_delay", type=float, default=0.2)
    p.add_argument("--min_identity", type=float, default=0.9)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)
    if "Entry_Isoform" not in loc_df.columns:
        log.error("loc_chrom missing 'Entry_Isoform' column")
        pd.DataFrame(columns=STRUCT_COLS).to_csv(outdir / "pdb_structures.tsv", sep="\t", index=False)
        pd.DataFrame(columns=MISSING_COLS).to_csv(outdir / "pdb_missing.tsv", sep="\t", index=False)
        return

    gene_to_rows, acc_to_seq, acc_main, gene_col = build_gene_lookup(loc_df)

    # Fetch structures for the main isoform of each gene (canonical UniProt acc)
    fetch_accs = [a for a, is_main in acc_main.items() if is_main] or list(acc_to_seq)
    log.info("Fetching PDB structures for %d accessions …", len(fetch_accs))

    try:
        from tqdm import tqdm as _tqdm
        acc_iter = _tqdm(fetch_accs, desc="PDB fetch", unit="protein")
    except ImportError:
        acc_iter = fetch_accs

    struct_rows, missing_rows = [], []

    for acc in acc_iter:
        ref_seq = acc_to_seq.get(acc, "")
        if not ref_seq:
            continue
        gene = next((g for (a, p, s, g) in
                     [r for rows in gene_to_rows.values() for r in rows] if a == acc), "")
        siblings = gene_to_rows.get(gene, [(acc, "", ref_seq, gene)])

        mappings = fetch_mappings(acc, args.request_delay)
        time.sleep(args.request_delay)

        for pdb_id, segs in mappings.items():
            resolution, method = fetch_summary(pdb_id, args.request_delay)
            coverage = fetch_coverage(pdb_id, args.request_delay)

            for seg in segs:
                unp_s, unp_e = seg["unp_start"], seg["unp_end"]
                if unp_s < 1 or unp_e < unp_s or unp_e > len(ref_seq):
                    # mapping outside this isoform's sequence — skip
                    continue
                region = ref_seq[unp_s - 1: unp_e]
                miss = missing_unp_ranges(seg, coverage.get(seg["struct_asym_id"], []))

                for (sib_acc, pid, sib_seq, _g) in siblings:
                    if not sib_seq or not pid:
                        continue
                    idx, ident = best_window(region, sib_seq, args.min_identity)
                    if idx == -1:
                        continue
                    prot_off = idx + 1 - unp_s        # prot = unp + prot_off
                    struct_rows.append({
                        "Protein_ID": pid, "Accession": acc, "pdb_id": pdb_id,
                        "chain_id": seg["chain_id"], "struct_asym_id": seg["struct_asym_id"],
                        "entity_id": seg["entity_id"],
                        "prot_start": unp_s + prot_off, "prot_end": unp_e + prot_off,
                        "unp_start": unp_s, "unp_end": unp_e,
                        "resolution": resolution, "experimental_method": method,
                    })
                    for (ms, me) in miss:
                        missing_rows.append({
                            "Protein_ID": pid, "Accession": acc, "pdb_id": pdb_id,
                            "chain_id": seg["chain_id"],
                            "prot_start": ms + prot_off, "prot_end": me + prot_off,
                            "unp_start": ms, "unp_end": me, "length": me - ms + 1,
                        })

    pd.DataFrame(struct_rows, columns=STRUCT_COLS).to_csv(
        outdir / "pdb_structures.tsv", sep="\t", index=False)
    pd.DataFrame(missing_rows, columns=MISSING_COLS).to_csv(
        outdir / "pdb_missing.tsv", sep="\t", index=False)
    log.info("PDB: %d structure-chain rows, %d missing-residue (disorder) rows",
             len(struct_rows), len(missing_rows))


if __name__ == "__main__":
    main()
