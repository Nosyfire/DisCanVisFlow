#!/usr/bin/env python3
"""
create_go_worker.py — Module 5f: GO Term Annotation

Maps Gene Ontology terms to GENCODE transcripts (Protein_ID) via UniProt
accession from the GOA human annotation file.

Parent terms are expanded so every ancestor is included (biological process,
molecular function, cellular component hierarchy).

Inputs
------
--loc_chrom    loc_chrom_with_names_isoforms_with_seq.tsv
--goa          goa_human.gaf.gz  (from geneontology.org)
--go_obo       go.obo            (GO ontology definition file)
--output_dir   output directory (default: .)

Outputs
-------
go_terms.tsv  — Protein_ID | Entry_Isoform | GO_Term | name | namespace | def | alt_id | is_a
"""

import argparse
import gzip
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GO OBO parser
# ---------------------------------------------------------------------------

def parse_go_obo(obo_path: str) -> dict:
    """Parse go.obo → {go_id: {'name': str, 'namespace': str, 'def': str, 'alt_id': list, 'is_a': list}}"""
    terms: dict = {}
    current: dict | None = None

    def _open(p):
        return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

    with _open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                current = {"alt_id": [], "is_a": []}
            elif line == "[Typedef]":
                current = None
            elif current is None:
                continue
            elif line.startswith("id: GO:"):
                current["id"] = line[4:]
                terms[current["id"]] = current
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("namespace: "):
                current["namespace"] = line[11:]
            elif line.startswith("def: "):
                # strip quotation marks and citation
                current["def"] = line[5:].split('"')[1] if '"' in line else line[5:]
            elif line.startswith("alt_id: "):
                current["alt_id"].append(line[8:])
            elif line.startswith("is_a: "):
                parent = line[6:].split(" !")[0].strip()
                current["is_a"].append(parent)

    log.info("Parsed GO OBO: %d terms", len(terms))
    return terms


def _all_ancestors(go_id: str, terms: dict, cache: dict | None = None) -> list[str]:
    """Collect go_id + all ancestor IDs (breadth-first)."""
    if cache is None:
        cache = {}
    if go_id in cache:
        return cache[go_id]
    result = [go_id]
    for parent in terms.get(go_id, {}).get("is_a", []):
        for anc in _all_ancestors(parent, terms, cache):
            if anc not in result:
                result.append(anc)
    cache[go_id] = result
    return result


# ---------------------------------------------------------------------------
# GOA parser
# ---------------------------------------------------------------------------

def parse_goa(goa_path: str, acc_set: set) -> dict:
    """
    Parse goa_human.gaf.gz → {uniprot_acc: [go_id, …]}
    Only accessions in acc_set are kept.
    """
    def _open(p):
        return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

    acc_to_goids: dict[str, list] = {}
    with _open(goa_path) as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 5:
                continue
            acc    = parts[1]
            go_id  = parts[4]
            if acc not in acc_set:
                continue
            acc_to_goids.setdefault(acc, [])
            if go_id not in acc_to_goids[acc]:
                acc_to_goids[acc].append(go_id)

    log.info("GOA: %d accessions with GO annotations", len(acc_to_goids))
    return acc_to_goids


# ---------------------------------------------------------------------------
# Build output table
# ---------------------------------------------------------------------------

def build_go_table(
    loc_df:        pd.DataFrame,
    acc_to_goids:  dict,
    terms:         dict,
    expand_parents: bool = True,
) -> pd.DataFrame:
    rows = []

    acc_col = "Entry_Isoform" if "Entry_Isoform" in loc_df.columns else "Accession"
    pid_col = "Protein_ID"

    ancestor_cache: dict = {}

    try:
        from tqdm import tqdm as _tqdm
        _iter = _tqdm(loc_df.iterrows(), total=len(loc_df), desc='GO mapping', unit='isoform')
    except ImportError:
        _iter = loc_df.iterrows()

    for _, row in _iter:
        acc = str(row.get(acc_col, "") or "")
        pid = str(row.get(pid_col, "") or "")
        if not acc or acc in ("nan", "") or not pid or pid in ("nan", ""):
            continue

        go_ids = acc_to_goids.get(acc, [])
        if not go_ids:
            continue

        all_ids = []
        if expand_parents:
            seen_ids: set = set()
            for gid in go_ids:
                for anc in _all_ancestors(gid, terms, ancestor_cache):
                    if anc not in seen_ids:
                        seen_ids.add(anc)
                        all_ids.append(anc)
        else:
            all_ids = go_ids

        for gid in all_ids:
            t = terms.get(gid, {})
            rows.append({
                "Protein_ID":    pid,
                "Entry_Isoform": acc,
                "GO_Term":       gid,
                "name":          t.get("name", ""),
                "namespace":     t.get("namespace", ""),
                "def":           t.get("def", ""),
                "alt_id":        str(t.get("alt_id", [])),
                "is_a":          str(t.get("is_a", [])),
            })

    df = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["Protein_ID", "Entry_Isoform", "GO_Term", "name", "namespace", "def", "alt_id", "is_a"])
    log.info("GO table: %d rows for %d unique proteins",
             len(df), df["Protein_ID"].nunique() if not df.empty else 0)
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Module 5f: GO Term annotation")
    p.add_argument("--loc_chrom",  required=True)
    p.add_argument("--goa",        required=True, help="goa_human.gaf or goa_human.gaf.gz")
    p.add_argument("--go_obo",     required=True, help="go.obo file")
    p.add_argument("--no_parents", action="store_true", default=False,
                   help="Do not expand to parent terms")
    p.add_argument("--output_dir", default=".")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Parsing GO OBO …")
    terms = parse_go_obo(args.go_obo)

    log.info("Loading loc_chrom …")
    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)

    acc_col = "Entry_Isoform" if "Entry_Isoform" in loc_df.columns else "Accession"
    acc_set = set(loc_df[acc_col].dropna().unique())
    log.info("Unique accessions in loc_chrom: %d", len(acc_set))

    log.info("Parsing GOA …")
    acc_to_goids = parse_goa(args.goa, acc_set)

    log.info("Building GO table …")
    df = build_go_table(loc_df, acc_to_goids, terms, expand_parents=not args.no_parents)
    df.to_csv(outdir / "go_terms.tsv", sep="\t", index=False)
    log.info("Done — %d GO annotation rows", len(df))


if __name__ == "__main__":
    main()
