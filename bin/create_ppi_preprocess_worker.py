#!/usr/bin/env python3
"""
Module 5j-prep — PPI raw-file preprocessing.

Converts raw IntAct MiTab2.7, BioGRID MiTab, and HIPPIE interaction files
into the standard Interaction_*.tsv format consumed by create_ppi_worker.py.

Filtering:
  - Human-only interactions (taxon:9606 for both interactors)
  - Self-interactions excluded (gene A == gene B)

Output columns match the pre-processed format expected by create_ppi_worker.py:
  Accession A | Accession B | ID Interactor A | ID Interactor B
  | Interaction Detection Methods | Publication Identifiers | Confidence Value

"Accession A/B" are UniProt base accessions (P04049, not P04049-2).
create_ppi_worker.py is extended to resolve these via loc_chrom.

Usage:
  create_ppi_preprocess_worker.py
      --intact   <intact_human.mitab or NO_FILE>
      --biogrid  <biogrid_human.mitab or NO_FILE>
      --hippie   <hippie_current.txt  or NO_FILE>
      --outdir   <output directory>
"""

import argparse
import gzip
import logging
import re
import sys
import zipfile
from io import TextIOWrapper
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_OUT_COLS = [
    "Accession A", "Accession B",
    "ID Interactor A", "ID Interactor B",
    "Interaction Detection Methods", "Publication Identifiers", "Confidence Value",
]


def _open_any(path: Path, prefer: str | None = None):
    """Return a text stream for .gz, .zip, or plain file.

    For a multi-member zip, open the first member whose name contains ``prefer``
    (case-insensitive) when given — e.g. the all-organism BioGRID archive bundles
    one MiTab per species, so we must pick the Homo_sapiens member rather than the
    alphabetically-first one. Falls back to the first member otherwise.
    """
    if not path.exists() or path.stat().st_size == 0:
        return None
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if path.suffix == ".zip":
        zf = zipfile.ZipFile(path)
        names = [n for n in zf.namelist() if not n.endswith("/")]
        if not names:
            return None
        chosen = names[0]
        if prefer:
            for n in names:
                if prefer.lower() in n.lower():
                    chosen = n
                    break
        return TextIOWrapper(zf.open(chosen), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _parse_uniprot(field: str) -> str:
    """Extract first UniProt base accession from e.g. 'uniprotkb:P04049-2'."""
    for part in str(field).split("|"):
        part = part.strip()
        if part.startswith("uniprotkb:"):
            acc = part[len("uniprotkb:"):]
            return acc.split("-")[0]
    return ""


def _is_human(taxon_field: str) -> bool:
    return "taxid:9606" in str(taxon_field)


def _pubmed_from_intact(pub_field: str) -> str:
    """Return pipe-joined pubmed:NNN entries from an IntAct publication field."""
    return "|".join(p for p in str(pub_field).split("|") if "pubmed:" in p.lower())


def parse_intact(path: Path) -> pd.DataFrame:
    """Parse IntAct MiTab2.7 file (plain, .gz, or .zip)."""
    fh = _open_any(path)
    if fh is None:
        return pd.DataFrame(columns=_OUT_COLS)

    rows = []
    try:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            if not (_is_human(parts[9]) and _is_human(parts[10])):
                continue
            acc_a = _parse_uniprot(parts[0])
            acc_b = _parse_uniprot(parts[1])
            if not acc_a or not acc_b or acc_a == acc_b:
                continue
            rows.append({
                "Accession A": acc_a,
                "Accession B": acc_b,
                "ID Interactor A": acc_a,
                "ID Interactor B": acc_b,
                "Interaction Detection Methods": parts[6].strip(),
                "Publication Identifiers": _pubmed_from_intact(parts[8]),
                "Confidence Value": _parse_intact_confidence(parts[14]),
            })
    finally:
        fh.close()

    df = pd.DataFrame(rows, columns=_OUT_COLS)
    log.info("IntAct: %d human interactions parsed", len(df))
    return df


def _parse_intact_confidence(conf_field: str) -> str:
    for part in str(conf_field).split("|"):
        m = re.search(r"miscore:([\d.]+)", part, re.IGNORECASE)
        if m:
            return m.group(1)
    return "0"


def parse_biogrid(path: Path) -> pd.DataFrame:
    """Parse BioGRID MiTab file (plain, .gz, or .zip).

    BioGRID MiTab uses Entrez Gene IDs in cols 1-2, Swiss-Prot in cols 24-25,
    gene symbols in cols 7-8, organism in cols 15-16.

    The download is the all-organism archive (one MiTab per species), so select
    the Homo_sapiens member; row-level taxon filtering below is the safety net.
    """
    fh = _open_any(path, prefer="Homo_sapiens")
    if fh is None:
        return pd.DataFrame(columns=_OUT_COLS)

    rows = []
    try:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 26:
                continue
            org_a = str(parts[15]).strip()
            org_b = str(parts[16]).strip()
            if org_a != "9606" or org_b != "9606":
                continue
            acc_a = _first_token(parts[24])
            acc_b = _first_token(parts[25])
            if not acc_a or not acc_b or acc_a == acc_b:
                continue
            pubmed = f"pubmed:{parts[14].strip()}" if parts[14].strip().isdigit() else parts[14].strip()
            rows.append({
                "Accession A": acc_a,
                "Accession B": acc_b,
                "ID Interactor A": acc_a,
                "ID Interactor B": acc_b,
                "Interaction Detection Methods": parts[11].strip() if len(parts) > 11 else "",
                "Publication Identifiers": pubmed,
                "Confidence Value": parts[18].strip() if len(parts) > 18 else "0",
            })
    finally:
        fh.close()

    df = pd.DataFrame(rows, columns=_OUT_COLS)
    log.info("BioGRID: %d human interactions parsed", len(df))
    return df


def _first_token(field: str) -> str:
    """First non-empty token from a pipe/bar-separated field."""
    for t in str(field).split("|"):
        t = t.strip().split("-")[0]
        if t and t != "-":
            return t
    return ""


def parse_hippie(path: Path) -> pd.DataFrame:
    """Parse HIPPIE custom format: GeneA UniprotA GeneB UniprotB score evidence."""
    fh = _open_any(path)
    if fh is None:
        return pd.DataFrame(columns=_OUT_COLS)

    rows = []
    try:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            acc_a = parts[1].strip().split("-")[0]
            acc_b = parts[3].strip().split("-")[0]
            if not acc_a or not acc_b or acc_a == acc_b:
                continue
            evidence = parts[5].strip()
            pmids = "|".join(
                f"pubmed:{m}"
                for m in re.findall(r"pmids?:(\d+)", evidence, re.IGNORECASE)
            )
            rows.append({
                "Accession A": acc_a,
                "Accession B": acc_b,
                "ID Interactor A": acc_a,
                "ID Interactor B": acc_b,
                "Interaction Detection Methods": "",
                "Publication Identifiers": pmids,
                "Confidence Value": parts[4].strip(),
            })
    finally:
        fh.close()

    df = pd.DataFrame(rows, columns=_OUT_COLS)
    log.info("HIPPIE: %d interactions parsed", len(df))
    return df


def main():
    p = argparse.ArgumentParser(description="Module 5j-prep — PPI raw preprocessing")
    p.add_argument("--intact",  required=True)
    p.add_argument("--biogrid", required=True)
    p.add_argument("--hippie",  required=True)
    p.add_argument("--outdir",  required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    no_file = Path("NO_FILE")

    for src_arg, parse_fn, out_name in [
        (args.intact,  parse_intact,  "Interaction_intact.tsv"),
        (args.biogrid, parse_biogrid, "Interaction_biogrid.tsv"),
        (args.hippie,  parse_hippie,  "Interaction_hippie.tsv"),
    ]:
        src = Path(src_arg)
        if src.name == "NO_FILE" or not src.exists() or src.stat().st_size == 0:
            pd.DataFrame(columns=_OUT_COLS).to_csv(outdir / out_name, sep="\t", index=False)
            continue
        df = parse_fn(src)
        df.to_csv(outdir / out_name, sep="\t", index=False)
        log.info("Written %d rows → %s", len(df), out_name)


if __name__ == "__main__":
    main()
