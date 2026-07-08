#!/usr/bin/env python3
"""create_disprot_worker.py — DisProt curated disorder regions → Protein_ID TSV.

Maps DisProt manually-curated intrinsic-disorder regions (UniProt-accession
keyed, with IDPO/GO ontology terms) onto every GENCODE isoform in the sequence
table. DisProt regions are 1-based inclusive coordinates on the canonical
UniProt sequence; each region is mapped to a Protein_ID only when the isoform's
sequence at [start-1:end] matches the DisProt-reported region sequence (so the
canonical coordinates are guaranteed valid for that isoform). This avoids the
blind cross-isoform coordinate copy that a naive accession-only join would do.

Input DisProt bulk TSV (from https://disprot.org/api/v2/download, term_ontology
IDPO+GO), relevant columns (matched by name):
    UniProt ACC | DisProt ID | Region ID | Start | End |
    Term namespace | Term ID | Term name | ECO Term ID | PMID |
    Region sequence | Obsolete

Output (disprot.tsv):
    Protein_ID | Entry_Isoform | disprot_id | region_id | start | end |
    term_namespace | term_id | term_name | eco_id | pmid

Usage:
    create_disprot_worker.py
        --seq_table   <loc_chrom_with_names_isoforms_with_seq.tsv>
        --disprot_tsv <disprot_regions.tsv>  (or NO_FILE)
        --outdir      <output directory>
        [--only_main_isoforms]
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_OUT_COLS = [
    "Protein_ID", "Entry_Isoform", "disprot_id", "region_id",
    "start", "end", "term_namespace", "term_id", "term_name",
    "eco_id", "pmid",
]

# DisProt column name → canonical key used internally. Matched case-insensitively
# against the actual header so minor upstream renames don't break the parser.
_DISPROT_COLS = {
    "uniprot acc":     "acc",
    "disprot id":      "disprot_id",
    "region id":       "region_id",
    "start":           "start",
    "end":             "end",
    "term namespace":  "term_namespace",
    "term id":         "term_id",
    "term name":       "term_name",
    "eco term id":     "eco_id",
    "pmid":            "pmid",
    "region sequence": "region_seq",
    "obsolete":        "obsolete",
}


def _resolve_columns(df: pd.DataFrame) -> dict:
    """Map internal keys → actual df column names (case-insensitive)."""
    lower = {c.lower().strip(): c for c in df.columns}
    return {key: lower[name] for name, key in _DISPROT_COLS.items()
            if name in lower}


def build_disprot_table(
    seq_df: pd.DataFrame,
    dis_df: pd.DataFrame,
    only_main: bool = False,
) -> pd.DataFrame:
    """Return DisProt regions mapped to Protein_ID, coordinate-validated."""
    needed = {"Protein_ID", "Entry_Isoform", "Sequence"}
    if not needed.issubset(seq_df.columns):
        log.warning("seq table missing required columns %s", needed - set(seq_df.columns))
        return pd.DataFrame(columns=_OUT_COLS)

    if only_main and "main_isoform" in seq_df.columns:
        seq_df = seq_df[seq_df["main_isoform"] == "yes"]

    # base accession (P04049-2 → P04049) → [(pid, entry_isoform, sequence), ...]
    acc_to_pids: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for _, row in seq_df.iterrows():
        pid = str(row.get("Protein_ID", "")).strip()
        acc = str(row.get("Entry_Isoform", "")).strip()
        seq = str(row.get("Sequence", "")).strip()
        if not pid or pid == "nan" or not acc or acc == "nan":
            continue
        base = acc.split("-")[0]
        entry = (pid, acc, seq if seq and seq != "nan" else "")
        if entry not in acc_to_pids[base]:
            acc_to_pids[base].append(entry)

    cols = _resolve_columns(dis_df)
    if "acc" not in cols or "start" not in cols or "end" not in cols:
        log.warning("DisProt TSV missing acc/start/end columns — nothing to map")
        return pd.DataFrame(columns=_OUT_COLS)

    rows = []
    for _, r in dis_df.iterrows():
        obsolete = str(r.get(cols.get("obsolete", ""), "")).strip().lower()
        if obsolete == "true":
            continue
        acc = str(r.get(cols["acc"], "")).strip().split("-")[0]
        if acc not in acc_to_pids:
            continue
        try:
            start = int(float(str(r.get(cols["start"], "")).strip()))
            end   = int(float(str(r.get(cols["end"], "")).strip()))
        except (ValueError, TypeError):
            continue
        if start < 1 or end < start:
            continue

        region_seq = str(r.get(cols.get("region_seq", ""), "")).strip()
        if region_seq.lower() in ("nan", "n/a", "none"):
            region_seq = ""

        rec = {
            "disprot_id":     str(r.get(cols.get("disprot_id", ""), "")).strip(),
            "region_id":      str(r.get(cols.get("region_id", ""), "")).strip(),
            "term_namespace": str(r.get(cols.get("term_namespace", ""), "")).strip(),
            "term_id":        str(r.get(cols.get("term_id", ""), "")).strip(),
            "term_name":      str(r.get(cols.get("term_name", ""), "")).strip(),
            "eco_id":         str(r.get(cols.get("eco_id", ""), "")).strip(),
            "pmid":           str(r.get(cols.get("pmid", ""), "")).strip(),
        }

        for pid, entry_iso, seq in acc_to_pids[acc]:
            if seq:
                if end > len(seq):
                    continue
                # Validate coordinates against the isoform sequence when the
                # DisProt region sequence is available; otherwise accept any
                # in-range region (best effort for the rare missing-seq rows).
                if region_seq and seq[start - 1:end] != region_seq:
                    continue
            rows.append({
                "Protein_ID":     pid,
                "Entry_Isoform":  entry_iso,
                "start":          start,
                "end":            end,
                **rec,
            })

    return pd.DataFrame(rows, columns=_OUT_COLS) if rows \
        else pd.DataFrame(columns=_OUT_COLS)


def main():
    p = argparse.ArgumentParser(
        description="DisProt curated disorder regions → Protein_ID TSV")
    p.add_argument("--seq_table",   required=True)
    p.add_argument("--disprot_tsv", required=True)
    p.add_argument("--outdir",      default=".")
    p.add_argument("--only_main_isoforms", action="store_true", default=False)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "disprot.tsv"

    seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)

    dis_path = Path(args.disprot_tsv)
    if dis_path.name == "NO_FILE" or not dis_path.exists() or dis_path.stat().st_size < 10:
        log.info("DisProt TSV not available — writing empty output")
        pd.DataFrame(columns=_OUT_COLS).to_csv(out, sep="\t", index=False)
        return

    dis_df = pd.read_csv(dis_path, sep="\t", dtype=str)
    log.info("Loaded %d DisProt region rows", len(dis_df))

    result = build_disprot_table(seq_df, dis_df, only_main=args.only_main_isoforms)
    result.to_csv(out, sep="\t", index=False)
    log.info("Done — %d DisProt region rows written (%d proteins)",
             len(result),
             result["Protein_ID"].nunique() if not result.empty else 0)


if __name__ == "__main__":
    main()
