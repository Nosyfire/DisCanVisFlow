#!/usr/bin/env python3
"""
create_blast_table_worker.py — Module 0→1 bridge

Parses the two reciprocal BLAST XML outputs and produces the
bestsequences.tsv (and optional allsequences.tsv + isoformssequences.tsv)
that create_id_map_worker.py consumes.

Legacy equivalents
------------------
  create_blast_table.py :: blasttable()   → parse_blast_xml()
  create_blast_table.py :: main()         → main()

Data flow
---------
  blast1_xml  (uniprotdb_gencode_query.xml)
      query = GENCODE protein FASTA headers
      hits  = UniProt entries
      → df1  [Gencode | Uniprot | alignmentpuntcuality | coverage]

  blast2_xml  (gencodedb_uniprot_query.xml)
      query = UniProt FASTA headers
      hits  = GENCODE protein FASTA headers
      → df2  [Uniprot | Gencode | alignmentpuntcuality | coverage]

  pd.merge(df1, df2, on=['Gencode', 'Uniprot'])
      → bestsequences.tsv   (reciprocal hits, coverage ≥ threshold)
      → allsequences.tsv    (all reciprocal hits, no filter)

  blast1_xml (all hits, extended columns)
      → isoformssequences.tsv

Output column contract (bestsequences.tsv / allsequences.tsv):
    Gencode | Uniprot | alignmentpuntcuality_x | coverage_x
    alignmentpuntcuality_y | coverage_y

Output column contract (isoformssequences.tsv):
    Gencode | Uniprot | alignmentpuntcuality_x | coverage_x
    alignmentpuntcuality_y | coverage_y
    identity | identity/aln_len | aln_len
    seq_len_x | seq_len_y | region_len_x | region_len_y
    Starx | Endx | Stary | Endy
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio.Blast import NCBIXML

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
COL_GENCODE   = "Gencode"
COL_UNIPROT   = "Uniprot"
COL_QUALITY   = "alignmentpuntcuality"
COL_COVERAGE  = "coverage"

# bestsequences.tsv / allsequences.tsv column order
BASE_COLS = [COL_GENCODE, COL_UNIPROT, COL_QUALITY, COL_COVERAGE]

# isoformssequences.tsv extended columns (mirroring legacy isoform_sequences=True)
ISO_EXTRA_COLS = [
    "identity", "identity/aln_len", "aln_len",
    "seq_len_x", "seq_len_y",
    "region_len_x", "region_len_y",
    "Starx", "Endx", "Stary", "Endy",
]


# ---------------------------------------------------------------------------
# XML → DataFrame
# ---------------------------------------------------------------------------

def _alignment_quality(
    identities: int,
    align_length: int,
    query_letters: int,
    query_start: int,
    sbjct_start: int,
) -> str:
    """
    Reproduce the legacy three-way quality classification:
      'identical' : perfect full-length hit (same start, all residues match)
      'aligned'   : partial or imperfect hit above coverage threshold
    """
    if (
        align_length == identities
        and query_letters == align_length
        and query_start == sbjct_start
    ):
        return "identical"
    return "aligned"


def parse_blast_xml(
    xml_path: str,
    query_col: str,
    hit_col: str,
    coverage_threshold: float = 0.0,
    include_iso_fields: bool = False,
) -> pd.DataFrame:
    """
    Parse a BLAST XML file (outfmt 5) and return a DataFrame.

    Parameters
    ----------
    xml_path : str
        Path to the BLAST XML output.
    query_col : str
        Column name to use for the query sequence title (e.g. 'Gencode').
    hit_col : str
        Column name to use for the hit sequence title   (e.g. 'Uniprot').
    coverage_threshold : float
        Minimum coverage (%) required to include a hit row. 0 = keep all.
    include_iso_fields : bool
        When True, add the extended HSP fields needed for isoformssequences.tsv.

    Notes
    -----
    In the BLAST XML produced without -parse_seqids, the UniProt header is
    stored in alignment.hit_def (e.g.
    'sp|Q8NH21|OR4F5_HUMAN Olfactory receptor 4F5 OS=Homo sapiens ...')
    while the GENCODE header is stored verbatim in blast_record.query.
    """
    records = []

    with open(xml_path) as fh:
        for blast_record in NCBIXML.parse(fh):
            query_title  = blast_record.query        # full FASTA header
            query_letters = blast_record.query_letters

            for alignment in blast_record.alignments:
                hit_title = alignment.hit_def        # full FASTA header of hit
                hit_len   = alignment.length

                for hsp in alignment.hsps:
                    coverage = round(
                        (hsp.identities / query_letters) * 100, 3
                    )
                    if coverage < coverage_threshold:
                        continue

                    quality = _alignment_quality(
                        identities=hsp.identities,
                        align_length=hsp.align_length,
                        query_letters=query_letters,
                        query_start=hsp.query_start,
                        sbjct_start=hsp.sbjct_start,
                    )

                    row: dict = {
                        query_col: query_title,
                        hit_col:   hit_title,
                        COL_QUALITY:  quality,
                        COL_COVERAGE: coverage,
                    }

                    if include_iso_fields:
                        row["identity"]        = hsp.identities
                        row["identity/aln_len"] = round(
                            hsp.identities / hsp.align_length * 100, 3
                        )
                        row["aln_len"]         = hsp.align_length
                        row["seq_len_x"]       = query_letters
                        row["seq_len_y"]       = hit_len
                        row["Starx"]           = hsp.sbjct_start
                        row["Endx"]            = hsp.sbjct_end
                        row["region_len_x"]    = hsp.sbjct_end - hsp.sbjct_start + 1
                        row["Stary"]           = hsp.query_start
                        row["Endy"]            = hsp.query_end
                        row["region_len_y"]    = hsp.query_end - hsp.query_start + 1

                    records.append(row)

    cols = BASE_COLS + (ISO_EXTRA_COLS if include_iso_fields else [])
    if not records:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(records, columns=cols)


# ---------------------------------------------------------------------------
# Reciprocal merge
# ---------------------------------------------------------------------------

def build_reciprocal_table(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
) -> pd.DataFrame:
    """
    Inner-join df1 (Gencode-query → Uniprot-hit) with
                df2 (Uniprot-query → Gencode-hit)
    on [Gencode, Uniprot] to keep only reciprocal best hits.

    Result columns:
        Gencode | Uniprot
        alignmentpuntcuality_x | coverage_x   (from df1, GENCODE queried vs UniProt)
        alignmentpuntcuality_y | coverage_y   (from df2, UniProt queried vs GENCODE)
    """
    merged = pd.merge(
        df1[[COL_GENCODE, COL_UNIPROT, COL_QUALITY, COL_COVERAGE]],
        df2[[COL_GENCODE, COL_UNIPROT, COL_QUALITY, COL_COVERAGE]],
        on=[COL_GENCODE, COL_UNIPROT],
        suffixes=("_x", "_y"),
    )
    return merged


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Parse reciprocal BLAST XML outputs and produce bestsequences.tsv, "
            "allsequences.tsv, and isoformssequences.tsv"
        )
    )
    p.add_argument(
        "--blast1_xml",
        required=True,
        help="uniprotdb_gencode_query.xml  (query=GENCODE, db=UniProt)",
    )
    p.add_argument(
        "--blast2_xml",
        required=True,
        help="gencodedb_uniprot_query.xml  (query=UniProt, db=GENCODE)",
    )
    p.add_argument(
        "--output_dir",
        default=".",
        help="Directory for output TSV files (default: current dir)",
    )
    p.add_argument(
        "--coverage",
        type=float,
        default=90.0,
        help="Minimum coverage %% for bestsequences.tsv (default: 90)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── blast1: query=GENCODE, db=UniProt ─────────────────────────────────
    log.info("Parsing blast1 (GENCODE → UniProt): %s", args.blast1_xml)
    df1_filtered = parse_blast_xml(
        args.blast1_xml,
        query_col=COL_GENCODE,
        hit_col=COL_UNIPROT,
        coverage_threshold=args.coverage,
    )
    df1_all = parse_blast_xml(
        args.blast1_xml,
        query_col=COL_GENCODE,
        hit_col=COL_UNIPROT,
        coverage_threshold=0.0,
    )
    df1_iso = parse_blast_xml(
        args.blast1_xml,
        query_col=COL_GENCODE,
        hit_col=COL_UNIPROT,
        coverage_threshold=0.0,
        include_iso_fields=True,
    )
    log.info("  blast1 hits (≥%.0f%% cov): %d  |  all: %d",
             args.coverage, len(df1_filtered), len(df1_all))

    # ── blast2: query=UniProt, db=GENCODE ─────────────────────────────────
    log.info("Parsing blast2 (UniProt → GENCODE): %s", args.blast2_xml)
    df2_filtered = parse_blast_xml(
        args.blast2_xml,
        query_col=COL_UNIPROT,
        hit_col=COL_GENCODE,
        coverage_threshold=args.coverage,
    )
    df2_all = parse_blast_xml(
        args.blast2_xml,
        query_col=COL_UNIPROT,
        hit_col=COL_GENCODE,
        coverage_threshold=0.0,
    )
    log.info("  blast2 hits (≥%.0f%% cov): %d  |  all: %d",
             args.coverage, len(df2_filtered), len(df2_all))

    # ── Reciprocal merge ──────────────────────────────────────────────────
    best = build_reciprocal_table(df1_filtered, df2_filtered)
    log.info("Reciprocal hits (bestsequences): %d", len(best))

    all_hits = build_reciprocal_table(df1_all, df2_all)
    log.info("All reciprocal hits:             %d", len(all_hits))

    # ── Write outputs ──────────────────────────────────────────────────────
    best_out = outdir / "bestsequences.tsv"
    best.to_csv(best_out, sep="\t", index=False, header=True)
    log.info("Written: %s", best_out)

    all_out = outdir / "allsequences.tsv"
    all_hits.to_csv(all_out, sep="\t", index=False, header=True)
    log.info("Written: %s", all_out)

    iso_out = outdir / "isoformssequences.tsv"
    df1_iso.to_csv(iso_out, sep="\t", index=False, header=True)
    log.info("Written: %s", iso_out)

    log.info("create_blast_table_worker.py complete.")


if __name__ == "__main__":
    main()
