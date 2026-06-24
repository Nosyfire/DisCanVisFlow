#!/usr/bin/env python3
"""
create_exon_worker.py — Exon Boundary Extraction

Parses combined_map.map (Module 3 output) and identifies exon boundaries
for each transcript by looking for >20 bp gaps between consecutive mapped
residues' genomic coordinates.

Algorithm (ported from legacy exon.py):
  For each protein block in combined_map.map:
    - Walk residues in order.
    - Skip residues with no genomic coordinate ('-').
    - When the absolute difference between consecutive genomic positions > 20 bp
      → new exon boundary.
  Output: "exon_index  aa_start  aa_end" lines per protein.

Inputs
------
  --combined_map  combined_map.map (Module 3 output)
  --loc_chrom     loc_chrom_with_names_isoforms_with_seq.tsv (Module 2 output)
  --output_dir    output directory (default: .)

Outputs
-------
  exon.tsv   — long format, one row per exon:
             Protein_ID | exon_number | total_exons | aa_start | aa_end |
             aa_length | genomic_start | genomic_end
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# combined_map.map parser
# ---------------------------------------------------------------------------

def parse_combined_map(map_path: str) -> dict:
    """
    Parse combined_map.map into:
        {transcript_id: {'chromosome': str, 'strand': str,
                         'residues': [{'num': int, 'aa': str, 'gene': [str,...]}]}}

    Header lines: # FASTA_header CHR STRAND START-END
    Residue lines: num aa cdna_0,1,2 codon aa chrom_0,1,2, chrom_codon chrom_aa
    """
    blocks: dict = {}
    current_key: str | None = None
    current_block: dict = {}

    with open(map_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.lstrip("#").split()
                if len(parts) < 3:
                    continue
                fasta_header = parts[0]
                chrom        = parts[1]
                strand       = parts[2]
                # GENCODE pipe-delimited header, e.g.
                #   ENST00000691899.1|ENSG...|OTT...|-|RAF1-262|RAF1|...
                # field 0 = transcript stable ID, field 4 = Protein_ID
                # (Gencode transcript name). loc_chrom keys on Protein_ID,
                # so register the block under every pipe field for lookup.
                header_fields = fasta_header.split("|")
                current_block = {
                    "chromosome": chrom,
                    "strand":     strand,
                    "residues":   [],
                }
                current_key = header_fields[0]
                for field in header_fields:
                    if field and field != "-":
                        blocks[field] = current_block
            else:
                if current_key is None:
                    continue
                cols = line.split()
                # combined_map.map residue lines have 8 whitespace-separated
                # columns; we only need up to col[5] (genomic positions).
                if len(cols) < 6:
                    continue
                try:
                    num    = int(cols[0])
                    aa     = cols[1]
                    # Genomic positions in col 5 → "pos0,pos1,pos2,"
                    gene   = cols[5].rstrip(",").split(",")
                except (ValueError, IndexError):
                    continue
                current_block["residues"].append({
                    "num":  num,
                    "aa":   aa,
                    "gene": gene,
                })

    log.info("Parsed %d transcript blocks from combined_map.map", len(blocks))
    return blocks


# ---------------------------------------------------------------------------
# Exon boundary calculation
# ---------------------------------------------------------------------------
EXON_GAP_BP = 20   # genomic gap threshold (bp) between consecutive residues


def calc_exons(ts_id: str, residues: list[dict]) -> list[dict]:
    """
    Split a transcript's CDS residues into exons at genomic gaps > EXON_GAP_BP.

    Returns a list of exon dicts (in transcript order), each with:
        aa_start, aa_end   1-based protein residue range covered by the exon
        gstart,  gend      genomic coordinate of the first/last CDS nucleotide
    Boundaries derived from combined_map.map per-residue genomic coordinates.
    """
    mapped = [r for r in residues if r["gene"] and r["gene"][0] != "-"]
    if not mapped:
        return []

    def _g(res, last=False):
        gl = res["gene"]
        return gl[-1] if last else gl[0]

    exons = []
    cur = {"beg": mapped[0]["num"], "end": mapped[0]["num"],
           "gstart": _g(mapped[0]), "gend": _g(mapped[0], last=True)}
    try:
        prev_g = int(_g(mapped[0]))
    except (ValueError, TypeError):
        prev_g = 0

    for res in mapped[1:]:
        try:
            g0 = int(_g(res))
        except (ValueError, TypeError):
            g0 = prev_g
        if abs(g0 - prev_g) > EXON_GAP_BP:
            exons.append(cur)
            cur = {"beg": res["num"], "end": res["num"],
                   "gstart": _g(res), "gend": _g(res, last=True)}
        else:
            cur["end"] = res["num"]
            cur["gend"] = _g(res, last=True)
        prev_g = g0

    exons.append(cur)
    return exons


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract exon boundaries from combined_map.map")
    p.add_argument("--combined_map", required=True,
                   help="combined_map.map (Module 3 output)")
    p.add_argument("--loc_chrom",    required=True,
                   help="loc_chrom_with_names_isoforms_with_seq.tsv (Module 2 output)")
    p.add_argument("--output_dir",   default=".",
                   help="Output directory (default: current dir)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)
    # Accept both 'Protein_ID' and 'transcript_stable_id' style id columns
    id_col = "Protein_ID" if "Protein_ID" in loc_df.columns else "transcript_stable_id"

    blocks = parse_combined_map(args.combined_map)

    rows = []
    skipped = 0
    for pid in tqdm(loc_df[id_col].dropna().unique(), desc="Exon calc"):
        # Try to find block by Protein_ID (e.g. 'RAF1-201') OR transcript id
        block = blocks.get(pid)
        if block is None:
            # Try matching via transcript_stable_id column
            matching = loc_df[loc_df[id_col] == pid]
            if not matching.empty and "transcript_stable_id" in matching.columns:
                ts_raw = matching.iloc[0]["transcript_stable_id"]
                ts_base = ts_raw.split(".")[0] if isinstance(ts_raw, str) else ts_raw
                block = blocks.get(ts_base) or blocks.get(ts_raw)
        if block is None:
            skipped += 1
            continue
        exons = calc_exons(pid, block["residues"])
        if not exons:
            continue
        total = len(exons)
        for i, ex in enumerate(exons, start=1):
            aa_start = ex["beg"] + 1          # map num is 0-based → 1-based
            aa_end   = ex["end"] + 1
            rows.append({
                "Protein_ID":   pid,
                "exon_number":  i,
                "total_exons":  total,
                "aa_start":     aa_start,
                "aa_end":       aa_end,
                "aa_length":    aa_end - aa_start + 1,
                "genomic_start": ex["gstart"],
                "genomic_end":   ex["gend"],
            })

    n_proteins = len({r["Protein_ID"] for r in rows})
    log.info("Exon boundaries computed: %d exons across %d proteins, %d skipped",
             len(rows), n_proteins, skipped)

    cols = ["Protein_ID", "exon_number", "total_exons", "aa_start", "aa_end",
            "aa_length", "genomic_start", "genomic_end"]
    exon_df = pd.DataFrame(rows, columns=cols)
    out_path = outdir / "exon.tsv"
    exon_df.to_csv(out_path, sep="\t", index=False)
    log.info("Written: %s", out_path)


if __name__ == "__main__":
    main()
