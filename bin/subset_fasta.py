#!/usr/bin/env python3
"""
subset_fasta.py — FASTA subsetting worker

Reads an input FASTA file (plain or gzip-compressed) and writes only those
sequences whose header line contains the given search string to the output file.

When --search is omitted or empty the entire input is copied unchanged
(pass-through mode), making this safe to call unconditionally in Nextflow.

Examples
--------
# Subset UniProt human proteome to RAF1 entry
subset_fasta.py \
    --input  UP000005640_9606.fasta \
    --output raf1_uniprot.fasta \
    --search GN=RAF1

# Subset GENCODE protein FASTA to all RAF1 transcripts  (pipe-field match)
subset_fasta.py \
    --input  gencode.v44.pc_translations.fa \
    --output raf1_gencode.fasta \
    --search '|RAF1|'

# Pass-through (no search term — copies everything)
subset_fasta.py --input all.fasta --output all_copy.fasta
"""

import argparse
import gzip
import logging
import sys
from pathlib import Path
from typing import Iterator

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FASTA I/O helpers
# ---------------------------------------------------------------------------

def open_fasta(path: str):
    """Return a text-mode file handle, transparently decompressing .gz files."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8")
    return open(p, "r", encoding="utf-8")


def fasta_records(fh) -> Iterator[tuple[str, list[str]]]:
    """
    Yield (header_line, sequence_lines) tuples.

    header_line includes the leading '>'.
    sequence_lines is the list of raw sequence lines (without newlines).
    """
    header: str | None = None
    seq_lines: list[str] = []

    for raw in fh:
        line = raw.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, seq_lines
            header = line
            seq_lines = []
        else:
            if line:
                seq_lines.append(line)

    if header is not None:
        yield header, seq_lines


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def build_patterns(gene_list: str, fasta_type: str) -> list[str]:
    """
    Build exact header search patterns from a comma-separated gene list.

    fasta_type 'uniprot' → "GN=GENE " (trailing space avoids GN=RAF1A)
    fasta_type 'gencode' or 'cdna' → "|GENE|" (pipe-delimited)
    """
    genes = [g.strip() for g in gene_list.split(",") if g.strip()]
    if fasta_type == "uniprot":
        return [f"GN={g} " for g in genes]
    else:
        return [f"|{g}|" for g in genes]


def subset_fasta(
    input_path: str,
    output_path: str,
    search: str,
    invert: bool = False,
    patterns: list[str] | None = None,
) -> tuple[int, int]:
    """
    Write matching sequences from *input_path* to *output_path*.

    Parameters
    ----------
    search : str
        Case-sensitive substring to look for in each sequence header.
        Empty string → include all sequences (pass-through).
        Ignored when *patterns* is provided.
    invert : bool
        If True, keep sequences whose header does NOT match.
    patterns : list[str] | None
        When provided, keep headers that contain ANY of these strings
        (OR logic).  Overrides *search*.

    Returns
    -------
    (total, written) — count of all sequences seen and sequences written.
    """
    if patterns:
        pass_through = False
        search_fn = lambda hdr: any(p in hdr for p in patterns)
    else:
        pass_through = not search
        search_fn = lambda hdr: search in hdr

    total = 0
    written = 0

    with open_fasta(input_path) as fh_in, \
         open(output_path, "w", encoding="utf-8") as fh_out:

        for header, seq_lines in fasta_records(fh_in):
            total += 1
            if pass_through:
                match = True
            else:
                match = search_fn(header)

            if invert:
                match = not match

            if match:
                fh_out.write(header + "\n")
                for sl in seq_lines:
                    fh_out.write(sl + "\n")
                written += 1

    return total, written


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Subset a FASTA file by a header search string"
    )
    p.add_argument(
        "--input",
        required=True,
        help="Input FASTA file (plain or .gz)",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output FASTA file (always plain text)",
    )
    p.add_argument(
        "--search",
        default="",
        help="Case-sensitive string to find in FASTA headers. "
             "Empty (default) = copy all sequences.",
    )
    p.add_argument(
        "--gene_list",
        default="",
        help="Comma-separated gene names (e.g. RAF1,TP53,BRAF). "
             "Overrides --search; builds exact patterns per --fasta_type.",
    )
    p.add_argument(
        "--fasta_type",
        default="gencode",
        choices=["uniprot", "gencode", "cdna"],
        help="Header format used when --gene_list is provided: "
             "'uniprot' → GN=GENE  patterns; "
             "'gencode'/'cdna' → |GENE| patterns.",
    )
    p.add_argument(
        "--invert",
        action="store_true",
        default=False,
        help="Invert match: keep sequences whose header does NOT contain --search",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    patterns: list[str] | None = None
    if args.gene_list:
        patterns = build_patterns(args.gene_list, args.fasta_type)
        log.info("Multi-gene subsetting '%s' → '%s'  (genes=%s, type=%s)",
                 args.input, args.output, args.gene_list, args.fasta_type)
        log.info("  Patterns: %s", patterns[:5])
    elif args.search:
        log.info("Subsetting '%s' → '%s'  (search='%s')",
                 args.input, args.output, args.search)
    else:
        log.info("Pass-through mode: copying '%s' → '%s'",
                 args.input, args.output)

    total, written = subset_fasta(
        input_path=args.input,
        output_path=args.output,
        search=args.search,
        invert=args.invert,
        patterns=patterns,
    )

    log.info("Done: %d / %d sequences written.", written, total)

    if (args.search or args.gene_list) and written == 0:
        log.warning(
            "No sequences matched in '%s'. "
            "Check your search term and FASTA headers.",
            args.input,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
