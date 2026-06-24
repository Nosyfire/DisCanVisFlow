#!/usr/bin/env python3
"""
create_genome_query_worker.py — Genome ↔ Protein reference tables.

Produces two query-friendly TSVs from combined_map.map so the DisCanVis2 web
layer can resolve genomic coordinates to protein level and, crucially, look up
or *construct* any point mutation.

Outputs
-------
  genome_protein_index.tsv
      One row per genomic nucleotide that is part of a CDS codon:
        chrom  gpos  strand  Protein_ID  prot_pos  codon_offset  aa  codon
      (prot_pos is 1-based.)

  genome_protein_mutations.tsv
      A *reference* table of EVERY possible single-nucleotide substitution in
      every CDS codon — not tied to any data source (ClinVar/TCGA/…).  This is
      the lookup that turns a genomic SNV (as found in a MAF/VCF) into its
      protein-level consequence, and conversely lets a desired protein mutation
      be realised as the genomic change(s) that produce it.
        chrom  gpos  strand  ref  alt  Protein_ID  prot_pos  codon_offset
        ref_codon  alt_codon  ref_aa  alt_aa  consequence  hgvs_g  hgvs_p
      ref/alt are on the genome (plus-strand) — i.e. VCF/MAF convention.

Inputs
------
  --combined_map  combined_map.map (Module 3 output)
  --outdir        output directory
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
_BASES = ("A", "C", "G", "T")

# Standard genetic code (DNA codons, T not U). '*' = stop.
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

INDEX_COLS = ["chrom", "gpos", "strand", "Protein_ID",
              "prot_pos", "codon_offset", "aa", "codon"]
MUT_COLS = ["chrom", "gpos", "strand", "ref", "alt", "Protein_ID", "prot_pos",
            "codon_offset", "ref_codon", "alt_codon", "ref_aa", "alt_aa",
            "consequence", "hgvs_g", "hgvs_p"]


def _translate(codon: str) -> str:
    return _CODON_TABLE.get(codon.upper(), "X")


def _consequence(ref_aa: str, alt_aa: str) -> str:
    if alt_aa == ref_aa:
        return "synonymous"
    if alt_aa == "*":
        return "nonsense"
    if ref_aa == "*":
        return "stop_loss"
    return "missense"


def parse_combined_map(map_path: str):
    """Yield per-codon records: (pid, chrom, strand, prot_pos, aa, codon, gpos_list)."""
    cur_pid = cur_chrom = cur_strand = None
    with open(map_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.lstrip("#").split()
                if len(parts) < 3:
                    cur_pid = None
                    continue
                header_fields = parts[0].split("|")
                cur_pid = next((f for f in header_fields if re.match(r".+-\d+$", f)),
                               header_fields[0])
                cur_chrom = parts[1]
                cur_strand = parts[2]
                continue
            if cur_pid is None:
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            try:
                prot_pos = int(cols[0]) + 1
                aa = cols[1]
                codon = cols[3]
                gpos_list = [g for g in cols[5].rstrip(",").split(",") if g and g != "-"]
            except (ValueError, IndexError):
                continue
            if len(codon) != 3 or len(gpos_list) != 3:
                continue
            yield cur_pid, cur_chrom, cur_strand, prot_pos, aa, codon, gpos_list


def build_tables(map_path: str):
    index_rows = []
    mut_rows = []
    proteins = set()

    for pid, chrom, strand, prot_pos, aa, codon, gpos_list in parse_combined_map(map_path):
        proteins.add(pid)
        minus = strand == "-"
        for offset in range(3):                       # 0,1,2 within codon (coding order)
            gpos = gpos_list[offset]
            coding_ref = codon[offset].upper()
            genomic_ref = _COMPLEMENT.get(coding_ref, coding_ref) if minus else coding_ref

            index_rows.append({
                "chrom": chrom, "gpos": gpos, "strand": strand,
                "Protein_ID": pid, "prot_pos": prot_pos,
                "codon_offset": offset + 1, "aa": aa, "codon": codon,
            })

            for galt in _BASES:                        # all 3 alternative genomic bases
                if galt == genomic_ref:
                    continue
                coding_alt = _COMPLEMENT[galt] if minus else galt
                alt_codon = codon[:offset] + coding_alt + codon[offset + 1:]
                alt_aa = _translate(alt_codon)
                mut_rows.append({
                    "chrom": chrom, "gpos": gpos, "strand": strand,
                    "ref": genomic_ref, "alt": galt,
                    "Protein_ID": pid, "prot_pos": prot_pos,
                    "codon_offset": offset + 1,
                    "ref_codon": codon.upper(), "alt_codon": alt_codon.upper(),
                    "ref_aa": aa, "alt_aa": alt_aa,
                    "consequence": _consequence(aa, alt_aa),
                    "hgvs_g": f"g.{gpos}{genomic_ref}>{galt}",
                    "hgvs_p": f"p.{aa}{prot_pos}{alt_aa}",
                })

    return index_rows, mut_rows, proteins


def main():
    ap = argparse.ArgumentParser(description="Build genome↔protein reference tables")
    ap.add_argument("--combined_map", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    index_rows, mut_rows, proteins = build_tables(args.combined_map)

    pd.DataFrame(index_rows, columns=INDEX_COLS).to_csv(
        outdir / "genome_protein_index.tsv", sep="\t", index=False)
    log.info("genome_protein_index.tsv: %d nucleotide rows (%d proteins)",
             len(index_rows), len(proteins))

    pd.DataFrame(mut_rows, columns=MUT_COLS).to_csv(
        outdir / "genome_protein_mutations.tsv", sep="\t", index=False)
    log.info("genome_protein_mutations.tsv: %d reference SNV rows (%d proteins)",
             len(mut_rows), len(proteins))


if __name__ == "__main__":
    main()
