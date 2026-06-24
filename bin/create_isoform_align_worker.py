#!/usr/bin/env python3
"""
create_isoform_align_worker.py — Insertion-free isoform alignment

For each gene, aligns every alternative isoform to the main (canonical)
isoform using global pairwise alignment (Needleman-Wunsch, BLOSUM62).
Positions where the *main* isoform has a gap in the alignment (i.e.,
insertions in the alt isoform relative to main) are stripped, producing
an "insertion-free" aligned sequence whose length equals the main isoform.

Output format mirrors GOPHER's insertion_free.tsv for orthologs:
    Protein_ID       — main isoform (e.g. RAF1-201)
    alt_Protein_ID   — alternative isoform (e.g. RAF1-205)
    gene             — gene name
    main_seq_len     — length of main isoform
    sequence         — alt isoform mapped to main-isoform positions
                       ('-' where the alt has a deletion relative to main)

Inputs
------
--seq_table       loc_chrom_with_names_isoforms_with_seq.tsv
                  (Module 2 output; used to get main isoform sequences)
--isoforms_table  loc_chrom_with_names_isoforms_only.tsv
                  (Module 2 optional output; alternative isoforms)
                  May be set to NO_FILE sentinel to skip alignment (only
                  self-alignment rows for each gene are written).
--outdir          output directory
"""

import argparse
import os
import sys

import pandas as pd
from tqdm import tqdm

try:
    from Bio.Align import PairwiseAligner
    from Bio.Align import substitution_matrices
    _BIOPYTHON_OK = True
except ImportError:
    _BIOPYTHON_OK = False


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def _build_aligner():
    """Return a configured Needleman-Wunsch aligner (BLOSUM62, affine gaps)."""
    if not _BIOPYTHON_OK:
        raise RuntimeError("Biopython not available; cannot perform alignment")
    aligner = PairwiseAligner()
    aligner.mode             = 'global'
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score   = -10
    aligner.extend_gap_score = -0.5
    return aligner


def insertion_free_align(main_seq: str, alt_seq: str, aligner) -> str:
    """
    Align alt_seq to main_seq (global) then strip all columns where
    main_seq has a gap character.  The resulting string has length
    len(main_seq): '-' at positions deleted in alt, residue otherwise.

    Returns the insertion-free alt sequence string.
    """
    if main_seq == alt_seq:
        return main_seq

    alignments = aligner.align(main_seq, alt_seq)
    best = next(iter(alignments))
    aligned_main, aligned_alt = best[0], best[1]

    # Strip positions where main has a gap (insertions in alt)
    result = []
    for m_aa, a_aa in zip(aligned_main, aligned_alt):
        if m_aa != '-':
            result.append(a_aa)

    aligned_str = ''.join(result)

    # Sanity check: must match main_seq length
    if len(aligned_str) != len(main_seq):
        # Fallback: return gaps of the right length
        return '-' * len(main_seq)

    return aligned_str


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Insertion-free isoform alignment against the main isoform'
    )
    parser.add_argument('--seq_table',      required=True,
                        help='loc_chrom_with_names_isoforms_with_seq.tsv')
    parser.add_argument('--isoforms_table', required=True,
                        help='loc_chrom_with_names_isoforms_only.tsv or NO_FILE')
    parser.add_argument('--outdir',         required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, 'isoform_alignment.tsv')

    # ── Load seq_table (main isoforms) ───────────────────────────────────────
    seq_df = pd.read_csv(args.seq_table, sep='\t', dtype=str)
    print(f'Seq table: {len(seq_df)} rows, columns: {list(seq_df.columns[:6])} ...')

    # Identify Protein_ID column (Transcript name or Protein_ID)
    pid_col = next(
        (c for c in ['Protein_ID', 'Transcript name'] if c in seq_df.columns),
        None
    )
    gene_col = next(
        (c for c in ['Gene_Uniprot', 'Gene_Gencode', 'Gene'] if c in seq_df.columns),
        None
    )

    if pid_col is None or 'Sequence' not in seq_df.columns:
        print('ERROR: seq_table missing Protein_ID/Transcript name or Sequence column')
        sys.exit(1)

    # Build main isoform lookup: Protein_ID → sequence, gene
    main_rows = seq_df
    if 'main_isoform' in seq_df.columns:
        main_rows = seq_df[seq_df['main_isoform'] == 'yes']
    # Fallback: use rows where Protein_ID == main_isoform_id
    elif 'main_isoform_id' in seq_df.columns:
        main_rows = seq_df[seq_df[pid_col] == seq_df['main_isoform_id']]

    main_seqs = {}        # Protein_ID → sequence
    main_gene  = {}       # Protein_ID → gene name
    for _, row in main_rows.iterrows():
        pid  = str(row[pid_col])
        seq  = str(row['Sequence']) if pd.notna(row['Sequence']) else ''
        gene = str(row[gene_col]) if gene_col and pd.notna(row[gene_col]) else ''
        if seq:
            main_seqs[pid] = seq
            main_gene[pid] = gene

    print(f'Main isoforms loaded: {len(main_seqs)}')

    # Only align to transcripts that were actually selected for the run (i.e. the
    # ones present in seq_table — for all_isoform_mapping these are the one-per-
    # UniProt-isoform winners). This avoids aligning to every GENCODE transcript.
    selected_pids = set(seq_df[pid_col].dropna().astype(str))
    print(f'Selected transcripts to align against: {len(selected_pids)}')

    # ── Load isoforms_only (alternative isoforms) ────────────────────────────
    no_file = (
        not os.path.isfile(args.isoforms_table)
        or os.path.basename(args.isoforms_table) == 'NO_FILE'
        or os.path.getsize(args.isoforms_table) == 0
    )

    alt_rows_by_main = {}  # main_Protein_ID → list of (alt_Protein_ID, alt_seq, gene)

    if not no_file:
        iso_df = pd.read_csv(args.isoforms_table, sep='\t', dtype=str)
        print(f'Isoforms table: {len(iso_df)} rows')

        alt_pid_col = next(
            (c for c in ['Protein_ID', 'Transcript name'] if c in iso_df.columns),
            None
        )
        alt_gene_col = next(
            (c for c in ['Gene_Uniprot', 'Gene_Gencode', 'Gene'] if c in iso_df.columns),
            None
        )

        if alt_pid_col and 'Sequence' in iso_df.columns and 'main_isoform_id' in iso_df.columns:
            for _, row in iso_df.iterrows():
                main_id = str(row['main_isoform_id']) if pd.notna(row['main_isoform_id']) else ''
                alt_id  = str(row[alt_pid_col]) if pd.notna(row[alt_pid_col]) else ''
                alt_seq = str(row['Sequence']) if pd.notna(row['Sequence']) else ''
                gene    = str(row[alt_gene_col]) if alt_gene_col and pd.notna(row[alt_gene_col]) else ''
                # Restrict to transcripts selected for the run (the UniProt-isoform
                # winners); skip every other GENCODE transcript.
                if alt_id not in selected_pids:
                    continue
                if main_id and alt_id and alt_seq and alt_seq != 'nan':
                    alt_rows_by_main.setdefault(main_id, []).append((alt_id, alt_seq, gene))
        else:
            print('WARNING: isoforms_table missing required columns — skipping alt isoforms')
    else:
        print('No isoforms_only table provided — only self-alignment rows will be written')

    # ── Align ─────────────────────────────────────────────────────────────────
    aligner = _build_aligner() if _BIOPYTHON_OK else None

    records = []

    for main_pid, main_seq in tqdm(main_seqs.items(), desc='Isoform alignment'):
        gene = main_gene.get(main_pid, '')
        main_len = len(main_seq)

        # Self-alignment row (main vs main)
        records.append({
            'Protein_ID':     main_pid,
            'alt_Protein_ID': main_pid,
            'gene':           gene,
            'main_seq_len':   main_len,
            'sequence':       main_seq,
        })

        alt_list = alt_rows_by_main.get(main_pid, [])
        for alt_pid, alt_seq, alt_gene in alt_list:
            if alt_pid == main_pid:
                continue  # already added self-alignment

            if aligner is None or not alt_seq:
                # No aligner: fill with gaps for positions, or skip
                aln_seq = '-' * main_len
            else:
                try:
                    aln_seq = insertion_free_align(main_seq, alt_seq, aligner)
                except Exception as exc:
                    print(f'WARNING: alignment failed for {alt_pid} vs {main_pid}: {exc}',
                          file=sys.stderr)
                    aln_seq = '-' * main_len

            records.append({
                'Protein_ID':     main_pid,
                'alt_Protein_ID': alt_pid,
                'gene':           gene or alt_gene,
                'main_seq_len':   main_len,
                'sequence':       aln_seq,
            })

    out_df = pd.DataFrame(records, columns=[
        'Protein_ID', 'alt_Protein_ID', 'gene', 'main_seq_len', 'sequence'
    ])
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f'Written {len(out_df)} rows → {out_path}')
    n_genes = out_df['Protein_ID'].nunique()
    n_alts  = len(out_df) - n_genes  # subtract self-alignment rows
    print(f'  {n_genes} genes, {n_alts} alternative isoform alignments')


if __name__ == '__main__':
    main()
