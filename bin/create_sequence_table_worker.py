#!/usr/bin/env python3
"""
create_sequence_table_worker.py — Module 2: Sequence Process

Merges BLAST ID-map results with GENCODE GTF annotations and protein
sequences to produce the canonical transcript/isoform table used by
downstream modules.

Inputs
------
--blast_best        bestmaps_blast_gene_transcript.tsv   (Module 1 output)
--blast_isoforms    blastmaps_isoforms.tsv               (Module 1, optional)
--gtf               GENCODE annotation GTF (.gz)
--translations      GENCODE pc_translations FASTA        (for sequences)
--output_dir        output directory
--cutoff            minimum coverage % for above_cutoff flag (default: 80)

Outputs (all tab-separated, in --output_dir)
-------------------------------------------
loc_chrom_with_names.tsv                 — all mapped transcripts
loc_chrom_with_names_isoforms_with_seq.tsv — above + protein sequences
                                             + main_isoform flag
"""

import argparse
import gzip
import io
import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def _open_maybe_gz(path: str):
    """Return a text-mode file handle, transparent to gzip."""
    if path.endswith('.gz'):
        return io.TextIOWrapper(gzip.open(path, 'rb'), encoding='utf-8')
    return open(path, 'r', encoding='utf-8')


def parse_gtf(gtf_path: str) -> pd.DataFrame:
    """
    Parse a GENCODE GTF (plain or .gz) and return a DataFrame of transcript-
    level annotations with isoform flag columns.

    Output columns:
        Transcript ID, Chromosome, gene_type,
        Ensembl_canonical, MANE_Select, appris_principal_1
    """
    rows = []
    attr_re_map = {
        'gene_id':         'gene_id',
        'gene_name':       'gene_name',
        'gene_type':       'gene_type',
        'transcript_id':   'Transcript ID',
        'transcript_name': 'transcript_name',
    }
    flag_tags = ('Ensembl_canonical', 'MANE_Select', 'appris_principal_1')

    with _open_maybe_gz(gtf_path) as fh:
        for line in fh:
            # Skip GTF header/comment lines
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature != 'transcript':
                continue

            chrom   = cols[0]
            info    = cols[8]

            row = {'Chromosome': chrom}
            for key, col_name in attr_re_map.items():
                import re
                m = re.search(rf'{key} "([^"]*)"', info)
                row[col_name] = m.group(1) if m else None

            for tag in flag_tags:
                row[tag] = 'yes' if f'tag "{tag}"' in info else 'no'

            rows.append(row)

    df = pd.DataFrame(rows)
    df.dropna(subset=['Transcript ID'], inplace=True)
    df.drop_duplicates(inplace=True)
    return df


# ---------------------------------------------------------------------------
# Translations FASTA  →  transcript_id → sequence
# ---------------------------------------------------------------------------

def load_translations(fasta_path: str) -> dict:
    """
    Parse the GENCODE pc_translations FASTA.

    Header format:
        >ENST00000423430.6|ENSG00000132155.11|...|RAF1-201|RAF1|648

    Returns a dict: transcript_id (ENST.version) → protein sequence string.
    """
    seqs = {}
    current_id = None
    current_seq_parts = []

    def _flush():
        if current_id is not None:
            seqs[current_id] = ''.join(current_seq_parts)

    opener = gzip.open if fasta_path.endswith('.gz') else open
    mode   = 'rt' if fasta_path.endswith('.gz') else 'r'

    with opener(fasta_path, mode) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                _flush()
                # Header can be ENST|ENSG|... (older GENCODE) or
                # ENSP|ENST|ENSG|... (newer GENCODE v44+).
                # We always use the ENST field as the key.
                parts = line[1:].split('|')
                enst  = next((p for p in parts if p.startswith('ENST')), parts[0])
                current_id = enst
                current_seq_parts = []
            else:
                current_seq_parts.append(line)
    _flush()
    return seqs


# ---------------------------------------------------------------------------
# Isoform identification  (ported from create_loc_chrom_with_names.py)
# ---------------------------------------------------------------------------

def isoform_identification(df: pd.DataFrame,
                           gene_col: str = 'Gene_Gencode') -> pd.DataFrame:
    """
    For each gene (Gene_Gencode), select one 'main isoform' using a
    priority ladder:
      1. identical + MANE_Select + Ensembl_canonical + appris_principal_1
      2. identical + any of (MANE | Ensembl | appris)
      3. identical (alone)
      4. MANE + Ensembl + appris (all three)
      5. any two of (MANE | Ensembl | appris)
      6. any one of (MANE | Ensembl | appris)
      7. highest coverage (fallback)

    Adds columns: main_isoform ('yes'/'no'), main_isoform_id
    """
    df = df.copy()
    df['main_isoform']    = 'no'
    df['main_isoform_id'] = ''

    main_isoform_ids: dict = {}

    for gene in tqdm(df[gene_col].unique(), desc='Main isoform selection'):
        sub = df[df[gene_col] == gene]

        if len(sub) == 1:
            idx = sub.index[0]
            df.loc[idx, 'main_isoform']    = 'yes'
            df.loc[idx, 'main_isoform_id'] = sub.at[idx, 'Transcript name']
            main_isoform_ids[gene]         = sub.at[idx, 'Transcript name']
            continue

        identical = sub['alignmentpuntcuality'] == 'identical'
        ensembl   = sub['Ensembl_canonical']   == 'yes'
        mane      = sub['MANE_Select']         == 'yes'
        appris    = sub['appris_principal_1']  == 'yes'

        cond_all  = identical & ensembl & mane & appris
        cond_id_c = identical & (ensembl | mane | appris)
        cond_id   = identical
        cond_abc  = ensembl & mane & appris
        cond_any2 = (ensembl.astype(int) + mane.astype(int) + appris.astype(int)) >= 2
        cond_any1 = ensembl | mane | appris

        if cond_all.any():
            chosen = sub[cond_all]
        elif cond_id_c.any():
            chosen = sub[cond_id_c]
        elif cond_id.any():
            chosen = sub[cond_id]
        elif cond_abc.any():
            chosen = sub[cond_abc]
        elif cond_any2.any():
            chosen = sub[cond_any2]
        elif cond_any1.any():
            chosen = sub[cond_any1]
        else:
            chosen = sub.sort_values('coverage', ascending=False).head(1)

        if not chosen.empty:
            winner = chosen.sort_values('coverage', ascending=False).iloc[0]
            df.loc[winner.name, 'main_isoform']    = 'yes'
            df.loc[winner.name, 'main_isoform_id'] = winner['Transcript name']
            main_isoform_ids[gene]                 = winner['Transcript name']

    # Propagate main_isoform_id to all isoforms within the same gene
    for gene, main_id in tqdm(main_isoform_ids.items(),
                               desc='Propagate main_isoform_id'):
        df.loc[df[gene_col] == gene, 'main_isoform_id'] = main_id

    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Module 2: Sequence Process — merge BLAST results with '
                    'GENCODE GTF annotations and protein sequences.'
    )
    parser.add_argument('--blast_best',     required=True,
                        help='bestmaps_blast_gene_transcript.tsv (Module 1 output)')
    parser.add_argument('--blast_isoforms', default=None,
                        help='blastmaps_isoforms.tsv (optional)')
    parser.add_argument('--gtf',            required=True,
                        help='GENCODE annotation GTF (.gz or plain)')
    parser.add_argument('--translations',   required=True,
                        help='GENCODE pc_translations FASTA (.gz or plain)')
    parser.add_argument('--output_dir',     required=True)
    parser.add_argument('--cutoff',         type=int, default=80,
                        help='Coverage %% threshold for above_cutoff flag')
    parser.add_argument('--mapping_mode',   default='main_isoform_mapping',
                        choices=['main_isoform_mapping', 'all_isoform_mapping'],
                        help='all_isoform_mapping keeps one winner transcript per '
                             'UniProt isoform (1:1); main_isoform_mapping keeps all '
                             'matched transcripts and flags one main per gene.')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── 1. Parse GTF ─────────────────────────────────────────────────────────
    print('Parsing GENCODE GTF ...')
    gtf_df = parse_gtf(args.gtf)
    gtf_df = gtf_df[['Transcript ID', 'Chromosome', 'gene_type',
                      'Ensembl_canonical', 'MANE_Select', 'appris_principal_1']]
    gtf_df.drop_duplicates(inplace=True)
    print(f'  GTF transcripts: {len(gtf_df):,}')

    # ── 2. Read BLAST best-hit table ──────────────────────────────────────────
    print('Reading BLAST best-hit table ...')
    blast_df = pd.read_csv(args.blast_best, sep='\t', header=0)
    print(f'  BLAST rows: {len(blast_df):,}   genes: {blast_df["Gene_Gencode"].nunique():,}')

    # ── 3. Merge on Transcript ID ─────────────────────────────────────────────
    merged = blast_df.merge(gtf_df, on='Transcript ID', how='left')

    # Save diagnostics for transcripts without full annotation
    no_annot = merged[merged.isna().any(axis=1)].copy()
    if not no_annot.empty:
        no_annot_path = os.path.join(args.output_dir, 'transcripts_missing_annotation.tsv')
        no_annot.to_csv(no_annot_path, sep='\t', index=False)
        print(f'  Transcripts without full annotation: {len(no_annot):,}  '
              f'(saved to {no_annot_path})')

    merged.dropna(inplace=True)
    print(f'  After dropna: {len(merged):,} rows, '
          f'{merged["Gene_Gencode"].nunique():,} genes')

    # ── 4. Isoform selection ──────────────────────────────────────────────────
    if args.mapping_mode == 'all_isoform_mapping':
        # One winner per *UniProt isoform*: group by Entry_Isoform, pick the best
        # representative transcript for each isoform, and KEEP ONLY those winners.
        # RAF1 (P04049 canonical + P04049-2) → exactly 2 rows.
        print('Selecting one winner per UniProt isoform (all_isoform_mapping) ...')
        merged = isoform_identification(merged, gene_col='Entry_Isoform')
        merged = merged[merged['main_isoform'] == 'yes'].copy()
        # The gene's main isoform is the winner of the canonical accession
        # (no '-N' suffix); alternative-isoform winners are kept but flagged 'no'.
        is_canon = ~merged['Entry_Isoform'].astype(str).str.contains('-', na=False)
        merged['main_isoform'] = np.where(is_canon, 'yes', 'no')
        canon_map = (merged[is_canon]
                     .drop_duplicates('Gene_Gencode')
                     .set_index('Gene_Gencode')['Transcript name'].to_dict())
        merged['main_isoform_id'] = merged['Gene_Gencode'].map(canon_map)
        merged['main_isoform_id'] = merged['main_isoform_id'].fillna(merged['Transcript name'])
        print(f'  Kept {len(merged):,} isoform winners '
              f'({merged["Entry_Isoform"].nunique():,} UniProt isoforms)')
    else:
        print('Selecting main isoforms ...')
        merged = isoform_identification(merged)
    merged['above_cutoff'] = merged['coverage'] >= args.cutoff

    # ── 5. Add protein sequences ──────────────────────────────────────────────
    print('Loading protein sequences from translations FASTA ...')
    seq_map = load_translations(args.translations)
    print(f'  Loaded {len(seq_map):,} sequences')

    merged['Sequence'] = merged['Transcript ID'].map(seq_map)
    missing_seq = merged['Sequence'].isna().sum()
    if missing_seq:
        print(f'  WARNING: {missing_seq:,} transcripts have no sequence '
              f'in translations FASTA')

    # ── 6. Write outputs ──────────────────────────────────────────────────────
    # Add Protein_ID (= Transcript name, e.g. "RAF1-201") and Gene aliases
    # so that downstream workers can use the same column names as the legacy pipeline.
    if 'Protein_ID' not in merged.columns and 'Transcript name' in merged.columns:
        merged.insert(merged.columns.get_loc('Transcript name') + 1,
                      'Protein_ID', merged['Transcript name'])
    if 'Gene' not in merged.columns:
        gene_src = next((c for c in ['Gene_Gencode', 'Gene_Uniprot'] if c in merged.columns), None)
        if gene_src:
            merged.insert(merged.columns.get_loc(gene_src) + 1, 'Gene', merged[gene_src])

    # Base table (no sequences)
    base_path = os.path.join(args.output_dir, 'loc_chrom_with_names.tsv')
    merged.drop(columns=['Sequence']).to_csv(base_path, sep='\t', index=False)
    print(f'Written: {base_path}')

    # Full table with sequences
    seq_path = os.path.join(args.output_dir,
                            'loc_chrom_with_names_isoforms_with_seq.tsv')
    merged.to_csv(seq_path, sep='\t', index=False)
    print(f'Written: {seq_path}')

    # Main-isoform-only table
    main_path = os.path.join(args.output_dir,
                             'loc_chrom_with_names_main_isoform.tsv')
    main_df = merged[merged['main_isoform'] == 'yes'].copy()
    main_df.to_csv(main_path, sep='\t', index=False)
    print(f'Written: {main_path}  ({len(main_df):,} rows)')

    # Optional isoforms table (if provided)
    if args.blast_isoforms:
        iso_df = pd.read_csv(args.blast_isoforms, sep='\t', header=0)
        iso_merged = iso_df.merge(gtf_df, on='Transcript ID', how='left')
        iso_merged['Sequence'] = iso_merged['Transcript ID'].map(seq_map)

        # Propagate main_isoform_id from main table
        main_id_map = main_df.set_index('Gene_Gencode')['main_isoform_id'].to_dict()
        iso_merged['main_isoform_id'] = iso_merged['Gene_Gencode'].map(main_id_map)

        iso_path = os.path.join(args.output_dir,
                                'loc_chrom_with_names_isoforms_only.tsv')
        iso_merged.to_csv(iso_path, sep='\t', index=False)
        print(f'Written: {iso_path}  ({len(iso_merged):,} rows)')

    print('Module 2 complete.')


if __name__ == '__main__':
    main()
