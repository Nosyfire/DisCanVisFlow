#!/usr/bin/env python3
"""
create_genome_map_worker.py — Module 3: Genome Mapping

For each transcript that has a BLAT alignment in the PSL file, builds a
per-residue coordinate map:

    protein AA index  →  cDNA codon positions  →  genomic positions

Output format (combined_map.map):
    One block per transcript:
        # Qname Tname strand Tstart-Tend
        0 M 0,1,2 ATG M    12600000,12600001,12600002, ATG M
        1 A 3,4,5 GCT A    12600003,12600004,12600005, GCT A
        ...

Ported from the legacy mapping_64.py + remap.py with these improvements:
  • BioPython SeqIO for FASTA access  (no custom index files needed)
  • Modern tblastn arguments         (BLAST+ 2.13+, not legacy bl2seq)
  • Streamlined PSL parsing          (column names hardcoded from standard PSL)

Inputs
------
--psl           combined_output.psl   (from BLAT_ALIGN)
--cdna_fasta    GENCODE pc_transcripts FASTA  (plain or .gz)
--prot_fasta    GENCODE pc_translations FASTA (plain or .gz)
--loc_chrom     loc_chrom_with_names.tsv      (from SEQUENCE_PROCESS)
--hg38_2bit     hg38.2bit
--output_dir    output directory
--num_processes number of parallel workers   (default: 4)
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Codon → amino acid table
# ---------------------------------------------------------------------------
CODON_TABLE = {
    'TTT':'F','TTC':'F',
    'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I',
    'ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'_','TAG':'_',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'_','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def codon2aa(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), 'X')


# ---------------------------------------------------------------------------
# PSL parsing
# ---------------------------------------------------------------------------
PSL_COLUMNS = [
    'matches','misMatches','repMatches','nCount',
    'qNumInsert','qBaseInsert','tNumInsert','tBaseInsert',
    'strand',
    'Qname','Qsize','Qstart','Qend',
    'Tname','Tsize','Tstart','Tend',
    'blockcount','blockSizes','qStarts','tStarts',
]

def load_psl(psl_path: str) -> pd.DataFrame:
    """Read a PSL file and return a DataFrame.  Skips BLAT header lines."""
    rows = []
    with open(psl_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#') or 'psLayout' in line \
                    or line.startswith('match') or line.startswith('-'):
                continue
            parts = line.split('\t')
            if len(parts) < 21:
                continue
            rows.append(parts[:21])
    df = pd.DataFrame(rows, columns=PSL_COLUMNS)
    for col in ['Qsize','Qstart','Qend','Tsize','Tstart','Tend','blockcount',
                'matches','misMatches']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['transcript_id'] = df['Qname'].str.split('|').str[0]
    return df


# ---------------------------------------------------------------------------
# FASTA loading (BioPython-free, memory-efficient dict)
# ---------------------------------------------------------------------------
def _open_fasta(path: str):
    if path.endswith('.gz'):
        import io
        return io.TextIOWrapper(gzip.open(path, 'rb'), encoding='utf-8')
    return open(path, 'r')


def load_fasta_dict(path: str, key_field: int = 0) -> dict:
    """
    Parse a GENCODE FASTA and return a dict:
        transcript_id  →  {'header': str, 'seq': str, 'fasta': str}

    key_field : pipe-split field index used as the lookup key.
                0  → cDNA (field 0 = ENST transcript ID)
                1  → protein (field 0 = ENSP protein ID; field 1 = ENST transcript ID)
    """
    seqs = {}
    current_key = None
    current_seq = []
    current_hdr = None

    def _flush():
        if current_key:
            seq = ''.join(current_seq).replace('*', '')
            seqs[current_key] = {
                'header': current_hdr,
                'seq':    seq,
                'fasta':  f'{current_hdr}\n{seq}\n',
            }

    with _open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                _flush()
                current_hdr = line
                current_key = line[1:].split('|')[key_field]
                current_seq = []
            else:
                current_seq.append(line)
    _flush()
    return seqs


# ---------------------------------------------------------------------------
# tblastn:  protein → cDNA alignment  (prot2cdna map)
# ---------------------------------------------------------------------------
def _parse_tblastn_output(blast_output: str,
                           prot_len: int) -> list:
    """
    Parse tblastn text output (format 0) to build a cdnamap list.

    Returns list of length prot_len where each element is either:
        None                          → residue not mapped
        {'0':ii,'1':ii+1,'2':ii+2, 'codon':str, 'aa':str}
    """
    # Find best-score hit
    hits = blast_output.split(' Score = ')
    best_hit = None
    best_score = -1.0
    score_re = re.compile(r'([\d.]+(?:e[+-]?\d+)?)\s+bits')
    for hit in hits[1:]:
        m = score_re.search(hit)
        if m:
            score = float(m.group(1))
            if score > best_score:
                best_score, best_hit = score, hit

    if best_hit is None:
        return [None] * prot_len

    hit_text = f' Score = {best_hit}'

    frame = ''
    Qbeg = Tbeg = 0
    seq_prot = seq_cdna = ''
    Tend = 0
    on = -1

    def _parse_al_line(line):
        # tblastn format 0 uses spaces (not a colon) before the position number:
        #   Query  1   MEHIQ...  60
        #   Sbjct  402 MEHIQ...  581
        # Older BLAST versions occasionally used a colon, so we allow ":?" here.
        m = re.match(r'(Sbjct|Query):?\s+(\d+)\s+(\S+(?:\s+\S+)*?)\s+(\d+)\s*$', line.strip())
        if m:
            return m.group(2), m.group(3), m.group(4)
        return None, '', None

    for line in hit_text.split('\n'):
        if 'Frame' in line:
            frame = line.split()[-1]
        if on > 1:
            break
        if 'Score' in line:
            on += 1
        if frame:
            if 'Sbjct' in line:
                start, seq, end = _parse_al_line(line)
                if start:
                    seq_prot += seq
                    if Qbeg == 0:
                        Qbeg = int(start)
                    Tend = int(end) if end else Tend
            if 'Query' in line:
                start, seq, end = _parse_al_line(line)
                if start:
                    seq_cdna += seq
                    if Tbeg == 0:
                        Tbeg = int(start)

    Qbeg -= 1   # 0-indexed cDNA nucleotide start
    Tbeg -= 1   # 0-indexed protein residue start

    # Walk through the alignment character-by-character, handling gaps.
    # seq_cdna  = Query  line = protein AA sequence (1 char per residue)
    # seq_prot  = Sbjct  line = translated cDNA AA sequence (1 char per codon)
    # Each non-gap cDNA character advances cdna_pos by 3 nucleotides.
    # Codon strings are filled later from the actual cDNA sequence in _worker.
    cdnamap  = [None] * prot_len
    prot_i   = Tbeg
    cdna_pos = Qbeg

    for p_aa, s_aa in zip(seq_cdna, seq_prot):
        p_gap = (p_aa == '-')
        s_gap = (s_aa == '-')
        if p_gap and s_gap:
            continue
        elif p_gap:
            cdna_pos += 3            # insertion in cDNA vs protein
        elif s_gap:
            prot_i   += 1            # deletion in cDNA (protein residue unmapped)
        else:
            if 0 <= prot_i < prot_len:
                cdnamap[prot_i] = {'0': cdna_pos, '1': cdna_pos + 1,
                                   '2': cdna_pos + 2,
                                   'codon': '---', 'aa': '?'}
            prot_i   += 1
            cdna_pos += 3

    return cdnamap


def map_prot_cdna(cdna_fasta_str: str,
                  prot_fasta_str: str,
                  prot_len: int) -> list:
    """
    Run tblastn to align a protein against a cDNA and build the cdnamap.
    Temporary files are used and cleaned up afterwards.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tf_cdna, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tf_prot, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tf_out:

        tf_cdna.write(cdna_fasta_str)
        tf_prot.write(prot_fasta_str)
        cdna_path = tf_cdna.name
        prot_path = tf_prot.name
        out_path  = tf_out.name

    try:
        cmd = [
            'tblastn',
            '-query',            prot_path,
            '-subject',          cdna_path,
            '-out',              out_path,
            '-outfmt',           '0',
            '-num_alignments',   '1',
            '-num_descriptions', '1',
            '-seg',              'no',
        ]
        subprocess.run(cmd, check=False, capture_output=True)

        with open(out_path) as fh:
            blast_text = fh.read()

        return _parse_tblastn_output(blast_text, prot_len)
    finally:
        for p in [cdna_path, prot_path, out_path]:
            try:
                os.remove(p)
            except FileNotFoundError:
                pass


# ---------------------------------------------------------------------------
# cDNA → genome alignment  (from PSL blocks)
# ---------------------------------------------------------------------------
def _int_list(s: str) -> list:
    return [int(x) for x in s.strip().rstrip(',').split(',') if x]


def ali_cdna_chrom(psl_row: dict,
                   cdna_seq: str,
                   hg38_2bit: str) -> tuple:
    """
    Given a PSL row (dict) for one transcript, extract the genomic sequence
    via twoBitToFa and build a positional map:
        gpos[cdna_pos]  → relative genomic index
    Returns (gpos, gene_seq, Qstart, Qend).
    """
    strand     = psl_row['strand']
    Tname      = psl_row['Tname']
    Tstart     = int(psl_row['Tstart'])
    Tend       = int(psl_row['Tend'])
    Qsize      = int(psl_row['Qsize'])
    Qstart     = int(psl_row['Qstart'])
    Qend       = int(psl_row['Qend'])
    blockcount = int(psl_row['blockcount'])
    blockSizes = _int_list(str(psl_row['blockSizes']))
    qStarts    = _int_list(str(psl_row['qStarts']))
    tStarts    = _int_list(str(psl_row['tStarts']))

    # Extract genomic sequence with twoBitToFa
    with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        cmd = [
            'twoBitToFa', hg38_2bit,
            f'-seq={Tname}',
            f'-start={Tstart}',
            f'-end={Tend}',
            tmp_path,
        ]
        subprocess.run(cmd, check=False, capture_output=True)

        gene_seq = ''
        with open(tmp_path) as fh:
            for line in fh:
                if not line.startswith('>'):
                    gene_seq += line.strip()
    finally:
        try:
            os.remove(tmp_path)
        except FileNotFoundError:
            pass

    if strand == '-':
        gene_seq = gene_seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))

    # Normalise tStarts relative to Tstart
    min_t = min(tStarts)
    tStarts = [ts - min_t for ts in tStarts]

    cdna = list(cdna_seq)
    gpos = [None] * len(cdna)

    if strand == '-':
        qStarts = [Qsize - qStarts[i] - blockSizes[i]
                   for i in range(blockcount)][::-1]
        tStarts = [tStarts[i] + blockSizes[i] - 1
                   for i in range(blockcount)][::-1]
        blockSizes = blockSizes[::-1]
        step = -1
    else:
        step = 1

    for i in range(blockcount):
        for j in range(blockSizes[i]):
            q_pos = qStarts[i] + j
            t_pos = tStarts[i] + step * j
            if 0 <= q_pos < len(gpos):
                gpos[q_pos] = t_pos

    out_gene = [None] * (Qend - Qstart)
    for i in range(Qstart, Qend):
        try:
            if gpos[i] is not None:
                out_gene[i] = gpos[i]
        except IndexError:
            break

    return out_gene, gene_seq, Qstart, Qend, Tstart


# ---------------------------------------------------------------------------
# Output text generation
# ---------------------------------------------------------------------------
def create_text(psl_row: dict,
                prot_seq: str,
                cdnamap: list,
                out_gene: list,
                gene_seq: str,
                Tstart: int) -> str:
    L    = len(prot_seq)
    prot = list(prot_seq)
    gene = list(gene_seq)

    # Build genemap: prot index → {0,1,2} cDNA positions → genomic offset
    genemap = [None] * L
    for i in range(L):
        if cdnamap[i] is not None:
            gm = {}
            for k_str in ('0', '1', '2'):
                cdna_pos = cdnamap[i].get(k_str)
                if cdna_pos is not None and cdna_pos < len(out_gene):
                    gm[k_str] = out_gene[cdna_pos]
                else:
                    gm[k_str] = None
            genemap[i] = gm

    Qname = psl_row['Qname']
    Tname = psl_row['Tname']
    strand = psl_row['strand']

    fout = (f'# {Qname} {Tname} {strand} '
            f'{psl_row["Tstart"]}-{psl_row["Tend"]}\n')

    for i in range(L):
        fout += f'{i} {prot[i]} '

        if cdnamap[i] is not None and '0' in cdnamap[i]:
            c = cdnamap[i]
            fout += f"{c['0']},{c['1']},{c['2']} {c['codon']} {c['aa']} "
        else:
            fout += '-,-,-, --- - '

        fout += '    '
        codon = ''
        badcodon = False

        for k in range(3):
            if genemap[i] and str(k) in genemap[i] \
                    and genemap[i][str(k)] is not None:
                gval = genemap[i][str(k)]
                fout += f'{gval + Tstart},'
                try:
                    codon += gene[gval]
                except IndexError:
                    codon += '-'
                    badcodon = True
            else:
                fout += '-,'
                codon += '-'
                badcodon = True

        fout += f' {codon} '
        if not badcodon:
            fout += f'{codon2aa(codon)} '
        else:
            fout += '- '
        fout += '\n'

    return fout


# ---------------------------------------------------------------------------
# Per-transcript worker  (called by ProcessPoolExecutor)
# ---------------------------------------------------------------------------
def _worker(args: tuple) -> str:
    (psl_row, cdna_rec, prot_rec, hg38_2bit) = args
    try:
        cdna_seq = cdna_rec['seq']
        prot_seq = prot_rec['seq']

        # 1. Protein → cDNA alignment (tblastn)
        cdnamap = map_prot_cdna(
            cdna_fasta_str=cdna_rec['fasta'],
            prot_fasta_str=prot_rec['fasta'],
            prot_len=len(prot_seq),
        )

        # Fill actual codon strings from the cDNA sequence using position info.
        # _parse_tblastn_output stores '---' placeholders; replace them here.
        for entry in cdnamap:
            if entry is not None:
                pos = entry['0']
                if 0 <= pos and pos + 3 <= len(cdna_seq):
                    codon = cdna_seq[pos:pos + 3]
                    entry['codon'] = codon
                    entry['aa']    = codon2aa(codon)

        # 2. cDNA → genome alignment (PSL blocks + twoBitToFa)
        out_gene, gene_seq, Qstart, Qend, Tstart = ali_cdna_chrom(
            psl_row, cdna_seq, hg38_2bit
        )

        return create_text(psl_row, prot_seq, cdnamap,
                           out_gene, gene_seq, Tstart)
    except Exception as exc:
        tb = traceback.format_exc()
        return f'# Error in {psl_row["Qname"]}: {exc}\n{tb}\n'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Module 3: Genome Mapping — '
                    'protein AA → cDNA → genomic coordinate map.'
    )
    parser.add_argument('--psl',           required=True,
                        help='BLAT output PSL file (combined_output.psl)')
    parser.add_argument('--cdna_fasta',    required=True,
                        help='GENCODE pc_transcripts FASTA (cDNA)')
    parser.add_argument('--prot_fasta',    required=True,
                        help='GENCODE pc_translations FASTA (protein)')
    parser.add_argument('--loc_chrom',     required=True,
                        help='loc_chrom_with_names.tsv  (from Module 2)')
    parser.add_argument('--hg38_2bit',     required=True,
                        help='Path to hg38.2bit genome file')
    parser.add_argument('--output_dir',    required=True)
    parser.add_argument('--num_processes', type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Load PSL ─────────────────────────────────────────────────────────────
    print('Loading PSL file ...')
    psl_df = load_psl(args.psl)
    print(f'  PSL alignments: {len(psl_df):,}')

    # ── Load loc_chrom to get transcript list + Chromosome map ───────────────
    print('Loading loc_chrom TSV ...')
    loc_df = pd.read_csv(args.loc_chrom, sep='\t')
    transcript_list = loc_df['Transcript ID'].unique().tolist()
    print(f'  Transcripts in loc_chrom: {len(transcript_list):,}')

    # Filter PSL to transcripts in our table
    psl_df = psl_df[psl_df['transcript_id'].isin(transcript_list)]
    print(f'  PSL rows after filter: {len(psl_df):,}')

    # Keep only best (highest matches) PSL row per transcript
    psl_df = psl_df.sort_values('matches', ascending=False) \
                   .drop_duplicates(subset='transcript_id')

    # ── Load FASTA dicts ──────────────────────────────────────────────────────
    print('Loading cDNA FASTA ...')
    cdna_dict = load_fasta_dict(args.cdna_fasta, key_field=0)
    print(f'  cDNA entries: {len(cdna_dict):,}')

    print('Loading protein FASTA ...')
    # GENCODE protein headers: >ENSP...|ENST...|ENSG...|...
    # field[1] is the transcript ID (ENST) used as the join key
    prot_dict = load_fasta_dict(args.prot_fasta, key_field=1)
    print(f'  Protein entries: {len(prot_dict):,}')

    # ── Build work list ───────────────────────────────────────────────────────
    work = []
    no_cdna = []
    no_prot = []
    no_psl  = []

    for tid in tqdm(transcript_list, desc='Preparing work'):
        psl_rows = psl_df[psl_df['transcript_id'] == tid]
        if psl_rows.empty:
            no_psl.append(tid)
            continue
        psl_row = psl_rows.iloc[0].to_dict()

        cdna_rec = cdna_dict.get(tid)
        prot_rec = prot_dict.get(tid)

        if cdna_rec is None:
            no_cdna.append(tid)
            continue
        if prot_rec is None:
            no_prot.append(tid)
            continue

        work.append((psl_row, cdna_rec, prot_rec, args.hg38_2bit))

    print(f'  Work items ready: {len(work):,}')
    print(f'  No PSL:  {len(no_psl):,}   No cDNA: {len(no_cdna):,}'
          f'   No prot: {len(no_prot):,}')

    # ── Process in parallel ───────────────────────────────────────────────────
    map_path   = os.path.join(args.output_dir, 'combined_map.map')
    error_path = os.path.join(args.output_dir, 'error_map.txt')

    mapped = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=args.num_processes) as pool, \
         open(map_path, 'w') as mf, \
         open(error_path, 'w') as ef:

        futures = {pool.submit(_worker, item): item for item in work}

        for future in tqdm(as_completed(futures),
                           total=len(futures),
                           desc='Genome mapping'):
            result = future.result()
            if '# Error in' in result:
                ef.write(result)
                errors += 1
            else:
                mf.write(result)
                mapped += 1
            mf.flush()
            ef.flush()

    print(f'Mapped: {mapped:,}   Errors: {errors:,}   Total: {len(work):,}')
    print(f'Written: {map_path}')


if __name__ == '__main__':
    main()
