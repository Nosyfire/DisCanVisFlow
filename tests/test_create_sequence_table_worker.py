"""
Tests for bin/create_sequence_table_worker.py (Module 2 — Sequence Process)
"""
import os
import sys
import gzip
import tempfile
import subprocess

import pandas as pd
import pytest

# Add bin/ to path for direct import
BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, BIN_DIR)

DUMMY_DIR = os.path.join(os.path.dirname(__file__), 'dummy_data')
DUMMY_GTF        = os.path.join(DUMMY_DIR, 'dummy_gencode.gtf.gz')
DUMMY_BESTMAPS   = os.path.join(DUMMY_DIR, 'dummy_bestmaps.tsv')
DUMMY_TRANS      = os.path.join(DUMMY_DIR, 'dummy_translations.fasta')


# ── Import module functions ───────────────────────────────────────────────────

from create_sequence_table_worker import (
    parse_gtf,
    load_translations,
    isoform_identification,
)


# ── GTF parsing ───────────────────────────────────────────────────────────────

class TestParseGtf:
    def test_returns_dataframe(self):
        df = parse_gtf(DUMMY_GTF)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        df = parse_gtf(DUMMY_GTF)
        for col in ['Transcript ID', 'Chromosome', 'gene_type',
                    'Ensembl_canonical', 'MANE_Select', 'appris_principal_1']:
            assert col in df.columns, f'Missing column: {col}'

    def test_transcript_count(self):
        df = parse_gtf(DUMMY_GTF)
        assert len(df) == 3, f'Expected 3 transcripts, got {len(df)}'

    def test_raf1_201_flags(self):
        df = parse_gtf(DUMMY_GTF)
        row = df[df['Transcript ID'] == 'ENST00000423430.6'].iloc[0]
        assert row['Ensembl_canonical'] == 'yes'
        assert row['MANE_Select'] == 'yes'
        assert row['appris_principal_1'] == 'yes'

    def test_raf1_202_no_flags(self):
        df = parse_gtf(DUMMY_GTF)
        row = df[df['Transcript ID'] == 'ENST00000398417.3'].iloc[0]
        assert row['Ensembl_canonical'] == 'no'
        assert row['MANE_Select'] == 'no'
        assert row['appris_principal_1'] == 'no'

    def test_chromosome_extracted(self):
        df = parse_gtf(DUMMY_GTF)
        raf1 = df[df['Transcript ID'] == 'ENST00000423430.6'].iloc[0]
        assert raf1['Chromosome'] == 'chr3'
        braf = df[df['Transcript ID'] == 'ENST00000288602.11'].iloc[0]
        assert braf['Chromosome'] == 'chr7'

    def test_plain_gtf(self, tmp_path):
        """Plain (non-gzip) GTF should also work."""
        plain = tmp_path / 'test.gtf'
        with gzip.open(DUMMY_GTF, 'rt') as src, open(plain, 'w') as dst:
            dst.write(src.read())
        df = parse_gtf(str(plain))
        assert len(df) == 3


# ── Translations FASTA ────────────────────────────────────────────────────────

class TestLoadTranslations:
    def test_returns_dict(self):
        seqs = load_translations(DUMMY_TRANS)
        assert isinstance(seqs, dict)

    def test_three_entries(self):
        seqs = load_translations(DUMMY_TRANS)
        assert len(seqs) == 3

    def test_keys_are_transcript_ids(self):
        seqs = load_translations(DUMMY_TRANS)
        assert 'ENST00000423430.6' in seqs
        assert 'ENST00000398417.3' in seqs
        assert 'ENST00000288602.11' in seqs

    def test_sequence_not_empty(self):
        seqs = load_translations(DUMMY_TRANS)
        assert len(seqs['ENST00000423430.6']) > 0

    def test_no_header_char_in_sequence(self):
        seqs = load_translations(DUMMY_TRANS)
        for seq in seqs.values():
            assert '>' not in seq


# ── Isoform identification ────────────────────────────────────────────────────

class TestIsoformIdentification:
    @pytest.fixture
    def merged_df(self):
        blast = pd.read_csv(DUMMY_BESTMAPS, sep='\t')
        gtf   = parse_gtf(DUMMY_GTF)[['Transcript ID', 'Chromosome', 'gene_type',
                                       'Ensembl_canonical', 'MANE_Select',
                                       'appris_principal_1']]
        return blast.merge(gtf, on='Transcript ID', how='left').dropna()

    def test_main_isoform_column_added(self, merged_df):
        result = isoform_identification(merged_df)
        assert 'main_isoform' in result.columns

    def test_main_isoform_id_column_added(self, merged_df):
        result = isoform_identification(merged_df)
        assert 'main_isoform_id' in result.columns

    def test_one_main_per_gene(self, merged_df):
        result = isoform_identification(merged_df)
        for gene, grp in result.groupby('Gene_Gencode'):
            n_main = (grp['main_isoform'] == 'yes').sum()
            assert n_main == 1, f'{gene}: expected 1 main isoform, got {n_main}'

    def test_raf1_201_is_main(self, merged_df):
        """RAF1-201 has identical + MANE + Ensembl + appris → must be main."""
        result = isoform_identification(merged_df)
        row = result[result['Transcript ID'] == 'ENST00000423430.6'].iloc[0]
        assert row['main_isoform'] == 'yes'

    def test_raf1_202_not_main(self, merged_df):
        result = isoform_identification(merged_df)
        row = result[result['Transcript ID'] == 'ENST00000398417.3'].iloc[0]
        assert row['main_isoform'] == 'no'

    def test_single_transcript_gene_is_main(self, merged_df):
        """BRAF has only one transcript → must be main."""
        result = isoform_identification(merged_df)
        row = result[result['Gene_Gencode'] == 'ENSG00000157764.13'].iloc[0]
        assert row['main_isoform'] == 'yes'


# ── End-to-end CLI invocation ─────────────────────────────────────────────────

class TestCLI:
    def test_cli_produces_output_files(self, tmp_path):
        worker = os.path.join(BIN_DIR, 'create_sequence_table_worker.py')
        result = subprocess.run(
            [sys.executable, worker,
             '--blast_best',   DUMMY_BESTMAPS,
             '--gtf',          DUMMY_GTF,
             '--translations', DUMMY_TRANS,
             '--output_dir',   str(tmp_path),
             '--cutoff',       '80'],
            capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr

        for fname in ['loc_chrom_with_names.tsv',
                      'loc_chrom_with_names_isoforms_with_seq.tsv',
                      'loc_chrom_with_names_main_isoform.tsv']:
            assert (tmp_path / fname).exists(), f'Missing output: {fname}'

    def test_output_has_chromosome_column(self, tmp_path):
        worker = os.path.join(BIN_DIR, 'create_sequence_table_worker.py')
        subprocess.run([sys.executable, worker,
                        '--blast_best',   DUMMY_BESTMAPS,
                        '--gtf',          DUMMY_GTF,
                        '--translations', DUMMY_TRANS,
                        '--output_dir',   str(tmp_path),
                        '--cutoff',       '80'],
                       capture_output=True, text=True)
        df = pd.read_csv(tmp_path / 'loc_chrom_with_names_isoforms_with_seq.tsv', sep='\t')
        assert 'Chromosome' in df.columns
        assert 'main_isoform' in df.columns
        assert 'Sequence' in df.columns


# ---------------------------------------------------------------------------
# all_isoform_mapping: one winner per UniProt isoform
# ---------------------------------------------------------------------------
import numpy as np  # noqa: E402


def _alliso_winners(df):
    """Mirror the worker's all_isoform_mapping selection: one winner per
    Entry_Isoform, canonical accession flagged as the gene main."""
    sel = isoform_identification(df.copy(), gene_col="Entry_Isoform")
    winners = sel[sel["main_isoform"] == "yes"].copy()
    is_canon = ~winners["Entry_Isoform"].astype(str).str.contains("-", na=False)
    winners["main_isoform"] = np.where(is_canon, "yes", "no")
    return winners


def _raf1_candidates():
    rows = []

    def add(tn, ei, punct, mane="no", ens="no", appr="no", cov=99.0):
        rows.append({"Transcript name": tn, "Protein_ID": tn, "Entry_Isoform": ei,
                     "Gene_Gencode": "RAF1", "alignmentpuntcuality": punct,
                     "MANE_Select": mane, "Ensembl_canonical": ens,
                     "appris_principal_1": appr, "coverage": cov})
    add("RAF1-201", "P04049", "identical", mane="yes", ens="yes", appr="yes", cov=100)
    add("RAF1-217", "P04049", "identical", cov=100)
    add("RAF1-262", "P04049", "identical", cov=100)
    add("RAF1-244", "P04049", "aligned", cov=95)
    add("RAF1-205", "P04049-2", "identical", cov=100)
    add("RAF1-219", "P04049-2", "aligned", cov=92)
    add("RAF1-253", "P04049-2", "aligned", cov=90)
    return pd.DataFrame(rows)


class TestAllIsoformWinners:
    def test_one_winner_per_uniprot_isoform(self):
        winners = _alliso_winners(_raf1_candidates())
        assert len(winners) == 2
        assert set(winners["Entry_Isoform"]) == {"P04049", "P04049-2"}

    def test_canonical_winner_is_mane_transcript(self):
        winners = _alliso_winners(_raf1_candidates())
        canon = winners[winners["Entry_Isoform"] == "P04049"].iloc[0]
        assert canon["Transcript name"] == "RAF1-201"
        assert canon["main_isoform"] == "yes"

    def test_alt_isoform_winner_is_identical_match(self):
        winners = _alliso_winners(_raf1_candidates())
        alt = winners[winners["Entry_Isoform"] == "P04049-2"].iloc[0]
        assert alt["Transcript name"] == "RAF1-205"   # the identical match
        assert alt["main_isoform"] == "no"            # not the gene canonical
