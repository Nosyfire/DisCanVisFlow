# SNP Polymorphisms

## Description

Population-level single-nucleotide polymorphisms (SNPs) at the protein position level, derived from NCBI dbSNP. Each residue position is annotated as `Common` (minor allele frequency ≥ 1% in at least one population) or `All` (any observed SNP regardless of frequency).

This is distinct from the UniProt natural variants / polymorphism annotation (Module 5g), which is fetched from the UniProt REST API and covers curated variants. This module covers the full positional dbSNP catalogue.

## Data source

- **Parameter:** `--snp_pos_tsv` pointing to a pre-processed file `polymorphism_pos.tsv`
- **Origin:** NCBI dbSNP, processed via the `DisCanVis_Data_Process` Perl pipeline (`positional_data_process/`)
- **Format:** `AccessionPosition` column encodes `Protein_ID|position`; `Polymorphism` column is `Common` or `All`
- **Update policy:** Static pre-processed file. Replace to update.

## Output file

`unmapped/annotations/snp_polymorphisms.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name (e.g. `RAF1-201`) |
| `Position` | Residue position (1-based) |
| `Polymorphism` | `Common` (MAF ≥ 1%) or `All` (any observed SNP) |

## Notes

- The input file is large (full proteome); the worker processes it in 500,000-row chunks using pandas.
- The `AccessionPosition` field is split on `|` to extract `Protein_ID` and `Position`; rows without this format are skipped.
- If `--snp_pos_tsv` is `NO_FILE`, not provided, or the file is empty, the output is an empty TSV and the pipeline does not fail.
- Only the three columns above are kept in the output; additional columns in the source file are dropped.
- Worker: `bin/create_snp_worker.py`
