# ClinVar Mutations

## Description

ClinVar mutations are disease-associated variants from NCBI ClinVar, mapped from genomic coordinates (hg38 VCF) to protein isoform positions using `combined_map.map`. The pipeline splits variants into four mutation types and applies isoform expansion so each genomic hit is translated to all isoforms of the same gene.

## Data source

- **Auto-download:** NCBI FTP (`ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`), cached via `storeDir` in `references/clinvar/`
- **Local override:** `--clinvar_vcf /path/to/clinvar.vcf.gz` skips the download
- **Real file location on this machine:** `/dlab/home/norbi/PycharmProjects/DisCanVis_Automated_Pipeline/data/raw/mutations/clinvar/v_2026_03_09/clinvar.vcf.gz`
- **VCF version used in testing:** 2026-03-09

## Output files

Under `results/{gene_dir}/mutations/ClinVar/`:

| File | Mutation type |
|------|--------------|
| `Missense_filter_mutations_mapped.tsv` | Missense substitutions |
| `Frameshift_filter_mutations_mapped.tsv` | Frameshift insertions/deletions |
| `Nonsense_filter_mutations_mapped.tsv` | Stop-gain (nonsense) mutations |
| `Indel_filter_mutations_mapped.tsv` | In-frame insertions/deletions |

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name (e.g. `RAF1-201`) |
| `Accession` | UniProt accession |
| `Mutation` | HGVS-style protein change (e.g. `p.Val600Glu`) |
| `Mutation Description` | Full variant description string from ClinVar |
| `Protein_position` | Residue position (1-based) in Gencode isoform |
| `Start_Position` | Genomic start position (hg38) |
| `Study Abbreviation` | `ClinVar` |
| `Study Name` | `ClinVar` |
| `Sample name` | ClinVar submission ID or blank |
| `ClinicalSignificance` | ClinVar pathogenicity classification |
| `CLNDN` | Disease name(s) from ClinVar |
| `isoform_transfer` | `True` when position was derived from another isoform via 3-AA context search |

## Notes

- Isoform expansion: each genomic variant is mapped to the primary isoform via `combined_map.map`, then propagated to all other isoforms of the same gene using a 3-amino-acid context substring search. Disable with `--no_isoform_expand`.
- Only ClinVar variants with `CLNSIG` containing "Pathogenic", "Likely pathogenic", or "Pathogenic/Likely pathogenic" are included by default (check `create_mutation_map_worker.py` for the current significance filter).
- Variants that cannot be mapped to a protein position (e.g. intronic, UTR, or no matching codon in `combined_map.map`) are silently skipped.
- Module 3 (Genome Mapping) must complete before this module can run.
- Worker: `bin/create_mutation_map_worker.py`
- The same worker handles MAF and generic VCF inputs via `--mutation_maf` / `--mutation_vcf`.
