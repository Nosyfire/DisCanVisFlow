# Cancer Drivers

## Description

Two complementary cancer driver gene annotations are provided:

1. **Cancer Gene Census (CGC)** — curated list of genes causally implicated in cancer (COSMIC), with tier 1 (strong evidence) and tier 2 (moderate evidence) classifications.
2. **Cosmic Compendium** — per-cohort driver gene scores from large-scale computational driver discovery across TCGA and other cancer cohorts.

Both tables are pre-processed and filtered to proteins in the current run.

## Data sources

### Cancer Gene Census

- **Parameter:** `--census_driver` pointing to a pre-processed file `census_driver.tsv`
- **Origin:** COSMIC Cancer Gene Census (https://cancer.sanger.ac.uk/census)
- **Coverage:** Tier 1 and tier 2 driver genes

### Cosmic Compendium

- **Parameter:** `--compendium_driver` pointing to a pre-processed file `compendium_driver.tsv`
- **Origin:** Cosmic Cancer Gene Compendium — per-cohort driver scores from Bailey et al. and other large-scale driver discovery studies

## Output files

`unmapped/annotations/census_driver.tsv`
`unmapped/annotations/compendium_driver.tsv`

## Output columns

### census_driver.tsv

The columns mirror the Cancer Gene Census export, filtered to `Protein_ID`:

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Gene Symbol` | HGNC gene symbol |
| `Tier` | CGC tier (1 = high confidence, 2 = moderate evidence) |
| `Hallmark` | Cancer hallmark(s) (e.g. `proliferation, evading growth suppressors`) |
| `Chr Band` | Chromosomal location |
| `Somatic` | Whether somatic mutations are reported (Yes/No) |
| `Germline` | Whether germline mutations are reported (Yes/No) |
| `Tumour Types(Somatic)` | Cancer types with somatic mutations |
| `Tumour Types(Germline)` | Cancer types with germline mutations |
| `Cancer Syndrome` | Associated cancer syndromes |
| `Molecular Genetics` | Dominant/recessive classification |
| `Role in Cancer` | Oncogene (OG), tumour suppressor (TSG), or fusion |

### compendium_driver.tsv

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Gene` | Gene symbol |
| `Cohort` | Cancer cohort (e.g. `TCGA-BRCA`) |
| `Score` | Per-cohort driver discovery score |
| `Rank` | Driver rank within cohort |

## Notes

- If either file is missing or not provided, the corresponding output is written as an empty TSV and the pipeline does not fail.
- These annotations describe the gene as a whole, not individual variants or positions.
- Worker: `bin/create_cancer_driver_worker.py`
