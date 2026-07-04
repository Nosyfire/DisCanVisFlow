# Disease Ontology — ClinVar and OMIM

## Description

Two complementary disease ontology annotations are provided:

1. **ClinVar disease ontology** (Module 8a): variants from ClinVar are linked to Disease Ontology terms and grouped into broad disease categories (`Final_Category`) used in the DisCanVis2 paper.
2. **OMIM disease ontology** (Module 8b): OMIM disease entries linked to Disease Ontology terms via DOID cross-references.

Both tables are filtered to proteins in the current run and keyed by `Protein_ID`.

---

## ClinVar disease ontology

### Data source

- **Parameter:** `--clinvar_disease` — pre-processed file `clinvar_table_with_ontology_only_disease_pathogen_annotate.tsv`
- **Parameter:** `--clinvar_category_tsv` — paper-derived category mapping `clinvar_diseases.tsv` (optional; maps disease names to `Final_Category`)
- **Origin:** ClinVar variant submissions + Disease Ontology (DOID) cross-references

### Output file

`final/annotations/clinvar_disease.tsv`

### Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Disease` | Disease name from ClinVar |
| `DOID` | Disease Ontology identifier |
| `DO Subset` | Disease Ontology subset classification |
| `synonyms` | Comma-separated disease synonyms |
| `PhenotypeList` | ClinVar phenotype list |
| `PhenotypeIDS` | ClinVar phenotype IDs (OMIM, MeSH, etc.) |
| `xref_source` | Cross-reference database source |
| `reference_id` | Reference identifier |
| `level1`–`level13` | Disease Ontology hierarchy levels (from root) |
| `Disordered` | Number of disordered-region mutations |
| `Ordered` | Number of ordered-region mutations |
| `Total Mutations` | Total variant count |
| `Disordered Percent` | Percentage of mutations in disordered regions |
| `Final_Category` | Broad disease class (see below) |

### Final_Category values

The `Final_Category` classification groups diseases into the categories used in the DisCanVis2 paper:

| Category | Examples |
|----------|---------|
| `Cancer` | Oncogenic variants, somatic cancer syndromes |
| `Neurodegenerative` | Parkinson's, Alzheimer's, ALS |
| `Cardiovascular/Hematopoietic` | Cardiomyopathies, blood disorders |
| `Metabolic` | Diabetes, lipid disorders |
| `Developmental` | Congenital and developmental syndromes |
| `Immunological` | Autoimmune, immunodeficiency |
| `Other` | Diseases not fitting other categories |

Category assignment requires `--clinvar_category_tsv`. Without it, `Final_Category` is left empty.

---

## OMIM disease ontology

### Data source

- **Parameter:** `--omim_table` — pre-processed file `omim_table_with_ontology_with_annotate.tsv`
- **Origin:** OMIM disease entries + Disease Ontology cross-references + disorder annotation

### Output file

`final/annotations/omim_disease.tsv`

### Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Disease` | Disease name from OMIM |
| `DOID` | Disease Ontology identifier |
| `DO Subset` | Disease Ontology subset |
| `synonyms` | Disease synonyms |
| `MIMID` | OMIM MIM number |
| `level1`–`level12` | DO hierarchy levels |
| `Disordered` | Mutations in disordered regions |
| `Ordered` | Mutations in ordered regions |
| `Total Mutations` | Total variant count |
| `Disordered Percent` | Fraction in disordered regions |
| `Name` | OMIM entry name |

## Notes

- Both tables are pre-processed outside this pipeline. The pipeline only filters them to the proteins in the current run.
- For proteins with no disease associations in either database the output file contains only the header.
- ClinVar disease (Module 8a) and OMIM (Module 8b) workers: `bin/create_clinvar_disease_worker.py` and `bin/create_omim_worker.py`.
