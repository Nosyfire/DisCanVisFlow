# Pathogenicity Scores (dbNSFP)

## Description

Per-variant pathogenicity predictor scores from 14 computational tools, compiled from dbNSFP (database for Nonsynonymous SNPs' Functional Predictions). The pipeline filters a pre-processed, Protein_ID-keyed variant table to the proteins in the current run.

## Data source

- **Parameter:** `--dbnsfp_tsv` pointing to a pre-processed file (`dbNSFP_custom/mapped_filtered_mutations.tsv`)
- **Origin:** dbNSFP v4+ (https://sites.google.com/site/jpopgen/dbNSFP), custom-processed to GENCODE Protein_ID coordinates
- **Update policy:** Static pre-processed file. Supply a new file to update.

## Output file

`final/annotations/pathogenicity_scores.tsv`

## Output columns

### Variant identifiers

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `chr` | Chromosome |
| `Start_Position` | Genomic start (hg38) |
| `End_Position` | Genomic end (hg38) |
| `Protein_position` | Residue position in isoform (e.g. `600`) |
| `aaref` | Reference amino acid |
| `aaalt` | Alternate amino acid |
| `aapos` | Numeric position |
| `ref` | Reference nucleotide |
| `alt` | Alternate nucleotide |
| `rs_dbSNP` | dbSNP rsID |

### Predictor scores

| Column | Tool | Description |
|--------|------|-------------|
| `AlphaMissense_score` | AlphaMissense | Structure-based missense pathogenicity (0–1) |
| `CADD_phred` | CADD | Combined Annotation Dependent Depletion, phred-scaled |
| `CADD_raw` | CADD | Raw CADD score |
| `ClinPred_score` | ClinPred | Clinical variant pathogenicity classifier |
| `ESM1b_score` | ESM-1b | Protein language model log-likelihood ratio |
| `EVE_score` | EVE | Evolutionary model of variant effect |
| `Polyphen2_HDIV_score` | PolyPhen-2 | HumDiv-trained pathogenicity score |
| `Polyphen2_HVAR_score` | PolyPhen-2 | HumVar-trained pathogenicity score |
| `PrimateAI_score` | PrimateAI | Primate deep learning pathogenicity |
| `SIFT_score` | SIFT | Sequence-based tolerance score (lower = more deleterious) |
| `VARITY_ER_LOO_score` | VARITY | Variant Annotation Retrieving Tool (ER-LOO) |
| `VARITY_ER_score` | VARITY | VARITY ER variant |
| `VARITY_R_LOO_score` | VARITY | VARITY R-LOO variant |
| `VARITY_R_score` | VARITY | VARITY R variant |
| `REVEL_score` | REVEL | Rare Exome Variant Ensemble Learner |
| `REVEL_rankscore` | REVEL | REVEL rank score |
| `gMVP_score` | gMVP | Genome-wide Missense Variant Pathogenicity |

## Notes

- Not all predictor columns are present for every variant; missing scores appear as empty cells.
- SIFT score is inverted relative to pathogenicity: lower scores (closer to 0) indicate more deleterious variants.
- CADD phred ≥ 20 is often used as a pathogenicity threshold (top 1% of deleterious variants in the genome).
- AlphaMissense also appears in `alphamissense.tsv` (Module 8d) but that file is GENCODE isoform–keyed; this file follows the dbNSFP canonical isoform mapping.
- Worker: `bin/create_pathogenicity_worker.py`
