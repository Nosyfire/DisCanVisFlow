# AlphaMissense

## Description

AlphaMissense (Google DeepMind, 2023) provides per-variant pathogenicity scores for all possible single amino acid substitutions in the human proteome, including all GENCODE isoforms. Unlike the canonical AlphaMissense release (which covers only canonical UniProt sequences), the pre-processed file used here is keyed by Gencode `Protein_ID` and covers all isoforms.

## Data source

- **Parameter:** `--alphamissense_tsv` pointing to a pre-processed file (`processed_alphamissense_results_hg38_mapping.tsv`)
- **Origin:** AlphaMissense proteome-wide predictions (Cheng et al., Science 2023), re-mapped to GENCODE isoforms
- **Genome build:** hg38 coordinates
- **Update policy:** Static pre-processed file keyed by Protein_ID. Supply a new file to update.

## Output file

`final/annotations/alphamissense.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name (e.g. `RAF1-201`) |
| `uniprot_id` | Corresponding UniProt accession |
| `protein_variant` | Amino acid substitution (e.g. `V600E`) |
| `pos` | Residue position (1-based) in the Gencode isoform |
| `am_pathogenicity` | AlphaMissense score (0–1; higher = more likely pathogenic) |
| `am_class` | Classification: `likely_pathogenic`, `ambiguous`, or `likely_benign` |

## Notes

- The pre-processed file is large (full proteome); the worker reads it in 500,000-row chunks and filters by `Protein_ID`.
- Each row represents one possible amino acid substitution at one position. Positions with no predicted substitution (e.g. already the wild-type) are absent.
- The `am_class` thresholds follow the original AlphaMissense paper: score < 0.34 → `likely_benign`; 0.34–0.564 → `ambiguous`; > 0.564 → `likely_pathogenic`.
- AlphaMissense is also included as a score column in `pathogenicity_scores.tsv` (Module 8f/dbNSFP), but that file uses the canonical isoform. This file covers all Gencode isoforms.
- Worker: `bin/create_alphamissense_worker.py`
