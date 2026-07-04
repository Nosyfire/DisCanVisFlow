# Protein-Protein Interactions (PPI)

## Description

Binary protein-protein interactions from three curated databases — IntAct, BioGRID, and HIPPIE — merged and deduplicated, with PubMed publication counts. The pipeline filters the merged table to proteins in the current run.

## Data sources

| Database | Parameter | Pre-processed file |
|----------|-----------|-------------------|
| IntAct | `--ppi_intact` | `Interaction_intact.tsv` |
| BioGRID | `--ppi_biogrid` | `Interaction_biogrid.tsv` |
| HIPPIE | `--ppi_hippie` | `Interaction_hippie.tsv` |

**Default paths:** `DisCanVis_Data_Process/Processed_Data/gencode_process/interactions/Interaction_{intact,biogrid,hippie}.tsv`

Each file uses the same format: `Accession A`, `Accession B`, `Publication Identifiers` (pipe-delimited PubMed IDs).

**Update policy:** Static pre-processed files. Replace to update.

## Output file

`final/annotations/interactions.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID_A` | Gencode transcript name of protein A |
| `Protein_ID_B` | UniProt accession or Gencode name of the interacting partner |
| `database` | Source database: `IntAct`, `BioGRID`, or `HIPPIE` |
| `number_of_pubmed` | Count of distinct PubMed IDs supporting the interaction |

## Notes

- All three source files are loaded and concatenated; rows with missing `Accession A` or `Accession B` columns are skipped.
- Duplicate interactions (same pair, same database) are deduplicated after PubMed counting.
- PubMed IDs are parsed from the `Publication Identifiers` field using `pubmed:(\d+)` regex; only distinct IDs are counted.
- If any of the three input files is missing or empty, it is silently skipped; at least one source must be present for a non-empty output.
- Protein_ID_A refers to the query protein in this run; Protein_ID_B is the interaction partner (may be a UniProt accession if not present in the current run's sequence table).
- Worker: `bin/create_ppi_worker.py`
