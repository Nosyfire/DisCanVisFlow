# MaveDB — Multiplexed Assay Variant Effects

## Description

Experimental variant-effect measurements from **MaveDB**, the database of
Multiplexed Assays of Variant Effect (MAVEs) — deep mutational scanning and
similar high-throughput functional assays. Each row is a measured protein-level
effect score for a specific variant, mapped to the isoform residue.

## Data source

- **Fetch:** `fetch_mavedb_worker.py` retrieves MaveDB score sets.
- **Origin:** [MaveDB](https://www.mavedb.org/).
- **Update policy:** Fetched/cached; refresh to update.

## Output file

`final/pathogenicity/mavedb.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Protein_position` | Residue position (1-based) in the isoform |
| `prot_expr` | Protein-level variant expression string |
| `score` | Measured functional/effect score from the assay |
| `mavedb_id` | MaveDB score-set identifier |
| `urn` | MaveDB URN (stable record identifier) |
| `gene_name` | HGNC gene symbol |
| `uniprot` | UniProt accession |
| `Transcript_ID` | Source transcript identifier |
| `is_double_mutant` | `True` if the measurement is for a double mutant |
| `mapping_type` | `direct` / `homology_similarity` |

## Notes

- Scores are assay-specific — their sign and scale depend on the source score set
  (`mavedb_id` / `urn`); consult the MaveDB record for interpretation.
- Complementary to the *computational* predictors in
  [dbNSFP](dbnsfp.md), [AlphaMissense](alphamissense.md), and
  [ProteinGym](proteingym.md).
- Workers: `bin/fetch_mavedb_worker.py`, `bin/create_mavedb_worker.py`.
