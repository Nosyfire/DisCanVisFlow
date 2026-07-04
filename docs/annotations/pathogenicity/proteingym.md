# ProteinGym — Deep Mutational Scanning Benchmarks

## Description

Variant-effect measurements from **ProteinGym**, a large curated collection of
deep mutational scanning (DMS) assays used to benchmark variant-effect
predictors. Each row is an experimental DMS score for a variant, mapped to the
isoform residue.

## Data source

- **Fetch:** `fetch_proteingym_worker.py` retrieves ProteinGym DMS substitution assays.
- **Origin:** [ProteinGym](https://proteingym.org/).
- **Update policy:** Fetched/cached; refresh to update.

## Output file

`final/pathogenicity/proteingym.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Protein_position` | Residue position (1-based) in the isoform |
| `protein_variant` | Variant string (e.g. `V600E`) |
| `DMS_score` | Continuous deep-mutational-scanning fitness/effect score |
| `DMS_score_bin` | Binarised effect (e.g. neutral vs deleterious) |
| `DMS_id` | ProteinGym assay identifier |
| `uniprot_id` | UniProt accession |
| `mapping_type` | `direct` / `homology_similarity` |

## Notes

- Like [MaveDB](mavedb.md), these are *experimental* measurements — the
  benchmark ground truth against which computational predictors
  ([dbNSFP](dbnsfp.md), [AlphaMissense](alphamissense.md)) are evaluated.
- `DMS_score` sign/scale is assay-specific (`DMS_id`); use `DMS_score_bin` for a
  normalised deleterious/neutral call.
- Workers: `bin/fetch_proteingym_worker.py`, `bin/create_proteingym_worker.py`.
