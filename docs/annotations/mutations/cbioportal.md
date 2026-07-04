# cBioPortal Mutations

## Description

Somatic cancer mutations from [cBioPortal](https://www.cbioportal.org/), mapped
from genomic coordinates to protein isoform positions using `combined_map.map`.
Like ClinVar, variants are split into four mutation types and each genomic hit is
expanded to all isoforms of the same gene.

## Data source

- **API mode (single gene):** `--fetch_cbioportal true` queries the public
  cBioPortal REST API — no study ID needed. Default in `test_one_protein`.
- **Bulk mode (full proteome):** `--fetch_cbioportal true --cbioportal_study <id>`
  downloads one cohort MAF (e.g. `tcga_pan_can_atlas_2018`).
- **Origin:** cBioPortal aggregated public studies.
- **Update policy:** Fetched at run time; API mode is always-current.

See [Configuration § `--fetch_cbioportal`](../../guide/configuration.md#--fetch_cbioportal--somatic-mutations-from-cbioportal).

## Output files

Under `results/<project>/final/mutations/CBioportal/`:

| File | Mutation type |
|------|--------------|
| `Missense_filter_mutations_mapped.tsv` | Missense substitutions |
| `Frameshift_filter_mutations_mapped.tsv` | Frameshift insertions/deletions |
| `Nonsense_filter_mutations_mapped.tsv` | Stop-gain (nonsense) mutations |
| `Indel_filter_mutations_mapped.tsv` | In-frame insertions/deletions |
| `mutation_stats.tsv` | Per-gene counts by type |

## Output columns

Same schema as [ClinVar](clinvar.md), with cancer-study provenance:

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Accession` | UniProt accession |
| `Gene` | HGNC gene symbol |
| `Mutation` | Protein change (e.g. `p.Val600Glu`) |
| `Protein_position` | Residue position in isoform |
| `Study Abbrevation` / `Study Name` | Source cBioPortal study |
| `Sample name` | Tumour sample identifier |
| `Start_Position` | Genomic start (hg38) |
| `isoform_mapped` | `True` when derived from another isoform via 3-AA context search |

## Notes

- cBioPortal is opt-in for full-proteome runs because per-gene API calls are slow
  (~1 min/gene); use `--cbioportal_study` for bulk cohorts.
- Isoform expansion is flagged in `isoform_mapped` (not homology transfer).
- Worker: `bin/create_mutation_map_worker.py`.
