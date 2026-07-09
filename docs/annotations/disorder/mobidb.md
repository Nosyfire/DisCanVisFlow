# MobiDB Disorder

## Description

MobiDB is a consensus database of intrinsic protein disorder that aggregates
experimental evidence (e.g. missing residues in X-ray structures, NMR
mobility), curated annotations, and predictor consensus. This track summarises
MobiDB disorder features per isoform.

## Data source

- **Fetch:** `FETCH_MOBIDB` downloads the bulk MobiDB table from the MobiDB API,
  cached via `storeDir` in `references/`.
- **Origin:** [MobiDB](https://mobidb.org/) (Piovesan et al.), mobidb.bio.unipd.it.
- **Update policy:** Always-current — re-fetched from the API (refresh with
  `bin/refresh_refs.sh mobidb`).

## Output file

`final/disorder/mobidb_disorder.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Entry_Isoform` | UniProt accession |
| `feature` | MobiDB disorder feature type (`curated-disorder-merge` or `homology-disorder-merge`) |
| `start_end` | Region span(s) as `start-end` (comma-separated for multiple regions) |
| `content_fraction` | Fraction of the sequence covered by the feature (0–1) |
| `content_count` | Number of residues covered |
| `length` | Sequence length used for the fraction |

## Notes

- MobiDB is one of the inputs to the [Combined disorder](disorder.md) track
  (MobiDB + AlphaFold pLDDT + IUPred3, with Pfam-domain exclusion).
- The output may be empty for isoforms absent from MobiDB.
- Worker: `bin/create_mobidb_worker.py` (Module 5o).
