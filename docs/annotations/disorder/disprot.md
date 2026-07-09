# DisProt Curated Disorder Regions

## Description

DisProt is the reference database of **manually curated** intrinsically
disordered regions (IDRs), each backed by experimental evidence and annotated
with Disorder Ontology (IDPO) and Gene Ontology (GO) function terms. Unlike the
predictor tracks, DisProt regions are literature-curated, not computed.

## Data source

- **Fetch:** `FETCH_DISPROT` downloads the current release from the DisProt API
  (`term_ontology=IDPO&term_ontology=GO`, TSV format), cached via `storeDir` in
  `references/disprot/`.
- **Origin:** [DisProt](https://disprot.org/) (Aspromonte et al.), disprot.org.
- **Update policy:** Always-current — re-fetched from the API (refresh with
  `bin/refresh_refs.sh disprot`).

## Output file

`final/disorder/disprot.tsv`

## Mapping logic

DisProt regions are keyed by UniProt accession with 1-based inclusive
coordinates on the **canonical** UniProt sequence. The worker maps each region
to every Gencode isoform of the same protein and **coordinate-validates** it:
when the region's own `Region sequence` is present, a region is emitted for an
isoform only if `sequence[start-1:end]` matches that region sequence exactly.
This avoids blindly copying a canonical-sequence region onto an isoform where
the coordinates no longer point at the same residues. Obsolete regions are
dropped.

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Entry_Isoform` | UniProt accession |
| `disprot_id` | DisProt entry ID (e.g. `DP00003`) |
| `region_id` | DisProt region ID |
| `start` | Region start (1-based, on the isoform sequence) |
| `end` | Region end (1-based, inclusive) |
| `term_namespace` | Ontology namespace (e.g. `Structural state`, `Disorder function`) |
| `term_id` | IDPO or GO term ID |
| `term_name` | Human-readable term name |
| `eco_id` | Evidence code (ECO) |
| `pmid` | Supporting publication |
| `dataset` | DisProt curation dataset (e.g. `Cancer-related proteins`); untagged regions are labelled `non-specific` |

## Notes

- Complements the predictor-based [disorder](disorder.md) tracks and the
  [MobiDB](mobidb.md) consensus with experimentally curated ground truth.
- The output may be empty for isoforms absent from DisProt.
- Worker: `bin/create_disprot_worker.py` (Module 5p). The `DISPROT_MAP` process
  runs with `--only_main_isoforms` in the pipeline.
