# ELM — Eukaryotic Linear Motifs

## Description

Short linear motifs (SLiMs) are compact functional sites in proteins, typically 3–10 residues long, that mediate protein-protein interactions, post-translational modifications, and targeting signals. The ELM database curates experimentally validated instances in eukaryotes.

The pipeline filters the ELM instance table to Homo sapiens proteins and maps them to all isoforms via the `TRANSCRIPT_MAP` step.

## Data source

- **File:** `legacy_data/elm/elm_instances-2023.tsv`
- **Origin:** ELM database (http://elm.eu.org), 2023 snapshot
- **Organism filter:** Homo sapiens only
- **Update policy:** Static local file; not auto-downloaded. Replace the file to upgrade.

## Output columns

`unmapped/annotations/elm.tsv` (UniProt-keyed):

| Column | Description |
|--------|-------------|
| `Accession` | UniProt accession |
| `ELMIdentifier` | ELM class name (e.g. `MOD_CDK_SPxK_1`) |
| `ProteinName` | UniProt entry name |
| `Start` | Motif start position (1-based, UniProt coordinates) |
| `End` | Motif end position (1-based, inclusive) |
| `Instances (Matched Sequence)` | Matched amino acid sequence |
| `Logic` | `true positive` / `false positive` / `unknown` |
| `Organism` | Source organism |

After `TRANSCRIPT_MAP`:

`mapped/annotations/elm.tsv` adds `Protein_ID` and `homology_transfer` columns. Positions are re-mapped to the Gencode transcript coordinate system.

## Notes

- The 2023 snapshot contains ~3,700 human instances across ~290 ELM classes.
- Only instances with `Logic = true positive` or `unknown` are typically used downstream; check the Django model for the current filter.
- Motif start/end are in UniProt sequence space; after transcript mapping they reflect the Gencode isoform position.
- DIBS (disordered binding sites) and MFIB (molecular function in intrinsically disordered) follow the same format as ELM but come from separate local files.
