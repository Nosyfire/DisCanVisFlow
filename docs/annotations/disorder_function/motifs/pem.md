# PEM — Predicted ELM Core Motifs

## Description

PEM (Predicted ELM Motifs) are short linear motif (SLiM) instances predicted by
scanning each isoform sequence against the ELM regular-expression library —
complementary to the experimentally annotated [ELM](elm.md) instances. The
"core motif" is the minimal residue span matched by the ELM class pattern.

## Data source

- **Patterns:** the ELM class definitions (see [ELM](elm.md)).
- **Input:** isoform sequences from the sequence table.
- **Update policy:** Predictions are recomputed each run from the sequences and the
  ELM class table.

## Output files

| File | Contents |
|------|----------|
| `final/annotations/pem_core_motifs.tsv` | Predicted core motifs on each isoform |
| `final/annotations/pem_core_motifs_mapped.tsv` | Isoform-transferred copy (when `--pem_transfer true`, the default) |

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `ELM_Accession` | ELM class accession |
| `ELMIdentifier` | ELM class name (e.g. `LIG_SH3_3`) |
| `ELMType` | Motif class (LIG / DOC / MOD / TRG / CLV / DEG) |
| `Start` | Motif start (1-based) |
| `End` | Motif end (1-based, inclusive) |
| `InstanceLogic` | Confidence flag from the pattern match |
| `References` | Supporting references from the ELM class |
| `Methods` | Detection methods from the ELM class |
| `PDB` | Associated PDB entries (where listed) |
| `Organism` | Source organism |
| `Found_Known` | `True` if the predicted motif overlaps a known ELM instance |

## Notes

- PEM is **prediction by pattern**; a hit is not experimental evidence. The
  `Found_Known` column flags overlap with a curated ELM instance.
- Isoform transfer for PEM uses sequence homology (`create_pem_transfer_worker.py`),
  distinct from the `TRANSCRIPT_MAP` path.
- Workers: `bin/create_pem_worker.py`, `bin/create_pem_transfer_worker.py`.
