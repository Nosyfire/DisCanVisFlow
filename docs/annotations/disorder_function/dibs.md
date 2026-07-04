# DIBS — Disordered Binding Sites

## Description

DIBS (Disordered Binding Site database) records intrinsically disordered protein
regions that fold upon binding to an ordered partner — one-to-one
disorder-to-order transitions that form well-defined complexes. Each entry marks
the disordered region on the UniProt sequence.

## Data source

- **File:** `legacy_data/dibs/` — curated DIBS instances (Homo sapiens filtered)
- **Origin:** [DIBS](http://dibs.enzim.ttk.mta.hu/) (Schad et al.)
- **Update policy:** Static local snapshot. Replace the directory contents to upgrade.

## Output file

`final/annotations/dibs.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `mapping_type` | `direct` / `homology_similarity` (see [Isoform mapping](../../pipeline/isoform_mapping.md)) |
| `homology_transfer` | `True` when transferred to a different isoform by sequence match |
| `homology_identity` | Sequence identity of the transfer (homology rows only) |
| `Entry_Isoform` | UniProt accession the region was curated against |
| `Accession` | UniProt accession |
| `Name` | DIBS entry name / partner description |
| `Start` | Region start (1-based, UniProt coordinates) |
| `End` | Region end (1-based, inclusive) |
| `Data` | Free-text metadata from the DIBS record |

## Notes

- DIBS regions describe **disorder-to-order** transitions (coupled folding and
  binding), complementary to [MFIB](mfib.md) (fuzzy/mutual folding) and
  [PhasePro](phasepro.md) (phase separation).
- Positions are curated in UniProt sequence space and converted to Gencode
  isoform coordinates by `TRANSCRIPT_MAP`.
- Worker: `bin/create_annotation_worker.py` (shared with ELM, MFIB, PhasePro,
  PTM, Pfam).
