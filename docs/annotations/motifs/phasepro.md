# PhasePro — Phase-Separation Drivers

## Description

PhasePro records protein regions that drive **liquid-liquid phase separation
(LLPS)** — the biomolecular condensate / membraneless organelle behaviour of
intrinsically disordered regions. Each entry marks the region responsible for
phase separation on the UniProt sequence.

## Data source

- **File:** `legacy_data/phasepro/` — curated PhasePro instances (Homo sapiens filtered)
- **Origin:** [PhasePro](https://phasepro.elte.hu/) (Mészáros et al.)
- **Update policy:** Static local snapshot. Replace the directory contents to upgrade.

## Output file

`final/annotations/phasepro.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Accession` | UniProt accession |
| `Entry_Isoform` | UniProt accession the region was curated against |
| `Name` | PhasePro entry name / condensate description |
| `Start` | Region start (1-based, UniProt coordinates) |
| `End` | Region end (1-based, inclusive) |
| `Data` | Free-text metadata from the PhasePro record |
| `homology_transfer` | `True` when transferred to a different isoform by sequence match |

## Notes

- PhasePro regions describe **phase-separation** propensity, complementary to the
  binding-coupled folding of [DIBS](dibs.md) and [MFIB](mfib.md).
- Positions are curated in UniProt sequence space and converted to Gencode
  isoform coordinates by `TRANSCRIPT_MAP`.
- Worker: `bin/create_annotation_worker.py`.
