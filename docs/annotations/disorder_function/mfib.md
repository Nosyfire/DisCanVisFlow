# MFIB — Mutual Folding Induced by Binding

## Description

MFIB (Mutual Folding Induced by Binding) records intrinsically disordered
regions that fold **only when binding to each other** — mutual synergistic
folding, where two or more disordered chains form a stable complex that neither
adopts alone. Each entry marks the disordered region on the UniProt sequence.

## Data source

- **File:** `legacy_data/mfib/` — curated MFIB instances (Homo sapiens filtered)
- **Origin:** [MFIB](http://mfib.enzim.ttk.mta.hu/) (Fichó et al.)
- **Update policy:** Static local snapshot. Replace the directory contents to upgrade.

## Output file

`final/annotations/mfib.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Accession` | UniProt accession |
| `Entry_Isoform` | UniProt accession the region was curated against |
| `Name` | MFIB entry name / complex description |
| `Start` | Region start (1-based, UniProt coordinates) |
| `End` | Region end (1-based, inclusive) |
| `Data` | Free-text metadata from the MFIB record |
| `homology_transfer` | `True` when transferred to a different isoform by sequence match |

## Notes

- MFIB captures **mutual** folding (disorder-to-order requiring a disordered
  partner), distinct from [DIBS](dibs.md) (disorder folds onto an *ordered*
  partner) and [PhasePro](phasepro.md) (phase separation).
- Positions are curated in UniProt sequence space and converted to Gencode
  isoform coordinates by `TRANSCRIPT_MAP`.
- Worker: `bin/create_annotation_worker.py`.
