# Low-Complexity Regions (LCR)

## Description

Low-complexity regions are stretches of biased amino-acid composition (e.g.
homopolymers, short-period repeats) that are frequently disordered and often
overlap linear motifs and phase-separating segments. This track marks the SEG
low-complexity intervals of every isoform.

## Data source

- **Computed:** `segmasker` (NCBI BLAST+ SEG algorithm) is run over each protein
  sequence; the masked intervals become the track.
- **Origin:** [NCBI SEG / BLAST+](https://www.ncbi.nlm.nih.gov/) (Wootton &
  Federhen).
- **Update policy:** Recomputed each run from the isoform sequences.

## Output file

`final/annotations/low_complexity.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `start` | Region start (1-based, inclusive) |
| `end` | Region end (1-based, inclusive) |
| `length` | Region length in residues |

## Notes

- If `segmasker` is not on `PATH`, the track is written empty (header only) and
  the run does not fail.
- Backbone parity with the legacy pipeline's LCR track.
- Worker: `bin/create_lcr_worker.py` (`LCR_MAP`, `modules/structure.nf`).
