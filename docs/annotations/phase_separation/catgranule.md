# catGRANULE

## Description

catGRANULE predicts the propensity of a protein to partition into
membraneless organelles / RNA granules through liquid-liquid phase separation
(LLPS), scoring per-residue contributions from disorder, RNA-binding, and
compositional features.

## Data source

- **Computed:** the catGRANULE algorithm is run over each isoform sequence.
- **Origin:** catGRANULE (Bolognesi et al. / Tartaglia lab).
- **Update policy:** Recomputed each run from the isoform sequences.

## Output file

`final/phase_separation/catgranule.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Position` | Residue position (1-based) |
| `catgranule_score` | Per-residue LLPS propensity contribution |
| `catgranule_total` | Whole-protein catGRANULE score (repeated per row) |

## Notes

- Complements [PLAAC](plaac.md), [PhasePro](../disorder_function/phasepro.md),
  and [FINCHES](../disorder_function/finches.md) in the phase-separation family.
- If the predictor environment/library is unavailable, the track is written
  empty (header only) and the run does not fail.
- Worker: `bin/create_catgranule_worker.py` (`CATGRANULE_MAP`,
  `modules/pathogenicity.nf`).
