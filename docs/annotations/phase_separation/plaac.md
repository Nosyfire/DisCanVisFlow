# PLAAC

## Description

PLAAC (Prion-Like Amino Acid Composition) scores each residue for prion-like
domain character using a hidden Markov model trained on known prion domains.
Prion-like domains are strongly associated with liquid-liquid phase separation
and aggregation.

## Data source

- **Computed:** the PLAAC Java tool is run over each isoform sequence.
- **Origin:** [PLAAC](http://plaac.wi.mit.edu/) (Lancaster et al.).
- **Update policy:** Recomputed each run from the isoform sequences.

## Output file

`final/phase_separation/plaac.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Position` | Residue position (1-based) |
| `plaac_score` | Per-residue prion-like propensity score |
| `in_PRD` | Whether the residue falls in a predicted prion-like domain (PRD) |

## Notes

- Complements [catGRANULE](catgranule.md),
  [PhasePro](../disorder_function/phasepro.md), and
  [FINCHES](../disorder_function/finches.md).
- If the PLAAC jar or Java is unavailable, the track is written empty (header
  only) and the run does not fail.
- Worker: `bin/create_plaac_worker.py` (`PLAAC_MAP`,
  `modules/pathogenicity.nf`).
