# FINCHES — LLPS Saturation Mutagenesis

## Description

FINCHES computes, for every possible single amino-acid substitution in a
sequence, how the mutation changes the region's **liquid-liquid phase separation
(LLPS)** tendency, expressed as a change in interaction energy (ε). It is an
in-silico saturation-mutagenesis scan of condensate-forming propensity.

## Data source

- **Predictor:** FINCHES, run on each isoform sequence.
- **Origin:** FINCHES (Ginell, Holehouse et al.).
- **Update policy:** Recomputed from sequences each run.
- **Off by default:** enable with `--skip_finches false`. Licensed CC BY-NC 4.0.

## Output file

`final/pathogenicity/finches_saturation.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Position` | Residue position (1-based) |
| `WT_AA` | Wild-type amino acid |
| `Mut_AA` | Substituted amino acid |
| `WT_Epsilon` | Interaction energy (ε) of the wild-type context |
| `Mut_Epsilon` | Interaction energy (ε) after the substitution |
| `Delta_Epsilon` | Δε = Mut − WT |

## Notes

- **Positive Δε** = the mutation *increases* LLPS tendency; **negative Δε** =
  *decreases* it.
- Complements the curated [PhasePro](phasepro.md) phase-separation
  regions with a per-mutation quantitative scan.
- Non-commercial licence (CC BY-NC 4.0); disabled by default for that reason.
- Worker: `bin/create_finches_worker.py` (Module 8h).
