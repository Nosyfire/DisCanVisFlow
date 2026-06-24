# Disorder Predictions

## Description

Intrinsic disorder predictions are computed per isoform using four tools: IUPred3, ANCHOR2, AIUPred, and AlphaFold pLDDT. These are then combined into a consensus disorder annotation (`CombinedDisorderNew`). All predictions run on every Gencode isoform (not just the canonical UniProt entry).

## Tools and data sources

| Tool | Library / API | Output |
|------|--------------|--------|
| IUPred3 | `External_Programs/iupred3` — `iupred3_lib.iupred(seq)[0]` | Per-residue disorder score (0–1) |
| ANCHOR2 | Same library — `iupred3_lib.anchor2(seq)` | Per-residue disordered binding region score |
| AIUPred disorder | `External_Programs/aiupred-caid3` — `init_models('disorder')` + `predict()` | Per-residue disorder score |
| AIUPred binding | `External_Programs/AIUPred` — `init_models('binding')` + `predict_binding()` | Per-residue binding score |
| AlphaFold pLDDT | EBI AlphaFold summary API (v6/v5/v4 fallback) | Per-residue pLDDT confidence score |
| Combined disorder | MobiDB + RSA from pLDDT + IUPred3 + Pfam exclusion | Binary disorder annotation |

## Output files

All files are placed under `unmapped/disorder/`:

| File | Description |
|------|-------------|
| `IUPredscores.tsv` | IUPred3 per-residue scores for all isoforms |
| `Anchorscores.tsv` | ANCHOR2 per-residue binding scores |
| `AIUPredscores.tsv` | AIUPred disorder scores |
| `AIUPredBinding.tsv` | AIUPred binding region scores |
| `AlphaFoldTable.tsv` | AlphaFold pLDDT scores and model metadata |
| `CombinedDisorderNew.tsv` | Consensus disorder regions (start/end, per isoform) |
| `CombinedDisorderNew_Pos.tsv` | Per-position binary disorder flag |

## Output columns (IUPredscores.tsv example)

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name (e.g. `RAF1-201`) |
| `Position` | Residue position (1-based) |
| `Residue` | Amino acid |
| `IUPredScore` | IUPred3 disorder score (0 = ordered, 1 = disordered) |

`CombinedDisorderNew_Pos.tsv` columns: `Protein_ID`, `Position`, `Residue`, `Disordered` (0/1), contributing scores.

## Notes

- IUPred3 is called in short mode (not `"long"`) for legacy parity with `DisCanVis_Data_Process`.
- AIUPred requires scipy and PyTorch; if the main conda environment lacks these, the worker falls back to a subprocess call via `params.aiupred_python` (path to a separate conda environment Python binary).
- AlphaFold isoform suffix (e.g. `-2`) is stripped before querying the EBI API; fallback tries v6 → v5 → v4 model versions.
- Combined disorder logic: a residue is called disordered when supported by MobiDB curated/homology evidence OR when pLDDT < 70 (RSA proxy) AND IUPred3 > 0.5, excluding residues inside Pfam domains.
- Regions shorter than 5 consecutive disordered residues are filtered out.
- `CombinedDisorderNew.tsv` and `CombinedDisorderNew_Pos.tsv` are also copied to `mapped/disorder/` as pass-through (they are already `Protein_ID`-keyed).
- Worker: `bin/create_disorder_worker.py`
