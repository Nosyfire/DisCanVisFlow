# Coiled Coils (DeepCoil)

## Description

Per-residue coiled-coil probability predicted by **DeepCoil**, a deep-learning
predictor of coiled-coil structural motifs (the α-helical bundles found in
transcription factors, motor proteins, and structural proteins).

## Data source

- **Predictor:** DeepCoil, run on each isoform sequence.
- **Environment:** DeepCoil needs its own conda env (`discanvis_deepcoil`,
  TensorFlow 2.x + PyTorch); set `deepcoil_python` in `local.config`
  (see [Installation](../../guide/installation.md#4-disorder-predictors--external-programs)).
- **Update policy:** Recomputed from sequences each run.

## Output files

| File | Contents |
|------|----------|
| `final/annotations/coiled_coils.tsv` | Per-residue coiled-coil probability array |
| `DeepCoil.tsv` | Raw DeepCoil output (when produced) |

## Output columns (`coiled_coils.tsv`)

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Prob_scores` | Comma-separated per-residue coiled-coil probability (one float per amino acid, 0–1) |

## Notes

- The score array is one value per residue in sequence order — the same layout as
  the [disorder](../disorder/disorder.md) and [conservation](../conservation/conservation.md) tracks.
- Skip with `--skip_coiledcoils true` (e.g. on CUDA 12+ hardware without the
  DeepCoil env). If skipped or the env is missing, the output is empty — the run
  does not crash.
- Worker: `bin/create_coiledcoils_worker.py` (Module 5i).
