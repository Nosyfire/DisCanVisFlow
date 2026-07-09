# DSSP Secondary Structure

## Description

Per-residue secondary structure and solvent accessibility assigned by **DSSP**
(`mkdssp`) from the AlphaFold model of each isoform. Provides the classic
8-state and collapsed 3-state secondary-structure calls plus relative solvent
accessibility (RSA) for every residue.

## Data source

- **Computed:** `mkdssp` is run over the AlphaFold PDB model fetched by the
  structure/disorder step.
- **Origin:** [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) (Kabsch & Sander);
  models from [AlphaFold DB](https://alphafold.ebi.ac.uk/).
- **Update policy:** Recomputed each run from the current AlphaFold model.

## Output file

`final/structure/dssp.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Position` | Residue position (1-based) |
| `aa` | Amino acid |
| `ss8` | 8-state DSSP secondary structure (`H G I E B T S -`) |
| `ss3` | Collapsed 3-state call (`H` helix / `E` strand / `C` coil) |
| `rsa` | Relative solvent accessibility (0–1) |

## Notes

- Complements the AlphaFold-derived [RSA](rsa.md) track and the
  [PDB coverage](pdb.md) track.
- If `mkdssp` is not on `PATH`, the track is written empty (header only) and the
  run does not fail.
- Worker: `bin/create_dssp_worker.py` (`DSSP_MAP`, `modules/structure.nf`).
