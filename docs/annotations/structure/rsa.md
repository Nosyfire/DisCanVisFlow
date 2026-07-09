# RSA & Position-Based Annotations

## Description

Two per-residue structural tracks:

- **RSA (relative solvent accessibility)** — how exposed each residue is to
  solvent, derived from the AlphaFold model. Low RSA = buried/core; high RSA =
  surface-exposed.
- **Position-based annotations** — a single wide table that gathers every
  per-residue score (pLDDT, RSA, disorder, conservation, Pfam membership) into
  one row per residue, for convenient position-level analysis and DB upload.

## Data source

- **RSA:** computed from the AlphaFold structure / pLDDT alongside the
  [disorder](../disorder/disorder.md) track.
- **Position table:** assembled from the disorder, conservation, and Pfam tracks
  after they are computed.
- **Update policy:** Recomputed from the per-residue tracks each run.

## Output files

| File | Location | Contents |
|------|----------|----------|
| `rsa_scores.tsv` | `final/structure/` | Per-residue RSA array |
| `position_based_annotations.tsv` | `final/position/` | One row per residue with all per-position features |

## Output columns

**`rsa_scores.tsv`**

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `rsascores` | Comma-separated per-residue RSA (one float per amino acid) |

**`position_based_annotations.tsv`** — one row per `(Protein_ID, position)`:

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `position` | Residue position (1-based) |
| `plddt` | AlphaFold pLDDT confidence |
| `rsa` | Relative solvent accessibility |
| `iupred` | IUPred3 disorder score |
| `edisorder` | AIUPred disorder score |
| `combineddisorder` | Combined disorder call |
| `phastCons` | phastCons conservation |
| `conservationGlobal` … `conservationViridiplantae` | GOPHER conservation per taxonomic level |
| `pfam` | Pfam domain membership at this residue |

## Notes

- The position table is the "long" per-residue join of tracks documented
  separately under [disorder/](../disorder/disorder.md) and
  [conservation/](../conservation/conservation.md).
- Workers: `bin/create_position_based_worker.py` (Module 5m); RSA is produced by
  `bin/create_disorder_worker.py`.
