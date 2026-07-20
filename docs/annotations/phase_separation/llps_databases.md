# LLPS Databases — Regions & Protein Information

## Description

Curated liquid–liquid phase-separation (LLPS) knowledge from three community
databases, mapped onto every curated isoform:

- **PhaSepDB** — manually curated phase-separation proteins with membraneless-
  organelle (MLO) localisation, material state, and phase-separating regions.
- **LLPSDB** — proteins shown to phase-separate *in vitro*, with their
  intrinsically disordered (IDR) and low-complexity (LCR) regions.
- **DisPhaseDB** — disease-related variations in LLPS proteins (variants land in
  the [LLPS variants](../mutations/llps_variants.md) track; its regions, when a
  dump is supplied, also appear here).

Regions and protein-level attributes are **annotations**; variants are handled
separately — see [LLPS variants](../mutations/llps_variants.md).

## Data sources

| Database | Endpoint | Notes |
|----------|----------|-------|
| PhaSepDB | `db2.phasep.pro/static/db/database/phasepdbv2_1_llps.xlsx` | v2.1 LLPS sheet, auto-fetched |
| LLPSDB | `bio-comp.org.cn/llpsdbv2/download/Phase_separation_Unambiguous.zip` | v2.0, `protein.xls`, auto-fetched |
| DisPhaseDB | `disphasedb.leloir.org.ar` | API-only SPA, frequently **offline** — supply `--disphasedb_path` |

Auto-fetched and cached under `references/llps/`. A dead endpoint degrades to a
skipped source rather than failing the run.

## Output files

`final/annotations/llps_regions.tsv`

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Gene` | HGNC gene symbol |
| `start`, `end` | Region bounds (1-based), residue-validated against the isoform |
| `region_type` | `PS_region` (PhaSepDB), `IDR` / `LCR` (LLPSDB) |
| `source` | `phasepdb` / `llpsdb` / `disphasedb` |
| `mapping_type` | `direct` (same UniProt protein) |

`final/annotations/llps_proteins.tsv`

| Column | Description |
|--------|-------------|
| `Protein_ID`, `Gene` | Isoform identity |
| `organelle_mlo` | Membraneless organelle / subcellular localisation |
| `material_state` | e.g. liquid, hydrogel, solid |
| `klass` | PhaSepDB phase-separation class |
| `source`, `mapping_type` | Provenance |

## Mapping

Regions and protein info are keyed by **UniProt accession** and attached to every
isoform in the proteome sharing that accession (isoform suffix stripped). Regions
are **residue-validated**: `start..end` must fit the isoform's own sequence, so
coordinates that fall off a shorter isoform are dropped rather than mis-placed.

## Notes

- LLPSDB's construct `Mutation` column uses opaque internal codes (`M14`) with no
  substitution, so LLPSDB contributes **regions and protein info only** — no
  variants. PhaSepDB's clean HGVS mutations do become variants.
- Complements the predicted phase-separation tracks
  [catGRANULE](catgranule.md) and [PLAAC](../../annotations/phase_separation/plaac.md)
  with curated/experimental evidence.
- Workers: `bin/parse_llps_sources.py` (normalise raw downloads),
  `bin/create_llps_regions_worker.py` (map). Module: `modules/llps.nf`.
