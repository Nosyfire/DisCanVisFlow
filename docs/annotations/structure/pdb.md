# PDB Structural Coverage

## Description

Which residues of each isoform are covered by an experimental PDB structure, and
which residues within a covered range are **unobserved** (missing from the
electron density — often flexible or disordered). Coverage is derived from the
SIFTS UniProt ↔ PDB residue mapping.

## Data source

- **Bulk mode (default, `pdb_bulk = true`):** SIFTS `uniprot_segments_observed.tsv.gz`,
  downloaded once and cached in `references/sifts/`.
- **API mode:** the PDBe REST API, queried per UniProt isoform (slower — see
  [Performance](../../guide/performance.md)).
- **Origin:** [PDBe](https://www.ebi.ac.uk/pdbe/) / [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/).
- **Update policy:** Bulk file refreshed on demand; always-current in API mode.

## Output files

| File | Contents |
|------|----------|
| `final/structure/pdb_structures.tsv` | PDB entries mapped to isoform residue ranges (coverage) |
| `final/structure/pdb_missing.tsv` | Residues within a mapped range that are unobserved in the structure |

## Output columns (`pdb_structures.tsv`)

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `pdb_id` | PDB accession (e.g. `1GUA`) |
| `chain` | Chain identifier |
| `start` | First isoform residue covered |
| `end` | Last isoform residue covered |

## Notes

- `pdb_missing.tsv` flags candidate disordered/flexible regions (mapped but not
  observed), complementary to the [disorder predictors](../disorder/disorder.md).
- Prefer `--pdb_bulk true` — it replaces N per-protein API calls with one SIFTS
  download.
- Worker: `bin/create_pdb_worker.py` (Module 5c).
