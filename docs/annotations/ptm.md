# PTM — Post-Translational Modification Sites

## Description

Post-translational modification (PTM) sites record experimentally characterized residues that are phosphorylated, ubiquitinated, acetylated, methylated, or otherwise modified. The pipeline merges two local databases: PTMdb and PhosphoSite.

## Data sources

- **PTMdb:** `legacy_data/ptm/ptmdb/` — curated PTM instances from the PTM database
- **PhosphoSite:** `legacy_data/ptm/ptmphs/` — phosphorylation, ubiquitination, acetylation sites from PhosphoSitePlus (static local snapshot)
- **Update policy:** Both are static local files. Replace the directory contents to upgrade.

## Output columns

`unmapped/annotations/ptm_merged.tsv` (UniProt-keyed):

| Column | Description |
|--------|-------------|
| `Accession` | UniProt accession |
| `Position` | Modified residue position (1-based, UniProt coordinates) |
| `Residue` | Amino acid single-letter code |
| `PTM_Type` | Modification type (e.g. `Phosphoserine`, `Ubiquitination`) |
| `Source` | `PTMdb` or `PhosphoSite` |
| `Evidence` | Experimental evidence code (where available) |

After `TRANSCRIPT_MAP`:

`mapped/annotations/ptm_merged.tsv` adds `Protein_ID` and `homology_transfer`.

## Notes

- Duplicate entries from the two sources are deduplicated on `(Accession, Position, PTM_Type)`.
- PhosphoSite covers phospho/ubiquitin/acetyl; PTMdb covers a broader set of modification types.
- PTM positions are in UniProt sequence space and converted to Gencode isoform coordinates after transcript mapping.
- The worker is `bin/create_annotation_worker.py` (shared with ELM, DIBS, MFIB, PhasePro, Pfam).
