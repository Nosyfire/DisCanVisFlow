# LLPS Variants

## Description

Amino-acid variants drawn from the LLPS databases, mapped residue-validated onto
every curated isoform:

- **DisPhaseDB** — disease-related variations in LLPS proteins (a merge of
  UniProt / ClinVar / COSMIC / DisGeNET), carried with their disease term(s).
  `variant_class = disease`.
- **PhaSepDB** — experimental construct mutations tested for effect on phase
  separation (clean HGVS, e.g. `p.V460P`). `variant_class = experimental`.

Regions and protein-level info from these databases are annotations instead —
see [LLPS databases](../phase_separation/llps_databases.md).

> **DisPhaseDB availability:** its server (`disphasedb.leloir.org.ar`) is an
> API-only SPA that is frequently offline and exposes no stable bulk-download
> URL. When it is unreachable the disease variants are simply absent; supply a
> dump with `--disphasedb_path /path/to/dump.csv` to include them. The parser
> reads any CSV/TSV with recognisable accession + variant columns.

## Output file

`final/mutations/llps_variants.tsv`

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Protein_position` | Residue position (1-based) |
| `WT_AA`, `Mut_AA` | Wild-type and substituted amino acid |
| `variant_class` | `disease` (DisPhaseDB) or `experimental` (PhaSepDB) |
| `disease_term` | Disease/phenotype term (disease variants) |
| `source_db` | `disphasedb` / `phasepdb` |
| `mapping_type` | `direct` |

## Mapping

Variants are keyed by **UniProt accession and residue identity**: a variant
(acc, pos, WT→Mut) is attached to an isoform only when that isoform's residue at
`pos` actually equals WT — the ProteinGym convention. This prevents mis-placing a
variant onto an isoform whose numbering has drifted (common where a database
annotated against a different canonical isoform). Variants that match no isoform
of their accession are dropped.

## Notes

- LLPSDB contributes no variants (its `Mutation` field is opaque internal codes).
- Complements the pipeline's other variant tracks (ClinVar, ProteinGym, MaveDB,
  dbNSFP) with an LLPS-specific, phase-separation-relevant set — useful for the
  `vep_benchmarking` project.
- Workers: `bin/parse_llps_sources.py` (normalise), `bin/create_llps_variants_worker.py`
  (map). Module: `modules/llps.nf`.
