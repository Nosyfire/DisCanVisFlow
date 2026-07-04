# TCGA Mutations

## Description

Somatic mutations from a TCGA (The Cancer Genome Atlas) MAF file, mapped from
genomic coordinates to protein isoform positions using `combined_map.map`.
Variants are split into four mutation types and expanded to all isoforms of the
gene.

## Data source

- **Input:** `--mutation_maf /path/to/tcga.maf --mutation_source TCGA`
- **Origin:** TCGA MAF (Mutation Annotation Format) files.
- **Update policy:** User-supplied file.

## Quality control

`--mutation_source TCGA` enables MAF-specific QC:

- Sample barcodes are truncated to 12 characters (patient-level).
- `--mutation_hypermutation_threshold 1500` drops hypermutated samples.
- `--no_hgvsp_validation` disables the reference-amino-acid check when set.

See [Configuration § Mutation inputs](../../guide/configuration.md#mutation-inputs-clinvar--maf--vcf).

## Output files

Under `results/<project>/final/mutations/TCGA/`:

| File | Mutation type |
|------|--------------|
| `Missense_filter_mutations_mapped.tsv` | Missense substitutions |
| `Frameshift_filter_mutations_mapped.tsv` | Frameshift insertions/deletions |
| `Nonsense_filter_mutations_mapped.tsv` | Stop-gain (nonsense) mutations |
| `Indel_filter_mutations_mapped.tsv` | In-frame insertions/deletions |
| `mutation_stats.tsv` | Per-gene counts by type |

## Output columns

Same schema as [ClinVar](clinvar.md) / [cBioPortal](cbioportal.md): `Protein_ID`,
`Accession`, `Gene`, `Mutation`, `Protein_position`, study/sample provenance,
`Start_Position`, and `isoform_mapped`.

## Notes

- `--clinvar_vcf`, `--mutation_maf`, and `--mutation_vcf` are mutually exclusive —
  one mutation source per run.
- Worker: `bin/create_mutation_map_worker.py`.
