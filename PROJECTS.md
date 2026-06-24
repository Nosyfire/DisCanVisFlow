# Predefined Project Runs

Each named run is a **Nextflow profile**. Outputs land under `results/<project>/`.
Combine the project profile with an environment profile (`conda` or `docker`):

```bash
nextflow run main.nf -profile <project>,conda -resume
```

Machine-specific reference/data paths are defined once in
[`conf/local_refs.config`](conf/local_refs.config) and shared by every project
profile (edit that file for a different installation).

| Profile | Scope | Mapping mode | Output |
|---------|-------|--------------|--------|
| `discanvis` | Full human proteome — DisCanVis2 DB update | `all_isoform_mapping` | `results/discanvis/` |
| `vep_benchmarking` | Full proteome, benchmark data only | `all_isoform_mapping` | `results/vep_benchmarking/` |
| `test_one_protein` | Single gene (default TP53) | `all_isoform_mapping` | `results/test_one_protein/` |
| `test_subset_of_protein` | Gene list (default TP53,RAF1,BRAF,KRAS,EGFR) | `all_isoform_mapping` | `results/test_subset_of_protein/` |

All four use **`all_isoform_mapping`**: each transcript is paired 1:1 to its best
UniProt isoform (one winner per UniProt isoform). Override with
`--mapping_mode main_isoform_mapping` if you want the canonical-only behaviour.

---

## `discanvis` — DisCanVis2 update (full proteome, all isoforms)

Produces the complete per-isoform annotation set for the DisCanVis2 web server.

```bash
nextflow run main.nf -profile discanvis,conda -resume
```

- Whole human proteome (`target_gene = null`), all UniProt isoforms.
- Every module runs (annotations, disorder, PDB, genome/mutation maps, drivers,
  conservation, pathogenicity, …).
- Output: `results/discanvis/final/…`

---

## `vep_benchmarking` — variant-effect-predictor benchmark set

Human proteome restricted to the data needed to benchmark variant effect
predictors. Keeps **ClinVar variants, MaveDB, ProteinGym, Combined Disorder,
conservation, dbNSFP/pathogenicity, AlphaMissense**; skips PDB, PPI, ScanSite,
coiled-coils, PEM, cancer drivers, OMIM, DepMap, and the TCGA/cBioPortal cohorts.

```bash
nextflow run main.nf -profile vep_benchmarking,conda -resume
```

- Output: `results/vep_benchmarking/final/…`
- Key benchmark tables:
  - `final/annotations/proteingym.tsv` — ProteinGym DMS scores + `DMS_score_bin`
  - `final/annotations/mavedb.tsv` — MaveDB single-mutant functional scores
  - `final/disorder/CombinedDisorderNew.tsv` — combined disorder
  - `final/pathogenicity/pathogenicity_scores.tsv` — dbNSFP predictor scores
  - `final/pathogenicity/alphamissense.tsv` — AlphaMissense
  - `final/conservation/…` — GOPHER + phastCons
  - `final/mutations/ClinVar/…` — mapped ClinVar variants

---

## `test_one_protein` — single-gene smoke test

```bash
# default gene is TP53
nextflow run main.nf -profile test_one_protein,conda -resume

# any other gene
nextflow run main.nf -profile test_one_protein,conda --target_gene RAF1 -resume
```

- Output: `results/test_one_protein/final/…`
- Fast (~10–20 min); good for validating changes end-to-end.

---

## `test_subset_of_protein` — small multi-gene test

```bash
# default subset: TP53,RAF1,BRAF,KRAS,EGFR
nextflow run main.nf -profile test_subset_of_protein,conda -resume

# custom subset (comma-separated, no spaces)
nextflow run main.nf -profile test_subset_of_protein,conda \
    --target_gene 'TP53,RAF1,PTEN' -resume
```

- Output: `results/test_subset_of_protein/final/…`

---

## Notes

- **`-resume`** reuses cached steps; each project has its own `workDir` so caches
  don't collide between projects.
- **SLURM**: add the `slurm` profile, e.g. `-profile discanvis,slurm,conda`.
- **Override anything** on the CLI: `--mapping_mode`, `--target_gene`,
  `--skip_pdb`, `--outdir`, individual data paths, etc.
- **Large runs** (`discanvis`, `vep_benchmarking`) hit external APIs
  (AlphaFold, optionally PDBe) and read large reference files (dbNSFP ~21 GB);
  expect multi-hour runtimes. `vep_benchmarking` skips the API-heavy structure
  modules to run faster.
- **`-stub`** validates the DAG without executing workers:
  `nextflow run main.nf -profile test_one_protein,conda -stub`.
