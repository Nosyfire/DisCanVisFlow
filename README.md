# DisCanVisFlow — Disease & Disorder Annotation for Human Protein Isoforms

> A Nextflow DSL2 pipeline that maps disease variants, functional annotations, and structural features onto every curated protein isoform in the human SwissProt proteome. Designed to power the DisCanVis2 web server, but fully usable as a standalone data-generation pipeline for any proteomics or structural-biology study.

---

## What it does

For each human protein (UniProt SwissProt × GENCODE):

1. **Maps each GENCODE transcript to its best UniProt isoform** via reciprocal BLASTP
2. **Produces a per-residue coordinate map** (protein position ↔ codon ↔ genomic position on hg38) via BLAT alignment
3. **Runs 20+ annotation modules** and delivers DB-ready TSVs:

| Category | Annotations |
|----------|-------------|
| Mutations | ClinVar (pathogenic/likely-pathogenic), TCGA MAF, cBioPortal MAF, custom VCF |
| Intrinsic disorder | IUPred3, ANCHOR2, AIUPred disorder, AIUPred-Binding, AlphaFold pLDDT, Combined disorder |
| SLiMs & PTMs | ELM motifs, DIBS, MFIB, PhasePro, PTMdb + PhosphoSite, Pfam domains, UniProt ROI/binding |
| Structure | PDB coverage + unobserved regions (structure-derived disorder) |
| Polymorphism | dbSNP 155 common SNPs with allele frequencies |
| Pathogenicity | dbNSFP (pre-mapped), AlphaMissense, MaveDB DMS scores, ProteinGym |
| Disease | ClinVar disease ontology (MONDO), OMIM disease + mutations |
| Interactions | IntAct, BioGRID, HIPPIE |
| Gene function | GO terms (GOA), ScanSite phospho motifs, PEM core motifs |
| Conservation | GOPHER multi-level, phastCons |
| Cancer | CGC census, Compendium, DepMap mutations |

All outputs are Protein_ID–keyed (GENCODE transcript name, e.g. `RAF1-201`) and published to `results/<project>/final/`.

---

## Quick start

### 1. Install dependencies

```bash
git clone https://github.com/Nosyfire/DisCanVisFlow.git
cd DisCanVisFlow

conda env create -f environment.yml
conda activate discanvis
```

### 2. Set up external programs (disorder predictors)

```bash
bash bin/setup_external_programs.sh
```

Clones AIUPred from GitHub, creates the `discanvis_aiupred` conda env, and installs `bigBedToBed`. IUPred3/ANCHOR2 requires free academic registration at [iupred2a.elte.hu/download](https://iupred2a.elte.hu/download) — extract into `External_Programs/iupred3/`.

With `--data discanvis_data` the `SETUP_DEPS` Nextflow process handles this automatically on first run.

### 3. Run the pipeline

```bash
conda activate discanvis

# Cellular vulnerability project — auto-downloads all references
nextflow run main.nf --project cellular_vulnerability --data discanvis_data --machine laptop \
    --description "Q4 2026 Turbine feature run" -resume

# Single-gene RAF1 test (uses pre-existing local paths, fastest)
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 -resume

# Validate the DAG without running anything
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 -stub
```

---

## Arguments

| Argument | Values | Default | Description |
|----------|--------|---------|-------------|
| `--project` | `cellular_vulnerability`, `full_discanvis`, `discanvis`, `vep_benchmarking`, `test_one_protein`, `test_subset` | `cellular_vulnerability` | Biological/annotation preset from `config/projects/` |
| `--data` | `local`, `discanvis_data` | `discanvis_data` | Reference source: `local` = pre-existing paths on this machine; `discanvis_data` = auto-download everything |
| `--machine` | `laptop`, `low`, `medium`, `hard`, `slurm` | `laptop` | Runtime/resource preset from `config/machines/` |
| `--env` | `conda`, `docker` | `conda` | Software environment |
| `--ram` | Nextflow memory string, e.g. `'4 GB'` | machine default | Override memory request |
| `--description` | any string | — | Written to `mapping_reports/mapping_summary.md` |
| `--target_gene` | HGNC symbol | project default | Override target gene for test/single-gene runs |
| `--skip_*` | `true`/`false` | project default | Disable or enable individual annotation tracks |

---

## Common run patterns

### Single gene — full annotation set

```bash
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 -resume
```

### Zero-config full run (downloads all references automatically)

```bash
nextflow run main.nf --project cellular_vulnerability --data discanvis_data --machine laptop \
    --description "Q4 2026 Turbine run" -resume
```

### Full human proteome — all tracks

```bash
nextflow run main.nf --project full_discanvis --data local --machine hard -resume
```

### HPC cluster

```bash
nextflow run main.nf --project full_discanvis --data local --machine slurm \
    --description "Full proteome DisCanVis2 update" -resume
```

### TCGA/cBioPortal MAF mutations

```bash
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene TP53 \
    --mutation_maf /path/to/tcga.maf --mutation_source TCGA -resume
```

### Custom VCF

```bash
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene TP53 \
    --mutation_vcf /path/to/variants.vcf.gz --mutation_source MyStudy -resume
```

### Skip individual annotation tracks

```bash
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 \
    --skip_alphafold true --skip_iupred true --skip_aiupred true -resume
```

---

## Output structure

```
results/<project>/
├── final/
│   ├── annotations/     ELM, DIBS, MFIB, PhasePro, PTM, Pfam, GO, polymorphism …
│   ├── disorder/        IUPredscores, AIUPredBinding, AlphaFoldTable, CombinedDisorder …
│   ├── genome/          combined_map.map, exon.tsv, genome_protein_index.tsv
│   ├── mutations/       ClinVar/, TCGA/, CBioportal/, DepMap/
│   ├── pathogenicity/   pathogenicity_scores.tsv, alphamissense.tsv, mavedb.tsv
│   ├── pdb/             pdb_structures.tsv, pdb_missing.tsv
│   ├── disease/         clinvar_disease.tsv, omim_disease.tsv
│   ├── drivers/         cancer_driver.tsv, census_driver.tsv
│   ├── position/        position_based_annotations.tsv, rsa_scores.tsv
│   └── sequence/        loc_chrom_with_names_isoforms_with_seq.tsv …
├── intermediate/        Staging TSVs (Entry_Isoform–keyed, before transcript mapping)
└── mapping_reports/
    ├── mapping_summary.md        Run metadata + provenance + run-wide coverage table
    └── mapping_coverage.tsv      Per-(Gene, annotation) coverage (full proteome runs)

work/
├── local/           Nextflow task cache for --data local runs
└── discanvis_data/  Nextflow task cache for --data discanvis_data runs

references/          Auto-downloaded reference data (storeDir cache, shared across all runs)
├── MANIFEST.tsv     Auto-generated: every cached file with size + modification date
└── ...
```

Every annotation output uses `Protein_ID` (GENCODE transcript name, e.g. `RAF1-201`) as its primary key. Annotations transferred from a main isoform to an alternative isoform via sequence homology are flagged `mapping_type=homology_similarity`.

---

## Cross-project data reuse

### Reuse pipeline cache across projects

Two runs with the **same `--data` flag share the same `work/` directory** (`work/local/` or `work/discanvis_data/`). Nextflow's `-resume` automatically reuses any task whose inputs are identical — so running `cellular_vulnerability` then `discanvis` (both `--data local`) reuses BLAST, genome mapping, and all FETCH steps without re-running them.

### Extract one gene from a completed full-proteome run

If you've already run the full proteome, extracting a single gene's data takes seconds:

```bash
python bin/extract_gene_from_results.py \
    --source results/discanvis \
    --gene   RAF1 \
    --out    results/discanvis_raf1

# Multiple genes
python bin/extract_gene_from_results.py \
    --source results/discanvis \
    --gene   RAF1,BRAF,KRAS \
    --out    results/discanvis_kinases
```

This filters all TSVs in `results/discanvis/final/` by Protein_ID prefix and writes filtered copies. No Nextflow, no recomputation.

---

## Reference data management

### Check what's cached

```bash
bin/refresh_refs.sh          # list all cached sources, sizes, and dates
python bin/generate_manifest.py --no_checksum   # write references/MANIFEST.tsv
```

### Force re-download of specific sources

```bash
bin/refresh_refs.sh clinvar              # delete cached ClinVar → re-download on next run
bin/refresh_refs.sh clinvar mobidb go   # multiple sources at once
bin/refresh_refs.sh all                  # refresh everything except hg38/dbsnp/alphafold
bin/refresh_refs.sh --force all         # refresh truly everything
```

Then re-run with `-resume` — only the deleted files are re-fetched.

### Two data modes

| Mode | When to use | Re-downloads? |
|------|-------------|--------------|
| `--data local` | Pre-existing paths on this machine; reproducible frozen snapshot | Never |
| `--data discanvis_data` | Zero-config; downloads any missing reference automatically | Only if missing (or deleted via `refresh_refs.sh`) |

---

## Environment setup details

### Disorder predictors

| Tool | Env | Notes |
|------|-----|-------|
| IUPred3 / ANCHOR2 | `discanvis` (main) | Academic licence — register at iupred2a.elte.hu; extract into `External_Programs/iupred3/` |
| AIUPred disorder | `discanvis_aiupred` | Auto-cloned from GitHub by `setup_external_programs.sh`; needs PyTorch ≥ 2.7.0+cu128 for NVIDIA Blackwell (SM 12.0) |
| AIUPred-Binding | `discanvis_aiupred` | Same env as above |
| DeepCoil | `discanvis_deepcoil` | TF 2.9 / CUDA 11 only — incompatible with CUDA 12+ Blackwell GPUs; use `--skip_coiledcoils true` on modern hardware |

### Licence-gated data (not redistributed)

| Track | How to enable |
|-------|--------------|
| IUPred3 / ANCHOR2 | Register at [iupred2a.elte.hu/download](https://iupred2a.elte.hu/download), extract into `External_Programs/iupred3/` |
| OMIM | Obtain OMIM API key; set `--omim_tsv` |
| dbNSFP | Pre-mapped TSV; set `--dbnsfp_tsv` |
| TCGA / cBioPortal | MAF files; set `--mutation_maf` |
| DepMap | Download `OmicsSomaticMutationsProfile.csv` from [DepMap downloads](https://depmap.org/portal/download/all/) and save as `references/depmap/OmicsSomaticMutationsProfile.csv`; or let auto-download try via `--data discanvis_data`. |
| GOPHER conservation | Local table; set `--gopher_conservation_table` |
| phastCons | Local bigWig dir; set `--phastcons_dir` |

---

## Running tests

```bash
conda activate discanvis
pytest tests/ -v                                          # all tests
pytest tests/test_create_disorder_worker.py -v           # single module
pytest tests/test_create_mutation_map_worker.py -v -k missense
```

---

## Configuration layout

All run configuration lives under `config/`:

| Folder | Meaning |
|--------|---------|
| `config/projects/` | Biological goal and annotation-track selection |
| `config/machines/` | Runtime resources: memory, CPUs, parallelism, executor |
| `config/data/` | Reference source paths/download behavior (`local` or `discanvis_data`) |
| `config/envs/` | Software environment: conda or docker |
| `config/gene_lists/` | Optional gene lists |

---

## Performance notes

Full-proteome benchmark on `gpu0.dlab.elte.hu` (64 CPUs): see `docs/performance_benchmark.md`.

Key findings:
- **Total wall time for single gene (RAF1)**: ~4 minutes (64 CPUs)
- **Top bottleneck**: `POLYMORPHISM_MAP` (chromosome-sweep optimized; was per-isoform bigBedToBed)
- **Second bottleneck**: `DBNSFP_MAP` — use `--dbnsfp_tsv` (pre-mapped) instead of `--dbnsfp_raw_dir` for full-proteome runs
- **PDB**: always use `--pdb_bulk true` (SIFTS join) — set by default in both `local.config` and `discanvis_data.config`

---

## Citation

If you use this pipeline, please cite the tools and databases listed in `CITATIONS.md`.
