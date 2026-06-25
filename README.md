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

# Create the conda environment (Python workers + Nextflow + UCSC tools)
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

Use direct Nextflow commands. Configuration is selected with separate axes:

```bash
conda activate discanvis

# Current target: cellular-vulnerability feature run on an 8 GB laptop
export NXF_OPTS='-Xms256m -Xmx1g'
nextflow run main.nf \
    --project cellular_vulnerability \
    --machine laptop \
    --description "Q4 2026 Turbine feature run" \
    -resume

# Validate the DAG without running anything
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 -stub

# Full DisCanVis update on a stronger machine
nextflow run main.nf --project full_discanvis --machine hard -resume

# Slurm
nextflow run main.nf --project full_discanvis --machine slurm -resume
```

---

## Arguments

| Argument | Values | Default | Description |
|----------|--------|---------|-------------|
| `--project` | `cellular_vulnerability`, `full_discanvis`, `discanvis`, `vep_benchmarking`, `test_one_protein`, `test_subset` | `cellular_vulnerability` | Biological/annotation track preset from `config/projects/` |
| `--machine` | `laptop`, `low`, `medium`, `hard`, `slurm` | `laptop` | Runtime/resource preset from `config/machines/` |
| `--data` | `discanvis_data`, `local` | `discanvis_data` | Reference source preset from `config/data/` |
| `--env` | `conda`, `docker`, `none` | `conda` | Software environment |
| `--ram` | Nextflow memory string, e.g. `'4 GB'` | machine default | Override the machine memory request |
| `--description` | any string | — | Written to `mapping_reports/mapping_summary.md` |
| `--target_gene` | HGNC symbol | project default | Override the target gene for test/single-gene runs |
| `--skip_*` | `true`/`false` | project default | Disable or enable individual annotation tracks |

Override any project setting directly: `nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene BRCA1 -resume`

---

## Common run patterns

### Single gene — full annotation set

```bash
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 -resume
```

### Single gene, subset of tracks

```bash
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene TP53 \
    --skip_pdb true --skip_conservation true --skip_ppi true -resume
```

### Full human proteome — all tracks

```bash
nextflow run main.nf --project full_discanvis --machine hard -resume
```

For an 8 GB laptop:

```bash
export NXF_OPTS='-Xms256m -Xmx1g'
nextflow run main.nf --project full_discanvis --machine laptop -resume
```

### Gene list from file

```bash
nextflow run main.nf --project full_discanvis --machine medium \
    --gene_list_file config/gene_lists/cellular_vulnerability.txt -resume
```

### Cellular vulnerability project

```bash
nextflow run main.nf --project cellular_vulnerability --machine laptop \
    --description "Q4 2026 Turbine feature run" -resume
```

### HPC cluster

```bash
nextflow run main.nf --project full_discanvis --machine slurm \
    --description "Full proteome DisCanVis2 update" -resume
```

### TCGA/cBioPortal MAF mutations

```bash
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene TP53 \
    --mutation_maf /path/to/tcga.maf --mutation_source TCGA -resume
```

### Custom VCF input

```bash
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene TP53 \
    --mutation_vcf /path/to/variants.vcf.gz --mutation_source MyStudy -resume
```

### Skip individual disorder predictors

```bash
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 \
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
```

Every annotation output uses `Protein_ID` (GENCODE transcript name, e.g. `RAF1-201`) as its primary key. Annotations transferred from a main isoform to an alternative isoform via sequence homology are flagged `mapping_type=homology_similarity`.

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
| DepMap | Download `OmicsSomaticMutations.csv` from [DepMap downloads](https://depmap.org/portal/download/all/) and save it as `references/depmap/OmicsSomaticMutations.csv`; the pipeline preflight confirms it and creates `references/depmap/depmap_mutations_raw.tsv` if needed. Set `--depmap_tsv` or `--depmap_raw_csv` only if using another path; use `--skip_depmap true` to omit it. |
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

## Configuration Layout

All run configuration lives under `config/`:

| Folder | Meaning |
|--------|---------|
| `config/projects/` | Biological goal and annotation-track selection |
| `config/machines/` | Runtime resources: memory, CPUs, parallelism, executor |
| `config/data/` | Reference source paths/download behavior |
| `config/envs/` | Software environment: conda, docker, or current shell |
| `config/gene_lists/` | Optional gene lists |

---

## Citation

If you use this pipeline, please cite the tools and databases listed in `CITATIONS.md`.
