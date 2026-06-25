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

Use the `./run` wrapper for a clean interface:

```bash
conda activate discanvis

# RAF1 single-gene test (~5-15 min) — validates the full DAG
./run --data local --target_gene RAF1

# Validate the DAG without running anything
./run --data local --target_gene RAF1 -stub

# Full human proteome, zero-config (all references auto-downloaded)
./run --data discanvis_data -resume

# Named project with description
./run --data discanvis_data --project cellular_vulnerability \
    --description "Q4 2026 Turbine feature run" -resume

# HPC cluster
./run --data discanvis_data --project discanvis --env slurm -resume
```

The wrapper resolves `--env/--data/--project` to Nextflow profiles and prints the resolved command. All other arguments (like `-resume`, `--target_gene`, `--skip_pdb`) pass through to Nextflow unchanged.

---

## Arguments

| Argument | Values | Default | Description |
|----------|--------|---------|-------------|
| `--env` | `conda`, `docker`, `slurm` | `conda` | Execution environment |
| `--data` | `local`, `discanvis_data` | — | Reference data source. `local` reads machine-specific paths from `conf/local_refs.config`; `discanvis_data` auto-downloads all references |
| `--project` | see below | — | Preset track selection + outdir |
| `--description` | any string | — | Written to `mapping_reports/mapping_summary.md` |

Any additional Nextflow argument (e.g. `-resume`, `-stub`, `--target_gene RAF1`) is passed through unchanged.

### Project presets

| `--project` | Scope | Key settings |
|-------------|-------|--------------|
| `test_one_protein` | Single gene (default: TP53) | Quick validation |
| `test_subset` | TP53, RAF1, BRAF, KRAS, EGFR | Regression testing |
| `discanvis` | Full proteome | All tracks (DisCanVis2 web server update) |
| `vep_benchmarking` | Full proteome | Mutations + pathogenicity + disorder only |
| `cellular_vulnerability` | Full proteome | ML feature set for Turbine cellular vulnerability model |

Override any project param directly: `./run --project test_one_protein --target_gene BRCA1`

---

## Common run patterns

### Single gene — full annotation set

```bash
./run --data local --target_gene RAF1 -resume
```

### Single gene, subset of tracks

```bash
./run --data local --target_gene TP53 \
    --skip_pdb true --skip_conservation true --skip_ppi true -resume
```

### Full human proteome — all tracks

```bash
./run --data discanvis_data \
    --scatter_chunks 20 --blat_chunks 16 -resume
```

### Gene list from file

```bash
./run --data discanvis_data \
    --gene_list_file projects/gene_lists/cellular_vulnerability.txt -resume
```

### Named project run

```bash
./run --data discanvis_data --project cellular_vulnerability \
    --description "Q4 2026 Turbine feature run" -resume
```

### HPC cluster

```bash
./run --data discanvis_data --project discanvis --env slurm \
    --description "Full proteome DisCanVis2 update" -resume
```

### TCGA/cBioPortal MAF mutations

```bash
./run --data local --target_gene TP53 \
    --mutation_maf /path/to/tcga.maf --mutation_source TCGA -resume
```

### Custom VCF input

```bash
./run --data local --target_gene TP53 \
    --mutation_vcf /path/to/variants.vcf.gz --mutation_source MyStudy -resume
```

### Skip individual disorder predictors

```bash
./run --data local --target_gene RAF1 \
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

## Advanced: direct Nextflow invocation

The `run` wrapper is syntactic sugar for Nextflow profiles. You can call Nextflow directly if preferred:

```bash
# Equivalent to: ./run --data local --target_gene RAF1
nextflow run main.nf -profile local,conda --target_gene RAF1

# Equivalent to: ./run --data discanvis_data --project cellular_vulnerability
nextflow run main.nf -profile discanvis_data,cellular_vulnerability,conda -resume
```

Profile load order: `data → project → env` (project settings override data defaults).

---

## Citation

If you use this pipeline, please cite the tools and databases listed in `CITATIONS.md`.
