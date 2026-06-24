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
# Clone the repo
git clone https://github.com/<org>/discanvisflow.git
cd discanvisflow

# Create the conda environment (Python workers + Nextflow + UCSC tools)
conda env create -f environment.yml
conda activate DisCanVis
```

### 2. Set up external programs (disorder predictors)

```bash
bash bin/setup_external_programs.sh
```

Clones AIUPred from GitHub, creates the `discanvis_aiupred` conda env, and installs `bigBedToBed`. IUPred3/ANCHOR2 requires free academic registration at [iupred2a.elte.hu/download](https://iupred2a.elte.hu/download) — extract into `External_Programs/iupred3/`.

With `-profile discanvis_data,conda` the `SETUP_DEPS` Nextflow process handles this automatically on first run (no manual setup needed).

### 3. Run the pipeline

```bash
conda activate DisCanVis

# RAF1 single-gene test (~5-15 min) — validates the full DAG
nextflow run main.nf -profile raf1,conda -resume

# Full human proteome with all tracks
nextflow run main.nf -profile discanvis_data,conda \
    --skip_coiledcoils true \
    -resume
```

All reference data (UniProt, GENCODE, AlphaFold, ClinVar, dbSNP, GO, MobiDB, SIFTS, AlphaMissense …) is auto-downloaded and `storeDir`-cached in `references/` — subsequent runs reuse the cache.

> **Track selection and gene subsets:** Use `python bin/new_project.py` to interactively generate a `projects/<name>.yaml` config with custom gene scope and annotation track selection. Then pass its parameters explicitly, e.g. `--target_gene "TP53,BRCA1" --skip_pdb false --skip_conservation true`.

---

## Common run patterns

### Single gene — full annotation set (~5-15 min)

```bash
nextflow run main.nf --target_gene RAF1 -profile raf1,conda -resume
```

### Single gene, subset of tracks

```bash
nextflow run main.nf --target_gene TP53 \
    --skip_pdb true --skip_conservation true --skip_ppi true \
    -profile raf1,conda -resume
```

### Full human proteome — all tracks (~4-8 h on GPU server)

```bash
nextflow run main.nf \
    --skip_coiledcoils true \
    --scatter_chunks 20 --blat_chunks 16 \
    -profile discanvis_data,conda -resume
```

### Gene list from file

```bash
nextflow run main.nf \
    --gene_list_file projects/gene_lists/cellular_vulnerability.txt \
    -profile discanvis_data,conda -resume
```

### VEP benchmarking — mutations + pathogenicity, full proteome

```bash
nextflow run main.nf -profile vep_benchmarking,conda -resume
```

### Validate the DAG without running anything

```bash
nextflow run main.nf -profile raf1,conda -stub
```

---

## Project config reference files

`projects/*.yaml` document the parameter sets for named studies. Use them as reference or copy + adapt:

| File | Purpose |
|------|---------|
| `projects/full_discanvis.yaml` | Full proteome, all tracks (DisCanVis2 update) |
| `projects/single_protein.yaml` | Single gene or small subset |
| `projects/vep_benchmarking.yaml` | Mutations + pathogenicity only |
| `projects/cellular_vulnerability.yaml` | Full proteome for Turbine cellular vulnerability model |

Generate a new config interactively:

```bash
python bin/new_project.py
```

Key parameters (pass directly on the command line or set in `nextflow.config`):

```
--target_gene "RAF1"             # or null (all), or "RAF1,TP53,BRCA1"
--gene_list_file path/to/genes.txt  # one gene per line, overrides --target_gene
--mapping_mode all_isoform_mapping  # or main_isoform_mapping
--scatter_chunks 20              # parallel DISORDER_MAP chunks (full proteome)
--skip_pdb true / false          # per-track switches for all 20+ tracks
--outdir results/my_study
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

## Profiles

Combine an **environment profile** with a **scope profile**:

```
-profile discanvis_data,conda   # zero-machine-config, all refs auto-download
-profile raf1,conda             # local machine paths, RAF1 only
-profile full,conda             # local machine paths, full proteome
-profile discanvis_data,slurm   # HPC cluster
```

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

These tracks are off by default until you supply the data:

| Track | How to enable |
|-------|--------------|
| IUPred3 / ANCHOR2 | Register at [iupred2a.elte.hu/download](https://iupred2a.elte.hu/download), extract into `External_Programs/iupred3/` |
| OMIM | Obtain OMIM API key; set `params.omim_tsv` |
| dbNSFP | Pre-mapped TSV; set `params.dbnsfp_tsv` |
| TCGA / cBioPortal | MAF files; set `params.mutation_maf` |
| DepMap | TSV; set `params.depmap_tsv` |
| GOPHER conservation | Local table; set `params.gopher_conservation_table` |
| phastCons | Local bigWig dir; set `params.phastcons_dir` |

---

## Running tests

```bash
conda activate DisCanVis
pytest tests/ -v                                          # all tests
pytest tests/test_create_disorder_worker.py -v           # single module
pytest tests/test_create_mutation_map_worker.py -v -k missense
```

---

## Docs

| File | Contents |
|------|---------|
| `PIPELINE_DESIGN.md` | Full pipeline overview: DAG, modules, design decisions, I/O |
| `PROJECTS.md` | Predefined project profiles and their purpose |
| `CITATIONS.md` | Citations for all data sources and tools |
| `docs/pipeline_overview.md` | DAG walkthrough and module descriptions |
| `docs/isoform_mapping.md` | Transcript-to-isoform mapping logic in detail |
| `docs/conservation_calculation.md` | GOPHER + phastCons conservation pipeline |
| `docs/annotations/` | Per-annotation-track data format documentation |

---

## Name

This pipeline is published as **DisCanVisFlow** — combining DisCanVis (the web server it was built to power) with Flow (the Nextflow-based workflow). The name reflects both the origin and the architecture.

---

## Citation

If you use this pipeline, please cite the tools and databases listed in `CITATIONS.md`.
