# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**DisCanVisFlow** — a Nextflow DSL2 pipeline that maps disease variants, functional annotations, and structural features onto every curated protein isoform in the human SwissProt proteome. Powers the DisCanVis2 web server (Django) and is usable as a standalone data-generation pipeline.

The pipeline consumes UniProt SwissProt + GENCODE FASTA/GTF references and produces per-residue coordinate maps and annotation TSVs ready for database upload or downstream ML pipelines.

## Environment Setup

```bash
# Conda env name is 'discanvis' (lowercase)
conda env create -f environment.yml
conda activate discanvis

# If the env was created before bioconda UCSC tools were added to environment.yml:
conda env update -n discanvis -f environment.yml --prune

# Disorder predictions (IUPred3/AIUPred) need scipy+torch in a separate env.
# With --data discanvis_data the SETUP_DEPS Nextflow process handles this on first run.
# For --data local runs: bash bin/setup_external_programs.sh
#
# IUPred3 / ANCHOR2: academic licence — register at https://iupred2a.elte.hu/download
# then extract into External_Programs/iupred3/

# bigBedToBed (polymorphism track): if 'which bigBedToBed' is empty:
#   curl -fsSL https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed \
#     -o "$CONDA_PREFIX/bin/bigBedToBed" && chmod +x "$CONDA_PREFIX/bin/bigBedToBed"
# A missing bigBedToBed no longer crashes the run — it just skips polymorphism.

# Java for Nextflow in NON-INTERACTIVE shells (cron/ssh/nohup/CI): the conda
# activate hook that puts Java on PATH may not run, so nextflow fails at launch.
# Export the env's JVM explicitly before `nextflow run`:
#   export JAVA_CMD="$CONDA_PREFIX/bin/java" JAVA_HOME="$CONDA_PREFIX"
#   export PATH="$CONDA_PREFIX/bin:$PATH"
```

## Running the Pipeline

Config is selected through four named axes — data source, project, machine, environment:

```bash
conda activate discanvis

# Single-gene RAF1 test — fastest (4 min on 64-CPU server, ~15-30 min on laptop)
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 -resume

# Validate DAG without running anything
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 -stub

# Cellular vulnerability run (auto-downloads all references)
nextflow run main.nf --project cellular_vulnerability --data discanvis_data --machine laptop \
    --description "Q4 2026 Turbine run" -resume

# Full DisCanVis2 update on the GPU server
nextflow run main.nf --project discanvis --data local --machine hard -resume

# Full proteome on SLURM cluster
nextflow run main.nf --project discanvis --data local --machine slurm -resume

# Include only specific annotation modules (preferred over stacking --skip flags)
# Example: RAF1 with cBioPortal + ClinVar mutations + AIUPred disorder/binding + ELM
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 \
    --modules mutations,disorder --fetch_cbioportal true --skip_iupred true -resume
# Available module names: mutations, disorder, mobidb, disprot, pdb, go, polymorphism, pem,
# coiledcoils, ppi, conservation, scansite, clinvar_disease, omim, cancer_drivers,
# alphamissense, depmap, mavedb, proteingym, dbnsfp, finches, lcr, dssp, catgranule, plaac
# ELM + Pfam + DIBS/MFIB/PhasePro/PTM always run as backbone regardless of --modules

# Skip individual predictors within a module
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 \
    --skip_alphafold true --skip_iupred true -resume

# Supply local ClinVar VCF instead of auto-download
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 \
    --clinvar_vcf /path/to/clinvar.vcf.gz -resume

# TCGA/cBioPortal MAF mutation input
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 \
    --mutation_maf /path/to/tcga.maf --mutation_source TCGA -resume

# Gene list from file
nextflow run main.nf --project cellular_vulnerability --data discanvis_data --machine laptop \
    --gene_list_file config/gene_lists/cellular_vulnerability.txt -resume
```

Config axes: `--data local|discanvis_data` + `--project <name>` + `--machine laptop|hard|slurm` + `--env conda|docker`

Docker:
```bash
docker build -t discanvis-pipeline:latest .
nextflow run main.nf --project test_one_protein --data local --machine hard --env docker --target_gene RAF1 -resume
```

## Reference Data Management

```bash
# List all cached references with sizes and dates
bin/refresh_refs.sh

# Force re-download of specific sources (then -resume to fetch only those)
bin/refresh_refs.sh clinvar
bin/refresh_refs.sh clinvar mobidb go
bin/refresh_refs.sh all          # everything except hg38/dbsnp/alphafold
bin/refresh_refs.sh --force all  # truly everything

# Generate MANIFEST.tsv (what's in references/, sizes, dates)
python bin/generate_manifest.py --no_checksum

# Extract one gene from a completed full-proteome run (no recomputation)
python bin/extract_gene_from_results.py --source results/discanvis --gene RAF1 --out results/discanvis_raf1
python bin/extract_gene_from_results.py --source results/discanvis --gene RAF1,BRAF,KRAS --out results/discanvis_kinases

# Docs-only: (re)generate enriched mapping reports for an existing results dir
# WITHOUT running the pipeline. Adds data-source versions (GENCODE/UniProt release
# + dates) and an "Input scale" section (# UniProt/GENCODE entries, direct vs
# curated-isoform mappings, genome-mapped isoforms). Use after extract_gene_*,
# or to backfill provenance onto an old run. Reference FASTAs auto-detected from
# config/data/local.config + references/uniprot/; cleans stale report artifacts.
bash bin/generate_docs.sh results/discanvis
bash bin/generate_docs.sh results/discanvis_raf1 --mapping-mode all_isoform_mapping
```

## Running Tests

```bash
# All tests from project root
pytest tests/ -v

# Single test file
pytest tests/test_create_mutation_map_worker.py -v

# Single test function
pytest tests/test_create_transcript_map_worker.py::TestBoundsCheck -v
```

Tests call `bin/*.py` scripts as subprocesses with dummy input. No Nextflow required.

## IDP Dataset Requests

When someone asks for IDP (intrinsically disordered protein) data, disorder annotations, or an annotated feature set for one or more proteins, follow this workflow:

### Step 1 — Check for existing full-proteome run
```bash
ls results/discanvis/final/disorder/ | head -3
```
If `CombinedDisorderNew.tsv` appears, extraction takes seconds. Skip straight to Step 2.

### Step 2 — Extract or run

**Existing run → extract (preferred):**
```bash
# Single gene
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis --gene RAF1 --out results/idp_RAF1

# Gene list (comma-separated or from file)
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis --gene RAF1,TP53,BRAF --out results/idp_custom

conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis --gene_list_file my_genes.txt --out results/idp_custom
```

**No existing run → pipeline:**
```bash
conda activate discanvis
# Single gene (~4–10 min on server)
nextflow run main.nf --project test_one_protein --data local --machine hard \
    --target_gene RAF1 -resume
# Gene list
nextflow run main.nf --project discanvis --data local --machine hard \
    --gene_list_file my_genes.txt -resume
# Full proteome (~24 h)
nextflow run main.nf --project discanvis --data local --machine hard -resume
```

### Step 3 — Propose, confirm, run
Show the exact command, ask "Shall I run this?", wait for confirmation, then execute.

### IDP-relevant outputs (all in `results/<project>/final/`)
| Directory | Contents |
|-----------|----------|
| `disorder/` | IUPred3, ANCHOR2, AIUPred, MobiDB, DisProt, CombinedDisorder |
| `annotations/` | ELM (+ classes, switches), DIBS, MFIB, PhasePro, PTM, Pfam, UniProt ROI/binding, PEM, GO, coiled-coils, PPI, low-complexity (SEG) |
| `sequence/` | Isoform table with sequences + genomic coordinates |
| `structure/` | AlphaFold pLDDT (`AlphaFoldTable.tsv`), RSA (`rsa_scores.tsv`), DSSP (`dssp.tsv`), PDB coverage + unobserved/disordered regions (`pdb_structures.tsv`, `pdb_missing.tsv`) |
| `phase_separation/` | catGRANULE (`catgranule.tsv`), PLAAC (`plaac.tsv`) |
| `position/` | Position-based annotations |

## Agentic Behavior & Autonomous TDD Protocol

You act as a fully autonomous developer. When assigned a task:

1. **Do not ask for permission to test.** Run bash commands, Python scripts, and `pytest` freely.
2. **Follow the Autonomous Loop:**
   - Write failing `pytest` cases first.
   - Write the implementation code.
   - AUTOMATICALLY run tests: `pytest tests/test_your_module.py -v`
3. **Self-Correction:** If tests fail, DO NOT stop. Analyze, fix, re-run immediately.
4. **When to stop:** Only report back when all tests are GREEN, or after 3 consecutive failures.
5. **Silence is golden:** Minimal explanations during the loop. Give a summary when done.

## Architecture

### Core Design Principles

- **Nextflow as orchestrator**: workflow logic, caching (`storeDir`), parallelism, and profiles
- **Python workers in `bin/`**: every compute step is a standalone `bin/*.py` script with `argparse`, callable independently for testing
- **Composable profiles**: combine an environment profile (`conda`/`docker`) with a run profile (`raf1`/`full`/`discanvis_data`/`slurm`)
- **`assets/NO_FILE`**: sentinel file used as a placeholder for optional Nextflow inputs
- **`Protein_ID` as primary key**: always `Gencode transcript name` (e.g. `RAF1-201`), not UniProt Entry_Isoform

### Pipeline DAG

```
UniProt FASTA + GENCODE FASTA
        │
   SUBSET_FASTA (exact gene match: |RAF1| / GN=RAF1 )
        │
MAKEBLASTDB × 2 ──► BLASTP × 2 (reciprocal) ──► MERGE_BLAST_HITS ──► ID_MAP
                                                                          │
                              ┌───────────────────────────────────────────┤
                              │                                           │
                       SEQUENCE_PROCESS ◄── GENCODE GTF              (blast maps)
                              │
         ┌────────────────────┼────────────────────────────────────────────────┐
         │                    │                                                │
    SUBSET_CDNA          BLAT_ALIGN        ANNOTATION_MAP ◄── ELM/DIBS/MFIB/PhasePro/PTM/Pfam
         │                    │            DISORDER_MAP ◄── IUPred3/ANCHOR2/AIUPred/AlphaFold
    GENOME_MAP ◄─────────────┘            PDB_MAP
    combined_map.map                       EXON_MAP ◄── GENOME_MAP
         │                                GO_MAP ◄── FETCH_GO (goa_human.gaf + go.obo)
    MUTATION_MAP ◄── ClinVar/MAF/VCF      POLYMORPHISM_MAP ◄── dbSnp155Common.bb
         │                                PEM_MAP + PEM_TRANSFER_MAP
    CLINVAR_DISEASE_BUILD ◄── MONDO OBO   COILEDCOILS_MAP, PPI_MAP, CONSERVATION_MAP
    + mapped mutations                    DBNSFP_MAP (raw dbNSFP 5.x/chr*.gz) / PATHOGENICITY_MAP (pre-mapped)
    DEPMAP_MAP ◄── DepMap TSV             TRANSCRIPT_MAP ◄── annotation + disorder
```

### Directory Structure

```
DisCanVisFlow/
├── work/
│   ├── local/           Nextflow task cache for --data local runs
│   ├── discanvis_data/  Nextflow task cache for --data discanvis_data runs
│   └── benchmark/       ad-hoc benchmark work dirs
├── references/          storeDir cache for all FETCH_* downloads (shared across all runs)
│   └── MANIFEST.tsv     auto-generated by bin/generate_manifest.py
├── results/<project>/
│   ├── intermediate/    Entry_Isoform-keyed staging inputs to TRANSCRIPT_MAP
│   └── final/           ALL DB-ready outputs (Protein_ID-keyed)
│       ├── sequence/, genome/, mutations/, annotations/, disorder/
│       ├── structure/ (AlphaFold pLDDT + RSA + DSSP + PDB), conservation/,
│       ├── position/, disease/, pathogenicity/, drivers/, phase_separation/
│       └── mapping_reports/
└── config/
    ├── data/            local.config | discanvis_data.config
    ├── projects/        cellular_vulnerability | discanvis | vep_benchmarking | test_* | full_discanvis
    ├── machines/        laptop | hard | medium | low | slurm
    └── envs/            conda | docker
```

**Cross-project data reuse**: same `--data` flag → same `work/<data>/` dir → Nextflow `-resume` shares all task cache automatically between projects. BLAST, GENOME_MAP, and FETCH_* steps computed for `cellular_vulnerability` are fully reused when running `discanvis` with the same `--data local`.

**Single-gene extraction from full run** (no recomputation):
```bash
python bin/extract_gene_from_results.py --source results/discanvis --gene RAF1 --out results/discanvis_raf1
```

### Module → File Mapping

| Module | Nextflow file | Python worker | Key output |
|--------|--------------|---------------|------------|
| 0 — FASTA / BLAST | `modules/blast_search.nf` | `create_blast_table_worker.py` | `bestsequences.tsv` |
| 1 — ID Map | `modules/blast_mapping.nf` | `create_id_map_worker.py` | `bestmaps_blast_gene_transcript.tsv` |
| 2 — Sequence Process | `modules/sequence_process.nf` | `create_sequence_table_worker.py` | `loc_chrom_with_names_isoforms_with_seq.tsv` |
| 3 — Genome Mapping | `modules/genome_mapping.nf` | `create_genome_map_worker.py` | `combined_map.map` |
| 4 — Mutation Mapping | `modules/mutation_mapping.nf` | `create_mutation_map_worker.py` | `Missense/Frameshift/Nonsense/Indel_filter_mutations_mapped.tsv` |
| 5a — Annotation | `modules/annotation_backbone.nf` | `create_annotation_worker.py` | `elm.tsv`, `dibs.tsv`, `mfib.tsv`, `phasepro.tsv`, `ptm_merged.tsv`, `pfam_domains.tsv`, `uniprot_roi.tsv`, `uniprot_binding.tsv` |
| 5b — Disorder | `modules/disorder.nf` | `create_disorder_worker.py` | `IUPredscores.tsv`, `AIUPredscores.tsv`, `AIUPredBinding.tsv`, `CombinedDisorderNew.tsv` (all → `final/disorder/`); `AlphaFoldTable.tsv` → `final/structure/` |
| 5c — PDB | `modules/structure.nf` | `create_pdb_worker.py` | `final/structure/pdb_structures.tsv`, `final/structure/pdb_missing.tsv` |
| 5d — Exon | `modules/structure.nf` | `create_exon_worker.py` | `exon.tsv` |
| 5e — Transcript Map | `modules/annotation_backbone.nf` | `create_transcript_map_worker.py` | Protein_ID-keyed mapped copies of all annotation TSVs |
| 5f — GO Terms | `modules/functional.nf` | `create_go_worker.py` | `go_terms.tsv` |
| 5g — Polymorphism | `modules/functional.nf` | `create_polymorphism_worker.py` | `polymorphism.tsv` (rsid + ref/alt + allele freq from dbSNP 155) |
| 5h — PEM | `modules/functional.nf` | `create_pem_worker.py`, `create_pem_transfer_worker.py` | `pem_core_motifs.tsv`, `pem_core_motifs_mapped.tsv` |
| 5i — Coiled Coils | `modules/functional.nf` | `create_coiledcoils_worker.py` | `coiled_coils.tsv`, `DeepCoil.tsv` |
| 5j — PPI | `modules/functional.nf` | `create_ppi_worker.py` | `interactions.tsv` |
| 5k — ScanSite | `modules/functional.nf` | `create_scansite_worker.py` | `scansite.tsv` |
| 5m — Position-Based | `modules/disorder.nf` | `create_position_based_worker.py` | `final/position/position_based_annotations.tsv`; `rsa_scores.tsv` → `final/structure/` |
| 5n — ELM Classes | `modules/annotation_backbone.nf` | `create_elm_class_worker.py` | `elm_classes.tsv` |
| 5o — MobiDB | `modules/disorder.nf` | `create_mobidb_worker.py` | `mobidb_disorder.tsv` |
| 5p — DisProt | `modules/disorder.nf` | `create_disprot_worker.py` | `final/disorder/disprot.tsv` (curated disorder regions, IDPO/GO terms + DisProt dataset; coordinate-validated per isoform) |
| 5q — ELM Switches | `modules/annotation_backbone.nf` | `create_elm_switches_worker.py` | `final/annotations/elmswitches_mapped.tsv` (ELM molecular switches) |
| 5r — LCR | `modules/structure.nf` | `create_lcr_worker.py` | `final/annotations/low_complexity.tsv` (SEG low-complexity regions via `segmasker`) |
| 5s — DSSP | `modules/structure.nf` | `create_dssp_worker.py` | `final/structure/dssp.tsv` (8/3-state secondary structure + RSA from AlphaFold model) |
| 7 — Conservation | `modules/functional.nf` | `create_conservation_worker.py` | `conservation_multiple_level.tsv`, `conservation_phastcons.tsv` |
| 8a — ClinVar disease | `modules/disease.nf` | `create_clinvar_disease_build_worker.py` | `final/disease/clinvar_disease.tsv` |
| 8f — dbNSFP / Pathogenicity | `modules/pathogenicity.nf` | `create_dbnsfp_map_worker.py` (raw) / `create_pathogenicity_worker.py` (pre-mapped) | `final/pathogenicity/dbnsfp_scores.tsv` (raw merged 5.x `.gz` or per-chr `chr*.gz` → `DBNSFP_MAP`; ~110 cols: 37 predictor scores+rankscores, CADD, conservation, gnomAD 4.1 joint AF) or `pathogenicity_scores.tsv` (pre-mapped → `PATHOGENICITY_MAP`) |
| 8g — ProteinGym | `modules/pathogenicity.nf` | `create_proteingym_worker.py` | `proteingym.tsv` |
| 8h — FINCHES | `modules/pathogenicity.nf` | `create_finches_worker.py` | `finches_saturation.tsv` (off by default; `--skip_finches false` to enable; CC BY-NC 4.0) |
| 8i — catGRANULE | `modules/pathogenicity.nf` | `create_catgranule_worker.py` | `final/phase_separation/catgranule.tsv` (per-residue LLPS propensity) |
| 8j — PLAAC | `modules/pathogenicity.nf` | `create_plaac_worker.py` | `final/phase_separation/plaac.tsv` (prion-like domain score) |
| Report | `modules/report.nf` | `create_mapping_report_worker.py` | `mapping_reports/` (runs last) |
| Scatter | `modules/annotation_backbone.nf` (`SPLIT_SEQ_TABLE`) | `split_seq_table.py` | N gene-balanced seq-table chunks (`--scatter_chunks N`) |
| Reference fetches | `modules/fetch_references.nf` | — | UniProt/GENCODE/ClinVar/GO/MobiDB/MONDO/AlphaMissense/IntAct/BioGRID/HIPPIE (cached via `storeDir`) |

### Key Conventions

- **`Protein_ID`**: Gencode transcript name (e.g. `RAF1-201`). Primary key for all mapped outputs.
- **Exact gene subsetting**: GENCODE headers use `|RAF1|` (pipe-delimited), UniProt uses `GN=RAF1 ` (trailing space) — both required to avoid `TRAF1`, `ZTRAF1`
- **`storeDir` caching**: Reference downloads cached in `references/`; if a storeDir file is 0 bytes (failed download), delete it and re-run. `-stub` writes to `references/_stub/` and never pollutes the real cache.
- **Mutation input is mutually exclusive**: `--clinvar_vcf` OR `--mutation_maf` OR `--mutation_vcf`, not combined
- **TCGA MAF QC**: `--mutation_source TCGA` truncates barcodes to 12 chars; `--mutation_hypermutation_threshold 1500` drops hot samples; `--no_hgvsp_validation` disables ref-AA check
- **ClinVar disease build**: when `hg38_2bit` + `clinvar_disease_from_mutations=true` + `mondo_obo` set, `CLINVAR_DISEASE_BUILD` runs from `MUTATION_MAP` outputs (not a static filter table)
- **dbNSFP dual mode**: `--dbnsfp_raw_dir` → `DBNSFP_MAP` (via `combined_map.map`); `--dbnsfp_tsv` → `PATHOGENICITY_MAP` (pre-mapped Protein_ID-keyed TSV). Mutually exclusive; raw takes priority. Raw mode accepts either a **single merged dbNSFP 5.x gzip** (e.g. `dbNSFP5.3.1a_grch38.gz`, ~50 GB — detected automatically, streamed once via an inverted `(chr,pos)` index built from `combined_map.map`, `pigz`-accelerated) or a **directory of legacy per-chr `chr*.gz`** files. Kept columns are pattern-selected (`select_keep_columns`): all `_score` + `_rankscore` predictors + CADD + GERP/phyloP/phastCons + `gnomAD4.1_joint_AF`/`_POPMAX_AF`. **gnomAD allele frequency comes from dbNSFP** here — no separate gnomAD fetch needed.
- **Polymorphisms (Module 5g)**: `--dbsnp_bb` (`dbSnp155Common.bb`) → `create_polymorphism_worker.py` runs `bigBedToBed` over each isoform's genomic region from `combined_map.map`, maps every SNV to a protein residue, emits `rsid + ref/alt + allele_frequency + Type` for all isoforms containing the codon.
- **PEM isoform transfer**: `--pem_transfer true` (default) writes `final/annotations/pem_core_motifs_mapped.tsv` via sequence homology
- **Module 3+ requires `hg38_2bit`**: Genome/Mutation/Exon/Polymorphism mapping is skipped when `params.hg38_2bit` is not set
- **Isoform expansion in mutations**: `MUTATION_MAP` translates each genomic hit to all isoforms of the same gene via 3-AA context substring search; disable with `--no_isoform_expand`. Flagged in the `isoform_mapped` column — NOT homology transfer.
- **Transcript mapping (annotations)**: `TRANSCRIPT_MAP` maps UniProt-keyed annotations to all isoforms. Same UniProt accession → `mapping_type=direct`. Different isoform → region transferred only if it aligns at ≥ `--homology_min_identity` (default 0.90), flagged `mapping_type=homology_similarity`.
- **Mapping mode (`--mapping_mode`, Modules 0–1)**: `main_isoform_mapping` (default) BLASTs against canonical Swiss-Prot only. `all_isoform_mapping` adds curated isoforms from `UP000005640_9606_additional.fasta`, then produces a 1:1 transcript → UniProt isoform assignment (winner selected by MANE/APPRIS/identical/coverage ladder).
- **`combined_map.map` format**: 8 tab-separated columns: `protein_pos  aa  nuc_pos  codon  aa  gpos_csv  codon  aa`
- **PPI auto-download**: when `ppi_intact/biogrid/hippie` are all null, `FETCH_INTACT/BIOGRID/HIPPIE` download raw MiTab files and `PPI_PREPROCESS` (`create_ppi_preprocess_worker.py`) builds the processed tables (cached in `references/ppi/processed/`).
- **`--gene_list_file`**: plain-text file (one HGNC name per line, `#` comments OK) read at workflow start; overrides `--target_gene`.

### Disorder Prediction

| Tool | Library path | API | Output column |
|------|-------------|-----|---------------|
| IUPred3 | `External_Programs/iupred3` | `iupred3_lib.iupred(seq)[0]` | `IUPredscores` |
| ANCHOR2 | same | `iupred3_lib.anchor2(seq)` | `AnchorScore` |
| AIUPred disorder | `External_Programs/aiupred-caid3` | `init_models('disorder')` + `predict()` | `AIUPredscores` |
| AIUPred-Binding | `External_Programs/AIUPred` | `init_models('binding')` + `predict_binding()` | `AIUPredBinding` |
| AlphaFold pLDDT | EBI API | summary API → current pdbUrl (v6/v5/v4 fallback); strip isoform suffix | `Plldtscores` |
| Combined disorder | — | MobiDB + RSA (pLDDT) + IUPred3 + Pfam exclusion | `CombinedDisorderNew.tsv` / `_Pos.tsv` |

If direct import fails, `create_disorder_worker.py` falls back to subprocess via `params.aiupred_python`.

### Reference Data Sources

| Source | How supplied |
|--------|-------------|
| ELM instances | `legacy_data/elm/elm_instances-2023.tsv` (Homo sapiens filtered) |
| ELM classes | `legacy_data/elm/elm_classes-2025.tsv` |
| DIBS / MFIB / PhasePro | `legacy_data/dibs/`, `legacy_data/mfib/`, `legacy_data/phasepro/` |
| PTM (PTMdb + PhosphoSitePlus) | `legacy_data/ptm/ptmdb/` + `legacy_data/ptm/ptmphs/` — **not in git** (licensed; provide manually) |
| Cancer drivers | `legacy_data/drivers/` (CGC census + Compendium) |
| MobiDB | `FETCH_MOBIDB` (API, cached via `storeDir`) |
| DisProt | `FETCH_DISPROT` (disprot.org API, `term_ontology=IDPO+GO`, cached in `references/disprot/`) |
| GO | `FETCH_GO` (`goa_human.gaf.gz` + `go.obo`, cached) |
| ClinVar | `FETCH_CLINVAR` (NCBI FTP, cached in `references/clinvar/`) |
| PPI | `FETCH_INTACT/BIOGRID/HIPPIE` → `PPI_PREPROCESS` (cached in `references/ppi/`) |
| Conservation (GOPHER) | `params.gopher_conservation_table` — external pre-computed table |
| Conservation (phastCons) | `params.phastcons_dir` — chr*.bw files; requires `bigWigToBedGraph` in PATH |
| hg38.2bit | `params.hg38_2bit` or `--fetch_hg38_2bit true` |

### Pending Modules (low priority)

| Module | Data Source | Status |
|--------|-------------|--------|
| FuzDrop LLPS probability | fuzdrop.bio.unipd.it API | pending |
| Complexity tracks (DUST/TRF) | local tools | pending (SEG done — see Module 5r LCR) |

Recently completed (previously pending): ELM Switches (Module 5q), SEG
low-complexity regions (Module 5r LCR), DSSP secondary structure (Module 5s),
catGRANULE + PLAAC phase-separation predictors (Modules 8i/8j).
