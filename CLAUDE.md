# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**DisCanVisFlow** — a Nextflow DSL2 pipeline that maps disease variants, functional annotations, and structural features onto every curated protein isoform in the human SwissProt proteome. Powers the DisCanVis2 web server (Django) and is usable as a standalone data-generation pipeline.

The pipeline consumes UniProt SwissProt + GENCODE FASTA/GTF references and produces per-residue coordinate maps and annotation TSVs ready for database upload or downstream ML pipelines.

## Environment Setup

```bash
# Conda env name is 'DisCanVis' (capital D and V)
conda env create -f environment.yml

# Activate before running pipeline or tests
conda activate DisCanVis

# If the env was created before bioconda UCSC tools were added to environment.yml,
# update it (installs blat, pslCDnaFilter, twoBitToFa, bigWigToBedGraph, nextflow):
conda env update -n DisCanVis -f environment.yml --prune

# Disorder predictions (IUPred3/AIUPred) need scipy+torch in a separate env.
# With -profile discanvis_data,conda the SETUP_DEPS Nextflow process runs once on
# first launch: clones AIUPred from GitHub, creates the discanvis_aiupred conda env,
# installs bigBedToBed, and auto-detects the aiupred_python path. No manual config
# needed. For raf1/full profiles, run `bash bin/setup_external_programs.sh` once.
#
# IUPred3 / ANCHOR2: academic licence (Dosztányi lab, ELTE) — NO redistribution.
# Cannot be bundled in the repo or on GitHub. Each user must register at:
#   https://iupred2a.elte.hu/download
# then extract into External_Programs/, OR set params.iupred3_url to a private URL
# and SETUP_DEPS will auto-download it. Without IUPred3 those tracks will be empty.

# bigBedToBed (polymorphism track): environment.yml lists ucsc-bigbedtobed, but the
# conda solver can be very slow. If `which bigBedToBed` is empty, drop the UCSC
# static binary straight into the env (instant):
#   curl -fsSL https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed \
#     -o "$CONDA_PREFIX/bin/bigBedToBed" && chmod +x "$CONDA_PREFIX/bin/bigBedToBed"
# A missing bigBedToBed no longer crashes the run — it just skips polymorphism.
```

## Running the Pipeline

```bash
conda activate DisCanVis

# RAF1 single-gene test — recommended first run (~5-15 min)
nextflow run main.nf -profile raf1,conda -resume

# Validate DAG without running anything
nextflow run main.nf -profile raf1,conda -stub

# Full human proteome (hours) — local machine paths
nextflow run main.nf -profile full,conda -resume

# Zero-config full run (all references auto-download)
nextflow run main.nf -profile discanvis_data,conda -resume
nextflow run main.nf -profile discanvis_data,conda --target_gene 'TP53,BRCA1' -resume

# Transcript→UniProt mapping mode (default: main_isoform_mapping = canonical only).
# all_isoform_mapping pairs each transcript to its best curated SwissProt isoform
# (needs --uniprot_isoform_fasta, set in raf1/full/discanvis_data profiles).
nextflow run main.nf -profile raf1,conda --mapping_mode all_isoform_mapping -resume

# Predefined project profiles
nextflow run main.nf -profile test_one_protein,conda -resume            # default TP53
nextflow run main.nf -profile test_one_protein,conda --target_gene RAF1 -resume
nextflow run main.nf -profile test_subset_of_protein,conda -resume      # TP53,RAF1,BRAF,KRAS,EGFR
nextflow run main.nf -profile vep_benchmarking,conda -resume            # mutations + pathogenicity, full proteome
nextflow run main.nf -profile discanvis,conda -resume                   # DisCanVis2 update (full proteome)

# Gene list from file
nextflow run main.nf --gene_list_file projects/gene_lists/cellular_vulnerability.txt \
    -profile discanvis_data,conda -resume

# Skip individual disorder predictors (for faster testing)
nextflow run main.nf -profile raf1,conda --skip_alphafold true --skip_iupred true --skip_aiupred true -resume

# Skip Pfam API lookup
nextflow run main.nf -profile raf1,conda --skip_pfam_api true -resume

# Supply local ClinVar VCF instead of auto-download
nextflow run main.nf -profile raf1,conda --clinvar_vcf /path/to/clinvar.vcf.gz -resume

# TCGA/cBioPortal MAF mutation input
nextflow run main.nf -profile raf1,conda \
    --mutation_maf /path/to/tcga.maf \
    --mutation_source TCGA -resume

# Generic VCF (non-ClinVar)
nextflow run main.nf -profile raf1,conda \
    --mutation_vcf /path/to/variants.vcf.gz \
    --mutation_source CustomStudy -resume
```

Profile combinations: `raf1,conda` | `full,conda` | `raf1,docker` | `full,docker` | `discanvis_data,conda` | `discanvis_data,slurm`

Docker image must be built once before using docker profiles:
```bash
docker build -t discanvis-pipeline:latest .
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
    + mapped mutations                    DBNSFP_MAP (raw chr*.gz) / PATHOGENICITY_MAP (pre-mapped)
    DEPMAP_MAP ◄── DepMap TSV             TRANSCRIPT_MAP ◄── annotation + disorder
```

### Output Directory Structure

```
results/<project>/
├── intermediate/             Entry_Isoform-keyed staging inputs to TRANSCRIPT_MAP
│   ├── annotations/          elm.tsv, dibs.tsv, mfib.tsv, phasepro.tsv,
│   │                         uniprot_roi.tsv, uniprot_binding.tsv, ptm_merged.tsv, pfam_domains.tsv
│   └── disorder/             mobidb_disorder.tsv
└── final/                    ALL DB-ready outputs (Protein_ID-keyed)
    ├── sequence/             isoform table with sequences, coordinates, MANE/APPRIS flags
    ├── genome/               combined_map.map, exon.tsv, genome_protein_index.tsv,
    │                         genome_protein_mutations.tsv (every possible SNV reference table)
    ├── mutations/            ClinVar/, TCGA/, CBioportal/, DepMap/ — per-source TSVs
    ├── annotations/          elm, dibs, mfib, phasepro, uniprot_roi, uniprot_binding,
    │                         ptm_merged, pfam_domains, go_terms, polymorphism,
    │                         pem_core_motifs, pem_core_motifs_mapped, coiled_coils, DeepCoil,
    │                         interactions, scansite, elm_classes, homology_similarity_manifest
    ├── disorder/             IUPredscores, Anchorscores, AIUPredscores, AIUPredBinding,
    │                         AlphaFoldTable, CombinedDisorderNew, CombinedDisorderNew_Pos, rsa_scores
    ├── pdb/                  pdb_structures.tsv (chain + region + resolution), pdb_missing.tsv
    ├── conservation/         conservation_multiple_level.tsv, conservation_phastcons.tsv
    ├── position/             position_based_annotations.tsv
    ├── disease/              clinvar_disease.tsv, clinvar_disease_mutations.tsv,
    │                         omim_disease.tsv, omim_mutations.tsv
    ├── pathogenicity/        pathogenicity_scores.tsv (dbNSFP), alphamissense.tsv,
    │                         mavedb.tsv, proteingym.tsv
    └── drivers/              cancer_driver.tsv, census_driver.tsv, compendium_driver.tsv
mapping_reports/
├── <GENE>_mapping_report.md  per-gene: isoform coverage, annotation source counts, mapping QC
├── mapping_summary.md        run command, tool versions, provenance, run-wide coverage
└── mapping_coverage.tsv      flat per-(Gene × annotation) table (full proteome runs)
```

Annotations transferred from a main isoform to an alternative isoform are flagged `mapping_type=homology_similarity` and collected in `final/annotations/homology_similarity_manifest.tsv`.

### Module → File Mapping

| Module | Nextflow file | Python worker | Key output |
|--------|--------------|---------------|------------|
| 0 — FASTA / BLAST | `modules/blast_search.nf` | `create_blast_table_worker.py` | `bestsequences.tsv` |
| 1 — ID Map | `modules/blast_mapping.nf` | `create_id_map_worker.py` | `bestmaps_blast_gene_transcript.tsv` |
| 2 — Sequence Process | `modules/sequence_process.nf` | `create_sequence_table_worker.py` | `loc_chrom_with_names_isoforms_with_seq.tsv` |
| 3 — Genome Mapping | `modules/genome_mapping.nf` | `create_genome_map_worker.py` | `combined_map.map` |
| 4 — Mutation Mapping | `modules/mutation_mapping.nf` | `create_mutation_map_worker.py` | `Missense/Frameshift/Nonsense/Indel_filter_mutations_mapped.tsv` |
| 5a — Annotation | `modules/annotation_mapping.nf` | `create_annotation_worker.py` | `elm.tsv`, `dibs.tsv`, `mfib.tsv`, `phasepro.tsv`, `ptm_merged.tsv`, `pfam_domains.tsv` |
| 5b — Disorder | `modules/annotation_mapping.nf` | `create_disorder_worker.py` | `IUPredscores.tsv`, `AIUPredscores.tsv`, `AIUPredBinding.tsv`, `AlphaFoldTable.tsv`, `CombinedDisorderNew.tsv` |
| 5c — PDB | `modules/annotation_mapping.nf` | `create_pdb_worker.py` | `pdb_structures.tsv`, `pdb_missing.tsv` |
| 5d — Exon | `modules/annotation_mapping.nf` | `create_exon_worker.py` | `exon.tsv` |
| 5e — Transcript Map | `modules/annotation_mapping.nf` | `create_transcript_map_worker.py` | Protein_ID-keyed mapped copies of all annotation TSVs |
| 5f — GO Terms | `modules/annotation_mapping.nf` | `create_go_worker.py` | `go_terms.tsv` |
| 5g — Polymorphism | `modules/annotation_mapping.nf` | `create_polymorphism_worker.py` | `polymorphism.tsv` (rsid + ref/alt + allele freq from dbSNP 155) |
| 5h — PEM | `modules/annotation_mapping.nf` | `create_pem_worker.py`, `create_pem_transfer_worker.py` | `pem_core_motifs.tsv`, `pem_core_motifs_mapped.tsv` |
| 5i — Coiled Coils | `modules/annotation_mapping.nf` | `create_coiledcoils_worker.py` | `coiled_coils.tsv`, `DeepCoil.tsv` |
| 5j — PPI | `modules/annotation_mapping.nf` | `create_ppi_worker.py` | `interactions.tsv` |
| 5k — ScanSite | `modules/annotation_mapping.nf` | `create_scansite_worker.py` | `scansite.tsv` |
| 5m — Position-Based | `modules/annotation_mapping.nf` | `create_position_based_worker.py` | `position_based_annotations.tsv`, `rsa_scores.tsv` |
| 5n — ELM Classes | `modules/annotation_mapping.nf` | `create_elm_class_worker.py` | `elm_classes.tsv` |
| 5o — MobiDB | `modules/annotation_mapping.nf` | `create_mobidb_worker.py` | `mobidb_disorder.tsv` |
| 7 — Conservation | `modules/annotation_mapping.nf` | `create_conservation_worker.py` | `conservation_multiple_level.tsv`, `conservation_phastcons.tsv` |
| 8a — ClinVar disease | `modules/annotation_mapping.nf` | `create_clinvar_disease_build_worker.py` | `final/disease/clinvar_disease.tsv` |
| 8f — Pathogenicity | `modules/annotation_mapping.nf` | `create_dbnsfp_map_worker.py` | `final/pathogenicity/pathogenicity_scores.tsv` |
| 8g — ProteinGym | `modules/annotation_mapping.nf` | `create_proteingym_worker.py` | `proteingym.tsv` |
| 8h — FINCHES | `modules/annotation_mapping.nf` | `create_finches_worker.py` | `finches_saturation.tsv` (off by default; `--skip_finches false` to enable; CC BY-NC 4.0) |
| Report | `modules/annotation_mapping.nf` | `create_mapping_report_worker.py` | `mapping_reports/` (runs last) |
| Scatter | `modules/annotation_mapping.nf` (`SPLIT_SEQ_TABLE`) | `split_seq_table.py` | N gene-balanced seq-table chunks (`--scatter_chunks N`) |
| Reference fetches | `modules/fetch_references.nf` | — | UniProt/GENCODE/ClinVar/GO/MobiDB/MONDO/AlphaMissense/IntAct/BioGRID/HIPPIE (cached via `storeDir`) |

### Key Conventions

- **`Protein_ID`**: Gencode transcript name (e.g. `RAF1-201`). Primary key for all mapped outputs.
- **Exact gene subsetting**: GENCODE headers use `|RAF1|` (pipe-delimited), UniProt uses `GN=RAF1 ` (trailing space) — both required to avoid `TRAF1`, `ZTRAF1`
- **`storeDir` caching**: Reference downloads cached in `references/`; if a storeDir file is 0 bytes (failed download), delete it and re-run. `-stub` writes to `references/_stub/` and never pollutes the real cache.
- **Mutation input is mutually exclusive**: `--clinvar_vcf` OR `--mutation_maf` OR `--mutation_vcf`, not combined
- **TCGA MAF QC**: `--mutation_source TCGA` truncates barcodes to 12 chars; `--mutation_hypermutation_threshold 1500` drops hot samples; `--no_hgvsp_validation` disables ref-AA check
- **ClinVar disease build**: when `hg38_2bit` + `clinvar_disease_from_mutations=true` + `mondo_obo` set, `CLINVAR_DISEASE_BUILD` runs from `MUTATION_MAP` outputs (not a static filter table)
- **dbNSFP dual mode**: `--dbnsfp_raw_dir` → `DBNSFP_MAP` (raw `chr*.gz` via `combined_map.map`); `--dbnsfp_tsv` → `PATHOGENICITY_MAP` (pre-mapped Protein_ID-keyed TSV). Mutually exclusive; raw takes priority.
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
| GO | `FETCH_GO` (`goa_human.gaf.gz` + `go.obo`, cached) |
| ClinVar | `FETCH_CLINVAR` (NCBI FTP, cached in `references/clinvar/`) |
| PPI | `FETCH_INTACT/BIOGRID/HIPPIE` → `PPI_PREPROCESS` (cached in `references/ppi/`) |
| Conservation (GOPHER) | `params.gopher_conservation_table` — external pre-computed table |
| Conservation (phastCons) | `params.phastcons_dir` — chr*.bw files; requires `bigWigToBedGraph` in PATH |
| hg38.2bit | `params.hg38_2bit` or `--fetch_hg38_2bit true` |

### Pending Modules (low priority)

| Module | Data Source | Status |
|--------|-------------|--------|
| ELM Switches | elm.eu.org/switches | pending |
| FuzDrop LLPS probability | fuzdrop.bio.unipd.it API | pending |
| Complexity tracks (SEG/DUST/TRF) | local tools | pending |
