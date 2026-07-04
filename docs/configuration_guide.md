# Configuration Guide

This guide explains every flag you will encounter when running DisCanVisFlow.
Most flags have sensible defaults — you only need to touch them when you want something non-standard.

---

## The four axes

Every run is shaped by four independent settings:

```
--project   what to compute and where to put results
--machine   how many CPUs/RAM to use (matches your hardware)
--data      where reference files come from
--env       conda (default) or docker
```

These map directly to config files:

| Axis | Flag | Config loaded |
|------|------|--------------|
| Project | `--project test_one_protein` | `config/projects/test_one_protein.config` |
| Machine | `--machine medium` | `config/machines/medium.config` |
| Data | `--data local` | `config/data/local.config` |
| Environment | `--env conda` | `config/envs/conda.config` |

---

## `--data` — where do reference files come from?

**Default: `discanvis_data` (auto-download everything)**

You almost never need to specify this. The default mode downloads UniProt, GENCODE, ClinVar,
GO, hg38, AlphaMissense, and all other open-access references automatically on first run
and caches them in `references/` via Nextflow's `storeDir`. Subsequent runs reuse the cache.

```bash
# These two commands are identical:
nextflow run Nosyfire/DisCanVisFlow --project test_one_protein --machine medium --target_gene RAF1 -resume
nextflow run Nosyfire/DisCanVisFlow --project test_one_protein --machine medium --target_gene RAF1 --data discanvis_data -resume
```

**When to use `--data local`:**

Switch to `--data local` only when you have pre-existing reference files and want to avoid
re-downloading them (e.g. frozen snapshots for reproducibility, or a shared server where
references live at a fixed path).

```bash
# 1. Copy the template (one time per machine)
cp config/data/local.config.template config/data/local.config

# 2. Fill in the paths to your local files
#    local.config is in .gitignore — it is machine-specific, not committed to git

# 3. Run with --data local
nextflow run main.nf --project discanvis --machine hard --data local -resume
```

**Summary:**

| `--data` value | Meaning | When to use |
|----------------|---------|-------------|
| `discanvis_data` | Auto-download all open-access references | **Default — fresh machine, first run, portable** |
| `local` | Read paths from `config/data/local.config` | Server with pre-existing files, frozen snapshots |

---

## `--machine` — how much CPU/RAM to claim

Matches your hardware. Affects parallelism (`blat_chunks`, `scatter_chunks`, `queueSize`)
and per-process memory limits. **Default: `laptop`** (safe but slow).

| `--machine` | Total RAM declared | CPUs | Parallel BLAT jobs | Use when |
|-------------|-------------------|------|--------------------|---------|
| `laptop` | 5 GB | 2 | 1 | 8 GB RAM laptop, WSL, very constrained |
| `low` | 32 GB | 6 | 4 | 32 GB cluster node or low-RAM workstation |
| `medium` | 64 GB | 16 | 8 | Workstation, cluster node with 64 GB+ |
| `hard` | 256 GB | 64 | 32 | Dedicated server with 256 GB+ RAM (e.g. gpu0) |
| `slurm` | cluster-managed | — | 32 | SLURM HPC cluster |

BLAT loads the hg38.2bit genome file (~4 GB) per process. With `--machine hard`,
32 parallel BLAT jobs = ~128 GB RAM. On a machine without that, use `medium` or `laptop`.

> **`--project test_one_protein` always overrides to `blat_chunks=1` regardless of `--machine`.**
> A single gene has ~10 transcripts — there is nothing to parallelize.

---

## `--project` — what to run and where results go

| `--project` | What it does | Output directory |
|-------------|-------------|-----------------|
| `test_one_protein` | Single-gene smoke test (default gene: TP53) | `results/test_one_protein/` |
| `discanvis` | Full human proteome, all annotation tracks (legacy preset) | `results/discanvis/` |
| `full_discanvis` | Full DisCanVis database update | `results/discanvis/` |
| `cellular_vulnerability` | Turbine ML feature set (gene-list driven) | `results/cellular_vulnerability/` |
| `vep_benchmarking` | Variant-effect-predictor benchmark set | `results/vep_benchmarking/` |
| `test_subset` | 5-gene regression set (TP53, RAF1, BRAF, KRAS, EGFR) | `results/test_subset/` |
| `individual_run_raf1` | Single-gene real run (RAF1) with blast/ + intermediate/ + reports/ | `results/individual_run_raf1/` |
| `individual_run_ctbnn1` | Single-gene real run (CTNNB1) | `results/individual_run_ctbnn1/` |
| `subset_run_test` | Sliced multi-gene subset (RAF1, CTNNB1, TP53, BRAF) | `results/subset_run_test/` |

Override the default gene with `--target_gene RAF1`, or supply a list with `--gene_list_file`.

**Defaults** (from `nextflow.config`, applied when a flag is omitted):
`--project cellular_vulnerability --data discanvis_data --machine laptop --env conda`.

Any project setting can be overridden directly on the CLI (`--target_gene`,
`--skip_pdb`, `--outdir`, `--ram '4 GB'`, individual data paths, …). `-resume`
reuses cached steps; `-stub` validates the DAG without executing workers.

---

## `--modules` — run only what you need

By default, all annotation modules run. Use `--modules` with a comma-separated list to
include only specific groups. The backbone (BLAST, SEQUENCE_PROCESS, ANNOTATION_MAP,
TRANSCRIPT_MAP, ELM, Pfam, DIBS/MFIB/PTM) always runs regardless of `--modules`.

```bash
# Only mutations + disorder:
--modules mutations,disorder

# Only PDB coverage + GO terms + PPI:
--modules pdb,go,ppi
```

Available names: `mutations`, `disorder`, `mobidb`, `pdb`, `go`, `polymorphism`,
`pem`, `coiledcoils`, `ppi`, `conservation`, `scansite`, `clinvar_disease`, `omim`,
`cancer_drivers`, `alphamissense`, `depmap`, `mavedb`, `proteingym`, `dbnsfp`, `finches`

---

## `--fetch_cbioportal` — somatic mutations from cBioPortal

**Default in `test_one_protein`: `true`** (enabled automatically for single-gene tests)
**Default in other projects: `false`** (opt-in for full proteome runs, where API per-gene would be slow)

ClinVar germline variants are always included in the `mutations` module.
cBioPortal somatic mutations are opt-in because:
- The REST API takes ~1 min per gene — fine for one gene, slow for 20,000
- For full-proteome runs, use `--cbioportal_study <id>` to download one bulk MAF instead

```bash
# Single gene: API mode, no study ID needed (default in test_one_protein)
nextflow run Nosyfire/DisCanVisFlow --project test_one_protein --target_gene RAF1 -resume

# Full proteome: bulk download from a specific cohort
nextflow run main.nf --project discanvis --machine hard \
    --fetch_cbioportal true --cbioportal_study tcga_pan_can_atlas_2018 -resume

# Disable cBioPortal entirely
nextflow run main.nf --project discanvis --machine hard \
    --fetch_cbioportal false -resume
```

---

## Skip flags — disable individual predictors within a module

These apply *within* a module (complement to `--modules`, which controls whole groups):

| Flag | What it skips | When useful |
|------|--------------|------------|
| `--skip_iupred true` | IUPred3 + ANCHOR2 | Keep only AIUPred; IUPred needs a separate conda env |
| `--skip_aiupred true` | AIUPred disorder + binding | Keep only IUPred3 |
| `--skip_alphafold true` | AlphaFold pLDDT fetch (EBI API) | Speed up rerun; use `--alphafold_precomputed_table` to reuse prior scores |
| `--skip_pdb true` | PDB structure mapping | Skip if you only need sequence/disorder |
| `--skip_conservation true` | GOPHER + phastCons | Needs external files; skip if not configured |
| `--skip_polymorphism true` | dbSNP 155 SNPs | Skip for pure IDP/disease analysis |
| `--skip_coiledcoils true` | DeepCoil predictions | Skip on CUDA 12+ hardware without DeepCoil env |

---

## Per-process resource usage (RAF1 single-gene, `--modules mutations,disorder`)

Measured with `-with-trace` on a 64-CPU server. All annotation processes are sub-second
and use < 500 MB. The two processes worth watching on memory-constrained machines:

| Process | Duration | Peak RAM | Notes |
|---------|----------|----------|-------|
| `SUBSET_*` (FASTA subsetting) | < 1 s | < 35 MB | trivial |
| `MAKEBLASTDB_*` | < 1 s | < 35 MB | trivial |
| `BLASTP_*` | < 200 ms | < 36 MB | single gene subset |
| `BLAT_ALIGN` (per chunk) | ~27 s | **4 GB** | loads full hg38.2bit — the memory bottleneck |
| `SEQUENCE_PROCESS` | 6.7 s | 310 MB | fine |
| `ANNOTATION_MAP` | 4.3 s | 459 MB | fine |
| `DISORDER_MAP` | 8.2 s | 844 MB | AIUPred model in memory |
| `MUTATION_MAP_CLINVAR` | 32 s | **5.1 GB** | scans VCF against genome map |
| `MUTATION_MAP_CBIOPORTAL` | 0.5 s | 24 MB | fast (small MAF for one gene) |
| All other processes | < 1 s | < 30 MB | trivial |

**Takeaway:** You need ~6 GB free RAM per BLAT chunk + ~5 GB for ClinVar mutation mapping.
With `blat_chunks=1` (enforced by `test_one_protein`), the peak footprint is ~6 GB total.
`--machine medium` (64 GB declared) is the right default for any machine with ≥ 16 GB free.

---

## Common flag combinations

The basic single-gene and full-proteome commands are in the
[README quick start](../README.md#quick-start). The combinations below show the
non-trivial variations.

```bash
# Fast focused run: only disorder + mutations (no PDB, GO, conservation, PPI...)
nextflow run main.nf \
    --project test_one_protein --machine medium --target_gene RAF1 \
    --modules mutations,disorder --skip_iupred true -resume

# Full proteome on local server (paths in local.config)
nextflow run main.nf \
    --project discanvis --machine hard --data local -resume

# Gene list, portable, medium server
nextflow run main.nf \
    --project discanvis --machine medium \
    --gene_list_file my_genes.txt -resume

# Validate pipeline wiring without running anything
nextflow run Nosyfire/DisCanVisFlow -latest \
    --project test_one_protein --target_gene RAF1 -stub
```

---

## Mutation inputs (ClinVar / MAF / VCF)

The `mutations` module accepts one mutation source. `--clinvar_vcf`,
`--mutation_maf`, and `--mutation_vcf` are **mutually exclusive**.

```bash
# ClinVar VCF override (instead of the auto-downloaded ClinVar)
nextflow run main.nf --project test_one_protein --machine hard --target_gene TP53 \
    --clinvar_vcf /path/to/clinvar.vcf.gz -resume

# TCGA / cBioPortal MAF
nextflow run main.nf --project test_one_protein --machine hard --target_gene TP53 \
    --mutation_maf /path/to/tcga.maf --mutation_source TCGA -resume

# Custom VCF
nextflow run main.nf --project test_one_protein --machine hard --target_gene TP53 \
    --mutation_vcf /path/to/variants.vcf.gz --mutation_source MyStudy -resume
```

For cBioPortal somatic mutations, see
[`--fetch_cbioportal`](#--fetch_cbioportal--somatic-mutations-from-cbioportal) above.

---

## Re-running disorder with pre-computed AlphaFold

Fetching AlphaFold pLDDT from the EBI API is the slow part of `DISORDER_MAP`
(~8 h for the full proteome). Reuse a prior run's scores to recompute
IUPred/AIUPred without re-fetching:

```bash
nextflow run main.nf --project discanvis --machine hard \
    --alphafold_precomputed_table results/discanvis/final/disorder/AlphaFoldTable.tsv \
    --skip_alphafold true \
    -resume
```

---

## Projects: primary run and derived projects

The canonical starting point is always **raw data + one full pipeline run**.
Secondary projects (subsets, alternative views) are **derived** from that run —
not computed separately.

### Primary run

```bash
nextflow run main.nf --project discanvis --machine hard -resume
```

Produces `results/discanvis/final/` with every annotation TSV for the full
proteome.

### Derived projects (no re-computation)

After `discanvis` completes, one script produces all derived directories:

```bash
bash bin/derive_projects_from_discanvis.sh            # default source: results/discanvis
bash bin/derive_projects_from_discanvis.sh results/discanvis
```

| Project | Method | Output |
|---------|--------|--------|
| `vep_benchmarking` | full copy | `results/vep_benchmarking/` |
| `cellular_vulnerability` | selective full-proteome copy (annotations, sequence, drivers, dbnsfp, combined disorder, alphamissense, DepMap) | `results/cellular_vulnerability/` |
| `test_subset` | 5-gene extraction (TP53, RAF1, BRAF, KRAS, EGFR) | `results/test_subset/` |

Ad-hoc extraction of any gene set from a completed run (filters all TSVs by
`Protein_ID` prefix — no pipeline re-run):

```bash
python bin/extract_gene_from_results.py \
    --source results/discanvis --gene RAF1,BRAF,KRAS --out results/kinase_subset

python bin/extract_gene_from_results.py \
    --source results/discanvis --gene_list_file config/gene_lists/my_genes.txt \
    --out results/my_gene_subset
```

### Re-run vs. re-derive

| Scenario | Action |
|----------|--------|
| Add a new annotation track | Re-run `discanvis` with `-resume`, then re-derive |
| Update a reference (e.g. new ClinVar) | `bin/refresh_refs.sh clinvar`, re-run with `-resume`, then re-derive |
| Add a gene to a derived subset | Re-derive from existing `discanvis` (no pipeline re-run) |

---

## SLURM cluster

```bash
nextflow run main.nf --project discanvis --machine slurm \
    --description "Full proteome — $(date +%Y-%m)" -resume
```

---

## Docker

```bash
docker build -t discanvis-pipeline:latest .
nextflow run main.nf --project test_one_protein --machine hard --env docker \
    --target_gene RAF1 -resume
```

---

## Configuration file layout

```
config/
├── projects/
│   ├── discanvis.config              Full proteome — all annotation tracks
│   ├── cellular_vulnerability.config Turbine ML features subset
│   ├── vep_benchmarking.config       VEP benchmark set
│   ├── test_one_protein.config       Single-gene smoke test
│   └── test_subset.config            Multi-gene regression (5 genes)
├── machines/
│   ├── hard.config                   Large server (64+ CPUs)
│   ├── medium.config                 Workstation (64 GB / 16 CPUs)
│   ├── low.config                    Low-RAM node (32 GB / 6 CPUs)
│   ├── laptop.config                 Memory-safe, 8 GB limit
│   └── slurm.config                  SLURM cluster
├── data/
│   ├── local.config                  Machine-specific paths (NOT in git)
│   ├── local.config.template         Template — copy and fill in paths
│   └── discanvis_data.config         Portable; auto-downloads references
└── envs/
    ├── conda.config
    └── docker.config
```
