# DisCanVisFlow — Disease & Disorder Annotation for Human Protein Isoforms

> A Nextflow DSL2 pipeline that maps disease variants, functional annotations, and structural features onto every curated protein isoform in the human SwissProt proteome. Designed to power the DisCanVis2 web server, but fully usable as a standalone data-generation pipeline for proteomics, structural biology, or ML feature pipelines.

---

## What it does

For each human protein (UniProt SwissProt × GENCODE):

1. **Maps each GENCODE transcript to its best UniProt isoform** via reciprocal BLASTP
2. **Builds a per-residue coordinate map** (protein pos ↔ codon ↔ hg38 position) via BLAT
3. **Runs 20+ annotation modules** — all DB-ready, Protein_ID-keyed TSVs:

| Category | Annotations |
|----------|-------------|
| Mutations | ClinVar (pathogenic/likely-pathogenic), TCGA MAF, cBioPortal MAF, custom VCF |
| Disorder | IUPred3, ANCHOR2, AIUPred disorder, AIUPred-Binding, AlphaFold pLDDT, Combined disorder |
| SLiMs & PTMs | ELM motifs, DIBS, MFIB, PhasePro, PTMdb, PhosphoSite, Pfam domains, UniProt ROI/binding |
| Structure | PDB coverage, unobserved regions, RSA scores |
| Polymorphism | dbSNP 155 common SNPs + allele frequencies |
| Pathogenicity | dbNSFP (raw chr*.gz or pre-mapped), AlphaMissense, MaveDB, ProteinGym |
| Disease | ClinVar disease ontology (MONDO), OMIM disease + mutations |
| Interactions | IntAct, BioGRID, HIPPIE |
| Gene function | GO terms (GOA), ScanSite phospho motifs, PEM core motifs |
| Conservation | GOPHER multi-level, phastCons per-residue |
| Cancer | CGC census, Compendium, DepMap somatic mutations |

All outputs use `Protein_ID` (GENCODE transcript name, e.g. `RAF1-201`) as the primary key and land in `results/<project>/final/`.

---

## Quick start

### 1. Clone and create conda environment

```bash
git clone https://github.com/Nosyfire/DisCanVisFlow
cd DisCanVisFlow

conda env create -f environment.yml
conda activate discanvis
```

### 2. Run a single-gene test (fastest, ~4 min on 64-CPU server)

> **Note:** Nextflow caches pipeline revisions locally. Use `-latest` to always pull the current version from GitHub — recommended on fresh machines and after updates.

**Full annotation run — all tracks:**

```bash
nextflow run Nosyfire/DisCanVisFlow -latest \
    --project test_one_protein \
    --data discanvis_data \
    --machine hard \
    --target_gene RAF1 \
    -resume
```

**Focused run — include only specific annotations (`--modules`):**

Use `--modules` to name exactly which annotation groups to run. Everything else is skipped.
The example below runs RAF1 with cBioPortal + ClinVar mutations, AIUPred disorder + binding prediction, and ELM motifs:

```bash
nextflow run Nosyfire/DisCanVisFlow -latest \
    --project test_one_protein \
    --data discanvis_data \
    --machine hard \
    --target_gene RAF1 \
    --modules mutations,disorder \
    --fetch_cbioportal true \
    --skip_iupred true \
    -resume
```

| Flag | Effect |
|------|--------|
| `--modules mutations,disorder` | Run only mutation mapping + disorder prediction; skip PDB, conservation, GO, PPI, etc. |
| `--fetch_cbioportal true` | Fetch RAF1 somatic mutations from cBioPortal across all public studies via REST API — no study ID required |
| `--skip_iupred true` | Within the disorder module, run AIUPred only — skip IUPred3/ANCHOR2 |

To use a specific cBioPortal study bundle instead of the API, add `--cbioportal_study <datahub_id>` (e.g. `tcga_pan_can_atlas_2018`). Study-bundle mode is better for full-proteome runs where fetching per-gene via API would be slow.

ELM motifs (`annotations/elm.tsv`) are always produced as part of the annotation backbone regardless of `--modules`.

`--data discanvis_data` downloads all references automatically (UniProt/GENCODE/ClinVar/GO/etc.) and caches them in `references/` for all future runs.

### 3. Validate the DAG without running

```bash
nextflow run Nosyfire/DisCanVisFlow --project test_one_protein --data local --machine laptop --target_gene RAF1 -stub
```

---

## New machine setup

### Option A — Portable (`--data discanvis_data`)

No config editing needed. All open-access references download on first run and are re-used via `storeDir` in `references/`. Disorder predictors and external programs must still be set up (see below).

```bash
nextflow run main.nf --project discanvis --data discanvis_data --machine hard -resume
```

### Option B — Local paths (`--data local`)

Use pre-existing local files to avoid downloading large references. Copy the template and fill in your paths:

```bash
cp config/data/local.config.template config/data/local.config
```

Then edit `config/data/local.config` with your machine's actual file paths for:
- UniProt FASTA, GENCODE FASTA/GTF
- hg38.2bit (for genome/mutation/polymorphism mapping)
- External_Programs directory (IUPred3, AIUPred, etc.)
- Optional: AlphaMissense, MaveDB, ProteinGym, DepMap, dbNSFP paths

> **Important**: `config/data/local.config` is machine-specific and is NOT committed to git. It is listed in `.gitignore`. Every collaborator maintains their own.

---

## External programs (disorder predictors)

Disorder prediction requires libraries not in conda. Run the setup script once per machine:

```bash
bash bin/setup_external_programs.sh
```

This clones AIUPred from GitHub, creates `discanvis_aiupred` and `discanvis_deepcoil` conda environments, and installs `bigBedToBed`. With `--data discanvis_data`, the `SETUP_DEPS` Nextflow process handles this automatically on first run.

### Manual setup

| Tool | Required | Setup |
|------|----------|-------|
| **IUPred3 / ANCHOR2** | For IUPred/ANCHOR scores | Register at [iupred2a.elte.hu/download](https://iupred2a.elte.hu/download) → extract into `External_Programs/iupred3/` |
| **AIUPred disorder** | For AIUPred scores | Cloned by `setup_external_programs.sh` → `External_Programs/aiupred-caid3/` |
| **AIUPred-Binding** | For binding-region scores | Cloned by `setup_external_programs.sh` → `External_Programs/AIUPred/` |
| **DeepCoil** | For coiled-coil predictions | Set `deepcoil_python` in local.config; skip with `--skip_coiledcoils true` on CUDA 12+ hardware |
| **bigBedToBed** | For polymorphism (dbSNP bigBed) | Installed via `bioconda::ucsc-bigbedtobed` in `environment.yml` — available automatically |
| **bigWigToBedGraph** | For phastCons conservation | Installed via `bioconda::ucsc-bigwigtobedgraph` in `environment.yml` — available automatically |

If a disorder predictor is unavailable, the scores for that predictor will be empty (the pipeline does not crash). Set `--skip_iupred true`, `--skip_aiupred true`, or `--skip_coiledcoils true` to explicitly skip.

### Python paths in local.config

When using `--data local`, specify the Python binary for each predictor:

```groovy
// config/data/local.config
params {
    ext_programs    = '/path/to/External_Programs'
    aiupred_python  = '/path/to/envs/aiupred/bin/python'     // has scipy + torch
    deepcoil_python = '/path/to/envs/discanvis_deepcoil/bin/python'
}
```

`bigWigToBedGraph` is found automatically via the conda env (installed as `ucsc-bigwigtobedgraph` from bioconda). No path override needed.

With `--data discanvis_data`, the `SETUP_DEPS` process auto-detects and writes the Python paths.

---

## Project derivation model

The canonical starting point is always **raw data + a full pipeline run**. Secondary projects (subsets, alternative views) are **derived** from the primary run — not re-run separately.

### Primary run: `discanvis`

Runs the full pipeline on the complete human proteome:

```bash
nextflow run main.nf --project discanvis --data local --machine hard -resume
```

Produces `results/discanvis/final/` with all annotation TSVs (~70 files, full proteome).

### Derived projects (no re-computation)

After `discanvis` completes, run one script to produce all derived directories:

```bash
bash bin/derive_projects_from_discanvis.sh
# or specify source:
bash bin/derive_projects_from_discanvis.sh results/discanvis
```

This generates:

| Project | Method | Output |
|---------|--------|--------|
| `vep_benchmarking` | `rsync` full copy | `results/vep_benchmarking/` |
| `cellular_vulnerability` | Selective full-proteome copy: annotations, sequence, drivers, dbnsfp, combined disorder, alphamissense, DepMap mutations | `results/cellular_vulnerability/` |
| `test_subset` | 5-gene extraction (TP53, RAF1, BRAF, KRAS, EGFR) | `results/test_subset/` |
| `raf1_example` | Single-gene extraction | `results/raf1_example/` |

Extraction uses `bin/extract_gene_from_results.py`, which filters all TSVs by `Protein_ID` prefix:

```bash
# Custom extraction
python bin/extract_gene_from_results.py \
    --source results/discanvis \
    --gene   RAF1,BRAF,KRAS \
    --out    results/kinase_subset

# From gene list file (custom subset example)
python bin/extract_gene_from_results.py \
    --source results/discanvis \
    --gene_list_file config/gene_lists/my_genes.txt \
    --out    results/my_gene_subset
```

### When to re-run the pipeline vs. re-derive

| Scenario | Action |
|----------|--------|
| Add a new annotation track | Re-run `discanvis` with `-resume` → re-derive |
| Update a reference dataset (e.g. new ClinVar) | `bin/refresh_refs.sh clinvar` → re-run with `-resume` → re-derive |
| Add a new gene to cellular_vulnerability | Re-derive from existing discanvis (no pipeline re-run) |
| Update disorder predictors (new IUPred version) | Re-run `discanvis --skip_alphafold true --alphafold_precomputed_table results/discanvis/final/disorder/AlphaFoldTable.tsv -resume` → re-derive |

---

## Common run commands

### Single-gene test

```bash
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 -resume
```

### Include only specific annotation modules (`--modules`)

Use `--modules` with a comma-separated list of module names to run only what you need.
When `--modules` is set, the backbone (BLAST, ID mapping, sequence processing, ANNOTATION_MAP)
always runs; only the named optional modules are added on top.

Available module names: `mutations`, `disorder`, `mobidb`, `pdb`, `go`, `polymorphism`,
`pem`, `coiledcoils`, `ppi`, `conservation`, `scansite`, `clinvar_disease`, `omim`,
`cancer_drivers`, `alphamissense`, `depmap`, `mavedb`, `proteingym`, `dbnsfp`, `finches`

```bash
# cBioPortal + ClinVar mutations + AIUPred disorder/binding + ELM (for RAF1)
# --fetch_cbioportal without --cbioportal_study uses the public API (no study ID needed)
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 \
    --modules mutations,disorder \
    --fetch_cbioportal true \
    --skip_iupred true \
    -resume

# Disorder + PDB unobserved regions only
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene TP53 \
    --modules disorder,pdb \
    -resume

# Mutations + GO terms + PPI
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene EGFR \
    --modules mutations,go,ppi \
    -resume
```

To skip individual predictors *within* a module (e.g. skip IUPred3 but keep AIUPred):

```bash
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene RAF1 \
    --skip_alphafold true --skip_iupred true -resume
```

### Re-run disorder only (with pre-computed pLDDT from a prior run)

```bash
nextflow run main.nf --project discanvis --data local --machine hard \
    --alphafold_precomputed_table results/discanvis/final/disorder/AlphaFoldTable.tsv \
    --skip_alphafold true \
    -resume
```

This avoids re-fetching AlphaFold from the EBI API (~8 hours for full proteome) while recomputing IUPred/AIUPred.

### Custom VCF / MAF mutations

```bash
# ClinVar VCF override
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene TP53 \
    --clinvar_vcf /path/to/clinvar.vcf.gz -resume

# TCGA MAF
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene TP53 \
    --mutation_maf /path/to/tcga.maf --mutation_source TCGA -resume

# Custom VCF
nextflow run main.nf --project test_one_protein --data local --machine hard --target_gene TP53 \
    --mutation_vcf /path/to/variants.vcf.gz --mutation_source MyStudy -resume
```

Note: `--clinvar_vcf`, `--mutation_maf`, and `--mutation_vcf` are mutually exclusive.

### Gene list from file

```bash
nextflow run main.nf --project discanvis --data local --machine hard \
    --gene_list_file config/gene_lists/my_genes.txt -resume
```

### SLURM cluster

```bash
nextflow run main.nf --project discanvis --data local --machine slurm \
    --description "Full proteome — $(date +%Y-%m)" -resume
```

### Docker

```bash
docker build -t discanvis-pipeline:latest .
nextflow run main.nf --project test_one_protein --data local --machine hard --env docker \
    --target_gene RAF1 -resume
```

---

## Reference data management

### List cached references

```bash
bin/refresh_refs.sh
python bin/generate_manifest.py --no_checksum   # writes references/MANIFEST.tsv
```

### Force re-download

```bash
bin/refresh_refs.sh clinvar             # ClinVar only
bin/refresh_refs.sh clinvar mobidb go   # multiple
bin/refresh_refs.sh all                  # everything except hg38/dbsnp/alphafold
bin/refresh_refs.sh --force all          # truly everything
```

Then `-resume` — only deleted files re-download.

### Two data modes

| Mode | Behavior | When to use |
|------|----------|-------------|
| `--data discanvis_data` | Downloads open references on demand | New machine, CI, portable |
| `--data local` | Reads pre-existing paths from `config/data/local.config` | Reproducibility with frozen snapshots |

---

## Output structure

```
results/<project>/
├── final/
│   ├── annotations/     elm.tsv, dibs.tsv, mfib.tsv, phasepro.tsv, ptm_merged.tsv,
│   │                    pfam_domains.tsv, go_terms.tsv, polymorphism.tsv,
│   │                    interactions.tsv, scansite.tsv, pem_core_motifs.tsv, coiled_coils.tsv …
│   ├── disorder/        IUPredscores.tsv, Anchorscores.tsv, AIUPredscores.tsv,
│   │                    AIUPredBinding.tsv, AlphaFoldTable.tsv, CombinedDisorderNew.tsv …
│   ├── genome/          combined_map.map, exon.tsv
│   ├── mutations/       Missense/Frameshift/Nonsense/Indel_filter_mutations_mapped.tsv …
│   ├── pathogenicity/   dbnsfp_scores.tsv, alphamissense.tsv, mavedb.tsv, proteingym.tsv
│   ├── pdb/             pdb_structures.tsv, pdb_missing.tsv
│   ├── disease/         clinvar_disease.tsv, omim_disease.tsv
│   ├── drivers/         cancer_driver.tsv, census_driver.tsv, compendium_driver.tsv
│   ├── conservation/    conservation_multiple_level.tsv, conservation_phastcons.tsv
│   ├── position/        position_based_annotations.tsv, rsa_scores.tsv
│   └── sequence/        loc_chrom_with_names_isoforms_with_seq.tsv, isoform_alignment.tsv …
├── intermediate/        Entry_Isoform-keyed staging TSVs (input to TRANSCRIPT_MAP)
└── mapping_reports/
    ├── mapping_summary.md        Run metadata, annotation sources, per-run coverage
    └── mapping_coverage.tsv      Per-(Gene × annotation) coverage matrix

work/local/              Nextflow task cache — --data local
work/discanvis_data/     Nextflow task cache — --data discanvis_data
references/              storeDir-cached downloads (shared across all runs)
```

**Cross-project cache sharing**: Two projects using the same `--data` flag share `work/<data>/`. Nextflow `-resume` automatically reuses tasks whose inputs are unchanged — so `cellular_vulnerability` and `discanvis` (both `--data local`) share BLAST, genome mapping, and reference downloads.

---

## Running tests

```bash
conda activate discanvis

# All tests
pytest tests/ -v

# Single module
pytest tests/test_create_disorder_worker.py -v

# Specific test function
pytest tests/test_create_mutation_map_worker.py::TestMissenseFilter -v
```

Tests call `bin/*.py` scripts as subprocesses with dummy data in `tests/dummy_data/`. No Nextflow required.

---

## Architecture: modular design

Every computation step is a standalone Python worker in `bin/`:

```
bin/create_disorder_worker.py      # IUPred3, AIUPred, AlphaFold pLDDT
bin/create_annotation_worker.py    # ELM, DIBS, MFIB, PTM, Pfam
bin/create_mutation_map_worker.py  # ClinVar/MAF/VCF → protein positions
bin/create_transcript_map_worker.py # UniProt-keyed → all isoforms
...
```

Each `bin/*.py` has `argparse` and can be called directly for debugging:

```bash
# Debug a single protein
python bin/create_disorder_worker.py \
    --loc_chrom results/discanvis/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv \
    --ext_programs /path/to/External_Programs \
    --aiupred_python /opt/anaconda3/envs/aiupred/bin/python \
    --output_dir /tmp/test_disorder \
    --skip_alphafold
```

Nextflow acts purely as orchestrator: it handles caching, parallelism, and data flow. The workers are tested independently via `pytest tests/`.

---

## Annotation sources — per-release update cadence

| Source | Update method | Freeze / always-current |
|--------|---------------|------------------------|
| UniProt SwissProt | `bin/refresh_refs.sh uniprot` | Frozen in `local.config` |
| GENCODE | `bin/refresh_refs.sh gencode` | Pinned to v44 by default |
| ClinVar | `bin/refresh_refs.sh clinvar` | Always-current via FETCH_CLINVAR |
| GO (GOA + OBO) | `bin/refresh_refs.sh go` | Always-current via FETCH_GO |
| MobiDB | `bin/refresh_refs.sh mobidb` | Always-current via FETCH_MOBIDB |
| ELM instances | `legacy_data/elm/elm_instances-2023.tsv` | Frozen 2023 snapshot |
| dbSNP bigBed | manual — see `bin/refresh_refs.sh dbsnp` | Large; rarely updated |
| AlphaMissense | `bin/refresh_refs.sh alphamissense` | v2023 frozen |
| dbNSFP | `--dbnsfp_raw_dir` or `--dbnsfp_tsv` | External; update manually |
| PPI (IntAct/BioGRID/HIPPIE) | `FETCH_INTACT/BIOGRID/HIPPIE + PPI_PREPROCESS` | Auto on first run |

---

## Configuration layout

```
config/
├── projects/
│   ├── discanvis.config            Full proteome — all annotation tracks
│   ├── cellular_vulnerability.config  Turbine ML features subset
│   ├── vep_benchmarking.config     VEP benchmark set
│   ├── test_one_protein.config     Single-gene smoke test
│   └── test_subset.config          Multi-gene regression (5 genes)
├── machines/
│   ├── hard.config                 Large server (64+ CPUs, 1+ TB RAM)
│   ├── medium.config               Workstation
│   ├── laptop.config               Memory-safe, 8 GB limit
│   └── slurm.config                SLURM cluster
├── data/
│   ├── local.config                Machine-specific paths (NOT in git)
│   ├── local.config.template       Template — copy and fill in paths
│   └── discanvis_data.config       Portable; auto-downloads references
└── envs/
    ├── conda.config
    └── docker.config
```

---

## Troubleshooting

### IUPredscores / AIUPredscores are empty (header only)

**Cause**: `aiupred_python` points to a non-existent or wrong Python binary.

**Fix**:
1. Verify the correct env: `conda run -n aiupred python -c "import iupred3_lib; print('OK')"`
2. Set it in `local.config`: `aiupred_python = '/opt/anaconda3/envs/aiupred/bin/python'`
3. Delete the cached work dirs for DISORDER_MAP (Nextflow cached the wrong results):
   ```bash
   find work/local -name ".command.sh" | xargs grep -l "create_disorder_worker" | \
       xargs -I{} dirname {} | xargs rm -rf
   ```
4. Re-run with `-resume`.

### coiled_coils.tsv is empty

Same root cause as above but for DeepCoil. Set `deepcoil_python` in `local.config`.

### pfam_domains.tsv is empty

**Cause**: `parse_uniprot_dat_worker.py` had wrong column indices for `protein2ipr.dat.gz`.
**Fix**: Already fixed in this codebase (commit `7fdfb24`). Delete the storeDir cache and re-run:
```bash
rm references/uniprot_parsed/pfam_domains.tsv
nextflow run main.nf ... -resume
```

### conservation_phastcons.tsv is empty

**Cause**: `bigWigToBedGraph` not in PATH.
**Fix**: Install from bioconda (already in `environment.yml`):
```bash
conda install -n discanvis -c bioconda ucsc-bigwigtobedgraph
```
Or if not using conda, set the full path in `local.config`: `bigwigtobedgraph = '/path/to/bigWigToBedGraph'`.

### Nextflow caches a task with wrong results

If a task produced an incorrect output (e.g. empty file) but exit code was 0, Nextflow `-resume` will not re-run it even after fixing the code.

**Fix**: Delete the specific work dir so Nextflow re-runs it:
```bash
# Find work dirs for a specific process
find work/local -name ".command.sh" | xargs grep -l "create_disorder_worker" | \
    xargs -I{} dirname {}
# Delete them
rm -rf work/local/XX/YYYYYYYY...
```

Then re-run with `-resume`.

### storeDir file is 0 bytes (failed download)

```bash
find references/ -empty -name "*.tsv" -o -empty -name "*.gz"
rm <empty-file>
nextflow run main.nf ... -resume
```

---

## Performance notes

Benchmarks on `gpu0.dlab.elte.hu` (64 CPUs, 1.4 TB RAM):

| Task | Time |
|------|------|
| Single gene (RAF1), all tracks | ~4 min |
| Full proteome (~20k genes), all tracks | ~24 h |
| BLAST (full proteome) | ~2 h |
| DISORDER_MAP (20 chunks, iupred+aiupred+alphafold) | ~8 h |
| DBNSFP_MAP (raw chr*.gz, 20 chunks) | ~3 h |

Key tips:
- Always use `--pdb_bulk true` (SIFTS join, ~10 min vs. 9 h for per-protein API)
- Use `--dbnsfp_tsv` (pre-mapped) instead of `--dbnsfp_raw_dir` for development runs
- Set `scatter_chunks=20` to parallelize DISORDER_MAP, DBNSFP_MAP, COILEDCOILS_MAP
- AlphaFold re-fetch takes ~25 min per chunk; use `--alphafold_precomputed_table` to skip if already done

---

## Citation

If you use this pipeline, please cite the tools and databases listed in `CITATIONS.md`.
