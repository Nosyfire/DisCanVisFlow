# Installation & Setup

Everything you need to get DisCanVisFlow running on a new machine: the conda
environment, the disorder-prediction helpers, and how references are supplied.

The [README](../README.md) quick start is enough for the common case (portable
mode auto-downloads open references). Read this page when you want local
reference files, the licensed disorder predictors, or a development setup.

---

## 1. Requirements

| Component | Notes |
|-----------|-------|
| [Nextflow](https://www.nextflow.io/) ≥ 24 | DSL2; installed into the conda env below |
| Java 11–21 | Provided by the conda env; see [Troubleshooting](troubleshooting.md#nextflow-cant-find-java) for non-interactive shells |
| conda / mamba | For the `discanvis` environment |
| Docker *(optional)* | Alternative to conda — see [Configuration guide § Docker](configuration_guide.md#docker) |

---

## 2. Conda environment

```bash
git clone https://github.com/Nosyfire/DisCanVisFlow
cd DisCanVisFlow

conda env create -f environment.yml
conda activate discanvis
```

If the environment was created before the bioconda UCSC tools were added:

```bash
conda env update -n discanvis -f environment.yml --prune
```

The environment includes Nextflow, BLAST+, BLAT, and the UCSC tools
(`bigBedToBed`, `bigWigToBedGraph`) used by the polymorphism and conservation
modules.

---

## 3. Reference data

References can be supplied two ways, selected with `--data`:

| Mode | Behavior | When to use |
|------|----------|-------------|
| `--data discanvis_data` *(default)* | Downloads open-access references on demand and caches them in `references/` via `storeDir` | New machine, CI, portable runs |
| `--data local` | Reads pre-existing paths from `config/data/local.config` | Reproducibility with frozen snapshots; shared server with references at a fixed path |

**Portable mode** needs no configuration — just run. Reference provenance,
versions, update cadence, and how to refresh the cache are documented in
[Reference data](reference_data.md).

**Local mode** — copy the template and fill in your paths:

```bash
cp config/data/local.config.template config/data/local.config
# then edit config/data/local.config with paths to:
#   UniProt FASTA, GENCODE FASTA/GTF, hg38.2bit,
#   External_Programs directory, and any optional
#   AlphaMissense / MaveDB / ProteinGym / DepMap / dbNSFP files
```

`config/data/local.config` is machine-specific and is **not** committed to git
(it is in `.gitignore`); every user maintains their own.

---

## 4. Disorder predictors & external programs

Disorder prediction relies on libraries that are not on conda. Run the setup
script once per machine:

```bash
bash bin/setup_external_programs.sh
```

This clones AIUPred, creates the `discanvis_aiupred` and `discanvis_deepcoil`
conda environments, and installs `bigBedToBed`. With `--data discanvis_data`,
the `SETUP_DEPS` Nextflow process performs this automatically on the first run.

### Programs and where they come from

| Tool | Needed for | How to obtain |
|------|-----------|---------------|
| **IUPred3 / ANCHOR2** | IUPred/ANCHOR scores | Academic licence — register at [iupred2a.elte.hu/download](https://iupred2a.elte.hu/download), extract into `External_Programs/iupred3/` |
| **AIUPred disorder** | AIUPred scores | Cloned by `setup_external_programs.sh` → `External_Programs/aiupred-caid3/` |
| **AIUPred-Binding** | Binding-region scores | Cloned by `setup_external_programs.sh` → `External_Programs/AIUPred/` |
| **DeepCoil** | Coiled-coil predictions | Set `deepcoil_python` in `local.config`; skip with `--skip_coiledcoils true` on CUDA 12+ hardware |
| **bigBedToBed** | Polymorphism (dbSNP bigBed) | Installed via `bioconda::ucsc-bigbedtobed` in `environment.yml` |
| **bigWigToBedGraph** | phastCons conservation | Installed via `bioconda::ucsc-bigwigtobedgraph` in `environment.yml` |

A missing predictor does **not** crash the run — its scores are simply left
empty. Skip predictors explicitly with `--skip_iupred true`,
`--skip_aiupred true`, or `--skip_coiledcoils true`.

### Python paths (local mode)

When using `--data local`, point each predictor at its Python interpreter:

```groovy
// config/data/local.config
params {
    ext_programs    = '/path/to/External_Programs'
    aiupred_python  = '/path/to/envs/aiupred/bin/python'      // needs scipy + torch
    deepcoil_python = '/path/to/envs/discanvis_deepcoil/bin/python'
}
```

`bigWigToBedGraph` is found automatically through the conda env — no override
needed. With `--data discanvis_data`, `SETUP_DEPS` auto-detects and writes these
paths for you.

---

## 5. Verify the install

Validate the workflow graph without running any computation:

```bash
nextflow run main.nf --project test_one_protein --target_gene RAF1 -stub
```

Then run the Python worker test suite (no Nextflow required — the workers are
called as subprocesses with dummy inputs):

```bash
conda activate discanvis
pytest tests/ -v
```

- Single file: `pytest tests/test_create_disorder_worker.py -v`
- Single test: `pytest tests/test_create_mutation_map_worker.py::TestMissenseFilter -v`

Once both pass, run a real single gene — see the
[README quick start](../README.md#quick-start).
