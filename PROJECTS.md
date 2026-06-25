# Run Configuration

DisCanVisFlow uses one configuration root: `config/`.

Runs are selected with separate Nextflow parameters:

```bash
nextflow run main.nf \
  --project cellular_vulnerability \
  --data discanvis_data \
  --machine laptop \
  -resume
```

## Axes

| Axis | Folder | Meaning |
|------|--------|---------|
| `--project` | `config/projects/` | Biological goal and annotation-track selection |
| `--data` | `config/data/` | Reference paths and download behavior |
| `--machine` | `config/machines/` | Runtime resources, parallelism, executor |
| `--env` | `config/envs/` | Software environment |

Defaults are set in `nextflow.config`:

```text
--project cellular_vulnerability
--data discanvis_data
--machine laptop
--env conda
```

## Current Main Run

For the Turbine cellular-vulnerability feature set on an 8 GB laptop:

```bash
export NXF_OPTS='-Xms256m -Xmx1g'
nextflow run main.nf \
  --project cellular_vulnerability \
  --machine laptop \
  --description "Q4 2026 Turbine feature run" \
  -resume
```

## Available Projects

| Project | Scope | Output |
|---------|-------|--------|
| `cellular_vulnerability` | Full proteome feature set for the Turbine model | `results/cellular_vulnerability/` |
| `full_discanvis` | Full DisCanVis database update | `results/discanvis/` |
| `discanvis` | Legacy full DisCanVis preset | `results/discanvis/` |
| `vep_benchmarking` | Variant-effect-predictor benchmark set | `results/vep_benchmarking/` |
| `test_one_protein` | Single-gene smoke test, default TP53 | `results/test_one_protein/` |
| `test_subset` | Small multi-gene regression set | `results/test_subset/` |

## Available Machines

| Machine | Intended use |
|---------|--------------|
| `laptop` | 8 GB laptop/WSL-safe mode; one 4 GB task at a time |
| `low` | Small local workstation |
| `medium` | Normal local workstation |
| `hard` | Large local workstation/server |
| `slurm` | Cluster execution |

Override memory explicitly with:

```bash
nextflow run main.nf --project cellular_vulnerability --machine laptop --ram '4 GB' -resume
```

## Data Sources

| Data | Meaning |
|------|---------|
| `discanvis_data` | Portable reference cache under `references/`; downloads open references when possible |
| `local` | Machine-specific local paths from `config/data/local_refs.config` |

## Notes

- `-resume` reuses cached steps.
- Override any project setting directly on the CLI: `--target_gene`, `--skip_pdb`, `--outdir`, individual data paths, etc.
- `-stub` validates the DAG without executing workers.
