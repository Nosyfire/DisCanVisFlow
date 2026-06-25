# Configuration

This is the only run-configuration tree.

Use separate Nextflow parameters for the three independent choices:

```bash
nextflow run main.nf \
  --project cellular_vulnerability \
  --data discanvis_data \
  --machine laptop \
  -resume
```

## Folders

| Folder | What belongs here |
|--------|-------------------|
| `projects/` | Biological goal, output directory, enabled/disabled annotation tracks |
| `data/` | Reference source policy: portable downloaded cache or machine-local paths |
| `machines/` | Runtime resources: memory, CPUs, parallelism, executor |
| `envs/` | Software environment: conda, docker, or current shell |
| `gene_lists/` | Optional project gene lists |

## Common Choices

```bash
# Current cellular-vulnerability run on an 8 GB laptop
export NXF_OPTS='-Xms256m -Xmx1g'
nextflow run main.nf --project cellular_vulnerability --machine laptop -resume

# Same biological project on a larger local machine
nextflow run main.nf --project cellular_vulnerability --machine hard -resume

# Explicit memory override
nextflow run main.nf --project cellular_vulnerability --machine laptop --ram '4 GB' -resume
```
