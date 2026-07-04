# Troubleshooting

Common failure modes and their fixes. Most empty-output problems trace back to
a missing external program or a stale Nextflow cache.

---

## Nextflow can't find Java

**Symptom**: `nextflow run ...` fails immediately with a Java error.

**Cause**: Nextflow needs a Java 11–21 runtime on `PATH`. In non-interactive
shells (cron, `ssh host 'cmd'`, CI, `nohup`) the conda activate hook that puts
Java on `PATH` does not run.

**Fix**: point Nextflow at the JVM inside the `discanvis` conda env explicitly:

```bash
conda activate discanvis
export JAVA_CMD="$CONDA_PREFIX/bin/java"
export JAVA_HOME="$CONDA_PREFIX"
export PATH="$CONDA_PREFIX/bin:$PATH"
nextflow run main.nf --project test_one_protein --target_gene RAF1 -resume
```

Interactive sessions that ran `conda activate discanvis` already have this — the
exports only matter for detached / scripted runs.

---

## IUPredscores / AIUPredscores are empty (header only)

**Cause**: `aiupred_python` points to a missing or wrong Python binary.

**Fix**:
1. Verify the env: `conda run -n discanvis_aiupred python -c "import iupred3_lib; print('OK')"`
2. Set it in `local.config`: `aiupred_python = '<env>/bin/python'`
3. Delete the cached `DISORDER_MAP` work dirs (Nextflow cached the empty result):
   ```bash
   find work -name ".command.sh" | xargs grep -l "create_disorder_worker" | \
       xargs -I{} dirname {} | xargs rm -rf
   ```
4. Re-run with `-resume`.

See [Installation § Disorder predictors](installation.md#4-disorder-predictors--external-programs).

---

## coiled_coils.tsv is empty

Same root cause as above but for DeepCoil. Set `deepcoil_python` in
`local.config`, or skip it with `--skip_coiledcoils true`.

---

## pfam_domains.tsv is empty

**Cause**: a stale `protein2ipr.dat.gz` parse cached under `references/`.

**Fix**: delete the cached parsed table and re-run:

```bash
rm references/uniprot_parsed/pfam_domains.tsv
nextflow run main.nf ... -resume
```

---

## conservation_phastcons.tsv is empty

**Cause**: `bigWigToBedGraph` is not on `PATH`.

**Fix**: install it from bioconda (already in `environment.yml`):

```bash
conda install -n discanvis -c bioconda ucsc-bigwigtobedgraph
```

Or set the full path in `local.config`: `bigwigtobedgraph = '/path/to/bigWigToBedGraph'`.

---

## Nextflow cached a task with wrong results

If a task produced an incorrect output (e.g. an empty file) but exited 0,
`-resume` will not re-run it even after you fix the code.

**Fix**: delete that specific work dir so Nextflow re-runs it:

```bash
# Find work dirs for a specific process
find work -name ".command.sh" | xargs grep -l "create_disorder_worker" | \
    xargs -I{} dirname {}
# Delete the offending one, then re-run with -resume
rm -rf work/<data>/XX/YYYYYYYY...
```

---

## storeDir file is 0 bytes (failed download)

A failed reference download can leave a 0-byte file that `storeDir` treats as
complete. Delete it and re-run:

```bash
find references/ -empty -name "*.tsv" -o -empty -name "*.gz"
rm <empty-file>
nextflow run main.nf ... -resume
```
