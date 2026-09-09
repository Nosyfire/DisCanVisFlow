---
name: discanvis-troubleshoot
description: >
  Diagnose DisCanVisFlow runs that failed, hung, or produced empty output. Use this skill
  when a pipeline run goes wrong in any way: "the run crashed", "nextflow can't find java",
  "my output TSV is empty / header only", "the process was killed", "it says exit status
  137", "conservation is empty", "IUPred produced nothing", "the download is 0 bytes",
  "it's been running for hours", "-resume isn't picking up my fix", "OutOfMemoryError",
  "storeDir file is corrupt", or when someone pastes a Nextflow error trace. Empty output
  is the important case — the pipeline is built to degrade rather than crash, so a missing
  track usually means a skipped or failed step rather than a genuine biological absence,
  and the run still reports success. Use it for diagnosing runs; to generate data use
  idp-dataset, to interpret existing data use discanvis-analysis.
---

# Diagnosing a DisCanVisFlow run

The pipeline is deliberately tolerant: a missing optional tool or an unreachable database
makes it skip a track rather than abort the run. That is good for finishing a 24-hour job
and bad for noticing problems, because **a run can report success and still be missing
tracks**. Empty output is therefore the most common symptom and almost never means the
biology is empty.

Work from evidence rather than guesses — Nextflow keeps the failing command, its output,
and its exit code on disk, so there is rarely a need to speculate.

## Step 1 — Establish the environment

```bash
python bin/find_discanvis.py
```

This alone explains a surprising share of failures: no conda environment, `nextflow` not
on `PATH`, or `local.config` absent while the command used `--data local`.

## Step 2 — Read the actual failure

Nextflow prints the work directory of the failing task. That directory is the diagnosis:

```bash
cd work/<data>/<hash-prefix>/<hash-rest>/
cat .command.sh    # exactly what ran
cat .command.err   # why it failed
cat .command.log   # combined output
cat .exitcode
```

Exit codes worth recognising immediately:

| Code | Meaning |
|---|---|
| 137 | Killed — out of memory, almost always. Not a bug in the worker. |
| 127 | Command not found — a tool is missing from the environment. |
| 1 | Ordinary error; read `.command.err`. |

To find the work directory again later:

```bash
grep -rn "process > <PROCESS_NAME>" .nextflow.log | tail -5
```

## Step 3 — Match the symptom

### Nextflow cannot find Java

Typical in cron, ssh, `nohup` and CI, where the conda activation hook that puts Java on
`PATH` never runs, so the launcher falls back to a stale interpreter path. Export the
environment's JVM explicitly:

```bash
export JAVA_CMD="$CONDA_PREFIX/bin/java" JAVA_HOME="$CONDA_PREFIX"
export PATH="$CONDA_PREFIX/bin:$PATH"
```

### A track is empty (header only)

Ask in this order — the first two account for most cases:

1. **Was the track even requested?** A `--modules` list excludes everything not named, and
   an unrecognised name is accepted silently. Check the run command before anything else.
2. **Was it skipped?** Search the run log for the process name. No mention means it never
   ran, which is a different problem from running and producing nothing.
3. **Is the required tool or input present?** Several tracks degrade quietly:
   `coiled_coils` needs DeepCoil, `low_complexity` needs `segmasker`, `dssp` needs
   `mkdssp` plus an AlphaFold model, `polymorphism` needs `bigBedToBed`, `conservation`
   needs external phastCons files, IUPred3 needs its licensed library in
   `External_Programs/`, and anything genome-anchored needs `params.hg38_2bit`.
4. **Did the input have the gene at all?** A curated database with no entry for the protein
   legitimately yields nothing — but confirm that before reporting it.

`docs/guide/troubleshooting.md` covers the specific empty-output cases in more detail.

### A downloaded reference is 0 bytes

`storeDir` caches by path, not by content, so a failed download is cached as a successful
one and every later run reuses the empty file. Delete it and re-run with `-resume`:

```bash
find references/ -size 0 -type f
rm references/<the-empty-file>
```

`bin/refresh_refs.sh <source>` does this properly for a whole source, and
`bin/refresh_refs.sh` with no arguments lists what is cached with sizes and dates — an
easy way to spot a truncated download.

### Out of memory (exit 137)

BLAT loads the ~4 GB genome per parallel job, so the peak is roughly
`blat_chunks × 4 GB`. `--machine hard` runs 32 of them, which needs about 128 GB. Drop to
`--machine medium` or `--machine laptop`. ClinVar mutation mapping separately needs around
5 GB. `docs/guide/performance.md` has the measured per-process figures.

### `-resume` ignores a fix

Nextflow caches on input hashes, so editing a worker invalidates only the tasks whose
inputs changed. If a task was cached with a wrong result, remove its work directory and
re-run with `-resume`. Note that `storeDir` outputs are *not* task cache — deleting the
work directory does not re-download a cached reference; delete the reference file itself.

### The run is slower than expected

Check `--machine` matches the hardware, and that `--modules` is not pulling in tracks the
user does not need. Fetching AlphaFold pLDDT across the proteome is inherently slow
(~8 h) — `--skip_alphafold true` with `--alphafold_precomputed_table` reuses a previous
run's scores.

## Step 4 — Report honestly

Say what failed, what the evidence was, and what is still missing. If a track cannot be
produced on this machine — a licensed predictor, an unreachable server, a missing genome
file — say so plainly and note what the run *did* produce, rather than leaving the user to
discover the gap later. A partially complete run is often perfectly usable, but only if
its gaps are stated.
