# Performance & Scaling

How long a run takes, which stages dominate, and the knobs that matter. Numbers
below are order-of-magnitude guidance from a 64-core / 256 GB server; your
absolute times depend on hardware, gene count, and which modules you enable.

For the exact timings of *your* run, read the Nextflow execution trace at
`results/<project>/reports/trace.tsv` and the per-process table in
`results/<project>/mapping_reports/mapping_summary.md`.

---

## Typical runtimes

| Scope | All tracks | Notes |
|-------|-----------|-------|
| Single gene (e.g. RAF1) | ~4 min | dominated by genome/mutation/pathogenicity mapping |
| Small gene list (10–100) | minutes–tens of minutes | scales roughly linearly with transcript count |
| Full proteome (~20k genes) | ~24 h | see stage breakdown below |

Full-proteome stage costs (all tracks, 20 scatter chunks):

| Stage | Approx. time | Why |
|-------|-------------|-----|
| BLAST (reciprocal, full proteome) | ~2 h | one-time; cached and reused across projects |
| `DISORDER_MAP` (IUPred + AIUPred + AlphaFold) | ~8 h | model inference per isoform + AlphaFold API fetch |
| `DBNSFP_MAP` (raw `chr*.gz`) | ~3 h | scans large per-chromosome files |
| Everything else | remainder | mostly sub-minute per gene, parallelized |

---

## What scales with proteome size

The genome-anchored and large-reference stages grow with the number of
transcripts and are where full-proteome time goes:

- **`DISORDER_MAP`** — per-isoform model inference plus AlphaFold pLDDT fetch.
- **`DBNSFP_MAP`** — scans raw dbNSFP `chr*.gz` unless a pre-mapped TSV is used.
- **`MUTATION_MAP`** — scans the mutation source (ClinVar VCF / MAF) against the
  genome map; peak RAM ~5 GB for ClinVar.
- **`POLYMORPHISM_MAP`** — extracts dbSNP intervals from the bigBed file.

BLAST, `GENOME_MAP`, and all `FETCH_*` downloads are computed once and reused by
every project that shares the same `--data` cache, so a second project is much
faster than the first.

---

## Tuning knobs

| Knob | Effect |
|------|--------|
| `--machine <laptop\|low\|medium\|hard\|slurm>` | Sets CPUs, RAM caps, and parallelism (`blat_chunks`, `scatter_chunks`, queue size) to match your hardware — see [Configuration guide § `--machine`](configuration.md#--machine--how-much-cpuram-to-claim) |
| `scatter_chunks` | Splits the sequence table into N gene-balanced chunks so `DISORDER_MAP`, `DBNSFP_MAP`, and `COILEDCOILS_MAP` run concurrently. Full proteome: 20 |
| `blat_chunks` | Parallel BLAT jobs; each loads the ~4 GB hg38.2bit, so cap at your CPU/RAM budget. `test_one_protein` forces 1 |
| `pdb_bulk = true` | Maps PDB coverage from one SIFTS download instead of per-protein PDBe API calls (~10 min vs ~9 h at full proteome). Default in both data configs |
| `--dbnsfp_tsv` | Uses a pre-mapped, Protein_ID-keyed dbNSFP TSV instead of scanning raw `chr*.gz` — seconds instead of hours. Build it once, reuse everywhere |
| `--alphafold_precomputed_table` | Reuses AlphaFold pLDDT scores from a prior run so `DISORDER_MAP` skips the ~8 h EBI fetch while still recomputing IUPred/AIUPred |

---

## Memory planning

The two processes worth watching on memory-constrained machines:

- **`BLAT_ALIGN`** — loads the full hg38.2bit (~4 GB) per chunk. Budget ~6 GB
  per parallel chunk.
- **`MUTATION_MAP_CLINVAR`** — ~5 GB peak while scanning the ClinVar VCF.

With `blat_chunks = 1` (enforced by `test_one_protein`) the peak footprint of a
single-gene run is ~6 GB. A machine with ≥ 16 GB free RAM comfortably runs
`--machine medium`. See
[Configuration guide § Per-process resource usage](configuration.md#per-process-resource-usage-raf1-single-gene---modules-mutationsdisorder)
for the detailed per-process table.
