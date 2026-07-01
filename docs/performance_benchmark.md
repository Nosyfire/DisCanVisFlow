# Performance Benchmark — DisCanVisFlow

**Run configuration**

| Field | Value |
|-------|-------|
| Gene | RAF1 (single-gene mode, `all_isoform_mapping`) |
| Project | `test_one_protein` |
| Data | `local` (all references pre-existing, zero downloads) |
| Machine | `hard` — 64 CPUs, 256 GB RAM cap, queueSize=24 |
| Server | `gpu0.dlab.elte.hu` — 128 physical cores, 1.8 TB RAM |
| Nextflow | 26.04.0 |
| Run name | `extravagant_kare` |
| Total wall time | **4m 8s** (63 tasks, 0 cached, 0 failed) |
| Trace | `results/test_one_protein/reports/trace.tsv` |

---

## Per-Process Timing (sorted by realtime)

| Process | Realtime | CPUs | Peak RSS | I/O read |
|---------|----------|------|----------|----------|
| POLYMORPHISM_MAP | **2m 34s** | 16 | 191 MB | 390 MB |
| DBNSFP_MAP | **2m 11s** | 24 | 81 MB | 2.1 GB |
| PDB_MAP | 44.3s | 16 | 87 MB | 80 MB |
| MUTATION_MAP_TCGA | 42.9s | 16 | 544 MB | 3.5 GB |
| BLAT_ALIGN ×16 chunks (parallel) | 38s wall* | 1 each | 4 GB each | 797 MB each |
| MUTATION_MAP_CLINVAR | 32.4s | 16 | **5.1 GB** | 205 MB |
| MUTATION_MAP_CBIOPORTAL | 26.0s | 16 | 211 MB | 1.6 GB |
| PPI_MAP | 22.7s | 2 | 215 MB | 70 MB |
| CONSERVATION_MAP | 9.1s | 16 | **1.7 GB** | 1.0 GB |
| SEQUENCE_PROCESS | 6.8s | 16 | 310 MB | 76 MB |
| ALPHAMISSENSE_MAP | 6.3s | 16 | 80 MB | **9.4 GB** |
| ANNOTATION_MAP | 3.5s | 16 | 460 MB | 269 MB |
| MAVEDB_MAP | 3.3s | 16 | 543 MB | 320 MB |
| CLINVAR_DISEASE_BUILD | 5.1s | 2 | 244 MB | 76 MB |
| All other tasks | < 2s each | — | < 50 MB | < 55 MB |

*16 BLAT chunks run in parallel on this server; longest chunk = 38s → 38s wall time.

---

## Critical Path

```
SUBSET (3.5s)
  → MAKEBLASTDB + BLASTP (7s)
  → MERGE/ID_MAP (5s)
  → SEQUENCE_PROCESS (8.9s)
    ├── [parallel annotation wave: MOBIDB, GO, PEM, SCANSITE, …]  ~5s
    ├── MAVEDB_MAP (3.3s) → ANNOTATION_MAP (3.5s)
    │     → PDB_MAP (44.3s)   ← bottleneck of annotation phase
    └── [BLAT_ALIGN ×16 parallel] (38s wall)
          → GENOME_MAP (0.5s)
  → MUTATION_MAP_CLINVAR (32.4s)
  → MUTATION_MAP_TCGA (44.9s)  ← mutation-phase bottleneck
  → CONSERVATION_MAP (9.1s)
  → DBNSFP_MAP (2m 11s)  ← serialised after GENOME_MAP
  → POLYMORPHISM_MAP (2m 34s)  ← longest single task
  → MAPPING_REPORT (0.9s)
```

The last two tasks on the critical path (**POLYMORPHISM_MAP** + **DBNSFP_MAP**) account for ~**4m 45s** of compute time on their own, yet the total run is only 4m 8s because they overlap in execution with several other late-stage tasks (TCGA, conservation). Eliminating or accelerating these two would cut wall time by roughly 50 % even for a single gene. For the full proteome they scale with the number of transcripts and become the dominant cost.

---

## Top Bottlenecks & Suggested Improvements

### 1. POLYMORPHISM_MAP — 2m 34s (single-gene)

**Root cause**: `bigBedToBed` is called once per isoform to extract SNPs from the full-genome `dbSnp155Common.bb` file. Each call loads nearly the entire bigBed index. With ~10 isoforms for RAF1 this is ~15s × 10 = serial; it scales to **hours** for the full proteome (≥20 000 isoforms).

**Fix — pre-extract all intervals in one pass**:
1. Add a `POLYMORPHISM_PREFETCH` process immediately after `GENOME_MAP`.
2. This process collects all isoform genomic intervals from `combined_map.map` into a single BED file, runs **one** `bigBedToBed -bed combined_intervals.bed` call, and emits the result.
3. `POLYMORPHISM_MAP` per-isoform workers filter from this pre-extracted BED rather than calling `bigBedToBed` individually.

**Expected speedup**: 10-100× depending on number of isoforms. Full-proteome runtime: hours → minutes.

---

### 2. DBNSFP_MAP — 2m 11s (single-gene); scales to hours for full proteome

**Root cause**: Scans all `chr*.gz` raw dbNSFP files (24 chromosomes, each ~10 GB compressed) looking for positions in `combined_map.map`. Each scatter chunk reads the entire directory.

**Fix A — use pre-mapped TSV** (already supported):
Switch `--dbnsfp_tsv` instead of `--dbnsfp_raw_dir`. This uses a pre-mapped Protein_ID-keyed TSV that filters to the target genes in seconds. Build the full-proteome TSV once and reuse across all project runs.

**Fix B — regional index** (for always-fresh approach):
Pre-extract the gene's chromosomal region from each chr*.gz using `tabix` or a region-filter script before the scatter phase. Feed the regional extract to each scatter chunk.

**Expected speedup**: Fix A is essentially instant (seconds) for any gene set.

---

### 3. PDB_MAP — 44s (single-gene); N × PDBe REST API calls

**Root cause**: Queries the PDBe REST API per UniProt isoform. Latency accumulates linearly.

**Fix — use SIFTS bulk download** (already implemented):
Pass `--pdb_bulk true`. This downloads the SIFTS `uniprot_segments_observed.tsv.gz` once (cached in `references/sifts/`), then maps all isoforms locally — no per-isoform HTTP calls.

`--pdb_bulk true` is already set in `discanvis_data.config`. It should also be the default in `local_refs.config` for any run that has internet access.

**Expected speedup**: REST latency eliminated; PDB_MAP drops to < 2s.

---

### 4. MUTATION_MAP_TCGA — 42.9s; reads 3.5 GB per scatter chunk

**Root cause**: The full TCGA MAF (3.5 GB) is loaded from disk for every scatter chunk. With `scatter_chunks=20` this means 20 × 3.5 GB = 70 GB of sequential I/O per full-proteome run.

**Fix — pre-filter MAF by chromosome / genomic region**:
Add a `MUTATION_PREFILTER` step that slices the MAF to only positions in the scatter chunk's genomic regions before dispatching to the worker.

Alternatively: index the MAF by chromosome using `tabix` and pass region-specific slices to each chunk.

**Expected speedup for full proteome**: 20×+ reduction in per-chunk I/O.

---

### 5. MUTATION_MAP_CLINVAR — 32.4s; 5.1 GB peak RSS

**Root cause**: Loads the full ClinVar VCF (205 MB on disk → 5.1 GB inflated in memory as a Python data structure) for every scatter chunk, then scans for matching positions.

**Fix — pre-filter with bcftools**:
Before `MUTATION_MAP`, run a `CLINVAR_REGION_EXTRACT` step using `bcftools view -R regions.bed clinvar.vcf.gz` to extract only positions relevant to the target genes/scatter chunk. Feed the small VCF to each worker.

**Expected speedup**: 5-10× per-chunk time reduction; peak RSS drops from 5.1 GB to < 100 MB.

---

### 6. ALPHAMISSENSE_MAP — reads 9.4 GB for a single gene

**Root cause**: The AlphaMissense TSV (`AlphaMissense_isoforms_hg38.tsv.gz`) is 2.6 GB compressed and covers all human transcripts. The worker reads it entirely to find RAF1 entries.

**Fix — pre-index by gene**:
Index the AlphaMissense file by gene (e.g. with `tabix` on a bgzip-compressed version keyed by HugoSymbol) so each worker only reads a kilobyte-range slice. Or build a per-gene lookup cache in `references/alphamissense/<gene>.tsv` on first access.

**Expected speedup**: I/O drops from 9.4 GB to < 1 MB per gene.

---

### 7. PPI_MAP — 22.7s on 2 CPUs (underallocated)

**Root cause**: PPI mapping is purely pandas-based (join + filter), but it is configured with `cpus = 2`. It holds 215 MB in memory. Increasing to 8 CPUs won't help pandas directly, but raising to `process_medium` allocation (16 CPUs) and using `--nthreads` inside the worker could.

**Fix**: Increase `PPI_MAP` label from `process_low` to `process_medium` in `nextflow.config`. Separately, rewrite the inner join to use DuckDB or polars for 3-5× speedup without needing extra CPUs.

---

### 8. BLAT_ALIGN — 38s wall time for RAF1 (already well-parallelized)

16 chunks run in parallel, each using 1 CPU and 4 GB RAM. For RAF1 this is good: 16 chunks × 1s/chunk on 16 cores = 38s wall.

For the full proteome, the GENCODE cDNA is ~170 MB. With `blat_chunks=16` each chunk is ~10 MB, completing in ~40s. With `blat_chunks=64` (and 64 CPUs) the wall time would drop to ~10s. BLAT is not the critical path bottleneck at full proteome scale.

**Suggestion**: Leave at 16 chunks for single-gene; increase `blat_chunks = 64` in `hard.config` for full-proteome runs.

---

## Summary Prioritization

| Priority | Bottleneck | Status | Estimated Full-Proteome Savings |
|----------|-----------|--------|-------------------------------|
| 🔴 Critical | POLYMORPHISM_MAP: chromosome-sweep bigBedToBed | ✅ **Implemented** | N calls → ≤24 calls (hours → minutes) |
| 🔴 Critical | DBNSFP_MAP: switch to pre-mapped TSV (`--dbnsfp_tsv`) | ⚠️ User choice | hours → seconds |
| 🟠 High | PDB_MAP: default to `pdb_bulk = true` | ✅ **Implemented** in `local.config` and `discanvis_data.config` | 44s → 2s |
| 🟠 High | AlphaMissense: per-gene pre-index | ❌ Pending | 9.4 GB/gene read → 1 MB/gene |
| 🟡 Medium | MUTATION_MAP_CLINVAR: bcftools region filter | ❌ Pending | 5.1 GB RSS → 100 MB |
| 🟡 Medium | MUTATION_MAP_TCGA: pre-filter MAF by region | ❌ Pending | 3.5 GB/chunk I/O → 10 MB/chunk |
| 🟢 Low | PPI_MAP: raise CPU allocation; explore polars | ❌ Pending | 22s → 5s |
| 🟢 Low | BLAT: increase to 64 chunks in hard.config | ❌ Pending | 38s → 10s |

### Implemented changes (this session)

**POLYMORPHISM_MAP** (`bin/create_polymorphism_worker.py`): replaced per-isoform `bigBedToBed` calls
with chromosome-level sweeps. Instead of one call per unique `(chrom, start, end)` region (= one per
isoform), the worker now computes the bounding box of all isoforms per chromosome and calls
`bigBedToBed` once per chromosome (at most 24 calls for the full human genome). The returned SNPs are
then filtered to each isoform's exact interval in memory using a range check before the `g2p` dict
lookup.

**PDB** (`config/data/local.config`, `config/data/discanvis_data.config`): `pdb_bulk = true` is now
the default in both data configs. This routes PDB annotation through the SIFTS bulk join path
(downloads `uniprot_segments_observed.tsv.gz` once) instead of per-isoform PDBe REST API calls.

**DBNSFP_MAP** for full proteome: raw-dir scan is still the default for fresh runs. For full-proteome
jobs, supply a pre-mapped TSV via `--dbnsfp_tsv` to bypass the 2+ minute per-scatter-chunk file scan.

The two critical implemented items together should reduce POLYMORPHISM_MAP from >2 min to seconds on
full-proteome runs, and PDB_MAP from 44s to ~2s.
