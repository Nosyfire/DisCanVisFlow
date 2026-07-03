# Design — Report provenance, base statistics & documentation overhaul

_Date: 2026-07-03 · Status: DRAFT (awaiting user approval)_

## Problem statement

Three issues surfaced after the Zenodo packaging run:

1. **cellular_vulnerability was copied, not pipelined.** The full-proteome content
   was `rsync`ed from `results/discanvis` instead of running
   `nextflow run --project cellular_vulnerability -resume`. Consequence: the
   `MAPPING_REPORT` process never fired, so the summary Markdown files under
   `results/cellular_vulnerability/` are **stale** — they describe the earlier
   partial run, not the merged full-proteome content now on disk.

2. **`--modules disorder,mutations` appeared to run everything.** Investigation
   (static trace of ~30 `mods == null || mods.contains(...)` gates in `main.nf`
   + a `-stub` DAG check) shows **gating works correctly**. What the user saw is
   the combination of (a) the always-on backbone (ELM, Pfam, DIBS, MFIB,
   PhasePro, PTM, AlphaFold — never gated, by design) and (b) `-resume` leaving
   previous `publishDir` outputs on disk, so the results tree *looks* complete
   even though the gated processes never executed. **Not a code bug** — a UX /
   observability gap.

3. **Generated reports lack base provenance and statistics.** `versions.txt`
   captures software versions (Nextflow, pipeline, python, pandas, blastp, blat)
   but the reports have **no data-source provenance** (GENCODE release, UniProt
   release, download dates) and **no base counts** (# UniProt entries,
   # GENCODE entries, # direct 1:1 mappings, per-annotation totals).

## Scope & sequencing

Three semi-independent workstreams. Recommended order:

- **Phase 1 — Report enrichment + standalone regenerator** (this is the core).
- **Phase 2 — Modules/rerun observability** (small, no behavioural change).
- **Phase 3 — Documentation overhaul** (README + diagrams + examples).

Phase 1 is designed so it **also fixes Problem 1**: a standalone regenerator run
against `results/cellular_vulnerability/final/` produces fresh, enriched reports
in seconds — no 24 h pipeline rerun.

---

## Phase 1 — Report enrichment

### Key architectural decision: standalone regeneration

The user has multi-hour existing runs (`discanvis`, `vep_benchmarking`,
`cellular_vulnerability`) that must **not** be recomputed. Therefore the
provenance/stats logic must be runnable **standalone against an existing
`final/` directory**, not only inline in the Nextflow `MAPPING_REPORT` process.

Chosen approach (**Option 1 + standalone entry point**):

- Extend `bin/create_mapping_report_worker.py` with the new sections **and** a
  `--standalone` mode: given `--final_dir`, `--gencode_fasta`, `--uniprot_fasta`
  (+ optional `--intermediate_dir`), it computes provenance/counts directly and
  writes the enriched reports — no Nextflow, no config, no recompute.
- `MAPPING_REPORT` in the pipeline calls the same worker, now passing the extra
  provenance args it already has access to (config params + launchDir). Both
  paths share one code path; standalone just fills the same args from the CLI.

Rejected alternatives:
- *Dedicated upstream PROVENANCE process*: more plumbing, and can't backfill
  existing runs.
- *Report-time only (no standalone)*: forces a pipeline rerun to fix cv — the
  exact thing the user wants to avoid.

### New report content

**A. Data provenance & versions table** (new, near the top of `mapping_summary.md`):

| Source | Release / Version | Dated | File |
|---|---|---|---|
| GENCODE | v44 (parsed from filename) | file mtime | `gencode.v44.pc_translations.fa` |
| UniProt SwissProt | release from `reldate.txt` if present, else `.dat` mtime | download date | `uniprot_sprot.dat.gz` |
| ClinVar | `##fileDate` from VCF header | — | `clinvar.vcf.gz` |
| AlphaFold DB | v6 (from tar filename) | mtime | `UP000005640_9606_HUMAN_v6.tar` |
| dbNSFP | version from dir name (already parsed) | — | — |
| dbSNP / MobiDB / GO / MONDO / AlphaMissense | filename/config + mtime | — | — |

Version extraction rules:
- GENCODE: regex `gencode\.(v\d+)` on `params.gencode_fasta`.
- UniProt: `references/uniprot/reldate.txt` if present; else the `DT` release
  line is absent in `.dat` records, so fall back to file mtime as "download date".
- ClinVar: `zcat … | grep -m1 '##fileDate'`.
- Everything else: filename regex + mtime, driven by a small registry.

**B. Input scale / base counts** (new "Input scale" section):

- UniProt SwissProt canonical entries: N (`grep -c '^>'` on the fasta; ~20,586)
- UniProt curated isoforms (if `all_isoform_mapping`): N (additional fasta count)
- GENCODE protein-coding entries: N (`grep -c '^>'` on gencode translations fa)
- Genes in run: N
- Transcripts / isoforms processed: N (main + non-main)
- **Direct 1:1 UniProt↔GENCODE mappings: N** (from the seq table — rows where a
  transcript maps to exactly one UniProt isoform with `mapping_type=direct` /
  matching accession; source: `loc_chrom_with_names.tsv` + `bestmaps_blast_gene_transcript.tsv`)
- Homology-transferred mappings: N
- Genome-mapped isoforms: N (already computed)

**C. Per-annotation base stats** (enrich the existing coverage table):

Existing columns: main iso w/ data, non-main iso w/ data, ann (main), ann (non-main).
Add:
- Source rows (before mapping) — from `intermediate/`
- Mapped rows (after mapping) — from `final/`
- Mapping success rate = mapped/source
- Proteome coverage % = isoforms-with-data / total isoforms

### Where the data is sourced (no recompute)

All derivable from files already on disk:
- Counts: `grep -c '^>'` on the UniProt/GENCODE fastas (paths from config).
- Direct mappings: `bestmaps_blast_gene_transcript.tsv` + seq table.
- Versions: filenames, `.dat`/tar mtimes, VCF `##fileDate`.

### Testing

- Unit tests in `tests/test_create_mapping_report_worker.py`: feed a tiny fixture
  `final/` tree + fake fastas → assert the provenance table, input-scale counts,
  and per-annotation stats render with correct numbers.
- Regenerate `results/discanvis` report standalone and eyeball the new sections
  against known truths (20,586 UniProt entries, GENCODE v44).

---

## Phase 2 — Modules / rerun observability (small)

No behavioural change to gating. Add clarity:

1. **Explicit log block** at workflow start when `--modules` is set:
   ```
   ▶ Modules requested : disorder, mutations
   ▶ Always-on backbone: elm, pfam, dibs, mfib, phasepro, ptm, alphafold
   ▶ NOT running       : go, ppi, pdb, polymorphism, scansite, coiledcoils, …
   ⚠ -resume note: previously published outputs for NOT-running modules remain
     on disk under results/<project>/final/ and are NOT regenerated this run.
   ```
2. Document the `-resume`/`publishDir` stale-file behaviour in README + CLAUDE.md.
3. *(Deferred / YAGNI)* a `--clean_publish` flag to wipe stale outputs — only if
   the user wants it; not in this phase.

---

## Phase 3 — Documentation overhaul

- **README restructure**: Quick-start → Concepts → Config axes → Modules →
  Provenance/Data sources (with the version table) → Troubleshooting.
- **Diagrams**: replace the ASCII DAG with a **Mermaid** flowchart (renders on
  GitHub) for the module map, plus a generated real DAG via
  `nextflow run … -with-dag docs/dag.mmd` committed as an artifact.
- **Data-sources section**: single canonical table of every source, its release,
  licence, and how it's supplied — mirrors the new report provenance table.
- **Worked examples**: verify and include the CLAUDE.md command examples as a
  copy-pasteable "recipes" section.

---

## Deliverables

1. `bin/create_mapping_report_worker.py` — new provenance/scale/per-annotation
   sections.
2. `bin/regenerate_reports.sh` (or `--standalone` mode) — backfill existing runs.
3. `tests/test_create_mapping_report_worker.py` — coverage for new sections.
4. `modules/annotation_mapping.nf` (`MAPPING_REPORT`) — pass new args.
5. `main.nf` — modules/backbone/NOT-running log block.
6. README + `docs/` — restructure, Mermaid diagrams, data-sources table, recipes.
7. Regenerated enriched reports for `discanvis`, `cellular_vulnerability`,
   `vep_benchmarking` (standalone, no recompute) — then re-archive for Zenodo.

## Out of scope

- Re-running any full-proteome pipeline for compute (only report regeneration).
- Changing module gating logic (verified correct).
- New annotation modules.
