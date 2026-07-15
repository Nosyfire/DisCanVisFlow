# Pathogenicity Scores (dbNSFP)

## Description

Per-variant pathogenicity, conservation, and population-frequency annotations
from **dbNSFP** (database for Nonsynonymous SNPs' Functional Predictions). Every
possible missense variant in the coding genome is annotated with the raw scores
and rank scores of **37 pathogenicity predictors**, CADD, several conservation
tracks, and gnomAD allele frequencies. The pipeline maps these variant-level
scores onto every curated isoform's residues via `combined_map.map`.

## Data source

dbNSFP has **two mutually-exclusive input modes** (raw takes priority):

1. **Raw single-file mode (primary) — `--dbnsfp_raw_dir`**
   Point at the merged dbNSFP academic release, e.g.
   `/dlab/home/norbi/data/dbNFSP/dbNSFP5.3.1a_grch38.gz` (a single ~50 GB
   bgzip-free gzip, 505 columns, GRCh38). The worker detects a self-describing
   merged `.gz`, builds an inverted `(chr, genomic_pos) → [(isoform, residue,
   aa)]` index from `combined_map.map`, and streams the file **once** (`pigz -dc`)
   with `aaref` validation against the mapped residue. A directory of legacy
   per-chromosome `chr*.gz` files (dbNSFP 4.x) is still supported via the same
   parameter + optional `--dbnsfp_bed_header`.
   → `DBNSFP_MAP` / `bin/create_dbnsfp_map_worker.py`

2. **Pre-mapped mode — `--dbnsfp_tsv`**
   A static, already Protein_ID-keyed TSV
   (`dbNSFP_custom/mapped_filtered_mutations.tsv`).
   → `PATHOGENICITY_MAP` / `bin/create_pathogenicity_worker.py`

- **Origin:** dbNSFP (https://sites.google.com/site/jpopgen/dbNSFP), Liu et al.
- **Update policy:** Static release file; supply a newer dbNSFP download to update.

## Output file

| Mode | File |
|------|------|
| Raw (`--dbnsfp_raw_dir`) | `final/pathogenicity/dbnsfp_scores.tsv` |
| Pre-mapped (`--dbnsfp_tsv`) | `final/pathogenicity/pathogenicity_scores.tsv` |

The raw single-file mode emits **110 columns** (2 identity + 3 genomic + 6
variant + 81 predictor score/rankscore + 16 conservation + 2 gnomAD AF), one row
per isoform residue × alternate allele. For the full proteome that is
**~100M rows / ~173 GiB uncompressed**, so it is delivered as a compressed,
randomly-sliceable bundle (see below) rather than a plain TSV.

### Where the processed output lives (for benchmarking)

The full-proteome (`--project discanvis`) processed dbNSFP is a **BGZF bundle**
under `results/discanvis/final/pathogenicity/`, sharing the `dbnsfp_scores` stem:

```
results/discanvis/final/pathogenicity/dbnsfp_scores.tsv.gz       # BGZF body (~12 GB), sorted by Protein_ID,Protein_position
results/discanvis/final/pathogenicity/dbnsfp_scores.tsv.gz.gzi   # bgzip random-access index
results/discanvis/final/pathogenicity/dbnsfp_scores.pidx         # Protein_ID <tab> offset <tab> length
results/discanvis/final/pathogenicity/dbnsfp_scores.header       # the single 110-column header line
```

Consume it in a benchmark **without decompressing the whole file** — slice the
isoform(s) you need (a slice reads a few KB):

```bash
# one isoform → stdout (with header)
bin/slice_dbnsfp.py --bgz results/discanvis/final/pathogenicity/dbnsfp_scores.tsv.gz --id RAF1-201

# a benchmark gene set → one TSV
bin/slice_dbnsfp.py --bgz results/discanvis/final/pathogenicity/dbnsfp_scores.tsv.gz \
    --id_file benchmark_genes_transcripts.txt --out bench_dbnsfp.tsv
```

`--id` takes GENCODE transcript names (`Protein_ID`, e.g. `RAF1-201`), the same
primary key used across every other `final/` table, so a benchmark row keyed by
`(Protein_ID, Protein_position)` joins directly against this file. The whole
uncompressed body can be streamed with `bgzip -dc dbnsfp_scores.tsv.gz` if a tool
needs the plain TSV.

### Mapping semantics (important)

Each isoform is mapped **directly** from its own genomic coordinates in
`combined_map.map` — there is **no homology/sequence transfer** between isoforms
in the merged path. combined_map already contains every curated isoform
independently genome-mapped, so a variant is attributed only to isoforms whose
own codon sits at that genomic position (a residue therefore spans **at most 3
distinct `Start_Position` values** — its 3 codon positions).

**On row counts / duplication.** A clean residue maps to 3 codon positions × 3
alternate alleles = up to 9 SNV rows. In practice a residue shows ~8–12 rows.
The cross-isoform coordinate inflation bug (a residue picking up dozens of
foreign `Start_Position`s via homology transfer) was fixed in `f627e27` — every
residue now correctly has ≤ 3 distinct `Start_Position`s (verified: `A1BG-201`
residue 469 = 12 rows / 3 positions; `RAF1-201` residues likewise). The residual
> 9 rows are **not a pipeline bug**: dbNSFP itself lists the *same* genomic SNV
under several transcript reading-frame annotations (e.g. the same `C>A` at one
position appears once with `aaref/aaalt = ./.` and once as `R>S`). Rows are kept
whenever the reference amino acid matches (or is `.`), so each is a legitimate
per-variant record carried verbatim from dbNSFP. If a benchmark wants exactly one
row per `(Start_Position, ref, alt)`, deduplicate on those three columns keeping
the row whose `aaref` equals the residue's amino acid.

### Packed, sliceable bundle

For the full proteome the raw output is post-processed by `bin/dbnsfp_pack.sh`
into four files sharing the `dbnsfp_scores` stem in `final/pathogenicity/`:

| File | Contents |
|------|----------|
| `dbnsfp_scores.tsv.gz` | BGZF (block-gzip) body, sorted by `Protein_ID` then `Protein_position` (~12 GB) |
| `dbnsfp_scores.tsv.gz.gzi` | bgzip random-access index |
| `dbnsfp_scores.pidx` | `Protein_ID  offset  length` — byte range of each isoform |
| `dbnsfp_scores.header` | the single header line |

Slice one (or several) isoforms without reading the whole file:

```bash
bin/slice_dbnsfp.py --bgz final/pathogenicity/dbnsfp_scores.tsv.gz --id RAF1-201
bin/slice_dbnsfp.py --bgz .../dbnsfp_scores.tsv.gz --id RAF1-201,BRAF-201 --out out.tsv
```

It looks up the `.pidx` and asks `bgzip -b <offset> -s <length>` to decompress
only that isoform's bytes (a few KB). Requires `bgzip` (htslib) on PATH.

## Output columns

> **Machine-readable column reference:** the full ordered
> `index → column → category` map is in
> [`dbnsfp_columns.tsv`](dbnsfp_columns.tsv) (repo-tracked), and the exact header
> line ships next to the data as
> `results/discanvis/final/pathogenicity/dbnsfp_scores.header`. Agents/tools should
> read one of those two rather than hard-coding column order.


### Variant identifiers

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Protein_position` | Residue position in isoform |
| `chr` | Chromosome (hg38) |
| `Start_Position` / `End_Position` | Genomic coordinates (hg38) |
| `ref` / `alt` | Reference / alternate nucleotide |
| `aaref` / `aaalt` | Reference / alternate amino acid |
| `aapos` | dbNSFP amino-acid position (per-transcript, `;`-separated) |
| `rs_dbSNP` | dbSNP rsID |

### Predictor scores (37 methods)

Each method contributes a raw `*_score` and a normalized `*_rankscore`
(0–1, higher = more deleterious), preserved verbatim from dbNSFP (per-transcript
values are `;`-separated; missing = `.`):

SIFT, SIFT4G, Polyphen2_HDIV, Polyphen2_HVAR, MutationTaster, MutationAssessor,
PROVEAN, VEST4, MetaSVM, MetaLR, MetaRNN, M-CAP, REVEL, MutPred2, MVP, gMVP,
MisFit_D, MisFit_S, MPC, PrimateAI, DEOGEN2, BayesDel_addAF, BayesDel_noAF,
ClinPred, LIST-S2, VARITY_R, VARITY_ER, VARITY_R_LOO, VARITY_ER_LOO, ESM1b,
AlphaMissense, PHACTboost, MutFormer, MutScore, popEVE, DANN, fathmm-XF_coding.

Plus `CADD_raw`, `CADD_raw_rankscore`, `CADD_phred`,
`Eigen-raw_coding_rankscore`, `Eigen-PC-raw_coding_rankscore`, and
`bStatistic_converted_rankscore` (rankscore-only columns from dbNSFP).

### Conservation

`GERP++_NR`, `GERP++_RS` (+rankscore), `GERP_92_mammals` (+rankscore),
`phyloP100way_vertebrate` / `phyloP470way_mammalian` / `phyloP17way_primate`
(+rankscores), `phastCons100way_vertebrate` / `phastCons470way_mammalian` /
`phastCons17way_primate` (+rankscores).

### Population allele frequency (gnomAD)

| Column | Description |
|--------|-------------|
| `gnomAD4.1_joint_AF` | gnomAD v4.1 joint (exomes+genomes) allele frequency |
| `gnomAD4.1_joint_POPMAX_AF` | gnomAD v4.1 joint popmax allele frequency |

gnomAD AF is taken **directly from dbNSFP** — the separate gnomAD server fetch
(`--fetch_gnomad_vcf` / `--gnomad_maf`) is not needed for these frequencies.

## Notes

- Missing scores appear as `.` (dbNSFP convention) or empty cells.
- SIFT/SIFT4G scores are inverted relative to pathogenicity (lower = more
  deleterious); use rankscores for cross-method comparison.
- CADD phred ≥ 20 ≈ top 1% most deleterious genome-wide.
- AlphaMissense also appears standalone in `alphamissense.tsv` (Module 8d),
  GENCODE isoform–keyed; here it follows the dbNSFP genomic mapping.
- Column selection lives in `select_keep_columns()`
  (`bin/create_dbnsfp_map_worker.py`) — a pattern rule over `_score` /
  `_rankscore` / CADD / GERP·phyloP·phastCons / gnomAD4.1 joint AF.
