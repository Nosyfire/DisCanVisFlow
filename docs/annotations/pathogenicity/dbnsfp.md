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

The raw single-file mode emits **~110 columns** (5 identity + 6 variant + all
predictor scores & rankscores + CADD + conservation + gnomAD AF). It is large
(tens of GB for the full proteome — one row per isoform residue × alternate allele).

## Output columns

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

Plus `CADD_raw`, `CADD_raw_rankscore`, `CADD_phred`.

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
