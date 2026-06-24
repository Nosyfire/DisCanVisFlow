# ScanSite 4.0 — Kinase and Phosphorylation Motifs

## Description

ScanSite 4.0 (MIT) predicts short sequence motifs that are likely to be phosphorylated by specific kinases or bound by specific domains. The pipeline uses pre-computed proteome-wide results at MAMMALIAN / High stringency by default.

## Modes

| Mode | When used | How |
|------|-----------|-----|
| **Pre-computed (default)** | `--scansite_tsv` points to a pre-processed TSV | Filter the file to proteins in this run |
| **Live API** | `--use_api` flag; `--scansite_tsv` is `NO_FILE` | Call `https://scansite4.mit.edu/webservice/proteinscan/` per protein |

Default stringency: **High** (MAMMALIAN motif class). Override with `--stringency Medium` or `--stringency Low`.

## Data source

- **Pre-computed file:** `params.scansite_tsv` — proteome-wide ScanSite 4.0 results, MAMMALIAN class, High stringency
- **Live API URL:** `https://scansite4.mit.edu/webservice/proteinscan/identifier={acc}/sequence={seq}/motifclass=MAMMALIAN/stringency={stringency}`
- **Update policy:** Replace the pre-computed TSV file to update. The API mode is always current but slow for large runs.

## Output file

`unmapped/annotations/scansite.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `motifName` | Full kinase/domain motif name (e.g. `CK2_PHOSPHO_1`) |
| `motifShortName` | Short motif identifier |
| `score` | ScanSite percentile score (lower = more significant match) |
| `site` | Phosphorylation or binding site position (1-based) |
| `siteSequence` | Matched peptide sequence (central residue is the site) |
| `Start` | Start position of the matched motif window (1-based) |
| `End` | End position of the matched motif window |

## Notes

- The `score` in ScanSite is a percentile: lower values indicate a better match. At High stringency the cutoff is typically 0.2% of all human sequences.
- The live API mode uses `ThreadPoolExecutor` with 10 workers by default (`--workers N` to change); it is slow for proteome-wide runs.
- In API mode, `Start` and `End` positions are derived by substring search of `siteSequence` in the full protein sequence.
- Position coordinates (`site`, `Start`, `End`) are in UniProt/input sequence space; they are mapped to Gencode isoform coordinates via `TRANSCRIPT_MAP`.
- Worker: `bin/create_scansite_worker.py`
