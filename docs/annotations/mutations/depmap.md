# DepMap Somatic Mutations

## Description

DepMap (Cancer Dependency Map) somatic mutation calls from cancer cell lines. The pipeline accepts either a normalized gene-keyed DepMap table produced by `bin/fetch_depmap_worker.py` or an older pre-processed Protein_ID-keyed table, then maps/filter it to the proteins in the current run.

## Data source

- **Normalized parameter:** `--depmap_tsv`, defaulting in the `discanvis_data` profile to `references/depmap/depmap_mutations_raw.tsv`
- **Raw-copy parameter:** `--depmap_raw_csv`, defaulting to `references/depmap/OmicsSomaticMutations.csv`
- **Origin:** DepMap portal, somatic mutation calls from CCLE/DepMap public releases
- **Manual download page:** <https://depmap.org/portal/download/all/>
- **Official file to choose:** `OmicsSomaticMutations.csv`
- **Update policy:** Manual refresh when DepMap changes the release.

The DepMap portal may put scripted downloads behind browser verification. For reproducible runs, download the official CSV manually before starting Nextflow:

```bash
# In a browser: download OmicsSomaticMutations.csv from:
# https://depmap.org/portal/download/all/
# Move/copy it to:
# references/depmap/OmicsSomaticMutations.csv
```

Then run the pipeline with `-resume`. At startup, the pipeline creates `references/depmap/` if needed and checks for either the normalized TSV or the raw CSV. If the normalized TSV exists, it prints a confirmation similar to:

```text
Manual reference confirmed: DepMap TSV -> /path/to/references/depmap/depmap_mutations_raw.tsv
```

If only the raw CSV exists, the pipeline prints a raw-file confirmation and runs a small normalization step before `DEPMAP_MAP`. If both files are missing, the pipeline stops before scheduling expensive tasks and prints a checklist of all missing manual references. To use a different location, pass `--depmap_raw_csv /path/to/OmicsSomaticMutations.csv` or `--depmap_tsv /path/to/depmap_mutations_raw.tsv`.

## Output file

`final/mutations/DepMap/depmap_mutations.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Chrom` | Chromosome |
| `Start_Position` | Genomic start (hg38) |
| `End_Position` | Genomic end (hg38) |
| `HugoSymbol` | HGNC gene symbol |
| `Protein_position` | Residue position in isoform |
| `HGVSp_Short` | Protein change in HGVS short notation (e.g. `p.V600E`) |
| `VariantType` | `SNP`, `INS`, `DEL`, etc. |
| `VariantInfo` | Functional class (Missense, Nonsense, etc.) |
| `DNAChange` | Nucleotide-level change |
| `ModelID` | DepMap cell line model ID (e.g. `ACH-000001`) |
| `Hotspot` | Boolean — whether the site is a known hotspot |
| `EntrezGeneID` | NCBI Entrez gene identifier |
| `Rescue` | Whether variant was rescued by a secondary event |
| `RescueReason` | Description of rescue event |

## Notes

- If DepMap is enabled and neither `--depmap_tsv` nor `--depmap_raw_csv` exists, the pipeline fails early with manual download instructions. Use `--skip_depmap true` to intentionally omit this track.
- Raw normalized DepMap mutations are gene-keyed and are expanded to matching isoforms by `create_depmap_worker.py`. Older Protein_ID-keyed preprocessed files are filtered directly to the current run.
- This annotation is distinct from ClinVar mutations: DepMap data is cancer cell line–specific, not clinical, and is not filtered by ClinVar significance.
- Worker: `bin/create_depmap_worker.py`
