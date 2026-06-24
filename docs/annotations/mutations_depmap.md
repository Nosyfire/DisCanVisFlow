# DepMap Somatic Mutations

## Description

DepMap (Cancer Dependency Map) somatic mutation calls from cancer cell lines. The pipeline filters a pre-processed, Protein_ID-keyed DepMap mutation table to the proteins in the current run.

## Data source

- **Parameter:** `--depmap_tsv` pointing to a pre-processed file (format: `mapped_filtered_mutations.tsv`)
- **Origin:** DepMap portal (depmap.org), somatic mutation calls from CCLE/DepMap public release
- **Pre-processing:** The raw DepMap mutation MAF is processed externally (outside this pipeline) to produce a Protein_ID-keyed TSV aligned to GENCODE isoforms.
- **Update policy:** Static pre-processed file; supply a new file to update.

## Output file

`unmapped/mutations/depmap_mutations.tsv` (also documented historically as `unmapped/annotations/` — use `mutations/` path)

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

- If `--depmap_tsv` is not provided or the file is empty/missing, the output is an empty TSV with the header columns only. The pipeline does not fail.
- DepMap mutations are kept at the isoform level (Protein_ID-keyed) and are not further expanded to additional isoforms.
- This annotation is distinct from ClinVar mutations: DepMap data is cancer cell line–specific, not clinical, and is not filtered by ClinVar significance.
- Worker: `bin/create_depmap_worker.py`
