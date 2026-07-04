# Conservation Scores

## Description

Two complementary conservation signals are computed: GOPHER trident scores across seven taxonomic levels (protein-level, from multiple sequence alignments) and phastCons scores from 100-vertebrate whole-genome alignments (nucleotide-level, requires `combined_map.map`).

## Data sources

### GOPHER conservation

- **File:** `params.gopher_conservation_table` (no default — point it at your pre-computed GOPHER conservation table; GOPHER is skipped when unset)
- **Method:** GOPHER trident algorithm applied to seven taxonomic levels
- **Taxonomic levels:** global, Mammalia, Vertebrata, Eukaryota, Eumetazoa, Opisthokonta, Viridiplantae

### phastCons 100-vertebrate

- **Directory:** `params.phastcons_dir` (no default — a directory of per-chromosome `chr*.bw` files; phastCons is skipped when unset)
- **Format:** UCSC bigWig files, one per chromosome
- **Tool required:** `bigWigToBedGraph` in PATH (installed via conda UCSC tools)
- **Coordinate source:** Genomic positions from `combined_map.map` (requires Module 3 to have run)

## Output files

Both files are placed under `final/conservation/`:

| File | Description |
|------|-------------|
| `conservation_multiple_level.tsv` | GOPHER trident scores per residue across 7 taxonomic levels |
| `conservation_phastcons.tsv` | phastCons per-residue scores (averaged over 3 codon nucleotides) |

## Output columns

### conservation_multiple_level.tsv

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Position` | Residue position (1-based) |
| `Residue` | Amino acid |
| `global` | GOPHER trident score — all species |
| `Mammalia` | Score restricted to mammals |
| `Vertebrata` | Score restricted to vertebrates |
| `Eukaryota` | Score restricted to eukaryotes |
| `Eumetazoa` | Score restricted to eumetazoa |
| `Opisthokonta` | Score restricted to opisthokonts |
| `Viridiplantae` | Score restricted to plants |

### conservation_phastcons.tsv

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Position` | Residue position (1-based) |
| `phastCons` | Average phastCons score across the three codon positions (0–1) |
| `chrom` | Chromosome |
| `gpos_csv` | Genomic positions (comma-separated, from combined_map.map) |

## Notes

- phastCons scores are extracted by converting the bigWig for the relevant chromosome to BedGraph, then looking up each genomic position. The three nucleotide positions of a codon are averaged.
- GOPHER scores are looked up by UniProt accession and then mapped to each Gencode isoform via the transcript map.
- Both outputs skip gracefully if the input files are not provided (`--skip_gopher`, `--skip_phastcons` flags).
- Worker: `bin/create_conservation_worker.py`
