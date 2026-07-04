# Exon Boundaries

## Description

The exon structure of each isoform mapped onto protein coordinates: which
residue range each exon encodes, and the genomic span of that exon. Useful for
relating protein features (domains, motifs, mutations) to exon architecture and
splice boundaries.

## Data source

- **Input:** `combined_map.map` (per-residue protein ↔ genomic coordinate map from
  `GENOME_MAP`) and the GENCODE GTF exon structure.
- **Requires:** genome mapping (`params.hg38_2bit` set) — this is a genome-anchored
  track.
- **Update policy:** Recomputed from the genome map each run.

## Output file

`final/genome/exon.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `exon_number` | Exon index within the transcript (1-based) |
| `total_exons` | Number of coding exons in the transcript |
| `aa_start` | First protein residue encoded by this exon |
| `aa_end` | Last protein residue encoded by this exon |
| `aa_length` | Residues encoded by this exon |
| `genomic_start` | Exon start (hg38) |
| `genomic_end` | Exon end (hg38) |

## Notes

- Exon boundaries in protein space can split a codon; `aa_start`/`aa_end` mark the
  residues wholly or partly encoded by the exon.
- Worker: `bin/create_exon_worker.py` (Module 5d), downstream of `GENOME_MAP`.
