# UniProt Sequence Features — Regions of Interest & Binding Sites

## Description

Sequence features curated in the UniProt flat file for each entry, split into two
tracks:

- **Regions of interest (ROI)** — `REGION` features: functionally or
  structurally noteworthy stretches (e.g. "Interaction with X",
  "Disordered", "Necessary for …").
- **Binding sites** — `BINDING` features: residues/ranges that bind a ligand,
  cofactor, metal, or nucleotide.

## Data source

- **File:** the UniProt Swiss-Prot flat file (`uniprot_sprot.dat.gz`), parsed by
  `PARSE_UNIPROT_DAT`.
- **Origin:** [UniProt](https://www.uniprot.org/) feature (FT) lines.
- **Update policy:** Fetched from the current UniProt release; cached in `references/`.

## Output files

| File | Feature |
|------|---------|
| `final/annotations/uniprot_roi.tsv` | `REGION` features |
| `final/annotations/uniprot_binding.tsv` | `BINDING` features |

## Output columns

Both files share the same schema:

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `mapping_type` | `direct` / `homology_similarity` |
| `homology_transfer` | `True` when transferred to a different isoform by sequence match |
| `homology_identity` | Sequence identity of the transfer (homology rows only) |
| `Accession` | UniProt accession |
| `Type` | Feature type (`REGION` or `BINDING`) |
| `Start` | Feature start (1-based, UniProt coordinates) |
| `End` | Feature end (1-based, inclusive) |
| `Note` | Feature description (e.g. ligand name, region role) |
| `Evidence` | UniProt evidence code (ECO) |
| `Ligand` | Bound ligand (binding features only; blank for regions) |

## Notes

- Positions are in UniProt sequence space and converted to Gencode isoform
  coordinates by `TRANSCRIPT_MAP`.
- Worker: `bin/create_annotation_worker.py` (features extracted during
  `PARSE_UNIPROT_DAT`).
