# Pfam Domains

## Description

Pfam domains are the protein family / domain assignments for each isoform —
the conserved structural and functional modules (kinase domains, SH2/SH3,
RRM, etc.). They are drawn from the InterPro → Pfam mapping keyed by UniProt
accession, with per-domain start/end envelope coordinates.

## Data source

- **File:** `protein2ipr.dat.gz` from InterPro (filtered to Pfam signatures during
  `PARSE_UNIPROT_DAT`), plus the UniProt flat file for accession → sequence.
- **Origin:** [Pfam](http://pfam.xfam.org/) / [InterPro](https://www.ebi.ac.uk/interpro/)
- **Update policy:** Fetched and parsed from the current InterPro release; cached in
  `references/`.

## Output file

`final/annotations/pfam_domains.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Accession` | UniProt accession |
| `hmm_acc` | Pfam accession (e.g. `PF00069`) |
| `hmm_name` | Pfam family name (e.g. `Pkinase`) |
| `type` | Pfam entry type (e.g. `Domain`, `Family`, `Repeat`) |
| `envelope_start` | Domain envelope start (1-based) |
| `envelope_end` | Domain envelope end (1-based, inclusive) |

## Notes

- Pfam domains are used by the [Combined disorder](../disorder/disorder.md) track
  to *exclude* structured domain regions from disorder calls.
- Coordinates are in the isoform's sequence space after transcript mapping.
- Parsing worker: `bin/parse_uniprot_dat_worker.py`; the InterPro column indices
  are Pfam-specific (see [Troubleshooting](../../guide/troubleshooting.md#pfam_domainstsv-is-empty)).
