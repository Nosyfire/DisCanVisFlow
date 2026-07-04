# GO — Gene Ontology Terms

## Description

Gene Ontology annotations for each isoform: the biological processes, molecular
functions, and cellular components associated with the protein, drawn from the
GOA (GO Annotation) human dataset and expanded with ontology metadata (name,
namespace, definition, parents) from the GO OBO file.

## Data source

- **Annotations:** `goa_human.gaf.gz` (GO Annotation file, keyed by UniProt accession).
- **Ontology:** `go.obo` (term definitions and hierarchy).
- **Fetch:** `FETCH_GO`, cached via `storeDir` in `references/`.
- **Origin:** [Gene Ontology](http://geneontology.org/).
- **Update policy:** Always-current — refresh with `bin/refresh_refs.sh go`.

## Output file

`final/annotations/go_terms.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Entry_Isoform` | UniProt accession |
| `GO_Term` | GO identifier (e.g. `GO:0004672`) |
| `name` | Term name (e.g. `protein kinase activity`) |
| `namespace` | `molecular_function` / `biological_process` / `cellular_component` |
| `def` | Term definition text |
| `alt_id` | Alternative GO IDs for the term |
| `is_a` | Parent term(s) in the ontology hierarchy |

## Notes

- Annotations are looked up by UniProt accession, so all Gencode isoforms of a
  gene inherit the same GO terms.
- Worker: `bin/create_go_worker.py` (Module 5f).
