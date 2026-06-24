# Annotation Overview — DisCanVis Pipeline

This index covers every annotation type produced by the pipeline. Annotations are divided into two groups based on whether they require genomic coordinates (from `combined_map.map`) or work directly on protein sequences and accessions.

---

## Genome-level annotations

These annotations require genomic coordinates derived from Module 3 (Genome Mapping). The pipeline must have been run with `params.hg38_2bit` set. Results are placed under `results/{gene_dir}/mutations/` or `results/{gene_dir}/unmapped/conservation/`.

| Annotation | Module | Output file(s) | Worker |
|------------|--------|----------------|--------|
| ClinVar disease variants | 4 | `mutations/ClinVar/Missense_filter_mutations_mapped.tsv` (+ Frameshift/Nonsense/Indel) | `create_mutation_map_worker.py` |
| DepMap somatic mutations | 8e | `unmapped/annotations/depmap_mutations.tsv` | `create_depmap_worker.py` |
| Exon boundaries | 5d | `genome/exon.tsv` | `create_exon_worker.py` |
| phastCons 100-vertebrate conservation | 7 | `unmapped/conservation/conservation_phastcons.tsv` | `create_conservation_worker.py` |

### Notes on genome-level data

- All genomic coordinates are hg38 (GRCh38).
- `combined_map.map` provides the protein-position ↔ genomic-position bridge used by mutation mapping, exon mapping, and phastCons.
- ClinVar auto-downloads from NCBI FTP on first run (cached via `storeDir`). Supply `--clinvar_vcf` to use a local copy.
- DepMap and phastCons require pre-processed local files (`--depmap_tsv`, `--phastcons_dir`).

---

## Protein-level annotations

These annotations operate on protein sequences and UniProt accessions. They do not depend on Module 3 and are computed in parallel with genome mapping.

### Linear motifs and post-translational modifications

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| ELM linear motifs | 5a | `unmapped/annotations/elm.tsv` | [elm.md](elm.md) |
| DIBS (disordered binding sites) | 5a | `unmapped/annotations/dibs.tsv` | |
| MFIB (molecular function in intrinsically disordered) | 5a | `unmapped/annotations/mfib.tsv` | |
| PhasePro (phase-separation drivers) | 5a | `unmapped/annotations/phasepro.tsv` | |
| PTM sites (PTMdb + PhosphoSite) | 5a | `unmapped/annotations/ptm_merged.tsv` | [ptm.md](ptm.md) |
| Pfam domains | 5a | `unmapped/annotations/pfam_domains.tsv` | |
| PEM Core Motifs | 5h | `unmapped/annotations/pem_core_motifs.tsv` | |
| ScanSite 4.0 kinase motifs | 5k | `unmapped/annotations/scansite.tsv` | [scansite.md](scansite.md) |
| SNP polymorphisms | 5l | `unmapped/annotations/snp_polymorphisms.tsv` | [snp_polymorphisms.md](snp_polymorphisms.md) |

### Disorder and structure

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| IUPred3 disorder + ANCHOR2 binding | 5b | `unmapped/disorder/IUPredscores.tsv`, `Anchorscores.tsv` | [disorder.md](disorder.md) |
| AIUPred disorder + binding | 5b | `unmapped/disorder/AIUPredscores.tsv`, `AIUPredBinding.tsv` | [disorder.md](disorder.md) |
| AlphaFold pLDDT | 5b | `unmapped/disorder/AlphaFoldTable.tsv` | [disorder.md](disorder.md) |
| Combined disorder (CombinedDisorderNew) | 5b | `unmapped/disorder/CombinedDisorderNew.tsv`, `CombinedDisorderNew_Pos.tsv` | [disorder.md](disorder.md) |
| PDB structures | 5c | `unmapped/pdb/pdb_structures.tsv`, `pdb_regions.tsv`, `pdb_disorder.tsv` | |
| Coiled coils (DeepCoil) | 5i | `unmapped/annotations/coiled_coils.tsv`, `DeepCoil.tsv` | |

### Functional annotation

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| GO terms | 5f | `unmapped/annotations/go_terms.tsv` | |
| UniProt natural variants (polymorphism) | 5g | `unmapped/annotations/polymorphism.tsv` | |
| Protein-protein interactions (PPI) | 5j | `unmapped/annotations/interactions.tsv` | [ppi.md](ppi.md) |

### Conservation

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| GOPHER 7-level conservation | 7 | `unmapped/conservation/conservation_multiple_level.tsv` | [conservation.md](conservation.md) |
| phastCons (also genome-level) | 7 | `unmapped/conservation/conservation_phastcons.tsv` | [conservation.md](conservation.md) |

### Pathogenicity and disease

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| AlphaMissense (GENCODE isoforms) | 8d | `unmapped/annotations/alphamissense.tsv` | [alphamissense.md](alphamissense.md) |
| dbNSFP pathogenicity scores | 8f | `unmapped/annotations/pathogenicity_scores.tsv` | [pathogenicity.md](pathogenicity.md) |
| ClinVar disease ontology + Final_Category | 8a | `unmapped/annotations/clinvar_disease.tsv` | [disease_ontology.md](disease_ontology.md) |
| OMIM disease ontology | 8b | `unmapped/annotations/omim_disease.tsv` | [disease_ontology.md](disease_ontology.md) |
| Cancer Gene Census (CGC) | 8c | `unmapped/annotations/census_driver.tsv` | [cancer_drivers.md](cancer_drivers.md) |
| Cosmic Compendium driver scores | 8c | `unmapped/annotations/compendium_driver.tsv` | [cancer_drivers.md](cancer_drivers.md) |

---

## Mapped vs. unmapped outputs

The `unmapped/` directory holds UniProt-accession-keyed raw outputs. The `TRANSCRIPT_MAP` step (Module 5e) translates these to Gencode transcript names (`Protein_ID`), placing results in `mapped/`. When an annotation is transferred to an isoform via 100% sequence substring match, the `homology_transfer` column is set to `True`.

Files in `unmapped/` must not have a `_mapped` suffix.

---

## Planned (not yet implemented)

| Annotation | Reason not yet done |
|------------|---------------------|
| RSA per-residue scores | Requires DSSP + AlphaFold PDB download pipeline |
| ELM Switches | Dataset not yet integrated |
| ScanSite live API mode | Pre-computed file preferred; API mode exists but not default |
