# Annotation Overview — DisCanVis Pipeline

This index covers every annotation type produced by the pipeline. The tables
below group annotations by whether they require genomic coordinates (from
`combined_map.map`) or work directly on protein sequences.

Per-track detail pages are organised into category folders that mirror the
`final/` output directories:

| Folder | Tracks |
|--------|--------|
| [`mutations/`](mutations/) | ClinVar, DepMap (cBioPortal / TCGA via config) |
| [`pathogenicity/`](pathogenicity/) | AlphaMissense, dbNSFP |
| [`disease/`](disease/) | ClinVar disease ontology, OMIM |
| [`drivers/`](drivers/) | Cancer Gene Census, Compendium |
| [`disorder/`](disorder/) | IUPred3, ANCHOR2, AIUPred, AlphaFold pLDDT, Combined disorder |
| [`motifs/`](motifs/) | ELM SLiMs, PTMs, ScanSite phospho motifs |
| [`interactions/`](interactions/) | IntAct / BioGRID / HIPPIE PPIs |
| [`conservation/`](conservation/) | GOPHER multi-level, phastCons |
| [`polymorphism/`](polymorphism/) | dbSNP 155 common SNPs |

---

## Genome-level annotations

These annotations require genomic coordinates derived from Module 3 (Genome Mapping). The pipeline must have been run with `params.hg38_2bit` set. Results are placed under `results/<project>/final/mutations/` or `results/<project>/final/conservation/`.

| Annotation | Module | Output file(s) | Worker |
|------------|--------|----------------|--------|
| ClinVar disease variants | 4 | `mutations/ClinVar/Missense_filter_mutations_mapped.tsv` (+ Frameshift/Nonsense/Indel) | `create_mutation_map_worker.py` |
| DepMap somatic mutations | 8e | `mutations/DepMap/depmap_mutations.tsv` | `create_depmap_worker.py` |
| Exon boundaries | 5d | `genome/exon.tsv` | `create_exon_worker.py` |
| phastCons 100-vertebrate conservation | 7 | `final/conservation/conservation_phastcons.tsv` | `create_conservation_worker.py` |

### Notes on genome-level data

- All genomic coordinates are hg38 (GRCh38).
- `combined_map.map` provides the protein-position ↔ genomic-position bridge used by mutation mapping, exon mapping, and phastCons.
- ClinVar auto-downloads from NCBI FTP on first run (cached via `storeDir`). Supply `--clinvar_vcf` to use a local copy.
- DepMap requires a local file when enabled. Download `OmicsSomaticMutations.csv` manually from <https://depmap.org/portal/download/all/> and save it as `references/depmap/OmicsSomaticMutations.csv`; the pipeline preflight confirms it and normalizes it to `references/depmap/depmap_mutations_raw.tsv` if needed. phastCons requires a local `--phastcons_dir` when conservation/phastCons are enabled.

---

## Protein-level annotations

These annotations operate on protein sequences and UniProt accessions. They do not depend on Module 3 and are computed in parallel with genome mapping.

### Linear motifs and post-translational modifications

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| ELM linear motifs | 5a | `final/annotations/elm.tsv` | [elm.md](motifs/elm.md) |
| DIBS (disordered binding sites) | 5a | `final/annotations/dibs.tsv` | |
| MFIB (molecular function in intrinsically disordered) | 5a | `final/annotations/mfib.tsv` | |
| PhasePro (phase-separation drivers) | 5a | `final/annotations/phasepro.tsv` | |
| PTM sites (PTMdb + PhosphoSite) | 5a | `final/annotations/ptm_merged.tsv` | [ptm.md](motifs/ptm.md) |
| Pfam domains | 5a | `final/annotations/pfam_domains.tsv` | |
| PEM Core Motifs | 5h | `final/annotations/pem_core_motifs.tsv` | |
| ScanSite 4.0 kinase motifs | 5k | `final/annotations/scansite.tsv` | [scansite.md](motifs/scansite.md) |
| SNP polymorphisms | 5l | `final/annotations/snp_polymorphisms.tsv` | [snp_polymorphisms.md](polymorphism/polymorphism.md) |

### Disorder and structure

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| IUPred3 disorder + ANCHOR2 binding | 5b | `final/disorder/IUPredscores.tsv`, `Anchorscores.tsv` | [disorder.md](disorder/disorder.md) |
| AIUPred disorder + binding | 5b | `final/disorder/AIUPredscores.tsv`, `AIUPredBinding.tsv` | [disorder.md](disorder/disorder.md) |
| AlphaFold pLDDT | 5b | `final/disorder/AlphaFoldTable.tsv` | [disorder.md](disorder/disorder.md) |
| Combined disorder (CombinedDisorderNew) | 5b | `final/disorder/CombinedDisorderNew.tsv`, `CombinedDisorderNew_Pos.tsv` | [disorder.md](disorder/disorder.md) |
| PDB structures | 5c | `final/pdb/pdb_structures.tsv`, `pdb_regions.tsv`, `pdb_disorder.tsv` | |
| Coiled coils (DeepCoil) | 5i | `final/annotations/coiled_coils.tsv`, `DeepCoil.tsv` | |

### Functional annotation

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| GO terms | 5f | `final/annotations/go_terms.tsv` | |
| UniProt natural variants (polymorphism) | 5g | `final/annotations/polymorphism.tsv` | |
| Protein-protein interactions (PPI) | 5j | `final/annotations/interactions.tsv` | [ppi.md](interactions/ppi.md) |

### Conservation

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| GOPHER 7-level conservation | 7 | `final/conservation/conservation_multiple_level.tsv` | [conservation.md](conservation/conservation.md) |
| phastCons (also genome-level) | 7 | `final/conservation/conservation_phastcons.tsv` | [conservation.md](conservation/conservation.md) |

### Pathogenicity and disease

| Annotation | Module | Output file | Detail |
|------------|--------|-------------|--------|
| AlphaMissense (GENCODE isoforms) | 8d | `final/annotations/alphamissense.tsv` | [alphamissense.md](pathogenicity/alphamissense.md) |
| dbNSFP pathogenicity scores | 8f | `final/annotations/pathogenicity_scores.tsv` | [pathogenicity.md](pathogenicity/dbnsfp.md) |
| ClinVar disease ontology + Final_Category | 8a | `final/annotations/clinvar_disease.tsv` | [disease_ontology.md](disease/disease_ontology.md) |
| OMIM disease ontology | 8b | `final/annotations/omim_disease.tsv` | [disease_ontology.md](disease/disease_ontology.md) |
| Cancer Gene Census (CGC) | 8c | `final/annotations/census_driver.tsv` | [cancer_drivers.md](drivers/cancer_drivers.md) |
| Cosmic Compendium driver scores | 8c | `final/annotations/compendium_driver.tsv` | [cancer_drivers.md](drivers/cancer_drivers.md) |

---

## Staging vs. final outputs

Annotations are first collected per UniProt isoform (`Entry_Isoform`) into
`results/<project>/intermediate/` (UniProt-accession-keyed staging TSVs). The
`TRANSCRIPT_MAP` step (Module 5e) then translates these to GENCODE transcript
names (`Protein_ID`) and writes the DB-ready copies to
`results/<project>/final/`. When a region is transferred to a different isoform
by sequence similarity (≥ `--homology_min_identity`, default 0.90), its row is
flagged `mapping_type=homology_similarity` (same UniProt accession →
`mapping_type=direct`).

See [Architecture § Key design decisions](../pipeline/architecture.md#key-design-decisions)
and [Isoform mapping](../pipeline/isoform_mapping.md) for the transfer logic.

---

## Planned (not yet implemented)

| Annotation | Reason not yet done |
|------------|---------------------|
| RSA per-residue scores | Requires DSSP + AlphaFold PDB download pipeline |
| ELM Switches | Dataset not yet integrated |
| ScanSite live API mode | Pre-computed file preferred; API mode exists but not default |
