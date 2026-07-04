# Annotation Reference — DisCanVisFlow

Every annotation track the pipeline produces, with a dedicated page describing
its biological meaning, data source, output file(s), and columns. Pages are
organised into category folders that mirror the `results/<project>/final/`
output directories.

All mapped outputs are keyed by `Protein_ID` (the GENCODE transcript name). How
UniProt-keyed annotations become isoform-keyed is described in
[Staging vs. final outputs](#staging-vs-final-outputs) below and in
[Isoform mapping](../pipeline/isoform_mapping.md).

---

## Mutations — [`mutations/`](mutations)

Somatic and germline variants mapped from genomic coordinates to isoform
residues (genome-anchored; require `params.hg38_2bit`).

| Track | Output | Page |
|-------|--------|------|
| ClinVar germline variants | `mutations/ClinVar/*_filter_mutations_mapped.tsv` | [clinvar.md](mutations/clinvar.md) |
| cBioPortal somatic mutations | `mutations/CBioportal/*` | [cbioportal.md](mutations/cbioportal.md) |
| TCGA somatic mutations | `mutations/TCGA/*` | [tcga.md](mutations/tcga.md) |
| DepMap somatic mutations | `mutations/DepMap/depmap_mutations.tsv` | [depmap.md](mutations/depmap.md) |

## Pathogenicity — [`pathogenicity/`](pathogenicity)

Predicted and experimentally measured variant effects.

| Track | Output | Page |
|-------|--------|------|
| AlphaMissense | `pathogenicity/alphamissense.tsv` | [alphamissense.md](pathogenicity/alphamissense.md) |
| dbNSFP (14 predictors) | `pathogenicity/dbnsfp_scores.tsv` | [dbnsfp.md](pathogenicity/dbnsfp.md) |
| MaveDB (MAVE assays) | `pathogenicity/mavedb.tsv` | [mavedb.md](pathogenicity/mavedb.md) |
| ProteinGym (DMS benchmarks) | `pathogenicity/proteingym.tsv` | [proteingym.md](pathogenicity/proteingym.md) |

## Disease — [`disease/`](disease)

| Track | Output | Page |
|-------|--------|------|
| ClinVar disease ontology (MONDO) | `disease/clinvar_disease.tsv` | [disease_ontology.md](disease/disease_ontology.md) |
| OMIM disease + mutations | `disease/omim_disease.tsv` | [disease_ontology.md](disease/disease_ontology.md) |

## Cancer drivers — [`drivers/`](drivers)

| Track | Output | Page |
|-------|--------|------|
| Cancer Gene Census (CGC) | `drivers/census_driver.tsv` | [cancer_drivers.md](drivers/cancer_drivers.md) |
| COSMIC Compendium driver scores | `drivers/compendium_driver.tsv` | [cancer_drivers.md](drivers/cancer_drivers.md) |

## Disorder — [`disorder/`](disorder)

Per-residue intrinsic disorder from multiple predictors.

| Track | Output | Page |
|-------|--------|------|
| IUPred3 + ANCHOR2 | `disorder/IUPredscores.tsv`, `Anchorscores.tsv` | [disorder.md](disorder/disorder.md) |
| AIUPred disorder + binding | `disorder/AIUPredscores.tsv`, `AIUPredBinding.tsv` | [disorder.md](disorder/disorder.md) |
| AlphaFold pLDDT | `disorder/AlphaFoldTable.tsv` | [disorder.md](disorder/disorder.md) |
| Combined disorder | `disorder/CombinedDisorderNew.tsv` | [disorder.md](disorder/disorder.md) |
| MobiDB consensus disorder | `disorder/mobidb_disorder.tsv` | [mobidb.md](disorder/mobidb.md) |

## Disorder-associated function — [`disorder_function/`](disorder_function)

Functional features that reside in — or act through — intrinsically disordered
regions. Short linear motifs (SLiMs) are grouped in a [`motifs/`](disorder_function/motifs)
sub-folder.

| Track | Output | Page |
|-------|--------|------|
| ELM linear motifs (SLiMs) | `annotations/elm.tsv` (+ `elm_classes`, `elmswitches`) | [motifs/elm.md](disorder_function/motifs/elm.md) |
| PEM predicted ELM motifs | `annotations/pem_core_motifs.tsv` | [motifs/pem.md](disorder_function/motifs/pem.md) |
| ScanSite phospho motifs | `annotations/scansite.tsv` | [motifs/scansite.md](disorder_function/motifs/scansite.md) |
| DIBS (disordered binding sites) | `annotations/dibs.tsv` | [dibs.md](disorder_function/dibs.md) |
| MFIB (mutual folding by binding) | `annotations/mfib.tsv` | [mfib.md](disorder_function/mfib.md) |
| PhasePro (phase-separation drivers) | `annotations/phasepro.tsv` | [phasepro.md](disorder_function/phasepro.md) |
| FINCHES (LLPS saturation mutagenesis) | `pathogenicity/finches_saturation.tsv` | [finches.md](disorder_function/finches.md) |

## Ordered-region function — [`order_function/`](order_function)

Functional features of structured / folded regions.

| Track | Output | Page |
|-------|--------|------|
| Pfam domains | `annotations/pfam_domains.tsv` | [pfam.md](order_function/pfam.md) |
| UniProt regions of interest & binding sites | `annotations/uniprot_roi.tsv`, `uniprot_binding.tsv` | [uniprot_features.md](order_function/uniprot_features.md) |

## Post-translational modifications — [`ptm/`](ptm)

| Track | Output | Page |
|-------|--------|------|
| PTM sites (PTMdb + PhosphoSite) | `annotations/ptm_merged.tsv` | [ptm.md](ptm/ptm.md) |

## Structure — [`structure/`](structure)

| Track | Output | Page |
|-------|--------|------|
| PDB coverage + unobserved regions | `pdb/pdb_structures.tsv`, `pdb_missing.tsv` | [pdb.md](structure/pdb.md) |
| Coiled coils (DeepCoil) | `annotations/coiled_coils.tsv` | [coiled_coils.md](structure/coiled_coils.md) |
| RSA & position-based annotations | `disorder/rsa_scores.tsv`, `position/position_based_annotations.tsv` | [rsa.md](structure/rsa.md) |

## Interactions — [`interactions/`](interactions)

| Track | Output | Page |
|-------|--------|------|
| Protein-protein interactions (IntAct/BioGRID/HIPPIE) | `annotations/interactions.tsv` | [ppi.md](interactions/ppi.md) |

## Conservation — [`conservation/`](conservation)

| Track | Output | Page |
|-------|--------|------|
| GOPHER 7-level conservation | `conservation/conservation_multiple_level.tsv` | [conservation.md](conservation/conservation.md) |
| phastCons 100-vertebrate | `conservation/conservation_phastcons.tsv` | [conservation.md](conservation/conservation.md) |

See also the [conservation method](../pipeline/conservation_method.md) deep-dive.

## Polymorphism — [`polymorphism/`](polymorphism)

| Track | Output | Page |
|-------|--------|------|
| dbSNP 155 common SNPs + allele frequencies | `annotations/polymorphism.tsv` | [polymorphism.md](polymorphism/polymorphism.md) |

## Function — [`function/`](function)

| Track | Output | Page |
|-------|--------|------|
| GO terms (GOA) | `annotations/go_terms.tsv` | [go.md](function/go.md) |
| Exon boundaries | `genome/exon.tsv` | [exon.md](function/exon.md) |

---

## Genome-anchored vs. protein-level

- **Genome-anchored** tracks (mutations, polymorphism, exon, phastCons) require
  `params.hg38_2bit` and use `combined_map.map` to bridge protein position ↔ hg38
  coordinate. They are skipped when genome mapping did not run.
- **Protein-level** tracks (motifs, disorder, domains, GO, PPI, GOPHER) operate on
  sequences/accessions and run in parallel with genome mapping.

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
