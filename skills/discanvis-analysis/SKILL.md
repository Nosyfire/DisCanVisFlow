---
name: discanvis-analysis
description: >
  Answer biological questions from DisCanVisFlow output that already exists — cross-referencing
  disorder, motifs, PTMs, domains, structure and variants to reach a conclusion rather than
  producing files. Use this skill when someone asks what the data means: "which region of TP53
  is most mutated", "are the ClinVar variants in disordered regions", "what motifs sit in the
  IDRs of CTNNB1", "is this protein a phase separator", "do the PTM sites overlap the binding
  regions", "compare disorder between these two isoforms", "which of my 200 genes have long
  IDRs", or any request to summarise, rank, plot or interpret annotation TSVs. It also covers
  loading these files correctly, which has real pitfalls — per-residue arrays exceed Python's
  default CSV field limit and isoform rows double-count if joined carelessly. Use it whenever
  the data is already on disk; if the data must be generated first, start with idp-dataset.
---

# Analysing DisCanVisFlow output

The outputs are per-residue tracks and region tables for every isoform of every human
protein. Most real questions are a **join across tracks** — variants against domains,
motifs against disorder, PTMs against binding regions — so the analysis is usually more
about lining up coordinates correctly than about any single file.

## Step 1 — Find the data

```bash
python bin/find_discanvis.py
```

It lists finished runs. Analysis needs one; if there are none, this is an
[idp-dataset](../idp-dataset/SKILL.md) job first.

`results/discanvis/` is the full proteome when present, so any gene is already in it. A
narrower directory like `results/idp_RAF1/` holds a previously extracted subset.

## Step 2 — Load the files without corrupting them

Three things routinely go wrong, and all three fail quietly rather than loudly:

```python
import csv
csv.field_size_limit(2**31 - 1)   # per-residue arrays blow the default limit
```

**Rows are per isoform, not per gene.** `Protein_ID` is a transcript name (`RAF1-201`).
A gene has several, and counting variants "in RAF1" by summing rows counts the same
genomic variant once per isoform. Decide whether the question is about one isoform or the
canonical one, and filter before aggregating.

**Per-residue arrays are comma-separated floats, position 1 first.** Their length equals
the sequence length; if a join makes those disagree, the join is wrong and every position
downstream is shifted.

Region coordinates are 1-based inclusive. `genome/combined_map.map` is the one 0-based file.

For a real analysis, prefer pandas and merge on `Protein_ID` explicitly rather than
positionally.

## Step 3 — Pick the files the question actually needs

| Question | Files |
|---|---|
| Where are the disordered regions? | `disorder/CombinedDisorderNew.tsv`, cross-checked against `structure/pdb_missing.tsv` |
| Which region is most mutated? | `mutations/ClinVar/Missense_filter_mutations_mapped.tsv` + `annotations/pfam_domains.tsv` + `annotations/uniprot_roi.tsv` |
| Are variants enriched in IDRs? | the same mutation file + `disorder/CombinedDisorderNew_Pos.tsv` |
| What motifs are in the IDRs? | `annotations/elm.tsv`, `dibs.tsv`, `mfib.tsv`, `pem_core_motifs.tsv` + disorder calls |
| Which PTMs sit in disordered regions? | `annotations/ptm_merged.tsv` + disorder calls |
| Is this a phase separator? | `annotations/phasepro.tsv`, `llps_regions.tsv`, `phase_separation/catgranule.tsv`, `plaac.tsv` |
| How confident is the structure? | `structure/AlphaFoldTable.tsv`, `dssp.tsv`, `rsa_scores.tsv`, `pdb_structures.tsv` |
| How pathogenic are these variants? | `pathogenicity/dbnsfp_scores.tsv`, `alphamissense.tsv` |
| Is the residue conserved? | `conservation/conservation_multiple_level.tsv`, `conservation_phastcons.tsv` |

Read [interpreting-scores.md](../idp-dataset/references/interpreting-scores.md) before
turning any number into a claim — it carries the thresholds (IUPred3 > 0.5, pLDDT < 50,
RSA > 0.25) and the traps. `docs/annotations/` has the per-column reference.

## Step 4 — Reason about evidence, not just numbers

This is where an analysis earns its keep, and where it most often goes wrong.

**Predictors that agree are often not independent.** IUPred3, AIUPred and ANCHOR2 share a
lineage; their agreement reflects that relationship, not accumulating evidence. Genuinely
independent support for disorder comes from DisProt curation, unresolved residues in
`pdb_missing.tsv`, or low pLDDT. Say which kind of evidence a conclusion rests on.

**Absence is not negative evidence.** A protein missing from DisProt is uncurated, not
ordered. A gene with no ClinVar variants may be understudied rather than tolerant. Empty
output frequently means a track was skipped for that run — check the track exists before
concluding it is empty.

**Say when a result rests on inference.** `mapping_type=homology_similarity` means an
annotation was transferred to a different isoform by sequence similarity, and
`isoform_mapped` means a variant was propagated by codon context. Both are legitimate and
both are weaker than a direct observation.

**Normalise before comparing proteins.** Long proteins accumulate more of everything.
Counts of variants, motifs or PTMs are only comparable per residue or per region.

## Step 5 — Report

State the answer first, then the evidence, then the caveat. A short table of the actual
numbers beats prose, and naming the files used lets the user check the work.

Where a chart genuinely helps — a positional track along the sequence, a distribution
across a gene set — matplotlib is in the conda environment. Keep it to cases where the
shape matters; for a handful of numbers a table is clearer.

If the question cannot be answered from the tracks present, say which track would be
needed and how to generate it (`--modules <name>`), rather than answering from a weaker
proxy.
