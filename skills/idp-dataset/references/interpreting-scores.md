# Reading DisCanVisFlow scores

`docs/annotations/` documents what every column *is*, track by track, and it is kept
current. This file covers what that reference deliberately leaves out: the thresholds and
judgement calls needed to turn a number into a statement about biology.

When the two disagree, `docs/annotations/` wins on column names and file paths; this file
wins on interpretation.

## Reading the files at all

Every output is a TSV under `results/<project>/final/`, keyed by `Protein_ID` (the GENCODE
transcript name, e.g. `RAF1-201`). One protein has several isoforms and therefore several
rows — decide early whether the question is about a specific isoform or the canonical one,
because summing across isoforms double-counts.

Per-residue tracks store the whole protein in one cell as comma-separated floats, in
sequence order, position 1 first. Two consequences worth internalising:

- Python's `csv` module refuses these by default. Call
  `csv.field_size_limit(2**31 - 1)` before reading, or you get a confusing
  `_csv.Error: field larger than field limit`.
- The array length equals the sequence length. If it does not, you are joining rows from
  different isoforms — a common and silent source of off-by-one position errors.

Coordinates in the annotation tracks are 1-based and inclusive. `combined_map.map` in
`genome/` is the exception: it is 0-based, and Module 4 normalises it internally.

## Disorder

| Track | Disordered when | Notes |
|---|---|---|
| IUPred3 (`IUPredscores`) | > 0.5 | The conventional cutoff. A stretch of >30 consecutive residues above it is a long IDR, which is the unit most IDP papers actually discuss. |
| AIUPred (`AIUPredscores`) | > 0.5 | Deep-learning successor; generally sharper boundaries than IUPred3. |
| AlphaFold pLDDT (`Plldtscores`) | < 50 | Strong disorder signal. 50–70 is genuinely ambiguous — flexible or low-confidence, not proof of either. Above 70 is a confident fold. |
| MobiDB | region-based | Consensus of predictors *and* curated evidence, so it is not an independent vote. |
| DisProt | region-based | Literature-curated with experimental backing. The closest thing to ground truth; absence means "not curated", never "ordered". |

`CombinedDisorderNew.tsv` already integrates MobiDB, pLDDT, IUPred3 and Pfam exclusion into
a binary call, filtered to runs of at least 5 residues. Prefer it over re-deriving a
consensus by hand, and note that it is already `Protein_ID`-keyed so it needs no mapping.

The single most common analytical mistake is treating agreement between IUPred3, AIUPred
and ANCHOR2 as independent confirmation. They share the same energy-estimation lineage;
they agree because they are related, not because the evidence is strong. Independent
support comes from DisProt, from unobserved regions in `structure/pdb_missing.tsv`, or
from low pLDDT.

## Binding within disorder

ANCHOR2 (`AnchorScore`) and AIUPred-Binding (`AIUPredBinding`) both predict regions that
fold on binding a structured partner. Peaks above 0.5 sitting inside an otherwise
disordered stretch are the interesting signal — a high score in an already-ordered region
usually means nothing. The two use different models, so here agreement *is* worth something.

## Structure

`rsa_scores.tsv` holds relative solvent accessibility. Above ~0.25 is surface-exposed,
below is buried; buried residues in a predicted-disordered region are a contradiction
worth investigating rather than reporting.

`dssp.tsv` gives 8-state and 3-state secondary structure from the AlphaFold model — so it
describes a *prediction*, not an experiment, and inherits the model's confidence.

`pdb_missing.tsv` lists residues in deposited structures that were not resolved. This is
experimental evidence of flexibility and is independent of every predictor above, which
makes it disproportionately valuable when a claim needs support.

## Motifs, PTMs, domains

ELM instances are curated and specific; PEM motifs are *predicted* and should be described
that way. Both matter mostly when they fall inside a disordered region — a linear motif
buried in a folded domain is usually not functional, and that cross-check (motif position
against `CombinedDisorderNew.tsv`) is the analysis users actually want.

Pfam domains are the natural negative control for disorder: the combined disorder track
already excludes them.

`mapping_type` appears on annotation tracks that were transferred between isoforms.
`direct` means the same UniProt accession. `homology_similarity` means the region was
carried to a different isoform because it aligned at or above the identity threshold
(default 0.90) — real, but inferred, and worth flagging when a result rests on it.

## Variants

`Mutation` is always the protein-level change (`V600E`). `Genomic_HGVS` carries the
source's genomic expression when one was derived. For frameshifts, `Mutation` is only a
position marker (`L10fs`) — the downstream consequence is not computed, so do not report
it as if it were.

`isoform_mapped` marks variants propagated to sibling isoforms by codon context. Numbering
legitimately differs between isoforms, so `V600E` on one transcript can be `V605E` on
another; both are correct for their own isoform.

ClinVar submission dates are keyed on the variant, not on `Protein_ID`. Join them via
`Genomic_HGVS`, falling back to `Mutation` when that column is empty.

## Phase separation and aggregation

catGRANULE and PLAAC are sequence-composition predictors — they answer "does this look
like known phase-separating protein", which is weaker than experimental evidence. PhasePro,
PhaSepDB and LLPSDB are curated; prefer them when they have an entry, and say which kind
of evidence a claim rests on.

Aggregation-prone regions come from the AGGRESCAN a3v scale with a threshold calibrated to
the measured proteome distribution, so scores are comparable across proteins in this
dataset but not directly against published AGGRESCAN numbers.

## Where to go next

- Column-level reference for any track: `docs/annotations/README.md` and the per-track page it links.
- How UniProt accessions become `Protein_ID`s: `docs/pipeline/isoform_mapping.md`.
- Conservation scoring details: `docs/pipeline/conservation_method.md`.
- Module structure and output layout: `docs/pipeline/architecture.md`.
