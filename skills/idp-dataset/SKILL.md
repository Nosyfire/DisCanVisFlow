---
name: idp-dataset
description: >
  Generate protein annotation datasets with DisCanVisFlow — per-residue disorder, linear
  motifs, PTMs, domains, structure, conservation, phase-separation and variant tracks for
  human genes. Use this skill whenever someone wants data for one or more proteins:
  "give me IDP data for RAF1", "I need disorder annotations for these 40 genes", "annotate
  my gene list", "build a feature matrix for my IDPs", "which regions of CTNNB1 are
  disordered", "get me ELM motifs and PTMs for TP53", "run the pipeline on BRAF", or a
  bare gene list with no verb at all. It also covers the setup question — "how do I get
  this running on this machine" — because the first thing this skill does is locate the
  pipeline and tell you what is missing. Use it even when the user never says DisCanVis,
  Nextflow, or pipeline; a request for protein annotation data is enough. For interpreting
  data that already exists, prefer discanvis-analysis; for a failing run, prefer
  discanvis-troubleshoot.
---

# Generating DisCanVisFlow datasets

The pipeline maps disease variants, functional annotations and structural features onto
every curated protein isoform in the human SwissProt proteome. This skill turns a
natural-language data request into the right command.

The single most valuable thing you can do here is **check for an existing run before
computing anything**. Slicing a finished full-proteome run takes seconds; recomputing the
same genes takes minutes to hours, and a fresh full proteome takes about a day. Users
rarely know a finished run is sitting there, so checking is on you, not them.

## Step 1 — Locate the pipeline

Run this before anything else. It answers where the checkout is, whether the conda
environment is built, and — crucially — which runs have already finished:

```bash
python "$(dirname "$0")/../../bin/find_discanvis.py"   # from inside the plugin
# or simply, when the repo is the working directory:
python bin/find_discanvis.py
```

Add `--json` when you want to branch on the result programmatically.

If it reports no checkout, the user needs one. Offer the three commands it prints
(clone, create the conda env, export `DISCANVIS_HOME`) and stop until that is done —
everything below depends on it.

If it reports the pipeline was found **in the plugin cache**, say so before running
anything long. That copy works, but a run writes tens of gigabytes of reference data into
a versioned plugin directory and the next plugin upgrade throws it away. A normal clone
with `DISCANVIS_HOME` pointing at it is what you want for real work.

Everything below assumes commands run from the checkout root.

## Step 2 — Understand the request

- **Which proteins?** One gene, a comma-separated list, a file of HGNC symbols, or the
  whole proteome. If given UniProt accessions instead of gene names, map them via
  `results/<project>/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv`.
- **Which tracks?** Everything (the default) or a subset. A subset is worth proposing when
  the user's question is narrow — a disorder question does not need dbNSFP or conservation,
  and dropping them is the difference between minutes and hours.
- **Where should it land?** Default to `results/idp_<gene>/`.

The pipeline covers the **human proteome only**. Say so plainly if asked for another
species rather than producing an empty result.

## Step 3 — Extract from an existing run, or run the pipeline

### Path A — a finished run exists (strongly preferred)

`find_discanvis.py` lists finished runs. A `discanvis` run is the full proteome, so any
gene the user asks for is already in it:

```bash
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis --gene RAF1 --out results/idp_RAF1

# several genes
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis --gene RAF1,TP53,BRAF --out results/idp_kinases

# from a file (one HGNC symbol per line, # comments allowed)
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis --gene_list_file my_genes.txt --out results/idp_custom
```

This filters every TSV by `Protein_ID` prefix. No recomputation, no conda/Nextflow risk.

### Path B — no finished run, or the user wants tracks the old run lacks

```bash
# one gene
nextflow run main.nf --project test_one_protein --target_gene RAF1 -resume

# a gene list
nextflow run main.nf --project discanvis --gene_list_file my_genes.txt -resume

# the whole proteome (~24 h on a big server)
nextflow run main.nf --project discanvis -resume
```

Add `-stub` to validate the workflow graph without computing anything — a cheap way to
confirm a command is well-formed before committing hours to it.

## Step 4 — Pick portable flags

**Default to `--data discanvis_data`, which is also the config default, so simply omit the
flag.** It downloads every open reference automatically and caches it. `--data local`
reads `config/data/local.config`, which is machine-specific and deliberately not in git —
on any machine that has not been set up by hand it will fail. Only reach for it when
`find_discanvis.py` reported `local.config: present`.

`--machine` should match the hardware, and getting it wrong is the most common cause of a
run dying on memory. BLAT loads a ~4 GB genome per parallel job:

| Hardware | Flag |
|---|---|
| 8 GB laptop, WSL | `--machine laptop` |
| 32 GB workstation | `--machine low` |
| 64 GB+ workstation | `--machine medium` |
| Dedicated 256 GB server | `--machine hard` |
| SLURM cluster | `--machine slurm` |

When you do not know the machine, `--machine medium` is the safe middle. You can check
with `nproc` and `free -g`.

## Step 5 — Narrow the tracks when it helps

`--modules` takes a comma-separated list and runs only those groups. The backbone (ELM,
Pfam, DIBS, MFIB, PhasePro, PTM) always runs and cannot be selected or excluded.

```bash
--modules mutations,disorder          # variants + disorder only
--modules disorder,disprot,mobidb     # predicted + curated disorder
--modules disprot,mobidb              # curated evidence only, no predictors
```

The authoritative list of module names, with what each one gates, is in
`docs/guide/configuration.md`. Read it rather than reciting names from memory — the set
grows, and a wrong name is accepted silently and simply produces nothing.

For finer control, `--skip_*` flags disable individual predictors inside a module
(`--skip_iupred`, `--skip_alphafold`, and so on); the full table is in the same document.

## Step 6 — Propose, confirm, run

Long runs are expensive and hard to undo, so state the plan before starting one:

1. Two or three lines: which proteins, extract or run, expected duration.
2. The exact command in a code block.
3. Ask whether to run it, and wait for a real answer.
4. Run it, then report where the output landed and roughly how much there is:
   ```bash
   wc -l results/<project>/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv
   ```

Extraction is fast and reversible, so for Path A a brief "extracting RAF1 from the
existing full run" and going ahead is fine. Reserve the confirmation step for actual
pipeline runs.

## What gets produced

Outputs land in `results/<project>/final/`, tab-separated, keyed by `Protein_ID` (the
GENCODE transcript name, e.g. `RAF1-201`). Per-residue scores are comma-separated arrays,
one value per residue in sequence order.

| Directory | Holds |
|---|---|
| `disorder/` | IUPred3, ANCHOR2, AIUPred, MobiDB, DisProt, combined disorder |
| `structure/` | AlphaFold pLDDT, RSA, DSSP, PDB coverage and unobserved regions |
| `annotations/` | ELM, DIBS, MFIB, PhasePro, PTM, Pfam, GO, PPI, coiled coils, low-complexity, aggregation-prone regions, LLPS regions, polymorphism |
| `mutations/` | ClinVar, TCGA, cBioPortal, DepMap, LLPS-associated variants |
| `pathogenicity/` | dbNSFP, AlphaMissense, MaveDB, ProteinGym |
| `phase_separation/` | catGRANULE, PLAAC |
| `disease/`, `drivers/` | ClinVar/MONDO and OMIM disease, cancer driver tables |
| `conservation/`, `sequence/`, `genome/`, `position/` | Conservation, isoform table, coordinate maps, position-based annotations |

`docs/annotations/README.md` indexes every track with a page per track explaining its
columns. Point users there, and read it yourself before explaining a column — it is kept
current, whereas anything restated here will drift.

For thresholds and how to read the scores, see
[references/interpreting-scores.md](references/interpreting-scores.md).
