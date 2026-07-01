---
name: idp-dataset
description: >
  Use this skill whenever someone asks for an IDP (intrinsically disordered protein) dataset,
  disorder-focused annotations, a feature matrix for disordered proteins, or wants to annotate
  a protein/gene list with disorder predictions and IDP biology tracks. Triggers on requests like
  "give me IDP data for RAF1", "I need disorder annotations for these proteins", "create an
  annotated IDP feature set", "generate a dataset for my IDPs", "map disorder for EGFR", or
  "I want per-residue IDP annotations for this gene list" — even if the user doesn't say
  DisCanVisFlow or pipeline explicitly. Always use this skill for DisCanVisFlow-generated IDP
  data requests, including single proteins, gene lists, and full-proteome jobs.
---

# IDP Dataset Generation — DisCanVisFlow

This skill translates a natural-language IDP data request into the right pipeline commands,
checks for existing completed runs first (fastest path), and either extracts or runs as needed.

## IDP-relevant outputs produced

All outputs land in `results/<project>/final/` as tab-separated TSVs keyed by `Protein_ID`
(GENCODE transcript name, e.g. `RAF1-201`):

| Directory | Contents |
|-----------|----------|
| `disorder/` | IUPredscores, AnchorScores, AIUPredscores, AIUPredBinding, AlphaFoldTable (pLDDT), MobiDB, CombinedDisorderNew |
| `annotations/` | ELM motifs, DIBS, MFIB, PhasePro, PTM, Pfam, PEM, GO terms, coiled-coils, PPI interactions, ScanSite |
| `sequence/` | Isoform table with full AA sequences, GENCODE coordinates, MANE/APPRIS flags |
| `position/` | RSA scores, position-based annotations |
| `pdb/` | PDB structure coverage, unobserved/disordered regions |
| `conservation/` | phastCons + GOPHER conservation (optional; skip for speed) |

---

## Step 1: Understand the request

Before generating commands, identify:
- **Proteins**: single gene name → comma-separated list → gene list file → full proteome
- **Format of input**: HGNC symbols (e.g. `RAF1`), UniProt accessions (convert to gene name first), or a `.txt` file
- **Output destination**: ask if not obvious (default: `results/idp_<label>/`)
- **Special needs**: specific annotation subsets, data format, downstream use (ML, web upload, exploration)

If the user gives UniProt accessions instead of gene names, map them to HGNC symbols first
(UniProt entry page or grep in `results/discanvis/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv`).

**The pipeline only covers the human proteome (UniProt SwissProt × GENCODE).**
Flag if the user asks for non-human proteins.

---

## Step 2: Check for an existing full-proteome run

Always check this first — extraction takes seconds vs. hours for a fresh run:

```bash
ls results/discanvis/final/disorder/ 2>/dev/null | head -3
```

If `CombinedDisorderNew.tsv` (or any TSV) appears, the full proteome run is complete and ready to extract from.

---

## Step 3: Choose the right approach

### Path A — Existing full-proteome run → extract (preferred)

Use `bin/extract_gene_from_results.py`, which filters every TSV by `Protein_ID` prefix:

```bash
# Single gene
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis \
    --gene RAF1 \
    --out results/idp_RAF1

# Comma-separated genes
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis \
    --gene RAF1,TP53,BRAF,KRAS,EGFR \
    --out results/idp_custom

# Gene list from file (one HGNC symbol per line; # comments OK)
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source results/discanvis \
    --gene_list_file /path/to/my_idps.txt \
    --out results/idp_my_dataset
```

For a full-proteome request, the data is already in `results/discanvis/` — just point there directly.

### Path B — No existing run → run the pipeline

```bash
conda activate discanvis

# Single gene (~4–10 min on 64-CPU server, ~30 min on laptop)
nextflow run main.nf \
    --project test_one_protein \
    --data local \
    --machine hard \
    --target_gene RAF1 \
    -resume

# Gene list from file
nextflow run main.nf \
    --project discanvis \
    --data local \
    --machine hard \
    --gene_list_file /path/to/my_idps.txt \
    -resume

# Full human proteome (~24 h on 64-CPU server)
nextflow run main.nf \
    --project discanvis \
    --data local \
    --machine hard \
    -resume
```

**Speed-up: skip non-IDP tracks** (add to any command above if user wants only disorder biology):
```
    --skip_polymorphism true     # dbSNP SNPs — not needed for IDP biology
    --skip_conservation true     # GOPHER phylogenetic conservation — optional
    --skip_mavedb true           # MaveDB DMS scores — optional
    --skip_proteingym true       # ProteinGym fitness scores — optional
    --skip_clinvar true          # ClinVar disease variants — optional for IDP focus
```

---

## Step 4: Pick the right machine and data flags

| Environment | `--machine` | `--data` |
|-------------|------------|---------|
| Server (64+ CPUs, this machine) | `hard` | `local` |
| Laptop / low RAM (≤16 GB) | `laptop` | `discanvis_data` |
| SLURM cluster | `slurm` | `local` |

Use `--data discanvis_data` when reference files aren't pre-downloaded; it auto-fetches everything.
Use `--data local` when references are already in place (faster, no downloads).

---

## Step 5: Propose → confirm → run

1. **Summarise**: 2–3 lines — which proteins, which approach (extract vs. pipeline run), expected time
2. **Show the exact command(s)** in a code block
3. **Ask**: "Shall I run this?" — wait for explicit confirmation before executing
4. **Run** via Bash tool; show output as it streams
5. **Report** where outputs landed and a quick count (how many genes, how many TSV files)

---

## Output summary to give the user

After the job completes, report:

```
results/<project>/final/
  disorder/         → IUPred3, ANCHOR2, AIUPred, AlphaFold pLDDT, CombinedDisorder, MobiDB
  annotations/      → ELM, DIBS, MFIB, PhasePro, PTM, Pfam, PEM, GO, coiled-coils, PPI
  sequence/         → isoform sequences + genomic coordinates
  position/         → RSA + position-based annotations
  pdb/              → PDB coverage + unobserved regions
```

Mention the protein/isoform count if you can get it:
```bash
wc -l results/<project>/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv
```
