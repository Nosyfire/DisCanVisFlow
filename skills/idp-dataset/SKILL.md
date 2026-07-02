---
name: idp-dataset
description: >
  Use this skill whenever someone asks for an IDP (intrinsically disordered protein) dataset,
  disorder-focused annotations, a feature matrix for disordered proteins, or wants to annotate
  a protein/gene list with disorder predictions and IDP biology tracks. Triggers on requests like
  "give me IDP data for RAF1", "I need disorder annotations for these proteins", "create an
  annotated IDP feature set", "generate a dataset for my IDPs", "map disorder for EGFR",
  "which regions of CTNNB1 are disordered?", "what ELM motifs does TP53 have?", or
  "I want per-residue IDP annotations for this gene list" — even if the user doesn't say
  DisCanVisFlow or pipeline explicitly. Always use this skill for DisCanVisFlow-generated IDP
  data requests, including single proteins, gene lists, and full-proteome jobs.
---

# IDP Dataset Generation — DisCanVisFlow

This skill translates a natural-language IDP data request into the right pipeline commands,
checks for existing completed runs first (fastest path), and either extracts or runs as needed.

## Documentation to consult

Before analysing data or explaining what an annotation means, read the relevant reference:

| Need | File to read |
|------|-------------|
| **What each annotation column means** | `skills/idp-dataset/references/annotations.md` |
| Pipeline modules, processes, output structure | `PIPELINE_DESIGN.md` |
| How UniProt ↔ GENCODE isoform mapping works | `docs/isoform_mapping.md` |
| Conservation score details (GOPHER / phastCons) | `docs/conservation_calculation.md` |
| Performance and bottleneck notes | `docs/performance_benchmark.md` |

Read `references/annotations.md` whenever the user asks what an annotation means, wants to
interpret scores, or asks about a specific track (disorder, ELM, PTM, RSA, etc.).

---

## IDP-relevant outputs produced

All outputs land in `results/<project>/final/` as tab-separated TSVs keyed by `Protein_ID`
(GENCODE transcript name, e.g. `RAF1-201`). Per-residue scores are comma-separated arrays
(one float per amino acid position, in sequence order).

| Directory | Key files | What they capture |
|-----------|-----------|------------------|
| `disorder/` | IUPredscores, AnchorScores, AIUPredscores, AIUPredBinding, AlphaFoldTable, CombinedDisorderNew | Per-residue disorder probability and binding-region predictions |
| `annotations/` | elm, dibs, mfib, phasepro, ptm_merged, pfam_domains, pem_core_motifs, mobidb_disorder, coiled_coils, uniprot_roi, uniprot_binding, go_terms, scansite, interactions | Functional sites, motifs, PTMs, domains, PPI |
| `sequence/` | loc_chrom_with_names_isoforms_with_seq.tsv | Isoform table: sequences, genomic coords, UniProt mapping |
| `position/` | rsa_scores, position_based_annotations | Relative solvent accessibility, secondary structure env |
| `pdb/` | pdb_structures, pdb_missing | PDB coverage + unobserved (candidate disordered) regions |
| `conservation/` | conservation_multiple_level, conservation_phastcons | Evolutionary conservation at 7 taxonomic levels + vertebrate genome |
| `mutations/ClinVar/` | Missense/Frameshift/Nonsense/Indel_filter_mutations_mapped | Clinical variants mapped to protein residues |

For full column-level descriptions, read `skills/idp-dataset/references/annotations.md`.

---

## Workflow

### Step 1: Understand the request

Determine:
- **What proteins?** Single gene name → comma-separated list → gene list file → full proteome
- **Input format?** HGNC gene symbols (e.g. `RAF1`) preferred; if UniProt accessions given, map them via `results/discanvis/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv`
- **Which data?** All IDP tracks (default) or a specific annotation subset?
- **Output destination?** Ask if not obvious (default: `results/idp_<gene>/`)
- **Analysis needed?** If the user wants insights (e.g. "which domain is most mutated"), plan to cross-reference mutation + domain files after extraction

The pipeline covers the **human proteome only** (UniProt SwissProt × GENCODE). Flag non-human requests.

### Step 2: Check for an existing full-proteome run

Always check first — extraction takes seconds vs. hours for a fresh run:

```bash
ls results/discanvis/final/disorder/ 2>/dev/null | head -3
```

If `CombinedDisorderNew.tsv` (or any TSV) appears, the full run is complete.

### Step 3: Choose the right approach

#### Path A — Existing full-proteome run → extract (preferred)

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

For a full-proteome request, `results/discanvis/` is already the dataset.

#### Path B — No existing run → run the pipeline

```bash
conda activate discanvis

# Single gene (~4–10 min on 64-CPU server)
nextflow run main.nf \
    --project test_one_protein \
    --data local \
    --machine hard \
    --target_gene RAF1 \
    -resume

# Gene list
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

**Include only specific annotation groups** with `--modules` (preferred over stacking `--skip_X` flags):

```bash
# Disorder + mutations only (skip PDB, conservation, GO, PPI, etc.)
nextflow run main.nf --project test_one_protein --data local --machine hard \
    --target_gene RAF1 \
    --modules mutations,disorder \
    --fetch_cbioportal true --cbioportal_study msk_impact_2017 \
    --skip_iupred true \
    -resume
```

Available module names: `mutations`, `disorder`, `mobidb`, `pdb`, `go`, `polymorphism`,
`pem`, `coiledcoils`, `ppi`, `conservation`, `scansite`, `clinvar_disease`, `omim`,
`cancer_drivers`, `alphamissense`, `depmap`, `mavedb`, `proteingym`, `dbnsfp`, `finches`

ELM + Pfam + DIBS/MFIB/PhasePro/PTM are backbone annotations — always produced regardless of `--modules`.

To skip individual predictors *within* a module (e.g. skip IUPred3 but keep AIUPred within `disorder`):
```
    --skip_iupred true
    --skip_alphafold true
    --skip_polymorphism true     # within polymorphism module
    --skip_conservation true     # within conservation module
```

### Step 4: Machine / data flags

| Environment | `--machine` | `--data` |
|-------------|------------|---------|
| Server (64+ CPUs) | `hard` | `local` (refs already present) |
| Laptop / low RAM | `laptop` | `discanvis_data` (auto-downloads) |
| SLURM cluster | `slurm` | `local` |

### Step 5: Propose → confirm → run

1. Summarise in 2–3 lines: which proteins, which approach, expected time
2. Show the exact command(s) in a code block
3. Ask "Shall I run this?" — wait for explicit confirmation
4. Run via Bash; stream output
5. Report output location and a quick count:
   ```bash
   wc -l results/<project>/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv
   ```

---

## Analysing the data after extraction

When the user wants insights (e.g. "which domain is most mutated"), follow this pattern:

1. **Read** `skills/idp-dataset/references/annotations.md` to understand the relevant columns
2. **Load** the appropriate TSVs with Python (`csv.field_size_limit(2**31-1)` first — some score arrays exceed the default limit)
3. **Cross-reference**: e.g. mutation positions vs. Pfam domain boundaries, or ANCHOR2 peaks vs. ELM motif positions
4. **Summarise** findings with a table + biological interpretation

Key files for common analyses:

| Question | Files to load |
|----------|--------------|
| Which region is most mutated? | `mutations/ClinVar/Missense_filter_mutations_mapped.tsv` + `annotations/pfam_domains.tsv` + `annotations/uniprot_roi.tsv` |
| Where are the disordered regions? | `disorder/CombinedDisorderNew.tsv` + `pdb/pdb_missing.tsv` |
| What SLiMs / binding sites are in the IDRs? | `annotations/elm.tsv`, `dibs.tsv`, `mfib.tsv`, `pem_core_motifs.tsv` |
| Which PTMs are in disordered regions? | `annotations/ptm_merged.tsv` vs `disorder/CombinedDisorderNew.tsv` |
| Is this protein a phase separator? | `annotations/phasepro.tsv`, `annotations/coiled_coils.tsv` |
| Binding sites vs. structural confidence | `disorder/AlphaFoldTable.tsv` vs `disorder/Anchorscores.tsv` |
