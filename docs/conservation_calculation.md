# Conservation Score Calculation

## Overview

The pipeline computes two independent, complementary conservation tracks per protein:

| Track | What it measures | Granularity | Source |
|-------|-----------------|-------------|--------|
| **GOPHER multi-level** | Evolutionary conservation via Trident score across 7 taxonomic levels | Per amino acid residue, per taxonomic level | Pre-computed from multi-species ortholog alignment |
| **phastCons** | Genome-level conservation from 100-vertebrate whole-genome alignment | Per amino acid residue (mean of 3 codon nucleotides) | UCSC hg38 phastCons bigWig files |

These two tracks are independent: GOPHER operates at the **protein sequence level** (UniProt-keyed, taxonomically stratified), while phastCons operates at the **genomic DNA level** (coordinate-keyed, all vertebrates). Together they give both evolutionary depth and genome-wide comparative context.

Both are produced by a single Nextflow process, `CONSERVATION_MAP` (Module 7), which calls `bin/create_conservation_worker.py`.

---

## 1. GOPHER Multi-Level Conservation

### What it measures

GOPHER (Gene Orthology-based Phylogenetic Homology and Evolutionary Reconstruction) detects orthologs across the QfO (Quest for Orthologs) reference proteome, builds multiple sequence alignments with MAFFT, removes insertion columns relative to the human sequence, and computes the **Trident score** per position. Trident scores are normalized to 0–1, where 1 = perfectly conserved.

The pipeline uses **seven taxonomic levels**, each reflecting a different evolutionary depth:

| Level | Approximate description |
|-------|------------------------|
| `global` | All QfO species (broadest comparison) |
| `Mammalia` | Mammals only |
| `Vertebrata` | All vertebrates |
| `Eukaryota` | All eukaryotes |
| `Eumetazoa` | All animals (excluding plants/fungi) |
| `Opisthokonta` | Animals + fungi |
| `Viridiplantae` | Green plants (useful for ancient conserved residues) |

### Computation pipeline

GOPHER runs are **pre-computed** for the full human proteome; the Nextflow module only performs a lookup. The conservation scores are stored in a flat TSV `conservation_table.tsv` supplied via `params.conservation_table`.

```mermaid
flowchart TD
    UNI[UniProt canonical sequences]
    UNI --> BLAST[BLAST against QfO reference proteome\n366K sequences across 66 species]
    BLAST --> GOPHER[GOPHER ortholog detection\nreciprocal best hit + synteny]
    GOPHER --> MAFFT[MAFFT multiple sequence alignment\nper gene orthogroup]
    MAFFT --> INS[insertion_free.py\nremove alignment columns\nwith insertions vs. human]
    INS --> TAX[taxonomy.py\nassign taxonomic level\nper ortholog sequence]
    TAX --> CALC[conservation_calculation.py\nTrident score per position\nper taxonomic level]
    CALC --> TABLE[conservation_table.tsv\nuniprot_acc | level | conservation_score]

    TABLE --> NF_IN[Nextflow CONSERVATION_MAP input\nparams.conservation_table]
    NF_IN --> STRIP[Strip isoform suffix\nP04049-2 → P04049]
    STRIP --> LOOKUP[Look up canonical accession\nin gopher_data dict]
    LOOKUP --> OUT[conservation_multiple_level.tsv\nProtein_ID | Entry_Isoform | level | conservationscores]
```

### Input table format

`conservation_table.tsv` has three columns:

| Column | Description |
|--------|-------------|
| `uniprot_acc` | UniProt canonical accession (no isoform suffix) |
| `level` | Taxonomic level string (e.g. `Mammalia`) |
| `conservation_score` | Comma-separated per-residue Trident scores (one per amino acid) |

### Nextflow lookup logic

`create_conservation_worker.py` strips isoform suffixes before lookup (e.g. `P04049-2` becomes `P04049`) so that all GENCODE transcripts of the same gene, even if matched to an isoform-suffixed accession, find the conservation data.

The output `conservation_multiple_level.tsv` has one row per `(Protein_ID, taxonomic_level)` pair:

```
Protein_ID    Entry_Isoform    level         conservationscores
RAF1-201      P04049           Mammalia      0.92,0.87,0.91,...
RAF1-201      P04049           Vertebrata    0.88,0.84,0.89,...
```

If a protein has no GOPHER data (e.g. it is absent from the QfO reference proteome), it is silently skipped. Skip via `--skip_gopher` flag or by passing `NO_FILE` as the `conservation_table` input.

---

## 2. phastCons Conservation

### What it measures

phastCons scores are derived from the **UCSC 100-vertebrate whole-genome alignment** against hg38. They represent the probability that each nucleotide position is under purifying selection, estimated by a phylogenetic hidden Markov model (phylo-HMM). Scores are 0–1.

For protein residues, the pipeline takes the **mean phastCons score of the three coding nucleotides** (the codon) as the per-residue score. This requires knowing the genomic coordinates of each codon, which come from `combined_map.map`.

### Data sources

phastCons bigWig files are stored locally at `/dlab/home/norbi/data/phastcons/` (per-chromosome files: `chr1.bw` ... `chrY.bw`). They are not downloaded by the pipeline; they must be present and the path set via `params.phastcons_dir`.

### Computation pipeline

```mermaid
flowchart TD
    BW[phastCons bigWig files\nchr1.bw ... chrY.bw\nparams.phastcons_dir]
    CMAP[combined_map.map\nfrom GENOME_MAP]
    BW --> PARSE_MAP
    CMAP --> PARSE_MAP[parse_combined_map\nextract genomic coords per residue\n3 nucleotide positions per amino acid]
    PARSE_MAP --> CHROM[For each Protein_ID:\nchrom, strand, residues list]
    CHROM --> BW2BG[bigWigToBedGraph\n-chrom=chrN -start=X -end=Y\noutput temp BedGraph]
    BW2BG --> PARSE_BG[Parse BedGraph\nbuild pos_score dict\ngenomic_position → phastCons_score]
    PARSE_BG --> MEAN[For each residue:\nmean score of 3 codon nucleotide positions\ngaps - → score 0.0]
    MEAN --> OUT[conservation_phastcons.tsv\nProtein_ID | Entry_Isoform | conservationscores]
```

### Strand handling

Genomic coordinates in `combined_map.map` are always on the **positive strand** (UCSC 0-based coordinates). The BLAT alignment and genome map worker handle reverse-complement conversion for minus-strand genes internally. The phastCons worker therefore reads coordinates directly from `combined_map.map` without needing to apply strand correction.

### Gap residues

Some residues in `combined_map.map` have `-` in their genomic coordinate fields (gaps introduced by the BLAT alignment or by UTR regions). For these residues, the phastCons score is set to `0.0`.

### Mitochondrial genes

Chromosome `chrM` genes are skipped entirely (phastCons bigWigs do not include `chrM`). These proteins will have an empty `conservationscores` field.

### Output format

`conservation_phastcons.tsv` has one row per `Protein_ID`:

```
Protein_ID    Entry_Isoform    conservationscores
RAF1-201      P04049           0.823,0.791,0.856,...
```

The `conservationscores` column contains a comma-separated list of floats, one per amino acid residue, in sequence order.

---

## 3. Combined Output Structure and DisCanVis2 Usage

### Files produced

| File | Location | Description |
|------|----------|-------------|
| `conservation_multiple_level.tsv` | `final/conservation/` | GOPHER scores, one row per Protein_ID × level |
| `conservation_phastcons.tsv` | `final/conservation/` | phastCons scores, one row per Protein_ID |

These files are already keyed by `Protein_ID` (GENCODE transcript name) and do not require the `TRANSCRIPT_MAP` transfer step; they land directly in `final/conservation/`.

### Django models

DisCanVis2 stores conservation data in two separate models:

**`Conservation_multiple_level`**
```python
# Fields: protein (FK), level (CharField), conservationscores (TextField)
# One instance per (protein, taxonomic_level) pair
# Displayed as 7 stacked per-residue heatmap tracks in the protein viewer
```

**`Conservation_phastCons`**
```python
# Fields: protein (FK), conservationscores (TextField)
# One instance per protein
# Displayed as a single per-residue heatmap track
```

Both are displayed as color-coded per-residue tracks in the DisCanVis2 protein sequence viewer, allowing researchers to identify residues that are conserved across different evolutionary depths and also under genome-level selection pressure.

---

## 4. Running Modes

```bash
# Full conservation (GOPHER + phastCons)
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 \
    --conservation_table /path/to/conservation_table.tsv \
    --phastcons_dir /dlab/home/norbi/data/phastcons/ \
    -resume

# GOPHER only (no phastCons bigWigs available)
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 \
    --conservation_table /path/to/conservation_table.tsv \
    --skip_phastcons -resume

# phastCons only (no GOPHER pre-computation)
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 \
    --phastcons_dir /dlab/home/norbi/data/phastcons/ \
    -resume
# (GOPHER is automatically skipped when conservation_table is not set)

# Skip both (fastest test runs)
nextflow run main.nf --project test_one_protein --data local --machine laptop --target_gene RAF1 \
    --skip_phastcons -resume
```

### Parameter reference

| Parameter | Description |
|-----------|-------------|
| `params.conservation_table` | Path to pre-computed GOPHER conservation TSV; if unset, GOPHER step is skipped |
| `params.phastcons_dir` | Directory with per-chromosome `.bw` files; if unset, phastCons step is skipped |
| `params.bigwigtobedgraph` | Path to `bigWigToBedGraph` binary (default: uses `$PATH`) |
| `--skip_gopher` | Force-skip GOPHER even if `conservation_table` is set |
| `--skip_phastcons` | Force-skip phastCons even if `phastcons_dir` is set |
