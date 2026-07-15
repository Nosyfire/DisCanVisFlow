# FINCHES — LLPS Saturation Mutagenesis

## Description

FINCHES computes, for every possible single amino-acid substitution in a
sequence, how the mutation changes the region's **liquid-liquid phase separation
(LLPS)** tendency, expressed as a change in interaction energy (ε). It is an
in-silico saturation-mutagenesis scan of condensate-forming propensity.

## Data source

- **Predictor:** FINCHES, run on each isoform sequence.
- **Origin:** FINCHES (Ginell, Holehouse et al.).
- **Update policy:** Recomputed from sequences each run.
- **Off by default:** enable with `--skip_finches false`. Licensed CC BY-NC 4.0.

## Output file

`final/pathogenicity/finches_saturation.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `Position` | Residue position (1-based) |
| `WT_AA` | Wild-type amino acid |
| `Mut_AA` | Substituted amino acid |
| `WT_Epsilon` | Interaction energy (ε) of the wild-type context |
| `Mut_Epsilon` | Interaction energy (ε) after the substitution |
| `Delta_Epsilon` | Δε = Mut − WT |

## Compute engines (`--finches_engine`)

The saturation scan touches every position × 19 substitutions, so speed matters
at proteome scale. Two engines produce **byte-identical** output:

| Engine | How | Cost / protein | Speed |
|--------|-----|----------------|-------|
| `incremental` (**default**) | rebuilds only the band of the L×L interaction matrix a point mutation actually changes | O(20·L²) | reference |
| `full` | one `calculate_epsilon_value(mut,mut)` per variant (rebuilds the whole matrix each time) | O(19·L³) | ~90× slower |

The FINCHES self-epsilon reduces to a normalised sum of an elementwise function of
the weighted matrix, `ε = (1/L)·Σ_ij h(w_ij)`, and a single substitution at
position *p* only perturbs row/col *p* plus the ±1 charge-window and any
aliphatic-cluster it changes. The incremental engine (`bin/finches_incremental.py`)
recomputes just that band and updates the cached WT sum. It is validated
bit-for-bit against `full` (max abs Δε ≈ 1e-14 across charged/aliphatic/terminal
test cases and real proteins). Benchmark: **RAF1-201 (648 aa) ≈ 32 min (`full`) →
≈ 22 s (`incremental`)** single-core.

Set the engine via `--finches_engine full` to force the reference path.

## Non-standard residues

The substitution alphabet is the standard 20, so a position is only mutated if
its wild-type residue is one of them. Residues outside that alphabet still count
as **sequence context** and contribute to ε.

**U (selenocysteine)** *is* parameterised by Mpipi. Selenoproteins are scored
normally, except that the U positions themselves are never mutated — so their
blocks hold `19 × (number of standard residues)` rows rather than `19 × length`.
The incremental engine's cached matrices cover the standard 20 only, so a
U-containing isoform falls back to the `full` engine and pays its ~90× cost;
with 25 such isoforms proteome-wide (the longest ~670 aa) this adds roughly half
an hour to the tail of a full run, in parallel with everything else.

**X (unknown residue)** has no parameters, so `calculate_epsilon_value` raises
and ε is undefined for the entire sequence — not just that site. Those isoforms
are skipped with a `WT epsilon failed … ('X')` warning and appear nowhere in the
output: **111 of 19,360** main isoforms on the current SwissProt/GENCODE
reference, mostly very short fragments (median ~21 aa).

The driver does not hardcode which residues are parameterised — it asks the
forcefield and skips a protein only if the reference path actually fails.
Isoforms longer than `--max_seq_len` (default 3000 aa) are skipped as well.

## Running / resuming

- Parallelism is **per protein** (`imap_unordered` over a length-sorted stream),
  so all workers stay saturated; a single **tqdm** progress bar reports
  proteins done / skipped / ETA.
- The run is **resumable**: `Protein_ID`s already present in
  `finches_saturation.tsv` are skipped and new results appended. Re-invoke the
  same command to continue an interrupted run (use `--no_resume` to recompute).
- `--order short_first` (default) front-loads short isoforms for fast early
  coverage; `--validate N` cross-checks the first *N* proteins against `full`.

## Notes

- **Positive Δε** = the mutation *increases* LLPS tendency; **negative Δε** =
  *decreases* it.
- Complements the curated [PhasePro](phasepro.md) phase-separation
  regions with a per-mutation quantitative scan.
- Non-commercial licence (CC BY-NC 4.0); disabled by default for that reason.
- Workers: `bin/create_finches_worker.py` (driver, Module 8h) +
  `bin/finches_incremental.py` (exact incremental engine).
