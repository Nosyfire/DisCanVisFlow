# Aggregation-prone regions (APR)

## Description

Contiguous stretches of a protein whose local amino-acid composition is intrinsically
aggregation-prone. Derived from sequence alone — no structure, no external reference,
no external binary.

## Data source

- **Computed:** a sliding window over each isoform sequence, scored with the
  **AGGRESCAN a3v** intrinsic aggregation-propensity scale. Positive a3v values are
  aggregation-prone.
- **Origin:** Conchillo-Solé O *et al.* (2007) *BMC Bioinformatics* 8:65.
- **Update policy:** Recomputed each run from the isoform sequences.

This is the *scale*, not the published AGGRESCAN HSA/NHSA hot-spot algorithm. That is a
deliberate simplification: it makes the track dependency-free and proteome-scale, at the
cost of not reproducing the AGGRESCAN server's exact hot-spot calls. Validation against
TANGO/Waltz on a narrower set is a separate, later step.

## Algorithm

1. Map every residue to its a3v value. Non-standard residues (X, U, B, Z, …) are NaN.
2. Take the centered sliding-window mean, window **5**, truncated at the termini. A window
   containing any non-standard residue is NaN and can never be part of a region.
3. An APR is a **maximal run of at least 5 consecutive positions** whose window mean is
   at or above the threshold.

### Threshold

In precedence order:

1. `--threshold VALUE` — used verbatim.
2. Input holds at least `--min_proteins_for_threshold` (default 1000) proteins — the
   threshold is the `--quantile` (default 0.90) of every window mean in the input.
3. Otherwise — the module constant `DEFAULT_A3V_THRESHOLD`, with a loud warning. This is
   the single-gene path: a run over one protein must not calibrate a quantile against
   itself.

`DEFAULT_A3V_THRESHOLD` is the measured q90 over the SwissProt human main-isoform
proteome. The worker always logs which branch it took and the value it used.

## Output file

`final/annotations/aggregation_prone_regions.tsv`

## Output columns

| Column | Description |
|--------|-------------|
| `Protein_ID` | Gencode transcript name |
| `start` | 1-based inclusive start (matching `low_complexity.tsv`) |
| `end` | 1-based inclusive end |
| `length` | `end - start + 1` |
| `mean_a3v` | Mean **raw** a3v of the segment's residues |
| `peak_a3v` | Maximum **window** mean inside the segment — the quantity the threshold acts on |

## Notes

- No per-residue score track is written. The continuous window-mean profile is
  recomputable from `final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv` with the
  same scale and window, so storing ~11M rows in the base data would be pure duplication.
- The disorder/order class of a region is deliberately absent: the disorder call is a
  downstream combination rule (pLDDT / IUPred3 / MobiDB / Pfam) and baking one version of
  it into the annotation would freeze it into the base data.
- Related but distinct: [PLAAC](../phase_separation/plaac.md) scores prion-like
  composition, [catGRANULE](../phase_separation/catgranule.md) scores LLPS propensity,
  and [low-complexity regions](../structure/lcr.md) mark compositional bias. None of them
  is an aggregation call.
- Unlike most prediction tracks there is no "missing dependency" degradation path — the
  worker is pure Python, so it either runs or the pipeline fails loudly.
- Worker: `bin/create_apr_worker.py` (`APR_MAP`, `modules/structure.nf`).

## Running it

```bash
# as part of a pipeline run
nextflow run main.nf --project test_one_protein --data local --machine hard \
    --target_gene RAF1 --modules apr -resume

# backfill into an existing results dir, no Nextflow rerun
bin/backfill_tracks.sh results/discanvis --tracks apr
```
