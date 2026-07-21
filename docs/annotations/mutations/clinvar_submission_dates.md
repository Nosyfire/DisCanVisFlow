# ClinVar Submission Dates

## Description

When was each ClinVar variant last evaluated by a submitter? The ClinVar VCF
carries no dates at all, so they come from `submission_summary.txt.gz` — one row
per submitted record (SCV) — aggregated to one row per variant.

## Output file

`final/mutations/ClinVar/clinvar_submission_dates.tsv`

| Column | Description |
|--------|-------------|
| `VariationID` | ClinVar VariationID (`ncbi.nlm.nih.gov/clinvar/<id>`) |
| `Mutation` | Genomic HGVS (`CLNHGVS`) — **the join key** |
| `Chromosome`, `Position`, `Ref`, `Alt` | VCF coordinates |
| `Gene` | Gene symbol from `GENEINFO` |
| `last_evaluated` | `max(DateLastEvaluated)` over all SCVs, ISO `YYYY-MM-DD` |
| `first_evaluated` | `min(DateLastEvaluated)` over all SCVs |
| `n_submissions` | Number of SCV records, **including undated ones** |
| `last_submitter` | Submitter of the SCV holding `last_evaluated` |
| `last_clinical_significance` | That SCV's germline classification |
| `last_review_status` | That SCV's review status |

## Joining

The table is **variant-keyed, not Protein_ID-keyed** — deliberately. Isoform
expansion means one variant appears on many `Protein_ID`s, so a per-isoform date
column would duplicate the same date millions of times. Instead join on
`Mutation`, which is byte-identical to the `Mutation` column of
`Missense/Indel/Nonsense/Frameshift_filter_mutations_mapped.tsv`:

```sql
SELECT m.Protein_ID, m.Protein_position, d.last_evaluated, d.n_submissions
FROM   missense_mapped m
JOIN   clinvar_submission_dates d USING (Mutation);
```

Coverage is complete: every distinct `Mutation` in the mapped ClinVar TSVs
resolves (verified 91,891 / 91,891 on a sample of the full-proteome run).

## Notes

- **Dateless submissions.** ClinVar writes a bare `-` when a submitter gave no
  evaluation date. Those rows still increment `n_submissions` but never move
  `last_evaluated` / `first_evaluated` — so `n_submissions` can exceed the number
  of dated records, and a variant can have submissions but no date at all.
- **Variants with no submissions** are kept with empty dates and
  `n_submissions=0`, so the table is a complete companion to the VCF rather than
  a filtered subset.
- **Release skew is harmless.** `submission_summary.txt.gz` is always the current
  NCBI release while `--clinvar_vcf` may be pinned to an older one. Because the
  join is on the stable `VariationID` (not on coordinates), submissions for
  variants absent from the VCF are simply ignored. Supply a matching file with
  `--clinvar_submission_summary` if you need exact release lock-step.
- **`last_evaluated` is the submitter's evaluation date, not the upload date.**
  ClinVar publishes no per-SCV upload timestamp; the evaluation date is the
  closest available proxy and is what the ClinVar web UI shows.

## Scale (2026-07 release, full proteome)

| Metric | Value |
|--------|-------|
| Variants in VCF | 4,397,693 |
| Submission rows read | 6,376,003 |
| Variants with a submission date | 4,124,269 (93.8 %) |
| Output size | 716 MB |
| Runtime / peak RSS | 79 s / 2.2 GB |
| Date range | 1965-01-01 → 2026-07-14 |

## Running

Enabled by default whenever ClinVar mutations run. To control it:

```bash
# skip it
nextflow run main.nf ... --skip_clinvar_dates true

# only this module
nextflow run main.nf ... --modules clinvar_dates

# supply a local submission_summary instead of downloading
nextflow run main.nf ... --clinvar_submission_summary /path/to/submission_summary.txt.gz
```

Worker: `bin/create_clinvar_dates_worker.py`.
Processes: `FETCH_CLINVAR_SUBMISSIONS` + `CLINVAR_DATES` in `modules/mutation_mapping.nf`.
