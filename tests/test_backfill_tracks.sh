#!/usr/bin/env bash
# Smoke test: backfill_tracks.sh with a tiny fake results dir, lcr track only.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

mkdir -p "$TMP/final/sequence"
printf 'Protein_ID\tmain_isoform\tSequence\n' > "$TMP/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv"
printf 'LCRX-201\tyes\tMKVLAAGDEFRHIKPWY%s\n' "$(printf 'S%.0s' {1..25})" \
    >> "$TMP/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv"

PYTHON="${PYTHON:-python}" bash "$ROOT/bin/backfill_tracks.sh" "$TMP" --tracks lcr
test -f "$TMP/final/annotations/low_complexity.tsv" || { echo "no LCR output"; exit 1; }
echo "OK: backfill lcr produced low_complexity.tsv"

# --- apr track ---------------------------------------------------------------
printf 'APRX-201\tyes\tDEKRNQDEKRNQVILFVILFVILFDEKRNQDEKRNQ\n' \
    >> "$TMP/final/sequence/loc_chrom_with_names_isoforms_with_seq.tsv"

PYTHON="${PYTHON:-python}" bash "$ROOT/bin/backfill_tracks.sh" "$TMP" --tracks apr
APR="$TMP/final/annotations/aggregation_prone_regions.tsv"
test -f "$APR" || { echo "no APR output"; exit 1; }
head -1 "$APR" | grep -q $'Protein_ID\tstart\tend\tlength\tmean_a3v\tpeak_a3v' \
    || { echo "bad APR header"; exit 1; }
grep -q '^APRX-201' "$APR" || { echo "no APR called on the hydrophobic test protein"; exit 1; }
echo "OK: backfill apr produced aggregation_prone_regions.tsv"
