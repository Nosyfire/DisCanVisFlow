#!/usr/bin/env bash
# copy_tracks_to_benchmark.sh — copy the five generated track files from a source
# results dir into a destination results dir at identical relative paths.
# Protein sets are identical between discanvis and vep_benchmarking, so no filter.
#
# Usage: bin/copy_tracks_to_benchmark.sh <src_results_dir> <dst_results_dir>
set -euo pipefail

SRC="${1:-}"; DST="${2:-}"
if [[ -z "$SRC" || -z "$DST" ]]; then
    echo "usage: $0 <src_results_dir> <dst_results_dir>" >&2
    exit 2
fi

RELPATHS=(
    "final/annotations/low_complexity.tsv"
    "final/structure/dssp.tsv"
    "final/phase_separation/catgranule.tsv"
    "final/phase_separation/plaac.tsv"
    "final/pathogenicity/finches_saturation.tsv"
)

copied=0
for rel in "${RELPATHS[@]}"; do
    s="$SRC/$rel"
    d="$DST/$rel"
    if [[ ! -f "$s" ]]; then
        echo "WARN: source missing, skipping: $s" >&2
        continue
    fi
    mkdir -p "$(dirname "$d")"
    cp -f "$s" "$d"
    echo "copied: $rel"
    copied=$((copied + 1))
done
echo "done: copied $copied/${#RELPATHS[@]} track files -> $DST"
