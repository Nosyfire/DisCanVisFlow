#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

mkdir -p "$TMP/src/final/annotations" "$TMP/src/final/structure" \
         "$TMP/src/final/phase_separation" "$TMP/src/final/pathogenicity"
echo x > "$TMP/src/final/annotations/low_complexity.tsv"
echo x > "$TMP/src/final/structure/dssp.tsv"
echo x > "$TMP/src/final/phase_separation/catgranule.tsv"
# leave plaac + finches absent → must warn+skip, not fail
bash "$ROOT/bin/copy_tracks_to_benchmark.sh" "$TMP/src" "$TMP/dst"
test -f "$TMP/dst/final/annotations/low_complexity.tsv"
test -f "$TMP/dst/final/structure/dssp.tsv"
test -f "$TMP/dst/final/phase_separation/catgranule.tsv"
test ! -f "$TMP/dst/final/phase_separation/plaac.tsv"
echo "OK: copy driver copied present files and skipped absent ones"
