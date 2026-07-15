#!/usr/bin/env bash
# dbnsfp_pack.sh — turn a large plain dbNSFP map TSV into a compact, slice-able
# bundle:
#   <outdir>/dbnsfp_scores.tsv.gz       BGZF (block-gzip) body, sorted by
#                                        Protein_ID then Protein_position
#   <outdir>/dbnsfp_scores.tsv.gz.gzi   bgzip random-access index
#   <outdir>/dbnsfp_scores.pidx         Protein_ID <tab> offset <tab> length
#   <outdir>/dbnsfp_scores.header       the single header line
#
# Slice one protein without reading the whole file:
#   bin/slice_dbnsfp.py --bgz <outdir>/dbnsfp_scores.tsv.gz --id RAF1-201
#
# Usage: bin/dbnsfp_pack.sh <input_map.tsv> <outdir>
# Env: SORT_BUF (default 800G), SORT_PAR (default 32), TMPDIR (default <outdir>/sorttmp)
set -euo pipefail

IN="${1:?usage: dbnsfp_pack.sh <input_map.tsv> <outdir>}"
OUTDIR="${2:?usage: dbnsfp_pack.sh <input_map.tsv> <outdir>}"
# $BGZIP overrides; require htslib >= 1.11 (older `bgzip -b` overflows >2GB offsets)
BGZIP="${BGZIP:-bgzip}"
command -v "$BGZIP" >/dev/null || { echo "bgzip not found (set \$BGZIP or install htslib)" >&2; exit 2; }
BGZ_VER="$("$BGZIP" --version | head -1 | grep -oE '[0-9]+\.[0-9]+' | head -1)"
BGZ_MAJ="${BGZ_VER%%.*}"; BGZ_MIN="${BGZ_VER##*.}"
if [ "$BGZ_MAJ" -lt 1 ] || { [ "$BGZ_MAJ" -eq 1 ] && [ "$BGZ_MIN" -lt 11 ]; }; then
  echo "WARN: bgzip $BGZ_VER is old — 'bgzip -b' overflows offsets > 2 GB. A full-proteome bundle needs htslib >= 1.11 (set \$BGZIP)." >&2
fi

mkdir -p "$OUTDIR"
export TMPDIR="${TMPDIR:-$OUTDIR/sorttmp}"; mkdir -p "$TMPDIR"
HDR="$OUTDIR/dbnsfp_scores.header"
SORTED="$OUTDIR/dbnsfp_scores.sorted.tsv"
PIDX="$OUTDIR/dbnsfp_scores.pidx"

echo "[$(date)] header + sort by Protein_ID,Protein_position …"
head -1 "$IN" > "$HDR"
# --compress-program keeps sort's spill files small on huge inputs
SORT_COMPRESS=""
command -v pigz >/dev/null && SORT_COMPRESS="--compress-program=pigz"
tail -n +2 "$IN" | LC_ALL=C sort -t$'\t' -k1,1 -k2,2n \
    -S "${SORT_BUF:-800G}" --parallel="${SORT_PAR:-32}" $SORT_COMPRESS > "$SORTED"

echo "[$(date)] building Protein_ID -> (offset,length) index …"
# offsets are into the UNCOMPRESSED sorted body (what bgzip -b addresses)
awk -F'\t' 'BEGIN{OFS="\t"; off=0}
{
  len=length($0)+1
  if($1!=cur){ if(cur!=""){print cur, start, off-start} cur=$1; start=off }
  off+=len
}
END{ if(cur!="") print cur, start, off-start }' "$SORTED" > "$PIDX"
echo "[$(date)] indexed $(wc -l < "$PIDX") proteins"

echo "[$(date)] bgzip (block-gzip) …"
"$BGZIP" -@ "${BGZIP_THREADS:-8}" -c "$SORTED" > "$OUTDIR/dbnsfp_scores.tsv.gz"
# build the random-access index explicitly: `-i -c` (stdout) leaves an empty .gzi
# on some htslib builds, whereas `bgzip -r` always reindexes the finished file.
echo "[$(date)] bgzip -r (random-access .gzi index) …"
"$BGZIP" -r "$OUTDIR/dbnsfp_scores.tsv.gz"
rm -f "$SORTED"
rmdir "$TMPDIR" 2>/dev/null || true

echo "[$(date)] DONE_PACK"
ls -la "$OUTDIR"/dbnsfp_scores.tsv.gz "$OUTDIR"/dbnsfp_scores.pidx "$OUTDIR"/dbnsfp_scores.header
