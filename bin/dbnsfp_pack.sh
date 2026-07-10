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
command -v bgzip >/dev/null || { echo "bgzip not found (install htslib)" >&2; exit 2; }

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

echo "[$(date)] bgzip (block-gzip + random-access index) …"
bgzip -i -I "$OUTDIR/dbnsfp_scores.tsv.gz.gzi" -c "$SORTED" > "$OUTDIR/dbnsfp_scores.tsv.gz"
rm -f "$SORTED"
rmdir "$TMPDIR" 2>/dev/null || true

echo "[$(date)] DONE_PACK"
ls -la "$OUTDIR"/dbnsfp_scores.tsv.gz "$OUTDIR"/dbnsfp_scores.pidx "$OUTDIR"/dbnsfp_scores.header
