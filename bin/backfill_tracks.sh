#!/usr/bin/env bash
# backfill_tracks.sh — run standalone annotation-track workers against an existing
# results dir's sequence table (no Nextflow rerun). Writes into <results_dir>/final/.
# Each track is routed to the interpreter that has its dependencies.
#
# Usage: bin/backfill_tracks.sh <results_dir> [--tracks lcr,dssp,catgranule,plaac,finches]
#
# Env overrides:
#   PYTHON            interpreter for lcr/dssp/catgranule/plaac  (default: python)
#   CATGRANULE_PYTHON python for the catGRANULE env  (default: /opt/anaconda3/envs/catgranule/bin/python)
#   CATGRANULE_LIB    catGRANULE 2.0 repo path        (default: /dlab/home/norbi/PycharmProjects/catGRANULE2.0)
#   FINCHES_PYTHON    python for the finches env      (default: /opt/anaconda3/envs/finches/bin/python)
#   PLAAC_JAR         PLAAC jar path (optional; worker default used if unset)
#   DISPROT_TSV       DisProt bulk TSV path           (default: references/disprot/disprot_regions.tsv)
set -euo pipefail

RESULTS_DIR="${1:-}"
if [[ -z "$RESULTS_DIR" ]]; then
    echo "usage: $0 <results_dir> [--tracks lcr,dssp,catgranule,plaac,finches]" >&2
    exit 2
fi
shift || true

TRACKS="lcr,dssp,catgranule,plaac"   # finches opt-in only
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tracks) TRACKS="$2"; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

PYTHON="${PYTHON:-python}"
CATGRANULE_PYTHON="${CATGRANULE_PYTHON:-/opt/anaconda3/envs/catgranule/bin/python}"
CATGRANULE_LIB="${CATGRANULE_LIB:-/dlab/home/norbi/PycharmProjects/catGRANULE2.0}"
FINCHES_PYTHON="${FINCHES_PYTHON:-/opt/anaconda3/envs/finches/bin/python}"
DISPROT_TSV="${DISPROT_TSV:-references/disprot/disprot_regions.tsv}"
MOBIDB_TSV="${MOBIDB_TSV:-references/mobidb/mobidb_human.tsv}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FINAL="$RESULTS_DIR/final"
SEQ="$FINAL/sequence/loc_chrom_with_names_isoforms_with_seq.tsv"
if [[ ! -f "$SEQ" ]]; then
    echo "sequence table not found: $SEQ" >&2
    exit 1
fi

run_track() {
    case "$1" in
        lcr)
            "$PYTHON" "$SCRIPT_DIR/create_lcr_worker.py" \
                --seq_table "$SEQ" --outdir "$FINAL/annotations" --only_main_isoforms ;;
        dssp)
            "$PYTHON" "$SCRIPT_DIR/create_dssp_worker.py" \
                --seq_table "$SEQ" --outdir "$FINAL/structure" --only_main_isoforms ;;
        catgranule)
            "$PYTHON" "$SCRIPT_DIR/create_catgranule_worker.py" \
                --seq_table "$SEQ" --outdir "$FINAL/phase_separation" --only_main_isoforms \
                --catgranule_python "$CATGRANULE_PYTHON" --catgranule_lib "$CATGRANULE_LIB" ;;
        plaac)
            "$PYTHON" "$SCRIPT_DIR/create_plaac_worker.py" \
                --seq_table "$SEQ" --outdir "$FINAL/phase_separation" --only_main_isoforms \
                ${PLAAC_JAR:+--plaac_jar "$PLAAC_JAR"} ;;
        finches)
            "$FINCHES_PYTHON" "$SCRIPT_DIR/create_finches_worker.py" \
                --loc_chrom "$SEQ" --output_dir "$FINAL/pathogenicity" --only_main_isoforms ;;
        disprot)
            [[ -f "$DISPROT_TSV" ]] || { echo "DisProt TSV not found: $DISPROT_TSV (set DISPROT_TSV)" >&2; return 1; }
            "$PYTHON" "$SCRIPT_DIR/create_disprot_worker.py" \
                --seq_table "$SEQ" --disprot_tsv "$DISPROT_TSV" \
                --outdir "$FINAL/disorder" --only_main_isoforms ;;
        mobidb)
            [[ -f "$MOBIDB_TSV" ]] || { echo "MobiDB TSV not found: $MOBIDB_TSV (set MOBIDB_TSV)" >&2; return 1; }
            "$PYTHON" "$SCRIPT_DIR/create_mobidb_worker.py" \
                --seq_table "$SEQ" --mobidb_tsv "$MOBIDB_TSV" \
                --outdir "$FINAL/disorder" ;;
        *) echo "unknown track: $1" >&2; return 1 ;;
    esac
}

IFS=',' read -ra SEL <<< "$TRACKS"
for t in "${SEL[@]}"; do
    echo ">>> backfilling track: $t"
    run_track "$t"
done
echo "backfill complete: $TRACKS -> $FINAL"
