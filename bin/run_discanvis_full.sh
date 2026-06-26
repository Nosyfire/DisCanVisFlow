#!/usr/bin/env bash
# bin/run_discanvis_full.sh
#
# Runs the full-proteome discanvis pipeline and then derives the other four
# project result directories from it automatically.
#
# Designed to be run in the background:
#   nohup bash bin/run_discanvis_full.sh &
#   echo $! > logs/discanvis_full_run.pid
#
# Monitor progress:
#   tail -f logs/discanvis_full_run.log

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$REPO_DIR/logs/discanvis_full_run.log"

cd "$REPO_DIR"

# Activate conda environment
CONDA_SH="$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh"
# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate discanvis

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "======================================================"
log "discanvis full-proteome run START"
log "REPO:    $REPO_DIR"
log "CONDA:   $CONDA_DEFAULT_ENV"
log "NF ver:  $(nextflow -version 2>&1 | grep 'version' | head -1 | tr -s ' ')"
log "======================================================"

# ── Phase 1: Full pipeline ─────────────────────────────────────────────────────
log "Phase 1: nextflow run main.nf --project discanvis --data local --machine hard -resume"

nextflow run main.nf \
    --project discanvis \
    --data local \
    --machine hard \
    -resume 2>&1 | tee -a "$LOG"
NF_STATUS=${PIPESTATUS[0]}

log "Phase 1 EXIT STATUS: $NF_STATUS"

if [[ $NF_STATUS -ne 0 ]]; then
    log "Pipeline failed — aborting derivations."
    exit "$NF_STATUS"
fi

log "======================================================"
log "Phase 2: Deriving project results from discanvis"
log "======================================================"

bash bin/derive_projects_from_discanvis.sh 2>&1 | tee -a "$LOG"

log "======================================================"
log "ALL DONE"
log "Results:"
for d in results/discanvis results/vep_benchmarking results/cellular_vulnerability results/test_subset results/raf1_example; do
    if [[ -d "$d/final" ]]; then
        n=$(find "$d/final" -name "*.tsv" | wc -l)
        log "  $d/  ($n TSV files)"
    fi
done
log "======================================================"
