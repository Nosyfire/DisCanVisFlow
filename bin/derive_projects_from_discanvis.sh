#!/usr/bin/env bash
# bin/derive_projects_from_discanvis.sh
#
# Run this AFTER the discanvis full-proteome pipeline completes.
# Derives four additional project result directories from the discanvis output
# without re-running the pipeline:
#
#   results/vep_benchmarking/     — full-proteome copy (rsync from discanvis)
#   results/cellular_vulnerability/ — 797-gene extraction from discanvis
#   results/test_subset/          — 5-gene extraction (TP53, RAF1, BRAF, KRAS, EGFR)
#   results/raf1_example/         — single-gene extraction (RAF1)
#
# Usage:
#   bash bin/derive_projects_from_discanvis.sh
#   bash bin/derive_projects_from_discanvis.sh --source results/discanvis   # override source

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SOURCE="${1:-results/discanvis}"
GENE_LIST="config/gene_lists/cellular_vulnerability.txt"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

cd "$REPO_DIR"

if [[ ! -d "$SOURCE/final" ]]; then
    echo "ERROR: $SOURCE/final not found. Run the discanvis pipeline first." >&2
    exit 1
fi

# ── 1. vep_benchmarking — full proteome, rsync from discanvis ─────────────────
log "=== vep_benchmarking: full-proteome rsync from $SOURCE ==="
mkdir -p results/vep_benchmarking
rsync -a --info=progress2 "$SOURCE/" results/vep_benchmarking/
log "  Done: results/vep_benchmarking/"

# ── 2. cellular_vulnerability — selective copy (full proteome, subset of data types) ─
# Includes: annotations, sequence, drivers, dbnsfp, alphamissense,
#           combined disorder only, DepMap mutations, intermediate, mapping_reports.
# Excludes: genome, conservation, pdb, position, disease, IUPred/ANCHOR/AIUPred
#           disorder scores, ClinVar/TCGA/cBioPortal mutations, mavedb, proteingym.
log "=== cellular_vulnerability: selective copy from $SOURCE ==="
DST=results/cellular_vulnerability
rm -rf "$DST"
mkdir -p "$DST/final"

for d in annotations sequence drivers dbnsfp; do
    mkdir -p "$DST/final/$d"
    rsync -a "$SOURCE/final/$d/" "$DST/final/$d/"
done

mkdir -p "$DST/final/disorder"
cp "$SOURCE/final/disorder/CombinedDisorderNew.tsv"     "$DST/final/disorder/"
cp "$SOURCE/final/disorder/CombinedDisorderNew_Pos.tsv" "$DST/final/disorder/"

mkdir -p "$DST/final/pathogenicity"
cp "$SOURCE/final/pathogenicity/alphamissense.tsv" "$DST/final/pathogenicity/"

mkdir -p "$DST/final/mutations/DepMap"
rsync -a "$SOURCE/final/mutations/DepMap/" "$DST/final/mutations/DepMap/"

rsync -a "$SOURCE/intermediate/" "$DST/intermediate/"
rsync -a "$SOURCE/final/mapping_reports/" "$DST/final/mapping_reports/" 2>/dev/null || true

log "  Done: results/cellular_vulnerability/"

# ── 3. test_subset — 5-gene extraction ────────────────────────────────────────
log "=== test_subset: extracting TP53,RAF1,BRAF,KRAS,EGFR from $SOURCE ==="
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source "$SOURCE" \
    --gene TP53,RAF1,BRAF,KRAS,EGFR \
    --out results/test_subset
log "  Done: results/test_subset/"

# ── 4. raf1_example — single-gene extraction ──────────────────────────────────
log "=== raf1_example: extracting RAF1 from $SOURCE ==="
conda run -n discanvis python bin/extract_gene_from_results.py \
    --source "$SOURCE" \
    --gene RAF1 \
    --out results/raf1_example
log "  Done: results/raf1_example/"

log "=== All projects derived from $SOURCE ==="
echo ""
echo "Results:"
for d in results/vep_benchmarking results/cellular_vulnerability results/test_subset results/raf1_example; do
    if [[ -d "$d/final" ]]; then
        n=$(find "$d/final" -name "*.tsv" | wc -l)
        echo "  $d/  ($n TSV files)"
    fi
done
