#!/usr/bin/env bash
#
# generate_docs.sh — regenerate the mapping-report Markdown for a results dir
# WITHOUT running the Nextflow pipeline. "Docs-only" mode.
#
# It drives bin/create_mapping_report_worker.py against an existing
# results/<project>/final/ tree, producing the enriched reports:
#   * mapping_summary.md   (reproducibility + data-source versions + input scale
#                           + annotation coverage)
#   * mapping_coverage.tsv (full-proteome runs) OR <GENE>_mapping_report.md
#                           (small runs / gene slices)
#
# Use it to:
#   * refresh reports after extracting a gene slice (extract_gene_from_results.py)
#   * backfill enriched provenance/base-stats onto an old run with no recompute
#
# Usage:
#   bin/generate_docs.sh <results_dir> [options]
#
#   <results_dir>   e.g. results/discanvis  (must contain final/sequence/)
#
# Options (all optional — sensible defaults are auto-detected):
#   --gencode  FASTA     GENCODE translations FASTA (entry count + version)
#   --uniprot  FASTA     UniProt SwissProt canonical FASTA (entry count)
#   --isoform  FASTA     UniProt curated-isoform FASTA (entry count)
#   --mapping-mode MODE  default: all_isoform_mapping
#   --outdir   DIR       default: <results_dir>/mapping_reports
#   --no-refs            skip reference FASTA counts (versions table omitted)
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORKER="${SCRIPT_DIR}/create_mapping_report_worker.py"

# ── parse args ───────────────────────────────────────────────────────────────
RESULTS_DIR="${1:-}"
if [[ -z "${RESULTS_DIR}" || "${RESULTS_DIR}" == -* ]]; then
    grep -E '^#( |$)' "$0" | sed -E 's/^# ?//' | head -40
    exit 1
fi
shift || true

GENCODE=""; UNIPROT=""; ISOFORM=""; MAPPING_MODE="all_isoform_mapping"
OUTDIR=""; USE_REFS=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gencode)      GENCODE="$2"; shift 2 ;;
        --uniprot)      UNIPROT="$2"; shift 2 ;;
        --isoform)      ISOFORM="$2"; shift 2 ;;
        --mapping-mode) MAPPING_MODE="$2"; shift 2 ;;
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --no-refs)      USE_REFS=0; shift ;;
        *) echo "[ERR] unknown option: $1" >&2; exit 1 ;;
    esac
done

RESULTS_DIR="$(cd "${RESULTS_DIR}" && pwd)"
# Relative label (repo-root-relative) for the publishable report — no abs paths.
REL_RESULTS="${RESULTS_DIR#${REPO_ROOT}/}"
FINAL_DIR="${RESULTS_DIR}/final"
INTER_DIR="${RESULTS_DIR}/intermediate"
[[ -d "${FINAL_DIR}/sequence" ]] || { echo "[ERR] no final/sequence/ under ${RESULTS_DIR}" >&2; exit 1; }
[[ -n "${OUTDIR}" ]] || OUTDIR="${RESULTS_DIR}/mapping_reports"

# ── locate the sequence table (pipeline's choice first, then fallbacks) ───────
SEQ_TABLE=""
for cand in loc_chrom_with_names_isoforms_with_seq.tsv \
            loc_chrom_with_names.tsv \
            loc_chrom_with_names_main_isoform.tsv; do
    if [[ -f "${FINAL_DIR}/sequence/${cand}" ]]; then
        SEQ_TABLE="${FINAL_DIR}/sequence/${cand}"; break
    fi
done
[[ -n "${SEQ_TABLE}" ]] || { echo "[ERR] no sequence table found in ${FINAL_DIR}/sequence/" >&2; exit 1; }

# ── auto-detect reference FASTAs (for entry counts + versions) ────────────────
if [[ "${USE_REFS}" -eq 1 ]]; then
    if [[ -z "${GENCODE}" ]]; then
        GENCODE="$(grep -oE "gencode_fasta[^=]*=[^']*'[^']*'" "${REPO_ROOT}/config/data/local.config" 2>/dev/null \
                    | grep -oE "/[^']*" | head -1 || true)"
    fi
    if [[ -z "${UNIPROT}" ]]; then
        for c in "${REPO_ROOT}/references/uniprot/uniprot_swissprot.fasta" \
                 "${REPO_ROOT}/references/uniprot/UP000005640_9606.fasta"; do
            [[ -r "$c" ]] && { UNIPROT="$c"; break; }
        done
    fi
    if [[ -z "${ISOFORM}" ]]; then
        c="${REPO_ROOT}/references/uniprot/UP000005640_9606_additional.fasta"
        [[ -r "$c" ]] && ISOFORM="$c"
    fi
fi

# ── pipeline version from git (best-effort) ───────────────────────────────────
PIPE_VER="$(git -C "${REPO_ROOT}" describe --tags --always 2>/dev/null || echo 'unknown')"

echo "=============================================================="
echo " Docs-only report generation"
echo "=============================================================="
echo " Results dir : ${RESULTS_DIR}"
echo " Seq table   : ${SEQ_TABLE##*/}"
echo " Output dir  : ${OUTDIR}"
echo " GENCODE fa  : ${GENCODE:-<none>}"
echo " UniProt fa  : ${UNIPROT:-<none>}"
echo " Isoform fa  : ${ISOFORM:-<none>}"
echo " Mapping mode: ${MAPPING_MODE}"
echo

ref_args=()
[[ -n "${GENCODE}" ]] && ref_args+=(--gencode_fasta "${GENCODE}")
[[ -n "${UNIPROT}" ]] && ref_args+=(--uniprot_fasta "${UNIPROT}")
[[ -n "${ISOFORM}" ]] && ref_args+=(--uniprot_isoform_fasta "${ISOFORM}")

inter_args=()
[[ -d "${INTER_DIR}" ]] && inter_args+=(--intermediate_dir "${INTER_DIR}")

mkdir -p "${OUTDIR}"
# Remove stale report artifacts so a small run doesn't leave behind a previous
# full-proteome mapping_coverage.tsv (and vice-versa). Only report outputs are
# touched — never source data.
rm -f "${OUTDIR}"/mapping_summary.md \
      "${OUTDIR}"/mapping_coverage.tsv \
      "${OUTDIR}"/*_mapping_report.md 2>/dev/null || true

python3 "${WORKER}" \
    --seq_table  "${SEQ_TABLE}" \
    --final_dir  "${FINAL_DIR}" \
    "${inter_args[@]}" \
    --outdir     "${OUTDIR}" \
    --mapping_mode "${MAPPING_MODE}" \
    --pipeline_version "${PIPE_VER}" \
    --command    "bin/generate_docs.sh ${REL_RESULTS} (docs-only, no pipeline run)" \
    "${ref_args[@]}"

echo
echo "[ok] Reports written to ${OUTDIR}"
ls -1 "${OUTDIR}"
