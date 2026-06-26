#!/usr/bin/env bash
# bin/refresh_refs.sh — force re-download of specific cached references.
#
# Usage:
#   bin/refresh_refs.sh                    # list all cached sources + sizes
#   bin/refresh_refs.sh clinvar            # refresh ClinVar only
#   bin/refresh_refs.sh clinvar mobidb go  # refresh multiple sources
#   bin/refresh_refs.sh all                # refresh everything (keeps hg38, dbsnp, alphafold by default)
#   bin/refresh_refs.sh --force all        # refresh truly everything
#
# After running, re-execute the pipeline with -resume — FETCH_* processes will
# re-download only the deleted files; everything else remains cached.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="$(cd "$SCRIPT_DIR/.." && pwd)/references"

# Sources and the files/dirs they own in references/
declare -A SOURCES=(
    [uniprot]="uniprot/uniprot_swissprot.fasta uniprot/UP000005640_9606_additional.fasta uniprot/uniprot_sprot.dat.gz uniprot_parsed/"
    [gencode]="gencode/"
    [clinvar]="clinvar/clinvar_grch38.vcf.gz"
    [mobidb]="mobidb/mobidb_human.tsv"
    [go]="go/goa_human.gaf.gz go/go.obo"
    [mondo]="mondo/mondo.obo"
    [alphamissense]="alphamissense/AlphaMissense_isoforms_hg38.tsv.gz alphamissense_parsed/"
    [ppi]="ppi/intact_human.mitab.zip ppi/biogrid_human.mitab.zip ppi/hippie_current.txt ppi/processed/"
    [sifts]="sifts/uniprot_segments_observed.tsv.gz"
    [mavedb]="mavedb/mave_single_mutant_protein.tsv"
    [proteingym]="proteingym/proteingym_substitutions.zip"
    [depmap]="depmap/OmicsSomaticMutationsProfile.csv depmap/OmicsSomaticMutations.csv depmap/depmap_mutations_raw.tsv"
    [omim]="omim/humsavar.txt omim/omim_disease.tsv"
    [interpro]="interpro/interpro_pfam.tsv"
    # Large / rarely changed — not included in 'all' unless --force
    [hg38]="hg38/hg38.2bit"
    [dbsnp]="dbsnp/dbSnp155Common.bb"
    [alphafold]="alphafold/"
)

# Sources skipped by 'all' unless --force is passed (too large / slow to re-download)
LARGE_SOURCES=(hg38 dbsnp alphafold)

usage() {
    echo "Usage: $0 [--force] [source...] | all | (no args = list)"
    echo ""
    echo "Available sources:"
    for src in $(echo "${!SOURCES[@]}" | tr ' ' '\n' | sort); do
        printf "  %-18s  %s\n" "$src" "${SOURCES[$src]}"
    done
    echo ""
    echo "  all   — refresh all sources except: ${LARGE_SOURCES[*]}"
    echo "  --force all — refresh everything"
}

list_cached() {
    echo "Cached references in $REF_DIR:"
    echo ""
    for src in $(echo "${!SOURCES[@]}" | tr ' ' '\n' | sort); do
        local found=0
        for rel in ${SOURCES[$src]}; do
            p="$REF_DIR/$rel"
            if [[ -e "$p" ]]; then
                sz=$(du -sh "$p" 2>/dev/null | cut -f1)
                ts=$(date -r "$p" '+%Y-%m-%d' 2>/dev/null || echo "?")
                printf "  %-18s  %-8s  %s  %s\n" "$src" "$sz" "$ts" "$rel"
                found=1
            fi
        done
        [[ $found -eq 0 ]] && printf "  %-18s  (not cached)\n" "$src"
    done
}

delete_source() {
    local src="$1"
    if [[ -z "${SOURCES[$src]+_}" ]]; then
        echo "ERROR: unknown source '$src'. Run without args to see available sources." >&2
        exit 1
    fi
    local deleted=0
    for rel in ${SOURCES[$src]}; do
        p="$REF_DIR/$rel"
        if [[ -e "$p" ]]; then
            echo "  deleting: $p"
            rm -rf "$p"
            deleted=1
        fi
    done
    [[ $deleted -eq 0 ]] && echo "  $src: nothing cached to delete"
}

# ── Main ─────────────────────────────────────────────────────────────────────

FORCE=0
TARGETS=()

for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) usage; exit 0 ;;
        *) TARGETS+=("$arg") ;;
    esac
done

if [[ ${#TARGETS[@]} -eq 0 ]]; then
    list_cached
    exit 0
fi

if [[ ${#TARGETS[@]} -eq 1 && "${TARGETS[0]}" == "all" ]]; then
    echo "Refreshing all sources (use --force all to include hg38, dbsnp, alphafold)..."
    for src in $(echo "${!SOURCES[@]}" | tr ' ' '\n' | sort); do
        skip=0
        if [[ $FORCE -eq 0 ]]; then
            for large in "${LARGE_SOURCES[@]}"; do
                [[ "$src" == "$large" ]] && skip=1 && break
            done
        fi
        [[ $skip -eq 1 ]] && echo "  skipping $src (large; use --force to include)" && continue
        delete_source "$src"
    done
else
    for src in "${TARGETS[@]}"; do
        delete_source "$src"
    done
fi

echo ""
echo "Done. Re-run your pipeline with -resume to re-download only deleted files."
