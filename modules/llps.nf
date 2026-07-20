/*
 * modules/llps.nf — LLPS database integration (PhaSepDB, LLPSDB, DisPhaseDB)
 *
 * FETCH_LLPS_SOURCES   downloads the three databases (each guarded so one dead
 *                      endpoint never fails the run) and normalises whatever
 *                      arrived into three flat TSVs via parse_llps_sources.py.
 *                      storeDir-cached in references/llps/.
 * LLPS_REGIONS_MAP     regions + protein-level info → final/annotations/
 * LLPS_VARIANTS_MAP    variants (residue-validated) → final/mutations/
 *
 * Confirmed bulk endpoints (2026-07):
 *   PhaSepDB  http://db2.phasep.pro/static/db/database/phasepdbv2_1_llps.xlsx
 *   LLPSDB    http://bio-comp.org.cn/llpsdbv2/download/Phase_separation_Unambiguous.zip
 *   DisPhaseDB  disphasedb.leloir.org.ar — API-only SPA, frequently offline;
 *               supply a dump via --disphasedb_path (or --disphasedb_url).
 */

process FETCH_LLPS_SOURCES {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/llps" : "${params.ref_dir}/llps" }

    output:
    path 'normalized', type: 'dir', emit: norm

    script:
    def do_phasepdb  = params.skip_phasepdb  ? false : true
    def do_llpsdb    = params.skip_llpsdb    ? false : true
    def do_disphase  = params.skip_disphasedb ? false : true
    def phasepdb_url = params.phasepdb_url ?: 'http://db2.phasep.pro/static/db/database/phasepdbv2_1_llps.xlsx'
    def llpsdb_url   = params.llpsdb_url   ?: 'http://bio-comp.org.cn/llpsdbv2/download/Phase_separation_Unambiguous.zip'
    def disphase_url = params.disphasedb_url ?: ''
    def disphase_path = params.disphasedb_path ?: ''
    """
    set +e
    mkdir -p raw normalized
    PARSE_ARGS=""

    # --- PhaSepDB -------------------------------------------------------------
    if [ "${do_phasepdb}" = "true" ]; then
        curl -fsSL --max-time 180 -o raw/phasepdb.xlsx '${phasepdb_url}' \\
            && PARSE_ARGS="\$PARSE_ARGS --phasepdb raw/phasepdb.xlsx" \\
            || echo "[WARN] PhaSepDB download failed — skipping this source"
    fi

    # --- LLPSDB ---------------------------------------------------------------
    if [ "${do_llpsdb}" = "true" ]; then
        if curl -fsSL --max-time 180 -o raw/llpsdb.zip '${llpsdb_url}'; then
            ( cd raw && unzip -o -q llpsdb.zip )
            PROT=\$(find raw -iname 'protein.xls' | head -1)
            [ -n "\$PROT" ] && PARSE_ARGS="\$PARSE_ARGS --llpsdb_protein \$PROT" \\
                            || echo "[WARN] LLPSDB protein.xls not found in archive"
        else
            echo "[WARN] LLPSDB download failed — skipping this source"
        fi
    fi

    # --- DisPhaseDB (manual path preferred; server usually offline) -----------
    if [ "${do_disphase}" = "true" ]; then
        if [ -n "${disphase_path}" ] && [ -f "${disphase_path}" ]; then
            cp "${disphase_path}" raw/disphasedb_dump
            PARSE_ARGS="\$PARSE_ARGS --disphasedb raw/disphasedb_dump"
        elif [ -n "${disphase_url}" ]; then
            curl -fsSL --max-time 300 -o raw/disphasedb_dump '${disphase_url}' \\
                && PARSE_ARGS="\$PARSE_ARGS --disphasedb raw/disphasedb_dump" \\
                || echo "[WARN] DisPhaseDB download failed — supply --disphasedb_path"
        else
            echo "[INFO] DisPhaseDB not fetched (server offline; no --disphasedb_path/url given)"
        fi
    fi

    set -e
    parse_llps_sources.py \$PARSE_ARGS --outdir normalized
    echo "Normalised LLPS sources:"; wc -l normalized/*.tsv
    """

    stub:
    """
    mkdir -p normalized
    printf 'acc\\tstart\\tend\\tregion_type\\tsource\\n'                       > normalized/llps_regions.tsv
    printf 'acc\\tgene\\torganelle_mlo\\tmaterial_state\\tklass\\tsource\\n'   > normalized/llps_proteins.tsv
    printf 'acc\\tpos\\twt_aa\\tmut_aa\\tvariant_class\\tdisease_term\\tsource\\n' > normalized/llps_variants.tsv
    """
}


process LLPS_REGIONS_MAP {
    tag { "llps_regions_map" }
    label 'process_medium'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path norm

    output:
    path "llps_regions.tsv",  emit: regions
    path "llps_proteins.tsv", emit: proteins

    script:
    """
    create_llps_regions_worker.py \\
        --regions  ${norm}/llps_regions.tsv \\
        --proteins ${norm}/llps_proteins.tsv \\
        --loc_chrom ${loc_chrom} \\
        --outdir   .
    """

    stub:
    """
    printf 'Protein_ID\\tGene\\tstart\\tend\\tregion_type\\tsource\\tmapping_type\\n' > llps_regions.tsv
    printf 'Protein_ID\\tGene\\torganelle_mlo\\tmaterial_state\\tklass\\tsource\\tmapping_type\\n' > llps_proteins.tsv
    """
}


process LLPS_VARIANTS_MAP {
    tag { "llps_variants_map" }
    label 'process_medium'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/mutations"
                                : "${params.outdir}/final/mutations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path norm

    output:
    path "llps_variants.tsv", emit: variants

    script:
    """
    create_llps_variants_worker.py \\
        --variants  ${norm}/llps_variants.tsv \\
        --loc_chrom ${loc_chrom} \\
        --outdir    .
    """

    stub:
    """
    printf 'Protein_ID\\tProtein_position\\tWT_AA\\tMut_AA\\tvariant_class\\tdisease_term\\tsource_db\\tmapping_type\\n' > llps_variants.tsv
    """
}
