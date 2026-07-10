/*
 * modules/pathogenicity.nf — Module 8 pathogenicity — AlphaMissense, dbNSFP, DepMap, MaveDB, ProteinGym, FINCHES
 *
 * Processes: PROTEINGYM_MAP, ALPHAMISSENSE_MAP, DEPMAP_MAP, DBNSFP_MAP, PATHOGENICITY_MAP, MAVEDB_MAP, FINCHES_MAP
 * (split out of the former annotation_mapping.nf monolith)
 */


// ──────────────────────────────────────────────────────────────────────────
// PROTEINGYM_MAP  — Module 8g: ProteinGym DMS scores + pathogenicity proxy
// ──────────────────────────────────────────────────────────────────────────
process PROTEINGYM_MAP {
    tag  { "proteingym" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom_seq
    path proteingym_tsv   // premapped: Protein_ID-keyed table; uniprot: raw UniProt-keyed table; or NO_FILE

    output:
    path "proteingym.tsv", emit: proteingym

    script:
    // uniprot mode = a freshly fetched UniProt-keyed raw table (FETCH_PROTEINGYM
    // or --proteingym_raw); otherwise filter the pre-mapped Protein_ID-keyed table.
    def uniprot_mode = (params.fetch_proteingym || params.proteingym_raw) as boolean
    def src_arg = uniprot_mode ? "--mapping_mode uniprot --proteingym_raw ${proteingym_tsv}"
                               : "--proteingym ${proteingym_tsv}"
    """
    create_proteingym_worker.py \\
        --seq_table   ${loc_chrom_seq} \\
        ${src_arg} \\
        --outdir      .
    """

    stub:
    """
    printf 'Protein_ID\\tProtein_position\\tprotein_variant\\tDMS_score\\tDMS_score_bin\\tDMS_id\\tuniprot_id\\tmapping_type\\n' > proteingym.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// ALPHAMISSENSE_MAP  — Module 8d: AlphaMissense GENCODE isoform scores
// ──────────────────────────────────────────────────────────────────────────
process ALPHAMISSENSE_MAP {
    tag  { "alphamissense" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path alphamissense_gz   // AlphaMissense_isoforms_hg38.tsv (decompressed by DECOMPRESS_ALPHAMISSENSE) or NO_FILE

    output:
    path "alphamissense.tsv", emit: alphamissense

    script:
    """
    create_alphamissense_worker.py \\
        --seq_table        ${loc_chrom} \\
        --alphamissense_gz ${alphamissense_gz} \\
        --outdir           .
    """

    stub:
    """
    echo -e "Protein_ID\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class" > alphamissense.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// DEPMAP_MAP  — Module 8e: DepMap cancer cell line somatic mutations
// ──────────────────────────────────────────────────────────────────────────
process DEPMAP_MAP {
    tag  { "depmap" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/mutations/DepMap"
                                : "${params.outdir}/final/mutations/DepMap" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path depmap_tsv   // mapped_filtered_mutations.tsv from DepMap or NO_FILE

    output:
    path "depmap_mutations.tsv", emit: depmap

    script:
    """
    create_depmap_worker.py \\
        --seq_table ${loc_chrom} \\
        --depmap_tsv ${depmap_tsv} \\
        --outdir    .
    """

    stub:
    """
    echo -e "Protein_ID\tHGVSp_Short\tModelID" > depmap_mutations.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// DBNSFP_MAP  — Module 8f: map raw dbNSFP (merged 5.x .gz or per-chr chr*.gz) via combined_map.map
// ──────────────────────────────────────────────────────────────────────────
process DBNSFP_MAP {
    tag  { "dbnsfp_map" }
    label 'process_high'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path combined_map
    path dbnsfp_raw_dir

    output:
    path "dbnsfp_scores.tsv", emit: scores

    script:
    def hdr_arg = params.dbnsfp_bed_header ? "--dbnsfp_bed_header ${params.dbnsfp_bed_header}" : ""
    def n_cpu   = task.cpus ?: 1
    """
    create_dbnsfp_map_worker.py \\
        --seq_table       ${loc_chrom} \\
        --combined_map    ${combined_map} \\
        --dbnsfp_raw_dir  ${dbnsfp_raw_dir} \\
        ${hdr_arg} \\
        --n_cpu           ${n_cpu} \\
        --outdir          .
    """

    stub:
    """
    echo -e "Protein_ID\tProtein_position\tAlphaMissense_score" > dbnsfp_scores.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// PATHOGENICITY_MAP  — Module 8f: dbNSFP pathogenicity predictor scores
// ──────────────────────────────────────────────────────────────────────────
process PATHOGENICITY_MAP {
    tag  { "pathogenicity" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path dbnsfp_tsv   // dbNSFP_custom/mapped_filtered_mutations.tsv or NO_FILE

    output:
    path "pathogenicity_scores.tsv", emit: scores

    script:
    """
    create_pathogenicity_worker.py \\
        --seq_table  ${loc_chrom} \\
        --dbnsfp_tsv ${dbnsfp_tsv} \\
        --outdir     .
    """

    stub:
    """
    echo -e "Protein_ID\tProtein_position\tAlphaMissense_score\tCADD_phred\tSIFT_score" > pathogenicity_scores.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// MAVEDB_MAP — MaveDB single-mutant functional scores (already Protein_ID-keyed)
// Filters the large pre-mapped MaveDB table to this run's isoforms.
// Python worker: create_mavedb_worker.py
// ──────────────────────────────────────────────────────────────────────────
process MAVEDB_MAP {
    tag  { "mavedb" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom_seq
    path mavedb_tsv      // premapped: Protein_ID-keyed table; uniprot: raw UniProt-keyed table; or NO_FILE

    output:
    path "mavedb.tsv", emit: mavedb

    script:
    // uniprot mode = a freshly fetched UniProt-keyed raw table (FETCH_MAVEDB or
    // --mavedb_raw); otherwise filter the pre-mapped Protein_ID-keyed table.
    def uniprot_mode = (params.fetch_mavedb || params.mavedb_raw) as boolean
    def src_arg = uniprot_mode ? "--mapping_mode uniprot --mavedb_raw ${mavedb_tsv}"
                               : "--mavedb ${mavedb_tsv}"
    """
    create_mavedb_worker.py \\
        --seq_table ${loc_chrom_seq} \\
        ${src_arg} \\
        --outdir    .
    """

    stub:
    """
    printf 'Protein_ID\\tProtein_position\\tprot_expr\\tscore\\tmavedb_id\\turn\\tgene_name\\tuniprot\\tTranscript_ID\\tis_double_mutant\\tmapping_type\\n' > mavedb.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// FINCHES_MAP  — Module 8h: site-wise Δε (LLPS-change score) via saturation
// mutagenesis using the FINCHES Mpipi_GGv1 force field.
//
// For every possible single-AA substitution at every position, computes:
//   Δε = epsilon(mut,mut) − epsilon(wt,wt)  (homotypic self-interaction energy)
// Positive Δε = mutation increases LLPS propensity.
//
// Output goes to final/pathogenicity/ alongside AlphaMissense / ProteinGym.
// ⚠  CC BY-NC 4.0 — non-commercial use only.
//
// Citation: Ginell et al. bioRxiv 2024.06.03.597104
// ──────────────────────────────────────────────────────────────────────────
process FINCHES_MAP {
    tag  { "finches_map" }
    label 'process_high'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom

    output:
    path "finches_saturation.tsv", emit: finches

    script:
    def finches_lib_arg = params.finches_lib         ? "--finches_lib ${params.finches_lib}"               : ''
    def n_cpu           = task.cpus                  ?: 1
    def only_main_arg   = params.only_main_isoforms  ? '--only_main_isoforms'                              : ''
    def max_len_arg     = params.finches_max_seq_len ? "--max_seq_len ${params.finches_max_seq_len}"        : '--max_seq_len 3000'
    // Use finches_python when set (needs finches + aiupred conda env); else fall
    // back to the script on PATH (assumes finches importable from default python3).
    def py              = params.finches_python       ?: null
    def invoke          = py ? "${py} ${baseDir}/bin/create_finches_worker.py" : "create_finches_worker.py"
    """
    ${invoke} \\
        --loc_chrom  ${loc_chrom} \\
        --output_dir . \\
        --n_cpu      ${n_cpu} \\
        ${finches_lib_arg} \\
        ${only_main_arg} \\
        ${max_len_arg}
    """

    stub:
    """
    printf 'Protein_ID\\tPosition\\tWT_AA\\tMut_AA\\tWT_Epsilon\\tMut_Epsilon\\tDelta_Epsilon\\n' > finches_saturation.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// CATGRANULE_MAP — catGRANULE 2.0 LLPS propensity (Monti et al.). Per-residue
// profile + per-protein RandomForest LLPS score. catGRANULE's deps live in a
// separate env, so the worker delegates to params.catgranule_python. Missing
// env/lib → empty track, never crashes. Worker: create_catgranule_worker.py.
// ──────────────────────────────────────────────────────────────────────────
process CATGRANULE_MAP {
    tag  { "catgranule_map" }
    label 'process_medium'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/phase_separation"
                                : "${params.outdir}/final/phase_separation" },
        mode: 'copy'
    )

    input:
    path loc_chrom

    output:
    path "catgranule.tsv", emit: catgranule

    script:
    def only_main_arg = params.only_main_isoforms ? '--only_main_isoforms' : ''
    def py_arg        = params.catgranule_python ? "--catgranule_python ${params.catgranule_python}" : ''
    def lib_arg       = params.catgranule_lib    ? "--catgranule_lib ${params.catgranule_lib}"       : ''
    """
    create_catgranule_worker.py \\
        --seq_table ${loc_chrom} \\
        --outdir    . \\
        ${py_arg} \\
        ${lib_arg} \\
        ${only_main_arg}
    """

    stub:
    """
    printf 'Protein_ID\\tPosition\\tcatgranule_score\\tcatgranule_total\\n' > catgranule.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// PLAAC_MAP — prion-like domain HMM (Lancaster 2014). Worker:
// create_plaac_worker.py. Missing jar/Java → empty track, never crashes.
// ──────────────────────────────────────────────────────────────────────────
process PLAAC_MAP {
    tag  { "plaac_map" }
    label 'process_medium'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/phase_separation"
                                : "${params.outdir}/final/phase_separation" },
        mode: 'copy'
    )

    input:
    path loc_chrom

    output:
    path "plaac.tsv", emit: plaac

    script:
    def only_main_arg = params.only_main_isoforms ? '--only_main_isoforms' : ''
    def jar_arg       = params.plaac_jar ? "--plaac_jar ${params.plaac_jar}" : ''
    """
    create_plaac_worker.py \\
        --seq_table ${loc_chrom} \\
        --outdir    . \\
        ${jar_arg} \\
        ${only_main_arg}
    """

    stub:
    """
    printf 'Protein_ID\\tPosition\\tplaac_score\\tin_PRD\\n' > plaac.tsv
    """
}
