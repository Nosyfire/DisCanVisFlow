/*
 * modules/disease.nf — Module 8 disease/drivers — ClinVar disease, OMIM, cancer drivers
 *
 * Processes: CLINVAR_DISEASE_BUILD, CLINVAR_DISEASE_MAP, OMIM_MAP, CANCER_DRIVER_MAP
 * (split out of the former annotation_mapping.nf monolith)
 */


// ──────────────────────────────────────────────────────────────────────────
// CLINVAR_DISEASE_BUILD  — Module 8a: build disease table from MONDO OBO + mutations
// ──────────────────────────────────────────────────────────────────────────
process CLINVAR_DISEASE_BUILD {
    tag  { "clinvar_disease_build" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disease"
                                : "${params.outdir}/final/disease" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path mondo_obo
    path missense,   stageAs: 'mutations/Missense_filter_mutations_mapped.tsv'
    path frameshift, stageAs: 'mutations/Frameshift_filter_mutations_mapped.tsv'
    path nonsense,   stageAs: 'mutations/Nonsense_filter_mutations_mapped.tsv'
    path indel,      stageAs: 'mutations/Indel_filter_mutations_mapped.tsv'

    output:
    path "clinvar_disease.tsv", emit: clinvar_disease
    path "clinvar_disease_mutations.tsv", emit: clinvar_disease_mutations

    script:
    """
    create_clinvar_disease_build_worker.py \\
        --seq_table     ${loc_chrom} \\
        --mondo_obo     ${mondo_obo} \\
        --mutation_dir  mutations \\
        --outdir        .
    """

    stub:
    """
    echo -e "Protein_ID\tDisease\tDOID\tFinal_Category" > clinvar_disease.tsv
    echo -e "Protein_ID\tMutation\tProtein_position\tDisease\tFinal_Category" > clinvar_disease_mutations.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// CLINVAR_DISEASE_MAP  — Module 8a: ClinVar disease ontology (filter fallback)
// ──────────────────────────────────────────────────────────────────────────
process CLINVAR_DISEASE_MAP {
    tag  { "clinvar_disease" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disease"
                                : "${params.outdir}/final/disease" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path clinvar_disease_tsv     // pre-processed ClinVar disease ontology TSV or NO_FILE
    path clinvar_category_tsv    // clinvar_diseases.tsv (paper disease categories) or NO_FILE

    output:
    path "clinvar_disease.tsv", emit: clinvar_disease

    script:
    """
    create_clinvar_disease_worker.py \\
        --seq_table            ${loc_chrom} \\
        --clinvar_disease      ${clinvar_disease_tsv} \\
        --clinvar_category_tsv ${clinvar_category_tsv} \\
        --outdir               .
    """

    stub:
    """
    echo -e "Protein_ID\tDisease\tDOID\tFinal_Category" > clinvar_disease.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// OMIM_MAP  — Module 8b: OMIM disease ontology
// ──────────────────────────────────────────────────────────────────────────
process OMIM_MAP {
    tag  { "omim_disease" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disease"
                                : "${params.outdir}/final/disease" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path omim_table
    path omim_mutations

    output:
    path "omim_disease.tsv", emit: omim_disease
    path "omim_mutations.tsv", emit: omim_mutations

    script:
    // raw mode: the first input is FETCH_OMIM's raw dir (genemap2.txt); parse it
    // directly. processed mode: pre-built disease/variant tables.
    def raw_mode = (params.fetch_omim || params.omim_raw_dir) as boolean
    def omim_mut_arg = (omim_mutations.name != 'NO_FILE') ? "--omim_mutations ${omim_mutations}" : ""
    def src_arg = raw_mode ? "--mapping_mode raw --omim_raw_dir ${omim_table}"
                           : "--omim_table ${omim_table} ${omim_mut_arg}"
    """
    create_omim_worker.py \\
        --seq_table  ${loc_chrom} \\
        ${src_arg} \\
        --outdir     .
    """

    stub:
    """
    echo -e "Protein_ID\tDisease\tMIMID" > omim_disease.tsv
    echo -e "Protein_ID\tProtein_position\taa_change\tDisease" > omim_mutations.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// CANCER_DRIVER_MAP  — Module 8c: Cancer Gene Census + Compendium
// ──────────────────────────────────────────────────────────────────────────
process CANCER_DRIVER_MAP {
    tag  { "cancer_drivers" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/drivers"
                                : "${params.outdir}/final/drivers" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    // Stage the source under a distinct name: the worker writes cancer_driver.tsv
    // as output, so a same-named input symlink would be clobbered (and would
    // truncate the vendored legacy file it points to).
    path(cancer_driver, stageAs: 'cancer_driver_src.tsv')   // combined legacy file or NO_FILE
    path census_roles                                       // gene-keyed CGC roles or NO_FILE
    path compendium_roles                                   // gene-keyed Compendium roles or NO_FILE

    output:
    path "cancer_driver.tsv",     emit: combined
    path "census_driver.tsv",     emit: census
    path "compendium_driver.tsv", emit: compendium

    script:
    """
    create_cancer_driver_worker.py \\
        --seq_table          ${loc_chrom} \\
        --cancer_driver      ${cancer_driver} \\
        --census_roles       ${census_roles} \\
        --compendium_roles   ${compendium_roles} \\
        --outdir             .
    """

    stub:
    """
    echo -e "Protein_ID\tCancer Driver\tRole in Cancer\tCompendium Role" > cancer_driver.tsv
    echo -e "Protein_ID\tGene\tTier\tRole in Cancer\tTumour Types(Somatic)\tTumour Types(Germline)" > census_driver.tsv
    echo -e "Protein_ID\tGene\tROLE\tCANCER_TYPE" > compendium_driver.tsv
    """
}
