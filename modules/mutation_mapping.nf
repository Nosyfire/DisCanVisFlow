/*
 * modules/mutation_mapping.nf  — Module 4: Mutation Mapping
 *
 * Processes
 * ─────────
 *  FETCH_CLINVAR      Download ClinVar VCF (GRCh38) from NCBI FTP (storeDir cached).
 *
 *  FETCH_CLINVAR_SUBMISSIONS
 *                     Download ClinVar submission_summary.txt.gz (storeDir cached).
 *
 *  CLINVAR_DATES      Aggregate per-SCV submission dates → one row per variant,
 *                     keyed on CLNHGVS (the mapped TSVs' `Mutation` column).
 *
 *  MUTATION_MAP       Map genomic mutations → protein positions on all isoforms.
 *                     Accepts ClinVar VCF, MAF, or generic VCF.
 *                     Outputs: Missense / Frameshift / Nonsense / Indel TSVs.
 *
 * Channel contracts
 * ─────────────────
 *  FETCH_CLINVAR.out.vcf          : path  (clinvar.vcf.gz)
 *  MUTATION_MAP.out.missense      : path
 *  MUTATION_MAP.out.frameshift    : path
 *  MUTATION_MAP.out.nonsense      : path
 *  MUTATION_MAP.out.indel         : path
 *  MUTATION_MAP.out.stats         : path
 *  FETCH_CLINVAR_SUBMISSIONS.out.submissions : path
 *  CLINVAR_DATES.out.dates        : path
 */

// ──────────────────────────────────────────────────────────────────────────
// FETCH_CLINVAR
// ──────────────────────────────────────────────────────────────────────────
process FETCH_CLINVAR {

    tag  { "clinvar_grch38" }
    label 'process_low'

    storeDir "${params.ref_dir}/clinvar"

    output:
    path "clinvar_grch38.vcf.gz", emit: vcf

    script:
    """
    wget -q -O clinvar_grch38.vcf.gz \\
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    """

    stub:
    """
    touch clinvar_grch38.vcf.gz
    """
}


// ──────────────────────────────────────────────────────────────────────────
// FETCH_CLINVAR_SUBMISSIONS
//   submission_summary.txt.gz — one row per submitted record (SCV), carrying the
//   per-submitter DateLastEvaluated the VCF does not have. Always the current
//   release; joined by VariationID, so skew against a pinned VCF is harmless.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_CLINVAR_SUBMISSIONS {

    tag  { "clinvar_submissions" }
    label 'process_low'

    storeDir "${params.ref_dir}/clinvar"

    output:
    path "submission_summary.txt.gz", emit: submissions

    script:
    """
    wget -q -O submission_summary.txt.gz \\
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
    """

    stub:
    """
    touch submission_summary.txt.gz
    """
}


// ──────────────────────────────────────────────────────────────────────────
// CLINVAR_DATES
//   One compact row per ClinVar variant with its submission-date summary.
//   Keyed on CLNHGVS == the `Mutation` column of the mapped mutation TSVs, so
//   it stays variant-keyed instead of being duplicated across every isoform.
// ──────────────────────────────────────────────────────────────────────────
process CLINVAR_DATES {

    tag  { "clinvar_dates" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir
            ? "${params.outdir}/${params.gene_dir}/final/mutations/ClinVar"
            : "${params.outdir}/final/mutations/ClinVar" },
        mode: 'copy'
    )

    input:
    path clinvar_vcf
    path submissions

    output:
    path "clinvar_submission_dates.tsv", emit: dates

    script:
    """
    create_clinvar_dates_worker.py \\
        --clinvar_vcf        ${clinvar_vcf} \\
        --submission_summary ${submissions} \\
        --output             clinvar_submission_dates.tsv
    """

    stub:
    """
    echo -e "VariationID\tMutation\tChromosome\tPosition\tRef\tAlt\tGene\tlast_evaluated\tfirst_evaluated\tn_submissions\tlast_submitter\tlast_clinical_significance\tlast_review_status" > clinvar_submission_dates.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// MUTATION_MAP
//   Maps mutations to protein positions on ALL isoforms of each gene.
//   Input format is controlled by the `input_format` val:
//     'clinvar_vcf' → --clinvar_vcf
//     'maf'         → --maf
//     'vcf'         → --vcf (generic VCF)
// ──────────────────────────────────────────────────────────────────────────
process MUTATION_MAP {

    tag  { "mutation_map_${source}" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir
            ? "${params.outdir}/${params.gene_dir}/final/mutations/${source}"
            : "${params.outdir}/final/mutations/${source}" },
        mode: 'copy'
    )

    input:
    path combined_map          // Module 3 output
    path loc_chrom             // Module 2 output
    path mutation_file         // ClinVar VCF, MAF, or generic VCF
    val  source                // 'ClinVar' | 'TCGA' | 'CBioportal' | …
    val  input_format          // 'clinvar_vcf' | 'maf' | 'vcf'

    output:
    path "Missense_filter_mutations_mapped.tsv",   emit: missense
    path "Frameshift_filter_mutations_mapped.tsv", emit: frameshift
    path "Nonsense_filter_mutations_mapped.tsv",   emit: nonsense
    path "Indel_filter_mutations_mapped.tsv",      emit: indel
    path "mutation_stats.tsv",                     emit: stats

    script:
    def mut_flag = input_format == 'clinvar_vcf' ? "--clinvar_vcf ${mutation_file}"
                 : input_format == 'maf'         ? "--maf ${mutation_file}"
                 :                                  "--vcf ${mutation_file}"
    def no_expand = params.no_isoform_expand ? "--no_isoform_expand" : ""
    def hgvsp_val = params.mutation_validate_hgvsp ? "" : "--no_hgvsp_validation"
    def hypermut  = params.mutation_hypermutation_threshold > 0
        ? "--hypermutation_threshold ${params.mutation_hypermutation_threshold}" : ""
    """
    create_mutation_map_worker.py \\
        --combined_map  ${combined_map} \\
        --loc_chrom     ${loc_chrom} \\
        ${mut_flag} \\
        --source        "${source}" \\
        --output_dir    . \\
        ${no_expand} \\
        ${hgvsp_val} \\
        ${hypermut}
    """

    stub:
    def header = "Protein_ID\tAccession\tGene\tMutation Description\tMutation\tProtein_position\tStudy Abbrevation\tStudy Name\tSample name\tStart_Position\thomology_transfer\tClinicalSignificance\tPhenotypeList\tPhenotypeIDS\tReviewStatus\tRCVaccession\tMONDO_ID\tMeSH_ID"
    """
    echo -e "${header}" > Missense_filter_mutations_mapped.tsv
    cp Missense_filter_mutations_mapped.tsv Frameshift_filter_mutations_mapped.tsv
    cp Missense_filter_mutations_mapped.tsv Nonsense_filter_mutations_mapped.tsv
    cp Missense_filter_mutations_mapped.tsv Indel_filter_mutations_mapped.tsv
    echo -e "source\ttotal_resolved\tmissense\tframeshift\tnonsense\tindel" > mutation_stats.tsv
    echo -e "${source}\t0\t0\t0\t0\t0" >> mutation_stats.tsv
    """
}
