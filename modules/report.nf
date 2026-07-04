/*
 * modules/report.nf — Mapping report — per-run annotation coverage audit (runs last)
 *
 * Processes: MAPPING_REPORT
 * (split out of the former annotation_mapping.nf monolith)
 */


// ──────────────────────────────────────────────────────────────────────────
// MAPPING_REPORT — per-gene before/after annotation mapping audit (Markdown)
// Explains, for one gene: annotation sources (local path / downloaded) and raw
// counts BEFORE mapping; the selected isoforms with their alignment + genomic
// locations; and per-isoform which annotations mapped (and which did NOT).
// Python worker: create_mapping_report_worker.py
// ──────────────────────────────────────────────────────────────────────────
process MAPPING_REPORT {
    tag  { "mapping_report" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/mapping_reports"
                                : "${params.outdir}/mapping_reports" },
        mode: 'copy'
    )

    input:
    path seq_table
    path gate_files          // collected final outputs — ordering only (run last)

    output:
    path "*.md",              emit: reports
    path "*.tsv", optional: true, emit: coverage_tsv

    script:
    // Absolute published locations the worker reads directly.
    def base = params.gene_dir ? "${params.outdir}/${params.gene_dir}"
                               : "${params.outdir}"
    def final_dir = "${workflow.launchDir}/${base}/final"
    def inter_dir = "${workflow.launchDir}/${base}/intermediate"
    // Source provenance: the worker has a built-in registry (relative repo paths
    // for legacy_data + descriptions for downloaded/computed/derived). Here we
    // override only the param-driven EXTERNAL file paths; the worker relativizes
    // any absolute path under --launch_dir.
    def ov = []
    if ( params.mavedb_tsv )     ov << "--source 'pathogenicity/mavedb.tsv=local|${params.mavedb_tsv}'"
    if ( params.proteingym_tsv ) ov << "--source 'pathogenicity/proteingym.tsv=local|${params.proteingym_tsv}'"
    if ( params.dbnsfp_tsv )     ov << "--source 'pathogenicity/dbnsfp_scores.tsv=local|${params.dbnsfp_tsv}'"
    // Extract dbNSFP version string from raw dir path (e.g. /path/to/dbNSFP4.8a → 4.8a)
    def dbnsfp_ver = params.dbnsfp_raw_dir
        ? (params.dbnsfp_raw_dir =~ /dbNSFP(\S+)/)[0]?.getAt(1) ?: "unknown"
        : (params.dbnsfp_tsv ? "pre-mapped" : "n/a")
    if ( params.omim_tsv ) {     ov << "--source 'disease/omim_disease.tsv=local|${params.omim_tsv}'"
                                 ov << "--source 'disease/omim_mutations.tsv=local|${params.omim_tsv}'" }
    if ( params.gopher_conservation_table ) ov << "--source 'conservation/conservation_multiple_level.tsv=local|${params.gopher_conservation_table}'"
    if ( params.phastcons_dir )  ov << "--source 'conservation/conservation_phastcons.tsv=local|${params.phastcons_dir}'"
    if ( params.pem_dataset )    ov << "--source 'annotations/pem_core_motifs.tsv=local|${params.pem_dataset}'"
    if ( params.depmap_tsv )     ov << "--source 'mutations/DepMap/depmap_mutations.tsv=local|${params.depmap_tsv}'"
    if ( params.skip_pfam_api )  ov << "--source 'annotations/pfam_domains.tsv=local|(Pfam API skipped)'"
    def src_args = ov.join(' ')
    // Reference FASTAs → data-source versions + input-scale base counts in the report.
    // reldate.txt is captured by FETCH_UNIPROT_FASTA into references/uniprot/;
    // the worker no-ops gracefully when the file is absent (e.g. --data local).
    def reldate_path = params.ref_dir ? "${params.ref_dir}/uniprot/uniprot_reldate.txt" : ""
    def ref_args = [
        params.gencode_fasta         ? "--gencode_fasta '${params.gencode_fasta}'"                 : "",
        params.uniprot_fasta         ? "--uniprot_fasta '${params.uniprot_fasta}'"                 : "",
        params.uniprot_isoform_fasta ? "--uniprot_isoform_fasta '${params.uniprot_isoform_fasta}'" : "",
        reldate_path                 ? "--uniprot_reldate '${reldate_path}'"                        : "",
    ].findAll { it }.join(' ')
    """
    # ── capture tool / data versions for reproducibility ──
    {
      echo "Nextflow: ${workflow.nextflow.version}"
      echo "Pipeline: ${workflow.manifest.version}"
      echo "Conda env: \$( (basename \"\${CONDA_PREFIX:-n/a}\") 2>/dev/null )"
      echo "Python: \$(python3 --version 2>&1 | awk '{print \$2}')"
      echo "pandas: \$(python3 -c 'import pandas;print(pandas.__version__)' 2>/dev/null)"
      echo "blastp: \$(blastp -version 2>/dev/null | head -1 | awk '{print \$2}')"
      echo "blat: \$(blat 2>&1 | head -1 | grep -oE '[0-9]+x[0-9]+' | head -1)"
      echo "Mapping mode: ${params.mapping_mode}"
    } > versions.txt 2>/dev/null || true

    # report_fix_2026-06-23: large runs emit mapping_coverage.tsv+mapping_summary.md
    # instead of one MD per gene (19k files was too slow for full-proteome runs)
    create_mapping_report_worker.py \\
        --seq_table        ${seq_table} \\
        --final_dir        '${final_dir}' \\
        --intermediate_dir '${inter_dir}' \\
        --mapping_mode     ${params.mapping_mode} \\
        --command          '${workflow.commandLine}' \\
        --pipeline_version '${workflow.manifest.version}' \\
        --nextflow_version '${workflow.nextflow.version}' \\
        --profile          'project=${params.project};data=${params.data};machine=${params.machine};env=${params.env}' \\
        --run_name         '${workflow.runName}' \\
        --start_time       '${workflow.start}' \\
        --work_dir         '${workflow.workDir}' \\
        --launch_dir       '${workflow.launchDir}' \\
        --versions_file    versions.txt \\
        --per_gene_md_threshold ${params.per_gene_md_threshold ?: 50} \\
        --dbnsfp_version   '${dbnsfp_ver}' \\
        ${ref_args} \\
        ${src_args} \\
        --outdir .
    """

    stub:
    """
    printf '# Mapping summary (stub)\\n' > mapping_summary.md
    touch mapping_coverage.tsv
    """
}
