/*
 * modules/sequence_process.nf — Module 2: Sequence Process
 *
 * Merges BLAST ID-map results with GENCODE GTF annotations and protein
 * sequences to produce the canonical transcript/isoform table.
 *
 * Input  channels
 * ───────────────
 *   blast_best    path   bestmaps_blast_gene_transcript.tsv  (Module 1)
 *   isoforms_tsv  path   blastmaps_isoforms.tsv              (Module 1)
 *   gtf           path   GENCODE annotation GTF (.gz)
 *   translations  path   GENCODE pc_translations FASTA
 *
 * Output channels
 * ───────────────
 *   loc_chrom       loc_chrom_with_names.tsv
 *   loc_chrom_seq   loc_chrom_with_names_isoforms_with_seq.tsv  ← main output
 *   main_isoform    loc_chrom_with_names_main_isoform.tsv
 */

process SEQUENCE_PROCESS {
    tag "cutoff_${params.seq_cutoff}"
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/sequence"
                                : "${params.outdir}/final/sequence" },
        mode: 'copy', pattern: '*.tsv'
    )

    input:
    path blast_best
    path isoforms_tsv
    path gtf
    path translations

    output:
    path 'loc_chrom_with_names.tsv',                   emit: loc_chrom
    path 'loc_chrom_with_names_isoforms_with_seq.tsv', emit: loc_chrom_seq
    path 'loc_chrom_with_names_main_isoform.tsv',      emit: main_isoform
    path 'loc_chrom_with_names_isoforms_only.tsv',     optional: true, emit: isoforms_only

    script:
    def iso_arg = (isoforms_tsv.name != 'NO_FILE') \
                  ? "--blast_isoforms ${isoforms_tsv}" \
                  : ''
    """
    create_sequence_table_worker.py \\
        --blast_best   ${blast_best} \\
        --gtf          ${gtf} \\
        --translations ${translations} \\
        --output_dir   . \\
        --cutoff       ${params.seq_cutoff} \\
        --mapping_mode ${params.mapping_mode ?: 'main_isoform_mapping'} \\
        ${iso_arg}
    """

    stub:
    """
    printf 'Transcript ID\tGene_Gencode\tGene_Uniprot\tAccession\tChromosome\tgene_type\tEnsembl_canonical\tMANE_Select\tappris_principal_1\tmain_isoform\tmain_isoform_id\tSequence\n' \\
        > loc_chrom_with_names_isoforms_with_seq.tsv
    printf 'ENST00000423430.6\tENSG00000132155.11\tRAF1\tP04049\tchr3\tprotein_coding\tyes\tyes\tyes\tyes\tRAF1-201\tMANTIQQFLK\n' \\
        >> loc_chrom_with_names_isoforms_with_seq.tsv

    cp loc_chrom_with_names_isoforms_with_seq.tsv loc_chrom_with_names.tsv
    cp loc_chrom_with_names_isoforms_with_seq.tsv loc_chrom_with_names_main_isoform.tsv
    """
}
