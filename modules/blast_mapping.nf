/*
 * modules/blast_mapping.nf
 *
 * Module 1 — Transcript → UniProt ID Mapping
 *
 * Takes a pre-computed reciprocal BLAST hit table (bestsequences.tsv) and
 * the optional full-isoform BLAST table, then runs create_id_map_worker.py
 * to produce:
 *
 *   bestmaps_blast_gene_transcript.tsv  — best single UniProt per transcript
 *   blastmaps_isoforms.tsv             — formatted isoform table (optional)
 *
 * ────────────────────────────────────────────────────────────────────────────
 * Input channels
 *   blast_tsv     : path  — reciprocal BLAST hit TSV (required)
 *   isoforms_tsv  : path  — full isoform BLAST TSV, OR assets/NO_FILE
 *
 * Output channels
 *   id_map        : path  — bestmaps_blast_gene_transcript.tsv
 *   isoforms      : path  — blastmaps_isoforms.tsv  (optional, may be absent)
 * ────────────────────────────────────────────────────────────────────────────
 */

process ID_MAP {

    tag { blast_tsv.simpleName }

    label 'process_low'

    publishDir(
        path:    "${params.outdir}/blast",
        mode:    'copy',
        pattern: "*.tsv",
    )

    input:
    path blast_tsv
    path isoforms_tsv   // pass assets/NO_FILE when no isoform table available

    output:
    path "bestmaps_blast_gene_transcript.tsv", emit: id_map
    path "blastmaps_isoforms.tsv",             emit: isoforms, optional: true

    script:
    /*
     * Only pass --isoforms_tsv when the caller supplied a real file.
     * The sentinel 'NO_FILE' has zero size; we skip it explicitly.
     */
    def isoforms_arg = (isoforms_tsv.name != 'NO_FILE' && isoforms_tsv.size() > 0)
                       ? "--isoforms_tsv ${isoforms_tsv}"
                       : ''
    """
    create_id_map_worker.py \\
        --blast_tsv     ${blast_tsv} \\
        --output_dir    . \\
        --database      ${params.blast_database} \\
        --coverage      ${params.blast_coverage} \\
        --mapping_mode  ${params.mapping_mode ?: 'main_isoform_mapping'} \\
        ${isoforms_arg}
    """

    stub:
    /*
     * Stub block — used with `nextflow run -stub` for fast DAG validation
     * without actually running Python.
     */
    """
    echo -e "Entry_Name\tGene_Uniprot\tGene_Gencode\tName\tTranscript name\ttranscript_stable_id\tTranscript ID\tEntry_Isoform\tDatabase\tcoverage_x\tcoverage_y\tcoverage\talignmentpuntcuality" \
        > bestmaps_blast_gene_transcript.tsv
    echo -e "STUB_HUMAN\tSTUB\tSTUB\tStub protein\tSTUB-201\tENST00000000001\tENST00000000001.1\tP00001\tUniprot/SWISSPROT\t100.0\t100.0\t100.0\tidentical" \
        >> bestmaps_blast_gene_transcript.tsv
    """
}
