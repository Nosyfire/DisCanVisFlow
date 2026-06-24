/*
 * modules/blast_search.nf
 *
 * Module 0 — Reciprocal BLAST-based transcript → UniProt mapping
 *
 * Processes (in execution order)
 * ─────────────────────────────────────────────────────────────────────────
 *  SUBSET_FASTA    — extract sequences whose header contains a target string
 *                    (pass-through when params.target_gene is empty)
 *
 *  MAKEBLASTDB     — format a FASTA as a BLAST protein database
 *
 *  BLASTP          — run blastp (outfmt 5 / XML) of query FASTA vs a DB
 *
 *  MERGE_BLAST_HITS — parse two reciprocal XML outputs with
 *                    create_blast_table_worker.py:
 *                    → bestsequences.tsv, allsequences.tsv,
 *                      isoformssequences.tsv
 *
 * Channel contracts
 * ─────────────────────────────────────────────────────────────────────────
 *  SUBSET_FASTA.out.subset      : path  (subsetted FASTA)
 *
 *  MAKEBLASTDB.out.blastdb      : tuple val(db_name), path(db_dir)
 *                                 db_dir contains all DB index files under
 *                                 ${db_name}_blastdb/${db_name}
 *
 *  BLASTP.out.blast_xml         : path  (BLAST XML result file)
 *
 *  MERGE_BLAST_HITS.out.bestsequences  : path  bestsequences.tsv
 *  MERGE_BLAST_HITS.out.allsequences   : path  allsequences.tsv
 *  MERGE_BLAST_HITS.out.isoforms       : path  isoformssequences.tsv
 */


// ──────────────────────────────────────────────────────────────────────────
// SUBSET_FASTA
//   Subset a FASTA to sequences whose header contains `search_term`.
//   When search_term is '' the input is copied unchanged (pass-through).
// ──────────────────────────────────────────────────────────────────────────
process SUBSET_FASTA {

    tag  { "${label}${search_term ? ':' + search_term : ':all'}" }
    label 'process_low'

    input:
    path  input_fasta
    val   search_term   // empty string = copy all sequences
    val   label         // 'uniprot' or 'gencode' (used for output file name)

    output:
    path "${label}_subset.fasta", emit: subset

    script:
    def search_arg = ''
    if (search_term) {
        if (search_term.contains(',')) {
            // Multi-gene mode: comma-separated list → build patterns per fasta_type
            def ftype = (label == 'uniprot') ? 'uniprot' : (label == 'cdna' ? 'cdna' : 'gencode')
            search_arg = "--gene_list '${search_term}' --fasta_type ${ftype}"
        } else {
            // Single-gene mode: build exact pattern inline (existing behaviour)
            // For GENCODE / cDNA (pipe-delimited headers):  match |GENE| exactly,
            // so that 'RAF1' does not hit 'TRAF1', 'ZTRAF1', etc.
            // For UniProt (space-delimited):  GN=GENE tag is already specific;
            // we append a trailing space to avoid prefix collisions (e.g. GN=RAF1A).
            def exact = (label == 'uniprot')
                ? "GN=${search_term} "
                : "|${search_term}|"
            search_arg = "--search '${exact}'"
        }
    }
    """
    subset_fasta.py \\
        --input  ${input_fasta} \\
        --output ${label}_subset.fasta \\
        ${search_arg}
    """

    stub:
    """
    cp ${input_fasta} ${label}_subset.fasta
    """
}


// ──────────────────────────────────────────────────────────────────────────
// MERGE_UNIPROT_ISOFORMS
//   For 'all_isoform_mapping' mode: append the curated Swiss-Prot isoform
//   sequences (sp| entries only, TrEMBL fragments excluded) from the UniProt
//   "additional" FASTA onto the canonical proteome FASTA, so the BLAST DB
//   contains alternative isoforms (e.g. P04049-2) and each GENCODE transcript
//   can be paired to its true isoform. A trailing space is appended to every
//   isoform header so the downstream `GN=<gene> ` exact-subset match works
//   (additional-FASTA headers end at GN=<gene> with no trailing space).
// ──────────────────────────────────────────────────────────────────────────
process MERGE_UNIPROT_ISOFORMS {
    tag  { 'merge_uniprot_isoforms' }
    label 'process_low'

    input:
    path canonical_fasta
    path isoform_fasta

    output:
    path 'uniprot_with_isoforms.fasta', emit: fasta

    script:
    """
    cat ${canonical_fasta} > uniprot_with_isoforms.fasta
    awk '
        /^>/ { keep = (\$0 ~ /^>sp\\|/); if (keep) print \$0 " "; next }
        keep { print }
    ' ${isoform_fasta} >> uniprot_with_isoforms.fasta
    """

    stub:
    """
    cp ${canonical_fasta} uniprot_with_isoforms.fasta
    """
}


// ──────────────────────────────────────────────────────────────────────────
// MAKEBLASTDB
//   Build a BLAST protein database from a FASTA file.
//   All DB index files are placed in a named subdirectory to avoid
//   filename collisions with the staged input FASTA.
// ──────────────────────────────────────────────────────────────────────────
process MAKEBLASTDB {

    tag  { db_name }
    label 'process_low'

    input:
    path fasta
    val  db_name   // 'uniprot' or 'gencode'

    output:
    tuple val(db_name), path("${db_name}_blastdb"), emit: blastdb

    script:
    // blast_bin may be '' (conda/docker PATH) or a dir path (legacy local install)
    def blast_prefix = params.blast_bin ? "${params.blast_bin}/" : ""
    """
    mkdir -p ${db_name}_blastdb
    ${blast_prefix}makeblastdb \\
        -in      ${fasta} \\
        -dbtype  prot \\
        -out     ${db_name}_blastdb/${db_name}
    """

    stub:
    """
    mkdir -p ${db_name}_blastdb
    touch ${db_name}_blastdb/${db_name}.phr \\
          ${db_name}_blastdb/${db_name}.pin \\
          ${db_name}_blastdb/${db_name}.psq
    """
}


// ──────────────────────────────────────────────────────────────────────────
// BLASTP
//   Run blastp of query_fasta against the BLAST database in db_dir.
//
//   Naming convention (mirrors legacy):
//     run_id = 'uniprotdb_gencode_query'  → query=GENCODE,  db=UniProt
//     run_id = 'gencodedb_uniprot_query'  → query=UniProt,  db=GENCODE
//
//   Output format: XML (-outfmt 5) — required by create_blast_table_worker.py
// ──────────────────────────────────────────────────────────────────────────
process BLASTP {

    tag  { run_id }
    label 'process_medium'

    input:
    path  query_fasta
    tuple val(db_name), path(db_dir)
    val   run_id          // descriptive label for output file

    output:
    path "${run_id}.xml", emit: blast_xml

    script:
    def blast_prefix = params.blast_bin ? "${params.blast_bin}/" : ""
    """
    ${blast_prefix}blastp \\
        -query        ${query_fasta} \\
        -db           ${db_dir}/${db_name} \\
        -out          ${run_id}.xml \\
        -outfmt       5 \\
        -best_hit_score_edge 0.05 \\
        -num_threads  ${task.cpus}
    """

    stub:
    // Minimal valid BLAST XML so downstream parsing doesn't crash
    """
    cat > ${run_id}.xml << 'XMLEOF'
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN"
  "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP stub</BlastOutput_version>
  <BlastOutput_reference>stub</BlastOutput_reference>
  <BlastOutput_db>stub_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>stub</BlastOutput_query-def>
  <BlastOutput_query-len>10</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>STUB_GENCODE|STUB_ENST|STUB_ENSG|STUB|STUB|STUB-201|STUB|10</Iteration_query-def>
      <Iteration_query-len>10</Iteration_query-len>
      <Iteration_hits/>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1</Statistics_db-num>
          <Statistics_db-len>10</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0</Statistics_kappa>
          <Statistics_lambda>0</Statistics_lambda>
          <Statistics_entropy>0</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
XMLEOF
    """
}


// ──────────────────────────────────────────────────────────────────────────
// MERGE_BLAST_HITS
//   Calls create_blast_table_worker.py to:
//     1. Parse both reciprocal BLAST XML outputs
//     2. Compute coverage + alignment quality per HSP
//     3. Inner-join on (Gencode, Uniprot) → reciprocal hits
//     4. Output bestsequences.tsv, allsequences.tsv, isoformssequences.tsv
//
//   Input naming contract:
//     blast1_xml = uniprotdb_gencode_query.xml  (query=GENCODE, db=UniProt)
//     blast2_xml = gencodedb_uniprot_query.xml  (query=UniProt, db=GENCODE)
// ──────────────────────────────────────────────────────────────────────────
process MERGE_BLAST_HITS {

    tag  { "coverage_${params.blast_coverage}" }
    label 'process_medium'

    publishDir(
        path:    "${params.outdir}/blast",
        mode:    'copy',
        pattern: "*.tsv",
    )

    input:
    path blast1_xml   // uniprotdb_gencode_query.xml
    path blast2_xml   // gencodedb_uniprot_query.xml

    output:
    path "bestsequences.tsv",    emit: bestsequences
    path "allsequences.tsv",     emit: allsequences
    path "isoformssequences.tsv", emit: isoforms

    script:
    """
    create_blast_table_worker.py \\
        --blast1_xml  ${blast1_xml} \\
        --blast2_xml  ${blast2_xml} \\
        --output_dir  . \\
        --coverage    ${params.blast_coverage}
    """

    stub:
    """
    printf 'Gencode\tUniprot\talignmentpuntcuality_x\tcoverage_x\talignmentpuntcuality_y\tcoverage_y\n' > bestsequences.tsv
    printf 'Gencode\tUniprot\talignmentpuntcuality_x\tcoverage_x\talignmentpuntcuality_y\tcoverage_y\n' > allsequences.tsv
    printf 'Gencode\tUniprot\talignmentpuntcuality\tcoverage\n' > isoformssequences.tsv
    """
}
