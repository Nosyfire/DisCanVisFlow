/*
 * modules/functional.nf — Module 5 functional — GO, polymorphism, PEM, coiled-coils, PPI, ScanSite, conservation
 *
 * Processes: FETCH_GO, GO_MAP, POLYMORPHISM_MAP, PEM_MAP, PEM_TRANSFER_MAP, COILEDCOILS_MAP, PPI_PREPROCESS, PPI_MAP, CONSERVATION_MAP, GOPHER_RECOMPUTE, SCANSITE_MAP
 * (split out of the former annotation_mapping.nf monolith)
 */



// ──────────────────────────────────────────────────────────────────────────
// FETCH_GO  — Download GOA human annotation + GO OBO (storeDir cached)
// ──────────────────────────────────────────────────────────────────────────
process FETCH_GO {
    tag  { "go_annotation" }
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/go" : "${params.ref_dir}/go" }

    output:
    path "goa_human.gaf.gz", emit: goa
    path "go.obo",           emit: obo

    script:
    """
    wget -q -O goa_human.gaf.gz \\
        https://current.geneontology.org/annotations/goa_human.gaf.gz
    wget -q -O go.obo \\
        https://current.geneontology.org/ontology/go.obo
    """

    stub:
    """
    touch goa_human.gaf.gz go.obo
    """
}


// ──────────────────────────────────────────────────────────────────────────
// GO_MAP  — Map GO terms to GENCODE Protein_IDs
// ──────────────────────────────────────────────────────────────────────────
process GO_MAP {
    tag  { "go_map" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path goa
    path go_obo

    output:
    path "go_terms.tsv", emit: go_terms

    script:
    """
    create_go_worker.py \\
        --loc_chrom  ${loc_chrom} \\
        --goa        ${goa} \\
        --go_obo     ${go_obo} \\
        --output_dir .
    """

    stub:
    """
    echo -e "Protein_ID\tEntry_Isoform\tGO_Term\tname\tnamespace\tdef\talt_id\tis_a" > go_terms.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// POLYMORPHISM_MAP  — Natural variant / SNP annotation (UniProt REST API)
// ──────────────────────────────────────────────────────────────────────────
process POLYMORPHISM_MAP {
    tag  { "polymorphism_map" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path combined_map                            // combined_map.map or NO_FILE
    path snp_common,  stageAs: 'snp_common/*'   // legacy common_poly.out or NO_FILE
    path snp_all,     stageAs: 'snp_all/*'      // legacy all_poly.out or NO_FILE
    path snp_pos_tsv, stageAs: 'snp_pos/*'      // polymorphism_pos.tsv or NO_FILE
    path dbsnp_bb                                // dbSnp*Common.bb bigBed or NO_FILE
    path dbsnp_maf_gz,  stageAs: 'dbsnp_maf/*'  // compact dbSNP MAF TSV (gzipped) or NO_FILE
    path gnomad_maf_gz, stageAs: 'gnomad_maf/*' // compact gnomAD MAF TSV (gzipped) or NO_FILE
    path setup_done, stageAs: 'setup.done'       // sentinel: ensures bigBedToBed is installed

    output:
    path "polymorphism.tsv", emit: polymorphism

    script:
    def ucsc_arg      = (params.ucsc_bin) ? "--ucsc_bin ${params.ucsc_bin}" : ""
    def maf_arg       = (dbsnp_maf_gz.name  != 'NO_FILE') ? "--dbsnp_maf ${dbsnp_maf_gz}"  : ""
    def gnomad_arg    = (gnomad_maf_gz.name != 'NO_FILE') ? "--gnomad_maf ${gnomad_maf_gz}" : ""
    def dbsnp_api_arg = params.fetch_dbsnp_api  ? "--use_dbsnp_api"  : ""
    def gnomad_api_arg= params.fetch_gnomad_api ? "--use_gnomad_api" : ""
    """
    create_polymorphism_worker.py \\
        --loc_chrom    ${loc_chrom} \\
        --combined_map ${combined_map} \\
        --snp_common   ${snp_common} \\
        --snp_all      ${snp_all} \\
        --snp_pos_tsv  ${snp_pos_tsv} \\
        --dbsnp_bb     ${dbsnp_bb} \\
        ${maf_arg} \\
        ${gnomad_arg} \\
        ${dbsnp_api_arg} \\
        ${gnomad_api_arg} \\
        ${ucsc_arg} \\
        --output_dir   .
    """

    stub:
    """
    echo -e "Protein_ID\tPosition\trsid\tref\talt\tallele_frequency\tType" > polymorphism.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// PEM_MAP  — Predicted ELM Motifs (HotspotPEM)
// ──────────────────────────────────────────────────────────────────────────
process PEM_MAP {
    tag  { "pem_map" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path pem_dataset

    output:
    path "pem_core_motifs.tsv", emit: pem

    script:
    """
    create_pem_worker.py \\
        --loc_chrom   ${loc_chrom} \\
        --pem_dataset ${pem_dataset} \\
        --output_dir  .
    """

    stub:
    """
    echo -e "Protein_ID\tELM_Accession\tELMIdentifier\tELMType\tStart\tEnd\tInstanceLogic\tReferences\tMethods\tPDB\tOrganism\tFound_Known" > pem_core_motifs.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// PEM_TRANSFER_MAP  — Map PEM motifs to alternative isoforms
// ──────────────────────────────────────────────────────────────────────────
process PEM_TRANSFER_MAP {
    tag  { "pem_transfer" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path pem_tsv

    output:
    path "pem_core_motifs_mapped.tsv", emit: pem_mapped

    script:
    """
    create_pem_transfer_worker.py \\
        --loc_chrom ${loc_chrom} \\
        --pem_tsv   ${pem_tsv} \\
        --outdir    .
    """

    stub:
    """
    echo -e "Protein_ID\tELM_Accession\tStart\tEnd\thomology_transfer" > pem_core_motifs_mapped.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// COILEDCOILS_MAP  — DeepCoil coiled-coil region prediction
// ──────────────────────────────────────────────────────────────────────────
process COILEDCOILS_MAP {
    tag  { "coiledcoils_map" }
    label 'process_high'

    // Per-chunk when scattering; a MERGE step publishes the combined table.
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy',
        enabled: ( ((params.scatter_chunks ?: 1) as Integer) <= 1 )
    )

    input:
    path loc_chrom
    path setup_done,       stageAs: 'setup.done'         // sentinel from SETUP_DEPS (ordering)
    path deepcoil_py_file, stageAs: 'deepcoil_py.txt'   // detected path from SETUP_DEPS

    output:
    path "coiled_coils.tsv", emit: coiled_coils

    script:
    def deepcoil_param = params.deepcoil_python ?: ""
    """
    _deepcoil_py="${deepcoil_param}"
    [[ -z "\${_deepcoil_py}" ]] && [[ -f deepcoil_py.txt ]] && \\
        _deepcoil_py="\$(cat deepcoil_py.txt | tr -d '[:space:]')"

    create_coiledcoils_worker.py \\
        --loc_chrom       ${loc_chrom} \\
        --deepcoil_python "\${_deepcoil_py:-python}" \\
        --n_cpu           ${task.cpus} \\
        --output_dir      .
    """

    stub:
    """
    echo -e "Protein_ID\tProb_scores" > coiled_coils.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// PPI_PREPROCESS  — Module 5j-prep: convert raw MiTab → standard Interaction_*.tsv
//
// Filters raw IntAct/BioGRID/HIPPIE downloads to human-only interactions and
// outputs the standard tab-separated format consumed by PPI_MAP.
// Outputs are cached in references/ppi/ via storeDir so the 300 MB downloads
// are parsed only once.
// ──────────────────────────────────────────────────────────────────────────
process PPI_PREPROCESS {
    tag  { "ppi_preprocess" }
    label 'process_medium'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/ppi/processed" : "${params.ref_dir}/ppi/processed" }

    input:
    path intact_raw
    path biogrid_raw
    path hippie_raw

    output:
    path "Interaction_intact.tsv",  emit: intact
    path "Interaction_biogrid.tsv", emit: biogrid
    path "Interaction_hippie.tsv",  emit: hippie

    script:
    """
    create_ppi_preprocess_worker.py \\
        --intact   ${intact_raw} \\
        --biogrid  ${biogrid_raw} \\
        --hippie   ${hippie_raw} \\
        --outdir   .
    """

    stub:
    """
    echo -e "Accession A\tAccession B\tID Interactor A\tID Interactor B\tInteraction Detection Methods\tPublication Identifiers\tConfidence Value" > Interaction_intact.tsv
    echo -e "Accession A\tAccession B\tID Interactor A\tID Interactor B\tInteraction Detection Methods\tPublication Identifiers\tConfidence Value" > Interaction_biogrid.tsv
    echo -e "Accession A\tAccession B\tID Interactor A\tID Interactor B\tInteraction Detection Methods\tPublication Identifiers\tConfidence Value" > Interaction_hippie.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// PPI_MAP  — Module 5j: Protein-Protein Interactions (BioGRID + IntAct + HIPPIE)
// ──────────────────────────────────────────────────────────────────────────
process PPI_MAP {
    tag  { "ppi_map" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path intact_file
    path biogrid_file
    path hippie_file

    output:
    path "interactions.tsv", emit: interactions

    script:
    """
    create_ppi_worker.py \\
        --seq_table ${loc_chrom} \\
        --intact    ${intact_file} \\
        --biogrid   ${biogrid_file} \\
        --hippie    ${hippie_file} \\
        --outdir    .
    """

    stub:
    """
    echo -e "Protein_ID_A\tProtein_ID_B\tdatabase\tnumber_of_pubmed" > interactions.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// CONSERVATION_MAP  — Module 7: GOPHER + phastCons per-residue conservation
// ──────────────────────────────────────────────────────────────────────────
process CONSERVATION_MAP {
    tag  { "conservation_map" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/conservation"
                                : "${params.outdir}/final/conservation" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path combined_map
    path conservation_table   // GOPHER conservation_table.tsv (or NO_FILE)

    output:
    path "conservation_multiple_level.tsv", emit: gopher
    path "conservation_phastcons.tsv",      emit: phastcons

    script:
    def skip_gopher  = (conservation_table.name == 'NO_FILE') ? '--skip_gopher'   : ''
    def gopher_arg   = (conservation_table.name != 'NO_FILE') ? "--conservation_table ${conservation_table}" : "--conservation_table /dev/null"
    def skip_pcons   = params.skip_phastcons                  ? '--skip_phastcons' : ''
    def phastcons_arg = params.phastcons_dir                  ? "--phastcons_dir ${params.phastcons_dir}" : '--skip_phastcons'
    def bw2bg_arg    = params.bigwigtobedgraph                ? "--bigwigtobedgraph ${params.bigwigtobedgraph}" : ''
    """
    create_conservation_worker.py \\
        --seq_table          ${loc_chrom} \\
        ${gopher_arg} \\
        --combined_map       ${combined_map} \\
        --outdir             . \\
        ${skip_gopher} \\
        ${phastcons_arg} \\
        ${bw2bg_arg} \\
        ${skip_pcons}
    """

    stub:
    """
    echo -e "Protein_ID\tEntry_Isoform\tlevel\tconservationscores" > conservation_multiple_level.tsv
    echo -e "Protein_ID\tEntry_Isoform\tconservationscores"        > conservation_phastcons.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// GOPHER_RECOMPUTE — recompute the GOPHER multi-level conservation table from
// orthologue alignments (instead of consuming a precomputed table). Slow but the
// goal calls for recomputation. Supply orthologue alignments via --gopher_aln_dir
// (a dir of <ACC>.orthaln.fas), or a command template (--gopher_cmd) that runs the
// real GOPHER/SLiMSuite to GENERATE them. Output feeds CONSERVATION_MAP.
// ──────────────────────────────────────────────────────────────────────────
process GOPHER_RECOMPUTE {
    tag  { "gopher_recompute" }
    label 'process_high'

    input:
    path loc_chrom
    path(aln_dir,   stageAs: 'gopher_aln/*')    // orthologue alignments or NO_FILE
    path(taxon_map, stageAs: 'gopher_taxon/*')  // species→levels TSV or NO_FILE

    output:
    path "conservation_table.tsv", emit: table

    script:
    def aln_arg   = (aln_dir.name   != 'NO_FILE') ? "--aln_dir ${aln_dir}"     : ''
    def taxon_arg = (taxon_map.name != 'NO_FILE') ? "--taxon_map ${taxon_map}" : ''
    def cmd_arg   = params.gopher_cmd ? "--gopher_cmd '${params.gopher_cmd}'"   : ''
    """
    run_gopher_worker.py \\
        --seq_table ${loc_chrom} \\
        --out conservation_table.tsv \\
        ${aln_arg} \\
        ${taxon_arg} \\
        ${cmd_arg}
    """

    stub:
    """ printf 'uniprot_acc\\tlevel\\tconservation_score\\n' > conservation_table.tsv """
}

// ──────────────────────────────────────────────────────────────────────────
// SCANSITE_MAP  — Module 5k: ScanSite 4.0 kinase / phospho motif predictions
// ──────────────────────────────────────────────────────────────────────────
process SCANSITE_MAP {
    tag  { "scansite" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path scansite_tsv   // pre-computed scansite.tsv or NO_FILE

    output:
    path "scansite.tsv", emit: scansite

    script:
    def use_api = params.scansite_use_api ? '--use_api' : ''
    def str_arg = params.scansite_stringency ? "--stringency ${params.scansite_stringency}" : ''
    """
    create_scansite_worker.py \\
        --seq_table    ${loc_chrom} \\
        --scansite_tsv ${scansite_tsv} \\
        --outdir       . \\
        ${use_api} ${str_arg}
    """

    stub:
    """
    echo -e "Protein_ID\tmotifName\tmotifShortName\tscore\tsite\tsiteSequence\tStart\tEnd" > scansite.tsv
    """
}
