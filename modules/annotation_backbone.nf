/*
 * modules/annotation_backbone.nf — Module 5 backbone — sequence-derived annotation + isoform transfer
 *
 * Processes: FETCH_ELM, ELM_SWITCHES_MAP, ELM_CLASS_MAP, PARSE_UNIPROT_DAT, ANNOTATION_MAP, SPLIT_SEQ_TABLE, TRANSCRIPT_MAP, ISOFORM_ALIGN_MAP, HOMOLOGY_MANIFEST
 * (split out of the former annotation_mapping.nf monolith)
 */


// ──────────────────────────────────────────────────────────────────────────
// FETCH_ELM  (fallback; normally legacy_data/elm/ is used directly)
// ──────────────────────────────────────────────────────────────────────────
process FETCH_ELM {
    tag  { "elm_instances" }
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/elm" : "${params.ref_dir}/elm" }

    output:
    path "elm_instances.tsv", emit: elm_tsv

    script:
    """
    wget -q -O elm_instances.tsv \\
        'https://elm.eu.org/instances.tsv?q=organism=9606&taxon=Homo+sapiens'
    """

    stub:
    """
    printf '#ELM stub\\n#ELM stub 2\\n#ELM stub 3\\n#ELM stub 4\\n#ELM stub 5\\n' > elm_instances.tsv
    printf '"Accession"\\t"ELMType"\\t"ELMIdentifier"\\t"ProteinName"\\t"Primary_Acc"\\t"Accessions"\\t"Start"\\t"End"\\t"References"\\t"Methods"\\t"InstanceLogic"\\t"PDB"\\t"Organism"\\n' >> elm_instances.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// ELM_SWITCHES_MAP  — Module 5p: ELM molecular switches → Elm_Switches TSV
//
// Maps the raw ELM switches dataset (elm.eu.org/switches.tsv) to all GENCODE
// isoforms via substring-based coordinate remapping.  Output: Protein_ID-keyed
// elmswitches_mapped.tsv for the Elm_Switches Django model.
// ──────────────────────────────────────────────────────────────────────────
process ELM_SWITCHES_MAP {
    tag  { "elm_switches_map" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path switches_tsv   // elmswitches.tsv from legacy_data or FETCH_ELM_SWITCHES (or NO_FILE)

    output:
    path "elmswitches_mapped.tsv", emit: elm_switches

    script:
    """
    create_elm_switches_worker.py \\
        --seq_table ${loc_chrom} \\
        --switches  ${switches_tsv} \\
        --outdir    .
    """

    stub:
    """
    echo -e "Protein_ID\tEntry_Isoform\thomology_transfer\tSwitch ID\tStatus\tInteraction ID\tIntramolecular\tID A\tBindingsite A ID\tBindingsite A Start\tBindingsite A End\tID B\tBindingsite B ID\tBindingsite B Start\tBindingsite B End\tAffected interactor\tSwitch type\tSwitch subtype\tSwitch mechanism\tSwitch direction\tSwitch outcome direction\tSwitch outcome\tModification\tModification sites\tModifying enzymes\tEffector\tCell cycle phase\tLocalisation\tPathway\tPMID" > elmswitches_mapped.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// ELM_CLASS_MAP  — Module 5n: ELM class definitions lookup table
//
// Parses elm_classes-*.tsv (shipped in legacy_data/elm/) into a flat TSV
// for ElmProteomeClassMatch Django model. Output is a per-run lookup table
// (not protein-specific) stored in mapped/annotations/.
// ──────────────────────────────────────────────────────────────────────────
process ELM_CLASS_MAP {
    tag  { "elm_class_map" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path elm_classes_tsv

    output:
    path "elm_classes.tsv", emit: elm_classes

    script:
    """
    create_elm_class_worker.py \\
        --elm_classes ${elm_classes_tsv} \\
        --outdir      .
    """

    stub:
    """
    echo -e "elm_accession\telm_identifier\tfunctional_site_name\tdescription\tregex\tprobability\tn_instances\tn_instances_in_pdb\telm_type" > elm_classes.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// PARSE_UNIPROT_DAT — extract feature TSVs from the Swiss-Prot flat file +
// InterPro protein2ipr file.  Runs once per accession set.
// Outputs consumed by ANNOTATION_MAP in place of per-protein REST API calls.
// Only invoked when flat-file mode is active (see main.nf bulk-mode selection).
// ──────────────────────────────────────────────────────────────────────────
process PARSE_UNIPROT_DAT {
    tag  { "parse_uniprot_dat" }
    label 'process_low'
    // No storeDir: output varies by which accessions are in loc_chrom (RAF1-only
    // vs full proteome).  A fixed storeDir path would reuse a single-gene result
    // for a full-proteome run.  The dat.gz downloads ARE storeDir-cached upstream
    // (FETCH_UNIPROT_SPROT_DAT / FETCH_INTERPRO_PFAM), so re-parsing is just a
    // fast streaming filter over already-local files.

    input:
    path uniprot_dat          // uniprot_sprot.dat.gz
    path interpro_pfam        // protein2ipr.dat.gz (or NO_FILE)
    path accessions_tsv       // loc_chrom TSV — for accession filter

    output:
    path 'uniprot_features.tsv', emit: features
    path 'pfam_domains.tsv',     emit: pfam

    script:
    def ipr_arg   = (interpro_pfam.name != 'NO_FILE') ? "--interpro_pfam ${interpro_pfam}" : ""
    def acc_col   = "Accession"   // column in loc_chrom TSV that holds UniProt accession
    """
    # Extract the accession list from the sequence table so we only scan the human entries
    python3 -c "
import pandas as pd, sys
df = pd.read_csv('${accessions_tsv}', sep='\\t', usecols=['Entry_Isoform'])
# Strip isoform suffix to get canonical accession (P04049-2 → P04049)
accs = df['Entry_Isoform'].dropna().str.split('-').str[0].unique()
print('\\n'.join(accs))
" > accessions.txt

    parse_uniprot_dat_worker.py \\
        --uniprot_dat   ${uniprot_dat} \\
        ${ipr_arg} \\
        --accessions    accessions.txt \\
        --outdir        .
    """

    stub:
    """
    printf 'Accession\\tType\\tStart\\tEnd\\tNote\\tEvidence\\tLigand\\n' > uniprot_features.tsv
    printf 'Accession\\thmm_acc\\thmm_name\\tstart\\tend\\ttype\\n'        > pfam_domains.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// ANNOTATION_MAP
// ──────────────────────────────────────────────────────────────────────────
process ANNOTATION_MAP {
    tag  { "annotation_map" }
    label 'process_medium'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/intermediate/annotations"
                                : "${params.outdir}/intermediate/annotations" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path elm_tsv,           stageAs: 'elm_in.tsv'
    path dibs_tsv,          stageAs: 'dibs_in'
    path mfib_tsv,          stageAs: 'mfib_in'
    path phasepro_tsv,      stageAs: 'phasepro_in'
    path uniprot_features,  stageAs: 'uniprot_features_in.tsv'   // from PARSE_UNIPROT_DAT or NO_FILE
    path pfam_dat,          stageAs: 'pfam_dat_in.tsv'           // from PARSE_UNIPROT_DAT or NO_FILE

    output:
    path "elm.tsv",               emit: elm_mapped
    path "dibs.tsv",              emit: dibs
    path "mfib.tsv",              emit: mfib
    path "phasepro.tsv",          emit: phasepro
    path "uniprot_roi.tsv",       emit: uniprot_roi
    path "uniprot_binding.tsv",   emit: uniprot_binding
    path "ptm_merged.tsv",        emit: ptm
    path "pfam_domains.tsv",      emit: pfam
    path "annotation_stats.tsv",  emit: stats

    script:
    // Pre-parsed bulk files take priority over per-protein REST API calls.
    // --uniprot_features_tsv and --pfam_tsv make the worker skip all API calls
    // and do local joins instead.  Fallback to the old REST API paths only when
    // both bulk files are NO_FILE (e.g., for small offline tests).
    def skip_uni  = (params.skip_uniprot_api || uniprot_features.name != 'NO_FILE') ? "--skip_uniprot" : ""
    def skip_pfam = (params.skip_pfam_api    || pfam_dat.name         != 'NO_FILE') ? "--skip_pfam"    : ""
    def uni_arg   = (uniprot_features.name   != 'NO_FILE') ? "--uniprot_features_tsv uniprot_features_in.tsv" : ""
    def pfam_arg  = (pfam_dat.name           != 'NO_FILE') ? "--pfam_tsv pfam_dat_in.tsv" : ""
    def elm_arg   = (elm_tsv.name      != 'NO_FILE') ? "--elm_tsv elm_in.tsv"        : ""
    def dibs_arg  = (dibs_tsv.name     != 'NO_FILE') ? "--dibs_tsv dibs_in"          : ""
    def mfib_arg  = (mfib_tsv.name     != 'NO_FILE') ? "--mfib_tsv mfib_in"          : ""
    def pp_arg    = (phasepro_tsv.name != 'NO_FILE') ? "--phasepro_tsv phasepro_in"  : ""
    def ptmdb_arg = params.legacy_ptmdb_dir  ? "--ptmdb_dir ${params.legacy_ptmdb_dir}"   : ""
    def ptmphs_arg= params.legacy_ptmphs_dir ? "--ptmphs_dir ${params.legacy_ptmphs_dir}" : ""
    """
    create_annotation_worker.py \\
        --loc_chrom     ${loc_chrom} \\
        ${elm_arg} \\
        ${dibs_arg} ${mfib_arg} ${pp_arg} \\
        ${ptmdb_arg} ${ptmphs_arg} \\
        ${uni_arg} ${pfam_arg} \\
        --output_dir    . \\
        --request_delay ${params.annotation_api_delay ?: 0.5} \\
        ${skip_uni} ${skip_pfam}
    """

    stub:
    """
    for f in elm.tsv dibs.tsv mfib.tsv phasepro.tsv \\
              uniprot_roi.tsv uniprot_binding.tsv ptm_merged.tsv \\
              pfam_domains.tsv annotation_stats.tsv; do
        touch "\$f"
    done
    """
}


// ──────────────────────────────────────────────────────────────────────────
// DISORDER_MAP  — local IUPred3/AIUPred libs + AlphaFold pLDDT API
// ──────────────────────────────────────────────────────────────────────────
// ──────────────────────────────────────────────────────────────────────────
// SPLIT_SEQ_TABLE  — split the sequence table into N gene-balanced chunks so the
// per-isoform heavy steps (DISORDER_MAP, COILEDCOILS_MAP) run as K concurrent
// tasks (maxForks / SLURM). Every isoform of a gene stays in one chunk.
// Only used when params.scatter_chunks > 1.
// ──────────────────────────────────────────────────────────────────────────
process SPLIT_SEQ_TABLE {
    tag  { "split_seq_table" }
    label 'process_low'

    input:
    path loc_chrom
    val  n_chunks

    output:
    path "chunk_*.tsv", emit: chunks

    script:
    """
    split_seq_table.py \\
        --loc_chrom ${loc_chrom} \\
        --n_chunks  ${n_chunks} \\
        --prefix    chunk_ \\
        --outdir    .
    """

    stub:
    """
    cp ${loc_chrom} chunk_001.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// TRANSCRIPT_MAP  — Map UniProt-keyed annotations → Gencode transcripts
//
// For each annotation keyed by UniProt accession (Entry_Isoform):
//   1. Find the canonical (main_isoform=yes) sequence
//   2. For each other transcript (isoform) of the same gene:
//      a. Direct match:   annotation region = exact substring → copy as-is
//      b. Homology:       100% sequence match of region → copy with homology_transfer=True
//      c. No match:       skip (unmappable for this transcript)
// ──────────────────────────────────────────────────────────────────────────
process TRANSCRIPT_MAP {
    tag  { "transcript_map" }
    label 'process_medium'
    // Annotation outputs → mapped/annotations/  (everything except disorder files)
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy',
        saveAs: { fn -> fn.startsWith("CombinedDisorder") ? null : fn }
    )
    // Disorder outputs → mapped/disorder/
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disorder"
                                : "${params.outdir}/final/disorder" },
        mode: 'copy',
        saveAs: { fn -> fn.startsWith("CombinedDisorder") ? fn : null }
    )

    input:
    path loc_chrom
    path elm_tsv,          stageAs: 'ann_elm.tsv'
    path dibs_tsv,         stageAs: 'ann_dibs.tsv'
    path mfib_tsv,         stageAs: 'ann_mfib.tsv'
    path phasepro_tsv,     stageAs: 'ann_phasepro.tsv'
    path uniprot_roi_tsv,  stageAs: 'ann_roi.tsv'
    path uniprot_bind_tsv, stageAs: 'ann_bind.tsv'
    path ptm_tsv,          stageAs: 'ann_ptm.tsv'
    path pfam_tsv,         stageAs: 'ann_pfam.tsv'
    path disorder_tsv,     stageAs: 'ann_disorder.tsv'
    path disorder_pos_tsv, stageAs: 'ann_disorder_pos.tsv'

    output:
    path "elm.tsv",              emit: elm
    path "dibs.tsv",             emit: dibs
    path "mfib.tsv",             emit: mfib
    path "phasepro.tsv",         emit: phasepro
    path "uniprot_roi.tsv",      emit: roi
    path "uniprot_binding.tsv",  emit: bind
    path "ptm_merged.tsv",       emit: ptm
    path "pfam_domains.tsv",     emit: pfam
    path "CombinedDisorderNew.tsv",     emit: disorder_regions
    path "CombinedDisorderNew_Pos.tsv", emit: disorder_pos
    path "transcript_map_stats.tsv",    emit: stats

    script:
    def only_main = params.only_main_isoforms ? "--only_main_isoforms" : ""
    """
    create_transcript_map_worker.py \\
        --loc_chrom    ${loc_chrom} \\
        --elm          ann_elm.tsv \\
        --dibs         ann_dibs.tsv \\
        --mfib         ann_mfib.tsv \\
        --phasepro     ann_phasepro.tsv \\
        --uniprot_roi  ann_roi.tsv \\
        --uniprot_bind ann_bind.tsv \\
        --ptm          ann_ptm.tsv \\
        --pfam         ann_pfam.tsv \\
        --disorder     ann_disorder.tsv \\
        --disorder_pos ann_disorder_pos.tsv \\
        --output_dir   . \\
        ${only_main}
    """

    stub:
    """
    for f in elm.tsv dibs.tsv mfib.tsv phasepro.tsv \\
              uniprot_roi.tsv uniprot_binding.tsv ptm_merged.tsv \\
              pfam_domains.tsv \\
              CombinedDisorderNew.tsv CombinedDisorderNew_Pos.tsv \\
              transcript_map_stats.tsv; do
        touch "\$f"
    done
    """
}


// ──────────────────────────────────────────────────────────────────────────
// ISOFORM_ALIGN_MAP  — Insertion-free isoform alignment
//
// Aligns every alternative isoform to the main (canonical) isoform using
// global pairwise alignment (Needleman-Wunsch, BLOSUM62).  Positions where
// the main isoform has a gap are stripped, yielding an insertion-free
// sequence of length == len(main isoform) — the same representation used
// by GOPHER for ortholog alignments.
//
// Output: mapped/sequence/isoform_alignment.tsv
//   Protein_ID  alt_Protein_ID  gene  main_seq_len  sequence
// ──────────────────────────────────────────────────────────────────────────
process ISOFORM_ALIGN_MAP {
    tag  { "isoform_align" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/sequence"
                                : "${params.outdir}/final/sequence" },
        mode: 'copy'
    )

    input:
    path loc_chrom_seq    // loc_chrom_with_names_isoforms_with_seq.tsv (main isoforms)
    path isoforms_only    // loc_chrom_with_names_isoforms_only.tsv (or NO_FILE)

    output:
    path "isoform_alignment.tsv", emit: isoform_alignment

    script:
    """
    create_isoform_align_worker.py \\
        --seq_table      ${loc_chrom_seq} \\
        --isoforms_table ${isoforms_only} \\
        --outdir         .
    """

    stub:
    """
    echo -e "Protein_ID\talt_Protein_ID\tgene\tmain_seq_len\tsequence" > isoform_alignment.tsv
    echo -e "RAF1-201\tRAF1-201\tRAF1\t648\tMANTIQQFLK..." >> isoform_alignment.tsv
    echo -e "RAF1-201\tRAF1-205\tRAF1\t648\tMANTIQQFLK..." >> isoform_alignment.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// HOMOLOGY_MANIFEST — audit table of homology-similarity transferred rows
// Scans mapped annotation TSVs carrying mapping_type / homology_transfer and
// records every row transferred from a main isoform onto an alternative one.
// Python worker: create_homology_manifest_worker.py
// ──────────────────────────────────────────────────────────────────────────
process HOMOLOGY_MANIFEST {
    tag  { "homology_manifest" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                : "${params.outdir}/final/annotations" },
        mode: 'copy'
    )

    input:
    path manifest_inputs   // collected mapped annotation TSVs

    output:
    path "homology_similarity_manifest.tsv", emit: manifest

    script:
    def in_args = (manifest_inputs instanceof List ? manifest_inputs : [manifest_inputs])
                      .collect { "${it}" }.join(' ')
    """
    create_homology_manifest_worker.py \\
        --inputs ${in_args} \\
        --outdir .
    """

    stub:
    """
    printf 'annotation\\tProtein_ID\\tsource_accession\\tidentifier\\tstart\\tend\\tposition\\tmapping_type\\n' > homology_similarity_manifest.tsv
    """
}
