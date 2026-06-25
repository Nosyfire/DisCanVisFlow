#!/usr/bin/env nextflow
/*
 * DisCanVis Nextflow Pipeline  (DSL2)
 *
 * Full BLAST → ID-map DAG:
 *
 *   [FTP or local] UniProt FASTA ──┬──► SUBSET_FASTA(UniProt) ──► MAKEBLASTDB(UniProt) ─┐
 *                                  │                                                      │
 *                                  └──────────────────────────────► BLASTP (uni→gen DB) ──┤
 *                                                                                          │
 *   [FTP or local] GENCODE FASTA ──┬──► SUBSET_FASTA(GENCODE) ──► MAKEBLASTDB(GENCODE) ──┤
 *                                  │                                                      │
 *                                  └──────────────────────────────► BLASTP (gen→uni DB) ──┘
 *                                                                           │
 *                                                              MERGE_BLAST_HITS
 *                                                                           │
 *                                                                        ID_MAP
 *                                                                           │
 *                                                   bestmaps_blast_gene_transcript.tsv
 *
 * FASTA sourcing (in priority order)
 * ───────────────────────────────────
 *   1. --uniprot_fasta / --gencode_fasta  CLI / config path  → used as-is
 *   2. params are null                                        → auto-download from FTP
 *      Files are cached in params.ref_dir via storeDir;
 *      subsequent runs skip the download entirely.
 *
 * Run modes
 * ─────────
 *   # Cellular-vulnerability project on an 8 GB laptop
 *   nextflow run main.nf --project cellular_vulnerability --machine laptop -resume
 *
 *   # Full DisCanVis update on a stronger workstation/server
 *   nextflow run main.nf --project full_discanvis --machine hard -resume
 *
 *   # Provide local FASTAs explicitly
 *   nextflow run main.nf \
 *       --uniprot_fasta /path/to/uniprot_swissprot.fasta \
 *       --gencode_fasta /path/to/gencode_pc_translations.fasta
 *
 *   # DAG-only stub validation (no BLAST, no download)
 *   nextflow run main.nf --project test_one_protein --data local --machine laptop -stub
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Module imports
// ---------------------------------------------------------------------------
include { SETUP_DEPS                 } from './modules/setup'
include { FETCH_UNIPROT_FASTA;
          FETCH_GENCODE_TRANSLATIONS;
          FETCH_GENCODE_TRANSCRIPTS;
          FETCH_GENCODE_GTF;
          FETCH_GENCODE_REFSEQ;
          FETCH_MONDO;
          FETCH_ALPHAMISSENSE;
          DECOMPRESS_ALPHAMISSENSE;
          FETCH_INTACT;
          FETCH_BIOGRID;
          FETCH_HIPPIE;
          FETCH_HG38_2BIT;
          FETCH_DBSNP_BB;
          FETCH_UNIPROT_ISOFORMS;
          FETCH_SIFTS;
          FETCH_ALPHAFOLD_BULK;
          FETCH_MAVEDB;
          FETCH_PROTEINGYM;
          FETCH_DEPMAP;
          NORMALISE_DEPMAP_MANUAL;
          FETCH_DBNSFP;
          FETCH_OMIM;
          FETCH_CBIOPORTAL;
          FETCH_ZENODO;
          FETCH_UNIPROT_SPROT_DAT;
          FETCH_INTERPRO_PFAM      } from './modules/fetch_references'

// Each process that is called more than once in the same workflow must be
// imported under a unique alias (DSL2 restriction).
include { SUBSET_FASTA  as SUBSET_UNIPROT;
          SUBSET_FASTA  as SUBSET_GENCODE;
          SUBSET_FASTA  as SUBSET_CDNA;
          MAKEBLASTDB   as MAKEBLASTDB_UNIPROT;
          MAKEBLASTDB   as MAKEBLASTDB_GENCODE;
          BLASTP        as BLASTP_GEN2UNI;
          BLASTP        as BLASTP_UNI2GEN;
          MERGE_UNIPROT_ISOFORMS;
          MERGE_BLAST_HITS         } from './modules/blast_search'

include { ID_MAP                   } from './modules/blast_mapping'
include { SEQUENCE_PROCESS         } from './modules/sequence_process'
include { SPLIT_CDNA_FASTA;
          BLAT_ALIGN;
          MERGE_BLAT_PSL;
          GENOME_MAP;
          GENOME_QUERY_MAP         } from './modules/genome_mapping'
include { FETCH_CLINVAR;
          MUTATION_MAP;
          MUTATION_MAP as MUTATION_MAP_CLINVAR;
          MUTATION_MAP as MUTATION_MAP_TCGA;
          MUTATION_MAP as MUTATION_MAP_CBIOPORTAL } from './modules/mutation_mapping'
include { FETCH_ELM;
          FETCH_MOBIDB;
          FETCH_GO;
          ELM_CLASS_MAP;
          ELM_SWITCHES_MAP;
          MOBIDB_MAP;
          ANNOTATION_MAP;
          SPLIT_SEQ_TABLE;
          DISORDER_MAP;
          GO_MAP;
          POLYMORPHISM_MAP;
          PEM_MAP;
          COILEDCOILS_MAP;
          PPI_PREPROCESS;
          PPI_MAP;
          CONSERVATION_MAP;
          GOPHER_RECOMPUTE;
          SCANSITE_MAP;
          CLINVAR_DISEASE_MAP;
          CLINVAR_DISEASE_BUILD;
          OMIM_MAP;
          CANCER_DRIVER_MAP;
          ALPHAMISSENSE_MAP;
          MAVEDB_MAP;
          PROTEINGYM_MAP;
          DEPMAP_MAP;
          DBNSFP_MAP;
          PATHOGENICITY_MAP;
          POSITION_BASED_MAP;
          ISOFORM_ALIGN_MAP;
          PEM_TRANSFER_MAP;
          PDB_MAP;
          PDB_BULK_MAP;
          FINCHES_MAP;
          EXON_MAP;
          TRANSCRIPT_MAP;
          HOMOLOGY_MANIFEST;
          MAPPING_REPORT;
          PARSE_UNIPROT_DAT;
          PARSE_ALPHAFOLD_PLDDT    } from './modules/annotation_mapping'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    // ── gene_list_file → target_gene resolution ───────────────────────────────
    // If --gene_list_file is set, read the file (one gene per line) and use it
    // as the target_gene filter, overriding any explicit --target_gene value.
    if ( params.gene_list_file ) {
        def glf = file(params.gene_list_file, checkIfExists: true)
        def genes_from_file = glf.readLines().findAll { it.trim() && !it.startsWith('#') }*.trim().join(',')
        params.target_gene = genes_from_file
        log.info "Gene list file : ${params.gene_list_file} → ${genes_from_file.split(',').size()} genes"
    }

    // ── Manual-reference preflight ────────────────────────────────────────────
    // Some public/licensed sources cannot always be fetched from a non-browser
    // process. Create conventional reference folders and stop early with one
    // actionable checklist when an enabled manual reference is missing.
    new File(params.ref_dir.toString()).mkdirs()
    [
        "${params.ref_dir}/depmap",
        "${params.ref_dir}/dbnsfp",
        "${params.ref_dir}/phastcons"
    ].each { new File(it.toString()).mkdirs() }

    def manual_missing = []
    def depmap_manual_ch = null
    def depmap_raw_csv = params.depmap_raw_csv ?: "${params.ref_dir}/depmap/OmicsSomaticMutations.csv"

    if ( !params.skip_depmap && !params.fetch_depmap ) {
        def depmap_tsv_f = params.depmap_tsv ? file(params.depmap_tsv) : null
        def depmap_raw_f = file(depmap_raw_csv)
        def have_depmap_tsv = depmap_tsv_f && depmap_tsv_f.exists() && depmap_tsv_f.toFile().length() > 0
        def have_depmap_raw = depmap_raw_f.exists() && depmap_raw_f.toFile().length() > 0

        if ( have_depmap_tsv ) {
            log.info "Manual reference confirmed: DepMap TSV → ${depmap_tsv_f}"
            depmap_manual_ch = Channel.value(depmap_tsv_f)
        } else if ( have_depmap_raw ) {
            log.info "Manual reference confirmed: DepMap raw CSV → ${depmap_raw_f}"
            log.info "DepMap normalized TSV missing; pipeline will create ${params.depmap_tsv}"
            depmap_manual_ch = NORMALISE_DEPMAP_MANUAL( Channel.value(depmap_raw_f) ).tsv
        } else {
            manual_missing << """
DepMap mutations
  Needed because: --skip_depmap false and --fetch_depmap false
  Download page: https://depmap.org/portal/download/all/
  File to choose: OmicsSomaticMutations.csv
  Copy to: ${depmap_raw_csv}
  Pipeline output after normalization: ${params.depmap_tsv ?: "${params.ref_dir}/depmap/depmap_mutations_raw.tsv"}
  To skip this annotation: --skip_depmap true
"""
        }
    }

    def genome_enabled = (params.hg38_2bit || params.fetch_hg38_2bit) as boolean

    if ( !params.skip_pathogenicity && !params.fetch_dbnsfp &&
         !params.dbnsfp_raw_dir && !params.dbnsfp_tsv ) {
        manual_missing << """
dbNSFP pathogenicity
  Needed because: --skip_pathogenicity false and --fetch_dbnsfp false
  Download/setup: obtain dbNSFP under its academic terms
  Copy raw chr*.gz files under: ${params.ref_dir}/dbnsfp/
  Then set: --dbnsfp_raw_dir ${params.ref_dir}/dbnsfp
  Alternative: set --dbnsfp_tsv /path/to/pre_mapped_dbnsfp.tsv
  To skip this annotation: --skip_pathogenicity true
"""
    }

    if ( genome_enabled && !params.skip_conservation && !params.skip_phastcons &&
         !params.phastcons_dir ) {
        manual_missing << """
phastCons conservation
  Needed because: genome mapping and conservation are enabled, but no --phastcons_dir was provided
  Copy hg38 per-chromosome bigWig files (chr1.bw ... chrY.bw) under: ${params.ref_dir}/phastcons/
  Then set: --phastcons_dir ${params.ref_dir}/phastcons
  To skip conservation entirely: --skip_conservation true
  To skip only phastCons: --skip_phastcons true
"""
    }

    if ( manual_missing ) {
        error """
Manual reference files are required before this run can continue.

The reference folders have been created under:
  ${params.ref_dir}

Missing manual references:
${manual_missing.join('\n')}
After copying/downloading the files, rerun the same command with -resume.
"""
    }

    // ── FASTA sourcing ───────────────────────────────────────────────────────
    // When fasta params are null, download from FTP (cached in ref_dir).
    // When paths are provided, use them directly.

    if ( params.uniprot_fasta ) {
        log.info "UniProt  FASTA : ${params.uniprot_fasta}  (local)"
        uniprot_ch = Channel.fromPath( params.uniprot_fasta, checkIfExists: true )
    } else {
        log.info "UniProt  FASTA : not provided → downloading SwissProt from FTP"
        log.info "  URL  : ${params.uniprot_ftp_url}"
        log.info "  Cache: ${params.ref_dir}/uniprot/uniprot_swissprot.fasta"
        uniprot_ch = FETCH_UNIPROT_FASTA().fasta
    }

    // ── Mapping mode ─────────────────────────────────────────────────────────
    // 'all_isoform_mapping' appends curated Swiss-Prot isoforms to the UniProt
    // BLAST DB so each GENCODE transcript can be paired to its true isoform.
    // 'main_isoform_mapping' (default) BLASTs against the canonical proteome only.
    if ( params.mapping_mode == 'all_isoform_mapping' ) {
        log.info "Mapping mode   : all_isoform_mapping (curated Swiss-Prot isoforms)"
        // local 'additional' FASTA if given, else download from UniProt FTP
        def iso_fa = params.uniprot_isoform_fasta
            ? Channel.fromPath( params.uniprot_isoform_fasta, checkIfExists: true )
            : FETCH_UNIPROT_ISOFORMS().fasta
        uniprot_ch = MERGE_UNIPROT_ISOFORMS( uniprot_ch, iso_fa ).fasta
    } else {
        log.info "Mapping mode   : main_isoform_mapping (canonical Swiss-Prot only)"
    }

    if ( params.gencode_fasta ) {
        log.info "GENCODE  FASTA : ${params.gencode_fasta}  (local)"
        gencode_ch = Channel.fromPath( params.gencode_fasta, checkIfExists: true )
    } else {
        log.info "GENCODE  FASTA : not provided → downloading translations from FTP"
        log.info "  URL  : ${params.gencode_translations_url}"
        log.info "  Cache: ${params.ref_dir}/gencode/gencode_pc_translations.fasta"
        gencode_ch = FETCH_GENCODE_TRANSLATIONS().fasta
    }

    // ── Pre-fetch Module 2 reference files ──────────────────────────────────
    // GTF, cDNA transcripts, and RefSeq metadata are needed for the Sequence
    // Process step (Module 2).  We fetch them now so the reference directory
    // is fully populated for downstream modules.
    // storeDir means this is a no-op when the files already exist.
    //
    // When local paths are provided via params, use those instead.

    if ( params.gencode_gtf ) {
        log.info "GENCODE  GTF   : ${params.gencode_gtf}  (local)"
        gencode_gtf_ch = Channel.fromPath( params.gencode_gtf, checkIfExists: true )
    } else {
        log.info "GENCODE  GTF   : not provided → downloading from FTP (cached for Module 2)"
        gencode_gtf_ch = FETCH_GENCODE_GTF().gtf
    }

    if ( params.gencode_transcripts ) {
        log.info "GENCODE  cDNA  : ${params.gencode_transcripts}  (local)"
        gencode_transcripts_ch = Channel.fromPath( params.gencode_transcripts, checkIfExists: true )
    } else {
        log.info "GENCODE  cDNA  : not provided → downloading pc_transcripts from FTP (cached for Module 2)"
        gencode_transcripts_ch = FETCH_GENCODE_TRANSCRIPTS().fasta
    }

    if ( params.gencode_refseq ) {
        log.info "GENCODE RefSeq : ${params.gencode_refseq}  (local)"
        gencode_refseq_ch = Channel.fromPath( params.gencode_refseq, checkIfExists: true )
    } else {
        log.info "GENCODE RefSeq : not provided → downloading metadata from FTP (cached for Module 2)"
        gencode_refseq_ch = FETCH_GENCODE_REFSEQ().tsv
    }

    // Reference files are ready in channels for downstream modules.
    // gencode_gtf_ch, gencode_transcripts_ch, gencode_refseq_ch are used by
    // SEQUENCE_PROCESS and GENOME_MAP below.

    // ── Target-gene subsetting ───────────────────────────────────────────────
    // Resolve the target-gene search term (empty string = pass-through mode).
    // SUBSET_FASTA always runs; when target == '' it copies the full FASTA.
    // --target_gene null on the CLI arrives as the string "null" — treat it as empty.
    def target = (params.target_gene && params.target_gene.toString() != 'null') ? params.target_gene.toString() : ''
    if ( target ) {
        log.info "Single-gene mode: subsetting FASTAs to '${target}'"
    }

    // ── Step 0a: Subset FASTAs ───────────────────────────────────────────────
    uni_sub = SUBSET_UNIPROT( uniprot_ch, target, 'uniprot' )
    gen_sub = SUBSET_GENCODE( gencode_ch, target, 'gencode' )

    // ── Step 0b: Build BLAST databases ───────────────────────────────────────
    uni_db = MAKEBLASTDB_UNIPROT( uni_sub.subset, 'uniprot' )
    gen_db = MAKEBLASTDB_GENCODE( gen_sub.subset, 'gencode' )

    // ── Step 0c: Reciprocal BLAST ────────────────────────────────────────────
    //  blast1: query = GENCODE proteins,  database = UniProt
    //          → uniprotdb_gencode_query.xml
    blast1_xml = BLASTP_GEN2UNI(
        gen_sub.subset,
        uni_db.blastdb,
        'uniprotdb_gencode_query'
    )
    //  blast2: query = UniProt proteins,  database = GENCODE
    //          → gencodedb_uniprot_query.xml
    blast2_xml = BLASTP_UNI2GEN(
        uni_sub.subset,
        gen_db.blastdb,
        'gencodedb_uniprot_query'
    )

    // ── Step 0d: Parse & merge reciprocal hits ────────────────────────────────
    merged = MERGE_BLAST_HITS(
        blast1_xml.blast_xml,   // uniprotdb_gencode_query.xml
        blast2_xml.blast_xml    // gencodedb_uniprot_query.xml
    )

    // ── Step 1: Transcript → UniProt ID mapping ───────────────────────────────
    ID_MAP(
        merged.bestsequences,
        merged.isoforms
    )

    // ── Step 2: Sequence Process ──────────────────────────────────────────────
    // The isoforms output from ID_MAP is optional; fall back to NO_FILE sentinel.
    // gen_sub.subset is used here (the subsetted translations FASTA) — for RAF1
    // mode it contains only RAF1 proteins (matching BLAST hits); for the full run
    // it is the complete translations FASTA (pass-through from SUBSET_GENCODE).
    isoforms_input_ch = ID_MAP.out.isoforms
        .ifEmpty( file("${projectDir}/assets/NO_FILE") )

    SEQUENCE_PROCESS(
        ID_MAP.out.id_map,      // bestmaps_blast_gene_transcript.tsv
        isoforms_input_ch,      // blastmaps_isoforms.tsv (or NO_FILE sentinel)
        gencode_gtf_ch,
        gen_sub.subset          // pc_translations FASTA (subsetted or full)
    )

    SEQUENCE_PROCESS.out.loc_chrom_seq.view { f ->
        "\n✔  Transcript sequence table: ${f}\n"
    }

    // ── Isoform alignment (insertion-free) ───────────────────────────────────
    // Runs immediately after SEQUENCE_PROCESS; no genome mapping required.
    def isoforms_only_ch = SEQUENCE_PROCESS.out.isoforms_only
        .ifEmpty( file("${projectDir}/assets/NO_FILE") )
    ISOFORM_ALIGN_MAP(
        SEQUENCE_PROCESS.out.loc_chrom_seq,
        isoforms_only_ch
    )
    ISOFORM_ALIGN_MAP.out.isoform_alignment.view { f ->
        "\n✔  Insertion-free isoform alignment: ${f}\n"
    }

    // ── Step 3: Genome Mapping ────────────────────────────────────────────────
    // Requires hg38.2bit: a local --hg38_2bit path, or fetched when --fetch_hg38_2bit.
    // Skipped entirely when neither is set.
    if ( genome_enabled ) {
        // Value channel: can be read by both BLAT_ALIGN and GENOME_MAP
        hg38_ch = params.hg38_2bit
            ? Channel.value( file(params.hg38_2bit, checkIfExists: true) )
            : FETCH_HG38_2BIT().twobit

        // Subset cDNA FASTA to target gene (same logic as for protein FASTA)
        cdna_sub = SUBSET_CDNA( gencode_transcripts_ch, target, 'cdna' )

        // BLAT: align cDNA against hg38 genome → PSL file.
        // blat_chunks > 1 → split cDNA FASTA into N chunks, run N BLAT jobs in
        // parallel, merge raw PSLs, then run pslCDnaFilter once on the merged file
        // (-bestOverlap needs all hits for a query to be visible together).
        // blat_chunks == 1 (default) → single BLAT job, still goes through
        // MERGE_BLAT_PSL so pslCDnaFilter always runs in one dedicated process.
        def n_blat   = Math.max(1, (params.blat_chunks ?: 1) as int)
        def blat_in  = ( n_blat > 1 )
            ? SPLIT_CDNA_FASTA( cdna_sub.subset, n_blat ).chunks.flatten()
            : cdna_sub.subset
        def raw_psls = BLAT_ALIGN( blat_in, hg38_ch ).psl.collect()
        blat_out     = MERGE_BLAT_PSL( raw_psls )

        // Genome coordinate map: PSL + FASTAs + hg38.2bit → combined_map.map
        GENOME_MAP(
            blat_out.psl,
            cdna_sub.subset,
            gen_sub.subset,         // protein translations (subsetted)
            SEQUENCE_PROCESS.out.loc_chrom,
            hg38_ch
        )

        GENOME_MAP.out.map_file.view { f ->
            "\n✔  Genome coordinate map: ${f}\n"
        }

        // ── Step 4: Mutation Mapping ──────────────────────────────────────────
        // ClinVar + TCGA + cBioPortal MAFs run in parallel (each → mutations/<source>/).
        // ClinVar outputs feed CLINVAR_DISEASE_BUILD via MUTATION_MAP_CLINVAR alias.

        def cancer_mut_specs = []
        if ( params.mutation_maf ) {
            cancer_mut_specs << [file(params.mutation_maf, checkIfExists: true),
                                 params.mutation_source ?: 'custom', 'maf']
        }
        if ( params.mutation_vcf ) {
            cancer_mut_specs << [file(params.mutation_vcf, checkIfExists: true),
                                 params.mutation_source ?: 'custom', 'vcf']
        }

        if ( params.tcga_maf ) {
            MUTATION_MAP_TCGA(
                GENOME_MAP.out.map_file,
                SEQUENCE_PROCESS.out.loc_chrom_seq,
                file(params.tcga_maf, checkIfExists: true),
                Channel.value('TCGA'),
                Channel.value('maf')
            )
            MUTATION_MAP_TCGA.out.stats.view { f ->
                "\n✔  Mutation mapping stats (TCGA): ${f}\n"
            }
        }

        // cBioPortal: open datahub study bundle (FETCH_CBIOPORTAL) or a local MAF.
        if ( params.fetch_cbioportal || params.cbioportal_maf ) {
            def cbio_maf = params.fetch_cbioportal
                ? FETCH_CBIOPORTAL().maf
                : file(params.cbioportal_maf, checkIfExists: true)
            MUTATION_MAP_CBIOPORTAL(
                GENOME_MAP.out.map_file,
                SEQUENCE_PROCESS.out.loc_chrom_seq,
                cbio_maf,
                Channel.value('CBioportal'),
                Channel.value('maf')
            )
            MUTATION_MAP_CBIOPORTAL.out.stats.view { f ->
                "\n✔  Mutation mapping stats (CBioportal): ${f}\n"
            }
        }

        if ( cancer_mut_specs ) {
            cancer_mut_ch = Channel.fromList(cancer_mut_specs).multiMap { spec ->
                mutation_file: spec[0]
                source:        spec[1]
                input_format:  spec[2]
            }
            MUTATION_MAP(
                GENOME_MAP.out.map_file,
                SEQUENCE_PROCESS.out.loc_chrom_seq,
                cancer_mut_ch.mutation_file,
                cancer_mut_ch.source,
                cancer_mut_ch.input_format
            )
            MUTATION_MAP.out.stats.view { f ->
                "\n✔  Mutation mapping stats: ${f}\n"
            }
        }

        if ( params.clinvar_vcf ) {
            MUTATION_MAP_CLINVAR(
                GENOME_MAP.out.map_file,
                SEQUENCE_PROCESS.out.loc_chrom_seq,
                file(params.clinvar_vcf, checkIfExists: true),
                Channel.value('ClinVar'),
                Channel.value('clinvar_vcf')
            )
            MUTATION_MAP_CLINVAR.out.stats.view { f ->
                "\n✔  Mutation mapping stats (ClinVar): ${f}\n"
            }
        } else if ( !params.tcga_maf && !params.cbioportal_maf && !cancer_mut_specs ) {
            clinvar_ch = FETCH_CLINVAR()
            MUTATION_MAP_CLINVAR(
                GENOME_MAP.out.map_file,
                SEQUENCE_PROCESS.out.loc_chrom_seq,
                clinvar_ch.vcf,
                Channel.value('ClinVar'),
                Channel.value('clinvar_vcf')
            )
            MUTATION_MAP_CLINVAR.out.stats.view { f ->
                "\n✔  Mutation mapping stats (ClinVar): ${f}\n"
            }
        }

        // ── Step 8a (build): ClinVar disease from MONDO OBO + ClinVar mutations only ─
        // mondo_obo: use local path if provided, else download via FETCH_MONDO (storeDir cached)
        if ( !params.skip_clinvar_disease && params.clinvar_disease_from_mutations
             && (params.clinvar_vcf || (!params.tcga_maf && !params.cbioportal_maf && !cancer_mut_specs)) ) {
            def mondo_ch = params.mondo_obo
                ? Channel.value( file(params.mondo_obo, checkIfExists: true) )
                : FETCH_MONDO().obo
            CLINVAR_DISEASE_BUILD(
                SEQUENCE_PROCESS.out.loc_chrom_seq,
                mondo_ch,
                MUTATION_MAP_CLINVAR.out.missense,
                MUTATION_MAP_CLINVAR.out.frameshift,
                MUTATION_MAP_CLINVAR.out.nonsense,
                MUTATION_MAP_CLINVAR.out.indel
            )
            CLINVAR_DISEASE_BUILD.out.clinvar_disease.view { f ->
                "\n✔  ClinVar disease (MONDO build): ${f}\n"
            }
            CLINVAR_DISEASE_BUILD.out.clinvar_disease_mutations.view { f ->
                "\n✔  ClinVar disease + mutations (merged): ${f}\n"
            }
        }

        // ── Step 8f: dbNSFP pathogenicity (raw map via combined_map.map) ─────────
        // Raw chr*.gz dir from FETCH_DBNSFP (academic download) or a local path.
        if ( !params.skip_pathogenicity && (params.dbnsfp_raw_dir || params.fetch_dbnsfp) ) {
            def dbnsfp_dir = params.fetch_dbnsfp
                ? FETCH_DBNSFP().dir
                : Channel.value( file(params.dbnsfp_raw_dir, checkIfExists: false) )
            DBNSFP_MAP(
                SEQUENCE_PROCESS.out.loc_chrom_seq,
                GENOME_MAP.out.map_file,
                dbnsfp_dir
            )
            DBNSFP_MAP.out.scores.view { f ->
                "\n✔  Pathogenicity scores (raw dbNSFP map): ${f}\n"
            }
        }

        // ── Genome ↔ Protein reference tables (per-residue index + all-SNV reference) ──
        GENOME_QUERY_MAP( GENOME_MAP.out.map_file )
        GENOME_QUERY_MAP.out.index.view { f ->
            "\n✔  Genome→protein index: ${f}\n"
        }
        GENOME_QUERY_MAP.out.mutations.view { f ->
            "\n✔  Genome→protein mutations: ${f}\n"
        }

    } else {
        log.warn "params.hg38_2bit not set — skipping Modules 3, 4 (Genome + Mutation Mapping)"
    }

    // ── Step 5: Annotation Mapping ────────────────────────────────────────────
    // ELM / DIBS / MFIB / PhasePro come from legacy_data/ (project-local)
    // MobiDB is downloaded from API (FETCH_MOBIDB) unless local file provided.

    no_file = file("${projectDir}/assets/NO_FILE")

    // combined_map.map is only produced when hg38_2bit is set (Module 3).
    // Downstream annotators that need genomic↔protein mapping fall back to NO_FILE.
    combined_map_ch = genome_enabled
        ? GENOME_MAP.out.map_file
        : Channel.value(no_file)

    // ELM: always from legacy_data (params.elm_tsv is set in config)
    elm_tsv_ch = params.elm_tsv
        ? Channel.value( file(params.elm_tsv, checkIfExists: true) )
        : Channel.value( no_file )

    // MobiDB: download via FETCH_MOBIDB unless local file provided
    mobidb_tsv_ch = params.mobidb_tsv
        ? Channel.value( file(params.mobidb_tsv, checkIfExists: true) )
        : FETCH_MOBIDB().mobidb_tsv
    mobidb_disorder_input_ch = params.mobidb_tsv
        ? mobidb_tsv_ch
        : mobidb_tsv_ch.first()

    // DIBS / MFIB / PhasePro from legacy_data/ (auto-set in config)
    dibs_ch     = params.dibs_tsv     ? Channel.value( file(params.dibs_tsv,     checkIfExists: true) ) : Channel.value( no_file )
    mfib_ch     = params.mfib_tsv     ? Channel.value( file(params.mfib_tsv,     checkIfExists: true) ) : Channel.value( no_file )
    phasepro_ch = params.phasepro_tsv ? Channel.value( file(params.phasepro_tsv, checkIfExists: true) ) : Channel.value( no_file )

    // ── Bulk FTP pre-parse: UniProt features + Pfam domains ─────────────────
    // Replace ~37k per-protein UniProt REST + InterPro REST calls in ANNOTATION_MAP
    // with a single streaming parse of two FTP flat files (storeDir-cached).
    // When params.uniprot_dat_gz / interpro_pfam_dat_gz are set, use local files;
    // otherwise auto-download from UniProt/EBI FTP (params.fetch_uniprot_dat true).
    def _uni_dat_file = params.uniprot_dat_gz
        ? Channel.value( file(params.uniprot_dat_gz) )
        : ( params.fetch_uniprot_dat != false ? FETCH_UNIPROT_SPROT_DAT().dat : Channel.value(no_file) )
    def _ipr_dat_file = params.interpro_pfam_dat_gz
        ? Channel.value( file(params.interpro_pfam_dat_gz) )
        : ( params.fetch_interpro_pfam != false ? FETCH_INTERPRO_PFAM().dat : Channel.value(no_file) )

    def uniprot_features_ch
    def pfam_bulk_ch
    if ( params.skip_uniprot_api && params.skip_pfam_api ) {
        // Both APIs disabled and no bulk file requested — pass NO_FILE
        uniprot_features_ch = Channel.value(no_file)
        pfam_bulk_ch        = Channel.value(no_file)
    } else {
        def parsed = PARSE_UNIPROT_DAT(
            _uni_dat_file,
            _ipr_dat_file,
            SEQUENCE_PROCESS.out.loc_chrom_seq
        )
        uniprot_features_ch = parsed.features
        pfam_bulk_ch        = parsed.pfam
    }

    ANNOTATION_MAP(
        SEQUENCE_PROCESS.out.loc_chrom_seq,
        elm_tsv_ch,
        dibs_ch,
        mfib_ch,
        phasepro_ch,
        uniprot_features_ch,
        pfam_bulk_ch
    )

    ANNOTATION_MAP.out.stats.view { f ->
        "\n✔  Annotation mapping stats: ${f}\n"
    }

    // ── Step 5o: MobiDB disorder features (MobiDBDisorder) ───────────────────
    MOBIDB_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq, mobidb_tsv_ch )
    MOBIDB_MAP.out.mobidb_disorder.view { f ->
        "\n✔  MobiDB disorder features: ${f}\n"
    }

    // ── Step 5n: ELM class lookup table (ElmProteomeClassMatch) ──────────────
    def elm_classes_f = params.elm_classes_tsv
        ? file(params.elm_classes_tsv, checkIfExists: true)
        : file("${projectDir}/legacy_data/elm/elm_classes-2025.tsv")
    ELM_CLASS_MAP( Channel.value(elm_classes_f) )
    ELM_CLASS_MAP.out.elm_classes.view { f ->
        "\n✔  ELM class definitions: ${f}\n"
    }

    // ── Step 5p: ELM molecular switches (Elm_Switches) ───────────────────────
    def elm_switches_f = params.elm_switches_tsv
        ? file(params.elm_switches_tsv, checkIfExists: true)
        : file("${projectDir}/legacy_data/elm/elmswitches-2023.tsv")
    ELM_SWITCHES_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq, Channel.value(elm_switches_f) )
    ELM_SWITCHES_MAP.out.elm_switches.view { f ->
        "\n✔  ELM switches mapped: ${f}\n"
    }

    // ── MaveDB single-mutant functional scores ──────────────────────────────
    // Source priority: fresh fetch (api.mavedb.org, uniprot mode) → --mavedb_raw
    // (uniprot mode) → --mavedb_tsv (pre-mapped, premapped mode).
    def mavedb_enabled = !params.skip_mavedb &&
                         (params.fetch_mavedb || params.mavedb_raw || params.mavedb_tsv)
    if ( mavedb_enabled ) {
        def mavedb_src = params.fetch_mavedb ? FETCH_MAVEDB().raw
                       : params.mavedb_raw   ? Channel.value( file(params.mavedb_raw, checkIfExists: true) )
                       :                        Channel.value( file(params.mavedb_tsv, checkIfExists: true) )
        MAVEDB_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq, mavedb_src )
        MAVEDB_MAP.out.mavedb.view { f ->
            "\n✔  MaveDB functional scores: ${f}\n"
        }
    }

    // ── ProteinGym DMS scores + pathogenicity proxy ─────────────────────────
    def proteingym_enabled = !params.skip_proteingym &&
                             (params.fetch_proteingym || params.proteingym_raw || params.proteingym_tsv)
    if ( proteingym_enabled ) {
        def proteingym_src = params.fetch_proteingym ? FETCH_PROTEINGYM().raw
                           : params.proteingym_raw   ? Channel.value( file(params.proteingym_raw, checkIfExists: true) )
                           :                            Channel.value( file(params.proteingym_tsv, checkIfExists: true) )
        PROTEINGYM_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq, proteingym_src )
        PROTEINGYM_MAP.out.proteingym.view { f ->
            "\n✔  ProteinGym DMS scores: ${f}\n"
        }
    }

    // ── One-time machine setup: External_Programs + bigBedToBed ──────────────
    // SETUP_DEPS clones AIUPred, creates conda envs, and installs bigBedToBed on
    // first run. Outputs are storeDir-cached so it is skipped on every subsequent
    // run. When params.auto_setup is false (default for raf1/full profiles), we
    // pass NO_FILE sentinels and the three dependent processes start immediately.
    def setup_done_ch
    def setup_aiupred_py_ch
    def setup_deepcoil_py_ch
    if ( params.auto_setup ) {
        def _setup = SETUP_DEPS()
        setup_done_ch        = _setup.done
        setup_aiupred_py_ch  = _setup.aiupred_python
        setup_deepcoil_py_ch = _setup.deepcoil_python
    } else {
        setup_done_ch        = Channel.value(no_file)
        setup_aiupred_py_ch  = Channel.value(no_file)
        setup_deepcoil_py_ch = Channel.value(no_file)
    }

    // ── Step 5b: Disorder Mapping ─────────────────────────────────────────────
    // Uses local IUPred3/AIUPred libs (ext_programs param) + AlphaFold API.
    // MobiDB comes from FETCH_MOBIDB; Pfam from ANNOTATION_MAP.
    // ── Scatter the per-isoform heavy steps (DISORDER, COILEDCOILS) ───────────
    def scatter_n = Math.max(1, (params.scatter_chunks ?: 1) as int)
    def split_chunks = ( scatter_n > 1 )
        ? SPLIT_SEQ_TABLE( SEQUENCE_PROCESS.out.loc_chrom_seq, scatter_n ).chunks
        : null
    def disorder_loc_ch = ( scatter_n > 1 ) ? split_chunks.flatten()
                                            : SEQUENCE_PROCESS.out.loc_chrom_seq

    // Pre-extract AlphaFold pLDDT from bulk tar → local dict lookup (no EBI API per protein)
    def af_plddt_ch
    def _have_af_bulk = !params.skip_alphafold && (params.alphafold_tar || params.fetch_alphafold_bulk)
    if ( _have_af_bulk ) {
        def af_tar_ch = params.alphafold_tar
            ? Channel.value( file(params.alphafold_tar) )
            : FETCH_ALPHAFOLD_BULK().tar
        af_plddt_ch = PARSE_ALPHAFOLD_PLDDT( af_tar_ch ).plddt
    } else {
        af_plddt_ch = Channel.value(no_file)
    }

    DISORDER_MAP(
        disorder_loc_ch,
        mobidb_disorder_input_ch,
        ANNOTATION_MAP.out.pfam.first(),
        af_plddt_ch.first(),
        setup_done_ch.first(),
        setup_aiupred_py_ch.first()
    )

    // Merge per-chunk disorder tables. Genes never split across chunks, so this
    // concatenates whole-gene blocks → functionally identical to the single task.
    def _dis_pub = params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disorder"
                                   : "${params.outdir}/final/disorder"
    def _mergeDis = { ch, fname -> ( scatter_n > 1
        ? ch.collectFile(name: fname, keepHeader: true, skip: 1, sort: true, storeDir: _dis_pub)
        : ch ).first() }
    def disorder_iupred  = _mergeDis.call( DISORDER_MAP.out.iupred,           'IUPredscores.tsv' )
    def disorder_anchor  = _mergeDis.call( DISORDER_MAP.out.anchor,           'Anchorscores.tsv' )
    def disorder_aiupred = _mergeDis.call( DISORDER_MAP.out.aiupred,          'AIUPredscores.tsv' )
    def disorder_aiubind = _mergeDis.call( DISORDER_MAP.out.aiupred_binding,  'AIUPredBinding.tsv' )
    def disorder_plddt   = _mergeDis.call( DISORDER_MAP.out.plddt,            'AlphaFoldTable.tsv' )
    def disorder_regions = _mergeDis.call( DISORDER_MAP.out.disorder_regions, 'CombinedDisorderNew.tsv' )
    def disorder_pos     = _mergeDis.call( DISORDER_MAP.out.disorder_pos,     'CombinedDisorderNew_Pos.tsv' )

    disorder_regions.view { f ->
        "\n✔  Combined disorder regions: ${f}\n"
    }

    // ── Step 5c: PDB structure mapping ───────────────────────────────────────
    // Two routes: PDB_MAP (PDBe API per protein, also yields pdb_missing) or, when
    // params.pdb_bulk, PDB_BULK_MAP (SIFTS bulk join — ~1000× faster, no
    // pdb_missing). SIFTS comes from a local --sifts_tsv or FETCH_SIFTS.
    if ( !params.skip_pdb ) {
        if ( params.pdb_bulk ) {
            def sifts_ch = params.sifts_tsv
                ? Channel.value( file(params.sifts_tsv, checkIfExists: false) )
                : FETCH_SIFTS().sifts
            PDB_BULK_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq, sifts_ch )
            PDB_BULK_MAP.out.structures.view { f ->
                "\n✔  PDB structures (SIFTS bulk): ${f}\n"
            }
        } else {
            PDB_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq )
            PDB_MAP.out.structures.view { f ->
                "\n✔  PDB structures (transcript region + chain): ${f}\n"
            }
            PDB_MAP.out.pdb_missing.view { f ->
                "\n✔  PDB missing residues (structure-derived disorder): ${f}\n"
            }
        }
    }

    // ── Step 5d: Exon boundaries ──────────────────────────────────────────────
    if ( genome_enabled ) {
        EXON_MAP(
            GENOME_MAP.out.map_file,
            SEQUENCE_PROCESS.out.loc_chrom_seq
        )
        EXON_MAP.out.exon.view { f ->
            "\n✔  Exon boundaries: ${f}\n"
        }
    }

    // ── Step 5f: GO Term annotation ───────────────────────────────────────────
    go_refs = FETCH_GO()
    GO_MAP(
        SEQUENCE_PROCESS.out.loc_chrom_seq,
        go_refs.goa,
        go_refs.obo
    )
    GO_MAP.out.go_terms.view { f ->
        "\n✔  GO terms: ${f}\n"
    }

    // ── Step 5g: Polymorphism — ALL polymorphisms + allele frequency ──────────
    // common_poly.out provides allele frequency / rsid (mapped via combined_map);
    // polymorphism_pos.tsv provides the comprehensive all-polymorphism set.
    // (Supersedes the former separate SNP_MAP / snp_polymorphisms.tsv output.)
    if ( !params.skip_polymorphism ) {
        def snp_common_f = params.snp_out_common ? file(params.snp_out_common, checkIfExists: false) : no_file
        def snp_all_f    = params.snp_out_all    ? file(params.snp_out_all,    checkIfExists: false) : no_file
        def snp_pos_f    = params.snp_pos_tsv    ? file(params.snp_pos_tsv,     checkIfExists: false) : no_file
        def dbsnp_bb_f   = params.dbsnp_bb       ? file(params.dbsnp_bb,        checkIfExists: false)
                           : ( params.fetch_dbsnp ? FETCH_DBSNP_BB().bb : no_file )
        POLYMORPHISM_MAP(
            SEQUENCE_PROCESS.out.loc_chrom_seq,
            combined_map_ch,
            snp_common_f,
            snp_all_f,
            snp_pos_f,
            dbsnp_bb_f,
            setup_done_ch.first()
        )
        POLYMORPHISM_MAP.out.polymorphism.view { f ->
            "\n✔  Polymorphism (all + allele frequency): ${f}\n"
        }
    }

    // ── Step 5h: PEM Core Motifs ──────────────────────────────────────────────
    pem_file_ch = params.pem_dataset
        ? Channel.value( file(params.pem_dataset, checkIfExists: false) )
        : Channel.value( file("${projectDir}/assets/NO_FILE") )
    PEM_MAP(
        SEQUENCE_PROCESS.out.loc_chrom_seq,
        pem_file_ch
    )
    PEM_MAP.out.pem.view { f ->
        "\n✔  PEM Core Motifs: ${f}\n"
    }

    if ( !params.skip_pem && params.pem_transfer ) {
        PEM_TRANSFER_MAP(
            SEQUENCE_PROCESS.out.loc_chrom_seq,
            PEM_MAP.out.pem
        )
        PEM_TRANSFER_MAP.out.pem_mapped.view { f ->
            "\n✔  PEM motifs (isoform transfer): ${f}\n"
        }
    }

    // ── Step 5i: CoiledCoils (DeepCoil) ──────────────────────────────────────
    // Scattered alongside disorder (reuses the same gene-balanced chunk list).
    def coiled_coils_ch = null
    if ( !params.skip_coiledcoils ) {
        def cc_loc_ch = ( scatter_n > 1 ) ? split_chunks.flatten()
                                          : SEQUENCE_PROCESS.out.loc_chrom_seq
        COILEDCOILS_MAP( cc_loc_ch, setup_done_ch.first(), setup_deepcoil_py_ch.first() )
        def _cc_pub = params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/annotations"
                                      : "${params.outdir}/final/annotations"
        coiled_coils_ch = ( scatter_n > 1
            ? COILEDCOILS_MAP.out.coiled_coils.collectFile(
                  name: 'coiled_coils.tsv', keepHeader: true, skip: 1, sort: true, storeDir: _cc_pub)
            : COILEDCOILS_MAP.out.coiled_coils ).first()
        coiled_coils_ch.view { f ->
            "\n✔  CoiledCoils: ${f}\n"
        }
    }

    // ── Step 5j: PPI (BioGRID + IntAct + HIPPIE) ─────────────────────────────
    // Pre-processed files (params.ppi_intact etc.) take priority.
    // When not set: download raw MiTab via FETCH_* and preprocess via PPI_PREPROCESS.
    if ( !params.skip_ppi ) {
        def intact_f  = params.ppi_intact  ? file(params.ppi_intact,  checkIfExists: false) : null
        def biogrid_f = params.ppi_biogrid ? file(params.ppi_biogrid, checkIfExists: false) : null
        def hippie_f  = params.ppi_hippie  ? file(params.ppi_hippie,  checkIfExists: false) : null

        def intact_ch
        def biogrid_ch
        def hippie_ch

        if ( intact_f && biogrid_f && hippie_f ) {
            intact_ch  = Channel.value(intact_f)
            biogrid_ch = Channel.value(biogrid_f)
            hippie_ch  = Channel.value(hippie_f)
        } else {
            // At least one pre-processed file is missing → download + preprocess
            log.info "PPI: pre-processed files not fully set → downloading raw PPI databases"
            def raw_intact  = intact_f  ? Channel.value(intact_f)  : FETCH_INTACT().zip
            def raw_biogrid = biogrid_f ? Channel.value(biogrid_f) : FETCH_BIOGRID().zip
            def raw_hippie  = hippie_f  ? Channel.value(hippie_f)  : FETCH_HIPPIE().txt
            def prep = PPI_PREPROCESS(raw_intact, raw_biogrid, raw_hippie)
            intact_ch  = prep.intact
            biogrid_ch = prep.biogrid
            hippie_ch  = prep.hippie
        }

        PPI_MAP(
            SEQUENCE_PROCESS.out.loc_chrom_seq,
            intact_ch,
            biogrid_ch,
            hippie_ch
        )
        PPI_MAP.out.interactions.view { f ->
            "\n✔  PPI interactions: ${f}\n"
        }
    }

    // ── Step 7: Conservation (GOPHER multi-level + phastCons) ────────────────
    // phastCons requires combined_map.map (GENOME_MAP); GOPHER runs without it.
    if ( !params.skip_conservation && genome_enabled ) {
        // GOPHER conservation: recompute from orthologue alignments (run_gopher)
        // or consume a precomputed table (gopher_conservation_table).
        def cons_table
        if ( params.run_gopher ) {
            def aln_d = params.gopher_aln_dir
                ? file(params.gopher_aln_dir, checkIfExists: false)
                : file("${projectDir}/assets/NO_FILE")
            def taxon_f = params.gopher_taxon_map
                ? file(params.gopher_taxon_map, checkIfExists: false)
                : file("${projectDir}/assets/NO_FILE")
            cons_table = GOPHER_RECOMPUTE(
                SEQUENCE_PROCESS.out.loc_chrom_seq, aln_d, taxon_f
            ).table
        } else {
            cons_table = params.gopher_conservation_table
                ? file(params.gopher_conservation_table, checkIfExists: false)
                : file("${projectDir}/assets/NO_FILE")
        }
        CONSERVATION_MAP(
            SEQUENCE_PROCESS.out.loc_chrom_seq,
            GENOME_MAP.out.map_file,
            cons_table
        )
        CONSERVATION_MAP.out.gopher.view { f ->
            "\n✔  Conservation (GOPHER): ${f}\n"
        }
        CONSERVATION_MAP.out.phastcons.view { f ->
            "\n✔  Conservation (phastCons): ${f}\n"
        }
    }

    // ── Step 5k: ScanSite 4.0 kinase motif predictions ───────────────────────
    if ( !params.skip_scansite ) {
        def scansite_f = params.scansite_tsv
            ? file(params.scansite_tsv, checkIfExists: false)
            : no_file
        SCANSITE_MAP(SEQUENCE_PROCESS.out.loc_chrom_seq, scansite_f)
        SCANSITE_MAP.out.scansite.view { f ->
            "\n✔  ScanSite motifs: ${f}\n"
        }
    }

    // ── Step 5l: (removed) Population SNP polymorphisms ──────────────────────
    // snp_polymorphisms.tsv is superseded by POLYMORPHISM_MAP, which now emits
    // ALL polymorphisms (from polymorphism_pos.tsv) enriched with allele
    // frequency (from common_poly.out) in a single polymorphism.tsv.

    // ── Step 8a: ClinVar Disease Ontology (filter fallback) ─────────────────
    // Skipped when CLINVAR_DISEASE_BUILD already ran (hg38 + clinvar_disease_from_mutations).
    if ( !params.skip_clinvar_disease
         && !(genome_enabled && params.clinvar_disease_from_mutations) ) {
        def clinvar_disease_f = params.clinvar_disease_tsv
            ? file(params.clinvar_disease_tsv, checkIfExists: false)
            : no_file
        def clinvar_cat_f = params.clinvar_category_tsv
            ? file(params.clinvar_category_tsv, checkIfExists: false)
            : no_file
        CLINVAR_DISEASE_MAP(SEQUENCE_PROCESS.out.loc_chrom_seq, clinvar_disease_f, clinvar_cat_f)
        CLINVAR_DISEASE_MAP.out.clinvar_disease.view { f ->
            "\n✔  ClinVar disease (filter fallback): ${f}\n"
        }
    }

    // ── Step 8b: OMIM Disease Ontology ───────────────────────────────────────
    // Raw path: FETCH_OMIM (key-gated download of genemap2.txt) → raw parse.
    // Processed path: pre-built disease/variant tables (--omim_tsv).
    if ( !params.skip_omim ) {
        def omim_f = params.fetch_omim
            ? FETCH_OMIM().dir
            : ( params.omim_raw_dir ? file(params.omim_raw_dir, checkIfExists: false)
              : ( params.omim_tsv   ? file(params.omim_tsv, checkIfExists: false) : no_file ) )
        def omim_mut_f = params.omim_mutations_tsv
            ? file(params.omim_mutations_tsv, checkIfExists: false)
            : no_file
        OMIM_MAP(SEQUENCE_PROCESS.out.loc_chrom_seq, omim_f, omim_mut_f)
        OMIM_MAP.out.omim_disease.view { f ->
            "\n✔  OMIM disease: ${f}\n"
        }
        OMIM_MAP.out.omim_mutations.view { f ->
            "\n✔  OMIM mutations: ${f}\n"
        }
    }

    // ── Step 8c: Cancer Gene Census + Compendium ─────────────────────────────
    if ( !params.skip_cancer_drivers ) {
        def cancer_driver_f = params.cancer_driver_tsv
            ? file(params.cancer_driver_tsv, checkIfExists: false)
            : no_file
        def census_roles_f = params.census_roles_tsv
            ? file(params.census_roles_tsv, checkIfExists: false)
            : no_file
        def compendium_roles_f = params.compendium_roles_tsv
            ? file(params.compendium_roles_tsv, checkIfExists: false)
            : no_file
        CANCER_DRIVER_MAP(SEQUENCE_PROCESS.out.loc_chrom_seq, cancer_driver_f,
                          census_roles_f, compendium_roles_f)
        CANCER_DRIVER_MAP.out.combined.view { f ->
            "\n✔  Cancer drivers (Census/Compendium): ${f}\n"
        }
    }

    // ── Step 8d: AlphaMissense GENCODE isoform scores ─────────────────────────
    // Decompress the .gz once (storeDir-cached); ALPHAMISSENSE_MAP reads the
    // plain TSV with the pandas C engine (much faster than gzip line-by-line).
    if ( !params.skip_alphamissense ) {
        def am_gz_ch = params.alphamissense_gz
            ? Channel.value( file(params.alphamissense_gz, checkIfExists: false) )
            : FETCH_ALPHAMISSENSE().gz
        def am_ch = DECOMPRESS_ALPHAMISSENSE(am_gz_ch).tsv
        ALPHAMISSENSE_MAP(SEQUENCE_PROCESS.out.loc_chrom_seq, am_ch)
        ALPHAMISSENSE_MAP.out.alphamissense.view { f ->
            "\n✔  AlphaMissense isoforms: ${f}\n"
        }
    }

    // ── Step 8e: DepMap cancer cell line mutations ────────────────────────────
    // Default source: local manual --depmap_tsv. FETCH_DEPMAP remains available
    // only when explicitly requested, because the DepMap portal may require
    // browser verification and block scripted catalogue access.
    if ( !params.skip_depmap ) {
        def depmap_src = params.fetch_depmap
            ? FETCH_DEPMAP().tsv
            : ( depmap_manual_ch ?: Channel.value(no_file) )
        DEPMAP_MAP(SEQUENCE_PROCESS.out.loc_chrom_seq, depmap_src)
        DEPMAP_MAP.out.depmap.view { f ->
            "\n✔  DepMap mutations: ${f}\n"
        }
    }

    // ── Step 8f: dbNSFP pathogenicity scores (pre-mapped fallback when no raw dir) ─
    if ( !params.skip_pathogenicity && !(genome_enabled && (params.dbnsfp_raw_dir || params.fetch_dbnsfp)) ) {
        def dbnsfp_f = params.dbnsfp_tsv
            ? file(params.dbnsfp_tsv, checkIfExists: false)
            : no_file
        PATHOGENICITY_MAP(SEQUENCE_PROCESS.out.loc_chrom_seq, dbnsfp_f)
        PATHOGENICITY_MAP.out.scores.view { f ->
            "\n✔  Pathogenicity scores (dbNSFP filter): ${f}\n"
        }
    }

    // ── Step 8h: FINCHES saturation mutagenesis (Δε LLPS-change score) ───────
    // Computes homotypic self-interaction energy change for every possible
    // single-AA substitution using the Mpipi_GGv1 force field.
    // ⚠  Non-commercial only (FINCHES CC BY-NC 4.0).
    if ( !params.skip_finches ) {
        FINCHES_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq )
        FINCHES_MAP.out.finches.view { f ->
            "\n✔  FINCHES Δε saturation scan: ${f}\n"
        }
    }

    // ── Step 5m: PositionBasedAnnotations + RSAscores ────────────────────────
    // Aggregates per-residue data: IUPred, pLDDT, CombinedDisorder, phastCons,
    // GOPHER conservation, and Pfam domain coverage into one row-per-position
    // table. RSA is derived from pLDDT as (100 - plddt) / 100.
    def pb_cons_gopher_ch = (!params.skip_conservation && genome_enabled)
        ? CONSERVATION_MAP.out.gopher
        : Channel.value(no_file)
    def pb_cons_pcons_ch  = (!params.skip_conservation && genome_enabled)
        ? CONSERVATION_MAP.out.phastcons
        : Channel.value(no_file)
    POSITION_BASED_MAP(
        SEQUENCE_PROCESS.out.loc_chrom_seq,
        disorder_iupred,
        disorder_plddt,
        disorder_pos,
        pb_cons_pcons_ch,
        pb_cons_gopher_ch,
        ANNOTATION_MAP.out.pfam,
    )
    POSITION_BASED_MAP.out.pos_annotations.view { f ->
        "\n✔  Position-based annotations: ${f}\n"
    }
    POSITION_BASED_MAP.out.rsa_scores.view { f ->
        "\n✔  RSA scores: ${f}\n"
    }

    // ── Step 5e: Transcript Mapping ───────────────────────────────────────────
    // Map all UniProt-keyed annotations onto each Gencode transcript.
    // Annotations that match 100% within a transcript → direct map.
    // Annotations that match 100% in a different isoform → homology_transfer.
    TRANSCRIPT_MAP(
        SEQUENCE_PROCESS.out.loc_chrom_seq,
        ANNOTATION_MAP.out.elm_mapped,
        ANNOTATION_MAP.out.dibs,
        ANNOTATION_MAP.out.mfib,
        ANNOTATION_MAP.out.phasepro,
        ANNOTATION_MAP.out.uniprot_roi,
        ANNOTATION_MAP.out.uniprot_binding,
        ANNOTATION_MAP.out.ptm,
        ANNOTATION_MAP.out.pfam,
        disorder_regions,   // CombinedDisorderNew.tsv (pass-through)
        disorder_pos        // CombinedDisorderNew_Pos.tsv (pass-through)
    )

    TRANSCRIPT_MAP.out.stats.view { f ->
        "\n✔  Transcript mapping stats: ${f}\n"
    }

    // ── Homology-similarity manifest (audit of main→alt isoform transfers) ───
    homology_manifest_inputs = TRANSCRIPT_MAP.out.elm
        .mix( TRANSCRIPT_MAP.out.dibs,
              TRANSCRIPT_MAP.out.mfib,
              TRANSCRIPT_MAP.out.phasepro,
              TRANSCRIPT_MAP.out.roi,
              TRANSCRIPT_MAP.out.bind,
              TRANSCRIPT_MAP.out.ptm,
              TRANSCRIPT_MAP.out.pfam,
              ELM_SWITCHES_MAP.out.elm_switches )
        .collect()
    HOMOLOGY_MANIFEST( homology_manifest_inputs )
    HOMOLOGY_MANIFEST.out.manifest.view { f ->
        "\n✔  Homology-similarity manifest: ${f}\n"
    }

    // ── Comprehensive per-protein mapping report + reproducible summary ──────
    // Reads the published final/ + intermediate/ directories, so it must run
    // LAST. We gate it on a broad bundle of the slowest / terminal outputs
    // (guarded by the same skip conditions that produce them) so every
    // publishDir copy has completed before the report scans the directory.
    report_gate = HOMOLOGY_MANIFEST.out.manifest
        .mix( TRANSCRIPT_MAP.out.stats,
              disorder_pos,
              POSITION_BASED_MAP.out.rsa_scores,
              GO_MAP.out.go_terms )
    // When scattering, the merge (collectFile) of the disorder tables that have no
    // downstream consumer must still be materialised so they publish; gating the
    // report on them also guarantees publish completes before the report scans.
    if ( scatter_n > 1 ) report_gate = report_gate.mix( disorder_anchor, disorder_aiupred, disorder_aiubind )
    if ( genome_enabled )                report_gate = report_gate.mix( GENOME_MAP.out.map_file, GENOME_QUERY_MAP.out.mutations )
    if ( !params.skip_pdb )                report_gate = report_gate.mix( params.pdb_bulk ? PDB_BULK_MAP.out.structures : PDB_MAP.out.structures )
    if ( !params.skip_ppi )                report_gate = report_gate.mix( PPI_MAP.out.interactions )
    if ( !params.skip_coiledcoils )        report_gate = report_gate.mix( coiled_coils_ch )
    if ( !params.skip_scansite )           report_gate = report_gate.mix( SCANSITE_MAP.out.scansite )
    if ( !params.skip_polymorphism )       report_gate = report_gate.mix( POLYMORPHISM_MAP.out.polymorphism )
    if ( !params.skip_conservation && genome_enabled ) report_gate = report_gate.mix( CONSERVATION_MAP.out.phastcons )
    if ( !params.skip_alphamissense )      report_gate = report_gate.mix( ALPHAMISSENSE_MAP.out.alphamissense )
    if ( mavedb_enabled )                  report_gate = report_gate.mix( MAVEDB_MAP.out.mavedb )
    if ( proteingym_enabled )              report_gate = report_gate.mix( PROTEINGYM_MAP.out.proteingym )
    if ( !params.skip_depmap )             report_gate = report_gate.mix( DEPMAP_MAP.out.depmap )
    if ( !params.skip_cancer_drivers )     report_gate = report_gate.mix( CANCER_DRIVER_MAP.out.combined )
    if ( !params.skip_omim )               report_gate = report_gate.mix( OMIM_MAP.out.omim_disease )
    if ( !params.skip_pathogenicity && !(genome_enabled && (params.dbnsfp_raw_dir || params.fetch_dbnsfp)) ) report_gate = report_gate.mix( PATHOGENICITY_MAP.out.scores )
    if ( !params.skip_finches )              report_gate = report_gate.mix( FINCHES_MAP.out.finches )

    MAPPING_REPORT(
        SEQUENCE_PROCESS.out.loc_chrom_seq,
        report_gate.collect()
    )
    MAPPING_REPORT.out.reports.view { f ->
        "\n✔  Mapping reports: ${f}\n"
    }

}
