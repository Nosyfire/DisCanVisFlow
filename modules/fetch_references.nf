/*
 * modules/fetch_references.nf
 *
 * Downloads all reference FASTA and annotation files from public FTP servers.
 *
 * storeDir semantics
 * ──────────────────
 * Output files are stored permanently in params.ref_dir (outside Nextflow's
 * work/ directory). If the expected output file already exists there, the
 * process is skipped entirely — no re-download on subsequent runs.
 * Delete the file from ref_dir to force a fresh download.
 *
 * Processes
 * ─────────
 *   FETCH_UNIPROT_FASTA        → uniprot_swissprot.fasta   (SwissProt >sp| only)
 *   FETCH_GENCODE_TRANSLATIONS → gencode_pc_translations.fasta
 *   FETCH_GENCODE_TRANSCRIPTS  → gencode_pc_transcripts.fasta  (Module 2)
 *   FETCH_GENCODE_GTF          → gencode_annotation.gtf.gz     (Module 2)
 *   FETCH_GENCODE_REFSEQ       → gencode_metadata_refseq.tsv   (Module 2)
 *   FETCH_MONDO                → mondo.obo   (Module 8a — ClinVar disease build)
 *   FETCH_ALPHAMISSENSE        → AlphaMissense_isoforms_hg38.tsv.gz  (Module 8d)
 *   FETCH_INTACT               → intact_human.mitab.zip  (Module 5j raw download)
 *   FETCH_BIOGRID              → biogrid_human.mitab.zip (Module 5j raw download)
 *   FETCH_HIPPIE               → hippie_current.txt      (Module 5j raw download)
 */


// ── UniProt SwissProt-only reference proteome ─────────────────────────────────
//
// Source: UP000005640 human reference proteome (contains both sp| and tr|).
// We stream-filter with awk to keep only reviewed SwissProt entries (>sp|).
// TrEMBL (>tr|) entries are discarded — they are not needed for BLAST ID map.
//
process FETCH_UNIPROT_FASTA {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/uniprot" : "${params.ref_dir}/uniprot" }

    output:
    path 'uniprot_swissprot.fasta', emit: fasta

    script:
    """
    echo "Downloading UniProt reference proteome (SwissProt filter)..."
    echo "  Source: ${params.uniprot_ftp_url}"

    curl -fsSL '${params.uniprot_ftp_url}' \\
        | zcat \\
        | awk '/^>/{keep=/^>sp\\|/} keep' \\
        > uniprot_swissprot.fasta

    n=\$(grep -c '^>' uniprot_swissprot.fasta)
    echo "Done — \${n} SwissProt entries written to uniprot_swissprot.fasta"
    """

    stub:
    """
    printf '>sp|P04049|RAF1_HUMAN Proto-oncogene serine/threonine-protein kinase Raf OS=Homo sapiens OX=9606 GN=RAF1 PE=1 SV=1\\nMANTIQQFLHR\\n' \\
        > uniprot_swissprot.fasta
    """
}


// ── GENCODE protein translations (all protein-coding isoforms) ────────────────
//
// Used in Module 0+1 (BLAST query / database).
// Contains all annotated protein-coding transcript translations (isoforms
// included) — no filtering; SUBSET_FASTA handles gene-level subsetting.
//
process FETCH_GENCODE_TRANSLATIONS {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/gencode" : "${params.ref_dir}/gencode" }

    output:
    path 'gencode_pc_translations.fasta', emit: fasta

    script:
    """
    echo "Downloading GENCODE protein translations..."
    echo "  Source: ${params.gencode_translations_url}"

    curl -fsSL '${params.gencode_translations_url}' \\
        | zcat > gencode_pc_translations.fasta

    n=\$(grep -c '^>' gencode_pc_translations.fasta)
    echo "Done — \${n} GENCODE translation entries written."
    """

    stub:
    """
    printf '>ENST00000423430.6|ENSG00000132155.11|OTTHUMT00000076045.3|ENST00000423430.6|RAF1-201|RAF1|648\\nMANTIQQFLHR\\n' \\
        > gencode_pc_translations.fasta
    """
}


// ── GENCODE protein-coding cDNA transcripts (Module 2 – Sequence Process) ─────
//
// Used in Module 2 to obtain cDNA sequences for genome mapping (BLAT).
//
process FETCH_GENCODE_TRANSCRIPTS {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/gencode" : "${params.ref_dir}/gencode" }

    output:
    path 'gencode_pc_transcripts.fasta', emit: fasta

    script:
    """
    echo "Downloading GENCODE cDNA transcripts..."
    echo "  Source: ${params.gencode_transcripts_url}"

    curl -fsSL '${params.gencode_transcripts_url}' \\
        | zcat > gencode_pc_transcripts.fasta

    n=\$(grep -c '^>' gencode_pc_transcripts.fasta)
    echo "Done — \${n} GENCODE transcript entries written."
    """

    stub:
    """
    printf '>ENST00000423430.6|ENSG00000132155.11|OTTHUMT00000076045.3|ENST00000423430.6|RAF1-201|RAF1|2583\\nATGGCGAATACGATG\\n' \\
        > gencode_pc_transcripts.fasta
    """
}


// ── GENCODE comprehensive annotation GTF (Module 2 – genomic coordinates) ─────
//
// Used in Module 2 to extract chromosome, strand, exon/CDS coordinates,
// and isoform flags (Ensembl_canonical, MANE_Select, appris_principal_1).
// Kept gzip-compressed; downstream Python workers handle gz input natively.
//
process FETCH_GENCODE_GTF {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/gencode" : "${params.ref_dir}/gencode" }

    output:
    path 'gencode_annotation.gtf.gz', emit: gtf

    script:
    """
    echo "Downloading GENCODE annotation GTF..."
    echo "  Source: ${params.gencode_gtf_url}"

    curl -fsSL '${params.gencode_gtf_url}' -o gencode_annotation.gtf.gz

    sz=\$(du -sh gencode_annotation.gtf.gz | cut -f1)
    echo "Done — GTF size: \${sz}"
    """

    stub:
    """
    printf '##description: stub GENCODE GTF\\nchr1\\tHAVANA\\tgene\\t1\\t1000\\t.\\t+\\t.\\tgene_id "ENSG00000132155.11"; gene_name "RAF1";\\n' \\
        | gzip > gencode_annotation.gtf.gz
    """
}


// ── GENCODE RefSeq cross-reference metadata (Module 2 – RefSeq IDs) ──────────
//
// Tab-separated: transcript_id <TAB> RefSeq_mRNA_id <TAB> RefSeq_protein_id
// Used in Module 2 to annotate transcripts with their RefSeq accessions.
//
process FETCH_GENCODE_REFSEQ {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/gencode" : "${params.ref_dir}/gencode" }

    output:
    path 'gencode_metadata_refseq.tsv', emit: tsv

    script:
    """
    echo "Downloading GENCODE RefSeq metadata..."
    echo "  Source: ${params.gencode_refseq_url}"

    curl -fsSL '${params.gencode_refseq_url}' \\
        | zcat > gencode_metadata_refseq.tsv

    n=\$(wc -l < gencode_metadata_refseq.tsv)
    echo "Done — \${n} RefSeq metadata rows written."
    """

    stub:
    """
    printf 'ENST00000423430.6\\tNM_002880.4\\tNP_002871.1\\n' > gencode_metadata_refseq.tsv
    """
}


// ── MONDO OBO ontology (Module 8a — ClinVar disease build) ───────────────────
//
// MONDO is used by create_clinvar_disease_build_worker.py to map ClinVar
// disease names → structured DO ontology hierarchy.
// Fallback: supply --mondo_obo <path> in nextflow.config or CLI.
//
process FETCH_MONDO {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/mondo" : "${params.ref_dir}/mondo" }

    output:
    path 'mondo.obo', emit: obo

    script:
    """
    echo "Downloading MONDO OBO ontology..."
    curl -fsSL \\
        'https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo' \\
        -o mondo.obo

    n=\$(grep -c '^\\[Term\\]' mondo.obo)
    echo "Done — \${n} MONDO terms written."
    """

    stub:
    """
    printf '[Term]\\nid: MONDO:0000001\\nname: disease or disorder\\n' > mondo.obo
    """
}


// ── AlphaMissense GENCODE isoform scores (Module 8d) ─────────────────────────
//
// Raw AlphaMissense per-residue pathogenicity scores for all GENCODE transcripts
// (hg38, ENST-keyed). ~1.1 GB gzipped; cached in ref_dir.
// Fallback: supply --alphamissense_gz <path> in nextflow.config or CLI.
//
process FETCH_ALPHAMISSENSE {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/alphamissense" : "${params.ref_dir}/alphamissense" }

    output:
    path 'AlphaMissense_isoforms_hg38.tsv.gz', emit: gz

    script:
    """
    echo "Downloading AlphaMissense isoforms (hg38)..."
    curl -fsSL \\
        'https://storage.googleapis.com/dm_alphamissense/AlphaMissense_isoforms_hg38.tsv.gz' \\
        -o AlphaMissense_isoforms_hg38.tsv.gz

    sz=\$(du -sh AlphaMissense_isoforms_hg38.tsv.gz | cut -f1)
    echo "Done — \${sz} written."
    """

    stub:
    """
    printf '#AlphaMissense stub\\ntranscript_id\\tprotein_variant\\tam_pathogenicity\\tam_class\\n' \\
        | gzip > AlphaMissense_isoforms_hg38.tsv.gz
    """
}

// Decompress the AlphaMissense .gz once and cache the plain TSV via storeDir.
// ALPHAMISSENSE_MAP reads the uncompressed file with the pandas C engine — much
// faster than Python's line-by-line gzip scanning over NFS (12 min → ~1 min).
process DECOMPRESS_ALPHAMISSENSE {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/alphamissense" : "${params.ref_dir}/alphamissense" }

    input:
    path gz_file

    output:
    path 'AlphaMissense_isoforms_hg38.tsv', emit: tsv

    script:
    """
    gzip -dc ${gz_file} > AlphaMissense_isoforms_hg38.tsv
    """

    stub:
    """
    printf '#CHROM\\tPOS\\tREF\\tALT\\tgenome\\ttranscript_id\\tprotein_variant\\tam_pathogenicity\\tam_class\\n' \
        > AlphaMissense_isoforms_hg38.tsv
    """
}


// ── IntAct human PPI (Module 5j raw source) ───────────────────────────────────
//
// Human-only IntAct interaction data in MiTab2.7 format.
// Converted to Interaction_intact.tsv by PPI_PREPROCESS.
// Fallback: supply pre-processed --ppi_intact <path> directly.
//
process FETCH_INTACT {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/ppi" : "${params.ref_dir}/ppi" }

    output:
    path 'intact_human.mitab.zip', emit: zip

    script:
    """
    echo "Downloading IntAct human interactions (MiTab2.7)..."
    # EBI retired the old intact-homo_sapiens.zip path; the per-species human
    # export now lives under psimitab/species/Human.zip.
    curl -fsSL \\
        'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/Human.zip' \\
        -o intact_human.mitab.zip

    sz=\$(du -sh intact_human.mitab.zip | cut -f1)
    echo "Done — \${sz} written."
    """

    stub:
    """
    printf 'uniprotkb:P11111\\tuniprotkb:P22222\\t-\\t-\\t-\\t-\\tpsi-mi:"MI:0018"\\t-\\tpubmed:12345678\\ttaxid:9606(human)\\ttaxid:9606(human)\\t-\\t-\\t-\\tintact-miscore:0.37\\n' \\
        | zip intact_human.mitab.zip /dev/stdin 2>/dev/null || \\
    printf 'uniprotkb:P11111\\tuniprotkb:P22222\\t-\\t-\\t-\\t-\\tpsi-mi:"MI:0018"\\t-\\tpubmed:12345678\\ttaxid:9606(human)\\ttaxid:9606(human)\\t-\\t-\\t-\\tintact-miscore:0.37\\n' \\
        > intact_human.mitab.zip
    """
}


// ── BioGRID human PPI (Module 5j raw source) ─────────────────────────────────
//
// Human BioGRID interactions in MiTab format.
//
process FETCH_BIOGRID {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/ppi" : "${params.ref_dir}/ppi" }

    output:
    path 'biogrid_human.mitab.zip', emit: zip

    script:
    """
    echo "Downloading BioGRID interactions (all-organism MiTab; the per-species"
    echo "Homo_sapiens-LATEST path was retired — PPI_PREPROCESS extracts the human"
    echo "member and taxon-filters)..."
    curl -fsSL \\
        'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip' \\
        -o biogrid_human.mitab.zip

    sz=\$(du -sh biogrid_human.mitab.zip | cut -f1)
    echo "Done — \${sz} written."
    """

    stub:
    """
    touch biogrid_human.mitab.zip
    """
}


// ── HIPPIE human PPI (Module 5j raw source) ───────────────────────────────────
//
// HIPPIE human interactome (all interactions are already human).
//
process FETCH_HIPPIE {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/ppi" : "${params.ref_dir}/ppi" }

    output:
    path 'hippie_current.txt', emit: txt

    script:
    """
    echo "Downloading HIPPIE interactome..."
    curl -fsSL \\
        'https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt' \\
        -o hippie_current.txt

    n=\$(wc -l < hippie_current.txt)
    echo "Done — \${n} interactions written."
    """

    stub:
    """
    printf 'GENEA\\tP11111\\tGENEB\\tP22222\\t0.63\\texperiments:1,pmids:12345678\\n' > hippie_current.txt
    """
}


// ──────────────────────────────────────────────────────────────────────────
// Reproducibility FETCH processes (verified URLs — see REF_FETCH_NOTES.md).
// All cache via storeDir so a clean machine downloads once; on a machine that
// already has the file, pre-seed ${params.ref_dir}/<sub>/ to skip the download.
// ──────────────────────────────────────────────────────────────────────────

// UCSC hg38 2bit genome (BLAT / genome mapping) — ~797 MiB
process FETCH_HG38_2BIT {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/hg38" : "${params.ref_dir}/hg38" }
    output:
    path 'hg38.2bit', emit: twobit
    script:
    """
    echo "Downloading UCSC hg38.2bit (~797 MiB)..."
    curl -fsSL 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit' -o hg38.2bit
    echo "Done — \$(du -sh hg38.2bit | cut -f1)."
    """
    stub:
    """ touch hg38.2bit """
}

// UCSC dbSnp155 Common track (polymorphisms) — ~1.83 GiB
process FETCH_DBSNP_BB {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/dbsnp" : "${params.ref_dir}/dbsnp" }
    output:
    path 'dbSnp155Common.bb', emit: bb
    script:
    """
    echo "Downloading UCSC dbSnp155Common.bb (~1.83 GiB)..."
    curl -fsSL 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb' -o dbSnp155Common.bb
    echo "Done — \$(du -sh dbSnp155Common.bb | cut -f1)."
    """
    stub:
    """ touch dbSnp155Common.bb """
}

// UniProt human proteome additional (isoform) FASTA — sp|Pxxxxx-N| headers — ~40 MiB gz
process FETCH_UNIPROT_ISOFORMS {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/uniprot" : "${params.ref_dir}/uniprot" }
    output:
    path 'UP000005640_9606_additional.fasta', emit: fasta
    script:
    """
    echo "Downloading UniProt UP000005640 additional (isoform) FASTA..."
    curl -fsSL 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.fasta.gz' -o add.fasta.gz
    gunzip -c add.fasta.gz > UP000005640_9606_additional.fasta
    rm -f add.fasta.gz
    echo "Done — \$(grep -c '^>' UP000005640_9606_additional.fasta) isoform sequences."
    """
    stub:
    """ printf '>sp|P04637-2|P53_HUMAN Isoform 2\\nMEEPQSDPSV\\n' > UP000005640_9606_additional.fasta """
}

// SIFTS UniProt↔PDB chain mapping (bulk PDB, replaces per-protein PDBe API) — ~6 MiB gz
process FETCH_SIFTS {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/sifts" : "${params.ref_dir}/sifts" }
    output:
    path 'pdb_chain_uniprot.tsv.gz', emit: sifts
    script:
    """
    echo "Downloading SIFTS pdb_chain_uniprot.tsv.gz (~6 MiB)..."
    curl -fsSL 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz' -o pdb_chain_uniprot.tsv.gz
    echo "Done — \$(du -sh pdb_chain_uniprot.tsv.gz | cut -f1)."
    """
    stub:
    """ printf 'PDB\\tCHAIN\\tSP_PRIMARY\\tRES_BEG\\tRES_END\\tPDB_BEG\\tPDB_END\\tSP_BEG\\tSP_END\\n' | gzip > pdb_chain_uniprot.tsv.gz """
}

// AlphaFold human-proteome bulk archive (pLDDT) — replaces per-accession API — ~4.82 GiB
process FETCH_ALPHAFOLD_BULK {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/alphafold" : "${params.ref_dir}/alphafold" }
    output:
    path 'UP000005640_9606_HUMAN_v6.tar', emit: tar
    script:
    """
    echo "Downloading AlphaFold human proteome v6 (~4.82 GiB)..."
    curl -fsSL 'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar' -o UP000005640_9606_HUMAN_v6.tar
    echo "Done — \$(du -sh UP000005640_9606_HUMAN_v6.tar | cut -f1)."
    """
    stub:
    """ tar -cf UP000005640_9606_HUMAN_v6.tar --files-from /dev/null """
}


// ──────────────────────────────────────────────────────────────────────────
// FETCH_ZENODO — pull the non-refreshing curated archive (ELM/DIBS/MFIB/
// PhasePro/PTM/drivers/PEM/…) from a Zenodo record and extract it. The extracted
// tree mirrors legacy_data/ so params can point at ${ref_dir}/zenodo/<...>.
// Set params.zenodo_url to the record's file download URL. storeDir-cached.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_ZENODO {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/zenodo" : "${params.ref_dir}/zenodo" }
    output:
    path 'discanvis_legacy', type: 'dir', emit: dir
    script:
    """
    echo "Downloading DisCanVis curated archive from Zenodo..."
    curl -fsSL '${params.zenodo_url}' -o archive.tar.gz
    mkdir -p discanvis_legacy
    tar -xzf archive.tar.gz -C discanvis_legacy --strip-components=0
    rm -f archive.tar.gz
    echo "Extracted: \$(find discanvis_legacy -type f | wc -l) files."
    """
    stub:
    """ mkdir -p discanvis_legacy """
}


// ──────────────────────────────────────────────────────────────────────────
// FETCH_MAVEDB — MaveDB single-mutant functional scores from the official
// api.mavedb.org. Emits a UniProt-keyed raw table (mavedb_raw.tsv) that
// MAVEDB_MAP fans out onto the run's isoforms (--mapping_mode uniprot).
// storeDir-cached so the (slow) API enumeration runs once.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_MAVEDB {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/mavedb" : "${params.ref_dir}/mavedb" }
    output:
    path 'mavedb_raw.tsv', emit: raw
    script:
    def lim = params.mavedb_max_sets ? "--max_sets ${params.mavedb_max_sets}" : ''
    """
    echo "Fetching MaveDB score sets from api.mavedb.org..."
    fetch_mavedb_worker.py --outdir . --cache_dir . ${lim}
    test -s mavedb_raw.tsv || printf 'uniprot\\tgene_name\\turn\\tmavedb_id\\tprot_expr\\tprotein_start\\tscore\\tis_double_mutant\\n' > mavedb_raw.tsv
    echo "Done — \$(wc -l < mavedb_raw.tsv) lines."
    """
    stub:
    """ printf 'uniprot\\tgene_name\\turn\\tmavedb_id\\tprot_expr\\tprotein_start\\tscore\\tis_double_mutant\\n' > mavedb_raw.tsv """
}

// ──────────────────────────────────────────────────────────────────────────
// FETCH_PROTEINGYM — ProteinGym substitution DMS assays. Emits a UniProt-keyed
// raw table (proteingym_raw.tsv) that PROTEINGYM_MAP fans out onto the run's
// isoforms (--mapping_mode uniprot). storeDir-cached.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_PROTEINGYM {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/proteingym" : "${params.ref_dir}/proteingym" }
    output:
    path 'proteingym_raw.tsv', emit: raw
    script:
    def url = params.proteingym_zip_url ? "--zip_url '${params.proteingym_zip_url}'" : ''
    def lim = params.proteingym_max_assays ? "--max_assays ${params.proteingym_max_assays}" : ''
    """
    echo "Fetching ProteinGym substitution DMS assays..."
    fetch_proteingym_worker.py --out proteingym_raw.tsv --cache_dir . ${url} ${lim}
    test -s proteingym_raw.tsv || printf 'uniprot\\tgene_name\\tDMS_id\\tprotein_variant\\tpos\\tDMS_score\\tDMS_score_bin\\n' > proteingym_raw.tsv
    echo "Done — \$(wc -l < proteingym_raw.tsv) lines."
    """
    stub:
    """ printf 'uniprot\\tgene_name\\tDMS_id\\tprotein_variant\\tpos\\tDMS_score\\tDMS_score_bin\\n' > proteingym_raw.tsv """
}


// ──────────────────────────────────────────────────────────────────────────
// FETCH_DEPMAP — DepMap cancer cell-line somatic mutations (OPEN data).
// Resolves the newest OmicsSomaticMutations.csv via the DepMap download
// catalogue (presigned URL) and normalises it to the columns DEPMAP_MAP wants.
// Override the release with --depmap_release or a direct file with --depmap_url.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_DEPMAP {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/depmap" : "${params.ref_dir}/depmap" }
    output:
    path 'depmap_mutations_raw.tsv', emit: tsv
    script:
    def rel = params.depmap_release ? "--release '${params.depmap_release}'" : ''
    def url = params.depmap_url     ? "--file_url '${params.depmap_url}'"     : ''
    """
    echo "Fetching DepMap somatic mutations (open data)..."
    fetch_depmap_worker.py --out depmap_mutations_raw.tsv --cache_dir . ${rel} ${url}
    echo "Done — \$(wc -l < depmap_mutations_raw.tsv) lines."
    """
    stub:
    """ printf 'HugoSymbol\\tProtein_position\\tHGVSp_Short\\tModelID\\tStart_Position\\tEntrezGeneID\\tHotspot\\n' > depmap_mutations_raw.tsv """
}

// ──────────────────────────────────────────────────────────────────────────
// FETCH_DBNSFP — raw dbNSFP academic distribution (per-chromosome variant gz).
// dbNSFP is free for academic use; supply the release zip URL via --dbnsfp_url
// (e.g. the dbNSFP S3/Drive link or a Zenodo mirror). Extracts the chr*.gz into
// a directory that DBNSFP_MAP consumes through --dbnsfp_raw_dir.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_DBNSFP {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/dbnsfp" : "${params.ref_dir}/dbnsfp" }
    output:
    path 'dbnsfp_raw', type: 'dir', emit: dir
    script:
    if( !params.dbnsfp_url )
        error "FETCH_DBNSFP needs --dbnsfp_url (academic dbNSFP zip / Zenodo mirror URL)."
    """
    echo "Downloading raw dbNSFP distribution..."
    curl -fsSL '${params.dbnsfp_url}' -o dbnsfp.zip
    mkdir -p dbnsfp_raw
    # dbNSFP ships dbNSFPx.y_variant.chr*.gz (+ tbi/readme); flatten into dbnsfp_raw/
    unzip -o -j dbnsfp.zip '*variant.chr*' -d dbnsfp_raw
    rm -f dbnsfp.zip
    echo "Extracted: \$(ls dbnsfp_raw | wc -l) files."
    """
    stub:
    """ mkdir -p dbnsfp_raw; printf '' | gzip > dbnsfp_raw/dbNSFP_variant.chr1.gz """
}

// ──────────────────────────────────────────────────────────────────────────
// FETCH_OMIM — raw OMIM download files (key-gated). OMIM grants per-user
// download keys; the files live at https://data.omim.org/downloads/<KEY>/.
// Supply the key via --omim_download_key (or a full mirror via --omim_url).
// Emits the raw genemap2.txt + morbidmap.txt + mimTitles.txt that a raw OMIM
// parse consumes (create_omim_worker.py --mapping_mode raw).
// ──────────────────────────────────────────────────────────────────────────
// humsavar.txt (UniProt, OPEN) is the primary protein-level disease-variant
// source (the data the user mapped from). genemap2.txt (OMIM key-gated) is an
// optional gene→phenotype disease-ontology enrichment.
process FETCH_OMIM {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/omim" : "${params.ref_dir}/omim" }
    output:
    path 'omim_raw', type: 'dir', emit: dir
    script:
    def humsavar_url = params.humsavar_url ?:
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/humsavar.txt"
    def base = params.omim_url ?: (params.omim_download_key ?
        "https://data.omim.org/downloads/${params.omim_download_key}" : "")
    """
    echo "Downloading UniProt humsavar.txt (open)..."
    mkdir -p omim_raw
    curl -fsSL '${humsavar_url}' -o omim_raw/humsavar.txt
    if [ -n "${base}" ]; then
        echo "Downloading OMIM genemap2 (key-gated)..."
        for f in genemap2.txt morbidmap.txt mimTitles.txt; do
            curl -fsSL "${base}/\$f" -o "omim_raw/\$f" || echo "[WARN] \$f not available"
        done
    fi
    echo "OMIM files: \$(ls omim_raw | wc -l) (\$(wc -l < omim_raw/humsavar.txt) humsavar lines)."
    """
    stub:
    """ mkdir -p omim_raw; printf '_________\\nKRAS P01116 VAR_1 p.Gly12Asp LP/P rs1 Noonan syndrome\\n' > omim_raw/humsavar.txt """
}

// ──────────────────────────────────────────────────────────────────────────
// FETCH_CBIOPORTAL — public cBioPortal study bundle from the datahub. Supply
// the study id via --cbioportal_study (e.g. 'msk_impact_2017') or a full URL
// via --cbioportal_url. Extracts data_mutations.txt (MAF) for MUTATION_MAP.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_CBIOPORTAL {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/cbioportal" : "${params.ref_dir}/cbioportal" }
    output:
    path 'cbioportal_mutations.maf', emit: maf
    script:
    if( !params.cbioportal_study && !params.cbioportal_url )
        error "FETCH_CBIOPORTAL needs --cbioportal_study (datahub id) or --cbioportal_url."
    def url = params.cbioportal_url ?: "https://cbioportal-datahub.s3.amazonaws.com/${params.cbioportal_study}.tar.gz"
    """
    echo "Downloading cBioPortal study bundle..."
    curl -fsSL '${url}' -o study.tar.gz
    mkdir -p study
    tar -xzf study.tar.gz -C study
    rm -f study.tar.gz
    # locate the mutations MAF (data_mutations.txt or data_mutations_extended.txt)
    src=\$(find study -type f \\( -name 'data_mutations.txt' -o -name 'data_mutations_extended.txt' \\) | head -1)
    test -n "\$src" || { echo "no data_mutations file in study" >&2; exit 1; }
    cp "\$src" cbioportal_mutations.maf
    echo "MAF: \$(wc -l < cbioportal_mutations.maf) lines from \$src."
    """
    stub:
    """ printf 'Hugo_Symbol\\tHGVSp_Short\\tTumor_Sample_Barcode\\n' > cbioportal_mutations.maf """
}


// ──────────────────────────────────────────────────────────────────────────
// FETCH_UNIPROT_SPROT_DAT — Swiss-Prot flat file (all features + Pfam DR refs)
// Replaces ~37 k per-protein UniProt REST API calls in ANNOTATION_MAP.
// The flat file contains every FT feature line (signal, transmembrane, binding,
// active site, region, etc.) for all Swiss-Prot entries in one ~700 MB download.
// storeDir-cached; PARSE_UNIPROT_DAT extracts the feature TSVs from this file.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_UNIPROT_SPROT_DAT {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/uniprot" : "${params.ref_dir}/uniprot" }

    output:
    path 'uniprot_sprot.dat.gz', emit: dat

    script:
    """
    echo "Downloading UniProt Swiss-Prot flat file (~700 MB compressed)..."
    curl -fsSL \\
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz' \\
        -o uniprot_sprot.dat.gz
    echo "Done — \$(du -sh uniprot_sprot.dat.gz | cut -f1)."
    """

    stub:
    """
    # Minimal stub with one human entry featuring SIGNAL + BINDING + Pfam DR
    python3 -c "
import gzip
s = (
    'ID   TEST_HUMAN   Reviewed;   100 AA.\\n'
    'AC   P12345;\\n'
    'OX   NCBI_TaxID=9606;\\n'
    'FT   SIGNAL          1..25\\n'
    'FT                   /evidence=\"ECO:0000255\"\\n'
    'FT   BINDING         70..70\\n'
    'FT                   /ligand=\"ATP\"\\n'
    'DR   Pfam; PF00001; 7tm_1; 1.\\n'
    '//\\n'
)
with gzip.open('uniprot_sprot.dat.gz', 'wt', encoding='latin-1') as f:
    f.write(s)
"
    """
}


// ──────────────────────────────────────────────────────────────────────────
// FETCH_INTERPRO_PFAM — InterPro protein→Pfam domain mapping with positions.
// Replaces ~37 k per-protein InterPro REST API calls in ANNOTATION_MAP.
// protein2ipr.dat.gz maps every UniProt accession to all InterPro member DB
// signatures (filtered to Pfam in PARSE_UNIPROT_DAT) including start/end.
// File size: ~1.2 GB compressed / ~9 GB uncompressed (all species). Filtered
// to our accession set during parsing so only ~20 k rows are kept.
// storeDir-cached.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_INTERPRO_PFAM {
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/interpro" : "${params.ref_dir}/interpro" }

    output:
    path 'protein2ipr.dat.gz', emit: dat

    script:
    """
    echo "Downloading InterPro protein2ipr.dat.gz (~1.2 GB compressed)..."
    curl -fsSL \\
        'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz' \\
        -o protein2ipr.dat.gz
    echo "Done — \$(du -sh protein2ipr.dat.gz | cut -f1)."
    """

    stub:
    """
    # Minimal stub: one Pfam domain row for P12345
    python3 -c "
import gzip
row = 'P12345\\tmd5\\t100\\tPfam\\tPF00001\\t7tm_1\\t1\\t96\\t1e-50\\tT\\t01-01-2024\\tIPR000001\\tTest\\n'
with gzip.open('protein2ipr.dat.gz', 'wt') as f:
    f.write(row)
"
    """
}
