/*
 * modules/annotation_mapping.nf  — Module 5: Annotation Mapping
 *
 * Output structure (when params.gene_dir is set, e.g. 'chr3/raf1'):
 *   results/{gene_dir}/unmapped/annotations/   ← Entry_Isoform-keyed (ELM/DIBS/PTM/etc.) → input to TRANSCRIPT_MAP
 *   results/{gene_dir}/unmapped/disorder/       ← per-residue disorder scores (Entry_Isoform-keyed) → input to TRANSCRIPT_MAP
 *   results/{gene_dir}/unmapped/pdb/            ← PDB structure mapping (UniProt-keyed)
 *   results/{gene_dir}/mapped/annotations/      ← Protein_ID-keyed: TRANSCRIPT_MAP outputs + GO/PPI/PEM/Scansite/SNP/Polymorphism/CoiledCoils
 *   results/{gene_dir}/mapped/disorder/         ← Protein_ID-keyed: CombinedDisorder pass-through from TRANSCRIPT_MAP
 *   results/{gene_dir}/mapped/conservation/     ← Protein_ID-keyed: GOPHER + phastCons
 *   results/{gene_dir}/mapped/disease/          ← Protein_ID-keyed: ClinVar + OMIM disease ontology
 *   results/{gene_dir}/mapped/drivers/          ← Protein_ID-keyed: CGC census + compendium
 *   results/{gene_dir}/mapped/pathogenicity/    ← Protein_ID-keyed: dbNSFP + AlphaMissense
 *   results/{gene_dir}/mapped/mutations/        ← Protein_ID-keyed: ClinVar/TCGA/MAF/DepMap mutations
 *
 * Processes
 * ─────────
 *  FETCH_ELM       (fallback downloader if legacy file not used)
 *  FETCH_MOBIDB    Download MobiDB curated disorder — cached in ref_dir.
 *  FETCH_GO        Download GOA human annotation + GO OBO — cached in ref_dir.
 *
 *  ANNOTATION_MAP  ELM + DIBS + MFIB + PhasePro (legacy_data/) + UniProt REST
 *                  + PTMdb + PhosphoSite (legacy_data/ptm/) + Pfam (InterPro)
 *                  + isoform annotation transfer (100% identity)
 *
 *  DISORDER_MAP    IUPred3 + ANCHOR2 + AIUPred + AIUPred-Binding (local libs)
 *                  + AlphaFold pLDDT (API) + combined disorder
 *
 *  GO_MAP          GO term annotation (GOA human + GO OBO)
 *
 *  POLYMORPHISM_MAP Natural variant / SNP annotation via UniProt REST API
 *
 *  PEM_MAP         Predicted ELM Motifs (HotspotPEM supplementary dataset)
 *  PEM_TRANSFER_MAP Map PEM motifs to alternative isoforms (sequence homology)
 *
 *  COILEDCOILS_MAP DeepCoil coiled-coil predictions
 *
 *  CONSERVATION_MAP GOPHER multi-level + phastCons per-residue scores (Module 7)
 *
 *  SCANSITE_MAP         ScanSite 4.0 kinase/phospho motifs (Module 5k)
 *  SNP_MAP              Population SNP polymorphisms (Module 5l)
 *  CLINVAR_DISEASE_MAP  ClinVar disease ontology (Module 8a, filter fallback)
 *  CLINVAR_DISEASE_BUILD Build ClinVar disease from MONDO OBO + MUTATION_MAP (Module 8a)
 *  OMIM_MAP             OMIM disease ontology (Module 8b)
 *  CANCER_DRIVER_MAP    CGC census + compendium (Module 8c)
 *  ALPHAMISSENSE_MAP    AlphaMissense isoform scores (Module 8d)
 *  DEPMAP_MAP           DepMap cancer cell line mutations (Module 8e)
 *  DBNSFP_MAP           Raw dbNSFP chr*.gz → Protein_ID via combined_map.map (Module 8f)
 *  PATHOGENICITY_MAP    Pre-mapped dbNSFP filter fallback (Module 8f)
 *
 *  PDB_MAP         PDBe REST API → per-residue + per-region ordered/disordered
 *
 *  EXON_MAP        Exon boundaries from combined_map.map (>20 bp gap)
 *
 *  TRANSCRIPT_MAP  Map UniProt-keyed annotations → each Gencode transcript
 *                  Marks homology_transfer=True where 100% sequence match
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
// MOBIDB_MAP  — Module 5o: MobiDB disorder features → MobiDBDisorder TSV
//
// Takes the bulk MobiDB TSV (acc | feature | source | start..end) and maps
// each feature to all GENCODE isoforms via Entry_Isoform, computing
// content_fraction, content_count, and total disordered length per feature.
// ──────────────────────────────────────────────────────────────────────────
process MOBIDB_MAP {
    tag  { "mobidb_map" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/intermediate/disorder"
                                : "${params.outdir}/intermediate/disorder" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path mobidb_tsv    // mobidb_human.tsv from FETCH_MOBIDB (or NO_FILE)

    output:
    path "mobidb_disorder.tsv", emit: mobidb_disorder

    script:
    """
    create_mobidb_worker.py \\
        --seq_table  ${loc_chrom} \\
        --mobidb_tsv ${mobidb_tsv} \\
        --outdir     .
    """

    stub:
    """
    echo -e "Protein_ID\tEntry_Isoform\tfeature\tstart_end\tcontent_fraction\tcontent_count\tlength" > mobidb_disorder.tsv
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
// FETCH_MOBIDB  — curated disorder for homo sapiens (separate projection)
// ──────────────────────────────────────────────────────────────────────────
process FETCH_MOBIDB {
    tag  { "mobidb_9606" }
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/mobidb" : "${params.ref_dir}/mobidb" }

    output:
    path "mobidb_human.tsv", emit: mobidb_tsv

    script:
    """
    python3 - << 'PYEOF'
import urllib.request, sys, time

downloads = [
    ('https://mobidb.bio.unipd.it/api/download?ncbi_taxon_id=9606&projection=curated-disorder-merge&format=tsv',
     'mobidb_curated.tsv'),
    ('https://mobidb.bio.unipd.it/api/download?ncbi_taxon_id=9606&projection=homology-disorder-merge&format=tsv',
     'mobidb_homol.tsv'),
]

for url, out in downloads:
    for attempt in range(1, 4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'DisCanVisFlow/0.6.0'})
            with urllib.request.urlopen(req, timeout=1800) as resp, open(out, 'wb') as fh:
                fh.write(resp.read())
            lines = sum(1 for _ in open(out))
            if lines < 2:
                print(f'ERROR: {out} has only {lines} line(s) — server returned no data', file=sys.stderr)
                sys.exit(1)
            print(f'{out}: {lines} lines OK')
            break
        except Exception as exc:
            print(f'Attempt {attempt}/3 failed for {out}: {exc}', file=sys.stderr)
            if attempt == 3:
                sys.exit(1)
            time.sleep(10)
PYEOF

    cat mobidb_curated.tsv mobidb_homol.tsv | sort -u > mobidb_human.tsv
    """

    stub:
    """
    printf 'acc\\tfeature\\tsource\\tstart..end\\n' > mobidb_human.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// PARSE_UNIPROT_DAT — extract feature TSVs from the Swiss-Prot flat file +
// InterPro protein2ipr file.  Runs once per dat.gz (storeDir-cached).
// Outputs consumed by ANNOTATION_MAP in place of per-protein REST API calls.
// ──────────────────────────────────────────────────────────────────────────
process PARSE_UNIPROT_DAT {
    tag  { "parse_uniprot_dat" }
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/uniprot_parsed"
                                : "${params.ref_dir}/uniprot_parsed" }

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
// PARSE_ALPHAFOLD_PLDDT — extract per-residue pLDDT scores from the EBI
// AlphaFold human-proteome bulk tar (UP000005640_9606_HUMAN_v*.tar).
// Replaces per-accession EBI API calls in DISORDER_MAP.
// The tar contains AF-{acc}-F1-confidences_v*.json.gz files with pLDDT lists.
// Output: alphafold_plddt.tsv (Accession, Plldtscores comma-separated).
// storeDir-cached.
// ──────────────────────────────────────────────────────────────────────────
process PARSE_ALPHAFOLD_PLDDT {
    tag  { "parse_alphafold_plddt" }
    label 'process_medium'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/alphafold_parsed"
                                : "${params.ref_dir}/alphafold_parsed" }

    input:
    path af_tar   // UP000005640_9606_HUMAN_v*.tar

    output:
    path 'alphafold_plddt.tsv', emit: plddt

    script:
    """
    python3 - <<'PYEOF'
import gzip, io, tarfile, re, csv, sys, json
from pathlib import Path

# v4 tar: AF-{acc}-F1-confidences_v4.json.gz   → parse JSON confidenceScore/pLDDT
# v5/v6 tar: AF-{acc}-F1-model_v6.pdb.gz       → parse B-factor (pLDDT) from CA atoms
PDB_RE   = re.compile(r'AF-([A-Z0-9]+(?:-\\d+)?)-F\\d+-model_v\\d+\\.pdb\\.gz')
JSON_RE  = re.compile(r'AF-([A-Z0-9]+(?:-\\d+)?)-F\\d+-confidences_v\\d+\\.json\\.gz')

def scores_from_pdb_gz(raw_bytes):
    # pLDDT = B-factor (col 60-66) of CA atoms
    text   = gzip.decompress(raw_bytes).decode('utf-8', errors='ignore')
    scores = []
    for line in text.splitlines():
        if len(line) < 66 or not line.startswith('ATOM'):
            continue
        if line[12:16].strip() != 'CA':
            continue
        try:
            scores.append(round(float(line[60:66]), 1))
        except ValueError:
            pass
    return scores

def scores_from_json_gz(raw_bytes):
    # AlphaFold v4 confidence JSON format
    data = json.loads(gzip.decompress(raw_bytes))
    return data.get('confidenceScore') or data.get('pLDDT') or []

out_path = Path('alphafold_plddt.tsv')
written  = 0

with open(out_path, 'w', newline='') as fout:
    writer = csv.writer(fout, delimiter='\\t')
    writer.writerow(['Accession', 'Plldtscores'])
    with tarfile.open('${af_tar}', 'r') as tar:
        for member in tar:
            name = member.name
            m_pdb  = PDB_RE.search(name)
            m_json = JSON_RE.search(name)
            if not (m_pdb or m_json):
                continue
            raw = tar.extractfile(member)
            if raw is None:
                continue
            raw_bytes = raw.read()
            acc    = (m_pdb or m_json).group(1)
            scores = (scores_from_pdb_gz(raw_bytes) if m_pdb
                      else scores_from_json_gz(raw_bytes))
            if scores:
                writer.writerow([acc, ','.join(f'{v:.1f}' for v in scores)])
                written += 1

print(f'alphafold_plddt.tsv: {written} proteins written', file=sys.stderr)
PYEOF
    """

    stub:
    """
    printf 'Accession\\tPlldtscores\\n'                           > alphafold_plddt.tsv
    printf 'P04637\\t90.0,88.5,92.1,87.3,95.0\\n' >> alphafold_plddt.tsv
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


process DISORDER_MAP {
    tag  { "disorder_map" }
    label 'process_medium'
    // When scattering (scatter_chunks > 1) this runs per-chunk and a MERGE step
    // publishes the combined tables, so per-task publishing is disabled.
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disorder"
                                : "${params.outdir}/final/disorder" },
        mode: 'copy',
        enabled: ( ((params.scatter_chunks ?: 1) as Integer) <= 1 )
    )

    input:
    path loc_chrom
    path mobidb_tsv,            stageAs: 'mobidb_dis.tsv'
    path pfam_tsv,              stageAs: 'pfam_dis.tsv'
    path alphafold_plddt,       stageAs: 'alphafold_plddt_in.tsv'     // from PARSE_ALPHAFOLD_PLDDT or NO_FILE
    path af_precomputed_table,  stageAs: 'af_precomputed_table.tsv'   // pre-computed Protein_ID|Plldtscores or NO_FILE
    path setup_done,            stageAs: 'setup.done'                  // sentinel from SETUP_DEPS (ordering)
    path aiupred_py_file,       stageAs: 'aiupred_py.txt'             // detected path from SETUP_DEPS

    output:
    path "IUPredscores.tsv",         emit: iupred
    path "Anchorscores.tsv",         emit: anchor
    path "AIUPredscores.tsv",        emit: aiupred
    path "AIUPredBinding.tsv",       emit: aiupred_binding
    path "AlphaFoldTable.tsv",       emit: plddt
    path "CombinedDisorderNew.tsv",  emit: disorder_regions
    path "CombinedDisorderNew_Pos.tsv", emit: disorder_pos

    script:
    def skip_af  = params.skip_alphafold ? "--skip_alphafold" : ""
    def skip_iu  = params.skip_iupred  ? "--skip_iupred"  : ""
    def skip_aiu = params.skip_aiupred ? "--skip_aiupred" : ""
    def mob_arg  = (mobidb_tsv.name != 'NO_FILE') ? "--mobidb_tsv mobidb_dis.tsv" : ""
    def pfam_arg = (pfam_tsv.name   != 'NO_FILE') ? "--pfam_tsv pfam_dis.tsv"     : ""
    // pLDDT source priority: bulk-tar parsed (Accession-keyed) > pre-computed merged (Protein_ID-keyed)
    def af_local = (alphafold_plddt.name      != 'NO_FILE') ? "--alphafold_plddt_tsv alphafold_plddt_in.tsv"
                 : (af_precomputed_table.name != 'NO_FILE') ? "--alphafold_plddt_tsv af_precomputed_table.tsv"
                 : ""
    def aiupred_param = params.aiupred_python ?: ""
    """
    # Resolve aiupred_python: explicit param > SETUP_DEPS file > empty (graceful skip)
    _aiupred_py="${aiupred_param}"
    [[ -z "\${_aiupred_py}" ]] && [[ -f aiupred_py.txt ]] && \\
        _aiupred_py="\$(cat aiupred_py.txt | tr -d '[:space:]')"

    create_disorder_worker.py \\
        --loc_chrom      ${loc_chrom} \\
        ${mob_arg} ${pfam_arg} \\
        ${af_local} \\
        --ext_programs   ${params.ext_programs} \\
        --aiupred_python "\${_aiupred_py:-python}" \\
        --output_dir     . \\
        --request_delay  ${params.disorder_api_delay ?: 0.5} \\
        ${skip_af} ${skip_iu} ${skip_aiu}
    """

    stub:
    """
    for f in IUPredscores.tsv Anchorscores.tsv AIUPredscores.tsv \\
              AIUPredBinding.tsv AlphaFoldTable.tsv \\
              CombinedDisorderNew.tsv CombinedDisorderNew_Pos.tsv; do
        touch "\$f"
    done
    """
}


// ──────────────────────────────────────────────────────────────────────────
// PDB_MAP  — PDBe REST API: metadata + per-residue + per-region
// ──────────────────────────────────────────────────────────────────────────
process PDB_MAP {
    tag  { "pdb_map" }
    label 'process_medium'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pdb"
                                : "${params.outdir}/final/pdb" },
        mode: 'copy'
    )

    input:
    path loc_chrom

    output:
    path "pdb_structures.tsv", emit: structures
    path "pdb_missing.tsv",    emit: pdb_missing

    script:
    """
    create_pdb_worker.py \\
        --loc_chrom     ${loc_chrom} \\
        --output_dir    . \\
        --request_delay ${params.annotation_api_delay ?: 0.5}
    """

    stub:
    """
    printf 'Protein_ID\\tAccession\\tpdb_id\\tchain_id\\tstruct_asym_id\\tentity_id\\tprot_start\\tprot_end\\tunp_start\\tunp_end\\tresolution\\texperimental_method\\n' > pdb_structures.tsv
    printf 'Protein_ID\\tAccession\\tpdb_id\\tchain_id\\tprot_start\\tprot_end\\tunp_start\\tunp_end\\tlength\\n' > pdb_missing.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// EXON_MAP  — exon boundaries from combined_map.map (>20 bp gap)
// ──────────────────────────────────────────────────────────────────────────
process EXON_MAP {
    tag  { "exon_map" }
    label 'process_low'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/genome"
                                : "${params.outdir}/final/genome" },
        mode: 'copy'
    )

    input:
    path combined_map
    path loc_chrom

    output:
    path "exon.tsv", emit: exon

    script:
    """
    create_exon_worker.py \\
        --combined_map  ${combined_map} \\
        --loc_chrom     ${loc_chrom} \\
        --output_dir    .
    """

    stub:
    """
    printf 'Protein_ID\\texon_number\\ttotal_exons\\taa_start\\taa_end\\taa_length\\tgenomic_start\\tgenomic_end\\n' > exon.tsv
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

// ──────────────────────────────────────────────────────────────────────────
// PROTEINGYM_MAP  — Module 8g: ProteinGym DMS scores + pathogenicity proxy
// ──────────────────────────────────────────────────────────────────────────
process PROTEINGYM_MAP {
    tag  { "proteingym" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom_seq
    path proteingym_tsv   // premapped: Protein_ID-keyed table; uniprot: raw UniProt-keyed table; or NO_FILE

    output:
    path "proteingym.tsv", emit: proteingym

    script:
    // uniprot mode = a freshly fetched UniProt-keyed raw table (FETCH_PROTEINGYM
    // or --proteingym_raw); otherwise filter the pre-mapped Protein_ID-keyed table.
    def uniprot_mode = (params.fetch_proteingym || params.proteingym_raw) as boolean
    def src_arg = uniprot_mode ? "--mapping_mode uniprot --proteingym_raw ${proteingym_tsv}"
                               : "--proteingym ${proteingym_tsv}"
    """
    create_proteingym_worker.py \\
        --seq_table   ${loc_chrom_seq} \\
        ${src_arg} \\
        --outdir      .
    """

    stub:
    """
    printf 'Protein_ID\\tProtein_position\\tprotein_variant\\tDMS_score\\tDMS_score_bin\\tDMS_id\\tuniprot_id\\tmapping_type\\n' > proteingym.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// SUMMARY_MAP  — Final annotation count summary
// ──────────────────────────────────────────────────────────────────────────
process SUMMARY_MAP {
    tag  { "annotation_summary" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}"
                                : "${params.outdir}" },
        mode: 'copy'
    )

    input:
    path results_dir   // the gene output directory

    output:
    path "annotation_summary.tsv", emit: summary

    script:
    def gene = params.gene_dir ? params.gene_dir.split("/")[-1].toUpperCase() : "proteome"
    """
    create_summary_worker.py \\
        --gene_name   ${gene} \\
        --results_dir ${results_dir} \\
        --outdir      .
    """

    stub:
    """
    echo -e "gene\tannotation_type\tcount\tnote" > annotation_summary.tsv
    """
}

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

// ──────────────────────────────────────────────────────────────────────────
// ALPHAMISSENSE_MAP  — Module 8d: AlphaMissense GENCODE isoform scores
// ──────────────────────────────────────────────────────────────────────────
process ALPHAMISSENSE_MAP {
    tag  { "alphamissense" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path alphamissense_gz   // AlphaMissense_isoforms_hg38.tsv (decompressed by DECOMPRESS_ALPHAMISSENSE) or NO_FILE

    output:
    path "alphamissense.tsv", emit: alphamissense

    script:
    """
    create_alphamissense_worker.py \\
        --seq_table        ${loc_chrom} \\
        --alphamissense_gz ${alphamissense_gz} \\
        --outdir           .
    """

    stub:
    """
    echo -e "Protein_ID\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class" > alphamissense.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// DEPMAP_MAP  — Module 8e: DepMap cancer cell line somatic mutations
// ──────────────────────────────────────────────────────────────────────────
process DEPMAP_MAP {
    tag  { "depmap" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/mutations/DepMap"
                                : "${params.outdir}/final/mutations/DepMap" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path depmap_tsv   // mapped_filtered_mutations.tsv from DepMap or NO_FILE

    output:
    path "depmap_mutations.tsv", emit: depmap

    script:
    """
    create_depmap_worker.py \\
        --seq_table ${loc_chrom} \\
        --depmap_tsv ${depmap_tsv} \\
        --outdir    .
    """

    stub:
    """
    echo -e "Protein_ID\tHGVSp_Short\tModelID" > depmap_mutations.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// DBNSFP_MAP  — Module 8f: map raw dbNSFP chr*.gz via combined_map.map
// ──────────────────────────────────────────────────────────────────────────
process DBNSFP_MAP {
    tag  { "dbnsfp_map" }
    label 'process_high'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path combined_map
    path dbnsfp_raw_dir

    output:
    path "dbnsfp_scores.tsv", emit: scores

    script:
    def hdr_arg = params.dbnsfp_bed_header ? "--dbnsfp_bed_header ${params.dbnsfp_bed_header}" : ""
    def n_cpu   = task.cpus ?: 1
    """
    create_dbnsfp_map_worker.py \\
        --seq_table       ${loc_chrom} \\
        --combined_map    ${combined_map} \\
        --dbnsfp_raw_dir  ${dbnsfp_raw_dir} \\
        ${hdr_arg} \\
        --n_cpu           ${n_cpu} \\
        --outdir          .
    """

    stub:
    """
    echo -e "Protein_ID\tProtein_position\tAlphaMissense_score" > dbnsfp_scores.tsv
    """
}

// ──────────────────────────────────────────────────────────────────────────
// PATHOGENICITY_MAP  — Module 8f: dbNSFP pathogenicity predictor scores
// ──────────────────────────────────────────────────────────────────────────
process PATHOGENICITY_MAP {
    tag  { "pathogenicity" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path dbnsfp_tsv   // dbNSFP_custom/mapped_filtered_mutations.tsv or NO_FILE

    output:
    path "pathogenicity_scores.tsv", emit: scores

    script:
    """
    create_pathogenicity_worker.py \\
        --seq_table  ${loc_chrom} \\
        --dbnsfp_tsv ${dbnsfp_tsv} \\
        --outdir     .
    """

    stub:
    """
    echo -e "Protein_ID\tProtein_position\tAlphaMissense_score\tCADD_phred\tSIFT_score" > pathogenicity_scores.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// POSITION_BASED_MAP  — Module 5m: per-residue PositionBasedAnnotations +
//                       RSAscores (RSA derived from pLDDT)
//
// Aggregates: IUPred3 · pLDDT · CombinedDisorder · GOPHER · phastCons · Pfam
// Output lives in mapped/ (Protein_ID-keyed).
// ──────────────────────────────────────────────────────────────────────────
process POSITION_BASED_MAP {
    tag  { "position_based_map" }
    label 'process_medium'

    // position_based_annotations.tsv → final/position ; rsa_scores.tsv → final/disorder
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/position"
                                : "${params.outdir}/final/position" },
        mode: 'copy',
        pattern: "position_based_annotations.tsv"
    )
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disorder"
                                : "${params.outdir}/final/disorder" },
        mode: 'copy',
        pattern: "rsa_scores.tsv"
    )

    input:
    path loc_chrom
    path iupred_tsv        // IUPredscores.tsv from DISORDER_MAP
    path plddt_tsv         // AlphaFoldTable.tsv from DISORDER_MAP
    path combined_pos_tsv  // CombinedDisorderNew_Pos.tsv from DISORDER_MAP (or TRANSCRIPT_MAP)
    path phastcons_tsv,    stageAs: 'phastcons/*'  // conservation_phastcons.tsv or NO_FILE
    path conservation_tsv, stageAs: 'conslevel/*'  // conservation_multiple_level.tsv or NO_FILE
    path pfam_tsv,         stageAs: 'pfam/*'        // pfam_domains.tsv from ANNOTATION_MAP or NO_FILE

    output:
    path "position_based_annotations.tsv", emit: pos_annotations
    path "rsa_scores.tsv",                 emit: rsa_scores

    script:
    """
    create_position_based_worker.py \\
        --seq_table          ${loc_chrom} \\
        --iupred_tsv         ${iupred_tsv} \\
        --plddt_tsv          ${plddt_tsv} \\
        --combined_pos_tsv   ${combined_pos_tsv} \\
        --phastcons_tsv      ${phastcons_tsv} \\
        --conservation_tsv   ${conservation_tsv} \\
        --pfam_tsv           ${pfam_tsv} \\
        --outdir             .
    """

    stub:
    """
    echo -e "Protein_ID\tposition\tplddt\trsa\tiupred\tedisorder\tcombineddisorder\tphastCons\tconservationGlobal\tconservationMammal\tconservationVertebrate\tconservationEukaryota\tconservationEumetazoa\tconservationOpisthokonta\tconservationViridiplantae\tpfam" > position_based_annotations.tsv
    echo -e "Protein_ID\trsascores" > rsa_scores.tsv
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
// MAVEDB_MAP — MaveDB single-mutant functional scores (already Protein_ID-keyed)
// Filters the large pre-mapped MaveDB table to this run's isoforms.
// Python worker: create_mavedb_worker.py
// ──────────────────────────────────────────────────────────────────────────
process MAVEDB_MAP {
    tag  { "mavedb" }
    label 'process_medium'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom_seq
    path mavedb_tsv      // premapped: Protein_ID-keyed table; uniprot: raw UniProt-keyed table; or NO_FILE

    output:
    path "mavedb.tsv", emit: mavedb

    script:
    // uniprot mode = a freshly fetched UniProt-keyed raw table (FETCH_MAVEDB or
    // --mavedb_raw); otherwise filter the pre-mapped Protein_ID-keyed table.
    def uniprot_mode = (params.fetch_mavedb || params.mavedb_raw) as boolean
    def src_arg = uniprot_mode ? "--mapping_mode uniprot --mavedb_raw ${mavedb_tsv}"
                               : "--mavedb ${mavedb_tsv}"
    """
    create_mavedb_worker.py \\
        --seq_table ${loc_chrom_seq} \\
        ${src_arg} \\
        --outdir    .
    """

    stub:
    """
    printf 'Protein_ID\\tProtein_position\\tprot_expr\\tscore\\tmavedb_id\\turn\\tgene_name\\tuniprot\\tTranscript_ID\\tis_double_mutant\\tmapping_type\\n' > mavedb.tsv
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
        ${src_args} \\
        --outdir .
    """

    stub:
    """
    printf '# Mapping summary (stub)\\n' > mapping_summary.md
    touch mapping_coverage.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// PDB_BULK_MAP — PDB structures via SIFTS bulk join (no per-protein PDBe API).
// ~1000× faster than PDB_MAP; reproduces the same pdb_structures.tsv coverage.
// pdb_missing (unobserved residues) needs the API, so it is emitted header-only
// here — run PDB_MAP too if you need pdb_missing. Worker: create_pdb_bulk_worker.py
// ──────────────────────────────────────────────────────────────────────────
process PDB_BULK_MAP {
    tag  { "pdb_bulk_map" }
    label 'process_low'
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pdb"
                                : "${params.outdir}/final/pdb" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path sifts_gz

    output:
    path "pdb_structures.tsv", emit: structures
    path "pdb_missing.tsv",    emit: pdb_missing

    script:
    """
    create_pdb_bulk_worker.py \\
        --seq_table ${loc_chrom} \\
        --sifts_tsv ${sifts_gz} \\
        --min_identity ${params.homology_min_identity ?: 0.9} \\
        --outdir    .
    """

    stub:
    """
    printf 'Protein_ID\\tAccession\\tpdb_id\\tchain_id\\tstruct_asym_id\\tentity_id\\tprot_start\\tprot_end\\tunp_start\\tunp_end\\tresolution\\texperimental_method\\n' > pdb_structures.tsv
    printf 'Protein_ID\\tAccession\\tpdb_id\\tchain_id\\tprot_start\\tprot_end\\tunp_start\\tunp_end\\tlength\\n' > pdb_missing.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// FINCHES_MAP  — Module 8h: site-wise Δε (LLPS-change score) via saturation
// mutagenesis using the FINCHES Mpipi_GGv1 force field.
//
// For every possible single-AA substitution at every position, computes:
//   Δε = epsilon(mut,mut) − epsilon(wt,wt)  (homotypic self-interaction energy)
// Positive Δε = mutation increases LLPS propensity.
//
// Output goes to final/pathogenicity/ alongside AlphaMissense / ProteinGym.
// ⚠  CC BY-NC 4.0 — non-commercial use only.
//
// Citation: Ginell et al. bioRxiv 2024.06.03.597104
// ──────────────────────────────────────────────────────────────────────────
process FINCHES_MAP {
    tag  { "finches_map" }
    label 'process_high'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/pathogenicity"
                                : "${params.outdir}/final/pathogenicity" },
        mode: 'copy'
    )

    input:
    path loc_chrom

    output:
    path "finches_saturation.tsv", emit: finches

    script:
    def finches_lib_arg = params.finches_lib         ? "--finches_lib ${params.finches_lib}"               : ''
    def n_cpu           = task.cpus                  ?: 1
    def only_main_arg   = params.only_main_isoforms  ? '--only_main_isoforms'                              : ''
    def max_len_arg     = params.finches_max_seq_len ? "--max_seq_len ${params.finches_max_seq_len}"        : '--max_seq_len 3000'
    // Use finches_python when set (needs finches + aiupred conda env); else fall
    // back to the script on PATH (assumes finches importable from default python3).
    def py              = params.finches_python       ?: null
    def invoke          = py ? "${py} ${baseDir}/bin/create_finches_worker.py" : "create_finches_worker.py"
    """
    ${invoke} \\
        --loc_chrom  ${loc_chrom} \\
        --output_dir . \\
        --n_cpu      ${n_cpu} \\
        ${finches_lib_arg} \\
        ${only_main_arg} \\
        ${max_len_arg}
    """

    stub:
    """
    printf 'Protein_ID\\tPosition\\tWT_AA\\tMut_AA\\tWT_Epsilon\\tMut_Epsilon\\tDelta_Epsilon\\n' > finches_saturation.tsv
    """
}
