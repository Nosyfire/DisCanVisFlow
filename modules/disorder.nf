/*
 * modules/disorder.nf — Module 5 disorder — IUPred/ANCHOR/AIUPred/AlphaFold/MobiDB + RSA/position
 *
 * Processes: MOBIDB_MAP, FETCH_MOBIDB, PARSE_ALPHAFOLD_PLDDT, DISORDER_MAP, POSITION_BASED_MAP
 * (split out of the former annotation_mapping.nf monolith)
 */



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
// DISPROT_MAP  — DisProt curated disorder regions → Protein_ID TSV.
//
// Maps DisProt manually-curated intrinsic-disorder regions (UniProt-accession
// keyed, IDPO/GO ontology terms) onto every GENCODE isoform, validating the
// canonical coordinates against each isoform's sequence.
// ──────────────────────────────────────────────────────────────────────────
process DISPROT_MAP {
    tag  { "disprot_map" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disorder"
                                : "${params.outdir}/final/disorder" },
        mode: 'copy'
    )

    input:
    path loc_chrom
    path disprot_tsv   // disprot_regions.tsv from FETCH_DISPROT (or NO_FILE)

    output:
    path "disprot.tsv", emit: disprot

    script:
    """
    create_disprot_worker.py \\
        --seq_table   ${loc_chrom} \\
        --disprot_tsv ${disprot_tsv} \\
        --outdir      . \\
        --only_main_isoforms
    """

    stub:
    """
    echo -e "Protein_ID\tEntry_Isoform\tdisprot_id\tregion_id\tstart\tend\tterm_namespace\tterm_id\tterm_name\teco_id\tpmid\tdataset" > disprot.tsv
    """
}


// ──────────────────────────────────────────────────────────────────────────
// FETCH_DISPROT  — bulk DisProt release (IDPO + GO term ontologies), cached.
// ──────────────────────────────────────────────────────────────────────────
process FETCH_DISPROT {
    tag  { "disprot_current" }
    label 'process_low'
    storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/disprot" : "${params.ref_dir}/disprot" }

    output:
    path "disprot_regions.tsv", emit: disprot_tsv

    script:
    """
    python3 - << 'PYEOF'
import urllib.request, sys, time

url = ('https://disprot.org/api/v2/download?format=tsv&release=current'
       '&term_ontology=IDPO&term_ontology=GO')
out = 'disprot_regions.tsv'
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
    """

    stub:
    """
    printf 'UniProt ACC\\tDisProt ID\\tRegion ID\\tStart\\tEnd\\tTerm namespace\\tTerm ID\\tTerm name\\tECO Term ID\\tPMID\\tRegion sequence\\tObsolete\\n' > disprot_regions.tsv
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


process DISORDER_MAP {
    tag  { "disorder_map" }
    label 'process_medium'
    // When scattering (scatter_chunks > 1) this runs per-chunk and a MERGE step
    // publishes the combined tables, so per-task publishing is disabled.
    // AlphaFoldTable.tsv is a structural feature → final/structure; all other
    // disorder tables → final/disorder.
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/disorder"
                                : "${params.outdir}/final/disorder" },
        mode: 'copy',
        enabled: ( ((params.scatter_chunks ?: 1) as Integer) <= 1 ),
        saveAs: { fn -> fn == "AlphaFoldTable.tsv" ? null : fn }
    )
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/structure"
                                : "${params.outdir}/final/structure" },
        mode: 'copy',
        enabled: ( ((params.scatter_chunks ?: 1) as Integer) <= 1 ),
        saveAs: { fn -> fn == "AlphaFoldTable.tsv" ? fn : null }
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
// POSITION_BASED_MAP  — Module 5m: per-residue PositionBasedAnnotations +
//                       RSAscores (RSA derived from pLDDT)
//
// Aggregates: IUPred3 · pLDDT · CombinedDisorder · GOPHER · phastCons · Pfam
// Output lives in mapped/ (Protein_ID-keyed).
// ──────────────────────────────────────────────────────────────────────────
process POSITION_BASED_MAP {
    tag  { "position_based_map" }
    label 'process_medium'

    // position_based_annotations.tsv → final/position ; rsa_scores.tsv → final/structure
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/position"
                                : "${params.outdir}/final/position" },
        mode: 'copy',
        pattern: "position_based_annotations.tsv"
    )
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/structure"
                                : "${params.outdir}/final/structure" },
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
