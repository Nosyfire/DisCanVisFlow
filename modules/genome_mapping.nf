/*
 * modules/genome_mapping.nf — Module 3: Genome Mapping
 *
 * Maps protein sequences to genomic coordinates via:
 *   1. BLAT_ALIGN  — aligns cDNA transcripts to hg38 genome → PSL file
 *   2. GENOME_MAP  — builds protein AA → cDNA → genomic coordinate map
 *
 * Dependencies (all in Docker image):
 *   blat, pslCDnaFilter, twoBitToFa, tblastn (BLAST+)
 *
 * Process: SUBSET_FASTA (imported in main.nf with alias SUBSET_CDNA)
 *   is applied to the cDNA FASTA before BLAT_ALIGN.
 */


// ── Process 0: SPLIT_CDNA_FASTA ──────────────────────────────────────────────
//
// Split the cDNA FASTA into N roughly equal chunks so BLAT_ALIGN can run in
// parallel.  Only used when params.blat_chunks > 1; main.nf passes the full
// FASTA directly when blat_chunks == 1 (default).
//
process SPLIT_CDNA_FASTA {
    tag  { "split_cdna_${n_chunks}" }
    label 'process_low'

    input:
    path  cdna_fasta
    val   n_chunks

    output:
    path 'chunk_*.fasta', emit: chunks

    script:
    """
    python3 - <<'PYEOF'
from pathlib import Path
seqs, header, buf = [], None, []
with open("${cdna_fasta}") as f:
    for line in f:
        if line.startswith(">"):
            if header is not None:
                seqs.append((header, "".join(buf)))
            header, buf = line, []
        else:
            buf.append(line)
    if header is not None:
        seqs.append((header, "".join(buf)))
n = int("${n_chunks}")
chunk_size = max(1, -(-len(seqs) // n))   # ceiling division
for i in range(n):
    chunk = seqs[i * chunk_size:(i + 1) * chunk_size]
    if not chunk:
        break
    with open(f"chunk_{i:04d}.fasta", "w") as out:
        for h, s in chunk:
            out.write(h + s)
PYEOF
    """

    stub:
    """
    cp ${cdna_fasta} chunk_0000.fasta
    """
}


// ── Process 1: BLAT_ALIGN ────────────────────────────────────────────────────
//
// Align one cDNA chunk against hg38 using BLAT.  Outputs a raw PSL file named
// after the input chunk (unique names allow safe collect() into MERGE_BLAT_PSL).
//
// pslCDnaFilter is intentionally NOT run here — MERGE_BLAT_PSL runs it once on
// the merged output so -bestOverlap can see all hits for every query.
//
// BLAT parameters:
//   -t=dna      target = genomic DNA (2bit file)
//   -q=rna      query  = spliced RNA / cDNA
//   -fine       fine-grained alignment (better for short exons)
//   -minIdentity= minimum percentage of bases that must match
//   -out=psl    output in PSL format
//
process BLAT_ALIGN {
    tag "blat_${cdna_fasta.simpleName}"
    // BLAT is single-threaded; claim only 1 CPU so N chunks can run in parallel.
    // Memory: hg38.2bit loads into ~4 GB RAM per process.
    cpus   1
    memory '5 GB'

    input:
    path cdna_fasta
    path hg38_2bit

    output:
    path "${cdna_fasta.simpleName}.psl", emit: psl

    script:
    def min_id = params.blat_min_identity ?: 95
    """
    echo "Running BLAT: ${cdna_fasta} vs ${hg38_2bit}"
    blat \\
        ${hg38_2bit} \\
        ${cdna_fasta} \\
        ${cdna_fasta.simpleName}.psl \\
        -t=dna \\
        -q=rna \\
        -fine \\
        -minIdentity=${min_id} \\
        -out=psl
    echo "Raw BLAT hits: \$(tail -n +6 ${cdna_fasta.simpleName}.psl | wc -l)"
    """

    stub:
    """
    printf 'psLayout version 3\\n\\nmatch\\t...\\n     \\t...\\n---\\n' > ${cdna_fasta.simpleName}.psl
    printf '0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t+\\tENST00000423430.6|ENSG00000132155.11|...\\t1944\\t0\\t1944\\tchr3\\t198295559\\t12600000\\t12800000\\t1\\t1944,\\t0,\\t12600000,\\n' \\
        >> ${cdna_fasta.simpleName}.psl
    """
}


// ── Process 1b: MERGE_BLAT_PSL ───────────────────────────────────────────────
//
// Merge all raw PSL chunks produced by BLAT_ALIGN (one per chunk), then run
// pslCDnaFilter on the combined file.  -bestOverlap needs to see all hits for
// every query, so filtering must happen after merging.
//
// pslCDnaFilter flags:
//   -minId        minimum fractional identity
//   -minCover     minimum coverage of query bases
//   -bestOverlap  keep only the best overlapping alignment per query
//
process MERGE_BLAT_PSL {
    tag  { "merge_blat_psl" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/genome"
                                : "${params.outdir}/final/genome" },
        mode: 'copy', pattern: '*.psl'
    )

    input:
    path psls   // collected list of per-chunk raw PSL files (all uniquely named)

    output:
    path 'combined_output.psl', emit: psl

    script:
    def min_id  = params.blat_min_identity ?: 95
    def min_cov = params.blat_min_coverage ?: 25
    """
    # Merge: PSL header is 5 lines; take it from the first file, skip it in the rest
    sorted_psls=(\$(ls *.psl | sort))
    head -5 "\${sorted_psls[0]}" > merged_raw.psl
    for f in "\${sorted_psls[@]}"; do
        tail -n +6 "\$f" >> merged_raw.psl
    done
    echo "Merged raw BLAT hits: \$(tail -n +6 merged_raw.psl | wc -l)"

    pslCDnaFilter \\
        -minId=0.${min_id} \\
        -minCover=0.${min_cov} \\
        -localNearBest=0.001 \\
        -minQSize=20 \\
        -minNonRepSize=16 \\
        -ignoreNs \\
        -bestOverlap \\
        merged_raw.psl \\
        combined_output.psl
    echo "BLAT alignments after filter: \$(grep -v '^#' combined_output.psl | wc -l)"
    """

    stub:
    """
    printf '0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t+\\tENST00000423430.6|ENSG00000132155.11|...\\t1944\\t0\\t1944\\tchr3\\t198295559\\t12600000\\t12800000\\t1\\t1944,\\t0,\\t12600000,\\n' \\
        > combined_output.psl
    """
}


// ── Process 2: GENOME_MAP ────────────────────────────────────────────────────
//
// For each transcript in the PSL file, build a per-residue coordinate map:
//   protein AA index → cDNA codon positions → genomic positions
//
// The map format (combined_map.map) is a text file with one block per
// transcript:
//   # Qname Tname strand Tstart-Tend
//   0 M 0,1,2 ATG M    12600000,12600001,12600002, ATG M
//   1 A 3,4,5 GCT A    ...
//   ...
//
// Python worker: create_genome_map_worker.py
//
process GENOME_MAP {
    tag "gmap_${psl.simpleName}"
    label 'process_high'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/genome"
                                : "${params.outdir}/final/genome" },
        mode: 'copy', pattern: '*.map'
    )
    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/genome"
                                : "${params.outdir}/final/genome" },
        mode: 'copy', pattern: '*.txt'
    )

    input:
    path psl
    path cdna_fasta
    path prot_fasta
    path loc_chrom_tsv
    path hg38_2bit

    output:
    path 'combined_map.map',  emit: map_file
    path 'error_map.txt',     optional: true, emit: errors

    script:
    def cpus = task.cpus ?: 4
    """
    create_genome_map_worker.py \\
        --psl           ${psl} \\
        --cdna_fasta    ${cdna_fasta} \\
        --prot_fasta    ${prot_fasta} \\
        --loc_chrom     ${loc_chrom_tsv} \\
        --hg38_2bit     ${hg38_2bit} \\
        --output_dir    . \\
        --num_processes ${cpus}
    """

    stub:
    """
    printf '# ENST00000423430.6|ENSG00000132155.11 chr3 + 12600000-12800000\\n' > combined_map.map
    printf '0 M 0,1,2 ATG M    12600000,12600001,12600002, ATG M\\n'            >> combined_map.map
    touch error_map.txt
    """
}

// ──────────────────────────────────────────────────────────────────────────
// GENOME_QUERY_MAP — genome ↔ protein reference tables (for the DisCanVis2 web layer)
//   genome_protein_index.tsv      per-nucleotide chrom/gpos ↔ Protein_ID/prot_pos
//   genome_protein_mutations.tsv  EVERY possible SNV (reference), genome ↔ protein
// Python worker: create_genome_query_worker.py
// ──────────────────────────────────────────────────────────────────────────
process GENOME_QUERY_MAP {
    tag  { "genome_query" }
    label 'process_low'

    publishDir(
        path: { params.gene_dir ? "${params.outdir}/${params.gene_dir}/final/genome"
                                : "${params.outdir}/final/genome" },
        mode: 'copy'
    )

    input:
    path combined_map

    output:
    path 'genome_protein_index.tsv',     emit: index
    path 'genome_protein_mutations.tsv', emit: mutations

    script:
    """
    create_genome_query_worker.py \\
        --combined_map ${combined_map} \\
        --outdir .
    """

    stub:
    """
    printf 'chrom\\tgpos\\tstrand\\tProtein_ID\\tprot_pos\\tcodon_offset\\taa\\tcodon\\n' > genome_protein_index.tsv
    printf 'chrom\\tgpos\\tstrand\\tref\\talt\\tProtein_ID\\tprot_pos\\tcodon_offset\\tref_codon\\talt_codon\\tref_aa\\talt_aa\\tconsequence\\thgvs_g\\thgvs_p\\n' > genome_protein_mutations.tsv
    """
}
