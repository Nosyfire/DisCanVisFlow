/*
 * modules/structure.nf — Module 5 structure — PDB coverage + exon boundaries
 *
 * Processes: PDB_MAP, EXON_MAP, PDB_BULK_MAP
 * (split out of the former annotation_mapping.nf monolith)
 */



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
