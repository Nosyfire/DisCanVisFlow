#!/usr/bin/env python3
"""
Interactive project config generator for the DisCanVis Nextflow pipeline.

Asks a series of questions and writes a ready-to-use projects/<name>.yaml
that can be passed directly to the pipeline:

    nextflow run main.nf --project projects/<name>.yaml -profile reproducible,conda -resume

Usage:
    python bin/new_project.py
    python bin/new_project.py --output projects/my_study.yaml   # skip the filename question
    python bin/new_project.py --non-interactive projects/my_study.yaml  # defaults only (for CI)
"""

import argparse
import sys
from pathlib import Path
from textwrap import dedent


# ── helpers ───────────────────────────────────────────────────────────────────

def ask(prompt, default=None, choices=None):
    suffix = f" [{default}]" if default is not None else ""
    if choices:
        suffix += f" ({'/'.join(choices)})"
    while True:
        raw = input(f"{prompt}{suffix}: ").strip()
        if not raw and default is not None:
            return default
        if choices and raw.lower() not in [c.lower() for c in choices]:
            print(f"  Please enter one of: {', '.join(choices)}")
            continue
        return raw


def ask_bool(prompt, default=True):
    d = "y" if default else "n"
    ans = ask(prompt, default=d, choices=["y", "n", "yes", "no"])
    return ans.lower() in ("y", "yes")


def ask_multiselect(prompt, options, defaults=None):
    """Show numbered options, let user enter comma-separated numbers or 'all'/'none'."""
    defaults = defaults or []
    print(f"\n{prompt}")
    for i, (key, desc) in enumerate(options, 1):
        mark = "*" if key in defaults else " "
        print(f"  [{mark}] {i:2d}. {key:25s}  {desc}")
    print("  Enter numbers (e.g. 1,3,5), 'all', 'none', or press Enter for defaults [*]")
    raw = input("  > ").strip()
    if not raw:
        return list(defaults)
    if raw.lower() == "all":
        return [k for k, _ in options]
    if raw.lower() == "none":
        return []
    selected = []
    for part in raw.split(","):
        part = part.strip()
        if part.isdigit():
            idx = int(part) - 1
            if 0 <= idx < len(options):
                selected.append(options[idx][0])
            else:
                print(f"  Warning: {part} out of range, skipped")
        else:
            print(f"  Warning: '{part}' not a number, skipped")
    return selected


# ── track catalogue ───────────────────────────────────────────────────────────

ALL_TRACKS = [
    ("mutations",       "Missense/frameshift/indel/nonsense mutations (ClinVar/MAF/VCF)"),
    ("clinvar_disease", "ClinVar disease ontology (MONDO) + disease-mutation links"),
    ("disorder",        "Intrinsic disorder: IUPred3 + ANCHOR2 + AIUPred + AlphaFold pLDDT"),
    ("slim_annotations","SLiMs & PTMs: ELM + DIBS + MFIB + PhasePro + PTM + Pfam + UniProt ROI"),
    ("pdb",             "PDB 3D structures (coverage, resolution, chain regions)"),
    ("go_terms",        "Gene Ontology (GO) annotations from GOA"),
    ("polymorphism",    "Common SNPs from dbSNP 155 with allele frequencies"),
    ("pathogenicity",   "Pathogenicity scores (dbNSFP pre-mapped table)"),
    ("alphamissense",   "AlphaMissense per-residue pathogenicity (Google/DeepMind)"),
    ("mavedb",          "MaveDB single-mutant functional scores (DMS assays)"),
    ("proteingym",      "ProteinGym deep mutational scanning scores"),
    ("cancer_drivers",  "Cancer driver classification (CGC census + Compendium)"),
    ("disease",         "OMIM disease ontology + OMIM mutations (requires OMIM key)"),
    ("depmap",          "DepMap cancer cell line mutations"),
    ("ppi",             "Protein–protein interactions (IntAct + BioGRID + HIPPIE)"),
    ("scansite",        "Phosphorylation and kinase motifs (MIT ScanSite 4.0)"),
    ("pem",             "PEM core motifs + isoform transfer"),
    ("conservation",    "Conservation scores (GOPHER multi-level, requires local table)"),
    ("coiled_coils",    "Coiled-coil predictions (DeepCoil — requires TF 2.9 / CUDA 11)"),
    ("finches",         "FINCHES Δε LLPS-change scores (CC BY-NC 4.0, off by default)"),
]

DEFAULT_TRACKS = {
    "full":        [k for k, _ in ALL_TRACKS if k not in ("conservation","coiled_coils","finches")],
    "mutations":   ["mutations", "clinvar_disease", "disease"],
    "disorder":    ["disorder", "slim_annotations", "pdb"],
    "pathogenicity": ["mutations", "pathogenicity", "alphamissense", "mavedb", "proteingym"],
    "minimal":     ["mutations", "disorder", "slim_annotations"],
}

PRESET_DESCS = {
    "full":         "All tracks (full DisCanVis update)",
    "mutations":    "Mutations + disease ontology only",
    "disorder":     "Disorder + SLiMs + PDB only",
    "pathogenicity":"Mutations + all pathogenicity predictors (VEP benchmarking)",
    "minimal":      "Mutations + disorder + SLiMs (quick single-protein run)",
    "custom":       "Let me choose individual tracks",
}

# ── gene-list format helpers ──────────────────────────────────────────────────

GENE_LIST_FORMATS = [
    ("gene_name",     "HGNC gene names (e.g. RAF1, TP53, BRCA1)"),
    ("uniprot_id",    "UniProt accession (e.g. P04049, P04637)"),
    ("uniprot_name",  "UniProt entry name (e.g. RAF1_HUMAN, P53_HUMAN)"),
    ("gencode_id",    "GENCODE gene ID (e.g. ENSG00000132155)"),
    ("transcript_id", "GENCODE transcript ID (e.g. ENST00000251849)"),
    ("protein_id",    "GENCODE transcript name (e.g. RAF1-201)"),
]

MUTATION_SOURCES = [
    ("ClinVar",    "NCBI ClinVar (auto-downloaded) — pathogenic/likely-pathogenic variants"),
    ("TCGA",       "TCGA MAF file (somatic mutations from cancer cohorts)"),
    ("CBioportal", "cBioPortal MAF file"),
    ("CustomVCF",  "Any VCF file (non-ClinVar)"),
]

# ── YAML writer ───────────────────────────────────────────────────────────────

def tracks_to_yaml(selected):
    lines = []
    for key, _ in ALL_TRACKS:
        val = "true" if key in selected else "false"
        lines.append(f"  {key}: {val}")
    return "\n".join(lines)


def write_config(cfg, out_path):
    target_line = f'target_genes: "{cfg["target_genes"]}"' if cfg["target_genes"] else "target_genes: null  # all SwissProt human proteins"

    mut_extra = ""
    if cfg["mutation_source"] == "ClinVar":
        mut_extra = "  # override with clinvar_vcf: /path/to/clinvar.vcf.gz"
    elif cfg["mutation_source"] in ("TCGA", "CBioportal"):
        mut_extra = f'\nmutation_maf: ""     # required: path to .maf file'
    elif cfg["mutation_source"] == "CustomVCF":
        mut_extra = f'\nmutation_vcf: ""     # required: path to .vcf.gz file'

    scatter = 1 if cfg["single_gene"] else 20
    blat    = 1 if cfg["single_gene"] else 16

    text = dedent(f"""\
        # DisCanVis Nextflow — project config
        # Generated by bin/new_project.py
        #
        # Run:
        #   nextflow run main.nf --project {out_path.name} -profile reproducible,conda -resume
        #
        # Override any param on the command line, e.g.:
        #   nextflow run main.nf --project {out_path.name} --target_gene TP53 -profile reproducible,conda -resume

        name: {cfg['name']}
        description: "{cfg['description']}"
        version: "{cfg['version']}"

        # ── Gene scope ──────────────────────────────────────────────────────────────
        {target_line}
        mapping_mode: {cfg['mapping_mode']}
        # To use a gene list file instead of --target_gene, set:
        # gene_list_file: "gene_lists/my_genes.txt"   (one gene per line, HGNC names)

        # ── Tracks ──────────────────────────────────────────────────────────────────
        tracks:
        {tracks_to_yaml(cfg['tracks'])}

        # ── Mutation source ──────────────────────────────────────────────────────────
        mutation_source: {cfg['mutation_source']}{mut_extra}

        # ── Compute ─────────────────────────────────────────────────────────────────
        scatter_chunks: {scatter}     # parallel DISORDER_MAP chunks (increase for full proteome on GPU)
        blat_chunks: {blat}

        # ── Output ──────────────────────────────────────────────────────────────────
        outdir: "{cfg['outdir']}"
        per_gene_md_threshold: 50    # write per-gene MD reports only for runs ≤ 50 genes
    """)
    out_path.write_text(text)
    print(f"\n✔  Written to: {out_path}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Interactive DisCanVis project config generator")
    ap.add_argument("--output", default=None, help="Output YAML path (skips filename question)")
    ap.add_argument("--non-interactive", action="store_true", help="Use defaults silently")
    args = ap.parse_args()

    ni = args.non_interactive

    print("\n" + "═" * 60)
    print("  DisCanVis Nextflow — new project config generator")
    print("═" * 60 + "\n")

    # ── basic metadata ────────────────────────────────────────────────────────
    name    = "my_project" if ni else ask("Project name (used as output dir suffix)", "my_project")
    desc    = ""           if ni else ask("Short description", "")
    version = "2026-q4"   if ni else ask("Version tag", "2026-q4")

    # ── gene scope ────────────────────────────────────────────────────────────
    print("\n── Gene scope ───────────────────────────────────────────")
    scope_choice = "full" if ni else ask(
        "Scope", default="full",
        choices=["full", "subset", "single"])

    target_genes = None
    single_gene = False

    if scope_choice == "full":
        target_genes = None
        print("  → Full human proteome (all 19,627 SwissProt proteins)")
    elif scope_choice == "single":
        single_gene = True
        if ni:
            target_genes = "RAF1"
        else:
            target_genes = ask("Gene name (HGNC)", "RAF1")
    else:
        if ni:
            target_genes = "TP53,BRCA1"
        else:
            print("\n  Gene list format:")
            fmt_choice = ask_multiselect("Which format is your gene list?",
                                         GENE_LIST_FORMATS, defaults=["gene_name"])
            fmt = fmt_choice[0] if fmt_choice else "gene_name"
            print(f"\n  Enter comma-separated values ({fmt}) or a file path (one per line):")
            raw = input("  > ").strip()
            if raw and Path(raw).exists():
                genes = [l.strip() for l in Path(raw).read_text().splitlines() if l.strip()]
                target_genes = ",".join(genes)
                print(f"  → Loaded {len(genes)} genes from {raw}")
            else:
                target_genes = raw
                print(f"  → {len(target_genes.split(','))} gene(s) selected")

    mapping_mode = "all_isoform_mapping" if ni else ask(
        "Mapping mode", default="all_isoform_mapping",
        choices=["all_isoform_mapping", "main_isoform_mapping"])

    # ── track selection ────────────────────────────────────────────────────────
    print("\n── Track selection ──────────────────────────────────────")
    print("  Presets:")
    for k, d in PRESET_DESCS.items():
        print(f"    {k:15s}  {d}")

    preset = "minimal" if (ni and single_gene) else ("full" if ni else
              ask("Preset", default="minimal" if single_gene else "full",
                  choices=list(PRESET_DESCS.keys())))

    if preset == "custom":
        default_sel = DEFAULT_TRACKS["full"] if not single_gene else DEFAULT_TRACKS["minimal"]
        selected_tracks = ask_multiselect(
            "Select tracks to enable:", ALL_TRACKS, defaults=default_sel)
    else:
        selected_tracks = list(DEFAULT_TRACKS[preset])
        print(f"  → {len(selected_tracks)} tracks enabled from preset '{preset}'")

    # ── mutation source ───────────────────────────────────────────────────────
    print("\n── Mutation source ──────────────────────────────────────")
    if "mutations" not in selected_tracks:
        mut_source = "ClinVar"
        print("  (mutations track disabled — source ignored)")
    else:
        if ni:
            mut_source = "ClinVar"
        else:
            print("  Options:")
            for i, (k, d) in enumerate(MUTATION_SOURCES, 1):
                print(f"    {i}. {k:12s}  {d}")
            raw = ask("Choice", default="1")
            idx = (int(raw) - 1) if raw.isdigit() else 0
            mut_source = MUTATION_SOURCES[min(idx, len(MUTATION_SOURCES)-1)][0]

    # ── output ────────────────────────────────────────────────────────────────
    outdir = f"results/{name}" if ni else ask("Output directory", f"results/{name}")

    # ── write ─────────────────────────────────────────────────────────────────
    projects_dir = Path(__file__).parent.parent / "projects"
    projects_dir.mkdir(exist_ok=True)

    if args.output:
        out_path = Path(args.output)
    else:
        default_fname = f"{name}.yaml"
        fname = default_fname if ni else ask(
            "Save config as", f"projects/{default_fname}")
        out_path = Path(fname) if fname.startswith("projects/") or "/" in fname \
                   else projects_dir / fname

    cfg = {
        "name":          name,
        "description":   desc,
        "version":       version,
        "target_genes":  target_genes,
        "single_gene":   single_gene,
        "mapping_mode":  mapping_mode,
        "tracks":        selected_tracks,
        "mutation_source": mut_source,
        "outdir":        outdir,
    }
    write_config(cfg, out_path)

    # ── run hint ──────────────────────────────────────────────────────────────
    print(f"""
Run command:
  conda activate discanvis
  nextflow run main.nf --project {out_path} \\
      --fetch_uniprot_dat --fetch_interpro_pfam --fetch_alphafold_bulk \\
      --skip_coiledcoils true \\
      -profile reproducible,conda -resume
""")


if __name__ == "__main__":
    main()
