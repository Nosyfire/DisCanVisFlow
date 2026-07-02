#!/usr/bin/env python3
"""Fetch somatic mutations from the cBioPortal public REST API.

Works without a specific study ID: queries all public mutation profiles
and returns mutations for the requested gene(s) as a MAF-like TSV.

Requires: requests (in discanvis conda env)
"""
import argparse
import csv
import sys
import time
import requests

BASE = "https://www.cbioportal.org/api"
CHUNK = 500  # max profile IDs per POST request


def get_entrez_ids(genes: list[str]) -> dict[str, int]:
    ids = {}
    for gene in genes:
        r = requests.get(f"{BASE}/genes/{gene}", timeout=30)
        if r.ok:
            ids[gene] = r.json()["entrezGeneId"]
        else:
            print(f"WARN: gene '{gene}' not found in cBioPortal ({r.status_code})", file=sys.stderr)
    return ids


def get_mutation_profile_ids() -> list[str]:
    r = requests.get(
        f"{BASE}/molecular-profiles",
        params={"projection": "ID", "molecularAlterationType": "MUTATION_EXTENDED"},
        timeout=120,
    )
    r.raise_for_status()
    return [p["molecularProfileId"] for p in r.json()]


def fetch_mutations(entrez_ids: list[int], profile_ids: list[str]) -> list[dict]:
    mutations = []
    total = len(profile_ids)
    for i in range(0, total, CHUNK):
        chunk = profile_ids[i : i + CHUNK]
        r = requests.post(
            f"{BASE}/mutations/fetch",
            json={"entrezGeneIds": entrez_ids, "molecularProfileIds": chunk},
            params={"projection": "SUMMARY"},
            timeout=120,
        )
        if r.ok:
            mutations.extend(r.json())
        else:
            print(f"WARN: chunk {i//CHUNK + 1} failed ({r.status_code})", file=sys.stderr)
        if i + CHUNK < total:
            time.sleep(0.15)
    return mutations


def write_maf(mutations: list[dict], entrez_map: dict[str, int], out_path: str) -> None:
    # reverse map: entrezGeneId → hugoSymbol
    id_to_symbol = {v: k for k, v in entrez_map.items()}
    cols = [
        "Hugo_Symbol",
        "Entrez_Gene_Id",
        "Variant_Classification",
        "Tumor_Sample_Barcode",
        "HGVSp_Short",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for m in mutations:
            eid = m.get("entrezGeneId")
            hugo = id_to_symbol.get(eid) or (m.get("gene") or {}).get("hugoGeneSymbol", "")
            pc = m.get("proteinChange", "")
            if pc and not pc.startswith("p."):
                pc = "p." + pc
            w.writerow({
                "Hugo_Symbol":          hugo,
                "Entrez_Gene_Id":       m.get("entrezGeneId", ""),
                "Variant_Classification": m.get("mutationType", ""),
                "Tumor_Sample_Barcode": m.get("sampleId", ""),
                "HGVSp_Short":          pc,
                "Chromosome":           m.get("chr", ""),
                "Start_Position":       m.get("startPosition", ""),
                "End_Position":         m.get("endPosition", ""),
                "Reference_Allele":     m.get("referenceAllele", ""),
                "Tumor_Seq_Allele2":    m.get("variantAllele", ""),
            })


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", help="Comma-separated HGNC gene symbols (e.g. RAF1,BRAF)")
    ap.add_argument("--gene_list_file", help="File with one HGNC symbol per line (# comments OK)")
    ap.add_argument("--out", required=True, help="Output MAF file path")
    args = ap.parse_args()

    if args.genes:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    elif args.gene_list_file:
        with open(args.gene_list_file) as fh:
            genes = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    else:
        ap.error("Provide --genes or --gene_list_file")

    print(f"Genes: {', '.join(genes)}", file=sys.stderr)

    entrez_map = get_entrez_ids(genes)
    if not entrez_map:
        sys.exit("ERROR: none of the requested genes found in cBioPortal")
    entrez_ids = list(entrez_map.values())
    print(f"Entrez IDs: {entrez_map}", file=sys.stderr)

    print("Fetching mutation profile IDs from all public studies...", file=sys.stderr)
    profile_ids = get_mutation_profile_ids()
    print(f"  {len(profile_ids)} mutation profiles found", file=sys.stderr)

    print("Fetching mutations (this may take a minute)...", file=sys.stderr)
    mutations = fetch_mutations(entrez_ids, profile_ids)
    print(f"  {len(mutations)} mutations retrieved", file=sys.stderr)

    write_maf(mutations, entrez_map, args.out)
    print(f"Written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
