"""
MONDO OBO disease categorization and ClinVar disease finalization rules.

Ported (simplified) from DisCanVis_Data_Process Complete_Ontology scripts and
IDP finalize_ontology_rule.py.
"""

from __future__ import annotations

import re
from collections import deque

try:
    import obonet
except ImportError:
    obonet = None

# MONDO ancestor name → organ-system categories (from 1_Categorize_ontologies.py)
CATEGORY_DICTIONARY_MONDO: dict[str, list[str]] = {
    "cardiovascular disorder": ["Cardiovascular/Hematopoietic"],
    "cardiogenetic disease": ["Cardiovascular/Hematopoietic"],
    "inherited hemoglobinopathy": ["Cardiovascular/Hematopoietic"],
    "inherited blood coagulation disorder": ["Cardiovascular/Hematopoietic"],
    "cancer or benign tumor": ["Cancer"],
    "prostate cancer, hereditary": ["Cancer"],
    "hereditary neoplastic syndrome": ["Cancer"],
    "endocrine system disorder": ["Endocrine"],
    "auditory system disorder": ["Neurodegenerative"],
    "inherited auditory system disease": ["Neurodegenerative"],
    "nervous system disorder": ["Neurodegenerative"],
    "disorder of visual system": ["Neurodegenerative"],
    "hereditary neurological disease": ["Neurodegenerative"],
    "hereditary dementia": ["Neurodegenerative"],
    "cognitive disorder": ["Neurodegenerative"],
    "syndromic intellectual disability": ["Neurodegenerative", "Developmental"],
    "complex neurodevelopmental disorder": ["Neurodegenerative", "Developmental"],
    "intellectual disability": ["Neurodegenerative", "Developmental"],
    "hematologic disorder": ["Cardiovascular/Hematopoietic"],
    "digestive system disorder": ["Gastrointestinal"],
    "immune system disorder": ["Immune"],
    "immune deficiency disease": ["Immune"],
    "immunodeficiency disease": ["Immune"],
    "immunodeficiency": ["Immune"],
    "integumentary system disorder": ["Integumentary"],
    "metabolic disease": ["Metabolic"],
    "musculoskeletal system disorder": ["Musculoskeletal"],
    "hereditary skeletal muscle disorder": ["Musculoskeletal"],
    "osteogenesis imperfecta": ["Musculoskeletal"],
    "skeletal dysplasia": ["Musculoskeletal"],
    "reproductive system disorder": ["Reproductive"],
    "respiratory system disorder": ["Respiratory"],
    "urinary system disorder": ["Urinary"],
    "disorder of development or morphogenesis": ["Developmental"],
    "developmental disorder of mental health": ["Developmental"],
}

DESIRED_CATEGORIES = [
    "Cancer", "Cardiovascular/Hematopoietic", "Developmental", "Endocrine",
    "Gastrointestinal", "Immune", "Integumentary", "Metabolic", "Musculoskeletal",
    "Neurodegenerative", "Reproductive", "Respiratory", "Urinary",
]

MONDO_ID_RE = re.compile(r"MONDO:(\d+)")


def load_mondo_graph(mondo_obo_path: str):
    if obonet is None:
        raise ImportError("obonet is required for MONDO disease categorization")
    return obonet.read_obo(mondo_obo_path)


def _ancestor_names(graph, node_id: str) -> set[str]:
    if node_id not in graph:
        return set()
    names: set[str] = set()
    queue = deque([node_id])
    seen = {node_id}
    while queue:
        cur = queue.popleft()
        name = graph.nodes[cur].get("name", "")
        if name:
            names.add(name.lower())
        for parent in graph.predecessors(cur):
            if parent not in seen:
                seen.add(parent)
                queue.append(parent)
    return names


def categorize_mondo_id(graph, mondo_id: str) -> str:
    """Return Final_Category for a MONDO ID using ancestor walk."""
    if not mondo_id.startswith("MONDO:"):
        mondo_id = f"MONDO:{mondo_id}"
    if mondo_id not in graph:
        return "Other"

    anc_names = _ancestor_names(graph, mondo_id)
    matched: set[str] = set()
    for anc in anc_names:
        for cat in CATEGORY_DICTIONARY_MONDO.get(anc, []):
            matched.add(cat)

    if not matched:
        return "Other"
    if "Cancer" in matched:
        return "Cancer"
    organ = [c for c in matched if c in DESIRED_CATEGORIES]
    if len(organ) == 1:
        return organ[0]
    if len(organ) > 1:
        return "Mixed"
    return "Other"


def extract_mondo_ids(phenotype_ids: str) -> list[str]:
    if not phenotype_ids or str(phenotype_ids) in ("nan", "-"):
        return []
    ids = []
    for m in MONDO_ID_RE.finditer(str(phenotype_ids)):
        ids.append(f"MONDO:{m.group(1)}")
    return ids


def finalize_disease_row(row: dict) -> dict:
    """Apply IDP/DisCanVis finalize rules on a disease summary row."""
    dg = str(row.get("disease_group", row.get("Disease", ""))).strip()
    fc = str(row.get("Final_Category", "Other")).strip()
    if fc in ("", "-", "nan"):
        fc = "Other"

    if dg.lower() in {"not provided", "not specified", "-", "nan", ""}:
        fc = "Unknown"
    elif dg == "Inborn genetic diseases":
        fc = "Inborn Genetic Diseases"

    row = dict(row)
    row["disease_group"] = dg
    row["Final_Category"] = _categorize_not_found(row, fc)
    if row.get("Developmental") in (True, "True", "true", "1"):
        if row["Final_Category"] == "Neurodegenerative":
            row["Final_Category"] = "Neurodevelopmental"
    return row


def _categorize_not_found(row: dict, fc: str) -> str:
    """Keyword heuristics for Other/Mixed rows (subset of finalize_ontology.py)."""
    if fc not in ("Other", "Mixed"):
        return fc
    dg = str(row.get("disease_group", "")).lower()
    checks = [
        (["cancer", "tumor", "malignant", "carcinoma", "sarcoma"], "Cancer"),
        (["cardiomyopathy", "thoracic aortic aneurysm", "long qt"], "Cardiovascular/Hematopoietic"),
        (["spastic paraplegia", " myopathy", "muscular dystrophy"], "Musculoskeletal"),
        (["retinitis pigmentosa", "deafness", "encephalopathy", "neurodegeneration"], "Neurodegenerative"),
        (["kidney disease", "renal disease"], "Urinary"),
        (["mitochondrial complex", "mucopolysaccharidosis"], "Metabolic"),
    ]
    for keywords, cat in checks:
        if any(k in dg for k in keywords):
            return cat
    return fc
