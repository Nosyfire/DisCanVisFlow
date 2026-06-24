"""Tests for create_clinvar_disease_build_worker.py and clinvar_disease_lib.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("obonet")

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from clinvar_disease_lib import categorize_mondo_id, extract_mondo_ids, finalize_disease_row

WORKER = Path(__file__).parent.parent / "bin" / "create_clinvar_disease_build_worker.py"
MONDO_OBO = Path(
    "/dlab/home/norbi/PycharmProjects/DisCanVis_Data_Process/"
    "External_Data/gencode_process/disease_ontology/mondo/mondo.obo"
)


class TestClinvarDiseaseLib:
    def test_extract_mondo_ids(self):
        ids = extract_mondo_ids("MONDO:MONDO:0018997,MeSH:D009634")
        assert "MONDO:0018997" in ids

    def test_finalize_unknown_group(self):
        row = finalize_disease_row({"disease_group": "not provided", "Final_Category": "Other"})
        assert row["Final_Category"] == "Unknown"

    @pytest.mark.skipif(not MONDO_OBO.exists(), reason="MONDO OBO not available")
    def test_categorize_noonan(self):
        from clinvar_disease_lib import load_mondo_graph
        g = load_mondo_graph(str(MONDO_OBO))
        fc = categorize_mondo_id(g, "MONDO:0018997")
        assert fc in ("Cardiovascular/Hematopoietic", "Developmental", "Mixed", "Other", "Neurodegenerative")


class TestBuildWorkerCLI:
    @pytest.mark.skipif(not MONDO_OBO.exists(), reason="MONDO OBO not available")
    def test_build_from_mutations(self, tmp_path):
        loc = tmp_path / "loc.tsv"
        pd.DataFrame({"Protein_ID": ["RAF1-201"]}).to_csv(loc, sep="\t", index=False)

        mut_dir = tmp_path / "mutations"
        mut_dir.mkdir()
        pd.DataFrame({
            "Protein_ID": ["RAF1-201"],
            "PhenotypeList": ["Noonan_syndrome"],
            "PhenotypeIDS": ["MONDO:MONDO:0018997"],
        }).to_csv(mut_dir / "Missense_filter_mutations_mapped.tsv", sep="\t", index=False)

        out = tmp_path / "out"
        out.mkdir()
        r = subprocess.run(
            [sys.executable, str(WORKER),
             "--seq_table", str(loc),
             "--mondo_obo", str(MONDO_OBO),
             "--mutation_dir", str(mut_dir),
             "--outdir", str(out)],
            capture_output=True, text=True,
        )
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(out / "clinvar_disease.tsv", sep="\t")
        assert len(df) >= 1
        assert "Final_Category" in df.columns
