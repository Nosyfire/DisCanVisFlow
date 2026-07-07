"""One-off: build a tiny valid AlphaFold-style mmCIF fixture for DSSP tests.

Downloads the real AlphaFold model for human hemoglobin subunit alpha
(UniProt P69905, AF-P69905-F1, pLDDT very-high across nearly the whole chain)
and trims it to the first N residues (default 8: MVLSPADK) while keeping every
mmCIF category `mkdssp` 4.6.1 needs to parse the file (entity / entity_poly /
entity_poly_seq / struct_asym / pdbx_poly_seq_scheme / atom_site).

The ModelArchive-only `_ma_*` categories and the `_audit_conform` pointer to
`mmcif_ma.dic` are dropped entirely: this build of libcifpp does not ship
that dictionary and mkdssp segfaults trying to load it. `_struct_ref*` /
`_struct_conf*` (full-length reference-sequence and precomputed secondary
-structure metadata) are also dropped since they describe the untrimmed
142-residue chain and are not needed by mkdssp (it recomputes SS itself).

This script requires network access and is NOT run by the test suite — it was
run once to produce the committed tests/fixtures/dssp/AF-P0TEST-F1.cif, which
the tests read as a fully offline static fixture (renamed to the synthetic
accession P0TEST expected by test_create_dssp_worker.py).
"""
from pathlib import Path

import requests

SRC_URL = "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif"
N_RES = 8  # keep residues 1..8 (sequence MVLSPADK)
OUT = Path(__file__).with_name("AF-P0TEST-F1.cif")


def _split_blocks(lines):
    """Split mmCIF body into '#'-delimited blocks (this file uses '#' as a
    separator after every category), returning list of line-lists."""
    blocks, cur = [], []
    for ln in lines:
        if ln.rstrip() == "#":
            if cur:
                blocks.append(cur)
                cur = []
        else:
            cur.append(ln)
    if cur:
        blocks.append(cur)
    return blocks


def _category_of(block):
    for ln in block:
        if ln.startswith("_"):
            return ln.split(".", 1)[0]
    return None


def _is_loop(block):
    return block[0].strip() == "loop_"


def _loop_headers(block):
    return [ln for ln in block[1:] if ln.startswith("_")]


def _loop_rows(block, headers):
    return block[1 + len(headers):]


def _trim_loop_by_col(block, colnames):
    headers = _loop_headers(block)
    cols = [h.split(".", 1)[1] for h in headers]
    rows = _loop_rows(block, headers)
    seqcol = next((c for c in colnames if c in cols), None)
    if seqcol is None:
        return block
    idx = cols.index(seqcol)
    kept = []
    for r in rows:
        toks = r.split()
        if len(toks) <= idx:
            continue
        try:
            seqid = int(toks[idx])
        except ValueError:
            continue
        if seqid <= N_RES:
            kept.append(r)
    return ["loop_"] + headers + kept


def main():
    text = requests.get(SRC_URL, timeout=60).text
    blocks = _split_blocks(text.splitlines())

    KEEP_ASIS = {"_entry", "_atom_type", "_struct_asym", "_software", "_entity"}
    DROP = {
        "_audit_author", "_audit_conform", "_pdbx_database_status",
        "_pdbx_audit_revision_details", "_pdbx_audit_revision_history",
        "_pdbx_data_usage", "_database_2", "_struct_ref", "_struct_ref_seq",
        "_struct_conf", "_struct_conf_type",
    }

    out_blocks = []
    for b in blocks:
        cat = _category_of(b)
        if cat is None:
            continue
        if cat.startswith("_ma_") or cat in DROP:
            continue
        if cat == "_entity_poly":
            seq = "MVLSPADK"[:N_RES]
            # Rebuild explicitly rather than filtering the original multi-line
            # semicolon-quoted FASTA block (avoids fragile line-skipping logic).
            new_block = [
                "_entity_poly.entity_id                    1",
                "_entity_poly.nstd_linkage                 no",
                "_entity_poly.nstd_monomer                 no",
                "_entity_poly.pdbx_seq_one_letter_code",
                ";" + seq,
                ";",
                "_entity_poly.pdbx_seq_one_letter_code_can",
                ";" + seq,
                ";",
                "_entity_poly.pdbx_strand_id               A",
                "_entity_poly.type                         polypeptide(L)",
            ]
            out_blocks.append(new_block)
        elif cat in ("_entity_poly_seq",):
            out_blocks.append(_trim_loop_by_col(b, ["num", "seq_id"]))
        elif cat in ("_pdbx_poly_seq_scheme",):
            out_blocks.append(_trim_loop_by_col(b, ["seq_id"]))
        elif cat == "_atom_site":
            out_blocks.append(_trim_loop_by_col(b, ["label_seq_id"]))
        elif cat in KEEP_ASIS:
            out_blocks.append(b)
        else:
            continue  # drop anything else not explicitly whitelisted

    out_lines = ["data_AFP0TEST", "#"]
    for b in out_blocks:
        out_lines.extend(b)
        out_lines.append("#")

    OUT.write_text("\n".join(out_lines) + "\n")
    print(f"wrote {OUT} ({len(out_lines)} lines)")


if __name__ == "__main__":
    main()
