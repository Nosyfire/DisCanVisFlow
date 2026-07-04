# Reference Data

Where every reference comes from, how current it is, and how to refresh the
local cache. For *how references are supplied to a run* (`--data` modes), see
[Installation § Reference data](installation.md#3-reference-data).

All downloads are cached under `references/` via Nextflow's `storeDir`, shared
across every project and run on the machine. The exact versions and entry counts
used by a completed run are recorded in
`results/<project>/mapping_reports/release.json` and the "Data source versions"
/ "Input scale" sections of `mapping_summary.md`.

---

## Sources & update cadence

| Source | Origin | Refresh | Freeze / currency |
|--------|--------|---------|-------------------|
| UniProt SwissProt | [ftp.uniprot.org](https://ftp.uniprot.org/pub/databases/uniprot/current_release/) | `bin/refresh_refs.sh uniprot` | Release captured in `release.json` |
| GENCODE | [gencodegenes.org](https://www.gencodegenes.org/human/) | `bin/refresh_refs.sh gencode` | Pinned to v44 by default |
| ClinVar | [NCBI ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/) | `bin/refresh_refs.sh clinvar` | Always-current via `FETCH_CLINVAR` |
| GO (GOA + OBO) | [geneontology.org](http://geneontology.org/) | `bin/refresh_refs.sh go` | Always-current via `FETCH_GO` |
| MobiDB | [mobidb.org](https://mobidb.org/) | `bin/refresh_refs.sh mobidb` | Always-current via `FETCH_MOBIDB` |
| ELM instances | [elm.eu.org](http://elm.eu.org/) | `legacy_data/elm/elm_instances-2023.tsv` | Frozen 2023 snapshot |
| dbSNP bigBed | [UCSC dbSnp155Common](https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/) | `bin/refresh_refs.sh dbsnp` (manual) | Large; rarely updated |
| AlphaMissense | [Zenodo 8208688](https://zenodo.org/records/8208688) | `bin/refresh_refs.sh alphamissense` | v2023 frozen |
| dbNSFP | [dbnsfp.org](https://www.dbnsfp.org/) | `--dbnsfp_raw_dir` or `--dbnsfp_tsv` | External; update manually |
| PPI (IntAct / BioGRID / HIPPIE) | [IntAct](https://www.ebi.ac.uk/intact/) · [BioGRID](https://thebiogrid.org/) · [HIPPIE](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/) | auto (`FETCH_INTACT/BIOGRID/HIPPIE` + `PPI_PREPROCESS`) | Auto on first run |

Licences and citation requirements for every source are in
[CITATIONS.md](../../CITATIONS.md).

---

## Managing the cache

List what is cached, and (re)generate the manifest:

```bash
bin/refresh_refs.sh                              # list cached references
python bin/generate_manifest.py --no_checksum    # writes references/MANIFEST.tsv
```

Force a re-download of specific sources, then `-resume` (only deleted files are
fetched again):

```bash
bin/refresh_refs.sh clinvar              # ClinVar only
bin/refresh_refs.sh clinvar mobidb go    # several
bin/refresh_refs.sh all                  # everything except hg38 / dbsnp / alphafold
bin/refresh_refs.sh --force all          # truly everything
```

If a cached file is 0 bytes (a failed download), delete it and re-run — see
[Troubleshooting § storeDir file is 0 bytes](troubleshooting.md#storedir-file-is-0-bytes-failed-download).
