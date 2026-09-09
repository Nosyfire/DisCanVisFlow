---
name: discanvis-new-track
description: >
  Add a new annotation track to the DisCanVisFlow pipeline end to end — worker, tests,
  Nextflow process, wiring, config parameters, and every place the track has to be
  registered so it is discoverable and documented. Use this skill whenever someone wants to
  extend the pipeline with a new data source or predictor: "add a track for X", "integrate
  this database", "wire up a new predictor", "we should pull in FuzDrop", "add DUST/TRF
  complexity", "can we include this new LLPS dataset", or when reviewing whether an
  existing track was fully integrated. The value is the registration checklist: a track
  wired only into main.nf runs correctly but stays invisible in --modules listings, the
  docs index, the refresh tooling and the citation file, and past additions have each
  missed a different subset of those. Use it for pipeline development work, not for
  generating or interpreting data.
---

# Adding an annotation track

A track is easy to make *work* and easy to leave half-registered. The compute side —
worker, process, wiring — is what makes it run. The registration side is what makes it
findable by the next person, and it is the half that has actually been missed.

That is not hypothetical. Of the two most recent tracks, one updated `CITATIONS.md` and the
annotation index but not `CLAUDE.md`; the other did the reverse. Neither touched
`docs/guide/configuration.md`, so their module names were silently absent from the list
users are told is authoritative, and a wrong `--modules` name is accepted without
complaint and simply produces nothing.

Work through both halves. Follow the repo's own TDD protocol in `CLAUDE.md`: failing test
first, then implementation, then run the tests without being asked.

## Step 1 — Locate the checkout

```bash
python bin/find_discanvis.py
```

Development needs a real clone, not the plugin cache.

## Step 2 — Decide the shape of the track

Answer these before writing code, because they determine every later step:

- **Name.** One lowercase word, used identically as the `--modules` name, the `--skip_<name>`
  flag, and the worker filename. Divergence here is a lasting source of confusion.
- **Input.** A downloaded reference (needs a `FETCH_*` process and a `storeDir` cache), a
  local file supplied by a parameter, or pure sequence (needs neither).
- **Keying.** Per-residue scores, or regions with start/end? Keyed by UniProt accession
  (needs `TRANSCRIPT_MAP` to reach `Protein_ID`) or already per transcript?
- **Output location.** Which `final/` subdirectory, matching the existing categories.
- **Genome dependency.** Anything using genomic coordinates needs `combined_map.map` and
  therefore `params.hg38_2bit`, and must be skipped when genome mapping did not run.

## Step 3 — Build the compute half

**Worker** — `bin/create_<name>_worker.py`, a standalone `argparse` script. Every worker
here is callable on its own precisely so it can be tested without Nextflow. Match the
structure of a neighbouring worker rather than inventing one.

**Tests** — `tests/test_create_<name>_worker.py`, invoking the worker as a subprocess with
dummy input, as the existing tests do. Write them before the implementation.

**Process** — add to the `modules/*.nf` file matching the category (`structure.nf`,
`disorder.nf`, `functional.nf`, `pathogenicity.nf`, …), or a new file for a genuinely new
group. Fetches use `storeDir` so downloads survive across runs and projects:

```groovy
storeDir { workflow.stubRun ? "${params.ref_dir}/_stub/<name>" : "${params.ref_dir}/<name>" }
```

A `stub:` block that touches the output files is what lets `-stub` validate the graph
without computing anything — include one.

**Wiring** — `include` the process in `main.nf`, then gate it. The established pattern,
which keeps `--modules` and the skip flag independent:

```groovy
if ( (mods == null || mods.contains('<name>')) && !params.skip_<name> ) {
    <NAME>_MAP( SEQUENCE_PROCESS.out.loc_chrom_seq )
    <NAME>_MAP.out.<name>.view { f -> "\n✔  <Track name>: ${f}\n" }
}
```

Add the same condition to the `report_gate` mix near the end of `main.nf`, or the mapping
report runs before your track finishes and reports it as absent.

**Config** — declare parameters in `nextflow.config`: the input path (`null` when a
`FETCH_*` supplies it) and `skip_<name> = false`. Add the name to the module-name comment
in the same file.

## Step 4 — Register it, or it stays invisible

Each of these is where a *different* kind of user looks for the track. Skipping one does
not break the pipeline, which is exactly why they get skipped.

| File | What to add | Who this fails if missed |
|---|---|---|
| `docs/annotations/<category>/<name>.md` | New page: what it measures, source, output, columns, caveats | Anyone asking what a column means |
| `docs/annotations/README.md` | Row in the category table linking the page | Anyone browsing available tracks |
| `docs/guide/configuration.md` | The `--modules` name, and the `--skip_<name>` row | Anyone trying to run a subset |
| `docs/pipeline/architecture.md` | Module table row, and the output tree | Anyone learning the pipeline structure |
| `CLAUDE.md` | Row in the Module → File Mapping table, and the `--modules` list | Every future Claude session |
| `CITATIONS.md` | Citation **and licence** for the tool or database | Anyone publishing, or checking commercial use |
| `bin/refresh_refs.sh` | A `SOURCES` entry, if the track downloads a reference | Anyone refreshing that reference |
| `docs/guide/reference_data.md` | Source row, if it downloads a reference | Anyone auditing provenance |
| `README.md` | The category table, if it is a new kind of annotation | Anyone evaluating the pipeline |

The licence column in `CITATIONS.md` matters more than it looks. Several sources are
non-commercial or registration-gated, and `CITATIONS.md § 0` is the summary users rely on
to decide what they may run commercially. A missing entry there is a real problem, not a
documentation nicety.

## Step 5 — Verify

Run the tests, then prove the wiring works from the outside:

```bash
pytest tests/test_create_<name>_worker.py -v

# the module runs when asked for by name, on its own
nextflow run main.nf --project test_one_protein --machine laptop \
    --target_gene RAF1 --modules <name> -stub

# and is absent when not asked for
nextflow run main.nf --project test_one_protein --machine laptop \
    --target_gene RAF1 --modules disorder -stub
```

`-stub` exercises the graph in seconds without computing anything, so both directions are
cheap to check. Grepping the run output for your process name is the actual evidence that
the gate works — a misspelled module name fails silently and looks exactly like a track
that legitimately did not run.

Finally, confirm the registration actually landed:

```bash
grep -rn "<name>" docs/guide/configuration.md docs/annotations/README.md CLAUDE.md CITATIONS.md
```

If any of those comes back empty, the track is not finished.
