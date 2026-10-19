# buildDB — versioned, provenance-first reference builder

Standalone tooling to build the mapping reference databases (`assets_v?`) for the
fragmentomics pipeline. It is **deliberately independent of the Snakemake
pipeline**: the pipeline only ever *consumes* prebuilt assets; this folder is how
those assets are *made*, with a full record of when, from where, and which
version each database came from.

## Why versioned

Each build is tied to a numbered script:

- `buildDB_v1.py` → produces `assets_v1/`
- to refresh later: copy to `buildDB_v2.py`, bump the URLs/versions in `SOURCES`,
  run it → `assets_v2/`. Old versions stay frozen and auditable.

## Layout

```
buildDB/
  buildDB_v1.py     # driver: pinned SOURCES + downloads + provenance (edit per version)
  builder_core.py   # the index-building steps (kept format-compatible with the pipeline)
  inputs/           # non-downloadable inputs you provide (see inputs/README.md)
  README.md         # this file
```

Output (created by running the script, at repo root, next to `assets/`):

```
assets_v1/
  univec/ ribospike/ rrna/ trna/ mirna/ pirna/ gtf/ genome/   # the indices
  human.bed
  human.dbpaths.csv           # mapping tiers; index paths relative to this folder (local and Docker)
  MANIFEST_human.json         # per-DB: url, version, timestamp, sha256, tool versions;
                              # sha256 and size of the dbpaths file, the BED and every index file
  PROVENANCE.md               # human-readable manifest
```

## Quickstart

```bash
cd buildDB

# 1) confirm every source URL is alive and show the pinned versions (no big downloads)
python buildDB_v1.py --species human --check-urls

# 2) put ribospike.fa in ./inputs (see inputs/README.md), then build
python buildDB_v1.py --species human --inputs ./inputs --output ../assets_v1 --threads 8
```

**Dependencies** — on PATH: `bowtie-build` (bowtie v1), `bowtie2-build`, `tar`,
`gunzip`. Python: `pandas`, `biopython`, `pysam`. (The `ribomarker` conda env
already has these.)

## Using assets_v1 with the pipeline

`assets_v1` is a drop-in for the existing `assets/`:

- **Local:** set `dbdir: ".../assets_v1"` in your local config.
- **Docker:** pass the folder to the run command with `--references` (it is mounted read-only;
  the relative `human.dbpaths.csv` works in both places).
- **Folders built before 2026-10-08** have `human.local.dbpaths.csv` (absolute paths) and
  `human.docker.dbpaths.csv` instead. Bring them up to date without rebuilding:
  `python buildDB/update_manifest.py <folder> --species human` (writes `human.dbpaths.csv` and the file
  checksums in the MANIFEST; `--check` later compares the files with it). Locally the older file is still
  accepted, with a warning; in Docker it is refused.

## Rebuilding one species in a shared folder (e.g. mouse in `docker/assets`)

Never build a second species into a folder that already holds one: `buildDB_v1.py` always rebuilds the
shared UniVec and RiboSpike indexes and writes unprefixed intermediates (`trna/trnas.fa`,
`trna/pretrnas.bed`, `rrna/rrna.species.fa`). Build into an empty folder and install from there:

```bash
python buildDB/buildDB_v1.py --species mouse --output ../assets_mouse_build     # empty folder, hours
python buildDB/install_species.py ../assets_mouse_build docker/assets --species mouse
python test/levels/check_mouse_references.py docker/assets                     # content + both MANIFESTs
bash test/levels/run_mouse_synthetic.sh                                         # synthetic mouse reads
```

`install_species.py` copies only the species' files (`mouse.*`, `mouse-*`, `mm39-*`, MANIFEST, PROVENANCE),
reports whether the new build's shared indexes equal the installed ones (the installed ones are kept),
moves the replaced files into `<build>/_replaced_mouse_<time>/` (nothing is deleted) and recomputes the
MANIFEST's file checksums from the installed files.

## Sources pinned for v1 (resolved 2026-07-02)

| DB | Version | Notes |
|----|---------|-------|
| genome / gtf | Ensembl **release-116** (GRCh38) | liveness-checked ✓ |
| tRNA | GtRNAdb2 **hg38** tarball | switched from a bundled tarball to a real URL |
| miRNA | **MirGeneDB 3.0** | |
| piRNA | **piRBase v3.0** gold set | newest; upstream MD5 recorded in the manifest |
| rRNA | mirtrace `rRNA_reference.fa` (**master**) | content is hashed in the manifest |
| UniVec | NCBI UniVec (current build) | |
| ribospike | local `inputs/ribospike.fa` | synthetic; not downloadable |

### Things to know (flagged honestly)

- **Ensembl FTP migration:** Ensembl is moving to a new platform/FTP in summer
  2026; 116 is the last release on the classic `ftp.ensembl.org`. Run
  `--check-urls` before a build; if the genome/gtf URLs 404, switch them to the
  new host (`ftp.ebi.ac.uk/pub/ensemblorganisms/...`) in `SOURCES`.
- **rRNA is pinned to a moving branch (`master`)**, not a commit. The manifest
  records the sha256 of the exact bytes downloaded, so a build is still
  reproducible by content; for strict URL reproducibility, pin a commit SHA.
- **Index tools match production, not the old script.** `univec` and `rrna` are
  built as **bowtie v1** indices (`univec_bt1`, `human-rrna_bt1`), everything
  else bowtie2 — verified against the current `assets/` `dbpaths`. (The old
  `buildDB.py` built univec/rrna with bowtie2 while tagging them `bowtie`, which
  would not have matched the pipeline.)
- **`mouse` is a placeholder** (`status: PLACEHOLDER-unverified`). Its URLs
  follow the human pattern but are unchecked — verify each (especially the
  GtRNAdb genome code `Mmusc39/mm39`) and set `status: "verified"` before use.
  The builder refuses to run a placeholder species.

## Important: this is a NEWER build than `assets/`

`assets_v1` uses **newer source versions** than the current `assets/` (Ensembl
116 vs 113, current MirGeneDB/piRBase, etc.). It is therefore **not byte-identical**
to `assets/`, and read counts produced with it will differ — that is the point.
What is guaranteed is **structural/format compatibility** (same index tools,
names, padding, BED6 columns, dbpaths layout). Before trusting `assets_v1`,
validate with a real pipeline run and sanity-check the outputs.
