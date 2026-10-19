# Changelog

All notable changes are listed here. Versions follow `MAJOR.MINOR.PATCH`; while the version is 0.x, a change that alters results (including a new reference build) raises the minor version, other changes the patch version. The git tag `v<version>` and the image tag `<version>` belong to the same commit.

## 0.1.0 – 2026-10-19

First public release.

- FASTQ → counts per gene, feature and RNA type, per-transcript coverage (`bedgraph/`, `bedgraph_full/`) and a technical QC report with consistency checks.
- One command, `run`, on a sample sheet (CSV: `sample`, `fastq_file`, optional `fastq_file2`, any further columns copied to the results); starter scripts for Linux, macOS and Windows.
- Container image with the conda environment of `docker/env.lock.txt` and the human and mouse references built in; `run --demo` checks an installation against synthetic reads with known results.
- Human and mouse references: see `docker/assets/PROVENANCE_human.md` and `PROVENANCE_mouse.md`. Mouse has been checked on synthetic reads only, not yet compared with results on real mouse data.
- Own references with `--references FOLDER`, also for other species (`--species NAME`); not tested beyond human, mouse and the demo references.
