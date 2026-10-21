# Changelog

All notable changes are listed here. Versions follow `MAJOR.MINOR.PATCH`; while the version is 0.x, a change that alters results (including a new reference build) raises the minor version, other changes the patch version. The git tag `v<version>` and the image tag `<version>` belong to the same commit.

## 0.2.0 – 2026-10-21

- Same results on every computer: overlapping features of equal rank are now ordered by a stable sort, so the feature listed first in the reference BED file wins on every computer. In version 0.1.0 the order of such ties depended on the computer's processor. Measured change against version 0.1.0: in two full samples 64 and 187 reads (0.005–0.006 %) are assigned to another feature, all on pre-mRNA intron features; the mapping statistics and the gene table are unchanged. Ties in output tables and top-100 lists (e.g. the read-sequence lists of single-end trimming) are now also in a fixed order.
- Documentation: what was tested for reproducibility, and on which systems (Windows on two computers, WSL2, macOS on an Intel Mac); the Docker Desktop licence terms; the harmless `x86-64-v3` warning.
- `buildDB/README.md` rewritten: which script version built the published references, how to rebuild them, and why a rebuild may differ from the published references.
- Zenodo DOI (all versions, 10.5281/zenodo.23475622) in the README and `CITATION.cff`.

## 0.1.0 – 2026-10-19

First public release.

- FASTQ → counts per gene, feature and RNA type, per-transcript coverage (`bedgraph/`, `bedgraph_full/`) and a technical QC report with consistency checks.
- One command, `run`, on a sample sheet (CSV: `sample`, `fastq_file`, optional `fastq_file2`, any further columns copied to the results); starter scripts for Linux, macOS and Windows.
- Container image with the conda environment of `docker/env.lock.txt` and the human and mouse references built in; `run --demo` checks an installation against synthetic reads with known results.
- Human and mouse references: see `docker/assets/PROVENANCE_human.md` and `PROVENANCE_mouse.md`. Mouse has been checked on synthetic reads only, not yet compared with results on real mouse data.
- Own references with `--references FOLDER`, also for other species (`--species NAME`); not tested beyond human, mouse and the demo references.
