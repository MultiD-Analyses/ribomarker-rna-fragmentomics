# buildDB: building the reference databases

These scripts build the references that the pipeline maps reads to: the Bowtie indexes, the BED file of features, and the list of mapping tiers in priority order (`<species>.dbpaths.csv`), together with a record of where every source came from (`MANIFEST_<species>.json`, `PROVENANCE_<species>.md`). The pipeline never runs them: it only uses finished references, built into the container image or given with `--references`.

Most users do not need this folder: the human and mouse references are built into the image. It is for rebuilding references from their sources, for building references of another species, and for reviewing how the published references were made.

## The published references

| Species | Built | Sources | Built with |
|---|---|---|---|
| human (GRCh38) | 2026-07-07 | Ensembl 116, GtRNAdb2 hg38, MirGeneDB 3.0, piRBase v3.0 (gold set), miRTrace rRNA, NCBI UniVec | an earlier version of `buildDB_v1.py` that was not archived (its checksum is in the MANIFEST) |
| mouse (GRCm39) | 2026-10-08 | Ensembl 116, GtRNAdb2 mm39, MirGeneDB 3.0, piRBase v3.0 (gold set), miRTrace rRNA, NCBI UniVec | `buildDB_v1.py` and `builder_core.py` as published, apart from comments |

Exact URLs, download times, tool versions and checksums: [`docker/assets/PROVENANCE_<species>.md`](../docker/assets/) and `MANIFEST_<species>.json`. Sources, citations and terms: [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).

## Rebuilding

```bash
conda activate ribomarker        # the environment of docker/env.lock.txt has every tool the build needs
cd buildDB
python buildDB_v1.py --species human --check-urls     # are all sources still online? (no large downloads)
python buildDB_v1.py --species human --inputs ./inputs --output ../my_references --threads 8
```

- The build downloads the genome and annotation and builds the indexes; it takes hours. The finished human and mouse folders hold about 8.6 GB together.
- [`inputs/`](inputs/README.md) holds the one file that cannot be downloaded: the RiboSpike spike-in sequences (`ribospike.fa`).
- Build each species into its own **empty** folder: the builder always rebuilds the indexes shared by all species (UniVec, RiboSpike) and writes intermediate files without a species prefix. `install_species.py` copies one species from such a folder into a folder that already holds another species, as `docker/assets` does.
- Tools used: `bowtie-build` (Bowtie 1, for UniVec and rRNA), `bowtie2-build` (all other tiers), `tar`, `gunzip`; Python packages `pandas`, `biopython`, `pysam`.

**A rebuild is not guaranteed to equal the published references.** The script downloads the files as they are today: NCBI UniVec has no version, the rRNA file comes from the `master` branch of miRTrace, and other sources may be updated at the same address. The new folder's MANIFEST records the size and checksum of every file it contains; compare them with `docker/assets/MANIFEST_<species>.json`: files with the same checksum are identical. Validate a changed build with a pipeline run before trusting its results.

`python update_manifest.py <folder> --species human --check` compares the files of a references folder with its own MANIFEST, e.g. after copying the folder.

## Using a built folder

Give the folder to the pipeline with `--references FOLDER`, with the starter, with `docker run` or without Docker; see the [user guide, *Own references*](../USER_GUIDE.md#own-references). The paths in `<species>.dbpaths.csv` are relative to the folder, so it can be moved or mounted anywhere.

Folder layout:

```
my_references/
  univec/ ribospike/ rrna/ trna/ mirna/ pirna/ gtf/ genome/   # the indexes
  human.bed                   # features: sequence name, start, end, feature, biotype, strand
  human.dbpaths.csv           # mapping tiers in priority order; index paths relative to this folder
  MANIFEST_human.json         # per source: URL, version, download time, sha256; tool versions;
                              # size and sha256 of the dbpaths file, the BED and every index file
  PROVENANCE_human.md         # the same, readable
```

## Things to know

- **Versions:** each build is tied to a numbered script. `buildDB_v1.py` stays as it is; a build with new sources goes into a copy, `buildDB_v2.py`, with the new URLs and versions in its `SOURCES`.
- **Ensembl's FTP is moving** to a new platform; release 116 is the last on the classic `ftp.ensembl.org`. If `--check-urls` reports the genome or annotation URLs as missing, switch them to the new host (`ftp.ebi.ac.uk/pub/ensemblorganisms/...`) in a new script version.
- **The rRNA file is pinned to a branch (`master`), not a commit.** The MANIFEST records the checksum of the exact file downloaded, so a build can be checked by content; for an exact URL, pin a commit.
- **Another species:** add an entry to `SOURCES` in a new script version, following the human and mouse entries, check every URL, and mark it `"verified"`; the builder refuses entries marked as placeholders. Run the pipeline with `--species NAME --references FOLDER`.
