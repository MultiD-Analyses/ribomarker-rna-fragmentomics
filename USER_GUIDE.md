# User guide

This guide describes the RiboMarker pipeline in full: installation, inputs, running, outputs, the QC report, how the pipeline works and its known limitations. For a short start, see the [README](README.md).

In the commands below, `IMAGE` stands for the pipeline's container image, `dockermultid/ribomarker-rna-fragmentomics:<version>` on Docker Hub; the version of each release is in [CHANGELOG.md](CHANGELOG.md).

## Contents

1. [What the pipeline does](#1-what-the-pipeline-does)
2. [Installation](#2-installation)
3. [Inputs](#3-inputs)
4. [Running](#4-running)
5. [Advanced settings](#5-advanced-settings)
6. [Outputs](#6-outputs)
7. [The QC report](#7-the-qc-report)
8. [How the pipeline works](#8-how-the-pipeline-works)
9. [Known limitations](#9-known-limitations)
10. [Running without Docker](#10-running-without-docker)
11. [Building the image](#11-building-the-image)
12. [Support, citing and licence](#12-support-citing-and-licence)
13. [Glossary](#13-glossary)

---

## 1. What the pipeline does

The pipeline processes small-RNA sequencing reads (single- or paired-end FASTQ files) from libraries prepared with the RiboMarker protocol of RealSeq Biosciences ([Clark et al. 2026](https://doi.org/10.1261/rna.080861.125)), for example from cell-free RNA in plasma. It writes:

- **counts** per gene, per reference feature and per RNA type (biotype);
- **per-transcript coverage** files (bedgraph): how many reads cover each position of each transcript;
- a **technical QC report**: every step's numbers, with consistency checks between the steps.

The pipeline stops there. It applies **no normalization** and no statistics; the gene counts are integer counts, suitable for count-based differential-expression tools.

Scope of this version:

- **Library protocol:** RiboMarker (RealSeq) only. The adapter sequences and the removal of one extra base at the insert's 5′ end are fixed for this protocol (§8).
- **Species:** human (default) and mouse (`--species mouse`). The references of both are built into the image (§8, *References*); other species with your own references (§5, *Own references*). For mouse, the mouse references were checked on synthetic reads cut from them (every read lands on its own reference), but no result on real mouse data has been compared yet; the QC report of a mouse run says so (§9).
- **Systems:** the image is built for x86-64 (Intel/AMD). It has been tested on Windows (PowerShell and WSL2) and on Linux under WSL2; not yet on a native Linux installation. macOS has not been tested; on Apple-silicon Macs Docker runs x86-64 images under emulation, which is slower.

## 2. Installation

### 2.1 Requirements

| | |
|---|---|
| Docker | Docker Desktop (Windows, macOS) or Docker Engine (Linux), installed and running |
| Memory | 12 GB available to Docker. The peak grows with sample size: with 8 threads, a paired-end run of 2 samples (0.30 GB of R1 files) peaked at 4.1 GB, one of 6 samples (1.72 GB, the largest sample 0.67 GB) at 8.1 GB. `run` warns when less than 12 GB is available. Docker Desktop gives containers part of the computer's memory: check *Settings → Resources*. |
| Disk | about 21 GB for the image, as Docker Desktop reports it (human and mouse references); per run about 3 × the size of the FASTQ files + 1 GB (checked before the run starts) |
| Time | about 2–3 hours per GB of R1 FASTQ files with 8 threads. Two measurements on one laptop: 37–41 min for 2 samples (0.30 GB of R1 files), and 4 h 41 min for 6 samples (1.72 GB) through the Windows starter. It depends on the computer and the samples; samples are processed one after another. |

### 2.2 Docker on each system

- **Linux:** install Docker Engine (or Docker Desktop) and start it. Your user must be allowed to run Docker (the `docker` group, or `sudo`).
- **Windows:** install Docker Desktop with the **WSL2 backend**. You can then use either **PowerShell** with Windows paths (`C:\data\...`), or a **WSL shell** (e.g. Ubuntu from the Microsoft Store, with Docker's WSL integration switched on under *Settings → Resources → WSL integration*) with Linux paths (`/mnt/c/data/...`). In the WSL shell, if commands fail with "No such file or directory" after a run, `cd` into your folder again by its full path.
- **macOS (not tested):** install and start Docker Desktop and use Terminal with paths such as `/Users/you/data`. Docker Desktop shares your home folder by default; other locations may need adding under *Settings → Resources → File sharing*.

### 2.3 The starter

The starter is a small script that writes the `docker run` command from your own paths. Download the repository (*Code → Download ZIP*, or a release's source archive) and unzip it; the starters are at the top level:

- `ribomarker` (bash): Linux, macOS, Windows WSL shell. If `./ribomarker` says *permission denied*, run `chmod +x ribomarker` once.
- `ribomarker.cmd` with `ribomarker_windows.ps1`: Windows PowerShell or command prompt. Keep the two files together; the `.cmd` file starts the PowerShell script without changing PowerShell's script setting (execution policy).

Run the commands from that folder (`./ribomarker`, `.\ribomarker`), or add the folder to your `PATH` and type `ribomarker`. Each starter pins the image version of its release; Docker downloads the image the first time it is used. To use another image, set the environment variable `RIBOMARKER_IMAGE`.

### 2.4 Test the installation

```bash
./ribomarker --demo results/demo            # Linux, macOS, Windows WSL shell
```
```powershell
.\ribomarker --demo results\demo            # Windows PowerShell
```

The demo runs synthetic paired-end reads whose results are known at every step (built-in, a few minutes) and compares every known number. It ends with `DEMO PASSED: all <n> known results reproduced`, or lists what differs. The results are in `results/demo/demo/`.

## 3. Inputs

The input is a **sample sheet**: a CSV file in the folder of the FASTQ files, with one row per sample.

```
sample,fastq_file,fastq_file2,age,group
Healthy_1,H1_R1.fastq.gz,H1_R2.fastq.gz,54,control
Patient_7,P7_R1.fastq.gz,P7_R2.fastq.gz,61,case
```

### 3.1 Columns

| Column | Content |
|---|---|
| `sample` (required) | The sample's name, used in all results: count tables, coverage and BAM files, QC report. |
| `fastq_file` (required) | The reads (R1) of the sample. |
| `fastq_file2` (optional) | The second reads (R2) of a paired-end sample. Filled in for every sample: a paired-end run; empty for every sample (or no column): a single-end run. A mix stops the run. For paired-end reads, R1 and R2 overlap (small-RNA inserts are shorter than the reads) and are merged into one insert sequence. |
| any other columns | E.g. clinical data. Not used by the pipeline; the sheet is copied unchanged into the results as `sample_sheet.csv`. |

The column names are written exactly as above, in any order. The sheet may be comma- or semicolon-separated and must be saved as UTF-8 text (in Excel: *CSV UTF-8*). Spaces around values are ignored, and so are empty rows.

### 3.2 Rules

- **One folder:** the sheet and the FASTQ files are in the same folder, and `fastq_file` / `fastq_file2` hold file names only, without a path. The folder is only read; nothing is written into it.
- **FASTQ files:** gzip-compressed, ending in `.fastq.gz` or `.fq.gz`; names of letters, digits, `.`, `_` and `-`, not starting with `-` (no spaces). Each file belongs to one sample and one read.
- **Sample names:** letters, digits, `_` and `-`, not starting with `-`, e.g. `Healthy_1`, `1_healthy`, `P-07`. Each name once; names that differ only in upper and lower case count as the same, because they become file names. Not allowed are names the pipeline would read as a number or a missing value (e.g. `001`, `1e5`, `NA`, `true`; use `S001`), and names containing `NEB` or `SMART` (the pipeline would treat the sample as another library kit).
- **One row per sample:** samples sequenced on several lanes are merged first, e.g. `zcat A_L001_R1.fastq.gz A_L002_R1.fastq.gz | gzip > A_R1.fastq.gz` (the same for R2).

## 4. Running

### 4.1 With the starter

```bash
./ribomarker SHEET OUTPUT [options]         # Linux, macOS, Windows WSL shell
.\ribomarker SHEET OUTPUT [options]         # Windows PowerShell
./ribomarker --help                         # all options
```

`SHEET` is your sample sheet (§3), `OUTPUT` the results folder (created if missing). The starter:

- mounts the sheet's folder read-only and `OUTPUT` writable into the container, and the `--config` file and the `--references` folder read-only;
- uses 8 threads unless `--threads N` is given;
- runs the container without network access (`--network none`), removes it afterwards (`--rm`), and lets Ctrl-C stop a run cleanly (`--init`; `-it` in a terminal);
- on Linux, macOS and WSL, runs as your user (`--user`), so the results belong to you;
- refuses an output folder inside the sheet's folder, and says so when Docker is not installed or not running.

### 4.2 With `docker run`

The same as the starter does, written by hand:

```bash
# Linux, macOS, Windows WSL shell; create the results folder first
mkdir -p /path/to/results
docker run --rm -it --init --network none --user "$(id -u):$(id -g)" \
  -v /path/to/fastq:/data/input:ro -v /path/to/results:/data/output \
  IMAGE run /data/input/samples.csv /data/output --threads 8

# Windows PowerShell (one line, no --user)
docker run --rm -it --init --network none -v C:\data\fastq:/data/input:ro -v C:\data\results:/data/output IMAGE run /data/input/samples.csv /data/output --threads 8

# the demo
docker run --rm -it --init --network none -v /path/to/demo_out:/data/output IMAGE run --demo /data/output
```

The image accepts only the `run` command; `docker run --rm IMAGE run --help` shows its options. Here `/path/to/fastq` is the folder with the sample sheet (`samples.csv`) and the FASTQ files. Files and folders given as options must be mounted into the container too, e.g. `-v /path/to/settings.yaml:/data/config/settings.yaml:ro --config /data/config/settings.yaml`.

### 4.3 Options

| Option | Meaning |
|---|---|
| `--threads N` | CPU cores to use (default: 8; at most the cores available) |
| `--run-id NAME` | name of the results folder and of the report title (default: `ribomarker`; letters, digits, `_`, `-`, starting with a letter) |
| `--config FILE` | advanced settings (§5) |
| `--species NAME` | the reference species: `human` (default) or `mouse` with the built-in references (§9); with `--references`, the species name of your own files (§5, *Own references*) |
| `--references FOLDER` | your own references instead of the built-in ones (§5, *Own references*) |
| `--unlock` | remove a lock left by a run that was stopped hard (only if no other run uses `OUTPUT`) |
| `--demo` | run the demo (§2.4) instead of your data: `run --demo OUTPUT` |

### 4.4 Checks before the run starts

`run` stops at once, with a message saying what to do, when:

| Problem | What to do |
|---|---|
| the sample sheet is not found, not UTF-8 text, or lacks the columns `sample` / `fastq_file` | check the path; save it as *CSV UTF-8*; name the columns as in §3.1 |
| a row without a sample name or without `fastq_file`; `fastq_file2` filled in for some samples only | complete the sheet |
| a sample name breaks the rules of §3.2, or is given twice | rename it in the sheet |
| a FASTQ file is not found in the sheet's folder, is given with a path, has an unusable name, is not `.fastq.gz` / `.fq.gz`, or is listed twice | fix the sheet, or move or rename the file |
| the `--config` file is not found | check the path |
| `--config` contains a setting not listed in §5, or a value of the wrong kind | see §5 |
| `--threads` is not a whole number of 1 or more | e.g. `--threads 4` |
| the references are incomplete, or a FASTQ file is damaged (`gzip -t` fails) | rebuild or copy again |
| not enough free disk space (3 × the FASTQ size + 1 GB) | free space, or use another output folder |
| without Docker: a folder path contains a space | use a path without spaces |

It warns, and starts, when less than 12 GB of memory is available and when the references have no `MANIFEST_<species>.json`.

### 4.5 The end of a run

The run ends with one line:

- `RESULTS VALID: all <n> consistency checks equal …`: no consistency check of the QC report (§7) is `differs` or `not evaluable`;
- `RESULTS NOT VALID: …`: the pipeline stopped at a step (the line names the log), or a consistency check is `differs` or `not evaluable` (the line names the QC report);
- `RESULTS NOT CHECKED`: the QC report was switched off (`qc_report: False`, §5); exit status 0.

Exit status: `0` results valid (or not checked); `1` the pipeline failed or the results are not valid; `2` refused before starting; `130` interrupted.

### 4.6 Stopping and resuming

- **Ctrl-C** (also `docker stop`, or closing the terminal) stops a run cleanly and marks it `interrupted`.
- **Resume** by repeating the same command: finished steps are not redone; `--threads` may differ.
- A **finished** run is refused, so it is never redone by accident. Use another `--run-id` or another output folder for a new run.
- If a run was stopped hard (e.g. the computer restarted) and the next attempt says the output folder is locked, add `--unlock`.

## 5. Advanced settings

`--config FILE` takes a YAML file with any of these settings; the defaults are the validated ones.

| Setting | Default | Effect |
|---|---|---|
| `remove_multimap` | `True` | Leave out reads that align equally well to several loci (multimappers). With `False`, each multimapper is counted at **one** of its loci, chosen by a deterministic but arbitrary rule (§9), and there is no gene table: the counts are in `*-readcounts.csv`, with multimappers as separate `<feature>_multimap` rows. |
| `remove_antisense` | `True` | Leave out reads on the strand opposite to their feature. With `False` there is no gene table, as for `remove_multimap`, and the strand-agnostic lookup of §9 matters more (the QC report warns). |
| `multimap` | `100` | The maximum number of alignments reported per read (Bowtie 2 `-k`). |
| `keep_types` | `all` | RNA types written to the coverage files: `all`, or a comma-separated list without spaces (e.g. `miRNA,tRNA`). |
| `count_ribospike` | `True` | Count the RiboSpike spike-in reads (QC report, *Spike-ins*). |
| `qc_report` | `True` | Write the QC report (§7). With `False`, the run ends with `RESULTS NOT CHECKED`. |

Example:

```yaml
remove_multimap: True
keep_types: miRNA,tRNA,piRNA
```

Write `True` / `False` and numbers without quotes (`"False"` in quotes would count as true; `run` refuses it). Other settings (folders, threads, paired- or single-end, trimming) are set by `run` itself and are refused in `--config`.

### Own references

`--references FOLDER` replaces the built-in references with a folder of your own, e.g. with another miRNA or piRNA database, an extra spike-in, or another priority order. The starter mounts the folder read-only. The file names carry the species name given with `--species`: letters, digits, `_` and `-`, e.g. `--species rat` with `rat.dbpaths.csv` and `rat.bed`. For a species other than human and mouse, `run` warns and the QC report skips two species-specific items: the check that the BED holds no contigs named for another species (marked *not applicable*) and the red-blood-cell miRNA ratio. The folder holds:

| File | Content |
|---|---|
| `<species>.dbpaths.csv` | The references (tiers), one per line, **in priority order** (§8): `aligner,index,retain`, e.g. `bowtie2,mirna/human-mirna,T`. `aligner`: `bowtie2` (up to `multimap` alignments per read), `bowtie` (Bowtie 1, one alignment) or `ribospikebt2` (the spike-in settings), as in the built-in file. `index`: the index prefix, relative to the folder. `retain`: `T` keeps and counts the reads that align there; `F` removes them without counting. |
| Bowtie / Bowtie 2 indexes | One per line of the dbpaths file (`bowtie-build` or `bowtie2-build` on a FASTA file). |
| `<species>.bed` | The features: tab-separated, no header, six columns: sequence name (as in the index), start, end, feature name, biotype, strand. The biotype becomes the RNA type in the results and decides which of overlapping features a read gets (§8, step 3); biotypes the pipeline does not know rank **above** all known ones. |
| `MANIFEST_<species>.json` | Optional. Without it, `run` warns that the QC report cannot identify the references. |

`run` checks before the start that every index and the BED file exist. Where to start: a copy of the built-in references (§10, step 3), or the small demo references in the image (`/pipeline/demo/references`), changed one reference at a time; to rebuild references from their original sources, see [`buildDB/`](buildDB/README.md).

Caveats: the pipeline has been tested only with the built-in and the demo references. The gene table (`*-genecounts.csv`) and parts of the QC report rely on the biotype names and tiers of the built-in references (§8, step 5), so check the report of the first runs with your own references.

## 6. Outputs

All results are in `OUTPUT/<run id>/` (default `OUTPUT/ribomarker/`); the tables are named `<run id>-….csv`. Sample names are those of the sample sheet (§3).

### 6.1 Counts and coverage

| File | Content |
|---|---|
| `*-genecounts.csv` | **Counts per gene** (genes × samples): mature-exonic reads, collapsed to one row per gene (§8, step 5). Unique and sense reads only (with the default settings). Not written when `remove_multimap` or `remove_antisense` is `False`. |
| `*-genelengths.csv` | One length per gene, with its biotype. |
| `*-readcounts.csv` | Counts per reference feature, before the collapse to genes. |
| `*-typecounts.csv` | Counts per RNA type. |
| `*-typelengths.csv` | Read-length distribution per RNA type (lengths above 150 nt are counted as 150, §9). |
| `*-featuretypes.csv` | Feature → RNA type, with the feature length. |
| `bedgraph/<sample>.bedgraph` | Per-feature coverage consistent with the counts (the same reads; each multimapper at its one assigned locus). Columns: `feature start end depth type`. |
| `bedgraph_full/<sample>.bedgraph` | Coverage from every equally best locus of every read, with `multimap_depth` (how much of the depth comes from multimappers) and `strand`. Columns: `feature start end depth multimap_depth strand type`. |
| `<sample>.bam` | The assigned alignments; the assigned feature is in the `CS` tag. |

### 6.2 Statistics and logs

| File | Content |
|---|---|
| `*-mappingstats.csv` | Mapped, unmapped, unique and multimapping reads per sample |
| `*-trimstats.csv` | Trimming statistics per sample (adapter dimers, trimmed reads, merged pairs) |
| `ribospike/` | Spike-in counts by length and type |
| `trimmed/logs/` | cutadapt and (paired-end) fastp reports |
| `sample_sheet.csv` | Your sample sheet, unchanged (with any extra columns) |
| `samples.csv`, `run_config.yaml` | What the pipeline was given: its own sample list and the settings |
| `run.log` | The full log of the run |
| `run_status.json` | Verdict, start and end times, pipeline version, code fingerprint and git commit; in Docker also the memory peak |

The output folder also holds `.snakemake/`, the bookkeeping of the workflow manager (Snakemake), which is needed to resume a run.

### 6.3 QC files (`qc/`)

| File | Content |
|---|---|
| `qc/report_technical.html` | **The QC report** (§7): one self-contained HTML file that opens offline in any browser |
| `qc/metrics_long.tsv` | Every QC number, one row per sample and metric (`run_id, sample, section, metric_id, value, unit, pipeline_version`) |
| `qc/qc_metrics.json` | The same numbers with denominators, shares and flags, plus the consistency checks, warnings and provenance |
| `qc/consistency_checks.tsv` | Every consistency check, one row per check and sample (`check, status, sample, subject, expected, observed, difference, identity`) |
| `qc/run_config.yaml` | The settings as the pipeline resolved them |
| `qc/mqc/` | Intermediate files for MultiQC |

Every metric (name, definition, unit, denominator, caveats) is defined in [`docker/qc/metrics.yaml`](docker/qc/metrics.yaml).

## 7. The QC report

The report describes; it does not judge. **No quality thresholds are applied:** no value is marked pass or fail, and every percentage names its denominator (e.g. "% of merged reads"). Sections, top to bottom:

1. **Run overview:** the number of consistency checks per status, the checks that are not `equal` (only if there are any) and any warnings.
2. **General statistics:** one row per sample with the headline numbers (columns from fastp and cutadapt are hidden; *Configure columns* shows them).
3. **Read fate before mapping:** where the input reads (pairs) went: adapter dimers, inserts of 1–14 nt, removed by fastp filters and not merged (both paired-end only), into mapping. The categories add up to the input.
4. **Reads per reference tier:** reads assigned in each mapping tier (§8), before multimappers and antisense reads are removed.
5. **Composition by RNA type:** kept reads per RNA type after that removal; further tabs without spike-ins, and without rRNA and spike-ins.
6. **Multimappers and antisense reads:** the mapped reads split into kept reads, reads removed as multimappers and reads removed as antisense (the last is derived, not counted directly); below it, per tier, the share of reads with a second alignment.
7. **Read length by RNA type:** all types together, then one tab per type.
8. **Spike-ins:** share of reads, fragments, length fidelity, end-base composition of exact-length reads (with its chance level), length offset and recovery per designed length. Without the amounts added, recovery can only be compared between samples and lengths.
9. **Detection and counts:** library size (with and without rRNA), genes detected with more than 0, 5 and 10 reads, the ratio of two red-blood-cell-enriched miRNAs (miR-451a / miR-23a-3p; descriptive, not comparable to qPCR cut-offs without calibration), and the 20 most abundant RNAs per sample, overall and per RNA type.

These are followed by MultiQC's own fastp (paired-end) and cutadapt sections, **Provenance** (library parameters, reference species check, reference build, settings) and **Software versions** (the pipeline version and every tool version).

Plots with one line per sample have no legend in MultiQC's own sections: hover over a line to see its sample, or use the toolbox on the right (*Highlight samples*, *Show / hide samples*).

Two kinds of marks need attention:

- **Consistency checks** compare the same quantity in the outputs of two steps (e.g. reads written by trimming = reads entering mapping). `equal` shows that read totals are conserved and the files are complete; it does not show that reads were classified correctly. `differs` means two outputs disagree, which points to a pipeline problem rather than a sample problem: please open an issue. `not evaluable` means a number the check needs could not be read; like `differs`, it makes the run `RESULTS NOT VALID`. `not applicable` is normal for checks that do not apply to the layout (e.g. the fastp check in single-end runs). Every check is in `qc/consistency_checks.tsv`.
- **`low n` / `not evaluable`:** fewer than 100 reads, or none, are behind a share, so it is not reliable. The report shows these flags in the spike-in and counts tables; all flags are in `qc/qc_metrics.json`.

### Sending the report

For help with a run, send `qc/report_technical.html` (and `qc/qc_metrics.json` if asked). Absolute paths are removed from the report. It contains sample names and FASTQ file names, so check these before sending it outside your organisation. It contains aggregate numbers, not read sequences (the pipeline does not switch on fastp's over-represented-sequence analysis).

## 8. How the pipeline works

The workflow is run by Snakemake; each step is a Python script in `docker/scripts/`. The steps, in order:

1. **Trimming.** cutadapt removes the 3′ adapter (`TGGAATTCTCGGGTGCCAAGG`; R2: `TGATCGTCGGACT`) and the one extra base that the library protocol adds at the insert's 5′ end, and drops inserts shorter than 15 nt (reported as adapter dimers and short inserts). For paired-end reads, fastp then merges R1 and R2 into one insert (overlap of at least 10 nt, at most 2 differences) and filters reads with its default criteria (low quality, too many N); pairs that cannot be merged are not used. No base-quality trimming is done.
2. **Tiered mapping.** Reads are aligned to a stack of references **in a fixed order**, and each read is assigned to the *first* reference it aligns to:
   `UniVec → RiboSpike → rRNA → tRNA → miRNA → miRNA star → miRNA precursor → piRNA → ncRNA → cDNA → genome`.
   Contaminants, spike-ins and rRNA are claimed first; reads that do not align go on to the next reference. UniVec and rRNA are aligned with Bowtie 1, reporting one alignment; the other tiers with Bowtie 2 (up to `multimap` alignments per read).
3. **Feature assignment.** For each read the pipeline resolves, in order: (a) the reference (the tier order above); (b) the overlapping feature, by a ranking of biotypes in which structural RNAs outrank others and host mRNA ranks lowest; (c) the locus of a multimapper, by a deterministic tie-break. Multimappers and antisense reads are then removed (§5).
4. **Counts.** Reads are counted per feature (`*-readcounts.csv`) and per RNA type (`*-typecounts.csv`).
5. **Gene counts.** The feature counts are collapsed to a gene-level, mature-exonic table: precursor and intronic features are dropped, spike-ins excluded, tRNA fragments summed to their gene, while near-identical rRNA copies and miRNA arms are kept separate.
6. **Coverage.** Per-feature coverage is written from the assigned alignments (`bedgraph/`) and from all equally best alignments (`bedgraph_full/`), see §6.1.
7. **QC.** A collector reads the outputs of all steps, writes every QC number with its definition and runs the consistency checks; MultiQC renders the report. This step only reads existing outputs.

### References

The image contains the human (GRCh38) and mouse (GRCm39) references with their Bowtie indexes: NCBI UniVec, the RiboSpike spike-in sequences, rRNA (from miRTrace, with sequences from NCBI and SILVA), tRNAs (GtRNAdb), miRNAs (MirGeneDB), piRNAs (piRBase), and ncRNA, cDNA and genome (Ensembl). The transcript references use **one representative (longest) transcript per gene**. The exact files, source versions and checksums are in `MANIFEST_<species>.json` and `PROVENANCE_<species>.md` in [`docker/assets/`](docker/assets/); the QC report shows the reference build. Sources, citations and terms: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The references were built with the scripts in [`buildDB/`](buildDB/).

### Versions and reproducibility

The pipeline version is in [`docker/VERSION`](docker/VERSION); the git tag `v<version>` and the image tag `<version>` are the same. Each run records the version, the git commit, a fingerprint of the code and the reference build in `run_status.json` and in the QC report. The same image (by digest) with the same references gives the same results; the image digest of each release is in its release notes.

## 9. Known limitations

Read these before interpreting results.

- **RNA type depends on the tier order.** A read is assigned to the first reference it aligns to (§8), so the boundaries between RNA types (e.g. piRNA, lncRNA, protein-coding) partly reflect the order of the references, not only biology.
- **Chance alignment of short reads.** Reads of 15–17 nt can align by chance to the large ncRNA and cDNA references before they reach the genome tier, and are then counted as lncRNA or mRNA. Under a uniform-random model, a 15-nt read has about an 11 % chance of an exact match somewhere in the 62.5 Mb cDNA reference, and nearly 100 % if one mismatch is allowed; at 20 nt about 0.01 % and 0.7 %. These are model estimates, not measurements. The QC report shows the 15–17 nt share per RNA type.
- **Strand-agnostic feature lookup in the genome tier.** Features of reads assigned in the genome tier are looked up without regard to strand, which affects the sense/antisense labels and counts of pre-mRNA and pre-tRNA features.
- **Equal-rank overlapping features.** When a read overlaps two features of the same rank, which one wins is not specified (an unstable sort). In our tests outputs were identical between runs and between local and Docker runs, and the order affected a few reads per sample; whether it is the same with other package versions or CPUs is not established.
- **Strand-specific libraries assumed.** The sense/antisense distinction assumes a strand-specific library.
- **Multimappers.** By default multimappers are left out of the counts. With `remove_multimap: False`, each multimapper is placed at one locus by a deterministic but arbitrary tie-break; no fractional or EM-based assignment is done. In `bedgraph_full/`, `multimap_depth` shows how much depth comes from multimappers.
- **UniVec and rRNA report one alignment** (Bowtie 1), so multimapping within these references is resolved by the aligner and does not appear in `multimap_depth` or in the multimapper counts.
- **One transcript per gene.** Isoform-level information (alternative starts and 3′ ends) is not represented.
- **Composition after removal.** The composition and detection numbers of the QC report are measured after multimappers and antisense reads are removed.
- **Read lengths above 150 nt** are counted as 150 nt in `*-typelengths.csv`.
- **Mouse not yet compared with real data.** For mouse, the mouse references were checked on synthetic reads cut from them (every read lands on its own reference), but no result on real mouse data has been compared yet; the QC report of a mouse run says so. Other species: only with your own references (§5, *Own references*), not tested by us.

## 10. Running without Docker

On Linux x86-64 (including WSL2 on Windows) the pipeline can run in a conda environment. Other systems are not supported for this route.

1. Install [Miniforge](https://github.com/conda-forge/miniforge) (conda with the conda-forge channel).
2. Download the repository and create the environment from the lock file, which lists the exact package builds used in the image:
   ```bash
   conda create -n ribomarker --file docker/env.lock.txt
   conda activate ribomarker
   ```
3. Copy the references out of the image (this needs Docker once, and the disk space of the references):
   ```bash
   docker create --name ribomarker-refs IMAGE
   docker cp ribomarker-refs:/pipeline/assets ./references
   docker rm ribomarker-refs
   ```
4. Test and run (paths must not contain spaces):
   ```bash
   bash docker/make_demo.sh                                 # builds the demo data once
   python docker/scripts/main.py run --demo results/demo
   python docker/scripts/main.py run /path/to/fastq/samples.csv /path/to/results --references ./references --threads 8
   ```

All options of §4.3 apply. Without `--references`, the pipeline looks for the references in `docker/assets/`.

## 11. Building the image

To build the image yourself, from the repository root with the conda environment of §10 active:

```bash
# references copied out of the image into ./references, as in §10 step 3
python docker/stage_references.py --assets references     # checks them against the MANIFESTs → staging/references
bash docker/make_demo.sh
docker build --provenance=false --build-context references=staging/references --build-context root=. \
  --build-arg GIT_COMMIT=$(git rev-parse HEAD) --build-arg VERSION=$(cat docker/VERSION) \
  -t ribomarker:local docker
```

The references can also be rebuilt from their sources with `buildDB/` (see its README), but the source files may have changed since, and the result is then not the published reference build.

## 12. Support, citing and licence

- **Questions and problems:** open an issue on GitHub, with the QC report if you can (§7, *Sending the report*). Contributions: [CONTRIBUTING.md](CONTRIBUTING.md).
- **Citing:** see [CITATION.cff](CITATION.cff) and the [README](README.md#citing); cite the reference databases as listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
- **Licence:** MIT ([LICENSE](LICENSE)); logos and names: [NOTICE](NOTICE); third-party software and data: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The image contains the same three files in `/pipeline/`, and the licence texts of its packages in `/opt/conda/envs/ribomarker/share/conda-licenses/`.

## 13. Glossary

- **cfRNA:** cell-free RNA, e.g. circulating in plasma.
- **Fragmentomics:** the study of which parts of an RNA are present as fragments, not only how much of it.
- **Biotype / RNA type:** the RNA class of a feature (miRNA, tRNA, rRNA, lncRNA, protein-coding, …).
- **Feature:** an annotated sequence in a reference (a transcript, a miRNA, a tRNA, an intron, …).
- **Tiered mapping:** alignment to an ordered stack of references; the first reference a read aligns to wins.
- **Sense / antisense:** a read on the same / opposite strand as its feature.
- **Multimapper:** a read that aligns equally well to several loci.
- **Mature-exonic:** restricted to the processed transcript (no intronic or precursor reads).
- **Coverage / bedgraph:** the number of reads covering each position of a feature.
- **Spike-in (RiboSpike):** synthetic control RNAs of known sequence and length added to the sample.
- **Consistency check:** a comparison of the same quantity in the outputs of two pipeline steps.
