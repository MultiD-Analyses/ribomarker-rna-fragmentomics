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

- **counts** per feature (the gene table: one row per gene, miRNA arm, tRNA family, …), per reference sequence and per RNA class;
- **per-transcript coverage** files (bedgraph): how many reads cover each position of each transcript;
- a **technical QC report**: every step's numbers, with consistency checks between the steps.

The pipeline stops there. It applies no normalization and no statistics; the counts of the gene table are integer counts, suitable for count-based differential-expression tools. When comparing samples, keep in mind that library sizes differ: the method paper subsampled all libraries of an experiment to the read count of the smallest one before comparing them (Clark et al. 2026, Methods, *Read processing and feature counting*); the pipeline does not, so subsample or normalize in your own analysis.

Scope of this version:

- **Library protocol:** RiboMarker (RealSeq) only. The adapter sequences and the removal of one extra base at the insert's 5′ end are fixed for this protocol (§8).
- **Species:** human (default) and mouse (`--species mouse`). The references of both are built into the image (§8, *References*); other species with your own references (§5, *Own references*). For mouse, the mouse references were checked on synthetic reads cut from them (every read lands on its own reference), but no result on real mouse data has been compared yet; the QC report of a mouse run says so (§9).
- **Systems:** the image is built for x86-64 (Intel/AMD). It has been tested on Windows (two computers; PowerShell and WSL2), on Linux under WSL2, and on macOS on an Intel Mac (§8, *Versions and reproducibility*). Not tested: a native Linux installation, and Apple-silicon Macs (M1 and later), where Docker runs x86-64 images under emulation, which is slower.

## 2. Installation

### 2.1 Requirements

| | |
|---|---|
| Docker | Docker Desktop (Windows, macOS) or Docker Engine (Linux), installed and running |
| Docker Desktop licence | Free for companies with fewer than 250 employees and less than US$10 million annual revenue, and for non-commercial open-source projects; other companies and public bodies need a paid Docker subscription ([Docker pricing FAQ](https://www.docker.com/pricing/faq/)). The condition is Docker's, for Docker Desktop; it does not concern the pipeline (MIT licence). |
| Memory | 12 GB available to Docker. The peak grows with sample size: with 8 threads, a paired-end run of 2 samples (0.30 GB of R1 files) peaked at 4.1 GB, one of 6 samples (1.72 GB, the largest sample 0.67 GB) at 8.1 GB. `run` warns when less than 12 GB is available. Docker Desktop gives containers part of the computer's memory: check *Settings → Resources*. |
| Disk | about 21 GB for the image, as Docker Desktop reports it (human and mouse references); per run about 3 × the size of the FASTQ files + 1 GB (checked before the run starts) |
| Time | about 2–3 hours per GB of R1 FASTQ files with 8 threads. Two measurements on one laptop: 37–41 min for 2 samples (0.30 GB of R1 files), and 4 h 41 min for 6 samples (1.72 GB) through the Windows starter. It depends on the computer and the samples; samples are processed one after another. |

### 2.2 Docker on each system

- **Linux:** install Docker Engine (or Docker Desktop) and start it. Your user must be allowed to run Docker (the `docker` group, or `sudo`).
- **Windows:** install Docker Desktop with the **WSL2 backend**. You can then use either **PowerShell** with Windows paths (`C:\data\...`), or a **WSL shell** (e.g. Ubuntu from the Microsoft Store, with Docker's WSL integration switched on under *Settings → Resources → WSL integration*) with Linux paths (`/mnt/c/data/...`). In the WSL shell, if commands fail with "No such file or directory" after a run, `cd` into your folder again by its full path.
- **macOS (tested on an Intel Mac):** install and start Docker Desktop and use Terminal with paths such as `/Users/you/data`. On Apple-silicon Macs (not tested), Docker runs the x86-64 image under emulation, which is slower. Give Docker at least 12 GB of memory (*Settings → Resources*). Docker Desktop shares your home folder by default; other locations may need adding under *Settings → Resources → File sharing*.

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
sample,fastq_file,fastq_file2,group,age,sex
Healthy_1,H1_R1.fastq.gz,H1_R2.fastq.gz,Healthy,25,M
Patient_1,P1_R1.fastq.gz,P1_R2.fastq.gz,Sick,27,F
```

### 3.1 Columns

| Column | Content |
|---|---|
| `sample` (required) | The sample's name, used in all results: count tables, coverage and BAM files, QC report. |
| `fastq_file` (required) | The reads (R1) of the sample. |
| `fastq_file2` (optional) | The second reads (R2) of a paired-end sample. Filled in for every sample: a paired-end run; empty for every sample (or no column): a single-end run. A mix stops the run. For paired-end reads, R1 and R2 overlap (small-RNA inserts are shorter than the reads) and are merged into one insert sequence. |
| any other columns | E.g. clinical data. Not used for processing. The sheet is copied unchanged into the results as `sample_sheet.csv`, and **all its columns are shown in the QC report** (section *Samples*): leave out what the report's readers should not see. In the report, *Colour samples by* colours the samples by any of these columns; `report_group` (§5) sets the first choice. |

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

Folder paths with spaces (e.g. `C:\Users\First Last\data`) work through the starters; tested with the Windows starter. Without Docker (§10) they do not.

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
| `--tiers FILE` | a tier list replacing the one of the references (§5, *Tier list and RNA class table*) |
| `--rna-classes FILE` | an RNA class table replacing the default one (§5, *Tier list and RNA class table*) |
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
| `report_group` names a column that is not in the sample sheet | name an existing column |
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

The output of the mapping step may contain `[WARNING] Failed to launch x86-64-v3 version, staying with default`: a program looks for a build optimised for newer processors and uses its standard build instead. It does not affect the results; it appeared on every computer we tested.

### 4.6 Stopping and resuming

- **Ctrl-C** (also `docker stop`, or closing the terminal) stops a run cleanly and marks it `interrupted`.
- **Resume** by repeating the same command: finished steps are not redone; `--threads` may differ.
- A **finished** run is refused, so it is never redone by accident. Use another `--run-id` or another output folder for a new run.
- If a run was stopped hard (e.g. the computer restarted) and the next attempt says the output folder is locked, add `--unlock`.
- A run started with version 0.2.x cannot be resumed with version 0.3.0 or later: start it again in an empty folder.

## 5. Advanced settings

`--config FILE` takes a YAML file with any of these settings; the defaults are the validated ones.

| Setting | Default | Effect |
|---|---|---|
| `remove_multimap` | `True` | Leave out reads that align equally well to several loci (multimappers). With `False`, each multimapper is counted at **one** of its loci, chosen by a deterministic but arbitrary rule (§9), and there is no gene table: the counts are in `*-readcounts.csv`, with multimappers as separate `<feature>_multimap` rows. |
| `remove_antisense` | `True` | Leave out reads on the strand opposite to their feature. With `False` there is no gene table, as for `remove_multimap`, and the strand-agnostic lookup of §9 matters more (the QC report warns). |
| `multimap` | `100` | The maximum number of alignments reported per read (Bowtie 2 `-k`). |
| `keep_classes` | `all` | RNA classes written to the coverage files: `all`, or a comma-separated list of class names without spaces (e.g. `miRNA,tRNA`). The old name `keep_types` is still accepted, with a warning. |
| `count_ribospike` | `True` | Count the RiboSpike spike-in reads (QC report, *Spike-ins*). |
| `qc_report` | `True` | Write the QC report (§7). With `False`, the run ends with `RESULTS NOT CHECKED`. |
| `report_group` | none | A column of the sample sheet: the column the QC report first colours the samples by (*Colour samples by* in section *Samples* changes it; up to 8 groups, with more the samples are not coloured). |
| `tiers`, `rna_classes` | none | The tier list and the RNA class table (§5, *Tier list and RNA class table*). With the starter or Docker, give them as `--tiers` and `--rna-classes`, which mount the files. |

Example:

```yaml
remove_multimap: True
keep_classes: miRNA,tRNA,piRNA
report_group: group
```

Write `True` / `False` and numbers without quotes (`"False"` in quotes would count as true; `run` refuses it). Other settings (folders, threads, paired- or single-end, trimming) are set by `run` itself and are refused in `--config`.

### Tier list and RNA class table

Two tables decide how an aligned read becomes an RNA class:

- **The tier list** names the references in their mapping order (§8, step 2): `<species>.dbpaths.csv` in the references folder, one line per reference, `aligner,index,retain` (format: *Own references* below).
- **The RNA class table** turns the labels of the reference BED (its fifth column, e.g. Ensembl's biotype) into RNA classes: `<species>.rna_classes.tsv`, tab-separated, one row per label (a label found in several references may have one row per reference; in this version these rows must agree on class and rank, since a read's class does not depend on its reference). Its columns:

| Column | Content |
|---|---|
| `source_label` | the label in the BED, or a label the pipeline sets itself (`unannotated`, `unmapped`, the spike-ins, …) |
| `tier` | the reference the label comes from (the index name without the species, e.g. `ncrna`), or `*` for all; tiers that the tier list does not have are named in a warning |
| `class_id` | the RNA class written into the result tables (no spaces) |
| `display_name` | the class's name in the QC report |
| `group` | its group in the report's tier and composition plots |
| `display_order` | the order of the classes in the report |
| `rank_in_tier` | which of overlapping features of one reference a read gets (higher wins) |
| `positions` | `yes`: reads are split by their place in the feature (`_5p`, `_3p`, `_internal`, `_whole`; §8, step 3) |
| `gene_table` | `yes`: in the gene table; `drop`: precursor or gene-span classes, left out; `exclude`: not a gene (spike-ins, vector contaminants, unannotated, unmapped) |
| `counts_as_rrna` | `yes`: left out of the "without rRNA" numbers of the QC report |
| `description` | one line for the report |

The built-in tables are in [`docker/tables/`](docker/tables/). Every run copies the two tables it used into `tables/` of its results folder (`rna_classes.tsv`, `tiers.csv`), and the QC report shows them (section *Tiers and RNA classes*; files and checksums under *Provenance*). To change them, start from these copies and give the changed file with `--tiers FILE` or `--rna-classes FILE`:

```bash
./ribomarker ~/data/fastq/samples.csv ~/data/results --rna-classes my_classes.tsv
```

Which class table is used (the first found): `--rna-classes` (or the setting `rna_classes`); a `<species>.rna_classes.tsv` in the references folder; the built-in table of the species; `docker/tables/default.rna_classes.tsv`. A table replaces the built-in one as a whole; rows for the labels the pipeline sets itself are added from the built-in table when missing. BED labels missing from the table become a class of the same name, with rank 10 (the highest built-in rank), the group "lncRNA and other ncRNA", the gene table, and a warning. The paths of a replaced tier list are read relative to the references folder.

Changing the tier order or the ranks changes which reads are assigned to which feature. The pipeline has been tested with the built-in tables only; check the QC report of the first runs with your own.

### Own references

`--references FOLDER` replaces the built-in references with a folder of your own, e.g. with another miRNA or piRNA database or an extra spike-in. The starter mounts the folder read-only. The file names carry the species name given with `--species`: letters, digits, `_` and `-`, e.g. `--species rat` with `rat.dbpaths.csv` and `rat.bed`. For a species other than human and mouse, `run` warns and the QC report skips two species-specific items: the check that the BED holds no contigs named for another species (marked *not applicable*) and the red-blood-cell miRNA ratio. The folder holds:

| File | Content |
|---|---|
| `<species>.dbpaths.csv` | The references (tiers), one per line, **in priority order** (§8): `aligner,index,retain`, e.g. `bowtie2,mirna/human-mirna,T`. `aligner`: `bowtie2` (up to `multimap` alignments per read), `bowtie` (Bowtie 1, one alignment) or `ribospikebt2` (the spike-in settings), as in the built-in file. `index`: the index prefix, relative to the folder. `retain`: `T` keeps and counts the reads that align there; `F` removes them without counting. |
| Bowtie / Bowtie 2 indexes | One per line of the dbpaths file (`bowtie-build` or `bowtie2-build` on a FASTA file). |
| `<species>.bed` | The reference sequences: tab-separated, no header, six columns: sequence name (as in the index), start, end, name, label, strand. The RNA class table turns the label into an RNA class and decides which of overlapping sequences a read gets (§8, step 3); labels missing from the table get a warning and rank 10 (above *Tier list and RNA class table*). |
| `<species>.rna_classes.tsv` | Optional: the RNA class table for these references (above). Without it, the built-in table of the species is used, or the default table for other species. |
| `MANIFEST_<species>.json` | Optional. Without it, `run` warns that the QC report cannot identify the references. |

`run` checks before the start that every index and the BED file exist. Where to start: a copy of the built-in references (§10, step 3), or the small demo references in the image (`/pipeline/demo/references`), changed one reference at a time; to rebuild references from their original sources, see [`buildDB/`](buildDB/README.md).

Caveats: the pipeline has been tested only with the built-in and the demo references. The gene table (`*-genecounts.csv`) and the QC report follow the RNA class table (§8, step 5): give your own labels rows in it, and check the report of the first runs with your own references.

## 6. Outputs

All results are in `OUTPUT/<run id>/` (default `OUTPUT/ribomarker/`); the tables are named `<run id>-….csv`. Sample names are those of the sample sheet (§3).

### 6.1 Counts and coverage

| File | Content |
|---|---|
| `*-genecounts.csv` | **The gene table: counts per feature** (features × samples; column `gene`): one row per gene, miRNA arm, tRNA family, …, from mature-exonic reads (§8, step 5). Unique and sense reads only (with the default settings). Not written when `remove_multimap` or `remove_antisense` is `False`. |
| `*-genelengths.csv` | The length (`length_nt`), RNA class (`rna_class`) and gene symbols (below) of every feature of the gene table. |
| `*-readcounts.csv` | Counts per reference sequence (column `feature`), before the collapse into the gene table. tRNA, snRNA and snoRNA sequences have one row per position (`_5p`, `_3p`, `_internal`, `_whole`). |
| `*-featureclasses.csv` | The RNA class (`rna_class`), length (`length_nt`) and gene symbols (below) of every reference sequence of `*-readcounts.csv`. |
| `*-featuredistinctseqs.csv` | The number of distinct read sequences per reference sequence. |
| `*-classcounts.csv` | Counts per RNA class (`rna_class`), including the spike-ins, vector contaminants (`UniVec`), `genome_unannotated` and `unmapped`. |
| `*-classlengths.csv` | Read-length distribution per RNA class, from 0 nt to the longest read of the run (at least 150 nt). |
| `*-classdistinctseqs.csv` | The number of distinct read sequences per RNA class. |
| `bedgraph/<sample>.bedgraph` | Per-sequence coverage consistent with the counts (the same reads; each multimapper at its one assigned locus). Columns: `feature start end depth rna_class`. |
| `bedgraph_full/<sample>.bedgraph` | Coverage from every equally best locus of every read, with `multimap_depth` (how much of the depth comes from multimappers) and `strand`. Columns: `feature start end depth multimap_depth strand rna_class`. |
| `<sample>.bam` | The assigned alignments; the assigned reference sequence is in the `CS` tag. |
| `tables/` | The tier list (`tiers.csv`) and the RNA class table (`rna_classes.tsv`) the run used (§5). |

**Gene symbols.** `*-featureclasses.csv` and `*-genelengths.csv` have the columns `hgnc_id` and `hgnc_symbol` (human) or `mgi_id` and `mgi_symbol` (mouse), and `mirbase_name` (the miRBase name of a miRNA arm, star arm or precursor, e.g. `hsa-miR-21-5p`). They come from `docker/tables/<species>.gene_symbols.tsv` (sources, versions and checksums in its `.provenance.md`), joined by Ensembl gene id (no second try by gene name, so a symbol is never taken from another gene of the same name); pre-mRNA through the gene names of the Ensembl GTF, miRNA arms and precursors through MirGeneDB's miRBase accessions, tRNAs through GtRNAdb ids. A feature that spans several genes (e.g. a tRNA family) lists them separated by `;`; in `*-genelengths.csv` a feature gets the first non-empty value of its reference sequences. How many have a value (built-in references): human protein_coding and pre-mRNA 97 %, snRNA > 99 %, miRNA arms and precursors 98 % (miRBase names: 96 % of the arms, 85 % of the star arms), tRNA 98 %, snoRNA 51 %, lncRNA 19 % (most human lncRNAs have no approved symbol); mouse protein_coding and pre-mRNA > 99 %, lncRNA 98 %, miRNA and tRNA none (MGI has no miRBase or GtRNAdb cross-references; the miRBase names are filled in, 98–100 %). piRNA, the rRNA copies of miRTrace, spike-ins and vector contaminants have no symbol. With your own references (§5) the columns stay empty, and the run warns, unless the table was built for your BED (`buildDB/make_gene_symbols.py`).

**Joining the tables.** The count tables hold counts only; the annotations are in separate tables with the same first column. Join `*-genecounts.csv` with `*-genelengths.csv` on `gene`, and `*-readcounts.csv` with `*-featureclasses.csv` on `feature`, e.g. in R `merge(counts, lengths, by = "gene")`, or in Python `counts.merge(lengths, on="gene")` (pandas).

### 6.2 Statistics and logs

| File | Content |
|---|---|
| `*-mappingstats.csv` | Mapped, unmapped, unique and multimapping reads per sample |
| `*-trimstats.csv` | Trimming statistics per sample (adapter dimers, trimmed reads, merged pairs) |
| `ribospike/` | Spike-in counts by length and RNA Type (`*-spikelengths.csv`) and the origin of every read: natural, spike-in or unmapped (`*-readorigin.csv`) |
| `trimmed/logs/` | cutadapt and (paired-end) fastp reports |
| `sample_sheet.csv` | Your sample sheet, unchanged (with any extra columns) |
| `samples.csv`, `run_config.yaml` | What the pipeline was given: its own sample list and the settings |
| `run.log` | The full log of the run |
| `run_status.json` | Verdict, start and end times, pipeline version, code fingerprint and git commit, the tier list and RNA class table used (`tables`); in Docker also the memory peak |

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

1. **Samples:** the columns of the sample sheet (the FASTQ file names are hidden at first; *Configure columns* shows them), and *Colour samples by*. Its first line gives the result of the consistency checks and the number of warnings.
2. **How to read this report:** the settings above the plots, the toolbox, and a glossary.
3. **General Statistics:** one row per sample with the headline numbers, including the mean read length of the natural reads; below the heading, the checks that are not `equal` and the warnings (only if there are any).
4. **Trimming and merging: read fate:** where the input reads (pairs) went: adapter dimers, inserts of 1–14 nt, removed by fastp filters and not merged (both paired-end only), into mapping. The categories add up to the input.
5. **Mapping per reference tier:** reads assigned in each mapping tier (§8), before multimappers and antisense reads are removed. The genome tier is split into the unannotated reads and the rest (pre-mRNA, pre-tRNA).
6. **Composition by RNA class:** kept reads per RNA class after that removal, grouped as in the RNA class table; further tabs without spike-ins, and without rRNA and spike-ins.
7. **Multimappers and antisense reads:** the mapped reads split into kept reads, reads removed as multimappers and reads removed as antisense (the last is derived, not counted directly); below it, per tier, the share of reads with a second alignment.
8. **Read length by RNA class**, in two modes (setting *Values*): **% of natural reads**, the reads of a class at a length as a share of all natural reads of the sample (all classes and lengths together add up to 100 %), and **% within class**, the length distribution of each class on its own. Natural reads are the kept reads without spike-ins and vector contaminants. The method paper shows such sequencing length profiles without stating their denominator (Clark et al. 2026); here it is the natural reads, counted after multimappers and antisense reads are removed (the default settings). One button per main class; the other classes are together under "Other RNAs" (the description lists them).
9. **Spike-ins:** share of reads, fragments, length fidelity, end-base composition of exact-length reads (with its chance level), length offset and recovery per designed length. Without the amounts added, recovery can only be compared between samples and lengths.
10. **Detection and counts:** library size (with and without rRNA), features detected with more than 0, 5 and 10 reads (with a plot per RNA class), the ratio of two red-blood-cell-enriched miRNAs (miR-451a / miR-23a-3p; descriptive, not comparable to qPCR cut-offs without calibration), and the 20 most abundant features per sample, overall and per RNA class.
11. **Tiers and RNA classes:** the tier list and the RNA class table of the run (§5).

These are followed by MultiQC's own fastp (paired-end) and cutadapt sections, **Provenance** (run, library parameters, reference species check, reference build, the tables, settings, and how to cite the pipeline and the method) and **Software versions** (the pipeline version and every tool version).

*Colour samples by*, above the sample table, colours the samples by the values of a sample-sheet column in every plot and table (at most 8 groups; `report_group`, §5, sets the first choice); a click on a colour in the key next to it changes that colour. The choice and the colours are not kept when the report is opened again, and it replaces colours set by hand in the toolbox's *Highlight* tab. Plots with one line per sample have no legend in MultiQC's own sections: hover over a line to see its sample, or use the toolbox on the right (*Highlight samples*, *Show / hide samples*).

Two kinds of marks need attention:

- **Consistency checks** compare the same quantity in the outputs of two steps (e.g. reads written by trimming = reads entering mapping). `equal` shows that read totals are conserved and the files are complete; it does not show that reads were classified correctly. `differs` means two outputs disagree, which points to a pipeline problem rather than a sample problem: please open an issue. `not evaluable` means a number the check needs could not be read; like `differs`, it makes the run `RESULTS NOT VALID`. `not applicable` is normal for checks that do not apply to the layout (e.g. the fastp check in single-end runs). Every check is in `qc/consistency_checks.tsv`.
- **`low n` / `not evaluable`:** fewer than 100 reads, or none, are behind a share, so it is not reliable. The report shows these flags in the spike-in and counts tables; all flags are in `qc/qc_metrics.json`.

### Sending the report

For help with a run, send `qc/report_technical.html` (and `qc/qc_metrics.json` if asked). Absolute paths are removed from the report. It contains the sample names, the FASTQ file names and **all columns of your sample sheet** (so does `qc/qc_metrics.json`), so check these before sending it outside your organisation. It contains aggregate numbers, not read sequences (the pipeline does not switch on fastp's over-represented-sequence analysis).

## 8. How the pipeline works

The workflow is run by Snakemake; each step is a Python script in `docker/scripts/`. The steps, in order:

1. **Trimming.** cutadapt removes the 3′ adapter (`TGGAATTCTCGGGTGCCAAGG`; R2: `TGATCGTCGGACT`) and the one extra base that the library protocol adds at the insert's 5′ end, and drops inserts shorter than 15 nt (reported as adapter dimers and short inserts). For paired-end reads, fastp then merges R1 and R2 into one insert (overlap of at least 10 nt, at most 2 differences) and filters reads with its default criteria (low quality, too many N); pairs that cannot be merged are not used. No base-quality trimming is done.
2. **Tiered mapping.** Reads are aligned to a stack of references **in a fixed order**, and each read is assigned to the *first* reference it aligns to:
   `UniVec → RiboSpike → rRNA → tRNA → miRNA → miRNA star → miRNA precursor → piRNA → ncRNA → cDNA → genome`.
   Contaminants, spike-ins and rRNA are claimed first; reads that do not align go on to the next reference. The order is the tier list (§5). The miRNA star arms are a reference of their own but are counted as the class miRNA. UniVec and rRNA are aligned with Bowtie 1, reporting one alignment; the other tiers with Bowtie 2 (up to `multimap` alignments per read).
3. **Assignment.** For each read the pipeline resolves, in order: (a) the reference (the tier order above); (b) the overlapping reference sequence, by the rank of its label (`rank_in_tier` of the class table), in which structural RNAs outrank others and protein_coding ranks lowest; (c) the locus of a multimapper, by a deterministic tie-break. Multimappers and antisense reads are then removed (§5). Reads of classes with `positions` (tRNA, snRNA, snoRNA) are split by their place in the sequence: `_5p` or `_3p` when the read lies within 10 nt of that end only, `_whole` when it lies within 10 nt of both ends, `_internal` otherwise; the 5′ end is on the sequence's strand. Except for `_whole`, this follows the method paper (Clark et al. 2026, Methods).
4. **Counts.** Reads are counted per reference sequence (`*-readcounts.csv`) and per RNA class (`*-classcounts.csv`).
5. **The gene table.** The counts per reference sequence are collapsed into a mature-exonic table of features (`*-genecounts.csv`): classes with `gene_table` `drop` (pre-miRNA, pre-mRNA, pre-tRNA) are left out, those with `exclude` (spike-ins, vector contaminants, unannotated, unmapped) are not genes, the position rows of tRNA, snRNA and snoRNA are summed per sequence, while near-identical rRNA copies and miRNA arms are kept separate.
6. **Coverage.** Per-sequence coverage is written from the assigned alignments (`bedgraph/`) and from all equally best alignments (`bedgraph_full/`), see §6.1.
7. **QC.** A collector reads the outputs of all steps, writes every QC number with its definition and runs the consistency checks; MultiQC renders the report. This step only reads existing outputs.

### References

The image contains the human (GRCh38) and mouse (GRCm39) references with their Bowtie indexes: NCBI UniVec, the RiboSpike spike-in sequences, rRNA (from miRTrace, with sequences from NCBI and SILVA), tRNAs (GtRNAdb), miRNAs (MirGeneDB), piRNAs (piRBase), and ncRNA, cDNA and genome (Ensembl). The transcript references use **one representative (longest) transcript per gene**. The exact files, source versions and checksums are in `MANIFEST_<species>.json` and `PROVENANCE_<species>.md` in [`docker/assets/`](docker/assets/); the QC report shows the reference build. Sources, citations and terms: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The references were built with the scripts in [`buildDB/`](buildDB/).

### Versions and reproducibility

The pipeline version is in [`docker/VERSION`](docker/VERSION); the git tag `v<version>` and the image tag `<version>` are the same. Each run records the version, the git commit, a fingerprint of the code and the reference build in `run_status.json` and in the QC report. What was tested: the same image (by digest), the same references and the same `--threads` gave identical results in repeated runs on one computer, between a local (conda) and a Docker run, and, from version 0.2.0, between two computers (a Windows laptop and an Intel Mac, two full samples: identical alignments, and identical tables, coverage files and assigned alignments of the read assignment). Compared were the count and statistics tables, the coverage files, `metrics_long.tsv` of the QC report and the records of the BAM files; not compared were logs, the HTML report (it holds the run date) and fastp's duplication estimate, which can differ in the sixth decimal between identical runs (the pipeline does not use it). Not tested: Apple-silicon Macs, and other thread counts than 8. Version 0.1.0 is reproducible on one computer, but between computers about 0.03 % of the reads could be assigned to another feature, all on pre-mRNA intron features, with an identical gene table (§9). The image digest of each release is in its release notes.

## 9. Known limitations

Read these before interpreting results.

- **RNA class depends on the tier order.** A read is assigned to the first reference it aligns to (§8), so the boundaries between RNA classes (e.g. piRNA, lncRNA, protein_coding) partly reflect the order of the references, not only biology.
- **Chance alignment of short reads.** Reads of 15–17 nt can align by chance to the large ncRNA and cDNA references before they reach the genome tier, and are then counted as lncRNA or protein_coding. Under a uniform-random model, a 15-nt read has about an 11 % chance of an exact match somewhere in the 62.5 Mb cDNA reference, and nearly 100 % if one mismatch is allowed; at 20 nt about 0.01 % and 0.7 %. These are model estimates, not measurements. `qc/metrics_long.tsv` has the number of 15–17 nt reads per RNA class (`assign.class.<class>.floor_reads`).
- **Strand-agnostic feature lookup in the genome tier.** Features of reads assigned in the genome tier are looked up without regard to strand, which affects the sense/antisense labels and counts of pre-mRNA and pre-tRNA features.
- **Equal-rank overlapping features.** When a read overlaps two features of the same rank, the feature listed first in the reference BED file wins (from version 0.2.0). The rule is technical, not biological: because the feature lookup in the genome tier ignores the strand (item above), it also decides whether such a read counts as sense (kept) or antisense (removed), so the later-listed feature loses these reads. Few reads are concerned: changing from the rule of version 0.1.0 moved 64 and 187 reads (0.005–0.006 %) of two full samples, all on pre-mRNA intron features, with an unchanged gene table. In version 0.1.0 the winner depended on the computer's processor (an unstable sort): between two computers about 0.03 % of the reads differed, all on pre-mRNA intron features, with an identical gene table.
- **Strand-specific libraries assumed.** The sense/antisense distinction assumes a strand-specific library.
- **Multimappers.** By default multimappers are left out of the counts. With `remove_multimap: False`, each multimapper is placed at one locus by a deterministic but arbitrary tie-break; no fractional or EM-based assignment is done. In `bedgraph_full/`, `multimap_depth` shows how much depth comes from multimappers.
- **UniVec and rRNA report one alignment** (Bowtie 1), so multimapping within these references is resolved by the aligner and does not appear in `multimap_depth` or in the multimapper counts.
- **One transcript per gene.** Isoform-level information (alternative starts and 3′ ends) is not represented.
- **Composition after removal.** The composition and detection numbers of the QC report are measured after multimappers and antisense reads are removed.
- **pre-mRNA names.** The pre-mRNA sequences are named `<gene>_intron`, but they span the whole gene on the genome: reads in exons that did not match the mature transcript count there too.
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
- **RNA class:** what a read is counted as (miRNA, tRNA, rRNA, lncRNA, protein_coding, …), from the RNA class table (§5).
- **RNA Type** (capital T): the combination of 5′ and 3′ end chemistry of an RNA, as in the method paper; here only for the spike-ins (RiboSpike-T1 to T4).
- **Label (biotype):** the fifth column of the reference BED, e.g. Ensembl's biotype; the RNA class table turns it into an RNA class.
- **Reference sequence:** one sequence of a reference BED (a transcript, a miRNA arm, a tRNA, a gene span, …); the rows of `*-readcounts.csv` and `*-featureclasses.csv` (column `feature`).
- **Feature:** a row of the gene table (`*-genecounts.csv`): a gene, a miRNA arm, a tRNA family, …
- **Tier:** one reference in the mapping order (§8, step 2).
- **Natural reads:** the kept reads without the spike-ins and the vector contaminants (UniVec).
- **Kept reads:** the mapped reads after multimappers and antisense reads are removed (default settings).
- **Tiered mapping:** alignment to an ordered stack of references; the first reference a read aligns to wins.
- **Sense / antisense:** a read on the same / opposite strand as its reference sequence.
- **Multimapper:** a read with two or more equally good best alignments.
- **Mature-exonic:** restricted to the processed transcript (no gene-span or precursor reads).
- **Coverage / bedgraph:** the number of reads covering each position of a reference sequence.
- **Unmapped:** a read that aligned to none of the references (the class `unmapped`); "unaligned" in the per-tier numbers means not aligned in that tier, so passed on to the next.
- **Spike-in (RiboSpike):** synthetic control RNAs of known sequence and length added to the sample.
- **Consistency check:** a comparison of the same quantity in the outputs of two pipeline steps.
