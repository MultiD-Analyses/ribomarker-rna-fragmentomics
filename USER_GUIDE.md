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
- **Species:** human (default) and mouse (`--species mouse`). The references of both are built into the image (§8, *References*); other species with your own references (§5.3). For mouse, the mouse references were checked on synthetic reads cut from them (every read lands on its own reference), but no result on real mouse data has been compared yet; the QC report of a mouse run says so (§9).
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
| any other columns | E.g. clinical data. Not used for processing. The sheet is copied unchanged into the results as `sample_sheet.csv`, and **all its columns are shown in the QC report** (section *Samples*): leave out what the report's readers should not see. In the report, *Colour samples by* colours the samples by any of these columns; `report_group` (§5.1) sets the first choice. |

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
| `--config FILE` | run settings (§5.1) |
| `--species NAME` | the reference species: `human` (default) or `mouse` with the built-in references (§9); with `--references`, the species name of your own files (§5.3) |
| `--references FOLDER` | your own references instead of the built-in ones (§5.3) |
| `--tiers FILE` | experimental: another order of the references (an order file; §5.2) |
| `--passes FILE` | experimental: several alignment passes with their own settings (§5.2) |
| `--rna-classes FILE` | an RNA class table replacing the default one (§5.4) |
| `--unlock` | remove a lock left by a run that was stopped hard (only if no other run uses `OUTPUT`) |
| `--demo` | run the demo (§2.4) instead of your data: `run --demo OUTPUT` |

### 4.4 Checks before the run starts

Before anything is computed, `run` checks the sample sheet, the FASTQ files (present and readable as gzip), the settings and the files given with `--config`, `--tiers`, `--passes` and `--rna-classes`, the references and the free disk space. If something is wrong it stops at once with exit status `2`, and its message says what is wrong and how to fix it. Some findings only give a warning and the run starts (for example less than 12 GB of memory, references without a `MANIFEST_<species>.json`, or references an order file leaves out); the warnings are printed before the run starts.

### 4.5 The end of a run

The run ends with one line:

- `RESULTS VALID: all <n> consistency checks equal …`: no consistency check of the QC report (§7) is `differs` or `not evaluable`;
- `RESULTS NOT VALID: …`: the pipeline stopped at a step (the line names the log), or a consistency check is `differs` or `not evaluable` (the line names the QC report);
- `RESULTS NOT CHECKED`: the QC report was switched off (`qc_report: False`, §5.1); exit status 0.

Exit status: `0` results valid (or not checked); `1` the pipeline failed or the results are not valid; `2` refused before starting; `130` interrupted.

The output of the mapping step may contain `[WARNING] Failed to launch x86-64-v3 version, staying with default`: a program looks for a build optimised for newer processors and uses its standard build instead. It does not affect the results; it appeared on every computer we tested.

### 4.6 Stopping and resuming

- **Ctrl-C** (also `docker stop`, or closing the terminal) stops a run cleanly and marks it `interrupted`.
- **Resume** by repeating the same command: finished steps are not redone; `--threads` may differ.
- A **finished** run is refused, so it is never redone by accident. Use another `--run-id` or another output folder for a new run.
- If a run was stopped hard (e.g. the computer restarted) and the next attempt says the output folder is locked, add `--unlock`.
- A run started with an earlier version cannot be resumed with this version (0.2.x with 0.3.0 or later, 0.3.x with 0.4.0): start it again in an empty folder.
- A run with `--tiers` or `--passes` is resumed only with the same files: an edited order or passes file is refused.

## 5. Advanced settings

Everything in this section is optional: the defaults are the validated settings. Four things can be changed, each with its own option:

| What you change | Option | Section |
|---|---|---|
| how reads are counted (multimappers, antisense reads) and what the report shows | `--config FILE` | §5.1 |
| the order of the references, and alignment passes with other settings (experimental) | `--tiers FILE`, `--passes FILE` | §5.2 |
| the references themselves | `--references FOLDER` | §5.3 |
| how the labels of the references become RNA classes | `--rna-classes FILE` | §5.4 |

**The order of the references and the rank of the RNA class table are two different things:**

- The **order of the references** (§5.2) decides **which reference** a read is assigned to. The references are tried one after another; the first one that aligns the read gets it, and the later ones never see it.
- The **rank** (`rank_in_tier`, §5.4) works **inside that one reference**: when a read overlaps several reference sequences there, or aligns equally well to several of them (a multimapper), the rank decides which one gets the read (for a multimapper after the strand: sense before antisense). A rank never moves a read to another reference.

For example, a read that aligns to a snoRNA in the reference `ncrna_snsno` and also to the genome is assigned to `ncrna_snsno`, because that reference comes first; a read that aligns only to the genome, at a pre-tRNA locus that lies inside the span of a gene, is counted as pre-tRNA (rank 4) and not as pre-mRNA (rank 1). The RNA class table has a column `tier` with the same reference names as the order file, but it only says which reference a label belongs to; it does not set the order.

### 5.1 Run settings

`--config FILE` takes a YAML file with any of these settings:

| Setting | Default | Effect |
|---|---|---|
| `remove_multimap` | `True` | Leave out reads that align equally well to several loci (multimappers). With `False`, each multimapper is counted at **one** of its loci, chosen by a deterministic but arbitrary rule (§9), and there is no gene table: the counts are in `*-readcounts.csv`, with multimappers as separate `<feature>_multimap` rows. |
| `remove_antisense` | `True` | Leave out reads on the strand opposite to their reference sequence. With `False` there is no gene table, as for `remove_multimap`; antisense reads keep their RNA class and are counted in `<feature>_antisense` rows of `*-readcounts.csv` (every class, since 0.4.0), and the strand-agnostic lookup of §9 matters more. Whatever the setting, `*-mappingstats.csv` counts the antisense reads (§6.2). |
| `multimap` | `100` | The maximum number of alignments reported per read (Bowtie 2 `-k`). |
| `keep_classes` | `all` | RNA classes written to the coverage files: `all`, or a comma-separated list of class names without spaces (e.g. `miRNA,tRNA`). The old name `keep_types` is still accepted, with a warning. |
| `count_ribospike` | `True` | Count the RiboSpike spike-in reads (QC report, *Spike-ins*). |
| `qc_report` | `True` | Write the QC report (§7). With `False`, the run ends with `RESULTS NOT CHECKED`. |
| `report_group` | none | A column of the sample sheet: the column the QC report first colours the samples by (*Colour samples by* in section *Samples* changes it; up to 8 groups, with more the samples are not coloured). |
| `tiers`, `passes`, `rna_classes` | none | The files of §5.2 and §5.4, without Docker only; with the starter or Docker give them as `--tiers`, `--passes` and `--rna-classes`, which mount the files. |

Example:

```yaml
remove_multimap: True
keep_classes: miRNA,tRNA,piRNA
report_group: group
```

Write `True` / `False` and numbers without quotes (`"False"` in quotes would count as true; `run` refuses it). Other settings (folders, threads, paired- or single-end, trimming) are set by `run` itself and are refused in `--config`. The settings `collapse`, `subsample_reads`, `trim_fastq`, `trim_only` and `dump_unmapped` were removed in 0.4.0: with their old default value (`False`, `0`, `True`, `False`, `False`) they are ignored with a warning, any other value stops the run.

### 5.2 Order of the references and passes (experimental)

By default each read is aligned once, through the 14 references in the order of the references folder (§8, step 2), each with its standard settings. Two optional files change this:

| You want | Give |
|---|---|
| another order of the references | `--tiers FILE` (an *order file*) |
| several alignment passes, each with its own settings | `--passes FILE` (a *passes file*) |
| both | both |

```bash
./ribomarker samples.csv results --tiers order.tsv                      # another order, one pass
./ribomarker samples.csv results --passes passes.tsv                    # the default order, several passes
./ribomarker samples.csv results --tiers order.tsv --passes passes.tsv  # both
```

Both are experimental: only the default settings are validated against known results (§9). Runs with either file are marked *Non-default settings (experimental)* at the top of the QC report (§7). A run with them is resumed only with the same files (§4.6).

#### The order file

A text file with the header `tier` and then one reference name per line, in the order to use. Every run writes the order it used to `tables/order.tsv` of its results folder; copy that file and change it. The default order (human and mouse):

```
tier
univec_bt1
ribospike
rrna_bt1
trnas
mirna
star
pre
pirna
ncrna_mt
ncrna_snsno
ncrna_other
ncrna_lnc
cdna
genome
```

For example, to align snRNA, snoRNA and scaRNA before piRNA, as in the method paper, move `ncrna_snsno` above `pirna` and `ncrna_mt` below it:

```
tier
univec_bt1
ribospike
rrna_bt1
trnas
mirna
star
pre
ncrna_snsno
pirna
ncrna_mt
ncrna_other
ncrna_lnc
cdna
genome
```

Rules:

- The names are the reference names of §8, step 2 (with your own references: the index name without the species, §5.3), each at most once.
- UniVec (`univec_bt1`) and the spike-ins (`ribospike`) come first, the genome last.
- A reference left out is not used (a warning). The spike-ins cannot be left out while the spike-in count is on (`count_ribospike`, §5.1).
- A tier list of version 0.3.0 (`aligner,index,retain`) is refused: write an order file.

#### The passes file

A pass is one walk through the references, in the order in effect, for the reads that are still unaligned; a read assigned in one pass is not aligned again. A passes file says which settings each reference gets in each pass. It is tab-separated, with a header line of six columns in this order, then one row per rule:

| Column | Values |
|---|---|
| `pass` | 1, 2, …: consecutive from 1 |
| `applies_to` | reference names separated by commas, or one group: `all`, `bowtie1`, `bowtie2` (the references of that aligner) or `other` (the references of that pass that no other row names). A group stands alone, not mixed with names. A reference is matched by at most one row per pass; a reference that no row of a pass matches is skipped in that pass. The spike-ins cannot be named |
| `max_mismatches` | `standard` (the default settings), or 0, 1, 2. Bowtie 1 references: at most this many mismatches over the whole read. Bowtie 2 references: a score budget of that many mismatches, in which a position on an N (in the read, or a read running onto the N padding at the end of a reference sequence) can take the place of a mismatch, and `2` also allows one gap of 1–2 bases; `0` is exact |
| `min_length`, `max_length` | read lengths in nt (at least 15; `max_length` at least `min_length`); other reads skip that reference in that pass and go on to the next one; reads that skip the last reference of the last pass are unmapped |
| `n_overhang` | only with `max_mismatches` 0, and for Bowtie 2 references only (ignored for Bowtie 1): at most 1 or 2 positions on an N and no mismatch; or `standard` |

Empty cells may be left empty or written `-` (= `standard`, no length limit). The order of the rows does not matter. Example: first exact alignments (at most one position on an N), then the default settings for the reads still unaligned:

```
pass	applies_to	max_mismatches	min_length	max_length	n_overhang
1	all	0	-	-	1
2	all	standard	-	-	standard
```

The spike-ins always run once, in pass 1, with their own settings. A budget is a limit, not a promise: Bowtie 2 finds an alignment through exact 20-nt pieces of the read (seeds), so a short read whose every 20-nt piece holds a mismatch or a gap may stay unaligned although it is within the budget (e.g. a 30-nt read with a gap in its middle is not aligned even with the standard settings). `run` warns when UniVec has no row in pass 1 (vector reads are then not removed first), when a reference runs in no pass, and when `max_mismatches` allows more than the standard settings (Bowtie 1: `2`; Bowtie 2: `2` when the row's `min_length` is under 19 nt, or not set).

#### The schedule

Before the run starts, `run` prints the *mapping schedule*: the references down, the passes across. Check it there. For the example passes file above:

```
Mapping schedule: 14 references, 2 passes; order: default (references folder); pass settings: passes file passes.tsv
  step  reference     aligner  pass 1               pass 2
     1  univec_bt1    bowtie   0 mismatches         standard
     2  ribospike     bowtie2  spike-in settings    -
     3  rrna_bt1      bowtie   0 mismatches         standard
     4  trnas         bowtie2  0 mismatches, <=1 N  standard
     5  mirna         bowtie2  0 mismatches, <=1 N  standard
     ...
    13  cdna          bowtie2  0 mismatches, <=1 N  standard
    14  genome        bowtie2  0 mismatches, <=1 N  standard, keeps unaligned
```

Each run saves the order (`tables/order.tsv`), the pass settings (`tables/passes.tsv`) and the schedule with the exact aligner options of every reference and pass (`tables/schedule.tsv`). With a passes file the BAM records carry the pass that aligned them (`XP` tag), and `qc/metrics_long.tsv` has the reads per reference and pass in place of `map.<index>.reads_in` and `.unaligned`: `map.p<pass>.<index>.arrived`, `.aligned`, `.unaligned`, `.skipped_window`, where `<index>` is the index name with the species (e.g. `map.p1.human-mirna.arrived`); the report's tier section shows the sums over the passes.

### 5.3 Own references

`--references FOLDER` replaces the built-in references with a folder of your own, e.g. with another miRNA or piRNA database or an extra spike-in. The starter mounts the folder read-only. The file names carry the species name given with `--species`: letters, digits, `_` and `-`, e.g. `--species rat` with `rat.dbpaths.csv` and `rat.bed`. For a species other than human and mouse, `run` warns and the QC report skips two species-specific items: the check that the BED holds no contigs named for another species (marked *not applicable*) and the red-blood-cell miRNA ratio. The folder holds:

| File | Content |
|---|---|
| `<species>.dbpaths.csv` | The references, one per line, **in their default order** (§5.2): `aligner,index,T`, e.g. `bowtie2,mirna/human-mirna,T`. See below. |
| Bowtie / Bowtie 2 indexes | One per line of the dbpaths file (`bowtie-build` or `bowtie2-build` on a FASTA file). |
| `<species>.bed` | The reference sequences: tab-separated, no header, six columns: sequence name (as in the index), start, end, name, label, strand. The RNA class table turns the label into an RNA class and a rank (§5.4); labels missing from the table get a warning and rank 10. |
| `<species>.rna_classes.tsv` | Optional: the RNA class table for these references (§5.4). Without it, the built-in table of the species is used, or the default table for other species. |
| `MANIFEST_<species>.json` | Optional. Without it, `run` warns that the QC report cannot identify the references. |

The three columns of `<species>.dbpaths.csv`:

- `aligner`: `bowtie2` (Bowtie 2, up to `multimap` alignments per read), `bowtie` (Bowtie 1, one alignment) or `ribospikebt2` (the spike-in settings), as in the built-in file. References with the aligner `ribospikebt2`, `ribospike` or `mind` are spike-ins: in an order file they come first, with UniVec; they run in pass 1 only and cannot be named in a passes file. Further forms (`bowtie1`, `mind`, `ribospike`, `bowtie2best`, `bowtie2sense`) exist for own references but are used by no built-in reference; an unknown aligner is refused.
- `index`: the index prefix, relative to the folder. The index name without the species is the reference name used in order and passes files and in the RNA class table (e.g. `mirna` for `mirna/human-mirna`).
- The third column must be `T` (every reference is kept and counted); `F`, which removed reads without counting them before 0.4.0, is refused.

`run` checks before the start that every index and the BED file exist. Where to start: a copy of the built-in references (§10, step 3), or the small demo references in the image (`/pipeline/demo/references`), changed one reference at a time; to rebuild references from their original sources, see [`buildDB/`](buildDB/README.md).

Caveats: the pipeline has been tested only with the built-in and the demo references. The gene table (`*-genecounts.csv`) and the QC report follow the RNA class table (§8, step 5): give your own labels rows in it, and check the report of the first runs with your own references.

### 5.4 RNA class table

The RNA class table turns the labels of the reference BED (its fifth column, e.g. Ensembl's biotype) into RNA classes, and gives each label the rank that decides between overlapping reference sequences of one reference (see the start of §5). It is a tab-separated file, `<species>.rna_classes.tsv`, with one row per label; a label found in several references may have one row per reference (in this version these rows must agree on class and rank, since a read's class does not depend on its reference). Its columns:

| Column | Content |
|---|---|
| `source_label` | the label in the BED, or a label the pipeline sets itself (`unannotated`, `unmapped`, the spike-ins, …) |
| `tier` | the reference the label comes from (its reference name, e.g. `ncrna_snsno`), or `*` for all. It does not set the order of the references (§5.2). References the run does not use are named in a warning (e.g. the old name `ncrna` of 0.3.0, or a reference an order file leaves out) |
| `class_id` | the RNA class written into the result tables (no spaces) |
| `display_name` | the class's name in the QC report |
| `group` | its group in the report's tier and composition plots |
| `display_order` | the order of the classes in the report |
| `rank_in_tier` | which of overlapping reference sequences of one reference a read gets, and at which of its equally good loci a multimapper is counted (after the strand: sense before antisense): the higher rank wins. In the built-in tables structural and small RNAs rank high and protein_coding and pre-mRNA lowest |
| `positions` | `yes`: reads are split by their place in the reference sequence (`_5p`, `_3p`, `_internal`, `_whole`; §8, step 3) |
| `gene_table` | `yes`: in the gene table; `drop`: precursor or gene-span classes, left out; `exclude`: not a gene (spike-ins, vector contaminants, unannotated, unmapped) |
| `counts_as_rrna` | `yes`: left out of the "without rRNA" numbers of the QC report |
| `description` | one line for the report |

The built-in tables are in [`docker/tables/`](docker/tables/). Every run copies the class table it used into `tables/rna_classes.tsv` of its results folder, and the QC report shows it (section *Tiers and RNA classes*; file and checksum under *Provenance*). To change it, start from this copy and give the changed file with `--rna-classes FILE`:

```bash
./ribomarker ~/data/fastq/samples.csv ~/data/results --rna-classes my_classes.tsv
```

Which class table is used (the first found): `--rna-classes` (or the setting `rna_classes`); a `<species>.rna_classes.tsv` in the references folder; the built-in table of the species; `docker/tables/default.rna_classes.tsv`. A table replaces the built-in one as a whole; rows for the labels the pipeline sets itself are added from the built-in table when missing. BED labels missing from the table become a class of the same name, with rank 10 (the highest built-in rank), the group "lncRNA and other ncRNA", the gene table, and a warning. Since 0.4.0 the table has no row `antisense` (antisense reads keep their class); a table copied from 0.3.0 that still has it gets a warning, and the row is not used.

Changing the ranks changes which reads are assigned to which reference sequence. The pipeline has been tested with the built-in tables only; check the QC report of the first runs with your own.

## 6. Outputs

All results are in `OUTPUT/<run id>/` (default `OUTPUT/ribomarker/`); the tables are named `<run id>-….csv`. Sample names are those of the sample sheet (§3).

### 6.1 Counts and coverage

| File | Content |
|---|---|
| `*-genecounts.csv` | **The gene table: counts per feature** (features × samples; column `gene`): one row per gene, miRNA arm, tRNA family, …, from mature-exonic reads (§8, step 5). Unique and sense reads only (with the default settings). Not written when `remove_multimap` or `remove_antisense` is `False`. |
| `*-genelengths.csv` | The length (`length_nt`), RNA class (`rna_class`) and gene symbols (below) of every feature of the gene table. |
| `*-readcounts.csv` | Counts per reference sequence (column `feature`), before the collapse into the gene table. tRNA, snRNA and snoRNA sequences have one row per position (`_5p`, `_3p`, `_internal`, `_whole`). With `remove_multimap` or `remove_antisense` `False`, multimappers and antisense reads have rows of their own (`_multimap`, `_antisense`). |
| `*-featureclasses.csv` | The RNA class (`rna_class`), length (`length_nt`) and gene symbols (below) of every reference sequence of `*-readcounts.csv`. |
| `*-featuredistinctseqs.csv` | The number of distinct read sequences per reference sequence. |
| `*-classcounts.csv` | Counts per RNA class (`rna_class`), including the spike-ins, vector contaminants (`UniVec`), `genome_unannotated` and `unmapped`. |
| `*-classlengths.csv` | Read-length distribution per RNA class, from 0 nt to the longest read of the run (at least 150 nt). |
| `*-classdistinctseqs.csv` | The number of distinct read sequences per RNA class. |
| `bedgraph/<sample>.bedgraph` | Per-sequence coverage consistent with the counts (the same reads; each multimapper at its one assigned locus). Columns: `feature start end depth rna_class`. |
| `bedgraph_full/<sample>.bedgraph` | Coverage from every equally best locus of every read, with `multimap_depth` (how much of the depth comes from multimappers) and `strand`. Columns: `feature start end depth multimap_depth strand rna_class`. |
| `alignments/<sample>.bam` | The assigned alignments; the assigned reference sequence is in the `CS` tag. |
| `tables/` | What the run used (§5.2, §5.4): the order of the references (`order.tsv`, an order file), their lines of `dbpaths.csv` in that order (`tiers.csv`), the mapping schedule (`schedule.tsv`: one row per reference and pass, with the aligner options), which tier each reference sequence belongs to (`contig_tiers.tsv`), the RNA class table (`rna_classes.tsv`) and, with `--passes`, the passes file (`passes.tsv`). |

**Gene symbols.** `*-featureclasses.csv` and `*-genelengths.csv` have the columns `hgnc_id` and `hgnc_symbol` (human) or `mgi_id` and `mgi_symbol` (mouse), and `mirbase_name` (the miRBase name of a miRNA arm, star arm or precursor, e.g. `hsa-miR-21-5p`). They come from `docker/tables/<species>.gene_symbols.tsv` (sources, versions and checksums in its `.provenance.md`), joined by Ensembl gene id (no second try by gene name, so a symbol is never taken from another gene of the same name); pre-mRNA through the gene names of the Ensembl GTF, miRNA arms and precursors through MirGeneDB's miRBase accessions, tRNAs through GtRNAdb ids. A feature that spans several genes (e.g. a tRNA family) lists them separated by `;`; in `*-genelengths.csv` a feature gets the first non-empty value of its reference sequences. How many have a value (built-in references): human protein_coding and pre-mRNA 97 %, snRNA > 99 %, miRNA arms and precursors 98 % (miRBase names: 96 % of the arms, 85 % of the star arms), tRNA 98 %, snoRNA 51 %, lncRNA 19 % (most human lncRNAs have no approved symbol); mouse protein_coding and pre-mRNA > 99 %, lncRNA 98 %, miRNA and tRNA none (MGI has no miRBase or GtRNAdb cross-references; the miRBase names are filled in, 98–100 %). piRNA, the rRNA copies of miRTrace, spike-ins and vector contaminants have no symbol. With your own references (§5.3) the columns stay empty, and the run warns, unless the table was built for your BED (`buildDB/make_gene_symbols.py`).

**Joining the tables.** The count tables hold counts only; the annotations are in separate tables with the same first column. Join `*-genecounts.csv` with `*-genelengths.csv` on `gene`, and `*-readcounts.csv` with `*-featureclasses.csv` on `feature`, e.g. in R `merge(counts, lengths, by = "gene")`, or in Python `counts.merge(lengths, on="gene")` (pandas).

### 6.2 Statistics and logs

| File | Content |
|---|---|
| `*-mappingstats.csv` | Per sample: `Mapped`, `Unmapped`, `Unique`, `Multi` (multimappers), `Antisense` (unique reads on the opposite strand of their feature), `AntisenseMulti` (multimappers whose assigned locus is antisense; since 0.4.0), `Total`. The antisense columns count the reads whatever `remove_antisense` says. |
| `*-trimstats.csv` | Trimming statistics per sample (adapter dimers, trimmed reads, merged pairs) |
| `ribospike/` | Spike-in counts by length and RNA Type (`*-spikelengths.csv`) and the origin of every read: natural, spike-in or unmapped (`*-readorigin.csv`) |
| `trimmed/logs/` | cutadapt and (paired-end) fastp reports |
| `alignments/` | The merged alignments of each sample (`<sample>.merged.bam`: every reported alignment of every reference and pass, with the unmapped reads), the assigned alignments (`<sample>.bam`, §6.1), the aligners' logs per reference and pass (`alignments/logs/<sample>.p<pass>.<index>.log`) and, per sample, which references ran in each pass and how many reads a length window let skip (`alignments/logs/<sample>.cells.tsv`). The alignments per reference and pass are merged and then deleted. |
| `sample_sheet.csv` | Your sample sheet, unchanged (with any extra columns) |
| `samples.csv`, `run_config.yaml` | What the pipeline was given: its own sample list and the settings |
| `run.log` | The full log of the run |
| `run_status.json` | Verdict, start and end times, pipeline version, code fingerprint and git commit, the order, pass settings and RNA class table used (`tables`); in Docker also the memory peak |

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

1. **Samples:** the columns of the sample sheet (the FASTQ file names are hidden at first; *Configure columns* shows them), and *Colour samples by*. Its first line gives the result of the consistency checks and the number of warnings; for a run with `--tiers` or `--passes` it starts with *Non-default settings (experimental)* and names them.
2. **How to read this report:** the settings above the plots, the toolbox, and a glossary.
3. **General Statistics:** one row per sample with the headline numbers, including the mean read length of the natural reads; below the heading, the checks that are not `equal` and the warnings (only if there are any).
4. **Trimming and merging: read fate:** where the input reads (pairs) went: adapter dimers, inserts of 1–14 nt, removed by fastp filters and not merged (both paired-end only), into mapping. The categories add up to the input.
5. **Mapping per reference tier:** reads assigned in each mapping tier (§8), before multimappers and antisense reads are removed; with `--passes`, summed over the passes. The genome tier is split into the unannotated reads and the rest (pre-mRNA, pre-tRNA).
6. **Composition by RNA class:** kept reads per RNA class after that removal, grouped as in the RNA class table; further tabs without spike-ins, and without rRNA and spike-ins.
7. **Multimappers and antisense reads:** the mapped reads split into kept reads, reads removed as multimappers and reads removed as antisense (with `remove_multimap` on, antisense multimappers count as multimappers), counted from `*-mappingstats.csv`; below it, per tier, the share of reads with a second alignment.
8. **Read length by RNA class**, in two modes (setting *Values*): **% of natural reads**, the reads of a class at a length as a share of all natural reads of the sample (all classes and lengths together add up to 100 %), and **% within class**, the length distribution of each class on its own. Natural reads are the kept reads without spike-ins and vector contaminants. The method paper shows such sequencing length profiles without stating their denominator (Clark et al. 2026); here it is the natural reads, counted after multimappers and antisense reads are removed (the default settings). One button per main class; the other classes are together under "Other RNAs" (the description lists them).
9. **Spike-ins:** share of reads, fragments, length fidelity, end-base composition of exact-length reads (with its chance level), length offset and recovery per designed length. Without the amounts added, recovery can only be compared between samples and lengths.
10. **Detection and counts:** library size (with and without rRNA), features detected with more than 0, 5 and 10 reads (with a plot per RNA class), the ratio of two red-blood-cell-enriched miRNAs (miR-451a / miR-23a-3p; descriptive, not comparable to qPCR cut-offs without calibration), and the 20 most abundant features per sample, overall and per RNA class.
11. **Tiers and RNA classes:** the order of the references (`tiers.csv`), the mapping schedule (references down, passes across) and the RNA class table of the run (§5.2, §5.4).

These are followed by MultiQC's own fastp (paired-end) and cutadapt sections, **Provenance** (run, library parameters, reference species check, reference build, the tables, settings, and how to cite the pipeline and the method) and **Software versions** (the pipeline version and every tool version).

*Colour samples by*, above the sample table, colours the samples by the values of a sample-sheet column in every plot and table (at most 8 groups; `report_group`, §5.1, sets the first choice); a click on a colour in the key next to it changes that colour. The choice and the colours are not kept when the report is opened again, and it replaces colours set by hand in the toolbox's *Highlight* tab. Plots with one line per sample have no legend in MultiQC's own sections: hover over a line to see its sample, or use the toolbox on the right (*Highlight samples*, *Show / hide samples*).

Two kinds of marks need attention:

- **Consistency checks** compare the same quantity in the outputs of two steps (e.g. reads written by trimming = reads entering mapping). `equal` shows that read totals are conserved and the files are complete; it does not show that reads were classified correctly. `differs` means two outputs disagree, which points to a pipeline problem rather than a sample problem: please open an issue. `not evaluable` means a number the check needs could not be read; like `differs`, it makes the run `RESULTS NOT VALID`. `not applicable` is normal for checks that do not apply to the layout (e.g. the fastp check in single-end runs). Every check is in `qc/consistency_checks.tsv`.
- **`low n` / `not evaluable`:** fewer than 100 reads, or none, are behind a share, so it is not reliable. The report shows these flags in the spike-in and counts tables; all flags are in `qc/qc_metrics.json`.

### Sending the report

For help with a run, send `qc/report_technical.html` (and `qc/qc_metrics.json` if asked). Absolute paths are removed from the report. It contains the sample names, the FASTQ file names and **all columns of your sample sheet** (so does `qc/qc_metrics.json`), so check these before sending it outside your organisation. It contains aggregate numbers, not read sequences (the pipeline does not switch on fastp's over-represented-sequence analysis).

## 8. How the pipeline works

The workflow is run by Snakemake; each step is a Python script in `docker/scripts/`. The steps, in order:

1. **Trimming.** cutadapt removes the 3′ adapter (`TGGAATTCTCGGGTGCCAAGG`; R2: `TGATCGTCGGACT`) and the one extra base that the library protocol adds at the insert's 5′ end, and drops inserts shorter than 15 nt (reported as adapter dimers and short inserts). For paired-end reads, fastp then merges R1 and R2 into one insert (overlap of at least 10 nt, at most 2 differences) and filters reads with its default criteria (low quality, too many N); pairs that cannot be merged are not used. No base-quality trimming is done.
2. **Tiered mapping.** Reads are aligned to a stack of references **in a fixed order**, and each read is assigned to the *first* reference it aligns to:
   `UniVec → RiboSpike → rRNA → tRNA → miRNA → miRNA star → miRNA precursor → piRNA → ncRNA (mitochondrial) → ncRNA (snRNA, snoRNA, scaRNA) → ncRNA (other) → lncRNA → cDNA → genome`.
   Contaminants, spike-ins and rRNA are claimed first; reads that do not align go on to the next reference. The reference names (used in order and passes files, and shown in `tables/order.tsv`) are `univec_bt1`, `ribospike`, `rrna_bt1`, `trnas`, `mirna`, `star`, `pre`, `pirna`, `ncrna_mt`, `ncrna_snsno`, `ncrna_other`, `ncrna_lnc`, `cdna`, `genome`. Since 0.4.0 the Ensembl ncRNA reference is four references (Mt_tRNA and Mt_rRNA; snRNA, snoRNA and scaRNA; misc_RNA and nuclear rRNA; lncRNA), so a read that aligns to two of them is assigned to the first instead of being a multimapper. The order is that of `<species>.dbpaths.csv`, or of an order file (§5.2); with a passes file the reads go through the references once per pass. The miRNA star arms are a reference of their own but are counted as the class miRNA. UniVec and rRNA are aligned with Bowtie 1, reporting one alignment; the other tiers with Bowtie 2 (up to `multimap` alignments per read).
3. **Assignment.** For each read the pipeline resolves, in order: (a) the reference (the tier order above); within it, the best-scoring alignments: a read with one best alignment is unique, a read with several is a multimapper; (b) the overlapping reference sequence, by the rank of its label (`rank_in_tier` of the class table), in which structural RNAs outrank others and protein_coding ranks lowest; a read overlaps a sequence when at least 2 of its bases lie in it; (c) the locus of a multimapper, by a deterministic tie-break: sense before antisense before no sequence, then the higher rank, then the position. A read is antisense when its (assigned) alignment lies on the opposite strand of its sequence; a read on the genome without a sequence is not antisense. Multimappers and antisense reads are then removed (§5.1); with `remove_antisense: False` an antisense read keeps its class. Reads of classes with `positions` (tRNA, snRNA, snoRNA) are split by their place in the sequence: `_5p` or `_3p` when the read lies within 10 nt of that end only, `_whole` when it lies within 10 nt of both ends, `_internal` otherwise; the 5′ end is on the sequence's strand. Except for `_whole`, this follows the method paper (Clark et al. 2026, Methods).
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

- **RNA class depends on the tier order.** A read is assigned to the first reference it aligns to (§8), so the boundaries between RNA classes (e.g. piRNA, lncRNA, protein_coding) partly reflect the order of the references, not only biology. The best alignment wins only within a reference: an earlier reference can win with a worse alignment than a later one (since the ncRNA split of 0.4.0 also between the four ncRNA references).
- **Order and passes are experimental.** Other orders (`--tiers`) and passes (`--passes`) are checked for consistency (the reads add up, every read is assigned once), but only the default settings are validated against known results; the QC report marks such runs.
- **Chance alignment of short reads.** Reads of 15–17 nt can align by chance to the large ncRNA and cDNA references before they reach the genome tier, and are then counted as lncRNA or protein_coding. Under a uniform-random model, a 15-nt read has about an 11 % chance of an exact match somewhere in the 62.5 Mb cDNA reference, and nearly 100 % if one mismatch is allowed; at 20 nt about 0.01 % and 0.7 %. These are model estimates, not measurements. `qc/metrics_long.tsv` has the number of 15–17 nt reads per RNA class (`assign.class.<class>.floor_reads`).
- **Strand-agnostic feature lookup in the genome tier.** Features of reads assigned in the genome tier are looked up without regard to strand, which affects the sense/antisense labels and counts of pre-mRNA and pre-tRNA features.
- **Equal-rank overlapping features.** When a read overlaps two features of the same rank, the feature listed first in the reference BED file wins (from version 0.2.0). The rule is technical, not biological: because the feature lookup in the genome tier ignores the strand (item above), it also decides whether such a read counts as sense (kept) or antisense (removed), so the later-listed feature loses these reads. Few reads are concerned: changing from the rule of version 0.1.0 moved 64 and 187 reads (0.005–0.006 %) of two full samples, all on pre-mRNA intron features, with an unchanged gene table. In version 0.1.0 the winner depended on the computer's processor (an unstable sort): between two computers about 0.03 % of the reads differed, all on pre-mRNA intron features, with an identical gene table.
- **Strand-specific libraries assumed.** The sense/antisense distinction assumes a strand-specific library.
- **Multimappers.** By default multimappers are left out of the counts. With `remove_multimap: False`, each multimapper is placed at one locus by a deterministic but arbitrary tie-break; no fractional or EM-based assignment is done. In `bedgraph_full/`, `multimap_depth` shows how much depth comes from multimappers.
- **UniVec and rRNA report one alignment** (Bowtie 1), so multimapping within these references is resolved by the aligner and does not appear in `multimap_depth` or in the multimapper counts. Which copy is reported is not guaranteed to be the best one (`-k 1` without `--best`), and can differ between settings.
- **Borderline reads at the ends of reference sequences.** A read one or two bases longer than a short sequence (e.g. a miRNA arm or a snoRNA) aligns only by paying for a position on the N padding, and whether Bowtie 2 finds such an alignment can depend on the reference searched (a few reads of SNORD14C were no longer aligned after the ncRNA split: at most 38 of 4.6 million).
- **One transcript per gene.** Isoform-level information (alternative starts and 3′ ends) is not represented.
- **Composition after removal.** The composition and detection numbers of the QC report are measured after multimappers and antisense reads are removed.
- **pre-mRNA names.** The pre-mRNA sequences are named `<gene>_intron`, but they span the whole gene on the genome: reads in exons that did not match the mature transcript count there too.
- **Mouse not yet compared with real data.** For mouse, the mouse references were checked on synthetic reads cut from them (every read lands on its own reference), but no result on real mouse data has been compared yet; the QC report of a mouse run says so. Other species: only with your own references (§5.3), not tested by us.

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
- **RNA class:** what a read is counted as (miRNA, tRNA, rRNA, lncRNA, protein_coding, …), from the RNA class table (§5.4).
- **RNA Type** (capital T): the combination of 5′ and 3′ end chemistry of an RNA, as in the method paper; here only for the spike-ins (RiboSpike-T1 to T4).
- **Label (biotype):** the fifth column of the reference BED, e.g. Ensembl's biotype; the RNA class table turns it into an RNA class.
- **Reference sequence:** one sequence of a reference BED (a transcript, a miRNA arm, a tRNA, a gene span, …); the rows of `*-readcounts.csv` and `*-featureclasses.csv` (column `feature`).
- **Feature:** a row of the gene table (`*-genecounts.csv`): a gene, a miRNA arm, a tRNA family, …
- **Tier (reference):** one reference in the mapping order (§8, step 2), named by its index without the species (e.g. `ncrna_snsno`).
- **Pass:** one round of alignment through the references, with its own settings (§5.2); by default there is one.
- **Natural reads:** the kept reads without the spike-ins and the vector contaminants (UniVec).
- **Kept reads:** the mapped reads after multimappers and antisense reads are removed (default settings).
- **Tiered mapping:** alignment to an ordered stack of references; the first reference a read aligns to wins.
- **Sense / antisense:** a read on the same / opposite strand as its reference sequence; for a multimapper, at its assigned locus.
- **Multimapper:** a read with two or more equally good best alignments.
- **Mature-exonic:** restricted to the processed transcript (no gene-span or precursor reads).
- **Coverage / bedgraph:** the number of reads covering each position of a reference sequence.
- **Unmapped:** a read that aligned to none of the references (the class `unmapped`); "unaligned" in the per-tier numbers means not aligned in that tier, so passed on to the next.
- **Spike-in (RiboSpike):** synthetic control RNAs of known sequence and length added to the sample.
- **Consistency check:** a comparison of the same quantity in the outputs of two pipeline steps.
