# RiboMarker: a pipeline for RNA fragmentomics

The RiboMarker pipeline turns small-RNA sequencing reads (FASTQ files) into:

- **counts** per gene and per RNA type,
- **per-transcript coverage** files (bedgraph),
- a **technical QC report** (one HTML file) with every step's numbers and consistency checks.

It is made for libraries prepared with the RiboMarker small-RNA protocol of RealSeq Biosciences and for **human** or **mouse** samples (for mouse, see the [known limitations](USER_GUIDE.md#9-known-limitations)). It runs in Docker; one command starts it ([installation](USER_GUIDE.md#2-installation)). The pipeline stops at counts and coverage: it applies no normalization and no statistics.

Full documentation: **[USER_GUIDE.md](USER_GUIDE.md)**.

## What the pipeline does

Every read gets **one RNA type**, by **prioritized (tiered) mapping**: the reads are aligned to the references below **in this order**, and each read is assigned to the **first** reference it aligns to. Only the reads that do not align go on to the next reference.

| Priority | Reference | Source |
|---:|---|---|
| 1 | UniVec: vector and adapter sequences (contaminants) | NCBI |
| 2 | RiboSpike spike-ins | RealSeq Biosciences |
| 3 | rRNA | miRTrace (NCBI, SILVA) |
| 4 | tRNA | GtRNAdb |
| 5 | mature miRNA | MirGeneDB |
| 6 | miRNA star | MirGeneDB |
| 7 | miRNA precursor | MirGeneDB |
| 8 | piRNA | piRBase |
| 9 | other non-coding RNA (e.g. lncRNA, snoRNA, snRNA) | Ensembl ncRNA |
| 10 | mRNA | Ensembl cDNA |
| 11 | genome (the rest, e.g. introns) | Ensembl |

So contaminants, spike-ins and rRNA are claimed first, the small, well-defined RNA classes come before the large transcript references, and the genome comes last. Within a reference, a read that overlaps several annotated features gets the feature of the highest-ranked biotype (structural RNAs such as tRNA and rRNA first, mRNA last). The RNA types therefore partly reflect this order, not only biology ([known limitations](USER_GUIDE.md#9-known-limitations)).

The steps of a run:

1. **Trimming:** the adapters and the extra 5′ base of the protocol are removed, inserts shorter than 15 nt are dropped, and paired-end reads are merged into one insert.
2. **Prioritized mapping** to the references above (Bowtie, Bowtie 2).
3. **Counting** per gene, per reference feature and per RNA type; reads that align equally well to several places (multimappers) and antisense reads are left out by default.
4. **Coverage:** per-transcript bedgraph files.
5. **QC:** a technical report with every step's numbers and consistency checks between the steps.

Details: [How the pipeline works](USER_GUIDE.md#8-how-the-pipeline-works).

## Quick start

**1. Install Docker** and start it. On Windows, use Docker Desktop with the WSL2 backend.

**2. Get the starter.** Download this repository (*Code → Download ZIP*, or the source archive of a release) and unzip it. The starter scripts are at its top level: `ribomarker` (Linux, macOS, Windows WSL shell) and `ribomarker.cmd` with `ribomarker_windows.ps1` (Windows PowerShell or command prompt; keep the two together). The starter writes the `docker run` command for you; Docker downloads the image the first time.

**3. Test the installation** (a few minutes; synthetic reads whose results are known):

```bash
./ribomarker --demo results/demo            # Linux, macOS, Windows WSL shell
```
```powershell
.\ribomarker --demo results\demo            # Windows PowerShell
```

It ends with `DEMO PASSED: all … known results reproduced`.

**4. Run your data.** Put the FASTQ files (`.fastq.gz`) into one folder, together with a sample sheet, a CSV file with one row per sample (more columns, e.g. clinical data, are allowed and copied to the results):

```
sample,fastq_file,fastq_file2
Healthy_1,H1_R1.fastq.gz,H1_R2.fastq.gz
Patient_7,P7_R1.fastq.gz,P7_R2.fastq.gz
```

Leave out `fastq_file2` for single-end reads. Then:

```bash
./ribomarker ~/data/fastq/samples.csv ~/data/results    # Linux, macOS, Windows WSL shell
```
```powershell
.\ribomarker C:\data\fastq\samples.csv C:\data\results  # Windows PowerShell
```

The FASTQ folder is only read; the results go into the second folder. The sample names of the sheet are used in all results ([rules](USER_GUIDE.md#3-inputs)). The run ends with one line, `RESULTS VALID` or `RESULTS NOT VALID`, and the path of the QC report. `./ribomarker --help` lists the options. If `./ribomarker` says *permission denied*, run `chmod +x ribomarker` once.

## What you get

In `results/ribomarker/` (the run id; change it with `--run-id`):

| File | Content |
|---|---|
| `qc/report_technical.html` | The technical QC report. Open this first. |
| `ribomarker-genecounts.csv` | Counts per gene and sample: integer counts, suitable for count-based differential-expression tools |
| `ribomarker-typecounts.csv`, `ribomarker-readcounts.csv` | Counts per RNA type and per reference feature |
| `bedgraph/`, `bedgraph_full/` | Per-transcript coverage per sample |
| `sample_sheet.csv` | Your sample sheet, unchanged |
| `run_status.json`, `run.log` | Verdict, times, pipeline version, and the full log |

All outputs are described in the [user guide](USER_GUIDE.md#6-outputs).

## Advanced settings

The defaults are the validated settings; most runs need no options. A few options go on the command line, for example:

```bash
./ribomarker ~/data/fastq/samples.csv ~/data/results --species mouse --threads 16 --config settings.yaml
```

| Option | Meaning |
|---|---|
| `--species mouse` | mouse references instead of human |
| `--threads N` | CPU cores to use (default 8) |
| `--run-id NAME` | name of the results folder (default `ribomarker`) |
| `--config FILE` | a YAML file with advanced settings |

The config file changes how reads are counted and what is written, for example which RNA types get coverage files:

```yaml
keep_types: miRNA,tRNA,piRNA
```

All options and settings, with their effects: [user guide §4.3 and §5](USER_GUIDE.md#5-advanced-settings). `./ribomarker --help` lists the options.

## Own references

The built-in references can be replaced by your own, for example another miRNA or piRNA database, an extra spike-in, references of another species (`--species rat` with your own `rat` files), or another priority order (the order is the line order of one small file):

```bash
./ribomarker ~/data/fastq/samples.csv ~/data/results --references ~/my_references
```

The folder holds the Bowtie indexes, a list of them in priority order, and a BED file with the features and their RNA types. The easiest start is a copy of the built-in references, changed one reference at a time. Format and caveats: [user guide, *Own references*](USER_GUIDE.md#own-references); building references from their original sources: [`buildDB/`](buildDB/README.md). Own references have not been validated by us: check the QC report of the run.

## Citing

If you use the pipeline, please cite this repository, with the version you used (shown in the QC report and in `run_status.json`):

> Smelik M. RiboMarker: a pipeline for RNA fragmentomics, version X.Y.Z. MultiD Analyses AB; 2026. https://github.com/MultiD-Analyses/ribomarker-rna-fragmentomics

The *Cite this repository* button on GitHub (from [CITATION.cff](CITATION.cff)) gives this citation in APA and BibTeX. Please also cite the method:

> Clark RC, Manning AC, Howard JM, Barberán-Soler S, Kazakov SA. A method of comprehensive sequencing analysis of the small RNA fragmentome (RiboMarker). *RNA* 2026;32(7):1156–1181. https://doi.org/10.1261/rna.080861.125

The reference databases built into the image have their own citations: see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Questions and problems

Please open an issue. For a problem with a run, attach `qc/report_technical.html` if you can (see [what it contains](USER_GUIDE.md#sending-the-report)). Contributions: [CONTRIBUTING.md](CONTRIBUTING.md).

## Licence

MIT, see [LICENSE](LICENSE). The logos and the names RiboMarker, RealSeq and MultiD are not covered by the licence ([NOTICE](NOTICE)). Third-party software and data: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
