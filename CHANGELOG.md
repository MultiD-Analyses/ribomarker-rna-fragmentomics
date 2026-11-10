# Changelog

All notable changes are listed here. Versions follow `MAJOR.MINOR.PATCH`; while the version is 0.x, a change that alters results (including a new reference build) raises the minor version, other changes the patch version. The git tag `v<version>` and the image tag `<version>` belong to the same commit.

## 0.4.0 – 2026-11-09

Results change through the ncRNA split (see *Changed results*); with the default settings nothing else changes. Two experimental controls are new: another order of the references (`--tiers`) and several alignment passes (`--passes`). Several settings and one output file are gone (see *Removed*). Runs started with 0.3.x cannot be resumed with 0.4.0; start them again in an empty folder.

### Changed results

- **The Ensembl ncRNA reference is four references**, aligned in this order where the one reference was: `ncrna_mt` (Mt_tRNA, Mt_rRNA), `ncrna_snsno` (snRNA, snoRNA, scaRNA), `ncrna_other` (misc_RNA, nuclear rRNA loci) and `ncrna_lnc` (lncRNA); the default order of the references is otherwise unchanged. Every record is unchanged; only the indexes are split. The effects, measured on two full paired-end samples (1.6 and 4.6 million reads):
  - a read that aligned equally well to records of two of the four references was a multimapper and was removed by default; it is now assigned to the first of them and kept: 2,830 and 8,648 reads, most of them in `ncrna_snsno` (in the larger sample 6,877, the rest in `ncrna_other` 1,510 and `ncrna_mt` 261);
  - the other way round, 149 and 312 kept reads are now removed: almost all of them former lncRNA reads (308 of the 312) whose alignments now fall first in an earlier ncRNA reference, where they are multimappers or antisense;
  - the best alignment wins only within a reference, so a read can now be assigned to an earlier ncRNA reference with a worse alignment than one in a later reference: 72 and 191 reads;
  - per class (kept reads, forecast from the 0.3.0 alignments; the observed moves in the larger sample are close to it): snRNA, snoRNA and scaRNA about +24 % and +32 %, misc_RNA and nuclear rRNA loci about ×3.5 to ×4 (from a small base), lncRNA −2 %;
  - in all, 28,138 and 70,255 reads (1.8 % and 1.5 %) changed their alignments or assignment, most of them reads that are removed before and after (21,512 and 51,853: multimappers with other loci reported); a few borderline reads at the end of SNORD14C are no longer aligned (up to 38 of 4.6 million).
- **With `remove_antisense: False`**, antisense reads keep their RNA class (before: the class `antisense`), and antisense reads of protein_coding and pre-mRNA features get `_antisense` rows in `*-readcounts.csv` and `*-featuredistinctseqs.csv`, like every other class (before: they were counted with the sense feature). The class `antisense` is gone from the class tables; a table copied from 0.3.0 that still has the row gets a warning. A class table of 0.3.0 that names the reference `ncrna` still works; the warning names the four new references. The tier metrics `map.<species>-ncrna.*` in `qc/metrics_long.tsv` are now `map.<species>-ncrna_mt.*`, `-ncrna_snsno.*`, `-ncrna_other.*` and `-ncrna_lnc.*`.
- **Snakefile with a hand-written configuration only:** a configuration without a `collapse` key collapsed identical reads before mapping; read collapsing is removed, so such runs change. `run` never collapsed.

### New

- **Order of the references (`--tiers FILE`, experimental):** an *order file*, one column `tier` with the reference names in the order to use (UniVec and the spike-ins first, the genome last; references left out are not used). The 0.3.0 tier-list format (`aligner,index,retain`) is refused. Every run writes the order it used to `tables/order.tsv`, to copy and change.
- **Passes (`--passes FILE`, experimental):** several alignment passes, each through the references in the order in effect, for the reads still unaligned, with rules per reference: `max_mismatches` (`standard`, 0, 1, 2), a read-length window (`min_length`, `max_length`) and `n_overhang` (at most 1 or 2 positions on an N, with 0 mismatches). The spike-ins run once, in pass 1, with their own settings. The BAM records carry the pass that aligned them (`XP` tag); `qc/metrics_long.tsv` has the reads per reference and pass (`map.p<pass>.<index>.arrived`, `.aligned`, `.unaligned`, `.skipped_window`, in place of `map.<index>.reads_in` and `.unaligned`), and the report's tier section shows the sums over the passes; `tables/passes.tsv` is a copy of the passes file. Idea and first implementation: Ruslan Klassen (pull request #1); 0.4.0 keeps reads under 20 nt in play instead of dropping them, and runs the spike-ins once.
- **The mapping schedule** (references down, passes across) is printed before every run, saved in `tables/schedule.tsv` (one row per reference and pass, with the aligner options) and shown in the QC report (*Tiers and RNA classes*). Runs with `--tiers` or `--passes` are marked as non-default at the top of the QC report; a resume with an edited order or passes file is refused.
- **`*-mappingstats.csv`** has two more columns, `Antisense` (unique reads on the opposite strand of their feature) and `AntisenseMulti` (multimappers whose assigned locus is antisense), counted whatever `remove_antisense` says. The reads removed as antisense are now counted from them, not derived as a difference, and a new consistency check, C11, compares the removed reads with the multimappers and antisense reads removed. New metrics `assign.antisense`, `assign.antisense_multi`.
- `alignments/logs/<sample>.cells.tsv`: which references ran in each pass (every run).

### Changed

- The aligners' logs carry the pass in their names (`alignments/logs/<sample>.p1.<index>.log`), also in default runs; so do the intermediate alignment files, which are deleted after the merge.
- The third column of `<species>.dbpaths.csv` (`retain`) must be `T`: every reference is kept and counted (`F`, which removed reads without counting them, is refused).
- QC report: the category "removed: antisense, maps once" is now "removed: antisense" (with `remove_multimap` on, antisense multimappers still count as multimappers).
- Own references: the Bowtie 1 commands `bowtie1`, `mind` and `ribospike` no longer pass `--no-unal`, which the Bowtie 1 of the image refuses, and read gzipped input; before, they could not run. No built-in reference uses them.
- QC report, own references with miND spike-ins: the read-length plot leaves the miND reads out of the natural reads, as the rest of the report does.

### Removed

- Settings: `collapse`, `subsample_reads`, `trim_fastq`, `trim_only`, `dump_unmapped`. With their old default value (`False`, `0`, `True`, `False`, `False`) they are ignored with a warning; any other value stops the run. To subsample, subsample the FASTQ files before the run.
- The file `<run>-ribospike.csv` at the top of the results folder, a copy of `ribospike/<run>-ribospike.csv`.
- seqtk and VSEARCH from the image (they served read collapsing and subsampling only).

### Other

- The software citation names Ruslan Klassen as co-author (README, CITATION.cff); the README has an *Acknowledgements* section.

## 0.3.0 – 2026-10-28

Results change (see *Changed results*), and many files, columns, RNA classes and settings have new names (see *Renamed*). Runs started with 0.2.x cannot be resumed with 0.3.0; start them again in an empty folder.

### Changed results

- **miRNA star arms** are counted as the class `miRNA` (before: their own class, `miRNA-star`). The feature names keep the arm (`_3p*`, `_5p*`), and the star arms are still mapped as a reference of their own, so no read is assigned differently.
- **snRNA and snoRNA** are split by position like tRNA: `_5p`, `_3p`, `_internal` and `_whole` rows in `*-readcounts.csv`, `*-featuredistinctseqs.csv` and `*-featureclasses.csv`, in the `CS` tag of their BAM records and in their `bedgraph/` tracks; the gene table sums them back to one row per gene. The 5′ end now follows the feature's strand (before: the reference coordinates, which swapped 5p and 3p for features on the minus strand; the built-in references have no tRNA, snRNA or snoRNA feature on the minus strand).
- **Vector contaminants** get their own class, `UniVec` (shown as "vector contaminants (UniVec)"); before, they were part of `unannotated`, which is now `genome_unannotated` and holds only reads on the genome without a feature.
- **Feature lengths of miRNA arms** (`length_nt` in `*-genelengths.csv` and `*-featureclasses.csv`) are filled in; before, they were empty. This includes arms that MirGeneDB merges into one feature (`A_5p|B_5p`) and, with `remove_antisense: False`, their `_antisense` rows.
- **`*-classlengths.csv`** has columns from 0 nt to the longest read of the run (at least 150 nt); before, they ended at 150 nt.
- `qc/metrics_long.tsv`: the new metric `assign.natural_mean_length` (mean read length of the natural reads), and the unit `features` instead of `genes` for the detection counts.
- **Gene symbols and miRBase names** in `*-featureclasses.csv` and `*-genelengths.csv`: `hgnc_id`, `hgnc_symbol` (human) or `mgi_id`, `mgi_symbol` (mouse), and `mirbase_name` (miRNA arms, star arms and precursors; miRBase 23.1). Built into `docker/tables/<species>.gene_symbols.tsv` from HGNC (2026-10-06), MGI, Ensembl 116, MirGeneDB 3.0 and miRBase, joined by Ensembl gene id (no second try by gene name), pre-mRNA through the gene names of the Ensembl GTF, miRNAs through MirGeneDB's miRBase accessions and tRNAs through GtRNAdb ids; several genes of one feature are separated by `;`. The columns stay empty, with a warning, when the references in use are not the ones the table was built for (BED checksum).

### New

- **The tier list and the RNA class table are files** (`docker/tables/` in the repository; `<run>/tables/` in every run). The class table gives each BED label its RNA class, the class's name, group and order in the report, its rank for overlapping features and whether it is in the gene table. Replace them with `--tiers FILE` and `--rna-classes FILE` (config keys `tiers`, `rna_classes`); unknown BED labels get a default row and a warning. In this version a label has one class and rank in every reference: a table whose rows of one label disagree is refused, and tiers that the tier list does not have are named in a warning. `run_status.json` records both tables (key `tables`), and the QC report shows them.
- **QC report:**
  - a *Samples* section with all columns of the sample sheet (on the config route: config key `sample_sheet`), with the result of the consistency checks on its first line;
  - *Colour samples by*: the reader chooses in the report which sample-sheet column colours the samples in every plot and table (at most 8 groups; a click on a colour in the key changes it); `report_group: COLUMN` sets the first choice;
  - a *How to read this report* section with a glossary;
  - read length by RNA class in two modes: % of natural reads (spike-ins, vector contaminants and unmapped reads left out) and % within class;
  - the mean read length of the natural reads in *General Statistics*;
  - one *Detection and counts* section with a plot of detected features per RNA class;
  - the tier list and the RNA class table of the run;
  - settings (Reads/Percentage, Linear/Log) shown the same way above every plot, a description under every plot, and the citations of the pipeline and the method.
- Section order of the report: *Samples*, *How to read this report*, *General Statistics*, then the steps of the pipeline. MultiQC's own General Statistics table is replaced by a section of the same name (the fastp and cutadapt numbers are in their own sections). The *Run overview* section is gone: the check summary is on the first line of *Samples*, the checks that are not `equal` and the warnings are under *General Statistics*, and run, mode and version are under *Provenance*. The spike-in length-offset plot shows percentages only.
- The gene-symbol table is listed under *Provenance* with whether it was applied (it is applied only to the references BED it was built for).
- `keep_classes` with an old class name (e.g. `miRNA-star`, `unannotated`) gives a warning and keeps the new class: `miRNA-star` now keeps every miRNA, since the star arms are part of miRNA.
- miND spike-ins (own references only) are in the report group *spike-ins* and are left out of the natural reads.

### Renamed

Scripts that read the result files need the new names. Of the old settings, only `keep_types` still works (with a warning). Own references with the old BED labels (`miRNA-star`, `miRNA-pre`, `RiboSpike-Short`) work unchanged: the RNA class table translates them.

| Kind | 0.2.x | 0.3.0 |
|---|---|---|
| File | `*-typecounts.csv` | `*-classcounts.csv` |
| File | `*-typelengths.csv` | `*-classlengths.csv` |
| File | `*-featuretypes.csv` | `*-featureclasses.csv` |
| File | `*-biotypedistinctseqs.csv` | `*-classdistinctseqs.csv` |
| File | `ribospike/*-ribotype.csv` | `ribospike/*-readorigin.csv` |
| File | `ribospike/*-ribolengths.csv` | `ribospike/*-spikelengths.csv` |
| Column | `classcounts`: first column (unnamed) | `rna_class` |
| Column | `classlengths`, `classdistinctseqs`, `genelengths`: `biotype` | `rna_class` |
| Column | `featureclasses`: `feature_type` | `rna_class` |
| Column | `featureclasses`: `gene_length`; `genelengths`: `length` | `length_nt` |
| Column | `readcounts`: first column (unnamed) | `feature` |
| Column | `mappingstats`, `trimstats`: first column (unnamed) | `sample` |
| Column | `*-ribospike.csv`: first column (unnamed); `ribospike/*-ribospike.csv`, `*-ribospikeends.csv`, `*-ribospikeends_exact.csv`: `ribospike` | `spike` |
| Column | `readorigin`: `biotype` | `read_origin` |
| Column | `spikelengths`: `biotype` | `spike_type` |
| RNA class | `miRNA-star` | `miRNA` |
| RNA class | `miRNA-pre` | `pre-miRNA` |
| RNA class | `unannotated` | `genome_unannotated`, and `UniVec` for vector contaminants |
| RNA class | `RiboSpike-Short` | `RiboSpike-Fragment` |
| Value | `Unmapped` (`readorigin`) | `unmapped` |
| Value | BAM tag `CS`: `RiboSpike-fragment` | `RiboSpike-Fragment` |
| Value | coverage files (`bedgraph/` column 5, `bedgraph_full/` column 7) | RNA class names as above |
| Metric | `assign.type.{type}.*` | `assign.class.{class}.*` |
| Metric | `counts.detected.{type}.gt{n}`, `counts.{type}.top20` | the class part renamed as above |
| Setting | `keep_types` | `keep_classes` |
| Report | "RNA type", "genes" (for rows of the gene table) | "RNA class", "features" |
| Report | "not aligned", "% unaligned" (reads aligned in no tier) | "unmapped", "% unmapped" |
| Report | groups "mRNA", "genome, precursors, introns" | "protein_coding", "genome: pre-mRNA, pre-tRNA" and "genome (unannotated)" |

The RNA classes `rRNA`, `Mt_rRNA`, `tRNA`, `Mt_tRNA`, `pre-tRNA`, `miRNA`, `piRNA`, `snRNA`, `snoRNA`, `scaRNA`, `misc_RNA`, `lncRNA`, `protein_coding`, `pre-mRNA`, `unmapped`, `RiboSpike-T1` to `T4` and `miND` keep their names. "RNA Type" (capital T) is kept for the end chemistry of the spike-ins, as in the method paper.

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
