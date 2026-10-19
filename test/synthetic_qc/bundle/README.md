# Synthetic RiboMarker test bundle

Small, frozen test data with known truth, for testing any small-RNA / fragmentomics pipeline on RiboMarker-type libraries. Created by `../make_bundle.sh`; provenance and checksums in `MANIFEST.txt`. Do not edit by hand: the level-A test (`test/levels/run_level_a.sh`) refuses to run if a file differs from the manifest.

## Contents

| Path | Content |
|---|---|
| `fastq/` | Three samples (S1, S2; S3 holds only the exploratory categories C32–C35), paired-end, 76 nt, gzip. `samples_pe.csv` / `samples_se.csv` are sample sheets for the MultiD pipeline; single-end mode uses the R1 files only |
| `references/` | The reference set the reads were designed against: one FASTA per mapping tier (`tiers.csv` gives the search order), `human.bed` with the annotation. Small-RNA and transcript sequences are padded with 10 N on both ends, as in the production assets; the genome tier is chromosome 21, positions 38,775,000–38,855,000 (0-based, Ensembl 116), named `21`, with BED coordinates relative to the window; two genome BED rows are synthetic (`SYNTH-SNORD-minus`, `SYNTH-LNC-edge`, see `../make_mini_assets.py`) |
| `truth/reads_truth.tsv` | **Pipeline-independent truth, one row per read pair** (columns below) |
| `truth/truth.tsv` | Expected values in the MultiD pipeline's own output tables, per mode and sample; table `qc` holds expected QC metric values (ids as in `docker/qc/metrics.yaml`), `tier_aligned_multi` reads with two or more reported alignments per tier, `exploratory` the S3 categories (reads per category, outcome not checked) |
| `truth/categories.tsv` | The 35 read categories and their counts per sample |

## Read layout

R1 = 1 extra base + insert + `TGGAATTCTCGGGTGCCAAGG…` (3′ adapter), R2 = reverse complement of the insert + `TGATCGTCGGACT…`, both 76 nt, filled with poly-G; all qualities Phred 40 except category C30 (Phred 10). So a correct pipeline removes the first base of R1 only, trims both adapters, and (paired-end) merges the mates back into the insert.

## `reads_truth.tsv` columns

| Column | Meaning |
|---|---|
| `read`, `sample`, `category` | Read name (as in the FASTQ), sample, category (C01–C35) |
| `biotype`, `feature` | What the insert is (e.g. `protein_coding`, `ETS2`); multimapper features are joined with `|` |
| `orientation` | `sense` / `antisense` relative to the feature; empty if not applicable |
| `insert_length`, `insert_sequence` | The true insert, 5′→3′ in RNA orientation |
| `source_contig`, `source_start` | Where the insert was taken from: contig in `references/` and 0-based start (padded contigs include the 10-N pad; genome positions are relative to the window) |
| `fate_paired_end`, `fate_single_end` | `lost:dimer` (empty insert), `lost:short` (< 15 nt), `lost:merge` (paired-end mates do not overlap), `lost:fastp` (paired-end pair removed by fastp's filters), `mapped`, `unmapped`, `exploratory` (outcome not defined) |
| `reference_tier`, `alignments` | First tier in `tiers.csv` order where the read aligns, and the number of reported alignments there (1; 2 for the multimapper C08, and for C19, whose second alignment has one mismatch; 0 if unmapped) |

Pipeline-specific choices (e.g. that multimappers and antisense reads are dropped from the MultiD counts) are in `truth.tsv`, not here.

## Sources

Sequences are taken from public databases via the MultiD `buildDB` assets: Ensembl 116 (GRCh38 genome, cDNA, ncRNA), GtRNAdb (tRNA), MirGeneDB 3.0 (miRNA), piRBase v3.0 (piRNA), mirtrace rRNA reference, NCBI UniVec, and the RealSeq RiboSpike design (`buildDB/inputs/ribospike.fa`). Sources, citations and redistribution terms: `THIRD_PARTY_NOTICES.md` at the repository root.
