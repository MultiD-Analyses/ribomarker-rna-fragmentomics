# Synthetic truth test

**Purpose.** Check that the pipeline's numbers are *correct*, not only consistent: a small set of synthetic reads in which every read has a known fate, run through the full pipeline (paired-end and single-end) and compared with the expected values. The expected values are derived from the read design only, never from a pipeline run. The same data are the pipeline's demo (`run --demo`).

## Making the test data

This folder **creates** the frozen test data (`bundle/`); it is rarely needed, because the bundle is part of the repository.

```bash
bash test/synthetic_qc/make_bundle.sh            # refuses if a bundle exists
bash test/synthetic_qc/make_bundle.sh --force    # replace it (new truth!)
```

Run from the repository root with the `ribomarker` conda environment active. It downloads chromosome 21 of Ensembl release 116 (~12 MB, the same release as the image's references) into `work/`, extracts the UniVec sequences from the reference index with `bowtie-inspect`, and writes `bundle/references`, `bundle/fastq`, `bundle/truth` and `bundle/MANIFEST.txt` (sources, pipeline commit, checksums).

## Files

| File | Role |
|---|---|
| `make_bundle.sh` | Runs the two generators below and writes `bundle/MANIFEST.txt` |
| `make_mini_assets.py` | Cuts the small reference set from `docker/assets`: a few records per tier, an 80-kb window of chromosome 21 around ETS2, matching BED rows, tier order (`tiers.csv`) |
| `make_reads.py` | Synthetic reads (fixed seed) and the truth files. The read categories are defined at the top of the file |
| `build_indexes.py` | Builds the Bowtie indexes of the small reference set (used by `docker/make_demo.sh`) |
| `bundle/` | The frozen test data; see `bundle/README.md` |
| `work/` | chr21 download and UniVec sequences. Not part of the repository |

## Read design

Read layout from the RiboMarker protocol (Appendix D, p. 20) and a check on raw reads: R1 = 1 extra base + insert + 3′ adapter …; R2 = reverse complement of the insert + R2 adapter …; 76 nt each, all bases Phred 40. The categories (`bundle/truth/categories.tsv`) cover adapter dimers, 10-nt inserts, a 150-nt insert whose mates cannot be merged, all 20 spike-ins, one fragment type per reference tier, a tRNA multimapper, an antisense mRNA fragment, intronic and intergenic genome fragments, random unmappable reads and 15-nt inserts, plus harder cases (C19–C31): near-multimappers, single mismatches, opposite-strand features, feature edges, spike-in offsets and fragments, and pairs that fastp removes. Where a fate also depends on aligner or fastp defaults, the category text in `make_reads.py` says so. Sample S3 holds only exploratory categories (C32–C35: 16- and 17-nt fragments): their outcome is reported, not checked.

Every insert that should map uniquely occurs **exactly once** (both strands, exact match) in the whole small reference set; the multimapper exactly twice; random reads not at all. Genome fragments avoid repeat-masked (lower-case) sequence and share no 12-mer with the transcript references. Trace any failure through `bundle/truth/reads_truth.tsv`.

## Not covered

Equal-rank overlapping features (order unspecified), spike-ins longer than designed, per-base depth (`bedgraph_full`), `collapse: True`, the NEB and SMART trimming branches, an empty sample.
