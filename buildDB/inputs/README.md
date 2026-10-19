# buildDB inputs — required local files

`buildDB_v1.py` downloads every reference it can. A few files **cannot** be
downloaded (they are synthetic or project-specific), so you must place them
here. The builder reads this folder via `--inputs` (default: `buildDB/inputs/`).

## Required (per build)

| File | What it is | Used for | Notes |
|------|-----------|----------|-------|
| `ribospike.fa` | The synthetic RiboSpike spike-in sequences (FASTA) | Builds the `ribospike` bowtie2 index (`ribospikebt2`) | Must be the **same** spike-in reference the wet-lab used. Not downloadable. You currently have this as `assets/ribospike.fa` in the repo — copy that file here (do **not** point the builder at `assets/`). |

That is the only required input for the default build.

## Optional — only if you enable the tRNA Sprinzl alignment

The tRNA Sprinzl alignment step (Infernal `cmalign` → `*-trnaalign.csv`) is
**OFF by default**, because that CSV is not part of the mapping `dbpaths` and the
pipeline does not consume it. If you turn it on, you also need:

| File | What it is | Notes |
|------|-----------|-------|
| `cm/trnamature-euk.cm` | Infernal covariance model (mature eukaryotic tRNA) | From the tRNAscan-SE / GtRNAdb model set. Requires `infernal` (`cmalign`) installed. |

(For non-eukaryotic clades you would instead provide `trnamature-arch.cm`,
`trnamature-bact.cm`, or `trnamature-mito.cm`.)

## What must NOT go here

Everything else (genome, GTF, tRNA tarball, miRNA, piRNA, rRNA, UniVec) is
downloaded by the script and recorded in `assets_v1/MANIFEST.json`. Do not stage
those manually — that would defeat the provenance tracking.
