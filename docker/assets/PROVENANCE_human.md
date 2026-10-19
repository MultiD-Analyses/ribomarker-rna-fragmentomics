# Provenance -- assets v1 (human)

Built: 2026-07-07T11:49:06Z  |  builder: buildDB_v1.py (sha256 `e51e855bb26b...`)

## Sources

| DB | Version | URL | Downloaded (UTC) | sha256 (first 12) |
|----|---------|-----|------------------|-------------------|
| genome | Ensembl release-116 / GRCh38 | https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz | 2026-07-07T11:51:39Z | `746dc3b40529` |
| gtf | Ensembl release-116 / GRCh38 | https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz | 2026-07-07T11:55:37Z | `ed992f0eac71` |
| trna | GtRNAdb2 hg38 | https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz | 2026-07-07T11:55:39Z | `55a354318748` |
| mirna | MirGeneDB 3.0 | https://mirgenedb.org/fasta/hsa?all=1 | 2026-07-07T11:55:39Z | `0cebd4662f6c` |
| pirna | piRBase v3.0 (gold) | http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/hsa.gold.fa.gz | 2026-07-07T11:55:46Z | `6a62c4fa96c3` |
| rrna | mirtrace rRNA_reference.fa (branch: master) | https://raw.githubusercontent.com/friedlanderlab/mirtrace/master/src/lib/inputs/rRNA_reference.fa | 2026-07-07T11:55:47Z | `4cb69b85228b` |
| univec | NCBI UniVec (current build) | https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec | 2026-07-07T11:55:50Z | `2fb50fd9749e` |

## Local inputs

- ribospike: `buildDB/inputs/ribospike.fa` (sha256 `25143e9dde9f`)

## Tools

- python: 3.13.13
- bowtie-build: bowtie-build version 1.0.0
- bowtie2-build: bowtie2-build-s version 2.5.5
