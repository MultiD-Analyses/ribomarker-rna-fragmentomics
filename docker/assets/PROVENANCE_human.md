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

## ncRNA split (release 0.4.0)

2026-11-03T11:44:40Z: `split_ncrna.py` (sha256 `4d426a3b234a`), bowtie2-build-s version 2.5.5. The Ensembl ncRNA reference `gtf/human.ncRNA.fa` (sha256 `dfcd632595d8`) cut into four, every record unchanged and in its order:

| Sub-reference | Records | Biotypes |
|---|---|---|
| `gtf/human-ncrna_mt` | 24 | Mt_tRNA 22, Mt_rRNA 2 |
| `gtf/human-ncrna_snsno` | 2720 | snRNA 1837, snoRNA 832, scaRNA 51 |
| `gtf/human-ncrna_other` | 1324 | misc_RNA 1282, rRNA 42 |
| `gtf/human-ncrna_lnc` | 34997 | lncRNA 34997 |
