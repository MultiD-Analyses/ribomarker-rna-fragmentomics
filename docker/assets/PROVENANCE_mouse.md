# Provenance -- assets v1 (mouse)

Built: 2026-10-08T15:48:57Z  |  builder: buildDB_v1.py (sha256 `1fc0dee301f6...`)

## Sources

| DB | Version | URL | Downloaded (UTC) | sha256 (first 12) |
|----|---------|-----|------------------|-------------------|
| genome | Ensembl release-116 / GRCm39 | https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz | 2026-10-08T15:51:58Z | `0669729d3276` |
| gtf | Ensembl release-116 / GRCm39 | https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz | 2026-10-08T15:53:45Z | `5c29fd9e3157` |
| trna | GtRNAdb2 mm39 | https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz | 2026-10-08T15:53:48Z | `fbcc8b2dd8a1` |
| mirna | MirGeneDB 3.0 | https://mirgenedb.org/fasta/mmu?all=1 | 2026-10-08T15:53:49Z | `57855e35a404` |
| pirna | piRBase v3.0 (gold) | http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/mmu.gold.fa.gz | 2026-10-08T15:54:03Z | `459f2228fd62` |
| rrna | mirtrace rRNA_reference.fa (branch: master) | https://raw.githubusercontent.com/friedlanderlab/mirtrace/master/src/lib/inputs/rRNA_reference.fa | 2026-10-08T15:54:04Z | `4cb69b85228b` |
| univec | NCBI UniVec (current build) | https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec | 2026-10-08T15:54:06Z | `2fb50fd9749e` |

## Local inputs

- ribospike: `buildDB/inputs/ribospike.fa` (sha256 `25143e9dde9f`)

## Tools

- python: 3.13.13
- bowtie-build: bowtie-build version 1.0.0
- bowtie2-build: bowtie2-build-s version 2.5.5

## ncRNA split (release 0.4.0)

2026-11-03T11:47:10Z: `split_ncrna.py` (sha256 `4d426a3b234a`), bowtie2-build-s version 2.5.5. The Ensembl ncRNA reference `gtf/mouse.ncRNA.fa` (sha256 `d9916f20a916`) cut into four, every record unchanged and in its order:

| Sub-reference | Records | Biotypes |
|---|---|---|
| `gtf/mouse-ncrna_mt` | 24 | Mt_tRNA 22, Mt_rRNA 2 |
| `gtf/mouse-ncrna_snsno` | 2911 | snRNA 1381, snoRNA 1484, scaRNA 46 |
| `gtf/mouse-ncrna_other` | 916 | misc_RNA 562, rRNA 354 |
| `gtf/mouse-ncrna_lnc` | 32879 | lncRNA 32879 |
