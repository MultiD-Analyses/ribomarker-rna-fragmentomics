# Third-party notices

The pipeline's own code is under the MIT License ([LICENSE](LICENSE)). The container image also contains third-party software and reference data, each under its own licence or terms, listed here. This file is informational; the licence texts of the software are in the image (see *Software*).

## Reference data

The human and mouse references in the image are built from the sources below (exact files, versions, download dates and checksums: `docker/assets/MANIFEST_<species>.json` and `docker/assets/PROVENANCE_<species>.md`). If you use results based on them, please cite the sources.

| Used for | Source | Terms | Please cite |
|---|---|---|---|
| Genome, ncRNA, cDNA, annotation | [Ensembl](https://www.ensembl.org), release 116 (GRCh38, GRCm39) | "Ensembl imposes no restrictions on access to, or use of, the data"; some included data may be subject to third-party constraints ([Ensembl disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html)) | Dyer SC et al. Ensembl 2025. *Nucleic Acids Res* 53(D1):D948–D957 (2025). https://doi.org/10.1093/nar/gkae1071 — and the release used (116) |
| Vector contamination | [NCBI UniVec](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) | See *NCBI* below | NCBI UniVec (download date in `PROVENANCE_<species>.md`) |
| rRNA | `rRNA_reference.fa` of [miRTrace](https://github.com/friedlanderlab/mirtrace), sequences curated from NCBI Nucleotide, SILVA and Ensembl | miRTrace is distributed under the GNU GPL v3; we redistribute this file under the same licence. SILVA data are under [CC BY 4.0](https://www.arb-silva.de/silva-license-information/) | Kang W et al. miRTrace reveals the organismal origins of microRNA sequencing data. *Genome Biol* 19:213 (2018). https://doi.org/10.1186/s13059-018-1588-9 ; Quast C et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Res* 41(D1):D590–D596 (2013). https://doi.org/10.1093/nar/gks1219 |
| tRNA | [GtRNAdb](https://gtrnadb.ucsc.edu) | Redistributed with the permission of the GtRNAdb authors | Chan PP, Lowe TM. GtRNAdb 2.0: an expanded database of transfer RNA genes identified in complete and draft genomes. *Nucleic Acids Res* 44(D1):D184–D189 (2016). https://doi.org/10.1093/nar/gkv1309 |
| miRNA | [MirGeneDB 3.0](https://mirgenedb.org) | Redistributed with the permission of the MirGeneDB authors | Clarke AW et al. MirGeneDB 3.0: improved taxonomic sampling, uniform nomenclature of novel conserved microRNA families and updated covariance models. *Nucleic Acids Res* 53(D1):D116–D128 (2025). https://doi.org/10.1093/nar/gkae1094 |
| piRNA | [piRBase](http://bigdata.ibp.ac.cn/piRBase/) v3.0 | Redistributed with the permission of the piRBase authors | Wang J et al. piRBase: integrating piRNA annotation in all aspects. *Nucleic Acids Res* 50(D1):D265–D272 (2022). https://doi.org/10.1093/nar/gkab1012 |
| Spike-in controls | RiboSpike sequences, RealSeq Biosciences | Published with the permission of RealSeq Biosciences | Clark RC et al. *RNA* 32(7):1156–1181 (2026). https://doi.org/10.1261/rna.080861.125 |
| Gene symbols, human (`docker/tables/human.gene_symbols.tsv`; release 0.3.0) | [HGNC](https://www.genenames.org) complete set, quarterly release (version and checksum in `docker/tables/human.gene_symbols.provenance.md`); joined through Ensembl, MirGeneDB and miRBase identifiers | "All data is released under the Creative Commons Public Domain (CC0) License" ([HGNC licence](https://www.genenames.org/about/license)) | Seal RL et al. Genenames.org: the HGNC resources in 2023. *Nucleic Acids Res* 51(D1):D1003–D1009 (2023). https://doi.org/10.1093/nar/gkac888 ; HGNC (RRID:SCR_002827) |
| Gene symbols, mouse (`docker/tables/mouse.gene_symbols.tsv`; release 0.3.0) | [MGI](https://www.informatics.jax.org), report `MRK_ENSEMBL.rpt` (download date and checksum in `docker/tables/mouse.gene_symbols.provenance.md`) | MGI data and annotations are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0), provided "as is" without warranty ([MGI copyright and licensing](https://www.informatics.jax.org/mgihome/other/copyright.shtml), checked 2026-10-27). Attribution: this row and the citation | Baldarelli RM, Smith CL, Ringwald M, Richardson JE, Bult CJ, Mouse Genome Informatics Group. Mouse Genome Informatics: an integrated knowledgebase system for the laboratory mouse. *Genetics* 227(1):iyae031 (2024). https://doi.org/10.1093/genetics/iyae031 ([MGI citing page](https://www.informatics.jax.org/mgihome/other/citation.shtml); data retrieved: see the provenance file) |
| miRNA names in the gene-symbol tables (`mirbase_name`, release 0.3.0; only accession → name pairs are used, no sequences are redistributed) | [miRBase](https://www.mirbase.org) release 23.1, `mature.fa`, `hairpin.fa` (download date and checksum in `docker/tables/<species>.gene_symbols.provenance.md`) | "miRBase is in the public domain. It is not copyrighted. You may freely modify, redistribute, or use it for any purpose" ([miRBase LICENSE](https://mirbase.org/download/LICENSE), release 23.1, checked 2026-10-27); miRBase asks to be cited, and the primary sources of the miRNA sequences too | Kozomara A, Birgaoanu M, Griffiths-Jones S. miRBase: from microRNA sequences to function. *Nucleic Acids Res* 47:D155–D162 (2019). https://doi.org/10.1093/nar/gky1141; and, as miRBase asks: Kozomara A, Griffiths-Jones S. *Nucleic Acids Res* 42:D68–D73 (2014); Kozomara A, Griffiths-Jones S. *Nucleic Acids Res* 39:D152–D157 (2011); Griffiths-Jones S, Saini HK, van Dongen S, Enright AJ. *Nucleic Acids Res* 36:D154–D158 (2008); Griffiths-Jones S, Grocock RJ, van Dongen S, Bateman A, Enright AJ. *Nucleic Acids Res* 34:D140–D144 (2006); Griffiths-Jones S. *Nucleic Acids Res* 32:D109–D111 (2004) |

### NCBI

From NCBI's [policies](https://www.ncbi.nlm.nih.gov/home/about/policies/): "NCBI itself places no restrictions on the use or distribution of the data contained therein." Some submitters may claim rights to their data, and NCBI cannot grant permission to use such data. "NCBI does not endorse or recommend any commercial products, processes, or services." The U.S. Government does not warrant or assume legal liability for the accuracy, completeness or usefulness of the data; see the policies page for the full disclaimer.

## Software

The image is based on [Miniforge](https://github.com/conda-forge/miniforge) and contains the conda packages listed, with exact builds, in [`docker/env.lock.txt`](docker/env.lock.txt) (from the conda-forge and bioconda channels). Each package's licence is recorded in the image under `/opt/conda/envs/ribomarker/conda-meta/` (field `license` of each package's JSON file), and the licence texts that each package provides are in `/opt/conda/envs/ribomarker/share/conda-licenses/<package>/`. A few packages do not include their licence text (bowtie, bedtools, dnaio, pysam and the GCC runtime libraries); their texts, taken from the upstream source at the packaged version, are in `conda-licenses/_added/` (in this repository: [`docker/licenses/`](docker/licenses/)). To read them: `docker run --rm --entrypoint sh IMAGE -c 'ls /opt/conda/envs/ribomarker/share/conda-licenses'`.

The programs and libraries the pipeline calls directly:

| Package | Version | Licence (as recorded by the conda package) | Project |
|---|---|---|---|
| Python | 3.13.13 | Python-2.0 | https://www.python.org |
| Snakemake (snakemake-minimal) | 9.20.0 | MIT | https://snakemake.github.io |
| cutadapt | 5.2 | MIT | https://github.com/marcelm/cutadapt |
| fastp | 1.3.3 | MIT | https://github.com/OpenGene/fastp |
| Bowtie | 1.0.0 | Artistic License (upstream `LICENSE` at v1.0.0; the conda package records GPLv3) | https://github.com/BenLangmead/bowtie |
| Bowtie 2 | 2.5.5 | GPL-3.0-or-later | https://github.com/BenLangmead/bowtie2 |
| SAMtools | 1.23.1 | MIT | https://github.com/samtools/samtools |
| seqtk | 1.5 | MIT | https://github.com/lh3/seqtk |
| VSEARCH | 2.31.0 | GPL-3.0-or-later OR BSD-2-Clause | https://github.com/torognes/vsearch |
| pigz | 2.8 | Zlib | https://zlib.net/pigz/ |
| MultiQC | 1.35 | GPL-3.0-or-later | https://github.com/MultiQC/MultiQC |
| pandas | 3.0.2 | BSD-3-Clause | https://pandas.pydata.org |
| NumPy | 2.4.3 | BSD-3-Clause | https://numpy.org |
| pysam | 0.24.0 | MIT | https://github.com/pysam-developers/pysam |
| PyYAML | 6.0.3 | MIT | https://pyyaml.org |
| tqdm | 4.67.3 | MPL-2.0 and MIT | https://github.com/tqdm/tqdm |

### Corresponding source

Some packages in the image are under licences that require their source code to be available:

- **GPL:** bowtie2 2.5.5 (GPL-3.0-or-later), multiqc 1.35 (GPL-3.0-or-later), vsearch 2.31.0 (GPL-3.0-or-later OR BSD-2-Clause), readline 8.3 (GPL-3.0-only), lzo 2.10 (GPL-2.0-or-later), ld_impl_linux-64 2.45.1 (GPL-3.0-only), perl 5.32.1 (GPL-1.0-or-later OR Artistic-1.0-Perl), freetype / libfreetype / libfreetype6 2.14.3 (GPL-2.0-only OR FTL), dbus 1.16.2 (AFL-2.1 OR GPL-2.0-or-later); the GCC runtime libraries libgcc, libgcc-ng, libgfortran, libgfortran5, libgomp, libstdcxx, libstdcxx-ng 15.2.0 (GPL-3.0-only WITH GCC-exception-3.1).
- **LGPL:** alsa-lib 1.2.15.3, cairo 1.18.4 (LGPL-2.1-only or MPL-1.1), graphite2 1.3.14, keyutils 1.6.3, libglib 2.88.1, libiconv 1.18, libntlm 1.8, libxcrypt 4.4.36, pyside6 6.11.0, qt6-main 6.11.0.

The pipeline runs the GPL programs it uses (Bowtie, Bowtie 2, MultiQC, VSEARCH) as separate programs; it does not include their code in its own. Their source code is available from the projects listed above and from the package recipes of [conda-forge](https://github.com/conda-forge) and [bioconda](https://github.com/bioconda/bioconda-recipes), which record the exact source archive and checksum of each build. In addition, for at least three years after each release, MultiD Analyses AB will provide the complete corresponding source code of any GPL- or LGPL-licensed package in that release's image, on request and free of charge: open an issue in this repository.
