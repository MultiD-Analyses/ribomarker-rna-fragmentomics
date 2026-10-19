#!/usr/bin/env python3
"""
builder_core.py -- reference-index building steps (frozen alongside buildDB_v1).

These functions turn already-downloaded reference files into the bowtie / bowtie2
indices + BED annotation that the mapping pipeline consumes. The logic is ported
from the original buildDB.py and kept close to it ON PURPOSE, so that assets_v?
stays FORMAT-compatible with the pipeline (identical N-padding, FASTA headers,
BED6 columns, index tools/names, and dbpaths layout). Only the surrounding
plumbing was cleaned up.

Design rules:
  * NO network here. Every function takes LOCAL input paths; the driver
    (buildDB_v1.py) is the only place that downloads, so provenance lives in one
    spot.
  * Each builder returns (dbrows, bedfiles):
        dbrows   : list of [map_cmd, rel_index_path, retain]
                   -> rows of {species}.dbpaths.csv (rel to the assets root)
        bedfiles : list of BED6 file paths to concatenate into {species}.bed
  * Index tools / names MATCH the production assets (verified 2026-07):
        univec    -> bowtie v1   name: univec_bt1        map_cmd 'bowtie'
        rrna      -> bowtie v1   name: {sp}-rrna_bt1      map_cmd 'bowtie'
        ribospike -> bowtie2     name: ribospike          map_cmd 'ribospikebt2'
        trna/mirna/pirna/cdna/ncrna/genome -> bowtie2     map_cmd 'bowtie2'

Requires on PATH: bowtie-build (v1), bowtie2-build, tar, gunzip.
Python deps: pandas, biopython (Bio.SeqIO), pysam.
"""

import os
import sys
import gzip
import re
import subprocess
import tarfile
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

NPAD = "N" * 10  # flanking pad added around every small-RNA reference sequence


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def run(cmd):
    """Run a shell command, raising on failure (unlike the old os.system)."""
    print(f"  $ {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True)


def _run_build(tool, fasta, index_prefix, logfile, extra=""):
    """Run an index builder; log verbose output but surface the tail on failure."""
    cmd = " ".join(x for x in [tool, extra, fasta, index_prefix] if x)
    print(f"  $ {cmd}", flush=True)
    with open(logfile, "a") as log:
        proc = subprocess.run(cmd, shell=True, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        tail = "".join(open(logfile).readlines()[-15:])
        raise RuntimeError(f"{tool} failed (exit {proc.returncode}); see {logfile}\n"
                           f"--- log tail ---\n{tail}")


def _bowtie1_build(fasta, index_prefix, threads, logfile):
    # bowtie v1's bowtie-build has NO threading option (unlike bowtie2-build).
    _run_build("bowtie-build", fasta, index_prefix, logfile)


def _bowtie2_build(fasta, index_prefix, threads, logfile):
    _run_build("bowtie2-build", fasta, index_prefix, logfile, extra=f"--threads {threads}")


def _only(names, wanted, what, source):
    """The one name that matches; an error if there is none or more than one."""
    found = [n for n in names if wanted(n)]
    if len(found) != 1:
        raise RuntimeError(f"expected exactly one {what} in {source}, found {found}")
    return found[0]


def _rel(out_root, path):
    """Path relative to the assets root, for the dbpaths file."""
    return os.path.relpath(path, out_root)


# --------------------------------------------------------------------------- #
# UniVec  (bowtie v1)
# --------------------------------------------------------------------------- #
def build_univec(univec_fa, out_root, threads):
    d = os.path.join(out_root, "univec")
    os.makedirs(d, exist_ok=True)
    index = os.path.join(d, "univec_bt1")
    _bowtie1_build(univec_fa, index, threads, os.path.join(d, "bowtie-build.log"))
    return [["bowtie", _rel(out_root, index), "T"]], []


# --------------------------------------------------------------------------- #
# RiboSpike  (bowtie2)  -- input is the synthetic spike-in FASTA (local)
# --------------------------------------------------------------------------- #
def build_ribospike(ribospike_fa, out_root, threads):
    d = os.path.join(out_root, "ribospike")
    os.makedirs(d, exist_ok=True)
    index = os.path.join(d, "ribospike")
    _bowtie2_build(ribospike_fa, index, threads, os.path.join(d, "bowtie-build.log"))
    return [["ribospikebt2", _rel(out_root, index), "T"]], []


# --------------------------------------------------------------------------- #
# rRNA  (mirtrace reference, filtered by species prefix)  (bowtie v1)
# --------------------------------------------------------------------------- #
def build_rrna(rrna_fa, species_prefix, out_root, threads, dbprefix):
    d = os.path.join(out_root, "rrna")
    os.makedirs(d, exist_ok=True)

    # keep only this species' rRNA (e.g. records starting 'Homo_sapiens_')
    species_fa = os.path.join(d, "rrna.species.fa")
    with open(species_fa, "w") as out:
        for rec in SeqIO.parse(rrna_fa, "fasta"):
            if rec.id.startswith(species_prefix):
                out.write(f">{rec.id}\n{rec.seq}\n")

    padded_fa = os.path.join(d, f"{dbprefix}.rrna.fa")
    bed_rows = []
    with open(padded_fa, "w") as f1:
        for rec in SeqIO.parse(species_fa, "fasta"):
            header = rec.id.strip().split()[0]
            seq = str(rec.seq).upper().replace("U", "T")
            f1.write(f">{header}\n{NPAD}{seq}{NPAD}\n")
            bed_rows.append([header, len(NPAD), len(seq) + len(NPAD), header, "rRNA", "+"])

    bed = os.path.join(d, f"{dbprefix}.rrna.bed")
    pd.DataFrame(bed_rows).drop_duplicates().to_csv(bed, sep="\t", header=False, index=False)

    index = os.path.join(d, f"{dbprefix}-rrna_bt1")
    _bowtie1_build(padded_fa, index, threads, os.path.join(d, "bowtie-build.log"))
    return [["bowtie", _rel(out_root, index), "T"]], [bed]


# --------------------------------------------------------------------------- #
# tRNA  (GtRNAdb tarball)  (bowtie2)
# --------------------------------------------------------------------------- #
def build_gtrnadb(trna_tar, out_root, threads, dbprefix, add_cca=True):
    d = os.path.join(out_root, "trna")
    os.makedirs(d, exist_ok=True)
    run(f"tar -xvf {trna_tar} -C {d}")

    # The FASTA and the BED are taken from THIS archive's file list, not from the folder: a folder that
    # already held another species' GtRNAdb files gave the mouse build the human tRNAs and pre-tRNA loci
    # (found 2026-10-07).
    with tarfile.open(trna_tar) as tar:
        names = [os.path.basename(m.name) for m in tar.getmembers() if m.isfile()]
    trna_fa = _only(names, lambda f: f.endswith("-tRNAs.fa") and not f.endswith("-mature-tRNAs.fa")
                    and not f.endswith("-filtered-tRNAs.fa"), "tRNA FASTA", trna_tar)
    trna_bed = _only(names, lambda f: f.endswith("-tRNAs.bed"), "tRNA BED", trna_tar)
    trna_fa = os.path.join(d, trna_fa)
    trna_bed = os.path.join(d, trna_bed)

    fastaD, fastabed, loci_list = {}, [], []
    for rec in SeqIO.parse(trna_fa, "fasta"):
        header = rec.id.strip().split()[0].split("_")[-1]
        species = rec.id.strip().split("_tRNA")[0]
        loci_list.append(header)
        parts = header.split("-")
        outheader = f"{species}_{parts[0]}-{parts[1]}-{parts[2]}-{parts[3]}"
        seq = str(rec.seq).upper().replace("U", "T")
        if add_cca and not header.startswith("tRNA-Lys-CTT-11"):
            seq = seq + "CCA"
        if outheader in fastaD:
            continue
        fastaD[outheader] = seq
        fastabed.append([outheader, len(NPAD), len(seq) + len(NPAD), outheader, "tRNA", "+"])

    # pre-tRNA loci (genomic), taken from the GtRNAdb BED
    loci = pd.read_csv(trna_bed, sep="\t", header=None)
    loci = loci[loci[3].isin(loci_list)].copy()
    loci["biotype"] = "pre-tRNA"
    loci = loci[[0, 1, 2, 3, "biotype", 5]]
    loci[0] = loci[0].str.replace("chr", "", regex=False)
    pretrnas_bed = os.path.join(d, "pretrnas.bed")
    loci.to_csv(pretrnas_bed, sep="\t", header=False, index=False)

    padded_fa = os.path.join(d, "trnas.fa")
    with open(padded_fa, "w") as f:
        for k, v in fastaD.items():
            f.write(f">{k}\n{NPAD}{v}{NPAD}\n")

    trnas_bed = os.path.join(d, "trnas.bed")
    pd.DataFrame(fastabed).drop_duplicates().to_csv(trnas_bed, sep="\t", header=False, index=False)

    index = os.path.join(d, f"{dbprefix}-trnas")
    _bowtie2_build(padded_fa, index, threads, os.path.join(d, "bowtie-build.log"))
    return [["bowtie2", _rel(out_root, index), "T"]], [trnas_bed, pretrnas_bed]


# --------------------------------------------------------------------------- #
# miRNA  (MirGeneDB)  -> mature / star / precursor  (bowtie2 x3)
# --------------------------------------------------------------------------- #
def build_mirna(mirna_fa, out_root, threads, dbprefix):
    d = os.path.join(out_root, "mirna")
    os.makedirs(d, exist_ok=True)

    mature, star, pre = defaultdict(list), defaultdict(list), defaultdict(list)
    for rec in SeqIO.parse(mirna_fa, "fasta"):
        header = rec.id.strip()
        seq = str(rec.seq).replace("U", "T")
        if header.endswith("_pre"):
            pre[seq].append(header)
        elif header.endswith("_3p*") or header.endswith("_5p*"):
            star[seq].append(header)
        elif header.endswith("3p") or header.endswith("5p"):
            mature[seq].append(header)

    def _write(group, biotype, fa_path):
        rows = []
        with open(fa_path, "w") as f:
            for seq, headers in group.items():
                h = "|".join(headers)
                f.write(f">{h}\n{NPAD}{seq}{NPAD}\n")
                rows.append([h, len(NPAD), len(seq) + len(NPAD), h, biotype, "+"])
        return pd.DataFrame(rows).drop_duplicates()

    mat_fa = os.path.join(d, f"{dbprefix}.mirna.fa")
    star_fa = os.path.join(d, f"{dbprefix}.star.fa")
    pre_fa = os.path.join(d, f"{dbprefix}.pre.fa")
    mat_bed = _write(mature, "miRNA", mat_fa)
    star_bed = _write(star, "miRNA-star", star_fa)
    pre_bed = _write(pre, "miRNA-pre", pre_fa)

    bed = os.path.join(d, f"{dbprefix}.mirna.bed")
    pd.concat([mat_bed, star_bed, pre_bed], axis=0).to_csv(bed, sep="\t", header=False, index=False)

    log = os.path.join(d, "bowtie-build.log")
    idx_mat = os.path.join(d, f"{dbprefix}-mirna")
    idx_star = os.path.join(d, f"{dbprefix}-star")
    idx_pre = os.path.join(d, f"{dbprefix}-pre")
    _bowtie2_build(mat_fa, idx_mat, threads, log)
    _bowtie2_build(star_fa, idx_star, threads, log)
    _bowtie2_build(pre_fa, idx_pre, threads, log)

    dbrows = [
        ["bowtie2", _rel(out_root, idx_mat), "T"],
        ["bowtie2", _rel(out_root, idx_star), "T"],
        ["bowtie2", _rel(out_root, idx_pre), "T"],
    ]
    return dbrows, [bed]


# --------------------------------------------------------------------------- #
# piRNA  (piRBase gold set)  (bowtie2)
# --------------------------------------------------------------------------- #
def build_pirbase(pirna_fa, out_root, threads, dbprefix):
    d = os.path.join(out_root, "pirna")
    os.makedirs(d, exist_ok=True)

    fastaD = defaultdict(list)
    for rec in SeqIO.parse(pirna_fa, "fasta"):
        header = rec.description.strip().split()[-1]
        if "tRNA" in header:
            continue
        fastaD[str(rec.seq)].append(header)

    padded_fa = os.path.join(d, f"{dbprefix}-pirna.fa")
    bed_rows = []
    with open(padded_fa, "w") as f:
        for seq, ids in fastaD.items():
            h = ids[0] if len(ids) == 1 else "|".join(ids)
            f.write(f">{h}\n{NPAD}{seq}{NPAD}\n")
            bed_rows.append([h, len(NPAD), len(seq) + len(NPAD), h, "piRNA", "+"])

    bed = os.path.join(d, f"{dbprefix}-pirna.bed")
    pd.DataFrame(bed_rows).drop_duplicates().to_csv(bed, sep="\t", header=False, index=False)

    index = os.path.join(d, f"{dbprefix}-pirna")
    _bowtie2_build(padded_fa, index, threads, os.path.join(d, "bowtie-build.log"))
    return [["bowtie2", _rel(out_root, index), "T"]], [bed]


# --------------------------------------------------------------------------- #
# transcriptome: GTF -> longest-per-gene BED12 -> cDNA + ncRNA fasta (bowtie2 x2)
# --------------------------------------------------------------------------- #
def _extract_fasta_from_bed12(bed12_file, genome_fa, out_fa, kind):
    """Reconstruct spliced transcript sequences from BED12 blocks + genome."""
    import pysam
    fasta = pysam.FastaFile(genome_fa)
    out_prefix = out_fa[:-3] if out_fa.endswith(".fa") else out_fa

    tx_bed, genome_bed = [], []
    with open(out_fa, "w") as out, open(bed12_file) as bed:
        for line in bed:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            chrom, start, end, name, biotype, strand, _s2, _e2, gene_name, \
                block_count, block_sizes, block_starts = fields[:12]
            start, end, block_count = int(start), int(end), int(block_count)
            sizes = [int(x) for x in block_sizes.split(",") if x != ""]
            starts = [int(x) for x in block_starts.split(",") if x != ""]

            seq = ""
            for i in range(block_count):
                bs = start + starts[i]
                seq += fasta.fetch(chrom, bs, bs + sizes[i])
            seq = seq.upper().replace("U", "T")
            if strand == "-":
                seq = seq.translate(str.maketrans("NATCGatcg", "NTAGCtagc"))[::-1]

            outname = f"{name}|{gene_name}"
            out.write(f">{outname}\n{NPAD}{seq}{NPAD}\n")
            tx_bed.append([outname, len(NPAD), len(seq) + len(NPAD), gene_name, biotype, "+"])
            genome_bed.append([chrom, start, end, gene_name + "_intron", biotype, strand])

    tx_bed_path = f"{out_prefix}.bed"
    pd.DataFrame(tx_bed).to_csv(tx_bed_path, sep="\t", index=False, header=False)

    gdf = pd.DataFrame(genome_bed, columns=["chrom", "start", "end", "name", "biotype", "strand"])
    genome_bed_path = f"{out_prefix}.genome.bed"
    if "cdna" in out_fa:
        gdf["biotype"] = "pre-mRNA"
    gdf.to_csv(genome_bed_path, sep="\t", index=False, header=False)
    return tx_bed_path, genome_bed_path


def build_transcriptome(gtf_path, genome_fa, out_root, threads, dbprefix):
    d = os.path.join(out_root, "gtf")
    os.makedirs(d, exist_ok=True)
    out_file = os.path.join(d, dbprefix)

    open_func = gzip.open if gtf_path.endswith(".gz") else open
    transcripts = defaultdict(lambda: {"exons": [], "strand": None, "chrom": None,
                                       "gene_id": None, "gene_name": None, "biotype": None})
    with open_func(gtf_path, "rt") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] not in ("transcript", "exon"):
                continue
            chrom, start, end, strand, attrs = fields[0], int(fields[3]) - 1, int(fields[4]), fields[6], fields[8]
            try:
                tid = re.search('transcript_id "(.*?)";', attrs).group(1)
                gid = re.search('gene_id "(.*?)";', attrs).group(1)
                gname = re.search('gene_name "(.*?)";', attrs).group(1) if "gene_name" in attrs else gid
                biotype = re.search('gene_biotype "(.*?)";', attrs).group(1)
            except AttributeError:
                continue
            if fields[2] == "transcript":
                transcripts[tid].update(chrom=chrom, strand=strand, gene_id=gid,
                                        gene_name=gname, biotype=biotype)
            else:
                transcripts[tid]["exons"].append((start, end))

    tmp_bed12 = f"{out_file}.tmp.bed12"
    with open(tmp_bed12, "w") as bed:
        for tid, info in transcripts.items():
            if not info["exons"]:
                continue
            exons = sorted(info["exons"])
            tx_start, tx_end = exons[0][0], exons[-1][1]
            sizes = [e - s for s, e in exons]
            starts = [s - tx_start for s, _ in exons]
            bed.write(f'{info["chrom"]}\t{tx_start}\t{tx_end}\t{tid}\t{info["biotype"]}\t{info["strand"]}\t'
                      f'{tx_start}\t{tx_end}\t{info["gene_name"]}\t{len(exons)}\t'
                      f'{",".join(map(str, sizes))},\t{",".join(map(str, starts))},\n')

    bed12 = pd.read_csv(tmp_bed12, sep="\t", header=None).sort_values(by=[0, 1, 2])
    unq = bed12.copy()
    unq["blen"] = unq[10].apply(lambda x: sum(int(i) for i in x.split(",") if i != ""))
    unq = unq.sort_values(by=[0, 1, 2, "blen"], ascending=[True, True, True, False])
    unq = unq.drop_duplicates(subset=[4, 8], keep="first").drop("blen", axis=1)  # longest per (biotype, gene)
    bed12.to_csv(f"{out_file}.bed12", sep="\t", index=False, header=False)
    unq.to_csv(f"{out_file}.longest.bed12", sep="\t", index=False, header=False)
    os.remove(tmp_bed12)

    unq = unq[~unq[4].isin(["miRNA", "tRNA"])]

    # protein_coding -> cDNA
    pc = unq[unq[4] == "protein_coding"]
    pc_bed12 = f"{out_file}.protein_coding.bed12"
    pc.to_csv(pc_bed12, sep="\t", index=False, header=False)
    cdna_bed, cdna_genome_bed = _extract_fasta_from_bed12(pc_bed12, genome_fa, f"{out_file}.cdna.fa", "cdna")
    os.remove(pc_bed12)

    # non-coding classes -> ncRNA
    nc_types = ["lncRNA", "snoRNA", "snRNA", "rRNA", "scaRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA"]
    nc = unq[unq[4].isin(nc_types)]
    nc_bed12 = f"{out_file}.ncRNA.bed12"
    nc.to_csv(nc_bed12, sep="\t", index=False, header=False)
    ncrna_bed, _ = _extract_fasta_from_bed12(nc_bed12, genome_fa, f"{out_file}.ncRNA.fa", "ncrna")

    log = f"{out_file}.bowtie-build.log"
    idx_cdna = os.path.join(d, f"{dbprefix}-cdna")
    idx_ncrna = os.path.join(d, f"{dbprefix}-ncrna")
    _bowtie2_build(f"{out_file}.cdna.fa", idx_cdna, threads, log)
    _bowtie2_build(f"{out_file}.ncRNA.fa", idx_ncrna, threads, log)

    # dbpaths ORDER (matches production): ncrna before cdna
    dbrows = [
        ["bowtie2", _rel(out_root, idx_ncrna), "T"],
        ["bowtie2", _rel(out_root, idx_cdna), "T"],
    ]
    return dbrows, [cdna_bed, ncrna_bed, cdna_genome_bed]


# --------------------------------------------------------------------------- #
# whole genome  (bowtie2)
# --------------------------------------------------------------------------- #
def build_genome(genome_fa, out_root, threads, dbprefix):
    d = os.path.join(out_root, "genome")
    os.makedirs(d, exist_ok=True)
    index = os.path.join(d, f"{dbprefix}-genome")
    _bowtie2_build(genome_fa, index, threads, os.path.join(d, "bowtie-build.log"))
    return [["bowtie2", _rel(out_root, index), "T"]], []


# --------------------------------------------------------------------------- #
# assemble the final dbpaths + concatenated BED
# --------------------------------------------------------------------------- #
def write_outputs(dbrows, bedfiles, out_root, species):
    """
    Write {species}.dbpaths.csv and the concatenated {species}.bed. dbrows carry
    paths RELATIVE to out_root, and the dbpaths file keeps them relative, so the
    folder works wherever it is (locally, mounted in Docker). Until
    2026-10-08 this wrote {species}.local.dbpaths.csv (absolute paths) and
    {species}.docker.dbpaths.csv (/pipeline/assets/ paths) instead.
    """
    dbpaths_csv = os.path.join(out_root, f"{species}.dbpaths.csv")
    with open(dbpaths_csv, "w") as f:
        f.write("\n".join(f"{map_cmd},{rel},{retain}" for map_cmd, rel, retain in dbrows) + "\n")

    bedfile = os.path.join(out_root, f"{species}.bed")
    with open(bedfile, "w") as out:
        for bed in bedfiles:
            with open(bed) as b:
                out.write(b.read())

    return dbpaths_csv, bedfile
