#!/usr/bin/env python3
"""
make_mini_assets.py -- cut a SMALL reference set for the synthetic test bundle.

What it does
  1. Takes a handful of records from the production reference FASTA files in
     docker/assets. These files are the buildDB outputs, so names, N-padding
     and sequences are identical to what the pipeline normally uses.
  2. Cuts an 80-kb window of chromosome 21 around ETS2 as the genome tier.
  3. Keeps the matching rows of docker/assets/human.bed; genome rows are
     shifted into window coordinates. Adds two synthetic genome features
     (SYNTHETIC_BED_ROWS) for the hard cases.
  4. Writes one FASTA per tier, human.bed and tiers.csv (mapping order) into --out.
No indexes and no machine-specific paths: build_indexes.py makes those at run time.
Nothing in buildDB/ or docker/assets/ is modified.

Usage (from the repository root):
  python test/synthetic_qc/make_mini_assets.py \
      --assets docker/assets --chr21 <chr21.fa[.gz]> --univec <univec.fa> \
      --out test/synthetic_qc/bundle/references
"""
import argparse
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

# Genome tier: window of chromosome 21 (0-based, half-open) around ETS2.
# ETS2 gene body in human.bed: 21:38805114-38824956 (+); no other annotation within 30 kb.
GENOME_CONTIG = "21"
WINDOW_START, WINDOW_END = 38_775_000, 38_855_000

# One row per tier, in the production mapping order.
# (map_cmd, folder, index name, source FASTA relative to --assets, records to keep)
TIERS = [
    ("bowtie",       "univec",    "univec_bt1",     "UNIVEC",               "first 3"),
    ("ribospikebt2", "ribospike", "ribospike",      "RIBOSPIKE",            "all"),
    ("bowtie",       "rrna",      "human-rrna_bt1", "rrna/human.rrna.fa",   ["Homo_sapiens_18S_1"]),
    ("bowtie2",      "trna",      "human-trnas",    "trna/trnas.fa",        ["Homo_sapiens_tRNA-Ala-CGC-3",
                                                                             "Homo_sapiens_tRNA-Leu-CAG-1",
                                                                             "Homo_sapiens_tRNA-Leu-CAG-2"]),
    ("bowtie2",      "mirna",     "human-mirna",    "mirna/human.mirna.fa", ["Hsa-Mir-21_5p"]),
    ("bowtie2",      "mirna",     "human-star",     "mirna/human.star.fa",  ["Hsa-Mir-21_3p*"]),
    ("bowtie2",      "mirna",     "human-pre",      "mirna/human.pre.fa",   ["Hsa-Mir-21_pre"]),
    ("bowtie2",      "pirna",     "human-pirna",    "pirna/human-pirna.fa", ["piR-hsa-9"]),
    ("bowtie2",      "gtf",       "human-ncrna",    "gtf/human.ncRNA.fa",   ["ENST00000786877|H19"]),
    ("bowtie2",      "gtf",       "human-cdna",     "gtf/human.cdna.fa",    ["ENST00000859953|ETS2",
                                                                             "ENST00001010609|BTG3",
                                                                             "ENST00001110240|MRPL39"]),
    ("bowtie2",      "genome",    "human-genome",   "CHR21",                "window"),
]

# Synthetic genome features (not in the production BED), in window coordinates, for the
# hard cases of make_reads.py (2026-10-05). Both lie in unmasked sequence.
#   SYNTH-SNORD-minus: a minus-strand snoRNA (rank 5) inside the ETS2 gene body (pre-mRNA, rank 1, +),
#                      for reads labelled by the higher-priority feature on the other strand (C22, C23)
#   SYNTH-LNC-edge:    an intergenic plus-strand lncRNA, for reads that only touch its start (C24, C25)
SYNTHETIC_BED_ROWS = [
    [GENOME_CONTIG, 38_000, 38_300, "SYNTH-SNORD-minus", "snoRNA", "-"],
    [GENOME_CONTIG, 60_000, 60_300, "SYNTH-LNC-edge", "lncRNA", "+"],
]


def read_fasta(path):
    """Return {name: sequence}; name = header up to the first whitespace. Reads .gz too."""
    opener = gzip.open if path.endswith(".gz") else open
    records, name, parts = {}, None, []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts)
    return records


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n{seq}\n")


def check_tier_order(assets):
    """The mini set must use exactly the production tiers, in the production order."""
    production = []
    dbpaths = os.path.join(assets, "human.dbpaths.csv")              # since 2026-10-08
    if not os.path.exists(dbpaths):
        dbpaths = os.path.join(assets, "human.local.dbpaths.csv")    # older references folders
    with open(dbpaths) as f:
        for line in f:
            if line.strip():
                map_cmd, path, _ = line.strip().split(",")
                production.append((map_cmd, os.path.basename(path)))
    mine = [(t[0], t[2]) for t in TIERS]
    if production != mine:
        sys.exit(f"ERROR: tier list differs from production dbpaths.\n production: {production}\n mini:       {mine}")


def select_records(source, keep, args):
    if source == "UNIVEC":
        records = read_fasta(args.univec)
    elif source == "RIBOSPIKE":
        records = read_fasta(os.path.join(REPO, "buildDB", "inputs", "ribospike.fa"))
    elif source == "CHR21":
        chrom = read_fasta(args.chr21)[GENOME_CONTIG]
        return {GENOME_CONTIG: chrom[WINDOW_START:WINDOW_END]}   # soft-masking (lower case) kept
    else:
        records = read_fasta(os.path.join(args.assets, source))

    if keep == "all":
        return records
    if keep == "first 3":
        return dict(list(records.items())[:3])
    missing = [k for k in keep if k not in records]
    if missing:
        sys.exit(f"ERROR: records not found in {source}: {missing}")
    return {k: records[k] for k in keep}


def mini_bed_rows(assets, kept_contigs):
    """Production BED rows for the kept contigs, plus genome rows inside the window (shifted)."""
    rows = []
    with open(os.path.join(assets, "human.bed")) as f:
        for line in f:
            chrom, start, end, feature, biotype, strand = line.rstrip("\n").split("\t")
            start, end = int(start), int(end)
            if chrom in kept_contigs and chrom != GENOME_CONTIG:
                rows.append([chrom, start, end, feature, biotype, strand])
            elif chrom == GENOME_CONTIG and end > WINDOW_START and start < WINDOW_END:
                new_start = max(start, WINDOW_START) - WINDOW_START
                new_end = min(end, WINDOW_END) - WINDOW_START
                rows.append([chrom, new_start, new_end, feature, biotype, strand])
    return rows + SYNTHETIC_BED_ROWS


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assets", default="docker/assets", help="production assets folder (read only)")
    ap.add_argument("--chr21", required=True, help="chromosome 21 FASTA (Ensembl 116, .fa or .fa.gz)")
    ap.add_argument("--univec", required=True, help="UniVec FASTA (bowtie-inspect of the production index)")
    ap.add_argument("--out", default="test/synthetic_qc/bundle/references")
    args = ap.parse_args()

    check_tier_order(args.assets)
    os.makedirs(args.out, exist_ok=True)

    tier_rows, kept_contigs = [], set()
    for map_cmd, folder, index_name, source, keep in TIERS:
        records = select_records(source, keep, args)
        kept_contigs.update(records)
        os.makedirs(os.path.join(args.out, folder), exist_ok=True)
        fasta = f"{folder}/{index_name}.fa"
        write_fasta(os.path.join(args.out, fasta), records)
        tier_rows.append(f"{map_cmd},{fasta},{index_name}")
        print(f"{index_name:16s} {len(records):3d} record(s)")

    with open(os.path.join(args.out, "tiers.csv"), "w") as f:     # mapping order
        f.write("map_cmd,fasta,index\n" + "\n".join(tier_rows) + "\n")
    with open(os.path.join(args.out, "human.bed"), "w") as f:
        for row in mini_bed_rows(args.assets, kept_contigs):
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"wrote FASTA files, tiers.csv and human.bed in {args.out}")


if __name__ == "__main__":
    main()
