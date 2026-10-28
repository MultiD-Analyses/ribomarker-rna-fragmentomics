#!/usr/bin/env python3
"""
make_reads.py -- synthetic RiboMarker reads with a known fate, plus the truth table.

Read layout (RiboMarker protocol, Appendix D p. 20, and a raw-read check on
a real sample, 2026-09-29):
  R1 = 1 extra base + insert + 3' adapter and the rest of the library (76 nt)
  R2 = reverse complement of the insert + R2 adapter and the rest (76 nt)
All bases have quality 'I' (Phred 40), so no read is removed for quality.

Every read belongs to one CATEGORY with a known fate (see the table below).
The expected numbers in truth.tsv are derived from these categories only,
never from running the pipeline.

Usage (from the repository root, after make_mini_assets.py):
  python test/synthetic_qc/make_reads.py --references test/synthetic_qc/bundle/references \
      --out test/synthetic_qc/bundle

Outputs in --out:
  fastq/<sample>_R1_001.fastq.gz, fastq/<sample>_R2_001.fastq.gz
  fastq/samples_pe.csv, fastq/samples_se.csv
  truth/reads_truth.tsv  every read: what it is and its expected fate (pipeline-independent)
  truth/truth.tsv        expected values in THIS pipeline's output tables, per mode and sample
  truth/categories.tsv   the read categories and their counts

Hard cases (added 2026-10-05): C19-C31 test the
classification rules on reads that are not exact unique matches (near-multimappers,
mismatches, opposite-strand features, feature edges, spike-in offsets and fragments,
fastp-filtered pairs). Their expected fate follows the documented rules
(see the user guide); where it also depends on
aligner or fastp defaults, the category text says so. Sample S3 holds only the
exploratory categories C32-C35: their outcome is reported, not checked.
"""
import argparse
import gzip
import os
import random
import re
import sys
from collections import defaultdict

SEED = 20260929
SAMPLES = ("S1", "S2", "S3")
EXPLORATORY_SAMPLES = ("S3",)   # outcome reported by check_truth.py, no expected values beyond trimming
READ_LEN = 76
MARGIN = 20   # fragments stay >= 20 nt away from contig ends (N pad of 10 + 2-nt clipping + safety)
K = 12        # genome fragments may not share any 12-mer with the transcript references (see Screen)
TRANSCRIPT_TIERS = ("human-ncrna", "human-cdna")

# Library sequence after the insert, as observed in raw reads (then poly-G on this instrument).
R1_TAIL = "TGGAATTCTCGGGTGCCAAGG" + "AACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG" + "A" * 10
R2_TAIL = "TGATCGTCGGACT" + "GTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT" + "A" * 10

# ---------------------------------------------------------------------------
# Read categories. n = reads per sample (S1, S2).
#   draw      how the insert is made:
#               ("empty",) | ("random", length) | ("spike",)
#               ("contig", index, contig, length, expected_occurrences, orientation)
#               ("whole", index, contig)          the full reference sequence
#               ("intron", length) | ("intergenic", length)   from the genome window
#   lost      removed before mapping: "dimer" (empty insert), "short" (< 15 nt),
#             "merge" (paired-end only: pair cannot be merged)
#   mapping   "unique" | "multi" | "unmapped" | "spike"
#   tier      index the read aligns to
#   rna_class row in *-classcounts.csv  (class id; None: removed by remove_multimap / remove_antisense)
#   gene      row in *-genecounts.csv   (None: not a mature gene)
# Optional (hard cases):
#   n           a third number = reads in S3 (missing = 0)
#   removed     "multimap" | "antisense": removed by policy (assign.removed_*)
#   alignments  reported alignments in the tier, if not given by `mapping`
#   mismatches  insert positions (0-based) changed to another base after drawing
#   leading_n   number of leading insert bases replaced by N
#   quality     Phred character for all bases (default "I" = Phred 40)
#   lost        also "fastp": paired-end pair removed by fastp's filters (single-end: kept)
#   exploratory True: S3 only, outcome reported, not checked
#   draw        also ("near_multi", index, contig, other contig, length): exact to contig,
#                   one mismatch to the other copy (same length, aligned position by position)
#               ("in_feature", feature, length, "plus"|"minus"): genome read inside a BED feature
#               ("edge", feature, overlap, length): genome read whose last `overlap` bases
#                   are the first bases of a plus-strand feature
#               ("spike", trim5, trim3): every spike-in with bases removed at its ends
#               ("spike_piece", length, spike length): internal piece of a spike-in core
# ---------------------------------------------------------------------------
CDNA_ETS2, CDNA_BTG3, CDNA_MRPL39 = "ENST00000859953|ETS2", "ENST00001010609|BTG3", "ENST00001110240|MRPL39"

CATEGORIES = [
    dict(id="C01", what="adapter dimer (empty insert)", n=(101, 83), draw=("empty",), lost="dimer"),
    dict(id="C02", what="10-nt insert (below the 15-nt minimum)", n=(53, 47), draw=("random", 10), lost="short"),
    dict(id="C03", what="150-nt BTG3 insert: paired-end mates do not overlap (not merged); single-end R1 maps",
         n=(37, 29), draw=("contig", "human-cdna", CDNA_BTG3, 150, 1, "sense"), lost="merge",
         mapping="unique", tier="human-cdna", rna_class="protein_coding", gene="BTG3"),
    dict(id="C04", what="spike-in, n reads for each of the 20 spikes", n=(7, 0), draw=("spike", 0, 0),
         mapping="spike", tier="ribospike"),
    dict(id="C05", what="UniVec fragment (UniVec contigs are not in the BED; class UniVec since 0.3.0)", n=(29, 31),
         draw=("contig", "univec_bt1", None, 30, 1, "sense"),
         mapping="unique", tier="univec_bt1", rna_class="UniVec", gene=None),
    dict(id="C06", what="rRNA fragment", n=(41, 43),
         draw=("contig", "human-rrna_bt1", "Homo_sapiens_18S_1", 30, 1, "sense"),
         mapping="unique", tier="human-rrna_bt1", rna_class="rRNA", gene="Homo_sapiens_18S_1"),
    dict(id="C07", what="tRNA fragment, unique", n=(47, 53),
         draw=("contig", "human-trnas", "Homo_sapiens_tRNA-Ala-CGC-3", 30, 1, "sense"),
         mapping="unique", tier="human-trnas", rna_class="tRNA", gene="Homo_sapiens_tRNA-Ala-CGC-3"),
    dict(id="C08", what="tRNA fragment shared by Leu-CAG-1 and -2 (multimapper, removed)", n=(31, 37),
         draw=("contig", "human-trnas", "Homo_sapiens_tRNA-Leu-CAG-1", 30, 2, "sense"),
         mapping="multi", removed="multimap", tier="human-trnas", rna_class=None, gene=None),
    dict(id="C09", what="mature miRNA, full length", n=(59, 61), draw=("whole", "human-mirna", "Hsa-Mir-21_5p"),
         mapping="unique", tier="human-mirna", rna_class="miRNA", gene="Hsa-Mir-21_5p"),
    dict(id="C10", what="piRNA, full length", n=(43, 41), draw=("whole", "human-pirna", "piR-hsa-9"),
         mapping="unique", tier="human-pirna", rna_class="piRNA", gene="piR-hsa-9"),
    dict(id="C11", what="lncRNA fragment", n=(61, 67),
         draw=("contig", "human-ncrna", "ENST00000786877|H19", 30, 1, "sense"),
         mapping="unique", tier="human-ncrna", rna_class="lncRNA", gene="H19"),
    dict(id="C12", what="mRNA fragment, ETS2", n=(67, 71), draw=("contig", "human-cdna", CDNA_ETS2, 30, 1, "sense"),
         mapping="unique", tier="human-cdna", rna_class="protein_coding", gene="ETS2"),
    dict(id="C13", what="mRNA fragment, BTG3", n=(71, 73), draw=("contig", "human-cdna", CDNA_BTG3, 30, 1, "sense"),
         mapping="unique", tier="human-cdna", rna_class="protein_coding", gene="BTG3"),
    dict(id="C14", what="mRNA fragment, ETS2, antisense (removed)", n=(19, 23),
         draw=("contig", "human-cdna", CDNA_ETS2, 30, 1, "antisense"),
         mapping="unique", removed="antisense", tier="human-cdna", rna_class=None, gene=None),
    dict(id="C15", what="ETS2 gene body, not in the cDNA (pre-mRNA)", n=(23, 29), draw=("intron", 30),
         mapping="unique", tier="human-genome", rna_class="pre-mRNA", gene=None),
    dict(id="C16", what="intergenic genome fragment", n=(17, 19), draw=("intergenic", 30),
         mapping="unique", tier="human-genome", rna_class="genome_unannotated", gene=None),
    dict(id="C17", what="random 30-nt sequence, absent from every reference", n=(13, 11), draw=("random", 30),
         mapping="unmapped", tier=None, rna_class="unmapped", gene=None),
    dict(id="C18", what="15-nt MRPL39 fragment (exactly the minimum length)", n=(11, 13),
         draw=("contig", "human-cdna", CDNA_MRPL39, 15, 1, "sense"),
         mapping="unique", tier="human-cdna", rna_class="protein_coding", gene="MRPL39"),
    # ---- hard cases (2026-10-05) ----
    dict(id="C19", what="tRNA fragment exact to Leu-CAG-1, 1 mismatch to Leu-CAG-2: unique by best score (kept), "
                        "two reported alignments (aligner-dependent)",
         n=(23, 19), draw=("near_multi", "human-trnas", "Homo_sapiens_tRNA-Leu-CAG-1", "Homo_sapiens_tRNA-Leu-CAG-2", 30),
         mapping="unique", alignments=2, tier="human-trnas", rna_class="tRNA", gene="Homo_sapiens_tRNA-Leu-CAG-1"),
    dict(id="C20", what="mRNA fragment, MRPL39, 1 mismatch (Bowtie 2 tier)", n=(17, 13),
         draw=("contig", "human-cdna", CDNA_MRPL39, 30, 1, "sense"), mismatches=(15,),
         mapping="unique", tier="human-cdna", rna_class="protein_coding", gene="MRPL39"),
    dict(id="C21", what="rRNA fragment, 1 mismatch (Bowtie 1 tier)", n=(13, 17),
         draw=("contig", "human-rrna_bt1", "Homo_sapiens_18S_1", 30, 1, "sense"), mismatches=(15,),
         mapping="unique", tier="human-rrna_bt1", rna_class="rRNA", gene="Homo_sapiens_18S_1"),
    dict(id="C22", what="ETS2 gene body, plus strand, inside a minus-strand snoRNA: the higher-rank snoRNA is chosen "
                        "whatever its strand (§6.1), so antisense (removed)",
         n=(19, 17), draw=("in_feature", "SYNTH-SNORD-minus", 30, "plus"), orientation="antisense",
         mapping="unique", removed="antisense", tier="human-genome", rna_class=None, gene=None),
    dict(id="C23", what="same snoRNA region, minus strand: sense to the snoRNA (kept)", n=(15, 21),
         draw=("in_feature", "SYNTH-SNORD-minus", 30, "minus"), orientation="sense",
         mapping="unique", tier="human-genome", rna_class="snoRNA", gene="SYNTH-SNORD-minus"),
    dict(id="C24", what="intergenic read whose last base is the first base of a lncRNA (1-nt overlap, "
                        "removed by the 2-nt clipping, §6.1): unannotated",
         n=(11, 7), draw=("edge", "SYNTH-LNC-edge", 1, 30), orientation="",
         mapping="unique", tier="human-genome", rna_class="genome_unannotated", gene=None),
    dict(id="C25", what="intergenic read overlapping the same lncRNA by 5 nt: assigned to it", n=(7, 11),
         draw=("edge", "SYNTH-LNC-edge", 5, 30), orientation="sense",
         mapping="unique", tier="human-genome", rna_class="lncRNA", gene="SYNTH-LNC-edge"),
    dict(id="C26", what="spike-in without its last base (offset -1: full-length, not exact), n per spike", n=(2, 0),
         draw=("spike", 0, 1), mapping="spike", tier="ribospike"),
    dict(id="C27", what="spike-in without its first and last base (offset -2: full-length, not exact), n per spike",
         n=(1, 0), draw=("spike", 1, 1), mapping="spike", tier="ribospike"),
    dict(id="C28", what="spike-in without its last 3 bases (offset -3: fragment), n per spike", n=(1, 0),
         draw=("spike", 0, 3), mapping="spike", tier="ribospike"),
    dict(id="C29", what="20-nt internal piece of a 60-nt spike-in core (fragment), n per 60-nt spike", n=(3, 0),
         draw=("spike_piece", 20, 60), mapping="spike", tier="ribospike"),
    dict(id="C30", what="ETS2 fragment, all bases Phred 10: removed by fastp's quality filter in paired-end mode "
                        "(fastp defaults); single-end has no quality filter and maps it",
         n=(13, 11), draw=("contig", "human-cdna", CDNA_ETS2, 30, 1, "sense"), quality="+", lost="fastp",
         mapping="unique", tier="human-cdna", rna_class="protein_coding", gene="ETS2"),
    dict(id="C31", what="60-nt ETS2 fragment, first 6 bases N: removed by fastp's N filter in paired-end mode "
                        "(fastp defaults); single-end maps it (6 N within Bowtie 2's default N ceiling, aligner-dependent)",
         n=(9, 7), draw=("contig", "human-cdna", CDNA_ETS2, 60, 1, "sense"), leading_n=6, lost="fastp",
         mapping="unique", tier="human-cdna", rna_class="protein_coding", gene="ETS2"),
    # ---- exploratory, sample S3 only: outcome reported, not checked (sub-seed mismatches) ----
    dict(id="C32", what="16-nt BTG3 fragment, exact (control)", n=(0, 0, 20), exploratory=True,
         draw=("contig", "human-cdna", CDNA_BTG3, 16, 1, "sense"), mapping="exploratory", gene="BTG3"),
    dict(id="C33", what="16-nt ETS2 fragment, 1 mismatch at position 9", n=(0, 0, 20), exploratory=True,
         draw=("contig", "human-cdna", CDNA_ETS2, 16, 1, "sense"), mismatches=(8,), mapping="exploratory", gene="ETS2"),
    dict(id="C34", what="17-nt MRPL39 fragment, 1 mismatch at position 9", n=(0, 0, 20), exploratory=True,
         draw=("contig", "human-cdna", CDNA_MRPL39, 17, 1, "sense"), mismatches=(8,), mapping="exploratory",
         gene="MRPL39"),
    dict(id="C35", what="17-nt H19 fragment, exact (control)", n=(0, 0, 20), exploratory=True,
         draw=("contig", "human-ncrna", "ENST00000786877|H19", 17, 1, "sense"), mapping="exploratory", gene="H19"),
]


def n_of(cat, sample_idx):
    """Reads of a category in one sample; a missing S3 number means 0."""
    return cat["n"][sample_idx] if sample_idx < len(cat["n"]) else 0

# What each read IS, independent of any pipeline (for reads_truth.tsv): (biotype, feature).
# feature None = the read's source contig.
IDENTITY = {
    "C01": ("adapter dimer", ""),            "C02": ("insert too short", ""),
    "C03": ("protein_coding", "BTG3"),       "C04": ("spike-in", None),
    "C05": ("vector (UniVec)", None),        "C06": ("rRNA", "Homo_sapiens_18S_1"),
    "C07": ("tRNA", "Homo_sapiens_tRNA-Ala-CGC-3"),
    "C08": ("tRNA", "Homo_sapiens_tRNA-Leu-CAG-1|Homo_sapiens_tRNA-Leu-CAG-2"),
    "C09": ("miRNA", "Hsa-Mir-21_5p"),       "C10": ("piRNA", "piR-hsa-9"),
    "C11": ("lncRNA", "H19"),                "C12": ("protein_coding", "ETS2"),
    "C13": ("protein_coding", "BTG3"),       "C14": ("protein_coding", "ETS2"),
    "C15": ("intron (gene body)", "ETS2"),   "C16": ("intergenic", ""),
    "C17": ("random sequence", ""),          "C18": ("protein_coding", "MRPL39"),
    "C19": ("tRNA", "Homo_sapiens_tRNA-Leu-CAG-1"), "C20": ("protein_coding", "MRPL39"),
    "C21": ("rRNA", "Homo_sapiens_18S_1"),
    "C22": ("intron (gene body) inside a synthetic minus-strand snoRNA", "ETS2|SYNTH-SNORD-minus"),
    "C23": ("snoRNA (synthetic feature)", "SYNTH-SNORD-minus"),
    "C24": ("intergenic, 1 nt into a synthetic lncRNA", ""),
    "C25": ("lncRNA (synthetic feature), 5 nt overlap", "SYNTH-LNC-edge"),
    "C26": ("spike-in", None), "C27": ("spike-in", None), "C28": ("spike-in", None), "C29": ("spike-in", None),
    "C30": ("protein_coding", "ETS2"),       "C31": ("protein_coding", "ETS2"),
    "C32": ("protein_coding", "BTG3"),       "C33": ("protein_coding", "ETS2"),
    "C34": ("protein_coding", "MRPL39"),     "C35": ("lncRNA", "H19"),
}
ALIGNMENTS = {"unique": 1, "multi": 2, "spike": 1, "unmapped": 0}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    records, name, parts = {}, None, []
    with open(path) as f:
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


def load_references(references):
    """{index name: {contig: sequence}} from the bundle's reference FASTA files, and the BED."""
    refs = {}
    with open(os.path.join(references, "tiers.csv")) as f:
        next(f)                                              # header: map_cmd,fasta,index
        for line in f:
            if line.strip():
                _, fasta, index = line.strip().split(",")
                refs[index] = read_fasta(os.path.join(references, fasta))
    bed = []
    with open(os.path.join(references, "human.bed")) as f:
        for line in f:
            chrom, start, end, feature, biotype, strand = line.rstrip("\n").split("\t")
            bed.append((chrom, int(start), int(end), feature, biotype, strand))
    return refs, bed


class Screen:
    """Counts exact occurrences of a sequence, on both strands, in all mini references.

    Genome fragments are also checked against the transcript tiers (ncRNA, cDNA),
    which are searched BEFORE the genome: a fragment that overlaps an exon by most
    of its length aligns there with a mismatch or two. The first synthetic run
    (2026-09-29) had one such read per sample (29 exon nt + 1 intron nt, NM 1).
    Any shared 12-mer rejects the fragment.
    """

    def __init__(self, refs):
        self.all_seq = "|".join(seq.upper() for tier in refs.values() for seq in tier.values())
        self.transcript_kmers = set()
        for tier in TRANSCRIPT_TIERS:
            for seq in refs[tier].values():
                for strand in (seq.upper(), revcomp(seq.upper())):
                    self.transcript_kmers.update(strand[i:i + K] for i in range(len(strand) - K + 1))

    def shares_transcript_kmer(self, seq):
        seq = seq.upper()
        return any(seq[i:i + K] in self.transcript_kmers for i in range(len(seq) - K + 1))

    def occurrences(self, seq):
        seq = seq.upper()
        n = self.all_seq.count(seq)
        rc = revcomp(seq)
        return n if rc == seq else n + self.all_seq.count(rc)


def draw_insert(cat, refs, bed, screen, rng):
    """Return (insert, source contig, start) for one read of this category."""
    kind = cat["draw"][0]
    if kind == "empty":
        return "", "-", 0
    if kind == "random":
        length = cat["draw"][1]
        for _ in range(1000):
            seq = random_seq(length, rng)
            if cat.get("mapping") is None or screen.occurrences(seq) == 0:
                return seq, "random", 0
    if kind == "whole":
        _, index, contig = cat["draw"]
        return refs[index][contig].strip("N"), contig, 10
    if kind == "contig":
        _, index, contig, length, expected, orientation = cat["draw"]
        contig = contig or next(iter(refs[index]))
        ref = refs[index][contig]
        for _ in range(1000):
            start = rng.randint(MARGIN, len(ref) - MARGIN - length)
            seq = ref[start:start + length].upper()
            parts = [seq] if length <= READ_LEN else [seq, seq[:70]]   # C03: also the part that maps single-end
            if "N" not in seq and all(screen.occurrences(p) == expected for p in parts):
                return (seq if orientation == "sense" else revcomp(seq)), contig, start
    if kind == "near_multi":
        _, index, contig, other, length = cat["draw"]
        ref, other_ref = refs[index][contig].upper(), refs[index][other].upper()
        starts = list(range(10, len(ref) - 10 - length + 1))          # inside the 10-nt N padding
        rng.shuffle(starts)
        for start in starts:
            seq = ref[start:start + length]
            differences = sum(a != b for a, b in zip(seq, other_ref[start:start + length]))
            if "N" not in seq and differences == 1 and screen.occurrences(seq) == 1:
                return seq, contig, start
    if kind in ("in_feature", "edge"):
        chrom, genome = next(iter(refs["human-genome"].items()))
        feature = [r for r in bed if r[3] == cat["draw"][1]][0]
        if kind == "in_feature":
            _, _, length, strand = cat["draw"]
            starts = [rng.randint(feature[1] + MARGIN, feature[2] - MARGIN - length) for _ in range(1000)]
        else:
            _, _, overlap, length = cat["draw"]
            strand = "plus"
            starts = [feature[1] + overlap - length]     # the read's last `overlap` bases are inside the feature
        for start in starts:
            seq = genome[start:start + length]
            if (seq.isupper() and "N" not in seq and screen.occurrences(seq) == 1
                    and not screen.shares_transcript_kmer(seq)):
                return (seq if strand == "plus" else revcomp(seq)), chrom, start
    if kind in ("intron", "intergenic"):
        length = cat["draw"][1]
        chrom, genome = next(iter(refs["human-genome"].items()))
        body = [r for r in bed if r[0] == chrom and r[3] == "ETS2_intron"][0]
        annotated = [(r[1], r[2]) for r in bed if r[0] == chrom]
        inside_body = [(r[1], r[2]) for r in bed if r[0] == chrom and r[3] != "ETS2_intron"]   # synthetic features
        for _ in range(1000):
            if kind == "intron":
                start = rng.randint(body[1] + MARGIN, body[2] - MARGIN - length)
                if any(start < e + MARGIN and start + length > s - MARGIN for s, e in inside_body):
                    continue
            else:
                start = rng.randint(MARGIN, len(genome) - MARGIN - length)
                if any(start < e + 1000 and start + length > s - 1000 for s, e in annotated):
                    continue
            seq = genome[start:start + length]
            if (seq.isupper() and "N" not in seq and screen.occurrences(seq) == 1   # upper case = not repeat-masked
                    and not screen.shares_transcript_kmer(seq)):
                return seq, chrom, start
    sys.exit(f"ERROR: could not draw a valid insert for {cat['id']} after 1000 tries")


def draw_read_insert(cat, refs, bed, screen, rng):
    """draw_insert plus the designed changes: mismatches (the changed read must match
    nothing exactly) and leading N."""
    for _ in range(1000):
        insert, source, start = draw_insert(cat, refs, bed, screen, rng)
        for pos in cat.get("mismatches", ()):
            insert = insert[:pos] + rng.choice([b for b in "ACGT" if b != insert[pos]]) + insert[pos + 1:]
        if cat.get("mismatches") and screen.occurrences(insert) > 0:
            continue
        if cat.get("leading_n"):
            insert = "N" * cat["leading_n"] + insert[cat["leading_n"]:]
        return insert, source, start
    sys.exit(f"ERROR: could not draw a valid changed insert for {cat['id']} after 1000 tries")


def spike_targets(cat, spikes):
    """The spikes a spike category draws from: all 20, or those of one length (spike_piece)."""
    if cat["draw"][0] == "spike_piece":
        return [sp for sp in spikes if sp[1] == cat["draw"][2]]
    return spikes


def spike_insert(cat, spike, screen, rng):
    """(insert, start in the spike contig) for one read of a spike category."""
    name, length, spike_type, core = spike
    if cat["draw"][0] == "spike_piece":
        piece = cat["draw"][1]
        for _ in range(1000):
            start = rng.randint(0, len(core) - piece)
            if screen.occurrences(core[start:start + piece]) == 1:
                return core[start:start + piece], start + 4
        sys.exit(f"ERROR: no unique {piece}-nt piece in {name}")
    _, trim5, trim3 = cat["draw"]
    full = random_seq(4, rng) + core + random_seq(4, rng)
    return full[trim5:len(full) - trim3], trim5


def spike_class(cat, length):
    """(full-length?, offset) of a spike read: full-length if at least designed length - 2
    pieces are fragments."""
    if cat["draw"][0] == "spike_piece":
        return False, None
    offset = -(cat["draw"][1] + cat["draw"][2])
    return offset >= -2, offset


def spike_list(refs):
    """[(contig name, designed length, type, core)] for the 20 spikes; the core excludes the 4N ends."""
    out = []
    for name, seq in refs["ribospike"].items():
        length, spike_type = re.search(r"SR_(\d+)_4N-T(\d+)", name).groups()
        out.append((name, int(length), f"T{spike_type}", seq.strip("N")))
    return out


def random_seq(length, rng):
    return "".join(rng.choice("ACGT") for _ in range(length))


def make_pair(insert, rng):
    r1 = rng.choice("ACGT") + insert + R1_TAIL
    r2 = revcomp(insert) + R2_TAIL
    return (r1 + "G" * READ_LEN)[:READ_LEN], (r2 + "G" * READ_LEN)[:READ_LEN]


def write_fastq_gz(path, records):
    """Deterministic gzip (mtime 0), so reruns give byte-identical files. records: (name, mate, seq, quality char)."""
    text = "".join(f"@{name} {mate}:N:0:TTAGGC\n{seq}\n+\n{q * len(seq)}\n" for name, mate, seq, q in records)
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(text.encode())


def expected_numbers(sample_idx, spikes):
    """All expected numbers for one sample, both modes: {(mode, table, key): (value, categories)}.
    Table "qc" holds QC metric ids (qc_metrics.json) that follow from the design;
    "tier_aligned_multi" reads with >= 2 reported alignments per tier; "exploratory"
    the S3 categories (key category:gene, value reads), whose outcome is only reported."""
    truth = defaultdict(lambda: [0, set()])
    exploratory_sample = SAMPLES[sample_idx] in EXPLORATORY_SAMPLES

    def add(mode, table, key, n, cat_id):
        truth[(mode, table, key)][0] += n
        truth[(mode, table, key)][1].add(cat_id)

    for mode in ("pe", "se"):                    # QC metrics that are 0 unless a category adds to them
        keys = [] if exploratory_sample else ["assign.removed_multimap", "assign.removed_antisense",
                                               "spike.fragment_reads"]
        if mode == "pe":
            keys += ["trim.fastp_filtered_pairs", "trim.fastp_low_quality_reads",
                     "trim.fastp_too_many_n_reads", "trim.unmerged_pairs"]
        for key in keys:
            add(mode, "qc", key, 0, "design")
    offsets = defaultdict(int)                   # (mode, designed length, offset): full-length spike reads

    for cat in CATEGORIES:
        n = n_of(cat, sample_idx)
        if n == 0:
            continue
        for mode in ("pe", "se"):
            lost = cat.get("lost")
            if lost in ("merge", "fastp") and mode == "se":
                lost = None
            if cat.get("mapping") == "spike":
                units = [(n, spike) for spike in spike_targets(cat, spikes)]   # n reads of every spike
            else:
                units = [(n, None)]
            total = sum(u[0] for u in units)
            add(mode, "input", "reads", total, cat["id"])
            # trimming (paired-end: trim log counts pairs < 15 nt as 'dimer_reads')
            if mode == "pe":
                if lost in ("dimer", "short"):
                    add(mode, "trimstats", "dimer_reads", total, cat["id"])
                else:
                    add(mode, "trimstats", "trimmed_reads", total, cat["id"])
                    if lost not in ("merge", "fastp"):
                        add(mode, "trimstats", "merged_reads", total, cat["id"])
                if lost == "merge":
                    add(mode, "qc", "trim.unmerged_pairs", total, cat["id"])
                if lost == "fastp":
                    add(mode, "qc", "trim.fastp_filtered_pairs", total, cat["id"])
                    if cat.get("quality"):               # every base below Phred 15
                        add(mode, "qc", "trim.fastp_low_quality_reads", 2 * total, cat["id"])   # R1 and R2
                    if cat.get("leading_n", 0) > 5:
                        add(mode, "qc", "trim.fastp_too_many_n_reads", 2 * total, cat["id"])
            else:
                if lost == "dimer":
                    add(mode, "trimstats", "dimer_reads", total, cat["id"])
                if lost in ("dimer", "short"):
                    add(mode, "trimstats", "<15nt_reads", total, cat["id"])
                else:
                    add(mode, "trimstats", "trimmed_reads", total, cat["id"])
            if lost:
                continue
            if cat.get("exploratory"):
                add(mode, "exploratory", f"{cat['id']}:{cat['gene']}", total, cat["id"])
                continue
            # mapping and assignment
            mapping = cat["mapping"]
            if mapping == "unmapped":
                add(mode, "mappingstats", "Unmapped", total, cat["id"])
                add(mode, "tier_aligned", "final_unaligned", total, cat["id"])
            else:
                add(mode, "mappingstats", "Mapped", total, cat["id"])
                add(mode, "tier_aligned", cat["tier"], total, cat["id"])
                if mapping == "unique":
                    add(mode, "mappingstats", "Unique", total, cat["id"])
                if mapping == "multi":
                    add(mode, "mappingstats", "Multi", total, cat["id"])
                if cat.get("alignments", ALIGNMENTS[mapping]) >= 2:
                    add(mode, "tier_aligned_multi", cat["tier"], total, cat["id"])
            add(mode, "mappingstats", "Total", total, cat["id"])
            if cat.get("removed"):
                add(mode, "qc", f"assign.removed_{cat['removed']}", total, cat["id"])
            if mapping == "spike":
                for n_spike, (name, length, spike_type, _) in units:
                    full_length, offset = spike_class(cat, length)
                    if not full_length:
                        add(mode, "classcounts", "RiboSpike-Fragment", n_spike, cat["id"])
                        add(mode, "qc", "spike.fragment_reads", n_spike, cat["id"])
                        continue
                    add(mode, "classcounts", f"RiboSpike-{spike_type}", n_spike, cat["id"])
                    add(mode, "ribospike", f"RiboSpike-{length}-{spike_type}", n_spike, cat["id"])
                    add(mode, "spikelengths", f"RiboSpike-{spike_type}:{length + offset}", n_spike, cat["id"])
                    add(mode, "qc", f"spike.{length}.reads", n_spike, cat["id"])
                    offsets[(mode, length, offset)] += n_spike
                    if offset == 0:              # designed length, whole spike: an exact read
                        add(mode, "qc", f"spike.{length}.exact_reads", n_spike, cat["id"])
                continue
            if cat.get("rna_class"):
                add(mode, "classcounts", cat["rna_class"], total, cat["id"])
            if cat.get("gene"):
                add(mode, "genecounts", cat["gene"], total, cat["id"])

    for mode, length in sorted({(m, l) for m, l, _ in offsets}):
        by_offset = {o: v for (m, l, o), v in offsets.items() if (m, l) == (mode, length)}
        modal = max(by_offset, key=by_offset.get)
        add(mode, "qc", f"spike.{length}.modal_offset", modal, "design")
        add(mode, "qc", f"spike.{length}.modal_offset_reads", by_offset[modal], "design")
    return truth


def fate(cat, mode):
    """Expected fate of a read of this category: lost:<stage>, mapped or unmapped."""
    lost = cat.get("lost")
    if lost in ("merge", "fastp") and mode == "se":
        lost = None
    if lost:
        return f"lost:{lost}"
    if cat.get("exploratory"):
        return "exploratory"
    return "unmapped" if cat["mapping"] == "unmapped" else "mapped"


def orientation(cat):
    if "orientation" in cat:
        return cat["orientation"]
    draw = cat["draw"]
    if draw[0] in ("empty", "random", "intergenic"):
        return ""
    return draw[5] if draw[0] == "contig" else "sense"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--references", default="test/synthetic_qc/bundle/references")
    ap.add_argument("--out", default="test/synthetic_qc/bundle")
    args = ap.parse_args()

    rng = random.Random(SEED)
    refs, bed = load_references(args.references)
    cat_by_id = {c["id"]: c for c in CATEGORIES}
    screen = Screen(refs)
    spikes = spike_list(refs)
    fastq_dir = os.path.join(args.out, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)

    key_rows, truth_rows = [], []
    for s_idx, sample in enumerate(SAMPLES):
        pairs = []
        for cat in CATEGORIES:
            n = n_of(cat, s_idx)
            jobs = [(cat["id"], None)] * n if cat.get("mapping") != "spike" else \
                   [(cat["id"], sp) for sp in spike_targets(cat, spikes) for _ in range(n)]
            for cat_id, spike in jobs:
                if spike:
                    insert, start = spike_insert(cat, spike, screen, rng)
                    source = spike[0]
                else:
                    insert, source, start = draw_read_insert(cat, refs, bed, screen, rng)
                pairs.append((cat_id, source, start, insert, make_pair(insert, rng)))
        rng.shuffle(pairs)

        r1, r2 = [], []
        for i, (cat_id, source, start, insert, (s1, s2)) in enumerate(pairs, 1):
            name = f"SYN_{sample}_{cat_id}_{i:05d}"
            cat = cat_by_id[cat_id]
            quality = cat.get("quality", "I")
            r1.append((name, 1, s1, quality))
            r2.append((name, 2, s2, quality))
            biotype, feature = IDENTITY[cat_id]
            key_rows.append([name, sample, cat_id, biotype, source if feature is None else feature,
                             orientation(cat), len(insert), insert, source, start,
                             fate(cat, "pe"), fate(cat, "se"), cat.get("tier") or "",
                             cat.get("alignments", ALIGNMENTS.get(cat.get("mapping"), ""))])
        write_fastq_gz(os.path.join(fastq_dir, f"{sample}_R1_001.fastq.gz"), r1)
        write_fastq_gz(os.path.join(fastq_dir, f"{sample}_R2_001.fastq.gz"), r2)

        for (mode, table, key), (value, cats) in sorted(expected_numbers(s_idx, spikes).items()):
            truth_rows.append([mode, sample, table, key, value, ",".join(sorted(cats))])
        print(f"{sample}: {len(pairs)} read pairs")

    with open(os.path.join(fastq_dir, "samples_pe.csv"), "w") as f:
        f.writelines(f"{s},synthetic,{s}_R1_001.fastq.gz,{s}_R2_001.fastq.gz\n" for s in SAMPLES)
    with open(os.path.join(fastq_dir, "samples_se.csv"), "w") as f:
        f.writelines(f"{s},synthetic,{s}_R1_001.fastq.gz\n" for s in SAMPLES)

    truth_dir = os.path.join(args.out, "truth")
    os.makedirs(truth_dir, exist_ok=True)

    def write_tsv(name, header, rows):
        with open(os.path.join(truth_dir, name), "w") as f:
            f.write("\t".join(header) + "\n")
            f.writelines("\t".join(str(x) for x in row) + "\n" for row in rows)

    write_tsv("truth.tsv", ["mode", "sample", "table", "key", "expected", "categories"], truth_rows)
    write_tsv("reads_truth.tsv", ["read", "sample", "category", "biotype", "feature", "orientation",
                                  "insert_length", "insert_sequence", "source_contig", "source_start",
                                  "fate_paired_end", "fate_single_end", "reference_tier", "alignments"], key_rows)
    write_tsv("categories.tsv", ["id", "what"] + [f"n_{s}" for s in SAMPLES] +
              ["lost", "mapping", "removed", "rna_class", "gene"],
              [[c["id"], c["what"]] + [n_of(c, i) for i in range(len(SAMPLES))] +
               [c.get("lost", ""), c.get("mapping", ""), c.get("removed", ""),
                c.get("rna_class") or "", c.get("gene") or ""] for c in CATEGORIES])
    print(f"wrote fastq/ and truth/ (reads_truth.tsv, truth.tsv with {len(truth_rows)} rows, categories.tsv) in {args.out}")


if __name__ == "__main__":
    main()
