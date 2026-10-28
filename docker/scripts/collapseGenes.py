#!/usr/bin/env python3
"""
collapseGenes.py -- collapse the per-feature counts matrix to a GENE-LEVEL,
mature-exonic counts matrix (integer counts, suitable for count-based differential-expression tools).

Rules (agreed 2026-07; since 0.3.0 read from the class table, column gene_table, see refTables.py):
  DROP  (not mature)      : gene_table = drop: pre-mRNA (the *_intron features), pre-tRNA, pre-miRNA
  EXCLUDE (not genes)     : gene_table = exclude: RiboSpike-* (spike-in controls), genome_unannotated, unmapped
  GENE KEY
    positions = yes       : tRNA / Mt_tRNA / snRNA / snoRNA: strip a trailing _5p/_3p/_internal/_whole (these are
                           parseMappings positional FRAGMENT tags) and sum the
                           fragments. Copies/isodecoders (…-1 vs …-14) stay separate.
    everything else kept  : the feature name AS-IS. So miRNA arms (…_5p/…_3p are
                           MirGeneDB arm names, NOT fragment tags), miRNA-star,
                           rRNA copies, and pipe-joined shared names all stay separate.
  Sum counts per gene key; one length_nt per gene via groupby(gene).max().

Assumes the input already reflects remove_multimap/remove_antisense (no _multimap/
_antisense features); it errors if any are present so the choice stays explicit.

Inputs
  --readcounts    {out}-readcounts.csv    (feature x sample)
  --featureclasses {out}-featureclasses.csv (feature, rna_class[, length_nt])
  --rna-classes   class table (default: docker/tables/default.rna_classes.tsv)
Outputs (next to --output prefix)
  {out}-genecounts.csv    (gene x sample)
  {out}-genelengths.csv   (gene, length_nt, rna_class[, the gene-symbol columns of featureclasses])
Since 0.3.0 the RNA classes are class ids of the class table (refTables.py); the columns were feature_type,
gene_length (featuretypes.csv) and length, biotype (genelengths.csv) before.
"""

import argparse
import re
import sys
import pandas as pd

import refTables

FRAGMENT_SUFFIX = re.compile(r"_(5p|3p|internal|whole)$")


def gene_key(feature, rna_class, position_classes):
    """Classes with positions (tRNA, Mt_tRNA, snRNA, snoRNA): drop the positional fragment tag (sum fragments).
    Else: as-is."""
    if rna_class in position_classes:
        return FRAGMENT_SUFFIX.sub("", feature)
    return feature


def main():
    ap = argparse.ArgumentParser(description="Collapse per-feature counts to gene level (design Z, 5a).")
    ap.add_argument("--readcounts", required=True, help="{out}-readcounts.csv (feature x sample)")
    ap.add_argument("--featureclasses", required=True,
                    help="{out}-featureclasses.csv (feature, rna_class[, length_nt])")
    ap.add_argument("--rna-classes", default=None, help="RNA class table (default: docker/tables/default.rna_classes.tsv)")
    ap.add_argument("--output", required=True, help="output prefix; writes {output}-genecounts.csv / -genelengths.csv")
    args = ap.parse_args()
    classes = refTables.ClassTable(args.rna_classes or refTables.default_table("default"), "default")

    counts = pd.read_csv(args.readcounts, index_col=0)              # index = feature, cols = samples
    ft = pd.read_csv(args.featureclasses).drop_duplicates("feature").set_index("feature")

    # every counted feature must have a known biotype
    biotype = ft["rna_class"].reindex(counts.index)
    missing = counts.index[biotype.isna()]
    if len(missing):
        sys.exit(f"ERROR: {len(missing)} features in readcounts have no RNA class in featureclasses "
                 f"(e.g. {list(missing)[:5]}). Refusing to guess.")

    # guard: input should already be remove-at-source (unique + sense)
    leaked = [f for f in counts.index if "_multimap" in f or "_antisense" in f]
    if leaked:
        sys.exit(f"ERROR: {len(leaked)} _multimap/_antisense features present (e.g. {leaked[:5]}). "
                 f"Gene-collapse expects remove_multimap/remove_antisense output. Aborting.")

    # classify each feature
    rule = biotype.map(classes.gene_table)           # yes / drop / exclude
    dropped = rule == "drop"
    excluded = rule == "exclude"
    keep = ~dropped & ~excluded
    if keep.sum() == 0:
        sys.exit("ERROR: no features left after drop/exclude — check inputs.")

    kept = counts[keep].copy()
    kept_bt = biotype[keep]
    position_classes = classes.position_classes()
    kept["gene"] = [gene_key(f, b, position_classes) for f, b in zip(kept.index, kept_bt)]

    sample_cols = [c for c in counts.columns]
    gene_counts = kept.groupby("gene")[sample_cols].sum()

    # one length + biotype per gene (length via max, as in 5.9)
    glen = ft["length_nt"].reindex(kept.index) if "length_nt" in ft.columns else pd.Series(index=kept.index, dtype=float)
    meta = pd.DataFrame({"gene": kept["gene"].values, "length_nt": glen.values, "rna_class": kept_bt.values})
    gene_lengths = meta.groupby("gene").agg(length_nt=("length_nt", "max"), rna_class=("rna_class", "first"))
    # gene-symbol columns of featureclasses (release 0.3.0, 3f): the first non-empty value of the gene's features
    for column in [c for c in ft.columns if c not in ("rna_class", "length_nt")]:
        values = pd.Series(ft[column].reindex(kept.index).fillna("").astype(str).values, index=kept["gene"].values)
        gene_lengths[column] = values[values != ""].groupby(level=0).first().reindex(gene_lengths.index).fillna("")
    gene_lengths = gene_lengths.reindex(gene_counts.index)

    # conservation: the collapse is a pure relabel+sum of the KEPT features
    in_sum = kept[sample_cols].to_numpy().sum()
    out_sum = gene_counts.to_numpy().sum()
    if in_sum != out_sum:
        sys.exit(f"ERROR: conservation failed ({in_sum} kept != {out_sum} collapsed).")

    gene_counts.to_csv(f"{args.output}-genecounts.csv")
    gene_lengths.to_csv(f"{args.output}-genelengths.csv")

    # accounting summary
    tot = counts[sample_cols].to_numpy().sum()
    print("-" * 55)
    print(f"gene-collapse: {len(counts)} features -> {len(gene_counts)} genes")
    print(f"  kept    : {int(keep.sum())} features, {in_sum:.0f} counts ({100*in_sum/tot:.1f}%)")
    print(f"  dropped : {int(dropped.sum())} features (gene_table = drop: pre-mRNA/pre-tRNA/pre-miRNA), "
          f"{counts[dropped][sample_cols].to_numpy().sum():.0f} counts")
    print(f"  excluded: {int(excluded.sum())} features (RiboSpike/unannotated/unmapped), "
          f"{counts[excluded][sample_cols].to_numpy().sum():.0f} counts")
    print(f"  conservation OK: {in_sum:.0f} == {out_sum:.0f}")
    print(f"  wrote {args.output}-genecounts.csv, {args.output}-genelengths.csv")
    print("-" * 55)


if __name__ == "__main__":
    main()
