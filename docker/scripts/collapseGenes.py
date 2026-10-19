#!/usr/bin/env python3
"""
collapseGenes.py -- collapse the per-feature counts matrix to a GENE-LEVEL,
mature-exonic counts matrix (integer counts, suitable for count-based differential-expression tools).

Rules (agreed 2026-07):
  DROP  (not mature)      : pre-mRNA (the *_intron features), pre-tRNA, miRNA-pre
  EXCLUDE (not genes)     : RiboSpike-* (spike-in controls), unannotated, unmapped
  GENE KEY
    tRNA / Mt_tRNA        : strip a trailing _5p/_3p/_internal/_whole (these are
                           parseMappings positional FRAGMENT tags) and sum the
                           fragments. Copies/isodecoders (…-1 vs …-14) stay separate.
    everything else kept  : the feature name AS-IS. So miRNA arms (…_5p/…_3p are
                           MirGeneDB arm names, NOT fragment tags), miRNA-star,
                           rRNA copies, and pipe-joined shared names all stay separate.
  Sum counts per gene key; one gene_length per gene via groupby(gene).max().

Assumes the input already reflects remove_multimap/remove_antisense (no _multimap/
_antisense features); it errors if any are present so the choice stays explicit.

Inputs
  --readcounts    {out}-readcounts.csv    (feature x sample)
  --featuretypes  {out}-featuretypes.csv  (feature, feature_type[, gene_length])
Outputs (next to --output prefix)
  {out}-genecounts.csv    (gene x sample)
  {out}-genelengths.csv   (gene, length, biotype)
"""

import argparse
import re
import sys
import pandas as pd

DROP_BIOTYPES = {"pre-mRNA", "pre-tRNA", "miRNA-pre"}   # precursor / intronic -> not mature
EXCLUDE_BIOTYPES = {"unannotated", "unmapped"}          # + anything starting "RiboSpike"
FRAGMENT_BIOTYPES = {"tRNA", "Mt_tRNA"}                  # only these carry positional tags
FRAGMENT_SUFFIX = re.compile(r"_(5p|3p|internal|whole)$")


def is_excluded(biotype):
    return biotype in EXCLUDE_BIOTYPES or str(biotype).startswith("RiboSpike")


def gene_key(feature, biotype):
    """tRNA/Mt_tRNA: drop the positional fragment tag (sum fragments). Else: as-is."""
    if biotype in FRAGMENT_BIOTYPES:
        return FRAGMENT_SUFFIX.sub("", feature)
    return feature


def main():
    ap = argparse.ArgumentParser(description="Collapse per-feature counts to gene level (design Z, 5a).")
    ap.add_argument("--readcounts", required=True, help="{out}-readcounts.csv (feature x sample)")
    ap.add_argument("--featuretypes", required=True, help="{out}-featuretypes.csv (feature, feature_type[, gene_length])")
    ap.add_argument("--output", required=True, help="output prefix; writes {output}-genecounts.csv / -genelengths.csv")
    args = ap.parse_args()

    counts = pd.read_csv(args.readcounts, index_col=0)              # index = feature, cols = samples
    ft = pd.read_csv(args.featuretypes).drop_duplicates("feature").set_index("feature")

    # every counted feature must have a known biotype
    biotype = ft["feature_type"].reindex(counts.index)
    missing = counts.index[biotype.isna()]
    if len(missing):
        sys.exit(f"ERROR: {len(missing)} features in readcounts have no biotype in featuretypes "
                 f"(e.g. {list(missing)[:5]}). Refusing to guess.")

    # guard: input should already be remove-at-source (unique + sense)
    leaked = [f for f in counts.index if "_multimap" in f or "_antisense" in f]
    if leaked:
        sys.exit(f"ERROR: {len(leaked)} _multimap/_antisense features present (e.g. {leaked[:5]}). "
                 f"Gene-collapse expects remove_multimap/remove_antisense output. Aborting.")

    # classify each feature
    dropped = biotype.isin(DROP_BIOTYPES)
    excluded = biotype.map(is_excluded)
    keep = ~dropped & ~excluded
    if keep.sum() == 0:
        sys.exit("ERROR: no features left after drop/exclude — check inputs.")

    kept = counts[keep].copy()
    kept_bt = biotype[keep]
    kept["gene"] = [gene_key(f, b) for f, b in zip(kept.index, kept_bt)]

    sample_cols = [c for c in counts.columns]
    gene_counts = kept.groupby("gene")[sample_cols].sum()

    # one length + biotype per gene (length via max, as in 5.9)
    glen = ft["gene_length"].reindex(kept.index) if "gene_length" in ft.columns else pd.Series(index=kept.index, dtype=float)
    meta = pd.DataFrame({"gene": kept["gene"].values, "length": glen.values, "biotype": kept_bt.values})
    gene_lengths = meta.groupby("gene").agg(length=("length", "max"), biotype=("biotype", "first"))
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
    print(f"  dropped : {int(dropped.sum())} features (pre-mRNA/pre-tRNA/miRNA-pre), "
          f"{counts[dropped][sample_cols].to_numpy().sum():.0f} counts")
    print(f"  excluded: {int(excluded.sum())} features (RiboSpike/unannotated/unmapped), "
          f"{counts[excluded][sample_cols].to_numpy().sum():.0f} counts")
    print(f"  conservation OK: {in_sum:.0f} == {out_sum:.0f}")
    print(f"  wrote {args.output}-genecounts.csv, {args.output}-genelengths.csv")
    print("-" * 55)


if __name__ == "__main__":
    main()
