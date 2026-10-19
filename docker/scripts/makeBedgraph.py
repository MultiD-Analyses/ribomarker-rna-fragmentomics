#!/usr/bin/env python3
# makeBedgraph.py
#
# Per-feature coverage of the parsed alignments, written as a feature-keyed bedgraph.
#
# Reads the PARSED {sample}.bam (one record per assigned read, carrying a 'CS' tag = the
# exact feature name) and writes:  feature  start  end  depth  type.
# Because parseMappings already applied remove_multimap / remove_antisense at the read
# level, this coverage matches the counts table and respects the flags. Each multimapper
# is drawn once, at its assigned locus.
#
# Feature names are the EXACT names from *-featuretypes.csv (e.g. ABALON_multimap,
# A1CF_intron) - no suffix stripping. ABALON and ABALON_multimap are separate tracks.

import argparse
import collections


def read_feature_types(feature_file):
    """{feature: type} and {feature: gene_length} keyed by the EXACT feature name
    (featuretypes columns: feature, feature_type[, gene_length])."""
    type_lookup, length_lookup = {}, {}
    with open(feature_file) as f:
        next(f, None)                                  # header
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split(",")
            type_lookup[parts[0]] = parts[1]
            if len(parts) > 2 and parts[2] != "":
                try:
                    length_lookup[parts[0]] = int(float(parts[2]))
                except ValueError:
                    pass
    return type_lookup, length_lookup


def _keep_set(keep_types):
    types = {t.strip() for t in keep_types.split(",")} if keep_types else {"all"}
    return ("all" in types or not types), types


def _depth_runs(counter):
    """Turn {pos: depth} into contiguous constant-depth runs (start, end, depth)."""
    out = []
    if not counter:
        return out
    xs = sorted(counter)
    s = prev = xs[0]
    d = counter[s]
    for x in xs[1:]:
        if x == prev + 1 and counter[x] == d:
            prev = x
        else:
            out.append((s, prev + 1, d))
            s = prev = x
            d = counter[x]
    out.append((s, prev + 1, d))
    return out


# --------------------------------------------------------------------- coverage
def step_coverage(in_bam, out_bedgraph, type_lookup, keep_types):
    import pysam
    keep_all, keep = _keep_set(keep_types)
    cov = collections.defaultdict(collections.Counter)     # feature -> {pos: depth}
    bam = pysam.AlignmentFile(in_bam, "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.reference_end is None:
            continue
        try:
            feat = r.get_tag("CS")
        except KeyError:
            continue
        typ = type_lookup.get(feat)
        if typ is None or not (keep_all or typ in keep):
            continue
        c = cov[feat]
        for p in range(r.reference_start, r.reference_end):
            c[p] += 1
    bam.close()
    with open(out_bedgraph, "w") as f:
        for feat in cov:
            typ = type_lookup.get(feat, "NA")
            for s, e, d in _depth_runs(cov[feat]):
                f.write(f"{feat}\t{s}\t{e}\t{d}\t{typ}\n")


def main():
    ap = argparse.ArgumentParser(description="Per-feature coverage of the parsed alignments (bedgraph).")
    ap.add_argument("--step", default="coverage", choices=["coverage"])
    ap.add_argument("--input", required=True, help="parsed {sample}.bam")
    ap.add_argument("--output", required=True, help="output bedgraph")
    ap.add_argument("--featuretypes", required=True, help="*-featuretypes.csv")
    ap.add_argument("--keep-types", default="all", help="comma-separated types to keep, or 'all'")
    args = ap.parse_args()

    type_lookup, _ = read_feature_types(args.featuretypes)
    step_coverage(args.input, args.output, type_lookup, args.keep_types)


if __name__ == "__main__":
    main()
