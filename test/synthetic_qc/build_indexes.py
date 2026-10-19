#!/usr/bin/env python3
"""
build_indexes.py -- turn the bundle's reference FASTA files into pipeline assets.

Builds the bowtie / bowtie2 indexes with buildDB's own functions
(buildDB/builder_core.py, imported unchanged) and writes human.bed and
human.dbpaths.csv (index paths relative to --out) into --out.
Run at test time, so indexes always match the local bowtie version.

Usage (from the repository root, pipeline conda environment active):
  python test/synthetic_qc/build_indexes.py \
      --references test/synthetic_qc/bundle/references --out test/synthetic_qc/work/assets_a
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "buildDB"))
import builder_core  # noqa: E402  -- buildDB code, used as-is


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--references", default="test/synthetic_qc/bundle/references")
    ap.add_argument("--out", default="test/synthetic_qc/work/assets_a")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    dbrows = []
    with open(os.path.join(args.references, "tiers.csv")) as f:
        for tier in csv.DictReader(f):                       # rows are in mapping order
            fasta = os.path.join(args.references, tier["fasta"])
            index_dir = os.path.join(args.out, os.path.dirname(tier["fasta"]))
            os.makedirs(index_dir, exist_ok=True)
            index = os.path.join(index_dir, tier["index"])
            log = index + ".build.log"
            if tier["map_cmd"] == "bowtie":
                builder_core._bowtie1_build(fasta, index, args.threads, log)
            else:
                builder_core._bowtie2_build(fasta, index, args.threads, log)
            dbrows.append([tier["map_cmd"], builder_core._rel(args.out, index), "T"])

    bed = os.path.join(args.references, "human.bed")
    dbpaths_csv, _ = builder_core.write_outputs(dbrows, [bed], args.out, "human")
    print(f"built {len(dbrows)} indexes; wrote {dbpaths_csv}")


if __name__ == "__main__":
    main()
