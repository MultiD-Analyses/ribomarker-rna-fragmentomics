#!/usr/bin/env bash
###############################################################################
# Demo data for `run --demo`: the frozen level A bundle (paired-end reads and
# their truth) and its small references with bowtie indexes and a MANIFEST, in docker/demo/. Not in
# git; the Dockerfile copies it into the image, and a local `run --demo` uses it from here.
# Run from anywhere in the repository, pipeline conda environment active, before `docker build`:
#
#   bash docker/make_demo.sh
###############################################################################
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."          # repository root
BUNDLE=test/synthetic_qc/bundle
DEMO=docker/demo
LOG=test/synthetic_qc/work/make_demo.log

# the bundle must be exactly the frozen one (the same check as the level A runner)
(cd "$BUNDLE" && sed '1,/^files (sha256):$/d' MANIFEST.txt | sha256sum --check --quiet) \
  || { echo "ERROR: bundle files differ from $BUNDLE/MANIFEST.txt"; exit 1; }

rm -rf "$DEMO"
mkdir -p "$DEMO/fastq" "$(dirname "$LOG")"
cp "$BUNDLE"/fastq/*_R[12]_001.fastq.gz "$DEMO/fastq/"
cp "$BUNDLE/truth/truth.tsv" "$DEMO/truth.tsv"
{ echo "sample,fastq_file,fastq_file2"                 # the sample sheet that run --demo reads
  for r1 in "$DEMO"/fastq/*_R1_001.fastq.gz; do
    s=$(basename "$r1" _R1_001.fastq.gz); echo "$s,${s}_R1_001.fastq.gz,${s}_R2_001.fastq.gz"
  done; } > "$DEMO/fastq/samples.csv"
python test/synthetic_qc/build_indexes.py --references "$BUNDLE/references" --out "$DEMO/references" > "$LOG" 2>&1 \
  || { echo "ERROR: index building failed, see $LOG"; exit 1; }
python buildDB/update_manifest.py "$DEMO/references" --species human >> "$LOG" 2>&1 \
  || { echo "ERROR: MANIFEST not written, see $LOG"; exit 1; }
echo "demo data in $DEMO ($(du -sh "$DEMO" | cut -f1)); try it: python docker/scripts/main.py run --demo results/demo"
