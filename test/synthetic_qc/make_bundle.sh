#!/usr/bin/env bash
###############################################################################
# Step 1 (rarely): create the frozen synthetic test bundle -- see README.md
#
#   bash test/synthetic_qc/make_bundle.sh            # refuses if a bundle exists
#   bash test/synthetic_qc/make_bundle.sh --force    # replace it (new truth!)
#
# Pipeline conda environment active (needs python and bowtie-inspect).
# Writes bundle/references, bundle/fastq, bundle/truth and bundle/MANIFEST.txt;
# bundle/README.md is kept. Downloads go to work/ (not committed).
###############################################################################
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."          # repository root
T=test/synthetic_qc
B=$T/bundle
W=$T/work

if [[ -e "$B/truth/truth.tsv" && "${1:-}" != "--force" ]]; then
  echo "A bundle already exists in $B. It is the frozen truth; use --force to replace it."
  exit 1
fi
command -v bowtie-inspect > /dev/null || { echo "ERROR: bowtie-inspect not found (activate the pipeline conda environment)"; exit 1; }
mkdir -p "$W"

# chromosome 21, same Ensembl release as the production assets
CHR21_URL=https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.21.fa.gz
CHR21=$W/Homo_sapiens.GRCh38.dna_sm.chromosome.21.fa.gz
[[ -s "$CHR21" ]] || curl -fL -o "$CHR21" "$CHR21_URL"

# UniVec sequences from the production bowtie index (no download)
[[ -s "$W/univec.fa" ]] || bowtie-inspect docker/assets/univec/univec_bt1 > "$W/univec.fa"

# generated parts only; bundle/README.md is kept
rm -rf "$B/references" "$B/fastq" "$B/truth" "$B/MANIFEST.txt"
python "$T/make_mini_assets.py" --assets docker/assets --chr21 "$CHR21" --univec "$W/univec.fa" --out "$B/references"
python "$T/make_reads.py" --references "$B/references" --out "$B"

# provenance + checksums (test/levels/run_level_a.sh verifies them before every test)
commit=$(git rev-parse HEAD 2>/dev/null || echo unavailable)
[[ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ]] && commit="$commit (with uncommitted changes)"
{
  echo "created_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "pipeline_commit: $commit"
  echo "python: $(python --version 2>&1)"
  echo "chr21_source: $CHR21_URL"
  echo "chr21_sha256: $(sha256sum "$CHR21" | cut -d' ' -f1)"
  echo "production_assets_manifest_sha256: $(sha256sum docker/assets/MANIFEST_human.json | cut -d' ' -f1)"
  echo "univec_source: bowtie-inspect docker/assets/univec/univec_bt1 (first 3 records)"
  echo "selection_and_seed: set in make_mini_assets.py (window, records) and make_reads.py (SEED, categories)"
  echo
  echo "files (sha256):"
  (cd "$B" && find references fastq truth -type f | sort | xargs sha256sum)
} > "$B/MANIFEST.txt"

echo "Bundle written to $B. Review it, then commit it."
