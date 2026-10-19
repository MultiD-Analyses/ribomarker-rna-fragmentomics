#!/bin/bash
# Container entrypoint: runs the pipeline's run command inside the ribomarker conda environment.
#
# Usage (from docker run; the starter scripts write this command):
#   docker run ... IMAGE run SHEET OUTPUT [options]
#   docker run ... IMAGE run --demo OUTPUT
#   docker run ... IMAGE run --help

set -euo pipefail

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate ribomarker

# The run command: docker run ... IMAGE run SHEET OUTPUT [options]. It works only in
# OUTPUT (Snakemake's working folder), so the image can be read-only and any --user works. exec, so
# that a Ctrl-C or `docker stop` reaches it; no --unlock here (a lock is only removed with run --unlock).
if [ "${1:-}" = "run" ]; then
  shift
  if [ ! -w "${HOME:-/}" ]; then
    export HOME=/tmp       # a --user without a home folder: Snakemake and MultiQC keep caches under $HOME
  fi
  exec python /pipeline/scripts/main.py run "$@"
fi

# Without arguments, or with --help: the run command's help. Anything else is refused.
case "${1:-}" in
  ""|-h|--help) exec python /pipeline/scripts/main.py run --help ;;
esac
echo "ERROR: unknown command '$1'. This image runs: run SHEET OUTPUT [options]   (help: run --help)" >&2
exit 2
