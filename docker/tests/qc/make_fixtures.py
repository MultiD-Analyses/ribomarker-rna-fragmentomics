#!/usr/bin/env python3
"""
make_fixtures.py -- copy the text outputs of the level A (synthetic) runs that
qcCollect.py reads into docker/tests/qc/fixtures/, for the QC tests.

Run once after a level A run (bash test/levels/run_level_a.sh), from the repository root:
    python docker/tests/qc/make_fixtures.py
Absolute paths in the copied files (map.log) are replaced by /repo. Existing fixture
files are only overwritten with --force; files that are no longer produced are listed,
not deleted.
"""
import argparse
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
FIXTURES = REPO / 'docker' / 'tests' / 'qc' / 'fixtures'
SHEETS = REPO / 'test' / 'synthetic_qc' / 'bundle' / 'fastq'

# files read by qcCollect.py, relative to the run folder; {run} = run folder name
FILES = [
    'samples.trimmed.csv',
    'logs/map.log',
    'trimmed/logs/*.trim.cutadapt.json',
    'trimmed/logs/*.dimer.cutadapt.json',
    'trimmed/logs/*.fastp.json',
    'alignments/logs/*.log',
    'metrics/*.log',
    'results/{run}-mappingstats.csv',
    'results/{run}-classcounts.csv',
    'results/{run}-classlengths.csv',
    'results/{run}-classdistinctseqs.csv',
    'results/{run}-genecounts.csv',
    'results/{run}-genelengths.csv',
    'ribospike/results/{run}-readorigin.csv',
    'ribospike/results/{run}-ribospike.csv',
    'ribospike/results/{run}-spikelengths.csv',
    'ribospike/results/{run}-ribospikeends.csv',
    'ribospike/results/{run}-ribospikeends_exact.csv',
]
ABSOLUTE_REPO_PATH = re.compile(r'/\S*/' + re.escape(REPO.name))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--force', action='store_true', help='overwrite existing fixture files')
    args = parser.parse_args()
    if FIXTURES.exists() and not args.force:
        sys.exit(f'{FIXTURES} exists; use --force to overwrite its files')
    before = {path for path in FIXTURES.rglob('*') if path.is_file() and '__pycache__' not in path.parts}
    written = set()

    for mode in ('pe', 'se'):
        run = f'levelA_{mode}'
        source = REPO / 'results' / run
        target = FIXTURES / run
        n_files = 0
        for pattern in FILES:
            for path in sorted(source.glob(pattern.format(run=run))):
                copy = target / path.relative_to(source)
                copy.parent.mkdir(parents=True, exist_ok=True)
                copy.write_text(ABSOLUTE_REPO_PATH.sub('/repo', path.read_text()))
                written.add(copy)
                n_files += 1
        shutil.copy(SHEETS / f'samples_{mode}.csv', FIXTURES / f'samples_{mode}.csv')
        written.add(FIXTURES / f'samples_{mode}.csv')
        print(f'{run}: {n_files} files copied to {target.relative_to(REPO)}')
    for stale in sorted(before - written):
        print(f'no longer produced (remove by hand): {stale.relative_to(REPO)}')


if __name__ == '__main__':
    main()
