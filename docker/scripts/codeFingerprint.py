#!/usr/bin/env python3
"""
codeFingerprint.py -- one sha256 over the pipeline code that determines a run's outputs (review 2026-10-07).

A git commit alone does not identify the code of a run: a run made with uncommitted changes records a
commit that does not contain them (it happened: the level B/C runs of 2026-10-06 recorded 25e649c, but
their code was committed later as 3f6c30e). The fingerprint is computed from the files themselves, so
equal fingerprints mean the same code. Files (relative to docker/), see CODE_FILES:
    Snakefile, VERSION, env.lock.txt, scripts/*.py, qc/*.yaml, qc/*.css, qc/*.js, qc/logos/*.png, tables/*.tsv (since 0.3.0)
Line endings are normalised (CRLF -> LF) in text files, so Windows and Linux checkouts agree. The
Docker image computes it at build time from the same files (Dockerfile, file /pipeline/CODE_FINGERPRINT).

    python docker/scripts/codeFingerprint.py [DOCKER_DIR]     print the fingerprint (default: this docker/)
    python docker/scripts/codeFingerprint.py --uncommitted    list code files with uncommitted changes
                                                              (git status); exit code 1 if there are any
"""
import hashlib
import subprocess
import sys
from pathlib import Path

CODE_FILES = ['Snakefile', 'VERSION', 'env.lock.txt', 'scripts/*.py', 'qc/*.yaml', 'qc/*.css', 'qc/*.js', 'qc/logos/*.png', 'tables/*.tsv']
DOCKER_DIR = Path(__file__).resolve().parents[1]


def code_files(docker_dir):
    """The files of CODE_FILES, sorted by their path relative to docker_dir; a missing single file is an error."""
    docker_dir = Path(docker_dir)
    found = []
    for pattern in CODE_FILES:
        matches = sorted(docker_dir.glob(pattern))
        if not matches and '*' not in pattern:
            raise FileNotFoundError(f'{pattern} not found in {docker_dir}')
        found += matches
    return sorted(found, key=lambda path: path.relative_to(docker_dir).as_posix())


def fingerprint(docker_dir=DOCKER_DIR):
    docker_dir = Path(docker_dir)
    total = hashlib.sha256()
    for path in code_files(docker_dir):
        data = path.read_bytes()
        if path.suffix != '.png':
            data = data.replace(b'\r\n', b'\n')
        total.update(path.relative_to(docker_dir).as_posix().encode() + b'\0' + hashlib.sha256(data).digest())
    return total.hexdigest()


def uncommitted(docker_dir=DOCKER_DIR):
    """Code files that differ from the last commit or are not committed (git status --porcelain)."""
    docker_dir = Path(docker_dir)
    pathspecs = [f'docker/{pattern}' for pattern in CODE_FILES]
    result = subprocess.run(['git', 'status', '--porcelain', '--', *pathspecs], cwd=docker_dir.parent,
                            capture_output=True, text=True, check=True)
    return [line[3:] for line in result.stdout.splitlines()]


if __name__ == '__main__':
    if sys.argv[1:] == ['--uncommitted']:
        try:
            changed = uncommitted()
        except (OSError, subprocess.CalledProcessError) as error:
            sys.exit(f'cannot check for uncommitted code (git): {error}')
        for name in changed:
            print(f'uncommitted: {name}')
        sys.exit(1 if changed else 0)
    print(fingerprint(sys.argv[1] if len(sys.argv) > 1 else DOCKER_DIR))
