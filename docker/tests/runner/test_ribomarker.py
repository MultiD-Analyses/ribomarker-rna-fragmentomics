"""The host-side starter `ribomarker` (bash): the docker run command it writes, with a
stand-in for docker that prints its arguments. Run from the repository root: python -m pytest docker/tests/runner
(ribomarker_windows.ps1 does the same in PowerShell; it is tested by hand on Windows.)"""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
SCRIPT = REPO / 'ribomarker'
# the image the starter uses by default (its name and version change only at a release)
IMAGE = re.search(r'^IMAGE=\$\{RIBOMARKER_IMAGE:-(\S+)\}$', SCRIPT.read_text(), re.M).group(1)


@pytest.fixture
def docker(tmp_path):
    """A stand-in docker: 'info' succeeds (or fails if DOCKER_DOWN is set); otherwise one argument per line."""
    fake = tmp_path / 'fake_docker'
    fake.write_text('#!/bin/bash\nif [[ $1 == info ]]; then [[ -z ${DOCKER_DOWN:-} ]]; exit; fi\n'
                    'printf "%s\\n" "$@"\n')
    fake.chmod(0o755)
    return fake


def ribomarker(docker, *args, cwd, **env):
    result = subprocess.run(['bash', str(SCRIPT), *map(str, args)], cwd=cwd, capture_output=True, text=True,
                            env=dict(os.environ, RIBOMARKER_DOCKER=str(docker), **env), stdin=subprocess.DEVNULL)
    return result.returncode, result.stdout.splitlines(), result.stderr


def mounts(lines):
    return [lines[i + 1] for i, a in enumerate(lines) if a == '-v']


def sheet(folder):
    """A folder with a sample sheet (the starter only checks that the sheet exists)."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'samples.csv').write_text('sample,fastq_file\nA,a.fastq.gz\n')
    return folder / 'samples.csv'


def test_run_mounts_the_sheet_folder_and_defaults_to_8_threads(tmp_path, docker):
    sheet(tmp_path / 'my fastq')                         # a space in the path
    code, lines, _ = ribomarker(docker, 'my fastq/samples.csv', 'out/results', cwd=tmp_path)
    assert code == 0
    assert lines[:6] == ['run', '--rm', '--init', '--network', 'none', '--user']
    assert '-it' not in lines                            # no terminal here
    assert mounts(lines) == [f'{tmp_path}/my fastq:/data/input:ro', f'{tmp_path}/out/results:/data/output']
    assert (tmp_path / 'out' / 'results').is_dir()     # made by the script, so it belongs to the user
    i = lines.index(IMAGE)
    assert lines[i:] == [IMAGE, 'run', '/data/input/samples.csv', '/data/output', '--threads', '8']


def test_options_are_passed_on_and_files_mounted(tmp_path, docker):
    sheet(tmp_path / 'in')
    (tmp_path / 'refs').mkdir()
    (tmp_path / 'settings.yaml').write_text('multimap: 50\n')
    code, lines, _ = ribomarker(docker, 'in/samples.csv', 'out', '--threads', '4', '--species', 'human',
                                '--run-id', 'r1', '--config', 'settings.yaml', '--references', 'refs', '--unlock',
                                cwd=tmp_path)
    assert code == 0
    assert f'{tmp_path}/settings.yaml:/data/config/settings.yaml:ro' in mounts(lines)
    assert f'{tmp_path}/refs:/data/references:ro' in mounts(lines)
    tail = lines[lines.index('run', 1):]
    assert tail[:5] == ['run', '/data/input/samples.csv', '/data/output', '--threads', '4']
    for pair in (['--species', 'human'], ['--run-id', 'r1'], ['--config', '/data/config/settings.yaml'],
                 ['--references', '/data/references']):
        assert any(tail[k:k + 2] == pair for k in range(len(tail))), pair
    assert '--unlock' in tail and tail.count('--threads') == 1


def test_demo_and_help(tmp_path, docker):
    code, lines, _ = ribomarker(docker, '--demo', 'demo', cwd=tmp_path)
    assert code == 0 and mounts(lines) == [f'{tmp_path}/demo:/data/output']
    assert lines[-5:] == ['run', '--demo', '/data/output', '--threads', '8']
    code, lines, _ = ribomarker(docker, '--help', cwd=tmp_path)
    assert code == 0 and lines[0].startswith('ribomarker SHEET OUTPUT')
    assert lines[-5:] == ['run', '--rm', IMAGE, 'run', '--help']


@pytest.mark.parametrize('args, message', [
    (['missing.csv', 'out'], 'sample sheet not found'),
    (['in', 'out'], 'sample sheet not found'),                 # a folder is not a sheet
    (['in/samples.csv', 'in/results'], 'must not be the folder of the sample sheet or lie inside it'),
    (['in/samples.csv', 'in'], 'must not be the folder of the sample sheet or lie inside it'),
    (['in/samples.csv'], 'usage: ribomarker SHEET OUTPUT'),
    (['in/samples.csv', 'out', '--threds', '4'], 'unknown option --threds'),
    (['in/samples.csv', 'out', '--metadata', 'm.csv'], 'unknown option --metadata'),
    (['in/samples.csv', 'out', '--config', 'none.yaml'], 'file not found'),
    (['in/samples.csv', 'out', '--threads'], 'needs a value'),
])
def test_mistakes_are_refused_before_docker_starts(tmp_path, docker, args, message):
    sheet(tmp_path / 'in')
    code, lines, err = ribomarker(docker, *args, cwd=tmp_path)
    assert code == 2 and message in err and lines == []
    assert not (tmp_path / 'in' / 'results').exists()     # nothing left behind in the input folder


def test_docker_not_running(tmp_path, docker):
    sheet(tmp_path / 'in')
    code, lines, err = ribomarker(docker, 'in/samples.csv', 'out', cwd=tmp_path, DOCKER_DOWN='1')
    assert code == 2 and 'Docker is not running' in err


def test_other_image(tmp_path, docker):
    sheet(tmp_path / 'in')
    code, lines, _ = ribomarker(docker, 'in/samples.csv', 'out', cwd=tmp_path, RIBOMARKER_IMAGE='ribomarker:test')
    assert code == 0 and 'ribomarker:test' in lines


def test_both_starters_use_the_same_image():
    ps1 = (REPO / 'ribomarker_windows.ps1').read_text()
    assert f"else {{ '{IMAGE}' }}" in ps1
    assert re.fullmatch(r'[a-z0-9._-]+/[a-z0-9._-]+:\d+\.\d+\.\d+', IMAGE)    # Docker Hub: NAMESPACE/NAME:X.Y.Z
