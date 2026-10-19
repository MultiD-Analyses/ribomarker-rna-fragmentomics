"""Shared helpers for the QC tests (run from the repository root: python -m pytest docker/tests/qc)."""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[3]
FIXTURES = Path(__file__).resolve().parent / 'fixtures'
COLLECTOR = REPO / 'docker' / 'scripts' / 'qcCollect.py'
DICTIONARY = yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8'))

sys.path.insert(0, str(REPO / 'docker' / 'scripts'))      # lets the tests import qcCollect


def prepare_run(tmp_path, mode):
    """Copy the level A fixture of one mode to tmp_path and write a config for it.
    Returns (run folder, config file)."""
    run = tmp_path / f'levelA_{mode}'
    shutil.copytree(FIXTURES / f'levelA_{mode}', run)
    shutil.copy(FIXTURES / f'samples_{mode}.csv', tmp_path / f'samples_{mode}.csv')
    config = {
        # synthetic references: BED for check C8, but no MANIFEST (one expected warning)
        'dbdir': str(REPO / 'test' / 'synthetic_qc' / 'bundle' / 'references') + '/',
        'input_dir': str(tmp_path) + '/',
        'output_prefix': str(run),
        'samples': f'samples_{mode}.csv',
        'species': 'human',
        'paired_end': mode == 'pe',
        'remove_multimap': True,
        'remove_antisense': True,
    }
    config_file = tmp_path / f'config_{mode}.yaml'
    config_file.write_text(yaml.safe_dump(config))
    return run, config_file


def run_collector(run, config_file):
    """Run qcCollect.py as the pipeline would; returns the parsed qc_metrics.json."""
    result = subprocess.run([sys.executable, str(COLLECTOR), '--output', str(run), '--config', str(config_file)],
                            cwd=REPO, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return json.loads((run / 'qc' / 'qc_metrics.json').read_text())


def values(report):
    """{(sample, metric id): value} for all scalar metrics and distributions of a report."""
    found = {}
    for sample, data in report['samples'].items():
        for record in data['metrics']:
            found[(sample, record['id'])] = record['value']
        for record in data['distributions']:
            found[(sample, record['id'])] = (record['x'], record['y'])
    return found


def template_pattern(template):
    """Regular expression that matches the concrete ids of a dictionary id with placeholders."""
    return re.compile('^' + re.sub(r'\\\{\w+\\\}', r'[^.]+', re.escape(template)) + '$')


@pytest.fixture(scope='session')
def baselines(tmp_path_factory):
    """{mode: report} of the collector run on the unchanged level A fixtures."""
    reports = {}
    for mode in ('pe', 'se'):
        run, config_file = prepare_run(tmp_path_factory.mktemp(f'baseline_{mode}'), mode)
        reports[mode] = run_collector(run, config_file)
    return reports


@pytest.fixture(params=['pe', 'se'])
def collected(request, baselines):
    """(mode, report) of the undamaged run; tests using it run once per mode."""
    return request.param, baselines[request.param]
