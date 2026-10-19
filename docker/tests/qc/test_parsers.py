"""Unit tests of the collector's small parsers and helpers."""
import shutil

import pytest

import codeFingerprint
import qcCollect
from conftest import REPO

BOWTIE2_LOG = """[WARNING] Failed to launch x86-64-v3 version, staying with default
300067 reads; of these:
  300067 (100.00%) were unpaired; of these:
    281908 (93.95%) aligned 0 times
    11996 (4.00%) aligned exactly 1 time
    6163 (2.05%) aligned >1 times
6.05% overall alignment rate
"""

BOWTIE1_LOG = """# reads processed: 1601301
# reads with at least one reported alignment: 17079 (1.07%)
# reads that failed to align: 1584222 (98.93%)
Reported 17079 alignments to 1 output stream(s)
"""


def test_bowtie2_log_ignores_warning_lines(tmp_path):
    path = tmp_path / 'tier.log'
    path.write_text(BOWTIE2_LOG)
    assert qcCollect.read_tier_log(path) == {'reads_in': 300067, 'aligned': 18159, 'unaligned': 281908,
                                             'multi': 6163}


def test_bowtie1_log(tmp_path):
    path = tmp_path / 'tier.log'
    path.write_text(BOWTIE1_LOG)
    assert qcCollect.read_tier_log(path) == {'reads_in': 1601301, 'aligned': 17079, 'unaligned': 1584222,
                                             'multi': None}


def test_bowtie2_log_without_reads(tmp_path):
    """A tier that receives no reads (all aligned earlier): Bowtie 2 prints only two lines
    (verbatim from level A sample S3, 2026-10-05)."""
    path = tmp_path / 'tier.log'
    path.write_text('[WARNING] Failed to launch x86-64-v3 version, staying with default\n'
                    '0 reads\n0.00% overall alignment rate\n')
    assert qcCollect.read_tier_log(path) == {'reads_in': 0, 'aligned': 0, 'unaligned': 0, 'multi': 0}


def test_tier_log_with_inconsistent_totals_is_rejected(tmp_path):
    path = tmp_path / 'tier.log'
    path.write_text(BOWTIE2_LOG.replace('281908', '281907'))
    with pytest.raises(ValueError):
        qcCollect.read_tier_log(path)


@pytest.mark.parametrize('text, expected', [
    ('/mnt/c/Users/someone/project/docker/assets/univec/univec_bt1 T', 'univec_bt1 T'),
    ('/home/someone/env/bin/bowtie2-build-s version 2.5.5', 'bowtie2-build-s version 2.5.5'),
    ('/tmp/run/levelA_pe/', 'levelA_pe'),
    ('https://ftp.ensembl.org/pub/release-116/gtf/x.gtf.gz', 'https://ftp.ensembl.org/pub/release-116/gtf/x.gtf.gz'),
    ('results/levelA_pe/trimmed', 'results/levelA_pe/trimmed'),
    ('ratio 5/3 and/or', 'ratio 5/3 and/or'),
])
def test_strip_paths(text, expected):
    assert qcCollect.strip_paths(text) == expected


def test_strip_paths_recurses():
    assert qcCollect.strip_paths({'a': ['/x/y/z'], 'b': 3}) == {'a': ['z'], 'b': 3}


@pytest.mark.parametrize('counts, expected', [([0, 3, 0], 1.0), ([1, 1, 0], 0.5), ([1, 0, 1], 1.0), ([0, 0, 5], 2.0)])
def test_histogram_median(counts, expected):
    assert qcCollect.histogram_median([0, 1, 2], counts) == expected


@pytest.mark.parametrize('text', ['1.5', 'abc', None])
def test_to_int_rejects_non_whole_numbers(text):
    with pytest.raises((ValueError, TypeError)):
        qcCollect.to_int(text)


def test_short_csv_row_is_rejected(tmp_path):
    path = tmp_path / 'table.csv'
    path.write_text('biotype,S1,S2\nrRNA,1,2\nmiRNA,3\n')
    with pytest.raises(ValueError):
        qcCollect.read_rows(path)


# ----------------------------------------------------------------------------- code fingerprint

def copy_code(tmp_path):
    """A copy of the code files of docker/ (codeFingerprint.CODE_FILES) in tmp_path/docker."""
    target = tmp_path / 'docker'
    for path in codeFingerprint.code_files(REPO / 'docker'):
        relative = path.relative_to(REPO / 'docker')
        (target / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(path, target / relative)
    return target


def test_code_fingerprint_follows_the_code(tmp_path):
    code = copy_code(tmp_path)
    original = codeFingerprint.fingerprint(code)
    assert original == codeFingerprint.fingerprint(REPO / 'docker')   # tests, docs and assets do not count
    script = code / 'scripts' / 'qcCollect.py'
    script.write_bytes(script.read_bytes().replace(b'\n', b'\r\n'))
    assert codeFingerprint.fingerprint(code) == original                # nor do line endings
    script.write_bytes(script.read_bytes() + b'# changed\r\n')
    assert codeFingerprint.fingerprint(code) != original                # a code change does


def test_collector_records_the_code_fingerprint(collected, tmp_path, monkeypatch):
    _, report = collected
    version = report['provenance']['prov.pipeline_version']
    assert version['code_fingerprint'] == codeFingerprint.fingerprint(REPO / 'docker')
    # Docker image: no Python sources next to VERSION, the fingerprint written at build time is read
    (tmp_path / 'VERSION').write_text('0.1.0\n')
    (tmp_path / 'CODE_FINGERPRINT').write_text('abc123\n')
    monkeypatch.setattr(qcCollect, 'VERSION_FILE', tmp_path / 'VERSION')
    assert qcCollect.read_version()['code_fingerprint'] == 'abc123'
