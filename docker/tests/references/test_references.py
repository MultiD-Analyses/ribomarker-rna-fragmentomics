"""References folder: relative dbpaths file and MANIFEST checksums.
Run from the repository root: python -m pytest docker/tests/references"""
import io
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'docker' / 'scripts'))
sys.path.insert(0, str(REPO / 'buildDB'))
import dbPaths          # noqa: E402
import update_manifest  # noqa: E402

UPDATE = REPO / 'buildDB' / 'update_manifest.py'


# ----------------------------------------------------------------------------- dbPaths.find_dbpaths
def test_new_dbpaths_file_is_preferred(tmp_path):
    (tmp_path / 'human.dbpaths.csv').write_text('bowtie2,genome/human-genome,T\n')
    (tmp_path / 'human.local.dbpaths.csv').write_text('bowtie2,/x/genome/human-genome,T\n')
    for in_container in (False, True):
        assert dbPaths.find_dbpaths(str(tmp_path), 'human', in_container) == str(tmp_path / 'human.dbpaths.csv')


def test_older_file_only_outside_a_container(tmp_path, capsys):
    (tmp_path / 'human.local.dbpaths.csv').write_text('bowtie2,/x/genome/human-genome,T\n')
    (tmp_path / 'human.docker.dbpaths.csv').write_text('bowtie2,/pipeline/assets/genome/human-genome,T\n')
    assert dbPaths.find_dbpaths(str(tmp_path), 'human', False) == str(tmp_path / 'human.local.dbpaths.csv')
    assert 'WARNING' in capsys.readouterr().out
    with pytest.raises(FileNotFoundError, match='update_manifest.py'):
        dbPaths.find_dbpaths(str(tmp_path), 'human', True)


def test_container_is_recognised_from_the_environment(tmp_path, monkeypatch):
    (tmp_path / 'human.local.dbpaths.csv').write_text('bowtie2,/x/genome/human-genome,T\n')
    monkeypatch.setenv(dbPaths.CONTAINER_ENV, '1')
    with pytest.raises(FileNotFoundError):
        dbPaths.find_dbpaths(str(tmp_path), 'human')
    monkeypatch.delenv(dbPaths.CONTAINER_ENV)
    assert dbPaths.find_dbpaths(str(tmp_path), 'human').endswith('human.local.dbpaths.csv')


def test_no_dbpaths_file_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        dbPaths.find_dbpaths(str(tmp_path), 'mouse', False)


# ----------------------------------------------------------------------------- dbPaths.read_dbpaths
def test_relative_paths_are_joined_with_the_folder_of_the_file(tmp_path):
    f = tmp_path / 'human.dbpaths.csv'
    f.write_text('bowtie,univec/univec_bt1,T\n\nbowtie2,/abs/genome/human-genome,T\n')
    assert dbPaths.read_dbpaths(str(f)) == [('bowtie', str(tmp_path / 'univec' / 'univec_bt1'), 'T'),
                                            ('bowtie2', '/abs/genome/human-genome', 'T')]


def test_relative_dbpaths_path_given_from_another_folder(tmp_path, monkeypatch):
    (tmp_path / 'refs').mkdir()
    (tmp_path / 'refs' / 'human.dbpaths.csv').write_text('bowtie2,trna/human-trnas,T\n')
    monkeypatch.chdir(tmp_path)
    assert dbPaths.read_dbpaths('refs/human.dbpaths.csv')[0][1] == str(tmp_path / 'refs' / 'trna' / 'human-trnas')


def test_a_replaced_tier_list_resolves_against_the_references_folder(tmp_path):
    """0.3.0 (--tiers): the replacement lies elsewhere, its relative index paths name tiers of the references."""
    (tmp_path / 'refs').mkdir()
    own = tmp_path / 'refs' / 'human.dbpaths.csv'
    own.write_text('bowtie,univec/univec_bt1,T\nbowtie2,trna/human-trnas,T\nbowtie2,genome/human-genome,T\n')
    mine = tmp_path / 'my_tiers.csv'
    mine.write_text('bowtie,univec/univec_bt1,T\nbowtie2,genome/human-genome,T\n')
    tiers = dbPaths.read_dbpaths(str(mine), str(tmp_path / 'refs'))
    assert [index for _, index, _ in tiers] == [str(tmp_path / 'refs' / 'univec' / 'univec_bt1'),
                                                str(tmp_path / 'refs' / 'genome' / 'human-genome')]
    assert dbPaths.check_tier_list(tiers, dbPaths.read_dbpaths(str(own)), 'human') == []


@pytest.mark.parametrize('lines, problem', [
    (['bowtie2,trna/human-trnas,T', 'bowtie,univec/univec_bt1,T'], 'must come first'),
    (['bowtie2,genome/human-genome,T', 'bowtie2,trna/human-trnas,T'], 'must be the last tier'),
    (['bowtie2,trna/human-trnas,T', 'bowtie2,trna/human-trnas,T'], 'listed twice'),
    (['bowtie2,mirna/human-mirna,T'], 'not a tier of the references folder'),
])
def test_a_replaced_tier_list_is_checked(tmp_path, lines, problem):
    own = tmp_path / 'human.dbpaths.csv'
    own.write_text('bowtie,univec/univec_bt1,T\nbowtie2,trna/human-trnas,T\nbowtie2,genome/human-genome,T\n')
    mine = tmp_path / 'mine.csv'
    mine.write_text('\n'.join(lines) + '\n')
    problems = dbPaths.check_tier_list(dbPaths.read_dbpaths(str(mine), str(tmp_path)),
                                       dbPaths.read_dbpaths(str(own)), 'human')
    assert any(problem in p for p in problems), problems


def test_tier_ids_are_index_basenames_without_the_species(tmp_path):
    """The tier ids of the built-in references, in mapping order."""
    lines = ['univec/univec_bt1', 'ribospike/ribospike', 'rrna/human-rrna_bt1', 'trna/human-trnas', 'mirna/human-mirna',
             'mirna/human-star', 'mirna/human-pre', 'pirna/human-pirna', 'gtf/human-ncrna', 'gtf/human-cdna',
             'genome/human-genome']
    assert [dbPaths.tier_id(index, 'human') for index in lines] == [
        'univec_bt1', 'ribospike', 'rrna_bt1', 'trnas', 'mirna', 'star', 'pre', 'pirna', 'ncrna', 'cdna', 'genome']


def test_mapreads_reads_tiers_as_absolute_paths(tmp_path, monkeypatch):
    import mapReads
    monkeypatch.setattr(mapReads.os, 'popen', lambda cmd: io.StringIO('tool version 1.0.0\n'))
    (tmp_path / 'samples.csv').write_text('S1,g,S1.fq.gz\n')
    (tmp_path / 'human.dbpaths.csv').write_text('bowtie,univec/univec_bt1,T\nbowtie2,genome/human-genome,T\n')
    args = SimpleNamespace(samples=str(tmp_path / 'samples.csv'), dbs=str(tmp_path / 'human.dbpaths.csv'),
                           threads=1, output=str(tmp_path / 'out'), multimap=100)
    m = mapReads.MapReads(args)
    assert list(m.dbs['map_cmd']) == ['bowtie', 'bowtie2']
    assert list(m.dbs['db_path']) == [str(tmp_path / 'univec' / 'univec_bt1'), str(tmp_path / 'genome' / 'human-genome')]
    assert list(m.dbs['retain']) == ['T', 'T']


# ----------------------------------------------------------------------------- builder_core.write_outputs
def test_builder_writes_one_relative_dbpaths_file(tmp_path):
    builder_core = pytest.importorskip('builder_core')
    bed = tmp_path / 'part.bed'
    bed.write_text('c\t0\t10\tf\trRNA\t+\n')
    out = tmp_path / 'assets'
    out.mkdir()
    rows = [['bowtie', 'univec/univec_bt1', 'T'], ['bowtie2', 'genome/human-genome', 'T']]
    dbpaths_csv, bedfile = builder_core.write_outputs(rows, [str(bed)], str(out), 'human')
    assert Path(dbpaths_csv).read_text() == 'bowtie,univec/univec_bt1,T\nbowtie2,genome/human-genome,T\n'
    assert sorted(p.name for p in out.iterdir()) == ['human.bed', 'human.dbpaths.csv']


# ----------------------------------------------------------------------------- update_manifest.py
def make_folder(root, local_root):
    """A references folder like today's docker/assets: two indexes, BED, MANIFEST from buildDB, and an
    older .local file whose absolute paths point to another folder (as mouse.local.dbpaths.csv does)."""
    root.mkdir()
    for index, ext in (('univec/univec_bt1', 'ebwt'), ('genome/human-genome', 'bt2')):
        (root / index).parent.mkdir(parents=True, exist_ok=True)
        for part in ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2'):
            (root / f'{index}{part}.{ext}').write_bytes(f'{index}{part}'.encode())
    (root / 'genome' / 'bowtie-build.log').write_text('log, not part of the index\n')
    (root / 'human.bed').write_text('c\t0\t10\tf\trRNA\t+\n')
    (root / 'human.local.dbpaths.csv').write_text(
        f'bowtie,{local_root}/univec/univec_bt1,T\nbowtie2,{local_root}/genome/human-genome,T\n')
    manifest = {'build': 'v1', 'species': 'human', 'sources': {'genome': {'url': 'u'}},
                'outputs': {'local_dbpaths': {'file': 'human.local.dbpaths.csv', 'sha256': 'x'},
                            'indices': [{'map_cmd': 'bowtie', 'path': 'univec/univec_bt1', 'retain': 'T'},
                                        {'map_cmd': 'bowtie2', 'path': 'genome/human-genome', 'retain': 'T'}]}}
    (root / 'MANIFEST_human.json').write_text(json.dumps(manifest))
    return root


def run_update(*args):
    return subprocess.run([sys.executable, str(UPDATE), *map(str, args)], capture_output=True, text=True)


def test_update_writes_relative_dbpaths_and_file_checksums(tmp_path):
    root = make_folder(tmp_path / 'assets', '/elsewhere/assets_v1')
    r = run_update(root, '--species', 'human')
    assert r.returncode == 0, r.stderr
    assert (root / 'human.dbpaths.csv').read_text() == 'bowtie,univec/univec_bt1,T\nbowtie2,genome/human-genome,T\n'
    m = json.loads((root / 'MANIFEST_human.json').read_text())
    assert m['sources'] == {'genome': {'url': 'u'}} and m['build'] == 'v1'          # other entries kept
    assert m['outputs_before_update']['local_dbpaths']['file'] == 'human.local.dbpaths.csv'
    files = [f['file'] for i in m['outputs']['indices'] for f in i['files']]
    assert len(files) == 12 and 'genome/bowtie-build.log' not in files
    first = m['outputs']['indices'][0]['files'][0]
    assert first['sha256'] == update_manifest.sha256_file(root / first['file'])
    assert m['outputs']['dbpaths']['file'] == 'human.dbpaths.csv'
    assert (root / 'human.local.dbpaths.csv').exists()                               # left in place
    assert run_update(root, '--species', 'human', '--check').returncode == 0


def test_check_finds_a_changed_file(tmp_path):
    root = make_folder(tmp_path / 'assets', str(tmp_path / 'assets'))
    assert run_update(root, '--species', 'human').returncode == 0
    (root / 'genome' / 'human-genome.3.bt2').write_bytes(b'changed')
    r = run_update(root, '--species', 'human', '--check')
    assert r.returncode == 1 and 'differs:  genome/human-genome.3.bt2' in r.stdout


def test_second_update_keeps_the_first_outputs(tmp_path):
    root = make_folder(tmp_path / 'assets', str(tmp_path / 'assets'))
    assert run_update(root, '--species', 'human').returncode == 0
    assert run_update(root, '--species', 'human').returncode == 0
    m = json.loads((root / 'MANIFEST_human.json').read_text())
    assert 'local_dbpaths' in m['outputs_before_update'] and 'outputs_before_update' not in m['outputs_before_update']


def test_older_file_with_other_tiers_is_refused(tmp_path):
    root = make_folder(tmp_path / 'assets', str(tmp_path / 'assets'))
    (root / 'human.local.dbpaths.csv').write_text(f'bowtie2,{root}/genome/human-genome,T\n')
    r = run_update(root, '--species', 'human')
    assert r.returncode != 0 and 'differ' in r.stderr
    assert not (root / 'human.dbpaths.csv').exists()


def test_without_manifest_paths_outside_the_folder_are_refused(tmp_path):
    root = make_folder(tmp_path / 'assets', '/elsewhere/assets_v1')
    (root / 'MANIFEST_human.json').unlink()
    r = run_update(root, '--species', 'human')
    assert r.returncode != 0 and 'outside' in r.stderr


def test_without_manifest_paths_inside_the_folder_are_made_relative(tmp_path):
    root = make_folder(tmp_path / 'assets', str(tmp_path / 'assets'))
    (root / 'MANIFEST_human.json').unlink()
    assert run_update(root, '--species', 'human').returncode == 0
    assert (root / 'human.dbpaths.csv').read_text() == 'bowtie,univec/univec_bt1,T\nbowtie2,genome/human-genome,T\n'
    assert 'note' in json.loads((root / 'MANIFEST_human.json').read_text())


def test_missing_index_files_are_an_error(tmp_path):
    root = make_folder(tmp_path / 'assets', str(tmp_path / 'assets'))
    for f in (root / 'univec').glob('univec_bt1.*'):
        f.unlink()
    r = run_update(root, '--species', 'human')
    assert r.returncode != 0 and 'no index files for univec/univec_bt1' in r.stderr


# ----------------------------------------------------------------------------- no older file names left
def test_pipeline_code_no_longer_reads_the_older_dbpaths_files():
    snakefile = (REPO / 'docker' / 'Snakefile').read_text()
    dockerfile = (REPO / 'docker' / 'Dockerfile').read_text()
    entrypoint = (REPO / 'docker' / 'run-pipeline.sh').read_text()
    assert 'dbpaths = find_dbpaths(dbdir, species)' in snakefile and 'local.dbpaths' not in snakefile
    assert 'docker.dbpaths' not in dockerfile and 'ENV RIBOMARKER_CONTAINER=1' in dockerfile
    assert 'COPY scripts/ /pipeline/scripts/' in dockerfile
    assert 'docker.dbpaths' not in entrypoint and 'ln -s' not in entrypoint


# ----------------------------------------------------------------------------- mouse references
def make_gtrnadb_tar(tmp_path, code, species_name, loci):
    """A GtRNAdb-like archive: <code>-tRNAs.fa (+ mature, filtered) and <code>-tRNAs.bed."""
    import tarfile
    src = tmp_path / f'src_{code}'
    src.mkdir()
    fa = ''.join(f'>{species_name}_tRNA-Ala-AGC-{i}-1 (x)\nGGGGGTATAGCTCAGTGGTAGAGCGCGTGCTTAGCATGC\n' for i in loci)
    (src / f'{code}-tRNAs.fa').write_text(fa)
    (src / f'{code}-mature-tRNAs.fa').write_text(fa)
    (src / f'{code}-filtered-tRNAs.fa').write_text(fa)
    (src / f'{code}-tRNAs.bed').write_text(
        ''.join(f'chr1\t{100 * i}\t{100 * i + 72}\ttRNA-Ala-AGC-{i}-1\t1000\t+\n' for i in loci))
    tar_path = tmp_path / f'{code}-tRNAs.tar.gz'
    with tarfile.open(tar_path, 'w:gz') as tar:
        for f in src.iterdir():
            tar.add(f, arcname=f.name)
    return tar_path


def test_trna_files_come_from_the_archive_not_from_the_folder(tmp_path, monkeypatch):
    """A trna/ folder that already holds hg38 files must not give the mouse build human tRNAs."""
    builder_core = pytest.importorskip('builder_core')
    monkeypatch.setattr(builder_core, '_bowtie2_build', lambda *a, **k: None)
    out = tmp_path / 'assets'
    human_tar = make_gtrnadb_tar(tmp_path, 'hg38', 'Homo_sapiens', [1, 2])
    mouse_tar = make_gtrnadb_tar(tmp_path, 'mm39', 'Mus_musculus', [7])
    builder_core.build_gtrnadb(str(human_tar), str(out), 1, 'human')          # hg38 files now in trna/
    rows, beds = builder_core.build_gtrnadb(str(mouse_tar), str(out), 1, 'mouse')
    assert rows == [['bowtie2', 'trna/mouse-trnas', 'T']]
    trna_bed, pre_bed = (Path(b).read_text() for b in beds)
    assert 'Mus_musculus_tRNA-Ala-AGC-7' in trna_bed and 'Homo_sapiens' not in trna_bed
    assert pre_bed == '1\t700\t772\ttRNA-Ala-AGC-7-1\tpre-tRNA\t+\n'


def make_build(root, species='mouse'):
    """A finished separate build: species files, shared indexes, unprefixed intermediates, MANIFEST."""
    files = {f'{species}.dbpaths.csv': 'bowtie,univec/univec_bt1,T\nbowtie2,trna/mouse-trnas,T\n',
             f'{species}.bed': 'Mus_musculus_tRNA-Ala-AGC-7\t10\t50\tx\ttRNA\t+\n',
             f'PROVENANCE_{species}.md': '# prov\n', 'trna/mm39-tRNAs.bed': 'bed\n',
             'trna/trnas.fa': 'unprefixed intermediate\n', 'trna/bowtie-build.log': 'log\n'}
    for part in ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2'):
        files[f'trna/mouse-trnas{part}.bt2'] = f'new mouse {part}'
        files[f'univec/univec_bt1{part}.ebwt'] = f'univec {part}'
    for rel, text in files.items():
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(text)
    manifest = {'species': species, 'sources': {'trna': {'url': 'https://x/Mmusc39/mm39-tRNAs.tar.gz'}}}
    (root / f'MANIFEST_{species}.json').write_text(json.dumps(manifest))
    return root


def test_install_copies_only_the_species_and_keeps_what_it_replaces(tmp_path):
    build = make_build(tmp_path / 'build')
    target = tmp_path / 'assets'
    (target / 'trna').mkdir(parents=True)
    (target / 'univec').mkdir()
    (target / 'trna' / 'trnas.fa').write_text('human intermediate\n')
    (target / 'trna' / 'mouse-trnas.1.bt2').write_text('old mouse')
    (target / 'human.bed').write_text('human\n')
    (target / 'mouse.local.dbpaths.csv').write_text('bowtie2,/old/trna/mouse-trnas,T\n')
    for part in ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2'):
        (target / 'univec' / f'univec_bt1{part}.ebwt').write_text(f'univec {part}')
    r = subprocess.run([sys.executable, str(REPO / 'buildDB' / 'install_species.py'), str(build), str(target),
                        '--species', 'mouse'], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (target / 'trna' / 'mouse-trnas.1.bt2').read_text() == 'new mouse .1'
    assert (target / 'trna' / 'trnas.fa').read_text() == 'human intermediate\n'     # unprefixed: untouched
    assert (target / 'human.bed').read_text() == 'human\n'
    assert (target / 'trna' / 'mm39-tRNAs.bed').exists() and not (target / 'mouse.local.dbpaths.csv').exists()
    kept = list(build.glob('_replaced_mouse_*'))
    assert len(kept) == 1 and (kept[0] / 'trna' / 'mouse-trnas.1.bt2').read_text() == 'old mouse'
    assert (kept[0] / 'mouse.local.dbpaths.csv').exists()
    assert 'univec/univec_bt1: same as installed' in r.stdout
    m = json.loads((target / 'MANIFEST_mouse.json').read_text())
    assert 'installed' in m and len(m['outputs']['indices']) == 2
    assert run_update(target, '--species', 'mouse', '--check').returncode == 0


def test_install_reports_a_different_shared_index_and_keeps_the_installed_one(tmp_path):
    build = make_build(tmp_path / 'build')
    (build / 'univec' / 'univec_bt1.3.ebwt').write_text('newer UniVec download')
    target = tmp_path / 'assets'
    (target / 'univec').mkdir(parents=True)
    for part in ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2'):
        (target / 'univec' / f'univec_bt1{part}.ebwt').write_text(f'univec {part}')
    r = subprocess.run([sys.executable, str(REPO / 'buildDB' / 'install_species.py'), str(build), str(target),
                        '--species', 'mouse'], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert 'univec/univec_bt1: DIFFERS' in r.stdout
    assert (target / 'univec' / 'univec_bt1.3.ebwt').read_text() == 'univec .3'
    assert run_update(target, '--species', 'mouse', '--check').returncode == 0

