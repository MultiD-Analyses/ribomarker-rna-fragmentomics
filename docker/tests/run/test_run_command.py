"""The run command's decisions, on small made-up folders; no pipeline run.
Run from the repository root: python -m pytest docker/tests/run"""
import gzip
import json
import os
import signal
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'docker' / 'scripts'))
import runCommand as rc  # noqa: E402
from runCommand import RunError  # noqa: E402


def fastq(folder, *names):
    """Small valid gzip FASTQ files."""
    folder.mkdir(parents=True, exist_ok=True)
    for name in names:
        with gzip.open(folder / name, 'wt') as f:
            f.write('@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n')
    return folder


def write_sheet(folder, text, name='samples.csv'):
    """A sample sheet in folder (the FASTQ files are made separately, with fastq())."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / name).write_text(text)
    return folder / name


# ----------------------------------------------------------------------------- the sample sheet
def test_paired_end_sheet_with_free_file_names(tmp_path):
    fastq(tmp_path / 'in', 'a_R1.fastq.gz', 'a_R2.fq.gz')
    sheet = write_sheet(tmp_path / 'in', 'sample,fastq_file,fastq_file2\nHealthy_1,a_R1.fastq.gz,a_R2.fq.gz\n')
    samples, mode = rc.read_sheet(sheet)
    assert samples == [('Healthy_1', 'a_R1.fastq.gz', 'a_R2.fq.gz')] and mode == 'pe'
    assert rc.sheet_rows(samples) == [['Healthy_1', 'RiboMarker', 'a_R1.fastq.gz', 'a_R2.fq.gz']]


@pytest.mark.parametrize('header, rows', [('sample,fastq_file', 'B,b.fastq.gz\nA,a.fastq.gz\n'),
                                          ('sample,fastq_file,fastq_file2', 'B,b.fastq.gz,\nA,a.fastq.gz,\n')])
def test_single_end_without_or_with_an_empty_fastq_file2_column_keeps_the_order(tmp_path, header, rows):
    fastq(tmp_path / 'in', 'a.fastq.gz', 'b.fastq.gz')
    samples, mode = rc.read_sheet(write_sheet(tmp_path / 'in', header + '\n' + rows))
    assert samples == [('B', 'b.fastq.gz', ''), ('A', 'a.fastq.gz', '')] and mode == 'se'
    assert rc.sheet_rows(samples)[0] == ['B', 'RiboMarker', 'b.fastq.gz']


def test_excel_sheet_semicolons_byte_order_mark_any_column_order_extra_columns(tmp_path):
    fastq(tmp_path / 'in', 'x.fq.gz')
    sheet = tmp_path / 'in' / 'samples.csv'
    sheet.write_bytes('\ufeffage;fastq_file;sample;group\r\n54; x.fq.gz ; P1 ;case\r\n;;;\r\n'.encode('utf-8'))
    assert rc.read_sheet(sheet) == ([('P1', 'x.fq.gz', '')], 'se')


@pytest.mark.parametrize('text, words', [
    ('Sample,fastq_file\nA,a.fastq.gz\n', 'missing: sample'),
    ('sample,fastq\nA,a.fastq.gz\n', 'missing: fastq_file'),
    ('sample,fastq_file,sample\nA,a.fastq.gz,B\n', 'given twice: sample'),
    ('sample,fastq_file\n', 'no rows with samples'),
])
def test_the_header_must_name_the_columns(tmp_path, text, words):
    fastq(tmp_path / 'in', 'a.fastq.gz')
    with pytest.raises(RunError, match=words):
        rc.read_sheet(write_sheet(tmp_path / 'in', text))


@pytest.mark.parametrize('name', ['Healthy_1', '1_healthy', '_x', 'S001', 'Nebraska', 'a-b_C'])
def test_usable_sample_names(tmp_path, name):
    fastq(tmp_path / 'in', 'a.fastq.gz')
    assert rc.read_sheet(write_sheet(tmp_path / 'in', f'sample,fastq_file\n{name},a.fastq.gz\n'))[0][0][0] == name


@pytest.mark.parametrize('name, words', [
    ('plasma 1', 'only letters, digits'), ('plasma.A', 'only letters, digits'), ('p(2)', 'only letters, digits'),
    ('-p1', 'must not start with'), ('001', 'numbers or missing values'), ('1e5', 'numbers or missing values'),
    ('NA', 'numbers or missing values'), ('true', 'numbers or missing values'), ('lib_NEB', 'NEB'), ('SMART1', 'SMART'),
])
def test_unusable_sample_names_are_refused(tmp_path, name, words):
    fastq(tmp_path / 'in', 'a.fastq.gz')
    with pytest.raises(RunError, match=words):
        rc.read_sheet(write_sheet(tmp_path / 'in', f'sample,fastq_file\n{name},a.fastq.gz\n'))


def test_sample_names_given_twice_also_in_another_case_are_refused(tmp_path):
    fastq(tmp_path / 'in', 'a.fastq.gz', 'b.fastq.gz')
    with pytest.raises(RunError, match='more than once(.|\n)*Ctrl / ctrl'):
        rc.read_sheet(write_sheet(tmp_path / 'in', 'sample,fastq_file\nCtrl,a.fastq.gz\nctrl,b.fastq.gz\n'))


def test_mixed_paired_and_single_end_is_refused(tmp_path):
    fastq(tmp_path / 'in', 'a1.fastq.gz', 'a2.fastq.gz', 'b1.fastq.gz')
    with pytest.raises(RunError, match='either paired-end(.|\n)*Without fastq_file2: B'):
        rc.read_sheet(write_sheet(tmp_path / 'in', 'sample,fastq_file,fastq_file2\n'
                                                   'A,a1.fastq.gz,a2.fastq.gz\nB,b1.fastq.gz,\n'))


@pytest.mark.parametrize('file, words', [
    ('sub/a.fastq.gz', 'without a path'), ('C:\\data\\a.fastq.gz', 'without a path'),
    ('my reads.fastq.gz', 'no spaces'), ('a.fastq', 'gzip-compressed'), ('missing.fastq.gz', 'not found in'),
])
def test_unusable_file_entries_are_refused(tmp_path, file, words):
    fastq(tmp_path / 'in', 'a.fastq')
    with pytest.raises(RunError, match=words):
        rc.read_sheet(write_sheet(tmp_path / 'in', f'sample,fastq_file\nA,{file}\n'))


def test_a_file_listed_twice_or_as_both_reads_is_refused(tmp_path):
    fastq(tmp_path / 'in', 'a.fastq.gz', 'b.fastq.gz')
    with pytest.raises(RunError, match='more than once(.|\n)*a.fastq.gz'):
        rc.read_sheet(write_sheet(tmp_path / 'in', 'sample,fastq_file,fastq_file2\nA,a.fastq.gz,a.fastq.gz\n'))
    with pytest.raises(RunError, match='more than once(.|\n)*b.fastq.gz'):
        rc.read_sheet(write_sheet(tmp_path / 'in', 'sample,fastq_file\nA,b.fastq.gz\nB,b.fastq.gz\n'))


def test_all_problems_are_reported_together_with_row_numbers(tmp_path):
    fastq(tmp_path / 'in', 'a.fastq.gz')
    with pytest.raises(RunError) as error:
        rc.read_sheet(write_sheet(tmp_path / 'in', 'sample,fastq_file\nA,a.fastq.gz\nB,\n,a.fastq.gz\nx y,c.fastq.gz\n'))
    text = str(error.value)
    for words in ('row 3: no fastq_file for B', 'row 4: no sample name', 'x y', 'c.fastq.gz'):
        assert words in text, words


def test_a_sheet_not_saved_as_utf8_is_refused_with_a_hint(tmp_path):
    fastq(tmp_path / 'in', 'a.fastq.gz')
    sheet = tmp_path / 'in' / 'samples.csv'
    sheet.write_bytes('sample;fastq_file;ort\nA;a.fastq.gz;Malmö\n'.encode('cp1252'))   # an Excel "CSV" export
    with pytest.raises(RunError, match='CSV UTF-8'):
        rc.read_sheet(sheet)


def test_a_missing_sheet_is_refused(tmp_path):
    with pytest.raises(RunError, match='sample sheet not found'):
        rc.read_sheet(tmp_path / 'none.csv')


# ----------------------------------------------------------------------------- configuration
def write_yaml(path, values):
    path.write_text(yaml.safe_dump(values))
    return path


def test_allowed_settings_are_taken(tmp_path):
    cfg = write_yaml(tmp_path / 'c.yaml', {'multimap': 50, 'remove_antisense': False, 'collapse': False})
    assert rc.user_config(cfg) == {'multimap': 50, 'remove_antisense': False, 'collapse': False}


@pytest.mark.parametrize('values, word', [
    ({'dbdir': '/x'}, 'set by run'), ({'threads': 4}, 'set by run'), ({'paired_end': True}, 'set by run'),
    ({'collapse': True}, 'multimapper'), ({'trim_fastq': False}, 'trimming'), ({'trim_only': True}, 'QC report'),
    ({'mapper': 'bowtie'}, 'not available under run'), ({'multimapp': 10}, 'not available under run'),
])
def test_settings_run_cannot_use_are_refused(tmp_path, values, word):
    with pytest.raises(RunError, match=word):
        rc.user_config(write_yaml(tmp_path / 'c.yaml', values))


@pytest.mark.parametrize('text, word', [
    ('remove_multimap: "False"', 'not True or False'), ('qc_report: 1', 'not True or False'),
    ('multimap: ten', 'whole number of 1'), ('multimap: 0', 'whole number of 1'), ('multimap: true', 'whole number of 1'),
    ('keep_types: miRNA, tRNA', 'without spaces'), ('keep_types: [miRNA]', 'without spaces'),
])
def test_values_of_allowed_settings_are_checked(tmp_path, text, word):
    (tmp_path / 'c.yaml').write_text(text + '\n')
    with pytest.raises(RunError, match=word):
        rc.user_config(tmp_path / 'c.yaml')


def test_valid_values_of_allowed_settings_are_taken(tmp_path):
    (tmp_path / 'c.yaml').write_text('remove_multimap: false\nqc_report: no\nmultimap: 1\n'
                                     'keep_types: miRNA,tRNA,Mt_tRNA,pre-mRNA\n')
    assert rc.user_config(tmp_path / 'c.yaml') == {'remove_multimap': False, 'qc_report': False, 'multimap': 1,
                                                   'keep_types': 'miRNA,tRNA,Mt_tRNA,pre-mRNA'}


def test_a_config_file_that_is_not_yaml_is_refused(tmp_path):
    (tmp_path / 'c.yaml').write_text('multimap: [1\n')
    with pytest.raises(RunError, match='cannot be read as a YAML file'):
        rc.user_config(tmp_path / 'c.yaml')


def test_pipeline_config_uses_absolute_folders_and_the_bare_run_id(tmp_path):
    c = rc.pipeline_config({'multimap': 50}, input_dir=Path('/data'), references=Path('/refs'), run_id='ribomarker',
                           sheet=Path('/out/ribomarker/samples.csv'), species='human', threads=8, mode='se',
                           title='RiboMarker run 2026-10-09')
    assert c['output_prefix'] == 'ribomarker' and c['input_dir'] == '/data/' and c['dbdir'] == '/refs/'
    assert c['samples'] == '/out/ribomarker/samples.csv' and c['paired_end'] is False
    assert c['collapse'] is False and c['remove_multimap'] is True and c['multimap'] == 50
    assert c['report_title'] == 'RiboMarker run 2026-10-09'


# ----------------------------------------------------------------------------- references
def references(folder, with_manifest=False):
    (folder / 'genome').mkdir(parents=True)
    (folder / 'genome' / 'human-genome.1.bt2').write_text('x')
    (folder / 'human.dbpaths.csv').write_text('bowtie2,genome/human-genome,T\n')
    (folder / 'human.bed').write_text('c\t1\t2\tf\tpre-mRNA\t+\n')
    if with_manifest:
        (folder / 'MANIFEST_human.json').write_text('{}')
    return folder


def test_references_without_manifest_give_a_warning(tmp_path):
    assert 'no MANIFEST_human.json' in rc.check_references(references(tmp_path / 'r'), 'human')[0]
    assert rc.check_references(references(tmp_path / 'r2', with_manifest=True), 'human') == []


def test_incomplete_references_are_refused(tmp_path):
    r = references(tmp_path / 'r')
    (r / 'genome' / 'human-genome.1.bt2').unlink()
    with pytest.raises(RunError, match='incomplete'):
        rc.check_references(r, 'human')
    with pytest.raises(RunError, match='mouse.dbpaths.csv not found'):
        rc.check_references(r, 'mouse')


def test_another_species_needs_own_references(tmp_path):
    one_sample(tmp_path)
    with pytest.raises(RunError, match='built-in references are human and mouse'):
        rc.prepare(args(tmp_path, species='rat', references=str(rc.REFERENCES)))
    with pytest.raises(RunError, match='--species'):
        rc.prepare(args(tmp_path, species='rat 2'))


def test_another_species_with_own_references(tmp_path):
    one_sample(tmp_path)
    refs = tmp_path / 'refs'
    for ending in ('dbpaths.csv', 'bed'):
        (refs / f'rat.{ending}').write_text((refs / f'human.{ending}').read_text())
    p = rc.prepare(args(tmp_path, species='rat'))
    assert yaml.safe_load(p['config_text'])['species'] == 'rat'
    assert any('tested with human and mouse only' in w for w in p['warnings'])


# ----------------------------------------------------------------------------- run folder: new, resume, refuse
def test_run_state(tmp_path):
    run = tmp_path / 'ribomarker'
    assert rc.run_state(run, 'S', 'threads: 8\nx: 1\n') == 'new'
    run.mkdir()
    with pytest.raises(RunError, match='not made by run'):
        rc.run_state(run, 'S', 'threads: 8\nx: 1\n')
    (run / 'samples.csv').write_text('S')
    (run / 'run_config.yaml').write_text('threads: 8\nx: 1\n')
    (run / 'run_status.json').write_text(json.dumps({'status': 'failed'}))
    assert rc.run_state(run, 'S', 'threads: 4\nx: 1\n') == 'resume'          # other threads: still resumable
    with pytest.raises(RunError, match='other samples or settings'):
        rc.run_state(run, 'S', 'threads: 8\nx: 2\n')
    (run / 'run_status.json').write_text(json.dumps({'status': 'finished', 'verdict': 'RESULTS VALID'}))
    with pytest.raises(RunError, match='finished run'):
        rc.run_state(run, 'S', 'threads: 8\nx: 1\n')


# ----------------------------------------------------------------------------- verdict
def qc_folder(run, statuses, warnings=()):
    (run / 'qc').mkdir(parents=True)
    rows = ''.join(f'C{i}\t{s}\tS1\t\t1\t1\t0\tx\n' for i, s in enumerate(statuses))
    (run / 'qc' / 'consistency_checks.tsv').write_text('check\tstatus\tsample\tsubject\texpected\tobserved\t'
                                                        'difference\tidentity\n' + rows)
    (run / 'qc' / 'qc_metrics.json').write_text(json.dumps({'warnings': [{'source': w} for w in warnings]}))


def test_verdict(tmp_path):
    qc_folder(tmp_path / 'a', ['equal', 'equal', 'not applicable'], ['x'])
    assert rc.verdict(tmp_path / 'a') == (True, 'RESULTS VALID: all 2 consistency checks equal '
                                                '(2 equal, 1 not applicable; 1 warning(s))')
    qc_folder(tmp_path / 'b', ['equal', 'differs', 'not evaluable'])
    valid, line = rc.verdict(tmp_path / 'b')
    assert not valid and line.startswith('RESULTS NOT VALID: 2 of 3 consistency checks failed')
    (tmp_path / 'c').mkdir()
    assert rc.verdict(tmp_path / 'c')[1].startswith('RESULTS NOT CHECKED')


# ----------------------------------------------------------------------------- prepare and main, without Snakemake
def one_sample(tmp_path, extra_header='', extra_value=''):
    """A single-end sample A with its sheet in tmp_path/in, and references in tmp_path/refs; returns the sheet."""
    fastq(tmp_path / 'in', 'A_R1.fastq.gz')
    references(tmp_path / 'refs', with_manifest=True)
    return write_sheet(tmp_path / 'in', 'sample,fastq_file' + extra_header + '\nA,A_R1.fastq.gz' + extra_value + '\n')


def args(tmp_path, **extra):
    base = dict(sheet=str(tmp_path / 'in' / 'samples.csv'), output=str(tmp_path / 'out'), species='human',
                threads=1, config=None, run_id=None, references=str(tmp_path / 'refs'), unlock=False)
    base.update(extra)
    return SimpleNamespace(**base)


def command(tmp_path, *more):
    return [str(tmp_path / 'in' / 'samples.csv'), str(tmp_path / 'out'), '--references', str(tmp_path / 'refs'),
            '--threads', '1', *more]


def test_prepare_a_new_run(tmp_path):
    fastq(tmp_path / 'in', 'A_R1.fastq.gz', 'A_R2.fastq.gz')
    write_sheet(tmp_path / 'in', 'sample,fastq_file,fastq_file2\nA,A_R1.fastq.gz,A_R2.fastq.gz\n')
    references(tmp_path / 'refs', with_manifest=True)
    p = rc.prepare(args(tmp_path))
    assert p['state'] == 'new' and p['mode'] == 'pe' and p['run_dir'] == (tmp_path / 'out' / 'ribomarker').resolve()
    assert p['input_dir'] == (tmp_path / 'in').resolve()
    assert p['title'] == f'RiboMarker run {p["first_start"][:10]}'
    config = yaml.safe_load(p['config_text'])
    assert config['output_prefix'] == 'ribomarker' and config['paired_end'] is True
    assert config['samples'] == str(p['run_dir'] / 'samples.csv')
    assert config['input_dir'] == str((tmp_path / 'in').resolve()) + '/'
    assert rc.prepare(args(tmp_path, run_id='plasma_batch1'))['title'] == 'plasma_batch1'


def test_prepare_refuses_a_missing_config_file(tmp_path):
    one_sample(tmp_path)
    with pytest.raises(RunError, match='--config: file not found'):
        rc.prepare(args(tmp_path, config=str(tmp_path / 'missing.yaml')))


def test_there_is_no_metadata_option_any_more(capsys):
    with pytest.raises(SystemExit):
        rc.parse_args(['in.csv', 'out', '--metadata', 'm.csv'])
    assert 'unrecognized arguments: --metadata' in capsys.readouterr().err


@pytest.mark.parametrize('value', ['0', '-2', 'eight'])
def test_threads_must_be_a_whole_number_of_1_or_more(value, capsys):
    with pytest.raises(SystemExit) as stop:
        rc.parse_args(['in.csv', 'out', '--threads', value])
    assert stop.value.code == 2 and 'whole number of 1 or more' in capsys.readouterr().err


def test_threads_default_to_8(tmp_path):
    assert rc.parse_args(['in.csv', 'out']).threads == 8


def test_a_finished_run_is_refused_before_the_gzip_check(tmp_path, monkeypatch):
    one_sample(tmp_path)
    run = tmp_path / 'out' / 'ribomarker'
    run.mkdir(parents=True)
    (run / 'run_status.json').write_text('{"status": "finished", "verdict": "RESULTS VALID"}')
    checked = []
    monkeypatch.setattr(rc, 'broken_gzip', lambda files, threads: checked.append(files) or [])
    with pytest.raises(RunError, match='holds a finished run'):
        rc.prepare(args(tmp_path))
    assert checked == []


def test_prepare_refuses_output_inside_input_spaces_and_damaged_files(tmp_path):
    one_sample(tmp_path)
    with pytest.raises(RunError, match='inside it'):
        rc.prepare(args(tmp_path, output=str(tmp_path / 'in' / 'results')))
    with pytest.raises(RunError, match='space'):
        rc.prepare(args(tmp_path, output=str(tmp_path / 'my results')))
    with pytest.raises(RunError, match='--run-id'):
        rc.prepare(args(tmp_path, run_id='2026 run'))
    (tmp_path / 'in' / 'A_R1.fastq.gz').write_bytes(gzip.compress(b'@r\nACGT\n+\nIIII\n')[:-6])
    with pytest.raises(RunError, match='damaged or incomplete'):
        rc.prepare(args(tmp_path))


def test_main_writes_the_run_folder_copies_the_sheet_and_gives_the_verdict(tmp_path, monkeypatch, capsys):
    sheet = one_sample(tmp_path, ',age', ',54')

    def fake_snakemake(p, unlock):                     # what a finished pipeline leaves behind
        qc_folder(p['run_dir'], ['equal', 'equal'])
        return 0
    monkeypatch.setattr(rc, 'run_snakemake', fake_snakemake)
    assert rc.main(command(tmp_path)) == 0
    run = tmp_path / 'out' / 'ribomarker'
    assert (run / 'sample_sheet.csv').read_text() == sheet.read_text() == 'sample,fastq_file,age\nA,A_R1.fastq.gz,54\n'
    assert (run / 'samples.csv').read_text() == 'A,RiboMarker,A_R1.fastq.gz\n'
    status = json.loads((run / 'run_status.json').read_text())
    assert status['status'] == 'finished' and status['verdict'].startswith('RESULTS VALID')
    version = status['pipeline_version']
    assert version['version'] == (REPO / 'docker' / 'VERSION').read_text().strip()
    assert len(version['code_fingerprint']) == 64 and version['git_commit']
    assert 'RESULTS VALID' in capsys.readouterr().out
    # the same command again: a finished run is refused (exit 2), not repeated
    assert rc.main(command(tmp_path)) == 2


def test_failed_pipeline_can_be_resumed_with_the_same_command(tmp_path, monkeypatch):
    one_sample(tmp_path)
    monkeypatch.setattr(rc, 'run_snakemake', lambda p, unlock: 1)
    assert rc.main(command(tmp_path)) == 1
    status = json.loads((tmp_path / 'out' / 'ribomarker' / 'run_status.json').read_text())
    assert status['status'] == 'failed' and 'stopped at a step' in status['verdict']
    first = status['first_started_utc']

    def finish(p, unlock):
        assert p['state'] == 'resume'
        qc_folder(p['run_dir'], ['equal'])
        return 0
    monkeypatch.setattr(rc, 'run_snakemake', finish)
    assert rc.main(command(tmp_path)) == 0
    assert json.loads((tmp_path / 'out' / 'ribomarker' / 'run_status.json').read_text())['first_started_utc'] == first


def test_another_sample_in_the_sheet_stops_the_resume(tmp_path, monkeypatch):
    one_sample(tmp_path)
    monkeypatch.setattr(rc, 'run_snakemake', lambda p, unlock: 1)
    assert rc.main(command(tmp_path)) == 1
    fastq(tmp_path / 'in', 'B_R1.fastq.gz')
    write_sheet(tmp_path / 'in', 'sample,fastq_file\nA,A_R1.fastq.gz\nB,B_R1.fastq.gz\n')
    assert rc.main(command(tmp_path)) == 2                # other samples or settings: refused, not resumed


def test_interrupted_run_is_marked_and_resumed_with_the_same_command(tmp_path, monkeypatch, capsys):
    one_sample(tmp_path)

    def ctrl_c(p, unlock):
        raise KeyboardInterrupt
    monkeypatch.setattr(rc, 'run_snakemake', ctrl_c)
    assert rc.main(command(tmp_path)) == 130
    assert 'INTERRUPTED: repeat the same command to resume' in capsys.readouterr().out
    status_file = tmp_path / 'out' / 'ribomarker' / 'run_status.json'
    assert json.loads(status_file.read_text())['status'] == 'interrupted'

    def finish(p, unlock):
        assert p['state'] == 'resume'
        qc_folder(p['run_dir'], ['equal'])
        return 0
    monkeypatch.setattr(rc, 'run_snakemake', finish)
    assert rc.main(command(tmp_path)) == 0
    assert json.loads(status_file.read_text())['status'] == 'finished'


@pytest.mark.parametrize('n_bytes, minutes, size, duration', [
    (3e5, 0.2, '0 MB', '1 min'), (5e7, 12.4, '50 MB', '12 min'), (3e8, 37.5, '0.3 GB', '38 min'),
    (4.2e9, 525, '4.2 GB', '8.8 h')])
def test_size_and_duration_are_readable(n_bytes, minutes, size, duration):
    assert rc.size_text(n_bytes) == size and rc.duration_text(minutes) == duration


def test_note_when_the_gene_table_is_not_made(tmp_path, monkeypatch, capsys):
    one_sample(tmp_path)
    cfg = write_yaml(tmp_path / 'c.yaml', {'remove_multimap': False})

    def finish(p, unlock):
        qc_folder(p['run_dir'], ['equal', 'not applicable'])
        return 0
    monkeypatch.setattr(rc, 'run_snakemake', finish)
    assert rc.main(command(tmp_path, '--config', str(cfg))) == 0
    assert 'note: no gene table' in capsys.readouterr().out


def test_a_signal_to_run_alone_reaches_snakemake(tmp_path):
    """As in Docker without a terminal: only run gets the Ctrl-C and passes it on to Snakemake."""
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    (bin_dir / 'snakemake').write_text('#!/bin/bash\ntrap \'echo "stopping on SIGINT"; exit 1\' INT\n'
                                       'echo started\nfor i in $(seq 100); do sleep 0.1; done\n')
    (bin_dir / 'snakemake').chmod(0o755)
    run_dir = tmp_path / 'out' / 'run'
    run_dir.mkdir(parents=True)
    driver = tmp_path / 'driver.py'
    driver.write_text(f'import sys\nfrom pathlib import Path\nsys.path.insert(0, {str(REPO / "docker" / "scripts")!r})\n'
                      'import runCommand as rc\n'
                      f'p = dict(run_dir=Path({str(run_dir)!r}), output=Path({str(run_dir.parent)!r}), threads=1)\n'
                      'try:\n    rc.run_snakemake(p, False)\nexcept KeyboardInterrupt:\n    print("interrupted")\n')
    env = dict(os.environ, PATH=f'{bin_dir}{os.pathsep}{os.environ["PATH"]}')
    child = subprocess.Popen([sys.executable, '-u', str(driver)], env=env, stdout=subprocess.PIPE, text=True)
    assert child.stdout.readline() == 'started\n'      # the fake Snakemake runs
    child.send_signal(signal.SIGINT)                    # to run only, not to its process group
    rest, _ = child.communicate(timeout=20)
    assert 'stopping on SIGINT' in rest and 'interrupted' in rest
    assert 'stopping on SIGINT' in (run_dir / 'run.log').read_text()


def test_container_memory_from_cgroup_v1_or_v2(tmp_path, monkeypatch):
    monkeypatch.delenv('RIBOMARKER_CONTAINER', raising=False)
    assert rc.container_memory_gb() is None                        # outside a container
    monkeypatch.setenv('RIBOMARKER_CONTAINER', '1')
    v1 = tmp_path / 'v1.stat'
    v1.write_text('cache 9000000000\nrss 1000\ntotal_cache 9000000000\ntotal_rss 3900000000\n')
    monkeypatch.setattr(rc, 'MEMORY_STAT', [(str(tmp_path / 'missing'), 'anon'), (str(v1), 'total_rss')])
    assert rc.container_memory_gb() == 3.9                          # the file cache is not counted
    monkeypatch.setattr(rc, 'MEMORY_STAT', [(str(tmp_path / 'missing'), 'anon')])
    assert rc.container_memory_gb() is None


def test_help_names_both_ways_to_run_and_the_exit_codes(capsys):
    with pytest.raises(SystemExit):
        rc.parse_args(['--help'])
    text = capsys.readouterr().out
    for expected in ('docker run', 'PowerShell', 'python docker/scripts/main.py run', 'fastq_file2',
                     '130 = interrupted', 'remove_antisense', 'run --demo OUTPUT'):
        assert expected in text, expected


def test_demo_runs_its_reads_and_compares_the_truth(tmp_path, monkeypatch, capsys):
    demo = tmp_path / 'demo'
    fastq(demo / 'fastq', 'S1_R1_001.fastq.gz', 'S1_R2_001.fastq.gz')
    write_sheet(demo / 'fastq', 'sample,fastq_file,fastq_file2\nS1,S1_R1_001.fastq.gz,S1_R2_001.fastq.gz\n')
    references(demo / 'references', with_manifest=True)
    (demo / 'truth.tsv').write_text('mode\n')
    monkeypatch.setattr(rc, 'DEMO', demo)

    def finish(p, unlock):
        assert p['run_id'] == 'demo' and p['input_dir'] == demo / 'fastq'
        qc_folder(p['run_dir'], ['equal'])
        return 0
    monkeypatch.setattr(rc, 'run_snakemake', finish)
    monkeypatch.setattr(rc.checkTruth, 'check', lambda truth, mode, outdir: ([['PASS'] * 7] * 3, 'table'))
    assert rc.main(['--demo', str(tmp_path / 'out'), '--threads', '1']) == 0
    assert 'DEMO PASSED: all 3 known results' in capsys.readouterr().out

    monkeypatch.setattr(rc.checkTruth, 'check', lambda truth, mode, outdir:
                        ([['PASS'] * 7, ['FAIL', 'genecounts', 'S1', 'G', 5, 4, 'c']], 'table'))
    assert rc.main(['--demo', str(tmp_path / 'out2'), '--threads', '1']) == 1
    assert 'DEMO FAILED: 1 of 2 known results differ' in capsys.readouterr().out
    assert rc.main(['--demo', str(tmp_path / 'out2'), '--threads', '1']) == 2     # finished: not repeated
    assert 'DEMO NOT STARTED' in capsys.readouterr().out


def test_demo_without_demo_data(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(rc, 'DEMO', tmp_path / 'none')
    assert rc.main(['--demo', str(tmp_path / 'out')]) == 2
    assert 'no demo data' in capsys.readouterr().out


def test_scripts_run_unbuffered(tmp_path, monkeypatch):
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    (bin_dir / 'snakemake').write_text('#!/bin/bash\necho "unbuffered=$PYTHONUNBUFFERED"\n')
    (bin_dir / 'snakemake').chmod(0o755)
    monkeypatch.setenv('PATH', f'{bin_dir}{os.pathsep}{os.environ["PATH"]}')
    run_dir = tmp_path / 'out' / 'run'
    run_dir.mkdir(parents=True)
    assert rc.run_snakemake(dict(run_dir=run_dir, output=run_dir.parent, threads=1), False) == 0
    assert 'unbuffered=1' in (run_dir / 'run.log').read_text()


def test_refusal_before_start_writes_nothing_into_the_run_folder(tmp_path, capsys):
    one_sample(tmp_path)
    write_sheet(tmp_path / 'in', 'sample,fastq_file\nA,A_R2.fastq.gz\n')
    assert rc.main(command(tmp_path)) == 2
    assert 'not found in' in capsys.readouterr().out
    assert not (tmp_path / 'out' / 'ribomarker').exists()


def test_main_dispatches_run(tmp_path):
    text = (REPO / 'docker' / 'scripts' / 'main.py').read_text()
    assert 'if len(sys.argv) > 1 and sys.argv[1] == "run":' in text and 'from runCommand import main' in text
