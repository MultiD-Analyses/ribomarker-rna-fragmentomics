"""Consistency checks C1-C10 and Kind A rules: all applicable checks are
'equal' on the synthetic data, and each check reports 'differs' when an input is changed
so that its identity no longer holds."""
import json
import re
import shutil

import pytest
import yaml

import qcChecks
import qcCollect
from conftest import REPO, prepare_run, run_collector


def statuses(report):
    """{(check id, sample, subject): status}"""
    return {(c['id'], c['sample'], c['subject']): c['status'] for c in report['checks']}


def test_all_checks_equal_on_synthetic_data(collected):
    mode, report = collected
    found = statuses(report)
    ids = {check_id for check_id, _, _ in found}
    assert ids == set(qcChecks.CHECKS)                       # every check is reported
    for (check_id, sample, subject), status in found.items():
        if mode == 'se' and check_id == 'C1':
            assert status == 'not applicable'
        else:
            assert status == 'equal', f'{check_id} {sample} {subject}: {status}'
    assert sum(1 for c in report['checks'] if c['id'] == 'C3') == len(report['samples']) * 10   # 11 tiers, 10 transitions per sample
    assert sum(1 for c in report['checks'] if c['id'] == 'C9') == len(qcCollect.SAMPLE_TABLES)


@pytest.mark.parametrize('mode', ['pe', 'se'])
def test_spike_checks_not_applicable_without_spike_counting(tmp_path, mode):
    run, config_file = prepare_run(tmp_path, mode)
    shutil.rmtree(run / 'ribospike')                       # count_ribospike: False makes no spike-in tables
    config = yaml.safe_load(config_file.read_text())
    config['count_ribospike'] = False
    config_file.write_text(yaml.safe_dump(config))
    report = run_collector(run, config_file)
    found = statuses(report)
    assert 'not evaluable' not in found.values() and 'differs' not in found.values()
    off = [c for c in report['checks'] if c['status'] == 'not applicable' and c['id'] in ('C9', 'C10')]
    assert len(off) == 3 + len(report['samples'])           # three ribospike tables, C10 per sample
    assert all(qcChecks.SPIKE_OFF in c['subject'] for c in off)
    assert not any('ribospike' in w['source'] for w in report['warnings'])
    assert not any(record['id'].startswith('spike.')
                   for data in report['samples'].values() for record in data['metrics'])


@pytest.mark.parametrize('setting', ['remove_multimap', 'remove_antisense'])
def test_gene_checks_not_applicable_without_gene_table(tmp_path, setting):
    run, config_file = prepare_run(tmp_path, 'pe')
    for path in run.rglob('*-gene*.csv'):                  # the Snakefile does not make the gene table then
        path.unlink()
    config = yaml.safe_load(config_file.read_text())
    config[setting] = False
    config_file.write_text(yaml.safe_dump(config))
    report = run_collector(run, config_file)
    found = statuses(report)
    assert 'not evaluable' not in found.values()
    off = [c for c in report['checks'] if c['status'] == 'not applicable' and c['id'] in ('C7', 'C9')]
    assert len(off) == len(report['samples']) + 1           # C7 per sample, C9 for genecounts
    assert all(qcChecks.GENES_OFF in c['subject'] for c in off)
    assert not any('genecounts' in w['source'] for w in report['warnings'])
    assert not any(record['id'].startswith('counts.')
                   for data in report['samples'].values() for record in data['metrics'])
    antisense_warning = any(w['source'] == 'config: remove_antisense' for w in report['warnings'])
    assert antisense_warning == (setting == 'remove_antisense')


def test_missing_spike_tables_are_not_evaluable_when_spike_counting_is_on(tmp_path):
    run, config_file = prepare_run(tmp_path, 'pe')
    shutil.rmtree(run / 'ribospike')
    found = statuses(run_collector(run, config_file))
    assert found[('C10', 'S1', None)] == 'not evaluable'


# ----------------------------------------------------------------------------- one change per check

def edit_json(path, change):
    data = json.loads(path.read_text())
    change(data)
    path.write_text(json.dumps(data))


def add_one_to_table_cell(path, row_name, column):
    """Add 1 to one cell of a CSV table (row by its first field, column by header name)."""
    lines = path.read_text().splitlines()
    header = lines[0].split(',')
    for i, line in enumerate(lines[1:], start=1):
        fields = line.split(',')
        if fields[0] == row_name:
            index = header.index(column)
            fields[index] = str(int(float(fields[index])) + 1)
            lines[i] = ','.join(fields)
            path.write_text('\n'.join(lines) + '\n')
            return
    raise AssertionError(f'row {row_name} not found in {path.name}')


def first_gene(path):
    return path.read_text().splitlines()[1].split(',')[0]


def one_more_read_in_tier_log(path):
    """Bowtie 2 log: one more read in and one more unaligned (the log stays self-consistent)."""
    text = path.read_text()
    text = re.sub(r'^(\d+) reads; of these:', lambda m: f'{int(m.group(1)) + 1} reads; of these:', text, flags=re.M)
    text = re.sub(r'^(\s+)(\d+)( \([\d.]+%\) aligned 0 times)', lambda m: f'{m.group(1)}{int(m.group(2)) + 1}{m.group(3)}',
                  text, flags=re.M)
    path.write_text(text)


def add_sample_row(path):
    path.write_text(path.read_text() + 'S9,1,1,1,1,1\n')


def change(run, mode):
    """{check id: function that breaks that check's identity in the run folder}"""
    results = run / 'results'
    table = lambda name: results / f'{run.name}-{name}.csv'
    def written():                          # reads written by trimming: fastp (pe) or cutadapt (se)
        if mode == 'pe':
            edit_json(run / 'trimmed/logs/S1.fastp.json',
                      lambda d: d['merged_and_filtered'].__setitem__('total_reads', d['merged_and_filtered']['total_reads'] + 1))
        else:
            edit_json(run / 'trimmed/logs/S1.trim.cutadapt.json',
                      lambda d: d['read_counts'].__setitem__('output', d['read_counts']['output'] + 1))

    return {
        'C1': lambda: edit_json(run / 'trimmed/logs/S1.fastp.json',
                                lambda d: d['read1_before_filtering'].__setitem__('total_reads', d['read1_before_filtering']['total_reads'] + 1)),
        'C2': written,
        'C3': lambda: one_more_read_in_tier_log(run / 'alignments/logs/S1.human-mirna.log'),
        'C4': lambda: add_one_to_table_cell(table('mappingstats'), 'S1', 'Unmapped'),
        'C5': lambda: add_one_to_table_cell(table('mappingstats'), 'S1', 'Mapped'),
        'C6': lambda: add_one_to_table_cell(run / 'metrics/S1.log', 'reads_total', 'S1'),
        'C7': lambda: add_one_to_table_cell(table('genecounts'), first_gene(table('genecounts')), 'S1'),
        'C9': lambda: add_sample_row(table('mappingstats')),
        'C10': lambda: add_one_to_table_cell(table('mappingstats'), 'S1', 'Unique'),
    }


CASES = [(mode, check_id) for mode in ('pe', 'se') for check_id in qcChecks.CHECKS
         if check_id != 'C8' and not (mode == 'se' and check_id == 'C1')]


@pytest.mark.parametrize('mode, check_id', CASES, ids=[f'{m}-{c}' for m, c in CASES])
def test_changed_input_makes_check_differ(tmp_path, mode, check_id):
    run, config_file = prepare_run(tmp_path, mode)
    change(run, mode)[check_id]()
    report = run_collector(run, config_file)
    assert any(c['id'] == check_id and c['status'] == 'differs' for c in report['checks']), \
        [c for c in report['checks'] if c['id'] == check_id]


def test_wrong_species_makes_c8_differ_with_warning(tmp_path):
    run, config_file = prepare_run(tmp_path, 'pe')
    config = yaml.safe_load(config_file.read_text())
    config['species'] = 'mouse'                      # human BED, mouse run
    (tmp_path / 'references').mkdir()
    (tmp_path / 'references' / 'mouse.bed').write_text(
        (REPO / 'test/synthetic_qc/bundle/references/human.bed').read_text())
    config['dbdir'] = str(tmp_path / 'references') + '/'
    config_file.write_text(yaml.safe_dump(config))
    report = run_collector(run, config_file)
    c8 = [c for c in report['checks'] if c['id'] == 'C8'][0]
    assert c8['status'] == 'differs' and c8['observed'] > 0
    assert any('check C8' in w['message'] for w in report['warnings'])


def test_c8_does_not_apply_to_a_species_without_known_prefixes():
    assert qcChecks.species_check(['Rattus_norvegicus_tRNA-Ala-1'], 'rat')['status'] == 'not applicable'
    assert qcChecks.species_check(None, 'human')['status'] == 'not evaluable'


def test_missing_input_makes_check_not_evaluable(tmp_path):
    run, config_file = prepare_run(tmp_path, 'pe')
    (run / 'metrics' / 'S1.log').unlink()
    report = run_collector(run, config_file)
    assert statuses(report)[('C6', 'S1', None)] == 'not evaluable'
    assert statuses(report)[('C6', 'S2', None)] == 'equal'


# ----------------------------------------------------------------------------- Kind A

def test_shares_carry_n_and_flags(collected):
    _, report = collected
    records = {r['id']: r for r in report['samples']['S1']['metrics']}
    # S1 spike reads: 20 spikes x (7 exact C04 + 2 at -1 nt C26 + 1 at -2 nt C27 + 1 fragment C28)
    # + 4 x 3 pieces (C29) = 232; full-length reads of SR_20_4N-T1: 7 + 2 + 1 = 10
    spike_type = records['spike.20_T1.reads']                  # 10 of 232 spike reads: n >= 100, no flag
    assert spike_type['n'] == records['spike.reads']['value'] == 232 and 'flag' not in spike_type
    assert spike_type['share'] == 10 / 232
    modal = records['spike.20.modal_offset_reads']              # 28 at offset 0 of 40 full-length 20-nt reads
    assert modal['n'] == 40 and modal['share'] == 28 / 40 and modal['flag'] == 'low n'
    total = records['map.univec_bt1.aligned']                   # denominator map.total_in > 100
    assert total['n'] == records['map.total_in']['value'] and 'flag' not in total
    s2 = {r['id']: r for r in report['samples']['S2']['metrics']}
    assert s2['spike.fragment_reads']['flag'] == 'not evaluable' and 'share' not in s2['spike.fragment_reads']



def test_distributions_carry_n(collected):
    _, report = collected
    for record in report['samples']['S1']['distributions']:
        assert record['n'] == sum(record['y'])


@pytest.mark.parametrize('n, flag', [(0, 'not evaluable'), (1, 'low n'), (99, 'low n'), (100, None)])
def test_n_flag(n, flag):
    assert qcChecks.n_flag(n) == flag
