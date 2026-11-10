"""Collector rules:
ids come from the dictionary, mode-specific metrics are absent in the other mode,
a missing or damaged input gives a warning and omitted metrics (never a changed
number or a zero), and the output is identical between runs."""
import csv
import json
import subprocess
import sys

import pytest

import qcCollect
from conftest import DICTIONARY, FIXTURES, REPO, prepare_run, run_collector, template_pattern, values

# Dictionary ids the synthetic fixture cannot produce, with the reason.
NOT_PRODUCED = {
    'prov.reference_build': 'the fixture has no reference MANIFEST (expected warning)',
    'counts.rbc_mirna_ratio': 'the fixture has no miR-23a-3p reads, so the ratio is not evaluable (omitted)',
    # release 0.4.0: cell metrics exist only in runs with --passes (the fixtures are default runs; see
    # test_parsers.py for a run with passes made from the fixture)
    'map.p{pass}.{tier}.arrived': 'only with --passes',
    'map.p{pass}.{tier}.aligned': 'only with --passes',
    'map.p{pass}.{tier}.unaligned': 'only with --passes',
    'map.p{pass}.{tier}.skipped_window': 'only with --passes',
}


def produced_ids(report):
    ids = set(report['provenance'])
    for data in report['samples'].values():
        ids.update(record['id'] for record in data['metrics'])
        ids.update(record['id'] for record in data['distributions'])
    return ids


def test_every_dictionary_id_of_the_mode_is_produced(collected):
    mode, report = collected
    ids = produced_ids(report)
    missing = [entry['id'] for entry in DICTIONARY['metrics']
               if mode in entry['modes'] and entry['id'] not in NOT_PRODUCED
               and not any(template_pattern(entry['id']).match(metric_id) for metric_id in ids)]
    assert not missing, f'dictionary ids not produced in mode {mode}: {missing}'


def test_other_mode_metrics_are_absent(collected):
    mode, report = collected
    ids = produced_ids(report)
    other_mode_only = [entry['id'] for entry in DICTIONARY['metrics'] if mode not in entry['modes']]
    present = [metric_id for metric_id in ids
               for template in other_mode_only if template_pattern(template).match(metric_id)]
    assert not present, f'metrics of the other mode present in mode {mode}: {present}'


def test_output_is_identical_between_runs(tmp_path):
    for mode in ('pe', 'se'):
        outputs = []
        for attempt in ('first', 'second'):
            work = tmp_path / f'{mode}_{attempt}' / 'work'      # same folder name: paths are reduced to it
            work.mkdir(parents=True)
            run, config_file = prepare_run(work, mode)
            run_collector(run, config_file)
            outputs.append([(run / 'qc' / name).read_bytes() for name in ('metrics_long.tsv', 'consistency_checks.tsv', 'qc_metrics.json')])
        assert outputs[0] == outputs[1]


# ----------------------------------------------------------------------------- damaged inputs

def delete(path):
    path.unlink()


def cut_in_half(path):
    text = path.read_text()
    path.write_text(text[:len(text) // 2])


def garble_second_line(path):
    lines = path.read_text().splitlines(keepends=True)
    lines[1] = 'garbage\n'
    path.write_text(''.join(lines))


def drop_unaligned_line(path):
    path.write_text(''.join(line for line in path.read_text().splitlines(keepends=True)
                            if 'aligned 0 times' not in line and 'failed to align' not in line))


# (file relative to the run folder, damage, modes); {run} = run folder name
DAMAGE = [
    ('trimmed/logs/S1.trim.cutadapt.json', delete, 'pe se'),
    ('trimmed/logs/S1.trim.cutadapt.json', cut_in_half, 'pe se'),
    ('trimmed/logs/S1.dimer.cutadapt.json', delete, 'pe se'),
    ('trimmed/logs/S1.fastp.json', delete, 'pe'),
    ('trimmed/logs/S1.fastp.json', cut_in_half, 'pe'),
    ('logs/map.log', delete, 'pe se'),
    ('alignments/logs/S1.human-cdna.log', delete, 'pe se'),
    ('alignments/logs/S1.human-cdna.log', drop_unaligned_line, 'pe se'),
    ('alignments/logs/S1.univec_bt1.log', drop_unaligned_line, 'pe se'),
    ('metrics/S1.log', delete, 'pe se'),
    ('results/{run}-mappingstats.csv', garble_second_line, 'pe se'),
    ('results/{run}-classcounts.csv', delete, 'pe se'),
    ('results/{run}-classcounts.csv', garble_second_line, 'pe se'),
    ('results/{run}-classlengths.csv', cut_in_half, 'pe se'),
    ('results/{run}-classdistinctseqs.csv', delete, 'pe se'),
    ('results/{run}-genecounts.csv', garble_second_line, 'pe se'),
    ('results/{run}-genelengths.csv', delete, 'pe se'),
    ('ribospike/results/{run}-readorigin.csv', delete, 'pe se'),
    ('ribospike/results/{run}-ribospike.csv', garble_second_line, 'pe se'),
    ('ribospike/results/{run}-spikelengths.csv', cut_in_half, 'pe se'),
    ('ribospike/results/{run}-ribospikeends.csv', garble_second_line, 'pe se'),
    ('ribospike/results/{run}-ribospikeends_exact.csv', delete, 'pe se'),
]
CASES = [(mode, name, damage) for name, damage, modes in DAMAGE for mode in modes.split()]


@pytest.mark.parametrize('mode, name, damage', CASES,
                         ids=[f'{mode}-{name.split("/")[-1]}-{damage.__name__}' for mode, name, damage in CASES])
def test_damaged_input_gives_warning_and_omits_metrics(tmp_path, baselines, mode, name, damage):
    baseline = baselines[mode]
    run, config_file = prepare_run(tmp_path, mode)
    damage(run / name.format(run=run.name))
    report = run_collector(run, config_file)

    expected, found = values(baseline), values(report)
    changed = {key: (expected.get(key), value) for key, value in found.items() if expected.get(key) != value}
    assert not changed, f'values that differ from the undamaged run: {changed}'
    assert len(found) < len(expected), 'no metric was omitted'
    new_warnings = [w for w in report['warnings'] if w not in baseline['warnings']]
    assert new_warnings, 'no warning'


@pytest.mark.parametrize('reads451, reads23a', [(30, 12), (5, 0)])
def test_rbc_mirna_ratio_is_a_plain_ratio(tmp_path, reads451, reads23a):
    """miR-451a / miR-23a-3p reads, without pseudocount; omitted when miR-23a-3p has no reads."""
    run, config_file = prepare_run(tmp_path, 'pe')
    for name, rows in (('genecounts', [f'Hsa-Mir-451_5p,{reads451},0,0,0', f'Hsa-Mir-23-P3_3p,{reads23a},0,0,0']),   # S1-S4
                       ('genelengths', ['Hsa-Mir-451_5p,22.0,miRNA', 'Hsa-Mir-23-P3_3p,21.0,miRNA'])):
        table = next((run / 'results').glob(f'*-{name}.csv'))
        table.write_text(table.read_text() + '\n'.join(rows) + '\n')
    found = values(run_collector(run, config_file))
    if reads23a:
        assert found[('S1', 'counts.rbc_mirna_ratio')] == reads451 / reads23a
    else:
        assert ('S1', 'counts.rbc_mirna_ratio') not in found


# ----------------------------------------------------------------------------- metrics read from one field

def test_single_field_metrics_equal_their_source(collected):
    """Metrics that copy one field of a tool report (or divide two fields), compared with that field in
    the fixture; the synthetic truth does not define these numbers."""
    mode, report = collected
    found = values(report)
    run = FIXTURES / f'levelA_{mode}'
    for sample in report['samples']:
        trim = json.loads((run / 'trimmed' / 'logs' / f'{sample}.trim.cutadapt.json').read_text())
        bases = trim['adapters_read1'][0]['three_prime_end']['adjacent_bases']
        expected = {'trim.r1_base_before_adapter': (['A', 'C', 'G', 'T', 'none'],
                                                    [bases[base] for base in ('A', 'C', 'G', 'T', '')])}
        if mode == 'pe':
            fastp = json.loads((run / 'trimmed' / 'logs' / f'{sample}.fastp.json').read_text())
            merged, filtering = fastp['merged_and_filtered'], fastp['filtering_result']
            histogram = fastp['insert_size']['histogram']
            expected.update({
                'trim.r1_with_adapter': trim['read_counts']['read1_with_adapter'],
                'trim.r2_with_adapter': trim['read_counts']['read2_with_adapter'],
                'trim.fastp_too_short_reads': filtering['too_short_reads'],
                'trim.fastp_corrected_reads': filtering['corrected_reads'],
                'trim.merged_q30_rate': merged['q30_bases'] / merged['total_bases'],
                'trim.insert_size_peak': fastp['insert_size']['peak'],
                'trim.insert_size_hist': (list(range(len(histogram))), histogram)})
        else:
            expected['trim.reads_with_adapter'] = trim['read_counts']['read1_with_adapter']
        with open(run / 'metrics' / f'{sample}.log', newline='') as f:
            metrics_log = {row[0]: row[1] for row in csv.reader(f) if len(row) >= 2}
        expected['assign.median_read_length'] = float(metrics_log['median_read_length'])
        for metric_id, value in expected.items():
            assert found[(sample, metric_id)] == value, f'{sample} {metric_id}'


def test_floor_reads_are_the_reads_of_15_to_17_nt(collected):
    """assign.class.{rna_class}.floor_reads recomputed from the fixture's length table; the exploratory
    sample S3 has the 16- and 17-nt reads (C32-C35) that the truth samples lack."""
    mode, report = collected
    found = values(report)
    with open(FIXTURES / f'levelA_{mode}' / 'results' / f'levelA_{mode}-classlengths.csv', newline='') as f:
        rows = list(csv.reader(f))
    lengths = [int(length) for length in rows[0][2:]]
    for sample, biotype, *counts in rows[1:]:
        if biotype != 'unmapped':
            expected = sum(int(float(n)) for length, n in zip(lengths, counts) if 15 <= length <= 17)
            assert found[(sample, f'assign.class.{biotype}.floor_reads')] == expected, f'{sample} {biotype}'
    assert any(found[(s, f'assign.class.{b}.floor_reads')] for s, b, *_ in rows[1:] if s == 'S3' and b != 'unmapped')


# ----------------------------------------------------------------------------- register A4, warnings, layout

def is_number(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def zero_sample(run, sample):
    """Make one sample look as if no read entered mapping: tier logs of a tier without reads, and 0 for
    every number of that sample in the result tables (the trimming reports are left as they are)."""
    for log in (run / 'alignments' / 'logs').glob(f'{sample}.*.log'):
        log.write_text('0 reads\n0.00% overall alignment rate\n')
    tables = [*(run / 'results').glob('*.csv'), *(run / 'ribospike' / 'results').glob('*.csv'),
              run / 'metrics' / f'{sample}.log']
    for table in tables:
        with open(table, newline='') as f:
            rows = list(csv.reader(f))
        if sample in rows[0]:                            # one column per sample
            column = rows[0].index(sample)
            for row in rows[1:]:
                row[column] = '0'
        else:                                            # one row per sample (and item)
            for row in rows[1:]:
                if row and row[0] == sample:
                    row[:] = ['0' if is_number(cell) else cell for cell in row]
        with open(table, 'w', newline='') as f:
            csv.writer(f).writerows(rows)


def test_sample_without_reads_is_not_evaluable(tmp_path, baselines):
    """Register A4: a sample with 0 reads entering mapping gives a warning, its shares are 'not evaluable'
    (no share, never a 0 % that looks like a measurement), no medians are written, the other samples keep
    their values, and the report is still made."""
    run, config_file = prepare_run(tmp_path, 'pe')
    zero_sample(run, 'S1')
    report = run_collector(run, config_file)
    assert any(w['sample'] == 'S1' and w['source'] == 'map.total_in' for w in report['warnings'])
    records = {r['id']: r for r in report['samples']['S1']['metrics']}
    for record in records.values():
        if record.get('n') == 0:
            assert record.get('flag') == 'not evaluable' and 'share' not in record, record['id']
    assert records['map.final_unaligned'].get('flag') == 'not evaluable'
    assert not [metric_id for metric_id in records if metric_id.endswith('.median_length')]
    expected, found = values(baselines['pe']), values(report)
    assert all(found[key] == value for key, value in expected.items() if key[0] == 'S2')
    result = subprocess.run([sys.executable, str(REPO / 'docker' / 'scripts' / 'qcReport.py'), '--output', str(run)],
                            cwd=REPO, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('settings, sample, expected', [
    ({'paired_end': False}, 'S1_NEB', ['sample name']),
    ({'paired_end': True}, 'S1_SMART', ['sample name']),
    ({'paired_end': True}, 'S1_NEB', []),             # the NEB branch exists only in single-end trimming
    ({'paired_end': True, 'species': 'mouse'}, 'S1', ['config: species']),
    ({'paired_end': True, 'species': 'human'}, 'S1', []),
    ({'paired_end': True, 'remove_antisense': False}, 'S1', []),      # 0.4.0: the tagging gap is closed
    ({'paired_end': True, 'remove_multimap': False}, 'S1', []),
])
def test_unsupported_settings_give_a_warning(settings, sample, expected):
    """Settings and trimming branches the QC numbers are not designed for."""
    c = qcCollect.Collector({'metrics': []}, 'pe' if settings['paired_end'] else 'se', [sample])
    qcCollect.warn_unsupported(c, settings, [sample])
    assert [w['source'] for w in c.warnings] == expected


def test_result_tables_found_for_an_absolute_output_prefix(tmp_path):
    """The Snakefile writes <prefix>/<prefix>-<name>.csv, so the absolute prefix /a/run puts the tables at
    /a/run/a/run-<name>.csv."""
    run, config_file = prepare_run(tmp_path, 'pe')       # the config's output_prefix is absolute
    nested = run / run.as_posix().lstrip('/')
    nested.parent.mkdir(parents=True)
    for table in (run / 'results').glob('*.csv'):
        table.rename(nested.parent / f"{nested.name}-{table.name.split('-', 1)[1]}")
    report = run_collector(run, config_file)
    assert [w['source'] for w in report['warnings']] == ['<dbdir>/MANIFEST_human.json']
    assert 'assign.mapped' in {r['id'] for r in report['samples']['S1']['metrics']}


# ----------------------------------------------------------------------------- release 0.4.0: removal split, C11
# Level B pe sample S3 (a default run before the antisense change of 0.4.0), and the kept reads each flag state gives
MAPPED, MULTI, ANTISENSE, ANTISENSE_MULTI = 152_000, 11_826, 1_565, 767


@pytest.mark.parametrize('flags, removed_multimap, removed_antisense', [
    ((True, True), MULTI, ANTISENSE),
    ((True, False), MULTI, 0),
    ((False, True), 0, ANTISENSE + ANTISENSE_MULTI),      # antisense multimappers are removed as antisense
    ((False, False), 0, 0),
])
def test_removed_reads_split_by_the_flags_and_c11(flags, removed_multimap, removed_antisense):
    import qcChecks
    c = qcCollect.Collector(DICTIONARY, 'pe', ['S1'])
    for metric_id, value in (('assign.mapped', MAPPED), ('assign.multi', MULTI), ('assign.antisense', ANTISENSE),
                             ('assign.antisense_multi', ANTISENSE_MULTI),
                             ('assign.kept', MAPPED - removed_multimap - removed_antisense)):
        c.put('S1', metric_id, value)
    qcCollect.derive_removed_by_policy(c, 'S1', {'remove_multimap': flags[0], 'remove_antisense': flags[1]})
    assert (c.get('S1', 'assign.removed_multimap'), c.get('S1', 'assign.removed_antisense')) == \
           (removed_multimap, removed_antisense)
    values = {metric_id: record['value'] for metric_id, record in c.metrics['S1'].items()}
    c11 = [x for x in qcChecks.sample_checks('S1', values, {}, 'pe', [], True, True, set()) if x['id'] == 'C11']
    assert [x['status'] for x in c11] == ['equal']


def test_runs_without_antisense_columns_get_no_split_and_c11_not_evaluable():
    import qcChecks
    c = qcCollect.Collector(DICTIONARY, 'pe', ['S1'])
    for metric_id, value in (('assign.mapped', MAPPED), ('assign.multi', MULTI), ('assign.kept', MAPPED - MULTI)):
        c.put('S1', metric_id, value)
    qcCollect.derive_removed_by_policy(c, 'S1', {})
    assert c.get('S1', 'assign.removed_by_policy') == MULTI
    assert c.get('S1', 'assign.removed_multimap') is None and c.get('S1', 'assign.removed_antisense') is None
    values = {metric_id: record['value'] for metric_id, record in c.metrics['S1'].items()}
    c11 = [x for x in qcChecks.sample_checks('S1', values, {}, 'pe', [], True, True, set()) if x['id'] == 'C11']
    assert [x['status'] for x in c11] == ['not evaluable']
