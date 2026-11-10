"""The passes file and the schedule with passes (docker/scripts/mapSchedule.py; release 0.4.0). The expected
schedules were written in the specification before the code.
Run from the repository root: python -m pytest docker/tests/references"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'docker' / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mapSchedule  # noqa: E402
from test_map_schedule import B, B_LAST, DEFAULT_NAMES, HUMAN_DBPATHS, PRINTOUT, R, U, write_all  # noqa: E402

DEMO = Path(__file__).resolve().parent / 'data'    # the test inputs, kept with the tests
HEADER = 'pass\tapplies_to\tmax_mismatches\tmin_length\tmax_length\tn_overhang\n'

# the kinds of options
U0 = '-S -k 1 -v 0 --chunkmbs 256'
B0N1 = '--sam-no-qname-trunc --no-unal --np 1 --very-sensitive --ignore-quals -k 100 --score-min C,-1'
B1 = '--sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals -k 100 --score-min C,-6'

PASS_TEST_PRINTOUT = """\
Mapping schedule: 14 references, 2 passes; order: default (references folder); pass settings: passes file test_passes.tsv
  step  reference     aligner  pass 1               pass 2
     1  univec_bt1    bowtie   0 mismatches         standard
     2  ribospike     bowtie2  spike-in settings    -
     3  rrna_bt1      bowtie   0 mismatches         standard
     4  trnas         bowtie2  0 mismatches, <=1 N  1 mismatch
     5  mirna         bowtie2  0 mismatches, <=1 N  1 mismatch
     6  star          bowtie2  0 mismatches, <=1 N  1 mismatch
     7  pre           bowtie2  0 mismatches, <=1 N  1 mismatch
     8  pirna         bowtie2  0 mismatches, <=1 N  1 mismatch
     9  ncrna_mt      bowtie2  0 mismatches, <=1 N  1 mismatch
    10  ncrna_snsno   bowtie2  0 mismatches, <=1 N  1 mismatch
    11  ncrna_other   bowtie2  0 mismatches, <=1 N  1 mismatch
    12  ncrna_lnc     bowtie2  0 mismatches, <=1 N  standard, 20+ nt
    13  cdna          bowtie2  0 mismatches, <=1 N  standard, 20+ nt
    14  genome        bowtie2  0 mismatches, <=1 N  standard, 20+ nt, keeps unaligned"""


@pytest.fixture
def refs(tmp_path):
    (tmp_path / 'human.dbpaths.csv').write_text(HUMAN_DBPATHS)
    return tmp_path


def resolve(refs, passes, order=None):
    return mapSchedule.resolve(str(refs / 'human.dbpaths.csv'), 'human', 100, order, True, passes)


def passes_file(tmp_path, *rows):
    path = tmp_path / 'passes.tsv'
    path.write_text(HEADER + ''.join('\t'.join(row) + '\n' for row in rows))
    return str(path)


def cells(rows):
    return [(r['pass'], r['step'], r['reference'], r['options'], r['min_length'], r['max_length'],
             r['max_mismatches'], r['n_overhang']) for r in rows]


# ----------------------------------------------------------------------------- §8.1 the identity passes file
def test_the_identity_passes_file_gives_the_default_schedule(refs, tmp_path):
    rows, left_out, problems, warnings = resolve(refs, str(DEMO / 'test_passes_identity.tsv'))
    assert (left_out, problems, warnings) == ([], [], [])
    default = write_all(resolve(refs, None)[0], refs, tmp_path / 'a')
    given = write_all(rows, refs, tmp_path / 'b')
    for name in ('schedule.tsv', 'order.tsv', 'tiers.csv'):
        assert (given / name).read_bytes() == (default / name).read_bytes(), name
    first, *rest = mapSchedule.render(rows, None, [], str(DEMO / 'test_passes_identity.tsv')).split('\n')
    assert first.endswith('pass settings: passes file test_passes_identity.tsv')
    assert rest == PRINTOUT.split('\n')[1:]


# ----------------------------------------------------------------------------- §8.2 the pass test file
def test_the_pass_test_file_gives_the_expected_27_cells(refs):
    rows, left_out, problems, warnings = resolve(refs, str(DEMO / 'test_passes.tsv'))
    assert (left_out, problems, warnings) == ([], [], [])
    bt2 = DEFAULT_NAMES[3:]                     # trnas ... genome
    expected = [('1', '1', 'univec_bt1', U0, '-', '-', '0', 'standard'),
                ('1', '2', 'ribospike', R, '-', '-', 'spike-in', 'spike-in'),
                ('1', '3', 'rrna_bt1', U0, '-', '-', '0', 'standard')]
    expected += [('1', str(4 + i), n, B0N1, '-', '-', '0', '1') for i, n in enumerate(bt2)]
    expected += [('2', '1', 'univec_bt1', U, '-', '-', 'standard', 'standard'),
                 ('2', '2', 'rrna_bt1', U, '-', '-', 'standard', 'standard')]
    expected += [('2', str(3 + i), n, B1, '-', '-', '1', 'standard') for i, n in enumerate(DEFAULT_NAMES[3:11])]
    expected += [('2', '11', 'ncrna_lnc', B, '20', '-', 'standard', 'standard'),
                 ('2', '12', 'cdna', B, '20', '-', 'standard', 'standard'),
                 ('2', '13', 'genome', B_LAST, '20', '-', 'standard', 'standard')]
    assert cells(rows) == expected


def test_the_pass_test_printout_is_the_expected_one(refs):
    rows = resolve(refs, str(DEMO / 'test_passes.tsv'))[0]
    assert mapSchedule.render(rows, None, [], str(DEMO / 'test_passes.tsv')) == PASS_TEST_PRINTOUT


def test_order_and_tiers_tables_list_each_reference_once(refs, tmp_path):
    rows = resolve(refs, str(DEMO / 'test_passes.tsv'))[0]
    out = write_all(rows, refs, tmp_path / 't')
    assert (out / 'order.tsv').read_text() == 'tier\n' + ''.join(f'{n}\n' for n in DEFAULT_NAMES)
    assert (out / 'tiers.csv').read_text() == HUMAN_DBPATHS


# ----------------------------------------------------------------------------- §7 from a row to a cell
@pytest.mark.parametrize('map_cmd, mm, n, expected', [
    ('bowtie2', '0', 'standard', B + ' --score-min C,0'),
    ('bowtie2', '1', 'standard', B1),
    ('bowtie2', '2', 'standard', B + ' --score-min C,-12'),
    ('bowtie2', '0', '2', B.replace('--np 5', '--np 1') + ' --score-min C,-2'),
    ('bowtie', '2', 'standard', '-S -k 1 -v 2 --chunkmbs 256'),
    ('bowtie1', '1', 'standard', '--fullref --best --strata -v 1 --chunkmbs 256 -m 100 -S'),    # no --no-unal (Bowtie 1.0.0)
    ('bowtie2', 'standard', 'standard', B),
])
def test_cell_options(map_cmd, mm, n, expected):
    assert ' '.join(mapSchedule.cell_options(map_cmd, mm, n, 100)) == expected


# ----------------------------------------------------------------------------- §6 the checks
@pytest.mark.parametrize('rows, problem', [
    ([('1', 'all', '3', '-', '-', '-')], 'max_mismatches'),
    ([('1', 'all', 'standard', '-', '-', '1')], 'needs max_mismatches 0'),
    ([('1', 'all', '0', '14', '-', '-')], 'min_length'),
    ([('1', 'all', '0', '30', '20', '-')], 'max_length is below min_length'),
    ([('0', 'all', '0', '-', '-', '-')], 'pass'),
    ([('1', 'all', '0', '-', '-', '-'), ('3', 'all', '0', '-', '-', '-')], 'consecutive from 1'),
    ([('1', 'all', '0', '-', '-', '-'), ('1', 'cdna', '1', '-', '-', '-')], 'cdna is already matched in pass 1'),
    ([('1', 'bowtie2', '0', '-', '-', '-'), ('1', 'trnas,mirna', '1', '-', '-', '-')], 'already matched'),
    ([('1', 'ribospike', '0', '-', '-', '-')], 'spike-in reference'),
    ([('1', 'mirtrons', '0', '-', '-', '-')], 'not in the order in effect'),
    ([('1', 'all,cdna', '0', '-', '-', '-')], 'must stand alone'),
    ([('1', 'other', '0', '-', '-', '-'), ('1', 'other', '1', '-', '-', '-')], 'two rows with other'),
])
def test_passes_files_are_checked(refs, tmp_path, rows, problem):
    out, _, problems, _ = resolve(refs, passes_file(tmp_path, *rows))
    assert out == [] and any(problem in p for p in problems), problems


@pytest.mark.parametrize('text, problem', [
    ('pass\tapplies_to\n1\tall\n', 'the first line must be the columns'),
    (HEADER, 'has no rows'),
    (HEADER + '1\tall\t0\n', '3 columns instead of 6'),
])
def test_unreadable_passes_files_are_refused(refs, tmp_path, text, problem):
    (tmp_path / 'p.tsv').write_text(text)
    _, _, problems, _ = resolve(refs, str(tmp_path / 'p.tsv'))
    assert len(problems) == 1 and problem in problems[0]


def test_a_name_left_out_of_the_order_is_refused(refs, tmp_path):
    order = tmp_path / 'order.tsv'
    order.write_text('tier\n' + '\n'.join(n for n in DEFAULT_NAMES if n != 'cdna') + '\n')
    _, _, problems, _ = resolve(refs, passes_file(tmp_path, ('1', 'cdna', '0', '-', '-', '-'),
                                                  ('1', 'other', '-', '-', '-', '-')), str(order))
    assert any('cdna is not in the order in effect' in p for p in problems)


@pytest.mark.parametrize('rows, warning', [
    ([('1', 'bowtie2', '0', '-', '-', '-')], 'univec_bt1 has no row in pass 1: vector reads are not removed first'),
    ([('1', 'bowtie2', '0', '-', '-', '-')], 'univec_bt1, rrna_bt1 run(s) in no pass'),
    ([('1', 'all', '2', '-', '-', '-')], 'line 2: max_mismatches 2 allows more than the standard settings for '
                                         + ', '.join(DEFAULT_NAMES[:1] + DEFAULT_NAMES[2:])),
    ([('1', 'all', '2', '19', '-', '-')], 'line 2: max_mismatches 2 allows more than the standard settings for '
                                          'univec_bt1, rrna_bt1'),          # Bowtie 1 always; Bowtie 2 below 19 nt
])
def test_passes_files_give_warnings(refs, tmp_path, rows, warning):
    out, _, problems, warnings = resolve(refs, passes_file(tmp_path, *rows))
    assert problems == [] and out
    assert f'passes.tsv: {warning}' in warnings, warnings


def test_one_mismatch_never_warns_and_two_warns_for_bowtie2_only_below_19_nt(refs, tmp_path):
    for min_length, warned in (('-', True), ('18', True), ('19', False)):
        _, _, _, warnings = resolve(refs, passes_file(tmp_path, ('1', 'bowtie2', '2', min_length, '-', '-'),
                                                      ('1', 'bowtie1', '1', '-', '-', '-')))
        assert any('allows more' in w for w in warnings) == warned, (min_length, warnings)


def test_a_reference_without_a_row_in_a_pass_has_no_cell_there(refs, tmp_path):
    rows, _, problems, _ = resolve(refs, passes_file(tmp_path, ('1', 'all', '0', '-', '-', '-'),
                                                     ('2', 'cdna,genome', '-', '-', '-', '-')))
    assert problems == []
    assert [(r['pass'], r['reference']) for r in rows if r['pass'] == '2'] == [('2', 'cdna'), ('2', 'genome')]
    assert '--no-unal' not in rows[-1]['options'] and all('--no-unal' in r['options'] for r in rows[3:-1])
    printout = mapSchedule.render(rows, None, [], 'p.tsv').split('\n')
    assert printout[3].split()[-1] == '-'                          # rrna_bt1 has no cell in pass 2


# ----------------------------------------------------------------------------- the order in effect with a partial pass 1
def test_a_partial_first_pass_keeps_the_order_in_effect(refs, tmp_path):
    """Pass 1 covers two references only: order.tsv, tiers.csv and the printout still follow the order in effect
    (the folder's order), not the order of the first cells (review of 2026-11-05), and order.tsv can be given
    back as an order file."""
    passes = passes_file(tmp_path, ('1', 'mirna,pirna', '0', '-', '-', '-'), ('2', 'all', 'standard', '-', '-', '-'))
    rows, _, problems, _ = resolve(refs, passes)
    assert problems == [] and [r['reference'] for r in rows[:3]] == ['ribospike', 'mirna', 'pirna']   # execution
    names = mapSchedule.order_in_effect(str(refs / 'human.dbpaths.csv'), 'human')
    out = tmp_path / 'tables'
    out.mkdir()
    mapSchedule.write_order(rows, out / 'order.tsv', names)
    mapSchedule.write_tiers(rows, str(refs / 'human.dbpaths.csv'), 'human', out / 'tiers.csv', names)
    assert (out / 'order.tsv').read_text() == 'tier\n' + ''.join(f'{n}\n' for n in DEFAULT_NAMES)
    assert (out / 'tiers.csv').read_text() == HUMAN_DBPATHS                     # byte copy: the default order
    printout = mapSchedule.render(rows, None, [], passes, names).split('\n')
    assert [line.split()[1] for line in printout[2:]] == DEFAULT_NAMES
    assert printout[2].split()[-2:] == ['-', 'standard'] and printout[5].split()[-2:] == ['-', 'standard']  # univec, trnas
    again, _, problems, _ = resolve(refs, passes, str(out / 'order.tsv'))       # order.tsv re-used as order file
    assert problems == [] and again == rows


def test_the_order_in_effect_of_an_order_file(refs, tmp_path):
    order = tmp_path / 'o.tsv'
    order.write_text('tier\n' + ''.join(f'{n}\n' for n in DEFAULT_NAMES[::-1][1:][::-1]))       # genome left out
    assert mapSchedule.order_in_effect(str(refs / 'human.dbpaths.csv'), 'human', str(order)) == DEFAULT_NAMES[:-1]
    rows = [{'reference': n} for n in ('mirna', 'univec_bt1', 'mirna')]
    assert mapSchedule.in_order(rows, DEFAULT_NAMES) == ['univec_bt1', 'mirna']
    assert mapSchedule.in_order(rows) == ['mirna', 'univec_bt1']                # without the order: first cells
