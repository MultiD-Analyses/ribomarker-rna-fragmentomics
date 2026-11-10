"""The order file and the mapping schedule (docker/scripts/mapSchedule.py; release 0.4.0). The expected schedules
were written in the specification before the code.
Run from the repository root: python -m pytest docker/tests/references"""
import hashlib
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'docker' / 'scripts'))
import dbPaths      # noqa: E402
import mapSchedule  # noqa: E402

# docker/assets/human.dbpaths.csv (not in git; its sha256 is in MANIFEST_human.json, checked below)
HUMAN_DBPATHS = ('bowtie,univec/univec_bt1,T\nribospikebt2,ribospike/ribospike,T\nbowtie,rrna/human-rrna_bt1,T\n'
                 'bowtie2,trna/human-trnas,T\nbowtie2,mirna/human-mirna,T\nbowtie2,mirna/human-star,T\n'
                 'bowtie2,mirna/human-pre,T\nbowtie2,pirna/human-pirna,T\nbowtie2,gtf/human-ncrna_mt,T\n'
                 'bowtie2,gtf/human-ncrna_snsno,T\nbowtie2,gtf/human-ncrna_other,T\nbowtie2,gtf/human-ncrna_lnc,T\n'
                 'bowtie2,gtf/human-cdna,T\nbowtie2,genome/human-genome,T\n')
DEFAULT_NAMES = ['univec_bt1', 'ribospike', 'rrna_bt1', 'trnas', 'mirna', 'star', 'pre', 'pirna', 'ncrna_mt',
                 'ncrna_snsno', 'ncrna_other', 'ncrna_lnc', 'cdna', 'genome']
# Test inputs kept with the tests (data/), so that the tests run in every copy of the repository
ORDER_DEFAULT = Path(__file__).resolve().parent / 'data' / 'test_order_default.tsv'
# The order as in the RealSeq paper: ncrna_snsno before pirna
REALSEQ_NAMES = ['univec_bt1', 'ribospike', 'rrna_bt1', 'trnas', 'mirna', 'star', 'pre', 'ncrna_snsno', 'pirna',
                 'ncrna_mt', 'ncrna_other', 'ncrna_lnc', 'cdna', 'genome']

# the options of the four kinds of row (multimap 100)
U = '-S -k 1 -n 1 -l 18 -e 80 --chunkmbs 256'
R = ('--sam-no-qname-trunc --no-unal --no-1mm-upfront --mp 6 --rdg 6,3 --rfg 6,3 --norc --n-ceil L,5,0.15 -k 1 -D 20 '
     '-R 2 -N 0 -L 10 -i S,1,0.40 --score-min L,0,-0.6 --np 1 --ignore-quals')
B = '--sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals -k 100'
B_LAST = '--sam-no-qname-trunc --np 5 --very-sensitive --ignore-quals -k 100'
KIND = {'univec_bt1': U, 'ribospike': R, 'rrna_bt1': U, 'genome': B_LAST}

PRINTOUT = """Mapping schedule: 14 references, 1 pass; order: default (references folder); pass settings: none
  step  reference     aligner  pass 1
     1  univec_bt1    bowtie   standard
     2  ribospike     bowtie2  spike-in settings
     3  rrna_bt1      bowtie   standard
     4  trnas         bowtie2  standard
     5  mirna         bowtie2  standard
     6  star          bowtie2  standard
     7  pre           bowtie2  standard
     8  pirna         bowtie2  standard
     9  ncrna_mt      bowtie2  standard
    10  ncrna_snsno   bowtie2  standard
    11  ncrna_other   bowtie2  standard
    12  ncrna_lnc     bowtie2  standard
    13  cdna          bowtie2  standard
    14  genome        bowtie2  standard, keeps unaligned"""


@pytest.fixture
def refs(tmp_path):
    (tmp_path / 'human.dbpaths.csv').write_text(HUMAN_DBPATHS)
    return tmp_path


def resolve(refs, order=None, count_ribospike=True):
    return mapSchedule.resolve(str(refs / 'human.dbpaths.csv'), 'human', 100, order, count_ribospike)[:3]  # passes: + warnings


def expected_rows(names):
    rows = []
    for step, name in enumerate(names, 1):
        map_cmd = {'univec_bt1': 'bowtie', 'ribospike': 'ribospikebt2', 'rrna_bt1': 'bowtie'}.get(name, 'bowtie2')
        spike = name == 'ribospike'
        rows.append({'pass': '1', 'step': str(step), 'reference': name,
                     'index_name': name if name in ('univec_bt1', 'ribospike') else f'human-{name}',
                     'map_cmd': map_cmd, 'aligner': 'bowtie' if map_cmd == 'bowtie' else 'bowtie2',
                     'min_length': '-', 'max_length': '-',
                     'max_mismatches': 'spike-in' if spike else 'standard',
                     'n_overhang': 'spike-in' if spike else 'standard', 'options': KIND.get(name, B)})
    return rows


def write_all(rows, refs, out):
    out.mkdir()
    mapSchedule.write_schedule(rows, out / 'schedule.tsv')
    mapSchedule.write_order(rows, out / 'order.tsv')
    mapSchedule.write_tiers(rows, str(refs / 'human.dbpaths.csv'), 'human', out / 'tiers.csv')
    return out


def test_the_embedded_dbpaths_file_is_that_of_the_human_references():
    manifest = json.loads((REPO / 'docker' / 'assets' / 'MANIFEST_human.json').read_text())
    assert hashlib.sha256(HUMAN_DBPATHS.encode()).hexdigest() == manifest['outputs']['dbpaths']['sha256']


# ----------------------------------------------------------------------------- default
def test_the_default_schedule_is_the_expected_one(refs, tmp_path):
    rows, left_out, problems = resolve(refs)
    assert (problems, left_out) == ([], [])
    assert rows == expected_rows(DEFAULT_NAMES)
    out = write_all(rows, refs, tmp_path / 'tables')
    assert mapSchedule.read_schedule(out / 'schedule.tsv') == rows
    assert (out / 'tiers.csv').read_bytes() == HUMAN_DBPATHS.encode()          # a byte copy, as before 0.4.0
    assert (out / 'order.tsv').read_text() == 'tier\n' + ''.join(f'{n}\n' for n in DEFAULT_NAMES)


def test_the_default_printout_is_the_expected_one(refs):
    rows, _, _ = resolve(refs)
    assert mapSchedule.render(rows) == PRINTOUT


# ----------------------------------------------------------------------------- §4.2: an order file in the default order
def test_an_order_file_in_the_default_order_gives_the_same_files(refs, tmp_path):
    default = write_all(resolve(refs)[0], refs, tmp_path / 'a')
    rows, left_out, problems = resolve(refs, str(ORDER_DEFAULT))
    assert (problems, left_out) == ([], [])
    given = write_all(rows, refs, tmp_path / 'b')
    for name in ('schedule.tsv', 'order.tsv', 'tiers.csv'):
        assert (given / name).read_bytes() == (default / name).read_bytes(), name
    assert ORDER_DEFAULT.read_text() == (default / 'order.tsv').read_text()     # the file is a default order.tsv
    first, *rest = mapSchedule.render(rows, str(ORDER_DEFAULT)).split('\n')
    assert first == ('Mapping schedule: 14 references, 1 pass; order: order file test_order_default.tsv; '
                     'pass settings: none')
    assert rest == PRINTOUT.split('\n')[1:]


# ----------------------------------------------------------------------------- §4.3: the order test file
def test_the_order_test_file_gives_the_expected_schedule(refs, tmp_path):
    names = DEFAULT_NAMES[:7] + ['ncrna_snsno', 'pirna', 'ncrna_mt'] + DEFAULT_NAMES[10:]
    order_realseq = tmp_path / 'order_realseq.tsv'
    order_realseq.write_text('tier\n' + ''.join(f'{n}\n' for n in REALSEQ_NAMES))
    assert REALSEQ_NAMES == names
    rows, left_out, problems = resolve(refs, str(order_realseq))
    assert (problems, left_out) == ([], [])
    assert rows == expected_rows(names)
    out = write_all(rows, refs, tmp_path / 'tables')
    assert (out / 'tiers.csv').read_text() == (
        'bowtie,univec/univec_bt1,T\nribospikebt2,ribospike/ribospike,T\nbowtie,rrna/human-rrna_bt1,T\n'
        'bowtie2,trna/human-trnas,T\nbowtie2,mirna/human-mirna,T\nbowtie2,mirna/human-star,T\n'
        'bowtie2,mirna/human-pre,T\nbowtie2,gtf/human-ncrna_snsno,T\nbowtie2,pirna/human-pirna,T\n'
        'bowtie2,gtf/human-ncrna_mt,T\nbowtie2,gtf/human-ncrna_other,T\nbowtie2,gtf/human-ncrna_lnc,T\n'
        'bowtie2,gtf/human-cdna,T\nbowtie2,genome/human-genome,T\n')
    assert (out / 'order.tsv').read_text() == order_realseq.read_text()
    # tables/tiers.csv is read back as a tier list (contig_tiers.tsv, recount_bam.py): the same indexes in the new order
    assert [dbPaths.tier_id(i, 'human') for _, i, _ in dbPaths.read_dbpaths(str(out / 'tiers.csv'), str(refs))] == names
    lines = mapSchedule.render(rows, str(order_realseq)).split('\n')
    assert lines[0].endswith('order: order file order_realseq.tsv; pass settings: none')
    assert [line.split()[1] for line in lines[2:]] == names


# ----------------------------------------------------------------------------- the checks of the order file
def order_file(tmp_path, names, header='tier'):
    path = tmp_path / 'order.tsv'
    path.write_text('\n'.join(([header] if header else []) + names) + '\n')
    return str(path)


@pytest.mark.parametrize('names, problem', [
    (['ribospike', 'univec_bt1'] + DEFAULT_NAMES[2:], None),                       # both first, any order: allowed
    (['univec_bt1', 'rrna_bt1', 'ribospike'] + DEFAULT_NAMES[3:], 'univec_bt1 and ribospike must come first'),
    (['trnas', 'univec_bt1', 'ribospike'] + DEFAULT_NAMES[4:], 'must come first'),
    (DEFAULT_NAMES[:-2] + ['genome', 'cdna'], 'genome must come last'),
    (DEFAULT_NAMES + ['trnas'], 'trnas is listed twice'),
    (DEFAULT_NAMES + ['mirtrons'], 'mirtrons is not a reference of the references folder'),
    (['univec_bt1'] + DEFAULT_NAMES[2:], 'ribospike is needed while count_ribospike is on'),
    ([], 'names no reference'),
])
def test_order_files_are_checked(refs, tmp_path, names, problem):
    rows, _, problems = resolve(refs, order_file(tmp_path, names))
    if problem is None:
        assert problems == [] and [r['reference'] for r in rows] == names
    else:
        assert any(problem in p for p in problems), problems
        assert rows == []


def test_references_left_out_are_not_used_with_a_warning(refs, tmp_path):
    names = [n for n in DEFAULT_NAMES if n not in ('star', 'cdna')]
    rows, left_out, problems = resolve(refs, order_file(tmp_path, names))
    assert problems == [] and left_out == ['star', 'cdna']
    assert [r['reference'] for r in rows] == names
    assert mapSchedule.render(rows, 'x.tsv', left_out).endswith('\n  not used: star, cdna')


def test_the_spike_ins_may_be_left_out_when_they_are_not_counted(refs, tmp_path):
    names = ['univec_bt1'] + DEFAULT_NAMES[2:]
    rows, left_out, problems = resolve(refs, order_file(tmp_path, names), count_ribospike=False)
    assert problems == [] and left_out == ['ribospike']


def test_the_last_reference_keeps_the_unaligned_reads_whatever_it_is(refs, tmp_path):
    rows, _, _ = resolve(refs, order_file(tmp_path, DEFAULT_NAMES[:-1]))          # genome left out: cdna is last
    assert rows[-1]['reference'] == 'cdna' and rows[-1]['options'] == B_LAST
    assert all('--no-unal' in r['options'] for r in rows[1:-1] if r['aligner'] == 'bowtie2')


@pytest.mark.parametrize('text, problem', [
    ('bowtie,univec/univec_bt1,T\nbowtie2,genome/human-genome,T\n', 'tier list of release 0.3.0'),
    ('univec_bt1\ngenome\n', 'the first line must be "tier"'),
])
def test_the_old_format_and_a_missing_header_are_refused_with_the_new_format(refs, tmp_path, text, problem):
    (tmp_path / 'mine.tsv').write_text(text)
    _, _, problems = resolve(refs, str(tmp_path / 'mine.tsv'))
    assert len(problems) == 1 and problem in problems[0] and 'tables/order.tsv' in problems[0]


def test_an_order_file_with_windows_line_ends_and_a_byte_order_mark_is_read(refs, tmp_path):
    (tmp_path / 'mine.tsv').write_bytes(b'\xef\xbb\xbftier\r\n' + '\r\n'.join(DEFAULT_NAMES).encode() + b'\r\n')
    rows, _, problems = resolve(refs, str(tmp_path / 'mine.tsv'))
    assert problems == [] and [r['reference'] for r in rows] == DEFAULT_NAMES


def test_a_reference_not_kept_in_dbpaths_is_refused(tmp_path):
    (tmp_path / 'human.dbpaths.csv').write_text(HUMAN_DBPATHS.replace('human-cdna,T', 'human-cdna,F'))
    _, _, problems = resolve(tmp_path)
    assert len(problems) == 1 and 'human-cdna has retain F' in problems[0]


def test_an_unknown_aligner_in_dbpaths_is_refused(tmp_path):
    (tmp_path / 'human.dbpaths.csv').write_text(HUMAN_DBPATHS.replace('bowtie2,gtf/human-cdna', 'bwa,gtf/human-cdna'))
    _, _, problems = resolve(tmp_path)
    assert len(problems) == 1 and 'unknown aligner bwa' in problems[0]


def test_mouse_references_get_the_same_names(tmp_path):
    (tmp_path / 'mouse.dbpaths.csv').write_text(HUMAN_DBPATHS.replace('human-', 'mouse-'))
    rows, _, problems, _ = mapSchedule.resolve(str(tmp_path / 'mouse.dbpaths.csv'), 'mouse', 100)
    assert problems == [] and [r['reference'] for r in rows] == DEFAULT_NAMES
    assert rows[-1]['index_name'] == 'mouse-genome'


# ----------------------------------------------------------------------------- the commands mapReads runs
# mapReads.py's command templates before 0.4.0 (mapping_cmds, as committed before the ncRNA split), with the file arguments filled in
BEFORE = {
    'bowtie': 'zcat -f IN | bowtie -S -k 1 -n 1 -l 18 -e 80 --chunkmbs 256 -p 8 --un UNAL IDX - SAM 2> LOG',
    'bowtie1': ('bowtie --fullref --no-unal --best --strata -n 1 -l 18 -e 80 --chunkmbs 256 -m 100 -p 8 --un UNAL '
                '-S IDX IN SAM 2> LOG'),
    'mind': 'bowtie --fullref --no-unal -k 1 -v 0 --best --norc -p 8 --un UNAL -S IDX IN SAM 2> LOG',
    'ribospike': 'bowtie --fullref --no-unal -k 1 -n 2 --best --norc -p 8 --un UNAL -S IDX IN SAM 2> LOG',
    'ribospikebt2': (f'bowtie2 {R} --threads 8 --un-gz UNAL -x IDX -U IN -S SAM 2> LOG'),
    'bowtie2': ('bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads 8 --un-gz UNAL '
                '-x IDX -k 100 -U IN -S SAM 2> LOG'),
    'bowtie2best': ('bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads 8 '
                    '--un-gz UNAL -x IDX -U IN -S SAM 2> LOG'),
    'bowtie2sense': ('bowtie2 --sam-no-qname-trunc --norc --no-unal --np 5 --very-sensitive --ignore-quals --threads 8 '
                     '--un-gz UNAL -x IDX -U IN -S SAM 2> LOG'),
}


# Changed on purpose (2026-11-05): bowtie1, mind and ribospike without --no-unal (the pinned Bowtie 1.0.0 refuses it)
# and with their input from zcat -f, like `bowtie`; mapSchedule.py OPTIONS, IO
NOW = dict(BEFORE, **{
    'bowtie1': ('zcat -f IN | bowtie --fullref --best --strata -n 1 -l 18 -e 80 --chunkmbs 256 -m 100 -p 8 --un UNAL '
                '-S IDX - SAM 2> LOG'),
    'mind': 'zcat -f IN | bowtie --fullref -k 1 -v 0 --best --norc -p 8 --un UNAL -S IDX - SAM 2> LOG',
    'ribospike': 'zcat -f IN | bowtie --fullref -k 1 -n 2 --best --norc -p 8 --un UNAL -S IDX - SAM 2> LOG',
})


def pairs(command):
    """The command as (option, value) pairs and positional words, in a set: the order of the options does not
    matter to the aligners [inferred; checked read by read], the order of the positional words does."""
    words = command.split()
    flags_with_value = {'-k', '-n', '-l', '-e', '--chunkmbs', '-p', '--un', '-m', '-v', '--threads', '--un-gz', '-x',
                        '-U', '--mp', '--rdg', '--rfg', '--n-ceil', '-D', '-R', '-N', '-L', '-i', '--score-min', '--np'}
    found, positional, i = set(), [], 0
    while i < len(words):
        if words[i] in flags_with_value:
            found.add((words[i], words[i + 1]))
            i += 2
        elif words[i].startswith('-'):
            found.add((words[i], None))
            i += 1
        else:
            positional.append(words[i])
            i += 1
    return found, positional


@pytest.mark.parametrize('map_cmd', sorted(BEFORE))
def test_each_command_has_the_options_and_files_it_had_before(map_cmd):
    """Every command as before 0.4.0, apart from the changes of NOW (the three Bowtie 1 forms)."""
    row = {'map_cmd': map_cmd, 'aligner': mapSchedule.aligner(map_cmd),
           'options': mapSchedule.OPTIONS[map_cmd].format(multimap=100)}
    assert pairs(mapSchedule.command(row, 'IDX', 'IN', 'UNAL', 'SAM', 'LOG', 8)) == pairs(NOW[map_cmd])
    if map_cmd in ('bowtie1', 'mind', 'ribospike'):
        before, now = pairs(BEFORE[map_cmd]), pairs(NOW[map_cmd])
        assert before[0] - now[0] == {('--no-unal', None)}
        assert now[0] - before[0] == {('-f', None), ('-', None)}          # zcat -f IN | bowtie ... IDX - SAM
        assert now[1][:2] == ['zcat', 'IN'] and 'IN' not in now[1][2:]


def test_the_bowtie1_command_of_the_default_references_is_unchanged_character_for_character():
    row = {'map_cmd': 'bowtie', 'aligner': 'bowtie', 'options': U}
    assert mapSchedule.command(row, 'IDX', 'IN', 'UNAL', 'SAM', 'LOG', 8) == BEFORE['bowtie']
