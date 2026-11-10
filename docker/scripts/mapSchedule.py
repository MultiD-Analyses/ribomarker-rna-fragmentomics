"""The order file, the passes file and the resolved mapping schedule (release 0.4.0).

Order file (run --tiers FILE, config key tiers): tab-separated, a header line `tier`, then one reference name per
line: the index name without the species prefix (gtf/human-ncrna_snsno -> ncrna_snsno, dbPaths.tier_id).

Passes file (run --passes FILE, config key passes; experimental): tab-separated, one row per rule, columns PASS_COLUMNS.
Each pass aligns the reads still unaligned, through the references in the order in effect; a row says which settings
apply to which references.

The schedule is one row per cell (one reference in one pass) in execution order. It holds the aligner options exactly
as mapReads.py runs them, without threads, paths and input/output arguments, which mapReads adds (IO below). It is
written to <run>/tables/schedule.tsv by the Snakefile (rule run_tables) and read by mapReads.py, qcCollect.py,
qcReport.py and compare_reads.py.

Used by runCommand.py (checks and printout before the start), the Snakefile, mapReads.py, qcCollect.py, qcReport.py.
"""
import csv
import os
import shutil

import dbPaths

SPIKE_CMDS = ('ribospike', 'ribospikebt2', 'mind')    # the spike-in references, by aligner
FIRST = 'univec_bt1'     # vector contaminants and the spike-ins are removed first
LAST = 'genome'          # the genome catches what no annotated reference took

# Aligner options per map_cmd, as in mapReads.py before 0.4.0 (mapping_cmds), without threads and file arguments.
# {multimap} is filled in from the configuration. bowtie1, mind and ribospike no longer pass --no-unal: the pinned
# Bowtie 1 (1.0.0) refuses it ("unrecognized option"), so these three never ran with it; the unaligned records of a
# Bowtie 1 cell are dropped when its SAM becomes a BAM, as for `bowtie` (found 2026-11-05 by test_budgets.py).
OPTIONS = {
    'bowtie': '-S -k 1 -n 1 -l 18 -e 80 --chunkmbs 256',
    'bowtie1': '--fullref --best --strata -n 1 -l 18 -e 80 --chunkmbs 256 -m {multimap} -S',
    'mind': '--fullref -k 1 -v 0 --best --norc -S',
    'ribospike': '--fullref -k 1 -n 2 --best --norc -S',
    'ribospikebt2': '--sam-no-qname-trunc --no-unal --no-1mm-upfront --mp 6 --rdg 6,3 --rfg 6,3 --norc '
                    '--n-ceil L,5,0.15 -k 1 -D 20 -R 2 -N 0 -L 10 -i S,1,0.40 --score-min L,0,-0.6 --np 1 '
                    '--ignore-quals',
    'bowtie2': '--sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals -k {multimap}',
    'bowtie2best': '--sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals',
    'bowtie2sense': '--sam-no-qname-trunc --norc --no-unal --np 5 --very-sensitive --ignore-quals',
}

# The command of one cell: the aligner, the schedule's options, then threads and files. Every Bowtie 1 command reads
# its input from zcat -f (gzipped or not), as the `bowtie` command did before 0.4.0: a Bowtie 1 cell after a Bowtie 2
# cell gets a .fq.gz file (2026-11-05; before, bowtie1, mind and ribospike read the file themselves).
IO = {
    'bowtie': 'zcat -f {fastq} | bowtie {options} -p {threads} --un {unal} {index} - {sam} 2> {log}',
    'bowtie1': 'zcat -f {fastq} | bowtie {options} -p {threads} --un {unal} {index} - {sam} 2> {log}',
    'bowtie2': 'bowtie2 {options} --threads {threads} --un-gz {unal} -x {index} -U {fastq} -S {sam} 2> {log}',
}

COLUMNS = ['pass', 'step', 'reference', 'index_name', 'map_cmd', 'aligner', 'min_length', 'max_length',
           'max_mismatches', 'n_overhang', 'options']

PASS_COLUMNS = ['pass', 'applies_to', 'max_mismatches', 'min_length', 'max_length', 'n_overhang']
GROUPS = ('all', 'bowtie1', 'bowtie2', 'other')
MIN_LENGTH = 15                     # the trimming's shortest insert; a window may not start below it
BOWTIE1_SEED = ('-n', '-l', '-e', '-v')   # Bowtie 1 options replaced by -v m


class OrderError(Exception):
    """A reason to refuse an order file, written for the user."""


class PassesError(Exception):
    """A reason to refuse a passes file, written for the user."""


def aligner(map_cmd):
    """'bowtie' or 'bowtie2'."""
    if map_cmd in dbPaths.BOWTIE1_CMDS:
        return 'bowtie'
    if map_cmd in dbPaths.BOWTIE2_CMDS:
        return 'bowtie2'
    raise ValueError(f'unknown map_cmd {map_cmd!r}')


# ----------------------------------------------------------------------------- the order file
def read_order(path):
    """The reference names of an order file, in order. OrderError for the 0.3.0 format or a missing header."""
    try:
        with open(path, encoding='utf-8-sig') as f:
            lines = [line.strip() for line in f if line.strip()]
    except UnicodeDecodeError:
        raise OrderError(f'{path} is not a text file (UTF-8)')
    new_format = ('the order file has one column: a header line "tier", then one reference name per line '
                  '(copy tables/order.tsv of any 0.4.0 run and change the order)')
    if any(',' in line for line in lines):
        raise OrderError(f'{path} is a tier list of release 0.3.0 (map_cmd,index,retain); {new_format}')
    if not lines or lines[0] != 'tier':
        raise OrderError(f'{path}: the first line must be "tier"; {new_format}')
    return lines[1:]


def check_order(names, own, species, count_ribospike=True):
    """(problems, left_out) of an order file against the references folder's tier list.
    names: as read_order returns them; own: [(map_cmd, index, retain)] as dbPaths.read_dbpaths returns them.
    left_out: the folder's references not named (they are not used; a warning)."""
    by_name = {dbPaths.tier_id(index, species): map_cmd for map_cmd, index, _ in own}
    if not names:
        return ['the order file names no reference'], []
    problems = [f'{n} is not a reference of the references folder (its references: {", ".join(by_name)})'
                for n in names if n not in by_name]
    problems += [f'{n} is listed twice' for n in sorted({n for n in names if names.count(n) > 1})]
    first = [n for n in names if n == FIRST or by_name.get(n) in SPIKE_CMDS]
    if names[:len(first)] != first:
        problems.append(f'{" and ".join(first)} must come first (they remove vector and spike-in reads)')
    if LAST in names and names[-1] != LAST:
        problems.append(f'{LAST} must come last')
    left_out = [n for n in by_name if n not in names]
    spikes_left_out = [n for n in left_out if by_name[n] in SPIKE_CMDS]
    if spikes_left_out and count_ribospike:
        problems.append(f'the spike-in reference {", ".join(spikes_left_out)} is needed while count_ribospike is on')
    return problems, left_out


def check_dbpaths(own, path):
    """dbpaths.csv keeps three columns; the aligner (map_cmd) must be known, and the third column (retain) must be
    T since 0.4.0. Returns problems."""
    problems = [f'{path}: {index} has the unknown aligner {map_cmd} (known: {", ".join(OPTIONS)})'
                for map_cmd, index, _ in own if map_cmd not in OPTIONS]
    return problems + [f'{path}: {index} has retain {retain}; since release 0.4.0 every reference is kept (T)'
                       for _, index, retain in own if retain != 'T']


# ----------------------------------------------------------------------------- the passes file
def read_passes(path):
    """The rows of a passes file: [{column: value}] with 'row' (its line number); empty cells become '-'.
    PassesError if the file cannot be read or its header is not PASS_COLUMNS."""
    try:
        with open(path, encoding='utf-8-sig', newline='') as f:
            lines = [(n, line.rstrip('\r\n')) for n, line in enumerate(f, 1) if line.strip()]
    except UnicodeDecodeError:
        raise PassesError(f'{path} is not a text file (UTF-8)')
    if not lines or [c.strip() for c in lines[0][1].split('\t')] != PASS_COLUMNS:
        raise PassesError(f'{path}: the first line must be the columns {", ".join(PASS_COLUMNS)}, separated by tabs')
    rows = []
    for n, line in lines[1:]:
        cells = [c.strip() or '-' for c in line.split('\t')]
        if len(cells) != len(PASS_COLUMNS):
            raise PassesError(f'{path} line {n}: {len(cells)} columns instead of {len(PASS_COLUMNS)}')
        rows.append(dict(zip(PASS_COLUMNS, cells), row=n))
    if not rows:
        raise PassesError(f'{path} has no rows')
    return rows


def row_problems(row):
    """Problems of one row's values (not of the names it applies to)."""
    where = f"line {row['row']}"
    problems = []
    if not row['pass'].isdigit() or int(row['pass']) < 1:
        problems.append(f"{where}: pass {row['pass']!r} is not a whole number of 1 or more")
    if row['max_mismatches'] not in ('-', 'standard', '0', '1', '2'):
        problems.append(f"{where}: max_mismatches {row['max_mismatches']!r} is not standard, 0, 1 or 2")
    if row['n_overhang'] not in ('-', 'standard', '1', '2'):
        problems.append(f"{where}: n_overhang {row['n_overhang']!r} is not standard, 1 or 2")
    elif row['n_overhang'] in ('1', '2') and row['max_mismatches'] != '0':
        problems.append(f"{where}: n_overhang {row['n_overhang']} needs max_mismatches 0")
    lengths = {}
    for column in ('min_length', 'max_length'):
        value = row[column]
        if value != '-' and not (value.isdigit() and int(value) >= MIN_LENGTH):
            problems.append(f'{where}: {column} {value!r} is not a whole number of {MIN_LENGTH} or more')
        elif value != '-':
            lengths[column] = int(value)
    if len(lengths) == 2 and lengths['max_length'] < lengths['min_length']:
        problems.append(f'{where}: max_length is below min_length')
    return problems


def match_passes(rows, names, map_cmd_of):
    """(matched, problems): matched = {pass: {reference: row}} for the references of the order in effect (names,
    without the spike-ins), by the rows' applies_to; problems for names not in the order, spike-ins, and references
    matched twice in a pass."""
    spikes = {n for n in names if map_cmd_of[n] in SPIKE_CMDS}
    usable = [n for n in names if n not in spikes]
    groups = {'all': usable,
              'bowtie1': [n for n in usable if map_cmd_of[n] in dbPaths.BOWTIE1_CMDS],
              'bowtie2': [n for n in usable if map_cmd_of[n] in dbPaths.BOWTIE2_CMDS]}
    matched, problems = {}, []
    others = [r for r in rows if r['applies_to'] == 'other']
    for p in sorted({r['pass'] for r in others}):
        lines = [str(r['row']) for r in others if r['pass'] == p]
        if len(lines) > 1:
            problems.append(f"lines {' and '.join(lines)}: pass {p} has two rows with other")
    for row in sorted(rows, key=lambda r: r['applies_to'] == 'other'):      # 'other' after the rest of its pass
        where = f"line {row['row']}"
        given = [t.strip() for t in row['applies_to'].split(',') if t.strip()]
        in_pass = matched.setdefault(int(row['pass']), {})
        if given == ['other']:
            targets = [n for n in usable if n not in in_pass]
        elif len(given) == 1 and given[0] in groups:
            targets = groups[given[0]]
        else:
            if any(g in GROUPS for g in given):
                problems.append(f'{where}: a group ({", ".join(GROUPS)}) must stand alone in applies_to')
                continue
            targets = []
            for g in given:
                if g in spikes:
                    problems.append(f'{where}: {g} is a spike-in reference; it runs in pass 1 with its own settings')
                elif g not in names:
                    problems.append(f'{where}: {g} is not in the order in effect ({", ".join(names)})')
                else:
                    targets.append(g)
        for n in targets:
            if n in in_pass:
                problems.append(f"{where}: {n} is already matched in pass {row['pass']} by line {in_pass[n]['row']}")
            else:
                in_pass[n] = row
    return matched, problems


def check_passes(rows, names, map_cmd_of):
    """(matched, problems, warnings) of the rows of a passes file against the order in effect."""
    problems = [p for row in rows for p in row_problems(row)]
    if problems:
        return {}, problems, []
    passes = sorted({int(r['pass']) for r in rows})
    if passes != list(range(1, len(passes) + 1)):
        return {}, [f'pass numbers must be consecutive from 1 (found: {", ".join(map(str, passes))})'], []
    matched, problems = match_passes(rows, names, map_cmd_of)
    warnings = []
    if FIRST in names and FIRST not in matched.get(1, {}):
        warnings.append(f'{FIRST} has no row in pass 1: vector reads are not removed first')
    usable = [n for n in names if map_cmd_of[n] not in SPIKE_CMDS]
    unused = [n for n in usable if not any(n in m for m in matched.values())]
    if unused:
        warnings.append(f'{", ".join(unused)} run(s) in no pass')
    for row in rows:
        targets = [n for p, m in matched.items() for n, r in m.items() if r is row]
        wider = [n for n in targets if allows_more_than_standard(map_cmd_of[n], row)]
        if wider:
            warnings.append(f"line {row['row']}: max_mismatches {row['max_mismatches']} allows more than the standard "
                            f"settings for {', '.join(wider)}")
    return matched, problems, warnings


def allows_more_than_standard(map_cmd, row):
    """Bowtie 2: m mismatches cost 6m, the standard budget is 0.6 + 0.6 x length (at the row's min_length, 15 if not
    set); Bowtie 1: -v 2 allows two mismatches in the 18-nt seed, the standard -n 1 one."""
    if row['max_mismatches'] not in ('1', '2'):
        return False
    m = int(row['max_mismatches'])
    if aligner(map_cmd) == 'bowtie':
        return m == 2
    length = int(row['min_length']) if row['min_length'] != '-' else MIN_LENGTH
    return 6 * m > 0.6 + 0.6 * length


# ----------------------------------------------------------------------------- cells
def cell_options(map_cmd, max_mismatches, n_overhang, multimap):
    """The options of one cell. max_mismatches/n_overhang: 'standard' or a number as text."""
    words = OPTIONS[map_cmd].format(multimap=multimap).split()
    if max_mismatches == 'standard':
        return words
    m = int(max_mismatches)
    if aligner(map_cmd) == 'bowtie':                     # -v m in the place of the seed options
        place = next(i for i, w in enumerate(words) if w in BOWTIE1_SEED)
        kept = [w for i, w in enumerate(words)
                if not (w in BOWTIE1_SEED or (i > 0 and words[i - 1] in BOWTIE1_SEED))]
        return kept[:place] + ['-v', str(m)] + kept[place:]
    if n_overhang != 'standard':                         # only with max_mismatches 0 (row_problems)
        words[words.index('--np') + 1] = '1'
        return words + ['--score-min', f'C,-{n_overhang}']
    return words + ['--score-min', 'C,0' if m == 0 else f'C,-{6 * m}']


def standard(value):
    return 'standard' if value in ('-', 'standard') else value


def build(own, species, multimap, names=None, matched=None):
    """The schedule rows (dicts with COLUMNS), in execution order. names: the order in effect (default: the folder's
    order); matched: {pass: {reference: passes-file row}} (default: one pass, every reference, standard settings).
    The spike-ins run in pass 1 only, with their own settings. The last cell keeps the unaligned reads."""
    by_name = {dbPaths.tier_id(index, species): (map_cmd, index) for map_cmd, index, _ in own}
    names = names if names is not None else list(by_name)
    if matched is None:
        matched = {1: {n: {'max_mismatches': '-', 'n_overhang': '-', 'min_length': '-', 'max_length': '-'}
                       for n in names}}
    rows = []
    for p in sorted(matched):
        step = 0
        for name in names:
            map_cmd, index = by_name[name]
            spike = map_cmd in SPIKE_CMDS
            if spike and p != 1 or not spike and name not in matched[p]:
                continue
            row = matched[p].get(name, {'max_mismatches': '-', 'n_overhang': '-', 'min_length': '-', 'max_length': '-'})
            mm = 'spike-in' if spike else standard(row['max_mismatches'])
            n = 'spike-in' if spike else ('standard' if aligner(map_cmd) == 'bowtie' else standard(row['n_overhang']))
            options = OPTIONS[map_cmd].format(multimap=multimap).split() if spike else \
                cell_options(map_cmd, mm, n, multimap)
            step += 1
            rows.append({'pass': str(p), 'step': str(step), 'reference': name, 'index_name': os.path.basename(index),
                         'map_cmd': map_cmd, 'aligner': aligner(map_cmd),
                         'min_length': '-' if spike else row['min_length'],
                         'max_length': '-' if spike else row['max_length'],
                         'max_mismatches': mm, 'n_overhang': n, 'options': options})
    if rows:                                   # the last cell keeps the unaligned reads (as before 0.4.0)
        rows[-1]['options'] = [o for o in rows[-1]['options'] if o != '--no-unal']
    for row in rows:
        row['options'] = ' '.join(row['options'])
    return rows


def write_schedule(rows, path):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_schedule(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def order_in_effect(dbpaths, species, order_file=None):
    """The reference names in the order in effect: the order file's, else the references folder's (dbpaths.csv).
    The schedule alone does not give it when pass 1 does not cover every reference."""
    if order_file:
        return read_order(order_file)
    return [dbPaths.tier_id(index, species) for _, index, _ in dbPaths.read_dbpaths(dbpaths)]


def in_order(rows, names=None):
    """The references of the schedule's cells, each once: in the order in effect (names, e.g. of tables/order.tsv)
    if given, else in the order of their first cell (the same whenever pass 1 has every reference)."""
    present = dict.fromkeys(r['reference'] for r in rows)
    return [n for n in names if n in present] if names else list(present)


def write_order(rows, path, names=None):
    """tables/order.tsv: the order in effect, as an order file (copy and edit it for another order). names: the order
    in effect (order_in_effect); without it, the references of the cells (in_order)."""
    names = names if names is not None else in_order(rows)
    with open(path, 'w') as f:
        f.write('tier\n' + ''.join(f'{n}\n' for n in names))


def write_tiers(rows, dbpaths, species, path, names=None):
    """tables/tiers.csv: the dbpaths.csv lines of the references in effect, verbatim, in the order in effect (names,
    as in write_order). In the folder's own order (all references, same order) it is a byte copy of dbpaths.csv, as
    before 0.4.0."""
    with open(dbpaths, newline='') as f:
        lines = [line for line in f if line.strip()]
    by_name = {dbPaths.tier_id(line.split(',')[1].strip(), species): line for line in lines}
    order = names if names is not None else in_order(rows)
    if order == list(by_name):
        shutil.copyfile(dbpaths, path)
        return
    with open(path, 'w', newline='') as f:
        f.write(''.join(by_name[name] for name in order))


# ----------------------------------------------------------------------------- the printout
def window_words(row):
    low, high = row['min_length'], row['max_length']
    if low != '-' and high != '-':
        return f'{low}-{high} nt'
    if low != '-':
        return f'{low}+ nt'
    return f'up to {high} nt' if high != '-' else ''


def cell_words(row, last):
    """A schedule cell in words (as printed before a run)."""
    mm, n = row['max_mismatches'], row['n_overhang']
    if mm == 'spike-in':
        words = ['spike-in settings']
    elif mm == 'standard':
        words = ['standard']
    else:
        words = [f"{mm} mismatch{'' if mm == '1' else 'es'}"]
        if n not in ('standard', 'spike-in'):
            words.append(f'<={n} N')
    if window_words(row):
        words.append(window_words(row))
    if last:
        words.append('keeps unaligned')
    return ', '.join(words)


def render(rows, order_file=None, left_out=(), passes_file=None, names=None):
    """The printout: references down (those with a cell, in the order in effect: names, see in_order), passes
    across. order_file/passes_file: the files' names, or None."""
    order = f'order file {os.path.basename(order_file)}' if order_file else 'default (references folder)'
    settings = f'passes file {os.path.basename(passes_file)}' if passes_file else 'none'
    names = in_order(rows, names)
    passes = sorted({int(r['pass']) for r in rows})
    aligner_of = {r['reference']: r['aligner'] for r in rows}
    words = {(r['reference'], int(r['pass'])): cell_words(r, i == len(rows) - 1) for i, r in enumerate(rows)}
    columns = [[words.get((n, p), '-') for n in names] for p in passes]
    widths = [max(len(f'pass {p}'), *(len(w) for w in column)) for p, column in zip(passes, columns)]
    width = 1 + max(len('reference'), *(len(n) for n in names))

    def line(first, cells):
        padded = [c.ljust(w) for c, w in zip(cells[:-1], widths)] + cells[-1:]
        return (first + '  '.join(padded)).rstrip()
    plural = '' if len(passes) == 1 else 'es'
    lines = [f'Mapping schedule: {len(names)} references, {len(passes)} pass{plural}; order: {order}; '
             f'pass settings: {settings}',
             line(f"  step  {'reference':<{width}}  aligner  ", [f'pass {p}' for p in passes])]
    for i, n in enumerate(names, 1):
        lines.append(line(f"  {i:>4}  {n:<{width}}  {aligner_of[n]:<7}  ", [column[i - 1] for column in columns]))
    if left_out:
        lines.append(f"  not used: {', '.join(left_out)}")
    return '\n'.join(lines)


# ----------------------------------------------------------------------------- the run's schedule
def resolve(dbpaths, species, multimap, order_file=None, count_ribospike=True, passes_file=None):
    """(rows, left_out, problems, warnings) for a run: the references folder's dbpaths.csv, the optional order file
    and the optional passes file. rows is empty when there are problems."""
    own = dbPaths.read_dbpaths(dbpaths)
    problems = check_dbpaths(own, dbpaths)
    names, left_out, warnings = None, [], []
    if order_file:
        try:
            names = read_order(order_file)
        except OrderError as error:
            return [], [], problems + [str(error)], []
        order_problems, left_out = check_order(names, own, species, count_ribospike)
        problems += order_problems
    if problems:
        return [], left_out, problems, []
    matched = None
    if passes_file:
        map_cmd_of = {dbPaths.tier_id(index, species): map_cmd for map_cmd, index, _ in own}
        in_effect = names if names is not None else list(map_cmd_of)
        try:
            matched, passes_problems, warnings = check_passes(read_passes(passes_file), in_effect, map_cmd_of)
        except PassesError as error:
            return [], left_out, [str(error)], []
        if passes_problems:
            return [], left_out, [f'{passes_file}: {p}' for p in passes_problems], []
        warnings = [f'{os.path.basename(passes_file)}: {w}' for w in warnings]
    return build(own, species, multimap, names, matched), left_out, [], warnings


def command(row, index, fastq, unal, sam, log, threads):
    """The shell command of one cell (mapReads.py). unal: the unaligned reads' file (.fq for Bowtie 1, .fq.gz for
    Bowtie 2)."""
    if row['map_cmd'] == 'bowtie':
        template = IO['bowtie']
    elif row['aligner'] == 'bowtie':
        template = IO['bowtie1']
    else:
        template = IO['bowtie2']
    return template.format(options=row['options'], threads=threads, index=index, fastq=fastq, unal=unal, sam=sam,
                           log=log)
