"""run: the whole pipeline on the samples of a sample sheet, with one command.

  python docker/scripts/main.py run SHEET OUTPUT [--species NAME] [--threads N]
                                    [--config FILE] [--run-id NAME] [--references FOLDER] [--unlock]

run reads the samples from SHEET (a CSV file in the folder of the FASTQ files), checks everything it can
before starting, writes the pipeline's sample sheet and the configuration into OUTPUT/<run id>/, runs the pipeline with OUTPUT as the working folder and
ends with one line: RESULTS VALID or RESULTS NOT VALID, and the path of the QC report.

Exit status: 0 = results valid; 1 = the pipeline failed or the results are not valid; 2 = refused
before starting (the message says why and what to do); 130 = interrupted (Ctrl-C).
The user-facing help text is HELP (run --help).

"""
import argparse
import csv
import datetime as dt
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

import checkTruth
import codeFingerprint
import dbPaths

PIPELINE_DIR = Path(os.environ.get('RIBOMARKER_PIPELINE_DIR', Path(__file__).resolve().parents[1]))
REFERENCES = Path(os.environ.get('RIBOMARKER_REFERENCES', PIPELINE_DIR / 'assets'))
DEMO = PIPELINE_DIR / 'demo'      # run --demo: fastq/, references/, truth.tsv, made by docker/make_demo.sh

NAME = re.compile(r'^[A-Za-z][A-Za-z0-9_-]*$')        # run ids and species names
BUILTIN_SPECIES = ('human', 'mouse')                  # the references built into the image
SAMPLE = re.compile(r'^[A-Za-z0-9_][A-Za-z0-9_-]*$')  # sample names: letters, digits, "_", "-"; no leading "-"
FILE = re.compile(r'^[A-Za-z0-9_.][A-Za-z0-9_.-]*$')  # FASTQ file names: also "."; no spaces, no paths
# sample names that the pipeline's pandas would read as a number, a missing value or True/False, not as text
NOT_TEXT = re.compile(r'^(\d+([eE]-?\d+)?|na|nan|null|none|inf|infinity|true|false)$', re.IGNORECASE)
FASTQ_ENDINGS = ('.fastq.gz', '.fq.gz')
COLUMNS = ('sample', 'fastq_file', 'fastq_file2')
GROUP = 'RiboMarker'                                    # the sample sheet's group column (not used)

DEFAULTS = {                                            # today's validated settings
    'remove_multimap': True, 'remove_antisense': True, 'collapse': False, 'keep_types': 'all',
    'qc_report': True, 'trim_fastq': True, 'count_ribospike': True,
    'multimap': 100, 'trim_only': False, 'subsample_reads': 0,
}
ALLOWED_KEYS = {'multimap', 'keep_types', 'remove_multimap', 'remove_antisense',
                'count_ribospike', 'qc_report'}         # --config under run
BOOLEAN_KEYS = {'remove_multimap', 'remove_antisense', 'count_ribospike', 'qc_report'}
TYPE_LIST = re.compile(r'^[A-Za-z0-9_.+-]+(,[A-Za-z0-9_.+-]+)*$')   # keep_types: all, or miRNA,tRNA,...
SET_BY_RUN = {'dbdir', 'input_dir', 'output_prefix', 'samples', 'host_input_dir', 'host_assets_dir',
              'species', 'threads', 'paired_end', 'report_title'}

DEFAULT_THREADS = 8          # the starters' default too; MINUTES_PER_R1_GB was measured with 8
MIN_MEMORY_GB = 12           # warning below this (8 threads: peak 4.1 GB with 2 samples, 8.1 GB with 6
                             # samples whose largest had 0.67 GB of R1; the peak grows with sample size)
MINUTES_PER_R1_GB = 165      # measured on a laptop, 8 threads: 37-41 min for 0.30 GB of R1 files (2 samples),
                             # 281 min for 1.72 GB (6 samples, Windows starter); the slower one is used
DISK_FACTOR = 3              # free space needed: 3 x input + 1 GB (level C: outputs 1.5 x input; peak not measured)


class RunError(Exception):
    """A reason to refuse the run, written for the user."""


def say(message=''):
    print(message, flush=True)


# ----------------------------------------------------------------------------- the sample sheet
def read_sheet(path):
    """The user's sample sheet: ([(sample, fastq_file, fastq_file2 or '')], 'pe' or 'se').
    One row per sample; the columns sample, fastq_file and (paired-end) fastq_file2 by name, in any order;
    other columns are left alone (the sheet is copied to the results). Comma or semicolon separated, UTF-8
    with or without a byte-order mark. The FASTQ files are in the sheet's folder, named without a path.
    All problems are reported together."""
    path = Path(path)
    if not path.is_file():
        raise RunError(f'sample sheet not found: {path}')
    try:
        text = path.read_text(encoding='utf-8-sig')
    except UnicodeDecodeError:
        raise RunError(f'{path} is not saved as UTF-8 text; save it again as "CSV UTF-8" (Excel: File > Save As)')
    lines = text.splitlines()
    first = lines[0] if lines else ''
    delimiter = ';' if first.count(';') > first.count(',') else ','
    rows = [(n, [c.strip() for c in row]) for n, row in enumerate(csv.reader(lines, delimiter=delimiter), 1)
            if any(c.strip() for c in row)]
    if not rows:
        raise RunError(f'{path} is empty')
    header = rows[0][1]
    missing = [c for c in COLUMNS[:2] if c not in header]
    twice = [c for c in COLUMNS if header.count(c) > 1]
    if missing or twice:
        raise RunError(f'{path}: the first row must name the columns sample, fastq_file and, for paired-end '
                       'reads, fastq_file2 (written exactly so, in any order); '
                       + (f'missing: {", ".join(missing)}' if missing else f'given twice: {", ".join(twice)}')
                       + f'\n  first row found: {", ".join(header)}')
    if len(rows) < 2:
        raise RunError(f'{path} has no rows with samples')
    index = {c: header.index(c) for c in COLUMNS if c in header}

    def cell(row, column):
        i = index.get(column)
        return row[i] if i is not None and i < len(row) else ''
    rows = [(n, cell(r, 'sample'), cell(r, 'fastq_file'), cell(r, 'fastq_file2')) for n, r in rows[1:]]

    problems = []
    for n, name, r1, _ in rows:
        if not name:
            problems.append(f'row {n}: no sample name')
        if not r1:
            problems.append(f'row {n}: no fastq_file' + (f' for {name}' if name else ''))
    names = [name for _, name, _, _ in rows if name]
    bad = [x for x in names if not SAMPLE.match(x)]
    if bad:
        problems.append('sample names may contain only letters, digits, "_" and "-", and must not start with "-": '
                        + ', '.join(bad))
    not_text = [x for x in names if SAMPLE.match(x) and NOT_TEXT.match(x)]
    if not_text:
        problems.append('these sample names would be read as numbers or missing values by the pipeline; '
                        'add a letter, e.g. "S001" instead of "001": ' + ', '.join(not_text))
    kits = [x for x in names if 'NEB' in x or 'SMART' in x]
    if kits:
        problems.append('sample names containing "NEB" or "SMART" make the pipeline treat the sample as another '
                        'library kit, so they cannot be used: ' + ', '.join(kits))
    same = {}
    for x in names:
        same.setdefault(x.lower(), []).append(x)
    twice = [' / '.join(v) for v in same.values() if len(v) > 1]
    if twice:
        problems.append('sample names given more than once (names that differ only in upper and lower case count '
                        'as the same, because they become file names): ' + '; '.join(twice))

    with_r2 = [n for n, _, _, r2 in rows if r2]
    if with_r2 and len(with_r2) != len(rows):
        without = [name or f'row {n}' for n, name, _, r2 in rows if not r2]
        problems.append('fastq_file2 is filled in for some samples but not for others; a run is either paired-end '
                        '(fastq_file2 for every sample) or single-end (for none). Without fastq_file2: '
                        + ', '.join(without))
    files = [f for _, _, r1, r2 in rows for f in (r1, r2) if f]
    with_path = [f for f in files if '/' in f or '\\' in f]
    if with_path:
        problems.append(f'fastq_file and fastq_file2 must be names of files in the folder of the sample sheet '
                        f'({path.parent}), without a path: ' + ', '.join(with_path))
    bad_name = [f for f in files if f not in with_path and not FILE.match(f)]
    if bad_name:
        problems.append('FASTQ file names may contain only letters, digits, ".", "_" and "-" (no spaces or '
                        'brackets) and must not start with "-": ' + ', '.join(bad_name))
    not_gz = [f for f in files if not f.lower().endswith(FASTQ_ENDINGS)]
    if not_gz:
        problems.append('FASTQ files must be gzip-compressed and end in .fastq.gz or .fq.gz: ' + ', '.join(not_gz))
    not_found = [f for f in files if f not in with_path and not (path.parent / f).is_file()]
    if not_found:
        problems.append(f'not found in {path.parent}: ' + ', '.join(not_found))
    repeated = sorted({f for f in files if files.count(f) > 1})
    if repeated:
        problems.append('files listed more than once (each file belongs to one sample and one read): '
                        + ', '.join(repeated))
    if problems:
        raise RunError(f'the sample sheet {path} cannot be used:\n  ' + '\n  '.join(problems))
    return [(name, r1, r2) for _, name, r1, r2 in rows], ('pe' if with_r2 else 'se')


def sheet_rows(samples):
    """Rows of the pipeline's own sample sheet, in the order of the user's sheet: sample, group, R1[, R2]."""
    return [[name, GROUP, r1] + ([r2] if r2 else []) for name, r1, r2 in samples]


# ----------------------------------------------------------------------------- configuration
def rejected_reason(key, value):
    """Why a --config key/value cannot be used under run, or None."""
    if key in SET_BY_RUN:
        return 'set by run itself (from the command line or the files)'
    if key == 'trim_fastq' and value is False:
        return 'run needs trimming'
    if key == 'trim_only' and value is True:
        return 'a trimming-only run has no QC report'
    if key == 'collapse' and value is True:
        return ('collapse: True maps each distinct sequence once and splits multimapper counts over the tied '
                'loci, so the QC read counts (checks C2-C4) are not read counts')
    if key in DEFAULTS and value == DEFAULTS[key] and type(value) is type(DEFAULTS[key]):
        return None                                   # the default value itself is harmless
    if key not in ALLOWED_KEYS:
        return 'not available under run (not tested); use the config-file route for it'
    return value_problem(key, value)


def value_problem(key, value):
    """Why the value of an allowed key cannot be used, or None. A quoted "False" would count as true."""
    is_whole_number = isinstance(value, int) and not isinstance(value, bool)
    if key in BOOLEAN_KEYS and not isinstance(value, bool):
        return f'{value!r} is not True or False (write it without quotes)'
    if key == 'multimap' and not (is_whole_number and value >= 1):
        return f'{value!r} is not a whole number of 1 or more'
    if key == 'keep_types' and not (isinstance(value, str) and TYPE_LIST.match(value)):
        return f'{value!r} is not "all" or RNA types separated by commas without spaces (e.g. miRNA,tRNA)'
    return None


def user_config(path):
    """The allowed settings of a --config file; RunError listing every key that cannot be used."""
    if path is None:
        return {}
    try:
        with open(path, encoding='utf-8') as f:
            values = yaml.safe_load(f) or {}
    except (UnicodeDecodeError, yaml.YAMLError) as error:
        raise RunError(f'{path} cannot be read as a YAML file of "key: value" lines ({type(error).__name__})')
    if not isinstance(values, dict):
        raise RunError(f'{path} is not a YAML file of "key: value" lines')
    problems = [f'  {key}: {reason}' for key, value in values.items()
                if (reason := rejected_reason(key, value))]
    if problems:
        raise RunError(f'{path} contains settings run cannot use:\n' + '\n'.join(problems))
    return values


def pipeline_config(settings, *, input_dir, references, run_id, sheet, species, threads, mode, title):
    """The configuration run gives to Snakemake. Folders absolute; output_prefix the bare run id."""
    config = dict(DEFAULTS)
    config.update(settings)
    config.update({'input_dir': f'{input_dir}/', 'dbdir': f'{references}/', 'output_prefix': run_id,
                   'samples': str(sheet), 'species': species, 'threads': threads,
                   'paired_end': mode == 'pe', 'report_title': title})
    return config


# ----------------------------------------------------------------------------- checks before starting
def check_references(references, species):
    """The tiers of the species exist; returns warnings (e.g. no MANIFEST)."""
    try:
        dbpaths = dbPaths.find_dbpaths(str(references), species)
    except FileNotFoundError as error:
        raise RunError(str(error))
    missing = [index for _, index, _ in dbPaths.read_dbpaths(dbpaths)
               if not any(os.path.exists(index + ext) for ext in ('.1.bt2', '.1.ebwt', '.1.bt2l'))]
    bed = Path(references) / f'{species}.bed'
    if missing or not bed.exists():
        raise RunError(f'the {species} references in {references} are incomplete; missing: '
                       + ', '.join(missing + ([] if bed.exists() else [str(bed)])))
    if not (Path(references) / f'MANIFEST_{species}.json').exists():
        return [f'no MANIFEST_{species}.json in {references}: the references cannot be identified in the report']
    return []


def broken_gzip(files, threads):
    """Files that fail `gzip -t` (truncated or damaged), checked in parallel."""
    def bad(path):
        return subprocess.run(['gzip', '-t', str(path)], capture_output=True).returncode != 0
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        return [f for f, is_bad in zip(files, pool.map(bad, files)) if is_bad]


def memory_gb():
    """Memory available to this process: the container's limit if there is one, else the machine's."""
    for path in ('/sys/fs/cgroup/memory.max', '/sys/fs/cgroup/memory/memory.limit_in_bytes'):
        try:
            value = Path(path).read_text().strip()
            if value.isdigit() and int(value) < 1 << 60:
                return int(value) / 1e9
        except OSError:
            pass
    try:
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemTotal:'):
                return int(line.split()[1]) * 1024 / 1e9
    except OSError:
        pass
    return None


def run_state(run_dir, sheet_text, config_text):
    """'new' or 'resume'; RunError if the folder holds a finished run or a different one."""
    status_file = run_dir / 'run_status.json'
    if not run_dir.exists():
        return 'new'
    if not status_file.exists():
        raise RunError(f'{run_dir} exists but was not made by run; choose another --run-id or output folder')
    status = json.loads(status_file.read_text())
    if status.get('status') == 'finished':
        raise RunError(f'{run_dir} holds a finished run ({status.get("verdict", "")}). '
                       'Choose another --run-id or output folder for a new run.')
    def same(name, text):
        path = run_dir / name
        return path.exists() and without_threads(path.read_text()) == without_threads(text)
    if not (same('samples.csv', sheet_text) and same('run_config.yaml', config_text)):
        raise RunError(f'{run_dir} holds an unfinished run of other samples or settings. Choose another '
                       '--run-id or output folder, or repeat the earlier command to resume that run.')
    return 'resume'


def without_threads(config_text):
    """Threads do not change the results, so a run may be resumed with another --threads."""
    return '\n'.join(line for line in config_text.splitlines() if not line.startswith('threads:'))


# ----------------------------------------------------------------------------- verdict
def verdict(run_dir):
    """(valid, one-line verdict) from qc/consistency_checks.tsv and the warnings of qc/qc_metrics.json."""
    checks_file = run_dir / 'qc' / 'consistency_checks.tsv'
    if not checks_file.exists():
        return True, 'RESULTS NOT CHECKED: the QC report is switched off (qc_report: False)'
    with open(checks_file, newline='') as f:
        statuses = [row['status'] for row in csv.DictReader(f, delimiter='\t')]
    counts = {s: statuses.count(s) for s in sorted(set(statuses))}
    warnings = json.loads((run_dir / 'qc' / 'qc_metrics.json').read_text()).get('warnings', [])
    bad = counts.get('differs', 0) + counts.get('not evaluable', 0)
    detail = ', '.join(f'{n} {s}' for s, n in counts.items()) + (f'; {len(warnings)} warning(s)' if warnings else '')
    report = run_dir / 'qc' / 'report_technical.html'
    if bad:
        return False, f'RESULTS NOT VALID: {bad} of {len(statuses)} consistency checks failed ({detail}) - see {report}'
    return True, f'RESULTS VALID: all {counts.get("equal", 0)} consistency checks equal ({detail})'


# ----------------------------------------------------------------------------- the command
HELP = """Runs the whole pipeline on the samples of a sample sheet and ends with one line: RESULTS VALID or
RESULTS NOT VALID, and where the QC report is.

In Docker:
  docker run --rm -it --init --network none --user "$(id -u):$(id -g)" \\
    -v /path/to/fastq:/data/input:ro -v /path/to/output:/data/output \\
    IMAGE run /data/input/SHEET /data/output [options]
  Windows PowerShell: the same on one line, without --user, with Windows paths
    (-v C:\\data\\fastq:/data/input:ro -v C:\\data\\results:/data/output).
Local (conda environment ribomarker active):
  python docker/scripts/main.py run SHEET OUTPUT [options]

SHEET is a CSV file in the folder of the FASTQ files, with one row per sample and the columns
sample (the name used in all results), fastq_file (the reads, R1) and, for paired-end reads,
fastq_file2 (R2); file names without a path. Further columns (e.g. clinical data) are copied to the
results unchanged. Comma or semicolon separated. The folder is only read; the results go into
OUTPUT/<run id>/.

Repeat the same command to resume an unfinished run; a finished run is refused. Ctrl-C stops
a run cleanly.

Test the installation in a few minutes (synthetic reads whose results are known):
  run --demo OUTPUT          (Docker: ... IMAGE run --demo /data/output)

IMAGE is the pipeline's container image; its name and version are in the README.

Exit status: 0 = results valid; 1 = the pipeline failed or the results are not valid;
2 = refused before starting (the message says why and what to do); 130 = interrupted.
"""


def positive_int(text):
    """--threads: a whole number of 1 or more (otherwise argparse refuses with exit status 2)."""
    if not text.isdigit() or int(text) < 1:
        raise argparse.ArgumentTypeError(f'{text!r} is not a whole number of 1 or more')
    return int(text)


def parse_args(argv):
    ap = argparse.ArgumentParser(prog='run', description=HELP,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('sheet', help='sample sheet (CSV) in the folder of the FASTQ files')
    ap.add_argument('output', help='folder for the results (created if missing)')
    ap.add_argument('--species', default='human',
                    help='reference species: human (default) or mouse with the built-in references; with '
                         '--references FOLDER, the species name of the files in that folder')
    ap.add_argument('--threads', type=positive_int, default=DEFAULT_THREADS,
                    help='cores to use (default: %(default)s; at most the cores available)')
    ap.add_argument('--config', help='advanced settings: a YAML file with only these keys: '
                    + ', '.join(sorted(ALLOWED_KEYS)))
    ap.add_argument('--run-id', help='name of the run folder and the report title (default: ribomarker)')
    ap.add_argument('--references', default=str(REFERENCES), help='references folder (default: %(default)s)')
    ap.add_argument('--unlock', action='store_true',
                    help='remove a lock left by an interrupted run (only if no other run uses OUTPUT)')
    return ap.parse_args(argv)


def utc_now():
    return dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def prepare(args):
    """All checks and files before Snakemake starts. Returns what the run needs, or raises RunError."""
    sheet_path, output = Path(args.sheet).resolve(), Path(args.output).resolve()
    references = Path(args.references).resolve()
    if args.config is not None and not Path(args.config).is_file():
        raise RunError(f'--config: file not found: {args.config}')
    samples, mode = read_sheet(sheet_path)
    input_dir = sheet_path.parent
    if output == input_dir or input_dir in output.parents:
        raise RunError('the output folder must not be the folder of the sample sheet or lie inside it')
    for label, path in (('sample sheet', input_dir), ('output', output), ('references', references)):
        if ' ' in str(path):
            raise RunError(f'the {label} folder path contains a space ({path}); the pipeline cannot handle '
                           'spaces in paths. Use a folder without spaces (in Docker this does not apply).')
    run_id = args.run_id or 'ribomarker'
    if not NAME.match(run_id):
        raise RunError(f'--run-id {run_id!r}: use letters, digits, "_" and "-" only, starting with a letter')
    if not NAME.match(args.species):
        raise RunError(f'--species {args.species!r}: use letters, digits, "_" and "-" only, starting with a letter')
    if args.species not in BUILTIN_SPECIES and references == REFERENCES.resolve():
        raise RunError(f'--species {args.species}: the built-in references are human and mouse; for another '
                       'species give your own references with --references FOLDER (user guide, "Own references")')

    settings = user_config(args.config)
    warnings = check_references(references, args.species)
    if args.species == 'mouse':
        warnings.append('mouse references not validated on mouse data: the numbers of this run '
                        'have not been compared with a known mouse result')
    elif args.species not in BUILTIN_SPECIES:
        warnings.append(f'species {args.species}: own references; the pipeline has been tested with human and '
                        'mouse only, and the species-specific QC items are skipped')

    threads = args.threads
    cores = os.cpu_count() or 1
    if threads > cores:
        warnings.append(f'{threads} threads asked for, but only {cores} cores are available; using {cores}')
        threads = cores
    memory = memory_gb()
    if memory is not None and memory < MIN_MEMORY_GB:
        warnings.append(f'only {memory:.1f} GB of memory available; large samples need about 8 GB at the peak '
                        f'and {MIN_MEMORY_GB} GB or more is recommended (Docker Desktop: Settings > Resources)')

    files = [input_dir / f for _, r1, r2 in samples for f in (r1, r2) if f]
    input_bytes = sum(f.stat().st_size for f in files)
    output.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(output).free
    needed = DISK_FACTOR * input_bytes + 1e9
    if free < needed:
        raise RunError(f'not enough free disk space in {output}: {free / 1e9:.1f} GB free, about '
                       f'{needed / 1e9:.1f} GB needed for {input_bytes / 1e9:.1f} GB of FASTQ files')
    run_dir = output / run_id
    status_file = run_dir / 'run_status.json'
    first_start = json.loads(status_file.read_text()).get('first_started_utc') if status_file.exists() else None
    first_start = first_start or utc_now()
    title = run_id if args.run_id else f'RiboMarker run {first_start[:10]}'
    sheet = run_dir / 'samples.csv'
    sheet_text = ''.join(','.join(row) + '\n' for row in sheet_rows(samples))
    config = pipeline_config(settings, input_dir=input_dir, references=references, run_id=run_id,
                             sheet=sheet, species=args.species, threads=threads, mode=mode, title=title)
    config_text = yaml.safe_dump(config, sort_keys=False)
    state = run_state(run_dir, sheet_text, config_text)   # before gzip -t: a finished run is refused at once
    say(f'checking {len(files)} FASTQ files (gzip -t) ...')
    broken = broken_gzip(files, threads)
    if broken:
        raise RunError('these FASTQ files are damaged or incomplete (gzip -t failed); copy them again:\n  '
                       + '\n  '.join(str(f) for f in broken))

    r1_bytes = sum((input_dir / r1).stat().st_size for _, r1, _ in samples)
    minutes = MINUTES_PER_R1_GB * r1_bytes / 1e9 * 8 / threads
    return dict(input_dir=input_dir, sheet_path=sheet_path, output=output, run_dir=run_dir, run_id=run_id,
                title=title, mode=mode, samples=samples, warnings=warnings, threads=threads, state=state,
                sheet_text=sheet_text, config_text=config_text, first_start=first_start,
                input_bytes=input_bytes, minutes=minutes)


def pipeline_version():
    """Version, code fingerprint and git commit of the pipeline, the same values as in the QC report:
    in the Docker image from the files written at build time, in a checkout computed from the code."""
    def text(name):
        path = PIPELINE_DIR / name
        return path.read_text().strip() if path.exists() else None
    try:
        fingerprint = text('CODE_FINGERPRINT') or codeFingerprint.fingerprint(PIPELINE_DIR)
    except OSError:
        fingerprint = 'unknown'
    commit = text('GIT_COMMIT')
    if commit is None:
        try:
            result = subprocess.run(['git', '-C', str(PIPELINE_DIR), 'rev-parse', 'HEAD'],
                                    capture_output=True, text=True)
            commit = result.stdout.strip() if result.returncode == 0 else 'unknown'
        except OSError:
            commit = 'unknown'
    return {'version': text('VERSION') or 'unknown', 'code_fingerprint': fingerprint, 'git_commit': commit}


def write_status(run_dir, **fields):
    path = run_dir / 'run_status.json'
    status = json.loads(path.read_text()) if path.exists() else {}
    status.update(fields)
    path.write_text(json.dumps(status, indent=2) + '\n')


def copy_output(process, log):
    for line in process.stdout:
        sys.stdout.write(line)
        log.write(line)


# The memory the container's processes hold, without the file cache (which the kernel can free):
# 'anon' in cgroup v2, 'total_rss' in cgroup v1 (Docker Desktop on WSL2, kernel 5.15: v1).
MEMORY_STAT = [('/sys/fs/cgroup/memory.stat', 'anon'), ('/sys/fs/cgroup/memory/memory.stat', 'total_rss')]
MEMORY_SAMPLE_SECONDS = 5


def container_memory_gb():
    """Memory the container's processes hold now, or None outside a container."""
    if not os.environ.get(dbPaths.CONTAINER_ENV):
        return None
    for path, key in MEMORY_STAT:
        try:
            for line in Path(path).read_text().splitlines():
                name, value = line.split()
                if name == key:
                    return int(value) / 1e9
        except (OSError, ValueError):
            continue
    return None


def watch_memory(process, peak):
    """Until the process ends: the highest container_memory_gb() in peak['gb'] (sampled every few seconds,
    enough for the plateaus of 1-2 min measured on level C)."""
    while process.poll() is None:
        now = container_memory_gb()
        if now is not None:
            peak['gb'] = max(peak.get('gb', 0.0), now)
        time.sleep(MEMORY_SAMPLE_SECONDS)


def size_text(n_bytes):
    return f'{n_bytes / 1e9:.1f} GB' if n_bytes >= 1e8 else f'{n_bytes / 1e6:.0f} MB'


def duration_text(minutes):
    return f'{minutes / 60:.1f} h' if minutes >= 60 else f'{max(1, round(minutes))} min'


def run_snakemake(p, unlock):
    """Snakemake in OUTPUT, output to the terminal and to <run>/run.log. Returns the exit status."""
    base = ['snakemake', '--snakefile', str(PIPELINE_DIR / 'Snakefile'),
            '--configfile', str(p['run_dir'] / 'run_config.yaml')]
    if shutil.which('snakemake') is None:
        raise RunError('snakemake not found: activate the pipeline environment (conda activate ribomarker)')
    if unlock:
        subprocess.run(base + ['--unlock'], cwd=p['output'], capture_output=True)
    log = p['run_dir'] / 'run.log'
    with open(log, 'a') as out:
        # Snakemake runs in its own process group, so a Ctrl-C reaches it exactly once, from here: from a
        # terminal, and also in Docker without one, where only the first process gets the signal.
        # Snakemake then stops its jobs and removes its lock; its last messages still reach run.log, and
        # main() marks the run as interrupted. Closing the terminal (SIGHUP) or `docker stop` (SIGTERM)
        # stop it the same way.
        process = subprocess.Popen(base + ['--cores', str(p['threads']), '--rerun-incomplete'],   # no --keep-incomplete
                                   cwd=p['output'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, errors='replace', start_new_session=True,
                                   env=dict(os.environ, PYTHONUNBUFFERED='1'))   # the scripts' lines reach run.log in order
        interrupted = []

        def stop(signum, frame):
            interrupted.append(signum)
            try:
                os.killpg(process.pid, signal.SIGINT)
            except ProcessLookupError:
                pass                                   # Snakemake has already ended
        previous = {s: signal.signal(s, stop) for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}
        peak = p.setdefault('memory_peak', {})
        threading.Thread(target=watch_memory, args=(process, peak), daemon=True).start()
        try:
            copy_output(process, out)
            status = process.wait()
        finally:
            for s, handler in previous.items():
                signal.signal(s, handler)
    if interrupted:
        raise KeyboardInterrupt
    if status != 0 and 'LockException' in log.read_text(errors='replace'):
        say(f'\nThe output folder is locked by an earlier run that was interrupted. If no other run is using '
            f'{p["output"]}, repeat the command with --unlock.')
    return status


def run_demo(argv):
    """run --demo OUTPUT [--threads N]: the synthetic level A reads (paired-end) and their small references
    through the whole pipeline, then every number with a known answer compared with it."""
    ap = argparse.ArgumentParser(prog='run --demo', description='Tests the installation in a few minutes: '
                                 'synthetic reads whose results are known go through the whole pipeline.')
    ap.add_argument('output', help='folder for the demo results (created if missing)')
    ap.add_argument('--threads', type=positive_int, default=DEFAULT_THREADS)
    args = ap.parse_args(argv)
    if not (DEMO / 'truth.tsv').exists():
        say(f'ERROR: no demo data in {DEMO} (docker/make_demo.sh makes it, before the image is built)')
        return 2
    status = main([str(DEMO / 'fastq' / 'samples.csv'), args.output, '--references', str(DEMO / 'references'),
                   '--threads', str(args.threads), '--run-id', 'demo'])
    if status == 2:
        say(f'DEMO NOT STARTED: see above (for another demo, use an empty output folder)')
        return 2
    if status != 0:
        say('DEMO FAILED: the pipeline did not end with valid results (see above)')
        return status
    results, table = checkTruth.check(DEMO / 'truth.tsv', 'pe', Path(args.output) / 'demo')
    passed = [r for r in results if r[0] == 'PASS']
    failed = [r for r in results if r[0] == 'FAIL']
    if failed:
        for r in failed[:10]:
            say(f'  {r[1]} {r[2]} {r[3]}: expected {r[4]}, observed {r[5]}')
        say(f'DEMO FAILED: {len(failed)} of {len(passed) + len(failed)} known results differ - see {table}')
        return 1
    say(f'DEMO PASSED: all {len(passed)} known results reproduced (details: {table})')
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ['--demo']:
        return run_demo(argv[1:])
    args = parse_args(argv)
    try:
        p = prepare(args)
    except RunError as error:
        say(f'ERROR: {error}')
        return 2
    for warning in p['warnings']:
        say(f'WARNING: {warning}')
    p['run_dir'].mkdir(parents=True, exist_ok=True)
    (p['run_dir'] / 'samples.csv').write_text(p['sheet_text'])
    (p['run_dir'] / 'run_config.yaml').write_text(p['config_text'])
    shutil.copyfile(p['sheet_path'], p['run_dir'] / 'sample_sheet.csv')        # the user's sheet, unchanged
    write_status(p['run_dir'], run_id=p['run_id'], title=p['title'], first_started_utc=p['first_start'],
                 started_utc=utc_now(), status='running', samples=len(p['samples']), mode=p['mode'],
                 pipeline_version=pipeline_version())
    say(f"{'resuming' if p['state'] == 'resume' else 'starting'} {p['title']}: {len(p['samples'])} samples "
        f"({'paired-end' if p['mode'] == 'pe' else 'single-end'}), {size_text(p['input_bytes'])} of FASTQ files, "
        f"{p['threads']} threads; roughly {duration_text(p['minutes'])} (keep the computer awake)")
    say(f"results: {p['run_dir']}")

    try:
        status = run_snakemake(p, args.unlock)
    except RunError as error:
        say(f'ERROR: {error}')
        write_status(p['run_dir'], status='failed', finished_utc=utc_now(), verdict=str(error))
        return 2
    except KeyboardInterrupt:
        write_status(p['run_dir'], status='interrupted', finished_utc=utc_now())
        say('\nINTERRUPTED: repeat the same command to resume '
            '(if it then says the output folder is locked, add --unlock)')
        return 130
    peak = p.get('memory_peak', {}).get('gb')
    if peak is not None:
        write_status(p['run_dir'], memory_peak_gb=round(peak, 2))
        say(f'memory peak of the container: {peak:.1f} GB (processes, without file cache)')
    if status != 0:
        line = f"RESULTS NOT VALID: the pipeline stopped at a step - see {p['run_dir'] / 'run.log'}"
        write_status(p['run_dir'], status='failed', finished_utc=utc_now(), verdict=line)
        say('\n' + line)
        return 1
    valid, line = verdict(p['run_dir'])
    write_status(p['run_dir'], status='finished', finished_utc=utc_now(), verdict=line)
    say('\n' + line)
    config = yaml.safe_load(p['config_text'])
    if not (config['remove_multimap'] and config['remove_antisense']):
        say(f"note: no gene table (remove_multimap / remove_antisense False); the counts are in "
            f"{p['run_dir'] / (p['run_id'] + '-readcounts.csv')}")
    say(f"report: {p['run_dir'] / 'qc' / 'report_technical.html'}")
    return 0 if valid else 1


if __name__ == '__main__':
    sys.exit(main())
