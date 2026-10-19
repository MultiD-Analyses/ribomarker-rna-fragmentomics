#!/usr/bin/env python3
"""
qcCollect.py -- collect the QC metrics of one pipeline run (and run the consistency checks).

Reads the files that the pipeline stages already write and produces
    <output>/qc/qc_metrics.json   metrics, distributions, provenance, consistency checks
                                  (qcChecks.py) and warnings
    <output>/qc/metrics_long.tsv  one row per numeric per-sample metric; no timestamps,
                                  so two runs can be compared byte for byte
    <output>/qc/consistency_checks.tsv  one row per consistency check (all statuses); the
                                  report shows only a summary unless a check is not equal
Metric ids, units and denominators come from docker/qc/metrics.yaml. Writing an id
that is not in the dictionary, or not defined for the run's mode, stops the script:
that is a programming error, not a data problem.

Rules:
  - a missing or unparseable input gives a warning and the metrics that depend on
    it are omitted, never set to 0;
  - a metric not defined for the run's mode (paired-end / single-end) is absent.

Usage (from the repository root):
    python docker/scripts/qcCollect.py --output results/levelA_pe \
        --config test/levels/work/config_levelA_pe.yaml
"""
import argparse
import ast
import csv
import json
import math
import os
import re
from pathlib import Path

import yaml

import codeFingerprint
import qcChecks

REPO = Path(__file__).resolve().parents[2]
DICTIONARY_FILE = REPO / 'docker' / 'qc' / 'metrics.yaml'
VERSION_FILE = REPO / 'docker' / 'VERSION'

SPIKE_OFFSETS = [-2, -1, 0, 1, 2]    # full-length spike reads are >= designed length - 2
FLOOR_LENGTHS = [15, 16, 17]         # assign.type.{biotype}.floor_reads
DETECTION_ABOVE = [0, 5, 10]         # counts.detected.*.gt{n}: genes with more than n counts
RRNA_BIOTYPES = {'rRNA', 'Mt_rRNA'}  # left out of counts.*_non_rrna
RBC_MIRNAS = {'human': ('Hsa-Mir-451_5p', 'Hsa-Mir-23-P3_3p'),
              'mouse': ('Mmu-Mir-451_5p', 'Mmu-Mir-23-P3_3p')}
NUMERIC_KINDS = {'count', 'rate', 'length', 'ratio'}   # written to metrics_long.tsv

# Errors that mean "this input is missing or cannot be read" (json errors are ValueErrors).
INPUT_ERRORS = (OSError, KeyError, ValueError, IndexError, TypeError, AttributeError)

# An absolute path (not part of a URL or of a relative path); group 1 = last component.
ABSOLUTE_PATH = re.compile(r'(?<![\w:/.~-])/(?:[^/\s\'"]+/)*([^/\s\'"]+)/?')


class DictionaryError(Exception):
    """The code writes an id that metrics.yaml does not define (for this mode)."""


# ----------------------------------------------------------------------------- helpers

def to_int(text):
    """A whole number from a table cell or JSON field ('12', '12.0', 12); anything else is an error."""
    if text is None:
        raise ValueError('value is missing (null)')
    number = float(text)
    if not number.is_integer():
        raise ValueError(f'not a whole number: {text!r}')
    return int(number)


def read_json(path):
    with open(path) as f:
        return json.load(f)


def read_rows(path):
    """All rows of a CSV file; every row must have as many fields as the header
    (a cut-off or damaged file raises instead of giving silently wrong numbers)."""
    with open(path, newline='') as f:
        rows = list(csv.reader(f))
    for number, row in enumerate(rows[1:], start=2):
        if row and len(row) != len(rows[0]):
            raise ValueError(f'line {number} has {len(row)} fields, the header has {len(rows[0])}')
    return rows


def sample_columns(path, sample):
    """Table with one row per item and one column per sample -> {item: whole number} for one sample."""
    rows = read_rows(path)
    column = rows[0].index(sample)                     # ValueError if the sample is not in the table
    return {row[0]: to_int(row[column]) for row in rows[1:]}


def sample_row(path, sample):
    """Table with one row per sample -> {column: text} for one sample."""
    rows = read_rows(path)
    for row in rows[1:]:
        if row[0] == sample:
            return dict(zip(rows[0], row))
    raise KeyError(f'sample {sample} not in {Path(path).name}')


def strip_paths(value):
    """Reduce absolute paths to their last component, recursively (privacy: no analysis paths in the outputs)."""
    if isinstance(value, dict):
        return {key: strip_paths(item) for key, item in value.items()}
    if isinstance(value, list):
        return [strip_paths(item) for item in value]
    if isinstance(value, str):
        return ABSOLUTE_PATH.sub(lambda match: match.group(1), value)
    return value


def histogram_median(x, y):
    """Median of the values x with counts y (mean of the two middle values when the total is even)."""
    total = sum(y)
    middle = []
    for rank in ((total - 1) // 2, total // 2):        # 0-based ranks of the middle value(s)
        cumulative = 0
        for value, count in zip(x, y):
            cumulative += count
            if cumulative > rank:
                middle.append(value)
                break
    return (middle[0] + middle[1]) / 2


def find_table(out, name):
    """Result table. The Snakefile writes <output_prefix>/<output_prefix>-<name>.csv (for
    output_prefix results/x: results/x/results/x-<name>.csv); older runs may have it directly
    in the run folder."""
    # out is the output_prefix as given in the config; an absolute prefix /a/run gives
    # /a/run//a/run-<name>.csv, i.e. /a/run/a/run-<name>.csv
    candidates = [out / f"{out.as_posix().lstrip('/')}-{name}.csv",
                  out / 'results' / f'{out.name}-{name}.csv', out / f'{out.name}-{name}.csv']
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f'{out.name}-{name}.csv not found in the run folder')


def spike_table(out, name):
    """countRiboSpike.py table: <output_prefix>/ribospike/<output_prefix>-<name>.csv."""
    candidates = [out / 'ribospike' / f"{out.as_posix().lstrip('/')}-{name}.csv",
                  out / 'ribospike' / 'results' / f'{out.name}-{name}.csv']
    for path in candidates:
        if path.exists():
            return path
    return candidates[-1]                         # missing: the reader reports it


# ----------------------------------------------------------------------------- collector

class Collector:
    """Holds the collected values; checks every id against the metric dictionary."""

    def __init__(self, dictionary, mode, samples):
        self.entries = {entry['id']: entry for entry in dictionary['metrics']}
        self.mode = mode
        self.samples = samples
        self.metrics = {sample: {} for sample in samples}        # sample -> {id: record}
        self.distributions = {sample: [] for sample in samples}
        self.provenance = {}
        self.warnings = []
        self.checks = []
        self.facts = {sample: {} for sample in samples}          # numbers read only for the checks
        self.section_kind = {}                                   # concrete id -> (section, kind)

    def entry(self, template):
        if template not in self.entries:
            raise DictionaryError(f'{template} is not in {DICTIONARY_FILE.name}')
        entry = self.entries[template]
        if self.mode not in entry['modes']:
            raise DictionaryError(f'{template} is not defined in mode {self.mode}')
        return entry

    def put(self, sample, template, value, n_reads=None, **names):
        """n_reads: reads a median, mode or deviation was computed from (Kind A rules, qcChecks.py)."""
        entry = self.entry(template)
        metric_id = template.format(**names)
        record = {'id': metric_id, 'value': value, 'unit': entry['unit']}
        if entry.get('denominator_id'):
            record['denominator_id'] = entry['denominator_id'].format(**names)
        if n_reads is not None:
            record['n'] = n_reads
        self.metrics[sample][metric_id] = record
        self.section_kind[metric_id] = (entry['section'], entry['kind'])

    def put_distribution(self, sample, template, x, y, **names):
        entry = self.entry(template)
        self.distributions[sample].append(
            {'id': template.format(**names), 'x': list(x), 'y': list(y), 'unit': entry['unit']})

    def put_run(self, template, value):
        self.entry(template)
        self.provenance[template] = strip_paths(value)

    def get(self, sample, metric_id):
        record = self.metrics[sample].get(metric_id)
        return None if record is None else record['value']

    def warn(self, sample, source, message):
        self.warnings.append({'sample': sample, 'source': source, 'message': strip_paths(message)})

    def attempt(self, sample, source, reader, *args):
        """Run one reader; if its input is missing or unreadable, record a warning instead.
        Readers parse their whole input before writing any metric, so a failure writes nothing."""
        try:
            reader(self, *args)
        except INPUT_ERRORS as error:
            self.warn(sample, source, f'{type(error).__name__}: {error}')


# ----------------------------------------------------------------------------- 5.3 trimming

def read_cutadapt_trim(c, out, sample):
    report = read_json(out / 'trimmed' / 'logs' / f'{sample}.trim.cutadapt.json')
    counts = report['read_counts']
    reads_in = to_int(counts['input'])
    r1_with_adapter = to_int(counts['read1_with_adapter'])
    too_short = to_int(counts['filtered']['too_short'])
    written = to_int(counts['output'])
    bases = report['adapters_read1'][0]['three_prime_end']['adjacent_bases']
    base_counts = [to_int(bases[base]) for base in ('A', 'C', 'G', 'T', '')]
    if c.mode == 'pe':
        r2_with_adapter = to_int(counts['read2_with_adapter'])
        c.put(sample, 'trim.pairs_in', reads_in)
        c.put(sample, 'trim.r1_with_adapter', r1_with_adapter)
        c.put(sample, 'trim.r2_with_adapter', r2_with_adapter)
        c.put(sample, 'trim.pairs_too_short', too_short)
        c.put(sample, 'trim.pairs_written', written)
    else:
        c.put(sample, 'trim.reads_in', reads_in)
        c.put(sample, 'trim.reads_with_adapter', r1_with_adapter)
        c.put(sample, 'trim.reads_too_short', too_short)
        c.put(sample, 'trim.reads_written', written)
    c.put_distribution(sample, 'trim.r1_base_before_adapter', ['A', 'C', 'G', 'T', 'none'], base_counts)


def read_cutadapt_dimer(c, out, sample):
    report = read_json(out / 'trimmed' / 'logs' / f'{sample}.dimer.cutadapt.json')
    dimers = to_int(report['read_counts']['filtered']['too_short'])
    if c.mode == 'pe':
        c.put(sample, 'trim.adapter_dimer_pairs', dimers)
    else:
        c.put(sample, 'trim.adapter_dimer_reads', dimers)


def read_lane_counts(c, out, sample):
    """Reads per lane of a sample whose lane files the Snakefile merged (logs/lanes/<sample>.csv)."""
    template = 'trim.{lane}.pairs_in' if c.mode == 'pe' else 'trim.{lane}.reads_in'
    for read, lane, reads in read_rows(out / 'logs' / 'lanes' / f'{sample}.csv')[1:]:
        if read == 'R1':
            c.put(sample, template, int(reads), lane=f'L{lane}')


def lane_samples(config):
    """Samples whose lanes were merged: those in the run command's lanes.csv (config 'lanes')."""
    if not config.get('lanes'):
        return set()
    return {row[0] for row in read_rows(Path(config['lanes']))[1:]}


def derive_short_inserts(c, sample):
    unit = 'pairs' if c.mode == 'pe' else 'reads'
    too_short = c.get(sample, f'trim.{unit}_too_short')
    dimers = c.get(sample, f'trim.adapter_dimer_{unit}')
    if too_short is not None and dimers is not None:
        c.put(sample, f'trim.short_insert_{unit}', too_short - dimers)


def read_fastp(c, out, sample):
    source = f'trimmed/logs/{sample}.fastp.json'
    report = read_json(out / source)
    result = report['filtering_result']
    merged = report['merged_and_filtered']
    reads_before = (to_int(report['read1_before_filtering']['total_reads'])
                    + to_int(report['read2_before_filtering']['total_reads']))
    passed = to_int(result['passed_filter_reads'])
    merged_reads = to_int(merged['total_reads'])
    merged_bases = to_int(merged['total_bases'])
    q30_bases = to_int(merged['q30_bases'])
    peak = to_int(report['insert_size']['peak'])
    histogram = [to_int(n) for n in report['insert_size']['histogram']]
    low_quality = to_int(result['low_quality_reads'])
    too_many_n = to_int(result['too_many_N_reads'])
    too_short = to_int(result['too_short_reads'])
    corrected = to_int(result['corrected_reads'])

    c.put(sample, 'trim.fastp_low_quality_reads', low_quality)
    c.put(sample, 'trim.fastp_too_many_n_reads', too_many_n)
    c.put(sample, 'trim.fastp_too_short_reads', too_short)
    c.put(sample, 'trim.fastp_corrected_reads', corrected)
    c.put(sample, 'trim.merged_reads', merged_reads)
    if merged_bases > 0:
        c.put(sample, 'trim.merged_q30_rate', q30_bases / merged_bases)
    c.put(sample, 'trim.insert_size_peak', peak, n_reads=sum(histogram))
    c.put_distribution(sample, 'trim.insert_size_hist', range(len(histogram)), histogram)
    c.facts[sample]['fastp_read1_in'] = to_int(report['read1_before_filtering']['total_reads'])   # check C1

    # Pairs removed by the filters vs pairs not merged: valid only if fastp removes
    # both reads of a pair together (metrics.yaml, trim.fastp_filtered_pairs).
    if passed % 2 == 0 and (reads_before - passed) % 2 == 0:
        c.put(sample, 'trim.fastp_filtered_pairs', (reads_before - passed) // 2)
        c.put(sample, 'trim.unmerged_pairs', passed // 2 - merged_reads)
    else:
        c.warn(sample, source, 'odd read counts: pairs were not filtered together; '
                               'trim.fastp_filtered_pairs and trim.unmerged_pairs omitted')


# ----------------------------------------------------------------------------- 5.4 mapping tiers

DATABASE_LINE = re.compile(r'^(\S+): (\S+) ([TF])$')        # "<tool>: <index path> <T|F>"
COMMAND_LINE = re.compile(r'^(\S+) : (.*)$')                 # "<tool> : <command>"
K_ONE = re.compile(r'(^|\s)-k\s+1(\s|$)')                    # bowtie option: report one alignment


def read_map_log(path):
    """Tier order, tools, tool versions and commands from logs/map.log (written by mapReads.py)."""
    tiers, commands, versions = [], {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            database = DATABASE_LINE.match(line)
            command = COMMAND_LINE.match(line)
            if database:
                tiers.append({'tier': os.path.basename(database.group(2)), 'tool': database.group(1),
                              'retain': database.group(3)})
            elif command:
                commands[command.group(1)] = command.group(2)
            elif line.startswith('tools: '):
                versions = dict(tool.split('=', 1) for tool in ast.literal_eval(line[len('tools: '):]))
    if not tiers:
        raise ValueError('no databases listed')
    for tier in tiers:
        tier['reports_multi'] = not K_ONE.search(commands.get(tier['tool'], ''))
    return {'tiers': tiers, 'commands': commands, 'versions': versions}


def log_number(pattern, text):
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        raise ValueError(f'line not found: {pattern}')
    return int(match.group(1))


def read_tier_log(path):
    """Summary of one Bowtie 1 or Bowtie 2 log. Unrelated lines (warnings) are ignored."""
    with open(path) as f:
        text = f.read()
    if '# reads processed:' in text:                                     # Bowtie 1
        reads_in = log_number(r'^# reads processed: (\d+)', text)
        aligned = log_number(r'^# reads with at least one reported alignment: (\d+)', text)
        unaligned = log_number(r'^# reads that failed to align: (\d+)', text)
        multi = None
    elif 'reads; of these:' in text:                                     # Bowtie 2, unpaired reads
        reads_in = log_number(r'^(\d+) reads; of these:', text)
        unaligned = log_number(r'^\s+(\d+) \([\d.]+%\) aligned 0 times', text)
        once = log_number(r'^\s+(\d+) \([\d.]+%\) aligned exactly 1 time', text)
        multi = log_number(r'^\s+(\d+) \([\d.]+%\) aligned >1 times', text)
        aligned = once + multi
    elif re.search(r'^0 reads$', text, re.MULTILINE):                   # Bowtie 2 given no reads: only
        reads_in = aligned = unaligned = multi = 0                       # '0 reads' and the overall rate
    else:                                                                # (level A sample S3, 2026-10-05)
        raise ValueError('neither a Bowtie 1 nor a Bowtie 2 summary')
    if aligned + unaligned != reads_in:
        raise ValueError(f'aligned + unaligned ({aligned} + {unaligned}) != reads in ({reads_in})')
    return {'reads_in': reads_in, 'aligned': aligned, 'unaligned': unaligned, 'multi': multi}


def collect_mapping(c, out, sample, tiers):
    parsed = {}
    for tier in tiers:
        name = tier['tier']
        source = f'alignments/logs/{sample}.{name}.log'
        try:
            parsed[name] = read_tier_log(out / source)
        except INPUT_ERRORS as error:
            c.warn(sample, source, f'{type(error).__name__}: {error}')
    for tier in tiers:
        name = tier['tier']
        if name not in parsed:
            continue
        log = parsed[name]
        c.put(sample, 'map.{tier}.reads_in', log['reads_in'], tier=name)
        c.put(sample, 'map.{tier}.aligned', log['aligned'], tier=name)
        c.put(sample, 'map.{tier}.unaligned', log['unaligned'], tier=name)
        if log['multi'] is not None and tier['reports_multi']:
            c.put(sample, 'map.{tier}.aligned_multi', log['multi'], tier=name)
    first, last = tiers[0]['tier'], tiers[-1]['tier']
    if first in parsed:
        c.put(sample, 'map.total_in', parsed[first]['reads_in'])
    if last in parsed:
        c.put(sample, 'map.final_unaligned', parsed[last]['unaligned'])


# ----------------------------------------------------------------------------- 5.5 assignment

def read_mappingstats(c, out, sample):
    row = sample_row(find_table(out, 'mappingstats'), sample)
    values = {column: to_int(row[column]) for column in ('Mapped', 'Unmapped', 'Unique', 'Multi')}
    c.put(sample, 'assign.mapped', values['Mapped'])
    c.put(sample, 'assign.unmapped', values['Unmapped'])
    c.put(sample, 'assign.unique', values['Unique'])
    c.put(sample, 'assign.multi', values['Multi'])


def read_typecounts(c, out, sample):
    counts = sample_columns(find_table(out, 'typecounts'), sample)
    kept = sum(n for biotype, n in counts.items() if biotype != 'unmapped')
    spikes = sum(n for biotype, n in counts.items() if biotype.startswith('RiboSpike'))
    # A category missing from this count table had no reads in any sample: 0 is a measurement here.
    unannotated = counts.get('unannotated', 0)
    c.put(sample, 'assign.kept', kept)
    c.put(sample, 'assign.annotated', kept - unannotated - spikes)
    for biotype, n in counts.items():
        if biotype != 'unmapped':
            c.put(sample, 'assign.type.{biotype}.reads', n, biotype=biotype)


def derive_removed_by_policy(c, sample, config):
    mapped, kept = c.get(sample, 'assign.mapped'), c.get(sample, 'assign.kept')
    if mapped is None or kept is None:
        return
    removed = mapped - kept
    c.put(sample, 'assign.removed_by_policy', removed)
    # Split by reason. With remove_multimap every Multi read is removed (antisense multimappers included);
    # the antisense reads then are the rest. mappingstats has no antisense column.
    remove_multi = config.get('remove_multimap', True)          # Snakefile defaults
    remove_antisense = config.get('remove_antisense', True)
    multi = c.get(sample, 'assign.multi')
    if multi is None:
        return
    removed_multi = multi if remove_multi else 0
    removed_antisense = removed - removed_multi if remove_antisense else 0
    if removed_antisense < 0 or (not remove_antisense and removed != removed_multi):
        c.warn(sample, 'assign.removed_by_policy', f'removed reads ({removed}) do not split into multimappers '
               f'({removed_multi}) and antisense reads; split not reported')
        return
    c.put(sample, 'assign.removed_multimap', removed_multi)
    c.put(sample, 'assign.removed_antisense', removed_antisense)


def read_typelengths(c, out, sample):
    rows = read_rows(find_table(out, 'typelengths'))
    lengths = [int(length) for length in rows[0][2:]]
    histograms = {row[1]: [to_int(n) for n in row[2:]] for row in rows[1:] if row[0] == sample}
    if not histograms:
        raise KeyError(f'sample {sample} not in typelengths')
    floor_positions = [lengths.index(length) for length in FLOOR_LENGTHS]
    for biotype, counts in histograms.items():
        if biotype == 'unmapped':
            continue
        c.put_distribution(sample, 'assign.type.{biotype}.length_hist', lengths, counts, biotype=biotype)
        if sum(counts) > 0:
            c.put(sample, 'assign.type.{biotype}.median_length', histogram_median(lengths, counts),
                  n_reads=sum(counts), biotype=biotype)
        c.put(sample, 'assign.type.{biotype}.floor_reads', sum(counts[i] for i in floor_positions),
              biotype=biotype)


def read_distinct_sequences(c, out, sample):
    counts = sample_columns(find_table(out, 'biotypedistinctseqs'), sample)
    for biotype, n in counts.items():
        if biotype != 'unmapped':
            c.put(sample, 'assign.type.{biotype}.distinct_sequences', n, biotype=biotype)


def read_sample_metrics_log(c, out, sample):
    values = dict(row[:2] for row in read_rows(out / 'metrics' / f'{sample}.log') if len(row) >= 2)
    median = float(values['median_read_length'])
    totals = {key: to_int(values[key]) for key in ('mapping_total', 'reads_total')}    # check C6
    kept, unmapped = c.get(sample, 'assign.kept'), c.get(sample, 'assign.unmapped')
    n = kept + unmapped if kept is not None and unmapped is not None else None   # the reads in the length tally
    c.put(sample, 'assign.median_read_length', median, n_reads=n)
    c.facts[sample].update(totals)


# ----------------------------------------------------------------------------- 5.6 spike-ins

SPIKE_NAME = re.compile(r'^RiboSpike-(\d+)-(T\d+)$')        # e.g. RiboSpike-20-T1


def read_spike_types(c, out, sample):
    counts = sample_columns(spike_table(out, 'ribotype'), sample)
    full_length = sum(n for name, n in counts.items() if re.fullmatch(r'RiboSpike-T\d+', name))
    fragments = counts.get('RiboSpike-Fragment', 0)     # count table: a missing row means no reads
    c.put(sample, 'spike.reads', full_length + fragments)
    c.put(sample, 'spike.fragment_reads', fragments)


def read_spike_counts(c, out, sample):
    counts = sample_columns(spike_table(out, 'ribospike'), sample)
    per_length = {}
    for name, n in counts.items():                    # first pass: parse all names before writing anything
        length, spike_type = SPIKE_NAME.match(name).groups()      # AttributeError if a name is unexpected
        per_length[int(length)] = per_length.get(int(length), 0) + n
    for name, n in counts.items():
        length, spike_type = SPIKE_NAME.match(name).groups()
        c.put(sample, 'spike.{len}_{type}.reads', n, len=length, type=spike_type)
    for length in sorted(per_length):
        c.put(sample, 'spike.{len}.reads', per_length[length], len=length)


def designed_spike_lengths(out):
    names = [row[0] for row in read_rows(spike_table(out, 'ribospike'))[1:]]
    return sorted({int(SPIKE_NAME.match(name).group(1)) for name in names})


def read_spike_lengths(c, out, sample):
    rows = read_rows(spike_table(out, 'ribolengths'))
    lengths = [int(length) for length in rows[0][2:]]
    by_length = {}                                      # observed length -> reads, summed over spike types
    for row in rows[1:]:
        if row[0] == sample and row[1].startswith('RiboSpike-T'):
            for length, n in zip(lengths, row[2:]):
                by_length[length] = by_length.get(length, 0) + to_int(n)
    if not by_length:
        raise KeyError(f'sample {sample} not in ribolengths')
    designed_lengths = designed_spike_lengths(out)
    for designed in designed_lengths:
        reads = [by_length.get(designed + offset, 0) for offset in SPIKE_OFFSETS]
        c.put_distribution(sample, 'spike.{len}.length_offset_hist', SPIKE_OFFSETS, reads, len=designed)
        if sum(reads) > 0:
            # most reads; ties go to the offset closest to 0
            modal = max(SPIKE_OFFSETS, key=lambda offset: (reads[SPIKE_OFFSETS.index(offset)], -abs(offset)))
            c.put(sample, 'spike.{len}.modal_offset', modal, n_reads=sum(reads), len=designed)
            c.put(sample, 'spike.{len}.modal_offset_reads', reads[SPIKE_OFFSETS.index(modal)], len=designed)


# Expected value of the largest of the 32 deviations |share - 0.25| (8 random positions x 4 bases)
# if the ends were uniform, in units of the SD of one share, sqrt(0.1875 / n). From a simulation
# (multinomial, 8 independent positions, 20,000 draws, seed 0; 2026-10-02): 2.31 at n = 100, 193,
# 1,000 and 2,037, 2.32 at n = 20,000.
END_BASE_CHANCE_FACTOR = 2.31


def end_base_deviations(out, sample, table):
    """{designed length: (largest |share - 0.25| over the 8 end positions and 4 bases, reads)} from a
    countRiboSpike end table; reads = the smallest number of reads at a position (all equal in practice)."""
    rows = read_rows(spike_table(out, table))
    header = rows[0]
    bases = ['A', 'C', 'G', 'T']
    totals = {}                                         # (designed length, end position) -> base counts
    for row in rows[1:]:
        record = dict(zip(header, row))
        if record['sample'] != sample:
            continue
        length = int(SPIKE_NAME.match(record['ribospike']).group(1))
        counts = totals.setdefault((length, record['end']), {base: 0 for base in bases + ['N']})
        for base in counts:
            counts[base] += to_int(record[base])
    result = {}
    for length in sorted({length for length, _ in totals}):
        deviations, position_reads = [], []
        for (spike_length, _), counts in totals.items():
            reads = sum(counts.values())
            if spike_length == length and reads > 0:
                deviations.append(max(abs(counts[base] / reads - 0.25) for base in bases))
                position_reads.append(reads)
        if deviations:
            result[length] = (max(deviations), min(position_reads))
    return result


def read_spike_ends(c, out, sample):
    """All full-length reads (kept for comparison; mixes length offset into the ends)."""
    for length, (deviation, reads) in end_base_deviations(out, sample, 'ribospikeends').items():
        c.put(sample, 'spike.{len}.end_base_deviation', deviation, n_reads=reads, len=length)


def read_spike_ends_exact(c, out, sample):
    """Reads of exactly the designed length aligned over the whole spike: their end bases are the 4N positions."""
    for length, (deviation, reads) in end_base_deviations(out, sample, 'ribospikeends_exact').items():
        c.put(sample, 'spike.{len}.exact_reads', reads, len=length)
        c.put(sample, 'spike.{len}.end_base_composition', deviation, n_reads=reads, len=length)
        c.put(sample, 'spike.{len}.end_base_chance', END_BASE_CHANCE_FACTOR * math.sqrt(0.1875 / reads),
              n_reads=reads, len=length)


# ----------------------------------------------------------------------------- 5.8 counts

def read_gene_counts(c, out, sample, species):
    counts = sample_columns(find_table(out, 'genecounts'), sample)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))   # ties: by gene name
    c.put(sample, 'counts.library_size', sum(counts.values()))
    for above in DETECTION_ABOVE:
        c.put(sample, 'counts.detected.gt{n}', sum(1 for n in counts.values() if n > above), n=above)
    c.put(sample, 'counts.top10_reads', sum(n for _, n in ranked[:10]))
    c.put(sample, 'counts.top20', [[gene, n] for gene, n in ranked[:20] if n > 0])
    if species in RBC_MIRNAS:
        mir451, mir23a = RBC_MIRNAS[species]
        # a miRNA missing from the gene table had no reads (count table): 0 is a measurement here
        reads451, reads23a = counts.get(mir451, 0), counts.get(mir23a, 0)
        if reads23a > 0:                     # with 0 miR-23a-3p reads the ratio is not evaluable: omitted
            c.put(sample, 'counts.rbc_mirna_ratio', reads451 / reads23a, n_reads=reads451 + reads23a)
    else:
        c.warn(sample, 'config: species', f'no miRNA names for species {species!r}; '
                                          'counts.rbc_mirna_ratio omitted')


def read_gene_biotypes(c, out, sample):
    counts = sample_columns(find_table(out, 'genecounts'), sample)
    rows = read_rows(find_table(out, 'genelengths'))
    column = rows[0].index('biotype')
    biotype_of = {row[0]: row[column] for row in rows[1:]}
    biotypes = sorted(set(biotype_of[gene] for gene in counts))    # KeyError if a gene has no biotype
    for biotype in biotypes:
        for above in DETECTION_ABOVE:
            detected = sum(1 for gene, n in counts.items() if biotype_of[gene] == biotype and n > above)
            c.put(sample, 'counts.detected.{biotype}.gt{n}', detected, biotype=biotype, n=above)
        ranked = sorted(((gene, n) for gene, n in counts.items() if biotype_of[gene] == biotype and n > 0),
                        key=lambda item: (-item[1], item[0]))                 # ties: by gene name
        c.put(sample, 'counts.{biotype}.top20', [[gene, n] for gene, n in ranked[:20]], biotype=biotype)
    non_rrna = sorted((n for gene, n in counts.items() if biotype_of[gene] not in RRNA_BIOTYPES), reverse=True)
    c.put(sample, 'counts.library_size_non_rrna', sum(non_rrna))
    c.put(sample, 'counts.top10_reads_non_rrna', sum(non_rrna[:10]))


# ----------------------------------------------------------------------------- 5.1 provenance

def read_version():
    """VERSION, the code fingerprint and the git commit. The fingerprint (codeFingerprint.py) is computed
    from the files when the Python sources are next to VERSION (a repository checkout), otherwise read
    from CODE_FINGERPRINT (Docker image, written at build time). The commit comes from .git if readable
    (git itself is not run), otherwise from a GIT_COMMIT file next to VERSION (Docker image)."""
    version = {'version': VERSION_FILE.read_text().strip()}
    code_dir = VERSION_FILE.parent
    if (code_dir / 'scripts' / 'qcCollect.py').exists():
        try:
            version['code_fingerprint'] = codeFingerprint.fingerprint(code_dir)
        except OSError as error:
            version['code_fingerprint'] = f'unavailable ({type(error).__name__})'
    elif (code_dir / 'CODE_FINGERPRINT').exists():
        version['code_fingerprint'] = (code_dir / 'CODE_FINGERPRINT').read_text().strip()
    head_file = REPO / '.git' / 'HEAD'
    if head_file.exists():
        head = head_file.read_text().strip()
        if head.startswith('ref: '):
            ref = head[len('ref: '):]
            ref_file = REPO / '.git' / ref
            if ref_file.exists():
                version['git_commit'] = ref_file.read_text().strip()
            elif (REPO / '.git' / 'packed-refs').exists():
                for line in (REPO / '.git' / 'packed-refs').read_text().splitlines():
                    if line.endswith(' ' + ref):
                        version['git_commit'] = line.split()[0]
        else:
            version['git_commit'] = head                    # detached HEAD
    elif (VERSION_FILE.parent / 'GIT_COMMIT').exists():
        version['git_commit'] = (VERSION_FILE.parent / 'GIT_COMMIT').read_text().strip()
    return version


def read_samples(sheet):
    """Sample sheet rows: sample, condition, FASTQ file names (without directories)."""
    rows = [row for row in read_rows(sheet) if row]
    return [[row[0], row[1]] + [os.path.basename(name) for name in row[2:]] for row in rows]


def collect_provenance(c, out, config, map_log):
    c.attempt(None, 'docker/VERSION', lambda c: c.put_run('prov.pipeline_version', read_version()))
    c.put_run('prov.config', config)
    manifest = Path(config['dbdir']) / f"MANIFEST_{config['species']}.json"
    c.attempt(None, f"<dbdir>/MANIFEST_{config['species']}.json",
              lambda c: c.put_run('prov.reference_build', read_json(manifest)))

    first = c.samples[0]
    tools = {}
    parameters = {}

    def from_cutadapt(c):
        for kind in ('trim', 'dimer'):
            report = read_json(out / 'trimmed' / 'logs' / f'{first}.{kind}.cutadapt.json')
            tools['cutadapt'] = report['cutadapt_version']
            parameters[f'cutadapt {kind}'] = ' '.join(report['command_line_arguments'])

    def from_fastp(c):
        tools['fastp'] = read_json(out / 'trimmed' / 'logs' / f'{first}.fastp.json')['summary']['fastp_version']

    c.attempt(first, f'trimmed/logs/{first}.trim/dimer.cutadapt.json', from_cutadapt)
    if c.mode == 'pe':
        c.attempt(first, f'trimmed/logs/{first}.fastp.json', from_fastp)
    if map_log is not None and map_log['versions']:
        tools.update(map_log['versions'])
    if tools:
        c.put_run('prov.tool_versions', tools)
    if parameters:
        # file names in the command line are replaced by <file>; only the options remain
        c.put_run('prov.library_parameters',
                  {kind: re.sub(r'\S*\.(fastq|fq|json)(\.gz)?(?=\s|$)', '<file>', line)
                   for kind, line in parameters.items()})


# ----------------------------------------------------------------------------- checks (qcChecks.py)

# Result tables and where they keep the sample names (check C9). Not included:
# *-ribospikeends.csv, which has rows only for samples with spike-in reads.
SAMPLE_TABLES = [
    ('*-mappingstats.csv', lambda out: find_table(out, 'mappingstats'), 'rows'),
    ('*-typecounts.csv', lambda out: find_table(out, 'typecounts'), 'columns'),
    ('*-typelengths.csv', lambda out: find_table(out, 'typelengths'), 'rows'),
    ('*-biotypedistinctseqs.csv', lambda out: find_table(out, 'biotypedistinctseqs'), 'columns'),
    ('*-genecounts.csv', lambda out: find_table(out, 'genecounts'), 'columns'),
    ('ribospike/results/*-ribotype.csv', lambda out: spike_table(out, 'ribotype'), 'columns'),
    ('ribospike/results/*-ribospike.csv', lambda out: spike_table(out, 'ribospike'), 'columns'),
    ('ribospike/results/*-ribolengths.csv', lambda out: spike_table(out, 'ribolengths'), 'rows'),
    ('samples.trimmed.csv', lambda out: out / 'samples.trimmed.csv', 'rows without header'),
]


def tables_not_made(config):
    """{table of SAMPLE_TABLES: reason} for the result tables the settings switch off."""
    skipped = {}
    if not config.get('count_ribospike', True):                              # Snakefile defaults
        skipped.update({name: qcChecks.SPIKE_OFF for name, _, _ in SAMPLE_TABLES if name.startswith('ribospike/')})
    if not gene_table_made(config):
        skipped['*-genecounts.csv'] = qcChecks.GENES_OFF
    return skipped


def gene_table_made(config):
    """The Snakefile makes the gene table only with both removal settings on (collapseGenes.py)."""
    return config.get('remove_multimap', True) and config.get('remove_antisense', True)


def table_sample_sets(out, skipped=()):
    """{table: set of sample names in it, or None if the table cannot be read} (check C9).
    Tables in skipped are not made under the run's settings and are not looked for."""
    sets = {}
    for name, path_of, layout in SAMPLE_TABLES:
        if name in skipped:
            continue
        try:
            rows = read_rows(path_of(out))
        except INPUT_ERRORS:
            sets[name] = None                   # the readers have already warned about this table
            continue
        if layout == 'columns':
            sets[name] = set(rows[0][1:])
        elif layout == 'rows':
            sets[name] = {row[0] for row in rows[1:] if row}
        else:
            sets[name] = {row[0] for row in rows if row}
    return sets


def read_bed_contigs(path):
    with open(path) as f:
        return [line.split('\t', 1)[0] for line in f if line.strip()]


def run_checks(c, out, config, map_log):
    """Consistency checks C1-C10 and the error / computability rules (registers A4, B3)."""
    tiers = map_log['tiers'] if map_log is not None else []
    spike_counted = config.get('count_ribospike', True)                      # Snakefile default
    genes_counted = gene_table_made(config)
    for sample in c.samples:
        values = {metric_id: record['value'] for metric_id, record in c.metrics[sample].items()}
        c.checks += qcChecks.sample_checks(sample, values, c.facts[sample], c.mode, tiers, spike_counted,
                                           genes_counted)
        if values.get('map.total_in') == 0:                                   # register A4
            c.warn(sample, 'map.total_in', 'no reads entered mapping; later sections are not evaluable')

    species = config['species']
    bed = f'<dbdir>/{species}.bed'
    contigs = None
    try:
        contigs = read_bed_contigs(Path(config['dbdir']) / f'{species}.bed')
    except INPUT_ERRORS as error:
        c.warn(None, bed, f'{type(error).__name__}: {error}')
    species_check = qcChecks.species_check(contigs, species)
    c.checks.append(species_check)
    c.put_run('prov.reference_species_check', {'status': species_check['status'], 'detail': species_check['subject']})
    if species_check['status'] == 'differs':                                  # register B3: a warning in version 1
        c.warn(None, bed, f"{species_check['observed']} contigs are named for another species (check C8)")

    skipped = tables_not_made(config)
    c.checks += qcChecks.sample_set_checks(table_sample_sets(out, skipped), c.samples)
    c.checks += [qcChecks.not_applicable('C9', None, f'{name}: {reason}') for name, reason in skipped.items()]


# ----------------------------------------------------------------------------- main

def warn_unsupported(c, config, samples):
    """Settings and trimming branches the QC numbers and checks are not designed for."""
    if config.get('collapse'):
        c.warn(None, 'config: collapse', 'collapse: True — the mapping tiers count distinct sequences, not reads, '
               'so the map.* metrics are not read counts and checks C2-C4 may differ for this reason alone')
    if (config.get('subsample_reads') or 0) > 0 and not config.get('collapse'):   # the Snakefile subsamples only without collapse
        c.warn(None, 'config: subsample_reads', f"subsample_reads: {config['subsample_reads']} — only a subsample of "
               'the trimmed reads was mapped, so checks C2 and C5 compare it with all trimmed reads and differ')
    if not config.get('remove_antisense', True):   # protein_coding / pre-mRNA antisense tagging gap (known limitation)
        c.warn(None, 'config: remove_antisense', 'remove_antisense: False — antisense reads of protein_coding and '
               'pre-mRNA features are not tagged _antisense: they are counted with the sense feature in '
               '*-readcounts.csv (antisense reads of the other biotypes are separate _antisense rows)')
    if config.get('species') == 'mouse':        # no mouse dataset has been checked yet
        c.warn(None, 'config: species', 'mouse references not validated on mouse data: their content was checked '
               '(no human contigs, mouse tRNA loci, synthetic reads), but no real mouse result has been compared')
    for sample in samples:
        if 'NEB' in sample and not config['paired_end']:   # the NEB branch exists only in single-end trimming
            c.warn(sample, 'sample name', 'trimming branch not supported by QC: the name contains "NEB", so '
                   'trimAdapters.py used the NEB adapter without -u 1 and wrote no cutadapt JSON (single-end)')
        if 'SMART' in sample:
            c.warn(sample, 'sample name', 'not supported by QC: the name contains "SMART", which makes '
                   'parseMappings.py flip the strand of every read')


def collect(out, config):
    dictionary = yaml.safe_load(DICTIONARY_FILE.read_text(encoding='utf-8'))
    mode = 'pe' if config['paired_end'] else 'se'
    sheet = Path(config['input_dir']) / config['samples']
    sample_rows = read_samples(sheet)                   # without a sample list there is no report
    samples = [row[0] for row in sample_rows]
    c = Collector(dictionary, mode, samples)
    c.put_run('prov.samples', sample_rows)

    warn_unsupported(c, config, samples)

    merged = lane_samples(config)
    map_log = None
    try:
        map_log = read_map_log(out / 'logs' / 'map.log')
    except INPUT_ERRORS as error:
        c.warn(None, 'logs/map.log', f'{type(error).__name__}: {error}; no tier metrics')

    for sample in samples:
        logs = f'trimmed/logs/{sample}'
        c.attempt(sample, f'{logs}.trim.cutadapt.json', read_cutadapt_trim, out, sample)
        c.attempt(sample, f'{logs}.dimer.cutadapt.json', read_cutadapt_dimer, out, sample)
        derive_short_inserts(c, sample)
        if sample in merged:
            c.attempt(sample, f'logs/lanes/{sample}.csv', read_lane_counts, out, sample)
        if mode == 'pe':
            c.attempt(sample, f'{logs}.fastp.json', read_fastp, out, sample)
        if map_log is not None:
            collect_mapping(c, out, sample, map_log['tiers'])
        c.attempt(sample, '*-mappingstats.csv', read_mappingstats, out, sample)
        c.attempt(sample, '*-typecounts.csv', read_typecounts, out, sample)
        derive_removed_by_policy(c, sample, config)
        c.attempt(sample, '*-typelengths.csv', read_typelengths, out, sample)
        c.attempt(sample, '*-biotypedistinctseqs.csv', read_distinct_sequences, out, sample)
        c.attempt(sample, f'metrics/{sample}.log', read_sample_metrics_log, out, sample)
        if config.get('count_ribospike', True):          # off: the spike-in tables are not made (no warnings)
            c.attempt(sample, 'ribospike/results/*-ribotype.csv', read_spike_types, out, sample)
            c.attempt(sample, 'ribospike/results/*-ribospike.csv', read_spike_counts, out, sample)
            c.attempt(sample, 'ribospike/results/*-ribolengths.csv', read_spike_lengths, out, sample)
            c.attempt(sample, 'ribospike/results/*-ribospikeends.csv', read_spike_ends, out, sample)
            c.attempt(sample, 'ribospike/results/*-ribospikeends_exact.csv', read_spike_ends_exact, out, sample)
        if gene_table_made(config):                       # off: no gene table, so no counts.* metrics (no warnings)
            c.attempt(sample, '*-genecounts.csv', read_gene_counts, out, sample, config['species'])
            c.attempt(sample, '*-genecounts.csv + *-genelengths.csv', read_gene_biotypes, out, sample)

    collect_provenance(c, out, config, map_log)
    run_checks(c, out, config, map_log)
    for sample in samples:                                  # register A1-A3: n, share and flag per record
        qcChecks.annotate_computability(c.metrics[sample], c.distributions[sample])
    return c


def write_outputs(c, out):
    version = c.provenance.get('prov.pipeline_version', {}).get('version', '')
    report = {
        'schema_version': '2.0',
        'run_id': out.name,
        'mode': c.mode,
        'pipeline_version': version,
        'provenance': c.provenance,
        'samples': {sample: {'metrics': list(c.metrics[sample].values()),
                             'distributions': c.distributions[sample]} for sample in c.samples},
        'check_definitions': qcChecks.CHECKS,
        'checks': c.checks,
        'warnings': c.warnings,
    }
    qc_dir = out / 'qc'
    qc_dir.mkdir(exist_ok=True)
    with open(qc_dir / 'qc_metrics.json', 'w') as f:
        json.dump(report, f, indent=1)
        f.write('\n')
    with open(qc_dir / 'metrics_long.tsv', 'w') as f:
        f.write('run_id\tsample\tsection\tmetric_id\tvalue\tunit\tpipeline_version\n')
        for sample in c.samples:
            for record in c.metrics[sample].values():
                section, kind = c.section_kind[record['id']]
                if kind in NUMERIC_KINDS:
                    f.write(f"{out.name}\t{sample}\t{section}\t{record['id']}\t{record['value']}\t"
                            f"{record['unit']}\t{version}\n")
    with open(qc_dir / 'consistency_checks.tsv', 'w') as f:
        f.write('check\tstatus\tsample\tsubject\texpected\tobserved\tdifference\tidentity\n')
        for check in c.checks:
            cells = [check['id'], check['status'], check['sample'] or 'run', check['subject'],
                     check['expected'], check['observed'], check['difference'], qcChecks.CHECKS.get(check['id'])]
            f.write('\t'.join('' if cell is None else str(cell) for cell in cells) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', required=True, help='run folder (the output_prefix of the run)')
    parser.add_argument('--config', required=True, help='the configuration file the run was made with')
    parser.add_argument('--qc-dir', help='folder with metrics.yaml (default: docker/qc of this repository; '
                                         'the Snakefile passes <Snakefile folder>/qc, also inside Docker)')
    parser.add_argument('--version-file', help='VERSION file (default: docker/VERSION; the Snakefile passes '
                                               '<Snakefile folder>/VERSION)')
    args = parser.parse_args()
    global DICTIONARY_FILE, VERSION_FILE        # inside the compiled Docker binary __file__ is not in the repository
    if args.qc_dir:
        DICTIONARY_FILE = Path(args.qc_dir) / 'metrics.yaml'
    if args.version_file:
        VERSION_FILE = Path(args.version_file)

    out = Path(args.output)
    with open(args.config) as f:
        config = yaml.safe_load(f)
    c = collect(out, config)
    write_outputs(c, out)

    n_metrics = sum(len(metrics) for metrics in c.metrics.values())
    print(f'qcCollect: {n_metrics} metrics for {len(c.samples)} samples written to {out / "qc"}; '
          f'{len(c.warnings)} warning(s)')
    for warning in c.warnings:
        print(f"  WARNING {warning['sample'] or 'run'}: {warning['source']}: {warning['message']}")
    statuses = [check['status'] for check in c.checks]
    print('checks: ' + ', '.join(f'{statuses.count(s)} {s}' for s in
                                 ('equal', 'differs', 'not evaluable', 'not applicable')))
    for check in c.checks:
        if check['status'] in ('differs', 'not evaluable'):
            print(f"  {check['status'].upper()} {check['id']} {check['sample'] or 'run'} {check['subject'] or ''}: "
                  f"expected {check['expected']}, observed {check['observed']}")


if __name__ == '__main__':
    main()
