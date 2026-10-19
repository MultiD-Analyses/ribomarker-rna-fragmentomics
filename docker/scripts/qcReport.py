#!/usr/bin/env python3
"""
qcReport.py -- technical QC report of one pipeline run.

Reads <output>/qc/qc_metrics.json (written by qcCollect.py), writes MultiQC custom-content
files to <output>/qc/mqc/ and runs MultiQC on them together with the cutadapt (trim step)
and fastp reports of each sample:
    <output>/qc/report_technical.html

All numbers in the RiboMarker sections come from qc_metrics.json: counts as collected
(the same values as metrics_long.tsv), percentages as the stored share x 100. Section
texts and column descriptions come from docker/qc/metrics.yaml.

Usage (from the repository root):
    python docker/scripts/qcReport.py --output results/levelC_pe
"""
import argparse
import html
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
DICTIONARY_FILE = REPO / 'docker' / 'qc' / 'metrics.yaml'
MULTIQC_CONFIG = REPO / 'docker' / 'qc' / 'multiqc_config.yaml'
PARENT = {'parent_id': 'ribomarker_qc', 'parent_name': 'RiboMarker QC',
          'parent_description': 'Sections computed by the pipeline\'s QC collector.'}
RRNA_BIOTYPES = {'rRNA', 'Mt_rRNA'}
ALL_TYPES = 'All RNA types'
SPIKE_OFFSETS = ['-2', '-1', '0', '1', '2']

# Colours (reference palette of the dataviz method; each set checked with its validator for
# colour-vision-deficiency and normal-vision separation of neighbouring stack segments).
# One colour per RNA entity, the same in the tier plot and the composition plot.
ENTITY_COLOURS = {'spike-ins': '#eb6834', 'rRNA': '#2a78d6', 'tRNA': '#e34948', 'miRNA': '#4a3aa7',
                  'piRNA': '#e87ba4', 'lncRNA and other ncRNA': '#008300', 'mRNA': '#eda100',
                  'genome, precursors, introns': '#1baf7a'}
NEUTRAL_DARK, NEUTRAL_LIGHT = '#52514e', '#c3c2b7'           # UniVec / unannotated; not aligned
FATE_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
OFFSET_COLOURS = ['#4a3aa7', '#eb6834', '#2a78d6', '#eda100', '#008300']   # -2 .. +2, 0 = blue
LENGTH_RAMP = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281']      # spike lengths, light -> dark

# Tiers grouped by the entity they hold (tier names are <species>-<database>, see logs/map.log).
TIER_GROUPS = [(r'univec', 'UniVec'), (r'spike', 'spike-ins'), (r'rrna', 'rRNA'), (r'trna', 'tRNA'),
               (r'mirna|star|-pre$', 'miRNA'), (r'pirna', 'piRNA'), (r'ncrna', 'lncRNA and other ncRNA'),
               (r'cdna', 'mRNA'), (r'genome', 'genome, precursors, introns')]
# RNA types (typecounts) grouped by the reference tier their features are in, so that the composition plot
# (after multimapper/antisense removal) uses the same groups as the tier plot (before removal); review
# 2026-10-07. Checked on the human assets: the features of each type lie in one reference, except 42 5S rRNA
# genes in the Ensembl ncRNA reference (type rRNA) and 'unannotated' (UniVec reads and genome reads without a
# feature; composition() moves the UniVec reads to the UniVec group). Unlisted types (lncRNA, snRNA, snoRNA,
# scaRNA, misc_RNA) are in the Ensembl ncRNA reference. Mt_rRNA and Mt_tRNA are in it too and stay with it
# for now: where to show them is an open decision.
TYPE_GROUPS = {'rRNA': 'rRNA', 'tRNA': 'tRNA', 'miRNA': 'miRNA', 'miRNA-star': 'miRNA', 'miRNA-pre': 'miRNA',
               'piRNA': 'piRNA', 'protein_coding': 'mRNA',
               'Mt_rRNA': 'lncRNA and other ncRNA', 'Mt_tRNA': 'lncRNA and other ncRNA',
               'pre-mRNA': 'genome, precursors, introns', 'pre-tRNA': 'genome, precursors, introns',
               'unannotated': 'genome, precursors, introns'}
# Groups in the order of the reference stack (logs/map.log), in both plots.
GROUP_ORDER = ['UniVec', 'spike-ins', 'rRNA', 'tRNA', 'miRNA', 'piRNA', 'lncRNA and other ncRNA', 'mRNA',
               'genome, precursors, introns']


def tier_group(tier):
    for pattern, group in TIER_GROUPS:
        if re.search(pattern, tier, re.IGNORECASE):
            return group
    return 'other tier'


def type_group(biotype):
    if biotype.startswith('RiboSpike'):
        return 'spike-ins'
    return TYPE_GROUPS.get(biotype, 'lncRNA and other ncRNA')


def colour(group):
    return ENTITY_COLOURS.get(group, NEUTRAL_DARK)


def run_tiers(run):
    """Tier names in tier order (from the collected map.{tier}.reads_in ids)."""
    return [metric_id.split('.')[1] for metric_id in run.ids(r'map\.[^.]+\.reads_in')]


def tier_groups(run, sample, tiers):
    """{group: reads aligned in its tiers} for one sample (before multimapper/antisense removal)."""
    row = {}
    for tier in tiers:
        aligned = run.value(sample, f'map.{tier}.aligned')
        if aligned is not None:
            row[tier_group(tier)] = row.get(tier_group(tier), 0) + aligned
    return row


# ----------------------------------------------------------------------------- data access

class Run:
    """qc_metrics.json of one run, with lookups by sample and metric id."""

    def __init__(self, report, dictionary):
        self.report = report
        self.mode = report['mode']
        self.samples = list(report['samples'])
        self.metrics = {s: {r['id']: r for r in d['metrics']} for s, d in report['samples'].items()}
        self.distributions = {s: {r['id']: r for r in d['distributions']} for s, d in report['samples'].items()}
        # placeholders such as {n} become named groups, so a concrete id can fill them in texts
        self.templates = [(re.compile('^' + re.sub(r'\\\{(\w+)\\\}', r'(?P<\1>[^.]+)', re.escape(e['id'])) + '$'), e)
                          for e in dictionary['metrics']]

    def value(self, sample, metric_id):
        record = self.metrics[sample].get(metric_id)
        return None if record is None else record['value']

    def percent(self, sample, metric_id):
        """share x 100, or None if the metric or its denominator is missing or 0."""
        record = self.metrics[sample].get(metric_id)
        return None if record is None or 'share' not in record else record['share'] * 100

    def entry(self, metric_id):
        """Dictionary entry of a concrete id (placeholders matched)."""
        return self.match(metric_id)[1]

    def match(self, metric_id):
        """(placeholder values, dictionary entry) of a concrete id, e.g. ({'n': '10'}, entry of counts.detected.gt{n})."""
        for pattern, entry in self.templates:
            found = pattern.match(metric_id)
            if found:
                return found.groupdict(), entry
        raise KeyError(metric_id)

    def ids(self, pattern):
        """Concrete metric ids (first sample that has any) matching a regular expression, in collector order."""
        found = []
        for sample in self.samples:
            for metric_id in self.metrics[sample]:
                if re.fullmatch(pattern, metric_id) and metric_id not in found:
                    found.append(metric_id)
        return found


def described(run, metric_id):
    """Help text of a metric: definition, denominator and caveat from metrics.yaml, placeholders filled in."""
    names, entry = run.match(metric_id)

    def filled(text):
        text = ' '.join(str(text).split())
        for name, value in names.items():
            text = text.replace('{' + name + '}', value)
        return text

    text = filled(entry['definition'])
    if entry.get('denominator_id'):
        denominator_id = filled(entry['denominator_id'])
        text += f' Percent of: {run.entry(denominator_id)["label"]} ({denominator_id}).'
    if entry.get('caveat'):
        text += ' Caveat: ' + filled(entry['caveat'])
    return text


def without_empty(data):
    """Drop samples without values (a missing metric is absent, never 0)."""
    return {sample: row for sample, row in data.items() if row}


# ----------------------------------------------------------------------------- sections

def general_stats(run):
    pe = run.mode == 'pe'
    columns = [   # (column, metric id, percent?, title, format)
        ('reads_in', 'trim.pairs_in' if pe else 'trim.reads_in', False, 'Pairs in' if pe else 'Reads in', '{:,.0f}'),
        ('short', 'trim.pairs_too_short' if pe else 'trim.reads_too_short', True,
         '% < 15 nt (of pairs in)' if pe else '% < 15 nt (of reads in)', '{:,.1f}'),
    ]
    if pe:   # titles name the denominator; the help text gives the metric ids
        columns.append(('merged', 'trim.merged_reads', True, '% merged (of pairs after cutadapt)', '{:,.1f}'))
    columns += [
        ('mapping_in', 'map.total_in', False, 'Reads to mapping', '{:,.0f}'),
        ('unaligned', 'map.final_unaligned', True, '% unaligned (of reads to mapping)', '{:,.1f}'),
        ('removed', 'assign.removed_by_policy', True, '% removed by policy (of mapped)', '{:,.1f}'),
        ('library', 'counts.library_size', False, 'Library size', '{:,.0f}'),
        ('genes10', 'counts.detected.gt10', False, 'Genes > 10', '{:,.0f}'),
    ]
    data = {s: {} for s in run.samples}
    headers = {}
    for column, metric_id, is_percent, title, fmt in columns:
        for s in run.samples:
            v = run.percent(s, metric_id) if is_percent else run.value(s, metric_id)
            if v is not None:
                data[s][column] = v
        headers[column] = {'title': title, 'description': described(run, metric_id), 'format': fmt,
                           'suffix': '%' if is_percent else '', 'scale': 'Blues'}   # one neutral scale: no traffic-light colours
        if is_percent:
            headers[column].update({'min': 0, 'max': 100})
    return {'id': 'ribomarker_general_stats', 'namespace': 'RiboMarker', 'plot_type': 'generalstats', 'headers': headers,
            'data': without_empty(data)}


def read_fate(run):
    pe = run.mode == 'pe'
    unit = 'pairs' if pe else 'reads'
    one, many = ('read pair', 'read pairs') if pe else ('read', 'reads')
    categories = [('trim.merged_reads', 'Merged: to mapping')] if pe else [('trim.reads_written', 'Written: to mapping')]
    categories += [(f'trim.adapter_dimer_{unit}', 'Adapter dimers (empty insert)'),
                   (f'trim.short_insert_{unit}', 'Inserts 1-14 nt')]
    if pe:
        categories += [('trim.fastp_filtered_pairs', 'Removed by fastp filters'),
                       ('trim.unmerged_pairs', 'Not merged')]
    data = {s: {cid: run.value(s, cid) for cid, _ in categories if run.value(s, cid) is not None}
            for s in run.samples}
    return {**PARENT, 'id': 'ribomarker_read_fate', 'section_name': 'Trimming and merging: read fate',
            'plot_type': 'bargraph',
            'categories': {cid: {'name': name, 'color': FATE_COLOURS[i]} for i, (cid, name) in enumerate(categories)},
            'pconfig': {'id': 'ribomarker_read_fate_plot', 'title': 'RiboMarker: read fate before mapping',
                        'ylab': many.capitalize(), 'sort_samples': False},
            'data': without_empty(data)}


def tier_fate(run):
    tiers = run_tiers(run)
    groups = []                                           # in tier order
    for tier in tiers:
        if tier_group(tier) not in groups:
            groups.append(tier_group(tier))
    data = {}
    for s in run.samples:
        row = tier_groups(run, s, tiers)
        if run.value(s, 'map.final_unaligned') is not None:
            row['not aligned'] = run.value(s, 'map.final_unaligned')
        data[s] = row
    categories = {g: {'name': g, 'color': NEUTRAL_DARK if g in ('UniVec', 'other tier') else colour(g)} for g in groups}
    categories['not aligned'] = {'name': 'not aligned', 'color': NEUTRAL_LIGHT}
    detail = html_table(['Sample'] + tiers + ['not aligned'],
                        [[html.escape(s)] + [fmt(run.value(s, f'map.{tier}.aligned')) for tier in tiers]
                         + [fmt(run.value(s, 'map.final_unaligned'))] for s in run.samples])
    return {**PARENT, 'id': 'ribomarker_tiers', 'section_name': 'Mapping per reference tier',
            'description': 'Reads aligned per tier, in tier order; a read aligned in one tier is not offered to '
                           'later tiers. Composition <b>before</b> multimapper/antisense removal. '
                           'Bowtie 1 tiers (UniVec, rRNA) report one alignment per read. The three miRNA tiers '
                           '(mature, star, precursor) are shown together. "lncRNA and other ncRNA" is the Ensembl '
                           'ncRNA tier (lncRNA, snRNA, snoRNA, scaRNA, misc_RNA, Mt_rRNA, Mt_tRNA); "genome, '
                           'precursors, introns" is the genome tier (pre-mRNA, pre-tRNA and reads without '
                           'annotation). The composition plot below uses the same groups.'
                           '<details><summary>Reads aligned in every tier</summary>' + detail + '</details>',
            'plot_type': 'bargraph', 'categories': categories,
            'pconfig': {'id': 'ribomarker_tiers_plot', 'title': 'RiboMarker: reads aligned per tier',
                        'ylab': 'Reads', 'sort_samples': False},
            'data': without_empty(data)}


def composition(run):
    biotypes = [metric_id.split('.')[2] for metric_id in run.ids(r'assign\.type\.[^.]+\.reads')]
    tiers = run_tiers(run)
    univec_tiers = [t for t in tiers if tier_group(t) == 'UniVec']
    grouped, not_split, larger = {}, [], []
    for s in run.samples:
        row = {}
        for b in biotypes:
            n = run.value(s, f'assign.type.{b}.reads')
            if n is not None:
                row[type_group(b)] = row.get(type_group(b), 0) + n
        # The pipeline counts UniVec reads as 'unannotated' (UniVec contigs have no features). They are never
        # removed (Bowtie 1 with -k 1: no multimappers; no features: no antisense), so the kept UniVec reads are
        # the reads aligned in the UniVec tier. If that does not fit, 'unannotated' stays in the genome group.
        univec = [run.value(s, f'map.{t}.aligned') for t in univec_tiers]
        unannotated = run.value(s, 'assign.type.unannotated.reads') or 0
        if univec and (None in univec or sum(univec) > unannotated):
            not_split.append(s)
        elif univec:
            row['UniVec'] = sum(univec)
            row['genome, precursors, introns'] = row.get('genome, precursors, introns', 0) - sum(univec)
        before = tier_groups(run, s, tiers)       # a group should not gain reads by the removal
        larger += [f'{s}: {g}' for g, n in row.items() if n > before.get(g, 0)]
        grouped[s] = row
    tabs = [('All RNA types', set()), ('Without spike-ins', {'spike-ins'}),
            ('Without rRNA and spike-ins', {'rRNA', 'spike-ins'})]       # (tab, groups left out)
    categories = {g: {'name': g, 'color': NEUTRAL_DARK if g == 'UniVec' else colour(g)} for g in GROUP_ORDER}
    detail = html_table(['Sample'] + biotypes,
                        [[html.escape(s)] + [fmt(run.value(s, f'assign.type.{b}.reads')) for b in biotypes]
                         for s in run.samples])
    notes = ''
    if not_split:
        notes += ('<p>Not split, because the UniVec tier holds more reads than "unannotated" (or is missing): '
                  + html.escape(', '.join(not_split)) + '. There all "unannotated" reads are in the genome group.</p>')
    if larger:
        notes += ('<p>Groups with more reads here than aligned in their tiers (an RNA type also comes from another '
                  'tier): ' + html.escape('; '.join(larger)) + '.</p>')
    return {**PARENT, 'id': 'ribomarker_composition', 'section_name': 'Composition by RNA type',
            'description': 'Kept reads per RNA type <b>after</b> multimapper/antisense removal, in the groups of '
                           'the tier plot above: each RNA type is shown with the tier its reads align in, so the '
                           'difference between the two plots is the reads removed in that group. The pipeline '
                           'counts UniVec reads as "unannotated"; here they are shown as UniVec (the reads aligned '
                           'in the UniVec tier) and the rest of "unannotated" is in the genome group.'
                           + notes +
                           '<details><summary>Reads of every RNA type</summary>' + detail + '</details>',
            'plot_type': 'bargraph',
            'categories': [{g: c for g, c in categories.items() if g not in left_out} for _, left_out in tabs],
            'pconfig': {'id': 'ribomarker_composition_plot', 'title': 'RiboMarker: composition by RNA type',
                        'ylab': 'Reads', 'sort_samples': False,
                        'data_labels': [{'name': tab, 'ylab': 'Reads'} for tab, _ in tabs]},
            'data': [without_empty({s: {g: n for g, n in row.items() if g not in left_out}
                                    for s, row in grouped.items()}) for _, left_out in tabs]}


def lengths(run):
    biotypes = [metric_id.split('.')[2] for metric_id in run.ids(r'assign\.type\.[^.]+\.reads')]
    totals = {b: sum((run.value(s, f'assign.type.{b}.reads') or 0) for s in run.samples) for b in biotypes}
    biotypes = sorted((b for b in biotypes if totals[b] > 0 and not b.startswith('RiboSpike')),
                      key=lambda b: -totals[b])
    datasets, labels, longest = [], [], 15
    pooled = {}                                   # 'All RNA types': the counts of the types below, summed per length
    for s in run.samples:
        counts = {}
        for b in biotypes:
            dist = run.distributions[s].get(f'assign.type.{b}.length_hist')
            for x, y in zip(dist['x'], dist['y']) if dist else ():
                counts[x] = counts.get(x, 0) + y
        n = sum(counts.values())
        if n > 0:
            pooled[s] = {x: 100 * y / n for x, y in counts.items() if x >= 15}
    if pooled:
        datasets.append(pooled)
        labels.append({'name': ALL_TYPES, 'ylab': '% of reads', 'xlab': 'Read length (nt)'})
    for b in biotypes:
        dataset = {}
        for s in run.samples:
            dist = run.distributions[s].get(f'assign.type.{b}.length_hist')
            if dist and dist['n'] > 0:
                dataset[s] = {x: 100 * y / dist['n'] for x, y in zip(dist['x'], dist['y']) if x >= 15}
                longest = max([longest] + [x for x, share in dataset[s].items() if share >= 0.1])
        if dataset:
            datasets.append(dataset)
            labels.append({'name': b, 'ylab': '% of reads of this type', 'xlab': 'Read length (nt)'})
    return {**PARENT, 'id': 'ribomarker_lengths', 'section_name': 'Read length by RNA type',
            'description': 'Share of the kept reads of one RNA type at each length (after multimapper/antisense '
                           'removal); lengths are integers, no smoothing. First tab: all RNA types together '
                           '(the kept reads of the types listed, without spike-ins and unmapped reads); then one '
                           'tab per RNA type. Types with fewer '
                           'than 100 reads in a sample are noisy here. Spike-in lengths '
                           'are in the spike-in sections.',
            'plot_type': 'linegraph',
            'pconfig': {'id': 'ribomarker_lengths_plot', 'title': 'RiboMarker: read length by RNA type',
                        'ylab': '% of reads of this type', 'xlab': 'Read length (nt)', 'ymin': 0,
                        'xmin': 15, 'xmax': longest + 5,
                        'smooth_points': None, 'data_labels': labels,
                        'showlegend': True},             # click a sample in the legend to hide or show it
            'data': datasets}


def multimapping(run):
    """What happened to the mapped reads (kept / removed as multimappers / removed as antisense), with the
    share of reads with a second (not necessarily equally good) alignment per tier in a table below."""
    parts = [('assign.kept', 'kept'), ('assign.removed_multimap', 'removed: multimapper'),
             ('assign.removed_antisense', 'removed: antisense, maps once')]
    data = {s: {name: run.value(s, metric_id) for metric_id, name in parts if run.value(s, metric_id) is not None}
            for s in run.samples}
    data = {s: row for s, row in data.items() if len(row) == len(parts)}
    if not data:
        return None
    categories = {'kept': {'name': 'kept', 'color': NEUTRAL_LIGHT},
                  'removed: multimapper': {'name': 'removed: multimapper', 'color': FATE_COLOURS[0]},
                  'removed: antisense, maps once': {'name': 'removed: antisense, maps once', 'color': FATE_COLOURS[1]}}
    tiers = [m.split('.')[1] for m in run.ids(r'map\.[^.]+\.aligned_multi')]
    rows = []
    for s in run.samples:
        row = [html.escape(s)]
        for tier in tiers:
            share = run.percent(s, f'map.{tier}.aligned_multi')
            row.append(flagged(run, s, f'map.{tier}.aligned_multi', '—' if share is None else f'{share:.1f} %'))
        rows.append(row)
    table = html_table(['Sample'] + [f'{t}' for t in tiers], rows)
    return {**PARENT, 'id': 'ribomarker_multimapping', 'section_name': 'Multimappers and antisense reads',
            'description': 'UniVec, rRNA '
                           'and spike-in tiers report one alignment per read, so their multimapping is not visible. '
                           '<details><summary>Reads with two or more reported alignments, per tier (% of the reads '
                           'aligned in the tier; Bowtie 2 "aligned &gt;1 times")</summary>This counts every read with a '
                           'second alignment, also when that alignment scores worse, so it can be much higher than the '
                           'removed multimappers above, which need equally good alignments. It shows which tiers hold '
                           'multi-copy or similar sequences.' + table + '</details>',
            'plot_type': 'bargraph', 'categories': categories,
            'pconfig': {'id': 'ribomarker_multimapping_plot', 'title': 'RiboMarker: multimappers and antisense reads',
                        'ylab': 'Mapped reads', 'sort_samples': False},
            'data': data}


def spike_offsets(run):
    lengths_ = sorted({int(m.split('.')[1]) for m in run.ids(r'spike\.\d+\.reads')})
    datasets, labels = [], []
    for length in lengths_:
        dataset = {}
        for s in run.samples:
            dist = run.distributions[s].get(f'spike.{length}.length_offset_hist')
            if dist and dist['n'] > 0:
                dataset[s] = {str(x): y for x, y in zip(dist['x'], dist['y'])}
        if dataset:
            datasets.append(dataset)
            labels.append({'name': f'{length} nt', 'ylab': 'Full-length spike reads'})
    if not datasets:
        return None
    return {**PARENT, 'id': 'ribomarker_spike_offsets', 'section_name': 'Spike-ins: length fidelity',
            'description': 'Observed minus designed length of the full-length spike-in reads (designed lengths '
                           'include the 4N ends), one tab per designed length. Expected at 0; shorter reads '
                           '(< -2) are counted as fragments.',
            'plot_type': 'bargraph',
            'categories': [{o: {'name': f'{int(o):+d} nt' if o != '0' else '0 nt', 'color': OFFSET_COLOURS[i]}
                            for i, o in enumerate(SPIKE_OFFSETS)}] * len(datasets),
            'pconfig': {'id': 'ribomarker_spike_offsets_plot', 'title': 'RiboMarker: spike-in length offset',
                        'ylab': 'Reads', 'sort_samples': False, 'cpswitch_c_active': False,
                        'data_labels': labels},
            'data': datasets}


def spike_recovery(run):
    ids = run.ids(r'spike\.\d+\.reads')
    data = {s: {m.split('.')[1] + ' nt': run.value(s, m) for m in ids if run.value(s, m) is not None}
            for s in run.samples}
    data = {s: row for s, row in data.items() if row and sum(row.values()) > 0}
    if not data:
        return None
    names = [m.split('.')[1] + ' nt' for m in ids]
    categories = {name: {'name': name, 'color': LENGTH_RAMP[i % len(LENGTH_RAMP)]} for i, name in enumerate(names)}
    return {**PARENT, 'id': 'ribomarker_spike_recovery', 'section_name': 'Spike-ins: recovery by length',
            'description': 'Full-length spike-in reads per designed length (all four types). Reads as capture '
                           'efficiency only if the spikes are added in equal amounts.',
            'plot_type': 'bargraph', 'categories': categories,
            'pconfig': {'id': 'ribomarker_spike_recovery_plot', 'title': 'RiboMarker: spike-in reads by length',
                        'ylab': 'Reads', 'sort_samples': False, 'cpswitch_c_active': False},
            'data': data}


def fmt(value, digits=0):
    if value is None:
        return '—'
    if isinstance(value, float) and digits:
        return f'{value:,.{digits}f}'
    return f'{value:,.0f}' if isinstance(value, (int, float)) else html.escape(str(value))


def flagged(run, sample, metric_id, text):
    record = run.metrics[sample].get(metric_id)
    if record and record.get('flag'):
        return f'{text} <span class="label label-warning">{html.escape(record["flag"])}</span>'
    return text


def html_table(headers, rows):
    head = ''.join(f'<th>{html.escape(h)}</th>' for h in headers)
    body = ''.join('<tr>' + ''.join(f'<td>{cell}</td>' for cell in row) + '</tr>' for row in rows)
    return f'<table class="table table-condensed table-hover"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def lane_table(run):
    """Reads per lane of the samples whose lane files were merged; no section without lanes."""
    ids = run.ids(r'trim\.L\d+\.(pairs|reads)_in')
    if not ids:
        return None
    headers = ['Sample'] + [metric_id.split('.')[1] for metric_id in ids]
    rows = [[html.escape(s)] + [fmt(run.value(s, metric_id)) for metric_id in ids] for s in run.samples]
    return {**PARENT, 'id': 'ribomarker_lane_table', 'section_name': 'Reads per lane',
            'description': 'Reads (paired-end: read pairs) of each sequencing lane, for the samples whose lane '
                           'files were merged before trimming; a sample\'s lanes add up to its reads in.',
            'plot_type': 'html', 'data': html_table(headers, rows)}


def spike_table(run):
    lengths_ = sorted({int(m.split('.')[1]) for m in run.ids(r'spike\.\d+\.reads')})
    rows = []
    for s in run.samples:
        row = [html.escape(s), fmt(run.value(s, 'spike.reads')),
               fmt(run.percent(s, 'spike.reads'), 2), fmt(run.percent(s, 'spike.fragment_reads'), 1)]
        for length in lengths_:
            at_mode = run.percent(s, f'spike.{length}.modal_offset_reads')
            text = '—' if at_mode is None else f'{run.value(s, f"spike.{length}.modal_offset")} ({at_mode:.0f} %)'
            row.append(flagged(run, s, f'spike.{length}.modal_offset_reads', text))
        for length in lengths_:
            composition, chance = (run.value(s, f'spike.{length}.end_base_{m}') for m in ('composition', 'chance'))
            text = '—' if composition is None else f'{composition:.3f} ({chance:.3f})'
            row.append(flagged(run, s, f'spike.{length}.end_base_composition', text))
        rows.append(row)
    headers = (['Sample', 'Spike-in reads', '% of reads into mapping', '% fragments (of spike-in reads)']
               + [f'Modal offset {n} nt (% of full-length reads at mode)' for n in lengths_]
               + [f'End-base composition {n} nt (chance level)' for n in lengths_])
    return {**PARENT, 'id': 'ribomarker_spike_table', 'section_name': 'Spike-ins: summary',
            'description': 'Modal length offset per designed length with the share of full-length reads at it. '
                           'End-base composition: in reads of exactly the designed length, the largest deviation '
                           'of any base share from 0.25 at the four random positions of each end; in brackets the '
                           'value expected if those bases were uniform, for the same number of reads.',
            'plot_type': 'html', 'data': html_table(headers, rows)}


def counts_table(run):
    rows = []
    for s in run.samples:
        mir = run.metrics[s].get('counts.rbc_mirna_ratio')
        ratio = '—' if mir is None else flagged(run, s, 'counts.rbc_mirna_ratio',
                                                f'{mir["value"]:.2f} (n = {mir.get("n", "—")})')
        rows.append([html.escape(s), fmt(run.value(s, 'counts.library_size')),
                     fmt(run.value(s, 'counts.library_size_non_rrna')),
                     fmt(run.value(s, 'counts.detected.gt0')), fmt(run.value(s, 'counts.detected.gt5')),
                     fmt(run.value(s, 'counts.detected.gt10')),
                     fmt(run.percent(s, 'counts.top10_reads'), 1), fmt(run.percent(s, 'counts.top10_reads_non_rrna'), 1),
                     ratio])
    headers = ['Sample', 'Library size', 'Without rRNA', 'Genes > 0', 'Genes > 5', 'Genes > 10',
               '% in top-10 genes', '% in top-10 genes without rRNA', 'miR-451a / miR-23a-3p (reads)']
    biotypes = sorted({m.split('.')[1] for m in run.ids(r'counts\.[^.]+\.top20')})
    per_type = ''.join(f'<details><summary>{html.escape(b)}</summary>{top_table(run, f"counts.{b}.top20")}</details>'
                       for b in biotypes)
    content = (html_table(headers, rows)
               + '<details><summary>20 most abundant RNAs per sample (reads)</summary>'
               + top_table(run, 'counts.top20') + '</details>'
               + '<details><summary>20 most abundant RNAs per sample, by RNA type (reads)</summary>'
               + per_type + '</details>')
    return {**PARENT, 'id': 'ribomarker_counts', 'section_name': 'Detection and counts',
            'description': 'Gene-level counts of every RNA type (mature exonic, after multimapper/antisense '
                           'removal). The miRNA ratio is a plain ratio of reads, descriptive only.',
            'plot_type': 'html', 'data': content}


def top_table(run, metric_id):
    """Rank x sample table of a top-20 list metric ('gene (reads)')."""
    rows = []
    for rank in range(20):
        row = [str(rank + 1)]
        for s in run.samples:
            top = run.value(s, metric_id) or []
            row.append(f'{html.escape(top[rank][0])} ({top[rank][1]:,})' if rank < len(top) else '')
        rows.append(row)
    return html_table(['Rank'] + run.samples, rows)


def overview(run):
    report = run.report
    prov = report['provenance']
    version = prov.get('prov.pipeline_version', {})
    statuses = [c['status'] for c in report.get('checks', [])]
    summary = ', '.join(f'{statuses.count(s)} {s}' for s in ('equal', 'differs', 'not evaluable', 'not applicable')
                        if statuses.count(s))
    problems = [c for c in report.get('checks', []) if c['status'] in ('differs', 'not evaluable')]
    parts = [html_table(['Item', 'Value'], [
        ['Run', html.escape(report['run_id'])],
        ['Mode', 'paired-end' if run.mode == 'pe' else 'single-end'],
        ['Pipeline version', html.escape(pipeline_label(version))],
        ['Samples', str(len(run.samples))],
        ['Consistency checks', html.escape(summary) + ' (every check: <code>qc/consistency_checks.tsv</code>)'],
        ['Warnings', str(len(report['warnings']))],
    ])]
    if problems:   # the only place checks are listed in the report; all of them are in qc/consistency_checks.tsv
        parts.append('<p><b>Checks that are not equal:</b></p>' + checks_table(report, problems))
    if report['warnings']:
        parts.append('<p><b>Warnings</b> (missing or unreadable inputs; the affected metrics are omitted):</p>'
                     + html_table(['Sample', 'Source', 'Message'],
                                  [[html.escape(w['sample'] or 'run'), html.escape(w['source']),
                                    html.escape(w['message'])] for w in report['warnings']]))
    return {**PARENT, 'id': 'ribomarker_overview', 'section_name': 'Run overview',
            'description': 'Run identity, consistency checks and warnings. In plots with one line per sample, '
                           'hover over a line to see its sample, or use the toolbox on the right (Highlight '
                           'Samples, Show / Hide Samples); in RiboMarker line plots, click a sample in the legend '
                           'to hide or show it.', 'plot_type': 'html',
            'data': ''.join(parts)}


def checks_table(report, checks):
    definitions = report.get('check_definitions', {})
    rows = [[html.escape(c['id']), html.escape(c['status']), html.escape(c['sample'] or 'run'),
             html.escape(c['subject'] or ''), html.escape(str(c['expected'])), html.escape(str(c['observed'])),
             html.escape(definitions.get(c['id'], ''))] for c in checks]
    return html_table(['Check', 'Status', 'Sample', 'Subject', 'Expected', 'Observed', 'Identity'], rows)


def provenance(run):
    prov = run.report['provenance']
    rows = []
    for key in ('prov.library_parameters', 'prov.reference_species_check'):   # tool versions: software_versions()
        value = prov.get(key)
        if isinstance(value, dict):
            value = '<br>'.join(f'{html.escape(str(k))}: {html.escape(str(v))}' for k, v in value.items())
        rows.append([html.escape(run.entry(key)['label']), value if value is not None else '—'])
    build = prov.get('prov.reference_build')
    if build:
        sources = build.get('sources', {})
        text = (f"{html.escape(str(build.get('species')))}, build {html.escape(str(build.get('build')))}, "
                f"built {html.escape(str(build.get('built_utc')))}<br>"
                + '<br>'.join(f'{html.escape(k)}: {html.escape(str(v.get("version")))}' for k, v in sources.items()))
        rows.append(['Reference build', text])
    else:
        rows.append(['Reference build', '— (no MANIFEST found)'])
    config = prov.get('prov.config', {})
    rows.append(['Configuration', '<br>'.join(f'{html.escape(str(k))}: {html.escape(str(v))}'
                                               for k, v in config.items())])
    # a section of its own (no parent), placed at the end just above MultiQC's Software Versions
    return {'id': 'ribomarker_provenance', 'section_name': 'Provenance',
            'description': 'Parameters and references used; the tool versions are listed under Software Versions '
                           'below. Absolute paths are reduced to their last component.', 'plot_type': 'html',
            'data': html_table(['Item', 'Value'], rows)}


def pipeline_label(version):
    """'0.1.0 (git 24c11cf8cc3f, code 1a2b3c4d5e6f)': version, commit and code fingerprint (codeFingerprint.py),
    each shortened to 12 characters; the fingerprint identifies the code also when the commit does not."""
    label = str(version.get('version', '—'))
    details = [f"{name} {str(version[key])[:12]}" for key, name in (('git_commit', 'git'), ('code_fingerprint', 'code'))
               if version.get(key)]
    return label + (f" ({', '.join(details)})" if details else '')


def software_versions(run):
    """{group: {software: version}} for MultiQC's Software Versions section (*_mqc_versions.yml)."""
    prov = run.report['provenance']
    pipeline = pipeline_label(prov.get('prov.pipeline_version', {}))
    # cutadapt and fastp are left out: MultiQC lists them itself, from their own reports
    tools = {str(k): str(v) for k, v in (prov.get('prov.tool_versions') or {}).items() if k not in ('cutadapt', 'fastp')}
    return {'RiboMarker pipeline': {'RiboMarker pipeline': pipeline, **tools}}


SECTIONS = [overview, general_stats, lane_table, read_fate, tier_fate, composition, multimapping, lengths,
            spike_table, spike_offsets, spike_recovery, counts_table, provenance]


# ----------------------------------------------------------------------------- main

def write_custom_content(run, mqc_dir):
    mqc_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for section in SECTIONS:
        content = section(run)
        if content is None:
            continue
        if 'pconfig' in content:   # no 'N bars' / 'N lines' subtitle under the plot titles
            content['pconfig']['series_label'] = False
        path = mqc_dir / f"{content['id']}_mqc.yaml"
        path.write_text(yaml.safe_dump(content, sort_keys=False, allow_unicode=True))
        written.append(path)
    versions = mqc_dir / 'ribomarker_mqc_versions.yml'
    versions.write_text(yaml.safe_dump(software_versions(run), sort_keys=False, allow_unicode=True))
    written.append(versions)
    return written


def tool_reports(out, samples, mode):
    """The trim-step cutadapt report and (paired-end) fastp report of every sample that has them."""
    logs = out / 'trimmed' / 'logs'
    names = [f'{s}.trim.cutadapt.json' for s in samples] + ([f'{s}.fastp.json' for s in samples] if mode == 'pe' else [])
    return [logs / name for name in names if (logs / name).exists()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', required=True, help='run folder that contains qc/qc_metrics.json')
    parser.add_argument('--qc-dir', help='folder with metrics.yaml and multiqc_config.yaml (default: docker/qc of '
                                         'this repository; the Snakefile passes <Snakefile folder>/qc)')
    parser.add_argument('--title', help='report title (default: the run id; the run command sets it)')
    args = parser.parse_args()
    global DICTIONARY_FILE, MULTIQC_CONFIG      # inside the compiled Docker binary __file__ is not in the repository
    if args.qc_dir:
        DICTIONARY_FILE = Path(args.qc_dir) / 'metrics.yaml'
        MULTIQC_CONFIG = Path(args.qc_dir) / 'multiqc_config.yaml'

    out = Path(args.output)
    qc_dir = out / 'qc'
    report = json.loads((qc_dir / 'qc_metrics.json').read_text())
    dictionary = yaml.safe_load(DICTIONARY_FILE.read_text(encoding='utf-8'))
    run = Run(report, dictionary)

    custom = write_custom_content(run, qc_dir / 'mqc')
    inputs = custom + tool_reports(out, run.samples, run.mode)
    # the multiqc command of the environment; python -m multiqc if it is not on PATH
    multiqc = [shutil.which('multiqc')] if shutil.which('multiqc') else [sys.executable, '-m', 'multiqc']
    command = multiqc + ['--force', '--quiet', '--no-data-dir',
               '--config', str(MULTIQC_CONFIG), '--title', args.title or report['run_id'],
               '--cl-config', 'custom_logo: ' + json.dumps(str((MULTIQC_CONFIG.parent / 'logos' / 'report_logos.png').resolve())),
               '--outdir', str(qc_dir), '--filename', 'report_technical.html'] + [str(p) for p in inputs]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f'MultiQC failed:\n{result.stderr}')
    size = (qc_dir / 'report_technical.html').stat().st_size
    print(f'qcReport: {qc_dir / "report_technical.html"} ({size / 1e6:.1f} MB; {sum(p.name.endswith("_mqc.yaml") for p in custom)} RiboMarker sections, '
          f'{len(inputs) - len(custom)} tool reports)')


if __name__ == '__main__':
    main()
