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
import base64
import html
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

import refTables
from qcCollect import RBC_MIRNAS as RBC_MIRNA_NAMES      # (miR-451a, miR-23a-3p) as MirGeneDB names per species

REPO = Path(__file__).resolve().parents[2]
DICTIONARY_FILE = REPO / 'docker' / 'qc' / 'metrics.yaml'
MULTIQC_CONFIG = REPO / 'docker' / 'qc' / 'multiqc_config.yaml'      # report.css is read from the same folder
PARENT = {'parent_id': 'ribomarker_qc', 'parent_name': 'RiboMarker QC',
          'parent_description': 'Sections computed by the pipeline\'s QC collector.',
          # the pipeline (Zenodo, all versions) and the method (Clark et al., RNA 2026): MultiQC lists them under
          # Tool Citations (toolbox); the full citations as in the README are under Provenance
          'doi': ['10.5281/zenodo.23475622', '10.1261/rna.080861.125'],
          'section_href': 'https://github.com/MultiD-Analyses/ribomarker-rna-fragmentomics'}   # like the tool sections
CITATIONS = [
    'Smelik M. RiboMarker: a pipeline for RNA fragmentomics, version {version}. MultiD Analyses AB; 2026. '
    'https://doi.org/10.5281/zenodo.23475622 (code: https://github.com/MultiD-Analyses/ribomarker-rna-fragmentomics)',
    'Clark RC, Manning AC, Howard JM, Barberán-Soler S, Kazakov SA. A method of comprehensive sequencing analysis of '
    'the small RNA fragmentome (RiboMarker). RNA 2026;32(7):1156–1181. https://doi.org/10.1261/rna.080861.125',
]
REPORT_JS = 'report.js'          # next to multiqc_config.yaml; settings above the plots, inlined into the guide
ALL_CLASSES = 'All RNA classes'
OTHER_RNAS = 'Other RNAs'
SPIKE_OFFSETS = ['-2', '-1', '0', '1', '2']

# Colours (reference palette of the dataviz method; each set checked with its validator for
# colour-vision-deficiency and normal-vision separation of neighbouring stack segments).
# One colour per RNA entity, the same in the tier plot and the composition plot.
# The two genome groups (RealSeq review #9, 0.3.0): the genome tier holds the reads of pre-mRNA and pre-tRNA and
# the reads without an annotated feature.
GENOME_PRE, GENOME_UNANNOTATED = 'genome: pre-mRNA, pre-tRNA', 'genome (unannotated)'
ENTITY_COLOURS = {'spike-ins': '#eb6834', 'rRNA': '#2a78d6', 'tRNA': '#e34948', 'miRNA': '#4a3aa7',
                  'piRNA': '#e87ba4', 'lncRNA and other ncRNA': '#008300', 'protein_coding': '#eda100',
                  GENOME_PRE: '#1baf7a',
                  # a neutral: no annotation. Checked with the dataviz validator (2026-10-26): normal-vision
                  # distance 16.0 to its neighbour GENOME_PRE, more to 'unmapped', 15.0 to UniVec (not adjacent)
                  GENOME_UNANNOTATED: '#7d7c75'}
NEUTRAL_DARK, NEUTRAL_LIGHT = '#52514e', '#c3c2b7'           # UniVec and unknown groups; unmapped
# Sample colours by report_group: the categorical slots of the dataviz reference palette in fixed order; with more
# than 8 groups the samples keep MultiQC's colours (a 9th hue would not be distinguishable).
SAMPLE_GROUP_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
FATE_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
OFFSET_COLOURS = ['#4a3aa7', '#eb6834', '#2a78d6', '#eda100', '#008300']   # -2 .. +2, 0 = blue
LENGTH_RAMP = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281']      # spike lengths, light -> dark

# Tiers grouped by the entity they hold (tier names are <species>-<database>, see logs/map.log).
TIER_GROUPS = [(r'univec', 'UniVec'), (r'spike', 'spike-ins'), (r'rrna', 'rRNA'), (r'trna', 'tRNA'),
               (r'mirna|star|-pre$', 'miRNA'), (r'pirna', 'piRNA'), (r'ncrna', 'lncRNA and other ncRNA'),
               (r'cdna', 'protein_coding'), (r'genome', 'genome')]      # 'genome': split by tier_groups()
# RNA classes grouped by the reference tier their features are in, so that the composition plot (after
# multimapper/antisense removal) uses the same groups as the tier plot (before removal); review 2026-10-07.
# Since 0.3.0 the groups come from the class table (column group; refTables.py), which keeps the grouping
# described here. Checked on the human assets: the features of each class lie in one reference, except 42
# Ensembl rRNA genes (36 5S, 6 5.8S) in the Ensembl ncRNA reference (class rRNA). UniVec reads are their own
# class since 0.3.0 (3g; group UniVec); the genome tier is split into its two groups (#9). Unlisted classes
# (lncRNA, snRNA, snoRNA, scaRNA, misc_RNA) are in the Ensembl ncRNA reference. Mt_rRNA and Mt_tRNA are in it
# too and stay with it for now: where to show them waits for RealSeq's answer.
# Groups in the order of the reference stack (logs/map.log), in both plots; groups of an own class table that
# are not listed here follow them.
GROUP_ORDER = ['UniVec', 'spike-ins', 'rRNA', 'tRNA', 'miRNA', 'piRNA', 'lncRNA and other ncRNA', 'protein_coding',
               GENOME_PRE, GENOME_UNANNOTATED]
# Group names of class tables copied into runs before the report changes of 0.3.0 (S4)
OLD_GROUPS = {'mRNA': 'protein_coding', 'genome, precursors, introns': GENOME_PRE,
              'not aligned': 'unmapped'}   # S5 review (Martin 2026-10-26): one word for reads aligned in no tier
# Read-length plots: classes with a button of their own (RealSeq review #13, Martin 2026-10-26); the reads of
# every other natural class are pooled under "Other RNAs"; genome (unannotated) comes last (#14).
OWN_LENGTH_BUTTON = {'rRNA', 'Mt_rRNA', 'tRNA', 'Mt_tRNA', 'miRNA', 'pre-miRNA', 'piRNA', 'snRNA', 'snoRNA', 'lncRNA',
                     'protein_coding', 'pre-mRNA', 'pre-tRNA'}


def terms(lead, items, after=''):
    """A description in the report's one structure (Martin 2026-10-26): a lead sentence, then one term per line
    (the words of the glossary in "How to read this report")."""
    return (f'<p>{lead}</p><ul class="rm-terms">' + ''.join(f'<li><b>{t}</b>: {text}</li>' for t, text in items)
            + '</ul>' + after)


GLOSSARY = [
    ('RNA class', 'rRNA, tRNA, miRNA, ...: what a read is counted as (the class table, in the section Tiers and RNA '
                  'classes)'),
    ('Tier', 'one reference in the mapping order (UniVec, spike-ins, rRNA, tRNA, miRNA, piRNA, Ensembl ncRNA, Ensembl '
             'cDNA, genome); a read aligned in one tier is not offered to the later tiers'),
    ('Multimapper', 'a read with two or more equally good best alignments; removed with remove_multimap'),
    ('Antisense read', 'a read that aligns once, on the opposite strand of its feature; removed with remove_antisense'),
    ('Kept reads', 'the mapped reads after the multimapper and antisense settings'),
    ('Natural reads', 'the kept reads without the spike-ins and the vector contaminants (UniVec)'),
    ('Feature', 'a row of the gene table: a gene, a miRNA arm, a tRNA family, ...'),
]


def tier_group(tier):
    for pattern, group in TIER_GROUPS:
        if re.search(pattern, tier, re.IGNORECASE):
            return group
    return 'other tier'


def default_classes():
    """The code's species-neutral class table (runs whose class table is not found)."""
    return refTables.ClassTable(refTables.default_table('default'), 'default')


def type_group(rna_class, classes=None):
    """Report group of an RNA class: the class table's column group (RiboSpike-*: spike-ins); group names of
    tables copied into runs before 0.3.0 S4 are translated (OLD_GROUPS)."""
    classes = classes or default_classes()
    group = classes.group(rna_class)
    if rna_class == unannotated_class(classes):
        return GENOME_UNANNOTATED if group in OLD_GROUPS else group
    return OLD_GROUPS.get(group, group)


def unannotated_class(classes=None):
    """Class id of the genome reads without a feature (genome_unannotated)."""
    return (classes or default_classes()).class_of('unannotated')


def natural(rna_class, classes=None):
    """True for the classes of natural RNA: every kept class except the spike-ins and UniVec (0.3.0 plan §8 item
    11; the collector's assign.natural_mean_length uses the same classes)."""
    return rna_class != (classes or default_classes()).class_of('UniVec') and not rna_class.startswith('RiboSpike')


def class_order(rna_classes, classes=None):
    """Classes in the class table's display order (then by name)."""
    classes = classes or default_classes()
    return sorted(rna_classes, key=lambda c: (classes.display_order(c), c))


def run_classes(out, report):
    """The class table of a run: the copy in <run>/tables/ (0.3.0), else the code's table for its species."""
    species = report['provenance'].get('prov.config', {}).get('species', 'default')
    copy = Path(out) / 'tables' / 'rna_classes.tsv'
    return refTables.ClassTable(str(copy if copy.exists() else refTables.default_table(species)), species)


def colour(group):
    return ENTITY_COLOURS.get(group, NEUTRAL_DARK)


def ordered_groups(groups):
    """GROUP_ORDER first, then other groups (own class tables) in the order given."""
    return [g for g in GROUP_ORDER if g in groups] + [g for g in groups if g not in GROUP_ORDER]


def run_tiers(run):
    """Tier names in tier order (from the collected map.{tier}.reads_in ids)."""
    return [metric_id.split('.')[1] for metric_id in run.ids(r'map\.[^.]+\.reads_in')]


def tier_groups(run, sample, tiers):
    """{group: reads aligned in its tiers} for one sample (before multimapper/antisense removal). The genome
    tier is split (#9): genome (unannotated) = the kept reads of the class genome_unannotated, which is known only
    after assignment; the rest of the genome reads, including the reads removed later, are pre-mRNA and pre-tRNA."""
    row = {}
    for tier in tiers:
        aligned = run.value(sample, f'map.{tier}.aligned')
        if aligned is not None:
            row[tier_group(tier)] = row.get(tier_group(tier), 0) + aligned
    if 'genome' in row:
        genome = row.pop('genome')
        unannotated = min(genome, run.value(sample, f'assign.class.{unannotated_class(run.classes)}.reads') or 0)
        row[GENOME_PRE] = genome - unannotated
        row[GENOME_UNANNOTATED] = unannotated
    return row


def with_percentages(data, categories, title, ylab, tabs=None):
    """A bar plot with MultiQC's percentage switch, shown by report.js as the setting Values: Reads | Percentage
    (RealSeq review #4; Martin 2026-10-26: the same control for every setting). data and categories: one dataset,
    or one per tab (tabs = their names). Percentages are of each sample's bar (MultiQC)."""
    plot = {'plot_type': 'bargraph', 'categories': categories, 'data': data,
            'pconfig': {'title': title, 'ylab': ylab, 'sort_samples': False, 'cpswitch': True,
                        'cpswitch_counts_label': 'Reads', 'cpswitch_percent_label': 'Percentage'}}
    if tabs:
        plot['pconfig']['data_labels'] = [{'name': tab, 'ylab': ylab} for tab in tabs]
    return plot


# ----------------------------------------------------------------------------- data access

class Run:
    """qc_metrics.json of one run, with lookups by sample and metric id."""

    def __init__(self, report, dictionary):
        self.report = report
        self.mode = report['mode']
        self.samples = list(report['samples'])
        self.metrics = {s: {r['id']: r for r in d['metrics']} for s, d in report['samples'].items()}
        self.distributions = {s: {r['id']: r for r in d['distributions']} for s, d in report['samples'].items()}
        self.classes = None             # the run's class table (main: run_classes); None: the code's default
        self.out = None                 # the run folder (main), for the tables shown in the report
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
        ('unaligned', 'map.final_unaligned', True, '% unmapped (of reads to mapping)', '{:,.1f}'),
        ('removed', 'assign.removed_by_policy', True, '% removed by policy (of mapped)', '{:,.1f}'),
        ('natural_length', 'assign.natural_mean_length', False, 'Mean length, natural reads (nt)', '{:,.1f}'),
        ('library', 'counts.library_size', False, 'Library size', '{:,.0f}'),
        ('genes10', 'counts.detected.gt10', False, 'Features > 10 reads', '{:,.0f}'),
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
    # A table section of its own (Martin 2026-10-26): MultiQC always places its General Statistics directly under
    # the header, and the order is Samples, How to read this report, General Statistics (multiqc_config.yaml).
    return {**PARENT, 'id': 'ribomarker_general_stats', 'section_name': 'General Statistics',
            'description': '<p>Headline numbers per sample; hover over a column header for its definition.</p>'
                           + f'<p>{status_note(run)}</p>',
            'plot_type': 'table', 'headers': headers,
            'pconfig': {'id': 'ribomarker_general_stats_table', 'title': 'RiboMarker: general statistics',
                        'col1_header': 'Sample', 'sort_rows': False},
            'data': without_empty(data)}


def read_fate(run):
    pe = run.mode == 'pe'
    unit = 'pairs' if pe else 'reads'
    many = 'Read pairs' if pe else 'Reads'
    categories = [('trim.merged_reads', 'Merged: to mapping')] if pe else [('trim.reads_written', 'Kept: to mapping')]
    categories += [(f'trim.adapter_dimer_{unit}', 'Removed: adapter dimers (empty insert)'),
                   (f'trim.short_insert_{unit}', 'Removed: inserts 1–14 nt')]
    if pe:
        categories += [('trim.fastp_filtered_pairs', 'Removed: fastp filters'),
                       ('trim.unmerged_pairs', 'Removed: not merged')]
    data = {s: {cid: run.value(s, cid) for cid, _ in categories if run.value(s, cid) is not None}
            for s in run.samples}
    if pe:
        description = terms('What happened to each read pair before mapping; the bars add up to the read pairs in.', [
            ('Merged: to mapping', 'fastp merged the two reads of the pair into one read of the whole insert; only '
                                   'these reads are mapped'),
            ('Removed: adapter dimers', 'cutadapt found no insert between the adapters'),
            ('Removed: inserts 1–14 nt', 'the insert is shorter than 15 nt (cutadapt)'),
            ('Removed: fastp filters', 'low quality, too many N or too short (fastp)'),
            ('Removed: not merged', 'the two reads do not overlap enough to be merged')])
    else:
        description = terms('What happened to each read before mapping; the bars add up to the reads in.', [
            ('Kept: to mapping', 'the reads that go on to mapping'),
            ('Removed: adapter dimers', 'cutadapt found no insert between the adapters'),
            ('Removed: inserts 1–14 nt', 'the insert is shorter than 15 nt (cutadapt)')])
    plot = with_percentages(without_empty(data),
                            {cid: {'name': name, 'color': FATE_COLOURS[i]} for i, (cid, name) in enumerate(categories)},
                            'RiboMarker: read fate before mapping', many)
    plot['pconfig']['id'] = 'ribomarker_read_fate_plot'
    return {**PARENT, 'id': 'ribomarker_read_fate', 'section_name': 'Trimming and merging: read fate',
            'description': description, **plot}


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
            row['unmapped'] = run.value(s, 'map.final_unaligned')
        data[s] = row
    groups = [x for g in groups for x in ([GENOME_PRE, GENOME_UNANNOTATED] if g == 'genome' else [g])]
    categories = {g: {'name': g, 'color': NEUTRAL_DARK if g in ('UniVec', 'other tier') else colour(g)} for g in groups}
    categories['unmapped'] = {'name': 'unmapped', 'color': NEUTRAL_LIGHT}
    detail = html_table(['Sample'] + tiers + ['unmapped'],
                        [[html.escape(s)] + [fmt(run.value(s, f'map.{tier}.aligned')) for tier in tiers]
                         + [fmt(run.value(s, 'map.final_unaligned'))] for s in run.samples])
    plot = with_percentages(without_empty(data), categories, 'RiboMarker: reads aligned per tier', 'Reads')
    plot['pconfig']['id'] = 'ribomarker_tiers_plot'
    return {**PARENT, 'id': 'ribomarker_tiers', 'section_name': 'Mapping per reference tier',
            'description': terms(
                'Reads aligned per tier, counted before the multimapper/antisense settings; the tiers are searched in '
                'order, so a read aligned in one tier is not offered to the later tiers. The composition plot below '
                'uses the same groups.', [
                    ('UniVec', 'vector contaminants (Bowtie 1: one alignment per read)'),
                    ('spike-ins', 'the RiboSpike sequences'),
                    ('rRNA', 'the miRTrace rRNA reference (Bowtie 1: one alignment per read)'),
                    ('tRNA', 'mature tRNAs (GtRNAdb)'),
                    ('miRNA', 'the three MirGeneDB tiers together: mature arms, star arms, precursors'),
                    ('piRNA', 'piRBase'),
                    ('lncRNA and other ncRNA', 'the Ensembl ncRNA tier: lncRNA, snRNA, snoRNA, scaRNA, misc_RNA, '
                                               'Mt_rRNA, Mt_tRNA (and 42 rRNA genes: 36 5S, 6 5.8S)'),
                    ('protein_coding', 'the Ensembl cDNA tier: one mature transcript per gene (the longest), i.e. its '
                                       'exons with the 5′ and 3′ UTRs; no introns'),
                    (GENOME_PRE, 'genome reads in the span of a protein-coding gene outside its mature transcript '
                                 '(mostly intronic) or at a tRNA locus, and the genome reads removed later as '
                                 'multimappers or antisense'),
                    (GENOME_UNANNOTATED, 'genome reads without an annotated feature (the kept reads of this class: '
                                         'known only after assignment)'),
                    ('unmapped', 'reads that aligned in no tier')],
                '<details><summary>Reads aligned in every tier</summary>' + detail + '</details>'),
            **plot}


def composition(run):
    biotypes = class_order([metric_id.split('.')[2] for metric_id in run.ids(r'assign\.class\.[^.]+\.reads')],
                           run.classes)
    tiers = run_tiers(run)
    grouped, larger = {}, []
    for s in run.samples:
        row = {}
        for b in biotypes:
            n = run.value(s, f'assign.class.{b}.reads')
            if n is not None:
                group = type_group(b, run.classes)
                row[group] = row.get(group, 0) + n
        before = tier_groups(run, s, tiers)       # a group should not gain reads by the removal
        larger += [f'{s}: {g}' for g, n in row.items() if n > before.get(g, 0)]
        grouped[s] = row
    tabs = [(ALL_CLASSES, set()), ('Without spike-ins', {'spike-ins'}),
            ('Without rRNA and spike-ins', {'rRNA', 'spike-ins'})]       # (tab, groups left out)
    present = ordered_groups(list(dict.fromkeys(g for row in grouped.values() for g in row)))
    categories = {g: {'name': g, 'color': NEUTRAL_DARK if g == 'UniVec' else colour(g)}
                  for g in ordered_groups(GROUP_ORDER + present)}
    detail = html_table(['Sample'] + [run.classes.display_name(b) if run.classes else b for b in biotypes],
                        [[html.escape(s)] + [fmt(run.value(s, f'assign.class.{b}.reads')) for b in biotypes]
                         for s in run.samples])
    notes = ''
    if larger:
        notes += ('<p>Groups with more reads here than aligned in their tiers (an RNA class also comes from another '
                  'tier): ' + html.escape('; '.join(larger)) + '.</p>')
    plot = with_percentages([without_empty({s: {g: n for g, n in row.items() if g not in left_out}
                                            for s, row in grouped.items()}) for _, left_out in tabs],
                            [{g: c for g, c in categories.items() if g not in left_out} for _, left_out in tabs],
                            'RiboMarker: composition by RNA class', 'Reads', tabs=[tab for tab, _ in tabs])
    plot['pconfig']['id'] = 'ribomarker_composition_plot'
    return {**PARENT, 'id': 'ribomarker_composition', 'section_name': 'Composition by RNA class',
            'description': terms(
                'Kept reads per RNA class, after the multimapper/antisense settings, in the groups of the tier plot '
                'above.', [
                    ('Groups', 'each RNA class is shown in the group of the tier its reads align in, so the difference '
                               'between the two plots is the reads removed in that group'),
                    ('Buttons', 'all RNA classes, or without the spike-ins (and rRNA)')],
                notes + '<details><summary>Reads of every RNA class</summary>' + detail + '</details>'),
            **plot}


def length_buttons(run):
    """[(button name, [classes])] of the read-length plots (RealSeq review #11-#16): all natural classes; the
    classes with a button of their own in the class table's display order (pre-miRNA follows miRNA); "Other RNAs"
    (every other natural class, listed in the description); genome (unannotated) last. Spike-ins and UniVec are
    not natural RNA and are left out; classes without reads have no button."""
    classes = run.classes or default_classes()
    present = [m.split('.')[2] for m in run.ids(r'assign\.class\.[^.]+\.reads')]
    present = [c for c in present if natural(c, classes)
               and sum((run.value(s, f'assign.class.{c}.reads') or 0) for s in run.samples) > 0]
    unannotated = unannotated_class(classes)
    own = class_order([c for c in present if c in OWN_LENGTH_BUTTON], classes)
    other = class_order([c for c in present if c not in OWN_LENGTH_BUTTON and c != unannotated], classes)
    buttons = [(ALL_CLASSES, class_order(present, classes))] if present else []
    buttons += [(classes.display_name(c), [c]) for c in own]
    if other:
        buttons.append((OTHER_RNAS, other))
    if unannotated in present:
        buttons.append((classes.display_name(unannotated), [unannotated]))
    return buttons


def length_counts(run, sample, members):
    """{length: reads} of the given classes for one sample, summed."""
    counts = {}
    for c in members:
        dist = run.distributions[sample].get(f'assign.class.{c}.length_hist')
        for x, y in zip(dist['x'], dist['y']) if dist else ():
            counts[x] = counts.get(x, 0) + y
    return counts


LENGTH_MODES = ['% of natural reads', '% within class']


def lengths(run):
    """Read length per RNA class in two modes (#11), one plot: the setting Values (report.js) chooses the mode, the
    buttons below it the class. '% of natural reads' = % of all natural reads of the sample (the paper's 'sequencing
    length profiles', whose denominator the paper does not state; the curves of all buttons but the first add up to
    the first); '% within class' = % of the reads of the button.
    The datasets are mode x button, in that order."""
    buttons = length_buttons(run)
    if not buttons:
        return None
    natural_total = {s: sum(length_counts(run, s, buttons[0][1]).values()) for s in run.samples}
    datasets, labels, longest = [], [], 15
    for mode in ('natural', 'class'):
        for name, members in buttons:
            dataset = {}
            for s in run.samples:
                counts = length_counts(run, s, members)
                n = natural_total[s] if mode == 'natural' else sum(counts.values())
                if n > 0 and sum(counts.values()) > 0:
                    dataset[s] = {x: 100 * y / n for x, y in counts.items() if x >= 15}
                    # the x axis ends where every class's own profile has fallen below 0.1 % (same in both modes)
                    longest = max([longest] + [x for x, y in counts.items() if 100 * y / sum(counts.values()) >= 0.1])
            datasets.append(dataset)
            labels.append({'name': name if mode == 'natural' else f'{name} (within class)',
                           'ylab': '% of natural reads' if mode == 'natural' else '% of the reads of the class',
                           'xlab': 'Read length (nt)'})
    classes = run.classes or default_classes()
    other = next((members for name, members in buttons if name == OTHER_RNAS), [])
    items = [('Values: % of natural reads', 'divides by all natural reads of the sample, so a class\'s curve shows '
                                           'how much of the sample it makes up at each length'),
             ('Values: % within class', 'divides by the reads of the class: its own length profile, whatever its share '
                                        'of the sample'),
             ('protein_coding', 'reads of the mature transcripts (exons with their UTRs)'),
             ('pre-mRNA', 'reads in the span of a protein-coding gene on the genome that do not match its mature '
                          'transcript (mostly intronic, also reads across an exon-intron boundary)')]
    if other:
        items.append((OTHER_RNAS, html.escape(', '.join(classes.display_name(c) for c in other))))
    modes = html.escape(json.dumps(LENGTH_MODES), quote=True)
    return {**PARENT, 'id': 'ribomarker_lengths', 'section_name': 'Read length by RNA class',
            'description': terms('Share of each read length in the kept reads, per RNA class. Lengths are integers, no '
                                 'smoothing; a class with fewer than 100 reads in a sample is noisy here. Spike-in '
                                 'lengths are in the spike-in sections.', items)
                           + f'<div class="rm-modes" data-anchor="ribomarker_lengths_plot" data-label="Values" '
                             f'data-modes="{modes}"></div>',
            'plot_type': 'linegraph',
            'pconfig': {'id': 'ribomarker_lengths_plot', 'title': 'RiboMarker: read length by RNA class',
                        'ylab': labels[0]['ylab'], 'xlab': 'Read length (nt)', 'ymin': 0,
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
    plot = with_percentages(data, categories, 'RiboMarker: multimappers and antisense reads', 'Mapped reads')
    plot['pconfig']['id'] = 'ribomarker_multimapping_plot'
    return {**PARENT, 'id': 'ribomarker_multimapping', 'section_name': 'Multimappers and antisense reads',
            'description': terms('The mapped reads: kept, or removed by the multimapper/antisense settings.', [
                ('removed: multimapper', 'two or more equally good best alignments (remove_multimap)'),
                ('removed: antisense, maps once', 'one alignment, on the opposite strand of its feature '
                                                  '(remove_antisense)'),
                ('UniVec, rRNA, spike-ins', 'these tiers report one alignment per read, so their reads are never '
                                            'multimappers; rRNA is most of the mapped reads')],
                '<details><summary>Reads with a second reported alignment, per tier (% of the reads aligned in the '
                'tier)</summary><p>Per Bowtie 2 tier, the reads with any second reported alignment, also one that '
                'scores worse (Bowtie 2 "aligned &gt;1 times"), as % of the reads aligned in that tier: it shows which '
                'tiers hold multi-copy or similar sequences. The percentages are much higher than the removed share in '
                'the plot mainly because they are per tier, while the plot divides by all mapped reads, most of them '
                'rRNA (never multimappers); also, a second alignment that scores worse counts here but does not make a '
                'read a multimapper.</p>' + table + '</details>'),
            **plot}


def spike_offsets(run):
    lengths_ = sorted({int(m.split('.')[1]) for m in run.ids(r'spike\.\d+\.reads')})
    datasets, labels = [], []
    for length in lengths_:
        dataset = {}
        for s in run.samples:
            dist = run.distributions[s].get(f'spike.{length}.length_offset_hist')
            if dist and dist['n'] > 0:
                dataset[s] = {str(x): 100 * y / dist['n'] for x, y in zip(dist['x'], dist['y'])}
        if dataset:
            datasets.append(dataset)
            labels.append({'name': f'{length} nt', 'ylab': '% of full-length spike reads', 'ymax': 100})
    if not datasets:
        return None
    return {**PARENT, 'id': 'ribomarker_spike_offsets', 'section_name': 'Spike-ins: length fidelity',
            'description': 'Observed minus designed length of the full-length spike-in reads (designed lengths '
                           'include the 4N ends), as % of the full-length reads of that length; one button per '
                           'designed length. Expected at 0; shorter reads (< -2) are counted as fragments. The read '
                           'numbers are in the spike-in summary above.',
            'plot_type': 'bargraph',
            'categories': [{o: {'name': f'{int(o):+d} nt' if o != '0' else '0 nt', 'color': OFFSET_COLOURS[i]}
                            for i, o in enumerate(SPIKE_OFFSETS)}] * len(datasets),
            'pconfig': {'id': 'ribomarker_spike_offsets_plot', 'title': 'RiboMarker: spike-in length offset',
                        'ylab': '% of full-length spike reads', 'sort_samples': False, 'cpswitch': False,
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
    plot = with_percentages(data, categories, 'RiboMarker: spike-in reads by length', 'Full-length spike reads')
    plot['pconfig']['id'] = 'ribomarker_spike_recovery_plot'
    return {**PARENT, 'id': 'ribomarker_spike_recovery', 'section_name': 'Spike-ins: recovery by length',
            'description': 'Full-length spike-in reads per designed length (all four types), as reads or as % of '
                           'the sample\'s full-length spike-in reads. Reads as capture efficiency only if the spikes '
                           'are added in equal amounts.',
            **plot}


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


DETECTION_TABS = [(0, '> 0 reads'), (5, '> 5 reads'), (10, '> 10 reads')]


def counts_table(run):
    """Detection and counts (RealSeq review #17, #18; Martin 2026-10-26: one section, 'features', the same
    thresholds in the plot and the table): detected features per RNA class and sample (plot), the per-sample table
    and the top lists below it, collapsed."""
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
    species = run.report.get('provenance', {}).get('prov.config', {}).get('species', 'human')
    mir451, mir23a = RBC_MIRNA_NAMES.get(species, RBC_MIRNA_NAMES['human'])
    headers = ['Sample', 'Library size', 'Without rRNA'] + [f'Features {tab}' for _, tab in DETECTION_TABS] + [
               '% in top-10 features', '% in top-10 features without rRNA',
               f'miR-451a / miR-23a-3p (MirGeneDB {mir451} / {mir23a}; reads)']
    classes = run.classes or default_classes()
    biotypes = class_order({m.split('.')[1] for m in run.ids(r'counts\.[^.]+\.top20')}, classes)
    per_type = ''.join(f'<details><summary>{html.escape(classes.display_name(b))}</summary>'
                       f'{top_table(run, f"counts.{b}.top20")}</details>' for b in biotypes)
    tables = ('<div class="rm-after-plot"><details><summary>Counts per sample (table)</summary>' + html_table(headers, rows) + '</details>'
              + '<details><summary>20 most abundant features per sample (reads)</summary>'
              + top_table(run, 'counts.top20') + '</details>'
              + '<details><summary>20 most abundant features per sample, by RNA class (reads)</summary>'
              + per_type + '</details></div>')                 # report.js moves it below the plot
    section = {**PARENT, 'id': 'ribomarker_counts', 'section_name': 'Detection and counts',
               'description': terms(
                   'Counts of the gene table (*-genecounts.csv: mature exonic reads after the multimapper/antisense '
                   'settings).', [
                       ('Feature', 'a row of the gene table: a gene, a miRNA arm, a tRNA family, ...'),
                       ('Plot', 'the features with more than 0, 5 or 10 reads per RNA class and sample: one group of '
                                'bars per class, one bar per sample'
                                + (', coloured by ' + html.escape(str(report_group(run))) if report_group(run) else '')
                                + '. Classes that are not in the gene table (pre-mRNA, pre-tRNA, pre-miRNA, '
                                'unannotated, spike-ins, UniVec) have no bars; Scale: Log makes the small classes '
                                'readable'),
                       ('Below the plot', 'collapsed: the counts per sample, and the most abundant features (as many '
                                          'ranks as there are features with reads)'),
                       ('miRNA ratio', 'miR-451a / miR-23a-3p, a marker of red-blood-cell contamination: a plain ratio '
                                       'of reads, descriptive only')]) + tables}
    found = class_order({m.split('.')[2] for m in run.ids(r'counts\.detected\.[^.]+\.gt\d+')}, classes)
    if not found:                     # no gene table (remove_multimap / remove_antisense off): the tables only
        return {**section, 'plot_type': 'html', 'data': ''}
    datasets, labels = [], []
    for above, tab in DETECTION_TABS:
        dataset = {}
        for c in found:
            row = {s: run.value(s, f'counts.detected.{c}.gt{above}') for s in run.samples}
            row = {s: n for s, n in row.items() if n is not None}
            if any(row.values()):
                dataset[classes.display_name(c)] = row
        datasets.append(dataset)
        labels.append({'name': tab, 'ylab': 'Features detected'})
    colours = sample_colours(run) or (dict(zip(run.samples, SAMPLE_GROUP_COLOURS))
                                      if len(run.samples) <= len(SAMPLE_GROUP_COLOURS) else {})
    categories = {s: {'name': s, **({'color': colours[s]} if s in colours else {})} for s in run.samples}
    return {**section, 'plot_type': 'bargraph', 'categories': [categories] * len(datasets),
            'pconfig': {'id': 'ribomarker_counts_plot', 'title': 'RiboMarker: detected features per RNA class',
                        'ylab': 'Features detected', 'sort_samples': False, 'cpswitch': False, 'stacking': 'group',
                        'logswitch': True, 'data_labels': labels},
            'data': datasets}


def top_table(run, metric_id):
    """Rank x sample table of a top-20 list metric ('gene (reads)'); only the ranks that some sample fills (#17)."""
    rows = []
    for rank in range(max([0] + [len(run.value(s, metric_id) or []) for s in run.samples])):
        row = [str(rank + 1)]
        for s in run.samples:
            top = run.value(s, metric_id) or []
            row.append(f'{html.escape(top[rank][0])} ({top[rank][1]:,})' if rank < len(top) else '')
        rows.append(row)
    return html_table(['Rank'] + run.samples, rows)


def check_summary(run):
    """'44 equal, 1 differs, ...': the consistency checks by status."""
    statuses = [c['status'] for c in run.report.get('checks', [])]
    return ', '.join(f'{statuses.count(s)} {s}' for s in ('equal', 'differs', 'not evaluable', 'not applicable')
                     if statuses.count(s))


def status_line(run):
    """One line at the top of the Samples section (Martin 2026-10-26): the checks and warnings at a glance; the
    details are under General Statistics."""
    report = run.report
    problems = sum(c['status'] in ('differs', 'not evaluable') for c in report.get('checks', []))
    warnings = len(report['warnings'])
    details = ' (details under General Statistics)' if problems or warnings else ''
    return (f'<b>Consistency checks:</b> {html.escape(check_summary(run))}{details if problems else ""}. '
            f'<b>Warnings:</b> {warnings or "none"}{details if warnings else ""}.')


def status_note(run):
    """Consistency checks and warnings in the General Statistics section (Martin 2026-10-26: instead of a Run
    overview section): the summary line, then the checks that are not equal and the warnings, if any."""
    report = run.report
    problems = [c for c in report.get('checks', []) if c['status'] in ('differs', 'not evaluable')]
    warnings = report['warnings']
    parts = [f'<b>Consistency checks:</b> {html.escape(check_summary(run))} (every check: '
             '<code>qc/consistency_checks.tsv</code>). '
             + ('<b>Warnings:</b> none.' if not warnings else f'<b>Warnings:</b> {len(warnings)}, listed below.')]
    if problems:   # the only place checks are listed in the report; all of them are in qc/consistency_checks.tsv
        parts.append('<p><b>Checks that are not equal:</b></p>' + checks_table(report, problems))
    if warnings:
        parts.append('<p><b>Warnings</b> (missing or unreadable inputs; the affected metrics are omitted):</p>'
                     + html_table(['Sample', 'Source', 'Message'],
                                  [[html.escape(w['sample'] or 'run'), html.escape(w['source']),
                                    html.escape(w['message'])] for w in warnings]))
    return ''.join(parts)


def checks_table(report, checks):
    definitions = report.get('check_definitions', {})
    rows = [[html.escape(c['id']), html.escape(c['status']), html.escape(c['sample'] or 'run'),
             html.escape(c['subject'] or ''), html.escape(str(c['expected'])), html.escape(str(c['observed'])),
             html.escape(definitions.get(c['id'], ''))] for c in checks]
    return html_table(['Check', 'Status', 'Sample', 'Subject', 'Expected', 'Observed', 'Identity'], rows)


def provenance(run):
    prov = run.report['provenance']
    report = run.report
    rows = [['Run', html.escape(report['run_id'])],
            ['Mode', 'paired-end' if run.mode == 'pe' else 'single-end'],
            ['Pipeline version', html.escape(pipeline_label(prov.get('prov.pipeline_version', {})))],
            ['Samples', str(len(run.samples))]]
    for key in ('prov.library_parameters', 'prov.reference_species_check'):   # tool versions: software_versions()
        value = prov.get(key)
        if isinstance(value, dict):
            value = '<br>'.join(f'{html.escape(str(k))}: {html.escape(str(v))}' for k, v in value.items())
        rows.append([html.escape(run.entry(key)['label']), value if value is not None else '—'])
    tables = prov.get('prov.class_tables')
    if tables:
        rows.append([html.escape(run.entry('prov.class_tables')['label']),
                     '<br>'.join(f"{html.escape(k)}: {html.escape(str(v.get('file')))}"
                                 f"{' (replaces the default)' if v.get('replaced') else ''}"
                                 f"{' (not applied: built for another references BED)' if v.get('applied') is False else ''}"
                                 f", sha256 {html.escape(str(v.get('sha256', ''))[:12])}" for k, v in tables.items())])
    build = prov.get('prov.reference_build')
    if build:
        sources = build.get('sources', {})
        text = (f"{html.escape(str(build.get('species')))}, build {html.escape(str(build.get('build')))}, "
                f"built {html.escape(str(build.get('built_utc')))}<br>"
                + '<br>'.join(f'{html.escape(k)}: {html.escape(str(v.get("version")))}' for k, v in sources.items()))
        rows.append(['Reference build', text])
    else:
        rows.append(['Reference build', '— (no MANIFEST found)'])
    for i, citation in enumerate(CITATIONS):
        rows.append(['Please cite' if i == 0 else '', html.escape(citation.format(
            version=prov.get('prov.pipeline_version', {}).get('version', '—')))])
    config = prov.get('prov.config', {})
    rows.append(['Configuration', '<br>'.join(f'{html.escape(str(k))}: {html.escape(str(v))}'
                                               for k, v in config.items())])
    # a section of its own (no parent), placed at the end just above MultiQC's Software Versions
    return {'id': 'ribomarker_provenance', 'section_name': 'Provenance',
            'description': 'The run, its parameters and references, and how to cite the pipeline; the tool versions '
                           'are listed under Software Versions below. Absolute paths are reduced to their last '
                           'component.', 'plot_type': 'html',
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



def guide(run):
    """How to read the report (RealSeq review #3)."""
    items = [
        ('Settings', 'blue, rounded, on their own line above a plot: they change how the plot shows its numbers. '
                     '<i>Values</i>: Reads or Percentage (% of each sample\'s bar), or % of natural reads or % within '
                     'class in the read-length plot; <i>Scale</i>: Linear or Log. The grey buttons below them choose '
                     'what is shown, e.g. one RNA class or a threshold'),
        ('Toolbox', 'the tab at the right edge: highlight, rename or hide samples in every plot, and export plots (PNG, '
                    'SVG, PDF) or their data. A plot can also be exported with the <i>Export</i> button above it'),
        ('Tables', '<i>Configure columns</i> above the General Statistics table and the sample table chooses the '
                   'columns; a click on a column header sorts by it'),
        ('Line plots', 'hover over a line to see its sample; click a sample in the legend to hide or show it'),
    ]
    if colour_selector(run):
        items.append(('Colours', '<i>Colour samples by</i>, above the sample table, chooses the sample-sheet column '
                                 'whose values colour the samples in every plot and table: their names in the bar plots '
                                 'and tables, their lines in the line plots (at most 8 groups; the configuration key '
                                 'report_group sets the first choice). A click on a colour in the key changes it. The toolbox\'s '
                                 'Highlight tab shows the colours'))
    return {**PARENT, 'id': 'ribomarker_guide', 'section_name': 'How to read this report',
            'description': 'Navigation and terms used in this report.', 'plot_type': 'html',
            'data': '<h4>Using the report</h4>' + terms('', items).replace('<p></p>', '')
                    + '<h4>Glossary</h4>' + terms('', GLOSSARY).replace('<p></p>', '')
                    + '<script>' + (MULTIQC_CONFIG.parent / REPORT_JS).read_text(encoding='utf-8') + '</script>'}


def samples_table(run):
    """The user's sample sheet (RealSeq review #1); without one, the samples and their FASTQ files."""
    sheet = run.report['provenance'].get('prov.sample_sheet')
    if sheet:
        columns = [c for c in sheet['columns'] if c != 'sample']
        sample_col = sheet['columns'].index('sample')
        data = {row[sample_col]: {c: row[sheet['columns'].index(c)] for c in columns} for row in sheet['rows']}
        source = f"the sample sheet ({html.escape(sheet['file'])})"
    else:
        columns = ['fastq_file', 'fastq_file2'] if run.mode == 'pe' else ['fastq_file']
        rows = run.report['provenance'].get('prov.samples') or []
        data = {row[0]: dict(zip(columns, row[2:])) for row in rows}
        source = 'the run\'s sample list (no sample sheet: on the configuration route give it with sample_sheet:)'
    if not data:
        return None
    headers = {c: {'title': c, 'scale': False, 'hidden': c in ('fastq_file', 'fastq_file2')} for c in columns}
    return {**PARENT, 'id': 'ribomarker_samples', 'section_name': 'Samples',
            'description': f'<p>{status_line(run)}</p>'
                           f'The samples of the run, from {source}, in its order. <i>Configure columns</i> chooses '
                           'the columns shown (the FASTQ file names are hidden at first).' + colour_selector(run),
            'plot_type': 'table', 'headers': headers,
            'pconfig': {'id': 'ribomarker_samples_table', 'title': 'RiboMarker: samples', 'col1_header': 'Sample',
                        'sort_rows': False, 'no_violin': True},
            'data': {s: row for s, row in data.items()}}


def colour_selector(run):
    """The 'Colour samples by' selector above the sample table (Martin 2026-10-26, S4c): report.js fills it and
    colours the samples by the chosen sample-sheet column in every plot and table through MultiQC's highlighting.
    Its data (samples, the values of every sheet column except the FASTQ files, the palette, the first choice =
    report_group) travel base64-encoded, so that MultiQC's Markdown reading of descriptions cannot change them.
    Empty without a sample sheet."""
    sheet = run.report['provenance'].get('prov.sample_sheet') or {}
    columns = [c for c in sheet.get('columns', []) if c not in ('sample', 'fastq_file', 'fastq_file2')]
    if not columns:
        return ''
    by_sample = {row[sheet['columns'].index('sample')]: row for row in sheet['rows']}
    samples = [s for s in run.samples if s in by_sample]
    values = {c: [by_sample[s][sheet['columns'].index(c)] or '(empty)' for s in samples] for c in columns}
    data = {'samples': samples, 'columns': values, 'palette': SAMPLE_GROUP_COLOURS,
            'first': report_group(run) or '', 'series_plots': ['ribomarker_counts_plot']}
    encoded = base64.b64encode(json.dumps(data).encode('utf-8')).decode('ascii')
    return f'<div class="rm-colour-by" data-config="{encoded}"></div>'


def report_group(run):
    """The sample-sheet column that groups the samples (configuration key report_group), or None."""
    return (run.report['provenance'].get('prov.sample_sheet') or {}).get('report_group')


def sample_group_of(run):
    """{sample: group} by report_group, in sample order; {} without one."""
    sheet = run.report['provenance'].get('prov.sample_sheet') or {}
    group = sheet.get('report_group')
    if not group:
        return {}
    columns = sheet['columns']
    return {row[columns.index('sample')]: row[columns.index(group)] or '(empty)' for row in sheet['rows']}


def sample_colours(run):
    """{sample: colour} by group; {} without groups or with more groups than colours."""
    groups = sample_group_of(run)
    names = list(dict.fromkeys(groups.values()))
    if not names or len(names) > len(SAMPLE_GROUP_COLOURS):
        return {}
    return {s: SAMPLE_GROUP_COLOURS[names.index(g)] for s, g in groups.items()}


def tables_section(run):
    """The tier list and the RNA class table of the run (0.3.0 plan, item 4: shown in the report)."""
    out = Path(run.out) if getattr(run, 'out', None) else None
    parts = []
    tiers_file = out / 'tables' / 'tiers.csv' if out else None
    if tiers_file and tiers_file.exists():
        rows = [line.split(',') for line in tiers_file.read_text().splitlines() if line.strip()]
        parts.append('<details><summary>Tier list (mapping order)</summary>'
                     + html_table(['#', 'Aligner', 'Index'],
                                  [[str(i + 1), html.escape(r[0]), html.escape(Path(r[1]).name if len(r) > 1 else '')]
                                   for i, r in enumerate(rows)])
                     + '</details>')
    classes = run.classes or default_classes()
    by_class = {}
    for row in classes.rows:
        by_class.setdefault(row['class_id'], []).append(row)
    order = class_order(by_class, classes)
    parts.append('<details><summary>RNA classes (class table)</summary>'
                 + html_table(['RNA class', 'Shown as', 'Group', 'From labels (tier)', 'In the gene table', 'Description'],
                              [[html.escape(c), html.escape(classes.display_name(c)), html.escape(type_group(c, classes)),
                                html.escape(', '.join(f"{r['source_label']} ({r['tier']})" for r in by_class[c])),
                                html.escape(by_class[c][0]['gene_table']), html.escape(classes.description(c))]
                               for c in order])
                 + '</details>')
    return {**PARENT, 'id': 'ribomarker_tables', 'section_name': 'Tiers and RNA classes',
            'description': 'The tier list and the RNA class table this run used (their files and checksums are under '
                           'Provenance). The class table gives each BED label its RNA class, the class\'s name in '
                           'the report, its group in the tier and composition plots, its order, and whether it is '
                           'in the gene table.',
            'plot_type': 'html', 'data': ''.join(parts)}


MULTIQC_VERSION = '1.35'   # the version report.js was written for (docker/env.lock.txt pins it)

SECTIONS = [samples_table, guide, general_stats, lane_table, read_fate, tier_fate, composition, multimapping, lengths,
            spike_table, spike_offsets, spike_recovery, counts_table, tables_section, provenance]


# ----------------------------------------------------------------------------- main

def write_custom_content(run, mqc_dir):
    mqc_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for section in SECTIONS:
        content = section(run)
        if content is None:
            continue
        if 'description' in content:   # MultiQC reads descriptions as Markdown: a '*' (as in *-genecounts.csv) would
            content['description'] = content['description'].replace('*', '&#42;')   # open an <em> that never closes
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
    run.classes = run_classes(out, report)
    run.out = out

    custom = write_custom_content(run, qc_dir / 'mqc')
    inputs = custom + tool_reports(out, run.samples, run.mode)
    # the multiqc command of the environment; python -m multiqc if it is not on PATH
    multiqc = [shutil.which('multiqc')] if shutil.which('multiqc') else [sys.executable, '-m', 'multiqc']
    command = multiqc + ['--force', '--quiet', '--no-data-dir',
               '--config', str(MULTIQC_CONFIG), '--title', args.title or report['run_id'],
               '--cl-config', 'custom_logo: ' + json.dumps(str((MULTIQC_CONFIG.parent / 'logos' / 'report_logos.png').resolve())),
               '--cl-config', 'custom_css_files: ' + json.dumps([str((MULTIQC_CONFIG.parent / 'report.css').resolve())]),
               '--outdir', str(qc_dir), '--filename', 'report_technical.html'] + [str(p) for p in inputs]
    # report.js uses MultiQC 1.35's report scripts (highlighting, plot data, colour format); another version may
    # break the settings line and the colour selector (critique 2026-10-27): say so
    version = subprocess.run(multiqc + ['--version'], capture_output=True, text=True).stdout.strip()
    if version and MULTIQC_VERSION not in version:
        print(f'qcReport: WARNING: {version}; the report scripts (docker/qc/report.js) were written for MultiQC '
              f'{MULTIQC_VERSION}: check the settings above the plots and "Colour samples by"', file=sys.stderr)
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f'MultiQC failed:\n{result.stderr}')
    size = (qc_dir / 'report_technical.html').stat().st_size
    print(f'qcReport: {qc_dir / "report_technical.html"} ({size / 1e6:.1f} MB; {sum(p.name.endswith("_mqc.yaml") for p in custom)} RiboMarker sections, '
          f'{len(inputs) - len(custom)} tool reports)')


if __name__ == '__main__':
    main()
