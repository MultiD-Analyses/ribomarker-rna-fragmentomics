"""The technical report: MultiQC renders it from the collector's output; the
numbers in the RiboMarker sections are the collected values; no absolute path reaches the HTML."""
import base64
import csv
import html
import json
import re
import subprocess
import sys

import pytest
import yaml

from conftest import REPO, prepare_run, run_collector

REPORTER = REPO / 'docker' / 'scripts' / 'qcReport.py'


@pytest.fixture(scope='session')
def reports(tmp_path_factory):
    """{mode: (run folder, qc_metrics.json)} after the collector and the reporter ran on the fixtures."""
    runs = {}
    for mode in ('pe', 'se'):
        run, config_file = prepare_run(tmp_path_factory.mktemp(f'report_{mode}'), mode)
        report = run_collector(run, config_file)
        result = subprocess.run([sys.executable, str(REPORTER), '--output', str(run)], cwd=REPO,
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        runs[mode] = (run, report)
    return runs


@pytest.fixture(params=['pe', 'se'])
def reported(request, reports):
    run, report = reports[request.param]
    long = {}
    with open(run / 'qc' / 'metrics_long.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            long[(row['sample'], row['metric_id'])] = float(row['value'])
    sections = {path.name[:-len('_mqc.yaml')]: yaml.safe_load(path.read_text())
                for path in (run / 'qc' / 'mqc').glob('*_mqc.yaml')}
    records = {s: {r['id']: r for r in d['metrics']} for s, d in report['samples'].items()}
    return request.param, run, report, long, sections, records


def test_html_exists_without_absolute_paths(reported):
    mode, run, _, _, _, _ = reported
    text = (run / 'qc' / 'report_technical.html').read_text()
    assert str(run.parent) not in text and str(REPO) not in text
    assert not re.search(r'/(home|mnt|Users|sessions)/', text)
    assert ('id="fastp"' in text) == (mode == 'pe')          # fastp module only in paired-end mode


def test_read_fate_equals_collected_counts(reported):
    mode, _, _, long, sections, _ = reported
    section = sections['ribomarker_read_fate']
    data = section['data']
    reads_in = 'trim.pairs_in' if mode == 'pe' else 'trim.reads_in'
    # #4: MultiQC's percentage switch, shown by report.js as "Values: Reads | Percentage"
    assert section['pconfig']['cpswitch'] is True and section['pconfig']['cpswitch_counts_label'] == 'Reads'
    assert list(data) == ['S1', 'S2', 'S3']                 # sample-sheet order
    for sample, row in data.items():
        for metric_id, value in row.items():
            assert value == long[(sample, metric_id)]
        assert sum(row.values()) == long[(sample, reads_in)]   # categories add up to the input
    names = [c['name'] for c in section['categories'].values()]
    assert names[0] == ('Merged: to mapping' if mode == 'pe' else 'Kept: to mapping')   # #5, #6
    assert all(name.startswith('Removed: ') for name in names[1:])
    assert 'only these reads are mapped' in section['description'] or 'go on to mapping' in section['description']  # #6, #7


def test_tiers_and_composition_add_up(reported):
    _, _, _, long, sections, records = reported
    for sample, row in sections['ribomarker_tiers']['data'].items():
        assert sum(row.values()) == long[(sample, 'map.total_in')]
    all_types, no_spikes, no_rrna_spikes = sections['ribomarker_composition']['data']
    for sample in all_types:
        assert sum(all_types[sample].values()) == long[(sample, 'assign.kept')]
        assert all_types[sample].get('UniVec') == long[(sample, 'map.univec_bt1.aligned')]
        # the rRNA group is the rRNA tier; Mt_rRNA is in the ncRNA group for now
        rrna = sum(r['value'] for i, r in records[sample].items() if i == 'assign.class.rRNA.reads')
        spikes = sum(r['value'] for i, r in records[sample].items()
                     if i.startswith('assign.class.RiboSpike') and i.endswith('.reads'))
        assert sum(no_spikes.get(sample, {}).values()) == long[(sample, 'assign.kept')] - spikes
        assert sum(no_rrna_spikes.get(sample, {}).values()) == long[(sample, 'assign.kept')] - rrna - spikes


# What each general-statistics column must show, written out by hand so that a title and the metric
# behind it are reviewed together (the earlier test accepted any metric value in any column, which let a
# column titled "Reads mapped" show the reads entering mapping; review 2026-10-02).
# column: (title, metric id, percent?) per mode
GENERAL_STATS_SPEC = {
    'pe': {'reads_in': ('Pairs in', 'trim.pairs_in', False),
           'short': ('% < 15 nt (of pairs in)', 'trim.pairs_too_short', True),
           'merged': ('% merged (of pairs after cutadapt)', 'trim.merged_reads', True),
           'mapping_in': ('Reads to mapping', 'map.total_in', False),
           'unaligned': ('% unmapped (of reads to mapping)', 'map.final_unaligned', True),
           'removed': ('% removed by policy (of mapped)', 'assign.removed_by_policy', True),
           'natural_length': ('Mean length, natural reads (nt)', 'assign.natural_mean_length', False),
           'library': ('Library size', 'counts.library_size', False),
           'genes10': ('Features > 10 reads', 'counts.detected.gt10', False)},
    'se': {'reads_in': ('Reads in', 'trim.reads_in', False),
           'short': ('% < 15 nt (of reads in)', 'trim.reads_too_short', True),
           'mapping_in': ('Reads to mapping', 'map.total_in', False),
           'unaligned': ('% unmapped (of reads to mapping)', 'map.final_unaligned', True),
           'removed': ('% removed by policy (of mapped)', 'assign.removed_by_policy', True),
           'natural_length': ('Mean length, natural reads (nt)', 'assign.natural_mean_length', False),
           'library': ('Library size', 'counts.library_size', False),
           'genes10': ('Features > 10 reads', 'counts.detected.gt10', False)},
}


def test_general_stats_columns_show_their_metric(reported):
    mode, _, _, _, sections, records = reported
    section = sections['ribomarker_general_stats']
    spec = GENERAL_STATS_SPEC[mode]
    assert list(section['headers']) == list(spec)
    for column, (title, metric_id, is_percent) in spec.items():
        header = section['headers'][column]
        assert header['title'] == title
        assert ('%' in title) == is_percent and (header['suffix'] == '%') == is_percent
        denominator_id = None
        for sample, row in section['data'].items():
            record = records[sample][metric_id]
            denominator_id = record.get('denominator_id')
            expected = 100 * record['share'] if is_percent else record['value']
            assert row[column] == pytest.approx(expected), f'{sample} {column}'
        if is_percent:   # the help text names the denominator the collector used
            assert f'({denominator_id}).' in header['description'], column
        if re.search(r'\bmapped\b', title.lower()):   # "mapped" (not "unmapped") in a title must mean assign.mapped
            assert 'assign.mapped' in (metric_id, denominator_id), column


def table_cells(html_text):
    """{row label: {column header: cell text without tags}} of the first HTML table."""
    first = html_text.split('</table>')[0]
    headers = [html.unescape(re.sub('<[^>]+>', '', h)) for h in re.findall(r'<th>(.*?)</th>', first)]
    rows = {}
    for row in re.findall(r'<tr>(.*?)</tr>', first.split('<tbody>')[1]):
        cells = [html.unescape(re.sub('<[^>]+>', ' ', c)).strip() for c in re.findall(r'<td>(.*?)</td>', row)]
        rows[cells[0]] = dict(zip(headers, cells))
    return rows


def number(cell):
    return float(cell.split()[0].replace(',', '').rstrip('%'))


def test_spike_and_counts_tables_show_their_metric(reported):
    _, _, _, _, sections, records = reported
    spike = table_cells(sections['ribomarker_spike_table']['data'])
    counts = table_cells(sections['ribomarker_counts']['description'])
    for sample, metrics in records.items():
        if metrics.get('spike.reads', {}).get('value'):
            row = spike[sample]
            assert number(row['Spike-in reads']) == metrics['spike.reads']['value']
            assert number(row['% of reads into mapping']) == round(100 * metrics['spike.reads']['share'], 2)
            assert number(row['% fragments (of spike-in reads)']) == round(100 * metrics['spike.fragment_reads']['share'], 1)
        row = counts[sample]
        for header, metric_id in [('Library size', 'counts.library_size'), ('Without rRNA', 'counts.library_size_non_rrna'),
                                  ('Features > 0 reads', 'counts.detected.gt0'), ('Features > 5 reads', 'counts.detected.gt5'),
                                  ('Features > 10 reads', 'counts.detected.gt10')]:
            assert number(row[header]) == metrics[metric_id]['value'], f'{sample} {header}'
        for header, metric_id in [('% in top-10 features', 'counts.top10_reads'),
                                  ('% in top-10 features without rRNA', 'counts.top10_reads_non_rrna')]:
            assert number(row[header]) == round(100 * metrics[metric_id]['share'], 1), f'{sample} {header}'


def test_spike_sections_equal_collected_values(reported):
    _, _, report, long, sections, _ = reported
    recovery = sections['ribomarker_spike_recovery']['data']
    for sample, row in recovery.items():
        for name, value in row.items():
            assert value == long[(sample, f"spike.{name.split()[0]}.reads")]
    for dataset in sections['ribomarker_spike_offsets']['data']:     # % of the full-length reads of a length
        for sample, row in dataset.items():
            assert sum(row.values()) == pytest.approx(100)
    assert 'S2' not in recovery and 'S3' not in recovery     # S2 and S3 have no spike-in reads


def length_hist(report, sample, rna_class):
    return next((d for d in report['samples'][sample]['distributions']
                 if d['id'] == f'assign.class.{rna_class}.length_hist'), None)


def button_classes(sections, report):
    """{button name: [classes]} of the read-length plots, from the class table (display names)."""
    import qcReport
    import refTables
    classes = refTables.ClassTable(refTables.default_table('human'), 'human')
    present = [d['id'].split('.')[2] for d in report['samples']['S1']['distributions']
               if d['id'].startswith('assign.class.') and d['id'].endswith('.length_hist')]
    natural = [c for c in present if qcReport.natural(c, classes) and c != 'unmapped'
               and any((length_hist(report, s, c) or {}).get('n') for s in report['samples'])]
    names = {classes.display_name(c): [c] for c in natural}
    names['All RNA classes'] = natural
    names['Other RNAs'] = [c for c in natural if c not in qcReport.OWN_LENGTH_BUTTON and c != 'genome_unannotated']
    return names


def length_labels(sections):
    """(button names, number of buttons) of the read-length plot; its datasets are mode x button."""
    labels = [label['name'] for label in sections['ribomarker_lengths']['pconfig']['data_labels']]
    n = len(labels) // 2
    assert labels[n:] == [f'{name} (within class)' for name in labels[:n]]
    return labels[:n], n


def test_length_buttons_follow_the_class_table(reported):
    """#12-#16: all classes first, then the classes in the class table's order (pre-miRNA right after miRNA),
    "Other RNAs", genome (unannotated) last; no spike-in or UniVec button."""
    import refTables
    _, _, report, _, sections, _ = reported
    classes = refTables.ClassTable(refTables.default_table('human'), 'human')
    labels, _ = length_labels(sections)
    assert labels[0] == 'All RNA classes' and labels[-1] == 'genome (unannotated)'
    assert not any(name.startswith('RiboSpike') or name.startswith('vector') for name in labels)
    middle = [n for n in labels[1:-1] if n != 'Other RNAs']
    orders = [classes.display_order(next(c for c in classes._by_class if classes.display_name(c) == n)) for n in middle]
    assert orders == sorted(orders)
    if 'pre-miRNA' in labels:
        assert labels.index('pre-miRNA') == labels.index('miRNA') + 1
    if 'Other RNAs' in labels:
        assert labels.index('Other RNAs') == len(labels) - 2
        assert 'Other RNAs' in sections['ribomarker_lengths']['description']


def test_length_curves_in_both_modes(reported):
    """#11: one plot; the setting Values (report.js, div.rm-modes) chooses % of natural reads (the sample's natural
    reads: every class but the spike-ins and UniVec) or % within the button's classes."""
    _, _, report, _, sections, _ = reported
    members = button_classes(sections, report)
    labels, n = length_labels(sections)
    section = sections['ribomarker_lengths']
    assert 'class="rm-modes" data-anchor="ribomarker_lengths_plot"' in section['description']
    for mode in ('natural', 'class'):
        datasets = section['data'][:n] if mode == 'natural' else section['data'][n:]
        for name, dataset in zip(labels, datasets):
            for sample, curve in dataset.items():
                dists = [length_hist(report, sample, c) for c in members[name]]
                dists = [d for d in dists if d]
                natural_n = sum(length_hist(report, sample, c)['n'] for c in members['All RNA classes']
                                if length_hist(report, sample, c))
                total = natural_n if mode == 'natural' else sum(d['n'] for d in dists)
                for x, share in curve.items():
                    count = sum(d['y'][d['x'].index(int(x))] for d in dists)
                    assert share == pytest.approx(100 * count / total), (mode, name, sample, x)
    # natural mode: the curves of the buttons after the first add up to the first
    first, rest = section['data'][0], section['data'][1:n]
    for sample, curve in first.items():
        for x, share in curve.items():
            assert sum(d.get(sample, {}).get(x, 0) for d in rest) == pytest.approx(share)


def test_natural_mean_length_from_the_class_lengths(reported):
    _, run, report, _, _, records = reported
    rows = list(csv.reader(open(next((run / 'results').glob('*-classlengths.csv')))))
    lengths = [int(x) for x in rows[0][2:]]
    for sample in report['samples']:
        total = reads = 0
        for row in rows[1:]:
            if row[0] == sample and row[1] not in ('unmapped', 'UniVec') and not row[1].startswith('RiboSpike'):
                counts = [int(float(n)) for n in row[2:]]
                total += sum(x * n for x, n in zip(lengths, counts))
                reads += sum(counts)
        assert records[sample]['assign.natural_mean_length']['value'] == pytest.approx(total / reads)
        assert records[sample]['assign.natural_mean_length']['n'] == reads


def test_multimapping_bars_add_up_to_mapped(reported):
    _, _, report, _, sections, _ = reported
    for sample, row in sections['ribomarker_multimapping']['data'].items():
        values = {r['id']: r['value'] for r in report['samples'][sample]['metrics']}
        assert sum(row.values()) == values['assign.mapped']
        assert row['removed: multimapper'] + row['removed: antisense, maps once'] == values['assign.removed_by_policy']


def test_table_columns_with_distinct_values():
    """The same column checks on made-up numbers that are all different, so that two swapped columns
    cannot pass by holding equal values (in the synthetic fixtures several counts coincide)."""
    import qcReport
    dictionary = yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8'))
    ids = ['trim.pairs_in', 'trim.pairs_too_short', 'trim.merged_reads', 'map.total_in', 'map.final_unaligned',
           'assign.mapped', 'assign.removed_by_policy', 'assign.natural_mean_length', 'counts.library_size',
           'counts.library_size_non_rrna',
           'counts.detected.gt0', 'counts.detected.gt5', 'counts.detected.gt10', 'counts.top10_reads',
           'counts.top10_reads_non_rrna', 'spike.reads', 'spike.fragment_reads', 'spike.20.reads']
    records = [{'id': metric_id, 'value': 1000 + 7 * i, 'share': (i + 1) / 97} for i, metric_id in enumerate(ids)]
    for record in records:   # denominators as the collector writes them
        entry = next(e for e in dictionary['metrics'] if e['id'] == record['id'] or
                     re.fullmatch(re.sub(r'\\\{\w+\\\}', r'[^.]+', re.escape(e['id'])), record['id']))
        if entry.get('denominator_id'):
            record['denominator_id'] = entry['denominator_id'].replace('{len}', '20')
    report = {'run_id': 'x', 'mode': 'pe', 'samples': {'X': {'metrics': records, 'distributions': []}},
              'provenance': {}, 'checks': [], 'warnings': []}
    run = qcReport.Run(report, dictionary)
    by_id = {r['id']: r for r in records}

    stats = qcReport.general_stats(run)
    for column, (title, metric_id, is_percent) in GENERAL_STATS_SPEC['pe'].items():
        record = by_id[metric_id]
        expected = 100 * record['share'] if is_percent else record['value']
        assert stats['data']['X'][column] == pytest.approx(expected), column

    counts = table_cells(qcReport.counts_table(run)['description'])['X']
    for header, metric_id in [('Library size', 'counts.library_size'), ('Without rRNA', 'counts.library_size_non_rrna'),
                              ('Features > 0 reads', 'counts.detected.gt0'), ('Features > 5 reads', 'counts.detected.gt5'),
                              ('Features > 10 reads', 'counts.detected.gt10')]:
        assert number(counts[header]) == by_id[metric_id]['value'], header
    for header, metric_id in [('% in top-10 features', 'counts.top10_reads'),
                              ('% in top-10 features without rRNA', 'counts.top10_reads_non_rrna')]:
        assert number(counts[header]) == round(100 * by_id[metric_id]['share'], 1), header

    spike = table_cells(qcReport.spike_table(run)['data'])['X']
    assert number(spike['Spike-in reads']) == by_id['spike.reads']['value']
    assert number(spike['% of reads into mapping']) == round(100 * by_id['spike.reads']['share'], 2)
    assert number(spike['% fragments (of spike-in reads)']) == round(100 * by_id['spike.fragment_reads']['share'], 1)


def test_checks_tsv_lists_every_check(reported):
    _, run, report, _, _, _ = reported
    with open(run / 'qc' / 'consistency_checks.tsv', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert len(rows) == len(report['checks'])
    for row, check in zip(rows, report['checks']):
        assert (row['check'], row['status'], row['sample']) == (check['id'], check['status'], check['sample'] or 'run')
        assert row['identity'] == report['check_definitions'][check['id']]


def test_only_non_equal_checks_under_the_general_statistics(reported):
    """No check section and no Run overview (Martin 2026-10-26): in the General Statistics section the summary
    line and the warnings ("none" if none); a check that differs is listed there, exactly that check."""
    import qcReport
    _, run_dir, report, _, sections, _ = reported
    assert 'ribomarker_checks' not in sections and 'ribomarker_overview' not in sections
    text = (run_dir / 'qc' / 'report_technical.html').read_text()
    note = text[text.index('id="ribomarker_general_stats"'):][:6000]
    assert 'Consistency checks:' in note and 'qc/consistency_checks.tsv' in note
    assert 'Checks that are not equal' not in note
    run = qcReport.Run(report, yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8')))
    assert ('Warnings:</b> none.' in qcReport.status_note(run)) == (not report['warnings'])
    damaged = dict(report, checks=[dict(c) for c in report['checks']])
    damaged['checks'][0]['status'] = 'differs'
    run = qcReport.Run(damaged, yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8')))
    overview = qcReport.status_note(run)
    listed = overview[overview.index('Checks that are not equal'):]
    listed = listed[:listed.index('</table>')]                               # the checks table only
    assert listed.count('<tr>') == 1 + 1                                    # header + the one check


def test_composition_tabs_leave_out_spike_ins_and_rrna(reported):
    _, _, _, _, sections, _ = reported
    section = sections['ribomarker_composition']
    labels = [d['name'] for d in section['pconfig']['data_labels']]
    assert labels == ['All RNA classes', 'Without spike-ins', 'Without rRNA and spike-ins']
    full, no_spikes, no_rrna_spikes = section['data']
    for sample, row in full.items():
        assert no_spikes.get(sample, {}) == {g: n for g, n in row.items() if g != 'spike-ins'}
        assert no_rrna_spikes.get(sample, {}) == {g: n for g, n in row.items() if g not in ('rRNA', 'spike-ins')}
    assert 'spike-ins' not in section['categories'][1] and 'rRNA' not in section['categories'][2]


def test_plots_without_series_subtitle_and_length_legend(reported):
    _, _, _, _, sections, _ = reported
    plots = [s for s in sections.values() if 'pconfig' in s]
    assert plots and all(s['pconfig']['series_label'] is False for s in plots)
    assert sections['ribomarker_lengths']['pconfig']['showlegend'] is True


def test_top20_per_rna_type(reported):
    """Recomputed from the fixture's gene table and gene biotypes: per RNA type, genes with counts > 0,
    highest first, ties by name; the overall list leaves out genes with 0 counts too."""
    _, run, _, _, sections, records = reported
    tables = {name: list(csv.reader(open(next((run / 'results').glob(f'*-{name}.csv')))))
              for name in ('genecounts', 'genelengths')}
    biotype_of = {row[0]: row[2] for row in tables['genelengths'][1:]}
    header = tables['genecounts'][0]
    for column, sample in enumerate(header[1:], start=1):
        counts = {row[0]: int(float(row[column])) for row in tables['genecounts'][1:]}
        for biotype in set(biotype_of[g] for g in counts):
            expected = sorted(((g, n) for g, n in counts.items() if biotype_of[g] == biotype and n > 0),
                              key=lambda item: (-item[1], item[0]))[:20]
            assert records[sample][f'counts.{biotype}.top20']['value'] == [list(item) for item in expected]
        assert all(n > 0 for _, n in records[sample]['counts.top20']['value'])
    content = sections['ribomarker_counts']['description']
    assert '20 most abundant features per sample (reads)' in content and 'by RNA class' in content


def test_detected_genes_are_counts_above_0_5_10(reported):
    """counts.detected.gt{n} recomputed from the fixture's gene table (S3 has genes with 0 counts and
    with 20, so > and >= give different numbers there)."""
    _, run, _, _, _, records = reported
    rows = list(csv.reader(open(next((run / 'results').glob('*-genecounts.csv')))))
    for column, sample in enumerate(rows[0][1:], start=1):
        counts = [int(float(row[column])) for row in rows[1:]]
        for above in (0, 5, 10):
            assert records[sample][f'counts.detected.gt{above}']['value'] == sum(n > above for n in counts)


def test_report_header_logo_versions_and_provenance(reported):
    """Report feedback of 2026-10-06: the title is the run id, no subtitle or intro text, the logos are
    embedded, fastp's last description sentence is gone, provenance is a section of its own directly
    above Software Versions and holds no tool versions; those are in the Software Versions table."""
    _, run, report, _, sections, _ = reported
    text = (run / 'qc' / 'report_technical.html').read_text()
    assert re.search(r'<h1[^>]*>\s*' + re.escape(report['run_id']) + r'\s*</h1>', text)
    assert 'RiboMarker QC: ' not in text and 'Technical QC of one' not in text
    assert 'class="custom_logo' in text and 'data:image/png;base64,' in text
    assert 'Notably, the latter' not in text
    provenance = sections['ribomarker_provenance']
    assert 'parent_id' not in provenance and 'Tool versions' not in provenance['data']
    order = re.findall(r'<h2 [^>]*id="([^"]+)"', text)
    assert order[-2:] == ['ribomarker_provenance', 'multiqc_software_versions']
    versions = yaml.safe_load((run / 'qc' / 'mqc' / 'ribomarker_mqc_versions.yml').read_text())
    tools = report['provenance']['prov.tool_versions']
    version = report['provenance']['prov.pipeline_version']
    details = [f"git {str(version['git_commit'])[:12]}"] if version.get('git_commit') else []   # Docker or a git checkout
    details.append(f"code {version['code_fingerprint'][:12]}")       # always: from the code files or the image
    pipeline = f"{version['version']} ({', '.join(details)})"
    assert html.escape(pipeline) in provenance['data']
    # #4 of Martin's review: the pipeline and the method under Provenance, as in the README, and in Tool Citations
    assert 'Please cite' in provenance['data'] and '10.1261/rna.080861.125' in provenance['data']
    assert 'data-doi="10.5281/zenodo.23475622"' in text and 'data-doi="10.1261/rna.080861.125"' in text
    assert versions['RiboMarker pipeline'] == {
        'RiboMarker pipeline': pipeline, **{k: str(v) for k, v in tools.items() if k not in ('cutadapt', 'fastp')}}


# ----------------------------------------------------------------------------- groups of the two plots

NC, PRE, UNANNOTATED = 'lncRNA and other ncRNA', 'genome: pre-mRNA, pre-tRNA', 'genome (unannotated)'
# Where every tier and RNA type of the human assets goes, written out by hand (review 2026-10-07): an RNA
# type is shown with the tier its features are in. Mt_rRNA and Mt_tRNA: provisional.
TIER_GROUP_SPEC = {'univec_bt1': 'UniVec', 'ribospike': 'spike-ins', 'human-rrna_bt1': 'rRNA', 'human-trnas': 'tRNA',
                   'human-mirna': 'miRNA', 'human-star': 'miRNA', 'human-pre': 'miRNA', 'human-pirna': 'piRNA',
                   'human-ncrna': NC, 'human-cdna': 'protein_coding', 'human-genome': 'genome'}   # split: tier_groups
TYPE_GROUP_SPEC = {'RiboSpike-T1': 'spike-ins', 'RiboSpike-Fragment': 'spike-ins', 'rRNA': 'rRNA', 'tRNA': 'tRNA',
                   'miRNA': 'miRNA', 'pre-miRNA': 'miRNA', 'piRNA': 'piRNA', 'UniVec': 'UniVec',
                   'lncRNA': NC, 'snRNA': NC, 'snoRNA': NC, 'scaRNA': NC, 'misc_RNA': NC, 'Mt_rRNA': NC, 'Mt_tRNA': NC,
                   'protein_coding': 'protein_coding', 'pre-mRNA': PRE, 'pre-tRNA': PRE, 'genome_unannotated': UNANNOTATED}


def test_tiers_and_rna_types_are_grouped_by_reference():
    import qcReport
    for tier, group in TIER_GROUP_SPEC.items():
        assert qcReport.tier_group(tier) == group, tier
    for biotype, group in TYPE_GROUP_SPEC.items():
        assert qcReport.type_group(biotype) == group, biotype


# one real sample (2026-10-06): reads aligned per tier (Bowtie logs) and kept reads per RNA
# class (classcounts; names translated to the class ids of 0.3.0); real data, in which every RNA class comes
# from one tier.
LEVEL_C_S3_ALIGNED = {'univec_bt1': 17079, 'ribospike': 15265, 'human-rrna_bt1': 1082639, 'human-trnas': 4611,
                      'human-mirna': 62827, 'human-star': 3133, 'human-pre': 1375, 'human-pirna': 46484,
                      'human-ncrna': 67821, 'human-cdna': 18159, 'human-genome': 41264}
# 0.3.0: the star arm is class miRNA (56,213 + 3,100, 3a) and the UniVec reads are their own class (17,079 of the
# 28,012 unannotated reads, 3g)
LEVEL_C_S3_KEPT = {'rRNA': 1082413, 'miRNA': 59313, 'UniVec': 17079, 'genome_unannotated': 10933, 'piRNA': 18726, 'lncRNA': 10187,
                   'protein_coding': 9544, 'RiboSpike-Fragment': 6772, 'snRNA': 5339, 'pre-mRNA': 4413,
                   'snoRNA': 2897, 'RiboSpike-T2': 2846, 'RiboSpike-T1': 2296, 'RiboSpike-T4': 1848, 'RiboSpike-T3': 1503,
                   'Mt_rRNA': 1437, 'tRNA': 1226, 'pre-miRNA': 522, 'Mt_tRNA': 422, 'misc_RNA': 222, 'scaRNA': 166,
                   'pre-tRNA': 35}
LEVEL_C_S3_REMOVED = 120518                       # assign.removed_by_policy


def level_c_run(kept=LEVEL_C_S3_KEPT):
    import qcReport
    metrics = [{'id': f'map.{tier}.reads_in', 'value': None} for tier in LEVEL_C_S3_ALIGNED]   # only the ids are used
    metrics += [{'id': f'map.{tier}.aligned', 'value': n} for tier, n in LEVEL_C_S3_ALIGNED.items()]
    metrics += [{'id': 'map.final_unaligned', 'value': 240644}]
    metrics += [{'id': f'assign.class.{biotype}.reads', 'value': n} for biotype, n in kept.items()]
    report = {'run_id': 'x', 'mode': 'pe', 'samples': {'S3': {'metrics': metrics, 'distributions': []}},
              'provenance': {}, 'checks': [], 'warnings': []}
    return qcReport.Run(report, yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8')))


def test_no_group_gains_reads_by_the_removal():
    """Same groups in both plots: per group, the kept reads (composition) are at most the aligned reads
    (tiers), and the differences add up to the reads removed by the multimapper/antisense settings."""
    import qcReport
    run = level_c_run()
    before = qcReport.tier_fate(run)['data']['S3']
    composition = qcReport.composition(run)
    after = composition['data'][0]['S3']
    for group, reads in after.items():
        assert reads <= before[group], group
    assert sum(before.values()) - before['unmapped'] - sum(after.values()) == LEVEL_C_S3_REMOVED
    assert after['UniVec'] == before['UniVec'] == 17079
    assert 'Groups with more reads here' not in composition['description']     # no note
    assert list(composition['categories'][0]) == qcReport.GROUP_ORDER
    assert before[UNANNOTATED] == after[UNANNOTATED] == 10933          # #9: the genome tier split by the class
    assert before[PRE] == 41264 - 10933


def test_univec_reads_are_their_own_group():
    """Since 0.3.0 (3g) the UniVec reads are the class UniVec: the composition shows them in the UniVec group
    from the class itself, and the genome group holds the genome reads without a feature."""
    import qcReport
    after = qcReport.composition(level_c_run())['data'][0]['S3']
    assert after['UniVec'] == 17079
    assert after[PRE] == 4413 + 35 and after[UNANNOTATED] == 10933     # pre-mRNA + pre-tRNA; genome_unannotated


def test_composition_names_groups_that_gain_reads(reported):
    """The synthetic bundle has snoRNA and lncRNA features in the genome region (C23, C25), so the ncRNA
    group holds more kept reads than the ncRNA tier aligned; the section names it."""
    _, _, _, _, sections, _ = reported
    description = sections['ribomarker_composition']['description']
    assert 'S1: lncRNA and other ncRNA' in description and 'S2: lncRNA and other ncRNA' in description


# ----------------------------------------------------------------------------- 0.3.0 S4 (RealSeq review)

def qcReport_palette():
    import qcReport
    return qcReport.SAMPLE_GROUP_COLOURS


def test_sample_table_from_the_users_sheet(reported):
    """#1: the user's sheet (semicolons, byte-order mark, own column order) in the run's sample order; the FASTQ
    columns hidden at first; a 'Colour samples by' selector with every other column (report.js colours the samples
    through MultiQC's highlighting; Martin 2026-10-26: chosen in the report, report_group is the first choice)."""
    _, run, report, _, sections, _ = reported
    table = sections['ribomarker_samples']
    assert list(table['data']) == ['S1', 'S2', 'S3']
    assert table['data']['S2'] == {'group': 'case', 'fastq_file': 'S2_R1_001.fastq.gz',
                                   'fastq_file2': 'S2_R2_001.fastq.gz', 'age': ''}
    assert table['headers']['fastq_file']['hidden'] and not table['headers']['age']['hidden']
    assert report['provenance']['prov.sample_sheet']['report_group'] == 'group'
    selector = re.search(r'<div class="rm-colour-by" data-config="([^"]+)"></div>', table['description'])
    config = json.loads(base64.b64decode(selector.group(1)))
    assert config['samples'] == ['S1', 'S2', 'S3'] and config['first'] == 'group'
    assert config['columns'] == {'group': ['control', 'case', 'control'], 'age': ['41', '(empty)', '38']}
    assert config['palette'] == qcReport_palette() and 'ribomarker_counts_plot' in config['series_plots']
    # the colours are set in the report only: no static line colours, legend groups or highlight patterns
    assert 'colors' not in sections['ribomarker_lengths']['pconfig']
    assert 'legend_groups' not in sections['ribomarker_lengths']['pconfig']
    assert 'sample_groups' not in sections['ribomarker_read_fate']['pconfig']
    text = (run / 'qc' / 'report_technical.html').read_text()
    assert '"highlight_patterns": []' in text                                   # none from the configuration
    assert 'initColourSelector' in text                                         # report.js inlined
    assert 'dataset-switch-group { flex-wrap: wrap; }' in text                  # report.css: buttons wrap


def test_without_a_sheet_the_samples_are_listed(tmp_path):
    """Configuration route without sample_sheet: the sample list (no groups, no warning)."""
    from conftest import prepare_run, run_collector
    import qcReport
    run, config_file = prepare_run(tmp_path, 'pe')
    config = yaml.safe_load(config_file.read_text())
    del config['sample_sheet'], config['report_group']
    config_file.write_text(yaml.safe_dump(config))
    report = run_collector(run, config_file)
    assert 'prov.sample_sheet' not in report['provenance'] and not report['warnings'][1:]
    r = qcReport.Run(report, yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8')))
    table = qcReport.samples_table(r)
    assert list(table['data']) == ['S1', 'S2', 'S3'] and table['data']['S1']['fastq_file'] == 'S1_R1_001.fastq.gz'
    assert 'rm-colour-by' not in table['description']                           # no sheet columns: no selector


def test_a_group_column_missing_from_the_sheet_is_a_warning(tmp_path):
    from conftest import prepare_run, run_collector
    run, config_file = prepare_run(tmp_path, 'se')
    config = yaml.safe_load(config_file.read_text())
    config['report_group'] = 'condition'
    config_file.write_text(yaml.safe_dump(config))
    report = run_collector(run, config_file)
    assert report['provenance']['prov.sample_sheet']['report_group'] is None
    assert any("'condition' is not in" in w['message'] for w in report['warnings'])


def test_detected_genes_per_class_plot(reported):
    """#18: one bar per sample and class, the collected counts.detected.{class}.gt{n}."""
    import refTables
    _, _, _, long, sections, _ = reported
    section = sections['ribomarker_counts']                      # one section with the counts (Martin 2026-10-26)
    classes = refTables.ClassTable(refTables.default_table('human'), 'human')
    assert [d['name'] for d in section['pconfig']['data_labels']] == ['> 0 reads', '> 5 reads', '> 10 reads']
    assert section['pconfig']['logswitch'] is True
    assert 'Features &gt; 5 reads' in section['description']     # the table uses the same thresholds as the plot
    for above, dataset in zip((0, 5, 10), section['data']):
        for name, row in dataset.items():
            rna_class = next(c for c in classes._by_class if classes.display_name(c) == name) \
                if name in {classes.display_name(c) for c in classes._by_class} else name
            for sample, n in row.items():
                assert n == long[(sample, f'counts.detected.{rna_class}.gt{above}')]


def test_top_lists_only_as_long_as_needed(reported):
    """#17: a top list has as many ranks as the longest list of a sample."""
    import qcReport
    _, _, report, _, _, records = reported
    run = qcReport.Run(report, yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8')))
    for metric_id in {i for r in records.values() for i in r if i.endswith('.top20')}:
        longest = max(len(records[s].get(metric_id, {}).get('value') or []) for s in records)
        assert qcReport.top_table(run, metric_id).count('<tr>') == 1 + longest, metric_id


def test_guide_tables_and_class_names(reported):
    """#3: a guide; the tier list and class table in the report; class names, not ids, in the buttons;
    #8: protein_coding in the tier legend, explained; #2: "RNA class" in the titles."""
    _, run, _, _, sections, _ = reported
    assert 'Toolbox' in sections['ribomarker_guide']['data'] and 'Configure columns' in sections['ribomarker_guide']['data']
    assert 'RNA classes (class table)' in sections['ribomarker_tables']['data']
    tiers = sections['ribomarker_tiers']
    assert 'protein_coding' in tiers['description'] and 'UTR' in tiers['description']
    assert 'mRNA' not in [c['name'] for c in tiers['categories'].values()]
    text = (run / 'qc' / 'report_technical.html').read_text()
    for title in ('Composition by RNA class', 'Read length by RNA class', 'Detection and counts'):
        assert title in text
    # settings: report.js inlined once; the native toggles it drives are in the report
    assert text.count('rm-settings-line') >= 2 and 'percent-switch' in text and 'log10-switch' in text
    assert 'RNA type' not in text.replace('RNA Type', '')


def test_samples_first_then_guide_then_general_statistics(reported):
    """Martin 2026-10-26: Samples, How to read this report, General Statistics, then the other sections; MultiQC's
    own General Statistics (fixed under the header) is turned off; the check status is one line in Samples."""
    _, run_dir, report, _, sections, _ = reported
    text = (run_dir / 'qc' / 'report_technical.html').read_text()
    positions = [text.index(f'id="{anchor}"') for anchor in
                 ('ribomarker_samples', 'ribomarker_guide', 'ribomarker_general_stats', 'ribomarker_read_fate')]
    assert positions == sorted(positions)
    assert 'id="general_stats"' not in text                # MultiQC's own section (templates/default/general_stats.html)
    assert sections['ribomarker_general_stats']['plot_type'] == 'table'
    line = sections['ribomarker_samples']['description']
    assert line.startswith('<p><b>Consistency checks:</b>')
    assert ('<b>Warnings:</b> none.' in line) == (not report['warnings'])


def test_natural_reads_leave_out_every_spike_in():
    """Natural reads: every kept class except the spike-ins (RiboSpike-* and miND; critique 2026-10-27: miND
    counted as natural before), UniVec and unmapped."""
    import types
    import qcCollect
    import refTables
    c = types.SimpleNamespace(classes=refTables.ClassTable(refTables.default_table('human'), 'human'))
    classes = ['miRNA', 'miND', 'RiboSpike-T1', 'RiboSpike-Fragment', 'UniVec', 'unmapped', 'genome_unannotated']
    assert qcCollect.natural_classes(c, classes) == ['miRNA', 'genome_unannotated']


def test_report_scripts_match_the_pinned_multiqc():
    """report.js uses MultiQC's report scripts (critique 2026-10-27): the version it was written for is the one
    the environment pins; a MultiQC update has to check report.js and change MULTIQC_VERSION."""
    import qcReport
    lock = (REPO / 'docker' / 'env.lock.txt').read_text()
    pinned = re.search(r'/multiqc-([0-9.]+)-', lock).group(1)
    assert pinned == qcReport.MULTIQC_VERSION
    script = (REPO / 'docker' / 'qc' / 'report.js').read_text()
    for name in ('apply_mqc_highlights', 'mqc_plots', 'make_colorsamples_filter'):
        assert name in script
