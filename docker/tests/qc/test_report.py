"""The technical report: MultiQC renders it from the collector's output; the
numbers in the RiboMarker sections are the collected values; no absolute path reaches the HTML."""
import csv
import html
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
    data = sections['ribomarker_read_fate']['data']
    reads_in = 'trim.pairs_in' if mode == 'pe' else 'trim.reads_in'
    assert list(data) == ['S1', 'S2', 'S3']                 # sample-sheet order
    for sample, row in data.items():
        for metric_id, value in row.items():
            assert value == long[(sample, metric_id)]
        assert sum(row.values()) == long[(sample, reads_in)]   # categories add up to the input


def test_tiers_and_composition_add_up(reported):
    _, _, _, long, sections, records = reported
    for sample, row in sections['ribomarker_tiers']['data'].items():
        assert sum(row.values()) == long[(sample, 'map.total_in')]
    all_types, no_spikes, no_rrna_spikes = sections['ribomarker_composition']['data']
    for sample in all_types:
        assert sum(all_types[sample].values()) == long[(sample, 'assign.kept')]
        assert all_types[sample].get('UniVec') == long[(sample, 'map.univec_bt1.aligned')]
        # the rRNA group is the rRNA tier; Mt_rRNA is in the ncRNA group for now
        rrna = sum(r['value'] for i, r in records[sample].items() if i == 'assign.type.rRNA.reads')
        spikes = sum(r['value'] for i, r in records[sample].items()
                     if i.startswith('assign.type.RiboSpike') and i.endswith('.reads'))
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
           'unaligned': ('% unaligned (of reads to mapping)', 'map.final_unaligned', True),
           'removed': ('% removed by policy (of mapped)', 'assign.removed_by_policy', True),
           'library': ('Library size', 'counts.library_size', False),
           'genes10': ('Genes > 10', 'counts.detected.gt10', False)},
    'se': {'reads_in': ('Reads in', 'trim.reads_in', False),
           'short': ('% < 15 nt (of reads in)', 'trim.reads_too_short', True),
           'mapping_in': ('Reads to mapping', 'map.total_in', False),
           'unaligned': ('% unaligned (of reads to mapping)', 'map.final_unaligned', True),
           'removed': ('% removed by policy (of mapped)', 'assign.removed_by_policy', True),
           'library': ('Library size', 'counts.library_size', False),
           'genes10': ('Genes > 10', 'counts.detected.gt10', False)},
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
        if 'mapped' in title.lower():   # "mapped" in a title must mean assign.mapped
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
    counts = table_cells(sections['ribomarker_counts']['data'])
    for sample, metrics in records.items():
        if metrics.get('spike.reads', {}).get('value'):
            row = spike[sample]
            assert number(row['Spike-in reads']) == metrics['spike.reads']['value']
            assert number(row['% of reads into mapping']) == round(100 * metrics['spike.reads']['share'], 2)
            assert number(row['% fragments (of spike-in reads)']) == round(100 * metrics['spike.fragment_reads']['share'], 1)
        row = counts[sample]
        for header, metric_id in [('Library size', 'counts.library_size'), ('Without rRNA', 'counts.library_size_non_rrna'),
                                  ('Genes > 0', 'counts.detected.gt0'), ('Genes > 5', 'counts.detected.gt5'),
                                  ('Genes > 10', 'counts.detected.gt10')]:
            assert number(row[header]) == metrics[metric_id]['value'], f'{sample} {header}'
        for header, metric_id in [('% in top-10 genes', 'counts.top10_reads'),
                                  ('% in top-10 genes without rRNA', 'counts.top10_reads_non_rrna')]:
            assert number(row[header]) == round(100 * metrics[metric_id]['share'], 1), f'{sample} {header}'


def test_spike_sections_equal_collected_values(reported):
    _, _, report, long, sections, _ = reported
    recovery = sections['ribomarker_spike_recovery']['data']
    for sample, row in recovery.items():
        for name, value in row.items():
            assert value == long[(sample, f"spike.{name.split()[0]}.reads")]
    for dataset in sections['ribomarker_spike_offsets']['data']:
        for sample, row in dataset.items():
            assert sum(row.values()) > 0
    assert 'S2' not in recovery and 'S3' not in recovery     # S2 and S3 have no spike-in reads


def test_length_curves_are_shares_of_the_type(reported):
    _, _, report, _, sections, _ = reported
    labels = [label['name'] for label in sections['ribomarker_lengths']['pconfig']['data_labels']]
    assert labels[0] == 'All RNA types'
    for biotype, dataset in zip(labels[1:], sections['ribomarker_lengths']['data'][1:]):
        for sample, curve in dataset.items():
            dist = next(d for d in report['samples'][sample]['distributions']
                        if d['id'] == f'assign.type.{biotype}.length_hist')
            for x, share in curve.items():
                assert share == pytest.approx(100 * dist['y'][dist['x'].index(int(x))] / dist['n'])


def test_all_types_curve_pools_the_type_counts(reported):
    _, _, report, _, sections, _ = reported
    labels = [label['name'] for label in sections['ribomarker_lengths']['pconfig']['data_labels']]
    pooled = sections['ribomarker_lengths']['data'][0]
    for sample, curve in pooled.items():
        dists = [d for d in report['samples'][sample]['distributions']
                 if d['id'] in {f'assign.type.{b}.length_hist' for b in labels[1:]}]
        total = sum(d['n'] for d in dists)
        for x, share in curve.items():
            count = sum(d['y'][d['x'].index(int(x))] for d in dists)
            assert share == pytest.approx(100 * count / total)


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
           'assign.mapped', 'assign.removed_by_policy', 'counts.library_size', 'counts.library_size_non_rrna',
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

    counts = table_cells(qcReport.counts_table(run)['data'])['X']
    for header, metric_id in [('Library size', 'counts.library_size'), ('Without rRNA', 'counts.library_size_non_rrna'),
                              ('Genes > 0', 'counts.detected.gt0'), ('Genes > 5', 'counts.detected.gt5'),
                              ('Genes > 10', 'counts.detected.gt10')]:
        assert number(counts[header]) == by_id[metric_id]['value'], header
    for header, metric_id in [('% in top-10 genes', 'counts.top10_reads'),
                              ('% in top-10 genes without rRNA', 'counts.top10_reads_non_rrna')]:
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


def test_only_non_equal_checks_in_the_overview(reported):
    """No separate check section. All checks equal on the fixtures: the overview has only the summary
    line. One check that differs: the overview lists exactly that check."""
    import qcReport
    _, _, report, _, sections, _ = reported
    assert 'ribomarker_checks' not in sections
    overview = sections['ribomarker_overview']['data']
    assert 'qc/consistency_checks.tsv' in overview and 'Checks that are not equal' not in overview
    damaged = dict(report, checks=[dict(c) for c in report['checks']])
    damaged['checks'][0]['status'] = 'differs'
    run = qcReport.Run(damaged, yaml.safe_load((REPO / 'docker' / 'qc' / 'metrics.yaml').read_text(encoding='utf-8')))
    overview = qcReport.overview(run)['data']
    listed = overview[overview.index('Checks that are not equal'):]
    listed = listed[:listed.index('</table>')]                               # the checks table only
    assert listed.count('<tr>') == 1 + 1                                    # header + the one check


def test_composition_tabs_leave_out_spike_ins_and_rrna(reported):
    _, _, _, _, sections, _ = reported
    section = sections['ribomarker_composition']
    labels = [d['name'] for d in section['pconfig']['data_labels']]
    assert labels == ['All RNA types', 'Without spike-ins', 'Without rRNA and spike-ins']
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
    content = sections['ribomarker_counts']['data']
    assert '20 most abundant RNAs per sample (reads)' in content and 'by RNA type' in content


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
    assert html.escape(pipeline) in sections['ribomarker_overview']['data']
    assert versions['RiboMarker pipeline'] == {
        'RiboMarker pipeline': pipeline, **{k: str(v) for k, v in tools.items() if k not in ('cutadapt', 'fastp')}}


# ----------------------------------------------------------------------------- groups of the two plots

NC, GENOME = 'lncRNA and other ncRNA', 'genome, precursors, introns'
# Where every tier and RNA type of the human assets goes, written out by hand (review 2026-10-07): an RNA
# type is shown with the tier its features are in. Mt_rRNA and Mt_tRNA: provisional.
TIER_GROUP_SPEC = {'univec_bt1': 'UniVec', 'ribospike': 'spike-ins', 'human-rrna_bt1': 'rRNA', 'human-trnas': 'tRNA',
                   'human-mirna': 'miRNA', 'human-star': 'miRNA', 'human-pre': 'miRNA', 'human-pirna': 'piRNA',
                   'human-ncrna': NC, 'human-cdna': 'mRNA', 'human-genome': GENOME}
TYPE_GROUP_SPEC = {'RiboSpike-T1': 'spike-ins', 'RiboSpike-Short': 'spike-ins', 'rRNA': 'rRNA', 'tRNA': 'tRNA',
                   'miRNA': 'miRNA', 'miRNA-star': 'miRNA', 'miRNA-pre': 'miRNA', 'piRNA': 'piRNA',
                   'lncRNA': NC, 'snRNA': NC, 'snoRNA': NC, 'scaRNA': NC, 'misc_RNA': NC, 'Mt_rRNA': NC, 'Mt_tRNA': NC,
                   'protein_coding': 'mRNA', 'pre-mRNA': GENOME, 'pre-tRNA': GENOME, 'unannotated': GENOME}


def test_tiers_and_rna_types_are_grouped_by_reference():
    import qcReport
    for tier, group in TIER_GROUP_SPEC.items():
        assert qcReport.tier_group(tier) == group, tier
    for biotype, group in TYPE_GROUP_SPEC.items():
        assert qcReport.type_group(biotype) == group, biotype


# one real sample (2026-10-06): reads aligned per tier (Bowtie logs) and kept reads per RNA
# type (typecounts); real data, in which every RNA type comes from one tier.
LEVEL_C_S3_ALIGNED = {'univec_bt1': 17079, 'ribospike': 15265, 'human-rrna_bt1': 1082639, 'human-trnas': 4611,
                      'human-mirna': 62827, 'human-star': 3133, 'human-pre': 1375, 'human-pirna': 46484,
                      'human-ncrna': 67821, 'human-cdna': 18159, 'human-genome': 41264}
LEVEL_C_S3_KEPT = {'rRNA': 1082413, 'miRNA': 56213, 'unannotated': 28012, 'piRNA': 18726, 'lncRNA': 10187,
                   'protein_coding': 9544, 'RiboSpike-Short': 6772, 'snRNA': 5339, 'pre-mRNA': 4413, 'miRNA-star': 3100,
                   'snoRNA': 2897, 'RiboSpike-T2': 2846, 'RiboSpike-T1': 2296, 'RiboSpike-T4': 1848, 'RiboSpike-T3': 1503,
                   'Mt_rRNA': 1437, 'tRNA': 1226, 'miRNA-pre': 522, 'Mt_tRNA': 422, 'misc_RNA': 222, 'scaRNA': 166,
                   'pre-tRNA': 35}
LEVEL_C_S3_REMOVED = 120518                       # assign.removed_by_policy


def level_c_run(kept=LEVEL_C_S3_KEPT):
    import qcReport
    metrics = [{'id': f'map.{tier}.reads_in', 'value': None} for tier in LEVEL_C_S3_ALIGNED]   # only the ids are used
    metrics += [{'id': f'map.{tier}.aligned', 'value': n} for tier, n in LEVEL_C_S3_ALIGNED.items()]
    metrics += [{'id': 'map.final_unaligned', 'value': 240644}]
    metrics += [{'id': f'assign.type.{biotype}.reads', 'value': n} for biotype, n in kept.items()]
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
    assert sum(before.values()) - before['not aligned'] - sum(after.values()) == LEVEL_C_S3_REMOVED
    assert after['UniVec'] == before['UniVec'] == 17079
    assert '<p>' not in composition['description']                 # no note
    assert list(composition['categories'][0]) == qcReport.GROUP_ORDER


def test_unannotated_not_split_when_univec_does_not_fit():
    """If the UniVec tier holds more reads than 'unannotated', the split cannot be right: all 'unannotated'
    reads stay in the genome group and the section says so."""
    import qcReport
    composition = qcReport.composition(level_c_run(dict(LEVEL_C_S3_KEPT, unannotated=100)))
    after = composition['data'][0]['S3']
    assert 'UniVec' not in after
    assert after[GENOME] == 4413 + 35 + 100                          # pre-mRNA + pre-tRNA + unannotated
    assert 'Not split' in composition['description'] and 'S3' in composition['description']


def test_composition_names_groups_that_gain_reads(reported):
    """The synthetic bundle has snoRNA and lncRNA features in the genome region (C23, C25), so the ncRNA
    group holds more kept reads than the ncRNA tier aligned; the section names it."""
    _, _, _, _, sections, _ = reported
    description = sections['ribomarker_composition']['description']
    assert 'S1: lncRNA and other ncRNA' in description and 'S2: lncRNA and other ncRNA' in description
