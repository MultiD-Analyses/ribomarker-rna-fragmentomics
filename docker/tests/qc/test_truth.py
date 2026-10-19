"""The collector on the synthetic level A data: every metric that the truth defines
must equal it exactly, in both modes."""
import csv
import statistics

from conftest import REPO, values

TRUTH_DIR = REPO / 'test' / 'synthetic_qc' / 'bundle' / 'truth'


def read_truth(mode, sample):
    """({(table, key): expected}, {category: reads}) for one mode and sample."""
    with open(TRUTH_DIR / 'truth.tsv', newline='') as f:
        truth = {(row['table'], row['key']): int(row['expected'])
                 for row in csv.DictReader(f, delimiter='\t') if row['mode'] == mode and row['sample'] == sample}
    with open(TRUTH_DIR / 'categories.tsv', newline='') as f:
        categories = {row['id']: int(row[f'n_{sample}']) for row in csv.DictReader(f, delimiter='\t')}
    return truth, categories


def expected_metrics(mode, sample, tiers):
    """Metric values that follow from the synthetic design (categories.tsv / truth.tsv).
    Rows of table "qc" in truth.tsv are QC metric ids with their expected value."""
    truth, n = read_truth(mode, sample)
    expected = {key: value for (table, key), value in truth.items() if table == 'qc'}
    if mode == 'pe':
        expected['trim.pairs_in'] = truth[('input', 'reads')]
        expected['trim.adapter_dimer_pairs'] = n['C01']                  # empty inserts
        expected['trim.short_insert_pairs'] = n['C02']                   # 10-nt inserts
        expected['trim.pairs_too_short'] = n['C01'] + n['C02']
        expected['trim.pairs_written'] = truth[('trimstats', 'trimmed_reads')]
        expected['trim.merged_reads'] = truth[('trimstats', 'merged_reads')]
        expected['map.total_in'] = truth[('trimstats', 'merged_reads')]
    else:
        expected['trim.reads_in'] = truth[('input', 'reads')]
        expected['trim.adapter_dimer_reads'] = n['C01']
        expected['trim.short_insert_reads'] = n['C02']
        expected['trim.reads_too_short'] = truth.get(('trimstats', '<15nt_reads'), 0)   # S3 has none
        expected['trim.reads_written'] = truth[('trimstats', 'trimmed_reads')]
        expected['map.total_in'] = truth[('trimstats', 'trimmed_reads')]
    if ('mappingstats', 'Total') not in truth:           # exploratory sample (S3): mapping is not checked
        return expected

    for tier in tiers:                                   # a tier without a truth row gets no reads
        expected[f'map.{tier}.aligned'] = truth.get(('tier_aligned', tier), 0)
    expected['map.final_unaligned'] = truth[('tier_aligned', 'final_unaligned')]

    for column in ('Mapped', 'Unmapped', 'Unique', 'Multi'):
        expected[f'assign.{column.lower()}'] = truth[('mappingstats', column)]
    type_counts = {key: value for (table, key), value in truth.items() if table == 'typecounts'}
    kept = sum(value for key, value in type_counts.items() if key != 'unmapped')
    expected['assign.kept'] = kept
    expected['assign.removed_by_policy'] = expected['assign.removed_multimap'] + expected['assign.removed_antisense']
    spikes = sum(value for key, value in type_counts.items() if key.startswith('RiboSpike'))
    expected['assign.annotated'] = kept - type_counts.get('unannotated', 0) - spikes
    for key, value in type_counts.items():
        if key != 'unmapped':
            expected[f'assign.type.{key}.reads'] = value

    expected['spike.reads'] = truth.get(('tier_aligned', 'ribospike'), 0)
    for (table, key), value in truth.items():
        if table == 'ribospike':
            _, length, spike_type = key.split('-')
            expected[f'spike.{length}_{spike_type}.reads'] = value

    gene_counts = [value for (table, _), value in truth.items() if table == 'genecounts']
    expected['counts.library_size'] = sum(gene_counts)
    expected['counts.detected.gt0'] = len(gene_counts)
    expected['counts.top10_reads'] = sum(sorted(gene_counts, reverse=True)[:10])
    with open(TRUTH_DIR / 'categories.tsv', newline='') as f:
        rrna_genes = {row['gene'] for row in csv.DictReader(f, delimiter='\t') if row['typecount'] == 'rRNA'}
    non_rrna = sorted((value for (table, key), value in truth.items() if table == 'genecounts' and key not in rrna_genes),
                      reverse=True)
    expected['counts.library_size_non_rrna'] = sum(non_rrna)
    expected['counts.top10_reads_non_rrna'] = sum(non_rrna[:10])
    return expected


def expected_aligned_multi(mode, sample, found):
    """map.{tier}.aligned_multi (Bowtie 2 tiers only, so only where the metric exists)."""
    truth, _ = read_truth(mode, sample)
    if ('mappingstats', 'Total') not in truth:
        return {}
    return {metric_id: truth.get(('tier_aligned_multi', metric_id.split('.')[1]), 0)
            for (s, metric_id) in found if s == sample and metric_id.endswith('.aligned_multi')}


def expected_end_base_composition(sample):
    """{length: (largest |base share - 0.25|, n)} over the 4 first and 4 last bases of the
    exact-length spike reads (insert length = designed length), from reads_truth.tsv."""
    shares = {}
    with open(TRUTH_DIR / 'reads_truth.tsv', newline='') as f:
        rows = [r for r in csv.DictReader(f, delimiter='\t')
                if r['sample'] == sample and r['source_contig'].startswith('SR_')]
    for length in (20, 30, 40, 50, 60):
        inserts = [r['insert_sequence'] for r in rows
                   if int(r['source_contig'].split('_')[1]) == length and int(r['insert_length']) == length]
        if not inserts:
            continue
        positions = [[s[i] for s in inserts] for i in list(range(4)) + list(range(length - 4, length))]
        largest = max(abs(column.count(base) / len(inserts) - 0.25) for column in positions for base in 'ACGT')
        shares[length] = (largest, len(inserts))
    return shares


def test_metrics_equal_truth(collected):
    mode, report = collected
    found = values(report)
    for sample in report['samples']:
        tiers = [record['id'].split('.')[1] for record in report['samples'][sample]['metrics']
                 if record['id'].startswith('map.') and record['id'].endswith('.reads_in')]
        assert len(tiers) == 11
        mismatches = []
        expected = expected_metrics(mode, sample, tiers) | expected_aligned_multi(mode, sample, found)
        for metric_id, value in expected.items():
            if found.get((sample, metric_id)) != value:
                mismatches.append(f'{sample} {metric_id}: expected {value}, found {found.get((sample, metric_id))}')
        assert not mismatches, '\n'.join(mismatches)


def test_end_base_composition_on_exact_reads(collected):
    """The exact-read end composition equals the one computed from the designed inserts
    (only spike reads of the designed length; the -1/-2 nt reads of C26/C27 are excluded),
    and the chance level follows from their number."""
    _, report = collected
    found = values(report)
    for sample in report['samples']:
        for length, (largest, n) in expected_end_base_composition(sample).items():
            assert abs(found[(sample, f'spike.{length}.end_base_composition')] - largest) < 1e-9
            assert found[(sample, f'spike.{length}.exact_reads')] == n
            assert abs(found[(sample, f'spike.{length}.end_base_chance')] - 2.31 * (0.1875 / n) ** 0.5) < 1e-12


def test_spike_offset_histograms_equal_truth(collected):
    """spike.{len}.length_offset_hist: full-length spike reads at observed - designed length -2..+2, summed
    over the four types, from the truth's ribolengths rows (C04: offset 0, C26: -1, C27: -2)."""
    mode, report = collected
    found = values(report)
    with open(TRUTH_DIR / 'truth.tsv', newline='') as f:
        rows = [row for row in csv.DictReader(f, delimiter='\t') if row['mode'] == mode and row['table'] == 'ribolengths']
    offsets = [-2, -1, 0, 1, 2]
    checked = 0
    for sample in report['samples']:
        reads_at = {}                                    # observed length -> reads, all spike types
        for row in rows:
            if row['sample'] == sample:
                length = int(row['key'].split(':')[1])
                reads_at[length] = reads_at.get(length, 0) + int(row['expected'])
        for designed in (20, 30, 40, 50, 60):
            expected = [reads_at.get(designed + offset, 0) for offset in offsets]
            assert found[(sample, f'spike.{designed}.length_offset_hist')] == (offsets, expected), f'{sample} {designed}'
            checked += sum(expected)
    assert checked > 0                                   # the truth has spike reads (S1)


def kept_reads(mode, sample):
    """{RNA type: [insert sequences]} of the reads the design keeps (mapped and not removed by the
    multimapper/antisense settings), from reads_truth.tsv and categories.tsv; spike-ins left out."""
    with open(TRUTH_DIR / 'categories.tsv', newline='') as f:
        categories = {row['id']: row for row in csv.DictReader(f, delimiter='\t')}
    fate = 'fate_paired_end' if mode == 'pe' else 'fate_single_end'
    kept = {}
    with open(TRUTH_DIR / 'reads_truth.tsv', newline='') as f:
        for read in csv.DictReader(f, delimiter='\t'):
            category = categories[read['category']]
            if (read['sample'] == sample and read[fate] == 'mapped' and not category['removed']
                    and category['typecount'] not in ('', 'unmapped')):
                kept.setdefault(category['typecount'], []).append(read['insert_sequence'])
    return kept


def test_length_and_sequence_metrics_equal_truth(collected):
    """Per RNA type: reads of 15-17 nt (floor_reads), median length and distinct sequences of the kept
    reads, computed from the designed inserts (samples with a mapping truth: S1, S2)."""
    mode, report = collected
    found = values(report)
    for sample in report['samples']:
        truth, _ = read_truth(mode, sample)
        if ('mappingstats', 'Total') not in truth:       # exploratory sample (S3)
            continue
        for biotype, sequences in kept_reads(mode, sample).items():
            lengths = [len(sequence) for sequence in sequences]
            assert len(sequences) == truth[('typecounts', biotype)], f'{sample} {biotype}: kept reads'
            assert found[(sample, f'assign.type.{biotype}.floor_reads')] == sum(15 <= n <= 17 for n in lengths), biotype
            assert found[(sample, f'assign.type.{biotype}.median_length')] == statistics.median(lengths), biotype
            assert found[(sample, f'assign.type.{biotype}.distinct_sequences')] == len(set(sequences)), biotype


def test_only_expected_warning(collected):
    """The fixture has no reference MANIFEST; nothing else may warn."""
    _, report = collected
    assert [w['source'] for w in report['warnings']] == ['<dbdir>/MANIFEST_human.json']
