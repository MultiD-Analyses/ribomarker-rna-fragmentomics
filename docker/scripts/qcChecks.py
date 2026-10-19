"""
qcChecks.py -- consistency checks and computability rules for qcCollect.py.

  Consistency checks C1-C10: exact
  identities between pipeline stages, tolerance 0 (register B1). Each check gets one
  status: 'equal', 'differs', 'not evaluable' (an input is missing) or
  'not applicable' (not defined for the run's mode).

  Kind A rules (section 9.2, A1-A3): every share gets its n (the denominator) and a
  flag 'not evaluable' (n = 0 or denominator missing) or 'low n' (n < 100); the same
  for distributions and for medians / modes / deviations computed from n reads.

Pure functions on values that qcCollect.py has already read; no file access here.
"""

LOW_N = 100          # register A1 (proposed default): below this a share is marked "low n"

CHECKS = {
    'C1': 'trim.pairs_written = fastp read-1 input count',
    'C2': 'reads written by trimming (pe: trim.merged_reads, se: trim.reads_written) = map.total_in',
    'C3': "map.{tier}.unaligned = the next tier's map.{tier}.reads_in",
    'C4': 'map.final_unaligned = assign.unmapped',
    'C5': 'assign.mapped + assign.unmapped = reads written by trimming',
    'C6': 'assign.removed_by_policy = mapping_total - reads_total in metrics/<sample>.log',
    'C7': 'counts.library_size = assign.kept - reads of pre-mRNA, pre-tRNA, miRNA-pre, '
          'unannotated and RiboSpike-* (typecounts)',
    'C8': 'the reference BED contains no contigs named for another species',
    'C9': 'every result table contains exactly the samples of the sample sheet',
    'C10': 'assign.unique + assign.multi + spike.reads = assign.mapped',
}

# Contig name prefixes per species in <species>.bed (tRNA/rRNA, miRNA, piRNA, Ensembl transcripts).
# Contigs without one of these prefixes (genome chromosomes) are not judged.
SPECIES_PREFIXES = {
    'human': ('Homo_sapiens_', 'Hsa-', 'piR-hsa-', 'ENST'),
    'mouse': ('Mus_musculus_', 'Mmu-', 'piR-mmu-', 'ENSMUST'),
}

NOT_GENES = ('pre-mRNA', 'pre-tRNA', 'miRNA-pre', 'unannotated')   # dropped by collapseGenes.py


def check(check_id, expected, observed, sample=None, subject=None):
    """One check record; 'not evaluable' if either side is missing."""
    if expected is None or observed is None:
        status = 'not evaluable'
    else:
        status = 'equal' if expected == observed else 'differs'
    numbers = isinstance(expected, int) and isinstance(observed, int)
    return {'id': check_id, 'sample': sample, 'subject': subject, 'expected': expected, 'observed': observed,
            'difference': observed - expected if numbers else None, 'status': status}


SPIKE_OFF = 'spike-in counting switched off (count_ribospike: False)'
GENES_OFF = 'no gene table (remove_multimap or remove_antisense False)'


def not_applicable(check_id, sample, reason):
    return {'id': check_id, 'sample': sample, 'subject': reason, 'expected': None, 'observed': None,
            'difference': None, 'status': 'not applicable'}


def add(*numbers):
    """Sum, or None if any part is missing."""
    return None if any(n is None for n in numbers) else sum(numbers)


def subtract(a, b):
    return None if a is None or b is None else a - b


def sample_checks(sample, values, facts, mode, tiers, spike_counted=True, genes_counted=True):
    """Checks C1-C7 and C10 for one sample. values = {metric id: value}; facts = numbers read
    for checks only (fastp read-1 input, metrics log totals); tiers = [{'tier', 'retain'}] in order;
    spike_counted = False when count_ribospike is off (no spike-in tables, so C10 is not defined);
    genes_counted = False when there is no gene table (C7 is not defined)."""
    get = values.get
    written = get('trim.merged_reads') if mode == 'pe' else get('trim.reads_written')
    checks = []

    if mode == 'pe':
        checks.append(check('C1', get('trim.pairs_written'), facts.get('fastp_read1_in'), sample))
    else:
        checks.append(not_applicable('C1', sample, 'no fastp step in single-end mode'))

    checks.append(check('C2', written, get('map.total_in'), sample))

    if len(tiers) < 2:
        checks.append(check('C3', None, None, sample, 'tier order unknown (logs/map.log)'))
    for upper, lower in zip(tiers, tiers[1:]):
        subject = f"{upper['tier']} -> {lower['tier']}"
        if upper['retain'] != 'T':       # reads aligned in this tier are passed on, too
            checks.append(not_applicable('C3', sample, f'{subject}: tier does not remove its reads'))
        else:
            checks.append(check('C3', get(f"map.{upper['tier']}.unaligned"), get(f"map.{lower['tier']}.reads_in"),
                                sample, subject))

    checks.append(check('C4', get('map.final_unaligned'), get('assign.unmapped'), sample))
    checks.append(check('C5', written, add(get('assign.mapped'), get('assign.unmapped')), sample))
    checks.append(check('C6', subtract(facts.get('mapping_total'), facts.get('reads_total')),
                        get('assign.removed_by_policy'), sample))

    # typecounts: a category without a row had no reads (0 is a measurement there)
    not_genes = sum(value for metric_id, value in values.items()
                    if metric_id.startswith('assign.type.') and metric_id.endswith('.reads')
                    and (metric_id.split('.')[2] in NOT_GENES or metric_id.split('.')[2].startswith('RiboSpike')))
    if genes_counted:
        checks.append(check('C7', subtract(get('assign.kept'), not_genes), get('counts.library_size'), sample))
    else:
        checks.append(not_applicable('C7', sample, GENES_OFF))

    if spike_counted:
        checks.append(check('C10', get('assign.mapped'),
                            add(get('assign.unique'), get('assign.multi'), get('spike.reads')), sample))
    else:
        checks.append(not_applicable('C10', sample, SPIKE_OFF))
    return checks


def species_check(contigs, species):
    """C8: number of BED contigs named for another species (expected 0)."""
    if species not in SPECIES_PREFIXES:
        return not_applicable('C8', None, f'no contig name prefixes known for species {species!r}')
    if contigs is None:
        return check('C8', 0, None, None, f'reference BED of species {species!r}')
    foreign = {}
    for other, prefixes in SPECIES_PREFIXES.items():
        if other == species:
            continue
        for prefix in prefixes:
            n = sum(1 for contig in contigs if contig.startswith(prefix))
            if n:
                foreign[prefix] = n
    subject = f'reference BED of species {species!r}'
    if foreign:
        subject += '; contigs of another species: ' + ', '.join(f'{p}* {n}' for p, n in foreign.items())
    return check('C8', 0, sum(foreign.values()), None, subject)


def sample_set_checks(sample_sets, samples):
    """C9: one check per result table; sample_sets = {table: samples found, or None if unreadable}."""
    return [check('C9', sorted(samples), None if found is None else sorted(found), None, table)
            for table, found in sample_sets.items()]


# ----------------------------------------------------------------------------- Kind A

def n_flag(n):
    if n == 0:
        return 'not evaluable'
    if n < LOW_N:
        return 'low n'
    return None


def annotate_computability(metrics, distributions):
    """Kind A rules for one sample: add n, share and flag to the records in place.
    metrics = {id: record}; distributions = [record]."""
    for record in metrics.values():
        denominator_id = record.get('denominator_id')
        if denominator_id:
            denominator = metrics.get(denominator_id, {}).get('value')
            if denominator is None:
                record['flag'] = 'not evaluable'          # denominator missing
                continue
            record['n'] = denominator
            if denominator > 0:
                record['share'] = record['value'] / denominator
        if record.get('n') is not None:
            flag = n_flag(record['n'])
            if flag:
                record['flag'] = flag
    for record in distributions:
        record['n'] = sum(record['y'])
        flag = n_flag(record['n'])
        if flag:
            record['flag'] = flag
