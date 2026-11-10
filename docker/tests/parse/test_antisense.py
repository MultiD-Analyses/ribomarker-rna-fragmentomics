"""Release 0.4.0: antisense reads in parseMappings. Every antisense
read keeps its class and gets _antisense (also protein_coding); mappingstats counts Antisense and AntisenseMulti
whatever the flags. Runs parseMappings on a small merged BAM; needs samtools.
Run from the repository root: python -m pytest docker/tests/parse"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DOCKER = Path(__file__).resolve().parents[2]
CONTIGS = {'ENST1': 1000, 'LNC': 1000, 'chrX': 1000}
BED = 'ENST1\t10\t990\tGENE1\tprotein_coding\t+\nLNC\t10\t990\tLNC1\tlncRNA\t+\n'      # chrX: no feature
READS = {                                   # name: [(contig, reverse)], all alignments with AS 0 (equal scores)
    'a_sense': [('LNC', False)],                            # Unique, sense
    'b_lnc_anti': [('LNC', True)],                          # Unique, antisense: lncRNA, LNC1_antisense
    'c_pc_anti': [('ENST1', True)],                         # Unique, antisense: protein_coding, GENE1_antisense
    'd_tie': [('ENST1', True), ('LNC', False)],             # Multi: sense before antisense -> LNC1, sense
    'e_anti': [('ENST1', True), ('LNC', True)],             # Multi: both antisense, lncRNA ranks higher
    'f_none': [('chrX', False), ('LNC', True)],             # Multi: antisense before no feature -> LNC1, antisense
}


def run_parse(tmp_path, flags):
    pysam = pytest.importorskip('pysam')
    pytest.importorskip('numba')                            # parseMappings imports numba
    if shutil.which('samtools') is None:
        pytest.skip('samtools not on PATH (parseMappings sorts and indexes with it)')
    run = tmp_path / 'R'
    (run / 'alignments').mkdir(parents=True)
    (run / 'logs').mkdir()
    (run / 'samples.csv').write_text('S1,g,S1.fq.gz\n')
    (tmp_path / 'x.bed').write_text(BED)
    header = {'HD': {'VN': '1.6', 'SO': 'queryname'}, 'SQ': [{'SN': c, 'LN': n} for c, n in CONTIGS.items()]}
    with pysam.AlignmentFile(str(run / 'alignments' / 'S1.merged.bam'), 'wb', header=header) as bam:
        for name, alignments in READS.items():
            for number, (contig, reverse) in enumerate(alignments):
                read = pysam.AlignedSegment(bam.header)
                read.query_name, read.query_sequence = name, 'ACGT' * 6
                read.flag = (16 if reverse else 0) | (256 if number else 0)
                read.reference_id, read.reference_start = list(CONTIGS).index(contig), 100
                read.mapping_quality, read.cigarstring = 255, '24M'
                read.query_qualities = pysam.qualitystring_to_array('I' * 24)
                read.set_tag('AS', 0)
                bam.write(read)
    command = [sys.executable, str(DOCKER / 'scripts' / 'parseMappings.py'), '--samples', 'R/samples.csv',
               '--threads', '1', '--output', 'R', '--bed', str(tmp_path / 'x.bed'), '--bampath', 'R/alignments']
    command += ['--remove-multimap'] * flags[0] + ['--remove-antisense'] * flags[1]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-2000:]
    table = lambda name: [line.split(',') for line in (run / f'R-{name}.csv').read_text().splitlines()]
    stats = table('mappingstats')
    classes = {row[0]: int(row[1]) for row in table('classcounts')[1:]}
    features = {row[0]: int(row[1]) for row in table('readcounts')[1:]}
    return dict(zip(stats[0], stats[1])), classes, features


FLAGS = [(False, False), (False, True), (True, False), (True, True)]


@pytest.mark.parametrize('flags', FLAGS)
def test_mappingstats_count_antisense_reads_in_every_flag_state(tmp_path, flags):
    stats, _, _ = run_parse(tmp_path, flags)
    assert list(stats) == ['sample', 'Mapped', 'Unmapped', 'Unique', 'Multi', 'Antisense', 'AntisenseMulti', 'Total']
    assert {k: stats[k] for k in ('Mapped', 'Unique', 'Multi', 'Antisense', 'AntisenseMulti')} == \
           {'Mapped': '6', 'Unique': '3', 'Multi': '3', 'Antisense': '2', 'AntisenseMulti': '2'}


def test_kept_antisense_reads_keep_their_class_and_get_the_suffix(tmp_path):
    _, classes, features = run_parse(tmp_path, (False, False))
    assert {k: v for k, v in classes.items() if v} == {'lncRNA': 5, 'protein_coding': 1}     # no class antisense
    assert {k: v for k, v in features.items() if v} == {
        'LNC1': 1, 'LNC1_antisense': 1, 'GENE1_antisense': 1,              # protein_coding tagged (the gap closed)
        'LNC1_multimap': 1, 'LNC1_antisense_multimap': 2}


@pytest.mark.parametrize('flags, kept', [
    ((False, True), {'lncRNA': 2}),                         # mixed state: antisense multimappers removed as antisense
    ((True, False), {'lncRNA': 2, 'protein_coding': 1}),
    ((True, True), {'lncRNA': 1}),
])
def test_removal_by_the_flags(tmp_path, flags, kept):
    _, classes, _ = run_parse(tmp_path, flags)
    assert {k: v for k, v in classes.items() if v} == kept
