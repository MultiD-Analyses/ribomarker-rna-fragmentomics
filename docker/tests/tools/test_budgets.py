"""The pass settings with the aligners themselves (release 0.4.0): which reads each value of
max_mismatches and n_overhang lets align, with the exact options of mapSchedule.cell_options. Written after the
review of 2026-11-05: before, only the option strings were tested.

Bowtie 2 (--ignore-quals): a mismatch costs 6, a position on an N 5 (--np 5; 1 with n_overhang), a gap 5 + 3 per base;
--score-min C,-6m: m = 0 exact, 1 one mismatch or one N, 2 two mismatches, a mismatch and an N, or a gap of 1-2 bases.
Bowtie 1 (-v m): at most m mismatches end to end. Reads of 50 nt whose first 20 nt are unchanged, so that the seed
search finds them: a budget is a limit, not a promise that the heuristics find every alignment within it (a 30-nt read
with a gap in its middle is not found even with the standard settings; test_a_gap_in_a_short_read_is_not_found).
Skipped (and so untested) without Bowtie 1 and 2; run them where the pipeline runs (WSL, conda environment).
Run from the repository root: python -m pytest docker/tests/tools"""
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'docker' / 'scripts'))
import mapSchedule  # noqa: E402

TOOLS = ('bowtie', 'bowtie-build', 'bowtie2', 'bowtie2-build')
MISSING = [t for t in TOOLS if shutil.which(t) is None]     # the tests that run an aligner are skipped without them

random.seed(4)
CONTIG = ''.join(random.choice('ACGT') for _ in range(400))      # one record; its end padded with N (below)
PAD = 'N' * 10
START, LENGTH = 100, 50                                          # reads come from CONTIG[100:150]


def change(sequence, position):
    """Another base at one position (a mismatch)."""
    return sequence[:position] + {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[sequence[position]] + sequence[position + 1:]


def read(mismatches=(), deleted=0, n=(), length=LENGTH, gap_at=32):
    """length nt of CONTIG from START: mismatches at the given read positions, `deleted` bases of the reference
    missing after read position gap_at (a gap in the read), N at the given positions."""
    sequence = CONTIG[START:START + gap_at] + CONTIG[START + gap_at + deleted:START + length + deleted]
    for position in mismatches:
        sequence = change(sequence, position)
    for position in n:
        sequence = sequence[:position] + 'N' + sequence[position + 1:]
    return sequence


READS = {
    'exact': read(),
    'mm1': read(mismatches=(25,)),
    'mm2': read(mismatches=(25, 45)),
    'mm3': read(mismatches=(25, 38, 45)),
    'gap1': read(deleted=1),
    'gap2': read(deleted=2),
    'mm1_gap1': read(mismatches=(25,), deleted=1),
    'n1': read(n=(40,)),
    'n2': read(n=(40, 41)),
    'mm1_n1': read(mismatches=(25,), n=(40,)),
    'overhang1': CONTIG[-(LENGTH - 1):] + 'A',     # the last base lies on the record's N padding
    'short_gap1': read(deleted=1, length=30, gap_at=15),
}


@pytest.fixture(scope='module')
def indexes(tmp_path_factory):
    if MISSING:
        pytest.skip(f'not installed: {", ".join(MISSING)}')
    folder = tmp_path_factory.mktemp('budget')
    (folder / 'ref.fa').write_text(f'>ref\n{CONTIG}{PAD}\n')
    fastq = folder / 'reads.fq'
    fastq.write_text(''.join(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n' for name, seq in READS.items()))
    for tool in ('bowtie-build', 'bowtie2-build'):
        subprocess.run([tool, '-q', str(folder / 'ref.fa'), str(folder / tool.split('-')[0])], check=True,
                       capture_output=True)
    return folder, fastq


def aligned(indexes, map_cmd, max_mismatches, n_overhang='standard'):
    """{read name: score} of the reads aligned with the cell's options (AS:i, or -NM:i for Bowtie 1)."""
    folder, fastq = indexes
    options = mapSchedule.cell_options(map_cmd, max_mismatches, n_overhang, 100)
    row = {'map_cmd': map_cmd, 'aligner': mapSchedule.aligner(map_cmd), 'options': ' '.join(options)}
    label = f'{map_cmd}_{max_mismatches}_{n_overhang}'
    sam = folder / f'{label}.sam'
    index = folder / ('bowtie' if row['aligner'] == 'bowtie' else 'bowtie2')
    command = mapSchedule.command(row, index, fastq, folder / f'{label}.unal.fq', sam, folder / f'{label}.log', 1)
    subprocess.run(command, shell=True, check=True)
    with pysam.AlignmentFile(str(sam), 'r') as f:
        return {r.query_name: (r.get_tag('AS') if r.has_tag('AS') else -r.get_tag('NM'))
                for r in f if not r.is_unmapped and not r.is_secondary}


# ----------------------------------------------------------------------------- Bowtie 2
@pytest.mark.parametrize('max_mismatches, n_overhang, expected', [
    ('0', 'standard', {'exact': 0}),
    ('1', 'standard', {'exact': 0, 'mm1': -6, 'n1': -5, 'overhang1': -5}),      # an N of the padding is an N
    ('2', 'standard', {'exact': 0, 'mm1': -6, 'mm2': -12, 'gap1': -8, 'gap2': -11, 'n1': -5, 'n2': -10, 'mm1_n1': -11,
                       'overhang1': -5}),
    ('0', '1', {'exact': 0, 'n1': -1, 'overhang1': -1}),
    ('0', '2', {'exact': 0, 'n1': -1, 'n2': -2, 'overhang1': -1}),
])
def test_bowtie2_budget(indexes, max_mismatches, n_overhang, expected):
    found = aligned(indexes, 'bowtie2', max_mismatches, n_overhang)
    assert found == expected


def test_a_gap_never_takes_the_place_of_a_mismatch(indexes):
    """max_mismatches 1 (C,-6): one gap costs at least 8, so no gapped read aligns; 2 (C,-12) admits a 1- or 2-nt
    gap but not a gap plus a mismatch (14)."""
    one, two = aligned(indexes, 'bowtie2', '1'), aligned(indexes, 'bowtie2', '2')
    assert not {'gap1', 'gap2', 'mm1_gap1'} & set(one)
    assert {'gap1', 'gap2'} <= set(two) and 'mm1_gap1' not in two


def test_bowtie2_standard_settings_are_unchanged(indexes):
    """standard: the pipeline's budget of 0.6 + 0.6 x 50 = 30.6 (three mismatches, a gap and a mismatch, ...)."""
    found = aligned(indexes, 'bowtie2', 'standard')
    assert set(READS) - {'short_gap1'} == set(found)


def test_a_gap_in_a_short_read_is_not_found(indexes):
    """A 30-nt read with a 1-nt gap after position 15 is within every budget from 2 on (8 points), but each 20-nt
    seed of Bowtie 2 --very-sensitive overlaps the gap, so the read is not found even with the standard settings:
    the documentation must not promise that `2` aligns every read with a short gap (finding of these tests)."""
    for max_mismatches in ('2', 'standard'):
        assert 'short_gap1' not in aligned(indexes, 'bowtie2', max_mismatches)


# ----------------------------------------------------------------------------- Bowtie 1
@pytest.mark.parametrize('map_cmd', ['bowtie', 'bowtie1'])
@pytest.mark.parametrize('max_mismatches, expected', [
    ('0', {'exact'}),
    ('1', {'exact', 'mm1', 'n1'}),                 # an N in the read is a mismatch
    ('2', {'exact', 'mm1', 'mm2', 'n1', 'n2', 'mm1_n1'}),
])
def test_bowtie1_v(indexes, map_cmd, max_mismatches, expected):
    """-v m: at most m mismatches end to end, no gaps; an alignment touching an N of the reference is invalid."""
    found = aligned(indexes, map_cmd, max_mismatches)
    assert set(found) == expected
    assert all(score == -sum(a != b for a, b in zip(READS[name], CONTIG[START:START + LENGTH]))
               for name, score in found.items())


@pytest.mark.parametrize('map_cmd', ['bowtie', 'bowtie1', 'mind', 'ribospike'])
def test_every_bowtie1_form_runs_on_gzipped_input(indexes, tmp_path, map_cmd):
    """Each Bowtie 1 form runs with the pinned Bowtie 1 (1.0.0 has no --no-unal; found by these tests 2026-11-05) on
    the .fq.gz input a Bowtie 1 cell gets after a Bowtie 2 cell, and aligns the exact read."""
    import gzip
    folder, fastq = indexes
    gz = tmp_path / 'reads.fq.gz'
    with gzip.open(gz, 'wt') as f:
        f.write(fastq.read_text())
    row = {'map_cmd': map_cmd, 'aligner': 'bowtie', 'options': mapSchedule.OPTIONS[map_cmd].format(multimap=100)}
    sam = tmp_path / 'out.sam'
    command = mapSchedule.command(row, folder / 'bowtie', gz, tmp_path / 'unal.fq', sam, tmp_path / 'log', 1)
    result = subprocess.run(command, shell=True)
    assert result.returncode == 0, (tmp_path / 'log').read_text()[-500:]
    with pysam.AlignmentFile(str(sam), 'r') as f:
        names = {r.query_name for r in f if not r.is_unmapped}
    assert 'exact' in names


def test_no_bowtie1_form_passes_no_unal():
    """Bowtie 1.0.0 (env.lock.txt) stops on --no-unal; checked here also where Bowtie is not installed."""
    for map_cmd in ('bowtie', 'bowtie1', 'mind', 'ribospike'):
        assert '--no-unal' not in mapSchedule.OPTIONS[map_cmd], map_cmd
