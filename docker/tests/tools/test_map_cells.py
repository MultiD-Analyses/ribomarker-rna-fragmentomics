"""mapReads.py's loop over the cells of a schedule with passes and windows (release 0.4.0),
with a fake aligner (the commands are parsed, no Bowtie is run). The comparison with the aligners' own unaligned
records is in test_unmapped_records.py (needs Bowtie 1 and 2).
Run from the repository root: python -m pytest docker/tests/tools"""
import gzip
import io
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'docker' / 'scripts'))
import dbPaths  # noqa: E402
import mapReads  # noqa: E402
import mapSchedule  # noqa: E402

DBPATHS = 'bowtie,univec/univec_bt1,T\nbowtie2,gtf/human-cdna,T\nbowtie2,genome/human-genome,T\n'
PASSES = ('pass\tapplies_to\tmax_mismatches\tmin_length\tmax_length\tn_overhang\n'
          '1\tall\t0\t-\t-\t-\n'
          '2\tcdna,genome\tstandard\t20\t-\t-\n')
READS = {                                        # name (with a comment, as fastp writes it): sequence
    'r_exact merged_1': 'A' * 18,               # aligns in pass 1, cdna
    'r_pass2 merged_2': 'C' * 25,               # aligns in pass 2, genome
    'r_short merged_3': 'G' * 17,               # unaligned in pass 1, under 20 nt: skips pass 2, unmapped (pysam)
    'r_long merged_4': 'T' * 30,                # aligns nowhere: unmapped (the genome cell's own record)
}
ALIGNS = {('1', 'human-cdna'): {'A' * 18}, ('2', 'human-genome'): {'C' * 25}}


def fake_aligner(cmd, calls, aligns):
    """Writes the SAM, the unaligned reads and the log of one cell from its command."""
    calls.append(cmd)
    if cmd.startswith('rm -f'):
        return
    fastq = re.search(r'(?:zcat -f|-U) (\S+)', cmd).group(1)
    sam = re.search(r'(\S+\.sam)', cmd).group(1)
    unal = re.search(r'--un(?:-gz)? (\S+)', cmd).group(1)
    pass_number, index_name = re.search(r'\.p(\d+)\.([^/]+)\.sam$', sam).groups()
    aligned = aligns.get((pass_number, index_name), set())
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'SN': index_name, 'LN': 1000}]}
    with pysam.AlignmentFile(sam, 'w', header=header) as out, mapReads.open_fastq(unal, 'wt') as left:
        for name, sequence, quality in mapReads.fastq_records(fastq):
            record = pysam.AlignedSegment(out.header)
            record.query_name, record.query_sequence = name, sequence
            record.query_qualities = pysam.qualitystring_to_array(quality)
            if sequence in aligned:
                record.flag, record.reference_id, record.reference_start = 0, 0, 10
                record.cigarstring = f'{len(sequence)}M'
                out.write(record)
            else:
                left.write(f'@{name}\n{sequence}\n+\n{quality}\n')
                if '--no-unal' not in cmd:           # the last cell writes its unaligned reads too
                    record.flag = 4
                    record.set_tag('YT', 'UU')
                    out.write(record)
    Path(re.search(r'2> (\S+)', cmd).group(1)).write_text('log\n')


def fake_merge(self, sample, bams):
    """samtools merge + sort -n, replaced by pysam (no samtools here)."""
    records, sq = [], {}
    for bam in bams:
        with pysam.AlignmentFile(bam, 'rb', check_sq=False) as f:
            for s in f.header.to_dict().get('SQ', []):
                sq[s['SN']] = s['LN']
            records += [r.to_string() for r in f]
    header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.0'}, 'SQ': [{'SN': n, 'LN': n_} for n, n_ in sq.items()]})
    path = f'{self.args.output}/alignments/{sample}.merged.bam'
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for line in sorted(records):
            out.write(pysam.AlignedSegment.fromstring(line, header))


def run_cells(tmp_path, monkeypatch, aligns):
    """Maps sample S1 (READS) through the schedule of DBPATHS and PASSES; returns (output folder, merged BAM records,
    commands run)."""
    (tmp_path / 'human.dbpaths.csv').write_text(DBPATHS)
    (tmp_path / 'passes.tsv').write_text(PASSES)
    rows, _, problems, _ = mapSchedule.resolve(str(tmp_path / 'human.dbpaths.csv'), 'human', 100, None, False,
                                               str(tmp_path / 'passes.tsv'))
    assert problems == []
    mapSchedule.write_schedule(rows, tmp_path / 'schedule.tsv')
    with gzip.open(tmp_path / 'S1.fq.gz', 'wt') as f:
        f.write(''.join(f'@{n}\n{s}\n+\n{"I" * len(s)}\n' for n, s in READS.items()))
    (tmp_path / 'samples.csv').write_text(f'S1,g,{tmp_path / "S1.fq.gz"}\n')
    out = tmp_path / 'out'
    for folder in ('logs', 'alignments/logs', 'alignments/unal'):
        (out / folder).mkdir(parents=True)
    calls = []
    monkeypatch.setattr(mapReads, 'run', lambda cmd: fake_aligner(cmd, calls, aligns))
    monkeypatch.setattr(mapReads.os, 'popen', lambda cmd: io.StringIO('tool version 1.0.0\n'))
    monkeypatch.setattr(mapReads.MapReads, 'merge_bams', fake_merge)
    args = SimpleNamespace(samples=str(tmp_path / 'samples.csv'), dbs=str(tmp_path / 'human.dbpaths.csv'),
                           schedule=str(tmp_path / 'schedule.tsv'), threads=1, output=str(out), xp=True)
    mapReads.MapReads(args).map_reads('S1')
    with pysam.AlignmentFile(str(out / 'alignments' / 'S1.merged.bam'), 'rb', check_sq=False) as bam:
        records = list(bam)
    return out, records, calls


@pytest.fixture
def mapped(tmp_path, monkeypatch):
    return run_cells(tmp_path, monkeypatch, ALIGNS)


def test_every_read_is_in_the_merged_bam_once(mapped):
    _, records, _ = mapped
    assert sorted(r.query_name for r in records) == sorted(READS)


def test_mapped_records_carry_their_pass_and_unmapped_ones_none(mapped):
    _, records, _ = mapped
    found = {r.query_name: (r.is_unmapped, r.get_tag('XP') if r.has_tag('XP') else None) for r in records}
    assert found == {'r_exact merged_1': (False, 1), 'r_pass2 merged_2': (False, 2),
                     'r_short merged_3': (True, None), 'r_long merged_4': (True, None)}


def test_a_read_skipped_by_the_last_cells_window_gets_the_bowtie2_unaligned_form(mapped):
    _, records, _ = mapped
    short = next(r for r in records if r.query_name == 'r_short merged_3')
    assert (short.flag, short.reference_name, short.get_tag('YT'), short.query_sequence) == (4, None, 'UU', 'G' * 17)


def test_the_cells_table_records_runs_and_window_skips(mapped):
    out, _, calls = mapped
    rows = (out / 'alignments' / 'logs' / 'S1.cells.tsv').read_text().splitlines()
    assert rows == ['pass\tindex_name\tran\tskipped_window',
                    '1\tunivec_bt1\t1\t0', '1\thuman-cdna\t1\t0', '1\thuman-genome\t1\t0',
                    '2\thuman-cdna\t1\t1', '2\thuman-genome\t1\t1']
    commands = [c for c in calls if not c.startswith('rm -f')]
    assert [re.search(r'\.p(\d)\.(\S+)\.sam', c).groups() for c in commands] == [
        ('1', 'univec_bt1'), ('1', 'human-cdna'), ('1', 'human-genome'), ('2', 'human-cdna'), ('2', 'human-genome')]
    assert '--score-min C,0' in commands[1] and '-v 0' in commands[0] and '--no-unal' not in commands[-1]


def test_a_cell_without_reads_is_not_run(tmp_path, monkeypatch):
    out, records, calls = run_cells(tmp_path, monkeypatch, {('1', 'univec_bt1'): set(READS.values())})
    rows = (out / 'alignments' / 'logs' / 'S1.cells.tsv').read_text().splitlines()[1:]
    assert [r.split('\t')[2] for r in rows] == ['1', '0', '0', '0', '0']       # only UniVec ran
    assert len([c for c in calls if not c.startswith('rm -f')]) == 1
    assert sorted(r.query_name for r in records if not r.is_unmapped) == sorted(READS)


def test_fastq_split_by_length_keeps_order_and_counts(tmp_path):
    source = tmp_path / 'in.fq'
    source.write_text(''.join(f'@r{n} x\n{"A" * n}\n+\n{"I" * n}\n' for n in (15, 20, 25, 30)))
    counts = mapReads.split_by_length(source, tmp_path / 'a.fq.gz', tmp_path / 'b.fq', '20', '25')
    assert counts == (2, 2)
    assert [r[0] for r in mapReads.fastq_records(tmp_path / 'a.fq.gz')] == ['r20 x', 'r25 x']
    assert [r[0] for r in mapReads.fastq_records(tmp_path / 'b.fq')] == ['r15 x', 'r30 x']
    mapReads.append_fastq(tmp_path / 'b.fq', tmp_path / 'a.fq.gz')               # a second gzip member
    assert [r[0] for r in mapReads.fastq_records(tmp_path / 'a.fq.gz')] == ['r20 x', 'r25 x', 'r15 x', 'r30 x']


def test_map_log_lists_the_references_in_the_order_in_effect(tmp_path, monkeypatch):
    """Pass 1 has the genome only: map.log's databases (read by qcCollect and compare_reads) follow tables/order.tsv,
    not the order of the first cells (review of 2026-11-05)."""
    (tmp_path / 'human.dbpaths.csv').write_text(DBPATHS)
    (tmp_path / 'passes.tsv').write_text(PASSES.split('\n')[0] + '\n1\tgenome\t0\t-\t-\t-\n2\tall\tstandard\t-\t-\t-\n')
    rows, _, problems, _ = mapSchedule.resolve(str(tmp_path / 'human.dbpaths.csv'), 'human', 100, None, False,
                                               str(tmp_path / 'passes.tsv'))
    assert problems == [] and [r['reference'] for r in rows] == ['genome', 'univec_bt1', 'cdna', 'genome']
    (tmp_path / 'samples.csv').write_text(f'S1,g,{tmp_path / "S1.fq.gz"}\n')
    monkeypatch.setattr(mapReads.os, 'popen', lambda cmd: io.StringIO('tool version 1.0.0\n'))    # no aligner run
    args = SimpleNamespace(samples=str(tmp_path / 'samples.csv'), dbs=str(tmp_path / 'human.dbpaths.csv'),
                           schedule=str(tmp_path / 'schedule.tsv'), threads=1, output=str(tmp_path / 'out'), xp=True)
    mapSchedule.write_schedule(rows, tmp_path / 'schedule.tsv')
    first_cells = [Path(p).name for p in mapReads.MapReads(args).dbs['db_path']]
    assert first_cells == ['human-genome', 'univec_bt1', 'human-cdna']               # runs without order.tsv
    names = mapSchedule.order_in_effect(str(tmp_path / 'human.dbpaths.csv'), 'human')
    mapSchedule.write_order(rows, tmp_path / 'order.tsv', names)                      # as rule run_tables writes it
    assert [Path(p).name for p in mapReads.MapReads(args).dbs['db_path']] == ['univec_bt1', 'human-cdna', 'human-genome']
