"""The ncRNA split (buildDB/split_ncrna.py; release 0.4.0): the split itself and its build checks.
Run from the repository root: python -m pytest docker/tests/references"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'buildDB'))
import split_ncrna  # noqa: E402

RECORDS = [('T1|LNCA', 'lncRNA', 'ACGT\nACGT\n'), ('T2|SNORD', 'snoRNA', 'GGGG\n'), ('T3|MTT', 'Mt_tRNA', 'TTTT\n'),
           ('T4|U6', 'snRNA', 'CCCC\nCC\n'), ('T5|MISC', 'misc_RNA', 'AAAA\n'), ('T6|LNCB', 'lncRNA', 'ACAC\n'),
           ('T7|RNA5S', 'rRNA', 'GAGA\n'), ('T8|SCA', 'scaRNA', 'TCTC\n'), ('T9|MTR', 'Mt_rRNA', 'AGAG\n')]


def write_inputs(tmp_path, records=RECORDS, bed_records=None):
    fasta, bed = tmp_path / 'x.ncRNA.fa', tmp_path / 'x.ncRNA.bed'
    fasta.write_text(''.join(f'>{name}\n{sequence}' for name, _, sequence in records))
    bed.write_text(''.join(f'{name}\t10\t20\t{name.split("|")[1]}\t{biotype}\t+\n'
                           for name, biotype, _ in (bed_records or records)))
    return fasta, bed


def out(tmp_path):
    return lambda sub: str(tmp_path / f'x.{sub}.fa')


def test_every_record_goes_unchanged_to_its_sub_reference_in_order(tmp_path):
    fasta, bed = write_inputs(tmp_path)
    counts = split_ncrna.split_fasta(str(fasta), str(bed), out(tmp_path))
    assert counts == {'ncrna_mt': {'Mt_tRNA': 1, 'Mt_rRNA': 1}, 'ncrna_snsno': {'snoRNA': 1, 'snRNA': 1, 'scaRNA': 1},
                      'ncrna_other': {'misc_RNA': 1, 'rRNA': 1}, 'ncrna_lnc': {'lncRNA': 2}}
    assert (tmp_path / 'x.ncrna_lnc.fa').read_text() == '>T1|LNCA\nACGT\nACGT\n>T6|LNCB\nACAC\n'
    assert (tmp_path / 'x.ncrna_snsno.fa').read_text() == '>T2|SNORD\nGGGG\n>T4|U6\nCCCC\nCC\n>T8|SCA\nTCTC\n'
    assert (tmp_path / 'x.ncrna_mt.fa').read_text() == '>T3|MTT\nTTTT\n>T9|MTR\nAGAG\n'


def test_the_order_of_the_sub_references_is_the_decided_one():
    assert [name for name, _ in split_ncrna.SUB_REFERENCES] == ['ncrna_mt', 'ncrna_snsno', 'ncrna_other', 'ncrna_lnc']


def test_a_biotype_outside_every_sub_reference_stops(tmp_path):
    fasta, bed = write_inputs(tmp_path, RECORDS + [('T10|VT', 'vault_RNA', 'ACGT\n')])
    with pytest.raises(SystemExit, match='vault_RNA'):
        split_ncrna.split_fasta(str(fasta), str(bed), out(tmp_path))


def test_records_without_a_bed_row_stop(tmp_path):
    fasta, bed = write_inputs(tmp_path, bed_records=RECORDS[:-1])
    with pytest.raises(SystemExit, match='differ'):
        split_ncrna.split_fasta(str(fasta), str(bed), out(tmp_path))


def test_the_build_check_finds_a_changed_or_lost_record(tmp_path):
    fasta, bed = write_inputs(tmp_path)
    counts = split_ncrna.split_fasta(str(fasta), str(bed), out(tmp_path))
    records = split_ncrna.read_records(str(fasta))
    biotype = split_ncrna.read_biotypes(str(bed))
    (tmp_path / 'x.ncrna_other.fa').write_text('>T5|MISC\nAAAT\n>T7|RNA5S\nGAGA\n')       # one base changed
    with pytest.raises(SystemExit, match='changed'):
        split_ncrna.check_split(records, biotype, out(tmp_path), counts)
    (tmp_path / 'x.ncrna_other.fa').write_text('>T5|MISC\nAAAA\n')                         # one record lost
    with pytest.raises(SystemExit, match='exactly once'):
        split_ncrna.check_split(records, biotype, out(tmp_path), counts)
    (tmp_path / 'x.ncrna_other.fa').write_text('>T7|RNA5S\nGAGA\n>T5|MISC\nAAAA\n')       # order changed
    with pytest.raises(SystemExit, match='order'):
        split_ncrna.check_split(records, biotype, out(tmp_path), counts)


def test_sequence_lengths_ignore_line_breaks():
    assert split_ncrna.sequence_length('>a\nACGT\nAC\n') == 6
