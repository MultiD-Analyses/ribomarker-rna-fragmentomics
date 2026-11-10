"""The unmapped records mapReads.py writes for reads that skip the last cell (release 0.4.0)
must look like the aligners' own unaligned records. Runs Bowtie 1 and 2 on a tiny index; skipped (and so untested)
where the aligners are not installed. Run from the repository root: python -m pytest docker/tests/tools"""
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'docker' / 'scripts'))
import mapReads  # noqa: E402
import mapSchedule  # noqa: E402

TOOLS = ('bowtie', 'bowtie-build', 'bowtie2', 'bowtie2-build')
NAME = 'A00123:8:H5:1:1101:1000:2000 merged_150_20'       # a fastp-merged read name, with its comment
SEQUENCE, QUALITY = 'ACGTTGCAAGCTTACGGATCC', 'IIIIIIIIIIIIIIIIIIII5'   # matches nothing in the reference


@pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS), reason='Bowtie 1 and 2 not installed')
@pytest.mark.parametrize('map_cmd', ['bowtie', 'bowtie2'])
def test_written_unmapped_record_equals_the_aligners_own(tmp_path, map_cmd):
    reference = tmp_path / 'ref.fa'
    reference.write_text('>ref\n' + 'GATTACA' * 30 + '\n')
    index = tmp_path / 'ref'
    subprocess.run([f'{"bowtie" if map_cmd == "bowtie" else "bowtie2"}-build', '-q', str(reference), str(index)],
                   check=True, capture_output=True)
    fastq = tmp_path / 'r.fq'
    fastq.write_text(f'@{NAME}\n{SEQUENCE}\n+\n{QUALITY}\n')
    options = mapSchedule.OPTIONS[map_cmd].format(multimap=100).replace('--no-unal ', '')    # the last cell's form
    row = {'map_cmd': map_cmd, 'aligner': mapSchedule.aligner(map_cmd), 'options': options}
    sam = tmp_path / 'own.sam'
    command = mapSchedule.command(row, index, fastq, tmp_path / 'unal.fq', sam, tmp_path / 'log', 1)
    subprocess.run(command, shell=True, check=True)
    with pysam.AlignmentFile(str(sam), 'r') as f:
        own = next(r for r in f if r.is_unmapped)
        written = mapReads.unmapped_record(f.header, NAME, SEQUENCE, QUALITY, row['aligner'])
    fields = lambda r: (r.query_name, r.flag, r.reference_id, r.reference_start, r.mapping_quality, r.cigarstring,
                        r.query_sequence, list(r.query_qualities), sorted(r.get_tags()))
    assert fields(written) == fields(own)
