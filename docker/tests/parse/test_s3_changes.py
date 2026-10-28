"""Release 0.3.0, step S3: the result changes of parseMappings and their helpers (3b positions, 3c read lengths,
3g UniVec contigs, 3h gene lengths). Run from the repository root: python -m pytest docker/tests/parse"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DOCKER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DOCKER / 'scripts'))
import dbPaths  # noqa: E402


def parse_mappings():
    pytest.importorskip('numba')                         # parseMappings imports numba
    import parseMappings
    return parseMappings


# ----------------------------------------------------------------------------- 3b
@pytest.mark.parametrize('read, strand, suffix', [
    ((100, 130), '+', '_5p'), ((170, 200), '+', '_3p'), ((140, 160), '+', '_internal'), ((105, 195), '+', '_whole'),
    ((100, 130), '-', '_3p'), ((170, 200), '-', '_5p'), ((140, 160), '-', '_internal'), ((105, 195), '-', '_whole'),
    ((110, 130), '+', '_5p'), ((111, 130), '+', '_internal'),       # within 10 nt of the start: start + 10 included
])
def test_position_follows_the_feature_strand(read, strand, suffix):
    """Feature 100-200 (0-based, end exclusive): 5p/3p by the feature's 5' end, on both strands."""
    assert parse_mappings().position_suffix(read[0], read[1], 100, 200, strand) == suffix


# ----------------------------------------------------------------------------- 3h
def test_gene_length_keeps_mirna_arm_names(tmp_path):
    """MirGeneDB arms are named ..._5p in the BED: their length is found; positional and policy suffixes of
    other features are removed only until the name is a BED feature."""
    import pandas as pd
    pm = parse_mappings()
    bed = tmp_path / 'x.bed'
    bed.write_text('Hsa-Mir-21\t10\t32\tHsa-Mir-21_5p\tmiRNA\t+\n'
                   'Homo_sapiens_tRNA-Ala\t10\t82\tHomo_sapiens_tRNA-Ala\ttRNA\t+\n'
                   'ENST1\t10\t160\tRNU2-1\tsnRNA\t+\n')
    features = pd.DataFrame({'feature': ['Hsa-Mir-21_5p', 'Hsa-Mir-21_5p_multimap', 'Homo_sapiens_tRNA-Ala_3p',
                                         'Homo_sapiens_tRNA-Ala_antisense_5p', 'RNU2-1_internal', 'nothing_5p'],
                             'rna_class': ['miRNA', 'miRNA', 'tRNA', 'antisense', 'snRNA', 'x']})
    lengths = pm.add_gene_length(features, str(bed))['length_nt'].tolist()
    assert lengths[:5] == [22, 22, 72, 72, 150]
    assert pd.isna(lengths[5])


# ----------------------------------------------------------------------------- 3g
def test_univec_contigs_from_the_contig_to_tier_map(tmp_path):
    path = tmp_path / 'contig_tiers.tsv'
    path.write_text('contig\ttier\ngnl|uv|X1.1:1-10\tunivec_bt1\nchr1\tgenome\nSR_20_4N-T1\tribospike\n')
    assert parse_mappings().read_univec_contigs(str(path)) == {'gnl|uv|X1.1:1-10'}


def test_contig_to_tier_map_from_the_index_names(tmp_path, monkeypatch):
    """bowtie-inspect -n / bowtie2-inspect -n print one name per line (with a description after a space); a
    stand-in tool prints the names of a fake index."""
    tools = tmp_path / 'bin'
    tools.mkdir()
    for tool in ('bowtie-inspect', 'bowtie2-inspect'):
        script = tools / tool
        script.write_text('#!/bin/sh\ncase "$2" in *univec*) printf "gnl|uv|A1.1:1-5 vector\\ngnl|uv|B2.1:1-9\\n";;'
                          ' *) printf "1 dna:chromosome\\n2\\n";; esac\n')
        script.chmod(0o755)
    monkeypatch.setenv('PATH', f'{tools}{os.pathsep}{os.environ["PATH"]}')
    tiers = [('bowtie', '/r/univec/univec_bt1', 'T'), ('bowtie2', '/r/genome/human-genome', 'T')]
    out = tmp_path / 'contig_tiers.tsv'
    dbPaths.write_contig_tiers(tiers, 'human', str(out))
    assert out.read_text() == ('contig\ttier\ngnl|uv|A1.1:1-5\tunivec_bt1\ngnl|uv|B2.1:1-9\tunivec_bt1\n'
                               '1\tgenome\n2\tgenome\n')
    with pytest.raises(ValueError, match='unknown map_cmd'):
        dbPaths.contig_names('blast', '/r/x')
    # bowtie1, mind and ribospike run with --fullref: SAM holds the whole name, so the map must too (critique
    # 2026-10-27; with the first word only, UniVec mapped with bowtie1 would give no UniVec class)
    assert dbPaths.contig_names('bowtie1', '/r/univec/univec_bt1') == ['gnl|uv|A1.1:1-5 vector', 'gnl|uv|B2.1:1-9']
    assert dbPaths.contig_names('bowtie', '/r/univec/univec_bt1') == ['gnl|uv|A1.1:1-5', 'gnl|uv|B2.1:1-9']


# ----------------------------------------------------------------------------- 3c
def test_read_length_table_has_no_upper_limit(tmp_path):
    """Reads of 151, 301 and 600 nt (before 0.3.0: columns 0-150, and reads over 300 nt stopped the step)."""
    pysam = pytest.importorskip('pysam')
    pm_path = DOCKER / 'scripts' / 'parseMappings.py'
    parse_mappings()
    if shutil.which('samtools') is None:
        pytest.skip('samtools not on PATH (parseMappings sorts and indexes with it)')
    run = tmp_path / 'R'
    (run / 'alignments').mkdir(parents=True)
    (run / 'logs').mkdir()
    (run / 'samples.csv').write_text('S1,g,S1.fq.gz\n')
    contig = 'ENST1'
    bed = tmp_path / 'x.bed'
    bed.write_text(f'{contig}\t10\t990\tGENE1\tlncRNA\t+\n')
    header = {'HD': {'VN': '1.6', 'SO': 'queryname'}, 'SQ': [{'SN': contig, 'LN': 1000}]}
    with pysam.AlignmentFile(str(run / 'alignments' / 'S1.merged.bam'), 'wb', header=header) as bam:
        for number, length in enumerate((151, 301, 600)):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f'r{number}'
            read.query_sequence = 'A' * length
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 20
            read.mapping_quality = 255
            read.cigarstring = f'{length}M'
            read.query_qualities = pysam.qualitystring_to_array('I' * length)
            read.set_tag('AS', 0)
            read.set_tag('XF', f'{contig}|1|Unique')
            bam.write(read)
    result = subprocess.run([sys.executable, str(pm_path), '--samples', 'R/samples.csv', '--threads', '1', '--output',
                             'R', '--bed', str(bed), '--bampath', 'R/alignments', '-q'], cwd=tmp_path,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-2000:]
    rows = (run / 'R-classlengths.csv').read_text().splitlines()
    header_cells = rows[0].split(',')
    assert header_cells[:3] == ['sample', 'rna_class', '0'] and header_cells[-1] == '600'
    counts = dict(zip(header_cells[2:], next(r for r in rows if ',lncRNA,' in r).split(',')[2:]))
    assert counts['151'] == counts['301'] == counts['600'] == '1'


# ----------------------------------------------------------------------------- 3f
def symbol_table(tmp_path, bed, bed_sha256=None):
    import hashlib
    sha = bed_sha256 or hashlib.sha256(bed.read_bytes()).hexdigest()
    path = tmp_path / 'human.gene_symbols.tsv'
    path.write_text(f'# made by a test\n# references_bed_sha256: {sha}\n'
                    'feature\thgnc_id\thgnc_symbol\tmirbase_name\n'
                    'Hsa-Mir-21_5p\tHGNC:31586\tMIR21\thsa-miR-21-5p\n'
                    'Homo_sapiens_tRNA-Ala\tHGNC:1;HGNC:2\tTRA-AGC1-1;TRA-AGC1-2\t\n')
    return path


def test_gene_symbols_join_by_base_feature(tmp_path, capsys):
    import pandas as pd
    pm = parse_mappings()
    bed = tmp_path / 'x.bed'
    bed.write_text('Hsa-Mir-21\t10\t32\tHsa-Mir-21_5p\tmiRNA\t+\n'
                   'Homo_sapiens_tRNA-Ala\t10\t82\tHomo_sapiens_tRNA-Ala\ttRNA\t+\n')
    features = pd.DataFrame({'feature': ['Hsa-Mir-21_5p_multimap', 'Homo_sapiens_tRNA-Ala_3p', 'other'],
                             'rna_class': ['miRNA', 'tRNA', 'x']})
    table = pm.add_gene_symbols(features, str(bed), str(symbol_table(tmp_path, bed)))
    assert list(table.columns) == ['feature', 'rna_class', 'hgnc_id', 'hgnc_symbol', 'mirbase_name']
    assert table['hgnc_symbol'].tolist() == ['MIR21', 'TRA-AGC1-1;TRA-AGC1-2', '']
    assert table['mirbase_name'].tolist() == ['hsa-miR-21-5p', '', '']
    other = pm.add_gene_symbols(features, str(bed), str(symbol_table(tmp_path, bed, bed_sha256='0' * 64)))
    assert other['hgnc_symbol'].tolist() == ['', '', ''] and 'another BED' in capsys.readouterr().err


def test_gene_table_takes_the_symbols_of_its_features(tmp_path):
    """genelengths gets the symbol columns of featureclasses: tRNA fragments are summed into one gene."""
    (tmp_path / 'r-readcounts.csv').write_text('feature,S1\nHomo_sapiens_tRNA-Ala_3p,4\nHomo_sapiens_tRNA-Ala_5p,1\n'
                                               'Hsa-Mir-21_5p,7\n')
    (tmp_path / 'r-featureclasses.csv').write_text(
        'feature,rna_class,length_nt,hgnc_id,hgnc_symbol,mirbase_name\n'
        'Homo_sapiens_tRNA-Ala_3p,tRNA,72.0,HGNC:1,TRA-AGC1-1,\nHomo_sapiens_tRNA-Ala_5p,tRNA,72.0,HGNC:1,TRA-AGC1-1,\n'
        'Hsa-Mir-21_5p,miRNA,22.0,HGNC:31586,MIR21,hsa-miR-21-5p\n')
    result = subprocess.run([sys.executable, str(DOCKER / 'scripts' / 'collapseGenes.py'), '--readcounts',
                             'r-readcounts.csv', '--featureclasses', 'r-featureclasses.csv', '--output', 'r'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'r-genelengths.csv').read_text().splitlines() == [
        'gene,length_nt,rna_class,hgnc_id,hgnc_symbol,mirbase_name',
        'Homo_sapiens_tRNA-Ala,72.0,tRNA,HGNC:1,TRA-AGC1-1,',
        'Hsa-Mir-21_5p,22.0,miRNA,HGNC:31586,MIR21,hsa-miR-21-5p']


def test_merged_mirgenedb_names_get_their_length():
    """MirGeneDB arms with identical sequences are one BED feature 'A_5p|B_5p' (22 of them in level B)."""
    pm = parse_mappings()
    known = {'Hsa-Let-7-P1d_5p|Hsa-Let-7-P2a1_5p', 'Hsa-Mir-21_5p', 'tRNA-X'}
    assert pm.base_feature('Hsa-Let-7-P1d_5p|Hsa-Let-7-P2a1_5p_multimap', known) == 'Hsa-Let-7-P1d_5p|Hsa-Let-7-P2a1_5p'
    assert pm.base_feature('Hsa-Mir-21_5p_antisense', known) == 'Hsa-Mir-21_5p'
    assert pm.base_feature('other|tRNA-X_3p', known) == 'tRNA-X'               # the part after the last '|'


# ----------------------------------------------------------------------------- 3f: join completeness
@pytest.mark.parametrize('species', ['human', 'mouse'])
def test_every_counted_name_leads_back_to_its_bed_feature(species):
    """Plan 3f / §8 (critique 2026-10-27: the planned join test was missing): the names a run writes into
    readcounts/featureclasses (BED feature + _5p/_3p/_internal/_whole for classes with positions, _multimap,
    _antisense, in either order) lead back to their BED feature, through which the annotation columns (length,
    gene symbols) are joined. Every feature of the installed references BED; skipped without it."""
    pm = parse_mappings()
    import refTables
    bed = DOCKER / 'assets' / f'{species}.bed'
    if not bed.exists():
        pytest.skip(f'{bed} not installed')
    labels = {}
    with open(bed) as f:
        for line in f:
            cells = line.rstrip('\n').split('\t')
            labels.setdefault(cells[3], cells[4])
    known = set(labels)
    positional = refTables.ClassTable(refTables.default_table(species), species).position_labels()
    wrong = []
    for feature, label in labels.items():
        names = [feature, feature + '_multimap', feature + '_antisense']
        if label in positional:
            for position in ('_5p', '_3p', '_internal', '_whole'):
                names += [feature + position, feature + position + '_multimap', feature + '_antisense' + position]
        wrong += [(name, pm.base_feature(name, known)) for name in names if pm.base_feature(name, known) != feature]
    assert not wrong, f'{len(wrong)} names do not lead back to their feature, e.g. {wrong[:5]}'
