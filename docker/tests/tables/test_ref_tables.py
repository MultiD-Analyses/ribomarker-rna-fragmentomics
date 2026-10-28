"""The RNA class tables (release 0.3.0): they reproduce the constants they replace (S1), their class ids are
those of the translation map (S2), and broken tables are refused."""
import os
import sys
from pathlib import Path

import pytest

DOCKER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DOCKER / 'scripts'))
import refTables  # noqa: E402

# The constants of the pipeline before 0.3.0, copied here as the reference (parseMappings.py genetype_hierarchy,
# collapseGenes.py DROP_BIOTYPES / EXCLUDE_BIOTYPES / FRAGMENT_BIOTYPES, qcCollect.py and qcReport.py
# RRNA_BIOTYPES, qcChecks.py NOT_GENES, qcReport.py TYPE_GROUPS with its default).
OLD_HIERARCHY = {
    'unannotated': 0, 'protein_coding': 1, 'pre-mRNA': 1, 'sinRNA': 1, 'scdRNA': 1, 'sutRNA': 1, 'piRNA': 2,
    'lncRNA': 2, 'SINE': 3, 'LINE': 3, 'pre-tRNA': 4, 'IG_V_pseudogene': 4, 'IG_V_gene': 5, 'IG_C_gene': 5,
    'IG_J_gene': 5, 'TR_V_gene': 5, 'TR_V_pseudogene': 4, 'TR_D_gene': 5, 'TR_J_gene': 5, 'TR_C_gene': 5,
    'IG_C_pseudogene': 4, 'TR_J_pseudogene': 4, 'IG_J_pseudogene': 4, 'IG_D_gene': 5, 'Y_RNA': 5,
    'vault_RNA': 5, 'scRNA': 5, 'scaRNA': 5, 'snoRNA': 5, 'snRNA': 5, 'pre-miRNA': 4, 'star-miRNA': 4,
    'Mt_rRNA': 6, 'miRNA': 6, 'rRNA': 7, 'Mt_tRNA': 8, 'tRNA': 9}
OLD_DROP = {'pre-mRNA', 'pre-tRNA', 'miRNA-pre'}
OLD_FRAGMENT = {'tRNA', 'Mt_tRNA'}
OLD_RRNA = {'rRNA', 'Mt_rRNA'}
OLD_NOT_GENES = {'pre-mRNA', 'pre-tRNA', 'miRNA-pre', 'unannotated'}
OLD_TYPE_GROUPS = {'rRNA': 'rRNA', 'tRNA': 'tRNA', 'miRNA': 'miRNA', 'miRNA-star': 'miRNA', 'miRNA-pre': 'miRNA',
                   'piRNA': 'piRNA', 'protein_coding': 'mRNA',
                   'Mt_rRNA': 'lncRNA and other ncRNA', 'Mt_tRNA': 'lncRNA and other ncRNA',
                   'pre-mRNA': 'genome, precursors, introns', 'pre-tRNA': 'genome, precursors, introns',
                   'unannotated': 'genome, precursors, introns'}
BED_LABELS = ['rRNA', 'Mt_rRNA', 'tRNA', 'Mt_tRNA', 'pre-tRNA', 'miRNA', 'miRNA-star', 'miRNA-pre', 'piRNA',
              'snRNA', 'snoRNA', 'scaRNA', 'misc_RNA', 'lncRNA', 'protein_coding', 'pre-mRNA']   # human = mouse


def old_rank(label):
    return OLD_HIERARCHY.get(label, 10)


# S4 (report, RealSeq review #8, #9): mRNA is shown as protein_coding, and the genome group is split in two
S4_GROUPS = {'mRNA': 'protein_coding', 'genome, precursors, introns': 'genome: pre-mRNA, pre-tRNA'}


def old_group(label):
    if label.startswith('RiboSpike') or label == 'miND':   # miND: spike-ins since the critique of 2026-10-27
        return 'spike-ins'
    if label == 'unannotated':
        return 'genome (unannotated)'
    group = OLD_TYPE_GROUPS.get(label, 'lncRNA and other ncRNA')
    return S4_GROUPS.get(group, group)


@pytest.fixture(params=['human', 'mouse', 'default'])
def table(request):
    path = DOCKER / 'tables' / f'{request.param}.rna_classes.tsv'
    return refTables.ClassTable(str(path), request.param)


def test_ranks_reproduce_the_old_hierarchy(table):
    # every label of the table, system labels included: a spike or miND label in a BED had rank 10 before 0.3.0
    ranks = table.rank_map(BED_LABELS + ['unannotated'])
    for label in set(BED_LABELS) | {'unannotated'} | set(OLD_HIERARCHY) | set(ranks):
        assert ranks[label] == old_rank(label), label


def test_gene_table_reproduces_collapse_genes(table):
    """collapseGenes sees class ids since S2: the rule of a label's class is the old rule of the label."""
    for label in BED_LABELS + ['unannotated', 'unmapped', 'RiboSpike-T1', 'RiboSpike-Short', 'miND']:
        old = 'drop' if label in OLD_DROP else (
            'exclude' if label in ('unannotated', 'unmapped') or label.startswith('RiboSpike') else 'yes')
        assert table.gene_table(table.class_of(label)) == old, label
    # 3b: snRNA and snoRNA reads get positions too
    assert table.position_labels() == OLD_FRAGMENT | {'snRNA', 'snoRNA'}
    assert table.position_classes() == OLD_FRAGMENT | {'snRNA', 'snoRNA'}


def test_qc_sets_reproduce_the_old_constants(table):
    assert table.rrna_classes() == OLD_RRNA
    not_genes = {c for c in table.not_gene_classes() if not c.startswith('RiboSpike')}
    # unmapped and Natural are not assign.class metrics, so C7 never sees them; UniVec: 3g
    assert not_genes - {'unmapped', 'Natural'} == {table.class_of(label) for label in OLD_NOT_GENES} | {'UniVec'}
    for label in BED_LABELS + ['unannotated', 'RiboSpike-T1', 'miND', 'antisense', 'something_new']:
        assert table.group(table.class_of(label)) == old_group(label), label


def translation_map():
    """{old: (new, display)} of the class rows of the hand-written translation map of release 0.3.0: step S2, then
    the later steps (S3a: the star arm is class miRNA) applied on top. The map is kept in the internal repository
    only, so the test is skipped in the public one."""
    path = DOCKER.parent / 'internal' / 'levels' / 'translation_map_0.3.0.tsv'
    if not path.exists():
        pytest.skip('translation map not in this copy of the repository')
    rows = [line.rstrip('\n').split('\t') for line in path.read_text(encoding='utf-8').splitlines(True)
            if line.strip() and not line.startswith('#')]
    header = rows[0]
    rows = [dict(zip(header, row)) for row in rows[1:]]
    classes = [r for r in rows if r['kind'] == 'class' and not r['old'].startswith('(')]
    result = {r['old']: (r['new'], r['display']) for r in classes if r['step'] == 'S2'}
    for r in classes:
        if r['step'] != 'S2':
            result[r['old']] = (r['new'], r['display'] or result[r['new']][1] if r['new'] in result else r['display'])
    return result


def test_class_ids_and_names_are_those_of_the_translation_map(table):
    for old, (new, display) in translation_map().items():
        if table.row(old) is None:                      # e.g. RiboSpike-Fragment: a countRiboSpike value only
            continue
        assert table.class_of(old) == new, old
        assert table.display_name(new) == display, old


def test_univec_and_star_arm(table):
    """3g: reads on UniVec contigs are the class UniVec; 3a: the star arm is the class miRNA (rank kept)."""
    assert table.class_of('UniVec') == 'UniVec' and table.display_name('UniVec') == 'vector contaminants (UniVec)'
    assert table.gene_table('UniVec') == 'exclude' and table.group('UniVec') == 'UniVec'
    assert table.class_of('miRNA-star') == 'miRNA' and table.rank_map()['miRNA-star'] == 10


def test_rows_of_one_class_must_agree(tmp_path):
    rows = ('rRNA\trrna_bt1\trRNA\trRNA\trRNA\t10\t7\tno\tyes\tyes\td\n'
            'rRNA\tncrna\trRNA\trRNA\tother group\t10\t7\tno\tyes\tyes\td\n')
    with pytest.raises(refTables.TableError, match='differ in group'):
        refTables.ClassTable(write(tmp_path, HEADER + rows), 'human')


def test_keep_values_are_translated(table):
    assert refTables.translate_keep('all', table) == 'all'
    assert refTables.translate_keep('miRNA,miRNA-pre,unannotated,pre-miRNA', table) == 'miRNA,pre-miRNA,genome_unannotated'
    assert refTables.translate_keep('tRNA,all', table) == 'all'


def test_every_bed_label_is_covered(table):
    assert table.check_labels(BED_LABELS) == []


def test_unknown_label_warns_and_gets_rank_10(table, capsys):
    assert table.check_labels(['brand_new_RNA']) == ['brand_new_RNA']
    assert 'brand_new_RNA' in capsys.readouterr().err
    assert table.rank_map(['brand_new_RNA'])['brand_new_RNA'] == 10


def write(tmp_path, text):
    path = tmp_path / 'classes.tsv'
    path.write_text(text, encoding='utf-8')
    return str(path)


HEADER = '\t'.join(refTables.COLUMNS) + '\n'


@pytest.mark.parametrize('row, message', [
    ('rRNA\t*\tr RNA\tr\tg\t1\t7\tno\tyes\tno\td\n', 'spaces'),
    ('rRNA\t*\trRNA\tr\tg\tx\t7\tno\tyes\tno\td\n', 'whole number'),
    ('rRNA\t*\trRNA\tr\tg\t1\t7\tmaybe\tyes\tno\td\n', 'yes or no'),
    ('rRNA\t*\trRNA\tr\tg\t1\t7\tno\tsometimes\tno\td\n', 'gene_table'),
    ('rRNA\t*\trRNA\tr\tg\t1\t7\tno\tyes\tno\td\nrRNA\t*\trRNA\tr\tg\t1\t7\tno\tyes\tno\td\n', 'twice'),
])
def test_broken_tables_are_refused(tmp_path, row, message):
    with pytest.raises(refTables.TableError, match=message):
        refTables.read_table(write(tmp_path, HEADER + row))


def test_wrong_columns_are_refused(tmp_path):
    with pytest.raises(refTables.TableError, match='columns'):
        refTables.read_table(write(tmp_path, 'label\tclass\nrRNA\trRNA\n'))


def test_replacement_gets_the_system_rows(tmp_path):
    path = write(tmp_path, '# my own table\n' + HEADER + 'rRNA\t*\trRNA\trRNA\trRNA\t1\t7\tno\tyes\tyes\tmine\n')
    table = refTables.ClassTable(path, 'human')
    assert table.row('rRNA')['description'] == 'mine'
    assert table.row('unmapped') is not None and table.row('RiboSpike-T1') is not None
    assert table.row('lncRNA') is None                      # a replacement replaces the whole table


def test_rows_of_one_label_must_agree(tmp_path):
    """In 0.3.0 the lookups do not use the tier (critique 2026-10-27): a label with another class or rank in
    another tier is refused instead of silently taking the first row."""
    other_class = ('rRNA\trrna_bt1\trRNA\tr\tg\t1\t7\tno\tyes\tyes\td\n'
                   'rRNA\tncrna\tlncRNA\tl\tg\t2\t7\tno\tyes\tno\td\n')
    with pytest.raises(refTables.TableError, match='label rRNA .* differ in class_id'):
        refTables.ClassTable(write(tmp_path, HEADER + other_class))
    other_rank = ('rRNA\trrna_bt1\trRNA\tr\tg\t1\t7\tno\tyes\tyes\td\n'
                  'rRNA\tncrna\trRNA\tr\tg\t1\t5\tno\tyes\tyes\td\n')
    with pytest.raises(refTables.TableError, match='differ in rank_in_tier'):
        refTables.ClassTable(write(tmp_path, HEADER + other_rank))
    table = refTables.ClassTable(str(DOCKER / 'tables' / 'human.rna_classes.tsv'), 'human')
    assert [r['tier'] for r in table.rows if r['source_label'] == 'rRNA'] == ['rrna_bt1', 'ncrna']   # these agree


def test_tiers_missing_from_the_tier_list_are_named():
    table = refTables.ClassTable(str(DOCKER / 'tables' / 'human.rna_classes.tsv'), 'human')
    built_in = ['univec_bt1', 'ribospike', 'rrna_bt1', 'trnas', 'mirna', 'star', 'pre', 'pirna', 'ncrna', 'cdna',
                'genome']
    assert table.unknown_tiers(built_in) == []
    assert table.unknown_tiers([t for t in built_in if t != 'star']) == ['star']


def test_find_table_precedence(tmp_path):
    explicit = write(tmp_path, HEADER)
    assert refTables.find_table('human', str(tmp_path), explicit) == explicit
    own = tmp_path / 'human.rna_classes.tsv'
    own.write_text(HEADER, encoding='utf-8')
    assert refTables.find_table('human', str(tmp_path)) == str(own)
    assert refTables.find_table('human', None).endswith(os.path.join('tables', 'human.rna_classes.tsv'))
    assert refTables.find_table('zebrafish', None).endswith('default.rna_classes.tsv')
    with pytest.raises(refTables.TableError, match='not found'):
        refTables.find_table('human', None, str(tmp_path / 'missing.tsv'))


def test_tier_id():
    assert refTables.tier_id('/x/genome/human-genome', 'human') == 'genome'
    assert refTables.tier_id('gtf/human-ncrna', 'human') == 'ncrna'
    assert refTables.tier_id('univec/univec_bt1', 'human') == 'univec_bt1'
    assert refTables.tier_id('rrna/mouse-rrna_bt1', 'mouse') == 'rrna_bt1'


def test_gene_symbols_are_applied_only_for_their_bed(tmp_path):
    """The provenance says whether the gene-symbol table was applied (critique 2026-10-27): only when it was
    built for the references BED of the run, as parseMappings checks."""
    bed = tmp_path / 'human.bed'
    bed.write_text('c\t0\t10\tA\tmiRNA\t+\n')
    table = tmp_path / 'human.gene_symbols.tsv'
    table.write_text(f'# references_bed_sha256: {refTables.sha256(str(bed))}\nfeature\thgnc_id\nA\tHGNC:1\n')
    assert refTables.gene_symbols_applied(str(table), str(bed))
    bed.write_text('c\t0\t10\tB\tmiRNA\t+\n')                    # another BED
    assert not refTables.gene_symbols_applied(str(table), str(bed))
    assert not refTables.gene_symbols_applied(str(table), str(tmp_path / 'missing.bed'))
