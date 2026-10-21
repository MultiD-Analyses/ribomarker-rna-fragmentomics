"""Same results on every computer (release 0.2.0): every sort in the pipeline scripts is stable, and the rank
sort of the BED keeps the file order for features of equal rank."""
import re
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / 'scripts'
SORT_CALL = re.compile(r'\b(sort_values|sort_index|argsort)\(|\bnp\.sort\(')


def test_every_sort_in_the_scripts_is_stable():
    unstable = []
    for path in sorted(SCRIPTS.glob('*.py')):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            code = line.split('#', 1)[0]
            if SORT_CALL.search(code) and 'kind=' not in code:
                unstable.append(f'{path.name}:{number}')
    assert unstable == [], 'sorts without kind=\'stable\' (their order of ties can differ between computers)'


def test_equal_ranks_keep_the_bed_order():
    pytest.importorskip('numba')
    pytest.importorskip('pysam')
    pd = pytest.importorskip('pandas')
    sys.path.insert(0, str(SCRIPTS))
    import parseMappings

    # 2,000 rows in four biotypes, as in the quicksort test that differed between two computers
    biotypes = ['tRNA', 'pre-mRNA', 'pre-mRNA', 'lncRNA'] * 500
    bed = pd.DataFrame({'chrom': 'c1', 'start': 0, 'end': 50, 'feature': [f'f{i}' for i in range(2000)],
                        'biotype': biotypes, 'strand': '+'})
    ranked = parseMappings.sort_by_rank(bed, {'tRNA': 9, 'lncRNA': 2, 'pre-mRNA': 1})

    assert ranked['biotype'].tolist() == ['tRNA'] * 500 + ['lncRNA'] * 500 + ['pre-mRNA'] * 1000
    for biotype, rows in ranked.groupby('biotype', sort=False):
        numbers = [int(name[1:]) for name in rows['feature']]
        assert numbers == sorted(numbers), f'{biotype}: equal ranks left the BED order'
    assert 'order' not in ranked.columns
