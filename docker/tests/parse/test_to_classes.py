"""parseMappings.to_classes (release 0.3.0): the per-label tallies of a sample become per-class tallies; labels of
one class are merged (needed from S3a on, when the miRNA star arm joins the class miRNA)."""
import sys
from pathlib import Path

import pytest

DOCKER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DOCKER / 'scripts'))


def test_labels_are_renamed_and_merged():
    pytest.importorskip('numba')                         # parseMappings imports numba
    import parseMappings
    class_of = {'unannotated': 'genome_unannotated', 'miRNA-star': 'miRNA', 'miRNA': 'miRNA'}.get
    counts = {'miRNA': 5, 'unannotated': 2, 'miRNA-star': 3}
    lengths = {'miRNA': {21: 4, 22: 1}, 'unannotated': {30: 2}, 'miRNA-star': {21: 1, 22: 2}}
    distinct = {'miRNA': {'ACGT': 4, 'ACGA': 1}, 'unannotated': {'TTTT': 2}, 'miRNA-star': {'ACGT': 1, 'GGGG': 2}}
    new_counts, new_lengths, new_distinct = parseMappings.to_classes(class_of, counts, lengths, distinct)
    assert new_counts == {'miRNA': 8, 'genome_unannotated': 2}
    assert list(new_counts) == ['miRNA', 'genome_unannotated']          # first-seen order kept (output row order)
    assert new_lengths == {'miRNA': {21: 5, 22: 3}, 'genome_unannotated': {30: 2}}
    assert new_distinct['miRNA'] == {'ACGT': 5, 'ACGA': 1, 'GGGG': 2}  # one sequence counted once per class
    assert len(new_distinct['miRNA']) == 3
