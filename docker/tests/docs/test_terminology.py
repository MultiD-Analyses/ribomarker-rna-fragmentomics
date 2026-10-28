"""Terminology of release 0.3.0 in the user documentation.

Since 0.3.0 the documentation says "RNA class" (not "RNA type"; "RNA Type" with a capital T is the end
chemistry of the method paper), "feature" for a row of the gene table (not "gene"), and uses the renamed
files, classes, metric ids and settings. The CHANGELOG is not checked: its rename table lists the old names.
A line that names an old version (0.2, "before 0.3.0") or an old name ("old name") may use old terms.
"""
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
DOCUMENTS = ['README.md', 'USER_GUIDE.md', 'CITATION.cff', 'buildDB/README.md',
             'docker/qc/metrics.yaml']   # labels and definitions: shown in the report (critique 2026-10-27)

OLD_TERMS = {
    'RNA type (lower-case t)': r'\bRNA types?\b',
    'biotype': r"(?i)(?<!Ensembl's )\bbiotypes?\b(?!\))",
    'old file name': r'typecounts|typelengths|featuretypes|biotypedistinctseqs|ribotype\.csv|ribolengths',
    'old RNA class': r'miRNA-(pre|star)\b|RiboSpike-Short',
    'old setting keep_types': r'\bkeep_types\b',
    'old metric id': r'assign\.type\.',
    '"genes" for rows of the gene table': r'(?i)\bgenes? detected|\bdetected genes|\bgenes? > ?\d|\bcounts per gene\b',
    'P(T,L) (not the paper\'s notation; v5.12)': r'P\(T, ?L\)',
}


def old_terms(line):
    """The old terms in one line, unless the line names an old version or an old name."""
    if '0.2' in line or 'old name' in line or 'before 0.3.0' in line:
        return []
    return [name for name, pattern in OLD_TERMS.items() if re.search(pattern, line)]


@pytest.mark.parametrize('document', DOCUMENTS)
def test_documentation_uses_the_terms_of_0_3_0(document):
    path = REPO / document
    if not path.exists():
        pytest.skip(f'{document} is not in this copy of the repository')
    found = [f'{document}:{n}: {", ".join(terms)}: {line.strip()[:100]}'
             for n, line in enumerate(path.read_text(encoding='utf-8').splitlines(), start=1)
             for terms in [old_terms(line)] if terms]
    assert not found, '\n'.join(found)


def test_the_check_finds_old_terms():
    assert old_terms('Counts per RNA type and per gene') == ['RNA type (lower-case t)']
    assert old_terms('see `*-typecounts.csv`') == ['old file name']
    assert old_terms('Genes detected with more than 10 reads') == ['"genes" for rows of the gene table']
    assert old_terms('In 0.2.x the file was `*-typecounts.csv`') == []
    assert old_terms('RNA Type (capital T): the end chemistry') == []
    assert old_terms("the label, e.g. Ensembl's biotype") == []
    assert old_terms('% of natural reads (Clark et al. 2026: P(T,L))') != []
    assert old_terms('RiboSpike-Short in typecounts.csv before 0.3.0') == []
