"""Configuration keys removed in release 0.4.0 (docker/scripts/removedKeys.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import removedKeys  # noqa: E402


def test_old_defaults_give_a_warning_and_no_error():
    old = {'collapse': False, 'subsample_reads': 0, 'trim_fastq': True, 'trim_only': False, 'dump_unmapped': False}
    warnings, errors = removedKeys.check_removed_keys(old)
    assert errors == [] and len(warnings) == len(old)
    assert all('obsolete since release 0.4.0 and ignored' in w for w in warnings)


def test_other_values_are_errors():
    warnings, errors = removedKeys.check_removed_keys({'collapse': True, 'subsample_reads': 1000, 'trim_fastq': False,
                                                       'trim_only': True, 'dump_unmapped': True})
    assert warnings == [] and len(errors) == 5 and all('removed in release 0.4.0' in e for e in errors)


def test_the_type_counts():
    """0 is not False and False is not 0: a YAML value of the wrong type is not the old default."""
    _, errors = removedKeys.check_removed_keys({'collapse': 0, 'subsample_reads': False})
    assert len(errors) == 2


def test_other_keys_are_not_checked():
    assert removedKeys.check_removed_keys({'multimap': 100, 'remove_multimap': True}) == ([], [])
