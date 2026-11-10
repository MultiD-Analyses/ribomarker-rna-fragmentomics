"""
Configuration keys removed in release 0.4.0 (cleanup), for both routes (Snakefile, run).

A removed key with its old default value is accepted with a warning ("obsolete, ignored"); any other value stops
the run, so a configuration that relied on the removed feature fails with a message instead of changing silently.
"""

REMOVED_KEYS = {   # key: (old default value, what to do instead)
    'collapse': (False, 'reads are no longer collapsed before mapping'),
    'subsample_reads': (0, 'subsample the FASTQ files before the run if needed'),
    'trim_fastq': (True, 'every run trims the reads'),
    'trim_only': (False, 'every run trims and maps the reads'),
    'dump_unmapped': (False, 'alignments/<sample>.merged.bam holds the unmapped reads'),
}


def is_old_default(key, value):
    """True if value is the old default of a removed key (True is not taken for 1, nor 0 for False)."""
    default = REMOVED_KEYS[key][0]
    return value == default and type(value) is type(default)


def check_removed_keys(config):
    """(warnings, errors) for the removed keys found in a configuration dictionary."""
    warnings, errors = [], []
    for key, (default, instead) in REMOVED_KEYS.items():
        if key not in config:
            continue
        if is_old_default(key, config[key]):
            warnings.append(f'{key}: {config[key]} is obsolete since release 0.4.0 and ignored; the key can go')
        else:
            errors.append(f'{key}: {config[key]!r}: removed in release 0.4.0 ({instead}); remove the key')
    return warnings, errors
