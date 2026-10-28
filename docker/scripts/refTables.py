"""The RNA class table: how the labels of a reference BED become RNA classes (release 0.3.0).

One table per species, `<species>.rna_classes.tsv`, tab-separated, '#' lines are comments, columns:

  source_label   the label in the 5th column of the references BED (or a label the pipeline sets itself:
                 UniVec, unannotated, unmapped, RiboSpike-*, miND, antisense, Natural, Unmapped)
  tier           the tier the label comes from (index basename without the species prefix, e.g. ncrna), or *
  class_id       the class written into output tables (no spaces)
  display_name   the name in the QC report and documentation
  group          the group in the QC report's composition and tier plots
  display_order  order of classes in the report (integer, ascending)
  rank_in_tier   which of overlapping features of one reference wins (higher wins; parseMappings)
  positions      yes: reads get _5p/_3p/_internal/_whole by their place in the feature (tRNA today)
  gene_table     yes: in the gene table; drop: precursor/intronic, left out; exclude: not a gene (spikes,
                 unannotated, unmapped)
  counts_as_rrna yes: left out of the "without rRNA" library size and top-10 share (QC)
  description    one line for the report and the glossary

Which table is used (first found): an explicit file (starter --rna-classes / config key rna_classes), the
references folder's own <species>.rna_classes.tsv, the pipeline's docker/tables/<species>.rna_classes.tsv,
docker/tables/default.rna_classes.tsv. A table replaces the default as a whole; the system rows
(SYSTEM_LABELS) are added from the default when a replacement lacks them.
BED labels missing from the table get class = label, the group UNKNOWN_GROUP, rank 10 and a warning: the rank
parseMappings gave unknown biotypes before 0.3.0.
"""
import csv
import hashlib
import os
import sys

import dbPaths

TABLES_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tables'))
COLUMNS = ['source_label', 'tier', 'class_id', 'display_name', 'group', 'display_order', 'rank_in_tier',
           'positions', 'gene_table', 'counts_as_rrna', 'description']
SYSTEM_LABELS = ('UniVec', 'unannotated', 'unmapped', 'RiboSpike-T1', 'RiboSpike-T2', 'RiboSpike-T3', 'RiboSpike-T4',
                 'RiboSpike-Short', 'miND', 'antisense', 'Natural', 'Unmapped')
UNKNOWN_RANK = 10            # rank of BED labels missing from the table (parseMappings before 0.3.0)
UNKNOWN_ORDER = 1000                       # report position of classes missing from the table (last)
UNKNOWN_GROUP = 'lncRNA and other ncRNA'   # composition group of unlisted labels (as before 0.3.0); their read-length
                                           # curves are in the report's "Other RNAs" (qcReport.OWN_LENGTH_BUTTON)
GENE_TABLE_VALUES = ('yes', 'drop', 'exclude')


class TableError(ValueError):
    """A class table that cannot be used."""


def read_table(path):
    """Rows of a class table as dicts (all values strings), checked for columns and values."""
    rows = []
    with open(path, newline='', encoding='utf-8') as f:
        lines = [line for line in f if line.strip() and not line.startswith('#')]
    reader = csv.DictReader(lines, delimiter='\t')
    if reader.fieldnames is None or [c.strip() for c in reader.fieldnames] != COLUMNS:
        raise TableError(f'{path}: columns must be {", ".join(COLUMNS)} (tab-separated), found '
                         f'{", ".join(reader.fieldnames or [])}')
    seen = set()
    for n, row in enumerate(reader, start=2):
        row = {k: (v or '').strip() for k, v in row.items()}
        where = f'{path}, row {n} ({row["source_label"]})'
        if not row['source_label'] or not row['tier'] or not row['class_id']:
            raise TableError(f'{where}: source_label, tier and class_id are required')
        if any(c.isspace() for c in row['class_id']):
            raise TableError(f'{where}: class_id must not contain spaces ("{row["class_id"]}")')
        for column in ('display_order', 'rank_in_tier'):
            try:
                row[column] = int(row[column])
            except ValueError:
                raise TableError(f'{where}: {column} must be a whole number, found "{row[column]}"')
        for column in ('positions', 'counts_as_rrna'):
            if row[column] not in ('yes', 'no'):
                raise TableError(f'{where}: {column} must be yes or no, found "{row[column]}"')
        if row['gene_table'] not in GENE_TABLE_VALUES:
            raise TableError(f'{where}: gene_table must be one of {", ".join(GENE_TABLE_VALUES)}')
        key = (row['source_label'], row['tier'])
        if key in seen:
            raise TableError(f'{where}: the pair ({key[0]}, {key[1]}) appears twice')
        seen.add(key)
        rows.append(row)
    return rows


def default_table(species):
    """The pipeline's own table for this species, or the species-neutral one."""
    own = os.path.join(TABLES_DIR, f'{species}.rna_classes.tsv')
    return own if os.path.exists(own) else os.path.join(TABLES_DIR, 'default.rna_classes.tsv')


def find_table(species, references_dir=None, explicit=None):
    """Path of the class table to use (precedence in the module text)."""
    if explicit:
        if not os.path.exists(explicit):
            raise TableError(f'class table not found: {explicit}')
        return explicit
    if references_dir:
        own = os.path.join(references_dir, f'{species}.rna_classes.tsv')
        if os.path.exists(own):
            return own
    return default_table(species)


CLASS_COLUMNS = ('display_name', 'group', 'display_order', 'positions', 'gene_table', 'counts_as_rrna')
LABEL_COLUMNS = ('class_id', 'rank_in_tier')   # one value per label (lookups do not use the tier in 0.3.0)


class ClassTable:
    """The class table of a run. The tier column records which reference a label comes from (or '*' for all).
    In 0.3.0 the pipeline does not use it for lookups: a label has one class, rank and position setting in every
    tier, so the rows of one label must agree on LABEL_COLUMNS (checked here; critique 2026-10-27). Tier-specific
    rows can come with the split of the Ensembl ncRNA reference (0.4.0). row(label, tier) is kept for that.

    Two views: parseMappings works with the BED labels (rank_map, position_labels, class_of) and writes class
    ids; every later step reads class ids (the class-level methods gene_table, position_classes, rrna_classes,
    not_gene_classes, group, display_name). Rows that share a class id must agree on CLASS_COLUMNS."""

    def __init__(self, path, species='human'):
        self.path = path
        self.rows = read_table(path)
        default = default_table(species)
        if os.path.abspath(path) != os.path.abspath(default):        # a replacement: add missing system rows
            present = {row['source_label'] for row in self.rows}
            self.rows += [row for row in read_table(default)
                          if row['source_label'] in SYSTEM_LABELS and row['source_label'] not in present]
        self.warnings = []
        self._by_label = {}
        self._by_class = {}
        for row in self.rows:
            self._by_label.setdefault(row['source_label'], []).append(row)
            first_of_label = self._by_label[row['source_label']][0]
            disagree = [c for c in LABEL_COLUMNS if row[c] != first_of_label[c]]
            if disagree:
                raise TableError(f'{path}: rows of label {row["source_label"]} (tiers {first_of_label["tier"]}, '
                                 f'{row["tier"]}) differ in {", ".join(disagree)}; in this version a label has '
                                 'one class and rank in every tier')
            first = self._by_class.setdefault(row['class_id'], row)
            different = [c for c in CLASS_COLUMNS if row[c] != first[c]]
            if different:
                raise TableError(f'{path}: rows of class {row["class_id"]} ({first["source_label"]}, '
                                 f'{row["source_label"]}) differ in {", ".join(different)}')

    def row(self, label, tier=None):
        """The row of a label (and tier), or None."""
        rows = self._by_label.get(label, [])
        if tier is not None:
            for row in rows:
                if row['tier'] == tier:
                    return row
            for row in rows:
                if row['tier'] == '*':
                    return row
            return None
        return rows[0] if rows else None

    def unknown_tiers(self, tier_ids):
        """Tiers named in the table that the run's tier list does not have ('*' excluded): a misspelt tier, or a
        table made for other references. Harmless in 0.3.0 (lookups do not use the tier), but worth a warning."""
        return sorted({row['tier'] for row in self.rows if row['tier'] != '*'} - set(tier_ids))

    def check_labels(self, labels, quiet=False):
        """Warn about labels (e.g. of a BED) that the table does not cover; returns them. quiet: no message
        (the caller reports them)."""
        missing = sorted({label for label in labels if label not in self._by_label})
        if missing and not quiet:
            message = (f'WARNING: {len(missing)} label(s) not in the class table {self.path}: {", ".join(missing)}; '
                       f'used as their own class, group "{UNKNOWN_GROUP}", rank {UNKNOWN_RANK}')
            self.warnings.append(message)
            print(message, file=sys.stderr)
        return missing

    # --- labels (parseMappings) ------------------------------------------------------------------------

    def rank_map(self, labels=()):
        """{label: rank} for parseMappings: every table label (first row) plus the given labels; labels
        missing from the table get UNKNOWN_RANK."""
        ranks = {label: rows[0]['rank_in_tier'] for label, rows in self._by_label.items()}
        for label in labels:
            ranks.setdefault(label, UNKNOWN_RANK)
        return ranks

    def position_labels(self):
        """Labels whose reads get _5p/_3p/_internal/_whole (parseMappings)."""
        return {label for label, rows in self._by_label.items() if rows[0]['positions'] == 'yes'}

    def class_of(self, label, tier=None):
        """Class id of a label; a label missing from the table is its own class. A class id given as label is
        returned as it is, so old and new names can be mixed (e.g. config key keep_classes)."""
        row = self.row(label, tier)
        if row is not None:
            return row['class_id']
        return label

    # --- classes (every step after parseMappings) -------------------------------------------------------

    def class_row(self, class_id):
        return self._by_class.get(class_id)

    def classes_where(self, column, value):
        return {class_id for class_id, row in self._by_class.items() if row[column] == value}

    def position_classes(self):
        """Classes whose features carry _5p/_3p/_internal/_whole (collapseGenes sums them back)."""
        return self.classes_where('positions', 'yes')

    def gene_table(self, class_id):
        """yes / drop / exclude for collapseGenes; RiboSpike* classes missing from the table are excluded,
        other missing classes are kept (as before 0.3.0)."""
        row = self.class_row(class_id)
        if row is not None:
            return row['gene_table']
        return 'exclude' if str(class_id).startswith('RiboSpike') else 'yes'

    def not_gene_classes(self):
        """Classes that are not in the gene table (QC check C7)."""
        return {class_id for class_id, row in self._by_class.items() if row['gene_table'] != 'yes'}

    def rrna_classes(self):
        return self.classes_where('counts_as_rrna', 'yes')

    def group(self, class_id):
        """Report group of a class; RiboSpike* classes missing from the table are spike-ins."""
        row = self.class_row(class_id)
        if row is not None:
            return row['group']
        return 'spike-ins' if str(class_id).startswith('RiboSpike') else UNKNOWN_GROUP

    def display_name(self, class_id):
        row = self.class_row(class_id)
        return row['display_name'] if row is not None else class_id

    def display_order(self, class_id):
        """Position of a class in the report (column display_order); classes missing from the table come last."""
        row = self.class_row(class_id)
        try:
            return int(row['display_order']) if row is not None else UNKNOWN_ORDER
        except ValueError:
            return UNKNOWN_ORDER

    def description(self, class_id):
        row = self.class_row(class_id)
        return row.get('description', '') if row is not None else ''


def gene_symbols_table(species):
    """docker/tables/<species>.gene_symbols.tsv (buildDB/make_gene_symbols.py, release 0.3.0, 3f), or None."""
    path = os.path.join(TABLES_DIR, f'{species}.gene_symbols.tsv')
    return path if os.path.exists(path) else None


def gene_symbols_applied(symbols_path, bed_path):
    """Whether parseMappings fills the gene-symbol columns: the table was built for this BED (the sha256 in its
    header, as parseMappings.add_gene_symbols checks). False when the BED is missing."""
    built_for = None
    with open(symbols_path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('# references_bed_sha256:'):
                built_for = line.split(':', 1)[1].strip()
                break
            if not line.startswith('#'):
                break
    return built_for is not None and os.path.exists(bed_path) and sha256(bed_path) == built_for


def renamed_keep(value, table):
    """The names of a coverage filter that translate_keep changes: [(old name, class id)]. The 0.3.0 plan (D2)
    accepts old class names with a warning; e.g. miRNA-star became part of miRNA, so it now keeps every miRNA."""
    names = [name.strip() for name in str(value).split(',') if name.strip() and name.strip() != 'all']
    return [(name, table.class_of(name)) for name in names if table.class_of(name) != name]


def translate_keep(value, table):
    """The coverage filter (config key keep_classes; keep_types before 0.3.0): 'all', or comma-separated class ids.
    Old names (BED labels and class ids before 0.3.0, e.g. miRNA-pre, unannotated) are translated with the class
    table; order kept, duplicates removed."""
    names = [name.strip() for name in str(value).split(',') if name.strip()]
    if not names or 'all' in names:
        return 'all'
    classes = []
    for name in names:
        rna_class = table.class_of(name)
        if rna_class not in classes:
            classes.append(rna_class)
    return ','.join(classes)


def load(species, references_dir=None, explicit=None):
    """The ClassTable to use for a run."""
    return ClassTable(find_table(species, references_dir, explicit), species)


tier_id = dbPaths.tier_id     # tier ids: index basename without the species prefix


# --- provenance -----------------------------------------------------------------------------------------

def sha256(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def provenance(source, default, copy=None):
    """What qc_metrics.json and run_status.json record about a table a run used: its file name, sha256 (of the
    copy in the run folder when there is one) and whether it replaced the default."""
    used = copy if copy is not None and os.path.exists(copy) else source
    return {'file': os.path.basename(str(source)), 'sha256': sha256(used),
            'replaced': os.path.abspath(str(source)) != os.path.abspath(str(default))}
