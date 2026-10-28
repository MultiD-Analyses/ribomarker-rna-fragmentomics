#!/usr/bin/env python3
"""
make_gene_symbols.py -- gene symbols and miRBase names for the features of a references BED (release 0.3.0, 3f).

Writes docker/tables/<species>.gene_symbols.tsv (one row per BED feature that has a symbol or a miRBase name) and
<species>.gene_symbols.provenance.md next to it. The pipeline adds the columns to *-featureclasses.csv and
*-genelengths.csv when the references in use have the BED this table was built for (sha256 in the header).

Columns: feature, then for human hgnc_id, hgnc_symbol, for mouse mgi_id, mgi_symbol, then mirbase_name. Several
values for one feature are separated by ';'.

Like buildDB_v1.py: a fixed, recorded set of downloads (SOURCES), each recorded with URL, version, download time
and the sha256 of the bytes; `--check-urls` checks them without downloading. A new set of sources goes into a copy
of this script (make_gene_symbols_v2.py), so a table can always be traced to the script that made it.

  python buildDB/make_gene_symbols.py --species human --check-urls
  python buildDB/make_gene_symbols.py --species human --work WORK --inspect         # download, show the formats
  python buildDB/make_gene_symbols.py --species human --work WORK --references docker/assets

How a feature is joined (the BED of buildDB: Ensembl contigs are named <transcript id>|<name>):
  Ensembl classes (protein_coding, lncRNA, snRNA, snoRNA, misc_RNA, scaRNA, Mt_rRNA, Mt_tRNA, ...):
      transcript id -> gene id (Ensembl GTF) -> HGNC / MGI by Ensembl gene id
  pre-mRNA (<gene name>_intron): gene name -> gene id of the protein_coding gene (GTF) -> HGNC / MGI
  miRNA arms (MirGeneDB, e.g. Hsa-Mir-21_5p, star arms ..._3p*): MirGeneDB GFF -> MIMAT -> miRBase mature name;
      the precursor (the name without the arm) -> MI -> HGNC (human; HGNC lists miRBase precursor accessions)
  pre-miRNA (<precursor>_pre): MirGeneDB GFF -> MI -> miRBase hairpin name; HGNC as above
  tRNA (<species>_tRNA-Ala-AGC-1, a GtRNAdb family) and pre-tRNA (tRNA-Ala-AGC-1-1, one locus): HGNC by GtRNAdb id
      (human only; MGI has no GtRNAdb cross-references, so mouse tRNA symbols stay empty)
  no symbol: piRNA, miRTrace rRNA copies, spike-ins, UniVec

FILE FORMATS, written 2026-10-23 without access to the sites; checked with --inspect on 2026-10-26 (HGNC quarterly
2026-10-06, MGI, MirGeneDB 3.0 GFF) and 2026-10-27 (miRBase 23.1):
  HGNC complete set: tab-separated with a header; columns hgnc_id, symbol, ensembl_gene_id, mirbase, gtrnadb;
      several values in one cell separated by '|'
  MGI MRK_ENSEMBL.rpt: tab-separated, no header; column 1 MGI id, 2 symbol, 6 Ensembl gene id (SOURCES 'columns')
  MirGeneDB GFF: GFF3; attribute ID = MirGeneDB name, Alias = miRBase accession (MI for precursors, MIMAT for arms)
  miRBase mature.fa / hairpin.fa: '>name accession species ...'
SOURCES stay 'PLACEHOLDER' until the URLs and formats are confirmed; the build refuses placeholders. The v1 sources
were verified on 2026-10-27 (miRBase: '>hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p', as read below).
"""
import argparse
import csv
import datetime
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from buildDB_v1 import SOURCES as REFERENCE_SOURCES, check_url, download, sha256_file, utc_now  # noqa: E402

VERSION = 'v1'
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES = os.path.join(REPO, 'docker', 'tables')

# miRBase serves only its current release; the release this builder was checked with (its README, 2026-10-27).
# A later release at the same URL changes the sha256 recorded with the download.
MIRBASE_VERSION = 'miRBase release 23.1 (Sep 2026; public domain)'

# Pinned sources. The GTF is the one the references were built from (buildDB_v1.py SOURCES).
SOURCES = {
    'human': {
        'status': 'verified 2026-10-27: URLs and formats checked with --inspect (HGNC/MGI, GTF, MirGeneDB 2026-10-26; miRBase 2026-10-27)',
        'symbols': 'hgnc',
        'downloads': {
            # the newest quarterly file (link from genenames.org/download/archive/quarterly/tsv/, Martin 2026-10-26)
            'hgnc': {'url': 'https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/'
                            'hgnc_complete_set_2026-10-06.tsv',
                     'version': 'HGNC complete set, quarterly 2026-10-06 (CC0)'},
            'gtf': REFERENCE_SOURCES['human']['downloads']['gtf'],
            'mirgenedb_gff': {'url': 'https://mirgenedb.org/gff/hsa?sort=pos&all=1', 'version': 'MirGeneDB 3.0'},
            'mirbase_mature': {'url': 'https://mirbase.org/download/mature.fa', 'version': MIRBASE_VERSION},
            'mirbase_hairpin': {'url': 'https://mirbase.org/download/hairpin.fa', 'version': MIRBASE_VERSION},
        },
    },
    'mouse': {
        'status': 'verified 2026-10-27: URLs and formats checked with --inspect (HGNC/MGI, GTF, MirGeneDB 2026-10-26; miRBase 2026-10-27)',
        'symbols': 'mgi',
        # columns checked 2026-10-23 against the MGI reports index: 1 MGI Marker Accession ID, 2 Marker Symbol,
        # 6 Ensembl Accession ID (0-based here); terms of use: CC BY 4.0 (THIRD_PARTY_NOTICES, checked 2026-10-27)
        'mgi_columns': {'id': 0, 'symbol': 1, 'ensembl_gene_id': 5},
        'downloads': {
            'mgi': {'url': 'https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt',
                    'version': 'MGI MRK_ENSEMBL.rpt (current; CC BY 4.0)'},
            'gtf': REFERENCE_SOURCES['mouse']['downloads']['gtf'],
            'mirgenedb_gff': {'url': 'https://mirgenedb.org/gff/mmu?sort=pos&all=1', 'version': 'MirGeneDB 3.0'},
            'mirbase_mature': {'url': 'https://mirbase.org/download/mature.fa', 'version': MIRBASE_VERSION},
            'mirbase_hairpin': {'url': 'https://mirbase.org/download/hairpin.fa', 'version': MIRBASE_VERSION},
        },
    },
}
SPECIES_PREFIX = {'human': 'Homo_sapiens_', 'mouse': 'Mus_musculus_'}
ENSEMBL_CONTIG = re.compile(r'^(ENS[A-Z]*T\d+(?:\.\d+)?)\|')
ENSEMBL_GENE = re.compile(r'^ENS[A-Z]*G\d+$')
MIRNA_LABELS = ('miRNA', 'miRNA-star')


# ----------------------------------------------------------------------------- readers
def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_bed(path):
    """[(contig, feature, label)] of a references BED."""
    with open(path) as f:
        return [(cells[0], cells[3], cells[4]) for cells in (line.rstrip('\n').split('\t') for line in f)
                if len(cells) >= 5]


def read_gtf(path):
    """{transcript id: gene id} and {gene name: {gene biotype: [gene ids]}} from an Ensembl GTF."""
    transcript_gene, by_name = {}, {}
    attribute = re.compile(r'(\w+) "([^"]*)"')
    with opener(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cells = line.rstrip('\n').split('\t')
            if len(cells) < 9 or cells[2] not in ('gene', 'transcript'):
                continue
            values = dict(attribute.findall(cells[8]))
            if cells[2] == 'transcript':
                transcript_gene[values['transcript_id']] = values['gene_id']
            elif 'gene_name' in values:
                ids = by_name.setdefault(values['gene_name'], {}).setdefault(values.get('gene_biotype', ''), [])
                ids.append(values['gene_id'])
    return transcript_gene, by_name


def read_hgnc(path):
    """HGNC complete set -> {'ensembl': {gene id: [(hgnc id, symbol)]}, 'mirbase': {MI: [...]}, 'gtrnadb': {...}}."""
    index = {'ensembl': {}, 'mirbase': {}, 'gtrnadb': {}}
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        missing = {'hgnc_id', 'symbol', 'ensembl_gene_id', 'mirbase', 'gtrnadb'} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f'ERROR: {path}: columns not found: {", ".join(sorted(missing))} (check the format with --inspect)')
        for row in reader:
            entry = (row['hgnc_id'], row['symbol'])
            for column, key in (('ensembl_gene_id', 'ensembl'), ('mirbase', 'mirbase'), ('gtrnadb', 'gtrnadb')):
                for value in filter(None, (v.strip() for v in (row[column] or '').split('|'))):
                    index[key].setdefault(value, []).append(entry)
    return index


def read_mgi(path, columns):
    """MGI MRK_ENSEMBL.rpt -> {'ensembl': {gene id: [(mgi id, symbol)]}}."""
    index = {'ensembl': {}, 'mirbase': {}, 'gtrnadb': {}}
    with open(path, encoding='utf-8') as f:
        for line in f:
            cells = line.rstrip('\n').split('\t')
            if len(cells) <= max(columns.values()) or not cells[columns['id']].startswith('MGI:'):
                continue
            for gene in filter(None, (v.strip() for v in cells[columns['ensembl_gene_id']].split())):
                index['ensembl'].setdefault(gene, []).append((cells[columns['id']], cells[columns['symbol']]))
    return index


def read_mirgenedb_gff(path):
    """{MirGeneDB name: miRBase accession} from MirGeneDB's GFF (attributes ID and Alias)."""
    accession = {}
    with opener(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cells = line.rstrip('\n').split('\t')
            if len(cells) < 9:
                continue
            values = dict(part.split('=', 1) for part in cells[8].split(';') if '=' in part)
            if 'ID' in values and values.get('Alias'):
                accession[values['ID']] = values['Alias'].split(',')[0]
    return accession


def read_mirbase_fasta(path):
    """{accession: name} from miRBase mature.fa / hairpin.fa ('>hsa-let-7a-5p MIMAT0000062 Homo sapiens ...')."""
    names = {}
    with opener(path) as f:
        for line in f:
            if line.startswith('>'):
                cells = line[1:].split()
                if len(cells) >= 2:
                    names[cells[1]] = cells[0]
    return names


# ----------------------------------------------------------------------------- joins
def joined(entries):
    """(ids, symbols) as ';'-separated text, in order of first appearance, duplicates removed."""
    seen = []
    for entry in entries:
        if entry not in seen:
            seen.append(entry)
    return ';'.join(e[0] for e in seen), ';'.join(e[1] for e in seen)


def precursor_of(feature):
    """MirGeneDB precursor of an arm (Hsa-Mir-21_5p, Hsa-Let-7-P1b_3p* -> ...) or of a _pre feature."""
    return re.sub(r'_(5p\*?|3p\*?|pre)$', '', feature)


def symbol_rows(species, bed, transcript_gene, genes_by_name, symbols, mirgenedb, mature, hairpin):
    """[(feature, symbol id, symbol, mirbase name)] for every BED feature with at least one value."""
    prefix = SPECIES_PREFIX.get(species, '')
    rows, seen = [], set()
    for contig, feature, label in bed:
        if feature in seen:
            continue
        seen.add(feature)
        entries, mirbase_name = [], ''
        ensembl = ENSEMBL_CONTIG.match(contig)
        if ensembl:
            gene = transcript_gene.get(ensembl.group(1)) or transcript_gene.get(ensembl.group(1).split('.')[0])
            entries = symbols['ensembl'].get(gene, []) if gene else []
        elif label == 'pre-mRNA' and feature.endswith('_intron'):
            name = feature[:-len('_intron')]
            # unnamed genes are named by their gene id in the BED (e.g. ENSG00000285629_intron)
            genes = [name] if ENSEMBL_GENE.match(name) else genes_by_name.get(name, {}).get('protein_coding', [])
            for gene in genes:
                entries += symbols['ensembl'].get(gene, [])
        elif label in MIRNA_LABELS or label == 'miRNA-pre':
            # MirGeneDB ids (--inspect 2026-10-26): precursors <name>_pre (Alias MI...), arms <name>_5p / _3p* (Alias
            # MIMAT...). Arms or precursors with identical sequences are one BED feature 'A_5p|B_5p|...'.
            names = []
            for part in feature.split('|'):
                precursor = mirgenedb.get(precursor_of(part) + '_pre')
                entries += symbols['mirbase'].get(precursor, []) if precursor else []
                if label == 'miRNA-pre':
                    names.append(hairpin.get(precursor, '') if precursor else '')
                else:
                    arm = mirgenedb.get(part)
                    names.append(mature.get(arm, '') if arm else '')
            mirbase_name = ';'.join(n for n in dict.fromkeys(names) if n)
        elif label == 'tRNA' and feature.startswith(prefix + 'tRNA-'):
            family = feature[len(prefix):] + '-'                 # tRNA-Ala-AGC-1 -> loci tRNA-Ala-AGC-1-1, -1-2, ...
            for locus in sorted(k for k in symbols['gtrnadb'] if k.startswith(family)):
                entries += symbols['gtrnadb'][locus]
        elif label == 'pre-tRNA':
            entries = symbols['gtrnadb'].get(feature, [])
        ids, names = joined(entries)
        if ids or mirbase_name:
            rows.append((feature, ids, names, mirbase_name))
    return rows


# ----------------------------------------------------------------------------- modes
def fetch(spec, work, skip_failed=False):
    """Download every source into WORK (reusing complete files already there); {name: (path, record)}.
    A download goes to <file>.part first and gets its name only when complete, so an interrupted download is
    never reused. skip_failed (--inspect): a source that cannot be downloaded is reported and left out."""
    os.makedirs(work, exist_ok=True)
    files = {}
    for name, source in spec['downloads'].items():
        ending = '.gz' if source['url'].split('?')[0].endswith('.gz') else ''
        path = os.path.join(work, f'{name}{ending}')
        if os.path.exists(path) and os.path.getsize(path) > 0:
            modified = datetime.datetime.fromtimestamp(os.path.getmtime(path), datetime.timezone.utc)
            record = {'downloaded_utc': modified.strftime('%Y-%m-%dT%H:%M:%SZ') + ' (file time; reused download)',
                      'download_bytes': os.path.getsize(path)}
        else:
            try:
                record = download(source['url'], path + '.part')
            except OSError as error:            # URLError, timeouts and HTTP errors are OSErrors
                if os.path.exists(path + '.part'):
                    os.remove(path + '.part')
                if not skip_failed:
                    sys.exit(f'ERROR: {name}: {source["url"]} cannot be downloaded ({error})')
                print(f'[SKIP] {name}: {source["url"]} cannot be downloaded ({error})')
                continue
            os.replace(path + '.part', path)
        record.update({'url': source['url'], 'version': source['version'], 'sha256': sha256_file(path)})
        files[name] = (path, record)
    return files


def inspect(files):
    for name, (path, record) in files.items():
        print(f'\n== {name}: {record["url"]}\n   {record["download_bytes"]} bytes, sha256 {record["sha256"][:12]}')
        with opener(path) as f:
            for number, line in enumerate(f):
                if number >= 6:
                    break
                print('   | ' + line.rstrip('\n')[:200])


def build(species, spec, files, references):
    if spec['status'].startswith('PLACEHOLDER'):
        sys.exit(f'ERROR: the {species} sources are placeholders ({spec["status"]}); confirm the URLs (--check-urls) '
                 'and formats (--inspect), then set status to "verified"')
    bed_path = os.path.join(references, f'{species}.bed')
    bed = read_bed(bed_path)
    transcript_gene, genes_by_name = read_gtf(files['gtf'][0])
    if spec['symbols'] == 'hgnc':
        symbols = read_hgnc(files['hgnc'][0])
    else:
        symbols = read_mgi(files['mgi'][0], spec['mgi_columns'])
    mirgenedb = read_mirgenedb_gff(files['mirgenedb_gff'][0])
    mature, hairpin = read_mirbase_fasta(files['mirbase_mature'][0]), read_mirbase_fasta(files['mirbase_hairpin'][0])
    rows = symbol_rows(species, bed, transcript_gene, genes_by_name, symbols, mirgenedb, mature, hairpin)

    key = spec['symbols']
    provenance = {'builder': os.path.basename(__file__), 'builder_version': VERSION,
                  'builder_sha256': sha256_file(os.path.abspath(__file__)), 'built_utc': utc_now(),
                  'species': species, 'references_bed': shown_path(bed_path),
                  'references_bed_sha256': sha256_file(bed_path),
                  'sources': {name: record for name, (_, record) in files.items()},
                  'features_with_values': len(rows), 'bed_features': len({feature for _, feature, _ in bed})}
    os.makedirs(TABLES, exist_ok=True)
    out = os.path.join(TABLES, f'{species}.gene_symbols.tsv')
    with open(out, 'w', newline='') as f:
        f.write(f'# Gene symbols and miRBase names of the features of {species}.bed (release 0.3.0, 3f); made by '
                f'buildDB/{provenance["builder"]} ({VERSION}).\n')
        f.write(f'# references_bed_sha256: {provenance["references_bed_sha256"]}\n')
        f.write(f'# built_utc: {provenance["built_utc"]}\n')
        for name, record in provenance['sources'].items():
            f.write(f'# source {name}: {record["version"]}; {record["url"]}; downloaded {record["downloaded_utc"]}; '
                    f'sha256 {record["sha256"]}\n')
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['feature', f'{key}_id', f'{key}_symbol', 'mirbase_name'])
        writer.writerows(rows)
    write_provenance_md(provenance, os.path.join(TABLES, f'{species}.gene_symbols.provenance.md'))
    print(f'{out}: {len(rows)} of {provenance["bed_features"]} features have a symbol or a miRBase name')
    return out


def write_provenance_md(p, path):
    lines = [f'# Provenance -- gene symbols ({p["species"]})', '',
             f'Built: {p["built_utc"]}  |  builder: {p["builder"]} {p["builder_version"]} '
             f'(sha256 `{p["builder_sha256"][:12]}...`)', '',
             f'References BED: `{p["references_bed"]}` (sha256 `{p["references_bed_sha256"][:12]}...`); '
             f'{p["features_with_values"]} of {p["bed_features"]} features have a value.', '',
             '## Sources', '',
             '| Source | Version | URL | Downloaded (UTC) | sha256 (first 12) |',
             '|--------|---------|-----|------------------|-------------------|']
    for name, s in p['sources'].items():
        lines.append(f'| {name} | {s["version"]} | {s["url"]} | {s["downloaded_utc"]} | `{s["sha256"][:12]}` |')
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    with open(path.replace('.md', '.json'), 'w') as f:
        json.dump(p, f, indent=2)


def shown_path(path):
    """A path for the provenance: relative to the repository, or only the file name when it lies outside (no user
    or computer names in published files; the release leak scan refuses them)."""
    relative = os.path.relpath(os.path.abspath(path), REPO)
    return os.path.basename(path) if relative.startswith('..') else relative.replace(os.sep, '/')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--species', required=True, choices=sorted(SOURCES))
    ap.add_argument('--check-urls', action='store_true', help='check that every source is online (no downloads)')
    ap.add_argument('--work', help='folder for the downloads, in a subfolder per species (reused when present)')
    ap.add_argument('--inspect', action='store_true', help='download and print the first lines of every source')
    ap.add_argument('--references', default=os.path.join(REPO, 'docker', 'assets'),
                    help='references folder with <species>.bed (default: %(default)s)')
    args = ap.parse_args(argv)
    spec = SOURCES[args.species]
    if args.check_urls:
        ok = True
        for name, source in spec['downloads'].items():
            reachable, status, size = check_url(source['url'])
            ok &= reachable
            print(f'[{"OK " if reachable else "DEAD"}] {name:16s} {source["version"]}\n       {source["url"]}  '
                  f'status={status} size={size}')
        return 0 if ok else 1
    if not args.work:
        ap.error('--work is required to download')
    # one subfolder per species: the downloads are named by source (gtf.gz, mirgenedb_gff, ...), so a shared folder
    # would give the second species the first one's files (2026-10-27)
    files = fetch(spec, os.path.join(args.work, args.species), skip_failed=args.inspect)
    if args.inspect:
        inspect(files)
        return 0
    build(args.species, spec, files, args.references)
    return 0


if __name__ == '__main__':
    sys.exit(main())
