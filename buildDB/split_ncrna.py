#!/usr/bin/env python3
"""
split_ncrna.py -- split the Ensembl ncRNA reference into four sub-references (release 0.4.0).

  python buildDB/split_ncrna.py docker/assets --species human \
         --saved-manifest <copy of the MANIFEST before the split>/MANIFEST_human.json [--threads 8]

The four sub-references, in this order at the position of the old ncrna tier (order decided by Martin):

  ncrna_mt      Mt_tRNA, Mt_rRNA
  ncrna_snsno   snRNA, snoRNA, scaRNA
  ncrna_other   misc_RNA, rRNA (nuclear rRNA gene loci)
  ncrna_lnc     lncRNA

On an existing references folder (no full rebuild):
  1. checks that the record names and lengths of gtf/<species>.ncRNA.fa are those of the old index
     (bowtie2-inspect -s gtf/<species>-ncrna);
  2. cuts gtf/<species>.<sub-reference>.fa from it: every record unchanged and in its order, its biotype
     taken from gtf/<species>.ncRNA.bed (column 5, by contig name in column 1);
  3. builds the four Bowtie 2 indexes gtf/<species>-<sub-reference>;
  4. replaces the ncrna row of <species>.dbpaths.csv by the four rows, then runs update_manifest.py;
  5. records the split (this script's sha256, bowtie2-build version, records per sub-reference, FASTA
     checksums) in MANIFEST_<species>.json and PROVENANCE_<species>.md;
  6. checks that every other file of the MANIFEST, the BED included, has the checksum of the saved 0.3.0
     MANIFEST (so the gene-symbol tables stay valid).
The old combined index gtf/<species>-ncrna is left in place; it is removed by hand once the split is checked.
builder_core.py uses split_fasta() for new builds. Needs bowtie2 and bowtie2-build (the pipeline environment).
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys

SUB_REFERENCES = [
    ('ncrna_mt', ('Mt_tRNA', 'Mt_rRNA')),
    ('ncrna_snsno', ('snRNA', 'snoRNA', 'scaRNA')),
    ('ncrna_other', ('misc_RNA', 'rRNA')),
    ('ncrna_lnc', ('lncRNA',)),
]
SUB_REFERENCE_OF = {biotype: name for name, biotypes in SUB_REFERENCES for biotype in biotypes}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_records(fasta):
    """[(name, record text)] in file order; the text is the header line and the sequence lines, unchanged."""
    records = []
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                records.append([line[1:].split()[0], line])
            elif records:
                records[-1][1] += line
            elif line.strip():
                raise SystemExit(f'ERROR: {fasta}: sequence before the first header')
    return [tuple(r) for r in records]


def sequence_length(text):
    return sum(len(line.strip()) for line in text.split('\n')[1:])


def read_biotypes(bed):
    """contig -> biotype, from the transcript BED (contig in column 1, biotype in column 5)."""
    biotype = {}
    with open(bed) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 5:
                biotype[fields[0]] = fields[4]
    return biotype


def split_fasta(fasta, bed, out_fasta):
    """Writes out_fasta(sub-reference) for each sub-reference; returns {sub-reference: {biotype: records}}.
    Stops if a record has no biotype or one outside every sub-reference, or if a build check fails."""
    records = read_records(fasta)
    biotype = read_biotypes(bed)
    names = [name for name, _ in records]
    if len(set(names)) != len(names):
        raise SystemExit(f'ERROR: {fasta}: record names are not unique')
    if set(names) != set(biotype):
        raise SystemExit(f'ERROR: the records of {fasta} and the contigs of {bed} differ')
    outside = sorted({biotype[n] for n in names} - set(SUB_REFERENCE_OF))
    if outside:
        raise SystemExit(f'ERROR: biotypes outside every sub-reference: {outside}')
    counts = {sub: {} for sub, _ in SUB_REFERENCES}
    for sub, _ in SUB_REFERENCES:
        with open(out_fasta(sub), 'w') as out:
            for name, text in records:
                if SUB_REFERENCE_OF[biotype[name]] == sub:
                    out.write(text)
                    counts[sub][biotype[name]] = counts[sub].get(biotype[name], 0) + 1
    check_split(records, biotype, out_fasta, counts)
    return counts


def check_split(records, biotype, out_fasta, counts):
    """Build checks: every old record in exactly one new FASTA, byte-identical and in the same
    relative order; each FASTA holds only its biotypes; the counts add up."""
    position = {name: i for i, (name, _) in enumerate(records)}
    text = dict(records)
    seen = []
    for sub, biotypes in SUB_REFERENCES:
        new = read_records(out_fasta(sub))
        order = [position[name] for name, _ in new]
        if order != sorted(order):
            raise SystemExit(f'ERROR: {out_fasta(sub)}: records not in their old order')
        for name, record in new:
            if record != text[name] or biotype[name] not in biotypes:
                raise SystemExit(f'ERROR: {out_fasta(sub)}: record {name} changed or of another biotype')
        seen += [name for name, _ in new]
    if sorted(seen) != sorted(position) or len(seen) != len(records):
        raise SystemExit('ERROR: the new FASTAs do not hold every old record exactly once')
    if sum(n for c in counts.values() for n in c.values()) != len(records):
        raise SystemExit('ERROR: the counts per biotype do not add up to the records')


def check_index(index, fasta):
    """Check: the old index holds the records of the FASTA (names and lengths, in order)."""
    out = subprocess.run(['bowtie2-inspect', '-s', index], check=True, capture_output=True, text=True).stdout
    indexed = [(f[1].split()[0], int(f[2])) for f in (line.split('\t') for line in out.splitlines())
               if f[0].startswith('Sequence-')]
    expected = [(name, sequence_length(text)) for name, text in read_records(fasta)]
    if indexed != expected:
        differ = sum(a != b for a, b in zip(indexed, expected)) + abs(len(indexed) - len(expected))
        raise SystemExit(f'ERROR: {index} and {fasta} differ ({len(indexed)} vs {len(expected)} records, '
                         f'{differ} names or lengths differ)')


def tool_version(tool):
    out = subprocess.run([tool, '--version'], check=True, capture_output=True, text=True).stdout
    return out.splitlines()[0].split('/')[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder', help='references folder (dbdir), e.g. docker/assets')
    ap.add_argument('--species', default='human')
    ap.add_argument('--saved-manifest', required=True, help='copy of the MANIFEST, made before the split')
    ap.add_argument('--threads', type=int, default=8)
    args = ap.parse_args()
    folder, species = args.folder, args.species
    gtf = os.path.join(folder, 'gtf')
    fasta, bed = os.path.join(gtf, f'{species}.ncRNA.fa'), os.path.join(gtf, f'{species}.ncRNA.bed')
    old_index = os.path.join(gtf, f'{species}-ncrna')
    dbpaths = os.path.join(folder, f'{species}.dbpaths.csv')
    old_row_path = f'gtf/{species}-ncrna'

    # refusals before anything is written
    for older in (f'{species}.local.dbpaths.csv', f'{species}.docker.dbpaths.csv'):
        if os.path.exists(os.path.join(folder, older)):
            raise SystemExit(f'ERROR: remove {os.path.join(folder, older)} first (an older tier list; '
                             'update_manifest.py would refuse the new tiers)')
    rows = [line.strip().split(',') for line in open(dbpaths) if line.strip()]
    if [r[1] for r in rows].count(old_row_path) != 1:
        raise SystemExit(f'ERROR: {dbpaths} has no single row {old_row_path} (already split?)')
    saved = json.load(open(args.saved_manifest))

    print(f'1. checking {old_index} against {fasta}')
    check_index(old_index, fasta)

    print('2. cutting the four FASTAs')
    def out_fasta(sub):
        return os.path.join(gtf, f'{species}.{sub}.fa')
    counts = split_fasta(fasta, bed, out_fasta)
    for sub, by_biotype in counts.items():
        print(f'   {sub}: {sum(by_biotype.values())} records {by_biotype}')

    print('3. building the four Bowtie 2 indexes')
    log = os.path.join(gtf, f'{species}.ncrna_split.bowtie-build.log')
    for sub, _ in SUB_REFERENCES:
        with open(log, 'a') as f:
            subprocess.run(['bowtie2-build', '--threads', str(args.threads), out_fasta(sub),
                            os.path.join(gtf, f'{species}-{sub}')], check=True, stdout=f, stderr=subprocess.STDOUT)

    print(f'4. {dbpaths}: the ncrna row replaced by the four sub-references; update_manifest.py')
    new_rows = []
    for row in rows:
        if row[1] == old_row_path:
            new_rows += [[row[0], f'gtf/{species}-{sub}', row[2]] for sub, _ in SUB_REFERENCES]
        else:
            new_rows.append(row)
    with open(dbpaths, 'w') as f:
        f.write(''.join(','.join(r) + '\n' for r in new_rows))
    here = os.path.dirname(os.path.abspath(__file__))
    subprocess.run([sys.executable, os.path.join(here, 'update_manifest.py'), folder, '--species', species], check=True)

    print('5. recording the split in the MANIFEST and PROVENANCE')
    manifest_path = os.path.join(folder, f'MANIFEST_{species}.json')
    manifest = json.load(open(manifest_path))
    manifest['ncrna_split'] = {
        'utc': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'script': os.path.basename(__file__), 'script_sha256': sha256_file(__file__),
        'bowtie2-build': tool_version('bowtie2-build'),
        'input_fasta': {'file': f'gtf/{species}.ncRNA.fa', 'sha256': sha256_file(fasta)},
        'sub_references': [{'name': sub, 'index': f'gtf/{species}-{sub}', 'biotypes': counts[sub],
                            'fasta': f'gtf/{species}.{sub}.fa', 'fasta_sha256': sha256_file(out_fasta(sub))}
                           for sub, _ in SUB_REFERENCES]}
    tmp = manifest_path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(manifest, f, indent=2)
    os.replace(tmp, manifest_path)
    split = manifest['ncrna_split']
    with open(os.path.join(folder, f'PROVENANCE_{species}.md'), 'a') as f:
        f.write(f"\n## ncRNA split (release 0.4.0)\n\n{split['utc']}: `{split['script']}` (sha256 "
                f"`{split['script_sha256'][:12]}`), {split['bowtie2-build']}. The Ensembl ncRNA reference "
                f"`gtf/{species}.ncRNA.fa` (sha256 `{split['input_fasta']['sha256'][:12]}`) cut into four, every "
                "record unchanged and in its order:\n\n| Sub-reference | Records | Biotypes |\n|---|---|---|\n")
        for sub in split['sub_references']:
            biotypes = ', '.join(f'{b} {n}' for b, n in sub['biotypes'].items())
            f.write(f"| `{sub['index']}` | {sum(sub['biotypes'].values())} | {biotypes} |\n")

    print('6. every other file against the saved MANIFEST')
    new_files = {e['file']: e['sha256'] for e in files_of(manifest['outputs'])}
    changed = [e['file'] for e in files_of(saved['outputs'])
               if not e['file'].startswith(f'{old_row_path}.') and e['file'] != f'{species}.dbpaths.csv'
               and new_files.get(e['file']) != e['sha256']]
    if changed:
        raise SystemExit(f'ERROR: files changed against the saved MANIFEST: {changed}')
    print(f'   {len(files_of(saved["outputs"]))} files of the saved MANIFEST: all others unchanged')
    print(f'done. The old index {old_index}.* stays until you remove it.')


def files_of(outputs):
    return [outputs['dbpaths'], outputs['bed']] + [f for i in outputs['indices'] for f in i['files']]


if __name__ == '__main__':
    main()
