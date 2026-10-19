#!/usr/bin/env python3
"""
update_manifest.py -- bring a references folder to the relative dbpaths file and record the checksum of
every reference file in its MANIFEST, without rebuilding anything.

  python buildDB/update_manifest.py docker/assets --species human            # update
  python buildDB/update_manifest.py docker/assets --species human --check    # only compare files with MANIFEST

Update:
  1. The tiers come from <species>.dbpaths.csv if it exists, else from the MANIFEST's list of indexes
     (relative paths, written by buildDB), else from the older <species>.local.dbpaths.csv (absolute
     paths, made relative; refused if an index lies outside the folder). If the older file exists, its
     tiers must be the same (order, command, index name, retain flag).
  2. Writes <species>.dbpaths.csv with paths relative to the folder.
  3. Rewrites "outputs" of MANIFEST_<species>.json: the dbpaths file, the BED and every index file with
     size and sha256. The previous "outputs" is kept once, as "outputs_before_update"; all other entries
     are kept unchanged.
The older <species>.local.dbpaths.csv / .docker.dbpaths.csv are left in place (the pipeline no longer
reads them); delete them by hand once the new file has been used. Hashing the human set (~9 GB) takes
a few minutes.

buildDB_v1.py uses outputs_entry() for new builds.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import sys

# the files of a bowtie (ebwt) or bowtie2 (bt2, bt2l) index
INDEX_SUFFIXES = [f'{part}.{ext}' for ext in ('ebwt', 'bt2', 'bt2l')
                  for part in ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2')]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def file_entry(folder, rel):
    path = os.path.join(folder, rel)
    return {'file': rel, 'bytes': os.path.getsize(path), 'sha256': sha256_file(path)}


def index_files(folder, index):
    """The files of one index (path relative to folder) that exist."""
    return [index + s for s in INDEX_SUFFIXES if os.path.exists(os.path.join(folder, index + s))]


def outputs_entry(folder, species, tiers):
    """MANIFEST "outputs": dbpaths file, BED and every index file, each with size and sha256.
    tiers = [(map_cmd, index relative to folder, retain)] in mapping order."""
    indices = []
    for map_cmd, index, retain in tiers:
        files = index_files(folder, index)
        if not files:
            sys.exit(f'ERROR: no index files for {index} in {folder}')
        indices.append({'map_cmd': map_cmd, 'path': index, 'retain': retain,
                        'files': [file_entry(folder, f) for f in files]})
    return {'dbpaths': file_entry(folder, f'{species}.dbpaths.csv'),
            'bed': file_entry(folder, f'{species}.bed'),
            'indices': indices}


def read_rows(path):
    with open(path) as f:
        return [tuple(field.strip() for field in line.split(',')) for line in f if line.strip()]


def relative(folder, index):
    """index relative to folder; an absolute index must lie inside the folder."""
    if not os.path.isabs(index):
        return index
    rel = os.path.relpath(index, os.path.abspath(folder))
    if rel.startswith('..'):
        sys.exit(f'ERROR: {index} lies outside {folder}; cannot make it relative')
    return rel


def find_tiers(folder, species, manifest):
    new = os.path.join(folder, f'{species}.dbpaths.csv')
    old = os.path.join(folder, f'{species}.local.dbpaths.csv')
    if os.path.exists(new):
        tiers, source = read_rows(new), new
    elif manifest and manifest.get('outputs', {}).get('indices'):
        tiers = [(i['map_cmd'], i['path'], i['retain']) for i in manifest['outputs']['indices']]
        source = 'the MANIFEST'
    elif os.path.exists(old):
        tiers, source = read_rows(old), old
    else:
        sys.exit(f'ERROR: no dbpaths file and no index list in the MANIFEST in {folder}')
    tiers = [(c, relative(folder, i), r) for c, i, r in tiers]
    if os.path.exists(old):                  # the older file must describe the same tiers
        def names(rows):
            return [(c, os.path.basename(i), r) for c, i, r in rows]
        if names(read_rows(old)) != names(tiers):
            sys.exit(f'ERROR: the tiers of {old} differ from those of {source}; fix by hand')
    return tiers, source


def check(folder, species, manifest):
    """Compare every file in the MANIFEST's outputs with the folder. Returns the number of problems."""
    outputs = manifest.get('outputs', {})
    entries = [outputs[k] for k in ('dbpaths', 'bed') if k in outputs]
    entries += [f for i in outputs.get('indices', []) for f in i.get('files', [])]
    if not entries or 'files' not in (outputs.get('indices') or [{}])[0]:
        print('the MANIFEST has no file checksums; run without --check first')
        return 1
    problems = 0
    for e in entries:
        path = os.path.join(folder, e['file'])
        if not os.path.exists(path):
            print(f'missing:  {e["file"]}'); problems += 1
        elif os.path.getsize(path) != e['bytes'] or sha256_file(path) != e['sha256']:
            print(f'differs:  {e["file"]}'); problems += 1
    print(f'{len(entries) - problems} of {len(entries)} files match MANIFEST_{species}.json')
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder', help='references folder (dbdir), e.g. docker/assets')
    ap.add_argument('--species', default='human')
    ap.add_argument('--check', action='store_true', help='only compare the files with the MANIFEST')
    args = ap.parse_args()

    manifest_path = os.path.join(args.folder, f'MANIFEST_{args.species}.json')
    manifest = None
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            manifest = json.load(f)
    if args.check:
        if manifest is None:
            sys.exit(f'ERROR: {manifest_path} not found')
        sys.exit(1 if check(args.folder, args.species, manifest) else 0)

    tiers, source = find_tiers(args.folder, args.species, manifest)
    dbpaths = os.path.join(args.folder, f'{args.species}.dbpaths.csv')
    with open(dbpaths, 'w') as f:
        f.write('\n'.join(','.join(t) for t in tiers) + '\n')
    print(f'wrote {dbpaths} ({len(tiers)} tiers, from {source})')

    if manifest is None:
        manifest = {'species': args.species,
                    'note': 'created by update_manifest.py: the references were not built by buildDB_v1.py'}
    if 'outputs' in manifest and 'outputs_before_update' not in manifest:
        manifest['outputs_before_update'] = manifest['outputs']
    print('hashing the reference files ...')
    manifest['outputs'] = outputs_entry(args.folder, args.species, tiers)
    manifest['outputs_updated'] = {
        'utc': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'script': os.path.basename(__file__), 'script_sha256': sha256_file(__file__)}
    tmp = manifest_path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(manifest, f, indent=2)
    os.replace(tmp, manifest_path)
    n = sum(len(i['files']) for i in manifest['outputs']['indices'])
    print(f'wrote {manifest_path} ({n} index files, dbpaths file and BED with sha256)')


if __name__ == '__main__':
    main()
