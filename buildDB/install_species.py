#!/usr/bin/env python3
"""
install_species.py -- copy one species' references from a separate build folder into a shared references
folder (e.g. docker/assets), without touching the other species.

  python buildDB/buildDB_v1.py --species mouse --output ../assets_mouse_build      # build in an EMPTY folder
  python buildDB/install_species.py ../assets_mouse_build docker/assets --species mouse

Copied: every file whose name starts with the species ("mouse.", "mouse-"), the GtRNAdb files of its
genome ("mm39-"), MANIFEST_<species>.json and PROVENANCE_<species>.md. Not copied: the indexes shared by
all species (UniVec, RiboSpike), and the unprefixed intermediate files (trnas.fa, pretrnas.bed, logs ...),
which the pipeline does not read. A shared index of the new build that differs from the installed one is
reported; the installed one is kept (the other species uses it).

Files of the species that are replaced, and its older dbpaths files, are moved (not deleted) into the
build folder, <build>/_replaced_<species>_<UTC time>/ (not into the references folder, which goes into
the Docker image). Afterwards the MANIFEST's "outputs" are recomputed from the
installed files (update_manifest.py), so `update_manifest.py <folder> --species <species> --check` passes.
"""
import argparse
import datetime as dt
import json
import os
import shutil
import sys

import update_manifest


def genome_code(manifest):
    """'mm39' from the GtRNAdb archive name in the MANIFEST (…/mm39-tRNAs.tar.gz)."""
    name = os.path.basename(manifest['sources']['trna']['url'])
    return name.split('-')[0]


def species_files(build, species, code):
    """Paths (relative to build) of the files to install."""
    prefixes = (f'{species}.', f'{species}-', f'{code}-')
    found = []
    for root, _, files in os.walk(build):
        if any(part.startswith('_') for part in os.path.relpath(root, build).split(os.sep)):
            continue                          # _raw/ downloads, _replaced_*/ from an earlier install
        for name in files:
            if name.startswith(prefixes) or name in (f'MANIFEST_{species}.json', f'PROVENANCE_{species}.md'):
                found.append(os.path.relpath(os.path.join(root, name), build))
    return sorted(found)


def shared_indexes(build, target, species, tiers):
    """Indexes without the species prefix (shared by all species): (index, same as installed?)."""
    report = []
    for _, index, _ in tiers:
        if os.path.basename(index).startswith(f'{species}-'):
            continue
        new = update_manifest.index_files(build, index)
        same = bool(new) and all(os.path.exists(os.path.join(target, f)) and
                                 update_manifest.sha256_file(os.path.join(build, f)) ==
                                 update_manifest.sha256_file(os.path.join(target, f)) for f in new)
        report.append((index, same))
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('build', help='folder of the separate build (buildDB_v1.py --output)')
    ap.add_argument('target', help='shared references folder, e.g. docker/assets')
    ap.add_argument('--species', required=True)
    args = ap.parse_args()
    sp = args.species

    with open(os.path.join(args.build, f'MANIFEST_{sp}.json')) as f:
        manifest = json.load(f)
    if manifest.get('species') != sp:
        sys.exit(f'ERROR: the MANIFEST in {args.build} is for {manifest.get("species")!r}, not {sp!r}')
    tiers = update_manifest.read_rows(os.path.join(args.build, f'{sp}.dbpaths.csv'))
    files = species_files(args.build, sp, genome_code(manifest))
    for _, index, _ in tiers:          # every species-specific index must be among the files
        if os.path.basename(index).startswith(f'{sp}-') and not update_manifest.index_files(args.build, index):
            sys.exit(f'ERROR: index {index} of {sp}.dbpaths.csv is missing in {args.build}')

    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    keep = os.path.join(args.build, f'_replaced_{sp}_{stamp}')     # outside the references folder
    old_dbpaths = [f'{sp}.local.dbpaths.csv', f'{sp}.docker.dbpaths.csv']
    moved = 0
    for rel in files + [f for f in old_dbpaths if os.path.exists(os.path.join(args.target, f))]:
        dest = os.path.join(args.target, rel)
        if os.path.exists(dest):
            os.makedirs(os.path.dirname(os.path.join(keep, rel)), exist_ok=True)
            shutil.move(dest, os.path.join(keep, rel))
            moved += 1
    for rel in files:
        os.makedirs(os.path.dirname(os.path.join(args.target, rel)), exist_ok=True)
        shutil.copy2(os.path.join(args.build, rel), os.path.join(args.target, rel))
    print(f'installed {len(files)} files of {sp} into {args.target}; '
          f'{moved} replaced files moved to {keep if moved else "(none)"}')

    for index, same in shared_indexes(args.build, args.target, sp, tiers):
        print(f'shared index {index}: ' + ('same as installed' if same else
              'DIFFERS from the installed one (newer download); the installed one is kept and used'))

    # outputs from the installed files, so that --check holds for what the pipeline will read
    target_manifest = os.path.join(args.target, f'MANIFEST_{sp}.json')
    with open(target_manifest) as f:
        installed = json.load(f)
    installed['outputs'] = update_manifest.outputs_entry(args.target, sp, tiers)
    installed['installed'] = {'utc': stamp, 'from': os.path.abspath(args.build),
                              'script': os.path.basename(__file__),
                              'script_sha256': update_manifest.sha256_file(__file__)}
    with open(target_manifest + '.tmp', 'w') as f:
        json.dump(installed, f, indent=2)
    os.replace(target_manifest + '.tmp', target_manifest)
    print(f'wrote {target_manifest}; check with: python buildDB/update_manifest.py {args.target} --species {sp} --check')


if __name__ == '__main__':
    main()
