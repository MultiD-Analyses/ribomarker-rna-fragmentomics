"""The dbpaths file of a references folder: the mapping tiers in mapping order.

<species>.dbpaths.csv has one tier per line, `map_cmd,index,retain`, e.g.
`bowtie2,genome/human-genome,T`. The index paths are relative to the folder of the file (buildDB writes
them so since 2026-10-08), so a references folder can be moved or mounted anywhere. Absolute paths
(older folders) are used as they are.

Used by the Snakefile (which file to read) and by mapReads.py (the tiers).
"""
import os

CONTAINER_ENV = 'RIBOMARKER_CONTAINER'   # set to 1 in the Docker image (docker/Dockerfile)


def find_dbpaths(dbdir, species, in_container=None):
    """Path of the dbpaths file to use for this references folder and species.

    <species>.dbpaths.csv if present. Otherwise, outside a container, the older
    <species>.local.dbpaths.csv (absolute paths) with a warning. In a container the older files are
    refused: <species>.docker.dbpaths.csv points at the image's own references, so a mounted folder
    would be ignored without notice."""
    if in_container is None:
        in_container = os.environ.get(CONTAINER_ENV) == '1'
    new = os.path.join(dbdir, f'{species}.dbpaths.csv')
    if os.path.exists(new):
        return new
    old = os.path.join(dbdir, f'{species}.local.dbpaths.csv')
    if os.path.exists(old) and not in_container:
        print(f'WARNING: {new} not found; using the older {old} (absolute paths). '
              f'Update the folder with: python buildDB/update_manifest.py {dbdir} --species {species}')
        return old
    raise FileNotFoundError(
        f'{new} not found. The references folder is missing or older than 2026-10-08; update it with '
        f'"python buildDB/update_manifest.py <folder> --species {species}" or rebuild it with buildDB.')


def read_dbpaths(path):
    """[(map_cmd, index, retain)] in mapping order, each index as an absolute path (a relative one is
    joined with the folder of the dbpaths file)."""
    folder = os.path.dirname(os.path.abspath(path))
    tiers = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            map_cmd, index, retain = [field.strip() for field in line.split(',')]
            tiers.append((map_cmd, os.path.join(folder, index), retain))   # join keeps an absolute index
    return tiers
