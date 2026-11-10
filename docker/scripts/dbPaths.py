"""The dbpaths file of a references folder: the mapping tiers in mapping order.

<species>.dbpaths.csv has one tier per line, `map_cmd,index,retain`, e.g.
`bowtie2,genome/human-genome,T`. The index paths are relative to the folder of the file (buildDB writes
them so since 2026-10-08), so a references folder can be moved or mounted anywhere. Absolute paths
(older folders) are used as they are.

Used by the Snakefile (which file to read), mapSchedule.py and mapReads.py (where the indexes are).

Since 0.4.0 a run can change the order of the tiers with an order file (starter --tiers, config key tiers), which
names the tiers by tier_id(); the order rules and the mapping schedule are in mapSchedule.py.
"""
import os
import subprocess

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


def read_dbpaths(path, base=None):
    """[(map_cmd, index, retain)] in mapping order, each index as an absolute path (a relative one is
    joined with base, by default the folder of the dbpaths file)."""
    folder = os.path.abspath(base) if base else os.path.dirname(os.path.abspath(path))
    tiers = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            map_cmd, index, retain = [field.strip() for field in line.split(',')]
            tiers.append((map_cmd, os.path.join(folder, index), retain))   # join keeps an absolute index
    return tiers


def tier_id(index, species):
    """Tier id of an index path: its basename without the species prefix (genome/human-genome -> genome)."""
    name = os.path.basename(str(index).rstrip('/'))
    prefix = f'{species}-'
    return name[len(prefix):] if name.startswith(prefix) else name


# map_cmd (mapReads.py) -> the tool that lists an index's sequence names
BOWTIE1_CMDS = ('bowtie', 'bowtie1', 'ribospike', 'mind')
BOWTIE2_CMDS = ('bowtie2', 'bowtie2best', 'bowtie2sense', 'ribospikebt2')
FULLREF_CMDS = ('bowtie1', 'mind', 'ribospike')   # mapReads runs these with bowtie --fullref: SAM has the whole name


def contig_names(map_cmd, index):
    """Reference sequence names of a Bowtie 1 / Bowtie 2 index (`bowtie-inspect -n`), as SAM writes them: up to the
    first whitespace, or the whole name for the commands run with --fullref (critique 2026-10-27: otherwise a tier
    list mapping UniVec with one of them would give no UniVec class)."""
    if map_cmd not in BOWTIE1_CMDS + BOWTIE2_CMDS:
        raise ValueError(f'unknown map_cmd {map_cmd!r} for {index}')
    tool = 'bowtie-inspect' if map_cmd in BOWTIE1_CMDS else 'bowtie2-inspect'
    result = subprocess.run([tool, '-n', index], capture_output=True, text=True, check=True)
    names = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    return names if map_cmd in FULLREF_CMDS else [name.split()[0] for name in names]


def write_contig_tiers(tiers, species, path):
    """contig<TAB>tier for every reference sequence of every tier (release 0.3.0): which tier a read aligned
    in, for reads without a feature (UniVec reads become the class UniVec). tiers: as read_dbpaths returns them."""
    with open(path, 'w') as f:
        f.write('contig\ttier\n')
        for map_cmd, index, _ in tiers:
            tier = tier_id(index, species)
            for name in contig_names(map_cmd, index):
                f.write(f'{name}\t{tier}\n')
