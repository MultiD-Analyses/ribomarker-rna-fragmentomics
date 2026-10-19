#!/usr/bin/env python3
"""Stages the reference files that go into the pipeline image (see the Dockerfile).

Only the files listed in MANIFEST_<species>.json (dbpaths, BED, index files), plus the MANIFEST and the
PROVENANCE, each checked against the size and sha256 in the MANIFEST, for each species. Genome index files
go to big/, everything else to rest/, so the Dockerfile can give each large file its own layer. Files that
both species list (UniVec, RiboSpike) must be identical and are staged once. Files are hard-linked when
possible (no extra disk space), otherwise copied. Run from the repository root:

  python docker/stage_references.py                 # docker/assets, human and mouse -> staging/references
  python docker/stage_references.py --no-verify     # sizes only, no sha256 (fast)
"""
import argparse
import hashlib
import json
import os
import shutil
import sys


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def listed_files(manifest):
    """[{'file', 'bytes', 'sha256'}, ...] of the MANIFEST's outputs."""
    out = manifest["outputs"]
    return [out["dbpaths"], out["bed"]] + [f for index in out["indices"] for f in index["files"]]


def place(src, dst):
    if os.path.exists(dst):                      # listed by both species (checked identical against both MANIFESTs)
        return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assets", default="docker/assets")
    ap.add_argument("--species", nargs="+", default=["human", "mouse"])
    ap.add_argument("--out", default="staging/references")
    ap.add_argument("--no-verify", action="store_true", help="check sizes only, not sha256")
    args = ap.parse_args()
    out = args.out
    if "staging" not in os.path.abspath(out):
        sys.exit(f"ERROR: refusing to replace {out}: the output folder must be a staging folder")
    if os.path.exists(out):
        shutil.rmtree(out)                       # links only; the files in docker/assets are not touched

    problems, staged = [], {}
    for species in args.species:
        with open(os.path.join(args.assets, f"MANIFEST_{species}.json")) as fh:
            manifest = json.load(fh)
        for entry in listed_files(manifest):
            src = os.path.join(args.assets, entry["file"])
            if not os.path.isfile(src):
                problems.append(f"missing: {entry['file']}")
                continue
            if os.path.getsize(src) != entry["bytes"] or (not args.no_verify and sha256(src) != entry["sha256"]):
                problems.append(f"differs from MANIFEST_{species}.json: {entry['file']}")
                continue
            if staged.get(entry["file"], entry["sha256"]) != entry["sha256"]:
                problems.append(f"listed by two species with different contents: {entry['file']}")
                continue
            part = "big" if entry["file"].startswith("genome/") else "rest"
            place(src, os.path.join(out, part, entry["file"]))
            staged[entry["file"]] = entry["sha256"]
        for name in (f"MANIFEST_{species}.json", f"PROVENANCE_{species}.md"):
            place(os.path.join(args.assets, name), os.path.join(out, "rest", name))

    if problems:
        print("\n".join(problems))
        sys.exit(f"ERROR: {len(problems)} file(s) do not match a MANIFEST; nothing should be built from {out}")
    total = sum(os.path.getsize(os.path.join(args.assets, f)) for f in staged)
    print(f"staged {len(staged)} files of {', '.join(args.species)} ({total / 1e9:.2f} GB) in {out} "
          f"({'sizes' if args.no_verify else 'sha256'} checked against the MANIFESTs)")


if __name__ == "__main__":
    main()
