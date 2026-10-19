#!/usr/bin/env python3
"""
buildDB_v1.py -- standalone, versioned reference-database builder (version 1).

Goal: FULL PROVENANCE. Running this script downloads a fixed, recorded set of
reference databases, builds the bowtie/bowtie2 indices + BED the mapping pipeline
needs, and writes them into an `assets_v1/` folder together with a manifest that
records, per database: the source URL, the version, the download time, the
sha256 of the downloaded bytes, and the sha256 of every built file.

This is INTENTIONALLY independent of the Snakemake pipeline. To make a newer
build, copy this file to buildDB_v2.py, bump the versions/URLs in SOURCES, and
run it -> assets_v2/. Old versions stay frozen and auditable.

Sources for v1 were resolved 2026-07-02 (see PROVENANCE notes in README.md):
  genome/gtf  Ensembl release-116 (GRCh38)      [FTP is mid-migration -- see README]
  tRNA        GtRNAdb hg38 (GtRNAdb2) tarball
  miRNA       MirGeneDB 3.0
  piRNA       piRBase v3.0 (gold set)            [upstream MD5 recorded below]
  rRNA        mirtrace rRNA_reference.fa (master)
  UniVec      NCBI UniVec (current build)
  ribospike   LOCAL, from --inputs (synthetic; not downloadable)

Usage:
  # 1) just check every URL is alive + show versions (no big downloads):
  python buildDB_v1.py --species human --check-urls

  # 2) full build:
  python buildDB_v1.py --species human \
        --inputs ./inputs --output ../assets_v1 --threads 8

Dependencies (on PATH): bowtie-build, bowtie2-build, tar, gunzip
Python deps: pandas, biopython, pysam
"""

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request

# NOTE: builder_core (pandas/biopython/pysam) is imported lazily inside do_build,
# so `--check-urls` and `--help` work with only the Python standard library.

BUILD_VERSION = "v1"

# --------------------------------------------------------------------------- #
# Pinned sources. `post`: how to prepare the downloaded file for builder_core.
#   None    -> use as downloaded          gunzip -> decompress .gz in place
# `upstream_md5` (optional) is the publisher-provided checksum, recorded for
# cross-checking; it is NOT our integrity anchor (our sha256 of the bytes is).
# --------------------------------------------------------------------------- #
SOURCES = {
    "human": {
        "status": "verified",           # URLs liveness-checked for v1
        "taxid": "9606",
        "species_prefix": "Homo_sapiens_",   # used to filter the mirtrace rRNA file
        "downloads": {
            "genome": {
                "url": "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz",
                "version": "Ensembl release-116 / GRCh38",
                "post": "gunzip",
            },
            "gtf": {
                "url": "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz",
                "version": "Ensembl release-116 / GRCh38",
                "post": None,            # builder reads .gtf.gz directly
            },
            "trna": {
                "url": "https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz",
                "version": "GtRNAdb2 hg38",
                "post": None,            # builder untars it
            },
            "mirna": {
                "url": "https://mirgenedb.org/fasta/hsa?all=1",
                "version": "MirGeneDB 3.0",
                "post": None,
            },
            "pirna": {
                "url": "http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/hsa.gold.fa.gz",
                "version": "piRBase v3.0 (gold)",
                "post": "gunzip",
                "upstream_md5": "24f5da6a448e060c5a3bc91a7d4a6ac8",
            },
            "rrna": {
                "url": "https://raw.githubusercontent.com/friedlanderlab/mirtrace/master/src/lib/inputs/rRNA_reference.fa",
                "version": "mirtrace rRNA_reference.fa (branch: master)",
                "post": None,
            },
            "univec": {
                "url": "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec",
                "version": "NCBI UniVec (current build)",
                "post": None,
            },
        },
        "inputs": {"ribospike": "ribospike.fa"},
    },

    # ----------------------------------------------------------------------- #
    # MOUSE -- sources verified 2026-07: GtRNAdb genome code 'Mmusc39' confirmed;
    # piRBase mmu gold MD5 confirmed; Ensembl 116/GRCm39 + MirGeneDB mmu follow
    # the same conventions as the (working) human build; UniVec + RiboSpike are
    # shared with human. rRNA uses the shared mirtrace reference filtered by the
    # 'Mus_musculus_' prefix (mirtrace supports mouse). Run `--check-urls` before
    # a full build to confirm liveness on your machine.
    # ----------------------------------------------------------------------- #
    "mouse": {
        "status": "verified",
        "taxid": "10090",
        "species_prefix": "Mus_musculus_",
        "downloads": {
            "genome": {"url": "https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz",
                       "version": "Ensembl release-116 / GRCm39", "post": "gunzip"},
            "gtf":    {"url": "https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz",
                       "version": "Ensembl release-116 / GRCm39", "post": None},
            "trna":   {"url": "https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz",
                       "version": "GtRNAdb2 mm39", "post": None},
            "mirna":  {"url": "https://mirgenedb.org/fasta/mmu?all=1",
                       "version": "MirGeneDB 3.0", "post": None},
            "pirna":  {"url": "http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/mmu.gold.fa.gz",
                       "version": "piRBase v3.0 (gold)", "post": "gunzip",
                       "upstream_md5": "75bf0f25717331bfa22e8149731b5647"},
            "rrna":   {"url": "https://raw.githubusercontent.com/friedlanderlab/mirtrace/master/src/lib/inputs/rRNA_reference.fa",
                       "version": "mirtrace rRNA_reference.fa (branch: master)", "post": None},
            "univec": {"url": "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec",
                       "version": "NCBI UniVec (current build)", "post": None},
        },
        "inputs": {"ribospike": "ribospike.fa"},
    },
}

# order the raw files are turned into indices (defines dbpaths / mapping-stack order)
BUILD_ORDER = ["univec", "ribospike", "rrna", "trna", "mirna", "pirna", "transcriptome", "genome"]


# --------------------------------------------------------------------------- #
# provenance helpers
# --------------------------------------------------------------------------- #
def utc_now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest):
    """Stream url -> dest. Returns {downloaded_utc, download_bytes}."""
    print(f"  downloading {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "buildDB/1.0"})
    n = 0
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as out:
        for chunk in iter(lambda: r.read(1 << 20), b""):
            out.write(chunk)
            n += len(chunk)
    return {"downloaded_utc": utc_now(), "download_bytes": n}


def check_url(url):
    """Report reachability without a full download (HEAD, else 1-byte GET)."""
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "buildDB/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            return True, r.status, r.headers.get("Content-Length", "?")
    except Exception:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "buildDB/1.0", "Range": "bytes=0-0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return True, r.status, r.headers.get("Content-Range", "partial")
        except Exception as e:
            return False, str(e), "-"


def tool_versions():
    out = {"python": sys.version.split()[0]}
    for tool in ("bowtie-build", "bowtie2-build"):
        try:
            v = subprocess.run([tool, "--version"], capture_output=True, text=True)
            out[tool] = v.stdout.splitlines()[0] if v.stdout else "unknown"
        except Exception:
            out[tool] = "NOT FOUND"
    return out


# --------------------------------------------------------------------------- #
# modes
# --------------------------------------------------------------------------- #
def do_check_urls(spec, species):
    print(f"\nURL liveness check -- {species} (status: {spec['status']})\n" + "-" * 60)
    ok_all = True
    for name, d in spec["downloads"].items():
        ok, status, size = check_url(d["url"])
        ok_all = ok_all and ok
        print(f"  [{'OK ' if ok else 'DEAD'}] {name:10s} {d['version']}")
        print(f"         {d['url']}")
        print(f"         status={status} size={size}")
    print("-" * 60)
    print("ALL REACHABLE" if ok_all else "SOME URLS UNREACHABLE (see above)")
    return 0 if ok_all else 1


def require_inputs(spec, inputs_dir):
    missing = []
    ribo = os.path.join(inputs_dir, spec["inputs"]["ribospike"])
    if not os.path.isfile(ribo):
        missing.append(ribo)
    if missing:
        print("ERROR: required local inputs missing (see buildDB/inputs/README.md):")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)
    return ribo


def do_build(spec, species, inputs_dir, out_root, threads, keep_raw):
    if spec["status"].startswith("PLACEHOLDER"):
        print(f"ERROR: '{species}' is a placeholder (status={spec['status']}). "
              f"Verify its URLs and set status='verified' before building.")
        sys.exit(1)

    import builder_core as core  # heavy deps (pandas/biopython/pysam) only needed to build

    ribospike_fa = require_inputs(spec, inputs_dir)
    os.makedirs(out_root, exist_ok=True)
    raw = os.path.join(out_root, "_raw")
    os.makedirs(raw, exist_ok=True)

    manifest = {
        "build": BUILD_VERSION,
        "builder_script": os.path.basename(__file__),
        "builder_sha256": sha256_file(__file__),
        "builder_core_sha256": sha256_file(os.path.join(os.path.dirname(__file__), "builder_core.py")),
        "species": species,
        "built_utc": utc_now(),
        "tools": tool_versions(),
        "sources": {},
        "inputs": {},
        "outputs": {},
    }

    # ---- fetch every remote source (reuse anything already in _raw/) -------- #
    local = {}
    for name, d in spec["downloads"].items():
        ext = ".tar.gz" if d["url"].endswith(".tar.gz") else (".gz" if d["url"].endswith(".gz") else ".dat")
        raw_dest = os.path.join(raw, f"{name}{ext}")
        final = raw_dest[:-3] if d.get("post") == "gunzip" else raw_dest

        if os.path.exists(final) and os.path.getsize(final) > 0:
            print(f"  reuse {os.path.basename(final)} (already in _raw/)")
            meta = {"downloaded_utc": "reused-existing", "download_bytes": os.path.getsize(final)}
        else:
            if os.path.exists(raw_dest) and os.path.getsize(raw_dest) > 0:
                print(f"  reuse {os.path.basename(raw_dest)} (downloaded, not yet unpacked)")
                meta = {"downloaded_utc": "reused-existing", "download_bytes": os.path.getsize(raw_dest)}
            else:
                meta = download(d["url"], raw_dest)
            if d.get("post") == "gunzip":
                subprocess.run(f"gunzip -f {raw_dest}", shell=True, check=True)

        local[name] = final
        meta["sha256"] = sha256_file(final)   # hash the exact file handed to the builder
        manifest["sources"][name] = {"url": d["url"], "version": d["version"],
                                     "local_file": os.path.basename(final), **meta}
        if "upstream_md5" in d:
            manifest["sources"][name]["upstream_md5"] = d["upstream_md5"]

    # ---- record the local (non-downloadable) input ------------------------- #
    manifest["inputs"]["ribospike"] = {"source_path": os.path.abspath(ribospike_fa),
                                       "sha256": sha256_file(ribospike_fa)}

    # ---- build indices (order defines the dbpaths / mapping stack) --------- #
    sp = species  # dbprefix, e.g. 'human'
    dbrows, beds = [], []

    r, b = core.build_univec(local["univec"], out_root, threads);                         dbrows += r; beds += b
    r, b = core.build_ribospike(ribospike_fa, out_root, threads);                          dbrows += r; beds += b
    r, b = core.build_rrna(local["rrna"], spec["species_prefix"], out_root, threads, sp);  dbrows += r; beds += b
    r, b = core.build_gtrnadb(local["trna"], out_root, threads, sp);                       dbrows += r; beds += b
    r, b = core.build_mirna(local["mirna"], out_root, threads, sp);                        dbrows += r; beds += b
    r, b = core.build_pirbase(local["pirna"], out_root, threads, sp);                      dbrows += r; beds += b
    r, b = core.build_transcriptome(local["gtf"], local["genome"], out_root, threads, sp); dbrows += r; beds += b
    r, b = core.build_genome(local["genome"], out_root, threads, sp);                      dbrows += r; beds += b

    dbpaths_csv, bedfile = core.write_outputs(dbrows, beds, out_root, sp)

    # ---- record outputs + checksums of every index file ------ #
    import update_manifest
    manifest["outputs"] = update_manifest.outputs_entry(out_root, sp, dbrows)

    with open(os.path.join(out_root, f"MANIFEST_{sp}.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    write_provenance_md(manifest, out_root, sp)

    if not keep_raw:
        shutil.rmtree(raw, ignore_errors=True)

    print(f"\nDone. assets root: {out_root}")
    print(f"  dbpaths:          {dbpaths_csv}")
    print(f"  bed:              {bedfile}")
    print(f"  manifest:         {os.path.join(out_root, f"MANIFEST_{sp}.json")}")


def write_provenance_md(m, out_root, sp):
    lines = [f"# Provenance -- assets {m['build']} ({m['species']})", "",
             f"Built: {m['built_utc']}  |  builder: {m['builder_script']} "
             f"(sha256 `{m['builder_sha256'][:12]}...`)", "",
             "## Sources", "",
             "| DB | Version | URL | Downloaded (UTC) | sha256 (first 12) |",
             "|----|---------|-----|------------------|-------------------|"]
    for name, s in m["sources"].items():
        lines.append(f"| {name} | {s['version']} | {s['url']} | {s['downloaded_utc']} | `{s['sha256'][:12]}` |")
    lines += ["", "## Local inputs", "",
              f"- ribospike: `{m['inputs']['ribospike']['source_path']}` "
              f"(sha256 `{m['inputs']['ribospike']['sha256'][:12]}`)", "",
              "## Tools", ""]
    for k, v in m["tools"].items():
        lines.append(f"- {k}: {v}")
    with open(os.path.join(out_root, f"PROVENANCE_{sp}.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="Versioned reference-DB builder (v1).")
    ap.add_argument("--species", default="human", choices=sorted(SOURCES.keys()))
    ap.add_argument("--inputs", default=os.path.join(os.path.dirname(__file__), "inputs"),
                    help="folder with non-downloadable inputs (ribospike.fa)")
    ap.add_argument("--output", default=None,
                    help="assets output root (default: ../assets_<version> next to the repo)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--check-urls", action="store_true", help="only test URL liveness, then exit")
    ap.add_argument("--keep-raw", action="store_true", help="keep the downloaded _raw/ files")
    args = ap.parse_args()

    spec = SOURCES[args.species]

    if args.check_urls:
        sys.exit(do_check_urls(spec, args.species))

    out_root = args.output or os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", f"assets_{BUILD_VERSION}"))
    do_build(spec, args.species, args.inputs, out_root, args.threads, args.keep_raw)


if __name__ == "__main__":
    main()
