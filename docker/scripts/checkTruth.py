"""
checkTruth.py -- compare a pipeline run on the synthetic reads (level A) with truth.tsv.
Used by `run --demo` and by the level A test.

Checked tables (found anywhere below --outdir):
  *-trimstats.csv     dimer_reads, trimmed_reads, merged_reads (PE) / <15nt_reads (SE)
  *-mappingstats.csv  Mapped, Unmapped, Unique, Multi, Total
  *-classcounts.csv   reads per RNA class (class ids; typecounts.csv before 0.3.0)
  *-genecounts.csv    reads per gene
  *-ribospike.csv     reads per spike
  *-spikelengths.csv  read length per spike type (ribolengths.csv before 0.3.0)
For the count tables, a row that is absent counts as 0, and any row with a
non-zero count that is NOT in the truth is reported as a failure too.
truth.tsv rows for other tables (input, tier_aligned, tier_aligned_multi, qc)
are for the QC collector (docker/tests/qc/test_truth.py) and are skipped here.
Exploratory samples (rows of table "exploratory", sample S3): only their trimming
is checked; for each exploratory category the reads counted for its gene are
printed (outcome reported, not checked; no "unexpected" check for these samples).

Usage (from the repository root):
  python test/levels/check_truth.py --truth test/synthetic_qc/bundle/truth/truth.tsv \
      --mode pe --outdir results/levelA_pe
Writes <outdir>/check_truth.tsv; exit code 1 if anything fails.
"""
import argparse
import csv
import glob
import os
import sys

COUNT_TABLES = ("classcounts", "genecounts", "ribospike", "spikelengths")
STAT_TABLES = ("trimstats", "mappingstats")


def find_one(outdir, pattern):
    hits = glob.glob(os.path.join(outdir, "**", pattern), recursive=True)
    if len(hits) > 1:   # parseMappings also writes a *-ribospike.csv; use the one from countRiboSpike
        hits = [h for h in hits if f"{os.sep}ribospike{os.sep}" in h]
    return hits[0] if len(hits) == 1 else None


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


def load_observed(outdir):
    """{(table, sample, key): value} from the pipeline outputs."""
    obs = {}
    files = {t: find_one(outdir, f"*-{t}.csv") for t in STAT_TABLES + COUNT_TABLES}
    for table, path in files.items():
        if path is None:
            print(f"WARNING: no unique *-{table}.csv below {outdir}")
            continue
        rows = read_csv(path)
        header = rows[0]
        if table in STAT_TABLES:                  # one row per sample, one column per statistic
            for row in rows[1:]:
                for col, value in zip(header[1:], row[1:]):
                    obs[(table, row[0], col)] = value
        elif table == "spikelengths":             # sample, spike_type, 0, 1, 2, ...
            for row in rows[1:]:
                sample, spike_type = row[0], row[1]
                for length, value in zip(header[2:], row[2:]):
                    obs[(table, sample, f"{spike_type}:{length}")] = value
        else:                                     # one row per feature, one column per sample
            for row in rows[1:]:
                for sample, value in zip(header[1:], row[1:]):
                    obs[(table, sample, row[0])] = value
    return obs


def as_number(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def check(truth_path, mode, outdir):
    """Compare the run in outdir with the truth of one mode; writes <outdir>/check_truth.tsv.
    Returns (results, path of check_truth.tsv); a result row is [status, table, sample, key, expected,
    observed, categories] with status PASS, FAIL or INFO."""
    with open(truth_path, newline="") as f:
        truth = [r for r in csv.DictReader(f, delimiter="\t") if r["mode"] == mode]
    observed = load_observed(outdir)

    results, expected_keys = [], set()
    for t in truth:
        if t["table"] not in STAT_TABLES + COUNT_TABLES:
            continue
        key = (t["table"], t["sample"], t["key"])
        expected_keys.add(key)
        raw = observed.get(key)
        value = as_number(raw) if raw is not None else (0 if t["table"] in COUNT_TABLES else None)
        status = "PASS" if value == int(t["expected"]) else "FAIL"
        results.append([status, t["table"], t["sample"], t["key"], t["expected"], value, t["categories"]])

    # non-zero rows in the count tables that the truth does not expect
    exploratory = [t for t in truth if t["table"] == "exploratory"]
    exploratory_samples = {t["sample"] for t in exploratory}
    samples = {t["sample"] for t in truth} - exploratory_samples
    for (table, sample, key), raw in observed.items():
        if table in COUNT_TABLES and sample in samples and (table, sample, key) not in expected_keys:
            value = as_number(raw)
            if value:
                results.append(["FAIL", table, sample, key, 0, value, "unexpected"])

    # exploratory categories: reads counted for the category's gene (each category has its own gene)
    for t in exploratory:
        category, gene = t["key"].split(":", 1)
        value = as_number(observed.get(("genecounts", t["sample"], gene))) or 0
        results.append(["INFO", "exploratory", t["sample"], t["key"], t["expected"], value, category])

    out = os.path.join(outdir, "check_truth.tsv")
    with open(out, "w") as f:
        f.write("status\ttable\tsample\tkey\texpected\tobserved\tcategories\n")
        f.writelines("\t".join(str(x) for x in r) + "\n" for r in results)
    return results, out


def report(results, mode, out):
    """Print the summary line, the failures and the exploratory outcomes; returns the exit code."""
    failures = [r for r in results if r[0] == "FAIL"]
    passes = [r for r in results if r[0] == "PASS"]
    print(f"\n=== {mode.upper()}: {len(passes)} PASS, {len(failures)} FAIL ===")
    for r in failures:
        print(f"  FAIL  {r[1]:13s} {r[2]:4s} {r[3]:32s} expected {r[4]!s:>6}  observed {r[5]!s:>6}  ({r[6]})")
    for r in results:
        if r[0] == "INFO":
            print(f"  exploratory {r[2]} {r[3]:16s} {r[5]} of {r[4]} reads counted for the gene")
    print(f"details: {out}")
    return 1 if failures else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default="test/synthetic_qc/bundle/truth/truth.tsv")
    ap.add_argument("--mode", required=True, choices=["pe", "se"])
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    results, out = check(args.truth, args.mode, args.outdir)
    return report(results, args.mode, out)


if __name__ == "__main__":
    sys.exit(main())
