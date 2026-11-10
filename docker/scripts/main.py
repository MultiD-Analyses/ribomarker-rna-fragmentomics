# main.py (the new single entry point)
import sys
import argparse

# Import each script’s main function
from trimAdapters import main as trimAdapters_main
from mapReads import main as mapReads_main
from parseMappings import main as parseMappings_main
from countRiboSpike import main as countRiboSpike_main
from makeBedgraph import main as makeBedgraph_main
from collapseGenes import main as collapseGenes_main
from qcCollect import main as qcCollect_main
from qcReport import main as qcReport_main
from toolRun import ToolError
from runCommand import main as run_main

# the pipeline scripts, by the name given with --script
SCRIPTS = {
    "trimAdapters": trimAdapters_main,
    "mapReads": mapReads_main,
    "parseMappings": parseMappings_main,
    "countRiboSpike": countRiboSpike_main,
    "makeBedgraph": makeBedgraph_main,
    "collapseGenes": collapseGenes_main,
    "qcCollect": qcCollect_main,
    "qcReport": qcReport_main,
}


def main():
    parser = argparse.ArgumentParser(
        description="Unified entry point for all pipeline scripts."
    )
    parser.add_argument(
        "--script",
        choices=list(SCRIPTS),
        required=True,
        help="Which pipeline script to run"
    )
    # Everything after `--` can be passed through to the called script
    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Arguments to forward to the chosen script"
    )

    args = parser.parse_args()

    # Rebuild sys.argv for the chosen script: the arguments after `--`
    script_argv = args.script_args[1:]
    print(f"Running {args.script} with arguments: {script_argv}")
    sys.argv = [sys.argv[0]] + script_argv
    SCRIPTS[args.script]()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":   # main.py run SHEET OUTPUT ...
        sys.exit(run_main(sys.argv[2:]))
    try:
        main()
    except ToolError as error:    # a failed tool: one clear message instead of a traceback
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
