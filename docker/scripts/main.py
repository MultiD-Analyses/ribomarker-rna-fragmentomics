# main.py (the new single entry point)
import sys
import argparse

# Import each script’s main function
from trimAdapters import main as trimAdapters_main
from collapseReads import main as collapseReads_main
from mapReads import main as mapReads_main
from parseMappings import main as parseMappings_main
from countRiboSpike import main as countRiboSpike_main
from makeBedgraph import main as makeBedgraph_main
from collapseGenes import main as collapseGenes_main
from qcCollect import main as qcCollect_main
from qcReport import main as qcReport_main
from toolRun import ToolError
from runCommand import main as run_main

#from sNAXreport import main as sNAXreport_main

def main():
    parser = argparse.ArgumentParser(
        description="Unified entry point for all pipeline scripts."
    )
    parser.add_argument(
        "--script",
        choices=["trimAdapters","collapseReads","mapReads","parseMappings","countRiboSpike","makeBedgraph","collapseGenes","qcCollect","qcReport"],
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

    # Rebuild sys.argv for the chosen script
    # We'll strip out the known arguments and pass the remainder to the sub-script.
    script_argv = args.script_args[1:]
    print(f"Running {args.script} with arguments: {script_argv}")

    if args.script == "trimAdapters":
        sys.argv = [sys.argv[0]] + script_argv
        trimAdapters_main()
    elif args.script == "collapseReads":
        sys.argv = [sys.argv[0]] + script_argv
        collapseReads_main()
    elif args.script == "mapReads":
        sys.argv = [sys.argv[0]] + script_argv
        mapReads_main()
    elif args.script == "parseMappings":
        sys.argv = [sys.argv[0]] + script_argv
        parseMappings_main()
    elif args.script == "countRiboSpike":
        sys.argv = [sys.argv[0]] + script_argv
        countRiboSpike_main()
    elif args.script == "makeBedgraph":
        sys.argv = [sys.argv[0]] + script_argv
        makeBedgraph_main()
    elif args.script == "collapseGenes":
        sys.argv = [sys.argv[0]] + script_argv
        collapseGenes_main()
    elif args.script == "qcCollect":
        sys.argv = [sys.argv[0]] + script_argv
        qcCollect_main()
    elif args.script == "qcReport":
        sys.argv = [sys.argv[0]] + script_argv
        qcReport_main()
    # elif args.script == "sNAXreport":
    #     sys.argv = [sys.argv[0]] + script_argv
    #     sNAXreport_main()

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":   # main.py run SHEET OUTPUT ...
        sys.exit(run_main(sys.argv[2:]))
    try:
        main()
    except ToolError as error:    # a failed tool: one clear message instead of a traceback
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
