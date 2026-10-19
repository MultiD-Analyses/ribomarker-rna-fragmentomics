"""Run an external tool and stop the step when it fails.

Until 2026-10-08 the scripts ran cutadapt, fastp, bowtie, bowtie2, samtools and seqtk with os.system()
and ignored the exit status, so a failed tool (out of memory, a missing or broken index, a full disk)
went unnoticed and the next step worked on incomplete files.
"""
import re
import subprocess


class ToolError(RuntimeError):
    """An external command exited with a non-zero status."""


def run(cmd):
    """Run a shell command and raise ToolError if it exits with a non-zero status.

    bash with pipefail, so that a failure anywhere in a pipe counts (os.system used sh, where only the
    last command of a pipe counts). Output goes to the same place as with os.system()."""
    status = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd]).returncode
    if status != 0:
        log = re.search(r'2>\s*(\S+)', cmd)                 # the tool's error log, if redirected
        where = f'; see {log.group(1)}' if log else ''
        raise ToolError(f'command failed (exit status {status}){where}:\n  {" ".join(cmd.split())}')
