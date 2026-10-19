"""Failed tools stop the pipeline.
Run from the repository root: python -m pytest docker/tests/tools"""
import io
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[3]
SCRIPTS = REPO / 'docker' / 'scripts'
sys.path.insert(0, str(SCRIPTS))
from toolRun import ToolError, run  # noqa: E402

PIPELINE_SCRIPTS = ['trimAdapters.py', 'collapseReads.py', 'mapReads.py', 'parseMappings.py', 'countRiboSpike.py']


def test_a_successful_command_passes(tmp_path):
    run(f'echo hello > {tmp_path}/out.txt')
    assert (tmp_path / 'out.txt').read_text() == 'hello\n'


def test_a_failed_command_raises_with_its_status():
    with pytest.raises(ToolError, match=r'exit status 3'):
        run('exit 3')


def test_a_missing_tool_raises():
    with pytest.raises(ToolError, match=r'exit status 127'):
        run('no_such_tool_ribomarker --version 2> /dev/null')


def test_a_failure_inside_a_pipe_counts():
    with pytest.raises(ToolError):
        run('false | cat')          # with sh (os.system) this was status 0


def test_the_message_names_the_error_log_and_the_command(tmp_path):
    log = tmp_path / 'tier.log'
    with pytest.raises(ToolError) as e:
        run(f'bowtie2_missing -x idx -U r.fq -S out.sam 2> {log}')
    assert f'see {log}' in str(e.value) and 'bowtie2_missing -x idx' in str(e.value)


def active_lines(name):
    """Code lines of a pipeline script, comments left out."""
    return [line for line in (SCRIPTS / name).read_text(encoding='utf-8').splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def test_no_unchecked_os_system_left():
    for name in PIPELINE_SCRIPTS:
        assert not [line for line in active_lines(name) if 'os.system(' in line], name


def test_error_branches_exit_with_a_non_zero_status():
    """sys.exit() without an argument exits with status 0, which Snakemake takes for success. The only one
    left is collapseReads' normal end after subsampling."""
    found = {name: sum(line.strip() == 'sys.exit()' for line in active_lines(name)) for name in PIPELINE_SCRIPTS}
    assert found == {'trimAdapters.py': 0, 'collapseReads.py': 1, 'mapReads.py': 0,
                     'parseMappings.py': 0, 'countRiboSpike.py': 0}


def test_cleanup_of_missing_files_does_not_fail():
    """rm became rm -f: a missing temporary file is not a tool failure."""
    for name in PIPELINE_SCRIPTS:
        assert not [line for line in active_lines(name) if re.search(r"run\(f?'rm (?!-f)", line)], name


def make_mapreads(tmp_path, monkeypatch, tiers):
    import mapReads
    monkeypatch.setattr(mapReads.os, 'popen', lambda cmd: io.StringIO('tool version 1.0.0\n'))
    (tmp_path / 'samples.csv').write_text('S1,g,S1.fq.gz\n')
    (tmp_path / 'human.dbpaths.csv').write_text(tiers)
    args = SimpleNamespace(samples=str(tmp_path / 'samples.csv'), dbs=str(tmp_path / 'human.dbpaths.csv'),
                           threads=1, output=str(tmp_path / 'out'), multimap=100)
    return mapReads.MapReads(args)


def test_missing_database_exits_with_status_1(tmp_path, monkeypatch):
    mr = make_mapreads(tmp_path, monkeypatch, 'bowtie2,genome/human-genome,T\n')
    with pytest.raises(SystemExit) as e:
        mr.verify_bowtie_database()
    assert e.value.code == 1


def test_failed_samtools_step_stops_mapping(tmp_path, monkeypatch):
    """samtools on a SAM file that does not exist (status 1, or 127 without samtools) stops the step
    before the SAM file would be deleted."""
    mr = make_mapreads(tmp_path, monkeypatch, 'bowtie2,genome/human-genome,T\n')
    with pytest.raises(ToolError):
        mr.sam_to_bam(str(tmp_path / 'missing.sam'))


def test_launchers_let_snakemake_delete_outputs_of_a_failed_step():
    """With --keep-incomplete a half-written per-sample file survives a failure, and the scripts
    skip samples whose output exists, so a rerun could reuse it."""
    for launcher in (REPO / 'docker' / 'run-pipeline.sh',):
        code = [line for line in launcher.read_text().splitlines() if not line.lstrip().startswith('#')]
        assert not [line for line in code if '--keep-incomplete' in line], launcher
