# Contributing

Thank you for your interest.

## Questions, problems, ideas

Please open an issue. For a problem with a run, include:

- the command you ran and the last lines of its output;
- the pipeline version (in the QC report under *Software versions*, or in `run_status.json`);
- `qc/report_technical.html` if you can share it (see the user guide, *Sending the report*).

## Pull requests

This repository is published from an internal repository, where development happens. Pull requests are welcome as proposals: we apply accepted changes in the internal repository, credit you in the `CHANGELOG.md` entry, and they appear here with the next release. The pull request itself is then closed. By opening a pull request you agree that your contribution is licensed under the MIT License of this project.

Before proposing a change, please run the unit tests and the demo:

```bash
conda activate ribomarker              # environment from docker/env.lock.txt, see USER_GUIDE.md §10
python -m pytest docker/tests
bash docker/make_demo.sh && python docker/scripts/main.py run --demo /tmp/demo_check
```

A change that alters results (counts, coverage or QC numbers) needs a reason and a description of the effect.
