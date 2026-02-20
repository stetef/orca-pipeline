# ORCA and CORVUS pipeline automation

This directory contains template inputs and helper scripts to run ORCA geometry optimizations and then prepare CORVUS/FEFF inputs. The typical flow is:

1) `prepare-orca.py` → run ORCA geometry optimization
2) `prepare-corvus.py` → convert ORCA output to FEFF inputs and set up CORVUS
3) `plot_xas.py` → plot XANES/EXAFS from completed CORVUS runs

## Scripts

`prepare-orca.py`
- Copies `orca-template.in` (or `orca-template-h-only.in`) and fills placeholders
- Cleans XYZ files and writes sidecar comments
- Generates ORCA input files per structure
- Generates `generated-<name>-orca.script` from `orca-qsub.script`
- Use `--dry-run` to skip `qsub` submission
- Key args: `path`, `--cys`, `--his`, `--out-dir`, `--dry-run`

`prepare-corvus.py`
- Requires ORCA `.hess` output
- Converts `.hess` → `.dym`, then runs `dym2feffinp` to build FEFF inputs
- Copies `corvus-template.in` and `corvus-qsub.script` (submit manually)
- Key args: `path`, `--out-dir`

`plot_xas.py`
- Takes a finished CORVUS run directory
- Plots XANES/EXAFS and converts chi(k) to R space via Larch
- Key args: `path`, `--out-dir`

## Templates
- `orca-template.in`: ORCA input template (placeholders filled by `prepare-orca.py`)
- `orca-template-h-only.in`: ORCA H-only input template
- `orca-qsub.script`: ORCA submission script template (filled by `prepare-orca.py`)
- `corvus-template.in`: CORVUS input template (filled by `prepare-corvus.py`)
- `corvus-qsub.script`: CORVUS submission script (copied by `prepare-corvus.py`)
- `corvus-wrapper-qsub.script`: CORVUS stage wrapper template (filled by `run-batch-pipeline.py`)
- `postprocess-qsub.script`: Batch postprocess template (filled by `run-batch-pipeline.py`)

## Quick examples
```bash
python prepare-orca.py /path/to/xyz --cys 3 --his 1 --out-dir /path/to/output
python prepare-corvus.py /path/to/orca/output --out-dir /path/to/output
python plot_xas.py /path/to/corvus/run --out-dir /path/to/plots
```