# ORCA and CORVUS pipeline automation

This directory contains template inputs and helper scripts to run ORCA geometry optimizations and then prepare CORVUS/FEFF inputs. The typical flow is:

1) `script-prepare-orca.py` → run ORCA geometry optimization
2) `script-prepare-corvus.py` → convert ORCA output to FEFF inputs and set up CORVUS
3) `script-process-feff-output.py` → postprocess FEFF output (plots + chi(R) + `dw.dat`)

## Scripts

`script-prepare-orca.py`
- Copies ORCA input templates from `orca-templates/` and fills placeholders
- Cleans XYZ files and writes sidecar comments
- Generates ORCA input files per structure
- Generates `generated-<name>-orca.script` from `orca-qsub.script`
- Use `--dry-run` to skip `qsub` submission
- Key args: `path`, `--out-dir`, mode flags, `--dry-run`
- Reads charge/multiplicity from XYZ header line 2 using:
	- `CHARGE_ROUNDED=<int>` or `ROUNDED_CHARGE=<int>`
	- `MULTIPLICITY=<int>`

`script-prepare-corvus.py`
- Requires ORCA `.hess` output
- Converts `.hess` → `.dym`, then runs `dym2feffinp` to build FEFF inputs
- Copies `corvus-template.in` and `corvus-qsub.script` (submit manually)
- Key args: `path`, `--out-dir`

`script-process-feff-output.py`
- Takes a finished CORVUS/FEFF run directory
- Plots XANES/EXAFS and converts chi(k) to R space via Larch
- Writes `dw.dat` and copies selected outputs into `output-<name>` directories
- Key args: `path`, `--out-dir`

Additional helper scripts with `script-` prefix are included for reporting and packaging:
- `script-count-imag-freq.py`
- `script-extract-orca-compute-times.py`
- `script-prepare-files-for-download.py`

## Template Inventory
ORCA input templates in `orca-templates/`:
- `orca-template-ca-fixed.in`
- `orca-template-h-only.in`
- `orca-template-single-point.in`
- `orca-template-no-constraints.in`
- `orca-template-backbone-charges.in`
- `orca-template-xtb-free.in`
- `orca-template-xtb-constrained.in`
- `orca-template-quick.in` (available, not selected by default CLI mode flags)
- `orca-template-ca-fixed-p450.in` (available, not selected by default CLI mode flags)

ORCA qsub template:
- `orca-qsub.script`

CORVUS input template:
- `corvus-template.in`

CORVUS qsub templates:
- `corvus-qsub.script`
- `corvus-wrapper-qsub.script`
- `postprocess-qsub.script`

## Quick examples
```bash
python prepare-orca.py /path/to/xyz --out-dir /path/to/output
python prepare-corvus.py /path/to/orca/output --out-dir /path/to/output
python script-process-feff-output.py /path/to/corvus/run
```