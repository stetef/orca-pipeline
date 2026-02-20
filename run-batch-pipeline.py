#!/usr/bin/env python3
"""Submit an end-to-end ORCA -> CORVUS -> postprocess PBS workflow for a batch.

Workflow overview (per structure ID):
1) Prepare ORCA inputs/scripts with prepare-orca.py (always called with --dry-run)
2) Submit ORCA qsub job
3) Submit dependent CORVUS qsub job (afterok on ORCA) that:
   - runs prepare-corvus.py inside the ORCA run directory
   - fails fast if <ID>.hess is missing (prepare-corvus enforces this)
   - executes generated corvus-qsub.script inline in the same allocated PBS job

Batch-level postprocess:
4) Submit one dependent postprocess qsub job (afterok on *all* CORVUS jobs) that runs:
   - extract-orca-compute-times.py
   - process-feff-output.py --recursive
   - prepare-files-for-download.py
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


JOB_ID_RE = re.compile(r"(?P<id>\d+)(?:\.[^\s]+)?")


@dataclass
class JobRecord:
    run_id: str
    run_dir: str
    orca_script: str
    orca_job_id: str
    orca_submitted_utc: str
    corvus_wrapper_script: str
    corvus_job_id: str
    corvus_submitted_utc: str


@dataclass
class BatchState:
    created_utc: str
    input_path: str
    output_root: str
    download_destination: str
    cys: int
    his: int
    h_only: bool
    postprocess_job_id: str | None
    runs: list[JobRecord]


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _check_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")


def _run_command(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def _parse_qsub_job_id(stdout_text: str) -> str:
    for line in stdout_text.splitlines():
        line = line.strip()
        if not line:
            continue
        match = JOB_ID_RE.match(line)
        if match:
            return match.group("id")
    raise ValueError(f"Unable to parse qsub job id from output: {stdout_text!r}")


def _discover_xyz_files(path_arg: Path) -> tuple[list[Path], Path]:
    if not path_arg.exists():
        raise FileNotFoundError(f"Path not found: {path_arg}")

    if path_arg.is_file():
        if path_arg.suffix.lower() != ".xyz":
            raise ValueError(f"Expected an XYZ file, got: {path_arg}")
        return [path_arg.resolve()], path_arg.resolve().parent

    xyz_files = sorted(path_arg.glob("*.xyz"))
    if not xyz_files:
        raise ValueError(f"No .xyz files found in directory: {path_arg}")
    return [p.resolve() for p in xyz_files], path_arg.resolve()


def _run_id_from_xyz(xyz_file: Path, h_only: bool) -> str:
    base = xyz_file.name.split(".")[0].split("_")[0]
    return f"{base}-H-only" if h_only else base


def _write_corvus_wrapper_script(script_path: Path, run_dir: Path, run_id: str, prepare_corvus_py: Path) -> None:
    template_path = Path(__file__).resolve().parent / "corvus-wrapper-qsub.script"
    if not template_path.exists():
        raise FileNotFoundError(f"Missing template: {template_path}")

    script = template_path.read_text(encoding="utf-8")
    script = script.replace("[RUN_DIR]", str(run_dir))
    script = script.replace("[RUN_ID]", run_id)
    script = script.replace("[PREP_CORVUS]", str(prepare_corvus_py))

    script_path.write_text(script if script.endswith("\n") else script + "\n", encoding="utf-8")
    script_path.chmod(0o755)


def _write_postprocess_script(
    script_path: Path,
    script_dir: Path,
    output_root: Path,
    download_destination: Path,
    skip_extract: bool,
    skip_process_feff: bool,
    skip_prepare_download: bool,
) -> None:
    extract_py = script_dir / "extract-orca-compute-times.py"
    process_feff_py = script_dir / "process-feff-output.py"
    prepare_download_py = script_dir / "prepare-files-for-download.py"

    template_path = script_dir / "postprocess-qsub.script"
    if not template_path.exists():
        raise FileNotFoundError(f"Missing template: {template_path}")

    extract_cmd = (
        f"python \"{extract_py}\" \"{output_root}\" --output-dir \"{output_root}\""
        if not skip_extract
        else "true"
    )
    process_feff_cmd = (
        f"python \"{process_feff_py}\" \"{output_root}\" --recursive"
        if not skip_process_feff
        else "true"
    )
    prepare_download_cmd = (
        f"python \"{prepare_download_py}\" \"{output_root}\" -d \"{download_destination}\""
        if not skip_prepare_download
        else "true"
    )

    script = template_path.read_text(encoding="utf-8")
    script = script.replace("[BATCH_NAME]", output_root.name)
    script = script.replace("[OUTPUT_ROOT]", str(output_root))
    script = script.replace("[EXTRACT_CMD]", extract_cmd)
    script = script.replace("[PROCESS_FEFF_CMD]", process_feff_cmd)
    script = script.replace("[PREPARE_DOWNLOAD_CMD]", prepare_download_cmd)

    script_path.write_text(script if script.endswith("\n") else script + "\n", encoding="utf-8")
    script_path.chmod(0o755)


def _submit_job(script_path: Path, cwd: Path, depend_afterok: Iterable[str] | None = None) -> str:
    cmd = ["qsub"]
    if depend_afterok:
        dep_expr = "afterok:" + ":".join(str(jobid) for jobid in depend_afterok)
        cmd.extend(["-W", f"depend={dep_expr}"])
    cmd.append(script_path.name)

    result = _run_command(cmd, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(
            f"qsub failed for {script_path} (cwd={cwd})\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return _parse_qsub_job_id(result.stdout)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Submit ORCA->CORVUS pipeline with PBS dependencies, then submit a final "
            "batch postprocess job after all CORVUS jobs succeed."
        )
    )
    parser.add_argument("path", type=Path, help="XYZ directory or single XYZ file")
    parser.add_argument("--cys", type=int, required=True, help="Number of cysteine ligands")
    parser.add_argument("--his", type=int, required=True, help="Number of histidine ligands")
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Batch output directory where per-ID run dirs are created",
    )
    parser.add_argument(
        "--H",
        action="store_true",
        help="Use H-only ORCA template (propagates to prepare-orca)",
    )
    parser.add_argument(
        "--download-destination",
        type=Path,
        default=Path("downloading-station"),
        help="Destination for prepare-files-for-download.py",
    )
    parser.add_argument(
        "--skip-extract",
        action="store_true",
        help="Skip extract-orca-compute-times.py in final postprocess job",
    )
    parser.add_argument(
        "--skip-process-feff",
        action="store_true",
        help="Skip process-feff-output.py in final postprocess job",
    )
    parser.add_argument(
        "--skip-prepare-download",
        action="store_true",
        help="Skip prepare-files-for-download.py in final postprocess job",
    )
    parser.add_argument(
        "--state-file",
        type=Path,
        default=None,
        help="Optional explicit path for pipeline state JSON",
    )
    parser.add_argument(
        "--no-submit",
        action="store_true",
        help="Generate scripts and state file only; do not qsub",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.cys + args.his != 4:
        raise SystemExit(f"ERROR: cys + his must equal 4 (got {args.cys + args.his})")

    if not args.skip_process_feff:
        # Keep explicit: process-feff-output imports numpy/matplotlib/larch at runtime.
        # This pre-check catches missing python early on head/login nodes.
        _check_executable("python")

    if not args.no_submit:
        _check_executable("qsub")

    script_dir = Path(__file__).resolve().parent
    prepare_orca_py = script_dir / "prepare-orca.py"
    prepare_corvus_py = script_dir / "prepare-corvus.py"

    if not prepare_orca_py.exists() or not prepare_corvus_py.exists():
        raise SystemExit("ERROR: Missing prepare-orca.py or prepare-corvus.py next to this script")

    output_root = args.out_dir.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    download_destination = args.download_destination.expanduser().resolve()

    xyz_files, _ = _discover_xyz_files(args.path.expanduser())

    prepare_cmd = [
        "python",
        str(prepare_orca_py),
        str(args.path.expanduser()),
        "--cys",
        str(args.cys),
        "--his",
        str(args.his),
        "--out-dir",
        str(output_root),
        "--dry-run",
    ]
    if args.H:
        prepare_cmd.append("--H")

    prep_result = _run_command(prepare_cmd)
    if prep_result.returncode != 0:
        raise RuntimeError(
            "prepare-orca.py failed:\n"
            f"stdout:\n{prep_result.stdout}\n"
            f"stderr:\n{prep_result.stderr}"
        )

    records: list[JobRecord] = []
    for xyz in xyz_files:
        run_id = _run_id_from_xyz(xyz, h_only=args.H)
        run_dir = output_root / run_id
        if not run_dir.is_dir():
            raise FileNotFoundError(f"Expected run directory not found: {run_dir}")

        orca_script = run_dir / f"generated-{run_id}-orca.script"
        if not orca_script.is_file():
            matches = sorted(run_dir.glob("generated-*-orca.script"))
            if len(matches) == 1:
                orca_script = matches[0]
            else:
                raise FileNotFoundError(
                    f"Could not locate generated ORCA qsub script in {run_dir}"
                )

        corvus_wrapper = run_dir / f"generated-{run_id}-corvus-wrapper.script"
        _write_corvus_wrapper_script(corvus_wrapper, run_dir, run_id, prepare_corvus_py)

        if args.no_submit:
            orca_job_id = "NO_SUBMIT"
            orca_submitted_utc = _utc_now_iso()
            corvus_job_id = "NO_SUBMIT"
            corvus_submitted_utc = _utc_now_iso()
        else:
            orca_submitted_utc = _utc_now_iso()
            orca_job_id = _submit_job(orca_script, cwd=run_dir)
            corvus_submitted_utc = _utc_now_iso()
            corvus_job_id = _submit_job(corvus_wrapper, cwd=run_dir, depend_afterok=[orca_job_id])

        records.append(
            JobRecord(
                run_id=run_id,
                run_dir=str(run_dir),
                orca_script=str(orca_script),
                orca_job_id=orca_job_id,
                orca_submitted_utc=orca_submitted_utc,
                corvus_wrapper_script=str(corvus_wrapper),
                corvus_job_id=corvus_job_id,
                corvus_submitted_utc=corvus_submitted_utc,
            )
        )

    postprocess_script = output_root / f"generated-postprocess-{output_root.name}.script"
    _write_postprocess_script(
        postprocess_script,
        script_dir,
        output_root,
        download_destination,
        skip_extract=args.skip_extract,
        skip_process_feff=args.skip_process_feff,
        skip_prepare_download=args.skip_prepare_download,
    )

    postprocess_job_id: str | None
    if args.no_submit:
        postprocess_job_id = "NO_SUBMIT"
    else:
        corvus_ids = [rec.corvus_job_id for rec in records]
        postprocess_job_id = _submit_job(postprocess_script, cwd=output_root, depend_afterok=corvus_ids)

    state_file = (
        args.state_file.expanduser().resolve()
        if args.state_file is not None
        else output_root / f"pipeline-state-{output_root.name}.json"
    )

    state = BatchState(
        created_utc=_utc_now_iso(),
        input_path=str(args.path.expanduser().resolve()),
        output_root=str(output_root),
        download_destination=str(download_destination),
        cys=args.cys,
        his=args.his,
        h_only=bool(args.H),
        postprocess_job_id=postprocess_job_id,
        runs=records,
    )
    state_file.write_text(
        json.dumps(
            {
                **asdict(state),
                "prepare_orca_stdout": prep_result.stdout,
                "prepare_orca_stderr": prep_result.stderr,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Prepared runs: {len(records)}")
    for rec in records:
        print(
            f"  {rec.run_id}: ORCA={rec.orca_job_id}, CORVUS={rec.corvus_job_id}"
        )
    print(f"Postprocess job: {postprocess_job_id}")
    print(f"State file: {state_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
