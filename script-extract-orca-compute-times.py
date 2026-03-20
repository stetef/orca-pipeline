#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


RUNTIME_LINE_RE = re.compile(r"^TOTAL RUN TIME:\s*(.+?)\s*$")
FINAL_GIBBS_RE = re.compile(
    r"^\s*Final Gibbs free energy\s*\.\.\.\s*([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+Eh\s*$"
)
TERMINATION_MARKER = "****ORCA TERMINATED NORMALLY****"


def extract_runtime_from_log(log_path: Path) -> str | None:
    """Return the last TOTAL RUN TIME value from a normally terminated ORCA log."""
    last_runtime = None
    saw_termination_marker = False

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if TERMINATION_MARKER in line:
                saw_termination_marker = True
                continue

            match = RUNTIME_LINE_RE.match(line)
            if saw_termination_marker and match:
                last_runtime = match.group(1)
                saw_termination_marker = False

    return last_runtime


def extract_final_gibbs_from_log(log_path: Path) -> str | None:
    """Return the last 'Final Gibbs free energy' value (in Eh) from an ORCA log."""
    last_final_gibbs = None

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            match = FINAL_GIBBS_RE.match(line)
            if match:
                last_final_gibbs = match.group(1)

    return last_final_gibbs


def find_orca_logs(parent_dir: Path) -> list[Path]:
    """Find all files matching *-orca.log anywhere below parent_dir."""
    return sorted(path for path in parent_dir.rglob("*-orca.log") if path.is_file())


def write_csv(output_path: Path, rows: list[tuple[str, str, str | None]]) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "total_run_time", "final_gibbs_free_energy_eh"])
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract ORCA TOTAL RUN TIME values from *-orca.log files under a parent directory "
            "and write them to a CSV."
        )
    )
    parser.add_argument("parent_dir", type=Path, help="Parent directory containing ORCA run directories")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory where the CSV is written (default: current working directory)",
    )

    args = parser.parse_args()
    parent_dir = args.parent_dir.resolve()

    if not parent_dir.is_dir():
        raise SystemExit(f"Error: '{parent_dir}' is not a directory.")

    log_files = find_orca_logs(parent_dir)
    rows: list[tuple[str, str, str | None]] = []

    for log_file in log_files:
        runtime = extract_runtime_from_log(log_file)
        if runtime is None:
            continue

        final_gibbs = extract_final_gibbs_from_log(log_file)

        name = log_file.name.removesuffix("-orca.log")
        rows.append((name, runtime, final_gibbs))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / f"{parent_dir.name}-orca-compute-times.csv"
    write_csv(output_path, rows)

    print(f"Scanned {len(log_files)} log file(s).")
    print(f"Wrote {len(rows)} row(s) to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())