#!/usr/bin/env python3
"""
Scans a family directory for imaginary frequency warnings in corvus-*.out files.
Writes results to a CSV in the parent directory.
"""

import os
import glob
import csv
import sys
from pathlib import Path


def count_imaginary_frequencies(family_dir: str) -> dict[str, int]:
    family_path = Path(family_dir).resolve()
    results = {}

    # Iterate over every subdirectory in the family dir
    for entry in sorted(family_path.iterdir()):
        if not entry.is_dir():
            continue

        subdir_name = entry.name

        # Find the working-* directory inside
        working_dirs = list(entry.glob(f"working-{subdir_name}"))
        if not working_dirs:
            # Try a looser glob in case naming differs slightly
            working_dirs = list(entry.glob("working-*"))

        if not working_dirs:
            print(f"  [SKIP] No working-* dir found in: {entry}")
            results[subdir_name] = None
            continue

        working_dir = working_dirs[0]

        # Find the corvus-*.out file
        out_files = list(working_dir.glob("corvus-*.out"))
        if not out_files:
            print(f"  [SKIP] No corvus-*.out file found in: {working_dir}")
            results[subdir_name] = None
            continue

        out_file = out_files[0]

        # Count occurrences of the target string
        target = "Found imaginary frequency with large weight"
        count = 0
        try:
            with open(out_file, "r", errors="replace") as f:
                for line in f:
                    count += line.count(target)
        except OSError as e:
            print(f"  [ERROR] Could not read {out_file}: {e}")
            results[subdir_name] = None
            continue

        print(f"  {subdir_name}: {count}")
        results[subdir_name] = count

    return results


def write_csv(family_dir: str, results: dict[str, int]) -> Path:
    family_path = Path(family_dir).resolve()
    csv_path = family_path / "imaginary_frequencies.csv"

    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["cluster", "imaginary_freq_count"])
        for name, count in results.items():
            writer.writerow([name, count if count is not None else "N/A"])

    return csv_path


def main():
    if len(sys.argv) < 2:
        print("Usage: python count-imag-freq.py <family_dir>")
        sys.exit(1)

    family_dir = sys.argv[1]

    if not os.path.isdir(family_dir):
        print(f"Error: '{family_dir}' is not a valid directory.")
        sys.exit(1)

    print(f"Scanning: {family_dir}\n")
    results = count_imaginary_frequencies(family_dir)

    csv_path = write_csv(family_dir, results)
    print(f"\nResults written to: {csv_path}")


if __name__ == "__main__":
    main()