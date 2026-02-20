#!/usr/bin/env python3
"""
Find and copy output* directories from each child directory under a parent directory.

Example:
    python prepare-files-for-download.py /path/to/parent -d ./downloaded_outputs
"""

import argparse
import shutil
import sys
from pathlib import Path


def find_output_dirs(parent_dir: Path):
    """Yield (child_dir, output_dir) for output* directories in each child of parent_dir."""
    for child_dir in sorted(parent_dir.iterdir()):
        if not child_dir.is_dir():
            continue

        for output_dir in sorted(child_dir.glob("output*")):
            if output_dir.is_dir():
                yield child_dir, output_dir


def copy_output_dirs(parent_dir: Path, destination_dir: Path, dry_run: bool = False):
    """Copy discovered output* directories into destination_dir.

    Destination naming:
        <child_dir_name>
    """
    destination_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    skipped = 0

    for child_dir, output_dir in find_output_dirs(parent_dir):
        target_name = child_dir.name
        target_dir = destination_dir / target_name

        if target_dir.exists():
            print(f"SKIP (exists): {target_dir}")
            skipped += 1
            continue

        print(f"COPY: {output_dir} -> {target_dir}")
        if not dry_run:
            shutil.copytree(output_dir, target_dir)
        copied += 1

    return copied, skipped


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Search each child directory in parent_dir for output* directories "
            "and copy them to a destination directory."
        )
    )
    parser.add_argument(
        "parent_dir",
        type=Path,
        help="Parent directory containing child directories to scan.",
    )
    parser.add_argument(
        "-d",
        "--destination",
        type=Path,
        default=Path("downloading-station"),
        help="Destination directory for copied output directories (default: ./downloading-station).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be copied without copying files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    parent_dir = args.parent_dir.expanduser().resolve()
    destination_dir = args.destination.expanduser().resolve()

    if not parent_dir.exists() or not parent_dir.is_dir():
        print(f"Error: parent_dir does not exist or is not a directory: {parent_dir}", file=sys.stderr)
        return 1

    copied, skipped = copy_output_dirs(parent_dir, destination_dir, dry_run=args.dry_run)
    print(f"Done. Copied: {copied}, Skipped: {skipped}")
    print(f"Destination: {destination_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
