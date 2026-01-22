#!/usr/bin/env python3
"""Clean XYZ for ORCA and extract trailing comments to a sidecar file.

Usage:
  xyz_clean_and_comments.py input.xyz [clean.xyz] [comments.txt]

- Removes anything after the 4th token per atom line.
- Preserves first two lines (natoms + comment) as-is.
- Writes sidecar comment lines as: Atom <index>: <comment>
  where <index> starts at 0 and follows atom order.
"""

from __future__ import annotations

import os
import stat
import sys
from pathlib import Path


def usage() -> None:
    print(
        "Usage: xyz_clean_and_comments.py input.xyz [clean.xyz] [comments.txt]",
        file=sys.stderr,
    )


def main() -> int:
    if len(sys.argv) < 2:
        usage()
        return 2

    input_path = Path(sys.argv[1])
    base = input_path.name.split(".")[0].split("_")[0]
    clean_path = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else input_path.with_name(f"{base}_clean.xyz")
    )
    comments_path = (
        Path(sys.argv[3])
        if len(sys.argv) > 3
        else input_path.with_name(f"{base}_comments.txt")
    )

    try:
        lines = input_path.read_text().splitlines()
    except FileNotFoundError:
        print(f"File not found: {input_path}", file=sys.stderr)
        return 1

    if len(lines) < 2:
        print("Input does not look like XYZ (missing header lines).", file=sys.stderr)
        return 1

    header = lines[:2]
    atom_lines = lines[2:]

    cleaned = [header[0], header[1]]
    comment_lines: list[str] = [header[1]]

    atom_index = 0
    for raw in atom_lines:
        if not raw.strip():
            continue

        # Split on # to separate inline comment.
        main_part, hash_part, comment = raw.partition("#")
        tokens = main_part.split()
        if len(tokens) < 4:
            print(f"Skipping invalid atom line: {raw}", file=sys.stderr)
            continue

        cleaned.append("{:<2} {:>12} {:>12} {:>12}".format(*tokens[:4]))

        if hash_part:
            comment_lines.append(f"Atom {atom_index}: # {comment.strip()}")
        else:
            comment_lines.append(f"Atom {atom_index}: (no comment)")

        atom_index += 1

    clean_path.write_text("\n".join(cleaned) + "\n")
    comments_path.write_text("\n".join(comment_lines) + "\n")

    # Ensure group read access for new files.
    for path in (clean_path, comments_path):
        mode = os.stat(path).st_mode
        os.chmod(path, mode | stat.S_IRGRP)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
