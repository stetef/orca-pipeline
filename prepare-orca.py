#!/usr/bin/env python3
"""
Process XYZ files for ORCA calculations with specified ligand composition.
"""

import argparse
import os
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path


def get_charge_multiplicity(cys, his):
    """
    Get charge and multiplicity based on cysteine and histidine counts.
    His = 54e, Cys = 34e, Zn = 28e, Water = 10e -> all even
    Thus no matter what, 2S + 1 = odd = 1
    """
    lookup = {
        (4, 0): (-2, 1),  # (cys, his) : (charge, mult)
        (3, 1): (-1, 1),
        (2, 2): (0, 1),
        (1, 3): (1, 1),
        (0, 4): (2, 1),
    }
    return lookup.get((cys, his), (None, None))


def extract_ca_atoms(comments_file):
    """Extract CA atom numbers from comments file."""
    ca_atoms = []
    
    if not os.path.exists(comments_file):
        print(f"Warning: Comments file not found: {comments_file}")
        return ca_atoms
    
    with open(comments_file, 'r') as f:
        for line in f:
            if 'ATOM=CA' in line:
                # Extract atom number from "Atom 18: " format
                parts = line.split()
                if len(parts) >= 2 and parts[0] == 'Atom':
                    atom_num = parts[1].rstrip(':')
                    ca_atoms.append(atom_num)
    
    return ca_atoms


def clean_xyz_and_comments(input_path, clean_path=None, comments_path=None):
    """Clean XYZ and extract trailing comments to a sidecar file."""
    input_path = Path(input_path)
    base = input_path.name.split(".")[0].split("_")[0]
    clean_path = (
        Path(clean_path)
        if clean_path is not None
        else input_path.with_name(f"{base}_clean.xyz")
    )
    comments_path = (
        Path(comments_path)
        if comments_path is not None
        else input_path.with_name(f"{base}_comments.txt")
    )

    try:
        lines = input_path.read_text().splitlines()
    except FileNotFoundError:
        print(f"File not found: {input_path}", file=sys.stderr)
        return None, None

    if len(lines) < 2:
        print("Input does not look like XYZ (missing header lines).", file=sys.stderr)
        return None, None

    header = lines[:2]
    atom_lines = lines[2:]

    cleaned = [header[0], header[1]]
    comment_lines = [header[1]]

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

    return clean_path, comments_path


def extract_nprocs(orcar_input_file):
    """Extract number of processors from ORCA input (PAL or %pal nprocs)."""
    try:
        lines = Path(orcar_input_file).read_text().splitlines()
    except FileNotFoundError:
        return 1

    pal_pattern = re.compile(r"\bPAL\s*([0-9]+)", re.IGNORECASE)
    nprocs_pattern = re.compile(r"\bnprocs\s*([0-9]+)", re.IGNORECASE)

    in_pal_block = False
    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith("!"):
            match = pal_pattern.search(line)
            if match:
                return int(match.group(1))

        if line.lower().startswith("%pal"):
            # Handle both styles:
            #   %pal nprocs 16
            #   %pal
            #     nprocs 16
            match = nprocs_pattern.search(line)
            if match:
                return int(match.group(1))
            in_pal_block = True
            continue

        if in_pal_block:
            match = nprocs_pattern.search(line)
            if match:
                return int(match.group(1))
            if line.lower() == "end":
                in_pal_block = False

    return 1


def generate_orca_qsub_script(template_dir, id_dir, input_filename, basename, nprocs):
    """Generate a PBS qsub script from template."""
    qsub_template = template_dir / "orca-qsub.script"
    if not qsub_template.exists():
        print(f"  Error: Qsub template not found: {qsub_template}")
        return None

    qsub_content = qsub_template.read_text()
    qsub_content = qsub_content.replace("[NODES]", str(nprocs or 1))
    qsub_content = qsub_content.replace("[BASENAME]", basename)
    qsub_content = qsub_content.replace("[INPUT_FILE]", input_filename)

    generated_qsub = id_dir / f"generated-{basename}-orca.script"
    generated_qsub.write_text(qsub_content if qsub_content.endswith("\n") else qsub_content + "\n")
    os.chmod(generated_qsub, 0o755)
    return generated_qsub


def process_xyz_file(xyz_file, cys, his, template_dir, output_root, dry_run, use_h_only):
    """Process a single XYZ file."""
    # Extract ID from filename
    filename = os.path.basename(xyz_file)
    # First split by ".", take first part, then split by "_", take first part
    id_name = filename.split(".")[0].split("_")[0]
    
    print(f"\nProcessing {filename} -> ID: {id_name}")

    output_base = f"{id_name}-H-only" if use_h_only else id_name
    
    # Create directory for this ID under output root
    id_dir = output_root / output_base
    id_dir.mkdir(exist_ok=True)
    print(f"  Created directory: {id_dir}")
    
    # Copy XYZ file to the directory (leave original in place)
    dest_xyz = id_dir / filename
    shutil.copy2(xyz_file, dest_xyz)
    os.chmod(dest_xyz, 0o644)
    print(f"  Copied {filename} to {id_dir}/")
    
    # Clean XYZ and write comments sidecar
    print("  Cleaning XYZ and extracting comments...")
    clean_path = id_dir / f"{output_base}_clean.xyz"
    comments_path = id_dir / f"{output_base}_comments.txt"
    clean_result, comments_result = clean_xyz_and_comments(
        dest_xyz,
        clean_path=clean_path,
        comments_path=comments_path,
    )
    if clean_result is None or comments_result is None:
        print("  Warning: Failed to clean XYZ or write comments file")
    
    # Get charge and multiplicity
    charge, multiplicity = get_charge_multiplicity(cys, his)
    
    # Copy and modify template
    template_file = template_dir / ("orca-template-h-only.in" if use_h_only else "orca-template.in")
    output_file = id_dir / f"{output_base}.in"
    
    if not template_file.exists():
        print(f"  Error: Template file not found: {template_file}")
        return
    
    with open(template_file, 'r') as f:
        template_content = f.read()
    
    # Replace placeholders
    template_content = template_content.replace('[CHARGE]', str(charge))
    template_content = template_content.replace('[MULTIPLICITY]', str(multiplicity))
    template_content = template_content.replace('[PDB_ID]', output_base)
    template_content = template_content.replace('[ID_DIR]', str(id_dir))
    
    ca_atoms = []
    if not use_h_only:
        # Extract CA atoms from comments file
        comments_file = id_dir / f"{output_base}_comments.txt"
        ca_atoms = extract_ca_atoms(comments_file)

        # Build CA constraint lines based on [CA_ATOM] placeholder
        ca_constraints = []
        for atom_num in ca_atoms:
            constraint_line = f"  {{C {atom_num} C}}  # Freeze all coordinates (X, Y, Z) of atom {atom_num}"
            ca_constraints.append(constraint_line)

        # Replace the single [CA_ATOM] line with multiple lines (one per CA atom)
        if '[CA_ATOM]' in template_content:
            if ca_constraints:
                template_lines = template_content.splitlines()
                replaced_lines = []
                replaced = False
                for line in template_lines:
                    if '[CA_ATOM]' in line and not replaced:
                        replaced_lines.extend(ca_constraints)
                        replaced = True
                    else:
                        replaced_lines.append(line)
                template_content = "\n".join(replaced_lines)
            else:
                print("  Warning: No CA atoms found to replace [CA_ATOM]")
        elif ca_constraints:
            print("  Warning: No [CA_ATOM] placeholder found in template")
            print(f"  Found {len(ca_constraints)} CA atoms to constrain")
    
    # Ensure final newline (some ORCA versions can misread the last line without it)
    if not template_content.endswith("\n"):
        template_content += "\n"

    # Write modified template
    with open(output_file, 'w') as f:
        f.write(template_content)
    os.chmod(output_file, 0o644)
    print(f"  Created input file: {output_file}")
    print(f"    CHARGE={charge}, MULTIPLICITY={multiplicity}")
    if use_h_only:
        print("    Only optimizing hydrogen atoms")
    else:
        print(f"    CA atoms to freeze: {len(ca_atoms)}")
    
    # Generate and (optionally) submit qsub script
    nprocs = extract_nprocs(output_file)
    generated_qsub = generate_orca_qsub_script(
        template_dir,
        id_dir,
        f"{output_base}.in",
        output_base,
        nprocs,
    )
    if generated_qsub is None:
        return

    print(f"  Generated qsub script: {generated_qsub.name}")
    if dry_run:
        print(f"  Dry run: generated {generated_qsub.name} (qsub skipped)")
    else:
        qsub_cmd = [
            "qsub",
            "-j",
            "oe",
            "-o",
            f"{generated_qsub.name}.out",
            generated_qsub.name,
        ]
        print("  Submitting with qsub...")
        result = subprocess.run(
            qsub_cmd,
            cwd=id_dir,
            capture_output=True,
            text=True,
        )
        if result.stdout:
            print(f"  qsub output:\n{result.stdout}")
        if result.stderr:
            print(f"  qsub stderr:\n{result.stderr}")
        if result.returncode != 0:
            print(f"  Warning: qsub submission failed (exit code {result.returncode})")


def main():
    parser = argparse.ArgumentParser(
        description='Process XYZ files for ORCA calculations with specified ligand composition.'
    )
    parser.add_argument('path', type=str, help='Directory containing XYZ files or a single XYZ file')
    parser.add_argument('--cys', type=int, required=True, help='Number of cysteines')
    parser.add_argument('--his', type=int, required=True, help='Number of histidines')
    parser.add_argument('--out-dir', type=str, default=None, help='Output directory for ID folders')
    parser.add_argument('--H', action='store_true', help='Use orca-template-h-only.in instead of orca-template.in')
    parser.add_argument('-n', '--dry-run', action='store_true', help='Generate qsub script but skip qsub submission')
    
    args = parser.parse_args()
    
    # Validate cys + his = 4
    if args.cys + args.his != 4:
        print(f"ERROR: Number of cysteines ({args.cys}) and histidines ({args.his}) must sum to 4!")
        print(f"Current sum: {args.cys + args.his}")
        sys.exit(1)
    
    # Validate charge/multiplicity combination exists
    charge, multiplicity = get_charge_multiplicity(args.cys, args.his)
    if charge is None:
        print(f"ERROR: Invalid combination of cys={args.cys} and his={args.his}")
        sys.exit(1)
    
    print(f"Configuration: {args.cys} Cys + {args.his} His")
    print(f"Charge: {charge}, Multiplicity: {multiplicity}")
    
    # Get template directory (assume script is in template directory)
    template_dir = Path(__file__).parent.absolute()
    
    # Determine output root (default to directory where the script is invoked)
    if args.out_dir:
        raw_out_dir = Path(args.out_dir)
        if not raw_out_dir.is_absolute() and str(raw_out_dir).startswith("home/"):
            fixed_out_dir = Path("/") / raw_out_dir
            print(
                "Warning: --out-dir looks like an absolute path missing a leading '/'. "
                f"Using: {fixed_out_dir}"
            )
            output_root = fixed_out_dir.resolve()
        else:
            output_root = raw_out_dir.resolve()
            if not raw_out_dir.is_absolute():
                print(
                    "Warning: --out-dir is relative. "
                    f"It resolves to: {output_root}"
                )
    else:
        output_root = Path.cwd().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    # Find all XYZ files in the specified directory
    target_path = Path(args.path)
    target_path_was_absolute = target_path.is_absolute()
    if not target_path.exists():
        print(f"ERROR: Path not found: {target_path}")
        sys.exit(1)

    target_path_resolved = target_path.resolve()

    if target_path.is_file():
        if target_path.suffix.lower() != ".xyz":
            print(f"ERROR: File is not an XYZ file: {target_path}")
            sys.exit(1)
        xyz_files = [target_path]
        input_base_dir = target_path_resolved.parent
    else:
        xyz_files = list(target_path.glob("*.xyz"))
        input_base_dir = target_path_resolved

        if not xyz_files:
            print(f"No XYZ files found in {target_path}")
            sys.exit(0)

    if input_base_dir != output_root:
        output_was_absolute = Path(args.out_dir).is_absolute() if args.out_dir else True
        print(
            "Warning: Input XYZ location and output directory differ: "
            f"input={input_base_dir} ({'absolute' if target_path_was_absolute else 'relative'} path), "
            f"output={output_root} ({'absolute' if output_was_absolute else 'relative'} path)"
        )
    
    print(f"\nFound {len(xyz_files)} XYZ file(s) to process")
    print(f"Output root: {output_root}")
    
    for xyz_file in xyz_files:
        process_xyz_file(
            xyz_file,
            args.cys,
            args.his,
            template_dir,
            output_root,
            args.dry_run,
            args.H,
        )
    
    print("\nProcessing complete!")


if __name__ == '__main__':
    main()
