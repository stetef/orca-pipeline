#!/usr/bin/env python3
"""
Process XYZ files for ORCA calculations with specified ligand composition.
"""

import argparse
import os
import sys
import shutil
import subprocess
from pathlib import Path


def get_charge_multiplicity(cys, his):
    """Get charge and multiplicity based on cysteine and histidine counts."""
    lookup = {
        (4, 0): (-2, 1),  # (cys, his) : (charge, mult)
        (3, 1): (-1, 2),
        (2, 2): (0, 1),
        (1, 3): (1, 2),
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


def process_xyz_file(xyz_file, cys, his, template_dir, output_root, dry_run):
    """Process a single XYZ file."""
    # Extract ID from filename
    filename = os.path.basename(xyz_file)
    # First split by ".", take first part, then split by "_", take first part
    id_name = filename.split(".")[0].split("_")[0]
    
    print(f"\nProcessing {filename} -> ID: {id_name}")
    
    # Create directory for this ID under output root
    id_dir = output_root / id_name
    id_dir.mkdir(exist_ok=True)
    print(f"  Created directory: {id_dir}")
    
    # Copy XYZ file to the directory (leave original in place)
    dest_xyz = id_dir / filename
    shutil.copy2(xyz_file, dest_xyz)
    os.chmod(dest_xyz, 0o644)
    print(f"  Copied {filename} to {id_dir}/")
    
    # Run xyz_clean_and_comments.py on the XYZ file
    xyz_script = template_dir / "xyz_clean_and_comments.py"
    if xyz_script.exists():
        print(f"  Running xyz_clean_and_comments.py...")
        result = subprocess.run(
            ["python3", str(xyz_script), str(dest_xyz)],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            print(f"  Warning: xyz_clean_and_comments.py returned error: {result.stderr}")
    else:
        print(f"  Warning: xyz_clean_and_comments.py not found at {xyz_script}")
    
    # Get charge and multiplicity
    charge, multiplicity = get_charge_multiplicity(cys, his)
    
    # Copy and modify template
    template_file = template_dir / "orca_template.in"
    output_file = id_dir / f"{id_name}.in"
    
    if not template_file.exists():
        print(f"  Error: Template file not found: {template_file}")
        return
    
    with open(template_file, 'r') as f:
        template_content = f.read()
    
    # Replace placeholders
    template_content = template_content.replace('[CHARGE]', str(charge))
    template_content = template_content.replace('[MULTIPLICITY]', str(multiplicity))
    template_content = template_content.replace('[PDB_ID]', id_name)
    template_content = template_content.replace('[ID_DIR]', str(id_dir))
    
    # Extract CA atoms from comments file
    comments_file = id_dir / f"{id_name}_comments.txt"
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
    print(f"    CA atoms to freeze: {len(ca_atoms)}")
    
    # Copy run_job.sh script
    run_script = template_dir / "run_job.sh"
    if run_script.exists():
        dest_script = id_dir / "run_job.sh"
        shutil.copy2(run_script, dest_script)
        os.chmod(dest_script, 0o755)
        print(f"  Copied run_job.sh to {id_dir}/")
        
        # Run job (optionally dry-run)
        run_args = ["./run_job.sh"]
        if dry_run:
            run_args.append("--dry-run")
        run_args.append(f"{id_name}.in")

        run_label = "dry-run" if dry_run else "run"
        print(f"  Running {run_label}...")
        result = subprocess.run(
            run_args,
            cwd=id_dir,
            capture_output=True,
            text=True
        )
        print(f"  {run_label.capitalize()} output:\n{result.stdout}")
        if result.stderr:
            print(f"  {run_label.capitalize()} stderr:\n{result.stderr}")
    else:
        print(f"  Warning: run_job.sh not found at {run_script}")


def main():
    parser = argparse.ArgumentParser(
        description='Process XYZ files for ORCA calculations with specified ligand composition.'
    )
    parser.add_argument('path', type=str, help='Directory containing XYZ files or a single XYZ file')
    parser.add_argument('--cys', type=int, required=True, help='Number of cysteines')
    parser.add_argument('--his', type=int, required=True, help='Number of histidines')
    parser.add_argument('--out-dir', type=str, default=None, help='Output directory for ID folders')
    parser.add_argument('-n', '--dry-run', action='store_true', help='Run run_job.sh with --dry-run')
    
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
    output_root = Path(args.out_dir).resolve() if args.out_dir else Path.cwd().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    # Find all XYZ files in the specified directory
    target_path = Path(args.path)
    if not target_path.exists():
        print(f"ERROR: Path not found: {target_path}")
        sys.exit(1)

    if target_path.is_file():
        if target_path.suffix.lower() != ".xyz":
            print(f"ERROR: File is not an XYZ file: {target_path}")
            sys.exit(1)
        xyz_files = [target_path]
    else:
        xyz_files = list(target_path.glob("*.xyz"))

        if not xyz_files:
            print(f"No XYZ files found in {target_path}")
            sys.exit(0)
    
    print(f"\nFound {len(xyz_files)} XYZ file(s) to process")
    print(f"Output root: {output_root}")
    
    for xyz_file in xyz_files:
        process_xyz_file(xyz_file, args.cys, args.his, template_dir, output_root, args.dry_run)
    
    print("\nProcessing complete!")


if __name__ == '__main__':
    main()
