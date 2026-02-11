#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np


def _resolve_dir(path_str: str) -> Path:
    run_dir = Path(path_str).expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Directory not found: {run_dir}")
    if not run_dir.is_dir():
        raise NotADirectoryError(f"Not a directory: {run_dir}")
    return run_dir


ANGSTROM_PER_BOHR = 0.52917724899

_ATOMIC_SYMBOLS = {
    "H": 1,
    "HE": 2,
    "LI": 3,
    "BE": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "NE": 10,
    "NA": 11,
    "MG": 12,
    "AL": 13,
    "SI": 14,
    "P": 15,
    "S": 16,
    "CL": 17,
    "AR": 18,
    "K": 19,
    "CA": 20,
    "SC": 21,
    "TI": 22,
    "V": 23,
    "CR": 24,
    "MN": 25,
    "FE": 26,
    "CO": 27,
    "NI": 28,
    "CU": 29,
    "ZN": 30,
    "GA": 31,
    "GE": 32,
    "AS": 33,
    "SE": 34,
    "BR": 35,
    "KR": 36,
    "RB": 37,
    "SR": 38,
    "Y": 39,
    "ZR": 40,
    "NB": 41,
    "MO": 42,
    "TC": 43,
    "RU": 44,
    "RH": 45,
    "PD": 46,
    "AG": 47,
    "CD": 48,
    "IN": 49,
    "SN": 50,
    "SB": 51,
    "TE": 52,
    "I": 53,
    "XE": 54,
}

_ATOMIC_MASSES_AMU = {
    1: 1.00794,
    2: 4.002602,
    3: 6.941,
    4: 9.012182,
    5: 10.811,
    6: 12.0107,
    7: 14.0067,
    8: 15.9994,
    9: 18.9984032,
    10: 20.1797,
    11: 22.98976928,
    12: 24.3050,
    13: 26.9815386,
    14: 28.0855,
    15: 30.973762,
    16: 32.065,
    17: 35.453,
    18: 39.948,
    19: 39.0983,
    20: 40.078,
    21: 44.955912,
    22: 47.867,
    23: 50.9415,
    24: 51.9961,
    25: 54.938045,
    26: 55.845,
    27: 58.933195,
    28: 58.6934,
    29: 63.546,
    30: 65.38,
    31: 69.723,
    32: 72.64,
    33: 74.92160,
    34: 78.96,
    35: 79.904,
    36: 83.798,
    37: 85.4678,
    38: 87.62,
    39: 88.90585,
    40: 91.224,
    41: 92.90638,
    42: 95.96,
    43: 98.0,
    44: 101.07,
    45: 102.90550,
    46: 106.42,
    47: 107.8682,
    48: 112.411,
    49: 114.818,
    50: 118.710,
    51: 121.760,
    52: 127.60,
    53: 126.90447,
    54: 131.293,
}


def _atomic_number_from_token(token: str) -> int:
    token = token.strip()
    if not token:
        raise ValueError("Empty atom token")
    if token.isdigit():
        return int(token)
    sym = token.upper()
    if sym not in _ATOMIC_SYMBOLS:
        raise ValueError(
            f"Unknown element symbol '{token}'. Extend _ATOMIC_SYMBOLS/_ATOMIC_MASSES_AMU."
        )
    return _ATOMIC_SYMBOLS[sym]


def _atomic_mass_amu(z: int) -> float:
    if z not in _ATOMIC_MASSES_AMU:
        raise ValueError(
            f"Unknown atomic mass for Z={z}. Extend _ATOMIC_MASSES_AMU for this element."
        )
    return float(_ATOMIC_MASSES_AMU[z])


def _read_orca_hessian(filename: Path) -> tuple[np.ndarray, int]:
    with open(filename, "r") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if "$HESSIAN" in line.upper():
            start = i + 1
            break

    if start is None:
        raise ValueError("HESSIAN section not found.")

    size = int(lines[start].strip())
    natoms = int(size / 3)
    hess_start = start + 1

    hessian = np.zeros((size, size))
    row = 0
    i = hess_start

    while row < size:
        header = lines[i].split()
        ncols = len(header)
        cols = [int(x) for x in header]
        i += 1
        for r in range(size):
            parts = lines[i + r].split()
            values = list(map(float, parts[1:1 + ncols]))
            for c, val in zip(cols, values):
                hessian[r, c] = val
        if c >= size - 1:
            break
        i += size
        row += 1

    return hessian, natoms


def _read_xyz(filename: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with open(filename, "r") as f:
        raw_lines = [ln.strip() for ln in f.readlines()]

    lines = [ln for ln in raw_lines if ln]
    if not lines:
        raise ValueError("Empty XYZ file")

    atom_lines = None
    try:
        n = int(lines[0].split()[0])
        if len(lines) < 2 + n:
            raise ValueError("XYZ file too short for declared atom count")
        atom_lines = lines[2:2 + n]
    except ValueError:
        atom_lines = [ln for ln in lines if len(ln.split()) >= 4]

    atomic_numbers = []
    masses_amu = []
    coords_bohr = []

    for ln in atom_lines:
        parts = ln.split()
        if len(parts) < 4:
            continue
        z = _atomic_number_from_token(parts[0])
        x_a, y_a, z_a = map(float, parts[1:4])
        atomic_numbers.append(z)
        masses_amu.append(_atomic_mass_amu(z))
        coords_bohr.append(
            [x_a / ANGSTROM_PER_BOHR, y_a / ANGSTROM_PER_BOHR, z_a / ANGSTROM_PER_BOHR]
        )

    if not atomic_numbers:
        raise ValueError("No atom records found in XYZ")

    return (
        np.array(atomic_numbers, dtype=int),
        np.array(masses_amu, dtype=float),
        np.array(coords_bohr, dtype=float),
    )


def _print_atom_pair_blocks(hessian: np.ndarray, natoms: int, stream=None) -> None:
    if stream is None:
        stream = sys.stdout

    for i in range(natoms):
        for j in range(natoms):
            print(i + 1, j + 1, file=stream)
            ii = 0
            while ii < 3:
                print(
                    " ".join(
                        f"{h:12.6e}" for h in hessian[3 * i + ii, 3 * j: 3 * j + 3]
                    ),
                    file=stream,
                )
                ii += 1


def _write_dym_file(
    dymout_filename: Path,
    atomic_numbers: np.ndarray,
    masses_amu: np.ndarray,
    coords_bohr: np.ndarray,
    hess: np.ndarray,
) -> None:
    """
    Write dynamic file from ORCA hessian and xyz files. Formatted as follows:
    Line 1 - dym_Type: Dynamical matrix file type (integer)
        This value is for future use. Set to 1 for now.
    Line 2 - nAt: Number of atoms (integer)
        Number of atoms in the system. (61 in this case)
    Lines 2..2+nAt - Atomic numbers (integer)
        Atomic numbers of atoms in the system.
    Lines 2+nAt+1..2+2*nAt - Atomic masses (real, in AMU)
        Atomic masses of the atoms in the system.
    Lines 2+2*nAt+1..2+3*nAt - Atomic coordinates (real, in Bohr)
        Cartesian coordinates ("x y z") of the atoms in the system.
        Directly from xyz file, but with conversion from angstrom to Bohr
    Lines 2+3*nAt+1..End - Dynamical matrix in atom pair block format (integer and
        real, see below,in atomic units)
        From `print_atom_pair_blocks` function
    """
    natoms = int(len(atomic_numbers))
    if hess.shape != (3 * natoms, 3 * natoms):
        raise ValueError(
            f"Hessian shape {hess.shape} does not match 3N x 3N with N={natoms}"
        )

    with open(dymout_filename, "w") as f:
        f.write("1\n")
        f.write(f"{natoms}\n")
        for z in atomic_numbers:
            f.write(f"{int(z)}\n")
        for m in masses_amu:
            f.write(f"{float(m):.10f}\n")
        for (x, y, z) in coords_bohr:
            f.write(f"{x:.10f} {y:.10f} {z:.10f}\n")
        _print_atom_pair_blocks(hess, natoms, stream=f)


def _copy_and_replace_qsub(template_path: Path, dest_path: Path, run_dir: Path, run_id: str) -> None:
    content = template_path.read_text(encoding="utf-8")
    content = content.replace("[directory]", f"{run_dir}/")
    content = content.replace("[ID]", run_id)
    dest_path.write_text(content, encoding="utf-8")


def _copy_and_replace_corvus(
    template_path: Path, dest_path: Path, run_dir: Path, run_id: str
) -> None:
    content = template_path.read_text(encoding="utf-8")
    content = content.replace("[DIRECTORY]", f"{run_dir}/")
    content = content.replace("[ID]", run_id)
    dest_path.write_text(content, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare Corvus run directory (dym + templates)."
    )
    parser.add_argument(
        "run_dir",
        help="Path to the run directory; its name is used as the run ID.",
    )
    args = parser.parse_args()

    run_dir = _resolve_dir(args.run_dir)
    run_id = run_dir.name

    script_dir = Path(__file__).resolve().parent
    corvus_template_path = script_dir / "corvus-template.in"
    qsub_template_path = script_dir / "corvus-qsub.script"

    if not corvus_template_path.exists():
        raise FileNotFoundError(f"Missing corvus-template.in at {corvus_template_path}")
    if not qsub_template_path.exists():
        raise FileNotFoundError(f"Missing corvus-qsub.script at {qsub_template_path}")

    hess_path = run_dir / f"{run_id}.hess"
    xyz_path = run_dir / f"{run_id}.xyz"
    dym_path = run_dir / f"{run_id}.dym"

    if not hess_path.exists():
        raise FileNotFoundError(f"Missing Hessian file: {hess_path}")
    if not xyz_path.exists():
        raise FileNotFoundError(f"Missing XYZ file: {xyz_path}")

    hess, natoms_hess = _read_orca_hessian(hess_path)
    atomic_numbers, masses_amu, coords_bohr = _read_xyz(xyz_path)
    if natoms_hess != len(atomic_numbers):
        raise ValueError(
            f"Atom count mismatch: Hessian implies {natoms_hess} atoms, XYZ has {len(atomic_numbers)}"
        )

    _write_dym_file(dym_path, atomic_numbers, masses_amu, coords_bohr, hess)

    zn_indices = np.where(atomic_numbers == _ATOMIC_SYMBOLS["ZN"])[0]
    if zn_indices.size == 0:
        raise ValueError("No Zn atoms found in XYZ; cannot center DYM on absorber.")

    zn_index_1based = int(zn_indices[0]) + 1
    feff_dym_path = run_dir / f"corvus-{run_id}.dym"
    dym2feffinp_cmd = [
        "/opt/feff10/bin/MPI/dym2feffinp",
        "--c",
        str(zn_index_1based),
        "--d",
        str(feff_dym_path.name),
        str(dym_path.name),
    ]
    subprocess.run(dym2feffinp_cmd, check=True, cwd=run_dir)

    corvus_in_dest = run_dir / f"corvus-{run_id}.in"
    _copy_and_replace_corvus(corvus_template_path, corvus_in_dest, run_dir, run_id)

    qsub_dest = run_dir / "corvus-qsub.script"
    _copy_and_replace_qsub(qsub_template_path, qsub_dest, run_dir, run_id)

    print(f"Job can be submitted with: qsub {run_dir}/corvus-qsub.script")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
