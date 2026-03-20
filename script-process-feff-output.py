#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 16,
    }
)


def load_feff_table(path: Path):
    data = np.genfromtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 6:
        raise ValueError(f"Expected at least 6 columns in {path}, got {data.shape[1]}")
    omega = data[:, 0]
    energy = data[:, 1]
    k = data[:, 2]
    mu = data[:, 3]
    mu0 = data[:, 4]
    chi = data[:, 5]
    return omega, energy, k, mu, mu0, chi


def load_chi_dat(path: Path):
    data = np.genfromtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"Expected at least 2 columns in {path}, got {data.shape[1]}")
    k = data[:, 0]
    chi = data[:, 1]
    return k, chi


def parse_xmu_nlegs2_entries(path: Path):
    in_table = False
    entries = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("#"):
                continue
            body = line.lstrip("#").strip()
            if "file" in body and "nlegs" in body and "reff" in body:
                in_table = True
                continue
            if not in_table or not body:
                continue

            fields = body.split()
            if len(fields) < 6:
                continue
            if not fields[0].isdigit():
                continue

            try:
                sig2_tot = float(fields[1])
                deg = int(float(fields[3]))
                nlegs = int(float(fields[-2]))
                reff = float(fields[-1])
            except ValueError:
                continue
            if nlegs == 2:
                entries.append({"reff": reff, "sig2_tot": sig2_tot, "deg": deg})

    if not entries:
        raise ValueError(f"No nlegs=2 entries found in {path}")
    return entries


def parse_feff_atoms(path: Path):
    atoms: List[Dict[str, float]] = []
    in_atoms = False

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            upper = stripped.upper()
            if upper.startswith("ATOMS"):
                in_atoms = True
                continue
            if not in_atoms:
                continue
            if upper == "END" or upper.startswith("END "):
                break

            fields = stripped.split()
            if len(fields) < 7:
                continue
            try:
                x = float(fields[0])
                y = float(fields[1])
                z = float(fields[2])
                atom_type = int(fields[3])
                symbol = fields[4]
                listed_distance = float(fields[5])
                atom_index = int(float(fields[6]))
            except ValueError:
                continue

            atoms.append(
                {
                    "index": atom_index,
                    "symbol": symbol,
                    "type": atom_type,
                    "x": x,
                    "y": y,
                    "z": z,
                    "distance": listed_distance,
                }
            )

    if not atoms:
        raise ValueError(f"No ATOMS block parsed in {path}")
    return atoms


def parse_xyz_atoms(path: Path):
    atoms: List[Dict[str, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        raw_lines = [line.strip() for line in handle if line.strip()]

    if not raw_lines:
        raise ValueError(f"Empty XYZ file: {path}")

    try:
        natoms = int(raw_lines[0].split()[0])
        atom_lines = raw_lines[2 : 2 + natoms]
    except (ValueError, IndexError):
        atom_lines = [line for line in raw_lines if len(line.split()) >= 4]

    for fields_line in atom_lines:
        fields = fields_line.split()
        if len(fields) < 4:
            continue
        symbol = fields[0]
        try:
            x = float(fields[1])
            y = float(fields[2])
            z = float(fields[3])
        except ValueError:
            continue
        atoms.append({"symbol": symbol, "x": x, "y": y, "z": z})

    if not atoms:
        raise ValueError(f"No atom coordinates parsed in XYZ file: {path}")
    return atoms


def parse_ca_indices_from_comments(path: Path):
    pattern = re.compile(r"^Atom\s+(\d+):.*\bATOM=CA\b", re.IGNORECASE)
    indices: List[int] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            match = pattern.match(line)
            if match:
                indices.append(int(match.group(1)))

    if not indices:
        raise ValueError(f"No ATOM=CA entries found in {path}")
    return indices


def find_first_zn_index(xyz_atoms):
    for i, atom in enumerate(xyz_atoms):
        if atom["symbol"].strip().upper() == "ZN":
            return i
    raise ValueError("No Zn atom found in XYZ file; cannot apply Zn-centering transform")


def find_latest_xyz_file(run_root: Path) -> Path:
    xyz_files = [path for path in run_root.glob("*.xyz") if path.is_file()]
    if not xyz_files:
        raise FileNotFoundError(f"No .xyz files found in {run_root}")

    # Prefer single-geometry xyz files over trajectory dumps (e.g. *_trj.xyz).
    # Trajectories often contain many frames, and parse_xyz_atoms reads only the first frame,
    # which may not match the FEFF geometry used for ATOMS matching.
    non_trj = [path for path in xyz_files if not path.stem.lower().endswith("_trj")]
    pool = non_trj if non_trj else xyz_files
    return max(pool, key=lambda path: path.stat().st_mtime)


def resolve_comments_file(run_root: Path, xyz_path: Path) -> Path:
    preferred = run_root / f"{xyz_path.stem}_comments.txt"
    if preferred.is_file():
        return preferred

    # Collect candidates from the run root first, then wider fallbacks.
    candidates = [path for path in run_root.glob("*_comments.txt") if path.is_file()]
    if not candidates:
        candidates = [path for path in run_root.rglob("*_comments.txt") if path.is_file()]
    if not candidates and run_root.parent.is_dir():
        candidates = [path for path in run_root.parent.glob("*_comments.txt") if path.is_file()]

    if not candidates:
        raise FileNotFoundError(
            "No *_comments.txt file found near "
            f"{run_root}; cannot determine CA atom indices"
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)


def match_ca_indices_to_feff_atoms(ca_indices, xyz_atoms, feff_atoms, tolerance=2e-3):
    zn_index = find_first_zn_index(xyz_atoms)
    zn_coords = np.array(
        [xyz_atoms[zn_index]["x"], xyz_atoms[zn_index]["y"], xyz_atoms[zn_index]["z"]],
        dtype=float,
    )

    matched = []
    used_feff_indices = set()
    for ca_idx in ca_indices:
        if ca_idx < 0 or ca_idx >= len(xyz_atoms):
            raise IndexError(
                f"CA atom index {ca_idx} from comments is out of range for XYZ with "
                f"{len(xyz_atoms)} atoms"
            )

        src_atom = xyz_atoms[ca_idx]
        src_symbol = src_atom["symbol"].strip().upper()
        src_coords = np.array([src_atom["x"], src_atom["y"], src_atom["z"]], dtype=float) - zn_coords

        symbol_candidates = [
            atom
            for atom in feff_atoms
            if atom["index"] not in used_feff_indices
            and atom["symbol"].strip().upper() == src_symbol
        ]
        candidates = symbol_candidates or [
            atom for atom in feff_atoms if atom["index"] not in used_feff_indices
        ]

        if not candidates:
            raise ValueError("No remaining FEFF ATOMS candidates for CA matching")

        best_atom = min(
            candidates,
            key=lambda atom: np.linalg.norm(
                src_coords - np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
            ),
        )
        best_err = float(
            np.linalg.norm(
                src_coords
                - np.array([best_atom["x"], best_atom["y"], best_atom["z"]], dtype=float)
            )
        )
        if best_err > tolerance:
            raise ValueError(
                f"CA atom index {ca_idx} failed coordinate match to FEFF ATOMS: "
                f"best error {best_err:.6f} A exceeds tolerance {tolerance:.6f} A"
            )

        atom_with_meta = dict(best_atom)
        atom_with_meta["source_xyz_index"] = ca_idx
        atom_with_meta["match_error"] = best_err
        matched.append(atom_with_meta)
        used_feff_indices.add(best_atom["index"])

    return matched


def has_atoms_block(path: Path) -> bool:
    if not path.is_file():
        return False
    try:
        with path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                if raw_line.strip().upper().startswith("ATOMS"):
                    return True
    except OSError:
        return False
    return False


def match_nlegs2_to_atoms(nlegs2_entries, atoms, tolerance=0.03):
    # Expand xmu rows by degeneracy in row order. Each expanded slot consumes one FEFF atom.
    slots = []
    for entry in nlegs2_entries:
        deg = max(1, int(entry.get("deg", 1)))
        for _ in range(deg):
            slots.append(
                {
                    "reff": float(entry["reff"]),
                    "sig2_tot": float(entry["sig2_tot"]),
                }
            )

    # Keep FEFF ATOMS order, excluding absorber at origin.
    ordered_atoms = [atom for atom in atoms if float(atom["distance"]) > 1e-8]
    used_atom_indices = set()

    matched = []
    for slot in slots:
        best_atom = None
        best_delta = None
        for atom in ordered_atoms:
            atom_idx = atom["index"]
            if atom_idx in used_atom_indices:
                continue
            delta = abs(float(atom["distance"]) - slot["reff"])
            if delta > tolerance:
                continue
            if best_delta is None or delta < best_delta:
                best_delta = delta
                best_atom = atom

        if best_atom is None:
            continue

        matched_atom = dict(best_atom)
        matched_atom["sig2_tot"] = slot["sig2_tot"]
        matched_atom["sig2_reff"] = slot["reff"]
        matched_atom["sig2_match_delta"] = float(best_delta)
        matched.append(matched_atom)
        used_atom_indices.add(best_atom["index"])

    if not matched:
        raise ValueError("Could not match nlegs=2 distances to any ATOMS entries")
    return matched


def find_sig2_for_distance(nlegs2_entries, distance, tolerance=0.03, require_within_tolerance=True):
    best_entry = None
    best_delta = None
    for entry in nlegs2_entries:
        delta = abs(float(distance) - float(entry["reff"]))
        if delta > tolerance:
            continue
        if best_delta is None or delta < best_delta:
            best_delta = delta
            best_entry = entry

    if best_entry is None and not require_within_tolerance:
        # Fall back to nearest available nlegs=2 path even when outside tolerance.
        if not nlegs2_entries:
            return None
        best_entry = min(nlegs2_entries, key=lambda item: abs(float(distance) - float(item["reff"])))
        best_delta = abs(float(distance) - float(best_entry["reff"]))

    if best_entry is None:
        return None
    return {
        "sig2_tot": float(best_entry["sig2_tot"]),
        "sig2_reff": float(best_entry["reff"]),
        "sig2_match_delta": float(best_delta),
    }


def pick_farthest_distinct_carbons(atoms, count=4, min_pair_separation=2.0):
    carbons = [atom for atom in atoms if atom["symbol"] == "C"]
    if len(carbons) < count:
        return sorted(carbons, key=lambda atom: atom["distance"], reverse=True)

    ordered = sorted(carbons, key=lambda atom: atom["distance"], reverse=True)
    selected = []
    for atom in ordered:
        if len(selected) == 0:
            selected.append(atom)
            if len(selected) == count:
                break
            continue

        coords = np.array([atom["x"], atom["y"], atom["z"]])
        too_close = False
        for prev in selected:
            prev_coords = np.array([prev["x"], prev["y"], prev["z"]])
            if np.linalg.norm(coords - prev_coords) < min_pair_separation:
                too_close = True
                break
        if not too_close:
            selected.append(atom)
            if len(selected) == count:
                break

    if len(selected) < count:
        selected = ordered[:count]
    return selected


def resolve_feff_inp_path(feff_dir: Path) -> Path:
    candidate = (feff_dir / "../../feff.inp").resolve()
    if candidate.exists():
        return candidate
    raise FileNotFoundError(f"Missing file: {candidate}")


def write_dw_dat(feff_dir: Path):
    xmu_path = feff_dir / "xmu.dat"
    if not xmu_path.exists():
        raise FileNotFoundError(f"Missing file: {xmu_path}")

    feff_path = resolve_feff_inp_path(feff_dir)
    run_root = feff_path.parent

    nlegs2_entries = parse_xmu_nlegs2_entries(xmu_path)
    atoms = parse_feff_atoms(feff_path)
    matched_atoms = match_nlegs2_to_atoms(nlegs2_entries, atoms)

    nearest_atoms = [atom for atom in matched_atoms if atom["distance"] > 1e-8]
    nearest_atoms = sorted(nearest_atoms, key=lambda atom: atom["distance"])[:4]

    latest_xyz = find_latest_xyz_file(run_root)
    comments_path = resolve_comments_file(run_root, latest_xyz)
    ca_indices = parse_ca_indices_from_comments(comments_path)
    xyz_atoms = parse_xyz_atoms(latest_xyz)
    ca_atoms_by_coords = match_ca_indices_to_feff_atoms(ca_indices, xyz_atoms, atoms)

    matched_by_index = {atom["index"]: atom for atom in matched_atoms}
    ca_atoms = []
    missing_sig2 = []
    ca_fallback_count = 0
    for ca_atom in ca_atoms_by_coords:
        idx = ca_atom["index"]
        if idx in matched_by_index:
            with_sig2 = dict(matched_by_index[idx])
        else:
            # Fallback: assign by closest nlegs=2 reff when atom-index pairing is ambiguous.
            fallback = find_sig2_for_distance(
                nlegs2_entries,
                ca_atom["distance"],
                tolerance=0.05,
                require_within_tolerance=True,
            )
            if fallback is None:
                fallback = find_sig2_for_distance(
                    nlegs2_entries,
                    ca_atom["distance"],
                    tolerance=0.05,
                    require_within_tolerance=False,
                )
            if fallback is None:
                missing_sig2.append(idx)
                continue
            with_sig2 = dict(ca_atom)
            with_sig2.update(fallback)
            with_sig2["sig2_extrapolated"] = fallback["sig2_match_delta"] > 0.05
            ca_fallback_count += 1

        with_sig2["source_xyz_index"] = ca_atom["source_xyz_index"]
        with_sig2["match_error"] = ca_atom["match_error"]
        ca_atoms.append(with_sig2)

    if missing_sig2:
        missing_str = ", ".join(str(val) for val in sorted(missing_sig2))
        raise ValueError(
            "Matched CA atoms in FEFF ATOMS are missing nlegs=2 sig2_tot assignments for "
            f"atom indices: {missing_str}"
        )

    out_path = feff_dir / "dw.dat"
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("# Derived from xmu.dat nlegs=2 reff distances matched to feff.inp ATOMS\n")
        handle.write("# Closest four atoms to Zn/origin\n")
        handle.write("group symbol x y z distance sig2_tot atom_index\n")
        for atom in nearest_atoms:
            handle.write(
                "nearest "
                f"{atom['symbol']} {atom['x']:.5f} {atom['y']:.5f} {atom['z']:.5f} "
                f"{atom['distance']:.5f} {atom['sig2_tot']:.5f} {atom['index']}\n"
            )

        handle.write(
            f"# CA atoms from {comments_path.name}, mapped from Zn-centered {latest_xyz.name} to FEFF ATOMS\n"
        )
        if ca_fallback_count:
            handle.write(
                "# note: some CA sig2_tot values used nearest available nlegs=2 reff "
                "because no close nlegs=2 path existed\n"
            )
        for atom in ca_atoms:
            handle.write(
                "ca "
                f"{atom['symbol']} {atom['x']:.5f} {atom['y']:.5f} {atom['z']:.5f} "
                f"{atom['distance']:.5f} {atom['sig2_tot']:.5f} {atom['index']}\n"
            )

    if ca_fallback_count:
        print(
            f"warning: assigned nearest nlegs=2 sig2_tot for {ca_fallback_count} CA atom(s) "
            f"in {feff_dir / 'dw.dat'} because close nlegs=2 matches were unavailable"
        )

    return out_path


def xftf_larch(k, chi, kmin, kmax, dk, kweight, kstep, rmax_out, window):
    from larch import Group
    from larch.xafs import xftf

    grp = Group()
    grp.k = k
    grp.chi = chi
    xftf(
        grp.k,
        grp.chi,
        kmin=kmin,
        kmax=kmax,
        dk=dk,
        kweight=kweight,
        kstep=kstep,
        rmax_out=rmax_out,
        window=window,
        group=grp,
    )
    return grp.r, grp.chir


def apply_plot_style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.tick_params(direction="in", width=2, length=8)


def run_for_feff_dir(feff_dir: Path, args: argparse.Namespace):
    xanes_path = feff_dir / "xanes_K.dat"
    exafs_path = feff_dir / "exafs_K.dat"

    if not xanes_path.exists():
        raise FileNotFoundError(f"Missing file: {xanes_path}")
    if not exafs_path.exists():
        raise FileNotFoundError(f"Missing file: {exafs_path}")

    x_omega, _, _, x_mu, _, _ = load_feff_table(xanes_path)
    _, _, ex_k, _, _, ex_chi = load_feff_table(exafs_path)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_omega, x_mu, lw=2)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"$\mu$")
    ax.set_title("XANES")
    apply_plot_style(ax)
    fig.tight_layout()
    xanes_png = feff_dir / "xanes_K.png"
    fig.savefig(xanes_png, dpi=300)
    if not args.show:
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(ex_k, ex_chi, lw=2)
    ax.set_xlabel(r"$k\ (1/\AA)$")
    ax.set_ylabel(r"$\chi(k)$")
    ax.set_title("EXAFS")
    apply_plot_style(ax)
    fig.tight_layout()
    exafs_png = feff_dir / "exafs_K.png"
    fig.savefig(exafs_png, dpi=300)
    if not args.show:
        plt.close(fig)

    dw_dat = write_dw_dat(feff_dir)

    saved_outputs = [xanes_png, exafs_png, dw_dat]
    if not args.skip_fft:
        chi_path = feff_dir / "chi.dat"
        if not chi_path.exists():
            raise FileNotFoundError(f"Missing file: {chi_path}")

        k, chi = load_chi_dat(chi_path)

        r, chir = xftf_larch(
            k,
            chi,
            kmin=args.kmin,
            kmax=args.kmax,
            dk=args.dk,
            kweight=args.kweight,
            kstep=args.kstep,
            rmax_out=args.rmax,
            window=args.window,
        )

        chir_mag = np.abs(chir)
        chir_re = chir.real
        chir_im = chir.imag

        out_dat = feff_dir / "chi_R.dat"
        header = "r  chir_mag  chir_re  chir_im"
        np.savetxt(out_dat, np.column_stack([r, chir_mag, chir_re, chir_im]), header=header)
        saved_outputs.append(out_dat)

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(r, chir_mag, lw=2)
        ax.set_xlabel(r"$R\ (\AA)$")
        ax.set_ylabel(r"$|\chi(R)|$")
        ax.set_title("EXAFS FT")
        apply_plot_style(ax)
        fig.tight_layout()

        out_png = feff_dir / "chi_R.png"
        fig.savefig(out_png, dpi=300)
        saved_outputs.append(out_png)
        if not args.show:
            plt.close(fig)

    if args.show:
        plt.show()
    else:
        for out_path in saved_outputs:
            print(f"Saved: {out_path}")


def is_feff_dir(path: Path) -> bool:
    return path.is_dir() and (path / "xanes_K.dat").is_file() and (path / "exafs_K.dat").is_file()


def find_feff_dir_in_tree(base: Path) -> Path | None:
    direct_candidates = [
        base,
        base / "Corvus3_cfavg_xas" / "Corvus1Zn_FEFF",
    ]
    for candidate in direct_candidates:
        if is_feff_dir(candidate):
            return candidate

    for child in base.iterdir() if base.is_dir() else []:
        if not child.is_dir():
            continue
        if child.name.startswith("working") and is_feff_dir(
            child / "Corvus3_cfavg_xas" / "Corvus1Zn_FEFF"
        ):
            return child / "Corvus3_cfavg_xas" / "Corvus1Zn_FEFF"
    return None


def has_working_output_pair(system_dir: Path) -> bool:
    name = system_dir.name
    return (system_dir / f"working-{name}").is_dir() and (system_dir / f"output-{name}").is_dir()


def is_process_target(system_dir: Path) -> bool:
    if not system_dir.is_dir():
        return False
    if has_working_output_pair(system_dir):
        return True
    return find_feff_dir_in_tree(system_dir) is not None


def resolve_system_targets(parent_dir: Path, recursive: bool) -> List[Path]:
    if not recursive:
        if not is_process_target(parent_dir):
            raise FileNotFoundError(
                f"No processing target found in {parent_dir}. "
                "Expected working/output directories or FEFF working files."
            )
        return [parent_dir]

    if is_process_target(parent_dir):
        return [parent_dir]

    targets = []
    for child in sorted(parent_dir.iterdir()):
        if child.is_dir() and is_process_target(child):
            targets.append(child)
    return targets


def move_unprocessed_contents_to_working(system_dir: Path, working_dir: Path, output_dir: Path):
    skip_names = {working_dir.name, output_dir.name}
    for entry in list(system_dir.iterdir()):
        if entry.name in skip_names:
            continue
        shutil.move(str(entry), str(working_dir / entry.name))


def copy_if_exists(src: Path, dst: Path, label: str):
    if src.is_file():
        shutil.copy2(src, dst)
        print(f"Copied {label}: {src} -> {dst}")
    else:
        print(f"warning: missing {label}: {src}")


def process_system_dir(system_dir: Path, args: argparse.Namespace):
    name = system_dir.name
    output_dir = system_dir / f"output-{name}"
    working_dir = system_dir / f"working-{name}"

    already_processed = output_dir.is_dir() and working_dir.is_dir()
    output_dir.mkdir(exist_ok=True)
    working_dir.mkdir(exist_ok=True)

    if not already_processed:
        move_unprocessed_contents_to_working(system_dir, working_dir, output_dir)

    feff_dir = find_feff_dir_in_tree(system_dir)
    if feff_dir is None:
        raise FileNotFoundError(
            f"Could not locate FEFF directory under {system_dir}. "
            "Expected xanes_K.dat and exafs_K.dat."
        )

    run_for_feff_dir(feff_dir, args)

    xyz_src_candidates = [working_dir / f"{name}.xyz", system_dir / f"{name}.xyz"]
    xyz_src = next((path for path in xyz_src_candidates if path.is_file()), xyz_src_candidates[0])

    chi_src = feff_dir / "chi_R.dat"
    dw_src = feff_dir / "dw.dat"
    exafs_src = feff_dir / "exafs_K.dat"

    copy_if_exists(xyz_src, output_dir / f"{name}.xyz", "xyz")
    copy_if_exists(chi_src, output_dir / f"chi-R-{name}.dat", "chi_R.dat")
    copy_if_exists(dw_src, output_dir / f"dw-{name}.dat", "dw.dat")
    copy_if_exists(exafs_src, output_dir / f"exafs-k-{name}.dat", "exafs_k.dat")


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Process FEFF output(s): generate XANES/EXAFS plots, chi_R.dat, dw.dat, "
            "and copy selected outputs to output-<name> directories."
        )
    )
    parser.add_argument(
        "parent_dir",
        type=Path,
        help=(
            "Directory to process. In recursive mode, this is scanned for child system directories "
            "when it is not itself a processing target."
        ),
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Process all matching child directories when parent_dir is not itself a target.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots after saving.",
    )
    parser.add_argument(
        "--skip-fft",
        action="store_true",
        help="Skip chi(k) to chi(R) Fourier transform.",
    )
    parser.add_argument("--kmin", type=float, default=3.0)
    parser.add_argument("--kmax", type=float, default=11.0)
    parser.add_argument("--dk", type=float, default=3.0)  # dk=1
    parser.add_argument("--kweight", type=int, default=2)
    parser.add_argument("--kstep", type=float, default=0.05)
    parser.add_argument("--rmax", type=float, default=6.0)
    parser.add_argument("--window", type=str, default="kaiser")  #window="hanning"
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    parent_dir = args.parent_dir
    if not parent_dir.is_dir():
        print(f"error: directory not found: {parent_dir}")
        return 2

    targets = resolve_system_targets(parent_dir, recursive=args.recursive)
    if not targets:
        print(
            f"error: no processable directories found under {parent_dir}. "
            "Expected working/output dirs or FEFF working files."
        )
        return 2

    failures = 0
    for target in targets:
        print(f"\nProcessing: {target}")
        try:
            process_system_dir(target, args)
        except Exception as exc:
            failures += 1
            print(f"error: failed processing {target}: {exc}")

    if failures:
        print(f"Completed with {failures} failure(s)")
        return 1
    print("Completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())