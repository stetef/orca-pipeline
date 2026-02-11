#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 16,
})


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


def xftf_larch(k, chi, kmin, kmax, dk, kweight, kstep, nfft, rmax_out):
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
        nfft=nfft,
        rmax_out=rmax_out,
        window="kaiser",
        group=grp,
    )
    return grp.r, grp.chir


def apply_plot_style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.tick_params(direction="in", width=2, length=8)


def main():
    parser = argparse.ArgumentParser(
        description="Plot xanes_K.dat and exafs_K.dat from a FEFF output directory."
    )
    parser.add_argument(
        "directory",
        type=Path,
        help="Directory containing xanes_K.dat and exafs_K.dat",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots after saving",
    )
    parser.add_argument(
        "--skip-fft",
        action="store_true",
        help="Skip chi(k) to chi(R) Fourier transform",
    )
    parser.add_argument("--kmin", type=float, default=2.0)
    parser.add_argument("--kmax", type=float, default=10.0)
    parser.add_argument("--dk", type=float, default=2.0)
    parser.add_argument("--kweight", type=int, default=3)
    parser.add_argument("--kstep", type=float, default=0.05)
    parser.add_argument("--nfft", type=int, default=2048)
    parser.add_argument("--rmax", type=float, default=10.0)
    args = parser.parse_args()

    directory = args.directory
    xanes_path = directory / "xanes_K.dat"
    exafs_path = directory / "exafs_K.dat"

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
    xanes_png = directory / "xanes_K.png"
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
    exafs_png = directory / "exafs_K.png"
    fig.savefig(exafs_png, dpi=300)
    if not args.show:
        plt.close(fig)

    saved_outputs = [xanes_png, exafs_png]
    if not args.skip_fft:
        chi_path = directory / "chi.dat"
        if not chi_path.exists():
            raise FileNotFoundError(f"Missing file: {chi_path}")

        k, chi = load_chi_dat(chi_path)
        k_shift = k - k.min()

        r, chir = xftf_larch(
            k_shift,
            chi,
            kmin=args.kmin,
            kmax=args.kmax,
            dk=args.dk,
            kweight=args.kweight,
            kstep=args.kstep,
            nfft=args.nfft,
            rmax_out=args.rmax,
        )

        chir_mag = np.abs(chir)
        chir_re = chir.real
        chir_im = chir.imag

        out_dat = directory / "chi_R.dat"
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

        out_png = directory / "chi_R.png"
        fig.savefig(out_png, dpi=300)
        saved_outputs.append(out_png)
        if not args.show:
            plt.close(fig)

    if args.show:
        plt.show()
    else:
        for out_path in saved_outputs:
            print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
