#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

METHODS = [
    ("res_jacobi.txt",   "Jacobi",             "tab:blue",   "-", 1.8),
    ("res_seidel.txt",   "Seidel",             "tab:orange", "-", 1.8),
    ("res_sor.txt",      "SOR (w = 1.5)",      "tab:green",  "-", 1.8),
    ("res_grad.txt",     "Gradient descent",   "tab:red",    "-", 1.8),
    ("res_minres.txt",   "Minimal residual",   "tab:purple", "-", 1.8),
    ("res_cg.txt",       "Conjugate gradient", "tab:brown",  "-", 1.8),
    ("res_bicgstab.txt", "BiCGStab",           "tab:pink",   "-", 1.8),
]

DATA_DIR = Path(".")


def load_history(path: Path):
    if not path.exists():
        print(f"file not found: {path}")
        return None
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1]


def plot_all(outfile="convergence_log.png", show=True):
    fig, ax = plt.subplots(figsize=(11, 7))
    loaded_any = False

    for fname, label, color, ls, lw in METHODS:
        loaded = load_history(DATA_DIR / fname)
        if loaded is None:
            continue
        it, res = loaded
        loaded_any = True
        ax.plot(it, res, ls, color=color, lw=lw, label=label)

    if not loaded_any:
        raise SystemExit("no res_*.txt files found")

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax.set_xlabel("iteration")
    ax.set_ylabel(r"residual $\|b - A x\|_2$  (log scale)")
    ax.set_title("Residual convergence")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend(loc="best", framealpha=0.9)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"saved: {outfile}")
    if show:
        plt.show()


def main():
    plot_all()


if __name__ == "__main__":
    main()
