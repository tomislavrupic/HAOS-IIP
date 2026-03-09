#!/usr/bin/env python3

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PLOTS = ROOT / "plots"
PLOTS.mkdir(exist_ok=True)
MPLCONFIG = ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    n_side = 5
    grid = np.linspace(0.0, 1.0, n_side)
    X, Y, Z = np.meshgrid(grid, grid, grid, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    h = 1.0 / (n_side - 1)
    epsilon_k = 1.5 * (h * h)

    diff = points[:, None, :] - points[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    cutoff = 2.5 * math.sqrt(epsilon_k)
    mask = (d2 > 0.0) & (d2 <= cutoff * cutoff)

    A = np.zeros_like(d2)
    A[mask] = np.exp(-d2[mask] / (2.0 * epsilon_k))
    D = np.diag(np.sum(A, axis=1))
    L = D - A

    evals, _ = np.linalg.eigh(L)
    evals = np.sort(evals)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(12), evals[:12], marker="o")
    ax.set_xlabel("eigen-index")
    ax.set_ylabel("eigenvalue")
    ax.set_title("starter scalar Laplacian spectrum")
    ax.grid(alpha=0.25)
    out = PLOTS / "starter_laplacian_spectrum.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print("epsilon_k =", epsilon_k)
    print("first eigenvalues =", [round(float(v), 6) for v in evals[:8]])
    print("plot =", out)


if __name__ == "__main__":
    main()
