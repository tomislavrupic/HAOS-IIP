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
    epsilon_k = h * h

    diff = points[:, None, :] - points[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    mask = (d2 > 0.0) & (d2 <= (1.05 * h) ** 2)

    A = np.zeros_like(d2)
    A[mask] = np.exp(-d2[mask] / (2.0 * epsilon_k))

    flux_per_plaquette = 0.2
    B = 2.0 * math.pi * flux_per_plaquette / (h * h)
    x_mid = 0.5 * (points[:, None, 0] + points[None, :, 0])
    theta = B * x_mid * diff[:, :, 1]
    theta = 0.5 * (theta - theta.T)
    U = np.exp(1j * theta)

    D = np.diag(np.sum(A, axis=1))
    L_theta = D - A * U
    evals, evecs = np.linalg.eigh(L_theta)
    order = np.argsort(evals.real)
    evals = evals[order]
    evecs = evecs[:, order]

    mode = evecs[:, 0]
    phase = np.angle(mode)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(points[:, 0], points[:, 1], points[:, 2], c=phase, cmap="twilight", s=45)
    ax.set_title("starter phase-dressed lowest mode")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    fig.colorbar(sc, ax=ax, shrink=0.7, label="phase")
    out = PLOTS / "starter_gauge_mode.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print("epsilon_k =", epsilon_k)
    print("first eigenvalues =", [round(float(v.real), 6) for v in evals[:8]])
    print("plot =", out)


if __name__ == "__main__":
    main()
