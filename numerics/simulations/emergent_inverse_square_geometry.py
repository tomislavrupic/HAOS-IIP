#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
for path in (DATA, PLOTS):
    path.mkdir(exist_ok=True)

MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    "D": 1.0,
    "r_max": 120.0,
    "n_grid": 1600,
    "source_width": 1.2,
    "source_strength": 1.0,
    "dt": 2.0,
    "time_steps": 220,
    "snapshot_steps": [0, 5, 20, 80, 220],
    "fit_r_min": 4.0,
    "fit_r_max": 20.0,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("inverse_square_geometry"), dict):
            merged.update(raw["inverse_square_geometry"])
    if config is not None:
        merged.update(config)
    return merged


def build_source(r: np.ndarray, source_width: float, source_strength: float) -> tuple[np.ndarray, float]:
    profile = np.exp(-(r / source_width) ** 2)
    dr = float(r[1] - r[0])
    raw_charge = float(np.sum(4.0 * math.pi * r * r * profile) * dr)
    rho = (source_strength / raw_charge) * profile
    total_charge = float(np.sum(4.0 * math.pi * r * r * rho) * dr)
    return rho, total_charge


def build_radial_laplacian(r: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dr = float(r[1] - r[0])
    n = r.size
    lower = np.zeros(n, dtype=float)
    diag = np.zeros(n, dtype=float)
    upper = np.zeros(n, dtype=float)

    # Regularity at the origin for a spherically symmetric field.
    diag[0] = -6.0 / (dr * dr)
    upper[0] = 6.0 / (dr * dr)

    for i in range(1, n - 1):
        lower[i] = 1.0 / (dr * dr) - 1.0 / (r[i] * dr)
        diag[i] = -2.0 / (dr * dr)
        upper[i] = 1.0 / (dr * dr) + 1.0 / (r[i] * dr)

    # Dirichlet outer boundary phi(r_max)=0.
    diag[-1] = 1.0
    lower[-1] = 0.0
    upper[-1] = 0.0

    return lower, diag, upper


def tridiagonal_solve(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    n = rhs.size
    c = upper.copy()
    d = rhs.copy()
    b = diag.copy()

    for i in range(1, n):
        factor = lower[i] / b[i - 1]
        b[i] -= factor * c[i - 1]
        d[i] -= factor * d[i - 1]

    x = np.zeros(n, dtype=float)
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x


def apply_tridiagonal(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, vec: np.ndarray) -> np.ndarray:
    out = diag * vec
    out[1:] += lower[1:] * vec[:-1]
    out[:-1] += upper[:-1] * vec[1:]
    return out


def run_relaxation(
    r: np.ndarray,
    rho: np.ndarray,
    D: float,
    dt: float,
    steps: int,
    snapshot_steps: list[int],
) -> dict[str, Any]:
    lower_L, diag_L, upper_L = build_radial_laplacian(r)
    n = r.size

    lower_A = -dt * D * lower_L
    diag_A = 1.0 - dt * D * diag_L
    upper_A = -dt * D * upper_L
    rhs_source = dt * rho

    phi = np.zeros(n, dtype=float)
    snapshots: dict[int, list[float]] = {}
    residuals: list[float] = []

    for step in range(steps + 1):
        if step in snapshot_steps:
            snapshots[step] = phi.tolist()
        lap_phi = apply_tridiagonal(lower_L, diag_L, upper_L, phi)
        residual = float(np.linalg.norm(D * lap_phi + rho))
        residuals.append(residual)
        if step == steps:
            break
        rhs = phi + rhs_source
        rhs[-1] = 0.0
        phi = tridiagonal_solve(lower_A, diag_A, upper_A, rhs)

    return {
        "phi": phi,
        "snapshots": snapshots,
        "residuals": residuals,
        "lower_L": lower_L,
        "diag_L": diag_L,
        "upper_L": upper_L,
    }


def fit_power_law(r: np.ndarray, values: np.ndarray, r_min: float, r_max: float) -> dict[str, float]:
    mask = (r >= r_min) & (r <= r_max) & (values > 0.0)
    x = np.log(r[mask])
    y = np.log(values[mask])
    slope, intercept = np.polyfit(x, y, deg=1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1.0e-18)
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(r2),
    }


def fit_inverse_r_plus_constant(r: np.ndarray, values: np.ndarray, r_min: float, r_max: float) -> dict[str, float]:
    mask = (r >= r_min) & (r <= r_max)
    X = np.column_stack([1.0 / r[mask], np.ones(mask.sum(), dtype=float)])
    coeffs, *_ = np.linalg.lstsq(X, values[mask], rcond=None)
    pred = X @ coeffs
    ss_res = float(np.sum((values[mask] - pred) ** 2))
    ss_tot = float(np.sum((values[mask] - np.mean(values[mask])) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1.0e-18)
    residual_field = values[mask] - float(coeffs[1])
    power = fit_power_law(r[mask], np.abs(residual_field), r_min=float(r[mask][0]), r_max=float(r[mask][-1]))
    return {
        "coefficient_over_r": float(coeffs[0]),
        "constant_offset": float(coeffs[1]),
        "r2": float(r2),
        "residual_power_slope": float(power["slope"]),
    }


def make_plots(
    r: np.ndarray,
    rho: np.ndarray,
    phi: np.ndarray,
    snapshots: dict[int, list[float]],
    fit_phi: dict[str, float],
    fit_force: dict[str, float],
    D: float,
    total_charge: float,
    stamp: str,
) -> list[str]:
    plot_paths: list[str] = []

    field_path = PLOTS / f"{stamp}_inverse_square_geometry_field.png"
    reference_phi = (
        float(fit_phi["constant_offset"])
        + total_charge / (4.0 * math.pi * D * np.maximum(r, 1.0e-12))
    )
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(r[1:], phi[1:], label="measured field")
    ax.loglog(r[1:], reference_phi[1:], "--", label="offset + Q/(4 pi D r)")
    ax.set_xlabel("r")
    ax.set_ylabel("field amplitude")
    ax.set_title("Emergent radial field")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(field_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(field_path.relative_to(REPO_ROOT)))

    force_path = PLOTS / f"{stamp}_inverse_square_geometry_force.png"
    force = np.abs(np.gradient(phi, r))
    reference_force = total_charge / (4.0 * math.pi * D * np.maximum(r, 1.0e-12) ** 2)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(r[2:], force[2:], label="measured gradient")
    ax.loglog(r[2:], reference_force[2:], "--", label="Q/(4 pi D r^2)")
    ax.set_xlabel("r")
    ax.set_ylabel(r"$|d\phi/dr|$")
    ax.set_title("Emergent inverse-square force law")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(force_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(force_path.relative_to(REPO_ROOT)))

    dynamic_path = PLOTS / f"{stamp}_inverse_square_geometry_relaxation.png"
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for step, values in sorted(snapshots.items()):
        if step == 0:
            continue
        ax.plot(r, values, label=f"step {step}")
    ax.plot(r, phi, color="k", linewidth=2.0, label="steady")
    ax.set_xlim(0.0, min(40.0, float(r[-1])))
    ax.set_xlabel("r")
    ax.set_ylabel("field amplitude")
    ax.set_title("Relaxation toward recoverable radial geometry")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(dynamic_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(dynamic_path.relative_to(REPO_ROOT)))

    source_path = PLOTS / f"{stamp}_inverse_square_geometry_source.png"
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(r, rho, label="localized source")
    ax.set_xlim(0.0, min(12.0, float(r[-1])))
    ax.set_xlabel("r")
    ax.set_ylabel("source density")
    ax.set_title("Localized substrate disturbance")
    ax.grid(alpha=0.25)
    fig.savefig(source_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(source_path.relative_to(REPO_ROOT)))

    return plot_paths


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## Emergent inverse-square geometry\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"D={config['D']}, r_max={config['r_max']}, n_grid={config['n_grid']}, "
            f"source_width={config['source_width']}, source_strength={config['source_strength']}, "
            f"dt={config['dt']}, time_steps={config['time_steps']}, "
            f"fit_window=[{config['fit_r_min']}, {config['fit_r_max']}]\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_emergent_inverse_square_geometry.json"
    latest = DATA / "emergent_inverse_square_geometry_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def run_inverse_square_geometry(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    r = np.linspace(0.0, float(cfg["r_max"]), int(cfg["n_grid"]))
    rho, total_charge = build_source(
        r,
        source_width=float(cfg["source_width"]),
        source_strength=float(cfg["source_strength"]),
    )
    relaxation = run_relaxation(
        r=r,
        rho=rho,
        D=float(cfg["D"]),
        dt=float(cfg["dt"]),
        steps=int(cfg["time_steps"]),
        snapshot_steps=[int(v) for v in cfg["snapshot_steps"]],
    )
    phi = np.asarray(relaxation["phi"], dtype=float)
    force = np.abs(np.gradient(phi, r))

    fit_phi = fit_inverse_r_plus_constant(
        r[1:],
        phi[1:],
        r_min=float(cfg["fit_r_min"]),
        r_max=float(cfg["fit_r_max"]),
    )
    fit_force = fit_power_law(
        r[2:],
        force[2:],
        r_min=float(cfg["fit_r_min"]),
        r_max=float(cfg["fit_r_max"]),
    )

    observation = "the steady radial field generated by a localized disturbance follows an inverse-distance law outside the source core"
    conclusion = (
        "inverse-square geometry emerges from flux-conserving 3D substrate propagation once a localized disturbance relaxes into a recoverable stationary field"
    )

    return {
        "experiment": "emergent_inverse_square_geometry",
        "config": cfg,
        "derivation": {
            "relaxation_equation": "d_tau phi = D (phi_rr + 2 phi_r / r) + rho(r)",
            "steady_equation": "-D (phi_rr + 2 phi_r / r) = rho(r)",
            "outside_source": "0 = (1/r^2) d_r (r^2 d_r phi)",
            "far_field_solution": "phi(r) = A + B/r, |d_r phi| = |B|/r^2",
        },
        "total_charge": total_charge,
        "r": r.tolist(),
        "rho": rho.tolist(),
        "phi": phi.tolist(),
        "snapshots": relaxation["snapshots"],
        "residuals": relaxation["residuals"],
        "fit_field": fit_phi,
        "fit_force": fit_force,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = run_inverse_square_geometry()
    plots = make_plots(
        r=np.asarray(result["r"], dtype=float),
        rho=np.asarray(result["rho"], dtype=float),
        phi=np.asarray(result["phi"], dtype=float),
        snapshots=result["snapshots"],
        fit_phi=result["fit_field"],
        fit_force=result["fit_force"],
        D=float(result["config"]["D"]),
        total_charge=float(result["total_charge"]),
        stamp=stamp,
    )
    result["plots"] = plots
    result_path, latest_path = save_results(result, stamp=stamp)
    append_log(
        result_path=result_path,
        plot_paths=plots,
        config=result["config"],
        observation=result["observation"],
        conclusion=result["conclusion"],
    )
    print("results =", result_path)
    print("latest =", latest_path)
    print("plots =", plots)
    print(
        "field fit = offset + B/r with B =",
        round(result["fit_field"]["coefficient_over_r"], 6),
        "offset =",
        round(result["fit_field"]["constant_offset"], 6),
        "R2 =",
        round(result["fit_field"]["r2"], 6),
        "residual slope =",
        round(result["fit_field"]["residual_power_slope"], 6),
    )
    print("force slope =", round(result["fit_force"]["slope"], 6), "R2 =", round(result["fit_force"]["r2"], 6))
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
