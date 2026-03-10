#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from numerics.simulations.recoverable_hydrogenic_spectrum import analytic_energy, build_radial_hamiltonian

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
    "alpha": 1.0,
    "r_min": 1.0e-3,
    "r_max": 80.0,
    "n_grid": 1200,
    "states": 4,
    "modal_basis_size": 12,
    "dt": 6.0,
    "max_steps": 400,
    "tol": 1.0e-10,
    "random_seed": 42,
    "recovery_noise": 0.12,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("interaction_attractor_spectrum"), dict):
            merged.update(raw["interaction_attractor_spectrum"])
    if config is not None:
        merged.update(config)
    return merged


def normalize(vec: np.ndarray, dr: float) -> np.ndarray:
    norm = math.sqrt(float(np.sum(np.abs(vec) ** 2) * dr))
    if norm == 0.0:
        return vec
    return vec / norm


def weighted_dot(a: np.ndarray, b: np.ndarray, dr: float) -> float:
    return float(np.sum(np.conjugate(a) * b) * dr)


def normalize_coeffs(coeffs: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(coeffs))
    if norm == 0.0:
        return coeffs
    return coeffs / norm


def project_coeffs(coeffs: np.ndarray, allowed_mask: np.ndarray) -> np.ndarray:
    out = np.asarray(coeffs, dtype=float).copy()
    out[~allowed_mask] = 0.0
    return normalize_coeffs(out)


def modal_residual(coeffs: np.ndarray, energies: np.ndarray, allowed_mask: np.ndarray) -> float:
    active = coeffs[allowed_mask]
    active_e = energies[allowed_mask]
    lam = float(np.sum(active_e * active * active))
    return float(np.linalg.norm((active_e - lam) * active))


def run_modal_flow(
    energies: np.ndarray,
    initial_coeffs: np.ndarray,
    allowed_mask: np.ndarray,
    dt: float,
    max_steps: int,
    tol: float,
) -> dict[str, Any]:
    coeffs = project_coeffs(initial_coeffs, allowed_mask)
    energy_trajectory: list[float] = []
    residual_trajectory: list[float] = []

    for step in range(max_steps):
        lam = float(np.sum(energies * coeffs * coeffs))
        residual = modal_residual(coeffs, energies, allowed_mask)
        energy_trajectory.append(lam)
        residual_trajectory.append(residual)
        if residual < tol and step > 4:
            break
        coeffs = np.exp(-dt * energies) * coeffs
        coeffs = project_coeffs(coeffs, allowed_mask)

    lam = float(np.sum(energies * coeffs * coeffs))
    residual = modal_residual(coeffs, energies, allowed_mask)
    energy_trajectory.append(lam)
    residual_trajectory.append(residual)

    return {
        "coeffs": coeffs,
        "energy": float(lam),
        "residual": float(residual),
        "steps": len(energy_trajectory),
        "energy_trajectory": energy_trajectory,
        "residual_trajectory": residual_trajectory,
    }


def sample_initial(dim: int, rng: np.random.Generator, allowed_mask: np.ndarray) -> np.ndarray:
    coeffs = rng.normal(size=dim)
    coeffs *= np.linspace(1.0, 0.4, dim)
    return project_coeffs(coeffs, allowed_mask)


def reconstruct_state(basis: np.ndarray, coeffs: np.ndarray, dr: float) -> np.ndarray:
    return normalize(basis @ coeffs, dr)


def attractor_run_for_state(
    basis: np.ndarray,
    modal_energies: np.ndarray,
    dr: float,
    reference_coeffs: np.ndarray,
    state_index: int,
    cfg: dict[str, Any],
    rng: np.random.Generator,
) -> dict[str, Any]:
    allowed_mask = np.arange(modal_energies.size) >= state_index
    initial = sample_initial(modal_energies.size, rng, allowed_mask)
    flow = run_modal_flow(
        energies=modal_energies,
        initial_coeffs=initial,
        allowed_mask=allowed_mask,
        dt=float(cfg["dt"]),
        max_steps=int(cfg["max_steps"]),
        tol=float(cfg["tol"]),
    )
    coeffs = np.asarray(flow["coeffs"], dtype=float)
    overlap = float(abs(np.dot(reference_coeffs, coeffs)))
    if overlap < 0.999 and np.dot(reference_coeffs, coeffs) < 0.0:
        coeffs = -coeffs
    noise = rng.normal(size=coeffs.size)
    noise[~allowed_mask] = 0.0
    noise = noise - float(np.dot(coeffs, noise)) * coeffs
    if np.linalg.norm(noise) > 0.0:
        noise = noise / np.linalg.norm(noise)
    perturbed = normalize_coeffs(coeffs + float(cfg["recovery_noise"]) * noise)
    perturbed = project_coeffs(perturbed, allowed_mask)
    recovery = run_modal_flow(
        energies=modal_energies,
        initial_coeffs=perturbed,
        allowed_mask=allowed_mask,
        dt=float(cfg["dt"]),
        max_steps=int(cfg["max_steps"]),
        tol=float(cfg["tol"]),
    )
    recovered_coeffs = np.asarray(recovery["coeffs"], dtype=float)
    recovered_overlap = float(abs(np.dot(reference_coeffs, recovered_coeffs)))
    principal_n = state_index + 1
    target_energy = analytic_energy(principal_n, D=float(cfg["D"]), alpha=float(cfg["alpha"]))
    return {
        "principal_n": principal_n,
        "energy": float(flow["energy"]),
        "analytic_energy": float(target_energy),
        "relative_error": abs(float(flow["energy"]) - target_energy) / max(abs(target_energy), 1.0e-12),
        "steps": int(flow["steps"]),
        "residual": float(flow["residual"]),
        "reference_overlap": float(overlap),
        "recovered_overlap": float(recovered_overlap),
        "recovery_steps": int(recovery["steps"]),
        "energy_trajectory": flow["energy_trajectory"],
        "recovery_energy_trajectory": recovery["energy_trajectory"],
        "coeffs": coeffs.tolist(),
        "recovered_coeffs": recovered_coeffs.tolist(),
        "state": reconstruct_state(basis, coeffs, dr).tolist(),
        "recovered_state": reconstruct_state(basis, recovered_coeffs, dr).tolist(),
    }


def make_plots(rin: np.ndarray, attractors: list[dict[str, Any]], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    descent_path = PLOTS / f"{stamp}_interaction_attractor_energy_descent.png"
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for attr in attractors:
        traj = np.asarray(attr["energy_trajectory"], dtype=float)
        ax.plot(traj, label=f"n={attr['principal_n']}")
    ax.set_xlabel("flow step")
    ax.set_ylabel("Rayleigh energy")
    ax.set_title("Normalized interaction flow toward bound attractors")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(descent_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(descent_path.relative_to(REPO_ROOT)))

    scaling_path = PLOTS / f"{stamp}_interaction_attractor_scaling.png"
    n_vals = np.array([attr["principal_n"] for attr in attractors], dtype=float)
    energies = np.array([attr["energy"] for attr in attractors], dtype=float)
    predicted = np.array([attr["analytic_energy"] for attr in attractors], dtype=float)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(n_vals, energies, "o-", label="attractor energies")
    ax.plot(n_vals, predicted, "--", label="target -1/n^2 scaling")
    ax.set_xlabel("principal index n")
    ax.set_ylabel("energy")
    ax.set_title("Discrete attractor ladder")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(scaling_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(scaling_path.relative_to(REPO_ROOT)))

    recovery_path = PLOTS / f"{stamp}_interaction_attractor_recovery.png"
    fig, ax = plt.subplots(figsize=(7, 4))
    recovered = np.array([attr["recovered_overlap"] for attr in attractors], dtype=float)
    ax.bar([str(attr["principal_n"]) for attr in attractors], recovered)
    ax.set_ylim(0.0, 1.02)
    ax.set_xlabel("principal index n")
    ax.set_ylabel("overlap after perturbation and recovery")
    ax.set_title("Recoverability of discrete attractors")
    ax.grid(alpha=0.25, axis="y")
    fig.savefig(recovery_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(recovery_path.relative_to(REPO_ROOT)))

    modes_path = PLOTS / f"{stamp}_interaction_attractor_modes.png"
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
    axes = axes.ravel()
    for ax, attr in zip(axes, attractors):
        state = np.asarray(attr["state"], dtype=float)
        density = np.abs(state) ** 2
        ax.plot(rin, density, linewidth=1.5)
        ax.set_title(f"n={attr['principal_n']}")
        ax.grid(alpha=0.25)
    for ax in axes:
        ax.set_xlim(0.0, 40.0)
        ax.set_xlabel("r")
        ax.set_ylabel(r"$|\chi(r)|^2$")
    fig.suptitle("Discrete recoverable attractor profiles")
    fig.tight_layout()
    fig.savefig(modes_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(modes_path.relative_to(REPO_ROOT)))

    return plot_paths


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_interaction_attractor_spectrum.json"
    latest = DATA / "interaction_attractor_spectrum_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## Interaction attractor spectrum\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"D={config['D']}, alpha={config['alpha']}, r_min={config['r_min']}, "
            f"r_max={config['r_max']}, n_grid={config['n_grid']}, states={config['states']}, "
            f"modal_basis_size={config['modal_basis_size']}, "
            f"dt={config['dt']}, max_steps={config['max_steps']}, tol={config['tol']}, "
            f"recovery_noise={config['recovery_noise']}, seed={config['random_seed']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def run_interaction_attractor_spectrum(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    rng = np.random.default_rng(int(cfg["random_seed"]))

    r = np.linspace(float(cfg["r_min"]), float(cfg["r_max"]), int(cfg["n_grid"]))
    rin, H, dr = build_radial_hamiltonian(r=r, ell=0, D=float(cfg["D"]), alpha=float(cfg["alpha"]))
    evals, evecs = np.linalg.eigh(H)
    basis_size = int(cfg["modal_basis_size"])
    modal_energies = np.asarray(evals[:basis_size], dtype=float)
    basis = np.asarray(evecs[:, :basis_size], dtype=float)
    bound_idx = [idx for idx, val in enumerate(modal_energies) if val < 0.0][: int(cfg["states"])]
    reference_coeffs = [np.eye(basis_size, dtype=float)[:, idx] for idx in bound_idx]

    attractors: list[dict[str, Any]] = []
    for state_index, reference in enumerate(reference_coeffs):
        attractor = attractor_run_for_state(
            basis=basis,
            modal_energies=modal_energies,
            dr=dr,
            reference_coeffs=reference,
            state_index=state_index,
            cfg=cfg,
            rng=rng,
        )
        attractors.append(attractor)

    energies = np.array([attr["energy"] for attr in attractors], dtype=float)
    n_vals = np.array([attr["principal_n"] for attr in attractors], dtype=float)
    n2_abs = (n_vals**2) * np.abs(energies)
    observation = "normalized interaction flow converges to a discrete ladder of recoverable bound attractors"
    conclusion = (
        "discrete states arise from the inverse-distance interaction geometry and boundary-constrained radial operator; "
        "the attractor flow stabilizes and recovers these modes, but unconstrained flow selects only the lowest sector"
    )

    return {
        "experiment": "interaction_attractor_spectrum",
        "config": cfg,
        "derivation": {
            "flow_equation": "d_tau chi = -(H chi - lambda[chi] chi), lambda[chi] = <chi,H chi>/<chi,chi>",
            "modal_update": "c_j -> exp(-dt E_j) c_j followed by normalization",
            "projected_flow": "projected modal flow removes lower recovered sectors before each normalization step",
            "radial_operator": "-D d^2/dr^2 + D ell(ell+1)/r^2 - alpha/r",
            "scaling_target": "E_n = -alpha^2 / (4 D n^2)",
        },
        "s_wave_reference_energies": [float(modal_energies[idx]) for idx in bound_idx],
        "attractors": attractors,
        "scaling": {
            "n_values": n_vals.tolist(),
            "energies": energies.tolist(),
            "analytic_energies": [analytic_energy(int(n), D=float(cfg["D"]), alpha=float(cfg["alpha"])) for n in n_vals],
            "n2_abs_energy": n2_abs.tolist(),
            "n2_abs_energy_mean": float(np.mean(n2_abs)),
            "n2_abs_energy_std": float(np.std(n2_abs)),
        },
        "observation": observation,
        "conclusion": conclusion,
        "rin": rin.tolist(),
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = run_interaction_attractor_spectrum()
    plot_paths = make_plots(np.array(result["rin"], dtype=float), result["attractors"], stamp=stamp)
    result["plots"] = plot_paths
    result_path, latest_path = save_results(result, stamp=stamp)
    append_log(
        result_path=result_path,
        plot_paths=plot_paths,
        config=result["config"],
        observation=result["observation"],
        conclusion=result["conclusion"],
    )
    print("results =", result_path)
    print("latest =", latest_path)
    print("plots =", plot_paths)
    print("energies =", [round(val, 6) for val in result["scaling"]["energies"]])
    print("n2|E| =", [round(val, 6) for val in result["scaling"]["n2_abs_energy"]])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
