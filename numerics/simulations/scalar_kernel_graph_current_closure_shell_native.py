#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from scalar_kernel_graph_geometry_bridge import build_cubic_points, heat_metrics, heat_trajectory
from scalar_kernel_graph_geometry_robustness import build_point_cloud_operator
from scalar_kernel_graph_metric_field import bulk_mask, coarse_grain_metric_field, local_metric_field

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Current_Closure_Shell_Native_v1.md"
for path in (DATA, PLOTS):
    path.mkdir(exist_ok=True)

MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    "n_sides": [13, 15, 17],
    "epsilon_coeff": 0.5,
    "cutoff_factor": 2.5,
    "heat_time_max": 0.03,
    "heat_time_steps": 61,
    "fit_window": [0.004, 0.02],
    "bulk_margin_layers": 2,
    "metric_radius_factor": 2.5,
    "shell_r_max": 0.45,
    "thresholds": {
        "max_median_relative_error": 0.10,
        "max_p90_relative_error": 0.30,
        "max_diffusivity_match_gap": 0.15,
        "max_shell_kappa_cv": 0.13,
        "max_refinement_profile_drift": 0.10,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_current_closure_shell_native"), dict):
            merged.update(
                {k: v for k, v in raw["scalar_kernel_graph_current_closure_shell_native"].items() if k != "thresholds"}
            )
            if isinstance(raw["scalar_kernel_graph_current_closure_shell_native"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_current_closure_shell_native"]["thresholds"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def build_bulk_shells(
    bulk_points: np.ndarray,
    bulk_metric: np.ndarray,
    source_point: np.ndarray,
    h: float,
    shell_r_max: float,
) -> tuple[np.ndarray, list[np.ndarray], np.ndarray, np.ndarray]:
    bulk_radii = np.linalg.norm(bulk_points - source_point, axis=1)
    labels = np.round(bulk_radii / max(h, 1.0e-12), 8) * h
    shell_radii: list[float] = []
    shell_masks: list[np.ndarray] = []
    shell_metric_rr: list[float] = []
    shell_counts: list[int] = []
    for label in np.unique(labels):
        mask_shell = np.isclose(labels, label)
        if not np.any(mask_shell):
            continue
        mean_r = float(np.mean(bulk_radii[mask_shell]))
        if mean_r <= 2.0 * h or mean_r > float(shell_r_max):
            continue
        radial_dirs = bulk_points[mask_shell] - source_point
        radial_norm = np.linalg.norm(radial_dirs, axis=1)
        radial_dirs = radial_dirs / np.maximum(radial_norm[:, None], 1.0e-18)
        shell_radii.append(mean_r)
        shell_masks.append(mask_shell)
        shell_metric_rr.append(
            float(np.mean(np.einsum("ni,nij,nj->n", radial_dirs, bulk_metric[mask_shell], radial_dirs)))
        )
        shell_counts.append(int(np.sum(mask_shell)))
    return (
        np.asarray(shell_radii, dtype=float),
        shell_masks,
        np.asarray(shell_metric_rr, dtype=float),
        np.asarray(shell_counts, dtype=float),
    )


def evaluate_case(
    n_side: int,
    epsilon_coeff: float,
    cutoff_factor: float,
    heat_time_max: float,
    heat_time_steps: int,
    fit_window: list[float],
    bulk_margin_layers: int,
    metric_radius_factor: float,
    shell_r_max: float,
) -> dict[str, Any]:
    points, _, source_idx, h = build_cubic_points(n_side)
    delta, epsilon_k = build_point_cloud_operator(
        points=points,
        h=h,
        epsilon_coeff=epsilon_coeff,
        cutoff_factor=cutoff_factor,
        kernel_family="gaussian_local",
    )
    times, trajectory = heat_trajectory(
        delta=delta,
        source_idx=source_idx,
        time_max=heat_time_max,
        time_steps=heat_time_steps,
    )
    heat = heat_metrics(
        points=points,
        source_idx=source_idx,
        times=times,
        trajectory=trajectory,
        fit_window=fit_window,
    )
    raw_metric = local_metric_field(points=points, delta=delta)
    mask = bulk_mask(points=points, h=h, bulk_margin_layers=bulk_margin_layers)
    bulk_points = points[mask]
    bulk_trajectory = trajectory[:, mask]
    bulk_metric = coarse_grain_metric_field(points=bulk_points, field=raw_metric[mask], radius=float(metric_radius_factor) * h)
    shell_radii, shell_masks, shell_metric_rr, shell_counts = build_bulk_shells(
        bulk_points=bulk_points,
        bulk_metric=bulk_metric,
        source_point=points[source_idx],
        h=h,
        shell_r_max=shell_r_max,
    )
    shell_mass = np.stack([np.sum(bulk_trajectory[:, mask_shell], axis=1) for mask_shell in shell_masks], axis=1)
    shell_density = shell_mass / (shell_counts[None, :] * (h**3))
    bulk_radii = np.linalg.norm(bulk_points - points[source_idx], axis=1)
    cumulative_mass = np.stack([np.sum(bulk_trajectory[:, bulk_radii <= shell_r], axis=1) for shell_r in shell_radii], axis=1)
    empirical_current = -np.gradient(cumulative_mass, times, axis=0)
    density_gradient = np.gradient(shell_density, shell_radii, axis=1)
    shell_area = 4.0 * np.pi * (shell_radii[None, :] ** 2) * shell_metric_rr[None, :]
    constitutive_basis = -(shell_area * density_gradient)

    time_mask = (times >= float(fit_window[0])) & (times <= float(fit_window[1]))
    basis_fit = constitutive_basis[time_mask].reshape(-1)
    current_fit = empirical_current[time_mask].reshape(-1)
    conductivity_fit = float(np.dot(basis_fit, current_fit) / max(np.dot(basis_fit, basis_fit), 1.0e-18))

    predicted_current = conductivity_fit * constitutive_basis[time_mask]
    denom = np.maximum(np.abs(empirical_current[time_mask]), 1.0e-12)
    relative_error = np.abs(empirical_current[time_mask] - predicted_current) / denom

    shell_ratio = np.divide(
        empirical_current[time_mask],
        constitutive_basis[time_mask],
        out=np.full_like(empirical_current[time_mask], np.nan),
        where=np.abs(constitutive_basis[time_mask]) > 1.0e-12,
    )
    shell_kappa_profile = np.nanmedian(shell_ratio, axis=0)
    shell_kappa_scale = max(float(np.nanmean(np.abs(shell_kappa_profile))), 1.0e-18)
    normalized_shell_kappa_profile = shell_kappa_profile / shell_kappa_scale

    heat_diffusivity = float(heat["effective_diffusivity"])
    conductivity_over_diffusivity = float(conductivity_fit / max(heat_diffusivity, 1.0e-18))
    return {
        "label": f"clean|n={n_side}",
        "n_side": int(n_side),
        "nodes": int(points.shape[0]),
        "lattice_spacing": float(h),
        "epsilon_coeff": float(epsilon_coeff),
        "epsilon_k": float(epsilon_k),
        "heat_diffusivity": heat_diffusivity,
        "conductivity_fit": conductivity_fit,
        "conductivity_over_diffusivity": conductivity_over_diffusivity,
        "diffusivity_match_gap": float(abs(conductivity_over_diffusivity - 1.0)),
        "median_relative_error": float(np.nanmedian(relative_error)),
        "mean_relative_error": float(np.nanmean(relative_error)),
        "p90_relative_error": float(np.nanquantile(relative_error, 0.9)),
        "shell_kappa_cv": float(np.nanstd(shell_kappa_profile) / shell_kappa_scale),
        "shell_radii": shell_radii.tolist(),
        "shell_metric_rr": shell_metric_rr.tolist(),
        "shell_kappa_profile": shell_kappa_profile.tolist(),
        "normalized_shell_kappa_profile": normalized_shell_kappa_profile.tolist(),
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        case["median_relative_error"] <= float(thresholds["max_median_relative_error"])
        and case["p90_relative_error"] <= float(thresholds["max_p90_relative_error"])
        and case["diffusivity_match_gap"] <= float(thresholds["max_diffusivity_match_gap"])
        and case["shell_kappa_cv"] <= float(thresholds["max_shell_kappa_cv"])
    )


def profile_drift(case_a: dict[str, Any], case_b: dict[str, Any]) -> float:
    profile_a = np.asarray(case_a["normalized_shell_kappa_profile"], dtype=float)
    profile_b = np.asarray(case_b["normalized_shell_kappa_profile"], dtype=float)
    min_len = min(len(profile_a), len(profile_b))
    return float(np.linalg.norm(profile_a[:min_len] - profile_b[:min_len]) / np.sqrt(min_len))


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    thresholds = cfg["thresholds"]
    cases = [
        evaluate_case(
            n_side=int(n_side),
            epsilon_coeff=float(cfg["epsilon_coeff"]),
            cutoff_factor=float(cfg["cutoff_factor"]),
            heat_time_max=float(cfg["heat_time_max"]),
            heat_time_steps=int(cfg["heat_time_steps"]),
            fit_window=[float(v) for v in cfg["fit_window"]],
            bulk_margin_layers=int(cfg["bulk_margin_layers"]),
            metric_radius_factor=float(cfg["metric_radius_factor"]),
            shell_r_max=float(cfg["shell_r_max"]),
        )
        for n_side in [int(v) for v in cfg["n_sides"]]
    ]
    for case in cases:
        case["pass"] = case_passes(case, thresholds)

    drifts = []
    for left, right in zip(cases[:-1], cases[1:]):
        drifts.append({"label": f"{left['n_side']}->{right['n_side']}", "drift": profile_drift(left, right)})

    max_profile_drift = max((row["drift"] for row in drifts), default=0.0)
    all_cases_pass = all(case["pass"] for case in cases)
    drift_pass = max_profile_drift <= float(thresholds["max_refinement_profile_drift"])
    weakest_case = max(
        cases,
        key=lambda item: (
            item["median_relative_error"],
            item["p90_relative_error"],
            item["diffusivity_match_gap"],
            item["shell_kappa_cv"],
        ),
    )

    if all_cases_pass and drift_pass:
        observation = (
            "on the clean scalar carrier, a shell-native transient current reconstruction now closes onto one bounded constitutive family: "
            "exact bulk shells plus density-normalized gradients recover a stable effective conductivity close to the heat diffusivity"
        )
        conclusion = (
            f"all {len(cases)} clean refinement cases pass, the maximum normalized shell-kappa drift is {max_profile_drift:.4f}, "
            f"and the weakest passing case is `{weakest_case['label']}` with kappa/D_eff {weakest_case['conductivity_over_diffusivity']:.4f}, "
            f"median relative error {weakest_case['median_relative_error']:.4f}, p90 error {weakest_case['p90_relative_error']:.4f}, "
            f"and shell-kappa CV {weakest_case['shell_kappa_cv']:.4f}; this supports a bounded shell-native transient current-closure read on the same scalar carrier"
        )
    else:
        failing_cases = [case["label"] for case in cases if not case["pass"]]
        observation = (
            "the shell-native transient current reconstruction improves the clean scalar-carrier read, but one bounded constitutive family is not yet stable enough under refinement"
        )
        conclusion = (
            f"open shell-native current-closure cases remain: failing_cases={failing_cases}, max_profile_drift={max_profile_drift:.4f}; "
            f"the weakest case is `{weakest_case['label']}` with kappa/D_eff {weakest_case['conductivity_over_diffusivity']:.4f}, "
            f"median relative error {weakest_case['median_relative_error']:.4f}, p90 error {weakest_case['p90_relative_error']:.4f}, "
            f"and shell-kappa CV {weakest_case['shell_kappa_cv']:.4f}; the correct boundary therefore stays short of a positive transient current-closure claim"
        )

    return {
        "experiment": "scalar_kernel_graph_current_closure_shell_native",
        "config": cfg,
        "cases": cases,
        "profile_drifts": drifts,
        "max_profile_drift": float(max_profile_drift),
        "observation": observation,
        "conclusion": conclusion,
    }


def make_plots(result: dict[str, Any], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_shell_native_summary.png"
    n_vals = [case["n_side"] for case in result["cases"]]
    median_err = [case["median_relative_error"] for case in result["cases"]]
    p90_err = [case["p90_relative_error"] for case in result["cases"]]
    ratio = [case["conductivity_over_diffusivity"] for case in result["cases"]]
    shell_cv = [case["shell_kappa_cv"] for case in result["cases"]]
    drift_vals = [row["drift"] for row in result["profile_drifts"]]
    drift_labels = [row["label"] for row in result["profile_drifts"]]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes[0, 0].plot(n_vals, median_err, marker="o", label="median relative error")
    axes[0, 0].plot(n_vals, p90_err, marker="s", label="p90 relative error")
    axes[0, 0].set_title("Shell-native current-closure error")
    axes[0, 0].set_xlabel("cubic side n")
    axes[0, 0].set_ylabel("relative error")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].plot(n_vals, ratio, marker="o", label=r"$\kappa_{\mathrm{fit}} / D_{\mathrm{eff}}$")
    axes[0, 1].axhline(1.0, color="k", linestyle="--", linewidth=1.0)
    axes[0, 1].set_title("Diffusivity match")
    axes[0, 1].set_xlabel("cubic side n")
    axes[0, 1].set_ylabel("ratio")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].plot(n_vals, shell_cv, marker="o", label="shell-kappa CV")
    axes[1, 0].set_title("Shellwise constitutive stability")
    axes[1, 0].set_xlabel("cubic side n")
    axes[1, 0].set_ylabel("CV")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].plot(range(len(drift_vals)), drift_vals, marker="o")
    axes[1, 1].set_xticks(range(len(drift_vals)), labels=drift_labels)
    axes[1, 1].set_title("Normalized shell-kappa drift")
    axes[1, 1].set_ylabel("drift")
    axes[1, 1].grid(alpha=0.25)

    fig.savefig(summary_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(summary_path.relative_to(REPO_ROOT)))

    detail_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_shell_native_profiles.png"
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for case in result["cases"]:
        radii = np.asarray(case["shell_radii"], dtype=float)
        profile = np.asarray(case["normalized_shell_kappa_profile"], dtype=float)
        ax.plot(radii, profile, marker="o", label=f"n={case['n_side']}")
    ax.axhline(1.0, color="k", linestyle="--", linewidth=1.0, label="flat constitutive profile")
    ax.set_title("Normalized shell-kappa profiles")
    ax.set_xlabel("shell radius")
    ax.set_ylabel("normalized kappa(r)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(detail_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(detail_path.relative_to(REPO_ROOT)))

    return plot_paths


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_current_closure_shell_native.json"
    latest = DATA / "scalar_kernel_graph_current_closure_shell_native_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    rows = "\n".join(
        f"| `{case['label']}` | `{case['heat_diffusivity']:.4f}` | `{case['conductivity_fit']:.4f}` | `{case['conductivity_over_diffusivity']:.4f}` | `{case['median_relative_error']:.4f}` | `{case['p90_relative_error']:.4f}` | `{case['shell_kappa_cv']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["cases"]
    )
    drift_rows = "\n".join(f"| `{row['label']}` | `{row['drift']:.4f}` |" for row in result["profile_drifts"])
    note = f"""# Scalar Kernel Graph Current Closure Shell-Native v1

## Purpose

Retest the open transient current-closure question on the same scalar carrier without changing the carrier, the `51.4` local metric field, or the clean refinement line.

The only change is the reconstruction:

- keep the same local scalar kernel graph
- keep the same coarse local metric field
- replace the nearest-shell profile fit with exact bulk shells on the clean scaffold
- read the constitutive law against shell density rather than raw node mass

This is a shell-native reconstruction audit, not a new substrate claim.

## Construction

For each clean cubic refinement:

1. evolve the same point-source heat trajectory;
2. build exact bulk shells from the clean radius labels;
3. convert shell mass to shell density via `rho = m_shell / (count * h^3)`;
4. compare empirical cumulative shell current against

$$
I_{{const}}(r,t) = - \\kappa \\; 4\\pi r^2 g_{{rr}}(r) \\; \\partial_r \\rho(r,t).
$$

The closure diagnostic is no longer the raw current profile alone. It is whether one bounded `\\kappa` family stays close to the heat diffusivity and remains stable across radius and refinement.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_current_closure_shell_native.py`
- config: `config.json -> scalar_kernel_graph_current_closure_shell_native`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_current_closure_shell_native_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Clean refinement scan

| case | `D_eff` | `kappa_fit` | `kappa_fit / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
{rows}

## Normalized shell-kappa drift

| pair | drift |
| --- | --- |
{drift_rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Current boundary

If this note is read positively, the correct bounded statement is:

> on the clean scalar carrier, a shell-native transient current reconstruction closes onto one bounded constitutive family whose fitted conductivity stays close to the heat diffusivity and whose shellwise profile remains refinement-stable.

This note still does **not** claim:

- disorder or kernel-family universality
- conserved-current closure
- curvature extraction
- spacetime or ontology
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph current closure shell-native\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, epsilon_coeff={config['epsilon_coeff']}, "
            f"cutoff_factor={config['cutoff_factor']}, fit_window={config['fit_window']}, "
            f"metric_radius_factor={config['metric_radius_factor']}, shell_r_max={config['shell_r_max']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = make_result()
    plots = make_plots(result, stamp)
    result_path, _ = save_results(result, stamp)
    write_note(result, result_path, plots)
    append_log(result_path, plots, result["config"], result["observation"], result["conclusion"])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
