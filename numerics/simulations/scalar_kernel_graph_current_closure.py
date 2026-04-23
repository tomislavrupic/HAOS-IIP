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
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Current_Closure_v1.md"
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
        "max_median_relative_error": 0.25,
        "max_p90_relative_error": 0.45,
        "max_refinement_profile_drift": 0.15,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_current_closure"), dict):
            merged.update({k: v for k, v in raw["scalar_kernel_graph_current_closure"].items() if k != "thresholds"})
            if isinstance(raw["scalar_kernel_graph_current_closure"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_current_closure"]["thresholds"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def build_clean_shell_radii(n_side: int, shell_r_max: float) -> tuple[np.ndarray, np.ndarray, int, float]:
    points, _, source_idx, h = build_cubic_points(n_side)
    source = points[source_idx]
    radii = np.linalg.norm(points - source, axis=1)
    labels = np.round(radii / max(h, 1.0e-12), 8) * h
    shell_radii: list[float] = []
    for label in np.unique(labels):
        mask = np.isclose(labels, label)
        mean_r = float(np.mean(radii[mask]))
        if mean_r <= 2.0 * h or mean_r > float(shell_r_max):
            continue
        shell_radii.append(mean_r)
    return points, np.asarray(sorted(shell_radii), dtype=float), source_idx, h


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
    points, shell_radii, source_idx, h = build_clean_shell_radii(n_side=n_side, shell_r_max=shell_r_max)
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

    source = points[source_idx]
    bulk_radii = np.linalg.norm(bulk_points - source, axis=1)
    nearest_shell = np.argmin(np.abs(bulk_radii[:, None] - shell_radii[None, :]), axis=1)

    shell_masks: list[np.ndarray] = []
    shell_metric_rr: list[float] = []
    kept_radii: list[float] = []
    for idx, shell_r in enumerate(shell_radii):
        mask_shell = nearest_shell == idx
        if not np.any(mask_shell):
            continue
        dirs = bulk_points[mask_shell] - source
        shell_norm = np.linalg.norm(dirs, axis=1)
        dirs = dirs / np.maximum(shell_norm[:, None], 1.0e-18)
        g_rr = float(np.mean(np.einsum("ni,nij,nj->n", dirs, bulk_metric[mask_shell], dirs)))
        shell_masks.append(mask_shell)
        shell_metric_rr.append(g_rr)
        kept_radii.append(float(shell_r))

    shell_radii = np.asarray(kept_radii, dtype=float)
    shell_metric_rr = np.asarray(shell_metric_rr, dtype=float)
    shell_u = np.stack([np.mean(bulk_trajectory[:, mask_shell], axis=1) for mask_shell in shell_masks], axis=1)
    cumulative_mass = np.stack([np.sum(bulk_trajectory[:, bulk_radii <= shell_r], axis=1) for shell_r in shell_radii], axis=1)

    dMdt = np.gradient(cumulative_mass, times, axis=0)
    dudr = np.gradient(shell_u, shell_radii, axis=1)

    area_metric = 4.0 * np.pi * (shell_radii[None, :] ** 2) * shell_metric_rr[None, :]
    constitutive_basis = -(area_metric * dudr)
    empirical_current = -dMdt

    time_mask = (times >= float(fit_window[0])) & (times <= float(fit_window[1]))
    basis_fit = constitutive_basis[time_mask].reshape(-1)
    current_fit = empirical_current[time_mask].reshape(-1)
    conductivity_fit = float(np.dot(basis_fit, current_fit) / max(np.dot(basis_fit, basis_fit), 1.0e-18))

    predicted_current = conductivity_fit * constitutive_basis[time_mask]
    denom = np.maximum(np.abs(empirical_current[time_mask]), 1.0e-12)
    relative_error = np.abs(empirical_current[time_mask] - predicted_current) / denom

    mean_empirical_profile = np.mean(np.abs(empirical_current[time_mask]), axis=0)
    normalized_profile = mean_empirical_profile / max(float(np.mean(mean_empirical_profile)), 1.0e-18)

    return {
        "label": f"clean|n={n_side}",
        "n_side": int(n_side),
        "nodes": int(points.shape[0]),
        "lattice_spacing": float(h),
        "epsilon_coeff": float(epsilon_coeff),
        "epsilon_k": float(epsilon_k),
        "heat_diffusivity": float(heat["effective_diffusivity"]),
        "conductivity_fit": conductivity_fit,
        "conductivity_over_diffusivity": float(conductivity_fit / max(float(heat["effective_diffusivity"]), 1.0e-18)),
        "median_relative_error": float(np.nanmedian(relative_error)),
        "mean_relative_error": float(np.nanmean(relative_error)),
        "p90_relative_error": float(np.nanquantile(relative_error, 0.9)),
        "shell_radii": shell_radii.tolist(),
        "shell_metric_rr": shell_metric_rr.tolist(),
        "normalized_current_profile": normalized_profile.tolist(),
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        case["median_relative_error"] <= float(thresholds["max_median_relative_error"])
        and case["p90_relative_error"] <= float(thresholds["max_p90_relative_error"])
    )


def profile_drift(case_a: dict[str, Any], case_b: dict[str, Any]) -> float:
    profile_a = np.asarray(case_a["normalized_current_profile"], dtype=float)
    profile_b = np.asarray(case_b["normalized_current_profile"], dtype=float)
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
    weakest_case = max(cases, key=lambda item: (item["median_relative_error"], item["p90_relative_error"]))

    if all_cases_pass and drift_pass:
        observation = (
            "on the clean scalar carrier, the same local metric field now supports a bounded transient response/current closure: "
            "a single constitutive coefficient matches cumulative shell-mass flow to the shell-gradient current with small residual and stable normalized profile under refinement"
        )
        conclusion = (
            f"all {len(cases)} clean refinement cases pass, the maximum normalized profile drift is {max_profile_drift:.4f}, "
            f"and the weakest passing case is `{weakest_case['label']}` with median relative error {weakest_case['median_relative_error']:.4f} and p90 error {weakest_case['p90_relative_error']:.4f}; "
            "this supports a bounded transient current-closure read on the same scalar carrier"
        )
    else:
        failing_cases = [case["label"] for case in cases if not case["pass"]]
        observation = (
            "the scalar carrier and local metric field now support a bounded inverse-square response law, but the stronger transient response/current closure remains open on the clean refinement scan"
        )
        conclusion = (
            f"open current-closure cases remain: failing_cases={failing_cases}, max_profile_drift={max_profile_drift:.4f}; "
            f"the weakest case is `{weakest_case['label']}` with median relative error {weakest_case['median_relative_error']:.4f}, "
            f"p90 error {weakest_case['p90_relative_error']:.4f}, and conductivity/diffusivity ratio {weakest_case['conductivity_over_diffusivity']:.4f}; "
            "the correct boundary therefore stays at the scalar carrier, local metric field, and shell-native inverse-square response law without yet promoting a positive transient current-closure claim"
        )

    return {
        "experiment": "scalar_kernel_graph_current_closure",
        "config": cfg,
        "cases": cases,
        "profile_drifts": drifts,
        "max_profile_drift": float(max_profile_drift),
        "observation": observation,
        "conclusion": conclusion,
    }


def make_plots(result: dict[str, Any], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_summary.png"
    n_vals = [case["n_side"] for case in result["cases"]]
    median_err = [case["median_relative_error"] for case in result["cases"]]
    p90_err = [case["p90_relative_error"] for case in result["cases"]]
    ratio = [case["conductivity_over_diffusivity"] for case in result["cases"]]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(n_vals, median_err, marker="o", label="median relative error")
    axes[0].plot(n_vals, p90_err, marker="s", label="p90 relative error")
    axes[0].set_title("Transient current-closure error")
    axes[0].set_xlabel("cubic side n")
    axes[0].set_ylabel("relative error")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].plot(n_vals, ratio, marker="o", label=r"$\kappa_{\mathrm{fit}} / D_{\mathrm{eff}}$")
    axes[1].set_title("Constitutive coefficient mismatch")
    axes[1].set_xlabel("cubic side n")
    axes[1].set_ylabel("ratio")
    axes[1].grid(alpha=0.25)
    axes[1].legend()

    fig.savefig(summary_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(summary_path.relative_to(REPO_ROOT)))

    detail_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_profiles.png"
    rep = result["cases"][-1]
    radii = np.asarray(rep["shell_radii"], dtype=float)
    profile = np.asarray(rep["normalized_current_profile"], dtype=float)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(radii, profile, marker="o", label="normalized empirical current profile")
    ax.axhline(1.0, color="k", linestyle="--", linewidth=1.0, label="flat profile")
    ax.set_title(f"Current profile shape (n={rep['n_side']})")
    ax.set_xlabel("shell radius")
    ax.set_ylabel("normalized shell current")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(detail_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(detail_path.relative_to(REPO_ROOT)))

    return plot_paths


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_current_closure.json"
    latest = DATA / "scalar_kernel_graph_current_closure_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    rows = "\n".join(
        f"| `{case['label']}` | `{case['heat_diffusivity']:.4f}` | `{case['conductivity_fit']:.4f}` | `{case['conductivity_over_diffusivity']:.4f}` | `{case['median_relative_error']:.4f}` | `{case['p90_relative_error']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["cases"]
    )
    drift_rows = "\n".join(f"| `{row['label']}` | `{row['drift']:.4f}` |" for row in result["profile_drifts"])
    note = f"""# Scalar Kernel Graph Current Closure v1

## Purpose

Push the scalar carrier from the static `52.2` inverse-square response law to the next stronger question:

> does the same scalar carrier support a transient response/current closure on the clean refinement line?

The carrier is unchanged:

- same local scalar kernel graph
- same `51.4` local metric field
- same bounded Euclidean shell reconstruction

The new ingredient is transient heat-flow current.

## Construction

For each clean cubic refinement:

1. evolve the point-source heat trajectory on the same scalar operator;
2. build shell-mean profiles on the clean shell scaffold;
3. infer empirical shell current from cumulative mass flow,
4. compare it against a constitutive shell current built from the same metric field and shell gradient:

$$
I_{{const}}(r,t) = - \\kappa \\; 4\\pi r^2 g_{{rr}}(r) \\; \\partial_r u(r,t)
$$

where `\\kappa` is fitted per refinement on the bounded time window.

This tests whether one transient current law closes on the same scalar carrier, not just whether the static Green response is inverse-square.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_current_closure.py`
- config: `config.json -> scalar_kernel_graph_current_closure`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_current_closure_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Clean refinement scan

| case | `D_eff` | `kappa_fit` | `kappa_fit / D_eff` | median rel. error | p90 rel. error | status |
| --- | --- | --- | --- | --- | --- | --- |
{rows}

## Normalized profile drift

| pair | drift |
| --- | --- |
{drift_rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Current boundary

If this note is read negatively, the correct bounded statement is:

> the same scalar carrier now supports the static `52.2` inverse-square response law, but transient constitutive current closure remains open even on the clean refinement scan.

So the next honest move after this note is not a stronger claim. It is:

1. improve the transient shell-current reconstruction,
2. understand the scaling of the fitted constitutive coefficient,
3. then only afterward widen back toward disorder or kernel-family robustness.
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph current closure\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, epsilon_coeff={config['epsilon_coeff']}, "
            f"cutoff_factor={config['cutoff_factor']}, fit_window={config['fit_window']}, "
            f"metric_radius_factor={config['metric_radius_factor']}\n"
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
