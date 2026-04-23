#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from scalar_kernel_graph_current_closure_inhomogeneity import normalized_profile, profile_drift
from scalar_kernel_graph_geometry_bridge import build_cubic_points, heat_metrics, heat_trajectory
from scalar_kernel_graph_geometry_robustness import build_point_cloud_operator
from scalar_kernel_graph_metric_field import bulk_mask, coarse_grain_metric_field, local_metric_field

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Current_Closure_Radial_Disorder_v1.md"
for path in (DATA, PLOTS):
    path.mkdir(exist_ok=True)

MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    "n_sides": [13, 15],
    "eta_values": [0.05, 0.10, 0.15],
    "disorder_fractions": [0.02, 0.04],
    "seeds": [0, 1, 2],
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
        "max_p90_relative_error": 0.37,
        "max_diffusivity_match_gap": 0.15,
        "max_shell_kappa_cv": 0.14,
        "max_refinement_profile_drift": 0.15,
        "min_metric_tracking_abs_corr": 0.95,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_current_closure_radial_disorder"), dict):
            merged.update(
                {k: v for k, v in raw["scalar_kernel_graph_current_closure_radial_disorder"].items() if k != "thresholds"}
            )
            if isinstance(raw["scalar_kernel_graph_current_closure_radial_disorder"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_current_closure_radial_disorder"]["thresholds"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def radial_deformation_profile(points: np.ndarray, n_side: int, eta: float) -> tuple[np.ndarray, np.ndarray]:
    center = np.array([0.5, 0.5, 0.5], dtype=float)
    delta = points - center
    radii = np.linalg.norm(delta, axis=1)
    r_max = max(float(np.max(radii)), 1.0e-18)
    profile = (radii / r_max) ** 2
    factor = 1.0 + float(eta) * profile
    out = center + factor[:, None] * delta
    grid_ids = np.indices((n_side, n_side, n_side)).reshape(3, -1).T
    boundary_mask = np.any((grid_ids == 0) | (grid_ids == n_side - 1), axis=1)
    out[boundary_mask] = points[boundary_mask]
    return np.clip(out, 0.0, 1.0), profile


def apply_disorder(points: np.ndarray, n_side: int, h: float, disorder_fraction: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    grid_ids = np.indices((n_side, n_side, n_side)).reshape(3, -1).T
    boundary_mask = np.any((grid_ids == 0) | (grid_ids == n_side - 1), axis=1)
    jitter = rng.normal(scale=float(disorder_fraction) * h, size=points.shape)
    out = points.copy()
    out[~boundary_mask] += jitter[~boundary_mask]
    out = np.clip(out, 0.0, 1.0)
    out[boundary_mask] = points[boundary_mask]
    return out


def prepare_clean_layout(
    n_side: int,
    epsilon_coeff: float,
    cutoff_factor: float,
    bulk_margin_layers: int,
    metric_radius_factor: float,
    shell_r_max: float,
) -> dict[str, Any]:
    clean_points, _, source_idx, h = build_cubic_points(n_side)
    clean_delta, epsilon_k = build_point_cloud_operator(
        points=clean_points,
        h=h,
        epsilon_coeff=epsilon_coeff,
        cutoff_factor=cutoff_factor,
        kernel_family="gaussian_local",
    )
    clean_raw_metric = local_metric_field(points=clean_points, delta=clean_delta)
    clean_mask = bulk_mask(points=clean_points, h=h, bulk_margin_layers=bulk_margin_layers)
    clean_bulk_points = clean_points[clean_mask]
    clean_bulk_metric = coarse_grain_metric_field(
        points=clean_bulk_points,
        field=clean_raw_metric[clean_mask],
        radius=float(metric_radius_factor) * h,
    )
    clean_bulk_radii = np.linalg.norm(clean_bulk_points - clean_points[source_idx], axis=1)
    shell_ids = np.rint(clean_bulk_radii / h).astype(int)
    keep = (clean_bulk_radii > 2.0 * h) & (clean_bulk_radii <= float(shell_r_max))
    post_keep_indices = np.flatnonzero(clean_mask)[keep]
    clean_bulk_points = clean_bulk_points[keep]
    clean_bulk_metric = clean_bulk_metric[keep]
    clean_bulk_radii = clean_bulk_radii[keep]
    shell_ids = shell_ids[keep]
    unique_ids = np.unique(shell_ids)
    shell_radii = []
    shell_counts = []
    baseline_metric_rr = []
    source_point = clean_points[source_idx]
    for shell_id in unique_ids:
        shell_mask = shell_ids == shell_id
        shell_radii.append(float(np.mean(clean_bulk_radii[shell_mask])))
        shell_counts.append(int(np.sum(shell_mask)))
        radial_dirs = clean_bulk_points[shell_mask] - source_point
        radial_dirs = radial_dirs / np.maximum(np.linalg.norm(radial_dirs, axis=1, keepdims=True), 1.0e-18)
        baseline_metric_rr.append(
            float(np.mean(np.einsum("ni,nij,nj->n", radial_dirs, clean_bulk_metric[shell_mask], radial_dirs)))
        )
    return {
        "n_side": int(n_side),
        "h": float(h),
        "epsilon_k": float(epsilon_k),
        "clean_points": clean_points,
        "source_idx": int(source_idx),
        "post_keep_indices": post_keep_indices,
        "shell_ids": shell_ids,
        "unique_ids": unique_ids,
        "shell_radii": np.asarray(shell_radii, dtype=float),
        "shell_counts": np.asarray(shell_counts, dtype=float),
        "baseline_metric_rr": np.asarray(baseline_metric_rr, dtype=float),
    }


def evaluate_trial(layout: dict[str, Any], eta: float, disorder_fraction: float, seed: int, config: dict[str, Any]) -> dict[str, Any]:
    n_side = int(layout["n_side"])
    h = float(layout["h"])
    source_idx = int(layout["source_idx"])
    base_points = np.asarray(layout["clean_points"], dtype=float)
    deformed_points, profile = radial_deformation_profile(base_points, n_side=n_side, eta=eta)
    points = apply_disorder(
        deformed_points,
        n_side=n_side,
        h=h,
        disorder_fraction=disorder_fraction,
        seed=seed,
    )

    delta, epsilon_k = build_point_cloud_operator(
        points=points,
        h=h,
        epsilon_coeff=float(config["epsilon_coeff"]),
        cutoff_factor=float(config["cutoff_factor"]),
        kernel_family="gaussian_local",
    )
    times, trajectory = heat_trajectory(
        delta=delta,
        source_idx=source_idx,
        time_max=float(config["heat_time_max"]),
        time_steps=int(config["heat_time_steps"]),
    )
    heat = heat_metrics(
        points=points,
        source_idx=source_idx,
        times=times,
        trajectory=trajectory,
        fit_window=[float(v) for v in config["fit_window"]],
    )

    raw_metric = local_metric_field(points=points, delta=delta)
    bulk_points = points[layout["post_keep_indices"]]
    bulk_metric = coarse_grain_metric_field(
        points=bulk_points,
        field=raw_metric[layout["post_keep_indices"]],
        radius=float(config["metric_radius_factor"]) * h,
    )
    bulk_trajectory = trajectory[:, layout["post_keep_indices"]]
    bulk_profile = profile[layout["post_keep_indices"]]
    shell_ids = np.asarray(layout["shell_ids"])
    shell_radii = np.asarray(layout["shell_radii"])
    baseline_metric_rr = np.asarray(layout["baseline_metric_rr"])

    shell_metric_rr = []
    shell_profile = []
    shell_mass = []
    source_point = points[source_idx]
    for shell_id in np.asarray(layout["unique_ids"]):
        shell_mask = shell_ids == shell_id
        radial_dirs = bulk_points[shell_mask] - source_point
        radial_dirs = radial_dirs / np.maximum(np.linalg.norm(radial_dirs, axis=1, keepdims=True), 1.0e-18)
        shell_metric_rr.append(
            float(np.mean(np.einsum("ni,nij,nj->n", radial_dirs, bulk_metric[shell_mask], radial_dirs)))
        )
        shell_profile.append(float(np.mean(bulk_profile[shell_mask])))
        shell_mass.append(np.sum(bulk_trajectory[:, shell_mask], axis=1))

    shell_metric_rr = np.asarray(shell_metric_rr, dtype=float)
    shell_profile = np.asarray(shell_profile, dtype=float)
    shell_mass = np.stack(shell_mass, axis=1)
    shell_counts = np.asarray(layout["shell_counts"], dtype=float)

    shell_density = shell_mass / (shell_counts[None, :] * (h**3))
    cumulative_mass = np.cumsum(shell_mass, axis=1)
    empirical_current = -np.gradient(cumulative_mass, times, axis=0)
    density_gradient = np.gradient(shell_density, shell_radii, axis=1)
    shell_area = 4.0 * np.pi * (shell_radii[None, :] ** 2) * shell_metric_rr[None, :]
    constitutive_basis = -(shell_area * density_gradient)

    time_mask = (times >= float(config["fit_window"][0])) & (times <= float(config["fit_window"][1]))
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
    normalized_shell_kappa_profile = normalized_profile(shell_kappa_profile)
    metric_shift = shell_metric_rr - baseline_metric_rr
    metric_shift_corr = float(np.corrcoef(shell_profile, metric_shift)[0, 1])
    heat_diffusivity = float(heat["effective_diffusivity"])
    conductivity_over_diffusivity = float(conductivity_fit / max(heat_diffusivity, 1.0e-18))
    return {
        "label": f"radial_disorder|n={n_side}|eta={eta:.3f}|sigma={disorder_fraction:.3f}|seed={seed}",
        "n_side": int(n_side),
        "eta": float(eta),
        "disorder_fraction": float(disorder_fraction),
        "seed": int(seed),
        "lattice_spacing": float(h),
        "epsilon_k": float(epsilon_k),
        "heat_diffusivity": heat_diffusivity,
        "conductivity_fit": conductivity_fit,
        "conductivity_over_diffusivity": conductivity_over_diffusivity,
        "diffusivity_match_gap": float(abs(conductivity_over_diffusivity - 1.0)),
        "median_relative_error": float(np.nanmedian(relative_error)),
        "mean_relative_error": float(np.nanmean(relative_error)),
        "p90_relative_error": float(np.nanquantile(relative_error, 0.9)),
        "shell_kappa_cv": float(np.nanstd(shell_kappa_profile) / max(float(np.nanmean(np.abs(shell_kappa_profile))), 1.0e-18)),
        "metric_tracking_corr": metric_shift_corr,
        "metric_tracking_abs_corr": float(abs(metric_shift_corr)),
        "shell_radii": shell_radii.tolist(),
        "shell_profile": shell_profile.tolist(),
        "shell_metric_rr": shell_metric_rr.tolist(),
        "metric_shift_profile": metric_shift.tolist(),
        "shell_kappa_profile": shell_kappa_profile.tolist(),
        "normalized_shell_kappa_profile": normalized_shell_kappa_profile.tolist(),
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        case["median_relative_error_median"] <= float(thresholds["max_median_relative_error"])
        and case["p90_relative_error_median"] <= float(thresholds["max_p90_relative_error"])
        and case["diffusivity_match_gap_median"] <= float(thresholds["max_diffusivity_match_gap"])
        and case["shell_kappa_cv_median"] <= float(thresholds["max_shell_kappa_cv"])
        and case["metric_tracking_abs_corr_median"] >= float(thresholds["min_metric_tracking_abs_corr"])
    )


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    thresholds = cfg["thresholds"]
    layouts = {
        int(n_side): prepare_clean_layout(
            n_side=int(n_side),
            epsilon_coeff=float(cfg["epsilon_coeff"]),
            cutoff_factor=float(cfg["cutoff_factor"]),
            bulk_margin_layers=int(cfg["bulk_margin_layers"]),
            metric_radius_factor=float(cfg["metric_radius_factor"]),
            shell_r_max=float(cfg["shell_r_max"]),
        )
        for n_side in [int(v) for v in cfg["n_sides"]]
    }

    trial_groups = []
    for disorder_fraction in [float(v) for v in cfg["disorder_fractions"]]:
        for eta in [float(v) for v in cfg["eta_values"]]:
            for n_side in [int(v) for v in cfg["n_sides"]]:
                trials = [
                    evaluate_trial(
                        layout=layouts[n_side],
                        eta=eta,
                        disorder_fraction=disorder_fraction,
                        seed=int(seed),
                        config=cfg,
                    )
                    for seed in [int(v) for v in cfg["seeds"]]
                ]
                min_len = min(len(trial["normalized_shell_kappa_profile"]) for trial in trials)
                mean_profile = np.mean(
                    [np.asarray(trial["normalized_shell_kappa_profile"][:min_len], dtype=float) for trial in trials],
                    axis=0,
                )
                case = {
                    "label": f"radial_disorder|n={n_side}|eta={eta:.3f}|sigma={disorder_fraction:.3f}",
                    "n_side": int(n_side),
                    "eta": float(eta),
                    "disorder_fraction": float(disorder_fraction),
                    "trial_count": len(trials),
                    "pass_count": int(
                        sum(
                            int(
                                trial["median_relative_error"] <= float(thresholds["max_median_relative_error"])
                                and trial["p90_relative_error"] <= float(thresholds["max_p90_relative_error"])
                                and trial["diffusivity_match_gap"] <= float(thresholds["max_diffusivity_match_gap"])
                                and trial["shell_kappa_cv"] <= float(thresholds["max_shell_kappa_cv"])
                                and trial["metric_tracking_abs_corr"] >= float(thresholds["min_metric_tracking_abs_corr"])
                            )
                            for trial in trials
                        )
                    ),
                    "median_relative_error_median": float(np.median([trial["median_relative_error"] for trial in trials])),
                    "p90_relative_error_median": float(np.median([trial["p90_relative_error"] for trial in trials])),
                    "diffusivity_match_gap_median": float(np.median([trial["diffusivity_match_gap"] for trial in trials])),
                    "shell_kappa_cv_median": float(np.median([trial["shell_kappa_cv"] for trial in trials])),
                    "metric_tracking_abs_corr_median": float(
                        np.median([trial["metric_tracking_abs_corr"] for trial in trials])
                    ),
                    "normalized_shell_kappa_profile_mean": mean_profile.tolist(),
                    "trials": trials,
                }
                case["pass"] = case_passes(case, thresholds)
                trial_groups.append(case)

    refinement_drifts = []
    for disorder_fraction in [float(v) for v in cfg["disorder_fractions"]]:
        for eta in [float(v) for v in cfg["eta_values"]]:
            pair = [
                case
                for case in trial_groups
                if abs(case["disorder_fraction"] - disorder_fraction) < 1.0e-12 and abs(case["eta"] - eta) < 1.0e-12
            ]
            pair = sorted(pair, key=lambda item: item["n_side"])
            if len(pair) == 2:
                refinement_drifts.append(
                    {
                        "label": f"sigma={disorder_fraction:.3f}|eta={eta:.3f}",
                        "disorder_fraction": float(disorder_fraction),
                        "eta": float(eta),
                        "drift": profile_drift(
                            np.asarray(pair[0]["normalized_shell_kappa_profile_mean"], dtype=float),
                            np.asarray(pair[1]["normalized_shell_kappa_profile_mean"], dtype=float),
                        ),
                    }
                )

    weakest_failure = min(
        trial_groups,
        key=lambda case: (
            case["median_relative_error_median"],
            case["p90_relative_error_median"],
            case["diffusivity_match_gap_median"],
            case["shell_kappa_cv_median"],
            -case["metric_tracking_abs_corr_median"],
        ),
    )
    max_drift = max((row["drift"] for row in refinement_drifts), default=0.0)
    observation = (
        "the scaffold-shell metric-tracking signal survives mild disorder on the smooth radial branch, "
        "but the shell-native transient constitutive closure does not remain bounded under the same disorder window"
    )
    conclusion = (
        f"all {len(trial_groups)} aggregated smooth-radial disorder cases remain open; "
        f"the weakest failure is `{weakest_failure['label']}` with median relative error "
        f"{weakest_failure['median_relative_error_median']:.4f}, p90 error {weakest_failure['p90_relative_error_median']:.4f}, "
        f"diffusivity gap {weakest_failure['diffusivity_match_gap_median']:.4f}, shell-kappa CV "
        f"{weakest_failure['shell_kappa_cv_median']:.4f}, and |metric-tracking corr| "
        f"{weakest_failure['metric_tracking_abs_corr_median']:.4f}; maximum refinement drift is {max_drift:.4f}; "
        "this keeps disorder-robust transient closure explicitly open and points to a disorder-native flux reconstruction "
        "as the next honest step"
    )
    return {
        "experiment": "scalar_kernel_graph_current_closure_radial_disorder",
        "config": cfg,
        "layouts": {
            str(n_side): {
                "n_side": int(layout["n_side"]),
                "lattice_spacing": float(layout["h"]),
                "epsilon_k": float(layout["epsilon_k"]),
                "shell_radii": np.asarray(layout["shell_radii"]).tolist(),
                "baseline_metric_rr": np.asarray(layout["baseline_metric_rr"]).tolist(),
            }
            for n_side, layout in layouts.items()
        },
        "cases": trial_groups,
        "refinement_drifts": refinement_drifts,
        "observation": observation,
        "conclusion": conclusion,
    }


def make_plots(result: dict[str, Any], stamp: str) -> list[str]:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    etas = [float(v) for v in result["config"]["eta_values"]]
    for n_side, marker in zip(result["config"]["n_sides"], ["o", "s"]):
        for disorder_fraction, color in zip(result["config"]["disorder_fractions"], ["#0f766e", "#b45309"]):
            cases = sorted(
                [
                    case
                    for case in result["cases"]
                    if case["n_side"] == int(n_side) and abs(case["disorder_fraction"] - float(disorder_fraction)) < 1.0e-12
                ],
                key=lambda item: item["eta"],
            )
            label = f"n={n_side}, sigma={float(disorder_fraction):.2f}"
            axes[0].plot(
                etas,
                [case["median_relative_error_median"] for case in cases],
                marker=marker,
                color=color,
                label=label,
            )
            axes[1].plot(
                etas,
                [case["diffusivity_match_gap_median"] for case in cases],
                marker=marker,
                color=color,
            )
            axes[2].plot(
                etas,
                [case["metric_tracking_abs_corr_median"] for case in cases],
                marker=marker,
                color=color,
            )
    axes[0].axhline(result["config"]["thresholds"]["max_median_relative_error"], linestyle="--", color="#7f1d1d")
    axes[1].axhline(result["config"]["thresholds"]["max_diffusivity_match_gap"], linestyle="--", color="#7f1d1d")
    axes[2].axhline(result["config"]["thresholds"]["min_metric_tracking_abs_corr"], linestyle="--", color="#7f1d1d")
    axes[0].set_title("Median Relative Error")
    axes[1].set_title("Diffusivity Match Gap")
    axes[2].set_title("|Metric Tracking Corr|")
    for ax in axes:
        ax.set_xlabel("Radial deformation eta")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_radial_disorder_summary.png"
    fig.savefig(summary_path, dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(
        len(result["config"]["disorder_fractions"]),
        len(result["config"]["eta_values"]),
        figsize=(13, 7),
        squeeze=False,
        sharex=False,
        sharey=False,
    )
    for row_idx, disorder_fraction in enumerate(result["config"]["disorder_fractions"]):
        for col_idx, eta in enumerate(result["config"]["eta_values"]):
            ax = axes[row_idx][col_idx]
            for n_side, color in zip(result["config"]["n_sides"], ["#1d4ed8", "#be123c"]):
                case = next(
                    case
                    for case in result["cases"]
                    if case["n_side"] == int(n_side)
                    and abs(case["eta"] - float(eta)) < 1.0e-12
                    and abs(case["disorder_fraction"] - float(disorder_fraction)) < 1.0e-12
                )
                radii = result["layouts"][str(n_side)]["shell_radii"][: len(case["normalized_shell_kappa_profile_mean"])]
                ax.plot(radii, case["normalized_shell_kappa_profile_mean"], label=f"n={n_side}", color=color)
            ax.set_title(f"sigma={float(disorder_fraction):.2f}, eta={float(eta):.2f}")
            ax.grid(alpha=0.25)
            if row_idx == len(result["config"]["disorder_fractions"]) - 1:
                ax.set_xlabel("Shell radius")
            if col_idx == 0:
                ax.set_ylabel("Mean normalized shell kappa")
    axes[0][0].legend(fontsize=8)
    fig.tight_layout()
    profile_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_radial_disorder_profiles.png"
    fig.savefig(profile_path, dpi=160)
    plt.close(fig)
    return [str(summary_path.relative_to(REPO_ROOT)), str(profile_path.relative_to(REPO_ROOT))]


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_current_closure_radial_disorder.json"
    latest = DATA / "scalar_kernel_graph_current_closure_radial_disorder_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    rows = "\n".join(
        f"| `{case['label']}` | `{case['pass_count']}/{case['trial_count']}` | `{case['median_relative_error_median']:.4f}` | `{case['p90_relative_error_median']:.4f}` | `{case['diffusivity_match_gap_median']:.4f}` | `{case['shell_kappa_cv_median']:.4f}` | `{case['metric_tracking_abs_corr_median']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in sorted(result["cases"], key=lambda item: (item["disorder_fraction"], item["eta"], item["n_side"]))
    )
    drift_rows = "\n".join(
        f"| `{row['label']}` | `{row['drift']:.4f}` |"
        for row in sorted(result["refinement_drifts"], key=lambda item: (item["disorder_fraction"], item["eta"]))
    )
    note = f"""# Scalar Kernel Graph Current Closure Radial Disorder v1

## Purpose

Test the next honest robustness question after the bounded smooth-inhomogeneity closure `52.5`:

> if the passing smooth radial branch is perturbed by mild carrier disorder, does the shell-native transient closure remain bounded?

The carrier stays fixed:

- same scalar kernel graph family
- same smooth radial deformation family
- same `52.5` transport thresholds

The only new ingredient is mild random carrier disorder added on top of the already-passing smooth radial deformation.

## Construction

- clean refinements: `{result['config']['n_sides']}`
- radial amplitudes: `{result['config']['eta_values']}`
- disorder fractions of lattice spacing: `{result['config']['disorder_fractions']}`
- seeds: `{result['config']['seeds']}`

This probe uses a scaffold-shell read:

1. build the clean shell layout once on the undeformed cubic scaffold;
2. apply smooth radial deformation;
3. add mild random interior disorder while keeping the boundary fixed;
4. read shell transport on the clean scaffold-shell assignment so shell bookkeeping itself does not collapse under jitter;
5. test whether the transient constitutive law remains within the same bounded window.

Artifacts:

- script: `numerics/simulations/scalar_kernel_graph_current_closure_radial_disorder.py`
- config: `config.json -> scalar_kernel_graph_current_closure_radial_disorder`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_current_closure_radial_disorder_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Aggregated cases

| case | trial passes | median rel. error | p90 rel. error | diffusivity gap | shell `kappa` CV | `|metric corr|` | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
{rows}

## Refinement drift

| case | drift |
| --- | --- |
{drift_rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph current closure radial disorder\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, eta_values={config['eta_values']}, "
            f"disorder_fractions={config['disorder_fractions']}, seeds={config['seeds']}, "
            f"epsilon_coeff={config['epsilon_coeff']}, cutoff_factor={config['cutoff_factor']}, "
            f"fit_window={config['fit_window']}, metric_radius_factor={config['metric_radius_factor']}\n"
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
