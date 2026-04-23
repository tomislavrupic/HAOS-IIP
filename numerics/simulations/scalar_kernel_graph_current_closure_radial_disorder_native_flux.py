#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from scalar_kernel_graph_current_closure_inhomogeneity import normalized_profile, profile_drift
from scalar_kernel_graph_current_closure_radial_disorder import apply_disorder, radial_deformation_profile
from scalar_kernel_graph_geometry_bridge import build_cubic_points, heat_metrics, heat_trajectory
from scalar_kernel_graph_geometry_robustness import build_point_cloud_operator
from scalar_kernel_graph_metric_field import bulk_mask, coarse_grain_metric_field, local_metric_field

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Current_Closure_Radial_Disorder_Native_Flux_v1.md"
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
    "fit_window": [0.006, 0.02],
    "bulk_margin_layers": 2,
    "metric_radius_factor": 2.5,
    "shell_r_max": 0.45,
    "gradient_mode": "local_linear",
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
        if isinstance(raw.get("scalar_kernel_graph_current_closure_radial_disorder_native_flux"), dict):
            merged.update(
                {
                    k: v
                    for k, v in raw["scalar_kernel_graph_current_closure_radial_disorder_native_flux"].items()
                    if k != "thresholds"
                }
            )
            if isinstance(raw["scalar_kernel_graph_current_closure_radial_disorder_native_flux"].get("thresholds"), dict):
                merged["thresholds"].update(
                    raw["scalar_kernel_graph_current_closure_radial_disorder_native_flux"]["thresholds"]
                )
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def local_linear_gradient(values: np.ndarray, radii: np.ndarray) -> np.ndarray:
    out = np.zeros_like(values)
    n_shells = values.shape[1]
    for shell_idx in range(n_shells):
        lo = max(0, shell_idx - 1)
        hi = min(n_shells, shell_idx + 2)
        x = radii[lo:hi]
        design = np.column_stack([x, np.ones_like(x)])
        for time_idx in range(values.shape[0]):
            y = values[time_idx, lo:hi]
            slope, _ = np.linalg.lstsq(design, y, rcond=None)[0]
            out[time_idx, shell_idx] = slope
    return out


def build_native_shells(
    bulk_points: np.ndarray,
    bulk_metric: np.ndarray,
    bulk_profile: np.ndarray,
    source_point: np.ndarray,
    h: float,
    shell_r_max: float,
) -> tuple[np.ndarray, list[np.ndarray], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    bulk_radii = np.linalg.norm(bulk_points - source_point, axis=1)
    labels = np.rint(bulk_radii / max(h, 1.0e-12)).astype(int)
    shell_radii: list[float] = []
    shell_masks: list[np.ndarray] = []
    shell_metric_rr: list[float] = []
    shell_counts: list[int] = []
    shell_profile: list[float] = []
    for label in np.unique(labels):
        mask_shell = labels == label
        if not np.any(mask_shell):
            continue
        mean_r = float(np.mean(bulk_radii[mask_shell]))
        if mean_r <= 2.0 * h or mean_r > float(shell_r_max):
            continue
        radial_dirs = bulk_points[mask_shell] - source_point
        radial_dirs = radial_dirs / np.maximum(np.linalg.norm(radial_dirs, axis=1, keepdims=True), 1.0e-18)
        shell_radii.append(mean_r)
        shell_masks.append(mask_shell)
        shell_metric_rr.append(
            float(np.mean(np.einsum("ni,nij,nj->n", radial_dirs, bulk_metric[mask_shell], radial_dirs)))
        )
        shell_counts.append(int(np.sum(mask_shell)))
        shell_profile.append(float(np.mean(bulk_profile[mask_shell])))
    return (
        np.asarray(shell_radii, dtype=float),
        shell_masks,
        np.asarray(shell_metric_rr, dtype=float),
        np.asarray(shell_counts, dtype=float),
        np.asarray(shell_profile, dtype=float),
        bulk_radii,
    )


def baseline_case(
    n_side: int,
    epsilon_coeff: float,
    cutoff_factor: float,
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
    raw_metric = local_metric_field(points=points, delta=delta)
    mask = bulk_mask(points=points, h=h, bulk_margin_layers=bulk_margin_layers)
    bulk_points = points[mask]
    bulk_metric = coarse_grain_metric_field(points=bulk_points, field=raw_metric[mask], radius=float(metric_radius_factor) * h)
    shell_radii, _, shell_metric_rr, shell_counts, _, _ = build_native_shells(
        bulk_points=bulk_points,
        bulk_metric=bulk_metric,
        bulk_profile=np.zeros(len(bulk_points), dtype=float),
        source_point=points[source_idx],
        h=h,
        shell_r_max=shell_r_max,
    )
    return {
        "n_side": int(n_side),
        "lattice_spacing": float(h),
        "epsilon_k": float(epsilon_k),
        "shell_radii": shell_radii.tolist(),
        "shell_metric_rr": shell_metric_rr.tolist(),
        "shell_counts": shell_counts.tolist(),
    }


def evaluate_trial(
    n_side: int,
    eta: float,
    disorder_fraction: float,
    seed: int,
    epsilon_coeff: float,
    cutoff_factor: float,
    heat_time_max: float,
    heat_time_steps: int,
    fit_window: list[float],
    bulk_margin_layers: int,
    metric_radius_factor: float,
    shell_r_max: float,
    gradient_mode: str,
    baseline_radii: np.ndarray,
    baseline_metric_rr: np.ndarray,
) -> dict[str, Any]:
    points_clean, _, source_idx, h = build_cubic_points(n_side)
    deformed_points, profile = radial_deformation_profile(points_clean, n_side=n_side, eta=eta)
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
    bulk_profile = profile[mask]
    bulk_metric = coarse_grain_metric_field(
        points=bulk_points,
        field=raw_metric[mask],
        radius=float(metric_radius_factor) * h,
    )
    shell_radii, shell_masks, shell_metric_rr, shell_counts, shell_profile, bulk_radii = build_native_shells(
        bulk_points=bulk_points,
        bulk_metric=bulk_metric,
        bulk_profile=bulk_profile,
        source_point=points[source_idx],
        h=h,
        shell_r_max=shell_r_max,
    )
    shell_mass = np.stack([np.sum(bulk_trajectory[:, mask_shell], axis=1) for mask_shell in shell_masks], axis=1)
    shell_density = shell_mass / (shell_counts[None, :] * (h**3))
    if gradient_mode == "local_linear":
        density_gradient = local_linear_gradient(shell_density, shell_radii)
    elif gradient_mode == "centered":
        density_gradient = np.gradient(shell_density, shell_radii, axis=1)
    else:
        raise ValueError(f"unknown gradient mode: {gradient_mode}")

    cumulative_mass = np.stack([np.sum(bulk_trajectory[:, bulk_radii <= shell_r], axis=1) for shell_r in shell_radii], axis=1)
    empirical_current = -np.gradient(cumulative_mass, times, axis=0)
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
    normalized_shell_kappa_profile = normalized_profile(shell_kappa_profile)
    baseline_interp = np.interp(shell_radii, baseline_radii, baseline_metric_rr)
    metric_shift = shell_metric_rr - baseline_interp
    metric_tracking_corr = float(np.corrcoef(shell_profile, metric_shift)[0, 1])

    heat_diffusivity = float(heat["effective_diffusivity"])
    conductivity_over_diffusivity = float(conductivity_fit / max(heat_diffusivity, 1.0e-18))
    return {
        "label": f"radial_native_flux|n={n_side}|eta={eta:.3f}|sigma={disorder_fraction:.3f}|seed={seed}",
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
        "metric_tracking_corr": metric_tracking_corr,
        "metric_tracking_abs_corr": float(abs(metric_tracking_corr)),
        "shell_radii": shell_radii.tolist(),
        "shell_profile": shell_profile.tolist(),
        "shell_metric_rr": shell_metric_rr.tolist(),
        "metric_shift_profile": metric_shift.tolist(),
        "shell_kappa_profile": shell_kappa_profile.tolist(),
        "normalized_shell_kappa_profile": normalized_shell_kappa_profile.tolist(),
    }


def trial_passes(trial: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        trial["median_relative_error"] <= float(thresholds["max_median_relative_error"])
        and trial["p90_relative_error"] <= float(thresholds["max_p90_relative_error"])
        and trial["diffusivity_match_gap"] <= float(thresholds["max_diffusivity_match_gap"])
        and trial["shell_kappa_cv"] <= float(thresholds["max_shell_kappa_cv"])
        and trial["metric_tracking_abs_corr"] >= float(thresholds["min_metric_tracking_abs_corr"])
    )


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
    baseline_cases = {
        int(n_side): baseline_case(
            n_side=int(n_side),
            epsilon_coeff=float(cfg["epsilon_coeff"]),
            cutoff_factor=float(cfg["cutoff_factor"]),
            bulk_margin_layers=int(cfg["bulk_margin_layers"]),
            metric_radius_factor=float(cfg["metric_radius_factor"]),
            shell_r_max=float(cfg["shell_r_max"]),
        )
        for n_side in [int(v) for v in cfg["n_sides"]]
    }

    cases = []
    all_trials = []
    for disorder_fraction in [float(v) for v in cfg["disorder_fractions"]]:
        for eta in [float(v) for v in cfg["eta_values"]]:
            for n_side in [int(v) for v in cfg["n_sides"]]:
                trials = [
                    evaluate_trial(
                        n_side=n_side,
                        eta=eta,
                        disorder_fraction=disorder_fraction,
                        seed=int(seed),
                        epsilon_coeff=float(cfg["epsilon_coeff"]),
                        cutoff_factor=float(cfg["cutoff_factor"]),
                        heat_time_max=float(cfg["heat_time_max"]),
                        heat_time_steps=int(cfg["heat_time_steps"]),
                        fit_window=[float(v) for v in cfg["fit_window"]],
                        bulk_margin_layers=int(cfg["bulk_margin_layers"]),
                        metric_radius_factor=float(cfg["metric_radius_factor"]),
                        shell_r_max=float(cfg["shell_r_max"]),
                        gradient_mode=str(cfg["gradient_mode"]),
                        baseline_radii=np.asarray(baseline_cases[n_side]["shell_radii"], dtype=float),
                        baseline_metric_rr=np.asarray(baseline_cases[n_side]["shell_metric_rr"], dtype=float),
                    )
                    for seed in [int(v) for v in cfg["seeds"]]
                ]
                for trial in trials:
                    trial["pass"] = trial_passes(trial, thresholds)
                all_trials.extend(trials)
                min_len = min(len(trial["normalized_shell_kappa_profile"]) for trial in trials)
                mean_profile = np.mean(
                    [np.asarray(trial["normalized_shell_kappa_profile"][:min_len], dtype=float) for trial in trials],
                    axis=0,
                )
                case = {
                    "label": f"radial_native_flux|n={n_side}|eta={eta:.3f}|sigma={disorder_fraction:.3f}",
                    "n_side": int(n_side),
                    "eta": float(eta),
                    "disorder_fraction": float(disorder_fraction),
                    "trial_count": len(trials),
                    "pass_count": int(sum(int(trial["pass"]) for trial in trials)),
                    "median_relative_error_median": float(np.median([trial["median_relative_error"] for trial in trials])),
                    "p90_relative_error_median": float(np.median([trial["p90_relative_error"] for trial in trials])),
                    "diffusivity_match_gap_median": float(np.median([trial["diffusivity_match_gap"] for trial in trials])),
                    "shell_kappa_cv_median": float(np.median([trial["shell_kappa_cv"] for trial in trials])),
                    "metric_tracking_abs_corr_median": float(
                        np.median([trial["metric_tracking_abs_corr"] for trial in trials])
                    ),
                    "conductivity_over_diffusivity_median": float(
                        np.median([trial["conductivity_over_diffusivity"] for trial in trials])
                    ),
                    "normalized_shell_kappa_profile_mean": mean_profile.tolist(),
                    "trials": trials,
                }
                case["pass"] = case_passes(case, thresholds)
                cases.append(case)

    refinement_drifts = []
    for disorder_fraction in [float(v) for v in cfg["disorder_fractions"]]:
        for eta in [float(v) for v in cfg["eta_values"]]:
            pair = [
                case
                for case in cases
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

    max_drift = max((row["drift"] for row in refinement_drifts), default=0.0)
    all_cases_pass = all(case["pass"] for case in cases)
    drift_pass = max_drift <= float(thresholds["max_refinement_profile_drift"])
    weakest_case = max(
        cases,
        key=lambda item: (
            item["median_relative_error_median"],
            item["p90_relative_error_median"],
            item["diffusivity_match_gap_median"],
            item["shell_kappa_cv_median"],
        ),
    )
    weakest_trial = max(
        all_trials,
        key=lambda item: (
            item["median_relative_error"],
            item["p90_relative_error"],
            item["diffusivity_match_gap"],
            item["shell_kappa_cv"],
        ),
    )
    min_pass_count = min((case["pass_count"] for case in cases), default=0)

    if all_cases_pass and drift_pass:
        observation = (
            "once the smooth-radial mild-disorder family is read with native bulk shells, interior-ball cumulative flux, "
            "and a delayed asymptotic fit window, the bounded constitutive closure returns on aggregated observables"
        )
        conclusion = (
            f"all {len(cases)} aggregated mild-disorder cases pass, with maximum refinement drift {max_drift:.4f}; "
            f"the weakest aggregated case is `{weakest_case['label']}` with kappa/D_eff "
            f"{weakest_case['conductivity_over_diffusivity_median']:.4f}, median relative error "
            f"{weakest_case['median_relative_error_median']:.4f}, p90 error {weakest_case['p90_relative_error_median']:.4f}, "
            f"shell-kappa CV {weakest_case['shell_kappa_cv_median']:.4f}, and |metric-tracking corr| "
            f"{weakest_case['metric_tracking_abs_corr_median']:.4f}; the hardest seed-level trial is "
            f"`{weakest_trial['label']}` with median relative error {weakest_trial['median_relative_error']:.4f}, "
            f"p90 error {weakest_trial['p90_relative_error']:.4f}, shell-kappa CV {weakest_trial['shell_kappa_cv']:.4f}, "
            f"and pass state {str(weakest_trial['pass']).upper()}; this supports a bounded disorder-native transient "
            f"closure on the smooth-radial branch, with the minimum seed-level pass count still only {min_pass_count}/3 "
            "in the hardest tested case"
        )
    else:
        observation = (
            "the disorder-native flux reconstruction restores part of the smooth-radial mild-disorder transport law, "
            "but the aggregated closure still remains open on the tested window"
        )
        conclusion = (
            f"aggregated passing cases = {sum(int(case['pass']) for case in cases)}/{len(cases)}, maximum refinement drift "
            f"{max_drift:.4f}; the weakest aggregated case is `{weakest_case['label']}` with median relative error "
            f"{weakest_case['median_relative_error_median']:.4f}, p90 error {weakest_case['p90_relative_error_median']:.4f}, "
            f"diffusivity gap {weakest_case['diffusivity_match_gap_median']:.4f}, and shell-kappa CV "
            f"{weakest_case['shell_kappa_cv_median']:.4f}; this keeps disorder-native transient closure explicitly open"
        )

    return {
        "experiment": "scalar_kernel_graph_current_closure_radial_disorder_native_flux",
        "config": cfg,
        "baseline_cases": baseline_cases,
        "cases": cases,
        "refinement_drifts": refinement_drifts,
        "max_refinement_profile_drift": float(max_drift),
        "observation": observation,
        "conclusion": conclusion,
    }


def make_plots(result: dict[str, Any], stamp: str) -> list[str]:
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
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
                [case["p90_relative_error_median"] for case in cases],
                marker=marker,
                color=color,
            )
            axes[2].plot(
                etas,
                [case["diffusivity_match_gap_median"] for case in cases],
                marker=marker,
                color=color,
            )
            axes[3].plot(
                etas,
                [case["metric_tracking_abs_corr_median"] for case in cases],
                marker=marker,
                color=color,
            )
    axes[0].axhline(result["config"]["thresholds"]["max_median_relative_error"], linestyle="--", color="#7f1d1d")
    axes[1].axhline(result["config"]["thresholds"]["max_p90_relative_error"], linestyle="--", color="#7f1d1d")
    axes[2].axhline(result["config"]["thresholds"]["max_diffusivity_match_gap"], linestyle="--", color="#7f1d1d")
    axes[3].axhline(result["config"]["thresholds"]["min_metric_tracking_abs_corr"], linestyle="--", color="#7f1d1d")
    axes[0].set_title("Median Relative Error")
    axes[1].set_title("P90 Relative Error")
    axes[2].set_title("Diffusivity Match Gap")
    axes[3].set_title("|Metric Tracking Corr|")
    for ax in axes:
        ax.set_xlabel("Radial deformation eta")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_radial_disorder_native_flux_summary.png"
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
                radii = result["baseline_cases"][int(n_side)]["shell_radii"][: len(case["normalized_shell_kappa_profile_mean"])]
                ax.plot(radii, case["normalized_shell_kappa_profile_mean"], label=f"n={n_side}", color=color)
            ax.set_title(f"sigma={float(disorder_fraction):.2f}, eta={float(eta):.2f}")
            ax.grid(alpha=0.25)
            if row_idx == len(result["config"]["disorder_fractions"]) - 1:
                ax.set_xlabel("Shell radius")
            if col_idx == 0:
                ax.set_ylabel("Mean normalized shell kappa")
    axes[0][0].legend(fontsize=8)
    fig.tight_layout()
    profile_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_radial_disorder_native_flux_profiles.png"
    fig.savefig(profile_path, dpi=160)
    plt.close(fig)
    return [str(summary_path.relative_to(REPO_ROOT)), str(profile_path.relative_to(REPO_ROOT))]


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_current_closure_radial_disorder_native_flux.json"
    latest = DATA / "scalar_kernel_graph_current_closure_radial_disorder_native_flux_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    rows = "\n".join(
        f"| `{case['label']}` | `{case['pass_count']}/{case['trial_count']}` | `{case['conductivity_over_diffusivity_median']:.4f}` | `{case['median_relative_error_median']:.4f}` | `{case['p90_relative_error_median']:.4f}` | `{case['shell_kappa_cv_median']:.4f}` | `{case['metric_tracking_abs_corr_median']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in sorted(result["cases"], key=lambda item: (item["disorder_fraction"], item["eta"], item["n_side"]))
    )
    drift_rows = "\n".join(
        f"| `{row['label']}` | `{row['drift']:.4f}` |"
        for row in sorted(result["refinement_drifts"], key=lambda item: (item["disorder_fraction"], item["eta"]))
    )
    note = f"""# Scalar Kernel Graph Current Closure Radial Disorder Native Flux v1

## Purpose

Test the missing transport bridge after the scaffold-shell disorder failure:

> if smooth radial mild-disorder cases are read with a disorder-native interior-ball flux instead of scaffold-shell bookkeeping, does the bounded constitutive law return?

The carrier stays fixed:

- same scalar kernel graph family
- same smooth radial deformation family
- same mild-disorder window
- same constitutive target

What changes is only the transport readout:

- native bulk shells built from the actual disordered radii
- interior-ball cumulative mass on the actual bulk nodes
- local-linear radial density gradient
- delayed asymptotic fit window `{result['config']['fit_window']}` to exclude the earliest source-core transient layer

## Construction

- clean refinements: `{result['config']['n_sides']}`
- radial amplitudes: `{result['config']['eta_values']}`
- disorder fractions of lattice spacing: `{result['config']['disorder_fractions']}`
- seeds: `{result['config']['seeds']}`
- gradient mode: `{result['config']['gradient_mode']}`

Artifacts:

- script: `numerics/simulations/scalar_kernel_graph_current_closure_radial_disorder_native_flux.py`
- config: `config.json -> scalar_kernel_graph_current_closure_radial_disorder_native_flux`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_current_closure_radial_disorder_native_flux_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Aggregated cases

| case | trial passes | `kappa / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | `|metric corr|` | status |
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
        handle.write("\n## scalar kernel graph current closure radial disorder native flux\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, eta_values={config['eta_values']}, "
            f"disorder_fractions={config['disorder_fractions']}, seeds={config['seeds']}, "
            f"fit_window={config['fit_window']}, gradient_mode={config['gradient_mode']}, "
            f"epsilon_coeff={config['epsilon_coeff']}, cutoff_factor={config['cutoff_factor']}, "
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
