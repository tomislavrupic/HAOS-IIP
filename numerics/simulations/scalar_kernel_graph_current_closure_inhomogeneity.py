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
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Current_Closure_Inhomogeneity_v1.md"
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
    "families": ["radial", "bump"],
    "eta_values": [0.05, 0.10, 0.15],
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

DEFORMATION_LABELS = {
    "radial": "smooth radial deformation",
    "bump": "localized scalar bump",
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_current_closure_inhomogeneity"), dict):
            merged.update(
                {k: v for k, v in raw["scalar_kernel_graph_current_closure_inhomogeneity"].items() if k != "thresholds"}
            )
            if isinstance(raw["scalar_kernel_graph_current_closure_inhomogeneity"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_current_closure_inhomogeneity"]["thresholds"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def deformation_profile(points: np.ndarray, n_side: int, family: str, eta: float) -> tuple[np.ndarray, np.ndarray]:
    center = np.array([0.5, 0.5, 0.5], dtype=float)
    delta = points - center
    radii = np.linalg.norm(delta, axis=1)
    grid_ids = np.indices((n_side, n_side, n_side)).reshape(3, -1).T
    boundary_mask = np.any((grid_ids == 0) | (grid_ids == n_side - 1), axis=1)

    if family == "radial":
        r_max = max(float(np.max(radii)), 1.0e-18)
        profile = (radii / r_max) ** 2
    elif family == "bump":
        sigma = 0.18
        profile = np.exp(-(radii**2) / (2.0 * sigma * sigma))
    else:
        raise ValueError(f"unknown deformation family: {family}")

    factor = 1.0 + float(eta) * profile
    out = center + factor[:, None] * delta
    out[boundary_mask] = points[boundary_mask]
    return np.clip(out, 0.0, 1.0), profile


def build_bulk_shells(
    bulk_points: np.ndarray,
    bulk_metric: np.ndarray,
    bulk_profile: np.ndarray,
    source_point: np.ndarray,
    h: float,
    shell_r_max: float,
) -> tuple[np.ndarray, list[np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    bulk_radii = np.linalg.norm(bulk_points - source_point, axis=1)
    labels = np.round(bulk_radii / max(h, 1.0e-12), 8) * h
    shell_radii: list[float] = []
    shell_masks: list[np.ndarray] = []
    shell_metric_rr: list[float] = []
    shell_profile: list[float] = []
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
        shell_profile.append(float(np.mean(bulk_profile[mask_shell])))
        shell_counts.append(int(np.sum(mask_shell)))
    return (
        np.asarray(shell_radii, dtype=float),
        shell_masks,
        np.asarray(shell_metric_rr, dtype=float),
        np.asarray(shell_profile, dtype=float),
        np.asarray(shell_counts, dtype=float),
    )


def profile_drift(profile_a: np.ndarray, profile_b: np.ndarray) -> float:
    min_len = min(len(profile_a), len(profile_b))
    return float(np.linalg.norm(profile_a[:min_len] - profile_b[:min_len]) / np.sqrt(min_len))


def normalized_profile(profile: np.ndarray) -> np.ndarray:
    scale = max(float(np.mean(np.abs(profile))), 1.0e-18)
    return profile / scale


def evaluate_baseline(
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
    shell_radii, _, shell_metric_rr, _, _ = build_bulk_shells(
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
    }


def evaluate_case(
    n_side: int,
    family: str,
    eta: float,
    epsilon_coeff: float,
    cutoff_factor: float,
    heat_time_max: float,
    heat_time_steps: int,
    fit_window: list[float],
    bulk_margin_layers: int,
    metric_radius_factor: float,
    shell_r_max: float,
    baseline_shell_metric_rr: np.ndarray,
) -> dict[str, Any]:
    base_points, _, source_idx, h = build_cubic_points(n_side)
    points, profile = deformation_profile(base_points, n_side=n_side, family=family, eta=eta)
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
    bulk_metric = coarse_grain_metric_field(points=bulk_points, field=raw_metric[mask], radius=float(metric_radius_factor) * h)
    shell_radii, shell_masks, shell_metric_rr, shell_profile, shell_counts = build_bulk_shells(
        bulk_points=bulk_points,
        bulk_metric=bulk_metric,
        bulk_profile=bulk_profile,
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
    normalized_shell_kappa_profile = normalized_profile(shell_kappa_profile)

    min_len = min(len(shell_metric_rr), len(baseline_shell_metric_rr))
    metric_shift = shell_metric_rr[:min_len] - baseline_shell_metric_rr[:min_len]
    metric_shift_corr = float(np.corrcoef(shell_profile[:min_len], metric_shift)[0, 1])
    normalized_metric_shift = normalized_profile(metric_shift)

    heat_diffusivity = float(heat["effective_diffusivity"])
    conductivity_over_diffusivity = float(conductivity_fit / max(heat_diffusivity, 1.0e-18))
    return {
        "label": f"{family}|n={n_side}|eta={eta:.3f}",
        "family": family,
        "family_label": DEFORMATION_LABELS[family],
        "eta": float(eta),
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
        "shell_kappa_cv": float(np.nanstd(shell_kappa_profile) / max(float(np.nanmean(np.abs(shell_kappa_profile))), 1.0e-18)),
        "metric_tracking_corr": metric_shift_corr,
        "metric_tracking_abs_corr": float(abs(metric_shift_corr)),
        "metric_shift_span": float(np.max(metric_shift) - np.min(metric_shift)),
        "shell_radii": shell_radii.tolist(),
        "shell_profile": shell_profile.tolist(),
        "shell_metric_rr": shell_metric_rr.tolist(),
        "metric_shift_profile": metric_shift.tolist(),
        "normalized_metric_shift_profile": normalized_metric_shift.tolist(),
        "shell_kappa_profile": shell_kappa_profile.tolist(),
        "normalized_shell_kappa_profile": normalized_shell_kappa_profile.tolist(),
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        case["median_relative_error"] <= float(thresholds["max_median_relative_error"])
        and case["p90_relative_error"] <= float(thresholds["max_p90_relative_error"])
        and case["diffusivity_match_gap"] <= float(thresholds["max_diffusivity_match_gap"])
        and case["shell_kappa_cv"] <= float(thresholds["max_shell_kappa_cv"])
        and case["metric_tracking_abs_corr"] >= float(thresholds["min_metric_tracking_abs_corr"])
    )


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    thresholds = cfg["thresholds"]

    baseline_cases = {
        int(n_side): evaluate_baseline(
            n_side=int(n_side),
            epsilon_coeff=float(cfg["epsilon_coeff"]),
            cutoff_factor=float(cfg["cutoff_factor"]),
            bulk_margin_layers=int(cfg["bulk_margin_layers"]),
            metric_radius_factor=float(cfg["metric_radius_factor"]),
            shell_r_max=float(cfg["shell_r_max"]),
        )
        for n_side in [int(v) for v in cfg["n_sides"]]
    }

    cases: list[dict[str, Any]] = []
    for family in [str(v) for v in cfg["families"]]:
        for eta in [float(v) for v in cfg["eta_values"]]:
            for n_side in [int(v) for v in cfg["n_sides"]]:
                case = evaluate_case(
                    n_side=n_side,
                    family=family,
                    eta=eta,
                    epsilon_coeff=float(cfg["epsilon_coeff"]),
                    cutoff_factor=float(cfg["cutoff_factor"]),
                    heat_time_max=float(cfg["heat_time_max"]),
                    heat_time_steps=int(cfg["heat_time_steps"]),
                    fit_window=[float(v) for v in cfg["fit_window"]],
                    bulk_margin_layers=int(cfg["bulk_margin_layers"]),
                    metric_radius_factor=float(cfg["metric_radius_factor"]),
                    shell_r_max=float(cfg["shell_r_max"]),
                    baseline_shell_metric_rr=np.asarray(baseline_cases[n_side]["shell_metric_rr"], dtype=float),
                )
                case["pass"] = case_passes(case, thresholds)
                cases.append(case)

    refinement_drifts = []
    for family in [str(v) for v in cfg["families"]]:
        for eta in [float(v) for v in cfg["eta_values"]]:
            pair = [case for case in cases if case["family"] == family and abs(case["eta"] - eta) < 1.0e-12]
            pair = sorted(pair, key=lambda item: item["n_side"])
            if len(pair) == 2:
                refinement_drifts.append(
                    {
                        "label": f"{family}|eta={eta:.3f}",
                        "family": family,
                        "eta": eta,
                        "drift": profile_drift(
                            np.asarray(pair[0]["normalized_shell_kappa_profile"], dtype=float),
                            np.asarray(pair[1]["normalized_shell_kappa_profile"], dtype=float),
                        ),
                    }
                )

    family_summaries = []
    max_drift = 0.0
    for family in [str(v) for v in cfg["families"]]:
        family_cases = [case for case in cases if case["family"] == family]
        family_drifts = [row for row in refinement_drifts if row["family"] == family]
        family_max_drift = max((row["drift"] for row in family_drifts), default=0.0)
        max_drift = max(max_drift, family_max_drift)
        all_cases_pass = all(case["pass"] for case in family_cases)
        drift_pass = family_max_drift <= float(thresholds["max_refinement_profile_drift"])
        family_summaries.append(
            {
                "family": family,
                "family_label": DEFORMATION_LABELS[family],
                "case_count": len(family_cases),
                "all_cases_pass": bool(all_cases_pass),
                "drift_pass": bool(drift_pass),
                "max_refinement_profile_drift": float(family_max_drift),
            }
        )

    radial_summary = next((row for row in family_summaries if row["family"] == "radial"), None)
    bump_summary = next((row for row in family_summaries if row["family"] == "bump"), None)
    weakest_passing = max(
        [case for case in cases if case["pass"]],
        key=lambda item: (
            item["median_relative_error"],
            item["p90_relative_error"],
            item["diffusivity_match_gap"],
            item["shell_kappa_cv"],
        ),
        default=None,
    )
    weakest_failing = max(
        [case for case in cases if not case["pass"]],
        key=lambda item: (
            item["median_relative_error"],
            item["p90_relative_error"],
            item["diffusivity_match_gap"],
            item["shell_kappa_cv"],
        ),
        default=None,
    )

    if radial_summary is not None and radial_summary["all_cases_pass"] and radial_summary["drift_pass"]:
        if bump_summary is not None and (not bump_summary["all_cases_pass"] or not bump_summary["drift_pass"]):
            observation = (
                "the extracted metric and the shell-native constitutive law track one another under smooth radial inhomogeneity on the same scalar carrier, "
                "while a more localized bump deformation remains too sharp for the same bounded transport closure"
            )
            conclusion = (
                f"all {radial_summary['case_count']} smooth radial cases pass, with maximum radial refinement drift {radial_summary['max_refinement_profile_drift']:.4f}; "
                f"the weakest radial passing case is `{weakest_passing['label']}` with kappa/D_eff {weakest_passing['conductivity_over_diffusivity']:.4f}, "
                f"median relative error {weakest_passing['median_relative_error']:.4f}, p90 error {weakest_passing['p90_relative_error']:.4f}, "
                f"and |metric-tracking corr| {weakest_passing['metric_tracking_abs_corr']:.4f}; "
                f"by contrast, the sharpest failing bump case is `{weakest_failing['label']}` with kappa/D_eff {weakest_failing['conductivity_over_diffusivity']:.4f}, "
                f"median relative error {weakest_failing['median_relative_error']:.4f}, p90 error {weakest_failing['p90_relative_error']:.4f}, "
                f"and shell-kappa CV {weakest_failing['shell_kappa_cv']:.4f}; "
                "this supports a bounded smooth-inhomogeneity closure on the scalar carrier, while keeping localized inhomogeneity explicitly open"
            )
        else:
            observation = (
                "the extracted metric and the shell-native constitutive law now track the designed deformation together across the tested inhomogeneity window"
            )
            conclusion = (
                f"all tested inhomogeneity cases pass, the maximum refinement drift is {max_drift:.4f}, "
                f"and the weakest passing case is `{weakest_passing['label']}` with kappa/D_eff {weakest_passing['conductivity_over_diffusivity']:.4f}, "
                f"median relative error {weakest_passing['median_relative_error']:.4f}, p90 error {weakest_passing['p90_relative_error']:.4f}, "
                f"and |metric-tracking corr| {weakest_passing['metric_tracking_abs_corr']:.4f}; "
                "this supports a bounded controlled-inhomogeneity closure on the scalar carrier"
            )
    else:
        observation = (
            "the current shell-native transport law does not yet stay stable enough when the scalar carrier is intentionally deformed, even though the extracted metric still responds coherently to the designed inhomogeneity"
        )
        conclusion = (
            f"open inhomogeneity cases remain, with maximum refinement drift {max_drift:.4f}; "
            f"the weakest failing case is `{weakest_failing['label']}` with kappa/D_eff {weakest_failing['conductivity_over_diffusivity']:.4f}, "
            f"median relative error {weakest_failing['median_relative_error']:.4f}, p90 error {weakest_failing['p90_relative_error']:.4f}, "
            f"and shell-kappa CV {weakest_failing['shell_kappa_cv']:.4f}; "
            "the correct boundary therefore stays at clean-line current closure rather than a positive inhomogeneity statement"
        )

    return {
        "experiment": "scalar_kernel_graph_current_closure_inhomogeneity",
        "config": cfg,
        "baseline_cases": baseline_cases,
        "cases": cases,
        "refinement_drifts": refinement_drifts,
        "family_summaries": family_summaries,
        "max_refinement_profile_drift": float(max_drift),
        "observation": observation,
        "conclusion": conclusion,
    }


def make_plots(result: dict[str, Any], stamp: str) -> list[str]:
    plot_paths: list[str] = []
    cases = result["cases"]

    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_inhomogeneity_summary.png"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for family in result["config"]["families"]:
        subset = [case for case in cases if case["family"] == family]
        for n_side in sorted({case["n_side"] for case in subset}):
            rows = sorted([case for case in subset if case["n_side"] == n_side], key=lambda item: item["eta"])
            eta = [case["eta"] for case in rows]
            label = f"{family}|n={n_side}"
            axes[0, 0].plot(eta, [case["conductivity_over_diffusivity"] for case in rows], marker="o", label=label)
            axes[0, 1].plot(eta, [case["p90_relative_error"] for case in rows], marker="o", label=label)
            axes[1, 0].plot(eta, [case["shell_kappa_cv"] for case in rows], marker="o", label=label)
            axes[1, 1].plot(eta, [case["metric_tracking_abs_corr"] for case in rows], marker="o", label=label)

    axes[0, 0].axhline(1.0, color="k", linestyle="--", linewidth=1.0)
    axes[0, 0].set_title("Conductivity / diffusivity")
    axes[0, 0].set_xlabel("deformation amplitude eta")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].set_title("p90 relative error")
    axes[0, 1].set_xlabel("deformation amplitude eta")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].set_title("Shell-kappa CV")
    axes[1, 0].set_xlabel("deformation amplitude eta")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].set_title("Metric-tracking |corr|")
    axes[1, 1].set_xlabel("deformation amplitude eta")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(fontsize=8)

    fig.savefig(summary_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(summary_path.relative_to(REPO_ROOT)))

    detail_path = PLOTS / f"{stamp}_scalar_kernel_graph_current_closure_inhomogeneity_profiles.png"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    selected = []
    for family in result["config"]["families"]:
        family_cases = [case for case in cases if case["family"] == family]
        if family_cases:
            selected.append(sorted(family_cases, key=lambda item: (item["n_side"], item["eta"]))[-1])

    for axis, case in zip(axes, selected):
        radii = np.asarray(case["shell_radii"], dtype=float)
        profile = np.asarray(case["shell_profile"], dtype=float)
        metric_shift = np.asarray(case["normalized_metric_shift_profile"], dtype=float)
        kappa_profile = np.asarray(case["normalized_shell_kappa_profile"], dtype=float)
        axis.plot(radii[: len(profile)], profile / max(float(np.mean(np.abs(profile))), 1.0e-18), marker="o", label="deformation profile")
        axis.plot(radii[: len(metric_shift)], metric_shift, marker="s", label="metric shift")
        axis.plot(radii[: len(kappa_profile)], kappa_profile, marker="^", label="shell kappa")
        axis.axhline(1.0, color="k", linestyle="--", linewidth=1.0, alpha=0.5)
        axis.set_title(case["label"])
        axis.set_xlabel("shell radius")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)

    fig.savefig(detail_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(detail_path.relative_to(REPO_ROOT)))

    return plot_paths


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_current_closure_inhomogeneity.json"
    latest = DATA / "scalar_kernel_graph_current_closure_inhomogeneity_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    sections = []
    for family in result["config"]["families"]:
        family_cases = sorted(
            [case for case in result["cases"] if case["family"] == family],
            key=lambda item: (item["n_side"], item["eta"]),
        )
        rows = "\n".join(
            f"| `{case['label']}` | `{case['conductivity_over_diffusivity']:.4f}` | `{case['median_relative_error']:.4f}` | `{case['p90_relative_error']:.4f}` | `{case['shell_kappa_cv']:.4f}` | `{case['metric_tracking_corr']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
            for case in family_cases
        )
        sections.append(
            f"""## {DEFORMATION_LABELS[family].title()}

| case | `kappa / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | metric corr | status |
| --- | --- | --- | --- | --- | --- | --- |
{rows}
"""
        )

    drift_rows = "\n".join(
        f"| `{row['label']}` | `{row['drift']:.4f}` |" for row in result["refinement_drifts"]
    )
    note = f"""# Scalar Kernel Graph Current Closure Inhomogeneity v1

## Purpose

Test the next honest follow-on after the clean `52.4` shell-native current closure:

> do the extracted metric and the shell-native transport law track a designed scalar inhomogeneity together on the same carrier?

The carrier stays fixed:

- same scalar kernel graph family
- same local metric-field construction
- same shell-native current law

The only new ingredient is a small deterministic deformation profile on the point cloud.

## Construction

- families: `{result['config']['families']}`
- amplitudes: `{result['config']['eta_values']}`
- clean refinements: `{result['config']['n_sides']}`

For each case:

1. deform the same clean carrier by one designed radial profile;
2. rebuild the scalar operator on the deformed point cloud;
3. extract the coarse local metric field;
4. test whether the shell-native transient current closure still fits one bounded constitutive family;
5. compare the shellwise metric shift against the imposed deformation profile.

Artifacts:

- script: `numerics/simulations/scalar_kernel_graph_current_closure_inhomogeneity.py`
- config: `config.json -> scalar_kernel_graph_current_closure_inhomogeneity`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_current_closure_inhomogeneity_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

{''.join(sections)}
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
        handle.write("\n## scalar kernel graph current closure inhomogeneity\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, families={config['families']}, eta_values={config['eta_values']}, "
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
