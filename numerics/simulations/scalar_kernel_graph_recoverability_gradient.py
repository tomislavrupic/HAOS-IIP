#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from scalar_kernel_graph_geometry_bridge import solve_green_field
from scalar_kernel_graph_geometry_robustness import (
    build_cubic_points,
    build_jittered_points,
    build_point_cloud_operator,
)
from scalar_kernel_graph_metric_field import bulk_mask, coarse_grain_metric_field, local_metric_field

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Recoverability_Gradient_v1.md"
for path in (DATA, PLOTS):
    path.mkdir(exist_ok=True)

MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    "epsilon_coeff": 0.5,
    "cutoff_factor": 2.5,
    "disorder_n_sides": [11, 13],
    "jitter_fractions": [0.0, 0.02, 0.04],
    "jitter_seed": 42,
    "kernel_n_sides": [11, 13],
    "kernel_families": ["gaussian_local", "gaussian_half", "inverse_quadratic"],
    "bulk_margin_layers": 2,
    "metric_radius_factor": 2.5,
    "gradient_radius_factor": 2.5,
    "fit_r_max": 0.55,
    "thresholds": {
        "min_radial_alignment": 0.94,
        "max_power_deviation": 0.18,
        "min_power_fit_r2": 0.95,
        "max_flux_constancy_cv": 0.12,
        "max_refinement_profile_drift": 0.08,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_recoverability_gradient"), dict):
            merged.update({k: v for k, v in raw["scalar_kernel_graph_recoverability_gradient"].items() if k != "thresholds"})
            if isinstance(raw["scalar_kernel_graph_recoverability_gradient"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_recoverability_gradient"]["thresholds"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def fit_power(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = (x > 0.0) & (y > 0.0)
    log_x = np.log(x[mask])
    log_y = np.log(y[mask])
    slope, intercept = np.polyfit(log_x, log_y, 1)
    pred = slope * log_x + intercept
    ss_res = float(np.sum((log_y - pred) ** 2))
    ss_tot = float(np.sum((log_y - np.mean(log_y)) ** 2))
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
    }


def local_gradient_field(points: np.ndarray, values: np.ndarray, radius: float) -> np.ndarray:
    tree = cKDTree(points)
    gradients = np.zeros((points.shape[0], 3), dtype=float)
    for idx, point in enumerate(points):
        neighborhood = [j for j in tree.query_ball_point(point, radius) if j != idx]
        if len(neighborhood) < 4:
            continue
        deltas = points[neighborhood] - point
        rhs = values[neighborhood] - values[idx]
        grad, *_ = np.linalg.lstsq(deltas, rhs, rcond=None)
        gradients[idx] = grad
    return gradients


def shell_profile(
    bulk_points: np.ndarray,
    source_point: np.ndarray,
    radial_flux: np.ndarray,
    h: float,
    fit_r_max: float,
) -> dict[str, Any]:
    radii = np.linalg.norm(bulk_points - source_point, axis=1)
    labels = np.round(radii / max(h, 1.0e-12), 8) * h
    unique = np.unique(labels)
    rows = []
    for radius in unique:
        mask = np.isclose(labels, radius)
        if np.sum(mask) == 0:
            continue
        mean_r = float(np.mean(radii[mask]))
        if mean_r <= 0.0 or mean_r > float(fit_r_max):
            continue
        mean_flux = float(np.mean(radial_flux[mask]))
        rows.append(
            {
                "radius": mean_r,
                "radial_flux": mean_flux,
                "scaled_flux": float((mean_r**2) * mean_flux),
                "count": int(np.sum(mask)),
            }
        )
    radii_fit = np.asarray([row["radius"] for row in rows], dtype=float)
    flux_fit = np.asarray([abs(row["radial_flux"]) for row in rows], dtype=float)
    scaled_flux = np.asarray([row["scaled_flux"] for row in rows], dtype=float)
    power = fit_power(radii_fit, flux_fit)
    return {
        "rows": rows,
        "power_slope": float(power["slope"]),
        "power_intercept": float(power["intercept"]),
        "power_fit_r2": float(power["r2"]),
        "flux_constancy_cv": float(np.std(scaled_flux) / max(np.mean(np.abs(scaled_flux)), 1.0e-18)),
    }


def evaluate_case(
    label: str,
    points: np.ndarray,
    source_idx: int,
    h: float,
    epsilon_coeff: float,
    cutoff_factor: float,
    kernel_family: str,
    bulk_margin_layers: int,
    metric_radius_factor: float,
    gradient_radius_factor: float,
    fit_r_max: float,
) -> dict[str, Any]:
    delta, epsilon_k = build_point_cloud_operator(
        points=points,
        h=h,
        epsilon_coeff=epsilon_coeff,
        cutoff_factor=cutoff_factor,
        kernel_family=kernel_family,
    )
    phi = solve_green_field(delta=delta, source_idx=source_idx)
    raw_metric = local_metric_field(points=points, delta=delta)
    mask = bulk_mask(points=points, h=h, bulk_margin_layers=bulk_margin_layers)
    bulk_points = points[mask]
    bulk_metric = coarse_grain_metric_field(points=bulk_points, field=raw_metric[mask], radius=float(metric_radius_factor) * h)
    bulk_phi = phi[mask]
    bulk_grad = local_gradient_field(points=bulk_points, values=bulk_phi, radius=float(gradient_radius_factor) * h)
    bulk_flux = -np.einsum("nij,nj->ni", bulk_metric, bulk_grad)
    source_point = points[source_idx]
    radii = np.linalg.norm(bulk_points - source_point, axis=1)
    radial_dirs = np.zeros_like(bulk_flux)
    nonzero = radii > 1.0e-12
    radial_dirs[nonzero] = (bulk_points[nonzero] - source_point) / radii[nonzero, None]
    radial_component = np.einsum("ni,ni->n", bulk_flux, radial_dirs)
    orientation = 1.0 if float(np.mean(radial_component[nonzero])) >= 0.0 else -1.0
    radial_component = orientation * radial_component
    flux_norm = np.linalg.norm(bulk_flux, axis=1)
    radial_alignment = float(
        np.mean(
            np.clip(
                radial_component[nonzero] / np.maximum(flux_norm[nonzero], 1.0e-18),
                0.0,
                1.0,
            )
        )
    )
    profile = shell_profile(
        bulk_points=bulk_points,
        source_point=source_point,
        radial_flux=radial_component,
        h=h,
        fit_r_max=fit_r_max,
    )
    return {
        "label": label,
        "points_count": int(points.shape[0]),
        "lattice_spacing": float(h),
        "epsilon_coeff": float(epsilon_coeff),
        "epsilon_k": float(epsilon_k),
        "kernel_family": kernel_family,
        "radial_alignment": radial_alignment,
        "profile": profile,
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        case["radial_alignment"] >= float(thresholds["min_radial_alignment"])
        and abs(case["profile"]["power_slope"] + 2.0) <= float(thresholds["max_power_deviation"])
        and case["profile"]["power_fit_r2"] >= float(thresholds["min_power_fit_r2"])
        and case["profile"]["flux_constancy_cv"] <= float(thresholds["max_flux_constancy_cv"])
    )


def profile_drift(case_a: dict[str, Any], case_b: dict[str, Any]) -> float:
    rows_a = case_a["profile"]["rows"]
    rows_b = case_b["profile"]["rows"]
    min_len = min(len(rows_a), len(rows_b))
    flux_a = np.asarray([row["scaled_flux"] for row in rows_a[:min_len]], dtype=float)
    flux_b = np.asarray([row["scaled_flux"] for row in rows_b[:min_len]], dtype=float)
    norm = max(float(np.mean(np.abs(flux_a))), 1.0e-18)
    return float(np.linalg.norm(flux_a - flux_b) / (math.sqrt(min_len) * norm))


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    thresholds = cfg["thresholds"]

    disorder_cases: list[dict[str, Any]] = []
    for n_side in [int(v) for v in cfg["disorder_n_sides"]]:
        for jitter_fraction in [float(v) for v in cfg["jitter_fractions"]]:
            points, source_idx, h = build_jittered_points(
                n_side=n_side,
                jitter_fraction=jitter_fraction,
                seed=int(cfg["jitter_seed"]),
            )
            case = evaluate_case(
                label=f"disorder|n={n_side}|jitter={jitter_fraction:.3f}",
                points=points,
                source_idx=source_idx,
                h=h,
                epsilon_coeff=float(cfg["epsilon_coeff"]),
                cutoff_factor=float(cfg["cutoff_factor"]),
                kernel_family="gaussian_local",
                bulk_margin_layers=int(cfg["bulk_margin_layers"]),
                metric_radius_factor=float(cfg["metric_radius_factor"]),
                gradient_radius_factor=float(cfg["gradient_radius_factor"]),
                fit_r_max=float(cfg["fit_r_max"]),
            )
            case["jitter_fraction"] = float(jitter_fraction)
            case["pass"] = case_passes(case, thresholds)
            disorder_cases.append(case)

    kernel_cases: list[dict[str, Any]] = []
    for n_side in [int(v) for v in cfg["kernel_n_sides"]]:
        points, source_idx, h = build_cubic_points(n_side)
        for family in [str(v) for v in cfg["kernel_families"]]:
            case = evaluate_case(
                label=f"kernel|n={n_side}|family={family}",
                points=points,
                source_idx=source_idx,
                h=h,
                epsilon_coeff=float(cfg["epsilon_coeff"]),
                cutoff_factor=float(cfg["cutoff_factor"]),
                kernel_family=family,
                bulk_margin_layers=int(cfg["bulk_margin_layers"]),
                metric_radius_factor=float(cfg["metric_radius_factor"]),
                gradient_radius_factor=float(cfg["gradient_radius_factor"]),
                fit_r_max=float(cfg["fit_r_max"]),
            )
            case["pass"] = case_passes(case, thresholds)
            kernel_cases.append(case)

    disorder_drifts = []
    for jitter_fraction in [float(v) for v in cfg["jitter_fractions"]]:
        pair = [case for case in disorder_cases if abs(case["jitter_fraction"] - jitter_fraction) < 1.0e-12]
        if len(pair) == 2:
            pair = sorted(pair, key=lambda item: item["points_count"])
            disorder_drifts.append({"label": f"disorder|jitter={jitter_fraction:.3f}", "drift": profile_drift(pair[0], pair[1])})

    kernel_drifts = []
    for family in [str(v) for v in cfg["kernel_families"]]:
        pair = [case for case in kernel_cases if case["kernel_family"] == family]
        if len(pair) == 2:
            pair = sorted(pair, key=lambda item: item["points_count"])
            kernel_drifts.append({"label": f"kernel|family={family}", "drift": profile_drift(pair[0], pair[1])})

    max_profile_drift = max([row["drift"] for row in disorder_drifts + kernel_drifts], default=0.0)
    all_cases_pass = all(case["pass"] for case in disorder_cases + kernel_cases)
    drift_pass = max_profile_drift <= float(thresholds["max_refinement_profile_drift"])
    weakest_case = max(
        disorder_cases + kernel_cases,
        key=lambda item: (
            abs(item["profile"]["power_slope"] + 2.0),
            item["profile"]["flux_constancy_cv"],
            1.0 - item["radial_alignment"],
        ),
    )

    if all_cases_pass and drift_pass:
        observation = (
            "the scalar Green potential now supports a stable effective recoverability-gradient field on the same validated geometry carrier: "
            "after local metric and gradient coarse-graining, the response remains strongly radial and its shell profile follows the expected inverse-square law on the tested window"
        )
        conclusion = (
            f"all {len(disorder_cases)} disorder cases and all {len(kernel_cases)} kernel-family cases pass the recoverability-gradient thresholds, and the maximum refinement-profile drift is {max_profile_drift:.4f}; "
            f"the weakest passing case is `{weakest_case['label']}` with radial alignment {weakest_case['radial_alignment']:.4f}, "
            f"power slope {weakest_case['profile']['power_slope']:.4f}, and flux constancy CV {weakest_case['profile']['flux_constancy_cv']:.4f}; "
            "this gives the first bounded operator-native force-like read on the scalar carrier without introducing a new substrate or a stochastic recovery model"
        )
    else:
        failing_cases = [case["label"] for case in disorder_cases + kernel_cases if not case["pass"]]
        observation = (
            "the scalar Green potential supports a geometry carrier and a local metric field, but the recoverability-gradient response does not yet close as one stable inverse-square family across the tested window"
        )
        conclusion = (
            f"open recoverability-gradient cases remain: failing_cases={failing_cases}, max_profile_drift={max_profile_drift:.4f}; "
            "the correct boundary therefore stays at the scalar geometry carrier and local metric-field statement without yet promoting a positive force-like response claim"
        )

    return {
        "experiment": "scalar_kernel_graph_recoverability_gradient",
        "config": cfg,
        "disorder_cases": disorder_cases,
        "kernel_cases": kernel_cases,
        "disorder_profile_drifts": disorder_drifts,
        "kernel_profile_drifts": kernel_drifts,
        "max_profile_drift": float(max_profile_drift),
        "observation": observation,
        "conclusion": conclusion,
    }


def make_summary_plot(disorder_cases: list[dict[str, Any]], kernel_cases: list[dict[str, Any]], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_recoverability_gradient_summary.png"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for n_side in sorted({int(case["label"].split("|")[1].split("=")[1]) for case in disorder_cases}):
        subset = sorted(
            [case for case in disorder_cases if int(case["label"].split("|")[1].split("=")[1]) == n_side],
            key=lambda item: float(item["jitter_fraction"]),
        )
        x = [float(case["jitter_fraction"]) for case in subset]
        axes[0, 0].plot(x, [case["radial_alignment"] for case in subset], marker="o", label=f"n={n_side}")
        axes[0, 1].plot(x, [case["profile"]["power_slope"] for case in subset], marker="o", label=f"n={n_side}")
        axes[1, 0].plot(x, [case["profile"]["flux_constancy_cv"] for case in subset], marker="o", label=f"n={n_side}")

    family_order = [str(v) for v in DEFAULT_CONFIG["kernel_families"]]
    for n_side in sorted({int(case["label"].split("|")[1].split("=")[1]) for case in kernel_cases}):
        subset = sorted(
            [case for case in kernel_cases if int(case["label"].split("|")[1].split("=")[1]) == n_side],
            key=lambda item: family_order.index(str(item["kernel_family"])),
        )
        x = range(len(subset))
        axes[1, 1].plot(x, [case["profile"]["power_fit_r2"] for case in subset], marker="o", label=f"n={n_side}")
        axes[1, 1].set_xticks(list(x), labels=[case["kernel_family"] for case in subset], rotation=20)

    axes[0, 0].set_title("Radial alignment")
    axes[0, 0].set_xlabel("jitter fraction")
    axes[0, 0].set_ylim(0.9, 1.001)
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].set_title("Power-law slope")
    axes[0, 1].set_xlabel("jitter fraction")
    axes[0, 1].axhline(-2.0, color="k", linestyle="--", linewidth=1.0)
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].set_title("Scaled-flux constancy CV")
    axes[1, 0].set_xlabel("jitter fraction")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].set_title("Power-law fit $R^2$")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend()

    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def make_profile_plot(result: dict[str, Any], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_recoverability_gradient_profiles.png"
    selected = [
        [case for case in result["disorder_cases"] if abs(case["jitter_fraction"] - 0.04) < 1.0e-12][-1],
        [case for case in result["kernel_cases"] if case["kernel_family"] == "gaussian_local"][-1],
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for axis, case in zip(axes, selected):
        rows = case["profile"]["rows"]
        radii = np.asarray([row["radius"] for row in rows], dtype=float)
        flux = np.asarray([abs(row["radial_flux"]) for row in rows], dtype=float)
        scaled = np.asarray([row["scaled_flux"] for row in rows], dtype=float)
        axis.loglog(radii, flux, "o", label="|radial flux|")
        axis.loglog(radii, np.exp(case["profile"]["power_intercept"]) * radii ** case["profile"]["power_slope"], "--", label="power fit")
        axis2 = axis.twinx()
        axis2.plot(radii, scaled, color="tab:red", alpha=0.5, label="$r^2 F_r$")
        axis.set_title(case["label"])
        axis.set_xlabel("r")
        axis.set_ylabel("|$F_r$|")
        axis2.set_ylabel("$r^2 F_r$")
        axis.grid(alpha=0.25)
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_recoverability_gradient.json"
    latest = DATA / "scalar_kernel_graph_recoverability_gradient_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    disorder_rows = "\n".join(
        f"| `{case['label']}` | `{case['radial_alignment']:.4f}` | `{case['profile']['power_slope']:.4f}` | `{case['profile']['power_fit_r2']:.4f}` | `{case['profile']['flux_constancy_cv']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["disorder_cases"]
    )
    kernel_rows = "\n".join(
        f"| `{case['label']}` | `{case['radial_alignment']:.4f}` | `{case['profile']['power_slope']:.4f}` | `{case['profile']['power_fit_r2']:.4f}` | `{case['profile']['flux_constancy_cv']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["kernel_cases"]
    )
    drift_rows = "\n".join(
        f"| `{row['label']}` | `{row['drift']:.4f}` |"
        for row in result["disorder_profile_drifts"] + result["kernel_profile_drifts"]
    )
    note = f"""# Scalar Kernel Graph Recoverability Gradient v1

## Purpose

Build the first bounded `52`-line force-like observable on the validated scalar geometry carrier without importing a new stochastic recovery model.

The deterministic proxy used here is simple:

- recoverability potential: the same mean-zero scalar Green field already validated on the carrier
- local metric: the coarse local metric field from `51.4`
- effective recoverability-gradient field:

$$
J_{{rec}} = - G \\nabla \\phi
$$

This is not yet called a fundamental force. It is the first bounded operator-native response field on the same carrier.

## Construction

- scalar carrier: same `3D` kernel-graph family
- metric radius factor: `{result['config']['metric_radius_factor']}`
- gradient radius factor: `{result['config']['gradient_radius_factor']}`
- bulk margin layers: `{result['config']['bulk_margin_layers']}`

The bounded test asks whether the coarse response field is:

1. strongly radial
2. close to inverse-square in shell profile
3. close to constant in `r^2 F_r`
4. stable under the same mild disorder and bounded kernel-family window

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_recoverability_gradient.py`
- config: `config.json -> scalar_kernel_graph_recoverability_gradient`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_recoverability_gradient_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Disorder pass

| case | radial alignment | power slope | power fit `R^2` | scaled-flux CV | status |
| --- | --- | --- | --- | --- | --- |
{disorder_rows}

## Kernel-family pass

| case | radial alignment | power slope | power fit `R^2` | scaled-flux CV | status |
| --- | --- | --- | --- | --- | --- |
{kernel_rows}

## Refinement profile drift

| pair | drift |
| --- | --- |
{drift_rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Current boundary

If this note is read positively, the correct bounded statement is:

> on the validated scalar geometry carrier, the Green recoverability potential supports a stable inverse-square-like effective response field under the tested mild disorder and bounded kernel-family window.

This note still does **not** claim:

- a universal force law
- a fundamental ontic force
- closure on arbitrary substrates or strong disorder
- coupling to currents, curvature, or spacetime
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph recoverability gradient\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"disorder_n_sides={config['disorder_n_sides']}, "
            f"jitter_fractions={config['jitter_fractions']}, "
            f"kernel_n_sides={config['kernel_n_sides']}, "
            f"kernel_families={config['kernel_families']}, "
            f"metric_radius_factor={config['metric_radius_factor']}, "
            f"gradient_radius_factor={config['gradient_radius_factor']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write(f"- Plots: `{plot_paths[0]}`, `{plot_paths[1]}`\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def main() -> None:
    result = make_result()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_paths = [
        make_summary_plot(result["disorder_cases"], result["kernel_cases"], stamp),
        make_profile_plot(result, stamp),
    ]
    result_path, latest_path = save_results(result, stamp)
    write_note(result, result_path=result_path, plot_paths=plot_paths)
    append_log(
        result_path=result_path,
        plot_paths=plot_paths,
        config=result["config"],
        observation=result["observation"],
        conclusion=result["conclusion"],
    )
    print(json.dumps({"result_path": result_path, "latest_path": latest_path, "plots": plot_paths}, indent=2))


if __name__ == "__main__":
    main()
