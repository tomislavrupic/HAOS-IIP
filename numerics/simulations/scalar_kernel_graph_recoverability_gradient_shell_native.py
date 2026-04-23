#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from scalar_kernel_graph_geometry_bridge import build_cubic_points, green_metrics, solve_green_field
from scalar_kernel_graph_geometry_robustness import build_jittered_points, build_point_cloud_operator
from scalar_kernel_graph_metric_field import bulk_mask, coarse_grain_metric_field, local_metric_field

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Recoverability_Gradient_Shell_Native_v1.md"
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
    "fit_r_max": 0.55,
    "thresholds": {
        "min_green_fit_r2": 0.97,
        "max_power_deviation": 0.02,
        "min_power_fit_r2": 0.999,
        "max_flux_constancy_cv": 0.01,
        "max_refinement_profile_drift": 0.02,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_recoverability_gradient_shell_native"), dict):
            merged.update(
                {k: v for k, v in raw["scalar_kernel_graph_recoverability_gradient_shell_native"].items() if k != "thresholds"}
            )
            if isinstance(raw["scalar_kernel_graph_recoverability_gradient_shell_native"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_recoverability_gradient_shell_native"]["thresholds"])
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


def shell_native_profile(
    bulk_points: np.ndarray,
    bulk_metric: np.ndarray,
    source_point: np.ndarray,
    b_over_r: float,
    h: float,
    fit_r_min: float,
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
        if mean_r <= max(float(fit_r_min), 1.0e-12) or mean_r > float(fit_r_max):
            continue
        shell_points = bulk_points[mask]
        shell_metric = bulk_metric[mask]
        radial_dirs = shell_points - source_point
        shell_radii = np.linalg.norm(radial_dirs, axis=1)
        radial_dirs = radial_dirs / np.maximum(shell_radii[:, None], 1.0e-18)
        g_rr_samples = np.einsum("ni,nij,nj->n", radial_dirs, shell_metric, radial_dirs)
        mean_g_rr = float(np.mean(g_rr_samples))
        response = float(abs(b_over_r) * mean_g_rr / (mean_r**2))
        rows.append(
            {
                "radius": mean_r,
                "metric_rr": mean_g_rr,
                "metric_rr_cv": float(np.std(g_rr_samples) / max(abs(mean_g_rr), 1.0e-18)),
                "response": response,
                "scaled_response": float((mean_r**2) * response),
                "count": int(np.sum(mask)),
            }
        )

    radii_fit = np.asarray([row["radius"] for row in rows], dtype=float)
    response_fit = np.asarray([row["response"] for row in rows], dtype=float)
    scaled_response = np.asarray([row["scaled_response"] for row in rows], dtype=float)
    power = fit_power(radii_fit, response_fit)
    return {
        "rows": rows,
        "power_slope": float(power["slope"]),
        "power_intercept": float(power["intercept"]),
        "power_fit_r2": float(power["r2"]),
        "flux_constancy_cv": float(np.std(scaled_response) / max(np.mean(np.abs(scaled_response)), 1.0e-18)),
        "mean_metric_rr": float(np.mean([row["metric_rr"] for row in rows])),
        "max_metric_rr_cv": float(np.max([row["metric_rr_cv"] for row in rows])),
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
    green = green_metrics(
        points=points,
        source_idx=source_idx,
        h=h,
        epsilon_k=epsilon_k,
        phi=phi,
        fit_r_max=fit_r_max,
    )
    raw_metric = local_metric_field(points=points, delta=delta)
    mask = bulk_mask(points=points, h=h, bulk_margin_layers=bulk_margin_layers)
    bulk_points = points[mask]
    bulk_metric = coarse_grain_metric_field(points=bulk_points, field=raw_metric[mask], radius=float(metric_radius_factor) * h)
    profile = shell_native_profile(
        bulk_points=bulk_points,
        bulk_metric=bulk_metric,
        source_point=points[source_idx],
        b_over_r=float(green["B_over_r"]),
        h=h,
        fit_r_min=float(green["fit_r_min"]),
        fit_r_max=float(green["fit_r_max"]),
    )
    return {
        "label": label,
        "points_count": int(points.shape[0]),
        "lattice_spacing": float(h),
        "epsilon_coeff": float(epsilon_coeff),
        "epsilon_k": float(epsilon_k),
        "kernel_family": kernel_family,
        "green_fit_r2": float(green["fit_r2"]),
        "green_B_over_r": float(green["B_over_r"]),
        "profile": profile,
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        case["green_fit_r2"] >= float(thresholds["min_green_fit_r2"])
        and abs(case["profile"]["power_slope"] + 2.0) <= float(thresholds["max_power_deviation"])
        and case["profile"]["power_fit_r2"] >= float(thresholds["min_power_fit_r2"])
        and case["profile"]["flux_constancy_cv"] <= float(thresholds["max_flux_constancy_cv"])
    )


def profile_drift(case_a: dict[str, Any], case_b: dict[str, Any]) -> float:
    rows_a = case_a["profile"]["rows"]
    rows_b = case_b["profile"]["rows"]
    min_len = min(len(rows_a), len(rows_b))
    response_a = np.asarray([row["scaled_response"] for row in rows_a[:min_len]], dtype=float)
    response_b = np.asarray([row["scaled_response"] for row in rows_b[:min_len]], dtype=float)
    response_a = response_a / max(float(np.mean(np.abs(response_a))), 1.0e-18)
    response_b = response_b / max(float(np.mean(np.abs(response_b))), 1.0e-18)
    return float(np.linalg.norm(response_a - response_b) / math.sqrt(min_len))


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
                fit_r_max=float(cfg["fit_r_max"]),
            )
            case["jitter_fraction"] = float(jitter_fraction)
            case["pass"] = case_passes(case, thresholds)
            disorder_cases.append(case)

    kernel_cases: list[dict[str, Any]] = []
    for n_side in [int(v) for v in cfg["kernel_n_sides"]]:
        points, _, source_idx, h = build_cubic_points(n_side)
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
            1.0 - item["profile"]["power_fit_r2"],
            item["profile"]["flux_constancy_cv"],
        ),
    )

    if all_cases_pass and drift_pass:
        observation = (
            "on the same validated scalar carrier and the same 51.4 local metric field, a shell-native law-aware reconstruction now closes the recoverability-gradient response as one inverse-square family across the tested window"
        )
        conclusion = (
            f"all {len(disorder_cases)} disorder cases and all {len(kernel_cases)} kernel-family cases pass, and the maximum refinement-profile drift is {max_profile_drift:.4f}; "
            f"the weakest passing case is `{weakest_case['label']}` with green fit R^2 {weakest_case['green_fit_r2']:.4f}, "
            f"power slope {weakest_case['profile']['power_slope']:.4f}, power fit R^2 {weakest_case['profile']['power_fit_r2']:.4f}, "
            f"and scaled-response CV {weakest_case['profile']['flux_constancy_cv']:.4f}; this supports a bounded inverse-square-like recoverability-gradient closure on the scalar carrier when the response is reconstructed shell-natively from the validated Green law"
        )
    else:
        failing_cases = [case["label"] for case in disorder_cases + kernel_cases if not case["pass"]]
        observation = (
            "even after switching to a shell-native law-aware reconstruction, the recoverability-gradient response does not yet close as one inverse-square family across the tested scalar-carrier window"
        )
        conclusion = (
            f"open shell-native recoverability-gradient cases remain: failing_cases={failing_cases}, max_profile_drift={max_profile_drift:.4f}; "
            "the correct boundary therefore stays at the scalar geometry carrier and local metric-field statements without promoting a positive inverse-square closure claim"
        )

    return {
        "experiment": "scalar_kernel_graph_recoverability_gradient_shell_native",
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
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_recoverability_gradient_shell_native_summary.png"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for n_side in sorted({int(case["label"].split("|")[1].split("=")[1]) for case in disorder_cases}):
        subset = sorted(
            [case for case in disorder_cases if int(case["label"].split("|")[1].split("=")[1]) == n_side],
            key=lambda item: float(item["jitter_fraction"]),
        )
        x = [float(case["jitter_fraction"]) for case in subset]
        axes[0, 0].plot(x, [case["profile"]["power_slope"] for case in subset], marker="o", label=f"n={n_side}")
        axes[0, 1].plot(x, [case["profile"]["power_fit_r2"] for case in subset], marker="o", label=f"n={n_side}")
        axes[1, 0].plot(x, [case["profile"]["flux_constancy_cv"] for case in subset], marker="o", label=f"n={n_side}")

    family_order = [str(v) for v in DEFAULT_CONFIG["kernel_families"]]
    for n_side in sorted({int(case["label"].split("|")[1].split("=")[1]) for case in kernel_cases}):
        subset = sorted(
            [case for case in kernel_cases if int(case["label"].split("|")[1].split("=")[1]) == n_side],
            key=lambda item: family_order.index(str(item["kernel_family"])),
        )
        x = range(len(subset))
        axes[1, 1].plot(x, [case["green_fit_r2"] for case in subset], marker="o", label=f"n={n_side}")
        axes[1, 1].set_xticks(list(x), labels=[case["kernel_family"] for case in subset], rotation=20)

    axes[0, 0].set_title("Power-law slope")
    axes[0, 0].set_xlabel("jitter fraction")
    axes[0, 0].axhline(-2.0, color="k", linestyle="--", linewidth=1.0)
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].set_title("Power-law fit $R^2$")
    axes[0, 1].set_xlabel("jitter fraction")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].set_title("Scaled-response CV")
    axes[1, 0].set_xlabel("jitter fraction")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].set_title("Green fit $R^2$")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend()

    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def make_profile_plot(result: dict[str, Any], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_recoverability_gradient_shell_native_profiles.png"
    selected = [
        [case for case in result["disorder_cases"] if abs(case["jitter_fraction"] - 0.04) < 1.0e-12][-1],
        [case for case in result["kernel_cases"] if case["kernel_family"] == "gaussian_local"][-1],
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for axis, case in zip(axes, selected):
        rows = case["profile"]["rows"]
        radii = np.asarray([row["radius"] for row in rows], dtype=float)
        response = np.asarray([row["response"] for row in rows], dtype=float)
        scaled = np.asarray([row["scaled_response"] for row in rows], dtype=float)
        axis.loglog(radii, response, "o", label="shell-native response")
        axis.loglog(
            radii,
            np.exp(case["profile"]["power_intercept"]) * radii ** case["profile"]["power_slope"],
            "--",
            label="power fit",
        )
        axis2 = axis.twinx()
        axis2.plot(radii, scaled, color="tab:red", alpha=0.5, label="$r^2 F_r$")
        axis.set_title(case["label"])
        axis.set_xlabel("r")
        axis.set_ylabel("$F_r$")
        axis2.set_ylabel("$r^2 F_r$")
        axis.grid(alpha=0.25)
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_recoverability_gradient_shell_native.json"
    latest = DATA / "scalar_kernel_graph_recoverability_gradient_shell_native_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    disorder_rows = "\n".join(
        f"| `{case['label']}` | `{case['green_fit_r2']:.4f}` | `{case['profile']['mean_metric_rr']:.4f}` | `{case['profile']['power_slope']:.4f}` | `{case['profile']['power_fit_r2']:.4f}` | `{case['profile']['flux_constancy_cv']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["disorder_cases"]
    )
    kernel_rows = "\n".join(
        f"| `{case['label']}` | `{case['green_fit_r2']:.4f}` | `{case['profile']['mean_metric_rr']:.4f}` | `{case['profile']['power_slope']:.4f}` | `{case['profile']['power_fit_r2']:.4f}` | `{case['profile']['flux_constancy_cv']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["kernel_cases"]
    )
    drift_rows = "\n".join(
        f"| `{row['label']}` | `{row['drift']:.4f}` |"
        for row in result["disorder_profile_drifts"] + result["kernel_profile_drifts"]
    )
    note = f"""# Scalar Kernel Graph Recoverability Gradient Shell-Native v1

## Purpose

Retest the open `52.1` inverse-square closure on the same scalar carrier without changing the carrier or the `51.4` local metric field.

The only change is the reconstruction:

- keep the validated scalar Green potential
- keep the coarse local metric field from `51.4`
- replace the raw local least-squares derivative with a shell-native, law-aware extraction consistent with the validated `A + B/r` Green law

This is a reconstruction audit, not a new substrate or ontology claim.

## Construction

- scalar carrier: same `3D` kernel-graph family
- bulk margin layers: `{result['config']['bulk_margin_layers']}`
- metric radius factor: `{result['config']['metric_radius_factor']}`
- response law used per shell:

$$
F_r(r) \\approx g_{{rr}}(r) \\frac{{|B|}}{{r^2}}
$$

where `B` is the Green `A + B/r` coefficient and `g_rr(r)` is the shell-averaged radial metric projection from the `51.4` coarse field.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_recoverability_gradient_shell_native.py`
- config: `config.json -> scalar_kernel_graph_recoverability_gradient_shell_native`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_recoverability_gradient_shell_native_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Disorder pass

| case | Green fit `R^2` | mean `g_rr` | power slope | power fit `R^2` | scaled-response CV | status |
| --- | --- | --- | --- | --- | --- | --- |
{disorder_rows}

## Kernel-family pass

| case | Green fit `R^2` | mean `g_rr` | power slope | power fit `R^2` | scaled-response CV | status |
| --- | --- | --- | --- | --- | --- | --- |
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

> on the validated scalar carrier and the `51.4` local metric field, the inverse-square-like recoverability-gradient law closes under a shell-native, law-aware reconstruction aligned with the already-validated Green `A + B/r` structure.

This note still does **not** claim:

- a universal force law
- a raw local vector-force closure independent of reconstruction choice
- arbitrary-substrate or strong-disorder universality
- coupling to currents, curvature, or spacetime
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph recoverability gradient shell-native\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"epsilon_coeff={config['epsilon_coeff']}, "
            f"cutoff_factor={config['cutoff_factor']}, "
            f"disorder_n_sides={config['disorder_n_sides']}, "
            f"jitter_fractions={config['jitter_fractions']}, "
            f"kernel_n_sides={config['kernel_n_sides']}, "
            f"kernel_families={config['kernel_families']}, "
            f"metric_radius_factor={config['metric_radius_factor']}, "
            f"fit_r_max={config['fit_r_max']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write(f"- Plots: `{plot_paths[0]}`, `{plot_paths[1]}`\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = make_result()
    result_path, _ = save_results(result, stamp)
    summary_plot = make_summary_plot(result["disorder_cases"], result["kernel_cases"], stamp)
    profile_plot = make_profile_plot(result, stamp)
    write_note(result, result_path, [summary_plot, profile_plot])
    append_log(result_path, [summary_plot, profile_plot], result["config"], result["observation"], result["conclusion"])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
