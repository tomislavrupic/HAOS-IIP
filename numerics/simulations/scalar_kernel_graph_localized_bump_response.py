#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from scalar_kernel_graph_current_closure_inhomogeneity import (
    case_passes,
    evaluate_baseline,
    evaluate_case,
    normalized_profile,
    profile_drift,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Localized_Bump_Response_v1.md"
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
    "weak_eta_values": [0.01, 0.02, 0.03, 0.04, 0.05],
    "stress_eta_values": [0.075, 0.10, 0.15],
    "epsilon_coeff": 0.5,
    "cutoff_factor": 2.5,
    "heat_time_max": 0.03,
    "heat_time_steps": 61,
    "fit_window": [0.010, 0.026],
    "bulk_margin_layers": 2,
    "metric_radius_factor": 2.5,
    "shell_r_max": 0.45,
    "bump_sigma": 0.18,
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
        if isinstance(raw.get("scalar_kernel_graph_localized_bump_response"), dict):
            merged.update(
                {
                    key: value
                    for key, value in raw["scalar_kernel_graph_localized_bump_response"].items()
                    if key != "thresholds"
                }
            )
            if isinstance(raw["scalar_kernel_graph_localized_bump_response"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_localized_bump_response"]["thresholds"])
    if config is not None:
        merged.update({key: value for key, value in config.items() if key != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def evaluate_regime(
    regime: str,
    eta_values: list[float],
    baseline_cases: dict[int, dict[str, Any]],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for eta in eta_values:
        for n_side in [int(value) for value in config["n_sides"]]:
            case = evaluate_case(
                n_side=n_side,
                family="bump",
                eta=float(eta),
                epsilon_coeff=float(config["epsilon_coeff"]),
                cutoff_factor=float(config["cutoff_factor"]),
                heat_time_max=float(config["heat_time_max"]),
                heat_time_steps=int(config["heat_time_steps"]),
                fit_window=[float(value) for value in config["fit_window"]],
                bulk_margin_layers=int(config["bulk_margin_layers"]),
                metric_radius_factor=float(config["metric_radius_factor"]),
                shell_r_max=float(config["shell_r_max"]),
                baseline_shell_metric_rr=np.asarray(baseline_cases[n_side]["shell_metric_rr"], dtype=float),
            )
            case["regime"] = regime
            case["label"] = f"{regime}_bump|n={n_side}|eta={float(eta):.3f}"
            case["pass"] = case_passes(case, config["thresholds"])
            cases.append(case)
    return cases


def refinement_drifts(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for regime in sorted({case["regime"] for case in cases}):
        regime_cases = [case for case in cases if case["regime"] == regime]
        for eta in sorted({float(case["eta"]) for case in regime_cases}):
            pair = sorted(
                [case for case in regime_cases if abs(float(case["eta"]) - eta) < 1.0e-12],
                key=lambda item: item["n_side"],
            )
            if len(pair) == 2:
                out.append(
                    {
                        "label": f"{regime}_bump|eta={eta:.3f}",
                        "regime": regime,
                        "eta": float(eta),
                        "drift": profile_drift(
                            np.asarray(pair[0]["normalized_shell_kappa_profile"], dtype=float),
                            np.asarray(pair[1]["normalized_shell_kappa_profile"], dtype=float),
                        ),
                    }
                )
    return out


def summarize_regime(regime: str, cases: list[dict[str, Any]], drifts: list[dict[str, Any]], thresholds: dict[str, float]) -> dict[str, Any]:
    regime_cases = [case for case in cases if case["regime"] == regime]
    regime_drifts = [row for row in drifts if row["regime"] == regime]
    max_drift = max((float(row["drift"]) for row in regime_drifts), default=0.0)
    return {
        "regime": regime,
        "case_count": len(regime_cases),
        "eta_values": sorted({float(case["eta"]) for case in regime_cases}),
        "all_cases_pass": bool(all(bool(case["pass"]) for case in regime_cases)),
        "drift_pass": bool(max_drift <= float(thresholds["max_refinement_profile_drift"])),
        "max_refinement_profile_drift": float(max_drift),
        "max_median_relative_error": float(max(case["median_relative_error"] for case in regime_cases)),
        "max_p90_relative_error": float(max(case["p90_relative_error"] for case in regime_cases)),
        "max_diffusivity_match_gap": float(max(case["diffusivity_match_gap"] for case in regime_cases)),
        "max_shell_kappa_cv": float(max(case["shell_kappa_cv"] for case in regime_cases)),
        "min_metric_tracking_abs_corr": float(min(case["metric_tracking_abs_corr"] for case in regime_cases)),
    }


def strongest_case(cases: list[dict[str, Any]], regime: str) -> dict[str, Any]:
    return max(
        [case for case in cases if case["regime"] == regime],
        key=lambda item: (
            item["median_relative_error"],
            item["p90_relative_error"],
            item["diffusivity_match_gap"],
            item["shell_kappa_cv"],
        ),
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
        for n_side in [int(value) for value in cfg["n_sides"]]
    }

    cases = []
    cases.extend(evaluate_regime("weak", [float(value) for value in cfg["weak_eta_values"]], baseline_cases, cfg))
    cases.extend(evaluate_regime("stress", [float(value) for value in cfg["stress_eta_values"]], baseline_cases, cfg))
    drifts = refinement_drifts(cases)
    summaries = [summarize_regime(regime, cases, drifts, thresholds) for regime in ("weak", "stress")]
    weak_summary = next(item for item in summaries if item["regime"] == "weak")
    stress_summary = next(item for item in summaries if item["regime"] == "stress")
    weak_edge = strongest_case(cases, "weak")
    stress_edge = strongest_case(cases, "stress")

    if weak_summary["all_cases_pass"] and weak_summary["drift_pass"]:
        observation = (
            "after excluding the earliest source-core transient layer, weak localized bump excitations close as a bounded "
            "metric-tracking transport response, while stronger localized bumps remain outside the same closure"
        )
        conclusion = (
            f"all {weak_summary['case_count']} weak localized bump cases pass for eta values {weak_summary['eta_values']}, "
            f"with maximum refinement drift {weak_summary['max_refinement_profile_drift']:.4f}; the weakest weak case is "
            f"`{weak_edge['label']}` with kappa/D_eff {weak_edge['conductivity_over_diffusivity']:.4f}, median relative error "
            f"{weak_edge['median_relative_error']:.4f}, p90 error {weak_edge['p90_relative_error']:.4f}, shell-kappa CV "
            f"{weak_edge['shell_kappa_cv']:.4f}, and |metric-tracking corr| {weak_edge['metric_tracking_abs_corr']:.4f}; "
            f"the stress window remains open with maximum drift {stress_summary['max_refinement_profile_drift']:.4f}, and the "
            f"hardest stress case is `{stress_edge['label']}` with median relative error {stress_edge['median_relative_error']:.4f}, "
            f"p90 error {stress_edge['p90_relative_error']:.4f}, and shell-kappa CV {stress_edge['shell_kappa_cv']:.4f}"
        )
    else:
        observation = (
            "the localized bump response improves after delaying the fit window, but a bounded weak localized-excitation "
            "closure is still not recovered on the tested scalar carrier"
        )
        conclusion = (
            f"weak localized bump cases passing = {sum(int(case['pass']) for case in cases if case['regime'] == 'weak')}/"
            f"{weak_summary['case_count']}, maximum weak drift {weak_summary['max_refinement_profile_drift']:.4f}; "
            f"the correct boundary remains open"
        )

    return {
        "experiment": "scalar_kernel_graph_localized_bump_response",
        "config": cfg,
        "baseline_cases": baseline_cases,
        "cases": cases,
        "refinement_drifts": drifts,
        "regime_summaries": summaries,
        "observation": observation,
        "conclusion": conclusion,
    }


def make_plots(result: dict[str, Any], stamp: str) -> list[str]:
    cases = result["cases"]
    thresholds = result["config"]["thresholds"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    colors = {"weak": "#0f766e", "stress": "#b45309"}
    markers = {13: "o", 15: "s"}
    for regime in ("weak", "stress"):
        for n_side in sorted({case["n_side"] for case in cases}):
            rows = sorted(
                [case for case in cases if case["regime"] == regime and case["n_side"] == n_side],
                key=lambda item: item["eta"],
            )
            label = f"{regime}|n={n_side}"
            axes[0, 0].plot([case["eta"] for case in rows], [case["median_relative_error"] for case in rows], marker=markers[n_side], color=colors[regime], label=label)
            axes[0, 1].plot([case["eta"] for case in rows], [case["p90_relative_error"] for case in rows], marker=markers[n_side], color=colors[regime], label=label)
            axes[1, 0].plot([case["eta"] for case in rows], [case["shell_kappa_cv"] for case in rows], marker=markers[n_side], color=colors[regime], label=label)
            axes[1, 1].plot([case["eta"] for case in rows], [case["conductivity_over_diffusivity"] for case in rows], marker=markers[n_side], color=colors[regime], label=label)

    axes[0, 0].axhline(thresholds["max_median_relative_error"], linestyle="--", color="#7f1d1d")
    axes[0, 1].axhline(thresholds["max_p90_relative_error"], linestyle="--", color="#7f1d1d")
    axes[1, 0].axhline(thresholds["max_shell_kappa_cv"], linestyle="--", color="#7f1d1d")
    axes[1, 1].axhline(1.0, linestyle="--", color="#111827")
    titles = ["Median relative error", "P90 relative error", "Shell-kappa CV", "Conductivity / diffusivity"]
    for axis, title in zip(axes.ravel(), titles):
        axis.set_title(title)
        axis.set_xlabel("localized bump amplitude eta")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_localized_bump_response_summary.png"
    fig.savefig(summary_path, dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for axis, regime in zip(axes, ("weak", "stress")):
        edge = strongest_case(cases, regime)
        radii = np.asarray(edge["shell_radii"], dtype=float)
        shell_profile = np.asarray(edge["shell_profile"], dtype=float)
        metric_shift = np.asarray(edge["normalized_metric_shift_profile"], dtype=float)
        kappa_profile = np.asarray(edge["normalized_shell_kappa_profile"], dtype=float)
        axis.plot(radii[: len(shell_profile)], normalized_profile(shell_profile), marker="o", label="bump profile")
        axis.plot(radii[: len(metric_shift)], metric_shift, marker="s", label="metric shift")
        axis.plot(radii[: len(kappa_profile)], kappa_profile, marker="^", label="shell kappa")
        axis.set_title(edge["label"])
        axis.set_xlabel("shell radius")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    profile_path = PLOTS / f"{stamp}_scalar_kernel_graph_localized_bump_response_profiles.png"
    fig.savefig(profile_path, dpi=170)
    plt.close(fig)

    return [str(summary_path.relative_to(REPO_ROOT)), str(profile_path.relative_to(REPO_ROOT))]


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_localized_bump_response.json"
    latest = DATA / "scalar_kernel_graph_localized_bump_response_latest.json"
    stamped.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    latest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, latest_path: str, plot_paths: list[str]) -> None:
    summary_rows = "\n".join(
        f"| `{row['regime']}` | `{row['eta_values']}` | `{row['all_cases_pass']}` | `{row['max_refinement_profile_drift']:.4f}` | `{row['max_median_relative_error']:.4f}` | `{row['max_p90_relative_error']:.4f}` | `{row['max_shell_kappa_cv']:.4f}` |"
        for row in result["regime_summaries"]
    )
    case_rows = "\n".join(
        f"| `{case['label']}` | `{case['conductivity_over_diffusivity']:.4f}` | `{case['median_relative_error']:.4f}` | `{case['p90_relative_error']:.4f}` | `{case['shell_kappa_cv']:.4f}` | `{case['metric_tracking_abs_corr']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in sorted(result["cases"], key=lambda item: (item["regime"], item["eta"], item["n_side"]))
    )
    drift_rows = "\n".join(
        f"| `{row['label']}` | `{row['drift']:.4f}` |" for row in result["refinement_drifts"]
    )
    note = f"""# Scalar Kernel Graph Localized Bump Response v1

## Purpose

Resolve the broad `52.5` localized-bump OPEN row into a sharper threshold statement.

The original bump test used the same early transient window as the smooth radial branch. That made the localized source-core layer dominate the shellwise constitutive fit and produced a maximum refinement drift of `0.5124`.

This follow-up keeps the same scalar carrier and the same bump shape, but excludes the earliest source-core transient layer with fit window `{result['config']['fit_window']}`.

## Construction

- bump sigma: `{result['config']['bump_sigma']}`
- refinements: `{result['config']['n_sides']}`
- weak eta values: `{result['config']['weak_eta_values']}`
- stress eta values: `{result['config']['stress_eta_values']}`
- thresholds: `{result['config']['thresholds']}`

## Regime Summary

| regime | eta values | all cases pass | max drift | max median error | max p90 error | max shell-kappa CV |
| --- | --- | --- | --- | --- | --- | --- |
{summary_rows}

## Case Table

| case | `kappa / D_eff` | median error | p90 error | shell-kappa CV | metric corr | status |
| --- | --- | --- | --- | --- | --- | --- |
{case_rows}

## Refinement Drift

| case | drift |
| --- | --- |
{drift_rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Authority

- script: `numerics/simulations/scalar_kernel_graph_localized_bump_response.py`
- result: `{result_path}`
- latest: `{latest_path}`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph localized bump response\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, weak_eta_values={config['weak_eta_values']}, "
            f"stress_eta_values={config['stress_eta_values']}, fit_window={config['fit_window']}, "
            f"bump_sigma={config['bump_sigma']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = make_result()
    plot_paths = make_plots(result, stamp)
    result_path, latest_path = save_results(result, stamp)
    write_note(result, result_path, latest_path, plot_paths)
    append_log(result_path, plot_paths, result["config"], result["observation"], result["conclusion"])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
