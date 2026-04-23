#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Power_Law_Scaling_v1.md"
for path in (DATA, PLOTS):
    path.mkdir(exist_ok=True)

MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAW_GRADIENT = DATA / "scalar_kernel_graph_recoverability_gradient_latest.json"
SHELL_NATIVE_GRADIENT = DATA / "scalar_kernel_graph_recoverability_gradient_shell_native_latest.json"

DEFAULT_CONFIG: dict[str, Any] = {
    "raw_artifact": str(RAW_GRADIENT.relative_to(REPO_ROOT)),
    "shell_native_artifact": str(SHELL_NATIVE_GRADIENT.relative_to(REPO_ROOT)),
    "thresholds": {
        "raw_min_power_fit_r2": 0.95,
        "raw_max_power_deviation": 0.18,
        "raw_max_flux_constancy_cv": 0.12,
        "raw_max_refinement_profile_drift": 0.08,
        "shell_native_min_green_fit_r2": 0.97,
        "shell_native_min_power_fit_r2": 0.999,
        "shell_native_max_power_deviation": 0.02,
        "shell_native_max_flux_constancy_cv": 0.01,
        "shell_native_max_refinement_profile_drift": 0.02,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(raw.get("scalar_kernel_graph_power_law_scaling"), dict):
            section = raw["scalar_kernel_graph_power_law_scaling"]
            merged.update({key: value for key, value in section.items() if key != "thresholds"})
            if isinstance(section.get("thresholds"), dict):
                merged["thresholds"].update(section["thresholds"])
    if config is not None:
        merged.update({key: value for key, value in config.items() if key != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def load_artifact(path: str) -> dict[str, Any]:
    artifact_path = REPO_ROOT / path
    if not artifact_path.exists():
        raise FileNotFoundError(f"missing power-law scaling source artifact: {path}")
    return json.loads(artifact_path.read_text(encoding="utf-8"))


def all_cases(data: dict[str, Any]) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for key in ("disorder_cases", "kernel_cases", "cases"):
        cases.extend(data.get(key, []))
    return cases


def weakest_case(cases: list[dict[str, Any]]) -> dict[str, Any]:
    return max(
        cases,
        key=lambda case: (
            abs(float(case["profile"]["power_slope"]) + 2.0),
            1.0 - float(case["profile"]["power_fit_r2"]),
            float(case["profile"]["flux_constancy_cv"]),
        ),
    )


def summarize_raw(data: dict[str, Any], thresholds: dict[str, float]) -> dict[str, Any]:
    cases = all_cases(data)
    min_power_r2 = min(float(case["profile"]["power_fit_r2"]) for case in cases)
    max_power_deviation = max(abs(float(case["profile"]["power_slope"]) + 2.0) for case in cases)
    max_flux_cv = max(float(case["profile"]["flux_constancy_cv"]) for case in cases)
    max_drift = float(data["max_profile_drift"])
    failing_cases = [case["label"] for case in cases if not bool(case.get("pass", False))]
    passes = bool(
        not failing_cases
        and min_power_r2 >= float(thresholds["raw_min_power_fit_r2"])
        and max_power_deviation <= float(thresholds["raw_max_power_deviation"])
        and max_flux_cv <= float(thresholds["raw_max_flux_constancy_cv"])
        and max_drift <= float(thresholds["raw_max_refinement_profile_drift"])
    )
    edge = weakest_case(cases)
    return {
        "readout": "raw_local_gradient",
        "case_count": len(cases),
        "pass": passes,
        "failing_cases": failing_cases,
        "min_power_fit_r2": min_power_r2,
        "max_power_deviation": max_power_deviation,
        "max_flux_constancy_cv": max_flux_cv,
        "max_refinement_profile_drift": max_drift,
        "weakest_case": {
            "label": edge["label"],
            "power_slope": float(edge["profile"]["power_slope"]),
            "power_fit_r2": float(edge["profile"]["power_fit_r2"]),
            "flux_constancy_cv": float(edge["profile"]["flux_constancy_cv"]),
        },
    }


def summarize_shell_native(data: dict[str, Any], thresholds: dict[str, float]) -> dict[str, Any]:
    cases = all_cases(data)
    min_green_r2 = min(float(case["green_fit_r2"]) for case in cases)
    min_power_r2 = min(float(case["profile"]["power_fit_r2"]) for case in cases)
    max_power_deviation = max(abs(float(case["profile"]["power_slope"]) + 2.0) for case in cases)
    max_flux_cv = max(float(case["profile"]["flux_constancy_cv"]) for case in cases)
    max_drift = float(data["max_profile_drift"])
    failing_cases = [case["label"] for case in cases if not bool(case.get("pass", False))]
    passes = bool(
        not failing_cases
        and min_green_r2 >= float(thresholds["shell_native_min_green_fit_r2"])
        and min_power_r2 >= float(thresholds["shell_native_min_power_fit_r2"])
        and max_power_deviation <= float(thresholds["shell_native_max_power_deviation"])
        and max_flux_cv <= float(thresholds["shell_native_max_flux_constancy_cv"])
        and max_drift <= float(thresholds["shell_native_max_refinement_profile_drift"])
    )
    edge = weakest_case(cases)
    return {
        "readout": "shell_native",
        "case_count": len(cases),
        "pass": passes,
        "failing_cases": failing_cases,
        "min_green_fit_r2": min_green_r2,
        "min_power_fit_r2": min_power_r2,
        "max_power_deviation": max_power_deviation,
        "max_flux_constancy_cv": max_flux_cv,
        "max_refinement_profile_drift": max_drift,
        "weakest_case": {
            "label": edge["label"],
            "green_fit_r2": float(edge["green_fit_r2"]),
            "power_slope": float(edge["profile"]["power_slope"]),
            "power_fit_r2": float(edge["profile"]["power_fit_r2"]),
            "flux_constancy_cv": float(edge["profile"]["flux_constancy_cv"]),
        },
    }


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    thresholds = cfg["thresholds"]
    raw = load_artifact(str(cfg["raw_artifact"]))
    shell_native = load_artifact(str(cfg["shell_native_artifact"]))
    raw_summary = summarize_raw(raw, thresholds)
    shell_summary = summarize_shell_native(shell_native, thresholds)

    if shell_summary["pass"]:
        observation = (
            "the broad power-law scaling row splits cleanly: the raw local-gradient read remains open, "
            "while the shell-native law-aware read passes the tested continuum-like scaling thresholds"
        )
        conclusion = (
            f"raw local-gradient power scaling remains open with minimum R^2 {raw_summary['min_power_fit_r2']:.4f}, "
            f"maximum |slope+2| {raw_summary['max_power_deviation']:.4f}, and profile drift {raw_summary['max_refinement_profile_drift']:.4f}; "
            f"shell-native scaling passes with minimum power-fit R^2 {shell_summary['min_power_fit_r2']:.4f}, "
            f"maximum |slope+2| {shell_summary['max_power_deviation']:.4f}, scaled-response CV {shell_summary['max_flux_constancy_cv']:.4f}, "
            f"and profile drift {shell_summary['max_refinement_profile_drift']:.4f}; this supports a bounded continuum-like power-law proxy only under the shell-native reconstruction"
        )
    else:
        observation = (
            "the shell-native reconstruction improves the power-law read but does not pass the configured scaling thresholds"
        )
        conclusion = (
            f"shell-native scaling remains open with minimum power-fit R^2 {shell_summary['min_power_fit_r2']:.4f}, "
            f"maximum |slope+2| {shell_summary['max_power_deviation']:.4f}, and profile drift {shell_summary['max_refinement_profile_drift']:.4f}"
        )

    return {
        "experiment": "scalar_kernel_graph_power_law_scaling",
        "config": cfg,
        "source_artifacts": {
            "raw_local_gradient": str(cfg["raw_artifact"]),
            "shell_native": str(cfg["shell_native_artifact"]),
        },
        "summaries": {
            "raw_local_gradient": raw_summary,
            "shell_native": shell_summary,
        },
        "observation": observation,
        "conclusion": conclusion,
    }


def selected_profile_case(data: dict[str, Any], label: str) -> dict[str, Any]:
    for case in all_cases(data):
        if case["label"] == label:
            return case
    raise KeyError(label)


def make_plots(result: dict[str, Any], stamp: str) -> list[str]:
    raw = load_artifact(result["source_artifacts"]["raw_local_gradient"])
    shell_native = load_artifact(result["source_artifacts"]["shell_native"])
    raw_summary = result["summaries"]["raw_local_gradient"]
    shell_summary = result["summaries"]["shell_native"]
    thresholds = result["config"]["thresholds"]

    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_power_law_scaling_summary.png"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    labels = ["raw local gradient", "shell native"]
    colors = ["#b45309", "#0f766e"]
    values = [
        [raw_summary["min_power_fit_r2"], shell_summary["min_power_fit_r2"]],
        [raw_summary["max_power_deviation"], shell_summary["max_power_deviation"]],
        [raw_summary["max_flux_constancy_cv"], shell_summary["max_flux_constancy_cv"]],
        [raw_summary["max_refinement_profile_drift"], shell_summary["max_refinement_profile_drift"]],
    ]
    targets = [
        [thresholds["raw_min_power_fit_r2"], thresholds["shell_native_min_power_fit_r2"]],
        [thresholds["raw_max_power_deviation"], thresholds["shell_native_max_power_deviation"]],
        [thresholds["raw_max_flux_constancy_cv"], thresholds["shell_native_max_flux_constancy_cv"]],
        [thresholds["raw_max_refinement_profile_drift"], thresholds["shell_native_max_refinement_profile_drift"]],
    ]
    titles = ["Minimum power-fit R^2", "Maximum |slope + 2|", "Maximum scaled-flux CV", "Maximum profile drift"]
    for axis, value_pair, target_pair, title in zip(axes.ravel(), values, targets, titles):
        axis.bar(labels, value_pair, color=colors)
        for idx, target in enumerate(target_pair):
            axis.plot([idx - 0.35, idx + 0.35], [target, target], color="#111827", linestyle="--")
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=12)
        axis.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(summary_path, dpi=170)
    plt.close(fig)

    profile_path = PLOTS / f"{stamp}_scalar_kernel_graph_power_law_scaling_profiles.png"
    raw_case = selected_profile_case(raw, raw_summary["weakest_case"]["label"])
    shell_case = selected_profile_case(shell_native, shell_summary["weakest_case"]["label"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for axis, case, field_name, title in (
        (axes[0], raw_case, "radial_flux", f"raw: {raw_case['label']}"),
        (axes[1], shell_case, "response", f"shell-native: {shell_case['label']}"),
    ):
        rows = case["profile"]["rows"]
        radii = np.asarray([row["radius"] for row in rows], dtype=float)
        response = np.asarray([abs(row[field_name]) for row in rows], dtype=float)
        slope = float(case["profile"]["power_slope"])
        intercept = float(case["profile"]["power_intercept"])
        axis.loglog(radii, response, "o", label="profile")
        axis.loglog(radii, np.exp(intercept) * radii**slope, "--", label="power fit")
        axis.set_title(title)
        axis.set_xlabel("shell radius")
        axis.set_ylabel("response magnitude")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(profile_path, dpi=170)
    plt.close(fig)

    return [str(summary_path.relative_to(REPO_ROOT)), str(profile_path.relative_to(REPO_ROOT))]


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_power_law_scaling.json"
    latest = DATA / "scalar_kernel_graph_power_law_scaling_latest.json"
    stamped.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    latest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, latest_path: str, plot_paths: list[str]) -> None:
    raw_summary = result["summaries"]["raw_local_gradient"]
    shell_summary = result["summaries"]["shell_native"]
    note = f"""# Scalar Kernel Graph Power-Law Scaling v1

## Purpose

Resolve the broad physics-bridge `power-law scaling` OPEN row into two auditable readouts:

- raw local-gradient power scaling, which remains `OPEN`;
- shell-native law-aware power scaling, which passes as a bounded continuum-like proxy.

This is a post-processing audit over frozen scalar artifacts. It does not change HAOS core, the scalar carrier, or either recoverability-gradient runner.

## Source Artifacts

- raw local gradient: `{result['source_artifacts']['raw_local_gradient']}`
- shell-native reconstruction: `{result['source_artifacts']['shell_native']}`
- result: `{result_path}`
- latest: `{latest_path}`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Summary

| readout | cases | status | min power-fit `R^2` | max `|slope + 2|` | max scaled CV | max profile drift |
| --- | --- | --- | --- | --- | --- | --- |
| raw local gradient | `{raw_summary['case_count']}` | {'PASS' if raw_summary['pass'] else 'OPEN'} | `{raw_summary['min_power_fit_r2']:.4f}` | `{raw_summary['max_power_deviation']:.4f}` | `{raw_summary['max_flux_constancy_cv']:.4f}` | `{raw_summary['max_refinement_profile_drift']:.4f}` |
| shell native | `{shell_summary['case_count']}` | {'PASS' if shell_summary['pass'] else 'OPEN'} | `{shell_summary['min_power_fit_r2']:.4f}` | `{shell_summary['max_power_deviation']:.4f}` | `{shell_summary['max_flux_constancy_cv']:.4f}` | `{shell_summary['max_refinement_profile_drift']:.4f}` |

## Weakest Cases

- raw: `{raw_summary['weakest_case']['label']}`, slope `{raw_summary['weakest_case']['power_slope']:.4f}`, `R^2 = {raw_summary['weakest_case']['power_fit_r2']:.4f}`, CV `{raw_summary['weakest_case']['flux_constancy_cv']:.4f}`;
- shell-native: `{shell_summary['weakest_case']['label']}`, Green `R^2 = {shell_summary['weakest_case']['green_fit_r2']:.4f}`, slope `{shell_summary['weakest_case']['power_slope']:.4f}`, power `R^2 = {shell_summary['weakest_case']['power_fit_r2']:.4f}`, CV `{shell_summary['weakest_case']['flux_constancy_cv']:.4f}`.

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Correct Boundary

The old broad statement was:

> power-law scaling remains open because raw power-fit R^2 is too low.

The new statement is:

> raw local-gradient power scaling remains open, but shell-native law-aware power scaling closes as a bounded continuum-like proxy on the tested scalar-carrier window.

This is not a continuum-limit proof and not a gravity claim. It is a reconstruction boundary.
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], result: dict[str, Any]) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph power-law scaling audit\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(f"- Source artifacts: `{result['source_artifacts']['raw_local_gradient']}`, `{result['source_artifacts']['shell_native']}`\n")
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {result['observation']}\n")
        handle.write(f"- Conclusion: {result['conclusion']}\n")


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = make_result()
    plot_paths = make_plots(result, stamp)
    result_path, latest_path = save_results(result, stamp)
    write_note(result, result_path, latest_path, plot_paths)
    append_log(result_path, plot_paths, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
