#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from scalar_kernel_graph_geometry_robustness import (
    build_cubic_points,
    build_jittered_points,
    build_point_cloud_operator,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Metric_Field_v1.md"
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
    "coarse_radius_factor": 2.5,
    "thresholds": {
        "min_bulk_count": 100,
        "min_min_eigenvalue": 0.80,
        "max_mean_anisotropy": 0.02,
        "max_mean_offdiag_ratio": 0.02,
        "max_spatial_trace_cv": 0.05,
        "max_normalized_mean_tensor_error": 0.02,
        "max_refinement_drift": 0.03,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_metric_field"), dict):
            merged.update({k: v for k, v in raw["scalar_kernel_graph_metric_field"].items() if k != "thresholds"})
            if isinstance(raw["scalar_kernel_graph_metric_field"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_metric_field"]["thresholds"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def local_metric_field(points: np.ndarray, delta: Any) -> np.ndarray:
    delta_csr = delta.tocsr()
    tensors = np.zeros((points.shape[0], 3, 3), dtype=float)
    indptr = delta_csr.indptr
    indices = delta_csr.indices
    data = delta_csr.data
    for node in range(points.shape[0]):
        start, end = indptr[node], indptr[node + 1]
        neighbors = indices[start:end]
        weights = data[start:end]
        mask = neighbors != node
        if not np.any(mask):
            continue
        deltas = points[neighbors[mask]] - points[node]
        row_weights = weights[mask]
        tensor = 0.5 * (deltas.T @ (row_weights[:, None] * deltas))
        tensors[node] = 0.5 * (tensor + tensor.T)
    return tensors


def bulk_mask(points: np.ndarray, h: float, bulk_margin_layers: int) -> np.ndarray:
    margin = float(bulk_margin_layers) * h - 1.0e-12
    return np.all((points >= margin) & (points <= 1.0 - margin), axis=1)


def normalize_tensor(tensor: np.ndarray) -> tuple[np.ndarray, float]:
    scale = float(np.trace(tensor) / 3.0)
    return tensor / max(scale, 1.0e-18), scale


def summarize_metric_samples(field: np.ndarray) -> dict[str, Any]:
    eigenvalues = np.linalg.eigvalsh(field)
    diagonals = np.diagonal(field, axis1=1, axis2=2)
    traces = np.sum(diagonals, axis=1)
    mean_scales = traces / 3.0
    mean_tensor = np.mean(field, axis=0)
    normalized_mean, mean_scale = normalize_tensor(mean_tensor)
    identity = np.eye(3, dtype=float)
    offdiag = field - np.einsum("nii,ij->nij", field, identity)
    node_anisotropy = (eigenvalues[:, -1] - eigenvalues[:, 0]) / np.maximum(mean_scales, 1.0e-18)
    node_offdiag_ratio = np.linalg.norm(offdiag, axis=(1, 2)) / np.maximum(mean_scales, 1.0e-18)
    return {
        "mean_tensor": mean_tensor.tolist(),
        "normalized_mean_tensor": normalized_mean.tolist(),
        "mean_scale": float(mean_scale),
        "mean_anisotropy": float(np.mean(node_anisotropy)),
        "max_anisotropy": float(np.max(node_anisotropy)),
        "mean_offdiag_ratio": float(np.mean(node_offdiag_ratio)),
        "max_offdiag_ratio": float(np.max(node_offdiag_ratio)),
        "min_eigenvalue": float(np.min(eigenvalues)),
        "max_eigenvalue": float(np.max(eigenvalues)),
        "mean_eigenvalues": np.mean(eigenvalues, axis=0).tolist(),
        "spatial_trace_cv": float(np.std(traces) / max(np.mean(traces), 1.0e-18)),
        "normalized_mean_tensor_error": float(np.linalg.norm(normalized_mean - identity)),
    }


def coarse_grain_metric_field(points: np.ndarray, field: np.ndarray, radius: float) -> np.ndarray:
    tree = cKDTree(points)
    smoothed = np.empty_like(field)
    for idx, point in enumerate(points):
        neighborhood = tree.query_ball_point(point, radius)
        smoothed[idx] = np.mean(field[neighborhood], axis=0)
    return smoothed


def metric_field_metrics(points: np.ndarray, h: float, delta: Any, bulk_margin_layers: int, coarse_radius_factor: float) -> dict[str, Any]:
    field = local_metric_field(points=points, delta=delta)
    mask = bulk_mask(points=points, h=h, bulk_margin_layers=bulk_margin_layers)
    bulk_points = points[mask]
    raw_bulk = field[mask]
    coarse_bulk = coarse_grain_metric_field(points=bulk_points, field=raw_bulk, radius=float(coarse_radius_factor) * h)
    raw_summary = summarize_metric_samples(raw_bulk)
    coarse_summary = summarize_metric_samples(coarse_bulk)
    return {
        "field": field.tolist(),
        "bulk_count": int(np.sum(mask)),
        "bulk_fraction": float(np.mean(mask)),
        "coarse_radius_factor": float(coarse_radius_factor),
        "raw": raw_summary,
        "coarse": coarse_summary,
    }


def evaluate_case(
    label: str,
    points: np.ndarray,
    h: float,
    epsilon_coeff: float,
    cutoff_factor: float,
    kernel_family: str,
    bulk_margin_layers: int,
    coarse_radius_factor: float,
) -> dict[str, Any]:
    delta, epsilon_k = build_point_cloud_operator(
        points=points,
        h=h,
        epsilon_coeff=epsilon_coeff,
        cutoff_factor=cutoff_factor,
        kernel_family=kernel_family,
    )
    metrics = metric_field_metrics(
        points=points,
        h=h,
        delta=delta,
        bulk_margin_layers=bulk_margin_layers,
        coarse_radius_factor=coarse_radius_factor,
    )
    return {
        "label": label,
        "points_count": int(points.shape[0]),
        "lattice_spacing": float(h),
        "epsilon_coeff": float(epsilon_coeff),
        "epsilon_k": float(epsilon_k),
        "kernel_family": kernel_family,
        "metric_field": metrics,
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    metrics = case["metric_field"]["coarse"]
    return bool(
        case["metric_field"]["bulk_count"] >= int(thresholds["min_bulk_count"])
        and metrics["min_eigenvalue"] >= float(thresholds["min_min_eigenvalue"])
        and metrics["mean_anisotropy"] <= float(thresholds["max_mean_anisotropy"])
        and metrics["mean_offdiag_ratio"] <= float(thresholds["max_mean_offdiag_ratio"])
        and metrics["spatial_trace_cv"] <= float(thresholds["max_spatial_trace_cv"])
        and metrics["normalized_mean_tensor_error"] <= float(thresholds["max_normalized_mean_tensor_error"])
    )


def refinement_drift(case_a: dict[str, Any], case_b: dict[str, Any]) -> float:
    tensor_a = np.asarray(case_a["metric_field"]["coarse"]["normalized_mean_tensor"], dtype=float)
    tensor_b = np.asarray(case_b["metric_field"]["coarse"]["normalized_mean_tensor"], dtype=float)
    return float(np.linalg.norm(tensor_a - tensor_b))


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    thresholds = cfg["thresholds"]

    disorder_cases: list[dict[str, Any]] = []
    for n_side in [int(v) for v in cfg["disorder_n_sides"]]:
        for jitter_fraction in [float(v) for v in cfg["jitter_fractions"]]:
            points, _, h = build_jittered_points(
                n_side=n_side,
                jitter_fraction=jitter_fraction,
                seed=int(cfg["jitter_seed"]),
            )
            case = evaluate_case(
                label=f"disorder|n={n_side}|jitter={jitter_fraction:.3f}",
                points=points,
                h=h,
                epsilon_coeff=float(cfg["epsilon_coeff"]),
                cutoff_factor=float(cfg["cutoff_factor"]),
                kernel_family="gaussian_local",
                bulk_margin_layers=int(cfg["bulk_margin_layers"]),
                coarse_radius_factor=float(cfg["coarse_radius_factor"]),
            )
            case["jitter_fraction"] = float(jitter_fraction)
            case["pass"] = case_passes(case, thresholds)
            disorder_cases.append(case)

    kernel_cases: list[dict[str, Any]] = []
    for n_side in [int(v) for v in cfg["kernel_n_sides"]]:
        points, _, h = build_cubic_points(n_side)
        for family in [str(v) for v in cfg["kernel_families"]]:
            case = evaluate_case(
                label=f"kernel|n={n_side}|family={family}",
                points=points,
                h=h,
                epsilon_coeff=float(cfg["epsilon_coeff"]),
                cutoff_factor=float(cfg["cutoff_factor"]),
                kernel_family=family,
                bulk_margin_layers=int(cfg["bulk_margin_layers"]),
                coarse_radius_factor=float(cfg["coarse_radius_factor"]),
            )
            case["pass"] = case_passes(case, thresholds)
            kernel_cases.append(case)

    disorder_drifts = []
    for jitter_fraction in [float(v) for v in cfg["jitter_fractions"]]:
        pair = [case for case in disorder_cases if abs(case["jitter_fraction"] - jitter_fraction) < 1.0e-12]
        if len(pair) == 2:
            pair = sorted(pair, key=lambda item: item["points_count"])
            disorder_drifts.append(
                {
                    "label": f"disorder|jitter={jitter_fraction:.3f}",
                    "drift": refinement_drift(pair[0], pair[1]),
                }
            )

    kernel_drifts = []
    for family in [str(v) for v in cfg["kernel_families"]]:
        pair = [case for case in kernel_cases if case["kernel_family"] == family]
        if len(pair) == 2:
            pair = sorted(pair, key=lambda item: item["points_count"])
            kernel_drifts.append(
                {
                    "label": f"kernel|family={family}",
                    "drift": refinement_drift(pair[0], pair[1]),
                }
            )

    max_refinement_drift = max([item["drift"] for item in disorder_drifts + kernel_drifts], default=0.0)
    all_cases_pass = all(case["pass"] for case in disorder_cases + kernel_cases)
    drift_pass = max_refinement_drift <= float(thresholds["max_refinement_drift"])

    weakest_case = max(
        disorder_cases + kernel_cases,
        key=lambda item: (
            item["metric_field"]["coarse"]["normalized_mean_tensor_error"],
            item["metric_field"]["coarse"]["mean_anisotropy"],
            item["metric_field"]["coarse"]["spatial_trace_cv"],
        ),
    )

    if all_cases_pass and drift_pass:
        observation = (
            "the same scalar carrier now supports a stable local metric-like tensor field: the operator-level row-local quadratic response stays positive, near-isotropic, "
            "and nearly constant across the tested bulk window under refinement, mild disorder, and bounded kernel-family variation"
        )
        conclusion = (
            f"all {len(disorder_cases)} disorder cases and all {len(kernel_cases)} kernel-family cases pass the local metric-field thresholds, and the maximum normalized refinement drift is {max_refinement_drift:.4f}; "
            f"the weakest passing case is `{weakest_case['label']}` with normalized mean-tensor error {weakest_case['metric_field']['coarse']['normalized_mean_tensor_error']:.4f}, "
            f"mean anisotropy {weakest_case['metric_field']['coarse']['mean_anisotropy']:.4f}, and spatial trace CV {weakest_case['metric_field']['coarse']['spatial_trace_cv']:.4f}; "
            "this means the scalar-carrier geometry bridge now extends from one global coordinate-stiffness read to a bounded stable local metric-field statement on the tested window"
        )
    else:
        failing_cases = [case["label"] for case in disorder_cases + kernel_cases if not case["pass"]]
        observation = (
            "the scalar carrier supports a positive global geometry read, but the local metric-field extraction still exposes instability or drift beyond the bounded thresholds"
        )
        conclusion = (
            f"open local metric-field cases remain: failing_cases={failing_cases}, max_refinement_drift={max_refinement_drift:.4f}; "
            "the correct boundary therefore stays at the global scalar-carrier geometry closure and robustness statement without yet promoting a positive local metric-field claim"
        )

    return {
        "experiment": "scalar_kernel_graph_metric_field",
        "config": cfg,
        "disorder_cases": disorder_cases,
        "kernel_cases": kernel_cases,
        "disorder_refinement_drifts": disorder_drifts,
        "kernel_refinement_drifts": kernel_drifts,
        "max_refinement_drift": float(max_refinement_drift),
        "observation": observation,
        "conclusion": conclusion,
    }


def make_summary_plot(disorder_cases: list[dict[str, Any]], kernel_cases: list[dict[str, Any]], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_metric_field_summary.png"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for n_side in sorted({int(case["label"].split("|")[1].split("=")[1]) for case in disorder_cases}):
        subset = sorted(
            [case for case in disorder_cases if int(case["label"].split("|")[1].split("=")[1]) == n_side],
            key=lambda item: float(item["jitter_fraction"]),
        )
        x = [float(case["jitter_fraction"]) for case in subset]
        axes[0, 0].plot(x, [case["metric_field"]["coarse"]["mean_anisotropy"] for case in subset], marker="o", label=f"n={n_side}")
        axes[0, 1].plot(x, [case["metric_field"]["coarse"]["normalized_mean_tensor_error"] for case in subset], marker="o", label=f"n={n_side}")
        axes[1, 0].plot(x, [case["metric_field"]["coarse"]["spatial_trace_cv"] for case in subset], marker="o", label=f"n={n_side}")

    family_order = [str(v) for v in DEFAULT_CONFIG["kernel_families"]]
    for n_side in sorted({int(case["label"].split("|")[1].split("=")[1]) for case in kernel_cases}):
        subset = sorted(
            [case for case in kernel_cases if int(case["label"].split("|")[1].split("=")[1]) == n_side],
            key=lambda item: family_order.index(str(item["kernel_family"])),
        )
        x = range(len(subset))
        axes[1, 1].plot(x, [case["metric_field"]["coarse"]["mean_offdiag_ratio"] for case in subset], marker="o", label=f"n={n_side}")
        axes[1, 1].set_xticks(list(x), labels=[case["kernel_family"] for case in subset], rotation=20)

    axes[0, 0].set_title("Local mean anisotropy")
    axes[0, 0].set_xlabel("jitter fraction")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].set_title("Normalized mean-tensor error")
    axes[0, 1].set_xlabel("jitter fraction")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].set_title("Spatial trace coefficient of variation")
    axes[1, 0].set_xlabel("jitter fraction")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].set_title("Local mean offdiag ratio")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend()

    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def make_tensor_plot(result: dict[str, Any], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_metric_field_tensors.png"
    selected = [
        [case for case in result["disorder_cases"] if abs(case["jitter_fraction"] - 0.04) < 1.0e-12][-1],
        [case for case in result["kernel_cases"] if case["kernel_family"] == "gaussian_local"][-1],
        [case for case in result["kernel_cases"] if case["kernel_family"] == "gaussian_half"][-1],
        [case for case in result["kernel_cases"] if case["kernel_family"] == "inverse_quadratic"][-1],
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    for axis, case in zip(axes.ravel(), selected):
        tensor = np.asarray(case["metric_field"]["coarse"]["normalized_mean_tensor"], dtype=float)
        im = axis.imshow(tensor, cmap="viridis", vmin=0.0, vmax=1.05)
        axis.set_title(case["label"].replace("|", "\n"))
        axis.set_xticks(range(3), labels=["x", "y", "z"])
        axis.set_yticks(range(3), labels=["x", "y", "z"])
        fig.colorbar(im, ax=axis, fraction=0.046, pad=0.04)
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_metric_field.json"
    latest = DATA / "scalar_kernel_graph_metric_field_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    disorder_rows = "\n".join(
        f"| `{case['label']}` | `{case['metric_field']['bulk_count']}` | `{case['metric_field']['coarse']['min_eigenvalue']:.4f}` | `{case['metric_field']['coarse']['mean_anisotropy']:.4f}` | `{case['metric_field']['coarse']['mean_offdiag_ratio']:.4f}` | `{case['metric_field']['coarse']['spatial_trace_cv']:.4f}` | `{case['metric_field']['coarse']['normalized_mean_tensor_error']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["disorder_cases"]
    )
    kernel_rows = "\n".join(
        f"| `{case['label']}` | `{case['metric_field']['bulk_count']}` | `{case['metric_field']['coarse']['min_eigenvalue']:.4f}` | `{case['metric_field']['coarse']['mean_anisotropy']:.4f}` | `{case['metric_field']['coarse']['mean_offdiag_ratio']:.4f}` | `{case['metric_field']['coarse']['spatial_trace_cv']:.4f}` | `{case['metric_field']['coarse']['normalized_mean_tensor_error']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["kernel_cases"]
    )
    drift_rows = "\n".join(
        f"| `{row['label']}` | `{row['drift']:.4f}` |"
        for row in result["disorder_refinement_drifts"] + result["kernel_refinement_drifts"]
    )
    note = f"""# Scalar Kernel Graph Metric Field v1

## Purpose

Upgrade the scalar-carrier geometry read from one global coordinate-stiffness tensor to a **local metric-like tensor field** on the same operator/substrate family.

The bounded question is:

> does the same scalar kernel-graph carrier support a positive, near-isotropic, spatially stable local metric-like field in the bulk, under the same mild disorder and bounded kernel-family window already used in the `51.3` robustness pass?

## Construction

The carrier stays fixed:

- scalar kernel graph on `3D` cubic point clouds
- local regime `c_epsilon = 0.5`
- same mild point-set disorder window
- same bounded kernel-family window

The local metric-field estimator is taken directly from the operator's row-local quadratic response. For each node `i`,

$$
G_{{ab}}(i) = \\frac{{1}}{{2}} \\sum_{{j \\neq i}} \\widehat{{L}}_{{ij}} (x_j^a - x_i^a)(x_j^b - x_i^b).
$$

This is the right bounded local estimator on the present carrier because it measures the operator's second-order response directly at each node, rather than importing a new geometric fitting rule.

The raw node-level field is kept as provenance, but the bounded physical read is a local coarse-grained field averaged over a bulk neighborhood radius of `{result['config']['coarse_radius_factor']} h`. This is the same kind of move the repository already used in other places when raw local readouts were too stencil-sensitive to carry the physical statement alone.

Bulk readout excludes a boundary collar of `{result['config']['bulk_margin_layers']}` lattice layers.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_metric_field.py`
- config: `config.json -> scalar_kernel_graph_metric_field`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_metric_field_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Disorder pass

| case | bulk nodes | min eigenvalue | mean anisotropy | mean offdiag ratio | spatial trace CV | normalized mean-tensor error | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
{disorder_rows}

## Kernel-family pass

| case | bulk nodes | min eigenvalue | mean anisotropy | mean offdiag ratio | spatial trace CV | normalized mean-tensor error | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
{kernel_rows}

## Refinement drift of the normalized mean tensor

| pair | drift |
| --- | --- |
{drift_rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Current boundary

If this note is read positively, the correct bounded statement is:

> the scalar-carrier geometry bridge now extends to a stable local bulk metric-like tensor field on the tested window.

This note still does **not** claim:

- curvature extraction
- current conservation
- broad universality beyond the tested mild carrier window
- ontology or spacetime
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## scalar kernel graph metric field\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"disorder_n_sides={config['disorder_n_sides']}, "
            f"jitter_fractions={config['jitter_fractions']}, "
            f"kernel_n_sides={config['kernel_n_sides']}, "
            f"kernel_families={config['kernel_families']}, "
            f"bulk_margin_layers={config['bulk_margin_layers']}, "
            f"coarse_radius_factor={config['coarse_radius_factor']}\n"
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
        make_tensor_plot(result, stamp),
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
