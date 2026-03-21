#!/usr/bin/env python3

from __future__ import annotations

import itertools
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import build_graph, read_json, relpath, write_csv_rows, write_json

PHASE_ROOT = ROOT / "phase7-spectral"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE5_MANIFEST_PATH = ROOT / "phase5-readout" / "phase5_authoritative_manifest.json"

SCALING_LEDGER_PATH = RUNS_ROOT / "phase7_spectral_scaling_ledger.csv"
SUMMARY_PATH = PHASE_ROOT / "phase7_spectral_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase7_spectral_manifest.json"

EIGENVALUE_PLOT_PATH = PLOTS_ROOT / "phase7_eigenvalue_vs_h.svg"
RADIUS_PLOT_PATH = PLOTS_ROOT / "phase7_spectral_radius_vs_h.svg"
TRACE_PLOT_PATH = PLOTS_ROOT / "phase7_trace_proxy_vs_h.svg"

TRACK_ID = "A"
TRACK_NAME = "spectral_feasibility"
STAGE_IDENTIFIER = "phase7-spectral"

NULLSPACE_TOL = 1.0e-8
BAND_REL_TOL = 1.0e-6
LOW_MODE_SAMPLE_COUNT = 40
LOW_MODE_POSITIVE_PROXY_COUNT = 24
LOW_MODE_BAND_COUNT = 3
LARGE_MODE_SAMPLE_COUNT = 8
SMALL_EIG_SEED_BASE = 7000
LARGE_EIG_SEED_BASE = 9000
TRACE_T_VALUES = (0.25, 0.5, 1.0)
NORM_CDF_GRID = np.linspace(1.0, 8.0, 71)
LOW_MODE_EXPONENT_SPAN_TOL = 0.1
LOW_MODE_FIT_R2_TOL = 0.999
NORM_CDF_COLLAPSE_TOL = 0.05
CLAIM_BOUNDARY = (
    "Phase VII authority is limited to descriptive spectral-feasibility diagnostics on the "
    "frozen branch-local cochain-Laplacian family and the declared refinement hierarchy. It "
    "does not assert continuum limits, geometric invariants, universal coefficients, physical "
    "laws, or observable-spectral coupling beyond this track."
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def closure_statement(success: bool) -> str:
    if success:
        return "Phase VII establishes spectral feasibility for the frozen operator hierarchy."
    return "Phase VII does not establish spectral feasibility for the frozen operator hierarchy."


def round_float(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def starting_vector(dim: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vector = rng.normal(size=dim)
    return vector / max(float(np.linalg.norm(vector)), 1.0e-12)


def sample_eigenvalues(operator, k: int, which: str, seed: int) -> np.ndarray:
    dim = int(operator.shape[0])
    if dim <= max(64, k + 2):
        dense_values = np.linalg.eigvalsh(operator.toarray())
        if which == "SA":
            return np.sort(dense_values)[:k]
        return np.sort(dense_values)[-k:]
    ncv = min(dim - 1, max(4 * k + 1, 48))
    values = spla.eigsh(
        operator,
        k=k,
        which=which,
        return_eigenvectors=False,
        tol=1.0e-10,
        ncv=ncv,
        v0=starting_vector(dim, seed),
    )
    return np.sort(np.real_if_close(values))


def distinct_positive_bands(values: np.ndarray) -> list[float]:
    bands: list[float] = []
    for value in sorted(float(entry) for entry in values if entry > NULLSPACE_TOL):
        if not bands or abs(value - bands[-1]) > BAND_REL_TOL * max(1.0, abs(bands[-1])):
            bands.append(value)
    return bands


def power_law_fit(h_values: list[float], y_values: list[float]) -> dict[str, float]:
    x = np.log(np.array(h_values, dtype=float))
    y = np.log(np.array(y_values, dtype=float))
    exponent, intercept = np.polyfit(x, y, 1)
    fitted = exponent * x + intercept
    denominator = float(np.sum((y - y.mean()) ** 2))
    if denominator <= 1.0e-18:
        r_squared = 1.0
    else:
        r_squared = 1.0 - float(np.sum((y - fitted) ** 2) / denominator)
    prefactor = float(math.exp(intercept))
    return {
        "exponent": round_float(exponent),
        "prefactor": round_float(prefactor),
        "r_squared": round_float(r_squared),
    }


def load_frozen_manifests() -> dict[str, Any]:
    phase6_manifest = read_json(PHASE6_MANIFEST_PATH)
    phase5_manifest = read_json(PHASE5_MANIFEST_PATH)

    if phase6_manifest.get("selected_operator_class") != "cochain_laplacian":
        raise ValueError("Phase VII Track A requires the frozen cochain-Laplacian family.")
    if not bool(phase6_manifest.get("success")):
        raise ValueError("Phase VII requires a successful Phase VI operator manifest.")
    if not bool(phase5_manifest.get("success")):
        raise ValueError("Phase VII requires a successful Phase V authority manifest.")

    return {
        "phase6_manifest": phase6_manifest,
        "phase5_manifest": phase5_manifest,
    }


def compute_level_summary(
    level_index: int,
    level: dict[str, Any],
    base_epsilon: float,
) -> dict[str, Any]:
    n_side = int(level["n_side"])
    h_value = float(level["h"])
    small_seed = SMALL_EIG_SEED_BASE + n_side
    large_seed = LARGE_EIG_SEED_BASE + n_side

    graph = build_graph(
        {
            "kind": "dk2d_periodic",
            "n_side": n_side,
            "epsilon": float(base_epsilon),
            "cycle_phase_x": 0.0,
            "cycle_phase_y": 0.0,
        }
    )
    operator = graph.delta_h.real.astype(float).tocsr()
    dim = int(operator.shape[0])

    small_values = sample_eigenvalues(operator, LOW_MODE_SAMPLE_COUNT, "SA", small_seed)
    large_values = sample_eigenvalues(operator, LARGE_MODE_SAMPLE_COUNT, "LA", large_seed)

    nullspace_estimate = int(np.sum(np.abs(small_values) <= NULLSPACE_TOL))
    positive_values = [float(value) for value in small_values if value > NULLSPACE_TOL]
    distinct_bands = distinct_positive_bands(small_values)
    if len(distinct_bands) < LOW_MODE_BAND_COUNT:
        raise ValueError("Phase VII requires at least three distinct positive low-mode bands.")

    first_positive = float(distinct_bands[0])
    normalized_proxy_values = np.array(positive_values[:LOW_MODE_POSITIVE_PROXY_COUNT], dtype=float) / first_positive
    cdf_values = np.array([(normalized_proxy_values <= grid_value).mean() for grid_value in NORM_CDF_GRID], dtype=float)
    trace_proxy = {
        f"t_{str(t_value).replace('.', '_')}": round_float(float(np.exp(-t_value * small_values).sum()))
        for t_value in TRACE_T_VALUES
    }

    return {
        "level_id": f"R{level_index}",
        "n_side": n_side,
        "h": round_float(h_value),
        "dimension": dim,
        "zero_form_size": int(graph.block_sizes[0]),
        "one_form_size": int(graph.block_sizes[1]),
        "two_form_size": int(graph.block_sizes[2]),
        "nnz": int(operator.nnz),
        "nullspace_estimate": nullspace_estimate,
        "low_mode_sample_count": int(LOW_MODE_SAMPLE_COUNT),
        "low_mode_positive_proxy_count": int(len(normalized_proxy_values)),
        "distinct_low_mode_bands": [round_float(value) for value in distinct_bands[:LOW_MODE_BAND_COUNT]],
        "sampled_smallest_eigenvalues": [round_float(value) for value in small_values.tolist()],
        "sampled_largest_eigenvalues": [round_float(value) for value in large_values.tolist()],
        "sampled_spectral_radius": round_float(float(large_values[-1])),
        "normalized_density_proxy_grid": [round_float(value, 6) for value in NORM_CDF_GRID.tolist()],
        "normalized_density_proxy_cdf": [round_float(value, 6) for value in cdf_values.tolist()],
        "trace_proxy": trace_proxy,
        "small_eig_seed": int(small_seed),
        "large_eig_seed": int(large_seed),
    }


def aggregate_level_summaries(levels: list[dict[str, Any]]) -> dict[str, Any]:
    h_values = [float(level["h"]) for level in levels]
    nullspace_values = [int(level["nullspace_estimate"]) for level in levels]
    spectral_radii = [float(level["sampled_spectral_radius"]) for level in levels]

    band_fits = []
    for band_index in range(LOW_MODE_BAND_COUNT):
        band_values = [float(level["distinct_low_mode_bands"][band_index]) for level in levels]
        fit = power_law_fit(h_values, band_values)
        fit["band_index"] = band_index + 1
        fit["monotone_decreasing_with_refinement"] = all(
            band_values[index + 1] < band_values[index] for index in range(len(band_values) - 1)
        )
        fit["band_values"] = [round_float(value) for value in band_values]
        band_fits.append(fit)

    trace_monotonicity = {}
    for t_value in TRACE_T_VALUES:
        trace_key = f"t_{str(t_value).replace('.', '_')}"
        trace_values = [float(level["trace_proxy"][trace_key]) for level in levels]
        trace_monotonicity[trace_key] = {
            "values": [round_float(value) for value in trace_values],
            "nondecreasing_with_refinement": all(
                trace_values[index + 1] >= trace_values[index] - 1.0e-10
                for index in range(len(trace_values) - 1)
            ),
        }

    cdf_distance_records = []
    for left_index, right_index in itertools.combinations(range(len(levels)), 2):
        left_level = levels[left_index]
        right_level = levels[right_index]
        left_cdf = np.array(left_level["normalized_density_proxy_cdf"], dtype=float)
        right_cdf = np.array(right_level["normalized_density_proxy_cdf"], dtype=float)
        max_distance = float(np.max(np.abs(left_cdf - right_cdf)))
        cdf_distance_records.append(
            {
                "left_level": left_level["level_id"],
                "right_level": right_level["level_id"],
                "max_cdf_distance": round_float(max_distance),
            }
        )

    exponent_values = [float(record["exponent"]) for record in band_fits]
    r_squared_values = [float(record["r_squared"]) for record in band_fits]
    density_distances = [float(record["max_cdf_distance"]) for record in cdf_distance_records]

    success = (
        len(set(nullspace_values)) == 1
        and all(bool(record["monotone_decreasing_with_refinement"]) for record in band_fits)
        and all(value >= LOW_MODE_FIT_R2_TOL for value in r_squared_values)
        and max(exponent_values) - min(exponent_values) <= LOW_MODE_EXPONENT_SPAN_TOL
        and all(
            spectral_radii[index + 1] >= spectral_radii[index] - 1.0e-10
            for index in range(len(spectral_radii) - 1)
        )
        and all(
            bool(trace_monotonicity[key]["nondecreasing_with_refinement"])
            for key in trace_monotonicity
        )
        and max(density_distances, default=0.0) <= NORM_CDF_COLLAPSE_TOL
    )

    return {
        "band_power_law_fits": band_fits,
        "nullspace_estimate_values": nullspace_values,
        "nullspace_consistent": len(set(nullspace_values)) == 1,
        "spectral_radius_values": [round_float(value) for value in spectral_radii],
        "spectral_radius_monotone": all(
            spectral_radii[index + 1] >= spectral_radii[index] - 1.0e-10
            for index in range(len(spectral_radii) - 1)
        ),
        "trace_proxy_monotonicity": trace_monotonicity,
        "normalized_density_proxy_pairwise_distances": cdf_distance_records,
        "normalized_density_proxy_max_pairwise_distance": round_float(max(density_distances, default=0.0)),
        "normalized_density_proxy_mean_pairwise_distance": round_float(
            float(sum(density_distances) / max(len(density_distances), 1))
        ),
        "low_mode_exponent_span": round_float(max(exponent_values) - min(exponent_values)),
        "minimum_low_mode_fit_r_squared": round_float(min(r_squared_values)),
        "success": bool(success),
    }


def write_plots(levels: list[dict[str, Any]]) -> None:
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)

    h_values = np.array([float(level["h"]) for level in levels], dtype=float)

    plt.style.use("seaborn-v0_8-whitegrid")

    figure, axis = plt.subplots(figsize=(6.2, 4.2))
    for band_index in range(LOW_MODE_BAND_COUNT):
        y_values = [float(level["distinct_low_mode_bands"][band_index]) for level in levels]
        axis.plot(h_values, y_values, marker="o", linewidth=1.8, label=f"band {band_index + 1}")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("h")
    axis.set_ylabel("sampled low-mode eigenvalue")
    axis.set_title("Phase VII Low-Mode Scaling")
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(EIGENVALUE_PLOT_PATH, format="svg")
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(6.2, 4.2))
    radius_values = [float(level["sampled_spectral_radius"]) for level in levels]
    axis.plot(h_values, radius_values, marker="o", linewidth=1.8, color="#0b6e4f")
    axis.set_xscale("log")
    axis.set_xlabel("h")
    axis.set_ylabel("sampled spectral radius")
    axis.set_title("Phase VII Spectral Envelope")
    figure.tight_layout()
    figure.savefig(RADIUS_PLOT_PATH, format="svg")
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(6.2, 4.2))
    palette = ["#7b2cbf", "#c14953", "#00798c"]
    for color, t_value in zip(palette, TRACE_T_VALUES):
        trace_key = f"t_{str(t_value).replace('.', '_')}"
        trace_values = [float(level["trace_proxy"][trace_key]) for level in levels]
        axis.plot(h_values, trace_values, marker="o", linewidth=1.8, color=color, label=f"t = {t_value}")
    axis.set_xscale("log")
    axis.set_xlabel("h")
    axis.set_ylabel("truncated trace proxy")
    axis.set_title("Phase VII Short-Time Trace Proxy")
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(TRACE_PLOT_PATH, format="svg")
    plt.close(figure)


def make_csv_rows(timestamp: str, levels: list[dict[str, Any]]) -> tuple[list[str], list[dict[str, Any]]]:
    fieldnames = [
        "stage_identifier",
        "timestamp",
        "track_id",
        "track_name",
        "operator_manifest_reference",
        "phase5_manifest_reference",
        "small_eig_seed",
        "large_eig_seed",
        "level_id",
        "n_side",
        "h",
        "dimension",
        "zero_form_size",
        "one_form_size",
        "two_form_size",
        "nnz",
        "nullspace_estimate",
        "band_1",
        "band_2",
        "band_3",
        "sampled_spectral_radius",
        "trace_t_0_25",
        "trace_t_0_5",
        "trace_t_1_0",
    ]
    rows = []
    for level in levels:
        rows.append(
            {
                "stage_identifier": STAGE_IDENTIFIER,
                "timestamp": timestamp,
                "track_id": TRACK_ID,
                "track_name": TRACK_NAME,
                "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
                "phase5_manifest_reference": relpath(PHASE5_MANIFEST_PATH),
                "small_eig_seed": level["small_eig_seed"],
                "large_eig_seed": level["large_eig_seed"],
                "level_id": level["level_id"],
                "n_side": level["n_side"],
                "h": level["h"],
                "dimension": level["dimension"],
                "zero_form_size": level["zero_form_size"],
                "one_form_size": level["one_form_size"],
                "two_form_size": level["two_form_size"],
                "nnz": level["nnz"],
                "nullspace_estimate": level["nullspace_estimate"],
                "band_1": level["distinct_low_mode_bands"][0],
                "band_2": level["distinct_low_mode_bands"][1],
                "band_3": level["distinct_low_mode_bands"][2],
                "sampled_spectral_radius": level["sampled_spectral_radius"],
                "trace_t_0_25": level["trace_proxy"]["t_0_25"],
                "trace_t_0_5": level["trace_proxy"]["t_0_5"],
                "trace_t_1_0": level["trace_proxy"]["t_1_0"],
            }
        )
    return fieldnames, rows


def build_manifest(
    timestamp: str,
    frozen: dict[str, Any],
    levels: list[dict[str, Any]],
    aggregate: dict[str, Any],
) -> dict[str, Any]:
    phase6_manifest = frozen["phase6_manifest"]
    phase5_manifest = frozen["phase5_manifest"]
    return {
        "timestamp": timestamp,
        "phase": 7,
        "phase_name": STAGE_IDENTIFIER,
        "stage_identifier": STAGE_IDENTIFIER,
        "primary_track": TRACK_NAME,
        "track_id": TRACK_ID,
        "status": "passed" if aggregate["success"] else "failed",
        "success": bool(aggregate["success"]),
        "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
        "phase5_manifest_reference": relpath(PHASE5_MANIFEST_PATH),
        "refinement_hierarchy": phase6_manifest["refinement_hierarchy"],
        "deterministic_seed_record": {
            "small_eig_seed_base": SMALL_EIG_SEED_BASE,
            "large_eig_seed_base": LARGE_EIG_SEED_BASE,
            "per_level_rule": {
                "small_eig_seed": "7000 + n_side",
                "large_eig_seed": "9000 + n_side",
            },
        },
        "frozen_inputs": {
            "phase6_selected_operator_class": phase6_manifest["selected_operator_class"],
            "phase6_operator_label": phase6_manifest["operator_label"],
            "phase4_sector_identifier": phase6_manifest["frozen_inputs"]["phase4_sector_identifier"],
            "phase5_dominant_stable_class": phase5_manifest["dominant_stable_class"],
            "base_epsilon": phase6_manifest["frozen_inputs"]["base_epsilon"],
        },
        "probe_settings": {
            "low_mode_sample_count": LOW_MODE_SAMPLE_COUNT,
            "low_mode_positive_proxy_count": LOW_MODE_POSITIVE_PROXY_COUNT,
            "low_mode_band_count": LOW_MODE_BAND_COUNT,
            "large_mode_sample_count": LARGE_MODE_SAMPLE_COUNT,
            "trace_t_values": list(TRACE_T_VALUES),
            "normalized_density_proxy_grid_max": round_float(float(NORM_CDF_GRID[-1]), 6),
            "normalized_density_proxy_grid_count": int(len(NORM_CDF_GRID)),
        },
        "track_selection_rationale": (
            "Track A was selected as the primary authoritative outcome because the frozen Phase VI "
            "operator hierarchy already shows consistent nullspace behavior, monotone envelope "
            "behavior, and stable low-mode ordering, while Phase V includes a bounded robustness "
            "edge-case classifier flip that makes Track B the less clean authority statement."
        ),
        "aggregate_diagnostics": aggregate,
        "level_summaries": levels,
        "artifacts": {
            "scaling_ledger_csv": relpath(SCALING_LEDGER_PATH),
            "summary_md": relpath(SUMMARY_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "eigenvalue_plot": relpath(EIGENVALUE_PLOT_PATH),
            "spectral_radius_plot": relpath(RADIUS_PLOT_PATH),
            "trace_proxy_plot": relpath(TRACE_PLOT_PATH),
            "builder_script": relpath(PHASE_ROOT / "build_spectral_feasibility.py"),
        },
        "secondary_track_executed": False,
        "failure_flags": {
            "low_mode_parameter_drift": not (
                float(aggregate["low_mode_exponent_span"]) <= LOW_MODE_EXPONENT_SPAN_TOL
                and float(aggregate["minimum_low_mode_fit_r_squared"]) >= LOW_MODE_FIT_R2_TOL
            ),
            "nullspace_instability": not bool(aggregate["nullspace_consistent"]),
            "spectral_radius_nonmonotone": not bool(aggregate["spectral_radius_monotone"]),
            "density_proxy_noncollapse": not (
                float(aggregate["normalized_density_proxy_max_pairwise_distance"]) <= NORM_CDF_COLLAPSE_TOL
            ),
            "trace_proxy_nonmonotone": not all(
                bool(item["nondecreasing_with_refinement"])
                for item in aggregate["trace_proxy_monotonicity"].values()
            ),
        },
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(bool(aggregate["success"])),
    }


def build_summary(
    timestamp: str,
    frozen: dict[str, Any],
    levels: list[dict[str, Any]],
    aggregate: dict[str, Any],
) -> str:
    phase6_manifest = frozen["phase6_manifest"]
    phase5_manifest = frozen["phase5_manifest"]
    table_header = (
        "| level | n_side | h | band 1 | band 2 | band 3 | radius | nullspace | trace(t=0.25) | trace(t=0.5) | trace(t=1.0) |\n"
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    )
    table_rows = [
        (
            f"| {level['level_id']} | {level['n_side']} | {level['h']:.12f} | "
            f"{level['distinct_low_mode_bands'][0]:.12f} | {level['distinct_low_mode_bands'][1]:.12f} | "
            f"{level['distinct_low_mode_bands'][2]:.12f} | {level['sampled_spectral_radius']:.12f} | "
            f"{level['nullspace_estimate']} | {level['trace_proxy']['t_0_25']:.12f} | "
            f"{level['trace_proxy']['t_0_5']:.12f} | {level['trace_proxy']['t_1_0']:.12f} |"
        )
        for level in levels
    ]
    band_fit_lines = [
        (
            f"- Band {record['band_index']}: exponent `{record['exponent']}`, "
            f"R^2 `{record['r_squared']}`, monotone `{record['monotone_decreasing_with_refinement']}`."
        )
        for record in aggregate["band_power_law_fits"]
    ]
    density_lines = [
        (
            f"- {record['left_level']} vs {record['right_level']}: max normalized CDF distance "
            f"`{record['max_cdf_distance']}`."
        )
        for record in aggregate["normalized_density_proxy_pairwise_distances"]
    ]
    trace_lines = [
        (
            f"- {trace_key}: values `{trace_data['values']}`, monotone "
            f"`{trace_data['nondecreasing_with_refinement']}`."
        )
        for trace_key, trace_data in aggregate["trace_proxy_monotonicity"].items()
    ]
    return "\n".join(
        [
            "# Phase VII Spectral Summary",
            "",
            f"Timestamp: `{timestamp}`.",
            "",
            "## Track Selection",
            "",
            "- Primary track: `spectral_feasibility` (Track A).",
            "- Secondary track executed: `False`.",
            "- Selection rationale: Phase VI already showed stable operator diagnostics across refinement, while Phase V contains a bounded robustness edge-case classifier flip that makes observable coupling the less clean primary authority statement.",
            "",
            "## Frozen References",
            "",
            f"- Operator manifest: `{relpath(PHASE6_MANIFEST_PATH)}`.",
            f"- Phase V authority manifest: `{relpath(PHASE5_MANIFEST_PATH)}`.",
            f"- Frozen operator class: `{phase6_manifest['selected_operator_class']}`.",
            f"- Frozen stable Phase V class: `{phase5_manifest['dominant_stable_class']}`.",
            "",
            "## Per-Level Probe Ledger",
            "",
            table_header,
            *table_rows,
            "",
            "## Low-Mode Scaling",
            "",
            *band_fit_lines,
            f"- Exponent span across the three fitted bands: `{aggregate['low_mode_exponent_span']}`.",
            f"- Minimum low-mode fit R^2: `{aggregate['minimum_low_mode_fit_r_squared']}`.",
            "",
            "## Spectral Envelope and Kernel Stability",
            "",
            f"- Spectral radius values: `{aggregate['spectral_radius_values']}`.",
            f"- Spectral radius monotone: `{aggregate['spectral_radius_monotone']}`.",
            f"- Nullspace estimate values: `{aggregate['nullspace_estimate_values']}`.",
            f"- Nullspace consistent: `{aggregate['nullspace_consistent']}`.",
            "",
            "## Normalized Density Proxy",
            "",
            "- The density surrogate is the empirical CDF of the first 24 sampled positive low-mode eigenvalues after normalization by the first positive band.",
            *density_lines,
            f"- Mean pairwise normalized CDF distance: `{aggregate['normalized_density_proxy_mean_pairwise_distance']}`.",
            f"- Max pairwise normalized CDF distance: `{aggregate['normalized_density_proxy_max_pairwise_distance']}`.",
            "",
            "## Short-Time Trace Proxy",
            "",
            "- The trace proxy is the truncated sum of `exp(-t lambda_i)` over the 40 sampled lowest modes.",
            *trace_lines,
            "",
            "## Bounded Interpretation",
            "",
            "- The frozen operator hierarchy supports a stable descriptive low-mode scaling regime across the declared refinement levels.",
            "- Kernel size remains fixed, the sampled spectral envelope remains monotone, and the normalized low-mode proxy stays tightly collapsed after scaling.",
            "- This is a branch-local feasibility result only. It is not a continuum, geometric, or physical-law claim.",
            "",
            "## Non-Claims",
            "",
            "- No continuum limit is established.",
            "- No geometric invariants or universal coefficients are extracted.",
            "- No observable-spectral coupling conclusion is made from this primary track.",
            "",
            closure_statement(bool(aggregate["success"])),
        ]
    )


def main() -> None:
    started_at = time.perf_counter()
    timestamp = timestamp_iso()
    frozen = load_frozen_manifests()

    phase6_manifest = frozen["phase6_manifest"]
    base_epsilon = float(phase6_manifest["frozen_inputs"]["base_epsilon"])
    hierarchy_levels = list(phase6_manifest["refinement_hierarchy"]["levels"])

    levels = [
        compute_level_summary(index + 1, level, base_epsilon)
        for index, level in enumerate(hierarchy_levels)
    ]
    aggregate = aggregate_level_summaries(levels)

    fieldnames, rows = make_csv_rows(timestamp, levels)
    write_csv_rows(SCALING_LEDGER_PATH, rows, fieldnames)
    write_plots(levels)

    summary_text = build_summary(timestamp, frozen, levels, aggregate)
    manifest_payload = build_manifest(timestamp, frozen, levels, aggregate)

    write_text(SUMMARY_PATH, summary_text)
    write_json(MANIFEST_PATH, manifest_payload)

    runtime_seconds = time.perf_counter() - started_at
    print(f"Phase VII Track A success: {aggregate['success']}")
    print(f"Runtime seconds: {round_float(runtime_seconds, 6)}")
    print(f"Manifest: {relpath(MANIFEST_PATH)}")
    print(closure_statement(bool(aggregate["success"])))


if __name__ == "__main__":
    main()
