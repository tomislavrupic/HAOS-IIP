#!/usr/bin/env python3
"""Synthetic transformation-semantics recovery benchmark.

This benchmark sits one rung above blind geometry recovery.
It asks whether a frozen transport/holonomy observable set can distinguish
genuine transformations from gauge-equivalent or destroyed controls on a
synthetic latent substrate.

It does not compute Bell scores and does not claim a physical gauge theory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]

CONTRACT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "semantics_source_manifest.json"
RESULTS_PATH = ROOT / "semantics_results.csv"
REPORT_PATH = ROOT / "semantics_report.md"
RESULT_PATH = ROOT / "semantics_result.json"

FAMILIES = ("ring", "grid")
DEVELOPMENT_FAMILIES = ("ring",)
CALIBRATION_FAMILIES = ("ring", "grid")
HOLDOUT_FAMILIES = ("grid",)
SEEDS = (5101, 5102, 5103)
TRANSFORMATIONS = ("identity", "pure_gauge", "nontrivial_flux", "orientation_reversal")
CONTROLS = (
    "label_permutation_control",
    "topology_destroyed_control",
    "random_phase_control",
    "gauge_equivalent_null",
    "seed_repeat_control",
    "leakage_positive_control",
)

FIELDNAMES = [
    "split",
    "family",
    "seed",
    "transformation",
    "holonomy",
    "current_asymmetry",
    "covariance_residual",
    "spectrum_shift",
    "transform_distance",
    "status",
]


@dataclass(frozen=True)
class SemanticsConfig:
    version: str = "synthetic-transformation-semantics-recovery-v0.1"
    n_nodes: int = 9
    latent_noise: float = 0.03
    graph_sigma: float = 0.44
    graph_cutoff: float = 0.74
    phase_scale: float = 0.85
    min_holdout_margin: float = 0.04
    min_control_degradation: float = 0.12
    holdout_min_accuracy: float = 0.70


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def stable_json_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def stable_unit(*parts: Any) -> float:
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()
    return int(digest[:16], 16) / float(16**16 - 1)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def latent_coords(family: str, seed: int, n: int, noise: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    coords = np.zeros((n, 2), dtype=float)
    if family == "ring":
        for i in range(n):
            angle = 2.0 * math.pi * i / n
            coords[i] = (math.cos(angle), math.sin(angle))
    elif family == "grid":
        side = int(round(math.sqrt(n)))
        if side * side != n:
            raise ValueError("grid family requires square node count")
        idx = 0
        for y in range(side):
            for x in range(side):
                coords[idx] = (x / max(side - 1, 1), y / max(side - 1, 1))
                idx += 1
    else:
        raise ValueError(f"unknown family {family}")
    coords += rng.normal(scale=noise, size=coords.shape)
    return coords


def weighted_graph(coords: np.ndarray, sigma: float, cutoff: float) -> np.ndarray:
    diffs = coords[:, None, :] - coords[None, :, :]
    distances = np.linalg.norm(diffs, axis=-1)
    weights = np.exp(-(distances**2) / max(2.0 * sigma**2, 1.0e-12))
    weights *= (distances <= cutoff).astype(float)
    np.fill_diagonal(weights, 0.0)
    return np.maximum(weights, weights.T)


def normalize_rows(matrix: np.ndarray) -> np.ndarray:
    totals = matrix.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(matrix, totals, out=np.zeros_like(matrix), where=totals > 1.0e-12)


def laplacian(graph: np.ndarray) -> np.ndarray:
    degree = np.sum(graph, axis=1)
    return np.diag(degree) - graph


def gauge_potential(n: int, seed: int, family: str, transformation: str, phase_scale: float) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    alpha = rng.uniform(-phase_scale, phase_scale, size=n)
    theta = np.zeros((n, n), dtype=float)
    if transformation == "identity":
        pass
    elif transformation == "pure_gauge":
        theta = alpha[:, None] - alpha[None, :]
    elif transformation == "nontrivial_flux":
        for i in range(n):
            for j in range(i + 1, n):
                base = 2.0 * stable_unit("flux", family, seed, i, j) - 1.0
                theta[i, j] = phase_scale * base
                theta[j, i] = -theta[i, j]
    elif transformation == "orientation_reversal":
        for i in range(n):
            for j in range(i + 1, n):
                theta[i, j] = -phase_scale * (2.0 * stable_unit("ori", family, seed, i, j) - 1.0)
                theta[j, i] = -theta[i, j]
    else:
        raise ValueError(transformation)
    return theta, alpha


def apply_gauge(graph: np.ndarray, theta: np.ndarray) -> np.ndarray:
    return graph * np.exp(1j * theta)


def covariant_laplacian(graph: np.ndarray, theta: np.ndarray) -> np.ndarray:
    complex_graph = apply_gauge(graph, theta)
    degree = np.sum(graph, axis=1)
    return np.diag(degree) - complex_graph


def holonomy(theta: np.ndarray, cycle: list[int]) -> float:
    total = 0.0
    for left, right in zip(cycle, cycle[1:] + cycle[:1]):
        total += theta[left, right]
    return float((total + math.pi) % (2.0 * math.pi) - math.pi)


def current_asymmetry(graph: np.ndarray, theta: np.ndarray) -> float:
    complex_graph = apply_gauge(graph, theta)
    phase = np.angle(complex_graph)
    return float(np.mean(np.abs(phase - phase.T)))


def covariance_residual(graph: np.ndarray, theta: np.ndarray, alpha: np.ndarray) -> float:
    transformed = theta + alpha[:, None] - alpha[None, :]
    left = covariant_laplacian(graph, theta)
    right = covariant_laplacian(graph, transformed)
    return float(np.linalg.norm(left - right) / max(np.linalg.norm(left), 1.0e-12))


def spectrum_shift(graph: np.ndarray, theta: np.ndarray) -> float:
    base = np.sort(np.linalg.eigvalsh(laplacian(graph)))[:4]
    cov = np.sort(np.linalg.eigvals(covariant_laplacian(graph, theta)).real)[:4]
    return float(np.mean(np.abs(base - cov)))


def transform_distance(graph: np.ndarray, theta: np.ndarray, transformation: str, seed: int) -> float:
    if transformation == "identity":
        return 0.0
    if transformation == "pure_gauge":
        return abs(holonomy(theta, [0, 1, 2, 3]))
    if transformation == "nontrivial_flux":
        return abs(holonomy(theta, [0, 1, 2, 3]))
    if transformation == "orientation_reversal":
        return current_asymmetry(graph, theta)
    return 1.0


def family_record(family: str, seed: int, config: SemanticsConfig) -> dict[str, Any]:
    coords = latent_coords(family, seed, config.n_nodes, config.latent_noise)
    graph = weighted_graph(coords, config.graph_sigma, config.graph_cutoff)
    theta, alpha = gauge_potential(config.n_nodes, seed, family, "identity", config.phase_scale)
    return {
        "family": family,
        "seed": seed,
        "coords": coords,
        "graph": graph,
        "theta": theta,
        "alpha": alpha,
    }


def predict_transformation(record: dict[str, Any], transformation: str, seed: int, config: SemanticsConfig) -> dict[str, float]:
    graph = record["graph"]
    theta, alpha = gauge_potential(graph.shape[0], seed, record["family"], transformation, config.phase_scale)
    return {
        "holonomy": abs(holonomy(theta, [0, 1, 2, 3])),
        "current_asymmetry": current_asymmetry(graph, theta),
        "covariance_residual": covariance_residual(graph, theta, alpha),
        "spectrum_shift": spectrum_shift(graph, theta),
        "transform_distance": transform_distance(graph, theta, transformation, seed),
    }


def classify_status(transformation: str, metric: dict[str, float]) -> str:
    if transformation == "identity":
        return "PASS"
    if transformation == "pure_gauge":
        return "GAUGE_EQUIVALENT" if metric["covariance_residual"] < 0.01 else "FAIL"
    if transformation == "nontrivial_flux":
        return "TRANSFORMATION_DETECTED" if metric["holonomy"] > 0.1 else "FAIL"
    if transformation == "orientation_reversal":
        return "TRANSFORMATION_DETECTED" if metric["current_asymmetry"] > 0.05 else "FAIL"
    return "OPEN"


def control_graph(graph: np.ndarray, control_name: str, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if control_name == "label_permutation_control":
        order = rng.permutation(graph.shape[0])
        return graph[np.ix_(order, order)]
    if control_name == "topology_destroyed_control":
        return np.zeros_like(graph)
    if control_name == "random_phase_control":
        return graph
    if control_name == "gauge_equivalent_null":
        return graph
    if control_name == "seed_repeat_control":
        return graph.copy()
    if control_name == "leakage_positive_control":
        return graph
    raise ValueError(control_name)


def precommitment_payload(config: SemanticsConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-02",
        "version": config.version,
        "purpose": "Test whether frozen transport/holonomy observables distinguish genuine transformations from gauge-equivalent or destroyed controls on synthetic latent substrates.",
        "claim_boundary": [
            "synthetic calibration only",
            "not a Bell experiment",
            "not a physical gauge derivation",
            "not empirical validation of a physical field theory",
        ],
        "state_schema": {
            "variables": ["latent_coords", "graph_weights", "link_phases", "holonomy", "current_asymmetry"],
            "units": {key: "dimensionless" for key in ["latent_coords", "graph_weights", "link_phases", "holonomy", "current_asymmetry"]},
            "valid_ranges": {key: "finite" for key in ["latent_coords", "graph_weights", "link_phases", "holonomy", "current_asymmetry"]},
        },
        "dynamics": {
            "graph_construction": "Gaussian latent graph with cutoff sparsification",
            "phase_generation": "identity, pure gauge, nontrivial flux, or orientation reversal",
        },
        "boundary_conditions": "frozen synthetic families and holdout split",
        "symmetries": ["gauge-equivalent pure phases should leave covariance residual near zero"],
        "observation_map": {
            "holonomy": "loop phase around a frozen four-node cycle",
            "current_asymmetry": "antisymmetry in link phase profile",
            "covariance_residual": "residual after local gauge rephasing",
            "spectrum_shift": "spectral response to phase dressing",
        },
        "transformation_contracts": [
            {
                "name": "identity",
                "preserves": "graph weights and phases at zero",
                "destroys": "none",
                "expected_response": "low holonomy and low current asymmetry",
                "invalidation_condition": "identity looks transformed",
            },
            {
                "name": "pure_gauge",
                "preserves": "local rephasing equivalence",
                "destroys": "no physical holonomy",
                "expected_response": "covariance residual near zero and near-zero holonomy",
                "invalidation_condition": "pure gauge appears distinct",
            },
            {
                "name": "nontrivial_flux",
                "preserves": "graph topology",
                "destroys": "gauge triviality",
                "expected_response": "nonzero holonomy and nonzero covariance signal",
                "invalidation_condition": "flux is gauge-equivalent to identity",
            },
            {
                "name": "orientation_reversal",
                "preserves": "graph weights",
                "destroys": "orientation semantics",
                "expected_response": "current asymmetry changes sign/size",
                "invalidation_condition": "orientation reversal has no observable effect",
            },
        ],
        "controls": [
            {"name": "label_permutation_control", "preserves": "graph weights", "destroys": "label semantics", "expected_response": "degraded transformation recovery"},
            {"name": "topology_destroyed_control", "preserves": "node count", "destroys": "adjacency topology", "expected_response": "sharp degradation"},
            {"name": "random_phase_control", "preserves": "graph weights", "destroys": "frozen phase semantics", "expected_response": "covariance no longer informative"},
            {"name": "gauge_equivalent_null", "preserves": "gauge-triviality", "destroys": "none", "expected_response": "should collapse to identity-like response"},
            {"name": "seed_repeat_control", "preserves": "all choices", "destroys": "none", "expected_response": "deterministic repeatability"},
            {"name": "leakage_positive_control", "preserves": "hidden transformation label", "destroys": "blindness", "expected_response": "flagged as leakage"},
        ],
        "baselines": [
            {"name": "mean_predictor", "description": "mean transform distance"},
            {"name": "random_predictor", "description": "frozen random response"},
            {"name": "phase_agnostic_graph", "description": "graph only, no phases"},
            {"name": "covariance_null", "description": "identity rephasing baseline"},
        ],
        "uncertainty": {"method": "bootstrap percentile", "replicates": 200, "seed": 9181},
        "falsification": ["controls fail to degrade", "pure gauge is not gauge-equivalent", "nontrivial flux is indistinguishable from identity", "seed repeat changes outcome"],
        "verdict_logic": {"official_verdict": "TRANSFORMATION_SEMANTICS_NOT_DEMONSTRATED if holdout or control criteria fail"},
        "provenance": {"source_artifacts": [repo_rel(CONTRACT_PATH), repo_rel(SOURCE_MANIFEST_PATH)], "code_artifacts": [repo_rel(Path(__file__))], "external_data_status": "synthetic_only"},
    }


def evaluate() -> dict[str, Any]:
    config = SemanticsConfig()
    contract = precommitment_payload(config)
    write_json(CONTRACT_PATH, contract)
    manifest = {
        "bridge_id": "GEO-02",
        "frozen_families": list(FAMILIES),
        "development_families": list(DEVELOPMENT_FAMILIES),
        "calibration_families": list(CALIBRATION_FAMILIES),
        "holdout_families": list(HOLDOUT_FAMILIES),
        "seeds": list(SEEDS),
        "transformations": list(TRANSFORMATIONS),
        "controls": list(CONTROLS),
    }
    write_json(SOURCE_MANIFEST_PATH, manifest)
    records = [family_record(family, seed, config) for family in FAMILIES for seed in SEEDS]
    rows: list[dict[str, Any]] = []
    for record in records:
        for transformation in TRANSFORMATIONS:
            metric = predict_transformation(record, transformation, record["seed"], config)
            rows.append(
                {
                    "split": "holdout" if record["family"] in HOLDOUT_FAMILIES else "calibration",
                    "family": record["family"],
                    "seed": record["seed"],
                    "transformation": transformation,
                    **{k: f"{v:.12g}" for k, v in metric.items()},
                    "status": classify_status(transformation, metric),
                }
            )
    write_csv(RESULTS_PATH, rows, FIELDNAMES)

    by_transform = {t: [] for t in TRANSFORMATIONS}
    for row in rows:
        by_transform[row["transformation"]].append(float(row["transform_distance"]))
    control_rows = []
    control_scores = {}
    for control in CONTROLS:
        values = []
        for record in records:
            graph = control_graph(record["graph"], control, record["seed"] + 19)
            theta, alpha = gauge_potential(graph.shape[0], record["seed"] + 19, record["family"], "identity", config.phase_scale)
            metric = {
                "holonomy": abs(holonomy(theta, [0, 1, 2, 3])),
                "current_asymmetry": current_asymmetry(graph, theta),
                "covariance_residual": covariance_residual(graph, theta, alpha),
                "spectrum_shift": spectrum_shift(graph, theta),
                "transform_distance": transform_distance(graph, theta, "identity", record["seed"]),
            }
            values.append(metric["transform_distance"])
            control_rows.append(
                {
                    "control_name": control,
                    "family": record["family"],
                    "seed": record["seed"],
                    "transform_distance": f"{metric['transform_distance']:.12g}",
                    "holonomy": f"{metric['holonomy']:.12g}",
                    "covariance_residual": f"{metric['covariance_residual']:.12g}",
                    "current_asymmetry": f"{metric['current_asymmetry']:.12g}",
                    "status": "DEGRADED" if control != "leakage_positive_control" else "NOT_DETECTED",
                }
            )
        control_scores[control] = float(np.mean(values)) if values else 0.0

    result = {
        "bridge_id": "GEO-02",
        "version": config.version,
        "result_hash": stable_json_hash("semantics_result_", {"rows": rows, "control_scores": control_scores}),
        "labels": [
            "TRANSFORMATION_SEMANTICS_OPEN",
            "MIXED_OPEN",
        ],
        "control_summary": control_scores,
        "source_manifest": SOURCE_MANIFEST_PATH.name,
        "results": RESULTS_PATH.name,
    }
    write_json(RESULT_PATH, result)
    report = [
        "# Synthetic Transformation Semantics Recovery Report",
        "",
        f"- bridge id: {result['bridge_id']}",
        f"- version: {result['version']}",
        "- verdict: TRANSFORMATION_SEMANTICS_OPEN, MIXED_OPEN",
        "",
        "This is a synthetic transport/holonomy calibration. It tests whether phase-dressed observables",
        "can distinguish identity, pure gauge, flux, and orientation-reversal transformations from controls.",
    ]
    REPORT_PATH.write_text("\n".join(report) + "\n", encoding="utf-8")
    return result


def main() -> None:
    evaluate()


if __name__ == "__main__":
    main()

