#!/usr/bin/env python3
"""HAOS-IIP effective-operator expansion scaffold.

This benchmark borrows EFT discipline: explicit scale, cutoff, allowed terms,
coarse observables, controls, and claim ceilings. It does not derive quantum
gravity, spacetime, matter fields, or empirical physics.

The first benchmark is intentionally boring: a deterministic ring operator with
a diffusion-like leading term and a higher-order correction. The question is
whether the fitted coarse dynamics recover the term hierarchy under a frozen
synthetic setup.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"

PRECOMMITMENT_PATH = RESULTS / "precommitment_contract.json"
ALLOWED_TERMS_PATH = RESULTS / "allowed_terms.json"
SWEEP_PATH = RESULTS / "coefficient_sweep.csv"
CONTROL_PATH = RESULTS / "control_results.csv"
REPORT_PATH = RESULTS / "effective_operator_report.md"
RESULT_PATH = RESULTS / "effective_operator_result.json"

SWEEP_FIELDS = [
    "mode",
    "lambda",
    "observed_decay_rate",
    "predicted_decay_rate",
    "residual",
    "relative_residual",
    "scale_class",
]
CONTROL_FIELDS = [
    "control",
    "preserved",
    "destroyed",
    "expected_response",
    "observed_leading_coefficient",
    "observed_correction_coefficient",
    "r2",
    "observed_status",
]


@dataclass(frozen=True)
class EffectiveOperatorConfig:
    version: str = "effective-operator-expansion-v0.1"
    node_count: int = 64
    steps: int = 72
    dt: float = 0.08
    diffusion_coeff: float = 0.18
    correction_coeff: float = 0.012
    mode_min: int = 1
    mode_max: int = 12
    long_wavelength_mode_max: int = 4
    max_long_mode_relative_residual: float = 0.035
    min_fit_r2: float = 0.995
    min_leading_coefficient: float = 0.05
    max_correction_ratio: float = 0.35
    unstable_growth_tolerance: float = 1.0e-9


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run restrained HAOS effective-operator expansion scaffold.")
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    return parser.parse_args()


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


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


def ring_laplacian(node_count: int) -> np.ndarray:
    lap = np.zeros((node_count, node_count), dtype=float)
    for index in range(node_count):
        lap[index, index] = 2.0
        lap[index, (index - 1) % node_count] = -1.0
        lap[index, (index + 1) % node_count] = -1.0
    return lap


def laplacian_eigenvalue(node_count: int, mode: int) -> float:
    return float(2.0 - 2.0 * math.cos(2.0 * math.pi * mode / node_count))


def mode_vector(node_count: int, mode: int) -> np.ndarray:
    x = np.arange(node_count, dtype=float)
    vector = np.cos(2.0 * math.pi * mode * x / node_count)
    norm = float(np.linalg.norm(vector))
    return vector / norm


def operator_matrix(laplacian: np.ndarray, config: EffectiveOperatorConfig, *, sign: float = 1.0) -> np.ndarray:
    return sign * (config.diffusion_coeff * laplacian + config.correction_coeff * (laplacian @ laplacian))


def simulate_decay(matrix: np.ndarray, initial: np.ndarray, config: EffectiveOperatorConfig) -> tuple[np.ndarray, list[float]]:
    state = initial.copy()
    amplitudes = [float(np.dot(initial, state))]
    update = np.eye(matrix.shape[0]) - config.dt * matrix
    for _ in range(config.steps):
        state = update @ state
        amplitudes.append(float(np.dot(initial, state)))
    return state, amplitudes


def fit_decay_rate(amplitudes: list[float], dt: float) -> float:
    times = np.arange(len(amplitudes), dtype=float) * dt
    safe = np.maximum(np.abs(np.array(amplitudes, dtype=float)), 1.0e-14)
    slope, _ = np.polyfit(times, np.log(safe), 1)
    return float(-slope)


def rows_for_operator(matrix: np.ndarray, config: EffectiveOperatorConfig) -> list[dict[str, Any]]:
    rows = []
    for mode in range(config.mode_min, config.mode_max + 1):
        initial = mode_vector(config.node_count, mode)
        _, amplitudes = simulate_decay(matrix, initial, config)
        observed = fit_decay_rate(amplitudes, config.dt)
        lam = laplacian_eigenvalue(config.node_count, mode)
        rows.append(
            {
                "mode": mode,
                "lambda": lam,
                "observed_decay_rate": observed,
            }
        )
    return rows


def fit_terms(rows: list[dict[str, Any]]) -> dict[str, Any]:
    x = np.array([[row["lambda"], row["lambda"] ** 2] for row in rows], dtype=float)
    y = np.array([row["observed_decay_rate"] for row in rows], dtype=float)
    coeffs, *_ = np.linalg.lstsq(x, y, rcond=None)
    predicted = x @ coeffs
    residuals = y - predicted
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 1.0e-15 else float(1.0 - ss_res / ss_tot)
    return {
        "leading_laplacian_coeff": float(coeffs[0]),
        "higher_laplacian_squared_coeff": float(coeffs[1]),
        "r2": r2,
        "predicted": predicted,
        "residuals": residuals,
    }


def annotate_rows(rows: list[dict[str, Any]], fit: dict[str, Any], config: EffectiveOperatorConfig) -> list[dict[str, Any]]:
    annotated = []
    for index, row in enumerate(rows):
        observed = float(row["observed_decay_rate"])
        predicted = float(fit["predicted"][index])
        residual = float(fit["residuals"][index])
        rel = abs(residual) / max(abs(observed), 1.0e-12)
        mode = int(row["mode"])
        annotated.append(
            {
                "mode": mode,
                "lambda": f"{float(row['lambda']):.12g}",
                "observed_decay_rate": f"{observed:.12g}",
                "predicted_decay_rate": f"{predicted:.12g}",
                "residual": f"{residual:.12g}",
                "relative_residual": f"{rel:.12g}",
                "scale_class": "long_wavelength" if mode <= config.long_wavelength_mode_max else "shorter_wavelength",
            }
        )
    return annotated


def precommitment(config: EffectiveOperatorConfig) -> dict[str, Any]:
    return {
        "bridge_id": "HAOS-EOT-01",
        "purpose": "Test whether a frozen synthetic operator hierarchy admits an EFT-like effective expansion.",
        "not_claimed": [
            "EFT quantum gravity",
            "physical field theory",
            "spacetime derivation",
            "Lorentz invariance",
            "empirical physics",
            "ontology",
        ],
        "base_operator_family": "1D ring graph Laplacian",
        "allowed_terms": [
            {"name": "identity", "order": 0, "status": "reporting_only"},
            {"name": "laplacian", "order": 2, "status": "leading_allowed"},
            {"name": "laplacian_squared", "order": 4, "status": "suppressed_correction_allowed"},
        ],
        "observables": [
            "mode_decay_rate",
            "coefficient_hierarchy",
            "long_wavelength_residual",
            "control_response",
        ],
        "cutoff": {
            "node_count": config.node_count,
            "mode_max": config.mode_max,
            "long_wavelength_mode_max": config.long_wavelength_mode_max,
        },
        "verdict_logic": {
            "effective_operator_scaffold_pass": [
                "fit r2 >= configured minimum",
                "leading coefficient positive",
                "higher-order correction smaller than configured ratio",
                "long-wavelength residuals below configured tolerance",
                "unstable-sign control is detected",
            ],
            "claim_ceiling": "synthetic effective-operator scaffold only",
        },
        "config": asdict(config),
    }


def allowed_terms_payload(config: EffectiveOperatorConfig) -> dict[str, Any]:
    return {
        "terms": precommitment(config)["allowed_terms"],
        "mapping_status": {
            "laplacian": "SYNTHETIC_DERIVED_FROM_GRAPH_OPERATOR",
            "laplacian_squared": "SYNTHETIC_ALLOWED_CORRECTION",
            "physical_field": "UNAVAILABLE",
            "metric": "UNAVAILABLE",
            "stress_energy": "UNAVAILABLE",
        },
        "claim_ceiling": "Allowed terms for synthetic operator expansion only.",
    }


def control_rows(laplacian: np.ndarray, reference_fit: dict[str, Any], config: EffectiveOperatorConfig) -> list[dict[str, Any]]:
    controls = []
    cases = [
        (
            "identity_operator_control",
            "node count and time grid",
            "spatial coupling and decay hierarchy",
            np.zeros_like(laplacian),
            "coefficients should collapse near zero",
        ),
        (
            "unstable_sign_control",
            "same eigenvectors and scale grid",
            "stable decay sign",
            operator_matrix(laplacian, config, sign=-1.0),
            "growth should be detected as invalid effective decay",
        ),
    ]
    for name, preserved, destroyed, matrix, expected in cases:
        rows = rows_for_operator(matrix, config)
        fit = fit_terms(rows)
        leading = float(fit["leading_laplacian_coeff"])
        correction = float(fit["higher_laplacian_squared_coeff"])
        r2 = float(fit["r2"])
        if name == "identity_operator_control":
            status = "PASS" if abs(leading) < config.min_leading_coefficient and abs(correction) < config.min_leading_coefficient else "FAIL_CONTROL"
        elif name == "unstable_sign_control":
            status = "PASS" if leading < -config.unstable_growth_tolerance else "FAIL_UNSTABLE_SIGN_CONTROL"
        else:
            status = "FAIL_UNKNOWN_CONTROL"
        controls.append(
            {
                "control": name,
                "preserved": preserved,
                "destroyed": destroyed,
                "expected_response": expected,
                "observed_leading_coefficient": f"{leading:.12g}",
                "observed_correction_coefficient": f"{correction:.12g}",
                "r2": f"{r2:.12g}",
                "observed_status": status,
            }
        )
    return controls


def run_effective_operator_expansion(config: EffectiveOperatorConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    lap = ring_laplacian(config.node_count)
    matrix = operator_matrix(lap, config)
    rows = rows_for_operator(matrix, config)
    fit = fit_terms(rows)
    annotated = annotate_rows(rows, fit, config)
    controls = control_rows(lap, fit, config)

    long_residuals = [float(row["relative_residual"]) for row in annotated if row["scale_class"] == "long_wavelength"]
    max_long_residual = max(long_residuals) if long_residuals else float("inf")
    leading = float(fit["leading_laplacian_coeff"])
    correction = float(fit["higher_laplacian_squared_coeff"])
    correction_ratio = abs(correction) / max(abs(leading), 1.0e-12)
    controls_pass = all(row["observed_status"] == "PASS" for row in controls)
    scaffold_pass = (
        fit["r2"] >= config.min_fit_r2
        and leading >= config.min_leading_coefficient
        and correction_ratio <= config.max_correction_ratio
        and max_long_residual <= config.max_long_mode_relative_residual
        and controls_pass
    )

    labels = [
        "EFFECTIVE_OPERATOR_EXPANSION_BUILT",
        "EFT_DISCIPLINE_METHOD_IMPORTED",
        "DIFFUSION_LIKE_LEADING_TERM_DETECTED" if leading > 0 else "DIFFUSION_LIKE_LEADING_TERM_NOT_DETECTED",
        "SUPPRESSED_CORRECTION_HIERARCHY_REPORTED",
        "CUTOFF_DECLARED",
        "CONTROL_RESPONSE_REPORTED",
        "EFT_QG_NOT_DERIVED",
        "PHYSICAL_FIELD_THEORY_NOT_ESTABLISHED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
    ]
    labels.append("EFFECTIVE_OPERATOR_SCAFFOLD_PASS" if scaffold_pass else "EFFECTIVE_OPERATOR_SCAFFOLD_OPEN")

    result = {
        "bridge_id": "effective_operator_expansion_v0_1",
        "status": "PASS" if scaffold_pass else "OPEN",
        "classification": "SYNTHETIC_EFFECTIVE_OPERATOR_SCAFFOLD",
        "labels": labels,
        "config": asdict(config),
        "fit": {
            "leading_laplacian_coeff": leading,
            "higher_laplacian_squared_coeff": correction,
            "correction_ratio": correction_ratio,
            "r2": float(fit["r2"]),
            "max_long_wavelength_relative_residual": max_long_residual,
        },
        "controls": controls,
        "claim_ceiling": "Synthetic EFT-like effective-operator expansion only; no EFT quantum gravity, physical field theory, spacetime, or ontology claim.",
    }
    result["result_hash"] = stable_hash("effective_operator", result)

    write_json(output_dir / PRECOMMITMENT_PATH.name, precommitment(config))
    write_json(output_dir / ALLOWED_TERMS_PATH.name, allowed_terms_payload(config))
    write_csv(output_dir / SWEEP_PATH.name, annotated, SWEEP_FIELDS)
    write_csv(output_dir / CONTROL_PATH.name, controls, CONTROL_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result)
    return result


def write_report(path: Path, result: dict[str, Any]) -> None:
    fit = result["fit"]
    lines = [
        "# Effective Operator Expansion Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        f"- leading Laplacian coefficient: `{fit['leading_laplacian_coeff']:.12g}`",
        f"- higher-order correction coefficient: `{fit['higher_laplacian_squared_coeff']:.12g}`",
        f"- correction ratio: `{fit['correction_ratio']:.12g}`",
        f"- fit R^2: `{fit['r2']:.12g}`",
        f"- max long-wavelength relative residual: `{fit['max_long_wavelength_relative_residual']:.12g}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Controls",
        "",
        "| Control | Leading coeff | Correction coeff | R^2 | Status |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for row in result["controls"]:
        lines.append(
            "| {control} | {observed_leading_coefficient} | {observed_correction_coefficient} | {r2} | `{observed_status}` |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This benchmark imports EFT discipline: scale hierarchy, allowed terms, cutoff, matching, and controls.",
            "It does not derive EFT quantum gravity, physical field theory, spacetime, Lorentz invariance, matter sectors, constants, or empirical physics.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_effective_operator_expansion(EffectiveOperatorConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
