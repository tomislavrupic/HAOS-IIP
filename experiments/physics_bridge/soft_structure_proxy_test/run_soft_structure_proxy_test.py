#!/usr/bin/env python3
"""Phase 59 toy soft-structure proxy test.

This sidecar does not test real celestial amplitudes or gravity. It creates a
known synthetic soft-pole target and asks whether soft-specific telemetry can
detect residue, pole, and factorization breakage better than generic smoothness
checks.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent

CSV_PATH = ROOT / "soft_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "soft_structure_report.md"
FIGURE_DIR = ROOT / "figures"

FIELDNAMES = [
    "sample",
    "sample_role",
    "control_type",
    "channel",
    "true_residue",
    "fitted_residue",
    "residue_relative_error",
    "fitted_pole",
    "residue_recoverability",
    "pole_location_score",
    "factorization_score",
    "generic_smoothness_score",
    "soft_specificity_score",
]


@dataclass(frozen=True)
class SoftConfig:
    seed: int = 5901
    omega_min: float = 0.002
    omega_max: float = 0.35
    omega_count: int = 80
    soft_fit_max_omega: float = 0.055
    pass_min_soft_score: float = 0.86
    pass_min_residue_recoverability: float = 0.92
    pass_min_pole_score: float = 0.88
    pass_min_factorization_score: float = 0.90
    pass_min_margin: float = 0.18
    pass_max_control_false_positive_rate: float = 0.0


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    sample_role: str
    control_type: str
    amplitudes: np.ndarray


@dataclass(frozen=True)
class SampleSummary:
    sample: str
    sample_role: str
    control_type: str
    residue_recoverability: float
    pole_location_score: float
    factorization_score: float
    generic_smoothness_score: float
    soft_specificity_score: float
    mean_abs_pole: float
    residue_error_norm: float


@dataclass(frozen=True)
class ChannelDiagnostic:
    sample: str
    sample_role: str
    control_type: str
    channel: str
    true_residue: float
    fitted_residue: float
    residue_relative_error: float
    fitted_pole: float
    residue_recoverability: float
    pole_location_score: float
    factorization_score: float
    generic_smoothness_score: float
    soft_specificity_score: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 59 toy soft-structure proxy test.")
    parser.add_argument("--seed", type=int, default=5901)
    parser.add_argument("--omega-count", type=int, default=80)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, JSON, markdown report, and figures.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> SoftConfig:
    return SoftConfig(seed=int(args.seed), omega_count=int(args.omega_count))


def channel_labels() -> list[str]:
    return ["12", "13", "14", "23", "24", "34"]


def hard_charges() -> np.ndarray:
    return np.array([1.0, 0.80, -0.60, 0.45], dtype=float)


def true_residues() -> np.ndarray:
    charges = hard_charges()
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    return np.array([charges[i] * charges[j] for i, j in pairs], dtype=float)


def omega_grid(config: SoftConfig) -> np.ndarray:
    return np.geomspace(config.omega_min, config.omega_max, config.omega_count)


def finite_terms(omega: np.ndarray) -> np.ndarray:
    terms: list[np.ndarray] = []
    for index in range(len(channel_labels())):
        offset = 0.08 * (index + 1)
        slope = (-1.0) ** index * (0.22 + 0.025 * index)
        curve = 0.12 * np.cos((index + 1) * omega)
        terms.append(offset + slope * omega + 0.18 * omega * omega + curve)
    return np.vstack(terms)


def amplitude_from_residue(
    omega: np.ndarray,
    residues: np.ndarray,
    pole_locations: np.ndarray,
    finite: np.ndarray,
) -> np.ndarray:
    amplitudes = np.zeros((len(residues), len(omega)), dtype=float)
    for index, residue in enumerate(residues):
        amplitudes[index] = residue / (omega - pole_locations[index]) + finite[index]
    return amplitudes


def build_samples(config: SoftConfig) -> tuple[np.ndarray, np.ndarray, list[SampleSpec]]:
    rng = np.random.default_rng(config.seed)
    omega = omega_grid(config)
    residues = true_residues()
    finite = finite_terms(omega)
    zeros = np.zeros(len(residues), dtype=float)
    target = amplitude_from_residue(omega, residues, zeros, finite)

    shuffled_by_omega = target.copy()
    for column in range(shuffled_by_omega.shape[1]):
        shuffled_by_omega[:, column] = shuffled_by_omega[rng.permutation(len(residues)), column]

    residue_scramble = residues[rng.permutation(len(residues))] * np.array([1.0, -1.0, 1.0, 1.0, -1.0, 1.0])
    residue_scramble *= np.linalg.norm(residues) / max(np.linalg.norm(residue_scramble), 1.0e-12)

    pole_shifts = np.array([-0.032, -0.025, -0.038, -0.021, -0.030, -0.027], dtype=float)
    factorization_break = residues + np.array([0.18, -0.11, 0.16, 0.10, -0.14, 0.12], dtype=float)
    factorization_break *= np.linalg.norm(residues) / max(np.linalg.norm(factorization_break), 1.0e-12)

    smooth_wrong = 0.65 * residues[:, None] / (omega[None, :] + 0.14) + finite * 1.15

    samples = [
        SampleSpec("soft_target", "target", "none", target),
        SampleSpec(
            "spectrum_preserving_shuffle",
            "negative_control",
            "per_frequency_channel_shuffle",
            shuffled_by_omega,
        ),
        SampleSpec(
            "residue_scramble_control",
            "negative_control",
            "wrong_residue_same_pole",
            amplitude_from_residue(omega, residue_scramble, zeros, finite),
        ),
        SampleSpec(
            "pole_location_drift_control",
            "negative_control",
            "shifted_soft_pole",
            amplitude_from_residue(omega, residues, pole_shifts, finite),
        ),
        SampleSpec(
            "factorization_break_control",
            "negative_control",
            "non_factorizing_residue",
            amplitude_from_residue(omega, factorization_break, zeros, finite),
        ),
        SampleSpec(
            "smooth_wrong_control",
            "negative_control",
            "smooth_no_soft_pole",
            smooth_wrong,
        ),
    ]
    return omega, residues, samples


def fit_soft_channel(
    omega: np.ndarray,
    amplitude: np.ndarray,
    soft_fit_max_omega: float,
) -> tuple[float, float, float]:
    mask = omega <= soft_fit_max_omega
    omega_fit = omega[mask]
    amplitude_fit = amplitude[mask]
    best_residue = 0.0
    best_pole = 0.0
    best_error = float("inf")
    for pole in np.linspace(-0.060, 0.001, 180):
        if np.min(np.abs(omega_fit - pole)) < 1.0e-5:
            continue
        design = np.column_stack([1.0 / (omega_fit - pole), np.ones_like(omega_fit), omega_fit])
        coefficients = np.linalg.lstsq(design, amplitude_fit, rcond=None)[0]
        residual = amplitude_fit - design @ coefficients
        error = float(np.mean(residual * residual))
        if error < best_error:
            best_error = error
            best_residue = float(coefficients[0])
            best_pole = float(pole)
    return best_residue, best_pole, best_error


def residue_recoverability(fitted: np.ndarray, truth: np.ndarray) -> tuple[float, float]:
    error_norm = float(np.linalg.norm(fitted - truth) / max(np.linalg.norm(truth), 1.0e-12))
    return float(np.clip(1.0 - error_norm / 0.20, 0.0, 1.0)), error_norm


def pole_location_score(poles: np.ndarray) -> tuple[float, float]:
    mean_abs = float(np.mean(np.abs(poles)))
    return float(np.clip(1.0 - mean_abs / 0.018, 0.0, 1.0)), mean_abs


def factorization_score(residues: np.ndarray) -> float:
    products = np.array(
        [
            residues[0] * residues[5],
            residues[1] * residues[4],
            residues[2] * residues[3],
        ],
        dtype=float,
    )
    scale = max(float(np.mean(np.abs(products))), 1.0e-12)
    relative_spread = float(np.std(products) / scale)
    return float(np.clip(1.0 - relative_spread / 0.12, 0.0, 1.0))


def generic_smoothness_score(omega: np.ndarray, amplitudes: np.ndarray) -> float:
    scaled = omega[None, :] * amplitudes
    second = np.diff(scaled, n=2, axis=1)
    roughness = float(np.mean(np.abs(second)))
    amplitude_scale = float(np.mean(np.abs(scaled))) + 1.0e-12
    normalized = roughness / amplitude_scale
    return float(np.clip(1.0 - normalized / 0.020, 0.0, 1.0))


def evaluate_sample(
    config: SoftConfig,
    omega: np.ndarray,
    truth: np.ndarray,
    sample: SampleSpec,
) -> tuple[SampleSummary, list[ChannelDiagnostic]]:
    fitted_residues: list[float] = []
    fitted_poles: list[float] = []
    for channel_index in range(sample.amplitudes.shape[0]):
        residue, pole, _ = fit_soft_channel(
            omega,
            sample.amplitudes[channel_index],
            config.soft_fit_max_omega,
        )
        fitted_residues.append(residue)
        fitted_poles.append(pole)

    fitted = np.array(fitted_residues, dtype=float)
    poles = np.array(fitted_poles, dtype=float)
    residue_score, residue_error = residue_recoverability(fitted, truth)
    pole_score, mean_abs_pole = pole_location_score(poles)
    factor_score = factorization_score(fitted)
    smooth_score = generic_smoothness_score(omega, sample.amplitudes)
    soft_score = float(
        0.43 * residue_score + 0.25 * pole_score + 0.27 * factor_score + 0.05 * smooth_score
    )

    summary = SampleSummary(
        sample=sample.sample,
        sample_role=sample.sample_role,
        control_type=sample.control_type,
        residue_recoverability=residue_score,
        pole_location_score=pole_score,
        factorization_score=factor_score,
        generic_smoothness_score=smooth_score,
        soft_specificity_score=soft_score,
        mean_abs_pole=mean_abs_pole,
        residue_error_norm=residue_error,
    )

    diagnostics: list[ChannelDiagnostic] = []
    for label, true_value, fitted_value, fitted_pole in zip(channel_labels(), truth, fitted, poles):
        relative_error = abs(float(fitted_value - true_value)) / max(abs(float(true_value)), 1.0e-12)
        diagnostics.append(
            ChannelDiagnostic(
                sample=sample.sample,
                sample_role=sample.sample_role,
                control_type=sample.control_type,
                channel=label,
                true_residue=float(true_value),
                fitted_residue=float(fitted_value),
                residue_relative_error=float(relative_error),
                fitted_pole=float(fitted_pole),
                residue_recoverability=residue_score,
                pole_location_score=pole_score,
                factorization_score=factor_score,
                generic_smoothness_score=smooth_score,
                soft_specificity_score=soft_score,
            )
        )
    return summary, diagnostics


def run_probe(config: SoftConfig) -> tuple[list[SampleSummary], list[ChannelDiagnostic], dict[str, Any]]:
    omega, truth, samples = build_samples(config)
    summaries: list[SampleSummary] = []
    diagnostics: list[ChannelDiagnostic] = []
    for sample in samples:
        summary, rows = evaluate_sample(config, omega, truth, sample)
        summaries.append(summary)
        diagnostics.extend(rows)
    status = build_status(config, summaries)
    return summaries, diagnostics, status


def build_status(config: SoftConfig, summaries: list[SampleSummary]) -> dict[str, Any]:
    target = next(row for row in summaries if row.sample_role == "target")
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    control_scores = np.array([row.soft_specificity_score for row in controls], dtype=float)
    control_threshold = target.soft_specificity_score - config.pass_min_margin
    control_false_positive_rate = float(np.mean(control_scores >= control_threshold))
    control_max = max(controls, key=lambda row: row.soft_specificity_score)
    score_margin = float(target.soft_specificity_score - control_max.soft_specificity_score)
    separation_z = float(
        (target.soft_specificity_score - float(np.mean(control_scores)))
        / max(float(np.std(control_scores)), 1.0e-9)
    )
    target_passes = bool(
        target.soft_specificity_score >= config.pass_min_soft_score
        and target.residue_recoverability >= config.pass_min_residue_recoverability
        and target.pole_location_score >= config.pass_min_pole_score
        and target.factorization_score >= config.pass_min_factorization_score
    )
    controls_pass = bool(
        control_false_positive_rate <= config.pass_max_control_false_positive_rate
        and score_margin >= config.pass_min_margin
        and separation_z >= 1.5
    )
    if target_passes and controls_pass:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif target.soft_specificity_score >= 0.72 and separation_z >= 1.0:
        bridge_status = "MARGINAL"
        failure_mode = "SOFT_CONTROL_SEPARATION_INCOMPLETE"
    else:
        bridge_status = "OPEN"
        failure_mode = "SOFT_STRUCTURE_NOT_DISTINCT"

    return {
        "phase": "59",
        "bridge_name": "soft_structure_proxy_test",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS requires a known toy soft target to recover residues, pole location, "
            "and factorization while controls preserving smoothness or spectra fail."
        ),
        "config": asdict(config),
        "target_soft_specificity_score": float(target.soft_specificity_score),
        "target_residue_recoverability": float(target.residue_recoverability),
        "target_pole_location_score": float(target.pole_location_score),
        "target_factorization_score": float(target.factorization_score),
        "target_generic_smoothness_score": float(target.generic_smoothness_score),
        "control_max_score": float(control_max.soft_specificity_score),
        "control_max_sample": control_max.sample,
        "control_false_positive_threshold": float(control_threshold),
        "control_false_positive_rate": control_false_positive_rate,
        "score_margin_over_best_control": score_margin,
        "separation_z": separation_z,
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "60_claim_gated_physics_bridge_update",
        "non_claims": [
            "no real gravitational soft theorem recovery claim",
            "no celestial amplitude or CFT correlator claim",
            "no S-matrix reconstruction claim",
            "toy pole/residue/factorization data only",
        ],
    }


def write_csv(path: Path, rows: list[ChannelDiagnostic]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def write_report(path: Path, summaries: list[SampleSummary], status: dict[str, Any]) -> None:
    summary_rows = "\n".join(
        "| {sample} | {sample_role} | {control_type} | {soft_specificity_score:.6f} | {residue_recoverability:.6f} | {pole_location_score:.6f} | {factorization_score:.6f} | {generic_smoothness_score:.6f} |".format(
            **asdict(row)
        )
        for row in summaries
    )
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    leakage_rows = "\n".join(
        f"- `{row.sample}`: soft={row.soft_specificity_score:.6f}, "
        f"smooth={row.generic_smoothness_score:.6f}, factorization={row.factorization_score:.6f}"
        for row in sorted(controls, key=lambda item: item.soft_specificity_score, reverse=True)
    )
    report = f"""# Phase 59 Soft-Structure Proxy Test

## Scope

This file is generated by `experiments/physics_bridge/soft_structure_proxy_test/run_soft_structure_proxy_test.py`.

This is a toy amplitude-like benchmark with known soft-pole, residue, and factorization structure. It does not use real scattering amplitudes and does not claim recovery of gravitational soft theorems, celestial amplitudes, BMS charges, Virasoro structure, S-matrix data, or gravitational memory.

The narrow question is:

> Can soft-specific telemetry detect known pole/residue/factorization breakage better than generic smoothness metrics?

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- target soft score: `{status['target_soft_specificity_score']:.6f}`
- target residue recoverability: `{status['target_residue_recoverability']:.6f}`
- target pole-location score: `{status['target_pole_location_score']:.6f}`
- target factorization score: `{status['target_factorization_score']:.6f}`
- target generic smoothness score: `{status['target_generic_smoothness_score']:.6f}`
- control max sample: `{status['control_max_sample']}`
- control max score: `{status['control_max_score']:.6f}`
- control false-positive threshold: `{status['control_false_positive_threshold']:.6f}`
- control false-positive rate: `{status['control_false_positive_rate']:.6f}`
- separation z: `{status['separation_z']:.6f}`
- margin over best control: `{status['score_margin_over_best_control']:.6f}`

## Summary Table

| sample | role | control | soft_score | residue | pole | factorization | smoothness |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
{summary_rows}

## Control Leakage Read

{leakage_rows}

## Interpretation

`PASS` means the toy target is separable from controls on residue recovery, pole location, and factorization. It does not mean HAOS-IIP recovers physical soft theorems.

`MARGINAL` means the toy target is detected but at least one control remains too competitive.

`OPEN` means the telemetry is mostly measuring smoothness rather than soft structure.

## Boundary

The decisive distinction is between generic smoothness and soft-specific structure. Several controls remain smooth, but fail residue, pole, or factorization checks. This is the intended toy result and must not be promoted into a real amplitude claim.

## Authority

- CSV: `experiments/physics_bridge/soft_structure_proxy_test/soft_diagnostics.csv`
- JSON: `experiments/physics_bridge/soft_structure_proxy_test/bridge_status.json`
- Figures: `experiments/physics_bridge/soft_structure_proxy_test/figures/`
- Generator: `experiments/physics_bridge/soft_structure_proxy_test/run_soft_structure_proxy_test.py`
"""
    path.write_text(report, encoding="utf-8")


def write_plots(output_dir: Path, summaries: list[SampleSummary], diagnostics: list[ChannelDiagnostic]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    labels = [row.sample for row in summaries]
    x = np.arange(len(summaries))

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(labels, [row.soft_specificity_score for row in summaries], color=["tab:green"] + ["tab:orange"] * (len(summaries) - 1))
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("soft specificity score")
    ax.set_title("Phase 59 toy soft-structure target vs controls")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "control_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.20
    ax.bar(x - 1.5 * width, [row.residue_recoverability for row in summaries], width, label="residue")
    ax.bar(x - 0.5 * width, [row.pole_location_score for row in summaries], width, label="pole")
    ax.bar(x + 0.5 * width, [row.factorization_score for row in summaries], width, label="factorization")
    ax.bar(x + 1.5 * width, [row.generic_smoothness_score for row in summaries], width, label="smoothness")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("score")
    ax.set_title("Soft-specific metrics vs generic smoothness")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "metric_breakdown.png", dpi=180)
    plt.close(fig)

    target_rows = [row for row in diagnostics if row.sample == "soft_target"]
    best_control = max(
        [row for row in summaries if row.sample_role == "negative_control"],
        key=lambda row: row.soft_specificity_score,
    )
    control_rows = [row for row in diagnostics if row.sample == best_control.sample]
    channel_x = np.arange(len(channel_labels()))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(channel_x - width / 2, [row.true_residue for row in target_rows], width, label="true residue")
    ax.bar(channel_x + width / 2, [row.fitted_residue for row in target_rows], width, label="target fitted")
    ax.scatter(channel_x, [row.fitted_residue for row in control_rows], color="tab:red", label=best_control.sample)
    ax.set_xticks(channel_x, channel_labels())
    ax.set_ylabel("residue")
    ax.set_title("Residue recovery and best-control contrast")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "residue_recovery.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        [row.generic_smoothness_score for row in summaries],
        [row.soft_specificity_score for row in summaries],
        c=["tab:green" if row.sample_role == "target" else "tab:orange" for row in summaries],
    )
    for row in summaries:
        ax.annotate(row.sample, (row.generic_smoothness_score, row.soft_specificity_score), fontsize=7)
    ax.set_xlim(0.0, 1.05)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("generic smoothness")
    ax.set_ylabel("soft specificity")
    ax.set_title("Smoothness is not soft structure")
    fig.tight_layout()
    fig.savefig(figure_dir / "smoothness_vs_soft_specificity.png", dpi=180)
    plt.close(fig)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def main() -> None:
    args = parse_args()
    config = make_config(args)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries, diagnostics, status = run_probe(config)
    write_csv(output_dir / CSV_PATH.name, diagnostics)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, summaries, status)
    write_plots(output_dir, summaries, diagnostics)

    print("Phase 59 soft-structure proxy test generated.")
    print(f"CSV: {relative(output_dir / CSV_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"Figures: {relative(output_dir / 'figures')}")
    print(
        f"bridge_status: {status['bridge_status']} ({status['failure_mode']}); "
        f"target={status['target_soft_specificity_score']:.6f}, "
        f"control_max={status['control_max_score']:.6f}"
    )


if __name__ == "__main__":
    main()
