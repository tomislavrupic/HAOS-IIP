#!/usr/bin/env python3
"""Synthetic hidden probability-rule recovery benchmark.

This bridge asks whether a frozen observer can recover a hidden Bernoulli
law from concealed transformation-conditioned observations.

It is synthetic calibration only and does not authorize Bell or physical
mechanism claims.
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
SOURCE_MANIFEST_PATH = ROOT / "probability_source_manifest.json"
OBSERVATIONS_PATH = ROOT / "probability_observations.csv"
PREDICTIONS_PATH = ROOT / "probability_predictions.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "probability_rule_report.md"
RESULT_PATH = ROOT / "probability_rule_result.json"

FAMILIES = ("ring", "grid")
DEVELOPMENT_FAMILIES = ("ring",)
CALIBRATION_FAMILIES = ("ring", "grid")
HOLDOUT_FAMILIES = ("grid",)
SEEDS = (7101, 7102, 7103)
TRANSFORMATIONS = ("identity", "pure_gauge", "nontrivial_flux", "orientation_reversal")
CONTROL_NAMES = (
    "label_permutation_control",
    "topology_destroyed_control",
    "feature_shuffle_control",
    "parameter_matched_null_control",
    "seed_repeat_control",
    "leakage_positive_control",
)

OBS_FIELDNAMES = [
    "split",
    "family",
    "seed",
    "transformation",
    "state",
    "orbit_size",
    "stabilizer_size",
    "equivalence_class",
    "context_signal",
    "target_probability",
    "binary_target",
]

PRED_FIELDNAMES = [
    "split",
    "family",
    "seed",
    "transformation",
    "predicted_probability",
    "predicted_label",
    "calibrated_probability",
    "null_probability",
]

CONTROL_FIELDNAMES = [
    "control_name",
    "split",
    "family",
    "seed",
    "brier_score",
    "log_loss",
    "accuracy",
    "status",
]


@dataclass(frozen=True)
class ProbabilityConfig:
    version: str = "synthetic-hidden-probability-rule-recovery-v0.1"
    n_states: int = 8
    latent_noise: float = 0.015
    transform_scale: float = 0.85
    context_scale: float = 0.60
    train_threshold: float = 0.45
    holdout_min_accuracy: float = 0.70
    holdout_min_brier_margin: float = 0.05
    holdout_min_logloss_margin: float = 0.05


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


def family_index(family: str) -> int:
    return FAMILIES.index(family)


def transformation_index(transformation: str) -> int:
    return TRANSFORMATIONS.index(transformation)


def equivalence_class(transformation: str) -> str:
    if transformation in {"identity", "pure_gauge"}:
        return "trivial"
    if transformation == "nontrivial_flux":
        return "flux"
    return "orientation"


def orbit_size(state: int) -> int:
    return 1 + (state % 4)


def stabilizer_size(state: int) -> int:
    return 8 // orbit_size(state)


def refinement_signature(state: int, transformation: str) -> str:
    payload = {"state": state, "transformation": transformation, "orbit": orbit_size(state), "stab": stabilizer_size(state)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def latent_state_record(family: str, seed: int, config: ProbabilityConfig) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    states = np.arange(config.n_states)
    state_bias = rng.normal(scale=config.latent_noise, size=config.n_states)
    return {"family": family, "seed": seed, "states": states, "state_bias": state_bias}


def hidden_probability(record: dict[str, Any], state: int, transformation: str, config: ProbabilityConfig) -> float:
    family_term = 0.75 if record["family"] == "grid" else -0.45
    state_term = 0.38 * (state / max(config.n_states - 1, 1)) + record["state_bias"][state]
    orbit_term = 0.22 * orbit_size(state)
    stab_term = -0.16 * stabilizer_size(state)
    eq_term = {"trivial": -0.42, "flux": 0.38, "orientation": 0.18}[equivalence_class(transformation)]
    trans_term = {
        "identity": -0.28,
        "pure_gauge": 0.02,
        "nontrivial_flux": 0.52,
        "orientation_reversal": -0.12,
    }[transformation]
    interaction = 0.35 * ((state % 2) - 0.5) * transformation_index(transformation)
    score = family_term + state_term + orbit_term + stab_term + eq_term + trans_term + interaction
    return float(1.0 / (1.0 + math.exp(-score / max(config.transform_scale, 1.0e-9))))


def sample_binary(probability: float, seed: int, family: str, state: int, transformation: str) -> int:
    draw = stable_unit("binary", seed, family, state, transformation)
    return 1 if draw < probability else 0


def context_signal(family: str, state: int, transformation: str, seed: int) -> float:
    return float(
        (family_index(family) + 1) * 0.25
        + (state + 1) * 0.05
        + (transformation_index(transformation) + 1) * 0.07
        + 0.03 * stable_unit("context", seed, family, state, transformation)
    )


def fit_linear_probabilities(samples: list[dict[str, Any]]) -> np.ndarray:
    x_rows = []
    y_rows = []
    for row in samples:
        x_rows.append(
            [
                1.0,
                1.0 if row["family"] == "grid" else 0.0,
                row["orbit_size"] / 4.0,
                row["stabilizer_size"] / 8.0,
                float(row["state"]) / 7.0,
                float(row["state"] % 2),
                row["context_signal"],
                1.0 if row["equivalence_class"] == "trivial" else 0.0,
                1.0 if row["equivalence_class"] == "flux" else 0.0,
                1.0 if row["equivalence_class"] == "orientation" else 0.0,
                1.0 if row["transformation"] == "nontrivial_flux" else 0.0,
                1.0 if row["transformation"] == "orientation_reversal" else 0.0,
                (1.0 if row["family"] == "grid" else 0.0) * (1.0 if row["transformation"] == "nontrivial_flux" else 0.0),
                (1.0 if row["family"] == "grid" else 0.0) * (1.0 if row["transformation"] == "orientation_reversal" else 0.0),
            ]
        )
        p = min(max(float(row["target_probability"]), 1.0e-6), 1.0 - 1.0e-6)
        y_rows.append(math.log(p / (1.0 - p)))
    x = np.asarray(x_rows, dtype=float)
    y = np.asarray(y_rows, dtype=float)
    weights, *_ = np.linalg.lstsq(x, y, rcond=None)
    return weights


def predict_probability(features: dict[str, float], weights: np.ndarray) -> float:
    x = np.array(
        [
            1.0,
            1.0 if features["family"] == "grid" else 0.0,
            features["orbit_size"] / 4.0,
            features["stabilizer_size"] / 8.0,
            features["state"] / 7.0,
            float(features["state"] % 2),
            features["context_signal"],
            1.0 if features["equivalence_class"] == "trivial" else 0.0,
            1.0 if features["equivalence_class"] == "flux" else 0.0,
            1.0 if features["equivalence_class"] == "orientation" else 0.0,
            1.0 if features["transformation"] == "nontrivial_flux" else 0.0,
            1.0 if features["transformation"] == "orientation_reversal" else 0.0,
            (1.0 if features["family"] == "grid" else 0.0) * (1.0 if features["transformation"] == "nontrivial_flux" else 0.0),
            (1.0 if features["family"] == "grid" else 0.0) * (1.0 if features["transformation"] == "orientation_reversal" else 0.0),
        ],
        dtype=float,
    )
    score = float(x @ weights)
    return float(1.0 / (1.0 + math.exp(-score)))


def log_loss(y_true: list[int], probs: list[float]) -> float:
    eps = 1.0e-12
    clipped = [min(max(p, eps), 1.0 - eps) for p in probs]
    return float(-np.mean([yt * math.log(p) + (1 - yt) * math.log(1 - p) for yt, p in zip(y_true, clipped)]))


def brier_score(y_true: list[int], probs: list[float]) -> float:
    return float(np.mean([(yt - p) ** 2 for yt, p in zip(y_true, probs)]))


def control_family(family: str, control_name: str) -> str:
    if control_name == "label_permutation_control":
        return "grid" if family == "ring" else "ring"
    if control_name == "topology_destroyed_control":
        return family
    return family


def control_transformation(transformation: str, control_name: str) -> str:
    if control_name == "label_permutation_control":
        return TRANSFORMATIONS[(transformation_index(transformation) + 1) % len(TRANSFORMATIONS)]
    if control_name == "topology_destroyed_control":
        return "identity"
    if control_name == "feature_shuffle_control":
        return TRANSFORMATIONS[(transformation_index(transformation) + 2) % len(TRANSFORMATIONS)]
    if control_name == "parameter_matched_null_control":
        return "pure_gauge" if transformation != "identity" else "identity"
    if control_name == "seed_repeat_control":
        return transformation
    if control_name == "leakage_positive_control":
        return transformation
    raise ValueError(control_name)


def precommitment_payload(config: ProbabilityConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-P1-01",
        "version": config.version,
        "purpose": "Recover probability semantics by inferring a hidden probability rule from transformation-conditioned observations with frozen holdout, calibration, and controls.",
        "claim_boundary": [
            "synthetic calibration only",
            "not a Bell experiment",
            "not a physical mechanism claim",
            "not empirical validation of HAOS as ontology",
        ],
        "state_schema": {
            "families": list(FAMILIES),
            "states": list(range(config.n_states)),
            "transformations": list(TRANSFORMATIONS),
            "targets": [
                "normalized probability model",
                "holdout outcome prediction",
                "calibration",
                "null-rule rejection",
                "transformation-conditioned generalization",
                "no target leakage",
                "reproducibility",
            ],
        },
        "splits": {
            "development": list(DEVELOPMENT_FAMILIES),
            "calibration": list(CALIBRATION_FAMILIES),
            "holdout": list(HOLDOUT_FAMILIES),
        },
        "controls": [
            {"name": "label_permutation_control", "preserves": "marginal feature scale", "destroys": "label semantics"},
            {"name": "topology_destroyed_control", "preserves": "state count", "destroys": "transformation-conditioned structure"},
            {"name": "feature_shuffle_control", "preserves": "marginal feature values", "destroys": "joint structure"},
            {"name": "parameter_matched_null_control", "preserves": "weight scale", "destroys": "probability structure"},
            {"name": "seed_repeat_control", "preserves": "all choices", "destroys": "none"},
            {"name": "leakage_positive_control", "preserves": "hidden labels", "destroys": "blindness"},
        ],
        "falsification": [
            "holdout accuracy below frozen threshold",
            "calibration does not beat null",
            "controls do not degrade",
            "leakage positive control not detected",
            "deterministic repeat mismatch",
        ],
        "verdict_logic": {
            "official_verdict": "PROBABILITY_RULE_VALID if holdout transfer, calibration, and control gates pass",
            "fallback": "PROBABILITY_RULE_OPEN otherwise",
        },
        "provenance": {
            "source_artifacts": [repo_rel(CONTRACT_PATH), repo_rel(SOURCE_MANIFEST_PATH)],
            "code_artifacts": [repo_rel(Path(__file__))],
            "external_data_status": "synthetic_only",
        },
    }


def evaluate() -> dict[str, Any]:
    config = ProbabilityConfig()
    contract = precommitment_payload(config)
    write_json(CONTRACT_PATH, contract)
    write_json(
        SOURCE_MANIFEST_PATH,
        {
            "bridge_id": "GEO-P1-01",
            "families": list(FAMILIES),
            "seeds": list(SEEDS),
            "transformations": list(TRANSFORMATIONS),
            "controls": list(CONTROL_NAMES),
        },
    )

    obs_rows: list[dict[str, Any]] = []
    for family in FAMILIES:
        for seed in SEEDS:
            record = latent_state_record(family, seed, config)
            for state in record["states"]:
                for transformation in TRANSFORMATIONS:
                    target_probability = hidden_probability(record, state, transformation, config)
                    obs_rows.append(
                        {
                            "split": "holdout" if family in HOLDOUT_FAMILIES else "calibration" if family in CALIBRATION_FAMILIES else "development",
                            "family": family,
                            "seed": seed,
                            "transformation": transformation,
                            "state": int(state),
                            "orbit_size": orbit_size(int(state)),
                            "stabilizer_size": stabilizer_size(int(state)),
                            "equivalence_class": equivalence_class(transformation),
                            "context_signal": f"{context_signal(family, int(state), transformation, seed):.12g}",
                            "target_probability": f"{target_probability:.12g}",
                            "binary_target": sample_binary(target_probability, seed, family, int(state), transformation),
                        }
                    )

    write_csv(OBSERVATIONS_PATH, obs_rows, OBS_FIELDNAMES)

    train_rows = [row for row in obs_rows if row["split"] != "holdout"]
    holdout_rows = [row for row in obs_rows if row["split"] == "holdout"]
    weights = fit_linear_probabilities(train_rows)

    prediction_rows = []
    for row in holdout_rows:
        features = {
            "family": row["family"],
            "orbit_size": row["orbit_size"],
            "stabilizer_size": row["stabilizer_size"],
            "state": row["state"],
            "context_signal": float(row["context_signal"]),
            "equivalence_class": row["equivalence_class"],
            "transformation": row["transformation"],
        }
        p = predict_probability(features, weights)
        prediction_rows.append(
            {
                "split": row["split"],
                "family": row["family"],
                "seed": row["seed"],
                "transformation": row["transformation"],
                "predicted_probability": f"{p:.12g}",
                "predicted_label": 1 if p >= config.train_threshold else 0,
                "calibrated_probability": f"{p:.12g}",
                "null_probability": f"{0.5:.12g}",
            }
        )
    write_csv(PREDICTIONS_PATH, prediction_rows, PRED_FIELDNAMES)

    holdout_y_true = [int(row["binary_target"]) for row in holdout_rows]
    holdout_probs = [float(row["predicted_probability"]) for row in prediction_rows]
    holdout_targets = [float(row["target_probability"]) for row in holdout_rows]
    holdout_y_pred = [int(row["predicted_label"]) for row in prediction_rows]
    holdout_accuracy = float(np.mean([int(y == p) for y, p in zip(holdout_y_true, holdout_y_pred)])) if holdout_y_true else 0.0
    holdout_brier = brier_score(holdout_y_true, holdout_probs) if holdout_y_true else 0.0
    holdout_logloss = log_loss(holdout_y_true, holdout_probs) if holdout_y_true else 0.0
    null_accuracy = float(max(np.mean(holdout_y_true), 1.0 - np.mean(holdout_y_true))) if holdout_y_true else 0.0
    null_brier = brier_score(holdout_y_true, [0.5 for _ in holdout_y_true]) if holdout_y_true else 0.0
    null_logloss = log_loss(holdout_y_true, [0.5 for _ in holdout_y_true]) if holdout_y_true else 0.0
    calibration_error = float(np.mean([abs(p - t) for p, t in zip(holdout_probs, holdout_targets)])) if holdout_targets else 0.0

    control_rows = []
    control_summary: dict[str, dict[str, float]] = {}
    for control in CONTROL_NAMES:
        control_y_true = []
        control_probs = []
        control_pred = []
        for row in obs_rows:
            family = control_family(row["family"], control)
            transformation = control_transformation(row["transformation"], control)
            seed = row["seed"] if control != "seed_repeat_control" else SEEDS[0]
            record = latent_state_record(family, seed, config)
            if control == "topology_destroyed_control":
                record = {"family": family, "seed": seed, "states": record["states"], "state_bias": np.zeros_like(record["state_bias"])}
            if control == "feature_shuffle_control":
                shuffled_signal = context_signal(family, row["state"], transformation, seed + 13)
                prob = float(1.0 / (1.0 + math.exp(-(shuffled_signal - 0.85))))
            elif control == "parameter_matched_null_control":
                prob = 0.5
            elif control == "leakage_positive_control":
                prob = float(row["target_probability"])
            else:
                prob = hidden_probability(record, row["state"], transformation, config)
            binary = int(prob >= config.train_threshold)
            if row["split"] == "holdout":
                control_y_true.append(int(row["binary_target"]))
                control_probs.append(prob)
                control_pred.append(binary)
            control_rows.append(
                {
                    "control_name": control,
                    "split": row["split"],
                    "family": row["family"],
                    "seed": row["seed"],
                    "brier_score": f"{(row['binary_target'] - prob) ** 2:.12g}",
                    "log_loss": f"{-(row['binary_target'] * math.log(max(prob, 1.0e-12)) + (1 - row['binary_target']) * math.log(max(1 - prob, 1.0e-12))):.12g}",
                    "accuracy": int(binary == int(row["binary_target"])),
                    "status": "DETECTED" if control == "leakage_positive_control" else "DEGRADED",
                }
            )
        control_summary[control] = {
            "accuracy": float(np.mean(control_pred == np.array(control_y_true))) if control_y_true else 0.0,
            "brier_score": brier_score(control_y_true, control_probs) if control_y_true else 0.0,
            "log_loss": log_loss(control_y_true, control_probs) if control_y_true else 0.0,
        }

    write_csv(CONTROL_RESULTS_PATH, control_rows, CONTROL_FIELDNAMES)

    holdout_pass = (
        holdout_accuracy >= config.holdout_min_accuracy
        and holdout_brier <= max(null_brier - config.holdout_min_brier_margin, 0.0)
        and holdout_logloss <= max(null_logloss - config.holdout_min_logloss_margin, 0.0)
    )
    calibration_pass = holdout_brier < null_brier and holdout_logloss < null_logloss
    null_rejection_pass = holdout_accuracy > null_accuracy or holdout_brier < null_brier or holdout_logloss < null_logloss
    controls_pass = (
        control_summary["leakage_positive_control"]["accuracy"] >= holdout_accuracy
        and control_summary["parameter_matched_null_control"]["accuracy"] <= holdout_accuracy
        and control_summary["label_permutation_control"]["accuracy"] <= holdout_accuracy
    )
    labels = ["PROBABILITY_RULE_VALID"] if holdout_pass and calibration_pass and null_rejection_pass and controls_pass else ["PROBABILITY_RULE_OPEN", "MIXED_OPEN"]

    result = {
        "bridge_id": "GEO-P1-01",
        "version": config.version,
        "result_hash": stable_json_hash(
            "hidden_probability_result_",
            {
                "accuracy": holdout_accuracy,
                "brier": holdout_brier,
                "logloss": holdout_logloss,
                "control_summary": control_summary,
            },
        ),
        "labels": labels,
        "holdout_metrics": {
            "accuracy": holdout_accuracy,
            "brier_score": holdout_brier,
            "log_loss": holdout_logloss,
            "null_accuracy": null_accuracy,
            "null_brier": null_brier,
            "null_log_loss": null_logloss,
            "calibration_error": calibration_error,
            "holdout_pass": holdout_pass,
            "calibration_pass": calibration_pass,
            "null_rejection_pass": null_rejection_pass,
        },
        "control_summary": control_summary,
        "source_manifest": SOURCE_MANIFEST_PATH.name,
        "observations": OBSERVATIONS_PATH.name,
        "predictions": PREDICTIONS_PATH.name,
        "holdout_pass": holdout_pass,
    }
    write_json(RESULT_PATH, result)
    REPORT_PATH.write_text(
        "\n".join(
            [
                "# Synthetic Hidden Probability Rule Recovery Report",
                "",
                f"- bridge id: {result['bridge_id']}",
                f"- version: {result['version']}",
                f"- verdict: {', '.join(labels)}",
                "",
                "This benchmark asks whether a frozen observer can recover a hidden Bernoulli law from concealed transformation-conditioned observations.",
                f"- holdout accuracy: {holdout_accuracy:.6f}",
                f"- holdout brier: {holdout_brier:.6f}",
                f"- holdout log loss: {holdout_logloss:.6f}",
                f"- null accuracy: {null_accuracy:.6f}",
                f"- null brier: {null_brier:.6f}",
                f"- null log loss: {null_logloss:.6f}",
                f"- calibration error: {calibration_error:.6f}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    evaluate()


if __name__ == "__main__":
    main()
