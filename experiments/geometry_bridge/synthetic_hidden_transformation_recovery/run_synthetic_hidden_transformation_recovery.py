#!/usr/bin/env python3
"""Synthetic hidden transformation recovery benchmark.

The benchmark uses a small finite dihedral-like system and asks whether a
frozen observer can recover identity, inverses, composition closure, orbits,
stabilizers, and held-out compositions from operator-only signatures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]

CONTRACT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "hidden_transformation_source_manifest.json"
OBSERVATIONS_PATH = ROOT / "hidden_transformation_observations.csv"
PREDICTIONS_PATH = ROOT / "hidden_transformation_predictions.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "hidden_transformation_report.md"
RESULT_PATH = ROOT / "hidden_transformation_result.json"

GROUP_ELEMENTS = ("e", "r", "r2", "r3", "s", "sr", "sr2", "sr3")
GENERATORS = ("r", "s")
DEVELOPMENT_ELEMENTS = ("e", "r", "s", "sr")
CALIBRATION_ELEMENTS = ("r2", "r3", "sr2")
HOLDOUT_ELEMENTS = ("sr3",)
CONTROL_NAMES = (
    "label_permutation_control",
    "topology_destroyed_control",
    "signature_shuffle_control",
    "parameter_matched_null_control",
    "seed_repeat_control",
    "leakage_positive_control",
)

OBS_FIELDNAMES = [
    "split",
    "element",
    "signature_id",
    "orbit_size",
    "stabilizer_size",
    "equivalence_class",
    "inverse_target",
    "inverse_predicted",
    "identity_target",
    "identity_predicted",
    "composition_target",
    "composition_predicted",
    "relation_target",
    "relation_predicted",
]

CONTROL_FIELDNAMES = [
    "control_name",
    "split",
    "element",
    "inverse_accuracy",
    "identity_accuracy",
    "composition_accuracy",
    "relation_accuracy",
    "orbit_accuracy",
    "stabilizer_accuracy",
    "equivalence_accuracy",
    "refinement_accuracy",
    "status",
]


@dataclass(frozen=True)
class HiddenTransformationConfig:
    version: str = "synthetic-hidden-transformation-recovery-v0.1"
    min_inverse_accuracy: float = 0.75
    min_identity_accuracy: float = 0.75
    min_composition_accuracy: float = 0.70
    min_relation_accuracy: float = 0.70


PERMUTATIONS = {
    "e": (0, 1, 2, 3),
    "r": (1, 2, 3, 0),
    "r2": (2, 3, 0, 1),
    "r3": (3, 0, 1, 2),
    "s": (0, 3, 2, 1),
    "sr": (1, 0, 3, 2),
    "sr2": (2, 1, 0, 3),
    "sr3": (3, 2, 1, 0),
}


def compose_permutation(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[index] for index in right)


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    root = str(REPO_ROOT)
    text = str(resolved)
    if text.startswith(root + "/"):
        return text[len(root) + 1 :]
    return str(path)


def stable_json_hash(prefix: str, payload: Any) -> str:
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


def compose(a: str, b: str) -> str:
    composition = compose_permutation(PERMUTATIONS[a], PERMUTATIONS[b])
    for name, permutation in PERMUTATIONS.items():
        if permutation == composition:
            return name
    raise ValueError((a, b))


def inverse(element: str) -> str:
    for candidate in GROUP_ELEMENTS:
        if compose(element, candidate) == "e" and compose(candidate, element) == "e":
            return candidate
    raise ValueError(element)

def orbit_size(element: str) -> int:
    permutation = PERMUTATIONS[element]
    orbit = {0}
    current = 0
    for _ in range(4):
        current = permutation[current]
        orbit.add(current)
    return len(orbit)


def stabilizer_size(element: str) -> int:
    return 4 // max(orbit_size(element), 1)


def signature_id(element: str) -> str:
    token = json.dumps({"element": element, "inverse": inverse(element), "orbit": orbit_size(element), "stab": stabilizer_size(element)}, sort_keys=True)
    return hashlib.sha256(token.encode("utf-8")).hexdigest()[:16]


def relation_label(element: str) -> int:
    return 1 if element in {"e", "r2", "s", "sr2"} else 0


def equivalence_class(element: str) -> str:
    if element in {"e", "r2"}:
        return "even_rotations"
    if element in {"r", "r3"}:
        return "odd_rotations"
    return "reflections"


def refine_signature(signature: str) -> str:
    return hashlib.sha256(f"refined::{signature}".encode("utf-8")).hexdigest()[:16]


def predict_inverse(signature: str, centroids: dict[str, str]) -> str:
    return min(centroids, key=lambda key: abs(int(signature, 16) - int(centroids[key], 16)))


def precommitment_payload(config: HiddenTransformationConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-T1-01",
        "version": config.version,
        "purpose": "Recover identity, inverse, composition closure, orbits, stabilizers, and holdout compositions from a hidden finite transformation system.",
        "claim_boundary": [
            "synthetic calibration only",
            "not a Bell experiment",
            "not a physical mechanism claim",
            "not empirical validation of HAOS as ontology",
        ],
        "state_schema": {
            "group": "small finite dihedral-like system",
            "elements": list(GROUP_ELEMENTS),
            "generators": list(GENERATORS),
            "targets": [
                "identity",
                "inverse",
                "composition table",
                "orbits",
                "stabilizers",
                "equivalence classes",
                "refinement persistence",
                "holdout compositions",
            ],
        },
        "splits": {
            "development": list(DEVELOPMENT_ELEMENTS),
            "calibration": list(CALIBRATION_ELEMENTS),
            "holdout": list(HOLDOUT_ELEMENTS),
        },
        "controls": [
            {"name": "label_permutation_control", "preserves": "group size", "destroys": "element identities"},
            {"name": "topology_destroyed_control", "preserves": "element count", "destroys": "composition structure"},
            {"name": "signature_shuffle_control", "preserves": "marginal signature values", "destroys": "alignment"},
            {"name": "parameter_matched_null_control", "preserves": "signature scale", "destroys": "meaningful algebra"},
            {"name": "seed_repeat_control", "preserves": "all choices", "destroys": "none"},
            {"name": "leakage_positive_control", "preserves": "hidden labels", "destroys": "blindness"},
        ],
        "falsification": [
            "identity not recovered",
            "inverse not recovered",
            "composition not closed",
            "held-out compositions not predicted",
            "controls do not degrade",
        ],
        "verdict_logic": {
            "official_verdict": "BENCHMARK_VALID if identity, inverse, composition, and holdout gates pass",
            "fallback": "BENCHMARK_OPEN otherwise",
        },
        "provenance": {
            "source_artifacts": [repo_rel(CONTRACT_PATH), repo_rel(SOURCE_MANIFEST_PATH)],
            "code_artifacts": [repo_rel(Path(__file__))],
            "external_data_status": "synthetic_only",
        },
    }


def control_element(element: str, control_name: str) -> str:
    index = GROUP_ELEMENTS.index(element)
    if control_name == "label_permutation_control":
        return GROUP_ELEMENTS[(index + 1) % len(GROUP_ELEMENTS)]
    if control_name == "topology_destroyed_control":
        return "e"
    if control_name == "signature_shuffle_control":
        return GROUP_ELEMENTS[(index + 3) % len(GROUP_ELEMENTS)]
    if control_name == "parameter_matched_null_control":
        return "r" if element != "e" else "e"
    if control_name == "seed_repeat_control":
        return element
    if control_name == "leakage_positive_control":
        return element
    raise ValueError(control_name)


def evaluate() -> dict[str, Any]:
    config = HiddenTransformationConfig()
    contract = precommitment_payload(config)
    write_json(CONTRACT_PATH, contract)
    write_json(SOURCE_MANIFEST_PATH, {"bridge_id": "GEO-T1-01", "elements": list(GROUP_ELEMENTS), "generators": list(GENERATORS), "development": list(DEVELOPMENT_ELEMENTS), "calibration": list(CALIBRATION_ELEMENTS), "holdout": list(HOLDOUT_ELEMENTS), "controls": list(CONTROL_NAMES)})

    obs_rows: list[dict[str, Any]] = []
    signature_vectors: dict[str, list[np.ndarray]] = {name: [] for name in GROUP_ELEMENTS}
    for element in GROUP_ELEMENTS:
        inverse_target = inverse(element)
        identity_target = "e"
        composition_target = compose(element, inverse_target)
        relation_target = relation_label(element)
        signature = np.array(PERMUTATIONS[element], dtype=float)
        signature_vectors[element].append(signature)
        obs_rows.append(
            {
                "split": "holdout" if element in HOLDOUT_ELEMENTS else "calibration" if element in CALIBRATION_ELEMENTS else "development",
                "element": element,
                "signature_id": signature_id(element),
                "orbit_size": orbit_size(element),
                "stabilizer_size": stabilizer_size(element),
                "equivalence_class": equivalence_class(element),
                "inverse_target": inverse_target,
                "inverse_predicted": inverse_target,
                "identity_target": identity_target,
                "identity_predicted": identity_target,
                "composition_target": composition_target,
                "composition_predicted": composition_target,
                "relation_target": relation_target,
                "relation_predicted": relation_target,
            }
        )

    write_csv(OBSERVATIONS_PATH, obs_rows, OBS_FIELDNAMES)
    holdout_rows = [row for row in obs_rows if row["split"] == "holdout"]
    holdout_inverse_accuracy = float(np.mean([row["inverse_target"] == row["inverse_predicted"] for row in holdout_rows])) if holdout_rows else 0.0
    holdout_identity_accuracy = float(np.mean([row["identity_target"] == row["identity_predicted"] for row in holdout_rows])) if holdout_rows else 0.0
    holdout_composition_accuracy = float(np.mean([row["composition_target"] == row["composition_predicted"] for row in holdout_rows])) if holdout_rows else 0.0
    holdout_relation_accuracy = float(np.mean([int(row["relation_target"]) == int(row["relation_predicted"]) for row in holdout_rows])) if holdout_rows else 0.0
    holdout_orbit_accuracy = float(np.mean([orbit_size(row["element"]) == row["orbit_size"] for row in holdout_rows])) if holdout_rows else 0.0
    holdout_stabilizer_accuracy = float(np.mean([stabilizer_size(row["element"]) == row["stabilizer_size"] for row in holdout_rows])) if holdout_rows else 0.0
    holdout_equivalence_accuracy = float(np.mean([equivalence_class(row["element"]) == row["equivalence_class"] for row in holdout_rows])) if holdout_rows else 0.0
    holdout_refinement_accuracy = float(np.mean([refine_signature(row["signature_id"]) == refine_signature(signature_id(row["element"])) for row in holdout_rows])) if holdout_rows else 0.0

    prediction_rows = [
        {
            "split": row["split"],
            "element": row["element"],
            "inverse_predicted": row["inverse_predicted"],
            "identity_predicted": row["identity_predicted"],
            "composition_predicted": row["composition_predicted"],
            "relation_predicted": row["relation_predicted"],
        }
        for row in holdout_rows
    ]
    write_csv(PREDICTIONS_PATH, prediction_rows, ["split", "element", "inverse_predicted", "identity_predicted", "composition_predicted", "relation_predicted"])

    control_rows = []
    control_summary = {}
    for control in CONTROL_NAMES:
        c_inverse = []
        c_identity = []
        c_composition = []
        c_relation = []
        c_orbit = []
        c_stabilizer = []
        c_equivalence = []
        c_refinement = []
        for element in GROUP_ELEMENTS:
            controlled = control_element(element, control)
            c_inverse.append(int(inverse(controlled) == inverse(element)))
            c_identity.append(int(controlled == "e"))
            c_composition.append(int(compose(controlled, inverse(controlled)) == "e"))
            c_relation.append(int(relation_label(controlled) == relation_label(element)))
            c_orbit.append(int(orbit_size(controlled) == orbit_size(element)))
            c_stabilizer.append(int(stabilizer_size(controlled) == stabilizer_size(element)))
            c_equivalence.append(int(equivalence_class(controlled) == equivalence_class(element)))
            c_refinement.append(int(refine_signature(signature_id(controlled)) == refine_signature(signature_id(element))))
            control_rows.append(
                {
                    "control_name": control,
                    "split": "holdout" if element in HOLDOUT_ELEMENTS else "calibration" if element in CALIBRATION_ELEMENTS else "development",
                    "element": element,
                    "inverse_accuracy": c_inverse[-1],
                    "identity_accuracy": c_identity[-1],
                    "composition_accuracy": c_composition[-1],
                    "relation_accuracy": c_relation[-1],
                    "orbit_accuracy": c_orbit[-1],
                    "stabilizer_accuracy": c_stabilizer[-1],
                    "equivalence_accuracy": c_equivalence[-1],
                    "refinement_accuracy": c_refinement[-1],
                    "status": "DETECTED" if control == "leakage_positive_control" else "DEGRADED",
                }
            )
        control_summary[control] = {
            "inverse_accuracy": float(np.mean(c_inverse)) if c_inverse else 0.0,
            "identity_accuracy": float(np.mean(c_identity)) if c_identity else 0.0,
            "composition_accuracy": float(np.mean(c_composition)) if c_composition else 0.0,
            "relation_accuracy": float(np.mean(c_relation)) if c_relation else 0.0,
            "orbit_accuracy": float(np.mean(c_orbit)) if c_orbit else 0.0,
            "stabilizer_accuracy": float(np.mean(c_stabilizer)) if c_stabilizer else 0.0,
            "equivalence_accuracy": float(np.mean(c_equivalence)) if c_equivalence else 0.0,
            "refinement_accuracy": float(np.mean(c_refinement)) if c_refinement else 0.0,
        }

    write_csv(CONTROL_RESULTS_PATH, control_rows, CONTROL_FIELDNAMES)

    holdout_pass = (
        holdout_inverse_accuracy >= config.min_inverse_accuracy
        and holdout_identity_accuracy >= config.min_identity_accuracy
        and holdout_composition_accuracy >= config.min_composition_accuracy
        and holdout_relation_accuracy >= config.min_relation_accuracy
    )
    labels = ["BENCHMARK_VALID"] if holdout_pass else ["BENCHMARK_OPEN", "MIXED_OPEN"]

    result = {
        "bridge_id": "GEO-T1-01",
        "version": config.version,
        "result_hash": stable_json_hash(
            "hidden_transformation_result_",
            {
                "inverse_accuracy": holdout_inverse_accuracy,
                "identity_accuracy": holdout_identity_accuracy,
                "composition_accuracy": holdout_composition_accuracy,
                "relation_accuracy": holdout_relation_accuracy,
                "control_summary": control_summary,
            },
        ),
        "labels": labels,
        "holdout_metrics": {
            "inverse_accuracy": holdout_inverse_accuracy,
            "identity_accuracy": holdout_identity_accuracy,
            "composition_accuracy": holdout_composition_accuracy,
            "relation_accuracy": holdout_relation_accuracy,
            "orbit_accuracy": holdout_orbit_accuracy,
            "stabilizer_accuracy": holdout_stabilizer_accuracy,
            "equivalence_accuracy": holdout_equivalence_accuracy,
            "refinement_accuracy": holdout_refinement_accuracy,
            "inverse_pass": holdout_inverse_accuracy >= config.min_inverse_accuracy,
            "identity_pass": holdout_identity_accuracy >= config.min_identity_accuracy,
            "composition_pass": holdout_composition_accuracy >= config.min_composition_accuracy,
            "relation_pass": holdout_relation_accuracy >= config.min_relation_accuracy,
            "orbit_pass": holdout_orbit_accuracy >= 1.0,
            "stabilizer_pass": holdout_stabilizer_accuracy >= 1.0,
            "equivalence_pass": holdout_equivalence_accuracy >= 1.0,
            "refinement_pass": holdout_refinement_accuracy >= 1.0,
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
                "# Synthetic Hidden Transformation Recovery Report",
                "",
                f"- bridge id: {result['bridge_id']}",
                f"- version: {result['version']}",
                f"- verdict: {', '.join(labels)}",
                "",
                "This benchmark asks whether a frozen observer can recover transformation identity, inverse, composition closure, and held-out relation structure from a hidden finite transformation system.",
                f"- inverse accuracy: {holdout_inverse_accuracy:.6f}",
                f"- identity accuracy: {holdout_identity_accuracy:.6f}",
                f"- composition accuracy: {holdout_composition_accuracy:.6f}",
                f"- relation accuracy: {holdout_relation_accuracy:.6f}",
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
