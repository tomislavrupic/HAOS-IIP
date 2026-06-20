#!/usr/bin/env python3
"""Monstrous moonshine arithmetic support diagnostic.

This sidecar records a small, auditable subset of Monstrous Moonshine:

- the 15 supersingular primes,
- the Monster order prime support and exponents,
- a few low j-coefficient decomposition witnesses,
- factor support of the embedded Monster irreducible dimensions,
- Gaussian-prime ramified / inert / split classes on the same prime support.

It is a diagnostic bridge from arithmetic invariants to HAOS-IIP-style
telemetry. It does not construct the Moonshine module, prove moonshine, or
establish physics.
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
REPO_ROOT = Path(__file__).resolve().parents[3]

PRECOMMITMENT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "source_manifest.json"
PRIME_TABLE_PATH = ROOT / "supersingular_prime_table.csv"
WITNESS_TABLE_PATH = ROOT / "moonshine_witnesses.csv"
COMPONENT_SCORES_PATH = ROOT / "component_scores.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "monstrous_moonshine_diagnostic_report.md"
RESULT_PATH = ROOT / "monstrous_moonshine_diagnostic_result.json"

SUPERSINGULAR_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)

MONSTER_ORDER_EXPONENTS = {
    2: 46,
    3: 20,
    5: 9,
    7: 6,
    11: 2,
    13: 3,
    17: 1,
    19: 1,
    23: 1,
    29: 1,
    31: 1,
    41: 1,
    47: 1,
    59: 1,
    71: 1,
}

IRREP_DIMS = {
    "trivial": 1,
    "dim_196883": 196883,
    "dim_21296876": 21296876,
    "dim_842609326": 842609326,
    "dim_18538750076": 18538750076,
    "dim_19360062527": 19360062527,
}

J_DECOMPOSITION_WITNESSES = (
    {
        "name": "j_q1",
        "coefficient": 196884,
        "terms": {"trivial": 1, "dim_196883": 1},
    },
    {
        "name": "j_q2",
        "coefficient": 21493760,
        "terms": {"trivial": 1, "dim_196883": 1, "dim_21296876": 1},
    },
    {
        "name": "j_q3",
        "coefficient": 864299970,
        "terms": {"trivial": 2, "dim_196883": 2, "dim_21296876": 1, "dim_842609326": 1},
    },
)

CONDITIONS = (
    "known_positive",
    "nonsupersingular_prime_control",
    "exponent_shuffle_control",
    "decomposition_break_control",
    "gaussian_class_swap_control",
)

COMPONENT_SCORE_FIELDS = ["condition", "component", "metric", "value"]
CONTROL_RESULT_FIELDS = [
    "condition",
    "expected_response",
    "support_distance",
    "exponent_distance",
    "decomposition_distance",
    "gaussian_class_distance",
    "observed_status",
]
PRIME_TABLE_FIELDS = ["prime", "monster_order_exponent", "gaussian_class"]
WITNESS_TABLE_FIELDS = ["name", "coefficient", "computed_sum", "residual", "terms"]


@dataclass(frozen=True)
class MoonshineDiagnosticConfig:
    version: str = "monstrous-moonshine-diagnostic-v0.1"
    invariance_max: float = 1.0e-12
    support_control_min: float = 0.10
    exponent_control_min: float = 0.20
    decomposition_control_min: float = 0.20
    gaussian_class_control_min: float = 0.20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Monstrous Moonshine arithmetic support diagnostic.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def stable_unit(*parts: Any) -> float:
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()
    return int(digest[:16], 16) / float(16**16 - 1)


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


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


def factor_integer(value: int) -> dict[int, int]:
    if value < 1:
        raise ValueError("factor_integer expects a positive integer")
    remaining = value
    factors: dict[int, int] = {}
    divisor = 2
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            remaining //= divisor
        divisor += 1 if divisor == 2 else 2
    if remaining > 1:
        factors[remaining] = factors.get(remaining, 0) + 1
    return factors


def gaussian_prime_class(prime: int) -> str:
    if prime == 2:
        return "ramified"
    if prime % 4 == 1:
        return "split"
    if prime % 4 == 3:
        return "inert"
    return "unknown"


def condition_payload(condition: str) -> dict[str, Any]:
    supersingular = list(SUPERSINGULAR_PRIMES)
    exponents = dict(MONSTER_ORDER_EXPONENTS)
    witnesses = json.loads(json.dumps(J_DECOMPOSITION_WITNESSES))
    class_override: dict[int, str] = {}

    if condition == "known_positive":
        pass
    elif condition == "nonsupersingular_prime_control":
        supersingular = [prime for prime in supersingular if prime != 71] + [73]
        exponents[73] = exponents.pop(71)
    elif condition == "exponent_shuffle_control":
        primes = list(SUPERSINGULAR_PRIMES)
        shuffled = sorted(MONSTER_ORDER_EXPONENTS.values(), key=lambda exponent: stable_unit("exponent_shuffle", exponent))
        exponents = dict(zip(primes, shuffled))
    elif condition == "decomposition_break_control":
        witnesses[2]["coefficient"] = witnesses[2]["coefficient"] + 1
    elif condition == "gaussian_class_swap_control":
        class_override[3] = "split"
        class_override[5] = "inert"
    else:
        raise ValueError(f"unknown condition {condition}")

    return {
        "supersingular_primes": supersingular,
        "monster_order_exponents": exponents,
        "witnesses": witnesses,
        "class_override": class_override,
    }


def support_component(payload: dict[str, Any]) -> dict[str, Any]:
    support = sorted(payload["monster_order_exponents"])
    supersingular = sorted(payload["supersingular_primes"])
    support_set = set(support)
    supersingular_set = set(supersingular)
    intersection = support_set & supersingular_set
    union = support_set | supersingular_set
    return {
        "monster_prime_support": support,
        "supersingular_primes": supersingular,
        "support_jaccard": len(intersection) / len(union) if union else 1.0,
        "missing_from_supersingular": sorted(support_set - supersingular_set),
        "extra_supersingular_candidates": sorted(supersingular_set - support_set),
        "support_count": len(support),
    }


def exponent_component(payload: dict[str, Any]) -> dict[str, Any]:
    exponents = payload["monster_order_exponents"]
    vector = [int(exponents.get(prime, 0)) for prime in SUPERSINGULAR_PRIMES]
    return {
        "exponent_vector": vector,
        "exponent_sum": sum(vector),
        "nonzero_exponent_count": sum(1 for value in vector if value > 0),
        "log10_order": sum(value * math.log10(prime) for prime, value in exponents.items()),
    }


def witness_sum(witness: dict[str, Any]) -> int:
    return int(sum(int(multiplier) * IRREP_DIMS[name] for name, multiplier in witness["terms"].items()))


def decomposition_component(payload: dict[str, Any]) -> dict[str, Any]:
    witnesses = payload["witnesses"]
    rows = []
    residuals = []
    exact_count = 0
    for witness in witnesses:
        computed = witness_sum(witness)
        residual = int(witness["coefficient"]) - computed
        residuals.append(abs(residual))
        if residual == 0:
            exact_count += 1
        rows.append(
            {
                "name": witness["name"],
                "coefficient": int(witness["coefficient"]),
                "computed_sum": computed,
                "residual": residual,
                "terms": witness["terms"],
            }
        )
    support = set(payload["supersingular_primes"])
    irrep_factor_rows = []
    support_good = 0
    nontrivial_dims = {key: value for key, value in IRREP_DIMS.items() if value > 1}
    for name, value in nontrivial_dims.items():
        factors = factor_integer(value)
        factor_support = set(factors)
        is_supported = factor_support <= support
        if is_supported:
            support_good += 1
        irrep_factor_rows.append(
            {
                "name": name,
                "dimension": value,
                "factorization": {str(prime): exponent for prime, exponent in factors.items()},
                "factor_support_subset": is_supported,
            }
        )
    return {
        "witnesses": rows,
        "residual_abs_sum": int(sum(residuals)),
        "exact_witness_count": exact_count,
        "witness_count": len(witnesses),
        "irrep_factor_support_coverage": support_good / len(nontrivial_dims) if nontrivial_dims else 0.0,
        "irrep_factor_rows": irrep_factor_rows,
    }


def gaussian_class_component(payload: dict[str, Any]) -> dict[str, Any]:
    overrides = payload["class_override"]
    classes = [overrides.get(prime, gaussian_prime_class(prime)) for prime in SUPERSINGULAR_PRIMES]
    class_codes = [{"ramified": 0, "split": 1, "inert": -1, "unknown": 99}[label] for label in classes]
    return {
        "class_vector": class_codes,
        "ramified_count": classes.count("ramified"),
        "split_count": classes.count("split"),
        "inert_count": classes.count("inert"),
        "unknown_count": classes.count("unknown"),
    }


def run_condition(condition: str) -> dict[str, Any]:
    payload = condition_payload(condition)
    return {
        "payload": payload,
        "components": {
            "support": support_component(payload),
            "exponent": exponent_component(payload),
            "decomposition": decomposition_component(payload),
            "gaussian_class": gaussian_class_component(payload),
        },
    }


def vector_distance(left: Iterable[float], right: Iterable[float]) -> float:
    a = np.array(list(left), dtype=float)
    b = np.array(list(right), dtype=float)
    denom = max(float(np.linalg.norm(a)), 1.0e-12)
    return float(np.linalg.norm(a - b) / denom)


def jaccard_distance(left: Iterable[int], right: Iterable[int]) -> float:
    left_set = set(left)
    right_set = set(right)
    union = left_set | right_set
    if not union:
        return 0.0
    return float(1.0 - len(left_set & right_set) / len(union))


def support_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    monster_drift = jaccard_distance(reference["monster_prime_support"], observed["monster_prime_support"])
    supersingular_drift = jaccard_distance(reference["supersingular_primes"], observed["supersingular_primes"])
    internal_alignment_drift = abs(float(reference["support_jaccard"]) - float(observed["support_jaccard"]))
    return float(monster_drift + supersingular_drift + internal_alignment_drift)


def exponent_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    return vector_distance(reference["exponent_vector"], observed["exponent_vector"])


def decomposition_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    fields = ["residual_abs_sum", "exact_witness_count", "irrep_factor_support_coverage"]
    return vector_distance([reference[field] for field in fields], [observed[field] for field in fields])


def gaussian_class_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    fields = ["ramified_count", "split_count", "inert_count", "unknown_count"]
    return vector_distance(reference["class_vector"] + [reference[field] for field in fields], observed["class_vector"] + [observed[field] for field in fields])


def expected_response(condition: str) -> str:
    return {
        "known_positive": "all component distances are zero",
        "nonsupersingular_prime_control": "prime support and irrep support coverage should degrade",
        "exponent_shuffle_control": "prime support remains but Monster-order exponent address moves",
        "decomposition_break_control": "j-coefficient decomposition witness should fail exactly",
        "gaussian_class_swap_control": "Gaussian ramified/inert/split class relation should move",
    }[condition]


def observed_status(condition: str, distances: dict[str, float], config: MoonshineDiagnosticConfig) -> str:
    if condition == "known_positive":
        return "PASS" if max(distances.values()) <= config.invariance_max else "FAIL_SELF_RECOVERY"
    if condition == "nonsupersingular_prime_control":
        return "PASS" if distances["support"] >= config.support_control_min else "FAIL_SUPPORT_CONTROL"
    if condition == "exponent_shuffle_control":
        return "PASS" if distances["exponent"] >= config.exponent_control_min else "FAIL_EXPONENT_CONTROL"
    if condition == "decomposition_break_control":
        return "PASS" if distances["decomposition"] >= config.decomposition_control_min else "FAIL_DECOMPOSITION_CONTROL"
    if condition == "gaussian_class_swap_control":
        return "PASS" if distances["gaussian_class"] >= config.gaussian_class_control_min else "FAIL_GAUSSIAN_CLASS_CONTROL"
    raise ValueError(f"unknown condition {condition}")


def component_rows(condition: str, components: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for component, metrics in components.items():
        for metric, value in metrics.items():
            if metric in {"witnesses", "irrep_factor_rows"}:
                continue
            if isinstance(value, list):
                for index, item in enumerate(value):
                    rows.append({"condition": condition, "component": component, "metric": f"{metric}_{index}", "value": str(item)})
            else:
                rows.append({"condition": condition, "component": component, "metric": metric, "value": str(value)})
    return rows


def prime_rows(components: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    exponents = components["exponent"]["exponent_vector"]
    classes = components["gaussian_class"]["class_vector"]
    label_by_code = {0: "ramified", 1: "split", -1: "inert", 99: "unknown"}
    return [
        {
            "prime": prime,
            "monster_order_exponent": exponent,
            "gaussian_class": label_by_code[class_code],
        }
        for prime, exponent, class_code in zip(SUPERSINGULAR_PRIMES, exponents, classes)
    ]


def witness_rows(components: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "name": row["name"],
            "coefficient": row["coefficient"],
            "computed_sum": row["computed_sum"],
            "residual": row["residual"],
            "terms": json.dumps(row["terms"], sort_keys=True),
        }
        for row in components["decomposition"]["witnesses"]
    ]


def source_manifest() -> dict[str, Any]:
    return {
        "source": "Embedded arithmetic constants for a bounded moonshine diagnostic.",
        "code": repo_rel(ROOT / "run_monstrous_moonshine_diagnostic.py"),
        "mathematical_objects": [
            "15 supersingular primes",
            "Monster group order prime factorization",
            "low j-coefficient decomposition witnesses",
            "small embedded Monster irreducible dimensions",
            "Gaussian-prime residue classes over the same prime support",
        ],
        "external_data": "none",
        "voa_construction": "not included",
        "monster_irrep_catalog": "not complete; witness subset only",
    }


def precommitment_contract(config: MoonshineDiagnosticConfig) -> dict[str, Any]:
    return {
        "name": "Monstrous Moonshine Arithmetic Diagnostic v0.1",
        "status": "PRECOMMITTED_ARITHMETIC_SUPPORT_DIAGNOSTIC",
        "purpose": "Test a small auditable moonshine support pattern as HAOS-IIP-style arithmetic telemetry.",
        "source_artifacts": source_manifest(),
        "conditions": {condition: expected_response(condition) for condition in CONDITIONS},
        "thresholds": asdict(config),
        "components": [
            "supersingular_prime_support",
            "monster_order_exponent_address",
            "j_coefficient_decomposition_witnesses",
            "gaussian_prime_class_bridge",
        ],
        "claim_ceiling": "ARITHMETIC_SUPPORT_DIAGNOSTIC",
        "non_claims": [
            "not a proof of Monstrous Moonshine",
            "not a Moonshine module construction",
            "not a complete Monster irrep catalog",
            "not a physical bridge",
            "not a continuum limit",
            "not a quantum, gravity, or field-theory derivation",
            "not empirical validation",
            "not a change to frozen HAOS-IIP phases",
        ],
    }


def run_diagnostic(config: MoonshineDiagnosticConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    reference = run_condition("known_positive")
    score_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    condition_results: dict[str, Any] = {}

    for condition in CONDITIONS:
        observed = run_condition(condition)
        components = observed["components"]
        distances = {
            "support": support_distance(reference["components"]["support"], components["support"]),
            "exponent": exponent_distance(reference["components"]["exponent"], components["exponent"]),
            "decomposition": decomposition_distance(reference["components"]["decomposition"], components["decomposition"]),
            "gaussian_class": gaussian_class_distance(reference["components"]["gaussian_class"], components["gaussian_class"]),
        }
        status = observed_status(condition, distances, config)
        condition_results[condition] = {
            "components": components,
            "distances": distances,
            "expected_response": expected_response(condition),
            "observed_status": status,
        }
        score_rows.extend(component_rows(condition, components))
        control_rows.append(
            {
                "condition": condition,
                "expected_response": expected_response(condition),
                "support_distance": f"{distances['support']:.12g}",
                "exponent_distance": f"{distances['exponent']:.12g}",
                "decomposition_distance": f"{distances['decomposition']:.12g}",
                "gaussian_class_distance": f"{distances['gaussian_class']:.12g}",
                "observed_status": status,
            }
        )

    controls_pass = all(result["observed_status"] == "PASS" for result in condition_results.values())
    labels = [
        "MOONSHINE_ARITHMETIC_SUPPORT_BUILT",
        "SUPERSINGULAR_PRIME_SUPPORT_AVAILABLE",
        "MONSTER_ORDER_EXPONENT_ADDRESS_AVAILABLE",
        "J_DECOMPOSITION_WITNESSES_AVAILABLE",
        "GAUSSIAN_CLASS_BRIDGE_AVAILABLE",
        "MOONSHINE_PROOF_NOT_ESTABLISHED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_ARITHMETIC_DIAGNOSTIC",
    ]
    labels.append("COMPONENT_CONTROLS_PASS" if controls_pass else "COMPONENT_CONTROLS_PARTIAL")

    result_payload = {
        "bridge_id": "monstrous_moonshine_diagnostic_v0_1",
        "status": "PASS" if controls_pass else "MIXED_OPEN",
        "classification": "ARITHMETIC_SUPPORT_DIAGNOSTIC",
        "labels": labels,
        "condition_results": condition_results,
        "claim_ceiling": "No moonshine proof, VOA construction, complete irrep catalog, or physics claim.",
    }
    result_payload["result_hash"] = stable_hash("moonshine_diag", result_payload)

    reference_components = reference["components"]
    write_json(output_dir / PRECOMMITMENT_PATH.name, precommitment_contract(config))
    write_json(output_dir / SOURCE_MANIFEST_PATH.name, source_manifest())
    write_csv(output_dir / PRIME_TABLE_PATH.name, prime_rows(reference_components), PRIME_TABLE_FIELDS)
    write_csv(output_dir / WITNESS_TABLE_PATH.name, witness_rows(reference_components), WITNESS_TABLE_FIELDS)
    write_csv(output_dir / COMPONENT_SCORES_PATH.name, score_rows, COMPONENT_SCORE_FIELDS)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_RESULT_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result_payload)
    write_report(output_dir / REPORT_PATH.name, result_payload, control_rows)
    return result_payload


def write_report(path: Path, result: dict[str, Any], control_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Monstrous Moonshine Arithmetic Diagnostic Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Control Results",
        "",
        "| Condition | Support | Exponent | Decomposition | Gaussian class | Status |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in control_rows:
        lines.append(
            "| {condition} | {support_distance} | {exponent_distance} | {decomposition_distance} | {gaussian_class_distance} | `{observed_status}` |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This sidecar records a small arithmetic support diagnostic for moonshine-like telemetry.",
            "It does not prove Monstrous Moonshine, construct the Moonshine module, include the complete Monster irrep catalog, or establish any physical bridge.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_diagnostic(MoonshineDiagnosticConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
