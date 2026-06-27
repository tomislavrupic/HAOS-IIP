#!/usr/bin/env python3
"""Pinned source validation for the moonshine arithmetic diagnostic.

This script checks the embedded constants used by GEO-MM-01 against pinned
source claims and local arithmetic identities. It prevents silent drift in the
constants, but it does not prove Monstrous Moonshine or certify external Lean
claims.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_betti_component_count import base_nodes, representative_for_prime  # noqa: E402
from run_monstrous_moonshine_diagnostic import (  # noqa: E402
    IRREP_DIMS,
    J_DECOMPOSITION_WITNESSES,
    MONSTER_ORDER_EXPONENTS,
    SUPERSINGULAR_PRIMES,
    factor_integer,
    gaussian_prime_class,
    stable_hash,
    witness_sum,
    write_json,
)


PINNED_MANIFEST_PATH = ROOT / "source_manifest_pinned.json"
REPORT_PATH = ROOT / "source_validation_report.md"
RESULT_PATH = ROOT / "arithmetic_source_validation_result.json"

RETRIEVAL_DATE = "2026-06-20"

PINNED_SUPERSINGULAR_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)
PINNED_MONSTER_ORDER_EXPONENTS = {
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
PINNED_J_COEFFICIENTS = {
    "j_q1": 196884,
    "j_q2": 21493760,
    "j_q3": 864299970,
}
PINNED_IRREP_DIMS = {
    "trivial": 1,
    "dim_196883": 196883,
    "dim_21296876": 21296876,
    "dim_842609326": 842609326,
    "dim_18538750076": 18538750076,
    "dim_19360062527": 19360062527,
}


@dataclass(frozen=True)
class SourceValidationConfig:
    version: str = "arithmetic-source-validation-v0.1"
    retrieval_date: str = RETRIEVAL_DATE
    claim_ceiling: str = "PINNED_SOURCE_VALIDATION_ONLY"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate pinned arithmetic constants for GEO-MM-01.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def value_hash(claim: str, value: Any) -> str:
    return stable_hash("source_value", {"claim": claim, "value": value})


def source_entries() -> list[dict[str, Any]]:
    return [
        {
            "claim": "15 supersingular primes in Moonshine usage",
            "source": "https://en.wikipedia.org/wiki/Supersingular_prime_(moonshine_theory)",
            "retrieval_date": RETRIEVAL_DATE,
            "hash_or_version": value_hash("supersingular_primes", list(PINNED_SUPERSINGULAR_PRIMES)),
            "used_in": "SUPERSINGULAR_PRIMES; Betti relation graph node support",
        },
        {
            "claim": "Monster group order prime exponents",
            "source": "ATLAS of Finite Groups; https://en.wikipedia.org/wiki/Monster_group",
            "retrieval_date": RETRIEVAL_DATE,
            "hash_or_version": value_hash("monster_order_exponents", PINNED_MONSTER_ORDER_EXPONENTS),
            "used_in": "MONSTER_ORDER_EXPONENTS; arithmetic support diagnostic",
        },
        {
            "claim": "first j-function coefficients used by the witness subset",
            "source": "https://en.wikipedia.org/wiki/J-invariant; https://en.wikipedia.org/wiki/Monstrous_moonshine",
            "retrieval_date": RETRIEVAL_DATE,
            "hash_or_version": value_hash("j_coefficients", PINNED_J_COEFFICIENTS),
            "used_in": "J_DECOMPOSITION_WITNESSES",
        },
        {
            "claim": "Monster irrep dimensions used by the embedded witness subset",
            "source": "ATLAS of Finite Groups; https://en.wikipedia.org/wiki/Monstrous_moonshine",
            "retrieval_date": RETRIEVAL_DATE,
            "hash_or_version": value_hash("irrep_dimensions", PINNED_IRREP_DIMS),
            "used_in": "IRREP_DIMS; j-coefficient witness sums",
        },
        {
            "claim": "Gaussian-prime residue class representatives",
            "source": "elementary Gaussian integer splitting law for rational primes",
            "retrieval_date": RETRIEVAL_DATE,
            "hash_or_version": value_hash(
                "gaussian_representatives",
                [{"prime": node["prime"], "x": node["x"], "y": node["y"], "class": node["gaussian_class"]} for node in base_nodes()],
            ),
            "used_in": "betti_relation_graph_nodes.csv; Betti relation graph construction",
        },
    ]


def validate_supersingular_primes() -> dict[str, Any]:
    observed = tuple(SUPERSINGULAR_PRIMES)
    expected = PINNED_SUPERSINGULAR_PRIMES
    return {
        "claim": "15 supersingular primes",
        "status": "PASS" if observed == expected else "FAIL",
        "expected": list(expected),
        "observed": list(observed),
    }


def validate_monster_order_exponents() -> dict[str, Any]:
    observed = {int(prime): int(exponent) for prime, exponent in MONSTER_ORDER_EXPONENTS.items()}
    expected = PINNED_MONSTER_ORDER_EXPONENTS
    return {
        "claim": "Monster order prime exponents",
        "status": "PASS" if observed == expected else "FAIL",
        "expected": expected,
        "observed": observed,
        "exponent_sum": sum(observed.values()),
    }


def validate_j_witnesses() -> dict[str, Any]:
    rows = []
    all_pass = True
    for witness in J_DECOMPOSITION_WITNESSES:
        name = witness["name"]
        computed = witness_sum(witness)
        expected_coefficient = PINNED_J_COEFFICIENTS[name]
        residual = int(witness["coefficient"]) - computed
        status = "PASS" if int(witness["coefficient"]) == expected_coefficient and residual == 0 else "FAIL"
        if status != "PASS":
            all_pass = False
        rows.append(
            {
                "name": name,
                "expected_coefficient": expected_coefficient,
                "observed_coefficient": int(witness["coefficient"]),
                "computed_sum": computed,
                "residual": residual,
                "status": status,
            }
        )
    return {"claim": "j-function witness coefficients and sums", "status": "PASS" if all_pass else "FAIL", "rows": rows}


def validate_irrep_dimensions() -> dict[str, Any]:
    observed = {name: int(value) for name, value in IRREP_DIMS.items()}
    support = set(SUPERSINGULAR_PRIMES)
    rows = []
    all_pass = observed == PINNED_IRREP_DIMS
    for name, value in observed.items():
        factors = factor_integer(value)
        support_subset = set(factors) <= support
        if value > 1 and not support_subset:
            all_pass = False
        rows.append(
            {
                "name": name,
                "dimension": value,
                "factorization": {str(prime): exponent for prime, exponent in factors.items()},
                "factor_support_subset": support_subset,
            }
        )
    return {
        "claim": "Monster irrep dimensions used by witness subset",
        "status": "PASS" if all_pass else "FAIL",
        "observed": observed,
        "rows": rows,
    }


def validate_gaussian_representatives() -> dict[str, Any]:
    rows = []
    all_pass = True
    for node in base_nodes():
        prime = int(node["prime"])
        x = int(node["x"])
        y = int(node["y"])
        label = node["gaussian_class"]
        expected_label = gaussian_prime_class(prime)
        if prime == 2:
            representative_ok = (x, y) == (1, 1)
        elif prime % 4 == 3:
            representative_ok = (x, y) == (prime, 0)
        else:
            representative_ok = x > 0 and y > 0 and x * x + y * y == prime
        round_trip = representative_for_prime(prime)
        status = "PASS" if label == expected_label and representative_ok and (x, y) == round_trip else "FAIL"
        if status != "PASS":
            all_pass = False
        rows.append(
            {
                "prime": prime,
                "x": x,
                "y": y,
                "gaussian_class": label,
                "expected_class": expected_label,
                "representative_ok": representative_ok,
                "status": status,
            }
        )
    return {"claim": "Gaussian-prime representatives", "status": "PASS" if all_pass else "FAIL", "rows": rows}


def run_source_validation(config: SourceValidationConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    validations = [
        validate_supersingular_primes(),
        validate_monster_order_exponents(),
        validate_j_witnesses(),
        validate_irrep_dimensions(),
        validate_gaussian_representatives(),
    ]
    all_pass = all(validation["status"] == "PASS" for validation in validations)
    labels = [
        "ARITHMETIC_SOURCE_VALIDATION_BUILT",
        "SOURCE_MANIFEST_PINNED_BUILT",
        "SUPERSINGULAR_PRIME_SOURCE_PINNED",
        "MONSTER_ORDER_EXPONENT_SOURCE_PINNED",
        "J_COEFFICIENT_SOURCE_PINNED",
        "IRREP_DIMENSION_SOURCE_PINNED",
        "GAUSSIAN_REPRESENTATIVE_SOURCE_PINNED",
        "MOONSHINE_PROOF_NOT_ESTABLISHED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_PINNED_SOURCE_VALIDATION",
    ]
    labels.append("PINNED_SOURCE_VALIDATION_PASS" if all_pass else "PINNED_SOURCE_VALIDATION_FAIL")

    manifest = {
        "validation_id": "arithmetic_source_validation_v0_1",
        "status": "PASS" if all_pass else "FAIL",
        "classification": "PINNED_SOURCE_VALIDATION_ONLY",
        "config": asdict(config),
        "source_entries": source_entries(),
        "claim_ceiling": "Pinned source validation and local arithmetic checks only; no Moonshine proof or physical bridge claim.",
    }
    manifest["manifest_hash"] = stable_hash("source_manifest_pinned", manifest)

    result = {
        "bridge_id": "arithmetic_source_validation_v0_1",
        "status": "PASS" if all_pass else "FAIL",
        "classification": "PINNED_SOURCE_VALIDATION_ONLY",
        "labels": labels,
        "manifest_hash": manifest["manifest_hash"],
        "validations": validations,
        "claim_ceiling": manifest["claim_ceiling"],
        "non_claims": [
            "not a proof of Monstrous Moonshine",
            "not external source certification beyond pinned references",
            "not a Lean theorem",
            "not a physical bridge",
        ],
    }
    result["result_hash"] = stable_hash("source_validation", result)

    write_json(output_dir / PINNED_MANIFEST_PATH.name, manifest)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result)
    return result


def write_report(path: Path, result: dict[str, Any]) -> None:
    lines = [
        "# Arithmetic Source Validation Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        f"- Pinned manifest hash: `{result['manifest_hash']}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Validation Rows",
        "",
        "| Claim | Status |",
        "| --- | --- |",
    ]
    for validation in result["validations"]:
        lines.append(f"| {validation['claim']} | `{validation['status']}` |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This report checks embedded arithmetic constants against pinned source claims and local arithmetic identities.",
            "It does not prove Monstrous Moonshine, certify the full external literature, add a Lean theorem, or establish a physical bridge.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_source_validation(SourceValidationConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
