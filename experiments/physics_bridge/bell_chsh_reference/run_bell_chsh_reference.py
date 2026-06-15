#!/usr/bin/env python3
"""Bell/CHSH computational reference probe.

This sidecar runs a narrow CHSH sanity check:

- deterministic local response functions obey the CHSH local bound;
- a standard singlet-reference correlation table violates that bound;
- seeded finite-sample telemetry reproduces those reference expectations.

It is not a laboratory Bell experiment, not a CST derivation, and not a claim
that HAOS-IIP recovers Bell correlations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[2]
CONTRACT_PATH = ROOT / "precommitment_contract.json"
CSV_PATH = ROOT / "bell_chsh_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "bell_chsh_reference_report.md"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)

FIELDNAMES = [
    "sample",
    "sample_role",
    "model_family",
    "data_kind",
    "setting_pair",
    "left_setting",
    "right_setting",
    "left_angle_rad",
    "right_angle_rad",
    "expectation",
    "standard_error",
    "trials",
    "left_marginal",
    "right_marginal",
    "chsh_s",
    "chsh_ci_low",
    "chsh_ci_high",
    "local_bound_margin",
    "status",
]

SETTING_PAIRS = (
    ("a", "b", 1.0),
    ("a", "b_prime", 1.0),
    ("a_prime", "b", 1.0),
    ("a_prime", "b_prime", -1.0),
)


@dataclass(frozen=True)
class BellConfig:
    seed: int = 6701
    trials_per_setting: int = 50_000
    ci_z: float = 1.96
    local_bound: float = 2.0
    deterministic_bound_epsilon: float = 1.0e-12
    pass_min_quantum_analytic_s: float = 2.80
    pass_min_quantum_sample_ci_low: float = 2.02
    pass_max_uncorrelated_abs_s: float = 0.12
    sample_local_warning_margin: float = 0.04


@dataclass(frozen=True)
class PairDiagnostic:
    sample: str
    sample_role: str
    model_family: str
    data_kind: str
    setting_pair: str
    left_setting: str
    right_setting: str
    left_angle_rad: float | None
    right_angle_rad: float | None
    expectation: float
    standard_error: float | None
    trials: int
    left_marginal: float | None
    right_marginal: float | None
    chsh_s: float | None
    chsh_ci_low: float | None
    chsh_ci_high: float | None
    local_bound_margin: float | None
    status: str


@dataclass(frozen=True)
class AggregateDiagnostic:
    sample: str
    sample_role: str
    model_family: str
    data_kind: str
    expectations: dict[str, float]
    standard_errors: dict[str, float]
    trials_per_setting: int
    chsh_s: float
    chsh_ci_low: float | None
    chsh_ci_high: float | None
    local_bound_margin: float
    max_abs_left_marginal: float
    max_abs_right_marginal: float
    status: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Bell/CHSH computational reference probe.")
    parser.add_argument("--seed", type=int, default=6701)
    parser.add_argument("--trials-per-setting", type=int, default=50_000)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated contract, CSV, status JSON, and markdown report.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> BellConfig:
    return BellConfig(seed=int(args.seed), trials_per_setting=int(args.trials_per_setting))


def stable_json_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def setting_angles() -> dict[str, float]:
    return {
        "a": 0.0,
        "a_prime": math.pi / 2.0,
        "b": math.pi / 4.0,
        "b_prime": -math.pi / 4.0,
    }


def precommitment_contract(config: BellConfig) -> dict[str, Any]:
    return {
        "name": "Bell CHSH computational reference probe",
        "status": "claim_gated_sidecar",
        "scope": {
            "implemented_fact": "Run a deterministic computational CHSH reference and classical controls.",
            "design_choice": "Use the standard CHSH expression |E(a,b)+E(a,b')+E(a',b)-E(a',b')|.",
            "heuristic": "Finite-sample rows are seeded telemetry around analytic reference expectations.",
            "unverified_hypothesis": "No CST or HAOS Bell-recovery hypothesis is tested here.",
        },
        "non_claims": [
            "not a laboratory Bell experiment",
            "not a loophole-free Bell test",
            "not a CST derivation of Bell correlations",
            "not evidence that HAOS-IIP replaces quantum mechanics",
            "not evidence for a Born-rule derivation",
        ],
        "settings": {
            "angles_radians": setting_angles(),
            "setting_pairs": [
                {"left": left, "right": right, "coefficient": coeff}
                for left, right, coeff in SETTING_PAIRS
            ],
        },
        "models": {
            "deterministic_local_bound": {
                "purpose": "Enumerate all deterministic local two-setting response functions.",
                "expected_response": "maximum CHSH S must be <= 2.",
            },
            "local_hidden_variable_reference": {
                "purpose": "Seeded classical control with shared hidden angle and local sign responses.",
                "expected_response": "sampled CHSH should not materially exceed the local bound.",
            },
            "quantum_singlet_reference": {
                "purpose": "Standard analytic reference E(theta_left, theta_right) = -cos(delta theta).",
                "expected_response": "analytic S should be near 2*sqrt(2).",
            },
            "uncorrelated_control": {
                "purpose": "Independent random outcomes with zero expected correlations.",
                "expected_response": "sampled CHSH should remain near zero.",
            },
        },
        "thresholds": {
            "local_bound": config.local_bound,
            "deterministic_bound_epsilon": config.deterministic_bound_epsilon,
            "pass_min_quantum_analytic_s": config.pass_min_quantum_analytic_s,
            "pass_min_quantum_sample_ci_low": config.pass_min_quantum_sample_ci_low,
            "pass_max_uncorrelated_abs_s": config.pass_max_uncorrelated_abs_s,
            "sample_local_warning_margin": config.sample_local_warning_margin,
        },
        "verdict_logic": {
            "reference_pass": [
                "deterministic local maximum <= local_bound + epsilon",
                "quantum analytic S >= pass_min_quantum_analytic_s",
                "quantum sampled lower CI >= pass_min_quantum_sample_ci_low",
                "local hidden-angle sampled S <= local_bound + sample_local_warning_margin",
                "abs(uncorrelated sampled S) <= pass_max_uncorrelated_abs_s",
            ],
            "always_report": [
                "CST_BELL_STATUS_OPEN",
                "HAOS_DERIVATION_NOT_TESTED",
            ],
        },
        "seed": config.seed,
        "trials_per_setting": config.trials_per_setting,
        "ci_z": config.ci_z,
    }


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def chsh_s(expectations: dict[str, float]) -> float:
    signed = 0.0
    for left, right, coefficient in SETTING_PAIRS:
        signed += coefficient * expectations[f"{left}_{right}"]
    return abs(float(signed))


def chsh_ci(expectations: dict[str, float], standard_errors: dict[str, float], z_value: float) -> tuple[float, float]:
    signed = 0.0
    variance = 0.0
    for left, right, coefficient in SETTING_PAIRS:
        key = f"{left}_{right}"
        signed += coefficient * expectations[key]
        variance += standard_errors[key] * standard_errors[key]
    s_value = abs(float(signed))
    radius = z_value * math.sqrt(max(variance, 0.0))
    return max(0.0, s_value - radius), s_value + radius


def quantum_expectations() -> dict[str, float]:
    angles = setting_angles()
    values: dict[str, float] = {}
    for left, right, _coefficient in SETTING_PAIRS:
        values[f"{left}_{right}"] = -math.cos(angles[left] - angles[right])
    return values


def enumerate_deterministic_local_bound() -> tuple[float, list[dict[str, Any]]]:
    response_values = (-1.0, 1.0)
    records: list[dict[str, Any]] = []
    max_s = 0.0
    for a in response_values:
        for a_prime in response_values:
            for b in response_values:
                for b_prime in response_values:
                    expectations = {
                        "a_b": a * b,
                        "a_b_prime": a * b_prime,
                        "a_prime_b": a_prime * b,
                        "a_prime_b_prime": a_prime * b_prime,
                    }
                    s_value = chsh_s(expectations)
                    max_s = max(max_s, s_value)
                    records.append(
                        {
                            "A(a)": a,
                            "A(a_prime)": a_prime,
                            "B(b)": b,
                            "B(b_prime)": b_prime,
                            "chsh_s": s_value,
                        }
                    )
    return float(max_s), records


def sign_nonzero(values: np.ndarray) -> np.ndarray:
    signs = np.sign(values)
    signs[signs == 0.0] = 1.0
    return signs.astype(float)


def sample_quantum_pair(expectation: float, trials: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    left = rng.choice(np.array([-1.0, 1.0]), size=trials)
    same_probability = float(np.clip((1.0 + expectation) / 2.0, 0.0, 1.0))
    same = rng.random(trials) < same_probability
    right = np.where(same, left, -left)
    return left.astype(float), right.astype(float)


def sample_local_hidden_pair(
    left_angle: float,
    right_angle: float,
    trials: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    hidden = rng.uniform(0.0, 2.0 * math.pi, size=trials)
    left = sign_nonzero(np.cos(hidden - left_angle))
    right = -sign_nonzero(np.cos(hidden - right_angle))
    return left, right


def sample_uncorrelated_pair(trials: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    values = np.array([-1.0, 1.0], dtype=float)
    return rng.choice(values, size=trials), rng.choice(values, size=trials)


def product_stats(left: np.ndarray, right: np.ndarray) -> tuple[float, float, float, float]:
    products = left * right
    expectation = float(np.mean(products))
    if len(products) <= 1:
        standard_error = 0.0
    else:
        standard_error = float(np.std(products, ddof=1) / math.sqrt(len(products)))
    return expectation, standard_error, float(np.mean(left)), float(np.mean(right))


def stable_sample_offset(sample: str) -> int:
    digest = hashlib.sha256(sample.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % 1_000_000


def aggregate_sample(
    sample: str,
    sample_role: str,
    model_family: str,
    data_kind: str,
    config: BellConfig,
    sampler: Any,
    analytic_expectations: dict[str, float] | None = None,
) -> tuple[AggregateDiagnostic, list[PairDiagnostic]]:
    angles = setting_angles()
    rng = np.random.default_rng(config.seed + stable_sample_offset(sample))
    expectations: dict[str, float] = {}
    standard_errors: dict[str, float] = {}
    left_marginals: dict[str, float] = {}
    right_marginals: dict[str, float] = {}
    pair_rows: list[PairDiagnostic] = []

    for left_setting, right_setting, _coefficient in SETTING_PAIRS:
        key = f"{left_setting}_{right_setting}"
        if analytic_expectations is not None:
            left, right = sampler(analytic_expectations[key], config.trials_per_setting, rng)
        else:
            left, right = sampler(angles[left_setting], angles[right_setting], config.trials_per_setting, rng)
        expectation, standard_error, left_marginal, right_marginal = product_stats(left, right)
        expectations[key] = expectation
        standard_errors[key] = standard_error
        left_marginals[key] = left_marginal
        right_marginals[key] = right_marginal
        pair_rows.append(
            PairDiagnostic(
                sample=sample,
                sample_role=sample_role,
                model_family=model_family,
                data_kind=data_kind,
                setting_pair=key,
                left_setting=left_setting,
                right_setting=right_setting,
                left_angle_rad=angles[left_setting],
                right_angle_rad=angles[right_setting],
                expectation=expectation,
                standard_error=standard_error,
                trials=config.trials_per_setting,
                left_marginal=left_marginal,
                right_marginal=right_marginal,
                chsh_s=None,
                chsh_ci_low=None,
                chsh_ci_high=None,
                local_bound_margin=None,
                status="PAIR_RECORDED",
            )
        )

    s_value = chsh_s(expectations)
    ci_low, ci_high = chsh_ci(expectations, standard_errors, config.ci_z)
    margin = s_value - config.local_bound
    if model_family == "quantum_singlet_reference":
        status = "REFERENCE_VIOLATION_REPRODUCED" if ci_low >= config.pass_min_quantum_sample_ci_low else "REFERENCE_WEAK_OR_UNSTABLE"
    elif model_family == "uncorrelated_control":
        status = "CONTROL_NEAR_ZERO" if abs(s_value) <= config.pass_max_uncorrelated_abs_s else "CONTROL_DRIFT"
    else:
        status = "CLASSICAL_CONTROL_WITHIN_WARNING_MARGIN" if s_value <= config.local_bound + config.sample_local_warning_margin else "CLASSICAL_CONTROL_SAMPLE_DRIFT"

    aggregate = AggregateDiagnostic(
        sample=sample,
        sample_role=sample_role,
        model_family=model_family,
        data_kind=data_kind,
        expectations=expectations,
        standard_errors=standard_errors,
        trials_per_setting=config.trials_per_setting,
        chsh_s=s_value,
        chsh_ci_low=ci_low,
        chsh_ci_high=ci_high,
        local_bound_margin=margin,
        max_abs_left_marginal=max(abs(value) for value in left_marginals.values()),
        max_abs_right_marginal=max(abs(value) for value in right_marginals.values()),
        status=status,
    )
    return aggregate, pair_rows


def aggregate_analytic(
    sample: str,
    sample_role: str,
    model_family: str,
    expectations: dict[str, float],
    config: BellConfig,
    status: str,
) -> tuple[AggregateDiagnostic, list[PairDiagnostic]]:
    angles = setting_angles()
    zero_errors = {key: 0.0 for key in expectations}
    s_value = chsh_s(expectations)
    pair_rows: list[PairDiagnostic] = []
    for left_setting, right_setting, _coefficient in SETTING_PAIRS:
        key = f"{left_setting}_{right_setting}"
        pair_rows.append(
            PairDiagnostic(
                sample=sample,
                sample_role=sample_role,
                model_family=model_family,
                data_kind="analytic",
                setting_pair=key,
                left_setting=left_setting,
                right_setting=right_setting,
                left_angle_rad=angles[left_setting],
                right_angle_rad=angles[right_setting],
                expectation=expectations[key],
                standard_error=None,
                trials=0,
                left_marginal=None,
                right_marginal=None,
                chsh_s=None,
                chsh_ci_low=None,
                chsh_ci_high=None,
                local_bound_margin=None,
                status="PAIR_RECORDED",
            )
        )
    aggregate = AggregateDiagnostic(
        sample=sample,
        sample_role=sample_role,
        model_family=model_family,
        data_kind="analytic",
        expectations=expectations,
        standard_errors=zero_errors,
        trials_per_setting=0,
        chsh_s=s_value,
        chsh_ci_low=None,
        chsh_ci_high=None,
        local_bound_margin=s_value - config.local_bound,
        max_abs_left_marginal=0.0,
        max_abs_right_marginal=0.0,
        status=status,
    )
    return aggregate, pair_rows


def aggregate_rows(aggregates: list[AggregateDiagnostic]) -> list[PairDiagnostic]:
    rows: list[PairDiagnostic] = []
    for aggregate in aggregates:
        rows.append(
            PairDiagnostic(
                sample=aggregate.sample,
                sample_role=aggregate.sample_role,
                model_family=aggregate.model_family,
                data_kind=aggregate.data_kind,
                setting_pair="aggregate",
                left_setting="",
                right_setting="",
                left_angle_rad=None,
                right_angle_rad=None,
                expectation=float("nan"),
                standard_error=None,
                trials=aggregate.trials_per_setting,
                left_marginal=aggregate.max_abs_left_marginal,
                right_marginal=aggregate.max_abs_right_marginal,
                chsh_s=aggregate.chsh_s,
                chsh_ci_low=aggregate.chsh_ci_low,
                chsh_ci_high=aggregate.chsh_ci_high,
                local_bound_margin=aggregate.local_bound_margin,
                status=aggregate.status,
            )
        )
    return rows


def row_dict(row: PairDiagnostic) -> dict[str, Any]:
    payload = asdict(row)
    for key, value in list(payload.items()):
        if isinstance(value, float):
            if math.isnan(value):
                payload[key] = ""
            else:
                payload[key] = f"{value:.12g}"
        elif value is None:
            payload[key] = ""
    return payload


def write_csv(path: Path, rows: list[PairDiagnostic]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row_dict(row))


def run_probe(config: BellConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = precommitment_contract(config)
    contract["contract_hash"] = stable_json_hash("bell_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)

    deterministic_max_s, deterministic_records = enumerate_deterministic_local_bound()
    deterministic_status = (
        "CLASSICAL_BOUND_REPRODUCED"
        if deterministic_max_s <= config.local_bound + config.deterministic_bound_epsilon
        else "CLASSICAL_BOUND_FAILURE"
    )
    deterministic_aggregate = AggregateDiagnostic(
        sample="deterministic_local_response_enumeration",
        sample_role="classical_bound",
        model_family="deterministic_local_bound",
        data_kind="exhaustive_enumeration",
        expectations={},
        standard_errors={},
        trials_per_setting=len(deterministic_records),
        chsh_s=deterministic_max_s,
        chsh_ci_low=None,
        chsh_ci_high=None,
        local_bound_margin=deterministic_max_s - config.local_bound,
        max_abs_left_marginal=0.0,
        max_abs_right_marginal=0.0,
        status=deterministic_status,
    )

    quantum_analytic_expectations = quantum_expectations()
    quantum_analytic_s = chsh_s(quantum_analytic_expectations)
    quantum_analytic_status = (
        "REFERENCE_ANALYTIC_VIOLATION_REPRODUCED"
        if quantum_analytic_s >= config.pass_min_quantum_analytic_s
        else "REFERENCE_ANALYTIC_WEAK"
    )
    quantum_analytic_aggregate, quantum_analytic_rows = aggregate_analytic(
        "quantum_singlet_reference_analytic",
        "reference",
        "quantum_singlet_reference",
        quantum_analytic_expectations,
        config,
        quantum_analytic_status,
    )

    quantum_sample_aggregate, quantum_sample_rows = aggregate_sample(
        "quantum_singlet_reference_sampled",
        "reference",
        "quantum_singlet_reference",
        "seeded_finite_sample",
        config,
        sample_quantum_pair,
        quantum_analytic_expectations,
    )
    local_sample_aggregate, local_sample_rows = aggregate_sample(
        "local_hidden_angle_control_sampled",
        "negative_control",
        "local_hidden_variable_reference",
        "seeded_finite_sample",
        config,
        sample_local_hidden_pair,
    )

    def uncorr_sampler(
        _left_angle: float,
        _right_angle: float,
        trials: int,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, np.ndarray]:
        return sample_uncorrelated_pair(trials, rng)

    uncorrelated_aggregate, uncorrelated_rows = aggregate_sample(
        "uncorrelated_outcome_control_sampled",
        "negative_control",
        "uncorrelated_control",
        "seeded_finite_sample",
        config,
        uncorr_sampler,
    )

    aggregates = [
        deterministic_aggregate,
        quantum_analytic_aggregate,
        quantum_sample_aggregate,
        local_sample_aggregate,
        uncorrelated_aggregate,
    ]
    rows = (
        aggregate_rows(aggregates)
        + quantum_analytic_rows
        + quantum_sample_rows
        + local_sample_rows
        + uncorrelated_rows
    )
    write_csv(output_dir / CSV_PATH.name, rows)

    reference_pass = (
        deterministic_status == "CLASSICAL_BOUND_REPRODUCED"
        and quantum_analytic_status == "REFERENCE_ANALYTIC_VIOLATION_REPRODUCED"
        and quantum_sample_aggregate.status == "REFERENCE_VIOLATION_REPRODUCED"
        and local_sample_aggregate.status == "CLASSICAL_CONTROL_WITHIN_WARNING_MARGIN"
        and uncorrelated_aggregate.status == "CONTROL_NEAR_ZERO"
    )
    labels = [
        "BELL_REFERENCE_SANITY_PASS" if reference_pass else "BELL_REFERENCE_SANITY_FAIL",
        "CST_BELL_STATUS_OPEN",
        "HAOS_DERIVATION_NOT_TESTED",
        "NO_PHYSICAL_EXPERIMENT",
    ]

    result: dict[str, Any] = {
        "labels": labels,
        "contract_hash": contract["contract_hash"],
        "config": asdict(config),
        "settings": setting_angles(),
        "chsh_convention": "abs(E(a,b)+E(a,b_prime)+E(a_prime,b)-E(a_prime,b_prime))",
        "deterministic_local_max_s": deterministic_max_s,
        "quantum_analytic_s": quantum_analytic_s,
        "aggregates": [asdict(item) for item in aggregates],
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "bell_chsh_diagnostics": repo_rel(output_dir / CSV_PATH.name),
            "bridge_status": repo_rel(output_dir / STATUS_PATH.name),
            "bell_chsh_reference_report": repo_rel(output_dir / REPORT_PATH.name),
        },
        "limitations": contract["non_claims"],
    }
    result["result_hash"] = stable_json_hash("bell_chsh_result", result)
    write_json(output_dir / STATUS_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, aggregates)
    return result


def write_report(path: Path, result: dict[str, Any], aggregates: list[AggregateDiagnostic]) -> None:
    by_sample = {item.sample: item for item in aggregates}
    quantum_sample = by_sample["quantum_singlet_reference_sampled"]
    local_sample = by_sample["local_hidden_angle_control_sampled"]
    uncorrelated = by_sample["uncorrelated_outcome_control_sampled"]
    lines = [
        "# Bell/CHSH Computational Reference Probe",
        "",
        "Implemented fact: this sidecar runs a deterministic CHSH reference calculation and seeded finite-sample controls.",
        "Design choice: CHSH uses `|E(a,b)+E(a,b')+E(a',b)-E(a',b')|` with the precommitted angle table.",
        "Heuristic: sampled rows are finite-sample telemetry around analytic reference or classical-control generators.",
        "Unverified hypothesis: no CST or HAOS Bell-recovery hypothesis is tested here.",
        "",
        "## Verdict Labels",
    ]
    lines.extend(f"- {label}" for label in result["labels"])
    lines.extend(
        [
            "",
            "## CHSH Results",
            f"- deterministic local response max S: `{result['deterministic_local_max_s']:.12f}`",
            f"- quantum analytic reference S: `{result['quantum_analytic_s']:.12f}`",
            (
                "- quantum sampled reference S: "
                f"`{quantum_sample.chsh_s:.12f}` "
                f"(95% CI `{quantum_sample.chsh_ci_low:.12f}` to `{quantum_sample.chsh_ci_high:.12f}`)"
            ),
            (
                "- local hidden-angle control sampled S: "
                f"`{local_sample.chsh_s:.12f}` "
                f"(95% CI `{local_sample.chsh_ci_low:.12f}` to `{local_sample.chsh_ci_high:.12f}`)"
            ),
            (
                "- uncorrelated control sampled S: "
                f"`{uncorrelated.chsh_s:.12f}` "
                f"(95% CI `{uncorrelated.chsh_ci_low:.12f}` to `{uncorrelated.chsh_ci_high:.12f}`)"
            ),
            "",
            "## Boundary",
            "- This is not a laboratory Bell experiment.",
            "- This does not test loophole closure, detector effects, locality constraints, or real apparatus data.",
            "- This does not derive Bell correlations from CST or HAOS-IIP.",
            "- This does not change the frozen CST v0.2.2 checkpoint.",
            "- The sampled local hidden-angle row is a finite-sample control; the exhaustive local enumeration is the actual local-bound check.",
            "- `BELL_REFERENCE_SANITY_PASS` means only that the computational reference harness behaves as expected.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_probe(make_config(args), args.output_dir)
    print(json.dumps({"labels": result["labels"], "result_hash": result["result_hash"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
