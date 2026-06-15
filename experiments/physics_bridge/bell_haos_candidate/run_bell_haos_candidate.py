#!/usr/bin/env python3
"""HAOS-IIP-to-Bell candidate bridge.

This bundle tests candidate-generated Bell/CHSH correlation tables against a
frozen scoreboard convention. It does not modify or import the frozen Bell
reference sidecar, and it does not attempt to derive quantum mechanics.

B0 defines the neutral adapter.
B1 runs a deliberately local sanity candidate.
B2 runs a HAOS-IIP-local recoverability candidate using frozen Phase 18 data.
B3 runs an explicit joint-closure-cost candidate after a frozen assumption
invoice. It is still not a derivation of quantum correlations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent
PHASE18_SOURCE = REPO_ROOT / "phase18-distance-surrogate" / "runs" / "phase18_shell_ordering_metrics.csv"
REFERENCE_SIDECAR = REPO_ROOT / "experiments" / "physics_bridge" / "bell_chsh_reference"

CONTRACT_PATH = ROOT / "precommitment_contract.json"
ASSUMPTION_LEDGER_PATH = ROOT / "assumption_ledger.json"
TRIALS_PATH = ROOT / "candidate_trials.csv"
CORRELATIONS_PATH = ROOT / "candidate_correlations.json"
NO_SIGNALLING_PATH = ROOT / "no_signalling_diagnostics.csv"
SETTING_INDEPENDENCE_PATH = ROOT / "setting_independence_diagnostics.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "bell_candidate_report.md"
RESULT_PATH = ROOT / "bell_candidate_result.json"
B3_CONTRACT_PATH = ROOT / "b3_precommitment_contract.json"
B3_REPORT_PATH = ROOT / "b3_precommitment_report.md"
B3_COST_AUDIT_PATH = ROOT / "b3_joint_cost_audit.csv"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)

SETTING_PAIRS = (("a0", "b0"), ("a0", "b1"), ("a1", "b0"), ("a1", "b1"))
PAIR_KEYS = ("E00", "E01", "E10", "E11")
CHSH_SIGNS = {"E00": 1.0, "E01": 1.0, "E10": 1.0, "E11": -1.0}
LOCAL_BOUND = 2.0
FORBIDDEN_TARGET_PATTERNS = (
    "quantum_expectations",
    "bell_chsh_reference",
    "bridge_status.json",
    "target_chsh",
    "target_s",
    "reference_outcomes",
    "post_selected_target",
    "2*sqrt(2)",
    "2 * sqrt(2)",
    "sqrt(2)",
    "-cos(",
    "cos(theta",
    "cos(2 * delta",
    "cos(2*delta",
)


TRIAL_FIELDNAMES = [
    "candidate_id",
    "run_id",
    "stage",
    "seed",
    "trial_index",
    "setting_pair",
    "alice_setting",
    "bob_setting",
    "source_id",
    "source_bucket",
    "source_artifact_hash",
    "alice_outcome",
    "bob_outcome",
    "product",
    "retained",
    "rejection_reason",
    "alice_local_score",
    "bob_local_score",
    "provenance_note",
]

NO_SIGNALLING_FIELDNAMES = [
    "candidate_id",
    "run_id",
    "party",
    "local_setting",
    "remote_setting_left",
    "remote_setting_right",
    "p_plus_left",
    "p_plus_right",
    "deviation",
    "tolerance",
    "status",
]

SETTING_INDEPENDENCE_FIELDNAMES = [
    "candidate_id",
    "run_id",
    "mutual_information_bits",
    "tolerance_bits",
    "source_bucket_count",
    "setting_pair_count",
    "status",
]

CONTROL_FIELDNAMES = [
    "candidate_id",
    "run_id",
    "stage",
    "sample_role",
    "S",
    "ci_low",
    "ci_high",
    "E00",
    "E01",
    "E10",
    "E11",
    "retained_trials",
    "rejected_trials",
    "retention_rate",
    "no_signalling_max_deviation",
    "setting_source_mutual_information_bits",
    "target_leakage_detected",
    "post_selection_detected",
    "verdict_labels",
]

B3_COST_AUDIT_FIELDNAMES = [
    "candidate_id",
    "run_id",
    "seed",
    "trial_index",
    "setting_pair",
    "source_id",
    "source_bucket",
    "base_instability",
    "closure_strength",
    "orientation_score",
    "preferred_product",
    "cost_product_plus",
    "cost_product_minus",
    "p_product_plus",
    "sampled_product",
    "source_features",
]


@dataclass(frozen=True)
class BellSetting:
    label: str
    side: str
    index: int
    perturbation_label: str


@dataclass(frozen=True)
class BellCandidateConfig:
    version: str = "bell-haos-candidate-v0.2-b3"
    seeds: tuple[int, ...] = (7201, 7202, 7203)
    trials_per_setting_pair: int = 384
    setting_pairs: tuple[tuple[str, str], ...] = SETTING_PAIRS
    local_bound: float = LOCAL_BOUND
    ci_z: float = 1.96
    finite_sample_tolerance: float = 0.04
    no_signalling_tolerance: float = 0.08
    setting_independence_mi_tolerance_bits: float = 0.002
    source_bucket_count: int = 16
    b2_source_artifact: str = "phase18-distance-surrogate/runs/phase18_shell_ordering_metrics.csv"
    b2_hierarchy_label: str = "frozen_branch"
    b2_source_limit: int = 6
    b2_threshold_a0: float = 0.70
    b2_threshold_a1: float = 0.76
    b2_threshold_b0: float = 0.68
    b2_threshold_b1: float = 0.74
    b3_beta: float = 1.35
    b3_strength_floor: float = 0.05
    b3_strength_ceiling: float = 0.95
    b3_joint_model: str = "phase18_chain_order_joint_closure_cost"


@dataclass(frozen=True)
class BellCandidateProvenance:
    code_path: str
    code_hash: str
    reference_sidecar_path: str
    reference_sidecar_hashes: dict[str, str]
    source_artifact: str
    source_artifact_hash: str
    config_hash: str


@dataclass(frozen=True)
class BellOutcome:
    alice: int
    bob: int
    alice_local_score: float
    bob_local_score: float
    retained: bool = True
    rejection_reason: str = ""


@dataclass(frozen=True)
class BellTrial:
    candidate_id: str
    run_id: str
    stage: str
    seed: int
    trial_index: int
    setting_pair: str
    alice_setting: str
    bob_setting: str
    source_id: str
    source_bucket: int
    source_artifact_hash: str
    outcome: BellOutcome
    provenance_note: str


@dataclass(frozen=True)
class BellCorrelationTable:
    candidate_id: str
    run_id: str
    stage: str
    sample_role: str
    correlations: dict[str, float]
    standard_errors: dict[str, float]
    trials_by_pair: dict[str, int]
    alice_marginals_by_pair: dict[str, float]
    bob_marginals_by_pair: dict[str, float]
    S: float
    ci_low: float
    ci_high: float
    retained_trials: int
    rejected_trials: int
    retention_rate: float


@dataclass(frozen=True)
class BellAssumptionLedger:
    candidate_id: str
    run_id: str
    stage: str
    factorization: str
    setting_independence: str
    parameter_independence: str
    outcome_independence: str
    forward_only_causality: str
    remote_setting_access: bool
    joint_setting_access_declared: bool
    post_selection_detected: bool
    target_leakage_detected: bool
    assumption_violations: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class BellCandidateResult:
    candidate_id: str
    run_id: str
    stage: str
    sample_role: str
    correlation_table: BellCorrelationTable
    assumption_ledger: BellAssumptionLedger
    no_signalling_max_deviation: float
    setting_source_mutual_information_bits: float
    verdict_labels: list[str]


@dataclass(frozen=True)
class SourceState:
    source_id: str
    source_bucket: int
    source_hash: str
    fields: dict[str, Any]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run isolated HAOS-IIP-to-Bell candidate bridge.")
    parser.add_argument("--trials-per-setting-pair", type=int, default=BellCandidateConfig.trials_per_setting_pair)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(BellCandidateConfig.seeds))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated contract, ledgers, CSVs, JSON, and report.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> BellCandidateConfig:
    return BellCandidateConfig(
        seeds=tuple(int(seed) for seed in args.seeds),
        trials_per_setting_pair=int(args.trials_per_setting_pair),
    )


def stable_json_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_unit_interval(*parts: Any) -> float:
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()
    return int(digest[:16], 16) / float(16**16 - 1)


def stable_sign(*parts: Any) -> int:
    return 1 if stable_unit_interval(*parts) >= 0.5 else -1


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def source_bucket(source_id: str, bucket_count: int) -> int:
    digest = hashlib.sha256(source_id.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % bucket_count


def pair_key(alice_setting: str, bob_setting: str) -> str:
    return {
        ("a0", "b0"): "E00",
        ("a0", "b1"): "E01",
        ("a1", "b0"): "E10",
        ("a1", "b1"): "E11",
    }[(alice_setting, bob_setting)]


def pair_label(alice_setting: str, bob_setting: str) -> str:
    return f"{alice_setting}_{bob_setting}"


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


def code_hash() -> str:
    return hash_file(Path(__file__).resolve())


def reference_sidecar_hashes() -> dict[str, str]:
    paths = [
        REFERENCE_SIDECAR / "run_bell_chsh_reference.py",
        REFERENCE_SIDECAR / "precommitment_contract.json",
        REFERENCE_SIDECAR / "bridge_status.json",
    ]
    return {str(path.relative_to(REPO_ROOT)): hash_file(path) for path in paths if path.exists()}


def read_phase18_sources(config: BellCandidateConfig) -> list[SourceState]:
    source_path = REPO_ROOT / config.b2_source_artifact
    artifact_hash = hash_file(source_path)
    states: list[SourceState] = []
    with source_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("hierarchy_label") != config.b2_hierarchy_label:
                continue
            fields = {
                "hierarchy_label": row["hierarchy_label"],
                "n_side": int(row["n_side"]),
                "seed": int(row["seed"]),
                "probe_name": row["probe_name"],
                "shell_order_score": float(row["shell_order_score"]),
                "max_shell_overlap_fraction": float(row["max_shell_overlap_fraction"]),
                "near_mean_arrival": float(row["near_mean_arrival"]),
                "far_mean_arrival": float(row["far_mean_arrival"]),
                "phase17_mismatch_rate": float(row["phase17_mismatch_rate"]),
                "chain_signature": row["chain_signature"],
            }
            source_id = stable_json_hash("phase18_source", fields)
            states.append(
                SourceState(
                    source_id=source_id,
                    source_bucket=source_bucket(source_id, config.source_bucket_count),
                    source_hash=artifact_hash,
                    fields=fields,
                )
            )
            if len(states) >= config.b2_source_limit:
                break
    if not states:
        raise ValueError(f"no B2 source states found in {source_path}")
    return states


def synthetic_source_state(seed: int, trial_index: int, config: BellCandidateConfig) -> SourceState:
    fields = {"seed": seed, "trial_index": trial_index, "kind": "balanced_synthetic_local_lambda"}
    source_id = stable_json_hash("synthetic_lambda", fields)
    return SourceState(
        source_id=source_id,
        source_bucket=source_bucket(source_id, config.source_bucket_count),
        source_hash=stable_json_hash("synthetic_source_family", {"seed": seed}),
        fields=fields,
    )


def local_sanity_outcome(source: SourceState, party: str, local_setting: str) -> tuple[int, float]:
    setting_index = 0 if local_setting.endswith("0") else 1
    score = stable_unit_interval("B1", source.source_id, party, local_setting)
    hidden_bit = stable_sign("B1_hidden", source.source_id, setting_index)
    party_flip = 1 if party == "alice" else stable_sign("B1_bob_flip", source.source_id, setting_index)
    return hidden_bit * party_flip, score


def haos_local_recoverability_outcome(source: SourceState, party: str, local_setting: str, config: BellCandidateConfig) -> tuple[int, float]:
    fields = source.fields
    shell_order = float(fields["shell_order_score"])
    overlap = float(fields["max_shell_overlap_fraction"])
    mismatch = float(fields["phase17_mismatch_rate"])
    near = float(fields["near_mean_arrival"])
    far = float(fields["far_mean_arrival"])
    arrival_gap = max(0.0, min(1.0, (far - near) / max(far, 1.0e-12)))
    chain_bias = 0.04 * stable_sign("chain", fields["chain_signature"], party, local_setting)
    setting_penalty = 0.03 if local_setting.endswith("1") else 0.0
    party_bias = 0.015 if party == "alice" else -0.015
    recoverability_proxy = shell_order - 0.55 * overlap - 0.75 * mismatch + 0.10 * arrival_gap
    score = recoverability_proxy + chain_bias + party_bias - setting_penalty
    thresholds = {
        ("alice", "a0"): config.b2_threshold_a0,
        ("alice", "a1"): config.b2_threshold_a1,
        ("bob", "b0"): config.b2_threshold_b0,
        ("bob", "b1"): config.b2_threshold_b1,
    }
    threshold = thresholds[(party, local_setting)]
    if abs(score - threshold) <= 1.0e-12:
        return stable_sign("B2_tie", source.source_id, party, local_setting), score
    return (1 if score > threshold else -1), score


def chain_rank_fraction(chain_signature: str, token: str) -> float:
    tokens = chain_signature.split(">")
    if len(tokens) <= 1 or token not in tokens:
        return 0.5
    return tokens.index(token) / float(len(tokens) - 1)


def b3_setting_tokens(setting: str) -> tuple[str, str]:
    tokens = {
        "a0": ("radius_half", "front_near"),
        "a1": ("spectral_half", "width_half"),
        "b0": ("dispersion_half", "front_near"),
        "b1": ("low_k_half", "front_far"),
    }
    return tokens[setting]


def b3_source_features(source: SourceState, config: BellCandidateConfig) -> dict[str, float]:
    fields = source.fields
    shell_order = float(fields["shell_order_score"])
    overlap = float(fields["max_shell_overlap_fraction"])
    mismatch = float(fields["phase17_mismatch_rate"])
    near = float(fields["near_mean_arrival"])
    far = float(fields["far_mean_arrival"])
    arrival_gap = clamp((far - near) / max(far, 1.0e-12), 0.0, 1.0)
    closure_strength = clamp(
        shell_order * (1.0 - overlap) * (1.0 - mismatch) + 0.25 * arrival_gap,
        config.b3_strength_floor,
        config.b3_strength_ceiling,
    )
    base_instability = clamp((1.0 - shell_order) + overlap + mismatch + (1.0 - arrival_gap), 0.0, 3.0)
    return {
        "shell_order_score": shell_order,
        "max_shell_overlap_fraction": overlap,
        "phase17_mismatch_rate": mismatch,
        "arrival_gap": arrival_gap,
        "closure_strength": closure_strength,
        "base_instability": base_instability,
    }


def b3_orientation_score(source: SourceState, alice_setting: str, bob_setting: str) -> float:
    chain_signature = str(source.fields["chain_signature"])
    alice_primary, alice_secondary = b3_setting_tokens(alice_setting)
    bob_primary, bob_secondary = b3_setting_tokens(bob_setting)
    primary_delta = chain_rank_fraction(chain_signature, alice_primary) - chain_rank_fraction(chain_signature, bob_primary)
    secondary_delta = chain_rank_fraction(chain_signature, alice_secondary) - chain_rank_fraction(chain_signature, bob_secondary)
    return primary_delta + 0.5 * secondary_delta


def b3_joint_closure_costs(
    source: SourceState,
    alice_setting: str,
    bob_setting: str,
    config: BellCandidateConfig,
) -> dict[str, float]:
    features = b3_source_features(source, config)
    orientation_score = b3_orientation_score(source, alice_setting, bob_setting)
    preferred_product = 1 if orientation_score >= 0.0 else -1
    signed_strength = preferred_product * features["closure_strength"]
    base = features["base_instability"]
    cost_plus = base + 0.5 * (1.0 - signed_strength)
    cost_minus = base + 0.5 * (1.0 + signed_strength)
    weight_plus = math.exp(-config.b3_beta * cost_plus)
    weight_minus = math.exp(-config.b3_beta * cost_minus)
    p_plus = weight_plus / max(weight_plus + weight_minus, 1.0e-300)
    return {
        **features,
        "orientation_score": orientation_score,
        "preferred_product": float(preferred_product),
        "cost_product_plus": cost_plus,
        "cost_product_minus": cost_minus,
        "p_product_plus": p_plus,
    }


def b3_sample_product(source: SourceState, alice_setting: str, bob_setting: str, seed: int, trial_index: int, config: BellCandidateConfig) -> tuple[int, dict[str, float]]:
    costs = b3_joint_closure_costs(source, alice_setting, bob_setting, config)
    draw = stable_unit_interval("B3_joint_product", seed, trial_index // 2, source.source_id, alice_setting, bob_setting)
    product = 1 if draw <= costs["p_product_plus"] else -1
    return product, costs


def independent_random_outcome(source: SourceState, party: str, local_setting: str) -> tuple[int, float]:
    score = stable_unit_interval("independent_random", source.source_id, party, local_setting)
    return (1 if score >= 0.5 else -1), score


def deterministic_lhv_outcome(source: SourceState, party: str, local_setting: str) -> tuple[int, float]:
    setting_index = 0 if local_setting.endswith("0") else 1
    score = stable_unit_interval("deterministic_lhv", source.source_id, party, setting_index)
    return stable_sign("deterministic_lhv_sign", source.source_id, party, setting_index), score


def shared_randomness_outcome(source: SourceState, party: str, local_setting: str) -> tuple[int, float]:
    setting_index = 0 if local_setting.endswith("0") else 1
    shared = stable_sign("shared_randomness", source.source_id)
    local_flip = stable_sign("shared_local_flip", party, setting_index)
    score = stable_unit_interval("shared_randomness_score", source.source_id, party, setting_index)
    return shared * local_flip, score


def detect_target_leakage(source_texts: Iterable[str]) -> tuple[bool, list[str]]:
    hits: list[str] = []
    normalized_texts = [text.replace(" ", "").lower() for text in source_texts]
    for pattern in FORBIDDEN_TARGET_PATTERNS:
        normalized_pattern = pattern.replace(" ", "").lower()
        if any(normalized_pattern in text for text in normalized_texts):
            hits.append(pattern)
    return bool(hits), sorted(set(hits))


def generator_sources(functions: Iterable[Callable[..., Any]]) -> list[str]:
    return [inspect.getsource(function) for function in functions]


def balanced_trials(
    candidate_id: str,
    stage: str,
    seed: int,
    sources: list[SourceState],
    config: BellCandidateConfig,
    outcome_function: Callable[..., tuple[int, float]],
    provenance_note: str,
    post_select: bool = False,
    leak_source_to_setting: bool = False,
) -> list[BellTrial]:
    trials: list[BellTrial] = []
    run_id = f"{candidate_id}_seed_{seed}"
    for pair_index, (alice_setting, bob_setting) in enumerate(config.setting_pairs):
        setting_pair = pair_label(alice_setting, bob_setting)
        chsh_sign = 1 if pair_key(alice_setting, bob_setting) != "E11" else -1
        for trial_index in range(config.trials_per_setting_pair):
            if leak_source_to_setting:
                source = sources[(trial_index + pair_index * 3) % len(sources)]
                source = SourceState(
                    source_id=f"{source.source_id}_{setting_pair}",
                    source_bucket=pair_index,
                    source_hash=source.source_hash,
                    fields={**source.fields, "leak_setting_pair": setting_pair},
                )
            else:
                source = sources[trial_index % len(sources)]
            if outcome_function is haos_local_recoverability_outcome:
                alice, alice_score = outcome_function(source, "alice", alice_setting, config)
                bob, bob_score = outcome_function(source, "bob", bob_setting, config)
            else:
                alice, alice_score = outcome_function(source, "alice", alice_setting)
                bob, bob_score = outcome_function(source, "bob", bob_setting)
            retained = True
            rejection_reason = ""
            if post_select and alice * bob != chsh_sign:
                retained = False
                rejection_reason = "post_selection_favored_chsh_sign"
            outcome = BellOutcome(
                alice=int(alice),
                bob=int(bob),
                alice_local_score=float(alice_score),
                bob_local_score=float(bob_score),
                retained=retained,
                rejection_reason=rejection_reason,
            )
            trials.append(
                BellTrial(
                    candidate_id=candidate_id,
                    run_id=run_id,
                    stage=stage,
                    seed=seed,
                    trial_index=trial_index,
                    setting_pair=setting_pair,
                    alice_setting=alice_setting,
                    bob_setting=bob_setting,
                    source_id=source.source_id,
                    source_bucket=source.source_bucket,
                    source_artifact_hash=source.source_hash,
                    outcome=outcome,
                    provenance_note=provenance_note,
                )
            )
    return trials


def b3_joint_closure_trials(
    seed: int,
    sources: list[SourceState],
    config: BellCandidateConfig,
) -> tuple[list[BellTrial], list[dict[str, Any]]]:
    candidate_id = "haos_joint_closure_candidate"
    run_id = f"{candidate_id}_seed_{seed}"
    trials: list[BellTrial] = []
    audit_rows: list[dict[str, Any]] = []
    for alice_setting, bob_setting in config.setting_pairs:
        setting_pair = pair_label(alice_setting, bob_setting)
        for half_index in range((config.trials_per_setting_pair + 1) // 2):
            first_trial_index = half_index * 2
            source = sources[half_index % len(sources)]
            product, costs = b3_sample_product(source, alice_setting, bob_setting, seed, first_trial_index, config)
            for offset, alice in enumerate((1, -1)):
                trial_index = first_trial_index + offset
                if trial_index >= config.trials_per_setting_pair:
                    break
                bob = product * alice
                outcome = BellOutcome(
                    alice=alice,
                    bob=bob,
                    alice_local_score=0.0,
                    bob_local_score=float(costs["p_product_plus"]),
                    retained=True,
                    rejection_reason="",
                )
                trials.append(
                    BellTrial(
                        candidate_id=candidate_id,
                        run_id=run_id,
                        stage="B3",
                        seed=seed,
                        trial_index=trial_index,
                        setting_pair=setting_pair,
                        alice_setting=alice_setting,
                        bob_setting=bob_setting,
                        source_id=source.source_id,
                        source_bucket=source.source_bucket,
                        source_artifact_hash=source.source_hash,
                        outcome=outcome,
                        provenance_note="B3 joint closure cost samples an outcome product before symmetric local registration",
                    )
                )
            audit_rows.append(
                {
                    "candidate_id": candidate_id,
                    "run_id": run_id,
                    "seed": seed,
                    "trial_index": first_trial_index,
                    "setting_pair": setting_pair,
                    "source_id": source.source_id,
                    "source_bucket": source.source_bucket,
                    "base_instability": f"{costs['base_instability']:.12g}",
                    "closure_strength": f"{costs['closure_strength']:.12g}",
                    "orientation_score": f"{costs['orientation_score']:.12g}",
                    "preferred_product": int(costs["preferred_product"]),
                    "cost_product_plus": f"{costs['cost_product_plus']:.12g}",
                    "cost_product_minus": f"{costs['cost_product_minus']:.12g}",
                    "p_product_plus": f"{costs['p_product_plus']:.12g}",
                    "sampled_product": product,
                    "source_features": json.dumps(
                        {
                            "shell_order_score": costs["shell_order_score"],
                            "max_shell_overlap_fraction": costs["max_shell_overlap_fraction"],
                            "phase17_mismatch_rate": costs["phase17_mismatch_rate"],
                            "arrival_gap": costs["arrival_gap"],
                        },
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                }
            )
    return trials, audit_rows


def label_permutation_control(base_trials: list[BellTrial], config: BellCandidateConfig) -> list[BellTrial]:
    mapping = {
        ("a0", "b0"): ("a0", "b1"),
        ("a0", "b1"): ("a1", "b0"),
        ("a1", "b0"): ("a1", "b1"),
        ("a1", "b1"): ("a0", "b0"),
    }
    remapped: list[BellTrial] = []
    for trial in base_trials:
        alice_setting, bob_setting = mapping[(trial.alice_setting, trial.bob_setting)]
        remapped.append(
            BellTrial(
                candidate_id="label_permutation_control",
                run_id=trial.run_id.replace(trial.candidate_id, "label_permutation_control"),
                stage="control",
                seed=trial.seed,
                trial_index=trial.trial_index,
                setting_pair=pair_label(alice_setting, bob_setting),
                alice_setting=alice_setting,
                bob_setting=bob_setting,
                source_id=trial.source_id,
                source_bucket=trial.source_bucket,
                source_artifact_hash=trial.source_artifact_hash,
                outcome=trial.outcome,
                provenance_note="control permutes labels after local outcomes are generated",
            )
        )
    return remapped


def setting_permutation_control(base_trials: list[BellTrial], config: BellCandidateConfig) -> list[BellTrial]:
    pairs = list(config.setting_pairs)
    remapped: list[BellTrial] = []
    for trial in base_trials:
        pair_index = (trial.trial_index + trial.seed) % len(pairs)
        alice_setting, bob_setting = pairs[pair_index]
        remapped.append(
            BellTrial(
                candidate_id="setting_permutation_control",
                run_id=trial.run_id.replace(trial.candidate_id, "setting_permutation_control"),
                stage="control",
                seed=trial.seed,
                trial_index=trial.trial_index,
                setting_pair=pair_label(alice_setting, bob_setting),
                alice_setting=alice_setting,
                bob_setting=bob_setting,
                source_id=trial.source_id,
                source_bucket=trial.source_bucket,
                source_artifact_hash=trial.source_artifact_hash,
                outcome=trial.outcome,
                provenance_note="control permutes setting labels independently of source generation",
            )
        )
    return remapped


def validate_binary_outcomes(trials: list[BellTrial]) -> None:
    for trial in trials:
        if trial.outcome.alice not in (-1, 1) or trial.outcome.bob not in (-1, 1):
            raise ValueError(f"non-binary outcome in {trial.run_id}")


def correlation_table(candidate_id: str, run_id: str, stage: str, sample_role: str, trials: list[BellTrial], config: BellCandidateConfig) -> BellCorrelationTable:
    retained = [trial for trial in trials if trial.outcome.retained]
    rejected = len(trials) - len(retained)
    correlations: dict[str, float] = {}
    standard_errors: dict[str, float] = {}
    trials_by_pair: dict[str, int] = {}
    alice_marginals: dict[str, float] = {}
    bob_marginals: dict[str, float] = {}
    for alice_setting, bob_setting in config.setting_pairs:
        key = pair_key(alice_setting, bob_setting)
        pair_trials = [
            trial for trial in retained
            if trial.alice_setting == alice_setting and trial.bob_setting == bob_setting
        ]
        products = [trial.outcome.alice * trial.outcome.bob for trial in pair_trials]
        n_value = len(products)
        trials_by_pair[key] = n_value
        correlations[key] = sum(products) / float(n_value) if n_value else 0.0
        if n_value > 1:
            variance = max(0.0, 1.0 - correlations[key] * correlations[key])
            standard_errors[key] = math.sqrt(variance / float(n_value))
        else:
            standard_errors[key] = 0.0
        alice_marginals[key] = sum(trial.outcome.alice for trial in pair_trials) / float(n_value) if n_value else 0.0
        bob_marginals[key] = sum(trial.outcome.bob for trial in pair_trials) / float(n_value) if n_value else 0.0
    signed = sum(CHSH_SIGNS[key] * correlations[key] for key in PAIR_KEYS)
    s_value = abs(signed)
    se_chsh = math.sqrt(sum(standard_errors[key] * standard_errors[key] for key in PAIR_KEYS))
    ci_low = max(0.0, s_value - config.ci_z * se_chsh)
    ci_high = s_value + config.ci_z * se_chsh
    return BellCorrelationTable(
        candidate_id=candidate_id,
        run_id=run_id,
        stage=stage,
        sample_role=sample_role,
        correlations=correlations,
        standard_errors=standard_errors,
        trials_by_pair=trials_by_pair,
        alice_marginals_by_pair=alice_marginals,
        bob_marginals_by_pair=bob_marginals,
        S=s_value,
        ci_low=ci_low,
        ci_high=ci_high,
        retained_trials=len(retained),
        rejected_trials=rejected,
        retention_rate=len(retained) / float(len(trials)) if trials else 0.0,
    )


def no_signalling_rows(candidate_id: str, run_id: str, trials: list[BellTrial], config: BellCandidateConfig) -> tuple[list[dict[str, Any]], float, str]:
    retained = [trial for trial in trials if trial.outcome.retained]
    rows: list[dict[str, Any]] = []
    deviations: list[float] = []
    for local_setting in ("a0", "a1"):
        buckets = []
        for remote_setting in ("b0", "b1"):
            subset = [
                trial for trial in retained
                if trial.alice_setting == local_setting and trial.bob_setting == remote_setting
            ]
            p_plus = sum(1 for trial in subset if trial.outcome.alice == 1) / float(len(subset)) if subset else 0.0
            buckets.append((remote_setting, p_plus))
        deviation = abs(buckets[0][1] - buckets[1][1])
        deviations.append(deviation)
        rows.append(
            {
                "candidate_id": candidate_id,
                "run_id": run_id,
                "party": "alice",
                "local_setting": local_setting,
                "remote_setting_left": buckets[0][0],
                "remote_setting_right": buckets[1][0],
                "p_plus_left": buckets[0][1],
                "p_plus_right": buckets[1][1],
                "deviation": deviation,
                "tolerance": config.no_signalling_tolerance,
                "status": "NO_SIGNALLING_PASS" if deviation <= config.no_signalling_tolerance else "NO_SIGNALLING_FAIL",
            }
        )
    for local_setting in ("b0", "b1"):
        buckets = []
        for remote_setting in ("a0", "a1"):
            subset = [
                trial for trial in retained
                if trial.bob_setting == local_setting and trial.alice_setting == remote_setting
            ]
            p_plus = sum(1 for trial in subset if trial.outcome.bob == 1) / float(len(subset)) if subset else 0.0
            buckets.append((remote_setting, p_plus))
        deviation = abs(buckets[0][1] - buckets[1][1])
        deviations.append(deviation)
        rows.append(
            {
                "candidate_id": candidate_id,
                "run_id": run_id,
                "party": "bob",
                "local_setting": local_setting,
                "remote_setting_left": buckets[0][0],
                "remote_setting_right": buckets[1][0],
                "p_plus_left": buckets[0][1],
                "p_plus_right": buckets[1][1],
                "deviation": deviation,
                "tolerance": config.no_signalling_tolerance,
                "status": "NO_SIGNALLING_PASS" if deviation <= config.no_signalling_tolerance else "NO_SIGNALLING_FAIL",
            }
        )
    max_deviation = max(deviations, default=0.0)
    status = "NO_SIGNALLING_PASS" if max_deviation <= config.no_signalling_tolerance else "NO_SIGNALLING_FAIL"
    return rows, max_deviation, status


def mutual_information_source_setting_bits(trials: list[BellTrial]) -> float:
    retained = [trial for trial in trials if trial.outcome.retained]
    if not retained:
        return 0.0
    total = float(len(retained))
    joint: dict[tuple[int, str], int] = {}
    source_counts: dict[int, int] = {}
    setting_counts: dict[str, int] = {}
    for trial in retained:
        source_counts[trial.source_bucket] = source_counts.get(trial.source_bucket, 0) + 1
        setting_counts[trial.setting_pair] = setting_counts.get(trial.setting_pair, 0) + 1
        key = (trial.source_bucket, trial.setting_pair)
        joint[key] = joint.get(key, 0) + 1
    mi = 0.0
    for (bucket, setting), count in joint.items():
        p_joint = count / total
        p_source = source_counts[bucket] / total
        p_setting = setting_counts[setting] / total
        mi += p_joint * math.log2(p_joint / max(p_source * p_setting, 1.0e-300))
    return max(0.0, mi)


def setting_independence_row(candidate_id: str, run_id: str, trials: list[BellTrial], config: BellCandidateConfig) -> tuple[dict[str, Any], float, str]:
    mi = mutual_information_source_setting_bits(trials)
    source_buckets = sorted({trial.source_bucket for trial in trials if trial.outcome.retained})
    setting_pairs = sorted({trial.setting_pair for trial in trials if trial.outcome.retained})
    status = "SETTING_INDEPENDENCE_PASS" if mi <= config.setting_independence_mi_tolerance_bits else "SETTING_INDEPENDENCE_FAIL"
    return (
        {
            "candidate_id": candidate_id,
            "run_id": run_id,
            "mutual_information_bits": mi,
            "tolerance_bits": config.setting_independence_mi_tolerance_bits,
            "source_bucket_count": len(source_buckets),
            "setting_pair_count": len(setting_pairs),
            "status": status,
        },
        mi,
        status,
    )


def classify_result(
    candidate_id: str,
    table: BellCorrelationTable,
    no_signalling_status: str,
    setting_independence_status: str,
    target_leakage: bool,
    post_selection: bool,
    remote_setting_access: bool,
    config: BellCandidateConfig,
) -> list[str]:
    labels: list[str] = []
    robust_violation = table.ci_low > config.local_bound
    overshoot = table.S > config.local_bound + config.finite_sample_tolerance

    if candidate_id == "local_sanity_candidate" and not robust_violation:
        labels.append("LOCAL_CANDIDATE_SANITY_PASS")
    if not robust_violation:
        labels.extend(["LOCAL_BOUND_RECOVERED", "CHSH_VIOLATION_NOT_DETECTED"])
    else:
        labels.append("CHSH_VIOLATION_DETECTED")
    if candidate_id in {"local_sanity_candidate", "haos_local_recoverability_candidate"} and robust_violation:
        labels.append("IMPLEMENTATION_LEAKAGE_SUSPECTED")
    elif candidate_id == "local_sanity_candidate" and overshoot:
        labels.append("IMPLEMENTATION_LEAKAGE_SUSPECTED")

    labels.append(no_signalling_status)
    labels.append(setting_independence_status)
    if post_selection:
        labels.append("POST_SELECTION_DETECTED")
    if target_leakage:
        labels.append("TARGET_LEAKAGE_DETECTED")
    if remote_setting_access:
        labels.append("IMPLEMENTATION_LEAKAGE_SUSPECTED")
    labels.extend(["CST_BELL_STATUS_OPEN", "HAOS_DERIVATION_NOT_ESTABLISHED", "NO_PHYSICAL_EXPERIMENT"])
    if candidate_id == "haos_local_recoverability_candidate":
        labels.append("MIXED / OPEN")
    return list(dict.fromkeys(labels))


def make_assumption_ledger(
    candidate_id: str,
    run_id: str,
    stage: str,
    no_signalling_status: str,
    setting_independence_status: str,
    target_leakage: bool,
    post_selection: bool,
    remote_setting_access: bool = False,
    joint_setting_access_declared: bool = False,
    factorization: str = "factorized_local_response",
    outcome_independence: str = "deterministic_given_source_and_local_setting",
    forward_only_causality: str = "forward_only_by_construction",
) -> BellAssumptionLedger:
    violations: list[str] = []
    if setting_independence_status != "SETTING_INDEPENDENCE_PASS":
        violations.append("MEASUREMENT_INDEPENDENCE_VIOLATED")
    if no_signalling_status != "NO_SIGNALLING_PASS":
        violations.append("NO_SIGNALLING_FAIL")
    if post_selection:
        violations.append("POST_SELECTION_DETECTED")
    if target_leakage:
        violations.append("TARGET_LEAKAGE_DETECTED")
    if remote_setting_access:
        violations.append("REMOTE_SETTING_ACCESS_DETECTED")
    return BellAssumptionLedger(
        candidate_id=candidate_id,
        run_id=run_id,
        stage=stage,
        factorization=factorization,
        setting_independence=setting_independence_status,
        parameter_independence=no_signalling_status,
        outcome_independence=outcome_independence,
        forward_only_causality=forward_only_causality,
        remote_setting_access=remote_setting_access,
        joint_setting_access_declared=joint_setting_access_declared,
        post_selection_detected=post_selection,
        target_leakage_detected=target_leakage,
        assumption_violations=violations,
    )


def evaluate_candidate(
    candidate_id: str,
    stage: str,
    sample_role: str,
    trials: list[BellTrial],
    config: BellCandidateConfig,
    target_leakage: bool = False,
    remote_setting_access: bool = False,
    joint_setting_access_declared: bool = False,
    factorization: str = "factorized_local_response",
    outcome_independence: str = "deterministic_given_source_and_local_setting",
    forward_only_causality: str = "forward_only_by_construction",
) -> BellCandidateResult:
    validate_binary_outcomes(trials)
    run_id = trials[0].run_id if trials else f"{candidate_id}_empty"
    table = correlation_table(candidate_id, run_id, stage, sample_role, trials, config)
    _ns_rows, max_ns, ns_status = no_signalling_rows(candidate_id, run_id, trials, config)
    _si_row, mi, si_status = setting_independence_row(candidate_id, run_id, trials, config)
    post_selection = table.rejected_trials > 0
    labels = classify_result(candidate_id, table, ns_status, si_status, target_leakage, post_selection, remote_setting_access, config)
    ledger = make_assumption_ledger(
        candidate_id,
        run_id,
        stage,
        ns_status,
        si_status,
        target_leakage,
        post_selection,
        remote_setting_access=remote_setting_access,
        joint_setting_access_declared=joint_setting_access_declared,
        factorization=factorization,
        outcome_independence=outcome_independence,
        forward_only_causality=forward_only_causality,
    )
    return BellCandidateResult(
        candidate_id=candidate_id,
        run_id=run_id,
        stage=stage,
        sample_role=sample_role,
        correlation_table=table,
        assumption_ledger=ledger,
        no_signalling_max_deviation=max_ns,
        setting_source_mutual_information_bits=mi,
        verdict_labels=labels,
    )


def trial_rows(trials: list[BellTrial]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for trial in trials:
        rows.append(
            {
                "candidate_id": trial.candidate_id,
                "run_id": trial.run_id,
                "stage": trial.stage,
                "seed": trial.seed,
                "trial_index": trial.trial_index,
                "setting_pair": trial.setting_pair,
                "alice_setting": trial.alice_setting,
                "bob_setting": trial.bob_setting,
                "source_id": trial.source_id,
                "source_bucket": trial.source_bucket,
                "source_artifact_hash": trial.source_artifact_hash,
                "alice_outcome": trial.outcome.alice,
                "bob_outcome": trial.outcome.bob,
                "product": trial.outcome.alice * trial.outcome.bob,
                "retained": trial.outcome.retained,
                "rejection_reason": trial.outcome.rejection_reason,
                "alice_local_score": f"{trial.outcome.alice_local_score:.12g}",
                "bob_local_score": f"{trial.outcome.bob_local_score:.12g}",
                "provenance_note": trial.provenance_note,
            }
        )
    return rows


def result_control_row(result: BellCandidateResult) -> dict[str, Any]:
    table = result.correlation_table
    return {
        "candidate_id": result.candidate_id,
        "run_id": result.run_id,
        "stage": result.stage,
        "sample_role": result.sample_role,
        "S": f"{table.S:.12g}",
        "ci_low": f"{table.ci_low:.12g}",
        "ci_high": f"{table.ci_high:.12g}",
        "E00": f"{table.correlations['E00']:.12g}",
        "E01": f"{table.correlations['E01']:.12g}",
        "E10": f"{table.correlations['E10']:.12g}",
        "E11": f"{table.correlations['E11']:.12g}",
        "retained_trials": table.retained_trials,
        "rejected_trials": table.rejected_trials,
        "retention_rate": f"{table.retention_rate:.12g}",
        "no_signalling_max_deviation": f"{result.no_signalling_max_deviation:.12g}",
        "setting_source_mutual_information_bits": f"{result.setting_source_mutual_information_bits:.12g}",
        "target_leakage_detected": result.assumption_ledger.target_leakage_detected,
        "post_selection_detected": result.assumption_ledger.post_selection_detected,
        "verdict_labels": ";".join(result.verdict_labels),
    }


def precommitment_contract(config: BellCandidateConfig, provenance: BellCandidateProvenance) -> dict[str, Any]:
    return {
        "name": "Bell HAOS-IIP candidate bridge",
        "status": "isolated_candidate_bundle",
        "scope": {
            "implemented_fact": "Build B0 adapter, B1 local sanity candidate, B2 HAOS-IIP-local recoverability candidate, and B3 explicit joint-closure-cost candidate.",
            "design_choice": "Use the frozen Bell/CHSH scoreboard convention only after candidate correlations are generated.",
            "heuristic": "B2 maps frozen Phase 18 source signatures through local perturbation-specific recoverability proxies.",
            "unverified_hypothesis": "B3 tests one non-factorizable operational candidate; it does not establish a HAOS-IIP Bell derivation.",
        },
        "non_claims": [
            "not a physical Bell experiment",
            "not a loophole-free Bell test",
            "not evidence that HAOS-IIP derives quantum mechanics",
            "not evidence that CST recovers Bell correlations",
            "not a Born-rule derivation",
            "not a license to import quantum correlations by construction",
        ],
        "frozen_reference_policy": {
            "reference_sidecar_path": str(REFERENCE_SIDECAR.relative_to(REPO_ROOT)),
            "used_as": "post_generation_scoreboard_convention_only",
            "candidate_imports_reference_runner": False,
            "candidate_reads_quantum_reference_table": False,
            "reference_hashes": provenance.reference_sidecar_hashes,
        },
        "settings": {
            "labels": {
                "alice": ["a0", "a1"],
                "bob": ["b0", "b1"],
            },
            "setting_pairs": [list(pair) for pair in config.setting_pairs],
            "generation_method": "balanced_complete_setting_grid_per_seed_independent_of_source_state",
        },
        "source_artifacts": {
            "B1": "synthetic balanced local lambda; no HAOS claim",
            "B2": config.b2_source_artifact,
            "B2_filter": {"hierarchy_label": config.b2_hierarchy_label, "source_limit": config.b2_source_limit},
            "source_artifact_hash": provenance.source_artifact_hash,
        },
        "candidate_stages": {
            "B0": "neutral adapter records and diagnostics",
            "B1": "local sanity candidate A=f(lambda,a), B=g(lambda,b)",
            "B2": "HAOS-IIP-local recoverability candidate A=f(lambda,a), B=g(lambda,b)",
            "B3": "joint closure candidate P(A,B|a,b,lambda) from an explicit Phase 18 chain-order closure cost",
        },
        "forbidden_candidate_access": [
            "2sqrt2 target constants",
            "cosine quantum correlation curves",
            "frozen quantum reference correlations",
            "target S values",
            "reference outcomes",
            "adaptive optimization against CHSH score",
            "outcome-dependent trial post-selection",
        ],
        "precommitted_B2_mapping": {
            "lambda": "shared frozen Phase 18 source signature row",
            "local_setting": "party-local perturbation label a0/a1 or b0/b1",
            "recoverability_proxy": "shell_order_score - 0.55*max_shell_overlap_fraction - 0.75*phase17_mismatch_rate + 0.10*arrival_gap plus local fixed perturbation bias",
            "binary_rule": "outcome is +1 when local recoverability proxy exceeds the party/setting threshold; otherwise -1",
            "tie_behavior": "stable hash sign from source_id, party, and setting",
            "thresholds": {
                "a0": config.b2_threshold_a0,
                "a1": config.b2_threshold_a1,
                "b0": config.b2_threshold_b0,
                "b1": config.b2_threshold_b1,
            },
        },
        "precommitted_B3_mapping": {
            "lambda": "shared frozen Phase 18 source signature row",
            "joint_probability_model": "P(A,B | a,b,lambda) = 0.5 * P(product=A*B | a,b,lambda)",
            "joint_cost": "P(product=q | a,b,lambda) proportional to exp(-beta * C(q,a,b,lambda))",
            "cost_definition": "C(q,a,b,lambda) = base_instability(lambda) + 0.5 * (1 - q * preferred_product(a,b,lambda) * closure_strength(lambda))",
            "closure_strength": "clamped shell_order_score*(1-max_shell_overlap_fraction)*(1-phase17_mismatch_rate) + 0.25*arrival_gap",
            "preferred_product": "sign of Phase 18 chain-order displacement between the declared Alice and Bob setting token pairs",
            "local_registration": "after product sampling, Alice is assigned symmetric +/- registrations and Bob receives product*Alice",
            "tie_behavior": "zero orientation prefers +1 by predeclared convention",
            "beta": config.b3_beta,
            "strength_floor": config.b3_strength_floor,
            "strength_ceiling": config.b3_strength_ceiling,
            "invoice": ["FACTORIZATION_RELAXED", "OUTCOME_INDEPENDENCE_RELAXED"],
        },
        "trial_policy": {
            "seeds": list(config.seeds),
            "trials_per_setting_pair": config.trials_per_setting_pair,
            "all_valid_trials_retained": True,
            "post_selection_policy": "detected and invalidated; positive control intentionally violates this",
        },
        "audit_thresholds": {
            "local_bound": config.local_bound,
            "ci_z": config.ci_z,
            "finite_sample_tolerance": config.finite_sample_tolerance,
            "no_signalling_tolerance": config.no_signalling_tolerance,
            "setting_independence_mi_tolerance_bits": config.setting_independence_mi_tolerance_bits,
        },
        "verdict_logic": {
            "structural_interest_requires": [
                "CI lower bound excludes 2",
                "NO_SIGNALLING_PASS",
                "SETTING_INDEPENDENCE_PASS or declared relaxation",
                "no target leakage",
                "no post-selection",
                "reproducibility across predeclared seeds",
                "matched local controls remain at or below the local bound",
            ],
            "B1_expected": [
                "LOCAL_CANDIDATE_SANITY_PASS",
                "LOCAL_BOUND_RECOVERED",
                "CHSH_VIOLATION_NOT_DETECTED",
            ],
            "B2_expected": [
                "LOCAL_BOUND_RECOVERED",
                "CHSH_VIOLATION_NOT_DETECTED",
                "CST_BELL_STATUS_OPEN",
            ],
        },
        "provenance": asdict(provenance),
    }


def b3_precommitment_contract(config: BellCandidateConfig, provenance: BellCandidateProvenance) -> dict[str, Any]:
    return {
        "name": "B3_PRECOMMITMENT_CONTRACT",
        "status": "AUTHORIZED_PRECOMMITMENT",
        "authorization_status": "IMPLEMENTATION_AUTHORIZED_WITH_EXPLICIT_INVOICE",
        "implementation_authorized": True,
        "created_after": {
            "B0_adapter": "PASS",
            "B1_local_sanity": "PASS",
            "B2_haos_local_recoverability": "BELL_LOCAL",
            "prior_B3_gate": "NONE_IDENTIFIED_IMPLEMENTATION_NOT_AUTHORIZED",
        },
        "scope": {
            "implemented_fact": "This contract authorizes one explicit B3 joint-closure-cost candidate before execution.",
            "design_choice": "Use frozen Phase 18 source rows and setting-token chain-order displacements to define a joint outcome-product cost.",
            "heuristic": "The cost is a HAOS-native engineering candidate, not a physical law and not a derivation of quantum mechanics.",
            "unverified_hypothesis": "The candidate may or may not produce Bell-nonlocal correlations; no positive result is assumed.",
        },
        "joint_probability_model": {
            "required_form": "P(A,B | a,b,lambda)",
            "provided": True,
            "binary_outcomes": [-1, 1],
            "required_setting_pairs": [list(pair) for pair in config.setting_pairs],
            "construction": {
                "product_distribution": "P(Q=q | a,b,lambda) = exp(-beta*C(q,a,b,lambda)) / sum_{r in {-1,+1}} exp(-beta*C(r,a,b,lambda))",
                "local_registration": "P(A=x,B=y | a,b,lambda) = 0.5 * P(Q=x*y | a,b,lambda)",
                "cost": "C(q,a,b,lambda) = base_instability(lambda) + 0.5*(1 - q*preferred_product(a,b,lambda)*closure_strength(lambda))",
                "base_instability": "clamp((1-shell_order_score)+max_shell_overlap_fraction+phase17_mismatch_rate+(1-arrival_gap),0,3)",
                "closure_strength": "clamp(shell_order_score*(1-max_shell_overlap_fraction)*(1-phase17_mismatch_rate)+0.25*arrival_gap,strength_floor,strength_ceiling)",
                "preferred_product": "sign of chain-order displacement between setting token pairs",
                "beta": config.b3_beta,
                "strength_floor": config.b3_strength_floor,
                "strength_ceiling": config.b3_strength_ceiling,
            },
        },
        "non_factorizability_statement": {
            "provided": True,
            "statement": "For nonzero closure_strength, P(A,B | a,b,lambda) generally differs from P(A | a,lambda)P(B | b,lambda) because the product Q=A*B is sampled from one joint distribution.",
            "factorization_failure_mechanism": "outcome-product coupling under joint closure cost",
            "not_claimed": [
                "not quantum mechanics",
                "not a Born-rule derivation",
                "not a loophole-free Bell experiment",
                "not evidence of physical nonlocality",
            ],
        },
        "bell_assumption_invoice": {
            "selected_invoice": ["FACTORIZATION_RELAXED", "OUTCOME_INDEPENDENCE_RELAXED"],
            "factorization": "RELAXED_BY_EXPLICIT_JOINT_OUTCOME_PRODUCT_DISTRIBUTION",
            "measurement_independence": "PRESERVED_AND_TESTED_BY_SOURCE_SETTING_MI",
            "parameter_independence": "NOT_RELAXED; TESTED_BY_NO_SIGNALLING_MARGINALS",
            "outcome_independence": "RELAXED_BY_PRODUCT_LEVEL_JOINT_DISTRIBUTION",
            "forward_only_causality": "NO_SPACETIME_CAUSAL_CLAIM; COMPUTATIONAL_BATCH_CANDIDATE_AFTER_SETTINGS_ASSIGNED",
        },
        "observable_consequences_predeclared": {
            "predicted_chsh_range": "unknown before execution; candidate can return S<=2",
            "marginal_behavior": "Alice and Bob marginals expected near zero by symmetric local registration",
            "no_signalling_behavior": "expected to pass finite-sample tolerance if implementation has no parameter leakage",
            "seed_stability": "reported across predeclared seeds; no seed retuning",
            "control_separation": "invalid leakage and post-selection controls must still be detected",
        },
        "falsification_conditions": [
            "NO_SIGNALLING_FAIL invalidates a standard Bell-nonlocal interpretation",
            "SETTING_INDEPENDENCE_FAIL means measurement independence is not preserved",
            "POST_SELECTION_DETECTED disqualifies the run",
            "TARGET_LEAKAGE_DETECTED disqualifies the run",
            "S<=2 or CI overlapping 2 means CHSH violation is not detected",
            "matched local controls exceeding the bound robustly indicates implementation leakage or scoreboard failure",
        ],
        "prohibited_shortcuts": [
            "no quantum target curve import",
            "no target constant import",
            "no frozen quantum reference correlations",
            "no reference outcomes",
            "no post-selection",
            "no setting-source leakage",
            "no adaptive tuning against S",
            "no hidden remote-setting access in B1/B2 local paths",
        ],
        "source_artifacts": {
            "B3": config.b2_source_artifact,
            "source_artifact_hash": provenance.source_artifact_hash,
        },
        "trial_policy": {
            "seeds": list(config.seeds),
            "trials_per_setting_pair": config.trials_per_setting_pair,
            "all_valid_trials_retained": True,
            "scoreboard_blinding": "B3 outcomes are generated before CHSH classification and without reading the frozen quantum reference table.",
        },
        "provenance": asdict(provenance),
    }


def build_provenance(config: BellCandidateConfig) -> BellCandidateProvenance:
    config_hash = stable_json_hash("bell_candidate_config", asdict(config))
    return BellCandidateProvenance(
        code_path=str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        code_hash=code_hash(),
        reference_sidecar_path=str(REFERENCE_SIDECAR.relative_to(REPO_ROOT)),
        reference_sidecar_hashes=reference_sidecar_hashes(),
        source_artifact=config.b2_source_artifact,
        source_artifact_hash=hash_file(REPO_ROOT / config.b2_source_artifact),
        config_hash=config_hash,
    )


def build_all_trials(config: BellCandidateConfig) -> tuple[list[BellTrial], dict[str, list[BellTrial]], dict[str, bool], list[dict[str, Any]]]:
    b2_sources = read_phase18_sources(config)
    all_trials: list[BellTrial] = []
    grouped: dict[str, list[BellTrial]] = {}
    b3_cost_rows: list[dict[str, Any]] = []

    leakage_hits: dict[str, bool] = {}
    generators = {
        "local_sanity_candidate": local_sanity_outcome,
        "haos_local_recoverability_candidate": haos_local_recoverability_outcome,
        "independent_random_control": independent_random_outcome,
        "deterministic_lhv_control": deterministic_lhv_outcome,
        "shared_randomness_local_control": shared_randomness_outcome,
    }
    for name, function in generators.items():
        leakage_hits[name] = detect_target_leakage(generator_sources([function]))[0]
    leakage_hits["haos_joint_closure_candidate"] = detect_target_leakage(
        generator_sources([
            b3_source_features,
            b3_orientation_score,
            b3_joint_closure_costs,
            b3_sample_product,
            b3_joint_closure_trials,
        ])
    )[0]

    for seed in config.seeds:
        synthetic_sources = [synthetic_source_state(seed, trial_index, config) for trial_index in range(config.trials_per_setting_pair)]
        b1_trials = balanced_trials(
            "local_sanity_candidate",
            "B1",
            seed,
            synthetic_sources,
            config,
            local_sanity_outcome,
            "B1 deterministic local code path; no remote setting access",
        )
        grouped.setdefault("local_sanity_candidate", []).extend(b1_trials)
        all_trials.extend(b1_trials)

        b2_trials = balanced_trials(
            "haos_local_recoverability_candidate",
            "B2",
            seed,
            b2_sources,
            config,
            haos_local_recoverability_outcome,
            "B2 local recoverability proxy from frozen Phase 18 source rows",
        )
        grouped.setdefault("haos_local_recoverability_candidate", []).extend(b2_trials)
        all_trials.extend(b2_trials)

        b3_trials, seed_b3_cost_rows = b3_joint_closure_trials(seed, b2_sources, config)
        grouped.setdefault("haos_joint_closure_candidate", []).extend(b3_trials)
        all_trials.extend(b3_trials)
        b3_cost_rows.extend(seed_b3_cost_rows)

        for candidate_id, function in [
            ("independent_random_control", independent_random_outcome),
            ("deterministic_lhv_control", deterministic_lhv_outcome),
            ("shared_randomness_local_control", shared_randomness_outcome),
        ]:
            trials = balanced_trials(
                candidate_id,
                "control",
                seed,
                synthetic_sources,
                config,
                function,
                f"{candidate_id} generated locally without target correlations",
            )
            grouped.setdefault(candidate_id, []).extend(trials)
            all_trials.extend(trials)

        leakage_trials = balanced_trials(
            "source_setting_leakage_positive_control",
            "invalid_control",
            seed,
            b2_sources,
            config,
            haos_local_recoverability_outcome,
            "positive control deliberately makes source bucket depend on setting pair",
            leak_source_to_setting=True,
        )
        grouped.setdefault("source_setting_leakage_positive_control", []).extend(leakage_trials)
        all_trials.extend(leakage_trials)

        post_selection_trials = balanced_trials(
            "post_selection_positive_control",
            "invalid_control",
            seed,
            synthetic_sources,
            config,
            deterministic_lhv_outcome,
            "positive control deliberately removes outcomes unfavorable to CHSH signs",
            post_select=True,
        )
        grouped.setdefault("post_selection_positive_control", []).extend(post_selection_trials)
        all_trials.extend(post_selection_trials)

    label_trials = label_permutation_control(grouped["haos_local_recoverability_candidate"], config)
    setting_trials = setting_permutation_control(grouped["haos_local_recoverability_candidate"], config)
    grouped["label_permutation_control"] = label_trials
    grouped["setting_permutation_control"] = setting_trials
    all_trials.extend(label_trials)
    all_trials.extend(setting_trials)
    leakage_hits["label_permutation_control"] = False
    leakage_hits["setting_permutation_control"] = False
    leakage_hits["source_setting_leakage_positive_control"] = leakage_hits["haos_local_recoverability_candidate"]
    leakage_hits["post_selection_positive_control"] = leakage_hits["deterministic_lhv_control"]
    return all_trials, grouped, leakage_hits, b3_cost_rows


def evaluate_all(grouped_trials: dict[str, list[BellTrial]], leakage_hits: dict[str, bool], config: BellCandidateConfig) -> list[BellCandidateResult]:
    results: list[BellCandidateResult] = []
    sample_roles = {
        "local_sanity_candidate": "candidate",
        "haos_local_recoverability_candidate": "candidate",
        "haos_joint_closure_candidate": "candidate",
        "independent_random_control": "control",
        "deterministic_lhv_control": "control",
        "shared_randomness_local_control": "control",
        "label_permutation_control": "control",
        "setting_permutation_control": "control",
        "source_setting_leakage_positive_control": "invalid_positive_control",
        "post_selection_positive_control": "invalid_positive_control",
    }
    for candidate_id, trials in grouped_trials.items():
        stage = trials[0].stage if trials else "unknown"
        factorization = "factorized_local_response"
        outcome_independence = "deterministic_given_source_and_local_setting"
        forward_only = "forward_only_by_construction"
        joint_setting_access_declared = False
        if candidate_id == "haos_joint_closure_candidate":
            factorization = "FACTORIZATION_RELAXED_BY_EXPLICIT_JOINT_CLOSURE_COST"
            outcome_independence = "OUTCOME_INDEPENDENCE_RELAXED_BY_PRODUCT_LEVEL_JOINT_DISTRIBUTION"
            forward_only = "COMPUTATIONAL_BATCH_MODEL_AFTER_SETTINGS_ASSIGNED_NO_SPACETIME_CAUSAL_CLAIM"
            joint_setting_access_declared = True
        results.append(
            evaluate_candidate(
                candidate_id,
                stage,
                sample_roles[candidate_id],
                trials,
                config,
                target_leakage=leakage_hits.get(candidate_id, False),
                joint_setting_access_declared=joint_setting_access_declared,
                factorization=factorization,
                outcome_independence=outcome_independence,
                forward_only_causality=forward_only,
            )
        )
    return results


def aggregate_no_signalling_rows(grouped_trials: dict[str, list[BellTrial]], config: BellCandidateConfig) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for candidate_id, trials in grouped_trials.items():
        candidate_rows, _max_ns, _status = no_signalling_rows(candidate_id, trials[0].run_id, trials, config)
        rows.extend(candidate_rows)
    return rows


def aggregate_setting_independence_rows(grouped_trials: dict[str, list[BellTrial]], config: BellCandidateConfig) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for candidate_id, trials in grouped_trials.items():
        row, _mi, _status = setting_independence_row(candidate_id, trials[0].run_id, trials, config)
        rows.append(row)
    return rows


def result_labels(results: list[BellCandidateResult]) -> list[str]:
    b1 = next(result for result in results if result.candidate_id == "local_sanity_candidate")
    b2 = next(result for result in results if result.candidate_id == "haos_local_recoverability_candidate")
    b3 = next(result for result in results if result.candidate_id == "haos_joint_closure_candidate")
    labels: list[str] = []
    labels.extend(b1.verdict_labels)
    labels.extend(label for label in b2.verdict_labels if label not in labels)
    labels.extend(label for label in b3.verdict_labels if label not in labels)
    invalid_positive = [result for result in results if result.sample_role == "invalid_positive_control"]
    if any(result.assumption_ledger.post_selection_detected for result in invalid_positive):
        labels.append("POST_SELECTION_DETECTED")
    if any(result.assumption_ledger.setting_independence == "SETTING_INDEPENDENCE_FAIL" for result in invalid_positive):
        labels.append("SETTING_INDEPENDENCE_FAIL")
    return list(dict.fromkeys(labels))


def seed_variance(results: list[BellCandidateResult], grouped_trials: dict[str, list[BellTrial]], config: BellCandidateConfig) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for candidate_id in ("local_sanity_candidate", "haos_local_recoverability_candidate", "haos_joint_closure_candidate"):
        s_values: list[float] = []
        for seed in config.seeds:
            seed_trials = [trial for trial in grouped_trials[candidate_id] if trial.seed == seed]
            table = correlation_table(candidate_id, f"{candidate_id}_seed_{seed}", grouped_trials[candidate_id][0].stage, "candidate", seed_trials, config)
            s_values.append(table.S)
        mean_s = sum(s_values) / float(len(s_values))
        variance = sum((value - mean_s) ** 2 for value in s_values) / float(len(s_values)) if s_values else 0.0
        output[candidate_id] = {
            "seeds": list(config.seeds),
            "S_by_seed": s_values,
            "mean_S": mean_s,
            "variance_S": variance,
        }
    return output


def write_report(path: Path, result: dict[str, Any], candidate_results: list[BellCandidateResult]) -> None:
    b1 = next(item for item in candidate_results if item.candidate_id == "local_sanity_candidate")
    b2 = next(item for item in candidate_results if item.candidate_id == "haos_local_recoverability_candidate")
    b3 = next(item for item in candidate_results if item.candidate_id == "haos_joint_closure_candidate")
    controls = [item for item in candidate_results if item.sample_role in {"control", "invalid_positive_control"}]
    lines = [
        "# Bell HAOS-IIP Candidate Bridge",
        "",
        "Implemented fact: B0 adapter, B1 local sanity candidate, B2 HAOS-IIP-local recoverability candidate, and B3 joint-closure-cost candidate were executed.",
        "Design choice: the frozen CHSH sidecar remains a scoreboard convention only; candidate generation does not import its quantum table.",
        "Heuristic: B2 uses frozen Phase 18 source rows and local recoverability proxy thresholds.",
        "Heuristic: B3 samples a joint outcome product from a precommitted Phase 18 chain-order closure cost.",
        "Unverified hypothesis: no HAOS-IIP Bell derivation is established by this run.",
        "",
        "## Verdict Labels",
    ]
    lines.extend(f"- {label}" for label in result["labels"])
    lines.extend(
        [
            "",
            "## B1 Local Sanity Candidate",
            f"- S: `{b1.correlation_table.S:.12f}`",
            f"- CI: `{b1.correlation_table.ci_low:.12f}` to `{b1.correlation_table.ci_high:.12f}`",
            f"- labels: `{'; '.join(b1.verdict_labels)}`",
            "",
            "## B2 HAOS-IIP Local Recoverability Candidate",
            f"- S: `{b2.correlation_table.S:.12f}`",
            f"- CI: `{b2.correlation_table.ci_low:.12f}` to `{b2.correlation_table.ci_high:.12f}`",
            f"- E00/E01/E10/E11: `{b2.correlation_table.correlations['E00']:.12f}`, `{b2.correlation_table.correlations['E01']:.12f}`, `{b2.correlation_table.correlations['E10']:.12f}`, `{b2.correlation_table.correlations['E11']:.12f}`",
            f"- no-signalling max deviation: `{b2.no_signalling_max_deviation:.12f}`",
            f"- source-setting MI bits: `{b2.setting_source_mutual_information_bits:.12f}`",
            f"- labels: `{'; '.join(b2.verdict_labels)}`",
            "",
            "## B3 HAOS-IIP Joint Closure Candidate",
            f"- S: `{b3.correlation_table.S:.12f}`",
            f"- CI: `{b3.correlation_table.ci_low:.12f}` to `{b3.correlation_table.ci_high:.12f}`",
            f"- E00/E01/E10/E11: `{b3.correlation_table.correlations['E00']:.12f}`, `{b3.correlation_table.correlations['E01']:.12f}`, `{b3.correlation_table.correlations['E10']:.12f}`, `{b3.correlation_table.correlations['E11']:.12f}`",
            f"- no-signalling max deviation: `{b3.no_signalling_max_deviation:.12f}`",
            f"- source-setting MI bits: `{b3.setting_source_mutual_information_bits:.12f}`",
            f"- assumption invoice: `{b3.assumption_ledger.factorization}; {b3.assumption_ledger.outcome_independence}`",
            f"- labels: `{'; '.join(b3.verdict_labels)}`",
            "",
            "## Controls",
        ]
    )
    for control in controls:
        lines.append(
            "- {candidate}: S `{s:.12f}`, CI `{lo:.12f}` to `{hi:.12f}`, retained `{retained}`, rejected `{rejected}`, labels `{labels}`".format(
                candidate=control.candidate_id,
                s=control.correlation_table.S,
                lo=control.correlation_table.ci_low,
                hi=control.correlation_table.ci_high,
                retained=control.correlation_table.retained_trials,
                rejected=control.correlation_table.rejected_trials,
                labels="; ".join(control.verdict_labels),
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "- This is not a physical Bell experiment.",
            "- This does not derive quantum mechanics, the Born rule, or Bell correlations from CST or HAOS-IIP.",
            "- B1 and B2 are factorized local models by construction.",
            "- B3 is a computational joint distribution candidate with factorization and outcome independence explicitly relaxed.",
            "- B3 makes no physical nonlocality claim and no spacetime causal claim.",
            "- A clean S <= 2 result is valid and expected for B1/B2, and remains a valid B3 outcome.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_b3_precommitment_report(path: Path, contract: dict[str, Any]) -> None:
    invoice = contract["bell_assumption_invoice"]
    model = contract["joint_probability_model"]["construction"]
    lines = [
        "# B3 Precommitment Contract",
        "",
        "Implemented fact: this report is generated before B3 scoring from the same contract used to authorize execution.",
        "Design choice: B3 tests one explicit joint-closure-cost candidate rather than generic global closure language.",
        "Heuristic: the cost uses frozen Phase 18 closure features and chain-order token displacements.",
        "Unverified hypothesis: this may or may not produce CHSH violation; no HAOS Bell derivation is established.",
        "",
        "## Authorization",
        f"- status: `{contract['status']}`",
        f"- authorization: `{contract['authorization_status']}`",
        f"- implementation authorized: `{contract['implementation_authorized']}`",
        f"- contract hash: `{contract['contract_hash']}`",
        "",
        "## Joint Probability Model",
        f"- product distribution: `{model['product_distribution']}`",
        f"- local registration: `{model['local_registration']}`",
        f"- cost: `{model['cost']}`",
        f"- preferred product: `{model['preferred_product']}`",
        "",
        "## Assumption Invoice",
        f"- selected invoice: `{'; '.join(invoice['selected_invoice'])}`",
        f"- factorization: `{invoice['factorization']}`",
        f"- measurement independence: `{invoice['measurement_independence']}`",
        f"- parameter independence: `{invoice['parameter_independence']}`",
        f"- outcome independence: `{invoice['outcome_independence']}`",
        f"- forward-only causality: `{invoice['forward_only_causality']}`",
        "",
        "## Falsification Conditions",
    ]
    lines.extend(f"- {item}" for item in contract["falsification_conditions"])
    lines.extend(
        [
            "",
            "## Boundary",
            "- This is not a physical Bell experiment.",
            "- This is not a loophole-free Bell test.",
            "- This is not a derivation of quantum mechanics or the Born rule.",
            "- A clean `S <= 2` remains a valid result.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_bridge(config: BellCandidateConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = build_provenance(config)
    contract = precommitment_contract(config, provenance)
    contract["contract_hash"] = stable_json_hash("bell_haos_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)
    b3_contract = b3_precommitment_contract(config, provenance)
    b3_contract["contract_hash"] = stable_json_hash("b3_precommitment_contract", b3_contract)
    write_json(output_dir / B3_CONTRACT_PATH.name, b3_contract)
    write_b3_precommitment_report(output_dir / B3_REPORT_PATH.name, b3_contract)

    all_trials, grouped_trials, leakage_hits, b3_cost_rows = build_all_trials(config)
    candidate_results = evaluate_all(grouped_trials, leakage_hits, config)
    assumption_ledgers = [asdict(result.assumption_ledger) for result in candidate_results]
    write_json(
        output_dir / ASSUMPTION_LEDGER_PATH.name,
        {
            "ledgers": assumption_ledgers,
            "B3": {
                "status": "EXECUTED_UNDER_PRECOMMITMENT",
                "contract_hash": b3_contract["contract_hash"],
                "invoice": ["FACTORIZATION_RELAXED", "OUTCOME_INDEPENDENCE_RELAXED"],
            },
        },
    )
    write_csv(output_dir / TRIALS_PATH.name, trial_rows(all_trials), TRIAL_FIELDNAMES)
    write_csv(output_dir / B3_COST_AUDIT_PATH.name, b3_cost_rows, B3_COST_AUDIT_FIELDNAMES)
    write_csv(output_dir / NO_SIGNALLING_PATH.name, aggregate_no_signalling_rows(grouped_trials, config), NO_SIGNALLING_FIELDNAMES)
    write_csv(output_dir / SETTING_INDEPENDENCE_PATH.name, aggregate_setting_independence_rows(grouped_trials, config), SETTING_INDEPENDENCE_FIELDNAMES)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, [result_control_row(result) for result in candidate_results], CONTROL_FIELDNAMES)

    result: dict[str, Any] = {
        "labels": result_labels(candidate_results),
        "contract_hash": contract["contract_hash"],
        "b3_contract_hash": b3_contract["contract_hash"],
        "config": asdict(config),
        "provenance": asdict(provenance),
        "candidate_results": [asdict(result) for result in candidate_results],
        "seed_variance": seed_variance(candidate_results, grouped_trials, config),
        "B3": {
            "status": "EXECUTED",
            "contract_hash": b3_contract["contract_hash"],
            "invoice": ["FACTORIZATION_RELAXED", "OUTCOME_INDEPENDENCE_RELAXED"],
            "physical_claim": "NONE",
            "derivation_claim": "NOT_ESTABLISHED",
        },
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "b3_precommitment_contract": repo_rel(output_dir / B3_CONTRACT_PATH.name),
            "b3_precommitment_report": repo_rel(output_dir / B3_REPORT_PATH.name),
            "assumption_ledger": repo_rel(output_dir / ASSUMPTION_LEDGER_PATH.name),
            "candidate_trials": repo_rel(output_dir / TRIALS_PATH.name),
            "candidate_correlations": repo_rel(output_dir / CORRELATIONS_PATH.name),
            "b3_joint_cost_audit": repo_rel(output_dir / B3_COST_AUDIT_PATH.name),
            "no_signalling_diagnostics": repo_rel(output_dir / NO_SIGNALLING_PATH.name),
            "setting_independence_diagnostics": repo_rel(output_dir / SETTING_INDEPENDENCE_PATH.name),
            "control_results": repo_rel(output_dir / CONTROL_RESULTS_PATH.name),
            "bell_candidate_report": repo_rel(output_dir / REPORT_PATH.name),
            "bell_candidate_result": repo_rel(output_dir / RESULT_PATH.name),
        },
    }
    hash_payload = {key: value for key, value in result.items() if key != "outputs"}
    result["result_hash"] = stable_json_hash("bell_haos_candidate_result", hash_payload)
    write_json(output_dir / CORRELATIONS_PATH.name, {
        "scoreboard_convention": "S = |E00 + E01 + E10 - E11|",
        "candidate_results": [asdict(result) for result in candidate_results],
    })
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, candidate_results)
    return result


def main() -> None:
    args = parse_args()
    result = run_bridge(make_config(args), args.output_dir)
    print(json.dumps({"labels": result["labels"], "result_hash": result["result_hash"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
