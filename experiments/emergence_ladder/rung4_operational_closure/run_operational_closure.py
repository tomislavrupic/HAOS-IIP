#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
CONTRACT = ROOT / "precommitment_contract.json"
RESULTS = ROOT / "results"


@dataclass(frozen=True)
class EventRun:
    run_id: str
    hierarchy: str
    n_side: int
    ensemble_size: int
    seed: int
    probe: str
    chain: tuple[str, ...]


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_runs(path: Path) -> list[EventRun]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    grouped: dict[tuple[str, int, int, int, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row["hierarchy_label"], int(row["n_side"]), int(row["ensemble_size"]), int(row["seed"]), row["probe_name"])
        grouped[key].append(row)
    runs = []
    for key, members in sorted(grouped.items()):
        signatures = {row["chain_signature"] for row in members}
        if len(signatures) != 1:
            raise ValueError(f"conflicting chain signatures for {key}")
        ranked = tuple(row["event_label"] for row in sorted(members, key=lambda row: int(row["event_rank"])))
        declared = tuple(next(iter(signatures)).split(">"))
        if ranked != declared or len(set(ranked)) != len(ranked):
            raise ValueError(f"event ranks disagree with chain signature for {key}")
        run_id = stable_hash("event_run", key)
        runs.append(EventRun(run_id, key[0], key[1], key[2], key[3], key[4], ranked))
    return runs


def train_counts(runs: Iterable[EventRun]) -> dict[str, dict[tuple[Any, ...], Counter[str]]]:
    internal: dict[tuple[Any, ...], Counter[str]] = defaultdict(Counter)
    external: dict[tuple[Any, ...], Counter[str]] = defaultdict(Counter)
    position: dict[tuple[Any, ...], Counter[str]] = defaultdict(Counter)
    marginal: Counter[str] = Counter()
    for run in runs:
        for index in range(len(run.chain) - 1):
            current = run.chain[index]
            target = run.chain[index + 1]
            internal[(run.hierarchy, run.probe, current)][target] += 1
            external[(run.hierarchy, run.probe, index)][target] += 1
            position[(index,)][target] += 1
            marginal[target] += 1
    return {"internal": internal, "external": external, "position": position, "marginal": {(): marginal}}


def choose(counter: Counter[str], remaining: tuple[str, ...]) -> str:
    return min(remaining, key=lambda label: (-counter.get(label, 0), label))


def prediction_for(model: str, counts: dict[str, dict[tuple[Any, ...], Counter[str]]], run: EventRun, index: int, *, ablate_context: bool = False) -> str:
    remaining = tuple(label for label in run.chain if label not in run.chain[: index + 1])
    if not remaining:
        raise ValueError("no candidate labels remain")
    if model == "internal":
        key = (run.hierarchy, run.probe, run.chain[index]) if not ablate_context else ("ABLATE",)
        counter = counts["internal"].get(key, counts["marginal"][()])
    elif model == "external":
        counter = counts["external"].get((run.hierarchy, run.probe, index), counts["marginal"][()])
    elif model == "position":
        counter = counts["position"].get((index,), counts["marginal"][()])
    elif model == "marginal":
        counter = counts["marginal"][()]
    else:
        raise ValueError(f"unknown model: {model}")
    return choose(counter, remaining)


def score_runs(runs: Iterable[EventRun], counts: dict[str, dict[tuple[Any, ...], Counter[str]]], *, ablate_context: bool = False, permute_targets: bool = False, seed: int = 4402) -> tuple[list[dict[str, Any]], dict[str, float]]:
    run_list = list(runs)
    rng = np.random.default_rng(seed)
    target_pool = [run.chain[index + 1] for run in run_list for index in range(len(run.chain) - 1)]
    if permute_targets:
        target_pool = list(np.array(target_pool, dtype=object)[rng.permutation(len(target_pool))])
    target_cursor = 0
    rows = []
    correct = Counter()
    total = Counter()
    for run in run_list:
        for index in range(len(run.chain) - 1):
            target = target_pool[target_cursor] if permute_targets else run.chain[index + 1]
            target_cursor += 1
            row = {
                "run_id": run.run_id,
                "hierarchy": run.hierarchy,
                "n_side": run.n_side,
                "probe": run.probe,
                "position": index,
                "current_event": run.chain[index],
                "target_event": target,
            }
            for model in ("internal", "external", "position", "marginal"):
                prediction = prediction_for(model, counts, run, index, ablate_context=ablate_context and model == "internal")
                row[f"{model}_prediction"] = prediction
                hit = int(prediction == target)
                row[f"{model}_correct"] = hit
                correct[model] += hit
                total[model] += 1
            rows.append(row)
    metrics = {f"{model}_accuracy": correct[model] / total[model] for model in total}
    metrics["best_external_accuracy"] = max(metrics["external_accuracy"], metrics["position_accuracy"], metrics["marginal_accuracy"])
    metrics["internal_advantage"] = metrics["internal_accuracy"] - metrics["best_external_accuracy"]
    return rows, metrics


def shuffled_training_runs(runs: Iterable[EventRun], seed: int) -> list[EventRun]:
    rng = np.random.default_rng(seed)
    shuffled = []
    changed = False
    for run in runs:
        chain = tuple(np.array(run.chain, dtype=object)[rng.permutation(len(run.chain))])
        changed = changed or chain != run.chain
        shuffled.append(EventRun(run.run_id, run.hierarchy, run.n_side, run.ensemble_size, run.seed, run.probe, chain))
    if not changed:
        raise ValueError("transition shuffle did not alter any chain")
    return shuffled


def run_level_advantages(rows: list[dict[str, Any]]) -> list[float]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row["run_id"]].append(row)
    advantages = []
    for members in grouped.values():
        internal = float(np.mean([row["internal_correct"] for row in members]))
        external_scores = [float(np.mean([row[f"{name}_correct"] for row in members])) for name in ("external", "position", "marginal")]
        advantages.append(internal - max(external_scores))
    return advantages


def bootstrap_ci(values: list[float], seed: int, resamples: int) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    if array.size < 2:
        raise ValueError("at least two run-level values are required")
    rng = np.random.default_rng(seed)
    means = np.empty(resamples, dtype=float)
    for index in range(resamples):
        means[index] = float(np.mean(rng.choice(array, size=len(array), replace=True)))
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def stratum_positive_fraction(rows: list[dict[str, Any]]) -> float:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(row["hierarchy"], row["probe"])].append(row)
    positive = 0
    for members in grouped.values():
        internal = np.mean([row["internal_correct"] for row in members])
        best_external = max(np.mean([row[f"{name}_correct"] for row in members]) for name in ("external", "position", "marginal"))
        positive += int(internal > best_external)
    return positive / len(grouped)


def classify(gates: dict[str, bool], instrument_valid: bool = True) -> str:
    if not instrument_valid:
        return "INSTRUMENT_INVALID"
    return "SUPPORTED" if all(gates.values()) else "NEGATIVE_RESULT"


def write_outputs(output_dir: Path, result: dict[str, Any], rows: list[dict[str, Any]], controls: list[dict[str, Any]], uncertainty: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "operational_closure_result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "operational_closure_uncertainty.json").write_text(json.dumps(uncertainty, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    prediction_fields = list(rows[0])
    with (output_dir / "operational_closure_predictions.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=prediction_fields)
        writer.writeheader()
        writer.writerows(rows)
    control_fields = list(controls[0])
    with (output_dir / "operational_closure_controls.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=control_fields)
        writer.writeheader()
        writer.writerows(controls)
    report = [
        "# EL-R4-OC-01 Operational Closure Result",
        "",
        f"- Classification: `{result['classification']}`",
        f"- Validation internal advantage: `{result['validation']['internal_advantage']:.6f}`",
        f"- Holdout internal advantage: `{result['holdout']['internal_advantage']:.6f}`",
        f"- Holdout 95% run-bootstrap CI: `{uncertainty['holdout_internal_advantage_ci_95']}`",
        f"- Positive holdout strata fraction: `{result['holdout_positive_strata_fraction']:.6f}`",
        f"- Result hash: `{result['result_hash']}`",
        "",
        "## Gate Results",
        "",
    ]
    report.extend(f"- `{name}`: `{passed}`" for name, passed in result["gates"].items())
    report.extend([
        "",
        "This experiment tests next-event predictive closure in stored Phase XVI event chains. It does not establish thermodynamic, biological, causal, cross-scale, or physical closure.",
        "",
    ])
    (output_dir / "operational_closure_report.md").write_text("\n".join(report), encoding="utf-8")


def run(output_dir: Path = RESULTS) -> dict[str, Any]:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN":
        raise ValueError("operational-closure contract is not frozen")
    source_path = REPO / contract["source"]["artifact"]
    runs = load_runs(source_path)
    by_level: dict[int, list[EventRun]] = defaultdict(list)
    for item in runs:
        by_level[item.n_side].append(item)
    construction_levels = set(contract["source"]["construction_levels"])
    construction = [item for item in runs if item.n_side in construction_levels]
    validation = by_level[contract["source"]["validation_level"]]
    holdout = by_level[contract["source"]["holdout_level"]]
    if len(holdout) < 30:
        raise ValueError("fewer than 30 holdout runs")

    counts = train_counts(construction)
    validation_rows, validation_metrics = score_runs(validation, counts)
    holdout_rows, holdout_metrics = score_runs(holdout, counts)
    ci = bootstrap_ci(run_level_advantages(holdout_rows), contract["seed_policy"]["bootstrap_seed"], contract["seed_policy"]["resamples"])
    positive_fraction = stratum_positive_fraction(holdout_rows)

    shuffled_counts = train_counts(shuffled_training_runs(construction, contract["seed_policy"]["control_seed"]))
    _, shuffled_metrics = score_runs(holdout, shuffled_counts)
    _, ablated_metrics = score_runs(holdout, counts, ablate_context=True)
    _, permuted_metrics = score_runs(holdout, counts, permute_targets=True, seed=contract["seed_policy"]["control_seed"])
    controls = [
        {
            "control": "transition_shuffled_training",
            "mechanism_valid": True,
            "internal_accuracy": shuffled_metrics["internal_accuracy"],
            "internal_advantage": shuffled_metrics["internal_advantage"],
            "degradation_from_target": holdout_metrics["internal_advantage"] - shuffled_metrics["internal_advantage"],
        },
        {
            "control": "context_ablation",
            "mechanism_valid": True,
            "internal_accuracy": ablated_metrics["internal_accuracy"],
            "internal_advantage": ablated_metrics["internal_advantage"],
            "degradation_from_target": holdout_metrics["internal_advantage"] - ablated_metrics["internal_advantage"],
        },
        {
            "control": "target_label_permutation",
            "mechanism_valid": True,
            "internal_accuracy": permuted_metrics["internal_accuracy"],
            "internal_advantage": permuted_metrics["internal_advantage"],
            "degradation_from_target": holdout_metrics["internal_advantage"] - permuted_metrics["internal_advantage"],
        },
    ]
    threshold = contract["success_threshold"]
    gates = {
        "validation_advantage": validation_metrics["internal_advantage"] >= threshold["validation_accuracy_advantage_min"],
        "holdout_advantage": holdout_metrics["internal_advantage"] >= threshold["holdout_accuracy_advantage_min"],
        "holdout_ci_excludes_zero": ci[0] > threshold["holdout_advantage_ci_lower_min"],
        "strata_consistency": positive_fraction >= threshold["holdout_positive_strata_fraction_min"],
        "transition_shuffle_degrades": controls[0]["degradation_from_target"] >= threshold["control_degradation_min"],
        "context_ablation_degrades": controls[1]["degradation_from_target"] >= threshold["control_degradation_min"],
    }
    classification = classify(gates)
    uncertainty = {
        "holdout_run_count": len(holdout),
        "holdout_internal_advantage_ci_95": list(ci),
        "resamples": contract["seed_policy"]["resamples"],
        "unit": "stored run key",
    }
    result = {
        "experiment_id": contract["experiment_id"],
        "classification": classification,
        "status": "PASS" if classification == "SUPPORTED" else "FAIL",
        "source_hash": source_hash(source_path),
        "contract_hash": stable_hash("el_r4_contract", contract),
        "run_counts": {str(level): len(items) for level, items in sorted(by_level.items())},
        "validation": validation_metrics,
        "holdout": holdout_metrics,
        "holdout_positive_strata_fraction": positive_fraction,
        "controls": controls,
        "gates": gates,
        "labels": [
            "RUNG_4_OPERATIONAL_CLOSURE_SUPPORTED" if classification == "SUPPORTED" else "RUNG_4_OPERATIONAL_CLOSURE_NOT_SUPPORTED",
            "INTERNAL_VS_EXTERNAL_PREDICTION_TESTED",
            "DESTRUCTIVE_CONTROLS_EXECUTED",
            "NO_CROSS_SCALE_OR_PHYSICS_CLAIM",
        ],
    }
    result["result_hash"] = stable_hash("el_r4_oc_01", result)
    write_outputs(output_dir, result, holdout_rows, controls, uncertainty)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Run EL-R4-OC-01 operational closure test.")
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    args = parser.parse_args()
    print(json.dumps(run(args.output_dir), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
