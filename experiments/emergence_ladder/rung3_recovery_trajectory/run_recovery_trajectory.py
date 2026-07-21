#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from haos_core import build_graph, run_transport


ROOT = Path(__file__).resolve().parent
CONTRACT = ROOT / "precommitment_contract.json"
RESULTS = ROOT / "results"


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def reference_state(graph: Any) -> np.ndarray:
    center = np.array([0.31, 0.43])
    delta = (graph.points - center + 0.5) % 1.0 - 0.5
    node = np.exp(-np.sum(delta**2, axis=1) / (2.0 * 0.12**2)).astype(complex)
    edge = np.asarray(graph.d0 @ node).ravel()
    edge *= 0.2 / max(float(np.linalg.norm(edge)), 1.0e-12)
    state = np.concatenate([node, edge, np.zeros(graph.block_sizes[2], dtype=complex)])
    return state / np.linalg.norm(state)


def orthogonal_perturbation(reference: np.ndarray, amplitude: float, seed: int) -> tuple[np.ndarray, float]:
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=len(reference)) + 1j * rng.normal(size=len(reference))
    noise -= reference * np.vdot(reference, noise)
    noise /= np.linalg.norm(noise)
    orthogonality = float(abs(np.vdot(reference, noise)))
    perturbed = reference + amplitude * noise
    return perturbed / np.linalg.norm(perturbed), orthogonality


def fidelity(left: np.ndarray, right: np.ndarray) -> float:
    denom = max(float(np.linalg.norm(left) * np.linalg.norm(right)), 1.0e-12)
    return float(abs(np.vdot(left, right)) ** 2 / denom**2)


def grade_fractions(state: np.ndarray, block_sizes: tuple[int, int, int]) -> np.ndarray:
    offsets = np.cumsum((0,) + block_sizes)
    energy = np.array([np.sum(np.abs(state[offsets[i] : offsets[i + 1]]) ** 2) for i in range(3)], dtype=float)
    return energy / max(float(np.sum(energy)), 1.0e-12)


def evolve(graph: Any, state: np.ndarray, config: dict[str, Any], *, beta: float, lock_alpha: float) -> list[np.ndarray]:
    states = run_transport(
        graph,
        {
            "mode": "grade_locked",
            "state0": state,
            "dt": config["dt"],
            "steps": config["steps"],
            "beta": beta,
            "lock_alpha": lock_alpha,
        },
    )
    if not all(np.isfinite(np.asarray(item)).all() for item in states):
        raise ValueError("non-finite trajectory value")
    return states


def trajectory_metrics(reference: list[np.ndarray], perturbed: list[np.ndarray], block_sizes: tuple[int, int, int]) -> dict[str, Any]:
    series = [fidelity(left, right) for left, right in zip(reference, perturbed)]
    initial = series[0]
    final = series[-1]
    target = initial + 0.9 * (1.0 - initial)
    recovery_time = next((index for index, value in enumerate(series) if value >= target), None)
    grade_error = float(np.linalg.norm(grade_fractions(reference[-1], block_sizes) - grade_fractions(perturbed[-1], block_sizes), ord=1))
    return {
        "initial_fidelity": initial,
        "final_fidelity": final,
        "recovery_gain": final - initial,
        "recovery_time_step": recovery_time,
        "grade_signature_error": grade_error,
        "series": series,
    }


def trivial_trajectory(initial: np.ndarray, attractor: np.ndarray, steps: int, rho: float = 0.12) -> list[np.ndarray]:
    states = [initial.copy()]
    state = initial.copy()
    for _ in range(steps):
        state = (1.0 - rho) * state + rho * attractor
        state /= np.linalg.norm(state)
        states.append(state.copy())
    return states


def bootstrap_ci(values: list[float], seed: int, resamples: int) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    samples = [float(np.mean(rng.choice(array, size=len(array), replace=True))) for _ in range(resamples)]
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def classify(gates: dict[str, bool], instrument_valid: bool = True) -> str:
    if not instrument_valid:
        return "INSTRUMENT_INVALID"
    return "SUPPORTED" if all(gates.values()) else "NEGATIVE_RESULT"


def run(output_dir: Path = RESULTS) -> dict[str, Any]:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN":
        raise ValueError("recovery-trajectory contract is not frozen")
    config = contract["config"]
    branch_graph = build_graph({"kind": "dk2d_periodic", "n_side": config["n_side"], "epsilon": config["epsilon"]})
    topology_graph = build_graph({
        "kind": "dk2d_periodic",
        "n_side": config["n_side"],
        "epsilon": config["epsilon"],
        "cycle_phase_x": float(np.pi),
    })
    if np.allclose(branch_graph.dirac_kahler.toarray(), topology_graph.dirac_kahler.toarray()):
        raise ValueError("topology control did not alter the frozen operator")
    reference = reference_state(branch_graph)
    topology_reference = reference_state(topology_graph)
    branch_reference_trajectory = evolve(branch_graph, reference, config, beta=config["beta"], lock_alpha=config["lock_alpha"])
    operator_reference_trajectory = evolve(branch_graph, reference, config, beta=0.0, lock_alpha=0.0)
    topology_reference_trajectory = evolve(topology_graph, topology_reference, config, beta=config["beta"], lock_alpha=config["lock_alpha"])

    zero_mode = np.zeros_like(reference)
    zero_mode[: branch_graph.block_sizes[0]] = 1.0
    zero_mode /= np.linalg.norm(zero_mode)
    rows = []
    for amplitude in contract["perturbation_amplitudes"]:
        for seed in contract["seed_policy"]["perturbation_seeds"]:
            perturbed, orthogonality = orthogonal_perturbation(reference, amplitude, seed)
            if orthogonality > 1.0e-10:
                raise ValueError("perturbation is not orthogonal to reference")
            topology_perturbed, _ = orthogonal_perturbation(topology_reference, amplitude, seed)
            branch = trajectory_metrics(
                branch_reference_trajectory,
                evolve(branch_graph, perturbed, config, beta=config["beta"], lock_alpha=config["lock_alpha"]),
                branch_graph.block_sizes,
            )
            operator_only = trajectory_metrics(
                operator_reference_trajectory,
                evolve(branch_graph, perturbed, config, beta=0.0, lock_alpha=0.0),
                branch_graph.block_sizes,
            )
            topology = trajectory_metrics(
                topology_reference_trajectory,
                evolve(topology_graph, topology_perturbed, config, beta=config["beta"], lock_alpha=config["lock_alpha"]),
                topology_graph.block_sizes,
            )
            passive = trajectory_metrics(
                [reference] * (config["steps"] + 1),
                [perturbed] * (config["steps"] + 1),
                branch_graph.block_sizes,
            )
            trivial_reference = trivial_trajectory(reference, zero_mode, config["steps"])
            trivial_perturbed = trivial_trajectory(perturbed, zero_mode, config["steps"])
            trivial = trajectory_metrics(trivial_reference, trivial_perturbed, branch_graph.block_sizes)
            trivial_identity = fidelity(trivial_perturbed[-1], branch_reference_trajectory[-1])
            rows.append({
                "amplitude": amplitude,
                "seed": seed,
                "initial_fidelity": branch["initial_fidelity"],
                "branch_final_fidelity": branch["final_fidelity"],
                "branch_recovery_gain": branch["recovery_gain"],
                "branch_recovery_time_step": branch["recovery_time_step"],
                "branch_grade_signature_error": branch["grade_signature_error"],
                "passive_recovery_gain": passive["recovery_gain"],
                "operator_only_recovery_gain": operator_only["recovery_gain"],
                "topology_recovery_gain": topology["recovery_gain"],
                "trivial_recovery_gain": trivial["recovery_gain"],
                "trivial_identity_to_branch": trivial_identity,
            })

    means = {key: float(np.mean([row[key] for row in rows])) for key in (
        "branch_recovery_gain",
        "passive_recovery_gain",
        "operator_only_recovery_gain",
        "topology_recovery_gain",
        "trivial_recovery_gain",
        "branch_grade_signature_error",
        "trivial_identity_to_branch",
    )}
    ci = bootstrap_ci([row["branch_recovery_gain"] for row in rows], contract["seed_policy"]["bootstrap_seed"], contract["seed_policy"]["resamples"])
    threshold = contract["success_threshold"]
    recovered_fraction = float(np.mean([row["branch_recovery_gain"] >= threshold["branch_gain_min"] for row in rows]))
    gates = {
        "branch_gain": means["branch_recovery_gain"] >= threshold["branch_gain_min"],
        "branch_gain_ci": ci[0] >= threshold["branch_gain_ci_lower_min"],
        "beats_passive": means["branch_recovery_gain"] - means["passive_recovery_gain"] >= threshold["branch_over_passive_min"],
        "beats_operator_only": means["branch_recovery_gain"] - means["operator_only_recovery_gain"] >= threshold["branch_over_operator_only_min"],
        "topology_specific": means["branch_recovery_gain"] - means["topology_recovery_gain"] >= threshold["topology_degradation_min"],
        "grade_signature_retained": means["branch_grade_signature_error"] <= threshold["grade_signature_error_max"],
        "run_fraction": recovered_fraction >= threshold["recovered_run_fraction_min"],
        "trivial_attractor_rejected": means["trivial_identity_to_branch"] <= threshold["trivial_identity_max"],
    }
    classification = classify(gates)
    controls = [
        {"control": "passive_relaxation", "mean_recovery_gain": means["passive_recovery_gain"], "mechanism_valid": True},
        {"control": "operator_only_filtering", "mean_recovery_gain": means["operator_only_recovery_gain"], "mechanism_valid": True},
        {"control": "topology_altered_cycle_phase", "mean_recovery_gain": means["topology_recovery_gain"], "mechanism_valid": True},
        {"control": "trivial_attractor", "mean_recovery_gain": means["trivial_recovery_gain"], "mechanism_valid": True},
    ]
    uncertainty = {
        "branch_recovery_gain_ci_95": list(ci),
        "independent_parameter_seed_rows": len(rows),
        "resamples": contract["seed_policy"]["resamples"],
    }
    result = {
        "experiment_id": contract["experiment_id"],
        "classification": classification,
        "status": "PASS" if classification == "SUPPORTED" else "FAIL",
        "contract_hash": stable_hash("el_r3_contract", contract),
        "means": means,
        "recovered_run_fraction": recovered_fraction,
        "gates": gates,
        "controls": controls,
        "uncertainty": uncertainty,
        "labels": [
            "RUNG_3_RECOVERABILITY_SUPPORTED" if classification == "SUPPORTED" else "RUNG_3_RECOVERABILITY_NOT_SUPPORTED",
            "PASSIVE_AND_TRIVIAL_CONTROLS_EXECUTED",
            "OPERATOR_ONLY_EXPLANATION_TESTED",
            "NO_CLOSURE_OR_PHYSICS_CLAIM",
        ],
    }
    result["result_hash"] = stable_hash("el_r3_rt_01", result)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "recovery_trajectory_result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "recovery_trajectory_uncertainty.json").write_text(json.dumps(uncertainty, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output_dir / "recovery_trajectory_rows.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (output_dir / "recovery_trajectory_controls.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(controls[0]))
        writer.writeheader()
        writer.writerows(controls)
    report = [
        "# EL-R3-RT-01 Recovery Trajectory Result",
        "",
        f"- Classification: `{classification}`",
        f"- Mean branch recovery gain: `{means['branch_recovery_gain']:.6f}`",
        f"- Mean operator-only gain: `{means['operator_only_recovery_gain']:.6f}`",
        f"- Branch gain 95% CI: `{list(ci)}`",
        f"- Result hash: `{result['result_hash']}`",
        "",
        "This test distinguishes post-perturbation trajectory recovery from persistence, passive relaxation, ordinary operator filtering, topology alteration, and trivial-attractor collapse. It makes no closure, cross-scale, continuum, or physical claim.",
        "",
    ]
    (output_dir / "recovery_trajectory_report.md").write_text("\n".join(report), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Run EL-R3-RT-01 recovery trajectory test.")
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    args = parser.parse_args()
    print(json.dumps(run(args.output_dir), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
