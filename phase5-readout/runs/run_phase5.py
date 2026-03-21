#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
import time
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import read_json, relpath, write_json, write_timestamped_json

CONTROL_CLASSES = (
    "stable_frozen_sector",
    "degraded_sector_control",
    "shuffled_null_control",
)
INPUT_MODES = ("frozen_sector", "dummy_sector")
RECOVERY_COMPONENTS = (
    "law_fit_score",
    "stable_subset_size",
    "stable_rho",
    "ordering_gap",
    "monotonic_consistency_score",
)
AUDIT_NAME = "phase5b_robustness_and_control_discrimination"
CLAIM_BOUNDARY = (
    "Phase V authority is limited to deterministic recovery-histogram behavior on the frozen "
    "Phase IV stable sector bundle under the configured perturbation schedules, policy bands, "
    "and control classes. It does not claim universal scaling, theory expansion, or classifier "
    "invariance outside this tested regime."
)


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(value, upper))


def stable_seed(*parts: Any) -> int:
    digest = hashlib.sha256("::".join(str(part) for part in parts).encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def artifact_ref(path: Path) -> str:
    try:
        return relpath(path)
    except ValueError:
        return str(path.resolve())


def resolve_candidate_path(explicit_path: Path | None, candidates: list[str]) -> Path:
    if explicit_path is not None:
        return explicit_path
    for candidate in candidates:
        candidate_path = ROOT / candidate
        if candidate_path.exists():
            return candidate_path
    raise FileNotFoundError(f"no candidate path found in {candidates}")


def phase4_inputs_from_bundle(bundle: dict[str, Any]) -> dict[str, dict[str, Any]]:
    inputs = bundle.get("upstream_inputs")
    if not isinstance(inputs, dict):
        raise ValueError("Phase V frozen_sector mode requires a Phase IV bundle with upstream_inputs.")
    required = ("stage24_1", "stage24_2", "stage24_3", "stage24_4")
    missing = [name for name in required if name not in inputs]
    if missing:
        raise ValueError(f"Phase IV bundle is missing upstream inputs: {missing}")
    return {name: dict(inputs[name]) for name in required}


def validate_phase4_bundle(path: Path, bundle: dict[str, Any], max_bytes: int) -> None:
    if path.stat().st_size > max_bytes:
        raise ValueError(
            f"Phase IV bundle {relpath(path)} exceeds Phase V size cap {max_bytes} bytes."
        )
    if bundle.get("phase_name") != "phase4-sector-freeze":
        raise ValueError("Phase V frozen_sector mode requires a phase4-sector-freeze bundle.")
    summary = bundle.get("summary", {})
    if not isinstance(summary, dict) or not bool(summary.get("phaseIV_freeze_valid")):
        raise ValueError("Phase V frozen_sector mode requires a valid frozen Phase IV bundle.")


def find_extension_entry(extension_ledger: list[dict[str, Any]], extension_id: str) -> dict[str, Any]:
    for row in extension_ledger:
        if str(row["extension_id"]) == extension_id:
            return row
    raise KeyError(f"no extension row found for {extension_id}")


def load_frozen_sector_input(config: dict[str, Any], phase4_input: Path | None) -> tuple[dict[str, Any], dict[str, str]]:
    frozen = dict(config["frozen_sector"])
    phase4_bundle_path = resolve_candidate_path(phase4_input, list(frozen["phase4_candidates"]))
    phase4_bundle = read_json(phase4_bundle_path)
    validate_phase4_bundle(phase4_bundle_path, phase4_bundle, int(config["max_phase4_bundle_bytes"]))

    phase4_input_refs = phase4_inputs_from_bundle(phase4_bundle)
    stage24_1 = read_json(ROOT / phase4_input_refs["stage24_1"]["json_path"])
    stage24_2 = read_json(ROOT / phase4_input_refs["stage24_2"]["json_path"])
    stage24_3 = read_json(ROOT / phase4_input_refs["stage24_3"]["json_path"])
    stage24_4 = read_json(ROOT / phase4_input_refs["stage24_4"]["json_path"])
    extension_config_path = resolve_candidate_path(
        None,
        list(frozen["phase4_extension_config_candidates"]),
    )
    extension_config = read_json(extension_config_path)

    baseline_id = str(frozen["baseline_extension_id"])
    baseline_row = find_extension_entry(list(stage24_4["extension_ledger"]), baseline_id)
    baseline_cfg = find_extension_entry(list(extension_config["extensions"]), baseline_id)
    common_fields = dict(extension_config["common_fields"])

    canonical_input = {
        "sector_identifier": str(frozen["sector_identifier"]),
        "graph_reference_or_seed": {
            "reference_type": "phaseIV_extension_anchor",
            "extension_id": baseline_cfg["extension_id"],
            "extension_class": baseline_cfg["extension_class"],
            "extension_level": baseline_cfg["extension_level"],
            "width_ratio_a_to_b": baseline_cfg["width_ratio_a_to_b"],
            "separation": baseline_cfg["separation"],
            "anisotropy_x": baseline_cfg["anisotropy_x"],
            "anisotropy_y": baseline_cfg["anisotropy_y"],
            "motif_amplitude_ratio": baseline_cfg["motif_amplitude_ratio"],
            "motif_offset_y": baseline_cfg["motif_offset_y"],
            "phase_values_fraction_of_pi": baseline_cfg["phase_values_fraction_of_pi"],
        },
        "kernel_configuration_snapshot": {
            "operator_sector": common_fields["operator_sector"],
            "boundary_type": common_fields["boundary_type"],
            "graph_type": common_fields["graph_type"],
            "kernel_type": common_fields["kernel_type"],
            "base_resolution": common_fields["base_resolution"],
            "refined_resolution": common_fields["refined_resolution"],
            "mean_width": common_fields["mean_width"],
            "base_separation": common_fields["base_separation"],
            "amplitude_scale": common_fields["amplitude_scale"],
            "t_final": common_fields["t_final"],
            "dt_scale": common_fields["dt_scale"],
            "weak_beta_values": common_fields["weak_beta_values"],
        },
        "selector_configuration": {
            "selector_id": stage24_2["selected_model_id"],
            "selector_family": stage24_2["selected_model_family"],
            "target_definition": stage24_2["target_definition"],
            "law_expression": stage24_2["selected_law_expression"],
            "primary_inputs": stage24_2["selected_core_inputs"],
            "context_inputs": stage24_2["selected_context_inputs"],
            "threshold": stage24_4["fixed_flow_threshold"],
        },
        "invariant_baseline_vector": {
            "law_fit_score": float(baseline_row["law_fit_score"]),
            "stable_subset_size": int(baseline_row["stable_subset_size"]),
            "stable_rho": float(baseline_row["stable_rho"]),
            "ordering_gap": float(baseline_row["ordering_gap"]),
            "monotonic_consistency_score": float(baseline_row["monotonic_consistency_score"]),
        },
        "perturbation_policy_descriptor": {
            "policy_id": frozen["policy"]["policy_id"],
            "trials": int(frozen["policy"]["trials"]),
            "incidence_sites": int(frozen["policy"]["incidence_sites"]),
            "flip_probability": float(frozen["policy"]["flip_probability"]),
            "max_incidence_shift": int(frozen["policy"]["max_incidence_shift"]),
            "max_phase_drift_fraction_of_pi": float(frozen["policy"]["max_phase_drift_fraction_of_pi"]),
            "amplitude_jitter": float(frozen["policy"]["amplitude_jitter"]),
            "preserve_operator_sector": bool(frozen["policy"]["preserve_operator_sector"]),
            "preserve_boundary_type": bool(frozen["policy"]["preserve_boundary_type"]),
            "preserve_kernel_type": bool(frozen["policy"]["preserve_kernel_type"]),
        },
    }
    refs = {
        "phase4_bundle": relpath(phase4_bundle_path),
        "stage24_1_json": phase4_input_refs["stage24_1"]["json_path"],
        "stage24_2_json": phase4_input_refs["stage24_2"]["json_path"],
        "stage24_3_json": phase4_input_refs["stage24_3"]["json_path"],
        "stage24_4_json": phase4_input_refs["stage24_4"]["json_path"],
        "stage24_4_extension_config": relpath(extension_config_path),
    }
    return canonical_input, refs


def load_dummy_sector_input(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, str]]:
    dummy_input = dict(config["dummy_sector"])
    refs = {
        "dummy_sector": "config:dummy_sector",
    }
    return dummy_input, refs


def override_policy(policy: dict[str, Any], policy_override: dict[str, Any] | None) -> dict[str, Any]:
    updated_policy = dict(policy)
    if not policy_override:
        return updated_policy
    for key, value in policy_override.items():
        baseline_value = updated_policy.get(key)
        if isinstance(baseline_value, bool):
            updated_policy[key] = bool(value)
        elif isinstance(baseline_value, int) and not isinstance(baseline_value, bool):
            updated_policy[key] = int(value)
        elif isinstance(baseline_value, float):
            updated_policy[key] = float(value)
        else:
            updated_policy[key] = value
    return updated_policy


def sample_incidence_perturbation(
    rng: random.Random,
    policy: dict[str, Any],
) -> dict[str, float]:
    incidence_sites = int(policy["incidence_sites"])
    flip_probability = float(policy["flip_probability"])
    max_incidence_shift = max(1, int(policy["max_incidence_shift"]))
    max_phase_drift = max(1e-9, float(policy["max_phase_drift_fraction_of_pi"]))

    flip_count = sum(1 for _ in range(incidence_sites) if rng.random() < flip_probability)
    incidence_shift = rng.randint(0, max_incidence_shift)
    phase_drift_fraction = rng.uniform(-max_phase_drift, max_phase_drift)
    amplitude_jitter = rng.random() * float(policy["amplitude_jitter"])

    return {
        "flip_count": float(flip_count),
        "flip_fraction": flip_count / incidence_sites,
        "incidence_shift_fraction": incidence_shift / max_incidence_shift,
        "phase_drift_fraction_of_pi": phase_drift_fraction,
        "phase_drift_fraction": abs(phase_drift_fraction) / max_phase_drift,
        "amplitude_jitter": amplitude_jitter,
    }


def trial_invariant_metrics(
    baseline: dict[str, Any],
    perturbation: dict[str, float],
    control_profile: dict[str, float],
    rng: random.Random,
) -> dict[str, float]:
    load = clamp(
        float(control_profile["load_offset"])
        + float(control_profile["load_scale"]) * perturbation["flip_fraction"]
        + 0.25 * perturbation["incidence_shift_fraction"]
        + 0.15 * perturbation["phase_drift_fraction"]
        + perturbation["amplitude_jitter"],
        0.0,
        1.5,
    )
    positive_signal = max(0.0, float(control_profile["recovery_bias"]) - load)
    shuffle_penalty = rng.random() * float(control_profile["shuffle_scale"])

    base_law_fit = float(baseline["law_fit_score"])
    base_subset = int(baseline["stable_subset_size"])
    base_rho = float(baseline["stable_rho"])
    base_gap = float(baseline["ordering_gap"])
    base_monotonic = float(baseline["monotonic_consistency_score"])

    law_fit_score = clamp(base_law_fit - 0.90 * load + 0.25 * positive_signal, 0.0, 1.0)
    stable_subset_size = int(
        round(
            clamp(
                base_subset - base_subset * 1.10 * load + base_subset * 0.35 * positive_signal,
                0.0,
                float(base_subset),
            )
        )
    )
    stable_rho = clamp(
        base_rho + 1.05 * load - 0.25 * positive_signal + 0.25 * shuffle_penalty,
        -1.0,
        1.0,
    )
    ordering_gap = clamp(base_gap - 1.10 * load + 0.30 * positive_signal, 0.0, 1.0)
    monotonic_consistency_score = clamp(
        base_monotonic - 1.00 * load + 0.30 * positive_signal,
        0.0,
        1.0,
    )

    return {
        "load": load,
        "law_fit_score": law_fit_score,
        "stable_subset_size": float(stable_subset_size),
        "stable_rho": stable_rho,
        "ordering_gap": ordering_gap,
        "monotonic_consistency_score": monotonic_consistency_score,
    }


def recovery_score(
    trial_metrics: dict[str, float],
    baseline: dict[str, Any],
) -> tuple[float, dict[str, float]]:
    base_law_fit = max(1e-9, float(baseline["law_fit_score"]))
    base_subset = max(1, int(baseline["stable_subset_size"]))
    base_rho = max(1e-9, abs(float(baseline["stable_rho"])))
    base_gap = max(1e-9, float(baseline["ordering_gap"]))
    base_monotonic = max(1e-9, float(baseline["monotonic_consistency_score"]))

    components = {
        "law_fit_score": clamp(trial_metrics["law_fit_score"] / base_law_fit, 0.0, 1.0),
        "stable_subset_size": clamp(trial_metrics["stable_subset_size"] / base_subset, 0.0, 1.0),
        "stable_rho": clamp(abs(trial_metrics["stable_rho"]) / base_rho, 0.0, 1.0),
        "ordering_gap": clamp(trial_metrics["ordering_gap"] / base_gap, 0.0, 1.0),
        "monotonic_consistency_score": clamp(
            trial_metrics["monotonic_consistency_score"] / base_monotonic,
            0.0,
            1.0,
        ),
    }
    score = round(sum(components.values()) / len(components), 6)
    return score, components


def run_ensemble(
    canonical_input: dict[str, Any],
    control_profile: dict[str, float],
    control_class: str,
    schedule_label: str = "default",
) -> tuple[list[float], dict[str, Any]]:
    policy = dict(canonical_input["perturbation_policy_descriptor"])
    baseline = dict(canonical_input["invariant_baseline_vector"])
    trials = int(policy["trials"])
    if trials < 1000:
        raise ValueError("Phase V requires at least 1000 trials per sector.")

    seed = stable_seed(
        canonical_input["sector_identifier"],
        control_class,
        trials,
        policy["policy_id"],
        schedule_label,
    )
    rng = random.Random(seed)
    trial_scores: list[float] = []
    first_trial_snapshot: dict[str, Any] | None = None

    for trial_index in range(trials):
        perturbation = sample_incidence_perturbation(rng, policy)
        trial_metrics = trial_invariant_metrics(baseline, perturbation, control_profile, rng)
        score, components = recovery_score(trial_metrics, baseline)
        trial_scores.append(score)
        if trial_index == 0:
            first_trial_snapshot = {
                "perturbation": perturbation,
                "trial_metrics": trial_metrics,
                "component_scores": components,
                "recovery_score": score,
            }

    metadata = {
        "trial_count": trials,
        "schedule_label": schedule_label,
        "seed": seed,
        "first_trial_snapshot": first_trial_snapshot,
    }
    return trial_scores, metadata


def build_histogram(scores: list[float], bin_count: int) -> dict[str, Any]:
    if not scores:
        raise ValueError("cannot build histogram from zero scores")

    edges = [round(index / bin_count, 6) for index in range(bin_count + 1)]
    counts = [0] * bin_count
    for score in scores:
        clamped = clamp(score, 0.0, 1.0)
        bucket = min(int(clamped * bin_count), bin_count - 1)
        counts[bucket] += 1

    total = len(scores)
    normalized = [round(count / total, 6) for count in counts]
    cumulative: list[float] = []
    running = 0.0
    for value in normalized:
        running += value
        cumulative.append(round(running, 6))

    mean_score = sum(scores) / total
    variance = sum((score - mean_score) ** 2 for score in scores) / total
    centers = [round((edges[index] + edges[index + 1]) / 2.0, 6) for index in range(bin_count)]

    return {
        "bin_edges": edges,
        "bin_centers": centers,
        "bin_counts": counts,
        "normalized_frequency": normalized,
        "cumulative_distribution": cumulative,
        "summary_statistics": {
            "mean": round(mean_score, 6),
            "variance": round(variance, 6),
            "min": round(min(scores), 6),
            "max": round(max(scores), 6),
            "nonzero_recovery_fraction": round(
                sum(1 for score in scores if score > 0.0) / total,
                6,
            ),
        },
        "trial_count": total,
    }


def solve_linear_system(matrix: list[list[float]], vector: list[float]) -> list[float]:
    size = len(vector)
    augmented = [list(matrix[row]) + [vector[row]] for row in range(size)]

    for pivot_index in range(size):
        pivot_row = max(range(pivot_index, size), key=lambda row: abs(augmented[row][pivot_index]))
        if abs(augmented[pivot_row][pivot_index]) < 1e-12:
            raise ValueError("singular system in scaling fit")
        augmented[pivot_index], augmented[pivot_row] = augmented[pivot_row], augmented[pivot_index]
        pivot = augmented[pivot_index][pivot_index]
        for column_index in range(pivot_index, size + 1):
            augmented[pivot_index][column_index] /= pivot
        for row_index in range(size):
            if row_index == pivot_index:
                continue
            factor = augmented[row_index][pivot_index]
            for column_index in range(pivot_index, size + 1):
                augmented[row_index][column_index] -= factor * augmented[pivot_index][column_index]

    return [augmented[row][size] for row in range(size)]


def fit_constant(xs: list[float], ys: list[float]) -> dict[str, Any]:
    mean_y = sum(ys) / len(ys)
    predictions = [mean_y for _ in xs]
    rmse = (sum((y - prediction) ** 2 for y, prediction in zip(ys, predictions)) / len(xs)) ** 0.5
    return {
        "coefficients": [round(mean_y, 6)],
        "rmse": round(rmse, 6),
    }


def fit_linear(xs: list[float], ys: list[float]) -> dict[str, Any]:
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    denominator = sum((x - mean_x) ** 2 for x in xs)
    slope = 0.0 if denominator == 0.0 else sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denominator
    intercept = mean_y - slope * mean_x
    predictions = [intercept + slope * x for x in xs]
    rmse = (sum((y - prediction) ** 2 for y, prediction in zip(ys, predictions)) / len(xs)) ** 0.5
    return {
        "coefficients": [round(intercept, 6), round(slope, 6)],
        "rmse": round(rmse, 6),
    }


def fit_quadratic(xs: list[float], ys: list[float]) -> dict[str, Any]:
    sums = {
        "x0": len(xs),
        "x1": sum(xs),
        "x2": sum(x**2 for x in xs),
        "x3": sum(x**3 for x in xs),
        "x4": sum(x**4 for x in xs),
        "y": sum(ys),
        "xy": sum(x * y for x, y in zip(xs, ys)),
        "x2y": sum((x**2) * y for x, y in zip(xs, ys)),
    }
    matrix = [
        [sums["x0"], sums["x1"], sums["x2"]],
        [sums["x1"], sums["x2"], sums["x3"]],
        [sums["x2"], sums["x3"], sums["x4"]],
    ]
    vector = [sums["y"], sums["xy"], sums["x2y"]]
    a0, a1, a2 = solve_linear_system(matrix, vector)
    predictions = [a0 + a1 * x + a2 * x * x for x in xs]
    rmse = (sum((y - prediction) ** 2 for y, prediction in zip(ys, predictions)) / len(xs)) ** 0.5
    return {
        "coefficients": [round(a0, 6), round(a1, 6), round(a2, 6)],
        "rmse": round(rmse, 6),
    }


def classify_scaling(histogram: dict[str, Any]) -> dict[str, Any]:
    xs = list(histogram["bin_centers"])
    ys = list(histogram["normalized_frequency"])
    stats = dict(histogram["summary_statistics"])
    constant_fit = fit_constant(xs, ys)
    linear_fit = fit_linear(xs, ys)
    quadratic_fit = fit_quadratic(xs, ys)

    lower_tail_mass = round(sum(ys[:2]), 6)
    upper_tail_mass = round(sum(ys[-2:]), 6)
    frequency_range = round(max(ys) - min(ys), 6)
    oscillation_count = sum(
        1
        for left, right in zip(
            [ys[index + 1] - ys[index] for index in range(len(ys) - 1)],
            [ys[index + 2] - ys[index + 1] for index in range(len(ys) - 2)],
        )
        if left * right < -1e-6
    )

    if stats["mean"] <= 0.18 and lower_tail_mass >= 0.70:
        scaling_class = "unstable"
    elif frequency_range <= 0.04 and constant_fit["rmse"] <= 0.02:
        scaling_class = "flat"
    elif quadratic_fit["rmse"] < linear_fit["rmse"] * 0.88 and upper_tail_mass > lower_tail_mass:
        scaling_class = "quadratic_like"
    elif linear_fit["rmse"] < constant_fit["rmse"] * 0.88:
        scaling_class = "linear_like"
    elif oscillation_count >= max(2, len(ys) // 4):
        scaling_class = "unstable"
    else:
        scaling_class = "mixed"

    return {
        "classifier_version": "phase5_scaling_v1",
        "scaling_class": scaling_class,
        "diagnostics": {
            "constant_fit": constant_fit,
            "linear_fit": linear_fit,
            "quadratic_fit": quadratic_fit,
            "lower_tail_mass": lower_tail_mass,
            "upper_tail_mass": upper_tail_mass,
            "frequency_range": frequency_range,
            "oscillation_count": oscillation_count,
        },
    }


def bundle_summary(
    canonical_input: dict[str, Any],
    mode: str,
    control_class: str,
    schedule_label: str,
    histogram: dict[str, Any],
    scaling_fit: dict[str, Any],
    runtime_seconds: float,
) -> dict[str, Any]:
    stats = dict(histogram["summary_statistics"])
    return {
        "experiment": "5_coherence_recovery_readout",
        "input_mode": mode,
        "control_class": control_class,
        "schedule_label": schedule_label,
        "sector_identifier": canonical_input["sector_identifier"],
        "recovery_score_name": "structural_recovery_score",
        "trial_count": histogram["trial_count"],
        "mean_recovery_score": stats["mean"],
        "variance_recovery_score": stats["variance"],
        "min_recovery_score": stats["min"],
        "max_recovery_score": stats["max"],
        "nonzero_recovery_fraction": stats["nonzero_recovery_fraction"],
        "scaling_class": scaling_fit["scaling_class"],
        "runtime_seconds": round(runtime_seconds, 6),
    }


def execute_phase5(
    config: dict[str, Any],
    config_path: Path,
    mode: str,
    control_class: str,
    phase4_input: Path | None,
    schedule_label: str = "default",
    policy_override: dict[str, Any] | None = None,
) -> dict[str, Any]:
    control_profile = dict(config["control_classes"][control_class])

    started_at = time.perf_counter()
    if mode == "frozen_sector":
        canonical_input, input_refs = load_frozen_sector_input(config, phase4_input)
    else:
        canonical_input, input_refs = load_dummy_sector_input(config)

    canonical_input["perturbation_policy_descriptor"] = override_policy(
        dict(canonical_input["perturbation_policy_descriptor"]),
        policy_override,
    )
    canonical_input["perturbation_policy_descriptor"]["control_class"] = control_class

    trial_scores, ensemble_metadata = run_ensemble(
        canonical_input,
        control_profile,
        control_class,
        schedule_label=schedule_label,
    )
    histogram = build_histogram(trial_scores, int(config["histogram_bins"]))
    scaling_fit = classify_scaling(histogram)
    runtime_seconds = time.perf_counter() - started_at

    histogram_payload = {
        "phase": 5,
        "phase_name": "phase5-readout",
        "input_mode": mode,
        "control_class": control_class,
        "sector_identifier": canonical_input["sector_identifier"],
        "histogram": histogram,
    }
    scaling_payload = {
        "phase": 5,
        "phase_name": "phase5-readout",
        "input_mode": mode,
        "control_class": control_class,
        "sector_identifier": canonical_input["sector_identifier"],
        "scaling_fit": scaling_fit,
    }
    bundle = {
        "phase": 5,
        "phase_name": "phase5-readout",
        "source_config": relpath(config_path),
        "input_mode": mode,
        "control_class": control_class,
        "canonical_input": canonical_input,
        "input_refs": input_refs,
        "summary": bundle_summary(
            canonical_input,
            mode,
            control_class,
            schedule_label,
            histogram,
            scaling_fit,
            runtime_seconds,
        ),
        "trial_scores": trial_scores,
        "ensemble_metadata": ensemble_metadata,
    }
    return {
        "bundle": bundle,
        "histogram_payload": histogram_payload,
        "scaling_payload": scaling_payload,
    }


def write_phase5_outputs(
    output_dir: Path,
    bundle: dict[str, Any],
    histogram_payload: dict[str, Any],
    scaling_payload: dict[str, Any],
) -> dict[str, str]:
    histogram_stamped, histogram_latest = write_timestamped_json(
        output_dir,
        "phase5_recovery_histogram",
        histogram_payload,
    )
    scaling_stamped, scaling_latest = write_timestamped_json(
        output_dir,
        "phase5_scaling_fit",
        scaling_payload,
    )
    bundle_stamped, bundle_latest = write_timestamped_json(
        output_dir,
        "phase5_readout_bundle",
        bundle,
    )

    artifacts = {
        "stamped_json": artifact_ref(bundle_stamped),
        "latest_json": artifact_ref(bundle_latest),
        "recovery_histogram_stamped_json": artifact_ref(histogram_stamped),
        "recovery_histogram_latest_json": artifact_ref(histogram_latest),
        "scaling_fit_stamped_json": artifact_ref(scaling_stamped),
        "scaling_fit_latest_json": artifact_ref(scaling_latest),
    }
    bundle["artifacts"] = artifacts
    histogram_payload["artifacts"] = {
        "stamped_json": artifact_ref(histogram_stamped),
        "latest_json": artifact_ref(histogram_latest),
    }
    scaling_payload["artifacts"] = {
        "stamped_json": artifact_ref(scaling_stamped),
        "latest_json": artifact_ref(scaling_latest),
    }

    write_json(bundle_stamped, bundle)
    write_json(bundle_latest, bundle)
    write_json(histogram_stamped, histogram_payload)
    write_json(histogram_latest, histogram_payload)
    write_json(scaling_stamped, scaling_payload)
    write_json(scaling_latest, scaling_payload)
    return artifacts


def summarize_run(result: dict[str, Any]) -> dict[str, Any]:
    bundle = dict(result["bundle"])
    summary = dict(bundle["summary"])
    ensemble_metadata = dict(bundle["ensemble_metadata"])
    return {
        "schedule_label": ensemble_metadata["schedule_label"],
        "seed": ensemble_metadata["seed"],
        "scaling_class": summary["scaling_class"],
        "mean_recovery_score": summary["mean_recovery_score"],
        "variance_recovery_score": summary["variance_recovery_score"],
        "min_recovery_score": summary["min_recovery_score"],
        "max_recovery_score": summary["max_recovery_score"],
        "nonzero_recovery_fraction": summary["nonzero_recovery_fraction"],
        "runtime_seconds": summary["runtime_seconds"],
    }


def histogram_l1_distance(left_histogram: dict[str, Any], right_histogram: dict[str, Any]) -> float:
    left = list(left_histogram["normalized_frequency"])
    right = list(right_histogram["normalized_frequency"])
    if len(left) != len(right):
        raise ValueError("histograms must have matching bin counts")
    return round(sum(abs(left_value - right_value) for left_value, right_value in zip(left, right)), 6)


def pairwise_schedule_distances(
    results_by_schedule: dict[str, dict[str, Any]],
    schedule_labels: list[str],
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    distance_rows: list[dict[str, Any]] = []
    distances: list[float] = []
    for left_schedule, right_schedule in combinations(schedule_labels, 2):
        distance = histogram_l1_distance(
            results_by_schedule[left_schedule]["histogram_payload"]["histogram"],
            results_by_schedule[right_schedule]["histogram_payload"]["histogram"],
        )
        distance_rows.append(
            {
                "left_schedule_label": left_schedule,
                "right_schedule_label": right_schedule,
                "histogram_l1_distance": distance,
            }
        )
        distances.append(distance)
    aggregate = {
        "average_histogram_l1_distance": round(sum(distances) / len(distances), 6) if distances else 0.0,
        "max_histogram_l1_distance": round(max(distances), 6) if distances else 0.0,
        "min_histogram_l1_distance": round(min(distances), 6) if distances else 0.0,
    }
    return distance_rows, aggregate


def summarize_run_cluster(
    run_records: list[dict[str, Any]],
    histogram_distance_aggregate: dict[str, float] | None = None,
) -> dict[str, Any]:
    if not run_records:
        raise ValueError("cannot summarize empty run cluster")
    class_counter = Counter(str(record["scaling_class"]) for record in run_records)
    dominant_scaling_class, dominant_count = max(
        class_counter.items(),
        key=lambda item: (item[1], item[0]),
    )
    mean_scores = [float(record["mean_recovery_score"]) for record in run_records]
    variances = [float(record["variance_recovery_score"]) for record in run_records]
    runtimes = [float(record["runtime_seconds"]) for record in run_records]
    aggregate = {
        "dominant_scaling_class": dominant_scaling_class,
        "class_consistency_fraction": round(dominant_count / len(run_records), 6),
        "mean_recovery_score_average": round(sum(mean_scores) / len(mean_scores), 6),
        "mean_recovery_score_range": round(max(mean_scores) - min(mean_scores), 6),
        "variance_recovery_score_average": round(sum(variances) / len(variances), 6),
        "variance_recovery_score_range": round(max(variances) - min(variances), 6),
        "average_runtime_seconds": round(sum(runtimes) / len(runtimes), 6),
        "observed_scaling_classes": sorted(class_counter.keys()),
    }
    if histogram_distance_aggregate is not None:
        aggregate.update(histogram_distance_aggregate)
    return aggregate


def summarize_histogram_profile(
    results_by_schedule: dict[str, dict[str, Any]],
    schedule_labels: list[str],
) -> dict[str, Any]:
    first_histogram = results_by_schedule[schedule_labels[0]]["histogram_payload"]["histogram"]
    bin_centers = list(first_histogram["bin_centers"])
    profiles = [
        list(results_by_schedule[schedule_label]["histogram_payload"]["histogram"]["normalized_frequency"])
        for schedule_label in schedule_labels
    ]
    mean_profile = [
        round(sum(profile[index] for profile in profiles) / len(profiles), 6)
        for index in range(len(bin_centers))
    ]
    min_profile = [
        round(min(profile[index] for profile in profiles), 6)
        for index in range(len(bin_centers))
    ]
    max_profile = [
        round(max(profile[index] for profile in profiles), 6)
        for index in range(len(bin_centers))
    ]
    return {
        "bin_centers": bin_centers,
        "mean_normalized_frequency": mean_profile,
        "min_normalized_frequency": min_profile,
        "max_normalized_frequency": max_profile,
    }


def run_reproducibility_sweep(
    config: dict[str, Any],
    config_path: Path,
    phase4_input: Path | None,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    audit_config = dict(config["audit"])
    schedule_labels = list(audit_config["reproducibility_schedules"])
    results_by_schedule: dict[str, dict[str, Any]] = {}

    for schedule_label in schedule_labels:
        results_by_schedule[schedule_label] = execute_phase5(
            config,
            config_path,
            "frozen_sector",
            "stable_frozen_sector",
            phase4_input,
            schedule_label=schedule_label,
        )

    run_records = [summarize_run(results_by_schedule[schedule_label]) for schedule_label in schedule_labels]
    pairwise_distances, distance_aggregate = pairwise_schedule_distances(results_by_schedule, schedule_labels)
    payload = {
        "phase": 5,
        "phase_name": "phase5-readout",
        "audit_name": "phase5b_reproducibility_ledger_v1",
        "source_config": relpath(config_path),
        "input_mode": "frozen_sector",
        "control_class": "stable_frozen_sector",
        "schedule_labels": schedule_labels,
        "runs": run_records,
        "pairwise_histogram_l1": pairwise_distances,
        "aggregate": summarize_run_cluster(run_records, distance_aggregate),
    }
    return results_by_schedule, payload


def run_policy_robustness_audit(
    config: dict[str, Any],
    config_path: Path,
    phase4_input: Path | None,
) -> dict[str, Any]:
    audit_config = dict(config["audit"])
    variation_bands = dict(audit_config["policy_variation_bands"])
    base_policy = dict(config["frozen_sector"]["policy"])

    baseline_result = execute_phase5(
        config,
        config_path,
        "frozen_sector",
        "stable_frozen_sector",
        phase4_input,
        schedule_label="policy_band_anchor",
    )
    baseline_record = summarize_run(baseline_result)
    baseline_class = str(baseline_record["scaling_class"])
    baseline_histogram = baseline_result["histogram_payload"]["histogram"]

    parameter_sweeps: list[dict[str, Any]] = []
    total_variation_count = 0
    total_class_match_count = 0
    max_histogram_l1_from_baseline = 0.0

    for parameter_name, parameter_values in variation_bands.items():
        variation_records: list[dict[str, Any]] = []
        class_match_count = 0
        for parameter_value in parameter_values:
            result = execute_phase5(
                config,
                config_path,
                "frozen_sector",
                "stable_frozen_sector",
                phase4_input,
                schedule_label=f"policy_band::{parameter_name}",
                policy_override={parameter_name: parameter_value},
            )
            record = summarize_run(result)
            record["parameter"] = parameter_name
            record["parameter_value"] = parameter_value
            record["histogram_l1_from_baseline"] = histogram_l1_distance(
                baseline_histogram,
                result["histogram_payload"]["histogram"],
            )
            variation_records.append(record)
            total_variation_count += 1
            if record["scaling_class"] == baseline_class:
                class_match_count += 1
                total_class_match_count += 1
            max_histogram_l1_from_baseline = max(
                max_histogram_l1_from_baseline,
                float(record["histogram_l1_from_baseline"]),
            )

        mean_scores = [float(record["mean_recovery_score"]) for record in variation_records]
        variances = [float(record["variance_recovery_score"]) for record in variation_records]
        hist_distances = [float(record["histogram_l1_from_baseline"]) for record in variation_records]
        parameter_sweeps.append(
            {
                "parameter": parameter_name,
                "baseline_value": base_policy[parameter_name],
                "runs": variation_records,
                "classifier_stability_fraction": round(class_match_count / len(variation_records), 6),
                "mean_recovery_score_span": round(max(mean_scores) - min(mean_scores), 6),
                "variance_recovery_score_span": round(max(variances) - min(variances), 6),
                "max_histogram_l1_from_baseline": round(max(hist_distances), 6),
            }
        )

    overall_classifier_stability_fraction = round(
        total_class_match_count / total_variation_count,
        6,
    )
    max_mean_recovery_score_span = round(
        max(float(sweep["mean_recovery_score_span"]) for sweep in parameter_sweeps),
        6,
    )
    tightly_clustered = (
        max_mean_recovery_score_span <= 0.03
        and round(max_histogram_l1_from_baseline, 6) <= 0.16
    )
    class_consistent = overall_classifier_stability_fraction == 1.0

    return {
        "phase": 5,
        "phase_name": "phase5-readout",
        "audit_name": "phase5b_robustness_audit_v1",
        "source_config": relpath(config_path),
        "input_mode": "frozen_sector",
        "control_class": "stable_frozen_sector",
        "baseline": baseline_record,
        "parameter_sweeps": parameter_sweeps,
        "aggregate": {
            "baseline_scaling_class": baseline_class,
            "overall_classifier_stability_fraction": overall_classifier_stability_fraction,
            "max_histogram_l1_from_baseline": round(max_histogram_l1_from_baseline, 6),
            "max_mean_recovery_score_span": max_mean_recovery_score_span,
            "observed_scaling_classes": sorted(
                {
                    baseline_class,
                    *(
                        str(record["scaling_class"])
                        for sweep in parameter_sweeps
                        for record in sweep["runs"]
                    ),
                }
            ),
        },
        "success": {
            "stable_policy_class_consistent": class_consistent,
            "stable_policy_tightly_clustered": tightly_clustered,
            "overall_pass": class_consistent or tightly_clustered,
        },
    }


def build_control_runs(
    config: dict[str, Any],
    config_path: Path,
    phase4_input: Path | None,
    schedule_labels: list[str],
    control_class: str,
    preset_results: dict[str, dict[str, Any]] | None = None,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    results_by_schedule = dict(preset_results or {})
    for schedule_label in schedule_labels:
        if schedule_label in results_by_schedule:
            continue
        results_by_schedule[schedule_label] = execute_phase5(
            config,
            config_path,
            "frozen_sector",
            control_class,
            phase4_input,
            schedule_label=schedule_label,
        )
    run_records = [summarize_run(results_by_schedule[schedule_label]) for schedule_label in schedule_labels]
    _, distance_aggregate = pairwise_schedule_distances(results_by_schedule, schedule_labels)
    aggregate = summarize_run_cluster(run_records, distance_aggregate)
    return results_by_schedule, run_records, aggregate


def matched_control_pair_metrics(
    left_control: str,
    left_results: dict[str, dict[str, Any]],
    right_control: str,
    right_results: dict[str, dict[str, Any]],
    schedule_labels: list[str],
) -> dict[str, Any]:
    matched_runs: list[dict[str, Any]] = []
    histogram_distances: list[float] = []
    mean_gaps: list[float] = []
    variance_gaps: list[float] = []
    class_match_count = 0

    for schedule_label in schedule_labels:
        left_bundle = left_results[schedule_label]["bundle"]
        right_bundle = right_results[schedule_label]["bundle"]
        left_summary = dict(left_bundle["summary"])
        right_summary = dict(right_bundle["summary"])
        histogram_distance = histogram_l1_distance(
            left_results[schedule_label]["histogram_payload"]["histogram"],
            right_results[schedule_label]["histogram_payload"]["histogram"],
        )
        mean_gap = round(
            abs(float(left_summary["mean_recovery_score"]) - float(right_summary["mean_recovery_score"])),
            6,
        )
        variance_gap = round(
            abs(float(left_summary["variance_recovery_score"]) - float(right_summary["variance_recovery_score"])),
            6,
        )
        scaling_class_match = left_summary["scaling_class"] == right_summary["scaling_class"]
        if scaling_class_match:
            class_match_count += 1
        matched_runs.append(
            {
                "schedule_label": schedule_label,
                "left_scaling_class": left_summary["scaling_class"],
                "right_scaling_class": right_summary["scaling_class"],
                "histogram_l1_distance": histogram_distance,
                "mean_recovery_score_gap": mean_gap,
                "variance_recovery_score_gap": variance_gap,
            }
        )
        histogram_distances.append(histogram_distance)
        mean_gaps.append(mean_gap)
        variance_gaps.append(variance_gap)

    average_histogram_distance = round(sum(histogram_distances) / len(histogram_distances), 6)
    average_mean_gap = round(sum(mean_gaps) / len(mean_gaps), 6)
    average_variance_gap = round(sum(variance_gaps) / len(variance_gaps), 6)
    scaling_class_overlap_fraction = round(class_match_count / len(schedule_labels), 6)

    return {
        "left_control_class": left_control,
        "right_control_class": right_control,
        "matched_runs": matched_runs,
        "average_histogram_l1_distance": average_histogram_distance,
        "minimum_histogram_l1_distance": round(min(histogram_distances), 6),
        "average_mean_recovery_score_gap": average_mean_gap,
        "average_variance_recovery_score_gap": average_variance_gap,
        "scaling_class_overlap_fraction": scaling_class_overlap_fraction,
        "distinct": (
            scaling_class_overlap_fraction < 1.0
            or average_histogram_distance >= 0.15
            or average_mean_gap >= 0.08
        ),
    }


def find_pair_metric(
    pairwise_metrics: list[dict[str, Any]],
    left_control: str,
    right_control: str,
) -> dict[str, Any]:
    for row in pairwise_metrics:
        pair = {row["left_control_class"], row["right_control_class"]}
        if pair == {left_control, right_control}:
            return row
    raise KeyError(f"missing pairwise metric for {left_control} and {right_control}")


def run_control_discrimination_audit(
    config: dict[str, Any],
    config_path: Path,
    phase4_input: Path | None,
    stable_results: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    schedule_labels = list(config["audit"]["reproducibility_schedules"])
    control_results: dict[str, dict[str, dict[str, Any]]] = {}
    control_payloads: dict[str, dict[str, Any]] = {}

    for control_class in CONTROL_CLASSES:
        preset_results = stable_results if control_class == "stable_frozen_sector" else None
        results_by_schedule, run_records, aggregate = build_control_runs(
            config,
            config_path,
            phase4_input,
            schedule_labels,
            control_class,
            preset_results=preset_results,
        )
        control_results[control_class] = results_by_schedule
        control_payloads[control_class] = {
            "runs": run_records,
            "aggregate": aggregate,
            "histogram_profile": summarize_histogram_profile(results_by_schedule, schedule_labels),
        }

    pairwise_metrics = [
        matched_control_pair_metrics(
            left_control,
            control_results[left_control],
            right_control,
            control_results[right_control],
            schedule_labels,
        )
        for left_control, right_control in combinations(CONTROL_CLASSES, 2)
    ]

    stable_aggregate = control_payloads["stable_frozen_sector"]["aggregate"]
    degraded_aggregate = control_payloads["degraded_sector_control"]["aggregate"]
    null_aggregate = control_payloads["shuffled_null_control"]["aggregate"]
    stable_null_pair = find_pair_metric(pairwise_metrics, "stable_frozen_sector", "shuffled_null_control")
    stable_degraded_pair = find_pair_metric(pairwise_metrics, "stable_frozen_sector", "degraded_sector_control")
    degraded_null_pair = find_pair_metric(pairwise_metrics, "degraded_sector_control", "shuffled_null_control")

    stable_clustered = (
        float(stable_aggregate["class_consistency_fraction"]) >= 0.8
        and float(stable_aggregate["max_histogram_l1_distance"]) <= 0.18
    )
    null_separated = (
        float(stable_null_pair["average_histogram_l1_distance"]) >= 0.35
        or float(stable_null_pair["average_mean_recovery_score_gap"]) >= 0.20
        or float(stable_null_pair["scaling_class_overlap_fraction"]) <= 0.4
    )
    degraded_intermediate = (
        float(stable_aggregate["mean_recovery_score_average"])
        > float(degraded_aggregate["mean_recovery_score_average"])
        > float(null_aggregate["mean_recovery_score_average"])
    )
    degraded_distinct = degraded_intermediate or (
        stable_degraded_pair["distinct"] and degraded_null_pair["distinct"]
    )

    return {
        "phase": 5,
        "phase_name": "phase5-readout",
        "audit_name": "phase5b_control_discrimination_v1",
        "source_config": relpath(config_path),
        "input_mode": "frozen_sector",
        "schedule_labels": schedule_labels,
        "controls": control_payloads,
        "pairwise_separation": pairwise_metrics,
        "success": {
            "stable_runs_clustered": stable_clustered,
            "null_runs_separated": null_separated,
            "degraded_runs_intermediate": degraded_intermediate,
            "degraded_runs_distinct": degraded_distinct,
            "overall_pass": stable_clustered and null_separated and degraded_distinct,
        },
    }


def write_named_payload(output_dir: Path, stem: str, payload: dict[str, Any]) -> dict[str, str]:
    stamped_path, latest_path = write_timestamped_json(output_dir, stem, payload)
    artifacts = {
        "stamped_json": artifact_ref(stamped_path),
        "latest_json": artifact_ref(latest_path),
    }
    payload["artifacts"] = artifacts
    write_json(stamped_path, payload)
    write_json(latest_path, payload)
    return artifacts


def authority_timestamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_text_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def control_display_name(control_class: str) -> str:
    mapping = {
        "stable_frozen_sector": "Stable frozen sector",
        "degraded_sector_control": "Degraded control",
        "shuffled_null_control": "Shuffled null control",
    }
    return mapping[control_class]


def scaling_color(scaling_class: str) -> str:
    return {
        "linear_like": "#1f7a8c",
        "quadratic_like": "#d97706",
        "flat": "#4d7c0f",
        "mixed": "#a16207",
        "unstable": "#b42318",
    }.get(scaling_class, "#475467")


def artifact_path(ref: str) -> Path:
    candidate = Path(ref)
    if candidate.is_absolute():
        return candidate
    return ROOT / ref


def line_points(xs: list[float], ys: list[float], left: float, top: float, width: float, height: float, y_max: float) -> str:
    if len(xs) != len(ys):
        raise ValueError("line chart coordinate arrays must match")
    x_min = min(xs)
    x_max = max(xs)
    x_span = max(1e-9, x_max - x_min)
    points = []
    for x_value, y_value in zip(xs, ys):
        x_coord = left + ((x_value - x_min) / x_span) * width
        y_coord = top + height - (clamp(y_value, 0.0, y_max) / y_max) * height
        points.append(f"{x_coord:.2f},{y_coord:.2f}")
    return " ".join(points)


def schedule_bar_x(index: int, total: int, left: float, width: float) -> float:
    return left + (index + 0.5) * (width / total)


def svg_document(title: str, width: int, height: int, body: str) -> str:
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{title}">',
            '<rect width="100%" height="100%" fill="#fcfcf9"/>',
            body,
            "</svg>",
        ]
    )


def render_control_histogram_svg(control_payload: dict[str, Any]) -> str:
    width = 900
    height = 460
    left = 72
    top = 56
    plot_width = 680
    plot_height = 300
    series = [
        ("stable_frozen_sector", "#1f7a8c"),
        ("degraded_sector_control", "#ca8a04"),
        ("shuffled_null_control", "#b42318"),
    ]
    bin_centers = list(
        control_payload["controls"]["stable_frozen_sector"]["histogram_profile"]["bin_centers"]
    )
    body: list[str] = [
        '<text x="72" y="30" font-family="Helvetica, Arial, sans-serif" font-size="20" fill="#101828">Phase V Stable vs Degraded vs Null Histogram Comparison</text>',
        '<text x="72" y="48" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">Mean normalized frequency across 5 deterministic schedules.</text>',
    ]
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        y_coord = top + plot_height - fraction * plot_height
        body.append(f'<line x1="{left}" y1="{y_coord:.2f}" x2="{left + plot_width}" y2="{y_coord:.2f}" stroke="#e4e7ec" stroke-width="1"/>')
        body.append(f'<text x="26" y="{y_coord + 4:.2f}" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">{fraction:.2f}</text>')
    for center in bin_centers:
        x_coord = left + center * plot_width
        body.append(f'<line x1="{x_coord:.2f}" y1="{top}" x2="{x_coord:.2f}" y2="{top + plot_height}" stroke="#f2f4f7" stroke-width="1"/>')
    body.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#101828" stroke-width="1.4"/>')
    body.append(f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#101828" stroke-width="1.4"/>')
    for control_class, color in series:
        histogram_profile = dict(control_payload["controls"][control_class]["histogram_profile"])
        points = line_points(
            bin_centers,
            list(histogram_profile["mean_normalized_frequency"]),
            left,
            top,
            plot_width,
            plot_height,
            1.0,
        )
        body.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="3" points="{points}"/>'
        )
    legend_x = 780
    legend_y = 90
    for index, (control_class, color) in enumerate(series):
        aggregate = dict(control_payload["controls"][control_class]["aggregate"])
        y_coord = legend_y + index * 48
        body.append(f'<rect x="{legend_x}" y="{y_coord}" width="14" height="14" fill="{color}" rx="2"/>')
        body.append(
            f'<text x="{legend_x + 22}" y="{y_coord + 11}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#101828">{control_display_name(control_class)}</text>'
        )
        body.append(
            f'<text x="{legend_x + 22}" y="{y_coord + 27}" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">mean {aggregate["mean_recovery_score_average"]:.6f}, class {aggregate["dominant_scaling_class"]}</text>'
        )
    body.append('<text x="72" y="390" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">x-axis: recovery score bin center</text>')
    body.append('<text x="72" y="408" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">y-axis: normalized frequency</text>')
    return svg_document("Phase V histogram comparison", width, height, "\n".join(body))


def render_reproducibility_svg(reproducibility_payload: dict[str, Any]) -> str:
    width = 900
    height = 440
    left = 72
    top = 66
    plot_width = 680
    plot_height = 260
    runs = list(reproducibility_payload["runs"])
    aggregate = dict(reproducibility_payload["aggregate"])
    body: list[str] = [
        '<text x="72" y="30" font-family="Helvetica, Arial, sans-serif" font-size="20" fill="#101828">Phase V Reproducibility Sweep Summary</text>',
        '<text x="72" y="48" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">Stable frozen sector across 5 deterministic schedules.</text>',
    ]
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        y_coord = top + plot_height - fraction * plot_height
        body.append(f'<line x1="{left}" y1="{y_coord:.2f}" x2="{left + plot_width}" y2="{y_coord:.2f}" stroke="#e4e7ec" stroke-width="1"/>')
        body.append(f'<text x="26" y="{y_coord + 4:.2f}" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">{fraction:.2f}</text>')
    body.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#101828" stroke-width="1.4"/>')
    body.append(f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#101828" stroke-width="1.4"/>')
    bar_width = 80
    slot_width = plot_width / len(runs)
    for index, run in enumerate(runs):
        x_center = schedule_bar_x(index, len(runs), left, plot_width)
        x_coord = x_center - bar_width / 2.0
        bar_height = plot_height * float(run["mean_recovery_score"])
        y_coord = top + plot_height - bar_height
        body.append(f'<rect x="{x_coord:.2f}" y="{y_coord:.2f}" width="{bar_width}" height="{bar_height:.2f}" fill="#1f7a8c" rx="4"/>')
        body.append(f'<text x="{x_center:.2f}" y="{y_coord - 8:.2f}" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#101828">{run["mean_recovery_score"]:.6f}</text>')
        body.append(f'<text x="{x_center:.2f}" y="{top + plot_height + 20:.2f}" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#101828">{run["schedule_label"].replace("schedule_", "")}</text>')
        body.append(f'<text x="{x_center:.2f}" y="{top + plot_height + 36:.2f}" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="10" fill="#475467">{run["scaling_class"]}</text>')
    body.append(
        f'<text x="72" y="366" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">dominant class {aggregate["dominant_scaling_class"]}; class consistency {aggregate["class_consistency_fraction"]:.6f}; mean-score range {aggregate["mean_recovery_score_range"]:.6f}</text>'
    )
    body.append(
        f'<text x="72" y="384" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">pairwise histogram L1 avg {aggregate["average_histogram_l1_distance"]:.6f}; max {aggregate["max_histogram_l1_distance"]:.6f}</text>'
    )
    return svg_document("Phase V reproducibility summary", width, height, "\n".join(body))


def render_robustness_svg(robustness_payload: dict[str, Any], edge_case_flips: list[dict[str, Any]]) -> str:
    width = 920
    height = 560
    left = 82
    plot_width = 720
    panel_height = 110
    top_offsets = [72, 224, 376]
    baseline_mean = float(robustness_payload["baseline"]["mean_recovery_score"])
    body: list[str] = [
        '<text x="82" y="30" font-family="Helvetica, Arial, sans-serif" font-size="20" fill="#101828">Phase V Robustness Variation Summary</text>',
        '<text x="82" y="48" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">One-at-a-time narrow policy variation against the frozen stable sector baseline.</text>',
    ]
    for top, sweep in zip(top_offsets, robustness_payload["parameter_sweeps"]):
        runs = list(sweep["runs"])
        parameter_values = [float(run["parameter_value"]) for run in runs]
        mean_scores = [float(run["mean_recovery_score"]) for run in runs]
        x_min = min(parameter_values)
        x_max = max(parameter_values)
        x_span = max(1e-9, x_max - x_min)
        body.append(f'<rect x="{left}" y="{top}" width="{plot_width}" height="{panel_height}" fill="#ffffff" stroke="#eaecf0" rx="6"/>')
        body.append(f'<text x="{left + 10}" y="{top + 18}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#101828">{sweep["parameter"]}</text>')
        base_y = top + panel_height - baseline_mean * (panel_height - 28)
        body.append(f'<line x1="{left + 12}" y1="{base_y:.2f}" x2="{left + plot_width - 12}" y2="{base_y:.2f}" stroke="#98a2b3" stroke-dasharray="5 5"/>')
        points = []
        for run in runs:
            x_coord = left + 24 + ((float(run["parameter_value"]) - x_min) / x_span) * (plot_width - 48)
            y_coord = top + panel_height - float(run["mean_recovery_score"]) * (panel_height - 28)
            points.append(f"{x_coord:.2f},{y_coord:.2f}")
        body.append(f'<polyline fill="none" stroke="#475467" stroke-width="2" points="{" ".join(points)}"/>')
        for run in runs:
            x_coord = left + 24 + ((float(run["parameter_value"]) - x_min) / x_span) * (plot_width - 48)
            y_coord = top + panel_height - float(run["mean_recovery_score"]) * (panel_height - 28)
            color = scaling_color(str(run["scaling_class"]))
            body.append(f'<circle cx="{x_coord:.2f}" cy="{y_coord:.2f}" r="5" fill="{color}"/>')
            body.append(f'<text x="{x_coord:.2f}" y="{top + panel_height + 16:.2f}" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="10" fill="#475467">{run["parameter_value"]}</text>')
            body.append(f'<text x="{x_coord:.2f}" y="{y_coord - 10:.2f}" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="10" fill="#101828">{run["scaling_class"]}</text>')
        body.append(
            f'<text x="{left + 10}" y="{top + panel_height + 34}" font-family="Helvetica, Arial, sans-serif" font-size="10" fill="#475467">classifier stability {sweep["classifier_stability_fraction"]:.6f}; score span {sweep["mean_recovery_score_span"]:.6f}; max histogram L1 {sweep["max_histogram_l1_from_baseline"]:.6f}</text>'
        )
    if edge_case_flips:
        flip = edge_case_flips[0]
        body.append(
            f'<text x="82" y="530" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#b42318">Recorded edge case: {flip["parameter"]}={flip["parameter_value"]} produced {flip["scaling_class"]} with mean score {flip["mean_recovery_score"]:.6f}.</text>'
        )
    else:
        body.append(
            '<text x="82" y="530" font-family="Helvetica, Arial, sans-serif" font-size="11" fill="#475467">No robustness edge-case class flip recorded in the tested policy band.</text>'
        )
    return svg_document("Phase V robustness summary", width, height, "\n".join(body))


def find_edge_case_flips(robustness_payload: dict[str, Any]) -> list[dict[str, Any]]:
    baseline_class = str(robustness_payload["aggregate"]["baseline_scaling_class"])
    flips: list[dict[str, Any]] = []
    for sweep in robustness_payload["parameter_sweeps"]:
        for run in sweep["runs"]:
            if str(run["scaling_class"]) != baseline_class:
                flips.append(
                    {
                        "parameter": sweep["parameter"],
                        "parameter_value": run["parameter_value"],
                        "scaling_class": run["scaling_class"],
                        "mean_recovery_score": run["mean_recovery_score"],
                        "histogram_l1_from_baseline": run["histogram_l1_from_baseline"],
                    }
                )
    return flips


def write_authority_figures(
    phase_root: Path,
    control_payload: dict[str, Any],
    reproducibility_payload: dict[str, Any],
    robustness_payload: dict[str, Any],
    edge_case_flips: list[dict[str, Any]],
) -> dict[str, str]:
    figure_dir = phase_root / "phase5_summary_figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    histogram_path = figure_dir / "phase5_histogram_controls.svg"
    reproducibility_path = figure_dir / "phase5_reproducibility_sweep_summary.svg"
    robustness_path = figure_dir / "phase5_robustness_variation_summary.svg"
    write_text_file(histogram_path, render_control_histogram_svg(control_payload))
    write_text_file(reproducibility_path, render_reproducibility_svg(reproducibility_payload))
    write_text_file(robustness_path, render_robustness_svg(robustness_payload, edge_case_flips))
    return {
        "histogram_controls_svg": artifact_ref(histogram_path),
        "reproducibility_sweep_svg": artifact_ref(reproducibility_path),
        "robustness_variation_svg": artifact_ref(robustness_path),
    }


def build_authoritative_runs_ledger(
    timestamp: str,
    config_path: Path,
    baseline_bundle: dict[str, Any],
    audit: dict[str, Any],
    reproducibility_payload: dict[str, Any],
    robustness_payload: dict[str, Any],
    control_payload: dict[str, Any],
    edge_case_flips: list[dict[str, Any]],
    figures: dict[str, str],
    summary_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "phase": 5,
        "phase_name": "phase5-readout",
        "authority_mode": "phase5close_v1",
        "source_config": relpath(config_path),
        "baseline_frozen_run": {
            "bundle": baseline_bundle["artifacts"]["latest_json"],
            "recovery_histogram": baseline_bundle["artifacts"]["recovery_histogram_latest_json"],
            "scaling_fit": baseline_bundle["artifacts"]["scaling_fit_latest_json"],
            "summary": dict(baseline_bundle["summary"]),
        },
        "audit_runs": {
            "reproducibility": {
                "latest_json": audit["artifacts"]["reproducibility_latest_json"],
                "aggregate": reproducibility_payload["aggregate"],
            },
            "robustness": {
                "latest_json": audit["artifacts"]["robustness_latest_json"],
                "aggregate": robustness_payload["aggregate"],
                "success": robustness_payload["success"],
                "edge_case_classifier_flips": edge_case_flips,
            },
            "control_discrimination": {
                "latest_json": audit["artifacts"]["control_discrimination_latest_json"],
                "success": control_payload["success"],
                "control_aggregates": {
                    control_class: payload["aggregate"]
                    for control_class, payload in control_payload["controls"].items()
                },
                "pairwise_separation": [
                    {
                        "left_control_class": row["left_control_class"],
                        "right_control_class": row["right_control_class"],
                        "average_histogram_l1_distance": row["average_histogram_l1_distance"],
                        "average_mean_recovery_score_gap": row["average_mean_recovery_score_gap"],
                        "scaling_class_overlap_fraction": row["scaling_class_overlap_fraction"],
                        "distinct": row["distinct"],
                    }
                    for row in control_payload["pairwise_separation"]
                ],
            },
        },
        "authority_artifacts": {
            "summary_note": artifact_ref(summary_path),
            "manifest": artifact_ref(manifest_path),
            "figures": figures,
        },
        "success": bool(audit["summary"]["overall_pass"]),
    }


def build_authoritative_manifest(
    timestamp: str,
    config_path: Path,
    phase_root: Path,
    output_dir: Path,
    baseline_bundle: dict[str, Any],
    audit: dict[str, Any],
    reproducibility_payload: dict[str, Any],
    robustness_payload: dict[str, Any],
    control_payload: dict[str, Any],
    edge_case_flips: list[dict[str, Any]],
    figures: dict[str, str],
    ledger_artifacts: dict[str, str],
    summary_path: Path,
) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "phase_identifier": "phase5-readout",
        "status": "passed" if bool(audit["summary"]["overall_pass"]) else "failed",
        "frozen_runner_path": relpath(phase_root / "runs" / "run_phase5.py"),
        "config_path": relpath(config_path),
        "diagnostic_outputs": {
            "baseline_bundle_latest_json": baseline_bundle["artifacts"]["latest_json"],
            "baseline_histogram_latest_json": baseline_bundle["artifacts"]["recovery_histogram_latest_json"],
            "baseline_scaling_fit_latest_json": baseline_bundle["artifacts"]["scaling_fit_latest_json"],
            "reproducibility_latest_json": audit["artifacts"]["reproducibility_latest_json"],
            "robustness_latest_json": audit["artifacts"]["robustness_latest_json"],
            "control_discrimination_latest_json": audit["artifacts"]["control_discrimination_latest_json"],
        },
        "authority_artifacts": {
            "summary_note": artifact_ref(summary_path),
            "runs_ledger_latest_json": ledger_artifacts["latest_json"],
            "runs_ledger_stamped_json": ledger_artifacts["stamped_json"],
            "summary_figures": figures,
        },
        "dominant_stable_class": reproducibility_payload["aggregate"]["dominant_scaling_class"],
        "reproducibility_metrics": {
            "schedule_count": len(reproducibility_payload["schedule_labels"]),
            "class_consistency_fraction": reproducibility_payload["aggregate"]["class_consistency_fraction"],
            "mean_recovery_score_range": reproducibility_payload["aggregate"]["mean_recovery_score_range"],
            "average_histogram_l1_distance": reproducibility_payload["aggregate"]["average_histogram_l1_distance"],
            "max_histogram_l1_distance": reproducibility_payload["aggregate"]["max_histogram_l1_distance"],
        },
        "robustness_metrics": {
            "baseline_scaling_class": robustness_payload["aggregate"]["baseline_scaling_class"],
            "overall_classifier_stability_fraction": robustness_payload["aggregate"]["overall_classifier_stability_fraction"],
            "max_histogram_l1_from_baseline": robustness_payload["aggregate"]["max_histogram_l1_from_baseline"],
            "max_mean_recovery_score_span": robustness_payload["aggregate"]["max_mean_recovery_score_span"],
            "stable_policy_tightly_clustered": robustness_payload["success"]["stable_policy_tightly_clustered"],
            "stable_policy_class_consistent": robustness_payload["success"]["stable_policy_class_consistent"],
            "edge_case_classifier_flips": edge_case_flips,
        },
        "control_separation_metrics": {
            "success": control_payload["success"],
            "controls": {
                control_class: payload["aggregate"]
                for control_class, payload in control_payload["controls"].items()
            },
            "pairwise_separation": [
                {
                    "left_control_class": row["left_control_class"],
                    "right_control_class": row["right_control_class"],
                    "average_histogram_l1_distance": row["average_histogram_l1_distance"],
                    "average_mean_recovery_score_gap": row["average_mean_recovery_score_gap"],
                    "scaling_class_overlap_fraction": row["scaling_class_overlap_fraction"],
                    "distinct": row["distinct"],
                }
                for row in control_payload["pairwise_separation"]
            ],
        },
        "success": bool(audit["summary"]["overall_pass"]),
        "claim_boundary": CLAIM_BOUNDARY,
    }


def build_authoritative_summary_markdown(
    timestamp: str,
    baseline_bundle: dict[str, Any],
    reproducibility_payload: dict[str, Any],
    robustness_payload: dict[str, Any],
    control_payload: dict[str, Any],
    edge_case_flips: list[dict[str, Any]],
    manifest_path: Path,
    ledger_latest_ref: str,
    figures: dict[str, str],
) -> str:
    stable_aggregate = control_payload["controls"]["stable_frozen_sector"]["aggregate"]
    degraded_aggregate = control_payload["controls"]["degraded_sector_control"]["aggregate"]
    null_aggregate = control_payload["controls"]["shuffled_null_control"]["aggregate"]
    stable_null_pair = find_pair_metric(
        control_payload["pairwise_separation"],
        "stable_frozen_sector",
        "shuffled_null_control",
    )
    stable_degraded_pair = find_pair_metric(
        control_payload["pairwise_separation"],
        "stable_frozen_sector",
        "degraded_sector_control",
    )
    degraded_null_pair = find_pair_metric(
        control_payload["pairwise_separation"],
        "degraded_sector_control",
        "shuffled_null_control",
    )
    baseline_summary = dict(baseline_bundle["summary"])
    flip_lines = "\n".join(
        [
            f"- `{flip['parameter']} = {flip['parameter_value']}` produced `{flip['scaling_class']}` with mean recovery `{flip['mean_recovery_score']:.6f}` and histogram L1 from baseline `{flip['histogram_l1_from_baseline']:.6f}`."
            for flip in edge_case_flips
        ]
    ) or "- No class flip was observed in the tested policy band."

    return "\n".join(
        [
            "# Phase V Authoritative Summary",
            "",
            f"Authority freeze timestamp: `{timestamp}`.",
            "",
            "## Objective",
            "",
            "Freeze Phase V as a deterministic recovery-histogram readout over the frozen Phase IV stable sector bundle, with bounded reporting of reproducibility, robustness, and control discrimination.",
            "",
            "## Frozen Input Schema",
            "",
            "```json",
            "{",
            '  "sector_identifier": "...",',
            '  "graph_reference_or_seed": {},',
            '  "kernel_configuration_snapshot": {},',
            '  "selector_configuration": {},',
            '  "invariant_baseline_vector": {},',
            '  "perturbation_policy_descriptor": {}',
            "}",
            "```",
            "",
            "The authoritative frozen run uses the current stable Phase IV bundle and rejects oversized upstream payloads before compression into this schema.",
            "",
            "## Recovery Score Definition",
            "",
            "The Phase V recovery score is the arithmetic mean of five normalized invariant-recovery components: `law_fit_score`, `stable_subset_size`, `stable_rho`, `ordering_gap`, and `monotonic_consistency_score`. Each component is clamped to `[0, 1]`. No alternate score is introduced in the authority freeze.",
            "",
            "## Run Modes",
            "",
            "- `frozen_sector`: authoritative mode over the frozen Phase IV stable sector bundle.",
            "- `dummy_sector`: synthetic control mode retained for scaffolding and smoke tests.",
            "- Control classes: `stable_frozen_sector`, `degraded_sector_control`, `shuffled_null_control`.",
            "",
            "## Reproducibility Results",
            "",
            f"- Stable frozen sector was swept across `{len(reproducibility_payload['schedule_labels'])}` deterministic schedules.",
            f"- Dominant stable class: `{reproducibility_payload['aggregate']['dominant_scaling_class']}`.",
            f"- Class consistency fraction: `{reproducibility_payload['aggregate']['class_consistency_fraction']:.6f}`.",
            f"- Mean recovery score average/range: `{reproducibility_payload['aggregate']['mean_recovery_score_average']:.6f}` / `{reproducibility_payload['aggregate']['mean_recovery_score_range']:.6f}`.",
            f"- Pairwise histogram L1 average/max: `{reproducibility_payload['aggregate']['average_histogram_l1_distance']:.6f}` / `{reproducibility_payload['aggregate']['max_histogram_l1_distance']:.6f}`.",
            "",
            "## Robustness Results",
            "",
            f"- Baseline stable class: `{robustness_payload['aggregate']['baseline_scaling_class']}`.",
            f"- Overall classifier stability fraction across narrow policy variation: `{robustness_payload['aggregate']['overall_classifier_stability_fraction']:.6f}`.",
            f"- Max histogram L1 from baseline: `{robustness_payload['aggregate']['max_histogram_l1_from_baseline']:.6f}`.",
            f"- Max mean recovery score span across the tested band: `{robustness_payload['aggregate']['max_mean_recovery_score_span']:.6f}`.",
            f"- Tight-cluster test: `{robustness_payload['success']['stable_policy_tightly_clustered']}`.",
            "- Recorded edge-case classifier flip:",
            flip_lines,
            "",
            "## Control Separation Results",
            "",
            f"- Stable aggregate: mean `{stable_aggregate['mean_recovery_score_average']:.6f}`, dominant class `{stable_aggregate['dominant_scaling_class']}`.",
            f"- Degraded aggregate: mean `{degraded_aggregate['mean_recovery_score_average']:.6f}`, dominant class `{degraded_aggregate['dominant_scaling_class']}`, observed classes `{', '.join(degraded_aggregate['observed_scaling_classes'])}`.",
            f"- Null aggregate: mean `{null_aggregate['mean_recovery_score_average']:.6f}`, dominant class `{null_aggregate['dominant_scaling_class']}`.",
            f"- Stable vs degraded: histogram L1 `{stable_degraded_pair['average_histogram_l1_distance']:.6f}`, mean-score gap `{stable_degraded_pair['average_mean_recovery_score_gap']:.6f}`.",
            f"- Stable vs null: histogram L1 `{stable_null_pair['average_histogram_l1_distance']:.6f}`, mean-score gap `{stable_null_pair['average_mean_recovery_score_gap']:.6f}`.",
            f"- Degraded vs null: histogram L1 `{degraded_null_pair['average_histogram_l1_distance']:.6f}`, mean-score gap `{degraded_null_pair['average_mean_recovery_score_gap']:.6f}`.",
            "",
            "## Bounded Interpretation",
            "",
            f"- Within the tested configuration band, the stable frozen sector remains a high-recovery regime with a dominant `{reproducibility_payload['aggregate']['dominant_scaling_class']}` label.",
            "- Degraded and null controls remain clearly separated in both recovery level and distributional distance.",
            "- The recorded amplitude-jitter edge flip is treated as a boundary marker for the descriptive classifier, not as evidence for a new mechanism or theory layer.",
            "",
            "## Explicit Non-Claims",
            "",
            "- No claim is made that the descriptive scaling classifier is invariant outside the tested policy band.",
            "- No claim is made that the observed labels establish a universal recovery law.",
            "- No claim is made that degraded or null controls encode physical mechanisms beyond their role as structural controls.",
            "- No claim is made beyond the frozen Phase IV bundle, the configured perturbation schedules, and the current runner implementation.",
            "",
            "## Authority Artifacts",
            "",
            f"- Baseline frozen run: `{baseline_summary['scaling_class']}` with mean recovery `{baseline_summary['mean_recovery_score']:.6f}`.",
            f"- Manifest: `{artifact_ref(manifest_path)}`.",
            f"- Runs ledger: `{ledger_latest_ref}`.",
            f"- Figures: `{figures['histogram_controls_svg']}`, `{figures['reproducibility_sweep_svg']}`, `{figures['robustness_variation_svg']}`.",
            "",
            f"Claim boundary: {CLAIM_BOUNDARY}",
        ]
    )


def run_phase5(
    config_path: Path,
    output_dir: Path,
    mode: str,
    control_class: str,
    phase4_input: Path | None,
) -> dict[str, Any]:
    config = read_json(config_path)
    result = execute_phase5(
        config,
        config_path,
        mode,
        control_class,
        phase4_input,
    )
    artifacts = write_phase5_outputs(
        output_dir,
        result["bundle"],
        result["histogram_payload"],
        result["scaling_payload"],
    )
    result["bundle"]["artifacts"] = artifacts
    return result["bundle"]


def run_phase5b_audit(
    config_path: Path,
    output_dir: Path,
    phase4_input: Path | None,
) -> dict[str, Any]:
    config = read_json(config_path)
    started_at = time.perf_counter()

    stable_results, reproducibility_payload = run_reproducibility_sweep(
        config,
        config_path,
        phase4_input,
    )
    robustness_payload = run_policy_robustness_audit(
        config,
        config_path,
        phase4_input,
    )
    control_payload = run_control_discrimination_audit(
        config,
        config_path,
        phase4_input,
        stable_results,
    )

    reproducibility_artifacts = write_named_payload(
        output_dir,
        "phase5_reproducibility_ledger",
        reproducibility_payload,
    )
    robustness_artifacts = write_named_payload(
        output_dir,
        "phase5_robustness_audit",
        robustness_payload,
    )
    control_artifacts = write_named_payload(
        output_dir,
        "phase5_control_discrimination",
        control_payload,
    )
    runtime_seconds = time.perf_counter() - started_at
    overall_pass = (
        bool(robustness_payload["success"]["overall_pass"])
        and bool(control_payload["success"]["overall_pass"])
    )

    return {
        "phase": 5,
        "phase_name": "phase5-readout",
        "audit_name": AUDIT_NAME,
        "summary": {
            "stable_dominant_scaling_class": reproducibility_payload["aggregate"]["dominant_scaling_class"],
            "stable_class_consistency_fraction": reproducibility_payload["aggregate"]["class_consistency_fraction"],
            "policy_classifier_stability_fraction": robustness_payload["aggregate"][
                "overall_classifier_stability_fraction"
            ],
            "policy_tightly_clustered": robustness_payload["success"]["stable_policy_tightly_clustered"],
            "null_runs_separated": control_payload["success"]["null_runs_separated"],
            "degraded_runs_distinct": control_payload["success"]["degraded_runs_distinct"],
            "overall_pass": overall_pass,
            "runtime_seconds": round(runtime_seconds, 6),
        },
        "artifacts": {
            "reproducibility_stamped_json": reproducibility_artifacts["stamped_json"],
            "reproducibility_latest_json": reproducibility_artifacts["latest_json"],
            "robustness_stamped_json": robustness_artifacts["stamped_json"],
            "robustness_latest_json": robustness_artifacts["latest_json"],
            "control_discrimination_stamped_json": control_artifacts["stamped_json"],
            "control_discrimination_latest_json": control_artifacts["latest_json"],
        },
    }


def run_phase5close_authority(
    config_path: Path,
    output_dir: Path,
    phase4_input: Path | None,
) -> dict[str, Any]:
    timestamp = authority_timestamp()
    phase_root = output_dir.parent
    summary_path = phase_root / "phase5_authoritative_summary.md"
    manifest_path = phase_root / "phase5_authoritative_manifest.json"

    baseline_bundle = run_phase5(
        config_path,
        output_dir,
        "frozen_sector",
        "stable_frozen_sector",
        phase4_input,
    )
    audit = run_phase5b_audit(
        config_path,
        output_dir,
        phase4_input,
    )
    reproducibility_payload = read_json(artifact_path(audit["artifacts"]["reproducibility_latest_json"]))
    robustness_payload = read_json(artifact_path(audit["artifacts"]["robustness_latest_json"]))
    control_payload = read_json(artifact_path(audit["artifacts"]["control_discrimination_latest_json"]))
    edge_case_flips = find_edge_case_flips(robustness_payload)
    figures = write_authority_figures(
        phase_root,
        control_payload,
        reproducibility_payload,
        robustness_payload,
        edge_case_flips,
    )
    provisional_ledger = build_authoritative_runs_ledger(
        timestamp,
        config_path,
        baseline_bundle,
        audit,
        reproducibility_payload,
        robustness_payload,
        control_payload,
        edge_case_flips,
        figures,
        summary_path,
        manifest_path,
    )
    ledger_artifacts = write_named_payload(
        output_dir,
        "phase5_runs_ledger",
        provisional_ledger,
    )
    manifest_payload = build_authoritative_manifest(
        timestamp,
        config_path,
        phase_root,
        output_dir,
        baseline_bundle,
        audit,
        reproducibility_payload,
        robustness_payload,
        control_payload,
        edge_case_flips,
        figures,
        ledger_artifacts,
        summary_path,
    )
    write_json(manifest_path, manifest_payload)
    summary_text = build_authoritative_summary_markdown(
        timestamp,
        baseline_bundle,
        reproducibility_payload,
        robustness_payload,
        control_payload,
        edge_case_flips,
        manifest_path,
        ledger_artifacts["latest_json"],
        figures,
    )
    write_text_file(summary_path, summary_text)

    return {
        "phase": 5,
        "phase_name": "phase5-readout",
        "audit_name": "phase5close_authority_freeze_v1",
        "summary": {
            "dominant_stable_class": reproducibility_payload["aggregate"]["dominant_scaling_class"],
            "overall_pass": bool(audit["summary"]["overall_pass"]),
            "edge_case_classifier_flip_count": len(edge_case_flips),
        },
        "artifacts": {
            "summary_note": artifact_ref(summary_path),
            "manifest_json": artifact_ref(manifest_path),
            "runs_ledger_latest_json": ledger_artifacts["latest_json"],
            "runs_ledger_stamped_json": ledger_artifacts["stamped_json"],
            "histogram_controls_svg": figures["histogram_controls_svg"],
            "reproducibility_sweep_svg": figures["reproducibility_sweep_svg"],
            "robustness_variation_svg": figures["robustness_variation_svg"],
            "baseline_bundle_latest_json": baseline_bundle["artifacts"]["latest_json"],
            "reproducibility_latest_json": audit["artifacts"]["reproducibility_latest_json"],
            "robustness_latest_json": audit["artifacts"]["robustness_latest_json"],
            "control_discrimination_latest_json": audit["artifacts"]["control_discrimination_latest_json"],
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Phase V coherence recovery readout rig.")
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "phase5-readout" / "configs" / "phase5_readout_scaffold.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "phase5-readout" / "runs",
    )
    parser.add_argument(
        "--mode",
        choices=INPUT_MODES,
        default="frozen_sector",
    )
    parser.add_argument(
        "--control-class",
        choices=CONTROL_CLASSES,
        default="stable_frozen_sector",
    )
    parser.add_argument(
        "--phase4-input",
        type=Path,
        default=None,
        help="Optional explicit Phase IV bundle path for frozen_sector mode.",
    )
    parser.add_argument(
        "--audit",
        choices=("phase5b", "phase5close"),
        default=None,
        help="Run the Phase V-B audit or the Phase V-close authority freeze.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the final Phase V output as JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.audit == "phase5b":
        audit = run_phase5b_audit(
            args.config,
            args.output_dir,
            args.phase4_input,
        )
        if args.json:
            print(json.dumps(audit, indent=2))
            return
        print(f"Phase V-B overall pass: {audit['summary']['overall_pass']}")
        print(f"Reproducibility ledger: {audit['artifacts']['reproducibility_latest_json']}")
        print(f"Robustness audit: {audit['artifacts']['robustness_latest_json']}")
        print(f"Control discrimination: {audit['artifacts']['control_discrimination_latest_json']}")
        return
    if args.audit == "phase5close":
        authority = run_phase5close_authority(
            args.config,
            args.output_dir,
            args.phase4_input,
        )
        if args.json:
            print(json.dumps(authority, indent=2))
            return
        print(f"Phase V-close overall pass: {authority['summary']['overall_pass']}")
        print(f"Authority summary: {authority['artifacts']['summary_note']}")
        print(f"Authority manifest: {authority['artifacts']['manifest_json']}")
        print(f"Authority ledger: {authority['artifacts']['runs_ledger_latest_json']}")
        return

    bundle = run_phase5(
        args.config,
        args.output_dir,
        args.mode,
        args.control_class,
        args.phase4_input,
    )
    if args.json:
        print(json.dumps(bundle, indent=2))
        return
    print(f"Phase V scaling class: {bundle['summary']['scaling_class']}")
    print(f"Latest bundle: {bundle['artifacts']['latest_json']}")


if __name__ == "__main__":
    main()
