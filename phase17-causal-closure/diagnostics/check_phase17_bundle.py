#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


PHASE_DIR = Path(__file__).resolve().parent.parent


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    manifest_path = PHASE_DIR / "phase17_manifest.json"
    summary_path = PHASE_DIR / "phase17_summary.md"
    runs_path = PHASE_DIR / "runs" / "phase17_runs.json"
    manifest = load_json(manifest_path)
    runs = load_json(runs_path)

    required_files = [
        manifest_path,
        summary_path,
        PHASE_DIR / "runs" / "phase17_influence_graph_ledger.csv",
        PHASE_DIR / "runs" / "phase17_dag_statistics.csv",
        PHASE_DIR / "runs" / "phase17_causal_distance_metrics.csv",
        PHASE_DIR / "runs" / "phase17_propagation_order_compatibility.csv",
        runs_path,
    ]

    success = True
    messages: list[str] = []

    for path in required_files:
        if not path.exists():
            success = False
            messages.append(f"missing_required_file:{path.name}")

    selected_slice = manifest.get("selected_slice", {})
    if selected_slice.get("probe_name") != "bias_onset":
        success = False
        messages.append("probe_name_must_remain_bias_onset")
    if selected_slice.get("refinements") != [60, 72, 84]:
        success = False
        messages.append("refinement_slice_must_remain_60_72_84")
    if selected_slice.get("ensemble_size") != 7:
        success = False
        messages.append("ensemble_size_must_remain_7")

    threshold_policy = runs.get("threshold_policy", {})
    expected_threshold_policy = {
        "radius_fraction": 0.5,
        "dispersion_fraction": 0.5,
        "spectral_fraction": 0.5,
        "low_k_fraction": 0.5,
        "width_fraction": 0.5,
    }
    if threshold_policy != expected_threshold_policy:
        success = False
        messages.append("threshold_policy_drifted_from_phase16")

    branch_flags = manifest.get("success_flags", {})
    control_flags = manifest.get("control_flags", {})
    manifest_success = bool(manifest.get("success"))
    expected_success = all(branch_flags.values()) and any(not value for value in control_flags.values())
    if manifest_success != expected_success:
        success = False
        messages.append("manifest_success_inconsistent_with_gate_logic")

    if manifest_success:
        for artifact_key in ("acyclicity_plot", "mismatch_plot", "depth_plot"):
            artifact_value = manifest["artifacts"].get(artifact_key)
            if not artifact_value:
                success = False
                messages.append(f"missing_plot_reference:{artifact_key}")
            elif not (PHASE_DIR.parent / artifact_value).exists():
                success = False
                messages.append(f"missing_plot_file:{artifact_key}")

    summary_text = summary_path.read_text(encoding="utf-8") if summary_path.exists() else ""
    if manifest.get("conclusion") and manifest["conclusion"] not in summary_text:
        success = False
        messages.append("summary_missing_conclusion")

    output = {
        "success": success,
        "manifest_success": manifest_success,
        "messages": messages,
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
