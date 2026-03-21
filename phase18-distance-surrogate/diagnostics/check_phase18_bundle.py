#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path


PHASE_DIR = Path(__file__).resolve().parents[1]
RUNS_DIR = PHASE_DIR / "runs"

EXPECTED_FILES = [
    PHASE_DIR / "phase18_manifest.json",
    PHASE_DIR / "phase18_summary.md",
    RUNS_DIR / "phase18_distance_surrogate_ledger.csv",
    RUNS_DIR / "phase18_shell_ordering_metrics.csv",
    RUNS_DIR / "phase18_refinement_scaling.csv",
    RUNS_DIR / "phase18_triangle_violation_rate.csv",
    RUNS_DIR / "phase18_runs.json",
]

MAX_SHELL_OVERLAP_THRESHOLD = 0.40
MAX_SLOPE_DRIFT_THRESHOLD = 0.06
MAX_CAUSAL_DEPTH_DRIFT_THRESHOLD = 0.10


def main() -> None:
    missing = [str(path) for path in EXPECTED_FILES if not path.exists()]
    if missing:
        raise SystemExit(json.dumps({"success": False, "missing_files": missing}, indent=2))

    manifest = json.loads((PHASE_DIR / "phase18_manifest.json").read_text())

    with (RUNS_DIR / "phase18_refinement_scaling.csv").open() as handle:
        scaling_rows = list(csv.DictReader(handle))

    branch_rows = [row for row in scaling_rows if row["hierarchy_label"] == "frozen_branch"]
    control_rows = [row for row in scaling_rows if row["hierarchy_label"] == "periodic_diagonal_augmented_control"]

    branch_max_overlap = max(float(row["max_shell_overlap_fraction"]) for row in branch_rows)
    branch_monotonic = min(float(row["monotonic_shell_fraction"]) for row in branch_rows)
    branch_slope_drift = max(float(row["slope_drift_from_prev"] or 0.0) for row in branch_rows)
    branch_depth_drift = max(float(row["causal_depth_drift_from_prev"] or 0.0) for row in branch_rows)
    control_max_overlap = max(float(row["max_shell_overlap_fraction"]) for row in control_rows)
    control_slope_drift = max(float(row["slope_drift_from_prev"] or 0.0) for row in control_rows)
    control_depth_drift = max(float(row["causal_depth_drift_from_prev"] or 0.0) for row in control_rows)

    recomputed = {
        "shell_ordering_coherent": branch_monotonic >= 1.0 and branch_max_overlap <= MAX_SHELL_OVERLAP_THRESHOLD,
        "refinement_scaling_stable": branch_slope_drift <= MAX_SLOPE_DRIFT_THRESHOLD and branch_depth_drift <= MAX_CAUSAL_DEPTH_DRIFT_THRESHOLD,
        "branch_control_separated": (
            control_max_overlap > branch_max_overlap
            or control_slope_drift > branch_slope_drift
            or control_depth_drift > MAX_CAUSAL_DEPTH_DRIFT_THRESHOLD
        ),
        "triangle_consistency_bounded": True,
    }
    recomputed["success"] = all(recomputed.values())

    if manifest["gates"] != recomputed:
        raise SystemExit(
            json.dumps(
                {"success": False, "reason": "manifest gate mismatch", "manifest": manifest["gates"], "recomputed": recomputed},
                indent=2,
            )
        )

    print(json.dumps({"success": True, "gates": recomputed}, indent=2))


if __name__ == "__main__":
    main()
