#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phase12-interactions"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

REQUIRED_PATHS = [
    RUNS_ROOT / "phase12_interaction_outcome_ledger.csv",
    RUNS_ROOT / "phase12_identity_metrics_ledger.csv",
    RUNS_ROOT / "phase12_interaction_thresholds.csv",
    RUNS_ROOT / "phase12_runs.json",
    PHASE_ROOT / "phase12_summary.md",
    PHASE_ROOT / "phase12_manifest.json",
    PLOTS_ROOT / "phase12_persistence_time_vs_separation.svg",
    PLOTS_ROOT / "phase12_identity_metric_vs_time.svg",
    PLOTS_ROOT / "phase12_interaction_regime_map.svg",
    PLOTS_ROOT / "phase12_threshold_scaling_vs_refinement.svg",
]

SUCCESS_SENTENCE = "Phase XII establishes two-mode interaction and identity-preservation feasibility for the frozen operator hierarchy."
FAILURE_SENTENCE = "Phase XII does not yet establish two-mode interaction and identity-preservation feasibility for the frozen operator hierarchy."


def main() -> None:
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED_PATHS if not path.exists()]
    success = False
    conclusion = "missing_artifacts"
    if not missing:
        manifest = json.loads((PHASE_ROOT / "phase12_manifest.json").read_text(encoding="utf-8"))
        summary = (PHASE_ROOT / "phase12_summary.md").read_text(encoding="utf-8").strip().splitlines()
        summary_closure = summary[-1] if summary else ""
        success = bool(manifest.get("success")) and summary_closure in {SUCCESS_SENTENCE, FAILURE_SENTENCE}
        conclusion = str(manifest.get("conclusion", summary_closure))
    print(
        json.dumps(
            {
                "missing": missing,
                "success": success,
                "conclusion": conclusion,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
