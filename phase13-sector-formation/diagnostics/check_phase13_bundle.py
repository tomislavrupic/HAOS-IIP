#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phase13-sector-formation"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

REQUIRED_PATHS = [
    RUNS_ROOT / "phase13_population_survival_ledger.csv",
    RUNS_ROOT / "phase13_spacing_statistics_ledger.csv",
    RUNS_ROOT / "phase13_cluster_metrics.csv",
    RUNS_ROOT / "phase13_spectral_ensemble_proxy.csv",
    RUNS_ROOT / "phase13_runs.json",
    PHASE_ROOT / "phase13_summary.md",
    PHASE_ROOT / "phase13_manifest.json",
    PLOTS_ROOT / "phase13_survival_fraction_vs_time.svg",
    PLOTS_ROOT / "phase13_pair_distance_histogram_vs_refinement.svg",
    PLOTS_ROOT / "phase13_nearest_neighbor_scale_vs_refinement.svg",
    PLOTS_ROOT / "phase13_cluster_count_vs_time.svg",
    PLOTS_ROOT / "phase13_ensemble_spectral_proxy_vs_population.svg",
]

SUCCESS_SENTENCE = "Phase XIII establishes multi-mode statistical sector formation feasibility for the frozen operator hierarchy."
FAILURE_SENTENCE = "Phase XIII does not yet establish multi-mode statistical sector formation feasibility for the frozen operator hierarchy."


def main() -> None:
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED_PATHS if not path.exists()]
    success = False
    conclusion = "missing_artifacts"
    if not missing:
        manifest = json.loads((PHASE_ROOT / "phase13_manifest.json").read_text(encoding="utf-8"))
        lines = (PHASE_ROOT / "phase13_summary.md").read_text(encoding="utf-8").strip().splitlines()
        summary_closure = lines[-1] if lines else ""
        success = bool(manifest.get("success")) and summary_closure in {SUCCESS_SENTENCE, FAILURE_SENTENCE}
        conclusion = str(manifest.get("conclusion", summary_closure))
    print(json.dumps({"missing": missing, "success": success, "conclusion": conclusion}, indent=2))


if __name__ == "__main__":
    main()
