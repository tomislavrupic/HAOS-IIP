#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phase11-protection"

REQUIRED = [
    PHASE_ROOT / "phase11_manifest.json",
    PHASE_ROOT / "phase11_summary.md",
    PHASE_ROOT / "runs" / "phase11_perturbation_survival_ledger.csv",
    PHASE_ROOT / "runs" / "phase11_failure_classification.csv",
    PHASE_ROOT / "runs" / "phase11_persistence_scaling.csv",
    PHASE_ROOT / "runs" / "phase11_runs.json",
    PHASE_ROOT / "plots" / "phase11_localization_width_vs_time.svg",
    PHASE_ROOT / "plots" / "phase11_survival_probability_vs_perturbation.svg",
    PHASE_ROOT / "plots" / "phase11_persistence_time_vs_refinement.svg",
    PHASE_ROOT / "plots" / "phase11_spectral_gap_shift_vs_defect_strength.svg"
]


def main() -> None:
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED if not path.exists()]
    manifest_path = PHASE_ROOT / "phase11_manifest.json"
    success = False
    conclusion = None
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        success = bool(manifest.get("success"))
        conclusion = manifest.get("conclusion")
    print(
        json.dumps(
            {
                "missing": missing,
                "success": bool(success and not missing),
                "conclusion": conclusion,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
