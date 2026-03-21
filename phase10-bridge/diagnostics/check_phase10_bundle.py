#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phase10-bridge"

REQUIRED = [
    PHASE_ROOT / "phase10_manifest.json",
    PHASE_ROOT / "phase10_summary.md",
    PHASE_ROOT / "runs" / "phase10_extension_ledger.csv",
    PHASE_ROOT / "runs" / "phase10_effective_scaling_ledger.csv",
    PHASE_ROOT / "runs" / "phase10_coarse_grain_summary.csv",
    PHASE_ROOT / "runs" / "phase10_runs.json",
    PHASE_ROOT / "plots" / "phase10_descriptor_vs_extended_refinement.svg",
    PHASE_ROOT / "plots" / "phase10_prediction_error_vs_refinement.svg",
    PHASE_ROOT / "plots" / "phase10_coarse_grain_spectral_comparison.svg",
    PHASE_ROOT / "plots" / "phase10_multi_descriptor_correlation.svg",
]


def main() -> None:
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED if not path.exists()]
    manifest_path = PHASE_ROOT / "phase10_manifest.json"
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
