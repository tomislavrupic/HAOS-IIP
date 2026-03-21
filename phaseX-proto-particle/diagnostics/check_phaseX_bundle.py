#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phaseX-proto-particle"


def main() -> None:
    required = [
        PHASE_ROOT / "runs" / "proto_particle_candidates.json",
        PHASE_ROOT / "runs" / "localization_metrics_ledger.csv",
        PHASE_ROOT / "runs" / "persistence_classification.csv",
        PHASE_ROOT / "runs" / "scaling_prediction_errors.csv",
        PHASE_ROOT / "phaseX_integrated_manifest.json",
        PHASE_ROOT / "phaseX_integrated_summary.md",
        PHASE_ROOT / "runs" / "phaseX_runs.json",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    manifest = json.loads((PHASE_ROOT / "phaseX_integrated_manifest.json").read_text(encoding="utf-8"))
    print(
        json.dumps(
            {
                "missing": missing,
                "success": bool(manifest.get("success")),
                "conclusion": manifest.get("conclusion"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
