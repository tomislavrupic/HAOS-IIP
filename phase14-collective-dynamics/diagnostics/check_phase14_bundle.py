#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phase14-collective-dynamics"

REQUIRED = [
    PHASE_ROOT / "phase14_manifest.json",
    PHASE_ROOT / "phase14_summary.md",
    PHASE_ROOT / "runs" / "phase14_transport_ledger.csv",
    PHASE_ROOT / "runs" / "phase14_relaxation_times.csv",
    PHASE_ROOT / "runs" / "phase14_fluctuation_spectra.csv",
    PHASE_ROOT / "runs" / "phase14_density_response.csv",
    PHASE_ROOT / "runs" / "phase14_equation_of_state_proxy.csv",
    PHASE_ROOT / "runs" / "phase14_runs.json",
]

SUCCESS_SENTENCE = "Phase XIV establishes collective dynamics feasibility for the frozen operator hierarchy."
FAILURE_SENTENCE = "Phase XIV does not yet establish collective dynamics feasibility for the frozen operator hierarchy."


def main() -> None:
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED if not path.exists()]
    manifest_path = PHASE_ROOT / "phase14_manifest.json"
    summary_path = PHASE_ROOT / "phase14_summary.md"
    success = False
    conclusion = FAILURE_SENTENCE
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        success = bool(manifest.get("success"))
        conclusion = str(manifest.get("conclusion", conclusion))
    if summary_path.exists():
        summary = summary_path.read_text(encoding="utf-8")
        if SUCCESS_SENTENCE not in summary and FAILURE_SENTENCE not in summary:
            conclusion = "summary_closure_missing"
            success = False
    payload = {
        "missing": missing,
        "success": bool(success and not missing and conclusion == SUCCESS_SENTENCE),
        "conclusion": conclusion,
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
