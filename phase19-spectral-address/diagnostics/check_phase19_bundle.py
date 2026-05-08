#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


PHASE_DIR = Path(__file__).resolve().parents[1]
RUNS_DIR = PHASE_DIR / "runs"

EXPECTED_FILES = [
    PHASE_DIR / "spectral_address_definition.md",
    PHASE_DIR / "slide_alignment.md",
    PHASE_DIR / "phase19_summary.md",
    PHASE_DIR / "phase19_manifest.json",
    RUNS_DIR / "phase19_runs.json",
]

REQUIRED_PHRASES = {
    "definition": [
        "stable, recoverable resonance pattern",
        "frozen branch-local cochain-Laplacian hierarchy",
        "persistent spectral signature",
        "interaction invariants",
        "recoverable cochain structures",
    ],
    "math": [
        "Laplacian eigenmodes",
        "spectral graph theory",
        "cochain complexes",
        "Kuramoto",
        "attractor",
    ],
}


def main() -> None:
    missing = [str(path) for path in EXPECTED_FILES if not path.exists()]
    if missing:
        raise SystemExit(json.dumps({"success": False, "missing_files": missing}, indent=2))

    text = "\n".join(path.read_text() for path in EXPECTED_FILES if path.suffix == ".md")
    missing_phrases = [
        phrase
        for group in REQUIRED_PHRASES.values()
        for phrase in group
        if phrase not in text
    ]
    if missing_phrases:
        raise SystemExit(json.dumps({"success": False, "missing_phrases": missing_phrases}, indent=2))

    manifest = json.loads((PHASE_DIR / "phase19_manifest.json").read_text())
    missing_inputs = [
        path
        for path in manifest.get("input_files", {}).values()
        if not Path(path).exists()
    ]
    if missing_inputs:
        raise SystemExit(json.dumps({"success": False, "missing_inputs": missing_inputs}, indent=2))

    gates = manifest.get("gates", {})
    failed_gates = [name for name, value in gates.items() if value is not True]
    if failed_gates:
        raise SystemExit(json.dumps({"success": False, "failed_gates": failed_gates}, indent=2))

    print(json.dumps({"success": True, "gates": gates}, indent=2))


if __name__ == "__main__":
    main()
