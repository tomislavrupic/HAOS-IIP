#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def detect_green_family(payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": "scalar_green_response",
        "operator_family": "kernel graph Laplacian",
        "substrate_family": "3D cubic kernel graph",
        "geometry_target": "inverse-distance field / effective 3D Green geometry",
        "dimension_hint": 3,
        "source_path": "data/kernel_graph_green_response_latest.json",
        "observation": str(payload.get("observation", "")),
        "conclusion": str(payload.get("conclusion", "")),
    }


def detect_trace_family(manifest: dict[str, Any]) -> dict[str, Any]:
    trace_method = str(manifest["trace_proxy_definition"]["label"])
    eval_method = str(manifest["trace_proxy_definition"]["evaluation_method"])
    return {
        "label": "heat_behavior",
        "operator_family": trace_method,
        "substrate_family": "periodic DK2D branch-local cochain-Laplacian hierarchy",
        "geometry_target": "short-time heat / trace asymptotic behavior",
        "dimension_hint": 2,
        "source_path": "phase8-trace/phase8_trace_manifest.json",
        "observation": eval_method,
        "conclusion": str(manifest.get("conclusion", "")),
    }


def detect_low_mode_family(manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": "low_mode_organization",
        "operator_family": "branch-local cochain-Laplacian spectral hierarchy",
        "substrate_family": "periodic DK2D branch-local cochain-Laplacian hierarchy",
        "geometry_target": "low-mode organization / spectral feasibility",
        "dimension_hint": 2,
        "source_path": "phase7-spectral/phase7_spectral_manifest.json",
        "observation": str(manifest.get("claim_boundary", "")),
        "conclusion": str(manifest.get("conclusion", "")),
    }


def detect_surrogate_family(manifest: dict[str, Any], refinement_rows: list[dict[str, str]]) -> dict[str, Any]:
    branch_rows = [row for row in refinement_rows if row["hierarchy_label"] == "frozen_branch"]
    slopes = [float(row["mean_arrival_vs_depth_slope"]) for row in branch_rows]
    return {
        "label": "shell_arrival_surrogate",
        "operator_family": "frozen branch-local cochain-Laplacian hierarchy with reused XV-XVII ledgers",
        "substrate_family": "periodic DK2D branch-local cochain-Laplacian hierarchy",
        "geometry_target": "shell-arrival / causal-depth distance surrogate",
        "dimension_hint": 2,
        "source_path": "phase18-distance-surrogate/phase18_manifest.json",
        "observation": f"branch slope band={min(slopes):.6f}-{max(slopes):.6f}" if slopes else "branch slope band unavailable",
        "conclusion": str(manifest.get("closure_statement", "")),
    }


def compatibility_matrix(families: list[dict[str, Any]]) -> dict[str, Any]:
    substrate_families = {item["substrate_family"] for item in families}
    operator_families = {item["operator_family"] for item in families}
    dimension_hints = {item["dimension_hint"] for item in families}
    shared_geometry_ready = len(substrate_families) == 1 and len(operator_families) == 1 and len(dimension_hints) == 1
    return {
        "unique_substrate_families": sorted(substrate_families),
        "unique_operator_families": sorted(operator_families),
        "unique_dimension_hints": sorted(dimension_hints),
        "shared_geometry_ready": bool(shared_geometry_ready),
    }


def make_result() -> dict[str, Any]:
    green_payload = read_json(REPO_ROOT / "data" / "kernel_graph_green_response_latest.json")
    phase7_manifest = read_json(REPO_ROOT / "phase7-spectral" / "phase7_spectral_manifest.json")
    phase8_manifest = read_json(REPO_ROOT / "phase8-trace" / "phase8_trace_manifest.json")
    phase18_manifest = read_json(REPO_ROOT / "phase18-distance-surrogate" / "phase18_manifest.json")
    phase18_refinement = read_csv_rows(REPO_ROOT / "phase18-distance-surrogate" / "runs" / "phase18_refinement_scaling.csv")

    families = [
        detect_green_family(green_payload),
        detect_trace_family(phase8_manifest),
        detect_low_mode_family(phase7_manifest),
        detect_surrogate_family(phase18_manifest, phase18_refinement),
    ]
    compatibility = compatibility_matrix(families)

    observation = (
        "the strict geometry-closure criterion is correct, but the current strongest ingredients still split into two incompatible authority families: "
        "a 3D cubic kernel-graph Green-response line and a 2D periodic branch-local torus line for heat, low-mode, and shell-arrival diagnostics"
    )
    if compatibility["shared_geometry_ready"]:
        conclusion = (
            "the current repo is already preflight-compatible for a single effective-geometry closure test: all selected ingredients live on one operator and substrate family"
        )
    else:
        conclusion = (
            "geometry closure is not yet preflight-ready: one effective geometry cannot be claimed until Green response, heat behavior, low-mode organization, "
            "and shell-arrival surrogates are rebuilt on the same operator/substrate family instead of being compared across the current 3D kernel-graph and 2D branch-local torus split"
        )

    return {
        "experiment": "geometry_closure_preflight",
        "families": families,
        "compatibility": compatibility,
        "observation": observation,
        "conclusion": conclusion,
    }


def save_result(result: dict[str, Any]) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = DATA / f"{stamp}_geometry_closure_preflight.json"
    path.write_text(json.dumps(result, indent=2))
    return path


def write_note(result: dict[str, Any], result_path: Path) -> Path:
    note_path = REPO_ROOT / "docs" / "notes" / "foundations" / "HAOS_IIP_Geometry_Closure_Preflight_v1.md"
    rows = "\n".join(
        f"| `{item['label']}` | {item['operator_family']} | {item['substrate_family']} | {item['dimension_hint']} | {item['geometry_target']} |"
        for item in result["families"]
    )
    note = f"""# HAOS-IIP Geometry Closure Preflight v1

## Purpose

Check whether the current strongest geometry-facing ingredients in the repo are even commensurate enough to support the strict geometry-closure success criterion:

> one effective geometry must organize Green response, heat behavior, shell-arrival slopes, and low-mode organization on the same operator/substrate family.

## Selected authority ingredients

| channel | operator family | substrate family | dimension hint | target |
| --- | --- | --- | ---: | --- |
{rows}

## Direct read

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Compatibility summary

- shared geometry ready: `{result['compatibility']['shared_geometry_ready']}`
- substrate families: `{result['compatibility']['unique_substrate_families']}`
- operator families: `{result['compatibility']['unique_operator_families']}`
- dimension hints: `{result['compatibility']['unique_dimension_hints']}`

## What this means

The next geometry-closure tranche should **not** claim success by stitching together the current best 3D Green-response result and the current best 2D branch-local heat / low-mode / shell-arrival stack.

The honest next move is a shared-family rebuild. In practice that means one of:

1. move Green response, heat behavior, low-mode organization, and shell-arrival diagnostics onto one common scalar kernel-graph family;
2. or rebuild the distance-surrogate / shell-arrival side on the same family as the operator side being used for Green and heat.

Only after that shared-family bridge exists does the strict success condition become a live theorem-like numerical target instead of a category mistake.

## Artifact

- result: `{result_path.relative_to(REPO_ROOT)}`
"""
    note_path.write_text(note, encoding="utf-8")
    return note_path


def append_log(result_path: Path, observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## geometry closure preflight\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write("- Config: authority_inputs=['kernel_graph_green_response_latest.json', 'phase7_spectral_manifest.json', 'phase8_trace_manifest.json', 'phase18_manifest.json']\n")
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def main() -> None:
    result = make_result()
    result_path = save_result(result)
    note_path = write_note(result, result_path)
    append_log(result_path, result["observation"], result["conclusion"])
    print(
        json.dumps(
            {
                "result_path": str(result_path),
                "note_path": str(note_path),
                "observation": result["observation"],
                "conclusion": result["conclusion"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
