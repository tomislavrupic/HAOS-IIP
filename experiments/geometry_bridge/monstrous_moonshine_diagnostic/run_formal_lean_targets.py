#!/usr/bin/env python3
"""Formal Lean target scaffold for the Betti graph diagnostic.

The repository does not currently contain a local Lean project. This runner
therefore records Lean theorem targets and required definitions without
claiming Lean certification. It links the targets to executable Python evidence
from the Betti vector diagnostic and keeps the claim ceiling explicit.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_betti_vector import BettiVectorConfig, run_betti_vector  # noqa: E402
from run_monstrous_moonshine_diagnostic import stable_hash, write_json  # noqa: E402


MANIFEST_PATH = ROOT / "formal_lean_targets_manifest.json"
REPORT_PATH = ROOT / "formal_lean_targets_report.md"
TARGET_PATH = ROOT / "lean_graph_invariance_targets.lean"


@dataclass(frozen=True)
class FormalLeanTargetConfig:
    version: str = "formal-lean-targets-v0.1"
    lean_project_present: bool = False
    lean_check_run: bool = False
    claim_ceiling: str = "FORMAL_TARGET_SCAFFOLD_ONLY"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write formal Lean target scaffold for Betti graph invariance.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def required_definitions() -> list[dict[str, str]]:
    return [
        {
            "name": "GaussianPrimeNode",
            "purpose": "arithmetic graph vertex carrying prime and Gaussian representative coordinates",
            "local_lean_status": "MISSING",
            "python_source": "run_betti_component_count.py:base_nodes",
        },
        {
            "name": "RelationGraph",
            "purpose": "finite undirected graph over declared Gaussian-prime representatives",
            "local_lean_status": "MISSING",
            "python_source": "run_betti_component_count.py:relation_edges",
        },
        {
            "name": "D4Action",
            "purpose": "square-lattice rotation/reflection action on Gaussian-prime coordinates",
            "local_lean_status": "MISSING",
            "python_source": "run_betti_component_count.py:transform_nodes",
        },
        {
            "name": "edgeRelation",
            "purpose": "Euclidean threshold relation used to build the graph",
            "local_lean_status": "MISSING",
            "python_source": "run_betti_component_count.py:relation_edges",
        },
        {
            "name": "componentCount",
            "purpose": "finite graph connected-component count",
            "local_lean_status": "MISSING",
            "python_source": "run_betti_component_count.py:component_count",
        },
        {
            "name": "bettiOne",
            "purpose": "finite graph cycle count E - V + componentCount",
            "local_lean_status": "MISSING",
            "python_source": "run_betti_vector.py:betti_signature",
        },
    ]


def theorem_ladder() -> list[dict[str, str]]:
    return [
        {
            "id": "L1",
            "name": "finite_graph_component_count_exists",
            "target_statement": "componentCount is defined for every finite RelationGraph",
            "status": "TARGET_ONLY_NOT_LEAN_CHECKED",
        },
        {
            "id": "L2",
            "name": "graph_iso_preserves_component_count",
            "target_statement": "graph isomorphism preserves componentCount",
            "status": "TARGET_ONLY_NOT_LEAN_CHECKED",
        },
        {
            "id": "L3",
            "name": "d4_action_induces_graph_iso",
            "target_statement": "D4Action induces graph isomorphism under the Euclidean threshold edgeRelation",
            "status": "TARGET_ONLY_NOT_LEAN_CHECKED",
        },
        {
            "id": "L4",
            "name": "d4_preserves_component_count",
            "target_statement": "componentCount (applyD4 s g) = componentCount g",
            "status": "TARGET_ONLY_NOT_LEAN_CHECKED",
        },
        {
            "id": "L5",
            "name": "graph_iso_preserves_betti_one",
            "target_statement": "graph isomorphism preserves bettiOne",
            "status": "TARGET_ONLY_NOT_LEAN_CHECKED",
        },
        {
            "id": "L6",
            "name": "d4_preserves_betti_vector",
            "target_statement": "D4Action preserves [Betti_0, Betti_1]",
            "status": "TARGET_ONLY_NOT_LEAN_CHECKED",
        },
    ]


def target_source() -> str:
    return """namespace HAOSIIP.BettiGraphTargets

/-!
Formal target scaffold only.

This repository does not currently include a Lean project, local graph
definitions, or checked proofs for these targets. The theorem statements below
are intentionally comments, not axioms, not placeholder proofs, and not
certification.

Required future definitions:
- GaussianPrimeNode
- RelationGraph
- D4Action
- edgeRelation
- componentCount
- bettiOne

Target ladder:
- theorem finite_graph_component_count_exists
- theorem graph_iso_preserves_component_count
- theorem d4_action_induces_graph_iso
- theorem d4_preserves_component_count
- theorem graph_iso_preserves_betti_one
- theorem d4_preserves_betti_vector
-/

def targetStatus : String := "TARGET_ONLY_NOT_LEAN_CHECKED"
def claimCeiling : String := "FORMAL_TARGET_SCAFFOLD_ONLY"

end HAOSIIP.BettiGraphTargets
"""


def run_formal_targets(config: FormalLeanTargetConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    betti_result = run_betti_vector(BettiVectorConfig(), output_dir)
    stable_conditions = ["known_positive", "d4_rotate_90", "d4_reflect_real", "isomorphism_relabel"]
    executable_evidence = {
        "betti_vector_result_hash": betti_result["result_hash"],
        "reference_signature": betti_result["reference_signature"],
        "stable_condition_status": {
            condition: betti_result["condition_results"][condition]["observed_status"] for condition in stable_conditions
        },
        "stable_condition_betti_vector_distance": {
            condition: betti_result["condition_results"][condition]["betti_vector_distance"] for condition in stable_conditions
        },
        "stable_condition_edge_distance": {
            condition: betti_result["condition_results"][condition]["edge_jaccard_distance"] for condition in stable_conditions
        },
    }

    labels = [
        "LEAN_TARGET_SCAFFOLD_BUILT",
        "GRAPH_ISO_COMPONENT_COUNT_TARGET_DECLARED",
        "D4_COMPONENT_COUNT_TARGET_DECLARED",
        "BETTI_VECTOR_INVARIANCE_TARGET_DECLARED",
        "LOCAL_LEAN_PROJECT_NOT_PRESENT",
        "LEAN_CHECK_NOT_RUN",
        "LEAN_THEOREM_NOT_INCLUDED",
        "MOONSHINE_PROOF_NOT_ESTABLISHED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_FORMAL_TARGET_ONLY",
    ]
    result = {
        "bridge_id": "formal_lean_targets_v0_1",
        "status": "OPEN",
        "classification": "FORMAL_TARGET_SCAFFOLD_ONLY",
        "labels": labels,
        "config": asdict(config),
        "required_definitions": required_definitions(),
        "theorem_ladder": theorem_ladder(),
        "executable_evidence": executable_evidence,
        "claim_ceiling": "Lean target scaffold only; no local Lean proof, Moonshine proof, or physical bridge claim.",
        "non_claims": [
            "not Lean-certified",
            "not a checked theorem inside this repository",
            "not a Moonshine proof",
            "not a Monster action",
            "not a physical bridge",
        ],
    }
    result["result_hash"] = stable_hash("formal_lean_targets", result)

    write_json(output_dir / MANIFEST_PATH.name, result)
    (output_dir / TARGET_PATH.name).write_text(target_source(), encoding="utf-8")
    write_report(output_dir / REPORT_PATH.name, result)
    return result


def write_report(path: Path, result: dict[str, Any]) -> None:
    lines = [
        "# Formal Lean Target Scaffold Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        f"- Betti evidence hash: `{result['executable_evidence']['betti_vector_result_hash']}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Required Local Lean Definitions",
        "",
        "| Name | Status | Executable source |",
        "| --- | --- | --- |",
    ]
    for definition in result["required_definitions"]:
        lines.append(
            f"| `{definition['name']}` | `{definition['local_lean_status']}` | `{definition['python_source']}` |"
        )
    lines.extend(["", "## Target Ladder", "", "| ID | Target | Status |", "| --- | --- | --- |"])
    for theorem in result["theorem_ladder"]:
        lines.append(f"| `{theorem['id']}` | `{theorem['name']}` | `{theorem['status']}` |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This is a formal target scaffold, not a Lean-certified theorem.",
            "The checked evidence remains executable Python diagnostics over the declared graph.",
            "The future Lean work should start with graph isomorphism preserving component count before any D4 specialization.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_formal_targets(FormalLeanTargetConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
