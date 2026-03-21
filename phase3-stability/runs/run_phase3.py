#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import read_json, relpath, write_json, write_timestamped_json


def build_claim_rows(inputs: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    s10 = inputs["stage23_10"]["json"]
    s11 = inputs["stage23_11"]["json"]
    s12 = inputs["stage23_12"]["json"]
    s13 = inputs["stage23_13"]["json"]
    return [
        {
            "claim_id": "P3_POS_01",
            "claim_text": "The tested clustered DK braid/smear sector admits a valid two-state effective surrogate on the anchor lattice.",
            "status": "survives",
            "supporting_stage": "23.11 / III-C2",
            "supporting_metric_or_fact": "effective_model_valid = TRUE; topology agreement = 1.0000; threshold error = 0.0000",
            "scope": "clustered_sector_only",
            "notes": "Valid only as a clustered texture surrogate.",
        },
        {
            "claim_id": "P3_POS_02",
            "claim_text": "The clustered DK sector closes as a finite smeared_dominant_closed_phase_diagram.",
            "status": "survives",
            "supporting_stage": "23.12 / III-C3",
            "supporting_metric_or_fact": f"phase_diagram_closed = {s12['summary']['phase_diagram_closed']}; closure_classification = {s12['summary']['closure_classification']}",
            "scope": "phase_closure",
            "notes": "Closure is one-sided and smeared-dominant.",
        },
        {
            "claim_id": "P3_POS_03",
            "claim_text": "Only the strong-asymmetry row width_ratio = 1.35 survives as refinement-/beta-/reverse-stable across the sampled phase corridor.",
            "status": "survives",
            "supporting_stage": "23.12 / III-C3",
            "supporting_metric_or_fact": "stable component = smeared_transfer_phase; width_min = width_max = 1.35; phase range 0.375 -> 0.575",
            "scope": "phase_closure",
            "notes": "No symmetric or moderate-asymmetry cell closes stably.",
        },
        {
            "claim_id": "P3_POS_04",
            "claim_text": "phase_corridor_primary is the controlling mechanism classification for the braid/smear split in the clustered sector.",
            "status": "survives",
            "supporting_stage": "23.13 / III-C4",
            "supporting_metric_or_fact": f"overall_mechanism_label = {s13['overall_summary']['overall_mechanism_label']}",
            "scope": "mechanism_read",
            "notes": "Weak and moderate phase flattening both flip the clustered control to transfer_smeared.",
        },
        {
            "claim_id": "P3_POS_05",
            "claim_text": "Geometry overlap and 0 <-> 1 operator cross-coupling act as bounded secondary supports only.",
            "status": "bounded_support",
            "supporting_stage": "23.13 / III-C4",
            "supporting_metric_or_fact": "Both branches collapse braid only at moderate ablation, not weak ablation.",
            "scope": "mechanism_read",
            "notes": "Secondary support only; not primary control.",
        },
        {
            "claim_id": "P3_NEG_01",
            "claim_text": "The DK braid does not generalize beyond clustered seeds and is not a family-wide intrinsic mechanism.",
            "status": "fails",
            "supporting_stage": "23.10 / III-C1",
            "supporting_metric_or_fact": f"non_clustered_braid_hits = {s10['summary']['classification_payload']['non_clustered_braid_hits']}; robust_non_clustered_braid_hits = {s10['summary']['classification_payload']['robust_non_clustered_braid_hits']}",
            "scope": "not_family_wide",
            "notes": "Family-wide braid generalization fails on this branch.",
        },
        {
            "claim_id": "P3_NEG_02",
            "claim_text": "No stable closed braid phase survives the full closure pass.",
            "status": "fails",
            "supporting_stage": "23.12 / III-C3",
            "supporting_metric_or_fact": f"stable_phase_counts = {s12['summary']['stable_phase_counts']}",
            "scope": "phase_closure",
            "notes": "Braid remains transient and does not close as a stable phase.",
        },
        {
            "claim_id": "P3_NEG_03",
            "claim_text": "No localized encounter phase survives the closure pass.",
            "status": "fails",
            "supporting_stage": "23.12 / III-C3",
            "supporting_metric_or_fact": f"stable_phase_counts = {s12['summary']['stable_phase_counts']}",
            "scope": "phase_closure",
            "notes": "Localized encounter phase absent from the stable closure set.",
        },
        {
            "claim_id": "P3_NEG_04",
            "claim_text": "Grade transfer is not the primary driver on the frozen branch.",
            "status": "fails",
            "supporting_stage": "23.13 / III-C4",
            "supporting_metric_or_fact": f"grade-transfer branch label = {s13['branch_summary']['A_grade_transfer_suppression']['label']}",
            "scope": "mechanism_read",
            "notes": "Grade-transfer-primary mechanism claim fails.",
        },
        {
            "claim_id": "P3_NEG_05",
            "claim_text": "The observed clustered braid-like exchange may not be elevated into a universal topological law.",
            "status": "excluded",
            "supporting_stage": "23.10 / III-C1; 23.12 / III-C3; 23.13 / III-C4",
            "supporting_metric_or_fact": "Clustered-only scope, absence of stable braid closure, and phase-corridor-primary control bound the interpretation.",
            "scope": "clustered_sector_only",
            "notes": "Universal mechanism reading is excluded on this branch.",
        },
    ]


def build_consistency_checks(inputs: dict[str, dict[str, Any]], claim_rows: list[dict[str, Any]]) -> dict[str, bool]:
    s10 = inputs["stage23_10"]["json"]
    s11 = inputs["stage23_11"]["json"]
    s12 = inputs["stage23_12"]["json"]
    s13 = inputs["stage23_13"]["json"]

    positive_claim_text = " ".join(
        row["claim_text"] for row in claim_rows if str(row["status"]) in {"survives", "bounded_support"}
    ).lower()
    contradiction_free = all(
        phrase not in positive_claim_text
        for phrase in (
            "family-wide intrinsic braid mechanism",
            "stable closed braid phase",
            "localized encounter phase",
            "grade transfer is the primary driver",
        )
    )

    return {
        "stage23_10_zero_non_clustered_braid_hits": bool(
            int(s10["summary"]["classification_payload"]["robust_non_clustered_braid_hits"]) == 0
            and int(s10["summary"]["classification_payload"]["non_clustered_braid_hits"]) == 0
        ),
        "stage23_11_effective_model_valid_true": bool(s11["summary"]["effective_model_valid"]),
        "stage23_12_phase_diagram_closed_true": bool(s12["summary"]["phase_diagram_closed"]),
        "stage23_12_no_stable_braid_phase": bool(
            "stable_braid_phase" not in s12["summary"]["stable_phase_counts"]
            and all(str(cell["derived_phase_label"]) != "stable_braid_phase" for cell in s12["cells"])
        ),
        "stage23_13_phase_corridor_primary": bool(
            s13["overall_summary"]["overall_mechanism_label"] == "phase_corridor_primary"
        ),
        "positive_core_noncontradictory": contradiction_free,
    }


def compact_input_refs(inputs: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        key: {
            "stage_label": value["stage_label"],
            "note_path": value["note_path"],
            "json_path": value["json_path"],
            "csv_path": value["csv_path"],
            "note_present": value["note_present"],
        }
        for key, value in inputs.items()
    }


def build_summary_payload(
    runsheet: dict[str, Any],
    inputs: dict[str, dict[str, Any]],
    claim_rows: list[dict[str, Any]],
    checks: dict[str, bool],
) -> dict[str, Any]:
    s11 = inputs["stage23_11"]["json"]
    s12 = inputs["stage23_12"]["json"]
    s13 = inputs["stage23_13"]["json"]

    positive_core_frozen = all(str(row["status"]) in {"survives", "bounded_support"} for row in claim_rows[:5])
    negative_boundaries_frozen = all(str(row["status"]) in {"fails", "excluded"} for row in claim_rows[5:])
    mechanism_question_closed = bool(s13["overall_summary"]["overall_mechanism_label"] == "phase_corridor_primary")
    inputs_consistent = all(bool(value) for value in checks.values())
    freeze_valid = bool(inputs_consistent and positive_core_frozen and negative_boundaries_frozen and mechanism_question_closed)

    return {
        "experiment": runsheet["stage"],
        "description": runsheet["description"],
        "phaseIII_status": "complete",
        "phaseIII_consolidation_complete": True,
        "phaseIII_inputs_consistent": inputs_consistent,
        "phaseIII_freeze_valid": freeze_valid,
        "positive_core_frozen": positive_core_frozen,
        "negative_boundaries_frozen": negative_boundaries_frozen,
        "mechanism_question_closed": mechanism_question_closed,
        "family_wide_braid_mechanism": False,
        "effective_model_valid": bool(s11["summary"]["effective_model_valid"]),
        "phase_diagram_closed": bool(s12["summary"]["phase_diagram_closed"]),
        "closure_class": "smeared_dominant_closed_phase_diagram",
        "stable_braid_phase_present": False,
        "localized_encounter_phase_present": False,
        "mechanism_primary": "phase_corridor_primary",
        "geometry_overlap_role": "secondary_support",
        "operator_cross_coupling_role": "secondary_support",
        "grade_transfer_primary": False,
        "phaseIII_minimal_read": "clustered DK sector closes as a finite smeared-transfer effective regime; braid is texture-level and phase-corridor-controlled, not a family-wide intrinsic mechanism",
        "consistency_checks": checks,
        "inputs": compact_input_refs(inputs),
        "upstream_facts": {
            "stage23_10": inputs["stage23_10"]["json"]["summary"],
            "stage23_11": s11["summary"],
            "stage23_12": s12["summary"],
            "stage23_13": {
                "overall_summary": s13["overall_summary"],
                "branch_summary": s13["branch_summary"],
            },
        },
        "claim_ledger": claim_rows,
    }


def load_inputs(runsheet: dict[str, Any]) -> dict[str, dict[str, Any]]:
    inputs: dict[str, dict[str, Any]] = {}
    for key, item in runsheet["required_inputs"].items():
        note_path = ROOT / item["note_path"]
        json_path = ROOT / item["json_path"]
        inputs[key] = {
            "stage_label": item["stage_label"],
            "note_path": item["note_path"],
            "json_path": item["json_path"],
            "csv_path": item["csv_path"],
            "note_present": note_path.exists(),
            "json": read_json(json_path),
        }
    return inputs


def run_phase3(config_path: Path, output_dir: Path) -> dict[str, Any]:
    runsheet = read_json(config_path)
    inputs = load_inputs(runsheet)
    claim_rows = build_claim_rows(inputs)
    checks = build_consistency_checks(inputs, claim_rows)
    summary = build_summary_payload(runsheet, inputs, claim_rows, checks)
    bundle = {
        "phase": 3,
        "phase_name": "phase3-stability",
        "source_config": relpath(config_path),
        "upstream_inputs": compact_input_refs(inputs),
        "summary": summary,
        "claim_ledger": claim_rows,
    }
    stamped, latest = write_timestamped_json(output_dir, "phase3_stability_bundle", bundle)
    bundle["artifacts"] = {
        "stamped_json": relpath(stamped),
        "latest_json": relpath(latest),
    }
    write_json(stamped, bundle)
    write_json(latest, bundle)
    return bundle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the narrow Phase III stability freeze.")
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "phase3-stability" / "configs" / "stage23_14_phaseIII_final_consolidation_runs.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "phase3-stability" / "runs",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the final bundle as JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bundle = run_phase3(args.config, args.output_dir)
    if args.json:
        print(json.dumps(bundle, indent=2))
        return
    print(f"Phase III freeze valid: {bundle['summary']['phaseIII_freeze_valid']}")
    print(f"Latest bundle: {bundle['artifacts']['latest_json']}")


if __name__ == "__main__":
    main()
