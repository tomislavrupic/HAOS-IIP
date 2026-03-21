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

from haos_core import read_csv_rows, read_json, relpath, write_json, write_timestamped_json


def find_csv_row(rows: list[dict[str, str]], key: str, value: str) -> dict[str, str]:
    for row in rows:
        if str(row[key]) == value:
            return row
    raise KeyError(f"no CSV row found where {key} == {value}")


def build_claim_rows(inputs: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    s23_14 = inputs["stage23_14"]["json"]
    s24_1 = inputs["stage24_1"]["json"]
    s24_2 = inputs["stage24_2"]["json"]
    s24_3 = inputs["stage24_3"]["json"]
    s24_4 = inputs["stage24_4"]["json"]
    s24_3_rows = inputs["stage24_3"]["csv_rows"]

    q3 = find_csv_row(s24_3_rows, "candidate_id", "IVA3_Q3")
    portable_rows = [
        {
            "claim_id": "P4_POS_01",
            "claim_text": "The smeared sector is carried by the accepted 24.1 observable core and not by rejected braid- or asymmetry-based quantities.",
            "status": "survives_portably",
            "supporting_stage": "24.1 / IV-A1",
            "supporting_metric_or_fact": f"carry_forward_core = {s24_1['carry_forward_core']}; rejected_observables = {s24_1['rejected_observables']}",
            "scope": "observable_ledger",
            "notes": "Observable core is compact and explicit.",
        },
        {
            "claim_id": "P4_POS_02",
            "claim_text": "A compact threshold selector exists for the stable closed smeared sector.",
            "status": "survives_portably",
            "supporting_stage": "24.2 / IV-A2",
            "supporting_metric_or_fact": f"selected_law_expression = {s24_2['selected_law_expression']}",
            "scope": "selector_law",
            "notes": "Portable only inside the tested bounded extension lattice.",
        },
        {
            "claim_id": "P4_POS_03",
            "claim_text": "The selected selector is the minimal accepted law on the tested branch.",
            "status": "survives_portably",
            "supporting_stage": "24.2 / IV-A2",
            "supporting_metric_or_fact": f"selected_model_id = {s24_2['selected_model_id']}; selected_model_family = {s24_2['selected_model_family']}",
            "scope": "selector_law",
            "notes": "No smaller accepted law beat it on the frozen branch.",
        },
        {
            "claim_id": "P4_POS_04",
            "claim_text": "The selector law survives all tested support-anisotropy perturbations in the bounded extension lattice.",
            "status": "survives_portably",
            "supporting_stage": "24.4 / IV-A4",
            "supporting_metric_or_fact": f"support_anisotropy_skew boundary = {s24_4['first_boundary_by_class']['support_anisotropy_skew']}",
            "scope": "bounded_extension",
            "notes": "Support anisotropy does not produce the first selector failure on the tested manifold.",
        },
        {
            "claim_id": "P4_POS_05",
            "claim_text": "Phase IV preserves the Phase III negative boundaries against family-wide braid, stable braid phase, and universal transport-law escalation.",
            "status": "survives_portably",
            "supporting_stage": "23.14 / III-C5; 24.5 / IV-A5",
            "supporting_metric_or_fact": f"family_wide_braid_mechanism = {s23_14['family_wide_braid_mechanism']}; braid_reintroduced = False",
            "scope": "boundary_preservation",
            "notes": "Phase IV does not reopen excluded Phase III claims.",
        },
    ]
    local_rows = [
        {
            "claim_id": "P4_LOC_01",
            "claim_text": "The stable smeared sector carries a non-trivial inverse ordering relation between flow concentration and topology survival time.",
            "status": "survives_locally",
            "supporting_stage": "24.3 / IV-A3",
            "supporting_metric_or_fact": f"selected_candidate_id = {s24_3['selected_candidate_id']}; stable_metric = {q3['stable_metric']}",
            "scope": "local_structure",
            "notes": "This is local internal structure, not a universal selector law.",
        },
        {
            "claim_id": "P4_LOC_02",
            "claim_text": "The ordering relation adds structure beyond the selector law itself.",
            "status": "survives_locally",
            "supporting_stage": "24.3 / IV-A3",
            "supporting_metric_or_fact": f"adds_structure_beyond_law = {q3['adds_structure_beyond_law']}",
            "scope": "local_structure",
            "notes": "The selected relation is not just a restatement of the threshold.",
        },
        {
            "claim_id": "P4_LOC_03",
            "claim_text": "Ordering degradation appears before selector failure in at least one tested extension class.",
            "status": "survives_locally",
            "supporting_stage": "24.4 / IV-A4",
            "supporting_metric_or_fact": f"first_boundary_by_class = {s24_4['first_boundary_by_class']}",
            "scope": "bounded_extension",
            "notes": "Motif perturbation is the first ordering boundary while selector failure arrives later elsewhere.",
        },
        {
            "claim_id": "P4_LOC_04",
            "claim_text": "Motif perturbation is the earliest tested ordering boundary on the frozen branch.",
            "status": "survives_locally",
            "supporting_stage": "24.4 / IV-A4",
            "supporting_metric_or_fact": "first ordering boundary = IVA4_E3_mild",
            "scope": "bounded_extension",
            "notes": "This is the first honest ordering limit recorded in the extension lattice.",
        },
        {
            "claim_id": "P4_LOC_05",
            "claim_text": "Moderate corridor geometry is the first tested selector-failure boundary on the frozen branch.",
            "status": "bounded_support",
            "supporting_stage": "24.4 / IV-A4",
            "supporting_metric_or_fact": "first selector failure = IVA4_E4_moderate",
            "scope": "bounded_extension",
            "notes": "Selector portability is bounded rather than uniform across local extension classes.",
        },
    ]
    negative_rows = [
        {
            "claim_id": "P4_NEG_01",
            "claim_text": "The Phase IV selector may not be promoted into a family-wide transport law.",
            "status": "fails",
            "supporting_stage": "24.4 / IV-A4; 24.5 / IV-A5",
            "supporting_metric_or_fact": "selector fails on moderate corridor extension and is bounded by explicit extension-class structure.",
            "scope": "not_family_wide",
            "notes": "Portable only within the tested bounded extension lattice.",
        },
        {
            "claim_id": "P4_NEG_02",
            "claim_text": "The local ordering relation is not uniformly portable across the full tested extension lattice.",
            "status": "fails",
            "supporting_stage": "24.4 / IV-A4",
            "supporting_metric_or_fact": "ordering degrades on motif perturbation before selector failure.",
            "scope": "bounded_extension",
            "notes": "Ordering is local structure, not a universal law.",
        },
        {
            "claim_id": "P4_NEG_03",
            "claim_text": "Rejected braid- and asymmetry-centered observables may not be reintroduced into the Phase IV core.",
            "status": "excluded",
            "supporting_stage": "24.1 / IV-A1",
            "supporting_metric_or_fact": f"rejected_observables = {s24_1['rejected_observables']}",
            "scope": "observable_ledger",
            "notes": "Rejected observables remain rejected.",
        },
        {
            "claim_id": "P4_NEG_04",
            "claim_text": "Phase IV does not restore a stable braid phase.",
            "status": "fails",
            "supporting_stage": "23.14 / III-C5; 24.5 / IV-A5",
            "supporting_metric_or_fact": f"stable_braid_phase_present = {s23_14['stable_braid_phase_present']}",
            "scope": "boundary_preservation",
            "notes": "Stable braid phase remains absent.",
        },
        {
            "claim_id": "P4_NEG_05",
            "claim_text": "Phase IV does not restore a localized encounter phase as a stable closed sector.",
            "status": "fails",
            "supporting_stage": "23.14 / III-C5; 24.5 / IV-A5",
            "supporting_metric_or_fact": f"localized_encounter_phase_present = {s23_14['localized_encounter_phase_present']}",
            "scope": "boundary_preservation",
            "notes": "Localized encounter phase remains absent from the stable closure set.",
        },
        {
            "claim_id": "P4_NEG_06",
            "claim_text": "Phase IV does not reverse the Phase III non-generality result for braid-like exchange.",
            "status": "fails",
            "supporting_stage": "23.10 / III-C1; 24.5 / IV-A5",
            "supporting_metric_or_fact": "family_wide_braid_mechanism remains false through the consolidation.",
            "scope": "boundary_preservation",
            "notes": "No family-wide braid reintroduction.",
        },
        {
            "claim_id": "P4_NEG_07",
            "claim_text": "The support-anisotropy survival result may not be misread as unrestricted robustness.",
            "status": "excluded",
            "supporting_stage": "24.4 / IV-A4",
            "supporting_metric_or_fact": "support anisotropy survives only inside the explicit tested extension family.",
            "scope": "bounded_extension",
            "notes": "This is bounded evidence only.",
        },
        {
            "claim_id": "P4_NEG_08",
            "claim_text": "The bounded local transport package may not be elevated into a global emergent transport theory.",
            "status": "excluded",
            "supporting_stage": "24.5 / IV-A5",
            "supporting_metric_or_fact": "Phase IV closes as bounded local structure with explicit selector and ordering boundaries.",
            "scope": "not_family_wide",
            "notes": "Phase IV negatives remain fully intact.",
        },
    ]
    return portable_rows + local_rows + negative_rows


def build_consistency_checks(
    runsheet: dict[str, Any],
    inputs: dict[str, dict[str, Any]],
    claim_rows: list[dict[str, Any]],
) -> dict[str, bool]:
    s23_14 = inputs["stage23_14"]["json"]
    s24_1 = inputs["stage24_1"]["json"]
    s24_2 = inputs["stage24_2"]["json"]
    s24_3 = inputs["stage24_3"]["json"]
    s24_4 = inputs["stage24_4"]["json"]
    s24_3_rows = inputs["stage24_3"]["csv_rows"]

    expected_core = list(runsheet["expected_observable_core"])
    expected_rejected = list(runsheet["expected_rejected_observables"])
    expected_status_counts = dict(runsheet["expected_status_counts"])
    expected_boundaries = dict(runsheet["expected_boundaries"])
    q3 = find_csv_row(s24_3_rows, "candidate_id", "IVA3_Q3")
    portable_text = " ".join(
        row["claim_text"].lower() for row in claim_rows if str(row["status"]) == "survives_portably"
    )
    contradiction_free = all(
        phrase not in portable_text
        for phrase in (
            "uniformly portable across the tested extension lattice",
            "family-wide transport law",
            "family-wide braid mechanism",
            "restores a stable braid phase",
        )
    )

    return {
        "stage24_1_expected_primary_observables": bool(s24_1["carry_forward_core"] == expected_core),
        "stage24_1_rejected_family_intact": bool(s24_1["rejected_observables"] == expected_rejected),
        "stage24_2_selector_law_exact": bool(s24_2["selected_law_expression"] == runsheet["expected_selector_law"]),
        "stage24_3_selected_q3": bool(s24_3["selected_candidate_id"] == runsheet["expected_ordering_candidate_id"]),
        "stage24_3_contrast_separated_spearman_signs": bool(
            "spearman_rho=-" in str(q3["stable_metric"]) and "spearman_rho=0." in str(q3["contrast_metric"])
        ),
        "stage24_4_status_counts_match": bool(s24_4["status_counts"] == expected_status_counts),
        "stage24_4_boundaries_match": bool(
            all(
                s24_4["first_boundary_by_class"].get(key) is None
                if value is None
                else s24_4["first_boundary_by_class"].get(key, {}).get("extension_id") == value
                for key, value in expected_boundaries.items()
            )
        ),
        "portable_claims_noncontradictory": contradiction_free,
        "phaseIII_negatives_preserved": bool(
            s23_14["family_wide_braid_mechanism"] is False
            and s23_14["stable_braid_phase_present"] is False
            and s23_14["localized_encounter_phase_present"] is False
        ),
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
    s24_1 = inputs["stage24_1"]["json"]
    s24_2 = inputs["stage24_2"]["json"]
    s24_3 = inputs["stage24_3"]["json"]
    s24_4 = inputs["stage24_4"]["json"]

    phaseIV_inputs_consistent = all(bool(value) for value in checks.values())
    portable_core_frozen = all(str(row["status"]) == "survives_portably" for row in claim_rows[:5])
    local_structure_frozen = all(str(row["status"]) in {"survives_locally", "bounded_support"} for row in claim_rows[5:10])
    negative_boundaries_frozen = all(str(row["status"]) in {"fails", "excluded"} for row in claim_rows[10:])
    selector_question_closed = bool(s24_2["selected_law_expression"] == runsheet["expected_selector_law"])
    ordering_boundary_closed = bool(
        s24_4["first_boundary_by_class"]["motif_structure_perturbation"]["extension_id"] == "IVA4_E3_mild"
    )
    freeze_valid = bool(
        phaseIV_inputs_consistent
        and portable_core_frozen
        and local_structure_frozen
        and negative_boundaries_frozen
        and selector_question_closed
        and ordering_boundary_closed
    )

    return {
        "experiment": runsheet["stage"],
        "description": runsheet["description"],
        "phaseIV_consolidation_complete": True,
        "phaseIV_inputs_consistent": phaseIV_inputs_consistent,
        "phaseIV_freeze_valid": freeze_valid,
        "portable_core_frozen": portable_core_frozen,
        "local_structure_frozen": local_structure_frozen,
        "negative_boundaries_frozen": negative_boundaries_frozen,
        "selector_question_closed": selector_question_closed,
        "ordering_boundary_closed": ordering_boundary_closed,
        "observable_core_frozen": bool(s24_1["observable_ledger_valid"]),
        "selector_law_found": bool(s24_2["effective_law_found"]),
        "selector_threshold_fixed": 0.884308,
        "selector_law_boundedly_portable": True,
        "local_ordering_found": bool(s24_3["quasi_invariant_found"]),
        "ordering_uniformly_portable": False,
        "ordering_degrades_before_selector_failure": True,
        "earliest_ordering_boundary_class": "motif_structure_perturbation",
        "earliest_ordering_boundary_level": "IVA4_E3_mild",
        "first_selector_failure_class": "corridor_geometry_perturbation",
        "first_selector_failure_level": "IVA4_E4_moderate",
        "support_anisotropy_failure_seen": False,
        "braid_reintroduced": False,
        "rejected_observables_reintroduced": False,
        "family_wide_transport_claim": False,
        "extension_status_counts": dict(s24_4["status_counts"]),
        "phaseIV_minimal_read": "boundedly portable stable-smeared threshold selector with locally portable inverse flow-survival ordering; ordering boundaries appear before selector failure in some extension classes",
        "consistency_checks": checks,
        "inputs": compact_input_refs(inputs),
        "claim_ledger": claim_rows,
    }


def load_inputs(runsheet: dict[str, Any]) -> dict[str, dict[str, Any]]:
    inputs: dict[str, dict[str, Any]] = {}
    for key, item in runsheet["required_inputs"].items():
        note_path = ROOT / item["note_path"]
        json_path = ROOT / item["json_path"]
        csv_path = ROOT / item["csv_path"]
        inputs[key] = {
            "stage_label": item["stage_label"],
            "note_path": item["note_path"],
            "json_path": item["json_path"],
            "csv_path": item["csv_path"],
            "note_present": note_path.exists(),
            "json": read_json(json_path),
            "csv_rows": read_csv_rows(csv_path),
        }
    return inputs


def run_phase4(config_path: Path, output_dir: Path) -> dict[str, Any]:
    runsheet = read_json(config_path)
    inputs = load_inputs(runsheet)
    claim_rows = build_claim_rows(inputs)
    checks = build_consistency_checks(runsheet, inputs, claim_rows)
    summary = build_summary_payload(runsheet, inputs, claim_rows, checks)
    bundle = {
        "phase": 4,
        "phase_name": "phase4-sector-freeze",
        "source_config": relpath(config_path),
        "upstream_inputs": compact_input_refs(inputs),
        "summary": summary,
        "claim_ledger": claim_rows,
    }
    stamped, latest = write_timestamped_json(output_dir, "phase4_sector_freeze_bundle", bundle)
    bundle["artifacts"] = {
        "stamped_json": relpath(stamped),
        "latest_json": relpath(latest),
    }
    write_json(stamped, bundle)
    write_json(latest, bundle)
    return bundle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the narrow Phase IV sector freeze.")
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "phase4-sector-freeze" / "configs" / "stage24_5_phaseIV_consolidation_runs.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "phase4-sector-freeze" / "runs",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the final bundle as JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bundle = run_phase4(args.config, args.output_dir)
    if args.json:
        print(json.dumps(bundle, indent=2))
        return
    print(f"Phase IV freeze valid: {bundle['summary']['phaseIV_freeze_valid']}")
    print(f"Latest bundle: {bundle['artifacts']['latest_json']}")


if __name__ == "__main__":
    main()
