#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = Path(__file__).resolve().parent / "results"

METRIC_FIELD = REPO_ROOT / "data" / "scalar_kernel_graph_metric_field_latest.json"
RECOVERABILITY_GRADIENT = REPO_ROOT / "data" / "scalar_kernel_graph_recoverability_gradient_latest.json"
RECOVERABILITY_GRADIENT_SHELL_NATIVE = (
    REPO_ROOT / "data" / "scalar_kernel_graph_recoverability_gradient_shell_native_latest.json"
)
CURRENT_CLOSURE = REPO_ROOT / "data" / "scalar_kernel_graph_current_closure_shell_native_latest.json"
INHOMOGENEITY_CLOSURE = REPO_ROOT / "data" / "scalar_kernel_graph_current_closure_inhomogeneity_latest.json"
DISORDER_NATIVE_FLUX = REPO_ROOT / "data" / "scalar_kernel_graph_current_closure_radial_disorder_native_flux_latest.json"

CSV_PATH = RESULTS_DIR / "physics_observables.csv"
JSON_PATH = RESULTS_DIR / "physics_observables.json"
SUMMARY_PATH = RESULTS_DIR / "physics_observables_summary.md"

FIELDNAMES = [
    "observable_id",
    "physics_analogue",
    "haos_source",
    "proxy_quantity",
    "value",
    "target",
    "status",
    "scope_note",
]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing frozen scalar artifact: {path.relative_to(REPO_ROOT)}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def artifact(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def fmt(value: float) -> str:
    return f"{value:.6f}"


def all_scalar_cases(data: dict[str, Any]) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for key in ("disorder_cases", "kernel_cases", "cases"):
        cases.extend(data.get(key, []))
    return cases


def row(
    observable_id: str,
    physics_analogue: str,
    haos_source: str,
    proxy_quantity: str,
    value: str,
    target: str,
    status: str,
    scope_note: str,
) -> dict[str, str]:
    return {
        "observable_id": observable_id,
        "physics_analogue": physics_analogue,
        "haos_source": haos_source,
        "proxy_quantity": proxy_quantity,
        "value": value,
        "target": target,
        "status": status,
        "scope_note": scope_note,
    }


def leq_status(value: float, threshold: float) -> str:
    return "PASS" if value <= threshold else "OPEN"


def geq_status(value: float, threshold: float) -> str:
    return "PASS" if value >= threshold else "OPEN"


def build_rows() -> tuple[list[dict[str, str]], dict[str, Any]]:
    metric = load_json(METRIC_FIELD)
    gradient = load_json(RECOVERABILITY_GRADIENT)
    shell_native_gradient = load_json(RECOVERABILITY_GRADIENT_SHELL_NATIVE)
    current = load_json(CURRENT_CLOSURE)
    inhomogeneity = load_json(INHOMOGENEITY_CLOSURE)
    disorder_flux = load_json(DISORDER_NATIVE_FLUX)

    rows: list[dict[str, str]] = []

    metric_thresholds = metric["config"]["thresholds"]
    metric_cases = all_scalar_cases(metric)
    min_metric_eigenvalue = min(
        float(case["metric_field"]["coarse"]["min_eigenvalue"]) for case in metric_cases
    )
    max_metric_anisotropy = max(
        float(case["metric_field"]["coarse"]["mean_anisotropy"]) for case in metric_cases
    )
    max_metric_drift = float(metric["max_refinement_drift"])
    rows.append(
        row(
            "metric_positive_bulk",
            "effective metric positivity proxy",
            artifact(METRIC_FIELD),
            "minimum coarse bulk metric eigenvalue",
            fmt(min_metric_eigenvalue),
            f">= {fmt(float(metric_thresholds['min_min_eigenvalue']))}",
            geq_status(min_metric_eigenvalue, float(metric_thresholds["min_min_eigenvalue"])),
            "Supports a bounded local metric-like tensor read, not a spacetime metric.",
        )
    )
    rows.append(
        row(
            "metric_anisotropy_bound",
            "local isotropy proxy",
            artifact(METRIC_FIELD),
            "maximum coarse mean anisotropy",
            fmt(max_metric_anisotropy),
            f"<= {fmt(float(metric_thresholds['max_mean_anisotropy']))}",
            leq_status(max_metric_anisotropy, float(metric_thresholds["max_mean_anisotropy"])),
            "Checks whether the tensor read stays close to the scalar-carrier isotropic baseline.",
        )
    )
    rows.append(
        row(
            "metric_refinement_stability",
            "continuum-limit sanity proxy",
            artifact(METRIC_FIELD),
            "maximum normalized metric refinement drift",
            fmt(max_metric_drift),
            f"<= {fmt(float(metric_thresholds['max_refinement_drift']))}",
            leq_status(max_metric_drift, float(metric_thresholds["max_refinement_drift"])),
            "A scaling sanity check only; it does not establish a continuum limit.",
        )
    )

    gradient_thresholds = gradient["config"]["thresholds"]
    gradient_cases = all_scalar_cases(gradient)
    min_radial_alignment = min(float(case["radial_alignment"]) for case in gradient_cases)
    min_power_r2 = min(float(case["profile"]["power_fit_r2"]) for case in gradient_cases)
    max_flux_cv = max(float(case["profile"]["flux_constancy_cv"]) for case in gradient_cases)
    gradient_drift = float(gradient["max_profile_drift"])
    rows.append(
        row(
            "recoverability_radial_alignment",
            "radial response-direction proxy",
            artifact(RECOVERABILITY_GRADIENT),
            "minimum radial alignment",
            fmt(min_radial_alignment),
            f">= {fmt(float(gradient_thresholds['min_radial_alignment']))}",
            geq_status(min_radial_alignment, float(gradient_thresholds["min_radial_alignment"])),
            "Directionality is strong, but direction alone is not a force law.",
        )
    )
    rows.append(
        row(
            "recoverability_inverse_square_residual",
            "raw inverse-square response-law proxy",
            artifact(RECOVERABILITY_GRADIENT),
            "worst flux constancy CV / profile drift",
            f"{fmt(max_flux_cv)} / {fmt(gradient_drift)}",
            f"<= {fmt(float(gradient_thresholds['max_flux_constancy_cv']))} / <= {fmt(float(gradient_thresholds['max_refinement_profile_drift']))}",
            "PASS"
            if max_flux_cv <= float(gradient_thresholds["max_flux_constancy_cv"])
            and gradient_drift <= float(gradient_thresholds["max_refinement_profile_drift"])
            else "OPEN",
            "Keeps the 52.1 boundary explicit: radial carrier response is present, inverse-square closure is not.",
        )
    )
    rows.append(
        row(
            "recoverability_power_fit_floor",
            "power-law scaling proxy",
            artifact(RECOVERABILITY_GRADIENT),
            "minimum power-fit R^2",
            fmt(min_power_r2),
            f">= {fmt(float(gradient_thresholds['min_power_fit_r2']))}",
            geq_status(min_power_r2, float(gradient_thresholds["min_power_fit_r2"])),
            "A deliberately hard guardrail against reading noisy response curves as laws.",
        )
    )

    shell_native_thresholds = shell_native_gradient["config"]["thresholds"]
    shell_native_cases = all_scalar_cases(shell_native_gradient)
    shell_native_min_green_r2 = min(float(case["green_fit_r2"]) for case in shell_native_cases)
    shell_native_max_power_deviation = max(
        abs(float(case["profile"]["power_slope"]) + 2.0) for case in shell_native_cases
    )
    shell_native_min_power_r2 = min(float(case["profile"]["power_fit_r2"]) for case in shell_native_cases)
    shell_native_max_flux_cv = max(float(case["profile"]["flux_constancy_cv"]) for case in shell_native_cases)
    shell_native_drift = float(shell_native_gradient["max_profile_drift"])
    shell_native_pass = all(
        [
            all(bool(case["pass"]) for case in shell_native_cases),
            shell_native_min_green_r2 >= float(shell_native_thresholds["min_green_fit_r2"]),
            shell_native_max_power_deviation <= float(shell_native_thresholds["max_power_deviation"]),
            shell_native_min_power_r2 >= float(shell_native_thresholds["min_power_fit_r2"]),
            shell_native_max_flux_cv <= float(shell_native_thresholds["max_flux_constancy_cv"]),
            shell_native_drift <= float(shell_native_thresholds["max_refinement_profile_drift"]),
        ]
    )
    rows.append(
        row(
            "shell_native_inverse_square_closure",
            "shell-native Newtonian-like response proxy",
            artifact(RECOVERABILITY_GRADIENT_SHELL_NATIVE),
            "min Green R2 / max |slope+2| / max scaled-response CV / profile drift",
            f"{fmt(shell_native_min_green_r2)} / {fmt(shell_native_max_power_deviation)} / {fmt(shell_native_max_flux_cv)} / {fmt(shell_native_drift)}",
            f">= {fmt(float(shell_native_thresholds['min_green_fit_r2']))} / <= {fmt(float(shell_native_thresholds['max_power_deviation']))} / <= {fmt(float(shell_native_thresholds['max_flux_constancy_cv']))} / <= {fmt(float(shell_native_thresholds['max_refinement_profile_drift']))}",
            "PASS" if shell_native_pass else "OPEN",
            "Closes the inverse-square proxy only after shell-native law-aware reconstruction; this is not a Newtonian gravity claim.",
        )
    )

    current_thresholds = current["config"]["thresholds"]
    current_cases = current["cases"]
    current_max_median_error = max(float(case["median_relative_error"]) for case in current_cases)
    current_max_p90_error = max(float(case["p90_relative_error"]) for case in current_cases)
    current_max_diffusivity_gap = max(float(case["diffusivity_match_gap"]) for case in current_cases)
    current_max_kappa_cv = max(float(case["shell_kappa_cv"]) for case in current_cases)
    current_max_drift = float(current["max_profile_drift"])
    current_pass = all(
        [
            current_max_median_error <= float(current_thresholds["max_median_relative_error"]),
            current_max_p90_error <= float(current_thresholds["max_p90_relative_error"]),
            current_max_diffusivity_gap <= float(current_thresholds["max_diffusivity_match_gap"]),
            current_max_kappa_cv <= float(current_thresholds["max_shell_kappa_cv"]),
            current_max_drift <= float(current_thresholds["max_refinement_profile_drift"]),
        ]
    )
    rows.append(
        row(
            "clean_current_closure",
            "conserved-current analogue proxy",
            artifact(CURRENT_CLOSURE),
            "max median error / max diffusivity gap / max profile drift",
            f"{fmt(current_max_median_error)} / {fmt(current_max_diffusivity_gap)} / {fmt(current_max_drift)}",
            f"<= {fmt(float(current_thresholds['max_median_relative_error']))} / <= {fmt(float(current_thresholds['max_diffusivity_match_gap']))} / <= {fmt(float(current_thresholds['max_refinement_profile_drift']))}",
            "PASS" if current_pass else "OPEN",
            "Supports bounded clean-line constitutive closure, not a conservation theorem.",
        )
    )
    rows.append(
        row(
            "clean_current_tail_error",
            "transport residual proxy",
            artifact(CURRENT_CLOSURE),
            "maximum p90 relative current error / shell-kappa CV",
            f"{fmt(current_max_p90_error)} / {fmt(current_max_kappa_cv)}",
            f"<= {fmt(float(current_thresholds['max_p90_relative_error']))} / <= {fmt(float(current_thresholds['max_shell_kappa_cv']))}",
            "PASS"
            if current_max_p90_error <= float(current_thresholds["max_p90_relative_error"])
            and current_max_kappa_cv <= float(current_thresholds["max_shell_kappa_cv"])
            else "OPEN",
            "Retains residual accounting rather than converting closure into exact equality.",
        )
    )

    inhomogeneity_thresholds = inhomogeneity["config"]["thresholds"]
    family_summaries = {item["family"]: item for item in inhomogeneity["family_summaries"]}
    radial_summary = family_summaries["radial"]
    bump_summary = family_summaries["bump"]
    radial_cases = [case for case in inhomogeneity["cases"] if case["family"] == "radial"]
    bump_cases = [case for case in inhomogeneity["cases"] if case["family"] == "bump"]
    radial_min_metric_corr = min(float(case["metric_tracking_abs_corr"]) for case in radial_cases)
    bump_min_metric_corr = min(float(case["metric_tracking_abs_corr"]) for case in bump_cases)
    rows.append(
        row(
            "smooth_inhomogeneity_metric_transport_tracking",
            "smooth weak-field co-deformation proxy",
            artifact(INHOMOGENEITY_CLOSURE),
            "radial max drift / minimum metric-tracking |corr|",
            f"{fmt(float(radial_summary['max_refinement_profile_drift']))} / {fmt(radial_min_metric_corr)}",
            f"<= {fmt(float(inhomogeneity_thresholds['max_refinement_profile_drift']))} / >= {fmt(float(inhomogeneity_thresholds['min_metric_tracking_abs_corr']))}",
            "PASS" if radial_summary["all_cases_pass"] and radial_summary["drift_pass"] else "OPEN",
            "The smooth branch co-deforms; this is still a scalar-carrier transport proxy, not a weak-field GR claim.",
        )
    )
    rows.append(
        row(
            "localized_bump_boundary",
            "localized matter-like excitation boundary",
            artifact(INHOMOGENEITY_CLOSURE),
            "bump max drift / minimum metric-tracking |corr|",
            f"{fmt(float(bump_summary['max_refinement_profile_drift']))} / {fmt(bump_min_metric_corr)}",
            f"<= {fmt(float(inhomogeneity_thresholds['max_refinement_profile_drift']))} / >= {fmt(float(inhomogeneity_thresholds['min_metric_tracking_abs_corr']))}",
            "PASS" if bump_summary["all_cases_pass"] and bump_summary["drift_pass"] else "OPEN",
            "This is the useful failure boundary: localized bump closure is not yet earned.",
        )
    )

    disorder_thresholds = disorder_flux["config"]["thresholds"]
    disorder_cases = disorder_flux["cases"]
    disorder_all_cases_pass = all(bool(case["pass"]) for case in disorder_cases)
    disorder_drift = float(disorder_flux["max_refinement_profile_drift"])
    disorder_max_median = max(float(case["median_relative_error_median"]) for case in disorder_cases)
    disorder_max_p90 = max(float(case["p90_relative_error_median"]) for case in disorder_cases)
    min_seed_pass = min(int(case["pass_count"]) for case in disorder_cases)
    max_seed_count = max(int(case["trial_count"]) for case in disorder_cases)
    rows.append(
        row(
            "disorder_native_aggregate_flux",
            "disorder-robust transport-law proxy",
            artifact(DISORDER_NATIVE_FLUX),
            "max aggregated median error / max aggregated p90 error / max drift",
            f"{fmt(disorder_max_median)} / {fmt(disorder_max_p90)} / {fmt(disorder_drift)}",
            f"<= {fmt(float(disorder_thresholds['max_median_relative_error']))} / <= {fmt(float(disorder_thresholds['max_p90_relative_error']))} / <= {fmt(float(disorder_thresholds['max_refinement_profile_drift']))}",
            "PASS"
            if disorder_all_cases_pass
            and disorder_drift <= float(disorder_thresholds["max_refinement_profile_drift"])
            else "OPEN",
            "Aggregated mild-disorder transport closes after native flux reconstruction.",
        )
    )
    rows.append(
        row(
            "disorder_native_seed_margin",
            "seed-universal disorder margin proxy",
            artifact(DISORDER_NATIVE_FLUX),
            "minimum seed-level pass count",
            f"{min_seed_pass}/{max_seed_count}",
            f"{max_seed_count}/{max_seed_count} for seed-universal closure",
            "PASS" if min_seed_pass == max_seed_count else "WATCH",
            "Aggregates pass, but one tested case still has a marginal seed-level failure.",
        )
    )

    metadata = {
        "program": "HAOS-IIP physics bridge 53.0",
        "claim_type": "post-processing proxy dictionary over frozen scalar artifacts",
        "non_claims": [
            "no modification of HAOS core",
            "no continuum ontology",
            "no Einstein-Hilbert recovery claim",
            "no conserved stress-energy theorem",
            "no arbitrary localized inhomogeneity closure",
        ],
        "source_artifacts": [
            artifact(METRIC_FIELD),
            artifact(RECOVERABILITY_GRADIENT),
            artifact(RECOVERABILITY_GRADIENT_SHELL_NATIVE),
            artifact(CURRENT_CLOSURE),
            artifact(INHOMOGENEITY_CLOSURE),
            artifact(DISORDER_NATIVE_FLUX),
        ],
        "status_counts": {
            "PASS": sum(1 for item in rows if item["status"] == "PASS"),
            "OPEN": sum(1 for item in rows if item["status"] == "OPEN"),
            "WATCH": sum(1 for item in rows if item["status"] == "WATCH"),
        },
    }
    return rows, metadata


def write_csv(rows: list[dict[str, str]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(rows: list[dict[str, str]], metadata: dict[str, Any]) -> None:
    payload = {
        "metadata": metadata,
        "observables": rows,
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_summary(rows: list[dict[str, str]], metadata: dict[str, Any]) -> None:
    table_rows = "\n".join(
        "| {observable_id} | {physics_analogue} | {value} | {target} | {status} |".format(**item)
        for item in rows
    )
    open_rows = [item for item in rows if item["status"] != "PASS"]
    open_list = "\n".join(
        f"- `{item['observable_id']}`: {item['scope_note']}" for item in open_rows
    )
    if not open_list:
        open_list = "- none"

    summary = f"""# HAOS-IIP Physics Bridge Observables 53.0

## Scope

This file is generated by `experiments/physics_bridge/physics_observables.py`.

It is a post-processing bridge over frozen scalar-carrier artifacts. It does not modify HAOS core, re-run the scalar simulations, or promote the bridge quantities into physical ontology.

## Bridge Table

| Observable | Physics-facing analogue | Value | Target | Status |
| --- | --- | --- | --- | --- |
{table_rows}

## Open Or Watched Boundaries

{open_list}

## Status Counts

- `PASS`: {metadata['status_counts']['PASS']}
- `OPEN`: {metadata['status_counts']['OPEN']}
- `WATCH`: {metadata['status_counts']['WATCH']}

## Authority

- CSV: `experiments/physics_bridge/results/physics_observables.csv`
- JSON: `experiments/physics_bridge/results/physics_observables.json`
- Generator: `experiments/physics_bridge/physics_observables.py`
"""
    SUMMARY_PATH.write_text(summary, encoding="utf-8")


def main() -> None:
    rows, metadata = build_rows()
    write_csv(rows)
    write_json(rows, metadata)
    write_summary(rows, metadata)
    print("Physics bridge observables generated.")
    print(f"CSV: {CSV_PATH.relative_to(REPO_ROOT)}")
    print(f"JSON: {JSON_PATH.relative_to(REPO_ROOT)}")
    print(f"Summary: {SUMMARY_PATH.relative_to(REPO_ROOT)}")
    print(
        "Statement: 53.0 is a physics-facing proxy dictionary over frozen scalar artifacts; "
        "it introduces no HAOS core changes and no continuum-ontology claim."
    )


if __name__ == "__main__":
    main()
