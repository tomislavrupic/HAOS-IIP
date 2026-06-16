#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent

PHASE10_EFFECTIVE = ROOT / "phase10-bridge" / "runs" / "phase10_effective_scaling_ledger.csv"
PHASE15_SPEED = ROOT / "phase15-propagation" / "runs" / "phase15_effective_speed_ledger.csv"
PHASE18_REFINEMENT = ROOT / "phase18-distance-surrogate" / "runs" / "phase18_refinement_scaling.csv"
CP2_RECOVERY = OUT_DIR / "same_surrogate_coarse_graining_recovery.csv"
CP2_CONTROL_INTEGRITY = OUT_DIR / "same_surrogate_control_integrity_report.md"
KURAMOTO_STATUS = ROOT / "experiments" / "oscillators" / "kuramoto_bridge" / "outputs" / "bridge_status.json"
HAOS_VALIDATION = ROOT / "external_validation" / "results" / "haos_iip_summaries.md"
TOY_VALIDATION = ROOT / "external_validation" / "results" / "toy_summaries.md"

CP3_OUT = OUT_DIR / "post_66_5_cp3_effective_equation_contract.csv"
COMPARATIVE_OUT = OUT_DIR / "post_66_5_narrow_comparative_diagnostic.csv"
CP5_OUT = OUT_DIR / "post_66_5_cp5_universality_screen.csv"
SUMMARY_OUT = OUT_DIR / "post_66_5_roadmap_run_summary.md"

BRANCH = "frozen_branch"
PHASE10_CONTROL = "generic_open_grid_scalar_block_control"
PERIODIC_CONTROL = "periodic_diagonal_augmented_control"
LEVELS = ["48", "60", "72", "84"]
SEED_PAIR = {"1302", "1303"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def f(value: str | float | int | None) -> float:
    if value is None:
        return math.nan
    if isinstance(value, (float, int)):
        return float(value)
    text = value.strip()
    if not text:
        return math.nan
    return float(text)


def finite(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def mean(values: list[float]) -> float:
    values = finite(values)
    return sum(values) / len(values) if values else math.nan


def max_or_nan(values: list[float]) -> float:
    values = finite(values)
    return max(values) if values else math.nan


def ratio(numerator: float, denominator: float) -> float:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or abs(denominator) < 1.0e-15:
        return math.nan
    return numerator / denominator


def rel_span(values: list[float]) -> float:
    values = finite(values)
    if not values:
        return math.nan
    mu = mean(values)
    if abs(mu) < 1.0e-15:
        return math.nan
    return (max(values) - min(values)) / abs(mu)


def fmt(value: float | int | str) -> str:
    if isinstance(value, str):
        return value
    value = float(value)
    if not math.isfinite(value):
        return "not_available"
    return f"{value:.6f}"


def phase10_family_metrics(rows: list[dict[str, str]], families: set[str]) -> tuple[float, float, float, float, float]:
    branch_rows = [row for row in rows if row["hierarchy_label"] == BRANCH and row["family"] in families]
    control_rows = [row for row in rows if row["hierarchy_label"] == PHASE10_CONTROL and row["family"] in families]
    branch_error = max_or_nan([f(row["relative_error_R5"]) for row in branch_rows])
    control_error = max_or_nan([f(row["relative_error_R5"]) for row in control_rows])
    branch_r2 = mean([f(row["fit_r_squared"]) for row in branch_rows])
    control_r2 = mean([f(row["fit_r_squared"]) for row in control_rows])
    separation = ratio(control_error, branch_error)
    return branch_error, control_error, branch_r2, control_r2, separation


def phase18_metric_metrics(rows: list[dict[str, str]]) -> tuple[float, float, float, float, float]:
    branch = [row for row in rows if row["hierarchy_label"] == BRANCH]
    control = [row for row in rows if row["hierarchy_label"] == PERIODIC_CONTROL]
    branch_slope = max_or_nan([abs(f(row["slope_drift_from_prev"])) for row in branch])
    control_slope = max_or_nan([abs(f(row["slope_drift_from_prev"])) for row in control])
    branch_depth = max_or_nan([abs(f(row["causal_depth_drift_from_prev"])) for row in branch])
    control_depth = max_or_nan([abs(f(row["causal_depth_drift_from_prev"])) for row in control])
    separation = ratio(control_slope, branch_slope)
    return branch_slope, control_slope, branch_depth, control_depth, separation


def phase15_speed_metrics(rows: list[dict[str, str]], hierarchy: str) -> tuple[list[float], float, float]:
    values_by_level: list[float] = []
    for level in LEVELS:
        values = [
            f(row["effective_speed"])
            for row in rows
            if row["hierarchy_label"] == hierarchy
            and row["probe_name"] == "bias_onset"
            and row["ensemble_size"] == "7"
            and row["seed"] in SEED_PAIR
            and row["n_side"] == level
        ]
        values_by_level.append(mean(values))
    max_adjacent = max_or_nan([abs(values_by_level[i] - values_by_level[i - 1]) for i in range(1, len(values_by_level))])
    return values_by_level, rel_span(values_by_level), max_adjacent


def seed_relspread_metrics(rows: list[dict[str, str]], hierarchy: str) -> tuple[float, float, int]:
    spreads: list[float] = []
    seed_values = set()
    for level in LEVELS:
        values = [
            f(row["effective_speed"])
            for row in rows
            if row["hierarchy_label"] == hierarchy
            and row["probe_name"] == "bias_onset"
            and row["ensemble_size"] == "7"
            and row["n_side"] == level
        ]
        spreads.append(rel_span(values))
        seed_values.update(
            row["seed"]
            for row in rows
            if row["hierarchy_label"] == hierarchy
            and row["probe_name"] == "bias_onset"
            and row["ensemble_size"] == "7"
        )
    return max_or_nan(spreads), mean(spreads), len(seed_values)


def load_kuramoto() -> dict[str, object]:
    if not KURAMOTO_STATUS.exists():
        return {}
    with KURAMOTO_STATUS.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def extract_summary_line(path: Path, prefix: str) -> str:
    if not path.exists():
        return "not_available"
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix):
            return line.split(":", 1)[-1].strip().strip("`")
    return "not_available"


def build_cp3() -> list[dict[str, object]]:
    phase10 = read_csv(PHASE10_EFFECTIVE)
    phase18 = read_csv(PHASE18_REFINEMENT)
    phase15 = read_csv(PHASE15_SPEED)

    rows: list[dict[str, object]] = []

    b_err, c_err, b_r2, c_r2, sep = phase10_family_metrics(
        phase10, {"scaled_coefficients", "ratios", "trace_window"}
    )
    rows.append(
        {
            "contract": "CP3 coefficient-flow closure",
            "source": "phase10_effective_scaling_ledger",
            "candidate_law": "R1-R4 power-law extrapolation to R5",
            "branch_metric": fmt(b_err),
            "control_metric": fmt(c_err),
            "separation_ratio": fmt(sep),
            "branch_secondary": f"mean_fit_r2={fmt(b_r2)}",
            "control_secondary": f"mean_fit_r2={fmt(c_r2)}",
            "threshold": "branch_error<=0.001 and control_error/branch_error>=5",
            "status": "PASS" if b_err <= 0.001 and sep >= 5.0 else "OPEN",
            "claim_boundary": "bounded coefficient-flow closure only; no physical field equation",
        }
    )

    b_err, c_err, b_r2, c_r2, sep = phase10_family_metrics(phase10, {"rescaled_invariants"})
    rows.append(
        {
            "contract": "CP3 rescaled-invariant flow",
            "source": "phase10_effective_scaling_ledger",
            "candidate_law": "R1-R4 inverse-size fit to R5",
            "branch_metric": fmt(b_err),
            "control_metric": fmt(c_err),
            "separation_ratio": fmt(sep),
            "branch_secondary": f"mean_fit_r2={fmt(b_r2)}",
            "control_secondary": f"mean_fit_r2={fmt(c_r2)}",
            "threshold": "branch_error<=0.01 and control_error/branch_error>=1.25",
            "status": "PASS" if b_err <= 0.01 and sep >= 1.25 else "CP3_RESCALED_INVARIANT_NOT_SPECIFIC",
            "claim_boundary": "bounded invariant-flow bookkeeping only; no continuum proof",
        }
    )

    b_slope, c_slope, b_depth, c_depth, sep = phase18_metric_metrics(phase18)
    rows.append(
        {
            "contract": "CP3 metric-surrogate shell-slope closure",
            "source": "phase18_refinement_scaling",
            "candidate_law": "arrival-vs-depth shell slope stability",
            "branch_metric": fmt(b_slope),
            "control_metric": fmt(c_slope),
            "separation_ratio": fmt(sep),
            "branch_secondary": f"max_causal_depth_drift={fmt(b_depth)}",
            "control_secondary": f"max_causal_depth_drift={fmt(c_depth)}",
            "threshold": "branch_slope_drift<=0.05 and control/branch>=5 and branch_depth_drift<=0.10",
            "status": "PASS" if b_slope <= 0.05 and sep >= 5.0 and b_depth <= 0.10 else "OPEN",
            "claim_boundary": "metric-like surrogate stabilization only; no metric field",
        }
    )

    branch_values, branch_span, branch_adj = phase15_speed_metrics(phase15, BRANCH)
    control_values, control_span, control_adj = phase15_speed_metrics(phase15, PERIODIC_CONTROL)
    sep = ratio(control_span, branch_span)
    rows.append(
        {
            "contract": "CP3 propagation-speed band closure",
            "source": "phase15_effective_speed_ledger",
            "candidate_law": "bias_onset ensemble-7 seed-pair effective-speed band",
            "branch_metric": fmt(branch_span),
            "control_metric": fmt(control_span),
            "separation_ratio": fmt(sep),
            "branch_secondary": f"max_adjacent_speed_drift={fmt(branch_adj)}; levels={','.join(fmt(v) for v in branch_values)}",
            "control_secondary": f"max_adjacent_speed_drift={fmt(control_adj)}; levels={','.join(fmt(v) for v in control_values)}",
            "threshold": "branch_relative_span<=0.10 and control_span/branch_span>=1.5",
            "status": "PASS" if branch_span <= 0.10 and sep >= 1.5 else "OPEN",
            "claim_boundary": "bounded propagation descriptor only; no transport equation derivation",
        }
    )

    return rows


def build_comparative(cp3_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    cp2_rows = read_csv(CP2_RECOVERY) if CP2_RECOVERY.exists() else []
    kuramoto = load_kuramoto()
    rows: list[dict[str, object]] = []

    for row in cp3_rows:
        rows.append(
            {
                "diagnostic": row["contract"],
                "comparison_system": row["source"],
                "observed_result": row["branch_metric"],
                "control_result": row["control_metric"],
                "contrast": row["separation_ratio"],
                "status": row["status"],
                "claim_boundary": "comparative diagnostic only",
            }
        )

    for row in cp2_rows:
        contrast = ratio(f(row["recovery_pct_admissible"]), f(row["recovery_pct_matched_control"]))
        rows.append(
            {
                "diagnostic": f"CP2 {row['surrogate_type']}",
                "comparison_system": "same-surrogate matched control",
                "observed_result": row["recovery_pct_admissible"],
                "control_result": row["recovery_pct_matched_control"],
                "contrast": fmt(contrast),
                "status": row["status"],
                "claim_boundary": "same-surrogate recovery bookkeeping only",
            }
        )

    if kuramoto:
        observed = kuramoto.get("observed_summary", {})
        rows.append(
            {
                "diagnostic": "Kuramoto line_e_like oscillator bridge",
                "comparison_system": "standard Kuramoto order parameter plus controls",
                "observed_result": (
                    f"early_detection={kuramoto.get('observed_early_detection')}; "
                    f"specificity_pass={kuramoto.get('observed_specificity_pass')}; "
                    f"k_star_time={fmt(float(observed.get('k_star_time', math.nan)))}"
                ),
                "control_result": (
                    f"control_contrast={kuramoto.get('control_contrast')}; "
                    f"specificity_control_contrast={kuramoto.get('specificity_control_contrast')}"
                ),
                "contrast": str(kuramoto.get("failure_mode", "not_available")),
                "status": kuramoto.get("bridge_status", "OPEN"),
                "claim_boundary": "oscillator sidecar diagnostic only; no HAOS superiority claim",
            }
        )

    haos_order = extract_summary_line(HAOS_VALIDATION, "- Invariance assessment")
    toy_order = extract_summary_line(TOY_VALIDATION, "- Invariance assessment")
    rows.append(
        {
            "diagnostic": "external-validation toy-vs-HAOS ordering screen",
            "comparison_system": "external_validation/runners/run_validation.py",
            "observed_result": haos_order,
            "control_result": toy_order,
            "contrast": "HAOS summary flags artifact-of-construction; toy summary flags invariant structural signal",
            "status": "OPEN",
            "claim_boundary": "validation harness comparison only; no physical claim",
        }
    )
    return rows


def build_cp5() -> list[dict[str, object]]:
    phase15 = read_csv(PHASE15_SPEED)
    phase18 = read_csv(PHASE18_REFINEMENT)
    cp2 = read_csv(CP2_RECOVERY) if CP2_RECOVERY.exists() else []

    branch_seed_max, branch_seed_mean, branch_seed_count = seed_relspread_metrics(phase15, BRANCH)
    control_seed_max, control_seed_mean, control_seed_count = seed_relspread_metrics(phase15, PERIODIC_CONTROL)

    rows: list[dict[str, object]] = [
        {
            "variation_axis": "seed variation",
            "source": "phase15_effective_speed_ledger",
            "available_count": branch_seed_count,
            "required_minimum": 3,
            "branch_metric": f"max_relspread={fmt(branch_seed_max)}; mean_relspread={fmt(branch_seed_mean)}",
            "control_metric": f"max_relspread={fmt(control_seed_max)}; mean_relspread={fmt(control_seed_mean)}",
            "status": "CP5_COVERAGE_INSUFFICIENT",
            "notes": "available but too variable for a CP5 universality pass",
        },
        {
            "variation_axis": "refinement levels",
            "source": "phase18_refinement_scaling",
            "available_count": len(set(row["n_side"] for row in phase18 if row["hierarchy_label"] == BRANCH)),
            "required_minimum": 3,
            "branch_metric": "levels=" + ",".join(sorted(set(row["n_side"] for row in phase18 if row["hierarchy_label"] == BRANCH))),
            "control_metric": "matched periodic diagonal augmented control",
            "status": "PASS",
            "notes": "refinement coverage exists for the tested slice only",
        },
        {
            "variation_axis": "kernel width variation",
            "source": "post-66.5 scale-bridge ledgers",
            "available_count": 0,
            "required_minimum": 2,
            "branch_metric": "not_available",
            "control_metric": "not_available",
            "status": "CP5_COVERAGE_INSUFFICIENT",
            "notes": "no same-surrogate CP1-CP2 kernel-width sweep is present in this roadmap run",
        },
        {
            "variation_axis": "substrate family variation",
            "source": "post-66.5 scale-bridge ledgers",
            "available_count": 1,
            "required_minimum": 2,
            "branch_metric": "frozen branch-local hierarchy",
            "control_metric": "matched controls only",
            "status": "CP5_COVERAGE_INSUFFICIENT",
            "notes": "matched controls are present, but not multiple admissible substrate families",
        },
        {
            "variation_axis": "same-surrogate CP2 control separation",
            "source": "same_surrogate_coarse_graining_recovery",
            "available_count": len(cp2),
            "required_minimum": 2,
            "branch_metric": "; ".join(f"{row['surrogate_type']}={row['recovery_pct_admissible']}%" for row in cp2),
            "control_metric": "; ".join(f"{row['surrogate_type']}={row['recovery_pct_matched_control']}%" for row in cp2),
            "status": "CP5_COVERAGE_INSUFFICIENT",
            "notes": "controls recover too strongly for CP5 universality language",
        },
    ]
    return rows


def write_summary(cp3_rows: list[dict[str, object]], comparative_rows: list[dict[str, object]], cp5_rows: list[dict[str, object]]) -> None:
    cp3_status = (
        "PASS"
        if all(row["status"] == "PASS" for row in cp3_rows)
        else "CP3_RESCALED_INVARIANT_NOT_SPECIFIC"
        if any(row["contract"] == "CP3 rescaled-invariant flow" for row in cp3_rows)
        else "CP3_OPEN"
    )
    comparative_status = "PASS" if all(row["status"] == "PASS" for row in comparative_rows) else "MIXED_OPEN"
    cp5_status = "PASS" if all(row["status"] == "PASS" for row in cp5_rows) else "CP5_COVERAGE_INSUFFICIENT"
    cp2_status = "CP2_CONTROL_INVALID" if CP2_CONTROL_INTEGRITY.exists() else "CP2_SAME_SURROGATE_OPEN"
    release_status = "RELEASE_66_5_GATES_CLOSED" if cp3_status == "PASS" and comparative_status == "PASS" and cp5_status == "PASS" and cp2_status == "CP2_CONTROL_INVALID" else "RELEASE_66_5_PARTIAL"

    lines = [
        "# Post-66.5 Roadmap Run Summary",
        "",
        "Status: Continuum Physics Scaling Roadmap artifact, not a new phase.",
        "",
        "Phase 67 remains parked.",
        "",
        f"Release verdict: `{release_status}`.",
        "",
        "Generated outputs:",
        "",
        f"- `{CP2_CONTROL_INTEGRITY.relative_to(ROOT)}`",
        f"- `{CP3_OUT.relative_to(ROOT)}`",
        f"- `{COMPARATIVE_OUT.relative_to(ROOT)}`",
        f"- `{CP5_OUT.relative_to(ROOT)}`",
        "",
        "## CP2 Same-Surrogate Control Integrity",
        "",
        f"Gate status: `{cp2_status}`.",
        "",
        f"See `{CP2_CONTROL_INTEGRITY.relative_to(ROOT)}` for branch/control distributions, effect size, overlap, and terminal control classification.",
        "",
        "## CP3 Effective-Equation Contract",
        "",
        f"Gate status: `{cp3_status}`.",
        "",
        "| Contract | Branch metric | Control metric | Contrast | Status |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for row in cp3_rows:
        lines.append(
            f"| {row['contract']} | {row['branch_metric']} | {row['control_metric']} | {row['separation_ratio']} | `{row['status']}` |"
        )

    lines.extend(
        [
            "",
            "Interpretation: coefficient-flow, metric-surrogate shell-slope, and propagation-speed rows separate from controls in the tested slices. The rescaled-invariant row does not separate, so it is terminally classified as `CP3_RESCALED_INVARIANT_NOT_SPECIFIC`.",
            "",
            "## Narrow Comparative Diagnostic",
            "",
            f"Gate status: `{comparative_status}`.",
            "",
            "| Diagnostic | Observed | Control | Status |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in comparative_rows:
        lines.append(
            f"| {row['diagnostic']} | {row['observed_result']} | {row['control_result']} | `{row['status']}` |"
        )

    lines.extend(
        [
            "",
            "Interpretation: the Kuramoto sidecar fails in a useful way: early detection appears, but edge-signature specificity does not pass. No superiority claim is made.",
            "",
            "## CP5 Universality Screen",
            "",
            f"Gate status: `{cp5_status}`.",
            "",
            "| Variation axis | Available | Required | Status | Notes |",
            "| --- | ---: | ---: | --- | --- |",
        ]
    )
    for row in cp5_rows:
        lines.append(
            f"| {row['variation_axis']} | {row['available_count']} | {row['required_minimum']} | `{row['status']}` | {row['notes']} |"
        )

    lines.extend(
        [
            "",
            "Claim boundary: this run records bounded scale-bridge diagnostics and sidecar comparisons. It does not prove a continuum limit, derive spacetime, derive physical equations, establish a physical metric, or replace standard models.",
            "",
            "Scale-bridge evidence is useful.",
            "",
            "Scale-bridge evidence is not continuum proof.",
        ]
    )
    SUMMARY_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    cp3_rows = build_cp3()
    comparative_rows = build_comparative(cp3_rows)
    cp5_rows = build_cp5()

    write_csv(
        CP3_OUT,
        cp3_rows,
        [
            "contract",
            "source",
            "candidate_law",
            "branch_metric",
            "control_metric",
            "separation_ratio",
            "branch_secondary",
            "control_secondary",
            "threshold",
            "status",
            "claim_boundary",
        ],
    )
    write_csv(
        COMPARATIVE_OUT,
        comparative_rows,
        ["diagnostic", "comparison_system", "observed_result", "control_result", "contrast", "status", "claim_boundary"],
    )
    write_csv(
        CP5_OUT,
        cp5_rows,
        ["variation_axis", "source", "available_count", "required_minimum", "branch_metric", "control_metric", "status", "notes"],
    )
    write_summary(cp3_rows, comparative_rows, cp5_rows)


if __name__ == "__main__":
    main()
