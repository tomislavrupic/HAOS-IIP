#!/usr/bin/env python3
"""Hydrogen spectra computational reference probe.

This sidecar runs a narrow gross-spectrum reference check:

- a declared hydrogen Rydberg-series catalog must reproduce its own line table;
- controls preserve line counts, ranges, ordering, or exact wavelength sets
  while breaking the transition law, absolute scale, identity, or selection rule;
- fine structure, Lamb shift, hyperfine structure, apparatus data, and CST/HAOS
  derivation remain outside this probe.

It is not a laboratory spectrum and not a claim that HAOS-IIP derives hydrogen.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[2]
CONTRACT_PATH = ROOT / "precommitment_contract.json"
CATALOG_PATH = ROOT / "transition_catalog.json"
CSV_PATH = ROOT / "hydrogen_spectra_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "hydrogen_spectra_reference_report.md"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)

FIELDNAMES = [
    "sample",
    "sample_role",
    "model_family",
    "data_kind",
    "series",
    "n_lower",
    "n_upper",
    "l_lower",
    "l_upper",
    "wavenumber_m_inv",
    "wavelength_nm",
    "reference_wavenumber_m_inv",
    "reference_wavelength_nm",
    "relative_wavenumber_error",
    "relative_wavelength_error",
    "fitted_rydberg_m_inv",
    "rydberg_relative_error",
    "transition_rms_relative_error",
    "series_law_rms_relative_error",
    "absolute_line_score",
    "rydberg_scale_score",
    "inverse_square_law_score",
    "line_order_score",
    "selection_rule_score",
    "hydrogen_reference_score",
    "status",
]

SERIES_SPECS = (
    ("lyman", 1, tuple(range(2, 7))),
    ("balmer", 2, tuple(range(3, 8))),
    ("paschen", 3, tuple(range(4, 9))),
)


@dataclass(frozen=True)
class HydrogenConfig:
    rydberg_h_m_inv: float = 10_967_758.340280
    pass_max_reference_relative_error: float = 1.0e-12
    pass_max_fit_rydberg_relative_error: float = 1.0e-12
    pass_max_series_law_rms_relative_error: float = 1.0e-12
    pass_min_line_order_score: float = 1.0
    pass_min_selection_rule_score: float = 1.0
    pass_min_reference_score: float = 0.999999
    control_false_positive_score_gate: float = 0.95
    line_score_error_scale: float = 0.02
    rydberg_score_error_scale: float = 0.02
    law_score_error_scale: float = 0.01
    scaled_rydberg_factor: float = 0.92


@dataclass(frozen=True)
class TransitionLine:
    series: str
    n_lower: int
    n_upper: int
    l_lower: int
    l_upper: int
    reference_wavenumber_m_inv: float
    reference_wavelength_nm: float
    wavenumber_m_inv: float
    wavelength_nm: float


@dataclass(frozen=True)
class SampleSummary:
    sample: str
    sample_role: str
    model_family: str
    transition_count: int
    fitted_rydberg_m_inv: float
    rydberg_relative_error: float
    transition_rms_relative_error: float
    max_transition_relative_error: float
    series_law_rms_relative_error: float
    absolute_line_score: float
    rydberg_scale_score: float
    inverse_square_law_score: float
    line_order_score: float
    selection_rule_score: float
    hydrogen_reference_score: float
    status: str


@dataclass(frozen=True)
class DiagnosticRow:
    sample: str
    sample_role: str
    model_family: str
    data_kind: str
    series: str
    n_lower: int | None
    n_upper: int | None
    l_lower: int | None
    l_upper: int | None
    wavenumber_m_inv: float | None
    wavelength_nm: float | None
    reference_wavenumber_m_inv: float | None
    reference_wavelength_nm: float | None
    relative_wavenumber_error: float | None
    relative_wavelength_error: float | None
    fitted_rydberg_m_inv: float | None
    rydberg_relative_error: float | None
    transition_rms_relative_error: float | None
    series_law_rms_relative_error: float | None
    absolute_line_score: float | None
    rydberg_scale_score: float | None
    inverse_square_law_score: float | None
    line_order_score: float | None
    selection_rule_score: float | None
    hydrogen_reference_score: float | None
    status: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run hydrogen spectra computational reference probe.")
    parser.add_argument("--rydberg-h", type=float, default=HydrogenConfig.rydberg_h_m_inv)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated contract, catalog, CSV, status JSON, and markdown report.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> HydrogenConfig:
    return HydrogenConfig(rydberg_h_m_inv=float(args.rydberg_h))


def stable_json_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def transition_delta(n_lower: int, n_upper: int) -> float:
    return 1.0 / float(n_lower * n_lower) - 1.0 / float(n_upper * n_upper)


def wavenumber_to_wavelength_nm(wavenumber_m_inv: float) -> float:
    return 1.0e9 / max(wavenumber_m_inv, 1.0e-300)


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def rms(values: list[float]) -> float:
    if not values:
        return 0.0
    return math.sqrt(sum(value * value for value in values) / float(len(values)))


def mean_abs(values: list[float]) -> float:
    if not values:
        return 0.0
    return sum(abs(value) for value in values) / float(len(values))


def build_reference_catalog(config: HydrogenConfig) -> list[TransitionLine]:
    lines: list[TransitionLine] = []
    for series, n_lower, upper_values in SERIES_SPECS:
        for n_upper in upper_values:
            wavenumber = config.rydberg_h_m_inv * transition_delta(n_lower, n_upper)
            lines.append(
                TransitionLine(
                    series=series,
                    n_lower=n_lower,
                    n_upper=n_upper,
                    l_lower=0,
                    l_upper=1,
                    reference_wavenumber_m_inv=wavenumber,
                    reference_wavelength_nm=wavenumber_to_wavelength_nm(wavenumber),
                    wavenumber_m_inv=wavenumber,
                    wavelength_nm=wavenumber_to_wavelength_nm(wavenumber),
                )
            )
    return lines


def replace_wavenumber(line: TransitionLine, wavenumber_m_inv: float) -> TransitionLine:
    return TransitionLine(
        series=line.series,
        n_lower=line.n_lower,
        n_upper=line.n_upper,
        l_lower=line.l_lower,
        l_upper=line.l_upper,
        reference_wavenumber_m_inv=line.reference_wavenumber_m_inv,
        reference_wavelength_nm=line.reference_wavelength_nm,
        wavenumber_m_inv=float(wavenumber_m_inv),
        wavelength_nm=wavenumber_to_wavelength_nm(float(wavenumber_m_inv)),
    )


def replace_l_values(line: TransitionLine, l_lower: int, l_upper: int) -> TransitionLine:
    return TransitionLine(
        series=line.series,
        n_lower=line.n_lower,
        n_upper=line.n_upper,
        l_lower=int(l_lower),
        l_upper=int(l_upper),
        reference_wavenumber_m_inv=line.reference_wavenumber_m_inv,
        reference_wavelength_nm=line.reference_wavelength_nm,
        wavenumber_m_inv=line.wavenumber_m_inv,
        wavelength_nm=line.wavelength_nm,
    )


def build_arithmetic_spacing_control(reference: list[TransitionLine]) -> list[TransitionLine]:
    control: list[TransitionLine] = []
    by_series = group_by_series(reference)
    for _series, lines in by_series.items():
        first = lines[0].reference_wavenumber_m_inv
        last = lines[-1].reference_wavenumber_m_inv
        count = len(lines)
        for index, line in enumerate(lines):
            fraction = 0.0 if count == 1 else float(index) / float(count - 1)
            control.append(replace_wavenumber(line, first + fraction * (last - first)))
    return control


def build_wrong_exponent_control(reference: list[TransitionLine]) -> list[TransitionLine]:
    wrong_basis = [1.0 / float(line.n_lower) - 1.0 / float(line.n_upper) for line in reference]
    references = [line.reference_wavenumber_m_inv for line in reference]
    numerator = sum(basis * value for basis, value in zip(wrong_basis, references))
    denominator = max(sum(basis * basis for basis in wrong_basis), 1.0e-300)
    scale = numerator / denominator
    return [replace_wavenumber(line, scale * basis) for line, basis in zip(reference, wrong_basis)]


def build_scaled_rydberg_control(reference: list[TransitionLine], config: HydrogenConfig) -> list[TransitionLine]:
    return [replace_wavenumber(line, config.scaled_rydberg_factor * line.reference_wavenumber_m_inv) for line in reference]


def build_reversed_transition_identity_control(reference: list[TransitionLine]) -> list[TransitionLine]:
    control: list[TransitionLine] = []
    by_series = group_by_series(reference)
    for _series, lines in by_series.items():
        reversed_wavenumbers = [line.reference_wavenumber_m_inv for line in reversed(lines)]
        for line, wavenumber in zip(lines, reversed_wavenumbers):
            control.append(replace_wavenumber(line, wavenumber))
    return control


def build_forbidden_selection_rule_control(reference: list[TransitionLine]) -> list[TransitionLine]:
    return [replace_l_values(line, 0, 0) for line in reference]


def group_by_series(lines: list[TransitionLine]) -> dict[str, list[TransitionLine]]:
    grouped: dict[str, list[TransitionLine]] = {}
    for line in lines:
        grouped.setdefault(line.series, []).append(line)
    for series_lines in grouped.values():
        series_lines.sort(key=lambda item: item.n_upper)
    return grouped


def fit_rydberg(lines: list[TransitionLine]) -> float:
    deltas = [transition_delta(line.n_lower, line.n_upper) for line in lines]
    values = [line.wavenumber_m_inv for line in lines]
    numerator = sum(delta * value for delta, value in zip(deltas, values))
    denominator = max(sum(delta * delta for delta in deltas), 1.0e-300)
    return numerator / denominator


def line_order_score(lines: list[TransitionLine]) -> float:
    checks = 0
    passed = 0
    for series_lines in group_by_series(lines).values():
        for left, right in zip(series_lines, series_lines[1:]):
            checks += 1
            if right.wavenumber_m_inv > left.wavenumber_m_inv:
                passed += 1
    if checks == 0:
        return 1.0
    return float(passed) / float(checks)


def selection_rule_score(lines: list[TransitionLine]) -> float:
    if not lines:
        return 0.0
    passed = sum(1 for line in lines if abs(line.l_upper - line.l_lower) == 1)
    return float(passed) / float(len(lines))


def summarize_sample(
    sample: str,
    sample_role: str,
    model_family: str,
    lines: list[TransitionLine],
    config: HydrogenConfig,
) -> SampleSummary:
    fitted_rydberg = fit_rydberg(lines)
    fitted_values = [
        fitted_rydberg * transition_delta(line.n_lower, line.n_upper)
        for line in lines
    ]
    observed_values = [line.wavenumber_m_inv for line in lines]
    reference_values = [line.reference_wavenumber_m_inv for line in lines]
    line_errors = [
        abs(observed - reference) / max(abs(reference), 1.0e-300)
        for observed, reference in zip(observed_values, reference_values)
    ]
    law_residuals = [
        (observed - fitted) / max(mean_abs(observed_values), 1.0e-300)
        for observed, fitted in zip(observed_values, fitted_values)
    ]
    transition_rms = rms(line_errors)
    series_law_rms = rms(law_residuals)
    rydberg_relative_error = abs(fitted_rydberg - config.rydberg_h_m_inv) / max(abs(config.rydberg_h_m_inv), 1.0e-300)
    absolute_line_score = clamp01(1.0 - transition_rms / config.line_score_error_scale)
    rydberg_scale_score = clamp01(1.0 - rydberg_relative_error / config.rydberg_score_error_scale)
    inverse_square_law_score = clamp01(1.0 - series_law_rms / config.law_score_error_scale)
    order_score = line_order_score(lines)
    rule_score = selection_rule_score(lines)
    reference_score = min(
        absolute_line_score,
        rydberg_scale_score,
        inverse_square_law_score,
        order_score,
        rule_score,
    )

    if sample_role == "reference":
        status = (
            "REFERENCE_RYDBERG_SERIES_REPRODUCED"
            if (
                max(line_errors, default=0.0) <= config.pass_max_reference_relative_error
                and rydberg_relative_error <= config.pass_max_fit_rydberg_relative_error
                and series_law_rms <= config.pass_max_series_law_rms_relative_error
                and order_score >= config.pass_min_line_order_score
                and rule_score >= config.pass_min_selection_rule_score
                and reference_score >= config.pass_min_reference_score
            )
            else "REFERENCE_FAILURE"
        )
    else:
        status = (
            "CONTROL_FALSE_POSITIVE"
            if reference_score >= config.control_false_positive_score_gate
            else "CONTROL_REJECTED"
        )

    return SampleSummary(
        sample=sample,
        sample_role=sample_role,
        model_family=model_family,
        transition_count=len(lines),
        fitted_rydberg_m_inv=fitted_rydberg,
        rydberg_relative_error=rydberg_relative_error,
        transition_rms_relative_error=transition_rms,
        max_transition_relative_error=max(line_errors, default=0.0),
        series_law_rms_relative_error=series_law_rms,
        absolute_line_score=absolute_line_score,
        rydberg_scale_score=rydberg_scale_score,
        inverse_square_law_score=inverse_square_law_score,
        line_order_score=order_score,
        selection_rule_score=rule_score,
        hydrogen_reference_score=reference_score,
        status=status,
    )


def make_rows(
    sample: str,
    sample_role: str,
    model_family: str,
    lines: list[TransitionLine],
    summary: SampleSummary,
) -> list[DiagnosticRow]:
    rows: list[DiagnosticRow] = [
        DiagnosticRow(
            sample=sample,
            sample_role=sample_role,
            model_family=model_family,
            data_kind="aggregate",
            series="all",
            n_lower=None,
            n_upper=None,
            l_lower=None,
            l_upper=None,
            wavenumber_m_inv=None,
            wavelength_nm=None,
            reference_wavenumber_m_inv=None,
            reference_wavelength_nm=None,
            relative_wavenumber_error=None,
            relative_wavelength_error=None,
            fitted_rydberg_m_inv=summary.fitted_rydberg_m_inv,
            rydberg_relative_error=summary.rydberg_relative_error,
            transition_rms_relative_error=summary.transition_rms_relative_error,
            series_law_rms_relative_error=summary.series_law_rms_relative_error,
            absolute_line_score=summary.absolute_line_score,
            rydberg_scale_score=summary.rydberg_scale_score,
            inverse_square_law_score=summary.inverse_square_law_score,
            line_order_score=summary.line_order_score,
            selection_rule_score=summary.selection_rule_score,
            hydrogen_reference_score=summary.hydrogen_reference_score,
            status=summary.status,
        )
    ]
    for line in lines:
        rel_wave = abs(line.wavelength_nm - line.reference_wavelength_nm) / max(abs(line.reference_wavelength_nm), 1.0e-300)
        rel_wavenumber = abs(line.wavenumber_m_inv - line.reference_wavenumber_m_inv) / max(abs(line.reference_wavenumber_m_inv), 1.0e-300)
        rows.append(
            DiagnosticRow(
                sample=sample,
                sample_role=sample_role,
                model_family=model_family,
                data_kind="line",
                series=line.series,
                n_lower=line.n_lower,
                n_upper=line.n_upper,
                l_lower=line.l_lower,
                l_upper=line.l_upper,
                wavenumber_m_inv=line.wavenumber_m_inv,
                wavelength_nm=line.wavelength_nm,
                reference_wavenumber_m_inv=line.reference_wavenumber_m_inv,
                reference_wavelength_nm=line.reference_wavelength_nm,
                relative_wavenumber_error=rel_wavenumber,
                relative_wavelength_error=rel_wave,
                fitted_rydberg_m_inv=None,
                rydberg_relative_error=None,
                transition_rms_relative_error=None,
                series_law_rms_relative_error=None,
                absolute_line_score=None,
                rydberg_scale_score=None,
                inverse_square_law_score=None,
                line_order_score=None,
                selection_rule_score=None,
                hydrogen_reference_score=None,
                status="LINE_RECORDED",
            )
        )
    return rows


def precommitment_contract(config: HydrogenConfig) -> dict[str, Any]:
    return {
        "name": "Hydrogen spectra computational reference probe",
        "status": "claim_gated_sidecar",
        "scope": {
            "implemented_fact": "Run a deterministic gross Rydberg-series reference and component controls.",
            "design_choice": "Use declared R_H and the gross line relation sigma = R_H*(1/n_lower^2 - 1/n_upper^2).",
            "heuristic": "Scores decompose absolute line match, fitted Rydberg scale, inverse-square law residual, line order, and delta_l selection-rule checks.",
            "unverified_hypothesis": "No CST or HAOS hydrogen-spectra recovery hypothesis is tested here.",
        },
        "non_claims": [
            "not a laboratory hydrogen spectrum",
            "not a fine-structure calculation",
            "not a Lamb-shift or hyperfine calculation",
            "not a QED derivation",
            "not a CST derivation of hydrogen spectra",
            "not evidence that HAOS-IIP recovers hydrogen spectra",
        ],
        "reference_constant": {
            "name": "declared hydrogen Rydberg constant",
            "symbol": "R_H",
            "value_m_inv": config.rydberg_h_m_inv,
            "note": "This sidecar uses a declared calibration constant; it does not estimate it from HAOS-IIP.",
        },
        "series_catalog": [
            {"series": series, "n_lower": n_lower, "n_upper_values": list(upper_values)}
            for series, n_lower, upper_values in SERIES_SPECS
        ],
        "control_parameters": {
            "scaled_rydberg_factor": config.scaled_rydberg_factor,
        },
        "controls": {
            "arithmetic_spacing_control": {
                "preserve": ["series labels", "line count", "first and last wavenumber per series", "monotonic order"],
                "destroy": ["inverse-square transition spacing"],
                "expected_response": "inverse_square_law_score decreases; absolute_line_score decreases for interior lines",
            },
            "wrong_exponent_control": {
                "preserve": ["transition labels", "monotonic line families", "rough spectral scale"],
                "destroy": ["1/n^2 energy law"],
                "expected_response": "inverse_square_law_score and absolute_line_score decrease",
            },
            "scaled_rydberg_control": {
                "preserve": ["inverse-square law", "line ordering", "selection labels"],
                "destroy": ["absolute Rydberg scale"],
                "expected_response": "rydberg_scale_score and absolute_line_score decrease while inverse_square_law_score remains high",
            },
            "reversed_transition_identity_control": {
                "preserve": ["exact wavenumber set within each series", "series labels", "selection labels"],
                "destroy": ["transition-to-line identity", "line ordering"],
                "expected_response": "line_order_score and absolute_line_score decrease",
            },
            "forbidden_selection_rule_control": {
                "preserve": ["gross wavelengths", "transition n labels", "inverse-square law"],
                "destroy": ["dipole delta_l selection rule"],
                "expected_response": "selection_rule_score decreases while spectral-line scores remain high",
            },
        },
        "thresholds": {
            "pass_max_reference_relative_error": config.pass_max_reference_relative_error,
            "pass_max_fit_rydberg_relative_error": config.pass_max_fit_rydberg_relative_error,
            "pass_max_series_law_rms_relative_error": config.pass_max_series_law_rms_relative_error,
            "pass_min_line_order_score": config.pass_min_line_order_score,
            "pass_min_selection_rule_score": config.pass_min_selection_rule_score,
            "pass_min_reference_score": config.pass_min_reference_score,
            "control_false_positive_score_gate": config.control_false_positive_score_gate,
            "line_score_error_scale": config.line_score_error_scale,
            "rydberg_score_error_scale": config.rydberg_score_error_scale,
            "law_score_error_scale": config.law_score_error_scale,
        },
        "verdict_logic": {
            "reference_pass": [
                "reference line max relative error <= pass_max_reference_relative_error",
                "fitted R_H relative error <= pass_max_fit_rydberg_relative_error",
                "series law RMS relative error <= pass_max_series_law_rms_relative_error",
                "line_order_score >= pass_min_line_order_score",
                "selection_rule_score >= pass_min_selection_rule_score",
                "all controls score < control_false_positive_score_gate",
            ],
            "always_report": [
                "CST_HYDROGEN_STATUS_OPEN",
                "HAOS_DERIVATION_NOT_TESTED",
                "FINE_STRUCTURE_OPEN_NOT_MODELED",
            ],
        },
        "future_candidate_rejection_conditions": [
            "candidate imports the Rydberg formula or declared R_H by construction",
            "candidate fits or tunes R_H after seeing the reference lines",
            "candidate passes gross wavelengths while failing declared selection-rule checks",
            "candidate passes only by using transition labels unavailable to controls",
            "candidate claims fine structure, Lamb shift, hyperfine structure, QED, or apparatus recovery from this gross-spectrum probe",
        ],
    }


def sample_definitions(reference: list[TransitionLine], config: HydrogenConfig) -> list[tuple[str, str, str, list[TransitionLine]]]:
    return [
        ("hydrogen_rydberg_reference", "reference", "gross_rydberg_series", reference),
        ("arithmetic_spacing_control", "negative_control", "line_count_range_control", build_arithmetic_spacing_control(reference)),
        ("wrong_exponent_control", "negative_control", "wrong_transition_law_control", build_wrong_exponent_control(reference)),
        ("scaled_rydberg_control", "negative_control", "scale_shift_control", build_scaled_rydberg_control(reference, config)),
        ("reversed_transition_identity_control", "negative_control", "identity_shuffle_control", build_reversed_transition_identity_control(reference)),
        ("forbidden_selection_rule_control", "negative_control", "selection_rule_control", build_forbidden_selection_rule_control(reference)),
    ]


def row_dict(row: DiagnosticRow) -> dict[str, Any]:
    payload = asdict(row)
    for key, value in list(payload.items()):
        if isinstance(value, float):
            payload[key] = f"{value:.12g}"
        elif value is None:
            payload[key] = ""
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[DiagnosticRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row_dict(row))


def run_probe(config: HydrogenConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = precommitment_contract(config)
    contract["contract_hash"] = stable_json_hash("hydrogen_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)

    reference = build_reference_catalog(config)
    write_json(output_dir / CATALOG_PATH.name, [asdict(line) for line in reference])

    summaries: list[SampleSummary] = []
    rows: list[DiagnosticRow] = []
    for sample, sample_role, model_family, lines in sample_definitions(reference, config):
        summary = summarize_sample(sample, sample_role, model_family, lines, config)
        summaries.append(summary)
        rows.extend(make_rows(sample, sample_role, model_family, lines, summary))
    write_csv(output_dir / CSV_PATH.name, rows)

    reference_summary = next(summary for summary in summaries if summary.sample_role == "reference")
    control_false_positive_count = sum(
        1 for summary in summaries if summary.sample_role == "negative_control" and summary.status == "CONTROL_FALSE_POSITIVE"
    )
    reference_pass = (
        reference_summary.status == "REFERENCE_RYDBERG_SERIES_REPRODUCED"
        and control_false_positive_count == 0
    )
    labels = [
        "HYDROGEN_REFERENCE_SANITY_PASS" if reference_pass else "HYDROGEN_REFERENCE_SANITY_FAIL",
        "CST_HYDROGEN_STATUS_OPEN",
        "HAOS_DERIVATION_NOT_TESTED",
        "FINE_STRUCTURE_OPEN_NOT_MODELED",
        "NO_PHYSICAL_EXPERIMENT",
    ]
    result: dict[str, Any] = {
        "labels": labels,
        "contract_hash": contract["contract_hash"],
        "config": asdict(config),
        "reference_constant": contract["reference_constant"],
        "summaries": [asdict(summary) for summary in summaries],
        "control_false_positive_count": control_false_positive_count,
        "fine_structure_status": "OPEN_NOT_MODELED",
        "cst_mechanism_slot": "EMPTY_NOT_TESTED",
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "transition_catalog": repo_rel(output_dir / CATALOG_PATH.name),
            "hydrogen_spectra_diagnostics": repo_rel(output_dir / CSV_PATH.name),
            "bridge_status": repo_rel(output_dir / STATUS_PATH.name),
            "hydrogen_spectra_reference_report": repo_rel(output_dir / REPORT_PATH.name),
        },
        "limitations": contract["non_claims"],
    }
    result["result_hash"] = stable_json_hash("hydrogen_spectra_result", result)
    write_json(output_dir / STATUS_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, summaries)
    return result


def write_report(path: Path, result: dict[str, Any], summaries: list[SampleSummary]) -> None:
    reference = next(summary for summary in summaries if summary.sample_role == "reference")
    controls = [summary for summary in summaries if summary.sample_role == "negative_control"]
    lines = [
        "# Hydrogen Spectra Computational Reference Probe",
        "",
        "Implemented fact: this sidecar runs a deterministic gross Rydberg-series reference and component controls.",
        "Design choice: gross line positions use the declared relation `sigma = R_H * (1/n_lower^2 - 1/n_upper^2)`.",
        "Heuristic: the reference score is the minimum of line match, fitted scale, inverse-square residual, line order, and selection-rule scores.",
        "Unverified hypothesis: no CST or HAOS hydrogen-spectra recovery hypothesis is tested here.",
        "",
        "## Verdict Labels",
    ]
    lines.extend(f"- {label}" for label in result["labels"])
    lines.extend(
        [
            "",
            "## Reference",
            f"- declared R_H: `{result['reference_constant']['value_m_inv']:.12f} m^-1`",
            f"- transition count: `{reference.transition_count}`",
            f"- fitted R_H: `{reference.fitted_rydberg_m_inv:.12f} m^-1`",
            f"- max transition relative error: `{reference.max_transition_relative_error:.12g}`",
            f"- reference score: `{reference.hydrogen_reference_score:.12f}`",
            "",
            "## Controls",
        ]
    )
    for control in controls:
        lines.append(
            "- {sample}: score `{score:.12f}`, line `{line:.12f}`, scale `{scale:.12f}`, law `{law:.12f}`, order `{order:.12f}`, selection `{selection:.12f}`, status `{status}`".format(
                sample=control.sample,
                score=control.hydrogen_reference_score,
                line=control.absolute_line_score,
                scale=control.rydberg_scale_score,
                law=control.inverse_square_law_score,
                order=control.line_order_score,
                selection=control.selection_rule_score,
                status=control.status,
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "- This is not a laboratory hydrogen spectrum.",
            "- This does not model fine structure, Lamb shift, hyperfine splitting, line intensities, or apparatus effects.",
            "- This does not derive hydrogen spectra from CST or HAOS-IIP.",
            "- This does not change the frozen CST v0.2.2 checkpoint.",
            "- `HYDROGEN_REFERENCE_SANITY_PASS` means only that the computational reference harness behaves as expected.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_probe(make_config(args), args.output_dir)
    print(json.dumps({"labels": result["labels"], "result_hash": result["result_hash"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
